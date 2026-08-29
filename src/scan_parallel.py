#!/usr/bin/env python3
"""scan_parallel.py -- exhaustive search for simultaneously blind traces.

Portable, parallel, resumable.  Intended to be run OFF this container, on a
multicore machine or cluster node.

    python3 scan_parallel.py --lo 3 --hi 1021022 --workers 64 --out results/

Dependencies: numpy, and the project's markoff_rigidity.py on the import path
(only canon / automorph / decomp are used, so the arithmetic conventions stay
identical to the published code).

Definitions (matching the paper):
  * GL-class      : narrow ideal class of D = T^2 - 4, folded by
                    <inverse, det -1 conjugation>.
  * invariant     : cycle type of the Vieta permutation of X(F_q) = F_q^3
                    induced by the automorph of the class (full space).
  * blind at S    : some pair of GL-classes has equal invariants for every q in S.
  * p_min(D)      : least p with {3,...,p} separating all pairs.

Sieve order.  Blindness at S requires collision at EVERY q in S, so the primes
may be applied in any order; cost is minimised by applying cheap, discriminating
primes first.  Measured at T~18000 (349 classes): q=3 splits into only 3 values
(useless as a filter), q=5 into 279, q=7 into 347, q=11 into 349.  Applying
5 -> 7 -> 11 -> 13 -> 3 evaluates the expensive primes on a handful of survivors
only, and is ~1.6x faster than the naive increasing order.

Checkpointing: one JSON per chunk in --out; existing chunk files are skipped, so
the job can be killed and restarted at will.
"""
import argparse, json, os, sys, time
from math import isqrt, gcd
from collections import Counter
from multiprocessing import Pool

import numpy as np
from markoff_rigidity import canon, automorph, decomp, _is_reduced, _rho

SIEVE_ORDER = [3, 5, 7, 11, 13, 17, 19, 23]   # q=3 and q=5 share one fused pass

# ---------------------------------------------------------------- forms
def _reduce(f, D):
    g, guard = f, 0
    while not _is_reduced(g, D):
        g = _rho(g, D); guard += 1
        if guard > 100000: raise RuntimeError('reduction loop')
    return g

def _egcd(a, b):
    x0, x1, y0, y1 = 1, 0, 0, 1
    while b:
        q, a, b = a // b, b, a % b
        x0, x1 = x1, x0 - q * x1
        y0, y1 = y1, y0 - q * y1
    return a, x0, y0

def _compose(f1, f2, D):
    a1, b1, c1 = f1; a2, b2, c2 = f2
    if a1 > a2: a1, b1, c1, a2, b2, c2 = a2, b2, c2, a1, b1, c1
    s = (b1 + b2) // 2; n = b2 - s
    if a2 % a1 == 0: y1, d = 0, a1
    else: d, y1, _ = _egcd(a2, a1)
    if s % d == 0: y2, x2, d1 = -1, 0, d
    else:
        d1, x2, v = _egcd(s, d); y2 = -v
    v1, v2 = a1 // d1, a2 // d1
    r = (y1 * y2 * n - x2 * c2) % v1
    b3 = b2 + 2 * v2 * r; a3 = v1 * v2
    return _reduce((a3, b3, (b3 * b3 - D) // (4 * a3)), D)

def _sqrt_mod_p(a, p):
    a %= p
    if a == 0: return 0
    if pow(a, (p - 1) // 2, p) != 1: return None
    if p % 4 == 3: return pow(a, (p + 1) // 4, p)
    q, s = p - 1, 0
    while q % 2 == 0: q //= 2; s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1: z += 1
    m, c, t, r = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
    while t != 1:
        i, t2 = 0, t
        while t2 != 1: t2 = t2 * t2 % p; i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c = i, b * b % p
        t = t * c % p; r = r * b % p
    return r

_PR = None
def _primes(n):
    global _PR
    if _PR is None or _PR[-1] < n:
        sieve = bytearray([1]) * (n + 1); sieve[0:2] = b'\x00\x00'
        for i in range(2, isqrt(n) + 1):
            if sieve[i]: sieve[i * i::i] = bytearray(len(sieve[i * i::i]))
        _PR = [i for i in range(n + 1) if sieve[i]]
    return [p for p in _PR if p <= n]

def class_reps(D, pbound=3000, gcap=200, idle_max=40):
    """Narrow class group of discriminant D, as canonical reduced reps.

    Coset construction: O(h) compositions.  NOTE two easily-missed points:
    (i) prime forms all have a>0 and generate only an index-2 subgroup in many
        cases -- the NEGATIVE principal form is added as a generator;
    (ii) truncating the generator list too early loses classes (24 was not
        enough; 200 with idle_max=40 was verified exhaustively for T <= 2000).
    """
    b0 = D % 2
    pr = canon(_reduce((1, b0, (b0 * b0 - D) // 4), D), D)
    gens = [_reduce((-1, b0, (D - b0 * b0) // 4), D)]
    for p in _primes(pbound):
        if D % p == 0: continue
        r = _sqrt_mod_p(D, p)
        if r is None: continue
        for t in (r, p - r):
            b = t if (t - D) % 2 == 0 else t + p
            if (b * b - D) % (4 * p): continue
            gens.append(_reduce((p, b, (b * b - D) // (4 * p)), D))
        if len(gens) >= gcap: break
    forms = [_reduce((1, b0, (b0 * b0 - D) // 4), D)]; seen = {pr}; idle = 0
    for g in gens:
        if canon(g, D) in seen:
            idle += 1
            if idle >= idle_max: break
            continue
        idle = 0
        base = list(forms); cur = g
        while canon(cur, D) not in seen:
            for f in base:
                h = _compose(f, cur, D); k = canon(h, D)
                if k not in seen: seen.add(k); forms.append(h)
            cur = _compose(cur, g, D)
    return sorted(seen)

def gl_reps(T, D, cls):
    def mc(M):
        (a, b), (c, d) = M
        return canon((c, d - a, -b), D)
    reps, seen = [], set()
    for c0 in cls:
        if c0 in seen: continue
        M = automorph(c0, T, 1); (a, b), (c, d) = M
        seen |= {mc(M), mc([[d, -b], [-c, a]]),
                 mc([[a, -b], [-c, d]]), mc([[d, b], [c, a]])}
        reps.append(c0)
    return reps

# ---------------------------------------------------------------- invariants
def _step(g, x, y, z, q):
    if g == 'T':  return x, z, (x * z - y) % q
    if g == 'Ti': return x, (x * y - z) % q, y
    if g in ('S', 'Si'): return y, x, (x * y - z) % q
    return x, y, (x * y - z) % q

_LET = {}
def _tables(q):
    if q not in _LET:
        i = np.arange(q ** 3)
        x, y, z = i // (q * q), (i // q) % q, i % q
        d = {}
        for g in ('T', 'Ti', 'S', 'J'):
            a, b, c = _step(g, x, y, z, q)
            d[g] = ((a * q + b) * q + c).astype(np.int32)
        d['Si'] = d['S']
        _LET[q] = d
    return _LET[q]

def cycle_type(word, q):
    L = _tables(q)
    p = np.arange(q ** 3, dtype=np.int32)
    for g in word: p = L[g][p]
    pl = p.tolist(); n = q ** 3; seen = bytearray(n); ct = Counter()
    for s in range(n):
        if seen[s]: continue
        c, t = 0, s
        while not seen[t]: seen[t] = 1; c += 1; t = pl[t]
        ct[c] += 1
    return tuple(sorted(ct.items()))

_FUSE = {}
def _fused35():
    """One index space holding X(F_3) and X(F_5): a single numpy op per letter
    yields BOTH cycle types, so the (nearly useless but definition-critical)
    prime 3 costs essentially nothing and p_min stays exact."""
    if not _FUSE:
        tot = 27 + 125; arrs = {}
        for g in ('T', 'Ti', 'S', 'J'):
            a = np.empty(tot, dtype=np.int32)
            for q, off in ((3, 0), (5, 27)):
                a[off:off + q ** 3] = _tables(q)[g] + off
            arrs[g] = a
        arrs['Si'] = arrs['S']
        _FUSE.update(arrs)
    return _FUSE

def cycle_types_35(word):
    L = _fused35()
    p = np.arange(152, dtype=np.int32)
    for g in word: p = L[g][p]
    pl = p.tolist(); out = []
    for q, off in ((3, 0), (5, 27)):
        n = q ** 3; seen = bytearray(n); ct = Counter()
        for s in range(n):
            if seen[s]: continue
            c, t = 0, s
            while not seen[t]: seen[t] = 1; c += 1; t = pl[off + t] - off
            ct[c] += 1
        out.append(tuple(sorted(ct.items())))
    return out

def analyse(T, order=SIEVE_ORDER):
    """Exact p_min (least p with {3,...,p} separating all pairs), and nGL."""
    D = T * T - 4
    reps = gl_reps(T, D, class_reps(D))
    n = len(reps)
    if n < 2: return 0, n
    W = [decomp(automorph(c, T, 1)) for c in reps]
    inv35 = [None] * n
    def inv(i, q):
        if q in (3, 5):
            if inv35[i] is None: inv35[i] = cycle_types_35(W[i])
            return inv35[i][0 if q == 3 else 1]
        return cycle_type(W[i], q)
    buckets = [list(range(n))]
    for q in order:
        nb = []
        for bk in buckets:
            if len(bk) == 1: continue
            sub = {}
            for i in bk: sub.setdefault(inv(i, q), []).append(i)
            nb += [v for v in sub.values() if len(v) > 1]
        buckets = nb
        if not buckets: return q, n
    return None, n

def _chunk(args):
    lo, hi = args
    out = {'lo': lo, 'hi': hi, 'dist': {}, 'blind': {}}
    for T in range(lo, hi + 1):
        p, n = analyse(T)
        if n < 2: continue
        k = str(p) if p else 'gt23'
        out['dist'][k] = out['dist'].get(k, 0) + 1
        for thr, name in ((7, 'S357'), (11, 'S35711'), (13, 'S3571113')):
            if p is None or p > thr:
                out['blind'].setdefault(name, []).append([T, p, n])
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--lo', type=int, default=3)
    ap.add_argument('--hi', type=int, required=True)
    ap.add_argument('--workers', type=int, default=os.cpu_count())
    ap.add_argument('--chunk', type=int, default=500)
    ap.add_argument('--out', default='results')
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    jobs = []
    for lo in range(a.lo, a.hi + 1, a.chunk):
        hi = min(lo + a.chunk - 1, a.hi)
        if not os.path.exists(f'{a.out}/chunk_{lo}_{hi}.json'):
            jobs.append((lo, hi))
    print(f"{len(jobs)} chunks to do on {a.workers} workers", flush=True)
    t0 = time.time()
    with Pool(a.workers) as pool:
        for i, r in enumerate(pool.imap_unordered(_chunk, jobs), 1):
            json.dump(r, open(f"{a.out}/chunk_{r['lo']}_{r['hi']}.json", 'w'))
            if r['blind'].get('S35711') or r['blind'].get('S3571113'):
                print(f"  *** {r['lo']}-{r['hi']}: {r['blind']}", flush=True)
            if i % 20 == 0:
                el = time.time() - t0
                print(f"  {i}/{len(jobs)} chunks, {el/60:.1f} min, "
                      f"ETA {el/i*(len(jobs)-i)/3600:.1f} h", flush=True)

if __name__ == '__main__':
    main()
