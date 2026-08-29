#!/usr/bin/env python3
"""tier1_scan.py -- headline verification under the IDENTIFIED original protocol.

Protocol (see window_protocol_IDENTIFIED_20260814.md):
  GL folding = narrow classes mod <inverse, det-1 conjugation>
  invariant  = full-space cycle type on X(F_q) = F_q^3   (cycle_type_on_Fq)
  pairs      = all pairs; p_min(D) = max over pairs of least witness
  primes     : two runs -- 'unram' (original) and 'all' (restriction lifted)

Fast class enumeration: for each b < sqrt(D) with b^2 = D mod 4, set
N = (D - b^2)/4 > 0; then c = -N/a, so a runs over (signed) divisors of N.
Validated against markoff_rigidity.classes().
"""
import os, sys
_R = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [_R, os.path.join(_R, "src"), os.path.dirname(os.path.abspath(__file__))]

import sys, json, time
from math import isqrt, gcd
from collections import Counter
import markoff_rigidity as mr
from markoff_rigidity import (canon, automorph, decomp, cycle_type_on_Fq,
                              kronecker_split, _is_reduced)

QLIST = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31]

def divisors(n):
    ds = []
    i = 1
    while i * i <= n:
        if n % i == 0:
            ds.append(i)
            if i != n // i: ds.append(n // i)
        i += 1
    return ds

def classes_fast(D):
    """All narrow ideal classes (canonical reduced reps), primitive forms."""
    seen = set(); reps = []
    r = isqrt(D)
    for b in range(1, r + 2):
        if (b * b - D) % 4: continue
        N = (D - b * b) // 4
        if N <= 0: continue
        for d in divisors(N):
            for a in (d, -d):
                c = -(N // d) if a > 0 else (N // d)
                f = (a, b, c)
                if a * c * 4 != b * b - D: continue
                if gcd(gcd(abs(a), abs(b)), abs(c)) != 1: continue
                if not _is_reduced(f, D): continue
                k = canon(f, D)
                if k not in seen:
                    seen.add(k); reps.append(k)
    return reps

def gl_reps(T, D, cls):
    def mcls(M):
        (a, b), (c, d) = M
        return canon((c, d - a, -b), D)
    reps = []; seen = set()
    for c0 in cls:
        if c0 in seen: continue
        M = automorph(c0, T, 1); (a, b), (c, d) = M
        seen |= {mcls(M), mcls([[d, -b], [-c, a]]),
                 mcls([[a, -b], [-c, d]]), mcls([[d, b], [c, a]])}
        reps.append(c0)
    return reps

def scan(t0, t1, mode):
    dist = Counter(); worst = (0, None); rows = {}
    for T in range(t0, t1 + 1):
        D = T * T - 4
        cls = classes_fast(D)
        reps = gl_reps(T, D, cls)
        n = len(reps)
        if n < 2:
            continue
        words = {c: decomp(automorph(c, T, 1)) for c in reps}
        qs = [q for q in QLIST if mode == 'all' or kronecker_split(D, q) != 'ram']
        cache = {}
        def inv(c, q):
            if (c, q) not in cache:
                cache[(c, q)] = cycle_type_on_Fq(words[c], q)
            return cache[(c, q)]
        pm = 0; failed = False
        for i in range(n):
            for j in range(i + 1, n):
                w = None
                for q in qs:
                    if inv(reps[i], q) != inv(reps[j], q):
                        w = q; break
                if w is None: failed = True; w = 0
                pm = max(pm, w)
        key = 'FAIL' if failed else pm
        dist[key] += 1
        rows[T] = {'nGL': n, 'p_min': key}
        if not failed and pm > worst[0]: worst = (pm, T)
        if T % 50 == 0:
            print(f"  ...T={T} (nGL={n}, p_min={key}) running dist={dict(dist)}",
                  flush=True)
    return dist, worst, rows

if __name__ == '__main__':
    t0, t1 = int(sys.argv[1]), int(sys.argv[2])
    # validation of fast enumeration
    for D in (221, 1517, 4620, 57596, 62496):
        a = sorted(classes_fast(D)); b = sorted(mr.classes(D))
        assert a == b, (D, len(a), len(b))
    print("classes_fast validated against markoff_rigidity.classes()", flush=True)
    out = {}
    for mode in ('unram', 'all'):
        st = time.time()
        d, w, rows = scan(t0, t1, mode)
        print(f"\n=== mode={mode}  T in [{t0},{t1}] ===")
        print("   distribution:", dict(sorted(d.items(), key=lambda kv: str(kv[0]))))
        print(f"   max p_min = {w[0]} at T={w[1]}   ({round(time.time()-st,1)}s)",
              flush=True)
        out[mode] = {'dist': {str(k): v for k, v in d.items()},
                     'max': w, 'rows': rows}
    json.dump(out, open(os.path.join(_R,'data',f'tier1_{t0}_{t1}.json'),'w'))
