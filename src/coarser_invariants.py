#!/usr/bin/env python3
"""Test candidate coarser invariants for the
frozen distribution {5:20, 7:18, 11:7, 13:1} on window 240<=T<=285.
Candidates (same GL folding as reproduce_window):
  A per-level cycle types            (reproduce baseline, known {5:37,7:9})
  B full-space cycle type (levels merged)
  C per-level fixed-point counts
  D total fixed-point count
  E per-level cycle types on the sign-quotient F_q^3 / K, K = Klein sign group
    {id, (-x,-y,z), (-x,y,-z), (x,-y,-z)}   (PSL-style coarsening)
For each: simultaneous p_min per D; distribution over the window.
"""
import os, sys
_R = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [_R, os.path.join(_R, "src"), os.path.dirname(os.path.abspath(__file__))]

import sys, json
from collections import Counter
from vieta import word_perm, level_ct, _letters
from markoff_rigidity import classes, canon, automorph, decomp

QLIST = [3, 5, 7, 11, 13, 17, 19, 23]

def gl_reps(T):
    D = T * T - 4
    cls = classes(D)
    def mcls(M):
        (a, b), (c, d) = M
        return canon((c, d - a, -b), D)
    reps = []; seen = set()
    for c0 in cls:
        if c0 in seen: continue
        M = automorph(c0, T, 1)
        (a, b), (c, d) = M
        orb = {mcls(M), mcls([[d, -b], [-c, a]]),
               mcls([[a, -b], [-c, d]]), mcls([[d, b], [c, a]])}
        seen |= orb
        reps.append(c0)
    return D, reps

def cycles_with_levels(perm, q):
    pts = _letters[('pts', q)]
    n = len(perm); seen = bytearray(n); out = []
    for s0 in range(n):
        if seen[s0]: continue
        c = 0; t = s0
        while True:
            seen[t] = 1; c += 1; t = perm[t]
            if t == s0: break
        x, y, z = pts[s0]
        out.append(((x*x + y*y + z*z - x*y*z) % q, c, s0))
    return out

def sign_orbit_perm(perm, q):
    """Induced permutation on F_q^3 / Klein sign group."""
    pts = _letters[('pts', q)]
    idx = {p: i for i, p in enumerate(pts)}
    def canon_pt(i):
        x, y, z = pts[i]
        cands = [(x, y, z), ((-x) % q, (-y) % q, z),
                 ((-x) % q, y, (-z) % q), (x, (-y) % q, (-z) % q)]
        return min(cands)
    orbit_rep = {}
    reps = []
    for i in range(len(pts)):
        cp = canon_pt(i)
        if cp not in orbit_rep:
            orbit_rep[cp] = len(reps); reps.append(cp)
    qperm = [None] * len(reps)
    for cp, oi in orbit_rep.items():
        j = perm[idx[cp]]
        qperm[oi] = orbit_rep[canon_pt(j)]
    return qperm, reps

def inv_all(word, q):
    perm = word_perm(word, q)
    cyc = cycles_with_levels(perm, q)
    # A: per-level cycle types
    lev = {}
    for k, c, _ in cyc: lev.setdefault(k, []).append(c)
    A = tuple(sorted((k, tuple(sorted(v))) for k, v in lev.items()))
    # B: full cycle type
    B = tuple(sorted(c for _, c, _ in cyc))
    # C: per-level fixed counts
    fc = Counter(k for k, c, _ in cyc if c == 1)
    C = tuple(sorted(fc.items()))
    # D: total fixed count
    Dv = sum(1 for _, c, _ in cyc if c == 1)
    # E: per-level cycle types on sign quotient
    qperm, reps = sign_orbit_perm(perm, q)
    n = len(qperm); seen = bytearray(n); levE = {}
    for s0 in range(n):
        if seen[s0]: continue
        c = 0; t = s0
        while True:
            seen[t] = 1; c += 1; t = qperm[t]
            if t == s0: break
        x, y, z = reps[s0]
        levE.setdefault((x*x + y*y + z*z - x*y*z) % q, []).append(c)
    E = tuple(sorted((k, tuple(sorted(v))) for k, v in levE.items()))
    return {'A': A, 'B': B, 'C': C, 'D': Dv, 'E': E}

def scan(t0, t1):
    dists = {k: Counter() for k in 'ABCDE'}
    per_T = {}
    for T in range(t0, t1 + 1):
        D, reps = gl_reps(T)
        n = len(reps)
        if n < 2: continue
        words = {c: decomp(automorph(c, T, 1)) for c in reps}
        cache = {}
        def get(c, q):
            if (c, q) not in cache: cache[(c, q)] = inv_all(words[c], q)
            return cache[(c, q)]
        row = {}
        for key in 'ABCDE':
            p = None
            for q in QLIST:
                vals = [get(c, q)[key] for c in reps]
                if len(set(vals)) == n:
                    p = q; break
            dists[key][p] += 1
            row[key] = p
        per_T[T] = row
        print(T, row, flush=True)
    for key in 'ABCDE':
        print(key, dict(dists[key]))
    print('frozen target: {5: 20, 7: 18, 11: 7, 13: 1}')
    json.dump({'per_T': per_T, 'dists': {k: dict(v) for k, v in dists.items()}},
              open(os.path.join(_R,'data',f'window_audit2.json'),'w'), default=str)

if __name__ == '__main__':
    scan(240, 285)
