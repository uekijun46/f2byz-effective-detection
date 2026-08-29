#!/usr/bin/env python3
"""hom_check.py -- independent verification, bypassing the Markoff model.

Instead of acting on the character variety F_q^3, act on the HONEST object:
   Hom(F_2, G) / Inn(G),   G = SL_2(F_q),
by the actual automorphism of F_2 (word in S,Si,T,Ti,J applied to a pair of
matrices), and compute the cycle type.  Cycle types of this action are what
counting Hom(Gamma, G wr Z/k) reads off.  No trace coordinates involved.

Question asked here: for a discriminant where the SMALLEST separating prime is
RAMIFIED, do the two classes really get different cycle types on Hom/Inn?
"""
import os, sys
_R = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [_R, os.path.join(_R, "src"), os.path.dirname(os.path.abspath(__file__))]

import sys
from collections import Counter
from markoff_rigidity import classes, canon, automorph, decomp, kronecker_split
from coarser_invariants import gl_reps, inv_all

# ---------------- SL_2(F_q) -------------------------------------------
def sl2(q):
    els = []
    for a in range(q):
        for b in range(q):
            for c in range(q):
                for d in range(q):
                    if (a*d - b*c) % q == 1:
                        els.append((a, b, c, d))
    return els

def mul(P, Q, q):
    a, b, c, d = P; e, f, g, h = Q
    return ((a*e + b*g) % q, (a*f + b*h) % q, (c*e + d*g) % q, (c*f + d*h) % q)

def inv(P, q):
    a, b, c, d = P
    return (d % q, (-b) % q, (-c) % q, a % q)   # det = 1

# ---------------- action of Aut(F_2) on pairs --------------------------
# conventions of markoff_rigidity.py:
#   S : a -> b,     b -> a^{-1}
#   Si: a -> b^{-1},b -> a
#   T : a -> a,     b -> ab
#   Ti: a -> a,     b -> a^{-1}b
#   J : a -> a,     b -> b^{-1}
def apply_letter(g, x, y, q):
    if g == 'S':  return (y, inv(x, q))
    if g == 'Si': return (inv(y, q), x)
    if g == 'T':  return (x, mul(x, y, q))
    if g == 'Ti': return (x, mul(inv(x, q), y, q))
    if g == 'J':  return (x, inv(y, q))
    raise ValueError(g)

def apply_word(w, x, y, q):
    for g in reversed(w):        # same order convention as word_perm
        x, y = apply_letter(g, x, y, q)
    return x, y

# ---------------- Hom(F_2,G)/Inn ---------------------------------------
def inn_orbits(G, q):
    idx = {P: i for i, P in enumerate(G)}
    pairs = [(x, y) for x in G for y in G]
    pidx = {p: i for i, p in enumerate(pairs)}
    rep_of = [-1] * len(pairs)
    reps = []
    for i, (x, y) in enumerate(pairs):
        if rep_of[i] >= 0: continue
        r = len(reps); reps.append((x, y))
        for g in G:
            gi = inv(g, q)
            xx = mul(mul(g, x, q), gi, q); yy = mul(mul(g, y, q), gi, q)
            rep_of[pidx[(xx, yy)]] = r
    return pairs, pidx, rep_of, reps

def cycle_type_on_hom(word, G, q, pairs, pidx, rep_of, reps):
    # induced permutation on Inn-orbits
    perm = [None] * len(reps)
    for r, (x, y) in enumerate(reps):
        xx, yy = apply_word(word, x, y, q)
        perm[r] = rep_of[pidx[(xx, yy)]]
    n = len(reps); seen = bytearray(n); cyc = []
    for s in range(n):
        if seen[s]: continue
        c = 0; t = s
        while True:
            seen[t] = 1; c += 1; t = perm[t]
            if t == s: break
        cyc.append(c)
    return tuple(sorted(cyc))

def run(T, q):
    D, reps_cls = gl_reps(T)
    G = sl2(q)
    pairs, pidx, rep_of, reps = inn_orbits(G, q)
    print(f"T={T} D={D} q={q} ({kronecker_split(D,q)}), |SL2(F_q)|={len(G)}, "
          f"|Hom(F_2,G)/Inn|={len(reps)}, #GL-classes={len(reps_cls)}")
    cts = []
    for c in reps_cls:
        w = decomp(automorph(c, T, 1))
        ct = cycle_type_on_hom(w, G, q, pairs, pidx, rep_of, reps)
        cts.append(ct)
        cc = Counter(ct)
        print(f"   class {c}: cycle type multiset {dict(sorted(cc.items()))}")
    print(f"   distinct cycle types: {len(set(cts))} of {len(cts)} classes "
          f"-> {'SEPARATES' if len(set(cts))==len(cts) else 'does NOT separate'}")
    # compare with character-variety computations
    cv = [inv_all(decomp(automorph(c, T, 1)), q) for c in reps_cls]
    print(f"   [char.variety] full-space CT distinct: {len({v['B'] for v in cv})}, "
          f"per-level CT distinct: {len({v['A'] for v in cv})}")
    return cts

if __name__ == '__main__':
    # T=243: ram = {5,7}; smallest separating prime under 'all primes' is 5 (RAMIFIED)
    run(243, 5)
    print()
    run(243, 11)   # smallest unramified separating prime, for comparison
