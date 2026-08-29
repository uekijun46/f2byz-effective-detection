#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify_blindness_law.py -- numerical verification of the Blindness law
(Theorem D / thm:family in "Effective detection of F_2-by-Z groups").

PROVISIONAL public version.  Stdlib only; reuses the conventions of
src/markoff_rigidity.py (letters S, Si, T, Ti, J and their Vieta action).

WHAT IS CHECKED
===============
Fix a prime q, a prime p0, and a parameter M.  Let
L(q) = q * lcm(q-1, q+1), a multiple of exp SL_2(F_q).  Put b = M*p0 and
T_trace = b + 2.  With R = T = [[1,1],[0,1]] and L = [[1,0],[1,1]] the two
matrices
    A1 = L^{-1}  R^{-b}  = [[1, -b ], [-1,  b+1]],
    A2 = L^{-p0} R^{-M}  = [[1, -M ], [-p0, b+1]]
share the characteristic polynomial x^2 - T_trace x + 1.  We verify:

 (1) the word identities above, exactly over Z, via word_to_mat;
 (2) Theorem D: if L(q) | M, the Vieta permutations of A1 and A2 on
     F_q^3 have identical cycle types on EVERY level set of
     kappa(x,y,z) = x^2 + y^2 + z^2 - x*y*z   (blindness at q);
 (3) a control run with L(q) not dividing M', reported for information
     (experimentally such pairs are detected; see Remark on the converse
     in the paper -- this direction is not asserted by Theorem D).

Words in the letters of markoff_rigidity.py:
    L^{-1} = S T Si   (verified),   R^{-1} = Ti.
so  A1 = [S,T,Si] + [Ti]*b,   A2 = [S] + [T]*p0 + [Si] + [Ti]*M.

USAGE
=====
    python3 verify_blindness_law.py                 # q=5, p0=7, M=L(5)=60
    python3 verify_blindness_law.py --q 7 --p0 11
    python3 verify_blindness_law.py --q 5 --M 120
Runtime O(q^3) per letter-permutation; q <= 13 takes seconds.
"""

import argparse, os, sys
from math import gcd
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "src"))
from markoff_rigidity import GENMAT, word_to_mat, _vstep  # noqa: E402

def lcm(a, b): return a // gcd(a, b) * b

# ---------- exact word identities ----------
def check_word_identities(M, p0):
    b = M * p0
    A1 = word_to_mat(["S", "T", "Si"] + ["Ti"] * b)
    A2 = word_to_mat(["S"] + ["T"] * p0 + ["Si"] + ["Ti"] * M)
    assert A1 == [[1, -b], [-1, b + 1]],  "word identity for A1 fails"
    assert A2 == [[1, -M], [-p0, b + 1]], "word identity for A2 fails"
    return b

# ---------- letter permutations on F_q^3, with exponent reduction ----------
def letter_perm(g, q):
    return {(x, y, z): _vstep(g, x, y, z, q)
            for x in range(q) for y in range(q) for z in range(q)}

def perm_order(p):
    idp = {k: k for k in p}
    cur, n = dict(p), 1
    while cur != idp:
        cur = {k: p[cur[k]] for k in cur}
        n += 1
    return n

def word_perm(word_runs, q):
    """word_runs: list of (letter, count); rightmost run acts FIRST
    (matching word_to_mat's left-to-right matrix product acting on column
    vectors is irrelevant here -- we only need one consistent convention
    for both words, which suffices for cycle-type comparison)."""
    cache = {}
    perm = None
    for g, k in reversed(word_runs):
        if g not in cache:
            cache[g] = (letter_perm(g, q), None)
        base, ordg = cache[g]
        if ordg is None:
            ordg = perm_order(base)
            cache[g] = (base, ordg)
        k %= ordg
        step = {t: t for t in base}
        for _ in range(k):
            step = {t: base[step[t]] for t in step}
        perm = step if perm is None else {t: step[perm[t]] for t in perm}
    return perm

def kappa(t, q):
    x, y, z = t
    return (x * x + y * y + z * z - x * y * z) % q

def cycle_types_by_level(perm, q):
    seen, levels = set(), {}
    for start in perm:
        if start in seen:
            continue
        n, t = 0, start
        while True:
            seen.add(t); n += 1; t = perm[t]
            if t == start:
                break
        levels.setdefault(kappa(start, q), []).append(n)
    return {lv: tuple(sorted(c)) for lv, c in levels.items()}

# ---------- runs ----------
def run(q, p0, M, expect_blind):
    b = check_word_identities(M, p0)
    pA = word_perm([("S", 1), ("T", 1), ("Si", 1), ("Ti", b)], q)
    pB = word_perm([("S", 1), ("T", p0), ("Si", 1), ("Ti", M)], q)
    same = cycle_types_by_level(pA, q) == cycle_types_by_level(pB, q)
    verdict = ("BLIND (equal cycle types on every level)" if same
               else "DETECTED (cycle types differ on some level)")
    print(f"  q={q}, p0={p0}, M={M}, b={b}, trace T={b+2}: {verdict}")
    if expect_blind is not None:
        assert same == expect_blind, "verdict contradicts expectation"
    return same

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--q",  type=int, default=5)
    ap.add_argument("--p0", type=int, default=7)
    ap.add_argument("--M",  type=int, default=None)
    a = ap.parse_args()
    Lq = a.q * lcm(a.q - 1, a.q + 1)
    M = a.M if a.M is not None else Lq
    print(f"[verify_blindness_law] q={a.q}, L(q)={Lq}")
    print("(1)+(2) Theorem D direction (L(q) | M => blindness is ASSERTED):")
    run(a.q, a.p0, M, expect_blind=True if M % Lq == 0 else None)
    print("(3) control, L(q) not dividing M' (informational; converse is")
    print("    experimental, see the paper's remark):")
    for Mc in (Lq // 2, Lq // 3, Lq // 4, Lq - 2):
        if Mc >= 2 and Mc % Lq != 0 and gcd(Mc, a.p0) == 1:
            run(a.q, a.p0, Mc, expect_blind=None)
            break
    print("all asserted checks passed.")

if __name__ == "__main__":
    main()
