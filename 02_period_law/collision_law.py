#!/usr/bin/env python3
"""collision_law.py -- reproduces the statistics of Section 7(b) on the
structure of colliding pairs.

For every trace T in the given range that is blind at {3,5,7} (i.e. some pair
of GL-classes survives all three primes), this reports
  * the sizes of the surviving buckets  (how many classes remain unseparated),
  * the class-group difference of each colliding pair, and whether it is
    ambiguous (2-torsion),
  * whether that ambiguous class is the one cut out by a factorization of T-2.

Reproduces: "all but two of the 46 blind traces T <= 4000 carry exactly one
undetected pair", "the two classes differ by an ambiguous class in 42 of the 52
pairs", "in 25 of those the ambiguous class comes from a factorization of T-2".

Usage:  python3 collision_law.py 3 4000
"""
import os, sys
_R = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [_R, os.path.join(_R, "src"), os.path.dirname(os.path.abspath(__file__))]

import sys
from collections import Counter

from scan_parallel import (class_reps, gl_reps, cycle_types_35, cycle_type,
                           _compose, _reduce)
from markoff_rigidity import canon, automorph, decomp


def inv_class(f, D):
    a, b, c = f
    return canon(_reduce((a, -b, c), D), D)


def mul(f, g, D):
    return canon(_compose(f, g, D), D)


def identity(D):
    b0 = D % 2
    return canon(_reduce((1, b0, (b0 * b0 - D) // 4), D), D)


def surviving_buckets(T):
    """Classes of D = T^2-4 left unseparated by q = 3, 5, 7."""
    D = T * T - 4
    reps = gl_reps(T, D, class_reps(D))
    n = len(reps)
    if n < 2:
        return D, reps, []
    W = [decomp(automorph(c, T, 1)) for c in reps]
    buckets = [list(range(n))]
    for q in (3, 5, 7):
        nb = []
        for bk in buckets:
            if len(bk) == 1:
                continue
            sub = {}
            for i in bk:
                v = (cycle_types_35(W[i])[0 if q == 3 else 1]
                     if q in (3, 5) else cycle_type(W[i], q))
                sub.setdefault(v, []).append(i)
            nb += [v for v in sub.values() if len(v) > 1]
        buckets = nb
        if not buckets:
            break
    return D, reps, [[reps[i] for i in b] for b in buckets]


def main(lo, hi):
    stats = Counter()
    sizes = Counter()
    for T in range(lo, hi + 1):
        D, reps, buckets = surviving_buckets(T)
        if not buckets:
            continue
        stats['blind traces'] += 1
        e = identity(D)
        for b in buckets:
            sizes[len(b)] += 1
            for i in range(len(b)):
                for j in range(i + 1, len(b)):
                    c1, c2 = tuple(b[i]), tuple(b[j])
                    stats['pairs'] += 1
                    # the GL-folding leaves the difference defined only up to
                    # inversion, so both products are candidates
                    cands = [mul(inv_class(c1, D), c2, D), mul(c1, c2, D)]
                    amb = [d for d in cands if mul(d, d, D) == e]
                    if not amb:
                        stats['difference not ambiguous'] += 1
                        continue
                    stats['difference ambiguous'] += 1
                    a, bb, c = amb[0]
                    if abs(a * c) == T - 2 or bb == T - 2:
                        stats['   ... cut out by a factorization of T-2'] += 1
        print(f"T={T:>6} nGL={len(reps):>4} surviving buckets "
              f"{[len(b) for b in buckets]}", flush=True)
    print()
    print("bucket sizes:", dict(sizes))
    for k, v in stats.items():
        print(f"{k}: {v}")


if __name__ == '__main__':
    main(int(sys.argv[1]), int(sys.argv[2]))
