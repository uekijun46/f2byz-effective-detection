#!/usr/bin/env python3
"""bigscan.py -- exhaustive record hunt for simultaneous blindness.

Uses fastscan.class_reps (coset construction, ~90x faster) and
fastscan.cycle_type_np (numpy Vieta permutation).

Records sought (draft §6 table):
  {3,5,7}-blind      least trace : claimed 1350   (minimal, exhaustive below)
  {3,5,7,11}-blind   least trace : claimed 18476  (least KNOWN)
  {3,5,7,11,13}-blind least trace: claimed 1021022 (least KNOWN)
"blind at S" = some genuine pair of that discriminant is separated by no q in S
             = p_min(D) > max(S).

Protocol: all primes (no ramification restriction), full-space cycle type,
GL folding by <inverse, det -1 conjugation>.  Progressive refinement:
classes are bucketed by the q=3 invariant, then only colliding buckets are
refined at 5, 7, 11, 13 -- so large primes are touched only when needed.
"""
import os, sys
_R = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [_R, os.path.join(_R, "src"), os.path.dirname(os.path.abspath(__file__))]

import sys, json, time
from fastscan import class_reps, cycle_type_np
from markoff_rigidity import canon, automorph, decomp

QS = [3, 5, 7, 11, 13, 17, 19, 23]

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

def p_min(T):
    """Least p such that {3,...,p} jointly separate all GL-classes; None if >23."""
    D = T * T - 4
    cls = class_reps(D)
    reps = gl_reps(T, D, cls)
    n = len(reps)
    if n < 2: return 0, n
    words = [decomp(automorph(c, T, 1)) for c in reps]
    # progressive refinement: buckets of indices sharing all invariants so far
    buckets = [list(range(n))]
    for q in QS:
        newb = []
        for bk in buckets:
            if len(bk) == 1: continue
            sub = {}
            for i in bk:
                sub.setdefault(cycle_type_np(words[i], q), []).append(i)
            for v in sub.values():
                if len(v) > 1: newb.append(v)
        buckets = newb
        if not buckets:
            return q, n
    return None, n          # not separated by any q <= 23

def scan(t0, t1, log_every=250):
    rec = {}          # threshold -> least T with p_min > threshold
    dist = {}
    t_start = time.time()
    for T in range(t0, t1 + 1):
        p, n = p_min(T)
        if n < 2: continue
        key = p if p else '>23'
        dist[key] = dist.get(key, 0) + 1
        for thr, name in ((7, '{3,5,7}'), (11, '{3,5,7,11}'), (13, '{3,5,7,11,13}')):
            if (p is None or p > thr) and name not in rec:
                rec[name] = T
                print(f"  *** first {name}-blind trace: T={T} "
                      f"(p_min={p}, nGL={n})", flush=True)
        if T % log_every == 0:
            el = time.time() - t_start
            print(f"  ...T={T}  ({el:.0f}s, {el/(T-t0+1)*1000:.1f} ms/T)  "
                  f"dist={dict(sorted(dist.items(), key=lambda kv: str(kv[0])))}",
                  flush=True)
    return rec, dist

if __name__ == '__main__':
    t0, t1 = int(sys.argv[1]), int(sys.argv[2])
    rec, dist = scan(t0, t1)
    print("\n=== scan", t0, "..", t1, "===")
    print("distribution:", dict(sorted(dist.items(), key=lambda kv: str(kv[0]))))
    print("records found:", rec)
    json.dump({'range': [t0, t1], 'dist': {str(k): v for k, v in dist.items()},
               'records': rec}, open(os.path.join(_R,'data',f'big_{t0}_{t1}.json'),'w'))
