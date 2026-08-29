"""Detection at ramified and at inert primes (Remark rem:ramdetect and the
failures at T = 102 and T = 258).
"""
import os, sys
_R = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [_R, os.path.join(_R, "src"), os.path.dirname(os.path.abspath(__file__))]

import sys
from collections import Counter
from scan_parallel import class_reps, gl_reps, cycle_types_35, cycle_type
from markoff_rigidity import automorph, decomp
from sympy import primerange
QS=[3,5,7,11,13]
lo,hi=int(sys.argv[1]),int(sys.argv[2])
diff=Counter(); dist_a=Counter(); dist_b=Counter(); ex=[]
for T in range(lo,hi+1):
    D=T*T-4
    reps=gl_reps(T,D,class_reps(D)); n=len(reps)
    if n<2: continue
    W=[decomp(automorph(c,T,1)) for c in reps]
    inv={}
    for q in QS:
        inv[q]=[cycle_types_35(w)[0 if q==3 else 1] if q in (3,5) else cycle_type(w,q) for w in W]
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    # (b) max over pairs of least detecting prime
    b=0
    for i,j in pairs:
        for q in QS:
            if inv[q][i]!=inv[q][j]: b=max(b,q); break
        else: b=99
    # (a) least single prime separating ALL pairs
    a=99
    for q in QS:
        if len(set(inv[q]))==n: a=q; break
    dist_a[a]+=1; dist_b[b]+=1
    if a!=b:
        diff[(b,a)]+=1
        if len(ex)<8: ex.append((T,n,b,a))
print(f"T in [{lo},{hi}]")
print(" (b) max-over-pairs  :", dict(sorted(dist_b.items())))
print(" (a) single prime    :", dict(sorted(dist_a.items())))
print(" disagreements (b,a):", dict(sorted(diff.items())), " total", sum(diff.values()))
print(" examples (T,nGL,b,a):", ex)
