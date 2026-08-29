"""The table of T_p: the least trace at which the two-block construction is
blind to every prime up to p (Corollaries cor:Mp and cor:order).
"""
import os, sys
_R = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [_R, os.path.join(_R, "src"), os.path.dirname(os.path.abspath(__file__))]

import numpy as np
from math import lcm, log10
from collections import Counter
from sympy import primerange
from scan_parallel import _tables
from markoff_rigidity import decomp
from twoblock import perm_of_word, power, ctype
def L_of(q): return lcm(q,q-1,q+1)
CACHE={}
def PLR(q):
    if q not in CACHE:
        CACHE[q]=(perm_of_word(decomp([[1,0],[1,1]]),q), perm_of_word(decomp([[1,1],[0,1]]),q))
    return CACHE[q]
def ct(u,v,q):
    PL,PR=PLR(q); Lq=L_of(q)
    return ctype(power(PL,u%Lq)[power(PR,v%Lq)])
print(f"{'p':>3} {'Lambda_p':>22} {'T_p':>24} {'log10 T':>8} {'witness':>8}")
Lam=1
for p in primerange(2,60):
    Lam=lcm(Lam,L_of(p))
    if p<3: continue
    B=2*(Lam+1); T=B+2
    u1,v1=1,B; u2,v2=Lam+1,2
    w=None
    for q in primerange(3,200):
        if ct(u1,v1,q)!=ct(u2,v2,q): w=q; break
    # sanity: blind at all q<=p
    ok=all(ct(u1,v1,q)==ct(u2,v2,q) for q in primerange(3,p+1))
    print(f"{p:>3} {Lam:>22} {T:>24} {log10(T):>8.1f} {w:>8}   blind<=p:{ok}")
