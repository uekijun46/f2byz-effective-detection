"""core.py -- reconstructed from psl_threshold_scan.py (which is self-contained).
Vieta action on F_q^3, letters S,T,Ti; L^{-1}=S T Si, R^{-1}=Ti.
compose(f,g): apply f first, then g.  word_perm applies letters left to right.
"""
import os, sys
_R = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [_R, os.path.join(_R, "src"), os.path.dirname(os.path.abspath(__file__))]

from math import gcd
from functools import reduce
from collections import Counter

def vstep(g,x,y,z,q):
    if g=='T':  return (x,z,(x*z-y)%q)
    if g=='Ti': return (x,(x*y-z)%q,y)
    if g in ('S','Si'): return (y,x,(x*y-z)%q)
    raise ValueError(g)

def pts(q): return [(x,y,z) for x in range(q) for y in range(q) for z in range(q)]

def letter_perm(g,q): return {p: vstep(g,*p,q) for p in pts(q)}

def compose(f,g):   # apply f first, then g   (point -> g(f(pt)))
    return {p: g[f[p]] for p in f}

def word_perm(word,q,cache={}):
    """word = list of letters, applied LEFT to RIGHT to a point."""
    perm = {p:p for p in pts(q)}
    for g in word:
        key=(g,q)
        if key not in cache: cache[key]=letter_perm(g,q)
        perm = compose(perm, cache[key])
    return perm

def perm_pow(base,k):
    res = {p:p for p in base}; b=dict(base)
    while k:
        if k&1: res = compose(res,b)
        b = compose(b,b); k>>=1
    return res

def perm_order(p):
    idp={k:k for k in p}; cur=dict(p); n=1
    while cur!=idp:
        cur=compose(cur,p); n+=1
    return n

def kappa(t,q):
    x,y,z=t; return (x*x+y*y+z*z-x*y*z)%q

def cyc_by_level(perm,q):
    seen=set(); lv={}
    for s in perm:
        if s in seen: continue
        n=0;t=s
        while True:
            seen.add(t);n+=1;t=perm[t]
            if t==s: break
        lv.setdefault(kappa(s,q),[]).append(n)
    return {k:tuple(sorted(v)) for k,v in lv.items()}

def cyc_total(perm):
    seen=set(); c=Counter()
    for s in perm:
        if s in seen: continue
        n=0;t=s
        while True:
            seen.add(t);n+=1;t=perm[t]
            if t==s: break
        c[n]+=1
    return tuple(sorted(c.items()))

def Lq(q): return q*(q*q-1)//gcd(q-1,q+1)


from math import gcd
def E(q): return Lq(q)//2
def signtwist(D,q):
    return all(any(((a*p[0])%q,(b*p[1])%q,(a*b*p[2])%q)==D[p] for a in(1,-1) for b in(1,-1)) for p in D)
