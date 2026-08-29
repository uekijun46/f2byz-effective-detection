#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""psl_threshold_scan.py -- PSL threshold experiments for Theorem D (2026-08-14).
Stdlib only.  Conventions match src/markoff_rigidity.py: letters S,T,Ti and
the Vieta action on F_q^3; L^{-1}=S T Si, R^{-1}=Ti.
A1 = L^{-1} R^{-b} (b = M p0),  A2 = L^{-p0} R^{-M}.
Outputs: for each q and each M | L(q), which p0 are blind."""
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
primes=[5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,73]
for q in (3,5,7):
    L=Lq(q); e=E(q); U=word_perm(['S','T','Si'],q)
    Rm=word_perm(['Ti'],q); oR=perm_order(Rm)
    print(f'=== q={q}  L={L} expPSL={e} ord(sigma_R)={oR} ord(sigma_L)={perm_order(U)}')
    for M in [d for d in range(2,L+1) if L%d==0]:
        Bl=[];Dt=[]
        JM=perm_pow(Rm, (-M)%oR)
        for p0 in primes:
            if gcd(p0,M)!=1: continue
            Jb=perm_pow(Rm, (-M*p0)%oR)
            A1=compose(U,Jb); A2=compose(perm_pow(U,p0),JM)
            (Bl if cyc_by_level(A1,q)==cyc_by_level(A2,q) else Dt).append(p0)
        g=0
        for p0 in Bl: g=gcd(g,p0-1)
        print(f' M={M:4d} expPSL|M?{"Y" if M%e==0 else "n"} signtwist(R^-M)?{"Y" if signtwist(JM,q) else "n"} blind@{len(Bl)}/{len(Bl)+len(Dt)} gcd(p0-1)={g} blind={Bl}')
