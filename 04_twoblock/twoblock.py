"""The two-block stratum: enumeration of the classes L^u R^v with uv = T - 2
and their stratified vocabulary (Table 4).
"""
import os, sys
_R = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [_R, os.path.join(_R, "src"), os.path.dirname(os.path.abspath(__file__))]

import numpy as np, sys, json
from math import lcm, gcd
from collections import Counter, defaultdict
from scan_parallel import _tables
from markoff_rigidity import decomp

def perm_of_word(word, q):
    Tb=_tables(q); p=np.arange(q**3,dtype=np.int32)
    for g in word: p=Tb[g][p]
    return p

def compose(p,r): return p[r]          # apply r then p?  fixed convention below

def power(p,n):
    res=np.arange(len(p),dtype=np.int32)
    b=p.copy()
    while n:
        if n&1: res=b[res]
        b=b[b]; n>>=1
    return res

def ctype(p):
    n=len(p); seen=bytearray(n); ct=Counter(); pl=p.tolist()
    for s in range(n):
        if seen[s]: continue
        c=0;t=s
        while not seen[t]: seen[t]=1;c+=1;t=pl[t]
        ct[c]+=1
    return tuple(sorted(ct.items()))

def L_of(q): return lcm(q, q-1, q+1)

def table(q):
    Lq=L_of(q)
    wR=decomp([[1,1],[0,1]]); wL=decomp([[1,0],[1,1]])
    PR=perm_of_word(wR,q); PL=perm_of_word(wL,q)
    # check orders
    def order(p):
        k=1; c=p.copy(); idn=np.arange(len(p),dtype=np.int32)
        while not np.array_equal(c,idn): c=p[c]; k+=1
        return k
    oR,oL=order(PR),order(PL)
    PRp=[power(PR,i) for i in range(Lq)]
    PLp=[power(PL,i) for i in range(Lq)]
    fib=defaultdict(list)
    for u in range(Lq):
        for v in range(Lq):
            fib[ctype(PLp[u][PRp[v]])].append((u,v))
    return Lq,oR,oL,fib

if __name__=='__main__':
    q=int(sys.argv[1])
    Lq,oR,oL,fib=table(q)
    sizes=Counter(len(v) for v in fib.values())
    print(f"q={q} L(q)={Lq}  ord sigma_R={oR} ord sigma_L={oL}")
    print(f"distinct cycle types among L^u R^v: {len(fib)} out of {Lq*Lq} residue pairs")
    print("fibre size distribution:", dict(sorted(sizes.items())))
    json.dump({str(k):v for k,v in fib.items()}, open(f'fib_{q}.json','w'))
