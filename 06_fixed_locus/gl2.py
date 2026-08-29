"""sigma_M(q) for an arbitrary M in GL_2(Z), by factoring M into the Nielsen
generators R, L, E = diag(1,-1)."""
import os, sys
_R = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [_R, os.path.join(_R, "src"), os.path.dirname(os.path.abspath(__file__))]

import numpy as np
from cyc import build, Lq

def factor(M):
    """Return a list of ('R',k)/('L',k)/('E',1) whose product is +-M."""
    a,b,c,d = M
    det = a*d-b*c
    assert det in (1,-1)
    tail = []
    if det == -1:                    # M = M' E with M' in SL_2 Z
        a,b,c,d = a,-b,c,-d
        tail = [('E',1)]
    ops = []                          # left multipliers U_1,...,U_n
    while c != 0:
        if a == 0:                    # then |c| = 1; make a nonzero
            a, b = a+c, b+d
            ops.append(('R', 1))
        else:
            k = -(a // c)             # row1 += k row2
            if k:
                a, b = a+k*c, b+k*d
                ops.append(('R', k))
            if c == 0:
                break
            if a == 0:
                continue
        k2 = -(c // a)                # row2 += k2 row1
        if k2:
            c, d = c+k2*a, d+k2*b
            ops.append(('L', k2))
    assert c == 0 and a*d == 1, (a,b,c,d)
    rem = [('R', a*b)]
    # M = U_1^{-1} ... U_n^{-1} * rem  (up to sign)
    word = [(t,-k) for (t,k) in ops] + rem + tail
    return word

class Sig:
    def __init__(self, q):
        self.q=q; self.n=q**3; self.Lq=Lq(q)
        R,L = build(q)
        self.R=R; self.L=L
        idx=np.arange(self.n); x=idx//(q*q); y=(idx//q)%q; z=idx%q
        self.E = (x*q*q + y*q + (x*y-z) % q).astype(np.int32)
        self.cacheR={}; self.cacheL={}
    def pw(self, P, k, cache):
        k %= self.Lq
        if k not in cache:
            arr=np.arange(self.n,dtype=np.int32)
            for _ in range(k): arr=P[arr]
            cache[k]=arr
        return cache[k]
    def perm_word(self, word):
        arr=np.arange(self.n,dtype=np.int32)
        for t,k in reversed(word):
            if t=='R': arr=self.pw(self.R,k,self.cacheR)[arr]
            elif t=='L': arr=self.pw(self.L,k,self.cacheL)[arr]
            else: arr=self.E[arr]
        return arr
    def ctype(self, M):
        p=self.perm_word(factor(M))
        n=len(p); seen=np.zeros(n,dtype=bool); out=[]
        for i in range(n):
            if not seen[i]:
                c=0;j=i
                while not seen[j]: seen[j]=True; j=p[j]; c+=1
                out.append(c)
        return tuple(sorted(out,reverse=True))

def check_factor(M):
    import numpy as np
    def mm(X,Y): return (X[0]*Y[0]+X[1]*Y[2], X[0]*Y[1]+X[1]*Y[3], X[2]*Y[0]+X[3]*Y[2], X[2]*Y[1]+X[3]*Y[3])
    P=(1,0,0,1)
    for t,k in factor(M):
        if t=='R': G=(1,k,0,1)
        elif t=='L': G=(1,0,k,1)
        else: G=(1,0,0,-1)
        P=mm(P,G)
    return P, M

if __name__=="__main__":
    for M in [(7,6,1,1),(7,3,2,1),(3,10,1,3),(3,5,2,3),(19,60,6,19),(19,30,12,19)]:
        P,_=check_factor(M)
        print(M, "->", P, "ok" if P==M or P==tuple(-x for x in M) else "MISMATCH")
