"""Test of a SHORTER proof of Prop 7.4, via F_2[sigma]-module algebra.
P = pointwise twists p_a (a: X -> Sigma, Sigma-invariant), sigma a := rho(a o U^{-1}).
  (i)   (JU)^k = p_{N_k c} U^k,  N_k = 1+sigma+...+sigma^{k-1}
  (ii)  p_a (p_c U^k) p_a = p_{c + (1+sigma^k)a} U^k
  (iii) k odd  =>  (N_k+1)c = (1+sigma) w,  w = sum_{i<(k-1)/2} sigma^{2i+1} c
  (iv)  gcd(k, ord U)=1 => 1+sigma = (1+sigma^k) Q,  Q = sum_{i<k'} sigma^{ki}, k'=k^{-1} mod n
  => a := Q w satisfies (1+sigma^k)a = (N_k+1)c, i.e. p_a conjugates JU^k to (JU)^k.
No use of rho^2=id, no parity case split, no descent to X/Sigma.
"""
import os, sys
_R = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [_R, os.path.join(_R, "src"), os.path.dirname(os.path.abspath(__file__))]
from signtwist_core import *
from math import gcd

def test(q,M,p0):
    L=Lq(q); U=word_perm(['S','T','Si'],q); Ui={v:k for k,v in U.items()}
    Rm=word_perm(['Ti'],q); J=perm_pow(Rm,(-M)%L)
    n=perm_order(U)
    SG=[(0,0),(0,1),(1,0),(1,1)]
    def s(d,t): return ((-1)**d[0]*t[0]%q,(-1)**d[1]*t[1]%q,(-1)**((d[0]+d[1])%2)*t[2]%q)
    rho={}
    for d in SG:
        for e in SG:
            if all(U[s(d,Ui[t])]==s(e,t) for t in U): rho[d]=e; break
    addv=lambda a,b:((a[0]+b[0])%2,(a[1]+b[1])%2)
    # coefficient c of J
    c={t:[d for d in SG if s(d,t)==J[t]][0] for t in U}
    def sig(a):                      # sigma a = rho(a o U^{-1})
        return {t: rho[a[Ui[t]]] for t in a}
    def addf(a,b): return {t: addv(a[t],b[t]) for t in a}
    def twist(a): return {t: s(a[t],t) for t in a}
    # (iii) w
    w={t:(0,0) for t in U}; cur=sig(c)                      # sigma^1 c
    for i in range((p0-1)//2):
        w=addf(w,cur); cur=sig(sig(cur))                    # sigma^{2i+1} c
    # (iv) a = Q w
    kp=pow(p0,-1,n)
    a={t:(0,0) for t in U}; cur=dict(w)
    for i in range(kp):
        a=addf(a,cur)
        for _ in range(p0): cur=sig(cur)                    # sigma^{k} step
    # check (1+sigma^k) a == (N_k+1) c
    sk=dict(a)
    for _ in range(p0): sk=sig(sk)
    lhs=addf(a,sk)
    Nc={t:(0,0) for t in U}; cur=dict(c)
    for j in range(p0):
        Nc=addf(Nc,cur); cur=sig(cur)
    rhs=addf(Nc,c)
    ok_alg = lhs==rhs
    # check the conjugation on permutations
    Uk=perm_pow(U,p0)
    W1=compose(U,J)          # J o U  (apply U then J)
    Wk=compose(Uk,J)
    pa=twist(a)
    lhsP=perm_pow(W1,p0)
    rhsP=compose(compose(pa,Wk),pa)
    ok_perm = lhsP==rhsP
    print(f'q={q} M={M} p0={p0}: ord U={n}  (1+s^k)a=(N_k+1)c: {ok_alg}   p_a(JU^k)p_a=(JU)^k: {ok_perm}'
          f'   cyc equal: {cyc_by_level(W1,q)==cyc_by_level(Wk,q)}')
    return ok_alg and ok_perm

allok=True
for q,Ms in ((3,[6,12,18]),(5,[30,60,90]),(7,[84,168]),(11,[330]),(13,[546])):
    for M in Ms:
        for p0 in (7,11,13):
            if gcd(p0,M)!=1: continue
            allok &= test(q,M,p0)
print('ALL OK' if allok else 'FAILURE')
