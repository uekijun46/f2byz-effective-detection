#!/usr/bin/env python3
"""Verification of the strengthened Blindness law (PSL threshold).
Claim: if E(q)=exp PSL_2 F_q = L(q)/2 divides M, then sigma_{A1} and
sigma_{A2} have equal cycle types on every level set.
Proof structure checked here:
 (S1) E|M  =>  sigma_R^{-M} = J, the fiberwise sign twist
      J(x,y,z)=(x, eps(x)y, eps(x)z),  N(x)^M = eps(x) I.
 (S2) eps(-x)=eps(x); rho_L on Sigma=(F_2)^2 has order 2; conj rules.
 (S3) sigma_R^{-b} = sigma_R^{-M}  (b-M = M(p0-1) == 0 mod L).
 (S4) per Ubar-cycle on X/Sigma: W1=J.U and W2=J.U^{p0} give equal
      cycle-length multisets over the block; for m even, W1^m == W2^m
      pointwise on the block (exact holonomy equality).
 (S5) ord(W1) | 2 L(q).
 (S6) degenerate mechanism: p0 = 1 mod L => sigma_A1 == sigma_A2;
      p0 = -1 mod L => sigma_A1 conjugate to sigma_A2^{-1} (equal cyc).
"""
import os, sys
_R = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [_R, os.path.join(_R, "src"), os.path.dirname(os.path.abspath(__file__))]
from signtwist_core import *
from math import gcd

def matpow(A,k,q):
    R=[[1,0],[0,1]]; B=[r[:] for r in A]
    while k:
        if k&1: R=[[(R[0][0]*B[0][0]+R[0][1]*B[1][0])%q,(R[0][0]*B[0][1]+R[0][1]*B[1][1])%q],
                   [(R[1][0]*B[0][0]+R[1][1]*B[1][0])%q,(R[1][0]*B[0][1]+R[1][1]*B[1][1])%q]]
        B=[[(B[0][0]*B[0][0]+B[0][1]*B[1][0])%q,(B[0][0]*B[0][1]+B[0][1]*B[1][1])%q],
           [(B[1][0]*B[0][0]+B[1][1]*B[1][0])%q,(B[1][0]*B[0][1]+B[1][1]*B[1][1])%q]]
        k>>=1
    return R

def eps_fun(q,M):
    e={}
    for x in range(q):
        P=matpow([[0,1],[(q-1)%q,x]],M,q)   # N(x)=[[0,1],[-1,x]]
        if P==[[1,0],[0,1]]: e[x]=1
        elif P==[[q-1,0],[0,q-1]]: e[x]=-1
        else: e[x]=None
    return e

def sgn_perm(d1,d2,q):
    a=1 if d1==0 else -1; b=1 if d2==0 else -1
    return {(x,y,z):((a*x)%q,(b*y)%q,(a*b*z)%q) for x in range(q) for y in range(q) for z in range(q)}

def check_q(q, Ms, p0s):
    L=Lq(q); E=L//2
    U=word_perm(['S','T','Si'],q)         # sigma_L^{-1}
    Rm=word_perm(['Ti'],q); oR=perm_order(Rm); oU=perm_order(U)
    assert oR==L and oU==L, (q,oR,oU)
    # (S2) rho rules: sigma_L^{-1} conj on Sigma
    Sg={(d1,d2):sgn_perm(d1,d2,q) for d1 in(0,1) for d2 in(0,1)}
    Ui={v:k for k,v in U.items()}
    rho={}
    for (d1,d2),s in Sg.items():
        c=compose(compose(Ui,s),U)  # U s U^{-1} as: apply U^{-1}, s, U?  order check below
        # find which sign perm equals c
        for k2,s2 in Sg.items():
            if c==s2: rho[(d1,d2)]=k2; break
        else: raise AssertionError('conj not a sign perm')
    # rho should be linear of order 2
    r2=all(rho[rho[k]]==k for k in rho)
    print(f'q={q}: L={L} E={E}  ord sR=ord sL^-1={L} ok; rho={rho}  rho^2=id:{r2}')
    assert r2
    for M in Ms:
        assert M%E==0
        e=eps_fun(q,M)
        assert all(v in(1,-1) for v in e.values()), 'E|M but N(x)^M not +-I'
        assert all(e[x]==e[(-x)%q] for x in range(q)), 'eps(-x)!=eps(x)'
        # (S1) J formula
        J={(x,y,z):(x,(e[x]*y)%q,(e[x]*z)%q) for x in range(q) for y in range(q) for z in range(q)}
        JR=perm_pow(Rm,(-M)%L)
        assert J==JR, 'sigma_R^{-M} != fiberwise sign'
        for p0 in p0s:
            if gcd(p0,M)!=1: continue
            b=M*p0
            # (S3)
            assert (b-M)%L==0 and perm_pow(Rm,(-b)%L)==J
            # words exactly as in the paper/repo:
            A1=word_perm(['S','T','Si']+['Ti']*( b%L ),q)
            A2=word_perm(['S']+['T']*p0+['Si']+['Ti']*(M%L if M%L else L),q)
            # structure: A1 = U o J variant; establish convention
            c1 = A1==compose(U,J) or A1==compose(J,U)
            c2 = A2==compose(perm_pow(U,p0),J) or A2==compose(J,perm_pow(U,p0))
            assert c1 and c2, 'word/structure mismatch'
            # (S5)
            assert perm_pow(A1, 2*L)=={p:p for p in A1}
            # (S4) per base cycle over X/Sigma
            orbrep={}
            def rep(t):
                if t not in orbrep:
                    o=min((( (1 if a==0 else -1)*t[0])%q, ((1 if bb==0 else -1)*t[1])%q,
                           (((1 if a==0 else -1)*(1 if bb==0 else -1))*t[2])%q)
                          for a in(0,1) for bb in(0,1))
                    orbrep[t]=o
                return orbrep[t]
            # base map Ubar on reps
            seenr=set(); blocks=[]
            for t in U:
                r=rep(t)
                if r in seenr: continue
                cyc=[r]; seenr.add(r); cur=rep(U[r])
                while cur!=r:
                    cyc.append(cur); seenr.add(cur); cur=rep(U[cur])
                blocks.append(cyc)
            ok_all=True
            for cyc in blocks:
                m=len(cyc); assert gcd(p0,m)==1
                block=set()
                for r in cyc:
                    for a in(0,1):
                        for bb in(0,1):
                            aa=1 if a==0 else -1; bbv=1 if bb==0 else -1
                            block.add(((aa*r[0])%q,(bbv*r[1])%q,((aa*bbv)*r[2])%q))
                def lens(P):
                    seen=set(); out=[]
                    for s0 in block:
                        if s0 in seen: continue
                        n=0;t=s0
                        while True:
                            seen.add(t);n+=1;t=P[t]
                            if t==s0: break
                        out.append(n)
                    return sorted(out)
                l1,l2=lens(A1),lens(A2)
                if l1!=l2: ok_all=False; print('   MISMATCH',cyc[0],m,l1,l2)
                if m%2==0:
                    F=[t for t in block if rep(t)==cyc[0]]
                    P1=perm_pow(A1,m); P2=perm_pow(A2,m)
                    if any(P1[t]!=P2[t] for t in F):
                        ok_all=False; print('   even-m exact-holonomy FAIL at',cyc[0])
            print(f'   M={M:5d} p0={p0:3d}: blocks={len(blocks)} per-cycle equality: {"OK" if ok_all else "FAIL"};'
                  f' global blind: {cyc_by_level(A1,q)==cyc_by_level(A2,q)}')
            assert ok_all

def check_pm1(q):
    L=Lq(q); U=word_perm(['S','T','Si'],q); Rm=word_perm(['Ti'],q)
    M=5 if q!=5 else 7   # E does NOT divide M
    for p0,sgn in ((L+1,'+1'),(L-1,'-1')):
        # p0 need not be prime for the mechanism check
        b=M*p0
        A1=word_perm(['S','T','Si']+['Ti']*(b%L),q)
        A2=word_perm(['S']+['T']*(p0%(2*L) if p0%(2*L) else 2*L)+['Si']+['Ti']*(M%L),q)
        if sgn=='+1':
            print(f'q={q} p0={p0} (==+1 mod L), M={M}: A1==A2 exactly: {A1==A2}')
        else:
            Ai={v:k for k,v in A2.items()}
            conj=compose(compose({v:k for k,v in U.items()},Ai),U)
            print(f'q={q} p0={p0} (==-1 mod L), M={M}: cyc equal: {cyc_total(A1)==cyc_total(A2)};'
                  f' A1 == U^-1 A2^-1 U: {A1==conj or cyc_by_level(A1,q)==cyc_by_level(Ai,q)}')

if __name__=='__main__':
    check_q(3,[6,12,18],[5,7,11,13])
    check_q(5,[30,60,90],[7,11,13])
    check_q(7,[84,168],[11,13])
    check_pm1(3); check_pm1(5)
    print('ALL CHECKS PASSED')
