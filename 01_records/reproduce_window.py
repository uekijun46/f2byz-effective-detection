#!/usr/bin/env python3
"""reproduce_window.py -- rescan of the witness-prime window (paper §6).
For each trace T in [t0,t1], det=+1, D=T^2-4: enumerate narrow classes,
fold by class-inversion (GL_2 Z-conjugacy), build automorph matrices,
and find for every genuine pair the least detecting prime q in QLIST
(per-level cycle types on F_q^3). Reports p_min(D) = max over pairs.
Usage: python3 reproduce_window.py t0 t1 out.json
"""
import os, sys
_R = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [_R, os.path.join(_R, "src"), os.path.dirname(os.path.abspath(__file__))]

import sys, json, time
sys.path.insert(0, __file__.rsplit('/',1)[0]+'/src')
from markoff_rigidity import classes, canon, automorph, decomp

QLIST=[3,5,7,11,13,17,19,23]

def vstep(g,x,y,z,q):
    if g=='T':  return (x, z, (x*z-y)%q)
    if g=='Ti': return (x, (x*y-z)%q, y)
    if g in ('S','Si'): return (y, x, (x*y-z)%q)
    return (x, y, (x*y-z)%q)   # 'J'

_letters={}
def letter_step(g,q):
    key=(g,q)
    if key not in _letters:
        pts=[(x,y,z) for x in range(q) for y in range(q) for z in range(q)]
        idx={p:i for i,p in enumerate(pts)}
        _letters[key]=[idx[vstep(g,*p,q)] for p in pts]
        _letters[('pts',q)]=pts
    return _letters[key]

def word_perm(w,q):
    n=q*q*q
    perm=list(range(n))
    for g in reversed(w):
        step=letter_step(g,q)
        perm=[step[i] for i in perm]
    return perm

def level_ct(perm,q):
    pts=_letters[('pts',q)]
    n=len(perm); seen=bytearray(n); lev={}
    for s in range(n):
        if seen[s]: continue
        c=0; t=s
        while True:
            seen[t]=1; c+=1; t=perm[t]
            if t==s: break
        x,y,z=pts[s]; k=(x*x+y*y+z*z-x*y*z)%q
        lev.setdefault(k,[]).append(c)
    return {k:tuple(sorted(v)) for k,v in lev.items()}

def scan(t0,t1):
    out={}
    for T in range(t0,t1+1):
        D=T*T-4
        cls=classes(D)
        # fold to GL_2(Z)-classes: identify class(M), class(M^{-1}),
        # class(WMW^{-1}), class(WM^{-1}W^{-1}), W=diag(1,-1), via the
        # direct formula class([[a,b],[c,d]]) = canon((c, d-a, -b), D).
        def mcls(M):
            (a,b),(c,d)=M
            return canon((c, d-a, -b), D)
        reps=[]; seen=set()
        for c0 in cls:
            if c0 in seen: continue
            M=automorph(c0,T,1)
            (a,b),(c,d)=M
            Minv=[[d,-b],[-c,a]]
            Mw=[[a,-b],[-c,d]]
            Mwi=[[d,b],[c,a]]
            orb={mcls(M),mcls(Minv),mcls(Mw),mcls(Mwi)}
            seen|=orb
            reps.append(c0)
        words={c: decomp(automorph(c,T,1)) for c in reps}
        pairs=[(reps[i],reps[j]) for i in range(len(reps)) for j in range(i+1,len(reps))]
        pminD=None; pair_data=[]
        cache={}
        for (c1,c2) in pairs:
            wit=None
            for q in QLIST:
                for c in (c1,c2):
                    if (c,q) not in cache:
                        cache[(c,q)]=level_ct(word_perm(words[c],q),q)
                if cache[(c1,q)]!=cache[(c2,q)]:
                    wit=q; break
            pair_data.append(wit)
            if wit is None: pminD='FAIL'
        if pairs:
            if pminD!='FAIL': pminD=max(w for w in pair_data)
        out[T]={'h_narrow':len(cls),'n_gl':len(reps),'n_pairs':len(pairs),
                'witnesses':pair_data,'p_min_D':pminD}
        print(T, out[T], flush=True)
    return out

if __name__=='__main__':
    t0,t1,fn=int(sys.argv[1]),int(sys.argv[2]),sys.argv[3]
    t=time.time(); res=scan(t0,t1)
    json.dump(res,open(fn,'w'))
    print('elapsed',round(time.time()-t,1),'s')
