"""Cycle type of a word in R,L on N_G = {N normal in F_2 : F_2/N = G}, for a
finite group G given as a list of permutations."""
import os, sys
_R = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [_R, os.path.join(_R, "src"), os.path.dirname(os.path.abspath(__file__))]

from itertools import permutations

def mul(p,q): return tuple(p[q[i]] for i in range(len(q)))
def inv(p):
    r=[0]*len(p)
    for i,v in enumerate(p): r[v]=i
    return tuple(r)

def epi_and_orbits(G, e):
    def closure(a,b):
        seen={e}; fr=[e]
        while fr:
            x=fr.pop()
            for g in (a,b,inv(a),inv(b)):
                y=mul(x,g)
                if y not in seen: seen.add(y); fr.append(y)
        return len(seen)
    epi=[(a,b) for a in G for b in G if closure(a,b)==len(G)]
    a0,b0=epi[0]; auts=[]
    for a1,b1 in epi:
        tab={e:e}; fr=[(e,e)]
        gens=[(a0,a1),(b0,b1),(inv(a0),inv(a1)),(inv(b0),inv(b1))]
        ok=True
        while fr and ok:
            x,y=fr.pop()
            for g,h in gens:
                nx,ny=mul(x,g),mul(y,h)
                if nx not in tab: tab[nx]=ny; fr.append((nx,ny))
                elif tab[nx]!=ny: ok=False; break
        if ok and len(tab)==len(G) and len(set(tab.values()))==len(G): auts.append(tab)
    reps={}; lab={}
    for p in epi:
        if p in lab: continue
        i=len(reps); reps[i]=p
        for t in auts: lab[(t[p[0]],t[p[1]])]=i
    return reps, lab

def ctype_level(reps, lab, word):
    def R(p): return (p[0], mul(p[0],p[1]))
    def L(p): return (mul(p[0],p[1]), p[1])
    def act(w,p):
        for i,a in enumerate(w):
            f = L if i%2==0 else R
            for _ in range(a): p=f(p)
        return p
    n=len(reps)
    perm=[lab[act(word,reps[i])] for i in range(n)]
    seen=[False]*n; out=[]
    for i in range(n):
        if not seen[i]:
            c=0;j=i
            while not seen[j]: seen[j]=True; j=perm[j]; c+=1
            out.append(c)
    return tuple(sorted(out))
