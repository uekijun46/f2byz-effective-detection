"""Action of Out(F_2) on N_G = {N normal in F_2 : F_2/N = G} for G = S_4, SL_2(F_3)."""
import os, sys
_R = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [_R, os.path.join(_R, "src"), os.path.dirname(os.path.abspath(__file__))]

from itertools import permutations, product

def sym_group(n):
    return list(permutations(range(n)))
def mul(p,q): return tuple(p[q[i]] for i in range(len(q)))
def inv(p):
    r=[0]*len(p)
    for i,v in enumerate(p): r[v]=i
    return tuple(r)

def level(G, mulf, invf, e):
    gen = [(a,b) for a in G for b in G]
    def closure(a,b):
        seen={e}; frontier=[e]
        while frontier:
            x=frontier.pop()
            for g in (a,b,invf(a),invf(b)):
                y=mulf(x,g)
                if y not in seen: seen.add(y); frontier.append(y)
        return seen
    epi=[(a,b) for a,b in gen if len(closure(a,b))==len(G)]
    # quotient by Aut(G); for S_4 and SL_2 F_3 use Inn + full Aut via conjugation-by-all-autos
    # we use the orbit of the full automorphism group computed as permutations of G
    return epi

def orbits_under_aut(epi, G, mulf, invf, e):
    # automorphisms: images of a fixed generating pair -> determined
    a0,b0 = epi[0]
    auts=[]
    for a1,b1 in epi:
        # build the map by word evaluation over the Cayley graph
        table={e:e}; frontier=[(e,e)]
        gens=[(a0,a1),(b0,b1),(invf(a0),invf(a1)),(invf(b0),invf(b1))]
        while frontier:
            x,y=frontier.pop()
            for g,h in gens:
                nx,ny=mulf(x,g),mulf(y,h)
                if nx not in table:
                    table[nx]=ny; frontier.append((nx,ny))
        if len(table)==len(G) and len(set(table.values()))==len(G):
            auts.append(table)
    reps={}; lab={}
    for p in epi:
        if p in lab: continue
        orb={(t[p[0]],t[p[1]]) for t in auts}
        i=len(reps); reps[i]=p
        for o in orb: lab[o]=i
    return reps, lab

def cycle_type(perm):
    n=len(perm); seen=[False]*n; out=[]
    for i in range(n):
        if not seen[i]:
            c=0; j=i
            while not seen[j]: seen[j]=True; j=perm[j]; c+=1
            out.append(c)
    return tuple(sorted(out, reverse=True))

def run(word):
    G=sym_group(4); e=tuple(range(4))
    epi=level(G,mul,inv,e)
    reps,lab=orbits_under_aut(epi,G,mul,inv,e)
    n=len(reps)
    # Nielsen moves on pairs:  R:(a,b)->(a,ab)   L:(a,b)->(ab,b)
    def R(p): return (p[0], mul(p[0],p[1]))
    def L(p): return (mul(p[0],p[1]), p[1])
    def act(word,p):
        # A = L^{a1} R^{a2} ... acting on the pair
        for i,a in enumerate(word):
            f = L if i%2==0 else R
            for _ in range(a): p=f(p)
        return p
    perm=[lab[act(word,reps[i])] for i in range(n)]
    return len(epi), n, cycle_type(perm)

if __name__=="__main__":
    for w in [(1,6),(2,3)]:
        print(w, run(w))
    for w in [(1,6),(2,3)]:
        print(w,'squared', run(w+w))
