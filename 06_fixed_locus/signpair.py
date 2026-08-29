"""Factorization of an element of GL_2(Z) into Nielsen generators, together
with its action on the pair (a, b) of generators of F_2.
"""
import os, sys
_R = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [_R, os.path.join(_R, "src"), os.path.dirname(os.path.abspath(__file__))]
import itertools
from sympy.combinatorics import Permutation, PermutationGroup
from sympy.combinatorics.named_groups import SymmetricGroup, AlternatingGroup, DihedralGroup, CyclicGroup, AbelianGroup

# ---- GL2Z -> Aut(F2) ---------------------------------------------------
# words: lists of ints, ±1 = a^{±1}, ±2 = b^{±1}
def reduce_word(w):
    out=[]
    for x in w:
        if out and out[-1]==-x: out.pop()
        else: out.append(x)
    return out
def inv_word(w): return [-x for x in reversed(w)]
def subst(w, wa, wb):
    out=[]
    for x in w:
        if x==1: out+=wa
        elif x==-1: out+=inv_word(wa)
        elif x==2: out+=wb
        else: out+=inv_word(wb)
    return reduce_word(out)
# elementary automorphisms and their abelianization matrices (columns = images of a,b)
# R: a->a, b->ab  ; abelianization: e1->e1, e2->e1+e2 : matrix [[1,1],[0,1]]
# L: a->ba? use a->ab? Let's define L: a->ab, b->b : e1->e1+e2 : [[1,0],[1,1]]
# I: a->a^{-1}, b->b : diag(-1,1)
# S: a->b, b->a : swap
GENS={'R':([1],[1,2]), 'L':([1,2],[2]), 'I':([-1],[2]), 'S':([2],[1])}
MATS={'R':((1,1),(0,1)), 'L':((1,0),(1,1)), 'I':((-1,0),(0,1)), 'S':((0,1),(1,0))}
def mmul(X,Y): return tuple(tuple(sum(X[i][k]*Y[k][j] for k in range(2)) for j in range(2)) for i in range(2))
def decompose(M):
    """return list of (g,k) with M = prod M_g^k (left to right)."""
    N=[list(r) for r in M]; ops=[]
    while N[1][0]!=0:
        a,c=N[0][0],N[1][0]
        if a!=0 and abs(a)>=abs(c):
            k=a//c; N[0]=[N[0][j]-k*N[1][j] for j in range(2)]; ops.append(('R',-k))
        elif a!=0:
            k=c//a; N[1]=[N[1][j]-k*N[0][j] for j in range(2)]; ops.append(('L',-k))
        else:
            N=[N[1],N[0]]; ops.append(('S',1))
    a,x=N[0]; d=N[1][1]
    assert abs(a)==1 and abs(d)==1
    k=x*d  # x/d
    if k: N[0]=[N[0][j]-k*N[1][j] for j in range(2)]; ops.append(('R',-k))
    if N[0][0]==-1: ops.append(('I',1)); N[0][0]=1
    if N[1][1]==-1: ops.append(('S',1)); ops.append(('I',1)); ops.append(('S',1)); N[1][1]=1
    return [(g,-k) for g,k in ops]
def aut_from_seq(seq):
    # automorphism phi with abelianization M = prod ; phi = phi_{g1}^{k1} ∘ ... ; check by abelianization
    # composition: (phi∘psi)(a)=phi(psi(a)) ; matrix of phi∘psi = M_phi M_psi (columns = images). So M = M_{g1}^{k1}...M_{gm}^{km} corresponds to phi = phi_{g1}^{k1}∘...∘phi_{gm}^{km}
    wa,wb=[1],[2]   # start with identity, build from the right: phi = phi_{g1}^{k1}∘(...) : apply innermost first
    # compute psi = phi_{gm}^{km}, then phi_{g(m-1)}^{k(m-1)}∘psi = substitute images of psi into phi_{g(m-1)}? 
    # (phi∘psi)(a) = phi(psi(a)) = subst(psi(a), phi(a), phi(b)).
    cur=([1],[2])
    for g,k in reversed(seq):
        ga,gb=GENS[g]
        if k<0:
            # inverse of generator: R^{-1}: a->a, b->a^{-1}b ; L^{-1}: a->ab^{-1}, b->b ; I,S self-inverse
            if g=='R': ga,gb=[1],[-1,2]
            elif g=='L': ga,gb=[1,-2],[2]
        for _ in range(abs(k)):
            cur=(subst(cur[0],ga,gb), subst(cur[1],ga,gb))
    return cur
def abel(w):
    return (sum(1 for x in w if x==1)-sum(1 for x in w if x==-1), sum(1 for x in w if x==2)-sum(1 for x in w if x==-2))
def aut_for_matrix(M):
    seq=decompose(M); phi=aut_from_seq(seq)
    ea,eb=abel(phi[0]),abel(phi[1])
    Mchk=((ea[0],eb[0]),(ea[1],eb[1]))
    assert Mchk==tuple(tuple(r) for r in M), (M,Mchk,seq)
    return phi
# ---- Hom counting -----------------------------------------------------
def evalw(w,g,h):
    r=None
    for x in w:
        e = g if x==1 else (~g if x==-1 else (h if x==2 else ~h))
        r = e if r is None else r*e
    return r if r is not None else g*~g
def hom_count(phi, G, twist=False):
    els=list(G.elements); n=len(els)
    total=0
    for g in els:
        for h in els:
            if twist: pa,pb=evalw(phi[0],~g,~h),evalw(phi[1],~g,~h)
            else: pa,pb=evalw(phi[0],g,h),evalw(phi[1],g,h)
            # count k with k g k^-1 = pa, k h k^-1 = pb
            total+=sum(1 for k in els if k*g*~k==pa and k*h*~k==pb)
    return total
if __name__=="__main__":
    A=[[11,-5],[-9,4]]   # class (-9,7,5) of D=229, det -1, tr 15
    phi=aut_for_matrix(A); print("phi:",phi)
    groups={'S3':SymmetricGroup(3),'D4':DihedralGroup(4),'Q8':None,'A4':AlternatingGroup(4),'D5':DihedralGroup(5),'D6':DihedralGroup(6),'S4':SymmetricGroup(4),
            'Z7xZ3':PermutationGroup([Permutation([1,2,3,4,5,6,0]),Permutation([0,2,4,6,1,3,5])]),
            'A5':AlternatingGroup(5),'SL2F3':PermutationGroup([Permutation(0,1,2,3)(4,5,6,7),Permutation(0,4,3,5)(1,7,2,6)]) }
    for name,G in groups.items():
        if G is None: continue
        c1=hom_count(phi,G); c2=hom_count(phi,G,twist=True)
        print(name, G.order(), c1, c2, "DIFF" if c1!=c2 else "")
