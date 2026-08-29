"""Fixed locus of sigma_M on the character variety, for M in GL_2(Z).

Recovers the defining ideal symbolically, writes a Singular session
fix_<name>.sing, and runs primdecGTZ on it. Requires Singular in the path.
"""
import os, sys
_R = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [_R, os.path.join(_R, "src"), os.path.dirname(os.path.abspath(__file__))]
import subprocess
from signpair import decompose
from sympy import symbols, expand, Poly
x,y,z=symbols('x y z')
def apply(seq):
    # moves for generators (as substitutions on (x,y,z)); powers negative -> use inverse move
    def R(p): return (p[0], p[2], expand(p[0]*p[2]-p[1]))
    def Rinv(p): return (p[0], expand(p[0]*p[1]-p[2]), p[1])   # b -> a^{-1} b : tr(a^{-1}b)=xy-z ; new (x, tr(a^{-1}b), tr(a a^{-1} b)=y)
    def L(p): return (p[2], p[1], expand(p[1]*p[2]-p[0]))
    def Linv(p): return (expand(p[0]*p[1]-p[2]), p[1], p[0])   # a -> a b^{-1}: tr = xy - z ; tr(a b^{-1} b)= x
    def S(p): return (p[1],p[0],p[2])
    def I(p): return (p[0],p[1],expand(p[0]*p[1]-p[2]))
    cur=(x,y,z)
    for g,k in seq:
        f={'R':R if k>0 else Rinv,'L':L if k>0 else Linv,'S':S,'I':I}[g]
        for _ in range(abs(k)): cur=f(cur)
    return cur
for name,M in ([] if __name__!="__main__" else [("A",[[7,6],[1,1]]),("B",[[7,3],[2,1]])]):
    seq=decompose(M); sig=apply(seq)
    # sanity: check kappa preserved
    kap=lambda p: expand(p[0]**2+p[1]**2+p[2]**2-p[0]*p[1]*p[2])
    assert kap(sig)==kap((x,y,z)), "kappa not preserved"
    eqs=[expand(sig[i]-v) for i,v in enumerate((x,y,z))]
    ideal=", ".join(str(e).replace('**','^') for e in eqs)
    sing=f'''LIB "primdec.lib";
ring r=0,(x,y,z),dp;
ideal I={ideal};
list L=primdecGTZ(I);
int i; for(i=1;i<=size(L);i++){{ print("component"); print(std(L[i][2])); print(dim(std(L[i][2]))); }}
quit;
'''
    open(f'fix_{name}.sing','w').write(sing)
    out=subprocess.run(['Singular','-q','fix_'+name+'.sing'],capture_output=True,text=True,timeout=250).stdout
    print("=====",name,seq); print(out)
