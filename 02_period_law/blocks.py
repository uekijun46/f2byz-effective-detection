"""Block structure of the automorph word: period, block count, and the
exponent vector, for the collision statistics of Table 3.
"""
import os, sys
_R = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [_R, os.path.join(_R, "src"), os.path.dirname(os.path.abspath(__file__))]

from math import isqrt
from markoff_rigidity import _cycle, _rho, canon

def cf_word(f, D):
    """partial quotients along the reduction cycle of the class of f.
    Returns the cyclic exponent sequence (the R,L block lengths)."""
    cyc = _cycle(f, D)
    qs = []
    for g in cyc:
        h = _rho(g, D)
        a, b, c = g; a2, b2, c2 = h
        q = (b + b2) // (2 * abs(c)) if c else 0
        qs.append(abs(q))
    return tuple(qs)

def canon_cyclic(t):
    n = len(t)
    return min(tuple(t[i:]+t[:i]) for i in range(n))

if __name__ == '__main__':
    # sanity: principal class of D=T^2-4 should be the two-block word (1, T-2)
    for T in (7, 12, 20, 1350):
        D = T*T-4
        b0 = D % 2
        print(T, "principal:", canon_cyclic(cf_word((1,b0,(b0*b0-D)//4), D))[:8],
              " len:", len(cf_word((1,b0,(b0*b0-D)//4), D)))
    # and the family L^u R^v: form (-u, T-2, v)
    T=2097; D=T*T-4
    for u in (1,5,419,2095):
        v=(T-2)//u
        print(f"T={T} u={u} v={v}:", canon_cyclic(cf_word((-u,T-2,v),D)))
