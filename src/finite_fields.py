"""Arithmetic of GF(p^n) and Vieta cycle types on X(F_{p^n}) for prime powers q.

Conventions match src/markoff_rigidity.py.
"""
import os
# Conventions:
#   letters 'T','Ti','S','Si','J';  S and Si act identically: (x,y,z)->(y,x,xy-z).
# Elements of GF(p^n) are encoded as integers 0..q-1 (base-p digit vectors,
# low digit = constant term).  MUL is a precomputed q x q table.

def _polmul(a, b, p):
    r = [0]*(len(a)+len(b)-1)
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                r[i+j] = (r[i+j] + ai*bj) % p
    return r

def _polmod(a, f, p):
    a = a[:]
    df, dm = len(f)-1, len(a)-1
    inv_lead = pow(f[-1], p-2, p)
    while dm >= df:
        c = a[dm]*inv_lead % p
        if c:
            for i in range(df+1):
                a[dm-df+i] = (a[dm-df+i] - c*f[i]) % p
        dm -= 1
        while dm >= 0 and a[dm] == 0 and dm >= df:
            dm -= 1
        a = a[:max(dm+1, df)]
    a = a[:df]
    return a + [0]*(df-len(a))

def _polgcd(a, b, p):
    a, b = a[:], b[:]
    def deg(u):
        d = len(u)-1
        while d >= 0 and u[d] == 0: d -= 1
        return d
    while deg(b) >= 0:
        db = deg(b)
        inv = pow(b[db], p-2, p)
        r = a[:]
        while deg(r) >= db:
            dr = deg(r)
            c = r[dr]*inv % p
            for i in range(db+1):
                r[dr-db+i] = (r[dr-db+i] - c*b[i]) % p
        a, b = b, r
    return a

def _xq_pow_mod(f, p, e):
    # x^(p^e) mod f, by repeated p-th powering
    df = len(f)-1
    r = [0,1] + [0]*(df-2) if df >= 2 else _polmod([0,1], f, p)
    r = _polmod(r + [0]*max(0, df-len(r)), f, p)
    for _ in range(e):
        s = [0]*df  # r^p mod f via binary pow
        base, ee = r[:], p
        s = [1] + [0]*(df-1)
        while ee:
            if ee & 1:
                s = _polmod(_polmul(s, base, p), f, p)
            base = _polmod(_polmul(base, base, p), f, p)
            ee >>= 1
        r = s
    return r

def find_irreducible(p, n):
    """Smallest monic irreducible of degree n over F_p (coeff list, low->high)."""
    if n == 1:
        return [0, 1]
    from itertools import product
    for tail in product(range(p), repeat=n):
        f = list(tail) + [1]
        if f[0] == 0:
            continue
        ok = True
        for e in range(1, n//2 + 1):
            xq = _xq_pow_mod(f, p, e)          # x^(p^e) mod f
            diff = xq[:]                        # x^(p^e) - x
            diff[1] = (diff[1] - 1) % p
            g = _polgcd(f, diff, p)
            dg = max((i for i, c in enumerate(g) if c), default=-1)
            if dg != 0:
                ok = False; break
        if ok:
            return f
    raise RuntimeError("no irreducible found")

def build_gf(p, n):
    """Return (q, ADD, SUB, MUL) tables for GF(p^n); elements 0..q-1."""
    q = p**n
    f = find_irreducible(p, n)
    def enc(v): return sum(c*p**i for i, c in enumerate(v))
    def dec(x): return [(x // p**i) % p for i in range(n)]
    ADD = [[0]*q for _ in range(q)]
    MUL = [[0]*q for _ in range(q)]
    for a in range(q):
        va = dec(a)
        for b in range(a, q):
            vb = dec(b)
            s = enc([(va[i]+vb[i]) % p for i in range(n)])
            ADD[a][b] = ADD[b][a] = s
            m = enc(_polmod(_polmul(va, vb, p), f, p))
            MUL[a][b] = MUL[b][a] = m
    NEG = [enc([(-c) % p for c in dec(a)]) for a in range(q)]
    SUB = [[ADD[a][NEG[b]] for b in range(q)] for a in range(q)]
    return q, ADD, SUB, MUL

def vstep_gf(g, x, y, z, SUB, MUL):
    if g == 'T':  return (x, z, SUB[MUL[x][z]][y])
    if g == 'Ti': return (x, SUB[MUL[x][y]][z], y)
    if g == 'S' or g == 'Si': return (y, x, SUB[MUL[x][y]][z])
    if g == 'J':  return (x, y, SUB[MUL[x][y]][z])
    raise ValueError(g)

def fixed_and_cycletype(word, p, n):
    """(fixed-point count, cycle type) of the Vieta permutation on X(GF(p^n))."""
    from collections import Counter
    q, ADD, SUB, MUL = build_gf(p, n)
    seen = bytearray(q*q*q)
    ct = Counter()
    for x in range(q):
        for y in range(q):
            base = (x*q + y)*q
            for z in range(q):
                if seen[base+z]:
                    continue
                L = 0; cx, cy, cz = x, y, z
                while True:
                    j = (cx*q + cy)*q + cz
                    if seen[j]:
                        break
                    seen[j] = 1; L += 1
                    for g in word:
                        cx, cy, cz = vstep_gf(g, cx, cy, cz, SUB, MUL)
                ct[L] += 1
    return dict(ct).get(1, 0), tuple(sorted(ct.items()))

if __name__ == "__main__":
    # Self-check against S8: central pair (D=40, T=38) has delta(3),delta(9),delta(27)=8,20,32.
    src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'markoff_rigidity.py')).read().replace(
        'if __name__ == "__main__":', 'if False:')
    exec(src)
    D = 40
    cl = classes(D); TU = fundamental_unit_norm1(D)
    wA, wB = (decomp(automorph(f, *TU)) for f in cl[:2])
    for nn in (1, 2, 3):
        fA, _ = fixed_and_cycletype(wA, 3, nn)
        fB, _ = fixed_and_cycletype(wB, 3, nn)
        print(f"n={nn}: q=3^{nn}  fix(A)={fA}  fix(B)={fB}  delta={fA-fB}")
