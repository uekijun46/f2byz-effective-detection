"""Cycle types of the Vieta permutation sigma_A(q) on X(F_q)=F_q^3."""
import numpy as np

def build(q):
    idx = np.arange(q**3)
    x = idx // (q*q); y = (idx // q) % q; z = idx % q
    # R:(x,y,z)->(x,z,xz-y)   L:(x,y,z)->(z,y,yz-x)
    R = (x*q*q + z*q + (x*z - y) % q)
    L = (z*q*q + y*q + (y*z - x) % q)
    return R.astype(np.int32), L.astype(np.int32)

def Lq(q): return q*(q*q-1)//2

def powers(P, n):
    idn = np.arange(len(P), dtype=np.int32)
    out = [idn]
    cur = idn
    for _ in range(n-1):
        cur = P[cur]
        out.append(cur)
    return out

class Q:
    def __init__(self, q):
        self.q = q; self.Lq = Lq(q)
        R, L = build(q)
        self.PR = powers(R, self.Lq)
        self.PL = powers(L, self.Lq)
    def perm(self, word):
        # word = (a1,a2,...,a_2k) -> sigma_L^{a1} o sigma_R^{a2} o ...
        arr = np.arange(self.q**3, dtype=np.int32)
        for i in range(len(word)-1, -1, -1):
            tab = self.PL if i % 2 == 0 else self.PR
            arr = tab[word[i] % self.Lq][arr]
        return arr
    def ctype(self, word):
        p = self.perm(word)
        n = len(p); seen = np.zeros(n, dtype=bool); out = []
        for i in range(n):
            if not seen[i]:
                c = 0; j = i
                while not seen[j]:
                    seen[j] = True; j = p[j]; c += 1
                out.append(c)
        out.sort()
        return tuple(out)
