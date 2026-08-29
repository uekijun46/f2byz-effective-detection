"""X_G = Hom(F_2,G)/Inn(G) with the Nielsen action; Hom counts and cycle types."""
import os, sys
_R = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [_R, os.path.join(_R, "src"), os.path.dirname(os.path.abspath(__file__))]

from levels import mul, inv

def build(G, e):
    idx = {g: i for i, g in enumerate(G)}
    pairs = [(a, b) for a in G for b in G]
    rep = {}
    lab = {}
    for p in pairs:
        if p in lab: continue
        i = len(rep); rep[i] = p
        for h in G:
            hi = inv(h)
            lab[(mul(mul(h, p[0]), hi), mul(mul(h, p[1]), hi))] = i
    return rep, lab

def perm(rep, lab, move):
    return [lab[move(rep[i])] for i in range(len(rep))]

def act_word(rep, lab, word):
    R = lambda p: (p[0], mul(p[0], p[1]))
    L = lambda p: (mul(p[0], p[1]), p[1])
    def f(p):
        for i, a in enumerate(word):
            g = L if i % 2 == 0 else R
            for _ in range(a): p = g(p)
        return p
    return [lab[f(rep[i])] for i in range(len(rep))]

def ctype(p):
    n = len(p); seen = [False]*n; out = []
    for i in range(n):
        if not seen[i]:
            c = 0; j = i
            while not seen[j]: seen[j] = True; j = p[j]; c += 1
            out.append(c)
    return tuple(sorted(out))

def order_of(p):
    n = len(p); idn = list(range(n)); cur = p[:]; k = 1
    while cur != idn:
        cur = [p[cur[i]] for i in range(n)]; k += 1
    return k
