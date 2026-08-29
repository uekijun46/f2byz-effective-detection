"""Vieta action of words in the Nielsen letters on X(F_q) = F_q^3.

The letters are 'T', 'Ti', 'S', 'Si', 'J'; S and Si act identically. A word is
a list of letters, applied from the right. These are the primitives shared by
the scans that work with per-level cycle types.
"""


def vstep(g, x, y, z, q):
    if g == 'T':  return (x, z, (x * z - y) % q)
    if g == 'Ti': return (x, (x * y - z) % q, y)
    if g in ('S', 'Si'): return (y, x, (x * y - z) % q)
    return (x, y, (x * y - z) % q)   # 'J'

_letters = {}
def letter_step(g, q):
    key = (g, q)
    if key not in _letters:
        pts = [(x, y, z) for x in range(q) for y in range(q) for z in range(q)]
        idx = {p: i for i, p in enumerate(pts)}
        _letters[key] = [idx[vstep(g, *p, q)] for p in pts]
        _letters[('pts', q)] = pts
    return _letters[key]

def word_perm(w, q):
    n = q * q * q
    perm = list(range(n))
    for g in reversed(w):
        step = letter_step(g, q)
        perm = [step[i] for i in perm]
    return perm

def level_ct(perm, q):
    pts = _letters[('pts', q)]
    n = len(perm); seen = bytearray(n); lev = {}
    for s0 in range(n):
        if seen[s0]: continue
        c = 0; t = s0
        while True:
            seen[t] = 1; c += 1; t = perm[t]
            if t == s0: break
        x, y, z = pts[s0]; kk = (x * x + y * y + z * z - x * y * z) % q
        lev.setdefault(kk, []).append(c)
