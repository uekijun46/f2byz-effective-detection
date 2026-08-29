"""Enumeration of GL_2(Z)-conjugacy classes (up to inversion) of Anosov matrices
with characteristic polynomial X^2 - T X + 1, via positive words in L and R."""
from collections import defaultdict


def canon(seq):
    n = len(seq)
    best = None
    for s in (seq, seq[::-1]):
        for i in range(n):
            r = s[i:] + s[:i]
            if best is None or r < best:
                best = r
    return best


def enum_words(Tmax):
    buckets = defaultdict(set)

    def rec(m, seq, nxt):
        p11, p12, p21, p22 = m
        a = 1
        while True:
            if nxt == 0:
                nm = (p11 + a * p12, p12, p21 + a * p22, p22)
            else:
                nm = (p11, a * p11 + p12, p21, a * p21 + p22)
            tr = nm[0] + nm[3]
            if tr > Tmax or a > Tmax:
                break
            ns = seq + (a,)
            if nxt == 1:
                buckets[tr].add(canon(ns))
            rec(nm, ns, 1 - nxt)
            a += 1

    rec((1, 0, 0, 1), tuple(), 0)
    return buckets


if __name__ == "__main__":
    import sys, time
    Tmax = int(sys.argv[1]) if len(sys.argv) > 1 else 60
    t = time.time()
    b = enum_words(Tmax)
    print("Tmax=%d traces=%d classes=%d %.1fs"
          % (Tmax, len(b), sum(len(v) for v in b.values()), time.time() - t))
    for T in (8, 38, 40, 102):
        if T <= Tmax:
            print(T, len(b[T]), sorted(b[T])[:6])
