"""Verification of prop:below24: every group of order less than 24 is
congruence, with N = exp G.

Let N_G = {N normal in F_2 : F_2/N = G}, which is Epi(F_2,G)/Aut(G), that is,
the set of generating pairs of G up to Aut(G). The group Out(F_2) = GL_2(Z)
acts on N_G by Nielsen moves. The group G is congruence of level N when the
kernel of this action contains the principal congruence subgroup Gamma(N),
equivalently when the action factors through GL_2(Z/N).

The test is exhaustive, not a sample. Traverse GL_2(Z/N) from the identity,
carrying with each matrix the permutation of N_G that a word reaching it
induces. If the traversal closes with no conflict, then the assignment f
satisfies f(Mg) = f(M) f(g) for every M and every generator g, so induction on
the length of a word in the generators shows that f is a homomorphism on
GL_2(Z/N); the action therefore factors, and Gamma(N) lies in the kernel. A
conflict is a proof that the action does not factor through GL_2(Z/N).

S_4 and SL_2(F_3) are run as controls: the paper asserts that they are not
congruence, so the traversal must report a conflict for them.

Usage:  python3 below24_check.py [--json OUT]
"""
import os, sys
_R = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [_R, os.path.join(_R, "src"), os.path.dirname(os.path.abspath(__file__))]

import json
from itertools import product
from groups_below_24 import catalogue, controls, minimal_generators, extend_hom


# ------------------------------------------------------------ N_G and Aut(G)

def automorphisms(G):
    gens = minimal_generators(G)
    orders = [G.order_of(g) for g in gens]
    cand = [[h for h in G.elements if G.order_of(h) == o] for o in orders]
    auts = []
    for images in product(*cand):
        f = extend_hom(G, G, gens, images)
        if f is not None:
            auts.append(f)
    return auts


def generating_pairs(G):
    n = len(G)
    return [(a, b) for a in G.elements for b in G.elements
            if len(G.closure((a, b))) == n]


def normal_subgroup_set(G):
    """N_G as a list of canonical representatives, together with the map from
    a generating pair to its index."""
    auts = automorphisms(G)
    pairs = generating_pairs(G)
    pos = {p: i for i, p in enumerate(pairs)}
    label = {}
    reps = []
    for p in pairs:
        if p in label:
            continue
        i = len(reps); reps.append(p)
        for f in auts:
            label[(f[p[0]], f[p[1]])] = i
    return reps, label, len(auts), len(pairs)


# ------------------------------------------------------------ Nielsen action

def nielsen_perms(G, reps, label):
    """Permutations of N_G induced by the three generators of GL_2(Z).

    R: (a,b) -> (a, ab)      L: (a,b) -> (ab, b)      E: (a,b) -> (a, b^{-1})
    """
    n = len(reps)
    def perm(move):
        return tuple(label[move(p)] for p in reps)
    R = perm(lambda p: (p[0], G.mul(p[0], p[1])))
    L = perm(lambda p: (G.mul(p[0], p[1]), p[1]))
    E = perm(lambda p: (p[0], G.inv(p[1])))
    return [R, L, E], n


MATS = [((1, 1), (0, 1)),      # R
        ((1, 0), (1, 1)),      # L
        ((1, 0), (0, -1))]     # E


def gl2_mod(N):
    """The elements of GL_2(Z/N), as tuples ((a,b),(c,d))."""
    units = {u for u in range(N) if _gcd(u, N) == 1}
    out = []
    for a, b, c, d in product(range(N), repeat=4):
        if (a * d - b * c) % N in units:
            out.append(((a, b), (c, d)))
    return out


def _gcd(a, b):
    while b:
        a, b = b, a % b
    return a


def matmul(M, K, N):
    (a, b), (c, d) = M
    (p, q), (r, s) = K
    return (((a * p + b * r) % N, (a * q + b * s) % N),
            ((c * p + d * r) % N, (c * q + d * s) % N))


def compose(sigma, tau):
    """The permutation i -> tau[sigma[i]], that is, sigma first."""
    return tuple(tau[i] for i in sigma)


def factors_through(perms, N, degree):
    """Traverse GL_2(Z/N) from the identity. Return (True, size) if the
    assignment is consistent everywhere, and (False, witness) otherwise."""
    ident = tuple(range(degree))
    gens = [tuple(((m[0][0]) % N, (m[0][1]) % N)) for m in MATS]
    G1 = [(((m[0][0]) % N, (m[0][1]) % N), ((m[1][0]) % N, (m[1][1]) % N))
          for m in MATS]
    start = ((1 % N, 0), (0, 1 % N))
    seen = {start: ident}
    frontier = [start]
    while frontier:
        nxt = []
        for M in frontier:
            sig = seen[M]
            for g, tau in zip(G1, perms):
                K = matmul(M, g, N)
                new = compose(sig, tau)
                if K in seen:
                    if seen[K] != new:
                        return False, (M, g)
                else:
                    seen[K] = new; nxt.append(K)
        frontier = nxt
    return True, len(seen)


# ------------------------------------------------------------ driver

def run_one(G):
    N = G.exponent()
    reps, label, naut, npairs = normal_subgroup_set(G)
    row = dict(name=G.name, order=len(G), exponent=N,
               n_aut=naut, n_gen_pairs=npairs, n_NG=len(reps))
    # cross-check: Aut(G) acts freely on generating pairs
    row["free_action_identity"] = (npairs == naut * len(reps))
    if len(reps) == 0:
        row["congruence"] = True
        row["note"] = "G is not 2-generated, so N_G is empty"
        return row
    perms, degree = nielsen_perms(G, reps, label)
    ok, info = factors_through(perms, N, degree)
    row["congruence"] = ok
    row["gl2_mod_N_visited"] = info if ok else None
    row["gl2_mod_N_order"] = len(gl2_mod(N)) if N <= 24 else None
    return row


def main(argv):
    out = []
    for G in catalogue():
        r = run_one(G); r["role"] = "claim"
        out.append(r)
        print("%-14s order %2d  exp %2d  |N_G| = %3d   congruence: %s"
              % (r["name"], r["order"], r["exponent"], r["n_NG"],
                 "yes" if r["congruence"] else "NO"))
    print("-" * 62)
    for G in controls():
        r = run_one(G); r["role"] = "control"
        out.append(r)
        print("%-14s order %2d  exp %2d  |N_G| = %3d   congruence: %s"
              % (r["name"], r["order"], r["exponent"], r["n_NG"],
                 "yes" if r["congruence"] else "NO"))
    claims = [r for r in out if r["role"] == "claim"]
    ctrls = [r for r in out if r["role"] == "control"]
    print("-" * 62)
    print("nonabelian groups of order below 24 tested:", len(claims))
    print("all congruence at N = exp G:", all(r["congruence"] for r in claims))
    print("controls not congruence at N = exp G:",
          all(not r["congruence"] for r in ctrls))
    print("Aut acts freely on generating pairs in every case:",
          all(r["free_action_identity"] for r in out))
    if "--json" in argv:
        path = argv[argv.index("--json") + 1]
        json.dump(out, open(path, "w"), indent=2)
        print("written:", path)


if __name__ == "__main__":
    main(sys.argv[1:])
