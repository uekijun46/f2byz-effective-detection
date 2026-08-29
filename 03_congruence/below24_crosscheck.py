"""Cross-check of the two implementations of the verification of prop:below24.

Runs `below24_check.py` and `below24_check_alt.py`, and compares every field
they both report: the exponent, the order of Aut(G), the number of generating
pairs, the size of N_G, and the verdict. For the groups that pass, it also
checks that the number of matrices visited by the traversal of the first
implementation equals the order of the image of GL_2(Z) in GL_2(Z/N) as
computed by the second, which the two obtain in unrelated ways: one by walking
the group, the other from the formula for |SL_2(Z/N)|.

Usage:  python3 below24_crosscheck.py
"""
import os, sys
_R = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [_R, os.path.join(_R, "src"), os.path.dirname(os.path.abspath(__file__))]

import below24_check as A
import below24_check_alt as B
from groups_below_24 import catalogue, controls, check_catalogue

FIELDS = ["order", "exponent", "n_aut", "n_gen_pairs", "n_NG", "congruence",
          "free_action_identity"]


def main():
    ok, report = check_catalogue()
    print("catalogue of the 25 nonabelian groups:",
          "verified" if ok else "FAILED")
    for line in report:
        if line.startswith("isomorphic"):
            print("   ", line)

    groups = [(G, "claim") for G in catalogue()] + \
             [(G, "control") for G in controls()]
    bad = []
    for G, role in groups:
        a = A.run_one(G)
        b = B.run_one(G)
        diff = [f for f in FIELDS if a.get(f) != b.get(f)]
        if diff:
            bad.append((G.name, diff))
        if a["congruence"] and a["n_NG"] > 0:
            if a["gl2_mod_N_visited"] != b["image_order"]:
                bad.append((G.name, ["visited vs image order"]))
        print("%-14s %-8s %s" % (G.name, role, "agree" if not diff else
                                 "DISAGREE on " + ", ".join(diff)))
    print("-" * 58)
    if bad:
        print("MISMATCHES:", bad)
        return 1
    claims = [A.run_one(G) for G, r in groups if r == "claim"]
    ctrls = [A.run_one(G) for G, r in groups if r == "control"]
    print("the two implementations agree on every field, for all %d groups"
          % len(groups))
    print("every nonabelian group of order below 24 is congruence at N = exp G:",
          all(r["congruence"] for r in claims))
    print("S_4 and SL_2(F_3) are not congruence at N = exp G:",
          all(not r["congruence"] for r in ctrls))
    return 0


if __name__ == "__main__":
    sys.exit(main())
