"""Primitivity of a positive word: whether its automorph matrix generates the
full multiplicator ring Z[X]/(X^2 - TX + 1) rather than a proper suborder.
"""
from math import gcd

def word_matrix(w):
    m = (1,0,0,1)
    for i,a in enumerate(w):
        p11,p12,p21,p22 = m
        if i % 2 == 0:   # L^a
            m = (p11+a*p12, p12, p21+a*p22, p22)
        else:            # R^a
            m = (p11, a*p11+p12, p21, a*p21+p22)
    return m

def content(w):
    """g such that the associated form has discriminant (T^2-4)/g^2."""
    p11,p12,p21,p22 = word_matrix(w)
    return gcd(gcd(p21, p22-p11), p12)
