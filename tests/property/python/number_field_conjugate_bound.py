"""Conjugate-count bound and its word-polynomial dependencies versus bundled PARI."""
from pari_qx_factor import pari_qx_factor


def nf_conjugate_bound(op, p, start, coefficients):
    operation = 16 if int(op) < 2 else 18 if int(op) == 2 else 17
    return pari_qx_factor(operation, p, start, 0, 0, 1, coefficients, [], [])
