"""Live-SageMath oracle for the ported hyperelliptic-curve surface."""

from sage.all import *


def _poly(K, coeffs):
    R = PolynomialRing(K, 'x')
    return R([K(c) for c in coeffs])


def _ints(xs):
    return ','.join(str(ZZ(x)) for x in xs)


def hyp_summary(p, f_coeffs, h_coeffs, extensions):
    K = GF(ZZ(p))
    f = _poly(K, f_coeffs)
    h = _poly(K, h_coeffs) if h_coeffs else 0
    H = HyperellipticCurve(f, h)
    counts = H.count_points(int(extensions))
    frob = H.frobenius_polynomial()
    return 'g=%s counts=%s frob=%s jac=%s' % (
        H.genus(),
        _ints(counts),
        _ints(frob.list()),
        frob(1),
    )


def hyp_cartier(p, f_coeffs):
    K = GF(ZZ(p))
    H = HyperellipticCurve(_poly(K, f_coeffs))
    cartier = ';'.join(_ints(row) for row in H.Cartier_matrix().rows())
    hasse_witt = ';'.join(_ints(row) for row in H.Hasse_Witt().rows())
    return 'C=%s HW=%s a=%s p=%s' % (
        cartier,
        hasse_witt,
        H.a_number(),
        H.p_rank(),
    )


FUNCTIONS = {
    'hyp_summary': hyp_summary,
    'hyp_cartier': hyp_cartier,
}
