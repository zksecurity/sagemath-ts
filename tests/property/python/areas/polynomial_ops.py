"""SageMath side of the ``polynomial_ops`` property-test area (ZZ[x] and GF(p)[x]).

Cases: tests/property/cases/polynomial_ops.cases.json
"""

from sage.all import *

from ._helpers import (
    poly_derivative,
    poly_eval,
    poly_factor,
    poly_gcd,
    poly_is_irreducible,
    poly_mul,
    poly_pow,
    poly_quo_rem,
    poly_roots,
)


def _make_poly_zz(coeffs):
    """Create a polynomial over ZZ from coefficient list [c0, c1, ..., cn]."""
    R = PolynomialRing(ZZ, 'x')
    return R(list(map(Integer, coeffs)))


def poly_gcd_zz(coeffs1, coeffs2):
    """Compute GCD of two polynomials over ZZ."""
    f = _make_poly_zz(coeffs1)
    g = _make_poly_zz(coeffs2)
    result = gcd(f, g)
    # Make content-free (primitive)
    if result != 0:
        content = result.content()
        if content != 0:
            result = result // content
        # Make leading coefficient positive
        if result.leading_coefficient() < 0:
            result = -result
    return [int(c) for c in result.list()] if result != 0 else [0]


def poly_eval_zz(coeffs, x):
    """Evaluate polynomial over ZZ at a point."""
    f = _make_poly_zz(coeffs)
    return int(f(Integer(x)))


FUNCTIONS = {
    'poly_gcd_zz': poly_gcd_zz,
    'poly_eval_zz': poly_eval_zz,
    'poly_gcd_ff': poly_gcd,
    'poly_eval_ff': poly_eval,
    'poly_factor_ff': poly_factor,
    'poly_roots_ff': poly_roots,
    'poly_derivative_ff': poly_derivative,
    'poly_is_irreducible_ff': poly_is_irreducible,
    'poly_mul_ff': poly_mul,
    'poly_quo_rem_ff': poly_quo_rem,
    'poly_pow_ff': poly_pow,
}
