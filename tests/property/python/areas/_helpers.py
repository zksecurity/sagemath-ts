"""Shared helpers used by more than one area module.

Only put something here when two or more areas genuinely need it (today: the
GF(p) polynomial helpers shared by ``polynomials`` and ``polynomial_ops``).
Area-private helpers belong in the area module itself so that agents adding
areas in parallel never touch the same file.
"""

from sage.all import *


def make_poly(p, coeffs):
    """Create a polynomial over GF(p) from coefficient list [c0, c1, ..., cn]."""
    F = GF(p)
    R = PolynomialRing(F, 'x')
    return R(list(map(F, coeffs)))


def poly_to_list(poly):
    """Convert polynomial to coefficient list."""
    if poly == 0:
        return []
    return [int(c) for c in poly.list()]


def poly_add(p, coeffs1, coeffs2):
    """Add two polynomials over GF(p)."""
    f = make_poly(p, coeffs1)
    g = make_poly(p, coeffs2)
    return poly_to_list(f + g)


def poly_mul(p, coeffs1, coeffs2):
    """Multiply two polynomials over GF(p)."""
    f = make_poly(p, coeffs1)
    g = make_poly(p, coeffs2)
    return poly_to_list(f * g)


def poly_quo_rem(p, coeffs1, coeffs2):
    """Compute quotient and remainder of polynomial division over GF(p)."""
    f = make_poly(p, coeffs1)
    g = make_poly(p, coeffs2)
    q, r = f.quo_rem(g)
    return (poly_to_list(q), poly_to_list(r))


def poly_mod(p, coeffs1, coeffs2):
    """Compute polynomial modulo over GF(p)."""
    f = make_poly(p, coeffs1)
    g = make_poly(p, coeffs2)
    return poly_to_list(f % g)


def poly_gcd(p, coeffs1, coeffs2):
    """Compute GCD of two polynomials over GF(p)."""
    f = make_poly(p, coeffs1)
    g = make_poly(p, coeffs2)
    result = gcd(f, g)
    # Make monic
    if result != 0:
        result = result.monic()
    return poly_to_list(result)


def poly_eval(p, coeffs, x):
    """Evaluate polynomial at a point over GF(p)."""
    f = make_poly(p, coeffs)
    F = GF(p)
    return int(f(F(x)))


def poly_factor(p, coeffs):
    """Factor polynomial over GF(p). Returns list of (factor_coeffs, multiplicity)."""
    f = make_poly(p, coeffs)
    factorization = f.factor()
    result = []
    for fac, mult in factorization:
        # Make factor monic
        monic_fac = fac.monic()
        result.append((poly_to_list(monic_fac), int(mult)))
    # Sort by degree then lexicographically
    result.sort(key=lambda x: (len(x[0]), x[0]))
    return result


def poly_derivative(p, coeffs):
    """Compute derivative of polynomial over GF(p)."""
    f = make_poly(p, coeffs)
    return poly_to_list(f.derivative())


def poly_is_irreducible(p, coeffs):
    """Test if polynomial is irreducible over GF(p)."""
    f = make_poly(p, coeffs)
    return f.is_irreducible()


def poly_roots(p, coeffs):
    """Find roots of polynomial over GF(p). Returns sorted list of roots."""
    f = make_poly(p, coeffs)
    roots = f.roots(multiplicities=False)
    return sorted([int(r) for r in roots])


def poly_pow(p, coeffs, n):
    """Compute polynomial power over GF(p)."""
    f = make_poly(p, coeffs)
    return poly_to_list(f ** int(n))
