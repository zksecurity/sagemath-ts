"""SageMath side of the ``elliptic_curves`` property-test area.

Cases: tests/property/cases/elliptic_curves.cases.json
"""

from sage.all import *


def _make_curve(p, a, b):
    """Create an elliptic curve y^2 = x^3 + ax + b over GF(p)."""
    F = GF(int(p))
    return EllipticCurve(F, [F(int(a)), F(int(b))])


def _make_point(E, x, y):
    """Create a point on elliptic curve E."""
    F = E.base_field()
    return E(F(int(x)), F(int(y)))


def _format_point(P):
    """Format an elliptic curve point as string."""
    if P.is_zero():
        return '(0 : 1 : 0)'
    x, y = P.xy()
    return f'({x} : {y} : 1)'


def ec_point_add(p, a, b, x1, y1, x2, y2):
    """Add two points on elliptic curve y^2 = x^3 + ax + b over GF(p)."""
    E = _make_curve(p, a, b)
    P1 = _make_point(E, x1, y1)
    P2 = _make_point(E, x2, y2)
    R = P1 + P2
    return _format_point(R)


def ec_point_add_identity(p, a, b, x, y):
    """Test P + O = P (adding identity)."""
    E = _make_curve(p, a, b)
    P = _make_point(E, x, y)
    O = E(0)
    R = P + O
    return _format_point(R)


def ec_scalar_mul(p, a, b, x, y, n):
    """Compute n*P on elliptic curve."""
    E = _make_curve(p, a, b)
    P = _make_point(E, x, y)
    R = int(n) * P
    return _format_point(R)


def ec_point_neg(p, a, b, x, y):
    """Compute -P on elliptic curve."""
    E = _make_curve(p, a, b)
    P = _make_point(E, x, y)
    R = -P
    return _format_point(R)


def ec_point_order(p, a, b, x, y):
    """Compute the order of point P."""
    E = _make_curve(p, a, b)
    P = _make_point(E, x, y)
    return str(P.order())


def ec_ellcard(p, a, b):
    """Compute the cardinality (order) of the elliptic curve group."""
    E = _make_curve(p, a, b)
    return str(E.cardinality())


def ec_discriminant(p, a, b):
    """Compute the discriminant of the elliptic curve."""
    E = _make_curve(p, a, b)
    return str(E.discriminant())


def ec_j_invariant(p, a, b):
    """Compute the j-invariant of the elliptic curve."""
    E = _make_curve(p, a, b)
    return str(E.j_invariant())


def ec_trace_of_frobenius(p, a, b):
    """Compute the trace of Frobenius."""
    E = _make_curve(p, a, b)
    return str(E.trace_of_frobenius())


def ec_is_supersingular(p, a, b):
    """Check if the curve is supersingular."""
    E = _make_curve(p, a, b)
    return E.is_supersingular()


def ec_is_ordinary(p, a, b):
    """Check if the curve is ordinary (not supersingular)."""
    E = _make_curve(p, a, b)
    return E.is_ordinary()


def ec_embedding_degree(p, a, b, n):
    """Compute the embedding degree for n-torsion."""
    E = _make_curve(p, a, b)
    # The embedding degree is the smallest k such that n | (p^k - 1)
    n = int(n)
    p_int = int(p)
    k = 1
    pk = p_int
    while (pk - 1) % n != 0:
        k += 1
        pk = pk * p_int
        if k > 1000:  # Safety limit
            return str(-1)
    return str(k)


def ec_division_polynomial(p, a, b, n):
    """Compute the n-th division polynomial."""
    E = _make_curve(p, a, b)
    psi = E.division_polynomial(int(n))
    # Format as comma-separated coefficients to match TypeScript
    coeffs = [str(int(c)) for c in psi.list()]
    return ','.join(coeffs)


def ec_isogeny_degree(p, a, b, deg):
    """Compute an isogeny of given degree from the curve (returns degree if valid)."""
    try:
        E = _make_curve(p, a, b)
        deg = int(deg)
        # Find a point of order deg for the isogeny kernel
        # For a degree-d isogeny, we need a subgroup of order d
        # Try to find such a point
        order = E.cardinality()
        if order % deg != 0:
            return 'no_isogeny'

        cofactor = order // deg
        for _ in range(100):
            P = E.random_point()
            Q = cofactor * P
            if not Q.is_zero() and (deg * Q).is_zero():
                # Found a point of order dividing deg
                # Create isogeny
                phi = E.isogeny(Q)
                return str(phi.degree())
        return 'no_point_found'
    except Exception as e:
        return f'error: {str(e)}'


def ec_isogeny_codomain_j(p, a, b, deg):
    """Compute j-invariant of codomain of isogeny of given degree."""
    try:
        E = _make_curve(p, a, b)
        deg = int(deg)
        order = E.cardinality()
        if order % deg != 0:
            return 'no_isogeny'

        cofactor = order // deg
        for _ in range(100):
            P = E.random_point()
            Q = cofactor * P
            if not Q.is_zero() and (deg * Q).is_zero():
                phi = E.isogeny(Q)
                return str(phi.codomain().j_invariant())
        return 'no_point_found'
    except Exception as e:
        return f'error: {str(e)}'


def ec_torsion_points_count(p, a, b):
    """Count the number of points on the curve (same as cardinality for finite fields)."""
    E = _make_curve(p, a, b)
    # For finite fields, all points are torsion
    return str(E.cardinality())


def ec_is_on_curve(p, a, b, x, y):
    """Check if point (x, y) is on the curve."""
    E = _make_curve(p, a, b)
    F = E.base_field()
    x_val = F(int(x))
    y_val = F(int(y))
    # Check: y^2 == x^3 + a*x + b
    lhs = y_val**2
    rhs = x_val**3 + F(int(a))*x_val + F(int(b))
    return lhs == rhs


def ec_lift_x(p, a, b, x):
    """Lift an x-coordinate to a point on the curve (if exists)."""
    E = _make_curve(p, a, b)
    F = E.base_field()
    x_val = F(int(x))
    try:
        P = E.lift_x(x_val)
        return _format_point(P)
    except ValueError:
        return 'no_point'


def ec_associativity_check(p, a, b, x1, y1, x2, y2, x3, y3):
    """Check (P + Q) + R == P + (Q + R)."""
    E = _make_curve(p, a, b)
    P = _make_point(E, x1, y1)
    Q = _make_point(E, x2, y2)
    R = _make_point(E, x3, y3)
    left = (P + Q) + R
    right = P + (Q + R)
    return left == right


def ec_inverse_check(p, a, b, x, y):
    """Check P + (-P) == O."""
    E = _make_curve(p, a, b)
    P = _make_point(E, x, y)
    R = P + (-P)
    return R.is_zero()


FUNCTIONS = {
    # Point operations
    'point_add': ec_point_add,
    'point_add_identity': ec_point_add_identity,
    'scalar_mul': ec_scalar_mul,
    'point_neg': ec_point_neg,
    'point_order': ec_point_order,
    # Curve properties
    'ellcard': ec_ellcard,
    'discriminant': ec_discriminant,
    'j_invariant': ec_j_invariant,
    'trace_of_frobenius': ec_trace_of_frobenius,
    'is_supersingular': ec_is_supersingular,
    'is_ordinary': ec_is_ordinary,
    'embedding_degree': ec_embedding_degree,
    # Division polynomials
    'division_polynomial': ec_division_polynomial,
    # Isogeny
    'isogeny_degree': ec_isogeny_degree,
    'isogeny_codomain_j': ec_isogeny_codomain_j,
    # Torsion
    'torsion_points_count': ec_torsion_points_count,
    # Point validation
    'is_on_curve': ec_is_on_curve,
    'lift_x': ec_lift_x,
    # Group law checks
    'associativity_check': ec_associativity_check,
    'inverse_check': ec_inverse_check,
}
