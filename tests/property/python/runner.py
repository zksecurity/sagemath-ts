#!/usr/bin/env sage
"""
Property test runner for SageMath.
Executes operations with deterministic seeding and outputs results.

This script reads test case definitions from stdin (JSON format) and
executes them in SageMath, outputting results in a format that can be
compared with the TypeScript implementation.

Usage:
    sage runner.py < test_cases.json
    sage runner.py --seed 42 < test_cases.json

Output format:
    {"function": "gcd", "args": [12, 8], "result": "4", "seed": 42}
"""

import json
import sys
import random as py_random
from sage.all import *
from sage.arith.misc import dedekind_psi


def set_seed(seed):
    """Set random seed for reproducibility across both SageMath and Python random."""
    set_random_seed(seed)
    py_random.seed(seed)


def generate_arg(generator_spec, seed):
    """
    Generate an argument based on the generator specification.

    Supported generators:
    - randomBigint(min, max): Random integer in range [min, max]
    - randomPrime(min, max): Random prime in range [min, max]
    - fixedValue(value): Fixed value
    - randomList(generator, length): List of random values

    Note: Uses Python's random module (Mersenne Twister) to match TypeScript.
    """
    if generator_spec.startswith('randomBigint('):
        # Parse: randomBigint(min, max)
        params = generator_spec[13:-1]  # Remove 'randomBigint(' and ')'
        parts = params.split(',')
        min_val = int(parts[0].strip())
        max_val = int(parts[1].strip())
        # Use Python's random.randint to match TypeScript's Mersenne Twister
        return Integer(py_random.randint(min_val, max_val))

    elif generator_spec.startswith('randomPrime('):
        # Parse: randomPrime(min, max)
        params = generator_spec[12:-1]
        parts = params.split(',')
        min_val = int(parts[0].strip())
        max_val = int(parts[1].strip())
        # Use Python's random with rejection sampling to match TypeScript
        max_attempts = 10000
        for _ in range(max_attempts):
            candidate = py_random.randint(min_val, max_val)
            if is_prime(candidate):
                return Integer(candidate)
        # Fallback: find next prime from min
        p = next_prime(min_val - 1)
        return p

    elif generator_spec.startswith('fixedValue('):
        # Parse: fixedValue(value) - value can be integer or array
        value_str = generator_spec[11:-1].strip()
        if value_str.startswith('[') and value_str.endswith(']'):
            # Parse array: [1, 2, 3]
            inner = value_str[1:-1].strip()
            if not inner:
                return []
            parts = inner.split(',')
            return [Integer(p.strip()) for p in parts]
        return Integer(value_str)

    elif generator_spec.startswith('randomList('):
        # Parse: randomList(generator, length)
        params = generator_spec[11:-1]
        # Find the last comma that separates length
        last_comma = params.rfind(',')
        inner_generator = params[:last_comma].strip()
        length = int(params[last_comma+1:].strip())
        return [generate_arg(inner_generator, seed + i) for i in range(length)]

    else:
        # Assume it's a literal value
        return Integer(generator_spec)


def format_result(result):
    """
    Format a SageMath result as a string for comparison.

    This ensures consistent string representation across Python and TypeScript.
    """
    if result is None:
        return 'null'
    if isinstance(result, bool):
        return 'True' if result else 'False'
    elif isinstance(result, Factorization):
        # Format: "2^2 * 3" or "1" for empty
        if len(result) == 0:
            return '1'
        parts = []
        for (p, e) in result:
            if e == 1:
                parts.append(str(p))
            else:
                parts.append(f'{p}^{e}')
        return ' * '.join(parts)
    elif isinstance(result, (list, tuple)):
        # Format lists/tuples to match Python's str() output
        if isinstance(result, tuple):
            return '(' + ', '.join(format_result(x) for x in result) + ')'
        else:
            return '[' + ', '.join(format_result(x) for x in result) + ']'
    elif hasattr(result, '__iter__') and not isinstance(result, str):
        # Other iterables
        return '[' + ', '.join(format_result(x) for x in result) + ']'
    else:
        return str(result)


# Polynomial helper functions for property tests
def _make_poly(p, coeffs):
    """Create a polynomial over GF(p) from coefficient list [c0, c1, ..., cn]."""
    F = GF(p)
    R = PolynomialRing(F, 'x')
    return R(list(map(F, coeffs)))


def _poly_to_list(poly):
    """Convert polynomial to coefficient list."""
    if poly == 0:
        return []
    return [int(c) for c in poly.list()]


def poly_add(p, coeffs1, coeffs2):
    """Add two polynomials over GF(p)."""
    f = _make_poly(p, coeffs1)
    g = _make_poly(p, coeffs2)
    return _poly_to_list(f + g)


def poly_mul(p, coeffs1, coeffs2):
    """Multiply two polynomials over GF(p)."""
    f = _make_poly(p, coeffs1)
    g = _make_poly(p, coeffs2)
    return _poly_to_list(f * g)


def poly_quo_rem(p, coeffs1, coeffs2):
    """Compute quotient and remainder of polynomial division over GF(p)."""
    f = _make_poly(p, coeffs1)
    g = _make_poly(p, coeffs2)
    q, r = f.quo_rem(g)
    return (_poly_to_list(q), _poly_to_list(r))


def poly_mod(p, coeffs1, coeffs2):
    """Compute polynomial modulo over GF(p)."""
    f = _make_poly(p, coeffs1)
    g = _make_poly(p, coeffs2)
    return _poly_to_list(f % g)


def poly_gcd(p, coeffs1, coeffs2):
    """Compute GCD of two polynomials over GF(p)."""
    f = _make_poly(p, coeffs1)
    g = _make_poly(p, coeffs2)
    result = gcd(f, g)
    # Make monic
    if result != 0:
        result = result.monic()
    return _poly_to_list(result)


def poly_eval(p, coeffs, x):
    """Evaluate polynomial at a point over GF(p)."""
    f = _make_poly(p, coeffs)
    F = GF(p)
    return int(f(F(x)))


def poly_factor(p, coeffs):
    """Factor polynomial over GF(p). Returns list of (factor_coeffs, multiplicity)."""
    f = _make_poly(p, coeffs)
    factorization = f.factor()
    result = []
    for fac, mult in factorization:
        # Make factor monic
        monic_fac = fac.monic()
        result.append((_poly_to_list(monic_fac), int(mult)))
    # Sort by degree then lexicographically
    result.sort(key=lambda x: (len(x[0]), x[0]))
    return result


def poly_derivative(p, coeffs):
    """Compute derivative of polynomial over GF(p)."""
    f = _make_poly(p, coeffs)
    return _poly_to_list(f.derivative())


def poly_is_irreducible(p, coeffs):
    """Test if polynomial is irreducible over GF(p)."""
    f = _make_poly(p, coeffs)
    return f.is_irreducible()


def poly_roots(p, coeffs):
    """Find roots of polynomial over GF(p). Returns sorted list of roots."""
    f = _make_poly(p, coeffs)
    roots = f.roots(multiplicities=False)
    return sorted([int(r) for r in roots])


def poly_pow(p, coeffs, n):
    """Compute polynomial power over GF(p)."""
    f = _make_poly(p, coeffs)
    return _poly_to_list(f ** int(n))


# ============================================================================
# Elliptic curve helper functions for property tests
# ============================================================================

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


# ============================================================================
# Matrix operations helper functions
# ============================================================================

def matrix_determinant_2x2(a, b, c, d):
    """Compute determinant of 2x2 matrix [[a,b],[c,d]]."""
    M = matrix(ZZ, [[a, b], [c, d]])
    return M.determinant()


def matrix_determinant_3x3(a11, a12, a13, a21, a22, a23, a31, a32, a33):
    """Compute determinant of 3x3 matrix."""
    M = matrix(ZZ, [[a11, a12, a13], [a21, a22, a23], [a31, a32, a33]])
    return M.determinant()


def matrix_determinant_4x4(*args):
    """Compute determinant of 4x4 matrix."""
    M = matrix(ZZ, 4, 4, list(args))
    return M.determinant()


def matrix_rank_2x3(a11, a12, a13, a21, a22, a23):
    """Compute rank of 2x3 matrix."""
    M = matrix(ZZ, [[a11, a12, a13], [a21, a22, a23]])
    return M.rank()


def matrix_rank_3x3(a11, a12, a13, a21, a22, a23, a31, a32, a33):
    """Compute rank of 3x3 matrix."""
    M = matrix(ZZ, [[a11, a12, a13], [a21, a22, a23], [a31, a32, a33]])
    return M.rank()


def matrix_hnf_2x2(a, b, c, d):
    """Compute HNF of 2x2 matrix, return as list of lists."""
    M = matrix(ZZ, [[a, b], [c, d]])
    H = M.hermite_form()
    return [[int(H[i,j]) for j in range(2)] for i in range(2)]


def matrix_hnf_3x3(*args):
    """Compute HNF of 3x3 matrix, return as list of lists."""
    M = matrix(ZZ, 3, 3, list(args))
    H = M.hermite_form()
    return [[int(H[i,j]) for j in range(3)] for i in range(3)]


def matrix_snf_2x2(a, b, c, d):
    """Compute SNF of 2x2 matrix, return diagonal elements."""
    M = matrix(ZZ, [[a, b], [c, d]])
    D, _, _ = M.smith_form()
    return [int(D[i,i]) for i in range(2)]


def matrix_snf_3x3(*args):
    """Compute SNF of 3x3 matrix, return diagonal elements."""
    M = matrix(ZZ, 3, 3, list(args))
    D, _, _ = M.smith_form()
    return [int(D[i,i]) for i in range(3)]


def matrix_lll_2x2(a, b, c, d):
    """Compute LLL-reduced basis of 2x2 matrix, return as list of lists."""
    M = matrix(ZZ, [[a, b], [c, d]])
    L = M.LLL()
    return [[int(L[i,j]) for j in range(2)] for i in range(2)]


def matrix_lll_3x3(*args):
    """Compute LLL-reduced basis of 3x3 matrix, return as list of lists."""
    M = matrix(ZZ, 3, 3, list(args))
    L = M.LLL()
    return [[int(L[i,j]) for j in range(3)] for i in range(3)]


def matrix_elementary_divisors_2x2(a, b, c, d):
    """Compute elementary divisors of 2x2 matrix."""
    M = matrix(ZZ, [[a, b], [c, d]])
    divs = M.elementary_divisors()
    return [int(d) for d in divs]


def matrix_elementary_divisors_3x3(*args):
    """Compute elementary divisors of 3x3 matrix."""
    M = matrix(ZZ, 3, 3, list(args))
    divs = M.elementary_divisors()
    return [int(d) for d in divs]


# ============================================================================
# Special arithmetic functions
# ============================================================================

def arith_binomial(n, k):
    """Compute binomial coefficient C(n,k)."""
    return binomial(n, k)


def arith_fibonacci(n):
    """Compute the n-th Fibonacci number."""
    return fibonacci(n)


def arith_lucas_number(n):
    """Compute the n-th Lucas number using lucas_number2."""
    return lucas_number2(n, 1, -1)


def arith_factorial(n):
    """Compute n!."""
    return factorial(n)


def arith_bernoulli_numerator(n):
    """Return the numerator of the n-th Bernoulli number."""
    return bernoulli(n).numerator()


def arith_bernoulli_denominator(n):
    """Return the denominator of the n-th Bernoulli number."""
    return bernoulli(n).denominator()


def arith_multinomial(*args):
    """Compute multinomial coefficient."""
    return multinomial(*args)


def arith_primitive_root(n):
    """Find a primitive root modulo n."""
    return primitive_root(n)


def arith_nth_prime(n):
    """Return the n-th prime number."""
    return nth_prime(n)


def arith_subfactorial(n):
    """Compute !n (number of derangements)."""
    if n == 0:
        return 1
    if n == 1:
        return 0
    # Use the formula: !n = (n-1) * (!(n-1) + !(n-2))
    a, b = 1, 0  # !0, !1
    for i in range(2, int(n) + 1):
        a, b = b, (i - 1) * (a + b)
    return b


def arith_carmichael_lambda(n):
    """Compute Carmichael's lambda function."""
    return carmichael_lambda(n)


def arith_dedekind_psi(n):
    """Compute Dedekind's psi function."""
    # Use the imported dedekind_psi from sage.arith.misc
    return int(dedekind_psi(n))


# ============================================================================
# Polynomial operations over integers
# ============================================================================

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


# ============================================================================
# Number field operations
# ============================================================================

def nf_quadratic_discriminant(d):
    """Compute the discriminant of Q(sqrt(d))."""
    K = QuadraticField(d)
    return K.discriminant()


def nf_degree(coeffs):
    """Compute the degree of a number field defined by polynomial with given coefficients."""
    R = PolynomialRing(QQ, 'x')
    f = R(list(map(QQ, coeffs)))
    K = NumberField(f, 'a')
    return K.degree()


def nf_quadratic_signature(d):
    """Compute the signature (r1, r2) of Q(sqrt(d))."""
    K = QuadraticField(d)
    return K.signature()


def nf_is_totally_real(d):
    """Check if Q(sqrt(d)) is totally real."""
    K = QuadraticField(d)
    return K.is_totally_real()


def nf_polynomial_discriminant(coeffs):
    """Compute the discriminant of a polynomial."""
    R = PolynomialRing(QQ, 'x')
    f = R(list(map(QQ, coeffs)))
    return f.discriminant()


def nf_fundamental_discriminant(d):
    """Compute the fundamental discriminant associated to d."""
    return fundamental_discriminant(d)


def execute_function(module, function_name, args):
    """
    Execute a SageMath function with the given arguments.

    Args:
        module: Module name (e.g., "arith", "rings.integer")
        function_name: Function name (e.g., "gcd", "is_prime")
        args: List of arguments

    Returns:
        The result of the function call
    """
    # Map module paths to actual SageMath functions
    function_map = {
        'arith': {
            'gcd': gcd,
            'lcm': lcm,
            'xgcd': xgcd,
            'factor': factor,
            'is_prime': is_prime,
            'is_prime_power': is_prime_power,
            'next_prime': next_prime,
            'previous_prime': previous_prime,
            'euler_phi': euler_phi,
            'radical': radical,
            'moebius': moebius,
            'kronecker_symbol': kronecker_symbol,
            'legendre_symbol': legendre_symbol,
            'jacobi_symbol': jacobi_symbol,
            'power_mod': power_mod,
            'inverse_mod': inverse_mod,
            'crt': crt,
            'isqrt': isqrt,
            'is_square': is_square,
            'is_squarefree': is_squarefree,
            'divisors': divisors,
            'number_of_divisors': number_of_divisors,
            'sigma': sigma,
            'prime_range': prime_range,
            'trial_division': trial_division,
            'squarefree_part': squarefree_part,
            'prime_factors': prime_factors,
            'valuation': valuation,
        },
        'finite_fields': {
            'ff_add': lambda p, a, b: GF(p)(a) + GF(p)(b),
            'ff_mul': lambda p, a, b: GF(p)(a) * GF(p)(b),
            'ff_inv': lambda p, a: GF(p)(a)**(-1),
            'ff_pow': lambda p, a, n: GF(p)(a)**n,
            'ff_pow_neg': lambda p, a, n: GF(p)(a)**n,
            'sqrt_mod': lambda a, p: (lambda F, x: (lambda r: p - r if r > (p-1)//2 else r)(Integer(x.sqrt())) if x.is_square() else None)(GF(p), GF(p)(a)),
            'sqrt_mod_p3mod4': lambda a, p: (lambda F, x: (lambda r: p - r if r > (p-1)//2 else r)(Integer(x.sqrt())) if x.is_square() else None)(GF(p), GF(p)(a)),
            'sqrt_mod_p5mod8': lambda a, p: (lambda F, x: (lambda r: p - r if r > (p-1)//2 else r)(Integer(x.sqrt())) if x.is_square() else None)(GF(p), GF(p)(a)),
            'sqrt_mod_general': lambda a, p: (lambda F, x: (lambda r: p - r if r > (p-1)//2 else r)(Integer(x.sqrt())) if x.is_square() else None)(GF(p), GF(p)(a)),
            'sqrt_mod_nonresidue': lambda a, p: (lambda F, x: (lambda r: p - r if r > (p-1)//2 else r)(Integer(x.sqrt())) if x.is_square() else None)(GF(p), GF(p)(a)),
            'primitive_root': primitive_root,
            'ff_multiplicative_generator': lambda p: Integer(GF(p).multiplicative_generator()),
            'discrete_log': lambda p, base, target: GF(p)(target).log(GF(p)(base)),
            'ff_ext_add': lambda p, n, a, b: (GF(p**n, 'a').fetch_int(Integer(a)) + GF(p**n, 'a').fetch_int(Integer(b))).integer_representation(),
            'ff_ext_mul': lambda p, n, a, b: (GF(p**n, 'a').fetch_int(Integer(a)) * GF(p**n, 'a').fetch_int(Integer(b))).integer_representation(),
            'ff_ext_inv': lambda p, n, a: (GF(p**n, 'a').fetch_int(Integer(a))**(-1)).integer_representation(),
            'ff_ext_pow': lambda p, n, a, e: (GF(p**n, 'a').fetch_int(Integer(a))**e).integer_representation(),
            'ff_ext_order': lambda p, n: p**n,
            'ff_ext_frobenius': lambda p, n, a: (GF(p**n, 'a').fetch_int(Integer(a))**p).integer_representation(),
            'ff_generator_order': lambda p, n: GF(p**n, 'a').multiplicative_generator().multiplicative_order(),
        },
        'polynomials': {
            'poly_add': poly_add,
            'poly_mul': poly_mul,
            'poly_quo_rem': poly_quo_rem,
            'poly_mod': poly_mod,
            'poly_gcd': poly_gcd,
            'poly_eval': poly_eval,
            'poly_factor': poly_factor,
            'poly_derivative': poly_derivative,
            'poly_is_irreducible': poly_is_irreducible,
            'poly_roots': poly_roots,
            'poly_pow': poly_pow,
        },
        'elliptic_curves': {
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
        },
        'matrix_ops': {
            # Determinant
            'determinant_2x2': matrix_determinant_2x2,
            'determinant_3x3': matrix_determinant_3x3,
            'determinant_4x4': matrix_determinant_4x4,
            # Rank
            'rank_2x3': matrix_rank_2x3,
            'rank_3x3': matrix_rank_3x3,
            # HNF
            'hnf_2x2': matrix_hnf_2x2,
            'hnf_3x3': matrix_hnf_3x3,
            # SNF
            'snf_2x2': matrix_snf_2x2,
            'snf_3x3': matrix_snf_3x3,
            # LLL
            'lll_2x2': matrix_lll_2x2,
            'lll_3x3': matrix_lll_3x3,
            # Elementary divisors
            'elementary_divisors_2x2': matrix_elementary_divisors_2x2,
            'elementary_divisors_3x3': matrix_elementary_divisors_3x3,
        },
        'arith_special': {
            'binomial': arith_binomial,
            'fibonacci': arith_fibonacci,
            'lucas_number': arith_lucas_number,
            'factorial': arith_factorial,
            'bernoulli_numerator': arith_bernoulli_numerator,
            'bernoulli_denominator': arith_bernoulli_denominator,
            'multinomial': arith_multinomial,
            'primitive_root': arith_primitive_root,
            'nth_prime': arith_nth_prime,
            'subfactorial': arith_subfactorial,
            'carmichael_lambda': arith_carmichael_lambda,
            'dedekind_psi': arith_dedekind_psi,
        },
        'polynomial_ops': {
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
        },
        'number_fields': {
            'quadratic_discriminant': nf_quadratic_discriminant,
            'nf_degree': nf_degree,
            'quadratic_signature': nf_quadratic_signature,
            'is_totally_real': nf_is_totally_real,
            'polynomial_discriminant': nf_polynomial_discriminant,
            'fundamental_discriminant_test': nf_fundamental_discriminant,
        },
    }

    if module not in function_map:
        raise ValueError(f"Unknown module: {module}")

    if function_name not in function_map[module]:
        raise ValueError(f"Unknown function: {module}.{function_name}")

    func = function_map[module][function_name]
    return func(*args)


def run_test_case(case, seed):
    """
    Run a single test case with the given seed.

    Args:
        case: Test case definition dict
        seed: Random seed to use

    Returns:
        Result dict with function, args, result, and seed
    """
    set_seed(seed)

    module = case.get('module', 'arith')
    function_name = case['function']
    arg_generators = case.get('argGenerators', [])

    # Generate arguments
    args = []
    for i, gen in enumerate(arg_generators):
        arg = generate_arg(gen, seed + i * 1000)
        args.append(arg)

    # Execute function
    try:
        result = execute_function(module, function_name, args)
        formatted_result = format_result(result)
        error = None
    except Exception as e:
        formatted_result = None
        error = str(e)

    return {
        'function': function_name,
        'args': [str(a) for a in args],
        'result': formatted_result,
        'error': error,
        'seed': seed,
    }


def run_test_suite(test_suite):
    """
    Run all test cases in a test suite.

    Args:
        test_suite: Dict with 'module' and 'cases' keys

    Returns:
        List of result dicts
    """
    module = test_suite.get('module', 'arith')
    cases = test_suite.get('cases', [])
    results = []

    for case in cases:
        case['module'] = module
        seeds = case.get('seeds', [42])

        for seed in seeds:
            result = run_test_case(case, seed)
            results.append(result)

    return results


def main():
    """Main entry point."""
    # Parse command line arguments
    seed_override = None
    for i, arg in enumerate(sys.argv[1:], 1):
        if arg == '--seed' and i < len(sys.argv) - 1:
            seed_override = int(sys.argv[i + 1])

    # Read test cases from stdin
    input_data = sys.stdin.read()

    if not input_data.strip():
        print("Error: No input provided. Pass test cases as JSON via stdin.", file=sys.stderr)
        sys.exit(1)

    try:
        test_suite = json.loads(input_data)
    except json.JSONDecodeError as e:
        print(f"Error parsing JSON: {e}", file=sys.stderr)
        sys.exit(1)

    # Override seeds if specified
    if seed_override is not None:
        for case in test_suite.get('cases', []):
            case['seeds'] = [seed_override]

    # Run tests
    results = run_test_suite(test_suite)

    # Output results as JSON
    print(json.dumps(results, indent=2))


if __name__ == '__main__':
    main()
