"""SageMath side of the ``small_roots`` property-test area.

Covers Coppersmith's ``small_roots``
(``sage/rings/polynomial/polynomial_modn_dense_ntl.pyx``), reached through the
``small_roots`` method of ``Zmod(N)[x]`` and ``GF(p)[x]`` polynomials, and the
``Integer ** RealNumber`` power that sets its default bound ``X``.

Cases: tests/property/cases/small_roots.cases.json
TypeScript counterpart: tests/property/typescript/areas/small_roots.ts

The shared generators produce integers only, so ``beta`` and ``epsilon`` are
indices into the tables below. KEEP THE TABLES IN SYNC with the TypeScript
counterpart. Roots are printed as their lifts, in the order Sage returns them
(the iteration order of the Python set ``small_roots`` builds).
"""

from sage.all import GF, RR, ZZ, Integer, PolynomialRing, Zmod, next_prime

#: ``beta`` values. The last two are rejected by ``small_roots``.
BETAS = [1.0, 0.5, 0.25, 0.75, 0.4, 0.9, 0.0, 1.5]

#: ``epsilon`` values; ``None`` selects the default ``beta / 8``.
EPSILONS = [None, 0.1, 0.05, -0.1, 0.0, 0.02]


def _fmt(roots):
    return '[' + ', '.join(str(r) for r in roots) + ']'


def _options(x_mode, x, beta_index, epsilon_index):
    options = {'beta': RR(BETAS[int(beta_index)])}
    if x_mode:
        options['X'] = Integer(x)
    epsilon = EPSILONS[int(epsilon_index)]
    if epsilon is not None:
        options['epsilon'] = RR(epsilon)
    return options


def _planted(R, roots, c, k):
    """``prod(x - r) * (x^k + sum_{i<k} c^(i+1) x^i)`` over ``R``."""
    x = R.gen()
    f = R.one()
    for r in roots:
        f *= x - r
    g = x ** int(k) + sum(R.base_ring()(c) ** (i + 1) * x ** i for i in range(int(k)))
    return f * g


def small_roots_monic(N, lower, x_mode, x, beta_index, epsilon_index):
    """``x^n + lower[n-1] x^(n-1) + ... + lower[0]`` over ``Zmod(N)``."""
    R = PolynomialRing(Zmod(N), 'x')
    f = R(list(lower) + [1])
    return _fmt(f.small_roots(**_options(x_mode, x, beta_index, epsilon_index)))


def small_roots_coeffs(N, coeffs, x_mode, x, beta_index, epsilon_index):
    R = PolynomialRing(Zmod(N), 'x')
    f = R(list(coeffs))
    return _fmt(f.small_roots(**_options(x_mode, x, beta_index, epsilon_index)))


def small_roots_planted(N, roots, c, k, beta_index, epsilon_index):
    R = PolynomialRing(Zmod(N), 'x')
    f = _planted(R, roots, c, k)
    return _fmt(f.small_roots(**_options(0, 0, beta_index, epsilon_index)))


def small_roots_field(a, roots, c, k):
    """Planted roots over the prime field ``GF(next_prime(a))``."""
    R = PolynomialRing(GF(next_prime(a)), 'x')
    return _fmt(_planted(R, roots, c, k).small_roots())


def small_roots_stereotyped(a, b, kbits, K):
    """RSA, e = 3, a ``kbits``-bit key ``K`` padded with 1s to 512 bits."""
    N = next_prime(a) * next_prime(b)
    ZmodN = Zmod(N)
    padding = ZZ(2) ** 512 - ZZ(2) ** int(kbits)
    C = ZmodN(padding + K) ** 3
    R = PolynomialRing(ZmodN, 'x')
    x = R.gen()
    f = (padding + x) ** 3 - C
    return _fmt(f.small_roots())


def small_roots_factor(a, b, hidden, noise):
    """Factor ``N = p q`` from ``q + noise`` with ``noise < 2^hidden``, beta = 0.5."""
    p = next_prime(a)
    q = next_prime(b)
    R = PolynomialRing(Zmod(p * q), 'x')
    f = R.gen() - (q + noise)
    return _fmt(f.small_roots(X=ZZ(2) ** int(hidden) - 1, beta=0.5))


def small_roots_method_error(kind):
    """Polynomial classes Sage does not give a ``small_roots`` method."""
    base = [GF(2), Zmod(2), GF(4, 'a'), GF(9, 'a')][int(kind)]
    return _fmt(PolynomialRing(base, 'x').gen().small_roots())


def integer_pow_rr(N, y_num, y_den):
    """``Integer ** RealNumber``, as the exact rational value of the RR result."""
    z = ZZ(N) ** (RR(y_num) / RR(y_den))
    return str(z.exact_rational())


FUNCTIONS = {
    'small_roots_monic': small_roots_monic,
    'small_roots_coeffs': small_roots_coeffs,
    'small_roots_planted': small_roots_planted,
    'small_roots_field': small_roots_field,
    'small_roots_stereotyped': small_roots_stereotyped,
    'small_roots_factor': small_roots_factor,
    'small_roots_method_error': small_roots_method_error,
    'integer_pow_rr': integer_pow_rr,
}
