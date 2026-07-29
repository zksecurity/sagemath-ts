"""SageMath side of the ``arith_special`` property-test area.

Cases: tests/property/cases/arith_special.cases.json
"""

from sage.all import *
from sage.arith.misc import dedekind_psi


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


FUNCTIONS = {
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
}
