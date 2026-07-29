"""SageMath side of the ``arith_extended`` property-test area.

Cases: tests/property/cases/arith_extended.cases.json
"""

from sage.all import bernoulli, binomial, fibonacci


def bernoulli_numerator(n):
    """Return the exact numerator of ``B_n``."""
    return bernoulli(n).numerator()


def bernoulli_denominator(n):
    """Return the positive exact denominator of ``B_n``."""
    return bernoulli(n).denominator()


FUNCTIONS = {
    'binomial': binomial,
    'fibonacci': fibonacci,
    'bernoulli_numerator': bernoulli_numerator,
    'bernoulli_denominator': bernoulli_denominator,
}
