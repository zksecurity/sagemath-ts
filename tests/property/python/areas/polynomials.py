"""SageMath side of the ``polynomials`` property-test area (GF(p)[x]).

Cases: tests/property/cases/polynomials.cases.json
"""

from sage.all import *

from ._helpers import (
    poly_add,
    poly_derivative,
    poly_eval,
    poly_factor,
    poly_gcd,
    poly_is_irreducible,
    poly_mod,
    poly_mul,
    poly_pow,
    poly_quo_rem,
    poly_roots,
)

FUNCTIONS = {
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
}
