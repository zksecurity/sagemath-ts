"""SageMath side of the ``arith`` property-test area.

Cases: tests/property/cases/arith.cases.json
"""

from sage.all import *

FUNCTIONS = {
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
}
