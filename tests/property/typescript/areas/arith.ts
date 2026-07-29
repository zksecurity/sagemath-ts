/**
 * sagemath-ts side of the `arith` property-test area.
 *
 * Cases: tests/property/cases/arith.cases.json
 * SageMath counterpart: tests/property/python/areas/arith.py
 */

import {
  arith,
  crt,
  divisors,
  euler_phi,
  factor,
  gcd,
  inverse_mod,
  is_prime,
  is_prime_power,
  is_square,
  is_squarefree,
  isqrt,
  jacobi_symbol,
  kronecker_symbol,
  lcm,
  legendre_symbol,
  moebius,
  next_prime,
  number_of_divisors,
  power_mod,
  previous_prime,
  prime_factors,
  prime_range,
  radical,
  sigma,
  squarefree_part,
  valuation,
  xgcd,
} from '../../../../packages/sagemath-ts/src/index.js';

export const functions = {
  gcd: (a: bigint, b: bigint) => gcd(a, b),
  lcm: (a: bigint, b: bigint) => lcm(a, b),
  xgcd: (a: bigint, b: bigint) => xgcd(a, b),
  factor: (n: bigint) => factor(n),
  is_prime: (n: bigint) => is_prime(n),
  is_prime_power: (n: bigint) => is_prime_power(n),
  next_prime: (n: bigint) => next_prime(n),
  previous_prime: (n: bigint) => previous_prime(n),
  euler_phi: (n: bigint) => euler_phi(n),
  radical: (n: bigint) => radical(n),
  moebius: (n: bigint) => moebius(n),
  kronecker_symbol: (a: bigint, n: bigint) => kronecker_symbol(a, n),
  legendre_symbol: (a: bigint, p: bigint) => legendre_symbol(a, p),
  jacobi_symbol: (a: bigint, n: bigint) => jacobi_symbol(a, n),
  power_mod: (a: bigint, n: bigint, m: bigint) => power_mod(a, n, m),
  inverse_mod: (a: bigint, m: bigint) => inverse_mod(a, m),
  crt: (a: bigint, b: bigint, m: bigint, n: bigint) => crt(a, b, m, n),
  isqrt: (n: bigint) => isqrt(n),
  is_square: (n: bigint) => is_square(n),
  is_squarefree: (n: bigint) => is_squarefree(n),
  divisors: (n: bigint) => divisors(n),
  number_of_divisors: (n: bigint) => number_of_divisors(n),
  sigma: (n: bigint, k?: bigint) => sigma(n, k),
  prime_range: (start: bigint, stop?: bigint) =>
    stop !== undefined ? prime_range(start, stop) : prime_range(start),
  trial_division: (n: bigint, bound?: bigint) =>
    bound !== undefined ? arith.trial_division(n, bound) : arith.trial_division(n),
  squarefree_part: (n: bigint) => squarefree_part(n),
  prime_factors: (n: bigint) => prime_factors(n),
  valuation: (n: bigint, p: bigint) => valuation(n, p),
};
