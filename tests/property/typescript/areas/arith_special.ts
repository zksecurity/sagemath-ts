/**
 * sagemath-ts side of the `arith_special` property-test area.
 *
 * Cases: tests/property/cases/arith_special.cases.json
 * SageMath counterpart: tests/property/python/areas/arith_special.py
 */

import {
  bernoulli,
  binomial,
  carmichael_lambda,
  dedekind_psi,
  factorial,
  fibonacci,
  lucas_number,
  multinomial,
  nth_prime,
  primitive_root,
  subfactorial,
} from '../../../../packages/sagemath-ts/src/arith/misc.js';

export const functions = {
  binomial: (n: bigint, k: bigint) => binomial(n, k),
  fibonacci: (n: bigint) => fibonacci(n),
  lucas_number: (n: bigint) => lucas_number(n),
  factorial: (n: bigint) => factorial(n),
  bernoulli_numerator: (n: bigint) => {
    const b = bernoulli(n);
    return b.numerator;
  },
  bernoulli_denominator: (n: bigint) => {
    const b = bernoulli(n);
    return b.denominator;
  },
  multinomial: (...args: bigint[]) => multinomial(...args),
  primitive_root: (n: bigint) => primitive_root(n),
  nth_prime: (n: bigint) => nth_prime(n),
  subfactorial: (n: bigint) => subfactorial(n),
  carmichael_lambda: (n: bigint) => carmichael_lambda(n),
  dedekind_psi: (n: bigint) => dedekind_psi(n),
};
