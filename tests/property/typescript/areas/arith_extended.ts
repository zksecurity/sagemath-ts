/**
 * sagemath-ts side of the `arith_extended` property-test area.
 *
 * Cases: tests/property/cases/arith_extended.cases.json
 * SageMath counterpart: tests/property/python/areas/arith_extended.py
 */

import { bernoulli, binomial, fibonacci } from '../../../../packages/sagemath-ts/src/arith/misc.js';

export const functions = {
  binomial: (n: bigint, k: bigint): bigint => binomial(n, k),
  fibonacci: (n: bigint): bigint => fibonacci(n),
  bernoulli_numerator: (n: bigint): bigint => bernoulli(n).numerator,
  bernoulli_denominator: (n: bigint): bigint => bernoulli(n).denominator,
};
