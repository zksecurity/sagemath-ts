/**
 * sagemath-ts side of the `polynomials` property-test area (GF(p)[x]).
 *
 * Cases: tests/property/cases/polynomials.cases.json
 * SageMath counterpart: tests/property/python/areas/polynomials.py
 */

import {
  polyAdd,
  polyDerivative,
  polyEval,
  polyFactor,
  polyGcd,
  polyIsIrreducible,
  polyMod,
  polyMul,
  polyPow,
  polyQuoRem,
  polyRoots,
} from './_helpers.js';

export const functions = {
  poly_add: (p: bigint, coeffs1: bigint[], coeffs2: bigint[]) => polyAdd(p, coeffs1, coeffs2),
  poly_mul: (p: bigint, coeffs1: bigint[], coeffs2: bigint[]) => polyMul(p, coeffs1, coeffs2),
  poly_quo_rem: (p: bigint, coeffs1: bigint[], coeffs2: bigint[]) =>
    polyQuoRem(p, coeffs1, coeffs2),
  poly_mod: (p: bigint, coeffs1: bigint[], coeffs2: bigint[]) => polyMod(p, coeffs1, coeffs2),
  poly_gcd: (p: bigint, coeffs1: bigint[], coeffs2: bigint[]) => polyGcd(p, coeffs1, coeffs2),
  poly_eval: (p: bigint, coeffs: bigint[], x: bigint) => polyEval(p, coeffs, x),
  poly_factor: (p: bigint, coeffs: bigint[]) => polyFactor(p, coeffs),
  poly_derivative: (p: bigint, coeffs: bigint[]) => polyDerivative(p, coeffs),
  poly_is_irreducible: (p: bigint, coeffs: bigint[]) => polyIsIrreducible(p, coeffs),
  poly_roots: (p: bigint, coeffs: bigint[]) => polyRoots(p, coeffs),
  poly_pow: (p: bigint, coeffs: bigint[], n: bigint) => polyPow(p, coeffs, n),
};
