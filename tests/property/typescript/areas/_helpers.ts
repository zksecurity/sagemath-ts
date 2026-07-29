/**
 * Shared helpers used by more than one area module.
 *
 * Only put something here when two or more areas genuinely need it (today: the
 * GF(p) polynomial helpers shared by `polynomials` and `polynomial_ops`).
 * Area-private helpers belong in the area module itself so that agents adding
 * areas in parallel never touch the same file.
 *
 * Mirrors tests/property/python/areas/_helpers.py.
 *
 * Files whose name starts with `_` are skipped by the runner's area discovery.
 */

import { GF } from '../../../../packages/sagemath-ts/src/index.js';
import type { FiniteFieldElement } from '../../../../packages/sagemath-ts/src/rings/finite_rings/index.js';
import {
  type Polynomial,
  PolynomialRingConstructor,
} from '../../../../packages/sagemath-ts/src/rings/polynomial/index.js';

/**
 * Create a polynomial over GF(p) from coefficient list [c0, c1, ..., cn].
 */
export function makePoly(p: bigint, coeffs: bigint[]): Polynomial<FiniteFieldElement> {
  const F = GF(p);
  const [R, _x] = PolynomialRingConstructor(F, 'x');
  const fieldCoeffs = coeffs.map((c) => F.__call__(c));
  return R.__call__(fieldCoeffs);
}

/**
 * Convert polynomial to coefficient list.
 */
export function polyToList(poly: Polynomial<FiniteFieldElement>): bigint[] {
  if (poly.isZero()) {
    return [];
  }
  return poly.coeffs.map((c) => c.value);
}

/**
 * Add two polynomials over GF(p).
 */
export function polyAdd(p: bigint, coeffs1: bigint[], coeffs2: bigint[]): bigint[] {
  const f = makePoly(p, coeffs1);
  const g = makePoly(p, coeffs2);
  return polyToList(f.add(g));
}

/**
 * Multiply two polynomials over GF(p).
 */
export function polyMul(p: bigint, coeffs1: bigint[], coeffs2: bigint[]): bigint[] {
  const f = makePoly(p, coeffs1);
  const g = makePoly(p, coeffs2);
  return polyToList(f.mul(g));
}

/**
 * Compute quotient and remainder of polynomial division over GF(p).
 */
export function polyQuoRem(p: bigint, coeffs1: bigint[], coeffs2: bigint[]): [bigint[], bigint[]] {
  const f = makePoly(p, coeffs1);
  const g = makePoly(p, coeffs2);
  const [q, r] = f.quo_rem(g);
  return [polyToList(q), polyToList(r)];
}

/**
 * Compute polynomial modulo over GF(p).
 */
export function polyMod(p: bigint, coeffs1: bigint[], coeffs2: bigint[]): bigint[] {
  const f = makePoly(p, coeffs1);
  const g = makePoly(p, coeffs2);
  return polyToList(f.mod(g));
}

/**
 * Compute GCD of two polynomials over GF(p).
 */
export function polyGcd(p: bigint, coeffs1: bigint[], coeffs2: bigint[]): bigint[] {
  const f = makePoly(p, coeffs1);
  const g = makePoly(p, coeffs2);
  const result = f.gcd(g);
  // gcd() already returns monic
  return polyToList(result);
}

/**
 * Evaluate polynomial at a point over GF(p).
 */
export function polyEval(p: bigint, coeffs: bigint[], x: bigint): bigint {
  const f = makePoly(p, coeffs);
  const F = GF(p);
  const xElem = F.__call__(x);
  return f.evaluate(xElem).value;
}

/**
 * Factor polynomial over GF(p). Returns list of (factor_coeffs, multiplicity).
 */
export function polyFactor(p: bigint, coeffs: bigint[]): Array<[bigint[], number]> {
  const f = makePoly(p, coeffs);
  const factorization = f.factor();
  const result: Array<[bigint[], number]> = [];

  for (const [fac, mult] of factorization) {
    // factor() already returns monic factors
    const facCoeffs = polyToList(fac);
    result.push([facCoeffs, mult]);
  }

  // Sort by degree then lexicographically (by coefficient list)
  result.sort((a, b) => {
    if (a[0].length !== b[0].length) {
      return a[0].length - b[0].length;
    }
    // Compare coefficient by coefficient
    for (let i = 0; i < a[0].length; i++) {
      if (a[0][i]! !== b[0][i]!) {
        return a[0][i]! < b[0][i]! ? -1 : 1;
      }
    }
    return 0;
  });

  return result;
}

/**
 * Compute derivative of polynomial over GF(p).
 */
export function polyDerivative(p: bigint, coeffs: bigint[]): bigint[] {
  const f = makePoly(p, coeffs);
  return polyToList(f.derivative());
}

/**
 * Test if polynomial is irreducible over GF(p).
 */
export function polyIsIrreducible(p: bigint, coeffs: bigint[]): boolean {
  const f = makePoly(p, coeffs);
  return f.is_irreducible();
}

/**
 * Find roots of polynomial over GF(p). Returns sorted list of roots.
 */
export function polyRoots(p: bigint, coeffs: bigint[]): bigint[] {
  const f = makePoly(p, coeffs);
  const roots = f.roots();
  // Extract just the roots (not multiplicities) and sort
  return roots.map(([r, _mult]) => r.value).sort((a, b) => (a < b ? -1 : a > b ? 1 : 0));
}

/**
 * Compute polynomial power over GF(p).
 */
export function polyPow(p: bigint, coeffs: bigint[], n: bigint): bigint[] {
  const f = makePoly(p, coeffs);
  return polyToList(f.pow(n));
}
