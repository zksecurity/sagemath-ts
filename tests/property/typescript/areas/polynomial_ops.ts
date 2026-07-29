/**
 * sagemath-ts side of the `polynomial_ops` property-test area (ZZ[x] and GF(p)[x]).
 *
 * Cases: tests/property/cases/polynomial_ops.cases.json
 * SageMath counterpart: tests/property/python/areas/polynomial_ops.py
 */

import { gcd } from '../../../../packages/sagemath-ts/src/index.js';
import {
  polyDerivative,
  polyEval,
  polyFactor,
  polyGcd,
  polyIsIrreducible,
  polyMul,
  polyPow,
  polyQuoRem,
  polyRoots,
} from './_helpers.js';

export const functions = {
  poly_gcd_zz: (coeffs1: bigint[], coeffs2: bigint[]) => {
    // Compute GCD of polynomials over ZZ using primitive polynomial algorithm
    // Helper to compute content (GCD of all coefficients)
    const content = (poly: bigint[]): bigint => {
      if (poly.length === 0) return 1n;
      let g = poly[0]!;
      for (let i = 1; i < poly.length; i++) {
        g = gcd(g, poly[i]!);
      }
      return g < 0n ? -g : g;
    };

    // Helper for polynomial pseudo-remainder
    const pseudoRemainder = (f: bigint[], g: bigint[]): bigint[] => {
      if (g.length === 0) throw new Error('division by zero');
      const m = f.length - 1;
      const n = g.length - 1;
      if (m < n) return [...f];

      const lc_g = g[n]!;
      const r = [...f];
      for (let i = m - n; i >= 0; i--) {
        if (r.length - 1 >= i + n) {
          const q = r[i + n]!;
          for (let j = 0; j <= n; j++) {
            r[i + j] = lc_g * r[i + j]! - q * g[j]!;
          }
        }
      }
      // Remove leading zeros
      while (r.length > 0 && r[r.length - 1] === 0n) {
        r.pop();
      }
      return r;
    };

    // Make polynomials primitive
    const c1 = content(coeffs1);
    const c2 = content(coeffs2);
    let f = coeffs1.map((c) => c / c1);
    let g = coeffs2.map((c) => c / c2);

    // Ensure f.degree >= g.degree
    if (f.length < g.length) {
      [f, g] = [g, f];
    }

    // Subresultant GCD algorithm simplified
    while (g.length > 0) {
      const r = pseudoRemainder(f, g);
      if (r.length === 0) break;
      const c = content(r);
      f = g;
      g = r.map((x) => x / c);
    }

    // Make primitive and positive leading coeff
    const result = g.length > 0 ? g : f;
    const cr = content(result);
    let finalResult = result.map((x) => x / cr);
    if (finalResult.length > 0 && finalResult[finalResult.length - 1]! < 0n) {
      finalResult = finalResult.map((x) => -x);
    }
    return finalResult.length > 0 ? finalResult : [1n];
  },
  poly_eval_zz: (coeffs: bigint[], x: bigint) => {
    // Evaluate polynomial over ZZ at a point using Horner's method
    if (coeffs.length === 0) return 0n;
    let result = coeffs[coeffs.length - 1]!;
    for (let i = coeffs.length - 2; i >= 0; i--) {
      result = result * x + coeffs[i]!;
    }
    return result;
  },
  poly_gcd_ff: (p: bigint, coeffs1: bigint[], coeffs2: bigint[]) => polyGcd(p, coeffs1, coeffs2),
  poly_eval_ff: (p: bigint, coeffs: bigint[], x: bigint) => polyEval(p, coeffs, x),
  poly_factor_ff: (p: bigint, coeffs: bigint[]) => polyFactor(p, coeffs),
  poly_roots_ff: (p: bigint, coeffs: bigint[]) => polyRoots(p, coeffs),
  poly_derivative_ff: (p: bigint, coeffs: bigint[]) => polyDerivative(p, coeffs),
  poly_is_irreducible_ff: (p: bigint, coeffs: bigint[]) => polyIsIrreducible(p, coeffs),
  poly_mul_ff: (p: bigint, coeffs1: bigint[], coeffs2: bigint[]) => polyMul(p, coeffs1, coeffs2),
  poly_quo_rem_ff: (p: bigint, coeffs1: bigint[], coeffs2: bigint[]) =>
    polyQuoRem(p, coeffs1, coeffs2),
  poly_pow_ff: (p: bigint, coeffs: bigint[], n: bigint) => polyPow(p, coeffs, n),
};
