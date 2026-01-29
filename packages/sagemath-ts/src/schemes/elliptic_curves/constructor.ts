/**
 * @module sage/schemes/elliptic_curves/constructor
 * @description Elliptic curve constructor function
 *
 * Port of: sage/schemes/elliptic_curves/constructor.py
 *
 * Provides the main `EllipticCurve()` factory function for creating elliptic curves
 * over various base rings/fields.
 */

import { ArithmeticError, ValueError } from '../../errors.js';
import { EllipticCurveGeneric } from './ell_generic.js';
import type { FieldElement, FieldParent } from './ell_point.js';

/**
 * Construct an elliptic curve.
 *
 * An elliptic curve is always specified by the coefficients of a
 * Weierstrass equation:
 *   y^2 + a1*x*y + a3*y = x^3 + a2*x^2 + a4*x + a6
 *
 * Input formats:
 * - `EllipticCurve(K, [a1, a2, a3, a4, a6])`: Full Weierstrass over field K
 * - `EllipticCurve(K, [a, b])`: Short form y^2 = x^3 + ax + b over K
 *
 * @example
 * ```typescript
 * import { GF } from '../rings/finite_rings/finite_field_constructor.js';
 * import { EllipticCurve } from './constructor.js';
 *
 * // Create a prime field
 * const F = GF(23n);
 *
 * // Short Weierstrass: y^2 = x^3 + x + 1
 * const E = EllipticCurve(F, [1n, 1n]);
 *
 * // Full Weierstrass: y^2 + xy + y = x^3 + 2x^2 + 3x + 4
 * const E2 = EllipticCurve(F, [1n, 2n, 1n, 3n, 4n]);
 *
 * // Check properties
 * console.log(E.discriminant());  // Non-zero for valid curve
 * console.log(E.j_invariant());   // The j-invariant
 *
 * // Create points
 * const P = E.point([F.__call__(0n), F.__call__(1n)]);
 * ```
 *
 * @param K - The base field
 * @param coeffs - Weierstrass coefficients: [a, b] for short form or [a1, a2, a3, a4, a6] for long form
 * @returns An elliptic curve over K
 */
export function EllipticCurve<F extends FieldElement>(
  K: FieldParent,
  coeffs:
    | [bigint | number, bigint | number]
    | [bigint | number, bigint | number, bigint | number, bigint | number, bigint | number]
): EllipticCurveGeneric<F>;

export function EllipticCurve<F extends FieldElement>(
  K: FieldParent,
  coeffs: [F, F] | [F, F, F, F, F]
): EllipticCurveGeneric<F>;

export function EllipticCurve<F extends FieldElement>(
  K: FieldParent,
  coeffs: unknown[]
): EllipticCurveGeneric<F> {
  // Convert coefficients to field elements if needed
  const fieldCoeffs = coeffs.map((c) => {
    if (typeof c === 'number' || typeof c === 'bigint') {
      return K.__call__(c) as F;
    }
    return c as F;
  });

  // Determine which form we have
  if (fieldCoeffs.length === 2) {
    // Short Weierstrass form: y^2 = x^3 + ax + b
    // Corresponds to a1 = a2 = a3 = 0, a4 = a, a6 = b
    const [a, b] = fieldCoeffs;
    const zero = K.zero() as F;
    const ainvs: [F, F, F, F, F] = [zero, zero, zero, a, b];
    return new EllipticCurveGeneric(K, ainvs);
  } else if (fieldCoeffs.length === 5) {
    // Long Weierstrass form: y^2 + a1*x*y + a3*y = x^3 + a2*x^2 + a4*x + a6
    const ainvs = fieldCoeffs as [F, F, F, F, F];
    return new EllipticCurveGeneric(K, ainvs);
  } else {
    throw new ValueError(
      `Invalid number of coefficients: ${coeffs.length}. ` +
        'Expected 2 (short form [a, b]) or 5 (long form [a1, a2, a3, a4, a6]).'
    );
  }
}

/**
 * Create an elliptic curve from a j-invariant.
 *
 * Returns an elliptic curve with the given j-invariant.
 *
 * Special cases:
 * - j = 0: Returns y^2 = x^3 + 1
 * - j = 1728: Returns y^2 = x^3 + x
 * - Otherwise: Returns y^2 = x^3 - 3jk*x - 2jk^2 where k = j - 1728
 *
 * @param K - The base field
 * @param j - The j-invariant
 * @returns An elliptic curve with j-invariant j
 *
 * @example
 * ```typescript
 * const F = GF(23n);
 * const E = EllipticCurve_from_j(F, 0n);  // j = 0
 * console.log(E.j_invariant());  // 0
 * ```
 */
export function EllipticCurve_from_j<F extends FieldElement>(
  K: FieldParent,
  j: bigint | number | F
): EllipticCurveGeneric<F> {
  const jElem = typeof j === 'number' || typeof j === 'bigint' ? (K.__call__(j) as F) : j;

  const zero = K.zero() as F;
  const one = K.one() as F;

  // Special case: j = 0
  if (jElem.isZero()) {
    // y^2 = x^3 + 1
    return EllipticCurve(K, [zero, one]);
  }

  // Special case: j = 1728
  const n1728 = K.__call__(1728) as F;
  if (jElem.eq(n1728)) {
    // y^2 = x^3 + x
    return EllipticCurve(K, [one, zero]);
  }

  // General case: y^2 = x^3 - 3jk*x - 2jk^2 where k = j - 1728
  const k = jElem.sub(n1728) as F;
  const n3 = K.__call__(3) as F;
  const n2 = K.__call__(2) as F;

  const a = jElem.mul(k).mul(n3).neg() as F; // -3jk
  const b = jElem.mul(k).mul(k).mul(n2).neg() as F; // -2jk^2

  return EllipticCurve(K, [a, b]);
}

/**
 * Create an elliptic curve from c4 and c6 invariants.
 *
 * Returns an elliptic curve with the given c-invariants.
 * The curve has the form: y^2 = x^3 - (c4/48)*x - (c6/864)
 *
 * @param K - The base field
 * @param c4 - The c4 invariant
 * @param c6 - The c6 invariant
 * @returns An elliptic curve with c-invariants c4, c6
 */
export function EllipticCurve_from_c4c6<F extends FieldElement>(
  K: FieldParent,
  c4: bigint | number | F,
  c6: bigint | number | F
): EllipticCurveGeneric<F> {
  const c4Elem = typeof c4 === 'number' || typeof c4 === 'bigint' ? (K.__call__(c4) as F) : c4;

  const c6Elem = typeof c6 === 'number' || typeof c6 === 'bigint' ? (K.__call__(c6) as F) : c6;

  const n48 = K.__call__(48) as F;
  const n864 = K.__call__(864) as F;

  // y^2 = x^3 - (c4/48)*x - (c6/864)
  const a = c4Elem.div(n48).neg() as F;
  const b = c6Elem.div(n864).neg() as F;

  return EllipticCurve(K, [a, b]);
}

// Re-export for convenience
export { EllipticCurveGeneric } from './ell_generic.js';
export { EllipticCurvePoint, type FieldElement, type FieldParent } from './ell_point.js';
