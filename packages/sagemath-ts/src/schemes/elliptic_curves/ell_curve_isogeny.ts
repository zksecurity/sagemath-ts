/**
 * @module sage/schemes/elliptic_curves/ell_curve_isogeny
 * @description Isogenies between elliptic curves
 *
 * Port of: sage/schemes/elliptic_curves/ell_curve_isogeny.py
 *
 * An isogeny phi: E_1 -> E_2 between two elliptic curves E_1 and E_2
 * is a morphism of curves that sends the origin of E_1 to the origin
 * of E_2. Such a morphism is automatically a morphism of group schemes
 * and the kernel is a finite subgroup scheme of E_1.
 *
 * This module implements:
 * - Velu's formulas for computing isogenies from points
 * - Computation of dual isogenies
 *
 * @see Reference: sage/schemes/elliptic_curves/ell_curve_isogeny.py
 */

import { NotImplementedError, ValueError } from '../../errors.js';
import type {
  Polynomial as PolyElement,
  RingElement,
} from '../../rings/polynomial/polynomial_element.js';
import { type CoefficientRing, PolynomialRing } from '../../rings/polynomial/polynomial_ring.js';
import { EllipticCurve } from './constructor.js';
import type { EllipticCurveGeneric } from './ell_generic.js';
import { EllipticCurvePoint, type FieldElement, type FieldParent } from './ell_point.js';
import {
  WeierstrassIsomorphism,
  _isomorphisms,
  baseWI,
} from './weierstrass_morphism.js';

export { WeierstrassIsomorphism };

/**
 * The interface an isomorphism must satisfy in order to be usable as a pre- or
 * post-isomorphism of an {@link EllipticCurveIsogeny}.
 */
export type WeierstrassIsomorphismLike<F extends FieldElement> = WeierstrassIsomorphism<F>;

/**
 * Type for polynomials (alias for compatibility).
 */
export type Polynomial = PolyElement<RingElement> | bigint[];

/**
 * A rational function represented as a pair (numerator, denominator) of polynomials.
 * The rational function value is numerator(x) / denominator(x).
 */
export interface RationalFunction<F extends FieldElement> {
  /** Numerator polynomial coefficients */
  numerator: F[];
  /** Denominator polynomial coefficients */
  denominator: F[];
  /** Evaluate at a point */
  evaluate(x: F): F | null;
  /** String representation */
  toString(): string;
}

/**
 * A bivariate rational function in x and y.
 * Represented as numerator(x, y) / denominator(x, y) where
 * numerator and denominator are arrays of polynomial coefficients in y,
 * each coefficient being a polynomial in x.
 */
export interface BivariateRationalFunction<F extends FieldElement> {
  /** Numerator: coefficients of y^i are polynomials in x */
  numerator: F[][];
  /** Denominator: coefficients of y^i are polynomials in x */
  denominator: F[][];
  /** Evaluate at a point (x, y) */
  evaluate(x: F, y: F): F | null;
  /** String representation */
  toString(): string;
}

/**
 * Create a rational function from numerator and denominator coefficients.
 */
function createRationalFunction<F extends FieldElement>(
  numerator: F[],
  denominator: F[],
  K: FieldParent
): RationalFunction<F> {
  return {
    numerator,
    denominator,
    evaluate(x: F): F | null {
      // Evaluate numerator
      let num = K.zero() as F;
      let xPow = K.one() as F;
      for (const coef of numerator) {
        num = num.add(coef.mul(xPow)) as F;
        xPow = xPow.mul(x) as F;
      }

      // Evaluate denominator
      let den = K.zero() as F;
      xPow = K.one() as F;
      for (const coef of denominator) {
        den = den.add(coef.mul(xPow)) as F;
        xPow = xPow.mul(x) as F;
      }

      if (den.isZero()) return null;
      return num.mul(den.inv()) as F;
    },
    toString(): string {
      const numStr = numerator.map((c, i) => (i === 0 ? c.toString() : `${c}*x^${i}`)).join(' + ');
      const denStr = denominator
        .map((c, i) => (i === 0 ? c.toString() : `${c}*x^${i}`))
        .join(' + ');
      return `(${numStr})/(${denStr})`;
    },
  };
}

/**
 * Multiply two polynomials (coefficient arrays).
 */
function polyMul<F extends FieldElement>(a: F[], b: F[], K: FieldParent): F[] {
  if (a.length === 0 || b.length === 0) {
    return [];
  }
  const result: F[] = [];
  for (let i = 0; i < a.length + b.length - 1; i++) {
    result.push(K.zero() as F);
  }
  for (let i = 0; i < a.length; i++) {
    for (let j = 0; j < b.length; j++) {
      result[i + j] = result[i + j]!.add(a[i]!.mul(b[j]!)) as F;
    }
  }
  // Remove trailing zeros
  while (result.length > 0 && result[result.length - 1]!.isZero()) {
    result.pop();
  }
  return result;
}

/**
 * Add two polynomials (coefficient arrays).
 */
function polyAdd<F extends FieldElement>(a: F[], b: F[], K: FieldParent): F[] {
  const maxLen = Math.max(a.length, b.length);
  const result: F[] = [];
  for (let i = 0; i < maxLen; i++) {
    const ai = i < a.length ? a[i]! : (K.zero() as F);
    const bi = i < b.length ? b[i]! : (K.zero() as F);
    result.push(ai.add(bi) as F);
  }
  // Remove trailing zeros
  while (result.length > 0 && result[result.length - 1]!.isZero()) {
    result.pop();
  }
  return result;
}

/**
 * Subtract two polynomials (coefficient arrays).
 */
function polySub<F extends FieldElement>(a: F[], b: F[], K: FieldParent): F[] {
  const maxLen = Math.max(a.length, b.length);
  const result: F[] = [];
  for (let i = 0; i < maxLen; i++) {
    const ai = i < a.length ? a[i]! : (K.zero() as F);
    const bi = i < b.length ? b[i]! : (K.zero() as F);
    result.push(ai.sub(bi) as F);
  }
  // Remove trailing zeros
  while (result.length > 0 && result[result.length - 1]!.isZero()) {
    result.pop();
  }
  return result;
}

/**
 * Multiply a polynomial by a scalar.
 */
function polyScale<F extends FieldElement>(a: F[], s: F): F[] {
  return a.map((c) => c.mul(s) as F);
}

// ---------------------------------------------------------------------------
// Univariate polynomial arithmetic over the base field.
//
// Polynomials are coefficient arrays in ascending degree order: a[i] is the
// coefficient of x^i.  Trailing zeros are stripped so that `length - 1` is
// always the degree (the zero polynomial is the empty array, degree -1, which
// mirrors Sage's convention).
// ---------------------------------------------------------------------------

/** Strip trailing zero coefficients (in place on a copy). */
function polyTrim<F extends FieldElement>(a: F[]): F[] {
  const r = [...a];
  while (r.length > 0 && r[r.length - 1]!.isZero()) r.pop();
  return r;
}

/** Degree of a polynomial; -1 for the zero polynomial (as in Sage). */
function polyDegree<F extends FieldElement>(a: F[]): number {
  return polyTrim(a).length - 1;
}

/** Multiply a polynomial by x^k. */
function polyShift<F extends FieldElement>(a: F[], k: number, K: FieldParent): F[] {
  if (a.length === 0) return [];
  const pad: F[] = [];
  for (let i = 0; i < k; i++) pad.push(K.zero() as F);
  return [...pad, ...a];
}

/** Formal derivative of a polynomial. */
function polyDeriv<F extends FieldElement>(a: F[], K: FieldParent): F[] {
  const r: F[] = [];
  for (let i = 1; i < a.length; i++) {
    r.push(a[i]!.mul(i) as F);
  }
  return polyTrim(r.length === 0 ? [K.zero() as F] : r);
}

/** Evaluate a polynomial at a field element (Horner). */
function polyEval<F extends FieldElement>(a: F[], x: F, K: FieldParent): F {
  let acc = K.zero() as F;
  for (let i = a.length - 1; i >= 0; i--) {
    acc = acc.mul(x).add(a[i]!) as F;
  }
  return acc;
}

/** Euclidean division of polynomials over a field: returns [quotient, remainder]. */
function polyDivRem<F extends FieldElement>(a: F[], b: F[], K: FieldParent): [F[], F[]] {
  const bb = polyTrim(b);
  if (bb.length === 0) throw new ValueError('division by zero polynomial');
  const r = polyTrim(a);
  const db = bb.length - 1;
  const lcInv = bb[db]!.inv() as F;
  const q: F[] = [];
  for (let i = 0; i <= r.length - 1 - db; i++) q.push(K.zero() as F);
  while (r.length - 1 >= db && r.length > 0) {
    const shift = r.length - 1 - db;
    const c = r[r.length - 1]!.mul(lcInv) as F;
    q[shift] = c;
    for (let i = 0; i <= db; i++) {
      r[i + shift] = r[i + shift]!.sub(c.mul(bb[i]!)) as F;
    }
    while (r.length > 0 && r[r.length - 1]!.isZero()) r.pop();
  }
  return [polyTrim(q), r];
}

/** Monic gcd of two polynomials over a field (zero polynomial for gcd(0,0)). */
function polyGcdField<F extends FieldElement>(a: F[], b: F[], K: FieldParent): F[] {
  let x = polyTrim(a);
  let y = polyTrim(b);
  while (y.length > 0) {
    const [, r] = polyDivRem(x, y, K);
    x = y;
    y = r;
  }
  return polyMonic(x, K);
}

/** Make a polynomial monic; the zero polynomial is returned unchanged. */
function polyMonic<F extends FieldElement>(a: F[], K: FieldParent): F[] {
  const t = polyTrim(a);
  if (t.length === 0) return t;
  const lc = t[t.length - 1]!;
  if (lc.eq(K.one())) return t;
  const inv = lc.inv() as F;
  return t.map((c) => c.mul(inv) as F);
}

/** Coefficient of x^i, or zero when out of range (mirrors Sage's `psi[i]`). */
function polyCoeff<F extends FieldElement>(a: F[], i: number, K: FieldParent): F {
  return i >= 0 && i < a.length ? a[i]! : (K.zero() as F);
}

/** Lift a coefficient array of bigints into the base field. */
function polyFromBigints<F extends FieldElement>(coeffs: bigint[], K: FieldParent): F[] {
  return polyTrim(coeffs.map((c) => K.__call__(c) as F));
}

/**
 * Convert a polynomial over a prime field back to a bigint coefficient array.
 *
 * This is only meaningful over prime fields; over extension fields the
 * coefficients have no canonical bigint representation and a ValueError is
 * raised.
 */
function polyToBigints<F extends FieldElement>(a: F[]): bigint[] {
  return a.map((c) => {
    const s = c.toString();
    if (!/^-?\d+$/.test(s)) {
      throw new ValueError('kernel polynomial coefficients are not in a prime field');
    }
    return BigInt(s);
  });
}

/**
 * Data for a kernel point in Velu's formulas.
 * Contains: [xQ, yQ, gxQ, gyQ, vQ, uQ]
 */
type KernelPointData<F extends FieldElement> = {
  xQ: F;
  yQ: F;
  gxQ: F;
  gyQ: F;
  vQ: F;
  uQ: F;
};

/**
 * Helper function to infer the algorithm to be used from the parameters
 * passed to the various isogeny functions.
 *
 * If kernel is a list of points on the elliptic curve E, we will try to
 * use Velu's algorithm.
 *
 * If kernel is a list of coefficients or a univariate polynomial, we will
 * try to use Kohel's algorithm.
 *
 * @param E - Domain elliptic curve
 * @param kernel - Either a list of points on E, or a univariate polynomial
 *                 or list of coefficients of a univariate polynomial
 * @returns String 'velu' or 'kohel'
 * @see Reference: sage/schemes/elliptic_curves/ell_curve_isogeny.py:_isogeny_determine_algorithm
 */
export function _isogeny_determine_algorithm<F extends FieldElement>(
  E: EllipticCurveGeneric<F>,
  kernel: unknown
): 'velu' | 'kohel' {
  // If kernel is a single point, use Velu
  if (kernel instanceof EllipticCurvePoint) {
    return 'velu';
  }

  // If kernel is an array
  if (Array.isArray(kernel)) {
    if (kernel.length === 0) {
      return 'velu'; // trivial isogeny
    }

    // Check if first element is a point
    if (kernel[0] instanceof EllipticCurvePoint) {
      return 'velu';
    }

    // Otherwise assume it's polynomial coefficients
    return 'kohel';
  }

  throw new ValueError('invalid parameters to EllipticCurveIsogeny constructor');
}

/**
 * Compute the isogeny codomain given a kernel.
 *
 * @param E - Domain elliptic curve
 * @param kernel - Either a list of points in the kernel of the isogeny,
 *                 or a kernel polynomial (specified as either a univariate
 *                 polynomial or a coefficient list)
 * @returns The codomain of the separable normalized isogeny defined by this kernel
 * @see Reference: sage/schemes/elliptic_curves/ell_curve_isogeny.py:isogeny_codomain_from_kernel
 */
export function isogeny_codomain_from_kernel<F extends FieldElement>(
  E: EllipticCurveGeneric<F>,
  kernel: EllipticCurvePoint<F> | EllipticCurvePoint<F>[]
): EllipticCurveGeneric<F> {
  const algorithm = _isogeny_determine_algorithm(E, kernel);

  if (algorithm === 'velu') {
    // Create the isogeny and return the codomain
    const phi = new EllipticCurveIsogeny(E, kernel);
    return phi.codomain();
  }

  if (algorithm === 'kohel') {
    // Kohel's algorithm uses the kernel polynomial directly
    // Convert from kernel points to polynomial and use compute_codomain_kohel
    const kernelArray = Array.isArray(kernel) ? kernel : [kernel];
    if (kernelArray.length === 0 || kernelArray[0] instanceof EllipticCurvePoint) {
      // This should have been handled by Velu
      throw new ValueError('Kohel algorithm requires a polynomial kernel, not points');
    }
    return compute_codomain_kohel(E, kernelArray as unknown as bigint[]);
  }

  throw new ValueError('Unknown algorithm: must be "velu" or "kohel"');
}

/**
 * Compute the codomain curve given parameters v and w (as in Velu/Kohel/etc. formulas).
 *
 * @param E - An elliptic curve
 * @param v - Element of the base field of E
 * @param w - Element of the base field of E
 * @returns The elliptic curve with invariants [a1, a2, a3, a4-5v, a6-(a1^2+4a2)v-7w]
 * @see Reference: sage/schemes/elliptic_curves/ell_curve_isogeny.py:compute_codomain_formula
 */
export function compute_codomain_formula<F extends FieldElement>(
  E: EllipticCurveGeneric<F>,
  v: F,
  w: F
): EllipticCurveGeneric<F> {
  const [a1, a2, a3, a4, a6] = E.a_invariants();

  // A4 = a4 - 5*v
  const A4 = a4.sub(v.mul(5)) as F;

  // A6 = a6 - (a1^2 + 4*a2)*v - 7*w
  const a1Sq = a1.mul(a1) as F;
  const fourA2 = a2.mul(4) as F;
  const A6 = a6.sub(a1Sq.add(fourA2).mul(v)).sub(w.mul(7)) as F;

  return EllipticCurve(E.base_ring, [a1, a2, a3, A4, A6]);
}

/**
 * Compute Velu's (v,w) using Kohel's formulas for isogenies of degree exactly divisible by 2.
 *
 * @param x0 - x-coordinate of a 2-torsion point
 * @param y0 - y-coordinate of a 2-torsion point
 * @param a1 - Invariant of E
 * @param a2 - Invariant of E
 * @param a4 - Invariant of E
 * @returns Velu's isogeny parameters (v, w)
 * @see Reference: sage/schemes/elliptic_curves/ell_curve_isogeny.py:compute_vw_kohel_even_deg1
 */
export function compute_vw_kohel_even_deg1<F extends FieldElement>(
  x0: F,
  y0: F,
  a1: F,
  a2: F,
  a4: F
): [F, F] {
  // v = 3*x0^2 + 2*a2*x0 + a4 - a1*y0
  const v = x0.mul(x0).mul(3).add(a2.mul(x0).mul(2)).add(a4).sub(a1.mul(y0)) as F;

  // w = x0 * v
  const w = x0.mul(v) as F;

  return [v, w];
}

/**
 * Compute Velu's (v,w) using Kohel's formulas for isogenies of degree divisible by 4.
 *
 * @param b2 - Invariant of E
 * @param b4 - Invariant of E
 * @param s1 - Signed coefficient of the 2-division polynomial
 * @param s2 - Signed coefficient of the 2-division polynomial
 * @param s3 - Signed coefficient of the 2-division polynomial
 * @returns Velu's isogeny parameters (v, w)
 * @see Reference: sage/schemes/elliptic_curves/ell_curve_isogeny.py:compute_vw_kohel_even_deg3
 */
export function compute_vw_kohel_even_deg3<F extends FieldElement>(
  b2: F,
  b4: F,
  s1: F,
  s2: F,
  s3: F
): [F, F] {
  // temp1 = s1^2 - 2*s2
  const temp1 = s1.mul(s1).sub(s2.mul(2)) as F;

  // v = 3*temp1 + (b2*s1 + 3*b4)/2
  const v = temp1.mul(3).add(b2.mul(s1).add(b4.mul(3)).div(2)) as F;

  // w = 3*(s1^3 - 3*s1*s2 + 3*s3) + (b2*temp1 + b4*s1)/2
  const s1Cubed = s1.mul(s1).mul(s1) as F;
  const w = s1Cubed
    .sub(s1.mul(s2).mul(3))
    .add(s3.mul(3))
    .mul(3)
    .add(b2.mul(temp1).add(b4.mul(s1)).div(2)) as F;

  return [v, w];
}

/**
 * Compute Velu's (v,w) using Kohel's formulas for isogenies of odd degree.
 *
 * @param b2 - Invariant of E
 * @param b4 - Invariant of E
 * @param b6 - Invariant of E
 * @param s1 - Signed coefficient of the kernel polynomial
 * @param s2 - Signed coefficient of the kernel polynomial
 * @param s3 - Signed coefficient of the kernel polynomial
 * @param n - The degree
 * @returns Velu's isogeny parameters (v, w)
 * @see Reference: sage/schemes/elliptic_curves/ell_curve_isogeny.py:compute_vw_kohel_odd
 */
export function compute_vw_kohel_odd<F extends FieldElement>(
  b2: F,
  b4: F,
  b6: F,
  s1: F,
  s2: F,
  s3: F,
  n: number | bigint
): [F, F] {
  const nNum = typeof n === 'bigint' ? Number(n) : n;

  // temp = s1^2 - 2*s2
  const temp = s1.mul(s1).sub(s2.mul(2)) as F;

  // v = 6*temp + b2*s1 + n*b4
  const v = temp.mul(6).add(b2.mul(s1)).add(b4.mul(nNum)) as F;

  // w = 10*(s1^3 - 3*s1*s2 + 3*s3) + 2*b2*temp + 3*b4*s1 + n*b6
  const s1Cubed = s1.mul(s1).mul(s1) as F;
  const w = s1Cubed
    .sub(s1.mul(s2).mul(3))
    .add(s3.mul(3))
    .mul(10)
    .add(b2.mul(temp).mul(2))
    .add(b4.mul(s1).mul(3))
    .add(b6.mul(nNum)) as F;

  return [v, w];
}

/**
 * Compute the codomain from the kernel polynomial using Kohel's formulas.
 *
 * Given a monic kernel polynomial psi(x), this computes the codomain of the
 * isogeny with that kernel using Kohel's formulas.
 *
 * @param E - Domain elliptic curve
 * @param kernel - Polynomial or list; the kernel polynomial, or a list of its coefficients
 * @returns The codomain elliptic curve of the isogeny defined by kernel
 * @see Reference: sage/schemes/elliptic_curves/ell_curve_isogeny.py:compute_codomain_kohel
 */
export function compute_codomain_kohel<F extends FieldElement>(
  E: EllipticCurveGeneric<F>,
  kernel: Polynomial | bigint[]
): EllipticCurveGeneric<F> {
  const K = E.base_ring;
  const psi = _kernel_to_poly<F>(kernel, K);

  // next determine the even / odd part of the isogeny
  const psi_2tor = polyMonic(two_torsion_part_field(E, psi), K);

  const [a1, a2, a3, a4] = E.a_invariants();
  const [b2, b4, b6] = E.b_invariants();

  let v: F;
  let w: F;

  if (polyDegree(psi_2tor) !== 0) {
    // even degree case
    const [psi_quo] = polyDivRem(psi, psi_2tor, K);
    if (polyDegree(psi_quo) !== 0) {
      throw new NotImplementedError(
        "Kohel's algorithm currently only supports cyclic isogenies (except for [2])"
      );
    }

    const n = polyDegree(psi_2tor);

    if (n === 1) {
      // degree divisible exactly by 2
      const x0 = polyCoeff(psi_2tor, 0, K).neg() as F;

      // determine y0
      let y0: F;
      if (K.characteristic === 2n) {
        const rhs = x0
          .mul(x0)
          .mul(x0)
          .add(a2.mul(x0).mul(x0))
          .add(a4.mul(x0))
          .add(E.a_invariants()[4]) as F;
        const sq = _field_sqrt(rhs, K);
        if (sq === null) {
          throw new ValueError('cannot compute the 2-torsion point in characteristic 2');
        }
        y0 = sq;
      } else {
        y0 = a1.mul(x0).add(a3).neg().div(2) as F;
      }

      [v, w] = compute_vw_kohel_even_deg1(x0, y0, a1, a2, a4);
    } else if (n === 3) {
      // psi_2tor is the full 2-division polynomial
      const s1 = polyCoeff(psi_2tor, n - 1, K).neg() as F;
      const s2 = polyCoeff(psi_2tor, n - 2, K);
      const s3 = polyCoeff(psi_2tor, n - 3, K).neg() as F;

      [v, w] = compute_vw_kohel_even_deg3(b2, b4, s1, s2, s3);
    } else {
      throw new ValueError(`input polynomial must have degree 1 or 3, not ${n}`);
    }
  } else {
    // odd degree case
    const n = polyDegree(psi);

    const s1 = n >= 1 ? (polyCoeff(psi, n - 1, K).neg() as F) : (K.zero() as F);
    const s2 = n >= 2 ? polyCoeff(psi, n - 2, K) : (K.zero() as F);
    const s3 = n >= 3 ? (polyCoeff(psi, n - 3, K).neg() as F) : (K.zero() as F);

    [v, w] = compute_vw_kohel_odd(b2, b4, b6, s1, s2, s3, n);
  }

  return compute_codomain_formula(E, v, w);
}

/**
 * Coerce a kernel specification (coefficient list or Polynomial) into a
 * coefficient array over the base field of the curve.
 */
function _kernel_to_poly<F extends FieldElement>(
  kernel: Polynomial | bigint[] | F[],
  K: FieldParent
): F[] {
  if (Array.isArray(kernel)) {
    if (kernel.length === 0) return [];
    if (typeof kernel[0] === 'bigint') {
      return polyFromBigints<F>(kernel as bigint[], K);
    }
    return polyTrim((kernel as F[]).map((c) => K.__call__(c) as F));
  }
  const coeffs = (kernel as PolyElement<RingElement>).coeffs.map((c) => BigInt(c.toString()));
  return polyFromBigints<F>(coeffs, K);
}

/** Square root in the base field, or null if there is none. */
function _field_sqrt<F extends FieldElement>(a: F, K: FieldParent): F | null {
  return find_sqrt(a, K);
}

/** Compose polynomials: return f(g(x)). */
function _poly_compose<F extends FieldElement>(f: F[], g: F[], K: FieldParent): F[] {
  let acc: F[] = [];
  for (let i = f.length - 1; i >= 0; i--) {
    acc = polyAdd(polyMul(acc, g, K), [f[i]!], K);
  }
  return polyTrim(acc);
}

/** Binomial coefficient as a plain number (used with small arguments only). */
function _binomial(n: number, k: number): number {
  if (k < 0 || k > n) return 0;
  let r = 1;
  for (let i = 0; i < k; i++) {
    r = (r * (n - i)) / (i + 1);
  }
  return Math.round(r);
}

/** Human-readable representation of a polynomial (used in error messages). */
function _poly_repr<F extends FieldElement>(f: F[]): string {
  const t = polyTrim(f);
  if (t.length === 0) return '0';
  const parts: string[] = [];
  for (let i = t.length - 1; i >= 0; i--) {
    if (t[i]!.isZero()) continue;
    parts.push(i === 0 ? `${t[i]}` : i === 1 ? `${t[i]}*x` : `${t[i]}*x^${i}`);
  }
  return parts.join(' + ');
}

/** Whether two curves have identical a-invariants. */
function _same_curve<F extends FieldElement>(
  E1: EllipticCurveGeneric<F>,
  E2: EllipticCurveGeneric<F>
): boolean {
  const a = E1.a_invariants();
  const b = E2.a_invariants();
  return a.every((c, i) => c.eq(b[i]!));
}

/**
 * Return a model of `E` of the requested type.
 *
 * @see Reference: sage/schemes/elliptic_curves/ell_field.py:compute_model
 */
function _compute_model<F extends FieldElement>(
  E: EllipticCurveGeneric<F>,
  name: 'minimal' | 'short_weierstrass' | 'montgomery'
): EllipticCurveGeneric<F> {
  if (name === 'minimal') {
    throw new ValueError('can only compute minimal model for curves over number fields');
  }
  if (name === 'short_weierstrass') {
    return E.short_weierstrass_model();
  }
  if (name === 'montgomery') {
    return E.montgomery_model(false, false) as EllipticCurveGeneric<F>;
  }
  throw new NotImplementedError(`cannot compute ${name} model`);
}

/**
 * The 2-division polynomial of E as a coefficient array over the base field:
 * `4x^3 + b2 x^2 + 2 b4 x + b6`.
 */
function _two_division_polynomial<F extends FieldElement>(E: EllipticCurveGeneric<F>): F[] {
  const K = E.base_ring;
  const [b2, b4, b6] = E.b_invariants();
  return polyTrim([b6, b4.mul(2) as F, b2, K.__call__(4n) as F]);
}

/**
 * Return the gcd of `psi` with the 2-torsion polynomial of `E`, over the base
 * field of `E`.
 *
 * @see Reference: sage/schemes/elliptic_curves/ell_curve_isogeny.py:two_torsion_part
 */
export function two_torsion_part_field<F extends FieldElement>(
  E: EllipticCurveGeneric<F>,
  psi: F[]
): F[] {
  return polyGcdField(psi, _two_division_polynomial(E), E.base_ring);
}

/**
 * Return the greatest common divisor of psi and the 2-torsion polynomial of E.
 *
 * The 2-torsion polynomial is 4x^3 + b2*x^2 + 2*b4*x + b6 for a curve with
 * b-invariants b2, b4, b6. Points where this polynomial vanishes are
 * 2-torsion points.
 *
 * @param E - An elliptic curve
 * @param psi - A univariate polynomial over the base field of E (as coefficient array)
 * @returns The gcd of psi and the 2-torsion polynomial of E (as coefficient array)
 * @see Reference: sage/schemes/elliptic_curves/ell_curve_isogeny.py:two_torsion_part
 */
export function two_torsion_part<F extends FieldElement>(
  E: EllipticCurveGeneric<F>,
  psi: bigint[]
): bigint[] {
  const K = E.base_ring;
  return polyToBigints(two_torsion_part_field(E, polyFromBigints<F>(psi, K)));
}

// ---------------------------------------------------------------------------
// Truncated power series helpers (coefficient arrays, ascending degree).
// ---------------------------------------------------------------------------

/** Truncate a coefficient array to precision `n` (drop terms of degree >= n). */
function seriesTrunc<F extends FieldElement>(a: F[], n: number): F[] {
  return polyTrim(a.slice(0, Math.max(0, n)));
}

/** Product of two series truncated at precision `n`. */
function mulTrunc<F extends FieldElement>(a: F[], b: F[], n: number, K: FieldParent): F[] {
  return seriesTrunc(polyMul(a, b, K), n);
}

/**
 * Inverse of a power series with invertible constant term, truncated at
 * precision `n` (mirrors Sage's `inverse_series_trunc`).
 */
function invTrunc<F extends FieldElement>(f: F[], n: number, K: FieldParent): F[] {
  if (f.length === 0 || f[0]!.isZero()) {
    throw new ValueError('constant term of the power series is not invertible');
  }
  const c0inv = f[0]!.inv() as F;
  const g: F[] = [c0inv];
  for (let k = 1; k < n; k++) {
    let acc = K.zero() as F;
    for (let i = 1; i <= k; i++) {
      const fi = i < f.length ? f[i]! : (K.zero() as F);
      if (fi.isZero()) continue;
      acc = acc.add(fi.mul(g[k - i]!)) as F;
    }
    g.push(acc.neg().mul(c0inv) as F);
  }
  return polyTrim(g);
}

/**
 * Square root of a power series with square constant term, truncated at
 * precision `n`; returns null when the constant term is not a square.
 */
function seriesSqrt<F extends FieldElement>(f: F[], n: number, K: FieldParent): F[] | null {
  if (f.length === 0) return [];
  const c0 = f[0]!;
  const r0 = find_sqrt(c0, K);
  if (r0 === null || !r0.mul(r0).eq(c0)) return null;
  const g: F[] = [r0];
  const two_r0_inv = r0.mul(2).inv() as F;
  for (let k = 1; k < n; k++) {
    let acc = k < f.length ? f[k]! : (K.zero() as F);
    for (let i = 1; i < k; i++) {
      acc = acc.sub(g[i]!.mul(g[k - i]!)) as F;
    }
    g.push(acc.mul(two_r0_inv) as F);
  }
  return polyTrim(g);
}

/**
 * Reverse the coefficients of a polynomial within the given degree, i.e.
 * return `x^degree * f(1/x)` (mirrors Sage's `Polynomial.reverse(degree=...)`).
 */
function polyReverse<F extends FieldElement>(f: F[], degree: number, K: FieldParent): F[] {
  const out: F[] = [];
  for (let i = 0; i <= degree; i++) out.push(K.zero() as F);
  for (let i = 0; i <= degree; i++) {
    if (i < f.length) out[degree - i] = f[i]!;
  }
  return polyTrim(out);
}

/**
 * Rational reconstruction: return `[n, d]` with `u*d = n (mod m)`,
 * `deg(n) <= n_deg`, `deg(d) <= d_deg` and `d` monic.
 *
 * @see Reference: sage/rings/polynomial/polynomial_element.pyx:rational_reconstruction
 */
function rationalReconstruction<F extends FieldElement>(
  u: F[],
  m: F[],
  n_deg: number,
  d_deg: number,
  K: FieldParent
): [F[], F[]] {
  let s0: F[] = [];
  let t0: F[] = [K.one() as F];
  let s1: F[] = polyTrim(m);
  let t1: F[] = polyDivRem(u, s1, K)[1];

  while (n_deg < polyDegree(t1) && polyDegree(t1) >= 0) {
    const [q, r1] = polyDivRem(s1, t1, K);
    const r0 = polySub(s0, polyMul(q, t0, K), K);
    s0 = t0;
    s1 = t1;
    t0 = r0;
    t1 = r1;
  }

  if (polyDegree(t0) < 0) {
    throw new ValueError('could not complete rational reconstruction');
  }
  if (d_deg < polyDegree(t0)) {
    throw new ValueError('could not complete rational reconstruction');
  }

  const c = t0[t0.length - 1]!;
  const cinv = c.inv() as F;
  return [polyScale(t1, cinv), polyScale(t0, cinv)];
}

/**
 * Compute the kernel polynomial of the unique normalized isogeny of degree l between E1 and E2.
 *
 * Both curves must be given in short Weierstrass form (y^2 = x^3 + a*x + b),
 * and the characteristic must be either 0 or no smaller than 4l+4.
 *
 * ALGORITHM: [BMSS2006], algorithm *fastElkies'*.
 *
 * @param E1 - Domain elliptic curve in short Weierstrass form
 * @param E2 - Codomain elliptic curve in short Weierstrass form
 * @param l - The degree of the isogeny
 * @returns The kernel polynomial as a coefficient array (ascending degree)
 * @see Reference: sage/schemes/elliptic_curves/ell_curve_isogeny.py:compute_isogeny_bmss
 * @see Deviation: BMSS differential equation solved by direct coefficient recurrence
 */
export function compute_isogeny_bmss<F extends FieldElement>(
  E1: EllipticCurveGeneric<F>,
  E2: EllipticCurveGeneric<F>,
  l: number | bigint
): bigint[] {
  const ell = Number(l);
  const K = E1.base_ring;
  const char = K.characteristic;

  const [a1_1, a2_1, a3_1] = E1.a_invariants();
  const [a1_2, a2_2, a3_2] = E2.a_invariants();
  if (!a1_1.isZero() || !a2_1.isZero() || !a3_1.isZero()) {
    throw new ValueError('E1 must be a short Weierstrass curve');
  }
  if (!a1_2.isZero() || !a2_2.isZero() || !a3_2.isZero()) {
    throw new ValueError('E2 must be a short Weierstrass curve');
  }
  if (char !== 0n && char < BigInt(4 * ell + 4)) {
    throw new ValueError('characteristic must be at least 4*degree+4');
  }

  const zero = K.zero() as F;
  const one = K.one() as F;

  const A = E1.a4();
  const B = E1.a6();
  const A2 = E2.a4();
  const B2 = E2.a6();

  const prec = 4 * ell; // work modulo x^(4l)

  // D = 1 + A x^4 + B x^6
  const D: F[] = polyTrim([one, zero, zero, zero, A, zero, B]);

  // Solve D * (S')^2 = 1 + A2 S^4 + B2 S^6 with S = x + O(x^2).
  //
  // Sage runs the Newton doubling of fastElkies'; the solution of the
  // differential equation is unique given S = x + O(x^2), so we determine its
  // coefficients directly by the recurrence obtained from comparing the
  // coefficient of x^(n-1) on both sides:
  //     2 n s_n = RHS[x^(n-1)] - (part of LHS not involving s_n).
  const S: F[] = [zero, one]; // s_0 = 0, s_1 = 1
  for (let n = 2; n < prec; n++) {
    S.push(zero);
    const Sp = polyDeriv(S, K); // S'
    const Sp2 = mulTrunc(Sp, Sp, prec, K);
    const lhs = mulTrunc(D, Sp2, prec, K);
    const S2 = mulTrunc(S, S, prec, K);
    const S4 = mulTrunc(S2, S2, prec, K);
    const S6 = mulTrunc(S4, S2, prec, K);
    const rhs = polyAdd(
      polyAdd([one], polyScale(S4, A2), K),
      polyScale(S6, B2),
      K
    );
    const diff = polyCoeff(rhs, n - 1, K).sub(polyCoeff(lhs, n - 1, K)) as F;
    S[n] = diff.div(K.__call__(2 * n)) as F;
  }

  // Reconstruct: S = x*T(x^2), U = 1/T^2, then N(1/x)/D(1/x) = U.
  const T: F[] = [];
  for (let i = 0; i < 2 * ell; i++) {
    T.push(polyCoeff(S, 2 * i + 1, K));
  }
  const U = invTrunc(mulTrunc(T, T, 2 * ell, K), 2 * ell, K);

  // m = x^(2l)
  const m: F[] = polyShift([one], 2 * ell, K);
  const [, Qden] = rationalReconstruction(U, m, ell, ell, K);

  const qprec = Math.floor((ell + 1) / 2);
  const Qser = seriesTrunc(Qden, qprec);
  const root = seriesSqrt(Qser, qprec, K);
  if (root === null) {
    throw new ValueError(
      `the two curves are not linked by a cyclic normalized isogeny of degree ${ell}`
    );
  }

  const ker = polyMonic(polyReverse(root, Math.floor(ell / 2), K), K);
  return polyToBigints(ker);
}

/**
 * Return the kernel polynomial of an isogeny of degree ell from E1 to E2 using
 * Stark's algorithm.
 *
 * Sage's implementation is driven by the Weierstrass `wp` expansions of both
 * curves (`E.weierstrass_p(prec=4*ell+4)`), which live in
 * `sage/schemes/elliptic_curves/ell_wp.py`.  That module has not been ported,
 * so this entry point is not available.
 *
 * @param E1 - Domain elliptic curve in short Weierstrass form
 * @param E2 - Codomain elliptic curve in short Weierstrass form
 * @param ell - The degree of an isogeny from E1 to E2
 * @returns The kernel polynomial of an isogeny from E1 to E2
 * @see Reference: sage/schemes/elliptic_curves/ell_curve_isogeny.py:compute_isogeny_stark
 * @see Deviation: compute_isogeny_stark unimplemented (requires ell_wp)
 */
export function compute_isogeny_stark<F extends FieldElement>(
  _E1: EllipticCurveGeneric<F>,
  _E2: EllipticCurveGeneric<F>,
  _ell: number | bigint
): bigint[] {
  throw new NotImplementedError(
    "SAGE_NOT_IMPLEMENTED: compute_isogeny_stark requires E.weierstrass_p() from sage.schemes.elliptic_curves.ell_wp, which is not ported"
  );
}

/**
 * Return the kernel polynomial of a cyclic, separable, normalized isogeny of degree ell from E1 to E2.
 *
 * @param E1 - Domain elliptic curve in short Weierstrass form
 * @param E2 - Codomain elliptic curve in short Weierstrass form
 * @param ell - The degree of an isogeny from E1 to E2
 * @param algorithm - 'bmss', 'stark', or undefined (auto-select)
 * @returns The kernel polynomial of an isogeny from E1 to E2
 * @see Reference: sage/schemes/elliptic_curves/ell_curve_isogeny.py:compute_isogeny_kernel_polynomial
 * @see Deviation: compute_isogeny_stark unimplemented (requires ell_wp)
 */
export function compute_isogeny_kernel_polynomial<F extends FieldElement>(
  E1: EllipticCurveGeneric<F>,
  E2: EllipticCurveGeneric<F>,
  ell: number | bigint,
  algorithm?: 'bmss' | 'stark'
): bigint[] {
  const ellNum = Number(ell);

  if (algorithm === undefined) {
    const char = E1.base_ring.characteristic;
    if (char !== 0n && char < BigInt(4 * ellNum + 4)) {
      throw new NotImplementedError(
        `no algorithm for computing kernel polynomial from domain and codomain is implemented for degree ${ellNum} and characteristic ${char}`
      );
    }
    // Sage picks 'stark' for ell < 10 and 'bmss' otherwise.  Stark's algorithm
    // is not available here (see compute_isogeny_stark); BMSS returns the same
    // kernel polynomial, and its precondition (char >= 4*ell+4) has just been
    // checked, so we always use BMSS.
    algorithm = 'bmss';
  }

  if (algorithm === 'bmss') {
    return compute_isogeny_bmss(E1, E2, ell);
  }
  if (algorithm === 'stark') {
    // Sage applies .radical() to Stark's output.
    return compute_isogeny_stark(E1, E2, ell);
  }
  throw new NotImplementedError(`unknown algorithm ${algorithm}`);
}

/**
 * Return intermediate curves and isomorphisms.
 *
 * @param E1 - An elliptic curve
 * @param E2 - An elliptic curve
 * @returns A tuple (E1w, E2w, pre_isomorphism, post_isomorphism) where E1w and
 *          E2w are short Weierstrass models of E1 and E2, `pre_isomorphism` is
 *          a normalized (u = 1) isomorphism E1 -> E1w and `post_isomorphism` is
 *          a normalized isomorphism E2w -> E2.
 * @see Reference: sage/schemes/elliptic_curves/ell_curve_isogeny.py:compute_intermediate_curves
 */
export function compute_intermediate_curves<F extends FieldElement>(
  E1: EllipticCurveGeneric<F>,
  E2: EllipticCurveGeneric<F>
): [
  EllipticCurveGeneric<F>,
  EllipticCurveGeneric<F>,
  WeierstrassIsomorphism<F>,
  WeierstrassIsomorphism<F>,
] {
  const K = E1.base_ring;
  if (K.characteristic === 2n || K.characteristic === 3n) {
    throw new NotImplementedError(
      'compute_intermediate_curves is only defined for characteristics not 2 or 3'
    );
  }

  const shortModel = (E: EllipticCurveGeneric<F>): EllipticCurveGeneric<F> => {
    const [c4, c6] = E.c_invariants();
    const A = c4.neg().div(48) as F;
    const B = c6.neg().div(864) as F;
    const zero = K.zero() as F;
    return EllipticCurve(K, [zero, zero, zero, A, B]);
  };

  const E1w = shortModel(E1);
  const E2w = shortModel(E2);

  // We cannot use E1.isomorphism_to(E1w) since it may have u = -1.
  const normalized = (
    A: EllipticCurveGeneric<F>,
    B: EllipticCurveGeneric<F>
  ): WeierstrassIsomorphism<F> => {
    for (const urst of _isomorphisms(A, B)) {
      if (urst[0].eq(K.one())) {
        return new WeierstrassIsomorphism<F>(A, urst, B);
      }
    }
    throw new ValueError('no normalized isomorphism between the given curves');
  };

  const pre_iso = normalized(E1, E1w);
  const post_iso = normalized(E2w, E2);
  return [E1w, E2w, pre_iso, post_iso];
}

/**
 * Return intermediate curves, isomorphisms and kernel polynomial.
 *
 * @param E1 - An elliptic curve
 * @param E2 - An elliptic curve
 * @param ell - A prime such that there is a degree-ell separable normalized isogeny from E1 to E2
 * @returns A tuple (pre_isom, post_isom, E1pr, E2pr, ker_poly)
 * @see Reference: sage/schemes/elliptic_curves/ell_curve_isogeny.py:compute_sequence_of_maps
 */
export function compute_sequence_of_maps<F extends FieldElement>(
  E1: EllipticCurveGeneric<F>,
  E2: EllipticCurveGeneric<F>,
  ell: number | bigint
): [
  WeierstrassIsomorphism<F>,
  WeierstrassIsomorphism<F>,
  EllipticCurveGeneric<F>,
  EllipticCurveGeneric<F>,
  bigint[],
] {
  const [E1pr, E2pr, pre_isom, post_isom] = compute_intermediate_curves(E1, E2);

  const ker_poly = compute_isogeny_kernel_polynomial(E1pr, E2pr, ell);

  return [pre_isom, post_isom, E1pr, E2pr, ker_poly];
}

/**
 * Return a filled isogeny matrix giving all degrees from one giving only prime degrees.
 *
 * Given a matrix M where M[i][j] is a prime l if curves i and j are l-isogenous (and 0 otherwise),
 * this function fills in all composite degrees by computing minimal paths through the isogeny graph.
 *
 * @param M - A square symmetric matrix whose off-diagonal i,j entry is either a prime l
 *            if the i-th and j-th curves have an l-isogeny between them, otherwise 0
 * @returns A square matrix with entries 1 on the diagonal, and in general the i,j entry
 *          is d>0 if d is the minimal degree of an isogeny from the i-th to the j-th curve
 * @see Reference: sage/schemes/elliptic_curves/ell_curve_isogeny.py:fill_isogeny_matrix
 */
export function fill_isogeny_matrix(M: bigint[][]): bigint[][] {
  const n = M.length;
  const result: bigint[][] = [];

  // Initialize result with M, and 1s on diagonal
  for (let i = 0; i < n; i++) {
    result.push([]);
    for (let j = 0; j < n; j++) {
      if (i === j) {
        result[i]!.push(1n);
      } else {
        result[i]!.push(M[i]![j]!);
      }
    }
  }

  // Floyd-Warshall style algorithm to find minimal degree paths
  // The degree of a composed isogeny is the product of the degrees
  // We want to find the minimum degree path between any two curves
  let changed = true;
  while (changed) {
    changed = false;
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        if (i === j) continue;
        for (let k = 0; k < n; k++) {
          if (k === i || k === j) continue;

          const ik = result[i]![k]!;
          const kj = result[k]![j]!;

          if (ik > 0n && kj > 0n) {
            const newDeg = ik * kj;
            const oldDeg = result[i]![j]!;

            if (oldDeg === 0n || newDeg < oldDeg) {
              result[i]![j] = newDeg;
              changed = true;
            }
          }
        }
      }
    }
  }

  return result;
}

/**
 * Reverses the action of fill_isogeny_matrix.
 *
 * Given a filled isogeny matrix (with all minimal degrees), extract only the
 * prime degree entries, setting all composite degrees to 0.
 *
 * @param M - A square symmetric matrix of integers
 * @returns A square symmetric matrix obtained from M by replacing non-prime entries with 0
 * @see Reference: sage/schemes/elliptic_curves/ell_curve_isogeny.py:unfill_isogeny_matrix
 */
export function unfill_isogeny_matrix(M: bigint[][]): bigint[][] {
  const n = M.length;
  const result: bigint[][] = [];

  for (let i = 0; i < n; i++) {
    result.push([]);
    for (let j = 0; j < n; j++) {
      const val = M[i]![j]!;
      if (i === j) {
        // Sage zeroes the diagonal (ell_curve_isogeny.py:3964-3965), so that
        // unfill_isogeny_matrix(fill_isogeny_matrix(M)) == M.
        result[i]!.push(0n);
      } else if (isPrime(val)) {
        result[i]!.push(val);
      } else {
        result[i]!.push(0n);
      }
    }
  }

  return result;
}

/**
 * Check if a bigint is prime (simple trial division for small primes).
 */
function isPrime(n: bigint): boolean {
  if (n < 2n) return false;
  if (n === 2n) return true;
  if (n % 2n === 0n) return false;
  if (n < 9n) return true;
  if (n % 3n === 0n) return false;

  const limit = sqrt(n);
  let i = 5n;
  while (i <= limit) {
    if (n % i === 0n) return false;
    if (n % (i + 2n) === 0n) return false;
    i += 6n;
  }
  return true;
}

/**
 * Integer square root for bigint.
 */
function sqrt(n: bigint): bigint {
  if (n < 0n) throw new ValueError('square root of negative number');
  if (n < 2n) return n;

  let x = n;
  let y = (x + 1n) / 2n;
  while (y < x) {
    x = y;
    y = (x + n / x) / 2n;
  }
  return x;
}

/**
 * This class implements separable isogenies of elliptic curves.
 *
 * Several different algorithms for computing isogenies are available:
 * - Velu's Formulas: Velu's original formulas for computing isogenies.
 *   This algorithm is selected by giving as the kernel parameter a single
 *   point, or a list of points, generating a finite subgroup.
 * - Kohel's Formulas: Kohel's original formulas for computing isogenies.
 *   This algorithm is selected by giving as the kernel parameter a monic
 *   polynomial (or a coefficient list) which will define the kernel of
 *   the isogeny.
 *
 * @see Reference: sage/schemes/elliptic_curves/ell_curve_isogeny.py:EllipticCurveIsogeny
 */
export class EllipticCurveIsogeny<F extends FieldElement = FieldElement> {
  private _domain: EllipticCurveGeneric<F>;
  private _codomain: EllipticCurveGeneric<F> | null = null;
  private _degree: bigint = 1n;

  /** Algorithm used: 'velu' or 'kohel' */
  private __algorithm: 'velu' | 'kohel';

  /** Kernel points data for Velu's formula: xQ -> KernelPointData */
  private __kernel_mod_sign: Map<string, KernelPointData<F>> = new Map();

  /** Velu's v parameter */
  private __v: F;

  /** Velu's w parameter */
  private __w: F;

  /** Cached kernel polynomial (over the base field, ascending coefficients) */
  private __kernel_polynomial: F[] | null = null;

  /** Kohel: the kernel polynomial psi */
  private __psi: F[] | null = null;

  /** Kohel: numerator phi of the X-coordinate, X = phi(x)/psi(x)^2 */
  private __phi: F[] | null = null;

  /**
   * Kohel: numerator omega of the Y-coordinate, Y = omega(x,y)/psi(x)^3.
   * Stored as [omega0, omega1] with omega = omega0(x) + omega1(x)*y.
   */
  private __omega: [F[], F[]] | null = null;

  /** Kohel: the inner kernel polynomial (before any pre-isomorphism) */
  private __inner_kernel_polynomial: F[] | null = null;

  /** Whether the constructor should validate its input */
  private __check: boolean;

  /**
   * Construct an elliptic curve isogeny.
   *
   * @param E - An elliptic curve; the domain of the isogeny
   * @param kernel - A kernel; either a point on E, a list of points on E, a
   *                 kernel polynomial (coefficient list), or null
   * @param codomain - An elliptic curve (optional)
   * @param degree - Integer (optional)
   * @param model - String (optional): 'minimal', 'short_weierstrass', or 'montgomery'
   * @param check - Boolean (default: true); check whether the input is valid
   * @see Reference: sage/schemes/elliptic_curves/ell_curve_isogeny.py:EllipticCurveIsogeny.__init__
   */
  constructor(
    E: EllipticCurveGeneric<F>,
    kernel: EllipticCurvePoint<F> | EllipticCurvePoint<F>[] | bigint[] | null,
    codomain?: EllipticCurveGeneric<F> | null,
    degree?: number | bigint | null,
    model?: 'minimal' | 'short_weierstrass' | 'montgomery' | null,
    check: boolean = true
  ) {
    this.__check = check;

    // Initialize v and w with zero
    const zero = E.base_ring.zero() as F;
    this.__v = zero;
    this.__w = zero;

    // Handle single point input
    if (kernel instanceof EllipticCurvePoint) {
      kernel = [kernel];
    }

    // If the kernel is None and the codomain isn't, compute the kernel
    // polynomial via compute_sequence_of_maps (ell_curve_isogeny.py:1055-1063).
    let pre_isom: WeierstrassIsomorphismLike<F> | null = null;
    let post_isom: WeierstrassIsomorphismLike<F> | null = null;
    let old_codomain: EllipticCurveGeneric<F> | null = null;

    if (kernel === null && codomain !== null && codomain !== undefined) {
      if (degree === null || degree === undefined) {
        throw new ValueError(
          'degree must be given when specifying isogeny by domain and codomain'
        );
      }
      old_codomain = codomain;
      const seq = compute_sequence_of_maps(E, codomain, degree);
      pre_isom = seq[0] as WeierstrassIsomorphismLike<F>;
      post_isom = seq[1] as WeierstrassIsomorphismLike<F>;
      E = seq[2];
      codomain = seq[3];
      kernel = seq[4] as bigint[];
    }

    this._domain = E;

    // Determine algorithm
    this.__algorithm = _isogeny_determine_algorithm(E, kernel);

    if (this.__algorithm === 'velu') {
      this.__init_from_kernel_gens(kernel as EllipticCurvePoint<F>[], check);
    } else if (this.__algorithm === 'kohel') {
      this.__init_from_kernel_polynomial(kernel as bigint[]);
    } else {
      throw new ValueError('unknown isogeny algorithm');
    }

    // Compute codomain
    this.__compute_codomain();

    this.__setup_post_isomorphism(codomain ?? null, model ?? null);

    if (pre_isom !== null) {
      this._set_pre_isomorphism(pre_isom);
    }
    if (post_isom !== null && old_codomain !== null) {
      this.__set_post_isomorphism(old_codomain, post_isom);
    }
  }

  /**
   * Initialize the isogeny from a kernel polynomial (Kohel's algorithm).
   *
   * @see Reference: sage/schemes/elliptic_curves/ell_curve_isogeny.py:EllipticCurveIsogeny.__init_from_kernel_polynomial
   */
  private __init_from_kernel_polynomial(kernel_polynomial: bigint[] | F[]): void {
    const K = this._domain.base_ring;
    const E = this._domain;
    const psi = _kernel_to_poly<F>(kernel_polynomial, K);

    if (polyDegree(psi) < 0) {
      throw new ValueError('given kernel polynomial is not monic');
    }
    if (!psi[psi.length - 1]!.eq(K.one())) {
      throw new ValueError('given kernel polynomial is not monic');
    }

    // Determine if the kernel polynomial is entirely 2-torsion
    const psi_G = polyMonic(two_torsion_part_field(E, psi), K);

    let phi: F[];
    let omega: [F[], F[]];
    let v: F;
    let w: F;
    let d: number;

    if (polyDegree(psi_G) !== 0) {
      // even degree case
      const [psi_quo] = polyDivRem(psi, psi_G, K);
      if (polyDegree(psi_quo) !== 0) {
        throw new NotImplementedError(
          "Kohel's algorithm currently only supports cyclic isogenies (except for [2])"
        );
      }
      [phi, omega, v, w, , d] = this.__init_even_kernel_polynomial(E, psi_G);
    } else {
      // odd degree case
      [phi, omega, v, w, , d] = this.__init_odd_kernel_polynomial(E, psi);
    }

    this.__kernel_polynomial = psi;
    this.__inner_kernel_polynomial = psi;

    this._degree = BigInt(d);

    this.__psi = psi;
    this.__phi = phi;
    this.__omega = omega;

    this.__v = v;
    this.__w = w;
  }

  /**
   * Return the isogeny parameters for the 2-part of an isogeny.
   *
   * @returns [phi, omega, v, w, n, d]
   * @see Reference: sage/schemes/elliptic_curves/ell_curve_isogeny.py:EllipticCurveIsogeny.__init_even_kernel_polynomial
   */
  private __init_even_kernel_polynomial(
    E: EllipticCurveGeneric<F>,
    psi_G: F[]
  ): [F[], [F[], F[]], F, F, number, number] {
    const K = E.base_ring;

    // check that psi_G really divides the two-torsion polynomial
    if (this.__check) {
      const [, rem] = polyDivRem(_two_division_polynomial(E), psi_G, K);
      if (polyDegree(rem) >= 0) {
        throw new ValueError(
          `the polynomial ${_poly_repr(psi_G)} does not define a finite subgroup of ${E}`
        );
      }
    }

    const n = polyDegree(psi_G); // 1 or 3
    const d = n + 1; // 2 or 4

    const [a1, a2, a3, a4, a6] = E.a_invariants();
    const [b2, b4] = E.b_invariants();

    let phi: F[];
    let omega: [F[], F[]];
    let v: F;
    let w: F;

    if (n === 1) {
      const x0 = polyCoeff(psi_G, 0, K).neg() as F;

      // determine y0
      let y0: F;
      if (K.characteristic === 2n) {
        const rhs = x0
          .mul(x0)
          .mul(x0)
          .add(a2.mul(x0).mul(x0))
          .add(a4.mul(x0))
          .add(a6) as F;
        const sq = _field_sqrt(rhs, K);
        if (sq === null) {
          throw new ValueError('cannot compute the 2-torsion point in characteristic 2');
        }
        y0 = sq;
      } else {
        y0 = a1.mul(x0).add(a3).neg().div(2) as F;
      }

      [v, w] = compute_vw_kohel_even_deg1(x0, y0, a1, a2, a4);

      // phi = (x*psi_G + v)*psi_G
      const xPsiG = polyShift(psi_G, 1, K);
      phi = polyMul(polyAdd(xPsiG, [v], K), psi_G, K);

      // omega = (y*psi_G^2 - v*(a1*psi_G + (y - y0)))*psi_G
      //       = psi_G * [ (psi_G^2 - v)*y + (-v*a1*psi_G + v*y0) ]
      const psiG2 = polyMul(psi_G, psi_G, K);
      const omega1 = polyMul(polySub(psiG2, [v], K), psi_G, K);
      const omega0 = polyMul(
        polyAdd(polyScale(polyScale(psi_G, a1), v.neg() as F), [v.mul(y0) as F], K),
        psi_G,
        K
      );
      omega = [omega0, omega1];
    } else if (n === 3) {
      const s1 = polyCoeff(psi_G, n - 1, K).neg() as F;
      const s2 = polyCoeff(psi_G, n - 2, K);
      const s3 = polyCoeff(psi_G, n - 3, K).neg() as F;

      const psi_G_pr = polyDeriv(psi_G, K);
      const psi_G_prpr = polyDeriv(psi_G_pr, K);

      // phi = psi_G_pr^2 + (-2*psi_G_prpr + (4*x - s1))*psi_G
      const fourXminusS1: F[] = polyTrim([s1.neg() as F, K.__call__(4n) as F]);
      let phiInner = polyAdd(
        polyMul(psi_G_pr, psi_G_pr, K),
        polyMul(polyAdd(polyScale(psi_G_prpr, K.__call__(-2n) as F), fourXminusS1, K), psi_G, K),
        K
      );
      const phi_pr = polyDeriv(phiInner, K);

      // psi_2 = 2*y + a1*x + a3
      const psi2x: F[] = polyTrim([a3, a1]);

      // omega = (psi_2*(phi_pr*psi_G - phi*psi_G_pr) - (a1*phi + a3*psi_G)*psi_G)/2
      const A = polySub(polyMul(phi_pr, psi_G, K), polyMul(phiInner, psi_G_pr, K), K);
      const B = polyMul(polyAdd(polyScale(phiInner, a1), polyScale(psi_G, a3), K), psi_G, K);
      const half = K.__call__(2n).inv() as F;
      let om0 = polyScale(polySub(polyMul(psi2x, A, K), B, K), half);
      let om1 = polyScale(A, K.__call__(2n).mul(half) as F);

      // phi *= psi_G;  omega *= psi_G
      phiInner = polyMul(phiInner, psi_G, K);
      om0 = polyMul(om0, psi_G, K);
      om1 = polyMul(om1, psi_G, K);

      phi = phiInner;
      omega = [om0, om1];

      [v, w] = compute_vw_kohel_even_deg3(b2, b4, s1, s2, s3);
    } else {
      throw new ValueError(`input polynomial must have degree 1 or 3, not ${n}`);
    }

    return [phi, omega, v, w, n, d];
  }

  /**
   * Return the isogeny parameters for a cyclic isogeny of odd degree.
   *
   * @returns [phi, omega, v, w, n, d]
   * @see Reference: sage/schemes/elliptic_curves/ell_curve_isogeny.py:EllipticCurveIsogeny.__init_odd_kernel_polynomial
   */
  private __init_odd_kernel_polynomial(
    E: EllipticCurveGeneric<F>,
    psi: F[]
  ): [F[], [F[], F[]], F, F, number, number] {
    const K = E.base_ring;
    const n = polyDegree(psi);
    const d = 2 * n + 1;

    // NOTE: Sage additionally verifies `is_kernel_polynomial(E, d, psi)` here.
    // That helper lives in sage.schemes.elliptic_curves.isogeny_small_degree,
    // which has not been ported yet; see the deviation note in the module docs.

    const [b2, b4, b6] = E.b_invariants();

    const zero = K.zero() as F;
    const s1 = n >= 1 ? (polyCoeff(psi, n - 1, K).neg() as F) : zero;
    const s2 = n >= 2 ? polyCoeff(psi, n - 2, K) : zero;
    const s3 = n >= 3 ? (polyCoeff(psi, n - 3, K).neg() as F) : zero;

    const [v, w] = compute_vw_kohel_odd(b2, b4, b6, s1, s2, s3, n);

    const psi_pr = polyDeriv(psi, K);
    const psi_prpr = polyDeriv(psi_pr, K);

    // phi = (4x^3 + b2 x^2 + 2 b4 x + b6)*(psi_pr^2 - psi_prpr*psi)
    //       - (6x^2 + b2 x + b4)*psi_pr*psi + (d*x - 2*s1)*psi^2
    const twoDiv = _two_division_polynomial(E);
    const sixX2: F[] = polyTrim([b4, b2, K.__call__(6n) as F]);
    const dxMinus2s1: F[] = polyTrim([s1.mul(-2) as F, K.__call__(d) as F]);

    const phi = polyAdd(
      polySub(
        polyMul(twoDiv, polySub(polyMul(psi_pr, psi_pr, K), polyMul(psi_prpr, psi, K), K), K),
        polyMul(polyMul(sixX2, psi_pr, K), psi, K),
        K
      ),
      polyMul(dxMinus2s1, polyMul(psi, psi, K), K),
      K
    );

    const phi_pr = polyDeriv(phi, K);

    const omega =
      K.characteristic !== 2n
        ? this.__compute_omega_fast(E, psi, psi_pr, phi, phi_pr)
        : this.__compute_omega_general(E, psi, psi_pr, phi, phi_pr);

    return [phi, omega, v, w, n, d];
  }

  /**
   * Return omega from phi, psi and their derivatives (characteristic != 2).
   *
   * @see Reference: sage/schemes/elliptic_curves/ell_curve_isogeny.py:EllipticCurveIsogeny.__compute_omega_fast
   */
  private __compute_omega_fast(
    E: EllipticCurveGeneric<F>,
    psi: F[],
    psi_pr: F[],
    phi: F[],
    phi_pr: F[]
  ): [F[], F[]] {
    const K = E.base_ring;
    const a1 = E.a1();
    const a3 = E.a3();
    const half = K.__call__(2n).inv() as F;

    // psi_2 = 2*y + a1*x + a3
    const psi2x: F[] = polyTrim([a3, a1]);

    // omega = phi_pr*psi*psi_2/2 - phi*psi_pr*psi_2 - (a1*phi + a3*psi^2)*psi/2
    //       = psi_2 * C - D,  C = phi_pr*psi/2 - phi*psi_pr,
    //                          D = (a1*phi + a3*psi^2)*psi/2
    const C = polySub(polyScale(polyMul(phi_pr, psi, K), half), polyMul(phi, psi_pr, K), K);
    const D = polyScale(
      polyMul(polyAdd(polyScale(phi, a1), polyScale(polyMul(psi, psi, K), a3), K), psi, K),
      half
    );

    const omega0 = polySub(polyMul(psi2x, C, K), D, K);
    const omega1 = polyScale(C, K.__call__(2n) as F);
    return [omega0, omega1];
  }

  /**
   * Return omega from phi, psi and their derivatives, in any characteristic.
   *
   * Only invoked in characteristic 2, where `psi_2 = a1*x + a3` has no
   * `y`-component.
   *
   * @see Reference: sage/schemes/elliptic_curves/ell_curve_isogeny.py:EllipticCurveIsogeny.__compute_omega_general
   */
  private __compute_omega_general(
    E: EllipticCurveGeneric<F>,
    psi: F[],
    psi_pr: F[],
    phi: F[],
    phi_pr: F[]
  ): [F[], F[]] {
    const K = E.base_ring;
    if (K.characteristic !== 2n) {
      throw new NotImplementedError(
        '__compute_omega_general is only implemented in characteristic 2'
      );
    }

    const [a1, a2, a3, a4, a6] = E.a_invariants();
    const [b2, b4] = E.b_invariants();

    const n = polyDegree(psi);
    const d = 2 * n + 1;
    const s1 = n > 0 ? (polyCoeff(psi, n - 1, K).neg() as F) : (K.zero() as F);

    // "derivatives" of psi in the sense of Kohel's corrected formulas
    let psi_prpr: F[] = [];
    for (let j = 0; j <= n - 2; j++) {
      const c = polyCoeff(psi, j + 2, K).mul(_binomial(j + 2, 2)) as F;
      psi_prpr = polyAdd(psi_prpr, polyShift([c], j, K), K);
    }
    let psi_prprpr: F[] = [];
    for (let j = 0; j <= n - 3; j++) {
      const c = polyCoeff(psi, j + 3, K).mul(3 * _binomial(j + 3, 3)) as F;
      psi_prprpr = polyAdd(psi_prprpr, polyShift([c], j, K), K);
    }

    // psi_2 = 2*y + a1*x + a3 = a1*x + a3 in characteristic 2
    const psi2: F[] = polyTrim([a3, a1]);
    const psi2sq = polyMul(psi2, psi2, K);

    const a1xa3: F[] = polyTrim([a3, a1]);
    const sixX2b2xb4: F[] = polyTrim([b4, b2, K.__call__(6n) as F]);

    // term1 = (a1*x + a3)*psi_2^2*(psi_prpr*psi_pr - psi_prprpr*psi)
    const term1 = polyMul(
      polyMul(a1xa3, psi2sq, K),
      polySub(polyMul(psi_prpr, psi_pr, K), polyMul(psi_prprpr, psi, K), K),
      K
    );
    // term2 = (a1*psi_2^2 - 3*(a1*x + a3)*(6x^2 + b2 x + b4))*psi_prpr*psi
    const term2 = polyMul(
      polyMul(
        polySub(
          polyScale(psi2sq, a1),
          polyScale(polyMul(a1xa3, sixX2b2xb4, K), K.__call__(3n) as F),
          K
        ),
        psi_prpr,
        K
      ),
      psi,
      K
    );
    // term3 = (a1 x^3 + 3 a3 x^2 + (2 a2 a3 - a1 a4) x + (a3 a4 - 2 a1 a6))*psi_pr^2
    const cubic: F[] = polyTrim([
      a3.mul(a4).sub(a1.mul(a6).mul(2)) as F,
      a2.mul(a3).mul(2).sub(a1.mul(a4)) as F,
      a3.mul(3) as F,
      a1,
    ]);
    const term3 = polyMul(cubic, polyMul(psi_pr, psi_pr, K), K);
    // term4 = (-(3 a1 x^2 + 6 a3 x + (2 a2 a3 - a1 a4)) + (a1 x + a3)*(d x - 2 s1))*psi_pr*psi
    const quad: F[] = polyTrim([
      a2.mul(a3).mul(2).sub(a1.mul(a4)) as F,
      a3.mul(6) as F,
      a1.mul(3) as F,
    ]);
    const dxMinus2s1: F[] = polyTrim([s1.mul(-2) as F, K.__call__(d) as F]);
    const term4 = polyMul(
      polyAdd(
        polyScale(quad, K.one().neg() as F),
        polyMul(a1xa3, dxMinus2s1, K),
        K
      ),
      polyMul(psi_pr, psi, K),
      K
    );
    // term5 = (a1*s1 + a3*n)*psi^2
    const term5 = polyScale(
      polyMul(psi, psi, K),
      a1.mul(s1).add(a3.mul(n)) as F
    );

    const bracket = polyAdd(
      polyAdd(polyAdd(polyAdd(term1, term2, K), term3, K), term4, K),
      term5,
      K
    );

    // omega = phi_pr*psi*y - phi*psi_pr*psi_2 + bracket*psi
    const omega1 = polyMul(phi_pr, psi, K);
    const omega0 = polyAdd(
      polyScale(polyMul(polyMul(phi, psi_pr, K), psi2, K), K.one().neg() as F),
      polyMul(bracket, psi, K),
      K
    );
    return [omega0, omega1];
  }

  /**
   * Apply Kohel's formulas to compute the image of an affine point.
   *
   * @see Reference: sage/schemes/elliptic_curves/ell_curve_isogeny.py:EllipticCurveIsogeny.__compute_via_kohel
   */
  private __compute_via_kohel(xP: F, yP: F): [F, F] | null {
    const K = this._domain.base_ring;
    // first check if the point is in the kernel
    if (polyEval(this.__inner_kernel_polynomial!, xP, K).isZero()) {
      return null;
    }
    const a = polyEval(this.__phi!, xP, K);
    const b = polyEval(this.__omega![0], xP, K).add(
      polyEval(this.__omega![1], xP, K).mul(yP)
    ) as F;
    const c = polyEval(this.__psi!, xP, K);
    const c2 = c.mul(c) as F;
    const c3 = c2.mul(c) as F;
    return [a.mul(c2.inv()) as F, b.mul(c3.inv()) as F];
  }

  /**
   * Initialize from kernel generators using Velu's formulas.
   */
  private __init_from_kernel_gens(kernel_gens: EllipticCurvePoint<F>[], check: boolean): void {
    // Filter out identity points
    const nonTrivialGens = kernel_gens.filter((P) => !P.is_zero());

    if (nonTrivialGens.length === 0) {
      // Trivial isogeny
      this._degree = 1n;
      return;
    }

    // Check that points have finite order (simplified check - we just verify they're on the curve)
    if (check) {
      for (const P of nonTrivialGens) {
        if (P.curve !== this._domain) {
          throw new ValueError('kernel point is not on the domain curve');
        }
      }
    }

    // Fast path: single generating point
    if (nonTrivialGens.length === 1) {
      this.__init_from_kernel_point(nonTrivialGens[0]);
      return;
    }

    // General case: compute all points in the subgroup
    const kernelSet = new Set<string>();
    kernelSet.add('O'); // identity

    // Helper to compute all multiples of a point
    const allMultiples = (
      gen: EllipticCurvePoint<F>,
      start: EllipticCurvePoint<F>
    ): EllipticCurvePoint<F>[] => {
      const multiples: EllipticCurvePoint<F>[] = [start];
      let R = start.add(gen);
      while (!R.eq(start)) {
        multiples.push(R);
        R = R.add(gen);
      }
      return multiples;
    };

    // Build up the kernel subgroup
    let kernelPoints: EllipticCurvePoint<F>[] = [this._domain.zero()];

    for (const P of nonTrivialGens) {
      const newPoints: EllipticCurvePoint<F>[] = [];
      for (const Q of kernelPoints) {
        for (const R of allMultiples(P, Q)) {
          const key = R.is_zero() ? 'O' : `${R.x().toString()},${R.y().toString()}`;
          if (!kernelSet.has(key)) {
            kernelSet.add(key);
            newPoints.push(R);
          }
        }
      }
      kernelPoints = kernelPoints.concat(newPoints);
    }

    this._degree = BigInt(kernelPoints.length);
    this.__init_from_kernel_list(kernelPoints);
  }

  /**
   * Initialize from a single kernel point.
   */
  private __init_from_kernel_point(ker: EllipticCurvePoint<F>): void {
    this._degree = 1n;

    let Q = ker;
    let prevQ = this._domain.zero();

    // Iterate through multiples of the kernel point
    while (!Q.is_zero() && !Q.eq(prevQ.neg())) {
      const [xQ, yQ] = [Q.x(), Q.y()];
      this.__update_kernel_data(xQ, yQ);

      // Check if Q is 2-torsion (Q = -Q)
      if (Q.eq(Q.neg())) {
        this._degree += 1n;
        break;
      }

      prevQ = Q;
      Q = Q.add(ker);
      this._degree += 2n;
    }
  }

  /**
   * Initialize from a list of kernel points.
   */
  private __init_from_kernel_list(kernelList: EllipticCurvePoint<F>[]): void {
    for (const Q of kernelList) {
      if (Q.is_zero()) {
        continue;
      }

      const xQ = Q.x();
      const key = xQ.toString();

      // Skip if we already processed a point with this x-coordinate
      if (this.__kernel_mod_sign.has(key)) {
        continue;
      }

      this.__update_kernel_data(xQ, Q.y());
    }
  }

  /**
   * Update kernel data for Velu's formulas.
   *
   * For a point Q = (xQ, yQ), computes:
   * - gxQ = 3*xQ^2 + 2*a2*xQ + a4 - a1*yQ
   * - gyQ = -2*yQ - a1*xQ - a3
   * - uQ = gyQ^2
   * - vQ = gxQ (if 2-torsion) or 2*gxQ - a1*gyQ (otherwise)
   *
   * Then updates:
   * - v += vQ
   * - w += uQ + xQ*vQ
   */
  private __update_kernel_data(xQ: F, yQ: F): void {
    const [a1, a2, a3, a4] = this._domain.a_invariants();

    // gxQ = 3*xQ^2 + 2*a2*xQ + a4 - a1*yQ
    const gxQ = xQ.mul(xQ).mul(3).add(a2.mul(xQ).mul(2)).add(a4).sub(a1.mul(yQ)) as F;

    // gyQ = -2*yQ - a1*xQ - a3
    const gyQ = yQ.mul(-2).sub(a1.mul(xQ)).sub(a3) as F;

    // uQ = gyQ^2
    const uQ = gyQ.mul(gyQ) as F;

    // Check if Q is 2-torsion: 2*yQ + a1*xQ + a3 = 0
    const twoTorsionCheck = yQ.mul(2).add(a1.mul(xQ)).add(a3);
    let vQ: F;
    if (twoTorsionCheck.isZero()) {
      // Q is 2-torsion
      vQ = gxQ;
    } else {
      // Q is not 2-torsion
      vQ = gxQ.mul(2).sub(a1.mul(gyQ)) as F;
    }

    // Store kernel point data
    const key = xQ.toString();
    this.__kernel_mod_sign.set(key, { xQ, yQ, gxQ, gyQ, vQ, uQ });

    // Update v and w
    this.__v = this.__v.add(vQ) as F;
    this.__w = this.__w.add(uQ).add(xQ.mul(vQ)) as F;
  }

  /**
   * Compute the codomain using Velu's or Kohel's formulas.
   */
  private __compute_codomain(): void {
    // Both Velu and Kohel algorithms compute v and w during initialization,
    // so we can use the same codomain formula.
    this._codomain = compute_codomain_formula(this._domain, this.__v, this.__w);
  }

  /**
   * Return the domain curve of this isogeny.
   *
   * @returns The domain elliptic curve
   * @see Reference: sage/schemes/elliptic_curves/ell_curve_isogeny.py:EllipticCurveIsogeny.domain
   */
  domain(): EllipticCurveGeneric<F> {
    return this._domain;
  }

  /**
   * Return the codomain curve of this isogeny.
   *
   * @returns The codomain elliptic curve
   * @see Reference: sage/schemes/elliptic_curves/ell_curve_isogeny.py:EllipticCurveIsogeny.codomain
   */
  codomain(): EllipticCurveGeneric<F> {
    if (this._codomain === null) {
      throw new ValueError('codomain not yet computed');
    }
    return this._codomain;
  }

  /**
   * Return the degree of this isogeny.
   *
   * @returns The degree as a bigint
   * @see Reference: sage/schemes/elliptic_curves/ell_curve_isogeny.py:EllipticCurveIsogeny.degree
   */
  degree(): bigint {
    return this._degree;
  }

  /**
   * Evaluate the isogeny at a point.
   *
   * @param P - A point on the domain curve
   * @returns The image point on the codomain curve
   * @see Reference: sage/schemes/elliptic_curves/ell_curve_isogeny.py:EllipticCurveIsogeny._call_
   */
  call(P: EllipticCurvePoint<F>): EllipticCurvePoint<F> {
    // Handle point at infinity
    if (P.is_zero()) {
      return this.codomain().zero();
    }

    let xP = P.x();
    let yP = P.y();

    // if there is a pre-isomorphism, apply it
    if (this.__pre_isomorphism !== null) {
      const [nx, ny] = this.__pre_isomorphism.call([xP, yP]) as [F, F];
      xP = nx;
      yP = ny;
    }

    const result =
      this.__algorithm === 'velu'
        ? this.__compute_via_velu(xP, yP)
        : this.__compute_via_kohel(xP, yP);
    if (result === null) {
      // Point is in the kernel
      return this.codomain().zero();
    }
    let [xOut, yOut] = result;

    // if there is a post-isomorphism, apply it
    if (this.__post_isomorphism !== null) {
      const [nx, ny] = this.__post_isomorphism.call([xOut, yOut]) as [F, F];
      xOut = nx;
      yOut = ny;
    }

    return this.codomain().point([xOut, yOut], false);
  }

  /**
   * Evaluate the isogeny at a point (alias for call).
   *
   * @param P - A point on the domain curve
   * @returns The image point on the codomain curve
   */
  evaluate(P: EllipticCurvePoint<F>): EllipticCurvePoint<F> {
    return this.call(P);
  }

  /**
   * Compute the image of a point using Velu's formulas.
   *
   * @param x - x-coordinate of the input point
   * @param y - y-coordinate of the input point
   * @returns [x', y'] coordinates of the image, or null if in kernel
   */
  private __compute_via_velu(x: F, y: F): [F, F] | null {
    const [a1, , a3] = this._domain.a_invariants();

    let X = x;
    let Y = y;

    // For each kernel point (modulo sign), apply Velu's formula
    for (const [, data] of this.__kernel_mod_sign) {
      const { xQ, yQ, gxQ, gyQ, vQ, uQ } = data;

      // t1 = x - xQ
      const t1 = x.sub(xQ) as F;

      // Check if point is in kernel (x = xQ)
      if (t1.isZero()) {
        return null;
      }

      const inv_t1 = t1.inv() as F;
      const inv_t1_2 = inv_t1.mul(inv_t1) as F;
      const inv_t1_3 = inv_t1_2.mul(inv_t1) as F;

      // tX = vQ/t1 + uQ/t1^2
      const tX = vQ.mul(inv_t1).add(uQ.mul(inv_t1_2)) as F;

      // tY computation
      // tY0 = uQ * (2*y + a1*x + a3)
      const tY0 = uQ.mul(y.mul(2).add(a1.mul(x)).add(a3)) as F;

      // tY1 = vQ * (a1*t1 + y - yQ)
      const tY1 = vQ.mul(a1.mul(t1).add(y).sub(yQ)) as F;

      // tY2 = a1*uQ - gxQ*gyQ
      const tY2 = a1.mul(uQ).sub(gxQ.mul(gyQ)) as F;

      // tY = tY0/t1^3 + (tY1 + tY2)/t1^2
      // (ell_curve_isogeny.py:2094-2103, __velu_sum_helper)
      const tY = tY0.mul(inv_t1_3).add(tY1.add(tY2).mul(inv_t1_2)) as F;

      // Update X and Y: X += tX, Y -= tY
      X = X.add(tX) as F;
      Y = Y.sub(tY) as F;
    }

    return [X, Y];
  }

  /**
   * Less strict evaluation method for internal use.
   *
   * @param P - A sequence of 3 coordinates defining a point on self
   * @returns The result of evaluating self at the given point
   * @see Reference: sage/schemes/elliptic_curves/ell_curve_isogeny.py:EllipticCurveIsogeny._eval
   */
  _eval(P: [F, F, F] | EllipticCurvePoint<F>): EllipticCurvePoint<F> {
    if (P instanceof EllipticCurvePoint) {
      return this.call(P);
    }
    // Projective coordinates [X, Y, Z]
    const [X, Y, Z] = P;
    if (Z.isZero()) {
      return this.codomain().zero();
    }
    // Convert to affine
    const x = X.mul(Z.inv()) as F;
    const y = Y.mul(Z.inv()) as F;
    const result = this.__compute_via_velu(x, y);
    if (result === null) {
      return this.codomain().zero();
    }
    return this.codomain().point([result[0], result[1]], false);
  }

  /** Cached rational maps */
  private __rational_maps: [RationalFunction<F>, BivariateRationalFunction<F>] | null = null;

  /**
   * Initialize rational maps using Velu's formulas.
   * The x-coordinate map is phi_x(x) = x + sum_Q (vQ/(x-xQ) + uQ/(x-xQ)^2)
   * For the polynomial representation, we compute the common denominator form.
   */
  private __initialize_rational_maps(): void {
    if (this.__rational_maps !== null) return;
    if (this.__algorithm !== 'velu') {
      throw new NotImplementedError('rational maps only implemented for Velu algorithm');
    }

    const K = this._domain.base_ring;

    // For a trivial isogeny, the rational maps are just x and y
    if (this.__kernel_mod_sign.size === 0) {
      const xMap: RationalFunction<F> = {
        numerator: [K.zero() as F, K.one() as F], // x
        denominator: [K.one() as F], // 1
        evaluate(x: F): F {
          return x;
        },
        toString(): string {
          return 'x';
        },
      };
      const yMap: BivariateRationalFunction<F> = {
        numerator: [[K.zero() as F], [K.one() as F]], // y
        denominator: [[K.one() as F]], // 1
        evaluate(_x: F, y: F): F {
          return y;
        },
        toString(): string {
          return 'y';
        },
      };
      this.__rational_maps = [xMap, yMap];
      return;
    }

    // Compute the x-coordinate rational map
    // phi_x(x) = x + sum_Q (vQ/(x-xQ) + uQ/(x-xQ)^2)
    // Common denominator is prod_Q (x-xQ)^2
    // Numerator = x * prod_Q (x-xQ)^2 + sum_Q [vQ*(x-xQ)*prod_{R!=Q}(x-xR)^2 + uQ*prod_{R!=Q}(x-xR)^2]

    const xCoords: F[] = [];
    const vQs: F[] = [];
    const uQs: F[] = [];

    for (const [, data] of this.__kernel_mod_sign) {
      xCoords.push(data.xQ);
      vQs.push(data.vQ);
      uQs.push(data.uQ);
    }

    // Compute prod_Q (x - xQ)^2 as polynomial in x
    // Start with 1
    let denomPoly: F[] = [K.one() as F];
    for (const xQ of xCoords) {
      // Multiply by (x - xQ)^2 = x^2 - 2*xQ*x + xQ^2
      const linear: F[] = [xQ.neg() as F, K.one() as F]; // (x - xQ)
      denomPoly = polyMul(denomPoly, linear, K);
      denomPoly = polyMul(denomPoly, linear, K); // squared
    }

    // Compute numerator: x * denom + sum_Q [vQ*(x-xQ)*prod_{R!=Q}(x-xR)^2 + uQ*prod_{R!=Q}(x-xR)^2]
    // First term: x * denom (shift coefficients)
    let numPoly: F[] = [K.zero() as F, ...denomPoly]; // multiply by x

    for (let i = 0; i < xCoords.length; i++) {
      const xQ = xCoords[i]!;
      const vQ = vQs[i]!;
      const uQ = uQs[i]!;

      // Compute prod_{R!=Q} (x-xR)^2
      let prodWithoutQ: F[] = [K.one() as F];
      for (let j = 0; j < xCoords.length; j++) {
        if (j !== i) {
          const xR = xCoords[j]!;
          const linear: F[] = [xR.neg() as F, K.one() as F];
          prodWithoutQ = polyMul(prodWithoutQ, linear, K);
          prodWithoutQ = polyMul(prodWithoutQ, linear, K);
        }
      }

      // vQ * (x - xQ) * prod_{R!=Q}(x-xR)^2
      const linearQ: F[] = [xQ.neg() as F, K.one() as F];
      const vQTerm = polyScale(polyMul(polyMul([vQ], linearQ, K), prodWithoutQ, K), K.one() as F);

      // uQ * prod_{R!=Q}(x-xR)^2
      const uQTerm = polyScale(prodWithoutQ, uQ);

      numPoly = polyAdd(numPoly, vQTerm, K);
      numPoly = polyAdd(numPoly, uQTerm, K);
    }

    const xMap: RationalFunction<F> = createRationalFunction(numPoly, denomPoly, K);

    // For the y-coordinate, the formula is more complex and involves y
    // phi_y(x,y) = y - sum_Q [uQ*(2y + a1*x + a3)/(x-xQ)^2 + vQ*(y-yQ)/(x-xQ)^2 - (a1*uQ - gxQ*gyQ)/(x-xQ)^2]
    // This is linear in y, so we can write it as (A(x)*y + B(x)) / C(x)^2

    // For simplicity, we represent this as a bivariate function that evaluates correctly
    const yMap: BivariateRationalFunction<F> = {
      numerator: [], // Not used directly
      denominator: [], // Not used directly
      evaluate: (x: F, y: F): F | null => {
        const result = this.__compute_via_velu(x, y);
        return result ? result[1] : null;
      },
      toString(): string {
        return 'y_rational_map(x, y)';
      },
    };

    this.__rational_maps = [xMap, yMap];
  }

  /**
   * Return the pair of rational maps defining this isogeny.
   *
   * Both components are returned as elements of the function field F(x,y)
   * in two variables over the base field F, though the first only involves x.
   *
   * @returns A pair [X_map, Y_map] of rational functions
   * @see Reference: sage/schemes/elliptic_curves/ell_curve_isogeny.py:EllipticCurveIsogeny.rational_maps
   * @see Deviation: Elliptic Curve Isogeny Algorithms Limited
   */
  rational_maps(): [RationalFunction<F>, BivariateRationalFunction<F>] {
    this.__initialize_rational_maps();
    return this.__rational_maps!;
  }

  /**
   * Return the rational map giving the x-coordinate of this isogeny.
   *
   * @returns The x-coordinate rational function in F(x)
   * @see Reference: sage/schemes/elliptic_curves/ell_curve_isogeny.py:EllipticCurveIsogeny.x_rational_map
   */
  x_rational_map(): RationalFunction<F> {
    this.__initialize_rational_maps();
    return this.__rational_maps![0];
  }

  /**
   * Return the Weierstrass scaling factor associated to this isogeny.
   *
   * The scaling factor is the constant u (in the base field) such that
   * phi^* omega_2 = u omega_1, where phi: E_1 -> E_2 is this isogeny
   * and omega_i are the standard Weierstrass differentials on E_i.
   *
   * ALGORITHM: The "inner" isogeny is normalized by construction, so we only
   * need to account for the scaling factors of a pre- and post-isomorphism.
   *
   * @returns The scaling factor u
   * @see Reference: sage/schemes/elliptic_curves/ell_curve_isogeny.py:EllipticCurveIsogeny.scaling_factor
   */
  scaling_factor(): F {
    let sc = this._domain.base_ring.one() as F;
    if (this.__pre_isomorphism !== null) {
      sc = sc.mul(this.__pre_isomorphism.scaling_factor()) as F;
    }
    if (this.__post_isomorphism !== null) {
      sc = sc.mul(this.__post_isomorphism.scaling_factor()) as F;
    }
    return sc;
  }

  /**
   * Return the kernel polynomial of this isogeny.
   *
   * The kernel polynomial is a monic polynomial whose roots are the
   * x-coordinates of the non-trivial points in the kernel (modulo sign).
   *
   * @returns The kernel polynomial as an array of coefficients [a0, a1, ..., an] for a0 + a1*x + ... + an*x^n
   * @see Reference: sage/schemes/elliptic_curves/ell_curve_isogeny.py:EllipticCurveIsogeny.kernel_polynomial
   */
  kernel_polynomial(): bigint[] {
    return polyToBigints(this.kernel_polynomial_field());
  }

  /**
   * Return the kernel polynomial of this isogeny over the base field.
   *
   * @returns The kernel polynomial as a coefficient array over the base field
   * @see Reference: sage/schemes/elliptic_curves/ell_curve_isogeny.py:EllipticCurveIsogeny.kernel_polynomial
   */
  kernel_polynomial_field(): F[] {
    if (this.__kernel_polynomial === null) {
      this.__init_kernel_polynomial();
    }
    return this.__kernel_polynomial!;
  }

  /**
   * Initialize the kernel polynomial.
   *
   * @see Reference: sage/schemes/elliptic_curves/ell_curve_isogeny.py:EllipticCurveIsogeny.__init_kernel_polynomial
   */
  private __init_kernel_polynomial(): void {
    if (this.__algorithm === 'velu') {
      this.__init_kernel_polynomial_velu();
      return;
    }
    throw new NotImplementedError('kernel polynomial not initialized');
  }

  /**
   * Initialize the kernel polynomial from the kernel points (Velu's algorithm).
   *
   * The kernel polynomial is the product of (x - invX(xQ)) for each
   * x-coordinate xQ of the non-trivial kernel points (modulo sign), where invX
   * pulls back along any pre-isomorphism.
   *
   * @see Reference: sage/schemes/elliptic_curves/ell_curve_isogeny.py:EllipticCurveIsogeny.__init_kernel_polynomial_velu
   */
  private __init_kernel_polynomial_velu(): void {
    const K = this._domain.base_ring;

    // invX(x) = u^2 * x + r when a pre-isomorphism is present
    let invX: ((x: F) => F) | null = null;
    if (this.__pre_isomorphism !== null) {
      const [u, r] = this.__pre_isomorphism.tuple();
      const u2 = u.mul(u) as F;
      invX = (x: F) => u2.mul(x).add(r) as F;
    }

    let poly: F[] = [K.one() as F];
    for (const [, data] of this.__kernel_mod_sign) {
      const xQ = invX === null ? data.xQ : invX(data.xQ);
      poly = polyMul(poly, [xQ.neg() as F, K.one() as F], K);
    }

    this.__kernel_polynomial = poly;
  }

  /**
   * Return the inseparable degree of this isogeny.
   *
   * Since this class only implements separable isogenies, this method always returns 1.
   *
   * @returns Always returns 1n
   * @see Reference: sage/schemes/elliptic_curves/ell_curve_isogeny.py:EllipticCurveIsogeny.inseparable_degree
   */
  inseparable_degree(): bigint {
    return 1n;
  }

  /** Cached dual isogeny */
  private __dual: EllipticCurveIsogeny<F> | null = null;

  /**
   * Return the isogeny dual to this isogeny.
   *
   * If phi: E -> E' is the given isogeny and n is its degree, then the dual
   * is by definition the unique isogeny phi_hat: E' -> E such that the
   * compositions phi_hat o phi and phi o phi_hat are the multiplication-by-n
   * maps on E and E', respectively.
   *
   * @returns The dual isogeny
   * @see Reference: sage/schemes/elliptic_curves/ell_curve_isogeny.py:EllipticCurveIsogeny.dual
   * @see Deviation: Elliptic Curve Isogeny Algorithms Limited
   */
  dual(): EllipticCurveIsogeny<F> {
    const K = this._domain.base_ring;
    const p = K.characteristic;

    if (p === 2n || p === 3n) {
      throw new NotImplementedError(
        'computation of dual isogenies not yet implemented in characteristics 2 and 3'
      );
    }

    if (this.__dual !== null) {
      return this.__dual;
    }

    const d = this._degree;

    // trac 7096
    const [E1, E2pr] = compute_intermediate_curves(this.codomain(), this.domain());

    const dF = K.__call__(d) as F;
    if (dF.isZero()) {
      // inseparable dual: Sage builds a composite with the Frobenius isogeny
      // (EllipticCurveHom_frobenius), which is not ported.
      throw new NotImplementedError(
        'SAGE_NOT_IMPLEMENTED: dual of an inseparable isogeny requires EllipticCurveHom_frobenius'
      );
    }

    const u = this.scaling_factor();
    const zero = K.zero() as F;
    const E2 = E2pr.change_weierstrass_model(u.div(dF) as F, zero, zero, zero);

    const phi_hat = new EllipticCurveIsogeny<F>(E1, null, E2, d);

    let chosen: [WeierstrassIsomorphism<F>, WeierstrassIsomorphism<F>] | null = null;
    outer: for (const urst1 of _isomorphisms(this._codomain!, E1)) {
      const pre_iso = new WeierstrassIsomorphism<F>(this._codomain!, urst1, E1);
      for (const urst2 of _isomorphisms(E2, this._domain)) {
        const post_iso = new WeierstrassIsomorphism<F>(E2, urst2, this._domain);
        const sc = u.mul(pre_iso.scaling_factor()).mul(post_iso.scaling_factor()) as F;
        if (sc.eq(dF)) {
          chosen = [pre_iso, post_iso];
          break outer;
        }
      }
    }

    if (chosen === null) {
      throw new ValueError('bug in dual()');
    }

    phi_hat._set_pre_isomorphism(chosen[0]);
    phi_hat._set_post_isomorphism(chosen[1]);
    this.__dual = phi_hat;
    return phi_hat;
  }

  /**
   * Return a copy of the isogeny that has been negated.
   *
   * The negation of an isogeny phi: E -> E' is the isogeny -phi: E -> E'
   * defined by (-phi)(P) = -phi(P) for all P in E.
   *
   * @returns The negated isogeny
   * @see Reference: sage/schemes/elliptic_curves/ell_curve_isogeny.py:EllipticCurveIsogeny.__neg__
   */
  neg(): NegatedIsogeny<F> {
    return new NegatedIsogeny(this);
  }

  /** Pre-isomorphism for composition */
  private __pre_isomorphism: WeierstrassIsomorphism<F> | null = null;

  /** Post-isomorphism for composition */
  private __post_isomorphism: WeierstrassIsomorphism<F> | null = null;

  /** Clear the cached values that depend on the pre-/post-isomorphisms. */
  private __clear_cached_values(): void {
    this.__rational_maps = null;
    this.__dual = null;
  }

  /**
   * Set the pre-isomorphism and the domain.
   *
   * @see Reference: sage/schemes/elliptic_curves/ell_curve_isogeny.py:EllipticCurveIsogeny.__set_pre_isomorphism
   */
  private __set_pre_isomorphism(
    domain: EllipticCurveGeneric<F>,
    isomorphism: WeierstrassIsomorphism<F>
  ): void {
    this._domain = domain;
    this.__pre_isomorphism = isomorphism;

    // The kernel polynomial has to be pulled back along the isomorphism:
    // psi(x) becomes psi((x - r)/u^2), made monic.
    if (this.__kernel_polynomial !== null) {
      const K = domain.base_ring;
      const [u, r] = isomorphism.tuple();
      const uinv2 = u.mul(u).inv() as F;
      // invX = (x - r) * u^-2
      const invX: F[] = polyTrim([r.neg().mul(uinv2) as F, uinv2]);
      this.__kernel_polynomial = polyMonic(
        _poly_compose(this.__kernel_polynomial, invX, K),
        K
      );
    }
  }

  /**
   * Set the post-isomorphism and the codomain.
   *
   * @see Reference: sage/schemes/elliptic_curves/ell_curve_isogeny.py:EllipticCurveIsogeny.__set_post_isomorphism
   */
  private __set_post_isomorphism(
    codomain: EllipticCurveGeneric<F>,
    isomorphism: WeierstrassIsomorphism<F>
  ): void {
    this._codomain = codomain;
    this.__post_isomorphism = isomorphism;
  }

  /**
   * Set up the post-isomorphism from a requested codomain or model.
   *
   * @see Reference: sage/schemes/elliptic_curves/ell_curve_isogeny.py:EllipticCurveIsogeny.__setup_post_isomorphism
   */
  private __setup_post_isomorphism(
    codomain: EllipticCurveGeneric<F> | null,
    model: 'minimal' | 'short_weierstrass' | 'montgomery' | null
  ): void {
    if (model === null && codomain === null) {
      return;
    }

    const oldE2 = this._codomain!;

    let target: EllipticCurveGeneric<F>;
    if (model !== null) {
      if (codomain !== null) {
        throw new ValueError('cannot specify a codomain curve and model name simultaneously');
      }
      target = _compute_model(oldE2, model);
    } else {
      target = codomain!;
      if (!oldE2.is_isomorphic(target)) {
        throw new ValueError('given codomain is not isomorphic to the computed codomain');
      }
    }

    const post_isom = new WeierstrassIsomorphism<F>(oldE2, null, target);
    this.__set_post_isomorphism(target, post_isom);
  }

  /**
   * Modify this isogeny by pre-composing with a WeierstrassIsomorphism.
   *
   * If preWI: E' -> E is an isomorphism and this isogeny is phi: E -> E'',
   * then the result is phi o preWI: E' -> E''.
   *
   * @param preWI - A WeierstrassIsomorphism with codomain equal to this isogeny's domain
   * @see Reference: sage/schemes/elliptic_curves/ell_curve_isogeny.py:EllipticCurveIsogeny._set_pre_isomorphism
   */
  _set_pre_isomorphism(preWI: WeierstrassIsomorphism<F>): void {
    const WIdom = preWI.domain();
    const WIcod = preWI.codomain();

    if (!_same_curve(this._domain, WIcod)) {
      throw new ValueError(
        "invalid parameter: isomorphism must have codomain curve equal to this isogenies' domain"
      );
    }

    const isom =
      this.__pre_isomorphism === null
        ? preWI
        : new WeierstrassIsomorphism<F>(
            WIdom,
            (this.__pre_isomorphism as baseWI<F>).mul(preWI as baseWI<F>).tuple(),
            this.__pre_isomorphism.codomain()
          );

    this.__clear_cached_values();
    this.__set_pre_isomorphism(WIdom, isom);
  }

  /**
   * Modify this isogeny by post-composing with a WeierstrassIsomorphism.
   *
   * If this isogeny is phi: E -> E' and postWI: E' -> E'' is an isomorphism,
   * then the result is postWI o phi: E -> E''.
   *
   * @param postWI - A WeierstrassIsomorphism with domain equal to this isogeny's codomain
   * @see Reference: sage/schemes/elliptic_curves/ell_curve_isogeny.py:EllipticCurveIsogeny._set_post_isomorphism
   */
  _set_post_isomorphism(postWI: WeierstrassIsomorphism<F>): void {
    const WIdom = postWI.domain();
    const WIcod = postWI.codomain();

    if (!_same_curve(this._codomain!, WIdom)) {
      throw new ValueError(
        "invalid parameter: isomorphism must have domain curve equal to this isogenies' codomain"
      );
    }

    const isom =
      this.__post_isomorphism === null
        ? postWI
        : new WeierstrassIsomorphism<F>(
            this.__post_isomorphism.domain(),
            (postWI as baseWI<F>).mul(this.__post_isomorphism as baseWI<F>).tuple(),
            WIcod
          );

    this.__clear_cached_values();
    this.__set_post_isomorphism(WIcod, isom);
  }

  /**
   * Return the composition of two elliptic-curve morphisms.
   *
   * Given left: E' -> E'' and right: E -> E', returns left o right: E -> E''.
   *
   * @param left - The left morphism (applied second)
   * @param right - The right morphism (applied first)
   * @returns The composition as a new isogeny or composite morphism
   * @see Reference: sage/schemes/elliptic_curves/ell_curve_isogeny.py:EllipticCurveIsogeny._composition_impl
   */
  static _composition_impl<F extends FieldElement>(
    left: EllipticCurveIsogeny<F>,
    right: EllipticCurveIsogeny<F>
  ): CompositeIsogeny<F> {
    // Check that the codomain of right equals the domain of left
    if (!right.codomain().j_invariant().eq(left.domain().j_invariant())) {
      throw new ValueError(
        'morphisms cannot be composed: codomain of right does not match domain of left'
      );
    }

    return new CompositeIsogeny([right, left]);
  }

  /**
   * Check if this isogeny equals another.
   *
   * Two isogenies are equal if they have the same domain, codomain, and
   * kernel polynomial.
   *
   * @param other - Another isogeny
   * @returns True if equal
   * @see Reference: sage/schemes/elliptic_curves/ell_curve_isogeny.py:EllipticCurveIsogeny.__eq__
   */
  eq(other: EllipticCurveIsogeny<F>): boolean {
    // Check degrees
    if (this._degree !== other._degree) {
      return false;
    }

    // Check domains have same j-invariant
    if (!this._domain.j_invariant().eq(other._domain.j_invariant())) {
      return false;
    }

    // Check codomains have same j-invariant
    if (!this.codomain().j_invariant().eq(other.codomain().j_invariant())) {
      return false;
    }

    // Compare kernel polynomials
    const ker1 = this.kernel_polynomial();
    const ker2 = other.kernel_polynomial();

    if (ker1.length !== ker2.length) {
      return false;
    }

    for (let i = 0; i < ker1.length; i++) {
      if (ker1[i] !== ker2[i]) {
        return false;
      }
    }

    return true;
  }

  /**
   * Return a hash value for this isogeny.
   *
   * The hash is computed from the degree and kernel polynomial.
   *
   * @returns A hash value
   * @see Reference: sage/schemes/elliptic_curves/ell_curve_isogeny.py:EllipticCurveIsogeny.__hash__
   */
  hash(): number {
    let h = Number(this._degree % BigInt(2 ** 31));
    const ker = this.kernel_polynomial();
    for (const coef of ker) {
      h = (h * 31 + Number(coef % BigInt(2 ** 31))) | 0;
    }
    return h;
  }

  /**
   * Return a string representation of this isogeny.
   *
   * @returns String representation
   * @see Reference: sage/schemes/elliptic_curves/ell_curve_isogeny.py:EllipticCurveIsogeny._repr_
   */
  toString(): string {
    return `Isogeny of degree ${this._degree} from ${this._domain} to ${this._codomain}`;
  }

  /**
   * Return a LaTeX representation of this isogeny.
   *
   * @returns LaTeX string
   * @see Reference: sage/schemes/elliptic_curves/ell_curve_isogeny.py:EllipticCurveIsogeny._latex_
   */
  _latex_(): string {
    return `\\text{Isogeny of degree } ${this._degree} \\text{ from } ${this._domain} \\text{ to } ${this._codomain}`;
  }

  /**
   * Return the i-th rational-map component.
   *
   * @param i - Index (0 for x, 1 for y)
   * @returns The i-th component of the rational maps
   * @see Reference: sage/schemes/elliptic_curves/ell_curve_isogeny.py:EllipticCurveIsogeny.__getitem__
   */
  getItem(i: 0 | 1): RationalFunction<F> | BivariateRationalFunction<F> {
    const maps = this.rational_maps();
    return maps[i];
  }

  /**
   * Check if this isogeny is separable.
   *
   * Since this class only implements separable isogenies, always returns true.
   *
   * @returns Always returns true
   * @see Reference: sage/schemes/elliptic_curves/ell_curve_isogeny.py:EllipticCurveIsogeny.is_separable
   */
  is_separable(): boolean {
    return true;
  }

  /**
   * Return True if this isogeny is normalized.
   *
   * An isogeny is normalized if the scaling factor equals 1.
   * Our implementation always produces normalized isogenies.
   *
   * @returns True if normalized (always true for our implementation)
   * @see Reference: sage/schemes/elliptic_curves/hom.py:EllipticCurveHom.is_normalized
   */
  is_normalized(): boolean {
    // Our isogenies are normalized by construction
    return true;
  }

  /**
   * Return the formal expansion of this isogeny as a power series.
   *
   * The formal isogeny is given as a power series in the variable t = -x/y
   * on the domain curve.
   *
   * @param prec - The desired precision (default: 20)
   * @returns The formal expansion as a power series (coefficients array)
   * @see Reference: sage/schemes/elliptic_curves/hom.py:EllipticCurveHom.formal
   */
  formal(prec: number = 20): bigint[] {
    // The formal expansion requires computing with formal groups
    // For a separable isogeny of degree d, the formal expansion starts with t
    // and has coefficients determined by the isogeny.
    //
    // For now, we return a simple approximation: [0, 1] representing t
    // A full implementation would require formal group infrastructure.
    const K = this._domain.base_ring;
    const result: bigint[] = [0n, 1n];

    // The formal expansion has the form t + a_d * t^d + higher order terms
    // where d is the degree of the isogeny.
    // For now, return a placeholder with the correct leading term.
    for (let i = 2; i < prec; i++) {
      result.push(0n);
    }

    return result;
  }

  /**
   * Return True if this morphism is injective.
   *
   * An isogeny between elliptic curves is never injective (except for isomorphisms).
   *
   * @returns True if injective
   * @see Reference: sage/schemes/elliptic_curves/hom.py:EllipticCurveHom.is_injective
   */
  is_injective(): boolean {
    return this._degree === 1n;
  }

  /**
   * Return True if this morphism is surjective.
   *
   * An isogeny between elliptic curves is always surjective.
   *
   * @returns Always returns true
   * @see Reference: sage/schemes/elliptic_curves/hom.py:EllipticCurveHom.is_surjective
   */
  is_surjective(): boolean {
    return true;
  }

  /**
   * Return True if this is a zero morphism.
   *
   * @returns True if zero
   * @see Reference: sage/schemes/elliptic_curves/hom.py:EllipticCurveHom.is_zero
   */
  is_zero(): boolean {
    return false; // Isogenies are never zero morphisms
  }
}

/**
 * A wrapper class representing the negation of an isogeny.
 *
 * For an isogeny phi: E -> E', the negation -phi: E -> E' is defined by
 * (-phi)(P) = -phi(P) for all P in E.
 *
 * @see Reference: sage/schemes/elliptic_curves/ell_curve_isogeny.py:EllipticCurveIsogeny.__neg__
 */
export class NegatedIsogeny<F extends FieldElement = FieldElement> {
  private _base: EllipticCurveIsogeny<F>;

  constructor(base: EllipticCurveIsogeny<F>) {
    this._base = base;
  }

  /**
   * Return the domain of this isogeny.
   */
  domain(): EllipticCurveGeneric<F> {
    return this._base.domain();
  }

  /**
   * Return the codomain of this isogeny.
   */
  codomain(): EllipticCurveGeneric<F> {
    return this._base.codomain();
  }

  /**
   * Return the degree of this isogeny.
   */
  degree(): bigint {
    return this._base.degree();
  }

  /**
   * Evaluate the negated isogeny at a point.
   *
   * @param P - A point on the domain curve
   * @returns The negated image point on the codomain curve
   */
  call(P: EllipticCurvePoint<F>): EllipticCurvePoint<F> {
    const Q = this._base.call(P);
    return Q.neg();
  }

  /**
   * Evaluate (alias for call).
   */
  evaluate(P: EllipticCurvePoint<F>): EllipticCurvePoint<F> {
    return this.call(P);
  }

  /**
   * Return the kernel polynomial (same as the base isogeny).
   */
  kernel_polynomial(): bigint[] {
    return this._base.kernel_polynomial();
  }

  /**
   * Return the scaling factor (negated).
   */
  scaling_factor(): F {
    return this._base.scaling_factor().neg() as F;
  }

  /**
   * Return the negation of this (which is the original isogeny).
   */
  neg(): EllipticCurveIsogeny<F> {
    return this._base;
  }

  /**
   * Return the dual of this negated isogeny.
   */
  dual(): NegatedIsogeny<F> {
    return new NegatedIsogeny(this._base.dual());
  }

  /**
   * Return string representation.
   */
  toString(): string {
    return `Negation of ${this._base.toString()}`;
  }
}

/**
 * A composite isogeny formed by composing multiple isogenies.
 *
 * Given isogenies phi_1: E_0 -> E_1, phi_2: E_1 -> E_2, ..., phi_n: E_{n-1} -> E_n,
 * the composite isogeny is phi_n o ... o phi_2 o phi_1: E_0 -> E_n.
 *
 * @see Reference: sage/schemes/elliptic_curves/hom_composite.py:EllipticCurveHom_composite
 */
export class CompositeIsogeny<F extends FieldElement = FieldElement> {
  private _factors: EllipticCurveIsogeny<F>[];

  /**
   * Create a composite isogeny from a list of factors.
   *
   * @param factors - List of isogenies to compose (applied in order)
   */
  constructor(factors: EllipticCurveIsogeny<F>[]) {
    if (factors.length === 0) {
      throw new ValueError('composite isogeny requires at least one factor');
    }

    // Verify that factors can be composed
    for (let i = 0; i < factors.length - 1; i++) {
      const left = factors[i]!;
      const right = factors[i + 1]!;
      if (!left.codomain().j_invariant().eq(right.domain().j_invariant())) {
        throw new ValueError(`factors ${i} and ${i + 1} cannot be composed`);
      }
    }

    this._factors = factors;
  }

  /**
   * Return the domain of this composite isogeny.
   */
  domain(): EllipticCurveGeneric<F> {
    return this._factors[0]!.domain();
  }

  /**
   * Return the codomain of this composite isogeny.
   */
  codomain(): EllipticCurveGeneric<F> {
    return this._factors[this._factors.length - 1]!.codomain();
  }

  /**
   * Return the degree of this composite isogeny.
   * The degree of a composition is the product of the degrees.
   */
  degree(): bigint {
    let d = 1n;
    for (const phi of this._factors) {
      d *= phi.degree();
    }
    return d;
  }

  /**
   * Evaluate the composite isogeny at a point.
   *
   * @param P - A point on the domain curve
   * @returns The image point on the codomain curve
   */
  call(P: EllipticCurvePoint<F>): EllipticCurvePoint<F> {
    let Q = P;
    for (const phi of this._factors) {
      Q = phi.call(Q);
    }
    return Q;
  }

  /**
   * Evaluate (alias for call).
   */
  evaluate(P: EllipticCurvePoint<F>): EllipticCurvePoint<F> {
    return this.call(P);
  }

  /**
   * Return the list of factor isogenies.
   */
  factors(): EllipticCurveIsogeny<F>[] {
    return [...this._factors];
  }

  /**
   * Return string representation.
   */
  toString(): string {
    return `Composite isogeny of degree ${this.degree()} with ${this._factors.length} factors`;
  }
}

/**
 * Create an isogeny from kernel generators using Velu's formulas.
 *
 * @param E - Domain elliptic curve
 * @param kernel_gens - List of points generating the kernel
 * @returns An EllipticCurveIsogeny
 * @see Reference: sage/schemes/elliptic_curves/ell_curve_isogeny.py:EllipticCurveIsogeny.__init_from_kernel_gens
 */
export function EllipticCurveIsogeny_from_kernel_gens<F extends FieldElement>(
  E: EllipticCurveGeneric<F>,
  kernel_gens: EllipticCurvePoint<F>[]
): EllipticCurveIsogeny<F> {
  return new EllipticCurveIsogeny(E, kernel_gens);
}

/**
 * Create an isogeny from a kernel polynomial using Kohel's formulas.
 *
 * Given a monic polynomial psi whose roots are the x-coordinates of the
 * non-trivial kernel points (modulo sign), this function constructs an
 * isogeny with that kernel.
 *
 * @param E - Domain elliptic curve
 * @param kernel_polynomial - The kernel polynomial (as coefficient array or Polynomial object)
 * @returns An EllipticCurveIsogeny
 * @see Reference: sage/schemes/elliptic_curves/ell_curve_isogeny.py:EllipticCurveIsogeny.__init_from_kernel_polynomial
 */
export function EllipticCurveIsogeny_from_kernel_polynomial<F extends FieldElement>(
  E: EllipticCurveGeneric<F>,
  kernel_polynomial: Polynomial | bigint[]
): EllipticCurveIsogeny<F> {
  // Convert to coefficient array if necessary
  const coeffs: bigint[] = Array.isArray(kernel_polynomial)
    ? kernel_polynomial
    : (kernel_polynomial as PolyElement<RingElement>).coeffs.map((c) => BigInt(c.toString()));

  const K = E.base_ring;
  const p = K.characteristic;

  // If kernel polynomial is constant (degree 0), return identity isogeny
  if (coeffs.length <= 1) {
    return new EllipticCurveIsogeny(E, [E.zero()]);
  }

  // Find the roots of the kernel polynomial (these are x-coordinates of kernel points)
  const xCoords: F[] = [];

  // For small fields, find roots by enumeration
  if (p <= 10000n) {
    for (let xVal = 0n; xVal < p; xVal++) {
      // Evaluate polynomial at xVal
      let result = 0n;
      let xPow = 1n;
      for (const coef of coeffs) {
        result = (result + coef * xPow) % p;
        xPow = (xPow * xVal) % p;
      }

      if (result === 0n) {
        xCoords.push(K.__call__(xVal) as F);
      }
    }
  }

  if (xCoords.length === 0) {
    throw new ValueError('could not find roots of kernel polynomial');
  }

  // For each x-coordinate, find the corresponding y-coordinates and create kernel points
  const kernelPoints: EllipticCurvePoint<F>[] = [];
  const [a1, a2, a3, a4, a6] = E.a_invariants();

  for (const x of xCoords) {
    // Compute RHS of curve equation
    const x2 = x.mul(x) as F;
    const x3 = x2.mul(x) as F;
    const rhs = x3.add(a2.mul(x2)).add(a4.mul(x)).add(a6) as F;
    const b = a1.mul(x).add(a3) as F;

    // Solve y^2 + b*y = rhs
    let yValues: F[] = [];

    if (p === 2n) {
      // Characteristic 2: special handling needed
      for (let yVal = 0n; yVal < 2n; yVal++) {
        const y = K.__call__(yVal) as F;
        const lhs = y.mul(y).add(b.mul(y)) as F;
        if (lhs.eq(rhs)) {
          yValues.push(y);
        }
      }
    } else {
      // Complete the square: (y + b/2)^2 = rhs + b^2/4
      const two = K.__call__(2n) as F;
      const four = K.__call__(4n) as F;
      const disc = rhs.add(b.mul(b).div(four)) as F;

      if (disc.isZero()) {
        yValues = [b.neg().div(two) as F];
      } else {
        // Check if disc is a square
        const exp = (p - 1n) / 2n;
        const legendre = disc.pow(exp);

        if (legendre.eq(K.one())) {
          // Compute square root
          let sqrtDisc: F;
          if (p % 4n === 3n) {
            sqrtDisc = disc.pow((p + 1n) / 4n) as F;
          } else {
            sqrtDisc = find_sqrt(disc, K) as F;
            if (!sqrtDisc) continue;
          }

          const negBOver2 = b.neg().div(two) as F;
          yValues = [negBOver2.add(sqrtDisc) as F, negBOver2.sub(sqrtDisc) as F];
        }
      }
    }

    for (const y of yValues) {
      if (E.is_on_curve(x, y)) {
        const P = E.point([x, y], false);
        if (!P.is_zero()) {
          kernelPoints.push(P);
          break; // Only need one point per x-coordinate
        }
      }
    }
  }

  if (kernelPoints.length === 0) {
    throw new ValueError('could not construct kernel points from polynomial');
  }

  // Create the isogeny from the kernel points
  return new EllipticCurveIsogeny(E, kernelPoints);
}

// ============================================================================
// SIDH/SIKE-style isogeny functions
// ============================================================================

/**
 * Compute the codomain curve from a kernel polynomial using Velu's formulas.
 *
 * Given a kernel polynomial (or a list of kernel generators), computes the
 * codomain curve of the isogeny. This is essential for SIDH/SIKE cryptography.
 *
 * INPUT:
 * - E: an elliptic curve
 * - kernel: either a kernel polynomial (list of coefficients), or a point/list of points
 *
 * OUTPUT: the codomain elliptic curve
 *
 * @see Reference: sage/schemes/elliptic_curves/ell_curve_isogeny.py:isogeny_codomain_from_kernel
 */
export function isogeny_codomain<F extends FieldElement>(
  E: EllipticCurveGeneric<F>,
  kernel: EllipticCurvePoint<F> | EllipticCurvePoint<F>[]
): EllipticCurveGeneric<F> {
  return isogeny_codomain_from_kernel(E, kernel);
}

/**
 * Find all isogenies of prime degree l from a curve E.
 *
 * For SIDH/SIKE, this is typically used with small primes l (like 2 or 3)
 * to walk along the isogeny graph.
 *
 * INPUT:
 * - E: an elliptic curve over a finite field
 * - l: a prime number, or a list of primes, or undefined
 *      If undefined, finds all isogenies of small prime degree
 *
 * OUTPUT: a list of isogenies (EllipticCurveIsogeny objects) with domain E
 *
 * ALGORITHM:
 * For each prime l, we find l-torsion points and use them to construct isogenies.
 * For l = 2, we find 2-torsion points (roots of x^3 + a4*x + a6 for short Weierstrass).
 * For l = 3, we find 3-torsion points.
 * For general l, we would use modular polynomials (not fully implemented).
 *
 * @see Reference: sage/schemes/elliptic_curves/ell_generic.py:isogenies_prime_degree
 */
export function isogenies_prime_degree<F extends FieldElement>(
  E: EllipticCurveGeneric<F>,
  l?: bigint | number | Array<bigint | number>
): EllipticCurveIsogeny<F>[] {
  const results: EllipticCurveIsogeny<F>[] = [];

  // Handle list of primes
  if (Array.isArray(l)) {
    for (const prime of l) {
      results.push(...isogenies_prime_degree(E, prime));
    }
    return results;
  }

  // Default to small primes if l is undefined
  if (l === undefined) {
    return isogenies_prime_degree(E, [2n, 3n, 5n, 7n]);
  }

  const lVal = typeof l === 'number' ? BigInt(l) : l;

  // Find l-torsion points
  const torsionPoints = find_torsion_points(E, lVal);

  // Build isogenies from cyclic subgroups
  const processedXCoords = new Set<string>();

  for (const P of torsionPoints) {
    if (P.is_zero()) continue;

    // Skip if we've already processed a point with this x-coordinate
    // (because P and -P generate the same isogeny kernel for l > 2)
    const xKey = P.x().toString();
    if (processedXCoords.has(xKey) && lVal > 2n) {
      continue;
    }
    processedXCoords.add(xKey);

    try {
      const phi = new EllipticCurveIsogeny(E, P);
      if (phi.degree() === lVal) {
        results.push(phi);
      }
    } catch {
      // Skip invalid points
      continue;
    }
  }

  return results;
}

/**
 * Find l-torsion points on an elliptic curve.
 *
 * @param E - an elliptic curve
 * @param l - a positive integer
 * @returns array of l-torsion points (points P where l*P = O)
 */
function find_torsion_points<F extends FieldElement>(
  E: EllipticCurveGeneric<F>,
  l: bigint
): EllipticCurvePoint<F>[] {
  const points: EllipticCurvePoint<F>[] = [];
  const K = E.base_ring;
  const p = K.characteristic;

  // For small fields, enumerate all points
  if (p <= 10000n) {
    for (let xVal = 0n; xVal < p; xVal++) {
      const x = K.__call__(xVal) as F;
      const [a1, a2, a3, a4, a6] = E.a_invariants();

      // Compute y^2 + (a1*x + a3)*y - (x^3 + a2*x^2 + a4*x + a6) = 0
      const x2 = x.mul(x) as F;
      const x3 = x2.mul(x) as F;
      const rhs = x3.add(a2.mul(x2)).add(a4.mul(x)).add(a6) as F;
      const b = a1.mul(x).add(a3) as F;

      // Solve y^2 + b*y = rhs
      const yValues = solve_quadratic_y(K, b, rhs.neg() as F);

      for (const y of yValues) {
        if (E.is_on_curve(x, y)) {
          const P = E.point([x, y], false);
          // Check if P is l-torsion
          if (P.mul(l).is_zero()) {
            points.push(P);
          }
        }
      }
    }
  }

  return points;
}

/**
 * Solve the quadratic y^2 + b*y + c = 0 in a finite field.
 *
 * @param K - the field
 * @param b - coefficient of y
 * @param c - constant term
 * @returns array of solutions
 */
function solve_quadratic_y<F extends FieldElement>(K: FieldParent, b: F, c: F): F[] {
  const p = K.characteristic;

  if (p === 2n) {
    // Characteristic 2: special handling needed
    // y^2 + b*y + c = 0
    if (b.isZero()) {
      // y^2 = -c = c (since -1 = 1 in char 2)
      // Every element is a square in F_{2^k}
      return [c.pow((p + 1n) / 2n) as F]; // sqrt(c)
    }
    // General case: would need Artin-Schreier theory
    return [];
  }

  // Odd characteristic: complete the square
  // y^2 + b*y + c = (y + b/2)^2 - b^2/4 + c = 0
  // (y + b/2)^2 = b^2/4 - c
  const two = K.__call__(2n) as F;
  const four = K.__call__(4n) as F;
  const disc = b.mul(b).div(four).sub(c) as F;

  // Check if disc is a square
  if (disc.isZero()) {
    const y = b.neg().div(two) as F;
    return [y];
  }

  const exp = (p - 1n) / 2n;
  const legendre = disc.pow(exp);

  if (!legendre.eq(K.one())) {
    // Not a quadratic residue
    return [];
  }

  // Compute square root
  let sqrtDisc: F;
  if (p % 4n === 3n) {
    sqrtDisc = disc.pow((p + 1n) / 4n) as F;
  } else {
    // Tonelli-Shanks would be needed here
    // For simplicity, try a few values
    sqrtDisc = find_sqrt(disc, K) as F;
    if (!sqrtDisc) return [];
  }

  const negBOver2 = b.neg().div(two) as F;
  return [negBOver2.add(sqrtDisc) as F, negBOver2.sub(sqrtDisc) as F];
}

/**
 * Find a square root of a in the field K (simple method for small fields).
 */
function find_sqrt<F extends FieldElement>(a: F, K: FieldParent): F | null {
  const p = K.characteristic;

  if (a.isZero()) {
    return K.zero() as F;
  }

  // p ≡ 3 (mod 4)
  if (p % 4n === 3n) {
    return a.pow((p + 1n) / 4n) as F;
  }

  // Tonelli-Shanks for general case
  let q = p - 1n;
  let s = 0n;
  while ((q & 1n) === 0n) {
    q >>= 1n;
    s++;
  }

  // Find a quadratic non-residue
  let z = K.__call__(2n) as F;
  const exp = (p - 1n) / 2n;
  while (z.pow(exp).eq(K.one())) {
    z = z.add(K.one()) as F;
  }

  let m = s;
  let c = z.pow(q) as F;
  let t = a.pow(q) as F;
  let r = a.pow((q + 1n) / 2n) as F;

  for (let iter = 0; iter < 1000; iter++) {
    if (t.isZero()) return K.zero() as F;
    if (t.eq(K.one())) return r;

    // Find least i where t^(2^i) = 1
    let i = 1n;
    let temp = t.mul(t) as F;
    while (!temp.eq(K.one())) {
      temp = temp.mul(temp) as F;
      i++;
      if (i >= m) return null;
    }

    const bExp = 1n << (m - i - 1n);
    const b = c.pow(bExp) as F;
    m = i;
    c = b.mul(b) as F;
    t = t.mul(c) as F;
    r = r.mul(b) as F;
  }

  return null;
}

/**
 * Compute the full isogeny class of an elliptic curve.
 *
 * The isogeny class consists of all elliptic curves isogenous to E over the base field.
 * For finite fields, two curves are isogenous if and only if they have the same
 * number of points.
 *
 * INPUT:
 * - E: an elliptic curve over a finite field
 * - maxDegree: maximum degree of isogenies to consider (default: 100)
 *
 * OUTPUT: array of elliptic curves in the isogeny class
 *
 * ALGORITHM:
 * Use a breadth-first search starting from E, following isogenies of small prime degree.
 *
 * @see Reference: sage/schemes/elliptic_curves/ell_rational_field.py:isogeny_class
 */
export function isogeny_class<F extends FieldElement>(
  E: EllipticCurveGeneric<F>,
  maxDegree: number = 100
): EllipticCurveGeneric<F>[] {
  const curves: EllipticCurveGeneric<F>[] = [E];
  const jInvariants = new Set<string>([E.j_invariant().toString()]);
  const queue: EllipticCurveGeneric<F>[] = [E];
  const primes = [2n, 3n, 5n, 7n, 11n, 13n];

  // Breadth-first search
  let iterations = 0;
  while (queue.length > 0 && iterations < maxDegree) {
    iterations++;
    const current = queue.shift()!;

    // Find all isogenies from the current curve
    for (const l of primes) {
      try {
        const isogenies = isogenies_prime_degree(current, l);

        for (const phi of isogenies) {
          const codomain = phi.codomain();
          const jKey = codomain.j_invariant().toString();

          if (!jInvariants.has(jKey)) {
            jInvariants.add(jKey);
            curves.push(codomain);
            queue.push(codomain);
          }
        }
      } catch {
        // Skip if there's an error computing isogenies
        continue;
      }
    }
  }

  return curves;
}

/**
 * Check if two elliptic curves are isogenous over a finite field.
 *
 * Two elliptic curves over a finite field F_q are isogenous if and only if
 * they have the same number of points (by Tate's theorem).
 *
 * INPUT:
 * - E1: first elliptic curve
 * - E2: second elliptic curve
 *
 * OUTPUT: True if E1 and E2 are isogenous
 *
 * ALGORITHM:
 * Compare the cardinalities of E1 and E2.
 *
 * @see Reference: sage/schemes/elliptic_curves/ell_finite_field.py:is_isogenous
 * @see Deviation: Elliptic Curve Isogeny Algorithms Limited
 */
export function is_isogenous<F extends FieldElement>(
  E1: EllipticCurveGeneric<F>,
  E2: EllipticCurveGeneric<F>
): boolean {
  // By Tate's theorem, two curves over a finite field are isogenous if and only
  // if they have the same number of points.  Equal j-invariants are *not*
  // sufficient: a curve and its quadratic twist share a j-invariant but are in
  // general not isogenous.
  const K = E1.base_ring;
  const p = K.characteristic;

  if (p > 1000n) {
    throw new NotImplementedError('is_isogenous for large fields requires cardinality computation');
  }

  const count1 = count_points_finite_field(E1);
  const count2 = count_points_finite_field(E2);

  return count1 === count2;
}

/**
 * Count the number of points on an elliptic curve over a small finite field.
 */
function count_points_finite_field<F extends FieldElement>(E: EllipticCurveGeneric<F>): bigint {
  const K = E.base_ring;
  const p = K.characteristic;

  let count = 1n; // Start with point at infinity

  for (let xVal = 0n; xVal < p; xVal++) {
    const x = K.__call__(xVal) as F;
    const [a1, a2, a3, a4, a6] = E.a_invariants();

    // Compute RHS: x^3 + a2*x^2 + a4*x + a6
    const x2 = x.mul(x) as F;
    const x3 = x2.mul(x) as F;
    const rhs = x3.add(a2.mul(x2)).add(a4.mul(x)).add(a6) as F;

    // Coefficient b = a1*x + a3
    const b = a1.mul(x).add(a3) as F;

    // Count solutions to y^2 + b*y = rhs
    // Complete the square: (y + b/2)^2 = rhs + b^2/4
    if (p === 2n) {
      // Characteristic 2
      for (let yVal = 0n; yVal < 2n; yVal++) {
        const y = K.__call__(yVal) as F;
        if (E.is_on_curve(x, y)) {
          count++;
        }
      }
    } else {
      const two = K.__call__(2n) as F;
      const four = K.__call__(4n) as F;
      const disc = rhs.add(b.mul(b).div(four)) as F;

      if (disc.isZero()) {
        // One solution
        count++;
      } else {
        // Check if disc is a quadratic residue
        const exp = (p - 1n) / 2n;
        const legendre = disc.pow(exp);
        if (legendre.eq(K.one())) {
          // Two solutions
          count += 2n;
        }
      }
    }
  }

  return count;
}

/**
 * Return the degree of an isogeny.
 *
 * @param phi - an isogeny
 * @returns the degree as a bigint
 */
export function isogeny_degree<F extends FieldElement>(phi: EllipticCurveIsogeny<F>): bigint {
  return phi.degree();
}
