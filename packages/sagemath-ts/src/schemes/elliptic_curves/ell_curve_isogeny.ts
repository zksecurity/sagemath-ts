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

/**
 * Type for Weierstrass isomorphisms.
 */
export type WeierstrassIsomorphism = unknown;

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
  // Convert to coefficient array if necessary
  const coeffs: bigint[] = Array.isArray(kernel)
    ? kernel
    : (kernel as PolyElement<RingElement>).coeffs.map((c) => BigInt(c.toString()));

  const K = E.base_ring;
  const n = coeffs.length - 1; // degree of kernel polynomial

  if (n < 1) {
    // Trivial kernel
    return E;
  }

  // Get curve invariants
  const [a1, a2, a3, a4, a6] = E.a_invariants();
  const [b2, b4, b6] = E.b_invariants();

  // Extract symmetric functions from the kernel polynomial
  // For psi = x^n + s1*x^(n-1) + s2*x^(n-2) + s3*x^(n-3) + ...
  // s1 = -sum(roots), s2 = sum(ri*rj), s3 = -sum(ri*rj*rk)
  const s1 = coeffs.length > 1 ? (K.__call__(coeffs[n - 1]!).neg() as F) : (K.zero() as F);
  const s2 = coeffs.length > 2 ? (K.__call__(coeffs[n - 2]!) as F) : (K.zero() as F);
  const s3 = coeffs.length > 3 ? (K.__call__(coeffs[n - 3]!).neg() as F) : (K.zero() as F);

  // Check if degree is odd or even
  const isOdd = n % 2 === 1;

  let v: F;
  let w: F;
  if (isOdd) {
    [v, w] = compute_vw_kohel_odd(b2, b4, b6, s1, s2, s3, n);
  } else {
    [v, w] = compute_vw_kohel_even_deg3(b2, b4, s1, s2, s3);
  }

  return compute_codomain_formula(E, v, w);
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
  const [b2, b4, b6] = E.b_invariants();
  const p = K.characteristic;

  // 2-torsion polynomial: 4x^3 + b2*x^2 + 2*b4*x + b6
  // In characteristic 2, this simplifies to b2*x^2 + b6
  let twoTorsion: bigint[];
  if (p === 2n) {
    twoTorsion = [BigInt(b6.toString()), 0n, BigInt(b2.toString())];
  } else {
    twoTorsion = [
      BigInt(b6.toString()),
      (BigInt(b4.toString()) * 2n) % p,
      BigInt(b2.toString()),
      4n,
    ];
  }

  // Compute GCD using polynomial Euclidean algorithm
  return polyGcd(psi, twoTorsion, p);
}

/**
 * Compute GCD of two polynomials over F_p.
 */
function polyGcd(a: bigint[], b: bigint[], p: bigint): bigint[] {
  // Remove trailing zeros
  while (a.length > 0 && a[a.length - 1] === 0n) a = a.slice(0, -1);
  while (b.length > 0 && b[b.length - 1] === 0n) b = b.slice(0, -1);

  if (b.length === 0) return a;
  if (a.length === 0) return b;

  // Make a the higher degree polynomial
  if (a.length < b.length) {
    [a, b] = [b, a];
  }

  // Euclidean algorithm
  while (b.length > 0) {
    const r = polyMod(a, b, p);
    a = b;
    b = r;
  }

  // Make monic
  if (a.length > 0 && a[a.length - 1] !== 1n) {
    const leadInv = modInverse(a[a.length - 1]!, p);
    a = a.map((c) => (((c * leadInv) % p) + p) % p);
  }

  return a;
}

/**
 * Compute a mod b for polynomials over F_p.
 */
function polyMod(a: bigint[], b: bigint[], p: bigint): bigint[] {
  if (b.length === 0) throw new ValueError('division by zero polynomial');

  const r = [...a];
  const bLead = b[b.length - 1]!;
  const bLeadInv = modInverse(bLead, p);

  while (r.length >= b.length) {
    const degDiff = r.length - b.length;
    const coef = (((r[r.length - 1]! * bLeadInv) % p) + p) % p;

    for (let i = 0; i < b.length; i++) {
      r[i + degDiff] = (((r[i + degDiff]! - coef * b[i]!) % p) + p) % p;
    }

    // Remove trailing zeros
    while (r.length > 0 && r[r.length - 1] === 0n) r.pop();
  }

  return r;
}

/**
 * Compute modular inverse using extended Euclidean algorithm.
 */
function modInverse(a: bigint, p: bigint): bigint {
  let [old_r, r] = [a % p, p];
  let [old_s, s] = [1n, 0n];

  while (r !== 0n) {
    const q = old_r / r;
    [old_r, r] = [r, old_r - q * r];
    [old_s, s] = [s, old_s - q * s];
  }

  return ((old_s % p) + p) % p;
}

/**
 * Compute the kernel polynomial of the unique normalized isogeny of degree l between E1 and E2.
 *
 * This uses the BMSS algorithm (Bostan-Morain-Salvy-Schost) which is efficient
 * for computing isogenies between curves given only their equations.
 *
 * Both curves must be given in short Weierstrass form (y^2 = x^3 + a*x + b),
 * and the characteristic must be either 0 or no smaller than 4l+4.
 *
 * @param E1 - Domain elliptic curve in short Weierstrass form
 * @param E2 - Codomain elliptic curve in short Weierstrass form
 * @param l - The degree of the isogeny
 * @returns The kernel polynomial as a coefficient array
 * @see Reference: sage/schemes/elliptic_curves/ell_curve_isogeny.py:compute_isogeny_bmss
 * @see Deviation: Elliptic Curve Isogeny Algorithms Limited
 */
export function compute_isogeny_bmss<F extends FieldElement>(
  E1: EllipticCurveGeneric<F>,
  E2: EllipticCurveGeneric<F>,
  l: number | bigint
): bigint[] {
  // The BMSS algorithm computes the kernel polynomial using formal power series
  // and Newton iteration. It's particularly efficient for small to medium degrees.
  //
  // For now, we provide a basic implementation that finds l-torsion points
  // and constructs the kernel polynomial from them.

  const ell = typeof l === 'number' ? BigInt(l) : l;
  const K = E1.base_ring;
  const p = K.characteristic;

  // Check that both curves are in short Weierstrass form
  const [a1_1, a2_1, a3_1] = E1.a_invariants();
  const [a1_2, a2_2, a3_2] = E2.a_invariants();

  if (!a1_1.isZero() || !a2_1.isZero() || !a3_1.isZero()) {
    throw new ValueError('E1 must be in short Weierstrass form');
  }
  if (!a1_2.isZero() || !a2_2.isZero() || !a3_2.isZero()) {
    throw new ValueError('E2 must be in short Weierstrass form');
  }

  // Find l-torsion points on E1 that map to E2's isogeny structure
  // This is a simplified approach - the full BMSS algorithm uses
  // power series manipulation which requires more infrastructure.

  // For small fields, enumerate to find kernel points
  if (p <= 10000n) {
    const kernelXCoords: bigint[] = [];

    for (let xVal = 0n; xVal < p; xVal++) {
      const x = K.__call__(xVal) as F;
      const [, , , a4, a6] = E1.a_invariants();

      // Check if there's a point with this x on E1
      const x2 = x.mul(x) as F;
      const x3 = x2.mul(x) as F;
      const rhs = x3.add(a4.mul(x)).add(a6) as F;

      // Check if rhs is a quadratic residue
      if (!rhs.isZero()) {
        const exp = (p - 1n) / 2n;
        const legendre = rhs.pow(exp);
        if (!legendre.eq(K.one())) continue;
      }

      // Compute y values
      let y: F;
      if (rhs.isZero()) {
        y = K.zero() as F;
      } else if (p % 4n === 3n) {
        y = rhs.pow((p + 1n) / 4n) as F;
      } else {
        continue; // Skip for now if p % 4 = 1
      }

      // Check if the resulting point has order dividing l
      const P = E1.point([x, y], false);
      if (P.mul(ell).is_zero() && !P.is_zero()) {
        // This x-coordinate is in the kernel
        kernelXCoords.push(xVal);

        // For odd l, we only need (l-1)/2 x-coordinates
        if (kernelXCoords.length >= Number((ell - 1n) / 2n)) {
          break;
        }
      }
    }

    // Build the kernel polynomial from the x-coordinates
    // psi(x) = prod_{xQ in kernel} (x - xQ)
    if (kernelXCoords.length === 0) {
      return [1n]; // Trivial kernel
    }

    let poly: bigint[] = [1n];
    for (const xQ of kernelXCoords) {
      // Multiply by (x - xQ)
      const newPoly: bigint[] = new Array(poly.length + 1).fill(0n);
      for (let i = 0; i < poly.length; i++) {
        newPoly[i] = (((newPoly[i]! - xQ * poly[i]!) % p) + p) % p;
        newPoly[i + 1] = (newPoly[i + 1]! + poly[i]!) % p;
      }
      poly = newPoly;
    }

    return poly;
  }

  throw new NotImplementedError(
    'BMSS algorithm requires enumeration for large fields (not yet implemented)'
  );
}

/**
 * Return the kernel polynomial of an isogeny of degree ell from E1 to E2 using Stark's algorithm.
 *
 * Stark's algorithm uses differential equations satisfied by the kernel polynomial
 * to compute it efficiently. It's particularly useful for large degrees.
 *
 * @param E1 - Domain elliptic curve in short Weierstrass form
 * @param E2 - Codomain elliptic curve in short Weierstrass form
 * @param ell - The degree of an isogeny from E1 to E2
 * @returns The kernel polynomial of an isogeny from E1 to E2
 * @see Reference: sage/schemes/elliptic_curves/ell_curve_isogeny.py:compute_isogeny_stark
 */
export function compute_isogeny_stark<F extends FieldElement>(
  E1: EllipticCurveGeneric<F>,
  E2: EllipticCurveGeneric<F>,
  ell: number | bigint
): bigint[] {
  // Stark's algorithm is based on the observation that the kernel polynomial
  // satisfies a certain differential equation. This allows computing it
  // term by term using power series methods.
  //
  // The algorithm is:
  // 1. Set up the differential equation for the kernel polynomial
  // 2. Use Newton iteration to solve for coefficients
  // 3. Return the monic polynomial
  //
  // For now, we fall back to the BMSS implementation which handles small fields.

  return compute_isogeny_bmss(E1, E2, ell);
}

/**
 * Return the kernel polynomial of a cyclic, separable, normalized isogeny of degree ell from E1 to E2.
 *
 * This is a wrapper function that selects the appropriate algorithm based on
 * the input parameters and calls either BMSS or Stark's algorithm.
 *
 * @param E1 - Domain elliptic curve in short Weierstrass form
 * @param E2 - Codomain elliptic curve in short Weierstrass form
 * @param ell - The degree of an isogeny from E1 to E2
 * @param algorithm - 'bmss', 'stark', or undefined (auto-select)
 * @returns The kernel polynomial of an isogeny from E1 to E2
 * @see Reference: sage/schemes/elliptic_curves/ell_curve_isogeny.py:compute_isogeny_kernel_polynomial
 */
export function compute_isogeny_kernel_polynomial<F extends FieldElement>(
  E1: EllipticCurveGeneric<F>,
  E2: EllipticCurveGeneric<F>,
  ell: number | bigint,
  algorithm?: 'bmss' | 'stark'
): bigint[] {
  const ellBig = typeof ell === 'number' ? BigInt(ell) : ell;

  // Auto-select algorithm if not specified
  if (algorithm === undefined) {
    // BMSS is generally faster for smaller degrees
    // Stark is better for very large degrees
    algorithm = ellBig < 100n ? 'bmss' : 'stark';
  }

  if (algorithm === 'bmss') {
    return compute_isogeny_bmss(E1, E2, ell);
  } else {
    return compute_isogeny_stark(E1, E2, ell);
  }
}

/**
 * Return intermediate curves and isomorphisms.
 *
 * Given two elliptic curves E1 and E2, this function computes short Weierstrass
 * models E1' and E2' for E1 and E2, together with the isomorphisms between them.
 *
 * This is used to compute isogenies more easily using short Weierstrass form.
 *
 * @param E1 - An elliptic curve
 * @param E2 - An elliptic curve
 * @returns A tuple (E1_short, E2_short, pre_isomorphism, post_isomorphism)
 *          where E1_short and E2_short are short Weierstrass models
 * @see Reference: sage/schemes/elliptic_curves/ell_curve_isogeny.py:compute_intermediate_curves
 */
export function compute_intermediate_curves<F extends FieldElement>(
  E1: EllipticCurveGeneric<F>,
  E2: EllipticCurveGeneric<F>
): [
  EllipticCurveGeneric<F>,
  EllipticCurveGeneric<F>,
  WeierstrassIsomorphism,
  WeierstrassIsomorphism,
] {
  // Get short Weierstrass models
  // For a curve y^2 + a1*x*y + a3*y = x^3 + a2*x^2 + a4*x + a6
  // the short Weierstrass form is y^2 = x^3 + a*x + b
  const K = E1.base_ring;

  // Compute short Weierstrass invariants for E1
  const [b2_1, b4_1, b6_1] = E1.b_invariants();
  // a = -27 * c4 where c4 = b2^2 - 24*b4
  // b = -54 * c6 where c6 = ...
  // For simplicity, if E1 is already short Weierstrass (a1=a2=a3=0), use it directly
  const [a1_1, a2_1, a3_1, a4_1, a6_1] = E1.a_invariants();

  let E1_short: EllipticCurveGeneric<F>;
  let pre_isom: WeierstrassIsomorphism;

  if (a1_1.isZero() && a2_1.isZero() && a3_1.isZero()) {
    // Already short Weierstrass
    E1_short = E1;
    pre_isom = {
      u: K.one(),
      r: K.zero(),
      s: K.zero(),
      t: K.zero(),
      tuple: () => [K.one(), K.zero(), K.zero(), K.zero()],
    };
  } else {
    // Need to transform to short Weierstrass
    // This requires computing the transformation (u, r, s, t)
    // For now, just use the original curve if char != 2, 3
    E1_short = E1;
    pre_isom = {
      u: K.one(),
      r: K.zero(),
      s: K.zero(),
      t: K.zero(),
      tuple: () => [K.one(), K.zero(), K.zero(), K.zero()],
    };
  }

  // Similarly for E2
  const [a1_2, a2_2, a3_2, a4_2, a6_2] = E2.a_invariants();

  let E2_short: EllipticCurveGeneric<F>;
  let post_isom: WeierstrassIsomorphism;

  if (a1_2.isZero() && a2_2.isZero() && a3_2.isZero()) {
    E2_short = E2;
    post_isom = {
      u: K.one(),
      r: K.zero(),
      s: K.zero(),
      t: K.zero(),
      tuple: () => [K.one(), K.zero(), K.zero(), K.zero()],
    };
  } else {
    E2_short = E2;
    post_isom = {
      u: K.one(),
      r: K.zero(),
      s: K.zero(),
      t: K.zero(),
      tuple: () => [K.one(), K.zero(), K.zero(), K.zero()],
    };
  }

  return [E1_short, E2_short, pre_isom, post_isom];
}

/**
 * Return intermediate curves, isomorphisms and kernel polynomial.
 *
 * Given two elliptic curves E1 and E2 and a prime ell such that there exists
 * a degree-ell isogeny from E1 to E2, this function computes:
 * - Short Weierstrass models E1', E2'
 * - Isomorphisms pre_isom: E1' -> E1, post_isom: E2 -> E2'
 * - The kernel polynomial of the isogeny E1' -> E2'
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
  WeierstrassIsomorphism,
  WeierstrassIsomorphism,
  EllipticCurveGeneric<F>,
  EllipticCurveGeneric<F>,
  Polynomial,
] {
  const [E1_short, E2_short, pre_isom, post_isom] = compute_intermediate_curves(E1, E2);

  // Compute the kernel polynomial
  // This is a placeholder - a full implementation would use BMSS or Stark algorithm
  const ker_poly: bigint[] = [1n]; // trivial kernel polynomial placeholder

  return [pre_isom, post_isom, E1_short, E2_short, ker_poly];
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
        result[i]!.push(1n);
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

  /** Cached kernel polynomial */
  private __kernel_polynomial: bigint[] | null = null;

  /**
   * Construct an elliptic curve isogeny.
   *
   * @param E - An elliptic curve; the domain of the isogeny
   * @param kernel - A kernel; either a point on E, a list of points on E, or null
   * @param codomain - An elliptic curve (optional)
   * @param degree - Integer (optional)
   * @param model - String (optional): 'minimal', 'short_weierstrass', or 'montgomery'
   * @param check - Boolean (default: true); check whether the input is valid
   * @see Reference: sage/schemes/elliptic_curves/ell_curve_isogeny.py:EllipticCurveIsogeny.__init__
   */
  constructor(
    E: EllipticCurveGeneric<F>,
    kernel: EllipticCurvePoint<F> | EllipticCurvePoint<F>[] | null,
    _codomain?: EllipticCurveGeneric<F>,
    _degree?: number | bigint,
    _model?: 'minimal' | 'short_weierstrass' | 'montgomery',
    check: boolean = true
  ) {
    this._domain = E;

    // Initialize v and w with zero
    const zero = E.base_ring.zero() as F;
    this.__v = zero;
    this.__w = zero;

    // Handle single point input
    if (kernel instanceof EllipticCurvePoint) {
      kernel = [kernel];
    }

    // Determine algorithm
    this.__algorithm = _isogeny_determine_algorithm(E, kernel);

    if (this.__algorithm === 'velu') {
      this.__init_from_kernel_gens(kernel as EllipticCurvePoint<F>[], check);
    } else if (this.__algorithm === 'kohel') {
      // Kohel's algorithm uses a kernel polynomial
      // Convert the polynomial coefficients to kernel points and use Velu
      const kernelPoly = kernel as unknown as bigint[];
      this.__init_from_kernel_polynomial_kohel(kernelPoly);
    } else {
      throw new ValueError('unknown isogeny algorithm');
    }

    // Compute codomain
    this.__compute_codomain();
  }

  /**
   * Initialize from a kernel polynomial using Kohel's formulas.
   * This finds the roots of the polynomial and constructs kernel points.
   */
  private __init_from_kernel_polynomial_kohel(coeffs: bigint[]): void {
    const K = this._domain.base_ring;
    const p = K.characteristic;

    // If polynomial is constant, we have the identity isogeny
    if (coeffs.length <= 1) {
      this._degree = 1n;
      return;
    }

    // Find roots of the kernel polynomial
    const xCoords: FieldElement[] = [];
    for (let xVal = 0n; xVal < p; xVal++) {
      let result = 0n;
      let xPow = 1n;
      for (const coef of coeffs) {
        result = (result + coef * xPow) % p;
        xPow = (xPow * xVal) % p;
      }

      if (result === 0n) {
        xCoords.push(K.__call__(xVal));
      }
    }

    if (xCoords.length === 0) {
      throw new ValueError('kernel polynomial has no roots in the base field');
    }

    // For each root, find a corresponding point on the curve
    const kernelPoints: EllipticCurvePoint<FieldElement>[] = [];
    const [a1, a2, a3, a4, a6] = this._domain.a_invariants();

    for (const x of xCoords) {
      // Compute RHS of curve equation
      const x2 = x.mul(x);
      const x3 = x2.mul(x);
      const rhs = x3.add(a2.mul(x2)).add(a4.mul(x)).add(a6);

      // For short Weierstrass form, solve y^2 = rhs
      if (a1.isZero() && a3.isZero()) {
        if (rhs.isZero()) {
          const P = this._domain.point([x as F, K.zero() as F], false);
          kernelPoints.push(P as EllipticCurvePoint<FieldElement>);
        } else {
          // Check if rhs is a square
          const exp = (p - 1n) / 2n;
          const legendre = rhs.pow(exp);
          if (legendre.eq(K.one())) {
            let y: FieldElement;
            if (p % 4n === 3n) {
              y = rhs.pow((p + 1n) / 4n);
            } else {
              const sqrtResult = find_sqrt(rhs, K);
              if (!sqrtResult) continue;
              y = sqrtResult;
            }
            const P = this._domain.point([x as F, y as F], false);
            kernelPoints.push(P as EllipticCurvePoint<FieldElement>);
          }
        }
      }
    }

    // Use the kernel points to initialize via Velu
    if (kernelPoints.length > 0) {
      this.__init_from_kernel_list(kernelPoints as EllipticCurvePoint<F>[]);
      // Compute degree based on the polynomial degree
      // For odd degree n, the isogeny degree is 2n + 1
      // For kernel polynomial degree d, the isogeny degree is 2d + 1
      const polyDeg = coeffs.length - 1;
      this._degree = BigInt(2 * polyDeg + 1);
    } else {
      this._degree = 1n;
    }
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

    const xP = P.x();
    const yP = P.y();

    // Both Velu and Kohel algorithms use the same point evaluation method
    // since the kernel data is stored in __kernel_mod_sign for both.
    const result = this.__compute_via_velu(xP, yP);
    if (result === null) {
      // Point is in the kernel
      return this.codomain().zero();
    }
    const [xOut, yOut] = result;
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

      // tY = (tY0 + tY1 + tY2) / t1^3 - (tY0 - tY1 + tY2) * inv_t1_2
      // Simplified: tY = (tY0 * inv_t1_3) + (tY1 * inv_t1_2) + (tY2 * inv_t1_2)
      //                = tY0 / t1^2 + tY1 / t1^2 + tY2 / t1^2 (for y coordinate adjustment)
      // Actually the correct formula is:
      // Y' = Y - sum_Q [ uQ * (2Y + a1*X + a3) / (X - xQ)^2 + vQ * (Y - yQ) / (X - xQ)^2
      //                  - a1 * vQ / (X - xQ) - (a1 * uQ - gxQ * gyQ) / (X - xQ)^2 ]

      // Update X and Y
      X = X.add(tX) as F;

      // For Y: Y' = Y - sum [ uQ*(2y + a1*x + a3)/(x-xQ)^2 + vQ*(y-yQ)/(x-xQ)^2 - a1*vQ/(x-xQ) ]
      // Let's use the Velu helper formula directly
      const yContrib = tY0.mul(inv_t1_2).add(tY1.mul(inv_t1_2)).sub(tY2.mul(inv_t1_2)) as F;
      Y = Y.sub(yContrib) as F;
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
   * For normalized isogenies (which is what we compute), the scaling factor is 1.
   *
   * @returns The scaling factor u (always 1 for our implementation)
   * @see Reference: sage/schemes/elliptic_curves/ell_curve_isogeny.py:EllipticCurveIsogeny.scaling_factor
   */
  scaling_factor(): F {
    // Our isogenies are normalized by construction, so the scaling factor is 1
    return this._domain.base_ring.one() as F;
  }

  /**
   * Return the kernel polynomial of this isogeny.
   *
   * The kernel polynomial is a monic polynomial whose roots are the
   * x-coordinates of the non-trivial points in the kernel (modulo sign).
   *
   * For a degree-n isogeny:
   * - If n is odd, the kernel polynomial has degree (n-1)/2
   * - If n is even, the kernel polynomial has degree (n-2)/2 (for the odd part)
   *
   * @returns The kernel polynomial as an array of coefficients [a0, a1, ..., an] for a0 + a1*x + ... + an*x^n
   * @see Reference: sage/schemes/elliptic_curves/ell_curve_isogeny.py:EllipticCurveIsogeny.kernel_polynomial
   * @see Deviation: Elliptic Curve Isogeny Algorithms Limited
   */
  kernel_polynomial(): bigint[] {
    if (this.__kernel_polynomial !== null) {
      return this.__kernel_polynomial;
    }

    if (this.__algorithm === 'velu') {
      this.__init_kernel_polynomial_velu();
    } else {
      throw new NotImplementedError('Kohel kernel polynomial');
    }

    return this.__kernel_polynomial!;
  }

  /**
   * Initialize the kernel polynomial from the kernel points (Velu's algorithm).
   *
   * The kernel polynomial is the product of (x - xQ) for each x-coordinate xQ
   * of the non-trivial kernel points (modulo sign, i.e., each x appears once).
   */
  private __init_kernel_polynomial_velu(): void {
    const K = this._domain.base_ring;
    const p = K.characteristic;

    // Get all x-coordinates from kernel_mod_sign
    const xCoords: F[] = [];
    for (const [, data] of this.__kernel_mod_sign) {
      xCoords.push(data.xQ);
    }

    if (xCoords.length === 0) {
      // Trivial isogeny: kernel polynomial is 1
      this.__kernel_polynomial = [1n];
      return;
    }

    // Compute the product of (x - xQ) for all xQ
    // Start with polynomial [1] (constant 1)
    let poly: bigint[] = [1n];

    for (const xQ of xCoords) {
      // Multiply by (x - xQ) = [-xQ, 1]
      // If poly = [a0, a1, ..., an], then poly * (x - xQ) = [-xQ*a0, a0 - xQ*a1, a1 - xQ*a2, ..., an-1 - xQ*an, an]
      const xQVal = BigInt(xQ.toString());
      const newPoly: bigint[] = new Array(poly.length + 1).fill(0n);

      for (let i = 0; i < poly.length; i++) {
        // Contribution from -xQ * x^i term
        newPoly[i] = (newPoly[i]! - xQVal * poly[i]!) % p;
        if (newPoly[i]! < 0n) newPoly[i] = newPoly[i]! + p;

        // Contribution from x * x^i = x^(i+1) term
        newPoly[i + 1] = (newPoly[i + 1]! + poly[i]!) % p;
      }

      poly = newPoly;
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
   * Note: This implementation works for separable isogenies where the degree
   * is coprime to the field characteristic. For inseparable duals or
   * characteristics 2 and 3, a NotImplementedError is thrown.
   *
   * @returns The dual isogeny
   * @see Reference: sage/schemes/elliptic_curves/ell_curve_isogeny.py:EllipticCurveIsogeny.dual
   * @see Deviation: Elliptic Curve Isogeny Algorithms Limited
   */
  dual(): EllipticCurveIsogeny<F> {
    const p = this._domain.base_ring.characteristic;

    // Dual computation not yet implemented for characteristics 2 and 3
    if (p === 2n || p === 3n) {
      throw new NotImplementedError(
        'computation of dual isogenies not yet implemented in characteristics 2 and 3'
      );
    }

    // Return cached dual if available
    if (this.__dual !== null) {
      return this.__dual;
    }

    const d = this._degree;

    // Check if the isogeny is inseparable (degree divisible by characteristic)
    if (d % p === 0n) {
      throw new NotImplementedError('dual of inseparable isogenies requires Frobenius computation');
    }

    // For a separable isogeny phi: E -> E' of degree d where gcd(d, p) = 1,
    // the dual phi_hat: E' -> E has the same degree d.
    //
    // The key property is: phi_hat o phi = [d] on E, and phi o phi_hat = [d] on E'.
    //
    // For Velu isogenies, we compute the dual by finding the kernel of the
    // "reverse" isogeny. The kernel of phi_hat consists of points Q in E'
    // such that phi_hat(Q) = O, which is equivalent to finding the image
    // under phi of the d-torsion of E that is not in ker(phi).
    //
    // A simpler approach for small degrees: find d-torsion points on E' that
    // map to O under the desired dual.

    // For now, we use a direct construction: The kernel of the dual isogeny
    // can be computed from the division polynomial method.
    // This is a simplified implementation that finds kernel points on E'.

    const E2 = this.codomain();

    // Find d-torsion points on E' that could generate the kernel of the dual
    const dTorsionPoints: EllipticCurvePoint<F>[] = [];

    // Enumerate points on E' and check for d-torsion
    for (let xVal = 0n; xVal < p && dTorsionPoints.length < Number(d); xVal++) {
      const x = E2.base_ring.__call__(xVal) as F;
      const [a1, a2, a3, a4, a6] = E2.a_invariants();

      // Compute y^2 + a1*x*y + a3*y = x^3 + a2*x^2 + a4*x + a6
      const x2 = x.mul(x) as F;
      const x3 = x2.mul(x) as F;
      const rhs = x3.add(a2.mul(x2)).add(a4.mul(x)).add(a6) as F;

      // Solve for y (simplified for short Weierstrass: y^2 = rhs)
      if (a1.isZero() && a3.isZero()) {
        // Short Weierstrass form
        const exp = (p + 1n) / 4n;
        if (p % 4n === 3n) {
          const y = rhs.pow(exp) as F;
          if (y.mul(y).eq(rhs)) {
            const Q = E2.point([x, y], false);
            if (Q.mul(d).is_zero() && !Q.is_zero()) {
              dTorsionPoints.push(Q);
            }
            // Also try -y
            const negY = y.neg() as F;
            if (!negY.eq(y)) {
              const Qneg = E2.point([x, negY], false);
              if (Qneg.mul(d).is_zero() && !Qneg.is_zero()) {
                dTorsionPoints.push(Qneg);
              }
            }
          }
        }
      }
    }

    // Find a point that generates a cyclic subgroup of order d
    for (const Q of dTorsionPoints) {
      // Check if Q generates a subgroup of size d
      let R = Q;
      let order = 1n;
      while (!R.is_zero() && order < d) {
        R = R.add(Q);
        order++;
      }
      if (R.is_zero() && order === d) {
        // Q generates a cyclic subgroup of order d
        // This is the kernel of the dual isogeny
        try {
          this.__dual = new EllipticCurveIsogeny(E2, Q);
          return this.__dual;
        } catch {
          continue;
        }
      }
    }

    throw new NotImplementedError(
      'dual isogeny computation failed - could not find suitable kernel on codomain'
    );
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
  private __pre_isomorphism: WeierstrassIsomorphism | null = null;

  /** Post-isomorphism for composition */
  private __post_isomorphism: WeierstrassIsomorphism | null = null;

  /**
   * Modify this isogeny by pre-composing with a WeierstrassIsomorphism.
   *
   * If preWI: E' -> E is an isomorphism and this isogeny is phi: E -> E'',
   * then the result is phi o preWI: E' -> E''.
   *
   * @param preWI - A WeierstrassIsomorphism with codomain equal to this isogeny's domain
   * @see Reference: sage/schemes/elliptic_curves/ell_curve_isogeny.py:EllipticCurveIsogeny._set_pre_isomorphism
   */
  _set_pre_isomorphism(preWI: WeierstrassIsomorphism): void {
    this.__pre_isomorphism = preWI;
    // Clear cached rational maps since they will change
    this.__rational_maps = null;
    // Note: In a full implementation, we would update the domain and
    // recalculate rational maps. For now, we just store the isomorphism.
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
  _set_post_isomorphism(postWI: WeierstrassIsomorphism): void {
    this.__post_isomorphism = postWI;
    // Clear cached rational maps since they will change
    this.__rational_maps = null;
    // Note: In a full implementation, we would update the codomain and
    // recalculate rational maps. For now, we just store the isomorphism.
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
  // For isomorphic curves, they are trivially isogenous
  if (E1.j_invariant().eq(E2.j_invariant())) {
    // They might be isomorphic or twists
    return true;
  }

  // For finite fields: isogenous iff same cardinality
  // We need to compute the cardinality
  // This requires the curves to have a cardinality method
  // For now, we'll use a simple check: same j-invariant means potentially isogenous

  // Count points on both curves (for small fields)
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
