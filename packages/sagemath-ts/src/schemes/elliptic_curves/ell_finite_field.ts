/**
 * @module sage/schemes/elliptic_curves/ell_finite_field
 * @description Elliptic curves over finite fields
 *
 * Port of: sage/schemes/elliptic_curves/ell_finite_field.py
 *
 * This module implements elliptic curve operations over finite fields.
 * Following SageMath's architecture, we delegate the core computations
 * to PARI/GP (via parigp-ts):
 * - Point order calculation -> ellorder
 * - Curve cardinality (order) -> ellcard
 * - Random point generation -> FpE_random
 * - Generator finding -> ellgenerators
 */

import {
  type EllipticCurveFp,
  // Point types and operations
  type EllipticPointFp,
  FpE_add,
  FpE_dbl,
  FpE_mul,
  FpE_neg,
  FpE_random,
  ell_is_inf,
  // Curve operations
  ellcard,
  elldivpol,
  ellgenerators,
  ellinf,
  // Curve initialization
  ellinit_Fp,
  ellisoncurve,
  elllift_x,
  ellorder,
  ellpoint,
  elltatepairing,
  // Pairings
  ellweilpairing,
  ellxn,
  trace_of_frobenius,
} from '@sagemath-ts/parigp-ts';
import { divisors, factor, gcd, is_prime, isqrt, lcm } from '../../arith/misc.js';
import {
  ArithmeticError,
  NotImplementedError,
  ValueError,
  ZeroDivisionError,
} from '../../errors.js';
import { has_order as generic_has_order } from '../../groups/generic.js';
import type {
  FiniteFieldElement,
  FiniteFieldPrime,
} from '../../rings/finite_rings/finite_field_prime.js';
import type { CoefficientRing, RingElement } from '../../rings/polynomial/polynomial_element.js';
import { PolynomialRing } from '../../rings/polynomial/polynomial_ring.js';
import { type IntegerLike, toBigInt } from '../../types/coercion.js';

/**
 * Type alias for field element (supporting both prime and extension fields)
 */
export type FieldElement = FiniteFieldElement;
export type BaseField = FiniteFieldPrime;

/**
 * Represents a point on an elliptic curve.
 *
 * Points can be either:
 * - The point at infinity (identity element)
 * - An affine point (x, y)
 *
 * Point operations delegate to PARI/GP functions via parigp-ts,
 * matching SageMath's architecture.
 */
export class EllipticCurvePoint {
  readonly curve: EllipticCurveFiniteField;
  readonly x: FieldElement | null;
  readonly y: FieldElement | null;
  readonly isInfinity: boolean;
  private _order: bigint | null = null;

  /**
   * Create a point on an elliptic curve.
   *
   * @param curve - The parent elliptic curve
   * @param x - x-coordinate (null for point at infinity)
   * @param y - y-coordinate (null for point at infinity)
   */
  constructor(curve: EllipticCurveFiniteField, x: FieldElement | null, y: FieldElement | null) {
    this.curve = curve;

    if (x === null && y === null) {
      this.isInfinity = true;
      this.x = null;
      this.y = null;
    } else if (x !== null && y !== null) {
      this.isInfinity = false;
      this.x = x;
      this.y = y;
    } else {
      throw new ValueError('Both x and y must be provided, or both must be null for infinity');
    }
  }

  /**
   * Create the point at infinity.
   */
  static infinity(curve: EllipticCurveFiniteField): EllipticCurvePoint {
    return new EllipticCurvePoint(curve, null, null);
  }

  /**
   * Check if this is the point at infinity.
   */
  isZero(): boolean {
    return this.isInfinity;
  }

  /**
   * Convert to PARI/GP point format.
   *
   * This is the equivalent of SageMath's __pari__() method.
   * Returns an EllipticPointFp that can be used with parigp-ts functions.
   */
  toPari(): EllipticPointFp {
    if (this.isInfinity) {
      return ellinf();
    }
    return ellpoint(this.x!.value, this.y!.value);
  }

  /**
   * Create an EllipticCurvePoint from a PARI point.
   */
  static fromPari(curve: EllipticCurveFiniteField, pariPoint: EllipticPointFp): EllipticCurvePoint {
    if (ell_is_inf(pariPoint)) {
      return EllipticCurvePoint.infinity(curve);
    }
    const x = curve.field.__call__(pariPoint.x!);
    const y = curve.field.__call__(pariPoint.y!);
    return new EllipticCurvePoint(curve, x, y);
  }

  /**
   * Check equality with another point.
   */
  eq(other: EllipticCurvePoint): boolean {
    if (this.isInfinity && other.isInfinity) {
      return true;
    }
    if (this.isInfinity || other.isInfinity) {
      return false;
    }
    return this.x!.eq(other.x!) && this.y!.eq(other.y!);
  }

  /**
   * Return the negation of this point.
   *
   * Delegates to PARI's FpE_neg.
   * For short Weierstrass form y^2 = x^3 + ax + b: -P = (x, -y)
   */
  neg(): EllipticCurvePoint {
    if (this.isInfinity) {
      return this;
    }
    const p = this.curve.field.characteristic;
    const result = FpE_neg(this.toPari(), p);
    return EllipticCurvePoint.fromPari(this.curve, result);
  }

  /**
   * Add this point to another point.
   *
   * Delegates to PARI's FpE_add.
   */
  add(other: EllipticCurvePoint): EllipticCurvePoint {
    const pariCurve = this.curve.toPari();
    const result = FpE_add(this.toPari(), other.toPari(), pariCurve.a4, pariCurve.p);
    return EllipticCurvePoint.fromPari(this.curve, result);
  }

  /**
   * Double this point.
   *
   * Delegates to PARI's FpE_dbl.
   */
  double(): EllipticCurvePoint {
    if (this.isInfinity) {
      return this;
    }
    const pariCurve = this.curve.toPari();
    const result = FpE_dbl(this.toPari(), pariCurve.a4, pariCurve.p);
    return EllipticCurvePoint.fromPari(this.curve, result);
  }

  /**
   * Subtract another point from this point.
   */
  sub(other: EllipticCurvePoint): EllipticCurvePoint {
    return this.add(other.neg());
  }

  /**
   * Scalar multiplication: compute [n]P
   *
   * Delegates to PARI's FpE_mul (double-and-add algorithm).
   */
  mul(n: bigint): EllipticCurvePoint {
    if (n === 0n || this.isInfinity) {
      return EllipticCurvePoint.infinity(this.curve);
    }
    const pariCurve = this.curve.toPari();
    const result = FpE_mul(this.toPari(), n, pariCurve.a4, pariCurve.p);
    return EllipticCurvePoint.fromPari(this.curve, result);
  }

  /**
   * Compute the order of this point.
   *
   * Delegates to PARI's ellorder.
   * SageMath: self.__pari__().ellorder()
   *
   * The order of a point P is the smallest positive integer n such that [n]P = O.
   * The order always divides the curve order.
   *
   * @see Deviation: PARI Elliptic Curve Advanced Algorithms Missing (parigp-ts)
   */
  order(): bigint {
    if (this._order !== null) {
      return this._order;
    }

    if (this.isInfinity) {
      this._order = 1n;
      return 1n;
    }

    // Delegate to PARI's ellorder
    const pariCurve = this.curve.toPari();
    const curveOrder = ellcard(pariCurve);
    this._order = ellorder(pariCurve, this.toPari(), curveOrder);
    return this._order;
  }

  /**
   * Set the cached order of this point.
   *
   * Use this when the order is known a priori to avoid computation.
   *
   * @param n - The order of the point
   * @param check - Whether to verify the order (default: true)
   */
  setOrder(n: bigint, check: boolean = true): void {
    if (check) {
      if (!this.mul(n).isZero()) {
        throw new ValueError(`${n} is not the order of this point`);
      }
      // Check that no smaller divisor works
      const divs = divisors(n);
      for (const d of divs) {
        if (d < n && this.mul(d).isZero()) {
          throw new ValueError(`${n} is not the order of this point (order is at most ${d})`);
        }
      }
    }
    this._order = n;
  }

  /**
   * Check if this point has exactly the given order n in the group.
   *
   * This verifies both:
   * 1. n*P = O (identity)
   * 2. For all prime divisors p of n: (n/p)*P != O
   *
   * The second condition ensures n is the exact order, not just a multiple.
   *
   * @param n - The proposed order
   * @returns true if the order of this point is exactly n
   */
  has_order(n: bigint | number): boolean {
    return generic_has_order(this, n, '+');
  }

  /**
   * Compute the Weil pairing of this point with another point Q.
   *
   * The Weil pairing is a bilinear, alternating, non-degenerate pairing
   * e_n: E[n] x E[n] -> mu_n where mu_n is the group of n-th roots of unity.
   *
   * Properties:
   * - Bilinearity: e_n(aP, bQ) = e_n(P, Q)^(ab)
   * - Alternating: e_n(P, P) = 1
   * - Non-degeneracy: if e_n(P, Q) = 1 for all Q in E[n], then P = O
   *
   * @param Q - Another n-torsion point on the same curve
   * @param n - Integer such that nP = nQ = O
   * @returns An n-th root of unity in the base field
   * @throws ValueError if points are not n-torsion
   *
   * @example
   * ```typescript
   * const E = EllipticCurve(GF(103n), [1n, 18n]);
   * const P = E.point(33n, 91n);
   * const Q = E.point(87n, 51n);
   * const n = 19n; // torsion order
   * const z = P.weil_pairing(Q, n); // n-th root of unity
   * ```
   */
  weil_pairing(Q: EllipticCurvePoint, n: bigint): FiniteFieldElement {
    // Delegate to PARI's ellweilpairing
    const pariCurve = this.curve.toPari();

    // Test if P, Q are both in E[n]
    if (!this.mul(n).isZero() || !Q.mul(n).isZero()) {
      throw new ValueError('points must both be n-torsion');
    }

    // Handle trivial cases
    if (this.isZero() || Q.isZero() || this.eq(Q)) {
      return this.curve.field.one();
    }

    // Call PARI's ellweilpairing
    const result = ellweilpairing(pariCurve, this.toPari(), Q.toPari(), n);
    return this.curve.field.__call__(result);
  }

  /**
   * Compute the Tate pairing of this n-torsion point with Q.
   *
   * The value returned is f_{n,P}(Q)^e where f_{n,P} is a function with
   * divisor n[P]-n[O]. This is the "modified Tate pairing".
   *
   * @param Q - Elliptic curve point on same curve as P
   * @param n - Positive integer, order of P
   * @param k - Positive integer, embedding degree (optional, computed if not given)
   * @returns An n-th root of unity in F_q^k
   *
   * @example
   * ```typescript
   * const p = 103n;
   * const E = EllipticCurve(GF(p), [1n, 18n]);
   * const P = E.point(33n, 91n);
   * const n = P.order(); // 19
   * const k = E.embedding_degree(n); // 6
   * const result = P.tate_pairing(Q, n, k);
   * ```
   */
  tate_pairing(Q: EllipticCurvePoint, n: bigint, k?: bigint | number): FiniteFieldElement {
    const pariCurve = this.curve.toPari();
    const p = this.curve.field.characteristic;

    // Compute embedding degree if not provided
    let kVal: bigint;
    if (k === undefined) {
      kVal = this.curve.embedding_degree(n);
    } else {
      kVal = typeof k === 'number' ? BigInt(k) : k;
    }

    // Check that P is n-torsion
    if (!this.mul(n).isZero()) {
      throw new ValueError('The point P must be n-torsion');
    }

    // Handle trivial cases
    if (this.isZero() || Q.isZero()) {
      return this.curve.field.one();
    }

    // Call PARI's elltatepairing (returns raw Miller function value)
    const millerValue = elltatepairing(pariCurve, this.toPari(), Q.toPari(), n);

    // Compute final exponentiation: (p^k - 1) / n
    let pk = 1n;
    for (let i = 0n; i < kVal; i++) {
      pk *= p;
    }
    const exp = (pk - 1n) / n;

    // Return millerValue^exp
    const result = modPow(millerValue, exp, p);
    return this.curve.field.__call__(result);
  }

  /**
   * Compute the Ate pairing of points P and Q.
   *
   * The Ate pairing is an optimized variant of the Tate pairing.
   * P should be in ker(pi - 1), Q should be in ker(pi - q).
   *
   * @param Q - Point of order n in ker(pi - q)
   * @param n - The order of P and Q
   * @param k - The embedding degree
   * @param t - The trace of Frobenius
   * @returns An n-th root of unity in F_q^k
   *
   * @example
   * ```typescript
   * const p = 7549n;
   * const E = EllipticCurve(GF(p), [0n, 1n]);
   * const P = E.point(3050n, 5371n);
   * const n = 157n;
   * const k = 6;
   * const t = E.trace_of_frobenius();
   * const result = P.ate_pairing(Q, n, k, t);
   * ```
   */
  ate_pairing(Q: EllipticCurvePoint, n: bigint, k: bigint | number, t: bigint): FiniteFieldElement {
    const pariCurve = this.curve.toPari();
    const p = this.curve.field.characteristic;
    const kVal = typeof k === 'number' ? BigInt(k) : k;

    // Handle trivial cases
    if (this.isZero() || Q.isZero()) {
      return this.curve.field.one();
    }

    // T = t - 1
    const T = t - 1n;
    const absT = T >= 0n ? T : -T;

    // Compute the Miller function f_{|T|,Q}(P)
    const millerValue = elltatepairing(pariCurve, Q.toPari(), this.toPari(), absT);

    // If T is negative, take reciprocal
    let adjustedMiller = millerValue;
    if (T < 0n) {
      adjustedMiller = modPow(millerValue, p - 2n, p); // modular inverse
    }

    // Compute final exponentiation: (p^k - 1) / n
    let pk = 1n;
    for (let i = 0n; i < kVal; i++) {
      pk *= p;
    }
    const exp = (pk - 1n) / n;

    // Return adjustedMiller^exp
    const result = modPow(adjustedMiller, exp, p);
    return this.curve.field.__call__(result);
  }

  /**
   * String representation.
   */
  toString(): string {
    if (this.isInfinity) {
      return '(0 : 1 : 0)';
    }
    return `(${this.x!.toString()} : ${this.y!.toString()} : 1)`;
  }

  /**
   * Repr for debugging.
   */
  repr(): string {
    return this.toString();
  }
}

/**
 * An elliptic curve over a finite field in short Weierstrass form.
 *
 * The curve is defined by: y^2 = x^3 + ax + b
 *
 * Following SageMath's architecture, the core computations are delegated
 * to PARI/GP (via parigp-ts):
 * - cardinality() -> ellcard
 * - gens() -> ellgenerators
 * - Point operations -> FpE_add, FpE_mul, etc.
 *
 * @example
 * ```typescript
 * const F = GF(101n);
 * const E = new EllipticCurveFiniteField(F, 2n, 3n);
 * console.log(E.cardinality());  // Delegates to PARI's ellcard
 *
 * const P = E.random_point();    // Delegates to PARI's FpE_random
 * console.log(P.order());        // Delegates to PARI's ellorder
 * ```
 */
export class EllipticCurveFiniteField {
  readonly field: BaseField;
  readonly a: FieldElement;
  readonly b: FieldElement;
  private _pariCurve: EllipticCurveFp | null = null;
  private _order: bigint | null = null;
  private _generators: EllipticCurvePoint[] | null = null;

  /**
   * Create an elliptic curve y^2 = x^3 + ax + b over a finite field.
   *
   * @param field - The base finite field
   * @param a - Coefficient a (or coefficient of x)
   * @param b - Coefficient b (constant term)
   * @param checkNonsingular - Whether to verify the curve is non-singular (default: true)
   */
  constructor(
    field: BaseField,
    a: bigint | number | FieldElement,
    b: bigint | number | FieldElement,
    checkNonsingular: boolean = true
  ) {
    this.field = field;

    // Convert coefficients to field elements
    if (typeof a === 'bigint' || typeof a === 'number') {
      this.a = field.__call__(a);
    } else {
      this.a = a;
    }

    if (typeof b === 'bigint' || typeof b === 'number') {
      this.b = field.__call__(b);
    } else {
      this.b = b;
    }

    // Check non-singularity: 4a^3 + 27b^2 != 0
    if (checkNonsingular) {
      const four = field.__call__(4n);
      const twentySeven = field.__call__(27n);
      const discriminant = four.mul(this.a.pow(3)).add(twentySeven.mul(this.b.pow(2)));

      if (discriminant.isZero()) {
        throw new ArithmeticError(
          `Curve y^2 = x^3 + ${this.a}*x + ${this.b} is singular (discriminant = 0)`
        );
      }
    }
  }

  /**
   * Convert to PARI/GP curve format.
   *
   * This is the equivalent of SageMath's __pari__() method.
   * Returns an EllipticCurveFp that can be used with parigp-ts functions.
   *
   * SageMath: self.__pari__() returns a PARI ell object
   */
  toPari(): EllipticCurveFp {
    if (this._pariCurve === null) {
      this._pariCurve = ellinit_Fp(this.a.value, this.b.value, this.field.characteristic);
    }
    return this._pariCurve;
  }

  /**
   * Alias for toPari() to match SageMath's naming convention.
   */
  __pari__(): EllipticCurveFp {
    return this.toPari();
  }

  /**
   * Return the characteristic of the base field.
   */
  get characteristic(): bigint {
    return this.field.characteristic;
  }

  /**
   * Return the discriminant of the curve.
   *
   * For y^2 = x^3 + ax + b: disc = -16(4a^3 + 27b^2)
   */
  discriminant(): FieldElement {
    const four = this.field.__call__(4n);
    const twentySeven = this.field.__call__(27n);
    const minusSixteen = this.field.__call__(-16n);
    return minusSixteen.mul(four.mul(this.a.pow(3)).add(twentySeven.mul(this.b.pow(2))));
  }

  /**
   * Return the j-invariant of the curve.
   *
   * j = 1728 * 4a^3 / (4a^3 + 27b^2)
   */
  j_invariant(): FieldElement {
    const four = this.field.__call__(4n);
    const twentySeven = this.field.__call__(27n);
    const val1728 = this.field.__call__(1728n);

    const fourACubed = four.mul(this.a.pow(3));
    const denominator = fourACubed.add(twentySeven.mul(this.b.pow(2)));

    if (denominator.isZero()) {
      throw new ZeroDivisionError('j-invariant undefined for singular curve');
    }

    return val1728.mul(fourACubed).div(denominator);
  }

  /**
   * Return the point at infinity (identity element).
   */
  zero(): EllipticCurvePoint {
    return EllipticCurvePoint.infinity(this);
  }

  /**
   * Create a point on this curve.
   *
   * @param x - x-coordinate
   * @param y - y-coordinate
   * @param check - Whether to verify the point is on the curve (default: true)
   */
  point(
    x: bigint | number | FieldElement,
    y: bigint | number | FieldElement,
    check?: boolean
  ): EllipticCurvePoint;
  point(
    coords: [bigint | number | FieldElement, bigint | number | FieldElement],
    check?: boolean
  ): EllipticCurvePoint;
  point(
    xOrCoords:
      | bigint
      | number
      | FieldElement
      | [bigint | number | FieldElement, bigint | number | FieldElement],
    yOrCheck?: bigint | number | FieldElement | boolean,
    check: boolean = true
  ): EllipticCurvePoint {
    let x: bigint | number | FieldElement;
    let y: bigint | number | FieldElement;
    let doCheck = check;

    if (Array.isArray(xOrCoords)) {
      [x, y] = xOrCoords;
      if (typeof yOrCheck === 'boolean') {
        doCheck = yOrCheck;
      }
    } else {
      x = xOrCoords;
      if (typeof yOrCheck === 'boolean' || yOrCheck === undefined) {
        throw new ValueError('y-coordinate must be provided');
      }
      y = yOrCheck;
    }

    // Convert to field elements
    const xElem = typeof x === 'bigint' || typeof x === 'number' ? this.field.__call__(x) : x;
    const yElem = typeof y === 'bigint' || typeof y === 'number' ? this.field.__call__(y) : y;

    if (doCheck && !this.is_on_curve(xElem, yElem)) {
      throw new ValueError(`Point (${xElem}, ${yElem}) is not on the curve`);
    }

    return new EllipticCurvePoint(this, xElem, yElem);
  }

  /**
   * Check if a point (x, y) is on the curve.
   *
   * Delegates to PARI's ellisoncurve.
   */
  is_on_curve(x: FieldElement, y: FieldElement): boolean {
    const pariCurve = this.toPari();
    const pariPoint = ellpoint(x.value, y.value);
    return ellisoncurve(pariCurve, pariPoint);
  }

  /**
   * Check if there is a point with the given x-coordinate.
   *
   * This checks if x^3 + ax + b is a quadratic residue.
   */
  is_x_coord(x: FieldElement | bigint | number): boolean {
    const xElem = typeof x === 'bigint' || typeof x === 'number' ? this.field.__call__(x) : x;

    // Compute y^2 = x^3 + ax + b
    const ySquared = xElem.pow(3).add(this.a.mul(xElem)).add(this.b);

    return ySquared.is_square();
  }

  /**
   * Find point(s) with a given x-coordinate.
   *
   * Delegates to PARI's elllift_x.
   *
   * @param x - The x-coordinate
   * @param all - If true, return all points (0, 1, or 2); if false, return one or throw
   * @returns Point or array of points
   *
   * @example
   * ```typescript
   * const E = new EllipticCurveFiniteField(GF(101n), 2n, 3n);
   * const P = E.lift_x(5n);           // Get one point with x=5
   * const pts = E.lift_x(5n, true);   // Get all points with x=5
   * ```
   */
  lift_x(x: FieldElement | bigint | number, all?: false): EllipticCurvePoint;
  lift_x(x: FieldElement | bigint | number, all: true): EllipticCurvePoint[];
  lift_x(
    x: FieldElement | bigint | number,
    all: boolean = false
  ): EllipticCurvePoint | EllipticCurvePoint[] {
    const xElem = typeof x === 'bigint' || typeof x === 'number' ? this.field.__call__(x) : x;
    const xVal = xElem.value;

    // Use PARI's elllift_x
    const pariCurve = this.toPari();
    const pariPoint = elllift_x(pariCurve, xVal);

    if (pariPoint === null) {
      if (all) {
        return [];
      }
      throw new ValueError(`No point with x-coordinate ${xElem} on ${this}`);
    }

    // We got one point, now construct both if needed
    const y = pariPoint.y!;
    const p = this.field.characteristic;
    const negY = (p - y) % p;

    const point1 = new EllipticCurvePoint(this, this.field.__call__(xVal), this.field.__call__(y));

    if (y === 0n || y === negY) {
      // Only one point (y = 0 or char = 2)
      return all ? [point1] : point1;
    }

    const point2 = new EllipticCurvePoint(
      this,
      this.field.__call__(xVal),
      this.field.__call__(negY)
    );

    // Sort for deterministic behavior (smaller y first)
    const points = y < negY ? [point1, point2] : [point2, point1];

    if (all) {
      return points;
    }
    return points[0]!;
  }

  /**
   * Return a random point on the curve.
   *
   * Delegates to PARI's FpE_random.
   * SageMath: Uses PARI's random_FpE internally
   *
   * @returns A random non-identity point on the curve
   */
  random_point(): EllipticCurvePoint {
    const pariCurve = this.toPari();
    const pariPoint = FpE_random(pariCurve);
    return EllipticCurvePoint.fromPari(this, pariPoint);
  }

  /**
   * Compute the cardinality (order) of the curve.
   *
   * Delegates to PARI's ellcard.
   * SageMath: return self.__pari__().ellcard()
   *
   * This is the number of points on E(F_q), including the point at infinity.
   *
   * @see Deviation: PARI Elliptic Curve Advanced Algorithms Missing (parigp-ts)
   */
  cardinality(): bigint {
    if (this._order !== null) {
      return this._order;
    }

    // Delegate to PARI's ellcard
    const pariCurve = this.toPari();
    this._order = ellcard(pariCurve);
    return this._order;
  }

  /**
   * Compute cardinality using PARI (explicit method).
   *
   * This is the equivalent of SageMath's cardinality_pari() method.
   * SageMath: return Integer(self.__pari__().ellcard())
   *
   * @see Deviation: PARI Elliptic Curve Advanced Algorithms Missing (parigp-ts)
   */
  cardinality_pari(): bigint {
    return ellcard(this.toPari());
  }

  /**
   * Set the cached curve order.
   *
   * Use this when the order is known a priori (e.g., for standard curves).
   *
   * @param n - The curve order
   * @param check - Whether to verify the order (default: true)
   */
  set_order(n: bigint, check: boolean = true, num_checks: number = 8): void {
    if (check && !this.has_order(n, num_checks)) {
      throw new ValueError(`${this} does not have order ${n}`);
    }

    this._order = n;
  }

  /**
   * Return True if the curve has order ``value``.
   *
   * INPUT:
   * - value: integer in the Hasse-Weil range for this curve
   * - num_checks: the number of times to check whether ``value`` times a
   *   random point on this curve equals the identity (default 8)
   *
   * @see Reference: sage/schemes/elliptic_curves/ell_finite_field.py:has_order
   */
  has_order(value: bigint, num_checks: number = 8): boolean {
    const q = this.field.order;

    // Hasse_bounds(q, 1): rq = isqrt(4*q); (q+1-rq, q+1+rq)
    const rq = isqrt(4n * q);
    const lower = q + 1n - rq;
    const upper = q + 1n + rq;
    if (value < lower || value > upper) {
      return false;
    }

    // For really small values, the random tests are too weak to detect wrong
    // orders, so we compute directly instead (Sage uses the same q <= 2^64
    // cutoff, see ell_finite_field.py:1320-1323).
    if (q <= 1n << 64n || this._order !== null) {
      return this.order() === value;
    }

    for (let i = 0; i < num_checks; i++) {
      let G = this.random_point();
      while (G.isZero()) {
        G = this.random_point();
      }
      if (!G.mul(value).isZero()) {
        return false;
      }
    }

    return true;
  }

  /**
   * Return the trace of Frobenius.
   *
   * Delegates to PARI's trace_of_frobenius.
   *
   * For a curve over F_q, the trace t satisfies #E = q + 1 - t.
   *
   * @see Deviation: PARI Elliptic Curve Advanced Algorithms Missing (parigp-ts)
   */
  trace_of_frobenius(): bigint {
    return trace_of_frobenius(this.toPari());
  }

  /**
   * Return generators of the group E(F_q).
   *
   * Delegates to PARI's ellgenerators.
   * SageMath: card, ords, pts = self.__pari__().ellgroup(flag=1)
   *
   * For prime-order curves, any non-identity point is a generator.
   * For curves with composite order, this returns a set of generators.
   *
   * @returns Array of generator points
   */
  gens(): EllipticCurvePoint[] {
    if (this._generators !== null) {
      return this._generators;
    }

    // Sage: card, ords, pts = self.__pari__().ellgroup(flag=1)
    const pariCurve = this.toPari();
    const pariGens = ellgenerators(pariCurve);

    const gens = pariGens.map((pg) => EllipticCurvePoint.fromPari(this, pg));

    // PARI documentation: "P is of order d_1", and the returned points
    // generate the whole group. Assert this, so that a broken ellgroup
    // surfaces here rather than silently corrupting abelian_group().
    const n = this.cardinality();
    let generated = 1n;
    for (const P of gens) {
      generated = lcm(generated, P.order());
    }
    if (gens.length === 1 && generated !== n) {
      throw new ArithmeticError(
        `gens(): the point returned by ellgroup generates a group of order ${generated}, not ${n}`
      );
    }
    if (gens.length === 2) {
      const [P, Q] = gens as [EllipticCurvePoint, EllipticCurvePoint];
      const n1 = P.order();
      if (n % n1 !== 0n || !Q.mul(n1).isZero()) {
        throw new ArithmeticError(
          'gens(): ellgroup returned generators that do not satisfy n1 | #E and n1*Q = O'
        );
      }
    }

    this._generators = gens;
    return this._generators;
  }

  /**
   * Return all rational points on the curve (including infinity).
   *
   * Warning: This enumerates all points, which is only practical for small fields.
   */
  points(): EllipticCurvePoint[] {
    const pts: EllipticCurvePoint[] = [this.zero()];

    for (const x of this.field) {
      const ySquared = x.pow(3).add(this.a.mul(x)).add(this.b);

      if (ySquared.isZero()) {
        pts.push(new EllipticCurvePoint(this, x, this.field.zero()));
      } else if (ySquared.is_square()) {
        const y = ySquared.sqrt();
        pts.push(new EllipticCurvePoint(this, x, y));
        if (!y.isZero()) {
          pts.push(new EllipticCurvePoint(this, x, y.neg()));
        }
      }
    }

    return pts;
  }

  /**
   * Alias for cardinality().
   */
  order(): bigint {
    return this.cardinality();
  }

  /**
   * Return the embedding degree of this curve with respect to n.
   *
   * The embedding degree is the smallest positive integer k such that
   * n divides q^k - 1, where q is the size of the base field.
   *
   * This is important for pairing-based cryptography as it determines
   * the extension field where the pairing takes values.
   *
   * @param n - A positive integer (typically the order of a subgroup)
   * @returns The embedding degree k
   *
   * @example
   * ```typescript
   * const E = EllipticCurve(GF(103n), [1n, 18n]);
   * const P = E.random_point();
   * const n = P.order();
   * const k = E.embedding_degree(n);
   * ```
   */
  embedding_degree(n: bigint): bigint {
    return embedding_degree(this, n);
  }

  /**
   * Return True if this elliptic curve is supersingular.
   *
   * An elliptic curve over a field of characteristic p is supersingular
   * if and only if its trace of Frobenius t satisfies t ≡ 0 (mod p).
   *
   * For curves over F_p, E is supersingular iff #E(F_p) = p + 1 (i.e., t = 0).
   *
   * @param proof - If True (default), return a proved result
   * @returns True if E is supersingular, False otherwise
   */
  is_supersingular(proof?: boolean): boolean {
    return is_supersingular(this, proof);
  }

  /**
   * Return True if this elliptic curve is ordinary.
   *
   * An elliptic curve is ordinary if and only if it is not supersingular.
   *
   * @param proof - If True (default), return a proved result
   * @returns True if E is ordinary, False otherwise
   */
  is_ordinary(proof?: boolean): boolean {
    return is_ordinary(this, proof);
  }

  /**
   * String representation.
   */
  toString(): string {
    return `Elliptic Curve defined by y^2 = x^3 + ${this.a}*x + ${this.b} over ${this.field}`;
  }
}

/**
 * Modular exponentiation: compute base^exp mod mod.
 */
function modPow(base: bigint, exp: bigint, mod: bigint): bigint {
  if (exp < 0n) {
    throw new ArithmeticError('Negative exponent not supported');
  }
  if (mod === 1n) return 0n;

  let result = 1n;
  base = ((base % mod) + mod) % mod;

  while (exp > 0n) {
    if ((exp & 1n) === 1n) {
      result = (result * base) % mod;
    }
    exp >>= 1n;
    base = (base * base) % mod;
  }

  return result;
}

/**
 * Create an elliptic curve over a prime field.
 *
 * @param field - The prime finite field
 * @param coeffs - Either [a, b] for y^2 = x^3 + ax + b, or [a1, a2, a3, a4, a6] for general form
 *
 * @see Deviation: Elliptic Curve Short Weierstrass Form Only
 *
 * @example
 * ```typescript
 * const F = GF(101n);
 * const E = EllipticCurve(F, [2n, 3n]);  // y^2 = x^3 + 2x + 3
 * ```
 */
export function EllipticCurve(
  field: BaseField,
  coeffs: [bigint | number, bigint | number]
): EllipticCurveFiniteField;
export function EllipticCurve(
  field: BaseField,
  coeffs:
    | [bigint | number, bigint | number]
    | [bigint | number, bigint | number, bigint | number, bigint | number, bigint | number]
): EllipticCurveFiniteField {
  if (coeffs.length === 2) {
    // Short Weierstrass form: y^2 = x^3 + ax + b
    return new EllipticCurveFiniteField(field, coeffs[0], coeffs[1]);
  } else if (coeffs.length === 5) {
    // General Weierstrass: y^2 + a1*x*y + a3*y = x^3 + a2*x^2 + a4*x + a6
    // For now, we only support char != 2, 3 and convert to short form
    const [a1, a2, a3, a4, a6] = coeffs.map((c) => field.__call__(c));

    if (field.characteristic === 2n) {
      throw new ValueError('General Weierstrass form in characteristic 2 not yet supported');
    }
    if (field.characteristic === 3n) {
      throw new ValueError('General Weierstrass form in characteristic 3 not yet supported');
    }

    // Convert to short form using standard transformation
    // b2 = a1^2 + 4*a2
    // b4 = a1*a3 + 2*a4
    // b6 = a3^2 + 4*a6
    // c4 = b2^2 - 24*b4
    // c6 = -b2^3 + 36*b2*b4 - 216*b6
    // a = -c4/48, b = -c6/864
    const four = field.__call__(4n);
    const two = field.__call__(2n);
    const twentyFour = field.__call__(24n);
    const thirtySix = field.__call__(36n);
    const twoSixteen = field.__call__(216n);
    const fortyEight = field.__call__(48n);
    const eightSixtyFour = field.__call__(864n);

    const b2 = a1.mul(a1).add(four.mul(a2));
    const b4 = a1.mul(a3).add(two.mul(a4));
    const b6 = a3.mul(a3).add(four.mul(a6));

    const c4 = b2.mul(b2).sub(twentyFour.mul(b4));
    const c6 = b2.pow(3).neg().add(thirtySix.mul(b2).mul(b4)).sub(twoSixteen.mul(b6));

    const aNew = c4.neg().div(fortyEight);
    const bNew = c6.neg().div(eightSixtyFour);

    return new EllipticCurveFiniteField(field, aNew, bNew);
  }

  throw new ValueError(`Invalid coefficients: expected 2 or 5 values, got ${coeffs.length}`);
}

// ============================================================================
// Stub functions for EllipticCurveFiniteField
// ============================================================================

/**
 * Represents the structure of an abelian group.
 */
export interface AbelianGroupStructure {
  /** The invariant factors [n1, n2, ...] where n_{i+1} | n_i */
  invariants: bigint[];
  /** The generators corresponding to each invariant */
  generators: EllipticCurvePoint[];
  /** The total order of the group */
  order: bigint;
}

/**
 * Return the abelian group structure of the group of points on this curve.
 *
 * OUTPUT: An AbelianGroupStructure object containing:
 * - invariants: the invariant factors [n1, n2] where n2 | n1
 * - generators: points [P1, P2] such that E(F_q) = <P1> x <P2>
 * - order: the total order of the group
 *
 * The group E(F_q) is isomorphic to Z/n1 x Z/n2 where n2 | n1.
 * For a cyclic group, n2 = 1 and we have just one generator.
 *
 * ALGORITHM:
 * 1. Get generators from PARI via gens()
 * 2. Compute orders of generators
 * 3. Adjust to get a proper basis (using discrete log if needed)
 *
 * @see Reference: sage/schemes/elliptic_curves/ell_finite_field.py:abelian_group
 */
export function abelian_group(E: EllipticCurveFiniteField): AbelianGroupStructure {
  let gens = E.gens();
  const n = E.cardinality();

  if (gens.length === 2) {
    // Direct port of ell_finite_field.py:963-981.
    const P = gens[0]!;
    let Q = gens[1]!;
    const n1 = P.order();
    const n2 = n / n1;

    // PARI should guarantee this.
    if (!Q.mul(n1).isZero()) {
      throw new ArithmeticError('abelian_group(): n1*Q is not the identity');
    }

    // k = n1.prime_to_m_part(n2); Q *= k  -- kill the part we do not need.
    const k = prime_to_m_part(n1, n2);
    Q = Q.mul(k);

    const nQ = n2 * order_from_multiple_point(Q.mul(n2), n1 / k / n2);

    const S = P.mul(n / nQ);
    const T = Q.mul(n2);
    S.setOrder(nQ / n2, false); // for .log()
    const x = discrete_log(S, T, nQ / n2);
    Q = Q.sub(P.mul((((x * (n1 / nQ)) % n1) + n1) % n1));

    // by construction
    if (!Q.mul(n2).isZero()) {
      throw new ArithmeticError('abelian_group(): basis correction failed');
    }
    Q.setOrder(n2, false);

    gens = [P, Q];

    const invariants = [n1, n2];
    if (invariants[0]! * invariants[1]! !== n) {
      throw new ArithmeticError(
        `abelian_group(): invariants ${invariants} do not multiply to the cardinality ${n}`
      );
    }
    return { invariants, generators: gens, order: n };
  }

  if (gens.length === 0) {
    // Trivial group (shouldn't happen for curves over F_q with q > 1)
    return {
      invariants: [],
      generators: [],
      order: 1n,
    };
  }

  // Cyclic group: gens() has already asserted that P generates the whole group.
  const P = gens[0]!;
  const n1 = P.order();
  return {
    invariants: [n1],
    generators: [P],
    order: n1,
  };
}

/**
 * Return the largest divisor of n that is coprime to m.
 *
 * @see Reference: sage/rings/integer.pyx:prime_to_m_part
 */
function prime_to_m_part(n: bigint, m: bigint): bigint {
  if (n === 0n) {
    throw new ArithmeticError('self must be nonzero');
  }
  let k = n < 0n ? -n : n;
  for (;;) {
    const g = gcd(k, m);
    if (g === 1n) {
      return k;
    }
    k /= g;
  }
}

/**
 * Return the exact order of the point P, given that ``m`` is a multiple of it.
 *
 * @see Reference: sage/groups/generic.py:order_from_multiple
 */
function order_from_multiple_point(P: EllipticCurvePoint, m: bigint): bigint {
  if (m <= 0n) {
    throw new ValueError('multiple must be positive');
  }
  let ord = m;
  for (const [q, e] of factor(m)) {
    for (let i = 0n; i < e; i++) {
      if (P.mul(ord / q).isZero()) {
        ord /= q;
      } else {
        break;
      }
    }
  }
  return ord;
}

/**
 * Compute the discrete logarithm of Q with respect to P.
 *
 * INPUT:
 * - P: a point on the curve (the base)
 * - Q: a point on the curve (target: we want n such that Q = nP)
 * - ord: (optional) the order of P
 * - operation: (optional) the group operation ('+' for additive)
 *
 * OUTPUT: an integer n such that Q = nP, or throws if no such n exists
 *
 * ALGORITHM:
 * Baby-step giant-step algorithm:
 * 1. Compute baby steps: 0*P, 1*P, 2*P, ..., m*P where m = ceil(sqrt(order))
 * 2. Compute giant steps: Q - j*m*P for j = 0, 1, 2, ...
 * 3. Find a collision: if i*P = Q - j*m*P, then n = i + j*m
 *
 * Time complexity: O(sqrt(order)) group operations and space
 *
 * @see Reference: sage/schemes/elliptic_curves/ell_finite_field.py (via discrete_log)
 */
export function discrete_log(
  P: EllipticCurvePoint,
  Q: EllipticCurvePoint,
  ord?: bigint,
  _operation?: string
): bigint {
  // Handle trivial cases
  if (Q.isZero()) {
    return 0n;
  }

  if (P.isZero()) {
    throw new ValueError('base point P cannot be the identity');
  }

  // Get the order of P
  const orderP = ord !== undefined ? ord : P.order();

  if (orderP <= 0n) {
    throw new ValueError('order must be positive');
  }

  // Quick check: verify Q is in the subgroup generated by P
  if (!Q.mul(orderP).isZero()) {
    throw new ValueError('Q is not in the subgroup generated by P');
  }

  // Baby-step giant-step algorithm
  const m = isqrt(orderP) + 1n;

  // Baby steps: compute iP for i = 0, 1, ..., m-1
  // Store in a map: point -> index
  const babySteps = new Map<string, bigint>();
  let current = P.curve.zero();

  for (let i = 0n; i < m; i++) {
    const key = pointToKey(current);
    babySteps.set(key, i);
    current = current.add(P);
  }

  // Giant steps: compute Q - j*m*P for j = 0, 1, ...
  const mP = P.mul(m); // m*P
  const negMP = mP.neg(); // -m*P for fast subtraction

  let giantStep = Q;

  for (let j = 0n; j <= m; j++) {
    const key = pointToKey(giantStep);
    const i = babySteps.get(key);

    if (i !== undefined) {
      // Found: i*P = Q - j*m*P, so Q = (i + j*m)*P
      const n = (i + j * m) % orderP;
      return n;
    }

    giantStep = giantStep.add(negMP);
  }

  throw new ValueError('discrete log does not exist');
}

/**
 * Convert a point to a string key for Map storage.
 */
function pointToKey(P: EllipticCurvePoint): string {
  if (P.isZero()) {
    return 'O';
  }
  return `${P.x!.value},${P.y!.value}`;
}

/**
 * Return the embedding degree of this curve with respect to n.
 *
 * The embedding degree is the smallest positive integer k such that
 * n divides q^k - 1, where q is the size of the base field.
 *
 * This is important for pairing-based cryptography:
 * - The Weil and Tate pairings map to the multiplicative group of F_{q^k}
 * - k determines the security level of pairing-based protocols
 *
 * INPUT:
 * - E: an elliptic curve over F_q
 * - n: a positive integer (typically the order of a subgroup)
 *
 * OUTPUT: the embedding degree k, i.e., the multiplicative order of q mod n
 *
 * ALGORITHM:
 * Find the smallest k > 0 such that q^k ≡ 1 (mod n).
 * This is the multiplicative order of q in (Z/nZ)*.
 *
 * @see Reference: sage/schemes/elliptic_curves/ell_finite_field.py:embedding_degree
 */
export function embedding_degree(E: EllipticCurveFiniteField, n: IntegerLike): bigint {
  const nBig = toBigInt(n);
  if (nBig <= 0n) {
    throw new ValueError('n must be positive');
  }

  if (nBig === 1n) {
    return 1n;
  }

  const q = E.field.characteristic; // For prime fields, q = p

  // Compute the multiplicative order of q mod nBig
  // i.e., smallest k > 0 such that q^k ≡ 1 (mod nBig)

  // First check if gcd(q, nBig) = 1
  const g = gcd(q, nBig);
  if (g !== 1n) {
    // nBig divides some power of q, or shares a factor
    // The embedding degree may not be well-defined in the usual sense
    throw new ValueError('n must be coprime to the field characteristic');
  }

  // Find the multiplicative order
  let qPow = q % nBig;
  let k = 1n;

  // Maximum possible order is phi(nBig) <= nBig
  while (k <= nBig) {
    if (qPow === 1n) {
      return k;
    }
    qPow = (qPow * q) % nBig;
    k++;
  }

  // Should not reach here if gcd(q, nBig) = 1
  throw new ValueError('embedding degree not found (should not happen)');
}

/**
 * Return True if this elliptic curve is supersingular.
 *
 * An elliptic curve over a field of characteristic p is supersingular
 * if and only if its trace of Frobenius t satisfies t ≡ 0 (mod p).
 *
 * Equivalently, E is supersingular iff:
 * - E has no points of order p over the algebraic closure
 * - The endomorphism ring is an order in a quaternion algebra
 * - The p-rank is 0
 *
 * For curves over F_p, E is supersingular iff #E(F_p) = p + 1 (i.e., t = 0).
 *
 * INPUT:
 * - E: an elliptic curve over a finite field
 * - proof: if True (default), return a proved result
 *
 * OUTPUT: True if E is supersingular, False otherwise
 *
 * ALGORITHM:
 * Compute the trace of Frobenius and check if it is divisible by p.
 *
 * @see Reference: sage/schemes/elliptic_curves/ell_finite_field.py:is_supersingular
 */
export function is_supersingular(E: EllipticCurveFiniteField, proof: boolean = true): boolean {
  return is_j_supersingular(E.j_invariant(), proof);
}

/**
 * Return True if ``j`` is a supersingular j-invariant.
 *
 * INPUT:
 * - j: a finite field element
 * - proof: if True (default), return a proved result. If False, a return
 *   value of False is certain but True may be based on a probabilistic test.
 *
 * ALGORITHM:
 * j = 0 is supersingular iff p = 3 or p = 2 (mod 3); j = 1728 is
 * supersingular iff p = 2 or p = 3 (mod 4); for p in {2,3,5,7,11} there are
 * no other supersingular invariants. Otherwise, over GF(p) a supersingular
 * curve has exactly p+1 points, so we test (p+1)*P = 0 for random points and
 * return False as soon as one fails. When ``proof`` is set we finish with the
 * trace-of-Frobenius test.
 *
 * @see Reference: sage/schemes/elliptic_curves/ell_finite_field.py:is_j_supersingular
 * @see Deviation: is_j_supersingular precomputed j-polynomials
 */
export function is_j_supersingular(j: FieldElement, proof: boolean = true): boolean {
  const K = j.parent as BaseField;
  const p = K.characteristic;

  if (j.isZero()) {
    return p === 3n || p % 3n === 2n;
  }

  if (j.sub(K.__call__(1728n)).isZero()) {
    return p === 2n || p % 4n === 3n;
  }

  // From now on we know that j != 0, 1728.
  if (p === 2n || p === 3n || p === 5n || p === 7n || p === 11n) {
    return false; // since j = 0, 1728 are the only s.s. invariants
  }

  // Over GF(p) the minimal polynomial of j has degree 1, so no extension is
  // needed. (Sage additionally consults a table of supersingular
  // j-polynomials here; the probabilistic test below plus the trace check
  // give the same answer whenever ``proof`` is set.)
  const k = j.sub(K.__call__(1728n));
  const a = j.mul(k).mul(K.__call__(-3n));
  const b = j.mul(k).mul(k).mul(K.__call__(-2n));
  const Ej = new EllipticCurveFiniteField(K, a, b);

  for (let i = 0; i < 10; i++) {
    const P = Ej.random_point();
    if (!P.mul(p + 1n).isZero()) {
      return false;
    }
  }

  // When proof is False we return True for any curve which passes the
  // probabilistic test.
  if (!proof) {
    return true;
  }

  return Ej.trace_of_frobenius() % p === 0n;
}

/**
 * Return True if this elliptic curve is ordinary.
 *
 * An elliptic curve is ordinary if and only if it is not supersingular.
 * Equivalently, E is ordinary iff its trace of Frobenius t satisfies gcd(t, p) = 1.
 *
 * For ordinary curves:
 * - E has points of order p over the algebraic closure
 * - The endomorphism ring is an order in an imaginary quadratic field
 * - The p-rank is 1
 *
 * INPUT:
 * - E: an elliptic curve over a finite field
 * - proof: if True (default), return a proved result
 *
 * OUTPUT: True if E is ordinary, False otherwise
 *
 * @see Reference: sage/schemes/elliptic_curves/ell_finite_field.py:is_ordinary
 */
export function is_ordinary(E: EllipticCurveFiniteField, proof: boolean = true): boolean {
  return !is_j_supersingular(E.j_invariant(), proof);
}

/**
 * Return True if this curve is isogenous to another curve.
 *
 * Two elliptic curves over a finite field are isogenous if and only if
 * they have the same number of points (Tate's theorem).
 *
 * INPUT:
 * - E: an elliptic curve over a finite field
 * - other: another elliptic curve over a finite field
 * - field: (optional) field over which to check (not yet supported)
 * - proof: if True (default), return a proved result
 *
 * OUTPUT: True if E is isogenous to other
 *
 * ALGORITHM:
 * Compare the cardinalities of E and other. By Tate's theorem,
 * two curves over a finite field are isogenous iff they have the same
 * number of points.
 *
 * @see Reference: sage/schemes/elliptic_curves/ell_finite_field.py:is_isogenous
 */
export function is_isogenous(
  E: EllipticCurveFiniteField,
  other: EllipticCurveFiniteField,
  _field?: unknown,
  _proof?: boolean
): boolean {
  // Quick check: same curve
  if (E === other) {
    return true;
  }

  // Check same base field characteristic
  if (E.field.characteristic !== other.field.characteristic) {
    throw new ValueError('The base fields must have the same characteristic.');
  }

  // Main check: compare cardinalities (Tate's theorem)
  return E.cardinality() === other.cardinality();
}

/**
 * Find a point of exact order n on the curve.
 *
 * INPUT:
 * - E: an elliptic curve
 * - n: a positive integer
 *
 * OUTPUT: a point P such that P.order() == n, or null if no such point exists
 *
 * ALGORITHM:
 * 1. Compute the curve order N
 * 2. If n does not divide N, return null
 * 3. Compute cofactor c = N / n
 * 4. Sample random points and multiply by cofactor until we find one of exact order n
 *
 * @see Reference: sage/schemes/elliptic_curves/ell_finite_field.py
 */
export function _find_point_of_order(
  E: EllipticCurveFiniteField,
  n: bigint
): EllipticCurvePoint | null {
  if (n <= 0n) {
    throw new ValueError('n must be positive');
  }

  if (n === 1n) {
    return E.zero();
  }

  const N = E.cardinality();

  // n must divide N
  if (N % n !== 0n) {
    return null;
  }

  const cofactor = N / n;

  // Try to find a point of exact order n
  // Maximum attempts to prevent infinite loops
  const maxAttempts = 1000;

  for (let attempt = 0; attempt < maxAttempts; attempt++) {
    const R = E.random_point();
    const P = R.mul(cofactor);

    if (P.isZero()) {
      continue;
    }

    // Check if P has exact order n
    // P has order dividing n, we need to verify it's exactly n
    const ordP = P.order();

    if (ordP === n) {
      return P;
    }

    // If ordP divides n but is not n, try to adjust
    // This can happen when gcd(cofactor, n) > 1
    if (n % ordP === 0n) {
      // P has order ordP, we need order n = ordP * (n/ordP)
      // Need to find another point
      continue;
    }
  }

  return null;
}

/**
 * Check if two points P and Q are independent in the n-torsion subgroup.
 *
 * Two points are independent if neither is a multiple of the other.
 * More precisely, P and Q are independent in E[n] if they generate
 * a group of order n^2.
 *
 * INPUT:
 * - P: a point of order dividing n
 * - Q: a point of order dividing n
 * - n: positive integer
 *
 * OUTPUT: true if P and Q are independent
 *
 * ALGORITHM:
 * Check if the Weil pairing e_n(P, Q) is a primitive n-th root of unity.
 * For finite fields, this reduces to checking that the order of the pairing
 * result is n.
 *
 * Simplified check: verify that Q is not in <P> by checking that
 * there's no k such that Q = kP.
 */
export function _is_independent(P: EllipticCurvePoint, Q: EllipticCurvePoint, n: bigint): boolean {
  if (P.isZero() || Q.isZero()) {
    return false;
  }

  // Check that Q is not a multiple of P
  // This is done by trying to compute discrete log
  // If it fails, they are independent

  try {
    // Try to express Q as a multiple of P
    // If we find k such that Q = kP, they are not independent
    discrete_log(P, Q, n);
    return false; // Q = kP for some k, so they're dependent
  } catch {
    // discrete_log threw an error, meaning Q is not in <P>
    return true;
  }
}

/**
 * Return a basis for the n-torsion subgroup of this curve.
 *
 * E[n] = {P : [n]P = O} is the n-torsion subgroup.
 *
 * INPUT:
 * - E: an elliptic curve over a finite field
 * - n: a positive integer
 *
 * OUTPUT:
 * - If E[n] is cyclic (of order n): returns [P] where P has order n
 * - If E[n] ≅ Z/n × Z/n: returns [P, Q] where P, Q are independent points of order n
 *
 * ALGORITHM:
 * 1. Find a point P of order n (using _find_point_of_order)
 * 2. If full n-torsion is rational (n^2 divides #E), find an independent point Q
 * 3. Return the basis [P] or [P, Q]
 *
 * If the full n-torsion is not rational over the base field, this function
 * throws an error. Use division_field to extend to a field where full
 * n-torsion is rational.
 *
 * @see Reference: sage/schemes/elliptic_curves/ell_finite_field.py:torsion_basis
 *
 * @example
 * ```typescript
 * const F = GF(62207n ** 2n);
 * const E = EllipticCurve(F, [1n, 0n]);
 * const [P, Q] = torsion_basis(E, 256n);  // 2^8-torsion basis
 * ```
 */
export function torsion_basis(
  E: EllipticCurveFiniteField,
  n: bigint | number
): EllipticCurvePoint[] {
  const nBig = BigInt(n);

  if (nBig <= 0n) {
    throw new ValueError('n must be positive');
  }

  // Sage: T = self.abelian_group().torsion_subgroup(n); the basis exists iff
  // T.invariants() == (n, n).
  const G = abelian_group(E);
  const n1 = G.invariants[0] ?? 1n;
  const n2 = G.invariants[1] ?? 1n;

  if (gcd(nBig, n1) !== nBig || gcd(nBig, n2) !== nBig) {
    throw new ValueError(`curve does not have full rational ${nBig}-torsion`);
  }

  const P = G.generators[0]!.mul(n1 / nBig);
  const Q = G.generators[1]!.mul(n2 / nBig);
  P.setOrder(nBig, false);
  Q.setOrder(nBig, false);
  return [P, Q];
}

/**
 * Return all points in the n-torsion subgroup E[n].
 *
 * E[n] = {P : [n]P = O} is the set of all n-torsion points.
 *
 * INPUT:
 * - E: an elliptic curve over a finite field
 * - n: a positive integer
 *
 * OUTPUT: array of all points P such that [n]P = O
 *
 * ALGORITHM:
 * 1. Get a basis for E[n] using torsion_basis
 * 2. Enumerate all linear combinations of basis elements
 *
 * @see Reference: sage/schemes/elliptic_curves/ell_torsion.py
 *
 * @example
 * ```typescript
 * const F = GF(101n);
 * const E = EllipticCurve(F, [1n, 1n]);
 * const torsion2 = torsion_subgroup(E, 2n);  // All 2-torsion points
 * ```
 */
export function torsion_subgroup(
  E: EllipticCurveFiniteField,
  n: bigint | number
): EllipticCurvePoint[] {
  const nBig = BigInt(n);

  if (nBig <= 0n) {
    throw new ValueError('n must be positive');
  }

  if (nBig === 1n) {
    return [E.zero()];
  }

  // E(F_q) = Z/n1 x Z/n2, so E[n] = Z/d1 x Z/d2 with di = gcd(n, ni),
  // generated by (n1/d1)*G1 and (n2/d2)*G2.
  const G = abelian_group(E);
  const n1 = G.invariants[0] ?? 1n;
  const n2 = G.invariants[1] ?? 1n;
  const d1 = gcd(nBig, n1);
  const d2 = gcd(nBig, n2);

  const points: EllipticCurvePoint[] = [];
  const P = d1 > 1n ? G.generators[0]!.mul(n1 / d1) : E.zero();
  const Q = d2 > 1n ? G.generators[1]!.mul(n2 / d2) : E.zero();

  for (let a = 0n; a < d1; a++) {
    const aP = P.mul(a);
    for (let b = 0n; b < d2; b++) {
      points.push(aP.add(Q.mul(b)));
    }
  }

  return points;
}

/**
 * Compute modular inverse of a mod m.
 * Returns null if inverse doesn't exist.
 */
function modInverse(a: bigint, m: bigint): bigint | null {
  if (m <= 0n) return null;

  a = ((a % m) + m) % m;

  if (a === 0n) return null;

  // Extended Euclidean algorithm
  let [old_r, r] = [a, m];
  let [old_s, s] = [1n, 0n];

  while (r !== 0n) {
    const quotient = old_r / r;
    [old_r, r] = [r, old_r - quotient * r];
    [old_s, s] = [s, old_s - quotient * s];
  }

  if (old_r !== 1n) {
    return null; // gcd(a, m) != 1, no inverse
  }

  return ((old_s % m) + m) % m;
}

/**
 * Return all points Q such that [n]Q = P.
 *
 * This finds all "n-th roots" of the point P, i.e., all solutions
 * to the equation [n]Q = P.
 *
 * INPUT:
 * - P: a point on an elliptic curve
 * - n: a positive integer
 *
 * OUTPUT: array of all points Q on E such that [n]Q = P
 *
 * ALGORITHM:
 * 1. If P = O, return E[n] (all n-torsion points)
 * 2. Otherwise, find one solution Q_0 (if it exists)
 * 3. All solutions are Q_0 + T for T in E[n]
 *
 * @see Reference: sage/schemes/elliptic_curves/ell_point.py:division_points
 *
 * @example
 * ```typescript
 * const F = GF(101n);
 * const E = EllipticCurve(F, [1n, 1n]);
 * const P = E.random_point();
 * const roots = division_points(P, 3n);  // All Q such that [3]Q = P
 * for (const Q of roots) {
 *   console.assert(Q.mul(3n).eq(P));
 * }
 * ```
 */
export function division_points(P: EllipticCurvePoint, n: bigint | number): EllipticCurvePoint[] {
  const nBig = BigInt(n);
  const E = P.curve;

  if (nBig <= 0n) {
    throw new ValueError('n must be positive');
  }

  if (nBig === 1n) {
    return [P];
  }

  // Compute a polynomial g whose roots are the possible x-coordinates of the
  // n-division points of P, exactly as Sage does in
  // ell_point.py:1516-1528 -- but with the underlying division polynomials
  // supplied by PARI (``elldivpol`` / ``ellxn``), which is what PARI itself
  // uses for ``elldivpol``.
  const pariCurve = E.toPari();
  const K = E.field;
  const ans: EllipticCurvePoint[] = [];

  let gCoeffs: bigint[];
  if (P.isZero()) {
    ans.push(P);
    gCoeffs = elldivpol(pariCurve, Number(nBig));
  } else {
    const [num, den] = ellxn(pariCurve, Number(nBig));
    const xP = P.x!.value;
    gCoeffs = polySubScaled(num, den, xP, K.characteristic);
  }

  const nP = P.neg();
  const P_is_2_torsion = P.eq(nP);

  for (const x of rootsOverPrimeField(gCoeffs, K)) {
    if (!E.is_x_coord(x)) {
      continue;
    }
    const Q = E.lift_x(x);
    const nQ = Q.neg();
    const mQ = Q.mul(nBig);
    if (P_is_2_torsion) {
      // If P == -P then Q works iff -Q works, so include both unless equal.
      if (mQ.eq(P)) {
        ans.push(Q);
        if (!nQ.eq(Q)) {
          ans.push(nQ);
        }
      }
    } else {
      // P is not 2-torsion, so at most one of Q, -Q works.
      if (mQ.eq(P)) {
        ans.push(Q);
      } else if (mQ.eq(nP)) {
        ans.push(nQ);
      }
    }
  }

  // Sage sorts the result by (Z, X, Y): the identity first, then by x then y.
  ans.sort((A, B) => {
    if (A.isZero() !== B.isZero()) return A.isZero() ? -1 : 1;
    if (A.isZero()) return 0;
    if (A.x!.value !== B.x!.value) return A.x!.value < B.x!.value ? -1 : 1;
    if (A.y!.value !== B.y!.value) return A.y!.value < B.y!.value ? -1 : 1;
    return 0;
  });

  return ans;
}

/**
 * Return the coefficients of ``num - c * den`` modulo p.
 */
function polySubScaled(num: bigint[], den: bigint[], c: bigint, p: bigint): bigint[] {
  const len = Math.max(num.length, den.length);
  const out: bigint[] = new Array(len).fill(0n);
  for (let i = 0; i < len; i++) {
    const a = num[i] ?? 0n;
    const b = den[i] ?? 0n;
    out[i] = (((a - c * b) % p) + p) % p;
  }
  while (out.length > 0 && out[out.length - 1] === 0n) {
    out.pop();
  }
  return out;
}

/**
 * Return the distinct roots in the prime field K of the polynomial given by
 * ``coeffs`` (constant term first).
 */
function rootsOverPrimeField(coeffs: bigint[], K: BaseField): FieldElement[] {
  if (coeffs.length === 0) {
    throw new ValueError('roots of zero polynomial are not defined');
  }
  const R = new PolynomialRing(K as unknown as CoefficientRing<RingElement>, 'x');
  const f = R.__call__(coeffs.map((c) => K.__call__(c) as unknown as RingElement));
  const rts = f.roots().map(([r]: [RingElement, number]) => r as unknown as FieldElement);
  rts.sort((a, b) => (a.value < b.value ? -1 : a.value > b.value ? 1 : 0));
  return rts;
}

/**
 * Return the multiplication-by-m map as an isogeny object.
 *
 * This returns the endomorphism [m]: E -> E that maps P to [m]P.
 * The result is an "isogeny" object with methods for applying the map.
 *
 * INPUT:
 * - E: an elliptic curve over a finite field
 * - m: a nonzero integer
 *
 * OUTPUT: a ScalarMultiplicationIsogeny object
 *
 * @see Reference: sage/schemes/elliptic_curves/ell_generic.py:scalar_multiplication
 * @see Reference: sage/schemes/elliptic_curves/hom_scalar.py
 *
 * @example
 * ```typescript
 * const F = GF(101n);
 * const E = EllipticCurve(F, [1n, 1n]);
 * const doubling = multiplication_by_m_isogeny(E, 2n);
 * const P = E.random_point();
 * console.assert(doubling.call(P).eq(P.mul(2n)));
 * ```
 */
export function multiplication_by_m_isogeny(
  E: EllipticCurveFiniteField,
  m: bigint | number
): ScalarMultiplicationIsogeny {
  const mBig = BigInt(m);

  if (mBig === 0n) {
    throw new ValueError('m must be nonzero');
  }

  return new ScalarMultiplicationIsogeny(E, mBig);
}

/**
 * Represents the multiplication-by-m map [m]: E -> E as an isogeny.
 *
 * This class implements the scalar multiplication endomorphism,
 * which maps P to [m]P for all points P on E.
 *
 * Properties:
 * - Domain and codomain are the same curve E
 * - Degree is m^2 (as a morphism)
 * - Kernel is E[m] (the m-torsion subgroup)
 */
export class ScalarMultiplicationIsogeny {
  readonly domain: EllipticCurveFiniteField;
  readonly codomain: EllipticCurveFiniteField;
  readonly scalar: bigint;

  constructor(E: EllipticCurveFiniteField, m: bigint) {
    this.domain = E;
    this.codomain = E; // Endomorphism, so domain = codomain
    this.scalar = m;
  }

  /**
   * Apply the isogeny to a point.
   *
   * @param P - A point on the domain curve
   * @returns [m]P on the codomain curve
   */
  call(P: EllipticCurvePoint): EllipticCurvePoint {
    return P.mul(this.scalar);
  }

  /**
   * Apply the isogeny (alias for call).
   */
  __call__(P: EllipticCurvePoint): EllipticCurvePoint {
    return this.call(P);
  }

  /**
   * Return the degree of the isogeny.
   *
   * For multiplication-by-m, the degree is m^2.
   */
  degree(): bigint {
    return this.scalar * this.scalar;
  }

  /**
   * Return the kernel of the isogeny.
   *
   * For multiplication-by-m, the kernel is E[m].
   */
  kernel(): EllipticCurvePoint[] {
    const m = this.scalar >= 0n ? this.scalar : -this.scalar;
    return torsion_subgroup(this.domain, m);
  }

  /**
   * Check if a point is in the kernel.
   */
  isInKernel(P: EllipticCurvePoint): boolean {
    return this.call(P).isZero();
  }

  /**
   * String representation.
   */
  toString(): string {
    return `Scalar-multiplication endomorphism [${this.scalar}] of ${this.domain}`;
  }
}

/**
 * Compute the number of points on the curve over extension fields.
 *
 * INPUT:
 * - n: compute cardinalities over F_q, F_q^2, ..., F_q^n (default 1)
 *
 * OUTPUT:
 * - If n=1, returns a single bigint (the cardinality over the base field)
 * - If n>1, returns a list [c_1, c_2, ..., c_n] where c_d is the cardinality
 *   over the extension of degree d
 *
 * The count over extension fields uses the recurrence relation based on
 * the Frobenius polynomial: if f(x) = x^2 - t*x + q with roots alpha, beta,
 * then #E(F_{q^n}) = q^n + 1 - (alpha^n + beta^n).
 *
 * @see Reference: sage/schemes/elliptic_curves/ell_finite_field.py:count_points
 */
export function count_points(E: EllipticCurveFiniteField, n: number = 1): bigint | bigint[] {
  if (n < 1) {
    throw new ValueError('n must be a positive integer');
  }

  if (n === 1) {
    return E.cardinality();
  }

  const q = E.field.order;
  const t = E.trace_of_frobenius();

  // Use recurrence for trace of Frobenius powers:
  // t_1 = t, t_{n+1} = t * t_n - q * t_{n-1}
  // where #E(F_{q^n}) = q^n + 1 - t_n
  const result: bigint[] = [];
  let tPrev = 2n; // t_0 = 2 (trace of identity)
  let tCurr = t; // t_1 = t

  let qPower = q;
  for (let i = 1; i <= n; i++) {
    const count = qPower + 1n - tCurr;
    result.push(count);

    if (i < n) {
      // Update for next iteration
      const tNext = t * tCurr - q * tPrev;
      tPrev = tCurr;
      tCurr = tNext;
      qPower *= q;
    }
  }

  return result;
}

/**
 * Simple representation of a polynomial over ZZ.
 * Coefficients are stored as coeffs[i] = coefficient of x^i.
 */
export interface ZZPolynomial {
  /** Coefficients of the polynomial, coeffs[i] is coefficient of x^i */
  coeffs: readonly bigint[];
  /** Evaluate the polynomial at x */
  evaluate(x: bigint): bigint;
  /** Return the discriminant of this quadratic polynomial */
  discriminant(): bigint;
  /** String representation */
  toString(): string;
}

/**
 * Create a ZZ polynomial from coefficients.
 */
function makeZZPolynomial(coeffs: bigint[]): ZZPolynomial {
  // Remove trailing zeros
  while (coeffs.length > 0 && coeffs[coeffs.length - 1] === 0n) {
    coeffs.pop();
  }

  return {
    coeffs: Object.freeze([...coeffs]),
    evaluate(x: bigint): bigint {
      if (coeffs.length === 0) return 0n;
      // Horner's method
      let result = coeffs[coeffs.length - 1]!;
      for (let i = coeffs.length - 2; i >= 0; i--) {
        result = result * x + coeffs[i]!;
      }
      return result;
    },
    discriminant(): bigint {
      // For ax^2 + bx + c, discriminant is b^2 - 4ac
      if (coeffs.length !== 3) {
        throw new ValueError('discriminant only implemented for quadratic polynomials');
      }
      const a = coeffs[2]!;
      const b = coeffs[1]!;
      const c = coeffs[0]!;
      return b * b - 4n * a * c;
    },
    toString(): string {
      if (coeffs.length === 0) return '0';
      const terms: string[] = [];
      for (let i = coeffs.length - 1; i >= 0; i--) {
        const c = coeffs[i]!;
        if (c === 0n) continue;

        let term: string;
        if (i === 0) {
          if (terms.length > 0 && c > 0n) {
            term = `+ ${c}`;
          } else if (terms.length > 0 && c < 0n) {
            term = `- ${-c}`;
          } else {
            term = c.toString();
          }
        } else if (i === 1) {
          if (terms.length === 0) {
            if (c === 1n) term = 'x';
            else if (c === -1n) term = '-x';
            else term = `${c}*x`;
          } else if (c === 1n) {
            term = '+ x';
          } else if (c === -1n) {
            term = '- x';
          } else if (c > 0n) {
            term = `+ ${c}*x`;
          } else {
            term = `- ${-c}*x`;
          }
        } else {
          if (terms.length === 0) {
            if (c === 1n) term = `x^${i}`;
            else if (c === -1n) term = `-x^${i}`;
            else term = `${c}*x^${i}`;
          } else if (c === 1n) {
            term = `+ x^${i}`;
          } else if (c === -1n) {
            term = `- x^${i}`;
          } else if (c > 0n) {
            term = `+ ${c}*x^${i}`;
          } else {
            term = `- ${-c}*x^${i}`;
          }
        }

        terms.push(term);
      }
      return terms.join(' ') || '0';
    },
  };
}

/**
 * Return the Frobenius polynomial of this curve.
 *
 * The Frobenius polynomial is x^2 - t*x + q where t is the trace
 * of Frobenius and q is the cardinality of the base field.
 *
 * OUTPUT: a polynomial x^2 - t*x + q
 *
 * @see Reference: sage/schemes/elliptic_curves/ell_finite_field.py:frobenius_polynomial
 */
export function frobenius_polynomial(E: EllipticCurveFiniteField): ZZPolynomial {
  const t = E.trace_of_frobenius();
  const q = E.field.order;
  // x^2 - t*x + q = [q, -t, 1] (coeffs of x^0, x^1, x^2)
  return makeZZPolynomial([q, -t, 1n]);
}

/**
 * Represents a quadratic order Z[pi] where pi satisfies a quadratic polynomial.
 * For supersingular curves where Frobenius is an integer, this is just Z.
 */
export interface QuadraticOrder {
  /** The defining polynomial of the order (x - a for Z, or x^2 - t*x + q for quadratic) */
  definingPolynomial: ZZPolynomial;
  /** The degree (1 for Z, 2 for quadratic) */
  degree: number;
  /** The conductor of the order (1 if maximal) */
  conductor?: bigint;
  /** The discriminant of the order */
  discriminant: bigint;
  /** String representation */
  toString(): string;
}

/**
 * Return the quadratic order Z[phi] where phi is the Frobenius endomorphism.
 *
 * For most ordinary curves, this is a quadratic order in an imaginary quadratic field.
 * For supersingular curves where Frobenius is an integer, this returns Z (represented
 * as a degree-1 order).
 *
 * OUTPUT: a quadratic order
 *
 * @see Reference: sage/schemes/elliptic_curves/ell_finite_field.py:frobenius_order
 */
export function frobenius_order(E: EllipticCurveFiniteField): QuadraticOrder {
  const t = E.trace_of_frobenius();
  const q = E.field.order;
  const D = t * t - 4n * q; // Frobenius discriminant

  // Check if the Frobenius polynomial is a perfect square (supersingular case)
  // This happens when D = 0, meaning t^2 = 4q, so t = +/-2*sqrt(q)
  if (D === 0n) {
    // Frobenius is an integer: pi = t/2
    const pi = t / 2n;
    return {
      definingPolynomial: makeZZPolynomial([-pi, 1n]), // x - pi
      degree: 1,
      discriminant: 0n,
      toString(): string {
        return `Order generated by [] in Number Field in pi with defining polynomial x + ${-pi}`;
      },
    };
  }

  // Generic case: the quadratic order Z[pi] of discriminant D = t^2 - 4q.
  //
  // Its conductor f is the largest positive integer with f^2 | D and
  // D/f^2 == 0 or 1 (mod 4) -- the latter condition is what makes D/f^2 a
  // discriminant. Stripping *every* square factor of |D| (which is what a
  // naive loop does) reports f = 2 for the already-fundamental D = -20.
  let f = 1n;
  for (const [pr, e] of factor(D < 0n ? -D : D)) {
    f *= pr ** (e / 2n);
  }
  // D = f^2 * d; if d == 2 or 3 (mod 4) then d is not a discriminant, and the
  // only prime that can cause this is 2.
  const mod4 = (x: bigint): bigint => ((x % 4n) + 4n) % 4n;
  while (f > 1n && f % 2n === 0n && mod4(D / (f * f)) > 1n) {
    f /= 2n;
  }
  if (mod4(D / (f * f)) > 1n) {
    throw new ArithmeticError(`frobenius_order: ${D} is not a discriminant`);
  }

  return {
    definingPolynomial: makeZZPolynomial([q, -t, 1n]), // x^2 - t*x + q
    degree: 2,
    conductor: f,
    discriminant: D,
    toString(): string {
      return `Order of conductor ${f} generated by pi in Number Field in pi with defining polynomial x^2 - ${t}*x + ${q}`;
    },
  };
}

/**
 * Represents the Frobenius endomorphism of an elliptic curve.
 *
 * The Frobenius endomorphism pi maps (x, y) -> (x^q, y^q) where q is the
 * cardinality of the base field.
 */
export class FrobeniusEndomorphism {
  readonly curve: EllipticCurveFiniteField;
  readonly degree: bigint;

  constructor(curve: EllipticCurveFiniteField) {
    this.curve = curve;
    this.degree = curve.field.order;
  }

  /**
   * Return the domain of this morphism.
   */
  domain(): EllipticCurveFiniteField {
    return this.curve;
  }

  /**
   * Return the codomain of this morphism.
   */
  codomain(): EllipticCurveFiniteField {
    return this.curve;
  }

  /**
   * Apply the Frobenius endomorphism to a point.
   *
   * For a point (x, y), returns (x^q, y^q).
   */
  call(P: EllipticCurvePoint): EllipticCurvePoint {
    if (P.isInfinity) {
      return EllipticCurvePoint.infinity(this.curve);
    }

    const q = this.degree;
    const xNew = P.x!.pow(q);
    const yNew = P.y!.pow(q);

    return new EllipticCurvePoint(this.curve, xNew, yNew);
  }

  /**
   * Return the degree of this endomorphism.
   */
  getDegree(): bigint {
    return this.degree;
  }

  /**
   * Check if this is a separable morphism.
   * The Frobenius is always inseparable (except in characteristic 0).
   */
  is_separable(): boolean {
    return false;
  }

  /**
   * Return the characteristic polynomial of the Frobenius.
   */
  characteristic_polynomial(): ZZPolynomial {
    return frobenius_polynomial(this.curve);
  }

  /**
   * String representation.
   */
  toString(): string {
    return `Frobenius endomorphism of degree ${this.degree}:\n  From: ${this.curve}\n  To:   ${this.curve}`;
  }
}

/**
 * Return the Frobenius endomorphism of this curve.
 *
 * The Frobenius endomorphism is the map pi: (x,y) -> (x^q, y^q) where q
 * is the cardinality of the base field.
 *
 * OUTPUT: the Frobenius endomorphism
 *
 * @see Reference: sage/schemes/elliptic_curves/ell_finite_field.py:frobenius_endomorphism
 */
export function frobenius_endomorphism(E: EllipticCurveFiniteField): FrobeniusEndomorphism {
  return new FrobeniusEndomorphism(E);
}

/**
 * Return the discriminant of the Frobenius endomorphism.
 *
 * The discriminant is t^2 - 4q where t is the trace of Frobenius
 * and q is the cardinality of the base field.
 *
 * OUTPUT: an integer D = t^2 - 4q
 *
 * @see Reference: sage/schemes/elliptic_curves/ell_finite_field.py:frobenius_discriminant
 */
export function frobenius_discriminant(E: EllipticCurveFiniteField): bigint {
  const t = E.trace_of_frobenius();
  const q = E.field.order;
  return t * t - 4n * q;
}

/**
 * Return the height above the floor in the ell-isogeny volcano.
 *
 * The ell-isogeny graph of elliptic curves over a finite field has a
 * "volcano" structure when organized by endomorphism ring. The "floor"
 * consists of curves with maximal endomorphism ring, and curves "above"
 * the floor have smaller endomorphism rings (larger conductors).
 *
 * INPUT:
 * - E: an elliptic curve over a finite field
 * - ell: a prime number (the isogeny degree)
 * - e: a positive integer (power of ell in the conductor)
 *
 * OUTPUT: the height h >= 0 such that E is at height h in the ell-volcano
 *
 * ALGORITHM:
 * The height is determined by the valuation of ell in the conductor of
 * the endomorphism ring. This requires computing the endomorphism ring
 * discriminant and factoring it.
 *
 * NOTE: This function requires computing the endomorphism ring, which
 * is a complex operation involving class field theory.
 *
 * @see Reference: sage/schemes/elliptic_curves/ell_finite_field.py:height_above_floor
 * @see [Koh1996] D. Kohel, "Endomorphism rings of elliptic curves over finite fields"
 */
export function height_above_floor(
  _E: EllipticCurveFiniteField,
  _ell: bigint | number,
  _e: number
): number {
  // This requires computing the endomorphism ring structure
  // which involves complex class field theory computations
  throw new NotImplementedError(
    'height_above_floor requires endomorphism ring computation. ' +
      'See endomorphism_order() for related functionality.'
  );
}

/**
 * Return the discriminant of the endomorphism ring given the class number.
 *
 * For an ordinary elliptic curve E over F_q, the endomorphism ring End(E)
 * is an order in an imaginary quadratic field. The discriminant D of this
 * order determines the endomorphism ring up to isomorphism.
 *
 * Given the class number h of the endomorphism ring, this function computes
 * the discriminant D.
 *
 * INPUT:
 * - E: an elliptic curve over a finite field
 * - h: the class number of the endomorphism ring
 *
 * OUTPUT: the discriminant D of the endomorphism ring
 *
 * ALGORITHM:
 * The class number h(D) can be computed from D using the class number formula.
 * This function inverts that relationship: given h, find D such that h(D) = h
 * and D is compatible with the Frobenius discriminant.
 *
 * NOTE: This is an advanced function that requires:
 * - Computing the Frobenius discriminant D_pi = t^2 - 4q
 * - Finding the fundamental discriminant D_0 dividing D_pi
 * - Determining the conductor f such that D = f^2 * D_0
 * - Verifying h(D) = h
 *
 * @see Reference: sage/schemes/elliptic_curves/ell_finite_field.py:endomorphism_discriminant_from_class_number
 */
export function endomorphism_discriminant_from_class_number(
  E: EllipticCurveFiniteField,
  _h: bigint | number
): bigint {
  // For ordinary curves, D_pi = t^2 - 4q where t is trace of Frobenius
  const t = E.trace_of_frobenius();
  const q = E.field.order;
  const D_pi = t * t - 4n * q;

  // The endomorphism ring discriminant D divides D_pi
  // and satisfies h(D) = h (the given class number)

  // Computing the exact discriminant from the class number requires
  // solving the inverse class number problem, which is computationally
  // intensive and requires class field theory.
  throw new NotImplementedError(
    'endomorphism_discriminant_from_class_number requires class number computation. ' +
      `Frobenius discriminant is ${D_pi}. Full implementation needs class field theory.`
  );
}

/**
 * Return the endomorphism ring of this curve as a quadratic order.
 *
 * For an ordinary elliptic curve E over a finite field F_q, the
 * endomorphism ring End(E) is an order O in an imaginary quadratic field K.
 * This order is determined by its discriminant D, and O = Z[(D + sqrt(D))/2].
 *
 * For supersingular curves, the endomorphism ring is a maximal order in
 * a quaternion algebra, which is more complex to represent.
 *
 * OUTPUT: a QuadraticOrder representing the endomorphism ring
 *
 * ALGORITHM:
 * 1. Compute the Frobenius polynomial x^2 - t*x + q
 * 2. The Frobenius generates an order Z[pi] with discriminant D_pi = t^2 - 4*q
 * 3. End(E) contains Z[pi] and may be larger (smaller discriminant)
 * 4. Use isogeny volcano structure to determine the exact order
 *
 * NOTE: This is a partial implementation that returns the Frobenius order.
 * The full endomorphism ring may be larger (when E is on the floor of the
 * isogeny volcano). Computing the exact endomorphism ring requires walking
 * the isogeny graph.
 *
 * @see Reference: sage/schemes/elliptic_curves/ell_finite_field.py:endomorphism_order
 * @see frobenius_order for the order generated by Frobenius
 */
export function endomorphism_order(E: EllipticCurveFiniteField): QuadraticOrder {
  // For supersingular curves, the endomorphism ring is a quaternion order
  if (is_supersingular(E)) {
    throw new NotImplementedError(
      'endomorphism_order for supersingular curves requires quaternion algebra infrastructure. ' +
        'The endomorphism ring is a maximal order in the quaternion algebra ramified at p and infinity.'
    );
  }

  // For ordinary curves, the endomorphism ring is an order in an imaginary quadratic field
  // The Frobenius order Z[pi] is always contained in End(E)
  // For curves on the floor of the isogeny volcano, End(E) = Z[pi]
  // For curves above the floor, End(E) may be strictly larger than Z[pi]

  // Return the Frobenius order as a lower bound
  // The actual endomorphism ring computation would require walking the isogeny graph
  const order = frobenius_order(E);

  // Note: This may not be the full endomorphism ring if E is above the floor
  // A complete implementation would use the algorithm from Kohel's thesis
  return order;
}

/**
 * Return the quadratic twist of this curve by D.
 *
 * For a curve y^2 = x^3 + ax + b over a field of characteristic != 2, 3,
 * the quadratic twist by D is y^2 = x^3 + aD^2*x + bD^3.
 *
 * If D is not specified, a non-square element is used.
 *
 * @param E - The elliptic curve
 * @param D - The twisting parameter (if not specified, a non-square is found)
 * @returns The quadratic twist of E by D
 *
 * @see Reference: sage/schemes/elliptic_curves/ell_field.py:quadratic_twist
 */
export function quadratic_twist(
  E: EllipticCurveFiniteField,
  D?: IntegerLike | FieldElement
): EllipticCurveFiniteField {
  const K = E.field;
  const p = K.characteristic;

  // For char 2 or 3, we don't support this yet
  if (p === 2n || p === 3n) {
    throw new NotImplementedError('quadratic_twist not implemented for characteristic 2 or 3');
  }

  let dElem: FieldElement;
  if (D === undefined) {
    // Find a non-square element
    const q = K.order;
    const q2 = (q - 1n) / 2n;
    dElem = K.__call__(K.gen().value);
    // Keep trying until we find a non-square
    let attempts = 0;
    while (dElem.isZero() || dElem.pow(q2).value === 1n) {
      // Try random elements
      dElem = K.__call__((BigInt(attempts + 2) * K.gen().value) % p);
      attempts++;
      if (attempts > 100) {
        // Fallback: find any non-square
        for (let i = 2n; i < p; i++) {
          const candidate = K.__call__(i);
          if (!candidate.isZero() && candidate.pow(q2).value !== 1n) {
            dElem = candidate;
            break;
          }
        }
        break;
      }
    }
  } else {
    // Check if D is a FieldElement (has isZero method) or an IntegerLike
    dElem =
      typeof D === 'object' && D !== null && 'isZero' in D
        ? (D as FieldElement)
        : K.__call__(toBigInt(D as IntegerLike));
  }

  if (dElem.isZero()) {
    throw new ValueError('twisting parameter D must be nonzero');
  }

  // For short Weierstrass form y^2 = x^3 + ax + b,
  // the quadratic twist by D is y^2 = x^3 + a*D^2*x + b*D^3
  const aNew = E.a.mul(dElem.pow(2));
  const bNew = E.b.mul(dElem.pow(3));

  return new EllipticCurveFiniteField(K, aNew, bNew);
}

/**
 * Return a complete list of pairwise nonisomorphic elliptic curves
 * with j-invariant 0 over the finite field K.
 *
 * For q = 1 (mod 6), there are 6 sextic twists.
 * For q = 5 (mod 6), there are 2 quadratic twists.
 *
 * @param K - The finite field
 * @returns List of curves with j-invariant 0
 *
 * @see Reference: sage/schemes/elliptic_curves/ell_finite_field.py:curves_with_j_0
 */
export function curves_with_j_0(K: BaseField): EllipticCurveFiniteField[] {
  const p = K.characteristic;

  if (p === 2n || p === 3n) {
    throw new NotImplementedError('curves_with_j_0 not implemented for characteristic 2 or 3');
  }

  const q = K.order;

  // For q = 2 (mod 3), we only have two quadratic twists
  // -3 is non-square in this case
  if (q % 3n === 2n) {
    return [
      new EllipticCurveFiniteField(K, 0n, 1n),
      new EllipticCurveFiniteField(K, 0n, K.__call__(-27n).value),
    ];
  }

  // For q = 1 (mod 3), we have genuine sextic twists
  // Find D generating K* mod 6th powers
  const q2 = (q - 1n) / 2n;
  const q3 = (q - 1n) / 3n;

  // Find a generator of K* mod 6th powers
  let D = K.__call__(2n);
  let attempts = 0;
  while (D.isZero() || D.pow(q2).value === 1n || D.pow(q3).value === 1n) {
    attempts++;
    D = K.__call__(BigInt(attempts + 2));
    if (attempts > 1000) {
      throw new ArithmeticError('Could not find generator for sextic twists');
    }
  }

  // Create 6 curves: y^2 = x^3 + D^i for i = 0, ..., 5
  const curves: EllipticCurveFiniteField[] = [];
  for (let i = 0; i < 6; i++) {
    const b = D.pow(BigInt(i));
    curves.push(new EllipticCurveFiniteField(K, 0n, b));
  }

  return curves;
}

/**
 * Return a complete list of pairwise nonisomorphic elliptic curves
 * with j-invariant 1728 over the finite field K.
 *
 * For q = 1 (mod 4), there are 4 quartic twists.
 * For q = 3 (mod 4), there are 2 quadratic twists.
 *
 * @param K - The finite field
 * @returns List of curves with j-invariant 1728
 *
 * @see Reference: sage/schemes/elliptic_curves/ell_finite_field.py:curves_with_j_1728
 */
export function curves_with_j_1728(K: BaseField): EllipticCurveFiniteField[] {
  const p = K.characteristic;

  if (p === 2n || p === 3n) {
    throw new NotImplementedError('curves_with_j_1728 not implemented for characteristic 2 or 3');
  }

  const q = K.order;

  // For q = 3 (mod 4), we only have two quadratic twists
  // -1 is non-square in this case
  if (q % 4n === 3n) {
    return [
      new EllipticCurveFiniteField(K, 1n, 0n),
      new EllipticCurveFiniteField(K, K.__call__(-1n).value, 0n),
    ];
  }

  // For q = 1 (mod 4), we have genuine quartic twists
  // Find D generating K* mod 4th powers
  const q2 = (q - 1n) / 2n;

  let D = K.__call__(2n);
  let attempts = 0;
  while (D.isZero() || D.pow(q2).value === 1n) {
    attempts++;
    D = K.__call__(BigInt(attempts + 2));
    if (attempts > 1000) {
      throw new ArithmeticError('Could not find generator for quartic twists');
    }
  }

  // Create 4 curves: y^2 = x^3 + D^i * x for i = 0, ..., 3
  const curves: EllipticCurveFiniteField[] = [];
  for (let i = 0; i < 4; i++) {
    const a = D.pow(BigInt(i));
    curves.push(new EllipticCurveFiniteField(K, a, 0n));
  }

  return curves;
}

/**
 * Return a list of all twists of this curve.
 *
 * A twist of E/k is an elliptic curve E' defined over k that is
 * isomorphic to E over the algebraic closure of k.
 *
 * Most elliptic curves over a finite field only admit a single
 * nontrivial twist (the quadratic twist); the only exceptions are
 * curves with j-invariant 0 or 1728.
 *
 * OUTPUT: a list of elliptic curves that are twists of self
 *
 * @see Reference: sage/schemes/elliptic_curves/ell_finite_field.py:twists
 */
export function twists(E: EllipticCurveFiniteField): EllipticCurveFiniteField[] {
  const K = E.field;
  const p = K.characteristic;
  const j = E.j_invariant();

  let allTwists: EllipticCurveFiniteField[] | null = null;
  if (j.isZero()) {
    allTwists = curves_with_j_0(K);
  } else if (j.eq(K.__call__(1728n))) {
    allTwists = curves_with_j_1728(K);
  }

  if (allTwists !== null) {
    // Sage (ell_finite_field.py:1939-1945) replaces the entry isomorphic to
    // self with twists[0] and puts self first. Note that upstream's ``break``
    // sits at the loop level, so only index 0 is ever examined; we reproduce
    // that here rather than inventing a different list.
    const t0 = allTwists[0]!;
    if (is_isomorphic(E, t0)) {
      allTwists[0] = E;
    }
    return allTwists;
  }

  // Now j is not 0 or 1728, and we only have a quadratic twist.
  if (p === 2n) {
    throw new NotImplementedError('twists not implemented for characteristic 2');
  }

  // Find a nonsquare D.
  const q = K.order;
  const q2 = (q - 1n) / 2n;
  let D = K.__call__(0n);
  for (let i = 2n; i < q; i++) {
    const cand = K.__call__(i);
    if (!cand.isZero() && cand.pow(q2).value !== 1n) {
      D = cand;
      break;
    }
  }
  if (D.isZero()) {
    throw new ArithmeticError('could not find a quadratic non-residue');
  }

  return [E, quadratic_twist(E, D)];
}

/**
 * Return whether the two short Weierstrass curves y^2 = x^3 + a*x + b are
 * isomorphic over their (common) base field.
 *
 * E and E' are isomorphic iff there exists u != 0 with a' = u^4*a, b' = u^6*b.
 *
 * @see Reference: sage/schemes/elliptic_curves/ell_generic.py:is_isomorphic
 */
export function is_isomorphic(E: EllipticCurveFiniteField, F: EllipticCurveFiniteField): boolean {
  const K = E.field;
  if (K.order !== F.field.order) {
    return false;
  }
  const q = K.order;

  const a = E.a;
  const b = E.b;
  const a2 = F.a;
  const b2 = F.b;

  // j = 0 (a = 0): need u^6 = b/b2.
  if (a.isZero() || a2.isZero()) {
    if (!a.isZero() || !a2.isZero()) {
      return false;
    }
    return hasRootOfUnityPower(b.div(b2), 6n, q);
  }
  // j = 1728 (b = 0): need u^4 = a/a2.
  if (b.isZero() || b2.isZero()) {
    if (!b.isZero() || !b2.isZero()) {
      return false;
    }
    return hasRootOfUnityPower(a.div(a2), 4n, q);
  }
  // Generic j: the only candidate is u^2 = (b*a2)/(b2*a); it must satisfy
  // u^4 = a/a2 and must itself be a square in K (otherwise u is only defined
  // over the quadratic extension, i.e. the curves are quadratic twists).
  const u2 = b.mul(a2).div(b2.mul(a));
  if (!u2.pow(2).eq(a.div(a2))) {
    return false;
  }
  return hasRootOfUnityPower(u2, 2n, q);
}

/**
 * Return whether ``c`` is an m-th power in GF(q)^*.
 */
function hasRootOfUnityPower(c: FieldElement, m: bigint, q: bigint): boolean {
  if (c.isZero()) {
    return false;
  }
  const d = gcd(m, q - 1n);
  return c.pow((q - 1n) / d).eq(c.parent.one() as FieldElement);
}

/**
 * Return the multiplication-by-p isogeny, where p is the characteristic.
 *
 * For an elliptic curve E over a finite field F_q of characteristic p,
 * the multiplication-by-p map [p]: E -> E is an isogeny of degree p^2.
 *
 * This isogeny factors as [p] = pi^ * pi where pi is the Frobenius
 * endomorphism of degree q and pi^ is its dual.
 *
 * OUTPUT: a ScalarMultiplicationIsogeny representing [p]: E -> E
 *
 * ALGORITHM:
 * In SageMath, this is computed as frob.dual() * frob where frob is the
 * Frobenius isogeny. Since we don't have full isogeny composition
 * infrastructure, we return a ScalarMultiplicationIsogeny which achieves
 * the same result functionally.
 *
 * @see Reference: sage/schemes/elliptic_curves/ell_finite_field.py:multiplication_by_p_isogeny
 *
 * @example
 * ```typescript
 * const p = 23n;
 * const K = GF(p ** 3n);
 * const E = EllipticCurve(K, [1n, 2n]);
 * const phi = multiplication_by_p_isogeny(E);
 * const P = E.random_point();
 * // phi(P) equals P.mul(p)
 * console.assert(phi.call(P).eq(P.mul(p)));
 * // Degree is p^2
 * console.assert(phi.degree() === p * p);
 * ```
 */
export function multiplication_by_p_isogeny(
  E: EllipticCurveFiniteField
): ScalarMultiplicationIsogeny {
  const p = E.field.characteristic;
  return new ScalarMultiplicationIsogeny(E, p);
}

/**
 * Return the j-invariants of curves that are l-isogenous to E.
 *
 * For an elliptic curve E with j-invariant j, this returns all j-invariants j'
 * such that there exists an isogeny of degree l from E to a curve with j-invariant j'.
 *
 * These are the j-invariants connected to j in the l-isogeny graph.
 *
 * The number of such j-invariants is:
 * - l+1 if j != 0, 1728 (for prime l)
 * - Varies for j = 0 or 1728 depending on the field
 *
 * ALGORITHM:
 * For SIDH/SIKE applications, this computes codomains of isogenies with
 * kernel generated by points of order l. For small l, we can enumerate
 * l-torsion points and compute the corresponding isogeny codomains.
 *
 * @param E - The elliptic curve
 * @param l - The prime degree of the isogeny
 * @returns Array of j-invariants of l-isogenous curves
 *
 * @see Reference: SIDH/SIKE key exchange uses this for isogeny graph walks
 */
export function j_invariant_neighbors(E: EllipticCurveFiniteField, l: bigint): FieldElement[] {
  if (l <= 1n) {
    throw new ValueError('l must be a prime > 1');
  }

  if (!is_prime(l)) {
    throw new ValueError('l must be prime');
  }

  const K = E.field;
  const p = K.characteristic;

  if (l === p) {
    // The Frobenius isogeny case
    // For a curve over F_p, the Frobenius is an isogeny of degree p
    // The codomain has the same j-invariant
    return [E.j_invariant()];
  }

  // For small l, we can find l-isogenous curves by computing isogenies
  // whose kernel is generated by a point of order l.
  //
  // The approach for SIDH/SIKE:
  // 1. Find all points of order l on E (over the base field or extension)
  // 2. For each cyclic subgroup of order l, compute the isogeny codomain
  // 3. Collect the distinct j-invariants

  const jInvariants: FieldElement[] = [];
  const seenJ = new Set<string>();

  // For prime l, the l-torsion E[l] is either:
  // - Trivial (if l | p-1 or l | p+1 depending on trace)
  // - Cyclic of order l
  // - Z/lZ x Z/lZ (full l-torsion)
  //
  // Over F_p, we can only access the F_p-rational l-torsion.

  // Try to find l-torsion points
  const N = E.cardinality();
  if (N % l !== 0n) {
    // No l-torsion over this field
    return [];
  }

  const cofactor = N / l;

  // Try random points to find l-torsion
  const maxAttempts = Number(l) * 10 + 100;
  const foundGenerators: EllipticCurvePoint[] = [];

  for (
    let attempt = 0;
    attempt < maxAttempts && foundGenerators.length < Number(l) + 1;
    attempt++
  ) {
    const P = E.random_point();
    const Q = P.mul(cofactor);

    if (!Q.isInfinity) {
      // Q has order dividing l, and l is prime, so order(Q) = l
      // Check if Q generates a new subgroup
      let isNew = true;
      for (const gen of foundGenerators) {
        // Check if Q is a multiple of gen
        for (let k = 1n; k < l; k++) {
          if (gen.mul(k).eq(Q)) {
            isNew = false;
            break;
          }
        }
        if (!isNew) break;
      }

      if (isNew) {
        foundGenerators.push(Q);

        // Compute the j-invariant of the isogeny codomain
        // Using Velu's formulas
        const jCodomain = compute_isogeny_j_invariant(E, Q, l);
        const jStr = jCodomain.value.toString();

        if (!seenJ.has(jStr)) {
          seenJ.add(jStr);
          jInvariants.push(jCodomain);
        }
      }
    }
  }

  return jInvariants;
}

/**
 * Compute the j-invariant of the codomain of an isogeny with given kernel generator.
 *
 * Using Velu's formulas: for an isogeny phi: E -> E' with kernel <P>,
 * we compute the curve E' and return its j-invariant.
 *
 * @param E - The source curve
 * @param P - A generator of the kernel (point of prime order l)
 * @param l - The order of P (degree of the isogeny)
 * @returns The j-invariant of the codomain curve
 */
function compute_isogeny_j_invariant(
  E: EllipticCurveFiniteField,
  P: EllipticCurvePoint,
  l: bigint
): FieldElement {
  const K = E.field;
  const a = E.a;
  const b = E.b;

  // Velu's formulas, exactly as in
  // ell_curve_isogeny.py:__update_kernel_data (lines 1949-1965):
  //
  //   gxQ = 3*xQ^2 + a4          (a1 = a2 = 0 here)
  //   gyQ = -2*yQ
  //   uQ  = gyQ^2
  //   vQ  = gxQ                  if Q is 2-torsion
  //       = 2*gxQ                otherwise
  //   v  += vQ ;  w += uQ + xQ*vQ
  //
  // summed over one representative of each pair {Q, -Q} in ker \ {O}.
  let v = K.__call__(0n);
  let w = K.__call__(0n);

  const three = K.__call__(3n);
  const two = K.__call__(2n);

  const halfL = l / 2n; // number of {Q,-Q} classes in <P> \ {O}
  for (let k = 1n; k <= halfL; k++) {
    const Q = P.mul(k);
    if (Q.isInfinity) continue;

    const xQ = Q.x!;
    const yQ = Q.y!;

    const gQx = three.mul(xQ.pow(2)).add(a);
    const gQy = two.neg().mul(yQ);
    const uQ = gQy.mul(gQy);

    // Q is 2-torsion iff 2*yQ == 0.
    const vQ = two.mul(yQ).isZero() ? gQx : two.mul(gQx);

    v = v.add(vQ);
    w = w.add(uQ.add(xQ.mul(vQ)));
  }

  // New curve coefficients (short Weierstrass form)
  // a' = a - 5*v
  // b' = b - 7*w
  const five = K.__call__(5n);
  const seven = K.__call__(7n);
  const aNew = a.sub(five.mul(v));
  const bNew = b.sub(seven.mul(w));

  // Compute j-invariant of the new curve
  // j = 1728 * 4a'^3 / (4a'^3 + 27b'^2)
  const four = K.__call__(4n);
  const twentySeven = K.__call__(27n);
  const val1728 = K.__call__(1728n);

  const fourANewCubed = four.mul(aNew.pow(3));
  const denominator = fourANewCubed.add(twentySeven.mul(bNew.pow(2)));

  if (denominator.isZero()) {
    // Singular curve - shouldn't happen for valid isogenies
    throw new ArithmeticError('Computed singular codomain curve');
  }

  return val1728.mul(fourANewCubed).div(denominator);
}
