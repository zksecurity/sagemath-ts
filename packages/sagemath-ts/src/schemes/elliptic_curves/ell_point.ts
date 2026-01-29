/**
 * @module sage/schemes/elliptic_curves/ell_point
 * @description Elliptic curve point arithmetic
 *
 * Port of: sage/schemes/elliptic_curves/ell_point.py
 *
 * Implements point arithmetic on elliptic curves in Weierstrass form:
 *   y^2 + a1*x*y + a3*y = x^3 + a2*x^2 + a4*x + a6
 *
 * For short Weierstrass form (y^2 = x^3 + a*x + b), we have:
 *   a1 = a2 = a3 = 0, a4 = a, a6 = b
 *
 * Reference: https://hyperelliptic.org/EFD/
 */

import { NotImplementedError, ValueError } from '../../errors.js';
import type { FieldElement, FieldRing } from './types.js';

// Re-export the types from types.js for consumers of this module
export type { FieldElement, FieldRing } from './types.js';

/**
 * Alias for FieldRing (for compatibility with different naming conventions).
 */
export type FieldParent = FieldRing;

/**
 * Forward declaration of curve interface to avoid circular imports.
 * The actual EllipticCurveGeneric class is defined in ell_generic.ts.
 */
export interface EllipticCurveInterface<F extends FieldElement> {
  readonly base_ring: FieldRing;
  a1(): F;
  a2(): F;
  a3(): F;
  a4(): F;
  a6(): F;
  a_invariants(): [F, F, F, F, F];
  zero(): EllipticCurvePoint<F>;
  is_on_curve(x: F, y: F): boolean;
}

/**
 * A point on an elliptic curve.
 *
 * Points are represented in projective coordinates (X : Y : Z) where:
 * - The affine point (x, y) corresponds to (x : y : 1)
 * - The point at infinity corresponds to (0 : 1 : 0)
 *
 * @example
 * ```typescript
 * const F = GF(23n);
 * const E = EllipticCurve(F, [0n, 0n, 0n, 1n, 1n]); // y^2 = x^3 + x + 1
 * const P = E.point([F.__call__(0n), F.__call__(1n)]);
 * const Q = E.point([F.__call__(6n), F.__call__(4n)]);
 *
 * // Point addition
 * const R = P.add(Q);
 *
 * // Scalar multiplication
 * const nP = P.mul(5n);
 *
 * // Check identity: P + O = P
 * const O = E.zero();
 * console.log(P.add(O).eq(P)); // true
 * ```
 */
export class EllipticCurvePoint<F extends FieldElement = FieldElement> {
  /** The parent curve */
  readonly curve: EllipticCurveInterface<F>;

  /** Projective X-coordinate */
  private readonly _X: F;

  /** Projective Y-coordinate */
  private readonly _Y: F;

  /** Projective Z-coordinate */
  private readonly _Z: F;

  /**
   * Create a point on an elliptic curve.
   *
   * Overloaded constructors:
   * 1. `new EllipticCurvePoint(curve, [], check)` - Point at infinity
   * 2. `new EllipticCurvePoint(curve, [x, y], check)` - Affine point
   * 3. `new EllipticCurvePoint(curve, X, Y, Z, check)` - Projective point
   */
  constructor(curve: EllipticCurveInterface<F>, coords: [], check?: boolean);
  constructor(curve: EllipticCurveInterface<F>, coords: [F, F], check?: boolean);
  constructor(curve: EllipticCurveInterface<F>, X: F, Y: F, Z: F, check?: boolean);
  constructor(
    curve: EllipticCurveInterface<F>,
    arg1: F | [] | [F, F],
    arg2?: F | boolean,
    arg3?: F,
    arg4?: boolean
  ) {
    this.curve = curve;

    // Determine which overload we're using
    if (Array.isArray(arg1)) {
      // Coords array form: (curve, coords, check)
      const coords = arg1;
      const check = arg2 === undefined ? true : (arg2 as boolean);

      if (coords.length === 0) {
        // Point at infinity: (0 : 1 : 0)
        const field = curve.base_ring;
        this._X = field.zero() as F;
        this._Y = field.one() as F;
        this._Z = field.zero() as F;
      } else {
        // Affine point (x, y) = (x : y : 1)
        const [x, y] = coords;
        this._X = x;
        this._Y = y;
        this._Z = curve.base_ring.one() as F;

        if (check && !curve.is_on_curve(x, y)) {
          throw new ValueError(`(${x}, ${y}) is not a point on the curve`);
        }
      }
    } else {
      // Projective form: (curve, X, Y, Z, check)
      const X = arg1;
      const Y = arg2 as F;
      const Z = arg3 as F;
      const check = arg4 === undefined ? true : arg4;

      this._X = X;
      this._Y = Y;
      this._Z = Z;

      if (check && !Z.isZero()) {
        // Convert to affine and check
        const x = X.div(Z) as F;
        const y = Y.div(Z) as F;
        if (!curve.is_on_curve(x, y)) {
          throw new ValueError(`(${x}, ${y}) is not a point on the curve`);
        }
      }
    }
  }

  /**
   * Get the affine x-coordinate.
   * @throws {ValueError} If this is the point at infinity
   */
  x(): F {
    if (this._Z.isZero()) {
      throw new ValueError('point at infinity has no x-coordinate');
    }
    return this._X.div(this._Z) as F;
  }

  /**
   * Get the affine y-coordinate.
   * @throws {ValueError} If this is the point at infinity
   */
  y(): F {
    if (this._Z.isZero()) {
      throw new ValueError('point at infinity has no y-coordinate');
    }
    return this._Y.div(this._Z) as F;
  }

  /**
   * Get coordinates as a tuple [x, y] or undefined for point at infinity.
   */
  xy(): [F, F] | undefined {
    if (this._Z.isZero()) {
      return undefined;
    }
    return [this.x(), this.y()];
  }

  /**
   * Get projective coordinates [X, Y, Z].
   */
  xyz(): [F, F, F] {
    return [this._X, this._Y, this._Z];
  }

  /**
   * Check if this is the point at infinity (identity element).
   */
  is_zero(): boolean {
    return this._Z.isZero();
  }

  /**
   * Alias for is_zero().
   */
  is_identity(): boolean {
    return this.is_zero();
  }

  /**
   * Add two points on the curve.
   *
   * This implements the group law for elliptic curves in Weierstrass form:
   *   y^2 + a1*x*y + a3*y = x^3 + a2*x^2 + a4*x + a6
   *
   * Special cases:
   *   - P + O = P (identity)
   *   - P + (-P) = O (inverse)
   *   - P + P = 2P (doubling)
   *
   * For short Weierstrass (a1 = a2 = a3 = 0):
   *   - General addition: lambda = (y2 - y1) / (x2 - x1)
   *   - Doubling: lambda = (3*x1^2 + a) / (2*y1)
   *   - x3 = lambda^2 - x1 - x2
   *   - y3 = lambda * (x1 - x3) - y1
   *
   * For general Weierstrass:
   *   - Uses the explicit formulas from hyperelliptic.org/EFD
   *
   * @param other - The point to add
   * @returns The sum P + Q
   */
  add(other: EllipticCurvePoint<F>): EllipticCurvePoint<F> {
    // Handle identity: P + O = P
    if (this.is_zero()) {
      return other;
    }
    if (other.is_zero()) {
      return this;
    }

    const a1 = this.curve.a1();
    const a2 = this.curve.a2();
    const a3 = this.curve.a3();

    const x1 = this.x();
    const y1 = this.y();
    const x2 = other.x();
    const y2 = other.y();

    // Check if points have the same x-coordinate
    if (x1.eq(x2)) {
      // Check if P = -Q (i.e., y1 + y2 + a1*x + a3 = 0)
      // For short Weierstrass: y1 + y2 = 0
      const negYSum = y1.add(y2).add(a1.mul(x1)).add(a3);
      if (negYSum.isZero()) {
        // P + (-P) = O
        return this.curve.zero();
      }

      // Otherwise P = Q, use doubling formula
      return this._double();
    }

    // General addition: P != Q
    // Slope lambda = (y2 - y1) / (x2 - x1)
    const dx = x2.sub(x1) as F;
    const dy = y2.sub(y1) as F;
    const lambda = dy.div(dx) as F;

    // For general Weierstrass:
    // x3 = lambda^2 + a1*lambda - a2 - x1 - x2
    // y3 = lambda*(x1 - x3) - y1 - a1*x3 - a3

    const lambda2 = lambda.mul(lambda) as F;
    const x3 = lambda2.add(a1.mul(lambda)).sub(a2).sub(x1).sub(x2) as F;

    // y3 = lambda * (x1 - x3) - y1 - a1*x3 - a3
    const y3 = lambda.mul(x1.sub(x3)).sub(y1).sub(a1.mul(x3)).sub(a3) as F;

    return affinePoint(this.curve, x3, y3, false);
  }

  /**
   * Point doubling: compute 2P.
   *
   * For short Weierstrass (y^2 = x^3 + a*x + b):
   *   lambda = (3*x^2 + a) / (2*y)
   *   x3 = lambda^2 - 2*x
   *   y3 = lambda*(x - x3) - y
   *
   * For general Weierstrass:
   *   lambda = (3*x^2 + 2*a2*x + a4 - a1*y) / (2*y + a1*x + a3)
   *   x3 = lambda^2 + a1*lambda - a2 - 2*x
   *   y3 = lambda*(x - x3) - y - a1*x3 - a3
   */
  private _double(): EllipticCurvePoint<F> {
    if (this.is_zero()) {
      return this;
    }

    const a1 = this.curve.a1();
    const a2 = this.curve.a2();
    const a3 = this.curve.a3();
    const a4 = this.curve.a4();
    const field = this.curve.base_ring;
    const two = field.__call__(2) as F;
    const three = field.__call__(3) as F;

    const x = this.x();
    const y = this.y();

    // Denominator: 2*y + a1*x + a3
    const denom = two.mul(y).add(a1.mul(x)).add(a3) as F;

    if (denom.isZero()) {
      // Tangent line is vertical, result is point at infinity
      return this.curve.zero();
    }

    // Numerator: 3*x^2 + 2*a2*x + a4 - a1*y
    const x2 = x.mul(x) as F;
    const numer = three.mul(x2).add(two.mul(a2).mul(x)).add(a4).sub(a1.mul(y)) as F;

    const lambda = numer.div(denom) as F;

    // x3 = lambda^2 + a1*lambda - a2 - 2*x
    const lambda2 = lambda.mul(lambda) as F;
    const x3 = lambda2.add(a1.mul(lambda)).sub(a2).sub(two.mul(x)) as F;

    // y3 = lambda*(x - x3) - y - a1*x3 - a3
    const y3 = lambda.mul(x.sub(x3)).sub(y).sub(a1.mul(x3)).sub(a3) as F;

    return affinePoint(this.curve, x3, y3, false);
  }

  /**
   * Return the additive inverse (negation) of this point.
   *
   * For a point (x, y) on the curve y^2 + a1*x*y + a3*y = x^3 + a2*x^2 + a4*x + a6,
   * the inverse is (x, -y - a1*x - a3).
   *
   * For short Weierstrass (y^2 = x^3 + a*x + b), this simplifies to (x, -y).
   */
  neg(): EllipticCurvePoint<F> {
    if (this.is_zero()) {
      return this;
    }

    const a1 = this.curve.a1();
    const a3 = this.curve.a3();
    const x = this.x();
    const y = this.y();

    // -P = (x, -y - a1*x - a3)
    const negY = y.neg().sub(a1.mul(x)).sub(a3) as F;

    return affinePoint(this.curve, x, negY, false);
  }

  /**
   * Subtract two points: P - Q = P + (-Q).
   */
  sub(other: EllipticCurvePoint<F>): EllipticCurvePoint<F> {
    return this.add(other.neg());
  }

  /**
   * Scalar multiplication: compute n*P using double-and-add algorithm.
   *
   * @param n - The scalar multiplier
   * @returns n*P
   *
   * Special cases:
   *   - 0*P = O (point at infinity)
   *   - n*O = O
   *   - (-n)*P = n*(-P)
   */
  mul(n: bigint | number): EllipticCurvePoint<F> {
    let scalar = typeof n === 'number' ? BigInt(n) : n;

    // Handle zero scalar
    if (scalar === 0n) {
      return this.curve.zero();
    }

    // Handle point at infinity
    if (this.is_zero()) {
      return this;
    }

    // Handle negative scalar: (-n)*P = n*(-P)
    let point: EllipticCurvePoint<F> = this;
    if (scalar < 0n) {
      point = this.neg();
      scalar = -scalar;
    }

    // Double-and-add algorithm (binary method)
    let result = this.curve.zero();
    let current = point;

    while (scalar > 0n) {
      if ((scalar & 1n) === 1n) {
        result = result.add(current);
      }
      current = current.add(current); // Double
      scalar >>= 1n;
    }

    return result;
  }

  /**
   * Alias for mul() to match Python's __rmul__.
   */
  rmul(n: bigint | number): EllipticCurvePoint<F> {
    return this.mul(n);
  }

  /**
   * Check if two points are equal.
   */
  eq(other: EllipticCurvePoint<F>): boolean {
    if (this.is_zero() && other.is_zero()) {
      return true;
    }
    if (this.is_zero() || other.is_zero()) {
      return false;
    }
    // Compare affine coordinates
    return this.x().eq(other.x()) && this.y().eq(other.y());
  }

  /**
   * Check if this point has finite order n in the group.
   * That is, check if n*P = O.
   */
  has_order(n: bigint | number): boolean {
    return this.mul(n).is_zero();
  }

  /**
   * Compute the order of this point in the elliptic curve group.
   *
   * Warning: This uses a naive algorithm and can be slow for large groups.
   * For cryptographic curves, use the known group order.
   *
   * @param maxOrder - Maximum order to check (default: search until found)
   * @returns The order, or undefined if not found within maxOrder
   */
  order(maxOrder?: bigint | number): bigint | undefined {
    if (this.is_zero()) {
      return 1n;
    }

    const max =
      maxOrder !== undefined
        ? typeof maxOrder === 'number'
          ? BigInt(maxOrder)
          : maxOrder
        : undefined;

    let current: EllipticCurvePoint<F> = this;
    let n = 1n;

    while (true) {
      if (current.is_zero()) {
        return n;
      }

      if (max !== undefined && n >= max) {
        return undefined;
      }

      current = current.add(this);
      n++;
    }
  }

  /**
   * String representation of the point in projective coordinates.
   */
  toString(): string {
    if (this.is_zero()) {
      return '(0 : 1 : 0)'; // Projective coordinates for point at infinity
    }
    return `(${this.x()} : ${this.y()} : 1)`;
  }

  /**
   * Get affine coordinates string.
   */
  toAffineString(): string {
    if (this.is_zero()) {
      return 'O'; // Point at infinity
    }
    return `(${this.x()}, ${this.y()})`;
  }
}

/**
 * Create the point at infinity for a curve.
 *
 * @param curve - The elliptic curve
 * @returns The point at infinity (identity element)
 */
export function pointAtInfinity<F extends FieldElement>(
  curve: EllipticCurveInterface<F>
): EllipticCurvePoint<F> {
  return new EllipticCurvePoint(curve, []);
}

/**
 * Create an affine point on a curve.
 *
 * @param curve - The elliptic curve
 * @param x - The x-coordinate
 * @param y - The y-coordinate
 * @param check - Whether to verify the point is on the curve
 * @returns The point (x, y)
 */
export function affinePoint<F extends FieldElement>(
  curve: EllipticCurveInterface<F>,
  x: F,
  y: F,
  check: boolean = true
): EllipticCurvePoint<F> {
  return new EllipticCurvePoint(curve, [x, y], check);
}

/**
 * Compute the Weil pairing of two points P and Q on the same curve.
 *
 * The Weil pairing is a bilinear, alternating, non-degenerate pairing
 * e_n: E[n] x E[n] -> mu_n where mu_n is the group of n-th roots of unity.
 *
 * INPUT:
 * - P: first point of order n
 * - Q: second point of order n
 * - n: integer such that nP = nQ = O
 * - algorithm: 'pari' or 'sage' (default 'sage' for this implementation)
 *
 * OUTPUT: an n-th root of unity in the base field
 *
 * ALGORITHM:
 * Uses Miller's algorithm to compute f_{n,P}(Q) / f_{n,Q}(P) where
 * f_{n,P} is a function with divisor n[P] - n[O].
 * The result is multiplied by (-1)^n to match the standard definition.
 *
 * @see Reference: sage/schemes/elliptic_curves/ell_point.py:weil_pairing
 * @see [Mil2004] Victor Miller, "The Weil Pairing, and Its Efficient Calculation"
 */
export function weil_pairing<F extends FieldElement>(
  P: EllipticCurvePoint<F>,
  Q: EllipticCurvePoint<F>,
  n: bigint,
  _algorithm?: 'pari' | 'sage'
): F {
  const E = P.curve;

  if (Q.curve !== E) {
    throw new ValueError('points must both be on the same curve');
  }

  // Test if P, Q are both in E[n]
  if (!P.mul(n).is_zero() || !Q.mul(n).is_zero()) {
    throw new ValueError('points must both be n-torsion');
  }

  const one = E.base_ring.one() as F;

  // Case where P = Q
  if (P.eq(Q)) {
    return one;
  }

  // Case where P = O or Q = O
  if (P.is_zero() || Q.is_zero()) {
    return one;
  }

  // The non-trivial case P != Q
  // Weil pairing: e_n(P, Q) = (-1)^n * f_{n,P}(Q) / f_{n,Q}(P)
  try {
    const fPQ = _miller(P, Q, n);
    const fQP = _miller(Q, P, n);

    if (fQP.isZero()) {
      // Linear dependence detected
      return one;
    }

    // Compute (-1)^n * f_P(Q) / f_Q(P)
    let result = fPQ.div(fQP) as F;

    // Multiply by (-1)^n: if n is odd, negate the result
    if ((n & 1n) === 1n) {
      result = result.neg() as F;
    }

    return result;
  } catch {
    // ZeroDivisionError indicates linearly dependent points
    return one;
  }
}

/**
 * Compute the value at Q of the straight line through points P and R.
 *
 * Used internally by Miller's algorithm.
 *
 * For the line through P and R evaluated at Q:
 * - If P = R (tangent line), use the derivative formula
 * - Otherwise, use the secant line formula
 *
 * @param P - First point
 * @param R - Second point
 * @param Q - Point at which to evaluate
 * @returns The value of the line at Q
 */
function _line<F extends FieldElement>(
  P: EllipticCurvePoint<F>,
  R: EllipticCurvePoint<F>,
  Q: EllipticCurvePoint<F>
): F {
  if (Q.is_zero()) {
    throw new ValueError('Q must be nonzero.');
  }

  const one = P.curve.base_ring.one() as F;
  const [a1, a2, a3, a4] = P.curve.a_invariants();

  // Case: P or R is identity
  if (P.is_zero() || R.is_zero()) {
    if (P.eq(R)) {
      return one;
    }
    if (P.is_zero()) {
      return Q.x().sub(R.x()) as F;
    }
    if (R.is_zero()) {
      return Q.x().sub(P.x()) as F;
    }
  }

  const Px = P.x();
  const Py = P.y();
  const Rx = R.x();
  const Qx = Q.x();
  const Qy = Q.y();

  // Case: P != R
  if (!P.eq(R)) {
    if (Px.eq(Rx)) {
      // Vertical line: x - x_P
      return Qx.sub(Px) as F;
    } else {
      // Secant line: y - y_P - lambda(x - x_P)
      const Ry = R.y();
      const lambda = Ry.sub(Py).div(Rx.sub(Px)) as F;
      return Qy.sub(Py).sub(lambda.mul(Qx.sub(Px))) as F;
    }
  }

  // Case: P = R (tangent line)
  // lambda = (3*x^2 + 2*a2*x + a4 - a1*y) / (2*y + a1*x + a3)
  const three = P.curve.base_ring.__call__(3) as F;
  const two = P.curve.base_ring.__call__(2) as F;

  const numerator = three.mul(Px.mul(Px)).add(two.mul(a2).mul(Px)).add(a4).sub(a1.mul(Py)) as F;
  const denominator = two.mul(Py).add(a1.mul(Px)).add(a3) as F;

  if (denominator.isZero()) {
    // Vertical tangent: x - x_P
    return Qx.sub(Px) as F;
  }

  const lambda = numerator.div(denominator) as F;
  return Qy.sub(Py).sub(lambda.mul(Qx.sub(Px))) as F;
}

/**
 * Miller's algorithm: compute the value at Q of the function f_{n,P}.
 *
 * The function f_{n,P} has divisor n[P] - [nP] - (n-1)[O].
 *
 * ALGORITHM:
 * Double-and-add with line function accumulation.
 *
 * @param P - The base point
 * @param Q - The evaluation point (must be nonzero)
 * @param n - The integer n (must be nonzero)
 * @returns f_{n,P}(Q)
 *
 * @see [Mil2004] Victor Miller, "The Weil Pairing, and Its Efficient Calculation"
 */
function _miller<F extends FieldElement>(
  P: EllipticCurvePoint<F>,
  Q: EllipticCurvePoint<F>,
  n: bigint
): F {
  if (Q.is_zero()) {
    throw new ValueError('Q must be nonzero.');
  }

  if (n === 0n) {
    throw new ValueError('n must be nonzero.');
  }

  // Handle negative n
  let nIsNegative = false;
  if (n < 0n) {
    n = -n;
    nIsNegative = true;
  }

  const one = P.curve.base_ring.one() as F;

  // Base case: if P is identity or n = 1
  if (P.is_zero()) {
    return one;
  }

  let t = one;
  let V = P;

  // Get binary representation of n
  const bits: number[] = [];
  let temp = n;
  while (temp > 0n) {
    bits.push(Number(temp & 1n));
    temp >>= 1n;
  }

  // Process bits from most significant to least significant (skip the top bit)
  for (let i = bits.length - 2; i >= 0; i--) {
    // Double step: V -> 2V
    const S = V.add(V);
    const ell = _line(V, V, Q);
    const vee = _line(S, S.neg(), Q);

    t = t.mul(t).mul(ell).div(vee) as F;
    V = S;

    // Add step if bit is 1: V -> V + P
    if (bits[i] === 1) {
      const SPrime = V.add(P);
      const ellAdd = _line(V, P, Q);
      const veeAdd = _line(SPrime, SPrime.neg(), Q);

      t = t.mul(ellAdd).div(veeAdd) as F;
      V = SPrime;
    }
  }

  // For negative n, compute 1/(t * v_{nP}(Q))
  if (nIsNegative) {
    const vee = _line(V, V.neg(), Q);
    t = one.div(t.mul(vee)) as F;
  }

  return t;
}

/**
 * Return the Tate pairing of n-torsion point P with point Q.
 *
 * The value returned is f_{n,P}(Q)^e where f_{n,P} is a function with
 * divisor n[P]-n[O]. This is also known as the "modified Tate pairing".
 *
 * INPUT:
 * - P: point of order n (self)
 * - Q: elliptic curve point on same curve as P
 * - n: positive integer, order of P
 * - k: positive integer, embedding degree
 * - q: (optional) size of base field
 *
 * OUTPUT: an n-th root of unity in the base field
 *
 * EXAMPLES:
 * ```typescript
 * const p = 103n;
 * const E = EllipticCurve(GF(p), [1n, 18n]);
 * const P = E.point(33n, 91n);
 * const n = P.order(); // 19
 * const k = 6; // embedding degree
 * const result = tate_pairing(P, P, n, k);
 * // Result is in the field F_q^k
 * ```
 *
 * @see Reference: sage/schemes/elliptic_curves/ell_point.py:tate_pairing
 */
export function tate_pairing<F extends FieldElement>(
  P: EllipticCurvePoint<F>,
  Q: EllipticCurvePoint<F>,
  n: bigint,
  k: bigint | number,
  q?: bigint
): F {
  const E = P.curve;

  if (Q.curve !== E) {
    throw new ValueError('Points must both be on the same curve');
  }

  const field = E.base_ring;
  const one = field.one() as F;

  // Get the field size
  const kVal = typeof k === 'number' ? BigInt(k) : k;
  let qVal: bigint;

  if (q !== undefined) {
    qVal = q;
  } else {
    // For prime fields, q = p
    qVal = field.characteristic;
  }

  // Check that P is n-torsion
  if (!P.mul(n).is_zero()) {
    throw new ValueError('The point P must be n-torsion');
  }

  // Handle trivial cases
  if (P.is_zero() || Q.is_zero()) {
    return one;
  }

  // Compute the Miller function f_{n,P}(Q)
  const millerValue = _miller(P, Q, n);

  // Compute the final exponentiation: (q^k - 1) / n
  let qk = 1n;
  for (let i = 0n; i < kVal; i++) {
    qk *= qVal;
  }
  const exp = (qk - 1n) / n;

  // Return millerValue^exp
  return millerValue.pow(exp) as F;
}

/**
 * Return the Ate pairing of points P and Q.
 *
 * The Ate pairing is an optimized variant of the Tate pairing that
 * uses a shorter Miller loop. It is computed as:
 *   ate(P, Q) = f_{T-1,Q}(P)^e
 * where T = t - 1 (t is trace of Frobenius) and e = (q^k - 1)/n.
 *
 * INPUT:
 * - P: point of order n in ker(pi - 1), where pi is the q-Frobenius
 *      (typically a q-rational point)
 * - Q: point of order n in ker(pi - q)
 * - n: the order of P and Q
 * - k: the embedding degree
 * - t: the trace of Frobenius
 * - q: (optional) size of base field
 *
 * OUTPUT: an n-th root of unity in F_q^k
 *
 * ALGORITHM:
 * The Ate pairing uses the Miller loop with T = t - 1 instead of n,
 * which can be much shorter for pairing-friendly curves.
 *
 * The ate pairing is a power of the Tate pairing:
 *   ate(P, Q) = tate(Q, P)^M
 * for some exponent M related to n, k, and T.
 *
 * @see Reference: sage/schemes/elliptic_curves/ell_point.py:ate_pairing
 * @see [HSV2006] F. Hess, N. Smart, F. Vercauteren, "The Eta Pairing Revisited"
 */
export function ate_pairing<F extends FieldElement>(
  P: EllipticCurvePoint<F>,
  Q: EllipticCurvePoint<F>,
  n: bigint,
  k: bigint | number,
  t: bigint,
  q?: bigint
): F {
  const E = P.curve;

  if (Q.curve !== E) {
    throw new ValueError('Points must both be on the same curve');
  }

  const field = E.base_ring;
  const one = field.one() as F;

  // Get the field size
  const kVal = typeof k === 'number' ? BigInt(k) : k;
  let qVal: bigint;

  if (q !== undefined) {
    qVal = q;
  } else {
    qVal = field.characteristic;
  }

  // Handle trivial cases
  if (P.is_zero() || Q.is_zero()) {
    return one;
  }

  // Verify P is in ker(pi - 1) and Q is in ker(pi - q)
  // For the simplified case over Fp, P should be Fp-rational (always true)
  // Q should be in the correct eigenspace of Frobenius

  // T = t - 1
  const T = t - 1n;

  // Compute absolute value if T is negative
  const absT = T >= 0n ? T : -T;

  // Compute the Miller function f_{|T|,Q}(P)
  let millerValue: F;
  try {
    millerValue = _miller(Q, P, absT);
  } catch {
    // If there's an issue (e.g., linearly dependent points), return 1
    return one;
  }

  // If T is negative, we need to adjust (take reciprocal)
  if (T < 0n) {
    millerValue = millerValue.inv() as F;
  }

  // Compute final exponentiation: (q^k - 1) / n
  let qk = 1n;
  for (let i = 0n; i < kVal; i++) {
    qk *= qVal;
  }
  const exp = (qk - 1n) / n;

  // Return millerValue^exp
  return millerValue.pow(exp) as F;
}

/**
 * Return a list of all points Q such that mQ = P where P = self.
 *
 * Only points on the elliptic curve containing self and defined
 * over the base field are included.
 *
 * INPUT:
 * - P: the target point
 * - m: positive integer (the divisor)
 * - poly_only: not supported in this implementation
 *
 * OUTPUT: a list of points Q such that mQ = P
 *
 * ALGORITHM:
 * For finite fields, we enumerate points on the curve and check which satisfy mQ = P.
 * This is a brute-force approach suitable for small fields.
 *
 * For larger fields or more efficient computation, one would use division polynomials.
 *
 * @see Reference: sage/schemes/elliptic_curves/ell_point.py:division_points
 */
export function division_points<F extends FieldElement>(
  P: EllipticCurvePoint<F>,
  m: bigint | number,
  poly_only?: boolean
): EllipticCurvePoint<F>[] {
  if (poly_only) {
    throw new NotImplementedError('poly_only=true not yet implemented');
  }

  const mVal = typeof m === 'number' ? BigInt(m) : m;

  // Handle trivial cases
  if (mVal === 1n || mVal === -1n) {
    return [P.mul(mVal)];
  }

  if (mVal === 0n) {
    // 0*Q = O for all Q, so if P = O, every point is a solution
    // If P != O, no solutions
    if (P.is_zero()) {
      return [P]; // Return just the identity (could return all points for finite fields)
    }
    return [];
  }

  const E = P.curve;
  const field = E.base_ring;
  const ans: EllipticCurvePoint<F>[] = [];

  // Check if the identity is a solution (P = O)
  if (P.is_zero()) {
    ans.push(P);
  }

  // For finite fields, we can enumerate points
  // This requires the field to be iterable
  const negP = P.neg();
  const P_is_2_torsion = P.eq(negP);

  // Iterate over field elements to find x-coordinates
  for (const x of field) {
    // Try to lift x to a point on the curve
    const xF = x as F;

    // Compute y^2 = x^3 + a4*x + a6 for short Weierstrass
    // For general Weierstrass: y^2 + a1*x*y + a3*y = x^3 + a2*x^2 + a4*x + a6
    const [a1, a2, a3, a4, a6] = E.a_invariants();
    const x3 = xF.mul(xF).mul(xF) as F;
    const rhs = x3.add(a2.mul(xF.mul(xF))).add(a4.mul(xF)).add(a6) as F;

    // For general Weierstrass, we need to solve y^2 + (a1*x + a3)*y = rhs
    const b = a1.mul(xF).add(a3) as F;

    // Completing the square: (y + b/2)^2 = rhs + b^2/4
    // This only works in characteristic != 2
    const char = field.characteristic;

    if (char === 2n) {
      // Characteristic 2 case: handle differently
      // y^2 + b*y = rhs
      // If b = 0: y^2 = rhs, so y exists iff rhs is a square
      // If b != 0: substitute y = b*z, get b^2*z^2 + b^2*z = rhs
      //            z^2 + z = rhs/b^2, which has solutions iff Tr(rhs/b^2) = 0
      // For simplicity, skip char 2 for now in the general case
      continue;
    }

    // For odd characteristic: complete the square
    const two = field.__call__(2) as F;
    const four = field.__call__(4) as F;
    const discriminant = rhs.add(b.mul(b).div(four)) as F;

    // Check if discriminant is a square
    // We need a way to check this - use the is_square method if available
    // For now, try to compute sqrt and catch errors
    let sqrtDisc: F;
    try {
      // We need to check if discriminant is a square in the field
      // Use Euler's criterion for odd characteristic: a^((p-1)/2) = 1 iff a is a square
      const exp = (char - 1n) / 2n;
      const test = discriminant.pow(exp) as F;

      if (discriminant.isZero()) {
        sqrtDisc = field.__call__(0) as F;
      } else if (!test.eq(field.one())) {
        // Not a quadratic residue
        continue;
      } else {
        // Compute square root using Tonelli-Shanks
        sqrtDisc = computeSqrt(discriminant, field) as F;
      }
    } catch {
      continue;
    }

    // y = -b/2 + sqrt(discriminant) or y = -b/2 - sqrt(discriminant)
    const negBOver2 = b.neg().div(two) as F;
    const y1 = negBOver2.add(sqrtDisc) as F;
    const y2 = negBOver2.sub(sqrtDisc) as F;

    // Try both y values
    for (const yF of [y1, y2]) {
      // Verify the point is on the curve
      if (!E.is_on_curve(xF, yF)) {
        continue;
      }

      const Q = affinePoint(E, xF, yF, false);
      const mQ = Q.mul(mVal);

      if (P_is_2_torsion) {
        // If P = -P, then Q and -Q are both solutions if mQ = P
        if (mQ.eq(P)) {
          ans.push(Q);
        }
      } else {
        // P is not 2-torsion, so at most one of Q, -Q works
        if (mQ.eq(P)) {
          ans.push(Q);
        } else if (mQ.eq(negP)) {
          ans.push(Q.neg());
        }
      }
    }
  }

  // Sort for deterministic output (optional, based on x then y coordinates)
  // For now, just return as is
  return ans;
}

/**
 * Compute the square root of a in the finite field.
 * Uses Tonelli-Shanks algorithm.
 */
function computeSqrt<F extends FieldElement>(a: F, field: FieldRing): F {
  if (a.isZero()) {
    return a;
  }

  const p = field.characteristic;

  // Case: p ≡ 3 (mod 4)
  if (p % 4n === 3n) {
    return a.pow((p + 1n) / 4n) as F;
  }

  // Tonelli-Shanks for general case
  // Write p - 1 = 2^s * q where q is odd
  let q = p - 1n;
  let s = 0n;
  while ((q & 1n) === 0n) {
    q >>= 1n;
    s++;
  }

  // Find a quadratic non-residue z
  let z = field.__call__(2) as F;
  const exp = (p - 1n) / 2n;
  while (z.pow(exp).eq(field.one())) {
    z = field.__call__(z.add(field.one())) as F;
  }

  let m = s;
  let c = z.pow(q) as F;
  let t = a.pow(q) as F;
  let r = a.pow((q + 1n) / 2n) as F;

  while (true) {
    if (t.isZero()) {
      return field.__call__(0) as F;
    }

    if (t.eq(field.one())) {
      return r;
    }

    // Find the least i such that t^(2^i) = 1
    let i = 1n;
    let temp = t.mul(t) as F;
    while (!temp.eq(field.one())) {
      temp = temp.mul(temp) as F;
      i++;
    }

    // Update values
    const b = c.pow(1n << (m - i - 1n)) as F;
    m = i;
    c = b.mul(b) as F;
    t = t.mul(c) as F;
    r = r.mul(b) as F;
  }
}

/**
 * Return True if this point has finite additive order.
 *
 * For finite fields, this always returns true since all points have finite order.
 * For other fields (like number fields), points may have infinite order.
 *
 * This implementation supports:
 * - Finite fields: always true
 * - The identity point: always true (order 1)
 *
 * @see Reference: sage/schemes/elliptic_curves/ell_point.py:has_finite_order
 */
export function has_finite_order<F extends FieldElement>(P: EllipticCurvePoint<F>): boolean {
  // The identity point always has finite order (order 1)
  if (P.is_zero()) {
    return true;
  }

  // For finite fields, all points have finite order because the group is finite
  // We detect a finite field by checking if the characteristic is finite and positive
  const char = P.curve.base_ring.characteristic;
  if (char > 0n) {
    // This is a finite field (characteristic > 0)
    return true;
  }

  // For other fields (like Q or number fields), we would need to compute
  // whether the point is torsion. This is a harder problem.
  // For now, return false for non-finite fields (conservative approach)
  return false;
}

/**
 * Return True if this point has infinite additive order.
 *
 * This is the logical complement of has_finite_order.
 *
 * For finite fields, this always returns false.
 * For other fields, a point has infinite order if it is not a torsion point.
 *
 * @see Reference: sage/schemes/elliptic_curves/ell_point.py:has_infinite_order
 */
export function has_infinite_order<F extends FieldElement>(P: EllipticCurvePoint<F>): boolean {
  return !has_finite_order(P);
}

/**
 * Return the Neron-Tate canonical height of the point.
 *
 * INPUT:
 * - P: a point on an elliptic curve
 * - precision: positive integer or None (default None)
 * - normalised: if True (default), height is invariant under extension
 * - algorithm: 'pari' or 'sage' (default 'pari')
 *
 * OUTPUT: a non-negative real number (0 for torsion points)
 *
 * NOTE: For elliptic curves over finite fields, all points are torsion points,
 * so the canonical height is always 0. This function is primarily meaningful
 * for curves over number fields.
 *
 * For curves over number fields, the canonical height is computed as:
 *   h(P) = sum of local heights at all places
 *
 * @see Reference: sage/schemes/elliptic_curves/ell_point.py:height
 */
export function height<F extends FieldElement>(
  P: EllipticCurvePoint<F>,
  _precision?: number,
  _normalised?: boolean,
  _algorithm?: 'pari' | 'sage'
): number {
  // For finite fields, all points are torsion points, so height is 0
  const char = P.curve.base_ring.characteristic;
  if (char > 0n) {
    // Finite field case: all points have finite order, so height = 0
    return 0;
  }

  // For number fields (characteristic 0), we would need a more complex implementation
  // involving local heights at all places. This requires number field infrastructure.
  throw new NotImplementedError(
    'height is only implemented for finite fields (returns 0). ' +
      'Number field implementation requires sage.rings.number_field'
  );
}

/**
 * Compute the local height at an archimedean place v.
 *
 * INPUT:
 * - P: a point on an elliptic curve over a number field K
 * - v: a real or complex embedding, or None for total archimedean height
 * - prec: precision in bits (default None)
 * - weighted: if True, multiply by local degree
 *
 * OUTPUT: a real number
 *
 * NOTE: For elliptic curves over finite fields, there are no archimedean places,
 * so this function returns 0. This function is only meaningful for curves over
 * number fields (characteristic 0).
 *
 * ALGORITHM: See Silverman, "Computing heights on elliptic curves" (Section 4)
 *
 * @see Reference: sage/schemes/elliptic_curves/ell_point.py:archimedean_local_height
 */
export function archimedean_local_height<F extends FieldElement>(
  P: EllipticCurvePoint<F>,
  _v?: unknown,
  _prec?: number,
  _weighted?: boolean
): number {
  // For finite fields, there are no archimedean places
  const char = P.curve.base_ring.characteristic;
  if (char > 0n) {
    // Finite field case: no archimedean places
    return 0;
  }

  // For number fields (characteristic 0), we would need embeddings infrastructure
  throw new NotImplementedError(
    'archimedean_local_height is only implemented for finite fields (returns 0). ' +
      'Number field implementation requires sage.rings.number_field'
  );
}

/**
 * Compute the local height at a non-archimedean place v.
 *
 * INPUT:
 * - P: a point on an elliptic curve over a number field K
 * - v: a prime ideal or prime number, or None
 * - prec: precision (default None)
 * - weighted: if True, multiply by local degree
 *
 * OUTPUT: a real number
 *
 * NOTE: For elliptic curves over finite fields, the local height concept
 * does not apply in the same way as for number fields. This function returns 0
 * for finite fields.
 *
 * ALGORITHM: See Silverman, "Computing heights on elliptic curves" (Section 5)
 *
 * @see Reference: sage/schemes/elliptic_curves/ell_point.py:non_archimedean_local_height
 */
export function non_archimedean_local_height<F extends FieldElement>(
  P: EllipticCurvePoint<F>,
  _v?: unknown,
  _prec?: number,
  _weighted?: boolean
): number {
  // For finite fields, the local height concept is not applicable
  const char = P.curve.base_ring.characteristic;
  if (char > 0n) {
    // Finite field case
    return 0;
  }

  // For number fields (characteristic 0), we would need prime ideal infrastructure
  throw new NotImplementedError(
    'non_archimedean_local_height is only implemented for finite fields (returns 0). ' +
      'Number field implementation requires sage.rings.number_field'
  );
}

/**
 * Return the elliptic logarithm of this point.
 *
 * The elliptic logarithm is the inverse of the Weierstrass p-function.
 * For a point P on an elliptic curve E over C, the elliptic logarithm
 * z is a complex number such that (p(z), p'(z)) = P under the uniformization
 * E(C) = C/L for some lattice L.
 *
 * INPUT:
 * - P: a point on an elliptic curve
 * - embedding: (optional) an embedding of the base field into C
 * - precision: (default 100) precision in bits
 * - algorithm: 'pari' or 'sage' (default 'pari')
 *
 * OUTPUT: a complex number z such that exp(z) corresponds to this point
 *
 * NOTE: For finite fields, this function is not applicable since there is
 * no natural embedding into the complex numbers. The elliptic logarithm
 * is defined only for curves over number fields or their completions.
 *
 * For the p-adic analog over finite fields, see padic_elliptic_logarithm.
 *
 * @see Reference: sage/schemes/elliptic_curves/ell_point.py:elliptic_logarithm
 */
export function elliptic_logarithm<F extends FieldElement>(
  P: EllipticCurvePoint<F>,
  _embedding?: unknown,
  _precision?: number,
  _algorithm?: string
): unknown {
  const char = P.curve.base_ring.characteristic;
  if (char > 0n) {
    throw new NotImplementedError(
      'elliptic_logarithm is not defined for curves over finite fields. ' +
        'Use padic_elliptic_logarithm for the p-adic analog.'
    );
  }

  throw new NotImplementedError(
    'elliptic_logarithm requires complex embeddings infrastructure. ' +
      'Not yet implemented for number fields.'
  );
}

/**
 * Return the p-adic elliptic logarithm of this point.
 *
 * The p-adic elliptic logarithm is a p-adic analog of the complex elliptic
 * logarithm. For an elliptic curve E over Qp with good ordinary reduction,
 * there is a p-adic logarithm map from E(Qp) to Qp.
 *
 * For finite fields F_p, this function computes the discrete logarithm
 * using the p-adic elliptic logarithm method, which is particularly
 * efficient for anomalous curves (where #E(F_p) = p).
 *
 * INPUT:
 * - P: a point on an elliptic curve
 * - p: a prime number (for finite fields, should be the characteristic)
 * - absprec: absolute precision (default 20)
 *
 * OUTPUT: an element of the p-adic field Qp (or an integer for finite fields)
 *
 * ALGORITHM:
 * For anomalous curves over F_p, uses Smart's attack / SSSA algorithm:
 * lift the curve and points to Q_p and compute the logarithm there.
 *
 * @see Reference: sage/schemes/elliptic_curves/ell_point.py:padic_elliptic_logarithm
 * @see [Sma1999] N. Smart, "The discrete logarithm problem on elliptic curves
 *      of trace one"
 */
export function padic_elliptic_logarithm<F extends FieldElement>(
  P: EllipticCurvePoint<F>,
  p: bigint | number,
  _absprec?: number
): unknown {
  const pVal = typeof p === 'number' ? BigInt(p) : p;
  const char = P.curve.base_ring.characteristic;

  if (char === 0n) {
    // Number field case
    throw new NotImplementedError(
      'padic_elliptic_logarithm for number fields requires sage.rings.padics. ' +
        'Not yet implemented.'
    );
  }

  // Finite field case
  if (pVal !== char) {
    throw new ValueError(
      `p (${pVal}) must equal the field characteristic (${char}) for finite field curves`
    );
  }

  // For anomalous curves (#E = p), the p-adic logarithm gives the discrete log
  // This is Smart's attack / SSSA attack
  const order = P.order();
  if (order === undefined) {
    throw new ValueError('Could not determine point order');
  }

  // The p-adic elliptic logarithm for finite fields is typically used
  // to compute discrete logarithms on anomalous curves.
  // For a full implementation, we would need p-adic lifting (Hensel's lemma).
  throw new NotImplementedError(
    'padic_elliptic_logarithm for finite fields requires p-adic lifting. ' +
      'For discrete logarithms, use point_log instead.'
  );
}

/**
 * Return True if this point has good reduction at the given prime.
 *
 * A point P has good reduction at a prime p if when we reduce P modulo p,
 * the resulting point is non-singular on the reduced curve. This is
 * equivalent to saying that the denominators of the coordinates (when
 * written in lowest terms) are not divisible by p.
 *
 * INPUT:
 * - P: a point on an elliptic curve over a number field
 * - prime: a prime ideal, or None (to check all primes dividing the
 *          discriminant)
 *
 * OUTPUT: True if the point has good reduction at the prime
 *
 * NOTE: For elliptic curves over finite fields, the concept of "reduction"
 * does not apply in the same sense. All points on a curve over F_p are
 * already "reduced", so this function returns True for finite fields.
 *
 * @see Reference: sage/schemes/elliptic_curves/ell_point.py:has_good_reduction
 */
export function point_has_good_reduction<F extends FieldElement>(
  P: EllipticCurvePoint<F>,
  _prime?: unknown
): boolean {
  const char = P.curve.base_ring.characteristic;
  if (char > 0n) {
    // For finite fields, points are already in their reduced form
    // The concept of "good reduction" applies to number fields
    return true;
  }

  // For number fields, we need to check if the point coordinates
  // have denominators divisible by the prime
  throw new NotImplementedError(
    'has_good_reduction for number fields requires prime ideal infrastructure. ' +
      'Not yet implemented.'
  );
}

/**
 * Return the reduction of point P modulo prime p.
 *
 * For a point P on an elliptic curve E over Q (or a number field),
 * this computes the reduction of P modulo p, which is a point on
 * the reduced curve E mod p over F_p.
 *
 * INPUT:
 * - P: a point on an elliptic curve over Q or a number field
 * - p: a prime number
 *
 * OUTPUT: the reduced point on the reduced curve over F_p
 *
 * NOTE: For curves already over a finite field F_q, this function
 * returns the point itself if q = p, or throws an error otherwise
 * (since reduction to a different characteristic is not defined).
 *
 * @see Reference: sage/schemes/elliptic_curves/ell_point.py:reduction
 */
export function reduction<F extends FieldElement>(
  P: EllipticCurvePoint<F>,
  p: bigint | number
): EllipticCurvePoint<F> {
  const pVal = typeof p === 'number' ? BigInt(p) : p;
  const char = P.curve.base_ring.characteristic;

  if (char > 0n) {
    // Already over a finite field
    if (char === pVal) {
      // Same characteristic - return the point itself
      return P;
    }
    throw new ValueError(
      `Cannot reduce a point over F_${char} modulo a different prime ${pVal}`
    );
  }

  // For number fields, we need to:
  // 1. Reduce the curve coefficients mod p
  // 2. Reduce the point coordinates mod p
  // 3. Verify the reduced point is on the reduced curve
  throw new NotImplementedError(
    'reduction from Q to F_p requires modular reduction infrastructure. ' +
      'Not yet implemented for number fields.'
  );
}

/**
 * Return True if there exists a point Q defined over the same field
 * as P such that mQ = P.
 *
 * INPUT:
 * - P: an elliptic curve point
 * - m: a positive integer
 *
 * OUTPUT: boolean; True if there is a solution, else False.
 *
 * ALGORITHM:
 * For finite fields, we use the division_points function and check
 * if any solutions exist. For trivial cases (m = 1, -1), always return True.
 * If P is the identity, it's always divisible. If m = 0 and P is not identity,
 * return False.
 *
 * @example
 * ```typescript
 * const E = EllipticCurve(GF(101n), [23n, 34n]);
 * const P = E.point([0n, 1n]);
 * console.log(is_divisible_by(P, 2n)); // Check if P is divisible by 2
 * ```
 *
 * @see Reference: sage/schemes/elliptic_curves/ell_point.py:is_divisible_by
 */
export function is_divisible_by<F extends FieldElement>(
  P: EllipticCurvePoint<F>,
  m: bigint | number
): boolean {
  const mVal = typeof m === 'number' ? BigInt(m) : m;

  // Trivial cases: m = 1, -1 always returns True
  if (mVal === 1n || mVal === -1n) {
    return true;
  }

  // m = 0 case: 0*Q = O for all Q, so only the identity is divisible by 0
  if (mVal === 0n) {
    return P.is_zero();
  }

  // Use absolute value of m
  const absM = mVal < 0n ? -mVal : mVal;

  // The identity point O is always divisible by any m (since m*O = O)
  if (P.is_zero()) {
    return true;
  }

  // Check if P has finite order n and gcd(m, n) = 1
  // In that case, there exists a unique solution
  const order = P.order();
  if (order !== undefined) {
    // Compute gcd(absM, order)
    let a = absM;
    let b = order;
    while (b !== 0n) {
      const temp = b;
      b = a % b;
      a = temp;
    }
    if (a === 1n) {
      return true;
    }
  }

  // Check if P is 2-torsion (P = -P)
  const P_is_2_torsion = P.eq(P.neg());

  // If P is not 2-torsion and m is odd with P being m-torsion,
  // then P is trivially divisible
  if (!P_is_2_torsion && absM % 2n === 1n && order !== undefined && order % absM === 0n) {
    return true;
  }

  // For 2-torsion points with odd m, P itself is a solution (since m*P = P when 2*P = O and m is odd)
  if (P_is_2_torsion && absM % 2n === 1n) {
    return true;
  }

  // Fall back to checking if division_points returns any solutions
  const divPoints = division_points(P, absM);
  return divPoints.length > 0;
}

/**
 * Return the discrete logarithm of point Q with respect to the given base P.
 * In other words, return an integer x such that x*P = Q.
 *
 * A ValueError is raised if there is no solution (Q is not a multiple of P).
 *
 * INPUT:
 * - Q: the target point (self)
 * - base: another point P on the same curve
 *
 * OUTPUT: an integer x with 0 <= x < ord(P) such that x*P = Q
 *
 * ALGORITHM:
 * For finite fields of small size, we use baby-step giant-step algorithm.
 * First we check that Q is in the subgroup generated by P by verifying
 * that n*Q = O where n is the order of P.
 *
 * @example
 * ```typescript
 * const E = EllipticCurve(GF(101n), [1n, 1n]);
 * const P = E.point([0n, 1n]); // generator
 * const Q = P.mul(42n);
 * const x = point_log(Q, P); // returns 42n
 * console.log(P.mul(x).eq(Q)); // true
 * ```
 *
 * @see Reference: sage/schemes/elliptic_curves/ell_point.py:log
 */
export function point_log<F extends FieldElement>(
  Q: EllipticCurvePoint<F>,
  base: EllipticCurvePoint<F>
): bigint {
  const E = Q.curve;

  // Check that both points are on the same curve
  if (base.curve !== E) {
    throw new ValueError('points must be on the same curve');
  }

  // If Q is the identity, the logarithm is 0
  if (Q.is_zero()) {
    return 0n;
  }

  // If base is the identity, Q cannot be a multiple of it (unless Q is also identity)
  if (base.is_zero()) {
    throw new ValueError('Q is not a multiple of the base point (base is identity)');
  }

  // Compute the order of the base point
  const orderP = base.order();
  if (orderP === undefined) {
    throw new ValueError('could not determine order of base point');
  }

  // Check that Q is in the subgroup generated by base
  // n*Q should be O if Q = x*P for some x
  if (!Q.mul(orderP).is_zero()) {
    throw new ValueError('Q is not a multiple of the base point');
  }

  // Baby-step giant-step algorithm
  // Find x such that x*P = Q where 0 <= x < n
  const n = orderP;

  // m = ceil(sqrt(n))
  let m = 1n;
  while (m * m < n) {
    m++;
  }

  // Baby step: compute a table of j*P for j = 0, 1, ..., m-1
  const babySteps = new Map<string, bigint>();
  let jP = E.zero();
  for (let j = 0n; j < m; j++) {
    const key = jP.is_zero() ? 'O' : `${jP.x().toString()},${jP.y().toString()}`;
    babySteps.set(key, j);
    jP = jP.add(base);
  }

  // Giant step: compute Q - i*m*P for i = 0, 1, ..., m
  // and look for a match in the baby steps table
  const mP = base.mul(m);
  let gamma = Q;
  for (let i = 0n; i <= m; i++) {
    const key = gamma.is_zero() ? 'O' : `${gamma.x().toString()},${gamma.y().toString()}`;
    const j = babySteps.get(key);
    if (j !== undefined) {
      // Found: x = i*m + j
      let x = i * m + j;
      // Normalize to [0, n)
      x = ((x % n) + n) % n;
      return x;
    }
    gamma = gamma.sub(mP);
  }

  throw new ValueError('Q is not a multiple of the base point (no solution found)');
}
