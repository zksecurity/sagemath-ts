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

import { NotImplementedError, TypeError as SageTypeError, ValueError } from '../../errors.js';
import {
  type OperationType,
  discrete_log as generic_discrete_log,
  has_order as generic_has_order,
  order_from_bounds,
} from '../../groups/generic.js';
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
          // `schemes/projective/projective_point.py`: SageMath raises a
          // TypeError naming the projective coordinates and the curve.
          throw new SageTypeError(
            `Coordinates [${x}, ${y}, ${curve.base_ring.one()}] do not define a point on ${curve}`
          );
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
          throw new SageTypeError(
            `Coordinates [${X}, ${Y}, ${Z}] do not define a point on ${curve}`
          );
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
   * Alias for is_zero() to satisfy AdditiveGroupElement interface.
   */
  isZero(): boolean {
    return this.is_zero();
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
    const two = field.__call__(2n) as F;
    const three = field.__call__(3n) as F;

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

  /** Cached order of this point */
  protected _order: bigint | undefined = undefined;

  /**
   * Compute the order of this point in the elliptic curve group.
   *
   * Uses the generic order_from_bounds algorithm which employs
   * baby-step giant-step (BSGS) for O(sqrt(n)) complexity.
   *
   * @param options - Configuration options
   * @param options.algorithm - Algorithm to use:
   *   - 'generic_small': Uses order_from_bounds with no bounds (gradually increases)
   *   - 'pari': Delegates to PARI (only for finite field subclasses)
   *   - 'hybrid': Combines generic_small with PARI when curve order is known
   * @returns The order of this point
   *
   * @example
   * ```typescript
   * const E = EllipticCurve(F, [a, b]);
   * const P = E.point([x, y]);
   * const ord = P.order(); // Uses BSGS algorithm
   * ```
   *
   * @see Reference: sage/schemes/elliptic_curves/ell_point.py:order
   */
  order(options?: { algorithm?: 'generic_small' | 'pari' | 'hybrid' }): bigint {
    // Return cached order if available
    if (this._order !== undefined) {
      return this._order;
    }

    // Identity has order 1
    if (this.is_zero()) {
      this._order = 1n;
      return 1n;
    }

    const algorithm = options?.algorithm ?? 'generic_small';

    if (algorithm === 'pari') {
      // PARI algorithm is only available for finite field subclasses
      throw new NotImplementedError(
        "algorithm 'pari' is only available for points on curves over finite fields"
      );
    }

    if (algorithm === 'generic_small' || algorithm === 'hybrid') {
      // Use order_from_bounds which employs BSGS for O(sqrt(n)) complexity
      // With no bounds provided, it will gradually increase the search range
      this._order = order_from_bounds(this, undefined, undefined, '+' as OperationType);
      return this._order;
    }

    throw new NotImplementedError(`algorithm '${algorithm}' not implemented`);
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
      // Verify n*P = O
      if (!this.mul(n).is_zero()) {
        throw new ValueError(`${n} is not the order of this point`);
      }
      // Verify no smaller divisor works (check prime divisors)
      // We would need factor() here but for now just do basic check
      if (n > 1n && this.mul(n / 2n).is_zero() && n % 2n === 0n) {
        throw new ValueError(`${n} is not the minimal order of this point`);
      }
    }
    this._order = n;
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
  const three = P.curve.base_ring.__call__(3n) as F;
  const two = P.curve.base_ring.__call__(2n) as F;

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
  const kVal = typeof k === 'number' ? BigInt(k) : k;

  // Sage requires a finite base field.
  const cardinality = fieldCardinality(field);
  if (cardinality === undefined) {
    throw new NotImplementedError(
      'Reduced Tate pairing is currently only implemented for finite fields'
    );
  }

  const p = field.characteristic;
  // K.degree(): the degree of K over its prime field.
  let d = 1n;
  {
    let acc = p;
    while (acc < cardinality) {
      acc *= p;
      d++;
    }
  }

  let qVal: bigint;
  if (q === undefined) {
    if (d === 1n) {
      qVal = cardinality;
    } else if (d === kVal) {
      qVal = p;
    } else {
      throw new ValueError(
        'Unexpected field degree: set keyword argument q equal to the size of the base field ' +
          `(big field is GF(q^${kVal})).`
      );
    }
  } else {
    // The user has supplied q, so we check here that it is a sensible value:
    // Mod(q, n)**k != 1
    qVal = q;
    if (modPow(((qVal % n) + n) % n, kVal, n) !== 1n % n) {
      throw new ValueError('n does not divide (q^k - 1) for the supplied value of q');
    }
  }

  if (!P.mul(n).is_zero()) {
    throw new ValueError('The point P must be n-torsion');
  }

  // The Miller value; a ZeroDivisionError propagates as in Sage's
  // EllipticCurve_finite_field branch (which calls PARI's elltatepairing).
  const ePQ = _miller(P, Q, n);

  const exp = (qVal ** kVal - 1n) / n;
  return ePQ.pow(exp) as F;
}

/**
 * Return the ate pairing of the n-torsion points P and Q.
 *
 * P must be GF(q)-rational (i.e. in ker(pi - 1)) and Q must lie in
 * ker(pi - q), where pi is the q-power Frobenius.
 *
 * INPUT:
 * - P: point of order n in ker(pi - 1)
 * - Q: point of order n in ker(pi - q)
 * - n: the order of P and Q
 * - k: the embedding degree
 * - t: the trace of Frobenius
 * - q: (optional) size of base field
 *
 * OUTPUT: an n-th root of unity in F_q^k
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
  const O = E.zero();

  if (Q.curve !== E) {
    throw new ValueError('Points must both be on the same curve');
  }

  const field = E.base_ring;
  const kVal = typeof k === 'number' ? BigInt(k) : k;

  // set q to be the order of the base field
  let qVal: bigint;
  if (q === undefined) {
    const cardinality = fieldCardinality(field);
    if (cardinality === undefined) {
      throw new ValueError(
        'Unexpected field degree: set keyword argument q equal to the size of the base field ' +
          `(big field is GF(q^${kVal})).`
      );
    }
    const p = field.characteristic;
    let d = 1n;
    let acc = p;
    while (acc < cardinality) {
      acc *= p;
      d++;
    }
    if (d === kVal) {
      qVal = p;
    } else {
      throw new ValueError(
        'Unexpected field degree: set keyword argument q equal to the size of the base field ' +
          `(big field is GF(q^${kVal})).`
      );
    }
  } else {
    qVal = q;
  }

  // check order of P
  if (!P.mul(n).eq(O)) {
    throw new ValueError(`This point ${P} is not of order n=${n}`);
  }

  // check for P in kernel pi - 1
  const piP = frobeniusImage(P, qVal);
  if (!piP.sub(P).eq(O)) {
    throw new ValueError(`This point ${P} is not in Ker(pi - 1)`);
  }

  // check for Q in kernel pi - q
  const piQ = frobeniusImage(Q, qVal);
  if (!piQ.sub(Q.mul(qVal)).eq(O)) {
    throw new ValueError(`Point ${Q} not in Ker(pi - q)`);
  }

  const T = t - 1n;
  // Sage passes the *signed* T to _miller_, which handles T < 0 by returning
  // 1/(v_{TQ} * f_{T,Q}) -- the vertical-line factor must not be dropped.
  let ret = _miller(Q, P, T);
  const e = (qVal ** kVal - 1n) / n;
  ret = ret.pow(e) as F;
  return ret;
}

/**
 * Return the cardinality of a finite field, or undefined if it is not finite.
 */
function fieldCardinality(field: FieldRing): bigint | undefined {
  const K = field as unknown as {
    cardinality?: () => bigint;
    order?: bigint | number;
    degree?: number;
    characteristic: bigint;
  };
  if (K.characteristic === 0n) {
    return undefined;
  }
  if (typeof K.cardinality === 'function') {
    return K.cardinality();
  }
  if (K.order !== undefined) {
    return typeof K.order === 'number' ? BigInt(K.order) : K.order;
  }
  if (K.degree !== undefined) {
    return K.characteristic ** BigInt(K.degree);
  }
  return K.characteristic;
}

/**
 * Apply the q-power Frobenius to a point: (x, y) -> (x^q, y^q).
 */
function frobeniusImage<F extends FieldElement>(
  P: EllipticCurvePoint<F>,
  q: bigint
): EllipticCurvePoint<F> {
  if (P.is_zero()) {
    return P;
  }
  return affinePoint(P.curve, P.x().pow(q) as F, P.y().pow(q) as F, false);
}

/**
 * Modular exponentiation on plain integers.
 */
function modPow(base: bigint, exp: bigint, mod: bigint): bigint {
  if (mod === 1n) return 0n;
  let result = 1n;
  let b = ((base % mod) + mod) % mod;
  let e = exp;
  while (e > 0n) {
    if (e & 1n) result = (result * b) % mod;
    b = (b * b) % mod;
    e >>= 1n;
  }
  return result;
}

/**
 * Structural view of the curve methods that ``division_points`` needs.
 *
 * ``ell_point.ts`` cannot import ``ell_generic.ts`` (circular dependency), so
 * the extra methods are accessed through this interface.
 */
interface CurveWithDivisionPolynomials<F extends FieldElement> extends EllipticCurveInterface<F> {
  division_polynomial(m: bigint | number): { roots(): Array<[unknown, number]> };
  _multiple_x_numerator(n: bigint | number): unknown;
  _multiple_x_denominator(n: bigint | number): unknown;
  is_x_coord(x: F): boolean;
  lift_x(x: F): EllipticCurvePoint<F>;
}

/**
 * Compare two points the way SageMath orders them: by (Z, X, Y), so that the
 * point at infinity is always the minimum.
 *
 * @see Reference: sage/schemes/elliptic_curves/ell_point.py:_richcmp_
 */
function comparePoints<F extends FieldElement>(
  P: EllipticCurvePoint<F>,
  Q: EllipticCurvePoint<F>
): number {
  const [Px, Py, Pz] = P.xyz();
  const [Qx, Qy, Qz] = Q.xyz();
  for (const [a, b] of [
    [Pz, Qz],
    [Px, Qx],
    [Py, Qy],
  ] as Array<[F, F]>) {
    const c = compareFieldElements(a, b);
    if (c !== 0) return c;
  }
  return 0;
}

/**
 * Compare two field elements the way SageMath orders them (integer
 * representative for prime fields, coefficient vector otherwise).
 */
function compareFieldElements(a: FieldElement, b: FieldElement): number {
  const av = (a as unknown as { value?: unknown }).value;
  const bv = (b as unknown as { value?: unknown }).value;
  if (typeof av === 'bigint' && typeof bv === 'bigint') {
    return av < bv ? -1 : av > bv ? 1 : 0;
  }
  const as = a.toString();
  const bs = b.toString();
  return as < bs ? -1 : as > bs ? 1 : 0;
}

/**
 * Return a list of all points Q such that mQ = P where P = self.
 *
 * Only points on the elliptic curve containing self and defined
 * over the base field are included.
 *
 * INPUT:
 * - P: the target point
 * - m: a nonzero integer
 * - poly_only: if true, return the polynomial whose roots are the
 *   x-coordinates of the solutions instead of the points
 *
 * OUTPUT: a sorted list of points Q such that mQ = P
 *
 * ALGORITHM:
 * Compute a polynomial g whose roots are exactly the possible x-coordinates
 * of the m-division points (via ``_multiple_x_numerator`` /
 * ``_multiple_x_denominator``, or the m-division polynomial when P is the
 * identity), take its **distinct** roots, lift each one to a single point Q,
 * and keep Q or -Q according to whether mQ equals P or -P.
 *
 * @see Reference: sage/schemes/elliptic_curves/ell_point.py:division_points
 */
export function division_points<F extends FieldElement>(
  P: EllipticCurvePoint<F>,
  m: bigint | number,
  poly_only?: false
): EllipticCurvePoint<F>[];
export function division_points<F extends FieldElement>(
  P: EllipticCurvePoint<F>,
  m: bigint | number,
  poly_only: true
): { roots(): Array<[unknown, number]> };
export function division_points<F extends FieldElement>(
  P: EllipticCurvePoint<F>,
  m: bigint | number,
  poly_only?: boolean
): EllipticCurvePoint<F>[] | { roots(): Array<[unknown, number]> } {
  const mVal = typeof m === 'number' ? BigInt(m) : m;

  // Check for trivial cases of m = 1, -1 and 0.
  if (mVal === 1n || mVal === -1n) {
    return [P.mul(mVal)];
  }
  if (mVal === 0n) {
    // then every point Q is a solution, but ...
    return P.is_zero() ? [P] : [];
  }

  const ans: EllipticCurvePoint<F>[] = [];

  const E = P.curve as CurveWithDivisionPolynomials<F>;
  if (
    typeof E.division_polynomial !== 'function' ||
    typeof E._multiple_x_numerator !== 'function' ||
    typeof E.lift_x !== 'function'
  ) {
    throw new NotImplementedError(
      'division_points requires a curve providing division_polynomial, ' +
        '_multiple_x_numerator/_multiple_x_denominator, is_x_coord and lift_x'
    );
  }

  const nP = P.neg();
  const P_is_2_torsion = P.eq(nP);

  // If self is 0, then self is a solution, and the correct poly is the m'th
  // division polynomial.
  let g: { roots(): Array<[unknown, number]> };
  if (P.is_zero()) {
    ans.push(P);
    g = E.division_polynomial(mVal < 0n ? -mVal : mVal);
  } else {
    // The poly g here is 0 at x(Q) iff x(m*Q) = x(P).
    const absM = mVal < 0n ? -mVal : mVal;
    const num = E._multiple_x_numerator(absM) as unknown as PolyLike;
    const den = E._multiple_x_denominator(absM) as unknown as PolyLike;
    g = num.sub(den.mul(den.parent.__call__(P.x()))) as unknown as {
      roots(): Array<[unknown, number]>;
    };

    // Sage additionally replaces g by its square root when 2*P = 0 (see
    // ell_point.py:1531-1557). That step only removes repeated factors, so
    // the *set* of roots -- all we use below -- is unchanged; we skip it.
  }

  if (poly_only) {
    return g;
  }

  for (const [xRoot] of g.roots()) {
    const x = xRoot as F;
    if (!E.is_x_coord(x)) {
      continue;
    }
    // Make a point on the curve with this x coordinate.
    const Q = E.lift_x(x);
    const nQ = Q.neg();
    const mQ = Q.mul(mVal);
    // if P == -P then Q works iff -Q works, so we include both unless they
    // are equal:
    if (P_is_2_torsion) {
      if (mQ.eq(P)) {
        ans.push(Q);
        if (!nQ.eq(Q)) {
          ans.push(nQ);
        }
      }
    } else {
      // P is not 2-torsion so at most one of Q, -Q works and we must try both:
      if (mQ.eq(P)) {
        ans.push(Q);
      } else if (mQ.eq(nP)) {
        ans.push(nQ);
      }
    }
  }

  // Finally, sort and return
  ans.sort(comparePoints);
  return ans;
}

/** Minimal structural view of a univariate polynomial. */
interface PolyLike {
  readonly parent: { __call__(x: unknown): PolyLike };
  sub(other: PolyLike): PolyLike;
  mul(other: PolyLike): PolyLike;
  roots(): Array<[unknown, number]>;
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
    throw new ValueError(`Cannot reduce a point over F_${char} modulo a different prime ${pVal}`);
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

  // Check for trivial cases of m = 1, -1 and 0.
  if (mVal === 1n || mVal === -1n) {
    return true;
  }
  if (mVal === 0n) {
    return P.is_zero(); // then m*self = self for all m!
  }
  const absM = mVal < 0n ? -mVal : mVal;

  // If P has finite order n and gcd(m, n) = 1 then the result is True.
  const n = P.order();
  let a = absM;
  let b = n;
  while (b !== 0n) {
    const t = b;
    b = a % b;
    a = t;
  }
  if (a === 1n) {
    return true;
  }

  const P_is_2_torsion = P.eq(P.neg());
  const g = division_points(P, absM, true);

  if (!P_is_2_torsion) {
    // In this case deg(g) = m^2, and each root in K lifts to two points
    // Q, -Q both in E(K), of which exactly one is a solution. So we just
    // check the existence of roots.
    return g.roots().length > 0;
  }

  // Now 2*P == 0.
  if (absM % 2n === 1n) {
    return true; // P itself is a solution when m is odd
  }

  // Now m is even and 2*P = 0. Roots of g in K may or may not lift to
  // solutions in E(K), so we fall back to the default.
  return division_points(P, absM).length > 0;
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
    throw new ValueError('not a point on the same curve');
  }

  // n = base.order()
  const n = base.order();

  // Sage: "if (hasattr(self, '_order') and not self._order.divides(n)) or n*self"
  if (!Q.mul(n).is_zero()) {
    throw new ValueError('ECDLog problem has no solution (order does not divide order of base)');
  }

  if (Q.is_zero()) {
    return 0n;
  }
  if (base.is_zero()) {
    // n == 1, and Q is nonzero: handled by the n*self check above.
    throw new ValueError('ECDLog problem has no solution (order does not divide order of base)');
  }

  // Sage rules out the case where Q lies outside <base> but still satisfies
  // n*Q = O, using the Weil pairing (ell_point.py:4640-4646).
  // (Sage's ``E._order.gcd(n**2) == n`` shortcut is only an optimisation.)
  const p = fieldCardinality(E.base_ring);
  if (p !== undefined && p === E.base_ring.characteristic && n === p) {
    // Anomalous case. Sage calls ``base.padic_elliptic_logarithm(self, p)``;
    // the port's p-adic logarithm is an unimplemented stub, so we fall through
    // to the generic algorithm, which returns the same value more slowly.
    // @see Deviation: anomalous ECDLog has no p-adic shortcut
  } else if (!weil_pairing(base, Q, n).eq(E.base_ring.one())) {
    throw new ValueError('ECDLog problem has no solution (non-trivial Weil pairing)');
  }

  // Sage delegates the actual logarithm to PARI's ``elllog``, which runs
  // Pohlig-Hellman with baby-step/giant-step inside each prime-power factor.
  // Our generic port of that algorithm lives in ``groups/generic.ts``.
  return generic_discrete_log(Q, base, n, '+' as OperationType);
}
