/**
 * @module sage/schemes/elliptic_curves/ell_generic
 * @description Base elliptic curve class over a general ring
 *
 * Port of: sage/schemes/elliptic_curves/ell_generic.py
 *
 * An elliptic curve is defined by a Weierstrass equation:
 *   y^2 + a1*x*y + a3*y = x^3 + a2*x^2 + a4*x + a6
 *
 * This is the "long Weierstrass form". The "short Weierstrass form"
 *   y^2 = x^3 + a*x + b
 * corresponds to a1 = a2 = a3 = 0, a4 = a, a6 = b.
 */

import { ArithmeticError, ValueError, ZeroDivisionError } from '../../errors.js';
import type { MPolynomial } from '../../rings/polynomial/multi_polynomial_element.js';
import { MPolynomialRing } from '../../rings/polynomial/multi_polynomial_ring.js';
import { Polynomial, type RingElement } from '../../rings/polynomial/polynomial_element.js';
import { type CoefficientRing, PolynomialRing } from '../../rings/polynomial/polynomial_ring.js';
import {
  type EllipticCurveInterface,
  type EllipticCurvePoint,
  type FieldElement,
  type FieldParent,
  affinePoint,
  pointAtInfinity,
} from './ell_point.js';

/**
 * Compute the GCD of two bigints.
 */
/**
 * Compare two field elements the way SageMath orders them, so that lists of
 * isomorphisms can be sorted exactly as ``sorted(...)`` does upstream.
 *
 * For prime fields, Sage compares the integer representatives; for extension
 * fields it compares the coefficient vectors. We use the numeric ``value``
 * when the element exposes one and fall back to the string form otherwise.
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

function gcd(a: bigint, b: bigint): bigint {
  a = a < 0n ? -a : a;
  b = b < 0n ? -b : b;
  while (b !== 0n) {
    const t = b;
    b = a % b;
    a = t;
  }
  return a;
}

/**
 * An elliptic curve over a ring/field.
 *
 * The curve is defined by the Weierstrass equation:
 *   y^2 + a1*x*y + a3*y = x^3 + a2*x^2 + a4*x + a6
 *
 * @template F - The type of field elements
 *
 * @example
 * ```typescript
 * const F = GF(23n);
 * const E = new EllipticCurveGeneric(F, [F(0), F(0), F(0), F(1), F(1)]);
 *
 * // Get invariants
 * console.log(E.discriminant());
 * console.log(E.j_invariant());
 *
 * // Create a point
 * const P = E.point([F(0), F(1)]);
 *
 * // Point arithmetic
 * const Q = P.add(P);  // 2P
 * const R = P.mul(5n); // 5P
 * ```
 */
import { NotImplementedError } from '../../errors.js';
import {
  EllipticCurveTorsionSubgroup,
  _p_primary_torsion_basis as _p_primary_torsion_basis_impl,
} from './ell_torsion.js';
import { EllipticCurveFormalGroup } from './formal_group.js';

/**
 * Forward declaration for EllipticCurveIsogeny to avoid circular dependency.
 * The actual implementation is in ell_curve_isogeny.ts.
 */
interface EllipticCurveIsogeny<F extends FieldElement> {
  domain(): EllipticCurveGeneric<F>;
  codomain(): EllipticCurveGeneric<F>;
  degree(): bigint;
  call(P: EllipticCurvePoint<F>): EllipticCurvePoint<F>;
  evaluate(P: EllipticCurvePoint<F>): EllipticCurvePoint<F>;
  kernel_polynomial(): bigint[];
  scaling_factor(): F;
  dual(): EllipticCurveIsogeny<F>;
  neg(): unknown;
  eq(other: EllipticCurveIsogeny<F>): boolean;
  is_separable(): boolean;
  is_normalized(): boolean;
  is_injective(): boolean;
  is_surjective(): boolean;
}

// Module-level variables for lazy loading to avoid circular dependency
let _isogenyModule: typeof import('./ell_curve_isogeny.js') | null = null;

async function getIsogenyModule() {
  if (_isogenyModule === null) {
    _isogenyModule = await import('./ell_curve_isogeny.js');
  }
  return _isogenyModule;
}

// Synchronous helper using dynamic import workaround
// Note: This creates an isogeny synchronously by constructing directly
function createIsogeny<F extends FieldElement>(
  E: EllipticCurveGeneric<F>,
  kernel: EllipticCurvePoint<F> | EllipticCurvePoint<F>[]
): EllipticCurveIsogeny<F> {
  // We need to construct the isogeny inline to avoid circular deps
  // This is a simplified inline implementation using Velu's formulas
  throw new NotImplementedError(
    'isogeny() requires importing from ell_curve_isogeny. ' +
      'Use: import { EllipticCurveIsogeny } from "./ell_curve_isogeny.js" and construct directly.'
  );
}

function isogenies_prime_degree_helper<F extends FieldElement>(
  _E: EllipticCurveGeneric<F>,
  _l?: bigint | number | Array<bigint | number>
): EllipticCurveIsogeny<F>[] {
  throw new NotImplementedError(
    'isogenies_prime_degree() requires importing from ell_curve_isogeny. ' +
      'Use: import { isogenies_prime_degree } from "./ell_curve_isogeny.js"'
  );
}

export class EllipticCurveGeneric<F extends FieldElement = FieldElement>
  implements EllipticCurveInterface<F>
{
  /** The base ring/field */
  readonly base_ring: FieldParent;

  /** Weierstrass coefficients [a1, a2, a3, a4, a6] */
  private readonly _ainvs: readonly [F, F, F, F, F];

  /** Cached b-invariants */
  private _binvs: [F, F, F, F] | null = null;

  /** Cached c-invariants */
  private _cinvs: [F, F] | null = null;

  /** Cached discriminant */
  private _discriminant: F | null = null;

  /** Cached j-invariant */
  private _j_invariant: F | null = null;

  /** The point at infinity */
  private _infinity: EllipticCurvePoint<F> | null = null;

  /**
   * Create an elliptic curve from Weierstrass coefficients.
   *
   * @param K - The base ring/field
   * @param ainvs - The Weierstrass coefficients [a1, a2, a3, a4, a6]
   */
  constructor(K: FieldParent, ainvs: readonly [F, F, F, F, F]) {
    this.base_ring = K;
    this._ainvs = ainvs;

    // Check that discriminant is non-zero
    const disc = this.discriminant();
    if (disc.isZero()) {
      throw new ArithmeticError(`${this._equation_string()} defines a singular curve`);
    }
  }

  /**
   * Get a1 coefficient.
   */
  a1(): F {
    return this._ainvs[0];
  }

  /**
   * Get a2 coefficient.
   */
  a2(): F {
    return this._ainvs[1];
  }

  /**
   * Get a3 coefficient.
   */
  a3(): F {
    return this._ainvs[2];
  }

  /**
   * Get a4 coefficient.
   */
  a4(): F {
    return this._ainvs[3];
  }

  /**
   * Get a6 coefficient.
   */
  a6(): F {
    return this._ainvs[4];
  }

  /**
   * Return the a-invariants [a1, a2, a3, a4, a6].
   */
  a_invariants(): [F, F, F, F, F] {
    return [this._ainvs[0], this._ainvs[1], this._ainvs[2], this._ainvs[3], this._ainvs[4]];
  }

  /**
   * Alias for a_invariants.
   */
  ainvs(): [F, F, F, F, F] {
    return this.a_invariants();
  }

  /**
   * Return the b-invariants [b2, b4, b6, b8].
   *
   * These are computed from the a-invariants:
   *   b2 = a1^2 + 4*a2
   *   b4 = a1*a3 + 2*a4
   *   b6 = a3^2 + 4*a6
   *   b8 = a1^2*a6 + 4*a2*a6 - a1*a3*a4 + a2*a3^2 - a4^2
   */
  b_invariants(): [F, F, F, F] {
    if (this._binvs !== null) {
      return this._binvs;
    }

    const [a1, a2, a3, a4, a6] = this._ainvs;
    const K = this.base_ring;
    const two = K.__call__(2n) as F;
    const four = K.__call__(4n) as F;

    // b2 = a1^2 + 4*a2
    const b2 = a1.mul(a1).add(a2.mul(four)) as F;

    // b4 = a1*a3 + 2*a4
    const b4 = a1.mul(a3).add(a4.mul(two)) as F;

    // b6 = a3^2 + 4*a6
    const b6 = a3.mul(a3).add(a6.mul(four)) as F;

    // b8 = a1^2*a6 + 4*a2*a6 - a1*a3*a4 + a2*a3^2 - a4^2
    const b8 = a1
      .mul(a1)
      .mul(a6)
      .add(a2.mul(a6).mul(four))
      .sub(a1.mul(a3).mul(a4))
      .add(a2.mul(a3).mul(a3))
      .sub(a4.mul(a4)) as F;

    this._binvs = [b2, b4, b6, b8];
    return this._binvs;
  }

  /**
   * Return b2.
   */
  b2(): F {
    return this.b_invariants()[0];
  }

  /**
   * Return b4.
   */
  b4(): F {
    return this.b_invariants()[1];
  }

  /**
   * Return b6.
   */
  b6(): F {
    return this.b_invariants()[2];
  }

  /**
   * Return b8.
   */
  b8(): F {
    return this.b_invariants()[3];
  }

  /**
   * Return the c-invariants [c4, c6].
   *
   * These are computed from the b-invariants:
   *   c4 = b2^2 - 24*b4
   *   c6 = -b2^3 + 36*b2*b4 - 216*b6
   */
  c_invariants(): [F, F] {
    if (this._cinvs !== null) {
      return this._cinvs;
    }

    const [b2, b4, b6] = this.b_invariants();
    const K = this.base_ring;
    const n24 = K.__call__(24n) as F;
    const n36 = K.__call__(36n) as F;
    const n216 = K.__call__(216n) as F;

    // c4 = b2^2 - 24*b4
    const c4 = b2.mul(b2).sub(b4.mul(n24)) as F;

    // c6 = -b2^3 + 36*b2*b4 - 216*b6
    const c6 = b2.mul(b2).mul(b2).neg().add(b2.mul(b4).mul(n36)).sub(b6.mul(n216)) as F;

    this._cinvs = [c4, c6];
    return this._cinvs;
  }

  /**
   * Return c4.
   */
  c4(): F {
    return this.c_invariants()[0];
  }

  /**
   * Return c6.
   */
  c6(): F {
    return this.c_invariants()[1];
  }

  /**
   * Return the discriminant of this elliptic curve.
   *
   * The discriminant is:
   *   Delta = -b2^2*b8 - 8*b4^3 - 27*b6^2 + 9*b2*b4*b6
   *
   * A curve is non-singular if and only if Delta != 0.
   */
  discriminant(): F {
    if (this._discriminant !== null) {
      return this._discriminant;
    }

    const [b2, b4, b6, b8] = this.b_invariants();
    const K = this.base_ring;
    const n8 = K.__call__(8n) as F;
    const n9 = K.__call__(9n) as F;
    const n27 = K.__call__(27n) as F;

    // Delta = -b2^2*b8 - 8*b4^3 - 27*b6^2 + 9*b2*b4*b6
    const disc = b2
      .mul(b2)
      .mul(b8)
      .neg()
      .sub(b4.mul(b4).mul(b4).mul(n8))
      .sub(b6.mul(b6).mul(n27))
      .add(b2.mul(b4).mul(b6).mul(n9)) as F;

    this._discriminant = disc;
    return disc;
  }

  /**
   * Return the j-invariant of this elliptic curve.
   *
   * The j-invariant is:
   *   j = c4^3 / Delta
   *
   * @throws {ZeroDivisionError} If the discriminant is zero (singular curve)
   */
  j_invariant(): F {
    if (this._j_invariant !== null) {
      return this._j_invariant;
    }

    const c4 = this.c4();
    const disc = this.discriminant();

    if (disc.isZero()) {
      throw new ZeroDivisionError('j-invariant is not defined for singular curves');
    }

    const j = c4.mul(c4).mul(c4).div(disc) as F;
    this._j_invariant = j;
    return j;
  }

  /**
   * Check if a point (x, y) is on this curve.
   *
   * A point (x, y) is on the curve if:
   *   y^2 + a1*x*y + a3*y = x^3 + a2*x^2 + a4*x + a6
   */
  is_on_curve(x: F, y: F): boolean {
    const [a1, a2, a3, a4, a6] = this._ainvs;

    // LHS: y^2 + a1*x*y + a3*y
    const lhs = y.mul(y).add(a1.mul(x).mul(y)).add(a3.mul(y));

    // RHS: x^3 + a2*x^2 + a4*x + a6
    const rhs = x.mul(x).mul(x).add(a2.mul(x).mul(x)).add(a4.mul(x)).add(a6);

    return lhs.eq(rhs);
  }

  /**
   * Return True if ``x`` is the x-coordinate of a rational point on this curve.
   *
   * @see Reference: sage/schemes/elliptic_curves/ell_generic.py:is_x_coord
   */
  is_x_coord(x: F | bigint | number): boolean {
    const K = this.base_ring;
    const xx = (typeof x === 'bigint' || typeof x === 'number' ? K.__call__(x) : x) as F;
    const [a1, a2, a3, a4, a6] = this._ainvs;
    const fx = xx.add(a2).mul(xx).add(a4).mul(xx).add(a6) as F;
    if (a1.isZero() && a3.isZero()) {
      return this._is_square(fx);
    }
    const b = a1.mul(xx).add(a3) as F;
    if (K.characteristic === 2n) {
      // Roots of y^2 + b*y - fx over K.
      return this._poly_roots([fx.neg() as F, b, K.one() as F]).length > 0;
    }
    const four = K.__call__(4n) as F;
    const D = b.mul(b).add(four.mul(fx)) as F;
    return this._is_square(D);
  }

  /**
   * Return one or all points with the given x-coordinate.
   *
   * @param x - the x-coordinate
   * @param all - if true return the (possibly empty) list of all such points
   *
   * @throws {ValueError} if ``all`` is false and there is no such point
   *
   * @see Reference: sage/schemes/elliptic_curves/ell_generic.py:lift_x
   */
  lift_x(x: F | bigint | number, all?: false): EllipticCurvePoint<F>;
  lift_x(x: F | bigint | number, all: true): EllipticCurvePoint<F>[];
  lift_x(
    x: F | bigint | number,
    all: boolean = false
  ): EllipticCurvePoint<F> | EllipticCurvePoint<F>[] {
    const K = this.base_ring;
    const xx = (typeof x === 'bigint' || typeof x === 'number' ? K.__call__(x) : x) as F;
    const [a1, a2, a3, a4, a6] = this._ainvs;
    const b = a1.mul(xx).add(a3) as F;
    const f = xx.add(a2).mul(xx).add(a4).mul(xx).add(a6) as F;

    let ys: F[];
    if (K.characteristic === 2n) {
      ys = this._poly_roots([f.neg() as F, b, K.one() as F]);
    } else {
      const two = K.__call__(2n) as F;
      const four = K.__call__(4n) as F;
      const D = b.mul(b).add(four.mul(f)) as F;
      ys = this._square_roots(D).map((d) => b.neg().add(d).div(two) as F);
    }

    // "ys.sort()  # ensure deterministic behavior"
    ys.sort((p, q) => compareFieldElements(p, q));
    // Remove duplicates (D == 0 yields the same y twice).
    const uniq: F[] = [];
    for (const y of ys) {
      if (uniq.length === 0 || !uniq[uniq.length - 1]!.eq(y)) {
        uniq.push(y);
      }
    }

    if (uniq.length > 0) {
      if (all) {
        return uniq.map((y) => this.point([xx, y], false));
      }
      return this.point([xx, uniq[0]!], false);
    }

    if (all) {
      return [];
    }
    throw new ValueError(`No point with x-coordinate ${xx} on ${this}`);
  }

  /**
   * Return the point at infinity (identity element of the group).
   */
  point_at_infinity(): EllipticCurvePoint<F> {
    if (this._infinity === null) {
      this._infinity = pointAtInfinity(this);
    }
    return this._infinity;
  }

  /**
   * Alias for point_at_infinity().
   */
  zero(): EllipticCurvePoint<F> {
    return this.point_at_infinity();
  }

  /**
   * Create a point on this curve from affine coordinates.
   *
   * @param coords - Coordinates [x, y] or [] for point at infinity
   * @param check - Whether to verify the point is on the curve (default: true)
   */
  point(coords: [F, F] | [], check: boolean = true): EllipticCurvePoint<F> {
    if (coords.length === 0) {
      return this.point_at_infinity();
    }
    const [x, y] = coords;
    return affinePoint(this, x, y, check);
  }

  /**
   * Check if this curve is in short Weierstrass form: y^2 = x^3 + a*x + b.
   */
  is_short_weierstrass(): boolean {
    const [a1, a2, a3] = this._ainvs;
    return a1.isZero() && a2.isZero() && a3.isZero();
  }

  /**
   * Return the string representation of the equation.
   */
  private _equation_string(): string {
    const [a1, a2, a3, a4, a6] = this._ainvs;

    let lhs = 'y^2';
    if (!a1.isZero()) {
      lhs += ` + ${a1}*x*y`;
    }
    if (!a3.isZero()) {
      lhs += ` + ${a3}*y`;
    }

    let rhs = 'x^3';
    if (!a2.isZero()) {
      rhs += ` + ${a2}*x^2`;
    }
    if (!a4.isZero()) {
      rhs += ` + ${a4}*x`;
    }
    if (!a6.isZero()) {
      rhs += ` + ${a6}`;
    }

    return `${lhs} = ${rhs}`;
  }

  /**
   * String representation.
   */
  toString(): string {
    return `Elliptic Curve defined by ${this._equation_string()} over ${this.base_ring}`;
  }

  /**
   * Return the n-th torsion (division) polynomial, without the 2-torsion factor
   * if n is even, as a polynomial in x.
   *
   * These are the polynomials g_n from [MT1991], but with sign flipped for even n
   * so that the leading coefficient is always positive.
   *
   * This is an internal method; users should use division_polynomial().
   *
   * INPUT:
   * - n: positive integer, or special values -1, -2
   *   -1 means B_6 = (2y + a1*x + a3)^2 = 4x^3 + b2*x^2 + 2*b4*x + b6
   *   -2 means B_6^2
   * - x: (optional) ring element to use as the x variable
   *
   * @see Reference: sage/schemes/elliptic_curves/ell_generic.py:division_polynomial_0
   */
  division_polynomial_0<T extends RingElement>(
    n: bigint | number | (bigint | number)[],
    x?: T | Polynomial<T>
  ): Polynomial<T> | Polynomial<T>[] {
    const K = this.base_ring;

    // Create polynomial ring if x is not provided
    let polyRing: PolynomialRing<T>;
    let xVar: Polynomial<T>;

    if (x === undefined) {
      // Create a polynomial ring over the base field
      // Need to cast K to CoefficientRing
      polyRing = new PolynomialRing(
        K as unknown as { zero(): T; one(): T; __call__(x: unknown): T },
        'x'
      );
      xVar = polyRing.gen();
    } else if (x instanceof Polynomial) {
      polyRing = x.parent as PolynomialRing<T>;
      xVar = x;
    } else {
      // x is a ring element - create polynomial ring and wrap it
      polyRing = new PolynomialRing(
        K as unknown as { zero(): T; one(): T; __call__(x: unknown): T },
        'x'
      );
      xVar = polyRing.__call__(x);
    }

    const [b2, b4, b6, b8] = this.b_invariants();

    // Cache for recursive computation
    const cache: Map<number, Polynomial<T>> = new Map();

    const poly = (n: number): Polynomial<T> => {
      const cached = cache.get(n);
      if (cached !== undefined) {
        return cached;
      }

      let ret: Polynomial<T>;

      if (n === -2) {
        // B_6^2
        ret = poly(-1).mul(poly(-1));
      } else if (n === -1) {
        // B_6 = 4x^3 + b2*x^2 + 2*b4*x + b6
        const four = polyRing.__call__(K.__call__(4n) as unknown as T);
        const two = K.__call__(2n) as F;
        const b2Poly = polyRing.__call__(b2 as unknown as T);
        const twoB4 = polyRing.__call__((b4 as FieldElement).mul(two) as unknown as T);
        const b6Poly = polyRing.__call__(b6 as unknown as T);

        ret = four
          .mul(xVar.pow(3))
          .add(b2Poly.mul(xVar.pow(2)))
          .add(twoB4.mul(xVar))
          .add(b6Poly);
      } else if (n <= 0) {
        throw new ValueError('n must be a positive integer (or -1 or -2)');
      } else if (n === 1 || n === 2) {
        ret = polyRing.one();
      } else if (n === 3) {
        // 3x^4 + b2*x^3 + 3*b4*x^2 + 3*b6*x + b8
        const three = polyRing.__call__(K.__call__(3n) as unknown as T);
        const threeF = K.__call__(3n) as F;
        const b2Poly = polyRing.__call__(b2 as unknown as T);
        const threeB4 = polyRing.__call__((b4 as FieldElement).mul(threeF) as unknown as T);
        const threeB6 = polyRing.__call__((b6 as FieldElement).mul(threeF) as unknown as T);
        const b8Poly = polyRing.__call__(b8 as unknown as T);

        ret = three
          .mul(xVar.pow(4))
          .add(b2Poly.mul(xVar.pow(3)))
          .add(threeB4.mul(xVar.pow(2)))
          .add(threeB6.mul(xVar))
          .add(b8Poly);
      } else if (n === 4) {
        // -B_6^2 + (6x^2 + b2*x + b4) * psi_3
        const negB6sq = poly(-2).neg();
        const six = polyRing.__call__(K.__call__(6n) as unknown as T);
        const b2Poly = polyRing.__call__(b2 as unknown as T);
        const b4Poly = polyRing.__call__(b4 as unknown as T);

        const factor = six.mul(xVar.pow(2)).add(b2Poly.mul(xVar)).add(b4Poly);

        ret = negB6sq.add(factor.mul(poly(3)));
      } else if (n % 2 === 0) {
        // Even n: psi_{m+1} * (psi_{m+3} * psi_m^2 - psi_{m-1} * psi_{m+2}^2)
        // where m = (n-2)/2
        const m = Math.floor((n - 2) / 2);
        ret = poly(m + 1).mul(
          poly(m + 3)
            .mul(poly(m).pow(2))
            .sub(poly(m - 1).mul(poly(m + 2).pow(2)))
        );
      } else {
        // Odd n
        const m = Math.floor((n - 1) / 2);
        if (m % 2 === 0) {
          // B_6^2 * psi_{m+2} * psi_m^3 - psi_{m-1} * psi_{m+1}^3
          ret = poly(-2)
            .mul(poly(m + 2))
            .mul(poly(m).pow(3))
            .sub(poly(m - 1).mul(poly(m + 1).pow(3)));
        } else {
          // psi_{m+2} * psi_m^3 - B_6^2 * psi_{m-1} * psi_{m+1}^3
          ret = poly(m + 2)
            .mul(poly(m).pow(3))
            .sub(
              poly(-2)
                .mul(poly(m - 1))
                .mul(poly(m + 1).pow(3))
            );
        }
      }

      cache.set(n, ret);
      return ret;
    };

    // Handle list input
    if (Array.isArray(n)) {
      return n.map((k) => poly(Number(k)));
    }

    return poly(Number(n));
  }

  /**
   * Return the m-th division polynomial of this elliptic curve evaluated at x.
   *
   * The division polynomial is cached if x is None.
   *
   * INPUT:
   * - m: positive integer
   * - x: (optional) ring element to use as the x variable
   * - two_torsion_multiplicity: 0, 1, or 2 (default 2)
   *   - If 0: For even m, omits factors whose roots are x-coords of 2-torsion points
   *   - If 2: For even m, includes a factor of degree 3 whose roots are x-coords of 2-torsion
   *   - If 1: For even m, includes factor 2y + a1*x + a3
   *
   * OUTPUT: a polynomial
   *
   * @see Reference: sage/schemes/elliptic_curves/ell_generic.py:division_polynomial
   */
  division_polynomial<T extends RingElement>(
    m: bigint | number,
    x?: T | Polynomial<T>,
    two_torsion_multiplicity?: 0 | 2
  ): Polynomial<T>;
  division_polynomial<T extends RingElement>(
    m: bigint | number,
    x: undefined,
    two_torsion_multiplicity: 1
  ): MPolynomial<T>;
  division_polynomial<T extends RingElement>(
    m: bigint | number,
    x: [T, T] | EllipticCurvePoint<F>,
    two_torsion_multiplicity: 1
  ): T;
  division_polynomial<T extends RingElement>(
    m: bigint | number,
    x?: T | Polynomial<T> | [T, T] | EllipticCurvePoint<F>,
    two_torsion_multiplicity: number = 2
  ): Polynomial<T> | MPolynomial<T> | T {
    if (![0, 1, 2].includes(two_torsion_multiplicity)) {
      throw new ValueError('two_torsion_multiplicity must be 0, 1, or 2');
    }

    let xy: [T, T] | undefined;
    if (x !== undefined && two_torsion_multiplicity === 1) {
      // Sage accepts a point (which is converted to its (x,y) tuple) or a
      // 2-tuple; anything else is an error.
      let cand: unknown = x;
      if (
        cand !== null &&
        typeof cand === 'object' &&
        typeof (cand as { xy?: unknown }).xy === 'function'
      ) {
        cand = (cand as EllipticCurvePoint<F>).xy();
      }
      if (!Array.isArray(cand) || cand.length !== 2) {
        throw new ValueError(
          'x should be a tuple of length 2 (or None) when two_torsion_multiplicity is 1'
        );
      }
      xy = cand as [T, T];
    }

    const mNum = Number(m);

    if (mNum <= 0) {
      throw new ValueError('m must be a positive integer');
    }

    if (two_torsion_multiplicity === 0) {
      return this.division_polynomial_0<T>(
        mNum,
        x as T | Polynomial<T> | undefined
      ) as Polynomial<T>;
    }

    if (two_torsion_multiplicity === 2) {
      const f = this.division_polynomial_0<T>(
        mNum,
        x as T | Polynomial<T> | undefined
      ) as Polynomial<T>;
      if (mNum % 2 === 0) {
        return f.mul(
          this.division_polynomial_0<T>(-1, x as T | Polynomial<T> | undefined) as Polynomial<T>
        );
      }
      return f;
    }

    // two_torsion_multiplicity === 1:
    //   f = psi_m, times (2*y + a1*x + a3) when m is even.
    const [a1, , a3] = this._ainvs;

    if (xy === undefined) {
      // Return a bivariate polynomial in x, y over the base ring.
      const R = new MPolynomialRing(this.base_ring as unknown as CoefficientRing<T>, ['x', 'y']);
      const X = R.gen(0);
      const Y = R.gen(1);

      // Lift the univariate psi_m(x) into R.
      const psi = this.division_polynomial_0<T>(mNum) as Polynomial<T>;
      let f = R.zero();
      for (let i = 0; i <= psi.degree(); i++) {
        const c = psi.getCoeff(i);
        if ((c as unknown as FieldElement).isZero()) continue;
        f = f.add(X.pow(i).scalarMul(c));
      }

      if (mNum % 2 === 0) {
        const two = this.base_ring.__call__(2n) as unknown as T;
        const lin = Y.scalarMul(two)
          .add(X.scalarMul(a1 as unknown as T))
          .add(R.__call__(a3));
        f = f.mul(lin);
      }
      return f as unknown as MPolynomial<T>;
    }

    // Evaluate the bivariate polynomial at the given (x, y).
    const [xv, yv] = xy;
    const f = this.division_polynomial_0<T>(mNum, xv) as Polynomial<T>;
    let val = f.getCoeff(0) as unknown as FieldElement;
    if (mNum % 2 === 0) {
      const two = this.base_ring.__call__(2n);
      const lin = two
        .mul(yv as unknown as FieldElement)
        .add(a1.mul(xv as unknown as FieldElement))
        .add(a3);
      val = val.mul(lin);
    }
    return val as unknown as T;
  }

  /**
   * Return the 2-division polynomial of this elliptic curve evaluated at x.
   *
   * The 2-division polynomial is 4x^3 + b2*x^2 + 2*b4*x + b6.
   * Its roots are the x-coordinates of the points of order 2.
   *
   * INPUT:
   * - x: (optional) ring element to use as the x variable
   *
   * OUTPUT: a polynomial in x
   *
   * @example
   * ```typescript
   * const E = EllipticCurve(GF(23n), [1n, 1n]);
   * const f = E.two_division_polynomial();
   * // f = 4*x^3 + 2*x + 4 for y^2 = x^3 + x + 1
   * ```
   *
   * @see Reference: sage/schemes/elliptic_curves/ell_generic.py:two_division_polynomial
   */
  two_division_polynomial<T extends RingElement>(x?: T | Polynomial<T>): Polynomial<T> {
    // The 2-division polynomial is division_polynomial_0(-1, x)
    // which computes B_6 = 4x^3 + b2*x^2 + 2*b4*x + b6
    return this.division_polynomial_0<T>(-1, x) as Polynomial<T>;
  }

  /** Cached formal group */
  private _formal_group: EllipticCurveFormalGroup | null = null;

  /**
   * Return the formal group associated to this elliptic curve.
   *
   * The formal group of an elliptic curve E is a one-dimensional
   * formal group law F(X, Y) = X + Y + higher order terms that encodes
   * the addition law on E near the origin.
   *
   * OUTPUT: an EllipticCurveFormalGroup object
   *
   * EXAMPLES:
   * ```typescript
   * const E = EllipticCurve(GF(23n), [1n, 1n]);
   * const F = E.formal_group();
   * const w = F.w(10);  // Power series w(t)
   * ```
   *
   * @see Reference: sage/schemes/elliptic_curves/ell_generic.py:formal_group
   */
  formal_group(): EllipticCurveFormalGroup {
    if (this._formal_group === null) {
      // Cast this curve to the interface expected by EllipticCurveFormalGroup
      this._formal_group = new EllipticCurveFormalGroup(
        this as unknown as Parameters<typeof EllipticCurveFormalGroup.prototype.constructor>[0]
      );
    }
    return this._formal_group;
  }

  /**
   * Alias for formal_group().
   *
   * @see Reference: sage/schemes/elliptic_curves/ell_generic.py:formal
   */
  formal(): EllipticCurveFormalGroup {
    return this.formal_group();
  }

  /**
   * Return an isogeny from this curve with the given kernel.
   *
   * INPUT:
   * - kernel: either a point on this curve, or a list of points generating the kernel
   * - codomain: (optional) the codomain curve (for verification)
   * - degree: (optional) the degree of the isogeny (for verification)
   *
   * OUTPUT: an EllipticCurveIsogeny from this curve
   *
   * The isogeny is computed using Velu's formulas when the kernel is given
   * as points.
   *
   * NOTE: Due to circular dependency issues, this method requires explicit import.
   * Use: `import { EllipticCurveIsogeny } from './ell_curve_isogeny.js'`
   * Then: `new EllipticCurveIsogeny(E, kernel)`
   *
   * @example
   * ```typescript
   * import { EllipticCurveIsogeny } from './ell_curve_isogeny.js';
   * const E = EllipticCurve(GF(23n), [1n, 1n]);
   * const P = E.point([...]);  // a point of finite order
   * const phi = new EllipticCurveIsogeny(E, P);
   * const E2 = phi.codomain();
   * ```
   *
   * @see Reference: sage/schemes/elliptic_curves/ell_generic.py:isogeny
   */
  isogeny(
    _kernel: EllipticCurvePoint<F> | EllipticCurvePoint<F>[],
    _codomain?: EllipticCurveGeneric<F>,
    _degree?: bigint | number
  ): EllipticCurveIsogeny<F> {
    return createIsogeny(this, _kernel);
  }

  /**
   * Return a list of all isogenies of a given prime degree with domain equal to this curve.
   *
   * INPUT:
   * - l: a prime number, or a list of primes, or undefined to use small primes [2, 3, 5, 7]
   *
   * OUTPUT: a list of EllipticCurveIsogeny objects
   *
   * ALGORITHM:
   * For each prime l, find all l-torsion points on the curve and construct
   * cyclic isogenies from them.
   *
   * NOTE: Due to circular dependency issues, this method requires explicit import.
   * Use: `import { isogenies_prime_degree } from './ell_curve_isogeny.js'`
   *
   * @example
   * ```typescript
   * import { isogenies_prime_degree } from './ell_curve_isogeny.js';
   * const E = EllipticCurve(GF(101n), [1n, 1n]);
   * const isogs = isogenies_prime_degree(E, 2n);
   * // isogs contains all degree-2 isogenies from E
   * ```
   *
   * @see Reference: sage/schemes/elliptic_curves/ell_generic.py:isogenies_prime_degree
   */
  isogenies_prime_degree(l?: bigint | number | Array<bigint | number>): EllipticCurveIsogeny<F>[] {
    return isogenies_prime_degree_helper(this, l);
  }

  /**
   * Return the torsion points on this elliptic curve.
   *
   * For curves over finite fields, this returns all points (since all points
   * have finite order). For small fields, this uses brute force enumeration.
   *
   * OUTPUT: a list of all torsion points on this curve
   *
   * @example
   * ```typescript
   * const E = EllipticCurve(GF(7n), [1n, 1n]);
   * const tors = E.torsion_points();
   * // tors contains all points on the curve over GF(7)
   * ```
   *
   * @see Reference: sage/schemes/elliptic_curves/ell_generic.py:torsion_points
   */
  torsion_points(): EllipticCurvePoint<F>[] {
    const K = this.base_ring;
    const p = K.characteristic;
    const results: EllipticCurvePoint<F>[] = [];

    // Always include the point at infinity
    results.push(this.point_at_infinity());

    // For finite fields, enumerate all points
    // This is only practical for small fields
    const maxFieldSize = 10000n;

    if (p > maxFieldSize) {
      throw new ValueError(
        `torsion_points() is only implemented for fields of size <= ${maxFieldSize}. ` +
          'Use _p_primary_torsion_basis for larger fields.'
      );
    }

    // Enumerate all x values and find corresponding y values
    for (let xVal = 0n; xVal < p; xVal++) {
      const x = K.__call__(xVal) as F;

      // Compute RHS: x^3 + a2*x^2 + a4*x + a6
      const [a1, a2, a3, a4, a6] = this._ainvs;
      const x2 = x.mul(x) as F;
      const x3 = x2.mul(x) as F;

      const rhs = x3.add(a2.mul(x2)).add(a4.mul(x)).add(a6) as F;

      // For general Weierstrass, we need to solve:
      // y^2 + (a1*x + a3)*y - RHS = 0
      // This is a quadratic in y: y^2 + By + C = 0 where B = a1*x + a3, C = -RHS
      const B = a1.mul(x).add(a3) as F;
      const C = rhs.neg() as F;

      // Solutions: y = (-B +/- sqrt(B^2 - 4C)) / 2
      // In characteristic 2, this is different

      if (p === 2n) {
        // Characteristic 2 case is special
        // y^2 + By + C = 0 => y^2 + By = -C
        // Try y = 0 and y = 1
        for (let yVal = 0n; yVal < 2n; yVal++) {
          const y = K.__call__(yVal) as F;
          if (this.is_on_curve(x, y)) {
            results.push(this.point([x, y], false));
          }
        }
      } else {
        // Characteristic != 2
        // Discriminant: B^2 - 4C = B^2 + 4*RHS
        const four = K.__call__(4n) as F;
        const two = K.__call__(2n) as F;
        const disc = B.mul(B).add(four.mul(rhs)) as F;

        // Check if discriminant is a quadratic residue
        const sqrtDisc = this._square_roots(disc);

        // Only use the first square root (if any) to avoid duplicates
        // since _square_roots returns both +s and -s
        if (sqrtDisc.length > 0) {
          const sd = sqrtDisc[0]!;
          // y = (-B + sd) / 2 and y = (-B - sd) / 2
          const y1 = B.neg().add(sd).div(two) as F;
          const y2 = B.neg().sub(sd).div(two) as F;

          if (this.is_on_curve(x, y1)) {
            results.push(this.point([x, y1], false));
          }

          // Only add y2 if it's different from y1
          if (!y1.eq(y2) && this.is_on_curve(x, y2)) {
            results.push(this.point([x, y2], false));
          }
        }
      }
    }

    return results;
  }

  /**
   * Return the torsion subgroup of this elliptic curve.
   *
   * OUTPUT: the torsion subgroup as an EllipticCurveTorsionSubgroup object
   *
   * For curves over finite fields, this returns the full group of points
   * (since all points have finite order over finite fields).
   *
   * @example
   * ```typescript
   * const E = EllipticCurve(GF(7n), [1n, 1n]);
   * const T = E.torsion_subgroup();
   * console.log(T.order());    // number of points
   * console.log(T.gens());     // generators
   * ```
   *
   * @see Reference: sage/schemes/elliptic_curves/ell_generic.py:torsion_subgroup
   */
  torsion_subgroup(): EllipticCurveTorsionSubgroup<F> {
    return new EllipticCurveTorsionSubgroup(this);
  }

  /**
   * Return True if this elliptic curve has good reduction at prime p.
   *
   * INPUT:
   * - p: a prime number, or None to test all primes dividing discriminant
   *
   * OUTPUT: True if E has good reduction at p
   *
   * NOTE: This method is only applicable to elliptic curves over number fields
   * (including Q). For curves over finite fields, reduction at primes does not
   * apply in the same way. The curve is defined over a finite field, not a
   * number field, so we check if the discriminant is non-zero modulo p.
   *
   * @example
   * ```typescript
   * // For a curve over Q (not implemented yet)
   * // const E = EllipticCurve(QQ, [0, 0, 0, -1, 0]);
   * // E.has_good_reduction(2n); // False (bad reduction at 2)
   * ```
   *
   * @see Reference: sage/schemes/elliptic_curves/ell_number_field.py:has_good_reduction
   */
  has_good_reduction(p?: bigint | number): boolean {
    const K = this.base_ring;
    const char = K.characteristic;

    // For finite fields, this concept doesn't apply directly
    if (char > 0n) {
      // If the curve is over Fp, we can check if it's non-singular (which it is)
      // For checking reduction at primes different from p, we need number field infrastructure
      if (p === undefined) {
        // Check if the discriminant is non-zero (already guaranteed by constructor)
        return !this.discriminant().isZero();
      }

      const pVal = BigInt(p);
      if (pVal === char) {
        // The curve is defined over Fp, so "reduction at p" is the curve itself
        // Good reduction means the discriminant is non-zero mod p
        return !this.discriminant().isZero();
      }

      throw new ValueError(
        `has_good_reduction at prime ${p} is not applicable for curves over F_${char}. ` +
          `The curve is already defined over a finite field of characteristic ${char}.`
      );
    }

    // For number fields (characteristic 0), requires Tate's algorithm
    throw new NotImplementedError(
      "has_good_reduction over number fields requires Tate's algorithm (local_data)"
    );
  }

  /**
   * Return True if this elliptic curve has bad reduction at prime p.
   *
   * INPUT:
   * - p: a prime number
   *
   * OUTPUT: True if E has bad reduction at p
   *
   * NOTE: This method is only applicable to elliptic curves over number fields
   * (including Q). For curves over finite fields, the concept of reduction at
   * primes different from the field characteristic does not apply.
   *
   * @example
   * ```typescript
   * // For a curve over Q (not implemented yet)
   * // const E = EllipticCurve(QQ, [0, 0, 0, -1, 0]);
   * // E.has_bad_reduction(2n); // True (bad reduction at 2)
   * ```
   *
   * @see Reference: sage/schemes/elliptic_curves/ell_number_field.py:has_bad_reduction
   */
  has_bad_reduction(p?: bigint | number): boolean {
    const K = this.base_ring;
    const char = K.characteristic;

    // For finite fields, this concept doesn't apply directly
    if (char > 0n) {
      if (p === undefined) {
        // Check if the discriminant is zero (which it never is for valid curves)
        return this.discriminant().isZero();
      }

      const pVal = BigInt(p);
      if (pVal === char) {
        // The curve is defined over Fp, so "reduction at p" is the curve itself
        // Bad reduction means the discriminant is zero mod p
        return this.discriminant().isZero();
      }

      throw new ValueError(
        `has_bad_reduction at prime ${p} is not applicable for curves over F_${char}. ` +
          `The curve is already defined over a finite field of characteristic ${char}.`
      );
    }

    // For number fields (characteristic 0), requires Tate's algorithm
    throw new NotImplementedError(
      "has_bad_reduction over number fields requires Tate's algorithm (local_data)"
    );
  }

  /**
   * Return the conductor of this elliptic curve.
   *
   * OUTPUT: the conductor as a positive integer (for curves over Q)
   *         or as a fractional ideal (for curves over number fields)
   *
   * The conductor measures the arithmetic complexity of the elliptic curve.
   * It is the product of local factors from primes of bad reduction.
   *
   * NOTE: The conductor is only defined for elliptic curves over number fields
   * (including Q). For curves over finite fields, this method throws an error.
   *
   * For curves over Q, the conductor equals the product of local conductors:
   *   N = prod_p p^{f_p}
   * where f_p is the conductor exponent at p, computed via Tate's algorithm.
   *
   * @example
   * ```typescript
   * // For a curve over Q (not fully implemented)
   * // const E = EllipticCurve(QQ, [0, -1, 1, 0, 0]);  // 11a1
   * // E.conductor();  // 11n
   * ```
   *
   * @see Reference: sage/schemes/elliptic_curves/ell_number_field.py:conductor
   */
  conductor(): bigint {
    const K = this.base_ring;
    const char = K.characteristic;

    // For finite fields, conductor is not defined
    if (char > 0n) {
      throw new ValueError(
        'conductor is not defined for elliptic curves over finite fields. ' +
          `The curve is defined over a field of characteristic ${char}.`
      );
    }

    // For number fields (characteristic 0), requires Tate's algorithm
    throw new NotImplementedError(
      "conductor over number fields requires local_data (Tate's algorithm)"
    );
  }

  /**
   * Return the local data at prime p.
   *
   * INPUT:
   * - p: a prime number, or None to return local data at all primes of bad reduction
   *
   * OUTPUT: local data object containing:
   *   - Kodaira symbol (I_n, II, III, IV, I_n*, II*, III*, IV*)
   *   - Conductor exponent (f_p)
   *   - Tamagawa number (c_p)
   *   - Minimal discriminant valuation
   *   - Reduction type (good, multiplicative, additive)
   *
   * NOTE: This method is only applicable to elliptic curves over number fields
   * (including Q). For curves over finite fields, local data does not apply.
   *
   * For curves over Q, local data at p is computed using PARI's ellglobalred
   * or Tate's algorithm.
   *
   * @example
   * ```typescript
   * // For a curve over Q (not fully implemented)
   * // const E = EllipticCurve(QQ, [0, -1, 1, 0, 0]);  // 11a1
   * // const ld = E.local_data(11n);
   * // console.log(ld.kodaira_symbol());  // I_5
   * // console.log(ld.tamagawa_number()); // 5
   * ```
   *
   * @see Reference: sage/schemes/elliptic_curves/ell_number_field.py:local_data
   * @see Reference: sage/schemes/elliptic_curves/ell_local_data.py:EllipticCurveLocalData
   */
  local_data(_p?: bigint | number): unknown {
    const K = this.base_ring;
    const char = K.characteristic;

    // For finite fields, local data is not defined
    if (char > 0n) {
      throw new ValueError(
        'local_data is not defined for elliptic curves over finite fields. ' +
          `The curve is defined over a field of characteristic ${char}.`
      );
    }

    // For number fields (characteristic 0), requires Tate's algorithm
    throw new NotImplementedError(
      "local_data over number fields requires Tate's algorithm. " +
        'This is implemented in sage/schemes/elliptic_curves/ell_local_data.py'
    );
  }

  /**
   * Extend the base field of this curve to R.
   *
   * INPUT:
   * - R: a ring or field that is an extension of the current base ring
   *
   * OUTPUT: an elliptic curve over R with the same a-invariants
   *
   * @see Reference: sage/schemes/elliptic_curves/ell_generic.py:base_extend
   */
  base_extend<G extends FieldElement>(R: FieldRing): EllipticCurveGeneric<G> {
    // Convert a-invariants to the new ring
    const [a1, a2, a3, a4, a6] = this.a_invariants();

    // Coerce the a-invariants into the new ring
    // This assumes the ring R has a __call__ method that can coerce elements
    const newA1 = R.__call__(a1.value ?? a1) as G;
    const newA2 = R.__call__(a2.value ?? a2) as G;
    const newA3 = R.__call__(a3.value ?? a3) as G;
    const newA4 = R.__call__(a4.value ?? a4) as G;
    const newA6 = R.__call__(a6.value ?? a6) as G;

    return new EllipticCurveGeneric(R, [newA1, newA2, newA3, newA4, newA6]);
  }

  /**
   * Return a new curve with base ring changed to R.
   *
   * This has the same effect as base_extend(R).
   *
   * INPUT:
   * - R: a ring
   *
   * OUTPUT: an elliptic curve over R
   *
   * @see Reference: sage/schemes/elliptic_curves/ell_generic.py:change_ring
   */
  change_ring<G extends FieldElement>(R: FieldRing): EllipticCurveGeneric<G> {
    return this.base_extend(R);
  }

  /**
   * Return a new Weierstrass model under the standard transformation (u,r,s,t).
   *
   * The transformation is: (x,y) -> (u^2*x' + r, u^3*y' + s*u^2*x' + t)
   *
   * INPUT:
   * - u, r, s, t: ring elements (u must be non-zero)
   *
   * OUTPUT: the transformed curve
   *
   * @see Reference: sage/schemes/elliptic_curves/ell_generic.py:change_weierstrass_model
   */
  change_weierstrass_model(u: F, r: F, s: F, t: F): EllipticCurveGeneric<F> {
    if (u.isZero()) {
      throw new ValueError('u must be non-zero for Weierstrass isomorphism');
    }

    const K = this.base_ring;

    // Get original a-invariants
    let [a1, a2, a3, a4, a6] = this.a_invariants();

    // Apply transformation formulas from SageMath's weierstrass_morphism.py
    // The transformation (u,r,s,t) maps (x,y) -> (u^2*x' + r, u^3*y' + s*u^2*x' + t)
    // New a-invariants are computed as:
    const two = K.__call__(2n) as F;
    const three = K.__call__(3n) as F;

    // a6' = a6 + r*(a4 + r*(a2 + r)) - t*(a3 + r*a1 + t)
    const ra2r = r.mul(a2.add(r) as F) as F;
    const a4_plus_ra2r = a4.add(ra2r) as F;
    const rPart = r.mul(a4_plus_ra2r) as F;
    const a3_plus_ra1_plus_t = a3.add(r.mul(a1) as F).add(t) as F;
    const tPart = t.mul(a3_plus_ra1_plus_t) as F;
    a6 = a6.add(rPart).sub(tPart) as F;

    // a4' = a4 - s*a3 + 2*r*a2 - (t + r*s)*a1 + 3*r^2 - 2*s*t
    const sa3 = s.mul(a3) as F;
    const twoRA2 = two.mul(r).mul(a2) as F;
    const t_plus_rs = t.add(r.mul(s) as F) as F;
    const t_plus_rs_a1 = t_plus_rs.mul(a1) as F;
    const threeRR = three.mul(r).mul(r) as F;
    const twoST = two.mul(s).mul(t) as F;
    a4 = a4.sub(sa3).add(twoRA2).sub(t_plus_rs_a1).add(threeRR).sub(twoST) as F;

    // a3' = a3 + r*a1 + 2*t
    const ra1 = r.mul(a1) as F;
    const twoT = two.mul(t) as F;
    a3 = a3.add(ra1).add(twoT) as F;

    // a2' = a2 - s*a1 + 3*r - s^2
    const sa1 = s.mul(a1) as F;
    const threeR = three.mul(r) as F;
    const ss = s.mul(s) as F;
    a2 = a2.sub(sa1).add(threeR).sub(ss) as F;

    // a1' = a1 + 2*s
    const twoS = two.mul(s) as F;
    a1 = a1.add(twoS) as F;

    // Divide by appropriate powers of u
    const u2 = u.mul(u) as F;
    const u3 = u2.mul(u) as F;
    const u4 = u2.mul(u2) as F;
    const u6 = u3.mul(u3) as F;

    // Use division if available, otherwise throw for non-fields
    const divU = (x: F): F => {
      if ('div' in x && typeof (x as unknown as { div: (y: F) => F }).div === 'function') {
        return (x as unknown as { div: (y: F) => F }).div(u);
      }
      // For fields, use inverse
      if ('inv' in u && typeof (u as unknown as { inv: () => F }).inv === 'function') {
        return x.mul((u as unknown as { inv: () => F }).inv()) as F;
      }
      throw new ValueError('Division not supported in base ring');
    };

    const divBy = (x: F, divisor: F): F => {
      if ('div' in x && typeof (x as unknown as { div: (y: F) => F }).div === 'function') {
        return (x as unknown as { div: (y: F) => F }).div(divisor);
      }
      if ('inv' in divisor && typeof (divisor as unknown as { inv: () => F }).inv === 'function') {
        return x.mul((divisor as unknown as { inv: () => F }).inv()) as F;
      }
      throw new ValueError('Division not supported in base ring');
    };

    const newA1 = divBy(a1, u);
    const newA2 = divBy(a2, u2);
    const newA3 = divBy(a3, u3);
    const newA4 = divBy(a4, u4);
    const newA6 = divBy(a6, u6);

    return new EllipticCurveGeneric(K, [newA1, newA2, newA3, newA4, newA6]);
  }

  /**
   * Apply an (r,s,t) transformation to this curve.
   *
   * This is a special case of change_weierstrass_model with u=1.
   *
   * INPUT:
   * - r, s, t: ring elements
   *
   * OUTPUT: the transformed curve
   *
   * @see Reference: sage/schemes/elliptic_curves/ell_generic.py:rst_transform
   */
  rst_transform(r: F, s: F, t: F): EllipticCurveGeneric<F> {
    return this.change_weierstrass_model(this.base_ring.one() as F, r, s, t);
  }

  /**
   * Return a short Weierstrass model for this curve.
   *
   * INPUT:
   * - complete_cube: boolean (default true)
   *   If true: Return a model of the form y^2 = x^3 + a*x + b
   *   If false: Return a model of the form y^2 = x^3 + ax^2 + bx + c
   *
   * OUTPUT: an elliptic curve in short Weierstrass form
   *
   * The characteristic must not be 2. In characteristic 3, it is only
   * possible to get the complete short form if b2 = 0.
   *
   * @example
   * ```typescript
   * const E = EllipticCurve(GF(23n), [1n, 2n, 3n, 4n, 5n]);
   * const F = E.short_weierstrass_model();
   * // F is now in the form y^2 = x^3 + a*x + b
   * console.log(E.is_isomorphic(F)); // true
   * ```
   *
   * @see Reference: sage/schemes/elliptic_curves/ell_generic.py:short_weierstrass_model
   */
  short_weierstrass_model(complete_cube: boolean = true): EllipticCurveGeneric<F> {
    const K = this.base_ring;
    const char = K.characteristic;

    // Any curve of the form y^2 = x^3 + .. is singular in characteristic 2
    if (char === 2n) {
      throw new ValueError(
        `short_weierstrass_model(): no short model for ${this} (characteristic is 2)`
      );
    }

    const [b2, b4, b6] = this.b_invariants();
    const [a1, a2, a3] = [this.a1(), this.a2(), this.a3()];

    // In characteristic 3, we can complete the square but can only complete
    // the cube if b2 is 0
    if (char === 3n) {
      if (complete_cube && !b2.isZero()) {
        throw new ValueError(
          `short_weierstrass_model(): no short model for ${this} (characteristic is 3)`
        );
      }
      // Return y^2 = x^3 + b2*x^2 + 8*b4*x + 16*b6
      const eight = K.__call__(8n) as F;
      const sixteen = K.__call__(16n) as F;
      const zero = K.zero() as F;

      return new EllipticCurveGeneric(K, [
        zero,
        b2 as F,
        zero,
        b4.mul(eight) as F,
        b6.mul(sixteen) as F,
      ]);
    }

    // Characteristic not 2 or 3
    if (complete_cube) {
      // Check if already in short Weierstrass form
      if (a1.isZero() && a2.isZero() && a3.isZero()) {
        return this;
      }

      if (b2.isZero()) {
        // If b2 = 0, use simpler transformation
        // Return y^2 = x^3 + 8*b4*x + 16*b6
        const eight = K.__call__(8n) as F;
        const sixteen = K.__call__(16n) as F;
        const zero = K.zero() as F;

        return new EllipticCurveGeneric(K, [
          zero,
          zero,
          zero,
          b4.mul(eight) as F,
          b6.mul(sixteen) as F,
        ]);
      } else {
        // General case: use c-invariants
        // Return y^2 = x^3 - 27*c4*x - 54*c6
        const [c4, c6] = this.c_invariants();
        const neg27 = K.__call__(-27n) as F;
        const neg54 = K.__call__(-54n) as F;
        const zero = K.zero() as F;

        return new EllipticCurveGeneric(K, [
          zero,
          zero,
          zero,
          c4.mul(neg27) as F,
          c6.mul(neg54) as F,
        ]);
      }
    } else {
      // complete_cube = false: return y^2 = x^3 + b2*x^2 + 8*b4*x + 16*b6
      if (a1.isZero() && a3.isZero()) {
        return this;
      }

      const eight = K.__call__(8n) as F;
      const sixteen = K.__call__(16n) as F;
      const zero = K.zero() as F;

      return new EllipticCurveGeneric(K, [
        zero,
        b2 as F,
        zero,
        b4.mul(eight) as F,
        b6.mul(sixteen) as F,
      ]);
    }
  }

  /**
   * Return the Montgomery model of this curve.
   *
   * A Montgomery curve has the form: By^2 = x^3 + Ax^2 + x
   * An untwisted Montgomery curve has B=1: y^2 = x^3 + Ax^2 + x
   *
   * INPUT:
   * - twisted: boolean (default: false); allow B != 1
   * - morphism: boolean (default: false); also return isomorphism
   *
   * OUTPUT:
   * - If twisted is false: an EllipticCurveGeneric in untwisted Montgomery form
   * - If morphism is true: a tuple [curve, isomorphism]
   *
   * Not all curves have a Montgomery model over their field of definition.
   *
   * ALGORITHM: [CS2018], Section 2.4
   *
   * @example
   * ```typescript
   * const F = GF(431n * 431n);  // Extension field
   * const E = EllipticCurve(F, [7n, 7n]);
   * const M = E.montgomery_model();  // y^2 = x^3 + A*x^2 + x
   * ```
   *
   * @see Reference: sage/schemes/elliptic_curves/ell_generic.py:montgomery_model
   */
  montgomery_model(
    twisted: boolean = false,
    morphism: boolean = false
  ): EllipticCurveGeneric<F> | [EllipticCurveGeneric<F>, [F, F, F, F]] {
    // Ew: y^2 = x^3 + a*x + b
    const Ew = this.short_weierstrass_model();
    const a = Ew.a4();
    const b = Ew.a6();
    const K = this.base_ring;
    const zero = K.zero() as F;
    const one = K.one() as F;
    const three = K.__call__(3n) as F;

    // sols = [(r, s) for r in P([b, a, 0, 1]).roots()
    //                for s in P([3*r^2 + a, 0, -1]).roots()]
    const sols: Array<[F, F]> = [];
    for (const r of this._poly_roots([b, a, zero, one])) {
      const c = three.mul(r).mul(r).add(a) as F;
      for (const s of this._poly_roots([c, zero, one.neg() as F])) {
        sols.push([r, s]);
      }
    }

    if (sols.length === 0) {
      throw new ValueError(`${this} has no Montgomery model`);
    }

    // "square s allows us to take B=1":
    //     r, s = max(sols, key=lambda t: t[1].is_square())
    // Python's max returns the *first* element attaining the maximum.
    let best = sols[0]!;
    for (const sol of sols) {
      if (this._is_square(sol[1])) {
        best = sol;
        break;
      }
    }
    const [r, s] = best;

    const A = three.mul(r).div(s) as F;
    const sIsSquare = this._is_square(s);
    const B: F = sIsSquare ? one : (one.div(s) as F);

    if (!twisted) {
      if (!B.eq(one)) {
        throw new ValueError(`${this} has no untwisted Montgomery model`);
      }
      const E = new EllipticCurveGeneric<F>(K, [zero, A, zero, one, zero]);
      if (morphism) {
        return [E, this.isomorphism_to(E)];
      }
      return E;
    }

    // Twisted case: B*y^2 = x^3 + A*x^2 + x is a plane cubic, not a
    // Weierstrass model, so it cannot be returned as an EllipticCurveGeneric.
    throw new NotImplementedError(
      'Twisted Montgomery models (B != 1) are not supported as EllipticCurve objects. ' +
        'They would need to be represented as ProjectivePlaneCurve.'
    );
  }

  /**
   * Return the cardinality of the base field, if it is finite.
   */
  private _field_order(): bigint {
    const K = this.base_ring as unknown as {
      cardinality?: () => bigint;
      order?: bigint | number;
      degree?: number;
      characteristic: bigint;
    };
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
   * Return whether ``x`` is a square in the base field (Euler's criterion).
   */
  private _is_square(x: F): boolean {
    if (x.isZero()) {
      return true;
    }
    const q = this._field_order();
    if (q % 2n === 0n) {
      // Every element of a field of characteristic 2 is a square.
      return true;
    }
    return x.pow((q - 1n) / 2n).eq(this.base_ring.one());
  }

  /**
   * Return the square roots of ``x`` in the base field (possibly empty).
   *
   * Uses Tonelli-Shanks over the full field cardinality; this is what
   * ``FiniteFieldElement.sqrt`` does in Sage and is O(log q) rather than a
   * root-finding call.
   */
  private _square_roots(x: F): F[] {
    const K = this.base_ring;
    if (x.isZero()) {
      return [K.zero() as F];
    }
    const q = this._field_order();
    if (q % 2n === 0n) {
      // Squaring is the Frobenius, hence a bijection: the unique square root
      // is x^(q/2).
      return [x.pow(q / 2n) as F];
    }
    if (!this._is_square(x)) {
      return [];
    }
    if (q % 4n === 3n) {
      const r = x.pow((q + 1n) / 4n) as F;
      return [r, r.neg() as F];
    }
    // Tonelli-Shanks: q - 1 = 2^s * m with m odd.
    let m = q - 1n;
    let s = 0n;
    while (m % 2n === 0n) {
      m /= 2n;
      s++;
    }
    // Find a non-residue.
    let z = K.__call__(2n) as F;
    let zi = 2n;
    while (this._is_square(z)) {
      zi++;
      z = K.__call__(zi) as F;
      if (zi > q) {
        throw new ArithmeticError('no quadratic non-residue found');
      }
    }
    let M = s;
    let c = z.pow(m) as F;
    let t = x.pow(m) as F;
    let r = x.pow((m + 1n) / 2n) as F;
    for (;;) {
      if (t.eq(K.one())) {
        return [r, r.neg() as F];
      }
      let i = 1n;
      let temp = t.mul(t) as F;
      while (!temp.eq(K.one())) {
        temp = temp.mul(temp) as F;
        i++;
      }
      const b = c.pow(1n << (M - i - 1n)) as F;
      M = i;
      c = b.mul(b) as F;
      t = t.mul(c) as F;
      r = r.mul(b) as F;
    }
  }

  /**
   * Return the scalar multiplication map [m].
   *
   * INPUT:
   * - m: an integer
   * - x_only: if True, return only the x-coordinate map (as phi_m / psi_m^2)
   *
   * OUTPUT: the multiplication-by-m map as rational functions
   *
   * For an elliptic curve E, the map [m]: E -> E is given by:
   *   x([m]P) = phi_m(x) / psi_m(x)^2
   *
   * where phi_m and psi_m are derived from the division polynomials.
   *
   * @example
   * ```typescript
   * const E = EllipticCurve(GF(23n), [1n, 1n]);
   * const [phiM, psiMSq] = E.multiplication_by_m(3n, true);
   * // x([3]P) = phiM(x(P)) / psiMSq(x(P))
   * ```
   *
   * @see Reference: sage/schemes/elliptic_curves/ell_generic.py:multiplication_by_m
   */
  multiplication_by_m<T extends RingElement>(
    m: bigint | number,
    x_only: boolean = false
  ): [Polynomial<T>, Polynomial<T>] | unknown {
    const mBig = typeof m === 'bigint' ? m : BigInt(m);

    if (mBig === 0n) {
      throw new ValueError('m must be a nonzero integer');
    }

    if (!x_only) {
      // Sage returns the pair of rational maps (mx, my); ``my`` needs the
      // bivariate rational function field, which the port does not have yet.
      throw new NotImplementedError(
        'multiplication_by_m with x_only=false requires bivariate polynomial support'
      );
    }

    const K = this.base_ring;
    const polyRing = new PolynomialRing(K as unknown as CoefficientRing<T>, 'x');
    const xPoly = polyRing.gen();

    // Special case of multiplication by +-1 is easy (ell_generic.py:2528-2538).
    if (mBig === 1n || mBig === -1n) {
      return [xPoly, polyRing.one()];
    }

    // The x-coordinate does not depend on the sign of m.
    const absM = mBig < 0n ? -mBig : mBig;
    return [this._multiple_x_numerator<T>(absM), this._multiple_x_denominator<T>(absM)];
  }

  /** Cache for _multiple_x_numerator (keyed by n, only when x is None). */
  private __mulxnums: Map<bigint, unknown> = new Map();

  /** Cache for _multiple_x_denominator (keyed by n, only when x is None). */
  private __mulxdens: Map<bigint, unknown> = new Map();

  /**
   * Return the numerator of the x-coordinate of the n-th multiple of a point,
   * using the division polynomials (without the 2-torsion factor).
   *
   * INPUT: n, x -- as described in division_polynomial_0
   *
   * If x is undefined, the result is cached.
   *
   * @see Reference: sage/schemes/elliptic_curves/ell_generic.py:_multiple_x_numerator
   */
  _multiple_x_numerator<T extends RingElement>(
    n: bigint | number,
    x?: T | Polynomial<T>
  ): Polynomial<T> {
    let nn = BigInt(n);
    if (nn < 0n) nn = -nn;
    if (nn === 0n) {
      throw new ValueError('n must be nonzero');
    }

    if (x === undefined) {
      const cached = this.__mulxnums.get(nn);
      if (cached !== undefined) {
        return cached as Polynomial<T>;
      }
    }

    const xx = this._x_as_polynomial<T>(x);

    if (nn === 1n) {
      if (x === undefined) {
        this.__mulxnums.set(nn, xx);
      }
      return xx;
    }

    const N = Number(nn);
    const polys = this.division_polynomial_0<T>([-2, -1, N - 1, N, N + 1], x) as Polynomial<T>[];

    let ret: Polynomial<T>;
    if (N % 2 === 0) {
      // xx * B_6 * psi_n^2 - psi_{n-1} * psi_{n+1}
      ret = xx.mul(polys[1]!).mul(polys[3]!.pow(2)).sub(polys[2]!.mul(polys[4]!));
    } else {
      // xx * psi_n^2 - B_6 * psi_{n-1} * psi_{n+1}
      ret = xx.mul(polys[3]!.pow(2)).sub(polys[1]!.mul(polys[2]!).mul(polys[4]!));
    }

    if (x === undefined) {
      this.__mulxnums.set(nn, ret);
    }
    return ret;
  }

  /**
   * Return the denominator of the x-coordinate of the n-th multiple of a
   * point, using the division polynomials.
   *
   * @see Reference: sage/schemes/elliptic_curves/ell_generic.py:_multiple_x_denominator
   */
  _multiple_x_denominator<T extends RingElement>(
    n: bigint | number,
    x?: T | Polynomial<T>
  ): Polynomial<T> {
    let nn = BigInt(n);
    if (nn < 0n) nn = -nn;
    if (nn === 0n) {
      throw new ValueError('n must be nonzero');
    }

    if (x === undefined) {
      const cached = this.__mulxdens.get(nn);
      if (cached !== undefined) {
        return cached as Polynomial<T>;
      }
    }

    const N = Number(nn);
    let ret = (this.division_polynomial_0<T>(N, x) as Polynomial<T>).pow(2);
    if (N % 2 === 0) {
      ret = ret.mul(this.division_polynomial_0<T>(-1, x) as Polynomial<T>);
    }

    if (x === undefined) {
      this.__mulxdens.set(nn, ret);
    }
    return ret;
  }

  /**
   * Build the polynomial standing for the "x variable" in the same ring that
   * division_polynomial_0 uses for the given ``x`` argument.
   */
  private _x_as_polynomial<T extends RingElement>(x?: T | Polynomial<T>): Polynomial<T> {
    if (x === undefined) {
      const polyRing = new PolynomialRing(this.base_ring as unknown as CoefficientRing<T>, 'x');
      return polyRing.gen() as unknown as Polynomial<T>;
    }
    if (x instanceof Polynomial) {
      return x;
    }
    const polyRing = new PolynomialRing(this.base_ring as unknown as CoefficientRing<T>, 'x');
    return polyRing.__call__(x) as unknown as Polynomial<T>;
  }

  /**
   * Return the scalar multiplication endomorphism [m].
   *
   * INPUT:
   * - m: an integer
   *
   * OUTPUT: the multiplication-by-m endomorphism as an object with:
   *   - domain(): the source curve
   *   - codomain(): the target curve (same as domain for endomorphism)
   *   - degree(): m^2 for the [m] map
   *   - call(P): compute [m]P
   *
   * @example
   * ```typescript
   * const E = EllipticCurve(GF(23n), [1n, 1n]);
   * const doubling = E.scalar_multiplication(2n);
   * const P = E.point([...]);
   * const twoP = doubling.call(P);  // Same as P.mul(2n)
   * ```
   *
   * @see Reference: sage/schemes/elliptic_curves/ell_generic.py:scalar_multiplication
   */
  scalar_multiplication(m: bigint | number): {
    domain: () => EllipticCurveGeneric<F>;
    codomain: () => EllipticCurveGeneric<F>;
    degree: () => bigint;
    call: (P: EllipticCurvePoint<F>) => EllipticCurvePoint<F>;
    toString: () => string;
  } {
    const mBig = typeof m === 'bigint' ? m : BigInt(m);

    return {
      domain: () => this,
      codomain: () => this,
      degree: () => mBig * mBig,
      call: (P: EllipticCurvePoint<F>) => P.mul(mBig),
      toString: () => `Scalar multiplication [${mBig}] endomorphism on ${this}`,
    };
  }

  /**
   * Return the n-th power Frobenius isogeny.
   *
   * For a curve E over F_q where q = p^k, the Frobenius endomorphism
   * pi: E -> E is defined by pi(x, y) = (x^p, y^p).
   *
   * INPUT:
   * - n: a positive integer (default 1)
   *
   * OUTPUT: the n-th power Frobenius isogeny as an object with:
   *   - domain(): the source curve
   *   - codomain(): the target curve (same curve for Frobenius)
   *   - degree(): p^n
   *   - call(P): compute pi^n(P) = (x^(p^n), y^(p^n))
   *
   * NOTE: For curves over prime fields F_p, the Frobenius is the identity
   * since x^p = x for all x in F_p.
   *
   * @example
   * ```typescript
   * const F = GF(23n * 23n);  // Extension field
   * const E = EllipticCurve(F, [1n, 1n]);
   * const frob = E.frobenius_isogeny();
   * const P = E.point([...]);
   * const frobP = frob.call(P);  // (x^23, y^23)
   * ```
   *
   * @see Reference: sage/schemes/elliptic_curves/ell_finite_field.py:frobenius_isogeny
   */
  frobenius_isogeny(n: number = 1): {
    domain: () => EllipticCurveGeneric<F>;
    codomain: () => EllipticCurveGeneric<F>;
    degree: () => bigint;
    call: (P: EllipticCurvePoint<F>) => EllipticCurvePoint<F>;
    is_separable: () => boolean;
    toString: () => string;
  } {
    const K = this.base_ring;
    const p = K.characteristic;

    if (p === 0n) {
      throw new ValueError('Frobenius isogeny is only defined for curves over finite fields');
    }
    const pn = p ** BigInt(n);

    return {
      domain: () => this,
      codomain: () => this,
      degree: () => pn,
      call: (P: EllipticCurvePoint<F>): EllipticCurvePoint<F> => {
        if (P.is_zero()) {
          return P;
        }
        // Frobenius: (x, y) -> (x^(p^n), y^(p^n))
        const x = P.x();
        const y = P.y();
        const xFrob = x.pow(pn) as F;
        const yFrob = y.pow(pn) as F;
        return this.point([xFrob, yFrob], false);
      },
      is_separable: () => false, // Frobenius is purely inseparable
      toString: () => `Frobenius isogeny pi^${n} on ${this}`,
    };
  }

  /**
   * Return the identity morphism on this curve.
   *
   * OUTPUT: the identity endomorphism as an object with:
   *   - domain(): the source curve
   *   - codomain(): the target curve (same as domain)
   *   - degree(): 1
   *   - call(P): return P
   *
   * @example
   * ```typescript
   * const E = EllipticCurve(GF(23n), [1n, 1n]);
   * const id = E.identity_morphism();
   * const P = E.point([...]);
   * const Q = id.call(P);  // Q equals P
   * ```
   *
   * @see Reference: sage/schemes/elliptic_curves/ell_generic.py:identity_morphism
   */
  identity_morphism(): {
    domain: () => EllipticCurveGeneric<F>;
    codomain: () => EllipticCurveGeneric<F>;
    degree: () => bigint;
    call: (P: EllipticCurvePoint<F>) => EllipticCurvePoint<F>;
    is_separable: () => boolean;
    is_injective: () => boolean;
    is_surjective: () => boolean;
    toString: () => string;
  } {
    return {
      domain: () => this,
      codomain: () => this,
      degree: () => 1n,
      call: (P: EllipticCurvePoint<F>) => P,
      is_separable: () => true,
      is_injective: () => true,
      is_surjective: () => true,
      toString: () => `Identity morphism on ${this}`,
    };
  }

  /**
   * Return an isomorphism from this curve to another curve.
   *
   * Given another Weierstrass model `other` of `self`, return an isomorphism
   * from `self` to `other` as a tuple [u, r, s, t] specifying the transformation:
   *
   *   (x, y) -> (u^2 * x' + r, u^3 * y' + s * u^2 * x' + t)
   *
   * INPUT:
   * - other: an elliptic curve isomorphic to self
   *
   * OUTPUT: An isomorphism [u, r, s, t] from self to other
   *
   * @throws {ValueError} If the curves are not isomorphic
   *
   * @example
   * ```typescript
   * const E = EllipticCurve(GF(23n), [1n, 1n]);
   * const F = E.short_weierstrass_model();
   * const iso = E.isomorphism_to(F);
   * // iso = [u, r, s, t] such that transforming E gives F
   * ```
   *
   * @see Reference: sage/schemes/elliptic_curves/ell_generic.py:isomorphism_to
   */
  isomorphism_to(other: EllipticCurveGeneric<F>): [F, F, F, F] {
    // Sage's ``isomorphism_to`` builds ``WeierstrassIsomorphism(self, None, other)``,
    // which takes ``next(_isomorphisms(E, F))`` -- the first tuple produced by the
    // *unsorted* generator (weierstrass_morphism.py:496-500).
    const isos = this._isomorphisms_unsorted(other);

    if (isos.length === 0) {
      throw new ValueError('elliptic curves not isomorphic');
    }

    return isos[0]!;
  }

  /**
   * Return the set of isomorphisms from this curve to itself (automorphisms).
   *
   * The identity and negation morphisms are guaranteed to appear as the first
   * and second entries of the returned list.
   *
   * OUTPUT: A list of [u, r, s, t] tuples representing automorphisms
   *
   * For a generic elliptic curve, the only automorphisms are:
   * - Identity: [1, 0, 0, 0]
   * - Negation: [-1, 0, 0, -a3] (maps P to -P)
   *
   * Curves with j = 0 have 6 automorphisms (over an algebraically closed field).
   * Curves with j = 1728 have 4 automorphisms (over an algebraically closed field).
   *
   * @example
   * ```typescript
   * const E = EllipticCurve(GF(23n), [1n, 1n]);
   * const auts = E.automorphisms();
   * // auts contains at least the identity and negation
   * ```
   *
   * @see Reference: sage/schemes/elliptic_curves/ell_generic.py:automorphisms
   */
  automorphisms(): Array<[F, F, F, F]> {
    return this._compute_isomorphisms(this);
  }

  /**
   * Return the set of isomorphisms from this curve to another (as a list).
   *
   * INPUT:
   * - other: another elliptic curve
   *
   * OUTPUT: A list of [u, r, s, t] tuples representing isomorphisms from
   *         self to other. Empty if the curves are not isomorphic.
   *
   * @example
   * ```typescript
   * const E = EllipticCurve(GF(23n), [1n, 1n]);
   * const F = E.short_weierstrass_model();
   * const isos = E.isomorphisms(F);
   * // isos contains all isomorphisms from E to F
   * ```
   *
   * @see Reference: sage/schemes/elliptic_curves/ell_generic.py:isomorphisms
   */
  isomorphisms(other: EllipticCurveGeneric<F>): Array<[F, F, F, F]> {
    return this._compute_isomorphisms(other);
  }

  /**
   * Internal method to compute all isomorphisms from self to other.
   *
   * An isomorphism is a tuple [u, r, s, t] with u != 0 such that the
   * transformation (x, y) -> ((x - r) / u^2, (y - s(x-r) - t) / u^3)
   * maps self to other.
   *
   * The inverse transformation is:
   *   (x', y') -> (u^2 * x' + r, u^3 * y' + s * u^2 * x' + t)
   */
  private _compute_isomorphisms(other: EllipticCurveGeneric<F>): Array<[F, F, F, F]> {
    const isos = this._isomorphisms_unsorted(other);

    // Sage returns ``sorted(...)`` of WeierstrassIsomorphism objects; the
    // ordering is given by ``WeierstrassIsomorphism._comparison_impl``'s
    // ``_sorting_key`` (weierstrass_morphism.py:568-574), which guarantees the
    // identity and the negation map come first.
    const a1 = this.a1();
    const a3 = this.a3();
    const one = this.base_ring.one() as F;
    const zero = this.base_ring.zero() as F;

    const negate = (v: [F, F, F, F]): [F, F, F, F] => {
      const [u, r, s, t] = v;
      return [u.neg() as F, r, s.neg().sub(a1) as F, t.neg().sub(a1.mul(r)).sub(a3) as F];
    };

    const isIdentity = (v: [F, F, F, F]): boolean =>
      v[0].eq(one) && v[1].isZero() && v[2].isZero() && v[3].isZero();

    const cmpTuple = (v: [F, F, F, F], w: [F, F, F, F]): number => {
      for (let i = 0; i < 4; i++) {
        const c = compareFieldElements(v[i]!, w[i]!);
        if (c !== 0) return c;
      }
      return 0;
    };

    const keyed = isos.map((v) => {
      const w = negate(v);
      const i = isIdentity(v) || isIdentity(w) ? 0 : 1;
      const j = v[0].eq(one) ? 0 : w[0].eq(one) ? 1 : 2;
      const mn = cmpTuple(v, w) <= 0 ? v : w;
      return { v, i, j, mn };
    });

    keyed.sort((A, B) => {
      if (A.i !== B.i) return A.i - B.i;
      const c = cmpTuple(A.mn, B.mn);
      if (c !== 0) return c;
      if (A.j !== B.j) return A.j - B.j;
      return cmpTuple(A.v, B.v);
    });

    void zero;
    return keyed.map((k) => k.v);
  }

  /**
   * Enumerate all isomorphisms (u, r, s, t) from this curve to ``other``.
   *
   * Direct port of ``sage.schemes.elliptic_curves.weierstrass_morphism._isomorphisms``.
   *
   * @see Reference: sage/schemes/elliptic_curves/weierstrass_morphism.py:_isomorphisms
   */
  private _isomorphisms_unsorted(other: EllipticCurveGeneric<F>): Array<[F, F, F, F]> {
    const K = this.base_ring;
    const result: Array<[F, F, F, F]> = [];

    const j = this.j_invariant();
    if (!j.eq(other.j_invariant())) {
      return result;
    }

    const [a1E, a2E, a3E, a4E, a6E] = this._ainvs;
    const [a1F, a2F, a3F, a4F, a6F] = other._ainvs;
    const char = K.characteristic;

    if (char === 2n) {
      if (j.isZero()) {
        // ulist = (x^3 - a3E/a3F).roots()
        const ulist = this._poly_roots([
          a3E.div(a3F).neg() as F,
          K.zero() as F,
          K.zero() as F,
          K.one() as F,
        ]);
        for (const u of ulist) {
          const u2 = u.mul(u) as F;
          const u4 = u2.mul(u2) as F;
          const u6 = u4.mul(u2) as F;
          // slist = (x^4 + a3E*x + (a2F^2 + a4F)*u^4 + a2E^2 + a4E).roots()
          const c0 = a2F.mul(a2F).add(a4F).mul(u4).add(a2E.mul(a2E)).add(a4E) as F;
          const slist = this._poly_roots([c0, a3E, K.zero() as F, K.zero() as F, K.one() as F]);
          for (const s of slist) {
            const r = s.mul(s).add(a2E).add(a2F.mul(u2)) as F;
            // tlist = (x^2 + a3E*x + r^3 + a2E*r^2 + a4E*r + a6E + a6F*u^6).roots()
            const d0 = r
              .mul(r)
              .mul(r)
              .add(a2E.mul(r).mul(r))
              .add(a4E.mul(r))
              .add(a6E)
              .add(a6F.mul(u6)) as F;
            const tlist = this._poly_roots([d0, a3E, K.one() as F]);
            for (const t of tlist) {
              result.push([u, r, s, t]);
            }
          }
        }
      } else {
        const u = a1E.div(a1F) as F;
        const u2 = u.mul(u) as F;
        const u3 = u2.mul(u) as F;
        const u4 = u2.mul(u2) as F;
        const r = a3E.add(a3F.mul(u3)).div(a1E) as F;
        // slist = (x^2 + a1E*x + r + a2E + a2F*u^2).roots()
        const slist = this._poly_roots([r.add(a2E).add(a2F.mul(u2)) as F, a1E, K.one() as F]);
        for (const s of slist) {
          const t = a4E
            .add(a4F.mul(u4))
            .add(s.mul(a3E))
            .add(r.mul(s).mul(a1E))
            .add(r.mul(r))
            .div(a1E) as F;
          result.push([u, r, s, t]);
        }
      }
      return result;
    }

    const [b2E, b4E, b6E] = this.b_invariants();
    const [b2F, b4F, b6F] = other.b_invariants();

    if (char === 3n) {
      if (j.isZero()) {
        // ulist = (x^4 - b4E/b4F).roots()
        const ulist = this._poly_roots([
          b4E.div(b4F).neg() as F,
          K.zero() as F,
          K.zero() as F,
          K.zero() as F,
          K.one() as F,
        ]);
        for (const u of ulist) {
          const u3 = u.mul(u).mul(u) as F;
          const u6 = u3.mul(u3) as F;
          const s = a1E.sub(a1F.mul(u)) as F;
          const t = a3E.sub(a3F.mul(u3)) as F;
          // rlist = (x^3 - b4E*x + b6E - b6F*u^6).roots()
          const rlist = this._poly_roots([
            b6E.sub(b6F.mul(u6)) as F,
            b4E.neg() as F,
            K.zero() as F,
            K.one() as F,
          ]);
          for (const r of rlist) {
            result.push([u, r, s, t.add(r.mul(a1E)) as F]);
          }
        }
      } else {
        // ulist = (x^2 - b2E/b2F).roots()
        const ulist = this._poly_roots([b2E.div(b2F).neg() as F, K.zero() as F, K.one() as F]);
        for (const u of ulist) {
          const u2 = u.mul(u) as F;
          const u3 = u2.mul(u) as F;
          const u4 = u2.mul(u2) as F;
          const r = b4F.mul(u4).sub(b4E).div(b2E) as F;
          const s = a1E.sub(a1F.mul(u)) as F;
          const t = a3E.sub(a3F.mul(u3)).add(a1E.mul(r)) as F;
          result.push([u, r, s, t]);
        }
      }
      return result;
    }

    // Now char != 2, 3.
    const [c4E, c6E] = this.c_invariants();
    const [c4F, c6F] = other.c_invariants();

    let m: number;
    let um: F;
    if (j.isZero()) {
      m = 6;
      um = c6E.div(c6F) as F;
    } else if (j.eq(K.__call__(1728n))) {
      m = 4;
      um = c4E.div(c4F) as F;
    } else {
      m = 2;
      um = c6E.mul(c4F).div(c6F.mul(c4E)) as F;
    }

    const coeffs: F[] = [um.neg() as F];
    for (let i = 1; i < m; i++) {
      coeffs.push(K.zero() as F);
    }
    coeffs.push(K.one() as F);

    const two = K.__call__(2n) as F;
    const three = K.__call__(3n) as F;

    for (const u of this._poly_roots(coeffs)) {
      const u2 = u.mul(u) as F;
      const u3 = u2.mul(u) as F;
      const s = a1F.mul(u).sub(a1E).div(two) as F;
      const r = a2F.mul(u2).add(a1E.mul(s)).add(s.mul(s)).sub(a2E).div(three) as F;
      const t = a3F.mul(u3).sub(a1E.mul(r)).sub(a3E).div(two) as F;
      result.push([u, r, s, t]);
    }

    return result;
  }

  /**
   * Return the roots in the base field of the polynomial whose coefficient
   * list (constant term first) is ``coeffs``.
   *
   * Mirrors Sage's use of ``(x**m - c).roots(multiplicities=False)``: a real
   * root finder, not a bounded brute-force search.
   */
  private _poly_roots(coeffs: F[]): F[] {
    const polyRing = new PolynomialRing(
      this.base_ring as unknown as CoefficientRing<RingElement>,
      'x'
    );
    const f = polyRing.__call__(coeffs as unknown as RingElement[]);
    // PARI's ``FpX_roots`` returns the roots sorted (see
    // reference/pari/src/basemath/factcyclo.c:491); sort so that the choice of
    // representative is deterministic and independent of the field size.
    const rts = f.roots().map(([r]: [RingElement, number]) => r as unknown as F);
    rts.sort((a, b) => compareFieldElements(a, b));
    return rts;
  }

  /**
   * Return True if this curve is isomorphic to another.
   *
   * Two elliptic curves over the same field are isomorphic if and only if
   * they have the same j-invariant (over an algebraically closed field).
   * Over non-algebraically closed fields, this may not be sufficient.
   *
   * INPUT:
   * - other: another elliptic curve
   *
   * OUTPUT: True if self is isomorphic to other over the base field
   *
   * @example
   * ```typescript
   * const E = EllipticCurve(GF(23n), [1n, 1n]);
   * const F = E.short_weierstrass_model();
   * console.log(E.is_isomorphic(F)); // true
   * ```
   *
   * @see Reference: sage/schemes/elliptic_curves/ell_generic.py:is_isomorphic
   */
  is_isomorphic(other: EllipticCurveGeneric<F>): boolean {
    // Quick check: j-invariants must match
    if (!this.j_invariant().eq(other.j_invariant())) {
      return false;
    }

    // Try to find an isomorphism
    const isos = this._compute_isomorphisms(other);
    return isos.length > 0;
  }

  /**
   * Return a pair of polynomials g(x), h(x) such that this elliptic curve
   * can be defined by the standard hyperelliptic equation y^2 + h(x)y = g(x).
   *
   * For the Weierstrass equation:
   *   y^2 + a1*x*y + a3*y = x^3 + a2*x^2 + a4*x + a6
   *
   * We have:
   *   g(x) = x^3 + a2*x^2 + a4*x + a6
   *   h(x) = a1*x + a3
   *
   * OUTPUT: a tuple [g, h] of polynomials
   *
   * @example
   * ```typescript
   * const E = EllipticCurve(GF(23n), [1n, 2n, 3n, 4n, 5n]);
   * const [g, h] = E.hyperelliptic_polynomials();
   * // g = x^3 + 2*x^2 + 4*x + 5
   * // h = x + 3
   * ```
   *
   * @see Reference: sage/schemes/elliptic_curves/ell_generic.py:hyperelliptic_polynomials
   */
  hyperelliptic_polynomials<T extends RingElement>(): [Polynomial<T>, Polynomial<T>] {
    const K = this.base_ring;
    const polyRing = new PolynomialRing(K as unknown as CoefficientRing<T>, 'x');

    const [a1, a2, a3, a4, a6] = this._ainvs;

    // g(x) = x^3 + a2*x^2 + a4*x + a6
    // Coefficients are [a6, a4, a2, 1] for constant, x, x^2, x^3
    const g = polyRing.__call__([a6, a4, a2, K.one()] as unknown as T[]);

    // h(x) = a1*x + a3
    // Coefficients are [a3, a1] for constant, x
    const h = polyRing.__call__([a3, a1] as unknown as T[]);

    return [g, h];
  }

  /**
   * Find a basis for the p-primary part of the torsion subgroup.
   *
   * The p-primary part consists of all points whose order is a power of p.
   * This function finds generators for this subgroup.
   *
   * INPUT:
   * - p: a prime number
   * - m: (optional) maximum power of p to consider
   *
   * OUTPUT: a list of [generator, order_exponent] pairs where order_exponent
   *         is the exponent e such that p^e is the order of the generator
   *
   * For finite fields, the p-primary part of E(F_q) is a subgroup of the form
   * Z/p^{e1} x Z/p^{e2} where e2 <= e1.
   *
   * @example
   * ```typescript
   * const E = EllipticCurve(GF(101n), [1n, 1n]);
   * const basis = E._p_primary_torsion_basis(2n);
   * // Returns generators for the 2-primary part
   * ```
   *
   * @see Reference: sage/schemes/elliptic_curves/ell_generic.py:_p_primary_torsion_basis
   */
  _p_primary_torsion_basis(
    p: bigint | number,
    m?: bigint | number
  ): Array<[EllipticCurvePoint<F>, number]> {
    return _p_primary_torsion_basis_impl(this, p, m);
  }
}
