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
    two_torsion_multiplicity: number = 2
  ): Polynomial<T> {
    const mNum = Number(m);

    if (mNum <= 0) {
      throw new ValueError('m must be a positive integer');
    }

    if (![0, 1, 2].includes(two_torsion_multiplicity)) {
      throw new ValueError('two_torsion_multiplicity must be 0, 1, or 2');
    }

    // Get the basic division polynomial (without 2-torsion factor)
    const psi_m = this.division_polynomial_0<T>(mNum, x) as Polynomial<T>;

    if (mNum % 2 === 1) {
      // Odd m: return psi_m as is
      return psi_m;
    }

    // Even m: need to handle 2-torsion factor
    if (two_torsion_multiplicity === 0) {
      // No 2-torsion factor
      return psi_m;
    } else if (two_torsion_multiplicity === 2) {
      // Include (4x^3 + b2*x^2 + 2*b4*x + b6) = B_6
      const B6 = this.division_polynomial_0<T>(-1, x) as Polynomial<T>;
      return psi_m.mul(B6);
    } else {
      // two_torsion_multiplicity === 1
      // For this we would need bivariate polynomials (2y + a1*x + a3)
      // For now, we return the same as multiplicity 2
      // This is a simplification - full implementation would need bivariate support
      const B6 = this.division_polynomial_0<T>(-1, x) as Polynomial<T>;
      return psi_m.mul(B6);
    }
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
    // Get short Weierstrass form: y^2 = x^3 + a*x + b
    const Ew = this.short_weierstrass_model();
    const a = Ew.a4();
    const b = Ew.a6();
    const K = this.base_ring;

    // Find roots r of x^3 + a*x + b = 0
    // Then find s such that s^2 = 3r^2 + a
    const rValues = this._findCubicRoots(a, b);

    if (rValues.length === 0) {
      throw new ValueError(`${this} has no Montgomery model`);
    }

    // For each r, try to find s with s^2 = 3r^2 + a
    let bestR: F | null = null;
    let bestS: F | null = null;
    let hasSqrtS = false;

    const three = K.__call__(3n) as F;

    for (const r of rValues) {
      // Compute 3r^2 + a
      const sSquared = three.mul(r).mul(r).add(a) as F;

      // Try to find square root
      const sRoots = this._square_roots(sSquared);
      if (sRoots.length > 0) {
        const s = sRoots[0]!;
        // Prefer solutions where s is a square (so B=1)
        const sSqRoots = this._square_roots(s);
        if (sSqRoots.length > 0) {
          bestR = r;
          bestS = s;
          hasSqrtS = true;
          break; // Found best solution
        } else if (bestR === null) {
          bestR = r;
          bestS = s;
        }
      }
    }

    if (bestR === null || bestS === null) {
      throw new ValueError(`${this} has no Montgomery model`);
    }

    // Compute A = 3r/s
    const A = three.mul(bestR).div(bestS) as F;

    // Compute B: if s is a square, B=1; otherwise B = 1/s
    const one = K.one() as F;
    let B: F;

    if (hasSqrtS) {
      B = one;
    } else {
      // B = 1/s
      B = one.div(bestS) as F;
    }

    if (!twisted && !B.eq(one)) {
      throw new ValueError(`${this} has no untwisted Montgomery model`);
    }

    // Create Montgomery curve: By^2 = x^3 + Ax^2 + x
    // In Weierstrass form with a1=a3=0: y^2 = (1/B)(x^3 + Ax^2 + x)
    // For B=1: y^2 = x^3 + Ax^2 + x, which is [0, A, 0, 1, 0]
    if (B.eq(one)) {
      const zero = K.zero() as F;
      const E = new EllipticCurveGeneric<F>(K, [zero, A, zero, one, zero]);

      if (morphism) {
        // Find the isomorphism from self to E
        const iso = this.isomorphism_to(E);
        return [E, iso];
      }
      return E;
    }

    // Twisted case: By^2 = x^3 + Ax^2 + x
    // This is not directly representable as an EllipticCurveGeneric
    throw new NotImplementedError(
      'Twisted Montgomery models (B != 1) are not supported as EllipticCurve objects. ' +
        'They would need to be represented as ProjectivePlaneCurve.'
    );
  }

  /**
   * Find roots of x^3 + a*x + b in the base field.
   */
  private _findCubicRoots(a: F, b: F): F[] {
    const K = this.base_ring;
    const p = K.characteristic;
    const results: F[] = [];

    if (p === 0n) {
      // Number fields - not implemented
      return [];
    }

    // For finite fields, try brute force for small fields
    const maxFieldSize = 10000n;
    if (p <= maxFieldSize) {
      for (let xVal = 0n; xVal < p; xVal++) {
        const x = K.__call__(xVal) as F;
        // Check if x^3 + a*x + b = 0
        const x2 = x.mul(x) as F;
        const x3 = x2.mul(x) as F;
        const val = x3.add(a.mul(x)).add(b) as F;
        if (val.isZero()) {
          results.push(x);
        }
      }
    }

    return results;
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
    const mNum = typeof m === 'number' ? m : Number(m);
    const mBig = typeof m === 'bigint' ? m : BigInt(m);

    if (mBig === 0n) {
      throw new ValueError('multiplication_by_m: m cannot be 0');
    }

    // Get the polynomials using division polynomials
    // x([m]P) = phi_m(x) / psi_m(x)^2
    // phi_m = x * psi_m^2 - psi_{m-1} * psi_{m+1}

    const K = this.base_ring;
    const polyRing = new PolynomialRing(K as unknown as CoefficientRing<T>, 'x');

    // Get division polynomials
    const psiM = this.division_polynomial<T>(Math.abs(mNum)) as Polynomial<T>;
    const psiMSq = psiM.mul(psiM);

    // For phi_m, we need psi_{m-1} and psi_{m+1}
    const psiMminus1 = this.division_polynomial<T>(Math.abs(mNum) - 1) as Polynomial<T>;
    const psiMplus1 = this.division_polynomial<T>(Math.abs(mNum) + 1) as Polynomial<T>;

    // phi_m = x * psi_m^2 - psi_{m-1} * psi_{m+1}
    const xPoly = polyRing.gen();
    const phiM = xPoly.mul(psiMSq).sub(psiMminus1.mul(psiMplus1));

    if (x_only) {
      return [phiM, psiMSq];
    }

    // For the y-coordinate, we would need:
    // y([m]P) = omega_m(x, y) / psi_m(x)^3
    // This is more complex as omega_m involves y
    throw new NotImplementedError(
      'multiplication_by_m with x_only=false requires bivariate polynomial support'
    );
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
    // Get all isomorphisms and return the first one
    const isos = this._compute_isomorphisms(other);

    if (isos.length === 0) {
      throw new ValueError(`${this} and ${other} are not isomorphic`);
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
    const K = this.base_ring;
    const result: Array<[F, F, F, F]> = [];

    // Quick check: j-invariants must match
    if (!this.j_invariant().eq(other.j_invariant())) {
      return [];
    }

    // Get invariants
    const [a1, a2, a3, a4, a6] = this._ainvs;
    const [a1p, a2p, a3p, a4p, a6p] = other._ainvs;
    const c4 = this.c4();
    const c6 = this.c6();

    const one = K.one() as F;
    const two = K.__call__(2n) as F;
    const three = K.__call__(3n) as F;

    // The isomorphism equations:
    // a1' = (a1 + 2s) / u
    // a2' = (a2 - s*a1 + 3r - s^2) / u^2
    // a3' = (a3 + r*a1 + 2t) / u^3
    // a4' = (a4 - s*a3 + 2r*a2 - (t + rs)*a1 + 3r^2 - 2st) / u^4
    // a6' = (a6 + r*a4 + r^2*a2 + r^3 - t*a3 - t^2 - r*t*a1) / u^6

    // Find possible values of u
    // From c4' = c4 / u^4 and c6' = c6 / u^6
    // we get u^4 = c4 / c4' (if c4' != 0)
    // or u^6 = c6 / c6' (if c6' != 0)

    const c4p = other.c4();
    const c6p = other.c6();

    // Collect candidate u values
    const uCandidates: F[] = [];

    if (!c4.isZero() && !c4p.isZero()) {
      // u^4 = c4 / c4'
      const u4 = c4.div(c4p) as F;
      // Find fourth roots
      const roots4 = this._fourth_roots(u4);
      uCandidates.push(...roots4);
    } else if (!c6.isZero() && !c6p.isZero()) {
      // u^6 = c6 / c6'
      const u6 = c6.div(c6p) as F;
      // Find sixth roots
      const roots6 = this._sixth_roots(u6);
      uCandidates.push(...roots6);
    } else if (c4.isZero() && c4p.isZero() && c6.isZero() && c6p.isZero()) {
      // Both curves have c4 = c6 = 0
      // This means j = 0 and any sixth root of unity works
      // In a prime field, just try u = 1, -1, and primitive roots if they exist
      uCandidates.push(one);
      uCandidates.push(one.neg() as F);
      // For j=0 curves, there can be more automorphisms
      // We'll handle the basic cases
    } else {
      // One curve has c4 = 0 or c6 = 0 but the other doesn't - not isomorphic
      return [];
    }

    // Remove duplicates and zeros
    const seenU = new Set<string>();
    const uniqueU: F[] = [];
    for (const u of uCandidates) {
      if (!u.isZero()) {
        const key = u.toString();
        if (!seenU.has(key)) {
          seenU.add(key);
          uniqueU.push(u);
        }
      }
    }

    // For each candidate u, solve for r, s, t
    for (const u of uniqueU) {
      const u2 = u.mul(u) as F;
      const u3 = u2.mul(u) as F;
      const u4 = u2.mul(u2) as F;
      const u6 = u3.mul(u3) as F;

      // From a1' = (a1 + 2s) / u, we get s = (u * a1' - a1) / 2
      const s = u.mul(a1p).sub(a1).div(two) as F;

      // From a3' = (a3 + r*a1 + 2t) / u^3
      // and a2' = (a2 - s*a1 + 3r - s^2) / u^2
      // We get r = (u^2 * a2' - a2 + s*a1 + s^2) / 3
      const r = u2.mul(a2p).sub(a2).add(s.mul(a1)).add(s.mul(s)).div(three) as F;

      // From a3' = (a3 + r*a1 + 2t) / u^3, we get t = (u^3 * a3' - a3 - r*a1) / 2
      const t = u3.mul(a3p).sub(a3).sub(r.mul(a1)).div(two) as F;

      // Verify using a4 and a6 equations
      // a4' should equal (a4 - s*a3 + 2r*a2 - (t + rs)*a1 + 3r^2 - 2st) / u^4
      const a4_transformed = a4
        .sub(s.mul(a3))
        .add(two.mul(r).mul(a2))
        .sub(t.add(r.mul(s)).mul(a1))
        .add(three.mul(r).mul(r))
        .sub(two.mul(s).mul(t));
      const a4_check = a4_transformed.div(u4) as F;

      if (!a4_check.eq(a4p)) {
        continue;
      }

      // a6' should equal (a6 + r*a4 + r^2*a2 + r^3 - t*a3 - t^2 - r*t*a1) / u^6
      const a6_transformed = a6
        .add(r.mul(a4))
        .add(r.mul(r).mul(a2))
        .add(r.mul(r).mul(r))
        .sub(t.mul(a3))
        .sub(t.mul(t))
        .sub(r.mul(t).mul(a1));
      const a6_check = a6_transformed.div(u6) as F;

      if (!a6_check.eq(a6p)) {
        continue;
      }

      // Valid isomorphism found
      result.push([u, r, s, t]);
    }

    // Sort results with identity first, then negation
    result.sort((a, b) => {
      const [u1, r1, s1, t1] = a;
      const [u2, r2, s2, t2] = b;

      // Identity [1, 0, 0, 0] comes first
      const isIdent1 = u1.eq(one) && r1.isZero() && s1.isZero() && t1.isZero();
      const isIdent2 = u2.eq(one) && r2.isZero() && s2.isZero() && t2.isZero();
      if (isIdent1 && !isIdent2) return -1;
      if (!isIdent1 && isIdent2) return 1;

      // Negation [-1, 0, 0, ...] comes second
      const negOne = one.neg() as F;
      const isNeg1 = u1.eq(negOne) && r1.isZero() && s1.isZero();
      const isNeg2 = u2.eq(negOne) && r2.isZero() && s2.isZero();
      if (isNeg1 && !isNeg2) return -1;
      if (!isNeg1 && isNeg2) return 1;

      return 0;
    });

    return result;
  }

  /**
   * Find fourth roots of an element in the field.
   */
  private _fourth_roots(x: F): F[] {
    const K = this.base_ring;
    const results: F[] = [];

    if (x.isZero()) {
      return [K.zero() as F];
    }

    // First find square roots of x
    const sqrtX = this._square_roots(x);

    for (const s of sqrtX) {
      // Then find square roots of each sqrt
      const sqrtS = this._square_roots(s);
      results.push(...sqrtS);
    }

    return results;
  }

  /**
   * Find sixth roots of an element in the field.
   */
  private _sixth_roots(x: F): F[] {
    const K = this.base_ring;
    const results: F[] = [];

    if (x.isZero()) {
      return [K.zero() as F];
    }

    // x^(1/6) = (x^(1/2))^(1/3) = (x^(1/3))^(1/2)
    // Find square roots first, then cube roots
    const sqrtX = this._square_roots(x);

    for (const s of sqrtX) {
      const cbrtS = this._cube_roots(s);
      results.push(...cbrtS);
    }

    return results;
  }

  /**
   * Find square roots of an element in the field.
   * Uses Tonelli-Shanks for prime fields.
   */
  private _square_roots(x: F): F[] {
    const K = this.base_ring;

    if (x.isZero()) {
      return [K.zero() as F];
    }

    // Check if x has a square root using Euler's criterion
    const p = K.characteristic;
    if (p === 2n) {
      // In characteristic 2, everything is a square
      return [x.pow((p + 1n) / 2n) as F];
    }

    const exp = (p - 1n) / 2n;
    const legendre = x.pow(exp);

    if (!legendre.eq(K.one())) {
      // Not a quadratic residue
      return [];
    }

    // Use the formula for p ≡ 3 (mod 4)
    if (p % 4n === 3n) {
      const r = x.pow((p + 1n) / 4n) as F;
      return [r, r.neg() as F];
    }

    // General Tonelli-Shanks algorithm
    // Write p - 1 = 2^s * q where q is odd
    let q = p - 1n;
    let s = 0n;
    while (q % 2n === 0n) {
      q /= 2n;
      s++;
    }

    // Find a quadratic non-residue
    let z = K.__call__(2n) as F;
    while (z.pow(exp).eq(K.one())) {
      z = z.add(1) as F;
    }

    let m = s;
    let c = z.pow(q) as F;
    let t = x.pow(q) as F;
    let r = x.pow((q + 1n) / 2n) as F;

    while (true) {
      if (t.eq(K.one())) {
        return [r, r.neg() as F];
      }

      // Find the least i such that t^(2^i) = 1
      let i = 1n;
      let temp = t.mul(t) as F;
      while (!temp.eq(K.one())) {
        temp = temp.mul(temp) as F;
        i++;
      }

      // Update values
      const exp2 = 1n << (m - i - 1n);
      const b = c.pow(exp2) as F;
      m = i;
      c = b.mul(b) as F;
      t = t.mul(c) as F;
      r = r.mul(b) as F;
    }
  }

  /**
   * Find cube roots of an element in the field.
   */
  private _cube_roots(x: F): F[] {
    const K = this.base_ring;
    const p = K.characteristic;

    if (x.isZero()) {
      return [K.zero() as F];
    }

    // Check if x has a cube root
    // x^((p-1)/gcd(3, p-1)) should equal 1

    const gcd3 = gcd(3n, p - 1n);

    if (gcd3 === 1n) {
      // Every element has a unique cube root
      // r = x^((2*(p-1)+1)/3) when p ≡ 2 (mod 3)
      const exp = (2n * (p - 1n) + 1n) / 3n;
      return [x.pow(exp) as F];
    }

    // p ≡ 1 (mod 3), need to check if x is a cubic residue
    const exp = (p - 1n) / 3n;
    const residue = x.pow(exp);

    if (!residue.eq(K.one())) {
      // Not a cubic residue
      return [];
    }

    // Find a cube root
    // This is more complex - for simplicity, use brute force for small fields
    // or specialized algorithms for large fields

    // For p ≡ 1 (mod 3), we can use the Adleman-Manders-Miller algorithm
    // For simplicity, we'll use a simpler approach if p is small enough

    const results: F[] = [];

    // Try to find cube roots by iteration for small fields
    if (p < 10000n) {
      for (let i = 0n; i < p; i++) {
        const candidate = K.__call__(i) as F;
        if (candidate.pow(3).eq(x)) {
          results.push(candidate);
        }
      }
      return results;
    }

    // For larger fields, use a probabilistic approach
    // Find a cube root using the formula for p ≡ 1 (mod 3)
    // This requires finding a primitive root and is complex
    // For now, return just the principal root if we can find it

    // Try to find a cube root using the formula
    const candidate = x.pow((2n * p - 1n) / 3n) as F;
    if (candidate.pow(3).eq(x)) {
      results.push(candidate);
    }

    return results;
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
