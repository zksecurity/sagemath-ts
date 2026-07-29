/**
 * @module sage/algebras/quatalg/quaternion_algebra_element
 * @description Elements of quaternion algebras
 *
 * Port of: sage/algebras/quatalg/quaternion_algebra_element.pyx
 * Reference: reference/sage/src/sage/algebras/quatalg/quaternion_algebra_element.pyx
 *
 * Only the rational case (`QuaternionAlgebraElement_rational_field`) is ported;
 * see the module docstring of `quaternion_algebra.ts` for the scope of this port.
 */

import { NotImplementedError, ValueError, ZeroDivisionError } from '../../errors.js';
import { Matrix } from '../../matrix/matrix_generic.js';
import {
  determinant as _determinant,
  inverse as _inverse,
  solve_left as _solve_left,
} from '../../matrix/matrix_operations.js';
import { PolynomialRing as _PolynomialRing } from '../../rings/polynomial/polynomial_ring.js';
import { Rational } from '../../rings/rational.js';
import { QQ } from '../../rings/rational_field.js';
import type { QuaternionAlgebra_ab } from './quaternion_algebra.js';

/** Anything that can be turned into a rational coefficient of a quaternion. */
export type QuaternionCoefficient = Rational | bigint | number | string;

/**
 * Anything the element constructor accepts:  a scalar (constant quaternion) or
 * a length-4 list/tuple of coefficients.
 *
 * @see Reference: quaternion_algebra_element.pyx:141 (to_quaternion)
 */
export type QuaternionLike =
  | QuaternionAlgebraElement_rational_field
  | QuaternionCoefficient
  | readonly QuaternionCoefficient[];

function toRational(x: QuaternionCoefficient): Rational {
  if (x instanceof Rational) return x;
  return QQ.__call__(x as never);
}

/**
 * Internal function used implicitly by all quaternion algebra printing.
 *
 * @see Reference: quaternion_algebra_element.pyx:165 (print_coeff)
 */
function print_coeff(y: Rational, i: string): string {
  if (y.isZero()) {
    return '';
  }
  if (y.eq(1n)) {
    return i;
  }
  if (y.eq(-1n)) {
    return `-${i}`;
  }
  // The base ring is QQ, whose elements print atomically, so we never need
  // the parenthesised '(%s)*%s' branch of the reference implementation.
  return `${y.toString()}*${i}`;
}

/**
 * A matrix over `QQ`.
 *
 * This is a structurally typed view of `sage/matrix/matrix_generic.Matrix`
 * specialised to `Rational` entries: the repo's `Matrix<R extends RingElement>`
 * constraint is not satisfied by `Rational` (whose `add` accepts
 * `Rational | IntegerLike` rather than exactly `this`), so the generic class
 * cannot be instantiated with `Rational` under `strict` type checking.  The
 * runtime objects created here *are* `Matrix` instances, so all functions of
 * `sage/matrix/matrix_operations` apply to them unchanged.
 *
 * @see Deviation: `Matrix<Rational>` façade
 */
export interface RationalMatrix {
  readonly nrows: number;
  readonly ncols: number;
  get(i: number, j: number): Rational;
  set(i: number, j: number, value: Rational): void;
  row(i: number): Rational[];
  rows(): Rational[][];
  columns(): Rational[][];
  mul(other: RationalMatrix): RationalMatrix;
  add(other: RationalMatrix): RationalMatrix;
  sub(other: RationalMatrix): RationalMatrix;
  transpose(): RationalMatrix;
  is_square(): boolean;
  copy(): RationalMatrix;
  eq(other: RationalMatrix): boolean;
  toString(): string;
}

/**
 * A univariate polynomial over `QQ`, structurally typed (see
 * {@link RationalMatrix} for why the generic class cannot be used directly).
 *
 * @see Deviation: `Matrix<Rational>` façade
 */
export interface RationalPolynomial {
  readonly coeffs: readonly Rational[];
  degree(): number;
  toString(): string;
}

const PolynomialRingQQ = _PolynomialRing as unknown as new (
  ring: unknown,
  name?: string
) => { __call__(x: Rational[]): RationalPolynomial };

/** Construct a matrix over `QQ` (see {@link RationalMatrix}). */
export const MatrixQQ = Matrix as unknown as new (
  ring: unknown,
  nrows: number,
  ncols: number,
  entries?: Rational[][] | Rational[] | ((i: number, j: number) => Rational)
) => RationalMatrix;

/** `determinant` of `sage/matrix/matrix_operations`, typed for `QQ`. */
export const determinantQQ = _determinant as unknown as (M: RationalMatrix) => Rational;

/** `inverse` of `sage/matrix/matrix_operations`, typed for `QQ`. */
export const inverseQQ = _inverse as unknown as (M: RationalMatrix) => RationalMatrix;

/** `solve_left` of `sage/matrix/matrix_operations`, typed for `QQ`. */
export const solveLeftQQ = _solve_left as unknown as (
  A: RationalMatrix,
  B: RationalMatrix,
  check?: boolean
) => RationalMatrix;

/**
 * An element of a rational quaternion algebra.
 *
 * The element `x + y*i + z*j + w*k` is stored by its four rational
 * coefficients.  SageMath stores `(x, y, z, w, d)` as integers with a common
 * denominator `d`; that is an internal optimisation with identical semantics.
 *
 * @see Reference: quaternion_algebra_element.pyx:883
 *   (QuaternionAlgebraElement_rational_field)
 */
export class QuaternionAlgebraElement_rational_field {
  readonly _parent: QuaternionAlgebra_ab;
  private readonly _x: Rational;
  private readonly _y: Rational;
  private readonly _z: Rational;
  private readonly _w: Rational;

  /**
   * @param parent - the ambient quaternion algebra
   * @param v - a scalar, or a list of four coefficients
   *
   * @see Reference: quaternion_algebra_element.pyx:1028 (__init__)
   */
  constructor(parent: QuaternionAlgebra_ab, v: QuaternionLike) {
    this._parent = parent;
    if (v instanceof QuaternionAlgebraElement_rational_field) {
      this._x = v._x;
      this._y = v._y;
      this._z = v._z;
      this._w = v._w;
      return;
    }
    if (Array.isArray(v)) {
      const t = v as readonly QuaternionCoefficient[];
      if (t.length !== 4) {
        throw new ValueError('quaternion element must have exactly 4 coefficients');
      }
      this._x = toRational(t[0] as QuaternionCoefficient);
      this._y = toRational(t[1] as QuaternionCoefficient);
      this._z = toRational(t[2] as QuaternionCoefficient);
      this._w = toRational(t[3] as QuaternionCoefficient);
      return;
    }
    this._x = toRational(v as QuaternionCoefficient);
    this._y = Rational.zero();
    this._z = Rational.zero();
    this._w = Rational.zero();
  }

  /** The ambient quaternion algebra. */
  parent(): QuaternionAlgebra_ab {
    return this._parent;
  }

  /** The base ring, always `QQ` for this class. */
  base_ring(): typeof QQ {
    return QQ;
  }

  /**
   * Coefficient number `i` (0-based, in the basis `1, i, j, k`).
   *
   * @see Reference: quaternion_algebra_element.pyx:1106 (__getitem__)
   */
  get(i: number): Rational {
    switch (i) {
      case 0:
        return this._x;
      case 1:
        return this._y;
      case 2:
        return this._z;
      case 3:
        return this._w;
      default:
        throw new IndexError('quaternion element index out of range');
    }
  }

  /**
   * Return 4-tuple of rational numbers which are the coefficients of this
   * quaternion.
   *
   * @see Reference: quaternion_algebra_element.pyx:1568 (coefficient_tuple)
   */
  coefficient_tuple(): [Rational, Rational, Rational, Rational] {
    return [this._x, this._y, this._z, this._w];
  }

  /** Alias for {@link coefficient_tuple} matching Sage's `list(theta)`. */
  list(): Rational[] {
    return [this._x, this._y, this._z, this._w];
  }

  /**
   * Return the least common multiple of the denominators of the coefficients.
   *
   * @see Reference: quaternion_algebra_element.pyx:1491 (denominator)
   */
  denominator(): bigint {
    let d = 1n;
    for (const c of this.list()) {
      const dc = c.denominator;
      d = (d / gcdBigInt(d, dc)) * dc;
    }
    return d;
  }

  /**
   * Return `[d, x, y, z, w]` with `self = (x + y*i + z*j + w*k)/d` and
   * `gcd(d, x, y, z, w) = 1`.
   *
   * @see Reference: quaternion_algebra_element.pyx:1516
   *   (denominator_and_integer_coefficient_tuple)
   */
  denominator_and_integer_coefficient_tuple(): [bigint, bigint, bigint, bigint, bigint] {
    const d = this.denominator();
    const t = this.list().map((c) => c.numerator * (d / c.denominator));
    return [d, t[0] as bigint, t[1] as bigint, t[2] as bigint, t[3] as bigint];
  }

  /**
   * Return the integer part of this quaternion, ignoring the common
   * denominator.
   *
   * @see Reference: quaternion_algebra_element.pyx:1544
   *   (integer_coefficient_tuple)
   */
  integer_coefficient_tuple(): [bigint, bigint, bigint, bigint] {
    const [, x, y, z, w] = this.denominator_and_integer_coefficient_tuple();
    return [x, y, z, w];
  }

  /**
   * Return ``True`` if this quaternion is constant, i.e., has no `i`, `j`,
   * or `k` term.
   *
   * @see Reference: quaternion_algebra_element.pyx:223 (is_constant)
   */
  is_constant(): boolean {
    return this._y.isZero() && this._z.isZero() && this._w.isZero();
  }

  /** Whether this quaternion is zero (Sage: `not bool(x)`). */
  is_zero(): boolean {
    return this._x.isZero() && this._y.isZero() && this._z.isZero() && this._w.isZero();
  }

  /** Alias of {@link is_zero} for the repo's `RingElement` interface. */
  isZero(): boolean {
    return this.is_zero();
  }

  /**
   * Equality of quaternions.
   *
   * @see Reference: quaternion_algebra_element.pyx:989 (_richcmp_)
   */
  eq(other: QuaternionAlgebraElement_rational_field | QuaternionCoefficient): boolean {
    const o =
      other instanceof QuaternionAlgebraElement_rational_field
        ? other
        : new QuaternionAlgebraElement_rational_field(this._parent, other);
    for (let i = 0; i < 4; i++) {
      if (!this.get(i).eq(o.get(i))) return false;
    }
    return true;
  }

  /** Sum of two quaternions. @see Reference: quaternion_algebra_element.pyx:1149 */
  add(other: QuaternionAlgebraElement_rational_field): QuaternionAlgebraElement_rational_field {
    return new QuaternionAlgebraElement_rational_field(this._parent, [
      this._x.add(other._x),
      this._y.add(other._y),
      this._z.add(other._z),
      this._w.add(other._w),
    ]);
  }

  /** Difference of two quaternions. @see Reference: quaternion_algebra_element.pyx:1205 */
  sub(other: QuaternionAlgebraElement_rational_field): QuaternionAlgebraElement_rational_field {
    return new QuaternionAlgebraElement_rational_field(this._parent, [
      this._x.sub(other._x),
      this._y.sub(other._y),
      this._z.sub(other._z),
      this._w.sub(other._w),
    ]);
  }

  /** Negation. */
  neg(): QuaternionAlgebraElement_rational_field {
    return new QuaternionAlgebraElement_rational_field(this._parent, [
      this._x.neg(),
      this._y.neg(),
      this._z.neg(),
      this._w.neg(),
    ]);
  }

  /**
   * Product of two quaternions.
   *
   * Uses the multiplication formulas of the reference implementation:
   * ```
   * x = x1*x2 + y1*y2*a + z1*z2*b - w1*w2*a*b
   * y = x1*y2 + y1*x2 - z1*w2*b + w1*z2*b
   * z = x1*z2 + y1*w2 + z1*x2 - w1*y2*a
   * w = x1*w2 + y1*z2 - z1*y2 + w1*x2
   * ```
   *
   * @see Reference: quaternion_algebra_element.pyx:826 (_mul_)
   */
  mul(other: QuaternionAlgebraElement_rational_field): QuaternionAlgebraElement_rational_field {
    const a = this._parent._a;
    const b = this._parent._b;
    const [x1, y1, z1, w1] = this.coefficient_tuple();
    const [x2, y2, z2, w2] = other.coefficient_tuple();

    const x = x1
      .mul(x2)
      .add(y1.mul(y2).mul(a))
      .add(z1.mul(z2).mul(b))
      .sub(w1.mul(w2).mul(a).mul(b));
    const y = x1.mul(y2).add(y1.mul(x2)).sub(z1.mul(w2).mul(b)).add(w1.mul(z2).mul(b));
    const z = x1.mul(z2).add(y1.mul(w2).mul(a)).add(z1.mul(x2)).sub(w1.mul(y2).mul(a));
    const w = x1.mul(w2).add(y1.mul(z2)).sub(z1.mul(y2)).add(w1.mul(x2));

    return new QuaternionAlgebraElement_rational_field(this._parent, [x, y, z, w]);
  }

  /**
   * Multiply by a rational scalar.
   *
   * @see Reference: quaternion_algebra_element.pyx:511 (_rmul_), :526 (_lmul_)
   */
  scalar_mul(c: QuaternionCoefficient): QuaternionAlgebraElement_rational_field {
    const r = toRational(c);
    return new QuaternionAlgebraElement_rational_field(this._parent, [
      this._x.mul(r),
      this._y.mul(r),
      this._z.mul(r),
      this._w.mul(r),
    ]);
  }

  /**
   * Return the conjugate `x - y*i - z*j - w*k` of `x + y*i + z*j + w*k`.
   *
   * @see Reference: quaternion_algebra_element.pyx:415 (conjugate)
   */
  conjugate(): QuaternionAlgebraElement_rational_field {
    return new QuaternionAlgebraElement_rational_field(this._parent, [
      this._x,
      this._y.neg(),
      this._z.neg(),
      this._w.neg(),
    ]);
  }

  /**
   * Return the reduced trace `2*x`.
   *
   * @see Reference: quaternion_algebra_element.pyx:439 (reduced_trace)
   */
  reduced_trace(): Rational {
    return this._x.mul(new Rational(2n));
  }

  /**
   * Return the reduced norm `x^2 - a*y^2 - b*z^2 + a*b*w^2`.
   *
   * @see Reference: quaternion_algebra_element.pyx:454 (reduced_norm)
   */
  reduced_norm(): Rational {
    const a = this._parent._a;
    const b = this._parent._b;
    const [x, y, z, w] = this.coefficient_tuple();
    return w.mul(w).mul(a).mul(b).sub(y.mul(y).mul(a)).sub(z.mul(z).mul(b)).add(x.mul(x));
  }

  /**
   * Return the inverse `~self = (1/N(self)) * conjugate(self)`.
   *
   * @see Reference: quaternion_algebra_element.pyx:472 (__invert__)
   */
  inverse(): QuaternionAlgebraElement_rational_field {
    const n = this.reduced_norm();
    if (n.isZero()) {
      throw new ZeroDivisionError('rational division by zero');
    }
    return this.conjugate().scalar_mul(n.inv());
  }

  /** Division: `self * ~right`. @see Reference: quaternion_algebra_element.pyx:541 (_div_) */
  div(right: QuaternionAlgebraElement_rational_field): QuaternionAlgebraElement_rational_field {
    return this.mul(right.inverse());
  }

  /**
   * Return the reduced characteristic polynomial `X^2 - t*X + n`.
   *
   * @see Reference: quaternion_algebra_element.pyx:558
   *   (reduced_characteristic_polynomial)
   */
  reduced_characteristic_polynomial(varName: string = 'x'): RationalPolynomial {
    const R = new PolynomialRingQQ(QQ, varName);
    return R.__call__([this.reduced_norm(), this.reduced_trace().neg(), new Rational(1n)]);
  }

  /**
   * Return the matrix of right (default) or left multiplication by this
   * element on the basis `1, i, j, k`.
   *
   * @see Reference: quaternion_algebra_element.pyx:586 (matrix)
   */
  matrix(action: 'left' | 'right' = 'right'): RationalMatrix {
    let v: Rational[][];
    if (action === 'right') {
      v = this._parent.basis().map((a) => a.mul(this).coefficient_tuple());
    } else if (action === 'left') {
      v = this._parent.basis().map((a) => this.mul(a).coefficient_tuple());
    } else {
      throw new ValueError("action must be either 'left' or 'right'");
    }
    return new MatrixQQ(QQ, 4, 4, v);
  }

  /**
   * Return `(self.conjugate() * right).reduced_trace()`.
   *
   * @see Reference: quaternion_algebra_element.pyx:660 (pair)
   */
  pair(right: QuaternionAlgebraElement_rational_field): Rational {
    return this.conjugate().mul(right).reduced_trace();
  }

  /**
   * Return the image of this quaternion under the morphism defined by
   * `im_gens`, where elements of the base ring are mapped by `base_map`.
   *
   * @see Reference: quaternion_algebra_element.pyx:685 (_im_gens_)
   */
  _im_gens_(
    codomain: QuaternionAlgebra_ab,
    im_gens: readonly QuaternionAlgebraElement_rational_field[],
    base_map?: ((c: Rational) => Rational) | null
  ): QuaternionAlgebraElement_rational_field {
    const f = base_map ?? ((c: Rational): Rational => c);
    const gens = [codomain.one(), ...im_gens];
    let s = codomain.zero();
    for (let i = 0; i < 4; i++) {
      s = s.add((gens[i] as QuaternionAlgebraElement_rational_field).scalar_mul(f(this.get(i))));
    }
    return s;
  }

  /**
   * Multiply this quaternion by an integer.
   *
   * @see Reference: quaternion_algebra_element.pyx:1599 (_multiply_by_integer)
   */
  _multiply_by_integer(n: bigint): QuaternionAlgebraElement_rational_field {
    return this.scalar_mul(n);
  }

  /**
   * Divide this quaternion by an integer.
   *
   * @see Reference: quaternion_algebra_element.pyx:1636 (_divide_by_integer)
   */
  _divide_by_integer(n: bigint): QuaternionAlgebraElement_rational_field {
    if (n === 0n) {
      throw new ZeroDivisionError('division by zero');
    }
    return this.scalar_mul(new Rational(1n, n));
  }

  /**
   * String representation, e.g. `-3/4 + i + j + k`.
   *
   * @see Reference: quaternion_algebra_element.pyx:338 (_do_print)
   */
  toString(): string {
    const [x, y, z, w] = this.coefficient_tuple();
    const [i, j, k] = this._parent.variable_names();
    const v: string[] = [];
    if (!x.isZero()) {
      v.push(x.toString());
    }
    for (const [coeff, name] of [
      [y, i],
      [z, j],
      [w, k],
    ] as [Rational, string][]) {
      const c = print_coeff(coeff, name);
      if (c) v.push(c);
    }
    if (v.length === 0) {
      return '0';
    }
    return v.join(' + ').replace(/\+ -/g, '- ');
  }
}

/** Alias: over `QQ` this is the only element class. */
export type QuaternionAlgebraElement = QuaternionAlgebraElement_rational_field;

/**
 * The generic and number-field element classes of SageMath are not ported.
 *
 * @see Reference: quaternion_algebra_element.pyx:710
 *   (QuaternionAlgebraElement_generic), :1671
 *   (QuaternionAlgebraElement_number_field)
 */
export function QuaternionAlgebraElement_generic(): never {
  throw new NotImplementedError(
    'SAGE_NOT_IMPLEMENTED: QuaternionAlgebraElement_generic (quaternion algebras over base rings other than QQ)'
  );
}

/** @see {@link QuaternionAlgebraElement_generic} */
export function QuaternionAlgebraElement_number_field(): never {
  throw new NotImplementedError(
    'SAGE_NOT_IMPLEMENTED: QuaternionAlgebraElement_number_field (quaternion algebras over number fields)'
  );
}

function gcdBigInt(a: bigint, b: bigint): bigint {
  let x = a < 0n ? -a : a;
  let y = b < 0n ? -b : b;
  while (y) {
    [x, y] = [y, x % y];
  }
  return x;
}

/** Thrown by {@link QuaternionAlgebraElement_rational_field.get} (Sage: IndexError). */
export class IndexError extends ValueError {
  override name = 'IndexError';
}
