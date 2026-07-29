/**
 * @module sage/quadratic_forms/quadratic_form
 * @description Quadratic forms in `n` variables over `ZZ` or `QQ`.
 *
 * Port of: `sage/quadratic_forms/quadratic_form.py`
 * Reference: `reference/sage/src/sage/quadratic_forms/quadratic_form.py`
 *
 * SageMath assembles {@link QuadraticForm} from a large number of
 * `quadratic_form__*.py` files which are imported *into the class body*.  We
 * mirror that layout: the core class lives here, the `p`-adic/real invariants
 * live in `quadratic_form__local_field_invariants.ts`, and the ternary
 * specialisation lives in `ternary_qf.ts`.
 *
 * Scope of this port (see the module notes at the bottom of the file for the
 * complete list of what is *not* implemented and throws
 * `NotImplementedError`):
 *
 * - base ring `ZZ` or `QQ` only (SageMath allows an arbitrary commutative ring)
 * - Gram / Hessian matrices, determinants, level, primitivity, evaluation
 * - direct sums, coefficient sums, change of ring, bilinear map
 * - rational diagonalisation, signature, Hasse invariants, (an)isotropy,
 *   definiteness (in `quadratic_form__local_field_invariants.ts`)
 * - rational isometry testing (Hasse--Minkowski)
 * - theta series via PARI's `qfrep` (delegated to `parigp-ts`)
 */

import { qfrep } from '@sagemath-ts/parigp-ts';
import { gcd as GCD, lcm as LCM } from '../arith/misc.js';
import { NotImplementedError, RuntimeError, ValueError } from '../errors.js';
import { Matrix } from '../matrix/matrix_generic.js';
import { adjugate, determinant, inverse } from '../matrix/matrix_operations.js';
import { type IntegerRing, ZZ } from '../rings/integer_ring.js';
import { Rational } from '../rings/rational.js';
import { QQ, type RationalField } from '../rings/rational_field.js';
import { type IntegerLike, type RationalLike, toBigInt, toRational } from '../types/coercion.js';
// Cyclic by design (see `require_local_field_invariants` at the bottom of the
// file); only ever dereferenced from inside a method body.
import * as LFI from './quadratic_form__local_field_invariants.js';

/**
 * The base rings supported by this port.
 *
 * @see Deviation: SageMath's `QuadraticForm` accepts any commutative ring.  We
 * support `ZZ` and `QQ`; everything else raises `NotImplementedError`.
 */
export type QFBaseRing = IntegerRing | RationalField;

/**
 * A matrix over `QQ`.
 *
 * This is a structurally typed view of `sage/matrix/matrix_generic.Matrix`
 * specialised to `Rational` entries: the repo's `Matrix<R extends RingElement>`
 * constraint is not satisfied by `Rational` (whose `add` accepts
 * `Rational | IntegerLike` rather than exactly `this`), so the generic class
 * cannot be instantiated with `Rational` under `strict` type checking.  The
 * runtime objects created here *are* `Matrix` instances, so all functions of
 * `sage/matrix/matrix_operations` apply to them unchanged.  This mirrors the
 * `RationalMatrix` façade already used by `sage/algebras/quatalg`.
 *
 * @see Deviation: `RationalMatrix` façade
 */
export interface RationalMatrix {
  readonly nrows: number;
  readonly ncols: number;
  get(i: number, j: number): Rational;
  set(i: number, j: number, value: Rational): void;
  row(i: number): Rational[];
  column(j: number): Rational[];
  rows(): Rational[][];
  columns(): Rational[][];
  list(): Rational[];
  mul(other: RationalMatrix): RationalMatrix;
  add(other: RationalMatrix): RationalMatrix;
  sub(other: RationalMatrix): RationalMatrix;
  scalar_mul(c: Rational): RationalMatrix;
  transpose(): RationalMatrix;
  is_square(): boolean;
  is_zero(): boolean;
  copy(): RationalMatrix;
  eq(other: RationalMatrix): boolean;
  toString(): string;
}

/** Construct a matrix over `QQ` (see {@link RationalMatrix}). */
const MatrixQQ = Matrix as unknown as new (
  ring: unknown,
  nrows: number,
  ncols: number,
  entries?: Rational[][] | Rational[] | ((i: number, j: number) => Rational)
) => RationalMatrix;

/** `determinant` of `sage/matrix/matrix_operations`, typed for `QQ`. */
export const determinantQQ = determinant as unknown as (M: RationalMatrix) => Rational;

/** `inverse` of `sage/matrix/matrix_operations`, typed for `QQ`. */
export const inverseQQ = inverse as unknown as (M: RationalMatrix) => RationalMatrix;

/** `adjugate` of `sage/matrix/matrix_operations`, typed for `QQ`. */
export const adjugateQQ = adjugate as unknown as (M: RationalMatrix) => RationalMatrix;

/** Runtime test for a {@link RationalMatrix} (i.e. a `Matrix` instance). */
export function isRationalMatrix(x: unknown): x is RationalMatrix {
  return x instanceof (Matrix as unknown as new (...args: never[]) => unknown);
}

/** Anything we accept where SageMath wants a matrix over the base ring. */
export type MatrixLike = RationalMatrix | RationalLike[][];

/** Anything we accept where SageMath wants a vector over the base ring. */
export type VectorLike = RationalLike[];

/* ------------------------------------------------------------------ */
/* Small helpers                                                       */
/* ------------------------------------------------------------------ */

/** Build an `n x m` matrix over `QQ` from a flat row-major list. */
export function matrixQQ(nrows: number, ncols: number, entries: Rational[]): RationalMatrix {
  return new MatrixQQ(QQ, nrows, ncols, entries);
}

/** The `n x n` identity matrix over `QQ`. */
export function identityQQ(n: number): RationalMatrix {
  return new MatrixQQ(QQ, n, n, (i, j) => (i === j ? Rational.one() : Rational.zero()));
}

/** Normalise a {@link MatrixLike} to a {@link RationalMatrix}. */
export function toMatrixQQ(M: MatrixLike): RationalMatrix {
  if (isRationalMatrix(M)) {
    return M;
  }
  const nrows = M.length;
  const ncols = nrows === 0 ? 0 : M[0]!.length;
  const entries: Rational[] = [];
  for (let i = 0; i < nrows; i++) {
    if (M[i]!.length !== ncols) {
      throw new ValueError('matrix rows must all have the same length');
    }
    for (let j = 0; j < ncols; j++) {
      entries.push(toRational(M[i]![j]!));
    }
  }
  return matrixQQ(nrows, ncols, entries);
}

/**
 * Coerce a matrix/coefficient *index* to a JavaScript number.
 *
 * Indices are positions, not arithmetic values, so (unlike the values handled
 * by `toBigInt`) a plain `number` is accepted here.
 */
function toIndex(x: number | IntegerLike): number {
  if (typeof x === 'number') {
    if (!Number.isInteger(x)) {
      throw new TypeError('index must be an integer');
    }
    return x;
  }
  return Number(toBigInt(x));
}

/** `true` when the value is one of the two supported base rings. */
function isQFBaseRing(x: unknown): x is QFBaseRing {
  return x === ZZ || x === QQ;
}

/**
 * Coerce a rational into the ring `R`, mirroring `R(x)`.
 *
 * @throws {TypeError} if `R` is `ZZ` and `x` is not an integer.
 */
function coerceToRing(R: QFBaseRing, x: Rational): Rational {
  if (R === ZZ && !x.isInteger()) {
    throw new TypeError(`no conversion of ${x.toString()} to an element of Integer Ring`);
  }
  return x;
}

/* ------------------------------------------------------------------ */
/* The QuadraticForm class                                             */
/* ------------------------------------------------------------------ */

/**
 * A quadratic form in `n` variables with coefficients in `R`.
 *
 * The form is `sum_{i <= j} a_{ij} x_i x_j`; the coefficients are stored in
 * upper-triangular reading order, exactly as in SageMath.
 *
 * @example
 * ```typescript
 * const Q = new QuadraticForm(ZZ, 2n, [1n, 3n, 5n]);
 * Q.toString();
 * // Quadratic form in 2 variables over Integer Ring with coefficients:
 * // [ 1 3 ]
 * // [ * 5 ]
 * ```
 *
 * @see Reference: sage/quadratic_forms/quadratic_form.py:159 (class QuadraticForm)
 */
export class QuadraticForm {
  private readonly _base_ring: QFBaseRing;
  private readonly _n: number;
  private _coeffs: Rational[];

  /* caches, mirroring SageMath's `self.__det`, `self.__level`, ... */
  private _det_cache?: Rational;
  private _level_cache?: bigint;
  private _definiteness_string?: string;
  private _rational_diagonal_cache?: [QuadraticForm, RationalMatrix];
  private _theta_vec?: bigint[];

  /**
   * `QuadraticForm(R, n, entries)`, `QuadraticForm(R, M)` or `QuadraticForm(M)`.
   *
   * @see Reference: sage/quadratic_forms/quadratic_form.py:492 (`__init__`)
   * @see Deviation: in the one-argument matrix form SageMath takes the base
   * ring from the matrix; a {@link RationalMatrix} carries no such information
   * here, so we infer `ZZ` when every entry is integral and `QQ` otherwise.
   */
  constructor(R: QFBaseRing, n: IntegerLike, entries?: RationalLike[]);
  constructor(R: QFBaseRing, M: MatrixLike);
  constructor(M: MatrixLike);
  constructor(R: QFBaseRing | MatrixLike, n?: IntegerLike | MatrixLike, entries?: RationalLike[]) {
    let M: RationalMatrix | null = null;
    let M_ring: QFBaseRing | null = null;

    if (isQFBaseRing(R)) {
      if (n !== undefined && (n instanceof Matrix || Array.isArray(n))) {
        // QuadraticForm(R, matrix)
        const A = toMatrixQQ(n as MatrixLike);
        if (!QuadraticForm._is_even_symmetric_matrix_(A, R)) {
          throw new TypeError('the matrix is not a symmetric with even diagonal defined over R');
        }
        M = A;
        M_ring = R;
      }
    } else if (R instanceof Matrix || Array.isArray(R)) {
      // QuadraticForm(matrix)
      const A = toMatrixQQ(R);
      // SageMath takes the base ring from the matrix.  A raw array of bigints
      // is a matrix over ZZ; anything with a non-integral entry is over QQ.
      const inferred: QFBaseRing = A.list().every((x) => x.isInteger()) ? ZZ : QQ;
      if (!QuadraticForm._is_even_symmetric_matrix_(A, inferred)) {
        throw new TypeError('the matrix is not a symmetric with even diagonal');
      }
      M = A;
      M_ring = inferred;
    } else {
      throw new TypeError('wrong input for QuadraticForm');
    }

    if (M !== null && M_ring !== null) {
      this._n = M.nrows;
      this._base_ring = M_ring;
      this._coeffs = [];
      for (let i = 0; i < M.nrows; i++) {
        for (let j = i; j < M.nrows; j++) {
          if (i === j) {
            this._coeffs.push(M.get(i, j).div(2n));
          } else {
            this._coeffs.push(M.get(i, j));
          }
        }
      }
      return;
    }

    // QuadraticForm(R, n, entries)
    const _R = R as QFBaseRing;
    const _n = toBigInt(n as IntegerLike);
    if (_n < 0n) {
      throw new ValueError(`the size must be a nonnegative integer, not ${_n}`);
    }
    const nn = Number(_n);
    const N = (nn * (nn + 1)) / 2;
    this._n = nn;
    this._base_ring = _R;
    this._coeffs = new Array<Rational>(N).fill(Rational.zero());

    if (entries !== undefined) {
      if (!Array.isArray(entries)) {
        throw new TypeError('entries must be an iterable');
      }
      if (entries.length !== N) {
        throw new TypeError(`the entries ${entries} must be a list of size n(n+1)/2`);
      }
      for (let i = 0; i < N; i++) {
        this._coeffs[i] = coerceToRing(_R, toRational(entries[i]!));
      }
    }
  }

  /**
   * Test whether `A` is symmetric, square, defined over `R` and has even
   * diagonal in `R`.
   *
   * @see Reference: sage/quadratic_forms/quadratic_form.py:1067
   */
  static _is_even_symmetric_matrix_(A: MatrixLike, R?: QFBaseRing): boolean {
    const M = toMatrixQQ(A);
    if (!M.is_square()) {
      return false;
    }
    const n = M.nrows;
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        if (!M.get(i, j).eq(M.get(j, i))) {
          return false;
        }
      }
    }
    const ring: QFBaseRing = R ?? (M.list().every((x) => x.isInteger()) ? ZZ : QQ);

    // Test that all entries coerce to R.
    if (ring === ZZ) {
      for (let i = 0; i < n; i++) {
        for (let j = i; j < n; j++) {
          if (!M.get(i, j).isInteger()) {
            return false;
          }
        }
      }
      // 2 is not a unit in ZZ, so the diagonal must be even.
      for (let i = 0; i < n; i++) {
        const d = M.get(i, i);
        if (d.denominator !== 1n || d.numerator % 2n !== 0n) {
          return false;
        }
      }
    }
    return true;
  }

  /* --------------------------------------------------------------- */
  /* Basic accessors                                                  */
  /* --------------------------------------------------------------- */

  /**
   * The number of variables.
   *
   * @see Reference: sage/quadratic_forms/quadratic_form.py:1425
   */
  dim(): bigint {
    return BigInt(this._n);
  }

  /** The number of variables, as a JavaScript number (internal convenience). */
  get n(): number {
    return this._n;
  }

  /**
   * The ring over which the form is defined.
   *
   * @see Reference: sage/quadratic_forms/quadratic_form.py:1444
   */
  base_ring(): QFBaseRing {
    return this._base_ring;
  }

  /**
   * The upper-triangular coefficients, read across the rows from the diagonal.
   *
   * @see Reference: sage/quadratic_forms/quadratic_form.py:1456
   */
  coefficients(): Rational[] {
    return this._coeffs.slice();
  }

  /**
   * The coefficient `a_{ij}` of `x_i x_j` (`Q[i, j]` in SageMath).
   *
   * @see Reference: sage/quadratic_forms/quadratic_form.py:765 (`__getitem__`)
   */
  get(i: number | IntegerLike, j: number | IntegerLike): Rational {
    let _i = toIndex(i);
    let _j = toIndex(j);
    if (_i > _j) {
      [_i, _j] = [_j, _i];
    }
    const idx = _i * this._n - (_i * (_i - 1)) / 2 + _j - _i;
    const v = this._coeffs[idx];
    if (v === undefined) {
      throw new ValueError(`index (${_i}, ${_j}) out of range`);
    }
    return v;
  }

  /**
   * Set the coefficient `a_{ij}` (`Q[i, j] = coeff` in SageMath).
   *
   * Mutates the form and invalidates the caches, exactly as SageMath's
   * `__setitem__` (which does *not* invalidate; see the deviation note).
   *
   * @see Reference: sage/quadratic_forms/quadratic_form.py:788 (`__setitem__`)
   * @see Deviation: we drop the memoised determinant/level/diagonalisation when
   * a coefficient is overwritten.  SageMath keeps the stale cache.
   */
  set(i: number | IntegerLike, j: number | IntegerLike, coeff: RationalLike): void {
    let _i = toIndex(i);
    let _j = toIndex(j);
    if (_i > _j) {
      [_i, _j] = [_j, _i];
    }
    const idx = _i * this._n - (_i * (_i - 1)) / 2 + _j - _i;
    if (idx < 0 || idx >= this._coeffs.length) {
      throw new ValueError(`index (${_i}, ${_j}) out of range`);
    }
    let value: Rational;
    try {
      value = coerceToRing(this._base_ring, toRational(coeff));
    } catch {
      throw new RuntimeError(
        'this coefficient cannot be coerced to an element of the base ring for the quadratic form'
      );
    }
    this._coeffs[idx] = value;
    this._det_cache = undefined;
    this._level_cache = undefined;
    this._definiteness_string = undefined;
    this._rational_diagonal_cache = undefined;
    this._theta_vec = undefined;
  }

  /**
   * The text representation used by SageMath.
   *
   * @see Reference: sage/quadratic_forms/quadratic_form.py:714 (`_repr_`)
   */
  _repr_(): string {
    const n = this._n;
    let out = `Quadratic form in ${n} variables over ${this._base_ring.toString()} with coefficients: \n`;
    for (let i = 0; i < n; i++) {
      if (i > 0) {
        out += '\n';
      }
      out += '[ ';
      for (let j = 0; j < n; j++) {
        if (i > j) {
          out += '* ';
        } else {
          out += `${this.get(i, j).toString()} `;
        }
      }
      out += ']';
    }
    return out;
  }

  toString(): string {
    return this._repr_();
  }

  /**
   * LaTeX representation.
   *
   * @see Reference: sage/quadratic_forms/quadratic_form.py:739 (`_latex_`)
   */
  _latex_(): string {
    const n = this._n;
    let out = `Quadratic form in ${n} variables over ${this._base_ring.toString()}`;
    out += ' with coefficients: \\newline';
    out += `\\left[ \\begin{array}{${'c'.repeat(n)}}`;
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        if (i > j) {
          out += ' * & ';
        } else {
          out += `${this.get(i, j).toString()} & `;
        }
      }
    }
    out += '\\end{array} \\right]';
    return out;
  }

  /**
   * Equality: same base ring and same coefficients.
   *
   * @see Reference: sage/quadratic_forms/quadratic_form.py:838 (`__eq__`)
   */
  equals(right: unknown): boolean {
    if (!(right instanceof QuadraticForm)) {
      return false;
    }
    if (this._base_ring !== right._base_ring) {
      return false;
    }
    if (this._coeffs.length !== right._coeffs.length) {
      return false;
    }
    return this._coeffs.every((c, i) => c.eq(right._coeffs[i]!));
  }

  /* --------------------------------------------------------------- */
  /* Sums                                                             */
  /* --------------------------------------------------------------- */

  /**
   * The direct sum (`Q + Q'` in SageMath).
   *
   * @see Reference: sage/quadratic_forms/quadratic_form.py:864 (`__add__`)
   */
  add(right: QuadraticForm): QuadraticForm {
    if (!(right instanceof QuadraticForm)) {
      throw new TypeError('cannot add these objects since they are not both quadratic forms');
    }
    if (this._base_ring !== right._base_ring) {
      throw new TypeError(
        'cannot add these since the quadratic forms do not have the same base rings'
      );
    }
    const Q = new QuadraticForm(this._base_ring, this.dim() + right.dim());
    const n = this._n;
    const m = right._n;
    for (let i = 0; i < n; i++) {
      for (let j = i; j < n; j++) {
        Q.set(i, j, this.get(i, j));
      }
    }
    for (let i = 0; i < m; i++) {
      for (let j = i; j < m; j++) {
        Q.set(n + i, n + j, right.get(i, j));
      }
    }
    return Q;
  }

  /**
   * The coefficient-wise sum of two forms of the same size.
   *
   * @see Reference: sage/quadratic_forms/quadratic_form.py:901
   */
  sum_by_coefficients_with(right: QuadraticForm): QuadraticForm {
    if (!(right instanceof QuadraticForm)) {
      throw new TypeError('cannot add these objects since they are not both quadratic forms');
    }
    if (this._n !== right._n) {
      throw new TypeError('cannot add these since the quadratic forms do not have the same sizes');
    }
    if (this._base_ring !== right._base_ring) {
      throw new TypeError(
        'cannot add these since the quadratic forms do not have the same base rings'
      );
    }
    return new QuadraticForm(
      this._base_ring,
      this.dim(),
      this._coeffs.map((c, i) => c.add(right._coeffs[i]!))
    );
  }

  /* --------------------------------------------------------------- */
  /* Evaluation                                                       */
  /* --------------------------------------------------------------- */

  /**
   * Evaluate the form on a vector (returning a scalar) or on a matrix
   * (returning the form `M^t Q M`).
   *
   * @see Reference: sage/quadratic_forms/quadratic_form.py:954 (`__call__`)
   * @see Reference: sage/quadratic_forms/quadratic_form__evaluate.pyx
   */
  __call__(v: VectorLike): Rational;
  __call__(v: RationalMatrix): QuadraticForm;
  __call__(v: VectorLike | MatrixLike): Rational | QuadraticForm;
  __call__(v: VectorLike | MatrixLike): Rational | QuadraticForm {
    const n = this._n;
    // Distinguish a vector (flat list) from a matrix (list of rows / Matrix).
    const isMatrix = v instanceof Matrix || (Array.isArray(v) && Array.isArray(v[0]));
    if (isMatrix) {
      const M = toMatrixQQ(v as MatrixLike);
      if (M.nrows !== n) {
        throw new TypeError(`the matrix must have ${n} rows`);
      }
      const m = M.ncols;
      const Q2 = new QuadraticForm(this._base_ring, BigInt(m));
      return QFEvaluateMatrix(this, M, Q2);
    }
    const vec = v as VectorLike;
    if (vec.length !== n) {
      throw new TypeError(`your vector needs to have length ${n}`);
    }
    return QFEvaluateVector(this, vec.map(toRational));
  }

  /**
   * The value of the associated bilinear map `B(v, w)`.
   *
   * `2 B(v, w) = Q(v + w) - Q(v) - Q(w)`.
   *
   * @see Reference: sage/quadratic_forms/quadratic_form.py:1666
   */
  bilinear_map(v: VectorLike, w: VectorLike): Rational {
    if (v.length !== this._n || w.length !== this._n) {
      throw new TypeError(`vectors must have length ${this._n}`);
    }
    // Both ZZ and QQ have characteristic 0, so SageMath's characteristic-2
    // guard can never fire here.
    const _v = v.map(toRational);
    const _w = w.map(toRational);
    const sum = _v.map((x, i) => x.add(_w[i]!));
    return QFEvaluateVector(this, sum)
      .sub(QFEvaluateVector(this, _v))
      .sub(QFEvaluateVector(this, _w))
      .div(2n);
  }

  /* --------------------------------------------------------------- */
  /* Matrices and determinants                                        */
  /* --------------------------------------------------------------- */

  /**
   * The Hessian matrix `A`, for which `Q(X) = (1/2) X^t A X`.
   *
   * @see Reference: sage/quadratic_forms/quadratic_form.py:1138
   * @see Deviation: the returned matrix always has `QQ` as its base ring, even
   * for forms over `ZZ` (SageMath returns a matrix over the form's base ring).
   * The entries are unchanged.
   */
  Hessian_matrix(): RationalMatrix {
    const n = this._n;
    const entries: Rational[] = [];
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        entries.push(i === j ? this.get(i, j).mul(2n) : this.get(i, j));
      }
    }
    return matrixQQ(n, n, entries);
  }

  /**
   * Alias for {@link Hessian_matrix}.
   *
   * @see Reference: sage/quadratic_forms/quadratic_form.py:1124
   */
  matrix(): RationalMatrix {
    return this.Hessian_matrix();
  }

  /**
   * A Gram matrix `A` with `Q(x) = x^t A x`, over the fraction field.
   *
   * @see Reference: sage/quadratic_forms/quadratic_form.py:1164
   */
  Gram_matrix_rational(): RationalMatrix {
    return this.matrix().scalar_mul(new Rational(1n, 2n));
  }

  /**
   * A Gram matrix over the base ring; raises if the form is not integral.
   *
   * @see Reference: sage/quadratic_forms/quadratic_form.py:1188
   */
  Gram_matrix(): RationalMatrix {
    const A = this.Gram_matrix_rational();
    if (this._base_ring === ZZ) {
      const n = this._n;
      for (let i = 0; i < n; i++) {
        for (let j = i; j < n; j++) {
          if (!A.get(i, j).isInteger()) {
            throw new TypeError('this form does not have an integral Gram matrix');
          }
        }
      }
    }
    return A;
  }

  /**
   * Whether the Gram matrix is integral over the base ring.
   *
   * @see Reference: sage/quadratic_forms/quadratic_form.py:1225
   */
  has_integral_Gram_matrix(): boolean {
    if (this._base_ring === QQ) {
      console.warn(
        'Warning -- A quadratic form over a field always has integral Gram matrix.  Do you really want to do this?!?'
      );
    }
    try {
      this.Gram_matrix();
      return true;
    } catch {
      return false;
    }
  }

  /**
   * `det(2Q)`: the determinant of the Hessian matrix.
   *
   * @see Reference: sage/quadratic_forms/quadratic_form.py:1469
   * @see Deviation: returns a {@link Rational} even over `ZZ` (SageMath returns
   * an element of the base ring).  Use `.numerator` when you need a `bigint`.
   */
  det(): Rational {
    if (this._det_cache === undefined) {
      this._det_cache = this._n === 0 ? Rational.one() : determinantQQ(this.matrix());
    }
    return this._det_cache;
  }

  /**
   * `det(Q) = det(2Q) / 2^n`.
   *
   * @see Reference: sage/quadratic_forms/quadratic_form.py:1497
   */
  Gram_det(): Rational {
    return this.det().div(new Rational(2n ** BigInt(this._n), 1n));
  }

  /* --------------------------------------------------------------- */
  /* Content, primitivity, level                                      */
  /* --------------------------------------------------------------- */

  /**
   * The gcd of the coefficients (over `ZZ` only).
   *
   * @see Reference: sage/quadratic_forms/quadratic_form.py:1256
   */
  gcd(): bigint {
    if (this._base_ring !== ZZ) {
      throw new TypeError('the given quadratic form must be defined over ZZ');
    }
    return GCD(this._coeffs.map((c) => c.numerator));
  }

  /**
   * Whether the (integer-valued) form is primitive.
   *
   * @see Reference: sage/quadratic_forms/quadratic_form.py:1367
   */
  is_primitive(): boolean {
    return this.gcd() === 1n;
  }

  /**
   * The primitive form obtained by dividing out the content.
   *
   * @see Reference: sage/quadratic_forms/quadratic_form.py:1385
   */
  primitive(): QuadraticForm {
    if (this._base_ring !== ZZ) {
      throw new TypeError('the given quadratic form must be defined over ZZ');
    }
    const g = this.gcd();
    return new QuadraticForm(
      ZZ,
      this.dim(),
      this._coeffs.map((x) => x.numerator / g)
    );
  }

  /**
   * The primitive adjoint form.
   *
   * @see Reference: sage/quadratic_forms/quadratic_form.py:1408
   */
  adjoint_primitive(): QuadraticForm {
    return new QuadraticForm(adjugateQQ(this.Hessian_matrix())).primitive();
  }

  /**
   * The level of the form over a PID.
   *
   * @see Reference: sage/quadratic_forms/quadratic_form.py:1558
   */
  level(): bigint {
    if (this._level_cache !== undefined) {
      return this._level_cache;
    }
    if (this._base_ring === QQ) {
      console.warn(
        'Warning -- The level of a quadratic form over a field is always 1.  Do you really want to do this?!?'
      );
    }
    let mat_inv: RationalMatrix;
    try {
      mat_inv = inverseQQ(this.matrix());
    } catch {
      throw new TypeError('the quadratic form is degenerate');
    }
    const inv_denoms: bigint[] = [];
    for (let i = 0; i < this._n; i++) {
      for (let j = i; j < this._n; j++) {
        if (i === j) {
          inv_denoms.push(mat_inv.get(i, j).div(2n).denominator);
        } else {
          inv_denoms.push(mat_inv.get(i, j).denominator);
        }
      }
    }
    let lvl = LCM(inv_denoms);
    if (this._base_ring === ZZ) {
      lvl = lvl < 0n ? -lvl : lvl;
    }
    this._level_cache = lvl;
    return lvl;
  }

  /**
   * Change the base ring of the form.
   *
   * @see Reference: sage/quadratic_forms/quadratic_form.py:1515
   */
  change_ring(R: QFBaseRing): QuadraticForm {
    if (!isQFBaseRing(R)) {
      throw new NotImplementedError(
        'SAGE_NOT_IMPLEMENTED: QuadraticForm.change_ring only supports ZZ and QQ'
      );
    }
    if (R === ZZ && this._base_ring === QQ) {
      throw new TypeError('there is no canonical coercion from Rational Field to R');
    }
    return new QuadraticForm(R, this.dim(), this.coefficients());
  }

  /* --------------------------------------------------------------- */
  /* Local field invariants (quadratic_form__local_field_invariants)  */
  /* --------------------------------------------------------------- */

  /**
   * @see {@link module:sage/quadratic_forms/quadratic_form__local_field_invariants}
   */
  rational_diagonal_form(): QuadraticForm;
  rational_diagonal_form(return_matrix: false): QuadraticForm;
  rational_diagonal_form(return_matrix: true): [QuadraticForm, RationalMatrix];
  rational_diagonal_form(return_matrix = false): QuadraticForm | [QuadraticForm, RationalMatrix] {
    const [Q, T] = this._rational_diagonal_form_and_transformation();
    // SageMath deep-copies the form because quadratic forms are mutable.
    const Qc = new QuadraticForm(Q.base_ring(), Q.dim(), Q.coefficients());
    return return_matrix ? [Qc, T] : Qc;
  }

  /** @internal cached `(D, T)` with `T^t * self.matrix() * T == D.matrix()`. */
  _rational_diagonal_form_and_transformation(): [QuadraticForm, RationalMatrix] {
    if (this._rational_diagonal_cache === undefined) {
      // Lazily required to avoid an import cycle at module-evaluation time.
      const mod = require_local_field_invariants();
      this._rational_diagonal_cache = mod._rational_diagonal_form_and_transformation(this);
    }
    return this._rational_diagonal_cache;
  }

  /** @see quadratic_form__local_field_invariants.signature_vector */
  signature_vector(): [bigint, bigint, bigint] {
    return require_local_field_invariants().signature_vector(this);
  }

  /** @see quadratic_form__local_field_invariants.signature */
  signature(): bigint {
    return require_local_field_invariants().signature(this);
  }

  /** @see quadratic_form__local_field_invariants.hasse_invariant */
  hasse_invariant(p: IntegerLike): bigint {
    return require_local_field_invariants().hasse_invariant(this, p);
  }

  /** @see quadratic_form__local_field_invariants.hasse_invariant__OMeara */
  hasse_invariant__OMeara(p: IntegerLike): bigint {
    return require_local_field_invariants().hasse_invariant__OMeara(this, p);
  }

  /** @see quadratic_form__local_field_invariants.is_hyperbolic */
  is_hyperbolic(p: IntegerLike): boolean {
    return require_local_field_invariants().is_hyperbolic(this, p);
  }

  /** @see quadratic_form__local_field_invariants.is_anisotropic */
  is_anisotropic(p: IntegerLike): boolean {
    return require_local_field_invariants().is_anisotropic(this, p);
  }

  /** @see quadratic_form__local_field_invariants.is_isotropic */
  is_isotropic(p: IntegerLike): boolean {
    return require_local_field_invariants().is_isotropic(this, p);
  }

  /** @see quadratic_form__local_field_invariants.anisotropic_primes */
  anisotropic_primes(): bigint[] {
    return require_local_field_invariants().anisotropic_primes(this);
  }

  /** @see quadratic_form__local_field_invariants.compute_definiteness */
  compute_definiteness(): void {
    this._definiteness_string = require_local_field_invariants().compute_definiteness_string(this);
  }

  /** @internal the cached definiteness string, computing it if necessary. */
  _definiteness(): string {
    if (this._definiteness_string === undefined) {
      this.compute_definiteness();
    }
    return this._definiteness_string as string;
  }

  /** @see quadratic_form__local_field_invariants.compute_definiteness_string_by_determinants */
  compute_definiteness_string_by_determinants(): string {
    return require_local_field_invariants().compute_definiteness_string_by_determinants(this);
  }

  /** @see Reference: sage/quadratic_forms/quadratic_form__local_field_invariants.py:924 */
  is_positive_definite(): boolean {
    const s = this._definiteness();
    return s === 'pos_def' || s === 'zero';
  }

  /** @see Reference: sage/quadratic_forms/quadratic_form__local_field_invariants.py:957 */
  is_negative_definite(): boolean {
    const s = this._definiteness();
    return s === 'neg_def' || s === 'zero';
  }

  /** @see Reference: sage/quadratic_forms/quadratic_form__local_field_invariants.py:991 */
  is_indefinite(): boolean {
    return this._definiteness() === 'indefinite';
  }

  /** @see Reference: sage/quadratic_forms/quadratic_form__local_field_invariants.py:1025 */
  is_definite(): boolean {
    const s = this._definiteness();
    return s === 'pos_def' || s === 'neg_def' || s === 'zero';
  }

  /* --------------------------------------------------------------- */
  /* Equivalence testing                                              */
  /* --------------------------------------------------------------- */

  /**
   * Determine whether two regular quadratic forms over `QQ` are isometric
   * over `QQ` (Hasse--Minkowski).
   *
   * @see Reference: sage/quadratic_forms/quadratic_form__equivalence_testing.py:is_rationally_isometric
   * @see Deviation: SageMath supports `QQ` and number fields; over `ZZ` it
   * raises `AttributeError` (`IntegerRing` has no `real_embeddings`).  We
   * raise `NotImplementedError` naming the missing support instead.
   */
  is_rationally_isometric(other: QuadraticForm): boolean {
    if (this.Gram_det().isZero() || other.Gram_det().isZero()) {
      throw new NotImplementedError('this only tests regular forms');
    }
    if (this._base_ring !== other._base_ring) {
      throw new TypeError('forms must have the same base ring.');
    }
    if (this._n !== other._n) {
      return false;
    }
    if (this._base_ring !== QQ) {
      throw new NotImplementedError(
        'SAGE_NOT_IMPLEMENTED: is_rationally_isometric over Integer Ring ' +
          '(SageMath only implements this for QQ and number fields; call ' +
          'change_ring(QQ) first)'
      );
    }
    if (!this.Gram_det().mul(other.Gram_det()).is_square()) {
      return false;
    }
    const L1 = this.Gram_det().support();
    const L2 = other.Gram_det().support();
    const primes = new Set<bigint>([...L1, ...L2]);
    for (const p of primes) {
      if (this.hasse_invariant(p) !== other.hasse_invariant(p)) {
        return false;
      }
    }
    return this.signature() === other.signature();
  }

  /* --------------------------------------------------------------- */
  /* Theta series (delegated to PARI's qfrep)                         */
  /* --------------------------------------------------------------- */

  /**
   * The list `[r(0), r(1), ..., r(B-1)]` of representation numbers.
   *
   * SageMath calls `pari(1).concat(self.__pari__().qfrep(B - 1, 1) * 2)`; we
   * delegate to `parigp-ts`'s port of `qfrep` in exactly the same way.
   *
   * @see Reference: sage/quadratic_forms/quadratic_form__ternary_Tornaria.py:representation_number_list
   */
  representation_number_list(B: IntegerLike): bigint[] {
    const _B = toBigInt(B);
    if (_B < 1n) {
      throw new ValueError('B must be positive');
    }
    const H = this.Hessian_matrix();
    const rows: bigint[][] = [];
    for (let i = 0; i < this._n; i++) {
      const row: bigint[] = [];
      for (let j = 0; j < this._n; j++) {
        const e = H.get(i, j);
        if (!e.isInteger()) {
          throw new TypeError('representation_number_list requires an integral Hessian matrix');
        }
        row.push(e.numerator);
      }
      rows.push(row);
    }
    const rep = qfrep(rows, _B - 1n, 1);
    return [1n, ...rep.map((x) => x * 2n)];
  }

  /**
   * The theta series coefficients up to `O(q^Max)`.
   *
   * @see Reference: sage/quadratic_forms/quadratic_form__theta.py:theta_by_pari
   * @see Deviation: SageMath returns a power series in `ZZ[[q]]` when
   * `var_str` is nonempty; we always return the coefficient vector (i.e.
   * SageMath's `var_str=''` behaviour), because that is the only shape the
   * result is consumed in here.
   */
  theta_by_pari(Max: IntegerLike): bigint[] {
    const M = toBigInt(Max);
    if (this._theta_vec !== undefined && BigInt(this._theta_vec.length) >= M) {
      return this._theta_vec.slice(0, Number(M));
    }
    const v = this.representation_number_list(M);
    this._theta_vec = v;
    return v.slice();
  }

  /**
   * Alias for {@link theta_by_pari}.
   *
   * @see Reference: sage/quadratic_forms/quadratic_form__theta.py:theta_series
   */
  theta_series(Max: IntegerLike = 10n): bigint[] {
    const M = toBigInt(Max);
    if (M < 0n) {
      throw new TypeError(`Max = ${Max} is not an integer >= 0 or an allowed string`);
    }
    return this.theta_by_pari(M);
  }

  /* --------------------------------------------------------------- */
  /* Honest stubs                                                     */
  /* --------------------------------------------------------------- */

  /**
   * @see Reference: sage/quadratic_forms/quadratic_form.py:1275
   * @throws {NotImplementedError} always
   */
  polynomial(_names = 'x'): never {
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: QuadraticForm.polynomial (needs the multivariate polynomial ring plumbing)'
    );
  }

  /**
   * @see Reference: sage/quadratic_forms/quadratic_form.py:1327
   * @throws {NotImplementedError} always
   */
  static from_polynomial(_poly: unknown): never {
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: QuadraticForm.from_polynomial (needs the multivariate polynomial ring plumbing)'
    );
  }

  /**
   * @see Reference: sage/quadratic_forms/quadratic_form__local_normal_form.py:local_normal_form
   * @throws {NotImplementedError} always
   */
  local_normal_form(_p: IntegerLike): never {
    throw new NotImplementedError('SAGE_NOT_IMPLEMENTED: QuadraticForm.local_normal_form');
  }

  /**
   * @see Reference: sage/quadratic_forms/quadratic_form__local_normal_form.py:jordan_blocks_by_scale_and_unimodular
   * @throws {NotImplementedError} always
   */
  jordan_blocks_by_scale_and_unimodular(_p: IntegerLike): never {
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: QuadraticForm.jordan_blocks_by_scale_and_unimodular ' +
        '(requires local_normal_form)'
    );
  }

  /**
   * @see Reference: sage/quadratic_forms/quadratic_form__genus.py:local_genus_symbol
   * @throws {NotImplementedError} always
   */
  local_genus_symbol(_p: IntegerLike): never {
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: QuadraticForm.local_genus_symbol (requires sage.quadratic_forms.genera.genus)'
    );
  }

  /**
   * @see Reference: sage/quadratic_forms/quadratic_form__genus.py:global_genus_symbol
   * @throws {NotImplementedError} always
   */
  global_genus_symbol(): never {
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: QuadraticForm.global_genus_symbol (requires sage.quadratic_forms.genera.genus)'
    );
  }

  /**
   * @see Reference: sage/quadratic_forms/quadratic_form__genus.py:CS_genus_symbol_list
   * @throws {NotImplementedError} always
   */
  CS_genus_symbol_list(): never {
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: QuadraticForm.CS_genus_symbol_list (requires sage.quadratic_forms.genera.genus)'
    );
  }

  /**
   * @see Reference: sage/quadratic_forms/quadratic_form__equivalence_testing.py:is_globally_equivalent_to
   * @throws {NotImplementedError} always
   */
  is_globally_equivalent_to(_other: QuadraticForm): never {
    throw new NotImplementedError(
      "SAGE_NOT_IMPLEMENTED: QuadraticForm.is_globally_equivalent_to (requires PARI's qfisominit/qfisom)"
    );
  }

  /**
   * @see Reference: sage/quadratic_forms/quadratic_form__equivalence_testing.py:is_locally_equivalent_to
   * @throws {NotImplementedError} always
   */
  is_locally_equivalent_to(_other: QuadraticForm): never {
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: QuadraticForm.is_locally_equivalent_to (requires the genus symbols)'
    );
  }

  /**
   * @see Reference: sage/quadratic_forms/qfsolve.py:solve
   * @throws {NotImplementedError} always
   */
  solve(_c: IntegerLike = 0n): never {
    throw new NotImplementedError(
      "SAGE_NOT_IMPLEMENTED: QuadraticForm.solve (requires PARI's qfsolve/qfparam)"
    );
  }

  /**
   * @see Reference: sage/quadratic_forms/quadratic_form__automorphisms.py:automorphisms
   * @throws {NotImplementedError} always
   */
  automorphisms(): never {
    throw new NotImplementedError(
      "SAGE_NOT_IMPLEMENTED: QuadraticForm.automorphisms (requires PARI's qfauto)"
    );
  }

  /**
   * @see Reference: sage/quadratic_forms/quadratic_form__automorphisms.py:number_of_automorphisms
   * @throws {NotImplementedError} always
   */
  number_of_automorphisms(): never {
    throw new NotImplementedError(
      "SAGE_NOT_IMPLEMENTED: QuadraticForm.number_of_automorphisms (requires PARI's qfauto)"
    );
  }

  /**
   * @see Reference: sage/quadratic_forms/quadratic_form__siegel_product.py:siegel_product
   * @throws {NotImplementedError} always
   */
  siegel_product(_u: IntegerLike): never {
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: QuadraticForm.siegel_product (requires the local density machinery)'
    );
  }

  /**
   * @see Reference: sage/quadratic_forms/quadratic_form__local_density_interfaces.py:local_density
   * @throws {NotImplementedError} always
   */
  local_density(_p: IntegerLike, _m: IntegerLike): never {
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: QuadraticForm.local_density (requires quadratic_form__count_local_2)'
    );
  }

  /**
   * @see Reference: sage/quadratic_forms/quadratic_form__mass.py
   * @throws {NotImplementedError} always
   */
  conway_mass(): never {
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: QuadraticForm.conway_mass (requires the Conway--Sloane mass machinery)'
    );
  }
}

/* ------------------------------------------------------------------ */
/* Evaluation helpers (sage/quadratic_forms/quadratic_form__evaluate)  */
/* ------------------------------------------------------------------ */

/**
 * Evaluate `Q` on the vector `v`.
 *
 * @see Reference: sage/quadratic_forms/quadratic_form__evaluate.pyx:QFEvaluateVector
 */
export function QFEvaluateVector(Q: QuadraticForm, v: Rational[]): Rational {
  const n = Q.n;
  let tmp = Rational.zero();
  for (let i = 0; i < n; i++) {
    for (let j = i; j < n; j++) {
      tmp = tmp.add(Q.get(i, j).mul(v[i]!).mul(v[j]!));
    }
  }
  return tmp;
}

/**
 * Evaluate `Q` on the matrix `M`, writing `M^t Q M` into `Q2`.
 *
 * @see Reference: sage/quadratic_forms/quadratic_form__evaluate.pyx:QFEvaluateMatrix
 */
export function QFEvaluateMatrix(
  Q: QuadraticForm,
  M: RationalMatrix,
  Q2: QuadraticForm
): QuadraticForm {
  const n = Q.n;
  const m = Q2.n;
  for (let k = 0; k < m; k++) {
    for (let l = k; l < m; l++) {
      let tmp = Rational.zero();
      if (k === l) {
        for (let i = 0; i < n; i++) {
          for (let j = i; j < n; j++) {
            tmp = tmp.add(Q.get(i, j).mul(M.get(i, k).mul(M.get(j, l))));
          }
        }
      } else {
        for (let i = 0; i < n; i++) {
          for (let j = i; j < n; j++) {
            tmp = tmp.add(
              Q.get(i, j).mul(
                M.get(i, k)
                  .mul(M.get(j, l))
                  .add(M.get(i, l).mul(M.get(j, k)))
              )
            );
          }
        }
      }
      Q2.set(k, l, tmp);
    }
  }
  return Q2;
}

/* ------------------------------------------------------------------ */
/* Constructors                                                        */
/* ------------------------------------------------------------------ */

/**
 * The diagonal quadratic form `sum diag[i] x_i^2` over `R`.
 *
 * @see Reference: sage/quadratic_forms/quadratic_form.py:1727
 */
export function DiagonalQuadraticForm(R: QFBaseRing, diag: RationalLike[]): QuadraticForm {
  const Q = new QuadraticForm(R, BigInt(diag.length));
  for (let i = 0; i < diag.length; i++) {
    Q.set(i, i, diag[i]!);
  }
  return Q;
}

/**
 * Return a rational quadratic form with the given invariants.
 *
 * @see Reference: sage/quadratic_forms/quadratic_form.py:48
 * @throws {NotImplementedError} always -- the rank-2 step of Kirschmer's
 * algorithm needs `QQ.hilbert_symbol_negative_at_S`, which this port does not
 * have yet.
 */
export function quadratic_form_from_invariants(
  _F: QFBaseRing,
  _rk: IntegerLike,
  _det: RationalLike,
  _P: IntegerLike[],
  _sminus: IntegerLike
): never {
  throw new NotImplementedError(
    'SAGE_NOT_IMPLEMENTED: quadratic_form_from_invariants ' +
      '(the rank-2 step needs QQ.hilbert_symbol_negative_at_S, ' +
      'sage/rings/rational_field.py, which is not ported yet)'
  );
}

/* ------------------------------------------------------------------ */
/* The local-field-invariants module                                   */
/* ------------------------------------------------------------------ */

/**
 * `quadratic_form__local_field_invariants.ts` imports {@link QuadraticForm},
 * so the dependency is cyclic -- exactly as in SageMath, where the module's
 * functions are imported *into* the `QuadraticForm` class body.  The cycle is
 * safe because neither module touches the other at module-evaluation time:
 * every use is inside a function body, and ESM hoists function declarations.
 */
function require_local_field_invariants(): typeof LFI {
  return LFI;
}
