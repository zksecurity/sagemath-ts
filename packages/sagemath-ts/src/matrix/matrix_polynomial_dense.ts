/**
 * @module sage/matrix/matrix_polynomial_dense
 * @description Dense matrices over univariate polynomials over fields
 *
 * Port of: `sage/matrix/matrix_polynomial_dense.pyx`
 *
 * SageMath implements these operations as methods of the Cython class
 * ``Matrix_polynomial_dense``.  This port has no such subclass -- a polynomial
 * matrix is simply a `Matrix<Polynomial<C>>` whose base ring is a
 * `PolynomialRing` -- so each method becomes a free function taking the matrix
 * as its first argument, exactly as the rest of `src/matrix/` does (compare
 * `matrix_operations.ts`, which turns `matrix2.pyx` methods into free
 * functions).  Method names, argument names, defaults, error messages and
 * outputs are otherwise identical to SageMath's.
 *
 * For a field `K`, we consider matrices over the univariate polynomial ring
 * `K[x]`.  They are often used to represent bases of `K[x]`-modules; see the
 * class description at `matrix_polynomial_dense.pyx:52-102` for the row-wise
 * versus column-wise conventions and for the definition of shifts.
 */

import { ValueError } from '../errors.js';
import type {
  CoefficientRing,
  Polynomial,
  PolynomialRingBase,
  RingElement,
} from '../rings/polynomial/polynomial_element.js';
import { Matrix, identity_matrix } from './matrix_generic.js';
import { determinant, inverse, rank } from './matrix_operations.js';

// ============================================================================
// Types and small helpers
// ============================================================================

/**
 * Interface for base-field elements: SageMath requires the base ring of the
 * polynomial ring to be a field for all algorithms in this module.
 */
export interface FieldElement extends RingElement {
  inv?(): this;
  inverse?(): this;
  div?(other: this): this;
}

/** A dense matrix over `K[x]`, i.e. SageMath's ``Matrix_polynomial_dense``. */
export type PolynomialMatrix<C extends FieldElement> = Matrix<Polynomial<C>>;

/** Multiplicative inverse of a base-field element. */
function _inv<C extends FieldElement>(c: C): C {
  if (typeof c.inv === 'function') return c.inv() as C;
  if (typeof c.inverse === 'function') return c.inverse() as C;
  throw new ValueError('base ring element does not support inversion');
}

/** Quotient `a / b` of two base-field elements. */
function _div<C extends FieldElement>(a: C, b: C): C {
  if (typeof a.div === 'function') return a.div(b) as C;
  return a.mul(_inv(b)) as C;
}

/** The polynomial ring `K[x]` a polynomial matrix lives over. */
function _pring<C extends FieldElement>(mat: PolynomialMatrix<C>): PolynomialRingBase<C> {
  return mat.base_ring as unknown as PolynomialRingBase<C>;
}

/** The base field `K`. */
function _field<C extends FieldElement>(mat: PolynomialMatrix<C>): CoefficientRing<C> {
  return _pring(mat).base_ring;
}

/** Build a polynomial from a coefficient list (index = degree). */
function _fromCoeffs<C extends FieldElement>(
  pR: PolynomialRingBase<C>,
  coeffs: C[]
): Polynomial<C> {
  return pR.__call__(coeffs);
}

/** The polynomial `c * x^d` with `c` in the base field and `d >= 0`. */
function _monomial<C extends FieldElement>(
  pR: PolynomialRingBase<C>,
  c: C,
  d: number
): Polynomial<C> {
  const coeffs: C[] = [];
  for (let i = 0; i < d; i++) coeffs.push(pR.base_ring.zero());
  coeffs.push(c);
  return _fromCoeffs(pR, coeffs);
}

/** Build a matrix over `ring` from a list of rows. */
function _fromRows<R extends RingElement>(
  ring: CoefficientRing<R>,
  m: number,
  n: number,
  rows: R[][]
): Matrix<R> {
  return new Matrix<R>(ring, m, n, rows);
}

/**
 * Mimic SageMath's ``matrix([[...] for i in range(m)])``: when the outer list
 * is empty the constructor produces the `0 x 0` matrix, not `0 x n`.  Used by
 * :func:`constant_matrix` and :func:`leading_matrix`, the two places where
 * upstream builds a matrix without passing explicit dimensions.
 */
function _fromRowsInferred<R extends RingElement>(
  ring: CoefficientRing<R>,
  m: number,
  n: number,
  rows: R[][]
): Matrix<R> {
  return m === 0 ? new Matrix<R>(ring, 0, 0) : new Matrix<R>(ring, m, n, rows);
}

/** `newRows[i] = oldRows[perm[i] - 1]`, SageMath's ``permute_rows``. */
function _permuteRows<R extends RingElement>(mat: Matrix<R>, perm: number[]): Matrix<R> {
  const rows: R[][] = perm.map((p) => mat.row(p - 1));
  return new Matrix<R>(mat.base_ring, mat.nrows, mat.ncols, rows);
}

/** The submatrix formed by the first `k` rows. */
function _firstRows<R extends RingElement>(mat: Matrix<R>, k: number): Matrix<R> {
  const rows: R[][] = [];
  for (let i = 0; i < k; i++) rows.push(mat.row(i));
  return new Matrix<R>(mat.base_ring, k, mat.ncols, rows);
}

/** The submatrix formed by the rows from index `k` onwards. */
function _rowsFrom<R extends RingElement>(mat: Matrix<R>, k: number): Matrix<R> {
  const rows: R[][] = [];
  for (let i = k; i < mat.nrows; i++) rows.push(mat.row(i));
  return new Matrix<R>(mat.base_ring, mat.nrows - k, mat.ncols, rows);
}

/** Vertical concatenation, SageMath's ``stack``. */
function _stack<R extends RingElement>(top: Matrix<R>, bottom: Matrix<R>): Matrix<R> {
  const rows: R[][] = [];
  for (let i = 0; i < top.nrows; i++) rows.push(top.row(i));
  for (let i = 0; i < bottom.nrows; i++) rows.push(bottom.row(i));
  return new Matrix<R>(top.base_ring, top.nrows + bottom.nrows, top.ncols, rows);
}

/** Horizontal concatenation, SageMath's ``augment``. */
function _augment<R extends RingElement>(left: Matrix<R>, right: Matrix<R>): Matrix<R> {
  const rows: R[][] = [];
  for (let i = 0; i < left.nrows; i++) rows.push([...left.row(i), ...right.row(i)]);
  return new Matrix<R>(left.base_ring, left.nrows, left.ncols + right.ncols, rows);
}

/** The submatrix formed by the given columns, SageMath's ``matrix_from_columns``. */
function _columnsOf<R extends RingElement>(mat: Matrix<R>, cols: number[]): Matrix<R> {
  const rows: R[][] = [];
  for (let i = 0; i < mat.nrows; i++) rows.push(cols.map((j) => mat.get(i, j)));
  return new Matrix<R>(mat.base_ring, mat.nrows, cols.length, rows);
}

/** Left-multiply a polynomial matrix by a constant matrix over the base field. */
function _constTimesPoly<C extends FieldElement>(
  cst: Matrix<C>,
  pol: PolynomialMatrix<C>
): PolynomialMatrix<C> {
  const pR = _pring(pol);
  const rows: Array<Array<Polynomial<C>>> = [];
  for (let i = 0; i < cst.nrows; i++) {
    const row: Array<Polynomial<C>> = [];
    for (let j = 0; j < pol.ncols; j++) {
      let s = pR.zero();
      for (let k = 0; k < cst.ncols; k++) {
        s = s.add(pol.get(k, j).scalar_mul(cst.get(i, k)));
      }
      row.push(s);
    }
    rows.push(row);
  }
  return _fromRows(pol.base_ring, cst.nrows, pol.ncols, rows);
}

/** Right-multiply a polynomial matrix by a constant matrix over the base field. */
function _polyTimesConst<C extends FieldElement>(
  pol: PolynomialMatrix<C>,
  cst: Matrix<C>
): PolynomialMatrix<C> {
  const pR = _pring(pol);
  const rows: Array<Array<Polynomial<C>>> = [];
  for (let i = 0; i < pol.nrows; i++) {
    const row: Array<Polynomial<C>> = [];
    for (let j = 0; j < cst.ncols; j++) {
      let s = pR.zero();
      for (let k = 0; k < pol.ncols; k++) {
        s = s.add(pol.get(i, k).scalar_mul(cst.get(k, j)));
      }
      row.push(s);
    }
    rows.push(row);
  }
  return _fromRows(pol.base_ring, pol.nrows, cst.ncols, rows);
}

/** Normalize an ``int or list`` argument into a list, as upstream does. */
function _asList(d: number | number[], len: number): number[] {
  return Array.isArray(d) ? [...d] : new Array<number>(len).fill(d);
}

// ============================================================================
// Basic degree / coefficient extraction
// ============================================================================

/**
 * Raise an exception if the ``shifts`` argument does not have the right length.
 *
 * For an `m x n` polynomial matrix, if working row-wise then ``shifts`` should
 * have `n` entries; if working column-wise, it should have `m` entries.
 *
 * @see Reference: sage/matrix/matrix_polynomial_dense.pyx:104-135
 */
export function _check_shift_dimension<C extends FieldElement>(
  mat: PolynomialMatrix<C>,
  shifts?: number[] | null,
  row_wise: boolean = true
): void {
  if (shifts !== undefined && shifts !== null && !row_wise && shifts.length !== mat.nrows) {
    throw new ValueError('shifts length should be the row dimension');
  }
  if (shifts !== undefined && shifts !== null && row_wise && shifts.length !== mat.ncols) {
    throw new ValueError('shifts length should be the column dimension');
  }
}

/**
 * Return the degree of this matrix: the maximum of the degrees of all its
 * entries.  The degree of the zero matrix (including empty matrices) is `-1`.
 *
 * @see Reference: sage/matrix/matrix_polynomial_dense.pyx:137-166
 */
export function degree<C extends FieldElement>(mat: PolynomialMatrix<C>): number {
  if (mat.nrows === 0 || mat.ncols === 0) return -1;
  let d = -1;
  for (let i = 0; i < mat.nrows; i++) {
    for (let j = 0; j < mat.ncols; j++) {
      const e = mat.get(i, j).degree();
      if (e > d) d = e;
    }
  }
  return d;
}

/**
 * Return the matrix of the (shifted) degrees in this matrix.
 *
 * The degree of the zero polynomial is `-1`; with shifts, the entry
 * corresponding to a zero polynomial is `min(shifts) - 1`.
 *
 * SageMath returns a matrix over `ZZ`; degrees are plain machine integers here
 * (as they already are in the return value of :func:`row_degrees`), so this
 * returns a nested array of numbers.
 *
 * @see Reference: sage/matrix/matrix_polynomial_dense.pyx:168-235
 */
export function degree_matrix<C extends FieldElement>(
  mat: PolynomialMatrix<C>,
  options?: { shifts?: number[] | null; row_wise?: boolean }
): number[][] {
  const shifts = options?.shifts ?? null;
  const row_wise = options?.row_wise ?? true;
  _check_shift_dimension(mat, shifts, row_wise);

  const out: number[][] = [];
  if (shifts === null) {
    for (let i = 0; i < mat.nrows; i++) {
      const row: number[] = [];
      for (let j = 0; j < mat.ncols; j++) row.push(mat.get(i, j).degree());
      out.push(row);
    }
    return out;
  }
  const zero_degree = Math.min(...shifts) - 1;
  for (let i = 0; i < mat.nrows; i++) {
    const row: number[] = [];
    for (let j = 0; j < mat.ncols; j++) {
      const e = mat.get(i, j);
      row.push(e.isZero() ? zero_degree : e.degree() + (row_wise ? shifts[j]! : shifts[i]!));
    }
    out.push(row);
  }
  return out;
}

/**
 * Return the constant coefficient of this matrix seen as a polynomial with
 * matrix coefficients; this is also this matrix evaluated at zero.
 *
 * @see Reference: sage/matrix/matrix_polynomial_dense.pyx:237-260
 */
export function constant_matrix<C extends FieldElement>(mat: PolynomialMatrix<C>): Matrix<C> {
  const K = _field(mat);
  const rows: C[][] = [];
  for (let i = 0; i < mat.nrows; i++) {
    const row: C[] = [];
    for (let j = 0; j < mat.ncols; j++) row.push(mat.get(i, j).getCoeff(0));
    rows.push(row);
  }
  return _fromRowsInferred(K, mat.nrows, mat.ncols, rows);
}

/**
 * Return whether this polynomial matrix is constant, that is, all its entries
 * are constant.
 *
 * @see Reference: sage/matrix/matrix_polynomial_dense.pyx:262-289
 */
export function is_constant<C extends FieldElement>(mat: PolynomialMatrix<C>): boolean {
  for (let j = 0; j < mat.ncols; j++) {
    for (let i = 0; i < mat.nrows; i++) {
      if (!mat.get(i, j).isConstant()) return false;
    }
  }
  return true;
}

/**
 * Return the constant matrix obtained by taking the coefficient of the entries
 * with degree specified by `d`.
 *
 * `d` is either an integer (same degree for all entries), or a list of length
 * the row dimension (``row_wise``, the default) or the column dimension.
 * Negative degrees give zero coefficients.
 *
 * @see Reference: sage/matrix/matrix_polynomial_dense.pyx:291-382
 */
export function coefficient_matrix<C extends FieldElement>(
  mat: PolynomialMatrix<C>,
  d: number | number[],
  options?: { row_wise?: boolean }
): Matrix<C> {
  const row_wise = options?.row_wise ?? true;
  const m = mat.nrows;
  const n = mat.ncols;
  const dd = _asList(d, row_wise ? m : n);

  if (row_wise && dd.length !== m) {
    throw new ValueError(
      'length of input degree list should be the row dimension of the input matrix'
    );
  } else if (!row_wise && dd.length !== n) {
    throw new ValueError(
      'length of input degree list should be the column dimension of the input matrix'
    );
  }

  const K = _field(mat);
  const rows: C[][] = [];
  for (let i = 0; i < m; i++) {
    const row: C[] = [];
    for (let j = 0; j < n; j++) {
      row.push(mat.get(i, j).getCoeff(row_wise ? dd[i]! : dd[j]!));
    }
    rows.push(row);
  }
  return _fromRows(K, m, n, rows);
}

/**
 * Return the matrix obtained by truncating all entries at the precisions
 * specified by `d` (an integer or a list).  The convention for univariate
 * polynomials is to take zero for the truncation for a negative `d`.
 *
 * @see Reference: sage/matrix/matrix_polynomial_dense.pyx:384-475
 */
export function truncate<C extends FieldElement>(
  mat: PolynomialMatrix<C>,
  d: number | number[],
  options?: { row_wise?: boolean }
): PolynomialMatrix<C> {
  const row_wise = options?.row_wise ?? true;
  const m = mat.nrows;
  const n = mat.ncols;
  const dd = _asList(d, row_wise ? m : n);

  if (row_wise && dd.length !== m) {
    throw new ValueError(
      'length of input precision list should be the row dimension of the input matrix'
    );
  } else if (!row_wise && dd.length !== n) {
    throw new ValueError(
      'length of input precision list should be the column dimension of the input matrix'
    );
  }

  const rows: Array<Array<Polynomial<C>>> = [];
  for (let i = 0; i < m; i++) {
    const row: Array<Polynomial<C>> = [];
    for (let j = 0; j < n; j++) row.push(mat.get(i, j).truncate(row_wise ? dd[i]! : dd[j]!));
    rows.push(row);
  }
  return _fromRows(mat.base_ring, m, n, rows);
}

/**
 * Return the matrix obtained by shifting all entries as specified by `d`
 * (an integer or a list): shifting by `d` means multiplying by the variable to
 * the power `d`; if `d` is negative then terms of negative degree after
 * shifting are discarded.
 *
 * @see Reference: sage/matrix/matrix_polynomial_dense.pyx:477-567
 */
export function shift<C extends FieldElement>(
  mat: PolynomialMatrix<C>,
  d: number | number[],
  options?: { row_wise?: boolean }
): PolynomialMatrix<C> {
  const row_wise = options?.row_wise ?? true;
  const m = mat.nrows;
  const n = mat.ncols;
  const dd = _asList(d, row_wise ? m : n);

  if (row_wise && dd.length !== m) {
    throw new ValueError(
      'length of input shift list should be the row dimension of the input matrix'
    );
  } else if (!row_wise && dd.length !== n) {
    throw new ValueError(
      'length of input shift list should be the column dimension of the input matrix'
    );
  }

  const rows: Array<Array<Polynomial<C>>> = [];
  for (let i = 0; i < m; i++) {
    const row: Array<Polynomial<C>> = [];
    for (let j = 0; j < n; j++) row.push(mat.get(i, j).shift(row_wise ? dd[i]! : dd[j]!));
    rows.push(row);
  }
  return _fromRows(mat.base_ring, m, n, rows);
}

/**
 * Return the matrix obtained by reversing all entries with respect to the
 * degree specified by ``degree``.
 *
 * If ``entry_wise`` is `true`, ``degree`` and ``row_wise`` are ignored and each
 * entry is reversed with respect to its own degree.  Otherwise ``degree`` may
 * be omitted (the matrix degree is used), an integer, or a list.
 *
 * @see Reference: sage/matrix/matrix_polynomial_dense.pyx:569-712
 */
export function reverse<C extends FieldElement>(
  mat: PolynomialMatrix<C>,
  options?: { degree?: number | number[] | null; row_wise?: boolean; entry_wise?: boolean }
): PolynomialMatrix<C> {
  const row_wise = options?.row_wise ?? true;
  const entry_wise = options?.entry_wise ?? false;
  const m = mat.nrows;
  const n = mat.ncols;

  if (entry_wise) {
    const rows: Array<Array<Polynomial<C>>> = [];
    for (let i = 0; i < m; i++) {
      const row: Array<Polynomial<C>> = [];
      for (let j = 0; j < n; j++) row.push(mat.get(i, j).reverse());
      rows.push(row);
    }
    return _fromRows(mat.base_ring, m, n, rows);
  }

  const deg: number | number[] = options?.degree ?? degree(mat);
  const dd = _asList(deg, row_wise ? m : n);

  if (row_wise && dd.length !== m) {
    throw new ValueError(
      'length of input degree list should be the row dimension of the input matrix'
    );
  } else if (!row_wise && dd.length !== n) {
    throw new ValueError(
      'length of input degree list should be the column dimension of the input matrix'
    );
  }

  // ``Polynomial.reverse`` rejects negative degrees
  // (polynomial_element.pyx:8137-8138); our Polynomial short-circuits on the
  // zero polynomial before that check, so validate here to keep the upstream
  // error for every entry.
  for (const e of dd) {
    if (e < 0) {
      throw new ValueError(`degree argument must be a nonnegative integer, got ${e}`);
    }
  }

  const rows: Array<Array<Polynomial<C>>> = [];
  for (let i = 0; i < m; i++) {
    const row: Array<Polynomial<C>> = [];
    for (let j = 0; j < n; j++) row.push(mat.get(i, j).reverse(row_wise ? dd[i]! : dd[j]!));
    rows.push(row);
  }
  return _fromRows(mat.base_ring, m, n, rows);
}

// ============================================================================
// Row / column degrees, leading matrix, leading positions
// ============================================================================

/**
 * Return the (shifted) row degrees of this matrix.
 *
 * Without shifts, `d_i = max_j deg(M[i,j])`, so `d_i = -1` for a zero row.
 * With shifts, `d_i = max_j (deg(M[i,j]) + s_j)`, and `d_i = min(shifts) - 1`
 * for a zero row.
 *
 * @see Reference: sage/matrix/matrix_polynomial_dense.pyx:1120-1197
 */
export function row_degrees<C extends FieldElement>(
  mat: PolynomialMatrix<C>,
  shifts?: number[] | null
): number[] {
  _check_shift_dimension(mat, shifts, true);
  if (mat.nrows === 0) return [];
  if (mat.ncols === 0) return new Array<number>(mat.nrows).fill(-1);
  if (shifts === undefined || shifts === null) {
    const out: number[] = [];
    for (let i = 0; i < mat.nrows; i++) {
      let d = Number.NEGATIVE_INFINITY;
      for (let j = 0; j < mat.ncols; j++) d = Math.max(d, mat.get(i, j).degree());
      out.push(d);
    }
    return out;
  }
  const zero_degree = Math.min(...shifts) - 1;
  const out: number[] = [];
  for (let i = 0; i < mat.nrows; i++) {
    let d = Number.NEGATIVE_INFINITY;
    for (let j = 0; j < mat.ncols; j++) {
      const e = mat.get(i, j);
      d = Math.max(d, e.isZero() ? zero_degree : e.degree() + shifts[j]!);
    }
    out.push(d);
  }
  return out;
}

/**
 * Return the (shifted) column degrees of this matrix.
 *
 * @see Reference: sage/matrix/matrix_polynomial_dense.pyx:1198-1266
 */
export function column_degrees<C extends FieldElement>(
  mat: PolynomialMatrix<C>,
  shifts?: number[] | null
): number[] {
  _check_shift_dimension(mat, shifts, false);
  if (mat.nrows === 0) return new Array<number>(mat.ncols).fill(-1);
  if (mat.ncols === 0) return [];
  if (shifts === undefined || shifts === null) {
    const out: number[] = [];
    for (let j = 0; j < mat.ncols; j++) {
      let d = Number.NEGATIVE_INFINITY;
      for (let i = 0; i < mat.nrows; i++) d = Math.max(d, mat.get(i, j).degree());
      out.push(d);
    }
    return out;
  }
  const zero_degree = Math.min(...shifts) - 1;
  const out: number[] = [];
  for (let j = 0; j < mat.ncols; j++) {
    let d = Number.NEGATIVE_INFINITY;
    for (let i = 0; i < mat.nrows; i++) {
      const e = mat.get(i, j);
      d = Math.max(d, e.isZero() ? zero_degree : e.degree() + shifts[i]!);
    }
    out.push(d);
  }
  return out;
}

/**
 * Return the (shifted) leading matrix of this matrix: the constant matrix
 * whose entry `(i,j)` is the coefficient of degree `d_i - s_j` of `M[i,j]`,
 * where `(d_i)` are the shifted row degrees (row-wise case).
 *
 * @see Reference: sage/matrix/matrix_polynomial_dense.pyx:1268-1362
 */
export function leading_matrix<C extends FieldElement>(
  mat: PolynomialMatrix<C>,
  options?: { shifts?: number[] | null; row_wise?: boolean }
): Matrix<C> {
  const shifts = options?.shifts ?? null;
  const row_wise = options?.row_wise ?? true;
  _check_shift_dimension(mat, shifts, row_wise);

  const K = _field(mat);
  const zero = K.zero();
  const rows: C[][] = [];

  if (row_wise) {
    const rdeg = row_degrees(mat, shifts);
    for (let i = 0; i < mat.nrows; i++) {
      const row: C[] = [];
      for (let j = 0; j < mat.ncols; j++) {
        const e = mat.get(i, j);
        const d = shifts === null ? e.degree() : e.degree() + shifts[j]!;
        row.push(d === rdeg[i] ? e.leading_coefficient() : zero);
      }
      rows.push(row);
    }
    return _fromRowsInferred(K, mat.nrows, mat.ncols, rows);
  }

  const cdeg = column_degrees(mat, shifts);
  for (let i = 0; i < mat.nrows; i++) {
    const row: C[] = [];
    for (let j = 0; j < mat.ncols; j++) {
      const e = mat.get(i, j);
      const d = shifts === null ? e.degree() : e.degree() + shifts[i]!;
      row.push(d === cdeg[j] ? e.leading_coefficient() : zero);
    }
    rows.push(row);
  }
  return _fromRowsInferred(K, mat.nrows, mat.ncols, rows);
}

/**
 * Assuming that this matrix is empty (`0 x n` or `m x 0`), return whether it is
 * in shifted Popov form.
 *
 * @see Reference: sage/matrix/matrix_polynomial_dense.pyx:1364-1415
 */
export function _is_empty_popov<C extends FieldElement>(
  mat: PolynomialMatrix<C>,
  row_wise: boolean = true,
  include_zero_vectors: boolean = true
): boolean {
  if (include_zero_vectors) return true;
  return row_wise ? mat.nrows === 0 : mat.ncols === 0;
}

/**
 * Return whether this matrix is in (shifted) reduced form: it has `k` nonzero
 * rows (resp. columns) with `k <= n` (resp. `<= m`) and its shifted leading
 * matrix has rank `k`.
 *
 * @see Reference: sage/matrix/matrix_polynomial_dense.pyx:1417-1495
 */
export function is_reduced<C extends FieldElement>(
  mat: PolynomialMatrix<C>,
  options?: { shifts?: number[] | null; row_wise?: boolean; include_zero_vectors?: boolean }
): boolean {
  const shifts = options?.shifts ?? null;
  const row_wise = options?.row_wise ?? true;
  const include_zero_vectors = options?.include_zero_vectors ?? true;
  _check_shift_dimension(mat, shifts, row_wise);
  if (mat.ncols === 0 || mat.nrows === 0) {
    return _is_empty_popov(mat, row_wise, include_zero_vectors);
  }

  let number_generators: number;
  if (include_zero_vectors) {
    number_generators = 0;
    if (row_wise) {
      for (let i = 0; i < mat.nrows; i++) {
        if (mat.row(i).some((p) => !p.isZero())) number_generators++;
      }
    } else {
      for (let j = 0; j < mat.ncols; j++) {
        if (mat.column(j).some((p) => !p.isZero())) number_generators++;
      }
    }
  } else {
    number_generators = row_wise ? mat.nrows : mat.ncols;
  }
  return number_generators === rank(leading_matrix(mat, { shifts, row_wise }));
}

/**
 * Return the (shifted) leading positions (pivot indices) and optionally the
 * (shifted) pivot degrees of this matrix.
 *
 * Working row-wise, the leading position of a row is the index `j` of the
 * rightmost nonzero entry `p_j` such that `deg(p_j) + s_j` equals the shifted
 * row degree of the row; its pivot degree is `deg(p_j)`.  For the zero row,
 * both are `-1`.
 *
 * @see Reference: sage/matrix/matrix_polynomial_dense.pyx:1497-1621
 */
export function leading_positions<C extends FieldElement>(
  mat: PolynomialMatrix<C>,
  options?: { shifts?: number[] | null; row_wise?: boolean; return_degree?: false }
): number[];
export function leading_positions<C extends FieldElement>(
  mat: PolynomialMatrix<C>,
  options: { shifts?: number[] | null; row_wise?: boolean; return_degree: true }
): [number[], number[]];
export function leading_positions<C extends FieldElement>(
  mat: PolynomialMatrix<C>,
  options?: { shifts?: number[] | null; row_wise?: boolean; return_degree?: boolean }
): number[] | [number[], number[]];
export function leading_positions<C extends FieldElement>(
  mat: PolynomialMatrix<C>,
  options?: { shifts?: number[] | null; row_wise?: boolean; return_degree?: boolean }
): number[] | [number[], number[]] {
  const shifts = options?.shifts ?? null;
  const row_wise = options?.row_wise ?? true;
  const return_degree = options?.return_degree ?? false;
  _check_shift_dimension(mat, shifts, row_wise);

  if (row_wise) {
    const rdeg = row_degrees(mat, shifts);
    const pivot_index: number[] = [];
    if (shifts === null) {
      for (let i = 0; i < mat.nrows; i++) {
        if (rdeg[i] === -1) {
          pivot_index.push(-1);
          continue;
        }
        let best = -1;
        for (let j = 0; j < mat.ncols; j++) {
          if (mat.get(i, j).degree() === rdeg[i]) best = j;
        }
        pivot_index.push(best);
      }
    } else {
      let zero_degree = -1;
      if (shifts.length > 0) zero_degree += Math.min(...shifts);
      for (let i = 0; i < mat.nrows; i++) {
        if (rdeg[i] === zero_degree) {
          pivot_index.push(-1);
          continue;
        }
        let best = -1;
        for (let j = 0; j < mat.ncols; j++) {
          const e = mat.get(i, j);
          if (!e.isZero() && e.degree() + shifts[j]! === rdeg[i]) best = j;
        }
        pivot_index.push(best);
      }
    }
    const pivot_degree = pivot_index.map((p, i) => (p === -1 ? -1 : mat.get(i, p).degree()));
    return return_degree ? [pivot_index, pivot_degree] : pivot_index;
  }

  // column-wise
  const cdeg = column_degrees(mat, shifts);
  const pivot_index: number[] = [];
  if (shifts === null) {
    for (let j = 0; j < mat.ncols; j++) {
      if (cdeg[j] === -1) {
        pivot_index.push(-1);
        continue;
      }
      let best = -1;
      for (let i = 0; i < mat.nrows; i++) {
        if (mat.get(i, j).degree() === cdeg[j]) best = i;
      }
      pivot_index.push(best);
    }
  } else {
    let zero_degree = -1;
    if (shifts.length > 0) zero_degree += Math.min(...shifts);
    for (let j = 0; j < mat.ncols; j++) {
      if (cdeg[j] === zero_degree) {
        pivot_index.push(-1);
        continue;
      }
      let best = -1;
      for (let i = 0; i < mat.nrows; i++) {
        const e = mat.get(i, j);
        if (!e.isZero() && e.degree() + shifts[i]! === cdeg[j]) best = i;
      }
      pivot_index.push(best);
    }
  }
  const pivot_degree = pivot_index.map((p, j) => (p === -1 ? -1 : mat.get(p, j).degree()));
  return return_degree ? [pivot_index, pivot_degree] : pivot_index;
}

// ============================================================================
// Form predicates
// ============================================================================

/**
 * Return whether this matrix is in (shifted) (ordered) weak Popov form: the
 * leading positions of its nonzero rows (resp. columns) are pairwise distinct,
 * and strictly increasing in the ``ordered`` variant.
 *
 * @see Reference: sage/matrix/matrix_polynomial_dense.pyx:1623-1752
 */
export function is_weak_popov<C extends FieldElement>(
  mat: PolynomialMatrix<C>,
  options?: {
    shifts?: number[] | null;
    row_wise?: boolean;
    ordered?: boolean;
    include_zero_vectors?: boolean;
  }
): boolean {
  const shifts = options?.shifts ?? null;
  const row_wise = options?.row_wise ?? true;
  const ordered = options?.ordered ?? false;
  const include_zero_vectors = options?.include_zero_vectors ?? true;
  _check_shift_dimension(mat, shifts, row_wise);
  if (mat.ncols === 0 || mat.nrows === 0) {
    return _is_empty_popov(mat, row_wise, include_zero_vectors);
  }

  let lpos = leading_positions(mat, { shifts, row_wise });
  // zero rows (resp. columns) get a leading position beyond every real one, so
  // that they sort to the bottom (resp. right)
  const pos_zero_vec = row_wise ? mat.ncols : mat.nrows;
  lpos = lpos.map((pos) => (pos >= 0 ? pos : pos_zero_vec + 1));
  if (!ordered) lpos = [...lpos].sort((a, b) => a - b);
  if (lpos[lpos.length - 1]! > pos_zero_vec && !include_zero_vectors) return false;
  for (let index = 0; index + 1 < lpos.length; index++) {
    const next_leading_position = lpos[index + 1]!;
    if (next_leading_position <= pos_zero_vec && next_leading_position <= lpos[index]!) {
      return false;
    }
  }
  return true;
}

/**
 * Return whether this matrix is in (shifted) Popov form: it is in ordered weak
 * Popov form (or weak Popov form if ``up_to_permutation``), each pivot entry is
 * monic, and each pivot has degree strictly larger than the other entries in
 * its column (resp. row).
 *
 * @see Reference: sage/matrix/matrix_polynomial_dense.pyx:1754-1891
 */
export function is_popov<C extends FieldElement>(
  mat: PolynomialMatrix<C>,
  options?: {
    shifts?: number[] | null;
    row_wise?: boolean;
    up_to_permutation?: boolean;
    include_zero_vectors?: boolean;
  }
): boolean {
  const shifts = options?.shifts ?? null;
  const row_wise = options?.row_wise ?? true;
  const up_to_permutation = options?.up_to_permutation ?? false;
  const include_zero_vectors = options?.include_zero_vectors ?? true;

  if (mat.ncols === 0 || mat.nrows === 0) {
    return _is_empty_popov(mat, row_wise, include_zero_vectors);
  }
  if (
    !is_weak_popov(mat, {
      shifts,
      row_wise,
      ordered: !up_to_permutation,
      include_zero_vectors,
    })
  ) {
    return false;
  }

  const [lpos, pivot_degree] = leading_positions(mat, {
    shifts,
    row_wise,
    return_degree: true,
  });
  for (let i = 0; i < lpos.length; i++) {
    const index = lpos[i]!;
    if (index >= 0) {
      if (row_wise) {
        if (!mat.get(i, index).is_monic()) return false;
        for (let k = 0; k < mat.nrows; k++) {
          if (k === i) continue;
          if (mat.get(k, index).degree() >= pivot_degree[i]!) return false;
        }
      } else {
        if (!mat.get(index, i).is_monic()) return false;
        for (let k = 0; k < mat.ncols; k++) {
          if (k === i) continue;
          if (mat.get(index, k).degree() >= pivot_degree[i]!) return false;
        }
      }
    }
  }
  return true;
}

/**
 * Return whether this matrix is in Hermite form: row echelon form with all
 * pivot entries monic and all entries above a pivot of degree less than the
 * pivot (row-wise, upper echelon).
 *
 * This is implemented exactly as upstream does, as a call to :func:`is_popov`
 * with the shift `((n-1)d, ..., 2d, d, 0)` where `d = degree(mat) + 1`.
 *
 * @see Reference: sage/matrix/matrix_polynomial_dense.pyx:1893-2008
 */
export function is_hermite<C extends FieldElement>(
  mat: PolynomialMatrix<C>,
  options?: { row_wise?: boolean; lower_echelon?: boolean; include_zero_vectors?: boolean }
): boolean {
  const row_wise = options?.row_wise ?? true;
  const lower_echelon = options?.lower_echelon ?? false;
  const include_zero_vectors = options?.include_zero_vectors ?? true;

  if (mat.ncols === 0 || mat.nrows === 0) {
    return _is_empty_popov(mat, row_wise, include_zero_vectors);
  }
  const d = degree(mat) + 1;
  const shift: number[] = row_wise
    ? Array.from({ length: mat.ncols }, (_, j) => j * d)
    : Array.from({ length: mat.nrows }, (_, j) => (mat.nrows - j) * d);
  if (!lower_echelon) shift.reverse();
  return is_popov(mat, { shifts: shift, row_wise, include_zero_vectors });
}

// ============================================================================
// Weak Popov form
// ============================================================================

/**
 * In-place Mulders-Storjohann transformation of ``M`` into weak Popov form.
 *
 * ``M`` is modified in place; the unimodular transformation is returned when
 * requested.  ``shifts``, when given, must be nonnegative.
 *
 * @see Reference: sage/matrix/matrix_polynomial_dense.pyx:2229-2343
 */
export function _weak_popov_form<C extends FieldElement>(
  M: PolynomialMatrix<C>,
  options?: { transformation?: boolean; shifts?: number[] | null }
): PolynomialMatrix<C> | undefined {
  const transformation = options?.transformation ?? false;
  const shifts = options?.shifts ?? null;

  const m = M.nrows;
  const n = M.ncols;
  const pR = _pring(M);

  let U: PolynomialMatrix<C> | undefined;
  if (transformation) {
    U = identity_matrix<Polynomial<C>>(M.base_ring, m);
  }

  // ``shifts`` is used as a Python truth value upstream: ``None`` and ``[]``
  // both disable the shift.
  const useShifts = shifts !== null && shifts.length > 0;

  /** Rightmost column reaching the (shifted) degree of row ``i``, and that degree. */
  const bestOfRow = (i: number): [number, number] => {
    let bestp = -1;
    let best = -1;
    for (let c = 0; c < n; c++) {
      let d = M.get(i, c).degree();
      if (useShifts && d >= 0) d += shifts![c]!;
      if (d >= best) {
        bestp = c;
        best = d;
      }
    }
    return [bestp, best];
  };

  const to_row: Array<Array<[number, number]>> = Array.from({ length: n }, () => []);
  const conflicts: number[] = [];
  for (let i = 0; i < m; i++) {
    const [bestp, best] = bestOfRow(i);
    if (best >= 0) {
      to_row[bestp]!.push([i, best]);
      if (to_row[bestp]!.length > 1) conflicts.push(bestp);
    }
  }

  while (conflicts.length > 0) {
    const c = conflicts.pop()!;
    const row = to_row[c]!;
    let [i, ideg] = row.pop()!;
    let [j, jdeg] = row.pop()!;

    if (jdeg > ideg) {
      [i, j] = [j, i];
      [ideg, jdeg] = [jdeg, ideg];
    }

    const coeff = _div(M.get(i, c).leading_coefficient(), M.get(j, c).leading_coefficient()).neg();
    const s = _monomial(pR, coeff as C, ideg - jdeg);

    for (let k = 0; k < n; k++) {
      M.set(i, k, M.get(i, k).add(s.mul(M.get(j, k))));
    }
    if (U !== undefined) {
      for (let k = 0; k < m; k++) {
        U.set(i, k, U.get(i, k).add(s.mul(U.get(j, k))));
      }
    }

    row.push([j, jdeg]);

    const [bestp, best] = bestOfRow(i);
    if (best >= 0) {
      to_row[bestp]!.push([i, best]);
      if (to_row[bestp]!.length > 1) conflicts.push(bestp);
    }
  }

  return U;
}

/**
 * Return a (shifted) (ordered) weak Popov form of this matrix.
 *
 * If the input matrix is `A`, a weak Popov form of `A` is any matrix `P` in
 * weak Popov form such that `UA = P` for some unimodular `U`; `U` is returned
 * as well when ``transformation`` is `true`.
 *
 * ALGORITHM: the Mulders-Storjohann algorithm of [MS2003], straightforwardly
 * extended to the case of shifted forms.
 *
 * @see Reference: sage/matrix/matrix_polynomial_dense.pyx:2010-2227
 */
export function weak_popov_form<C extends FieldElement>(
  mat: PolynomialMatrix<C>,
  options?: {
    transformation?: false;
    shifts?: number[] | null;
    row_wise?: boolean;
    ordered?: boolean;
    include_zero_vectors?: boolean;
  }
): PolynomialMatrix<C>;
export function weak_popov_form<C extends FieldElement>(
  mat: PolynomialMatrix<C>,
  options: {
    transformation: true;
    shifts?: number[] | null;
    row_wise?: boolean;
    ordered?: boolean;
    include_zero_vectors?: boolean;
  }
): [PolynomialMatrix<C>, PolynomialMatrix<C>];
export function weak_popov_form<C extends FieldElement>(
  mat: PolynomialMatrix<C>,
  options?: {
    transformation?: boolean;
    shifts?: number[] | null;
    row_wise?: boolean;
    ordered?: boolean;
    include_zero_vectors?: boolean;
  }
): PolynomialMatrix<C> | [PolynomialMatrix<C>, PolynomialMatrix<C>];
export function weak_popov_form<C extends FieldElement>(
  mat: PolynomialMatrix<C>,
  options?: {
    transformation?: boolean;
    shifts?: number[] | null;
    row_wise?: boolean;
    ordered?: boolean;
    include_zero_vectors?: boolean;
  }
): PolynomialMatrix<C> | [PolynomialMatrix<C>, PolynomialMatrix<C>] {
  const transformation = options?.transformation ?? false;
  const shifts = options?.shifts ?? null;
  const row_wise = options?.row_wise ?? true;
  const ordered = options?.ordered ?? false;
  const include_zero_vectors = options?.include_zero_vectors ?? true;

  // if column-wise, call the algorithm on the transpose
  if (!row_wise) {
    const W = weak_popov_form(mat.transpose(), {
      transformation,
      shifts,
      row_wise: true,
      ordered,
      include_zero_vectors,
    });
    if (transformation) {
      const [P, U] = W as [PolynomialMatrix<C>, PolynomialMatrix<C>];
      return [P.transpose(), U.transpose()];
    }
    return (W as PolynomialMatrix<C>).transpose();
  }

  const m = mat.nrows;
  _check_shift_dimension(mat, shifts, true);
  // make the shift nonnegative, as required by the main call ``_weak_popov_form``
  let nonnegative_shifts: number[] | null;
  if (shifts === null) {
    nonnegative_shifts = null;
  } else if (shifts.length === 0) {
    // upstream would raise from ``min([])``; this can only happen for a matrix
    // with zero columns, where the shift is irrelevant
    nonnegative_shifts = [];
  } else {
    const min_shifts = Math.min(...shifts);
    nonnegative_shifts = shifts.map((s) => s - min_shifts);
  }

  let M = mat.copy();
  let U = _weak_popov_form(M, { transformation, shifts: nonnegative_shifts });

  // move zero rows to the bottom of the matrix
  const zero_rows: number[] = [];
  const nonzero_rows: number[] = [];
  for (let i = 0; i < m; i++) {
    if (M.row(i).every((p) => p.isZero())) zero_rows.push(i + 1);
    else nonzero_rows.push(i + 1);
  }
  const gather = [...nonzero_rows, ...zero_rows];
  M = _permuteRows(M, gather);
  if (U !== undefined) U = _permuteRows(U, gather);

  // remove zero rows, if asked to (the corresponding rows of U are kept)
  const nnzr = m - zero_rows.length;
  if (!include_zero_vectors) M = _firstRows(M, nnzr);

  // order rows by increasing leading positions, if asked to
  if (ordered) {
    const lpos = leading_positions(_firstRows(M, nnzr), { shifts: nonnegative_shifts });
    if (include_zero_vectors) {
      for (let i = 0; i < m - nnzr; i++) lpos.push(m);
    }
    let row_permutation = _sortPermutation(lpos);
    M = _permuteRows(M, row_permutation);
    if (U !== undefined) {
      if (!include_zero_vectors) {
        for (let i = 0; i < m - nnzr; i++) lpos.push(m);
        row_permutation = _sortPermutation(lpos);
      }
      U = _permuteRows(U, row_permutation);
    }
  }

  return transformation ? [M, U!] : M;
}

/**
 * The one-based permutation sorting ``keys`` in increasing order, ties broken
 * by index -- SageMath's ``Permutation([elt[1] for elt in sorted(...)])``.
 */
function _sortPermutation(keys: number[]): number[] {
  return keys
    .map((k, i) => [k, i + 1] as [number, number])
    .sort((a, b) => (a[0] !== b[0] ? a[0] - b[0] : a[1] - b[1]))
    .map((t) => t[1]);
}

/**
 * Return a row reduced form of this matrix (a column reduced form when
 * ``row_wise`` is `false`).  Upstream this is a direct call to
 * :func:`weak_popov_form`.
 *
 * @see Reference: sage/matrix/matrix_polynomial_dense.pyx:2542-2669
 */
export function reduced_form<C extends FieldElement>(
  mat: PolynomialMatrix<C>,
  options?: {
    transformation?: false;
    shifts?: number[] | null;
    row_wise?: boolean;
    include_zero_vectors?: boolean;
  }
): PolynomialMatrix<C>;
export function reduced_form<C extends FieldElement>(
  mat: PolynomialMatrix<C>,
  options: {
    transformation: true;
    shifts?: number[] | null;
    row_wise?: boolean;
    include_zero_vectors?: boolean;
  }
): [PolynomialMatrix<C>, PolynomialMatrix<C>];
export function reduced_form<C extends FieldElement>(
  mat: PolynomialMatrix<C>,
  options?: {
    transformation?: boolean;
    shifts?: number[] | null;
    row_wise?: boolean;
    include_zero_vectors?: boolean;
  }
): PolynomialMatrix<C> | [PolynomialMatrix<C>, PolynomialMatrix<C>] {
  return weak_popov_form(mat, {
    transformation: options?.transformation ?? false,
    shifts: options?.shifts ?? null,
    row_wise: options?.row_wise ?? true,
    ordered: false,
    include_zero_vectors: options?.include_zero_vectors ?? true,
  });
}

// ============================================================================
// Popov form
// ============================================================================

/**
 * Return the (shifted) Popov form of this matrix: the unique matrix `P` in
 * (shifted) Popov form such that `UA = P` for some unimodular `U`.
 *
 * ALGORITHM: the Mulders-Storjohann algorithm of [MS2003] for transforming a
 * weak Popov form into Popov form, extended to shifted forms.
 *
 * @see Reference: sage/matrix/matrix_polynomial_dense.pyx:2345-2540
 */
export function popov_form<C extends FieldElement>(
  mat: PolynomialMatrix<C>,
  options?: {
    transformation?: false;
    shifts?: number[] | null;
    row_wise?: boolean;
    include_zero_vectors?: boolean;
  }
): PolynomialMatrix<C>;
export function popov_form<C extends FieldElement>(
  mat: PolynomialMatrix<C>,
  options: {
    transformation: true;
    shifts?: number[] | null;
    row_wise?: boolean;
    include_zero_vectors?: boolean;
  }
): [PolynomialMatrix<C>, PolynomialMatrix<C>];
export function popov_form<C extends FieldElement>(
  mat: PolynomialMatrix<C>,
  options?: {
    transformation?: boolean;
    shifts?: number[] | null;
    row_wise?: boolean;
    include_zero_vectors?: boolean;
  }
): PolynomialMatrix<C> | [PolynomialMatrix<C>, PolynomialMatrix<C>];
export function popov_form<C extends FieldElement>(
  mat: PolynomialMatrix<C>,
  options?: {
    transformation?: boolean;
    shifts?: number[] | null;
    row_wise?: boolean;
    include_zero_vectors?: boolean;
  }
): PolynomialMatrix<C> | [PolynomialMatrix<C>, PolynomialMatrix<C>] {
  const transformation = options?.transformation ?? false;
  const shifts = options?.shifts ?? null;
  const row_wise = options?.row_wise ?? true;
  const include_zero_vectors = options?.include_zero_vectors ?? true;

  // if column-wise, call the algorithm on the transpose
  if (!row_wise) {
    const res = popov_form(mat.transpose(), {
      transformation,
      shifts,
      row_wise: true,
      include_zero_vectors,
    });
    if (transformation) {
      const [P, U] = res as [PolynomialMatrix<C>, PolynomialMatrix<C>];
      return [P.transpose(), U.transpose()];
    }
    return (res as PolynomialMatrix<C>).transpose();
  }

  const nrows_zero = mat.nrows;

  // compute row-wise weak Popov form: non-ordered (rows are ordered below
  // anyway) and without zero rows (re-inserted later if asked to)
  const WP = weak_popov_form(mat, {
    transformation,
    shifts,
    row_wise: true,
    ordered: false,
    include_zero_vectors: false,
  });
  let P: PolynomialMatrix<C>;
  let UU: PolynomialMatrix<C> | undefined;
  if (transformation) {
    const pair = WP as [PolynomialMatrix<C>, PolynomialMatrix<C>];
    P = pair[0].copy();
    UU = pair[1];
  } else {
    P = (WP as PolynomialMatrix<C>).copy();
  }
  const m = P.nrows;
  let U: PolynomialMatrix<C> | undefined;
  if (transformation) U = _firstRows(UU!, m).copy();

  // compute leading positions and shifted row degrees
  let [lpos, rdeg] = leading_positions(P, { shifts, row_wise: true, return_degree: true });
  if (shifts !== null) {
    rdeg = rdeg.map((d, i) => d + shifts[lpos[i]!]!);
  }

  // 1/ transform P into ascending order ([MS2003, p394]; P has no zero rows)
  const sorted_rdeg_lpos = rdeg
    .map((d, i) => [d, lpos[i]!, i + 1] as [number, number, number])
    .sort((a, b) => (a[0] !== b[0] ? a[0] - b[0] : a[1] !== b[1] ? a[1] - b[1] : a[2] - b[2]));
  rdeg = sorted_rdeg_lpos.map((t) => t[0]);
  lpos = sorted_rdeg_lpos.map((t) => t[1]);
  let row_permutation = sorted_rdeg_lpos.map((t) => t[2]);
  P = _permuteRows(P, row_permutation);
  if (U !== undefined) U = _permuteRows(U, row_permutation);

  // 2/ ensure all pivots are monic
  for (let i = 0; i < m; i++) {
    const inv_lc = _inv(P.get(i, lpos[i]!).leading_coefficient());
    for (let j = 0; j < P.ncols; j++) P.set(i, j, P.get(i, j).scalar_mul(inv_lc));
    if (U !== undefined) {
      for (let j = 0; j < U.ncols; j++) U.set(i, j, U.get(i, j).scalar_mul(inv_lc));
    }
  }

  // 3/ reduce degrees as much as possible, row by row
  for (let i = 1; i < m; i++) {
    let delta = 0;
    while (delta >= 0) {
      // see [MS2003, Algo. PopovForm, p396]
      delta = -1;
      let j = -1;
      for (let k = 0; k < i; k++) {
        const diff = P.get(i, lpos[k]!).degree() - P.get(k, lpos[k]!).degree();
        if (diff > delta) {
          delta = diff;
          j = k;
        }
      }
      if (delta >= 0) {
        // the leading coefficient of P[j, lpos[j]] is 1
        const c = P.get(i, lpos[j]!).leading_coefficient().neg() as C;
        for (let k = 0; k < P.ncols; k++) {
          P.set(i, k, P.get(i, k).add(P.get(j, k).shift(delta).scalar_mul(c)));
        }
        if (U !== undefined) {
          for (let k = 0; k < U.ncols; k++) {
            U.set(i, k, U.get(i, k).add(U.get(j, k).shift(delta).scalar_mul(c)));
          }
        }
      }
    }
  }

  // 4/ transform so as to have increasing leading positions
  row_permutation = _sortPermutation(lpos);
  P = _permuteRows(P, row_permutation);
  if (U !== undefined) U = _permuteRows(U, row_permutation);

  // reinsert zero rows: in U in all cases, in P if asked to
  if (U !== undefined) U = _stack(U, _rowsFrom(UU!, m));
  if (include_zero_vectors) {
    P = _stack(P, new Matrix<Polynomial<C>>(mat.base_ring, nrows_zero - m, mat.ncols));
  }
  return transformation ? [P, U!] : P;
}

// ============================================================================
// Hermite form
// ============================================================================

/**
 * Transform ``A`` in place into Hermite normal form over a Euclidean domain,
 * optionally returning the transformation matrix.
 *
 * @see Reference: sage/matrix/matrix2.pyx:17113-17240 (``_hermite_form_euclidean``)
 */
function _hermite_form_euclidean<C extends FieldElement>(
  A: PolynomialMatrix<C>,
  transformation: boolean,
  normalization: (p: Polynomial<C>) => C
): PolynomialMatrix<C> | undefined {
  const m = A.nrows;
  const n = A.ncols;

  let i = 0;
  let j = 0;

  let U: PolynomialMatrix<C> | undefined;
  if (transformation) U = identity_matrix<Polynomial<C>>(A.base_ring, m);

  const pivot_cols: number[] = [];
  while (j < n) {
    let k = i;
    while (k < m && A.get(k, j).isZero()) k += 1; // first nonzero entry
    if (k < m) {
      let l = k + 1;
      while (l < m) {
        while (l < m && A.get(l, j).isZero()) l += 1; // nonzero entry below
        if (l >= m) break;

        const a = A.get(k, j);
        const b = A.get(l, j);
        const [d, p, q] = a.xgcd(b); // p*a + q*b = d = gcd(a,b)
        const e = a.quo_rem(d)[0];
        const f = b.quo_rem(d)[0];

        for (let c = j; c < n; c++) {
          const Akc = A.get(k, c);
          const Alc = A.get(l, c);
          A.set(k, c, p.mul(Akc).add(q.mul(Alc)));
          A.set(l, c, f.neg().mul(Akc).add(e.mul(Alc)));
        }
        if (U !== undefined) {
          for (let c = 0; c < m; c++) {
            const Ukc = U.get(k, c);
            const Ulc = U.get(l, c);
            U.set(k, c, p.mul(Ukc).add(q.mul(Ulc)));
            U.set(l, c, f.neg().mul(Ukc).add(e.mul(Ulc)));
          }
        }
      }
      if (i !== k) {
        _swapRows(A, i, k);
        if (U !== undefined) _swapRows(U, i, k);
      }
      pivot_cols.push(j);
      i += 1;
    }
    j += 1;
  }

  // reduce entries above pivots
  for (let r = 0; r < pivot_cols.length; r++) {
    const c0 = pivot_cols[r]!;
    let pivot = A.get(r, c0);

    // possibly normalize the pivot
    const coeff = normalization(pivot);
    for (let c = c0; c < n; c++) A.set(r, c, A.get(r, c).scalar_mul(coeff));
    if (U !== undefined) {
      for (let c = 0; c < m; c++) U.set(r, c, U.get(r, c).scalar_mul(coeff));
    }

    pivot = A.get(r, c0);
    for (let k = 0; k < r; k++) {
      const q = A.get(k, c0).quo_rem(pivot)[0].neg();
      if (!q.isZero()) {
        for (let c = c0; c < n; c++) A.set(k, c, A.get(k, c).add(q.mul(A.get(r, c))));
        if (U !== undefined) {
          for (let c = 0; c < m; c++) U.set(k, c, U.get(k, c).add(q.mul(U.get(r, c))));
        }
      }
    }
  }

  return U;
}

/** Swap two rows of a matrix in place. */
function _swapRows<R extends RingElement>(mat: Matrix<R>, r1: number, r2: number): void {
  for (let c = 0; c < mat.ncols; c++) {
    const t = mat.get(r1, c);
    mat.set(r1, c, mat.get(r2, c));
    mat.set(r2, c, t);
  }
}

/**
 * Return the Hermite form `H` of this matrix `A`: the unique matrix in Hermite
 * form such that `UA = H` for some unimodular `U`.
 *
 * @see Reference: sage/matrix/matrix_polynomial_dense.pyx:2671-2743
 */
export function hermite_form<C extends FieldElement>(
  mat: PolynomialMatrix<C>,
  options?: { include_zero_rows?: boolean; transformation?: false }
): PolynomialMatrix<C>;
export function hermite_form<C extends FieldElement>(
  mat: PolynomialMatrix<C>,
  options: { include_zero_rows?: boolean; transformation: true }
): [PolynomialMatrix<C>, PolynomialMatrix<C>];
export function hermite_form<C extends FieldElement>(
  mat: PolynomialMatrix<C>,
  options?: { include_zero_rows?: boolean; transformation?: boolean }
): PolynomialMatrix<C> | [PolynomialMatrix<C>, PolynomialMatrix<C>];
export function hermite_form<C extends FieldElement>(
  mat: PolynomialMatrix<C>,
  options?: { include_zero_rows?: boolean; transformation?: boolean }
): PolynomialMatrix<C> | [PolynomialMatrix<C>, PolynomialMatrix<C>] {
  const include_zero_rows = options?.include_zero_rows ?? true;
  const transformation = options?.transformation ?? false;

  let A = mat.copy();
  let U = _hermite_form_euclidean(A, transformation, (p) => _inv(p.leading_coefficient()));
  if (!include_zero_rows) {
    let i = A.nrows - 1;
    while (i >= 0 && A.row(i).every((p) => p.isZero())) i -= 1;
    A = _firstRows(A, i + 1);
    if (U !== undefined) U = _firstRows(U, i + 1);
  }
  return transformation ? [A, U!] : A;
}

// ============================================================================
// Minimal approximant bases
// ============================================================================

/** Evaluate every entry of a polynomial matrix at a base-field element. */
function _evaluateAt<C extends FieldElement>(mat: PolynomialMatrix<C>, x: C): Matrix<C> {
  const K = _field(mat);
  const rows: C[][] = [];
  for (let i = 0; i < mat.nrows; i++) {
    const row: C[] = [];
    for (let j = 0; j < mat.ncols; j++) row.push(mat.get(i, j).evaluate(x));
    rows.push(row);
  }
  return _fromRows(K, mat.nrows, mat.ncols, rows);
}

/**
 * Return whether ``mat`` is an approximant basis in ``shifts``-ordered weak
 * Popov form for the polynomial matrix ``pmat`` at order ``order``.
 *
 * ALGORITHM: verification that the matrix is formed by approximants is done via
 * a truncated matrix product; verification that the matrix is square,
 * nonsingular and in shifted weak Popov form is done via :func:`is_weak_popov`;
 * verification that the matrix generates the module of approximants is done via
 * the characterization in Theorem 2.1 of [GN2018].
 *
 * @see Reference: sage/matrix/matrix_polynomial_dense.pyx:3392-3590
 */
export function is_minimal_approximant_basis<C extends FieldElement>(
  mat: PolynomialMatrix<C>,
  pmat: PolynomialMatrix<C>,
  order: number | number[],
  options?: { shifts?: number[] | null; row_wise?: boolean; normal_form?: boolean }
): boolean {
  const row_wise = options?.row_wise ?? true;
  const normal_form = options?.normal_form ?? false;
  let shifts = options?.shifts ?? null;

  const m = pmat.nrows;
  const n = pmat.ncols;

  // set default shifts / check shifts dimension
  if (shifts === null) {
    shifts = new Array<number>(row_wise ? m : n).fill(0);
  } else if (row_wise && shifts.length !== m) {
    throw new ValueError('shifts length should be the row dimension of the input matrix');
  } else if (!row_wise && shifts.length !== n) {
    throw new ValueError('shifts length should be the column dimension of the input matrix');
  }

  // set default order / check order dimension
  const ord = _asList(order, row_wise ? n : m);
  if (row_wise && ord.length !== n) {
    throw new ValueError('order length should be the column dimension of the input matrix');
  } else if (!row_wise && ord.length !== m) {
    throw new ValueError('order length should be the row dimension of the input matrix');
  }

  // raise an error if mat does not have the right dimension
  if (row_wise && mat.ncols !== m) {
    throw new ValueError('column dimension should be the row dimension of the input matrix');
  } else if (!row_wise && mat.nrows !== n) {
    throw new ValueError('row dimension should be the column dimension of the input matrix');
  }

  // check square
  if (!mat.is_square()) return false;
  // check nonsingular and shifts-(ordered weak) Popov form
  if (
    normal_form &&
    !is_popov(mat, {
      shifts,
      row_wise,
      up_to_permutation: false,
      include_zero_vectors: false,
    })
  ) {
    return false;
  }
  if (
    !normal_form &&
    !is_weak_popov(mat, { shifts, row_wise, ordered: true, include_zero_vectors: false })
  ) {
    return false;
  }

  const pR = _pring(mat);
  const K = pR.base_ring;
  const X = pR.gen();

  if (row_wise) {
    // check that mat * pmat is 0 bmod x^order, and compute the certificate
    // matrix ``cert_mat``, the constant term of (mat * pmat) * x^(-order)
    const residual = mat.mul(pmat);
    if (!truncate(residual, ord, { row_wise: false }).is_zero()) return false;
    const cert_mat = coefficient_matrix(residual, ord, { row_wise: false });

    // 1/ determinant of mat should be a monomial c*x^d, d = sum of pivot degrees
    let d = 0;
    for (let i = 0; i < m; i++) d += mat.get(i, i).degree();
    const detOne = determinant(_evaluateAt(mat, K.one()));
    if (!determinant(mat).eq(X.pow(d).scalar_mul(detOne))) return false;
    // 2/ the m x (m+n) constant matrix [mat(0) | cert_mat] should have rank m
    if (rank(_augment(constant_matrix(mat), cert_mat)) < m) return false;
  } else {
    // check that pmat * mat is 0 bmod x^order
    const residual = pmat.mul(mat);
    if (!truncate(residual, ord).is_zero()) return false;
    const cert_mat = coefficient_matrix(residual, ord);

    let d = 0;
    for (let i = 0; i < n; i++) d += mat.get(i, i).degree();
    const detOne = determinant(_evaluateAt(mat, K.one()));
    if (!determinant(mat).eq(X.pow(d).scalar_mul(detOne))) return false;
    // the (m+n) x n constant matrix [mat(0).T | cert_mat.T].T should have rank n
    if (rank(_stack(constant_matrix(mat), cert_mat)) < n) return false;
  }

  return true;
}

/**
 * Return a ``shifts``-ordered weak Popov approximant basis for this polynomial
 * matrix at order ``order``, together with its ``shifts``-row degrees.
 *
 * The output basis is considered row-wise; the input dimensions are supposed to
 * be sound.
 *
 * ALGORITHM: inspired from the iterative algorithms described in [VBB1992] and
 * [BL1994].
 *
 * @see Reference: sage/matrix/matrix_polynomial_dense.pyx:3772-3919
 */
export function _approximant_basis_iterative<C extends FieldElement>(
  mat: PolynomialMatrix<C>,
  order: number[],
  shifts: number[]
): [PolynomialMatrix<C>, number[]] {
  const m = mat.nrows;
  const n = mat.ncols;

  // 'rem_order': the orders that remain to be dealt with
  // 'rem_index': indices of orders that remain to be dealt with
  const rem_order: number[] = [];
  const rem_index: number[] = [];
  for (let j = 0; j < n; j++) {
    if (order[j]! > 0) {
      rem_order.push(order[j]!);
      rem_index.push(j);
    }
  }

  // initialization of the residuals (= input mat, without the columns whose
  // order is zero) and of the approximant basis (= identity matrix)
  const appbas = identity_matrix<Polynomial<C>>(mat.base_ring, m);
  let residuals = _columnsOf(mat, rem_index);

  // throughout the algorithm, 'rdeg' is the shifts-row degree of 'appbas'
  const rdeg = [...shifts];

  while (rem_order.length > 0) {
    // choice for the next coefficient to be dealt with: first of the largest
    // entries in order
    const max_rem_order = Math.max(...rem_order);
    let j = 0;
    for (let ind = 0; ind < rem_order.length; ind++) {
      if (rem_order[ind] === max_rem_order) {
        j = ind;
        break;
      }
    }
    const d = order[rem_index[j]!]! - rem_order[j]!;

    // coefficient = the coefficient of degree d of column j of the residual
    const coefficient: C[] = [];
    for (let i = 0; i < m; i++) coefficient.push(residuals.get(i, j).getCoeff(d));

    // Lambda: rows [i] with nonzero coefficient[i]
    // pi: index of the first row with smallest shift, among those in Lambda
    const Lambda: number[] = [];
    let pi = -1;
    for (let i = 0; i < m; i++) {
      if (!coefficient[i]!.isZero()) {
        Lambda.push(i);
        if (pi < 0 || rdeg[i]! < rdeg[pi]!) pi = i;
      }
    }
    if (Lambda.length > 0) {
      // update all rows in Lambda--{pi}
      const idx = Lambda.indexOf(pi);
      Lambda.splice(idx, 1);
      for (const row of Lambda) {
        const scalar = _div(coefficient[row]!, coefficient[pi]!).neg() as C;
        for (let jj = 0; jj < m; jj++) {
          appbas.set(row, jj, appbas.get(row, jj).add(appbas.get(pi, jj).scalar_mul(scalar)));
        }
        for (let jj = 0; jj < residuals.ncols; jj++) {
          residuals.set(
            row,
            jj,
            residuals.get(row, jj).add(residuals.get(pi, jj).scalar_mul(scalar))
          );
        }
      }
      // update row pi: multiply by x
      rdeg[pi] = rdeg[pi]! + 1;
      for (let jj = 0; jj < m; jj++) appbas.set(pi, jj, appbas.get(pi, jj).shift(1));
      for (let jj = 0; jj < residuals.ncols; jj++) {
        residuals.set(pi, jj, residuals.get(pi, jj).shift(1));
      }
    }

    // decrement rem_order[j], unless there is no more work to do in this column
    if (rem_order[j] === 1) {
      const keep: number[] = [];
      for (let c = 0; c < residuals.ncols; c++) if (c !== j) keep.push(c);
      residuals = _columnsOf(residuals, keep);
      rem_order.splice(j, 1);
      rem_index.splice(j, 1);
    } else {
      rem_order[j] = rem_order[j]! - 1;
    }
  }
  return [appbas, rdeg];
}

/**
 * Return an approximant basis in ``shifts``-ordered weak Popov form for this
 * polynomial matrix at order ``order``.
 *
 * Working row-wise, an approximant basis for `F` at order `(d_0,...,d_{n-1})`
 * is a polynomial matrix whose rows form a basis of the module of polynomial
 * row vectors `p` such that the column `j` of `pF` has valuation at least
 * `d_j`.  If ``normal_form`` is `true` the output is furthermore in
 * ``shifts``-Popov form.
 *
 * ALGORITHM: inspired from the iterative algorithms of [VBB1992] and [BL1994];
 * the normal form relies on Lemmas 3.3 and 4.1 in [JNSV2016].
 *
 * @see Reference: sage/matrix/matrix_polynomial_dense.pyx:3593-3770
 */
export function minimal_approximant_basis<C extends FieldElement>(
  mat: PolynomialMatrix<C>,
  order: number | number[],
  options?: { shifts?: number[] | null; row_wise?: boolean; normal_form?: boolean }
): PolynomialMatrix<C> {
  const row_wise = options?.row_wise ?? true;
  const normal_form = options?.normal_form ?? false;
  let shifts = options?.shifts ?? null;

  const m = mat.nrows;
  const n = mat.ncols;

  // set default shifts / check shifts dimension
  if (shifts === null) {
    shifts = new Array<number>(row_wise ? m : n).fill(0);
  } else if (row_wise && shifts.length !== m) {
    throw new ValueError('shifts length should be the row dimension');
  } else if (!row_wise && shifts.length !== n) {
    throw new ValueError('shifts length should be the column dimension');
  }

  // set default order / check order dimension
  const ord = _asList(order, row_wise ? n : m);
  if (row_wise && ord.length !== n) {
    throw new ValueError('order length should be the column dimension');
  } else if (!row_wise && ord.length !== m) {
    throw new ValueError('order length should be the row dimension');
  }

  if (row_wise) {
    let [P, rdeg] = _approximant_basis_iterative(mat, ord, shifts);
    if (normal_form) {
      // the list "- pivot degree": -deg(P[i,i]) = shifts[i] - rdeg[i]
      const degree_shifts = shifts.map((s, i) => s - rdeg[i]!);
      [P, rdeg] = _approximant_basis_iterative(mat, ord, degree_shifts);
      const lmat = leading_matrix(P, { shifts: degree_shifts });
      P = _constTimesPoly(inverse(lmat), P);
    }
    return P;
  }

  const T = mat.transpose();
  let [P, rdeg] = _approximant_basis_iterative(T, ord, shifts);
  if (normal_form) {
    const degree_shifts = shifts.map((s, i) => s - rdeg[i]!);
    [P, rdeg] = _approximant_basis_iterative(T, ord, degree_shifts);
    let PT = P.transpose();
    const lmat = leading_matrix(PT, { shifts: degree_shifts, row_wise: false });
    PT = _polyTimesConst(PT, inverse(lmat));
    return PT;
  }
  return P.transpose();
}
