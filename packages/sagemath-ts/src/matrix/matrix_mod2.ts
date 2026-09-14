import { float as python_float, type FloatInput } from '../types/python_float.js';
import { str as matrix_str, type MatrixRepresentation } from './matrix0.js';
import { entries as matrix_entries } from './args.js';
import { PrimeField } from '../rings/finite_rings/finite_field_extension.js';
import { current_randstate } from '../misc/randstate.js';
import {
  mzd_pluq,
  _mzd_ple_russian,
  _mzd_pluq_russian,
  _mzd_ple_naive,
  _mzd_pluq_naive,
} from '@sagemath-ts/m4ri-ts';
import {
  mzd_echelonize,
  mzd_echelonize_m4ri,
  mzd_echelonize_pluq,
  mzd_ple,
} from '@sagemath-ts/m4ri-ts';
/**
 * @module sage/matrix/matrix_mod2
 * @description GF(2) matrix operations - important for cryptography
 *
 * Port of: sage/matrix/matrix_mod2_dense.pyx
 */

import {
  ArithmeticError,
  IndexError,
  NotImplementedError,
  ValueError,
  ZeroDivisionError,
} from '../errors.js';
import { RuntimeError, OverflowError } from '../errors.js';
import { type IntegerLike, toBigInt } from '../types/coercion.js';
import {
  mzd_init,
  mzd_mul,
  mzd_add,
  mzd_concat,
  mzd_transpose,
  mzd_inv_m4ri,
  mzd_solve_left,
  mzd_kernel_left_pluq,
  mzd_mul_m4rm,
  mzd_mul_naive,
  mzd_submatrix,
  mzd_density,
  type mzd_t,
} from '@sagemath-ts/m4ri-ts';
import { Integer } from '../rings/integer_ring.js';
import { Rational } from '../rings/rational.js';
import { create_RealNumber, type RealLiteral } from '../rings/real_mpfr.js';

type BinaryMatrixIndex = number | IntegerLike | Rational | boolean;
type BinarySubdivisionLine = BinaryMatrixIndex | string | null;
let binaryBaseRing: PrimeField | undefined;
function binaryEntry(value: unknown): number {
  binaryBaseRing ??= new PrimeField(2n);
  return Number(binaryBaseRing.__call__(value).value);
}

/** matrix0.pyx indexing and matrix_mod2_dense.pyx Py_ssize_t argument conversion.
 * Integer-valued JS numbers retain the port's existing array-index convention.
 * @see Deviation: Binary Matrix Indices and Slices
 */
function binaryMatrixIndex(value: BinaryMatrixIndex, kind: 'get' | 'set' | 'index'): bigint {
  let integer: bigint;
  if (typeof value === 'bigint') integer = value;
  else if (value instanceof Integer) integer = value.value;
  else if (value instanceof Rational) {
    // matrix0.__getitem__ assigns to C int (__int__, truncating); the
    // Py_ssize_t setters/row/submatrix use Rational.__index__ instead.
    if (kind !== 'get' && value.denominator !== 1n)
      throw new TypeError(`unable to convert rational ${value} to an integer`);
    integer = value.trunc();
  } else if (typeof value === 'boolean') integer = value ? 1n : 0n;
  else if (typeof value === 'number' && Number.isInteger(value)) integer = BigInt(value);
  else {
    if (kind === 'get') throw new TypeError('index must be an integer');
    if (kind === 'set')
      throw new TypeError(
        'index must be an integer or slice or a tuple/list of integers and slices'
      );
    const name =
      value == null
        ? 'NoneType'
        : typeof value === 'number'
          ? 'float'
          : typeof value === 'string'
            ? 'str'
            : 'object';
    throw new TypeError(`'${name}' object cannot be interpreted as an integer`);
  }
  if (integer < -(1n << 63n) || integer >= 1n << 63n)
    throw new OverflowError(
      `Python int too large to convert to C ${kind === 'get' ? 'long' : 'ssize_t'}`
    );
  if (kind === 'get' && (integer < -(1n << 31n) || integer >= 1n << 31n))
    throw new OverflowError('value too large to convert to int');
  return integer;
}

/**
 * A dense matrix over GF(2).
 *
 * This is a specialized matrix class for matrices over the field with two
 * elements, which is extremely important for cryptography and coding theory.
 *
 * Entries are stored as bits (0 or 1).
 *
 * @see Reference: sage/matrix/matrix_mod2_dense.pyx:Matrix_mod2_dense
 */
export class Matrix_mod2_dense {
  readonly nrows: number;
  readonly ncols: number;
  private _entries: number[][];
  private _immutable = false;
  private _subdivisions?: [bigint[], bigint[]];
  private _cache = new Map<string, unknown>();
  /** Original matrix0.pyx mutability flag; cached echelon forms are immutable. */
  set_immutable(): void {
    this._immutable = true;
  }
  is_immutable(): boolean {
    return this._immutable;
  }
  is_mutable(): boolean {
    return !this._immutable;
  }
  /** matrix2.pyx:subdivide; metadata changes invalidate caches before coercion.
   * @see Deviation: Binary Matrix Formatting and Subdivisions
   */
  subdivide(
    row_lines?: BinarySubdivisionLine | BinarySubdivisionLine[],
    col_lines?: BinarySubdivisionLine | BinarySubdivisionLine[]
  ): void {
    this.check_mutability();
    const rows =
      row_lines == null ? [] : Array.isArray(row_lines) ? row_lines.slice() : [row_lines];
    const cols =
      col_lines == null ? [] : Array.isArray(col_lines) ? col_lines.slice() : [col_lines];
    if (!rows.length && !cols.length) {
      this._subdivisions = undefined;
      return;
    }
    const sort = (left: unknown, right: unknown): number => {
      const name = (x: unknown) =>
        x == null
          ? 'NoneType'
          : typeof x === 'string'
            ? 'str'
            : typeof x === 'boolean'
              ? 'bool'
              : typeof x === 'number'
                ? 'int'
                : x instanceof Integer
                  ? 'sage.rings.integer.Integer'
                  : x instanceof Rational
                    ? 'sage.rings.rational.Rational'
                    : 'int';
      if (
        left == null ||
        right == null ||
        (typeof left === 'string') !== (typeof right === 'string')
      )
        throw new TypeError(
          `'<' not supported between instances of '${name(left)}' and '${name(right)}'`
        );
      if (typeof left === 'string' && typeof right === 'string')
        return left < right ? -1 : left > right ? 1 : 0;
      if (typeof left === 'number' && !Number.isFinite(left)) {
        if (Number.isNaN(left)) return 0;
        return left === right ? 0 : left < 0 ? -1 : 1;
      }
      if (typeof right === 'number' && !Number.isFinite(right)) {
        if (Number.isNaN(right)) return 0;
        return right < 0 ? 1 : -1;
      }
      const fraction = (x: unknown): [bigint, bigint] => {
        if (x instanceof Rational) return [x.numerator, x.denominator];
        if (x instanceof Integer) return [x.value, 1n];
        if (typeof x !== 'number' || Number.isInteger(x)) return [BigInt(x as bigint), 1n];
        // Exact IEEE input conversion for sorting; ZZ validation follows both sorts.
        const view = new DataView(new ArrayBuffer(8));
        view.setFloat64(0, x);
        const bits = view.getBigUint64(0),
          exponent = Number((bits >> 52n) & 2047n),
          shift = exponent ? exponent - 1075 : -1074;
        const mantissa =
          ((bits & ((1n << 52n) - 1n)) | (exponent ? 1n << 52n : 0n)) * (bits >> 63n ? -1n : 1n);
        return shift >= 0 ? [mantissa << BigInt(shift), 1n] : [mantissa, 1n << BigInt(-shift)];
      };
      const [a, b] = fraction(left),
        [c, d] = fraction(right),
        delta = a * d - c * b;
      return delta < 0n ? -1 : delta > 0n ? 1 : 0;
    };
    rows.sort(sort);
    cols.sort(sort);
    const r = rows.map((x) => new Integer(x).value),
      c = cols.map((x) => new Integer(x).value);
    this._subdivisions = [r, c];
  }
  /** Original visible subdivision lines, copied independently of stored metadata. */
  subdivisions(): [bigint[], bigint[]] {
    return this._subdivisions
      ? [this._subdivisions[0].slice(), this._subdivisions[1].slice()]
      : [[], []];
  }
  get_subdivisions(): [bigint[], bigint[]] {
    return this.subdivisions();
  }
  /** Clear computed properties without changing matrix mutability. */
  _clear_cache(): void {
    this._cache.clear();
  }
  private check_mutability(capital = false): void {
    if (this._immutable)
      throw new ValueError(
        (capital ? 'Matrix' : 'matrix') +
          ' is immutable; please change a copy instead (i.e., use copy(M) to change a copy of M).'
      );
    this._clear_cache();
  }
  private packed(): mzd_t {
    return mzd_init(
      this.nrows,
      this.ncols,
      this._entries.map((row) => row.reduce((v, b, j) => v | (BigInt(b) << BigInt(j)), 0n))
    );
  }
  private loadPacked(a: mzd_t): void {
    this._entries = a.rows.map((row) =>
      Array.from({ length: a.ncols }, (_, j) => Number((row >> BigInt(j)) & 1n))
    );
  }
  private currentPivots(): number[] {
    return this._entries.flatMap((row) => {
      const j = row.indexOf(1);
      return j < 0 ? [] : [j];
    });
  }

  /**
   * Create a GF(2) matrix.
   *
   * @param nrows - Number of rows
   * @param ncols - Number of columns
   * @param entries - Optional scalar, flat array or nested row array, converted through GF(2)
   * @see Deviation: Binary Constructor and Row Ownership
   */
  constructor(rows: BinaryMatrixIndex, cols: BinaryMatrixIndex, entries?: unknown) {
    // MatrixSpace.__classcall__: convert columns before rows, then validate rows first.
    const dimension = (value: BinaryMatrixIndex): bigint => {
      if (value instanceof Integer) return value.value;
      if (value instanceof Rational) return value.trunc();
      if (typeof value === 'number') {
        if (Number.isNaN(value)) throw new ValueError('cannot convert float NaN to integer');
        if (!Number.isFinite(value))
          throw new OverflowError('cannot convert float infinity to integer');
        return BigInt(Math.trunc(value));
      }
      return BigInt(value);
    };
    const c = dimension(cols),
      r = dimension(rows);
    if (r < 0n) throw new ArithmeticError('nrows must be nonnegative');
    if (c < 0n) throw new ArithmeticError('ncols must be nonnegative');
    if (r >= 1n << 63n || c >= 1n << 63n)
      throw new OverflowError('number of rows and columns may be at most 9223372036854775807');
    if (r >= 1n << 31n)
      throw new OverflowError(
        `matrices with ${r} rows over Finite Field of size 2 are not supported`
      );
    if (c >= 1n << 31n)
      throw new OverflowError(
        `matrices with ${c} columns over Finite Field of size 2 are not supported`
      );
    const nrows = (this.nrows = Number(r)),
      ncols = (this.ncols = Number(c));

    // Initialize entries array
    this._entries = [];
    for (let i = 0; i < nrows; i++) {
      this._entries.push([]);
      for (let j = 0; j < ncols; j++) {
        this._entries[i]!.push(0);
      }
    }
    for (const [i, j, value] of matrix_entries(nrows, ncols, entries, binaryEntry))
      this._entries[i]![j] = value;
  }

  /**
   * Get entry at (i, j).
   *
   * @param row - Row index
   * @param col - Column index
   * @returns 0 or 1
   * @see Deviation: Binary Matrix Indices and Slices
   */
  get(row: BinaryMatrixIndex, col: BinaryMatrixIndex): number {
    let i = Number(binaryMatrixIndex(row, 'get'));
    if (i < 0) i += this.nrows;
    if (i < 0 || i >= this.nrows) throw new IndexError('matrix index out of range');
    let j = Number(binaryMatrixIndex(col, 'get'));
    if (j < 0) j += this.ncols;
    if (j < 0 || j >= this.ncols) {
      throw new IndexError('matrix index out of range');
    }
    return this._entries[i]![j]!;
  }

  /**
   * Set entry at (i, j).
   *
   * @param row - Row index
   * @param col - Column index
   * @param value - A value coercible to GF(2)
   * @see Deviation: Binary Constructor and Row Ownership
   * @see Deviation: Binary Matrix Indices and Slices
   */
  set(row: BinaryMatrixIndex, col: BinaryMatrixIndex, value: unknown): void {
    this.check_mutability();
    let j = Number(binaryMatrixIndex(col, 'set'));
    if (j < 0) j += this.ncols;
    if (j < 0 || j >= this.ncols) throw new IndexError('index out of range');
    let i = Number(binaryMatrixIndex(row, 'set'));
    if (i < 0) i += this.nrows;
    if (i < 0 || i >= this.nrows) {
      throw new IndexError('index out of range');
    }
    this._entries[i]![j] = binaryEntry(value);
  }

  /**
   * Return string representation with custom mapping.
   *
   * @param rep_mapping - Mapping for values
   * @param zero - String for zero
   * @param plus_one - String for one
   * @param minus_one - String for minus one (same as one in GF(2))
   * @returns String representation
   * @see Reference: sage/matrix/matrix_mod2_dense.pyx:str
   * @see Deviation: Binary Matrix Formatting and Subdivisions
   */
  str(
    rep_mapping?: MatrixRepresentation | null,
    zero?: string | null,
    plus_one?: string | null,
    minus_one?: string | null
  ): string {
    const [rowDivisions, colDivisions] = this.subdivisions();
    let generic = rep_mapping != null || zero != null || plus_one != null;
    let previous = 0n;
    for (const col of colDivisions) {
      if (col === previous || col <= 0n || col >= BigInt(this.ncols)) generic = true;
      previous = col;
    }
    return matrix_str(
      {
        nrows: this.nrows,
        ncols: this.ncols,
        entries: () => (generic ? this.list() : this._entries.flat()),
        rowDivisions,
        colDivisions,
        binaryFast: !generic,
      },
      rep_mapping,
      zero,
      plus_one
    );
  }

  /**
   * Return a specific row.
   *
   * @param index - Row index
   * @param from_list - Return the shared immutable cached row when true
   * @see Deviation: Binary Constructor and Row Ownership
   * @returns The row as an array of 0s and 1s
   * @see Reference: sage/matrix/matrix_mod2_dense.pyx:row
   * @see Deviation: Binary Matrix Indices and Slices
   */
  row(index: BinaryMatrixIndex, from_list?: boolean): number[] {
    let i = Number(binaryMatrixIndex(index, 'index'));
    if (this.nrows === 0) throw new IndexError('matrix has no rows');
    if (i < 0) {
      i = i + this.nrows;
    }
    if (i < 0 || i >= this.nrows) {
      throw new IndexError('row index out of range');
    }
    if (from_list) {
      // matrix1.rows/dense_rows cache immutable vectors; default row extraction
      // remains independent and mutable. Matrix mutations clear both caches.
      let rows = this._cache.get('rows') as number[][] | undefined;
      if (rows === undefined) {
        rows = Array.from(
          { length: this.nrows },
          (_, j) =>
            new Proxy(Object.freeze(this.row(j)) as unknown as number[], {
              set() {
                throw new ValueError(
                  'vector is immutable; please change a copy instead (use copy())'
                );
              },
            })
        );
        this._cache.set('dense_rows', rows);
        this._cache.set('rows', rows);
      }
      return rows[i]!;
    }
    return this.submatrix(i, 0, 1, this.ncols)._entries[0]!;
  }

  /**
   * Return all columns.
   *
   * @param copy - Whether to copy
   * @returns List of columns (each column is an array of 0s and 1s)
   * @see Reference: sage/matrix/matrix_mod2_dense.pyx:columns
   * @see Deviation: Binary Basic Operations and Column Cache
   */
  columns(copy = true): number[][] {
    let columns = this._cache.get('columns') as number[][] | undefined;
    if (columns === undefined) {
      columns = this.transpose()._entries.map(
        (row) =>
          new Proxy(Object.freeze(row) as unknown as number[], {
            set() {
              throw new ValueError(
                'vector is immutable; please change a copy instead (use copy())'
              );
            },
          })
      );
      this._cache.set('columns', columns);
    }
    return copy ? columns.slice() : columns;
  }

  /**
   * Matrix addition (XOR in GF(2)).
   *
   * @param other - Matrix to add
   * @returns Sum of matrices
   */
  add(other: Matrix_mod2_dense): Matrix_mod2_dense {
    return this.binarySum(other, '+');
  }

  private binarySum(other: Matrix_mod2_dense, operation: '+' | '-'): Matrix_mod2_dense {
    if (this.nrows !== other.nrows || this.ncols !== other.ncols) {
      const parent = (A: Matrix_mod2_dense) =>
        `Full MatrixSpace of ${A.nrows} by ${A.ncols} dense matrices over Finite Field of size 2`;
      throw new TypeError(
        `unsupported operand parent(s) for ${operation}: '${parent(this)}' and '${parent(other)}'`
      );
    }
    const result = new Matrix_mod2_dense(this.nrows, this.ncols);
    if (this.nrows && this.ncols) result.loadPacked(mzd_add(this.packed(), other.packed()));
    return result;
  }

  /**
   * Matrix subtraction (same as addition in GF(2)).
   *
   * @param other - Matrix to subtract
   * @returns Difference of matrices
   */
  sub(other: Matrix_mod2_dense): Matrix_mod2_dense {
    return this.binarySum(other, '-');
  }

  /**
   * Classical matrix multiplication.
   *
   * @param right - The right operand
   * @returns The product
   * @see Reference: sage/matrix/matrix_mod2_dense.pyx:_multiply_classical
   * @see Deviation: Specialized Binary Matrix Multiplication
   */
  _multiply_classical(right: Matrix_mod2_dense): Matrix_mod2_dense {
    if (this.nrows === 0 || this.ncols === 0 || right.ncols === 0)
      return new Matrix_mod2_dense(this.nrows, right.ncols);
    if (this.ncols !== right.nrows)
      throw new ArithmeticError(
        `cannot multiply ${this.nrows}x${this.ncols} matrix by ${right.nrows}x${right.ncols} matrix`
      );
    return this._m4ri_product(right, 2, 0);
  }

  /**
   * Matrix multiplication through M4RI Strassen-Winograd and M4RM.
   * @see Deviation: Specialized Binary Matrix Multiplication
   *
   * @param other - Matrix to multiply
   * @returns Product of matrices
   */
  mul(other: Matrix_mod2_dense): Matrix_mod2_dense {
    if (this.ncols !== other.nrows)
      throw new TypeError(
        `unsupported operand parent(s) for *: 'Full MatrixSpace of ${this.nrows} by ${this.ncols} dense matrices over Finite Field of size 2' and 'Full MatrixSpace of ${other.nrows} by ${other.ncols} dense matrices over Finite Field of size 2'`
      );
    return this._multiply_strassen(other, 0n);
  }

  /** Original M4RM entry point. @see Deviation: Specialized Binary Matrix Multiplication */
  _multiply_m4rm(right: Matrix_mod2_dense, k: IntegerLike): Matrix_mod2_dense {
    const value = toBigInt(k);
    if (value < -(1n << 63n) || value >= 1n << 63n)
      throw new OverflowError('Python int too large to convert to C long');
    if (value < -(1n << 31n) || value >= 1n << 31n)
      throw new OverflowError('value too large to convert to int');
    if (this.ncols !== right.nrows) throw new ArithmeticError('left ncols must match right nrows');
    return this._m4ri_product(right, 1, Number(value));
  }

  /** Original Strassen entry point. @see Deviation: Specialized Binary Matrix Multiplication */
  _multiply_strassen(right: Matrix_mod2_dense, cutoff: IntegerLike): Matrix_mod2_dense {
    const value = toBigInt(cutoff);
    if (value < -(1n << 63n) || value >= 1n << 63n)
      throw new OverflowError('Python int too large to convert to C long');
    if (value < -(1n << 31n) || value >= 1n << 31n)
      throw new OverflowError('value too large to convert to int');
    if (this.ncols !== right.nrows) throw new ArithmeticError('left ncols must match right nrows');
    return this._m4ri_product(right, 0, Number(value));
  }

  private _m4ri_product(
    right: Matrix_mod2_dense,
    method: number,
    parameter: number
  ): Matrix_mod2_dense {
    if (this.nrows === 0 || this.ncols === 0 || right.ncols === 0)
      return new Matrix_mod2_dense(this.nrows, right.ncols);
    const pack = (matrix: Matrix_mod2_dense): mzd_t =>
      mzd_init(
        matrix.nrows,
        matrix.ncols,
        matrix._entries.map((row) => {
          const words: string[] = [];
          for (let j = 0; j < row.length; j += 64) {
            let word = 0n;
            for (let t = 0; t < 64 && j + t < row.length; t++)
              word |= BigInt(row[j + t]!) << BigInt(t);
            words.push(word.toString(16).padStart(16, '0'));
          }
          return BigInt('0x' + words.reverse().join(''));
        })
      );
    const a = pack(this),
      b = this === right ? a : pack(right);
    let c: mzd_t;
    try {
      c =
        method === 0
          ? mzd_mul(a, b, parameter)
          : method === 1
            ? mzd_mul_m4rm(a, b, parameter)
            : mzd_mul_naive(a, b);
    } catch (error) {
      if (error instanceof Error && error.name === 'RuntimeError')
        throw new RuntimeError(error.message);
      throw error;
    }
    const width = Math.ceil(c.ncols / 64);
    const entries = c.rows.map((row) => {
      const hex = row.toString(16).padStart(width * 16, '0'),
        out: number[] = [];
      for (let j = 0; j < width; j++) {
        const end = hex.length - j * 16,
          word = BigInt('0x' + hex.slice(end - 16, end));
        for (let t = 0; t < 64 && out.length < c.ncols; t++)
          out.push(Number((word >> BigInt(t)) & 1n));
      }
      return out;
    });
    return new Matrix_mod2_dense(c.nrows, c.ncols, entries);
  }

  /**
   * Return the negation (same as self in GF(2)).
   *
   * @returns The negation (a copy of self)
   * @see Reference: sage/matrix/matrix_mod2_dense.pyx:__neg__
   */
  neg(): Matrix_mod2_dense {
    return this.copy();
  }

  /**
   * Return the inverse.
   *
   * @returns The inverse
   * @throws {ZeroDivisionError} If the matrix is not invertible
   * @see Reference: sage/matrix/matrix_mod2_dense.pyx:__invert__
   * @see Deviation: Native Binary Inverse Adapter
   */
  inverse(): Matrix_mod2_dense {
    if (this.nrows !== this.ncols) {
      throw new ArithmeticError('self must be a square matrix');
    }

    const n = this.nrows;

    if (n === 0) {
      return this.copy();
    }

    if (this.rank() !== n) {
      throw new ZeroDivisionError('Matrix does not have full rank.');
    }

    const packed = mzd_inv_m4ri(this.packed());
    const result = new Matrix_mod2_dense(n, n);
    result.loadPacked(packed);
    return result;
  }

  /**
   * Return a copy.
   *
   * @returns A copy
   * @see Reference: sage/matrix/matrix_mod2_dense.pyx:__copy__
   */
  copy(): Matrix_mod2_dense {
    const result = new Matrix_mod2_dense(this.nrows, this.ncols);
    for (let i = 0; i < this.nrows; i++) {
      for (let j = 0; j < this.ncols; j++) {
        result._entries[i]![j] = this._entries[i]![j]!;
      }
    }
    if (this._subdivisions) result.subdivide(...this.subdivisions());
    return result;
  }

  /**
   * Return list of entries.
   *
   * @returns List of entries in row-major order
   * @see Reference: sage/matrix/matrix_mod2_dense.pyx:_list
   */
  list(): number[] {
    const result: number[] = [];
    for (let i = 0; i < this.nrows; i++) {
      for (let j = 0; j < this.ncols; j++) {
        result.push(this._entries[i]![j]!);
      }
    }
    return result;
  }

  /**
   * Put matrix in echelon form.
   *
   * @param algorithm - 'heuristic', 'pluq', 'm4ri', or 'classical'
   * @param cutoff - Cutoff for algorithm
   * @param reduced - Whether to compute reduced form (default: true)
   * @see Reference: sage/matrix/matrix_mod2_dense.pyx:echelonize
   */
  /** @see Deviation: Binary Elimination and Cached Forms */
  echelonize(
    algorithm = 'heuristic',
    cutoff?: number,
    reduced = true,
    options?: { k?: IntegerLike }
  ): void | this {
    if (!this.nrows || !this.ncols) {
      this._cache.set('in_echelon_form', true);
      this._cache.set('rank', 0);
      this._cache.set('pivots', []);
      return this;
    }
    if (this._cache.has('in_echelon_form')) return;
    if (algorithm === 'linbox') throw new NotImplementedError('');
    if (!['heuristic', 'm4ri', 'pluq', 'classical'].includes(algorithm))
      throw new ValueError(`no algorithm '${algorithm}'`);
    this.check_mutability();
    if (algorithm === 'classical') {
      let pivot = 0;
      for (let c = 0; c < this.ncols && pivot < this.nrows; c++) {
        let found = pivot;
        while (found < this.nrows && !this._entries[found]![c]) found++;
        if (found === this.nrows) continue;
        [this._entries[pivot], this._entries[found]] = [
          this._entries[found]!,
          this._entries[pivot]!,
        ];
        for (let i = 0; i < this.nrows; i++)
          if (i !== pivot && this._entries[i]![c])
            for (let j = c; j < this.ncols; j++) this._entries[i]![j] ^= this._entries[pivot]![j]!;
        pivot++;
      }
      this._cache.set('echelon_form', this);
    } else {
      let k = 0;
      if (algorithm === 'm4ri' && options?.k !== undefined) {
        const value = toBigInt(options.k);
        if (value < -(1n << 63n) || value >= 1n << 63n)
          throw new OverflowError('Python int too large to convert to C long');
        if (value < -(1n << 31n) || value >= 1n << 31n)
          throw new OverflowError('value too large to convert to int');
        k = Number(value);
        if (k < 1 || k > 16) throw new RuntimeError('k must be between 1 and 16');
      }
      const a = this.packed(),
        result =
          algorithm === 'm4ri'
            ? mzd_echelonize_m4ri(a, reduced, k)
            : algorithm === 'pluq'
              ? mzd_echelonize_pluq(a, reduced)
              : mzd_echelonize(a, reduced);
      this.loadPacked(result.matrix);
      this._cache.set('rank', result.rank);
    }
    this._cache.set('in_echelon_form', true);
    this._cache.set('pivots', this.currentPivots());
  }

  /**
   * Return the echelon form (without modifying self).
   *
   * @param algorithm - Algorithm to use
   * @param cutoff - Cutoff parameter
   * @param reduced - Whether to compute reduced form
   * @returns The echelon form
   */
  /** @see Deviation: Binary Elimination and Cached Forms */
  echelon_form(
    algorithm = 'default',
    cutoff?: number,
    reduced = true,
    options?: { k?: IntegerLike }
  ): Matrix_mod2_dense {
    const cached = this._cache.get('echelon_form');
    if (cached) return cached as Matrix_mod2_dense;
    const result = this.copy();
    result.echelonize(algorithm === 'default' ? 'heuristic' : algorithm, cutoff, reduced, options);
    result.set_immutable();
    this._cache.set('echelon_form', result);
    this._cache.set('pivots', result.pivots());
    return result;
  }

  /**
   * Return pivot columns.
   *
   * @returns List of pivot column indices
   * @see Reference: sage/matrix/matrix_mod2_dense.pyx:_pivots
   */
  pivots(): number[] {
    const cached = this._cache.get('pivots');
    if (cached) return [...(cached as number[])];
    const result = this.echelon_form().pivots();
    this._cache.set('pivots', result);
    return [...result];
  }

  /**
   * Randomize entries.
   *
   * @param density - Proportion selected for updates (default: 1.0 for uniform random)
   * @see Deviation: Binary Seeded Randomization and Factories
   * @see Deviation: Binary Density Input Conversion
   * @param nonzero - Whether all entries should be nonzero (all ones)
   * @see Reference: sage/matrix/matrix_mod2_dense.pyx:randomize
   */
  randomize(density?: FloatInput, nonzero?: boolean): void {
    if (!this.nrows || !this.ncols) return;
    let d = density === undefined ? 1.0 : python_float(density);
    if (d <= 0) return;
    if (d > 1) d = 1;
    this.check_mutability();
    const state = current_randstate();
    if (nonzero) {
      for (let i = 0; i < this.nrows; i++) {
        for (let j = 0; j < this.ncols; j++) {
          if (state.c_rand_double() <= d) {
            this._entries[i]![j] = 1;
          }
        }
      }
    } else if (d === 1) {
      // Original M4RI words consume two 32-bit GMP draws, including padding.
      // m4ri_swap_bits places the high draw's most-significant bit in column 0.
      for (let i = 0; i < this.nrows; i++) {
        for (let start = 0; start < this.ncols; start += 64) {
          const low = state.random_bits(32),
            high = state.random_bits(32);
          const word = low | (high << 32n);
          for (let j = 0; j < Math.min(64, this.ncols - start); j++) {
            this._entries[i]![start + j] = Number((word >> BigInt(63 - j)) & 1n);
          }
        }
      }
    } else {
      // The sparse path chooses positions with replacement, not Bernoulli trials.
      if (Number.isNaN(d)) throw new ValueError('cannot convert float NaN to integer');
      const perRow = Math.trunc(d * this.ncols);
      for (let i = 0; i < this.nrows; i++) {
        for (let j = 0; j < perRow; j++) {
          const column = state.c_random() % this.ncols;
          this._entries[i]![column] = state.c_random() % 2;
        }
      }
    }
  }

  /**
   * Return the determinant.
   *
   * @returns The determinant (0 or 1)
   * @see Reference: sage/matrix/matrix_mod2_dense.pyx:determinant
   */
  determinant(): number {
    if (this.nrows !== this.ncols) {
      throw new ValueError('self must be a square matrix');
    }

    // Over GF(2), det = 1 iff matrix has full rank
    return this.rank() === this.nrows ? 1 : 0;
  }

  /**
   * Return the transpose.
   *
   * @returns The transpose
   * @see Reference: sage/matrix/matrix_mod2_dense.pyx:transpose
   * @see Deviation: Binary Basic Operations and Column Cache
   */
  transpose(): Matrix_mod2_dense {
    const result = new Matrix_mod2_dense(this.ncols, this.nrows);
    if (this.nrows && this.ncols) result.loadPacked(mzd_transpose(this.packed()));
    if (this._subdivisions) {
      const [r, c] = this.subdivisions();
      result.subdivide(c, r);
    }
    return result;
  }

  /**
   * Augment with another matrix.
   *
   * @param right - Matrix to augment with
   * @param subdivide - Whether to subdivide
   * @returns Augmented matrix
   * @see Reference: sage/matrix/matrix_mod2_dense.pyx:augment
   * @see Deviation: Binary Basic Operations and Column Cache
   */
  augment(right: Matrix_mod2_dense, subdivide?: boolean): Matrix_mod2_dense {
    if (this.nrows !== right.nrows) throw new TypeError('Both numbers of rows must match.');
    if (!this.ncols) return right.copy();
    if (!right.ncols) return this.copy();
    const result = new Matrix_mod2_dense(this.nrows, this.ncols + right.ncols);
    if (this.nrows) result.loadPacked(mzd_concat(this.packed(), right.packed()));
    if (subdivide && this.nrows) {
      const [lr, lc] = this.subdivisions(),
        [rr, rc] = right.subdivisions();
      result.subdivide(lr.length === rr.length && lr.every((v, i) => v === rr[i]) ? lr : null, [
        ...lc,
        BigInt(this.ncols),
        ...rc.map((v) => v + BigInt(this.ncols)),
      ]);
    }
    return result;
  }

  /**
   * Return a submatrix.
   *
   * @param row - Starting row (default: 0)
   * @param col - Starting column (default: 0)
   * @param nrows - Number of rows (default: remaining rows)
   * @param ncols - Number of columns (default: remaining columns)
   * @returns The submatrix
   * @see Reference: sage/matrix/matrix_mod2_dense.pyx:submatrix
   * @see Deviation: Binary Matrix Indices and Slices
   */
  submatrix(
    row: BinaryMatrixIndex = 0,
    col: BinaryMatrixIndex = 0,
    nrows: BinaryMatrixIndex = -1,
    ncols: BinaryMatrixIndex = -1
  ): Matrix_mod2_dense {
    const startRow = binaryMatrixIndex(row, 'index');
    const startCol = binaryMatrixIndex(col, 'index');
    let numRows = binaryMatrixIndex(nrows, 'index');
    let numCols = binaryMatrixIndex(ncols, 'index');
    if (numRows < 0n) numRows = BigInt(this.nrows) > startRow ? BigInt(this.nrows) - startRow : 0n;
    if (numCols < 0n) numCols = BigInt(this.ncols) > startCol ? BigInt(this.ncols) - startCol : 0n;
    const highr = startRow + numRows,
      highc = startCol + numCols;
    if (startRow < 0) throw new TypeError(`Expected row >= 0, but got ${startRow} instead.`);
    if (startCol < 0) throw new TypeError(`Expected col >= 0, but got ${startCol} instead.`);
    if (highc > this.ncols)
      throw new TypeError(
        `Expected highc <= self.ncols(), but got ${highc} > ${this.ncols} instead.`
      );
    if (highr > this.nrows)
      throw new TypeError(
        `Expected highr <= self.nrows(), but got ${highr} > ${this.nrows} instead.`
      );

    const r0 = Number(startRow),
      c0 = Number(startCol),
      r = Number(numRows),
      c = Number(numCols);
    if (r === 0 || c === 0) return new Matrix_mod2_dense(r, c);
    // Pack only source words intersecting the requested slice. Packing the
    // whole source would make tiny slices linear in the whole matrix size.
    const firstWord = c0 - (c0 % 64),
      lastWord = Math.min(this.ncols, Math.ceil(Number(highc) / 64) * 64);
    const source = mzd_init(
      r,
      lastWord - firstWord,
      this._entries
        .slice(r0, r0 + r)
        .map((row) => BigInt('0b' + row.slice(firstWord, lastWord).reverse().join('')))
    );
    const result = mzd_submatrix(source, 0, c0 - firstWord, r, c0 - firstWord + c);
    return new Matrix_mod2_dense(
      r,
      c,
      result.rows.map((row) => Array.from(row.toString(2).padStart(c, '0')).reverse().map(Number))
    );
  }

  /**
   * Return the density (proportion of nonzero entries).
   *
   * @param approx - Whether to return approximate value
   * @returns The density of ones
   * @see Reference: sage/matrix/matrix_mod2_dense.pyx:density
   * @see Deviation: Binary Matrix Density
   */
  density(approx?: false): Rational | bigint;
  density(approx: true): RealLiteral;
  density(approx?: boolean): Rational | bigint | RealLiteral;
  density(approx?: boolean): Rational | bigint | RealLiteral {
    if (approx) {
      const packed = mzd_init(
        this.nrows,
        this.ncols,
        this._entries.map((row) =>
          this.ncols ? BigInt('0b' + row.slice().reverse().join('')) : 0n
        )
      );
      try {
        return create_RealNumber(mzd_density(packed, 1));
      } catch (error) {
        if ((error as Error).name === 'RuntimeError')
          throw new RuntimeError((error as Error).message);
        throw error;
      }
    }
    if (this.nrows === 0 || this.ncols === 0) {
      return 0n;
    }

    let count = 0n;
    for (let i = 0; i < this.nrows; i++) {
      for (let j = 0; j < this.ncols; j++) {
        if (this._entries[i]![j] === 1) {
          count++;
        }
      }
    }

    return new Rational(count, BigInt(this.nrows) * BigInt(this.ncols));
  }

  /**
   * Return the rank.
   *
   * @param algorithm - 'ple' or 'm4ri'
   * @returns The rank
   * @see Reference: sage/matrix/matrix_mod2_dense.pyx:rank
   */
  rank(algorithm = 'ple'): number {
    const cached = this._cache.get('rank');
    if (cached !== undefined) return cached as number;
    if (!this.nrows || !this.ncols) return 0;
    if (algorithm !== 'ple' && algorithm !== 'm4ri')
      throw new ValueError(`Algorithm '${algorithm}' unknown.`);
    const a = this.packed(),
      rank = (algorithm === 'ple' ? mzd_ple(a) : mzd_echelonize_m4ri(a, false)).rank;
    this._cache.set('rank', rank);
    return rank;
  }

  /**
   * Solve a linear system.
   *
   * @param B - Right-hand side
   * @param check - Whether to check the solution
   * @returns The solution
   * @see Reference: sage/matrix/matrix_mod2_dense.pyx:_solve_right_general
   * @see Deviation: Binary Linear Solve and Kernel Adapters
   */
  solve_right(B: Matrix_mod2_dense, check = true): Matrix_mod2_dense {
    if (this.nrows !== B.nrows)
      throw new ValueError('number of rows of self must equal number of rows of right-hand side');
    if (this.nrows === this.ncols && this.rank() === this.nrows) {
      const augmented = this.augment(B);
      augmented.echelonize();
      return augmented.submatrix(0, this.ncols, this.ncols, B.ncols);
    }
    const result = new Matrix_mod2_dense(this.ncols, B.ncols);
    if (!this.ncols || !B.ncols) {
      if (check && B._entries.some((row) => row.some(Boolean)))
        throw new ValueError('matrix equation has no solutions');
      return result;
    }
    const count = Math.max(this.nrows, this.ncols),
      a = this.packed(),
      b = B.packed();
    const solved = mzd_solve_left(
      mzd_init(count, this.ncols, [...a.rows, ...Array<bigint>(count - this.nrows).fill(0n)]),
      mzd_init(count, B.ncols, [...b.rows, ...Array<bigint>(count - B.nrows).fill(0n)]),
      0,
      check
    );
    if (solved.status !== 0) throw new ValueError('matrix equation has no solutions');
    result.loadPacked(mzd_init(this.ncols, B.ncols, solved.rhs.rows.slice(0, this.ncols)));
    return result;
  }

  /**
   * Return a matrix whose rows form a basis for the right kernel.
   *
   * @returns Kernel matrix
   * @see Reference: sage/matrix/matrix_mod2_dense.pyx:_right_kernel_matrix
   * @see Deviation: Binary Linear Solve and Kernel Adapters
   */
  right_kernel_matrix(options?: {
    basis?: 'default' | 'echelon' | 'pivot' | 'computed' | 'LLL';
    algorithm?: 'default' | 'generic' | 'pluq' | 'flint' | 'linbox' | 'pari' | 'padic';
    proof?: boolean | null;
  }): Matrix_mod2_dense {
    const algorithm = options?.algorithm ?? 'default';
    if (!['default', 'generic', 'flint', 'linbox', 'pari', 'padic', 'pluq'].includes(algorithm))
      throw new ValueError(`matrix kernel algorithm '${algorithm}' not recognized`);
    if (algorithm === 'padic' || algorithm === 'flint')
      throw new ValueError(
        `'${algorithm}' matrix kernel algorithm only available over the rationals and the integers, not over Finite Field of size 2`
      );
    if (algorithm === 'linbox')
      throw new ValueError(
        "'linbox' matrix kernel algorithm only available over the rationals, not over Finite Field of size 2"
      );
    if (algorithm === 'pari')
      throw new ValueError(
        "'pari' matrix kernel algorithm only available over non-trivial number fields and the integers, not over Finite Field of size 2"
      );
    const requested = options?.basis ?? 'default',
      basis = requested === 'default' ? 'echelon' : requested;
    if (!['computed', 'echelon', 'pivot', 'LLL'].includes(basis))
      throw new ValueError(`matrix kernel basis format '${basis}' not recognized`);
    if (basis === 'LLL')
      throw new ValueError(
        'LLL-reduced basis only available over the integers, not over Finite Field of size 2'
      );
    const proof: unknown = options?.proof;
    if (proof !== undefined && proof !== null) {
      if (![true, false, 0, 1, 0n, 1n].includes(proof as any))
        throw new ValueError(`'proof' must be one of True, False or None, not ${proof}`);
      throw new ValueError("'proof' flag only valid for matrices over the integers");
    }
    if (!this.ncols) return new Matrix_mod2_dense(0, 0);
    if (!this.nrows) {
      const identity = identity_matrix_gf2(this.ncols);
      identity.set_immutable();
      return identity;
    }
    if (algorithm === 'generic') return this.genericKernel(basis);
    const packed = mzd_kernel_left_pluq(this.packed()).kernel;
    let result: Matrix_mod2_dense;
    if (packed === null) result = new Matrix_mod2_dense(0, this.ncols);
    else {
      const columns = new Matrix_mod2_dense(packed.nrows, packed.ncols);
      columns.loadPacked(packed);
      result = columns.transpose();
    }
    if (basis === 'echelon') result.echelonize();
    else if (basis === 'pivot') {
      const pivots = new Set(this.pivots());
      const free = Array.from({ length: this.ncols }, (_, i) => i).filter((i) => !pivots.has(i));
      const C = new Matrix_mod2_dense(
        result.nrows,
        free.length,
        result._entries.map((row) => free.map((j) => row[j]!))
      );
      result = C.inverse().mul(result);
    }
    return result;
  }

  private genericKernel(basis: 'echelon' | 'computed' | 'pivot'): Matrix_mod2_dense {
    const m = this.nrows;
    const n = this.ncols;

    // Compute the reduced row echelon form
    const echelon = this.echelon_form();

    // Find pivot columns
    const pivotCols: number[] = [];
    let pivotRow = 0;
    for (let j = 0; j < n && pivotRow < m; j++) {
      if (echelon._entries[pivotRow]![j] === 1) {
        pivotCols.push(j);
        pivotRow++;
      }
    }

    // The free columns give the kernel basis
    const freeCols: number[] = [];
    for (let j = 0; j < n; j++) {
      if (!pivotCols.includes(j)) {
        freeCols.push(j);
      }
    }

    const kernelDim = freeCols.length;
    if (kernelDim === 0) {
      return new Matrix_mod2_dense(0, n);
    }

    // Build kernel basis
    const kernel = new Matrix_mod2_dense(kernelDim, n);

    for (let i = 0; i < kernelDim; i++) {
      const freeCol = freeCols[i]!;
      kernel._entries[i]![freeCol] = 1;

      // Fill in the pivot columns using the echelon form
      for (let r = 0; r < pivotCols.length; r++) {
        const pivotCol = pivotCols[r]!;
        // In GF(2), negation is identity
        kernel._entries[i]![pivotCol] = echelon._entries[r]![freeCol]!;
      }
    }

    // `matrix2.pyx:4266 right_kernel_matrix` defaults to `basis='echelon'`
    // over a field, and GF(2) is one.  Above we built the 'computed'/'pivot'
    // basis (they coincide in characteristic 2, where negation is the
    // identity); returning that unechelonized was the defect.
    if (basis === 'echelon') {
      kernel.echelonize();
    }

    return kernel;
  }

  /**
   * Return a doubly lexical ordering of the matrix.
   *
   * A doubly lexical ordering of a matrix is an ordering of the rows
   * and of the columns of the matrix so that both the rows and the
   * columns, as vectors, are lexically increasing. See [Lub1987].
   * A lexical ordering of vectors is the standard dictionary ordering,
   * except that vectors will be read from highest to lowest coordinate.
   * Thus row vectors will be compared from right to left, and column
   * vectors from bottom to top.
   *
   * @param inplace - Whether to modify in place (default: false)
   * @returns A pair [row_ordering, col_ordering] representing permutations (1-indexed)
   * @see Reference: sage/matrix/matrix_mod2_dense.pyx:doubly_lexical_ordering
   */
  doubly_lexical_ordering(inplace?: boolean): [number[], number[]] {
    if (inplace && this._immutable)
      throw new TypeError(
        'this matrix is immutable; use inplace=False or apply to a mutable copy.'
      );
    if (this.nrows === 0 || this.ncols === 0) {
      // Identity permutations for empty dimensions
      const rowPerm = Array.from({ length: this.nrows }, (_, i) => i + 1);
      const colPerm = Array.from({ length: this.ncols }, (_, i) => i + 1);
      return [rowPerm, colPerm];
    }

    // Track partition boundaries between rows
    const partitionRows: boolean[] = new Array(this.nrows - 1).fill(false);
    let partitionNum = 1;

    // Track row and column swaps (1-indexed as in SageMath)
    const rowSwapped: number[] = Array.from({ length: this.nrows }, (_, i) => i + 1);
    const colSwapped: number[] = Array.from({ length: this.ncols }, (_, i) => i + 1);

    // Work on a copy or self
    const A = inplace ? this : this.copy();

    // Process columns from right to left
    for (let i = this.ncols; i >= 1; i--) {
      // Count 1s for each partition and column
      // count1[col][partitionIdx] = count of 1s in that column within that partition
      const count1: number[][] = [];
      for (let col = 0; col < i; col++) {
        count1.push(new Array(partitionNum).fill(0));
      }

      for (let col = 0; col < i; col++) {
        let partitionIdx = 0;
        // Iterate rows from bottom to top
        for (let row = this.nrows - 1; row >= 1; row--) {
          count1[col]![partitionIdx] += A._entries[row]![col]!;
          if (partitionRows[row - 1]) {
            partitionIdx++;
          }
        }
        // Special case for row 0
        count1[col]![partitionIdx] += A._entries[0]![col]!;
      }

      // Find the column with lexicographically largest count vector.
      // Sage computes ``max((c, i) for i, c in enumerate(count1))``: ties in the
      // count vector are broken towards the *largest* column index.
      let largestCol = 0;
      let largestCount = count1[0]!;
      for (let col = 1; col < i; col++) {
        const currentCount = count1[col]!;
        // Compare lexicographically; `cmp > 0` means strictly larger,
        // `cmp === 0` means equal (and then the larger index wins).
        let cmp = 0;
        for (let p = 0; p < partitionNum; p++) {
          if (currentCount[p]! > largestCount[p]!) {
            cmp = 1;
            break;
          } else if (currentCount[p]! < largestCount[p]!) {
            cmp = -1;
            break;
          }
        }
        if (cmp >= 0) {
          largestCol = col;
          largestCount = currentCount;
        }
      }

      // Refine partitions and reorder rows according to values in largestCol
      let partitionStart = 0;
      for (let pIdx = 0; pIdx < partitionNum; pIdx++) {
        // Find partition end
        let partitionEnd = partitionStart;
        while (partitionEnd < this.nrows - 1 && !partitionRows[partitionEnd]) {
          partitionEnd++;
        }

        // Sort: move rows with 0 to top, rows with 1 to bottom within partition
        let rowStart = partitionStart;
        let rowEnd = partitionEnd;
        while (rowStart < rowEnd) {
          // Find first row with 1 from top
          while (rowStart < rowEnd && A._entries[rowStart]![largestCol] === 0) {
            rowStart++;
          }
          // Find first row with 0 from bottom
          while (rowStart < rowEnd && A._entries[rowEnd]![largestCol] === 1) {
            rowEnd--;
          }
          if (rowStart < rowEnd) {
            // Swap rows
            [A._entries[rowStart], A._entries[rowEnd]] = [
              A._entries[rowEnd]!,
              A._entries[rowStart]!,
            ];
            [rowSwapped[rowStart], rowSwapped[rowEnd]] = [
              rowSwapped[rowEnd]!,
              rowSwapped[rowStart]!,
            ];
          }
        }

        partitionStart = partitionEnd + 1;
      }

      // Create new partition boundaries based on transitions in largestCol
      for (let row = 0; row < this.nrows - 1; row++) {
        if (A._entries[row]![largestCol] !== A._entries[row + 1]![largestCol]) {
          if (!partitionRows[row]) {
            partitionRows[row] = true;
            partitionNum++;
          }
        }
      }

      // Swap the largest column to position i-1
      for (const row of A._entries) [row[largestCol], row[i - 1]] = [row[i - 1]!, row[largestCol]!];
      [colSwapped[largestCol], colSwapped[i - 1]] = [colSwapped[i - 1]!, colSwapped[largestCol]!];
    }

    return [rowSwapped, colSwapped];
  }

  /**
   * Check if the matrix is Gamma-free.
   *
   * A matrix is Gamma-free if it does not contain a 2x2 submatrix of the form:
   *   [1 1]
   *   [1 0]
   *
   * @param certificate - Whether to return a certificate
   * @returns True if Gamma-free, or [true/false, [r1, c1, r2, c2] | null] if certificate=true
   * @see Reference: sage/matrix/matrix_mod2_dense.pyx:is_Gamma_free
   */
  is_Gamma_free(
    certificate?: boolean
  ): boolean | [boolean, [number, number, number, number] | null] {
    // For each 1 entry, find the next 1 in the same row and the next 1 in the same column,
    // and check if the 2x2 submatrix forms a Gamma pattern
    for (let i = 0; i < this.nrows; i++) {
      let j = 0;
      while (j < this.ncols) {
        if (this._entries[i]![j] === 1) {
          // Find the next 1 in the row
          let jRight = j + 1;
          while (jRight < this.ncols && this._entries[i]![jRight] !== 1) {
            jRight++;
          }

          if (jRight < this.ncols) {
            // Found A[i][j] = 1 and A[i][jRight] = 1
            // Now find the next 1 in column j
            let iBottom = i + 1;
            while (iBottom < this.nrows && this._entries[iBottom]![j] !== 1) {
              iBottom++;
            }

            if (iBottom < this.nrows && this._entries[iBottom]![jRight] === 0) {
              // Found Gamma pattern:
              // A[i][j] = 1, A[i][jRight] = 1
              // A[iBottom][j] = 1, A[iBottom][jRight] = 0
              if (certificate) {
                return [false, [i, j, iBottom, jRight]];
              }
              return false;
            }
          }
          j = jRight;
        } else {
          j++;
        }
      }
    }

    if (certificate) {
      return [true, null];
    }
    return true;
  }

  /**
   * Swap two rows.
   *
   * @param i - First row index
   * @param j - Second row index
   * @see Reference: sage/matrix/matrix_mod2_dense.pyx:swap_rows_c
   * @see Deviation: Binary Matrix Indices and Slices
   */
  swap_rows(first: BinaryMatrixIndex, second: BinaryMatrixIndex): void {
    const a = binaryMatrixIndex(first, 'index'),
      b = binaryMatrixIndex(second, 'index');
    this.check_mutability(true);
    if (a < 0n || a >= BigInt(this.nrows) || b < 0n || b >= BigInt(this.nrows)) {
      throw new IndexError('matrix row index out of range');
    }
    const i = Number(a),
      j = Number(b);
    if (i !== j) {
      [this._entries[i], this._entries[j]] = [this._entries[j]!, this._entries[i]!];
    }
  }

  /**
   * Swap two columns.
   *
   * @param i - First column index
   * @param j - Second column index
   * @see Reference: sage/matrix/matrix_mod2_dense.pyx:swap_columns_c
   * @see Deviation: Binary Matrix Indices and Slices
   */
  swap_columns(first: BinaryMatrixIndex, second: BinaryMatrixIndex): void {
    const a = binaryMatrixIndex(first, 'index'),
      b = binaryMatrixIndex(second, 'index');
    this.check_mutability(true);
    if (a < 0n || a >= BigInt(this.ncols) || b < 0n || b >= BigInt(this.ncols)) {
      throw new IndexError('matrix column index out of range');
    }
    const i = Number(a),
      j = Number(b);
    if (i !== j) {
      for (let row = 0; row < this.nrows; row++) {
        [this._entries[row]![i], this._entries[row]![j]] = [
          this._entries[row]![j]!,
          this._entries[row]![i]!,
        ];
      }
    }
  }

  /**
   * Permute rows according to a permutation.
   *
   * @param perm - Permutation (1-indexed, as in SageMath)
   * @see Reference: sage/matrix/matrix0.pyx:permute_rows
   * @see Deviation: Binary Matrix Permutation Image Lists
   */
  permute_rows(perm: number[]): void {
    this.permute(perm, false);
  }

  /** Apply a 1-based permutation image list to columns.
   * @see Reference: sage/matrix/matrix0.pyx:permute_columns
   * @see Deviation: Binary Matrix Permutation Image Lists
   */
  permute_columns(perm: number[]): void {
    this.permute(perm, true);
  }

  private permute(perm: number[], columns: boolean): void {
    if (
      new Set(perm).size !== perm.length ||
      perm.some((x) => !Number.isInteger(x) || x < 1 || x > perm.length)
    )
      throw new ValueError('invalid data to initialize a permutation');
    this.check_mutability();
    const seen = new Set<number>();
    for (let start = 0; start < perm.length; start++) {
      if (seen.has(start)) continue;
      const cycle: number[] = [];
      let next = start;
      do {
        cycle.push(next);
        seen.add(next);
        next = perm[next]! - 1;
      } while (next !== start);
      if (cycle.length === 1) continue;
      cycle.reverse();
      const first = cycle[0]!;
      for (const other of cycle) {
        // Preserve native cycle order, including changes before a later bounds error.
        if (columns) this.swap_columns(first, other);
        else this.swap_rows(first, other);
      }
    }
  }

  /**
   * Check equality with another matrix.
   *
   * @param other - Matrix to compare
   * @returns True if equal
   */
  eq(other: Matrix_mod2_dense): boolean {
    if (this.nrows !== other.nrows || this.ncols !== other.ncols) {
      return false;
    }

    for (let i = 0; i < this.nrows; i++) {
      for (let j = 0; j < this.ncols; j++) {
        if (this._entries[i]![j] !== other._entries[i]![j]) {
          return false;
        }
      }
    }
    return true;
  }

  /**
   * Return string representation.
   *
   * @returns String representation
   */
  toString(): string {
    return this.str();
  }
}

// ============================================================================
// Module-level functions
// ============================================================================

/**
 * Create a GF(2) matrix from a PNG image data.
 *
 * NOTE: This is a TypeScript adaptation. The original SageMath function uses
 * libgd to read PNG files directly from disk. In TypeScript/JavaScript, we
 * instead work with image data that has already been loaded.
 *
 * The function expects raw pixel data where each pixel is either 0 (black=1)
 * or 255 (white=0). This matches SageMath's convention where black pixels
 * represent 1s and white pixels represent 0s.
 *
 * @param width - Image width in pixels
 * @param height - Image height in pixels
 * @param pixels - Grayscale pixel values (0-255), row-major order
 * @returns Matrix from image data
 * @see Reference: sage/matrix/matrix_mod2_dense.pyx:from_png
 */
export function from_png_data(
  width: number,
  height: number,
  pixels: Uint8Array | number[]
): Matrix_mod2_dense {
  if (pixels.length !== width * height) {
    throw new ValueError(`Expected ${width * height} pixels, got ${pixels.length}`);
  }

  const A = new Matrix_mod2_dense(height, width);

  for (let i = 0; i < height; i++) {
    for (let j = 0; j < width; j++) {
      // Black (0) -> 1, White (255) -> 0
      // Use threshold of 128 for grayscale images
      const pixelValue = pixels[i * width + j]!;
      A.set(i, j, pixelValue < 128 ? 1 : 0);
    }
  }

  return A;
}

/**
 * Create a GF(2) matrix from a PNG file.
 *
 * NOTE: File I/O is not available in all JavaScript environments.
 * This function is provided for API compatibility but will throw in
 * environments without file system access. Use from_png_data() for
 * a more portable alternative.
 *
 * @param filename - Path to PNG file
 * @returns Matrix from PNG
 * @throws {NotImplementedError} In environments without file system access
 * @see Reference: sage/matrix/matrix_mod2_dense.pyx:from_png
 */
export function from_png(filename: string): Matrix_mod2_dense {
  // This would require a PNG decoder library and file system access.
  // In Node.js, one could use 'sharp' or 'pngjs' packages.
  // For browser environments, one would use Canvas API.
  throw new NotImplementedError(
    'from_png requires environment-specific PNG decoding. Use from_png_data() with pre-loaded image data instead.'
  );
}

/**
 * Convert a GF(2) matrix to grayscale pixel data suitable for PNG encoding.
 *
 * NOTE: This is a TypeScript adaptation. The original SageMath function uses
 * libgd to write PNG files directly to disk. In TypeScript/JavaScript, we
 * instead return raw pixel data that can be encoded using environment-specific
 * PNG libraries.
 *
 * Returns grayscale pixel values where 1s become black (0) and 0s become white (255).
 *
 * @param A - The matrix
 * @returns Object with width, height, and pixels (Uint8Array in row-major order)
 * @see Reference: sage/matrix/matrix_mod2_dense.pyx:to_png
 * @see Deviation: GF(2) Matrix PNG Functions
 */
export function to_png_data(A: Matrix_mod2_dense): {
  width: number;
  height: number;
  pixels: Uint8Array;
} {
  const r = A.nrows;
  const c = A.ncols;

  if (r === 0 || c === 0) {
    throw new TypeError(`cannot write image with dimensions ${c} x ${r}`);
  }

  const pixels = new Uint8Array(r * c);

  for (let i = 0; i < r; i++) {
    for (let j = 0; j < c; j++) {
      // 1 -> black (0), 0 -> white (255)
      pixels[i * c + j] = A.get(i, j) === 1 ? 0 : 255;
    }
  }

  return { width: c, height: r, pixels };
}

/**
 * Write a GF(2) matrix to a PNG file.
 *
 * NOTE: File I/O is not available in all JavaScript environments.
 * This function is provided for API compatibility but will throw in
 * environments without file system access. Use to_png_data() for
 * a more portable alternative.
 *
 * @param A - The matrix
 * @param filename - Path to output file
 * @throws {NotImplementedError} In environments without file system access
 * @see Reference: sage/matrix/matrix_mod2_dense.pyx:to_png
 * @see Deviation: GF(2) Matrix PNG Functions
 */
export function to_png(A: Matrix_mod2_dense, filename: string): void {
  // This would require a PNG encoder library and file system access.
  // In Node.js, one could use 'sharp' or 'pngjs' packages.
  // For browser environments, one would use Canvas API.
  throw new NotImplementedError(
    'to_png requires environment-specific PNG encoding. Use to_png_data() to get pixel data for encoding.'
  );
}

/**
 * Compute PLUQ factorization of A.
 *
 * Returns a tuple (LU, P, Q) where:
 * - LU is a matrix containing both L (lower, strictly below the diagonal) and
 *   U (upper, on and above the diagonal) in packed form, exactly as M4RI
 *   stores them in place over A
 * - P is the row permutation
 * - Q is the column permutation
 *
 * As in M4RI (`mzp_t`), `P` and `Q` are **transposition lists**, not the
 * composed permutation: they are applied as `for i in 0..len-1: swap(i, X[i])`.
 *
 * The factorization satisfies: A = P * L * U * Q.
 *
 * @param A - The matrix
 * @param algorithm - 'standard', 'mmpf', or 'naive' (default: 'standard')
 * @param param - Algorithm parameter (default: 0)
 * @returns Tuple [LU, P, Q]
 * @see Deviation: Standalone Binary Factorizations
 * @see Reference: sage/matrix/matrix_mod2_dense.pyx:2707 (pluq), M4RI mzd_pluq
 */
export function pluq(
  A: Matrix_mod2_dense,
  algorithm = 'standard',
  param: BinaryFactorizationParameter = 0
): [Matrix_mod2_dense, number[], number[]] {
  return binaryFactorization(A, algorithm, param, true);
}

/**
 * Compute PLE factorization of A: `A = P * L * E` with `P` a permutation,
 * `L` m x r unit lower triangular and `E` an r x n matrix in row echelon form.
 *
 * Returns a tuple (LU, P, Q) where:
 * - LU packs L and E in place over A exactly the way M4RI's `mzd_ple` does:
 *   `LU[i][j] = L[i][j]` for `j < i` (the multipliers, **compacted to the
 *   left**), `LU[i][i] = 1` for `i < r` (the unit diagonal / pivot), and
 *   `LU[i][j] = E[i][j]` for `j > Q[i]` (the free part of the echelon row, at
 *   its original column). Everything else is zero. Applying the column
 *   permutation Q to this layout yields exactly the PLUQ layout.
 * - P is the row permutation, as a transposition list
 * - Q[i] is the pivot column of row i for `i < r`, and `i` otherwise (which is
 *   the same thing as M4RI's transposition list for the column permutation)
 *
 * @param A - The matrix
 * @param algorithm - 'standard', 'russian', or 'naive' (default: 'standard')
 * @param param - Algorithm parameter (default: 0)
 * @returns Tuple [LU, P, Q]
 * @see Deviation: Standalone Binary Factorizations
 * @see Reference: sage/matrix/matrix_mod2_dense.pyx:2771 (ple), M4RI mzd_ple
 */
export function ple(
  A: Matrix_mod2_dense,
  algorithm = 'standard',
  param: BinaryFactorizationParameter = 0
): [Matrix_mod2_dense, number[], number[]] {
  return binaryFactorization(A, algorithm, param, false);
}
type BinaryFactorizationParameter = number | IntegerLike | Rational | boolean;
/** Sage's C-int wrapper conversion occurs before algorithm dispatch. */
function binaryFactorization(
  A: Matrix_mod2_dense,
  algorithm: string,
  param: BinaryFactorizationParameter,
  pluq: boolean
): [Matrix_mod2_dense, number[], number[]] {
  let value: bigint;
  if (typeof param === 'number') {
    if (Number.isNaN(param)) throw new ValueError('cannot convert float NaN to integer');
    if (!Number.isFinite(param))
      throw new OverflowError('cannot convert float infinity to integer');
    value = BigInt(Math.trunc(param));
  } else if (typeof param === 'bigint' || typeof param === 'boolean') value = BigInt(param);
  else if (param instanceof Integer) value = param.value;
  else if (param instanceof Rational) value = param.trunc();
  else throw new TypeError('an integer is required');
  if (value < -(1n << 63n) || value >= 1n << 63n)
    throw new OverflowError('Python int too large to convert to C long');
  if (value < -(1n << 31n) || value >= 1n << 31n)
    throw new OverflowError('value too large to convert to int');
  const k = Number(value),
    russian = pluq ? 'mmpf' : 'russian';
  if (!['standard', russian, 'naive'].includes(algorithm))
    throw new ValueError(`Algorithm '${algorithm}' unknown.`);
  const rows = Array.from({ length: A.nrows }, (_, i) =>
      A.row(i).reduce((v, b, j) => v | (BigInt(b) << BigInt(j)), 0n)
    ),
    a = mzd_init(A.nrows, A.ncols, rows);
  const result =
    algorithm === 'standard'
      ? (pluq ? mzd_pluq : mzd_ple)(a, k)
      : algorithm === 'naive'
        ? (pluq ? _mzd_pluq_naive : _mzd_ple_naive)(a)
        : (pluq ? _mzd_pluq_russian : _mzd_ple_russian)(a, k);
  const B = new Matrix_mod2_dense(
    A.nrows,
    A.ncols,
    result.matrix.rows.map((row) =>
      Array.from({ length: A.ncols }, (_, j) => Number((row >> BigInt(j)) & 1n))
    )
  );
  return [B, result.P, result.Q];
}

// ============================================================================
// Factory functions
// ============================================================================

/**
 * Create a zero matrix over GF(2).
 *
 * @param nrows - Number of rows
 * @param ncols - Number of columns (default: nrows)
 * @returns Zero matrix
 * @see Deviation: Binary Seeded Randomization and Factories
 */
export function zero_matrix_gf2(
  nrows: BinaryMatrixIndex,
  ncols?: BinaryMatrixIndex
): Matrix_mod2_dense {
  return new Matrix_mod2_dense(nrows, ncols ?? nrows);
}

/**
 * Create an identity matrix over GF(2).
 *
 * @param n - Size
 * @returns Identity matrix
 * @see Deviation: Binary Seeded Randomization and Factories
 */
export function identity_matrix_gf2(n: BinaryMatrixIndex): Matrix_mod2_dense {
  return new Matrix_mod2_dense(n, n, 1);
}

/**
 * Create a random matrix over GF(2).
 *
 * @param nrows - Number of rows
 * @param ncols - Number of columns (default: nrows)
 * @param density - Explicit density selects nonzero entries; omitted density is uniform
 * @returns Random matrix
 * @see Deviation: Binary Seeded Randomization and Factories
 */
export function random_matrix_gf2(
  nrows: BinaryMatrixIndex,
  ncols?: BinaryMatrixIndex,
  density?: FloatInput
): Matrix_mod2_dense {
  const result = new Matrix_mod2_dense(nrows, ncols ?? nrows);
  result.randomize(density ?? 1.0, density != null);
  return result;
}

/**
 * Create a matrix over GF(2) from entries.
 *
 * @param entries - Nested coefficient arrays, converted through GF(2)
 * @see Deviation: Binary Seeded Randomization and Factories
 * @returns The matrix
 */
export function matrix_gf2_from_entries(entries: unknown[][]): Matrix_mod2_dense {
  if (entries.length === 0) {
    return new Matrix_mod2_dense(0, 0);
  }

  const nrows = entries.length;
  const ncols = entries[0]!.length;

  // Validate row lengths
  for (let i = 1; i < nrows; i++) {
    if (entries[i]!.length !== ncols) {
      throw new ValueError(
        `inconsistent number of columns: should be ${ncols} but got ${entries[i]!.length}`
      );
    }
  }

  return new Matrix_mod2_dense(nrows, ncols, entries);
}
