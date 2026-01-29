/**
 * @module sage/matrix/matrix_generic
 * @description Generic matrix class over arbitrary rings
 *
 * Port of: sage/matrix/matrix_generic_dense.py
 */

import { ArithmeticError, ValueError } from '../errors.js';
import type { CoefficientRing, RingElement } from '../rings/polynomial/polynomial_element.js';

/**
 * A matrix over a ring R.
 *
 * Matrices are stored as a 2D array of entries. All matrix operations
 * are performed using the ring operations of the base ring.
 *
 * @example
 * ```typescript
 * const F7 = GF(7n);
 * const MS = MatrixSpace(F7, 2, 2);
 * const A = MS([[1, 2], [3, 4]]);
 * const B = MS([[5, 6], [7, 8]]);
 * console.log(A.add(B));  // Matrix over GF(7)
 * console.log(A.mul(B));  // Matrix multiplication
 * ```
 */
export class Matrix<R extends RingElement> {
  private _entries: R[][];
  readonly nrows: number;
  readonly ncols: number;
  readonly base_ring: CoefficientRing<R>;

  /**
   * Create a matrix with the given entries.
   *
   * @param base_ring - The ring of coefficients
   * @param nrows - Number of rows
   * @param ncols - Number of columns
   * @param entries - 2D array of entries, or flat array, or a function (i, j) => R
   */
  constructor(
    base_ring: CoefficientRing<R>,
    nrows: number,
    ncols: number,
    entries?: R[][] | R[] | ((i: number, j: number) => R)
  ) {
    if (nrows < 0 || ncols < 0) {
      throw new ValueError('matrix dimensions must be non-negative');
    }

    this.base_ring = base_ring;
    this.nrows = nrows;
    this.ncols = ncols;

    // Initialize the entries array
    this._entries = [];
    for (let i = 0; i < nrows; i++) {
      this._entries.push([]);
      for (let j = 0; j < ncols; j++) {
        this._entries[i]!.push(base_ring.zero());
      }
    }

    // Fill in entries if provided
    if (entries !== undefined) {
      if (typeof entries === 'function') {
        // entries is a function (i, j) => R
        for (let i = 0; i < nrows; i++) {
          for (let j = 0; j < ncols; j++) {
            this._entries[i]![j] = entries(i, j);
          }
        }
      } else if (entries.length > 0 && Array.isArray(entries[0])) {
        // entries is a 2D array
        const entries2D = entries as R[][];
        if (entries2D.length !== nrows) {
          throw new ValueError(`expected ${nrows} rows, got ${entries2D.length}`);
        }
        for (let i = 0; i < nrows; i++) {
          if (entries2D[i]!.length !== ncols) {
            throw new ValueError(
              `expected ${ncols} columns in row ${i}, got ${entries2D[i]!.length}`
            );
          }
          for (let j = 0; j < ncols; j++) {
            this._entries[i]![j] = entries2D[i]![j]!;
          }
        }
      } else {
        // entries is a flat array
        const flat = entries as R[];
        if (flat.length !== nrows * ncols) {
          throw new ValueError(`expected ${nrows * ncols} entries, got ${flat.length}`);
        }
        for (let i = 0; i < nrows; i++) {
          for (let j = 0; j < ncols; j++) {
            this._entries[i]![j] = flat[i * ncols + j]!;
          }
        }
      }
    }
  }

  /**
   * Get the entry at position (i, j).
   *
   * @param i - Row index (0-based)
   * @param j - Column index (0-based)
   */
  get(i: number, j: number): R {
    if (i < 0 || i >= this.nrows || j < 0 || j >= this.ncols) {
      throw new ValueError(
        `index (${i}, ${j}) out of bounds for ${this.nrows}x${this.ncols} matrix`
      );
    }
    return this._entries[i]![j]!;
  }

  /**
   * Set the entry at position (i, j).
   *
   * @param i - Row index (0-based)
   * @param j - Column index (0-based)
   * @param value - The new value
   */
  set(i: number, j: number, value: R): void {
    if (i < 0 || i >= this.nrows || j < 0 || j >= this.ncols) {
      throw new ValueError(
        `index (${i}, ${j}) out of bounds for ${this.nrows}x${this.ncols} matrix`
      );
    }
    this._entries[i]![j] = value;
  }

  /**
   * Get a copy of row i as an array.
   *
   * @param i - Row index (0-based)
   */
  row(i: number): R[] {
    if (i < 0 || i >= this.nrows) {
      throw new ValueError(`row index ${i} out of bounds for ${this.nrows}x${this.ncols} matrix`);
    }
    return [...this._entries[i]!];
  }

  /**
   * Get a copy of column j as an array.
   *
   * @param j - Column index (0-based)
   */
  column(j: number): R[] {
    if (j < 0 || j >= this.ncols) {
      throw new ValueError(
        `column index ${j} out of bounds for ${this.nrows}x${this.ncols} matrix`
      );
    }
    const col: R[] = [];
    for (let i = 0; i < this.nrows; i++) {
      col.push(this._entries[i]![j]!);
    }
    return col;
  }

  /**
   * Return a list of all rows.
   */
  rows(): R[][] {
    return this._entries.map((row) => [...row]);
  }

  /**
   * Return a list of all columns.
   */
  columns(): R[][] {
    const cols: R[][] = [];
    for (let j = 0; j < this.ncols; j++) {
      cols.push(this.column(j));
    }
    return cols;
  }

  /**
   * Return the entries as a flat list (row-major order).
   */
  list(): R[] {
    const result: R[] = [];
    for (let i = 0; i < this.nrows; i++) {
      for (let j = 0; j < this.ncols; j++) {
        result.push(this._entries[i]![j]!);
      }
    }
    return result;
  }

  /**
   * Check if this is a square matrix.
   */
  is_square(): boolean {
    return this.nrows === this.ncols;
  }

  /**
   * Check if this is the zero matrix.
   */
  is_zero(): boolean {
    for (let i = 0; i < this.nrows; i++) {
      for (let j = 0; j < this.ncols; j++) {
        if (!this._entries[i]![j]!.isZero()) {
          return false;
        }
      }
    }
    return true;
  }

  /**
   * Add two matrices.
   */
  add(other: Matrix<R>): Matrix<R> {
    if (this.nrows !== other.nrows || this.ncols !== other.ncols) {
      throw new ArithmeticError(
        `cannot add ${this.nrows}x${this.ncols} matrix to ${other.nrows}x${other.ncols} matrix`
      );
    }

    const result = new Matrix<R>(this.base_ring, this.nrows, this.ncols);
    for (let i = 0; i < this.nrows; i++) {
      for (let j = 0; j < this.ncols; j++) {
        result._entries[i]![j] = this._entries[i]![j]!.add(other._entries[i]![j]!) as R;
      }
    }
    return result;
  }

  /**
   * Subtract two matrices.
   */
  sub(other: Matrix<R>): Matrix<R> {
    if (this.nrows !== other.nrows || this.ncols !== other.ncols) {
      throw new ArithmeticError(
        `cannot subtract ${other.nrows}x${other.ncols} matrix from ${this.nrows}x${this.ncols} matrix`
      );
    }

    const result = new Matrix<R>(this.base_ring, this.nrows, this.ncols);
    for (let i = 0; i < this.nrows; i++) {
      for (let j = 0; j < this.ncols; j++) {
        result._entries[i]![j] = this._entries[i]![j]!.sub(other._entries[i]![j]!) as R;
      }
    }
    return result;
  }

  /**
   * Return the negation of this matrix.
   */
  neg(): Matrix<R> {
    const result = new Matrix<R>(this.base_ring, this.nrows, this.ncols);
    for (let i = 0; i < this.nrows; i++) {
      for (let j = 0; j < this.ncols; j++) {
        result._entries[i]![j] = this._entries[i]![j]!.neg() as R;
      }
    }
    return result;
  }

  /**
   * Multiply this matrix by another matrix.
   */
  mul(other: Matrix<R>): Matrix<R> {
    if (this.ncols !== other.nrows) {
      throw new ArithmeticError(
        `cannot multiply ${this.nrows}x${this.ncols} matrix by ${other.nrows}x${other.ncols} matrix`
      );
    }

    const result = new Matrix<R>(this.base_ring, this.nrows, other.ncols);
    for (let i = 0; i < this.nrows; i++) {
      for (let j = 0; j < other.ncols; j++) {
        let sum = this.base_ring.zero();
        for (let k = 0; k < this.ncols; k++) {
          sum = sum.add(this._entries[i]![k]!.mul(other._entries[k]![j]!) as R) as R;
        }
        result._entries[i]![j] = sum;
      }
    }
    return result;
  }

  /**
   * Multiply this matrix by a scalar.
   */
  scalar_mul(c: R): Matrix<R> {
    const result = new Matrix<R>(this.base_ring, this.nrows, this.ncols);
    for (let i = 0; i < this.nrows; i++) {
      for (let j = 0; j < this.ncols; j++) {
        result._entries[i]![j] = this._entries[i]![j]!.mul(c) as R;
      }
    }
    return result;
  }

  /**
   * Return the transpose of this matrix.
   */
  transpose(): Matrix<R> {
    const result = new Matrix<R>(this.base_ring, this.ncols, this.nrows);
    for (let i = 0; i < this.nrows; i++) {
      for (let j = 0; j < this.ncols; j++) {
        result._entries[j]![i] = this._entries[i]![j]!;
      }
    }
    return result;
  }

  /**
   * Return the trace of this matrix (sum of diagonal elements).
   *
   * @throws {ArithmeticError} If the matrix is not square
   */
  trace(): R {
    if (!this.is_square()) {
      throw new ArithmeticError('trace is only defined for square matrices');
    }

    let sum = this.base_ring.zero();
    for (let i = 0; i < this.nrows; i++) {
      sum = sum.add(this._entries[i]![i]!) as R;
    }
    return sum;
  }

  /**
   * Return a deep copy of this matrix.
   */
  copy(): Matrix<R> {
    const result = new Matrix<R>(this.base_ring, this.nrows, this.ncols);
    for (let i = 0; i < this.nrows; i++) {
      for (let j = 0; j < this.ncols; j++) {
        result._entries[i]![j] = this._entries[i]![j]!;
      }
    }
    return result;
  }

  /**
   * Check equality with another matrix.
   */
  eq(other: Matrix<R>): boolean {
    if (this.nrows !== other.nrows || this.ncols !== other.ncols) {
      return false;
    }

    for (let i = 0; i < this.nrows; i++) {
      for (let j = 0; j < this.ncols; j++) {
        if (!this._entries[i]![j]!.eq(other._entries[i]![j]!)) {
          return false;
        }
      }
    }
    return true;
  }

  /**
   * Return the diagonal elements as an array.
   */
  diagonal(): R[] {
    const n = Math.min(this.nrows, this.ncols);
    const result: R[] = [];
    for (let i = 0; i < n; i++) {
      result.push(this._entries[i]![i]!);
    }
    return result;
  }

  /**
   * Compute the power of this matrix.
   *
   * @param n - The exponent (non-negative integer)
   * @throws {ArithmeticError} If the matrix is not square
   */
  pow(n: number | bigint): Matrix<R> {
    if (!this.is_square()) {
      throw new ArithmeticError('matrix power is only defined for square matrices');
    }

    let exp = typeof n === 'bigint' ? n : BigInt(n);
    if (exp < 0n) {
      throw new ValueError('negative exponent not supported for generic matrices');
    }

    if (exp === 0n) {
      return identity_matrix(this.base_ring, this.nrows);
    }

    // Binary exponentiation
    let result = identity_matrix(this.base_ring, this.nrows);
    let base: Matrix<R> = this.copy();

    while (exp > 0n) {
      if ((exp & 1n) === 1n) {
        result = result.mul(base);
      }
      base = base.mul(base);
      exp >>= 1n;
    }

    return result;
  }

  /**
   * Return a string representation of this matrix.
   */
  toString(): string {
    if (this.nrows === 0 || this.ncols === 0) {
      return `${this.nrows} x ${this.ncols} empty matrix`;
    }

    // Convert all entries to strings and find max width per column
    const strings: string[][] = [];
    const widths: number[] = new Array(this.ncols).fill(0);

    for (let i = 0; i < this.nrows; i++) {
      strings.push([]);
      for (let j = 0; j < this.ncols; j++) {
        const s = this._entries[i]![j]!.toString();
        strings[i]!.push(s);
        widths[j] = Math.max(widths[j]!, s.length);
      }
    }

    // Build the output string
    const lines: string[] = [];
    for (let i = 0; i < this.nrows; i++) {
      const row = strings[i]!.map((s, j) => s.padStart(widths[j]!));
      lines.push('[' + row.join(' ') + ']');
    }

    return lines.join('\n');
  }
}

/**
 * Create a zero matrix over the given ring.
 *
 * @param ring - The base ring
 * @param nrows - Number of rows
 * @param ncols - Number of columns (defaults to nrows for square matrix)
 */
export function zero_matrix<R extends RingElement>(
  ring: CoefficientRing<R>,
  nrows: number,
  ncols?: number
): Matrix<R> {
  return new Matrix<R>(ring, nrows, ncols ?? nrows);
}

/**
 * Create an identity matrix over the given ring.
 *
 * @param ring - The base ring
 * @param n - The size of the matrix
 */
export function identity_matrix<R extends RingElement>(
  ring: CoefficientRing<R>,
  n: number
): Matrix<R> {
  const result = new Matrix<R>(ring, n, n);
  for (let i = 0; i < n; i++) {
    result.set(i, i, ring.one());
  }
  return result;
}

/**
 * Create a diagonal matrix from a list of diagonal elements.
 *
 * @param ring - The base ring
 * @param entries - The diagonal entries
 */
export function diagonal_matrix<R extends RingElement>(
  ring: CoefficientRing<R>,
  entries: R[]
): Matrix<R> {
  const n = entries.length;
  const result = new Matrix<R>(ring, n, n);
  for (let i = 0; i < n; i++) {
    result.set(i, i, entries[i]!);
  }
  return result;
}

/**
 * Create a scalar matrix (diagonal matrix with the same value on the diagonal).
 *
 * @param ring - The base ring
 * @param n - The size of the matrix
 * @param value - The scalar value
 */
export function scalar_matrix<R extends RingElement>(
  ring: CoefficientRing<R>,
  n: number,
  value: R
): Matrix<R> {
  const result = new Matrix<R>(ring, n, n);
  for (let i = 0; i < n; i++) {
    result.set(i, i, value);
  }
  return result;
}
