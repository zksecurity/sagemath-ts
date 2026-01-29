/**
 * @module sage/matrix/matrix_space
 * @description Matrix space factory for creating matrices over rings
 *
 * Port of: sage/matrix/matrix_space.py
 */

import { ValueError } from '../errors.js';
import type { CoefficientRing, RingElement } from '../rings/polynomial/polynomial_element.js';
import { Matrix, identity_matrix, zero_matrix } from './matrix_generic.js';

/**
 * A matrix space is the set of all nrows x ncols matrices over a ring R.
 *
 * This class serves as a factory for creating matrices with a fixed base ring
 * and dimensions.
 *
 * @example
 * ```typescript
 * const F7 = GF(7n);
 * const MS = MatrixSpace(F7, 2, 3);
 *
 * // Create a matrix from a 2D array
 * const A = MS([[1, 2, 3], [4, 5, 6]]);
 *
 * // Create zero and identity matrices
 * const Z = MS.zero();
 *
 * // For square matrices
 * const MS2 = MatrixSpace(F7, 3, 3);
 * const I = MS2.identity();
 * ```
 */
export class MatrixSpaceClass<R extends RingElement> {
  readonly base_ring: CoefficientRing<R>;
  readonly nrows: number;
  readonly ncols: number;

  /**
   * Create a matrix space.
   *
   * @param base_ring - The ring of matrix entries
   * @param nrows - Number of rows
   * @param ncols - Number of columns (defaults to nrows for square matrices)
   */
  constructor(base_ring: CoefficientRing<R>, nrows: number, ncols?: number) {
    if (nrows < 0) {
      throw new ValueError('number of rows must be non-negative');
    }

    this.base_ring = base_ring;
    this.nrows = nrows;
    this.ncols = ncols ?? nrows;

    if (this.ncols < 0) {
      throw new ValueError('number of columns must be non-negative');
    }
  }

  /**
   * Create a matrix from entries.
   *
   * @param entries - Can be:
   *   - A 2D array of ring elements or coercible values
   *   - A flat array of ring elements or coercible values
   *   - A single value to fill the entire matrix
   *   - undefined for the zero matrix
   */
  __call__(entries?: unknown[][] | unknown[] | unknown): Matrix<R> {
    if (entries === undefined) {
      return this.zero();
    }

    // Single value - fill entire matrix
    if (!Array.isArray(entries)) {
      const val = this.base_ring.__call__(entries);
      const result = new Matrix<R>(this.base_ring, this.nrows, this.ncols);
      for (let i = 0; i < this.nrows; i++) {
        for (let j = 0; j < this.ncols; j++) {
          result.set(i, j, val);
        }
      }
      return result;
    }

    // 2D array
    if (entries.length > 0 && Array.isArray(entries[0])) {
      const entries2D = entries as unknown[][];
      if (entries2D.length !== this.nrows) {
        throw new ValueError(`expected ${this.nrows} rows, got ${entries2D.length}`);
      }

      const coerced: R[][] = [];
      for (let i = 0; i < this.nrows; i++) {
        if (entries2D[i]!.length !== this.ncols) {
          throw new ValueError(
            `expected ${this.ncols} columns in row ${i}, got ${entries2D[i]!.length}`
          );
        }
        coerced.push([]);
        for (let j = 0; j < this.ncols; j++) {
          coerced[i]!.push(this.base_ring.__call__(entries2D[i]![j]));
        }
      }

      return new Matrix<R>(this.base_ring, this.nrows, this.ncols, coerced);
    }

    // Flat array
    const flat = entries as unknown[];
    if (flat.length !== this.nrows * this.ncols) {
      throw new ValueError(`expected ${this.nrows * this.ncols} entries, got ${flat.length}`);
    }

    const coerced: R[] = flat.map((x) => this.base_ring.__call__(x));
    return new Matrix<R>(this.base_ring, this.nrows, this.ncols, coerced);
  }

  /**
   * Return the zero matrix.
   */
  zero(): Matrix<R> {
    return zero_matrix(this.base_ring, this.nrows, this.ncols);
  }

  /**
   * Return the identity matrix.
   *
   * @throws {ValueError} If the matrix space is not square
   */
  identity(): Matrix<R> {
    if (this.nrows !== this.ncols) {
      throw new ValueError(
        `identity matrix only exists for square matrix spaces (got ${this.nrows}x${this.ncols})`
      );
    }
    return identity_matrix(this.base_ring, this.nrows);
  }

  /**
   * Alias for identity().
   */
  one(): Matrix<R> {
    return this.identity();
  }

  /**
   * Check if this matrix space is square.
   */
  is_square(): boolean {
    return this.nrows === this.ncols;
  }

  /**
   * Return the dimension of this matrix space as a vector space.
   */
  dimension(): number {
    return this.nrows * this.ncols;
  }

  /**
   * Return the dimensions as a tuple.
   */
  dims(): [number, number] {
    return [this.nrows, this.ncols];
  }

  /**
   * Return a random matrix.
   *
   * This requires the base ring to have a random_element method.
   */
  random_element(): Matrix<R> {
    const ring = this.base_ring as CoefficientRing<R> & { random_element?: () => R };
    if (typeof ring.random_element !== 'function') {
      throw new ValueError('base ring does not support random element generation');
    }

    const result = new Matrix<R>(this.base_ring, this.nrows, this.ncols);
    for (let i = 0; i < this.nrows; i++) {
      for (let j = 0; j < this.ncols; j++) {
        result.set(i, j, ring.random_element());
      }
    }
    return result;
  }

  /**
   * String representation.
   */
  toString(): string {
    const density = 'dense';
    return `Full MatrixSpace of ${this.nrows} by ${this.ncols} ${density} matrices over ${this.base_ring}`;
  }
}

/**
 * Factory function to create a matrix space.
 *
 * @param ring - The base ring for matrix entries
 * @param nrows - Number of rows
 * @param ncols - Number of columns (defaults to nrows for square matrices)
 *
 * @example
 * ```typescript
 * const F7 = GF(7n);
 * const MS = MatrixSpace(F7, 2, 3);
 * const A = MS([[1, 2, 3], [4, 5, 6]]);
 * ```
 */
export function MatrixSpace<R extends RingElement>(
  ring: CoefficientRing<R>,
  nrows: number,
  ncols?: number
): MatrixSpaceClass<R> {
  return new MatrixSpaceClass(ring, nrows, ncols);
}

/**
 * Create a matrix from a ring and entries.
 *
 * This is a convenience function that creates a matrix directly without
 * explicitly constructing a matrix space.
 *
 * @param ring - The base ring for matrix entries
 * @param entries - 2D array of entries
 *
 * @example
 * ```typescript
 * const F7 = GF(7n);
 * const A = MatrixFromEntries(F7, [[1, 2], [3, 4]]);
 * ```
 */
export function MatrixFromEntries<R extends RingElement>(
  ring: CoefficientRing<R>,
  entries: unknown[][]
): Matrix<R> {
  if (entries.length === 0) {
    return new Matrix<R>(ring, 0, 0);
  }

  const nrows = entries.length;
  const ncols = entries[0]!.length;

  // Validate all rows have same length
  for (let i = 1; i < nrows; i++) {
    if (entries[i]!.length !== ncols) {
      throw new ValueError(
        `inconsistent row lengths: row 0 has ${ncols} entries, row ${i} has ${entries[i]!.length} entries`
      );
    }
  }

  const MS = MatrixSpace(ring, nrows, ncols);
  return MS.__call__(entries);
}

/**
 * Alias for MatrixFromEntries.
 */
export const matrix = MatrixFromEntries;
