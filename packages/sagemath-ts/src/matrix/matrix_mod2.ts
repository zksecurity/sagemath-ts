/**
 * @module sage/matrix/matrix_mod2
 * @description GF(2) matrix operations - important for cryptography
 *
 * Port of: sage/matrix/matrix_mod2_dense.pyx
 */

import { ArithmeticError, NotImplementedError, ValueError, ZeroDivisionError } from '../errors.js';

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

  /**
   * Create a GF(2) matrix.
   *
   * @param nrows - Number of rows
   * @param ncols - Number of columns
   * @param entries - Optional 2D array of entries (0 or 1)
   */
  constructor(nrows: number, ncols: number, entries?: (number | boolean)[][]) {
    this.nrows = nrows;
    this.ncols = ncols;

    // Initialize entries array
    this._entries = [];
    for (let i = 0; i < nrows; i++) {
      this._entries.push([]);
      for (let j = 0; j < ncols; j++) {
        if (entries !== undefined && entries[i] !== undefined && entries[i][j] !== undefined) {
          const val = entries[i][j];
          this._entries[i]!.push(typeof val === 'boolean' ? (val ? 1 : 0) : val! & 1);
        } else {
          this._entries[i]!.push(0);
        }
      }
    }
  }

  /**
   * Get entry at (i, j).
   *
   * @param i - Row index
   * @param j - Column index
   * @returns 0 or 1
   */
  get(i: number, j: number): number {
    if (i < 0 || i >= this.nrows || j < 0 || j >= this.ncols) {
      throw new ValueError(`index out of bounds: (${i}, ${j}) for ${this.nrows}x${this.ncols} matrix`);
    }
    return this._entries[i]![j]!;
  }

  /**
   * Set entry at (i, j).
   *
   * @param i - Row index
   * @param j - Column index
   * @param value - 0 or 1 (or any value that will be reduced mod 2)
   */
  set(i: number, j: number, value: number | boolean): void {
    if (i < 0 || i >= this.nrows || j < 0 || j >= this.ncols) {
      throw new ValueError(`index out of bounds: (${i}, ${j}) for ${this.nrows}x${this.ncols} matrix`);
    }
    this._entries[i]![j] = typeof value === 'boolean' ? (value ? 1 : 0) : value & 1;
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
   */
  str(
    rep_mapping?: Record<number, string>,
    zero?: string,
    plus_one?: string,
    minus_one?: string
  ): string {
    if (this.nrows === 0 || this.ncols === 0) {
      return '[]';
    }

    const zeroStr = rep_mapping?.[0] ?? zero ?? '0';
    const oneStr = rep_mapping?.[1] ?? plus_one ?? '1';

    const lines: string[] = [];
    for (let i = 0; i < this.nrows; i++) {
      const row: string[] = [];
      for (let j = 0; j < this.ncols; j++) {
        row.push(this._entries[i]![j] === 0 ? zeroStr : oneStr);
      }
      lines.push('[' + row.join(' ') + ']');
    }

    return lines.join('\n');
  }

  /**
   * Return a specific row.
   *
   * @param i - Row index
   * @param from_list - Whether to return from list
   * @returns The row as an array of 0s and 1s
   * @see Reference: sage/matrix/matrix_mod2_dense.pyx:row
   */
  row(i: number, from_list?: boolean): number[] {
    if (i < 0) {
      i = i + this.nrows;
    }
    if (i < 0 || i >= this.nrows) {
      throw new ValueError(`row index ${i} out of range`);
    }
    return [...this._entries[i]!];
  }

  /**
   * Return all columns.
   *
   * @param copy - Whether to copy
   * @returns List of columns (each column is an array of 0s and 1s)
   * @see Reference: sage/matrix/matrix_mod2_dense.pyx:columns
   */
  columns(copy?: boolean): number[][] {
    const result: number[][] = [];
    for (let j = 0; j < this.ncols; j++) {
      const col: number[] = [];
      for (let i = 0; i < this.nrows; i++) {
        col.push(this._entries[i]![j]!);
      }
      result.push(col);
    }
    return result;
  }

  /**
   * Matrix addition (XOR in GF(2)).
   *
   * @param other - Matrix to add
   * @returns Sum of matrices
   */
  add(other: Matrix_mod2_dense): Matrix_mod2_dense {
    if (this.nrows !== other.nrows || this.ncols !== other.ncols) {
      throw new ArithmeticError(
        `cannot add ${this.nrows}x${this.ncols} matrix to ${other.nrows}x${other.ncols} matrix`
      );
    }

    const result = new Matrix_mod2_dense(this.nrows, this.ncols);
    for (let i = 0; i < this.nrows; i++) {
      for (let j = 0; j < this.ncols; j++) {
        result._entries[i]![j] = this._entries[i]![j]! ^ other._entries[i]![j]!;
      }
    }
    return result;
  }

  /**
   * Matrix subtraction (same as addition in GF(2)).
   *
   * @param other - Matrix to subtract
   * @returns Difference of matrices
   */
  sub(other: Matrix_mod2_dense): Matrix_mod2_dense {
    return this.add(other);
  }

  /**
   * Classical matrix multiplication.
   *
   * @param right - The right operand
   * @returns The product
   * @see Reference: sage/matrix/matrix_mod2_dense.pyx:_multiply_classical
   */
  _multiply_classical(right: Matrix_mod2_dense): Matrix_mod2_dense {
    if (this.ncols !== right.nrows) {
      throw new ArithmeticError(
        `cannot multiply ${this.nrows}x${this.ncols} matrix by ${right.nrows}x${right.ncols} matrix`
      );
    }

    const result = new Matrix_mod2_dense(this.nrows, right.ncols);
    for (let i = 0; i < this.nrows; i++) {
      for (let j = 0; j < right.ncols; j++) {
        let sum = 0;
        for (let k = 0; k < this.ncols; k++) {
          sum ^= this._entries[i]![k]! & right._entries[k]![j]!;
        }
        result._entries[i]![j] = sum;
      }
    }
    return result;
  }

  /**
   * Matrix multiplication (alias for _multiply_classical).
   *
   * @param other - Matrix to multiply
   * @returns Product of matrices
   */
  mul(other: Matrix_mod2_dense): Matrix_mod2_dense {
    return this._multiply_classical(other);
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
   */
  inverse(): Matrix_mod2_dense {
    if (this.nrows !== this.ncols) {
      throw new ArithmeticError('inverse is only defined for square matrices');
    }

    const n = this.nrows;

    if (n === 0) {
      return this.copy();
    }

    if (this.rank() !== n) {
      throw new ZeroDivisionError('matrix is not invertible');
    }

    // Augment with identity matrix
    const augmented = new Matrix_mod2_dense(n, 2 * n);
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        augmented._entries[i]![j] = this._entries[i]![j]!;
      }
      augmented._entries[i]![n + i] = 1;
    }

    // Gaussian elimination
    for (let k = 0; k < n; k++) {
      // Find pivot
      let pivotRow = -1;
      for (let i = k; i < n; i++) {
        if (augmented._entries[i]![k] === 1) {
          pivotRow = i;
          break;
        }
      }

      if (pivotRow === -1) {
        throw new ZeroDivisionError('matrix is not invertible');
      }

      // Swap rows
      if (pivotRow !== k) {
        [augmented._entries[k], augmented._entries[pivotRow]] = [augmented._entries[pivotRow]!, augmented._entries[k]!];
      }

      // Eliminate other rows (no scaling needed in GF(2) since pivot is always 1)
      for (let i = 0; i < n; i++) {
        if (i !== k && augmented._entries[i]![k] === 1) {
          for (let j = 0; j < 2 * n; j++) {
            augmented._entries[i]![j] ^= augmented._entries[k]![j]!;
          }
        }
      }
    }

    // Extract inverse
    const result = new Matrix_mod2_dense(n, n);
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        result._entries[i]![j] = augmented._entries[i]![n + j]!;
      }
    }

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
   * @param algorithm - 'heuristic', 'mmpf', 'ple', 'm4ri', or 'classical'
   * @param cutoff - Cutoff for algorithm
   * @param reduced - Whether to compute reduced form (default: true)
   * @see Reference: sage/matrix/matrix_mod2_dense.pyx:echelonize
   */
  echelonize(algorithm?: string, cutoff?: number, reduced?: boolean): void {
    const doReduced = reduced !== false;

    let pivotRow = 0;

    for (let col = 0; col < this.ncols && pivotRow < this.nrows; col++) {
      // Find pivot in this column
      let found = -1;
      for (let i = pivotRow; i < this.nrows; i++) {
        if (this._entries[i]![col] === 1) {
          found = i;
          break;
        }
      }

      if (found === -1) {
        continue;
      }

      // Swap rows
      if (found !== pivotRow) {
        [this._entries[pivotRow], this._entries[found]] = [this._entries[found]!, this._entries[pivotRow]!];
      }

      // Eliminate (in GF(2), no scaling needed)
      const startRow = doReduced ? 0 : pivotRow + 1;
      for (let i = startRow; i < this.nrows; i++) {
        if (i !== pivotRow && this._entries[i]![col] === 1) {
          for (let j = 0; j < this.ncols; j++) {
            this._entries[i]![j] ^= this._entries[pivotRow]![j]!;
          }
        }
      }

      pivotRow++;
    }
  }

  /**
   * Return the echelon form (without modifying self).
   *
   * @param algorithm - Algorithm to use
   * @param cutoff - Cutoff parameter
   * @param reduced - Whether to compute reduced form
   * @returns The echelon form
   */
  echelon_form(algorithm?: string, cutoff?: number, reduced?: boolean): Matrix_mod2_dense {
    const result = this.copy();
    result.echelonize(algorithm, cutoff, reduced);
    return result;
  }

  /**
   * Return pivot columns.
   *
   * @returns List of pivot column indices
   * @see Reference: sage/matrix/matrix_mod2_dense.pyx:_pivots
   */
  pivots(): number[] {
    const echelon = this.echelon_form();
    const result: number[] = [];
    let pivotRow = 0;

    for (let col = 0; col < echelon.ncols && pivotRow < echelon.nrows; col++) {
      if (echelon._entries[pivotRow]![col] === 1) {
        result.push(col);
        pivotRow++;
      }
    }

    return result;
  }

  /**
   * Randomize entries.
   *
   * @param density - Density of ones (default: 1.0 for uniform random)
   * @param nonzero - Whether all entries should be nonzero (all ones)
   * @see Reference: sage/matrix/matrix_mod2_dense.pyx:randomize
   */
  randomize(density?: number, nonzero?: boolean): void {
    const d = density ?? 1.0;

    if (nonzero) {
      // All entries should be 1
      for (let i = 0; i < this.nrows; i++) {
        for (let j = 0; j < this.ncols; j++) {
          if (Math.random() <= d) {
            this._entries[i]![j] = 1;
          }
        }
      }
    } else {
      // Uniform random for selected entries
      if (d >= 1.0) {
        for (let i = 0; i < this.nrows; i++) {
          for (let j = 0; j < this.ncols; j++) {
            this._entries[i]![j] = Math.random() < 0.5 ? 0 : 1;
          }
        }
      } else {
        for (let i = 0; i < this.nrows; i++) {
          for (let j = 0; j < this.ncols; j++) {
            if (Math.random() <= d) {
              this._entries[i]![j] = Math.random() < 0.5 ? 0 : 1;
            }
          }
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
      throw new ArithmeticError('determinant is only defined for square matrices');
    }

    // Over GF(2), det = 1 iff matrix has full rank
    return this.rank() === this.nrows ? 1 : 0;
  }

  /**
   * Return the transpose.
   *
   * @returns The transpose
   * @see Reference: sage/matrix/matrix_mod2_dense.pyx:transpose
   */
  transpose(): Matrix_mod2_dense {
    const result = new Matrix_mod2_dense(this.ncols, this.nrows);
    for (let i = 0; i < this.nrows; i++) {
      for (let j = 0; j < this.ncols; j++) {
        result._entries[j]![i] = this._entries[i]![j]!;
      }
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
   */
  augment(right: Matrix_mod2_dense, subdivide?: boolean): Matrix_mod2_dense {
    if (this.nrows !== right.nrows) {
      throw new ArithmeticError('matrices must have the same number of rows');
    }

    const result = new Matrix_mod2_dense(this.nrows, this.ncols + right.ncols);
    for (let i = 0; i < this.nrows; i++) {
      for (let j = 0; j < this.ncols; j++) {
        result._entries[i]![j] = this._entries[i]![j]!;
      }
      for (let j = 0; j < right.ncols; j++) {
        result._entries[i]![this.ncols + j] = right._entries[i]![j]!;
      }
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
   */
  submatrix(row?: number, col?: number, nrows?: number, ncols?: number): Matrix_mod2_dense {
    const startRow = row ?? 0;
    const startCol = col ?? 0;
    const numRows = nrows ?? this.nrows - startRow;
    const numCols = ncols ?? this.ncols - startCol;

    if (startRow < 0 || startCol < 0 || startRow + numRows > this.nrows || startCol + numCols > this.ncols) {
      throw new ValueError('submatrix indices out of range');
    }

    const result = new Matrix_mod2_dense(numRows, numCols);
    for (let i = 0; i < numRows; i++) {
      for (let j = 0; j < numCols; j++) {
        result._entries[i]![j] = this._entries[startRow + i]![startCol + j]!;
      }
    }

    return result;
  }

  /**
   * Return the density (proportion of nonzero entries).
   *
   * @param approx - Whether to return approximate value
   * @returns The density of ones
   * @see Reference: sage/matrix/matrix_mod2_dense.pyx:density
   */
  density(approx?: boolean): number {
    if (this.nrows === 0 || this.ncols === 0) {
      return 0;
    }

    let count = 0;
    for (let i = 0; i < this.nrows; i++) {
      for (let j = 0; j < this.ncols; j++) {
        if (this._entries[i]![j] === 1) {
          count++;
        }
      }
    }

    return count / (this.nrows * this.ncols);
  }

  /**
   * Return the rank.
   *
   * @param algorithm - 'ple' or other
   * @returns The rank
   * @see Reference: sage/matrix/matrix_mod2_dense.pyx:rank
   */
  rank(algorithm?: string): number {
    // Compute echelon form and count non-zero rows
    const echelon = this.echelon_form();
    let r = 0;
    for (let i = 0; i < echelon.nrows; i++) {
      let isZeroRow = true;
      for (let j = 0; j < echelon.ncols; j++) {
        if (echelon._entries[i]![j] === 1) {
          isZeroRow = false;
          break;
        }
      }
      if (!isZeroRow) {
        r++;
      }
    }
    return r;
  }

  /**
   * Solve a linear system.
   *
   * @param B - Right-hand side
   * @param check - Whether to check the solution
   * @returns The solution
   * @see Reference: sage/matrix/matrix_mod2_dense.pyx:_solve_right_general
   */
  solve_right(B: Matrix_mod2_dense, check?: boolean): Matrix_mod2_dense {
    if (this.nrows !== this.ncols) {
      throw new ArithmeticError('solve_right requires a square matrix');
    }
    if (this.nrows !== B.nrows) {
      throw new ArithmeticError('incompatible matrix dimensions');
    }

    // Solve by computing inverse and multiplying
    const inv = this.inverse();
    const result = inv.mul(B);

    if (check) {
      const product = this.mul(result);
      if (!product.eq(B)) {
        throw new ArithmeticError('no solution exists');
      }
    }

    return result;
  }

  /**
   * Return a matrix whose rows form a basis for the right kernel.
   *
   * @returns Kernel matrix
   * @see Reference: sage/matrix/matrix_mod2_dense.pyx:_right_kernel_matrix
   */
  right_kernel_matrix(): Matrix_mod2_dense {
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

      // Find the column with lexicographically largest count vector
      let largestCol = 0;
      let largestCount = count1[0]!;
      for (let col = 1; col < i; col++) {
        const currentCount = count1[col]!;
        // Compare lexicographically
        for (let p = 0; p < partitionNum; p++) {
          if (currentCount[p]! > largestCount[p]!) {
            largestCol = col;
            largestCount = currentCount;
            break;
          } else if (currentCount[p]! < largestCount[p]!) {
            break;
          }
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
            A.swap_rows(rowStart, rowEnd);
            [rowSwapped[rowStart], rowSwapped[rowEnd]] = [rowSwapped[rowEnd]!, rowSwapped[rowStart]!];
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
      A.swap_columns(largestCol, i - 1);
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
  is_Gamma_free(certificate?: boolean): boolean | [boolean, [number, number, number, number] | null] {
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
   */
  swap_rows(i: number, j: number): void {
    if (i < 0 || i >= this.nrows || j < 0 || j >= this.nrows) {
      throw new ValueError(`row index out of bounds`);
    }
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
   */
  swap_columns(i: number, j: number): void {
    if (i < 0 || i >= this.ncols || j < 0 || j >= this.ncols) {
      throw new ValueError(`column index out of bounds`);
    }
    if (i !== j) {
      for (let row = 0; row < this.nrows; row++) {
        [this._entries[row]![i], this._entries[row]![j]] = [this._entries[row]![j]!, this._entries[row]![i]!];
      }
    }
  }

  /**
   * Permute rows according to a permutation.
   *
   * @param perm - Permutation (1-indexed, as in SageMath)
   * @see Reference: sage/matrix/matrix2.pyx:permute_rows
   */
  permute_rows(perm: number[]): void {
    if (perm.length !== this.nrows) {
      throw new ValueError('permutation length must match number of rows');
    }
    // Convert 1-indexed permutation to 0-indexed
    const newEntries: number[][] = [];
    for (let i = 0; i < this.nrows; i++) {
      newEntries.push(this._entries[perm[i]! - 1]!);
    }
    this._entries = newEntries;
  }

  /**
   * Permute columns according to a permutation.
   *
   * @param perm - Permutation (1-indexed, as in SageMath)
   * @see Reference: sage/matrix/matrix2.pyx:permute_columns
   */
  permute_columns(perm: number[]): void {
    if (perm.length !== this.ncols) {
      throw new ValueError('permutation length must match number of columns');
    }
    // Convert 1-indexed permutation to 0-indexed
    for (let row = 0; row < this.nrows; row++) {
      const newRow: number[] = [];
      for (let j = 0; j < this.ncols; j++) {
        newRow.push(this._entries[row]![perm[j]! - 1]!);
      }
      this._entries[row] = newRow;
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
 */
export function to_png_data(A: Matrix_mod2_dense): {
  width: number;
  height: number;
  pixels: Uint8Array;
} {
  const r = A.nrows;
  const c = A.ncols;

  if (r === 0 || c === 0) {
    throw new TypeError(`cannot create image with dimensions ${c} x ${r}`);
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
 * - LU is a matrix containing both L (lower) and U (upper) in packed form
 * - P is the row permutation (list of indices)
 * - Q is the column permutation (list of indices)
 *
 * The factorization satisfies: P * A * Q^T = L * U
 *
 * @param A - The matrix
 * @param algorithm - 'standard', 'mmpf', or 'naive' (default: 'standard')
 * @param param - Algorithm parameter (default: 0)
 * @returns Tuple [LU, P, Q]
 * @see Reference: sage/matrix/matrix_mod2_dense.pyx:pluq
 */
export function pluq(
  A: Matrix_mod2_dense,
  algorithm?: string,
  param?: number
): [Matrix_mod2_dense, number[], number[]] {
  const algo = algorithm ?? 'standard';
  if (algo !== 'standard' && algo !== 'mmpf' && algo !== 'naive') {
    throw new ValueError(`Algorithm '${algo}' unknown.`);
  }

  // Work on a copy
  const B = A.copy();
  const m = A.nrows;
  const n = A.ncols;

  // Initialize permutations as identity
  const P: number[] = Array.from({ length: m }, (_, i) => i);
  const Q: number[] = Array.from({ length: n }, (_, i) => i);

  // Gaussian elimination with row and column pivoting
  let pivotRow = 0;
  let pivotCol = 0;

  while (pivotRow < m && pivotCol < n) {
    // Find pivot (first 1 in submatrix starting at (pivotRow, pivotCol))
    let foundRow = -1;
    let foundCol = -1;

    outer: for (let j = pivotCol; j < n; j++) {
      for (let i = pivotRow; i < m; i++) {
        if (B._entries[i]![j] === 1) {
          foundRow = i;
          foundCol = j;
          break outer;
        }
      }
    }

    if (foundRow === -1) {
      // No pivot found, we're done
      break;
    }

    // Swap rows
    if (foundRow !== pivotRow) {
      B.swap_rows(pivotRow, foundRow);
      [P[pivotRow], P[foundRow]] = [P[foundRow]!, P[pivotRow]!];
    }

    // Swap columns
    if (foundCol !== pivotCol) {
      B.swap_columns(pivotCol, foundCol);
      [Q[pivotCol], Q[foundCol]] = [Q[foundCol]!, Q[pivotCol]!];
    }

    // Eliminate: for each row below pivot, if it has a 1 in pivot column, XOR with pivot row
    for (let i = pivotRow + 1; i < m; i++) {
      if (B._entries[i]![pivotCol] === 1) {
        // Store the multiplier in L part (below diagonal)
        // In GF(2), the multiplier is always 1 when we eliminate
        // The elimination: row[i] = row[i] XOR row[pivotRow]
        for (let j = pivotCol + 1; j < n; j++) {
          B._entries[i]![j] ^= B._entries[pivotRow]![j]!;
        }
        // The entry at (i, pivotCol) becomes part of L (it stays 1 to represent the multiplier)
        // Actually in PLUQ, we zero it out in U and store 1 implicitly in L
        // But in the packed format, we keep the lower part as-is
      }
    }

    pivotRow++;
    pivotCol++;
  }

  return [B, P, Q];
}

/**
 * Compute PLE factorization of A.
 *
 * Returns a tuple (LU, P, Q) where:
 * - LU is a matrix containing both L (lower) and E (echelon form) in packed form
 * - P is the row permutation (list of indices)
 * - Q is the pivot column positions (list of indices)
 *
 * The factorization is similar to PLUQ but with a different column strategy:
 * instead of column pivoting, it tracks which columns become pivot columns.
 *
 * @param A - The matrix
 * @param algorithm - 'standard', 'russian', or 'naive' (default: 'standard')
 * @param param - Algorithm parameter (default: 0)
 * @returns Tuple [LU, P, Q]
 * @see Reference: sage/matrix/matrix_mod2_dense.pyx:ple
 */
export function ple(
  A: Matrix_mod2_dense,
  algorithm?: string,
  param?: number
): [Matrix_mod2_dense, number[], number[]] {
  const algo = algorithm ?? 'standard';
  if (algo !== 'standard' && algo !== 'russian' && algo !== 'naive') {
    throw new ValueError(`Algorithm '${algo}' unknown.`);
  }

  // Work on a copy
  const B = A.copy();
  const m = A.nrows;
  const n = A.ncols;

  // Initialize row permutation as identity
  const P: number[] = Array.from({ length: m }, (_, i) => i);
  // Q will track pivot columns
  const Q: number[] = Array.from({ length: n }, (_, i) => i);

  let pivotRow = 0;

  for (let col = 0; col < n && pivotRow < m; col++) {
    // Find pivot in this column
    let foundRow = -1;
    for (let i = pivotRow; i < m; i++) {
      if (B._entries[i]![col] === 1) {
        foundRow = i;
        break;
      }
    }

    if (foundRow === -1) {
      // No pivot in this column, continue to next column
      // Update Q to reflect this column maps to itself
      continue;
    }

    // Swap rows
    if (foundRow !== pivotRow) {
      B.swap_rows(pivotRow, foundRow);
      [P[pivotRow], P[foundRow]] = [P[foundRow]!, P[pivotRow]!];
    }

    // Record pivot column position
    Q[pivotRow] = col;

    // Eliminate: for each row below pivot with a 1 in this column, XOR with pivot row
    for (let i = pivotRow + 1; i < m; i++) {
      if (B._entries[i]![col] === 1) {
        for (let j = col + 1; j < n; j++) {
          B._entries[i]![j] ^= B._entries[pivotRow]![j]!;
        }
        // Keep the 1 in the lower part to represent L
      }
    }

    pivotRow++;
  }

  return [B, P, Q];
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
 */
export function zero_matrix_gf2(nrows: number, ncols?: number): Matrix_mod2_dense {
  return new Matrix_mod2_dense(nrows, ncols ?? nrows);
}

/**
 * Create an identity matrix over GF(2).
 *
 * @param n - Size
 * @returns Identity matrix
 */
export function identity_matrix_gf2(n: number): Matrix_mod2_dense {
  const result = new Matrix_mod2_dense(n, n);
  for (let i = 0; i < n; i++) {
    result.set(i, i, 1);
  }
  return result;
}

/**
 * Create a random matrix over GF(2).
 *
 * @param nrows - Number of rows
 * @param ncols - Number of columns (default: nrows)
 * @param density - Density of ones (default: 0.5)
 * @returns Random matrix
 */
export function random_matrix_gf2(nrows: number, ncols?: number, density?: number): Matrix_mod2_dense {
  const result = new Matrix_mod2_dense(nrows, ncols ?? nrows);
  result.randomize(density ?? 1.0);
  return result;
}

/**
 * Create a matrix over GF(2) from entries.
 *
 * @param entries - 2D array of entries (0 or 1)
 * @returns The matrix
 */
export function matrix_gf2_from_entries(entries: (number | boolean)[][]): Matrix_mod2_dense {
  if (entries.length === 0) {
    return new Matrix_mod2_dense(0, 0);
  }

  const nrows = entries.length;
  const ncols = entries[0]!.length;

  // Validate row lengths
  for (let i = 1; i < nrows; i++) {
    if (entries[i]!.length !== ncols) {
      throw new ValueError(
        `inconsistent row lengths: row 0 has ${ncols} entries, row ${i} has ${entries[i]!.length} entries`
      );
    }
  }

  return new Matrix_mod2_dense(nrows, ncols, entries);
}
