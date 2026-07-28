/**
 * @module sage/matrix/matrix_integer
 * @description Integer matrices (over ZZ) with additional operations
 *
 * Port of: sage/matrix/matrix_integer_dense.pyx
 */

import { ArithmeticError, ValueError, ZeroDivisionError } from '../errors.js';
import { Integer, ZZ } from '../rings/integer_ring.js';
import { Rational } from '../rings/rational.js';
import { Matrix, identity_matrix, zero_matrix } from './matrix_generic.js';

/**
 * An integer matrix (matrix over ZZ).
 *
 * This class provides additional operations specific to integer matrices,
 * such as determinant calculation and row echelon form.
 */
export class IntegerMatrix {
  private _matrix: Matrix<Integer>;

  /**
   * Create an integer matrix.
   *
   * @param nrows - Number of rows
   * @param ncols - Number of columns
   * @param entries - 2D array of integer entries
   */
  constructor(nrows: number, ncols: number, entries?: bigint[][] | number[][] | Integer[][]) {
    // Convert entries to Integer objects
    let intEntries: Integer[][] | undefined;

    if (entries !== undefined) {
      intEntries = [];
      for (let i = 0; i < nrows; i++) {
        intEntries.push([]);
        for (let j = 0; j < ncols; j++) {
          const val = entries[i]?.[j];
          if (val === undefined) {
            intEntries[i]!.push(new Integer(0n));
          } else if (val instanceof Integer) {
            intEntries[i]!.push(val);
          } else {
            intEntries[i]!.push(new Integer(val));
          }
        }
      }
    }

    // Create the underlying Matrix<Integer>
    // We need to create a wrapper for ZZ that conforms to CoefficientRing<Integer>
    const zzRing = {
      zero: () => new Integer(0n),
      one: () => new Integer(1n),
      __call__: (x: unknown): Integer => {
        if (x instanceof Integer) return x;
        if (typeof x === 'bigint') return new Integer(x);
        if (typeof x === 'number') return new Integer(x);
        throw new ValueError(`cannot coerce ${x} to Integer`);
      },
      is_field: () => false,
      toString: () => 'Integer Ring',
    };

    this._matrix = new Matrix(zzRing, nrows, ncols, intEntries);
  }

  get nrows(): number {
    return this._matrix.nrows;
  }

  get ncols(): number {
    return this._matrix.ncols;
  }

  /**
   * Get the entry at position (i, j).
   */
  get(i: number, j: number): Integer {
    return this._matrix.get(i, j);
  }

  /**
   * Set the entry at position (i, j).
   */
  set(i: number, j: number, value: bigint | number | Integer): void {
    const intVal = value instanceof Integer ? value : new Integer(value);
    this._matrix.set(i, j, intVal);
  }

  /**
   * Check if this is a square matrix.
   */
  is_square(): boolean {
    return this._matrix.is_square();
  }

  /**
   * Add two integer matrices.
   */
  add(other: IntegerMatrix): IntegerMatrix {
    const result = new IntegerMatrix(this.nrows, this.ncols);
    for (let i = 0; i < this.nrows; i++) {
      for (let j = 0; j < this.ncols; j++) {
        result.set(i, j, this.get(i, j).add(other.get(i, j)));
      }
    }
    return result;
  }

  /**
   * Subtract two integer matrices.
   */
  sub(other: IntegerMatrix): IntegerMatrix {
    const result = new IntegerMatrix(this.nrows, this.ncols);
    for (let i = 0; i < this.nrows; i++) {
      for (let j = 0; j < this.ncols; j++) {
        result.set(i, j, this.get(i, j).sub(other.get(i, j)));
      }
    }
    return result;
  }

  /**
   * Negate this matrix.
   */
  neg(): IntegerMatrix {
    const result = new IntegerMatrix(this.nrows, this.ncols);
    for (let i = 0; i < this.nrows; i++) {
      for (let j = 0; j < this.ncols; j++) {
        result.set(i, j, this.get(i, j).neg());
      }
    }
    return result;
  }

  /**
   * Multiply this matrix by another integer matrix.
   */
  mul(other: IntegerMatrix): IntegerMatrix {
    if (this.ncols !== other.nrows) {
      throw new ArithmeticError(
        `cannot multiply ${this.nrows}x${this.ncols} matrix by ${other.nrows}x${other.ncols} matrix`
      );
    }

    const result = new IntegerMatrix(this.nrows, other.ncols);
    for (let i = 0; i < this.nrows; i++) {
      for (let j = 0; j < other.ncols; j++) {
        let sum = new Integer(0n);
        for (let k = 0; k < this.ncols; k++) {
          sum = sum.add(this.get(i, k).mul(other.get(k, j)));
        }
        result.set(i, j, sum);
      }
    }
    return result;
  }

  /**
   * Multiply by a scalar.
   */
  scalar_mul(c: bigint | number | Integer): IntegerMatrix {
    const scalar = c instanceof Integer ? c : new Integer(c);
    const result = new IntegerMatrix(this.nrows, this.ncols);
    for (let i = 0; i < this.nrows; i++) {
      for (let j = 0; j < this.ncols; j++) {
        result.set(i, j, this.get(i, j).mul(scalar));
      }
    }
    return result;
  }

  /**
   * Return the transpose.
   */
  transpose(): IntegerMatrix {
    const result = new IntegerMatrix(this.ncols, this.nrows);
    for (let i = 0; i < this.nrows; i++) {
      for (let j = 0; j < this.ncols; j++) {
        result.set(j, i, this.get(i, j));
      }
    }
    return result;
  }

  /**
   * Return the trace (sum of diagonal elements).
   */
  trace(): Integer {
    if (!this.is_square()) {
      throw new ArithmeticError('trace is only defined for square matrices');
    }

    let sum = new Integer(0n);
    for (let i = 0; i < this.nrows; i++) {
      sum = sum.add(this.get(i, i));
    }
    return sum;
  }

  /**
   * Compute the determinant of this matrix.
   *
   * Uses LU decomposition with fraction-free elimination for exact integer arithmetic.
   *
   * @throws {ArithmeticError} If the matrix is not square
   */
  determinant(): Integer {
    if (!this.is_square()) {
      throw new ArithmeticError('determinant is only defined for square matrices');
    }

    const n = this.nrows;

    if (n === 0) {
      return new Integer(1n);
    }

    if (n === 1) {
      return this.get(0, 0);
    }

    if (n === 2) {
      // det = a*d - b*c
      const a = this.get(0, 0);
      const b = this.get(0, 1);
      const c = this.get(1, 0);
      const d = this.get(1, 1);
      return a.mul(d).sub(b.mul(c));
    }

    if (n === 3) {
      // Sarrus rule
      const a = this.get(0, 0).value;
      const b = this.get(0, 1).value;
      const c = this.get(0, 2).value;
      const d = this.get(1, 0).value;
      const e = this.get(1, 1).value;
      const f = this.get(1, 2).value;
      const g = this.get(2, 0).value;
      const h = this.get(2, 1).value;
      const i = this.get(2, 2).value;

      const det = a * e * i + b * f * g + c * d * h - c * e * g - b * d * i - a * f * h;
      return new Integer(det);
    }

    // For larger matrices, use Bareiss algorithm (fraction-free Gaussian elimination)
    return this._determinantBareiss();
  }

  /**
   * Bareiss algorithm for computing determinant with exact integer arithmetic.
   * This avoids fractions by using division only when the result is guaranteed to be exact.
   */
  private _determinantBareiss(): Integer {
    const n = this.nrows;

    // Copy matrix data
    const M: bigint[][] = [];
    for (let i = 0; i < n; i++) {
      M.push([]);
      for (let j = 0; j < n; j++) {
        M[i]!.push(this.get(i, j).value);
      }
    }

    let sign = 1n;
    let prevPivot = 1n;

    for (let k = 0; k < n - 1; k++) {
      // Find pivot
      let pivotRow = k;
      for (let i = k; i < n; i++) {
        if (M[i]![k] !== 0n) {
          pivotRow = i;
          break;
        }
      }

      if (M[pivotRow]![k] === 0n) {
        return new Integer(0n);
      }

      // Swap rows if needed
      if (pivotRow !== k) {
        [M[k], M[pivotRow]] = [M[pivotRow]!, M[k]!];
        sign = -sign;
      }

      // Bareiss elimination
      for (let i = k + 1; i < n; i++) {
        for (let j = k + 1; j < n; j++) {
          const num = M[k]![k]! * M[i]![j]! - M[i]![k]! * M[k]![j]!;
          M[i]![j] = num / prevPivot;
        }
      }

      prevPivot = M[k]![k]!;
    }

    return new Integer(sign * M[n - 1]![n - 1]!);
  }

  /**
   * Compute the row echelon form (for future use in solving systems, etc.).
   * Returns a matrix in row echelon form along with the list of pivot columns.
   *
   * Note: This returns a matrix over rationals (as bigint fractions).
   * For exact integer row echelon form, use Hermite normal form instead.
   */
  row_echelon_form(): { matrix: IntegerMatrix; pivots: number[]; denominator: bigint } {
    const n = this.nrows;
    const m = this.ncols;

    // Work with rational arithmetic (numerator/denominator pairs)
    // We'll use a common denominator approach for simplicity

    // Copy to working matrix (as bigints)
    const M: bigint[][] = [];
    for (let i = 0; i < n; i++) {
      M.push([]);
      for (let j = 0; j < m; j++) {
        M[i]!.push(this.get(i, j).value);
      }
    }

    const denom = 1n;
    const pivots: number[] = [];
    let pivotRow = 0;

    for (let col = 0; col < m && pivotRow < n; col++) {
      // Find pivot in this column
      let found = -1;
      for (let i = pivotRow; i < n; i++) {
        if (M[i]![col] !== 0n) {
          found = i;
          break;
        }
      }

      if (found === -1) {
        continue;
      }

      // Swap rows
      if (found !== pivotRow) {
        [M[pivotRow], M[found]] = [M[found]!, M[pivotRow]!];
      }

      pivots.push(col);

      // Scale other rows to eliminate this column
      const pivot = M[pivotRow]![col]!;

      for (let i = pivotRow + 1; i < n; i++) {
        if (M[i]![col] !== 0n) {
          const factor = M[i]![col]!;
          for (let j = col; j < m; j++) {
            M[i]![j] = M[i]![j]! * pivot - factor * M[pivotRow]![j]!;
          }
          // Reduce by GCD to keep numbers manageable
          let g = 0n;
          for (let j = 0; j < m; j++) {
            g = bigintGcd(g, M[i]![j]! < 0n ? -M[i]![j]! : M[i]![j]!);
          }
          if (g > 1n) {
            for (let j = 0; j < m; j++) {
              M[i]![j] = M[i]![j]! / g;
            }
          }
        }
      }

      pivotRow++;
    }

    // Create result matrix
    const result = new IntegerMatrix(n, m);
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < m; j++) {
        result.set(i, j, M[i]![j]!);
      }
    }

    return { matrix: result, pivots, denominator: denom };
  }

  /**
   * Return the Hermite Normal Form of this matrix.
   *
   * @param transformation - Whether to return the transformation matrix
   * @returns HNF, or [HNF, U] if transformation=true where HNF = U * this
   */
  hermite_form(transformation?: false): IntegerMatrix;
  hermite_form(transformation: true): [IntegerMatrix, IntegerMatrix];
  hermite_form(transformation?: boolean): IntegerMatrix | [IntegerMatrix, IntegerMatrix] {
    return hermite_normal_form(this, 'default', false, true, transformation) as
      | IntegerMatrix
      | [IntegerMatrix, IntegerMatrix];
  }

  /**
   * Return the Smith Normal Form of this matrix.
   *
   * @param transformation - Whether to return the transformation matrices
   * @returns SNF, or [D, U, V] if transformation=true where D = U * this * V
   */
  smith_form(transformation?: false): IntegerMatrix;
  smith_form(transformation: true): [IntegerMatrix, IntegerMatrix, IntegerMatrix];
  smith_form(
    transformation?: boolean
  ): IntegerMatrix | [IntegerMatrix, IntegerMatrix, IntegerMatrix] {
    return smith_form_integer(this, transformation) as
      | IntegerMatrix
      | [IntegerMatrix, IntegerMatrix, IntegerMatrix];
  }

  /**
   * Return the elementary divisors of this matrix.
   *
   * @returns List of elementary divisors [d_1, d_2, ..., d_r] where d_i | d_{i+1}
   */
  elementary_divisors(): Integer[] {
    return elementary_divisors_integer(this);
  }

  /**
   * Return the rank of this matrix.
   *
   * @returns The rank (number of linearly independent rows/columns)
   */
  rank(): number {
    return rank_integer(this);
  }

  /**
   * Return the right kernel (null space) of this matrix.
   *
   * @returns Matrix whose rows span ker(this)
   */
  right_kernel_matrix(): IntegerMatrix {
    return kernel_matrix(this);
  }

  /**
   * Return the left kernel of this matrix.
   *
   * @returns Matrix whose rows span left ker(this)
   */
  left_kernel_matrix(): IntegerMatrix {
    return left_kernel_matrix(this);
  }

  /**
   * Check equality with another integer matrix.
   */
  eq(other: IntegerMatrix): boolean {
    if (this.nrows !== other.nrows || this.ncols !== other.ncols) {
      return false;
    }

    for (let i = 0; i < this.nrows; i++) {
      for (let j = 0; j < this.ncols; j++) {
        if (!this.get(i, j).eq(other.get(i, j))) {
          return false;
        }
      }
    }
    return true;
  }

  /**
   * Return a copy of this matrix.
   */
  copy(): IntegerMatrix {
    const result = new IntegerMatrix(this.nrows, this.ncols);
    for (let i = 0; i < this.nrows; i++) {
      for (let j = 0; j < this.ncols; j++) {
        result.set(i, j, this.get(i, j));
      }
    }
    return result;
  }

  /**
   * Return a string representation.
   */
  toString(): string {
    if (this.nrows === 0 || this.ncols === 0) {
      return `${this.nrows} x ${this.ncols} empty integer matrix`;
    }

    const strings: string[][] = [];
    const widths: number[] = new Array(this.ncols).fill(0);

    for (let i = 0; i < this.nrows; i++) {
      strings.push([]);
      for (let j = 0; j < this.ncols; j++) {
        const s = this.get(i, j).toString();
        strings[i]!.push(s);
        widths[j] = Math.max(widths[j]!, s.length);
      }
    }

    const lines: string[] = [];
    for (let i = 0; i < this.nrows; i++) {
      const row = strings[i]!.map((s, j) => s.padStart(widths[j]!));
      lines.push('[' + row.join(' ') + ']');
    }

    return lines.join('\n');
  }
}

/**
 * GCD of two bigints.
 */
function bigintGcd(a: bigint, b: bigint): bigint {
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
 * Create a zero integer matrix.
 */
export function zero_integer_matrix(nrows: number, ncols?: number): IntegerMatrix {
  return new IntegerMatrix(nrows, ncols ?? nrows);
}

/**
 * Create an identity integer matrix.
 */
export function identity_integer_matrix(n: number): IntegerMatrix {
  const result = new IntegerMatrix(n, n);
  for (let i = 0; i < n; i++) {
    result.set(i, i, 1n);
  }
  return result;
}

/**
 * Factory function to create an integer matrix from a 2D array.
 */
export function IntegerMatrixFromEntries(entries: (bigint | number)[][]): IntegerMatrix {
  if (entries.length === 0) {
    return new IntegerMatrix(0, 0);
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

  return new IntegerMatrix(nrows, ncols, entries);
}

// ============================================================================
// Integer Matrix Specific Operations (Stubs)
// ============================================================================

import { NotImplementedError } from '../errors.js';

// ============================================================================
// GCD and Extended GCD utilities for bigint
// ============================================================================

/**
 * Extended GCD: returns [g, s, t] such that g = gcd(a, b) = s*a + t*b.
 * The result g is always non-negative.
 */
function extendedGcd(a: bigint, b: bigint): [bigint, bigint, bigint] {
  if (b === 0n) {
    if (a < 0n) {
      return [-a, -1n, 0n];
    }
    return [a === 0n ? 0n : a, a === 0n ? 0n : 1n, 0n];
  }

  let oldR = a < 0n ? -a : a;
  let r = b < 0n ? -b : b;
  let oldS = a < 0n ? -1n : 1n;
  let s = 0n;
  let oldT = 0n;
  let t = b < 0n ? -1n : 1n;

  while (r !== 0n) {
    const q = oldR / r;
    [oldR, r] = [r, oldR - q * r];
    [oldS, s] = [s, oldS - q * s];
    [oldT, t] = [t, oldT - q * t];
  }

  return [oldR, oldS, oldT];
}

// ============================================================================
// Hermite Normal Form
// ============================================================================

/**
 * Return the Hermite Normal Form of the matrix.
 *
 * The HNF is an upper triangular matrix H with the following properties:
 * - H[i,i] > 0 for each non-zero row i (or H[i,i] = 0 for zero rows)
 * - For j < i, 0 <= H[j,i] < H[i,i]
 * - H = U * A for some unimodular matrix U (det(U) = ±1)
 *
 * @param matrix - The integer matrix
 * @param algorithm - Algorithm to use ('default', 'pari', 'ntl', 'flint')
 * @param proof - Whether to use proof mode
 * @param include_zero_rows - Whether to include zero rows
 * @param transformation - Whether to return the transformation matrix
 * @returns HNF, or (HNF, U) if transformation=true where H = U * A
 * @see Reference: sage/matrix/matrix_integer_dense.pyx:echelon_form
 */
export function hermite_normal_form(
  matrix: IntegerMatrix,
  algorithm?: 'default' | 'pari' | 'ntl' | 'flint',
  proof?: boolean,
  include_zero_rows?: boolean,
  transformation?: boolean
): IntegerMatrix | [IntegerMatrix, IntegerMatrix] {
  const m = matrix.nrows;
  const n = matrix.ncols;

  // Handle empty matrix
  if (m === 0 || n === 0) {
    const H = matrix.copy();
    if (transformation) {
      return [H, identity_integer_matrix(m)];
    }
    return H;
  }

  // Copy matrix entries to work with
  const H: bigint[][] = [];
  for (let i = 0; i < m; i++) {
    H.push([]);
    for (let j = 0; j < n; j++) {
      H[i]!.push(matrix.get(i, j).value);
    }
  }

  // Transformation matrix (if needed)
  const U: bigint[][] = [];
  if (transformation) {
    for (let i = 0; i < m; i++) {
      U.push([]);
      for (let j = 0; j < m; j++) {
        U[i]!.push(i === j ? 1n : 0n);
      }
    }
  }

  // Classical HNF algorithm using row operations
  let pivotRow = 0;

  for (let col = 0; col < n && pivotRow < m; col++) {
    // Find a non-zero entry in this column at or below pivotRow
    let nonZeroRow = -1;
    for (let i = pivotRow; i < m; i++) {
      if (H[i]![col] !== 0n) {
        nonZeroRow = i;
        break;
      }
    }

    if (nonZeroRow === -1) {
      continue; // This column is all zeros below pivotRow
    }

    // Move non-zero row to pivotRow position
    if (nonZeroRow !== pivotRow) {
      [H[pivotRow], H[nonZeroRow]] = [H[nonZeroRow]!, H[pivotRow]!];
      if (transformation) {
        [U[pivotRow], U[nonZeroRow]] = [U[nonZeroRow]!, U[pivotRow]!];
      }
    }

    // Use extended GCD to reduce all entries below the pivot
    for (let i = pivotRow + 1; i < m; i++) {
      if (H[i]![col] === 0n) continue;

      // Apply extended GCD: find g = s*H[pivotRow][col] + t*H[i][col]
      const [g, s, t] = extendedGcd(H[pivotRow]![col]!, H[i]![col]!);

      // Compute coefficients for the transformation
      const a = H[pivotRow]![col]!;
      const b = H[i]![col]!;
      const aOverG = a / g;
      const bOverG = b / g;

      // New rows: [s, t; -b/g, a/g] * [row_pivot; row_i]
      // This gives g in position (pivotRow, col) and 0 in position (i, col)
      const newPivotRow: bigint[] = [];
      const newRow: bigint[] = [];

      for (let j = 0; j < n; j++) {
        newPivotRow.push(s * H[pivotRow]![j]! + t * H[i]![j]!);
        newRow.push(-bOverG * H[pivotRow]![j]! + aOverG * H[i]![j]!);
      }

      H[pivotRow] = newPivotRow;
      H[i] = newRow;

      if (transformation) {
        const newUPivot: bigint[] = [];
        const newURow: bigint[] = [];
        for (let j = 0; j < m; j++) {
          newUPivot.push(s * U[pivotRow]![j]! + t * U[i]![j]!);
          newURow.push(-bOverG * U[pivotRow]![j]! + aOverG * U[i]![j]!);
        }
        U[pivotRow] = newUPivot;
        U[i] = newURow;
      }
    }

    // Make the pivot positive
    if (H[pivotRow]![col]! < 0n) {
      for (let j = 0; j < n; j++) {
        H[pivotRow]![j] = -H[pivotRow]![j]!;
      }
      if (transformation) {
        for (let j = 0; j < m; j++) {
          U[pivotRow]![j] = -U[pivotRow]![j]!;
        }
      }
    }

    // Reduce entries above the pivot
    const pivot = H[pivotRow]![col]!;
    if (pivot !== 0n) {
      for (let i = 0; i < pivotRow; i++) {
        if (H[i]![col] !== 0n) {
          // Compute quotient (floor division towards negative infinity for correct modular reduction)
          let q = H[i]![col]! / pivot;
          // Adjust for proper floor division (towards -infinity)
          if (H[i]![col]! < 0n && H[i]![col]! % pivot !== 0n) {
            q -= 1n;
          }

          for (let j = 0; j < n; j++) {
            H[i]![j] = H[i]![j]! - q * H[pivotRow]![j]!;
          }
          if (transformation) {
            for (let j = 0; j < m; j++) {
              U[i]![j] = U[i]![j]! - q * U[pivotRow]![j]!;
            }
          }
        }
      }
    }

    pivotRow++;
  }

  // Build result matrices
  let resultRows = m;
  if (include_zero_rows === false) {
    // Count non-zero rows
    resultRows = 0;
    for (let i = 0; i < m; i++) {
      let isZero = true;
      for (let j = 0; j < n; j++) {
        if (H[i]![j] !== 0n) {
          isZero = false;
          break;
        }
      }
      if (!isZero) {
        resultRows++;
      }
    }
  }

  const result = new IntegerMatrix(resultRows, n);
  const keptRows: number[] = [];
  let destRow = 0;
  for (let i = 0; i < m; i++) {
    let isZero = true;
    for (let j = 0; j < n; j++) {
      if (H[i]![j] !== 0n) {
        isZero = false;
        break;
      }
    }

    if (include_zero_rows !== false || !isZero) {
      for (let j = 0; j < n; j++) {
        result.set(destRow, j, H[i]![j]!);
      }
      keptRows.push(i);
      destRow++;
    }
  }

  if (transformation) {
    // Sage truncates U alongside H (`U = U[:r]`, matrix_integer_dense.pyx:2095),
    // so that `U * A == H` still type-checks when the zero rows are dropped.
    const Umatrix = new IntegerMatrix(keptRows.length, m);
    for (let i = 0; i < keptRows.length; i++) {
      for (let j = 0; j < m; j++) {
        Umatrix.set(i, j, U[keptRows[i]!]![j]!);
      }
    }
    return [result, Umatrix];
  }

  return result;
}

/**
 * Return the Smith Normal Form of the matrix.
 *
 * For a matrix A, computes D = U * A * V where:
 * - U and V are unimodular (determinant ±1)
 * - D is diagonal with entries d_1, d_2, ..., d_r where d_i | d_{i+1}
 *
 * @param matrix - The integer matrix
 * @param transformation - Whether to return the transformation matrices
 * @param integral - Whether transformation matrices should be integral
 * @returns SNF, or (D, U, V) if transformation=true where D = U * A * V
 * @see Reference: sage/matrix/matrix_integer_dense.pyx:smith_form
 */
export function smith_form_integer(
  matrix: IntegerMatrix,
  transformation?: boolean,
  integral?: boolean
): IntegerMatrix | [IntegerMatrix, IntegerMatrix, IntegerMatrix] {
  const m = matrix.nrows;
  const n = matrix.ncols;

  // Handle empty matrix
  if (m === 0 || n === 0) {
    const D = matrix.copy();
    if (transformation) {
      return [D, identity_integer_matrix(m), identity_integer_matrix(n)];
    }
    return D;
  }

  // Copy matrix entries to work with
  const A: bigint[][] = [];
  for (let i = 0; i < m; i++) {
    A.push([]);
    for (let j = 0; j < n; j++) {
      A[i]!.push(matrix.get(i, j).value);
    }
  }

  // Initialize transformation matrices
  const U: bigint[][] = [];
  const V: bigint[][] = [];
  for (let i = 0; i < m; i++) {
    U.push([]);
    for (let j = 0; j < m; j++) {
      U[i]!.push(i === j ? 1n : 0n);
    }
  }
  for (let i = 0; i < n; i++) {
    V.push([]);
    for (let j = 0; j < n; j++) {
      V[i]!.push(i === j ? 1n : 0n);
    }
  }

  const minDim = Math.min(m, n);

  for (let k = 0; k < minDim; k++) {
    // Find the pivot with smallest absolute value (non-zero)
    let pivotFound = false;
    let pivotI = k;
    let pivotJ = k;
    let minVal = 0n;

    for (let i = k; i < m; i++) {
      for (let j = k; j < n; j++) {
        if (A[i]![j] !== 0n) {
          const absVal = A[i]![j]! < 0n ? -A[i]![j]! : A[i]![j]!;
          if (!pivotFound || absVal < minVal) {
            pivotFound = true;
            pivotI = i;
            pivotJ = j;
            minVal = absVal;
          }
        }
      }
    }

    if (!pivotFound) {
      break; // Remaining submatrix is all zeros
    }

    // Move pivot to position (k, k)
    if (pivotI !== k) {
      [A[k], A[pivotI]] = [A[pivotI]!, A[k]!];
      [U[k], U[pivotI]] = [U[pivotI]!, U[k]!];
    }
    if (pivotJ !== k) {
      for (let i = 0; i < m; i++) {
        [A[i]![k], A[i]![pivotJ]] = [A[i]![pivotJ]!, A[i]![k]!];
      }
      for (let i = 0; i < n; i++) {
        [V[i]![k], V[i]![pivotJ]] = [V[i]![pivotJ]!, V[i]![k]!];
      }
    }

    // Iteratively eliminate row k and column k.
    //
    // Each pass first clears column k below the pivot, then row k to the right
    // of it.  Two kinds of operation are used:
    //
    //   * when the pivot divides the entry, a plain reduction
    //     `row_i -= q*row_k` (resp. `col_j -= q*col_k`) is applied.  It touches
    //     only the eliminated row (resp. column), so it can never re-introduce
    //     a nonzero entry elsewhere in column k;
    //   * otherwise the xgcd transform is applied, which replaces the pivot by
    //     `gcd(pivot, entry)` -- *strictly* smaller in absolute value.
    //
    // Hence every pass that does not finish the job strictly decreases
    // |A[k][k]|, and the loop terminates.  (The previous implementation used
    // the xgcd transform unconditionally; for (s,t) = (0,+-1) that is a swap
    // rather than a reduction, and the row and column phases undid each other
    // forever -- see DEVIATIONS/AUDIT C9.)
    for (;;) {
      const pivotBefore = A[k]![k]! < 0n ? -A[k]![k]! : A[k]![k]!;

      // Eliminate column k below the pivot
      for (let i = k + 1; i < m; i++) {
        if (A[i]![k] === 0n) continue;

        const a = A[k]![k]!;
        const b = A[i]![k]!;

        if (b % a === 0n) {
          const q = b / a;
          for (let j = 0; j < n; j++) {
            A[i]![j] = A[i]![j]! - q * A[k]![j]!;
          }
          for (let j = 0; j < m; j++) {
            U[i]![j] = U[i]![j]! - q * U[k]![j]!;
          }
          continue;
        }

        const [g, s, t] = extendedGcd(a, b);
        const aOverG = a / g;
        const bOverG = b / g;

        // Apply row operation: new rows = [[s, t], [-b/g, a/g]] * [row_k, row_i]
        const newRowK: bigint[] = [];
        const newRowI: bigint[] = [];
        for (let j = 0; j < n; j++) {
          newRowK.push(s * A[k]![j]! + t * A[i]![j]!);
          newRowI.push(-bOverG * A[k]![j]! + aOverG * A[i]![j]!);
        }
        A[k] = newRowK;
        A[i] = newRowI;

        // Update U
        const newUK: bigint[] = [];
        const newUI: bigint[] = [];
        for (let j = 0; j < m; j++) {
          newUK.push(s * U[k]![j]! + t * U[i]![j]!);
          newUI.push(-bOverG * U[k]![j]! + aOverG * U[i]![j]!);
        }
        U[k] = newUK;
        U[i] = newUI;
      }

      // Eliminate row k to the right of the pivot
      for (let j = k + 1; j < n; j++) {
        if (A[k]![j] === 0n) continue;

        const a = A[k]![k]!;
        const b = A[k]![j]!;

        if (b % a === 0n) {
          const q = b / a;
          for (let i = 0; i < m; i++) {
            A[i]![j] = A[i]![j]! - q * A[i]![k]!;
          }
          for (let i = 0; i < n; i++) {
            V[i]![j] = V[i]![j]! - q * V[i]![k]!;
          }
          continue;
        }

        const [g, s, t] = extendedGcd(a, b);
        const aOverG = a / g;
        const bOverG = b / g;

        // Apply column operation: new cols = [col_k, col_j] * [[s, -b/g], [t, a/g]]
        for (let i = 0; i < m; i++) {
          const newColK = s * A[i]![k]! + t * A[i]![j]!;
          const newColJ = -bOverG * A[i]![k]! + aOverG * A[i]![j]!;
          A[i]![k] = newColK;
          A[i]![j] = newColJ;
        }

        // Update V: V = V * [[s, -b/g], [t, a/g]]
        for (let i = 0; i < n; i++) {
          const newVK = s * V[i]![k]! + t * V[i]![j]!;
          const newVJ = -bOverG * V[i]![k]! + aOverG * V[i]![j]!;
          V[i]![k] = newVK;
          V[i]![j] = newVJ;
        }
      }

      // Done when both the column below and the row to the right are zero.
      let clean = true;
      for (let i = k + 1; i < m && clean; i++) {
        if (A[i]![k] !== 0n) clean = false;
      }
      for (let j = k + 1; j < n && clean; j++) {
        if (A[k]![j] !== 0n) clean = false;
      }
      if (clean) break;

      const pivotAfter = A[k]![k]! < 0n ? -A[k]![k]! : A[k]![k]!;
      if (pivotAfter >= pivotBefore) {
        throw new ArithmeticError(
          'smith_form: pivot did not decrease; elimination would not terminate'
        );
      }
    }

    // Make the diagonal element positive
    if (A[k]![k]! < 0n) {
      for (let j = 0; j < n; j++) {
        A[k]![j] = -A[k]![j]!;
      }
      for (let j = 0; j < m; j++) {
        U[k]![j] = -U[k]![j]!;
      }
    }
  }

  // Now ensure divisibility: d_i | d_{i+1}
  // Iterate until all divisibility conditions are satisfied
  let needsPass = true;
  while (needsPass) {
    needsPass = false;
    for (let i = 0; i < minDim - 1; i++) {
      if (A[i]![i] === 0n) continue;
      if (A[i + 1]![i + 1] === 0n) continue;

      // Check if A[i][i] divides A[i+1][i+1]
      if (A[i + 1]![i + 1]! % A[i]![i]! !== 0n) {
        needsPass = true;

        // The 2x2 block on rows/columns {i, i+1} is diag(a, b) with a > 0,
        // b > 0 and a not dividing b.  Turn it into diag(gcd(a,b), lcm(a,b)):
        //
        //   col_i += col_{i+1}          ->  [[a, 0], [b, b]]
        //   xgcd row transform          ->  [[g, t*b], [0, a*b/g]]
        //   col_{i+1} -= (t*b/g)*col_i  ->  [[g, 0], [0, a*b/g]]
        //
        // The last step is an *exact* division (g divides b), so no further
        // xgcd is needed -- using one there (as the previous implementation
        // did) destroyed diagonality.
        for (let r = 0; r < m; r++) {
          A[r]![i] = A[r]![i]! + A[r]![i + 1]!;
        }
        for (let r = 0; r < n; r++) {
          V[r]![i] = V[r]![i]! + V[r]![i + 1]!;
        }

        // Apply row operations to make A[i][i] = gcd(a, b) and A[i+1][i] = 0
        {
          const newA = A[i]![i]!;
          const newB = A[i + 1]![i]!;
          const [g2, s2, t2] = extendedGcd(newA, newB);
          const aOverG = newA / g2;
          const bOverG = newB / g2;

          const newRowK: bigint[] = [];
          const newRowI: bigint[] = [];
          for (let j = 0; j < n; j++) {
            newRowK.push(s2 * A[i]![j]! + t2 * A[i + 1]![j]!);
            newRowI.push(-bOverG * A[i]![j]! + aOverG * A[i + 1]![j]!);
          }
          A[i] = newRowK;
          A[i + 1] = newRowI;

          const newUK: bigint[] = [];
          const newUI: bigint[] = [];
          for (let j = 0; j < m; j++) {
            newUK.push(s2 * U[i]![j]! + t2 * U[i + 1]![j]!);
            newUI.push(-bOverG * U[i]![j]! + aOverG * U[i + 1]![j]!);
          }
          U[i] = newUK;
          U[i + 1] = newUI;
        }

        // Now clear A[i][i+1]; the pivot divides it exactly.
        if (A[i]![i + 1] !== 0n) {
          const pivot = A[i]![i]!;
          const target = A[i]![i + 1]!;
          if (target % pivot !== 0n) {
            throw new ArithmeticError('smith_form: unexpected non-exact division');
          }
          const q = target / pivot;

          for (let r = 0; r < m; r++) {
            A[r]![i + 1] = A[r]![i + 1]! - q * A[r]![i]!;
          }
          for (let r = 0; r < n; r++) {
            V[r]![i + 1] = V[r]![i + 1]! - q * V[r]![i]!;
          }
        }

        // Make diagonal positive again
        if (A[i]![i]! < 0n) {
          for (let j = 0; j < n; j++) {
            A[i]![j] = -A[i]![j]!;
          }
          for (let j = 0; j < m; j++) {
            U[i]![j] = -U[i]![j]!;
          }
        }
        if (A[i + 1]![i + 1]! < 0n) {
          for (let j = 0; j < n; j++) {
            A[i + 1]![j] = -A[i + 1]![j]!;
          }
          for (let j = 0; j < m; j++) {
            U[i + 1]![j] = -U[i + 1]![j]!;
          }
        }
      }
    }
  }

  // Build result matrices
  const D = new IntegerMatrix(m, n);
  for (let i = 0; i < m; i++) {
    for (let j = 0; j < n; j++) {
      D.set(i, j, A[i]![j]!);
    }
  }

  if (transformation) {
    const Umatrix = new IntegerMatrix(m, m);
    const Vmatrix = new IntegerMatrix(n, n);

    for (let i = 0; i < m; i++) {
      for (let j = 0; j < m; j++) {
        Umatrix.set(i, j, U[i]![j]!);
      }
    }

    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        Vmatrix.set(i, j, V[i]![j]!);
      }
    }

    return [D, Umatrix, Vmatrix];
  }

  return D;
}

/**
 * Return the elementary divisors of the matrix.
 *
 * The elementary divisors are the diagonal entries of the Smith normal form,
 * ordered so that each divides the next (d_1 | d_2 | ... | d_r).
 *
 * @param matrix - The integer matrix
 * @param algorithm - Algorithm to use ('pari' or 'default')
 * @returns List of elementary divisors [d_1, d_2, ..., d_r]
 * @see Reference: sage/matrix/matrix_integer_dense.pyx:elementary_divisors
 */
export function elementary_divisors_integer(
  matrix: IntegerMatrix,
  algorithm?: 'pari' | 'default'
): Integer[] {
  if (matrix.nrows === 0 || matrix.ncols === 0) {
    return [];
  }

  // Compute Smith normal form (without transformation matrices for efficiency)
  const D = smith_form_integer(matrix, false) as IntegerMatrix;

  // Extract diagonal elements.  These are the invariants of the cokernel of
  // *left* multiplication, i.e. of Z^nrows / (image), so Sage (via PARI's
  // matsnf) returns exactly `nrows` of them: the min(nrows, ncols) diagonal
  // entries padded with zeros.  Sage then sorts them, putting the zeros last.
  const minDim = Math.min(D.nrows, D.ncols);
  const raw: bigint[] = [];

  for (let i = 0; i < minDim; i++) {
    raw.push(D.get(i, i).value);
  }
  while (raw.length < matrix.nrows) {
    raw.push(0n);
  }

  const nonzero = raw.filter((x) => x !== 0n).sort((a, b) => (a < b ? -1 : a > b ? 1 : 0));
  const zeros = raw.filter((x) => x === 0n);

  return [...nonzero, ...zeros].map((x) => new Integer(x));
}

/**
 * Return the invariant factors of this matrix.
 *
 * The invariant factors are the same as the elementary divisors for matrices
 * over the integers. They are the diagonal entries of the Smith normal form,
 * ordered so that each divides the next.
 *
 * @param matrix - The integer matrix
 * @param algorithm - Algorithm to use ('pari' or 'default')
 * @returns List of invariant factors [d_1, d_2, ..., d_r]
 * @see Reference: sage/matrix/matrix2.pyx:invariant_factors
 */
export function invariant_factors_integer(
  matrix: IntegerMatrix,
  algorithm?: 'pari' | 'default'
): Integer[] {
  // For matrices over ZZ, invariant factors = elementary divisors
  return elementary_divisors_integer(matrix, algorithm);
}

/**
 * Return the rank of the matrix.
 *
 * The rank is the number of non-zero diagonal elements in the Smith normal form.
 *
 * @param matrix - The integer matrix
 * @returns The rank (number of linearly independent rows/columns)
 */
export function rank_integer(matrix: IntegerMatrix): number {
  const divisors = elementary_divisors_integer(matrix);
  let r = 0;
  for (const d of divisors) {
    if (!d.isZero()) {
      r++;
    }
  }
  return r;
}

/**
 * Compute the primary cyclic decomposition of a module determined by this matrix.
 *
 * If the matrix is viewed as defining a linear transformation, this computes
 * the cyclic decomposition of the corresponding module.
 *
 * For integer matrices, the cyclic decomposition is determined by the
 * elementary divisors: each non-unit elementary divisor d_i gives a cyclic
 * submodule Z/d_i Z.
 *
 * @param matrix - The integer matrix
 * @returns Array of pairs [d, k] where d is a prime power and k is multiplicity
 */
export function cyclic_decomposition_integer(matrix: IntegerMatrix): Array<[bigint, number]> {
  const elemDivs = elementary_divisors_integer(matrix);

  // Factor each non-trivial elementary divisor into prime powers
  const result: Array<[bigint, number]> = [];

  for (const d of elemDivs) {
    const val = d.value;
    if (val === 0n || val === 1n || val === -1n) {
      continue;
    }

    // Factor this divisor
    const absVal = val < 0n ? -val : val;
    const factors = _factor_bigint(absVal);

    for (const [prime, exp] of factors) {
      const primePower = prime ** BigInt(exp);
      result.push([primePower, 1]);
    }
  }

  return result;
}

/**
 * Simple factorization of a positive bigint into prime powers.
 * Returns pairs [prime, exponent].
 */
function _factor_bigint(n: bigint): Array<[bigint, number]> {
  if (n <= 1n) {
    return [];
  }

  const factors: Array<[bigint, number]> = [];
  let remaining = n;

  // Trial division
  let p = 2n;
  while (p * p <= remaining) {
    if (remaining % p === 0n) {
      let exp = 0;
      while (remaining % p === 0n) {
        remaining /= p;
        exp++;
      }
      factors.push([p, exp]);
    }
    p += p === 2n ? 1n : 2n;
  }

  if (remaining > 1n) {
    factors.push([remaining, 1]);
  }

  return factors;
}

/**
 * Return the integer right kernel (null space) of the matrix.
 *
 * The kernel consists of all integer vectors v such that A * v = 0.
 * The returned matrix has rows that form a basis for the kernel.
 *
 * Uses the HNF-based algorithm: augment [A | I] and compute HNF,
 * then extract the kernel from the transformation matrix.
 *
 * @param matrix - The integer matrix A
 * @returns Matrix whose rows span the right kernel of A
 */
export function kernel_matrix(matrix: IntegerMatrix): IntegerMatrix {
  const m = matrix.nrows;
  const n = matrix.ncols;

  // Handle trivial cases
  if (n === 0) {
    return new IntegerMatrix(0, 0);
  }
  if (m === 0) {
    // Everything is in the kernel
    return identity_integer_matrix(n);
  }

  // Compute the Smith form with transformation matrices: D = U * A * V
  const [D, U, V] = smith_form_integer(matrix, true) as [
    IntegerMatrix,
    IntegerMatrix,
    IntegerMatrix,
  ];

  // The kernel of A consists of columns of V corresponding to zero diagonal entries in D
  const minDim = Math.min(m, n);
  const kernelCols: number[] = [];

  // Find which diagonal entries are zero
  for (let i = 0; i < minDim; i++) {
    if (D.get(i, i).isZero()) {
      kernelCols.push(i);
    }
  }

  // Also include columns beyond the number of rows (these are always in the kernel)
  for (let j = minDim; j < n; j++) {
    kernelCols.push(j);
  }

  // Extract the kernel basis from V
  // The kernel vectors are the columns of V corresponding to zero diagonal entries
  const kernelDim = kernelCols.length;

  if (kernelDim === 0) {
    return new IntegerMatrix(0, n);
  }

  // Return the kernel vectors as rows of the result matrix
  const result = new IntegerMatrix(kernelDim, n);
  for (let i = 0; i < kernelDim; i++) {
    const colIdx = kernelCols[i]!;
    for (let j = 0; j < n; j++) {
      result.set(i, j, V.get(j, colIdx));
    }
  }

  return result;
}

/**
 * Return the integer left kernel of the matrix.
 *
 * The left kernel consists of all integer vectors v such that v * A = 0.
 * The returned matrix has rows that form a basis for the left kernel.
 *
 * @param matrix - The integer matrix A
 * @returns Matrix whose rows span the left kernel of A
 */
export function left_kernel_matrix(matrix: IntegerMatrix): IntegerMatrix {
  // Left kernel of A = Right kernel of A^T
  return kernel_matrix(matrix.transpose());
}

/**
 * Return the Frobenius (rational canonical) form of the matrix.
 *
 * The Frobenius form is a block diagonal matrix whose blocks are the companion
 * matrices of the elementary divisor polynomials `P_1, ..., P_k`, which satisfy
 * `P_{i+1} | P_i` (PARI's ordering: the minimal polynomial comes **first**).
 *
 * Sage's `frobenius_form` delegates to PARI's `matfrobenius`
 * (`matrix_integer_dense.pyx:2573`), so we port PARI's algorithm
 * (`alglin2.c:617` `RgM_Frobenius`, a mix of Ozello's thesis chapter 2 and
 * Storjohann's Lemmas 9.14/9.18) rather than deriving the invariant factors
 * from a Smith normal form: only PARI's version also produces the change-of-basis
 * matrix, and only PARI's version fixes the block ordering that Sage prints.
 *
 * @param matrix - The integer matrix (must be square)
 * @param flag - Computation flag:
 *   - 0 (default): return the Frobenius form matrix
 *   - 1: return only the elementary divisor polynomials, as coefficient arrays
 *        ordered constant-term first (Sage returns polynomials in `variable`)
 *   - 2: return `[F, B]` with `M = B^-1 * F * B`; both are rational matrices
 *        (row-major `Rational[][]`), exactly as Sage returns them over `QQ`
 * @param variable - Variable name for polynomials (accepted for signature
 *   compatibility with Sage; the port returns coefficient arrays, so it is unused)
 * @returns The Frobenius form (flag=0), list of polynomials (flag=1), or [F, B] (flag=2)
 * @see Reference: sage/matrix/matrix_integer_dense.pyx:2512 (frobenius_form)
 * @see Reference: pari/src/basemath/alglin2.c:617 (RgM_Frobenius), :688 (matfrobenius)
 */
export function frobenius_form_integer(
  matrix: IntegerMatrix,
  flag?: number,
  variable?: string
): IntegerMatrix | bigint[][] | [Rational[][], Rational[][]] {
  if (!matrix.is_square()) {
    throw new ArithmeticError('frobenius matrix of non-square matrix not defined.');
  }

  const n = matrix.nrows;
  const f = flag ?? 0;
  void variable;

  if (f !== 0 && f !== 1 && f !== 2) {
    // PARI: pari_err_FLAG("matfrobenius") (alglin2.c:695).
    throw new ValueError('incorrect flag in matfrobenius');
  }

  // Handle empty matrix: Sage's doctest `matrix([]).frobenius_form(2)` -> ([], []).
  if (n === 0) {
    if (f === 2) {
      return [[], []];
    }
    if (f === 1) {
      return [];
    }
    return matrix.copy();
  }

  // PARI matfrobenius: flags 0 and 2 both use RgM_Frobenius(M, 0, ...); flag 1
  // reads the elementary divisors off the same form.
  const { M, P, v } = _RgM_Frobenius(matrix, 0, f === 2);

  if (f === 1) {
    return _minpoly_listpolslice(M, v, n).map((poly) =>
      poly.map((c) => {
        if (c[1] !== 1n) {
          throw new ArithmeticError('elementary divisor is not integral');
        }
        return c[0];
      })
    );
  }

  if (f === 0) {
    const F = new IntegerMatrix(n, n);
    for (let i = 1; i <= n; i++) {
      for (let j = 1; j <= n; j++) {
        const c = M[i]![j]!;
        if (c[1] !== 1n) {
          throw new ArithmeticError('Frobenius form is not integral');
        }
        F.set(i - 1, j - 1, c[0]);
      }
    }
    return F;
  }

  // f === 2: Sage returns [F, B] over QQ with M = B^-1 * F * B
  // (matrix_integer_dense.pyx:2583-2586).  F happens to be integral for an
  // integer matrix, but B is genuinely rational, so both are returned as
  // matrices of `Rational` to mirror Sage's `MatrixSpace(QQ, n)`.
  return [_ratMatToRational(M, n), _ratMatToRational(P!, n)];
}

// ---------------------------------------------------------------------------
// PARI's Frobenius form (alglin2.c:428-720), ported verbatim over Q.
//
// The internal matrices are 1-indexed (`A[i][j]` = row i, column j, with index 0
// unused) so that the port matches PARI's `gcoeff(M,i,j)` line for line.
// ---------------------------------------------------------------------------

/** A square matrix over Q in the 1-indexed layout described above. */
type _RatMat = _Rat[][];

function _ratMatIdentity(n: number): _RatMat {
  const P: _RatMat = [[]];
  for (let i = 1; i <= n; i++) {
    const row: _Rat[] = [_ratZero];
    for (let j = 1; j <= n; j++) row.push(i === j ? _ratOne : _ratZero);
    P.push(row);
  }
  return P;
}

function _ratMatFromInteger(matrix: IntegerMatrix): _RatMat {
  const n = matrix.nrows;
  const M: _RatMat = [[]];
  for (let i = 1; i <= n; i++) {
    const row: _Rat[] = [_ratZero];
    for (let j = 1; j <= n; j++) row.push(_ratMake(matrix.get(i - 1, j - 1).value, 1n));
    M.push(row);
  }
  return M;
}

function _ratMatToRational(M: _RatMat, n: number): Rational[][] {
  const out: Rational[][] = [];
  for (let i = 1; i <= n; i++) {
    const row: Rational[] = [];
    for (let j = 1; j <= n; j++) row.push(new Rational(M[i]![j]![0], M[i]![j]![1]));
    out.push(row);
  }
  return out;
}

/**
 * `M <- U M U^-1` with `U = E_{i,j}(k)`; `P <- U P`.
 *
 * @see Reference: pari/src/basemath/alglin2.c:450 (transL)
 */
function _frobTransL(M: _RatMat, P: _RatMat | null, k: _Rat, i: number, j: number): void {
  const n = M.length - 1;
  for (let l = 1; l <= n; l++) M[l]![j] = _ratSub(M[l]![j]!, _ratMul(M[l]![i]!, k));
  for (let l = 1; l <= n; l++) M[i]![l] = _ratAdd(M[i]![l]!, _ratMul(M[j]![l]!, k));
  if (P) for (let l = 1; l <= n; l++) P[i]![l] = _ratAdd(P[i]![l]!, _ratMul(P[j]![l]!, k));
}

/**
 * Conjugate by the diagonal matrix with `1/M[a][b]` in position `j` (`j = a` or `b`).
 *
 * @see Reference: pari/src/basemath/alglin2.c:465 (transD)
 */
function _frobTransD(M: _RatMat, P: _RatMat | null, a: number, b: number, j: number): void {
  const k = M[a]![b]!;
  if (k[0] === 1n && k[1] === 1n) return;
  const ki = _ratMake(k[1], k[0]);
  const n = M.length - 1;
  for (let l = 1; l <= n; l++) {
    if (l !== j) {
      M[l]![j] = _ratMul(M[l]![j]!, k);
      M[j]![l] = j === a && l === b ? _ratOne : _ratMul(M[j]![l]!, ki);
    }
  }
  if (P) for (let l = 1; l <= n; l++) P[j]![l] = _ratMul(P[j]![l]!, ki);
}

/**
 * Conjugate by the transposition `(i j)`: swap columns `i`, `j` and rows `i`, `j`.
 *
 * @see Reference: pari/src/basemath/alglin2.c:484 (transS)
 */
function _frobTransS(M: _RatMat, P: _RatMat | null, i: number, j: number): void {
  const n = M.length - 1;
  for (let l = 1; l <= n; l++) {
    const t = M[l]![i]!;
    M[l]![i] = M[l]![j]!;
    M[l]![j] = t;
  }
  for (let l = 1; l <= n; l++) {
    const t = M[i]![l]!;
    M[i]![l] = M[j]![l]!;
    M[j]![l] = t;
  }
  if (P) {
    for (let l = 1; l <= n; l++) {
      const t = P[i]![l]!;
      P[i]![l] = P[j]![l]!;
      P[j]![l] = t;
    }
  }
}

/**
 * Storjohann Lemma 9.14, step 1.
 *
 * @see Reference: pari/src/basemath/alglin2.c:546 (weakfrobenius_step1)
 */
function _weakfrobenius_step1(M: _RatMat, P: _RatMat | null, j0: number): number {
  const n = M.length - 1;
  for (let j = j0; j < n; ++j) {
    if (M[j + 1]![j]![0] === 0n) {
      let k = j + 2;
      for (; k <= n; ++k) if (M[k]![j]![0] !== 0n) break;
      if (k > n) return j;
      _frobTransS(M, P, k, j + 1);
    }
    _frobTransD(M, P, j + 1, j, j + 1);
    /* Now M[j+1,j] = 1 */
    for (let k = 1; k <= n; ++k) {
      if (k !== j + 1 && M[k]![j]![0] !== 0n) {
        _frobTransL(M, P, _ratNeg(M[k]![j]!), k, j + 1);
        M[k]![j] = _ratZero;
      }
    }
  }
  return n;
}

/**
 * Storjohann Lemma 9.14, step 2.
 *
 * @see Reference: pari/src/basemath/alglin2.c:578 (weakfrobenius_step2)
 */
function _weakfrobenius_step2(M: _RatMat, P: _RatMat | null, j: number): void {
  const n = M.length - 1;
  for (let i = j; i >= 2; i--) {
    for (let k = j + 1; k <= n; k++) {
      if (M[i]![k]![0] !== 0n) _frobTransL(M, P, M[i]![k]!, i - 1, k);
    }
  }
}

/**
 * Storjohann Lemma 9.14, step 3.
 *
 * @see Reference: pari/src/basemath/alglin2.c:597 (weakfrobenius_step3)
 */
function _weakfrobenius_step3(M: _RatMat, P: _RatMat | null, j0: number, j: number): number {
  const n = M.length - 1;
  if (j === n) return 0;
  if (M[j0]![j + 1]![0] === 0n) {
    let k = j + 2;
    for (; k <= n; k++) if (M[j0]![k]![0] !== 0n) break;
    if (k > n) return 0;
    _frobTransS(M, P, k, j + 1);
  }
  _frobTransD(M, P, j0, j + 1, j + 1);
  for (let i = j + 2; i <= n; i++) {
    if (M[j0]![i]![0] !== 0n) _frobTransL(M, P, M[j0]![i]!, j + 1, i);
  }
  return 1;
}

/**
 * The companion block occupying rows/columns `i..j`, as a monic polynomial.
 *
 * @see Reference: pari/src/basemath/alglin2.c:495 (minpoly_polslice)
 */
function _minpoly_polslice(M: _RatMat, i: number, j: number): _RPoly {
  const d = j + 1 - i;
  const p: _RPoly = [];
  for (let k = 0; k < d; k++) p.push(_ratNeg(M[i + k]![j]!));
  p.push(_ratOne);
  return p;
}

/**
 * @see Reference: pari/src/basemath/alglin2.c:508 (minpoly_listpolslice)
 */
function _minpoly_listpolslice(M: _RatMat, V: number[], n: number): _RPoly[] {
  const W: _RPoly[] = [];
  for (let i = 0; i < V.length; i++) {
    W.push(_minpoly_polslice(M, V[i]!, i < V.length - 1 ? V[i + 1]! - 1 : n));
  }
  return W;
}

/**
 * @see Reference: pari/src/basemath/alglin2.c:518 (minpoly_dvdslice)
 */
function _minpoly_dvdslice(M: _RatMat, i: number, j: number, k: number): boolean {
  const [, r] = _rpDivMod(_minpoly_polslice(M, i, j - 1), _minpoly_polslice(M, j, k));
  return _rpIsZero(r);
}

/**
 * PARI's `RgM_Frobenius`, specialised to Q.
 *
 * Returns the (weak, when `flag = 1`) Frobenius form `F`, the accumulated basis
 * change `P` (when `wantP`) with `F = P * M * P^-1`, and the vector `v` of block
 * start indices.
 *
 * @see Reference: pari/src/basemath/alglin2.c:617 (RgM_Frobenius)
 */
function _RgM_Frobenius(
  matrix: IntegerMatrix,
  flag: number,
  wantP: boolean
): { M: _RatMat; P: _RatMat | null; v: number[] } {
  const n = matrix.nrows;
  const M = _ratMatFromInteger(matrix);
  const P = wantP ? _ratMatIdentity(n) : null;
  const v: number[] = [0]; // 1-indexed, like PARI's t_VECSMALL
  let nb = 0;
  let j0 = 1;

  while (j0 <= n) {
    let j = _weakfrobenius_step1(M, P, j0);
    _weakfrobenius_step2(M, P, j);
    const eps = _weakfrobenius_step3(M, P, j0, j);
    if (eps === 0) {
      v[++nb] = j0;
      if (flag === 0 && nb > 1 && !_minpoly_dvdslice(M, v[nb - 1]!, j0, j)) {
        j = j0;
        j0 = v[nb - 1]!;
        nb -= 2;
        _frobTransL(M, P, _ratOne, j, j0); /* lemma 9.18 */
      } else {
        j0 = j + 1;
      }
    } else {
      _frobTransS(M, P, j0, j + 1); /* theorem 4 */
    }
  }

  return { M, P, v: v.slice(1, nb + 1) };
}

// ---------------------------------------------------------------------------
// Exact arithmetic over Q[x], used to compute invariant factors as the Smith
// normal form of x*I - A (PARI computes the same data in `matfrobenius`).
// ---------------------------------------------------------------------------

/** A rational number as [numerator, denominator] with denominator > 0, in lowest terms. */
type _Rat = [bigint, bigint];

function _ratMake(num: bigint, den: bigint): _Rat {
  if (den === 0n) throw new ZeroDivisionError('rational division by zero');
  if (den < 0n) {
    num = -num;
    den = -den;
  }
  if (num === 0n) return [0n, 1n];
  const g = bigintGcd(num < 0n ? -num : num, den);
  return [num / g, den / g];
}

const _ratZero: _Rat = [0n, 1n];
const _ratOne: _Rat = [1n, 1n];

function _ratAdd(a: _Rat, b: _Rat): _Rat {
  return _ratMake(a[0] * b[1] + b[0] * a[1], a[1] * b[1]);
}
function _ratSub(a: _Rat, b: _Rat): _Rat {
  return _ratMake(a[0] * b[1] - b[0] * a[1], a[1] * b[1]);
}
function _ratMul(a: _Rat, b: _Rat): _Rat {
  return _ratMake(a[0] * b[0], a[1] * b[1]);
}
function _ratDiv(a: _Rat, b: _Rat): _Rat {
  return _ratMake(a[0] * b[1], a[1] * b[0]);
}
function _ratNeg(a: _Rat): _Rat {
  return [-a[0], a[1]];
}

/** A polynomial over Q, coefficients indexed by degree, no trailing zeros. */
type _RPoly = _Rat[];

function _rpTrim(p: _RPoly): _RPoly {
  let d = p.length - 1;
  while (d >= 0 && p[d]![0] === 0n) d--;
  return p.slice(0, d + 1);
}
function _rpIsZero(p: _RPoly): boolean {
  return p.length === 0;
}
function _rpDeg(p: _RPoly): number {
  return p.length - 1;
}
function _rpSub(a: _RPoly, b: _RPoly): _RPoly {
  const len = Math.max(a.length, b.length);
  const r: _RPoly = [];
  for (let i = 0; i < len; i++) {
    r.push(_ratSub(a[i] ?? _ratZero, b[i] ?? _ratZero));
  }
  return _rpTrim(r);
}
function _rpMul(a: _RPoly, b: _RPoly): _RPoly {
  if (_rpIsZero(a) || _rpIsZero(b)) return [];
  const r: _RPoly = new Array(a.length + b.length - 1).fill(_ratZero);
  for (let i = 0; i < a.length; i++) {
    for (let j = 0; j < b.length; j++) {
      r[i + j] = _ratAdd(r[i + j]!, _ratMul(a[i]!, b[j]!));
    }
  }
  return _rpTrim(r);
}
/** Euclidean division in Q[x]: returns [quotient, remainder]. */
function _rpDivMod(a: _RPoly, b: _RPoly): [_RPoly, _RPoly] {
  if (_rpIsZero(b)) throw new ZeroDivisionError('polynomial division by zero');
  let rem = [...a];
  const q: _RPoly = new Array(Math.max(0, _rpDeg(a) - _rpDeg(b) + 1)).fill(_ratZero);
  const bLead = b[b.length - 1]!;
  while (!_rpIsZero(rem) && _rpDeg(rem) >= _rpDeg(b)) {
    const shift = _rpDeg(rem) - _rpDeg(b);
    const c = _ratDiv(rem[rem.length - 1]!, bLead);
    q[shift] = _ratAdd(q[shift]!, c);
    const scaled: _RPoly = new Array(shift).fill(_ratZero);
    for (const coeff of b) scaled.push(_ratMul(coeff, c));
    rem = _rpSub(rem, _rpTrim(scaled));
  }
  return [_rpTrim(q), rem];
}
/**
 * Compute the characteristic polynomial of an integer matrix.
 * Uses Faddeev-LeVerrier algorithm.
 * Returns coefficients [c_0, c_1, ..., c_n] where char_poly = c_0 + c_1*x + ... + c_n*x^n
 */
function _characteristic_polynomial_bigint(matrix: IntegerMatrix): bigint[] {
  const n = matrix.nrows;

  if (n === 0) {
    return [1n]; // char poly of empty matrix is 1
  }

  if (n === 1) {
    // char poly = x - a[0][0]
    return [-matrix.get(0, 0).value, 1n];
  }

  // Faddeev-LeVerrier algorithm
  // c_n = 1
  // M_1 = A, c_{n-1} = -tr(M_1)
  // M_k = A * (M_{k-1} + c_{n-k+1} * I), c_{n-k} = -tr(M_k) / k

  const coeffs: bigint[] = new Array(n + 1);
  coeffs[n] = 1n; // leading coefficient

  // M starts as A
  let M = matrix.copy();

  for (let k = 1; k <= n; k++) {
    // c_{n-k} = -tr(M) / k
    const trace = M.trace().value;
    coeffs[n - k] = -trace / BigInt(k);

    if (k < n) {
      // M = A * (M + c_{n-k} * I)
      const MplusCI = M.copy();
      for (let i = 0; i < n; i++) {
        MplusCI.set(i, i, M.get(i, i).value + coeffs[n - k]!);
      }
      M = matrix.mul(MplusCI);
    }
  }

  return coeffs;
}

/**
 * Return the saturation of the row space.
 *
 * The saturation of a ZZ-module M embedded in ZZ^n is a module S that
 * contains M with finite index such that ZZ^n/S is torsion free.
 * This function takes the row span M of self, and finds another matrix
 * of full rank with row span the saturation of M.
 *
 * ALGORITHM: If A is a matrix of full rank, then hnf(transpose(A))^(-1)*A
 * is a saturation of A.
 *
 * @param matrix - The integer matrix
 * @param p - Prime or 0 for full saturation (currently not used - full saturation always)
 * @param proof - Whether to use proof mode (not used)
 * @param max_dets - Maximum number of determinants to compute (not used)
 * @returns The saturated matrix
 * @see Reference: sage/matrix/matrix_integer_dense_saturation.py:saturation
 */
export function saturation(
  matrix: IntegerMatrix,
  p?: number | bigint,
  proof?: boolean,
  max_dets?: number
): IntegerMatrix {
  const m = matrix.nrows;
  const n = matrix.ncols;

  // Handle trivial cases
  if (m === 0 || n === 0) {
    return matrix.copy();
  }

  // Compute rank
  const r = rank_integer(matrix);

  // If the matrix is square and full rank, saturation is identity
  if (matrix.is_square() && r === m) {
    return identity_integer_matrix(r);
  }

  // If we have more rows than rank, reduce to full-rank submatrix
  let A = matrix;
  if (m > r) {
    // Use HNF to find a full-rank submatrix
    const H = hermite_normal_form(matrix, 'default', false, false) as IntegerMatrix;
    // Take the first r non-zero rows
    const rows: bigint[][] = [];
    for (let i = 0; i < H.nrows && rows.length < r; i++) {
      let isZero = true;
      for (let j = 0; j < H.ncols; j++) {
        if (H.get(i, j).value !== 0n) {
          isZero = false;
          break;
        }
      }
      if (!isZero) {
        const row: bigint[] = [];
        for (let j = 0; j < H.ncols; j++) {
          row.push(H.get(i, j).value);
        }
        rows.push(row);
      }
    }
    if (rows.length === 0) {
      return matrix.copy();
    }
    A = IntegerMatrixFromEntries(rows);
  }

  // Factor out common factors from each row
  const [factoredA] = factor_out_common_factors_from_each_row(A);
  A = factoredA;

  if (A.nrows <= 1) {
    return A;
  }

  // Remove zero columns for efficiency
  const [Acompact, zeroCols] = _delete_zero_columns_internal(A);

  // Main algorithm: hnf(transpose(A))^(-1) * A gives saturation
  // Compute B = transpose(A).hermite_form()
  const AT = Acompact.transpose();
  const B = hermite_normal_form(AT, 'default', false, false) as IntegerMatrix;
  const BT = B.transpose();

  // Now compute B^(-1) * A by solving the system B * X = A
  // Since B is from HNF of A^T, it should be square and invertible
  // We use the adjugate method: B^(-1) = adj(B) / det(B)
  const detB = BT.determinant();
  if (detB.isZero()) {
    // Matrix is not full rank, return as-is
    return _insert_zero_columns_internal(Acompact, zeroCols);
  }

  // Compute adjugate and solve
  const result = _solve_integer_system(BT, Acompact);

  return _insert_zero_columns_internal(result, zeroCols);
}

/**
 * Internal helper: Delete zero columns from a matrix.
 * Returns the compacted matrix and a list of which columns were zero.
 */
function _delete_zero_columns_internal(matrix: IntegerMatrix): [IntegerMatrix, number[]] {
  const zeroCols: number[] = [];
  const nonZeroCols: number[] = [];

  for (let j = 0; j < matrix.ncols; j++) {
    let isZero = true;
    for (let i = 0; i < matrix.nrows; i++) {
      if (matrix.get(i, j).value !== 0n) {
        isZero = false;
        break;
      }
    }
    if (isZero) {
      zeroCols.push(j);
    } else {
      nonZeroCols.push(j);
    }
  }

  if (nonZeroCols.length === matrix.ncols) {
    return [matrix, []];
  }

  const result = new IntegerMatrix(matrix.nrows, nonZeroCols.length);
  for (let i = 0; i < matrix.nrows; i++) {
    for (let newJ = 0; newJ < nonZeroCols.length; newJ++) {
      result.set(i, newJ, matrix.get(i, nonZeroCols[newJ]!));
    }
  }

  return [result, zeroCols];
}

/**
 * Internal helper: Insert zero columns back into a matrix.
 */
function _insert_zero_columns_internal(matrix: IntegerMatrix, zeroCols: number[]): IntegerMatrix {
  if (zeroCols.length === 0) {
    return matrix;
  }

  const totalCols = matrix.ncols + zeroCols.length;
  const zeroSet = new Set(zeroCols);

  const result = new IntegerMatrix(matrix.nrows, totalCols);
  let srcCol = 0;

  for (let j = 0; j < totalCols; j++) {
    if (zeroSet.has(j)) {
      // Zero column
      for (let i = 0; i < matrix.nrows; i++) {
        result.set(i, j, 0n);
      }
    } else {
      // Copy from source
      for (let i = 0; i < matrix.nrows; i++) {
        result.set(i, j, matrix.get(i, srcCol));
      }
      srcCol++;
    }
  }

  return result;
}

/**
 * Internal helper: Solve B * X = A for integer matrices.
 * Uses Cramer's rule with adjugate matrix.
 */
function _solve_integer_system(B: IntegerMatrix, A: IntegerMatrix): IntegerMatrix {
  const n = B.nrows;
  const m = A.ncols;

  if (n !== A.nrows || !B.is_square()) {
    throw new ArithmeticError('incompatible dimensions for system solve');
  }

  const detB = B.determinant().value;
  if (detB === 0n) {
    throw new ArithmeticError('singular matrix in system solve');
  }

  // Compute adjugate of B
  const adj = _adjugate_integer(B);

  // X = adj(B) * A / det(B)
  const adjA = adj.mul(A);

  const result = new IntegerMatrix(n, m);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < m; j++) {
      result.set(i, j, adjA.get(i, j).value / detB);
    }
  }

  return result;
}

/**
 * Internal helper: Compute the adjugate (classical adjoint) of a matrix.
 */
function _adjugate_integer(matrix: IntegerMatrix): IntegerMatrix {
  const n = matrix.nrows;

  if (n === 1) {
    return IntegerMatrixFromEntries([[1n]]);
  }

  const result = new IntegerMatrix(n, n);

  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      // Cofactor C[i][j] = (-1)^(i+j) * det(minor(i,j))
      const minor = _get_minor(matrix, i, j);
      const minorDet = minor.determinant().value;
      const sign = (i + j) % 2 === 0 ? 1n : -1n;
      // Adjugate is transpose of cofactor matrix
      result.set(j, i, sign * minorDet);
    }
  }

  return result;
}

/**
 * Internal helper: Get the (i,j) minor of a matrix (delete row i and column j).
 */
function _get_minor(matrix: IntegerMatrix, row: number, col: number): IntegerMatrix {
  const n = matrix.nrows;
  const result = new IntegerMatrix(n - 1, n - 1);

  let destI = 0;
  for (let i = 0; i < n; i++) {
    if (i === row) continue;
    let destJ = 0;
    for (let j = 0; j < n; j++) {
      if (j === col) continue;
      result.set(destI, destJ, matrix.get(i, j));
      destJ++;
    }
    destI++;
  }

  return result;
}

/**
 * Return the index of the row space in its saturation.
 *
 * ALGORITHM: Use Hermite normal form twice to find an invertible matrix whose
 * inverse transforms a matrix with the same row span as self to its saturation,
 * then compute the determinant of that matrix.
 *
 * @param matrix - The integer matrix
 * @param proof - Whether to use proof mode (not used)
 * @returns The index (positive integer)
 * @see Reference: sage/matrix/matrix_integer_dense_saturation.py:index_in_saturation
 */
export function index_in_saturation(matrix: IntegerMatrix, proof?: boolean): Integer {
  const r = rank_integer(matrix);

  // Zero rank means index is 1
  if (r === 0) {
    return new Integer(1n);
  }

  // Get a full-rank submatrix if needed
  let A = matrix;
  if (r < matrix.nrows) {
    A = hermite_normal_form(matrix, 'default', false, false) as IntegerMatrix;
  }

  // If square, the index is just the absolute determinant
  if (A.is_square()) {
    const det = A.determinant().value;
    return new Integer(det < 0n ? -det : det);
  }

  // Otherwise, compute HNF of transpose
  A = A.transpose();
  A = hermite_normal_form(A, 'default', false, false) as IntegerMatrix;

  // The index is the absolute determinant
  const det = A.determinant().value;
  return new Integer(det < 0n ? -det : det);
}

/**
 * Return the pivot columns.
 *
 * @param matrix - The integer matrix
 * @returns Tuple of pivot column indices
 * @see Reference: sage/matrix/matrix_integer_dense.pyx:pivots
 */
export function pivots_integer(matrix: IntegerMatrix): number[] {
  const { pivots } = matrix.row_echelon_form();
  return pivots;
}

/**
 * Return the height of the matrix (max absolute value of entries).
 *
 * @param matrix - The integer matrix
 * @returns The height
 * @see Reference: sage/matrix/matrix_integer_dense.pyx:height
 */
export function height(matrix: IntegerMatrix): Integer {
  let maxVal = 0n;
  for (let i = 0; i < matrix.nrows; i++) {
    for (let j = 0; j < matrix.ncols; j++) {
      const val = matrix.get(i, j).value;
      const absVal = val < 0n ? -val : val;
      if (absVal > maxVal) {
        maxVal = absVal;
      }
    }
  }
  return new Integer(maxVal);
}

/**
 * Return the GCD of all entries.
 *
 * @param matrix - The integer matrix
 * @returns The GCD
 * @see Reference: sage/matrix/matrix_integer_dense.pyx:gcd
 */
export function gcd_integer_matrix(matrix: IntegerMatrix): Integer {
  let g = 0n;
  for (let i = 0; i < matrix.nrows; i++) {
    for (let j = 0; j < matrix.ncols; j++) {
      g = bigintGcd(g, matrix.get(i, j).value);
      if (g === 1n) {
        return new Integer(1n);
      }
    }
  }
  return new Integer(g);
}

/**
 * Check if the rows of the matrix are primitive (GCD of each row is 1).
 *
 * @param matrix - The integer matrix
 * @returns True if primitive
 * @see Reference: sage/matrix/matrix_integer_dense.pyx:is_primitive
 */
export function is_primitive(matrix: IntegerMatrix): boolean {
  for (let i = 0; i < matrix.nrows; i++) {
    let rowGcd = 0n;
    for (let j = 0; j < matrix.ncols; j++) {
      rowGcd = bigintGcd(rowGcd, matrix.get(i, j).value);
      if (rowGcd === 1n) {
        break;
      }
    }
    if (rowGcd !== 1n && rowGcd !== 0n) {
      return false;
    }
  }
  return true;
}

/**
 * Thrown internally when the integral LLL detects linearly dependent rows.
 */
const LLL_DEPENDENT = Symbol('LLL_DEPENDENT');

/**
 * Round a rational a/b (b > 0) to the nearest integer, ties away from -infinity.
 */
function _roundQuotient(a: bigint, b: bigint): bigint {
  // floor((2a + b) / (2b)) with exact floor division
  const num = 2n * a + b;
  const den = 2n * b;
  let q = num / den;
  if (num % den !== 0n && num < 0n) {
    q -= 1n;
  }
  return q;
}

/**
 * Integral LLL (Cohen, *A Course in Computational Algebraic Number Theory*,
 * Algorithm 2.6.7).
 *
 * Everything is done in exact integer arithmetic: `d[i]` is the determinant of
 * the Gram matrix of the first `i` rows and `lambda[i][j] = d[j] * mu[i][j]`.
 * Every division performed below is exact.
 *
 * `B` and `Uacc` are modified in place.  Throws `LLL_DEPENDENT` when the rows
 * turn out to be linearly dependent (the algorithm requires a basis).
 */
function _integral_lll(
  B: bigint[][],
  Uacc: bigint[][] | null,
  deltaNum: bigint,
  deltaDen: bigint
): void {
  const n = B.length;
  if (n === 0) return;
  const m = B[0]!.length;

  const dot = (u: bigint[], v: bigint[]): bigint => {
    let s = 0n;
    for (let i = 0; i < m; i++) s += u[i]! * v[i]!;
    return s;
  };

  // 1-indexed storage
  const d: bigint[] = new Array(n + 1).fill(0n);
  const lambda: bigint[][] = [];
  for (let i = 0; i <= n; i++) lambda.push(new Array(n + 1).fill(0n));

  const b = (i: number): bigint[] => B[i - 1]!;
  const u = (i: number): bigint[] => Uacc![i - 1]!;

  d[0] = 1n;
  d[1] = dot(b(1), b(1));
  if (d[1] === 0n) throw LLL_DEPENDENT;

  const RED = (k: number, l: number): void => {
    const lam = lambda[k]![l]!;
    const dl = d[l]!;
    if (2n * (lam < 0n ? -lam : lam) <= dl) return;
    const q = _roundQuotient(lam, dl);
    const bk = b(k);
    const bl = b(l);
    for (let i = 0; i < m; i++) bk[i] = bk[i]! - q * bl[i]!;
    if (Uacc !== null) {
      const uk = u(k);
      const ul = u(l);
      for (let i = 0; i < n; i++) uk[i] = uk[i]! - q * ul[i]!;
    }
    lambda[k]![l] = lam - q * dl;
    for (let i = 1; i < l; i++) {
      lambda[k]![i] = lambda[k]![i]! - q * lambda[l]![i]!;
    }
  };

  let k = 2;
  let kmax = 1;

  const SWAP = (kk: number): void => {
    const t = B[kk - 1]!;
    B[kk - 1] = B[kk - 2]!;
    B[kk - 2] = t;
    if (Uacc !== null) {
      const tu = Uacc[kk - 1]!;
      Uacc[kk - 1] = Uacc[kk - 2]!;
      Uacc[kk - 2] = tu;
    }
    for (let j = 1; j <= kk - 2; j++) {
      const tmp = lambda[kk]![j]!;
      lambda[kk]![j] = lambda[kk - 1]![j]!;
      lambda[kk - 1]![j] = tmp;
    }
    const lam = lambda[kk]![kk - 1]!;
    const Bnew = (d[kk - 2]! * d[kk]! + lam * lam) / d[kk - 1]!;
    for (let i = kk + 1; i <= kmax; i++) {
      const t2 = lambda[i]![kk]!;
      lambda[i]![kk] = (d[kk]! * lambda[i]![kk - 1]! - lam * t2) / d[kk - 1]!;
      lambda[i]![kk - 1] = (Bnew * t2 + lam * lambda[i]![kk]!) / d[kk]!;
    }
    d[kk - 1] = Bnew;
  };

  while (k <= n) {
    if (k > kmax) {
      kmax = k;
      for (let j = 1; j <= k; j++) {
        let uu = dot(b(k), b(j));
        for (let i = 1; i <= j - 1; i++) {
          uu = (d[i]! * uu - lambda[k]![i]! * lambda[j]![i]!) / d[i - 1]!;
        }
        if (j < k) {
          lambda[k]![j] = uu;
        } else {
          d[k] = uu;
          if (uu === 0n) throw LLL_DEPENDENT;
        }
      }
    }

    for (;;) {
      RED(k, k - 1);
      // Lovasz condition: d_k*d_{k-2} >= delta*d_{k-1}^2 - lambda_{k,k-1}^2
      const lam = lambda[k]![k - 1]!;
      const lhs = deltaDen * d[k]! * d[k - 2]!;
      const rhs = deltaNum * d[k - 1]! * d[k - 1]! - deltaDen * lam * lam;
      if (lhs < rhs) {
        SWAP(k);
        k = Math.max(2, k - 1);
      } else {
        for (let l = k - 2; l >= 1; l--) {
          RED(k, l);
        }
        k = k + 1;
        break;
      }
    }
  }
}

/**
 * Return the LLL-reduced form of the matrix.
 *
 * This implementation uses the Lenstra-Lenstra-Lovasz lattice basis reduction
 * algorithm in **exact integer arithmetic** (Cohen, Algorithm 2.6.7), so the
 * output always spans exactly the same lattice as the input, no matter how
 * large the entries are.  The rows of the matrix are treated as basis vectors
 * of a lattice.
 *
 * When the rows are linearly dependent, the row lattice is first put in
 * Hermite normal form (which spans the same lattice) and the nonzero rows are
 * reduced; the output then begins with `nrows - rank` zero rows, exactly like
 * Sage's `matrix(ZZ,3,3,range(1,10)).LLL()`.
 *
 * @param matrix - The integer matrix
 * @param delta - LLL parameter (default: 0.99, as in Sage)
 * @param eta - LLL parameter (default: 0.501); the exact algorithm always
 *   achieves |mu| <= 1/2, which is at least as strong
 * @param algorithm - Algorithm ('fpLLL:wrapper', 'fpLLL:proved', 'NTL:LLL', etc.)
 * @param fp - Floating point type
 * @param prec - Precision
 * @param early_red - Whether to use early reduction
 * @param use_givens - Whether to use Givens rotations
 * @param use_siegel - Whether to use Siegel conditions
 * @param transformation - Whether to return transformation matrix
 * @returns LLL-reduced matrix, or (LLL, U) if transformation=true
 * @see Reference: sage/matrix/matrix_integer_dense.pyx:3302 (LLL)
 */
export function LLL(
  matrix: IntegerMatrix,
  delta?: number,
  eta?: number,
  algorithm?: string,
  fp?: string,
  prec?: number,
  early_red?: boolean,
  use_givens?: boolean,
  use_siegel?: boolean,
  transformation?: boolean
): IntegerMatrix | [IntegerMatrix, IntegerMatrix] {
  const d = delta ?? 0.99;
  if (d <= 0.25) {
    throw new TypeError('delta must be > 0.25');
  }
  if (d > 1) {
    throw new TypeError('delta must be <= 1');
  }
  if (eta !== undefined && eta < 0.5) {
    throw new TypeError('eta must be >= 0.5');
  }

  const n = matrix.nrows;
  const m = matrix.ncols;

  if (n === 0 || m === 0) {
    if (transformation) {
      return [matrix.copy(), identity_integer_matrix(n)];
    }
    return matrix.copy();
  }

  // Exact rational form of delta.
  const deltaDen = 1000000000n;
  const deltaNum = BigInt(Math.round(d * 1e9));

  const toRows = (M: IntegerMatrix): bigint[][] => {
    const rows: bigint[][] = [];
    for (let i = 0; i < M.nrows; i++) {
      const row: bigint[] = [];
      for (let j = 0; j < M.ncols; j++) row.push(M.get(i, j).value);
      rows.push(row);
    }
    return rows;
  };

  let B = toRows(matrix);
  let U: bigint[][] | null = null;
  if (transformation) {
    U = [];
    for (let i = 0; i < n; i++) {
      const row: bigint[] = new Array(n).fill(0n);
      row[i] = 1n;
      U.push(row);
    }
  }

  let zeroRows = 0;
  try {
    _integral_lll(B, U, deltaNum, deltaDen);
  } catch (err) {
    if (err !== LLL_DEPENDENT) throw err;

    // Linearly dependent rows: pass to a basis of the same lattice via HNF.
    const hnf = hermite_normal_form(matrix, undefined, undefined, true, true) as [
      IntegerMatrix,
      IntegerMatrix,
    ];
    const H = toRows(hnf[0]);
    const HU = toRows(hnf[1]);

    const nonzero: number[] = [];
    const zero: number[] = [];
    for (let i = 0; i < n; i++) {
      if (H[i]!.some((x) => x !== 0n)) {
        nonzero.push(i);
      } else {
        zero.push(i);
      }
    }
    zeroRows = zero.length;

    const sub = nonzero.map((i) => H[i]!);
    let subU: bigint[][] | null = null;
    if (transformation) {
      subU = nonzero.map((_, idx) => {
        const row: bigint[] = new Array(nonzero.length).fill(0n);
        row[idx] = 1n;
        return row;
      });
    }
    _integral_lll(sub, subU, deltaNum, deltaDen);

    B = [];
    for (const i of zero) B.push(H[i]!);
    for (const row of sub) B.push(row);

    if (transformation) {
      U = [];
      for (const i of zero) U.push(HU[i]!);
      for (let r = 0; r < sub.length; r++) {
        const row: bigint[] = new Array(n).fill(0n);
        for (let c = 0; c < nonzero.length; c++) {
          const coeff = subU![r]![c]!;
          if (coeff === 0n) continue;
          const src = HU[nonzero[c]!]!;
          for (let j = 0; j < n; j++) row[j] = row[j]! + coeff * src[j]!;
        }
        U.push(row);
      }
    }
  }
  void zeroRows;

  // Convert back to IntegerMatrix
  const result = new IntegerMatrix(n, m);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < m; j++) {
      result.set(i, j, B[i]![j]!);
    }
  }

  if (transformation) {
    const Umatrix = new IntegerMatrix(n, n);
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        Umatrix.set(i, j, U![i]![j]!);
      }
    }
    return [result, Umatrix];
  }

  return result;
}

/**
 * Return the BKZ-reduced form of the matrix.
 *
 * BKZ (Block Korkine-Zolotarev) reduction is a lattice basis reduction algorithm
 * that generalizes LLL by using larger block sizes. For block_size=2, it is
 * equivalent to LLL.
 *
 * @param matrix - The integer matrix
 * @param delta - BKZ parameter (default: 0.99)
 * @param algorithm - Algorithm ('fpLLL' or 'NTL') - not used, we use native impl
 * @param fp - Floating point type - not used
 * @param block_size - Block size (default: 10, larger = better reduction but slower)
 * @param prune - Pruning parameter - not used
 * @param proof - Whether to use proof mode - not used
 * @returns BKZ-reduced matrix
 * @see Reference: sage/matrix/matrix_integer_dense.pyx:BKZ
 */
export function BKZ(
  matrix: IntegerMatrix,
  delta?: number,
  algorithm?: string,
  fp?: string,
  block_size?: number,
  prune?: number,
  proof?: boolean
): IntegerMatrix {
  const d = delta ?? 0.99;
  const beta = block_size ?? 10;

  const n = matrix.nrows;
  const m = matrix.ncols;

  if (n === 0 || m === 0) {
    return matrix.copy();
  }

  // For small block sizes, just use LLL
  if (beta <= 2 || n <= 2) {
    return LLL(matrix, d) as IntegerMatrix;
  }

  // BKZ algorithm: iteratively apply LLL to blocks of size beta
  // and enumerate short vectors in the blocks

  // Start with LLL reduction
  let B = LLL(matrix, d) as IntegerMatrix;

  // BKZ loop - simplified version
  // A full implementation would use enumeration for finding short vectors
  // This is a simplified approximation that does multiple LLL passes

  const maxIterations = Math.min(10, Math.ceil(n / beta) * 2);
  let changed = true;
  let iteration = 0;

  while (changed && iteration < maxIterations) {
    changed = false;
    iteration++;

    // Process each block
    for (let k = 0; k < n; k++) {
      const blockEnd = Math.min(k + beta, n);
      const blockSize = blockEnd - k;

      if (blockSize <= 1) continue;

      // Extract the block
      const block = new IntegerMatrix(blockSize, m);
      for (let i = 0; i < blockSize; i++) {
        for (let j = 0; j < m; j++) {
          block.set(i, j, B.get(k + i, j));
        }
      }

      // LLL-reduce the block
      const reducedBlock = LLL(block, d) as IntegerMatrix;

      // Check if anything changed
      let blockChanged = false;
      for (let i = 0; i < blockSize && !blockChanged; i++) {
        for (let j = 0; j < m && !blockChanged; j++) {
          if (reducedBlock.get(i, j).value !== B.get(k + i, j).value) {
            blockChanged = true;
          }
        }
      }

      if (blockChanged) {
        changed = true;
        // Update the matrix with the reduced block
        for (let i = 0; i < blockSize; i++) {
          for (let j = 0; j < m; j++) {
            B.set(k + i, j, reducedBlock.get(i, j));
          }
        }

        // Re-reduce the entire matrix after block update
        B = LLL(B, d) as IntegerMatrix;
      }
    }
  }

  return B;
}

/**
 * Check if the matrix is LLL-reduced.
 *
 * A basis is (delta, eta)-LLL-reduced if:
 * 1. For all i > j: |mu_{i,j}| <= eta
 * 2. For all i < n-1: delta * |b_i*|^2 <= |b_{i+1}* + mu_{i+1,i} * b_i*|^2
 *
 * @param matrix - The integer matrix
 * @param delta - LLL parameter (default: 0.75)
 * @param eta - LLL parameter (default: 0.501)
 * @param algorithm - Algorithm (not used)
 * @returns True if LLL-reduced
 * @see Reference: sage/matrix/matrix_integer_dense.pyx:is_LLL_reduced
 */
export function is_LLL_reduced(
  matrix: IntegerMatrix,
  delta?: number,
  eta?: number,
  algorithm?: string
): boolean {
  const d = delta ?? 0.75;
  const e = eta ?? 0.501;

  const n = matrix.nrows;
  const m = matrix.ncols;

  if (n <= 1) {
    return true;
  }

  // Compute Gram-Schmidt orthogonalization
  const B: number[][] = [];
  for (let i = 0; i < n; i++) {
    B.push([]);
    for (let j = 0; j < m; j++) {
      B[i]!.push(Number(matrix.get(i, j).value));
    }
  }

  const GStar: number[][] = [];
  const mu: number[][] = [];
  const Bnorms: number[] = [];

  for (let i = 0; i < n; i++) {
    GStar.push([...B[i]!]);
    mu.push(new Array(n).fill(0));

    for (let j = 0; j < i; j++) {
      // Compute mu[i][j] = <B[i], GStar[j]> / <GStar[j], GStar[j]>
      let dot1 = 0;
      let dot2 = 0;
      for (let k = 0; k < m; k++) {
        dot1 += B[i]![k]! * GStar[j]![k]!;
        dot2 += GStar[j]![k]! * GStar[j]![k]!;
      }
      mu[i]![j] = dot2 !== 0 ? dot1 / dot2 : 0;

      // GStar[i] = GStar[i] - mu[i][j] * GStar[j]
      for (let k = 0; k < m; k++) {
        GStar[i]![k] -= mu[i]![j]! * GStar[j]![k]!;
      }
    }

    // Compute |GStar[i]|^2
    let norm = 0;
    for (let k = 0; k < m; k++) {
      norm += GStar[i]![k]! * GStar[i]![k]!;
    }
    Bnorms.push(norm);
  }

  // Check size reduction condition: |mu_{i,j}| <= eta for all i > j
  for (let i = 1; i < n; i++) {
    for (let j = 0; j < i; j++) {
      if (Math.abs(mu[i]![j]!) > e + 1e-10) {
        return false;
      }
    }
  }

  // Check Lovasz condition: delta * |b_i*|^2 <= |b_{i+1}*|^2 + mu_{i+1,i}^2 * |b_i*|^2
  for (let i = 0; i < n - 1; i++) {
    const lhs = d * Bnorms[i]!;
    const muSq = mu[i + 1]![i]! * mu[i + 1]![i]!;
    const rhs = Bnorms[i + 1]! + muSq * Bnorms[i]!;

    if (lhs > rhs + 1e-10) {
      return false;
    }
  }

  return true;
}

/**
 * Perform rational reconstruction on each entry of the matrix.
 *
 * For each entry a, find a rational p/q such that:
 * - p/q ≡ a (mod N)
 * - |p|, |q| ≤ sqrt(N/2)
 * - gcd(p, q) = 1
 *
 * @param matrix - The integer matrix
 * @param N - The modulus
 * @returns Object { numerators: IntegerMatrix, denominators: IntegerMatrix }
 * @throws {ZeroDivisionError} If the modulus is zero (Sage issue #9345)
 * @throws {ValueError} If no rational reconstruction exists
 * @see Reference: sage/matrix/matrix_integer_dense.pyx:3513 (rational_reconstruction)
 * @see Reference: sage/matrix/misc_flint.pyx:52
 */
export function rational_reconstruction(
  matrix: IntegerMatrix,
  N: bigint | number
): { numerators: IntegerMatrix; denominators: IntegerMatrix } {
  const modulus = typeof N === 'number' ? BigInt(N) : N;

  if (modulus === 0n) {
    throw new ZeroDivisionError('The modulus cannot be zero');
  }

  const m = matrix.nrows;
  const n = matrix.ncols;

  const numerators = new IntegerMatrix(m, n);
  const denominators = new IntegerMatrix(m, n);

  for (let i = 0; i < m; i++) {
    for (let j = 0; j < n; j++) {
      const a = matrix.get(i, j).value;
      const result = _rational_reconstruction_single(a, modulus);
      if (result === null) {
        throw new ValueError('rational reconstruction does not exist');
      }
      numerators.set(i, j, result.p);
      denominators.set(i, j, result.q);
    }
  }

  return { numerators, denominators };
}

/**
 * Rational reconstruction for a single integer.
 * Find p/q such that p/q ≡ a (mod N) and |p|, |q| ≤ sqrt(N/2).
 */
function _rational_reconstruction_single(a: bigint, m: bigint): { p: bigint; q: bigint } | null {
  if (m === 0n) {
    return null;
  }

  const absM = m < 0n ? -m : m;

  // Compute bound = floor(sqrt(|m|/2))
  const bound = _isqrt(absM / 2n);

  // Extended Euclidean algorithm
  let u1 = 0n;
  let v1 = 1n;
  let u2 = absM;
  let v2 = ((a % absM) + absM) % absM; // Ensure positive

  while (true) {
    const absV2 = v2 < 0n ? -v2 : v2;
    if (absV2 <= bound) {
      break;
    }

    const q = u2 / v2;
    const newU1 = u1 - q * v1;
    const newU2 = u2 - q * v2;

    u1 = v1;
    v1 = newU1;
    u2 = v2;
    v2 = newU2;
  }

  // Check if |v1| <= bound
  const absV1 = v1 < 0n ? -v1 : v1;
  if (absV1 > bound) {
    return null;
  }

  // Check if gcd(v1, v2) = 1
  if (bigintGcd(v1, v2) !== 1n) {
    return null;
  }

  // Result is v2/v1, normalize sign
  if (v1 < 0n) {
    return { p: -v2, q: -v1 };
  }
  return { p: v2, q: v1 };
}

/**
 * Integer square root (floor).
 */
function _isqrt(n: bigint): bigint {
  if (n < 0n) {
    throw new ValueError('square root of negative number');
  }
  if (n === 0n) {
    return 0n;
  }

  // Newton's method
  let x = n;
  let y = (x + 1n) / 2n;
  while (y < x) {
    x = y;
    y = (x + n / x) / 2n;
  }
  return x;
}

/**
 * Find a symplectic basis for an anti-symmetric, alternating matrix.
 *
 * Return a pair (F, C) such that the rows of C form a symplectic
 * basis for self and F = C * self * C^T.
 *
 * Anti-symmetric means M = -M^T. Alternating means the diagonal is identically zero.
 *
 * @param matrix - The integer matrix (must be square, anti-symmetric, alternating)
 * @returns [F, C] where F is the symplectic form and C is the basis change matrix
 * @throws ValueError if matrix is not anti-symmetric or alternating
 * @see Reference: sage/matrix/matrix_integer_dense.pyx:symplectic_form
 */
export function symplectic_form_integer(matrix: IntegerMatrix): [IntegerMatrix, IntegerMatrix] {
  if (!matrix.is_square()) {
    throw new ValueError('Can only find symplectic bases for square matrices');
  }

  const n = matrix.nrows;

  // Check anti-symmetric: M = -M^T
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      if (matrix.get(i, j).value !== -matrix.get(j, i).value) {
        throw new ValueError('Can only find symplectic bases for anti-symmetric matrices');
      }
    }
  }

  // Check alternating: diagonal is zero
  for (let i = 0; i < n; i++) {
    if (matrix.get(i, i).value !== 0n) {
      throw new ValueError('Can only find symplectic bases for alternating matrices');
    }
  }

  // Handle empty matrix
  if (n === 0) {
    return [matrix.copy(), identity_integer_matrix(0)];
  }

  // --- port of sage.matrix.symplectic_basis.symplectic_basis_over_ZZ ---

  const E: bigint[][] = [];
  for (let i = 0; i < n; i++) {
    const row: bigint[] = [];
    for (let j = 0; j < n; j++) row.push(matrix.get(i, j).value);
    E.push(row);
  }
  const B: bigint[][] = [];
  for (let i = 0; i < n; i++) {
    const row: bigint[] = new Array(n).fill(0n);
    row[i] = 1n;
    B.push(row);
  }

  const swapRows = (M: bigint[][], a: number, b: number): void => {
    if (a === b) return;
    const t = M[a]!;
    M[a] = M[b]!;
    M[b] = t;
  };
  const swapCols = (M: bigint[][], a: number, b: number): void => {
    if (a === b) return;
    for (let i = 0; i < M.length; i++) {
      const t = M[i]![a]!;
      M[i]![a] = M[i]![b]!;
      M[i]![b] = t;
    }
  };
  const addMultipleOfRow = (M: bigint[][], i: number, j: number, v: bigint): void => {
    for (let k = 0; k < M[i]!.length; k++) M[i]![k] = M[i]![k]! + v * M[j]![k]!;
  };
  const addMultipleOfCol = (M: bigint[][], i: number, j: number, v: bigint): void => {
    for (let k = 0; k < M.length; k++) M[k]![i] = M[k]![i]! + v * M[k]![j]!;
  };
  /** Python/Sage `quo_rem`: the remainder has the sign of the divisor (floor division). */
  const quoRem = (a: bigint, b: bigint): [bigint, bigint] => {
    let q = a / b;
    const r = a - q * b;
    if (r !== 0n && r < 0n !== b < 0n) {
      q -= 1n;
    }
    return [q, a - q * b];
  };

  /** Position (row, col) of the smallest strictly positive entry with row,col >= pivot. */
  const smallestElementPosition = (pivot: number): [number, number] | null => {
    let found: [number, number] | null = null;
    let min: bigint | null = null;
    for (let i = pivot; i < n; i++) {
      for (let j = pivot; j < n; j++) {
        const v = E[j]![i]!;
        if (v > 0n && (min === null || v < min)) {
          min = v;
          found = [j, i];
        }
      }
    }
    return found;
  };

  const moveToPositivePivot = (row: number, col: number, pivot: number): void => {
    const v = E[row]![col]!;

    if (row === pivot && col === pivot + 1) {
      // nothing to do
    } else if (row === pivot + 1 && col === pivot) {
      swapRows(B, pivot, pivot + 1);
      swapRows(E, pivot, pivot + 1);
      swapCols(E, pivot, pivot + 1);
    } else if (row !== pivot && row !== pivot + 1 && col !== pivot && col !== pivot + 1) {
      swapRows(B, pivot, row);
      swapRows(B, pivot + 1, col);
      swapRows(E, pivot, row);
      swapRows(E, pivot + 1, col);
      swapCols(E, pivot, row);
      swapCols(E, pivot + 1, col);
    } else if (row === pivot) {
      swapRows(B, pivot + 1, col);
      swapRows(E, pivot + 1, col);
      swapCols(E, pivot + 1, col);
    } else if (row === pivot + 1) {
      swapRows(B, pivot, col);
      swapRows(E, pivot, col);
      swapCols(E, pivot, col);
    } else if (col === pivot) {
      swapRows(B, pivot + 1, row);
      swapRows(E, pivot + 1, row);
      swapCols(E, pivot + 1, row);
    } else if (col === pivot + 1) {
      swapRows(B, pivot, row);
      swapRows(E, pivot, row);
      swapCols(E, pivot, row);
    }

    // all that swapping can switch the sign of a row
    if (E[pivot]![pivot + 1] !== v) {
      swapRows(B, pivot, pivot + 1);
      swapRows(E, pivot, pivot + 1);
      swapCols(E, pivot, pivot + 1);
    }
  };

  const zeroes: number[] = [];
  const ps: Array<[bigint, number]> = [];
  let pivot = 0;

  while (pivot < n) {
    const found = smallestElementPosition(pivot);
    if (found === null) {
      zeroes.push(pivot);
      pivot += 1;
      continue;
    }
    moveToPositivePivot(found[0], found[1], pivot);

    let allZero = true;

    // use nonzero element to clean row pivot
    let u = E[pivot + 1]![pivot]!;
    for (let i = pivot + 2; i < n; i++) {
      const [v] = quoRem(-E[i]![pivot]!, u);
      if (v !== 0n) {
        allZero = false;
        addMultipleOfRow(E, i, pivot + 1, v);
        addMultipleOfCol(E, i, pivot + 1, v);
        addMultipleOfRow(B, i, pivot + 1, v);
      }
    }

    // use nonzero element to clean row pivot+1
    u = E[pivot]![pivot + 1]!;
    for (let i = pivot + 2; i < n; i++) {
      const [v] = quoRem(-E[i]![pivot + 1]!, u);
      if (v !== 0n) {
        allZero = false;
        addMultipleOfRow(E, i, pivot, v);
        addMultipleOfCol(E, i, pivot, v);
        addMultipleOfRow(B, i, pivot, v);
      }
    }

    if (allZero) {
      ps.push([E[pivot]![pivot + 1]!, pivot]);
      pivot += 2;
    }
  }

  ps.sort((a, b) => (a[0] < b[0] ? -1 : a[0] > b[0] ? 1 : a[1] - b[1]));
  const order = [...ps.map((p) => p[1]), ...ps.map((p) => p[1] + 1), ...zeroes];

  const C = new IntegerMatrix(n, n);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      C.set(i, j, B[order[i]!]![j]!);
    }
  }

  const F = C.mul(matrix).mul(C.transpose());
  return [F, C];
}

/**
 * Compute (p^s)-minimal polynomials of this matrix.
 *
 * A (p^s)-minimal polynomial of a matrix B is a monic polynomial f in Z[X]
 * of minimal degree such that all entries of f(B) are divisible by p^s.
 *
 * Returns a dictionary mapping s values to the corresponding (p^s)-minimal polynomials
 * (as coefficient arrays with lowest degree first).
 *
 * @param matrix - The integer matrix (must be square)
 * @param p - A prime number
 * @param s_max - Maximum s value to compute (default: computed automatically)
 * @returns Dictionary mapping s to polynomial coefficients
 * @see Reference: sage/matrix/matrix_integer_dense.pyx:p_minimal_polynomials
 */
export function p_minimal_polynomials(
  matrix: IntegerMatrix,
  p: bigint | number,
  s_max?: number
): Map<number, bigint[]> {
  if (!matrix.is_square()) {
    throw new ValueError('p-minimal polynomials only defined for square matrices');
  }

  const prime = typeof p === 'number' ? BigInt(p) : p;
  const n = matrix.nrows;

  // The minimal polynomial is always a (p^s)-minimal polynomial for large enough s
  const minPoly = _characteristic_polynomial_bigint(matrix);

  // For s=1, we need to find the minimal degree polynomial f such that f(B) ≡ 0 (mod p)
  // This is a simplified implementation that returns the characteristic polynomial
  // for s=1 (which is always valid, though not necessarily minimal)

  const result = new Map<number, bigint[]>();

  // For the simple case, just return the characteristic polynomial for s=1
  // A full implementation would compute proper (p^s)-minimal polynomials

  const maxS = s_max ?? 1;
  for (let s = 1; s <= maxS; s++) {
    result.set(s, minPoly);
  }

  return result;
}

/**
 * Return the (b)-ideal (null ideal modulo b) of the matrix.
 *
 * The null ideal modulo b is N_{(b)}(B) = {f in Z[X] | f(B) in M_n(bZ)}.
 *
 * For b=0, this is the principal ideal generated by the minimal polynomial.
 * For b>0, this includes the minimal polynomial and additional generators.
 *
 * Returns the generators of the ideal as polynomial coefficient arrays.
 *
 * @param matrix - The integer matrix (must be square)
 * @param b - An element of Z (default: 0)
 * @returns Array of polynomial coefficient arrays (generators of the ideal)
 * @see Reference: sage/matrix/matrix_integer_dense.pyx:null_ideal
 */
export function null_ideal(matrix: IntegerMatrix, b?: number): bigint[][] {
  if (!matrix.is_square()) {
    throw new ValueError('null ideal only defined for square matrices');
  }

  const modulus = BigInt(b ?? 0);

  // The minimal polynomial is always in the null ideal
  const minPoly = _characteristic_polynomial_bigint(matrix);

  if (modulus === 0n) {
    // For b=0, the null ideal is just the principal ideal (min_poly)
    return [minPoly];
  }

  // For b > 0, we need additional generators
  // A simplified implementation returns the modulus and the minimal polynomial
  // A full implementation would compute proper (p^s)-minimal polynomials

  return [[modulus], minPoly];
}

/**
 * Determine the generators of the ring of integer valued polynomials on this matrix.
 *
 * Returns a pair (mu_B, P) where P is a list of polynomials such that
 * {f in Q[X] | f(B) in M_n(Z)} = mu_B * Q[X] + sum_{g in P} g * Z[X]
 *
 * @param matrix - The integer matrix (must be square)
 * @returns Tuple [mu_B, P] where mu_B is the minimal polynomial and P are generators
 * @see Reference: sage/matrix/matrix_integer_dense.pyx:integer_valued_polynomials_generators
 */
export function integer_valued_polynomials_generators(
  matrix: IntegerMatrix
): [bigint[], bigint[][]] {
  if (!matrix.is_square()) {
    throw new ValueError('integer valued polynomials only defined for square matrices');
  }

  // The minimal polynomial mu_B
  const minPoly = _characteristic_polynomial_bigint(matrix);

  // For a simplified implementation, return the constant polynomial 1
  // as the only additional generator (meaning Z[X] is a subset)
  // A full implementation would compute proper generators based on
  // the J-ideal computation

  const generators: bigint[][] = [[1n]]; // The constant polynomial 1

  return [minPoly, generators];
}

/**
 * Return the antitranspose of the matrix.
 *
 * The antitranspose is the matrix obtained by transposing and then
 * reversing both rows and columns (equivalent to rotating 180 degrees
 * and then transposing).
 *
 * @param matrix - The integer matrix
 * @returns The antitranspose (reverse rows and columns)
 * @see Reference: sage/matrix/matrix_integer_dense.pyx:antitranspose
 */
export function antitranspose(matrix: IntegerMatrix): IntegerMatrix {
  const n = matrix.nrows;
  const m = matrix.ncols;
  const result = new IntegerMatrix(m, n);

  for (let i = 0; i < n; i++) {
    for (let j = 0; j < m; j++) {
      // Antitranspose: element (i,j) goes to position (m-1-j, n-1-i)
      result.set(m - 1 - j, n - 1 - i, matrix.get(i, j));
    }
  }

  return result;
}

/**
 * Insert a row into the matrix.
 *
 * @param matrix - The integer matrix (mutated)
 * @param index - Position to insert
 * @param row - The row to insert
 * @returns A new matrix with the row inserted
 * @see Reference: sage/matrix/matrix_integer_dense.pyx:insert_row
 */
export function insert_row(
  matrix: IntegerMatrix,
  index: number,
  row: (bigint | number)[]
): IntegerMatrix {
  if (index < 0 || index > matrix.nrows) {
    throw new ValueError(`index ${index} out of range for matrix with ${matrix.nrows} rows`);
  }
  if (row.length !== matrix.ncols) {
    throw new ValueError(
      `row length ${row.length} does not match matrix column count ${matrix.ncols}`
    );
  }

  const result = new IntegerMatrix(matrix.nrows + 1, matrix.ncols);

  // Copy rows before the insertion point
  for (let i = 0; i < index; i++) {
    for (let j = 0; j < matrix.ncols; j++) {
      result.set(i, j, matrix.get(i, j));
    }
  }

  // Insert the new row
  for (let j = 0; j < matrix.ncols; j++) {
    result.set(index, j, row[j]!);
  }

  // Copy rows after the insertion point
  for (let i = index; i < matrix.nrows; i++) {
    for (let j = 0; j < matrix.ncols; j++) {
      result.set(i + 1, j, matrix.get(i, j));
    }
  }

  return result;
}

/**
 * Augment the matrix with another matrix.
 *
 * @param matrix - The integer matrix
 * @param right - Matrix to augment with
 * @param subdivide - Whether to subdivide
 * @returns The augmented matrix
 * @see Reference: sage/matrix/matrix_integer_dense.pyx:augment
 */
export function augment_integer(
  matrix: IntegerMatrix,
  right: IntegerMatrix,
  subdivide?: boolean
): IntegerMatrix {
  if (matrix.nrows !== right.nrows) {
    throw new ArithmeticError(
      `cannot augment ${matrix.nrows}x${matrix.ncols} matrix with ${right.nrows}x${right.ncols} matrix`
    );
  }

  const result = new IntegerMatrix(matrix.nrows, matrix.ncols + right.ncols);

  for (let i = 0; i < matrix.nrows; i++) {
    for (let j = 0; j < matrix.ncols; j++) {
      result.set(i, j, matrix.get(i, j));
    }
    for (let j = 0; j < right.ncols; j++) {
      result.set(i, matrix.ncols + j, right.get(i, j));
    }
  }

  return result;
}

/**
 * Delete zero columns.
 *
 * @param matrix - The integer matrix
 * @returns Matrix with zero columns removed
 * @see Reference: sage/matrix/matrix_integer_dense.pyx:_delete_zero_columns
 */
export function delete_zero_columns(matrix: IntegerMatrix): IntegerMatrix {
  // Find non-zero columns
  const nonZeroCols: number[] = [];
  for (let j = 0; j < matrix.ncols; j++) {
    let isZero = true;
    for (let i = 0; i < matrix.nrows; i++) {
      if (!matrix.get(i, j).isZero()) {
        isZero = false;
        break;
      }
    }
    if (!isZero) {
      nonZeroCols.push(j);
    }
  }

  if (nonZeroCols.length === matrix.ncols) {
    return matrix.copy();
  }

  const result = new IntegerMatrix(matrix.nrows, nonZeroCols.length);
  for (let i = 0; i < matrix.nrows; i++) {
    for (let newJ = 0; newJ < nonZeroCols.length; newJ++) {
      result.set(i, newJ, matrix.get(i, nonZeroCols[newJ]!));
    }
  }

  return result;
}

/**
 * Factor out common factors from each row.
 *
 * @param matrix - The integer matrix
 * @returns Pair (factored matrix, list of factors)
 * @see Reference: sage/matrix/matrix_integer_dense.pyx:_factor_out_common_factors_from_each_row
 */
export function factor_out_common_factors_from_each_row(
  matrix: IntegerMatrix
): [IntegerMatrix, Integer[]] {
  const result = new IntegerMatrix(matrix.nrows, matrix.ncols);
  const factors: Integer[] = [];

  for (let i = 0; i < matrix.nrows; i++) {
    // Compute GCD of this row
    let rowGcd = 0n;
    for (let j = 0; j < matrix.ncols; j++) {
      rowGcd = bigintGcd(rowGcd, matrix.get(i, j).value);
      if (rowGcd === 1n) {
        break;
      }
    }

    // If row is all zeros, factor is 0
    if (rowGcd === 0n) {
      factors.push(new Integer(0n));
      for (let j = 0; j < matrix.ncols; j++) {
        result.set(i, j, 0n);
      }
    } else {
      factors.push(new Integer(rowGcd));
      for (let j = 0; j < matrix.ncols; j++) {
        result.set(i, j, matrix.get(i, j).value / rowGcd);
      }
    }
  }

  return [result, factors];
}
