/**
 * @module sage/matrix/matrix_modn
 * @description Z/nZ matrix operations - important for cryptography
 *
 * Port of: sage/matrix/matrix_modn_dense_double.pyx, sage/matrix/matrix_modn_dense_float.pyx
 */

import { ArithmeticError, NotImplementedError, ValueError, ZeroDivisionError } from '../errors.js';
import { gcd, xgcd, inverse_mod } from '../arith/misc.js';

/**
 * Compute a mod n, ensuring the result is in [0, n).
 */
function mod(a: bigint, n: bigint): bigint {
  const result = a % n;
  return result < 0n ? result + n : result;
}

/**
 * A dense matrix over Z/nZ where n fits in machine arithmetic.
 *
 * This is important for cryptographic applications involving modular arithmetic.
 *
 * @see Reference: sage/matrix/matrix_modn_dense_double.pyx:Matrix_modn_dense_double
 */
export class Matrix_modn_dense {
  readonly nrows: number;
  readonly ncols: number;
  readonly modulus: bigint;
  private _entries: bigint[][];

  /**
   * Create a matrix over Z/nZ.
   *
   * @param nrows - Number of rows
   * @param ncols - Number of columns
   * @param modulus - The modulus n
   * @param entries - Optional 2D array of entries
   */
  constructor(
    nrows: number,
    ncols: number,
    modulus: bigint,
    entries?: (bigint | number)[][]
  ) {
    if (modulus <= 0n) {
      throw new ValueError('modulus must be positive');
    }

    this.nrows = nrows;
    this.ncols = ncols;
    this.modulus = modulus;

    // Initialize entries array
    this._entries = [];
    for (let i = 0; i < nrows; i++) {
      this._entries.push([]);
      for (let j = 0; j < ncols; j++) {
        if (entries !== undefined && entries[i] !== undefined && entries[i][j] !== undefined) {
          const val = typeof entries[i][j] === 'number' ? BigInt(entries[i][j]!) : entries[i][j]!;
          this._entries[i]!.push(mod(val as bigint, modulus));
        } else {
          this._entries[i]!.push(0n);
        }
      }
    }
  }

  /**
   * Get entry at (i, j).
   *
   * @param i - Row index
   * @param j - Column index
   * @returns The entry modulo n
   */
  get(i: number, j: number): bigint {
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
   * @param value - The new value
   */
  set(i: number, j: number, value: bigint | number): void {
    if (i < 0 || i >= this.nrows || j < 0 || j >= this.ncols) {
      throw new ValueError(`index out of bounds: (${i}, ${j}) for ${this.nrows}x${this.ncols} matrix`);
    }
    const val = typeof value === 'number' ? BigInt(value) : value;
    this._entries[i]![j] = mod(val, this.modulus);
  }

  /**
   * Matrix addition.
   *
   * @param other - Matrix to add
   * @returns Sum of matrices
   */
  add(other: Matrix_modn_dense): Matrix_modn_dense {
    if (this.nrows !== other.nrows || this.ncols !== other.ncols) {
      throw new ArithmeticError(
        `cannot add ${this.nrows}x${this.ncols} matrix to ${other.nrows}x${other.ncols} matrix`
      );
    }
    if (this.modulus !== other.modulus) {
      throw new ArithmeticError('matrices must have the same modulus');
    }

    const result = new Matrix_modn_dense(this.nrows, this.ncols, this.modulus);
    for (let i = 0; i < this.nrows; i++) {
      for (let j = 0; j < this.ncols; j++) {
        result._entries[i]![j] = mod(this._entries[i]![j]! + other._entries[i]![j]!, this.modulus);
      }
    }
    return result;
  }

  /**
   * Matrix subtraction.
   *
   * @param other - Matrix to subtract
   * @returns Difference of matrices
   */
  sub(other: Matrix_modn_dense): Matrix_modn_dense {
    if (this.nrows !== other.nrows || this.ncols !== other.ncols) {
      throw new ArithmeticError(
        `cannot subtract ${other.nrows}x${other.ncols} matrix from ${this.nrows}x${this.ncols} matrix`
      );
    }
    if (this.modulus !== other.modulus) {
      throw new ArithmeticError('matrices must have the same modulus');
    }

    const result = new Matrix_modn_dense(this.nrows, this.ncols, this.modulus);
    for (let i = 0; i < this.nrows; i++) {
      for (let j = 0; j < this.ncols; j++) {
        result._entries[i]![j] = mod(this._entries[i]![j]! - other._entries[i]![j]!, this.modulus);
      }
    }
    return result;
  }

  /**
   * Matrix negation.
   *
   * @returns Negation of matrix
   */
  neg(): Matrix_modn_dense {
    const result = new Matrix_modn_dense(this.nrows, this.ncols, this.modulus);
    for (let i = 0; i < this.nrows; i++) {
      for (let j = 0; j < this.ncols; j++) {
        const val = this._entries[i]![j]!;
        result._entries[i]![j] = val === 0n ? 0n : this.modulus - val;
      }
    }
    return result;
  }

  /**
   * Matrix multiplication.
   *
   * @param other - Matrix to multiply
   * @returns Product of matrices
   */
  mul(other: Matrix_modn_dense): Matrix_modn_dense {
    if (this.ncols !== other.nrows) {
      throw new ArithmeticError(
        `cannot multiply ${this.nrows}x${this.ncols} matrix by ${other.nrows}x${other.ncols} matrix`
      );
    }
    if (this.modulus !== other.modulus) {
      throw new ArithmeticError('matrices must have the same modulus');
    }

    const result = new Matrix_modn_dense(this.nrows, other.ncols, this.modulus);
    for (let i = 0; i < this.nrows; i++) {
      for (let j = 0; j < other.ncols; j++) {
        let sum = 0n;
        for (let k = 0; k < this.ncols; k++) {
          sum = mod(sum + this._entries[i]![k]! * other._entries[k]![j]!, this.modulus);
        }
        result._entries[i]![j] = sum;
      }
    }
    return result;
  }

  /**
   * Scalar multiplication.
   *
   * @param scalar - Scalar to multiply by
   * @returns Scaled matrix
   */
  scalar_mul(scalar: bigint | number): Matrix_modn_dense {
    const s = mod(typeof scalar === 'number' ? BigInt(scalar) : scalar, this.modulus);
    const result = new Matrix_modn_dense(this.nrows, this.ncols, this.modulus);
    for (let i = 0; i < this.nrows; i++) {
      for (let j = 0; j < this.ncols; j++) {
        result._entries[i]![j] = mod(s * this._entries[i]![j]!, this.modulus);
      }
    }
    return result;
  }

  /**
   * Return the transpose.
   *
   * @returns The transpose
   */
  transpose(): Matrix_modn_dense {
    const result = new Matrix_modn_dense(this.ncols, this.nrows, this.modulus);
    for (let i = 0; i < this.nrows; i++) {
      for (let j = 0; j < this.ncols; j++) {
        result._entries[j]![i] = this._entries[i]![j]!;
      }
    }
    return result;
  }

  /**
   * Return the determinant.
   *
   * Uses LU decomposition. Only works when the modulus is prime.
   *
   * @returns The determinant modulo n
   */
  determinant(): bigint {
    if (this.nrows !== this.ncols) {
      throw new ArithmeticError('determinant is only defined for square matrices');
    }

    const n = this.nrows;

    if (n === 0) {
      return 1n;
    }

    if (n === 1) {
      return this._entries[0]![0]!;
    }

    if (n === 2) {
      const a = this._entries[0]![0]!;
      const b = this._entries[0]![1]!;
      const c = this._entries[1]![0]!;
      const d = this._entries[1]![1]!;
      return mod(a * d - b * c, this.modulus);
    }

    // Gaussian elimination with partial pivoting
    const M: bigint[][] = [];
    for (let i = 0; i < n; i++) {
      M.push([...this._entries[i]!]);
    }

    let det = 1n;
    let sign = 1n;

    for (let k = 0; k < n; k++) {
      // Find pivot
      let pivotRow = -1;
      for (let i = k; i < n; i++) {
        if (M[i]![k] !== 0n) {
          pivotRow = i;
          break;
        }
      }

      if (pivotRow === -1) {
        return 0n; // Matrix is singular
      }

      // Swap rows if needed
      if (pivotRow !== k) {
        [M[k], M[pivotRow]] = [M[pivotRow]!, M[k]!];
        sign = -sign;
      }

      const pivot = M[k]![k]!;
      det = mod(det * pivot, this.modulus);

      // Find inverse of pivot (requires modulus to be prime for general case)
      const [g, s] = xgcd(pivot, this.modulus);
      if (g !== 1n) {
        // If pivot is not invertible, we need a different approach
        // For now, throw an error for non-prime moduli where pivot is not invertible
        throw new ArithmeticError('determinant computation requires prime modulus when pivot is not invertible');
      }
      const pivotInv = mod(s, this.modulus);

      // Eliminate column k below diagonal
      for (let i = k + 1; i < n; i++) {
        if (M[i]![k] !== 0n) {
          const factor = mod(M[i]![k]! * pivotInv, this.modulus);
          for (let j = k; j < n; j++) {
            M[i]![j] = mod(M[i]![j]! - factor * M[k]![j]!, this.modulus);
          }
        }
      }
    }

    return mod(sign * det, this.modulus);
  }

  /**
   * Return the rank.
   *
   * @returns The rank
   */
  rank(): number {
    // Compute row echelon form and count non-zero rows
    const echelon = this.echelon_form();
    let r = 0;
    for (let i = 0; i < echelon.nrows; i++) {
      let isZeroRow = true;
      for (let j = 0; j < echelon.ncols; j++) {
        if (echelon._entries[i]![j] !== 0n) {
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
   * Put matrix in echelon form (in place).
   */
  echelonize(): void {
    const n = this.nrows;
    const m = this.ncols;

    let pivotRow = 0;

    for (let col = 0; col < m && pivotRow < n; col++) {
      // Find pivot in this column
      let found = -1;
      for (let i = pivotRow; i < n; i++) {
        if (this._entries[i]![col] !== 0n) {
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

      // Scale pivot row to have leading 1
      const pivot = this._entries[pivotRow]![col]!;
      const [g, s] = xgcd(pivot, this.modulus);
      if (g !== 1n) {
        // Pivot is not invertible - skip to next column for non-prime moduli
        pivotRow++;
        continue;
      }
      const pivotInv = mod(s, this.modulus);

      for (let j = col; j < m; j++) {
        this._entries[pivotRow]![j] = mod(this._entries[pivotRow]![j]! * pivotInv, this.modulus);
      }

      // Eliminate other rows
      for (let i = 0; i < n; i++) {
        if (i !== pivotRow && this._entries[i]![col] !== 0n) {
          const factor = this._entries[i]![col]!;
          for (let j = col; j < m; j++) {
            this._entries[i]![j] = mod(this._entries[i]![j]! - factor * this._entries[pivotRow]![j]!, this.modulus);
          }
        }
      }

      pivotRow++;
    }
  }

  /**
   * Return the echelon form.
   *
   * @returns The echelon form
   */
  echelon_form(): Matrix_modn_dense {
    const result = this.copy();
    result.echelonize();
    return result;
  }

  /**
   * Return the inverse.
   *
   * @returns The inverse
   */
  inverse(): Matrix_modn_dense {
    if (this.nrows !== this.ncols) {
      throw new ArithmeticError('inverse is only defined for square matrices');
    }

    const n = this.nrows;

    if (n === 0) {
      return this.copy();
    }

    // Augment with identity matrix
    const augmented = new Matrix_modn_dense(n, 2 * n, this.modulus);
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        augmented._entries[i]![j] = this._entries[i]![j]!;
      }
      augmented._entries[i]![n + i] = 1n;
    }

    // Gaussian elimination
    for (let k = 0; k < n; k++) {
      // Find pivot
      let pivotRow = -1;
      for (let i = k; i < n; i++) {
        if (augmented._entries[i]![k] !== 0n) {
          const [g] = xgcd(augmented._entries[i]![k]!, this.modulus);
          if (g === 1n) {
            pivotRow = i;
            break;
          }
        }
      }

      if (pivotRow === -1) {
        throw new ZeroDivisionError('matrix is not invertible');
      }

      // Swap rows
      if (pivotRow !== k) {
        [augmented._entries[k], augmented._entries[pivotRow]] = [augmented._entries[pivotRow]!, augmented._entries[k]!];
      }

      // Scale pivot row
      const pivot = augmented._entries[k]![k]!;
      const [, s] = xgcd(pivot, this.modulus);
      const pivotInv = mod(s, this.modulus);

      for (let j = 0; j < 2 * n; j++) {
        augmented._entries[k]![j] = mod(augmented._entries[k]![j]! * pivotInv, this.modulus);
      }

      // Eliminate other rows
      for (let i = 0; i < n; i++) {
        if (i !== k && augmented._entries[i]![k] !== 0n) {
          const factor = augmented._entries[i]![k]!;
          for (let j = 0; j < 2 * n; j++) {
            augmented._entries[i]![j] = mod(
              augmented._entries[i]![j]! - factor * augmented._entries[k]![j]!,
              this.modulus
            );
          }
        }
      }
    }

    // Extract inverse
    const result = new Matrix_modn_dense(n, n, this.modulus);
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        result._entries[i]![j] = augmented._entries[i]![n + j]!;
      }
    }

    return result;
  }

  /**
   * Return the characteristic polynomial.
   *
   * @param variable - Variable name
   * @returns The characteristic polynomial coefficients (from constant term to leading)
   */
  charpoly(variable?: string): bigint[] {
    if (this.nrows !== this.ncols) {
      throw new ArithmeticError('characteristic polynomial is only defined for square matrices');
    }

    const n = this.nrows;

    if (n === 0) {
      return [1n]; // det(xI - A) = 1 for 0x0 matrix
    }

    // Use Faddeev-LeVerrier algorithm
    // The characteristic polynomial is det(xI - A) = x^n - c_{n-1} x^{n-1} - ... - c_0
    // where the c_i are computed iteratively

    const coeffs: bigint[] = new Array(n + 1).fill(0n);
    coeffs[n] = 1n; // Leading coefficient

    let B = this.copy();
    let trace = 0n;

    for (let i = 0; i < n; i++) {
      trace = 0n;
      for (let j = 0; j < n; j++) {
        trace = mod(trace + B._entries[j]![j]!, this.modulus);
      }

      // c_{n-1-i} = trace / (i + 1)
      const [g, s] = xgcd(BigInt(i + 1), this.modulus);
      if (g !== 1n) {
        throw new ArithmeticError('characteristic polynomial computation requires prime modulus');
      }
      const divInv = mod(s, this.modulus);
      coeffs[n - 1 - i] = mod(-trace * divInv, this.modulus);

      if (i < n - 1) {
        // B = A * (B - c_{n-1-i} * I)
        const nextB = new Matrix_modn_dense(n, n, this.modulus);

        // First compute B - c * I
        const c = mod(-coeffs[n - 1 - i]!, this.modulus);
        for (let r = 0; r < n; r++) {
          for (let cc = 0; cc < n; cc++) {
            if (r === cc) {
              nextB._entries[r]![cc] = mod(B._entries[r]![cc]! + c, this.modulus);
            } else {
              nextB._entries[r]![cc] = B._entries[r]![cc]!;
            }
          }
        }

        // Then multiply by A
        B = this.mul(nextB);
      }
    }

    return coeffs;
  }

  /**
   * Return the minimal polynomial.
   *
   * Uses the Krylov sequence approach: compute v, Av, A^2v, ... until linear
   * dependence is found. The minimal polynomial is the polynomial of least
   * degree that annihilates the matrix.
   *
   * @param variable - Variable name
   * @returns The minimal polynomial coefficients (from constant term to leading)
   * @see Reference: sage/matrix/matrix2.pyx:minpoly
   */
  minpoly(variable?: string): bigint[] {
    if (this.nrows !== this.ncols) {
      throw new ArithmeticError('minimal polynomial is only defined for square matrices');
    }

    const n = this.nrows;

    if (n === 0) {
      return [1n]; // Minimal polynomial of 0x0 matrix is 1
    }

    // Get characteristic polynomial first
    const charpoly = this.charpoly(variable);

    // Helper: compute LCM of two monic polynomials over Z/pZ
    // We use the formula: LCM(f, g) = f * g / GCD(f, g)
    const polyGCD = (f: bigint[], g: bigint[]): bigint[] => {
      // Euclidean algorithm for polynomials
      while (g.length > 0 && !(g.length === 1 && g[0] === 0n)) {
        const [_, r] = polyDivMod(f, g);
        f = g;
        g = r;
      }
      // Make monic
      if (f.length === 0) return [1n];
      const lead = f[f.length - 1]!;
      if (lead === 0n) return [1n];
      const leadInv = inverse_mod(lead, this.modulus);
      return f.map(c => mod(c * leadInv, this.modulus));
    };

    const polyDivMod = (
      dividend: bigint[],
      divisor: bigint[]
    ): [bigint[], bigint[]] => {
      if (divisor.length === 0 || (divisor.length === 1 && divisor[0] === 0n)) {
        throw new ArithmeticError('division by zero polynomial');
      }

      const result: bigint[] = [];
      let remainder = [...dividend];

      const divisorLead = divisor[divisor.length - 1]!;
      const divisorLeadInv = inverse_mod(divisorLead, this.modulus);

      while (remainder.length >= divisor.length) {
        const coeff = mod(remainder[remainder.length - 1]! * divisorLeadInv, this.modulus);
        const degDiff = remainder.length - divisor.length;

        result.unshift(coeff);

        for (let i = 0; i < divisor.length; i++) {
          remainder[degDiff + i] = mod(
            remainder[degDiff + i]! - coeff * divisor[i]!,
            this.modulus
          );
        }

        // Remove leading zeros
        while (remainder.length > 0 && remainder[remainder.length - 1] === 0n) {
          remainder.pop();
        }
      }

      if (result.length === 0) {
        result.push(0n);
      }

      if (remainder.length === 0) {
        remainder.push(0n);
      }

      return [result, remainder];
    };

    const polyMul = (f: bigint[], g: bigint[]): bigint[] => {
      if (f.length === 0 || g.length === 0) return [0n];
      const result: bigint[] = new Array(f.length + g.length - 1).fill(0n);
      for (let i = 0; i < f.length; i++) {
        for (let j = 0; j < g.length; j++) {
          result[i + j] = mod(result[i + j]! + f[i]! * g[j]!, this.modulus);
        }
      }
      // Remove leading zeros
      while (result.length > 1 && result[result.length - 1] === 0n) {
        result.pop();
      }
      return result;
    };

    const polyLCM = (f: bigint[], g: bigint[]): bigint[] => {
      const gcd = polyGCD(f, g);
      const [fDivGcd] = polyDivMod(f, gcd);
      return polyMul(fDivGcd, g);
    };

    // Compute minimal polynomial using Krylov sequences with multiple starting vectors
    // The minimal polynomial is the LCM of the minimal polynomials for each standard basis vector
    let minPoly: bigint[] = [1n]; // Start with 1

    for (let startIdx = 0; startIdx < n; startIdx++) {
      // Use standard basis vector e_startIdx
      let v: bigint[] = new Array(n).fill(0n);
      v[startIdx] = 1n;

      // Compute the Krylov sequence: v, Av, A^2v, ...
      const krylov: bigint[][] = [v.slice()];
      for (let k = 0; k < n; k++) {
        // Multiply by A: Av
        const Av: bigint[] = new Array(n).fill(0n);
        for (let i = 0; i < n; i++) {
          let sum = 0n;
          for (let j = 0; j < n; j++) {
            sum = mod(sum + this._entries[i]![j]! * v[j]!, this.modulus);
          }
          Av[i] = sum;
        }
        krylov.push(Av.slice());
        v = Av;
      }

      // Find the first linear dependence by row-reducing the Krylov matrix
      // The columns are the Krylov vectors
      const krylovMatrix = new Matrix_modn_dense(n, krylov.length, this.modulus);
      for (let j = 0; j < krylov.length; j++) {
        for (let i = 0; i < n; i++) {
          krylovMatrix._entries[i]![j] = krylov[j]![i]!;
        }
      }

      // Compute RREF to find pivot columns
      const echelon = krylovMatrix.echelon_form();

      // Find pivot columns
      const pivotCols: number[] = [];
      let pivotRow = 0;
      for (let col = 0; col < krylov.length && pivotRow < n; col++) {
        if (echelon._entries[pivotRow]![col] !== 0n) {
          pivotCols.push(col);
          pivotRow++;
        }
      }

      // Find the first non-pivot column - this gives the minimal degree for this vector
      let localDegree = n + 1;
      for (let j = 0; j < krylov.length; j++) {
        if (!pivotCols.includes(j)) {
          localDegree = j;
          break;
        }
      }

      if (localDegree <= n) {
        // Extract the coefficients from the RREF
        // Column localDegree is a linear combination of previous columns
        const localCoeffs: bigint[] = new Array(localDegree + 1).fill(0n);
        for (let i = 0; i < localDegree; i++) {
          // Find the row where column i is a pivot
          for (let r = 0; r < pivotCols.length; r++) {
            if (pivotCols[r] === i) {
              localCoeffs[i] = mod(-echelon._entries[r]![localDegree]!, this.modulus);
              break;
            }
          }
        }
        localCoeffs[localDegree] = 1n; // Leading coefficient is 1

        // Compute LCM with current minimal polynomial
        minPoly = polyLCM(minPoly, localCoeffs);
      }

      // Early exit: if minPoly has degree n, it equals charpoly
      if (minPoly.length - 1 >= n) {
        return charpoly;
      }
    }

    return minPoly;
  }

  /**
   * Solve a linear system A * X = B.
   *
   * @param B - Right-hand side
   * @param check - Whether to verify the solution
   * @returns The solution X
   */
  solve_right(B: Matrix_modn_dense, check?: boolean): Matrix_modn_dense {
    if (this.nrows !== this.ncols) {
      throw new ArithmeticError('solve_right requires a square matrix');
    }
    if (this.nrows !== B.nrows) {
      throw new ArithmeticError('incompatible matrix dimensions');
    }
    if (this.modulus !== B.modulus) {
      throw new ArithmeticError('matrices must have the same modulus');
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
   */
  right_kernel_matrix(): Matrix_modn_dense {
    const m = this.nrows;
    const n = this.ncols;

    // Compute the reduced row echelon form of self^T
    // The kernel of A is the same as the cokernel of A^T
    const At = this.transpose();
    const echelon = At.echelon_form();

    // Find pivot columns
    const pivotCols: number[] = [];
    let pivotRow = 0;
    for (let j = 0; j < At.ncols && pivotRow < At.nrows; j++) {
      if (echelon._entries[pivotRow]![j] !== 0n) {
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
      return new Matrix_modn_dense(0, n, this.modulus);
    }

    // Build kernel basis
    const kernel = new Matrix_modn_dense(kernelDim, n, this.modulus);

    for (let i = 0; i < kernelDim; i++) {
      const freeCol = freeCols[i]!;
      kernel._entries[i]![freeCol] = 1n;

      // Fill in the pivot columns
      for (let r = 0; r < pivotCols.length; r++) {
        const pivotCol = pivotCols[r]!;
        // Find the value that makes the row sum to zero
        kernel._entries[i]![pivotCol] = mod(-echelon._entries[r]![freeCol]!, this.modulus);
      }
    }

    return kernel;
  }

  /**
   * Return a copy.
   *
   * @returns A copy
   */
  copy(): Matrix_modn_dense {
    const result = new Matrix_modn_dense(this.nrows, this.ncols, this.modulus);
    for (let i = 0; i < this.nrows; i++) {
      for (let j = 0; j < this.ncols; j++) {
        result._entries[i]![j] = this._entries[i]![j]!;
      }
    }
    return result;
  }

  /**
   * Check equality.
   *
   * @param other - Matrix to compare
   * @returns True if equal
   */
  eq(other: Matrix_modn_dense): boolean {
    if (this.nrows !== other.nrows || this.ncols !== other.ncols) {
      return false;
    }
    if (this.modulus !== other.modulus) {
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
    if (this.nrows === 0 || this.ncols === 0) {
      return `${this.nrows} x ${this.ncols} empty matrix over Z/${this.modulus}Z`;
    }

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

    const lines: string[] = [];
    for (let i = 0; i < this.nrows; i++) {
      const row = strings[i]!.map((s, j) => s.padStart(widths[j]!));
      lines.push('[' + row.join(' ') + ']');
    }

    return lines.join('\n');
  }
}

// ============================================================================
// Factory functions
// ============================================================================

/**
 * Create a zero matrix over Z/nZ.
 *
 * @param modulus - The modulus
 * @param nrows - Number of rows
 * @param ncols - Number of columns (default: nrows)
 * @returns Zero matrix
 */
export function zero_matrix_modn(modulus: bigint, nrows: number, ncols?: number): Matrix_modn_dense {
  return new Matrix_modn_dense(nrows, ncols ?? nrows, modulus);
}

/**
 * Create an identity matrix over Z/nZ.
 *
 * @param modulus - The modulus
 * @param n - Size
 * @returns Identity matrix
 */
export function identity_matrix_modn(modulus: bigint, n: number): Matrix_modn_dense {
  const result = new Matrix_modn_dense(n, n, modulus);
  for (let i = 0; i < n; i++) {
    result.set(i, i, 1n);
  }
  return result;
}

/**
 * Create a random matrix over Z/nZ.
 *
 * @param modulus - The modulus
 * @param nrows - Number of rows
 * @param ncols - Number of columns (default: nrows)
 * @returns Random matrix
 */
export function random_matrix_modn(
  modulus: bigint,
  nrows: number,
  ncols?: number
): Matrix_modn_dense {
  const nc = ncols ?? nrows;
  const result = new Matrix_modn_dense(nrows, nc, modulus);

  for (let i = 0; i < nrows; i++) {
    for (let j = 0; j < nc; j++) {
      // Generate a random value in [0, modulus)
      // Using Math.random() for simplicity; in production, use a cryptographically secure RNG
      const randVal = BigInt(Math.floor(Math.random() * Number(modulus)));
      result.set(i, j, randVal);
    }
  }

  return result;
}

/**
 * Create a matrix over Z/nZ from entries.
 *
 * @param modulus - The modulus
 * @param entries - 2D array of entries
 * @returns The matrix
 */
export function matrix_modn_from_entries(
  modulus: bigint,
  entries: (bigint | number)[][]
): Matrix_modn_dense {
  if (entries.length === 0) {
    return new Matrix_modn_dense(0, 0, modulus);
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

  return new Matrix_modn_dense(nrows, ncols, modulus, entries);
}
