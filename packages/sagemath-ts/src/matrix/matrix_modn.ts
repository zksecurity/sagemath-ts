/**
 * @module sage/matrix/matrix_modn
 * @description Z/nZ matrix operations - important for cryptography
 *
 * Port of: sage/matrix/matrix_modn_dense_double.pyx, sage/matrix/matrix_modn_dense_float.pyx
 */

import { matkermod } from '@sagemath-ts/parigp-ts';
import { inverse_mod, is_prime, xgcd } from '../arith/misc.js';
import { ArithmeticError, NotImplementedError, ValueError, ZeroDivisionError } from '../errors.js';

/**
 * Compute a mod n, ensuring the result is in [0, n).
 */
function mod(a: bigint, n: bigint): bigint {
  const result = a % n;
  return result < 0n ? result + n : result;
}

/**
 * Lift a residue in [0, n) to the symmetric representative in (-n/2, n/2].
 *
 * Mirrors `IntegerMod.lift_centered`, which Sage's generic determinant uses
 * before computing over ZZ (`matrix2.pyx:2394`).
 */
function liftCentered(a: bigint, n: bigint): bigint {
  return 2n * a > n ? a - n : a;
}

/**
 * Exact determinant over ZZ of a small integer matrix using the
 * fraction-free Bareiss algorithm (no floating point, no division that is
 * not exact).
 */
function determinantZZ(M: bigint[][]): bigint {
  const n = M.length;
  if (n === 0) return 1n;
  const A = M.map((row) => [...row]);
  let sign = 1n;
  let prev = 1n;

  for (let k = 0; k < n - 1; k++) {
    if (A[k]![k] === 0n) {
      let swapRow = -1;
      for (let i = k + 1; i < n; i++) {
        if (A[i]![k] !== 0n) {
          swapRow = i;
          break;
        }
      }
      if (swapRow === -1) {
        return 0n;
      }
      [A[k], A[swapRow]] = [A[swapRow]!, A[k]!];
      sign = -sign;
    }
    for (let i = k + 1; i < n; i++) {
      for (let j = k + 1; j < n; j++) {
        A[i]![j] = (A[i]![j]! * A[k]![k]! - A[i]![k]! * A[k]![j]!) / prev;
      }
      A[i]![k] = 0n;
    }
    prev = A[k]![k]!;
  }

  return sign * A[n - 1]![n - 1]!;
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
  constructor(nrows: number, ncols: number, modulus: bigint, entries?: (bigint | number)[][]) {
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
      throw new ValueError(
        `index out of bounds: (${i}, ${j}) for ${this.nrows}x${this.ncols} matrix`
      );
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
      throw new ValueError(
        `index out of bounds: (${i}, ${j}) for ${this.nrows}x${this.ncols} matrix`
      );
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
   * Follows Sage: for `n <= 3` the naive formula is used, otherwise the
   * matrix is lifted to ZZ with centered representatives and the exact
   * integer determinant is reduced modulo `n`.  This works for *every*
   * modulus, prime or composite -- no division by a possibly non-invertible
   * pivot is ever performed.
   *
   * @returns The determinant modulo n
   * @see Reference: sage/matrix/matrix_modn_dense_template.pxi:2406 (determinant)
   * @see Reference: sage/matrix/matrix2.pyx:2366-2398 (generic determinant, Z/nZ branch)
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

    if (n === 3) {
      const e = this._entries;
      const d =
        e[0]![0]! * (e[1]![1]! * e[2]![2]! - e[1]![2]! * e[2]![1]!) -
        e[1]![0]! * (e[0]![1]! * e[2]![2]! - e[0]![2]! * e[2]![1]!) +
        e[2]![0]! * (e[0]![1]! * e[1]![2]! - e[0]![2]! * e[1]![1]!);
      return mod(d, this.modulus);
    }

    // Lift to ZZ (centered representatives) and compute the exact determinant.
    const M: bigint[][] = [];
    for (let i = 0; i < n; i++) {
      M.push(this._entries[i]!.map((x) => liftCentered(x, this.modulus)));
    }

    return mod(determinantZZ(M), this.modulus);
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
   *
   * Sage only implements the echelon form over Z/nZ when n is prime (the ring
   * is then a field); for composite moduli it raises `NotImplementedError`.
   *
   * @see Reference: sage/matrix/matrix_modn_dense_template.pxi:1632 (echelonize)
   */
  echelonize(): void {
    if (!is_prime(this.modulus)) {
      throw new NotImplementedError(
        `Echelon form not implemented over 'Ring of integers modulo ${this.modulus}'.`
      );
    }

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
        [this._entries[pivotRow], this._entries[found]] = [
          this._entries[found]!,
          this._entries[pivotRow]!,
        ];
      }

      // Scale pivot row to have leading 1.  The modulus is prime here, so a
      // nonzero pivot is always invertible.
      const pivot = this._entries[pivotRow]![col]!;
      const pivotInv = inverse_mod(pivot, this.modulus);

      for (let j = col; j < m; j++) {
        this._entries[pivotRow]![j] = mod(this._entries[pivotRow]![j]! * pivotInv, this.modulus);
      }

      // Eliminate other rows
      for (let i = 0; i < n; i++) {
        if (i !== pivotRow && this._entries[i]![col] !== 0n) {
          const factor = this._entries[i]![col]!;
          for (let j = col; j < m; j++) {
            this._entries[i]![j] = mod(
              this._entries[i]![j]! - factor * this._entries[pivotRow]![j]!,
              this.modulus
            );
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
   * Return the pivot columns of the (reduced row) echelon form of this matrix.
   *
   * @returns The indices of the pivot columns, in increasing order
   * @see Reference: sage/matrix/matrix_modn_dense_template.pxi:2019 (pivots)
   */
  pivots(): number[] {
    const echelon = this.echelon_form();
    const result: number[] = [];
    let row = 0;
    for (let j = 0; j < echelon.ncols && row < echelon.nrows; j++) {
      if (echelon._entries[row]![j] !== 0n) {
        result.push(j);
        row++;
      }
    }
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
        [augmented._entries[k], augmented._entries[pivotRow]] = [
          augmented._entries[pivotRow]!,
          augmented._entries[k]!,
        ];
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
   * Return the characteristic polynomial `det(x*I - self)`.
   *
   * Uses the division-free algorithm of Seifullin that Sage falls back on
   * whenever LinBox is unavailable (p = 2 or a composite modulus).  Being
   * division free it is valid over *every* Z/nZ; in particular it does not
   * need to divide by `i + 1`, which is what made the previous
   * Faddeev-LeVerrier implementation fail for small primes.
   *
   * @param variable - Variable name (accepted for signature compatibility)
   * @returns The characteristic polynomial coefficients (from constant term to leading)
   * @see Reference: sage/matrix/matrix_modn_dense_template.pxi:1443 (charpoly)
   * @see Reference: sage/matrix/matrix2.pyx:3342 (_charpoly_df)
   */
  charpoly(variable?: string): bigint[] {
    if (this.nrows !== this.ncols) {
      throw new ArithmeticError('characteristic polynomial is only defined for square matrices');
    }

    const n = this.nrows;
    const p = this.modulus;

    if (n === 0) {
      return [1n]; // det(xI - A) = 1 for 0x0 matrix
    }

    const M = this._entries;

    // In the notation of Algorithm 3.1 of Seifullin (as ported in
    // matrix2.pyx:_charpoly_df):
    //   F[p] is the coefficient of x^{n-p-1} of the characteristic polynomial,
    //   a[p] is a vector of length n, A[p] a scalar.
    const F: bigint[] = new Array(n).fill(0n);
    const A: bigint[] = new Array(n).fill(0n);
    const a: bigint[][] = [];
    for (let i = 0; i < Math.max(n - 1, 1); i++) {
      a.push(new Array(n).fill(0n));
    }

    F[0] = mod(-M[0]![0]!, p);

    for (let t = 1; t < n; t++) {
      // a(0, t) := M(<=t, t)
      for (let i = 0; i <= t; i++) {
        a[0]![i] = M[i]![t]!;
      }

      A[0] = M[t]![t]!;

      for (let q = 1; q < t; q++) {
        // a(q, t) := M[<=t, <=t] * a(q-1, t)
        for (let i = 0; i <= t; i++) {
          let s = 0n;
          for (let j = 0; j <= t; j++) {
            s = mod(s + M[i]![j]! * a[q - 1]![j]!, p);
          }
          a[q]![i] = s;
        }
        A[q] = a[q]![t]!;
      }

      // A[t] := M[t, <=t] * a(t-1, t)
      let s = 0n;
      for (let j = 0; j <= t; j++) {
        s = mod(s + M[t]![j]! * a[t - 1]![j]!, p);
      }
      A[t] = s;

      for (let q = 0; q <= t; q++) {
        let acc = F[q]!;
        for (let k = 0; k < q; k++) {
          acc = mod(acc - A[k]! * F[q - k - 1]!, p);
        }
        F[q] = mod(acc - A[q]!, p);
      }
    }

    // f = x^n + sum_{p} F[p] * x^{n-p-1}; return constant term first.
    const coeffs: bigint[] = new Array(n + 1).fill(0n);
    coeffs[n] = 1n;
    for (let i = 0; i < n; i++) {
      coeffs[n - 1 - i] = F[i]!;
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
      return f.map((c) => mod(c * leadInv, this.modulus));
    };

    const polyDivMod = (dividend: bigint[], divisor: bigint[]): [bigint[], bigint[]] => {
      if (divisor.length === 0 || (divisor.length === 1 && divisor[0] === 0n)) {
        throw new ArithmeticError('division by zero polynomial');
      }

      const result: bigint[] = [];
      const remainder = [...dividend];

      const divisorLead = divisor[divisor.length - 1]!;
      const divisorLeadInv = inverse_mod(divisorLead, this.modulus);

      while (remainder.length >= divisor.length) {
        const coeff = mod(remainder[remainder.length - 1]! * divisorLeadInv, this.modulus);
        const degDiff = remainder.length - divisor.length;

        result.unshift(coeff);

        for (let i = 0; i < divisor.length; i++) {
          remainder[degDiff + i] = mod(remainder[degDiff + i]! - coeff * divisor[i]!, this.modulus);
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
   * Return a matrix whose rows form a basis for the right kernel of `self`,
   * i.e. a matrix `X` with `self * X.transpose() == 0`.
   *
   * The construction follows Sage exactly: `self` is echelonized, the
   * non-pivot columns index the basis vectors, and the resulting matrix is
   * echelonized again (the `'echelon'` default).
   *
   * If the modulus is composite the keyword arguments are ignored and the
   * computation is delegated to PARI's `matkermod`, exactly as Sage does
   * (`matrix_modn_dense_template.pxi:2136` falls back to
   * `Matrix_dense.right_kernel_matrix`, which reaches
   * `matrix2.pyx:4233 _right_kernel_matrix_over_integer_mod_ring`).
   *
   * @param options - `basis` is one of `'echelon'` (default), `'pivot'` or `'computed'`
   * @returns Kernel matrix
   * @see Reference: sage/matrix/matrix_modn_dense_template.pxi:2072 (right_kernel_matrix)
   * @see Reference: sage/matrix/matrix2.pyx:4233 (_right_kernel_matrix_over_integer_mod_ring)
   */
  right_kernel_matrix(options?: {
    basis?: 'echelon' | 'pivot' | 'computed';
  }): Matrix_modn_dense {
    if (!is_prime(this.modulus)) {
      // Composite modulus: Sage's echelon_form raises NotImplementedError and the
      // template falls back to Matrix_dense.right_kernel_matrix(self) — called with
      // *no* keyword arguments, so `basis` defaults to 'computed' and the result is
      // returned exactly as PARI's matkermod produced it (matrix2.pyx:4931).
      return this._right_kernel_matrix_over_integer_mod_ring();
    }

    const basis = options?.basis ?? 'echelon';
    if (basis !== 'echelon' && basis !== 'pivot' && basis !== 'computed') {
      throw new ValueError('matrix kernel basis format not recognized');
    }

    const ncols = this.ncols;

    // Echelonize self (NOT the transpose) and read off its pivots.
    const echelon = this.echelon_form();
    const pivots: number[] = [];
    let row = 0;
    for (let j = 0; j < ncols && row < echelon.nrows; j++) {
      if (echelon._entries[row]![j] !== 0n) {
        pivots.push(j);
        row++;
      }
    }
    const r = pivots.length;

    const nonpivots: number[] = [];
    for (let j = 0; j < ncols; j++) {
      if (!pivots.includes(j)) {
        nonpivots.push(j);
      }
    }

    const M = new Matrix_modn_dense(ncols - r, ncols, this.modulus);
    // 'computed' basis (as returned by Sage): free coordinate is -1 and the
    // pivot coordinates are the echelon entries themselves.
    const pm1 = this.modulus - 1n;
    for (let i = 0; i < ncols - r; i++) {
      M._entries[i]![nonpivots[i]!] = pm1;
      for (let j = 0; j < r; j++) {
        M._entries[i]![pivots[j]!] = echelon._entries[j]![nonpivots[i]!]!;
      }
    }

    if (basis === 'computed') {
      return M;
    }
    const P = M.neg();
    if (basis === 'pivot') {
      return P;
    }
    P.echelonize();
    return P;
  }

  /**
   * Return a matrix whose rows are a basis for the right kernel of `self` over
   * `Z/nZ`, as computed by PARI's `matkermod`.
   *
   * This is Sage's `_right_kernel_matrix_over_integer_mod_ring`, which lifts the
   * matrix to `ZZ`, calls `matkermod(n)` and reads the returned PARI *columns* as
   * the rows of the result.  The basis is returned unchanged ('computed' format):
   * a composite `Z/nZ` is not a field, so Sage's default basis format is
   * `'computed'` and no echelonization is performed.
   *
   * Deviation: the image is requested from `matkermod` even though it is
   * discarded.  PARI's `matkermod` (`bb_hnf.c:1049`) shortcuts a tall matrix via
   * `shallowtrans(matimagemod(shallowtrans(A), d))` when no image is wanted and
   * `m > 2n`; a `t_MAT` with zero columns does not record its row count, so for
   * `A == 0 (mod d)` — whose image is empty — the transpose collapses to `0 x 0`
   * and PARI reports an *empty* kernel instead of the whole module.  Asking for
   * the image disables that branch, so the kernel is always the true one.
   * Verified against SageMath 10.3 over 2400 small cases: the two agree
   * everywhere except exactly the zero matrices with `m > 2n`, where SageMath
   * returns a basis that does not span the kernel.
   *
   * @returns Matrix `X` with `self * X.transpose() == 0` over `Z/nZ`
   * @see Reference: sage/matrix/matrix2.pyx:4233 (_right_kernel_matrix_over_integer_mod_ring)
   * @see Reference: pari/src/basemath/bb_hnf.c:1036 (matkermod)
   */
  _right_kernel_matrix_over_integer_mod_ring(): Matrix_modn_dense {
    // PARI's t_MAT is column-major; build it directly so that a matrix with zero
    // rows still carries its column count (`zm_from_rows` cannot know it).
    const A: bigint[][] = [];
    for (let j = 0; j < this.ncols; j++) {
      const col: bigint[] = new Array(this.nrows);
      for (let i = 0; i < this.nrows; i++) col[i] = this._entries[i]![j]!;
      A.push(col);
    }

    const { ker } = matkermod(A, this.modulus, /* wantIm = */ true);

    // Each PARI column is one basis vector; they become the rows of the result.
    const M = new Matrix_modn_dense(ker.length, this.ncols, this.modulus);
    for (let i = 0; i < ker.length; i++) {
      const v = ker[i]!;
      for (let j = 0; j < this.ncols; j++) {
        M._entries[i]![j] = mod(v[j]!, this.modulus);
      }
    }
    return M;
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
export function zero_matrix_modn(
  modulus: bigint,
  nrows: number,
  ncols?: number
): Matrix_modn_dense {
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
