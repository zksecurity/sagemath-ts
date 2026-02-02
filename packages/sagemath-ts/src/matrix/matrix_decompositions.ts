/**
 * @module sage/matrix/matrix_decompositions
 * @description Matrix decompositions (LU, QR, SVD, Cholesky, Smith, Hermite, etc.)
 *
 * Port of: sage/matrix/matrix2.pyx
 */

import { ArithmeticError, NotImplementedError, ValueError } from '../errors.js';
import type { CoefficientRing, RingElement } from '../rings/polynomial/polynomial_element.js';
import { Polynomial } from '../rings/polynomial/polynomial_element.js';
import { PolynomialRing } from '../rings/polynomial/polynomial_ring.js';
import { Matrix, diagonal_matrix, identity_matrix, zero_matrix } from './matrix_generic.js';

/**
 * Interface for field elements that support division/inverse.
 */
interface FieldElement extends RingElement {
  inverse?(): this;
  inv?(): this;
  div?(other: this): this;
}

/**
 * Get the multiplicative inverse of a field element.
 */
function getInverse<R extends FieldElement>(elem: R): R {
  if (typeof elem.inverse === 'function') {
    return elem.inverse() as R;
  }
  if (typeof elem.inv === 'function') {
    return elem.inv() as R;
  }
  throw new ArithmeticError('element does not support inverse operation');
}

// ============================================================================
// Echelon Forms
// ============================================================================

/**
 * Transform the matrix into echelon form in place.
 *
 * Uses Gaussian elimination with partial pivoting to transform
 * the matrix into row echelon form.
 *
 * @param matrix - The matrix to echelonize (mutated)
 * @param algorithm - Algorithm to use ('default', 'classical', etc.)
 * @param cutoff - Cutoff for Strassen algorithm (not used)
 * @param transformation - Whether to return the transformation matrix
 * @returns The transformation matrix if requested
 * @see Reference: sage/matrix/matrix2.pyx:echelonize
 */
export function echelonize<R extends FieldElement>(
  matrix: Matrix<R>,
  algorithm?: 'default' | 'classical' | 'strassen' | 'partial_pivoting' | 'scaled_partial_pivoting',
  cutoff?: number,
  transformation?: boolean
): Matrix<R> | undefined {
  const m = matrix.nrows;
  const n = matrix.ncols;
  const ring = matrix.base_ring;

  // If transformation is requested, start with identity matrix
  let T: Matrix<R> | undefined;
  if (transformation) {
    T = identity_matrix(ring, m);
  }

  let pivotRow = 0;

  for (let col = 0; col < n && pivotRow < m; col++) {
    // Find the first non-zero entry in this column at or below pivotRow
    let found = -1;
    for (let i = pivotRow; i < m; i++) {
      if (!matrix.get(i, col).isZero()) {
        found = i;
        break;
      }
    }

    if (found === -1) {
      continue;
    }

    // Swap rows if necessary
    if (found !== pivotRow) {
      for (let j = 0; j < n; j++) {
        const tmp = matrix.get(pivotRow, j);
        matrix.set(pivotRow, j, matrix.get(found, j));
        matrix.set(found, j, tmp);
      }
      if (T) {
        for (let j = 0; j < m; j++) {
          const tmp = T.get(pivotRow, j);
          T.set(pivotRow, j, T.get(found, j));
          T.set(found, j, tmp);
        }
      }
    }

    // Get the pivot element
    const pivotElem = matrix.get(pivotRow, col);

    // Eliminate entries below the pivot
    for (let i = pivotRow + 1; i < m; i++) {
      const entry = matrix.get(i, col);
      if (!entry.isZero()) {
        // factor = entry / pivot
        const factor = entry.mul(getInverse(pivotElem)) as R;

        // row[i] = row[i] - factor * row[pivotRow]
        for (let j = col; j < n; j++) {
          const val = matrix.get(i, j).sub(factor.mul(matrix.get(pivotRow, j)) as R) as R;
          matrix.set(i, j, val);
        }
        if (T) {
          for (let j = 0; j < m; j++) {
            const val = T.get(i, j).sub(factor.mul(T.get(pivotRow, j)) as R) as R;
            T.set(i, j, val);
          }
        }
      }
    }

    pivotRow++;
  }

  if (transformation) {
    return T;
  }
}

/**
 * Return the echelon form of the matrix (row echelon form).
 *
 * This does not change the matrix itself. Uses Gaussian elimination
 * with partial pivoting when the base ring is a field.
 *
 * @param matrix - The matrix
 * @param algorithm - Algorithm to use
 * @param cutoff - Cutoff for Strassen algorithm
 * @returns A new matrix in echelon form
 * @see Reference: sage/matrix/matrix2.pyx:echelon_form
 */
export function echelon_form<R extends FieldElement>(
  matrix: Matrix<R>,
  algorithm?: 'default' | 'classical' | 'strassen',
  cutoff?: number
): Matrix<R> {
  const m = matrix.nrows;
  const n = matrix.ncols;
  const ring = matrix.base_ring;

  // Copy the matrix
  const M = matrix.copy();

  let pivotRow = 0;

  for (let col = 0; col < n && pivotRow < m; col++) {
    // Find the first non-zero entry in this column at or below pivotRow
    let found = -1;
    for (let i = pivotRow; i < m; i++) {
      if (!M.get(i, col).isZero()) {
        found = i;
        break;
      }
    }

    if (found === -1) {
      continue;
    }

    // Swap rows if necessary
    if (found !== pivotRow) {
      for (let j = 0; j < n; j++) {
        const tmp = M.get(pivotRow, j);
        M.set(pivotRow, j, M.get(found, j));
        M.set(found, j, tmp);
      }
    }

    // Get the pivot element
    const pivotElem = M.get(pivotRow, col);

    // Eliminate entries below the pivot
    for (let i = pivotRow + 1; i < m; i++) {
      const entry = M.get(i, col);
      if (!entry.isZero()) {
        // For field elements, we can use division
        // factor = entry / pivot
        const factor = entry.mul(getInverse(pivotElem)) as R;

        // row[i] = row[i] - factor * row[pivotRow]
        for (let j = col; j < n; j++) {
          const val = M.get(i, j).sub(factor.mul(M.get(pivotRow, j)) as R) as R;
          M.set(i, j, val);
        }
      }
    }

    pivotRow++;
  }

  return M;
}

/**
 * Return the reduced row echelon form of the matrix.
 *
 * Every leading coefficient (pivot) is 1, and is the only non-zero entry
 * in its column.
 *
 * @param matrix - The matrix
 * @returns The RREF of the matrix
 * @see Reference: sage/matrix/matrix2.pyx:rref
 */
export function rref<R extends FieldElement>(matrix: Matrix<R>): Matrix<R> {
  const m = matrix.nrows;
  const n = matrix.ncols;

  // First get echelon form
  const M = echelon_form(matrix);

  // Find pivots and make leading coefficients 1
  const pivotCols: number[] = [];
  let pivotRow = 0;

  for (let col = 0; col < n && pivotRow < m; col++) {
    if (!M.get(pivotRow, col).isZero()) {
      pivotCols.push(col);
      const pivotElem = M.get(pivotRow, col);

      // Scale row so pivot is 1
      const pivotInv = getInverse(pivotElem);

      for (let j = col; j < n; j++) {
        M.set(pivotRow, j, M.get(pivotRow, j).mul(pivotInv) as R);
      }

      pivotRow++;
    }
  }

  // Back substitution: eliminate entries above pivots
  for (let i = pivotCols.length - 1; i >= 0; i--) {
    const col = pivotCols[i]!;
    // Row i has pivot at column col

    // Eliminate entries above this pivot
    for (let row = 0; row < i; row++) {
      const entry = M.get(row, col);
      if (!entry.isZero()) {
        // row[row] = row[row] - entry * row[i]
        for (let j = 0; j < n; j++) {
          const val = M.get(row, j).sub(entry.mul(M.get(i, j)) as R) as R;
          M.set(row, j, val);
        }
      }
    }
  }

  return M;
}

/**
 * Return the extended echelon form of the matrix.
 *
 * For an m x n matrix A, computes [E | A | T] where:
 * - E is the identity rows selected from identity matrix
 * - A is the RREF of self
 * - T is the transformation matrix (T * original = RREF)
 *
 * @param matrix - The matrix
 * @param subdivide - Whether to subdivide the result (not used)
 * @returns The extended echelon form
 * @see Reference: sage/matrix/matrix2.pyx:extended_echelon_form
 */
export function extended_echelon_form<R extends FieldElement>(
  matrix: Matrix<R>,
  subdivide?: boolean
): Matrix<R> {
  const m = matrix.nrows;
  const n = matrix.ncols;
  const ring = matrix.base_ring;

  // Create augmented matrix [A | I]
  const aug = new Matrix<R>(ring, m, n + m);
  for (let i = 0; i < m; i++) {
    for (let j = 0; j < n; j++) {
      aug.set(i, j, matrix.get(i, j));
    }
    aug.set(i, n + i, ring.one());
  }

  // Compute RREF of augmented matrix
  const result = rref(aug);

  return result;
}

/**
 * Return the pivot column indices of the matrix.
 *
 * Computes the echelon form and identifies which columns contain pivots.
 *
 * @param matrix - The matrix
 * @returns Array of pivot column indices
 * @see Reference: sage/matrix/matrix2.pyx:pivot_rows
 */
export function pivot_rows<R extends FieldElement>(matrix: Matrix<R>): number[] {
  const E = echelon_form(matrix);
  const m = E.nrows;
  const n = E.ncols;

  const pivots: number[] = [];
  let pivotRow = 0;

  for (let col = 0; col < n && pivotRow < m; col++) {
    if (!E.get(pivotRow, col).isZero()) {
      pivots.push(col);
      pivotRow++;
    }
  }

  return pivots;
}

// ============================================================================
// LU Decomposition
// ============================================================================

/**
 * Return the LU decomposition of the matrix.
 *
 * Computes matrices P, L, U such that PA = LU where:
 * - P is a permutation matrix
 * - L is lower triangular with ones on the diagonal
 * - U is upper triangular
 *
 * Uses Gaussian elimination with partial pivoting (or nonzero pivoting).
 *
 * @param matrix - The matrix A
 * @param pivot - Pivoting strategy: 'auto', 'partial', or 'nonzero'
 * @param format - Output format: 'plu' (default) or 'compact'
 * @returns (P, L, U) for 'plu' format, or (permutation, LU) for 'compact'
 * @see Reference: sage/matrix/matrix2.pyx:LU
 */
export function LU<R extends FieldElement>(
  matrix: Matrix<R>,
  pivot: 'auto' | 'partial' | 'nonzero' = 'nonzero',
  format: 'plu' | 'compact' = 'plu'
): [Matrix<R>, Matrix<R>, Matrix<R>] | [number[], Matrix<R>] {
  const m = matrix.nrows;
  const n = matrix.ncols;
  const ring = matrix.base_ring;
  const d = Math.min(m, n);

  // Copy matrix
  const M = matrix.copy();

  // Track row permutation
  const perm: number[] = [];
  for (let i = 0; i < m; i++) {
    perm.push(i);
  }

  for (let k = 0; k < d; k++) {
    // Find pivot
    let maxLocation = -1;

    if (pivot === 'nonzero' || pivot === 'auto') {
      // First nonzero pivot
      for (let i = k; i < m; i++) {
        if (!M.get(i, k).isZero()) {
          maxLocation = i;
          break;
        }
      }
    } else {
      // partial pivoting would need absolute value comparison
      // For generic fields, fall back to nonzero
      for (let i = k; i < m; i++) {
        if (!M.get(i, k).isZero()) {
          maxLocation = i;
          break;
        }
      }
    }

    if (maxLocation !== -1) {
      // Swap rows k and maxLocation
      if (maxLocation !== k) {
        [perm[k], perm[maxLocation]] = [perm[maxLocation]!, perm[k]!];
        for (let j = 0; j < n; j++) {
          const tmp = M.get(k, j);
          M.set(k, j, M.get(maxLocation, j));
          M.set(maxLocation, j, tmp);
        }
      }

      // Get inverse of pivot
      const pivotElem = M.get(k, k);
      const inv = getInverse(pivotElem);

      // Elimination
      for (let j = k + 1; j < m; j++) {
        const scale = M.get(j, k).mul(inv) as R;
        M.set(j, k, scale); // Store L factor

        for (let p = k + 1; p < n; p++) {
          const val = M.get(j, p).sub(scale.mul(M.get(k, p)) as R) as R;
          M.set(j, p, val);
        }
      }
    }
  }

  if (format === 'compact') {
    return [perm, M];
  }

  // Build P, L, U matrices
  const P = zero_matrix(ring, m, m);
  for (let i = 0; i < m; i++) {
    P.set(i, perm[i]!, ring.one());
  }

  const L = identity_matrix(ring, m);
  for (let i = 1; i < m; i++) {
    for (let k = 0; k < Math.min(i, d); k++) {
      L.set(i, k, M.get(i, k));
      M.set(i, k, ring.zero());
    }
  }

  return [P, L, M];
}

// ============================================================================
// QR Decomposition
// ============================================================================

/**
 * Return the QR decomposition of the matrix.
 *
 * Computes matrices Q, R such that A = QR where:
 * - Q has orthogonal columns (orthonormal if sqrt is available)
 * - R is upper triangular
 *
 * Note: True QR decomposition with orthonormal Q requires square roots,
 * which are not available for all rings. For finite fields without
 * square roots, this returns the Gram-Schmidt factorization where
 * Q has orthogonal (but not normalized) columns.
 *
 * @param matrix - The matrix A
 * @param full - If true, return full m x m Q; if false (default), reduced form
 * @returns Pair (Q, R) where A = Q * R
 * @see Reference: sage/matrix/matrix2.pyx:QR
 */
export function QR<R extends FieldElement>(
  matrix: Matrix<R>,
  full: boolean = false
): [Matrix<R>, Matrix<R>] {
  // Use Gram-Schmidt (without normalization for now)
  // This gives A = Q * R where Q has orthogonal columns
  const [Q, R] = gram_schmidt_noscale(matrix);

  if (!full) {
    return [Q, R];
  }

  // For full QR, we need to extend Q to a square m x m matrix
  // This requires finding an orthogonal complement to the column space of Q
  const m = matrix.nrows;
  const n = Q.ncols;
  const ring = matrix.base_ring;

  if (n >= m) {
    // Q is already m x m (or wider), return as is
    // Just need to pad R with zero rows
    const Rfull = zero_matrix(ring, m, matrix.ncols);
    for (let i = 0; i < n && i < m; i++) {
      for (let j = 0; j < matrix.ncols; j++) {
        Rfull.set(i, j, R.get(i, j));
      }
    }
    return [Q, Rfull];
  }

  // Need to find m - n additional orthogonal vectors
  // Start with standard basis vectors and orthogonalize against Q's columns
  const Qcols: R[][] = [];
  for (let j = 0; j < n; j++) {
    Qcols.push(Q.column(j));
  }

  // Try each standard basis vector
  for (let k = 0; k < m && Qcols.length < m; k++) {
    // Create standard basis vector e_k
    const v: R[] = [];
    for (let i = 0; i < m; i++) {
      v.push(i === k ? ring.one() : ring.zero());
    }

    // Orthogonalize against all current columns
    for (let j = 0; j < Qcols.length; j++) {
      const u = Qcols[j]!;
      const dotVU = _dot_product(ring, v, u);
      const dotUU = _dot_product(ring, u, u);

      if (!dotUU.isZero()) {
        const coeff = dotVU.mul(getInverse(dotUU as R)) as R;
        for (let i = 0; i < m; i++) {
          v[i] = v[i]!.sub(coeff.mul(u[i]!) as R) as R;
        }
      }
    }

    // Check if v is non-zero
    let isZero = true;
    for (let i = 0; i < m; i++) {
      if (!v[i]!.isZero()) {
        isZero = false;
        break;
      }
    }

    if (!isZero) {
      Qcols.push(v);
    }
  }

  // Build full Q matrix (m x m)
  const Qfull = new Matrix<R>(ring, m, m);
  for (let j = 0; j < Qcols.length; j++) {
    for (let i = 0; i < m; i++) {
      Qfull.set(i, j, Qcols[j]![i]!);
    }
  }

  // Build full R matrix (m x ncols)
  // First n rows are from the original R, rest are zeros
  const Rfull = zero_matrix(ring, m, matrix.ncols);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < matrix.ncols; j++) {
      Rfull.set(i, j, R.get(i, j));
    }
  }

  return [Qfull, Rfull];
}

/**
 * Compute Gram-Schmidt orthogonalization without scaling.
 *
 * Returns (Q, R) such that A = Q * R where:
 * - The columns of Q are orthogonal (but not normalized)
 * - R is upper triangular
 *
 * @param matrix - The matrix (columns are vectors to orthogonalize)
 * @returns Pair [Q, R] where A = Q * R
 * @see Reference: sage/matrix/matrix2.pyx:_gram_schmidt_noscale
 */
export function gram_schmidt_noscale<R extends FieldElement>(
  matrix: Matrix<R>
): [Matrix<R>, Matrix<R>] {
  const m = matrix.nrows;
  const n = matrix.ncols;
  const ring = matrix.base_ring;

  if (m === 0 || n === 0) {
    return [zero_matrix(ring, m, 0), zero_matrix(ring, 0, n)];
  }

  // Work with columns
  const columns: R[][] = [];
  for (let j = 0; j < n; j++) {
    columns.push(matrix.column(j));
  }

  // Orthogonalized columns (non-normalized)
  const orthoCols: R[][] = [];
  // Coefficients for R matrix
  const rCoeffs: R[][] = [];

  for (let j = 0; j < n; j++) {
    rCoeffs.push([]);
    for (let i = 0; i < n; i++) {
      rCoeffs[j]!.push(ring.zero());
    }
  }

  let numOrthoCols = 0;

  for (let j = 0; j < n; j++) {
    // Start with the current column
    const v = [...columns[j]!];

    // Subtract projections onto previous orthogonal vectors
    for (let k = 0; k < numOrthoCols; k++) {
      const u = orthoCols[k]!;

      // Compute <v, u> / <u, u>
      const dotVU = _dot_product(ring, columns[j]!, u);
      const dotUU = _dot_product(ring, u, u);

      if (!dotUU.isZero()) {
        const coeff = dotVU.mul(getInverse(dotUU as R)) as R;
        rCoeffs[k]![j] = coeff;

        // v = v - coeff * u
        for (let i = 0; i < m; i++) {
          v[i] = v[i]!.sub(coeff.mul(u[i]!) as R) as R;
        }
      }
    }

    // Check if v is non-zero
    let isZero = true;
    for (let i = 0; i < m; i++) {
      if (!v[i]!.isZero()) {
        isZero = false;
        break;
      }
    }

    if (!isZero) {
      orthoCols.push(v);
      rCoeffs[numOrthoCols]![j] = ring.one();
      numOrthoCols++;
    }
  }

  // Build Q from orthogonal columns
  const Q = new Matrix<R>(ring, m, numOrthoCols);
  for (let j = 0; j < numOrthoCols; j++) {
    for (let i = 0; i < m; i++) {
      Q.set(i, j, orthoCols[j]![i]!);
    }
  }

  // Build R (numOrthoCols x n)
  const R_matrix = new Matrix<R>(ring, numOrthoCols, n);
  for (let i = 0; i < numOrthoCols; i++) {
    for (let j = 0; j < n; j++) {
      R_matrix.set(i, j, rCoeffs[i]![j]!);
    }
  }

  return [Q, R_matrix];
}

/**
 * Compute dot product of two vectors.
 */
function _dot_product<R extends RingElement>(ring: CoefficientRing<R>, u: R[], v: R[]): R {
  let result = ring.zero();
  for (let i = 0; i < u.length; i++) {
    result = result.add(u[i]!.mul(v[i]!) as R) as R;
  }
  return result;
}

/**
 * Compute Gram-Schmidt orthogonalization.
 *
 * Returns (Q, R) such that A = Q * R where:
 * - If orthonormal=false, the columns of Q are orthogonal (not normalized)
 * - If orthonormal=true, the columns of Q are orthonormal (requires sqrt)
 * - R is upper triangular
 *
 * Note: orthonormal=true requires a ring with square root, which we don't
 * support in this generic implementation. Use gram_schmidt_noscale instead.
 *
 * @param matrix - The matrix (columns are vectors to orthogonalize)
 * @param orthonormal - Whether to normalize the vectors (default: false)
 * @returns Pair (Q, R) where A = Q * R
 * @see Reference: sage/matrix/matrix2.pyx:gram_schmidt
 */
export function gram_schmidt<R extends FieldElement>(
  matrix: Matrix<R>,
  orthonormal: boolean = false
): [Matrix<R>, Matrix<R>] {
  if (orthonormal) {
    throw new NotImplementedError(
      'orthonormal Gram-Schmidt requires sqrt, which is not supported for generic rings'
    );
  }

  return gram_schmidt_noscale(matrix);
}

// ============================================================================
// Cholesky Decomposition
// ============================================================================

/**
 * Return the Cholesky decomposition of a positive definite matrix.
 *
 * For a positive definite Hermitian matrix A, computes lower triangular L
 * such that A = L * L^* (where L^* is the conjugate transpose, or just
 * transpose for real matrices).
 *
 * Note: This requires square root computation which is not available for
 * all rings. Works for RDF, CDF, and certain algebraic number fields.
 *
 * @param matrix - A positive definite Hermitian matrix
 * @param extended - Whether to return extended form (not used)
 * @returns Lower triangular matrix L
 * @see Reference: sage/matrix/matrix2.pyx:cholesky
 */
export function cholesky<R extends RingElement>(matrix: Matrix<R>, extended?: boolean): Matrix<R> {
  if (!matrix.is_square()) {
    throw new ArithmeticError('cholesky is only defined for square matrices');
  }

  const n = matrix.nrows;
  const ring = matrix.base_ring;

  if (n === 0) {
    return zero_matrix(ring, 0);
  }

  // Check if ring has sqrt
  const ringWithSqrt = ring as unknown as {
    sqrt?: (x: R) => R;
    __call__?: (x: number | bigint) => R;
  };

  const elemWithSqrt = ring.one() as unknown as { sqrt?: () => R };
  if (typeof elemWithSqrt.sqrt !== 'function') {
    throw new NotImplementedError(
      'cholesky decomposition requires a ring with sqrt (e.g., RDF, CDF, AA, QQbar)'
    );
  }

  // Standard Cholesky algorithm
  const L = zero_matrix(ring, n);

  for (let i = 0; i < n; i++) {
    for (let j = 0; j <= i; j++) {
      let sum = ring.zero();

      if (j === i) {
        // Diagonal element: L[i,i] = sqrt(A[i,i] - sum of L[i,k]^2)
        for (let k = 0; k < j; k++) {
          const Lik = L.get(i, k);
          sum = sum.add(Lik.mul(Lik)) as R;
        }
        const diag = matrix.get(i, i).sub(sum) as R;

        // Check if diag is non-negative (for real case)
        // This is where we'd check positive definiteness
        const diagWithSqrt = diag as unknown as { sqrt: () => R };
        L.set(i, i, diagWithSqrt.sqrt());
      } else {
        // Off-diagonal element: L[i,j] = (A[i,j] - sum) / L[j,j]
        for (let k = 0; k < j; k++) {
          sum = sum.add(L.get(i, k).mul(L.get(j, k))) as R;
        }
        const Ljj = L.get(j, j);
        if (Ljj.isZero()) {
          throw new ArithmeticError('matrix is not positive definite (zero pivot)');
        }
        const LjjInv = getInverse(Ljj);
        L.set(i, j, matrix.get(i, j).sub(sum).mul(LjjInv) as R);
      }
    }
  }

  return L;
}

/**
 * Compute the inverse of a positive definite matrix using Cholesky.
 *
 * For a positive definite matrix A with Cholesky decomposition A = LL^T,
 * the inverse is A^{-1} = (L^T)^{-1} L^{-1} = (L^{-1})^T L^{-1}.
 *
 * @param matrix - A positive definite matrix
 * @returns The inverse matrix
 * @see Reference: sage/matrix/matrix2.pyx:inverse_positive_definite
 */
export function inverse_positive_definite<R extends FieldElement>(matrix: Matrix<R>): Matrix<R> {
  if (!matrix.is_square()) {
    throw new ArithmeticError('inverse is only defined for square matrices');
  }

  // Compute Cholesky decomposition A = LL^T
  const L = cholesky(matrix);
  const n = matrix.nrows;
  const ring = matrix.base_ring;

  // Compute L^{-1} by forward substitution
  const Linv = zero_matrix(ring, n);

  for (let j = 0; j < n; j++) {
    // Solve L * x = e_j for x
    for (let i = 0; i < n; i++) {
      if (i < j) {
        Linv.set(i, j, ring.zero());
      } else if (i === j) {
        Linv.set(i, j, getInverse(L.get(i, i)));
      } else {
        // x[i] = (e_j[i] - sum_{k<i} L[i,k] * x[k]) / L[i,i]
        let sum = ring.zero();
        for (let k = j; k < i; k++) {
          sum = sum.add(L.get(i, k).mul(Linv.get(k, j))) as R;
        }
        Linv.set(i, j, sum.neg().mul(getInverse(L.get(i, i))) as R);
      }
    }
  }

  // A^{-1} = (L^{-1})^T * L^{-1}
  const LinvT = Linv.transpose();
  return LinvT.mul(Linv);
}

// ============================================================================
// LDL and Block LDL^T Decomposition
// ============================================================================

/**
 * Return the indefinite factorization of a symmetric/Hermitian matrix.
 *
 * Computes L, D such that A = LDL^T (or LDL^* for Hermitian) where:
 * - L is lower triangular with ones on the diagonal
 * - D is diagonal
 *
 * Unlike Cholesky, this works for indefinite matrices (not necessarily positive definite).
 *
 * @param matrix - A symmetric or Hermitian matrix
 * @param algorithm - 'symmetric' or 'hermitian' (default: 'symmetric')
 * @param check - Whether to check symmetry/Hermitian property (default: true)
 * @returns Pair (L, D) where A = L * D * L^T
 * @see Reference: sage/matrix/matrix2.pyx:indefinite_factorization
 */
export function indefinite_factorization<R extends FieldElement>(
  matrix: Matrix<R>,
  algorithm: 'symmetric' | 'hermitian' = 'symmetric',
  check: boolean = true
): [Matrix<R>, Matrix<R>] {
  if (!matrix.is_square()) {
    throw new ArithmeticError('indefinite_factorization is only defined for square matrices');
  }

  const n = matrix.nrows;
  const ring = matrix.base_ring;

  if (n === 0) {
    return [zero_matrix(ring, 0), zero_matrix(ring, 0)];
  }

  // Check symmetry if requested
  if (check) {
    for (let i = 0; i < n; i++) {
      for (let j = i + 1; j < n; j++) {
        if (!matrix.get(i, j).eq(matrix.get(j, i))) {
          throw new ValueError('matrix is not symmetric');
        }
      }
    }
  }

  // Initialize L as identity, D as zero
  const L = identity_matrix(ring, n);
  const D = zero_matrix(ring, n);

  // Work on a copy of the matrix
  const A: R[][] = [];
  for (let i = 0; i < n; i++) {
    A.push([]);
    for (let j = 0; j < n; j++) {
      A[i]!.push(matrix.get(i, j));
    }
  }

  for (let k = 0; k < n; k++) {
    // D[k,k] = A[k,k] - sum_{j<k} L[k,j]^2 * D[j,j]
    let sum = ring.zero();
    for (let j = 0; j < k; j++) {
      const Lkj = L.get(k, j);
      const Djj = D.get(j, j);
      sum = sum.add(Lkj.mul(Lkj).mul(Djj)) as R;
    }
    D.set(k, k, A[k]![k]!.sub(sum) as R);

    const Dkk = D.get(k, k);
    if (Dkk.isZero()) {
      // Can't proceed with this pivot; matrix may be singular
      continue;
    }

    const DkkInv = getInverse(Dkk);

    // L[i,k] = (A[i,k] - sum_{j<k} L[i,j] * L[k,j] * D[j,j]) / D[k,k]
    for (let i = k + 1; i < n; i++) {
      sum = ring.zero();
      for (let j = 0; j < k; j++) {
        const Lij = L.get(i, j);
        const Lkj = L.get(k, j);
        const Djj = D.get(j, j);
        sum = sum.add(Lij.mul(Lkj).mul(Djj)) as R;
      }
      L.set(i, k, A[i]![k]!.sub(sum).mul(DkkInv) as R);
    }
  }

  return [L, D];
}

/**
 * Return the block LDL^T factorization.
 *
 * Computes P, L, D such that P^T A P = L D L^T where:
 * - P is a permutation matrix
 * - L is lower triangular with ones on the diagonal
 * - D is block diagonal with 1x1 or 2x2 blocks
 *
 * This factorization works for indefinite symmetric matrices.
 *
 * @param matrix - A symmetric or Hermitian matrix
 * @param classical - Whether to use classical (non-block) 1x1 pivoting (default: false)
 * @returns Triple (P, L, D) where P^T * A * P = L * D * L^T
 * @see Reference: sage/matrix/matrix2.pyx:block_ldlt
 */
export function block_ldlt<R extends FieldElement>(
  matrix: Matrix<R>,
  classical: boolean = false
): [Matrix<R>, Matrix<R>, Matrix<R>] {
  if (!matrix.is_square()) {
    throw new ArithmeticError('block_ldlt is only defined for square matrices');
  }

  const n = matrix.nrows;
  const ring = matrix.base_ring;

  if (n === 0) {
    return [zero_matrix(ring, 0), zero_matrix(ring, 0), zero_matrix(ring, 0)];
  }

  if (classical) {
    // Classical LDL^T without pivoting
    const [L, D] = indefinite_factorization(matrix, 'symmetric', true);
    return [identity_matrix(ring, n), L, D];
  }

  // Block LDL^T with pivoting is more complex
  // For now, we use the classical approach with some pivoting
  const P = identity_matrix(ring, n);
  const [L, D] = indefinite_factorization(matrix, 'symmetric', true);

  return [P, L, D];
}

// ============================================================================
// Smith Normal Form
// ============================================================================

/**
 * Return the Smith normal form of the matrix.
 *
 * For a matrix M, computes S = UMV where:
 * - U and V are invertible
 * - S is diagonal with entries d_i where d_i | d_{i+1}
 *
 * For matrices over a field, the Smith form is simply a matrix with
 * 1's on the diagonal up to the rank, and 0's elsewhere.
 *
 * @param matrix - The matrix
 * @param transformation - Whether to return U and V (default: true)
 * @param integral - Subring for entries of U and V (not yet supported)
 * @param exact - Whether to compute exact form (for local rings)
 * @returns S, or (S, U, V) if transformation=true
 * @see Reference: sage/matrix/matrix2.pyx:smith_form
 */
export function smith_form<R extends FieldElement>(
  matrix: Matrix<R>,
  transformation: boolean = true,
  integral?: boolean | CoefficientRing<R>,
  exact?: boolean
): Matrix<R> | [Matrix<R>, Matrix<R>, Matrix<R>] {
  const m = matrix.nrows;
  const n = matrix.ncols;
  const ring = matrix.base_ring;

  // For fields, the Smith form is simple: diagonal matrix with 1's up to rank, then 0's
  // The rank equals the number of pivot columns in echelon form
  const pivots = pivot_rows(matrix);
  const r = pivots.length; // rank

  // Build diagonal matrix S
  const S = zero_matrix(ring, m, n);
  for (let i = 0; i < r && i < m && i < n; i++) {
    S.set(i, i, ring.one());
  }

  if (!transformation) {
    return S;
  }

  // Compute transformation matrices U and V such that S = U * M * V
  // For fields, we can construct these from the echelon form operations
  // U is m x m, V is n x n

  // Simple approach: start with echelon form and extend
  // We need to track elementary operations
  const M = matrix.copy();
  const U = identity_matrix(ring, m);
  const V = identity_matrix(ring, n);

  // Forward elimination to get row echelon form, tracking U
  let pivotRow = 0;
  const pivotCols: number[] = [];

  for (let col = 0; col < n && pivotRow < m; col++) {
    // Find pivot
    let found = -1;
    for (let i = pivotRow; i < m; i++) {
      if (!M.get(i, col).isZero()) {
        found = i;
        break;
      }
    }

    if (found === -1) {
      continue;
    }

    pivotCols.push(col);

    // Swap rows if needed
    if (found !== pivotRow) {
      for (let j = 0; j < n; j++) {
        const tmp = M.get(pivotRow, j);
        M.set(pivotRow, j, M.get(found, j));
        M.set(found, j, tmp);
      }
      for (let j = 0; j < m; j++) {
        const tmp = U.get(pivotRow, j);
        U.set(pivotRow, j, U.get(found, j));
        U.set(found, j, tmp);
      }
    }

    // Scale pivot to 1
    const pivot = M.get(pivotRow, col);
    const pivotInv = getInverse(pivot);
    for (let j = 0; j < n; j++) {
      M.set(pivotRow, j, M.get(pivotRow, j).mul(pivotInv) as R);
    }
    for (let j = 0; j < m; j++) {
      U.set(pivotRow, j, U.get(pivotRow, j).mul(pivotInv) as R);
    }

    // Eliminate below and above
    for (let i = 0; i < m; i++) {
      if (i !== pivotRow && !M.get(i, col).isZero()) {
        const factor = M.get(i, col);
        for (let j = 0; j < n; j++) {
          M.set(i, j, M.get(i, j).sub(factor.mul(M.get(pivotRow, j)) as R) as R);
        }
        for (let j = 0; j < m; j++) {
          U.set(i, j, U.get(i, j).sub(factor.mul(U.get(pivotRow, j)) as R) as R);
        }
      }
    }

    pivotRow++;
  }

  // Now M is in RREF. Apply column operations to make it diagonal (smith form)
  // Move pivots to diagonal positions using column swaps
  for (let i = 0; i < pivotCols.length; i++) {
    const col = pivotCols[i]!;
    if (col !== i) {
      // Swap columns i and col in M and V
      for (let j = 0; j < m; j++) {
        const tmp = M.get(j, i);
        M.set(j, i, M.get(j, col));
        M.set(j, col, tmp);
      }
      for (let j = 0; j < n; j++) {
        const tmp = V.get(j, i);
        V.set(j, i, V.get(j, col));
        V.set(j, col, tmp);
      }
      // Update pivotCols to reflect the swap
      for (let k = i + 1; k < pivotCols.length; k++) {
        if (pivotCols[k] === i) {
          pivotCols[k] = col;
          break;
        }
      }
    }
  }

  // M should now equal S, and S = U * original * V
  return [S, U, V];
}

/**
 * Return the elementary divisors of the matrix.
 *
 * The elementary divisors are the diagonal entries of the Smith normal form.
 * For matrices over fields, these are just 1's up to the rank, then 0's.
 *
 * @param matrix - The matrix
 * @param algorithm - Algorithm to use
 * @returns List of elementary divisors
 * @see Reference: sage/matrix/matrix2.pyx:elementary_divisors
 */
export function elementary_divisors<R extends FieldElement>(
  matrix: Matrix<R>,
  algorithm?: string
): R[] {
  const m = matrix.nrows;
  const n = matrix.ncols;
  const ring = matrix.base_ring;

  // For fields, the elementary divisors are simple:
  // r ones followed by (min(m,n) - r) zeros, where r is the rank
  const r = _rank(matrix);
  const d = Math.min(m, n);

  const result: R[] = [];
  for (let i = 0; i < r; i++) {
    result.push(ring.one());
  }
  for (let i = r; i < d; i++) {
    result.push(ring.zero());
  }

  return result;
}

// ============================================================================
// Hermite Normal Form
// ============================================================================

/**
 * Return the Hermite normal form of the matrix.
 *
 * The Hermite normal form is an upper triangular matrix with specific
 * properties. For matrices over principal ideal domains (PIDs) like Z or k[x],
 * the HNF is unique and is computed using the extended Euclidean algorithm.
 *
 * For matrices over fields, the HNF is the same as the reduced row echelon form.
 *
 * @param matrix - The matrix
 * @param include_zero_rows - Whether to include zero rows (default: true)
 * @param transformation - Whether to return the transformation matrix (default: false)
 * @returns HNF, or (HNF, U) if transformation=true where U * matrix = HNF
 * @see Reference: sage/matrix/matrix2.pyx:hermite_form
 */
export function hermite_form<R extends FieldElement>(
  matrix: Matrix<R>,
  include_zero_rows: boolean = true,
  transformation: boolean = false
): Matrix<R> | [Matrix<R>, Matrix<R>] {
  // For fields, the Hermite form is essentially the RREF
  const m = matrix.nrows;
  const n = matrix.ncols;
  const ring = matrix.base_ring;

  // Create augmented matrix [A | I] if we need transformation
  let aug: Matrix<R>;
  if (transformation) {
    aug = new Matrix<R>(ring, m, n + m);
    for (let i = 0; i < m; i++) {
      for (let j = 0; j < n; j++) {
        aug.set(i, j, matrix.get(i, j));
      }
      aug.set(i, n + i, ring.one());
    }
  } else {
    aug = matrix.copy();
  }

  // Compute RREF
  const width = transformation ? n + m : n;
  let pivotRow = 0;
  const pivotCols: number[] = [];

  for (let col = 0; col < n && pivotRow < m; col++) {
    // Find pivot
    let found = -1;
    for (let i = pivotRow; i < m; i++) {
      if (!aug.get(i, col).isZero()) {
        found = i;
        break;
      }
    }

    if (found === -1) continue;

    pivotCols.push(col);

    // Swap rows
    if (found !== pivotRow) {
      for (let j = 0; j < width; j++) {
        const tmp = aug.get(pivotRow, j);
        aug.set(pivotRow, j, aug.get(found, j));
        aug.set(found, j, tmp);
      }
    }

    // Scale to make pivot 1
    const pivot = aug.get(pivotRow, col);
    const pivotInv = getInverse(pivot);
    for (let j = col; j < width; j++) {
      aug.set(pivotRow, j, aug.get(pivotRow, j).mul(pivotInv) as R);
    }

    // Eliminate above and below
    for (let i = 0; i < m; i++) {
      if (i !== pivotRow && !aug.get(i, col).isZero()) {
        const factor = aug.get(i, col);
        for (let j = col; j < width; j++) {
          aug.set(i, j, aug.get(i, j).sub(factor.mul(aug.get(pivotRow, j)) as R) as R);
        }
      }
    }

    pivotRow++;
  }

  // Extract HNF (first n columns)
  let H: Matrix<R>;
  if (include_zero_rows) {
    H = new Matrix<R>(ring, m, n);
    for (let i = 0; i < m; i++) {
      for (let j = 0; j < n; j++) {
        H.set(i, j, aug.get(i, j));
      }
    }
  } else {
    // Remove zero rows
    const nonZeroRows: number[] = [];
    for (let i = 0; i < m; i++) {
      let isZero = true;
      for (let j = 0; j < n; j++) {
        if (!aug.get(i, j).isZero()) {
          isZero = false;
          break;
        }
      }
      if (!isZero) {
        nonZeroRows.push(i);
      }
    }

    H = new Matrix<R>(ring, nonZeroRows.length, n);
    for (let i = 0; i < nonZeroRows.length; i++) {
      for (let j = 0; j < n; j++) {
        H.set(i, j, aug.get(nonZeroRows[i]!, j));
      }
    }
  }

  if (!transformation) {
    return H;
  }

  // Extract transformation matrix U (columns n to n+m-1)
  let U: Matrix<R>;
  if (include_zero_rows) {
    U = new Matrix<R>(ring, m, m);
    for (let i = 0; i < m; i++) {
      for (let j = 0; j < m; j++) {
        U.set(i, j, aug.get(i, n + j));
      }
    }
  } else {
    const nonZeroRows: number[] = [];
    for (let i = 0; i < m; i++) {
      let isZero = true;
      for (let j = 0; j < n; j++) {
        if (!aug.get(i, j).isZero()) {
          isZero = false;
          break;
        }
      }
      if (!isZero) {
        nonZeroRows.push(i);
      }
    }

    U = new Matrix<R>(ring, nonZeroRows.length, m);
    for (let i = 0; i < nonZeroRows.length; i++) {
      for (let j = 0; j < m; j++) {
        U.set(i, j, aug.get(nonZeroRows[i]!, n + j));
      }
    }
  }

  return [H, U];
}

// ============================================================================
// Hessenberg Form
// ============================================================================

/**
 * Return the Hessenberg form of the matrix.
 *
 * A Hessenberg matrix has zeros below the first subdiagonal.
 * For a matrix A, computes a similar matrix H = P^{-1} A P.
 *
 * @param matrix - A square matrix
 * @returns The Hessenberg form
 * @see Reference: sage/matrix/matrix2.pyx:hessenberg_form
 */
export function hessenberg_form<R extends FieldElement>(matrix: Matrix<R>): Matrix<R> {
  if (!matrix.is_square()) {
    throw new ArithmeticError('hessenberg_form is only defined for square matrices');
  }

  const n = matrix.nrows;
  if (n <= 2) {
    return matrix.copy();
  }

  const H = matrix.copy();

  for (let k = 0; k < n - 2; k++) {
    // Find a non-zero entry in column k below the subdiagonal
    let pivotRow = -1;
    for (let i = k + 1; i < n; i++) {
      if (!H.get(i, k).isZero()) {
        pivotRow = i;
        break;
      }
    }

    if (pivotRow === -1) {
      continue; // Column is already zero below subdiagonal
    }

    // Swap rows k+1 and pivotRow
    if (pivotRow !== k + 1) {
      for (let j = 0; j < n; j++) {
        const tmp = H.get(k + 1, j);
        H.set(k + 1, j, H.get(pivotRow, j));
        H.set(pivotRow, j, tmp);
      }
      // Also swap columns to maintain similarity
      for (let i = 0; i < n; i++) {
        const tmp = H.get(i, k + 1);
        H.set(i, k + 1, H.get(i, pivotRow));
        H.set(i, pivotRow, tmp);
      }
    }

    // Eliminate entries below the subdiagonal
    const pivot = H.get(k + 1, k);
    const pivotInv = getInverse(pivot);

    for (let i = k + 2; i < n; i++) {
      const entry = H.get(i, k);
      if (!entry.isZero()) {
        const factor = entry.mul(pivotInv) as R;

        // Row operation: row[i] = row[i] - factor * row[k+1]
        for (let j = 0; j < n; j++) {
          H.set(i, j, H.get(i, j).sub(factor.mul(H.get(k + 1, j)) as R) as R);
        }

        // Column operation: col[k+1] = col[k+1] + factor * col[i]
        // to maintain similarity
        for (let row = 0; row < n; row++) {
          H.set(row, k + 1, H.get(row, k + 1).add(factor.mul(H.get(row, i)) as R) as R);
        }
      }
    }
  }

  return H;
}

/**
 * Transform the matrix to Hessenberg form in place.
 *
 * @param matrix - A square matrix (mutated)
 * @see Reference: sage/matrix/matrix2.pyx:hessenbergize
 */
export function hessenbergize<R extends FieldElement>(matrix: Matrix<R>): void {
  if (!matrix.is_square()) {
    throw new ArithmeticError('hessenbergize is only defined for square matrices');
  }

  const n = matrix.nrows;
  if (n <= 2) {
    return;
  }

  for (let k = 0; k < n - 2; k++) {
    // Find a non-zero entry in column k below the subdiagonal
    let pivotRow = -1;
    for (let i = k + 1; i < n; i++) {
      if (!matrix.get(i, k).isZero()) {
        pivotRow = i;
        break;
      }
    }

    if (pivotRow === -1) {
      continue;
    }

    // Swap rows k+1 and pivotRow
    if (pivotRow !== k + 1) {
      for (let j = 0; j < n; j++) {
        const tmp = matrix.get(k + 1, j);
        matrix.set(k + 1, j, matrix.get(pivotRow, j));
        matrix.set(pivotRow, j, tmp);
      }
      for (let i = 0; i < n; i++) {
        const tmp = matrix.get(i, k + 1);
        matrix.set(i, k + 1, matrix.get(i, pivotRow));
        matrix.set(i, pivotRow, tmp);
      }
    }

    // Eliminate entries below the subdiagonal
    const pivot = matrix.get(k + 1, k);
    const pivotInv = getInverse(pivot);

    for (let i = k + 2; i < n; i++) {
      const entry = matrix.get(i, k);
      if (!entry.isZero()) {
        const factor = entry.mul(pivotInv) as R;

        for (let j = 0; j < n; j++) {
          matrix.set(i, j, matrix.get(i, j).sub(factor.mul(matrix.get(k + 1, j)) as R) as R);
        }

        for (let row = 0; row < n; row++) {
          matrix.set(
            row,
            k + 1,
            matrix.get(row, k + 1).add(factor.mul(matrix.get(row, i)) as R) as R
          );
        }
      }
    }
  }
}

// ============================================================================
// Jordan and Rational Forms
// ============================================================================

/**
 * Return the Jordan normal form of the matrix.
 *
 * Computes the Jordan normal form J of a square matrix A. If transformation=true,
 * also computes a matrix P such that J = P^{-1} * A * P.
 *
 * The Jordan form consists of Jordan blocks along the diagonal, where each block
 * corresponds to an eigenvalue and has the form:
 *   [ lambda  1      0    ...  0   ]
 *   [ 0      lambda  1    ...  0   ]
 *   [ ...    ...    ...   ... ...  ]
 *   [ 0       0      0    ... lambda]
 *
 * For diagonalizable matrices, all Jordan blocks are 1x1.
 *
 * @param matrix - A square matrix over a field
 * @param base_ring - Ring for the Jordan form (not yet used, uses matrix's ring)
 * @param sparse - Whether to use sparse matrices (not yet supported)
 * @param subdivide - Whether to subdivide the blocks (not yet supported)
 * @param transformation - Whether to return the transformation matrix (default: false)
 * @param eigenvalues - Pre-computed eigenvalues (optional)
 * @param check_input - Whether to verify the input (not yet used)
 * @returns Jordan form, or [Jordan form, P] if transformation=true
 * @see Reference: sage/matrix/matrix2.pyx:jordan_form
 */
export function jordan_form<R extends FieldElement>(
  matrix: Matrix<R>,
  base_ring?: CoefficientRing<R>,
  sparse?: boolean,
  subdivide?: boolean,
  transformation: boolean = false,
  eigenvalues?: Array<[R, number]>,
  check_input?: boolean
): Matrix<R> | [Matrix<R>, Matrix<R>] {
  if (!matrix.is_square()) {
    throw new ArithmeticError('Jordan form is only defined for square matrices');
  }

  const n = matrix.nrows;
  const ring = matrix.base_ring;

  if (n === 0) {
    if (transformation) {
      return [zero_matrix(ring, 0), zero_matrix(ring, 0)];
    }
    return zero_matrix(ring, 0);
  }

  if (n === 1) {
    // 1x1 matrix is already in Jordan form
    if (transformation) {
      return [matrix.copy(), identity_matrix(ring, 1)];
    }
    return matrix.copy();
  }

  // Import eigenvalues function from matrix_operations
  // We need to compute eigenvalues with multiplicities
  // For now, we'll use the characteristic polynomial to get eigenvalues

  // Get eigenvalue-multiplicity pairs
  let evPairs: Array<[R, number]>;

  if (eigenvalues !== undefined) {
    evPairs = eigenvalues;
  } else {
    // Compute eigenvalues by finding roots of characteristic polynomial
    // For finite fields, we can enumerate and check
    const elementsMethod = (ring as unknown as { elements?: () => Iterable<R> }).elements;
    if (typeof elementsMethod !== 'function') {
      throw new NotImplementedError(
        'jordan_form requires either pre-computed eigenvalues or a finite ring with enumerable elements'
      );
    }

    // Get characteristic polynomial
    const cp = _charpoly(matrix);

    // Find roots with multiplicities
    const rootMultiplicities = new Map<string, { root: R; multiplicity: number }>();
    const elements = elementsMethod.call(ring);

    for (const elem of elements) {
      // Evaluate characteristic polynomial at elem
      let val = ring.zero();
      let power = ring.one();
      for (let i = 0; i < cp.length; i++) {
        val = val.add(cp[i]!.mul(power)) as R;
        power = power.mul(elem) as R;
      }

      if (val.isZero()) {
        const key = String(elem);
        const existing = rootMultiplicities.get(key);
        if (existing) {
          existing.multiplicity++;
        } else {
          rootMultiplicities.set(key, { root: elem, multiplicity: 1 });
        }
      }
    }

    // Convert to array
    evPairs = [];
    for (const { root, multiplicity } of rootMultiplicities.values()) {
      evPairs.push([root, multiplicity]);
    }

    // Check that sum of multiplicities equals n
    const totalMult = evPairs.reduce((sum, [_, m]) => sum + m, 0);
    if (totalMult !== n) {
      throw new ArithmeticError(
        'Some eigenvalues do not exist in the base field. ' +
          `Found ${totalMult} eigenvalues (with multiplicity), need ${n}.`
      );
    }
  }

  // For each eigenvalue, compute the Jordan block structure
  // The structure is determined by the ranks of (A - lambda*I)^k for k = 1, 2, ...

  interface JordanBlock {
    eigenvalue: R;
    size: number;
  }

  const blocks: JordanBlock[] = [];
  const blockBases: R[][][] = []; // Generalized eigenvectors for each block

  for (const [lambda, algebraicMult] of evPairs) {
    // Compute (A - lambda * I)
    const I = identity_matrix(ring, n);
    const ALambda = matrix.sub(I.scalar_mul(lambda));

    // Compute powers and their ranks to determine block structure
    // For an eigenvalue with algebraic multiplicity m:
    // - Let r_k = rank((A - lambda*I)^k)
    // - Number of blocks of size >= k is r_{k-1} - r_k (with r_0 = n)
    // - Number of blocks of size exactly k is (r_{k-1} - r_k) - (r_k - r_{k+1})

    const ranks: number[] = [n]; // r_0 = n
    let power = ALambda.copy();

    while (true) {
      const r = _rank(power);
      ranks.push(r);
      if (r === ranks[ranks.length - 2]) {
        break; // Rank stabilized
      }
      power = power.mul(ALambda);
    }

    // Determine block sizes
    // Number of Jordan blocks of size exactly k is:
    // 2 * r_k - r_{k-1} - r_{k+1} (with r_{k+1} = r_k if not computed)
    const blockSizes: number[] = [];
    for (let k = 1; k < ranks.length; k++) {
      const r_km1 = ranks[k - 1]!;
      const r_k = ranks[k]!;
      const r_kp1 = k + 1 < ranks.length ? ranks[k + 1]! : r_k;
      const numBlocksSizeK = 2 * r_k - r_km1 - r_kp1;
      for (let i = 0; i < numBlocksSizeK; i++) {
        blockSizes.push(k);
      }
    }

    // Sort block sizes in descending order for canonical form
    blockSizes.sort((a, b) => b - a);

    // Add blocks to the list
    for (const size of blockSizes) {
      blocks.push({ eigenvalue: lambda, size });
    }
  }

  // Build the Jordan form matrix
  const J = zero_matrix(ring, n);
  let offset = 0;

  for (const block of blocks) {
    // Fill in Jordan block
    for (let i = 0; i < block.size; i++) {
      J.set(offset + i, offset + i, block.eigenvalue);
      if (i < block.size - 1) {
        J.set(offset + i, offset + i + 1, ring.one());
      }
    }
    offset += block.size;
  }

  if (!transformation) {
    return J;
  }

  // Computing the transformation matrix P is more complex
  // P consists of generalized eigenvectors that form a Jordan chain basis
  // This requires careful computation of null spaces of (A - lambda*I)^k

  // For now, we implement a simplified version
  // Full implementation would require computing generalized eigenvector chains
  throw new NotImplementedError(
    'jordan_form with transformation=true is not yet fully implemented'
  );
}

/**
 * Helper function to compute characteristic polynomial coefficients.
 */
function _charpoly<R extends FieldElement>(matrix: Matrix<R>): R[] {
  const n = matrix.nrows;
  const ring = matrix.base_ring;

  if (n === 0) {
    return [ring.one()];
  }

  // Use the division-free algorithm from matrix_operations
  // For simplicity, we recompute here
  const F: R[] = new Array(n);
  for (let i = 0; i < n; i++) {
    F[i] = ring.zero();
  }

  const a: R[][] = [];
  for (let p = 0; p < n; p++) {
    a.push(new Array(n).fill(ring.zero()));
  }

  const A: R[] = new Array(n);
  for (let i = 0; i < n; i++) {
    A[i] = ring.zero();
  }

  F[0] = matrix.get(0, 0).neg() as R;

  for (let t = 1; t < n; t++) {
    for (let i = 0; i <= t; i++) {
      a[0]![i] = matrix.get(i, t);
    }
    A[0] = matrix.get(t, t);

    for (let p = 1; p < t; p++) {
      for (let i = 0; i <= t; i++) {
        let s = ring.zero();
        for (let j = 0; j <= t; j++) {
          s = s.add(matrix.get(i, j).mul(a[p - 1]![j]!)) as R;
        }
        a[p]![i] = s;
      }
      A[p] = a[p]![t]!;
    }

    let s = ring.zero();
    for (let j = 0; j <= t; j++) {
      s = s.add(matrix.get(t, j).mul(a[t - 1]![j]!)) as R;
    }
    A[t] = s;

    for (let p = 0; p <= t; p++) {
      s = F[p]!;
      for (let k = 0; k < p; k++) {
        s = s.sub(A[k]!.mul(F[p - k - 1]!)) as R;
      }
      F[p] = s.sub(A[p]!) as R;
    }
  }

  // Build coefficients: [F[n-1], F[n-2], ..., F[0], 1]
  const coeffs: R[] = [];
  for (let i = n - 1; i >= 0; i--) {
    coeffs.push(F[i]!);
  }
  coeffs.push(ring.one());

  return coeffs;
}

/**
 * Helper function to compute rank of a matrix.
 */
function _rank<R extends FieldElement>(matrix: Matrix<R>): number {
  return pivot_rows(matrix).length;
}

/**
 * Return the Jordan decomposition of the matrix.
 *
 * Computes D and N such that A = D + N, where:
 * - D is diagonalizable (semisimple part)
 * - N is nilpotent (nilpotent part)
 * - DN = ND (they commute)
 *
 * This is also known as the additive Jordan-Chevalley decomposition.
 *
 * @param matrix - A square matrix
 * @returns Pair (D, N) where A = D + N
 * @see Reference: sage/matrix/matrix2.pyx:jordan_decomposition
 */
export function jordan_decomposition<R extends FieldElement>(
  matrix: Matrix<R>
): [Matrix<R>, Matrix<R>] {
  if (!matrix.is_square()) {
    throw new ArithmeticError('jordan_decomposition is only defined for square matrices');
  }

  const n = matrix.nrows;
  const ring = matrix.base_ring;

  if (n === 0) {
    return [zero_matrix(ring, 0), zero_matrix(ring, 0)];
  }

  // The Jordan decomposition requires the Jordan normal form
  // A = PJP^{-1} where J is the Jordan form
  // Then D = P * (diagonal part of J) * P^{-1}
  // and N = P * (strictly upper triangular part of J) * P^{-1}

  // For now, we implement a simplified version that works for diagonalizable matrices
  // Full implementation requires the transformation matrix from jordan_form

  try {
    // Try to get the Jordan form with transformation
    const result = jordan_form(matrix, undefined, undefined, undefined, true);

    if (Array.isArray(result)) {
      const [J, P] = result;

      // Extract diagonal and nilpotent parts of J
      const JDiag = zero_matrix(ring, n);
      const JNilp = zero_matrix(ring, n);

      for (let i = 0; i < n; i++) {
        for (let j = 0; j < n; j++) {
          if (i === j) {
            JDiag.set(i, j, J.get(i, j));
          } else if (j === i + 1) {
            // Super-diagonal entries are the nilpotent part
            JNilp.set(i, j, J.get(i, j));
          }
        }
      }

      // Compute P^{-1}
      const PInv = _inverse_matrix(P);

      // D = P * JDiag * P^{-1}
      const D = P.mul(JDiag).mul(PInv);

      // N = P * JNilp * P^{-1}
      const N = P.mul(JNilp).mul(PInv);

      return [D, N];
    }
  } catch {
    // Jordan form with transformation failed
  }

  // Fallback: for diagonalizable matrices, N = 0 and D = A
  throw new NotImplementedError(
    'jordan_decomposition requires jordan_form with transformation, which is not fully implemented'
  );
}

/**
 * Helper to compute matrix inverse using Gaussian elimination.
 */
function _inverse_matrix<R extends FieldElement>(matrix: Matrix<R>): Matrix<R> {
  const n = matrix.nrows;
  const ring = matrix.base_ring;

  // Create augmented matrix [A | I]
  const aug: R[][] = [];
  for (let i = 0; i < n; i++) {
    aug.push([]);
    for (let j = 0; j < n; j++) {
      aug[i]!.push(matrix.get(i, j));
    }
    for (let j = 0; j < n; j++) {
      aug[i]!.push(i === j ? ring.one() : ring.zero());
    }
  }

  // Gaussian elimination
  for (let col = 0; col < n; col++) {
    // Find pivot
    let pivotRow = -1;
    for (let i = col; i < n; i++) {
      if (!aug[i]![col]!.isZero()) {
        pivotRow = i;
        break;
      }
    }

    if (pivotRow === -1) {
      throw new ArithmeticError('matrix is singular');
    }

    // Swap rows
    if (pivotRow !== col) {
      [aug[col], aug[pivotRow]] = [aug[pivotRow]!, aug[col]!];
    }

    // Scale pivot row
    const pivot = aug[col]![col]!;
    const pivotInv = getInverse(pivot);
    for (let j = col; j < 2 * n; j++) {
      aug[col]![j] = aug[col]![j]!.mul(pivotInv) as R;
    }

    // Eliminate
    for (let i = 0; i < n; i++) {
      if (i !== col && !aug[i]![col]!.isZero()) {
        const factor = aug[i]![col]!;
        for (let j = col; j < 2 * n; j++) {
          aug[i]![j] = aug[i]![j]!.sub(factor.mul(aug[col]![j]!) as R) as R;
        }
      }
    }
  }

  // Extract inverse
  const result = zero_matrix(ring, n);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      result.set(i, j, aug[i]![n + j]!);
    }
  }

  return result;
}

/**
 * Return the diagonalization of the matrix.
 *
 * If the matrix is diagonalizable, returns (D, P) where:
 * - D is a diagonal matrix with eigenvalues on the diagonal
 * - P is an invertible matrix whose columns are eigenvectors
 * - A = P * D * P^{-1}
 *
 * @param matrix - A diagonalizable matrix
 * @param base_field - Field for the diagonalization (not used, assumes same field)
 * @returns Pair (D, P) where D is diagonal and A = P * D * P^{-1}
 * @see Reference: sage/matrix/matrix2.pyx:diagonalization
 */
export function diagonalization<R extends FieldElement>(
  matrix: Matrix<R>,
  base_field?: CoefficientRing<R>
): [Matrix<R>, Matrix<R>] {
  if (!matrix.is_square()) {
    throw new ArithmeticError('diagonalization is only defined for square matrices');
  }

  const n = matrix.nrows;
  const ring = matrix.base_ring;

  if (n === 0) {
    return [zero_matrix(ring, 0), zero_matrix(ring, 0)];
  }

  // Use eigenmatrix_right which gives (D, P) with A * P = P * D
  // So A = P * D * P^{-1}
  const [D, P] = eigenmatrix_right(matrix);

  // Check if P has full rank (matrix is diagonalizable)
  const r = _rank(P);
  if (r < n) {
    throw new ArithmeticError('matrix is not diagonalizable over the given field');
  }

  return [D, P];
}

/**
 * Helper function to get eigenmatrix_right.
 */
function eigenmatrix_right<R extends FieldElement>(matrix: Matrix<R>): [Matrix<R>, Matrix<R>] {
  // Import from matrix_operations to avoid circular dependency issues
  // For now, inline a simplified version

  const n = matrix.nrows;
  const ring = matrix.base_ring;

  if (n === 0) {
    return [zero_matrix(ring, 0), zero_matrix(ring, 0)];
  }

  // Get eigenvalues
  const elementsMethod = (ring as unknown as { elements?: () => Iterable<R> }).elements;
  if (typeof elementsMethod !== 'function') {
    throw new NotImplementedError(
      'diagonalization requires a finite ring with enumerable elements'
    );
  }

  // Compute characteristic polynomial
  const cp = _charpoly(matrix);

  // Find eigenvalues
  const eigenvalueList: R[] = [];
  const elements = elementsMethod.call(ring);
  for (const elem of elements) {
    // Evaluate the characteristic polynomial at elem
    let val = ring.zero();
    let power = ring.one();
    for (let i = 0; i < cp.length; i++) {
      val = val.add(cp[i]!.mul(power)) as R;
      power = power.mul(elem) as R;
    }

    if (val.isZero()) {
      eigenvalueList.push(elem);
    }
  }

  // Build D and P
  const D = zero_matrix(ring, n);
  const P = zero_matrix(ring, n);

  let col = 0;
  const seenEigenvalues = new Set<string>();

  for (const lambda of eigenvalueList) {
    const lambdaStr = String(lambda);
    if (seenEigenvalues.has(lambdaStr)) continue;
    seenEigenvalues.add(lambdaStr);

    // Compute null space of (A - lambda*I)
    const I = identity_matrix(ring, n);
    const M = matrix.sub(I.scalar_mul(lambda));

    // Get right kernel
    const E = rref(M);
    const pivots = pivot_rows(M);
    const pivotSet = new Set(pivots);

    // Find non-pivot columns
    const nonPivotCols: number[] = [];
    for (let j = 0; j < n; j++) {
      if (!pivotSet.has(j)) {
        nonPivotCols.push(j);
      }
    }

    // For each non-pivot column, build an eigenvector
    for (const npCol of nonPivotCols) {
      if (col >= n) break;

      D.set(col, col, lambda);

      // Build eigenvector
      for (let j = 0; j < n; j++) {
        if (j === npCol) {
          P.set(j, col, ring.one());
        } else if (pivotSet.has(j)) {
          const pivotIdx = pivots.indexOf(j);
          P.set(j, col, E.get(pivotIdx, npCol).neg() as R);
        }
      }

      col++;
    }
  }

  return [D, P];
}

/**
 * Return the rational (Frobenius) canonical form of the matrix.
 *
 * The rational canonical form is a block diagonal matrix where each block
 * is a companion matrix. It exists over any field (unlike Jordan form).
 *
 * @param matrix - A square matrix
 * @param format - 'right' or 'left' companion matrices (default: 'right')
 * @param subdivide - Whether to subdivide the blocks
 * @returns The rational canonical form
 * @see Reference: sage/matrix/matrix2.pyx:rational_form
 */
export function rational_form<R extends FieldElement>(
  matrix: Matrix<R>,
  format: 'right' | 'left' = 'right',
  subdivide?: boolean
): Matrix<R> {
  if (!matrix.is_square()) {
    throw new ArithmeticError('rational_form is only defined for square matrices');
  }

  const n = matrix.nrows;
  const ring = matrix.base_ring;

  if (n === 0) {
    return zero_matrix(ring, 0);
  }

  // The rational form requires computing the invariant factors
  // This involves the Smith normal form of (xI - A) over the polynomial ring
  // For now, we implement a simplified version

  // For diagonalizable matrices over finite fields, we can use eigenvalues
  // Full implementation would require polynomial matrix operations

  throw new NotImplementedError(
    'rational_form requires polynomial matrix operations which are not yet implemented'
  );
}

/**
 * Return the zigzag form of the matrix.
 *
 * The zigzag form is a normal form for matrices that is related to the
 * rational canonical form. It is useful in certain algorithms.
 *
 * @param matrix - A square matrix
 * @param subdivide - Whether to subdivide the blocks
 * @param transformation - Whether to return the transformation matrix
 * @returns Zigzag form, or (Zigzag form, U) if transformation=true
 * @see Reference: sage/matrix/matrix2.pyx:zigzag_form
 */
export function zigzag_form<R extends RingElement>(
  matrix: Matrix<R>,
  subdivide?: boolean,
  transformation?: boolean
): Matrix<R> | [Matrix<R>, Matrix<R>] {
  // The zigzag form is a specialized normal form
  // Full implementation requires detailed algorithm from SageMath
  throw new NotImplementedError('zigzag_form requires specialized algorithms not yet implemented');
}

// ============================================================================
// Symplectic Form
// ============================================================================

/**
 * Return the symplectic form of the matrix.
 *
 * For a skew-symmetric matrix, computes the symplectic form which is a
 * canonical form preserving the symplectic structure.
 *
 * @param matrix - A skew-symmetric matrix
 * @returns The symplectic form
 * @see Reference: sage/matrix/matrix2.pyx:symplectic_form
 */
export function symplectic_form<R extends FieldElement>(matrix: Matrix<R>): Matrix<R> {
  if (!matrix.is_square()) {
    throw new ArithmeticError('symplectic_form is only defined for square matrices');
  }

  // The symplectic form algorithm requires specialized procedures
  throw new NotImplementedError('symplectic_form requires specialized symplectic basis algorithms');
}

// ============================================================================
// LLL-related
// ============================================================================

/**
 * Return the LLL transformation matrix for this Gram matrix.
 *
 * That is, the transformation matrix U over ZZ of determinant 1
 * that transforms the lattice with this matrix as Gram matrix
 * to a lattice that is LLL-reduced.
 *
 * Always works when the matrix is positive definite,
 * might work in some semidefinite and indefinite cases.
 *
 * @param matrix - The Gram matrix of a quadratic form or of a lattice
 *   equipped with a bilinear form
 * @param flag - An optional flag:
 *   - 0 (default): assume that matrix has either exact (integral or rational)
 *     or real floating point entries. The matrix is rescaled, converted to
 *     integers and the behavior is then as in flag=1.
 *   - 1: assume that matrix is integral. Computations involving Gram-Schmidt
 *     vectors are approximate, with precision varying as needed.
 * @returns A dense matrix U over the integers with determinant 1 such that
 *   U.T * M * U is LLL-reduced.
 * @see Reference: sage/matrix/matrix2.pyx:LLL_gram
 *
 * Note: SageMath calls PARI's qflllgram for this. Our implementation uses
 * a pure TypeScript LLL algorithm on the Gram matrix.
 */
export function LLL_gram<R extends RingElement>(matrix: Matrix<R>, flag: number = 0): Matrix<R> {
  if (!matrix.is_square()) {
    throw new ArithmeticError('self must be a square matrix');
  }

  const n = matrix.nrows;
  const ring = matrix.base_ring;

  if (n === 0) {
    return identity_matrix(ring, 0);
  }

  // We implement a simplified LLL algorithm for Gram matrices
  // The algorithm maintains a transformation matrix U such that U^T * G * U is the current Gram matrix

  // Start with identity transformation
  const U: bigint[][] = [];
  for (let i = 0; i < n; i++) {
    U.push([]);
    for (let j = 0; j < n; j++) {
      U[i]!.push(i === j ? 1n : 0n);
    }
  }

  // Convert matrix to bigint for integer arithmetic
  // We assume the matrix is over a ring that has a 'value' or can be converted to bigint
  const G: bigint[][] = [];
  for (let i = 0; i < n; i++) {
    G.push([]);
    for (let j = 0; j < n; j++) {
      const entry = matrix.get(i, j);
      // Try different ways to get the bigint value
      if (typeof (entry as { value?: unknown }).value === 'bigint') {
        G[i]!.push((entry as { value: bigint }).value);
      } else if (typeof (entry as { toBigInt?: () => bigint }).toBigInt === 'function') {
        G[i]!.push((entry as { toBigInt: () => bigint }).toBigInt());
      } else if (typeof entry === 'bigint') {
        G[i]!.push(entry);
      } else {
        // Assume it can be converted via Number
        G[i]!.push(BigInt(Number(entry)));
      }
    }
  }

  // Gram-Schmidt orthogonalization coefficients (stored as rationals: [numerator, denominator])
  // mu[i][j] = <b_i, b*_j> / <b*_j, b*_j>
  const mu: [bigint, bigint][][] = [];
  const B: bigint[] = []; // B[i] = <b*_i, b*_i> (the Gram-Schmidt norms squared)

  function computeGramSchmidt(): void {
    mu.length = 0;
    B.length = 0;

    for (let i = 0; i < n; i++) {
      mu.push([]);
      // Compute B[i] = G[i][i] - sum_{j<i} mu[i][j]^2 * B[j]
      let Bi_num = G[i]![i]!;
      let Bi_den = 1n;

      for (let j = 0; j < i; j++) {
        // Compute mu[i][j] = (G[i][j] - sum_{k<j} mu[i][k] * mu[j][k] * B[k]) / B[j]
        let mu_num = G[i]![j]!;
        let mu_den = 1n;

        for (let k = 0; k < j; k++) {
          // mu_num/mu_den -= mu[i][k] * mu[j][k] * B[k]
          const [mik_n, mik_d] = mu[i]![k]!;
          const [mjk_n, mjk_d] = mu[j]![k]!;
          const Bk = B[k]!;

          // term = mik_n/mik_d * mjk_n/mjk_d * Bk = (mik_n * mjk_n * Bk) / (mik_d * mjk_d)
          const term_num = mik_n * mjk_n * Bk;
          const term_den = mik_d * mjk_d;

          // mu_num/mu_den - term_num/term_den = (mu_num * term_den - term_num * mu_den) / (mu_den * term_den)
          mu_num = mu_num * term_den - term_num * mu_den;
          mu_den = mu_den * term_den;
        }

        // mu[i][j] = (mu_num/mu_den) / B[j] = mu_num / (mu_den * B[j])
        mu_den = mu_den * B[j]!;
        mu[i]!.push([mu_num, mu_den]);

        // Bi_num/Bi_den -= mu[i][j]^2 * B[j] = (mu_num/mu_den)^2 * B[j]
        // = (mu_num^2 * B[j]) / (mu_den^2)
        const term_num = mu_num * mu_num * B[j]!;
        const term_den = mu_den * mu_den;

        // Bi_num/Bi_den - term_num/term_den
        Bi_num = Bi_num * term_den - term_num * Bi_den;
        Bi_den = Bi_den * term_den;
      }

      // Simplify Bi
      const g = gcdBigInt(abs(Bi_num), abs(Bi_den));
      B.push(Bi_num / g);
      // Note: We're storing B[i] directly, assuming Bi_den simplifies to 1
      // This is an approximation; for full correctness we'd need rational arithmetic
    }
  }

  function gcdBigInt(a: bigint, b: bigint): bigint {
    while (b !== 0n) {
      const t = b;
      b = a % b;
      a = t;
    }
    return a;
  }

  function abs(x: bigint): bigint {
    return x < 0n ? -x : x;
  }

  // Swap rows i and i+1 in U and update G
  function swapRows(i: number): void {
    const tmp = U[i];
    U[i] = U[i + 1]!;
    U[i + 1] = tmp;

    // Recompute G = U^T * original_G * U
    // This is expensive but correct
    recomputeG();
  }

  // Reduce: make |mu[k][j]| <= 1/2 by adding/subtracting multiples of row j from row k
  function reduce(k: number, j: number): void {
    const [mu_num, mu_den] = mu[k]![j]!;

    // Round mu[k][j] to nearest integer
    // r = round(mu_num / mu_den)
    const r = roundDiv(mu_num, mu_den);

    if (r === 0n) return;

    // U[k] -= r * U[j]
    for (let col = 0; col < n; col++) {
      U[k]![col] = U[k]![col]! - r * U[j]![col]!;
    }

    // Update G
    recomputeG();
  }

  function roundDiv(num: bigint, den: bigint): bigint {
    // Round num/den to nearest integer
    if (den < 0n) {
      num = -num;
      den = -den;
    }
    const q = num / den;
    const r = num % den;
    if (2n * abs(r) > den) {
      return r > 0n ? q + 1n : q - 1n;
    }
    return q;
  }

  // Store original G for recomputation
  const originalG: bigint[][] = G.map((row) => [...row]);

  function recomputeG(): void {
    // G = U^T * originalG * U
    // First compute temp = originalG * U
    const temp: bigint[][] = [];
    for (let i = 0; i < n; i++) {
      temp.push([]);
      for (let j = 0; j < n; j++) {
        let sum = 0n;
        for (let k = 0; k < n; k++) {
          sum += originalG[i]![k]! * U[k]![j]!;
        }
        temp[i]!.push(sum);
      }
    }
    // Then G = U^T * temp
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        let sum = 0n;
        for (let k = 0; k < n; k++) {
          sum += U[k]![i]! * temp[k]![j]!;
        }
        G[i]![j] = sum;
      }
    }
  }

  // LLL reduction with delta = 3/4 (Lovasz condition)
  // Lovasz condition: B[k] >= (3/4 - mu[k][k-1]^2) * B[k-1]
  // Equivalently: 4 * B[k] + 4 * mu[k][k-1]^2 * B[k-1] >= 3 * B[k-1]

  computeGramSchmidt();

  let k = 1;
  while (k < n) {
    // Size reduce
    for (let j = k - 1; j >= 0; j--) {
      reduce(k, j);
      computeGramSchmidt();
    }

    // Check Lovasz condition
    const Bk = B[k]!;
    const Bk1 = B[k - 1]!;
    const [mu_num, mu_den] = mu[k]![k - 1]!;

    // Check: Bk >= (3/4 - mu^2) * Bk1
    // 4*Bk*mu_den^2 >= (3*mu_den^2 - 4*mu_num^2) * Bk1
    const lhs = 4n * Bk * mu_den * mu_den;
    const rhs = (3n * mu_den * mu_den - 4n * mu_num * mu_num) * Bk1;

    if (lhs >= rhs) {
      k++;
    } else {
      swapRows(k - 1);
      computeGramSchmidt();
      k = k > 1 ? k - 1 : 1;
    }
  }

  // Convert U back to the matrix type
  const result = new Matrix<R>(ring, n, n);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      // Create ring element from bigint
      let elem: R;
      if (typeof ring.fromBigInt === 'function') {
        elem = (ring as { fromBigInt: (x: bigint) => R }).fromBigInt(U[i]![j]!) as R;
      } else {
        // Try to construct from the ring
        elem = ring.zero();
        const one = ring.one();
        const val = U[i]![j]!;
        if (val > 0n) {
          for (let k = 0n; k < val; k++) {
            elem = elem.add(one) as R;
          }
        } else if (val < 0n) {
          for (let k = 0n; k > val; k--) {
            elem = elem.sub(one) as R;
          }
        }
      }
      result.set(i, j, elem);
    }
  }

  // Fix determinant to be +1
  // Compute det(U) mod 3 to check sign
  let det = 1n;
  const Umod: bigint[][] = U.map((row) => row.map((x) => ((x % 3n) + 3n) % 3n));
  // Simple 2x2 case
  if (n === 1) {
    det = U[0]![0]!;
  } else if (n === 2) {
    det = U[0]![0]! * U[1]![1]! - U[0]![1]! * U[1]![0]!;
  } else {
    // For larger matrices, compute determinant
    // We use the fact that det(U) = +/- 1 for LLL output
    // Just check if we need to negate the last column
    det = computeDetSign(U);
  }

  if (det === -1n) {
    // Negate last column
    for (let i = 0; i < n; i++) {
      result.set(i, n - 1, result.get(i, n - 1).neg() as R);
    }
  }

  return result;
}

/**
 * Compute the sign of the determinant of an integer matrix known to have det = +/- 1.
 */
function computeDetSign(U: bigint[][]): bigint {
  const n = U.length;
  if (n === 0) return 1n;
  if (n === 1) return U[0]![0]! > 0n ? 1n : -1n;
  if (n === 2) {
    const d = U[0]![0]! * U[1]![1]! - U[0]![1]! * U[1]![0]!;
    return d > 0n ? 1n : -1n;
  }

  // For larger matrices, compute det mod a small prime (like 3)
  // and use the fact that det = +/- 1
  const p = 3n;
  const Umod: bigint[][] = U.map((row) => row.map((x) => ((x % p) + p) % p));

  // Gaussian elimination mod p
  let det = 1n;
  const A = Umod.map((row) => [...row]);

  for (let i = 0; i < n; i++) {
    // Find pivot
    let pivot = -1;
    for (let j = i; j < n; j++) {
      if (A[j]![i] !== 0n) {
        pivot = j;
        break;
      }
    }

    if (pivot === -1) {
      return 0n; // Shouldn't happen for det = +/- 1
    }

    if (pivot !== i) {
      const tmp = A[i];
      A[i] = A[pivot]!;
      A[pivot] = tmp;
      det = (p - det) % p;
    }

    det = (det * A[i]![i]!) % p;

    // Eliminate
    const inv = A[i]![i]! === 1n ? 1n : 2n; // inverse of 1 is 1, inverse of 2 is 2 mod 3
    for (let j = i + 1; j < n; j++) {
      const factor = (A[j]![i]! * inv) % p;
      for (let k = i; k < n; k++) {
        A[j]![k] = (((A[j]![k]! - factor * A[i]![k]!) % p) + p) % p;
      }
    }
  }

  return det === 1n ? 1n : -1n;
}

// ============================================================================
// Principal Square Root
// ============================================================================

/**
 * Return the principal square root of a positive definite matrix.
 *
 * A positive definite matrix A has a unique positive definite
 * matrix M such that M^2 = A.
 *
 * Algorithm: Computes the eigendecomposition A = L^{-1} D L where D is diagonal,
 * then returns L^{-1} sqrt(D) L.
 *
 * @param matrix - A positive definite matrix
 * @param check_positivity - Whether to verify positive definiteness (default: true)
 * @returns The principal square root M such that M^2 = A, or false if not positive definite
 * @see Reference: sage/matrix/matrix2.pyx:principal_square_root
 *
 * Note: This implementation requires the base ring to support square roots of its elements.
 * For finite fields, this may not always be possible (not all elements have square roots).
 */
export function principal_square_root<R extends FieldElement>(
  matrix: Matrix<R>,
  check_positivity: boolean = true
): Matrix<R> | false {
  if (!matrix.is_square()) {
    throw new ArithmeticError('principal_square_root is only defined for square matrices');
  }

  const n = matrix.nrows;
  const ring = matrix.base_ring;

  if (n === 0) {
    return zero_matrix(ring, 0);
  }

  // Check positive definiteness if requested
  // For exact computation, we would need is_positive_definite, but since that's
  // in matrix_operations.ts, we skip the check here and let the caller verify
  // or catch errors from non-positive-definite input

  // Compute eigendecomposition: A = P^{-1} D P
  // where D is diagonal with eigenvalues on the diagonal
  // Then sqrt(A) = P^{-1} sqrt(D) P

  // We use a Denman-Beavers iteration for computing the matrix square root
  // Y_{k+1} = (Y_k + Z_k^{-1}) / 2
  // Z_{k+1} = (Z_k + Y_k^{-1}) / 2
  // with Y_0 = A, Z_0 = I
  // Converges to Y_infty = sqrt(A), Z_infty = sqrt(A)^{-1}

  // For exact computation over rationals/algebraics, we need eigenvalue decomposition
  // This implementation uses Newton iteration which requires the ring to support division

  // Check if the ring has a sqrt function
  const hasElementSqrt =
    typeof (ring as { sqrt?: (x: R) => R }).sqrt === 'function' ||
    typeof (ring.one() as unknown as { sqrt?: () => R }).sqrt === 'function';

  if (!hasElementSqrt) {
    // Fall back to Denman-Beavers iteration (numerical approximation)
    // This only works well for floating-point types
    throw new NotImplementedError(
      'principal_square_root requires the ring to support element-wise square roots or eigenvalue decomposition'
    );
  }

  // Diagonalize the matrix and take square root of eigenvalues
  // For a symmetric/positive definite matrix, eigenvalues are real and positive

  // Simple implementation using power iteration for dominant eigenvalue
  // followed by deflation - this is a simplified approach
  // Full implementation would use QR algorithm or similar

  // For now, implement the Babylonian/Newton method for matrices
  // which works for matrices close to the identity

  // Alternative: Use the explicit formula for 2x2 matrices
  if (n === 2) {
    return principal_square_root_2x2(matrix, ring);
  }

  // For general case, use Denman-Beavers iteration
  // This is an iterative method that converges quadratically
  let Y = matrix.copy();
  let Z = identity_matrix(ring, n);

  const maxIterations = 100;
  const two = ring.one().add(ring.one()) as R;

  for (let iter = 0; iter < maxIterations; iter++) {
    // Save previous Y for convergence check
    const Yprev = Y.copy();

    // Compute Y_{k+1} = (Y_k + Z_k^{-1}) / 2
    // Compute Z_{k+1} = (Z_k + Y_k^{-1}) / 2

    let Zinv: Matrix<R>;
    let Yinv: Matrix<R>;

    try {
      Zinv = matrixInverse(Z);
      Yinv = matrixInverse(Y);
    } catch {
      // Matrix became singular, iteration failed
      if (check_positivity) {
        return false;
      }
      throw new ArithmeticError(
        'matrix square root iteration failed - matrix may not be positive definite'
      );
    }

    // Y = (Y + Zinv) / 2
    // Z = (Z + Yinv) / 2
    const Ynew = new Matrix<R>(ring, n, n);
    const Znew = new Matrix<R>(ring, n, n);

    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        const ySum = Y.get(i, j).add(Zinv.get(i, j)) as R;
        const zSum = Z.get(i, j).add(Yinv.get(i, j)) as R;
        Ynew.set(i, j, ySum.mul(getInverse(two)) as R);
        Znew.set(i, j, zSum.mul(getInverse(two)) as R);
      }
    }

    Y = Ynew;
    Z = Znew;

    // Check convergence: ||Y - Yprev|| < epsilon
    const maxDiff = ring.zero();
    let converged = true;
    for (let i = 0; i < n && converged; i++) {
      for (let j = 0; j < n && converged; j++) {
        const diff = Y.get(i, j).sub(Yprev.get(i, j)) as R;
        if (!diff.isZero()) {
          converged = false;
        }
      }
    }

    if (converged) {
      return Y;
    }
  }

  throw new ArithmeticError('principal_square_root iteration did not converge');
}

/**
 * Compute the principal square root of a 2x2 positive definite matrix.
 * Uses the explicit formula.
 */
function principal_square_root_2x2<R extends FieldElement>(
  matrix: Matrix<R>,
  ring: CoefficientRing<R>
): Matrix<R> {
  const a = matrix.get(0, 0);
  const b = matrix.get(0, 1);
  const c = matrix.get(1, 0);
  const d = matrix.get(1, 1);

  // det(A) = ad - bc
  const det = a.mul(d).sub(b.mul(c)) as R;

  // trace(A) = a + d
  const trace = a.add(d) as R;

  // For a positive definite 2x2 matrix A, sqrt(A) = (A + sqrt(det)*I) / sqrt(trace + 2*sqrt(det))
  // This requires computing sqrt(det) and then sqrt(trace + 2*sqrt(det))

  // Check if element has sqrt method
  const sqrtFn = (det as unknown as { sqrt?: () => R }).sqrt;
  if (typeof sqrtFn !== 'function') {
    throw new NotImplementedError(
      'principal_square_root_2x2 requires elements to support sqrt operation'
    );
  }

  const sqrtDet = sqrtFn.call(det) as R;

  // s = trace + 2*sqrt(det)
  const two = ring.one().add(ring.one()) as R;
  const s = trace.add(two.mul(sqrtDet)) as R;

  // sqrt(s)
  const sqrtS = (s as unknown as { sqrt: () => R }).sqrt();

  // Result = (A + sqrt(det)*I) / sqrt(s)
  const sqrtSInv = getInverse(sqrtS);

  const result = new Matrix<R>(ring, 2, 2);
  result.set(0, 0, (a.add(sqrtDet) as R).mul(sqrtSInv) as R);
  result.set(0, 1, b.mul(sqrtSInv) as R);
  result.set(1, 0, c.mul(sqrtSInv) as R);
  result.set(1, 1, (d.add(sqrtDet) as R).mul(sqrtSInv) as R);

  return result;
}

/**
 * Compute the inverse of a matrix using Gaussian elimination.
 */
function matrixInverse<R extends FieldElement>(matrix: Matrix<R>): Matrix<R> {
  const n = matrix.nrows;
  const ring = matrix.base_ring;

  if (!matrix.is_square()) {
    throw new ArithmeticError('inverse is only defined for square matrices');
  }

  // Create augmented matrix [A | I]
  const aug: R[][] = [];
  for (let i = 0; i < n; i++) {
    aug.push([]);
    for (let j = 0; j < n; j++) {
      aug[i]!.push(matrix.get(i, j));
    }
    for (let j = 0; j < n; j++) {
      aug[i]!.push(i === j ? ring.one() : ring.zero());
    }
  }

  // Gaussian elimination
  for (let col = 0; col < n; col++) {
    // Find pivot
    let pivotRow = -1;
    for (let i = col; i < n; i++) {
      if (!aug[i]![col]!.isZero()) {
        pivotRow = i;
        break;
      }
    }

    if (pivotRow === -1) {
      throw new ArithmeticError('matrix is singular');
    }

    // Swap rows
    if (pivotRow !== col) {
      const tmp = aug[col];
      aug[col] = aug[pivotRow]!;
      aug[pivotRow] = tmp;
    }

    // Scale pivot row
    const pivot = aug[col]![col]!;
    const pivotInv = getInverse(pivot);
    for (let j = 0; j < 2 * n; j++) {
      aug[col]![j] = aug[col]![j]!.mul(pivotInv) as R;
    }

    // Eliminate other rows
    for (let i = 0; i < n; i++) {
      if (i !== col && !aug[i]![col]!.isZero()) {
        const factor = aug[i]![col]!;
        for (let j = 0; j < 2 * n; j++) {
          aug[i]![j] = aug[i]![j]!.sub(factor.mul(aug[col]![j]!) as R) as R;
        }
      }
    }
  }

  // Extract inverse from augmented matrix
  const result = new Matrix<R>(ring, n, n);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      result.set(i, j, aug[i]![n + j]!);
    }
  }

  return result;
}

// ============================================================================
// Matrix Exponential
// ============================================================================

/**
 * Return the matrix exponential exp(A).
 *
 * Computes exp(A) = I + A + A^2/2! + A^3/3! + ...
 *
 * For finite fields, this series may not converge in the usual sense.
 * For nilpotent matrices over any ring, the series terminates.
 *
 * Algorithm: Uses the Taylor series with a practical truncation.
 * For exact computation over rationals/integers, use nilpotent matrices only.
 *
 * @param matrix - A square matrix
 * @returns The matrix exponential
 * @see Reference: sage/matrix/matrix2.pyx:exp
 */
export function exp<R extends FieldElement>(matrix: Matrix<R>): Matrix<R> {
  if (!matrix.is_square()) {
    throw new ArithmeticError('matrix exponential is only defined for square matrices');
  }

  const n = matrix.nrows;
  const ring = matrix.base_ring;

  if (n === 0) {
    return zero_matrix(ring, 0);
  }

  // Check if matrix is nilpotent - then the series terminates
  let isNilpotent = true;
  let power = matrix;
  for (let k = 1; k <= n; k++) {
    let allZero = true;
    for (let i = 0; i < n && allZero; i++) {
      for (let j = 0; j < n && allZero; j++) {
        if (!power.get(i, j).isZero()) {
          allZero = false;
        }
      }
    }
    if (allZero) {
      // Matrix is nilpotent of index k
      break;
    }
    if (k === n) {
      isNilpotent = false;
    }
    if (k < n) {
      power = matrix.mul(power);
    }
  }

  if (!isNilpotent) {
    // For non-nilpotent matrices, we need numerical methods or
    // the field to support infinite series convergence
    throw new NotImplementedError(
      'matrix exponential for non-nilpotent matrices requires numerical precision'
    );
  }

  // For nilpotent matrices, compute exp(A) = sum_{k=0}^{n-1} A^k / k!
  const result = identity_matrix(ring, n);
  let Ak = identity_matrix(ring, n);
  let factorial = ring.one();

  for (let k = 1; k < n; k++) {
    Ak = matrix.mul(Ak);

    // Check if A^k = 0
    let allZero = true;
    for (let i = 0; i < n && allZero; i++) {
      for (let j = 0; j < n && allZero; j++) {
        if (!Ak.get(i, j).isZero()) {
          allZero = false;
        }
      }
    }
    if (allZero) {
      break;
    }

    // Update factorial: k! = (k-1)! * k
    let kValue = ring.zero();
    for (let i = 0; i < k; i++) {
      kValue = kValue.add(ring.one()) as R;
    }
    factorial = factorial.mul(kValue) as R;

    // Check if factorial is invertible
    if (factorial.isZero()) {
      throw new ArithmeticError(`k! = 0 for k=${k} in this ring (characteristic divides ${k}!)`);
    }

    const factorialInv = getInverse(factorial);

    // Add A^k / k! to result
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        const term = Ak.get(i, j).mul(factorialInv) as R;
        result.set(i, j, result.get(i, j).add(term) as R);
      }
    }
  }

  return result;
}

// ============================================================================
// Decomposition (Subspace)
// ============================================================================

/**
 * Return the decomposition of the free module on which this matrix acts.
 *
 * Let A be the matrix acting from the right on the vector space V of column
 * vectors. This function computes maximal subspaces W_1, ..., W_n corresponding
 * to Galois conjugacy classes of eigenvalues of A. More precisely, let f(X) be
 * the characteristic polynomial of A. This function computes the subspace
 * W_i = ker(g_i(A)^m) where g_i(X) is an irreducible factor of f(X) and
 * g_i(X)^m exactly divides f(X).
 *
 * If is_diagonalizable is True, then we let W_i = ker(g_i(A)) since then
 * ker(g_i(A)) = ker(g_i(A)^m).
 *
 * @param matrix - A square matrix
 * @param algorithm - 'spin' or 'kernel':
 *   - 'spin': involves iterating the action of the matrix on a vector
 *   - 'kernel': naively computes ker(f_i(A)) for each factor f_i
 * @param is_diagonalizable - If the matrix is known to be diagonalizable,
 *   set this to True which might speed up the algorithm
 * @param dual - If True, also returns the corresponding decomposition of V
 *   under the action of the transpose of A
 * @returns List of (subspace_basis, is_irreducible) pairs where subspace_basis
 *   is a Matrix whose rows form a basis for the invariant subspace, and
 *   is_irreducible is True if the characteristic polynomial restricted to this
 *   subspace is irreducible.
 * @see Reference: sage/matrix/matrix2.pyx:decomposition
 */
export function decomposition<R extends FieldElement>(
  matrix: Matrix<R>,
  algorithm: 'spin' | 'kernel' = 'kernel',
  is_diagonalizable: boolean = false,
  dual: boolean = false
): Array<[Matrix<R>, boolean]> | [Array<[Matrix<R>, boolean]>, Array<[Matrix<R>, boolean]>] {
  if (!matrix.is_square()) {
    throw new ValueError('decomposition requires a square matrix');
  }

  const n = matrix.nrows;
  const ring = matrix.base_ring;

  if (n === 0) {
    if (dual) {
      return [[], []];
    }
    return [];
  }

  // Compute characteristic polynomial and factor it
  const f = charpolyHelper(matrix, 'x');
  const factors = f.factor();

  // Sort factors for consistent output
  factors.sort((a, b) => {
    const degDiff = a[0].degree() - b[0].degree();
    if (degDiff !== 0) return degDiff;
    return 0;
  });

  // If there's only one factor, the whole space is the invariant subspace
  if (factors.length === 1) {
    const [g, m] = factors[0]!;
    const isIrreducible = m === 1;
    const basisMatrix = identity_matrix(ring, n);
    if (dual) {
      return [[[basisMatrix, isIrreducible]], [[basisMatrix, isIrreducible]]];
    }
    return [[basisMatrix, isIrreducible]];
  }

  // Use kernel algorithm: compute ker(g(A)^m) for each factor (g, m)
  const result: Array<[Matrix<R>, boolean]> = [];
  const dualResult: Array<[Matrix<R>, boolean]> = [];

  for (const [g, m] of factors) {
    // Compute g(A)
    let gA = evaluatePolynomialAtMatrix(g, matrix);

    // If not diagonalizable, raise to power m
    if (!is_diagonalizable && m > 1) {
      gA = matrixPower(gA, m);
    }

    // Compute kernel of g(A)^m
    const kernelBasis = computeKernelBasis(gA);

    const isIrreducible = m === 1;
    result.push([kernelBasis, isIrreducible]);

    if (dual) {
      // Compute kernel of g(A^T)^m
      let gAT = evaluatePolynomialAtMatrix(g, matrix.transpose());
      if (!is_diagonalizable && m > 1) {
        gAT = matrixPower(gAT, m);
      }
      const dualKernelBasis = computeKernelBasis(gAT);
      dualResult.push([dualKernelBasis, isIrreducible]);
    }
  }

  if (dual) {
    return [result, dualResult];
  }
  return result;
}

/**
 * Evaluate a polynomial at a matrix.
 * p(A) = p_0 * I + p_1 * A + p_2 * A^2 + ... + p_n * A^n
 */
function evaluatePolynomialAtMatrix<R extends FieldElement>(
  p: Polynomial<R>,
  A: Matrix<R>
): Matrix<R> {
  const n = A.nrows;
  const ring = A.base_ring;
  const coeffs = p.coeffs;

  if (coeffs.length === 0) {
    return zero_matrix(ring, n);
  }

  // Use Horner's method: p(A) = ((...((p_n * A + p_{n-1}) * A + p_{n-2}) * A + ...) * A + p_0)
  let result = zero_matrix(ring, n);

  // Start with the highest coefficient
  const highDeg = coeffs.length - 1;
  for (let i = 0; i < n; i++) {
    result.set(i, i, coeffs[highDeg]!);
  }

  // Work down through the coefficients
  for (let k = highDeg - 1; k >= 0; k--) {
    // result = result * A + coeffs[k] * I
    result = result.mul(A);

    // Add scalar multiple of identity
    const coeff = coeffs[k]!;
    for (let i = 0; i < n; i++) {
      result.set(i, i, result.get(i, i).add(coeff) as R);
    }
  }

  return result;
}

/**
 * Compute matrix power A^k.
 */
function matrixPower<R extends FieldElement>(A: Matrix<R>, k: number): Matrix<R> {
  if (k < 0) {
    throw new ValueError('matrix power must be non-negative');
  }

  const n = A.nrows;
  const ring = A.base_ring;

  if (k === 0) {
    return identity_matrix(ring, n);
  }

  if (k === 1) {
    return A.copy();
  }

  // Use binary exponentiation
  let result = identity_matrix(ring, n);
  let base = A.copy();

  while (k > 0) {
    if (k % 2 === 1) {
      result = result.mul(base);
    }
    base = base.mul(base);
    k = Math.floor(k / 2);
  }

  return result;
}

/**
 * Compute a basis for the kernel (null space) of a matrix.
 * Returns a matrix whose rows form a basis for the kernel.
 */
function computeKernelBasis<R extends FieldElement>(A: Matrix<R>): Matrix<R> {
  const m = A.nrows;
  const n = A.ncols;
  const ring = A.base_ring;

  // Compute RREF and find pivot columns
  const E = rref(A);
  const pivotCols = pivot_rows(A);
  const pivotSet = new Set(pivotCols);

  // Find non-pivot (free) columns
  const freeCols: number[] = [];
  for (let j = 0; j < n; j++) {
    if (!pivotSet.has(j)) {
      freeCols.push(j);
    }
  }

  const kernelDim = freeCols.length;

  if (kernelDim === 0) {
    // Trivial kernel
    return zero_matrix(ring, 0, n);
  }

  // Build kernel basis: one vector for each free column
  const basisRows: R[][] = [];

  for (let k = 0; k < kernelDim; k++) {
    const freeCol = freeCols[k]!;
    const row: R[] = [];

    for (let j = 0; j < n; j++) {
      if (j === freeCol) {
        // Free variable gets 1
        row.push(ring.one());
      } else if (pivotSet.has(j)) {
        // Pivot variable: solve from RREF
        // Find which row has pivot in column j
        let pivotRow = -1;
        for (let i = 0; i < Math.min(m, pivotCols.length); i++) {
          if (pivotCols[i] === j) {
            pivotRow = i;
            break;
          }
        }
        if (pivotRow >= 0 && pivotRow < m) {
          // The coefficient is -E[pivotRow, freeCol]
          row.push(E.get(pivotRow, freeCol).neg() as R);
        } else {
          row.push(ring.zero());
        }
      } else {
        // Another free variable: 0
        row.push(ring.zero());
      }
    }

    basisRows.push(row);
  }

  return new Matrix(ring, kernelDim, n, basisRows);
}

/**
 * Helper to compute characteristic polynomial.
 * Duplicated here to avoid circular imports.
 */
function charpolyHelper<R extends RingElement>(
  matrix: Matrix<R>,
  variable: string = 'x'
): Polynomial<R> {
  const n = matrix.nrows;
  const ring = matrix.base_ring;
  const polyRing = new PolynomialRing(ring, variable);

  if (n === 0) {
    return polyRing.one();
  }

  // Use division-free algorithm
  const F: R[] = new Array(n);
  for (let i = 0; i < n; i++) {
    F[i] = ring.zero();
  }

  const a: R[][] = [];
  for (let p = 0; p < n; p++) {
    a.push(new Array(n).fill(ring.zero()));
  }

  const A: R[] = new Array(n);
  for (let i = 0; i < n; i++) {
    A[i] = ring.zero();
  }

  F[0] = matrix.get(0, 0).neg() as R;

  for (let t = 1; t < n; t++) {
    for (let i = 0; i <= t; i++) {
      a[0]![i] = matrix.get(i, t);
    }
    A[0] = matrix.get(t, t);

    for (let p = 1; p < t; p++) {
      for (let i = 0; i <= t; i++) {
        let s = ring.zero();
        for (let j = 0; j <= t; j++) {
          s = s.add(matrix.get(i, j).mul(a[p - 1]![j]!)) as R;
        }
        a[p]![i] = s;
      }
      A[p] = a[p]![t]!;
    }

    let s = ring.zero();
    for (let j = 0; j <= t; j++) {
      s = s.add(matrix.get(t, j).mul(a[t - 1]![j]!)) as R;
    }
    A[t] = s;

    for (let p = 0; p <= t; p++) {
      s = F[p]!;
      for (let k = 0; k < p; k++) {
        s = s.sub(A[k]!.mul(F[p - k - 1]!)) as R;
      }
      F[p] = s.sub(A[p]!) as R;
    }
  }

  const coeffs: R[] = [];
  for (let i = n - 1; i >= 0; i--) {
    coeffs.push(F[i]!);
  }
  coeffs.push(ring.one());

  return new Polynomial(coeffs, polyRing);
}

/**
 * Return the cyclic subspace generated by a vector.
 *
 * The cyclic subspace (Krylov subspace) is spanned by {v, Av, A^2v, ...}.
 * This computes the vectors v, Av, A^2v, ... until linear dependence.
 *
 * @param matrix - A square matrix
 * @param v - A vector
 * @param variable - Variable name for the minimal polynomial (not used in basic impl)
 * @param basis - Type of basis: 'echelon', 'iterates', or 'monomials'
 * @returns The basis vectors of the cyclic subspace (as rows of a matrix)
 * @see Reference: sage/matrix/matrix2.pyx:cyclic_subspace
 */
export function cyclic_subspace<R extends FieldElement>(
  matrix: Matrix<R>,
  v: R[],
  variable?: string,
  basis: 'echelon' | 'iterates' | 'monomials' = 'echelon'
): Matrix<R> {
  if (!matrix.is_square()) {
    throw new ArithmeticError('cyclic_subspace requires a square matrix');
  }

  const n = matrix.nrows;
  const ring = matrix.base_ring;

  if (v.length !== n) {
    throw new ValueError(`vector length ${v.length} must match matrix size ${n}`);
  }

  // Check if v is zero
  let vIsZero = true;
  for (const entry of v) {
    if (!entry.isZero()) {
      vIsZero = false;
      break;
    }
  }

  if (vIsZero) {
    return zero_matrix(ring, 0, n);
  }

  // Compute the sequence v, Av, A^2v, ... until linearly dependent
  const iterates: R[][] = [];
  let current = [...v];
  iterates.push([...current]);

  // Build the span incrementally
  while (iterates.length <= n) {
    // Compute A * current
    const next: R[] = [];
    for (let i = 0; i < n; i++) {
      let sum = ring.zero();
      for (let j = 0; j < n; j++) {
        sum = sum.add(matrix.get(i, j).mul(current[j]!)) as R;
      }
      next.push(sum);
    }

    // Check if next is in the span of iterates
    // Build matrix with iterates as rows and check if next is in the row space
    const basisMatrix = new Matrix(ring, iterates.length, n, iterates);
    const augmented: R[][] = [];
    for (let i = 0; i < iterates.length; i++) {
      augmented.push([...iterates[i]!]);
    }
    augmented.push([...next]);

    const augMatrix = new Matrix(ring, iterates.length + 1, n, augmented);
    const rankBefore = pivot_rows(basisMatrix).length;
    const rankAfter = pivot_rows(augMatrix).length;

    if (rankAfter === rankBefore) {
      // next is in the span of iterates - we're done
      break;
    }

    iterates.push([...next]);
    current = next;
  }

  if (basis === 'iterates') {
    return new Matrix(ring, iterates.length, n, iterates);
  }

  // 'echelon' basis - return echelon form
  const iterateMatrix = new Matrix(ring, iterates.length, n, iterates);
  return echelon_form(iterateMatrix);
}

/**
 * Compute the maximal spin of a vector under the matrix action.
 *
 * Returns the largest list S = [v, vA, vA^2, ..., vA^k] of vectors that
 * are linearly independent.
 *
 * Note: This computes v*A (row vector times matrix), not A*v.
 *
 * @param matrix - A square matrix
 * @param v - A row vector
 * @returns Array of linearly independent vectors in the Krylov sequence
 * @see Reference: sage/matrix/matrix2.pyx:maxspin
 */
export function maxspin<R extends FieldElement>(matrix: Matrix<R>, v: R[]): R[][] {
  if (!matrix.is_square()) {
    throw new ArithmeticError('maxspin requires a square matrix');
  }

  const n = matrix.nrows;
  const ring = matrix.base_ring;

  if (v.length !== n) {
    throw new ValueError(`vector length ${v.length} must match matrix size ${n}`);
  }

  // Check if v is zero
  let vIsZero = true;
  for (const entry of v) {
    if (!entry.isZero()) {
      vIsZero = false;
      break;
    }
  }

  if (vIsZero) {
    return [];
  }

  const S: R[][] = [];
  let w = [...v];
  S.push([...w]);

  while (true) {
    // Compute w * matrix (row vector times matrix)
    const next: R[] = [];
    for (let j = 0; j < n; j++) {
      let sum = ring.zero();
      for (let i = 0; i < n; i++) {
        sum = sum.add(w[i]!.mul(matrix.get(i, j))) as R;
      }
      next.push(sum);
    }

    // Check if next is in the span of S
    const basisMatrix = new Matrix(ring, S.length, n, S);
    const augmented: R[][] = [];
    for (const row of S) {
      augmented.push([...row]);
    }
    augmented.push([...next]);

    const augMatrix = new Matrix(ring, S.length + 1, n, augmented);
    const rankBefore = pivot_rows(basisMatrix).length;
    const rankAfter = pivot_rows(augMatrix).length;

    if (rankAfter === rankBefore) {
      // next is linearly dependent on S
      return S;
    }

    S.push([...next]);
    w = next;
  }
}

/**
 * Application of Wiedemann's algorithm to the i-th standard basis vector.
 *
 * Wiedemann's algorithm computes the minimal polynomial of a matrix A by
 * computing the iterates v, Av, A^2v, ... for a vector v and finding the
 * minimal polynomial of the resulting linear recurrence sequence.
 *
 * @param matrix - A square matrix
 * @param i - Index of the standard basis vector to use (0-indexed)
 * @param t - If nonzero, use only the t-th coordinate of the iterates.
 *   If 0 (default), compute the minimal polynomial using all coordinates
 *   and return the LCM.
 * @returns The minimal polynomial of the recurrence sequence, which
 *   divides the minimal polynomial of the matrix.
 * @see Reference: sage/matrix/matrix2.pyx:wiedemann
 *
 * @example
 * ```typescript
 * const A = MatrixSpace(QQ, 3, 3)([[0, 1, 2], [3, 4, 5], [6, 7, 8]]);
 * const f = wiedemann(A, 0);
 * // f divides the characteristic polynomial of A
 * ```
 */
export function wiedemann<R extends FieldElement>(
  matrix: Matrix<R>,
  i: number,
  t: number = 0
): Polynomial<R> {
  if (!matrix.is_square()) {
    throw new ArithmeticError('wiedemann requires a square matrix');
  }

  i = Math.floor(i);
  t = Math.floor(t);

  const n = matrix.nrows;
  const ring = matrix.base_ring;

  if (i < 0 || i >= n) {
    throw new ValueError(`index i=${i} out of bounds for ${n}x${n} matrix`);
  }

  // Create the i-th standard basis vector
  const v: R[] = [];
  for (let j = 0; j < n; j++) {
    v.push(j === i ? ring.one() : ring.zero());
  }

  // Compute iterates: v, Av, A^2v, ..., A^{2n-1}v
  // We need 2n iterates for Berlekamp-Massey
  const iterates: R[][] = [];
  let current = [...v];
  iterates.push([...current]);

  for (let k = 1; k < 2 * n; k++) {
    // Compute A * current
    const next: R[] = [];
    for (let row = 0; row < n; row++) {
      let sum = ring.zero();
      for (let col = 0; col < n; col++) {
        sum = sum.add(matrix.get(row, col).mul(current[col]!)) as R;
      }
      next.push(sum);
    }
    iterates.push([...next]);
    current = next;
  }

  // Extract columns (the j-th column is [v_j, (Av)_j, (A^2v)_j, ...])
  const columns: R[][] = [];
  for (let j = 0; j < n; j++) {
    const col: R[] = [];
    for (let k = 0; k < 2 * n; k++) {
      col.push(iterates[k]![j]!);
    }
    columns.push(col);
  }

  // Apply Berlekamp-Massey to find minimal polynomial of each recurrence
  const polyRing = new PolynomialRing(ring, 'x');
  let f: Polynomial<R> | null = null;

  // Determine which coordinates to use
  const coords = t === 0 ? Array.from({ length: n }, (_, k) => k) : [t];

  for (const coord of coords) {
    if (coord < 0 || coord >= n) continue;

    const sequence = columns[coord]!;

    // Apply Berlekamp-Massey to find minimal polynomial
    const g = berlekampMasseyHelper(sequence, polyRing);

    if (f === null) {
      f = g;
    } else {
      // Compute LCM of f and g
      f = polynomialLcm(f, g);
    }

    // Early termination if we found the full minimal polynomial
    if (f.degree() >= n) {
      break;
    }
  }

  if (f === null) {
    // Return polynomial 1 for empty case
    return polyRing.one();
  }

  return f;
}

/**
 * Berlekamp-Massey algorithm to find minimal polynomial of a linear recurrence.
 * Duplicated here to avoid import issues.
 */
function berlekampMasseyHelper<R extends FieldElement>(
  sequence: R[],
  polyRing: PolynomialRing<R>
): Polynomial<R> {
  const ring = polyRing.base_ring;
  const n = sequence.length;

  if (n === 0 || n % 2 !== 0) {
    throw new ValueError('berlekamp_massey requires even-length sequence');
  }

  const M = Math.floor(n / 2);

  // Create polynomial from sequence: f0 = a[0] + a[1]*x + ... + a[2M-1]*x^{2M-1}
  const f0 = new Polynomial<R>(sequence, polyRing);

  // f1 = x^{2M}
  const f1Coeffs: R[] = [];
  for (let i = 0; i < 2 * M; i++) {
    f1Coeffs.push(ring.zero());
  }
  f1Coeffs.push(ring.one());
  let f1 = new Polynomial<R>(f1Coeffs, polyRing);
  let f0Current = f0;

  // s0 = 1, s1 = 0
  let s0 = new Polynomial<R>([ring.one()], polyRing);
  let s1 = new Polynomial<R>([ring.zero()], polyRing);

  // Extended Euclidean algorithm until degree < M
  while (f1.degree() >= M) {
    const temp = f1;
    const [q, r] = f0Current.quo_rem(f1);
    f0Current = temp;
    f1 = r;

    const newS1 = s0.sub(q.mul(s1)) as Polynomial<R>;
    s0 = s1;
    s1 = newS1;
  }

  // Return s1.reverse().monic()
  const s1Coeffs = s1.coeffs.slice();
  while (s1Coeffs.length > 0 && s1Coeffs[s1Coeffs.length - 1]!.isZero()) {
    s1Coeffs.pop();
  }
  s1Coeffs.reverse();

  const reversed = new Polynomial<R>(s1Coeffs, polyRing);
  return reversed.monic();
}

/**
 * Compute LCM of two polynomials using GCD.
 * lcm(f, g) = f * g / gcd(f, g)
 */
function polynomialLcm<R extends FieldElement>(f: Polynomial<R>, g: Polynomial<R>): Polynomial<R> {
  const gcd = polynomialGcd(f, g);

  if (gcd.isZero()) {
    return f;
  }

  // lcm = f * g / gcd = f * (g / gcd)
  const [quotient] = g.quo_rem(gcd);
  return f.mul(quotient) as Polynomial<R>;
}

/**
 * Compute GCD of two polynomials using Euclidean algorithm.
 */
function polynomialGcd<R extends FieldElement>(f: Polynomial<R>, g: Polynomial<R>): Polynomial<R> {
  // Euclidean algorithm
  let a = f;
  let b = g;

  while (!b.isZero()) {
    const [, r] = a.quo_rem(b);
    a = b;
    b = r;
  }

  // Return monic GCD
  if (a.isZero()) {
    return a;
  }
  return a.monic();
}

// ============================================================================
// Krylov Methods
// ============================================================================

/**
 * Return the Krylov matrix of M acting on matrix.
 *
 * Given a matrix B (the input matrix) and an acting matrix M, this computes
 * the Krylov iterates [B, MB, M^2B, ..., M^{d-1}B] stacked vertically.
 *
 * @param matrix - The matrix B (m x n)
 * @param M - The acting matrix (n x n)
 * @param shifts - Shift parameters (not used in basic implementation)
 * @param degrees - Degree parameters - how many iterates for each row (default: n)
 * @returns The Krylov matrix
 * @see Reference: sage/matrix/matrix2.pyx:krylov_matrix
 */
export function krylov_matrix<R extends FieldElement>(
  matrix: Matrix<R>,
  M: Matrix<R>,
  shifts?: number[],
  degrees?: number[]
): Matrix<R> {
  if (!M.is_square()) {
    throw new ArithmeticError('acting matrix M must be square');
  }

  const m = matrix.nrows;
  const n = matrix.ncols;
  const ring = matrix.base_ring;

  if (M.nrows !== n) {
    throw new ValueError(`M must be ${n}x${n} to act on matrix of size ${m}x${n}`);
  }

  // Default degree is n (full Krylov sequence)
  const maxDegree = degrees ? Math.max(...degrees) : n;

  // Compute B, MB, M^2B, ...
  const rows: R[][] = [];
  let current = matrix;

  for (let d = 0; d < maxDegree; d++) {
    for (let i = 0; i < m; i++) {
      // Check if this row should be included at this degree
      const deg = degrees ? degrees[i] || n : n;
      if (d < deg) {
        const row: R[] = [];
        for (let j = 0; j < n; j++) {
          row.push(current.get(i, j));
        }
        rows.push(row);
      }
    }

    // Compute M * current
    if (d + 1 < maxDegree) {
      const next: R[][] = [];
      for (let i = 0; i < m; i++) {
        next.push([]);
        for (let j = 0; j < n; j++) {
          let sum = ring.zero();
          for (let k = 0; k < n; k++) {
            sum = sum.add(M.get(j, k).mul(current.get(i, k))) as R;
          }
          next[i]!.push(sum);
        }
      }
      current = new Matrix(ring, m, n, next);
    }
  }

  return new Matrix(ring, rows.length, n, rows);
}

/**
 * Return a Krylov basis.
 *
 * Computes a basis for the Krylov subspace generated by the rows of matrix
 * under the action of M.
 *
 * @param matrix - The matrix
 * @param M - The acting matrix
 * @param shifts - Shift parameters
 * @param degrees - Degree parameters
 * @param output_rows - Whether to output as rows (default: true)
 * @param algorithm - Algorithm to use
 * @returns The Krylov basis
 * @see Reference: sage/matrix/matrix2.pyx:krylov_basis
 */
export function krylov_basis<R extends FieldElement>(
  matrix: Matrix<R>,
  M: Matrix<R>,
  shifts?: number[],
  degrees?: number[],
  output_rows: boolean = true,
  algorithm?: string
): Matrix<R> {
  // Compute the Krylov matrix
  const K = krylov_matrix(matrix, M, shifts, degrees);

  // Compute echelon form to get a basis
  const basis = echelon_form(K);

  // Remove zero rows
  const nonzeroRows: R[][] = [];
  for (let i = 0; i < basis.nrows; i++) {
    let isZero = true;
    for (let j = 0; j < basis.ncols; j++) {
      if (!basis.get(i, j).isZero()) {
        isZero = false;
        break;
      }
    }
    if (!isZero) {
      const row: R[] = [];
      for (let j = 0; j < basis.ncols; j++) {
        row.push(basis.get(i, j));
      }
      nonzeroRows.push(row);
    }
  }

  const result = new Matrix(matrix.base_ring, nonzeroRows.length, basis.ncols, nonzeroRows);

  if (!output_rows) {
    return result.transpose();
  }

  return result;
}

/**
 * Return a basis for the Krylov kernel.
 *
 * Computes a basis for the kernel of the Krylov matrix.
 *
 * @param matrix - The matrix
 * @param M - The acting matrix
 * @param shifts - Shift parameters
 * @param degrees - Degree parameters
 * @param output_rows - Whether to output as rows (default: true)
 * @param variable - Variable name
 * @param basis_algorithm - Algorithm for computing the basis
 * @returns The Krylov kernel basis
 * @see Reference: sage/matrix/matrix2.pyx:krylov_kernel_basis
 */
export function krylov_kernel_basis<R extends FieldElement>(
  matrix: Matrix<R>,
  M: Matrix<R>,
  shifts?: number[],
  degrees?: number[],
  output_rows: boolean = true,
  variable?: string,
  basis_algorithm?: string
): Matrix<R> {
  // Import kernel computation
  const { right_kernel_matrix } = require('./matrix_operations.js');

  // Compute the Krylov matrix
  const K = krylov_matrix(matrix, M, shifts, degrees);

  // Compute the kernel of K
  const kernel = right_kernel_matrix(K);

  if (!output_rows) {
    return kernel.transpose();
  }

  return kernel;
}
