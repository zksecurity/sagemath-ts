/**
 * @module sage/matrix/matrix_decompositions
 * @description Matrix decompositions (LU, QR, SVD, Cholesky, Smith, Hermite, etc.)
 *
 * Port of: sage/matrix/matrix2.pyx
 */

import { ArithmeticError, NotImplementedError, ValueError, ZeroDivisionError } from '../errors.js';
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
 * Reduce ``M`` to reduced row echelon form in place, optionally applying the
 * same row operations to ``T``.
 *
 * Over a field SageMath's echelon form *is* the reduced row echelon form
 * (``matrix2.pyx:echelon_form``, doctest over GF(19) gives ``[1 0 18 / 0 1 2]``),
 * so pivots are scaled to one and every entry above a pivot is cleared.
 *
 * @param M - matrix to reduce (mutated)
 * @param T - optional matrix receiving the same row operations (mutated)
 * @returns the list of pivot column indices
 */
function _echelonize_in_place<R extends FieldElement>(M: Matrix<R>, T?: Matrix<R>): number[] {
  const m = M.nrows;
  const n = M.ncols;

  const pivotCols: number[] = [];
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
      if (T) {
        for (let j = 0; j < T.ncols; j++) {
          const tmp = T.get(pivotRow, j);
          T.set(pivotRow, j, T.get(found, j));
          T.set(found, j, tmp);
        }
      }
    }

    // Scale the pivot row so that the pivot equals one
    const pivotInv = getInverse(M.get(pivotRow, col));
    for (let j = 0; j < n; j++) {
      M.set(pivotRow, j, M.get(pivotRow, j).mul(pivotInv) as R);
    }
    if (T) {
      for (let j = 0; j < T.ncols; j++) {
        T.set(pivotRow, j, T.get(pivotRow, j).mul(pivotInv) as R);
      }
    }

    // Clear the rest of the column, both below *and* above the pivot
    for (let i = 0; i < m; i++) {
      if (i === pivotRow) {
        continue;
      }
      const entry = M.get(i, col);
      if (entry.isZero()) {
        continue;
      }
      for (let j = 0; j < n; j++) {
        const val = M.get(i, j).sub(entry.mul(M.get(pivotRow, j)) as R) as R;
        M.set(i, j, val);
      }
      if (T) {
        for (let j = 0; j < T.ncols; j++) {
          const val = T.get(i, j).sub(entry.mul(T.get(pivotRow, j)) as R) as R;
          T.set(i, j, val);
        }
      }
    }

    pivotCols.push(col);
    pivotRow++;
  }

  return pivotCols;
}

/**
 * Transform the matrix into (reduced row) echelon form in place.
 *
 * Over a field SageMath's echelon form is the *reduced* row echelon form.
 *
 * @param matrix - The matrix to echelonize (mutated)
 * @param algorithm - Algorithm to use ('default', 'classical', etc.)
 * @param cutoff - Cutoff for Strassen algorithm (not used)
 * @param transformation - Whether to return the transformation matrix
 * @returns The transformation matrix T with T*A = E if requested
 * @see Reference: sage/matrix/matrix2.pyx:echelonize
 */
export function echelonize<R extends FieldElement>(
  matrix: Matrix<R>,
  algorithm?: 'default' | 'classical' | 'strassen' | 'partial_pivoting' | 'scaled_partial_pivoting',
  cutoff?: number,
  transformation?: boolean
): Matrix<R> | undefined {
  const ring = matrix.base_ring;

  // If transformation is requested, start with identity matrix
  let T: Matrix<R> | undefined;
  if (transformation) {
    T = identity_matrix(ring, matrix.nrows);
  }

  _echelonize_in_place(matrix, T);

  if (transformation) {
    return T;
  }
}

/**
 * Return the echelon form of the matrix.
 *
 * Over a field this is the *reduced* row echelon form: every pivot is one and
 * is the only non-zero entry of its column (``matrix2.pyx:echelon_form``).
 *
 * This does not change the matrix itself.
 *
 * @param matrix - The matrix
 * @param algorithm - Algorithm to use
 * @param cutoff - Cutoff for Strassen algorithm
 * @returns A new matrix in reduced row echelon form
 * @see Reference: sage/matrix/matrix2.pyx:echelon_form
 */
export function echelon_form<R extends FieldElement>(
  matrix: Matrix<R>,
  algorithm?: 'default' | 'classical' | 'strassen',
  cutoff?: number
): Matrix<R> {
  const M = matrix.copy();
  _echelonize_in_place(M);
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
  return echelon_form(matrix);
}

/**
 * Return the extended echelon form of the matrix.
 *
 * For an m x n matrix A, computes [E | T] where E is the RREF of A and T is
 * the transformation matrix (T * original = E).
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
 * These are the indices of the columns of the echelon form that contain a
 * leading one.
 *
 * @param matrix - The matrix
 * @returns Array of pivot column indices
 * @see Reference: sage/matrix/matrix2.pyx:pivots
 */
export function pivots<R extends FieldElement>(matrix: Matrix<R>): number[] {
  const M = matrix.copy();
  return _echelonize_in_place(M);
}

/**
 * Return the pivot row positions of the matrix.
 *
 * These are a topmost subset of the rows that span the row space and are
 * linearly independent.  SageMath computes this as ``self.transpose().pivots()``
 * (``matrix2.pyx:1014``); e.g. ``matrix(QQ,3,3,[0,0,0,1,2,3,2,4,6]).pivot_rows()``
 * is ``(1,)``.
 *
 * @param matrix - The matrix
 * @returns Array of pivot row indices
 * @see Reference: sage/matrix/matrix2.pyx:pivot_rows
 */
export function pivot_rows<R extends FieldElement>(matrix: Matrix<R>): number[] {
  return pivots(matrix.transpose());
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

  // Build P, L, U matrices.
  // Sage builds P as Permutation([perm[i]+1 ...]).to_matrix(), which places the
  // one at row perm[i], column i (sage/combinat/permutation.py:1373), so that
  // P*L*U == A.
  const P = zero_matrix(ring, m, m);
  for (let i = 0; i < m; i++) {
    P.set(perm[i]!, i, ring.one());
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
 * @param full - If true (default, as in Sage), return the full m x m Q;
 *   if false, the reduced form
 * @returns Pair (Q, R) where A = Q * R
 * @see Reference: sage/matrix/matrix2.pyx:QR
 * @see Deviation: Sage's Q is unitary (columns scaled by 1/sqrt(<v,v>)); we do
 *   not scale, so our Q has orthogonal but not orthonormal columns and our R
 *   differs from Sage's by the corresponding diagonal factor.  Scaling requires
 *   square roots in the base field, which generic rings do not provide.
 */
export function QR<R extends FieldElement>(
  matrix: Matrix<R>,
  full: boolean = true
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
 * Perform Gram-Schmidt orthogonalization on the **rows** of the matrix.
 *
 * Returns a pair (G, M) such that, writing A for ``matrix``:
 * - ``A = M * G``
 * - the rows of G are an orthogonal set spanning the row space of A
 * - ``G * G.transpose()`` is diagonal
 * - M is full rank with zeros above the diagonal
 *
 * Zero vectors arising from linear dependence are dropped, so G has exactly
 * ``rank(A)`` rows.
 *
 * Sage implements this as ``self.transpose()._gram_schmidt_noscale()`` followed
 * by transposing both results (``matrix2.pyx:11806``).
 *
 * @param matrix - The matrix whose rows are to be orthogonalized
 * @param orthonormal - Whether to normalize the vectors (default: false);
 *   requires square roots in the base field and is not supported here
 * @returns Pair (G, M) with A = M * G
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

  const [Q, Rm] = gram_schmidt_noscale(matrix.transpose());
  return [Q.transpose(), Rm.transpose()];
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
 * Return the conjugate of a ring element, or the element itself when the ring
 * has no conjugation (as is the case for real/finite fields).
 */
function _conjugate<R extends RingElement>(x: R): R {
  const c = (x as unknown as { conjugate?: () => R }).conjugate;
  if (typeof c === 'function') {
    return c.call(x);
  }
  return x;
}

/**
 * Utility function that decomposes a symmetric or Hermitian matrix into a unit
 * lower-triangular matrix ``L`` and the list ``d`` of diagonal entries, so that
 * ``A == L*diagonal_matrix(d)*L.transpose()`` (resp. ``L.conjugate_transpose()``).
 *
 * Returns ``false`` together with the size of the singular leading principal
 * submatrix when a zero pivot is met, exactly as Sage's
 * ``_indefinite_factorization`` does.
 *
 * @see Reference: sage/matrix/matrix2.pyx:_indefinite_factorization
 */
function _indefinite_factorization<R extends FieldElement>(
  matrix: Matrix<R>,
  algorithm: 'symmetric' | 'hermitian',
  check: boolean
): [Matrix<R>, R[]] | [false, number] {
  const m = matrix.nrows;
  const ring = matrix.base_ring;

  if (!matrix.is_square()) {
    throw new ValueError(`matrix must be square, not ${matrix.nrows} x ${matrix.ncols}`);
  }
  if (algorithm !== 'symmetric' && algorithm !== 'hermitian') {
    throw new ValueError(`'algorithm' must be 'symmetric' or 'hermitian', not ${algorithm}`);
  }

  const conjugate = algorithm === 'hermitian';

  if (check) {
    for (let i = 0; i < m; i++) {
      for (let j = i + 1; j < m; j++) {
        const upper = conjugate ? _conjugate(matrix.get(j, i)) : matrix.get(j, i);
        if (!matrix.get(i, j).eq(upper)) {
          if (conjugate) {
            throw new ValueError('matrix is not hermitian');
          }
          throw new ValueError("matrix is not symmetric (maybe try the 'hermitian' keyword)");
        }
      }
    }
  }

  const L = matrix.copy();
  const d: R[] = [];
  const dInv: R[] = [];

  for (let i = 0; i < m; i++) {
    for (let j = 0; j <= i; j++) {
      let t = L.get(i, j);
      for (let k = 0; k < j; k++) {
        const factor = conjugate ? _conjugate(L.get(j, k)) : L.get(j, k);
        t = t.sub(L.get(k, i).mul(factor) as R) as R;
      }
      if (i === j) {
        if (t.isZero()) {
          return [false, i + 1];
        }
        d.push(t);
        dInv.push(getInverse(t));
        L.set(i, i, ring.one());
      } else {
        L.set(j, i, t);
        L.set(i, j, dInv[j]!.mul(t) as R);
      }
    }
  }

  // Triangularize output matrix
  for (let i = 0; i < m; i++) {
    for (let j = i + 1; j < m; j++) {
      L.set(i, j, ring.zero());
    }
  }

  return [L, d];
}

/**
 * Decompose a symmetric or Hermitian matrix into a unit lower-triangular matrix
 * and a diagonal.
 *
 * Computes L and d such that ``A = L*diagonal_matrix(d)*L.transpose()`` (or the
 * conjugate transpose when ``algorithm='hermitian'``).
 *
 * Unlike Cholesky, this works for indefinite matrices, but it fails — and Sage
 * raises a :class:`ValueError` — as soon as a leading principal submatrix is
 * singular.
 *
 * @param matrix - A symmetric or Hermitian matrix
 * @param algorithm - 'symmetric' or 'hermitian' (default: 'symmetric')
 * @param check - Whether to check symmetry/Hermitian property (default: true)
 * @returns Pair (L, d) where d is the diagonal of D as a vector
 * @see Reference: sage/matrix/matrix2.pyx:indefinite_factorization
 * @see Deviation: Sage returns ``d`` as a Sage vector; we return a plain array,
 *   the port's representation of a vector.
 */
export function indefinite_factorization<R extends FieldElement>(
  matrix: Matrix<R>,
  algorithm: 'symmetric' | 'hermitian' = 'symmetric',
  check: boolean = true
): [Matrix<R>, R[]] {
  const result = _indefinite_factorization(matrix, algorithm, check);
  if (result[0] === false) {
    const k = result[1];
    throw new ValueError(
      `${k}x${k} leading principal submatrix is singular, so cannot create indefinite factorization`
    );
  }
  return result as [Matrix<R>, R[]];
}

/**
 * Convert ``|x|`` into a JavaScript number for the Bunch-Kaufman pivot
 * comparisons, or return ``undefined`` when the base ring carries no absolute
 * value (e.g. a finite field).
 *
 * Sage itself performs these comparisons in C ``double`` arithmetic
 * (``matrix2.pyx:15191``), so no exactness is lost by doing the same here: the
 * comparisons only select which rows/columns get swapped.
 */
function _bk_abs<R extends RingElement>(x: R): number | undefined {
  const absFn = (x as unknown as { abs?: () => unknown }).abs;
  if (typeof absFn !== 'function') {
    return undefined;
  }
  const a = absFn.call(x) as { toNumber?: () => number; value?: unknown };
  if (typeof a.toNumber === 'function') {
    const v = a.toNumber();
    return Number.isFinite(v) ? v : undefined;
  }
  if (typeof a.value === 'bigint') {
    return Number(a.value);
  }
  const v = Number(String(a));
  return Number.isFinite(v) ? v : undefined;
}

/**
 * Perform the 1x1 pivot update of the Bunch-Kaufman algorithm at position
 * ``(k,k)``, overwriting ``A`` in place.
 *
 * @see Reference: sage/matrix/matrix2.pyx:_block_ldlt_pivot1x1
 */
function _block_ldlt_pivot1x1<R extends FieldElement>(A: Matrix<R>, k: number): void {
  const n = A.nrows;
  const pivot = A.get(k, k);
  if (pivot.isZero()) {
    throw new ZeroDivisionError('zero pivot');
  }
  const pivotInv = getInverse(pivot);

  for (let i = 0; i < n - k - 1; i++) {
    for (let j = 0; j <= i; j++) {
      const val = A.get(k + 1 + i, k + 1 + j).sub(
        A.get(k + 1 + i, k).mul(A.get(k, k + 1 + j)).mul(pivotInv) as R
      ) as R;
      A.set(k + 1 + i, k + 1 + j, val);
      A.set(k + 1 + j, k + 1 + i, _conjugate(val));
    }
  }

  for (let i = 0; i < n - k - 1; i++) {
    A.set(k + i + 1, k, A.get(k + i + 1, k).mul(pivotInv) as R);
  }
}

/**
 * Swap rows ``i`` and ``j`` of ``A`` in place.
 */
function _swap_rows<R extends RingElement>(A: Matrix<R>, i: number, j: number): void {
  if (i === j) return;
  for (let c = 0; c < A.ncols; c++) {
    const t = A.get(i, c);
    A.set(i, c, A.get(j, c));
    A.set(j, c, t);
  }
}

/**
 * Swap columns ``i`` and ``j`` of ``A`` in place.
 */
function _swap_columns<R extends RingElement>(A: Matrix<R>, i: number, j: number): void {
  if (i === j) return;
  for (let r = 0; r < A.nrows; r++) {
    const t = A.get(r, i);
    A.set(r, i, A.get(r, j));
    A.set(r, j, t);
  }
}

/**
 * The user-unfriendly block-LDL^T factorization: returns the permutation as an
 * array, a matrix whose lower-triangular part is L, and the list of diagonal
 * blocks (each 1x1 or 2x2).
 *
 * @see Reference: sage/matrix/matrix2.pyx:_block_ldlt
 */
function _block_ldlt<R extends FieldElement>(
  matrix: Matrix<R>,
  classical: boolean
): [number[], Matrix<R>, Matrix<R>[]] {
  const ring = matrix.base_ring;
  const n = matrix.nrows;
  const A = matrix.copy();

  const p: number[] = [];
  for (let i = 0; i < n; i++) {
    p.push(i);
  }

  const d: Matrix<R>[] = [];

  // The magic constant (1 + sqrt(17))/8 used by Bunch-Kaufman.
  const alpha = 0.6403882032022076;

  let k = 0;
  while (k < n) {
    const A_kk = A.get(k, k);

    if (k === n - 1) {
      d.push(new Matrix<R>(ring, 1, 1, [[A_kk]]));
      k += 1;
      continue;
    }

    if (classical) {
      // Back door giving the standard non-block, non-pivoting LDL^T.
      if (A_kk.isZero()) {
        throw new ValueError('matrix has no classical LDL^T factorization');
      }
      d.push(new Matrix<R>(ring, 1, 1, [[A_kk]]));
      _block_ldlt_pivot1x1(A, k);
      k += 1;
      continue;
    }

    // Step (1) of Higham / Step (1) of Bunch and Kaufman: largest subdiagonal
    // entry (in magnitude) of column k.
    let omega_1 = 0;
    let r = -1;
    let exact = false; // true when the base ring has no absolute value
    for (let i = k + 1; i < n; i++) {
      const a = _bk_abs(A.get(i, k));
      if (a === undefined) {
        exact = true;
        break;
      }
      if (a > omega_1) {
        omega_1 = a;
        r = i;
      }
    }

    if (exact) {
      // No absolute value on the base ring: use an exact "first nonzero"
      // pivoting rule.  Any choice of nonzero pivot yields a valid
      // factorization; only numerical stability (irrelevant here) is lost.
      if (!A_kk.isZero()) {
        d.push(new Matrix<R>(ring, 1, 1, [[A_kk]]));
        _block_ldlt_pivot1x1(A, k);
        k += 1;
        continue;
      }
      let rr = -1;
      for (let i = k + 1; i < n; i++) {
        if (!A.get(i, k).isZero()) {
          rr = i;
          break;
        }
      }
      if (rr === -1) {
        // Column k is zero below the diagonal, and A_kk is zero too.
        d.push(new Matrix<R>(ring, 1, 1, [[A_kk]]));
        k += 1;
        continue;
      }
      if (!A.get(rr, rr).isZero()) {
        d.push(new Matrix<R>(ring, 1, 1, [[A.get(rr, rr)]]));
        _swap_columns(A, k, rr);
        _swap_rows(A, k, rr);
        const t = p[k]!;
        p[k] = p[rr]!;
        p[rr] = t;
        _block_ldlt_pivot1x1(A, k);
        k += 1;
        continue;
      }
      k = _block_ldlt_pivot2x2(A, p, d, k, rr);
      continue;
    }

    if (omega_1 === 0) {
      // A looks like [[a, 0], [0, B]]: record the 1x1 pivot and move on.
      d.push(new Matrix<R>(ring, 1, 1, [[A_kk]]));
      k += 1;
      continue;
    }

    const abs_A_kk = _bk_abs(A_kk)!;
    if (abs_A_kk > alpha * omega_1) {
      // First case of Higham's Step (1) / B&K's Step (2): 1x1 pivot in place.
      d.push(new Matrix<R>(ring, 1, 1, [[A_kk]]));
      _block_ldlt_pivot1x1(A, k);
      k += 1;
      continue;
    }

    // B&K's Step (3): largest off-diagonal entry (in magnitude) of column r.
    let omega_r = 0;
    for (let j = k; j < r; j++) {
      const a = _bk_abs(A.get(r, j))!;
      if (a > omega_r) {
        omega_r = a;
      }
    }

    if (abs_A_kk * omega_r >= alpha * omega_1 * omega_1) {
      // Higham's Step (2) / B&K's Step (4).
      d.push(new Matrix<R>(ring, 1, 1, [[A_kk]]));
      _block_ldlt_pivot1x1(A, k);
      k += 1;
      continue;
    }

    const A_rr = A.get(r, r);
    if (_bk_abs(A_rr)! > alpha * omega_r) {
      // Higham's Step (3) / B&K's Step (5): 1x1 pivot after swapping k and r.
      d.push(new Matrix<R>(ring, 1, 1, [[A_rr]]));
      _swap_columns(A, k, r);
      _swap_rows(A, k, r);
      const t = p[k]!;
      p[k] = p[r]!;
      p[r] = t;
      _block_ldlt_pivot1x1(A, k);
      k += 1;
      continue;
    }

    // Higham's Step (4) / B&K's Step (6): a 2x2 pivot.
    k = _block_ldlt_pivot2x2(A, p, d, k, r);
  }

  for (let i = 0; i < n; i++) {
    A.set(i, i, ring.one());
  }

  return [p, A, d];
}

/**
 * Perform the 2x2 pivot step of the Bunch-Kaufman algorithm, overwriting ``A``,
 * ``p`` and ``d`` in place.  Returns the next value of ``k``.
 *
 * @see Reference: sage/matrix/matrix2.pyx:_block_ldlt
 */
function _block_ldlt_pivot2x2<R extends FieldElement>(
  A: Matrix<R>,
  p: number[],
  d: Matrix<R>[],
  k: number,
  r: number
): number {
  const ring = A.base_ring;
  const n = A.nrows;

  _swap_columns(A, k + 1, r);
  _swap_rows(A, k + 1, r);
  const t = p[k + 1]!;
  p[k + 1] = p[r]!;
  p[r] = t;

  // The top-left 2x2 submatrix starting at (k,k) is the pivot.
  const e00 = A.get(k, k);
  const e01 = A.get(k, k + 1);
  const e10 = A.get(k + 1, k);
  const e11 = A.get(k + 1, k + 1);
  d.push(
    new Matrix<R>(ring, 2, 2, [
      [e00, e01],
      [e10, e11],
    ])
  );

  // X = C * E^{-1} where C = A[k+2:, k:k+2].
  const det = e00.mul(e11).sub(e01.mul(e10) as R) as R;
  if (det.isZero()) {
    throw new ValueError('matrix has no block LDL^T factorization: singular 2x2 pivot');
  }
  const detInv = getInverse(det);
  const i00 = e11.mul(detInv) as R;
  const i01 = e01.neg().mul(detInv) as R;
  const i10 = e10.neg().mul(detInv) as R;
  const i11 = e00.mul(detInv) as R;

  const rows = n - k - 2;
  const X: R[][] = [];
  for (let i = 0; i < rows; i++) {
    const c0 = A.get(k + 2 + i, k);
    const c1 = A.get(k + 2 + i, k + 1);
    X.push([
      c0.mul(i00).add(c1.mul(i10)) as R,
      c0.mul(i01).add(c1.mul(i11)) as R,
    ]);
  }

  // schur_complement = B - X * C^*
  const schur: R[][] = [];
  for (let i = 0; i < rows; i++) {
    schur.push([]);
    for (let j = 0; j < rows; j++) {
      const c0 = _conjugate(A.get(k + 2 + j, k));
      const c1 = _conjugate(A.get(k + 2 + j, k + 1));
      const prod = X[i]![0]!.mul(c0).add(X[i]![1]!.mul(c1)) as R;
      schur[i]!.push(A.get(k + 2 + i, k + 2 + j).sub(prod) as R);
    }
  }

  for (let i = 0; i < rows; i++) {
    for (let j = 0; j <= i; j++) {
      A.set(k + 2 + i, k + 2 + j, schur[i]![j]!);
      A.set(k + 2 + j, k + 2 + i, schur[j]![i]!);
    }
  }

  A.set(k + 1, k, ring.zero());
  for (let i = 0; i < rows; i++) {
    for (let j = 0; j < 2; j++) {
      A.set(k + i + 2, k + j, X[i]![j]!);
    }
  }

  return k + 2;
}

/**
 * Compute a block-LDL^T factorization of a Hermitian matrix.
 *
 * Returns a triple (P, L, D) such that ``A == P*L*D*L^*·P^T``, equivalently
 * ``P.transpose()*A*P == L*D*L.transpose()``, where
 * - P is a permutation matrix,
 * - L is unit lower-triangular,
 * - D is block diagonal with blocks of size one or two.
 *
 * With ``classical=true`` the permutation matrix is the identity and all blocks
 * are 1x1; a :class:`ValueError` is raised when no classical factorization
 * exists.
 *
 * ALGORITHM: "Algorithm A" of Bunch and Kaufman.
 *
 * @param matrix - A symmetric or Hermitian matrix
 * @param classical - Whether to force the classical non-block factorization
 * @returns Triple (P, L, D)
 * @see Reference: sage/matrix/matrix2.pyx:block_ldlt
 * @see Deviation: over base rings without an absolute value (finite fields, say)
 *   the Bunch-Kaufman magnitude comparisons are meaningless, so we fall back to
 *   an exact "first nonzero pivot" rule.  The factorization identity still
 *   holds; only the choice of permutation may differ from Sage's.
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

  const [p, L, d] = _block_ldlt(matrix, classical);

  // P[i,j] == 1 exactly when p[j] == i
  const P = zero_matrix(ring, n);
  for (let j = 0; j < n; j++) {
    P.set(p[j]!, j, ring.one());
  }

  // D is the block-diagonal matrix built from the blocks in d
  const D = zero_matrix(ring, n);
  let offset = 0;
  for (const block of d) {
    for (let i = 0; i < block.nrows; i++) {
      for (let j = 0; j < block.ncols; j++) {
        D.set(offset + i, offset + j, block.get(i, j));
      }
    }
    offset += block.nrows;
  }

  // Overwrite the strict upper-triangular part of L, which still holds scratch.
  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      L.set(i, j, ring.zero());
    }
  }

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

  // Row-reduce to RREF, recording the row operations in U so that U*matrix == S.
  const S = matrix.copy();
  const U = identity_matrix(ring, m);
  const V = identity_matrix(ring, n);

  const pivotCols = _echelonize_in_place(S, U);
  const r = pivotCols.length; // rank

  if (!transformation) {
    // Over a field the Smith form is diagonal with r ones followed by zeros.
    const D = zero_matrix(ring, m, n);
    for (let i = 0; i < r; i++) {
      D.set(i, i, ring.one());
    }
    return D;
  }

  // S is in RREF: row i has a leading one at column pivotCols[i], and possibly
  // non-zero entries in the non-pivot columns.  Clear those with column
  // operations, recording each one in V so that S == U*matrix*V throughout.
  const pivotSet = new Set(pivotCols);
  for (let c = 0; c < n; c++) {
    if (pivotSet.has(c)) {
      continue;
    }
    for (let i = 0; i < r; i++) {
      const f = S.get(i, c);
      if (f.isZero()) {
        continue;
      }
      const pc = pivotCols[i]!;
      // column_c -= f * column_{pc}
      for (let row = 0; row < m; row++) {
        S.set(row, c, S.get(row, c).sub(f.mul(S.get(row, pc)) as R) as R);
      }
      for (let row = 0; row < n; row++) {
        V.set(row, c, V.get(row, c).sub(f.mul(V.get(row, pc)) as R) as R);
      }
    }
  }

  // Move the pivot columns onto the diagonal with column swaps.
  for (let i = 0; i < r; i++) {
    const c = pivotCols[i]!;
    if (c === i) {
      continue;
    }
    for (let row = 0; row < m; row++) {
      const t = S.get(row, i);
      S.set(row, i, S.get(row, c));
      S.set(row, c, t);
    }
    for (let row = 0; row < n; row++) {
      const t = V.get(row, i);
      V.set(row, i, V.get(row, c));
      V.set(row, c, t);
    }
    // Record the swap so that later pivots are still found in the right place.
    for (let k = i + 1; k < r; k++) {
      if (pivotCols[k] === i) {
        pivotCols[k] = c;
        break;
      }
    }
    pivotCols[i] = i;
  }

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

  // Compute the eigenvalues of the matrix, with multiplicities: ``evPairs`` is
  // a list of pairs, each first entry a root of the characteristic polynomial
  // and each second entry the corresponding multiplicity.
  let evPairs: Array<[R, number]>;

  if (eigenvalues !== undefined) {
    evPairs = eigenvalues;
    if (check_input) {
      // The provided eigenvalues must reproduce the characteristic polynomial.
      const cp = _charpoly(matrix);
      let prod: R[] = [ring.one()];
      for (const [z, i] of evPairs) {
        for (let t = 0; t < i; t++) {
          prod = _multiply_by_linear(prod, z);
        }
      }
      let equal = prod.length === cp.length;
      if (equal) {
        for (let i = 0; i < cp.length; i++) {
          if (!cp[i]!.eq(prod[i]!)) {
            equal = false;
            break;
          }
        }
      }
      if (!equal) {
        throw new ValueError('The provided list of eigenvalues is not correct.');
      }
    }
  } else {
    // Compute eigenvalues as the roots of the characteristic polynomial, with
    // multiplicities obtained by repeated division by (x - lambda).
    const elementsMethod = (ring as unknown as { elements?: () => Iterable<R> }).elements;
    if (typeof elementsMethod !== 'function') {
      throw new NotImplementedError(
        'jordan_form requires either pre-computed eigenvalues or a finite ring with enumerable elements'
      );
    }

    let cp = _charpoly(matrix);
    evPairs = [];

    for (const elem of elementsMethod.call(ring)) {
      let multiplicity = 0;
      while (cp.length > 1 && _evaluate_poly(cp, elem).isZero()) {
        cp = _divide_by_linear(cp, elem);
        multiplicity++;
      }
      if (multiplicity > 0) {
        evPairs.push([elem, multiplicity]);
      }
    }
  }

  // Check that the sum of the multiplicities equals n
  const totalMult = evPairs.reduce((sum, [, mult]) => sum + mult, 0);
  if (totalMult < n) {
    throw new ArithmeticError(`Some eigenvalue does not exist in ${String(ring)}.`);
  }

  // Compute the block information.  ``blocks`` is a list of pairs, each first
  // entry a root and each second entry the size of a block.  Note that in
  // general there is more than one block per eigenvalue.
  interface JordanBlock {
    eigenvalue: R;
    size: number;
  }

  const blocks: JordanBlock[] = [];

  for (const [lambda, mult] of evPairs) {
    if (mult === 1) {
      blocks.push({ eigenvalue: lambda, size: 1 });
      continue;
    }

    const I = identity_matrix(ring, n);
    const B = matrix.sub(I.scalar_mul(lambda));
    let C = B;
    const ranks: number[] = [n, _rank(C)];
    let i = 0;
    while (ranks[i]! > ranks[i + 1]! && ranks[i + 1]! > n - mult) {
      C = B.mul(C);
      ranks.push(_rank(C));
      i += 1;
    }

    // The diagram is a partition; its conjugate lists the Jordan block sizes,
    // in decreasing order.
    const diagram: number[] = [];
    for (let j = 0; j + 1 < ranks.length; j++) {
      const v = ranks[j]! - ranks[j + 1]!;
      if (v > 0) {
        diagram.push(v);
      }
    }

    for (const size of _conjugate_partition(diagram)) {
      blocks.push({ eigenvalue: lambda, size });
    }
  }

  // Build the Jordan form matrix
  const J = zero_matrix(ring, n);
  let offset = 0;

  for (const block of blocks) {
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

  // Computing the transformation matrix P requires generalized eigenvector
  // chains, which are not implemented yet.
  throw new NotImplementedError(
    'jordan_form with transformation=true is not yet fully implemented'
  );
}

/**
 * Evaluate a polynomial given by ascending coefficients at ``x``.
 */
function _evaluate_poly<R extends FieldElement>(coeffs: R[], x: R): R {
  if (coeffs.length === 0) {
    throw new ValueError('empty polynomial');
  }
  let val = coeffs[coeffs.length - 1]!;
  for (let i = coeffs.length - 2; i >= 0; i--) {
    val = val.mul(x).add(coeffs[i]!) as R;
  }
  return val;
}

/**
 * Divide a polynomial given by ascending coefficients by ``x - root`` (exact
 * division by synthetic division; the remainder is discarded).
 */
function _divide_by_linear<R extends FieldElement>(coeffs: R[], root: R): R[] {
  const d = coeffs.length - 1;
  const q: R[] = new Array(d);
  let carry = coeffs[d]!;
  for (let i = d - 1; i >= 0; i--) {
    q[i] = carry;
    carry = coeffs[i]!.add(carry.mul(root)) as R;
  }
  return q;
}

/**
 * Multiply a polynomial given by ascending coefficients by ``x - root``.
 */
function _multiply_by_linear<R extends FieldElement>(coeffs: R[], root: R): R[] {
  const zero = coeffs[0]!.sub(coeffs[0]!) as R;
  const out: R[] = new Array(coeffs.length + 1);
  for (let i = 0; i < out.length; i++) {
    out[i] = zero;
  }
  for (let i = 0; i < coeffs.length; i++) {
    out[i + 1] = out[i + 1]!.add(coeffs[i]!) as R;
    out[i] = out[i]!.sub(coeffs[i]!.mul(root) as R) as R;
  }
  return out;
}

/**
 * Return the conjugate of a partition given in weakly decreasing order.
 */
function _conjugate_partition(part: number[]): number[] {
  if (part.length === 0) {
    return [];
  }
  const maxPart = part[0]!;
  const conj: number[] = [];
  for (let j = 1; j <= maxPart; j++) {
    let count = 0;
    for (const p of part) {
      if (p >= j) {
        count++;
      }
    }
    conj.push(count);
  }
  return conj;
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
  return pivots(matrix).length;
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
    const pivotCols = pivots(M);
    const pivotSet = new Set(pivotCols);

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
          const pivotIdx = pivotCols.indexOf(j);
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

  // Convert the Gram matrix to bigint for exact integer arithmetic.
  const G0: bigint[][] = [];
  for (let i = 0; i < n; i++) {
    G0.push([]);
    for (let j = 0; j < n; j++) {
      G0[i]!.push(_to_bigint(matrix.get(i, j)));
    }
  }

  // U holds the current basis in its COLUMNS: b_j = sum_i U[i][j] * e_i, so the
  // current Gram matrix is G = U^T * G0 * U.
  const U: bigint[][] = [];
  for (let i = 0; i < n; i++) {
    U.push([]);
    for (let j = 0; j < n; j++) {
      U[i]!.push(i === j ? 1n : 0n);
    }
  }

  const G: bigint[][] = [];
  for (let i = 0; i < n; i++) {
    G.push(new Array<bigint>(n).fill(0n));
  }

  function recomputeG(): void {
    // temp = G0 * U
    const temp: bigint[][] = [];
    for (let i = 0; i < n; i++) {
      temp.push([]);
      for (let j = 0; j < n; j++) {
        let sum = 0n;
        for (let k = 0; k < n; k++) {
          sum += G0[i]![k]! * U[k]![j]!;
        }
        temp[i]!.push(sum);
      }
    }
    // G = U^T * temp
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

  // Exact rational Gram-Schmidt data derived from G:
  //   mu[i][j] = <b_i, b*_j> / <b*_j, b*_j>,  B[i] = <b*_i, b*_i>
  const mu: Rat[][] = [];
  const B: Rat[] = [];

  function computeGramSchmidt(): void {
    mu.length = 0;
    B.length = 0;
    for (let i = 0; i < n; i++) {
      mu.push([]);
      let Bi = ratFromBigInt(G[i]![i]!);
      for (let j = 0; j < i; j++) {
        let mij = ratFromBigInt(G[i]![j]!);
        for (let k = 0; k < j; k++) {
          mij = ratSub(mij, ratMul(ratMul(mu[i]![k]!, mu[j]![k]!), B[k]!));
        }
        if (ratIsZero(B[j]!)) {
          throw new ValueError(
            'qflllgram did not return a square matrix, perhaps the matrix is not positive definite'
          );
        }
        mij = ratDiv(mij, B[j]!);
        mu[i]!.push(mij);
        Bi = ratSub(Bi, ratMul(ratMul(mij, mij), B[j]!));
      }
      B.push(Bi);
    }
  }

  // RED(k, l): size-reduce column k against column l.
  function reduce(k: number, l: number): void {
    const q = ratRoundToNearest(mu[k]![l]!);
    if (q === 0n) {
      return;
    }
    for (let row = 0; row < n; row++) {
      U[row]![k] = U[row]![k]! - q * U[row]![l]!;
    }
    recomputeG();
    computeGramSchmidt();
  }

  function swapColumns(i: number, j: number): void {
    for (let row = 0; row < n; row++) {
      const t = U[row]![i]!;
      U[row]![i] = U[row]![j]!;
      U[row]![j] = t;
    }
    recomputeG();
    computeGramSchmidt();
  }

  recomputeG();
  computeGramSchmidt();

  if (ratIsZero(B[0]!)) {
    // PARI's qflllgram bails out on such Gram matrices.
    throw new ValueError(
      'qflllgram did not return a square matrix, perhaps the matrix is not positive definite'
    );
  }

  // LLL with delta = 3/4.
  let k = 1;
  let guard = 0;
  const maxSteps = 1000000;
  while (k < n) {
    if (++guard > maxSteps) {
      throw new ArithmeticError('infinite loop while LLL-reducing the Gram matrix');
    }

    reduce(k, k - 1);

    // Lovasz condition: B[k] >= (3/4 - mu[k][k-1]^2) * B[k-1]
    const mk = mu[k]![k - 1]!;
    const bound = ratMul(ratSub(ratFromFraction(3n, 4n), ratMul(mk, mk)), B[k - 1]!);
    if (ratCompare(B[k]!, bound) < 0) {
      swapColumns(k - 1, k);
      k = k > 1 ? k - 1 : 1;
    } else {
      for (let l = k - 2; l >= 0; l--) {
        reduce(k, l);
      }
      k += 1;
    }
  }

  // PARI's qflllgram returns a transformation of determinant 1.  Negating a
  // column leaves U^T*G0*U unchanged and flips the sign of the determinant.
  if (computeDetSign(U) === -1n) {
    for (let i = 0; i < n; i++) {
      U[i]![n - 1] = -U[i]![n - 1]!;
    }
  }

  // Convert U back to a matrix over the base ring.
  const result = new Matrix<R>(ring, n, n);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      result.set(i, j, _from_bigint(ring, U[i]![j]!));
    }
  }

  return result;
}

/** Exact rational number as a normalized [numerator, denominator] pair. */
type Rat = [bigint, bigint];

function ratGcd(a: bigint, b: bigint): bigint {
  a = a < 0n ? -a : a;
  b = b < 0n ? -b : b;
  while (b !== 0n) {
    const t = b;
    b = a % b;
    a = t;
  }
  return a;
}

function ratNormalize(num: bigint, den: bigint): Rat {
  if (den === 0n) {
    throw new ZeroDivisionError('rational division by zero');
  }
  if (den < 0n) {
    num = -num;
    den = -den;
  }
  if (num === 0n) {
    return [0n, 1n];
  }
  const g = ratGcd(num, den);
  return [num / g, den / g];
}

function ratFromBigInt(x: bigint): Rat {
  return [x, 1n];
}

function ratFromFraction(num: bigint, den: bigint): Rat {
  return ratNormalize(num, den);
}

function ratIsZero(x: Rat): boolean {
  return x[0] === 0n;
}

function ratSub(a: Rat, b: Rat): Rat {
  return ratNormalize(a[0] * b[1] - b[0] * a[1], a[1] * b[1]);
}

function ratMul(a: Rat, b: Rat): Rat {
  return ratNormalize(a[0] * b[0], a[1] * b[1]);
}

function ratDiv(a: Rat, b: Rat): Rat {
  if (b[0] === 0n) {
    throw new ZeroDivisionError('rational division by zero');
  }
  return ratNormalize(a[0] * b[1], a[1] * b[0]);
}

function ratCompare(a: Rat, b: Rat): number {
  const lhs = a[0] * b[1];
  const rhs = b[0] * a[1];
  if (lhs < rhs) return -1;
  if (lhs > rhs) return 1;
  return 0;
}

/** Round a rational to the nearest integer, halves away from zero. */
function ratRoundToNearest(x: Rat): bigint {
  const [num, den] = x;
  const q = num / den;
  const r = num % den;
  const twice = 2n * (r < 0n ? -r : r);
  if (twice >= den) {
    return r > 0n ? q + 1n : q - 1n;
  }
  return q;
}

/** Extract a bigint from a ring element in as many ways as we can. */
function _to_bigint(entry: unknown): bigint {
  if (typeof entry === 'bigint') {
    return entry;
  }
  const withValue = entry as { value?: unknown; toBigInt?: () => bigint };
  if (typeof withValue.value === 'bigint') {
    return withValue.value;
  }
  if (typeof withValue.toBigInt === 'function') {
    return withValue.toBigInt();
  }
  return BigInt(Number(entry));
}

/** Build a ring element from a bigint. */
function _from_bigint<R extends RingElement>(ring: CoefficientRing<R>, val: bigint): R {
  const withFrom = ring as unknown as {
    fromBigInt?: (x: bigint) => R;
    __call__?: (x: bigint) => R;
  };
  if (typeof withFrom.fromBigInt === 'function') {
    return withFrom.fromBigInt(val);
  }
  if (typeof withFrom.__call__ === 'function') {
    return withFrom.__call__(val);
  }
  let elem = ring.zero();
  const one = ring.one();
  if (val > 0n) {
    for (let k = 0n; k < val; k++) {
      elem = elem.add(one) as R;
    }
  } else if (val < 0n) {
    for (let k = 0n; k > val; k--) {
      elem = elem.sub(one) as R;
    }
  }
  return elem;
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
 * A positive definite matrix A has a unique positive definite matrix M such
 * that M^2 = A.
 *
 * ALGORITHM: exactly Sage's — diagonalize, take the element-wise square roots of
 * the eigenvalues, and conjugate back:
 * ``d, L = self.eigenmatrix_left(); return L.inverse()*diagonal_matrix([sqrt(a) for a in d.diagonal()])*L``.
 *
 * @param matrix - A positive definite matrix
 * @param check_positivity - Whether to verify positive definiteness (default: true)
 * @returns The principal square root M such that M^2 = A, or false if not positive definite
 * @see Reference: sage/matrix/matrix2.pyx:principal_square_root
 * @see Deviation: Sage first calls ``is_positive_definite()`` and returns
 *   ``False`` when the matrix is not positive definite.  Positive definiteness
 *   is meaningless over the base rings this generic implementation supports
 *   (finite fields have no ordering), so the check is skipped.  Likewise, over a
 *   finite field "the" square root of an eigenvalue is only defined up to sign,
 *   so the result is *a* square root, not the principal one.
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

  // Diagonalize: A*P == P*D, i.e. A == P*D*P^{-1}.
  const [D, P] = eigenmatrix_right(matrix);

  if (_rank(P) < n) {
    throw new ArithmeticError(
      'principal_square_root requires the matrix to be diagonalizable over its base ring'
    );
  }

  // sqrt(D), entry-wise
  const sqrtD = zero_matrix(ring, n);
  for (let i = 0; i < n; i++) {
    const entry = D.get(i, i) as unknown as { sqrt?: () => R };
    if (typeof entry.sqrt !== 'function') {
      throw new NotImplementedError(
        'principal_square_root requires the base ring to support element-wise square roots'
      );
    }
    sqrtD.set(i, i, entry.sqrt());
  }

  // A^{1/2} = P * sqrt(D) * P^{-1}
  return P.mul(sqrtD).mul(matrixInverse(P));
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

    // Sage appends ``B.kernel()`` for the primal decomposition, and Sage's
    // ``kernel`` is the *left* kernel (matrix2.pyx:5503), i.e. the right kernel
    // of B^T.  The dual uses ``B.transpose().kernel()``, the left kernel of B^T,
    // which is the right kernel of B.
    const kernelBasis = computeKernelBasis(gA.transpose());

    const isIrreducible = m === 1;
    result.push([kernelBasis, isIrreducible]);

    if (dual) {
      const dualKernelBasis = computeKernelBasis(gA);
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
  const pivotCols = pivots(A);
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
    const rankBefore = pivots(basisMatrix).length;
    const rankAfter = pivots(augMatrix).length;

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
    const rankBefore = pivots(basisMatrix).length;
    const rankAfter = pivots(augMatrix).length;

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
 * Normalize the ``shifts``/``degrees`` arguments of the Krylov routines.
 */
function _krylov_normalize_args<R extends FieldElement>(
  matrix: Matrix<R>,
  M: Matrix<R>,
  shifts?: number[],
  degrees?: number | number[]
): [number[], number[]] {
  const m = matrix.nrows;
  const n = matrix.ncols;

  if (M.nrows !== n || M.ncols !== n) {
    throw new ValueError('M does not have correct dimensions');
  }

  let sh: number[];
  if (shifts === undefined) {
    sh = new Array<number>(m).fill(0);
  } else {
    if (shifts.length !== m) {
      throw new ValueError('x must be a list of the right length');
    }
    sh = shifts.slice();
  }

  let deg: number[];
  if (degrees === undefined) {
    deg = new Array<number>(m).fill(n);
  } else if (typeof degrees === 'number') {
    deg = new Array<number>(m).fill(degrees);
  } else {
    if (degrees.length !== m) {
      throw new ValueError('x must be a list of the right length');
    }
    deg = degrees.slice();
  }

  if (m > 0 && Math.min(...deg) < 0) {
    throw new ValueError('degrees must not contain a negative bound');
  }

  return [sh, deg];
}

/**
 * Return the (row, exponent) coordinates of the Krylov matrix rows, sorted by
 * the priority given by ``shifts``.
 *
 * Each entry is a triple ``[c, d, i]``: the row `E_c M^d` of the Krylov matrix,
 * whose position *before* sorting is ``i``.  Sorting is by the lexicographic
 * order on ``(shifts[c] + d, c)``.
 *
 * @see Reference: sage/matrix/matrix2.pyx:_krylov_row_coordinates
 */
function _krylov_row_coordinates(
  m: number,
  shifts: number[],
  degrees: number[],
  row_pairs?: Array<[number, number]>
): Array<[number, number, number]> {
  const blocks = degrees.length === 0 ? 1 : Math.max(...degrees) + 1;

  let pairs: Array<[number, number]>;
  if (row_pairs === undefined) {
    pairs = [];
    for (let j = 0; j < blocks; j++) {
      for (let i = 0; i < m; i++) {
        if (j <= degrees[i]!) {
          pairs.push([i, j]);
        }
      }
    }
  } else {
    pairs = row_pairs.filter((row) => row[1] <= degrees[row[0]]!);
  }

  const rows: Array<[number, number, number]> = pairs.map((row, i) => [row[0], row[1], i]);
  rows.sort((a, b) => {
    const ka = shifts[a[0]]! + a[1];
    const kb = shifts[b[0]]! + b[1];
    if (ka !== kb) return ka - kb;
    return a[0] - b[0];
  });
  return rows;
}

/**
 * Return the Krylov matrix built from the rows of ``matrix`` and using the
 * multiplication matrix ``M``.
 *
 * Writing E for ``matrix`` (an m x n matrix with rows E_0, ..., E_{m-1}) and M
 * for the n x n acting matrix, the Krylov matrix stacks the iterates
 * ``E_i * M^j`` for all ``0 <= i < m`` and ``0 <= j <= degrees[i]``
 * (**inclusive**), so it has ``d_0 + ... + d_{m-1} + m`` rows.  The rows are
 * ordered by the priority defined by ``shifts``: ascending in
 * ``(shifts[i] + j, i)``.
 *
 * @param matrix - The matrix E (m x n)
 * @param M - The acting matrix (n x n)
 * @param shifts - Row priority shifts (default: all zero)
 * @param degrees - Degree bounds, an integer or a list of m integers
 *   (default: n for every row)
 * @returns The Krylov matrix
 * @see Reference: sage/matrix/matrix2.pyx:krylov_matrix
 */
export function krylov_matrix<R extends FieldElement>(
  matrix: Matrix<R>,
  M: Matrix<R>,
  shifts?: number[],
  degrees?: number | number[]
): Matrix<R> {
  if (!M.is_square()) {
    throw new ArithmeticError('acting matrix M must be square');
  }

  const m = matrix.nrows;
  const n = matrix.ncols;
  const ring = matrix.base_ring;

  const [sh, deg] = _krylov_normalize_args(matrix, M, shifts, degrees);

  const maxDegree = m === 0 ? 0 : Math.max(...deg);

  // iterates[j][i] = row i of E * M^j
  const iterates: R[][][] = [];
  let current: R[][] = [];
  for (let i = 0; i < m; i++) {
    current.push(matrix.row(i));
  }
  iterates.push(current);

  for (let d = 1; d <= maxDegree; d++) {
    const next: R[][] = [];
    for (let i = 0; i < m; i++) {
      const row: R[] = [];
      for (let j = 0; j < n; j++) {
        let sum = ring.zero();
        for (let k = 0; k < n; k++) {
          sum = sum.add(current[i]![k]!.mul(M.get(k, j))) as R;
        }
        row.push(sum);
      }
      next.push(row);
    }
    iterates.push(next);
    current = next;
  }

  const coords = _krylov_row_coordinates(m, sh, deg);
  const rows: R[][] = coords.map(([c, d]) => iterates[d]![c]!);

  return new Matrix(ring, rows.length, n, rows);
}

/**
 * Return the matrix formed by stacking the first linearly independent rows of
 * the Krylov matrix, in the order defined by ``shifts``.
 *
 * That is, ``K.matrix_from_rows(K.pivot_rows())`` where ``K`` is
 * ``krylov_matrix(matrix, M, shifts, degrees)``.  Together with the basis, Sage
 * returns for each row its position in ``K`` and the pair ``(i, j)`` with the
 * row equal to ``E_i * M^j``; that information is returned when
 * ``output_rows`` is true.
 *
 * @param matrix - The matrix E
 * @param M - The acting matrix
 * @param shifts - Row priority shifts
 * @param degrees - Degree bounds
 * @param output_rows - Whether to also return the row coordinates (default: true)
 * @param algorithm - Algorithm to use ('naive' or 'elimination'); both produce
 *   the same result here
 * @returns The Krylov basis, together with the row coordinates when requested
 * @see Reference: sage/matrix/matrix2.pyx:krylov_basis
 */
export function krylov_basis<R extends FieldElement>(
  matrix: Matrix<R>,
  M: Matrix<R>,
  shifts?: number[],
  degrees?: number | number[],
  output_rows: boolean = true,
  algorithm?: string
): Matrix<R> | [Matrix<R>, Array<[number, number, number]>] {
  if (algorithm !== undefined && algorithm !== 'naive' && algorithm !== 'elimination') {
    throw new ValueError('algorithm must be one of None, "naive" or "elimination"');
  }

  const m = matrix.nrows;
  const ring = matrix.base_ring;
  const [sh, deg] = _krylov_normalize_args(matrix, M, shifts, degrees);

  const K = krylov_matrix(matrix, M, sh, deg);

  // The first linearly independent rows of K (its row rank profile).
  const rowProfile = pivot_rows(K);

  const rows: R[][] = rowProfile.map((i) => K.row(i));
  const kmat = new Matrix(ring, rows.length, K.ncols, rows);

  if (!output_rows) {
    return kmat;
  }

  const coords = _krylov_row_coordinates(m, sh, deg);
  const profile: Array<[number, number, number]> = rowProfile.map((i) => [
    coords[i]![0],
    coords[i]![1],
    i,
  ]);

  return [kmat, profile];
}

/**
 * Return a basis in canonical form for the left kernel of the Krylov matrix of
 * ``(matrix, M)`` with rows ordered according to ``shifts``.
 *
 * Following Sage, let ``B`` be the Krylov basis computed by
 * :func:`krylov_basis` with the same parameters, and let
 * ``[delta_0, ..., delta_{m-1}]`` be the exponents of first linear dependency
 * for each row (``delta_i = 0`` when row ``i`` never appears in ``B``, else one
 * more than the largest exponent appearing).  The result is a basis of the left
 * kernel of the Krylov matrix built from ``matrix`` and ``M`` with degree bounds
 * ``delta``.  It has ``m`` rows and ``m + rank(B)`` columns.
 *
 * @param matrix - The matrix E
 * @param M - The acting matrix
 * @param shifts - Row priority shifts
 * @param degrees - Degree bounds
 * @param output_rows - Whether to also return the row coordinates (default: true)
 * @param variable - Variable name for the polynomial-matrix representation
 *   (not supported)
 * @param basis_algorithm - Algorithm for computing the Krylov basis
 * @returns The Krylov kernel basis, with the row coordinates when requested
 * @see Reference: sage/matrix/matrix2.pyx:krylov_kernel_basis
 */
export function krylov_kernel_basis<R extends FieldElement>(
  matrix: Matrix<R>,
  M: Matrix<R>,
  shifts?: number[],
  degrees?: number | number[],
  output_rows: boolean = true,
  variable?: string,
  basis_algorithm?: string
): Matrix<R> | [Matrix<R>, Array<[number, number, number]>] {
  if (variable !== undefined) {
    throw new NotImplementedError(
      'the polynomial matrix representation of krylov_kernel_basis is not implemented'
    );
  }

  const m = matrix.nrows;
  const [sh, deg] = _krylov_normalize_args(matrix, M, shifts, degrees);

  const [, profile] = krylov_basis(matrix, M, sh, deg, true, basis_algorithm) as [
    Matrix<R>,
    Array<[number, number, number]>,
  ];

  // delta_i = 1 + (largest exponent of row i selected for the basis), or 0
  const delta = new Array<number>(m).fill(0);
  for (const [i, j] of profile) {
    delta[i] = Math.max(delta[i]!, j + 1);
  }

  const A = krylov_matrix(matrix, M, sh, delta);

  // The left kernel of A is the right kernel of A^T.
  const kernel = computeKernelBasis(A.transpose());

  if (!output_rows) {
    return kernel;
  }

  return [kernel, _krylov_row_coordinates(m, sh, delta)];
}
