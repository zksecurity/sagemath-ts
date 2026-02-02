/**
 * @module sage/matrix/matrix_decompositions_additions
 * @description Additional matrix decompositions for real matrices using IEEE 754 doubles.
 *
 * This module provides SVD, QR (Householder), and LU decompositions specifically
 * for real matrices represented as 2D number arrays.
 *
 * @see Deviation: Uses IEEE 754 double precision instead of arbitrary precision.
 *   See DEVIATIONS.md for details on SVD/Real Matrix Decompositions.
 *
 * Port of: sage/matrix/matrix_double_dense.pyx
 */

import { ArithmeticError } from '../errors.js';

// ============================================================================
// SVD (Singular Value Decomposition) for Real Matrices
// ============================================================================

/**
 * Result type for SVD decomposition.
 */
export interface SVDResult {
  /** Left singular vectors (m x m orthogonal matrix) */
  U: number[][];
  /** Singular values (as vector, in descending order) */
  S: number[];
  /** Right singular vectors (n x n orthogonal matrix) */
  V: number[][];
}

/**
 * Return the Singular Value Decomposition of a real matrix.
 *
 * Computes matrices U, S, V such that A = U * diag(S) * V^T where:
 * - U is an m x m orthogonal matrix (columns are left singular vectors)
 * - S is a vector of singular values in descending order
 * - V is an n x n orthogonal matrix (columns are right singular vectors)
 *
 * This implementation uses the Jacobi one-sided SVD algorithm, which
 * is numerically stable and works well for dense matrices.
 *
 * @param A - A real matrix as a 2D array of numbers (IEEE 754 doubles)
 * @param options - Optional settings:
 *   - tol: tolerance for convergence (default: 1e-14)
 *   - maxIterations: maximum number of sweeps (default: 100)
 * @returns SVD decomposition {U, S, V}
 * @see Deviation: Uses IEEE 754 double precision instead of arbitrary precision.
 *   See DEVIATIONS.md for details.
 * @see Reference: sage/matrix/matrix_double_dense.pyx:SVD
 *
 * @example
 * ```typescript
 * const A = [[1, 2], [3, 4], [5, 6]];
 * const { U, S, V } = SVD_double(A);
 * // U is 3x3 orthogonal
 * // S contains singular values [9.508..., 0.773...]
 * // V is 2x2 orthogonal
 * // A ≈ U * diag(S) * V^T
 * ```
 */
export function SVD_double(
  A: number[][],
  options?: { tol?: number; maxIterations?: number }
): SVDResult {
  const tol = options?.tol ?? 1e-14;
  const maxIterations = options?.maxIterations ?? 100;

  const m = A.length;
  if (m === 0) {
    return { U: [], S: [], V: [] };
  }

  const n = A[0]!.length;
  if (n === 0) {
    return { U: Array.from({ length: m }, () => []), S: [], V: [] };
  }

  // Handle case where m < n by computing SVD of A^T and swapping U and V
  if (m < n) {
    const AT = transpose_double(A);
    const { U: V, S, V: U } = SVD_double(AT, options);
    return { U, S, V };
  }

  // Copy A to working matrix B (m x n)
  const B: number[][] = A.map((row) => [...row]);

  // Initialize V as identity (n x n)
  const V: number[][] = [];
  for (let i = 0; i < n; i++) {
    V.push([]);
    for (let j = 0; j < n; j++) {
      V[i]!.push(i === j ? 1 : 0);
    }
  }

  // Jacobi one-sided SVD algorithm
  // Apply Givens rotations to make B^T * B diagonal
  let converged = false;

  for (let sweep = 0; sweep < maxIterations && !converged; sweep++) {
    converged = true;

    for (let i = 0; i < n - 1; i++) {
      for (let j = i + 1; j < n; j++) {
        // Compute elements of B^T * B for columns i and j
        let a_ii = 0;
        let a_jj = 0;
        let a_ij = 0;

        for (let k = 0; k < m; k++) {
          a_ii += B[k]![i]! * B[k]![i]!;
          a_jj += B[k]![j]! * B[k]![j]!;
          a_ij += B[k]![i]! * B[k]![j]!;
        }

        // Check if off-diagonal element is significant
        if (Math.abs(a_ij) > tol * Math.sqrt(a_ii * a_jj)) {
          converged = false;

          // Compute Jacobi rotation angle
          const tau = (a_jj - a_ii) / (2 * a_ij);
          let t: number;
          if (tau >= 0) {
            t = 1 / (tau + Math.sqrt(1 + tau * tau));
          } else {
            t = -1 / (-tau + Math.sqrt(1 + tau * tau));
          }

          const c = 1 / Math.sqrt(1 + t * t);
          const s = t * c;

          // Apply rotation to B (columns i and j)
          for (let k = 0; k < m; k++) {
            const bi = B[k]![i]!;
            const bj = B[k]![j]!;
            B[k]![i] = c * bi - s * bj;
            B[k]![j] = s * bi + c * bj;
          }

          // Apply rotation to V (columns i and j)
          for (let k = 0; k < n; k++) {
            const vi = V[k]![i]!;
            const vj = V[k]![j]!;
            V[k]![i] = c * vi - s * vj;
            V[k]![j] = s * vi + c * vj;
          }
        }
      }
    }
  }

  // Singular values are the norms of columns of B
  const S: number[] = [];

  for (let j = 0; j < n; j++) {
    let norm = 0;
    for (let k = 0; k < m; k++) {
      norm += B[k]![j]! * B[k]![j]!;
    }
    S.push(Math.sqrt(norm));
  }

  // Build U by normalizing columns of B
  const U: number[][] = [];
  for (let i = 0; i < m; i++) {
    U.push([]);
    for (let j = 0; j < n; j++) {
      if (S[j]! > tol) {
        U[i]!.push(B[i]![j]! / S[j]!);
      } else {
        U[i]!.push(0);
      }
    }
  }

  // Extend U to m x m by adding orthogonal vectors using Gram-Schmidt
  for (let k = n; k < m; k++) {
    const v: number[] = [];
    for (let i = 0; i < m; i++) {
      v.push(i === k ? 1 : 0);
    }

    for (let j = 0; j < k; j++) {
      let dot = 0;
      for (let i = 0; i < m; i++) {
        dot += v[i]! * U[i]![j]!;
      }
      for (let i = 0; i < m; i++) {
        v[i]! -= dot * U[i]![j]!;
      }
    }

    let norm = 0;
    for (let i = 0; i < m; i++) {
      norm += v[i]! * v[i]!;
    }
    norm = Math.sqrt(norm);

    if (norm > tol) {
      for (let i = 0; i < m; i++) {
        v[i]! /= norm;
      }
    }

    for (let i = 0; i < m; i++) {
      U[i]!.push(v[i]!);
    }
  }

  // Sort singular values in descending order
  const indices = Array.from({ length: n }, (_, i) => i);
  indices.sort((a, b) => S[b]! - S[a]!);

  const sortedS: number[] = [];
  const sortedU: number[][] = U.map(() => []);
  const sortedV: number[][] = V.map(() => []);

  for (let j = 0; j < n; j++) {
    const idx = indices[j]!;
    sortedS.push(S[idx]!);

    for (let i = 0; i < m; i++) {
      sortedU[i]!.push(U[i]![idx]!);
    }
    for (let i = 0; i < n; i++) {
      sortedV[i]!.push(V[i]![idx]!);
    }
  }

  for (let j = n; j < m; j++) {
    for (let i = 0; i < m; i++) {
      sortedU[i]!.push(U[i]![j]!);
    }
  }

  return { U: sortedU, S: sortedS, V: sortedV };
}

/**
 * Transpose a 2D number array.
 */
function transpose_double(A: number[][]): number[][] {
  if (A.length === 0) return [];
  const m = A.length;
  const n = A[0]!.length;
  const AT: number[][] = [];
  for (let j = 0; j < n; j++) {
    AT.push([]);
    for (let i = 0; i < m; i++) {
      AT[j]!.push(A[i]![j]!);
    }
  }
  return AT;
}

/**
 * Reconstruct matrix A from its SVD: A = U * diag(S) * V^T
 *
 * Useful for verifying the SVD result.
 *
 * @param svd - The SVD result {U, S, V}
 * @returns The reconstructed matrix
 */
export function SVD_reconstruct(svd: SVDResult): number[][] {
  const { U, S, V } = svd;
  const m = U.length;
  if (m === 0) return [];
  const n = V.length;

  const A: number[][] = [];
  for (let i = 0; i < m; i++) {
    A.push([]);
    for (let j = 0; j < n; j++) {
      let sum = 0;
      for (let k = 0; k < S.length; k++) {
        sum += U[i]![k]! * S[k]! * V[j]![k]!;
      }
      A[i]!.push(sum);
    }
  }
  return A;
}

/**
 * Compute the pseudoinverse of a matrix using its SVD.
 *
 * The pseudoinverse A^+ is computed as V * diag(1/S) * U^T,
 * where small singular values (below tolerance) are treated as zero.
 *
 * @param A - A real matrix as a 2D array of numbers
 * @param tol - Tolerance for small singular values (default: 1e-14)
 * @returns The pseudoinverse matrix
 */
export function pseudoinverse_SVD(A: number[][], tol: number = 1e-14): number[][] {
  const { U, S, V } = SVD_double(A, { tol });
  const m = U.length;
  const n = V.length;
  const numSingular = S.length;

  // A^+ = V * diag(1/S) * U^T
  // Result is n x m
  const Aplus: number[][] = [];
  for (let i = 0; i < n; i++) {
    Aplus.push([]);
    for (let j = 0; j < m; j++) {
      let sum = 0;
      for (let k = 0; k < numSingular; k++) {
        if (S[k]! > tol) {
          sum += V[i]![k]! * (1 / S[k]!) * U[j]![k]!;
        }
      }
      Aplus[i]!.push(sum);
    }
  }
  return Aplus;
}

/**
 * Compute the matrix rank using SVD.
 *
 * The rank is the number of singular values above the tolerance.
 *
 * @param A - A real matrix as a 2D array of numbers
 * @param tol - Tolerance for small singular values (default: 1e-14)
 * @returns The numerical rank
 */
export function rank_SVD(A: number[][], tol: number = 1e-14): number {
  const { S } = SVD_double(A, { tol });
  let rank = 0;
  for (const s of S) {
    if (s > tol) rank++;
  }
  return rank;
}

/**
 * Compute the condition number of a matrix using SVD.
 *
 * The condition number is sigma_max / sigma_min, where sigma_max and sigma_min
 * are the largest and smallest singular values.
 *
 * @param A - A real matrix as a 2D array of numbers
 * @param tol - Tolerance for small singular values (default: 1e-14)
 * @returns The condition number (Infinity if matrix is rank-deficient)
 */
export function condition_number_SVD(A: number[][], tol: number = 1e-14): number {
  const { S } = SVD_double(A, { tol });
  if (S.length === 0) return Number.POSITIVE_INFINITY;
  const maxS = S[0]!;
  const minS = S[S.length - 1]!;
  if (minS < tol) return Number.POSITIVE_INFINITY;
  return maxS / minS;
}

/**
 * Compute the Frobenius norm of a matrix using SVD.
 *
 * The Frobenius norm is sqrt(sum of squared singular values).
 *
 * @param A - A real matrix as a 2D array of numbers
 * @returns The Frobenius norm
 */
export function frobenius_norm_SVD(A: number[][]): number {
  const { S } = SVD_double(A);
  let sum = 0;
  for (const s of S) {
    sum += s * s;
  }
  return Math.sqrt(sum);
}

/**
 * Compute the 2-norm (spectral norm) of a matrix using SVD.
 *
 * The 2-norm is the largest singular value.
 *
 * @param A - A real matrix as a 2D array of numbers
 * @returns The spectral norm
 */
export function spectral_norm_SVD(A: number[][]): number {
  const { S } = SVD_double(A);
  return S.length > 0 ? S[0]! : 0;
}

/**
 * Low-rank approximation of a matrix using truncated SVD.
 *
 * Computes the best rank-k approximation in the Frobenius norm.
 *
 * @param A - A real matrix as a 2D array of numbers
 * @param k - Target rank
 * @returns The rank-k approximation
 */
export function low_rank_approx_SVD(A: number[][], k: number): number[][] {
  const { U, S, V } = SVD_double(A);
  const m = U.length;
  const n = V.length;
  const kActual = Math.min(k, S.length);

  const Ak: number[][] = [];
  for (let i = 0; i < m; i++) {
    Ak.push([]);
    for (let j = 0; j < n; j++) {
      let sum = 0;
      for (let l = 0; l < kActual; l++) {
        sum += U[i]![l]! * S[l]! * V[j]![l]!;
      }
      Ak[i]!.push(sum);
    }
  }
  return Ak;
}

// ============================================================================
// QR Decomposition for Real Matrices with Householder Reflections
// ============================================================================

/**
 * Result type for QR decomposition.
 */
export interface QRResult {
  /** Orthogonal matrix Q (m x m for full, m x k for reduced) */
  Q: number[][];
  /** Upper triangular matrix R (m x n for full, k x n for reduced) */
  R: number[][];
}

/**
 * Return the QR decomposition of a real matrix using Householder reflections.
 *
 * Computes matrices Q, R such that A = Q * R where:
 * - Q is an m x m orthogonal matrix (or m x k for reduced)
 * - R is an m x n upper triangular matrix (or k x n for reduced)
 *
 * This implementation uses Householder reflections, which is numerically
 * more stable than Gram-Schmidt orthogonalization.
 *
 * @param A - A real matrix as a 2D array of numbers (IEEE 754 doubles)
 * @param reduced - If true, return reduced QR (Q is m x k, R is k x n where k = min(m,n))
 * @returns QR decomposition {Q, R}
 * @see Deviation: Uses IEEE 754 double precision instead of arbitrary precision.
 *   See DEVIATIONS.md for details.
 * @see Reference: sage/matrix/matrix_double_dense.pyx:QR
 *
 * @example
 * ```typescript
 * const A = [[12, -51, 4], [6, 167, -68], [-4, 24, -41]];
 * const { Q, R } = QR_double(A);
 * // Q is 3x3 orthogonal: Q * Q^T = I
 * // R is 3x3 upper triangular
 * // A = Q * R
 * ```
 */
export function QR_double(A: number[][], reduced: boolean = false): QRResult {
  const m = A.length;
  if (m === 0) {
    return { Q: [], R: [] };
  }

  const n = A[0]!.length;
  if (n === 0) {
    return { Q: Array.from({ length: m }, () => []), R: A.map((row) => [...row]) };
  }

  // Copy A to working matrix R
  const R: number[][] = A.map((row) => [...row]);

  // Initialize Q as identity (m x m)
  const Q: number[][] = [];
  for (let i = 0; i < m; i++) {
    Q.push([]);
    for (let j = 0; j < m; j++) {
      Q[i]!.push(i === j ? 1 : 0);
    }
  }

  const k = Math.min(m, n);

  for (let col = 0; col < k; col++) {
    // Extract column vector from R[col:m, col]
    const x: number[] = [];
    for (let i = col; i < m; i++) {
      x.push(R[i]![col]!);
    }

    // Compute norm of x
    let normX = 0;
    for (let i = 0; i < x.length; i++) {
      normX += x[i]! * x[i]!;
    }
    normX = Math.sqrt(normX);

    if (normX < 1e-15) continue;

    // Compute Householder vector v = x - sign(x[0]) * ||x|| * e_1
    const sign = x[0]! >= 0 ? 1 : -1;
    const v = [...x];
    v[0] = v[0]! + sign * normX;

    // Compute beta = 2 / (v^T * v)
    let vTv = 0;
    for (let i = 0; i < v.length; i++) {
      vTv += v[i]! * v[i]!;
    }
    if (vTv < 1e-30) continue;
    const beta = 2 / vTv;

    // Apply Householder reflection to R[col:m, col:n]
    // R = R - beta * v * (v^T * R)
    for (let j = col; j < n; j++) {
      let vTRj = 0;
      for (let i = 0; i < v.length; i++) {
        vTRj += v[i]! * R[col + i]![j]!;
      }
      for (let i = 0; i < v.length; i++) {
        R[col + i]![j] = R[col + i]![j]! - beta * v[i]! * vTRj;
      }
    }

    // Apply Householder reflection to Q[:, col:m]
    // Q = Q - beta * (Q * v) * v^T
    for (let i = 0; i < m; i++) {
      let Qiv = 0;
      for (let j = 0; j < v.length; j++) {
        Qiv += Q[i]![col + j]! * v[j]!;
      }
      for (let j = 0; j < v.length; j++) {
        Q[i]![col + j] = Q[i]![col + j]! - beta * Qiv * v[j]!;
      }
    }
  }

  // Clean up small values in R below diagonal
  for (let i = 0; i < m; i++) {
    for (let j = 0; j < Math.min(i, n); j++) {
      if (Math.abs(R[i]![j]!) < 1e-14) {
        R[i]![j] = 0;
      }
    }
  }

  if (reduced) {
    // Return reduced QR: Q is m x k, R is k x n
    const reducedQ: number[][] = [];
    for (let i = 0; i < m; i++) {
      reducedQ.push([]);
      for (let j = 0; j < k; j++) {
        reducedQ[i]!.push(Q[i]![j]!);
      }
    }

    const reducedR: number[][] = [];
    for (let i = 0; i < k; i++) {
      reducedR.push([]);
      for (let j = 0; j < n; j++) {
        reducedR[i]!.push(R[i]![j]!);
      }
    }

    return { Q: reducedQ, R: reducedR };
  }

  return { Q, R };
}

/**
 * Solve a linear system Ax = b using QR decomposition.
 *
 * More numerically stable than direct methods for ill-conditioned systems.
 *
 * @param A - Coefficient matrix (m x n)
 * @param b - Right-hand side vector (length m)
 * @returns Solution vector x (length n), or least-squares solution if m > n
 */
export function solve_QR(A: number[][], b: number[]): number[] {
  const { Q, R } = QR_double(A);
  const m = Q.length;
  const n = R[0]?.length || 0;

  // Compute Q^T * b
  const QTb: number[] = [];
  for (let j = 0; j < m; j++) {
    let sum = 0;
    for (let i = 0; i < m; i++) {
      sum += Q[i]![j]! * b[i]!;
    }
    QTb.push(sum);
  }

  // Back substitution on R * x = Q^T * b
  const x: number[] = new Array(n).fill(0);
  const k = Math.min(m, n);

  for (let i = k - 1; i >= 0; i--) {
    let sum = QTb[i]!;
    for (let j = i + 1; j < n; j++) {
      sum -= R[i]![j]! * x[j]!;
    }
    if (Math.abs(R[i]![i]!) > 1e-14) {
      x[i] = sum / R[i]![i]!;
    } else {
      x[i] = 0;
    }
  }

  return x;
}

// ============================================================================
// LU Decomposition for Real Matrices with Partial Pivoting
// ============================================================================

/**
 * Result type for LU decomposition with partial pivoting.
 */
export interface LUResult {
  /** Permutation vector (P[i] = j means row i of PA is row j of A) */
  P: number[];
  /** Lower triangular matrix with 1s on diagonal (m x min(m,n)) */
  L: number[][];
  /** Upper triangular matrix (min(m,n) x n) */
  U: number[][];
}

/**
 * Return the LU decomposition of a real matrix with partial pivoting.
 *
 * Computes P, L, U such that PA = LU where:
 * - P is represented as a permutation vector
 * - L is lower triangular with 1s on the diagonal
 * - U is upper triangular
 *
 * Uses partial pivoting for numerical stability.
 *
 * @param A - A real matrix as a 2D array of numbers (IEEE 754 doubles)
 * @returns LU decomposition {P, L, U}
 * @see Deviation: Uses IEEE 754 double precision instead of arbitrary precision.
 *   See DEVIATIONS.md for details.
 * @see Reference: sage/matrix/matrix2.pyx:LU
 *
 * @example
 * ```typescript
 * const A = [[2, 1, 1], [4, 3, 3], [8, 7, 9]];
 * const { P, L, U } = LU_double(A);
 * // L is lower triangular with 1s on diagonal
 * // U is upper triangular
 * // PA = LU where P permutes rows
 * ```
 */
export function LU_double(A: number[][]): LUResult {
  const m = A.length;
  if (m === 0) {
    return { P: [], L: [], U: [] };
  }

  const n = A[0]!.length;
  if (n === 0) {
    return { P: Array.from({ length: m }, (_, i) => i), L: A.map(() => []), U: [] };
  }

  // Copy A to working matrix U
  const U: number[][] = A.map((row) => [...row]);

  // Initialize L as identity (m x min(m,n))
  const k = Math.min(m, n);
  const L: number[][] = [];
  for (let i = 0; i < m; i++) {
    L.push([]);
    for (let j = 0; j < k; j++) {
      L[i]!.push(i === j ? 1 : 0);
    }
  }

  // Initialize permutation as identity
  const P: number[] = Array.from({ length: m }, (_, i) => i);

  for (let col = 0; col < k; col++) {
    // Partial pivoting: find the row with largest absolute value in column col
    let maxVal = Math.abs(U[col]![col]!);
    let maxRow = col;

    for (let i = col + 1; i < m; i++) {
      const absVal = Math.abs(U[i]![col]!);
      if (absVal > maxVal) {
        maxVal = absVal;
        maxRow = i;
      }
    }

    // Swap rows in U, L, and P
    if (maxRow !== col) {
      [U[col], U[maxRow]] = [U[maxRow]!, U[col]!];
      [P[col], P[maxRow]] = [P[maxRow]!, P[col]!];

      // Swap L rows up to column col-1
      for (let j = 0; j < col; j++) {
        [L[col]![j], L[maxRow]![j]] = [L[maxRow]![j]!, L[col]![j]!];
      }
    }

    const pivot = U[col]![col]!;
    if (Math.abs(pivot) < 1e-15) {
      // Skip this column if pivot is too small
      continue;
    }

    // Elimination
    for (let i = col + 1; i < m; i++) {
      const factor = U[i]![col]! / pivot;
      L[i]![col] = factor;
      U[i]![col] = 0;

      for (let j = col + 1; j < n; j++) {
        U[i]![j] = U[i]![j]! - factor * U[col]![j]!;
      }
    }
  }

  return { P, L, U };
}

/**
 * Solve a linear system Ax = b using LU decomposition.
 *
 * First solves Ly = Pb (forward substitution), then Ux = y (back substitution).
 *
 * @param A - Coefficient matrix (n x n, must be square)
 * @param b - Right-hand side vector (length n)
 * @returns Solution vector x (length n)
 */
export function solve_LU(A: number[][], b: number[]): number[] {
  const { P, L, U } = LU_double(A);
  const n = A.length;
  if (n === 0) return [];

  // Apply permutation to b: Pb
  const Pb: number[] = [];
  for (let i = 0; i < n; i++) {
    Pb.push(b[P[i]!]!);
  }

  // Forward substitution: Ly = Pb
  const y: number[] = new Array(n).fill(0);
  for (let i = 0; i < n; i++) {
    let sum = Pb[i]!;
    for (let j = 0; j < i; j++) {
      sum -= L[i]![j]! * y[j]!;
    }
    y[i] = sum; // L[i][i] = 1
  }

  // Back substitution: Ux = y
  const x: number[] = new Array(n).fill(0);
  for (let i = n - 1; i >= 0; i--) {
    let sum = y[i]!;
    for (let j = i + 1; j < n; j++) {
      sum -= U[i]![j]! * x[j]!;
    }
    if (Math.abs(U[i]![i]!) > 1e-14) {
      x[i] = sum / U[i]![i]!;
    } else {
      x[i] = 0;
    }
  }

  return x;
}

/**
 * Compute the determinant using LU decomposition.
 *
 * det(A) = det(P)^{-1} * det(L) * det(U) = (-1)^{swaps} * prod(U[i][i])
 *
 * @param A - A square matrix
 * @returns The determinant
 */
export function det_LU(A: number[][]): number {
  const n = A.length;
  if (n === 0) return 1;
  if (A[0]!.length !== n) {
    throw new ArithmeticError('determinant requires a square matrix');
  }

  const { P, U } = LU_double(A);

  // Count inversions in P (parity of permutation)
  let swaps = 0;
  const visited = new Array(n).fill(false);
  for (let i = 0; i < n; i++) {
    if (!visited[i]) {
      let j = i;
      let cycleLength = 0;
      while (!visited[j]) {
        visited[j] = true;
        j = P[j]!;
        cycleLength++;
      }
      if (cycleLength > 1) {
        swaps += cycleLength - 1;
      }
    }
  }

  // det = (-1)^swaps * prod(U[i][i])
  let det = swaps % 2 === 0 ? 1 : -1;
  for (let i = 0; i < n; i++) {
    det *= U[i]![i]!;
  }

  return det;
}

/**
 * Compute the inverse of a matrix using LU decomposition.
 *
 * @param A - A square invertible matrix
 * @returns The inverse matrix
 */
export function inverse_LU(A: number[][]): number[][] {
  const n = A.length;
  if (n === 0) return [];
  if (A[0]!.length !== n) {
    throw new ArithmeticError('inverse requires a square matrix');
  }

  // Solve A * X = I column by column
  const result: number[][] = [];
  for (let i = 0; i < n; i++) {
    result.push([]);
  }

  for (let col = 0; col < n; col++) {
    // Create e_col (standard basis vector)
    const e: number[] = new Array(n).fill(0);
    e[col] = 1;

    // Solve A * x = e_col
    const x = solve_LU(A, e);

    // Store as column of result
    for (let i = 0; i < n; i++) {
      result[i]!.push(x[i]!);
    }
  }

  return result;
}
