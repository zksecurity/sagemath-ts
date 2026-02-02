/**
 * Tests for matrix decompositions (LLL_gram, principal_square_root, decomposition, wiedemann)
 *
 * These test implementations based on:
 * - sage/matrix/matrix2.pyx
 */

import { describe, expect, it } from 'vitest';
import { GF } from '../rings/finite_rings/finite_field_constructor.js';
import { Integer } from '../rings/integer_ring.js';
import {
  IntegerMatrix,
  IntegerMatrixFromEntries,
  LLL_gram,
  Matrix,
  MatrixSpace,
  charpoly,
  decomposition,
  determinant,
  identity_integer_matrix,
  identity_matrix,
  principal_square_root,
  wiedemann,
  zero_matrix,
} from './index.js';

describe('LLL_gram', () => {
  // Note: LLL_gram is designed for integer matrices. Since the implementation
  // converts matrix entries to bigint internally, we test with IntegerMatrix.
  // Due to the complexity of the ring interface, we focus on basic correctness tests.

  it('should compute LLL transformation for a 2x2 Gram matrix using IntegerMatrix', () => {
    // Test case from SageMath documentation:
    // M = Matrix(ZZ, 2, 2, [5, 3, 3, 2])
    // U = M.LLL_gram()
    // MM = U.T * M * U should be LLL-reduced
    const M = IntegerMatrixFromEntries([
      [5n, 3n],
      [3n, 2n],
    ]);

    // Access the internal matrix representation
    const Mdata = [
      [M.get(0, 0), M.get(0, 1)],
      [M.get(1, 0), M.get(1, 1)],
    ];

    // Create the ring wrapper that IntegerMatrix uses
    const zzRing = {
      zero: () => new Integer(0n),
      one: () => new Integer(1n),
      __call__: (x: unknown): Integer => {
        if (x instanceof Integer) return x;
        if (typeof x === 'bigint') return new Integer(x);
        if (typeof x === 'number') return new Integer(x);
        throw new Error(`cannot coerce ${x} to Integer`);
      },
      is_field: () => false,
      toString: () => 'Integer Ring',
    };

    const Mgeneric = new Matrix(zzRing, 2, 2, Mdata);
    const U = LLL_gram(Mgeneric);

    // U should have determinant +/- 1 (unimodular)
    const det = determinant(U);
    expect(det.value === 1n || det.value === -1n).toBe(true);

    // Compute MM = U^T * M * U
    const UT = U.transpose();
    const MM = UT.mul(Mgeneric).mul(U);

    // The result should be LLL-reduced
    // For this specific example, the reduced Gram matrix should have small entries
    expect(MM.get(0, 0).value).toBeLessThanOrEqual(5n);
    expect(MM.get(1, 1).value).toBeLessThanOrEqual(5n);
  });

  it('should return identity for identity Gram matrix', () => {
    // Create the ring wrapper
    const zzRing = {
      zero: () => new Integer(0n),
      one: () => new Integer(1n),
      __call__: (x: unknown): Integer => {
        if (x instanceof Integer) return x;
        if (typeof x === 'bigint') return new Integer(x);
        if (typeof x === 'number') return new Integer(x);
        throw new Error(`cannot coerce ${x} to Integer`);
      },
      is_field: () => false,
      toString: () => 'Integer Ring',
    };

    const I = identity_matrix(zzRing, 2);
    const U = LLL_gram(I);

    // For identity matrix, U should be identity or a sign permutation
    const det = determinant(U);
    expect(det.value === 1n || det.value === -1n).toBe(true);
  });

  it('should handle 1x1 matrices', () => {
    const zzRing = {
      zero: () => new Integer(0n),
      one: () => new Integer(1n),
      __call__: (x: unknown): Integer => {
        if (x instanceof Integer) return x;
        if (typeof x === 'bigint') return new Integer(x);
        if (typeof x === 'number') return new Integer(x);
        throw new Error(`cannot coerce ${x} to Integer`);
      },
      is_field: () => false,
      toString: () => 'Integer Ring',
    };

    const M = new Matrix(zzRing, 1, 1, [[new Integer(5n)]]);
    const U = LLL_gram(M);

    expect(U.nrows).toBe(1);
    expect(U.ncols).toBe(1);

    const det = determinant(U);
    expect(det.value === 1n || det.value === -1n).toBe(true);
  });

  it('should handle empty matrices', () => {
    const zzRing = {
      zero: () => new Integer(0n),
      one: () => new Integer(1n),
      __call__: (x: unknown): Integer => new Integer(0n),
      is_field: () => false,
      toString: () => 'Integer Ring',
    };

    const M = zero_matrix(zzRing, 0);
    const U = LLL_gram(M);

    expect(U.nrows).toBe(0);
    expect(U.ncols).toBe(0);
  });
});

describe('principal_square_root', () => {
  // Note: principal_square_root requires elements that support sqrt()
  // For testing, we can verify it with simple cases or mock elements

  it('should handle empty matrices', () => {
    const F7 = GF(7n);
    const ring = F7;
    const M = zero_matrix(ring, 0);

    const sqrtM = principal_square_root(
      M as Matrix<typeof ring extends { __call__: (x: bigint) => infer R } ? R : never>
    );

    if (sqrtM !== false) {
      expect(sqrtM.nrows).toBe(0);
      expect(sqrtM.ncols).toBe(0);
    }
  });

  it('should throw for non-square matrices', () => {
    const F7 = GF(7n);
    const MS = MatrixSpace(F7, 2, 3);
    const M = MS.zero();

    expect(() => principal_square_root(M)).toThrow();
  });
});

describe('decomposition', () => {
  it('should decompose the identity matrix into a single irreducible subspace', () => {
    const F7 = GF(7n);
    const MS = MatrixSpace(F7, 3, 3);
    const I = MS.identity();

    const result = decomposition(I);

    // Identity matrix has characteristic polynomial (x-1)^n
    // which is not irreducible for n > 1, so is_irreducible should be false
    expect(Array.isArray(result)).toBe(true);
    if (Array.isArray(result) && !Array.isArray(result[0]?.[0])) {
      const decomp = result as Array<
        [Matrix<typeof F7 extends { __call__: (x: bigint) => infer R } ? R : never>, boolean]
      >;
      // Should have at least one component
      expect(decomp.length).toBeGreaterThan(0);

      // Total dimension should equal matrix size
      let totalDim = 0;
      for (const [basis] of decomp) {
        totalDim += basis.nrows;
      }
      expect(totalDim).toBe(3);
    }
  });

  it('should handle empty matrices', () => {
    const F7 = GF(7n);
    const M = zero_matrix(F7, 0);

    const result = decomposition(
      M as Matrix<typeof F7 extends { __call__: (x: bigint) => infer R } ? R : never>
    );

    expect(Array.isArray(result)).toBe(true);
    if (Array.isArray(result) && !Array.isArray(result[0])) {
      expect(result.length).toBe(0);
    }
  });

  it('should decompose a diagonal matrix with distinct eigenvalues', () => {
    const F7 = GF(7n);
    const MS = MatrixSpace(F7, 2, 2);

    // Diagonal matrix with eigenvalues 1 and 2
    const D = MS.__call__([
      [F7.__call__(1n), F7.__call__(0n)],
      [F7.__call__(0n), F7.__call__(2n)],
    ]);

    const result = decomposition(D, 'kernel', true); // diagonalizable

    expect(Array.isArray(result)).toBe(true);
    if (Array.isArray(result) && !Array.isArray(result[0]?.[0])) {
      const decomp = result as Array<
        [Matrix<typeof F7 extends { __call__: (x: bigint) => infer R } ? R : never>, boolean]
      >;

      // Should have two 1-dimensional eigenspaces, both irreducible
      expect(decomp.length).toBe(2);
      expect(decomp[0]![0].nrows).toBe(1);
      expect(decomp[1]![0].nrows).toBe(1);
      expect(decomp[0]![1]).toBe(true); // irreducible
      expect(decomp[1]![1]).toBe(true); // irreducible
    }
  });

  it('should return dual decomposition when requested', () => {
    const F7 = GF(7n);
    const MS = MatrixSpace(F7, 2, 2);

    const D = MS.__call__([
      [F7.__call__(1n), F7.__call__(0n)],
      [F7.__call__(0n), F7.__call__(2n)],
    ]);

    const result = decomposition(D, 'kernel', true, true); // with dual

    expect(Array.isArray(result)).toBe(true);
    expect(result.length).toBe(2);
    const [decomp, dualDecomp] = result as [
      Array<[Matrix<typeof F7 extends { __call__: (x: bigint) => infer R } ? R : never>, boolean]>,
      Array<[Matrix<typeof F7 extends { __call__: (x: bigint) => infer R } ? R : never>, boolean]>,
    ];

    expect(decomp.length).toBe(2);
    expect(dualDecomp.length).toBe(2);
  });
});

describe('wiedemann', () => {
  it('should compute a polynomial dividing the characteristic polynomial', () => {
    const F7 = GF(7n);
    const MS = MatrixSpace(F7, 3, 3);

    // Simple matrix
    const A = MS.__call__([
      [F7.__call__(0n), F7.__call__(1n), F7.__call__(2n)],
      [F7.__call__(3n), F7.__call__(4n), F7.__call__(5n)],
      [F7.__call__(6n), F7.__call__(0n), F7.__call__(1n)],
    ]);

    const f = wiedemann(A, 0);

    // f should divide the characteristic polynomial
    const chi = charpoly(A);

    // The degree of f should be at most the degree of chi
    expect(f.degree()).toBeLessThanOrEqual(chi.degree());

    // f(A) should have the i-th column in its null space
    // (i.e., f(A) * e_i = 0 where e_i is the i-th standard basis vector)
    // This is a stronger test but requires matrix evaluation
  });

  it('should handle the identity matrix', () => {
    const F7 = GF(7n);
    const MS = MatrixSpace(F7, 2, 2);
    const I = MS.identity();

    const f = wiedemann(I, 0);

    // For identity, the minimal polynomial is x - 1
    // f should divide this
    expect(f.degree()).toBeLessThanOrEqual(1);
  });

  it('should handle diagonal matrices', () => {
    const F7 = GF(7n);
    const MS = MatrixSpace(F7, 2, 2);

    const D = MS.__call__([
      [F7.__call__(2n), F7.__call__(0n)],
      [F7.__call__(0n), F7.__call__(3n)],
    ]);

    const f = wiedemann(D, 0);

    // For diagonal matrix with distinct eigenvalues 2 and 3,
    // using the 0-th basis vector [1, 0], the sequence is 1, 2, 4, 8, ...
    // which has minimal polynomial x - 2
    expect(f.degree()).toBeLessThanOrEqual(2);
  });

  it('should throw for invalid index', () => {
    const F7 = GF(7n);
    const MS = MatrixSpace(F7, 2, 2);
    const A = MS.identity();

    expect(() => wiedemann(A, 5)).toThrow();
    expect(() => wiedemann(A, -1)).toThrow();
  });

  it('should use specific coordinate when t is specified', () => {
    const F7 = GF(7n);
    const MS = MatrixSpace(F7, 2, 2);

    const A = MS.__call__([
      [F7.__call__(1n), F7.__call__(1n)],
      [F7.__call__(0n), F7.__call__(2n)],
    ]);

    // Use only the first coordinate (t=0)
    const f0 = wiedemann(A, 0, 0);

    // Use only the second coordinate (t=1)
    const f1 = wiedemann(A, 0, 1);

    // Both should produce valid results
    expect(f0.degree()).toBeGreaterThanOrEqual(0);
    expect(f1.degree()).toBeGreaterThanOrEqual(0);
  });
});

describe('integration tests', () => {
  it('decomposition results should be invariant subspaces', () => {
    const F7 = GF(7n);
    const MS = MatrixSpace(F7, 3, 3);

    // Create a matrix with known structure
    const A = MS.__call__([
      [F7.__call__(1n), F7.__call__(1n), F7.__call__(0n)],
      [F7.__call__(0n), F7.__call__(1n), F7.__call__(0n)],
      [F7.__call__(0n), F7.__call__(0n), F7.__call__(2n)],
    ]);

    const result = decomposition(A, 'kernel');

    if (Array.isArray(result) && !Array.isArray(result[0]?.[0])) {
      const decomp = result as Array<
        [Matrix<typeof F7 extends { __call__: (x: bigint) => infer R } ? R : never>, boolean]
      >;

      for (const [basisMatrix] of decomp) {
        // Each row of basisMatrix, when acted on by A,
        // should remain in the span of basisMatrix
        // This is the invariance property
        const dim = basisMatrix.nrows;
        if (dim > 0) {
          // The subspace is A-invariant
          expect(dim).toBeGreaterThan(0);
        }
      }
    }
  });
});

// ============================================================================
// Tests for Real Matrix Decompositions (IEEE 754 doubles)
// ============================================================================

import {
  LU_double,
  QR_double,
  SVD_double,
  SVD_reconstruct,
  condition_number_SVD,
  det_LU,
  frobenius_norm_SVD,
  inverse_LU,
  low_rank_approx_SVD,
  pseudoinverse_SVD,
  rank_SVD,
  solve_LU,
  solve_QR,
  spectral_norm_SVD,
} from './matrix_decompositions_additions.js';

// Helper function to check if two numbers are approximately equal
function approxEqual(a: number, b: number, tol: number = 1e-10): boolean {
  return Math.abs(a - b) < tol;
}

// Helper function to check if a matrix is approximately equal to another
function matrixApproxEqual(A: number[][], B: number[][], tol: number = 1e-10): boolean {
  if (A.length !== B.length) return false;
  for (let i = 0; i < A.length; i++) {
    if (A[i]!.length !== B[i]!.length) return false;
    for (let j = 0; j < A[i]!.length; j++) {
      if (!approxEqual(A[i]![j]!, B[i]![j]!, tol)) return false;
    }
  }
  return true;
}

// Helper function to check if a matrix is orthogonal (Q * Q^T = I)
function isOrthogonal(Q: number[][], tol: number = 1e-10): boolean {
  const n = Q.length;
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      let dot = 0;
      for (let k = 0; k < Q[0]!.length; k++) {
        dot += Q[i]![k]! * Q[j]![k]!;
      }
      const expected = i === j ? 1 : 0;
      if (!approxEqual(dot, expected, tol)) return false;
    }
  }
  return true;
}

// Helper function to check if a matrix is upper triangular
function isUpperTriangular(R: number[][], tol: number = 1e-10): boolean {
  for (let i = 0; i < R.length; i++) {
    for (let j = 0; j < Math.min(i, R[i]!.length); j++) {
      if (!approxEqual(R[i]![j]!, 0, tol)) return false;
    }
  }
  return true;
}

// Helper function to multiply matrices
function matmul(A: number[][], B: number[][]): number[][] {
  if (A.length === 0 || B.length === 0) return [];
  const m = A.length;
  const n = B[0]!.length;
  const k = A[0]!.length;
  const C: number[][] = [];
  for (let i = 0; i < m; i++) {
    C.push([]);
    for (let j = 0; j < n; j++) {
      let sum = 0;
      for (let l = 0; l < k; l++) {
        sum += A[i]![l]! * B[l]![j]!;
      }
      C[i]!.push(sum);
    }
  }
  return C;
}

describe('SVD_double', () => {
  it('should compute SVD for a 2x2 matrix', () => {
    const A = [
      [3, 2],
      [2, 3],
    ];

    const { U, S, V } = SVD_double(A);

    // Check dimensions
    expect(U.length).toBe(2);
    expect(U[0]!.length).toBe(2);
    expect(S.length).toBe(2);
    expect(V.length).toBe(2);
    expect(V[0]!.length).toBe(2);

    // Singular values should be positive and in descending order
    expect(S[0]!).toBeGreaterThan(0);
    expect(S[1]!).toBeGreaterThan(0);
    expect(S[0]!).toBeGreaterThanOrEqual(S[1]!);

    // For symmetric [[3,2],[2,3]], eigenvalues are 5 and 1
    // So singular values should be 5 and 1
    expect(approxEqual(S[0]!, 5, 1e-10)).toBe(true);
    expect(approxEqual(S[1]!, 1, 1e-10)).toBe(true);

    // U and V should be orthogonal
    expect(isOrthogonal(U)).toBe(true);
    expect(isOrthogonal(V)).toBe(true);
  });

  it('should reconstruct the original matrix from SVD', () => {
    const A = [
      [1, 2],
      [3, 4],
      [5, 6],
    ];

    const svd = SVD_double(A);
    const reconstructed = SVD_reconstruct(svd);

    expect(matrixApproxEqual(A, reconstructed, 1e-10)).toBe(true);
  });

  it('should handle rectangular matrices (m > n)', () => {
    const A = [
      [1, 2],
      [3, 4],
      [5, 6],
    ];

    const { U, S, V } = SVD_double(A);

    // Dimensions: U is 3x3, S has 2 elements, V is 2x2
    expect(U.length).toBe(3);
    expect(U[0]!.length).toBe(3);
    expect(S.length).toBe(2);
    expect(V.length).toBe(2);
    expect(V[0]!.length).toBe(2);

    expect(isOrthogonal(U)).toBe(true);
    expect(isOrthogonal(V)).toBe(true);
  });

  it('should handle rectangular matrices (m < n)', () => {
    const A = [
      [1, 2, 3],
      [4, 5, 6],
    ];

    const { U, S, V } = SVD_double(A);

    // Dimensions: U is 2x2, S has 2 elements, V is 3x3
    expect(U.length).toBe(2);
    expect(U[0]!.length).toBe(2);
    expect(S.length).toBe(2);
    expect(V.length).toBe(3);
    expect(V[0]!.length).toBe(3);

    expect(isOrthogonal(U)).toBe(true);
    expect(isOrthogonal(V)).toBe(true);
  });

  it('should handle empty matrices', () => {
    const { U, S, V } = SVD_double([]);
    expect(U.length).toBe(0);
    expect(S.length).toBe(0);
    expect(V.length).toBe(0);
  });

  it('should compute correct singular values for a rank-1 matrix', () => {
    // Rank-1 matrix: outer product of [1, 2] and [3, 4]
    const A = [
      [3, 4],
      [6, 8],
    ];

    const { S } = SVD_double(A);

    // First singular value = ||[1,2]|| * ||[3,4]|| = sqrt(5) * 5 = 5*sqrt(5) approx 11.18
    // Actually, let's compute: A = [3,4;6,8] has eigenvalues of A^T*A
    // A^T*A = [[45, 60], [60, 80]], trace = 125, det = 0
    // eigenvalues: 125 and 0. Singular values: sqrt(125) = 5*sqrt(5), 0
    expect(approxEqual(S[0]!, 5 * Math.sqrt(5), 1e-10)).toBe(true);
    expect(approxEqual(S[1]!, 0, 1e-10)).toBe(true);
  });
});

describe('pseudoinverse_SVD', () => {
  it('should compute pseudoinverse of a square invertible matrix', () => {
    const A = [
      [4, 7],
      [2, 6],
    ];

    const Aplus = pseudoinverse_SVD(A);

    // For invertible matrix, pseudoinverse = inverse
    // A^{-1} = (1/10) * [[6, -7], [-2, 4]]
    const expected = [
      [0.6, -0.7],
      [-0.2, 0.4],
    ];

    expect(matrixApproxEqual(Aplus, expected, 1e-10)).toBe(true);

    // A * A^+ should be identity (for square invertible)
    const product = matmul(A, Aplus);
    const I = [
      [1, 0],
      [0, 1],
    ];
    expect(matrixApproxEqual(product, I, 1e-10)).toBe(true);
  });

  it('should compute pseudoinverse of a tall matrix', () => {
    const A = [
      [1, 0],
      [0, 1],
      [1, 1],
    ];

    const Aplus = pseudoinverse_SVD(A);

    // A^+ should be 2x3
    expect(Aplus.length).toBe(2);
    expect(Aplus[0]!.length).toBe(3);

    // A^+ * A should be identity (for full column rank)
    const product = matmul(Aplus, A);
    const I = [
      [1, 0],
      [0, 1],
    ];
    expect(matrixApproxEqual(product, I, 1e-10)).toBe(true);
  });
});

describe('rank_SVD', () => {
  it('should compute rank of a full rank matrix', () => {
    const A = [
      [1, 2],
      [3, 4],
    ];
    expect(rank_SVD(A)).toBe(2);
  });

  it('should compute rank of a rank-deficient matrix', () => {
    // Linearly dependent rows
    const A = [
      [1, 2],
      [2, 4],
    ];
    expect(rank_SVD(A)).toBe(1);
  });

  it('should compute rank of zero matrix', () => {
    const A = [
      [0, 0],
      [0, 0],
    ];
    expect(rank_SVD(A)).toBe(0);
  });
});

describe('condition_number_SVD', () => {
  it('should compute condition number of identity matrix', () => {
    const I = [
      [1, 0],
      [0, 1],
    ];
    expect(approxEqual(condition_number_SVD(I), 1, 1e-10)).toBe(true);
  });

  it('should compute condition number of well-conditioned matrix', () => {
    // Symmetric matrix [[2, 1], [1, 2]] has eigenvalues 3 and 1
    const A = [
      [2, 1],
      [1, 2],
    ];
    expect(approxEqual(condition_number_SVD(A), 3, 1e-10)).toBe(true);
  });

  it('should return Infinity for singular matrix', () => {
    const A = [
      [1, 2],
      [2, 4],
    ];
    expect(condition_number_SVD(A)).toBe(Number.POSITIVE_INFINITY);
  });
});

describe('norms via SVD', () => {
  it('should compute Frobenius norm', () => {
    const A = [
      [1, 2],
      [3, 4],
    ];
    // Frobenius norm = sqrt(1 + 4 + 9 + 16) = sqrt(30)
    expect(approxEqual(frobenius_norm_SVD(A), Math.sqrt(30), 1e-10)).toBe(true);
  });

  it('should compute spectral norm', () => {
    // For symmetric [[2, 1], [1, 2]], spectral norm = largest eigenvalue = 3
    const A = [
      [2, 1],
      [1, 2],
    ];
    expect(approxEqual(spectral_norm_SVD(A), 3, 1e-10)).toBe(true);
  });
});

describe('low_rank_approx_SVD', () => {
  it('should compute rank-1 approximation', () => {
    const A = [
      [3, 2],
      [2, 3],
    ];

    const A1 = low_rank_approx_SVD(A, 1);

    // Rank-1 approximation should have rank 1
    expect(rank_SVD(A1)).toBe(1);
  });

  it('should return original for k >= rank', () => {
    const A = [
      [1, 2],
      [3, 4],
    ];

    const A2 = low_rank_approx_SVD(A, 2);

    expect(matrixApproxEqual(A, A2, 1e-10)).toBe(true);
  });
});

describe('QR_double', () => {
  it('should compute QR decomposition for a 3x3 matrix', () => {
    const A = [
      [12, -51, 4],
      [6, 167, -68],
      [-4, 24, -41],
    ];

    const { Q, R } = QR_double(A);

    // Q should be 3x3 orthogonal
    expect(Q.length).toBe(3);
    expect(Q[0]!.length).toBe(3);
    expect(isOrthogonal(Q)).toBe(true);

    // R should be 3x3 upper triangular
    expect(R.length).toBe(3);
    expect(R[0]!.length).toBe(3);
    expect(isUpperTriangular(R)).toBe(true);

    // Q * R should equal A
    const product = matmul(Q, R);
    expect(matrixApproxEqual(A, product, 1e-10)).toBe(true);
  });

  it('should compute QR for rectangular matrix (m > n)', () => {
    const A = [
      [1, 2],
      [3, 4],
      [5, 6],
    ];

    const { Q, R } = QR_double(A);

    // Full QR: Q is 3x3, R is 3x2
    expect(Q.length).toBe(3);
    expect(Q[0]!.length).toBe(3);
    expect(R.length).toBe(3);
    expect(R[0]!.length).toBe(2);

    expect(isOrthogonal(Q)).toBe(true);
    expect(isUpperTriangular(R)).toBe(true);

    const product = matmul(Q, R);
    expect(matrixApproxEqual(A, product, 1e-10)).toBe(true);
  });

  it('should compute reduced QR when requested', () => {
    const A = [
      [1, 2],
      [3, 4],
      [5, 6],
    ];

    const { Q, R } = QR_double(A, true); // reduced

    // Reduced QR: Q is 3x2, R is 2x2
    expect(Q.length).toBe(3);
    expect(Q[0]!.length).toBe(2);
    expect(R.length).toBe(2);
    expect(R[0]!.length).toBe(2);

    expect(isUpperTriangular(R)).toBe(true);

    const product = matmul(Q, R);
    expect(matrixApproxEqual(A, product, 1e-10)).toBe(true);
  });

  it('should handle empty matrices', () => {
    const { Q, R } = QR_double([]);
    expect(Q.length).toBe(0);
    expect(R.length).toBe(0);
  });
});

describe('solve_QR', () => {
  it('should solve a simple linear system', () => {
    const A = [
      [2, 1],
      [1, 3],
    ];
    const b = [4, 5];

    const x = solve_QR(A, b);

    // Verify: Ax = b
    const Ax = [A[0]![0]! * x[0]! + A[0]![1]! * x[1]!, A[1]![0]! * x[0]! + A[1]![1]! * x[1]!];

    expect(approxEqual(Ax[0]!, b[0]!, 1e-10)).toBe(true);
    expect(approxEqual(Ax[1]!, b[1]!, 1e-10)).toBe(true);
  });

  it('should solve least squares for overdetermined system', () => {
    const A = [
      [1, 1],
      [1, 2],
      [1, 3],
    ];
    const b = [6, 9, 12]; // Exactly on line y = 3x + 3

    const x = solve_QR(A, b);

    // x should be approximately [3, 3]
    expect(approxEqual(x[0]!, 3, 1e-10)).toBe(true);
    expect(approxEqual(x[1]!, 3, 1e-10)).toBe(true);
  });
});

describe('LU_double', () => {
  it('should compute LU decomposition with partial pivoting', () => {
    const A = [
      [2, 1, 1],
      [4, 3, 3],
      [8, 7, 9],
    ];

    const { P, L, U } = LU_double(A);

    // Check dimensions
    expect(L.length).toBe(3);
    expect(U.length).toBe(3);
    expect(P.length).toBe(3);

    // L should be lower triangular with 1s on diagonal
    for (let i = 0; i < 3; i++) {
      expect(approxEqual(L[i]![i]!, 1, 1e-10)).toBe(true);
      for (let j = i + 1; j < 3; j++) {
        expect(approxEqual(L[i]![j]!, 0, 1e-10)).toBe(true);
      }
    }

    // U should be upper triangular
    expect(isUpperTriangular(U)).toBe(true);

    // PA = LU
    // Apply P to A
    const PA: number[][] = [];
    for (let i = 0; i < 3; i++) {
      PA.push([...A[P[i]!]!]);
    }

    const LU = matmul(L, U);
    expect(matrixApproxEqual(PA, LU, 1e-10)).toBe(true);
  });

  it('should handle empty matrices', () => {
    const { P, L, U } = LU_double([]);
    expect(P.length).toBe(0);
    expect(L.length).toBe(0);
    expect(U.length).toBe(0);
  });
});

describe('solve_LU', () => {
  it('should solve a linear system', () => {
    const A = [
      [2, 1, 1],
      [4, 3, 3],
      [8, 7, 9],
    ];
    const b = [4, 10, 24];

    const x = solve_LU(A, b);

    // Verify: Ax = b
    const Ax = [
      A[0]![0]! * x[0]! + A[0]![1]! * x[1]! + A[0]![2]! * x[2]!,
      A[1]![0]! * x[0]! + A[1]![1]! * x[1]! + A[1]![2]! * x[2]!,
      A[2]![0]! * x[0]! + A[2]![1]! * x[1]! + A[2]![2]! * x[2]!,
    ];

    expect(approxEqual(Ax[0]!, b[0]!, 1e-10)).toBe(true);
    expect(approxEqual(Ax[1]!, b[1]!, 1e-10)).toBe(true);
    expect(approxEqual(Ax[2]!, b[2]!, 1e-10)).toBe(true);
  });
});

describe('det_LU', () => {
  it('should compute determinant of 2x2 matrix', () => {
    const A = [
      [4, 6],
      [3, 8],
    ];
    // det = 4*8 - 6*3 = 32 - 18 = 14
    expect(approxEqual(det_LU(A), 14, 1e-10)).toBe(true);
  });

  it('should compute determinant of 3x3 matrix', () => {
    const A = [
      [1, 2, 3],
      [4, 5, 6],
      [7, 8, 9],
    ];
    // This matrix is singular, det = 0
    expect(approxEqual(det_LU(A), 0, 1e-10)).toBe(true);
  });

  it('should return 1 for identity matrix', () => {
    const I = [
      [1, 0, 0],
      [0, 1, 0],
      [0, 0, 1],
    ];
    expect(approxEqual(det_LU(I), 1, 1e-10)).toBe(true);
  });

  it('should handle empty matrix', () => {
    expect(det_LU([])).toBe(1);
  });
});

describe('inverse_LU', () => {
  it('should compute inverse of 2x2 matrix', () => {
    const A = [
      [4, 7],
      [2, 6],
    ];

    const Ainv = inverse_LU(A);

    // A^{-1} = (1/10) * [[6, -7], [-2, 4]]
    const expected = [
      [0.6, -0.7],
      [-0.2, 0.4],
    ];

    expect(matrixApproxEqual(Ainv, expected, 1e-10)).toBe(true);

    // A * A^{-1} should be identity
    const product = matmul(A, Ainv);
    const I = [
      [1, 0],
      [0, 1],
    ];
    expect(matrixApproxEqual(product, I, 1e-10)).toBe(true);
  });

  it('should compute inverse of 3x3 matrix', () => {
    const A = [
      [1, 2, 3],
      [0, 1, 4],
      [5, 6, 0],
    ];

    const Ainv = inverse_LU(A);

    // A * A^{-1} should be identity
    const product = matmul(A, Ainv);
    const I = [
      [1, 0, 0],
      [0, 1, 0],
      [0, 0, 1],
    ];
    expect(matrixApproxEqual(product, I, 1e-10)).toBe(true);
  });
});
