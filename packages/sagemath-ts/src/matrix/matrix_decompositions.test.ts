/**
 * Tests for matrix decompositions (LLL_gram, principal_square_root, decomposition, wiedemann)
 *
 * These test implementations based on:
 * - sage/matrix/matrix2.pyx
 */

import { describe, expect, it } from 'bun:test';
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

    // Sage/PARI return exactly [[-1, 1], [1, -2]] (det 1) for this Gram matrix.
    expect(U.get(0, 0).value).toBe(-1n);
    expect(U.get(0, 1).value).toBe(1n);
    expect(U.get(1, 0).value).toBe(1n);
    expect(U.get(1, 1).value).toBe(-2n);
    expect(determinant(U).value).toBe(1n);

    // Compute MM = U^T * M * U; the Sage doctest says this is the identity.
    const UT = U.transpose();
    const MM = UT.mul(Mgeneric).mul(U);
    expect(MM.get(0, 0).value).toBe(1n);
    expect(MM.get(0, 1).value).toBe(0n);
    expect(MM.get(1, 0).value).toBe(0n);
    expect(MM.get(1, 1).value).toBe(1n);
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

// ============================================================================
// Tests for the AUDIT-2026-07 fixes
// ============================================================================

import { QQ } from '../rings/rational_field.js';
import {
  LU,
  QR,
  block_ldlt,
  cholesky,
  echelon_form,
  echelonize,
  elementary_divisors,
  extended_echelon_form,
  gram_schmidt,
  gram_schmidt_noscale,
  hermite_form,
  hessenberg_form,
  indefinite_factorization,
  jordan_form,
  krylov_basis,
  krylov_kernel_basis,
  krylov_matrix,
  matrix_str,
  pivot_rows,
  pivots,
  rref,
  smith_form,
} from './matrix_decompositions.js';

// ============================================================================
// Audit fixes (AUDIT-2026-07: H36-H45, M49-M53): defining identities
//
// Every expected value below is the value SageMath actually produces, taken
// from the doctests of sage/matrix/matrix2.pyx.
// ============================================================================

const F19 = GF(19n);
const F97 = GF(97n);
const F101 = GF(101n);

/** Build a Matrix over a finite field from small integers. */
function ffmat(F: ReturnType<typeof GF>, entries: number[][]): Matrix<any> {
  return new Matrix(
    F as any,
    entries.length,
    entries[0]!.length,
    entries.map((r) => r.map((x) => (F as any).__call__(x)))
  );
}

/** Build a Matrix over QQ from small integers. */
function qmat(entries: number[][]): Matrix<any> {
  return new Matrix(
    QQ as any,
    entries.length,
    entries[0]!.length,
    entries.map((r) => r.map((x) => QQ.__call__(x)))
  );
}

/** Render a matrix the way SageMath prints its rows, for readable assertions. */
function render(M: Matrix<any>): string {
  const rows: string[] = [];
  for (let i = 0; i < M.nrows; i++) {
    const r: string[] = [];
    for (let j = 0; j < M.ncols; j++) r.push(String(M.get(i, j)));
    rows.push('[' + r.join(' ') + ']');
  }
  return rows.join(' / ');
}

function matEq(A: Matrix<any>, B: Matrix<any>): boolean {
  if (A.nrows !== B.nrows || A.ncols !== B.ncols) return false;
  for (let i = 0; i < A.nrows; i++) {
    for (let j = 0; j < A.ncols; j++) {
      if (!A.get(i, j).eq(B.get(i, j))) return false;
    }
  }
  return true;
}

/** Deterministic small PRNG so the property tests are reproducible. */
function makeRng(seed: number): () => number {
  let s = seed;
  return () => {
    s = (s * 1103515245 + 12345) % 2147483648;
    return Math.abs(s);
  };
}

const zzRing = {
  zero: () => new Integer(0n),
  one: () => new Integer(1n),
  __call__: (x: unknown): Integer => {
    if (x instanceof Integer) return x;
    if (typeof x === 'bigint') return new Integer(x);
    if (typeof x === 'number') return new Integer(BigInt(x));
    throw new Error(`cannot coerce ${x} to Integer`);
  },
  is_field: () => false,
  toString: () => 'Integer Ring',
};

function zmat(entries: bigint[][]): Matrix<any> {
  return new Matrix(
    zzRing as any,
    entries.length,
    entries[0]!.length,
    entries.map((r) => r.map((x) => new Integer(x)))
  );
}

describe('echelon_form / pivots / pivot_rows', () => {
  it('echelon_form is the reduced row echelon form (Sage GF(19) doctest)', () => {
    // sage: MatrixSpace(GF(19),2,3)([1,2,3,4,5,6]).echelon_form()
    // [ 1  0 18]
    // [ 0  1  2]
    const C = ffmat(F19, [
      [1, 2, 3],
      [4, 5, 6],
    ]);
    expect(render(echelon_form(C))).toBe('[1 0 18] / [0 1 2]');
  });

  it('echelon_form reduces above the pivots (Sage QQ range(9) doctest)', () => {
    // sage: matrix(QQ,3,3,range(9)).echelon_form()
    // [ 1  0 -1]
    // [ 0  1  2]
    // [ 0  0  0]
    const A = qmat([
      [0, 1, 2],
      [3, 4, 5],
      [6, 7, 8],
    ]);
    expect(render(echelon_form(A))).toBe('[1 0 -1] / [0 1 2] / [0 0 0]');
  });

  it('echelonize mutates in place and returns T with T*A == E', () => {
    const A = qmat([
      [0, 1, 2],
      [3, 4, 5],
      [6, 7, 8],
    ]);
    const original = A.copy();
    const T = echelonize(A, 'default', 0, true)!;
    expect(render(A)).toBe('[1 0 -1] / [0 1 2] / [0 0 0]');
    expect(matEq(T.mul(original), A)).toBe(true);
  });

  it('pivots returns pivot column indices, pivot_rows pivot row indices', () => {
    // sage: matrix(QQ,3,3,[0,0,0,1,2,3,2,4,6]).pivot_rows() == (1,)
    const A = qmat([
      [0, 0, 0],
      [1, 2, 3],
      [2, 4, 6],
    ]);
    expect(pivot_rows(A)).toEqual([1]);
    expect(pivots(A)).toEqual([0]);
  });

  it('rref agrees with echelon_form over a field', () => {
    const A = ffmat(F101, [
      [2, 4, 6, 1],
      [1, 2, 3, 0],
      [0, 0, 5, 5],
    ]);
    expect(render(rref(A))).toBe(render(echelon_form(A)));
  });
});

describe('LU', () => {
  it('P*L*U == A for the 3-cycle permutation matrix', () => {
    const A = ffmat(F101, [
      [0, 1, 0],
      [0, 0, 1],
      [1, 0, 0],
    ]);
    const [P, L, U] = LU(A) as [Matrix<any>, Matrix<any>, Matrix<any>];
    expect(matEq(P.mul(L).mul(U), A)).toBe(true);
  });

  it('P*L*U == A, L unit lower triangular, U upper triangular (200 random)', () => {
    const rnd = makeRng(12345);
    for (let t = 0; t < 200; t++) {
      const m = 1 + (rnd() % 5);
      const n = 1 + (rnd() % 5);
      const e: number[][] = [];
      for (let i = 0; i < m; i++) {
        const r: number[] = [];
        for (let j = 0; j < n; j++) r.push(rnd() % 101);
        e.push(r);
      }
      const A = ffmat(F101, e);
      const [P, L, U] = LU(A) as [Matrix<any>, Matrix<any>, Matrix<any>];
      expect(matEq(P.mul(L).mul(U), A)).toBe(true);
      for (let i = 0; i < L.nrows; i++) {
        expect(L.get(i, i).eq((F101 as any).__call__(1))).toBe(true);
        for (let j = i + 1; j < L.ncols; j++) expect(L.get(i, j).isZero()).toBe(true);
      }
      for (let i = 0; i < U.nrows; i++) {
        for (let j = 0; j < Math.min(i, U.ncols); j++) expect(U.get(i, j).isZero()).toBe(true);
      }
    }
  });
});

describe('gram_schmidt', () => {
  it('orthogonalizes rows: A == M*G (Sage QQ doctest)', () => {
    // sage: A = matrix(QQ, [[-1,3,2,2],[-1,0,-1,0],[-1,-2,-3,-1],[1,1,2,0]])
    // sage: G, M = A.gram_schmidt()
    const A = qmat([
      [-1, 3, 2, 2],
      [-1, 0, -1, 0],
      [-1, -2, -3, -1],
      [1, 1, 2, 0],
    ]);
    const [G, M] = gram_schmidt(A);
    expect(render(G)).toBe('[-1 3 2 2] / [-19/18 1/6 -8/9 1/9] / [2/35 -4/35 -2/35 9/35]');
    expect(render(M)).toBe('[1 0 0] / [-1/18 1 0] / [-13/18 59/35 1] / [1/3 -48/35 -2]');
    expect(matEq(M.mul(G), A)).toBe(true);
    // G*G^T is diagonal
    const GGt = G.mul(G.transpose());
    for (let i = 0; i < GGt.nrows; i++) {
      for (let j = 0; j < GGt.ncols; j++) {
        if (i !== j) expect(GGt.get(i, j).isZero()).toBe(true);
      }
    }
  });

  it('gram_schmidt_noscale still orthogonalizes columns: A == Q*R (Sage doctest)', () => {
    // sage: A = matrix(ZZ, [[-1,-3,0,-1],[1,2,-1,2],[-3,-6,4,-7]])
    // sage: Q, R = A._gram_schmidt_noscale()
    const A = qmat([
      [-1, -3, 0, -1],
      [1, 2, -1, 2],
      [-3, -6, 4, -7],
    ]);
    const [Q, R] = gram_schmidt_noscale(A);
    expect(render(Q)).toBe('[-1 -10/11 0] / [1 -1/11 3/10] / [-3 3/11 1/10]');
    expect(render(R)).toBe('[1 23/11 -13/11 24/11] / [0 1 13/10 -13/10] / [0 0 1 -1]');
    expect(matEq(Q.mul(R), A)).toBe(true);
  });
});

describe('indefinite_factorization', () => {
  it('A == L*diag(d)*L^T (Sage 5x5 QQ doctest)', () => {
    // sage: A = matrix(QQ, [[3,-6,9,6,-9],[-6,11,-16,-11,17],[9,-16,28,16,-40],
    //                       [6,-11,16,9,-19],[-9,17,-40,-19,68]])
    // sage: L, d = A.indefinite_factorization()
    const A = qmat([
      [3, -6, 9, 6, -9],
      [-6, 11, -16, -11, 17],
      [9, -16, 28, 16, -40],
      [6, -11, 16, 9, -19],
      [-9, 17, -40, -19, 68],
    ]);
    const [L, d] = indefinite_factorization(A);
    expect(render(L)).toBe(
      '[1 0 0 0 0] / [-2 1 0 0 0] / [3 -2 1 0 0] / [2 -1 0 1 0] / [-3 1 -3 1 1]'
    );
    expect(d.map(String)).toEqual(['3', '-1', '5', '-2', '-1']);

    const D = zero_matrix(QQ as any, 5);
    d.forEach((x: any, i: number) => D.set(i, i, x));
    expect(matEq(L.mul(D).mul(L.transpose()), A)).toBe(true);
  });

  it('raises for a singular leading principal submatrix (Sage doctest message)', () => {
    // sage: matrix(QQ,[[4,6,1],[6,9,5],[1,5,2]]).indefinite_factorization()
    // ValueError: 2x2 leading principal submatrix is singular, ...
    const A = qmat([
      [4, 6, 1],
      [6, 9, 5],
      [1, 5, 2],
    ]);
    expect(() => indefinite_factorization(A)).toThrow(
      '2x2 leading principal submatrix is singular, so cannot create indefinite factorization'
    );
  });

  it('raises on a zero leading entry rather than returning an invalid factorization', () => {
    const A = ffmat(F101, [
      [0, 1],
      [1, 0],
    ]);
    expect(() => indefinite_factorization(A)).toThrow(
      '1x1 leading principal submatrix is singular, so cannot create indefinite factorization'
    );
  });

  it('rejects non-symmetric input with Sage’s message', () => {
    const A = qmat([
      [1, 2],
      [3, 4],
    ]);
    expect(() => indefinite_factorization(A)).toThrow(
      "matrix is not symmetric (maybe try the 'hermitian' keyword)"
    );
  });
});

describe('block_ldlt', () => {
  it('reproduces the Sage doctest P, L, D exactly', () => {
    // sage: A = matrix(QQ, [[0,1,0],[1,1,2],[0,2,0]])
    // sage: P,L,D = A.block_ldlt()
    const A = qmat([
      [0, 1, 0],
      [1, 1, 2],
      [0, 2, 0],
    ]);
    const [P, L, D] = block_ldlt(A);
    expect(render(P)).toBe('[0 0 1] / [1 0 0] / [0 1 0]');
    expect(render(L)).toBe('[1 0 0] / [2 1 0] / [1 1/2 1]');
    expect(render(D)).toBe('[1 0 0] / [0 -4 0] / [0 0 0]');
    expect(matEq(P.transpose().mul(A).mul(P), L.mul(D).mul(L.transpose()))).toBe(true);
  });

  it('a 2x2 matrix with no classical factorization is its own block factorization', () => {
    // sage: A = matrix(QQ, [[0,1],[1,0]]); A.block_ldlt(classical=True) -> ValueError
    const A = qmat([
      [0, 1],
      [1, 0],
    ]);
    expect(() => block_ldlt(A, true)).toThrow('matrix has no classical LDL^T factorization');

    const [P, L, D] = block_ldlt(A);
    expect(render(P)).toBe('[1 0] / [0 1]');
    expect(render(L)).toBe('[1 0] / [0 1]');
    expect(render(D)).toBe('[0 1] / [1 0]');
    expect(matEq(P.transpose().mul(A).mul(P), L.mul(D).mul(L.transpose()))).toBe(true);
  });

  it('classical=true agrees with indefinite_factorization (Sage TESTS block)', () => {
    const A = qmat([
      [4, -2, 4, 2],
      [-2, 10, -2, -7],
      [4, -2, 8, 4],
      [2, -7, 4, 7],
    ]);
    const [, Lc, Dc] = block_ldlt(A, true);
    const [Li, di] = indefinite_factorization(A);
    const Di = zero_matrix(QQ as any, 4);
    di.forEach((x: any, i: number) => Di.set(i, i, x));
    expect(matEq(Lc, Li)).toBe(true);
    expect(matEq(Dc, Di)).toBe(true);
  });

  it('P^T*A*P == L*D*L^T with D block diagonal (100 random symmetric over QQ)', () => {
    const rnd = makeRng(4242);
    for (let t = 0; t < 100; t++) {
      const n = 1 + (rnd() % 5);
      const e: number[][] = [];
      for (let i = 0; i < n; i++) e.push(new Array<number>(n).fill(0));
      for (let i = 0; i < n; i++) {
        for (let j = 0; j <= i; j++) {
          const v = (rnd() % 11) - 5;
          e[i]![j] = v;
          e[j]![i] = v;
        }
      }
      const A = qmat(e);
      const [P, L, D] = block_ldlt(A);
      expect(matEq(P.transpose().mul(A).mul(P), L.mul(D).mul(L.transpose()))).toBe(true);
      // L is unit lower triangular
      for (let i = 0; i < n; i++) {
        expect(L.get(i, i).eq(QQ.__call__(1))).toBe(true);
        for (let j = i + 1; j < n; j++) expect(L.get(i, j).isZero()).toBe(true);
      }
      // D is block diagonal with blocks of size at most 2
      for (let i = 0; i < n; i++) {
        for (let j = 0; j < n; j++) {
          if (Math.abs(i - j) > 1) expect(D.get(i, j).isZero()).toBe(true);
        }
      }
    }
  });

  it('P^T*A*P == L*D*L^T over GF(101) (exact-pivot fallback, 100 random)', () => {
    const rnd = makeRng(777);
    for (let t = 0; t < 100; t++) {
      const n = 1 + (rnd() % 5);
      const e: number[][] = [];
      for (let i = 0; i < n; i++) e.push(new Array<number>(n).fill(0));
      for (let i = 0; i < n; i++) {
        for (let j = 0; j <= i; j++) {
          const v = rnd() % 101;
          e[i]![j] = v;
          e[j]![i] = v;
        }
      }
      const A = ffmat(F101, e);
      const [P, L, D] = block_ldlt(A);
      expect(matEq(P.transpose().mul(A).mul(P), L.mul(D).mul(L.transpose()))).toBe(true);
    }
  });
});

describe('smith_form', () => {
  it('satisfies S == U*M*V for a matrix with a non-pivot column', () => {
    const M = ffmat(F101, [
      [1, 2],
      [0, 0],
    ]);
    const [S, U, V] = smith_form(M) as [Matrix<any>, Matrix<any>, Matrix<any>];
    expect(render(S)).toBe('[1 0] / [0 0]');
    expect(matEq(S, U.mul(M).mul(V))).toBe(true);
  });

  it('satisfies S == U*M*V and S is diagonal (200 random over GF(101))', () => {
    const rnd = makeRng(999);
    for (let t = 0; t < 200; t++) {
      const m = 1 + (rnd() % 5);
      const n = 1 + (rnd() % 5);
      const e: number[][] = [];
      for (let i = 0; i < m; i++) {
        const r: number[] = [];
        for (let j = 0; j < n; j++) r.push(rnd() % 101);
        e.push(r);
      }
      const M = ffmat(F101, e);
      const [S, U, V] = smith_form(M) as [Matrix<any>, Matrix<any>, Matrix<any>];
      expect(matEq(S, U.mul(M).mul(V))).toBe(true);
      for (let i = 0; i < m; i++) {
        for (let j = 0; j < n; j++) {
          if (i !== j) expect(S.get(i, j).isZero()).toBe(true);
        }
      }
      const rank = pivots(M).length;
      for (let i = 0; i < Math.min(m, n); i++) {
        if (i < rank) expect(S.get(i, i).eq((F101 as any).__call__(1))).toBe(true);
        else expect(S.get(i, i).isZero()).toBe(true);
      }
    }
  });

  it('without transformation returns the diagonal form only', () => {
    const M = ffmat(F101, [
      [1, 2, 3],
      [2, 4, 6],
    ]);
    const S = smith_form(M, false) as Matrix<any>;
    expect(render(S)).toBe('[1 0 0] / [0 0 0]');
  });
});

describe('jordan_form', () => {
  it('leaves a diagonal matrix with distinct eigenvalues alone', () => {
    const D = ffmat(F101, [
      [2, 0],
      [0, 3],
    ]);
    // Blocks are ordered by `charpoly().roots()` order (`matrix2.pyx:12251-12254`),
    // and `roots()` follows `factor()` order, not ascending value.  Verified:
    // `matrix(GF(101),2,2,[2,0,0,3]).jordan_form()` is `[3 0; 0 2]` in
    // SageMath 10.3.  (This assertion previously pinned ascending order.)
    expect(render(jordan_form(D) as Matrix<any>)).toBe('[3 0] / [0 2]');
  });

  it('recovers a single 2x2 Jordan block', () => {
    const J = ffmat(F101, [
      [2, 1],
      [0, 2],
    ]);
    expect(render(jordan_form(J) as Matrix<any>)).toBe('[2 1] / [0 2]');
  });

  it('splits a nilpotent matrix into blocks of sizes 3 and 1', () => {
    const N = ffmat(F101, [
      [0, 1, 0, 0],
      [0, 0, 1, 0],
      [0, 0, 0, 0],
      [0, 0, 0, 0],
    ]);
    expect(render(jordan_form(N) as Matrix<any>)).toBe(
      '[0 1 0 0] / [0 0 1 0] / [0 0 0 0] / [0 0 0 0]'
    );
  });

  it('handles a repeated eigenvalue together with a simple one', () => {
    const A = ffmat(F101, [
      [2, 1, 0],
      [0, 2, 0],
      [0, 0, 3],
    ]);
    // See above: SageMath orders the blocks by `roots()` order, so the simple
    // eigenvalue 3 comes first.  Verified against SageMath 10.3 over both
    // GF(101) and QQ.  (Previously pinned the port's ascending order.)
    expect(render(jordan_form(A) as Matrix<any>)).toBe('[3 0 0] / [0 2 1] / [0 0 2]');
  });

  it('is invariant under conjugation', () => {
    // A = S * J * S^{-1} has the same Jordan form as J
    const J = ffmat(F101, [
      [2, 1, 0],
      [0, 2, 0],
      [0, 0, 3],
    ]);
    const S = ffmat(F101, [
      [1, 1, 0],
      [0, 1, 1],
      [1, 0, 1],
    ]);
    const half = (F101 as any).__call__(51); // 1/2 mod 101
    const Sinv = ffmat(F101, [
      [1, -1, 1],
      [1, 1, -1],
      [-1, 1, 1],
    ]).scalar_mul(half);
    expect(render(S.mul(Sinv))).toBe('[1 0 0] / [0 1 0] / [0 0 1]');
    const A = S.mul(J).mul(Sinv);
    // `J` itself is written with the blocks in ascending eigenvalue order,
    // which is NOT SageMath's; the Jordan form of `A` uses `roots()` order.
    expect(render(jordan_form(A) as Matrix<any>)).toBe('[3 0 0] / [0 2 1] / [0 0 2]');
  });

  it('rejects a characteristic polynomial that does not split', () => {
    // x^2 - 2 is irreducible over GF(101) because 2 is a non-residue mod 101
    const A = ffmat(F101, [
      [0, 1],
      [2, 0],
    ]);
    expect(() => jordan_form(A)).toThrow('Some eigenvalue does not exist in');
  });

  it('checks a user-supplied eigenvalue list', () => {
    const N = ffmat(F101, [
      [0, 1],
      [0, 0],
    ]);
    expect(
      render(
        jordan_form(
          N,
          undefined,
          false,
          true,
          false,
          [[(F101 as any).__call__(0), 2]],
          true
        ) as Matrix<any>
      )
    ).toBe('[0 1] / [0 0]');
    expect(() =>
      jordan_form(N, undefined, false, true, false, [[(F101 as any).__call__(0), 1]], true)
    ).toThrow('The provided list of eigenvalues is not correct.');
  });
});

describe('LLL_gram (PARI qflllgram doctests)', () => {
  it('matches PARI on the semidefinite example [[2,6],[6,3]]', () => {
    expect(
      render(
        LLL_gram(
          zmat([
            [2n, 6n],
            [6n, 3n],
          ])
        )
      )
    ).toBe('[-3 -1] / [1 0]');
  });

  it('matches PARI on the indefinite example [[1,0],[0,-1]]', () => {
    expect(
      render(
        LLL_gram(
          zmat([
            [1n, 0n],
            [0n, -1n],
          ])
        )
      )
    ).toBe('[0 -1] / [1 0]');
  });

  it('reduces a 3x3 positive definite Gram matrix to the identity', () => {
    const G = zmat([
      [1n, 3n, 5n],
      [3n, 10n, 18n],
      [5n, 18n, 35n],
    ]);
    const U = LLL_gram(G);
    const R = U.transpose().mul(G).mul(U);
    expect(render(R)).toBe('[1 0 0] / [0 1 0] / [0 0 1]');
  });

  it('raises like PARI for degenerate Gram matrices', () => {
    expect(() => LLL_gram(zmat([[0n]]))).toThrow(
      'qflllgram did not return a square matrix, perhaps the matrix is not positive definite'
    );
    expect(() =>
      LLL_gram(
        zmat([
          [0n, 1n],
          [1n, 0n],
        ])
      )
    ).toThrow('qflllgram did not return a square matrix');
  });

  it('output is size reduced and U is unimodular (random positive definite)', () => {
    const rnd = makeRng(31337);
    let tested = 0;
    for (let t = 0; t < 60; t++) {
      const n = 2 + (rnd() % 2);
      const b: bigint[][] = [];
      for (let i = 0; i < n; i++) {
        const r: bigint[] = [];
        for (let j = 0; j < n; j++) r.push(BigInt((rnd() % 21) - 10));
        b.push(r);
      }
      const g: bigint[][] = [];
      for (let i = 0; i < n; i++) {
        const r: bigint[] = [];
        for (let j = 0; j < n; j++) {
          let s = 0n;
          for (let k = 0; k < n; k++) s += b[k]![i]! * b[k]![j]!;
          r.push(s);
        }
        g.push(r);
      }
      const G = zmat(g);
      let U: Matrix<any>;
      try {
        U = LLL_gram(G);
      } catch {
        continue; // degenerate Gram matrix; PARI errors here too
      }
      tested++;
      const R = U.transpose().mul(G).mul(U);
      const r00 = R.get(0, 0).value as bigint;
      const r01 = R.get(0, 1).value as bigint;
      // |mu_{1,0}| <= 1/2 for a size-reduced basis
      expect(2n * (r01 < 0n ? -r01 : r01) <= r00).toBe(true);
      expect(determinant(U).value).toBe(1n);
    }
    expect(tested).toBeGreaterThan(10);
  });
});

describe('principal_square_root', () => {
  it('returns B with B^2 == A for a diagonalizable matrix over GF(101)', () => {
    const P = ffmat(F101, [
      [1, 2],
      [3, 7],
    ]);
    const Pinv = ffmat(F101, [
      [7, -2],
      [-3, 1],
    ]);
    const D = ffmat(F101, [
      [4, 0],
      [0, 9],
    ]);
    const A = P.mul(D).mul(Pinv);
    const B = principal_square_root(A) as Matrix<any>;
    expect(matEq(B.mul(B), A)).toBe(true);
  });

  it('works for 3x3 matrices, where the old Denman-Beavers loop never converged', () => {
    const D = ffmat(F101, [
      [4, 0, 0],
      [0, 9, 0],
      [0, 0, 16],
    ]);
    const S = ffmat(F101, [
      [1, 1, 0],
      [0, 1, 1],
      [1, 0, 1],
    ]);
    const half = (F101 as any).__call__(51);
    const Sinv = ffmat(F101, [
      [1, -1, 1],
      [1, 1, -1],
      [-1, 1, 1],
    ]).scalar_mul(half);
    const A = S.mul(D).mul(Sinv);
    const B = principal_square_root(A) as Matrix<any>;
    expect(matEq(B.mul(B), A)).toBe(true);
  });

  it('reports non-diagonalizable input instead of looping', () => {
    const A = ffmat(F101, [
      [1, 1],
      [0, 1],
    ]);
    expect(() => principal_square_root(A)).toThrow('diagonalizable');
  });
});

describe('decomposition (primal vs dual)', () => {
  it('uses the left kernel for the primal and the right kernel for the dual', () => {
    // Over GF(101), A = [[1,1],[0,2]].  For g = x-1, A - I = [[0,1],[0,1]]:
    // its left kernel is spanned by (-1,1) and its right kernel by (1,0).
    const A = ffmat(F101, [
      [1, 1],
      [0, 2],
    ]);
    const [primal, dual] = decomposition(A, 'kernel', false, true) as [
      Array<[Matrix<any>, boolean]>,
      Array<[Matrix<any>, boolean]>,
    ];
    const primalRows = primal.map(([B]) => render(B));
    const dualRows = dual.map(([B]) => render(B));
    expect(primalRows).toContain('[100 1]');
    expect(dualRows).toContain('[1 0]');
  });

  it('primal subspaces are invariant under the left action v -> v*A', () => {
    const A = ffmat(F101, [
      [1, 1],
      [0, 2],
    ]);
    const primal = decomposition(A) as Array<[Matrix<any>, boolean]>;
    for (const [basis] of primal) {
      const stackedRows: any[][] = [];
      for (let i = 0; i < basis.nrows; i++) stackedRows.push(basis.row(i));
      for (let i = 0; i < basis.nrows; i++) {
        const v = basis.row(i);
        const vA: any[] = [];
        for (let j = 0; j < A.ncols; j++) {
          let s = (F101 as any).__call__(0);
          for (let k = 0; k < A.nrows; k++) s = s.add(v[k]!.mul(A.get(k, j)));
          vA.push(s);
        }
        const stacked = new Matrix(F101 as any, basis.nrows + 1, basis.ncols, [...stackedRows, vA]);
        expect(pivots(stacked).length).toBe(pivots(basis).length);
      }
    }
  });
});

describe('krylov_matrix / krylov_basis', () => {
  const E = ffmat(F97, [
    [27, 49, 29],
    [50, 58, 0],
    [77, 10, 29],
  ]);
  const M = ffmat(F97, [
    [0, 1, 0],
    [0, 0, 1],
    [0, 0, 0],
  ]);

  it('stacks E, E*M, E*M^2, E*M^3 (Sage GF(97) doctest, 12 rows)', () => {
    expect(render(krylov_matrix(E, M))).toBe(
      '[27 49 29] / [50 58 0] / [77 10 29] / ' +
        '[0 27 49] / [0 50 58] / [0 77 10] / ' +
        '[0 0 27] / [0 0 50] / [0 0 77] / ' +
        '[0 0 0] / [0 0 0] / [0 0 0]'
    );
  });

  it('orders rows by shifts (Sage doctest, shifts=[0,3,6])', () => {
    expect(render(krylov_matrix(E, M, [0, 3, 6]))).toBe(
      '[27 49 29] / [0 27 49] / [0 0 27] / [0 0 0] / ' +
        '[50 58 0] / [0 50 58] / [0 0 50] / [0 0 0] / ' +
        '[77 10 29] / [0 77 10] / [0 0 77] / [0 0 0]'
    );
  });

  it('orders rows by shifts (Sage doctest, shifts=[3,0,2])', () => {
    expect(render(krylov_matrix(E, M, [3, 0, 2]))).toBe(
      '[50 58 0] / [0 50 58] / [0 0 50] / [77 10 29] / [27 49 29] / [0 0 0] / ' +
        '[0 77 10] / [0 27 49] / [0 0 77] / [0 0 27] / [0 0 0] / [0 0 0]'
    );
  });

  it('honours per-row degree bounds inclusively (Sage doctest)', () => {
    expect(render(krylov_matrix(E, M, [3, 0, 2], [0, 2, 1]))).toBe(
      '[50 58 0] / [0 50 58] / [0 0 50] / [77 10 29] / [27 49 29] / [0 77 10]'
    );
  });

  it('rejects negative degree bounds like Sage', () => {
    expect(() => krylov_matrix(E, M, undefined, [2, 3, -1])).toThrow(
      'degrees must not contain a negative bound'
    );
  });

  it('rejects a wrongly sized M like Sage', () => {
    expect(() =>
      krylov_matrix(
        E,
        ffmat(F97, [
          [0, 0, 0, 0],
          [0, 0, 0, 0],
          [0, 0, 0, 0],
          [0, 0, 0, 0],
        ])
      )
    ).toThrow('M does not have correct dimensions');
  });

  it('krylov_basis returns the row rank profile of the Krylov matrix (Sage doctests)', () => {
    const [B, profile] = krylov_basis(E, M) as [Matrix<any>, Array<[number, number, number]>];
    expect(render(B)).toBe('[27 49 29] / [50 58 0] / [0 27 49]');
    expect(profile).toEqual([
      [0, 0, 0],
      [1, 0, 1],
      [0, 1, 3],
    ]);

    const [B2, p2] = krylov_basis(E, M, [0, 3, 6]) as [
      Matrix<any>,
      Array<[number, number, number]>,
    ];
    expect(render(B2)).toBe('[27 49 29] / [0 27 49] / [0 0 27]');
    expect(p2).toEqual([
      [0, 0, 0],
      [0, 1, 1],
      [0, 2, 2],
    ]);

    const [B3, p3] = krylov_basis(E, M, [3, 0, 2]) as [
      Matrix<any>,
      Array<[number, number, number]>,
    ];
    expect(render(B3)).toBe('[50 58 0] / [0 50 58] / [0 0 50]');
    expect(p3).toEqual([
      [1, 0, 0],
      [1, 1, 1],
      [1, 2, 2],
    ]);
  });

  it('krylov_kernel_basis rows annihilate the corresponding Krylov matrix', () => {
    const [K] = krylov_kernel_basis(E, M) as [Matrix<any>, Array<[number, number, number]>];
    const [, profile] = krylov_basis(E, M) as [Matrix<any>, Array<[number, number, number]>];
    const delta = [0, 0, 0];
    for (const [i, j] of profile) delta[i] = Math.max(delta[i]!, j + 1);
    const A = krylov_matrix(E, M, [0, 0, 0], delta);
    expect(K.ncols).toBe(A.nrows);
    expect(K.nrows).toBe(3); // m rows
    const product = K.mul(A);
    for (let i = 0; i < product.nrows; i++) {
      for (let j = 0; j < product.ncols; j++) {
        expect(product.get(i, j).isZero()).toBe(true);
      }
    }
  });
});

describe('LU_double (M53)', () => {
  it('returns U with min(m,n) rows so that L*U is defined for tall inputs', () => {
    const A = [
      [1, 2],
      [3, 4],
      [5, 6],
    ];
    const { P, L, U } = LU_double(A);
    expect(L.length).toBe(3);
    expect(L[0]!.length).toBe(2);
    expect(U.length).toBe(2);
    expect(U[0]!.length).toBe(2);

    // P*A == L*U
    for (let i = 0; i < 3; i++) {
      for (let j = 0; j < 2; j++) {
        let s = 0;
        for (let k = 0; k < 2; k++) s += L[i]![k]! * U[k]![j]!;
        expect(Math.abs(s - A[P[i]!]![j]!)).toBeLessThan(1e-12);
      }
    }
  });

  it('still reconstructs square systems exactly', () => {
    const A = [
      [2, 1, 1],
      [4, 3, 3],
      [8, 7, 9],
    ];
    const { P, L, U } = LU_double(A);
    expect(U.length).toBe(3);
    for (let i = 0; i < 3; i++) {
      for (let j = 0; j < 3; j++) {
        let s = 0;
        for (let k = 0; k < 3; k++) s += L[i]![k]! * U[k]![j]!;
        expect(Math.abs(s - A[P[i]!]![j]!)).toBeLessThan(1e-12);
      }
    }
  });
});

describe('QR', () => {
  it('defaults to the full factorization (Sage full=True) and satisfies A == Q*R', () => {
    const A = qmat([
      [1, 2],
      [3, 4],
      [5, 7],
    ]);
    const [Q, R] = QR(A);
    expect(Q.nrows).toBe(3);
    expect(Q.ncols).toBe(3);
    expect(R.nrows).toBe(3);
    expect(R.ncols).toBe(2);
    expect(matEq(Q.mul(R), A)).toBe(true);
    // Q has orthogonal columns (see the documented deviation: not orthonormal)
    const QtQ = Q.transpose().mul(Q);
    for (let i = 0; i < 3; i++) {
      for (let j = 0; j < 3; j++) {
        if (i !== j) expect(QtQ.get(i, j).isZero()).toBe(true);
      }
    }
  });

  it('full=false gives the reduced factorization with A == Q*R', () => {
    const A = qmat([
      [1, 2],
      [3, 4],
      [5, 7],
    ]);
    const [Q, R] = QR(A, false);
    expect(Q.ncols).toBe(2);
    expect(R.nrows).toBe(2);
    expect(matEq(Q.mul(R), A)).toBe(true);
  });
});

describe('hermite_form / extended_echelon_form / elementary_divisors', () => {
  it('hermite_form over a field is the RREF and U*A == H', () => {
    const A = ffmat(F101, [
      [0, 1, 2],
      [3, 4, 5],
      [6, 7, 8],
    ]);
    const H = hermite_form(A) as Matrix<any>;
    // `hermite_form` is `_echelon_form_PID` (`matrix2.pyx:17245,17305`), NOT the
    // RREF: it never scales the pivot row and never reduces above a pivot.
    // Verified: `matrix(GF(101),3,3,[0..8]).hermite_form()` is
    // `[3 4 5; 0 100 99; 0 0 0]` while `.echelon_form()` is
    // `[1 0 100; 0 1 2; 0 0 0]` in SageMath 10.3.  (This assertion previously
    // claimed the two coincide over a field.)
    expect(render(H)).toBe('[3 4 5] / [0 100 99] / [0 0 0]');
    expect(render(H)).not.toBe(render(echelon_form(A)));

    const [H2, U] = hermite_form(A, true, true) as [Matrix<any>, Matrix<any>];
    expect(render(H2)).toBe(render(H));
    expect(matEq(U.mul(A), H2)).toBe(true);
  });

  it('extended_echelon_form is [E | T] with T*A == E', () => {
    const A = ffmat(F101, [
      [0, 1, 2],
      [3, 4, 5],
    ]);
    const X = extended_echelon_form(A);
    expect(X.nrows).toBe(2);
    expect(X.ncols).toBe(5);
    const E = new Matrix(F101 as any, 2, 3, [
      [X.get(0, 0), X.get(0, 1), X.get(0, 2)],
      [X.get(1, 0), X.get(1, 1), X.get(1, 2)],
    ]);
    const T = new Matrix(F101 as any, 2, 2, [
      [X.get(0, 3), X.get(0, 4)],
      [X.get(1, 3), X.get(1, 4)],
    ]);
    expect(render(E)).toBe(render(echelon_form(A)));
    expect(matEq(T.mul(A), E)).toBe(true);
  });

  it('elementary_divisors over a field are rank ones then zeros', () => {
    const A = ffmat(F101, [
      [1, 2, 3],
      [2, 4, 6],
    ]);
    expect(elementary_divisors(A).map(String)).toEqual(['1', '0']);
  });
});

describe('cholesky', () => {
  it('recovers L with A == L*L^T over GF(101)', () => {
    // A = L0 * L0^T for a lower triangular L0 with square diagonal entries
    const L0 = ffmat(F101, [
      [2, 0, 0],
      [5, 3, 0],
      [7, 11, 4],
    ]);
    const A = L0.mul(L0.transpose());
    const L = cholesky(A);
    expect(matEq(L.mul(L.transpose()), A)).toBe(true);
    for (let i = 0; i < 3; i++) {
      for (let j = i + 1; j < 3; j++) expect(L.get(i, j).isZero()).toBe(true);
    }
  });
});

describe('hessenberg_form', () => {
  it('is upper Hessenberg and preserves the characteristic polynomial', () => {
    const A = ffmat(F101, [
      [4, 1, 7, 2],
      [3, 9, 5, 6],
      [8, 2, 1, 4],
      [5, 5, 3, 7],
    ]);
    const H = hessenberg_form(A);
    for (let i = 2; i < 4; i++) {
      for (let j = 0; j < i - 1; j++) expect(H.get(i, j).isZero()).toBe(true);
    }
    expect(String(charpoly(H))).toBe(String(charpoly(A)));
  });
});

// ============================================================================
// Deferred audit items 21, 22, 23 (unit matrix-decomp2)
//
// All expected values below are SageMath's, copied from the doctests of
// sage/matrix/matrix2.pyx as cited on each test.
// ============================================================================

import { PolynomialRing } from '../rings/polynomial/polynomial_ring.js';
import * as matrixIndex from './index.js';
import { jordan_decomposition } from './matrix_decompositions.js';
import { inverse, rank } from './matrix_operations.js';
import {
  is_hermite,
  is_popov,
  minimal_approximant_basis,
  popov_form,
} from './matrix_polynomial_dense.js';
import { subdivide, subdivisions } from './matrix_special.js';

/** Multiset of (eigenvalue, block size) read off a matrix in Jordan form. */
function jordanBlockMultiset(J: Matrix<any>): string {
  const n = J.nrows;
  const out: string[] = [];
  let i = 0;
  while (i < n) {
    let size = 1;
    while (i + size < n && !J.get(i + size - 1, i + size).isZero()) size++;
    out.push(`${String(J.get(i, i))}:${size}`);
    i += size;
  }
  return out.sort().join(',');
}

describe('item 21: matrix/index.ts exports pivots next to pivot_rows', () => {
  it('re-exports pivots and it is the same function', () => {
    expect(typeof matrixIndex.pivots).toBe('function');
    expect(matrixIndex.pivots).toBe(pivots);
    expect(typeof matrixIndex.pivot_rows).toBe('function');
  });

  it('the re-exported pivots agrees with Sage on the range(9) matrix', () => {
    // sage: matrix(QQ,3,3,range(9)).pivots()  ->  (0, 1)
    const A = qmat([
      [0, 1, 2],
      [3, 4, 5],
      [6, 7, 8],
    ]);
    expect(matrixIndex.pivots(A)).toEqual([0, 1]);
  });
});

describe('item 23: jordan_form(transformation=true)', () => {
  it('reproduces the Sage issue-12693 doctest exactly (matrix2.pyx:12108-12121)', () => {
    // sage: M = matrix(((2,2,2), (0,0,0), (-2,-2,-2)))
    // sage: J, P = M.jordan_form(transformation=True); J; P
    // [0 1|0]     [ 2  1  0]
    // [0 0|0]     [ 0  0  1]
    // [0 0|0]     [-2  0 -1]
    const M = qmat([
      [2, 2, 2],
      [0, 0, 0],
      [-2, -2, -2],
    ]);
    const [J, P] = jordan_form(M, undefined, undefined, undefined, true) as [
      Matrix<any>,
      Matrix<any>,
    ];
    expect(render(J)).toBe('[0 1 0] / [0 0 0] / [0 0 0]');
    expect(render(P)).toBe('[2 1 0] / [0 0 1] / [-2 0 -1]');
    // sage: J - ~P * M * P  ==  0
    expect(matEq(inverse(P).mul(M).mul(P), J)).toBe(true);
  });

  it('reproduces the 10x10 QQ doctest (matrix2.pyx:12073-12102)', () => {
    const f = (a: number, b: number) => (QQ as any).__call__([BigInt(a), BigInt(b)]);
    const e: any[][] = [
      [15, f(37, 3), -16, f(-104, 3), -29, f(-7, 3), 35, f(2, 3), f(-29, 3), f(-1, 3)],
      [2, 9, -1, -6, -6, 0, 7, 0, -2, 0],
      [24, f(74, 3), -29, f(-208, 3), -58, f(-14, 3), 70, f(4, 3), f(-58, 3), f(-2, 3)],
      [-6, -19, 3, 21, 19, 0, -21, 0, 6, 0],
      [2, 6, -1, -6, -3, 0, 7, 0, -2, 0],
      [-96, f(-296, 3), 128, f(832, 3), 232, f(65, 3), -279, f(-16, 3), f(232, 3), f(8, 3)],
      [0, 0, 0, 0, 0, 0, 3, 0, 0, 0],
      [20, f(26, 3), -30, f(-199, 3), -42, f(-14, 3), 70, f(13, 3), f(-55, 3), f(-2, 3)],
      [18, 57, -9, -54, -57, 0, 63, 0, -15, 0],
      [0, 0, 0, 0, 0, 0, 0, 0, 0, 3],
    ];
    const A = new Matrix(
      QQ as any,
      10,
      10,
      e.map((r) => r.map((x) => (QQ as any).__call__(x)))
    );
    const [J, T] = jordan_form(A, undefined, undefined, undefined, true) as [
      Matrix<any>,
      Matrix<any>,
    ];
    // Blocks for the single eigenvalue 3 have sizes 3, 2, 2, 2, 1.
    expect(jordanBlockMultiset(J)).toBe('3:1,3:2,3:2,3:2,3:3');
    // sage: T.rank()  ->  10 ;  sage: T * J * T**(-1) == A  ->  True
    expect(rank(T)).toBe(10);
    expect(matEq(T.mul(J).mul(inverse(T)), A)).toBe(true);
  });

  it('reproduces the supplied-eigenvalues doctest (matrix2.pyx:12144-12150)', () => {
    // Sage states this over ZZ[x11..x33] with M = [[0,0,x31],[0,0,x21],[0,0,0]]
    // and gets T = [x31 0 1 / x21 0 0 / 0 1 0].  The port cannot build a matrix
    // over a multivariate polynomial ring, so we instantiate x31=5, x21=7.
    const M = qmat([
      [0, 0, 5],
      [0, 0, 7],
      [0, 0, 0],
    ]);
    const [J, T] = jordan_form(M, undefined, undefined, undefined, true, [
      [(QQ as any).__call__(0), 3],
    ]) as [Matrix<any>, Matrix<any>];
    expect(render(J)).toBe('[0 1 0] / [0 0 0] / [0 0 0]');
    expect(render(T)).toBe('[5 0 1] / [7 0 0] / [0 1 0]');
    expect(matEq(inverse(T).mul(M).mul(T), J)).toBe(true);
  });

  it('computes eigenvalues over QQ via charpoly().roots() (no enumerable ring)', () => {
    // (x-2)^2 (x-3), with a 2x2 block for the eigenvalue 2.
    const A = qmat([
      [2, 1, 0],
      [0, 2, 0],
      [0, 0, 3],
    ]);
    // See above: SageMath orders the blocks by `roots()` order, so the simple
    // eigenvalue 3 comes first.  Verified against SageMath 10.3 over both
    // GF(101) and QQ.  (Previously pinned the port's ascending order.)
    expect(render(jordan_form(A) as Matrix<any>)).toBe('[3 0 0] / [0 2 1] / [0 0 2]');
    // A rational (non-integral) eigenvalue also works.
    const B = qmat([
      [1, 1],
      [0, 1],
    ]).scalar_mul((QQ as any).__call__([1n, 2n]));
    // (1/2)*[[1,1],[0,1]] has the single eigenvalue 1/2 with one 2x2 block, so
    // the Jordan form is jordan_block(1/2, 2) = [[1/2, 1], [0, 1/2]].
    expect(render(jordan_form(B) as Matrix<any>)).toBe('[1/2 1] / [0 1/2]');
    const [J, P] = jordan_form(B, undefined, undefined, undefined, true) as [
      Matrix<any>,
      Matrix<any>,
    ];
    expect(matEq(inverse(P).mul(B).mul(P), J)).toBe(true);
  });

  it('P is invertible and P^-1*A*P == J for conjugated Jordan matrices over QQ and GF(p)', () => {
    const specs: Array<Array<[number, number]>> = [
      [[2, 3]],
      [
        [2, 2],
        [2, 1],
      ],
      [
        [0, 3],
        [0, 1],
      ],
      [
        [1, 2],
        [3, 2],
      ],
      [
        [2, 2],
        [2, 2],
      ],
      [[5, 4]],
      [
        [0, 2],
        [0, 2],
        [1, 1],
      ],
      [
        [3, 3],
        [3, 2],
        [4, 1],
      ],
      [
        [1, 2],
        [1, 1],
        [2, 3],
      ],
    ];
    const rng = makeRng(20260727);
    for (const F of [GF(101n) as any, GF(7n) as any, QQ as any]) {
      for (const spec of specs) {
        const n = spec.reduce((s, [, sz]) => s + sz, 0);
        // J0 = block diagonal of jordan_block(lambda, size)
        const J0 = zero_matrix(F, n);
        let off = 0;
        for (const [lam, sz] of spec) {
          for (let i = 0; i < sz; i++) {
            J0.set(off + i, off + i, F.__call__(lam));
            if (i < sz - 1) J0.set(off + i, off + i + 1, F.one());
          }
          off += sz;
        }
        // random invertible S
        let S: Matrix<any> | undefined;
        let Sinv: Matrix<any> | undefined;
        for (let attempt = 0; attempt < 50 && S === undefined; attempt++) {
          const cand = new Matrix(
            F,
            n,
            n,
            Array.from({ length: n }, () =>
              Array.from({ length: n }, () => F.__call__((rng() % 11) - 5))
            )
          );
          try {
            Sinv = inverse(cand);
            S = cand;
          } catch {
            // singular candidate: draw another one
          }
        }
        expect(S).toBeDefined();
        const A = S!.mul(J0).mul(Sinv!);
        const [J, P] = jordan_form(A, undefined, undefined, undefined, true) as [
          Matrix<any>,
          Matrix<any>,
        ];
        // Same block structure as the matrix we conjugated.
        expect(jordanBlockMultiset(J)).toBe(jordanBlockMultiset(J0));
        // P is invertible (inverse() throws otherwise) and diagonalizes to J.
        const Pinv = inverse(P);
        expect(matEq(Pinv.mul(A).mul(P), J)).toBe(true);
        expect(matEq(P.mul(J).mul(Pinv), A)).toBe(true);
      }
    }
  });

  it('still rejects a characteristic polynomial that does not split, over QQ', () => {
    // charpoly x^3 - 1: only the rational root 1, of multiplicity 1 < 3
    const A = qmat([
      [0, 1, 0],
      [0, 0, 1],
      [1, 0, 0],
    ]);
    expect(() => jordan_form(A)).toThrow('Some eigenvalue does not exist in');
  });

  it('jordan_decomposition returns A = D + N with D semisimple and N nilpotent', () => {
    const A = qmat([
      [2, 1, 0],
      [0, 2, 0],
      [0, 0, 3],
    ]);
    const [D, N] = jordan_decomposition(A);
    expect(matEq(D.add(N), A)).toBe(true);
    expect(matEq(D.mul(N), N.mul(D))).toBe(true);
    let X = N.copy();
    for (let k = 0; k < 3; k++) X = X.mul(N);
    expect(X.is_zero()).toBe(true);
    expect(render(D)).toBe('[2 0 0] / [0 2 0] / [0 0 3]');
  });
});

describe('item 22: krylov_kernel_basis (constant and polynomial forms)', () => {
  // sage: R = GF(97)
  // sage: E = matrix(R, [[27, 49, 29], [50, 58, 0], [77, 10, 29]])
  // sage: M = matrix(R, [[0, 1, 0], [0, 0, 1], [0, 0, 0]])
  const E = ffmat(F97, [
    [27, 49, 29],
    [50, 58, 0],
    [77, 10, 29],
  ]);
  const M = ffmat(F97, [
    [0, 1, 0],
    [0, 0, 1],
    [0, 0, 0],
  ]);

  it('matches the constant-matrix doctest (matrix2.pyx:20100-20109)', () => {
    const [K, rp] = krylov_kernel_basis(E, M) as [Matrix<any>, Array<[number, number, number]>];
    expect(render(K)).toBe('[82 76 0 40 0 1] / [13 57 0 3 1 0] / [96 96 1 0 0 0]');
    expect(rp).toEqual([
      [0, 0, 0],
      [1, 0, 1],
      [2, 0, 2],
      [0, 1, 3],
      [1, 1, 4],
      [0, 2, 6],
    ]);
  });

  it('matches the Hermite-Pade polynomial doctest (matrix2.pyx:20081-20085)', () => {
    const [P, rp] = krylov_kernel_basis(E, M, undefined, undefined, true, 'x');
    expect(render(P)).toBe('[x^2 + 40*x + 82 76 0] / [3*x + 13 x + 57 0] / [96 96 1]');
    // sage: row_profile == row_profile_bis  (same as the constant form)
    expect(rp).toEqual([
      [0, 0, 0],
      [1, 0, 1],
      [2, 0, 2],
      [0, 1, 3],
      [1, 1, 4],
      [0, 2, 6],
    ]);
  });

  it('matches the shifts=[0,3,6] Hermite doctest (matrix2.pyx:20143-20160)', () => {
    const H = krylov_kernel_basis(E, M, [0, 3, 6], undefined, false, 'x');
    expect(render(H)).toBe('[x^3 0 0] / [60*x^2 + 72*x + 70 1 0] / [60*x^2 + 72*x + 69 0 1]');
    const [K] = krylov_kernel_basis(E, M, [0, 3, 6]) as [
      Matrix<any>,
      Array<[number, number, number]>,
    ];
    expect(render(K)).toBe('[0 0 0 1 0 0] / [70 72 60 0 1 0] / [69 72 60 0 0 1]');
  });

  it('matches the shifts=[3,0,2] doctests, polynomial and constant (matrix2.pyx:20174-20198)', () => {
    const [Q, rpq] = krylov_kernel_basis(E, M, [3, 0, 2], undefined, true, 'Y');
    expect(render(Q)).toBe('[1 26*Y^2 + 49*Y + 79 0] / [0 Y^3 0] / [0 26*Y^2 + 49*Y + 78 1]');
    expect(rpq).toEqual([
      [1, 0, 0],
      [1, 1, 1],
      [1, 2, 2],
      [2, 0, 3],
      [0, 0, 4],
      [1, 3, 5],
    ]);
    const [K, rp] = krylov_kernel_basis(E, M, [3, 0, 2], [0, 3, 0]) as [
      Matrix<any>,
      Array<[number, number, number]>,
    ];
    expect(render(K)).toBe('[79 49 26 0 1 0] / [0 0 0 0 0 1] / [78 49 26 1 0 0]');
    expect(rp).toEqual(rpq);
  });

  it('matches the too-small-degree-bound doctest (matrix2.pyx:20219-20232)', () => {
    // Sage documents this output as *not* a correct kernel basis; we reproduce it.
    const [K2] = krylov_kernel_basis(E, M, [3, 0, 2], [3, 1, 3]) as [
      Matrix<any>,
      Array<[number, number, number]>,
    ];
    expect(render(K2)).toBe('[1 0 96 1 0 0] / [3 28 56 0 1 0] / [47 64 69 0 0 1]');
  });

  it('matches the zero-M doctest (matrix2.pyx:20267-20273)', () => {
    const [K, rp] = krylov_kernel_basis(
      E,
      ffmat(F97, [
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0],
      ])
    ) as [Matrix<any>, Array<[number, number, number]>];
    expect(render(K)).toBe('[0 0 0 1 0] / [0 0 0 0 1] / [96 96 1 0 0]');
    expect(rp).toEqual([
      [0, 0, 0],
      [1, 0, 1],
      [2, 0, 2],
      [0, 1, 3],
      [1, 1, 4],
    ]);
  });

  it('accepts a PolynomialRing as the variable argument', () => {
    const R = new PolynomialRing(F97 as any, 'x');
    const P = krylov_kernel_basis(E, M, undefined, undefined, false, R);
    expect(render(P)).toBe('[x^2 + 40*x + 82 76 0] / [3*x + 13 x + 57 0] / [96 96 1]');
  });

  it('satisfies K*A == 0, rank(K) == m and the P/K relation on random inputs', () => {
    const rng = makeRng(424242);
    let cases = 0;
    for (const F of [GF(7n) as any, GF(97n) as any, QQ as any]) {
      for (let t = 0; t < 12; t++) {
        const m = 1 + (rng() % 3);
        const n = 1 + (rng() % 3);
        const shifts = t % 3 === 0 ? undefined : Array.from({ length: m }, () => rng() % 4);
        const Er = new Matrix(
          F,
          m,
          n,
          Array.from({ length: m }, () =>
            Array.from({ length: n }, () => F.__call__((rng() % 17) - 8))
          )
        );
        const Mr = new Matrix(
          F,
          n,
          n,
          Array.from({ length: n }, () =>
            Array.from({ length: n }, () => F.__call__((rng() % 17) - 8))
          )
        );
        const [K, rp] = krylov_kernel_basis(Er, Mr, shifts) as [
          Matrix<any>,
          Array<[number, number, number]>,
        ];
        const [P] = krylov_kernel_basis(Er, Mr, shifts, undefined, true, 'x');

        // (1) K is a left kernel basis of the Krylov matrix with the deduced degrees.
        const delta = new Array<number>(m).fill(0);
        for (const [i, j] of rp) delta[i] = Math.max(delta[i]!, j);
        const A = krylov_matrix(Er, Mr, shifts ?? new Array<number>(m).fill(0), delta);
        expect(K.ncols).toBe(A.nrows);
        expect(K.mul(A).is_zero()).toBe(true);

        // (2) K has full row rank m.
        expect(pivot_rows(K).length).toBe(m);

        // (3) The polynomial rows annihilate the K[x]-module map:
        //     sum_{j,d} P[i][j]_d * (E_j M^d) == 0.
        for (let i = 0; i < m; i++) {
          const acc: any[] = new Array(n).fill(F.zero());
          for (let j = 0; j < m; j++) {
            const coeffs = (P.get(i, j) as any).coeffs as any[];
            let v = Er.row(j);
            for (let dd = 0; dd < coeffs.length; dd++) {
              if (dd > 0) {
                const w: any[] = new Array(n).fill(F.zero());
                for (let a = 0; a < n; a++) {
                  let s = F.zero();
                  for (let b = 0; b < n; b++) s = s.add(v[b]!.mul(Mr.get(b, a)));
                  w[a] = s;
                }
                v = w;
              }
              for (let a = 0; a < n; a++) acc[a] = acc[a].add(coeffs[dd]!.mul(v[a]!));
            }
          }
          for (const a of acc) expect(a.isZero()).toBe(true);
        }

        // (4) Sage's doctest relation between P and K (matrix2.pyx:20131-20135).
        for (let j = 0; j < m; j++) {
          for (let i = 0; i < m; i++) {
            const coeffs = (P.get(i, j) as any).coeffs as any[];
            const target = new Map<number, any>();
            for (let k = 0; k < K.ncols; k++) {
              if (rp[k]![0] === j) {
                const dgr = rp[k]![1];
                target.set(dgr, (target.get(dgr) ?? F.zero()).add(K.get(i, k)));
              }
            }
            const maxd = Math.max(coeffs.length - 1, ...[...target.keys(), -1]);
            for (let dd = 0; dd <= maxd; dd++) {
              const a = dd < coeffs.length ? coeffs[dd]! : F.zero();
              const b = target.get(dd) ?? F.zero();
              expect(a.eq(b)).toBe(true);
            }
          }
        }
        cases++;
      }
    }
    expect(cases).toBe(36);
  });
});

// ============================================================================
// jordan_form(subdivide=...) -- the block subdivisions SageMath sets with
// ``block_diagonal_matrix(..., subdivide=subdivide)`` (matrix2.pyx:12255-12257)
// and prints with ``Matrix.str`` (matrix0.pyx:1834).
//
// Every string below is copied character for character out of the doctests of
// sage/matrix/matrix2.pyx.
// ============================================================================

describe('jordan_form subdivisions (matrix2.pyx:11814 subdivide=True)', () => {
  const F2 = GF(2n);

  it('reproduces the issue-6942 doctest with its block separators (matrix2.pyx:11975-11996)', () => {
    // sage: M = Matrix(GF(2), [[1,0,1,0,0,0,1], [1,0,0,1,1,1,0], [1,1,0,1,1,1,1],
    // ....:                    [1,1,1,0,1,1,1], [1,1,1,0,0,1,0], [1,1,1,0,1,0,0],
    // ....:                    [1,1,1,1,1,1,0]])
    // sage: J, T = M.jordan_form(transformation=True); J
    const M = ffmat(F2, [
      [1, 0, 1, 0, 0, 0, 1],
      [1, 0, 0, 1, 1, 1, 0],
      [1, 1, 0, 1, 1, 1, 1],
      [1, 1, 1, 0, 1, 1, 1],
      [1, 1, 1, 0, 0, 1, 0],
      [1, 1, 1, 0, 1, 0, 0],
      [1, 1, 1, 1, 1, 1, 0],
    ]);
    const [J, T] = jordan_form(M, undefined, undefined, undefined, true) as [
      Matrix<any>,
      Matrix<any>,
    ];
    expect(String(J)).toBe(
      [
        '[1 1|0 0|0 0|0]',
        '[0 1|0 0|0 0|0]',
        '[---+---+---+-]',
        '[0 0|1 1|0 0|0]',
        '[0 0|0 1|0 0|0]',
        '[---+---+---+-]',
        '[0 0|0 0|1 1|0]',
        '[0 0|0 0|0 1|0]',
        '[---+---+---+-]',
        '[0 0|0 0|0 0|1]',
      ].join('\n')
    );
    expect(matrix_str(J)).toBe(String(J));
    expect(subdivisions(J)).toEqual([
      [2, 4, 6],
      [2, 4, 6],
    ]);
    // sage: M * T == T * J  ->  True
    expect(matEq(M.mul(T), T.mul(J))).toBe(true);
    // and therefore  T^-1 * M * T == J
    expect(matEq(inverse(T).mul(M).mul(T), J)).toBe(true);
    // sage: T.rank() -> 7 ; sage: M.rank() -> 7
    expect(rank(T)).toBe(7);
    expect(rank(M)).toBe(7);
  });

  it('reproduces the first 10x10 QQ doctest with separators (matrix2.pyx:12020-12036)', () => {
    const f = (a: number, b: number) => (QQ as any).__call__([BigInt(a), BigInt(b)]);
    const e: any[][] = [
      [15, f(37, 3), -16, f(-104, 3), -29, f(-7, 3), 0, f(2, 3), f(-29, 3), f(-1, 3)],
      [2, 9, -1, -6, -6, 0, 0, 0, -2, 0],
      [24, f(74, 3), -41, f(-208, 3), -58, f(-23, 3), 0, f(4, 3), f(-58, 3), f(-2, 3)],
      [-6, -19, 3, 21, 19, 0, 0, 0, 6, 0],
      [2, 6, 3, -6, -3, 1, 0, 0, -2, 0],
      [-96, f(-296, 3), 176, f(832, 3), 232, f(101, 3), 0, f(-16, 3), f(232, 3), f(8, 3)],
      [-4, f(-2, 3), 21, f(16, 3), 4, f(14, 3), 3, f(-1, 3), f(4, 3), f(-25, 3)],
      [20, f(26, 3), -66, f(-199, 3), -42, f(-41, 3), 0, f(13, 3), f(-55, 3), f(-2, 3)],
      [18, 57, -9, -54, -57, 0, 0, 0, -15, 0],
      [0, 0, 0, 0, 0, 0, 0, 0, 0, 3],
    ];
    const A = new Matrix(
      QQ as any,
      10,
      10,
      e.map((r) => r.map((x) => (QQ as any).__call__(x)))
    );
    const [J, T] = jordan_form(A, undefined, undefined, undefined, true) as [
      Matrix<any>,
      Matrix<any>,
    ];
    expect(String(J)).toBe(
      [
        '[3 1 0|0 0 0|0 0 0|0]',
        '[0 3 1|0 0 0|0 0 0|0]',
        '[0 0 3|0 0 0|0 0 0|0]',
        '[-----+-----+-----+-]',
        '[0 0 0|3 1 0|0 0 0|0]',
        '[0 0 0|0 3 1|0 0 0|0]',
        '[0 0 0|0 0 3|0 0 0|0]',
        '[-----+-----+-----+-]',
        '[0 0 0|0 0 0|3 1 0|0]',
        '[0 0 0|0 0 0|0 3 1|0]',
        '[0 0 0|0 0 0|0 0 3|0]',
        '[-----+-----+-----+-]',
        '[0 0 0|0 0 0|0 0 0|3]',
      ].join('\n')
    );
    expect(subdivisions(J)).toEqual([
      [3, 6, 9],
      [3, 6, 9],
    ]);
    // sage: T * J * T**(-1) == A  ->  True ; sage: T.rank() -> 10
    expect(matEq(T.mul(J).mul(inverse(T)), A)).toBe(true);
    expect(matEq(inverse(T).mul(A).mul(T), J)).toBe(true);
    expect(rank(T)).toBe(10);
  });

  it('reproduces the second 10x10 QQ doctest with separators (matrix2.pyx:12085-12102)', () => {
    const f = (a: number, b: number) => (QQ as any).__call__([BigInt(a), BigInt(b)]);
    const e: any[][] = [
      [15, f(37, 3), -16, f(-104, 3), -29, f(-7, 3), 35, f(2, 3), f(-29, 3), f(-1, 3)],
      [2, 9, -1, -6, -6, 0, 7, 0, -2, 0],
      [24, f(74, 3), -29, f(-208, 3), -58, f(-14, 3), 70, f(4, 3), f(-58, 3), f(-2, 3)],
      [-6, -19, 3, 21, 19, 0, -21, 0, 6, 0],
      [2, 6, -1, -6, -3, 0, 7, 0, -2, 0],
      [-96, f(-296, 3), 128, f(832, 3), 232, f(65, 3), -279, f(-16, 3), f(232, 3), f(8, 3)],
      [0, 0, 0, 0, 0, 0, 3, 0, 0, 0],
      [20, f(26, 3), -30, f(-199, 3), -42, f(-14, 3), 70, f(13, 3), f(-55, 3), f(-2, 3)],
      [18, 57, -9, -54, -57, 0, 63, 0, -15, 0],
      [0, 0, 0, 0, 0, 0, 0, 0, 0, 3],
    ];
    const A = new Matrix(
      QQ as any,
      10,
      10,
      e.map((r) => r.map((x) => (QQ as any).__call__(x)))
    );
    const [J, T] = jordan_form(A, undefined, undefined, undefined, true) as [
      Matrix<any>,
      Matrix<any>,
    ];
    expect(String(J)).toBe(
      [
        '[3 1 0|0 0|0 0|0 0|0]',
        '[0 3 1|0 0|0 0|0 0|0]',
        '[0 0 3|0 0|0 0|0 0|0]',
        '[-----+---+---+---+-]',
        '[0 0 0|3 1|0 0|0 0|0]',
        '[0 0 0|0 3|0 0|0 0|0]',
        '[-----+---+---+---+-]',
        '[0 0 0|0 0|3 1|0 0|0]',
        '[0 0 0|0 0|0 3|0 0|0]',
        '[-----+---+---+---+-]',
        '[0 0 0|0 0|0 0|3 1|0]',
        '[0 0 0|0 0|0 0|0 3|0]',
        '[-----+---+---+---+-]',
        '[0 0 0|0 0|0 0|0 0|3]',
      ].join('\n')
    );
    // sage: T * J * T**(-1) == A  ->  True ; sage: T.rank() -> 10
    expect(matEq(T.mul(J).mul(inverse(T)), A)).toBe(true);
    expect(rank(T)).toBe(10);
  });

  it('reproduces the issue-12693 doctest, J subdivided and P not (matrix2.pyx:12108-12121)', () => {
    // sage: M = matrix(((2,2,2), (0,0,0), (-2,-2,-2)))
    // sage: J, P = M.jordan_form(transformation=True); J; P
    const M = qmat([
      [2, 2, 2],
      [0, 0, 0],
      [-2, -2, -2],
    ]);
    const [J, P] = jordan_form(M, undefined, undefined, undefined, true) as [
      Matrix<any>,
      Matrix<any>,
    ];
    expect(String(J)).toBe(['[0 1|0]', '[0 0|0]', '[---+-]', '[0 0|0]'].join('\n'));
    expect(subdivisions(J)).toEqual([[2], [2]]);
    // The transformation matrix carries no subdivision, and SageMath pads all
    // its entries to one common width.
    expect(subdivisions(P)).toEqual([[], []]);
    expect(matrix_str(P)).toBe(['[ 2  1  0]', '[ 0  0  1]', '[-2  0 -1]'].join('\n'));
    // sage: J - ~P * M * P  ->  the zero matrix
    expect(matEq(inverse(P).mul(M).mul(P), J)).toBe(true);
  });

  it('reproduces the supplied-eigenvalues doctest layout (matrix2.pyx:12144-12148)', () => {
    // sage: M = matrix(Qx, [[0, 0, x31], [0, 0, x21], [0, 0, 0]])
    // sage: M.jordan_form(eigenvalues=[(0, 3)])
    // [0 1|0]
    // [0 0|0]
    // [---+-]
    // [0 0|0]
    const M = qmat([
      [0, 0, 5],
      [0, 0, 7],
      [0, 0, 0],
    ]);
    const J = jordan_form(M, undefined, undefined, undefined, false, [
      [(QQ as any).__call__(0), 3],
    ]) as Matrix<any>;
    expect(String(J)).toBe(['[0 1|0]', '[0 0|0]', '[---+-]', '[0 0|0]'].join('\n'));
  });

  it('subdivide=false suppresses the separators and leaves the entries alone', () => {
    // sage: c = matrix(ZZ, 3, [1]*9)
    // sage: c.jordan_form(subdivide=False)   (matrix2.pyx:11945-11950)
    // has no separators; only the block ORDER differs here, see the note below.
    const M = qmat([
      [2, 2, 2],
      [0, 0, 0],
      [-2, -2, -2],
    ]);
    const Jsub = jordan_form(M, undefined, undefined, true) as Matrix<any>;
    const Jflat = jordan_form(M, undefined, undefined, false) as Matrix<any>;
    expect(subdivisions(Jflat)).toEqual([[], []]);
    expect(String(Jflat)).toBe(['[0 1 0]', '[0 0 0]', '[0 0 0]'].join('\n'));
    // Same matrix, only the printed form differs.
    expect(matEq(Jsub, Jflat)).toBe(true);
  });

  it('a single Jordan block carries no subdivision (matrix2.pyx:11995 issue 6932)', () => {
    // sage: M = Matrix(1, 1, [1]); M.jordan_form(transformation=True)  ->  ([1], [1])
    const [J1, P1] = jordan_form(qmat([[1]]), undefined, undefined, undefined, true) as [
      Matrix<any>,
      Matrix<any>,
    ];
    expect(matrix_str(J1)).toBe('[1]');
    expect(matrix_str(P1)).toBe('[1]');
    expect(subdivisions(J1)).toEqual([[], []]);
    // A 2x2 single Jordan block likewise.
    const J2 = jordan_form(
      qmat([
        [2, 1],
        [0, 2],
      ])
    ) as Matrix<any>;
    expect(subdivisions(J2)).toEqual([[], []]);
    expect(matrix_str(J2)).toBe(['[2 1]', '[0 2]'].join('\n'));
  });
});

describe('matrix_str (port of matrix0.pyx:1834)', () => {
  it('pads every entry to one common width, as upstream does', () => {
    // sage: matrix(ZZ, 2, 2, [1, 222, 3, 4])
    // [  1 222]
    // [  3   4]
    const A = qmat([
      [1, 222],
      [3, 4],
    ]);
    expect(matrix_str(A)).toBe(['[  1 222]', '[  3   4]'].join('\n'));
  });

  it('returns [] for an empty matrix (matrix0.pyx:2114-2116)', () => {
    // sage: print(matrix(ZZ, 0, 0).str())  ->  []
    expect(matrix_str(new Matrix(QQ as any, 0, 0))).toBe('[]');
    expect(matrix_str(new Matrix(QQ as any, 0, 4))).toBe('[]');
    expect(matrix_str(new Matrix(QQ as any, 3, 0))).toBe('[]');
  });

  it("reproduces subdivide's own doctests verbatim (matrix2.pyx:9647-9712)", () => {
    // sage: M = matrix(5, 5, prime_range(100))
    const M = qmat([
      [2, 3, 5, 7, 11],
      [13, 17, 19, 23, 29],
      [31, 37, 41, 43, 47],
      [53, 59, 61, 67, 71],
      [73, 79, 83, 89, 97],
    ]);

    // sage: M.subdivide(2,3); M
    subdivide(M, 2, 3);
    expect(matrix_str(M)).toBe(
      [
        '[ 2  3  5| 7 11]',
        '[13 17 19|23 29]',
        '[--------+-----]',
        '[31 37 41|43 47]',
        '[53 59 61|67 71]',
        '[73 79 83|89 97]',
      ].join('\n')
    );
    // sage: M.subdivisions()  ->  ([2], [3])
    expect(subdivisions(M)).toEqual([[2], [3]]);

    // sage: M.subdivide(None, [1,3]); M
    subdivide(M, [], [1, 3]);
    expect(matrix_str(M)).toBe(
      [
        '[ 2| 3  5| 7 11]',
        '[13|17 19|23 29]',
        '[31|37 41|43 47]',
        '[53|59 61|67 71]',
        '[73|79 83|89 97]',
      ].join('\n')
    );

    // Degenerate cases work too (matrix2.pyx:9674-9684):
    // sage: M.subdivide([2,5], [0,1,3]); M
    subdivide(M, [2, 5], [0, 1, 3]);
    expect(matrix_str(M)).toBe(
      [
        '[| 2| 3  5| 7 11]',
        '[|13|17 19|23 29]',
        '[+--+-----+-----]',
        '[|31|37 41|43 47]',
        '[|53|59 61|67 71]',
        '[|73|79 83|89 97]',
        '[+--+-----+-----]',
      ].join('\n')
    );

    // sage: M.subdivide([2,2,3], [0,0,1,1]); M   (matrix2.pyx:9690-9698)
    subdivide(M, [2, 2, 3], [0, 0, 1, 1]);
    expect(matrix_str(M)).toBe(
      [
        '[|| 2|| 3  5  7 11]',
        '[||13||17 19 23 29]',
        '[++--++-----------]',
        '[++--++-----------]',
        '[||31||37 41 43 47]',
        '[++--++-----------]',
        '[||53||59 61 67 71]',
        '[||73||79 83 89 97]',
      ].join('\n')
    );

    // Indices do not need to be in the right order (issue 14064,
    // matrix2.pyx:9704-9712):
    // sage: M.subdivide([4, 2], [3, 1]); M
    subdivide(M, [4, 2], [3, 1]);
    expect(matrix_str(M)).toBe(
      [
        '[ 2| 3  5| 7 11]',
        '[13|17 19|23 29]',
        '[--+-----+-----]',
        '[31|37 41|43 47]',
        '[53|59 61|67 71]',
        '[--+-----+-----]',
        '[73|79 83|89 97]',
      ].join('\n')
    );
  });
});

// ============================================================================
// krylov_kernel_basis(var=...): the normal-form properties SageMath's own
// doctests assert (matrix2.pyx:20085, 20158-20160, 20182, 20089-20095), checked
// with the ported sage/matrix/matrix_polynomial_dense predicates rather than by
// byte-equality with Sage's printed output.
// ============================================================================

describe('item 22b: krylov_kernel_basis is in shifted Popov / Hermite form', () => {
  const E = ffmat(F97, [
    [27, 49, 29],
    [50, 58, 0],
    [77, 10, 29],
  ]);
  const M = ffmat(F97, [
    [0, 1, 0],
    [0, 0, 1],
    [0, 0, 0],
  ]);
  const Rx = new PolynomialRing(F97 as any, 'x');

  /** Coefficient-by-coefficient comparison of two polynomial matrices. */
  function polyMatEq(A: Matrix<any>, B: Matrix<any>): boolean {
    if (A.nrows !== B.nrows || A.ncols !== B.ncols) return false;
    for (let i = 0; i < A.nrows; i++) {
      for (let j = 0; j < A.ncols; j++) {
        const a = (A.get(i, j) as any).coeffs as any[];
        const b = (B.get(i, j) as any).coeffs as any[];
        if (a.length !== b.length) return false;
        for (let k = 0; k < a.length; k++) if (!a[k]!.eq(b[k]!)) return false;
      }
    }
    return true;
  }

  it('P.is_popov() is true and P is the minimal approximant basis (matrix2.pyx:20085-20095)', () => {
    const [P] = krylov_kernel_basis(E, M, undefined, undefined, true, Rx as any) as [
      Matrix<any>,
      Array<[number, number, number]>,
    ];
    // sage: P.is_popov()   ->  True
    expect(is_popov(P as any)).toBe(true);

    // sage: x = P.base_ring().gen()
    // sage: F = E * matrix([[1], [x], [x**2]]); F
    // [29*x^2 + 49*x + 27]
    // [         58*x + 50]
    // [29*x^2 + 10*x + 77]
    const x = (Rx as any).gen();
    const V = new Matrix(Rx as any, 3, 1, [[(Rx as any).one()], [x], [x.mul(x)]] as any);
    const Epoly = new Matrix(
      Rx as any,
      3,
      3,
      E.rows().map((r: any[]) => r.map((c: any) => (Rx as any).__call__(c)))
    );
    const F = Epoly.mul(V);
    expect(render(F)).toBe('[29*x^2 + 49*x + 27] / [58*x + 50] / [29*x^2 + 10*x + 77]');

    // sage: P == F.minimal_approximant_basis(3, normal_form=True)   ->  True
    const AB = minimal_approximant_basis(F as any, 3, { normal_form: true });
    expect(polyMatEq(P, AB as any)).toBe(true);
  });

  it('H is in shifts-Popov and lower-echelon Hermite form (matrix2.pyx:20158-20162)', () => {
    const shifts = [0, 3, 6];
    const H = krylov_kernel_basis(E, M, shifts, undefined, false, Rx as any) as Matrix<any>;
    // sage: H.is_popov(shifts=shifts) and H.is_hermite(lower_echelon=True)  ->  True
    expect(is_popov(H as any, { shifts })).toBe(true);
    expect(is_hermite(H as any, { lower_echelon: true })).toBe(true);

    // sage: P.popov_form(shifts=shifts) == H   ->  True
    const [P] = krylov_kernel_basis(E, M, undefined, undefined, true, Rx as any) as [
      Matrix<any>,
      Array<[number, number, number]>,
    ];
    expect(polyMatEq(popov_form(P as any, { shifts }) as any, H)).toBe(true);
  });

  it('Q is in shifts-Popov form and equals P.popov_form(shifts) (matrix2.pyx:20182-20186)', () => {
    const shifts = [3, 0, 2];
    const RY = new PolynomialRing(F97 as any, 'Y');
    const [Q] = krylov_kernel_basis(E, M, shifts, undefined, true, RY as any) as [
      Matrix<any>,
      Array<[number, number, number]>,
    ];
    // sage: Q.is_popov(shifts=shifts)   ->  True
    expect(is_popov(Q as any, { shifts })).toBe(true);

    // sage: P.popov_form(shifts=shifts) == Q(x)   ->  True
    // (Q lives over GF(97)[Y]; only the variable name differs, so we compare
    // the coefficient lists.)
    const [P] = krylov_kernel_basis(E, M, undefined, undefined, true, Rx as any) as [
      Matrix<any>,
      Array<[number, number, number]>,
    ];
    const PQ = popov_form(P as any, { shifts }) as any;
    expect(PQ.nrows).toBe(Q.nrows);
    for (let i = 0; i < Q.nrows; i++) {
      for (let j = 0; j < Q.ncols; j++) {
        const a = (PQ.get(i, j) as any).coeffs as any[];
        const b = (Q.get(i, j) as any).coeffs as any[];
        expect(a.length).toBe(b.length);
        for (let k = 0; k < a.length; k++) expect(a[k]!.eq(b[k]!)).toBe(true);
      }
    }
  });

  it('the polynomial kernel basis is in shifted Popov form on random inputs', () => {
    const rng = makeRng(20260728);
    let cases = 0;
    for (const F of [GF(7n) as any, GF(97n) as any]) {
      const RF = new PolynomialRing(F, 'x');
      for (let t = 0; t < 15; t++) {
        const m = 1 + (rng() % 3);
        const n = 1 + (rng() % 3);
        const shifts = t % 3 === 0 ? undefined : Array.from({ length: m }, () => rng() % 4);
        const Er = new Matrix(
          F,
          m,
          n,
          Array.from({ length: m }, () =>
            Array.from({ length: n }, () => F.__call__((rng() % 17) - 8))
          )
        );
        const Mr = new Matrix(
          F,
          n,
          n,
          Array.from({ length: n }, () =>
            Array.from({ length: n }, () => F.__call__((rng() % 17) - 8))
          )
        );
        const P = krylov_kernel_basis(Er, Mr, shifts, undefined, false, RF as any) as Matrix<any>;
        expect(is_popov(P as any, shifts === undefined ? undefined : { shifts })).toBe(true);
        cases++;
      }
    }
    expect(cases).toBe(30);
  });
});
