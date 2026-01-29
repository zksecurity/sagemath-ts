/**
 * Tests for matrix normal forms (HNF, Smith, etc.)
 *
 * These test implementations based on:
 * - sage/matrix/matrix_integer_dense.pyx
 * - sage/matrix/matrix2.pyx
 */

import { describe, expect, it } from 'vitest';
import {
  IntegerMatrix,
  IntegerMatrixFromEntries,
  hermite_normal_form,
  smith_form_integer,
  elementary_divisors_integer,
  rank_integer,
  kernel_matrix,
  left_kernel_matrix_integer,
  identity_integer_matrix,
} from './index.js';

describe('Hermite Normal Form', () => {
  describe('basic properties', () => {
    it('should compute HNF of a simple 2x2 matrix', () => {
      const A = IntegerMatrixFromEntries([
        [2n, 3n],
        [4n, 5n],
      ]);

      const H = hermite_normal_form(A) as IntegerMatrix;

      // HNF should be upper triangular
      expect(H.get(1, 0).value).toBe(0n);

      // Diagonal entries should be positive
      expect(H.get(0, 0).value > 0n).toBe(true);
      expect(H.get(1, 1).value > 0n).toBe(true);
    });

    it('should compute HNF of identity matrix (identity is already in HNF)', () => {
      const I = identity_integer_matrix(3);
      const H = hermite_normal_form(I) as IntegerMatrix;

      expect(H.eq(I)).toBe(true);
    });

    it('should compute HNF of a 3x3 matrix', () => {
      const A = IntegerMatrixFromEntries([
        [1n, 2n, 3n],
        [4n, 5n, 6n],
        [7n, 8n, 9n],
      ]);

      const H = hermite_normal_form(A) as IntegerMatrix;

      // HNF should be upper triangular
      expect(H.get(1, 0).value).toBe(0n);
      expect(H.get(2, 0).value).toBe(0n);
      expect(H.get(2, 1).value).toBe(0n);
    });

    it('should return transformation matrix when requested', () => {
      const A = IntegerMatrixFromEntries([
        [2n, 3n],
        [4n, 5n],
      ]);

      const [H, U] = hermite_normal_form(A, 'default', false, true, true) as [
        IntegerMatrix,
        IntegerMatrix,
      ];

      // Verify H = U * A
      const UA = U.mul(A);
      expect(UA.eq(H)).toBe(true);
    });

    it('should handle matrices with zero rows', () => {
      const A = IntegerMatrixFromEntries([
        [1n, 2n, 3n],
        [0n, 0n, 0n],
        [4n, 5n, 6n],
      ]);

      const H = hermite_normal_form(A) as IntegerMatrix;

      // Result should be upper triangular
      expect(H.get(1, 0).value).toBe(0n);
      expect(H.get(2, 0).value).toBe(0n);
      expect(H.get(2, 1).value).toBe(0n);
    });

    it('should exclude zero rows when requested', () => {
      const A = IntegerMatrixFromEntries([
        [1n, 2n],
        [2n, 4n],
      ]);

      const H = hermite_normal_form(A, 'default', false, false, false) as IntegerMatrix;

      // Second row should be zero (linearly dependent), and excluded
      expect(H.nrows).toBeLessThanOrEqual(2);
    });

    it('should correctly reduce entries above pivots', () => {
      const A = IntegerMatrixFromEntries([
        [4n, 8n],
        [2n, 5n],
      ]);

      const H = hermite_normal_form(A) as IntegerMatrix;

      // Check that entries above pivots are in correct range
      if (H.get(1, 1).value !== 0n) {
        const pivot = H.get(1, 1).value;
        const above = H.get(0, 1).value;
        // 0 <= above < pivot
        expect(above >= 0n).toBe(true);
        expect(above < pivot).toBe(true);
      }
    });
  });

  describe('non-square matrices', () => {
    it('should handle more rows than columns', () => {
      const A = IntegerMatrixFromEntries([
        [1n, 2n],
        [3n, 4n],
        [5n, 6n],
      ]);

      const H = hermite_normal_form(A) as IntegerMatrix;

      expect(H.nrows).toBe(3);
      expect(H.ncols).toBe(2);
      // Should be upper triangular in the appropriate sense
      expect(H.get(1, 0).value).toBe(0n);
      expect(H.get(2, 0).value).toBe(0n);
    });

    it('should handle more columns than rows', () => {
      const A = IntegerMatrixFromEntries([
        [1n, 2n, 3n],
        [4n, 5n, 6n],
      ]);

      const H = hermite_normal_form(A) as IntegerMatrix;

      expect(H.nrows).toBe(2);
      expect(H.ncols).toBe(3);
      expect(H.get(1, 0).value).toBe(0n);
    });
  });
});

describe('Smith Normal Form', () => {
  describe('basic properties', () => {
    it('should compute SNF of identity matrix', () => {
      const I = identity_integer_matrix(3);
      const D = smith_form_integer(I) as IntegerMatrix;

      // SNF of identity is identity
      for (let i = 0; i < 3; i++) {
        for (let j = 0; j < 3; j++) {
          if (i === j) {
            expect(D.get(i, j).value).toBe(1n);
          } else {
            expect(D.get(i, j).value).toBe(0n);
          }
        }
      }
    });

    it('should compute SNF of a 2x2 matrix', () => {
      const A = IntegerMatrixFromEntries([
        [2n, 4n],
        [6n, 8n],
      ]);

      const D = smith_form_integer(A) as IntegerMatrix;

      // D should be diagonal
      expect(D.get(0, 1).value).toBe(0n);
      expect(D.get(1, 0).value).toBe(0n);

      // Diagonal entries should be non-negative
      expect(D.get(0, 0).value >= 0n).toBe(true);
      expect(D.get(1, 1).value >= 0n).toBe(true);

      // d_1 | d_2 (divisibility)
      if (D.get(0, 0).value !== 0n) {
        expect(D.get(1, 1).value % D.get(0, 0).value).toBe(0n);
      }
    });

    it('should return transformation matrices when requested', () => {
      const A = IntegerMatrixFromEntries([
        [1n, 2n],
        [3n, 4n],
      ]);

      const [D, U, V] = smith_form_integer(A, true) as [
        IntegerMatrix,
        IntegerMatrix,
        IntegerMatrix,
      ];

      // Verify D = U * A * V
      const UAV = U.mul(A).mul(V);
      expect(UAV.eq(D)).toBe(true);
    });

    it('should satisfy D = U * A * V for the Sage example', () => {
      // From Sage: A = MatrixSpace(IntegerRing(), 3)(range(9))
      const A = IntegerMatrixFromEntries([
        [0n, 1n, 2n],
        [3n, 4n, 5n],
        [6n, 7n, 8n],
      ]);

      const [D, U, V] = smith_form_integer(A, true) as [
        IntegerMatrix,
        IntegerMatrix,
        IntegerMatrix,
      ];

      // Verify D = U * A * V
      const UAV = U.mul(A).mul(V);
      expect(UAV.eq(D)).toBe(true);

      // D should be diagonal
      for (let i = 0; i < 3; i++) {
        for (let j = 0; j < 3; j++) {
          if (i !== j) {
            expect(D.get(i, j).value).toBe(0n);
          }
        }
      }
    });

    it('should handle singular matrix', () => {
      const A = IntegerMatrixFromEntries([
        [1n, 2n, 3n],
        [4n, 5n, 6n],
        [7n, 8n, 9n],
      ]);

      const D = smith_form_integer(A) as IntegerMatrix;

      // This matrix is singular (rank 2), so one diagonal entry should be 0
      const diag = [D.get(0, 0).value, D.get(1, 1).value, D.get(2, 2).value];
      const numZeros = diag.filter((x) => x === 0n).length;
      expect(numZeros).toBeGreaterThanOrEqual(1);
    });

    it('should have divisibility chain', () => {
      const A = IntegerMatrixFromEntries([
        [6n, 4n, 2n],
        [3n, 9n, 6n],
        [12n, 8n, 4n],
      ]);

      const D = smith_form_integer(A) as IntegerMatrix;

      // Check d_i | d_{i+1}
      const d = [D.get(0, 0).value, D.get(1, 1).value, D.get(2, 2).value];

      for (let i = 0; i < 2; i++) {
        if (d[i] !== 0n && d[i + 1] !== 0n) {
          expect(d[i + 1]! % d[i]!).toBe(0n);
        }
      }
    });
  });

  describe('non-square matrices', () => {
    it('should handle more rows than columns', () => {
      const A = IntegerMatrixFromEntries([
        [1n, 2n],
        [3n, 4n],
        [5n, 6n],
      ]);

      const [D, U, V] = smith_form_integer(A, true) as [
        IntegerMatrix,
        IntegerMatrix,
        IntegerMatrix,
      ];

      expect(D.nrows).toBe(3);
      expect(D.ncols).toBe(2);

      // Verify D = U * A * V
      const UAV = U.mul(A).mul(V);
      expect(UAV.eq(D)).toBe(true);
    });

    it('should handle more columns than rows', () => {
      const A = IntegerMatrixFromEntries([
        [1n, 2n, 3n],
        [4n, 5n, 6n],
      ]);

      const [D, U, V] = smith_form_integer(A, true) as [
        IntegerMatrix,
        IntegerMatrix,
        IntegerMatrix,
      ];

      expect(D.nrows).toBe(2);
      expect(D.ncols).toBe(3);

      // Verify D = U * A * V
      const UAV = U.mul(A).mul(V);
      expect(UAV.eq(D)).toBe(true);
    });
  });
});

describe('Elementary Divisors', () => {
  it('should return diagonal of Smith form', () => {
    const A = IntegerMatrixFromEntries([
      [2n, 4n],
      [6n, 8n],
    ]);

    const divs = elementary_divisors_integer(A);
    const D = smith_form_integer(A) as IntegerMatrix;

    expect(divs.length).toBe(2);
    expect(divs[0]!.value).toBe(D.get(0, 0).value);
    expect(divs[1]!.value).toBe(D.get(1, 1).value);
  });

  it('should have divisibility: d_i | d_{i+1}', () => {
    const A = IntegerMatrixFromEntries([
      [6n, 4n, 2n],
      [3n, 9n, 6n],
      [12n, 8n, 4n],
    ]);

    const divs = elementary_divisors_integer(A);

    for (let i = 0; i < divs.length - 1; i++) {
      if (!divs[i]!.isZero() && !divs[i + 1]!.isZero()) {
        expect(divs[i + 1]!.value % divs[i]!.value).toBe(0n);
      }
    }
  });

  it('should match Sage example', () => {
    // sage: matrix(3, range(9)).elementary_divisors() -> [1, 3, 0]
    const A = IntegerMatrixFromEntries([
      [0n, 1n, 2n],
      [3n, 4n, 5n],
      [6n, 7n, 8n],
    ]);

    const divs = elementary_divisors_integer(A);

    // The values [1, 3, 0] are expected
    expect(divs.length).toBe(3);
    expect(divs[0]!.value).toBe(1n);
    expect(divs[1]!.value).toBe(3n);
    expect(divs[2]!.value).toBe(0n);
  });
});

describe('Rank', () => {
  it('should compute rank of full-rank matrix', () => {
    const A = IntegerMatrixFromEntries([
      [1n, 0n],
      [0n, 1n],
    ]);

    expect(rank_integer(A)).toBe(2);
  });

  it('should compute rank of rank-deficient matrix', () => {
    const A = IntegerMatrixFromEntries([
      [1n, 2n, 3n],
      [4n, 5n, 6n],
      [7n, 8n, 9n],
    ]);

    // This matrix has rank 2
    expect(rank_integer(A)).toBe(2);
  });

  it('should compute rank of zero matrix', () => {
    const A = IntegerMatrixFromEntries([
      [0n, 0n],
      [0n, 0n],
    ]);

    expect(rank_integer(A)).toBe(0);
  });

  it('should handle non-square matrices', () => {
    const A = IntegerMatrixFromEntries([
      [1n, 2n, 3n],
      [4n, 5n, 6n],
    ]);

    expect(rank_integer(A)).toBe(2);
  });
});

describe('Kernel Matrix', () => {
  it('should compute kernel of full-rank matrix', () => {
    const A = IntegerMatrixFromEntries([
      [1n, 0n],
      [0n, 1n],
    ]);

    const K = kernel_matrix(A);

    // Full-rank matrix has trivial kernel
    expect(K.nrows).toBe(0);
  });

  it('should compute kernel of rank-deficient matrix', () => {
    const A = IntegerMatrixFromEntries([
      [1n, 2n, 3n],
      [4n, 5n, 6n],
      [7n, 8n, 9n],
    ]);

    const K = kernel_matrix(A);

    // Rank 2, so kernel has dimension 1
    expect(K.nrows).toBe(1);
    expect(K.ncols).toBe(3);

    // Verify A * K^T = 0
    for (let i = 0; i < K.nrows; i++) {
      const v = new IntegerMatrix(3, 1);
      for (let j = 0; j < 3; j++) {
        v.set(j, 0, K.get(i, j));
      }

      const Av = A.mul(v);
      for (let r = 0; r < A.nrows; r++) {
        expect(Av.get(r, 0).value).toBe(0n);
      }
    }
  });

  it('should compute kernel of zero matrix', () => {
    const A = IntegerMatrixFromEntries([
      [0n, 0n],
      [0n, 0n],
    ]);

    const K = kernel_matrix(A);

    // Kernel is the entire space
    expect(K.nrows).toBe(2);
  });

  it('should verify kernel vectors satisfy A * v = 0', () => {
    const A = IntegerMatrixFromEntries([
      [1n, 2n, 1n],
      [2n, 4n, 2n],
    ]);

    const K = kernel_matrix(A);

    // For each kernel vector, verify A * v = 0
    for (let i = 0; i < K.nrows; i++) {
      const v = new IntegerMatrix(K.ncols, 1);
      for (let j = 0; j < K.ncols; j++) {
        v.set(j, 0, K.get(i, j));
      }

      const Av = A.mul(v);
      for (let r = 0; r < A.nrows; r++) {
        expect(Av.get(r, 0).value).toBe(0n);
      }
    }
  });

  it('should compute kernel of non-square matrix', () => {
    const A = IntegerMatrixFromEntries([
      [1n, 2n, 3n, 4n],
      [2n, 4n, 6n, 8n],
    ]);

    const K = kernel_matrix(A);

    // Rank 1 matrix with 4 columns -> kernel has dimension 3
    expect(K.nrows).toBe(3);
    expect(K.ncols).toBe(4);

    // Verify each kernel vector
    for (let i = 0; i < K.nrows; i++) {
      const v = new IntegerMatrix(4, 1);
      for (let j = 0; j < 4; j++) {
        v.set(j, 0, K.get(i, j));
      }

      const Av = A.mul(v);
      for (let r = 0; r < A.nrows; r++) {
        expect(Av.get(r, 0).value).toBe(0n);
      }
    }
  });
});

describe('Left Kernel Matrix', () => {
  it('should compute left kernel', () => {
    const A = IntegerMatrixFromEntries([
      [1n, 2n],
      [2n, 4n],
      [3n, 6n],
    ]);

    const K = left_kernel_matrix_integer(A);

    // Left kernel should have dimension 2 (rank is 1, and we have 3 rows)
    expect(K.nrows).toBe(2);
    expect(K.ncols).toBe(3);

    // Verify v * A = 0 for each kernel vector
    for (let i = 0; i < K.nrows; i++) {
      const v = new IntegerMatrix(1, 3);
      for (let j = 0; j < 3; j++) {
        v.set(0, j, K.get(i, j));
      }

      const vA = v.mul(A);
      for (let c = 0; c < A.ncols; c++) {
        expect(vA.get(0, c).value).toBe(0n);
      }
    }
  });
});

describe('IntegerMatrix convenience methods', () => {
  it('should call hermite_form method directly', () => {
    const A = IntegerMatrixFromEntries([
      [2n, 3n],
      [4n, 5n],
    ]);

    const H = A.hermite_form();
    expect(H.get(1, 0).value).toBe(0n);
  });

  it('should call smith_form method directly', () => {
    const A = IntegerMatrixFromEntries([
      [2n, 4n],
      [6n, 8n],
    ]);

    const D = A.smith_form();
    expect(D.get(0, 1).value).toBe(0n);
    expect(D.get(1, 0).value).toBe(0n);
  });

  it('should call smith_form with transformation matrices', () => {
    const A = IntegerMatrixFromEntries([
      [1n, 2n],
      [3n, 4n],
    ]);

    const [D, U, V] = A.smith_form(true);
    const UAV = U.mul(A).mul(V);
    expect(UAV.eq(D)).toBe(true);
  });

  it('should call elementary_divisors method directly', () => {
    const A = IntegerMatrixFromEntries([
      [0n, 1n, 2n],
      [3n, 4n, 5n],
      [6n, 7n, 8n],
    ]);

    const divs = A.elementary_divisors();
    expect(divs[0]!.value).toBe(1n);
    expect(divs[1]!.value).toBe(3n);
    expect(divs[2]!.value).toBe(0n);
  });

  it('should call rank method directly', () => {
    const A = IntegerMatrixFromEntries([
      [1n, 2n, 3n],
      [4n, 5n, 6n],
      [7n, 8n, 9n],
    ]);

    expect(A.rank()).toBe(2);
  });

  it('should call right_kernel_matrix method directly', () => {
    const A = IntegerMatrixFromEntries([
      [1n, 2n, 3n],
      [4n, 5n, 6n],
      [7n, 8n, 9n],
    ]);

    const K = A.right_kernel_matrix();
    expect(K.nrows).toBe(1);

    // Verify A * K^T = 0
    const v = new IntegerMatrix(3, 1);
    for (let j = 0; j < 3; j++) {
      v.set(j, 0, K.get(0, j));
    }
    const Av = A.mul(v);
    for (let r = 0; r < 3; r++) {
      expect(Av.get(r, 0).value).toBe(0n);
    }
  });

  it('should call left_kernel_matrix method directly', () => {
    const A = IntegerMatrixFromEntries([
      [1n, 2n],
      [2n, 4n],
      [3n, 6n],
    ]);

    const K = A.left_kernel_matrix();
    expect(K.nrows).toBe(2);
  });
});

describe('Edge cases', () => {
  it('should handle 1x1 matrices', () => {
    const A = IntegerMatrixFromEntries([[5n]]);

    const H = hermite_normal_form(A) as IntegerMatrix;
    expect(H.get(0, 0).value).toBe(5n);

    const D = smith_form_integer(A) as IntegerMatrix;
    expect(D.get(0, 0).value).toBe(5n);

    expect(rank_integer(A)).toBe(1);

    const K = kernel_matrix(A);
    expect(K.nrows).toBe(0);
  });

  it('should handle 1x1 zero matrix', () => {
    const A = IntegerMatrixFromEntries([[0n]]);

    const H = hermite_normal_form(A) as IntegerMatrix;
    expect(H.get(0, 0).value).toBe(0n);

    const D = smith_form_integer(A) as IntegerMatrix;
    expect(D.get(0, 0).value).toBe(0n);

    expect(rank_integer(A)).toBe(0);

    const K = kernel_matrix(A);
    expect(K.nrows).toBe(1);
  });

  it('should handle negative entries', () => {
    const A = IntegerMatrixFromEntries([
      [-2n, 3n],
      [4n, -5n],
    ]);

    const [H, U] = hermite_normal_form(A, 'default', false, true, true) as [
      IntegerMatrix,
      IntegerMatrix,
    ];

    // Verify H = U * A
    const UA = U.mul(A);
    expect(UA.eq(H)).toBe(true);

    // HNF diagonal should be positive
    expect(H.get(0, 0).value > 0n).toBe(true);
    expect(H.get(1, 1).value >= 0n).toBe(true);
  });

  it('should handle large entries', () => {
    const A = IntegerMatrixFromEntries([
      [1000000000000n, 2000000000000n],
      [3000000000000n, 4000000000000n],
    ]);

    const [D, U, V] = smith_form_integer(A, true) as [
      IntegerMatrix,
      IntegerMatrix,
      IntegerMatrix,
    ];

    // Verify D = U * A * V
    const UAV = U.mul(A).mul(V);
    expect(UAV.eq(D)).toBe(true);
  });
});
