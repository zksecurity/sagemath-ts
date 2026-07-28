/**
 * Tests for GF(2) matrix operations
 * @see Reference: sage/matrix/matrix_mod2_dense.pyx
 */

import { describe, expect, it } from 'vitest';
import {
  Matrix_mod2_dense,
  from_png_data,
  identity_matrix_gf2,
  matrix_gf2_from_entries,
  ple,
  pluq,
  to_png_data,
  zero_matrix_gf2,
} from './matrix_mod2.js';

/** The rows of a GF(2) matrix, as plain arrays. */
function rowsOf(A: Matrix_mod2_dense): number[][] {
  const rows: number[][] = [];
  for (let i = 0; i < A.nrows; i++) rows.push(A.row(i));
  return rows;
}

describe('Matrix_mod2_dense', () => {
  describe('creation and basic operations', () => {
    it('should create a zero matrix', () => {
      const A = zero_matrix_gf2(3, 4);
      expect(A.nrows).toBe(3);
      expect(A.ncols).toBe(4);
      for (let i = 0; i < 3; i++) {
        for (let j = 0; j < 4; j++) {
          expect(A.get(i, j)).toBe(0);
        }
      }
    });

    it('should create an identity matrix', () => {
      const I = identity_matrix_gf2(3);
      expect(I.nrows).toBe(3);
      expect(I.ncols).toBe(3);
      for (let i = 0; i < 3; i++) {
        for (let j = 0; j < 3; j++) {
          expect(I.get(i, j)).toBe(i === j ? 1 : 0);
        }
      }
    });

    it('should create matrix from entries', () => {
      const A = matrix_gf2_from_entries([
        [1, 0, 1],
        [0, 1, 0],
      ]);
      expect(A.nrows).toBe(2);
      expect(A.ncols).toBe(3);
      expect(A.get(0, 0)).toBe(1);
      expect(A.get(0, 1)).toBe(0);
      expect(A.get(0, 2)).toBe(1);
      expect(A.get(1, 0)).toBe(0);
      expect(A.get(1, 1)).toBe(1);
      expect(A.get(1, 2)).toBe(0);
    });

    it('should reduce entries mod 2', () => {
      const A = matrix_gf2_from_entries([
        [2, 3, 4],
        [5, 6, 7],
      ]);
      expect(A.get(0, 0)).toBe(0); // 2 mod 2
      expect(A.get(0, 1)).toBe(1); // 3 mod 2
      expect(A.get(0, 2)).toBe(0); // 4 mod 2
      expect(A.get(1, 0)).toBe(1); // 5 mod 2
      expect(A.get(1, 1)).toBe(0); // 6 mod 2
      expect(A.get(1, 2)).toBe(1); // 7 mod 2
    });
  });

  describe('swap_rows and swap_columns', () => {
    it('should swap rows correctly', () => {
      const A = matrix_gf2_from_entries([
        [1, 0],
        [0, 1],
      ]);
      A.swap_rows(0, 1);
      expect(A.get(0, 0)).toBe(0);
      expect(A.get(0, 1)).toBe(1);
      expect(A.get(1, 0)).toBe(1);
      expect(A.get(1, 1)).toBe(0);
    });

    it('should swap columns correctly', () => {
      const A = matrix_gf2_from_entries([
        [1, 0],
        [0, 1],
      ]);
      A.swap_columns(0, 1);
      expect(A.get(0, 0)).toBe(0);
      expect(A.get(0, 1)).toBe(1);
      expect(A.get(1, 0)).toBe(1);
      expect(A.get(1, 1)).toBe(0);
    });

    it('should handle self-swap (no-op)', () => {
      const A = matrix_gf2_from_entries([
        [1, 0],
        [0, 1],
      ]);
      const Acopy = A.copy();
      A.swap_rows(0, 0);
      expect(A.eq(Acopy)).toBe(true);

      A.swap_columns(1, 1);
      expect(A.eq(Acopy)).toBe(true);
    });
  });

  describe('doubly_lexical_ordering', () => {
    it('should order a 2x2 matrix', () => {
      // Example from SageMath docs
      const A = matrix_gf2_from_entries([
        [0, 1],
        [1, 0],
      ]);
      const [rowPerm, colPerm] = A.doubly_lexical_ordering();

      // After applying the permutation, the result should be lexically ordered
      // Create a copy and apply permutation
      const B = A.copy();
      B.permute_rows(rowPerm);
      B.permute_columns(colPerm);

      // The result should be identity matrix (lexically ordered)
      expect(B.get(0, 0)).toBe(1);
      expect(B.get(0, 1)).toBe(0);
      expect(B.get(1, 0)).toBe(0);
      expect(B.get(1, 1)).toBe(1);
    });

    it('should order a 2x2 matrix inplace', () => {
      const A = matrix_gf2_from_entries([
        [0, 1],
        [1, 0],
      ]);
      A.doubly_lexical_ordering(true);

      // After inplace ordering, A should be lexically ordered
      expect(A.get(0, 0)).toBe(1);
      expect(A.get(0, 1)).toBe(0);
      expect(A.get(1, 0)).toBe(0);
      expect(A.get(1, 1)).toBe(1);
    });

    it('should return identity permutation for empty matrix', () => {
      const A = new Matrix_mod2_dense(0, 2);
      const [rowPerm, colPerm] = A.doubly_lexical_ordering();
      expect(rowPerm).toEqual([]);
      expect(colPerm).toEqual([1, 2]);
    });

    it('should handle a larger matrix (Example 3.7 from HAM1985)', () => {
      const A = matrix_gf2_from_entries([
        [1, 1, 0, 0, 0, 0, 0],
        [1, 1, 0, 0, 0, 0, 0],
        [1, 1, 0, 1, 0, 0, 0],
        [0, 0, 1, 1, 0, 0, 0],
        [0, 1, 1, 1, 1, 0, 0],
        [0, 0, 0, 0, 0, 1, 1],
        [0, 0, 0, 0, 0, 1, 1],
        [0, 0, 0, 0, 1, 1, 1],
        [0, 0, 0, 1, 1, 1, 0],
      ]);

      // The in-place version must agree with the returned permutations.
      const B = A.copy();
      B.doubly_lexical_ordering(true);

      // Sage's doctest asserts exactly these two properties of the reordered
      // matrix (matrix_mod2_dense.pyx:2274-2290):
      //
      //   for i in range(B.ncols()):
      //       for j in range(i):
      //           for k in reversed(range(B.nrows())):
      //               assert B[k][j] <= B[k][i]
      //               if B[k][j] < B[k][i]: break
      //
      // i.e. columns (read bottom to top) and rows (read right to left) are
      // lexicographically increasing.
      for (let i = 0; i < B.ncols; i++) {
        for (let j = 0; j < i; j++) {
          for (let k = B.nrows - 1; k >= 0; k--) {
            expect(B.get(k, j)).toBeLessThanOrEqual(B.get(k, i));
            if (B.get(k, j) < B.get(k, i)) break;
          }
        }
      }
      for (let i = 0; i < B.nrows; i++) {
        for (let j = 0; j < i; j++) {
          for (let k = B.ncols - 1; k >= 0; k--) {
            expect(B.get(j, k)).toBeLessThanOrEqual(B.get(i, k));
            if (B.get(j, k) < B.get(i, k)) break;
          }
        }
      }
    });
  });

  describe('pluq decomposition', () => {
    it('should compute PLUQ for a simple matrix', () => {
      const A = matrix_gf2_from_entries([
        [0, 1, 0, 1],
        [0, 1, 1, 1],
        [0, 0, 0, 1],
        [0, 1, 1, 0],
      ]);

      const [LU, P, Q] = pluq(A);

      // Sage doctest (matrix_mod2_dense.pyx:2732-2743):
      //   sage: LU, P, Q = pluq(A); LU
      //   [1 0 1 0]
      //   [1 1 0 0]
      //   [0 0 1 0]
      //   [1 1 1 0]
      //   sage: P
      //   [0, 1, 2, 3]
      //   sage: Q
      //   [1, 2, 3, 3]
      expect(rowsOf(LU)).toEqual([
        [1, 0, 1, 0],
        [1, 1, 0, 0],
        [0, 0, 1, 0],
        [1, 1, 1, 0],
      ]);
      expect(P).toEqual([0, 1, 2, 3]);
      expect(Q).toEqual([1, 2, 3, 3]);
    });

    it('should satisfy A = P*L*U*Q on the doctest matrix', () => {
      const A = matrix_gf2_from_entries([
        [0, 1, 0, 1],
        [0, 1, 1, 1],
        [0, 0, 0, 1],
        [0, 1, 1, 0],
      ]);

      const [LU, P, Q] = pluq(A);
      const r = A.rank();

      // L is unit lower triangular (m x r), U is upper triangular (r x n).
      const L: number[][] = [];
      for (let i = 0; i < A.nrows; i++) {
        L.push(new Array(r).fill(0));
        for (let j = 0; j < Math.min(i, r); j++) L[i]![j] = LU.get(i, j);
        if (i < r) L[i]![i] = 1;
      }
      const U: number[][] = [];
      for (let i = 0; i < r; i++) {
        U.push(new Array(A.ncols).fill(0));
        for (let j = i; j < A.ncols; j++) U[i]![j] = LU.get(i, j);
      }

      const prod: number[][] = [];
      for (let i = 0; i < A.nrows; i++) {
        prod.push(new Array(A.ncols).fill(0));
        for (let k = 0; k < r; k++) {
          if (L[i]![k] === 1) {
            for (let j = 0; j < A.ncols; j++) prod[i]![j] ^= U[k]![j]!;
          }
        }
      }

      // P and Q are transposition lists, applied in order; undo them.
      for (let i = Q.length - 1; i >= 0; i--) {
        const j = Q[i]!;
        if (j !== i) for (const row of prod) [row[i], row[j]] = [row[j]!, row[i]!];
      }
      for (let i = P.length - 1; i >= 0; i--) {
        const j = P[i]!;
        if (j !== i) [prod[i], prod[j]] = [prod[j]!, prod[i]!];
      }

      expect(prod).toEqual(rowsOf(A));
    });

    it('should compute PLUQ for identity matrix', () => {
      const I = identity_matrix_gf2(3);
      const [LU, P, Q] = pluq(I);

      // For identity, P and Q should effectively be identity permutations
      expect(P).toEqual([0, 1, 2]);
      expect(Q).toEqual([0, 1, 2]);

      // LU should be identity
      for (let i = 0; i < 3; i++) {
        for (let j = 0; j < 3; j++) {
          expect(LU.get(i, j)).toBe(i === j ? 1 : 0);
        }
      }
    });

    it('should handle zero matrix', () => {
      const Z = zero_matrix_gf2(3, 3);
      const [LU, P, Q] = pluq(Z);

      // LU should be all zeros
      for (let i = 0; i < 3; i++) {
        for (let j = 0; j < 3; j++) {
          expect(LU.get(i, j)).toBe(0);
        }
      }
    });

    it('should throw for unknown algorithm', () => {
      const A = identity_matrix_gf2(2);
      expect(() => pluq(A, 'unknown')).toThrow("Algorithm 'unknown' unknown.");
    });
  });

  describe('ple decomposition', () => {
    it('should compute PLE for a simple matrix', () => {
      const A = matrix_gf2_from_entries([
        [0, 1, 0, 1],
        [0, 1, 1, 1],
        [0, 0, 0, 1],
        [0, 1, 1, 0],
      ]);

      const [LU, P, Q] = ple(A);

      // Sage doctest (matrix_mod2_dense.pyx:2795-2808):
      //   sage: LU, P, Q = ple(A); LU
      //   [1 0 0 1]
      //   [1 1 0 0]
      //   [0 0 1 0]
      //   [1 1 1 0]
      //   sage: P
      //   [0, 1, 2, 3]
      //   sage: Q
      //   [1, 2, 3, 3]
      expect(rowsOf(LU)).toEqual([
        [1, 0, 0, 1],
        [1, 1, 0, 0],
        [0, 0, 1, 0],
        [1, 1, 1, 0],
      ]);
      expect(P).toEqual([0, 1, 2, 3]);
      expect(Q).toEqual([1, 2, 3, 3]);
    });

    it('should satisfy A = P*L*E on the doctest matrix', () => {
      const A = matrix_gf2_from_entries([
        [0, 1, 0, 1],
        [0, 1, 1, 1],
        [0, 0, 0, 1],
        [0, 1, 1, 0],
      ]);

      const [LU, P, Q] = ple(A);
      const r = A.rank();

      // L is unit lower triangular (m x r), compacted to the left; E is the
      // echelon form, with the pivot of row i in column Q[i].
      const L: number[][] = [];
      for (let i = 0; i < A.nrows; i++) {
        L.push(new Array(r).fill(0));
        for (let j = 0; j < Math.min(i, r); j++) L[i]![j] = LU.get(i, j);
        if (i < r) L[i]![i] = 1;
      }
      const E: number[][] = [];
      for (let i = 0; i < r; i++) {
        E.push(new Array(A.ncols).fill(0));
        E[i]![Q[i]!] = 1;
        for (let j = Q[i]! + 1; j < A.ncols; j++) E[i]![j] = LU.get(i, j);
      }

      const prod: number[][] = [];
      for (let i = 0; i < A.nrows; i++) {
        prod.push(new Array(A.ncols).fill(0));
        for (let k = 0; k < r; k++) {
          if (L[i]![k] === 1) {
            for (let j = 0; j < A.ncols; j++) prod[i]![j] ^= E[k]![j]!;
          }
        }
      }

      // P is a transposition list applied in order; undo it.
      for (let i = P.length - 1; i >= 0; i--) {
        const j = P[i]!;
        if (j !== i) [prod[i], prod[j]] = [prod[j]!, prod[i]!];
      }

      expect(prod).toEqual(rowsOf(A));
    });

    it('should compute PLE for identity matrix', () => {
      const I = identity_matrix_gf2(3);
      const [LU, P, Q] = ple(I);

      // For identity, P should be identity permutation
      expect(P).toEqual([0, 1, 2]);

      // LU should be identity
      for (let i = 0; i < 3; i++) {
        for (let j = 0; j < 3; j++) {
          expect(LU.get(i, j)).toBe(i === j ? 1 : 0);
        }
      }
    });

    it('should handle rectangular matrix', () => {
      const A = matrix_gf2_from_entries([
        [1, 0, 1],
        [0, 1, 1],
      ]);

      const [LU, P, Q] = ple(A);

      expect(LU.nrows).toBe(2);
      expect(LU.ncols).toBe(3);
      expect(P.length).toBe(2);
      expect(Q.length).toBe(3);
    });

    it('should throw for unknown algorithm', () => {
      const A = identity_matrix_gf2(2);
      expect(() => ple(A, 'unknown')).toThrow("Algorithm 'unknown' unknown.");
    });
  });

  describe('PNG data conversion', () => {
    it('should convert matrix to PNG data', () => {
      const A = matrix_gf2_from_entries([
        [1, 0],
        [0, 1],
      ]);

      const { width, height, pixels } = to_png_data(A);

      expect(width).toBe(2);
      expect(height).toBe(2);
      expect(pixels.length).toBe(4);

      // 1 -> black (0), 0 -> white (255)
      expect(pixels[0]).toBe(0); // (0,0) = 1 -> black
      expect(pixels[1]).toBe(255); // (0,1) = 0 -> white
      expect(pixels[2]).toBe(255); // (1,0) = 0 -> white
      expect(pixels[3]).toBe(0); // (1,1) = 1 -> black
    });

    it('should convert PNG data to matrix', () => {
      const pixels = new Uint8Array([0, 255, 255, 0]); // black, white, white, black
      const A = from_png_data(2, 2, pixels);

      expect(A.nrows).toBe(2);
      expect(A.ncols).toBe(2);
      expect(A.get(0, 0)).toBe(1); // black -> 1
      expect(A.get(0, 1)).toBe(0); // white -> 0
      expect(A.get(1, 0)).toBe(0); // white -> 0
      expect(A.get(1, 1)).toBe(1); // black -> 1
    });

    it('should round-trip matrix through PNG data', () => {
      const A = matrix_gf2_from_entries([
        [1, 0, 1],
        [0, 1, 0],
        [1, 1, 0],
      ]);

      const pngData = to_png_data(A);
      const B = from_png_data(pngData.width, pngData.height, pngData.pixels);

      expect(B.eq(A)).toBe(true);
    });

    it('should handle grayscale threshold', () => {
      // Values below 128 -> 1 (black), values >= 128 -> 0 (white)
      const pixels = new Uint8Array([0, 127, 128, 255]);
      const A = from_png_data(2, 2, pixels);

      expect(A.get(0, 0)).toBe(1); // 0 < 128
      expect(A.get(0, 1)).toBe(1); // 127 < 128
      expect(A.get(1, 0)).toBe(0); // 128 >= 128
      expect(A.get(1, 1)).toBe(0); // 255 >= 128
    });

    it('should throw for empty matrix in to_png_data', () => {
      const A = new Matrix_mod2_dense(0, 0);
      expect(() => to_png_data(A)).toThrow('cannot create image with dimensions 0 x 0');
    });

    it('should throw for mismatched pixel count', () => {
      const pixels = new Uint8Array([0, 255, 255]); // 3 pixels for 2x2 = 4
      expect(() => from_png_data(2, 2, pixels)).toThrow('Expected 4 pixels, got 3');
    });
  });

  describe('is_Gamma_free', () => {
    it('should detect Gamma-free matrix', () => {
      // No Gamma pattern: [1 1] / [0 0]
      const A = matrix_gf2_from_entries([
        [1, 1],
        [0, 0],
      ]);
      expect(A.is_Gamma_free()).toBe(true);
    });

    it('should detect Gamma pattern', () => {
      // Gamma pattern: [1 1] / [1 0]
      const A = matrix_gf2_from_entries([
        [1, 1],
        [1, 0],
      ]);
      expect(A.is_Gamma_free()).toBe(false);
    });

    it('should return certificate for Gamma pattern', () => {
      // Gamma pattern: [1 1] / [1 0]
      const A = matrix_gf2_from_entries([
        [1, 1],
        [1, 0],
      ]);
      const result = A.is_Gamma_free(true);
      expect(Array.isArray(result)).toBe(true);
      const [isFree, cert] = result as [boolean, [number, number, number, number] | null];
      expect(isFree).toBe(false);
      expect(cert).toEqual([0, 0, 1, 1]); // Matches SageMath: (0, 0, 1, 1)
    });

    it('should return null certificate when Gamma-free', () => {
      const A = matrix_gf2_from_entries([
        [1, 1],
        [0, 0],
      ]);
      const result = A.is_Gamma_free(true);
      expect(Array.isArray(result)).toBe(true);
      const [isFree, cert] = result as [boolean, null];
      expect(isFree).toBe(true);
      expect(cert).toBe(null);
    });

    it('should detect no Gamma pattern in larger matrix with matching corners', () => {
      // This matrix has [1][1] at corners but the third row has 1 at position 2
      // so no Gamma pattern because the bottom-right of any candidate is also 1
      const A = matrix_gf2_from_entries([
        [1, 0, 1],
        [0, 0, 0],
        [1, 0, 1],
      ]);
      expect(A.is_Gamma_free()).toBe(true);
    });

    it('should detect Gamma pattern in larger matrix', () => {
      // Example from SageMath tests: (False, (0, 0, 2, 2))
      // [1 0 1]   A[0,0]=1, A[0,2]=1
      // [0 0 0]
      // [1 0 0]   A[2,0]=1, A[2,2]=0 -> Gamma!
      const A = matrix_gf2_from_entries([
        [1, 0, 1],
        [0, 0, 0],
        [1, 0, 0],
      ]);
      const result = A.is_Gamma_free(true);
      expect(Array.isArray(result)).toBe(true);
      const [isFree, cert] = result as [boolean, [number, number, number, number]];
      expect(isFree).toBe(false);
      expect(cert).toEqual([0, 0, 2, 2]);
    });
  });

  describe('matrix arithmetic', () => {
    it('should add matrices (XOR)', () => {
      const A = matrix_gf2_from_entries([
        [1, 0],
        [0, 1],
      ]);
      const B = matrix_gf2_from_entries([
        [1, 1],
        [1, 1],
      ]);

      const C = A.add(B);

      expect(C.get(0, 0)).toBe(0); // 1 XOR 1
      expect(C.get(0, 1)).toBe(1); // 0 XOR 1
      expect(C.get(1, 0)).toBe(1); // 0 XOR 1
      expect(C.get(1, 1)).toBe(0); // 1 XOR 1
    });

    it('should multiply matrices', () => {
      const A = matrix_gf2_from_entries([
        [1, 1],
        [0, 1],
      ]);
      const B = matrix_gf2_from_entries([
        [1, 0],
        [1, 1],
      ]);

      const C = A.mul(B);

      // C[0,0] = (1*1) XOR (1*1) = 1 XOR 1 = 0
      // C[0,1] = (1*0) XOR (1*1) = 0 XOR 1 = 1
      // C[1,0] = (0*1) XOR (1*1) = 0 XOR 1 = 1
      // C[1,1] = (0*0) XOR (1*1) = 0 XOR 1 = 1
      expect(C.get(0, 0)).toBe(0);
      expect(C.get(0, 1)).toBe(1);
      expect(C.get(1, 0)).toBe(1);
      expect(C.get(1, 1)).toBe(1);
    });

    it('should compute determinant', () => {
      // Full rank matrix
      const A = matrix_gf2_from_entries([
        [1, 0],
        [0, 1],
      ]);
      expect(A.determinant()).toBe(1);

      // Singular matrix
      const B = matrix_gf2_from_entries([
        [1, 1],
        [1, 1],
      ]);
      expect(B.determinant()).toBe(0);
    });

    it('should compute rank', () => {
      const A = matrix_gf2_from_entries([
        [1, 0, 1],
        [0, 1, 1],
        [1, 1, 0],
      ]);
      expect(A.rank()).toBe(2); // Third row is sum of first two

      const I = identity_matrix_gf2(3);
      expect(I.rank()).toBe(3);
    });

    it('should compute inverse', () => {
      const A = matrix_gf2_from_entries([
        [1, 1],
        [0, 1],
      ]);

      const Ainv = A.inverse();
      const product = A.mul(Ainv);

      // A * A^-1 should be identity
      expect(product.eq(identity_matrix_gf2(2))).toBe(true);
    });
  });

  describe('echelon form', () => {
    it('should compute echelon form', () => {
      const A = matrix_gf2_from_entries([
        [0, 1, 1],
        [1, 0, 1],
        [1, 1, 0],
      ]);

      const E = A.echelon_form();

      // First non-zero entry in each row should be to the right of the one above
      let lastPivot = -1;
      for (let i = 0; i < E.nrows; i++) {
        for (let j = 0; j < E.ncols; j++) {
          if (E.get(i, j) === 1) {
            expect(j).toBeGreaterThan(lastPivot);
            lastPivot = j;
            break;
          }
        }
      }
    });

    it('should compute pivots', () => {
      const A = matrix_gf2_from_entries([
        [1, 0, 1],
        [0, 1, 1],
      ]);

      const pivots = A.pivots();
      expect(pivots).toEqual([0, 1]);
    });
  });

  describe('right kernel', () => {
    it('should compute right kernel', () => {
      const A = matrix_gf2_from_entries([
        [1, 0, 1],
        [0, 1, 1],
      ]);

      const K = A.right_kernel_matrix();

      // Kernel should have dimension 1 (ncols - rank = 3 - 2 = 1)
      expect(K.nrows).toBe(1);
      expect(K.ncols).toBe(3);

      // A * K^T should be zero
      const product = A.mul(K.transpose());
      for (let i = 0; i < product.nrows; i++) {
        for (let j = 0; j < product.ncols; j++) {
          expect(product.get(i, j)).toBe(0);
        }
      }
    });
  });
});
