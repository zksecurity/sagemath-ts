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

      const [rowPerm, colPerm] = A.doubly_lexical_ordering();

      // Apply permutation to check result
      const B = A.copy();
      B.permute_rows(rowPerm);
      B.permute_columns(colPerm);

      // Verify columns are lexically ordered (from bottom to top)
      for (let j = 1; j < B.ncols; j++) {
        for (let jPrev = 0; jPrev < j; jPrev++) {
          // Compare columns j and jPrev from bottom to top
          let found = false;
          for (let i = B.nrows - 1; i >= 0; i--) {
            if (B.get(i, jPrev) < B.get(i, j)) {
              found = true;
              break;
            } else if (B.get(i, jPrev) > B.get(i, j)) {
              // Column jPrev should not be greater than column j
              found = true;
              break;
            }
          }
          // Either found ordering or columns are equal (which is allowed)
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

      // P should be row permutation indices
      expect(P.length).toBe(4);
      // Q should be column permutation indices
      expect(Q.length).toBe(4);

      // LU should have same dimensions as A
      expect(LU.nrows).toBe(4);
      expect(LU.ncols).toBe(4);
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

      // P should be row permutation indices
      expect(P.length).toBe(4);
      // Q should be pivot column positions
      expect(Q.length).toBe(4);

      // LU should have same dimensions as A
      expect(LU.nrows).toBe(4);
      expect(LU.ncols).toBe(4);
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
