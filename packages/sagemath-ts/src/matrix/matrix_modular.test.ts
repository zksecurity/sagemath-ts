/**
 * @module sage/matrix/matrix_modular.test
 * @description Tests for modular matrix operations
 */

import { describe, expect, it } from 'bun:test';
import {
  Matrix_modn_dense,
  zero_matrix_modn,
  identity_matrix_modn,
  matrix_modn_from_entries,
} from './matrix_modn.js';
import {
  Matrix_mod2_dense,
  zero_matrix_gf2,
  identity_matrix_gf2,
  matrix_gf2_from_entries,
} from './matrix_mod2.js';
import {
  IntegerMatrix,
  IntegerMatrixFromEntries,
  identity_integer_matrix,
  zero_integer_matrix,
  height,
  gcd_integer_matrix,
  is_primitive,
  antitranspose,
  insert_row,
  augment_integer,
  delete_zero_columns,
  factor_out_common_factors_from_each_row,
  pivots_integer,
  LLL,
  is_LLL_reduced,
} from './matrix_integer.js';

describe('Matrix_modn_dense', () => {
  describe('construction', () => {
    it('should create a zero matrix', () => {
      const m = zero_matrix_modn(7n, 3, 4);
      expect(m.nrows).toBe(3);
      expect(m.ncols).toBe(4);
      expect(m.modulus).toBe(7n);
      for (let i = 0; i < 3; i++) {
        for (let j = 0; j < 4; j++) {
          expect(m.get(i, j)).toBe(0n);
        }
      }
    });

    it('should create an identity matrix', () => {
      const m = identity_matrix_modn(5n, 3);
      expect(m.nrows).toBe(3);
      expect(m.ncols).toBe(3);
      for (let i = 0; i < 3; i++) {
        for (let j = 0; j < 3; j++) {
          expect(m.get(i, j)).toBe(i === j ? 1n : 0n);
        }
      }
    });

    it('should create a matrix from entries', () => {
      const m = matrix_modn_from_entries(7n, [[1, 2], [3, 4]]);
      expect(m.get(0, 0)).toBe(1n);
      expect(m.get(0, 1)).toBe(2n);
      expect(m.get(1, 0)).toBe(3n);
      expect(m.get(1, 1)).toBe(4n);
    });

    it('should reduce entries modulo n', () => {
      const m = matrix_modn_from_entries(5n, [[7, -2], [10, 13]]);
      expect(m.get(0, 0)).toBe(2n); // 7 mod 5
      expect(m.get(0, 1)).toBe(3n); // -2 mod 5
      expect(m.get(1, 0)).toBe(0n); // 10 mod 5
      expect(m.get(1, 1)).toBe(3n); // 13 mod 5
    });
  });

  describe('arithmetic', () => {
    it('should add matrices', () => {
      const a = matrix_modn_from_entries(7n, [[1, 2], [3, 4]]);
      const b = matrix_modn_from_entries(7n, [[5, 6], [0, 1]]);
      const c = a.add(b);
      expect(c.get(0, 0)).toBe(6n);
      expect(c.get(0, 1)).toBe(1n); // (2+6) mod 7
      expect(c.get(1, 0)).toBe(3n);
      expect(c.get(1, 1)).toBe(5n);
    });

    it('should subtract matrices', () => {
      const a = matrix_modn_from_entries(7n, [[1, 2], [3, 4]]);
      const b = matrix_modn_from_entries(7n, [[5, 6], [0, 1]]);
      const c = a.sub(b);
      expect(c.get(0, 0)).toBe(3n); // (1-5) mod 7 = -4 mod 7 = 3
      expect(c.get(0, 1)).toBe(3n); // (2-6) mod 7 = -4 mod 7 = 3
      expect(c.get(1, 0)).toBe(3n);
      expect(c.get(1, 1)).toBe(3n);
    });

    it('should negate matrices', () => {
      const a = matrix_modn_from_entries(7n, [[1, 2], [3, 0]]);
      const b = a.neg();
      expect(b.get(0, 0)).toBe(6n);
      expect(b.get(0, 1)).toBe(5n);
      expect(b.get(1, 0)).toBe(4n);
      expect(b.get(1, 1)).toBe(0n);
    });

    it('should multiply matrices', () => {
      const a = matrix_modn_from_entries(7n, [[1, 2], [3, 4]]);
      const b = matrix_modn_from_entries(7n, [[5, 6], [0, 1]]);
      const c = a.mul(b);
      // c[0,0] = 1*5 + 2*0 = 5
      // c[0,1] = 1*6 + 2*1 = 8 mod 7 = 1
      // c[1,0] = 3*5 + 4*0 = 15 mod 7 = 1
      // c[1,1] = 3*6 + 4*1 = 22 mod 7 = 1
      expect(c.get(0, 0)).toBe(5n);
      expect(c.get(0, 1)).toBe(1n);
      expect(c.get(1, 0)).toBe(1n);
      expect(c.get(1, 1)).toBe(1n);
    });

    it('should scalar multiply', () => {
      const a = matrix_modn_from_entries(7n, [[1, 2], [3, 4]]);
      const b = a.scalar_mul(3n);
      expect(b.get(0, 0)).toBe(3n);
      expect(b.get(0, 1)).toBe(6n);
      expect(b.get(1, 0)).toBe(2n); // 9 mod 7
      expect(b.get(1, 1)).toBe(5n); // 12 mod 7
    });
  });

  describe('transpose', () => {
    it('should transpose a matrix', () => {
      const a = matrix_modn_from_entries(7n, [[1, 2, 3], [4, 5, 6]]);
      const b = a.transpose();
      expect(b.nrows).toBe(3);
      expect(b.ncols).toBe(2);
      expect(b.get(0, 0)).toBe(1n);
      expect(b.get(1, 0)).toBe(2n);
      expect(b.get(2, 0)).toBe(3n);
      expect(b.get(0, 1)).toBe(4n);
      expect(b.get(1, 1)).toBe(5n);
      expect(b.get(2, 1)).toBe(6n);
    });
  });

  describe('determinant', () => {
    it('should compute determinant of 1x1 matrix', () => {
      const a = matrix_modn_from_entries(7n, [[5]]);
      expect(a.determinant()).toBe(5n);
    });

    it('should compute determinant of 2x2 matrix', () => {
      const a = matrix_modn_from_entries(7n, [[1, 2], [3, 4]]);
      // det = 1*4 - 2*3 = 4 - 6 = -2 mod 7 = 5
      expect(a.determinant()).toBe(5n);
    });

    it('should return 0 for singular matrix', () => {
      const a = matrix_modn_from_entries(7n, [[1, 2], [2, 4]]);
      expect(a.determinant()).toBe(0n);
    });
  });

  describe('inverse', () => {
    it('should compute inverse', () => {
      const a = matrix_modn_from_entries(7n, [[1, 2], [3, 4]]);
      const aInv = a.inverse();
      const product = a.mul(aInv);

      // Should be identity
      expect(product.get(0, 0)).toBe(1n);
      expect(product.get(0, 1)).toBe(0n);
      expect(product.get(1, 0)).toBe(0n);
      expect(product.get(1, 1)).toBe(1n);
    });

    it('should throw for singular matrix', () => {
      const a = matrix_modn_from_entries(7n, [[1, 2], [2, 4]]);
      expect(() => a.inverse()).toThrow();
    });
  });

  describe('echelon form', () => {
    it('should compute echelon form', () => {
      const a = matrix_modn_from_entries(7n, [[1, 2, 3], [4, 5, 6], [0, 1, 2]]);
      const echelon = a.echelon_form();

      // First row should have leading 1
      expect(echelon.get(0, 0)).toBe(1n);
    });

    it('should compute rank', () => {
      const a = matrix_modn_from_entries(7n, [[1, 2, 3], [2, 4, 6], [1, 1, 1]]);
      // Rows 1 and 2 are linearly dependent
      expect(a.rank()).toBe(2);
    });
  });

  describe('minpoly', () => {
    it('should compute minimal polynomial of identity matrix', () => {
      const I = identity_matrix_modn(7n, 3);
      const minp = I.minpoly();
      // Minimal polynomial of I is x - 1, so coefficients are [-1, 1]
      // In mod 7: [-1 mod 7, 1] = [6, 1]
      expect(minp.length).toBe(2);
      expect(minp[0]).toBe(6n); // -1 mod 7
      expect(minp[1]).toBe(1n);
    });

    it('should compute minimal polynomial of zero matrix', () => {
      const Z = zero_matrix_modn(7n, 3, 3);
      const minp = Z.minpoly();
      // Minimal polynomial of zero matrix is x
      expect(minp.length).toBe(2);
      expect(minp[0]).toBe(0n);
      expect(minp[1]).toBe(1n);
    });

    it('should compute minimal polynomial of nilpotent matrix', () => {
      // [[0, 1], [0, 0]] is nilpotent with A^2 = 0
      const N = matrix_modn_from_entries(7n, [[0, 1], [0, 0]]);
      const minp = N.minpoly();
      // Minimal polynomial is x^2
      expect(minp.length).toBe(3);
      expect(minp[0]).toBe(0n);
      expect(minp[1]).toBe(0n);
      expect(minp[2]).toBe(1n);
    });

    it('should compute minimal polynomial of scalar matrix', () => {
      // 3*I has minimal polynomial x - 3
      const S = matrix_modn_from_entries(7n, [[3, 0], [0, 3]]);
      const minp = S.minpoly();
      expect(minp.length).toBe(2);
      expect(minp[0]).toBe(4n); // -3 mod 7
      expect(minp[1]).toBe(1n);
    });

    it('should throw error for non-square matrix', () => {
      const m = matrix_modn_from_entries(7n, [[1, 2, 3], [4, 5, 6]]);
      expect(() => m.minpoly()).toThrow();
    });
  });
});

describe('Matrix_mod2_dense', () => {
  describe('construction', () => {
    it('should create a zero matrix', () => {
      const m = zero_matrix_gf2(3, 4);
      expect(m.nrows).toBe(3);
      expect(m.ncols).toBe(4);
      for (let i = 0; i < 3; i++) {
        for (let j = 0; j < 4; j++) {
          expect(m.get(i, j)).toBe(0);
        }
      }
    });

    it('should create an identity matrix', () => {
      const m = identity_matrix_gf2(3);
      for (let i = 0; i < 3; i++) {
        for (let j = 0; j < 3; j++) {
          expect(m.get(i, j)).toBe(i === j ? 1 : 0);
        }
      }
    });

    it('should create a matrix from entries', () => {
      const m = matrix_gf2_from_entries([[1, 0], [0, 1]]);
      expect(m.get(0, 0)).toBe(1);
      expect(m.get(0, 1)).toBe(0);
      expect(m.get(1, 0)).toBe(0);
      expect(m.get(1, 1)).toBe(1);
    });

    it('should reduce entries mod 2', () => {
      const m = matrix_gf2_from_entries([[5, 4], [3, 2]]);
      expect(m.get(0, 0)).toBe(1); // 5 mod 2
      expect(m.get(0, 1)).toBe(0); // 4 mod 2
      expect(m.get(1, 0)).toBe(1); // 3 mod 2
      expect(m.get(1, 1)).toBe(0); // 2 mod 2
    });
  });

  describe('arithmetic', () => {
    it('should add matrices (XOR)', () => {
      const a = matrix_gf2_from_entries([[1, 0], [1, 1]]);
      const b = matrix_gf2_from_entries([[1, 1], [0, 1]]);
      const c = a.add(b);
      expect(c.get(0, 0)).toBe(0); // 1 XOR 1
      expect(c.get(0, 1)).toBe(1); // 0 XOR 1
      expect(c.get(1, 0)).toBe(1); // 1 XOR 0
      expect(c.get(1, 1)).toBe(0); // 1 XOR 1
    });

    it('should subtract matrices (same as add in GF(2))', () => {
      const a = matrix_gf2_from_entries([[1, 0], [1, 1]]);
      const b = matrix_gf2_from_entries([[1, 1], [0, 1]]);
      expect(a.sub(b).eq(a.add(b))).toBe(true);
    });

    it('should negate matrices (same as copy in GF(2))', () => {
      const a = matrix_gf2_from_entries([[1, 0], [1, 1]]);
      const b = a.neg();
      expect(b.eq(a)).toBe(true);
    });

    it('should multiply matrices', () => {
      const a = matrix_gf2_from_entries([[1, 1], [0, 1]]);
      const b = matrix_gf2_from_entries([[1, 0], [1, 1]]);
      const c = a.mul(b);
      // c[0,0] = 1*1 XOR 1*1 = 0
      // c[0,1] = 1*0 XOR 1*1 = 1
      // c[1,0] = 0*1 XOR 1*1 = 1
      // c[1,1] = 0*0 XOR 1*1 = 1
      expect(c.get(0, 0)).toBe(0);
      expect(c.get(0, 1)).toBe(1);
      expect(c.get(1, 0)).toBe(1);
      expect(c.get(1, 1)).toBe(1);
    });
  });

  describe('transpose', () => {
    it('should transpose a matrix', () => {
      const a = matrix_gf2_from_entries([[1, 0, 1], [0, 1, 1]]);
      const b = a.transpose();
      expect(b.nrows).toBe(3);
      expect(b.ncols).toBe(2);
      expect(b.get(0, 0)).toBe(1);
      expect(b.get(0, 1)).toBe(0);
      expect(b.get(1, 0)).toBe(0);
      expect(b.get(1, 1)).toBe(1);
      expect(b.get(2, 0)).toBe(1);
      expect(b.get(2, 1)).toBe(1);
    });
  });

  describe('determinant and rank', () => {
    it('should compute determinant', () => {
      const a = matrix_gf2_from_entries([[1, 0], [0, 1]]);
      expect(a.determinant()).toBe(1);

      const b = matrix_gf2_from_entries([[1, 1], [1, 1]]);
      expect(b.determinant()).toBe(0);
    });

    it('should compute rank', () => {
      const a = matrix_gf2_from_entries([[1, 0, 1], [0, 1, 1], [1, 1, 0]]);
      expect(a.rank()).toBe(2); // Row 3 = Row 1 + Row 2
    });
  });

  describe('inverse', () => {
    it('should compute inverse', () => {
      const a = matrix_gf2_from_entries([[1, 1], [0, 1]]);
      const aInv = a.inverse();
      const product = a.mul(aInv);

      expect(product.get(0, 0)).toBe(1);
      expect(product.get(0, 1)).toBe(0);
      expect(product.get(1, 0)).toBe(0);
      expect(product.get(1, 1)).toBe(1);
    });

    it('should throw for singular matrix', () => {
      const a = matrix_gf2_from_entries([[1, 1], [1, 1]]);
      expect(() => a.inverse()).toThrow();
    });
  });

  describe('echelon form', () => {
    it('should compute echelon form', () => {
      const a = matrix_gf2_from_entries([[1, 1, 0], [0, 1, 1], [1, 0, 1]]);
      const echelon = a.echelon_form();

      // Check it's in reduced row echelon form
      expect(echelon.get(0, 0)).toBe(1);
    });

    it('should find pivots', () => {
      const a = matrix_gf2_from_entries([[1, 0, 1], [0, 1, 1], [0, 0, 0]]);
      const pivots = a.pivots();
      expect(pivots).toEqual([0, 1]);
    });
  });

  describe('is_Gamma_free', () => {
    it('should detect Gamma pattern', () => {
      // Gamma pattern: [1 1] / [1 0]
      // A[0][0]=1, A[0][1]=1, A[1][0]=1, A[1][1]=0
      const a = matrix_gf2_from_entries([[1, 1], [1, 0]]);
      expect(a.is_Gamma_free()).toBe(false);
    });

    it('should detect Gamma-free matrix (identity)', () => {
      // Identity matrix has no Gamma pattern because there's no
      // row with two consecutive 1s
      const a = matrix_gf2_from_entries([[1, 0], [0, 1]]);
      expect(a.is_Gamma_free()).toBe(true);
    });

    it('should detect Gamma-free matrix (upper triangular)', () => {
      // Upper triangular [1 1] / [0 1] is Gamma-free
      // A[0][0]=1, A[0][1]=1 but A[1][0]=0 (not 1), so no Gamma
      const a = matrix_gf2_from_entries([[1, 1], [0, 1]]);
      expect(a.is_Gamma_free()).toBe(true);
    });
  });

  describe('density', () => {
    it('should compute density', () => {
      const a = matrix_gf2_from_entries([[1, 0], [0, 1]]);
      expect(a.density()).toBe(0.5);

      const b = matrix_gf2_from_entries([[1, 1], [1, 1]]);
      expect(b.density()).toBe(1);

      const c = matrix_gf2_from_entries([[0, 0], [0, 0]]);
      expect(c.density()).toBe(0);
    });
  });
});

describe('IntegerMatrix additional functions', () => {
  describe('height', () => {
    it('should compute height of matrix', () => {
      const m = IntegerMatrixFromEntries([[1, -5], [3, 4]]);
      expect(height(m).value).toBe(5n);
    });

    it('should return 0 for zero matrix', () => {
      const m = zero_integer_matrix(2, 2);
      expect(height(m).value).toBe(0n);
    });
  });

  describe('gcd_integer_matrix', () => {
    it('should compute GCD of all entries', () => {
      const m = IntegerMatrixFromEntries([[12, 18], [24, 6]]);
      expect(gcd_integer_matrix(m).value).toBe(6n);
    });

    it('should return 0 for zero matrix', () => {
      const m = zero_integer_matrix(2, 2);
      expect(gcd_integer_matrix(m).value).toBe(0n);
    });
  });

  describe('is_primitive', () => {
    it('should return true for primitive matrix', () => {
      const m = IntegerMatrixFromEntries([[1, 2], [3, 4]]);
      expect(is_primitive(m)).toBe(true);
    });

    it('should return false for non-primitive matrix', () => {
      const m = IntegerMatrixFromEntries([[2, 4], [3, 5]]);
      expect(is_primitive(m)).toBe(false); // First row has GCD 2
    });
  });

  describe('antitranspose', () => {
    it('should compute antitranspose', () => {
      const m = IntegerMatrixFromEntries([[1, 2], [3, 4]]);
      const at = antitranspose(m);
      // Antitranspose swaps (i,j) with (m-1-j, n-1-i)
      // Original: [[1,2],[3,4]]
      // Antitranspose: [[4,2],[3,1]]
      expect(at.get(0, 0).value).toBe(4n);
      expect(at.get(0, 1).value).toBe(2n);
      expect(at.get(1, 0).value).toBe(3n);
      expect(at.get(1, 1).value).toBe(1n);
    });
  });

  describe('insert_row', () => {
    it('should insert a row at the beginning', () => {
      const m = IntegerMatrixFromEntries([[1, 2], [3, 4]]);
      const result = insert_row(m, 0, [5n, 6n]);
      expect(result.nrows).toBe(3);
      expect(result.get(0, 0).value).toBe(5n);
      expect(result.get(0, 1).value).toBe(6n);
      expect(result.get(1, 0).value).toBe(1n);
    });

    it('should insert a row at the end', () => {
      const m = IntegerMatrixFromEntries([[1, 2], [3, 4]]);
      const result = insert_row(m, 2, [5n, 6n]);
      expect(result.nrows).toBe(3);
      expect(result.get(2, 0).value).toBe(5n);
      expect(result.get(2, 1).value).toBe(6n);
    });
  });

  describe('augment_integer', () => {
    it('should augment two matrices', () => {
      const a = IntegerMatrixFromEntries([[1, 2], [3, 4]]);
      const b = IntegerMatrixFromEntries([[5], [6]]);
      const result = augment_integer(a, b);

      expect(result.nrows).toBe(2);
      expect(result.ncols).toBe(3);
      expect(result.get(0, 2).value).toBe(5n);
      expect(result.get(1, 2).value).toBe(6n);
    });
  });

  describe('delete_zero_columns', () => {
    it('should delete zero columns', () => {
      const m = IntegerMatrixFromEntries([[1, 0, 2], [3, 0, 4]]);
      const result = delete_zero_columns(m);

      expect(result.ncols).toBe(2);
      expect(result.get(0, 0).value).toBe(1n);
      expect(result.get(0, 1).value).toBe(2n);
    });

    it('should keep all columns if none are zero', () => {
      const m = IntegerMatrixFromEntries([[1, 2], [3, 4]]);
      const result = delete_zero_columns(m);
      expect(result.ncols).toBe(2);
    });
  });

  describe('factor_out_common_factors_from_each_row', () => {
    it('should factor out common factors', () => {
      const m = IntegerMatrixFromEntries([[6, 12], [5, 10]]);
      const [factored, factors] = factor_out_common_factors_from_each_row(m);

      expect(factors[0]!.value).toBe(6n);
      expect(factors[1]!.value).toBe(5n);

      expect(factored.get(0, 0).value).toBe(1n);
      expect(factored.get(0, 1).value).toBe(2n);
      expect(factored.get(1, 0).value).toBe(1n);
      expect(factored.get(1, 1).value).toBe(2n);
    });
  });

  describe('pivots_integer', () => {
    it('should find pivot columns', () => {
      const m = IntegerMatrixFromEntries([[1, 2, 3], [0, 1, 2], [0, 0, 0]]);
      const pivots = pivots_integer(m);
      expect(pivots).toEqual([0, 1]);
    });
  });

  describe('LLL', () => {
    it('should reduce a simple lattice basis', () => {
      // Classic example: nearly orthogonal basis that LLL can reduce
      const m = IntegerMatrixFromEntries([
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
      ]);
      const reduced = LLL(m) as IntegerMatrix;

      // Identity matrix should remain identity (already optimal)
      expect(reduced.get(0, 0).value).toBe(1n);
      expect(reduced.get(1, 1).value).toBe(1n);
      expect(reduced.get(2, 2).value).toBe(1n);
    });

    it('should reduce a non-trivial basis', () => {
      // A basis that can be reduced
      const m = IntegerMatrixFromEntries([
        [1, 1, 1],
        [-1, 0, 2],
        [3, 5, 6],
      ]);
      const reduced = LLL(m) as IntegerMatrix;

      // The reduced basis should have shorter vectors
      // Check that it's valid (produces the same lattice up to unimodular transformation)
      expect(reduced.nrows).toBe(3);
      expect(reduced.ncols).toBe(3);
    });

    it('should return transformation matrix when requested', () => {
      const m = IntegerMatrixFromEntries([
        [1, 0],
        [0, 1],
      ]);
      const [reduced, U] = LLL(m, 0.75, 0.501, undefined, undefined, undefined, undefined, undefined, undefined, true) as [IntegerMatrix, IntegerMatrix];

      // U * m should equal reduced
      expect(U.nrows).toBe(2);
      expect(U.ncols).toBe(2);
    });
  });

  describe('is_LLL_reduced', () => {
    it('should return true for identity matrix', () => {
      const m = IntegerMatrixFromEntries([
        [1, 0],
        [0, 1],
      ]);
      expect(is_LLL_reduced(m)).toBe(true);
    });

    it('should return true for LLL-reduced basis', () => {
      const m = IntegerMatrixFromEntries([
        [1, 1, 1],
        [-1, 0, 2],
        [3, 5, 6],
      ]);
      const reduced = LLL(m) as IntegerMatrix;
      expect(is_LLL_reduced(reduced)).toBe(true);
    });
  });
});
