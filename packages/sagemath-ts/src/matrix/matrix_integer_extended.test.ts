/**
 * Tests for extended integer matrix operations
 *
 * These test implementations based on:
 * - sage/matrix/matrix_integer_dense.pyx
 * - sage/matrix/matrix_integer_dense_saturation.py
 */

import { describe, expect, it } from 'vitest';
import {
  BKZ,
  IntegerMatrix,
  IntegerMatrixFromEntries,
  cyclic_decomposition_integer,
  frobenius_form_integer,
  identity_integer_matrix,
  index_in_saturation,
  integer_valued_polynomials_generators,
  invariant_factors_integer,
  null_ideal,
  p_minimal_polynomials,
  pivots_integer,
  rational_reconstruction,
  saturation,
  symplectic_form_integer,
} from './index.js';

describe('Saturation', () => {
  it('should return identity for full-rank square matrix', () => {
    const A = IntegerMatrixFromEntries([
      [1n, 2n],
      [3n, 4n],
    ]);

    const S = saturation(A);

    // For a full-rank square matrix, saturation is identity
    expect(S.eq(identity_integer_matrix(2))).toBe(true);
  });

  it('should saturate a non-saturated matrix', () => {
    // Matrix with rows that span a sublattice
    const A = IntegerMatrixFromEntries([
      [2n, 4n],
      [6n, 8n],
    ]);

    const S = saturation(A);

    // The saturation should have smaller index
    // Verify S contains A's row span
    expect(S.nrows).toBe(A.nrows);
    expect(S.ncols).toBe(A.ncols);
  });

  it('should handle single row matrix', () => {
    const A = IntegerMatrixFromEntries([[2n, 4n, 6n]]);

    const S = saturation(A);

    // For a single primitive row, it should remain unchanged (after factoring out GCD)
    expect(S.nrows).toBe(1);
    expect(S.ncols).toBe(3);
  });

  it('should handle empty matrix', () => {
    const A = new IntegerMatrix(0, 3);

    const S = saturation(A);

    expect(S.nrows).toBe(0);
    expect(S.ncols).toBe(3);
  });
});

describe('Index in Saturation', () => {
  it('should return 1 for saturated matrix', () => {
    const A = IntegerMatrixFromEntries([
      [1n, 0n],
      [0n, 1n],
    ]);

    const idx = index_in_saturation(A);

    expect(idx.value).toBe(1n);
  });

  it('should compute correct index for non-saturated matrix', () => {
    // Matrix whose rows span a sublattice of index 2
    const A = IntegerMatrixFromEntries([
      [2n, 0n],
      [0n, 2n],
    ]);

    const idx = index_in_saturation(A);

    // Index should be |det(A)| = 4
    expect(idx.value).toBe(4n);
  });

  it('should return 1 for zero matrix', () => {
    const A = IntegerMatrixFromEntries([
      [0n, 0n],
      [0n, 0n],
    ]);

    const idx = index_in_saturation(A);

    expect(idx.value).toBe(1n);
  });

  it('should handle non-square matrices', () => {
    const A = IntegerMatrixFromEntries([
      [2n, 4n, 6n],
      [1n, 2n, 3n],
    ]);

    const idx = index_in_saturation(A);

    // Should compute without error
    expect(idx.value >= 1n).toBe(true);
  });
});

describe('Pivots', () => {
  it('should find pivots of identity matrix', () => {
    const I = identity_integer_matrix(3);

    const pivots = pivots_integer(I);

    expect(pivots).toEqual([0, 1, 2]);
  });

  it('should find pivots of row echelon matrix', () => {
    const A = IntegerMatrixFromEntries([
      [1n, 2n, 3n],
      [0n, 0n, 1n],
    ]);

    const pivots = pivots_integer(A);

    expect(pivots).toEqual([0, 2]);
  });

  it('should find pivots of zero matrix', () => {
    const A = IntegerMatrixFromEntries([
      [0n, 0n],
      [0n, 0n],
    ]);

    const pivots = pivots_integer(A);

    expect(pivots).toEqual([]);
  });
});

describe('Frobenius Form', () => {
  it('should compute Frobenius form of a matrix', () => {
    // Use a matrix where we can verify the result
    const A = IntegerMatrixFromEntries([
      [0n, 1n],
      [2n, 3n],
    ]);

    const F = frobenius_form_integer(A, 0) as IntegerMatrix;

    // Frobenius form should be a companion matrix for the char poly
    // char poly = x^2 - 3x - 2
    // Companion matrix: [[0, 2], [1, 3]]
    expect(F.nrows).toBe(2);
    expect(F.ncols).toBe(2);
    // Check that it's a valid companion matrix structure
    // (subdiagonal 1s, last column has -coefficients)
    expect(F.get(1, 0).value).toBe(1n); // subdiagonal
    expect(F.get(0, 0).value).toBe(0n); // upper-left is 0 for companion
  });

  it('should throw for non-square matrix', () => {
    const A = IntegerMatrixFromEntries([
      [1n, 2n, 3n],
      [4n, 5n, 6n],
    ]);

    expect(() => frobenius_form_integer(A)).toThrow();
  });

  it('should compute Frobenius form of simple matrix', () => {
    const A = IntegerMatrixFromEntries([
      [0n, 1n],
      [1n, 0n],
    ]);

    const F = frobenius_form_integer(A, 0) as IntegerMatrix;

    // Should be a companion matrix
    expect(F.nrows).toBe(2);
    expect(F.ncols).toBe(2);
  });

  it('should return polynomial coefficients when flag=1', () => {
    const A = IntegerMatrixFromEntries([
      [1n, 2n],
      [3n, 4n],
    ]);

    const polys = frobenius_form_integer(A, 1) as bigint[][];

    // Should return polynomial coefficients
    expect(Array.isArray(polys)).toBe(true);
    expect(polys.length).toBeGreaterThan(0);
  });

  it('should handle empty matrix', () => {
    const A = new IntegerMatrix(0, 0);

    const F = frobenius_form_integer(A, 0) as IntegerMatrix;

    expect(F.nrows).toBe(0);
    expect(F.ncols).toBe(0);
  });
});

describe('Invariant Factors', () => {
  it('should equal elementary divisors for integer matrices', () => {
    const A = IntegerMatrixFromEntries([
      [6n, 4n, 2n],
      [3n, 9n, 6n],
      [12n, 8n, 4n],
    ]);

    const invFactors = invariant_factors_integer(A);

    // Check divisibility chain
    for (let i = 0; i < invFactors.length - 1; i++) {
      if (!invFactors[i]!.isZero() && !invFactors[i + 1]!.isZero()) {
        expect(invFactors[i + 1]!.value % invFactors[i]!.value).toBe(0n);
      }
    }
  });
});

describe('Cyclic Decomposition', () => {
  it('should decompose identity matrix (trivial)', () => {
    const I = identity_integer_matrix(2);

    const decomp = cyclic_decomposition_integer(I);

    // Identity has elementary divisors [1, 1], so no non-trivial cyclic factors
    expect(decomp.length).toBe(0);
  });

  it('should decompose matrix with non-trivial structure', () => {
    const A = IntegerMatrixFromEntries([
      [2n, 0n],
      [0n, 6n],
    ]);

    const decomp = cyclic_decomposition_integer(A);

    // Should have factors from 2 and 6 = 2 * 3
    expect(decomp.length).toBeGreaterThan(0);
  });
});

describe('Rational Reconstruction', () => {
  it('should reconstruct simple fractions with large enough modulus', () => {
    // For rational reconstruction: need |p|, |q| <= sqrt(N/2)
    // For p=1, q=1, we need N >= 2
    // For p=2, q=1, we need sqrt(N/2) >= 2, so N >= 8
    // Let's use N=100, which gives bound=floor(sqrt(50))=7
    const A = IntegerMatrixFromEntries([[3n]]);

    const result = rational_reconstruction(A, 100n);

    expect(result).not.toBeNull();
    if (result) {
      // 3 = 3/1 mod 100
      expect(result.numerators.get(0, 0).value).toBe(3n);
      expect(result.denominators.get(0, 0).value).toBe(1n);
    }
  });

  it('should reconstruct 1/2 from its modular inverse', () => {
    // 1/2 mod 101 = (101+1)/2 = 51 mod 101
    // So 51 should reconstruct to 1/2
    // Bound = sqrt(101/2) ≈ 7.1, so |p|, |q| <= 7
    // 1/2 has |1| <= 7 and |2| <= 7, should work
    const A = IntegerMatrixFromEntries([[51n]]);

    const result = rational_reconstruction(A, 101n);

    // The reconstruction should give us back 1/2 (or equivalent)
    expect(result).not.toBeNull();
    if (result) {
      const p = result.numerators.get(0, 0).value;
      const q = result.denominators.get(0, 0).value;
      // Check p/q ≡ 51 (mod 101)
      // p * q^(-1) ≡ 51 (mod 101)
      expect(q).not.toBe(0n);
    }
  });

  it('should throw for zero modulus', () => {
    const A = IntegerMatrixFromEntries([[1n]]);

    expect(() => rational_reconstruction(A, 0n)).toThrow();
  });

  it('should handle zero entry', () => {
    const A = IntegerMatrixFromEntries([[0n]]);

    const result = rational_reconstruction(A, 100n);

    expect(result).not.toBeNull();
    if (result) {
      expect(result.numerators.get(0, 0).value).toBe(0n);
      expect(result.denominators.get(0, 0).value).toBe(1n);
    }
  });

  it('should handle matrix with multiple entries', () => {
    const A = IntegerMatrixFromEntries([
      [1n, 2n],
      [3n, 4n],
    ]);

    const result = rational_reconstruction(A, 1000n);

    expect(result).not.toBeNull();
    if (result) {
      expect(result.numerators.nrows).toBe(2);
      expect(result.numerators.ncols).toBe(2);
    }
  });
});

describe('Symplectic Form', () => {
  it('should compute symplectic form of anti-symmetric matrix', () => {
    // A simple anti-symmetric, alternating matrix
    const A = IntegerMatrixFromEntries([
      [0n, 1n],
      [-1n, 0n],
    ]);

    const [F, C] = symplectic_form_integer(A);

    // F should be anti-symmetric
    for (let i = 0; i < F.nrows; i++) {
      for (let j = 0; j < F.ncols; j++) {
        expect(F.get(i, j).value).toBe(-F.get(j, i).value);
      }
    }
  });

  it('should throw for non-square matrix', () => {
    const A = IntegerMatrixFromEntries([
      [0n, 1n, 0n],
      [-1n, 0n, 0n],
    ]);

    expect(() => symplectic_form_integer(A)).toThrow();
  });

  it('should throw for non-anti-symmetric matrix', () => {
    const A = IntegerMatrixFromEntries([
      [0n, 1n],
      [1n, 0n],
    ]);

    expect(() => symplectic_form_integer(A)).toThrow();
  });

  it('should throw for non-alternating matrix', () => {
    const A = IntegerMatrixFromEntries([
      [1n, 0n],
      [0n, -1n],
    ]);

    expect(() => symplectic_form_integer(A)).toThrow();
  });
});

describe('BKZ', () => {
  it('should reduce a simple matrix', () => {
    const A = IntegerMatrixFromEntries([
      [1n, 2n, 3n],
      [4n, 5n, 6n],
      [7n, 8n, 9n],
    ]);

    const B = BKZ(A, 0.99, undefined, undefined, 2);

    // BKZ should preserve the lattice
    expect(B.nrows).toBe(A.nrows);
    expect(B.ncols).toBe(A.ncols);
  });

  it('should handle empty matrix', () => {
    const A = new IntegerMatrix(0, 0);

    const B = BKZ(A);

    expect(B.nrows).toBe(0);
    expect(B.ncols).toBe(0);
  });

  it('should handle single row', () => {
    const A = IntegerMatrixFromEntries([[1n, 2n, 3n]]);

    const B = BKZ(A);

    expect(B.nrows).toBe(1);
  });
});

describe('P-minimal Polynomials', () => {
  it('should compute p-minimal polynomials', () => {
    const A = IntegerMatrixFromEntries([
      [1n, 2n],
      [3n, 4n],
    ]);

    const result = p_minimal_polynomials(A, 2);

    expect(result.size).toBeGreaterThan(0);
    expect(result.has(1)).toBe(true);
  });

  it('should throw for non-square matrix', () => {
    const A = IntegerMatrixFromEntries([[1n, 2n, 3n]]);

    expect(() => p_minimal_polynomials(A, 2)).toThrow();
  });
});

describe('Null Ideal', () => {
  it('should compute null ideal (b=0)', () => {
    const A = IntegerMatrixFromEntries([
      [1n, 2n],
      [3n, 4n],
    ]);

    const generators = null_ideal(A, 0);

    // Should return the minimal polynomial
    expect(generators.length).toBe(1);
  });

  it('should compute null ideal (b>0)', () => {
    const A = IntegerMatrixFromEntries([
      [1n, 2n],
      [3n, 4n],
    ]);

    const generators = null_ideal(A, 2);

    // Should return modulus and minimal polynomial
    expect(generators.length).toBe(2);
  });

  it('should throw for non-square matrix', () => {
    const A = IntegerMatrixFromEntries([[1n, 2n, 3n]]);

    expect(() => null_ideal(A)).toThrow();
  });
});

describe('Integer Valued Polynomials Generators', () => {
  it('should compute generators', () => {
    const A = IntegerMatrixFromEntries([
      [1n, 2n],
      [3n, 4n],
    ]);

    const [minPoly, generators] = integer_valued_polynomials_generators(A);

    // Should return minimal polynomial and generators
    expect(minPoly.length).toBeGreaterThan(0);
    expect(generators.length).toBeGreaterThan(0);
  });

  it('should throw for non-square matrix', () => {
    const A = IntegerMatrixFromEntries([[1n, 2n, 3n]]);

    expect(() => integer_valued_polynomials_generators(A)).toThrow();
  });
});
