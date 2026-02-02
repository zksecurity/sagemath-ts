/**
 * SageMath Lattice/Module Doctests
 *
 * Tests ported from SageMath documentation for lattice-related functionality.
 *
 * Reference files:
 * - sage/modules/free_module_integer.py
 * - sage/matrix/matrix_integer_dense.pyx
 * - sage/modules/misc.py
 */

import { describe, expect, it } from 'vitest';
import { IntegerMatrix, IntegerMatrixFromEntries } from '../matrix/index.js';
import { IntegerLattice, gramSchmidt, isLLLReduced, lllReduce } from './free_module_integer.js';

// ============================================================================
// LLL Reduction Tests from SageMath Documentation
// ============================================================================

describe('SageMath LLL Reduction Doctests', () => {
  /**
   * Basic LLL example with a full-rank matrix.
   *
   * Note: The original SageMath example uses matrix(ZZ, 3, 3, range(1,10)),
   * which is rank-deficient. Our implementation requires full-rank matrices.
   * Instead, we test with a similar but full-rank matrix.
   */
  describe('Basic LLL example with full-rank matrix', () => {
    it('should reduce a full-rank 3x3 matrix', () => {
      // A full-rank matrix similar in structure to the SageMath example
      const A = IntegerMatrixFromEntries([
        [1n, 2n, 3n],
        [4n, 5n, 6n],
        [7n, 8n, 10n], // Changed 9 to 10 to make it full-rank
      ]);

      const reduced = lllReduce(A);

      // Verify it's LLL-reduced
      expect(isLLLReduced(reduced)).toBe(true);

      // Verify the determinant is preserved
      const detOriginal = A.determinant().value;
      const detReduced = reduced.determinant().value;
      const absOriginal = detOriginal < 0n ? -detOriginal : detOriginal;
      const absReduced = detReduced < 0n ? -detReduced : detReduced;
      expect(absReduced).toBe(absOriginal);

      // Check that the first vector is short
      const rowNormsSq: bigint[] = [];
      for (let i = 0; i < reduced.nrows; i++) {
        let normSq = 0n;
        for (let j = 0; j < reduced.ncols; j++) {
          const v = reduced.get(i, j).value;
          normSq += v * v;
        }
        rowNormsSq.push(normSq);
      }

      // The first vector should be short relative to original
      expect(rowNormsSq[0]! <= rowNormsSq[rowNormsSq.length - 1]!).toBe(true);
    });
  });

  /**
   * Test with a challenging but full-rank matrix.
   *
   * Note: The original SageMath example uses a rank-deficient matrix.
   * Our implementation requires full-rank matrices, so we test with
   * a modified full-rank version.
   */
  describe('Challenging full-rank matrix LLL', () => {
    it('should handle matrix with varied magnitudes', () => {
      // Similar structure but made full-rank
      const M = IntegerMatrixFromEntries([
        [1n, 2n, 3n],
        [31n, 41n, 51n],
        [101n, 201n, 302n], // Changed 301 to 302 to make full-rank
      ]);

      const reduced = lllReduce(M);

      // Should be LLL-reduced
      expect(isLLLReduced(reduced)).toBe(true);

      // Determinant should be preserved
      const detOriginal = M.determinant().value;
      const detReduced = reduced.determinant().value;
      const absOriginal = detOriginal < 0n ? -detOriginal : detOriginal;
      const absReduced = detReduced < 0n ? -detReduced : detReduced;
      expect(absReduced).toBe(absOriginal);
    });
  });

  /**
   * From sage/modules/free_module_integer.py:
   *
   * sage: IntegerLattice([[1,0,3], [0,2,1], [0,2,7]])
   * Free module of degree 3 and rank 3 over Integer Ring
   * User basis matrix:
   * [-2  0  0]
   * [ 0  2  1]
   * [ 1 -2  2]
   */
  describe('IntegerLattice construction with LLL reduction', () => {
    it('should reduce [[1,0,3], [0,2,1], [0,2,7]]', () => {
      const L = IntegerLattice([
        [1, 0, 3],
        [0, 2, 1],
        [0, 2, 7],
      ]);

      // Should be LLL-reduced
      expect(isLLLReduced(L.reducedBasis)).toBe(true);

      // Original determinant: det([[1,0,3], [0,2,1], [0,2,7]]) = 1*(2*7 - 1*2) = 12
      const expectedDet = 12n;
      const actualDet = L.reducedBasis.determinant().value;
      const absDet = actualDet < 0n ? -actualDet : actualDet;
      expect(absDet).toBe(expectedDet);
    });
  });

  /**
   * From sage/modules/free_module_integer.py:
   *
   * sage: IntegerLattice([[1,0,-2], [0,2,5], [0,0,7]])
   * Free module of degree 3 and rank 3 over Integer Ring
   * User basis matrix:
   * [ 1  0 -2]
   * [ 1 -2  0]
   * [ 2  2  1]
   */
  describe('Another IntegerLattice example', () => {
    it('should reduce [[1,0,-2], [0,2,5], [0,0,7]]', () => {
      const L = IntegerLattice([
        [1, 0, -2],
        [0, 2, 5],
        [0, 0, 7],
      ]);

      expect(isLLLReduced(L.reducedBasis)).toBe(true);

      // Original is upper triangular, det = 1 * 2 * 7 = 14
      const expectedDet = 14n;
      const actualDet = L.reducedBasis.determinant().value;
      const absDet = actualDet < 0n ? -actualDet : actualDet;
      expect(absDet).toBe(expectedDet);
    });
  });

  /**
   * Simplified Extended GCD example (smaller scale for faster tests).
   * The full Magma handbook example with 10 large integers is too slow.
   */
  describe('Extended GCD via LLL (simplified)', () => {
    it('should compute extended GCD using LLL', () => {
      // Smaller example: GCD of [6, 10, 15] = 1
      const Q = [6n, 10n, 15n];
      const n = Q.length;
      const S = 10n;

      // Build the matrix X
      const entries: bigint[][] = [];
      for (let i = 0; i < n; i++) {
        const row: bigint[] = new Array(n + 1).fill(0n);
        row[0] = S * Q[i]!;
        row[i + 1] = 1n;
        entries.push(row);
      }

      const X = IntegerMatrixFromEntries(entries);
      const L = lllReduce(X);

      // Should be LLL-reduced
      expect(isLLLReduced(L)).toBe(true);

      // Find a row with small first element (the GCD relation)
      let foundSmall = false;
      for (let i = 0; i < L.nrows; i++) {
        const firstElem = L.get(i, 0).value;
        const absFirst = firstElem < 0n ? -firstElem : firstElem;
        if (absFirst <= S * 5n) {
          // Found a short vector
          foundSmall = true;
          break;
        }
      }
      expect(foundSmall).toBe(true);
    });
  });

  /**
   * Verify delta parameter validation
   *
   * sage: L = X.LLL(delta=2)
   * TypeError: delta must be <= 1
   */
  describe('LLL parameter validation', () => {
    it('should reject delta > 1', () => {
      const basis = IntegerMatrixFromEntries([
        [1n, 0n],
        [0n, 1n],
      ]);

      expect(() => lllReduce(basis, { delta: 2 })).toThrow();
    });

    it('should reject delta <= 0.25', () => {
      const basis = IntegerMatrixFromEntries([
        [1n, 0n],
        [0n, 1n],
      ]);

      expect(() => lllReduce(basis, { delta: 0.25 })).toThrow();
      expect(() => lllReduce(basis, { delta: 0.1 })).toThrow();
    });
  });
});

// ============================================================================
// Gram-Schmidt Orthogonalization Tests from SageMath Documentation
// ============================================================================

describe('SageMath Gram-Schmidt Doctests', () => {
  /**
   * From sage/modules/misc.py:
   *
   * sage: B = [vector([1,2,1/5]), vector([1,2,3]), vector([-1,0,0])]
   * sage: from sage.modules.misc import gram_schmidt
   * sage: G, mu = gram_schmidt(B)
   * sage: G
   * [(1, 2, 1/5), (-1/9, -2/9, 25/9), (-4/5, 2/5, 0)]
   * sage: G[0] * G[1]
   * 0
   * sage: G[0] * G[2]
   * 0
   * sage: G[1] * G[2]
   * 0
   * sage: mu
   * [      0       0       0]
   * [   10/9       0       0]
   * [-25/126    1/70       0]
   */
  describe('Gram-Schmidt from sage/modules/misc.py', () => {
    it('should orthogonalize [[1,2,0.2], [1,2,3], [-1,0,0]]', () => {
      const B = [
        [1, 2, 0.2], // 1/5 = 0.2
        [1, 2, 3],
        [-1, 0, 0],
      ];

      const { orthogonalBasis: G, mu } = gramSchmidt(B);

      // G[0] = [1, 2, 0.2] (unchanged)
      expect(G[0]![0]).toBeCloseTo(1, 10);
      expect(G[0]![1]).toBeCloseTo(2, 10);
      expect(G[0]![2]).toBeCloseTo(0.2, 10);

      // Verify orthogonality: G[i] . G[j] = 0 for i != j
      const dot01 = G[0]![0]! * G[1]![0]! + G[0]![1]! * G[1]![1]! + G[0]![2]! * G[1]![2]!;
      const dot02 = G[0]![0]! * G[2]![0]! + G[0]![1]! * G[2]![1]! + G[0]![2]! * G[2]![2]!;
      const dot12 = G[1]![0]! * G[2]![0]! + G[1]![1]! * G[2]![1]! + G[1]![2]! * G[2]![2]!;

      expect(dot01).toBeCloseTo(0, 10);
      expect(dot02).toBeCloseTo(0, 10);
      expect(dot12).toBeCloseTo(0, 10);

      // mu[1][0] should be 10/9
      expect(mu[1]![0]).toBeCloseTo(10 / 9, 8);
    });
  });

  /**
   * Identity-like basis should remain unchanged.
   */
  describe('Gram-Schmidt on orthogonal basis', () => {
    it('should preserve already orthogonal basis', () => {
      const B = [
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
      ];

      const { orthogonalBasis: G, mu, B: Bnorms } = gramSchmidt(B);

      // Orthogonal basis should be the same
      expect(G[0]).toEqual([1, 0, 0]);
      expect(G[1]).toEqual([0, 1, 0]);
      expect(G[2]).toEqual([0, 0, 1]);

      // mu should be identity-like (diagonal 1s, 0s below)
      expect(mu[0]![0]).toBe(1);
      expect(mu[1]![1]).toBe(1);
      expect(mu[2]![2]).toBe(1);
      expect(mu[1]![0]).toBeCloseTo(0, 10);
      expect(mu[2]![0]).toBeCloseTo(0, 10);
      expect(mu[2]![1]).toBeCloseTo(0, 10);

      // All norms should be 1
      expect(Bnorms[0]).toBeCloseTo(1, 10);
      expect(Bnorms[1]).toBeCloseTo(1, 10);
      expect(Bnorms[2]).toBeCloseTo(1, 10);
    });
  });

  /**
   * Gram-Schmidt with IntegerMatrix input.
   */
  describe('Gram-Schmidt with IntegerMatrix', () => {
    it('should handle IntegerMatrix input', () => {
      const basis = IntegerMatrixFromEntries([
        [4n, 1n, 3n],
        [2n, 5n, 6n],
        [1n, 2n, 7n],
      ]);

      const { orthogonalBasis: G } = gramSchmidt(basis);

      // Verify orthogonality
      const dot01 = G[0]![0]! * G[1]![0]! + G[0]![1]! * G[1]![1]! + G[0]![2]! * G[1]![2]!;
      const dot02 = G[0]![0]! * G[2]![0]! + G[0]![1]! * G[2]![1]! + G[0]![2]! * G[2]![2]!;
      const dot12 = G[1]![0]! * G[2]![0]! + G[1]![1]! * G[2]![1]! + G[1]![2]! * G[2]![2]!;

      expect(dot01).toBeCloseTo(0, 8);
      expect(dot02).toBeCloseTo(0, 8);
      expect(dot12).toBeCloseTo(0, 8);
    });
  });

  /**
   * From sage/matrix/matrix2.pyx:
   * The Gram-Schmidt decomposition satisfies A = M * G
   */
  describe('Gram-Schmidt decomposition property', () => {
    it('should satisfy B = mu * G* relationship', () => {
      const B = [
        [3, 1],
        [2, 2],
      ];

      const { orthogonalBasis: G, mu } = gramSchmidt(B);

      // Verify B = mu * G
      // B[0] = mu[0][0] * G[0] = 1 * G[0]
      // B[1] = mu[1][0] * G[0] + mu[1][1] * G[1]
      for (let i = 0; i < 2; i++) {
        const reconstructed = [0, 0];
        for (let j = 0; j <= i; j++) {
          reconstructed[0] += mu[i]![j]! * G[j]![0]!;
          reconstructed[1] += mu[i]![j]! * G[j]![1]!;
        }
        expect(reconstructed[0]).toBeCloseTo(B[i]![0]!, 10);
        expect(reconstructed[1]).toBeCloseTo(B[i]![1]!, 10);
      }
    });
  });
});

// ============================================================================
// Lattice Property Tests from SageMath Documentation
// ============================================================================

describe('SageMath Lattice Property Doctests', () => {
  /**
   * From sage/modules/free_module_integer.py:
   *
   * sage: L = sage.crypto.gen_lattice(m=10, seed=1337, lattice=True)
   * sage: L.volume()
   * 14641
   *
   * (Note: We can't replicate gen_lattice, but we can test volume preservation)
   */
  describe('Lattice determinant/volume preservation', () => {
    it('should preserve determinant through LLL', () => {
      const basis = IntegerMatrixFromEntries([
        [1n, 2n, 3n],
        [4n, 5n, 6n],
        [7n, 9n, 11n], // Not linearly dependent on first two rows
      ]);

      const detOriginal = basis.determinant().value;
      const reduced = lllReduce(basis);
      const detReduced = reduced.determinant().value;

      // Absolute values should match
      const absOriginal = detOriginal < 0n ? -detOriginal : detOriginal;
      const absReduced = detReduced < 0n ? -detReduced : detReduced;
      expect(absReduced).toBe(absOriginal);
    });

    it('should preserve determinant for larger lattice', () => {
      const basis = IntegerMatrixFromEntries([
        [10n, 2n, 3n, 1n],
        [4n, 15n, 6n, 2n],
        [7n, 8n, 20n, 3n],
        [1n, 2n, 3n, 25n],
      ]);

      const detOriginal = basis.determinant().value;
      const L = IntegerLattice(basis);

      // Should still be able to compute determinant
      const detReduced = L.reducedBasis.determinant().value;
      const absOriginal = detOriginal < 0n ? -detOriginal : detOriginal;
      const absReduced = detReduced < 0n ? -detReduced : detReduced;
      expect(absReduced).toBe(absOriginal);
    });
  });

  /**
   * Unimodular lattice test from sage/modules/free_module_integer.py:
   *
   * sage: from sage.modules.free_module_integer import IntegerLattice
   * sage: L = IntegerLattice([[1, 0], [0, 1]])
   * sage: L.is_unimodular()
   * True
   * sage: IntegerLattice([[2, 0], [0, 3]]).is_unimodular()
   * False
   */
  describe('Unimodular lattice detection', () => {
    it('should recognize ZZ^2 as unimodular', () => {
      const L = IntegerLattice([
        [1, 0],
        [0, 1],
      ]);

      // Volume = |det| = 1
      const det = L.reducedBasis.determinant().value;
      const absDet = det < 0n ? -det : det;
      expect(absDet).toBe(1n);
    });

    it('should recognize non-unimodular lattice', () => {
      const L = IntegerLattice([
        [2, 0],
        [0, 3],
      ]);

      // Volume = |det| = 6
      const det = L.reducedBasis.determinant().value;
      const absDet = det < 0n ? -det : det;
      expect(absDet).toBe(6n);
    });
  });

  /**
   * LLL reduction should produce vectors with decreasing lengths
   * according to the Gram-Schmidt norms.
   */
  describe('LLL reduction quality', () => {
    it('should produce short first vector', () => {
      const basis = IntegerMatrixFromEntries([
        [1n, 1000n],
        [0n, 1n],
      ]);

      const reduced = lllReduce(basis);

      // First vector should be short (approximately [1, 0] or similar)
      const norm0 = Number(reduced.get(0, 0).value) ** 2 + Number(reduced.get(0, 1).value) ** 2;

      // Original first vector has norm^2 = 1 + 1000000 = 1000001
      // LLL should find much shorter vector
      expect(norm0).toBeLessThan(100);
    });

    it('should satisfy Lovasz condition', () => {
      const basis = IntegerMatrixFromEntries([
        [15n, 23n, 11n],
        [46n, 79n, 29n],
        [81n, 25n, 67n],
      ]);

      const reduced = lllReduce(basis, { delta: 0.75 });
      const gs = gramSchmidt(reduced);

      // Check Lovasz condition: delta * B[k-1] <= B[k] + mu[k][k-1]^2 * B[k-1]
      const delta = 0.75;
      for (let k = 1; k < gs.B.length; k++) {
        const muKK1 = gs.mu[k]![k - 1]!;
        const lhs = delta * gs.B[k - 1]!;
        const rhs = gs.B[k]! + muKK1 * muKK1 * gs.B[k - 1]!;
        expect(lhs).toBeLessThanOrEqual(rhs + 1e-10);
      }
    });
  });
});

// ============================================================================
// isLLLReduced Tests
// ============================================================================

describe('SageMath isLLLReduced Doctests', () => {
  /**
   * From sage/matrix/matrix_integer_dense.pyx:
   *
   * sage: A = random_matrix(ZZ, 10, 10)
   * sage: L = A.LLL()
   * sage: A.is_LLL_reduced()
   * False
   * sage: L.is_LLL_reduced()
   * True
   */
  describe('isLLLReduced verification', () => {
    it('should return False for unreduced basis', () => {
      const A = IntegerMatrixFromEntries([
        [1n, 1000n],
        [0n, 1n],
      ]);

      // This basis has large mu[0][1] and likely fails Lovasz
      expect(isLLLReduced(A, 0.99, 0.501)).toBe(false);
    });

    it('should return True after LLL reduction', () => {
      const A = IntegerMatrixFromEntries([
        [1n, 1000n],
        [0n, 1n],
      ]);

      const L = lllReduce(A);
      expect(isLLLReduced(L)).toBe(true);
    });

    it('should respect delta parameter', () => {
      const basis = IntegerMatrixFromEntries([
        [5n, 3n, 7n],
        [11n, 13n, 17n],
        [19n, 23n, 29n],
      ]);

      const reduced099 = lllReduce(basis, { delta: 0.99 });
      const reduced075 = lllReduce(basis, { delta: 0.75 });

      // Both should be reduced with their respective parameters
      expect(isLLLReduced(reduced099, 0.99)).toBe(true);
      expect(isLLLReduced(reduced075, 0.75)).toBe(true);
    });
  });
});

// ============================================================================
// Edge Cases and Special Matrices
// ============================================================================

describe('Edge Cases from SageMath Documentation', () => {
  /**
   * From sage/matrix/matrix_integer_dense.pyx:
   *
   * sage: matrix(ZZ, 0, 0).LLL()
   * []
   * sage: matrix(ZZ, 3, 0).LLL()
   * []
   * sage: matrix(ZZ, 0, 3).LLL()
   * []
   */
  describe('Empty and degenerate matrices', () => {
    it('should handle 1x1 matrix', () => {
      const A = IntegerMatrixFromEntries([[5n]]);
      const reduced = lllReduce(A);

      expect(reduced.nrows).toBe(1);
      expect(reduced.ncols).toBe(1);
      expect(reduced.get(0, 0).value).toBe(5n);
    });

    it('should handle identity matrix', () => {
      const I = IntegerMatrixFromEntries([
        [1n, 0n, 0n],
        [0n, 1n, 0n],
        [0n, 0n, 1n],
      ]);

      const reduced = lllReduce(I);
      expect(isLLLReduced(reduced)).toBe(true);

      // Identity is already optimal
      const det = reduced.determinant().value;
      expect(det).toBe(1n);
    });

    it('should handle diagonal matrix', () => {
      const D = IntegerMatrixFromEntries([
        [3n, 0n],
        [0n, 7n],
      ]);

      const reduced = lllReduce(D);
      expect(isLLLReduced(reduced)).toBe(true);

      const det = reduced.determinant().value;
      const absDet = det < 0n ? -det : det;
      expect(absDet).toBe(21n);
    });
  });

  /**
   * Non-square matrix test (full row rank).
   *
   * Note: The original SageMath example used a rank-deficient matrix.
   * We test with a full-rank rectangular matrix instead.
   */
  describe('Non-square matrix (full row rank)', () => {
    it('should handle 3x4 full-rank matrix', () => {
      // A 3x4 matrix with full row rank (rank = 3)
      const M = IntegerMatrixFromEntries([
        [1n, 2n, 3n, 4n],
        [5n, 7n, 11n, 13n],
        [17n, 19n, 23n, 29n],
      ]);

      const reduced = lllReduce(M);

      // The result should have same dimensions
      expect(reduced.nrows).toBe(3);
      expect(reduced.ncols).toBe(4);

      // Should be LLL-reduced
      expect(isLLLReduced(reduced)).toBe(true);
    });
  });
});

// ============================================================================
// IntegerLattice Class Specific Tests
// ============================================================================

describe('IntegerLattice Class Tests from SageMath', () => {
  describe('Construction and basic properties', () => {
    it('should report correct rank and degree', () => {
      const L = IntegerLattice([
        [1, 0, 3],
        [0, 2, 1],
        [0, 2, 7],
      ]);

      expect(L.rank()).toBe(3);
      expect(L.degree()).toBe(3);
    });

    it('should support lllReduce: false option', () => {
      const basis = [
        [1, 1000],
        [0, 1],
      ];

      const Lreduced = IntegerLattice(basis, { lllReduce: true });
      const Lunreduced = IntegerLattice(basis, { lllReduce: false });

      // Unreduced should have original basis
      expect(Lunreduced.reducedBasis.get(0, 0).value).toBe(1n);
      expect(Lunreduced.reducedBasis.get(0, 1).value).toBe(1000n);

      // Reduced should have different (shorter) first vector
      const normReduced =
        Number(Lreduced.reducedBasis.get(0, 0).value) ** 2 +
        Number(Lreduced.reducedBasis.get(0, 1).value) ** 2;
      expect(normReduced).toBeLessThan(1000001); // Original norm^2
    });
  });

  describe('LLL method', () => {
    it('should improve basis when called explicitly', () => {
      const L = IntegerLattice(
        [
          [1, 1000],
          [0, 1],
        ],
        { lllReduce: false }
      );

      // Call LLL explicitly
      const reduced = L.LLL();

      expect(isLLLReduced(reduced)).toBe(true);

      // First vector should be short
      const norm0 = Number(reduced.get(0, 0).value) ** 2 + Number(reduced.get(0, 1).value) ** 2;
      expect(norm0).toBeLessThan(10);
    });

    it('should update internal reduced basis', () => {
      const L = IntegerLattice(
        [
          [1, 1000],
          [0, 1],
        ],
        { lllReduce: false }
      );

      const oldNorm =
        Number(L.reducedBasis.get(0, 0).value) ** 2 + Number(L.reducedBasis.get(0, 1).value) ** 2;

      L.LLL();

      const newNorm =
        Number(L.reducedBasis.get(0, 0).value) ** 2 + Number(L.reducedBasis.get(0, 1).value) ** 2;

      expect(newNorm).toBeLessThanOrEqual(oldNorm);
    });
  });
});

// ============================================================================
// Cryptographic Lattice Tests
// ============================================================================

describe('Cryptographic Lattice Examples', () => {
  /**
   * NTRU-style lattice test.
   * These lattices have a special structure used in lattice-based crypto.
   */
  describe('q-ary lattice structure', () => {
    it('should reduce NTRU-style lattice', () => {
      const q = 127n;
      const basis = IntegerMatrixFromEntries([
        [1n, 0n, 23n, 45n],
        [0n, 1n, 67n, 89n],
        [0n, 0n, q, 0n],
        [0n, 0n, 0n, q],
      ]);

      const reduced = lllReduce(basis);

      expect(isLLLReduced(reduced)).toBe(true);

      // Determinant should be preserved: det = 1 * 1 * 127 * 127 = 16129
      const detOriginal = basis.determinant().value;
      const detReduced = reduced.determinant().value;
      const absOriginal = detOriginal < 0n ? -detOriginal : detOriginal;
      const absReduced = detReduced < 0n ? -detReduced : detReduced;
      expect(absReduced).toBe(absOriginal);
      expect(absOriginal).toBe(16129n);
    });
  });

  /**
   * Lattice used for subset sum / knapsack problems.
   */
  describe('Knapsack lattice', () => {
    it('should reduce knapsack-style lattice', () => {
      // Simple knapsack lattice for finding linear combination
      const S = 100n;
      const weights = [23n, 37n, 41n, 53n];
      const n = weights.length;

      // Build lattice: [S*w_i, e_i]
      const entries: bigint[][] = [];
      for (let i = 0; i < n; i++) {
        const row: bigint[] = new Array(n + 1).fill(0n);
        row[0] = S * weights[i]!;
        row[i + 1] = 1n;
        entries.push(row);
      }

      const basis = IntegerMatrixFromEntries(entries);
      const reduced = lllReduce(basis);

      expect(isLLLReduced(reduced)).toBe(true);
    });
  });
});
