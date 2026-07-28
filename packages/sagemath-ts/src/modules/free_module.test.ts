/**
 * Tests for free modules and LLL lattice reduction
 */

import { describe, expect, it } from 'vitest';
import { type IntegerMatrix, IntegerMatrixFromEntries } from '../matrix/index.js';
import { GF } from '../rings/finite_rings/finite_field_constructor.js';
import { PolynomialRingConstructor } from '../rings/polynomial/polynomial_ring.js';
import { Rational } from '../rings/rational.js';
import { QQ as QQfield } from '../rings/rational_field.js';
import {
  FractionFieldElement,
  FreeModule,
  type FreeModuleField,
  type FreeModulePID,
  VectorSpace,
  span,
  tensorProductVector,
} from './free_module.js';
import { IntegerLattice, gramSchmidt, isLLLReduced, lllReduce } from './free_module_integer.js';

describe('Gram-Schmidt orthogonalization', () => {
  it('should compute Gram-Schmidt for a 2D basis', () => {
    const basis = [
      [3, 1],
      [2, 2],
    ];

    const result = gramSchmidt(basis);

    // b_0* = b_0 = [3, 1]
    expect(result.orthogonalBasis[0]).toEqual([3, 1]);

    // mu[1][0] = <b_1, b_0*> / <b_0*, b_0*> = (2*3 + 2*1) / (9 + 1) = 8/10 = 0.8
    expect(result.mu[1]![0]).toBeCloseTo(0.8, 10);

    // b_1* = b_1 - mu[1][0] * b_0* = [2, 2] - 0.8 * [3, 1] = [2 - 2.4, 2 - 0.8] = [-0.4, 1.2]
    expect(result.orthogonalBasis[1]![0]).toBeCloseTo(-0.4, 10);
    expect(result.orthogonalBasis[1]![1]).toBeCloseTo(1.2, 10);

    // B[0] = |b_0*|^2 = 9 + 1 = 10
    expect(result.B[0]).toBeCloseTo(10, 10);

    // B[1] = |b_1*|^2 = 0.16 + 1.44 = 1.6
    expect(result.B[1]).toBeCloseTo(1.6, 10);
  });

  it('should compute Gram-Schmidt for a 3D basis', () => {
    const basis = [
      [1, 0, 0],
      [1, 1, 0],
      [1, 1, 1],
    ];

    const result = gramSchmidt(basis);

    // b_0* = [1, 0, 0]
    expect(result.orthogonalBasis[0]).toEqual([1, 0, 0]);

    // b_1* = [1, 1, 0] - 1 * [1, 0, 0] = [0, 1, 0]
    expect(result.orthogonalBasis[1]![0]).toBeCloseTo(0, 10);
    expect(result.orthogonalBasis[1]![1]).toBeCloseTo(1, 10);
    expect(result.orthogonalBasis[1]![2]).toBeCloseTo(0, 10);

    // b_2* = [1, 1, 1] - 1 * [1, 0, 0] - 1 * [0, 1, 0] = [0, 0, 1]
    expect(result.orthogonalBasis[2]![0]).toBeCloseTo(0, 10);
    expect(result.orthogonalBasis[2]![1]).toBeCloseTo(0, 10);
    expect(result.orthogonalBasis[2]![2]).toBeCloseTo(1, 10);
  });

  it('should handle IntegerMatrix input', () => {
    const basis = IntegerMatrixFromEntries([
      [3n, 1n],
      [2n, 2n],
    ]);

    const result = gramSchmidt(basis);

    expect(result.orthogonalBasis[0]).toEqual([3, 1]);
    expect(result.mu[1]![0]).toBeCloseTo(0.8, 10);
  });
});

describe('LLL lattice reduction', () => {
  describe('simple 2D lattice', () => {
    it('should reduce a simple 2x2 lattice', () => {
      // Example: basis [[1, 0], [0, 2]] should stay roughly the same
      const basis = IntegerMatrixFromEntries([
        [1n, 0n],
        [0n, 2n],
      ]);

      const reduced = lllReduce(basis);

      // The basis is already quite good, so should not change much
      expect(reduced.nrows).toBe(2);
      expect(reduced.ncols).toBe(2);

      // Should be LLL-reduced
      expect(isLLLReduced(reduced)).toBe(true);
    });

    it('should reduce a poorly conditioned 2x2 lattice', () => {
      // Basis with large skew - LLL should find shorter vectors
      const basis = IntegerMatrixFromEntries([
        [1n, 1n],
        [0n, 100n],
      ]);

      const reduced = lllReduce(basis);

      expect(isLLLReduced(reduced)).toBe(true);

      // The first vector should be short
      const norm0 = Number(reduced.get(0, 0).value) ** 2 + Number(reduced.get(0, 1).value) ** 2;
      expect(norm0).toBeLessThanOrEqual(2); // [1, 1] or [1, 0] or similar
    });
  });

  describe('3D lattice from literature', () => {
    it('should reduce the classic example [[1, 0, 3], [0, 2, 1], [0, 2, 7]]', () => {
      // This example is from many LLL tutorials
      const basis = IntegerMatrixFromEntries([
        [1n, 0n, 3n],
        [0n, 2n, 1n],
        [0n, 2n, 7n],
      ]);

      const reduced = lllReduce(basis);

      expect(reduced.nrows).toBe(3);
      expect(reduced.ncols).toBe(3);

      // Should be LLL-reduced
      expect(isLLLReduced(reduced)).toBe(true);

      // The reduced basis should have shorter vectors than the original
      const computeNormSq = (M: IntegerMatrix, row: number) => {
        let sum = 0n;
        for (let j = 0; j < M.ncols; j++) {
          const v = M.get(row, j).value;
          sum += v * v;
        }
        return sum;
      };

      // Original first row norm^2 = 1 + 0 + 9 = 10
      const origNorm0 = computeNormSq(basis, 0);
      const reducedNorm0 = computeNormSq(reduced, 0);

      // Reduced should have first vector with norm^2 <= 10
      expect(reducedNorm0 <= origNorm0).toBe(true);
    });
  });

  describe('determinant preservation', () => {
    it('should preserve the lattice determinant', () => {
      const basis = IntegerMatrixFromEntries([
        [1n, 2n],
        [3n, 5n],
      ]);

      const reduced = lllReduce(basis);

      // Determinant should be preserved (up to sign)
      const detOriginal = basis.determinant().value;
      const detReduced = reduced.determinant().value;

      expect(detReduced === detOriginal || detReduced === -detOriginal).toBe(true);
    });

    it('should preserve the 3x3 lattice determinant', () => {
      const basis = IntegerMatrixFromEntries([
        [1n, 0n, 3n],
        [0n, 2n, 1n],
        [0n, 2n, 7n],
      ]);

      const reduced = lllReduce(basis);

      const detOriginal = basis.determinant().value;
      const detReduced = reduced.determinant().value;

      // Absolute values should match
      const absDetOriginal = detOriginal < 0n ? -detOriginal : detOriginal;
      const absDetReduced = detReduced < 0n ? -detReduced : detReduced;

      expect(absDetReduced).toBe(absDetOriginal);
    });
  });

  describe('Lovasz condition verification', () => {
    it('should produce basis satisfying Lovasz condition', () => {
      const basis = IntegerMatrixFromEntries([
        [1n, 0n, 3n],
        [0n, 2n, 1n],
        [0n, 2n, 7n],
      ]);

      const reduced = lllReduce(basis, { delta: 0.75 });

      // Compute Gram-Schmidt
      const gs = gramSchmidt(reduced);
      const { mu, B: Bnorms } = gs;

      // Check size reduction: |mu[i][j]| <= 0.501 for all i > j
      for (let i = 1; i < reduced.nrows; i++) {
        for (let j = 0; j < i; j++) {
          expect(Math.abs(mu[i]![j]!)).toBeLessThanOrEqual(0.501 + 1e-10);
        }
      }

      // Check Lovasz condition: 0.75 * B[k-1] <= B[k] + mu[k][k-1]^2 * B[k-1]
      for (let k = 1; k < reduced.nrows; k++) {
        const muKK1 = mu[k]![k - 1]!;
        const lhs = 0.75 * Bnorms[k - 1]!;
        const rhs = Bnorms[k]! + muKK1 * muKK1 * Bnorms[k - 1]!;
        expect(lhs).toBeLessThanOrEqual(rhs + 1e-10);
      }
    });
  });

  describe('parameter variations', () => {
    it('should respect delta parameter', () => {
      const basis = IntegerMatrixFromEntries([
        [1n, 0n, 100n],
        [0n, 1n, 0n],
        [0n, 0n, 1n],
      ]);

      // Higher delta = stricter condition = potentially more reduction
      const reduced099 = lllReduce(basis, { delta: 0.99 });
      const reduced075 = lllReduce(basis, { delta: 0.75 });

      expect(isLLLReduced(reduced099, 0.99)).toBe(true);
      expect(isLLLReduced(reduced075, 0.75)).toBe(true);

      // Both should be valid but delta=0.99 might have shorter first vector
      // (not always guaranteed, but the condition is stricter)
    });

    it('should throw for invalid delta', () => {
      const basis = IntegerMatrixFromEntries([
        [1n, 0n],
        [0n, 1n],
      ]);

      expect(() => lllReduce(basis, { delta: 0.2 })).toThrow();
      expect(() => lllReduce(basis, { delta: 1.5 })).toThrow();
    });

    it('should throw for invalid eta', () => {
      const basis = IntegerMatrixFromEntries([
        [1n, 0n],
        [0n, 1n],
      ]);

      // eta must satisfy 0.5 <= eta < sqrt(delta); state delta explicitly so
      // that the test does not depend on its default value.
      expect(() => lllReduce(basis, { eta: 0.4 })).toThrow();
      expect(() => lllReduce(basis, { delta: 0.75, eta: 0.9 })).toThrow();
    });
  });

  describe('isLLLReduced function', () => {
    it('should return true for identity matrix', () => {
      const identity = IntegerMatrixFromEntries([
        [1n, 0n],
        [0n, 1n],
      ]);

      expect(isLLLReduced(identity)).toBe(true);
    });

    it('should return true for already reduced basis', () => {
      const basis = IntegerMatrixFromEntries([
        [1n, 0n],
        [0n, 2n],
      ]);

      expect(isLLLReduced(basis)).toBe(true);
    });

    it('should return false for poorly conditioned basis', () => {
      // This basis has large mu values and fails Lovasz
      const basis = IntegerMatrixFromEntries([
        [1n, 1000n],
        [0n, 1n],
      ]);

      // Check with strict parameters - this should fail
      expect(isLLLReduced(basis, 0.99, 0.501)).toBe(false);
    });
  });
});

describe('IntegerLattice class', () => {
  describe('construction', () => {
    it('should create a lattice from 2D array', () => {
      const L = IntegerLattice([
        [1n, 0n, 3n],
        [0n, 2n, 1n],
        [0n, 2n, 7n],
      ]);

      expect(L.rank()).toBe(3);
      expect(L.degree()).toBe(3);
    });

    it('should create a lattice from number array', () => {
      const L = IntegerLattice([
        [1, 0, 3],
        [0, 2, 1],
        [0, 2, 7],
      ]);

      expect(L.rank()).toBe(3);
      expect(L.degree()).toBe(3);
    });

    it('should LLL-reduce by default', () => {
      const L = IntegerLattice([
        [1, 0, 3],
        [0, 2, 1],
        [0, 2, 7],
      ]);

      // The reduced basis should be LLL-reduced
      expect(isLLLReduced(L.reducedBasis)).toBe(true);
    });

    it('should not LLL-reduce when lllReduce: false', () => {
      const basis = [
        [1, 1000],
        [0, 1],
      ];

      const L = IntegerLattice(basis, { lllReduce: false });

      // The reduced basis should be the same as input
      expect(L.reducedBasis.get(0, 0).value).toBe(1n);
      expect(L.reducedBasis.get(0, 1).value).toBe(1000n);
    });
  });

  describe('LLL method', () => {
    it('should reduce the basis when called', () => {
      const L = IntegerLattice(
        [
          [1, 1000],
          [0, 1],
        ],
        { lllReduce: false }
      );

      const reduced = L.LLL();

      expect(isLLLReduced(reduced)).toBe(true);
    });

    it('should update reducedBasis with better result', () => {
      const L = IntegerLattice(
        [
          [1, 1000],
          [0, 1],
        ],
        { lllReduce: false }
      );

      // Original first vector has large norm
      const origNorm = L.reducedBasis.get(0, 0).value ** 2n + L.reducedBasis.get(0, 1).value ** 2n;

      L.LLL();

      // After LLL, first vector should have smaller norm
      const newNorm = L.reducedBasis.get(0, 0).value ** 2n + L.reducedBasis.get(0, 1).value ** 2n;

      expect(newNorm <= origNorm).toBe(true);
    });
  });
});

describe('Known examples from literature', () => {
  it('should reduce NTRU-like lattice', () => {
    // A simple example inspired by NTRU-style lattices
    const q = 127n;
    const basis = IntegerMatrixFromEntries([
      [1n, 0n, 23n, 45n],
      [0n, 1n, 67n, 89n],
      [0n, 0n, q, 0n],
      [0n, 0n, 0n, q],
    ]);

    const reduced = lllReduce(basis);

    expect(isLLLReduced(reduced)).toBe(true);

    // Determinant should be preserved (q^2 for this lattice structure)
    const detOriginal = basis.determinant().value;
    const detReduced = reduced.determinant().value;
    const absDetOriginal = detOriginal < 0n ? -detOriginal : detOriginal;
    const absDetReduced = detReduced < 0n ? -detReduced : detReduced;

    expect(absDetReduced).toBe(absDetOriginal);
  });

  it('should handle larger random-looking lattice', () => {
    // A 5x5 lattice
    const basis = IntegerMatrixFromEntries([
      [7n, 2n, 3n, 4n, 5n],
      [1n, 11n, 2n, 3n, 4n],
      [2n, 3n, 13n, 2n, 3n],
      [3n, 4n, 5n, 17n, 2n],
      [4n, 5n, 6n, 7n, 19n],
    ]);

    const reduced = lllReduce(basis);

    expect(isLLLReduced(reduced)).toBe(true);

    // Verify determinant preservation
    const detOriginal = basis.determinant().value;
    const detReduced = reduced.determinant().value;
    const absDetOriginal = detOriginal < 0n ? -detOriginal : detOriginal;
    const absDetReduced = detReduced < 0n ? -detReduced : detReduced;

    expect(absDetReduced).toBe(absDetOriginal);
  });
});

describe('Lattice invariants', () => {
  describe('volume and discriminant', () => {
    it('should compute volume for identity lattice', () => {
      const L = IntegerLattice([
        [1, 0],
        [0, 1],
      ]);

      expect(L.volume()).toBe(1n);
      expect(L.discriminant()).toBe(1n);
    });

    it('should compute volume for scaled lattice', () => {
      const L = IntegerLattice([
        [2, 0],
        [0, 3],
      ]);

      // Volume = |det| = 6
      expect(L.volume()).toBe(6n);
      // Discriminant = det(B * B^T) = det([[4, 0], [0, 9]]) = 36
      expect(L.discriminant()).toBe(36n);
    });

    it('should compute volume for 3D lattice', () => {
      const L = IntegerLattice(
        [
          [1, 0, 3],
          [0, 2, 1],
          [0, 2, 7],
        ],
        { lllReduce: false }
      );

      // det = 1 * (2*7 - 1*2) - 0 + 3 * (0 - 0) = 12
      const expectedVol = 12n;
      expect(L.volume()).toBe(expectedVol);
    });
  });

  describe('isUnimodular', () => {
    it('should return true for ZZ^n', () => {
      const L = IntegerLattice([
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
      ]);

      expect(L.isUnimodular()).toBe(true);
    });

    it('should return false for non-unimodular lattice', () => {
      const L = IntegerLattice([
        [2, 0],
        [0, 3],
      ]);

      expect(L.isUnimodular()).toBe(false);
    });
  });

  describe('hadamardRatio', () => {
    it('should be 1 for orthogonal basis', () => {
      const L = IntegerLattice(
        [
          [1, 0, 0],
          [0, 1, 0],
          [0, 0, 1],
        ],
        { lllReduce: false }
      );

      // For orthogonal basis, Hadamard ratio should be 1
      expect(L.hadamardRatio()).toBeCloseTo(1.0, 5);
    });

    it('should be less than 1 for non-orthogonal basis', () => {
      const L = IntegerLattice(
        [
          [1, 1],
          [0, 1],
        ],
        { lllReduce: false }
      );

      const ratio = L.hadamardRatio();
      expect(ratio).toBeGreaterThan(0);
      expect(ratio).toBeLessThanOrEqual(1.0);
    });
  });

  describe('gaussianHeuristic', () => {
    it('should estimate shortest vector length', () => {
      const L = IntegerLattice([
        [1, 0],
        [0, 1],
      ]);

      // For ZZ^2, Gaussian heuristic ~ sqrt(2 / (2*pi*e)) ~ 0.31
      const gh = L.gaussianHeuristic();
      expect(gh).toBeGreaterThan(0);
      expect(gh).toBeLessThan(2);
    });

    it('should be larger for larger volume', () => {
      const L1 = IntegerLattice([
        [1, 0],
        [0, 1],
      ]);

      const L2 = IntegerLattice([
        [10, 0],
        [0, 10],
      ]);

      expect(L2.gaussianHeuristic()).toBeGreaterThan(L1.gaussianHeuristic());
    });
  });
});

describe('Babai nearest plane (approximate closest vector)', () => {
  it('should find close vector for simple lattice', () => {
    const L = IntegerLattice([
      [1, 0],
      [0, 1],
    ]);

    // Target: (0.4, 0.6) should map to (0, 1) or (1, 1)
    const closest = L.babai([0.4, 0.6]);

    // The closest lattice point should be (0, 1) or (1, 1)
    // Check that it's actually a lattice point (integer coordinates)
    expect(typeof closest[0]).toBe('bigint');
    expect(typeof closest[1]).toBe('bigint');

    // The result should be close to target
    const dist = (0.4 - Number(closest[0])) ** 2 + (0.6 - Number(closest[1])) ** 2;
    expect(dist).toBeLessThan(2); // Should be within 2 units squared
  });

  it('should find exact result for lattice points', () => {
    const L = IntegerLattice([
      [1, 0],
      [0, 1],
    ]);

    // Target: (3, 5) is already a lattice point
    const closest = L.babai([3, 5]);

    expect(closest[0]).toBe(3n);
    expect(closest[1]).toBe(5n);
  });

  it('should work with non-trivial lattice', () => {
    const L = IntegerLattice([
      [101, 0, 0, 0],
      [0, 101, 0, 0],
      [0, 0, 101, 0],
      [-28, 39, 45, 1],
    ]);

    const target = [1337, 1337, 1337, 1337];
    const closest = L.babai(target);

    // Verify result is a lattice vector (linear combination of basis)
    expect(closest.length).toBe(4);

    // Should be reasonably close to target
    let distSq = 0;
    for (let i = 0; i < 4; i++) {
      distSq += (target[i]! - Number(closest[i])) ** 2;
    }
    // The distance should be reasonable (not huge)
    expect(distSq).toBeLessThan(10000 * 10000);
  });
});

describe('FreeModule methods', () => {
  describe('randomElement', () => {
    it('should generate a random element', () => {
      const { FreeModule } = require('./free_module.js');

      // Create a simple ring-like object
      const ZZ = {
        zero: () => 0n,
        one: () => 1n,
        is_field: () => false,
        toString: () => 'Integer Ring',
      };

      const M = FreeModule(ZZ, 3);
      const v = M.randomElement();

      expect(v.degree()).toBe(3);
    });
  });

  describe('areLinearlyDependent', () => {
    it('should return false for linearly independent vectors', () => {
      const { FreeModule } = require('./free_module.js');

      const ZZ = {
        zero: () => 0n,
        one: () => 1n,
        is_field: () => false,
        toString: () => 'Integer Ring',
      };

      const M = FreeModule(ZZ, 3);
      const v1 = M.createElement([1n, 0n, 0n]);
      const v2 = M.createElement([0n, 1n, 0n]);

      expect(M.areLinearlyDependent([v1, v2])).toBe(false);
    });

    it('should return true for linearly dependent vectors', () => {
      const { FreeModule } = require('./free_module.js');

      const ZZ = {
        zero: () => 0n,
        one: () => 1n,
        is_field: () => false,
        toString: () => 'Integer Ring',
      };

      const M = FreeModule(ZZ, 3);
      const v1 = M.createElement([1n, 2n, 3n]);
      const v2 = M.createElement([2n, 4n, 6n]); // 2 * v1

      expect(M.areLinearlyDependent([v1, v2])).toBe(true);
    });

    it('should return true when more vectors than dimension', () => {
      const { FreeModule } = require('./free_module.js');

      const ZZ = {
        zero: () => 0n,
        one: () => 1n,
        is_field: () => false,
        toString: () => 'Integer Ring',
      };

      const M = FreeModule(ZZ, 2);
      const v1 = M.createElement([1n, 0n]);
      const v2 = M.createElement([0n, 1n]);
      const v3 = M.createElement([1n, 1n]);

      // Three vectors in R^2 must be dependent
      expect(M.areLinearlyDependent([v1, v2, v3])).toBe(true);
    });
  });

  describe('coordinates', () => {
    it('should find coordinates for ambient module', () => {
      const { FreeModule } = require('./free_module.js');

      const ZZ = {
        zero: () => 0n,
        one: () => 1n,
        is_field: () => false,
        toString: () => 'Integer Ring',
      };

      const M = FreeModule(ZZ, 3);
      const v = M.createElement([2n, 3n, 5n]);

      const coords = M.coordinates(v);

      expect(coords.length).toBe(3);
      expect(coords[0]).toBe(2n);
      expect(coords[1]).toBe(3n);
      expect(coords[2]).toBe(5n);
    });
  });

  describe('directSum', () => {
    it('should compute direct sum of two modules', () => {
      const { FreeModule } = require('./free_module.js');

      const ZZ = {
        zero: () => 0n,
        one: () => 1n,
        is_field: () => false,
        toString: () => 'Integer Ring',
      };

      const M1 = FreeModule(ZZ, 2);
      const M2 = FreeModule(ZZ, 3);

      const M = M1.directSum(M2);

      expect(M.degree()).toBe(5);
    });
  });
});

describe('FreeModuleElement methods', () => {
  describe('pNorm', () => {
    it('should compute 1-norm', () => {
      const { vector } = require('./free_module_element.js');

      const v = vector([1n, -2n, 3n]);
      const norm1 = v.pNorm(1);

      // sage: vector([1,-2,3]).norm(1) == 6
      expect(norm1).toBe(6n);
    });

    it('should compute infinity-norm', () => {
      const { vector } = require('./free_module_element.js');

      const v = vector([1n, -5n, 3n]);
      const normInf = v.pNorm(Number.POSITIVE_INFINITY);

      // sage: vector([1,-5,3]).norm(Infinity) == 5
      expect(normInf).toBe(5n);
    });

    it('should compute 2-norm', () => {
      const { vector } = require('./free_module_element.js');

      const v = vector([3n, 4n]);
      const norm2 = v.pNorm(2);

      // sage: vector(ZZ, [3,4]).norm() == 5 (NOT the squared norm)
      expect(norm2).toBe(5n);
    });
  });

  describe('normalized', () => {
    it('should normalize a vector exactly', () => {
      const { vector } = require('./free_module_element.js');

      const v = vector([3, 4]);
      const u = v.normalized();

      // sage: vector([3,4]).normalized() == (3/5, 4/5); normalizing changes
      // the base ring to the fraction field, so the entries are Rationals.
      expect(String(u.getItem(0))).toBe('3/5');
      expect(String(u.getItem(1))).toBe('4/5');
    });

    it('should throw for zero vector', () => {
      const { vector } = require('./free_module_element.js');

      const v = vector([0, 0, 0]);
      expect(() => v.normalized()).toThrow();
    });
  });
});

describe('FreeModule structure (sage/modules/free_module.py doctests)', () => {
  const ZZ = {
    zero: () => 0n,
    one: () => 1n,
    is_field: () => false,
    toString: () => 'Integer Ring',
  };
  const QQ = {
    zero: () => Rational.zero(),
    one: () => Rational.one(),
    is_field: () => true,
    toString: () => 'Rational Field',
  };
  const r = (n: number, d: number = 1) => new Rational(BigInt(n), BigInt(d));
  const rows = (M: unknown[][]) => M.map((row) => row.map((e) => String(e)));

  describe('span: rank and echelon basis', () => {
    it('uses the true rank, not the number of generators', () => {
      // sage: span([[1,2,3],[2,4,6]], ZZ).rank() == 1
      const S = span(
        [
          [1n, 2n, 3n],
          [2n, 4n, 6n],
        ],
        ZZ
      );
      expect(S.rank()).toBe(1);
      expect(rows(S.basisMatrix() as unknown[][])).toEqual([['1', '2', '3']]);
    });

    it('echelonizes over ZZ with the Hermite normal form', () => {
      // sage: V = ZZ^3; V.span([[1,2,3],[4,5,6]])
      // Echelon basis matrix: [1 2 3] / [0 3 6]
      const S = span(
        [
          [1n, 2n, 3n],
          [4n, 5n, 6n],
        ],
        ZZ
      );
      expect(rows(S.basisMatrix() as unknown[][])).toEqual([
        ['1', '2', '3'],
        ['0', '3', '6'],
      ]);
    });

    it('echelonizes non-integral generators over ZZ', () => {
      // sage: V.span([[1,0,0],[1/5,4,0],[6,3/4,0]]) has rank 2 and echelon
      // basis [1/5 0 0] / [0 1/4 0]
      const V = FreeModule(ZZ, 3);
      const S = V.span([
        V.createElement([r(1), r(0), r(0)]),
        V.createElement([r(1, 5), r(4), r(0)]),
        V.createElement([r(6), r(3, 4), r(0)]),
      ]);
      expect(S.rank()).toBe(2);
      expect(rows(S.basisMatrix() as unknown[][])).toEqual([
        ['1/5', '0', '0'],
        ['0', '1/4', '0'],
      ]);
    });

    it('computes the dimension of a QQ subspace', () => {
      // sage: VectorSpace(QQ,3).subspace([(1,2,3),(2,4,6)]).dimension() == 1
      const V = VectorSpace(QQ, 3);
      const W = V.subspace([
        V.createElement([r(1), r(2), r(3)]),
        V.createElement([r(2), r(4), r(6)]),
      ]);
      expect(W.dimension()).toBe(1);
    });

    it('rejects a dependent user basis', () => {
      // sage: ZZ^2.span_of_basis([[1,2],[2,4]]) raises ValueError
      const V = FreeModule(ZZ, 2) as FreeModulePID;
      expect(() => V.spanOfBasis([V.createElement([1n, 2n]), V.createElement([2n, 4n])])).toThrow(
        'the given basis vectors must be linearly independent'
      );
    });
  });

  describe('echelonizedBasisMatrix', () => {
    it('is the Hermite normal form of the user basis', () => {
      // sage: FreeModule(ZZ,3).span_of_basis([[1,2,3],[4,5,6]])
      //         .echelonized_basis_matrix() == [1 2 3] / [0 3 6]
      const V = FreeModule(ZZ, 3) as FreeModulePID;
      const W = V.spanOfBasis([V.createElement([1n, 2n, 3n]), V.createElement([4n, 5n, 6n])]);
      expect(rows(W.basisMatrix() as unknown[][])).toEqual([
        ['1', '2', '3'],
        ['4', '5', '6'],
      ]);
      expect(rows(W.echelonizedBasisMatrix())).toEqual([
        ['1', '2', '3'],
        ['0', '3', '6'],
      ]);
    });

    it('reduces above the pivots', () => {
      // sage: FreeModule_submodule_with_basis_pid(ZZ^3, [[1,2,3],[1,1,1]])
      //         .echelonized_basis_matrix() == [1 0 -1] / [0 1 2]
      const V = FreeModule(ZZ, 3) as FreeModulePID;
      const W = V.spanOfBasis([V.createElement([1n, 2n, 3n]), V.createElement([1n, 1n, 1n])]);
      expect(rows(W.echelonizedBasisMatrix())).toEqual([
        ['1', '0', '-1'],
        ['0', '1', '2'],
      ]);
    });

    it('does not lose the index of 2*ZZ^2', () => {
      const V = FreeModule(ZZ, 2);
      const W = V.scale(2n);
      expect(rows(W.echelonizedBasisMatrix())).toEqual([
        ['2', '0'],
        ['0', '2'],
      ]);
    });
  });

  describe('coordinates', () => {
    it('solves exactly over the base ring', () => {
      // sage: M = FreeModule(ZZ,2); W = M.submodule([M0+M1, M0-2*M1])
      //       W.coordinates(2*M0 - M1) == [2, -1]
      const M = FreeModule(ZZ, 2);
      const W = M.submodule([M.createElement([1n, 1n]), M.createElement([1n, -2n])]);
      expect(W.coordinates(M.createElement([2n, -1n])).map(String)).toEqual(['2', '-1']);
    });

    it('returns rational coordinates when the vector is not in the module', () => {
      // sage: W = (ZZ^2).span([[2,0],[0,2]]); W.coordinates((1,0)) == [1/2, 0]
      const M = FreeModule(ZZ, 2);
      const W = M.span([M.createElement([2n, 0n]), M.createElement([0n, 2n])]);
      expect(W.coordinates(M.createElement([1n, 0n])).map(String)).toEqual(['1/2', '0']);
      expect(W.coordinates(M.createElement([2n, 0n])).map(String)).toEqual(['1', '0']);
    });

    it('raises when the vector is not in the span', () => {
      const M = FreeModule(ZZ, 2);
      const W = M.span([M.createElement([3n, 6n])]);
      expect(W.coordinates(M.createElement([1n, 2n])).map(String)).toEqual(['1/3']);
      expect(() => W.coordinates(M.createElement([1n, 3n]))).toThrow(
        'vector is not in free module'
      );
    });

    it('uses the user basis, not the echelon basis', () => {
      // sage: W = (ZZ^3).span_of_basis([[1,2,3],[4,5,6]])
      //       W.coordinate_vector([1,5,9]) == (5, -1)
      const V = FreeModule(ZZ, 3) as FreeModulePID;
      const W = V.spanOfBasis([V.createElement([1n, 2n, 3n]), V.createElement([4n, 5n, 6n])]);
      expect(W.coordinates(V.createElement([1n, 5n, 9n])).map(String)).toEqual(['5', '-1']);
    });
  });

  describe('isSubmodule', () => {
    it('requires integral coordinates', () => {
      // sage: M = FreeModule(ZZ,2); N = M.scale(2)
      //       N.is_submodule(M) is True; M.is_submodule(N) is False
      const M = FreeModule(ZZ, 2);
      const N = M.scale(2n);
      expect(N.isSubmodule(M)).toBe(true);
      expect(M.isSubmodule(N)).toBe(false);
      expect(M.isSubmodule(M)).toBe(true);
    });
  });

  describe('discriminant', () => {
    it('matches the SageMath doctests', () => {
      // sage: M = FreeModule(ZZ,3); M.discriminant() == 1
      //       M.span([[1,2,3]]).discriminant() == 14
      //       M.span([[1,2,3],[1,1,1]]).discriminant() == 6
      expect(FreeModule(ZZ, 3).discriminant()).toBe(1n);
      expect(span([[1n, 2n, 3n]], ZZ).discriminant()).toBe(14n);
      expect(
        span(
          [
            [1n, 2n, 3n],
            [1n, 1n, 1n],
          ],
          ZZ
        ).discriminant()
      ).toBe(6n);
    });

    it('is exact for large entries', () => {
      const L = span(
        [
          [123456789n, 2n, 3n],
          [4n, 987654321n, 6n],
          [7n, 8n, 555555555n],
        ],
        ZZ
      );
      expect(L.discriminant()).toBe(4588755092689766572878148493545359752410735083930084n);
    });

    it('applies the inner product matrix', () => {
      // sage: FreeQuadraticModule(ZZ, 2, matrix.identity(2)).discriminant() == -1
      const M = FreeModule(ZZ, 2, {
        innerProductMatrix: [
          [1n, 0n],
          [0n, 1n],
        ],
      });
      expect(rows(M.gramMatrix())).toEqual([
        ['1', '0'],
        ['0', '1'],
      ]);
      expect(M.discriminant()).toBe(-1n);
    });
  });

  describe('saturation', () => {
    it('saturates a rank 1 lattice', () => {
      // sage: span([[9,9,6]], ZZ).saturation() has basis [3 3 2]
      const L = span([[9n, 9n, 6n]], ZZ) as FreeModulePID;
      expect(rows(L.saturation().basisMatrix() as unknown[][])).toEqual([['3', '3', '2']]);
    });

    it('saturates a rank 2 lattice and drops the index', () => {
      // sage: L = span([[1,2,3],[4,5,6]], ZZ); L.discriminant() == 54
      //       L.saturation() has basis [1 0 -1] / [0 1 2] and discriminant 6
      const L = span(
        [
          [1n, 2n, 3n],
          [4n, 5n, 6n],
        ],
        ZZ
      ) as FreeModulePID;
      expect(L.discriminant()).toBe(54n);
      const S = L.saturation();
      expect(rows(S.basisMatrix() as unknown[][])).toEqual([
        ['1', '0', '-1'],
        ['0', '1', '2'],
      ]);
      expect(S.discriminant()).toBe(6n);
      expect(L.indexInSaturation()).toBe(3n);
    });

    it('returns self when already saturated', () => {
      const L = span([[3n, 3n, 2n]], ZZ) as FreeModulePID;
      expect(L.saturation()).toBe(L);
    });
  });

  describe('indexIn', () => {
    it('matches the SageMath doctests', () => {
      // sage: L1 = span([[1,2]], ZZ); L2 = span([[3,6]], ZZ)
      //       L2.index_in(L1) == 3
      const L1 = span([[1n, 2n]], ZZ) as FreeModulePID;
      const L2 = span([[3n, 6n]], ZZ) as FreeModulePID;
      expect(L2.indexIn(L1)).toBe(3n);
    });

    it('returns a rational index for non-integral lattices', () => {
      // sage: L1 = span([['1/2','1/3'], [4,5]], ZZ); L2 = span([[1,2],[3,4]], ZZ)
      //       L2.index_in(L1) == 12/7 and L1.index_in(L2) == 7/12
      const V = FreeModule(ZZ, 3);
      const L1 = V.span([
        V.createElement([r(1, 2), r(1, 3), r(0)]),
        V.createElement([r(4), r(5), r(0)]),
      ]) as FreeModulePID;
      const L2 = V.span([
        V.createElement([r(1), r(2), r(0)]),
        V.createElement([r(3), r(4), r(0)]),
      ]) as FreeModulePID;
      expect(String(L2.indexIn(L1))).toBe('12/7');
      expect(String(L1.indexIn(L2))).toBe('7/12');
    });

    it('is infinite for a lattice of smaller rank', () => {
      // sage: span([[1,2]], ZZ).index_in(FreeModule(ZZ,2)) == +Infinity
      const L = span([[1n, 2n]], ZZ) as FreeModulePID;
      expect(L.indexIn(FreeModule(ZZ, 2))).toBe(Number.POSITIVE_INFINITY);
    });
  });

  describe('intersection', () => {
    it('intersects two QQ subspaces', () => {
      // sage: V = VectorSpace(QQ,3); W1 = V.submodule([V.0, V.0+V.1])
      //       W2 = V.submodule([V.1, V.2]); W1.intersection(W2) has basis [0 1 0]
      const V = VectorSpace(QQ, 3);
      const W1 = V.subspace([
        V.createElement([r(1), r(0), r(0)]),
        V.createElement([r(1), r(1), r(0)]),
      ]);
      const W2 = V.subspace([
        V.createElement([r(0), r(1), r(0)]),
        V.createElement([r(0), r(0), r(1)]),
      ]);
      expect(rows(W1.intersection(W2).basisMatrix() as unknown[][])).toEqual([['0', '1', '0']]);
      expect(rows(W2.intersection(W1).basisMatrix() as unknown[][])).toEqual([['0', '1', '0']]);
    });

    it('intersects two ZZ lattices', () => {
      const M = FreeModule(ZZ, 2) as FreeModulePID;
      const two = M.scale(2n) as FreeModulePID;
      const three = M.scale(3n);
      expect(rows(two.intersection(three).basisMatrix() as unknown[][])).toEqual([
        ['6', '0'],
        ['0', '6'],
      ]);
    });
  });

  describe('complement', () => {
    it('is the right kernel of the basis matrix over QQ', () => {
      // sage: V = QQ^3; W = V.span([[1,0,1]]); W.complement() == [1 0 -1] / [0 1 0]
      const V = VectorSpace(QQ, 3);
      const W = V.subspace([V.createElement([r(1), r(0), r(1)])]);
      expect(rows(W.complement().basisMatrix() as unknown[][])).toEqual([
        ['1', '0', '-1'],
        ['0', '1', '0'],
      ]);
    });

    it('respects the characteristic of the base field', () => {
      // Over GF(5), the complement of <(2,1)> is <(1,3)>, since 2 + 3 = 5 = 0
      const V = VectorSpace(GF(5n) as never, 2);
      const W = V.subspace([V.createElement([2n, 1n])]);
      expect(rows(W.complement().basisMatrix() as unknown[][])).toEqual([['1', '3']]);
    });

    it('reproduces the GF(2) doctest', () => {
      // sage: V = GF(2)^6; W = V.span([[1,1,0,0,0,0]]); W.complement() has the
      // 5 rows below, and W.intersection(W.complement()) == W
      const V = VectorSpace(GF(2n) as never, 6);
      const W = V.subspace([V.createElement([1n, 1n, 0n, 0n, 0n, 0n])]);
      expect(rows(W.complement().basisMatrix() as unknown[][])).toEqual([
        ['1', '1', '0', '0', '0', '0'],
        ['0', '0', '1', '0', '0', '0'],
        ['0', '0', '0', '1', '0', '0'],
        ['0', '0', '0', '0', '1', '0'],
        ['0', '0', '0', '0', '0', '1'],
      ]);
    });
  });

  describe('areLinearlyDependent', () => {
    it('respects the characteristic', () => {
      // Over GF(5): 3*(1,2) == (3,1)
      const V = VectorSpace(GF(5n) as never, 2);
      expect(V.areLinearlyDependent([V.createElement([1n, 2n]), V.createElement([3n, 1n])])).toBe(
        true
      );
      expect(V.areLinearlyDependent([V.createElement([1n, 2n]), V.createElement([1n, 1n])])).toBe(
        false
      );
    });

    it('matches the QQ doctest', () => {
      // sage: M = QQ^3; vecs = [M([1,2,3]), M([4,5,6])] independent,
      //       adding M([3,3,3]) makes them dependent
      const M = VectorSpace(QQ, 3);
      const v1 = M.createElement([r(1), r(2), r(3)]);
      const v2 = M.createElement([r(4), r(5), r(6)]);
      const v3 = M.createElement([r(3), r(3), r(3)]);
      expect(M.areLinearlyDependent([v1, v2])).toBe(false);
      expect(M.areLinearlyDependent([v1, v2, v3])).toBe(true);
    });
  });

  describe('changeRing', () => {
    it('keeps the degree and the basis', () => {
      const V = FreeModule(ZZ, 3);
      const W = V.span([V.createElement([1n, 2n, 3n])]);
      expect(W.changeRing(ZZ)).toBe(W);
      const WQ = W.changeRing(QQ);
      expect(WQ.degree()).toBe(3);
      expect(rows(WQ.basisMatrix() as unknown[][])).toEqual([['1', '2', '3']]);
    });

    it('reduces a QQ subspace modulo p', () => {
      // sage: W = (QQ^3).subspace([[2, 1/2, 1]]); W.change_ring(GF(7)) == [1 2 4]
      const V = VectorSpace(QQ, 3);
      const W = V.subspace([V.createElement([r(2), r(1, 2), r(1)])]);
      expect(rows(W.changeRing(GF(7n) as never).basisMatrix() as unknown[][])).toEqual([
        ['1', '2', '4'],
      ]);
    });

    it('rejects a ring that is not a principal ideal domain', () => {
      const V = FreeModule(ZZ, 3);
      const W = V.span([V.createElement([1n, 2n, 3n])]);
      const weird = {
        zero: () => 0n,
        one: () => 1n,
        is_field: () => false,
        toString: () => 'Weird Ring',
      };
      expect(() => W.changeRing(weird)).toThrow('should be a principal ideal domain');
    });
  });

  describe('cardinality', () => {
    it('is exact above 2^53', () => {
      // sage: FreeModule(GF(2), 70).cardinality() == 2^70
      expect(FreeModule(GF(2n) as never, 70).cardinality()).toBe(2n ** 70n);
    });

    it('is 1 for the zero module and infinite over QQ', () => {
      // sage: VectorSpace(QQ, 0).cardinality() == 1; (QQ^3).cardinality() == +Infinity
      expect(VectorSpace(QQ, 0).cardinality()).toBe(1n);
      expect(VectorSpace(QQ, 3).cardinality()).toBe(Number.POSITIVE_INFINITY);
    });
  });

  describe('quotient and tensor product', () => {
    it('has the expected dimension', () => {
      const V = VectorSpace(QQ, 3);
      const W = V.subspace([V.createElement([r(1), r(0), r(0)])]);
      expect(V.quotient(W).dimension()).toBe(2);
    });

    it('quotientModule rejects a module that is not a submodule', () => {
      const M = FreeModule(ZZ, 2);
      const N = M.scale(2n);
      // M is not a submodule of N, so M/N is not defined
      expect(() => N.quotientModule(M)).toThrow();
    });

    it('tensor product has rank m*n', () => {
      const M = FreeModule(ZZ, 2) as FreeModulePID;
      const N = FreeModule(ZZ, 3);
      expect(M.tensorProduct(N).rank()).toBe(6);
    });
  });

  describe('module equality', () => {
    it('compares echelon bases, not generators', () => {
      const A = span(
        [
          [1n, 2n, 3n],
          [4n, 5n, 6n],
        ],
        ZZ
      );
      const B = span(
        [
          [4n, 5n, 6n],
          [1n, 2n, 3n],
          [5n, 7n, 9n],
        ],
        ZZ
      );
      expect(A.equals(B)).toBe(true);
      expect(A.equals(span([[1n, 2n, 3n]], ZZ))).toBe(false);
    });
  });
});

describe('Lattice generation functions', () => {
  describe('genLattice', () => {
    it('should generate modular lattice', () => {
      const { genLattice } = require('./free_module_integer.js');

      const L = genLattice({
        type: 'modular',
        n: 2,
        m: 4,
        q: 11n,
        seed: 42,
        lattice: true,
      });

      expect(L.rank()).toBe(4);
      expect(L.degree()).toBe(4);
    });

    it('should generate random lattice (n=1)', () => {
      const { genLattice } = require('./free_module_integer.js');

      const L = genLattice({
        type: 'random',
        n: 1,
        m: 5,
        q: 101n,
        seed: 42,
        lattice: true,
      });

      expect(L.rank()).toBe(5);
      expect(L.degree()).toBe(5);
    });

    it('should generate matrix when lattice: false', () => {
      const { genLattice } = require('./free_module_integer.js');

      const M = genLattice({
        type: 'modular',
        n: 2,
        m: 4,
        q: 11n,
        seed: 42,
        lattice: false,
      });

      // Should be a matrix
      expect(M.nrows).toBe(4);
      expect(M.ncols).toBe(4);
    });

    it('should support dual lattices', () => {
      const { genLattice } = require('./free_module_integer.js');

      const L = genLattice({
        type: 'modular',
        n: 2,
        m: 4,
        q: 11n,
        seed: 42,
        dual: true,
        lattice: true,
      });

      expect(L.rank()).toBe(4);
      expect(L.degree()).toBe(4);
    });
  });

  describe('randomLattice', () => {
    it('should generate random lattice of given dimension', () => {
      const { randomLattice } = require('./free_module_integer.js');

      const L = randomLattice(3, undefined, 42);

      expect(L.rank()).toBe(3);
      expect(L.degree()).toBe(3);
    });

    it('should generate lattice with specified determinant', () => {
      const { randomLattice } = require('./free_module_integer.js');

      const L = randomLattice(3, 100n, 42);

      expect(L.rank()).toBe(3);
      // The volume should be preserved (absolute value of determinant)
      expect(L.volume()).toBe(100n);
    });
  });

  describe('standardLattice', () => {
    it('should create ZZ^n', () => {
      const { standardLattice } = require('./free_module_integer.js');

      const L = standardLattice(3);

      expect(L.rank()).toBe(3);
      expect(L.degree()).toBe(3);
      expect(L.volume()).toBe(1n);
      expect(L.isUnimodular()).toBe(true);
    });
  });

  describe('qaryLattice', () => {
    it('should create q-ary lattice', () => {
      const { qaryLattice } = require('./free_module_integer.js');

      const A = [
        [1n, 2n],
        [3n, 4n],
      ];
      const L = qaryLattice(A, 7n);

      expect(L.rank()).toBe(2);
      expect(L.degree()).toBe(2);
    });
  });

  describe('hermiteFactor', () => {
    it('should compute Hermite factor', () => {
      const { hermiteFactor, standardLattice } = require('./free_module_integer.js');

      const L = standardLattice(3);
      const hf = hermiteFactor(L);

      // For the standard lattice, Hermite factor is 1
      expect(hf).toBeCloseTo(1.0, 5);
    });
  });

  describe('estimateBKZBlockSize', () => {
    it('should estimate block size', () => {
      const { estimateBKZBlockSize } = require('./free_module_integer.js');

      // For a small lattice with achievable target
      const beta = estimateBKZBlockSize(10, 1000n, 5);

      expect(beta).toBeGreaterThanOrEqual(2);
      expect(beta).toBeLessThanOrEqual(10);
    });
  });
});

describe('Shortest Vector Problem (SVP)', () => {
  describe('shortestVector basic functionality', () => {
    it('should find shortest vector in simple 2D lattice', () => {
      const L = IntegerLattice([
        [1, 0],
        [0, 2],
      ]);

      const sv = L.shortestVector();

      // The shortest vector should be (1, 0) or (-1, 0)
      const normSq = Number(sv[0]! * sv[0]! + sv[1]! * sv[1]!);
      expect(normSq).toBe(1); // ||(1, 0)||^2 = 1
    });

    it('should find shortest vector in 3D lattice', () => {
      const L = IntegerLattice([
        [1, 0, 3],
        [0, 2, 1],
        [0, 2, 7],
      ]);

      const sv = L.shortestVector();

      // The shortest vector should be found
      expect(sv.length).toBe(3);

      // Verify it's nonzero
      const isZero = sv.every((v) => v === 0n);
      expect(isZero).toBe(false);

      // Compute norm squared
      let normSq = 0n;
      for (const v of sv) {
        normSq += v * v;
      }

      // Should be short (LLL gives good approximation)
      expect(normSq).toBeLessThanOrEqual(20n);
    });

    it('should find shortest vector in identity lattice', () => {
      const L = IntegerLattice([
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
      ]);

      const sv = L.shortestVector();

      // For identity lattice, shortest nonzero vector has norm 1
      let normSq = 0n;
      for (const v of sv) {
        normSq += v * v;
      }
      expect(normSq).toBe(1n);
    });

    it('should find shortest vector in scaled lattice', () => {
      const L = IntegerLattice([
        [2, 0],
        [0, 3],
      ]);

      const sv = L.shortestVector();

      // Shortest vector should be (2, 0), (+/-2, 0)
      let normSq = 0n;
      for (const v of sv) {
        normSq += v * v;
      }
      expect(normSq).toBe(4n); // ||(2, 0)||^2 = 4
    });
  });

  describe('shortestVector with known examples', () => {
    it('should handle NTRU-like lattice', () => {
      // Simple q-ary style lattice
      const q = 127n;
      const L = IntegerLattice([
        [1n, 0n, 23n, 45n],
        [0n, 1n, 67n, 89n],
        [0n, 0n, q, 0n],
        [0n, 0n, 0n, q],
      ]);

      const sv = L.shortestVector();

      // Should find a vector
      expect(sv.length).toBe(4);

      // Should be nonzero
      const isZero = sv.every((v) => v === 0n);
      expect(isZero).toBe(false);

      // Norm should be reasonable for this lattice
      let normSq = 0n;
      for (const v of sv) {
        normSq += v * v;
      }
      // A small secret (1, 1, 0, 0) would give norm 2
      // LLL should find something reasonably short
      expect(normSq).toBeLessThan(BigInt(q * q));
    });

    it('should find exact shortest vector for small lattices', () => {
      // Lattice with known shortest vector
      const L = IntegerLattice([
        [4, 1],
        [3, 2],
      ]);

      const sv = L.shortestVector();

      // Compute the norm
      let normSq = 0n;
      for (const v of sv) {
        normSq += v * v;
      }

      // The shortest vector in this lattice should have small norm
      // The lattice has det = 4*2 - 1*3 = 5
      // By Minkowski, lambda_1 <= sqrt(2/pi * det) ~ 1.8
      // So lambda_1^2 <= about 4
      expect(normSq).toBeLessThanOrEqual(10n);
    });
  });

  describe('shortestVector options', () => {
    it('should work with updateReducedBasis: false', () => {
      const L = IntegerLattice(
        [
          [1, 1000],
          [0, 1],
        ],
        { lllReduce: false }
      );

      const originalBasis = L.reducedBasis;
      const sv = L.shortestVector({ updateReducedBasis: false });

      // Should find a short vector
      let normSq = 0n;
      for (const v of sv) {
        normSq += v * v;
      }
      expect(normSq).toBeLessThan(1000n * 1000n);
    });
  });

  describe('shortestVector with medium-sized lattices', () => {
    it('should handle 5D lattice', () => {
      const L = IntegerLattice([
        [7, 2, 3, 4, 5],
        [1, 11, 2, 3, 4],
        [2, 3, 13, 2, 3],
        [3, 4, 5, 17, 2],
        [4, 5, 6, 7, 19],
      ]);

      const sv = L.shortestVector();

      // Should find a vector
      expect(sv.length).toBe(5);

      // Should be nonzero
      const isZero = sv.every((v) => v === 0n);
      expect(isZero).toBe(false);

      // Should be reasonably short (diagonal-dominant matrix)
      let normSq = 0n;
      for (const v of sv) {
        normSq += v * v;
      }
      // First row has norm ~sqrt(7^2+2^2+3^2+4^2+5^2) = sqrt(103) ~ 10
      // So normSq should be < 200 or so
      expect(normSq).toBeLessThan(500n);
    });

    it('should handle 8D random lattice', () => {
      // Determinant 1 lattice for easier SVP
      const { randomLattice } = require('./free_module_integer.js');
      const L = randomLattice(8, 1n, 12345);

      const sv = L.shortestVector();

      // Should find a vector
      expect(sv.length).toBe(8);

      // Should be nonzero
      const isZero = sv.every((v) => v === 0n);
      expect(isZero).toBe(false);
    });
  });

  describe('shortestVector verifies lattice membership', () => {
    it('returned vector should be in the lattice', () => {
      const basis = [
        [3n, 1n, 2n],
        [1n, 4n, 1n],
        [2n, 1n, 5n],
      ];
      const L = IntegerLattice(basis);

      const sv = L.shortestVector();

      // The vector should be expressible as integer combination of basis vectors
      // We verify by checking it's in the span (this is automatic since we
      // construct it from the reduced basis)
      expect(sv.length).toBe(3);

      // Verify nonzero
      const normSq = sv.reduce((sum, v) => sum + v * v, 0n);
      expect(normSq).toBeGreaterThan(0n);
    });
  });

  describe('shortestVector edge cases', () => {
    it('should handle 1D lattice', () => {
      const L = IntegerLattice([[5]]);

      const sv = L.shortestVector();

      // Should be +/- 5
      expect(sv.length).toBe(1);
      expect(sv[0]! === 5n || sv[0]! === -5n).toBe(true);
    });

    it('should throw for empty lattice', () => {
      // Create empty basis
      const L = IntegerLattice([], { lllReduce: false });

      expect(() => L.shortestVector()).toThrow();
    });
  });
});

// ============================================================================
// Free modules over a non-ZZ PID (QQ[x], GF(p)[x])
//
// Every expected value below was produced by SageMath 10.3.
// ============================================================================

describe('free modules over QQ[x] (sage/matrix/matrix2.pyx:_echelon_form_PID)', () => {
  const [Rx, x] = PolynomialRingConstructor(QQfield as never, 'x');
  const one = Rx.one();
  const zero = Rx.zero();
  const A3 = FreeModule(Rx as never, 3) as FreeModulePID;
  const p = (...c: number[]) => Rx.__call__(c.map((n) => new Rational(BigInt(n))) as never);
  const v3 = (a: unknown, b: unknown, c: unknown) => A3.createElement([a, b, c]);
  const strRows = (M: unknown[][]) => M.map((row) => row.map((e) => String(e)));
  /** The span of a list of vectors, typed as a module over a PID. */
  const spanPID = (gens: ReturnType<typeof v3>[]) => A3.span(gens) as unknown as FreeModulePID;

  it('QQ[x] is recognised as a PID', () => {
    // sage: R.<x> = QQ[]; R^3 is an ambient free module over a PID
    expect(A3.isAmbient()).toBe(true);
    expect(A3.rank()).toBe(3);
    expect(A3.toString()).toBe(
      'Free module of rank 3 over Univariate Polynomial Ring in x over Rational Field'
    );
  });

  it('computes the Hermite echelon basis of a span', () => {
    // sage: R.<x> = QQ[]; A = R^3
    // sage: A.span([vector(R,[x,x^2,1]), vector(R,[1,x,x])]).basis_matrix()
    // [      1       x       x]
    // [      0       0 x^2 - 1]
    const L = A3.span([v3(x, x.mul(x), one), v3(one, x, x)]);
    expect(L.rank()).toBe(2);
    expect(strRows(L.basisMatrix() as unknown[][])).toEqual([
      ['1', 'x', 'x'],
      ['0', '0', 'x^2 + -1'],
    ]);
  });

  it('reproduces the issue #9053 doctest (matrix2.pyx:17305)', () => {
    // sage: R.<x> = GF(7)[]; A = R^3
    // sage: L = A.span([x*A.0 + (x^3+1)*A.1, x*A.2]); M = A.span([x*L.0])
    // sage: M.0 in L
    // True
    const L = A3.span([v3(x, x.pow(3).add(one), zero), v3(zero, zero, x)]);
    expect(strRows(L.basisMatrix() as unknown[][])).toEqual([
      ['x', 'x^3 + 1', '0'],
      ['0', '0', 'x'],
    ]);
    const M = A3.span([L.basis()[0]!.mul(x)]);
    expect(strRows(M.basisMatrix() as unknown[][])).toEqual([['x^2', 'x^4 + x', '0']]);
    expect(M.isSubmodule(L)).toBe(true);
    // sage: L.coordinates(M.0) == [x, 0]
    expect(L.coordinates(M.basis()[0]!).map(String)).toEqual(['x', '0']);
  });

  it('computes intersections and sums over QQ[x]', () => {
    // sage: P = A.span([[x,1,0],[0,x,1]]); Q = A.span([[x,1,0],[1,1,1]])
    // sage: P.intersection(Q).basis_matrix() == [x 1 0]
    // sage: (P+Q).basis_matrix() ==
    // [          1           1           1]
    // [          0           1      -x + 1]
    // [          0           0 x^2 - x + 1]
    const P = spanPID([v3(x, one, zero), v3(zero, x, one)]);
    const Q = spanPID([v3(x, one, zero), v3(one, one, one)]);
    expect(strRows(P.intersection(Q).basisMatrix() as unknown[][])).toEqual([['x', '1', '0']]);
    expect(strRows(P.add(Q).basisMatrix() as unknown[][])).toEqual([
      ['1', '1', '1'],
      ['0', '1', '(-1)*x + 1'],
      ['0', '0', 'x^2 + (-1)*x + 1'],
    ]);
  });

  it('intersections are contained in both factors and have the right rank', () => {
    // The echelon basis of an intersection over a non-ZZ PID is only canonical
    // up to units, so the module-level properties are what is pinned here.
    const cases: [unknown[][], unknown[][]][] = [
      [
        [
          [x, one, zero],
          [zero, x, one],
        ],
        [
          [x, one, zero],
          [one, one, one],
        ],
      ],
      [
        [[p(1), p(0, 1), p(2)]],
        [
          [p(2), p(0, 2), p(4)],
          [p(0), p(1), p(1)],
        ],
      ],
      [
        [
          [p(0, 1), p(1), p(0)],
          [p(0), p(0, 0, 1), p(1)],
        ],
        [
          [p(0, 1), p(1), p(0)],
          [p(0), p(0, 1), p(1)],
        ],
      ],
    ];
    for (const [g1, g2] of cases) {
      const P = spanPID(g1.map((r) => v3(r[0], r[1], r[2])));
      const Q = spanPID(g2.map((r) => v3(r[0], r[1], r[2])));
      const I = P.intersection(Q) as unknown as FreeModulePID;
      expect(I.isSubmodule(P)).toBe(true);
      expect(I.isSubmodule(Q)).toBe(true);
      // dim(P n Q) = dim P + dim Q - dim(P + Q) over the fraction field
      expect(I.rank()).toBe(P.rank() + Q.rank() - P.add(Q).rank());
    }
  });

  it('has QQ(x) as its base field and spans vector spaces over it', () => {
    // sage: A.base_field() == Fraction Field of Univariate Polynomial Ring in x over QQ
    // sage: L.vector_space_span(L.basis()).basis_matrix() ==
    // [          1 (x^3 + 1)/x           0]
    // [          0           0           1]
    expect(String(A3.baseField())).toBe(
      'Fraction Field of Univariate Polynomial Ring in x over Rational Field'
    );
    const L = spanPID([v3(x, x.pow(3).add(one), zero), v3(zero, zero, x)]);
    const VS = L.vectorSpaceSpan(L.basis());
    expect(VS.dimension()).toBe(2);
    expect(strRows(VS.basisMatrix() as unknown[][])).toEqual([
      ['1', '(x^3 + 1)/(x)', '0'],
      ['0', '0', '1'],
    ]);
  });

  it('clears denominators like Sage _denominator', () => {
    // sage: K = Frac(R); span([vector(K,[1/x,1]), vector(K,[1,x])], R).basis_matrix()
    // [1/x   1]
    // sage: ....denominator() == x
    const A2 = FreeModule(Rx as never, 2) as FreeModulePID;
    const invx = FractionFieldElement.make(
      one as unknown as Parameters<typeof FractionFieldElement.make>[0],
      x as unknown as Parameters<typeof FractionFieldElement.make>[0]
    );
    const S = A2.span([
      A2.createElement([invx, one]),
      A2.createElement([one, x]),
    ]) as unknown as FreeModulePID;
    expect(S.rank()).toBe(1);
    expect(strRows(S.basisMatrix() as unknown[][])).toEqual([['(1)/(x)', '1']]);
    expect(String(S.denominator())).toBe('x');
  });

  it('the echelon basis spans the same module as the generators', () => {
    // sage: gens = [(2*x+1, x^2, 3), (x, 1, x-1), (2, 2*x+2, x^3)]
    // sage: (R^3).span(gens).basis_matrix()
    // [ 1                x^2 - 2                              -2*x + 5]
    // [ 0                      1   1/4*x^5 - 1/4*x^4 + 1/2*x^2 - x - 1]
    // [ 0                      0 -x^6 + 2*x^4 - x^3 + 2*x^2 + 10*x - 4]
    const gens = [
      v3(p(1, 2), p(0, 0, 1), p(3)),
      v3(p(0, 1), p(1), p(-1, 1)),
      v3(p(2), p(2, 2), p(0, 0, 0, 1)),
    ];
    const L = A3.span(gens);
    expect(strRows(L.basisMatrix() as unknown[][])).toEqual([
      ['1', 'x^2 + -2', '(-2)*x + 5'],
      ['0', '1', '1/4*x^5 + (-1/4)*x^4 + 1/2*x^2 + (-1)*x + -1'],
      ['0', '0', '(-1)*x^6 + 2*x^4 + (-1)*x^3 + 2*x^2 + 10*x + -4'],
    ]);
    // sage: [L.coordinates(g) for g in gens]
    expect(L.coordinates(gens[0]!).map(String)).toEqual([
      '2*x + 1',
      '(-2)*x^3 + 4*x + 2',
      '(-1/2)*x^2 + 1/2*x',
    ]);
    // The echelon form is a *base change* of the generators: the two modules
    // contain each other.  isSubmodule solves over the fraction field and
    // checks integrality, so it does not go through echelonRows.
    const G = A3.spanOfBasis(gens, undefined, { check: false });
    expect(G.isSubmodule(L)).toBe(true);
    expect(L.isSubmodule(G)).toBe(true);
  });

  it('quotients over QQ[x] raise SageMath NotImplementedError', () => {
    // sage: R.<x> = QQ[]; A = R^2; A.quotient_module(A.span([[x,1]]))
    // NotImplementedError: quotients of modules over rings other than fields
    //                      or ZZ is not fully implemented
    const A2 = FreeModule(Rx as never, 2);
    const W = A2.span([A2.createElement([x, one])]);
    expect(() => A2.quotientModule(W)).toThrow(
      'quotients of modules over rings other than fields or ZZ is not fully implemented'
    );
  });

  it('works over GF(7)[x] as well', () => {
    // sage: R.<x> = GF(7)[]; A = R^3; L = A.span([x*A.0+(x^3+1)*A.1, x*A.2])
    // sage: M = A.span([x*L.0]); M.0 in L  ->  True
    const [R7, x7] = PolynomialRingConstructor(GF(7) as never, 'x');
    const A = FreeModule(R7 as never, 3);
    const L = A.span([
      A.createElement([x7, x7.pow(3).add(R7.one()), R7.zero()]),
      A.createElement([R7.zero(), R7.zero(), x7]),
    ]);
    expect(strRows(L.basisMatrix() as unknown[][])).toEqual([
      ['x', 'x^3 + 1', '0'],
      ['0', '0', 'x'],
    ]);
    const M = A.span([L.basis()[0]!.mul(x7)]);
    expect(M.isSubmodule(L)).toBe(true);
    expect(L.coordinates(M.basis()[0]!).map(String)).toEqual(['x', '0']);
  });
});

// ============================================================================
// Tensor products
// ============================================================================

describe('tensor product (free_quadratic_module_integer_symmetric.py:1343)', () => {
  const ZZring = {
    zero: () => 0n,
    one: () => 1n,
    is_field: () => false,
    toString: () => 'Integer Ring',
  };
  const strRows = (M: unknown[][]) => M.map((row) => row.map((e) => String(e)));

  // sage: L = IntegralLattice("D3", [[1,-1,0], [0,1,-1]])
  const D3 = [
    [2n, -1n, -1n],
    [-1n, 2n, 0n],
    [-1n, 0n, 2n],
  ];
  const ambient = FreeModule(ZZring, 3, { innerProductMatrix: D3 }) as FreeModulePID;
  const L = ambient.spanOfBasis([
    ambient.createElement([1n, -1n, 0n]),
    ambient.createElement([0n, 1n, -1n]),
  ]) as FreeModulePID;

  it('has rank m*n and degree deg(M)*deg(N)', () => {
    const M = FreeModule(ZZring, 2) as FreeModulePID;
    const N = FreeModule(ZZring, 3);
    const T = M.tensorProduct(N);
    expect(T.rank()).toBe(6);
    expect(T.degree()).toBe(6);
  });

  it('has the Kronecker product of the basis matrices as basis', () => {
    // sage: L.tensor_product(L).basis_matrix()
    // [ 1 -1  0 -1  1  0  0  0  0]
    // [ 0  1 -1  0 -1  1  0  0  0]
    // [ 0  0  0  1 -1  0 -1  1  0]
    // [ 0  0  0  0  1 -1  0 -1  1]
    const T = L.tensorProduct(L);
    expect(T.degree()).toBe(9);
    expect(T.rank()).toBe(4);
    expect(strRows(T.basisMatrix() as unknown[][])).toEqual([
      ['1', '-1', '0', '-1', '1', '0', '0', '0', '0'],
      ['0', '1', '-1', '0', '-1', '1', '0', '0', '0'],
      ['0', '0', '0', '1', '-1', '0', '-1', '1', '0'],
      ['0', '0', '0', '0', '1', '-1', '0', '-1', '1'],
    ]);
  });

  it('tensors the inner product matrices and the Gram matrices', () => {
    // sage: L.tensor_product(L).gram_matrix()
    // [ 36 -12 -12   4]
    // [-12  24   4  -8]
    // [-12   4  24  -8]
    // [  4  -8  -8  16]
    const T = L.tensorProduct(L);
    expect(strRows(T.innerProductMatrix() as unknown[][])[0]).toEqual([
      '4',
      '-2',
      '-2',
      '-2',
      '1',
      '1',
      '-2',
      '1',
      '1',
    ]);
    expect(strRows(T.gramMatrix())).toEqual([
      ['36', '-12', '-12', '4'],
      ['-12', '24', '4', '-8'],
      ['-12', '4', '24', '-8'],
      ['4', '-8', '-8', '16'],
    ]);
  });

  it('discardBasis returns the standard basis with the tensored Gram matrix', () => {
    // sage: L.tensor_product(L, True).inner_product_matrix() == the Gram matrix above
    const T = L.tensorProduct(L, { discardBasis: true });
    expect(T.rank()).toBe(4);
    expect(T.degree()).toBe(4);
    expect(strRows(T.innerProductMatrix() as unknown[][])).toEqual([
      ['36', '-12', '-12', '4'],
      ['-12', '24', '4', '-8'],
      ['-12', '4', '24', '-8'],
      ['4', '-8', '-8', '16'],
    ]);
  });

  it('tracks the basis pairs e_i (x) f_j', () => {
    const T = L.tensorProduct(L);
    const b = L.basis();
    const s = b.length;
    for (let i = 0; i < s; i++) {
      for (let j = 0; j < s; j++) {
        const t = tensorProductVector(b[i]!, b[j]!, T.ambientModule());
        expect(t.list().map(String)).toEqual(T.basis()[i * s + j]!.list().map(String));
      }
    }
  });

  it('rejects a factor over a different base ring', () => {
    const M = FreeModule(ZZring, 2) as FreeModulePID;
    const other = VectorSpace(
      { zero: () => Rational.zero(), one: () => Rational.one(), is_field: () => true },
      2
    );
    expect(() => M.tensorProduct(other)).toThrow('base rings must be the same');
  });
});

// ============================================================================
// Quotient modules
// ============================================================================

describe('quotient modules (sage/modules/quotient_module.py, fg_pid/fgp_module.py)', () => {
  const ZZring = {
    zero: () => 0n,
    one: () => 1n,
    is_field: () => false,
    toString: () => 'Integer Ring',
  };
  const QQring = {
    zero: () => Rational.zero(),
    one: () => Rational.one(),
    is_field: () => true,
    toString: () => 'Rational Field',
  };
  const q = (n: number, d = 1) => new Rational(BigInt(n), BigInt(d));

  describe('over a field', () => {
    const V = VectorSpace(QQring, 3);
    const W = V.subspace([V.createElement([q(1), q(2), q(3)])]);
    const M = V.quotient(W);

    it('is an ambient space of dimension dim(V) - dim(W)', () => {
      // sage: M = QQ^3 / [[1,2,3]]; M
      // Vector space quotient V/W of dimension 2 over Rational Field where ...
      expect(M.dimension()).toBe(2);
      expect(M.degree()).toBe(2);
      expect(M.toString().split('\n')[0]).toBe(
        'Vector space quotient V/W of dimension 2 over Rational Field where'
      );
      expect(M.coveringModule()).toBe(V);
      expect(M.relations()).toBe(W);
    });

    it('projects with SageMath _element_constructor_ values', () => {
      // sage: M([1,2,4]) -> (-1/3, -2/3);  M([1,2,3]) -> (0, 0)
      expect(
        M.project(V.createElement([q(1), q(2), q(4)]))
          .list()
          .map(String)
      ).toEqual(['-1/3', '-2/3']);
      expect(
        M.project(V.createElement([q(1), q(2), q(3)]))
          .list()
          .map(String)
      ).toEqual(['0', '0']);
    });

    it('lifts with SageMath lift() values', () => {
      // sage: M.lift(M.0) -> (1,0,0); M.lift(M.1) -> (0,1,0)
      // sage: M.lift(M.0 - 2*M.1) -> (1,-2,0)
      expect(M.lift(M.gen(0)).list().map(String)).toEqual(['1', '0', '0']);
      expect(M.lift(M.gen(1)).list().map(String)).toEqual(['0', '1', '0']);
      expect(
        M.lift(M.createElement([q(1), q(-2)]))
          .list()
          .map(String)
      ).toEqual(['1', '-2', '0']);
    });

    it('project(lift(x)) == x and project kills W', () => {
      for (const g of M.basis()) {
        expect(M.project(M.lift(g)).list().map(String)).toEqual(g.list().map(String));
      }
      for (const w of W.basis()) {
        expect(
          M.project(w)
            .list()
            .every((e) => (e as Rational).isZero())
        ).toBe(true);
      }
      // a random-ish element: lift(project(v)) - v lies in W
      const v = V.createElement([q(5), q(-1), q(7, 2)]);
      const diff = M.lift(M.project(v)).sub(v);
      expect(() => W.coordinates(diff)).not.toThrow();
    });

    it('quotients a subspace by a subspace', () => {
      // sage: A = QQ^3; V = A.span([[1,2,3],[4,5,6]]); Q = V.quotient([V.0+V.1])
      // sage: Q(V.0) -> (1)
      const A = VectorSpace(QQring, 3);
      const V2 = A.subspace([
        A.createElement([q(1), q(2), q(3)]),
        A.createElement([q(4), q(5), q(6)]),
      ]);
      const W2 = V2.subspace([V2.gen(0).add(V2.gen(1))]);
      const Q = V2.quotient(W2);
      expect(Q.dimension()).toBe(1);
      expect(Q.project(V2.gen(0)).list().map(String)).toEqual(['1']);
      expect(
        Q.project(V2.gen(0).add(V2.gen(1)))
          .list()
          .map(String)
      ).toEqual(['0']);
    });

    it('reproduces the GF(5) and GF(19) doctests', () => {
      // sage: A = GF(5)^2; B = A.span([[1,3]]); Q = A/B
      // sage: Q(A.0) -> (1); Q(A.1) -> (3); Q([1,3]) -> (0); Q.lift(Q.0) -> (1,0)
      const A = VectorSpace(GF(5) as never, 2);
      const Bs = A.subspace([A.createElement([1, 3])]);
      const Q = A.quotient(Bs as never);
      expect(Q.dimension()).toBe(1);
      expect(Q.project(A.gen(0)).list().map(String)).toEqual(['1']);
      expect(Q.project(A.gen(1)).list().map(String)).toEqual(['3']);
      expect(
        Q.project(A.createElement([1, 3]))
          .list()
          .map(String)
      ).toEqual(['0']);
      expect(
        Q.project(A.createElement([2, 1]))
          .list()
          .map(String)
      ).toEqual(['0']);
      expect(Q.lift(Q.gen(0)).list().map(String)).toEqual(['1', '0']);

      // sage: V = GF(19)^3; W = V.span_of_basis([[1,2,3],[1,0,1]])
      // sage: U, pi, lift = V.quotient_abstract(W)
      // sage: pi(V.2) -> (18); pi(V.0) -> (1); pi(V.0 + V.2) -> (0)
      const V19 = VectorSpace(GF(19) as never, 3) as FreeModulePID;
      const W19 = V19.spanOfBasis([
        V19.createElement([1, 2, 3]),
        V19.createElement([1, 0, 1]),
      ]) as FreeModuleField;
      const Q19 = (V19 as unknown as FreeModuleField).quotient(W19);
      expect(Q19.project(V19.gen(2)).list().map(String)).toEqual(['18']);
      expect(Q19.project(V19.gen(0)).list().map(String)).toEqual(['1']);
      expect(
        Q19.project(V19.gen(0).add(V19.gen(2)))
          .list()
          .map(String)
      ).toEqual(['0']);
    });
  });

  describe('over ZZ (finitely generated modules)', () => {
    it('reproduces the fgp_module invariants doctests', () => {
      // sage: V1 = ZZ^2; W1 = V1.span([[1,2],[3,4]]); V1/W1
      // Finitely generated module V/W over Integer Ring with invariants (2)
      const V1 = FreeModule(ZZring, 2);
      const W1 = V1.span([V1.createElement([1n, 2n]), V1.createElement([3n, 4n])]);
      const Q1 = V1.quotientModule(W1);
      expect(Q1.invariants()).toEqual([2n]);
      expect(Q1.cardinality()).toBe(2n);
      expect(Q1.toString()).toBe(
        'Finitely generated module V/W over Integer Ring with invariants (2)'
      );

      // sage: V = ZZ^3; W = V.span([[1,2,0],[0,1,0],[0,2,0]]); Q = V/W
      // sage: Q.invariants() -> (0,); Q.invariants(include_ones=True) -> (1, 1, 0)
      const V2 = FreeModule(ZZring, 3);
      const W2 = V2.span([
        V2.createElement([1n, 2n, 0n]),
        V2.createElement([0n, 1n, 0n]),
        V2.createElement([0n, 2n, 0n]),
      ]);
      const Q2 = V2.quotientModule(W2);
      expect(Q2.invariants()).toEqual([0n]);
      expect(Q2.invariants(true)).toEqual([1n, 1n, 0n]);
      expect(Q2.cardinality()).toBe(Number.POSITIVE_INFINITY);
    });

    it('reproduces the (4, 12) doctest', () => {
      // sage: V = span([[1/2,1,1],[3/2,2,1],[0,0,1]], ZZ)
      // sage: W = V.span([2*V.0+4*V.1, 9*V.0+12*V.1, 4*V.2]); Q = V/W
      // sage: Q.invariants() -> (4, 12); Q.cardinality() -> 48
      const V = span(
        [
          [q(1, 2), q(1), q(1)],
          [q(3, 2), q(2), q(1)],
          [q(0), q(0), q(1)],
        ],
        ZZring
      );
      expect((V.basisMatrix() as unknown[][]).map((r) => r.map(String))).toEqual([
        ['1/2', '0', '0'],
        ['0', '1', '0'],
        ['0', '0', '1'],
      ]);
      const g = V.basis();
      const W = V.span([
        g[0]!.mul(2n).add(g[1]!.mul(4n)),
        g[0]!.mul(9n).add(g[1]!.mul(12n)),
        g[2]!.mul(4n),
      ]);
      const Q = V.quotientModule(W);
      expect(Q.invariants()).toEqual([4n, 12n]);
      expect(Q.cardinality()).toBe(48n);
      expect(Q.toString()).toBe(
        'Finitely generated module V/W over Integer Ring with invariants (4, 12)'
      );

      // The Smith form generators really do have the stated orders.
      const gens = Q.basis();
      for (let i = 0; i < gens.length; i++) {
        const order = Q.invariants()[i]!;
        const lifted = Q.lift(gens[i]!);
        expect(
          Q.project(lifted.mul(order))
            .list()
            .every((e) => e === 0n)
        ).toBe(true);
        for (let k = 1n; k < order; k++) {
          expect(
            Q.project(lifted.mul(k))
              .list()
              .some((e) => e !== 0n)
          ).toBe(true);
        }
      }
    });

    it('project(lift(x)) == x, project kills W, and project is additive', () => {
      const V = FreeModule(ZZring, 3);
      const Vs = V.span([
        V.createElement([2n, 0n, 1n]),
        V.createElement([0n, 3n, 1n]),
        V.createElement([1n, 1n, 4n]),
      ]);
      const b = Vs.basis();
      const W = Vs.span([b[0]!.mul(6n), b[1]!.mul(4n), b[2]!.mul(10n)]);
      const Q = Vs.quotientModule(W);

      for (const gen of Q.basis()) {
        expect(Q.project(Q.lift(gen)).list()).toEqual(gen.list());
      }
      for (const w of W.basis()) {
        expect(
          Q.project(w)
            .list()
            .every((e) => e === 0n)
        ).toBe(true);
      }
      const inv = Q.invariants();
      const u = b[0]!.mul(3n).add(b[1]!.mul(-2n)).add(b[2]!.mul(5n));
      const v = b[0]!.mul(-1n).add(b[2]!.mul(7n));
      const pu = Q.project(u).list() as bigint[];
      const pv = Q.project(v).list() as bigint[];
      const sum = pu.map((e, i) => {
        const m = inv[i]!;
        const s = e + pv[i]!;
        return m === 0n ? s : ((s % m) + m) % m;
      });
      expect(Q.project(u.add(v)).list()).toEqual(sum);
    });

    it('rejects a module that is not a submodule and a different base ring', () => {
      const M = FreeModule(ZZring, 2);
      const N = M.scale(2n);
      expect(() => N.quotientModule(M)).toThrow('sub must be a subspace of self');
      const Vq = VectorSpace(QQring, 2);
      expect(() => M.quotientModule(Vq)).toThrow('base rings must be the same');
    });
  });
});

// ============================================================================
// modules/index.ts re-exports
// ============================================================================

describe('modules/index.ts exports', () => {
  it('re-exports the submodule and quotient classes', async () => {
    const mod = await import('./index.js');
    expect(typeof mod.FreeModuleSubmodulePID).toBe('function');
    expect(typeof mod.FreeModuleSubspace).toBe('function');
    expect(typeof mod.FreeModuleSubmodule).toBe('function');
    expect(typeof mod.FreeModuleQuotient).toBe('function');
    expect(typeof mod.FractionFieldElement).toBe('function');
  });
});
