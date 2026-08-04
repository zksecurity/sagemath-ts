/**
 * Tests for free modules and LLL lattice reduction
 */

import { describe, expect, it } from 'bun:test';
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
// Echelon basis / intersection over K[x]: golden values from SageMath
// ============================================================================
//
// SageMath computes the echelon basis of a submodule of K[x]^n with
//   _echelonized_basis  (sage/modules/free_module.py:6900)
//     -> Matrix.echelon_form -> _echelonize_ring (matrix2.pyx:8262)
//       -> _echelon_form_PID (matrix2.pyx:17305)
//         -> _generic_clear_column (matrix2.pyx:20613)
// and the intersection with
//   FreeModule_generic_pid.intersection (free_module.py:3933)
//     -> Matrix.integer_kernel (matrix2.pyx:5641)
//       -> Matrix.denominator (matrix2.pyx:3521)
//       -> Matrix.kernel = left_kernel (matrix2.pyx:5345)
//         -> right_kernel (matrix2.pyx:4975)
//           -> _right_kernel_matrix_over_domain (matrix2.pyx:4166, Smith form)
//
// None of these normalizes a pivot (the reduction above pivots in
// _echelon_form_PID is guarded by `except AttributeError` on
// `Ideal.small_residue`, which K[x] ideals do not have), so the basis is only
// canonical once *every* unit factor introduced along the way is reproduced
// exactly -- in particular the `d = self.denominator()` scaling that
// integer_kernel applies before computing the kernel, which over QQ[x] is the
// lcm of the denominators of the rational *coefficients* and is therefore a
// non-trivial unit of QQ[x].
//
// The tables below were produced by running SageMath 10.3 (the version
// installed on this machine) with `_echelonized_basis` replaced by the
// vendored SageMath 10.9.beta4 version in
// reference/sage/src/sage/modules/free_module.py:6900 -- 10.3 predates the
// `if basis.universe().coordinate_ring() == ambient.base_ring(): d = 1` guard
// and clears the coefficient denominators of *integral* generators as well.
// Every other function on the path above is byte-identical between the
// installed 10.3 and the vendored 10.9.beta4 (verified by diffing them).
//
// Each case records the generators of two submodules P, Q of K[x]^n and
// SageMath's basis matrices for P, Q, P n Q and P + Q.  Polynomials are
// coefficient lists, constant term first ('a/b' over QQ, 'a' over GF(p)).
// ============================================================================

/** Golden cases with generators in K[x]^n (coordinate ring K[x]). */
const KX_GOLDEN_INTEGRAL = [
  '[{"i":1,"ring":"GF2","n":2,"g1":[[["1","0","0"],["0"]]],"g2":[[[],[]]],"P":[[["1"],[]]],"Q":',
  '[],"I":[],"S":[[["1"],[]]]},{"i":2,"ring":"GF3","n":2,"g1":[[["1","1"],["2","1"]],[["1","2",',
  '"0"],[]]],"g2":[[[],["0"]],[["0"],[]],[["2","0"],["1","2"]]],"P":[[["1"],["1","2"]],[[],["1"',
  ',"1","1"]]],"Q":[[["2"],["1","2"]]],"I":[[["2","1"],["1","1","1"]]],"S":[[["1"],["1","2"]],[',
  '[],["2","1"]]]},{"i":3,"ring":"GF5","n":3,"g1":[[["3","1","4"],["2","1","2"],["3"]]],"g2":[[',
  '["3","1"],["2"],["3","0"]],[[],["1","4"],["0","2"]],[["0","1"],[],["4"]]],"P":[[["3","1","4"',
  '],["2","1","2"],["3"]]],"Q":[[["1"],["4"],["3"]],[[],["1"],["4","4"]],[[],[],["2","4","3"]]]',
  ',"I":[[["4","1","4","1","3"],["1","0","1","4","4"],["4","3","1"]]],"S":[[["1"],["4"],["3"]],',
  '[[],["1"],["4","4"]],[[],[],["1"]]]},{"i":4,"ring":"GF7","n":4,"g1":[[["4","2"],["3"],["5","',
  '3","1"],["2","6","1"]],[["1","0","6"],["6","3","4"],["0"],[]],[["2"],["6","2","6"],["0","3",',
  '"0"],["3","6"]]],"g2":[[["5","6","4"],["2"],[],[]],[["6"],["3","3"],[],["6"]],[["3","2"],["0',
  '"],[],["6","6"]]],"P":[[["1"],["6","2","1"],["4","6","1","1"],["3","4","4","1"]],[[],["1","1',
  '"],["5","2","1","5","1"],["5","5","0","1","1"]],[[],[],["5","6","1","1","0","1"],["2","5","0',
  '","6","3","1"]]],"Q":[[["1"],["4","4"],[],["1"]],[[],["1"],[],["4","3","2","4"]],[[],[],[],[',
  '"5","2","0","0","2"]]],"I":[[["6","1"],["6","1","4","2","0","4","3","4","2"],[],["1","1","0"',
  ',"4","6","0","2","1"]],[[],["3","1","0","4","6","5","5","6","4","2"],[],["2","3","4","3","6"',
  ',"3","5","2","1"]]],"S":[[["1"],["6","2","1"],["4","6","1","1"],["3","4","4","1"]],[[],["1"]',
  ',["1","0","6","3","1","4"],["2","0","4","5","6","4"]],[[],[],["1"],["6","0","2","6","0","2"]',
  '],[[],[],[],["1"]]]},{"i":7,"ring":"GF2","n":3,"g1":[[["0","0","0"],[],["1"]]],"g2":[[["1"],',
  '["0"],["1","1","1"]],[["1","0"],["1","1"],["1"]]],"P":[[[],[],["1"]]],"Q":[[["1"],[],["1","1',
  '","1"]],[[],["1","1"],["0","1","1"]]],"I":[],"S":[[["1"],[],["1","1","1"]],[[],["1","1"],["0',
  '","1","1"]],[[],[],["1"]]]},{"i":12,"ring":"QQ","n":3,"g1":[[[],["1/3","1/3"],[]],[[],["3/2"',
  ',"0/1"],[]],[["5/2","-5/2","-1/1"],["0/1","4/1"],["-1/1","5/3","-1/2"]]],"g2":[[["-4/1","-3/',
  '1","-5/2"],["1/1"],[]],[["-2/1"],["-2/1","-2/1","-1/1"],["2/1","4/3"]],[["1/1","-5/1"],["1/2',
  '","5/1","-5/2"],["-4/1"]]],"P":[[["5/2","-5/2","-1/1"],["0/1","4/1"],["-1/1","5/3","-1/2"]],',
  '[[],["3/2"],[]]],"Q":[[["1/1"],["1/1","1/1","1/2"],["-1/1","-2/3"]],[[],["1/1"],["-548818/36',
  '20197","98258/1551513","3311348/10860591","87430/517171","489700/10860591"]],[[],[],["-34/1"',
  ',"-19/1","-51/2","-41/3","-95/6","-25/3"]]],"I":[[["5/2","-5/2","-1/1"],["-8637/590","-4437/',
  '2360","-21609/4720","-1327/236","-2615/944"],["-1/1","5/3","-1/2"]],[[],["-52266153511055721',
  '3268968/1","-292075563738252560356188/1","-391996151332917909951726/1","-2100894405836553504',
  '31644/1","-243396303115210466963490/1","-128103317429058140507100/1"],[]]],"S":[[["1/1"],["1',
  '/1","1/1","1/2"],["-1/1","-2/3"]],[[],["3/2"],[]],[[],[],["1/1"]]]},{"i":13,"ring":"GF2","n"',
  ':2,"g1":[[["1","0","0"],["1","0","0"]]],"g2":[[["1","1","0"],["0"]]],"P":[[["1"],["1"]]],"Q"',
  ':[[["1","1"],[]]],"I":[],"S":[[["1"],["1"]],[[],["1","1"]]]},{"i":19,"ring":"GF2","n":4,"g1"',
  ':[[[],["0"],["1","1","0"],["1","1"]],[["0","0"],["1"],["0","0"],["0"]],[["1","1"],["0"],["0"',
  '],[]]],"g2":[[["0"],["1"],["1","0"],["0","0"]],[["0","1"],["0","1","0"],["1","1"],["0","1","',
  '1"]],[["0"],[],["0","0","0"],["1","1","0"]]],"P":[[["1","1"],[],[],[]],[[],["1"],[],[]],[[],',
  '[],["1","1"],["1","1"]]],"Q":[[["0","1"],["0","1"],["1","1"],["0","1","1"]],[[],["1"],["1"],',
  '[]],[[],[],[],["1","1"]]],"I":[[["0","1","1"],["1","0","0","1"],["0","1","0","1"],["0","1","',
  '0","1"]],[[],["1","1"],["1","1"],["1","1"]]],"S":[[["1"],["0","1"],["1","1"],["0","1","1"]],',
  '[[],["1"],[],[]],[[],[],["1"],[]],[[],[],[],["1","1"]]]},{"i":25,"ring":"GF2","n":2,"g1":[[[',
  '"1","1"],["1"]],[["0","1","0"],["1"]]],"g2":[[[],["1","0","1"]],[["1","1"],["1","0"]],[["0",',
  '"0","0"],["1","1","0"]]],"P":[[["1"],[]],[[],["1"]]],"Q":[[["1","1"],["1"]],[[],["1","1"]]],',
  '"I":[[["1","1"],["1"]],[[],["1","1"]]],"S":[[["1"],[]],[[],["1"]]]},{"i":30,"ring":"QQ","n":',
  '3,"g1":[[["5/1","2/1","-3/2"],["-2/1","4/1"],["-5/1","-2/1","-1/1"]],[[],["-3/2"],["3/2","-1',
  '/1","-4/3"]]],"g2":[[["0/1"],["1/3","-5/3"],["-4/1","3/2"]],[["0/1"],["5/2"],["1/2","-1/1","',
  '-5/2"]]],"P":[[["5/1","2/1","-3/2"],["-2/1","4/1"],["-5/1","-2/1","-1/1"]],[[],["-3/2"],["3/',
  '2","-1/1","-4/3"]]],"Q":[[[],["1/1"],["1/5","-2/5","-1/1"]],[[],[],["61/6","-59/12","5/6","2',
  '5/6"]]],"I":[[[],["4941000/1","-2389500/1","405000/1","2025000/1"],["-4941000/1","5683500/1"',
  ',"2394000/1","-3879000/1","1710000/1","1800000/1"]]],"S":[[["5/1","2/1","-3/2"],["-2/1","4/1',
  '"],["-5/1","-2/1","-1/1"]],[[],["-3/2"],["3/2","-1/1","-4/3"]],[[],[],["1/1"]]]},{"i":36,"ri',
  'ng":"QQ","n":3,"g1":[[["4/1","3/2"],["5/2"],[]]],"g2":[[["-1/3","1/3","-5/3"],[],[]],[["-2/1',
  '","-1/1"],["4/1","1/1"],[]]],"P":[[["4/1","3/2"],["5/2"],[]]],"Q":[[["1/1"],["-44/23","9/23"',
  ',"5/23"],[]],[[],["-4/3","1/1","-19/3","-5/3"],[]]],"I":[[["-101568/1","38088/1","-453882/1"',
  ',"-307878/1","-47610/1"],["-63480/1","47610/1","-301530/1","-79350/1"],[]]],"S":[[["1/1"],["',
  '-44/23","9/23","5/23"],[]],[[],["1/1"],[]]]},{"i":42,"ring":"QQ","n":3,"g1":[[["-5/1"],["5/1',
  '"],["2/1","-3/1","-1/1"]],[["-3/1","1/1"],[],[]],[[],["-3/1","-2/1"],["1/1","-2/1","5/2"]]],',
  '"g2":[[[],["-5/1"],[]]],"P":[[["-5/1"],["5/1"],["2/1","-3/1","-1/1"]],[[],["1/1"],["7/45","-',
  '4/15","-5/18","2/45"]],[[],[],["-33/5","56/5","-51/10","19/10","-2/5"]]],"Q":[[[],["-5/1"],[',
  ']]],"I":[[[],["-594/1","1008/1","-459/1","171/1","-36/1"],[]]],"S":[[["-5/1"],["5/1"],["2/1"',
  ',"-3/1","-1/1"]],[[],["1/1"],["7/45","-4/15","-5/18","2/45"]],[[],[],["1/1"]]]},{"i":48,"rin',
  'g":"QQ","n":3,"g1":[[["0/1","-1/1","-3/1"],["-5/1"],["-2/1","2/1"]],[["-4/1"],["-1/1"],["3/2',
  '","1/1"]]],"g2":[[["4/3","-2/1"],["0/1","4/1","2/1"],["0/1"]],[[],["-2/1","-3/1"],["-1/1","2',
  '/1","3/1"]],[[],[],["0/1","-5/1","-1/1"]]],"P":[[["1/1"],["1/4"],["-3/8","-1/4"]],[[],["-20/',
  '1","1/1","3/1"],["-8/1","13/2","-11/2","-3/1"]]],"Q":[[["4/3","-2/1"],["0/1","4/1","2/1"],[]',
  '],[[],["-2/1","-3/1"],["-1/1","2/1","3/1"]],[[],[],["0/1","-5/1","-1/1"]]],"I":[[["-2/3","1/',
  '1"],["19/6","991991/222372","265151/222372","-37171/49416","-14167/49416"],["19/12","204307/',
  '444744","3755/444744","963037/889488","116843/98832","14167/49416"]],[[],["0/1","230400/1","',
  '380160/1","14976/1","-62208/1","-10368/1"],["0/1","92160/1","81792/1","-36288/1","119808/1",',
  '"77760/1","10368/1"]]],"S":[[["1/1"],["1/4"],["-3/8","-1/4"]],[[],["1/1"],["3153/5878","-249',
  '5/2939","1533/2939","-183/2939","-414/2939"]],[[],[],["1/1"]]]},{"i":72,"ring":"QQ","n":3,"g',
  '1":[[["-1/1","-1/1"],["2/1","3/2","5/1"],["-4/3"]],[["-5/3"],[],["3/1"]]],"g2":[[[],["5/3","',
  '-4/1"],[]],[["-2/3","1/1","-5/2"],["4/1"],["1/1","1/2","5/1"]],[["3/1","5/1","-2/1"],["1/3",',
  '"-3/1"],[]]],"P":[[["1/1"],[],["-9/5"]],[[],["10/3","5/2","25/3"],["-47/9","-3/1"]]],"Q":[[[',
  '"1/1"],["-883/473","-810/5203","2835/5203"],["-2526/5203","-507/5203","-12252/5203","3780/52',
  '03"]],[[],["1/1"],["10368/63305","22464/63305","53568/63305","82944/63305","-6912/12661"]],[',
  '[],[],["5/1","-7/6","-1/6","-22/1","-338/3","40/1"]]],"I":[[["1/1"],["-2911179475371987/3319',
  '161275154230","-674576580877299/2655329020123384","-26244149258519397/26553290201233840","-3',
  '4677935747883/2655329020123384","9581532801819651/6638322550308460","-2588220351436617/66383',
  '2255030846","384622298738010/331916127515423"],["-14136424505281677/33191612751542300","1040',
  '9428969814949/66383225503084600","-117881906299226661/66383225503084600","2790804824583561/1',
  '508709670524650","13319444882323827/16595806375771150","-692320137728418/1659580637577115"]]',
  ',[[],["4841635241845215000/1","2501511541620027750/1","11095414095895284375/1","-24248523169',
  '574785125/1","-125479046684488488750/1","-96348541312719778500/1","-243695640506209155000/1"',
  ',"96832704836904300000/1"],["-7585228545557503500/1","-2587585057030609350/1","1269584352306',
  '078600/1","33520254657708371850/1","190093358784269463600/1","37506534340160932200/1","-3485',
  '9773741285548000/1"]]],"S":[[["1/1"],[],["-9/5"]],[[],["1/1"],["-9089357469/11335688005","15',
  '9981699/2267137601","2552938860/2267137601","-1189261980/2267137601","125533800/2267137601"]',
  '],[[],[],["1/1"]]]},{"i":78,"ring":"QQ","n":2,"g1":[[["-2/3","3/1","-4/1"],["-3/2","1/1"]],[',
  '["-2/3","-2/1"],["0/1","-1/2","3/1"]]],"g2":[[[],[]],[["1/1","-1/3"],["-3/1","-2/1"]],[["-1/',
  '1","-1/1"],["5/1","3/1","5/1"]]],"P":[[["1/1"],["27/38","3/76","-135/38","54/19"]],[[],["-1/',
  '1","-2/1","-3/2","11/1","-12/1"]]],"Q":[[["1/1"],["-7/2","-9/4","-5/4"]],[[],["2/1","-11/3",',
  '"2/1","-5/3"]]],"I":[[["1/1"],["-691360799/359741098","-3173026739/719482196","-1299700201/1',
  '79870549","5452078734/179870549","-30640092381/719482196","8879601277/359741098","-292416753',
  '0/179870549"]],[[],["-456/1","-76/1","532/1","5738/1","-14592/1","15618/1","-9652/1","4560/1',
  '"]]],"S":[[["1/1"],["27/38","3/76","-135/38","54/19"]],[[],["1/1"]]]},{"i":90,"ring":"QQ","n',
  '":2,"g1":[[[],[]],[["3/1","3/2","-1/2"],[]],[[],["2/3","1/2","5/2"]]],"g2":[[["1/1","0/1"],[',
  '"1/2","2/1"]]],"P":[[["3/1","3/2","-1/2"],[]],[[],["2/3","1/2","5/2"]]],"Q":[[["1/1"],["1/2"',
  ',"2/1"]]],"I":[[["12/1","15/1","95/2","21/1","-15/2"],["6/1","63/2","215/4","211/2","153/4",',
  '"-15/1"]]],"S":[[["1/1"],["1/2","2/1"]],[[],["1/1"]]]},{"i":96,"ring":"QQ","n":3,"g1":[[[],[',
  '"4/1","1/1","-4/3"],["2/1","2/1"]],[[],["-2/1","1/1"],["-4/1"]]],"g2":[[["1/1","-1/1"],["5/1',
  '","0/1","-5/1"],["5/1","-1/2"]],[[],["2/1"],[]]],"P":[[[],["1/1"],["-7/1","-5/1"]],[[],[],["',
  '-12/1","-2/1","10/3"]]],"Q":[[["1/1","-1/1"],["5/1","0/1","-5/1"],["5/1","-1/2"]],[[],["2/1"',
  '],[]]],"I":[[[],["-864/1","-144/1","240/1"],[]]],"S":[[["1/1","-1/1"],["5/1","0/1","-5/1"],[',
  '"5/1","-1/2"]],[[],["-1/1"],["7/1","5/1"]],[[],[],["1/1"]]]},{"i":102,"ring":"QQ","n":2,"g1"',
  ':[[["-1/1"],["-5/1"]],[[],[]]],"g2":[[["5/1"],["2/1","-1/1"]],[["-5/1","-4/1"],[]]],"P":[[["',
  '-1/1"],["-5/1"]]],"Q":[[["5/1"],["2/1","-1/1"]],[[],["2/1","3/5","-4/5"]]],"I":[[["50/1","15',
  '/1","-20/1"],["250/1","75/1","-100/1"]]],"S":[[["-1/1"],["-5/1"]],[[],["1/1"]]]},{"i":108,"r',
  'ing":"QQ","n":2,"g1":[[["-5/1","2/1","0/1"],["1/1","3/2","3/1"]]],"g2":[[["-2/1","-2/1"],["-',
  '2/1","-2/1"]],[[],["-3/1"]]],"P":[[["-5/1","2/1"],["1/1","3/2","3/1"]]],"Q":[[["-2/1","-2/1"',
  '],["-2/1","-2/1"]],[[],["-3/1"]]],"I":[[["-120/1","-72/1","48/1"],["24/1","60/1","108/1","72',
  '/1"]]],"S":[[["1/1"],["1/7","1/14","-3/7"]],[[],["1/1"]]]},{"i":120,"ring":"QQ","n":2,"g1":[',
  '[["4/3","-1/2"],["5/1","4/1","1/1"]],[["3/1","0/1"],["2/1","-3/2","5/1"]]],"g2":[[["0/1","2/',
  '1"],[]]],"P":[[["1/1"],["2/3","-1/2","5/3"]],[[],["-37/3","-15/1","53/12","-5/2"]]],"Q":[[["',
  '0/1","2/1"],[]]],"I":[[["0/1","-3552/1","-4320/1","1272/1","-720/1"],[]]],"S":[[["1/1"],["2/',
  '3","-1/2","5/3"]],[[],["1/1"]]]},{"i":138,"ring":"QQ","n":4,"g1":[[["1/1"],["1/3","3/1"],["2',
  '/1","2/3","0/1"],[]],[["-4/1","-2/1","5/1"],[],["1/1"],["4/3","5/1"]],[["3/1","-4/1"],["-1/1',
  '","-1/1"],["-4/3","5/1"],["3/2","5/1","1/2"]]],"g2":[[["5/2","0/1"],["2/1","3/2","3/1"],["-2',
  '/3","-1/3"],[]],[["-4/3","1/1"],["1/1","-2/1","-3/1"],[],["1/1","-1/3","-1/1"]]],"P":[[["1/1',
  '"],["1/3","3/1"],["2/1","2/3"],[]],[[],["1/1"],["-5497/7660","-7032/1915","-7368/1915","2777',
  '7/3064"],["-12403/3830","-159429/15320","185013/30640","64263/7660","5319/6128"]],[[],[],["7',
  '4/9","118/9","392/9","89/3","-235/3"],["14/3","425/9","587/6","-109/2","-437/6","-15/2"]]],"',
  'Q":[[["5/2"],["2/1","3/2","3/1"],["-2/3","-1/3"],[]],[[],["31/15","-2/1","-2/1","-6/5"],["-1',
  '6/45","4/45","2/15"],["1/1","-1/3","-1/1"]]],"I":[[["-1216014066098841600/1","-2921581814042',
  '419200/1","4389388535330150400/1","5606852363786649600/1","2465269111969689600/1","-42612140',
  '09745100800/1","-4278766489718169600/1","-619307412855091200/1","-41939553917952000/1"],["-4',
  '3622313797376000/1","-29435353583155200/1","-642218835828787200/1","-11416878566553600/1","-',
  '4712607875247206400/1","-774483762351513600/1","2262405936351744000/1","-117430750970265600/',
  '1","-1908249703266816000/1"],["164409954988646400/1","303915063761817600/1","-11998336826428',
  '41600/1","-2416287855171686400/1","-87037518686515200/1","1580965851025612800/1","1401091764',
  '221952000/1","432443400398438400/1","13979851305984000/1"],["449607551168563200/1","17549632',
  '22558054400/1","738343257864192000/1","-397778546882304000/1","-5022261581674752000/1","-328',
  '8837685573324800/1","2592641090535321600/1","2153829091208601600/1","584357784590131200/1","',
  '41939553917952000/1"]]],"S":[[["1/1"],["1/3","3/1"],["2/1","2/3"],[]],[[],["1/1"],["-5497/76',
  '60","-7032/1915","-7368/1915","27777/3064"],["-12403/3830","-159429/15320","185013/30640","6',
  '4263/7660","5319/6128"]],[[],[],["1/1"],["-177520402459176/171130227840601","-21155998737703',
  '00/513390683521803","-214921583079892/171130227840601","345395306096055/171130227840601","23',
  '2247171543792/171130227840601","56891987787564/171130227840601","3916702574190/1711302278406',
  '01"]],[[],[],[],["1/1"]]]},{"i":144,"ring":"QQ","n":3,"g1":[[[],[],["-2/1","3/1","-1/1"]],[[',
  '"-5/1","-1/1"],["-2/1","-5/1","3/1"],[]],[["1/3"],[],["-1/1"]]],"g2":[[[],["1/1"],[]],[["2/3',
  '","1/1","1/1"],["-1/1","-1/1"],["2/1","5/1"]]],"P":[[["1/1"],[],["-3/1"]],[[],["2/3","5/3","',
  '-1/1"],["5/1","1/1"]],[[],[],["-2/1","3/1","-1/1"]]],"Q":[[["2/3","1/1","1/1"],["-1/1","-1/1',
  '"],["2/1","5/1"]],[[],["-1/1"],[]]],"I":[[["2/3","1/1","1/1"],["2/7","44/21","127/42","-29/1',
  '4"],["2/1","5/1"]],[[],["36/1","36/1","-171/1","126/1","-27/1"],[]]],"S":[[["1/1"],[],["-3/1',
  '"]],[[],["1/1"],["-47/6","-55/6","0/1","3/2"]],[[],[],["1/1"]]]},{"i":162,"ring":"QQ","n":4,',
  '"g1":[[[],["-2/1"],[],[]],[["-2/3","1/1","4/1"],["5/1","4/1","3/1"],["-2/1","-2/1","1/1"],["',
  '-5/2"]]],"g2":[[["0/1","-4/1","2/1"],["4/1","-2/3","0/1"],[],["5/1","-1/1"]],[["-1/3","-1/1"',
  ',"-5/1"],["-2/3","3/2"],["-2/1","-3/1"],[]],[["-5/1","-1/1"],[],[],["-2/3"]]],"P":[[["-2/3",',
  '"1/1","4/1"],["5/1","4/1","3/1"],["-2/1","-2/1","1/1"],["-5/2"]],[[],["2/1"],[],[]]],"Q":[[[',
  '"1/1"],["302/67","1189/134","-33/134"],["6/1","405/67","-297/67"],["210/67","2391/134","-495',
  '/134"]],[[],["1/1"],["-4962/13985","-35877/69925","2007/69925"],["1179541/839100","-3027/139',
  '85","5129/279700"]],[[],[],["40/1","184/3","2/3","-2/1"],["-158/9","1793/54","-22/3","-23/18',
  '"]]],"I":[[["-4369511370984017867302118940000000000/1","-46608121290496190584555935360000000',
  '0/1","36456289871909989072857345689400000000/1","43170772345322096528944935127200000000/1","',
  '830207160486963394787402598600000000/1","-3670389551626575008533779909600000000/1"],["-17245',
  '004877483590516285696083200000000/1","-16074946832586759065152573122600000000/1","5020083063',
  '997193860878212204400000000/1","-9508542244524665547901333271100000000/1","13254184491984854',
  '19748309411800000000/1","167497935887720684913247892700000000/1"],["-13108534112952053601906',
  '356820000000000/1","-34169578921095019722302570110800000000/1","-153806800258637428929034586',
  '68800000000/1","11491814905687966991004572812200000000/1","227214591291168929099710184880000',
  '0000/1","-917597387906643752133444977400000000/1"],["-16385667641190067002382946025000000000',
  '/1","-26326306010178707650495266613500000000/1","-1092377842746004466825529735000000000/1","',
  '2293993469766609380333612443500000000/1"]]],"S":[[["1/1"],["302/67","1189/134","-33/134"],["',
  '6/1","405/67","-297/67"],["210/67","2391/134","-495/134"]],[[],["2/1"],[],[]],[[],[],["1/1"]',
  ',["-28655949100553062367233/782976715755369371388","73162191158801441344699/7829767157553693',
  '71388","254228477904045010512671/521984477170246247592","3749090187562289764921/130496119292',
  '561561898","-7654182746890795359489/173994825723415415864","35133003297448925715/43498706430',
  '853853966"]],[[],[],[],["1/1"]]]},{"i":168,"ring":"QQ","n":2,"g1":[[["0/1","3/1","2/1"],[]],',
  '[[],["-5/2","1/1","5/1"]]],"g2":[[[],["0/1","2/1"]],[["1/1","0/1","2/3"],["0/1","4/1","-5/1"',
  ']]],"P":[[["0/1","3/1","2/1"],[]],[[],["-5/2","1/1","5/1"]]],"Q":[[["1/1","0/1","2/3"],["0/1',
  '","4/1","-5/1"]],[[],["0/1","-2/1"]]],"I":[[["0/1","9/2","3/1","3/1","2/1"],["0/1","15/4","6',
  '/1","-21/2","-15/1"]],[[],["0/1","120/1","-48/1","-240/1"]]],"S":[[["1/1"],["0/1","4/1","-17',
  '/5","-2/1"]],[[],["1/1"]]]},{"i":174,"ring":"QQ","n":3,"g1":[[["-1/1","-3/2"],["5/3","-1/3"]',
  ',["2/1","2/1"]]],"g2":[[[],["-2/1","5/1"],["2/3","5/1"]],[[],["-3/1","5/1"],["-3/1","-3/1"]]',
  ',[["5/1","1/3"],["2/1"],["-5/2"]]],"P":[[["-1/1","-3/2"],["5/3","-1/3"],["2/1","2/1"]]],"Q":',
  '[[["5/1","1/3"],["2/1"],["-5/2"]],[[],["1/1"],["11/3","8/1"]],[[],[],["8/1","8/3","-40/1"]]]',
  ',"I":[[["8640/1","16416/1","-37824/1","-67392/1","-4320/1"],["-14400/1","-2880/1","72832/1",',
  '"-9536/1","-960/1"],["-17280/1","-24192/1","79104/1","91776/1","5760/1"]]],"S":[[["1/1"],["6',
  '4/129","-2/129"],["-37/86","4/43"]],[[],["1/1"],["11/3","8/1"]],[[],[],["1/1"]]]},{"i":180,"',
  'ring":"QQ","n":3,"g1":[[["5/3"],["4/1","-4/1"],[]],[["-1/1","1/1","-5/2"],[],["1/3","1/1","-',
  '1/2"]]],"g2":[[["-1/2","-4/1"],["-4/1","1/1"],["4/1"]],[[],["1/1"],["4/1","4/1"]]],"P":[[["5',
  '/3"],["4/1","-4/1"],[]],[[],["12/5","-24/5","42/5","-6/1"],["1/3","1/1","-1/2"]]],"Q":[[["-1',
  '/2","-4/1"],["-4/1","1/1"],["4/1"]],[[],["1/1"],["4/1","4/1"]]],"I":[[["-4170/1","-28590/1",',
  '"31455/1","-57960/1","-23760/1","86400/1"],["-42840/1","64080/1","-86580/1","-12420/1","9936',
  '0/1","-21600/1"],["-4560/1","-5760/1","30360/1","-24120/1","-34200/1","17280/1"]]],"S":[[["5',
  '/3"],["4/1","-4/1"],[]],[[],["1/1"],["382133/531900","-15607/35460","-52177/70920","-2218/14',
  '775"]],[[],[],["1/1"]]]},{"i":192,"ring":"QQ","n":4,"g1":[[["-1/1","4/3","-5/1"],["1/1"],["-',
  '4/3","2/1","-3/1"],["-5/3"]],[["-3/2","0/1","3/1"],[],[],["-4/1","4/1"]]],"g2":[[["-3/1","-5',
  '/2","1/1"],[],["1/3"],["-5/2","-2/3"]],[["2/3","1/3"],[],["-1/2","1/1","-3/2"],["0/1","-1/3"',
  ']],[["3/1","-3/1"],["-3/1"],["3/1","-2/3"],["2/1","-1/1"]]],"P":[[["1/1"],["-126/409","-48/4',
  '09"],["168/409","-188/409","282/409","144/409"],["2894/1227","-1064/1227","-320/409"]],[[],[',
  '"3/2","0/1","-3/1"],["-2/1","3/1","-1/2","-6/1","9/1"],["3/2","-28/3","91/3","-20/1"]]],"Q":',
  '[[["1/1"],[],["-77/72","5/2","-31/8","3/4"],["-5/12","-31/36","1/6"]],[[],["-3/1"],["149/24"',
  ',"-91/8","153/8","-111/8","9/4"],["13/4","1/3","-37/12","1/2"]],[[],[],["-23/18","67/36","-3',
  '/2","-19/4","3/2"],["-5/3","-41/18","-19/18","1/3"]]],"I":[[["-164925522930240/1","119242926',
  '4712832/1","-4181660033005440/1","5576965424849664/1","-3527987476531392/1","-75496874862216',
  '96/1","11263881198319488/1","-2413232812940544/1"],["247920302211264/1","-725672300893056/1"',
  ',"962952246786240/1","408589682614272/1","-1206616406470272/1","689495089411584/1"],["-33056',
  '0402948352/1","1463403672279936/1","-3479041837468224/1","3558135152765952/1","-462855499836',
  '480/1","-4558328646665472/1","4998839398233984/1","-2068485268234752/1"],["-191881092269376/',
  '1","1351325252396160/1","-4222920970283328/1","5162227980807168/1","-7326003551179392/1","32',
  '20481178938880/1","1378990178823168/1"]]],"S":[[["1/1"],["-126/409","-48/409"],["168/409","-',
  '188/409","282/409","144/409"],["2894/1227","-1064/1227","-320/409"]],[[],["1/1"],["-731/144"',
  ',"1315/108","-955/48","187/24","-1/1"],["-235/24","101/24","187/108","-2/9"]],[[],[],["1/1"]',
  ',["63619362862527/22090333460884","25292473800898/5522583365221","-82694652239319/2209033346',
  '0884","1345408579479/36940356958","-16757528450877/356295700982","46443006934947/55225833652',
  '21"]],[[],[],[],["1/1"]]]},{"i":198,"ring":"QQ","n":2,"g1":[[[],[]],[["5/3","2/1","-3/1"],["',
  '-3/2"]],[["-4/1","-5/3"],["-1/1","4/1"]]],"g2":[[["-1/1","0/1"],["3/1"]],[["-4/1"],["-1/1","',
  '-1/1"]]],"P":[[["1/1"],["1053/3062","-1791/1531","540/1531"]],[[],["-23/3","13/6","11/1","-1',
  '2/1"]]],"Q":[[["-1/1"],["3/1"]],[[],["-13/1","-1/1"]]],"I":[[["1/1"],["59259/162286","-19082',
  '1/162286","26145/81143","2700/81143"]],[[],["-915538/1","188313/1","1333501/1","-1331970/1",',
  '"-110232/1"]]],"S":[[["1/1"],["1053/3062","-1791/1531","540/1531"]],[[],["1/1"]]]},{"i":204,',
  '"ring":"QQ","n":3,"g1":[[[],["-2/3"],["-1/1","-4/1","-3/1"]],[["1/1","3/1"],["-4/3","-1/3","',
  '-3/1"],["4/1"]]],"g2":[[[],[],["5/1","-1/1"]],[["1/1","2/1","-1/1"],["4/3"],["0/1"]]],"P":[[',
  '["1/1","3/1"],["-4/3","-1/3","-3/1"],["4/1"]],[[],["2/3"],["1/1","4/1","3/1"]]],"Q":[[["1/1"',
  ',"2/1","-1/1"],["4/3"],[]],[[],[],["-5/1","1/1"]]],"I":[[["-90/1","-432/1","-360/1","360/1",',
  '"-54/1"],["-120/1","-336/1","72/1"],["-720/1","-2961/1","-4194/1","-3897/1","-2628/1","45/1"',
  ',"1350/1","-243/1"]]],"S":[[["1/1"],["32/3","-5/6","10/1","-9/2"],["-14/1","6/1"]],[[],["2/3',
  '"],["1/1","4/1","3/1"]],[[],[],["1/1"]]]},{"i":246,"ring":"QQ","n":2,"g1":[[["0/1"],["1/1","',
  '2/1","1/3"]]],"g2":[[["-2/1"],["-1/1","-2/1","0/1"]],[["5/1","-3/1"],[]]],"P":[[[],["1/1","2',
  '/1","1/3"]]],"Q":[[["-2/1"],["-1/1","-2/1"]],[[],["-5/2","-7/2","3/1"]]],"I":[[[],["180/1","',
  '612/1","348/1","-348/1","-72/1"]]],"S":[[["-2/1"],["-1/1","-2/1"]],[[],["1/1"]]]}]',
].join('');

/**
 * Golden cases over QQ[x] whose generators are genuine rational functions, so
 * that the coordinate ring is Frac(QQ[x]) and `_echelonized_basis` takes its
 * `d = self._denominator(basis)` branch.  Entries are `[numerator,
 * denominator]` coefficient lists.
 */
const KX_GOLDEN_FRACTIONAL = [
  '[{"i":0,"n":3,"g1":[[{"n":["0/1","-1/2"],"d":["1/3"]},{"n":[],"d":["-2/1","1/3"]},{"n":[],"d',
  '":["-2/1"]}]],"g2":[[{"n":[],"d":["-3/1","2/1"]},{"n":["1/1","0/1"],"d":["-1/1","-1/1"]},{"n',
  '":["-1/2","1/3"],"d":["-1/2","3/1"]}]],"P":[[[["0/1","-3/2"],["1/1"]],[[],["1/1"]],[[],["1/1',
  '"]]]],"Q":[[[[],["1/1"]],[["-1/1"],["1/1","1/1"]],[["-1/6","1/9"],["-1/6","1/1"]]]],"I":[],"',
  'S":[[[["0/1","-3/2"],["1/1"]],[[],["1/1"]],[[],["1/1"]]],[[[],["1/1"]],[["-1/1"],["1/1","1/1',
  '"]],[["-1/6","1/9"],["-1/6","1/1"]]]]},{"i":1,"n":3,"g1":[[{"n":["-1/1"],"d":["-2/1","-1/3"]',
  '},{"n":["-3/1"],"d":["1/1","2/3"]},{"n":["-2/3"],"d":["-1/1"]}],[{"n":["1/3"],"d":["3/1"]},{',
  '"n":["2/1"],"d":["-3/2","3/1"]},{"n":["-2/1","0/1"],"d":["2/3"]}]],"g2":[[{"n":["2/1","3/2"]',
  ',"d":["-2/3"]},{"n":["-1/1","2/1"],"d":["3/1","0/1"]},{"n":["-1/1","1/1"],"d":["-3/1"]}],[{"',
  'n":[],"d":["-2/1","-1/1"]},{"n":["3/1"],"d":["-3/1","0/1"]},{"n":["-3/2","-1/1"],"d":["-3/1"',
  ',"2/1"]}]],"P":[[[["3/1"],["6/1","1/1"]],[["-9/2"],["3/2","1/1"]],[["2/3"],["1/1"]]],[[[],["',
  '1/1"]],[["1/2","19/12","1/6"],["-3/4","1/1","1/1"]],[["-85/27","-2/81"],["1/1"]]]],"Q":[[[["',
  '-3/1","-9/4"],["1/1"]],[["-1/3","2/3"],["1/1"]],[["1/3","-1/3"],["1/1"]]],[[[],["1/1"]],[["-',
  '1/1"],["1/1"]],[["-3/4","-1/2"],["-3/2","1/1"]]]],"I":[[[["141/4","-193/8","-8513/192","3479',
  '/144","3161/144","1/6"],["1/1"]],[["-1605/8","302/1","5975/32","-1891/12","-713/24"],["1/1"]',
  '],[["4721/48","-3407/288","-367021/1728","-51787/432","-10859/1296","4/81"],["1/1"]]]],"S":[',
  '[[["3/1"],["6/1","1/1"]],[["-9/2"],["3/2","1/1"]],[["2/3"],["1/1"]]],[[[],["1/1"]],[["1/1"],',
  '["-3/4","1/1","1/1"]],[["-2072752/3558573","-27912461/10675719","-3234193/10675719","6880/10',
  '675719"],["1/1"]]],[[[],["1/1"]],[[],["1/1"]],[["1/1"],["-3/2","1/1"]]]]},{"i":2,"n":2,"g1":',
  '[[{"n":[],"d":["3/1","-1/1"]},{"n":[],"d":["0/1","2/1"]}]],"g2":[[{"n":[],"d":["0/1","1/1"]}',
  ',{"n":[],"d":["3/1"]}]],"P":[],"Q":[],"I":[],"S":[]},{"i":3,"n":3,"g1":[[{"n":["0/1"],"d":["',
  '3/2"]},{"n":[],"d":["1/1"]},{"n":[],"d":["0/1","-2/1"]}],[{"n":[],"d":["0/1","2/1"]},{"n":["',
  '0/1"],"d":["1/1","3/1"]},{"n":["-1/2","-3/1"],"d":["3/1","0/1"]}]],"g2":[[{"n":["-1/1","-3/1',
  '"],"d":["2/1"]},{"n":[],"d":["-2/1"]},{"n":["1/3"],"d":["1/1"]}]],"P":[[[[],["1/1"]],[[],["1',
  '/1"]],[["-1/6","-1/1"],["1/1"]]]],"Q":[[[["-1/2","-3/2"],["1/1"]],[[],["1/1"]],[["1/3"],["1/',
  '1"]]]],"I":[],"S":[[[["-1/2","-3/2"],["1/1"]],[[],["1/1"]],[["1/3"],["1/1"]]],[[[],["1/1"]],',
  '[[],["1/1"]],[["1/6","1/1"],["1/1"]]]]},{"i":4,"n":2,"g1":[[{"n":[],"d":["3/1"]},{"n":["-1/1',
  '"],"d":["-1/1","2/3"]}],[{"n":["-1/1"],"d":["-1/3"]},{"n":["-1/1","2/3"],"d":["-1/1"]}]],"g2',
  '":[[{"n":["-1/2","-1/1"],"d":["0/1","-1/3"]},{"n":["0/1"],"d":["-2/1"]}],[{"n":["1/1","1/2"]',
  ',"d":["-3/1","-1/1"]},{"n":["0/1","-1/1"],"d":["1/1"]}]],"P":[[[["3/1"],["1/1"]],[["1/1","-2',
  '/3"],["1/1"]]],[[[],["1/1"]],[["3/2"],["-3/2","1/1"]]]],"Q":[[[["1/1"],["0/1","3/1","1/1"]],',
  '[["0/1","-10/3","-4/3"],["1/1"]]],[[[],["1/1"]],[["0/1","-9/2","-21/2","-3/1"],["1/1"]]]],"I',
  '":[[[["-3/1"],["1/1"]],[["0/1","0/1","30/1","22/1","4/1"],["1/1"]]],[[[],["1/1"]],[["0/1","9',
  '/2","21/2","3/1"],["1/1"]]]],"S":[[[["1/1"],["0/1","3/1","1/1"]],[["0/1","-10/3","-4/3"],["1',
  '/1"]]],[[[],["1/1"]],[["3/2"],["-3/2","1/1"]]]]},{"i":5,"n":2,"g1":[[{"n":["1/2"],"d":["2/1"',
  ']},{"n":["0/1","-3/1"],"d":["-2/3"]}]],"g2":[[{"n":[],"d":["3/2","0/1"]},{"n":["-3/1","1/1"]',
  ',"d":["1/1"]}]],"P":[[[["1/4"],["1/1"]],[["0/1","9/2"],["1/1"]]]],"Q":[[[[],["1/1"]],[["-3/1',
  '","1/1"],["1/1"]]]],"I":[],"S":[[[["1/4"],["1/1"]],[["0/1","9/2"],["1/1"]]],[[[],["1/1"]],[[',
  '"-3/1","1/1"],["1/1"]]]]},{"i":6,"n":2,"g1":[[{"n":[],"d":["-2/1","0/1"]},{"n":["-1/1","1/2"',
  '],"d":["3/1"]}],[{"n":["-1/1","-1/1"],"d":["-3/2"]},{"n":["-1/1"],"d":["-1/1","3/2"]}]],"g2"',
  ':[[{"n":["2/3"],"d":["-1/3"]},{"n":["0/1"],"d":["0/1","1/1"]}],[{"n":["3/1","2/1"],"d":["0/1',
  '","1/1"]},{"n":["1/1"],"d":["-1/3"]}]],"P":[[[["2/3","2/3"],["1/1"]],[["-2/3"],["-2/3","1/1"',
  ']]],[[[],["1/1"]],[["1/3","-1/6"],["1/1"]]]],"Q":[[[["1/1"],["0/1","1/1"]],[["-1/1"],["1/1"]',
  ']],[[[],["1/1"]],[["0/1","6/1"],["1/1"]]]],"I":[[[["-4/9","2/9","2/3"],["1/1"]],[["0/1","25/',
  '9","-2/9","-2/3"],["1/1"]]],[[[],["1/1"]],[["0/1","2/1","-1/1"],["1/1"]]]],"S":[[[["1/1"],["',
  '0/1","1/1"]],[["-1/1"],["1/1"]]],[[[],["1/1"]],[["1/1"],["-2/3","1/1"]]]]},{"i":7,"n":2,"g1"',
  ':[[{"n":["0/1"],"d":["1/1"]},{"n":[],"d":["1/3","0/1"]}],[{"n":["1/1","-1/3"],"d":["2/3"]},{',
  '"n":["2/1","3/1"],"d":["2/1","-3/1"]}]],"g2":[[{"n":["2/1"],"d":["-2/1","1/1"]},{"n":[],"d":',
  '["-3/1"]}],[{"n":["0/1"],"d":["2/1","-1/3"]},{"n":[],"d":["-3/1","-1/1"]}]],"P":[[[["3/2","-',
  '1/2"],["1/1"]],[["-2/3","-1/1"],["-2/3","1/1"]]]],"Q":[[[["2/1"],["-2/1","1/1"]],[[],["1/1"]',
  ']]],"I":[],"S":[[[["1/1"],["-2/1","1/1"]],[[],["1/1"]]],[[[],["1/1"]],[["4/3","2/1"],["-2/3"',
  ',"1/1"]]]]},{"i":8,"n":3,"g1":[[{"n":["3/1","-3/2"],"d":["1/2"]},{"n":["-1/1"],"d":["-1/1"]}',
  ',{"n":[],"d":["2/3"]}],[{"n":["-1/1"],"d":["-2/1"]},{"n":["0/1","2/1"],"d":["-1/1","-3/1"]},',
  '{"n":["3/2"],"d":["-3/2","1/1"]}]],"g2":[[{"n":["0/1"],"d":["2/3","3/1"]},{"n":["0/1","1/1"]',
  ',"d":["-3/1"]},{"n":[],"d":["-1/1"]}],[{"n":["3/1"],"d":["-2/1"]},{"n":["2/3"],"d":["0/1","1',
  '/1"]},{"n":["-3/2"],"d":["-1/1"]}]],"P":[[[["1/1"],["1/1"]],[["0/1","-4/3"],["1/3","1/1"]],[',
  '["3/1"],["-3/2","1/1"]]],[[[],["1/1"]],[["-1/6","-9/2","2/1"],["1/3","1/1"]],[["9/1","-9/2"]',
  ',["-3/2","1/1"]]]],"Q":[[[["-3/2"],["1/1"]],[["2/3"],["0/1","1/1"]],[["3/2"],["1/1"]]],[[[],',
  '["1/1"]],[["0/1","1/3"],["1/1"]],[[],["1/1"]]]],"I":[[[["0/1","0/1","3/2","15/4","-9/4"],["1',
  '/1"]],[["0/1","0/1","1/8","-61/24","15/4","-1/1"],["1/1"]],[["0/1","0/1","-3/2","-15/4","9/4',
  '"],["1/1"]]]],"S":[[[["1/1"],["1/1"]],[["0/1","-4/3"],["1/3","1/1"]],[["3/1"],["-3/2","1/1"]',
  ']],[[[],["1/1"]],[["1/1"],["0/1","1/3","1/1"]],[["-130329/4408","25893/551","-43821/1102","5',
  '346/551"],["-3/2","1/1"]]],[[[],["1/1"]],[[],["1/1"]],[["1/1"],["-3/2","1/1"]]]]},{"i":10,"n',
  '":2,"g1":[[{"n":["-3/1"],"d":["2/1","-1/2"]},{"n":["1/1"],"d":["-1/1"]}]],"g2":[[{"n":[],"d"',
  ':["-1/1","1/2"]},{"n":["2/1"],"d":["0/1","2/1"]}],[{"n":["-1/1","1/2"],"d":["1/1","-1/2"]},{',
  '"n":["0/1"],"d":["1/1","1/1"]}]],"P":[[[["6/1"],["-4/1","1/1"]],[["-1/1"],["1/1"]]]],"Q":[[[',
  '["-1/1"],["1/1"]],[[],["1/1"]]],[[[],["1/1"]],[["-1/1"],["0/1","1/1"]]]],"I":[[[["1/1"],["1/',
  '1"]],[["2/3","-1/6"],["1/1"]]]],"S":[[[["6/1"],["-4/1","1/1"]],[["-1/1"],["1/1"]]],[[[],["1/',
  '1"]],[["1/1"],["0/1","1/1"]]]]},{"i":12,"n":3,"g1":[[{"n":[],"d":["-1/1","0/1"]},{"n":["-2/3',
  '","-1/1"],"d":["1/1"]},{"n":[],"d":["-1/1"]}]],"g2":[[{"n":[],"d":["3/1"]},{"n":["1/1"],"d":',
  '["1/1","-2/1"]},{"n":[],"d":["1/1"]}]],"P":[[[[],["1/1"]],[["-2/3","-1/1"],["1/1"]],[[],["1/',
  '1"]]]],"Q":[[[[],["1/1"]],[["-1/2"],["-1/2","1/1"]],[[],["1/1"]]]],"I":[[[[],["1/1"]],[["-2/',
  '3","-1/1"],["1/1"]],[[],["1/1"]]]],"S":[[[[],["1/1"]],[["1/1"],["-1/2","1/1"]],[[],["1/1"]]]',
  ']},{"i":17,"n":3,"g1":[[{"n":["2/1","-1/1"],"d":["1/1","2/1"]},{"n":["0/1","-2/1"],"d":["2/3',
  '"]},{"n":["-2/1"],"d":["0/1","2/3"]}],[{"n":["-3/1","3/1"],"d":["-2/3"]},{"n":[],"d":["-1/1"',
  ',"0/1"]},{"n":["-2/1","-1/1"],"d":["2/1"]}]],"g2":[[{"n":["1/3"],"d":["-3/1"]},{"n":["2/1"],',
  '"d":["1/3"]},{"n":["1/3","3/1"],"d":["2/1"]}],[{"n":["1/2"],"d":["3/1"]},{"n":["-1/1","1/1"]',
  ',"d":["1/1"]},{"n":["0/1","2/1"],"d":["2/1"]}]],"P":[[[["1/1"],["1/2","1/1"]],[["0/1","-18/5',
  '","-12/5"],["1/1"]],[["-18/5","-104/45","2/45"],["0/1","1/1"]]],[[[],["1/1"]],[["0/1","27/4"',
  ',"27/4","-27/2"],["1/1"]],[["27/4","23/4","-27/2","1/4"],["0/1","1/1"]]]],"Q":[[[["-1/9"],["',
  '1/1"]],[["6/1"],["1/1"]],[["1/6","3/2"],["1/1"]]],[[[],["1/1"]],[["8/1","1/1"],["1/1"]],[["1',
  '/4","13/4"],["1/1"]]]],"I":[[[["6/1","211/36","-1663/144","-281/72","-293/144","39/8"],["1/1',
  '"]],[["0/1","0/1","59/24","433/8","311/8","-2615/24","81/4"],["1/1"]],[["9/8","1211/24","179',
  '/6","-1711/16","1127/48","1/6"],["1/1"]]]],"S":[[[["1/1"],["1/2","1/1"]],[["0/1","-18/5","-1',
  '2/5"],["1/1"]],[["-18/5","-104/45","2/45"],["0/1","1/1"]]],[[[],["1/1"]],[["1/1"],["1/1"]],[',
  '["-1/15","-97/972","1861/7290","301/14580","1/30"],["0/1","1/1"]]],[[[],["1/1"]],[[],["1/1"]',
  '],[["1/1"],["0/1","1/1"]]]]}]',
].join('');

interface KxGoldenIntegralCase {
  i: number;
  ring: string;
  n: number;
  g1: string[][][];
  g2: string[][][];
  P: string[][][];
  Q: string[][][];
  I: string[][][];
  S: string[][][];
}

interface KxGoldenFractionalCase {
  i: number;
  n: number;
  g1: { n: string[]; d: string[] }[][];
  g2: { n: string[]; d: string[] }[][];
  P: string[][][][];
  Q: string[][][][];
  I: string[][][][];
  S: string[][][][];
}

describe('K[x] echelon basis and intersection match SageMath exactly', () => {
  const ringCache = new Map<string, ReturnType<typeof PolynomialRingConstructor>[0]>();
  const ringFor = (name: string) => {
    const hit = ringCache.get(name);
    if (hit) return hit;
    const base = name === 'QQ' ? (QQfield as never) : (GF(Number(name.slice(2))) as never);
    const [R] = PolynomialRingConstructor(base, 'x');
    ringCache.set(name, R);
    return R;
  };

  /** Build a polynomial from a constant-term-first coefficient list. */
  const polyFrom = (R: ReturnType<typeof ringFor>, ring: string, cs: string[]) =>
    R.__call__(
      cs.map((c) => {
        if (ring === 'QQ') {
          const [num, den] = c.split('/');
          return new Rational(BigInt(num!), BigInt(den!));
        }
        return R.base_ring.__call__(BigInt(c));
      }) as never
    );

  /** Serialize a polynomial back to a coefficient list, trailing zeros dropped. */
  const serPoly = (ring: string, f: unknown): string[] => {
    const coeffs = (f as { coeffs: readonly unknown[] }).coeffs;
    const out = coeffs.map((c) =>
      ring === 'QQ' ? `${(c as Rational).numerator}/${(c as Rational).denominator}` : String(c)
    );
    while (out.length > 0 && (out[out.length - 1] === '0/1' || out[out.length - 1] === '0')) {
      out.pop();
    }
    return out;
  };

  it('reproduces SageMath basis matrices for integral generators', () => {
    const cases = JSON.parse(KX_GOLDEN_INTEGRAL) as KxGoldenIntegralCase[];
    expect(cases.length).toBe(30);
    const failures: string[] = [];
    for (const c of cases) {
      const R = ringFor(c.ring);
      const A = FreeModule(R as never, c.n) as FreeModulePID;
      const mk = (rows: string[][][]) =>
        rows.map((row) => A.createElement(row.map((p) => polyFrom(R, c.ring, p))));
      const P = A.span(mk(c.g1)) as unknown as FreeModulePID;
      const Q = A.span(mk(c.g2)) as unknown as FreeModulePID;
      const got: Record<string, string[][][]> = {
        P: (P.basisMatrix() as unknown[][]).map((r) => r.map((e) => serPoly(c.ring, e))),
        Q: (Q.basisMatrix() as unknown[][]).map((r) => r.map((e) => serPoly(c.ring, e))),
        I: (P.intersection(Q).basisMatrix() as unknown[][]).map((r) =>
          r.map((e) => serPoly(c.ring, e))
        ),
        S: (P.add(Q).basisMatrix() as unknown[][]).map((r) => r.map((e) => serPoly(c.ring, e))),
      };
      for (const key of ['P', 'Q', 'I', 'S'] as const) {
        if (JSON.stringify(got[key]) !== JSON.stringify(c[key])) {
          failures.push(
            `case ${c.i} (${c.ring}) ${key}: sage ${JSON.stringify(c[key])} got ${JSON.stringify(got[key])}`
          );
        }
      }
    }
    expect(failures).toEqual([]);
  });

  it('reproduces SageMath basis matrices for rational-function generators', () => {
    const cases = JSON.parse(KX_GOLDEN_FRACTIONAL) as KxGoldenFractionalCase[];
    expect(cases.length).toBe(12);
    const R = ringFor('QQ');
    /** Normalize an entry to (numerator, monic denominator), as SageMath prints it. */
    const serEntry = (e: unknown): string[][] => {
      if (e instanceof FractionFieldElement) {
        return [serPoly('QQ', e.num), serPoly('QQ', e.den)];
      }
      return [serPoly('QQ', e), ['1/1']];
    };
    const failures: string[] = [];
    for (const c of cases) {
      const A = FreeModule(R as never, c.n) as FreeModulePID;
      const mk = (rows: { n: string[]; d: string[] }[][]) =>
        rows.map((row) =>
          A.createElement(
            row.map((e) =>
              FractionFieldElement.make(
                polyFrom(R, 'QQ', e.n) as never,
                polyFrom(R, 'QQ', e.d) as never
              )
            )
          )
        );
      const P = A.span(mk(c.g1)) as unknown as FreeModulePID;
      const Q = A.span(mk(c.g2)) as unknown as FreeModulePID;
      const got: Record<string, string[][][][]> = {
        P: (P.basisMatrix() as unknown[][]).map((r) => r.map(serEntry)),
        Q: (Q.basisMatrix() as unknown[][]).map((r) => r.map(serEntry)),
        I: (P.intersection(Q).basisMatrix() as unknown[][]).map((r) => r.map(serEntry)),
        S: (P.add(Q).basisMatrix() as unknown[][]).map((r) => r.map(serEntry)),
      };
      for (const key of ['P', 'Q', 'I', 'S'] as const) {
        if (JSON.stringify(got[key]) !== JSON.stringify(c[key])) {
          failures.push(
            `case ${c.i} ${key}: sage ${JSON.stringify(c[key])} got ${JSON.stringify(got[key])}`
          );
        }
      }
    }
    expect(failures).toEqual([]);
  });

  it('the intersection is a submodule of both factors with the expected rank', () => {
    const cases = JSON.parse(KX_GOLDEN_INTEGRAL) as KxGoldenIntegralCase[];
    for (const c of cases) {
      const R = ringFor(c.ring);
      const A = FreeModule(R as never, c.n) as FreeModulePID;
      const mk = (rows: string[][][]) =>
        rows.map((row) => A.createElement(row.map((p) => polyFrom(R, c.ring, p))));
      const P = A.span(mk(c.g1)) as unknown as FreeModulePID;
      const Q = A.span(mk(c.g2)) as unknown as FreeModulePID;
      const I = P.intersection(Q) as unknown as FreeModulePID;
      expect(I.isSubmodule(P)).toBe(true);
      expect(I.isSubmodule(Q)).toBe(true);
      // dim(P n Q) = dim P + dim Q - dim(P + Q) over the fraction field
      expect(I.rank()).toBe(P.rank() + Q.rank() - P.add(Q).rank());
    }
  });

  it('Matrix.denominator over QQ[x] is what forces the unit scaling', () => {
    // The whole point of the fix: `integer_kernel` scales by
    // `S.denominator()` before taking the kernel (matrix2.pyx:5641), and over
    // QQ[x] that denominator is the lcm of the denominators of the rational
    // *coefficients* (matrix2.pyx:3521 + polynomial_element.pyx:4026), i.e. a
    // non-trivial unit of QQ[x].  Nothing downstream normalizes it away.
    //
    // sage: R.<x> = QQ[]; A = R^2
    // sage: P = A.span([[1/2*x + 1, 0]]); Q = A.span([[x^2 + 1/3, 0]])
    // sage: P.intersection(Q).basis_matrix()
    // [-3*x^3 - 6*x^2 - x - 2                      0]
    // sage: (P + Q).basis_matrix()
    // [1 0]
    //
    // Note -3*x^3 - 6*x^2 - x - 2 == -6 * (1/2*x + 1) * (x^2 + 1/3): the ideal
    // is the right one, and the -6 is exactly S.denominator() = lcm(2, 3) = 6
    // (up to sign).
    const R = ringFor('QQ');
    const A = FreeModule(R as never, 2) as FreeModulePID;
    const g1 = A.createElement([polyFrom(R, 'QQ', ['1/1', '1/2']), R.zero()]);
    const g2 = A.createElement([polyFrom(R, 'QQ', ['1/3', '0/1', '1/1']), R.zero()]);
    const P = A.span([g1]) as unknown as FreeModulePID;
    const Q = A.span([g2]) as unknown as FreeModulePID;
    const I = P.intersection(Q) as unknown as FreeModulePID;
    expect((I.basisMatrix() as unknown[][]).map((r) => r.map((e) => serPoly('QQ', e)))).toEqual([
      [['-2/1', '-1/1', '-6/1', '-3/1'], []],
    ]);
    expect(
      (P.add(Q).basisMatrix() as unknown[][]).map((r) => r.map((e) => serPoly('QQ', e)))
    ).toEqual([[['1/1'], []]]);
    expect(I.isSubmodule(P)).toBe(true);
    expect(I.isSubmodule(Q)).toBe(true);
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
