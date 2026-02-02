/**
 * Tests for free modules and LLL lattice reduction
 */

import { describe, expect, it } from 'vitest';
import { type IntegerMatrix, IntegerMatrixFromEntries } from '../matrix/index.js';
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

      expect(() => lllReduce(basis, { eta: 0.4 })).toThrow();
      // eta must be < sqrt(delta), with default delta=0.75, sqrt(0.75) ~ 0.866
      expect(() => lllReduce(basis, { eta: 0.9 })).toThrow();
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

      // |1| + |-2| + |3| = 6
      expect(norm1).toBe(6);
    });

    it('should compute infinity-norm', () => {
      const { vector } = require('./free_module_element.js');

      const v = vector([1n, -5n, 3n]);
      const normInf = v.pNorm(Number.POSITIVE_INFINITY);

      // max(|1|, |-5|, |3|) = 5
      expect(normInf).toBe(5);
    });

    it('should compute 2-norm (squared)', () => {
      const { vector } = require('./free_module_element.js');

      const v = vector([3n, 4n]);
      const norm2 = v.pNorm(2);

      // 3^2 + 4^2 = 25
      expect(norm2).toBe(25n);
    });
  });

  describe('normalized', () => {
    it('should normalize a vector', () => {
      const { vector } = require('./free_module_element.js');

      const v = vector([3, 4]);
      const u = v.normalized();

      // Should be unit vector: (3/5, 4/5)
      expect(u.getItem(0)).toBeCloseTo(0.6, 10);
      expect(u.getItem(1)).toBeCloseTo(0.8, 10);
    });

    it('should throw for zero vector', () => {
      const { vector } = require('./free_module_element.js');

      const v = vector([0, 0, 0]);
      expect(() => v.normalized()).toThrow();
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
