/**
 * Tests for closest vector problem (CVP) functions
 *
 * Tests closestVector() and approximateClosestVector() methods on IntegerLattice.
 */

import { describe, expect, it } from 'vitest';
import { IntegerLattice, gramSchmidt, lllReduce, isLLLReduced } from './free_module_integer.js';
import { IntegerMatrixFromEntries } from '../matrix/index.js';

describe('closestVector (exact CVP for small lattices)', () => {
  it('should find closest vector in ZZ^2 - SageMath docstring example', () => {
    // From SageMath: L.closest_vector((-6, 5/3)) should return (-6, 2)
    const L = IntegerLattice([[1, 0], [0, 1]]);

    const closest = L.closestVector([-6, 5 / 3]);

    expect(closest[0]).toBe(-6n);
    expect(closest[1]).toBe(2n);
  });

  it('should return exact lattice point when target is a lattice point', () => {
    const L = IntegerLattice([[1, 0], [0, 1]]);

    const closest = L.closestVector([3, 5]);

    expect(closest[0]).toBe(3n);
    expect(closest[1]).toBe(5n);
  });

  it('should handle negative coordinates', () => {
    const L = IntegerLattice([[1, 0], [0, 1]]);

    const closest = L.closestVector([-3.4, -2.7]);

    expect(closest[0]).toBe(-3n);
    expect(closest[1]).toBe(-3n);
  });

  it('should work with scaled lattice', () => {
    const L = IntegerLattice([[2, 0], [0, 3]]);

    // Target (3, 5) should map to closest point in 2Z x 3Z
    const closest = L.closestVector([3, 5]);

    // Check that it's a valid lattice point
    expect(closest[0]! % 2n).toBe(0n);
    expect(closest[1]! % 3n).toBe(0n);

    // And close to target
    const dist = Math.sqrt(
      Math.pow(3 - Number(closest[0]), 2) + Math.pow(5 - Number(closest[1]), 2)
    );
    expect(dist).toBeLessThan(5);
  });

  it('should handle 3D lattice', () => {
    const L = IntegerLattice([
      [1, 0, 0],
      [0, 1, 0],
      [0, 0, 1],
    ]);

    const closest = L.closestVector([1.2, 2.7, 3.1]);

    expect(closest[0]).toBe(1n);
    expect(closest[1]).toBe(3n);
    expect(closest[2]).toBe(3n);
  });

  it('should verify result is close to target', () => {
    const L = IntegerLattice([
      [1, 0, 3],
      [0, 2, 1],
      [0, 2, 7],
    ]);

    const target = [5.5, 3.3, 10.1];
    const closest = L.closestVector(target);

    // The closest vector should be reasonably close to target
    let distSq = 0;
    for (let i = 0; i < 3; i++) {
      distSq += Math.pow(target[i]! - Number(closest[i]!), 2);
    }
    expect(Math.sqrt(distSq)).toBeLessThan(20);
  });

  it('should handle zero target', () => {
    const L = IntegerLattice([[1, 0], [0, 1]]);

    const closest = L.closestVector([0, 0]);

    expect(closest[0]).toBe(0n);
    expect(closest[1]).toBe(0n);
  });

  it('should handle large coordinates', () => {
    const L = IntegerLattice([[1, 0], [0, 1]]);

    const closest = L.closestVector([1000000.3, 2000000.7]);

    expect(closest[0]).toBe(1000000n);
    expect(closest[1]).toBe(2000001n);
  });

  it('should throw for wrong dimension', () => {
    const L = IntegerLattice([[1, 0], [0, 1]]);

    expect(() => L.closestVector([1, 2, 3])).toThrow();
  });
});

describe('approximateClosestVector algorithms', () => {
  describe('nearest_plane algorithm (default)', () => {
    it('should approximate closest vector', () => {
      const L = IntegerLattice([[1, 0], [0, 1]]);

      const approx = L.approximateClosestVector([0.4, 0.6], { algorithm: 'nearest_plane' });

      expect(typeof approx[0]).toBe('bigint');
      expect(typeof approx[1]).toBe('bigint');
    });

    it('should match babai() alias', () => {
      const L = IntegerLattice([[1, 0], [0, 1]]);
      const target = [3.2, 4.7];

      const approx = L.approximateClosestVector(target, { algorithm: 'nearest_plane' });
      const babai = L.babai(target, { algorithm: 'nearest_plane' });

      expect(approx[0]).toBe(babai[0]);
      expect(approx[1]).toBe(babai[1]);
    });

    it('should work with non-trivial lattice', () => {
      const L = IntegerLattice([
        [101, 0, 0, 0],
        [0, 101, 0, 0],
        [0, 0, 101, 0],
        [-28, 39, 45, 1],
      ]);

      const target = [1337, 1337, 1337, 1337];
      const approx = L.approximateClosestVector(target, { algorithm: 'nearest_plane' });

      expect(approx.length).toBe(4);

      // Should be reasonably close to target
      let distSq = 0;
      for (let i = 0; i < 4; i++) {
        distSq += Math.pow(target[i]! - Number(approx[i]!), 2);
      }
      expect(Math.sqrt(distSq)).toBeLessThan(300);
    });
  });

  describe('rounding_off algorithm', () => {
    it('should approximate closest vector', () => {
      const L = IntegerLattice([[1, 0], [0, 1]]);

      const approx = L.approximateClosestVector([0.4, 0.6], { algorithm: 'rounding_off' });

      expect(typeof approx[0]).toBe('bigint');
      expect(typeof approx[1]).toBe('bigint');
    });

    it('should give reasonable approximation for simple lattice', () => {
      const L = IntegerLattice([[1, 0], [0, 1]]);

      const approx = L.approximateClosestVector([3, 5], { algorithm: 'rounding_off' });

      expect(approx[0]).toBe(3n);
      expect(approx[1]).toBe(5n);
    });
  });

  describe('embedding algorithm', () => {
    it('should approximate closest vector', () => {
      const L = IntegerLattice([[1, 0], [0, 1]]);

      const approx = L.approximateClosestVector([0.4, 0.6], { algorithm: 'embedding' });

      expect(typeof approx[0]).toBe('bigint');
      expect(typeof approx[1]).toBe('bigint');
    });

    it('should work for lattice points', () => {
      const L = IntegerLattice([[1, 0], [0, 1]]);

      const approx = L.approximateClosestVector([3, 5], { algorithm: 'embedding' });

      expect(approx[0]).toBe(3n);
      expect(approx[1]).toBe(5n);
    });
  });

  describe('algorithm comparison', () => {
    it('all algorithms should produce valid lattice vectors', () => {
      const L = IntegerLattice([
        [2, 0],
        [0, 3],
      ]);

      const target = [5, 7];

      const nearPlane = L.approximateClosestVector(target, { algorithm: 'nearest_plane' });
      const rounding = L.approximateClosestVector(target, { algorithm: 'rounding_off' });
      const embedding = L.approximateClosestVector(target, { algorithm: 'embedding' });

      // All should be valid lattice vectors (multiples of 2 and 3)
      for (const result of [nearPlane, rounding, embedding]) {
        expect(result[0]! % 2n).toBe(0n);
        expect(result[1]! % 3n).toBe(0n);
      }
    });

    it('results should be reasonably similar', () => {
      const L = IntegerLattice([[1, 0], [0, 1]]);
      const target = [10.5, 20.3];

      const nearPlane = L.approximateClosestVector(target, { algorithm: 'nearest_plane' });
      const rounding = L.approximateClosestVector(target, { algorithm: 'rounding_off' });
      const embedding = L.approximateClosestVector(target, { algorithm: 'embedding' });

      // All should give results within 2 units of each other
      const maxDiff = 2n;
      expect(Math.abs(Number(nearPlane[0]! - rounding[0]!))).toBeLessThanOrEqual(Number(maxDiff));
      expect(Math.abs(Number(nearPlane[1]! - rounding[1]!))).toBeLessThanOrEqual(Number(maxDiff));
      expect(Math.abs(Number(nearPlane[0]! - embedding[0]!))).toBeLessThanOrEqual(Number(maxDiff));
      expect(Math.abs(Number(nearPlane[1]! - embedding[1]!))).toBeLessThanOrEqual(Number(maxDiff));
    });
  });
});

describe('CVP quality metrics', () => {
  it('exact CVP should be optimal for small lattices', () => {
    const L = IntegerLattice([[1, 0], [0, 1]]);

    const target = [0.3, 0.3];
    const closest = L.closestVector(target);

    // The closest point to (0.3, 0.3) in ZZ^2 is (0, 0)
    expect(closest[0]).toBe(0n);
    expect(closest[1]).toBe(0n);
  });

  it('approximate CVP should be within expected bound', () => {
    // For Babai's nearest plane, the approximation factor is at most 2^(n/2)
    const L = IntegerLattice([[1, 0], [0, 1]]);

    const target = [0.6, 0.6];
    const approx = L.approximateClosestVector(target);

    // The optimal closest vector is (1, 1) with distance sqrt(0.16 + 0.16) = 0.566
    const dist = Math.sqrt(
      Math.pow(0.6 - Number(approx[0]), 2) + Math.pow(0.6 - Number(approx[1]), 2)
    );

    // Babai should give us a good approximation
    expect(dist).toBeLessThan(2);
  });
});

describe('Babai alias', () => {
  it('babai() should be an alias for approximateClosestVector', () => {
    const L = IntegerLattice([[1, 0], [0, 1]]);
    const target = [3.2, 4.7];

    const approx = L.approximateClosestVector(target);
    const babai = L.babai(target);

    expect(approx[0]).toBe(babai[0]);
    expect(approx[1]).toBe(babai[1]);
  });
});
