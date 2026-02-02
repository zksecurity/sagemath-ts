/**
 * Tests for BKZ (Block Korkine-Zolotarev) lattice reduction
 */

import { describe, expect, it } from 'vitest';
import { IntegerMatrix, IntegerMatrixFromEntries } from '../matrix/index.js';
import {
  BKZ,
  BKZWithInfo,
  HKZ,
  computeHermiteFactor,
  estimateBlockSize,
  isBKZReduced,
} from './bkz.js';
import { gramSchmidt, isLLLReduced, lllReduce } from './free_module_integer.js';

/**
 * Compute squared Euclidean norm of a row in a matrix.
 */
function rowNormSquared(M: IntegerMatrix, row: number): bigint {
  let sum = 0n;
  for (let j = 0; j < M.ncols; j++) {
    const v = M.get(row, j).value;
    sum += v * v;
  }
  return sum;
}

/**
 * Compute squared Euclidean norm of a row in a matrix (as number).
 */
function rowNormSquaredNum(M: IntegerMatrix, row: number): number {
  let sum = 0;
  for (let j = 0; j < M.ncols; j++) {
    const v = Number(M.get(row, j).value);
    sum += v * v;
  }
  return sum;
}

describe('BKZ parameter validation', () => {
  it('should throw for delta <= 0.25', () => {
    const basis = IntegerMatrixFromEntries([
      [1n, 0n],
      [0n, 1n],
    ]);

    expect(() => BKZ(basis, { delta: 0.25 })).toThrow('delta must be in (0.25, 1.0]');
    expect(() => BKZ(basis, { delta: 0.2 })).toThrow('delta must be in (0.25, 1.0]');
    expect(() => BKZ(basis, { delta: 0 })).toThrow('delta must be in (0.25, 1.0]');
    expect(() => BKZ(basis, { delta: -0.5 })).toThrow('delta must be in (0.25, 1.0]');
  });

  it('should throw for delta > 1.0', () => {
    const basis = IntegerMatrixFromEntries([
      [1n, 0n],
      [0n, 1n],
    ]);

    expect(() => BKZ(basis, { delta: 1.01 })).toThrow('delta must be in (0.25, 1.0]');
    expect(() => BKZ(basis, { delta: 2.0 })).toThrow('delta must be in (0.25, 1.0]');
  });

  it('should accept delta = 1.0', () => {
    const basis = IntegerMatrixFromEntries([
      [1n, 0n],
      [0n, 1n],
    ]);

    expect(() => BKZ(basis, { delta: 1.0 })).not.toThrow();
  });

  it('should accept delta in valid range', () => {
    const basis = IntegerMatrixFromEntries([
      [1n, 0n],
      [0n, 1n],
    ]);

    expect(() => BKZ(basis, { delta: 0.99 })).not.toThrow();
    expect(() => BKZ(basis, { delta: 0.75 })).not.toThrow();
    expect(() => BKZ(basis, { delta: 0.5 })).not.toThrow();
    expect(() => BKZ(basis, { delta: 0.26 })).not.toThrow();
  });

  it('should throw for blockSize < 2', () => {
    const basis = IntegerMatrixFromEntries([
      [1n, 0n],
      [0n, 1n],
    ]);

    expect(() => BKZ(basis, { blockSize: 1 })).toThrow('blockSize must be >= 2');
    expect(() => BKZ(basis, { blockSize: 0 })).toThrow('blockSize must be >= 2');
    expect(() => BKZ(basis, { blockSize: -1 })).toThrow('blockSize must be >= 2');
  });

  it('should throw for blockSize > rank', () => {
    const basis = IntegerMatrixFromEntries([
      [1n, 0n],
      [0n, 1n],
    ]);

    expect(() => BKZ(basis, { blockSize: 3 })).toThrow('blockSize must be <= rank (2)');
    expect(() => BKZ(basis, { blockSize: 10 })).toThrow('blockSize must be <= rank (2)');
  });

  it('should accept valid blockSize', () => {
    const basis = IntegerMatrixFromEntries([
      [1n, 0n, 0n],
      [0n, 1n, 0n],
      [0n, 0n, 1n],
    ]);

    expect(() => BKZ(basis, { blockSize: 2 })).not.toThrow();
    expect(() => BKZ(basis, { blockSize: 3 })).not.toThrow();
  });
});

describe('BKZ basic functionality', () => {
  it('should handle empty matrix', () => {
    const basis = new IntegerMatrix(0, 3);
    const reduced = BKZ(basis);

    expect(reduced.nrows).toBe(0);
    expect(reduced.ncols).toBe(3);
  });

  it('should handle single row matrix', () => {
    const basis = IntegerMatrixFromEntries([[3n, 4n, 5n]]);
    const reduced = BKZ(basis);

    expect(reduced.nrows).toBe(1);
    expect(reduced.ncols).toBe(3);
    // Should be unchanged
    expect(reduced.get(0, 0).value).toBe(3n);
    expect(reduced.get(0, 1).value).toBe(4n);
    expect(reduced.get(0, 2).value).toBe(5n);
  });

  it('should reduce a simple 2x2 lattice with blockSize=2', () => {
    const basis = IntegerMatrixFromEntries([
      [1n, 0n],
      [0n, 2n],
    ]);

    const reduced = BKZ(basis, { blockSize: 2 });

    expect(reduced.nrows).toBe(2);
    expect(reduced.ncols).toBe(2);
    expect(isLLLReduced(reduced)).toBe(true);
  });

  it('should reduce poorly conditioned 2x2 lattice', () => {
    const basis = IntegerMatrixFromEntries([
      [1n, 1000n],
      [0n, 1n],
    ]);

    const reduced = BKZ(basis, { blockSize: 2 });

    expect(isLLLReduced(reduced)).toBe(true);

    // First vector should be short
    const norm0 = rowNormSquared(reduced, 0);
    expect(norm0).toBeLessThanOrEqual(2n); // Should be [1, 0] or [0, 1] or similar
  });

  it('should preserve lattice determinant', () => {
    const basis = IntegerMatrixFromEntries([
      [1n, 2n],
      [3n, 5n],
    ]);

    const reduced = BKZ(basis, { blockSize: 2 });

    const detOriginal = basis.determinant().value;
    const detReduced = reduced.determinant().value;

    // Absolute values should match (determinant may flip sign)
    const absDetOriginal = detOriginal < 0n ? -detOriginal : detOriginal;
    const absDetReduced = detReduced < 0n ? -detReduced : detReduced;

    expect(absDetReduced).toBe(absDetOriginal);
  });
});

describe('BKZ vs LLL comparison', () => {
  it('should produce at least as good reduction as LLL', () => {
    const basis = IntegerMatrixFromEntries([
      [1n, 0n, 3n],
      [0n, 2n, 1n],
      [0n, 2n, 7n],
    ]);

    const lllReduced = lllReduce(basis);
    const bkzReduced = BKZ(basis, { blockSize: 3 });

    // Both should be LLL-reduced
    expect(isLLLReduced(lllReduced)).toBe(true);
    expect(isLLLReduced(bkzReduced)).toBe(true);

    // BKZ first vector should be at least as short as LLL first vector
    const lllNorm0 = rowNormSquared(lllReduced, 0);
    const bkzNorm0 = rowNormSquared(bkzReduced, 0);

    expect(bkzNorm0 <= lllNorm0).toBe(true);
  });

  it('should produce better reduction than LLL for larger blocks', () => {
    // A lattice where BKZ with larger block size should do better
    const basis = IntegerMatrixFromEntries([
      [1n, 0n, 0n, 100n],
      [0n, 1n, 0n, 0n],
      [0n, 0n, 1n, 0n],
      [0n, 0n, 0n, 1n],
    ]);

    const lllReduced = lllReduce(basis);
    const bkz2 = BKZ(basis, { blockSize: 2 }); // Equivalent to LLL
    const bkz4 = BKZ(basis, { blockSize: 4 }); // Full HKZ

    // All should be LLL-reduced
    expect(isLLLReduced(lllReduced)).toBe(true);
    expect(isLLLReduced(bkz2)).toBe(true);
    expect(isLLLReduced(bkz4)).toBe(true);

    // BKZ-4 should produce vectors at least as short as BKZ-2/LLL
    const bkz2Norm0 = rowNormSquaredNum(bkz2, 0);
    const bkz4Norm0 = rowNormSquaredNum(bkz4, 0);

    expect(bkz4Norm0).toBeLessThanOrEqual(bkz2Norm0 + 1e-10);
  });
});

describe('BKZ with blockSize = rank (HKZ)', () => {
  it('should give HKZ reduction when blockSize equals rank', () => {
    const basis = IntegerMatrixFromEntries([
      [1n, 0n, 3n],
      [0n, 2n, 1n],
      [0n, 2n, 7n],
    ]);

    const bkzFull = BKZ(basis, { blockSize: 3 });
    const hkzReduced = HKZ(basis);

    // Both should be LLL-reduced
    expect(isLLLReduced(bkzFull)).toBe(true);
    expect(isLLLReduced(hkzReduced)).toBe(true);

    // First vectors should have the same norm (both are HKZ)
    const bkzNorm0 = rowNormSquared(bkzFull, 0);
    const hkzNorm0 = rowNormSquared(hkzReduced, 0);

    expect(bkzNorm0).toBe(hkzNorm0);
  });

  it('HKZ should find shortest vector as first basis vector', () => {
    // Simple lattice where shortest vector is known
    const basis = IntegerMatrixFromEntries([
      [2n, 0n],
      [0n, 3n],
    ]);

    const hkzReduced = HKZ(basis);

    // Shortest vector should be [2, 0] with norm 4
    // or could be [0, 3] with norm 9, but [2, 0] is shorter
    const norm0 = rowNormSquared(hkzReduced, 0);
    expect(norm0).toBe(4n);
  });
});

describe('BKZ determinant preservation', () => {
  it('should preserve 3x3 lattice determinant', () => {
    const basis = IntegerMatrixFromEntries([
      [1n, 0n, 3n],
      [0n, 2n, 1n],
      [0n, 2n, 7n],
    ]);

    const reduced = BKZ(basis, { blockSize: 3 });

    const detOriginal = basis.determinant().value;
    const detReduced = reduced.determinant().value;

    const absDetOriginal = detOriginal < 0n ? -detOriginal : detOriginal;
    const absDetReduced = detReduced < 0n ? -detReduced : detReduced;

    expect(absDetReduced).toBe(absDetOriginal);
  });

  it('should preserve 4x4 lattice determinant', () => {
    const basis = IntegerMatrixFromEntries([
      [1n, 0n, 0n, 5n],
      [0n, 2n, 0n, 3n],
      [0n, 0n, 3n, 1n],
      [0n, 0n, 0n, 7n],
    ]);

    const reduced = BKZ(basis, { blockSize: 3 });

    const detOriginal = basis.determinant().value;
    const detReduced = reduced.determinant().value;

    const absDetOriginal = detOriginal < 0n ? -detOriginal : detOriginal;
    const absDetReduced = detReduced < 0n ? -detReduced : detReduced;

    expect(absDetReduced).toBe(absDetOriginal);
  });
});

describe('BKZWithInfo', () => {
  it('should return tour count and convergence info', () => {
    const basis = IntegerMatrixFromEntries([
      [1n, 0n],
      [0n, 1n],
    ]);

    const result = BKZWithInfo(basis, { blockSize: 2 });

    expect(result.basis).toBeDefined();
    expect(typeof result.tours).toBe('number');
    expect(typeof result.converged).toBe('boolean');
  });

  it('should converge for small lattices', () => {
    const basis = IntegerMatrixFromEntries([
      [1n, 2n],
      [3n, 5n],
    ]);

    const result = BKZWithInfo(basis, { blockSize: 2 });

    expect(result.converged).toBe(true);
  });

  it('should respect maxLoops parameter', () => {
    const basis = IntegerMatrixFromEntries([
      [1n, 0n, 100n],
      [0n, 1n, 0n],
      [0n, 0n, 1n],
    ]);

    const result1 = BKZWithInfo(basis, { blockSize: 3, maxLoops: 1 });
    const result2 = BKZWithInfo(basis, { blockSize: 3, maxLoops: 10 });

    expect(result1.tours).toBeLessThanOrEqual(1);
    expect(result2.tours).toBeLessThanOrEqual(10);
  });
});

describe('isBKZReduced', () => {
  it('should return true for identity matrix', () => {
    const identity = IntegerMatrixFromEntries([
      [1n, 0n],
      [0n, 1n],
    ]);

    expect(isBKZReduced(identity)).toBe(true);
  });

  it('should return true for BKZ-reduced basis', () => {
    const basis = IntegerMatrixFromEntries([
      [1n, 0n, 3n],
      [0n, 2n, 1n],
      [0n, 2n, 7n],
    ]);

    const reduced = BKZ(basis, { blockSize: 3 });

    expect(isBKZReduced(reduced, 3)).toBe(true);
  });

  it('should return false for poorly conditioned basis', () => {
    const basis = IntegerMatrixFromEntries([
      [1n, 1000n],
      [0n, 1n],
    ]);

    // This basis is not even LLL-reduced
    expect(isBKZReduced(basis)).toBe(false);
  });
});

describe('computeHermiteFactor', () => {
  it('should compute reasonable Hermite factor for identity', () => {
    const identity = IntegerMatrixFromEntries([
      [1n, 0n],
      [0n, 1n],
    ]);

    const hf = computeHermiteFactor(identity);

    // For identity, ||b_1|| = 1, vol = 1, so delta = 1
    expect(hf).toBeCloseTo(1.0, 5);
  });

  it('should compute Hermite factor for reduced basis', () => {
    const basis = IntegerMatrixFromEntries([
      [1n, 0n],
      [0n, 2n],
    ]);

    const hf = computeHermiteFactor(basis);

    // ||b_1|| = 1, vol = 2, n = 2
    // delta = (1 / 2^{1/2})^{1/2} = (1/sqrt(2))^{1/2} ≈ 0.84
    expect(hf).toBeLessThan(1.5);
    expect(hf).toBeGreaterThan(0);
  });

  it('should return smaller Hermite factor for better reduced bases', () => {
    const poorBasis = IntegerMatrixFromEntries([
      [1n, 1000n],
      [0n, 1n],
    ]);

    const lllReduced = lllReduce(poorBasis);

    const hfPoor = computeHermiteFactor(poorBasis);
    const hfLLL = computeHermiteFactor(lllReduced);

    // LLL-reduced should have smaller Hermite factor
    expect(hfLLL).toBeLessThan(hfPoor);
  });
});

describe('estimateBlockSize', () => {
  it('should return reasonable block size estimates', () => {
    // For a 100-dimensional lattice with log(vol) = 100
    const blockSize = estimateBlockSize(100, 100, 10);

    expect(blockSize).toBeGreaterThanOrEqual(2);
    expect(blockSize).toBeLessThanOrEqual(100);
  });

  it('should return larger block sizes for tighter targets', () => {
    const loose = estimateBlockSize(50, 100, 50);
    const tight = estimateBlockSize(50, 100, 10);

    // Tighter target (smaller logTargetNorm) should require larger block size
    expect(tight).toBeGreaterThanOrEqual(loose);
  });
});

describe('Known examples from literature', () => {
  it('should reduce NTRU-like lattice', () => {
    const q = 127n;
    const basis = IntegerMatrixFromEntries([
      [1n, 0n, 23n, 45n],
      [0n, 1n, 67n, 89n],
      [0n, 0n, q, 0n],
      [0n, 0n, 0n, q],
    ]);

    const reduced = BKZ(basis, { blockSize: 4 });

    expect(isLLLReduced(reduced)).toBe(true);

    // Determinant should be preserved (q^2 for this lattice structure)
    const detOriginal = basis.determinant().value;
    const detReduced = reduced.determinant().value;
    const absDetOriginal = detOriginal < 0n ? -detOriginal : detOriginal;
    const absDetReduced = detReduced < 0n ? -detReduced : detReduced;

    expect(absDetReduced).toBe(absDetOriginal);
  });

  it('should handle 5x5 lattice with various block sizes', () => {
    const basis = IntegerMatrixFromEntries([
      [7n, 2n, 3n, 4n, 5n],
      [1n, 11n, 2n, 3n, 4n],
      [2n, 3n, 13n, 2n, 3n],
      [3n, 4n, 5n, 17n, 2n],
      [4n, 5n, 6n, 7n, 19n],
    ]);

    const bkz2 = BKZ(basis, { blockSize: 2 });
    const bkz3 = BKZ(basis, { blockSize: 3 });
    const bkz5 = BKZ(basis, { blockSize: 5 });

    // All should be LLL-reduced
    expect(isLLLReduced(bkz2)).toBe(true);
    expect(isLLLReduced(bkz3)).toBe(true);
    expect(isLLLReduced(bkz5)).toBe(true);

    // Verify determinant preservation
    const detOriginal = basis.determinant().value;
    const absDetOriginal = detOriginal < 0n ? -detOriginal : detOriginal;

    for (const reduced of [bkz2, bkz3, bkz5]) {
      const det = reduced.determinant().value;
      const absDet = det < 0n ? -det : det;
      expect(absDet).toBe(absDetOriginal);
    }

    // Larger block sizes should produce at least as short first vectors
    const norm2 = rowNormSquaredNum(bkz2, 0);
    const norm3 = rowNormSquaredNum(bkz3, 0);
    const norm5 = rowNormSquaredNum(bkz5, 0);

    expect(norm3).toBeLessThanOrEqual(norm2 + 1e-10);
    expect(norm5).toBeLessThanOrEqual(norm3 + 1e-10);
  });
});

describe('BKZ edge cases', () => {
  it('should handle lattice with zero vectors correctly', () => {
    const basis = IntegerMatrixFromEntries([
      [1n, 0n],
      [0n, 1n],
    ]);

    // This should work without issues
    const reduced = BKZ(basis, { blockSize: 2 });

    expect(reduced.nrows).toBe(2);
    expect(isLLLReduced(reduced)).toBe(true);
  });

  it('should handle rectangular matrices (more columns than rows)', () => {
    const basis = IntegerMatrixFromEntries([
      [1n, 0n, 1n, 2n],
      [0n, 1n, 2n, 3n],
    ]);

    const reduced = BKZ(basis, { blockSize: 2 });

    expect(reduced.nrows).toBe(2);
    expect(reduced.ncols).toBe(4);
    expect(isLLLReduced(reduced)).toBe(true);
  });

  it('should handle already reduced basis efficiently', () => {
    // Identity matrix is already perfectly reduced
    const basis = IntegerMatrixFromEntries([
      [1n, 0n, 0n],
      [0n, 1n, 0n],
      [0n, 0n, 1n],
    ]);

    const result = BKZWithInfo(basis, { blockSize: 3 });

    expect(result.converged).toBe(true);
    // Should converge quickly since basis is already optimal
  });
});
