/**
 * Tests for MatrixSpace.
 *
 * Port of the doctests in sage/matrix/matrix_space.py and sage/matrix/args.pyx.
 */

import { describe, expect, it } from 'vitest';
import { GF } from '../rings/finite_rings/finite_field_constructor.js';
import { QQ } from '../rings/rational_field.js';
import { MatrixFromEntries, MatrixSpace, matrix } from './matrix_space.js';

describe('MatrixSpace construction', () => {
  const F7 = GF(7n);

  it('should default to a square space', () => {
    const MS = MatrixSpace(F7, 3);
    expect(MS.dims()).toEqual([3, 3]);
    expect(MS.is_square()).toBe(true);
    expect(MS.dimension()).toBe(9);
  });

  it('should reject negative dimensions', () => {
    expect(() => MatrixSpace(F7, -1)).toThrow(/non-negative/);
    expect(() => MatrixSpace(F7, 2, -3)).toThrow(/non-negative/);
  });

  it('should print like Sage', () => {
    expect(MatrixSpace(F7, 2, 3).toString()).toBe(
      `Full MatrixSpace of 2 by 3 dense matrices over ${F7}`
    );
  });
});

describe('MatrixSpace.__call__', () => {
  const F7 = GF(7n);

  it('should return the zero matrix for no argument', () => {
    expect(MatrixSpace(F7, 2, 3).__call__().is_zero()).toBe(true);
  });

  it('should build a scalar matrix from a single nonzero value', () => {
    // matrix/args.pyx: MA_ENTRIES_SCALAR puts the value on the diagonal only.
    const S = MatrixSpace(F7, 3, 3).__call__(5);
    for (let i = 0; i < 3; i++) {
      for (let j = 0; j < 3; j++) {
        expect(S.get(i, j).value).toBe(i === j ? 5n : 0n);
      }
    }
    expect(S.eq(MatrixSpace(F7, 3, 3).identity().scalar_mul(F7.__call__(5n)))).toBe(true);
  });

  it('should reject a nonzero scalar for a non-square space', () => {
    // sage: MatrixSpace(ZZ, 2, 3)(1)
    // TypeError: nonzero scalar matrix must be square
    expect(() => MatrixSpace(F7, 2, 3).__call__(5)).toThrow('nonzero scalar matrix must be square');
  });

  it('should accept the scalar 0 for a non-square space', () => {
    const Z = MatrixSpace(F7, 2, 3).__call__(0);
    expect(Z.nrows).toBe(2);
    expect(Z.ncols).toBe(3);
    expect(Z.is_zero()).toBe(true);
  });

  it('should build from a 2D array', () => {
    const A = MatrixSpace(F7, 2, 3).__call__([
      [1, 2, 3],
      [4, 5, 6],
    ]);
    expect(A.list().map((x) => x.value)).toEqual([1n, 2n, 3n, 4n, 5n, 6n]);
  });

  it('should build from a flat array', () => {
    const A = MatrixSpace(F7, 2, 3).__call__([1, 2, 3, 4, 5, 6]);
    expect(A.get(1, 0).value).toBe(4n);
  });

  it('should reject the wrong number of entries', () => {
    expect(() => MatrixSpace(F7, 2, 3).__call__([1, 2, 3])).toThrow(/expected 6 entries/);
    expect(() =>
      MatrixSpace(F7, 2, 3).__call__([
        [1, 2, 3],
        [4, 5],
      ])
    ).toThrow(/expected 3 columns/);
    expect(() => MatrixSpace(F7, 2, 3).__call__([[1, 2, 3]])).toThrow(/expected 2 rows/);
  });

  it('should coerce entries into the base ring', () => {
    const A = MatrixSpace(F7, 1, 2).__call__([[8, -1]]);
    expect(A.get(0, 0).value).toBe(1n);
    expect(A.get(0, 1).value).toBe(6n);
  });
});

describe('MatrixSpace helpers', () => {
  const F7 = GF(7n);

  it('should give zero and identity matrices', () => {
    expect(MatrixSpace(F7, 2, 2).zero().is_zero()).toBe(true);
    expect(
      MatrixSpace(F7, 2, 2)
        .one()
        .eq(MatrixSpace(F7, 2, 2).identity())
    ).toBe(true);
    expect(() => MatrixSpace(F7, 2, 3).identity()).toThrow(/square/);
  });

  it('should produce random matrices of the right shape', () => {
    const A = MatrixSpace(F7, 3, 4).random_element();
    expect(A.nrows).toBe(3);
    expect(A.ncols).toBe(4);
  });
});

describe('MatrixFromEntries', () => {
  it('should infer the dimensions', () => {
    const A = MatrixFromEntries(QQ, [
      [1, 2, 3],
      [4, 5, 6],
    ]);
    expect(A.nrows).toBe(2);
    expect(A.ncols).toBe(3);
    expect(A.get(1, 2).eq(6)).toBe(true);
  });

  it('should reject ragged input', () => {
    expect(() => MatrixFromEntries(QQ, [[1, 2], [3]])).toThrow(/inconsistent row lengths/);
  });

  it('should handle the empty matrix', () => {
    const A = MatrixFromEntries(QQ, []);
    expect(A.nrows).toBe(0);
    expect(A.ncols).toBe(0);
  });

  it('should be aliased as matrix', () => {
    expect(matrix).toBe(MatrixFromEntries);
  });
});
