/**
 * Tests for the generic dense matrix class.
 *
 * Port of the doctests in sage/matrix/matrix0.pyx and
 * sage/matrix/matrix_generic_dense.py.
 */

import { describe, expect, it } from 'vitest';
import { GF } from '../rings/finite_rings/finite_field_constructor.js';
import { QQ } from '../rings/rational_field.js';
import {
  Matrix,
  diagonal_matrix,
  identity_matrix,
  scalar_matrix,
  zero_matrix,
} from './matrix_generic.js';
import { MatrixSpace } from './matrix_space.js';

describe('Matrix.pow', () => {
  const F7 = GF(7n);
  const MS2 = MatrixSpace(F7, 2, 2);

  it('should return the identity for exponent 0', () => {
    // Sage follows Python's convention 0^0 = 1:
    //   sage: Matrix([[1,0],[0,0]])^0 == identity_matrix(2)
    const A = MS2.__call__([
      [1, 0],
      [0, 0],
    ]);
    expect(A.pow(0).eq(MS2.identity())).toBe(true);
  });

  it('should compute positive powers by binary exponentiation', () => {
    const A = MS2.__call__([
      [1, 2],
      [3, 4],
    ]);
    expect(A.pow(3).eq(A.mul(A).mul(A))).toBe(true);
    expect(A.pow(5).eq(A.mul(A).mul(A).mul(A).mul(A))).toBe(true);
  });

  it('should compute negative powers as powers of the inverse', () => {
    // matrix0.pyx:__pow__ delegates to generic_power, which inverts first.
    const A = MS2.__call__([
      [1, 2],
      [3, 4],
    ]);
    const Ainv = A.pow(-1);
    expect(A.mul(Ainv).eq(MS2.identity())).toBe(true);
    expect(A.pow(-2).eq(Ainv.mul(Ainv))).toBe(true);
    expect(A.pow(-3).mul(A.pow(3)).eq(MS2.identity())).toBe(true);
  });

  it('should match the Sage doctest A * A^(-1) == 1 over QQ', () => {
    // sage: MS = MatrixSpace(QQ, 3, 3)
    // sage: A = MS([0, 0, 1, 1, 0, '-2/11', 0, 1, '-3/11'])
    // sage: A * A^(-1) == 1
    // True
    const MS3 = MatrixSpace(QQ, 3, 3);
    const A = MS3.__call__([
      [0, 0, 1],
      [1, 0, '-2/11'],
      [0, 1, '-3/11'],
    ]);
    expect(A.mul(A.pow(-1)).eq(MS3.identity())).toBe(true);
  });

  it('should raise when a negative power is asked of a singular matrix', () => {
    const S = MS2.__call__([
      [1, 2],
      [2, 4],
    ]);
    expect(() => S.pow(-1)).toThrow(/singular/);
  });

  it('should raise for non-square matrices', () => {
    expect(() => MatrixSpace(F7, 2, 3).zero().pow(2)).toThrow('self must be a square matrix');
  });
});

describe('Matrix basics', () => {
  const F7 = GF(7n);

  it('should build from a 2D array, a flat array and a function', () => {
    const fromRows = new Matrix(F7, 2, 2, [
      [F7.__call__(1n), F7.__call__(2n)],
      [F7.__call__(3n), F7.__call__(4n)],
    ]);
    const fromFlat = new Matrix(F7, 2, 2, [
      F7.__call__(1n),
      F7.__call__(2n),
      F7.__call__(3n),
      F7.__call__(4n),
    ]);
    const fromFn = new Matrix(F7, 2, 2, (i, j) => F7.__call__(BigInt(2 * i + j + 1)));
    expect(fromRows.eq(fromFlat)).toBe(true);
    expect(fromRows.eq(fromFn)).toBe(true);
  });

  it('should reject inconsistent shapes', () => {
    expect(() => new Matrix(F7, 2, 2, [F7.__call__(1n)])).toThrow(/expected 4 entries/);
    expect(() => new Matrix(F7, -1, 2)).toThrow(/non-negative/);
  });

  it('should add, subtract, negate, multiply and transpose', () => {
    const MS = MatrixSpace(F7, 2, 2);
    const A = MS.__call__([
      [1, 2],
      [3, 4],
    ]);
    const B = MS.__call__([
      [5, 6],
      [0, 1],
    ]);
    expect(A.add(B).sub(B).eq(A)).toBe(true);
    expect(A.add(A.neg()).is_zero()).toBe(true);
    expect(A.mul(MS.identity()).eq(A)).toBe(true);
    expect(A.transpose().transpose().eq(A)).toBe(true);
    expect(A.mul(B).transpose().eq(B.transpose().mul(A.transpose()))).toBe(true);
  });

  it('should reject mismatched shapes in arithmetic', () => {
    const A = MatrixSpace(F7, 2, 2).zero();
    const C = MatrixSpace(F7, 3, 3).zero();
    expect(() => A.add(C)).toThrow(/cannot add/);
    expect(() => A.sub(C)).toThrow(/cannot subtract/);
    expect(() => A.mul(C)).toThrow(/cannot multiply/);
  });

  it('should compute the trace and the diagonal', () => {
    const A = MatrixSpace(F7, 3, 3).__call__([
      [1, 2, 3],
      [4, 5, 6],
      [0, 1, 2],
    ]);
    expect(A.trace().eq(1)).toBe(true); // 1 + 5 + 2 = 8 = 1 mod 7
    expect(A.diagonal().map((d) => d.value)).toEqual([1n, 5n, 2n]);
    expect(() => MatrixSpace(F7, 2, 3).zero().trace()).toThrow(/square/);
  });

  it('should expose rows, columns and a flat list', () => {
    const A = MatrixSpace(F7, 2, 3).__call__([
      [1, 2, 3],
      [4, 5, 6],
    ]);
    expect(A.row(1).map((x) => x.value)).toEqual([4n, 5n, 6n]);
    expect(A.column(2).map((x) => x.value)).toEqual([3n, 6n]);
    expect(A.list().map((x) => x.value)).toEqual([1n, 2n, 3n, 4n, 5n, 6n]);
    expect(A.rows().length).toBe(2);
    expect(A.columns().length).toBe(3);
    expect(() => A.row(2)).toThrow(/out of bounds/);
    expect(() => A.column(3)).toThrow(/out of bounds/);
  });

  it('should copy without aliasing', () => {
    const A = MatrixSpace(F7, 2, 2).__call__([
      [1, 2],
      [3, 4],
    ]);
    const B = A.copy();
    B.set(0, 0, F7.__call__(6n));
    expect(A.get(0, 0).eq(1)).toBe(true);
    expect(B.get(0, 0).eq(6)).toBe(true);
  });
});

describe('matrix constructors', () => {
  const F7 = GF(7n);

  it('should build zero, identity, diagonal and scalar matrices', () => {
    expect(zero_matrix(F7, 2, 3).is_zero()).toBe(true);
    expect(zero_matrix(F7, 3).nrows).toBe(3);
    expect(zero_matrix(F7, 3).ncols).toBe(3);

    const I = identity_matrix(F7, 3);
    expect(I.diagonal().every((d) => d.eq(1))).toBe(true);

    const D = diagonal_matrix(F7, [F7.__call__(1n), F7.__call__(2n), F7.__call__(3n)]);
    expect(D.diagonal().map((d) => d.value)).toEqual([1n, 2n, 3n]);
    expect(D.get(0, 1).isZero()).toBe(true);

    const S = scalar_matrix(F7, 3, F7.__call__(5n));
    expect(S.eq(I.scalar_mul(F7.__call__(5n)))).toBe(true);
  });
});
