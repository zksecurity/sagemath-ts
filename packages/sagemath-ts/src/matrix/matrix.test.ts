/**
 * Tests for matrix operations
 */

import { describe, expect, it } from 'vitest';
import { GF } from '../rings/finite_rings/finite_field_constructor.js';
import { Integer } from '../rings/integer_ring.js';
import {
  IntegerMatrix,
  IntegerMatrixFromEntries,
  Matrix,
  MatrixFromEntries,
  MatrixSpace,
  diagonal_matrix,
  identity_integer_matrix,
  identity_matrix,
  scalar_matrix,
  zero_integer_matrix,
  zero_matrix,
} from './index.js';

describe('Matrix over GF(p)', () => {
  const F7 = GF(7n);
  const MS2x2 = MatrixSpace(F7, 2, 2);
  const MS2x3 = MatrixSpace(F7, 2, 3);
  const MS3x2 = MatrixSpace(F7, 3, 2);

  describe('creation and access', () => {
    it('should create a matrix from 2D array', () => {
      const A = MS2x2.__call__([
        [1, 2],
        [3, 4],
      ]);

      expect(A.get(0, 0).value).toBe(1n);
      expect(A.get(0, 1).value).toBe(2n);
      expect(A.get(1, 0).value).toBe(3n);
      expect(A.get(1, 1).value).toBe(4n);
    });

    it('should create a matrix from flat array', () => {
      const A = MS2x2.__call__([1, 2, 3, 4]);

      expect(A.get(0, 0).value).toBe(1n);
      expect(A.get(0, 1).value).toBe(2n);
      expect(A.get(1, 0).value).toBe(3n);
      expect(A.get(1, 1).value).toBe(4n);
    });

    it('should reduce entries modulo p', () => {
      const A = MS2x2.__call__([
        [8, 9],
        [10, 14],
      ]);

      expect(A.get(0, 0).value).toBe(1n); // 8 mod 7 = 1
      expect(A.get(0, 1).value).toBe(2n); // 9 mod 7 = 2
      expect(A.get(1, 0).value).toBe(3n); // 10 mod 7 = 3
      expect(A.get(1, 1).value).toBe(0n); // 14 mod 7 = 0
    });

    it('should get and set entries', () => {
      const A = MS2x2.__call__();
      A.set(0, 1, F7.__call__(5n));
      A.set(1, 0, F7.__call__(3n));

      expect(A.get(0, 0).value).toBe(0n);
      expect(A.get(0, 1).value).toBe(5n);
      expect(A.get(1, 0).value).toBe(3n);
      expect(A.get(1, 1).value).toBe(0n);
    });

    it('should get rows and columns', () => {
      const A = MS2x3.__call__([
        [1, 2, 3],
        [4, 5, 6],
      ]);

      const row0 = A.row(0);
      expect(row0.map((e) => e.value)).toEqual([1n, 2n, 3n]);

      const col1 = A.column(1);
      expect(col1.map((e) => e.value)).toEqual([2n, 5n]);
    });

    it('should throw on out-of-bounds access', () => {
      const A = MS2x2.__call__([
        [1, 2],
        [3, 4],
      ]);
      expect(() => A.get(2, 0)).toThrow();
      expect(() => A.get(0, 2)).toThrow();
      expect(() => A.get(-1, 0)).toThrow();
    });
  });

  describe('zero and identity matrices', () => {
    it('should create zero matrix', () => {
      const Z = MS2x3.zero();

      expect(Z.nrows).toBe(2);
      expect(Z.ncols).toBe(3);
      expect(Z.is_zero()).toBe(true);

      for (let i = 0; i < 2; i++) {
        for (let j = 0; j < 3; j++) {
          expect(Z.get(i, j).value).toBe(0n);
        }
      }
    });

    it('should create identity matrix', () => {
      const I = MS2x2.identity();

      expect(I.get(0, 0).value).toBe(1n);
      expect(I.get(0, 1).value).toBe(0n);
      expect(I.get(1, 0).value).toBe(0n);
      expect(I.get(1, 1).value).toBe(1n);
    });

    it('should throw when creating identity for non-square space', () => {
      expect(() => MS2x3.identity()).toThrow();
    });

    it('should create diagonal matrix', () => {
      const D = diagonal_matrix(F7, [F7.__call__(1n), F7.__call__(2n), F7.__call__(3n)]);

      expect(D.nrows).toBe(3);
      expect(D.ncols).toBe(3);
      expect(D.get(0, 0).value).toBe(1n);
      expect(D.get(1, 1).value).toBe(2n);
      expect(D.get(2, 2).value).toBe(3n);
      expect(D.get(0, 1).value).toBe(0n);
    });

    it('should create scalar matrix', () => {
      const S = scalar_matrix(F7, 3, F7.__call__(5n));

      expect(S.get(0, 0).value).toBe(5n);
      expect(S.get(1, 1).value).toBe(5n);
      expect(S.get(2, 2).value).toBe(5n);
      expect(S.get(0, 1).value).toBe(0n);
    });
  });

  describe('addition and subtraction', () => {
    it('should add two matrices', () => {
      const A = MS2x2.__call__([
        [1, 2],
        [3, 4],
      ]);
      const B = MS2x2.__call__([
        [5, 6],
        [6, 5],
      ]);

      const C = A.add(B);

      expect(C.get(0, 0).value).toBe(6n); // 1+5
      expect(C.get(0, 1).value).toBe(1n); // 2+6 = 8 ≡ 1 mod 7
      expect(C.get(1, 0).value).toBe(2n); // 3+6 = 9 ≡ 2 mod 7
      expect(C.get(1, 1).value).toBe(2n); // 4+5 = 9 ≡ 2 mod 7
    });

    it('should subtract two matrices', () => {
      const A = MS2x2.__call__([
        [5, 6],
        [3, 4],
      ]);
      const B = MS2x2.__call__([
        [1, 2],
        [4, 1],
      ]);

      const C = A.sub(B);

      expect(C.get(0, 0).value).toBe(4n); // 5-1
      expect(C.get(0, 1).value).toBe(4n); // 6-2
      expect(C.get(1, 0).value).toBe(6n); // 3-4 = -1 ≡ 6 mod 7
      expect(C.get(1, 1).value).toBe(3n); // 4-1
    });

    it('should negate a matrix', () => {
      const A = MS2x2.__call__([
        [1, 2],
        [3, 4],
      ]);

      const B = A.neg();

      expect(B.get(0, 0).value).toBe(6n); // -1 ≡ 6 mod 7
      expect(B.get(0, 1).value).toBe(5n); // -2 ≡ 5 mod 7
      expect(B.get(1, 0).value).toBe(4n); // -3 ≡ 4 mod 7
      expect(B.get(1, 1).value).toBe(3n); // -4 ≡ 3 mod 7
    });

    it('should throw when adding matrices of different sizes', () => {
      const A = MS2x2.__call__([
        [1, 2],
        [3, 4],
      ]);
      const B = MS2x3.__call__([
        [1, 2, 3],
        [4, 5, 6],
      ]);
      expect(() => A.add(B as any)).toThrow();
    });
  });

  describe('matrix multiplication', () => {
    it('should multiply two 2x2 matrices', () => {
      const A = MS2x2.__call__([
        [1, 2],
        [3, 4],
      ]);
      const B = MS2x2.__call__([
        [5, 6],
        [0, 1],
      ]);

      const C = A.mul(B);

      // C[0,0] = 1*5 + 2*0 = 5
      // C[0,1] = 1*6 + 2*1 = 8 ≡ 1 mod 7
      // C[1,0] = 3*5 + 4*0 = 15 ≡ 1 mod 7
      // C[1,1] = 3*6 + 4*1 = 22 ≡ 1 mod 7
      expect(C.get(0, 0).value).toBe(5n);
      expect(C.get(0, 1).value).toBe(1n);
      expect(C.get(1, 0).value).toBe(1n);
      expect(C.get(1, 1).value).toBe(1n);
    });

    it('should multiply 2x3 by 3x2 to get 2x2', () => {
      const A = MS2x3.__call__([
        [1, 2, 3],
        [4, 5, 6],
      ]);
      const B = MS3x2.__call__([
        [1, 2],
        [3, 4],
        [5, 6],
      ]);

      const C = A.mul(B);

      expect(C.nrows).toBe(2);
      expect(C.ncols).toBe(2);

      // C[0,0] = 1*1 + 2*3 + 3*5 = 1 + 6 + 15 = 22 ≡ 1 mod 7
      // C[0,1] = 1*2 + 2*4 + 3*6 = 2 + 8 + 18 = 28 ≡ 0 mod 7
      expect(C.get(0, 0).value).toBe(1n);
      expect(C.get(0, 1).value).toBe(0n);
    });

    it('should multiply by identity to get same matrix', () => {
      const A = MS2x2.__call__([
        [1, 2],
        [3, 4],
      ]);
      const I = MS2x2.identity();

      const AI = A.mul(I);
      const IA = I.mul(A);

      expect(AI.eq(A)).toBe(true);
      expect(IA.eq(A)).toBe(true);
    });

    it('should multiply by zero to get zero', () => {
      const A = MS2x2.__call__([
        [1, 2],
        [3, 4],
      ]);
      const Z = MS2x2.zero();

      const AZ = A.mul(Z);
      const ZA = Z.mul(A);

      expect(AZ.is_zero()).toBe(true);
      expect(ZA.is_zero()).toBe(true);
    });

    it('should throw when dimensions are incompatible', () => {
      const A = MS2x2.__call__([
        [1, 2],
        [3, 4],
      ]);
      const B = MS2x2.__call__([
        [1, 2],
        [3, 4],
      ]);
      // 2x2 * 2x2 should work
      expect(() => A.mul(B)).not.toThrow();

      // 2x3 * 2x2 should fail
      const C = MS2x3.__call__([
        [1, 2, 3],
        [4, 5, 6],
      ]);
      expect(() => C.mul(A)).toThrow();
    });

    it('should scalar multiply', () => {
      const A = MS2x2.__call__([
        [1, 2],
        [3, 4],
      ]);

      const B = A.scalar_mul(F7.__call__(3n));

      expect(B.get(0, 0).value).toBe(3n); // 1*3
      expect(B.get(0, 1).value).toBe(6n); // 2*3
      expect(B.get(1, 0).value).toBe(2n); // 3*3 = 9 ≡ 2 mod 7
      expect(B.get(1, 1).value).toBe(5n); // 4*3 = 12 ≡ 5 mod 7
    });
  });

  describe('transpose', () => {
    it('should transpose a square matrix', () => {
      const A = MS2x2.__call__([
        [1, 2],
        [3, 4],
      ]);

      const At = A.transpose();

      expect(At.get(0, 0).value).toBe(1n);
      expect(At.get(0, 1).value).toBe(3n);
      expect(At.get(1, 0).value).toBe(2n);
      expect(At.get(1, 1).value).toBe(4n);
    });

    it('should transpose a non-square matrix', () => {
      const A = MS2x3.__call__([
        [1, 2, 3],
        [4, 5, 6],
      ]);

      const At = A.transpose();

      expect(At.nrows).toBe(3);
      expect(At.ncols).toBe(2);
      expect(At.get(0, 0).value).toBe(1n);
      expect(At.get(0, 1).value).toBe(4n);
      expect(At.get(1, 0).value).toBe(2n);
      expect(At.get(1, 1).value).toBe(5n);
      expect(At.get(2, 0).value).toBe(3n);
      expect(At.get(2, 1).value).toBe(6n);
    });

    it('should satisfy (A^T)^T = A', () => {
      const A = MS2x3.__call__([
        [1, 2, 3],
        [4, 5, 6],
      ]);

      const Att = A.transpose().transpose();
      expect(Att.eq(A)).toBe(true);
    });
  });

  describe('trace', () => {
    it('should compute trace of square matrix', () => {
      const A = MS2x2.__call__([
        [1, 2],
        [3, 4],
      ]);

      const tr = A.trace();
      expect(tr.value).toBe(5n); // 1 + 4 = 5
    });

    it('should compute trace with reduction mod p', () => {
      const A = MS2x2.__call__([
        [5, 2],
        [3, 6],
      ]);

      const tr = A.trace();
      expect(tr.value).toBe(4n); // 5 + 6 = 11 ≡ 4 mod 7
    });

    it('should throw for non-square matrix', () => {
      const A = MS2x3.__call__([
        [1, 2, 3],
        [4, 5, 6],
      ]);

      expect(() => A.trace()).toThrow();
    });
  });

  describe('matrix power', () => {
    it('should compute A^0 = I', () => {
      const A = MS2x2.__call__([
        [1, 2],
        [3, 4],
      ]);

      const A0 = A.pow(0);
      expect(A0.eq(MS2x2.identity())).toBe(true);
    });

    it('should compute A^1 = A', () => {
      const A = MS2x2.__call__([
        [1, 2],
        [3, 4],
      ]);

      const A1 = A.pow(1);
      expect(A1.eq(A)).toBe(true);
    });

    it('should compute A^2 = A*A', () => {
      const A = MS2x2.__call__([
        [1, 2],
        [3, 4],
      ]);

      const A2 = A.pow(2);
      const AA = A.mul(A);
      expect(A2.eq(AA)).toBe(true);
    });

    it('should throw for non-square matrix', () => {
      const A = MS2x3.__call__([
        [1, 2, 3],
        [4, 5, 6],
      ]);

      expect(() => A.pow(2)).toThrow();
    });
  });

  describe('equality', () => {
    it('should detect equal matrices', () => {
      const A = MS2x2.__call__([
        [1, 2],
        [3, 4],
      ]);
      const B = MS2x2.__call__([
        [1, 2],
        [3, 4],
      ]);

      expect(A.eq(B)).toBe(true);
    });

    it('should detect unequal matrices', () => {
      const A = MS2x2.__call__([
        [1, 2],
        [3, 4],
      ]);
      const B = MS2x2.__call__([
        [1, 2],
        [3, 5],
      ]);

      expect(A.eq(B)).toBe(false);
    });

    it('should handle equality with reduction', () => {
      const A = MS2x2.__call__([
        [1, 2],
        [3, 4],
      ]);
      const B = MS2x2.__call__([
        [8, 9],
        [10, 11],
      ]); // 8≡1, 9≡2, 10≡3, 11≡4 mod 7

      expect(A.eq(B)).toBe(true);
    });
  });
});

describe('IntegerMatrix', () => {
  describe('creation and access', () => {
    it('should create an integer matrix', () => {
      const A = new IntegerMatrix(2, 2, [
        [1n, 2n],
        [3n, 4n],
      ]);

      expect(A.get(0, 0).value).toBe(1n);
      expect(A.get(0, 1).value).toBe(2n);
      expect(A.get(1, 0).value).toBe(3n);
      expect(A.get(1, 1).value).toBe(4n);
    });

    it('should create from number array', () => {
      const A = new IntegerMatrix(2, 2, [
        [1, 2],
        [3, 4],
      ]);

      expect(A.get(0, 0).value).toBe(1n);
      expect(A.get(1, 1).value).toBe(4n);
    });

    it('should create using factory function', () => {
      const A = IntegerMatrixFromEntries([
        [1n, 2n],
        [3n, 4n],
      ]);

      expect(A.nrows).toBe(2);
      expect(A.ncols).toBe(2);
      expect(A.get(0, 0).value).toBe(1n);
    });
  });

  describe('arithmetic operations', () => {
    it('should add matrices', () => {
      const A = new IntegerMatrix(2, 2, [
        [1n, 2n],
        [3n, 4n],
      ]);
      const B = new IntegerMatrix(2, 2, [
        [5n, 6n],
        [7n, 8n],
      ]);

      const C = A.add(B);

      expect(C.get(0, 0).value).toBe(6n);
      expect(C.get(0, 1).value).toBe(8n);
      expect(C.get(1, 0).value).toBe(10n);
      expect(C.get(1, 1).value).toBe(12n);
    });

    it('should subtract matrices', () => {
      const A = new IntegerMatrix(2, 2, [
        [5n, 6n],
        [7n, 8n],
      ]);
      const B = new IntegerMatrix(2, 2, [
        [1n, 2n],
        [3n, 4n],
      ]);

      const C = A.sub(B);

      expect(C.get(0, 0).value).toBe(4n);
      expect(C.get(0, 1).value).toBe(4n);
      expect(C.get(1, 0).value).toBe(4n);
      expect(C.get(1, 1).value).toBe(4n);
    });

    it('should multiply matrices', () => {
      const A = new IntegerMatrix(2, 2, [
        [1n, 2n],
        [3n, 4n],
      ]);
      const B = new IntegerMatrix(2, 2, [
        [5n, 6n],
        [7n, 8n],
      ]);

      const C = A.mul(B);

      // C[0,0] = 1*5 + 2*7 = 19
      // C[0,1] = 1*6 + 2*8 = 22
      // C[1,0] = 3*5 + 4*7 = 43
      // C[1,1] = 3*6 + 4*8 = 50
      expect(C.get(0, 0).value).toBe(19n);
      expect(C.get(0, 1).value).toBe(22n);
      expect(C.get(1, 0).value).toBe(43n);
      expect(C.get(1, 1).value).toBe(50n);
    });

    it('should scalar multiply', () => {
      const A = new IntegerMatrix(2, 2, [
        [1n, 2n],
        [3n, 4n],
      ]);
      const B = A.scalar_mul(3n);

      expect(B.get(0, 0).value).toBe(3n);
      expect(B.get(0, 1).value).toBe(6n);
      expect(B.get(1, 0).value).toBe(9n);
      expect(B.get(1, 1).value).toBe(12n);
    });
  });

  describe('determinant', () => {
    it('should compute determinant of 1x1 matrix', () => {
      const A = new IntegerMatrix(1, 1, [[5n]]);
      expect(A.determinant().value).toBe(5n);
    });

    it('should compute determinant of 2x2 matrix', () => {
      const A = new IntegerMatrix(2, 2, [
        [1n, 2n],
        [3n, 4n],
      ]);

      // det = 1*4 - 2*3 = 4 - 6 = -2
      expect(A.determinant().value).toBe(-2n);
    });

    it('should compute determinant of 2x2 matrix (another example)', () => {
      const A = new IntegerMatrix(2, 2, [
        [3n, 8n],
        [4n, 6n],
      ]);

      // det = 3*6 - 8*4 = 18 - 32 = -14
      expect(A.determinant().value).toBe(-14n);
    });

    it('should compute determinant of 3x3 matrix', () => {
      const A = new IntegerMatrix(3, 3, [
        [1n, 2n, 3n],
        [4n, 5n, 6n],
        [7n, 8n, 9n],
      ]);

      // This is a singular matrix (rows are linearly dependent)
      expect(A.determinant().value).toBe(0n);
    });

    it('should compute determinant of 3x3 matrix (non-singular)', () => {
      const A = new IntegerMatrix(3, 3, [
        [1n, 2n, 3n],
        [0n, 1n, 4n],
        [5n, 6n, 0n],
      ]);

      // det = 1*(1*0 - 4*6) - 2*(0*0 - 4*5) + 3*(0*6 - 1*5)
      //     = 1*(-24) - 2*(-20) + 3*(-5)
      //     = -24 + 40 - 15
      //     = 1
      expect(A.determinant().value).toBe(1n);
    });

    it('should compute determinant of identity matrix', () => {
      const I3 = identity_integer_matrix(3);
      expect(I3.determinant().value).toBe(1n);

      const I4 = identity_integer_matrix(4);
      expect(I4.determinant().value).toBe(1n);
    });

    it('should compute determinant of 4x4 matrix', () => {
      const A = new IntegerMatrix(4, 4, [
        [1n, 2n, 3n, 4n],
        [5n, 6n, 7n, 8n],
        [9n, 10n, 11n, 12n],
        [13n, 14n, 15n, 16n],
      ]);

      // This matrix has linearly dependent rows
      expect(A.determinant().value).toBe(0n);
    });

    it('should compute determinant of 4x4 non-singular matrix', () => {
      const A = new IntegerMatrix(4, 4, [
        [1n, 0n, 0n, 0n],
        [0n, 2n, 0n, 0n],
        [0n, 0n, 3n, 0n],
        [0n, 0n, 0n, 4n],
      ]);

      // Diagonal matrix: det = 1*2*3*4 = 24
      expect(A.determinant().value).toBe(24n);
    });

    it('should throw for non-square matrix', () => {
      const A = new IntegerMatrix(2, 3, [
        [1n, 2n, 3n],
        [4n, 5n, 6n],
      ]);
      expect(() => A.determinant()).toThrow();
    });
  });

  describe('trace', () => {
    it('should compute trace', () => {
      const A = new IntegerMatrix(3, 3, [
        [1n, 2n, 3n],
        [4n, 5n, 6n],
        [7n, 8n, 9n],
      ]);

      expect(A.trace().value).toBe(15n); // 1 + 5 + 9
    });
  });

  describe('transpose', () => {
    it('should transpose a matrix', () => {
      const A = new IntegerMatrix(2, 3, [
        [1n, 2n, 3n],
        [4n, 5n, 6n],
      ]);

      const At = A.transpose();

      expect(At.nrows).toBe(3);
      expect(At.ncols).toBe(2);
      expect(At.get(0, 0).value).toBe(1n);
      expect(At.get(0, 1).value).toBe(4n);
      expect(At.get(1, 0).value).toBe(2n);
      expect(At.get(2, 1).value).toBe(6n);
    });
  });

  describe('zero and identity matrices', () => {
    it('should create zero matrix', () => {
      const Z = zero_integer_matrix(3, 2);

      expect(Z.nrows).toBe(3);
      expect(Z.ncols).toBe(2);

      for (let i = 0; i < 3; i++) {
        for (let j = 0; j < 2; j++) {
          expect(Z.get(i, j).value).toBe(0n);
        }
      }
    });

    it('should create identity matrix', () => {
      const I = identity_integer_matrix(3);

      expect(I.nrows).toBe(3);
      expect(I.ncols).toBe(3);
      expect(I.get(0, 0).value).toBe(1n);
      expect(I.get(1, 1).value).toBe(1n);
      expect(I.get(2, 2).value).toBe(1n);
      expect(I.get(0, 1).value).toBe(0n);
    });
  });
});

describe('Matrix over larger finite fields', () => {
  const F101 = GF(101n);
  const MS = MatrixSpace(F101, 3, 3);

  it('should work with GF(101)', () => {
    const A = MS.__call__([
      [50, 60, 70],
      [80, 90, 100],
      [10, 20, 30],
    ]);

    const B = MS.__call__([
      [1, 2, 3],
      [4, 5, 6],
      [7, 8, 9],
    ]);

    const C = A.add(B);

    expect(C.get(0, 0).value).toBe(51n);
    expect(C.get(0, 1).value).toBe(62n);
    expect(C.get(0, 2).value).toBe(73n);
    expect(C.get(1, 0).value).toBe(84n);
    expect(C.get(1, 1).value).toBe(95n);
    expect(C.get(1, 2).value).toBe(5n); // 100 + 6 = 106 ≡ 5 mod 101
  });

  it('should handle large matrix multiplication', () => {
    const A = MS.__call__([
      [10, 20, 30],
      [40, 50, 60],
      [70, 80, 90],
    ]);

    const I = MS.identity();
    const AI = A.mul(I);

    expect(AI.eq(A)).toBe(true);
  });
});

describe('MatrixSpace utilities', () => {
  const F5 = GF(5n);

  it('should report correct dimension', () => {
    const MS23 = MatrixSpace(F5, 2, 3);
    expect(MS23.dimension()).toBe(6);

    const MS44 = MatrixSpace(F5, 4, 4);
    expect(MS44.dimension()).toBe(16);
  });

  it('should report correct dims', () => {
    const MS23 = MatrixSpace(F5, 2, 3);
    expect(MS23.dims()).toEqual([2, 3]);
  });

  it('should check if square', () => {
    const MS22 = MatrixSpace(F5, 2, 2);
    const MS23 = MatrixSpace(F5, 2, 3);

    expect(MS22.is_square()).toBe(true);
    expect(MS23.is_square()).toBe(false);
  });

  it('should generate random matrices', () => {
    const MS = MatrixSpace(F5, 2, 2);
    const R = MS.random_element();

    expect(R.nrows).toBe(2);
    expect(R.ncols).toBe(2);

    // Check all entries are in range [0, 5)
    for (let i = 0; i < 2; i++) {
      for (let j = 0; j < 2; j++) {
        expect(R.get(i, j).value >= 0n).toBe(true);
        expect(R.get(i, j).value < 5n).toBe(true);
      }
    }
  });
});

describe('MatrixFromEntries helper', () => {
  const F7 = GF(7n);

  it('should create matrix from entries', () => {
    const A = MatrixFromEntries(F7, [
      [1, 2, 3],
      [4, 5, 6],
    ]);

    expect(A.nrows).toBe(2);
    expect(A.ncols).toBe(3);
    expect(A.get(0, 0).value).toBe(1n);
    expect(A.get(1, 2).value).toBe(6n);
  });

  it('should handle empty matrix', () => {
    const A = MatrixFromEntries(F7, []);

    expect(A.nrows).toBe(0);
    expect(A.ncols).toBe(0);
  });

  it('should throw on inconsistent row lengths', () => {
    expect(() =>
      MatrixFromEntries(F7, [
        [1, 2, 3],
        [4, 5], // Missing element
      ])
    ).toThrow();
  });
});

describe('Matrix string representation', () => {
  const F7 = GF(7n);
  const MS = MatrixSpace(F7, 2, 2);

  it('should produce readable string output', () => {
    const A = MS.__call__([
      [1, 2],
      [3, 4],
    ]);

    const str = A.toString();
    expect(str).toContain('1');
    expect(str).toContain('2');
    expect(str).toContain('3');
    expect(str).toContain('4');
  });
});
