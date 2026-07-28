/**
 * Tests for matrix operations (determinant, permanent, minors, etc.)
 */

import { describe, expect, it } from 'vitest';
import { GF } from '../rings/finite_rings/finite_field_constructor.js';
import { Zmod } from '../rings/finite_rings/integer_mod_ring.js';
import { QQ } from '../rings/rational_field.js';
import {
  type Matrix,
  MatrixSpace,
  QR,
  adjugate,
  augment,
  charpoly,
  column_space,
  companion_matrix,
  conjugate,
  conjugate_transpose,
  density,
  determinant,
  direct_sum,
  eigenmatrix_right,
  eigenspaces_right,
  eigenvalues,
  eigenvectors_left,
  eigenvectors_right,
  elementary_matrix,
  gram_schmidt_noscale,
  identity_matrix,
  image,
  inverse,
  is_diagonal,
  is_diagonalizable,
  is_nilpotent,
  is_normal,
  is_one,
  is_positive_definite,
  is_positive_semidefinite,
  is_scalar,
  is_semisimple,
  is_similar,
  is_triangular,
  is_unitary,
  kernel,
  left_kernel,
  left_nullity,
  minors,
  minpoly,
  norm,
  permanent,
  permanental_minor,
  pfaffian,
  pseudoinverse,
  quantum_determinant,
  rank,
  right_kernel,
  right_kernel_matrix,
  right_nullity,
  row_space,
  solve_left,
  solve_right,
  stack,
  zero_matrix,
} from './index.js';
import { change_ring, is_hermitian } from './matrix_operations.js';

describe('permanent', () => {
  const F7 = GF(7n);
  const MS2x2 = MatrixSpace(F7, 2, 2);
  const MS3x3 = MatrixSpace(F7, 3, 3);

  it('should compute permanent of 1x1 matrix', () => {
    const A = MS2x2.__call__([
      [5, 0],
      [0, 0],
    ]);
    // Create a 1x1 submatrix
    const B = MatrixSpace(F7, 1, 1).__call__([[5]]);
    expect(permanent(B).value).toBe(5n);
  });

  it('should compute permanent of 2x2 matrix', () => {
    const A = MS2x2.__call__([
      [1, 2],
      [3, 4],
    ]);

    // perm = 1*4 + 2*3 = 4 + 6 = 10 = 3 mod 7
    expect(permanent(A).value).toBe(3n);
  });

  it('should compute permanent of identity matrix', () => {
    const I = MS2x2.identity();
    // perm(I) = 1*1 + 0*0 = 1
    expect(permanent(I).value).toBe(1n);
  });

  it('should compute permanent of ones matrix', () => {
    const A = MS2x2.__call__([
      [1, 1],
      [1, 1],
    ]);
    // perm = 1*1 + 1*1 = 2
    expect(permanent(A).value).toBe(2n);
  });

  it('should give same result with definition algorithm', () => {
    const A = MS2x2.__call__([
      [2, 3],
      [4, 5],
    ]);

    const permRyser = permanent(A, 'Ryser');
    const permDef = permanent(A, 'definition');

    expect(permRyser.eq(permDef)).toBe(true);
  });
});

describe('permanental_minor', () => {
  const F7 = GF(7n);
  const MS2x2 = MatrixSpace(F7, 2, 2);

  it('should compute 0-th permanental minor as 1', () => {
    const A = MS2x2.__call__([
      [1, 2],
      [3, 4],
    ]);

    expect(permanental_minor(A, 0).value).toBe(1n);
  });

  it('should compute 1st permanental minor as sum of entries', () => {
    const A = MS2x2.__call__([
      [1, 2],
      [3, 4],
    ]);

    // Sum of all 1x1 permanents (which are just the entries)
    // 1 + 2 + 3 + 4 = 10 = 3 mod 7
    expect(permanental_minor(A, 1).value).toBe(3n);
  });

  it('should compute 2nd permanental minor as the permanent', () => {
    const A = MS2x2.__call__([
      [1, 2],
      [3, 4],
    ]);

    const perm = permanent(A);
    const pm2 = permanental_minor(A, 2);

    expect(pm2.eq(perm)).toBe(true);
  });
});

describe('minors', () => {
  const F7 = GF(7n);
  const MS2x2 = MatrixSpace(F7, 2, 2);
  const MS3x3 = MatrixSpace(F7, 3, 3);

  it('should return all 1x1 minors (entries)', () => {
    const A = MS2x2.__call__([
      [1, 2],
      [3, 4],
    ]);

    const minorsList = minors(A, 1);
    expect(minorsList.length).toBe(4);

    const values = minorsList.map((m) => m.value);
    expect(values).toContain(1n);
    expect(values).toContain(2n);
    expect(values).toContain(3n);
    expect(values).toContain(4n);
  });

  it('should return determinant for full-size minor', () => {
    const A = MS2x2.__call__([
      [1, 2],
      [3, 4],
    ]);

    const minorsList = minors(A, 2);
    expect(minorsList.length).toBe(1);
    expect(minorsList[0]!.eq(determinant(A))).toBe(true);
  });
});

describe('adjugate', () => {
  const F7 = GF(7n);
  const MS2x2 = MatrixSpace(F7, 2, 2);

  it('should compute adjugate of 2x2 matrix', () => {
    const A = MS2x2.__call__([
      [1, 2],
      [3, 4],
    ]);

    // adj([[a,b],[c,d]]) = [[d,-b],[-c,a]]
    // = [[4, -2], [-3, 1]] = [[4, 5], [4, 1]] mod 7
    const adj = adjugate(A);

    expect(adj.get(0, 0).value).toBe(4n);
    expect(adj.get(0, 1).value).toBe(5n); // -2 mod 7
    expect(adj.get(1, 0).value).toBe(4n); // -3 mod 7
    expect(adj.get(1, 1).value).toBe(1n);
  });

  it('should satisfy adj(A) * A = det(A) * I for invertible matrix', () => {
    const A = MS2x2.__call__([
      [1, 2],
      [3, 4],
    ]);

    const adj = adjugate(A);
    const detA = determinant(A);
    const product = adj.mul(A);
    const detI = MS2x2.identity().scalar_mul(detA);

    expect(product.eq(detI)).toBe(true);
  });
});

describe('is_one', () => {
  const F7 = GF(7n);
  const MS2x2 = MatrixSpace(F7, 2, 2);

  it('should return true for identity matrix', () => {
    const I = MS2x2.identity();
    expect(is_one(I)).toBe(true);
  });

  it('should return false for non-identity matrix', () => {
    const A = MS2x2.__call__([
      [1, 2],
      [0, 1],
    ]);
    expect(is_one(A)).toBe(false);
  });

  it('should return false for non-square matrix', () => {
    const A = MatrixSpace(F7, 2, 3).__call__([
      [1, 0, 0],
      [0, 1, 0],
    ]);
    expect(is_one(A)).toBe(false);
  });
});

describe('is_scalar', () => {
  const F7 = GF(7n);
  const MS2x2 = MatrixSpace(F7, 2, 2);

  it('should return true for scalar matrix', () => {
    const A = MS2x2.__call__([
      [3, 0],
      [0, 3],
    ]);
    expect(is_scalar(A)).toBe(true);
  });

  it('should return true for identity matrix', () => {
    const I = MS2x2.identity();
    expect(is_scalar(I)).toBe(true);
  });

  it('should return true for zero matrix', () => {
    const Z = MS2x2.zero();
    expect(is_scalar(Z)).toBe(true);
  });

  it('should return false for non-scalar matrix', () => {
    const A = MS2x2.__call__([
      [1, 2],
      [0, 1],
    ]);
    expect(is_scalar(A)).toBe(false);
  });

  it('should check against specific scalar', () => {
    const A = MS2x2.__call__([
      [3, 0],
      [0, 3],
    ]);
    expect(is_scalar(A, F7.__call__(3n))).toBe(true);
    expect(is_scalar(A, F7.__call__(2n))).toBe(false);
  });
});

describe('is_diagonal', () => {
  const F7 = GF(7n);
  const MS2x2 = MatrixSpace(F7, 2, 2);

  it('should return true for diagonal matrix', () => {
    const A = MS2x2.__call__([
      [1, 0],
      [0, 2],
    ]);
    expect(is_diagonal(A)).toBe(true);
  });

  it('should return true for identity matrix', () => {
    const I = MS2x2.identity();
    expect(is_diagonal(I)).toBe(true);
  });

  it('should return true for zero matrix', () => {
    const Z = MS2x2.zero();
    expect(is_diagonal(Z)).toBe(true);
  });

  it('should return false for non-diagonal matrix', () => {
    const A = MS2x2.__call__([
      [1, 2],
      [0, 3],
    ]);
    expect(is_diagonal(A)).toBe(false);
  });
});

describe('is_triangular', () => {
  const F7 = GF(7n);
  const MS3x3 = MatrixSpace(F7, 3, 3);

  it('should return true for upper triangular matrix', () => {
    const A = MS3x3.__call__([
      [1, 2, 3],
      [0, 4, 5],
      [0, 0, 6],
    ]);
    expect(is_triangular(A, 'upper')).toBe(true);
    expect(is_triangular(A, 'lower')).toBe(false);
    expect(is_triangular(A)).toBe(true);
  });

  it('should return true for lower triangular matrix', () => {
    const A = MS3x3.__call__([
      [1, 0, 0],
      [2, 3, 0],
      [4, 5, 6],
    ]);
    expect(is_triangular(A, 'lower')).toBe(true);
    expect(is_triangular(A, 'upper')).toBe(false);
    expect(is_triangular(A)).toBe(true);
  });

  it('should return true for diagonal matrix (both upper and lower)', () => {
    const A = MS3x3.__call__([
      [1, 0, 0],
      [0, 2, 0],
      [0, 0, 3],
    ]);
    expect(is_triangular(A, 'upper')).toBe(true);
    expect(is_triangular(A, 'lower')).toBe(true);
    expect(is_triangular(A)).toBe(true);
  });

  it('should return false for non-triangular matrix', () => {
    const A = MS3x3.__call__([
      [1, 2, 3],
      [4, 5, 6],
      [1, 2, 3],
    ]);
    expect(is_triangular(A)).toBe(false);
  });
});

describe('density', () => {
  const F7 = GF(7n);
  const MS2x2 = MatrixSpace(F7, 2, 2);
  const MS3x3 = MatrixSpace(F7, 3, 3);

  // Sage's density() returns an exact rational, not a float:
  //   sage: matrix(QQ, 3,3, [0,1,2,3,0,0,6,7,8]).density()
  //   2/3
  it('should return 0 for zero matrix', () => {
    const Z = MS2x2.zero();
    expect(density(Z).toString()).toBe('0');
  });

  it('should return 1 for matrix with no zeros', () => {
    const A = MS2x2.__call__([
      [1, 2],
      [3, 4],
    ]);
    expect(density(A).toString()).toBe('1');
  });

  it('should return correct fraction for sparse matrix', () => {
    const A = MS2x2.__call__([
      [1, 0],
      [0, 0],
    ]);
    expect(density(A).toString()).toBe('1/4');
  });

  it('should return correct fraction for identity matrix', () => {
    const I = MS3x3.identity();
    // 3 non-zeros out of 9 entries; Sage gives the exact 1/3
    expect(density(I).toString()).toBe('1/3');
  });

  it('should match the Sage doctest value 2/3', () => {
    // sage: matrix(QQ, 3,3, [0,1,2,3,0,0,6,7,8]).density() -> 2/3
    const A = MatrixSpace(QQ, 3, 3).__call__([
      [0, 1, 2],
      [3, 0, 0],
      [6, 7, 8],
    ]);
    expect(density(A).toString()).toBe('2/3');
  });
});

describe('charpoly', () => {
  const F7 = GF(7n);
  const MS2x2 = MatrixSpace(F7, 2, 2);
  const MS3x3 = MatrixSpace(F7, 3, 3);

  it('should compute characteristic polynomial of identity matrix', () => {
    const I = MS2x2.identity();
    const cp = charpoly(I);
    // charpoly(I) = det(xI - I) = (x-1)^2 = x^2 - 2x + 1
    // Coefficients: [1, -2, 1] (constant, x, x^2)
    // In F7: [1, 5, 1] since -2 mod 7 = 5
    expect(cp.coeffs.map((c) => c.value)).toEqual([1n, 5n, 1n]);
  });

  it('should compute characteristic polynomial of 2x2 matrix', () => {
    const A = MS2x2.__call__([
      [1, 2],
      [3, 4],
    ]);
    const cp = charpoly(A);
    // charpoly = x^2 - (a+d)x + (ad-bc)
    // = x^2 - 5x + (4-6) = x^2 - 5x - 2
    // In F7: x^2 + 2x + 5 (since -5 mod 7 = 2 and -2 mod 7 = 5)
    expect(cp.coeffs.map((c) => c.value)).toEqual([5n, 2n, 1n]);
  });

  it('should compute characteristic polynomial of zero matrix', () => {
    const Z = MS2x2.zero();
    const cp = charpoly(Z);
    // charpoly(0) = det(xI) = x^2
    expect(cp.coeffs.map((c) => c.value)).toEqual([0n, 0n, 1n]);
  });

  it('should compute characteristic polynomial of 3x3 matrix', () => {
    const A = MS3x3.__call__([
      [1, 0, 0],
      [0, 2, 0],
      [0, 0, 3],
    ]);
    const cp = charpoly(A);
    // Diagonal matrix: charpoly = (x-1)(x-2)(x-3)
    // = x^3 - 6x^2 + 11x - 6
    // In F7: x^3 + x^2 + 4x + 1 (since -6 mod 7 = 1, 11 mod 7 = 4)
    expect(cp.coeffs.map((c) => c.value)).toEqual([1n, 4n, 1n, 1n]);
  });

  it('should satisfy Cayley-Hamilton theorem', () => {
    // A matrix satisfies its own characteristic polynomial
    const A = MS2x2.__call__([
      [1, 2],
      [3, 4],
    ]);
    const cp = charpoly(A);

    // Evaluate charpoly at A: c0*I + c1*A + c2*A^2
    const I = MS2x2.identity();
    const A2 = A.mul(A);

    // cp.coeffs[0] + cp.coeffs[1]*A + cp.coeffs[2]*A^2 should be zero
    const c0 = cp.coeffs[0]!;
    const c1 = cp.coeffs[1]!;
    const c2 = cp.coeffs[2]!;

    const result = I.scalar_mul(c0).add(A.scalar_mul(c1)).add(A2.scalar_mul(c2));

    // Result should be the zero matrix
    expect(result.eq(MS2x2.zero())).toBe(true);
  });

  it('should compute characteristic polynomial with custom variable name', () => {
    const A = MS2x2.__call__([
      [1, 0],
      [0, 1],
    ]);
    const cp = charpoly(A, 't');
    expect(cp.parent.variable_name).toBe('t');
  });
});

describe('determinant via charpoly', () => {
  const F7 = GF(7n);
  const MS2x2 = MatrixSpace(F7, 2, 2);

  it('determinant equals (-1)^n * constant term of charpoly', () => {
    const A = MS2x2.__call__([
      [1, 2],
      [3, 4],
    ]);
    const det = determinant(A);
    const cp = charpoly(A);
    // det(A) = (-1)^n * cp(0) = (-1)^2 * constant_term = constant_term
    // For n=2, det = constant_term
    expect(det.eq(cp.coeffs[0]!)).toBe(true);
  });
});

describe('QR decomposition', () => {
  const F7 = GF(7n);
  const MS2x2 = MatrixSpace(F7, 2, 2);
  const MS3x3 = MatrixSpace(F7, 3, 3);
  const MS3x2 = MatrixSpace(F7, 3, 2);

  it('should factorize A = Q*R for square matrix', () => {
    const A = MS2x2.__call__([
      [1, 2],
      [3, 4],
    ]);
    const [Q, R] = QR(A);

    // Q * R should equal A
    const product = Q.mul(R);
    expect(product.eq(A)).toBe(true);
  });

  it('should give upper triangular R', () => {
    const A = MS2x2.__call__([
      [1, 2],
      [3, 4],
    ]);
    const [_Q, R] = QR(A);

    // R should be upper triangular (lower entries should be 0)
    for (let i = 0; i < R.nrows; i++) {
      for (let j = 0; j < i && j < R.ncols; j++) {
        expect(R.get(i, j).isZero()).toBe(true);
      }
    }
  });

  it('should produce orthogonal columns in Q', () => {
    const A = MS2x2.__call__([
      [1, 2],
      [3, 4],
    ]);
    const [Q, _R] = QR(A);

    // Check Q columns are orthogonal: col_i . col_j = 0 for i != j
    if (Q.ncols >= 2) {
      const col0 = Q.column(0);
      const col1 = Q.column(1);

      let dot = F7.zero();
      for (let i = 0; i < col0.length; i++) {
        dot = dot.add(col0[i]!.mul(col1[i]!));
      }
      expect(dot.isZero()).toBe(true);
    }
  });

  it('should handle identity matrix', () => {
    const I = MS2x2.identity();
    const [Q, R] = QR(I);

    // Q * R should equal I
    expect(Q.mul(R).eq(I)).toBe(true);
  });

  it('should handle rectangular matrix (tall)', () => {
    const A = MS3x2.__call__([
      [1, 2],
      [3, 4],
      [5, 6],
    ]);
    const [Q, R] = QR(A);

    // Q * R should equal A
    expect(Q.mul(R).eq(A)).toBe(true);
  });
});

describe('gram_schmidt_noscale', () => {
  const F7 = GF(7n);
  const MS3x3 = MatrixSpace(F7, 3, 3);

  it('should produce orthogonal columns', () => {
    const A = MS3x3.__call__([
      [1, 2, 3],
      [4, 5, 6],
      [0, 1, 2],
    ]);
    const [Q, _R] = gram_schmidt_noscale(A);

    // Check all pairs of columns are orthogonal
    for (let i = 0; i < Q.ncols; i++) {
      for (let j = i + 1; j < Q.ncols; j++) {
        const col_i = Q.column(i);
        const col_j = Q.column(j);

        let dot = F7.zero();
        for (let k = 0; k < col_i.length; k++) {
          dot = dot.add(col_i[k]!.mul(col_j[k]!));
        }
        expect(dot.isZero()).toBe(true);
      }
    }
  });

  it('should satisfy A = Q * R', () => {
    const A = MS3x3.__call__([
      [1, 2, 3],
      [4, 5, 6],
      [0, 1, 2],
    ]);
    const [Q, R] = gram_schmidt_noscale(A);

    expect(Q.mul(R).eq(A)).toBe(true);
  });
});

describe('is_nilpotent', () => {
  const F7 = GF(7n);
  const MS2x2 = MatrixSpace(F7, 2, 2);
  const MS3x3 = MatrixSpace(F7, 3, 3);

  it('should return true for zero matrix', () => {
    const Z = MS2x2.zero();
    expect(is_nilpotent(Z)).toBe(true);
  });

  it('should return false for identity matrix', () => {
    const I = MS2x2.identity();
    expect(is_nilpotent(I)).toBe(false);
  });

  it('should return true for strictly upper triangular matrix', () => {
    const A = MS2x2.__call__([
      [0, 1],
      [0, 0],
    ]);
    expect(is_nilpotent(A)).toBe(true);
  });

  it('should return false for non-nilpotent matrix', () => {
    const A = MS2x2.__call__([
      [1, 2],
      [3, 4],
    ]);
    expect(is_nilpotent(A)).toBe(false);
  });

  it('should return true for larger nilpotent matrix', () => {
    // Strictly upper triangular matrices are nilpotent
    const A = MS3x3.__call__([
      [0, 1, 2],
      [0, 0, 3],
      [0, 0, 0],
    ]);
    expect(is_nilpotent(A)).toBe(true);
  });
});

describe('elementary_matrix', () => {
  const F7 = GF(7n);

  it('should create row swap matrix', () => {
    const E = elementary_matrix(F7, 3, { row1: 0, row2: 2 });

    // Swap rows 0 and 2
    expect(E.get(0, 0).value).toBe(0n);
    expect(E.get(0, 2).value).toBe(1n);
    expect(E.get(1, 1).value).toBe(1n);
    expect(E.get(2, 0).value).toBe(1n);
    expect(E.get(2, 2).value).toBe(0n);
  });

  it('should create row scaling matrix', () => {
    const E = elementary_matrix(F7, 3, { row1: 1, scale: F7.__call__(3n) });

    // Scale row 1 by 3
    expect(E.get(0, 0).value).toBe(1n);
    expect(E.get(1, 1).value).toBe(3n);
    expect(E.get(2, 2).value).toBe(1n);
  });

  it('should create row addition matrix', () => {
    const E = elementary_matrix(F7, 3, { row1: 0, row2: 2, scale: F7.__call__(5n) });

    // Add 5 * row 2 to row 0
    expect(E.get(0, 0).value).toBe(1n);
    expect(E.get(0, 2).value).toBe(5n);
    expect(E.get(1, 1).value).toBe(1n);
    expect(E.get(2, 2).value).toBe(1n);
  });

  it('should perform row operation when multiplied', () => {
    const MS3x3 = MatrixSpace(F7, 3, 3);
    const A = MS3x3.__call__([
      [1, 2, 3],
      [4, 5, 6],
      [0, 1, 2],
    ]);

    // Row swap: swap rows 0 and 1
    const E = elementary_matrix(F7, 3, { row1: 0, row2: 1 });
    const result = E.mul(A);

    expect(result.get(0, 0).value).toBe(4n);
    expect(result.get(0, 1).value).toBe(5n);
    expect(result.get(0, 2).value).toBe(6n);
    expect(result.get(1, 0).value).toBe(1n);
    expect(result.get(1, 1).value).toBe(2n);
    expect(result.get(1, 2).value).toBe(3n);
  });
});

describe('companion_matrix', () => {
  const F7 = GF(7n);

  // Sage (matrix/special.py:companion_matrix) takes the *full* coefficient
  // list, low degree first, including the leading 1, and puts the
  // *negatives* of the coefficients along the indicated border:
  //   sage: companion_matrix([-2,-3,-4,-5,-6,1], format='right')
  //   [0 0 0 0 2]
  //   [1 0 0 0 3]
  //   [0 1 0 0 4]
  //   [0 0 1 0 5]
  //   [0 0 0 1 6]
  it('should create companion matrix with correct structure', () => {
    const coeffs = [-2, -3, -4, -5, -6, 1].map((c) => F7.__call__(BigInt(c)));
    const C = companion_matrix(F7, coeffs);

    expect(C.nrows).toBe(5);
    expect(C.ncols).toBe(5);
    for (let i = 0; i < 5; i++) {
      for (let j = 0; j < 4; j++) {
        expect(C.get(i, j).value).toBe(i === j + 1 ? 1n : 0n);
      }
    }
    // last column holds -(-2), -(-3), ... = 2, 3, 4, 5, 6
    expect([0, 1, 2, 3, 4].map((i) => C.get(i, 4).value)).toEqual([2n, 3n, 4n, 5n, 6n]);
  });

  it('should reject a non-monic coefficient list', () => {
    expect(() => companion_matrix(F7, [F7.__call__(4n), F7.__call__(3n), F7.__call__(2n)])).toThrow(
      /must be monic/
    );
  });

  it('should have characteristic polynomial matching input coefficients', () => {
    // The companion matrix of a monic polynomial has that polynomial as its
    // characteristic polynomial.
    const coeffs = [2, 3, 1].map((c) => F7.__call__(BigInt(c)));
    const C = companion_matrix(F7, coeffs);
    const cp = charpoly(C);

    expect(cp.coeffs.length).toBe(3);
    expect(cp.coeffs[0]!.value).toBe(2n);
    expect(cp.coeffs[1]!.value).toBe(3n);
    expect(cp.coeffs[2]!.value).toBe(1n);
  });
});

describe('eigenvalues', () => {
  const F7 = GF(7n);
  const MS2x2 = MatrixSpace(F7, 2, 2);
  const MS3x3 = MatrixSpace(F7, 3, 3);

  it('should find eigenvalues of diagonal matrix', () => {
    const A = MS2x2.__call__([
      [2, 0],
      [0, 5],
    ]);
    const evs = eigenvalues(A);

    // Should find 2 and 5
    expect(evs.length).toBe(2);
    const values = evs.map((e) => e.value);
    expect(values).toContain(2n);
    expect(values).toContain(5n);
  });

  it('should find eigenvalue of scalar matrix', () => {
    const A = MS2x2.__call__([
      [3, 0],
      [0, 3],
    ]);
    const evs = eigenvalues(A);

    // Should find 3 (possibly twice if we count multiplicities)
    const values = evs.map((e) => e.value);
    expect(values).toContain(3n);
  });

  it('should find eigenvalues of identity matrix', () => {
    const I = MS2x2.identity();
    const evs = eigenvalues(I);

    // Identity matrix has eigenvalue 1
    expect(evs.some((e) => e.value === 1n)).toBe(true);
  });

  it('should find no eigenvalues for zero matrix except 0', () => {
    const Z = MS2x2.zero();
    const evs = eigenvalues(Z);

    // Zero matrix has eigenvalue 0
    expect(evs.every((e) => e.value === 0n)).toBe(true);
  });

  it('should work with 3x3 diagonal matrix', () => {
    const A = MS3x3.__call__([
      [1, 0, 0],
      [0, 2, 0],
      [0, 0, 4],
    ]);
    const evs = eigenvalues(A);

    const values = evs.map((e) => e.value);
    expect(values).toContain(1n);
    expect(values).toContain(2n);
    expect(values).toContain(4n);
  });
});

describe('eigenvectors_right', () => {
  const F7 = GF(7n);
  const MS2x2 = MatrixSpace(F7, 2, 2);
  const MS3x3 = MatrixSpace(F7, 3, 3);

  it('should find eigenvectors of diagonal matrix', () => {
    const A = MS2x2.__call__([
      [2, 0],
      [0, 5],
    ]);
    const result = eigenvectors_right(A);

    // Should have 2 eigenvalue-eigenvector pairs
    expect(result.length).toBe(2);

    // Verify A * v = lambda * v for each eigenvector
    for (const [lambda, vecs, _mult] of result) {
      for (const v of vecs) {
        // Compute A * v
        const Av = [
          A.get(0, 0).mul(v[0]!).add(A.get(0, 1).mul(v[1]!)),
          A.get(1, 0).mul(v[0]!).add(A.get(1, 1).mul(v[1]!)),
        ];

        // Check A * v = lambda * v
        expect(Av[0]!.eq(lambda.mul(v[0]!))).toBe(true);
        expect(Av[1]!.eq(lambda.mul(v[1]!))).toBe(true);
      }
    }
  });

  it('should find eigenvectors of identity matrix', () => {
    const I = MS2x2.identity();
    const result = eigenvectors_right(I);

    // Identity has eigenvalue 1 with algebraic multiplicity 2
    expect(result.length).toBeGreaterThan(0);
    expect(result[0]![0].value).toBe(1n); // eigenvalue is 1
  });

  it('should find eigenvectors of scalar matrix', () => {
    const A = MS2x2.__call__([
      [4, 0],
      [0, 4],
    ]);
    const result = eigenvectors_right(A);

    // Should have eigenvalue 4
    expect(result.some(([ev, _vecs, _m]) => ev.value === 4n)).toBe(true);
  });
});

describe('eigenvectors_left', () => {
  const F7 = GF(7n);
  const MS2x2 = MatrixSpace(F7, 2, 2);

  it('should find left eigenvectors of diagonal matrix', () => {
    const A = MS2x2.__call__([
      [2, 0],
      [0, 5],
    ]);
    const result = eigenvectors_left(A);

    // Should have 2 eigenvalue-eigenvector pairs
    expect(result.length).toBe(2);

    // Verify v * A = lambda * v for each left eigenvector
    for (const [lambda, vecs, _mult] of result) {
      for (const v of vecs) {
        // Compute v * A (v is a row vector)
        const vA = [
          v[0]!.mul(A.get(0, 0)).add(v[1]!.mul(A.get(1, 0))),
          v[0]!.mul(A.get(0, 1)).add(v[1]!.mul(A.get(1, 1))),
        ];

        // Check v * A = lambda * v
        expect(vA[0]!.eq(lambda.mul(v[0]!))).toBe(true);
        expect(vA[1]!.eq(lambda.mul(v[1]!))).toBe(true);
      }
    }
  });
});

describe('eigenspaces_right', () => {
  const F7 = GF(7n);
  const MS2x2 = MatrixSpace(F7, 2, 2);

  it('should return eigenspace objects', () => {
    const A = MS2x2.__call__([
      [2, 0],
      [0, 5],
    ]);
    const spaces = eigenspaces_right(A);

    expect(spaces.length).toBe(2);

    for (const space of spaces) {
      expect(space.eigenvalue).toBeDefined();
      expect(space.basis).toBeDefined();
      expect(space.algebraicMultiplicity).toBeGreaterThan(0);
      expect(space.geometricMultiplicity).toBeGreaterThan(0);
    }
  });
});

describe('eigenmatrix_right', () => {
  const F7 = GF(7n);
  const MS2x2 = MatrixSpace(F7, 2, 2);

  it('should return diagonal and eigenvector matrices', () => {
    const A = MS2x2.__call__([
      [2, 0],
      [0, 5],
    ]);
    const [D, P] = eigenmatrix_right(A);

    // D should be diagonal
    expect(is_diagonal(D)).toBe(true);

    // P * A = P * D for diagonal matrices
    // Since A is already diagonal, this is simpler
    expect(D.nrows).toBe(2);
    expect(D.ncols).toBe(2);
    expect(P.nrows).toBe(2);
    expect(P.ncols).toBe(2);
  });
});

describe('is_unitary', () => {
  const F7 = GF(7n);
  const MS2x2 = MatrixSpace(F7, 2, 2);

  it('should return true for identity matrix', () => {
    const I = MS2x2.identity();
    expect(is_unitary(I)).toBe(true);
  });

  it('should return false for non-unitary matrix', () => {
    const A = MS2x2.__call__([
      [1, 2],
      [3, 4],
    ]);
    expect(is_unitary(A)).toBe(false);
  });

  it('should return false for non-square matrix', () => {
    const A = MatrixSpace(F7, 2, 3).__call__([
      [1, 0, 0],
      [0, 1, 0],
    ]);
    expect(is_unitary(A)).toBe(false);
  });
});

describe('is_normal', () => {
  const F7 = GF(7n);
  const MS2x2 = MatrixSpace(F7, 2, 2);

  it('should return true for diagonal matrix', () => {
    const A = MS2x2.__call__([
      [2, 0],
      [0, 5],
    ]);
    expect(is_normal(A)).toBe(true);
  });

  it('should return true for symmetric matrix', () => {
    const A = MS2x2.__call__([
      [1, 2],
      [2, 3],
    ]);
    expect(is_normal(A)).toBe(true);
  });

  it('should return true for identity matrix', () => {
    const I = MS2x2.identity();
    expect(is_normal(I)).toBe(true);
  });

  it('should return true for zero matrix', () => {
    const Z = MS2x2.zero();
    expect(is_normal(Z)).toBe(true);
  });
});

describe('is_positive_definite', () => {
  // Sage convention (matrix2.pyx:_is_positive_definite_or_semidefinite):
  // definiteness only makes sense over a subring of the reals/complexes, so a
  // finite field raises ValueError rather than reporting a boolean.
  const F7 = GF(7n);
  const MS2x2 = MatrixSpace(F7, 2, 2);
  const qq = (rows: (number | string)[][]) =>
    MatrixSpace(QQ, rows.length, rows[0]!.length).__call__(rows);

  it('should return true for the identity matrix over QQ', () => {
    expect(
      is_positive_definite(
        qq([
          [1, 0],
          [0, 1],
        ])
      )
    ).toBe(true);
  });

  it('should return false for the zero matrix over QQ', () => {
    expect(
      is_positive_definite(
        qq([
          [0, 0],
          [0, 0],
        ])
      )
    ).toBe(false);
  });

  it('should return false for a non-symmetric matrix', () => {
    expect(
      is_positive_definite(
        qq([
          [1, 2],
          [3, 4],
        ])
      )
    ).toBe(false);
  });

  it('should return false for a singular symmetric matrix', () => {
    expect(
      is_positive_definite(
        qq([
          [1, 1],
          [1, 1],
        ])
      )
    ).toBe(false);
  });

  it('should return false for negative definite matrices', () => {
    // leading principal minors -1 and -2 < 0
    expect(is_positive_definite(qq([[-1]]))).toBe(false);
    expect(
      is_positive_definite(
        qq([
          [-2, -1],
          [-1, -2],
        ])
      )
    ).toBe(false);
  });

  it('should match the Sage doctest for a 4x4 positive definite matrix', () => {
    // sage: A = matrix(QQ, [[4,-2,4,2],[-2,10,-2,-7],[4,-2,8,4],[2,-7,4,7]])
    // sage: A.is_positive_definite()
    // True   (leading principal minors 4, 36, 144, 144)
    const A = qq([
      [4, -2, 4, 2],
      [-2, 10, -2, -7],
      [4, -2, 8, 4],
      [2, -7, 4, 7],
    ]);
    expect(is_positive_definite(A)).toBe(true);
  });

  it('should match the Sage doctests for indefinite matrices', () => {
    // sage: matrix(QQ, [[3,-6,9,6,-9],[-6,11,-16,-11,17],[9,-16,28,16,-40],
    //                   [6,-11,16,9,-19],[-9,17,-40,-19,68]]).is_positive_definite()
    // False
    expect(
      is_positive_definite(
        qq([
          [3, -6, 9, 6, -9],
          [-6, 11, -16, -11, 17],
          [9, -16, 28, 16, -40],
          [6, -11, 16, 9, -19],
          [-9, 17, -40, -19, 68],
        ])
      )
    ).toBe(false);
    // sage: matrix(QQ, [[21,15,12,-2],[15,12,9,6],[12,9,7,3],[-2,6,3,8]])
    //         .is_positive_definite()
    // False
    expect(
      is_positive_definite(
        qq([
          [21, 15, 12, -2],
          [15, 12, 9, 6],
          [12, 9, 7, 3],
          [-2, 6, 3, 8],
        ])
      )
    ).toBe(false);
  });

  it('should raise over a finite field', () => {
    expect(() => is_positive_definite(MS2x2.identity())).toThrow(
      /as a subring of the real or complex numbers/
    );
  });
});

describe('is_positive_semidefinite', () => {
  const qq = (rows: (number | string)[][]) =>
    MatrixSpace(QQ, rows.length, rows[0]!.length).__call__(rows);

  it('should match the Sage doctests', () => {
    // sage: matrix(QQ, [[1,1],[1,1]]).is_positive_semidefinite()  -> True
    expect(
      is_positive_semidefinite(
        qq([
          [1, 1],
          [1, 1],
        ])
      )
    ).toBe(true);
    // sage: matrix(QQ, [[0,1],[1,0]]).is_positive_semidefinite()  -> False
    expect(
      is_positive_semidefinite(
        qq([
          [0, 1],
          [1, 0],
        ])
      )
    ).toBe(false);
    // sage: matrix(QQ, [[2,1],[0,0]]).is_positive_semidefinite()  -> False (not Hermitian)
    expect(
      is_positive_semidefinite(
        qq([
          [2, 1],
          [0, 0],
        ])
      )
    ).toBe(false);
    // sage: matrix(QQ, 0).is_positive_semidefinite()              -> True (vacuous)
    expect(is_positive_semidefinite(MatrixSpace(QQ, 0, 0).zero())).toBe(true);
  });

  it('should return false for a negative definite matrix', () => {
    expect(is_positive_semidefinite(qq([[-1]]))).toBe(false);
  });

  it('should agree with is_positive_definite on strictly definite matrices', () => {
    const A = qq([
      [4, -2, 4, 2],
      [-2, 10, -2, -7],
      [4, -2, 8, 4],
      [2, -7, 4, 7],
    ]);
    expect(is_positive_definite(A)).toBe(true);
    expect(is_positive_semidefinite(A)).toBe(true);
  });

  it('should raise over a finite field', () => {
    expect(() => is_positive_semidefinite(MatrixSpace(GF(7n), 1, 1).identity())).toThrow(
      /as a subring of the real or complex numbers/
    );
  });
});

describe('is_similar', () => {
  const F7 = GF(7n);
  const MS2x2 = MatrixSpace(F7, 2, 2);

  it('should return true for identical matrices', () => {
    const A = MS2x2.__call__([
      [1, 2],
      [3, 4],
    ]);
    expect(is_similar(A, A)).toBe(true);
  });

  it('should return true for identical diagonal matrices', () => {
    // Two identical diagonal matrices are trivially similar
    const A = MS2x2.__call__([
      [2, 0],
      [0, 5],
    ]);
    const B = MS2x2.__call__([
      [2, 0],
      [0, 5],
    ]);
    expect(is_similar(A, B)).toBe(true);
  });

  it('should return false for matrices with different characteristic polynomials', () => {
    const A = MS2x2.__call__([
      [1, 0],
      [0, 2],
    ]);
    const B = MS2x2.__call__([
      [1, 0],
      [0, 3],
    ]);
    expect(is_similar(A, B)).toBe(false);
  });
});

describe('conjugate and conjugate_transpose', () => {
  const F7 = GF(7n);
  const MS2x2 = MatrixSpace(F7, 2, 2);

  it('conjugate should return copy for real fields', () => {
    const A = MS2x2.__call__([
      [1, 2],
      [3, 4],
    ]);
    const conj = conjugate(A);
    expect(conj.eq(A)).toBe(true);
  });

  it('conjugate_transpose should equal transpose for real fields', () => {
    const A = MS2x2.__call__([
      [1, 2],
      [3, 4],
    ]);
    const ct = conjugate_transpose(A);
    const t = A.transpose();
    expect(ct.eq(t)).toBe(true);
  });
});

describe('norm', () => {
  // Sage's norm is `self.apply_map(abs, R=RDF)` followed by a row/column max
  // or a square root, and returns an RDF number (matrix2.pyx:16379).
  const qq = (rows: number[][]) => MatrixSpace(QQ, rows.length, rows[0]!.length).__call__(rows);

  it('should return zero for the zero matrix', () => {
    expect(
      norm(
        qq([
          [0, 0],
          [0, 0],
        ]),
        'frob'
      )
    ).toBe(0);
  });

  it('should compute the three norms of [[1,2],[3,4]]', () => {
    const A = qq([
      [1, 2],
      [3, 4],
    ]);
    expect(norm(A, 1)).toBe(6); // largest column sum |1|+|3|=4, |2|+|4|=6
    expect(norm(A, Number.POSITIVE_INFINITY)).toBe(7); // largest row sum 3+4
    expect(norm(A, 'frob')).toBeCloseTo(Math.sqrt(30), 12); // 5.477225575051661
  });

  it('should match the Sage doctest for norm(1) and norm(Infinity)', () => {
    // sage: A = matrix(ZZ, [[1,2,4,3], [-1,0,3,-10]])
    // sage: A.norm(1)          -> 13.0
    // sage: A.norm(Infinity)   -> 14.0
    const A = qq([
      [1, 2, 4, 3],
      [-1, 0, 3, -10],
    ]);
    expect(norm(A, 1)).toBe(13);
    expect(norm(A, Number.POSITIVE_INFINITY)).toBe(14);
  });

  it('should satisfy norm(A, Infinity) === norm(A^T, 1)', () => {
    const A = qq([
      [3, -1, 0],
      [2, 5, -7],
    ]);
    expect(norm(A, Number.POSITIVE_INFINITY)).toBe(norm(A.transpose(), 1));
  });

  it('should return 0 for an empty matrix', () => {
    expect(norm(MatrixSpace(QQ, 0, 0).zero(), 1)).toBe(0);
  });

  it('should reject matrices whose entries have no absolute value', () => {
    const A = MatrixSpace(GF(7n), 2, 2).__call__([
      [1, 2],
      [3, 4],
    ]);
    expect(() => norm(A, 1)).toThrow(/bad operand type for abs/);
  });

  // This test used to assert that the default (p = 2) norm threw, because the
  // port had no SVD.  It is now computed exactly (see the
  // "norm(A, 2) — largest singular value" block below), so the assertion is
  // replaced by Sage's own value rather than being deleted.
  it('should default to p = 2, the largest singular value', () => {
    // sage: matrix(ZZ, [[1, 2], [3, 4]], sparse=True).norm()
    // 5.464985704219043
    const A = qq([
      [1, 2],
      [3, 4],
    ]);
    expect(norm(A)).toBe(5.464985704219043);
    expect(norm(A, 2)).toBe(norm(A));
  });
});

describe('direct_sum', () => {
  const F7 = GF(7n);
  const MS2x2 = MatrixSpace(F7, 2, 2);
  const MS3x3 = MatrixSpace(F7, 3, 3);

  it('should create block diagonal matrix from two matrices', () => {
    const A = MS2x2.__call__([
      [1, 2],
      [3, 4],
    ]);
    const B = MS2x2.__call__([
      [5, 6],
      [0, 1],
    ]);
    const C = direct_sum(A, B);

    expect(C.nrows).toBe(4);
    expect(C.ncols).toBe(4);

    // Check upper-left block (A)
    expect(C.get(0, 0).value).toBe(1n);
    expect(C.get(0, 1).value).toBe(2n);
    expect(C.get(1, 0).value).toBe(3n);
    expect(C.get(1, 1).value).toBe(4n);

    // Check lower-right block (B)
    expect(C.get(2, 2).value).toBe(5n);
    expect(C.get(2, 3).value).toBe(6n);
    expect(C.get(3, 2).value).toBe(0n);
    expect(C.get(3, 3).value).toBe(1n);

    // Check off-diagonal blocks are zero
    expect(C.get(0, 2).value).toBe(0n);
    expect(C.get(0, 3).value).toBe(0n);
    expect(C.get(1, 2).value).toBe(0n);
    expect(C.get(1, 3).value).toBe(0n);
    expect(C.get(2, 0).value).toBe(0n);
    expect(C.get(2, 1).value).toBe(0n);
    expect(C.get(3, 0).value).toBe(0n);
    expect(C.get(3, 1).value).toBe(0n);
  });

  it('should handle matrices of different sizes', () => {
    const A = MS2x2.__call__([
      [1, 0],
      [0, 1],
    ]);
    const B = MS3x3.__call__([
      [2, 0, 0],
      [0, 2, 0],
      [0, 0, 2],
    ]);
    const C = direct_sum(A, B);

    expect(C.nrows).toBe(5);
    expect(C.ncols).toBe(5);
    expect(C.get(0, 0).value).toBe(1n);
    expect(C.get(2, 2).value).toBe(2n);
    expect(C.get(4, 4).value).toBe(2n);
  });
});

describe('stack', () => {
  const F7 = GF(7n);
  const MS2x3 = MatrixSpace(F7, 2, 3);
  const MS3x3 = MatrixSpace(F7, 3, 3);

  it('should stack two matrices vertically', () => {
    const A = MS2x3.__call__([
      [1, 2, 3],
      [4, 5, 6],
    ]);
    const B = MS2x3.__call__([
      [0, 1, 0],
      [1, 0, 1],
    ]);
    const C = stack(A, B);

    expect(C.nrows).toBe(4);
    expect(C.ncols).toBe(3);

    // Check first block (A)
    expect(C.get(0, 0).value).toBe(1n);
    expect(C.get(0, 2).value).toBe(3n);
    expect(C.get(1, 1).value).toBe(5n);

    // Check second block (B)
    expect(C.get(2, 0).value).toBe(0n);
    expect(C.get(2, 1).value).toBe(1n);
    expect(C.get(3, 2).value).toBe(1n);
  });

  it('should stack multiple matrices', () => {
    const A = MS2x3.__call__([
      [1, 1, 1],
      [1, 1, 1],
    ]);
    const B = MS2x3.__call__([
      [2, 2, 2],
      [2, 2, 2],
    ]);
    const D = MS2x3.__call__([
      [3, 3, 3],
      [3, 3, 3],
    ]);
    const C = stack(A, [B, D]);

    expect(C.nrows).toBe(6);
    expect(C.ncols).toBe(3);
    expect(C.get(0, 0).value).toBe(1n);
    expect(C.get(2, 0).value).toBe(2n);
    expect(C.get(4, 0).value).toBe(3n);
  });

  it('should throw error for incompatible column counts', () => {
    const A = MS2x3.__call__([
      [1, 2, 3],
      [4, 5, 6],
    ]);
    const B = MatrixSpace(F7, 2, 2).__call__([
      [1, 2],
      [3, 4],
    ]);
    expect(() => stack(A, B)).toThrow();
  });
});

describe('augment', () => {
  const F7 = GF(7n);
  const MS2x2 = MatrixSpace(F7, 2, 2);
  const MS2x3 = MatrixSpace(F7, 2, 3);

  it('should augment two matrices horizontally', () => {
    const A = MS2x2.__call__([
      [1, 2],
      [3, 4],
    ]);
    const B = MS2x2.__call__([
      [5, 6],
      [0, 1],
    ]);
    const C = augment(A, B);

    expect(C.nrows).toBe(2);
    expect(C.ncols).toBe(4);

    // Check first block (A)
    expect(C.get(0, 0).value).toBe(1n);
    expect(C.get(0, 1).value).toBe(2n);
    expect(C.get(1, 0).value).toBe(3n);
    expect(C.get(1, 1).value).toBe(4n);

    // Check second block (B)
    expect(C.get(0, 2).value).toBe(5n);
    expect(C.get(0, 3).value).toBe(6n);
    expect(C.get(1, 2).value).toBe(0n);
    expect(C.get(1, 3).value).toBe(1n);
  });

  it('should augment multiple matrices', () => {
    const A = MS2x2.__call__([
      [1, 1],
      [1, 1],
    ]);
    const B = MS2x2.__call__([
      [2, 2],
      [2, 2],
    ]);
    const D = MS2x2.__call__([
      [3, 3],
      [3, 3],
    ]);
    const C = augment(A, [B, D]);

    expect(C.nrows).toBe(2);
    expect(C.ncols).toBe(6);
    expect(C.get(0, 0).value).toBe(1n);
    expect(C.get(0, 2).value).toBe(2n);
    expect(C.get(0, 4).value).toBe(3n);
  });

  it('should throw error for incompatible row counts', () => {
    const A = MS2x2.__call__([
      [1, 2],
      [3, 4],
    ]);
    const B = MS2x3.__call__([
      [1, 2, 3],
      [4, 5, 6],
    ]);
    // This is compatible
    const C = augment(A, B);
    expect(C.ncols).toBe(5);

    // But this should throw
    const D = MatrixSpace(F7, 3, 2).__call__([
      [1, 2],
      [3, 4],
      [5, 6],
    ]);
    expect(() => augment(A, D)).toThrow();
  });
});

describe('right_kernel', () => {
  const F7 = GF(7n);
  const MS2x3 = MatrixSpace(F7, 2, 3);
  const MS3x3 = MatrixSpace(F7, 3, 3);

  it('should return a vector space', () => {
    const A = MS2x3.__call__([
      [1, 2, 3],
      [4, 5, 6],
    ]);
    const K = right_kernel(A);

    // Should be a module/vector space
    expect(K).toBeDefined();
    expect(typeof K.degree).toBe('function');
    expect(K.degree()).toBe(3); // Kernel vectors have length = ncols
  });

  it('should contain vectors that multiply to zero', () => {
    const A = MS2x3.__call__([
      [1, 2, 3],
      [4, 5, 6],
    ]);
    const K = right_kernel(A);
    const basis = K.gens();

    // Each basis vector v should satisfy A * v = 0
    for (const v of basis) {
      const entries = v.list();
      // Compute A * v
      for (let i = 0; i < A.nrows; i++) {
        let sum = F7.zero();
        for (let j = 0; j < A.ncols; j++) {
          sum = sum.add(A.get(i, j).mul(entries[j]!));
        }
        expect(sum.isZero()).toBe(true);
      }
    }
  });

  it('should have correct dimension for full rank matrix', () => {
    // Full rank 2x3 matrix: kernel dimension = 3 - 2 = 1
    const A = MS2x3.__call__([
      [1, 0, 1],
      [0, 1, 1],
    ]);
    const K = right_kernel(A);

    expect(K.rank()).toBe(1);
  });

  it('should have trivial kernel for full column rank matrix', () => {
    // 3x3 identity has trivial kernel
    const I = MS3x3.identity();
    const K = right_kernel(I);

    expect(K.rank()).toBe(0);
  });

  it('should have full kernel for zero matrix', () => {
    const Z = MS2x3.zero();
    const K = right_kernel(Z);

    // Kernel of zero matrix is entire ambient space
    expect(K.rank()).toBe(3);
  });
});

describe('left_kernel', () => {
  const F7 = GF(7n);
  const MS3x2 = MatrixSpace(F7, 3, 2);
  const MS3x3 = MatrixSpace(F7, 3, 3);

  it('should return a vector space', () => {
    const A = MS3x2.__call__([
      [1, 2],
      [3, 4],
      [5, 6],
    ]);
    const K = left_kernel(A);

    expect(K).toBeDefined();
    expect(K.degree()).toBe(3); // Kernel vectors have length = nrows
  });

  it('should contain vectors that multiply to zero from left', () => {
    const A = MS3x2.__call__([
      [1, 2],
      [3, 4],
      [5, 6],
    ]);
    const K = left_kernel(A);
    const basis = K.gens();

    // Each basis vector v should satisfy v * A = 0
    for (const v of basis) {
      const entries = v.list();
      // Compute v * A
      for (let j = 0; j < A.ncols; j++) {
        let sum = F7.zero();
        for (let i = 0; i < A.nrows; i++) {
          sum = sum.add(entries[i]!.mul(A.get(i, j)));
        }
        expect(sum.isZero()).toBe(true);
      }
    }
  });

  it('should have correct dimension', () => {
    // 3x2 matrix with rank 2: left kernel dimension = 3 - 2 = 1
    const A = MS3x2.__call__([
      [1, 0],
      [0, 1],
      [1, 1],
    ]);
    const K = left_kernel(A);

    expect(K.rank()).toBe(1);
  });

  it('should have trivial kernel for full row rank matrix', () => {
    const I = MS3x3.identity();
    const K = left_kernel(I);

    expect(K.rank()).toBe(0);
  });
});

describe('image', () => {
  const F7 = GF(7n);
  const MS3x3 = MatrixSpace(F7, 3, 3);

  it('should return the row space', () => {
    const A = MS3x3.__call__([
      [1, 2, 3],
      [4, 5, 6],
      [5, 0, 2],
    ]);
    const img = image(A);
    const rs = row_space(A);

    // image and row_space should be the same
    expect(img.degree()).toBe(rs.degree());
    expect(img.rank()).toBe(rs.rank());
  });

  it('should have correct dimension', () => {
    const A = MS3x3.__call__([
      [1, 0, 0],
      [0, 1, 0],
      [0, 0, 0],
    ]);
    const img = image(A);

    expect(img.rank()).toBe(2);
  });
});

describe('row_space', () => {
  const F7 = GF(7n);
  const MS2x3 = MatrixSpace(F7, 2, 3);
  const MS3x3 = MatrixSpace(F7, 3, 3);

  it('should return a vector space with correct degree', () => {
    const A = MS2x3.__call__([
      [1, 2, 3],
      [4, 5, 6],
    ]);
    const V = row_space(A);

    expect(V.degree()).toBe(3); // Degree = number of columns
  });

  it('should have dimension equal to rank', () => {
    const A = MS3x3.__call__([
      [1, 0, 0],
      [0, 1, 0],
      [1, 1, 0],
    ]);
    const V = row_space(A);

    expect(V.rank()).toBe(2); // rank of A is 2
  });

  it('should span all rows', () => {
    const A = MS2x3.__call__([
      [1, 0, 1],
      [0, 1, 1],
    ]);
    const V = row_space(A);

    // Should have rank 2 since rows are independent
    expect(V.rank()).toBe(2);
  });

  it('should return zero space for zero matrix', () => {
    const Z = MS2x3.zero();
    const V = row_space(Z);

    expect(V.rank()).toBe(0);
  });

  it('should return full space for identity matrix', () => {
    const I = MS3x3.identity();
    const V = row_space(I);

    expect(V.rank()).toBe(3);
    expect(V.degree()).toBe(3);
  });
});

describe('column_space', () => {
  const F7 = GF(7n);
  const MS3x2 = MatrixSpace(F7, 3, 2);
  const MS3x3 = MatrixSpace(F7, 3, 3);

  it('should return a vector space with correct degree', () => {
    const A = MS3x2.__call__([
      [1, 2],
      [3, 4],
      [5, 6],
    ]);
    const V = column_space(A);

    expect(V.degree()).toBe(3); // Degree = number of rows
  });

  it('should have dimension equal to rank', () => {
    const A = MS3x3.__call__([
      [1, 0, 1],
      [0, 1, 1],
      [0, 0, 0],
    ]);
    const V = column_space(A);

    expect(V.rank()).toBe(2); // rank of A is 2
  });

  it('should equal row space of transpose', () => {
    const A = MS3x2.__call__([
      [1, 2],
      [3, 4],
      [5, 6],
    ]);
    const colSpace = column_space(A);
    const rowSpaceT = row_space(A.transpose());

    expect(colSpace.degree()).toBe(rowSpaceT.degree());
    expect(colSpace.rank()).toBe(rowSpaceT.rank());
  });

  it('should return zero space for zero matrix', () => {
    const Z = MS3x2.zero();
    const V = column_space(Z);

    expect(V.rank()).toBe(0);
  });

  it('should return full space for identity matrix', () => {
    const I = MS3x3.identity();
    const V = column_space(I);

    expect(V.rank()).toBe(3);
    expect(V.degree()).toBe(3);
  });
});

describe('kernel (alias for left_kernel)', () => {
  const F7 = GF(7n);
  const MS3x3 = MatrixSpace(F7, 3, 3);

  it('should be the same as left_kernel', () => {
    const A = MS3x3.__call__([
      [1, 2, 3],
      [4, 5, 6],
      [7, 8, 9],
    ]);
    const K1 = kernel(A);
    const K2 = left_kernel(A);

    expect(K1.degree()).toBe(K2.degree());
    expect(K1.rank()).toBe(K2.rank());
  });
});

// ============================================================================
// Regression tests for the 2026-07 audit (C8, H46-H51, M55-M58, L42, L43, M56)
// ============================================================================

/** Evaluate a polynomial at a square matrix (Horner). */
function polyAtMatrix<R extends { add(o: R): R; mul(o: R): R; isZero(): boolean }>(
  poly: { degree(): number; getCoeff(i: number): R },
  A: Matrix<R>
): Matrix<R> {
  const n = A.nrows;
  const d = poly.degree();
  if (d < 0) return zero_matrix(A.base_ring, n);
  let res = identity_matrix(A.base_ring, n).scalar_mul(poly.getCoeff(d));
  for (let i = d - 1; i >= 0; i--) {
    res = res.mul(A).add(identity_matrix(A.base_ring, n).scalar_mul(poly.getCoeff(i)));
  }
  return res;
}

/** Small deterministic LCG so the sweeps below are reproducible. */
function makeRng(seed: number): (m: number) => number {
  let s = seed;
  return (m: number) => {
    s = (s * 1103515245 + 12345) % 2147483648;
    return s % m;
  };
}

describe('minpoly (C8)', () => {
  const F7 = GF(7n);

  it('should return the minimal polynomial of diag(1,2), not of e_0', () => {
    // Sage: matrix(GF(7), [[1,0],[0,2]]).minpoly() == x^2 + 4*x + 2
    const A = MatrixSpace(F7, 2, 2).__call__([
      [1, 0],
      [0, 2],
    ]);
    const mp = minpoly(A);
    expect(mp.coeffs.map((c) => c.value)).toEqual([2n, 4n, 1n]);
    expect(polyAtMatrix(mp, A).is_zero()).toBe(true);
  });

  it('should return (x-1)^2 for the Jordan block [[1,1],[0,1]]', () => {
    const A = MatrixSpace(F7, 2, 2).__call__([
      [1, 1],
      [0, 1],
    ]);
    // (x-1)^2 = x^2 - 2x + 1 = x^2 + 5x + 1 over GF(7)
    expect(minpoly(A).coeffs.map((c) => c.value)).toEqual([1n, 5n, 1n]);
  });

  it('should return (x-1)(x-2)(x-3) for diag(1,2,3)', () => {
    const A = MatrixSpace(F7, 3, 3).__call__([
      [1, 0, 0],
      [0, 2, 0],
      [0, 0, 3],
    ]);
    // x^3 - 6x^2 + 11x - 6 = x^3 + x^2 + 4x + 1 over GF(7)
    expect(minpoly(A).coeffs.map((c) => c.value)).toEqual([1n, 4n, 1n, 1n]);
  });

  it('should annihilate the matrix for 50 random matrices over GF(7)', () => {
    const rnd = makeRng(20260727);
    for (let t = 0; t < 50; t++) {
      const n = 2 + rnd(3);
      const rows: number[][] = [];
      for (let i = 0; i < n; i++) {
        const r: number[] = [];
        for (let j = 0; j < n; j++) r.push(rnd(7));
        rows.push(r);
      }
      const A = MatrixSpace(F7, n, n).__call__(rows);
      expect(polyAtMatrix(minpoly(A), A).is_zero()).toBe(true);
    }
  });

  it('should annihilate the matrix over QQ (Cayley-Hamilton)', () => {
    const A = MatrixSpace(QQ, 3, 3).__call__([
      [1, 2, 3],
      [4, 5, 6],
      [7, 8, 10],
    ]);
    expect(polyAtMatrix(minpoly(A), A).is_zero()).toBe(true);
    expect(polyAtMatrix(charpoly(A), A).is_zero()).toBe(true);
  });

  it('should equal the charpoly when the charpoly is squarefree', () => {
    const A = MatrixSpace(F7, 2, 2).__call__([
      [1, 2],
      [3, 4],
    ]);
    expect(minpoly(A).toString()).toBe(charpoly(A).toString());
  });
});

describe('determinant over non-fields (H51)', () => {
  const R8 = Zmod(8);

  it('should compute a 4x4 determinant over Z/8 without dividing', () => {
    // det over ZZ of [[1,2,3,4],[5,6,7,0],[1,3,5,7],[2,4,6,1]] is -240; -240 mod 8 = 0
    const A = MatrixSpace(R8, 4, 4).__call__([
      [1, 2, 3, 4],
      [5, 6, 7, 0],
      [1, 3, 5, 7],
      [2, 4, 6, 1],
    ]);
    expect(determinant(A).eq(0)).toBe(true);
  });

  it('should compute a 4x4 determinant over Z/8 with a unit result', () => {
    // det over ZZ of [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,7]] = 7
    const A = MatrixSpace(R8, 4, 4).__call__([
      [1, 0, 0, 0],
      [0, 1, 0, 0],
      [0, 0, 1, 0],
      [0, 0, 0, 7],
    ]);
    expect(determinant(A).eq(7)).toBe(true);
  });

  it('should satisfy adjugate(A)*A == det(A)*I over Z/8 (M57)', () => {
    const A = MatrixSpace(R8, 4, 4).__call__([
      [1, 2, 3, 4],
      [5, 6, 7, 0],
      [1, 3, 5, 7],
      [2, 4, 6, 1],
    ]);
    const d = determinant(A);
    expect(adjugate(A).mul(A).eq(identity_matrix(R8, 4).scalar_mul(d))).toBe(true);
    expect(A.mul(adjugate(A)).eq(identity_matrix(R8, 4).scalar_mul(d))).toBe(true);
  });

  it('should compute 4x4 minors over Z/8', () => {
    const A = MatrixSpace(R8, 4, 4).__call__([
      [1, 2, 3, 4],
      [5, 6, 7, 0],
      [1, 3, 5, 7],
      [2, 4, 6, 1],
    ]);
    expect(minors(A, 4).length).toBe(1);
  });

  it('should agree with the hessenberg algorithm over QQ', () => {
    const A = MatrixSpace(QQ, 4, 4).__call__([
      [2, -1, 0, 3],
      [1, 4, -2, 0],
      [0, 5, 1, -1],
      [3, 0, 2, 6],
    ]);
    expect(determinant(A, 'hessenberg').eq(determinant(A, 'df'))).toBe(true);
    expect(charpoly(A, 'x', 'hessenberg').toString()).toBe(charpoly(A, 'x', 'df').toString());
  });
});

describe('adjugate (M57)', () => {
  it('should match the Sage doctest for [[5,2],[3,4]]', () => {
    // sage: Matrix(ZZ,2,2,[5,2,3,4]).adjugate() == [[4,-2],[-3,5]]
    const A = MatrixSpace(QQ, 2, 2).__call__([
      [5, 2],
      [3, 4],
    ]);
    const N = adjugate(A);
    expect(N.get(0, 0).eq(4)).toBe(true);
    expect(N.get(0, 1).eq(-2)).toBe(true);
    expect(N.get(1, 0).eq(-3)).toBe(true);
    expect(N.get(1, 1).eq(5)).toBe(true);
  });

  it('should return [1] for the 1x1 matrix [2] (Sage doctest)', () => {
    const A = MatrixSpace(QQ, 1, 1).__call__([[2]]);
    expect(adjugate(A).get(0, 0).eq(1)).toBe(true);
  });

  it('should reject non-square matrices', () => {
    expect(() => adjugate(MatrixSpace(QQ, 2, 3).zero())).toThrow(/must be a square matrix/);
  });
});

describe('right_kernel_matrix (M55)', () => {
  const F11 = GF(11n);

  it('should satisfy A*K^T == 0 and have the right dimension (300 random)', () => {
    const rnd = makeRng(4242);
    for (let t = 0; t < 300; t++) {
      const m = 1 + rnd(4);
      const n = 1 + rnd(4);
      const rows: number[][] = [];
      for (let i = 0; i < m; i++) {
        const r: number[] = [];
        for (let j = 0; j < n; j++) r.push(rnd(11));
        rows.push(r);
      }
      const A = MatrixSpace(F11, m, n).__call__(rows);
      const K = right_kernel_matrix(A);
      expect(K.nrows).toBe(n - rank(A));
      if (K.nrows > 0) {
        expect(A.mul(K.transpose()).is_zero()).toBe(true);
      }
    }
  });

  // Regression: right_kernel_matrix needs pivot COLUMNS. It briefly used
  // pivot_rows (= pivots of the transpose), which agrees with the column pivots
  // whenever they form an initial segment [0..r-1] — which is the case for
  // almost every small dense random matrix, so the 300-random sweep above did
  // not catch it. These matrices have a gap in the pivot columns, where the two
  // genuinely differ.
  it('should be correct when the pivot columns are not an initial segment', () => {
    const cases: number[][][] = [
      [
        [1, 2, 0],
        [0, 0, 1],
        [0, 0, 0],
      ], // pivot cols [0,2], pivot rows [0,1]
      [
        [0, 1, 0, 3],
        [0, 0, 0, 1],
        [0, 0, 0, 0],
      ], // pivot cols [1,3]
      [
        [1, 0, 5, 0],
        [0, 0, 0, 1],
      ], // pivot cols [0,3]
    ];
    for (const rows of cases) {
      const A = MatrixSpace(F11, rows.length, rows[0]!.length).__call__(rows);
      const K = right_kernel_matrix(A);
      expect(K.nrows).toBe(rows[0]!.length - rank(A));
      expect(A.mul(K.transpose()).is_zero()).toBe(true);
    }
  });

  it('should return an echelonized basis by default', () => {
    // Sage's default basis over a field is 'echelon'.
    const A = MatrixSpace(F11, 1, 3).__call__([[1, 2, 3]]);
    const K = right_kernel_matrix(A);
    expect(K.nrows).toBe(2);
    // echelon: leading entries are 1 and strictly increasing, and are the only
    // nonzero entry in their column
    expect(K.get(0, 0).eq(1)).toBe(true);
    expect(K.get(1, 0).eq(0)).toBe(true);
    expect(K.get(1, 1).eq(1)).toBe(true);
    expect(K.get(0, 1).eq(0)).toBe(true);
  });

  it("should return the unechelonized pivot basis for basis: 'pivot'", () => {
    const A = MatrixSpace(F11, 1, 3).__call__([[1, 2, 3]]);
    const K = right_kernel_matrix(A, { basis: 'pivot' });
    expect(K.nrows).toBe(2);
    expect(A.mul(K.transpose()).is_zero()).toBe(true);
  });

  it('should reject an unknown basis format', () => {
    const A = MatrixSpace(F11, 1, 2).__call__([[1, 1]]);
    // biome-ignore lint/suspicious/noExplicitAny: testing an invalid option
    expect(() => right_kernel_matrix(A, { basis: 'bogus' as any })).toThrow(/not recognized/);
  });
});

describe('rank, nullity, solve and inverse (M59)', () => {
  const F7 = GF(7n);
  const A = MatrixSpace(F7, 3, 3).__call__([
    [2, 1, 1],
    [1, 3, 2],
    [1, 0, 0],
  ]);

  it('should compute rank and nullities', () => {
    expect(rank(A)).toBe(3);
    expect(left_nullity(A)).toBe(0);
    expect(right_nullity(A)).toBe(0);
    const S = MatrixSpace(F7, 2, 3).__call__([
      [1, 2, 3],
      [2, 4, 6],
    ]);
    expect(rank(S)).toBe(1);
    expect(left_nullity(S)).toBe(1);
    expect(right_nullity(S)).toBe(2);
  });

  it('should invert a matrix', () => {
    expect(A.mul(inverse(A)).eq(MatrixSpace(F7, 3, 3).identity())).toBe(true);
    expect(
      inverse(A)
        .mul(A)
        .eq(MatrixSpace(F7, 3, 3).identity())
    ).toBe(true);
  });

  it('should raise when inverting a singular matrix', () => {
    const S = MatrixSpace(F7, 2, 2).__call__([
      [1, 2],
      [2, 4],
    ]);
    expect(() => inverse(S)).toThrow(/singular/);
  });

  it('should solve A*X = B on the right', () => {
    const B = MatrixSpace(F7, 3, 1).__call__([[1], [2], [3]]);
    const X = solve_right(A, B);
    expect(A.mul(X).eq(B)).toBe(true);
  });

  it('should solve X*A = B on the left', () => {
    const B = MatrixSpace(F7, 1, 3).__call__([[1, 2, 3]]);
    const X = solve_left(A, B);
    expect(X.mul(A).eq(B)).toBe(true);
  });

  it('should raise for an inconsistent system', () => {
    const S = MatrixSpace(F7, 2, 2).__call__([
      [1, 2],
      [2, 4],
    ]);
    const B = MatrixSpace(F7, 2, 1).__call__([[1], [0]]);
    expect(() => solve_right(S, B)).toThrow(/no solutions/);
  });

  it('should give the inverse as pseudoinverse for invertible matrices', () => {
    expect(pseudoinverse(A).eq(inverse(A))).toBe(true);
  });
});

describe('pfaffian and quantum_determinant (M59)', () => {
  const F7 = GF(7n);

  it('should satisfy pf(A)^2 == det(A)', () => {
    const A = MatrixSpace(F7, 4, 4).__call__([
      [0, 1, 2, 3],
      [-1, 0, 4, 5],
      [-2, -4, 0, 6],
      [-3, -5, -6, 0],
    ]);
    const pf = pfaffian(A);
    expect(pf.mul(pf).eq(determinant(A))).toBe(true);
  });

  it('should reject a non skew-symmetric matrix', () => {
    const A = MatrixSpace(F7, 2, 2).__call__([
      [0, 1],
      [1, 0],
    ]);
    expect(() => pfaffian(A)).toThrow(/skew-symmetric/);
  });

  it('should reduce to det at q=1 and to the permanent at q=-1', () => {
    const A = MatrixSpace(F7, 3, 3).__call__([
      [1, 2, 3],
      [4, 5, 6],
      [0, 1, 2],
    ]);
    expect(quantum_determinant(A, F7.__call__(1n)).eq(determinant(A))).toBe(true);
    expect(quantum_determinant(A, F7.__call__(-1n)).eq(permanent(A))).toBe(true);
  });
});

describe('is_diagonalizable and is_semisimple (H46)', () => {
  const F7 = GF(7n);

  it('should return false for a non-trivial Jordan block', () => {
    const J = MatrixSpace(F7, 2, 2).__call__([
      [1, 1],
      [0, 1],
    ]);
    expect(is_diagonalizable(J)).toBe(false);
    expect(is_semisimple(J)).toBe(false);
  });

  it('should return true for a diagonal matrix with distinct entries', () => {
    const D = MatrixSpace(F7, 2, 2).__call__([
      [1, 0],
      [0, 2],
    ]);
    expect(is_diagonalizable(D)).toBe(true);
    expect(is_semisimple(D)).toBe(true);
  });

  it('should match the Sage doctest for [[0,-1],[1,0]] over QQ', () => {
    // sage: A = matrix([[0, -1], [1, 0]])
    // sage: A.is_semisimple()                 -> True
    // sage: A.change_ring(QQ).is_diagonalizable() -> False
    const A = MatrixSpace(QQ, 2, 2).__call__([
      [0, -1],
      [1, 0],
    ]);
    expect(is_semisimple(A)).toBe(true);
    expect(is_diagonalizable(A)).toBe(false);
  });

  it('should raise a TypeError for non-square matrices', () => {
    expect(() => is_diagonalizable(MatrixSpace(F7, 2, 3).zero())).toThrow('not a square matrix');
  });
});

describe('is_similar (H47)', () => {
  const F7 = GF(7n);

  it('should return false for [[1,1],[0,1]] and the identity', () => {
    const J = MatrixSpace(F7, 2, 2).__call__([
      [1, 1],
      [0, 1],
    ]);
    expect(is_similar(J, MatrixSpace(F7, 2, 2).identity())).toBe(false);
  });

  it('should distinguish nilpotent matrices with equal charpoly and minpoly', () => {
    // Both are 6x6 nilpotent with charpoly x^6 and minpoly x^3, but the Jordan
    // block structures are {3,3} and {3,2,1}, so they are not similar.
    const block = (sizes: number[]) => {
      const n = sizes.reduce((a, b) => a + b, 0);
      const rows: number[][] = Array.from({ length: n }, () => new Array<number>(n).fill(0));
      let off = 0;
      for (const s of sizes) {
        for (let i = 0; i < s - 1; i++) rows[off + i]![off + i + 1] = 1;
        off += s;
      }
      return MatrixSpace(F7, n, n).__call__(rows);
    };
    const X = block([3, 3]);
    const Y = block([3, 2, 1]);
    expect(charpoly(X).toString()).toBe(charpoly(Y).toString());
    expect(minpoly(X).toString()).toBe(minpoly(Y).toString());
    expect(is_similar(X, Y)).toBe(false);
    expect(is_similar(X, block([3, 3]))).toBe(true);
  });

  it('should be invariant under conjugation', () => {
    const A = MatrixSpace(F7, 3, 3).__call__([
      [1, 2, 3],
      [0, 4, 5],
      [0, 0, 6],
    ]);
    const P = MatrixSpace(F7, 3, 3).__call__([
      [1, 1, 0],
      [0, 1, 1],
      [1, 0, 1],
    ]);
    const B = inverse(P).mul(A).mul(P);
    expect(is_similar(A, B)).toBe(true);
  });

  it('should raise for non-square or mismatched sizes', () => {
    expect(() => is_similar(MatrixSpace(F7, 2, 3).zero(), MatrixSpace(F7, 2, 3).zero())).toThrow(
      'similarity only makes sense for square matrices'
    );
    expect(() => is_similar(MatrixSpace(F7, 2, 2).zero(), MatrixSpace(F7, 3, 3).zero())).toThrow(
      'matrices do not have the same size'
    );
  });
});

describe('is_hermitian', () => {
  it('should be true exactly for symmetric matrices over a real field', () => {
    expect(
      is_hermitian(
        MatrixSpace(QQ, 2, 2).__call__([
          [1, 2],
          [2, 3],
        ])
      )
    ).toBe(true);
    expect(
      is_hermitian(
        MatrixSpace(QQ, 2, 2).__call__([
          [1, 2],
          [3, 4],
        ])
      )
    ).toBe(false);
    expect(is_hermitian(MatrixSpace(QQ, 2, 3).zero())).toBe(false);
  });
});

describe('eigenvalues extend option (M58)', () => {
  const F7 = GF(7n);

  it('should raise with the default extend=true when the charpoly does not split', () => {
    // x^2 + 1 is irreducible over GF(7)
    const A = MatrixSpace(F7, 2, 2).__call__([
      [0, 6],
      [1, 0],
    ]);
    expect(charpoly(A).coeffs.map((c) => c.value)).toEqual([1n, 0n, 1n]);
    expect(() => eigenvalues(A)).toThrow(/algebraic closure/);
  });

  it('should return only the base-field eigenvalues with extend=false', () => {
    const A = MatrixSpace(F7, 2, 2).__call__([
      [0, 6],
      [1, 0],
    ]);
    expect(eigenvalues(A, false)).toEqual([]);
  });

  it('should return all eigenvalues when the charpoly splits', () => {
    const A = MatrixSpace(F7, 2, 2).__call__([
      [2, 0],
      [0, 5],
    ]);
    expect(
      eigenvalues(A)
        .map((e) => e.value)
        .sort()
    ).toEqual([2n, 5n]);
  });
});

// ============================================================================
// Deferred audit items 16-20 (unit matrix-ops2)
// ============================================================================

describe('minpoly over QQ, non-squarefree branch (deferred 16)', () => {
  const qq = (rows: number[][]) => MatrixSpace(QQ, rows.length, rows[0]!.length).__call__(rows);

  it("should match Sage's minimal_polynomial doctest for matrix(QQ,4,4,range(16))", () => {
    // sage: a = matrix(QQ, 4, 4, range(16))
    // sage: a.minimal_polynomial('z')   ->  z^3 - 30*z^2 - 80*z
    // sage: a.minpoly()                 ->  x^3 - 30*x^2 - 80*x
    // The charpoly x^4 - 30x^3 - 80x^2 is *not* squarefree, so this exercises
    // the radical + kernel-dimension branch of matrix2.pyx:3110-3128, which was
    // unreachable until Polynomial.factor() worked over QQ.
    const A = qq([
      [0, 1, 2, 3],
      [4, 5, 6, 7],
      [8, 9, 10, 11],
      [12, 13, 14, 15],
    ]);
    expect(charpoly(A).toString()).toBe('x^4 + (-30)*x^3 + (-80)*x^2');
    expect(minpoly(A).toString()).toBe('x^3 + (-30)*x^2 + (-80)*x');
    expect(minpoly(A, 'z').toString()).toBe('z^3 + (-30)*z^2 + (-80)*z');
    expect(polyAtMatrix(minpoly(A), A).is_zero()).toBe(true);
    // and it is a *proper* divisor of the charpoly, i.e. the branch really ran
    expect(minpoly(A).degree()).toBe(3);
  });

  it('should give the minpoly of a repeated irreducible quadratic factor', () => {
    // C is the companion matrix of x^2 + 1.
    // diag(C, C)      has charpoly (x^2+1)^2 and minpoly x^2+1
    // [[C, I], [0, C]] has charpoly (x^2+1)^2 and minpoly (x^2+1)^2
    const X = qq([
      [0, -1, 0, 0],
      [1, 0, 0, 0],
      [0, 0, 0, -1],
      [0, 0, 1, 0],
    ]);
    const Y = qq([
      [0, -1, 1, 0],
      [1, 0, 0, 1],
      [0, 0, 0, -1],
      [0, 0, 1, 0],
    ]);
    expect(charpoly(X).toString()).toBe('x^4 + 2*x^2 + 1');
    expect(charpoly(Y).toString()).toBe('x^4 + 2*x^2 + 1');
    expect(minpoly(X).toString()).toBe('x^2 + 1');
    expect(minpoly(Y).toString()).toBe('x^4 + 2*x^2 + 1');
    expect(polyAtMatrix(minpoly(X), X).is_zero()).toBe(true);
    expect(polyAtMatrix(minpoly(Y), Y).is_zero()).toBe(true);
  });

  it('should annihilate the matrix and divide the charpoly on 60 random QQ matrices', () => {
    const rnd = makeRng(20260728);
    for (let t = 0; t < 60; t++) {
      const n = 2 + rnd(3);
      const rows: number[][] = [];
      for (let i = 0; i < n; i++) {
        const r: number[] = [];
        for (let j = 0; j < n; j++) r.push(rnd(9) - 4);
        rows.push(r);
      }
      const A = qq(rows);
      const mp = minpoly(A);
      expect(polyAtMatrix(mp, A).is_zero()).toBe(true);
      expect(mp.degree()).toBeLessThanOrEqual(n);
      // charpoly is divisible by minpoly
      expect(charpoly(A).mod(mp).isZero()).toBe(true);
    }
  });

  it('should decide semisimplicity from a non-squarefree minpoly over QQ', () => {
    // sage: matrix(QQ, [[0,-1],[1,0]]).is_semisimple() -> True
    expect(
      is_semisimple(
        qq([
          [0, -1],
          [1, 0],
        ])
      )
    ).toBe(true);
    // diag(C,C): minpoly x^2+1, squarefree
    expect(
      is_semisimple(
        qq([
          [0, -1, 0, 0],
          [1, 0, 0, 0],
          [0, 0, 0, -1],
          [0, 0, 1, 0],
        ])
      )
    ).toBe(true);
    // [[C,I],[0,C]]: minpoly (x^2+1)^2, not squarefree
    expect(
      is_semisimple(
        qq([
          [0, -1, 1, 0],
          [1, 0, 0, 1],
          [0, 0, 0, -1],
          [0, 0, 1, 0],
        ])
      )
    ).toBe(false);
  });

  it('should separate the two (x^2+1)^2 similarity classes over QQ', () => {
    const X = qq([
      [0, -1, 0, 0],
      [1, 0, 0, 0],
      [0, 0, 0, -1],
      [0, 0, 1, 0],
    ]);
    const Y = qq([
      [0, -1, 1, 0],
      [1, 0, 0, 1],
      [0, 0, 0, -1],
      [0, 0, 1, 0],
    ]);
    expect(is_similar(X, Y)).toBe(false);
    const P = qq([
      [1, 1, 0, 0],
      [0, 1, 1, 0],
      [0, 0, 1, 1],
      [0, 0, 0, 1],
    ]);
    expect(is_similar(Y, inverse(P).mul(Y).mul(P))).toBe(true);
  });
});

describe('norm(A, 2) — largest singular value (deferred 17)', () => {
  const qq = (rows: number[][]) => MatrixSpace(QQ, rows.length, rows[0]!.length).__call__(rows);

  it("should match Sage's sparse-euclidean doctest for [[1,2],[3,4]]", () => {
    // sage: matrix(ZZ, [[1, 2], [3, 4]], sparse=True).norm()
    // 5.464985704219043      (matrix2.pyx:16456)
    expect(
      norm(
        qq([
          [1, 2],
          [3, 4],
        ])
      )
    ).toBe(5.464985704219043);
  });

  it("should match Sage's doctest Id.norm(2) == 1.0", () => {
    // sage: Id = identity_matrix(12); Id.norm(2) -> 1.0
    expect(norm(MatrixSpace(QQ, 12, 12).identity())).toBe(1);
  });

  it("should match Sage's RR doctest matrix(RR,2,2,[13,-4,-4,7]).norm()", () => {
    // sage: A.norm()   # rel tol 2e-16  ->  14.999999999999998
    // The eigenvalues of the symmetric matrix are 15 and 5, so the exact
    // spectral norm is 15; Sage's numerical SVD returns it to 2e-16.
    expect(
      norm(
        qq([
          [13, -4],
          [-4, 7],
        ])
      )
    ).toBe(15);
  });

  it('should be the largest |eigenvalue| for a symmetric matrix', () => {
    // diag(3,3) has a repeated singular value: the charpoly of A^T A is
    // (x-9)^2, so this exercises the squarefree-part reduction.
    expect(
      norm(
        qq([
          [3, 0],
          [0, 3],
        ])
      )
    ).toBe(3);
    expect(
      norm(
        qq([
          [0, 2],
          [2, 0],
        ])
      )
    ).toBe(2);
  });

  it('should handle exact non-integral entries and tiny values', () => {
    const half = QQ.__call__([1n, 2n]);
    const third = QQ.__call__([1n, 3n]);
    const zero = QQ.zero();
    const A = MatrixSpace(QQ, 2, 2).__call__([
      [half, zero],
      [zero, third],
    ]);
    expect(norm(A)).toBe(0.5);
    const tiny = MatrixSpace(QQ, 1, 1).__call__([[QQ.__call__([1n, 10n ** 20n])]]);
    expect(norm(tiny)).toBe(1e-20);
  });

  it('should return 0 for the zero matrix and for empty matrices', () => {
    expect(norm(MatrixSpace(QQ, 3, 3).zero())).toBe(0);
    expect(norm(MatrixSpace(QQ, 0, 0).zero())).toBe(0);
  });

  it('should agree with a numerical Jacobi eigenvalue oracle on 60 random matrices', () => {
    // Independent oracle: the largest eigenvalue of A^T A by the cyclic Jacobi
    // rotation method in double precision.
    const jacobiMaxEig = (M: number[][]): number => {
      const n = M.length;
      const A = M.map((r) => r.slice());
      for (let sweep = 0; sweep < 100; sweep++) {
        let off = 0;
        for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) off += A[i]![j]! * A[i]![j]!;
        if (off < 1e-28) break;
        for (let p = 0; p < n; p++) {
          for (let q = p + 1; q < n; q++) {
            if (Math.abs(A[p]![q]!) < 1e-300) continue;
            const theta = (A[q]![q]! - A[p]![p]!) / (2 * A[p]![q]!);
            const t = (theta >= 0 ? 1 : -1) / (Math.abs(theta) + Math.sqrt(theta * theta + 1));
            const c = 1 / Math.sqrt(t * t + 1);
            const s = t * c;
            for (let k = 0; k < n; k++) {
              const akp = A[k]![p]!;
              const akq = A[k]![q]!;
              A[k]![p] = c * akp - s * akq;
              A[k]![q] = s * akp + c * akq;
            }
            for (let k = 0; k < n; k++) {
              const apk = A[p]![k]!;
              const aqk = A[q]![k]!;
              A[p]![k] = c * apk - s * aqk;
              A[q]![k] = s * apk + c * aqk;
            }
          }
        }
      }
      let mx = Number.NEGATIVE_INFINITY;
      for (let i = 0; i < n; i++) mx = Math.max(mx, A[i]![i]!);
      return mx;
    };

    const rnd = makeRng(20260729);
    for (let t = 0; t < 60; t++) {
      const m = 1 + rnd(4);
      const n = 1 + rnd(4);
      const rows: number[][] = [];
      for (let i = 0; i < m; i++) {
        const r: number[] = [];
        for (let j = 0; j < n; j++) r.push(rnd(21) - 10);
        rows.push(r);
      }
      const got = norm(qq(rows));
      const AtA: number[][] = [];
      for (let i = 0; i < n; i++) {
        const r: number[] = [];
        for (let j = 0; j < n; j++) {
          let s = 0;
          for (let k = 0; k < m; k++) s += rows[k]![i]! * rows[k]![j]!;
          r.push(s);
        }
        AtA.push(r);
      }
      const expected = Math.sqrt(Math.max(0, jacobiMaxEig(AtA)));
      if (expected === 0) {
        expect(got).toBe(0);
      } else {
        expect(Math.abs(got - expected) / expected).toBeLessThan(1e-12);
      }
    }
  });

  it('should reject base rings with no embedding into the complex numbers', () => {
    const A = MatrixSpace(GF(7n), 2, 2).__call__([
      [1, 2],
      [3, 4],
    ]);
    expect(() => norm(A, 2)).toThrow(/no canonical coercion/);
  });
});

describe('is_similar with transformation (deferred 18)', () => {
  const qq = (rows: number[][]) => MatrixSpace(QQ, rows.length, rows[0]!.length).__call__(rows);

  it("should reproduce Sage's is_similar transformation doctest over QQ", () => {
    // sage: A = matrix(ZZ, [[-5, 2, -11], [-6, 7, -42], [0, 1, -6]])
    // sage: B = matrix(ZZ, [[ 1, 12,  3], [-1, -6, -1], [0,  6,  1]])
    // sage: _, T = A.is_similar(B, transformation=True)
    // sage: T.change_ring(QQ)
    // [   1    0    0]
    // [-2/3  1/6 -5/6]
    // [ 2/3    0 -1/3]
    // sage: A == T.inverse()*B*T
    // True
    const A = qq([
      [-5, 2, -11],
      [-6, 7, -42],
      [0, 1, -6],
    ]);
    const B = qq([
      [1, 12, 3],
      [-1, -6, -1],
      [0, 6, 1],
    ]);
    expect(is_similar(A, B)).toBe(true);
    const [ok, T] = is_similar(A, B, true) as [boolean, Matrix<never>];
    expect(ok).toBe(true);
    expect(T).not.toBeNull();
    expect(inverse(T).mul(B).mul(T).eq(A)).toBe(true);
    // entry-for-entry equal to Sage's T.change_ring(QQ)
    const expected = [
      ['1', '0', '0'],
      ['-2/3', '1/6', '-5/6'],
      ['2/3', '0', '-1/3'],
    ];
    for (let i = 0; i < 3; i++) {
      for (let j = 0; j < 3; j++) {
        expect(T.get(i, j).toString()).toBe(expected[i]![j]!);
      }
    }
  });

  it('should return (false, null) when the matrices are not similar', () => {
    const F7 = GF(7n);
    const J = MatrixSpace(F7, 2, 2).__call__([
      [1, 1],
      [0, 1],
    ]);
    expect(is_similar(J, MatrixSpace(F7, 2, 2).identity(), true)).toEqual([false, null]);
  });

  it('should return a verified transformation over GF(7), GF(2) and QQ (150 conjugate pairs)', () => {
    const rnd = makeRng(20260730);
    for (const [ring, q] of [
      [GF(7n), 7],
      [GF(2n), 2],
      [QQ, 9],
    ] as const) {
      for (let t = 0; t < 50; t++) {
        const n = 1 + rnd(4);
        const rows: number[][] = [];
        for (let i = 0; i < n; i++) {
          const r: number[] = [];
          for (let j = 0; j < n; j++) r.push(rnd(q));
          rows.push(r);
        }
        // half the sweep gets a repeated eigenvalue and a non-trivial
        // Jordan structure, which is what makes the transformation hard
        if (t % 2 === 0) {
          for (let i = 0; i < n; i++) for (let j = 0; j < i; j++) rows[i]![j] = 0;
          for (let i = 0; i < n; i++) rows[i]![i] = rows[0]![0]!;
        }
        // biome-ignore lint/suspicious/noExplicitAny: heterogeneous rings in one sweep
        const A = MatrixSpace(ring as any, n, n).__call__(rows);
        let P: typeof A | null = null;
        for (let tries = 0; tries < 60; tries++) {
          const pr: number[][] = [];
          for (let i = 0; i < n; i++) {
            const r: number[] = [];
            for (let j = 0; j < n; j++) r.push(rnd(q));
            pr.push(r);
          }
          // biome-ignore lint/suspicious/noExplicitAny: heterogeneous rings in one sweep
          const cand = MatrixSpace(ring as any, n, n).__call__(pr);
          if (rank(cand) === n) {
            P = cand;
            break;
          }
        }
        if (P === null) continue;
        const B = inverse(P).mul(A).mul(P);
        const [ok, T] = is_similar(A, B, true) as [boolean, typeof A];
        expect(ok).toBe(true);
        expect(inverse(T).mul(B).mul(T).eq(A)).toBe(true);
      }
    }
  });

  it('should agree with brute force over GF(2), 2x2, on every ordered pair', () => {
    // Exhaustive oracle: A ~ B iff some invertible P has P^-1 B P == A.
    const F2 = GF(2n);
    const MS = MatrixSpace(F2, 2, 2);
    const mats: Array<ReturnType<typeof MS.__call__>> = [];
    for (let code = 0; code < 16; code++) {
      mats.push(
        MS.__call__([
          [code & 1, (code >> 1) & 1],
          [(code >> 2) & 1, (code >> 3) & 1],
        ])
      );
    }
    const invertible = mats.filter((m) => rank(m) === 2);
    for (const A of mats) {
      for (const B of mats) {
        let brute = false;
        for (const P of invertible) {
          if (inverse(P).mul(B).mul(P).eq(A)) {
            brute = true;
            break;
          }
        }
        expect(is_similar(A, B)).toBe(brute);
        if (brute) {
          const [ok, T] = is_similar(A, B, true) as [boolean, typeof A];
          expect(ok).toBe(true);
          expect(inverse(T).mul(B).mul(T).eq(A)).toBe(true);
        }
      }
    }
  });
});

describe('change_ring and is_diagonalizable(base_field) (deferred 19)', () => {
  const qq = (rows: number[][]) => MatrixSpace(QQ, rows.length, rows[0]!.length).__call__(rows);

  it('should reduce a rational matrix modulo p', () => {
    // Sage's own doctest (matrix0.pyx:1679) is over GF(25):
    //   sage: matrix(QQ, 2, 2, [1/2, 1/3, 1/3, 1/4]).change_ring(GF(25,'a'))
    //   [3 2]
    //   [2 4]
    // GF(p^k) is not available in this port, so we run the same matrix over
    // GF(7): 1/2 = 4, 1/3 = 5, 1/4 = 2 (mod 7).
    const A = MatrixSpace(QQ, 2, 2).__call__([
      [QQ.__call__([1n, 2n]), QQ.__call__([1n, 3n])],
      [QQ.__call__([1n, 3n]), QQ.__call__([1n, 4n])],
    ]);
    const B = change_ring(A, GF(7n));
    expect(B.get(0, 0).value).toBe(4n);
    expect(B.get(0, 1).value).toBe(5n);
    expect(B.get(1, 0).value).toBe(5n);
    expect(B.get(1, 1).value).toBe(2n);
  });

  it('should map Z/8 onto Z/4 but refuse Z/8 -> GF(7)', () => {
    const A = MatrixSpace(Zmod(8n), 2, 2).__call__([
      [1, 2],
      [3, 4],
    ]);
    const B = change_ring(A, Zmod(4n));
    expect([B.get(0, 0).value, B.get(0, 1).value, B.get(1, 0).value, B.get(1, 1).value]).toEqual([
      1n,
      2n,
      3n,
      0n,
    ]);
    // there is no ring homomorphism Z/8 -> GF(7)
    expect(() => change_ring(A, GF(7n))).toThrow(/unable to coerce/);
  });

  it('should return a copy when the ring is unchanged', () => {
    const A = qq([
      [1, 2],
      [3, 4],
    ]);
    const B = change_ring(A, QQ);
    expect(B).not.toBe(A);
    expect(B.eq(A)).toBe(true);
  });

  it('should test diagonalizability over the requested base field', () => {
    // x^2 + 1 is irreducible over QQ and over GF(7), but splits as
    // (x-2)(x-3) over GF(5).
    const A = qq([
      [0, -1],
      [1, 0],
    ]);
    expect(is_diagonalizable(A)).toBe(false);
    // biome-ignore lint/suspicious/noExplicitAny: base-changing to another ring
    expect(is_diagonalizable(A, GF(5n) as any)).toBe(true);
    // biome-ignore lint/suspicious/noExplicitAny: base-changing to another ring
    expect(is_diagonalizable(A, GF(7n) as any)).toBe(false);
  });

  it('should notice that base-changing can destroy diagonalizability', () => {
    // Distinct eigenvalues 1 and 6 over QQ, but 6 == 1 mod 5, and mod 5 the
    // matrix is a non-trivial Jordan block.
    const A = qq([
      [1, 1],
      [0, 6],
    ]);
    expect(is_diagonalizable(A)).toBe(true);
    // biome-ignore lint/suspicious/noExplicitAny: base-changing to another ring
    expect(is_diagonalizable(A, GF(5n) as any)).toBe(false);
    // biome-ignore lint/suspicious/noExplicitAny: base-changing to another ring
    expect(is_diagonalizable(A, GF(7n) as any)).toBe(true);
  });
});

describe('right_kernel_matrix uses pivot columns (deferred 20)', () => {
  it('should satisfy A*K^T == 0 over QQ including gapped pivot columns', () => {
    const cases: number[][][] = [
      [
        [1, 2, 0],
        [0, 0, 1],
        [0, 0, 0],
      ],
      [
        [0, 1, 0, 3],
        [0, 0, 0, 1],
        [0, 0, 0, 0],
      ],
      [
        [1, 0, 5, 0],
        [0, 0, 0, 1],
      ],
      [
        [2, 4, 6, 8],
        [1, 2, 3, 4],
      ],
    ];
    for (const rows of cases) {
      const A = MatrixSpace(QQ, rows.length, rows[0]!.length).__call__(rows);
      const K = right_kernel_matrix(A);
      expect(K.nrows).toBe(rows[0]!.length - rank(A));
      if (K.nrows > 0) expect(A.mul(K.transpose()).is_zero()).toBe(true);
    }
  });

  it('should satisfy A*K^T == 0 on 150 random QQ matrices', () => {
    const rnd = makeRng(20260731);
    for (let t = 0; t < 150; t++) {
      const m = 1 + rnd(4);
      const n = 1 + rnd(4);
      const rows: number[][] = [];
      for (let i = 0; i < m; i++) {
        const r: number[] = [];
        for (let j = 0; j < n; j++) r.push(rnd(7) - 3);
        rows.push(r);
      }
      const A = MatrixSpace(QQ, m, n).__call__(rows);
      const K = right_kernel_matrix(A);
      expect(K.nrows).toBe(n - rank(A));
      if (K.nrows > 0) expect(A.mul(K.transpose()).is_zero()).toBe(true);
    }
  });
});
