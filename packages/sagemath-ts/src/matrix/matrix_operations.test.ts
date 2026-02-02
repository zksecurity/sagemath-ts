/**
 * Tests for matrix operations (determinant, permanent, minors, etc.)
 */

import { describe, expect, it } from 'vitest';
import { GF } from '../rings/finite_rings/finite_field_constructor.js';
import {
  Matrix,
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
  is_diagonal,
  is_nilpotent,
  is_normal,
  is_one,
  is_positive_definite,
  is_scalar,
  is_similar,
  is_triangular,
  is_unitary,
  kernel,
  left_kernel,
  minors,
  norm,
  permanent,
  permanental_minor,
  right_kernel,
  row_space,
  stack,
  zero_matrix,
} from './index.js';

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

  it('should return 0 for zero matrix', () => {
    const Z = MS2x2.zero();
    expect(density(Z)).toBe(0);
  });

  it('should return 1 for matrix with no zeros', () => {
    const A = MS2x2.__call__([
      [1, 2],
      [3, 4],
    ]);
    expect(density(A)).toBe(1);
  });

  it('should return correct fraction for sparse matrix', () => {
    const A = MS2x2.__call__([
      [1, 0],
      [0, 0],
    ]);
    expect(density(A)).toBe(0.25);
  });

  it('should return correct fraction for identity matrix', () => {
    const I = MS3x3.identity();
    // 3 non-zeros out of 9 entries
    expect(density(I)).toBeCloseTo(1 / 3);
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

  it('should create companion matrix with correct structure', () => {
    // For polynomial x^3 - 2x^2 - 3x - 4, coefficients are [4, 3, 2]
    // (we need coefficients with opposite signs since charpoly is x^n - ...)
    const C = companion_matrix(F7, [F7.__call__(4n), F7.__call__(3n), F7.__call__(2n)]);

    // Right companion form:
    // [0 0 4]
    // [1 0 3]
    // [0 1 2]
    expect(C.get(0, 0).value).toBe(0n);
    expect(C.get(0, 2).value).toBe(4n);
    expect(C.get(1, 0).value).toBe(1n);
    expect(C.get(1, 2).value).toBe(3n);
    expect(C.get(2, 1).value).toBe(1n);
    expect(C.get(2, 2).value).toBe(2n);
  });

  it('should have characteristic polynomial matching input coefficients', () => {
    // The companion matrix for coefficients [c_0, c_1] produces charpoly:
    // det(xI - C) = x^2 - c_1*x - c_0
    //
    // Let's use coefficients [2, 3] in F7
    // Then charpoly should be x^2 - 3x - 2 = x^2 + 4x + 5 (in F7, since -3 = 4, -2 = 5)
    const C = companion_matrix(F7, [F7.__call__(2n), F7.__call__(3n)]);
    const cp = charpoly(C);

    // The charpoly of this companion matrix
    // Companion matrix is:
    // [0 2]
    // [1 3]
    // det(xI - C) = det([x -2; -1 x-3]) = x(x-3) - (-2)(-1) = x^2 - 3x - 2
    // In F7: x^2 + 4x + 5
    expect(cp.coeffs[0]!.value).toBe(5n); // constant term = -2 = 5 mod 7
    expect(cp.coeffs[1]!.value).toBe(4n); // x coefficient = -3 = 4 mod 7
    expect(cp.coeffs[2]!.value).toBe(1n); // x^2 coefficient = 1
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
  const F7 = GF(7n);
  const MS2x2 = MatrixSpace(F7, 2, 2);

  it('should return true for identity matrix', () => {
    const I = MS2x2.identity();
    expect(is_positive_definite(I)).toBe(true);
  });

  it('should return false for zero matrix', () => {
    const Z = MS2x2.zero();
    expect(is_positive_definite(Z)).toBe(false);
  });

  it('should return false for non-symmetric matrix', () => {
    const A = MS2x2.__call__([
      [1, 2],
      [3, 4],
    ]);
    expect(is_positive_definite(A)).toBe(false);
  });

  it('should return false for singular symmetric matrix', () => {
    const A = MS2x2.__call__([
      [1, 1],
      [1, 1],
    ]);
    expect(is_positive_definite(A)).toBe(false);
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
  const F7 = GF(7n);
  const MS2x2 = MatrixSpace(F7, 2, 2);

  it('should return zero for zero matrix', () => {
    const Z = MS2x2.zero();
    expect(norm(Z, 'frob').isZero()).toBe(true);
  });

  it('should compute sum of squares for Frobenius norm', () => {
    const A = MS2x2.__call__([
      [1, 2],
      [3, 4],
    ]);
    // Sum of squares: 1 + 4 + 9 + 16 = 30 = 2 mod 7
    const n = norm(A, 'frob');
    expect(n.value).toBe(2n);
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
