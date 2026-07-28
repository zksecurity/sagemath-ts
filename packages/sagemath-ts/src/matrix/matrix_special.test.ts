/**
 * Tests for matrix_special functions
 */

import { describe, expect, it } from 'vitest';
import { GF } from '../rings/finite_rings/finite_field_constructor.js';
import { PolynomialRing } from '../rings/polynomial/polynomial_ring.js';
import type { Rational } from '../rings/rational.js';
import { QQ } from '../rings/rational_field.js';
import { Matrix } from './matrix_generic.js';
import { determinant, rank as matrix_rank } from './matrix_operations.js';
import {
  apply_map,
  apply_morphism,
  as_bipartite_graph,
  automorphisms_of_rows_and_columns,
  berlekamp_massey,
  block_diagonal_matrix,
  block_matrix,
  circulant,
  column_matrix,
  companion_matrix,
  denominator,
  derivative,
  elementary_matrix,
  elementwise_product,
  find,
  get_bandwidth,
  hadamard_bound,
  hankel,
  hilbert,
  is_permutation_of,
  jordan_block,
  lehmer,
  numerical_approx,
  ones_matrix,
  permutation_normal_form,
  prod_of_row_sums,
  random_diagonalizable_matrix,
  random_echelonizable_matrix,
  random_unimodular_matrix,
  randomize,
  rook_vector,
  set_block,
  subdivide,
  subdivision,
  subdivisions,
  subs,
  tensor_product,
  toeplitz,
  vandermonde,
} from './matrix_special.js';

describe('matrix_special - Subdivisions', () => {
  const F7 = GF(7n);

  describe('subdivide', () => {
    it('should set row and column subdivisions', () => {
      const M = new Matrix(F7, 5, 5, (i, j) => F7.__call__(BigInt(i * 5 + j)));

      subdivide(M, 2, 3);

      const [rows, cols] = subdivisions(M);
      expect(rows).toEqual([2]);
      expect(cols).toEqual([3]);
    });

    it('should handle multiple subdivision lines', () => {
      const M = new Matrix(F7, 5, 5, (i, j) => F7.__call__(BigInt(i * 5 + j)));

      subdivide(M, [1, 3], [2, 4]);

      const [rows, cols] = subdivisions(M);
      expect(rows).toEqual([1, 3]);
      expect(cols).toEqual([2, 4]);
    });

    it('should sort subdivision lines', () => {
      const M = new Matrix(F7, 5, 5, (i, j) => F7.__call__(BigInt(i * 5 + j)));

      subdivide(M, [4, 2], [3, 1]);

      const [rows, cols] = subdivisions(M);
      expect(rows).toEqual([2, 4]);
      expect(cols).toEqual([1, 3]);
    });

    it('should clear subdivisions with empty arrays', () => {
      const M = new Matrix(F7, 3, 3, (i, j) => F7.__call__(BigInt(i * 3 + j)));

      subdivide(M, 1, 1);
      expect(subdivisions(M)).toEqual([[1], [1]]);

      subdivide(M, [], []);
      expect(subdivisions(M)).toEqual([[], []]);
    });

    it('should handle tuple argument', () => {
      const M = new Matrix(F7, 4, 4, (i, j) => F7.__call__(BigInt(i * 4 + j)));

      subdivide(M, [
        [1, 2],
        [2, 3],
      ] as [number[], number[]]);

      const [rows, cols] = subdivisions(M);
      expect(rows).toEqual([1, 2]);
      expect(cols).toEqual([2, 3]);
    });
  });

  describe('subdivision', () => {
    it('should extract subdivision blocks', () => {
      const M = new Matrix(F7, 5, 5, (i, j) => F7.__call__(BigInt(i * 5 + j)));

      subdivide(M, 2, 3);

      // Block (0, 0): rows 0-1, cols 0-2
      const block00 = subdivision(M, 0, 0);
      expect(block00.nrows).toBe(2);
      expect(block00.ncols).toBe(3);
      expect(block00.get(0, 0).value).toBe(0n);
      // row 1, col 2 is index 1*5+2=7, which is 0 mod 7
      expect(block00.get(1, 2).value).toBe(0n);

      // Block (0, 1): rows 0-1, cols 3-4
      const block01 = subdivision(M, 0, 1);
      expect(block01.nrows).toBe(2);
      expect(block01.ncols).toBe(2);
      expect(block01.get(0, 0).value).toBe(3n); // row 0, col 3

      // Block (1, 0): rows 2-4, cols 0-2
      const block10 = subdivision(M, 1, 0);
      expect(block10.nrows).toBe(3);
      expect(block10.ncols).toBe(3);
      // row 2, col 0 is index 2*5+0=10, which is 3 mod 7
      expect(block10.get(0, 0).value).toBe(3n);
    });

    it('should return entire matrix for (0,0) without subdivisions', () => {
      const M = new Matrix(F7, 3, 3, (i, j) => F7.__call__(BigInt(i * 3 + j)));

      const block = subdivision(M, 0, 0);
      expect(block.nrows).toBe(3);
      expect(block.ncols).toBe(3);
      expect(block.eq(M)).toBe(true);
    });

    it('should throw for out-of-bounds indices', () => {
      const M = new Matrix(F7, 5, 5, (i, j) => F7.__call__(BigInt(i * 5 + j)));
      subdivide(M, 2, 3);

      expect(() => subdivision(M, 2, 0)).toThrow();
      expect(() => subdivision(M, 0, 2)).toThrow();
      expect(() => subdivision(M, -1, 0)).toThrow();
    });
  });

  describe('subdivisions', () => {
    it('should return empty arrays for unsubdivided matrix', () => {
      const M = new Matrix(F7, 3, 3, (i, j) => F7.__call__(BigInt(i * 3 + j)));
      expect(subdivisions(M)).toEqual([[], []]);
    });

    it('should return current subdivisions', () => {
      const M = new Matrix(F7, 5, 5, (i, j) => F7.__call__(BigInt(i * 5 + j)));
      subdivide(M, [2], [3]);

      const [rows, cols] = subdivisions(M);
      expect(rows).toEqual([2]);
      expect(cols).toEqual([3]);
    });
  });
});

describe('matrix_special - Bipartite Graph', () => {
  const F7 = GF(7n);

  describe('as_bipartite_graph', () => {
    it('should create bipartite graph representation', () => {
      const M = new Matrix(F7, 2, 3, [
        [F7.__call__(1n), F7.__call__(2n), F7.__call__(3n)],
        [F7.__call__(4n), F7.__call__(5n), F7.__call__(6n)],
      ]);

      const graph = as_bipartite_graph(M);

      // Left vertices: 1, 2 (rows)
      expect(graph.left).toEqual([1, 2]);

      // Right vertices: 3, 4, 5 (columns, offset by nrows)
      expect(graph.right).toEqual([3, 4, 5]);

      // Should have 6 edges (2 rows * 3 cols)
      expect(graph.edges.length).toBe(6);

      // Check specific edges
      // Edge from row 0 to col 0: (1, 3, M[0,0])
      const edge00 = graph.edges.find((e) => e[0] === 1 && e[1] === 3);
      expect(edge00).toBeDefined();
      expect(edge00![2].value).toBe(1n);
    });

    it('should handle empty matrix', () => {
      const M = new Matrix(F7, 0, 0);

      const graph = as_bipartite_graph(M);

      expect(graph.left).toEqual([]);
      expect(graph.right).toEqual([]);
      expect(graph.edges).toEqual([]);
    });
  });
});

describe('matrix_special - Automorphisms', () => {
  const F7 = GF(7n);

  describe('automorphisms_of_rows_and_columns', () => {
    it('should find identity automorphism', () => {
      const M = new Matrix(F7, 2, 2, [
        [F7.__call__(1n), F7.__call__(2n)],
        [F7.__call__(3n), F7.__call__(4n)],
      ]);

      const auts = automorphisms_of_rows_and_columns(M);

      // At minimum, identity should be present
      expect(auts.length).toBeGreaterThanOrEqual(1);

      // Check identity is present
      const hasIdentity = auts.some(
        ([rowPerm, colPerm]) => rowPerm.every((v, i) => v === i) && colPerm.every((v, i) => v === i)
      );
      expect(hasIdentity).toBe(true);
    });

    it('should find row swap automorphism for symmetric rows', () => {
      // Matrix with identical rows
      const M = new Matrix(F7, 2, 2, [
        [F7.__call__(1n), F7.__call__(2n)],
        [F7.__call__(1n), F7.__call__(2n)],
      ]);

      const auts = automorphisms_of_rows_and_columns(M);

      // Should include row swap [1, 0], [0, 1]
      expect(auts.length).toBeGreaterThan(1);
    });

    it('should work with empty matrix', () => {
      const M = new Matrix(F7, 0, 0);

      const auts = automorphisms_of_rows_and_columns(M);

      expect(auts).toEqual([[[], []]]);
    });
  });
});

describe('matrix_special - Permutation Normal Form', () => {
  const F7 = GF(7n);

  describe('permutation_normal_form', () => {
    it('should return lexicographically maximal permutation', () => {
      const M = new Matrix(F7, 2, 2, [
        [F7.__call__(1n), F7.__call__(2n)],
        [F7.__call__(3n), F7.__call__(4n)],
      ]);

      const normal = permutation_normal_form(M) as Matrix<
        typeof M extends Matrix<infer R> ? R : never
      >;

      // The maximal permutation should have largest elements first
      // 4, 3 in first row, 2, 1 in second (or some permutation giving max lex order)
      expect(normal.get(0, 0).value).toBeGreaterThanOrEqual(normal.get(1, 0).value);
    });

    it('should return permutation with check=true', () => {
      const M = new Matrix(F7, 2, 2, [
        [F7.__call__(1n), F7.__call__(2n)],
        [F7.__call__(3n), F7.__call__(4n)],
      ]);

      const result = permutation_normal_form(M, true) as [
        Matrix<typeof M extends Matrix<infer R> ? R : never>,
        [number[], number[]],
      ];

      expect(Array.isArray(result)).toBe(true);
      expect(result[0]).toBeInstanceOf(Matrix);
      expect(Array.isArray(result[1])).toBe(true);
    });

    it('should handle 3x3 matrix', () => {
      const M = new Matrix(F7, 3, 3, [
        [F7.__call__(0n), F7.__call__(0n), F7.__call__(1n)],
        [F7.__call__(1n), F7.__call__(0n), F7.__call__(2n)],
        [F7.__call__(0n), F7.__call__(0n), F7.__call__(0n)],
      ]);

      const normal = permutation_normal_form(M) as Matrix<
        typeof M extends Matrix<infer R> ? R : never
      >;

      // Result should be well-defined
      expect(normal.nrows).toBe(3);
      expect(normal.ncols).toBe(3);
    });
  });

  describe('is_permutation_of', () => {
    it('should return true for same matrix', () => {
      const M = new Matrix(F7, 2, 2, [
        [F7.__call__(1n), F7.__call__(2n)],
        [F7.__call__(3n), F7.__call__(4n)],
      ]);

      expect(is_permutation_of(M, M)).toBe(true);
    });

    it('should return true for row permutation', () => {
      const M = new Matrix(F7, 2, 2, [
        [F7.__call__(1n), F7.__call__(2n)],
        [F7.__call__(3n), F7.__call__(4n)],
      ]);

      const N = new Matrix(F7, 2, 2, [
        [F7.__call__(3n), F7.__call__(4n)],
        [F7.__call__(1n), F7.__call__(2n)],
      ]);

      expect(is_permutation_of(M, N)).toBe(true);
    });

    it('should return true for column permutation', () => {
      const M = new Matrix(F7, 2, 2, [
        [F7.__call__(1n), F7.__call__(2n)],
        [F7.__call__(3n), F7.__call__(4n)],
      ]);

      const N = new Matrix(F7, 2, 2, [
        [F7.__call__(2n), F7.__call__(1n)],
        [F7.__call__(4n), F7.__call__(3n)],
      ]);

      expect(is_permutation_of(M, N)).toBe(true);
    });

    it('should return false for different matrices', () => {
      const M = new Matrix(F7, 2, 2, [
        [F7.__call__(1n), F7.__call__(2n)],
        [F7.__call__(3n), F7.__call__(4n)],
      ]);

      const N = new Matrix(F7, 2, 2, [
        [F7.__call__(1n), F7.__call__(2n)],
        [F7.__call__(3n), F7.__call__(5n)], // Different value
      ]);

      expect(is_permutation_of(M, N)).toBe(false);
    });

    it('should return false for different dimensions', () => {
      const M = new Matrix(F7, 2, 2, [
        [F7.__call__(1n), F7.__call__(2n)],
        [F7.__call__(3n), F7.__call__(4n)],
      ]);

      const N = new Matrix(F7, 2, 3, [
        [F7.__call__(1n), F7.__call__(2n), F7.__call__(3n)],
        [F7.__call__(4n), F7.__call__(5n), F7.__call__(6n)],
      ]);

      expect(is_permutation_of(M, N)).toBe(false);
    });

    it('should return permutation with check=true', () => {
      const M = new Matrix(F7, 2, 2, [
        [F7.__call__(1n), F7.__call__(2n)],
        [F7.__call__(3n), F7.__call__(4n)],
      ]);

      const N = new Matrix(F7, 2, 2, [
        [F7.__call__(3n), F7.__call__(4n)],
        [F7.__call__(1n), F7.__call__(2n)],
      ]);

      const result = is_permutation_of(M, N, true) as [boolean, [number[], number[]] | null];

      expect(result[0]).toBe(true);
    });
  });
});

describe('matrix_special - Block Matrices', () => {
  const F7 = GF(7n);

  describe('block_matrix', () => {
    it('should construct 2x2 block matrix', () => {
      const A = new Matrix(F7, 2, 2, [
        [F7.__call__(1n), F7.__call__(2n)],
        [F7.__call__(3n), F7.__call__(4n)],
      ]);

      const B = new Matrix(F7, 2, 2, [
        [F7.__call__(5n), F7.__call__(6n)],
        [F7.__call__(0n), F7.__call__(1n)],
      ]);

      const result = block_matrix(F7, [
        [A, B],
        [B, A],
      ]);

      expect(result.nrows).toBe(4);
      expect(result.ncols).toBe(4);

      // Check top-left block
      expect(result.get(0, 0).value).toBe(1n);
      expect(result.get(0, 1).value).toBe(2n);

      // Check top-right block
      expect(result.get(0, 2).value).toBe(5n);
      expect(result.get(0, 3).value).toBe(6n);

      // Check bottom-left block
      expect(result.get(2, 0).value).toBe(5n);
      expect(result.get(2, 1).value).toBe(6n);

      // Check bottom-right block
      expect(result.get(2, 2).value).toBe(1n);
      expect(result.get(2, 3).value).toBe(2n);
    });

    it('should handle zero blocks', () => {
      const A = new Matrix(F7, 2, 2, [
        [F7.__call__(1n), F7.__call__(2n)],
        [F7.__call__(3n), F7.__call__(4n)],
      ]);

      const result = block_matrix(F7, [
        [A, 0],
        [0, A],
      ]);

      expect(result.nrows).toBe(4);
      expect(result.ncols).toBe(4);

      // Zero blocks should be zero
      expect(result.get(0, 2).value).toBe(0n);
      expect(result.get(2, 0).value).toBe(0n);
    });
  });

  describe('block_diagonal_matrix', () => {
    it('should construct block diagonal matrix', () => {
      const A = new Matrix(F7, 2, 2, [
        [F7.__call__(1n), F7.__call__(2n)],
        [F7.__call__(3n), F7.__call__(4n)],
      ]);

      const B = new Matrix(F7, 2, 2, [
        [F7.__call__(5n), F7.__call__(6n)],
        [F7.__call__(0n), F7.__call__(1n)],
      ]);

      const result = block_diagonal_matrix(F7, A, B);

      expect(result.nrows).toBe(4);
      expect(result.ncols).toBe(4);

      // A block
      expect(result.get(0, 0).value).toBe(1n);
      expect(result.get(1, 1).value).toBe(4n);

      // B block
      expect(result.get(2, 2).value).toBe(5n);
      expect(result.get(3, 3).value).toBe(1n);

      // Off-diagonal blocks should be zero
      expect(result.get(0, 2).value).toBe(0n);
      expect(result.get(2, 0).value).toBe(0n);
    });

    it('should handle single block', () => {
      const A = new Matrix(F7, 2, 2, [
        [F7.__call__(1n), F7.__call__(2n)],
        [F7.__call__(3n), F7.__call__(4n)],
      ]);

      const result = block_diagonal_matrix(F7, A);

      expect(result.eq(A)).toBe(true);
    });

    it('should handle empty input', () => {
      const result = block_diagonal_matrix(F7);

      expect(result.nrows).toBe(0);
      expect(result.ncols).toBe(0);
    });
  });
});

describe('matrix_special - Jordan Block', () => {
  const F7 = GF(7n);

  describe('jordan_block', () => {
    it('should create Jordan block of size 3', () => {
      const lambda = F7.__call__(2n);
      const J = jordan_block(F7, lambda, 3);

      expect(J.nrows).toBe(3);
      expect(J.ncols).toBe(3);

      // Diagonal entries
      expect(J.get(0, 0).value).toBe(2n);
      expect(J.get(1, 1).value).toBe(2n);
      expect(J.get(2, 2).value).toBe(2n);

      // Superdiagonal entries
      expect(J.get(0, 1).value).toBe(1n);
      expect(J.get(1, 2).value).toBe(1n);

      // Other entries
      expect(J.get(1, 0).value).toBe(0n);
      expect(J.get(2, 1).value).toBe(0n);
      expect(J.get(0, 2).value).toBe(0n);
    });

    it('should create 1x1 Jordan block', () => {
      const lambda = F7.__call__(5n);
      const J = jordan_block(F7, lambda, 1);

      expect(J.nrows).toBe(1);
      expect(J.ncols).toBe(1);
      expect(J.get(0, 0).value).toBe(5n);
    });

    it('should create empty Jordan block', () => {
      const lambda = F7.__call__(3n);
      const J = jordan_block(F7, lambda, 0);

      expect(J.nrows).toBe(0);
      expect(J.ncols).toBe(0);
    });
  });
});

describe('matrix_special - Companion Matrix', () => {
  const F7 = GF(7n);

  describe('companion_matrix', () => {
    // SageMath takes ALL coefficients of a monic polynomial, low degree first,
    // and puts their NEGATIVES on the indicated border.  The previous version of
    // these tests pinned the port's own convention (leading coefficient omitted,
    // coefficients not negated, 'left' laid out as Sage's 'top'); the values
    // below are the ones `sage: companion_matrix(...)` prints.
    // Reference: sage/matrix/special.py:2301 (doctest with poly = [-2,-3,-4,-5,-6,1]).
    const poly = [-2, -3, -4, -5, -6, 1].map((c) => F7.__call__(BigInt(((c % 7) + 7) % 7)));

    const asValues = (M: Matrix<ReturnType<typeof F7.__call__>>) =>
      M.rows().map((row) => row.map((x) => Number(x.value)));

    it('should create the right companion matrix', () => {
      // Sage over ZZ:            over GF(7) the negatives -(-2) = 2 etc. are unchanged
      // [0 0 0 0 2]
      // [1 0 0 0 3]
      // [0 1 0 0 4]
      // [0 0 1 0 5]
      // [0 0 0 1 6]
      expect(asValues(companion_matrix(F7, poly, 'right'))).toEqual([
        [0, 0, 0, 0, 2],
        [1, 0, 0, 0, 3],
        [0, 1, 0, 0, 4],
        [0, 0, 1, 0, 5],
        [0, 0, 0, 1, 6],
      ]);
    });

    it('should create the left companion matrix', () => {
      expect(asValues(companion_matrix(F7, poly, 'left'))).toEqual([
        [6, 1, 0, 0, 0],
        [5, 0, 1, 0, 0],
        [4, 0, 0, 1, 0],
        [3, 0, 0, 0, 1],
        [2, 0, 0, 0, 0],
      ]);
    });

    it('should create the bottom companion matrix', () => {
      expect(asValues(companion_matrix(F7, poly, 'bottom'))).toEqual([
        [0, 1, 0, 0, 0],
        [0, 0, 1, 0, 0],
        [0, 0, 0, 1, 0],
        [0, 0, 0, 0, 1],
        [2, 3, 4, 5, 6],
      ]);
    });

    it('should create the top companion matrix', () => {
      expect(asValues(companion_matrix(F7, poly, 'top'))).toEqual([
        [6, 5, 4, 3, 2],
        [1, 0, 0, 0, 0],
        [0, 1, 0, 0, 0],
        [0, 0, 1, 0, 0],
        [0, 0, 0, 1, 0],
      ]);
    });

    it('should build a 0x0 matrix from the constant polynomial 1', () => {
      const C = companion_matrix(F7, [F7.__call__(1n)]);

      expect(C.nrows).toBe(0);
      expect(C.ncols).toBe(0);
    });

    it('should reject an empty coefficient list', () => {
      expect(() => companion_matrix(F7, [])).toThrow(
        'polynomial cannot be specified by an empty list'
      );
    });

    it('should reject a non-monic polynomial', () => {
      expect(() =>
        companion_matrix(
          F7,
          [2, 3, 5].map((c) => F7.__call__(BigInt(c)))
        )
      ).toThrow('must be monic');
    });

    it('should reject an unknown format', () => {
      expect(() => companion_matrix(F7, poly, 'junk' as unknown as 'right')).toThrow(
        "format must be 'right', 'left', 'top' or 'bottom', not junk"
      );
    });
  });
});

describe('matrix_special - Structured Matrices', () => {
  const F7 = GF(7n);

  describe('circulant', () => {
    it('should create circulant matrix', () => {
      const v = [F7.__call__(1n), F7.__call__(2n), F7.__call__(3n)];
      const C = circulant(F7, v);

      expect(C.nrows).toBe(3);
      expect(C.ncols).toBe(3);

      // First row is v
      expect(C.get(0, 0).value).toBe(1n);
      expect(C.get(0, 1).value).toBe(2n);
      expect(C.get(0, 2).value).toBe(3n);

      // Each row is cyclic shift of previous
      expect(C.get(1, 0).value).toBe(3n);
      expect(C.get(1, 1).value).toBe(1n);
      expect(C.get(1, 2).value).toBe(2n);

      expect(C.get(2, 0).value).toBe(2n);
      expect(C.get(2, 1).value).toBe(3n);
      expect(C.get(2, 2).value).toBe(1n);
    });
  });

  describe('toeplitz', () => {
    // In SageMath `r` is the first row *counting from the second column*, so the
    // matrix has len(c) rows and len(r)+1 columns.  The earlier version of this
    // test treated r[0] as the diagonal entry, which is the port's old (wrong)
    // convention.  Reference: sage/matrix/special.py:3801.
    const asValues = (M: Matrix<ReturnType<typeof F7.__call__>>) =>
      M.rows().map((row) => row.map((x) => Number(x.value)));

    it('should reproduce the SageMath doctest matrix.toeplitz([1..4], [5..6])', () => {
      const c = [1, 2, 3, 4].map((x) => F7.__call__(BigInt(x)));
      const r = [5, 6].map((x) => F7.__call__(BigInt(x)));

      const T = toeplitz(F7, c, r);

      expect(T.nrows).toBe(4);
      expect(T.ncols).toBe(3);
      expect(asValues(T)).toEqual([
        [1, 5, 6],
        [2, 1, 5],
        [3, 2, 1],
        [4, 3, 2],
      ]);
    });

    it('should reproduce the boundary value problem doctest', () => {
      // matrix.toeplitz([-2, 1, 0, 0], [1, 0, 0]) over GF(7): -2 = 5
      const c = [5, 1, 0, 0].map((x) => F7.__call__(BigInt(x)));
      const r = [1, 0, 0].map((x) => F7.__call__(BigInt(x)));

      expect(asValues(toeplitz(F7, c, r))).toEqual([
        [5, 1, 0, 0],
        [1, 5, 1, 0],
        [0, 1, 5, 1],
        [0, 0, 1, 5],
      ]);
    });
  });

  describe('hankel', () => {
    // In SageMath `r` is the last row *from the second to the last column*, so
    // the matrix has len(c) rows and len(r)+1 columns; omitting r gives a square
    // matrix filled with zeros below the anti-diagonal.  The earlier version of
    // this test pinned the port's old convention, in which r[0] was ignored.
    // Reference: sage/matrix/special.py:3847.
    const asValues = (M: Matrix<ReturnType<typeof F7.__call__>>) =>
      M.rows().map((row) => row.map((x) => Number(x.value)));

    it('should reproduce the SageMath doctest matrix.hankel([1, 2, 3])', () => {
      const c = [1, 2, 3].map((x) => F7.__call__(BigInt(x)));

      const H = hankel(F7, c);

      expect(H.nrows).toBe(3);
      expect(H.ncols).toBe(3);
      expect(asValues(H)).toEqual([
        [1, 2, 3],
        [2, 3, 0],
        [3, 0, 0],
      ]);
    });

    it('should reproduce the SageMath doctest matrix.hankel([1..3], [7..10])', () => {
      const c = [1, 2, 3].map((x) => F7.__call__(BigInt(x)));
      const r = [7, 8, 9, 10].map((x) => F7.__call__(BigInt(x)));

      const H = hankel(F7, c, r);

      expect(H.nrows).toBe(3);
      expect(H.ncols).toBe(5);
      // over GF(7): 7 = 0, 8 = 1, 9 = 2, 10 = 3
      expect(asValues(H)).toEqual([
        [1, 2, 3, 0, 1],
        [2, 3, 0, 1, 2],
        [3, 0, 1, 2, 3],
      ]);
    });
  });

  describe('vandermonde', () => {
    it('should create Vandermonde matrix', () => {
      const v = [F7.__call__(1n), F7.__call__(2n), F7.__call__(3n)];
      const V = vandermonde(F7, v);

      expect(V.nrows).toBe(3);
      expect(V.ncols).toBe(3);

      // Row i: [1, v[i], v[i]^2, ..., v[i]^(n-1)]
      // Row 0: [1, 1, 1]
      expect(V.get(0, 0).value).toBe(1n);
      expect(V.get(0, 1).value).toBe(1n);
      expect(V.get(0, 2).value).toBe(1n);

      // Row 1: [1, 2, 4]
      expect(V.get(1, 0).value).toBe(1n);
      expect(V.get(1, 1).value).toBe(2n);
      expect(V.get(1, 2).value).toBe(4n);

      // Row 2: [1, 3, 9] = [1, 3, 2] mod 7
      expect(V.get(2, 0).value).toBe(1n);
      expect(V.get(2, 1).value).toBe(3n);
      expect(V.get(2, 2).value).toBe(2n); // 9 mod 7 = 2
    });
  });
});

describe('matrix_special - Utility Functions', () => {
  const F7 = GF(7n);

  describe('ones_matrix', () => {
    it('should create matrix of all ones', () => {
      const M = ones_matrix(F7, 2, 3);

      expect(M.nrows).toBe(2);
      expect(M.ncols).toBe(3);

      for (let i = 0; i < 2; i++) {
        for (let j = 0; j < 3; j++) {
          expect(M.get(i, j).value).toBe(1n);
        }
      }
    });

    it('should default to square matrix', () => {
      const M = ones_matrix(F7, 3);

      expect(M.nrows).toBe(3);
      expect(M.ncols).toBe(3);
    });
  });

  describe('column_matrix', () => {
    it('should create matrix from column vectors', () => {
      const col1 = [F7.__call__(1n), F7.__call__(2n)];
      const col2 = [F7.__call__(3n), F7.__call__(4n)];

      const M = column_matrix(F7, [col1, col2]);

      expect(M.nrows).toBe(2);
      expect(M.ncols).toBe(2);

      expect(M.get(0, 0).value).toBe(1n);
      expect(M.get(1, 0).value).toBe(2n);
      expect(M.get(0, 1).value).toBe(3n);
      expect(M.get(1, 1).value).toBe(4n);
    });
  });

  describe('elementary_matrix', () => {
    it('should create row swap elementary matrix', () => {
      const E = elementary_matrix(F7, 3, { row1: 0, row2: 2 });

      expect(E.nrows).toBe(3);
      expect(E.ncols).toBe(3);

      // Should swap rows 0 and 2
      expect(E.get(0, 0).value).toBe(0n);
      expect(E.get(0, 2).value).toBe(1n);
      expect(E.get(2, 0).value).toBe(1n);
      expect(E.get(2, 2).value).toBe(0n);
      expect(E.get(1, 1).value).toBe(1n);
    });

    it('should create row addition elementary matrix', () => {
      const E = elementary_matrix(F7, 3, { row1: 0, row2: 1, scale: F7.__call__(2n) });

      expect(E.nrows).toBe(3);
      expect(E.ncols).toBe(3);

      // Should add 2 * row1 to row0
      expect(E.get(0, 0).value).toBe(1n);
      expect(E.get(0, 1).value).toBe(2n);
      expect(E.get(1, 1).value).toBe(1n);
    });

    it('should create row scaling elementary matrix', () => {
      const E = elementary_matrix(F7, 3, { row1: 1, scale: F7.__call__(3n) });

      expect(E.nrows).toBe(3);
      expect(E.ncols).toBe(3);

      // Should scale row1 by 3
      expect(E.get(0, 0).value).toBe(1n);
      expect(E.get(1, 1).value).toBe(3n);
      expect(E.get(2, 2).value).toBe(1n);
    });

    it('should create column operation matrices as the transpose', () => {
      const row = elementary_matrix(F7, 4, { row1: 3, row2: 2, scale: F7.__call__(3n) });
      const col = elementary_matrix(F7, 4, { col1: 3, col2: 2, scale: F7.__call__(3n) });

      expect(col.eq(row.transpose())).toBe(true);
      expect(col.get(2, 3).value).toBe(3n);
    });

    it('should reject adding a multiple of a row/column to itself', () => {
      // Sage: ValueError: cannot add a multiple of a row to itself
      expect(() => elementary_matrix(F7, 5, { row1: 3, row2: 3, scale: F7.__call__(2n) })).toThrow(
        'cannot add a multiple of a row to itself'
      );
      expect(() => elementary_matrix(F7, 5, { col1: 3, col2: 3, scale: F7.__call__(2n) })).toThrow(
        'cannot add a multiple of a column to itself'
      );
    });

    it('should reject a zero scale', () => {
      // Sage: ValueError: scale parameter of row of elementary matrix must be nonzero
      expect(() => elementary_matrix(F7, 5, { row1: 3, scale: F7.__call__(0n) })).toThrow(
        'scale parameter of row of elementary matrix must be nonzero'
      );
      expect(() => elementary_matrix(F7, 5, { col1: 3, scale: F7.__call__(0n) })).toThrow(
        'scale parameter of column of elementary matrix must be nonzero'
      );
    });

    it('should validate the remaining arguments as SageMath does', () => {
      expect(() => elementary_matrix(F7, 0, { row1: 0, scale: F7.__call__(1n) })).toThrow(
        'size of elementary matrix must be 1 or greater, not 0'
      );
      expect(() => elementary_matrix(F7, 5, { scale: F7.__call__(1n) })).toThrow(
        'row1 or col1 must be specified'
      );
      expect(() => elementary_matrix(F7, 5, { row1: 1, col1: 2 })).toThrow(
        'cannot specify both row1 and col1'
      );
      expect(() => elementary_matrix(F7, 5, { row1: 7, scale: F7.__call__(1n) })).toThrow(
        'row of elementary matrix must be positive and smaller than 5, not 7'
      );
      expect(() => elementary_matrix(F7, 5, { row1: 2 })).toThrow(
        'insufficient parameters provided to construct elementary matrix'
      );
    });
  });

  describe('tensor_product', () => {
    it('should compute tensor product', () => {
      const A = new Matrix(F7, 2, 2, [
        [F7.__call__(1n), F7.__call__(2n)],
        [F7.__call__(3n), F7.__call__(4n)],
      ]);

      const B = new Matrix(F7, 2, 2, [
        [F7.__call__(5n), F7.__call__(6n)],
        [F7.__call__(0n), F7.__call__(1n)],
      ]);

      const T = tensor_product(A, B);

      expect(T.nrows).toBe(4);
      expect(T.ncols).toBe(4);

      // Top-left 2x2 block: 1 * B
      expect(T.get(0, 0).value).toBe(5n);
      expect(T.get(0, 1).value).toBe(6n);
      expect(T.get(1, 0).value).toBe(0n);
      expect(T.get(1, 1).value).toBe(1n);

      // Top-right 2x2 block: 2 * B
      expect(T.get(0, 2).value).toBe(3n); // 2*5 = 10 mod 7 = 3
      expect(T.get(0, 3).value).toBe(5n); // 2*6 = 12 mod 7 = 5
    });
  });

  describe('elementwise_product', () => {
    it('should compute Hadamard product', () => {
      const A = new Matrix(F7, 2, 2, [
        [F7.__call__(1n), F7.__call__(2n)],
        [F7.__call__(3n), F7.__call__(4n)],
      ]);

      const B = new Matrix(F7, 2, 2, [
        [F7.__call__(5n), F7.__call__(6n)],
        [F7.__call__(0n), F7.__call__(2n)],
      ]);

      const H = elementwise_product(A, B);

      expect(H.nrows).toBe(2);
      expect(H.ncols).toBe(2);

      expect(H.get(0, 0).value).toBe(5n); // 1*5
      expect(H.get(0, 1).value).toBe(5n); // 2*6 = 12 mod 7 = 5
      expect(H.get(1, 0).value).toBe(0n); // 3*0
      expect(H.get(1, 1).value).toBe(1n); // 4*2 = 8 mod 7 = 1
    });

    it('should throw for different dimensions', () => {
      const A = new Matrix(F7, 2, 2, [
        [F7.__call__(1n), F7.__call__(2n)],
        [F7.__call__(3n), F7.__call__(4n)],
      ]);

      const B = new Matrix(F7, 2, 3, [
        [F7.__call__(5n), F7.__call__(6n), F7.__call__(0n)],
        [F7.__call__(0n), F7.__call__(2n), F7.__call__(1n)],
      ]);

      expect(() => elementwise_product(A, B)).toThrow();
    });
  });

  describe('set_block', () => {
    it('should set a block in the matrix', () => {
      const M = new Matrix(F7, 4, 4, (i, j) => F7.__call__(0n));

      const block = new Matrix(F7, 2, 2, [
        [F7.__call__(1n), F7.__call__(2n)],
        [F7.__call__(3n), F7.__call__(4n)],
      ]);

      set_block(M, 1, 1, block);

      expect(M.get(1, 1).value).toBe(1n);
      expect(M.get(1, 2).value).toBe(2n);
      expect(M.get(2, 1).value).toBe(3n);
      expect(M.get(2, 2).value).toBe(4n);

      // Other entries should be unchanged
      expect(M.get(0, 0).value).toBe(0n);
      expect(M.get(3, 3).value).toBe(0n);
    });

    it('should throw for out-of-bounds block', () => {
      const M = new Matrix(F7, 3, 3, (i, j) => F7.__call__(0n));

      const block = new Matrix(F7, 2, 2, [
        [F7.__call__(1n), F7.__call__(2n)],
        [F7.__call__(3n), F7.__call__(4n)],
      ]);

      expect(() => set_block(M, 2, 2, block)).toThrow();
    });
  });

  describe('find', () => {
    it('should find entries matching condition', () => {
      const M = new Matrix(F7, 3, 3, (i, j) => F7.__call__(BigInt(i * 3 + j)));

      // Find entries > 4 (values mod 7: 0,1,2,3,4,5,6,0,1 -> only 5 and 6 are > 4)
      const matching = find(M, (x) => x.value > 4n) as Array<
        typeof M extends Matrix<infer R> ? R : never
      >;

      expect(matching.length).toBe(2); // 5 at (1,2) and 6 at (2,0)
    });

    it('should return indices when requested', () => {
      const M = new Matrix(F7, 3, 3, (i, j) => F7.__call__(BigInt(i * 3 + j)));

      // Find positions of entries > 4 (5 at (1,2) and 6 at (2,0))
      const positions = find(M, (x) => x.value > 4n, true) as Array<[number, number]>;

      expect(positions.length).toBe(2);
      expect(positions).toContainEqual([1, 2]); // value 5
      expect(positions).toContainEqual([2, 0]); // value 6
    });
  });

  describe('apply_map', () => {
    it('should apply function to each entry', () => {
      const M = new Matrix(F7, 2, 2, [
        [F7.__call__(1n), F7.__call__(2n)],
        [F7.__call__(3n), F7.__call__(4n)],
      ]);

      const doubled = apply_map(M, (x) => x.mul(F7.__call__(2n)), F7);

      expect(doubled.get(0, 0).value).toBe(2n);
      expect(doubled.get(0, 1).value).toBe(4n);
      expect(doubled.get(1, 0).value).toBe(6n);
      expect(doubled.get(1, 1).value).toBe(1n); // 8 mod 7 = 1
    });
  });

  describe('get_bandwidth', () => {
    it('should compute bandwidth of diagonal matrix', () => {
      const M = new Matrix(F7, 3, 3, (i, j) => (i === j ? F7.__call__(1n) : F7.__call__(0n)));

      const [lower, upper] = get_bandwidth(M);

      expect(lower).toBe(0);
      expect(upper).toBe(0);
    });

    it('should compute bandwidth of tridiagonal matrix', () => {
      const M = new Matrix(F7, 4, 4, (i, j) =>
        Math.abs(i - j) <= 1 ? F7.__call__(1n) : F7.__call__(0n)
      );

      const [lower, upper] = get_bandwidth(M);

      expect(lower).toBe(1);
      expect(upper).toBe(1);
    });

    it('should compute bandwidth of upper triangular matrix', () => {
      const M = new Matrix(F7, 3, 3, (i, j) => (i <= j ? F7.__call__(1n) : F7.__call__(0n)));

      const [lower, upper] = get_bandwidth(M);

      expect(lower).toBe(0);
      expect(upper).toBe(2);
    });
  });
});

// ===========================================================================
// Coverage added by the 2026-07 audit (finding M68): these exports previously
// had no test at all, which is how H52 (32-bit rook masks) and H53 (always
// throwing berlekamp_massey) survived.
// ===========================================================================

describe('matrix_special - rook_vector', () => {
  const asNumbers = (v: Rational[]) => v.map((x) => Number(x.toString()));
  const mkQ = (rows: number[][]) =>
    new Matrix(QQ, rows.length, rows[0]!.length, (i, j) => QQ.__call__(BigInt(rows[i]![j]!)));

  it('should reproduce the SageMath doctest for ones_matrix(8, 8)', () => {
    // sage: ones_matrix(8,8).rook_vector()
    // [1, 64, 1568, 18816, 117600, 376320, 564480, 322560, 40320]
    // This one exercises the automatic use_complement path and therefore the
    // inclusion-exclusion coefficients, which must be exact integers.
    expect(asNumbers(rook_vector(ones_matrix(QQ, 8, 8)))).toEqual([
      1, 64, 1568, 18816, 117600, 376320, 564480, 322560, 40320,
    ]);
  });

  it('should reproduce the SageMath doctest for the 3x6 band matrix', () => {
    const A = mkQ([
      [1, 1, 1, 1, 0, 0],
      [0, 1, 1, 1, 1, 0],
      [0, 0, 1, 1, 1, 1],
    ]);
    expect(asNumbers(rook_vector(A))).toEqual([1, 12, 40, 36]);
  });

  it('should agree across algorithms on the SageMath 4x4 doctest', () => {
    const A = mkQ([
      [1, 0, 0, 1],
      [0, 1, 1, 0],
      [0, 1, 1, 0],
      [1, 0, 0, 1],
    ]);
    expect(asNumbers(rook_vector(A, 'ButeraPernici'))).toEqual([1, 8, 20, 16, 4]);
    expect(asNumbers(rook_vector(A, 'Ryser'))).toEqual([1, 8, 20, 16, 4]);
  });

  it('should handle more than 32 columns (bignum Grassmann masks)', () => {
    // Regression test for H52: with 32-bit shifts, column 32 aliased column 0
    // and this returned [1, 2, 0].
    const A = new Matrix(QQ, 2, 33, (_i, j) =>
      j === 0 || j === 32 ? QQ.__call__(1n) : QQ.__call__(0n)
    );
    expect(asNumbers(rook_vector(A, 'ButeraPernici', false, false))).toEqual([1, 4, 2]);
    expect(asNumbers(rook_vector(A, 'Ryser', false, false))).toEqual([1, 4, 2]);
    expect(asNumbers(rook_vector(A, 'naive', false, false))).toEqual([1, 4, 2]);

    // all-ones 2 x 40 (previously returned the 2 x 32 answer [1, 64, 992])
    expect(asNumbers(rook_vector(ones_matrix(QQ, 2, 40), 'ButeraPernici', false, false))).toEqual([
      1, 80, 1560,
    ]);
    expect(asNumbers(rook_vector(ones_matrix(QQ, 2, 40), 'Ryser', false, false))).toEqual([
      1, 80, 1560,
    ]);
  });

  it('should count the derangements of 21 elements', () => {
    // sage: A = identity_matrix(21); A.rook_vector(complement=True)[-1]
    // 18795307255050944540
    // Regression test for M66: the inclusion-exclusion coefficients must be
    // exact BigInts and the scalar multiplication must be double-and-add.
    const I21 = new Matrix(QQ, 21, 21, (i, j) => (i === j ? QQ.__call__(1n) : QQ.__call__(0n)));
    const v = rook_vector(I21, 'ButeraPernici', true);
    expect(v[v.length - 1]!.toString()).toBe('18795307255050944540');
  });

  it('should reproduce the complement doctest for the 4x5 matrix', () => {
    const D = new Matrix(QQ, 4, 5, (i, j) =>
      i === j && j < 4 ? QQ.__call__(1n) : QQ.__call__(0n)
    );
    for (const algorithm of ['Ryser', 'ButeraPernici'] as const) {
      expect(asNumbers(rook_vector(D, algorithm, true, true))).toEqual([1, 16, 78, 128, 53]);
      expect(asNumbers(rook_vector(D, algorithm, true, false))).toEqual([1, 16, 78, 128, 53]);
    }
  });

  it('should reject a complement request on a non 0-1 matrix', () => {
    const A = mkQ([
      [1, 5],
      [0, 1],
    ]);
    expect(() => rook_vector(A, 'ButeraPernici', true)).toThrow(
      "coefficients must be zero or one, but we have '5' in position (0,1)."
    );
  });

  it('should return [1, 0, 0] for a zero matrix', () => {
    expect(asNumbers(rook_vector(new Matrix(QQ, 2, 2), 'ButeraPernici'))).toEqual([1, 0, 0]);
    expect(asNumbers(rook_vector(new Matrix(QQ, 2, 2), 'Ryser'))).toEqual([1, 0, 0]);
  });
});

describe('matrix_special - berlekamp_massey', () => {
  it('should reproduce the SageMath doctest over GF(7)', () => {
    // sage: berlekamp_massey([GF(7)(1), 19, 1, 19])
    // x^2 + 6
    const F7 = GF(7n);
    const g = berlekamp_massey([
      F7.__call__(1n),
      F7.__call__(19n),
      F7.__call__(1n),
      F7.__call__(19n),
    ]);
    expect(g.coeffs.map((c) => c.toString())).toEqual(['6', '0', '1']);
  });

  it('should reproduce the SageMath doctest over the rationals', () => {
    // sage: berlekamp_massey([1,2,1,2,1,2])
    // x^2 - 1
    expect(berlekamp_massey([1, 2, 1, 2, 1, 2]).coeffs.map((c) => c.toString())).toEqual([
      '-1',
      '0',
      '1',
    ]);

    // sage: berlekamp_massey([2,2,1,2,1,191,393,132])
    // x^4 - 36727/11711*x^3 + 34213/5019*x^2 + 7024942/35133*x - 335813/1673
    expect(
      berlekamp_massey([2, 2, 1, 2, 1, 191, 393, 132]).coeffs.map((c) => c.toString())
    ).toEqual(['-335813/1673', '7024942/35133', '34213/5019', '-36727/11711', '1']);

    // sage: berlekamp_massey(prime_range(2, 38))
    // x^6 - 14/9*x^5 - 7/9*x^4 + 157/54*x^3 - 25/27*x^2 - 73/18*x + 37/9
    expect(
      berlekamp_massey([2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]).coeffs.map((c) => c.toString())
    ).toEqual(['37/9', '-73/18', '-25/27', '157/54', '-7/9', '-14/9', '1']);
  });

  it('should reject an odd number of terms', () => {
    expect(() => berlekamp_massey([1, 2, 5])).toThrow('argument must have an even number of terms');
  });
});

describe('matrix_special - permutation normal form (exact)', () => {
  const mkQ = (rows: number[][]) =>
    new Matrix(QQ, rows.length, rows[0]!.length, (i, j) => QQ.__call__(BigInt(rows[i]![j]!)));
  const asNumbers = (M: Matrix<Rational>) =>
    M.rows().map((row) => row.map((x) => Number(x.toString())));

  it('should reproduce the SageMath doctests', () => {
    expect(
      asNumbers(
        permutation_normal_form(
          mkQ([
            [0, 0, 1],
            [1, 0, 2],
            [0, 0, 0],
          ])
        ) as Matrix<Rational>
      )
    ).toEqual([
      [2, 1, 0],
      [1, 0, 0],
      [0, 0, 0],
    ]);

    expect(
      asNumbers(
        permutation_normal_form(
          mkQ([
            [-1, 3],
            [-1, 5],
            [2, 4],
          ])
        ) as Matrix<Rational>
      )
    ).toEqual([
      [5, -1],
      [4, 2],
      [3, -1],
    ]);

    expect(
      asNumbers(
        permutation_normal_form(
          mkQ([
            [3, 4, 5],
            [3, 4, 5],
            [3, 5, 4],
            [2, 0, 1],
          ])
        ) as Matrix<Rational>
      )
    ).toEqual([
      [5, 4, 3],
      [5, 4, 3],
      [4, 5, 3],
      [1, 0, 2],
    ]);
  });

  it('should compare finite field elements by value, not by their printed form', () => {
    // Regression test for H55: "10" < "9" as strings, so the old code returned
    // the matrix unchanged here.
    const F101 = GF(101n);
    const vals = [
      [9, 0],
      [0, 10],
    ];
    const M = new Matrix(F101, 2, 2, (i, j) => F101.__call__(BigInt(vals[i]![j]!)));
    const nf = permutation_normal_form(M) as Matrix<ReturnType<typeof F101.__call__>>;
    expect(nf.rows().map((r) => r.map((x) => Number(x.value)))).toEqual([
      [10, 0],
      [0, 9],
    ]);
  });

  it('should agree with brute force on random 4x3 matrices over GF(11)', () => {
    const F11 = GF(11n);
    let rng = 987654321;
    const rnd = (n: number) => {
      rng = (rng * 1103515245 + 12345) % 2147483648;
      return rng % n;
    };
    const permsOf = (n: number): number[][] => {
      if (n === 0) return [[]];
      const out: number[][] = [];
      for (const p of permsOf(n - 1)) {
        for (let i = 0; i <= p.length; i++) {
          out.push([...p.slice(0, i), n - 1, ...p.slice(i)]);
        }
      }
      return out;
    };
    const rowPerms = permsOf(4);
    const colPerms = permsOf(3);
    const key = (rows: number[][]) =>
      rows.map((r) => r.map((x) => String(x).padStart(3, '0')).join(',')).join(';');

    for (let t = 0; t < 40; t++) {
      const A = new Matrix(F11, 4, 3, () => F11.__call__(BigInt(rnd(11))));
      let best: string | null = null;
      for (const rp of rowPerms) {
        for (const cp of colPerms) {
          const k = key(rp.map((r) => cp.map((c) => Number(A.get(r, c).value))));
          if (best === null || k > best) best = k;
        }
      }
      const nf = permutation_normal_form(A) as Matrix<ReturnType<typeof F11.__call__>>;
      expect(key(nf.rows().map((r) => r.map((x) => Number(x.value))))).toBe(best);
    }
  });

  it('should return a permutation that realises the normal form with check=true', () => {
    const M = mkQ([
      [-1, 3],
      [-1, 5],
      [2, 4],
    ]);
    const [nf, [rowPerm, colPerm]] = permutation_normal_form(M, true) as [
      Matrix<Rational>,
      [number[], number[]],
    ];
    const rebuilt = new Matrix(QQ, 3, 2, (i, j) => M.get(rowPerm[i]!, colPerm[j]!));
    expect(rebuilt.eq(nf)).toBe(true);
  });

  it('should detect permutations of 7x3 matrices over GF(11)', () => {
    // Regression test for H54: above 6x6 the old code compared a non-canonical
    // heuristic normal form and produced 220 false negatives in 300 trials.
    const F11 = GF(11n);
    let rng = 123456789;
    const rnd = (n: number) => {
      rng = (rng * 1103515245 + 12345) % 2147483648;
      return rng % n;
    };
    for (let trial = 0; trial < 100; trial++) {
      const A = new Matrix(F11, 7, 3, () => F11.__call__(BigInt(rnd(11))));
      const rowPerm = [0, 1, 2, 3, 4, 5, 6];
      for (let i = 6; i > 0; i--) {
        const j = rnd(i + 1);
        [rowPerm[i], rowPerm[j]] = [rowPerm[j]!, rowPerm[i]!];
      }
      const colPerm = [0, 1, 2];
      for (let i = 2; i > 0; i--) {
        const j = rnd(i + 1);
        [colPerm[i], colPerm[j]] = [colPerm[j]!, colPerm[i]!];
      }
      const B = new Matrix(F11, 7, 3, (i, j) => A.get(rowPerm[i]!, colPerm[j]!));
      expect(is_permutation_of(A, B)).toBe(true);

      // and a single perturbed entry must be rejected
      const C = B.copy();
      const i0 = rnd(7);
      const j0 = rnd(3);
      C.set(i0, j0, C.get(i0, j0).add(F11.__call__(1n)));
      expect(is_permutation_of(A, C)).toBe(false);
    }
  });

  it('should reproduce the SageMath is_permutation_of doctests', () => {
    const M = mkQ([
      [1, 2, 3],
      [3, 5, 3],
      [2, 6, 4],
    ]);
    expect(
      is_permutation_of(
        M,
        mkQ([
          [1, 2, 3],
          [2, 6, 4],
          [3, 5, 3],
        ])
      )
    ).toBe(true);
    expect(
      is_permutation_of(
        M,
        mkQ([
          [1, 2, 3],
          [4, 5, 6],
          [7, 8, 9],
        ])
      )
    ).toBe(false);
    expect(
      is_permutation_of(
        M,
        mkQ([
          [1, 2],
          [3, 4],
        ])
      )
    ).toBe(false);

    const N = mkQ([
      [3, 5, 3],
      [2, 6, 4],
      [1, 2, 3],
    ]);
    const [truth, perm] = is_permutation_of(M, N, true) as [boolean, [number[], number[]]];
    expect(truth).toBe(true);
    const rebuilt = new Matrix(QQ, 3, 3, (i, j) => M.get(perm[0][i]!, perm[1][j]!));
    expect(rebuilt.eq(N)).toBe(true);
  });
});

describe('matrix_special - random matrices', () => {
  it('should return determinant one unimodular matrices', () => {
    // Regression test for M61: the old implementation scaled rows by random
    // units and swapped rows, so the determinant was an arbitrary unit.
    const F101 = GF(101n);
    for (let t = 0; t < 20; t++) {
      expect(determinant(random_unimodular_matrix(F101, 5)).value).toBe(1n);
    }
    for (let t = 0; t < 10; t++) {
      expect(determinant(random_unimodular_matrix(QQ, 4)).toString()).toBe('1');
    }
  });

  it('should produce echelonizable matrices of the requested rank', () => {
    const F101 = GF(101n);
    for (const rank of [0, 1, 2, 3]) {
      const A = random_echelonizable_matrix(F101, 4, 5, rank);
      expect(A.nrows).toBe(4);
      expect(A.ncols).toBe(5);
      expect(matrix_rank(A)).toBe(rank);
    }
  });

  it('should refuse the unimplemented upper_bound size control', () => {
    expect(() => random_echelonizable_matrix(QQ, 3, 3, 3, 50)).toThrow('upper_bound');
  });

  it('should produce a genuinely non-diagonal diagonalizable matrix', () => {
    // Regression test for M60: the capability probe used to look for a method
    // named `inverse` (every element class provides `inv`), so the P*D*P^-1
    // code was dead and a plain diagonal matrix was always returned.
    const F3 = GF(3n);
    let offDiagonal = 0;
    for (let t = 0; t < 20; t++) {
      const M = random_diagonalizable_matrix(F3, 3);
      for (let i = 0; i < 3; i++) {
        for (let j = 0; j < 3; j++) {
          if (i !== j && !M.get(i, j).isZero()) offDiagonal++;
        }
      }
    }
    expect(offDiagonal).toBeGreaterThan(0);
  });

  it('should honour the requested eigenvalues and multiplicities', () => {
    for (let t = 0; t < 10; t++) {
      const M = random_diagonalizable_matrix(
        QQ,
        6,
        [-12, 4, 6].map((x) => QQ.__call__(BigInt(x))),
        [2, 3, 1]
      );
      // trace = sum of the eigenvalues with multiplicity
      expect(M.trace().toString()).toBe(String(-12 * 2 + 4 * 3 + 6 * 1));
      for (const e of [-12, 4, 6]) {
        const shifted = new Matrix(QQ, 6, 6, (i, j) =>
          i === j ? M.get(i, j).sub(QQ.__call__(BigInt(e))) : M.get(i, j)
        );
        expect(determinant(shifted).isZero()).toBe(true);
      }
    }
  });

  it('should validate the eigenvalue arguments as SageMath does', () => {
    expect(() => random_diagonalizable_matrix(QQ, 3, [QQ.__call__(1n)])).toThrow(
      'the list of eigenvalues must have a list of dimensions corresponding to each eigenvalue.'
    );
    expect(() => random_diagonalizable_matrix(QQ, 3, undefined, [1, 2])).toThrow(
      'the list of dimensions must have a list of corresponding eigenvalues.'
    );
    expect(() =>
      random_diagonalizable_matrix(QQ, 3, [QQ.__call__(1n), QQ.__call__(2n)], [1, 1])
    ).toThrow('the size of the matrix must equal the sum of the dimensions.');
    expect(() =>
      random_diagonalizable_matrix(QQ, 3, [QQ.__call__(1n), QQ.__call__(2n)], [3, 0])
    ).toThrow('eigenspaces must have a dimension of at least 1.');
  });

  it('should randomize a matrix in place', () => {
    const F101 = GF(101n);
    const M = new Matrix(F101, 4, 4);
    randomize(M, 1.0, true);
    for (let i = 0; i < 4; i++) {
      for (let j = 0; j < 4; j++) {
        expect(M.get(i, j).isZero()).toBe(false);
      }
    }
  });
});

describe('matrix_special - miscellaneous constructors and operations', () => {
  it('should build the Hilbert matrix (SageMath doctest)', () => {
    // sage: matrix.hilbert(5)
    // [  1 1/2 1/3 1/4 1/5]  ...
    const H = hilbert(5, QQ);
    expect(H.rows().map((r) => r.map((x) => x.toString()))).toEqual([
      ['1', '1/2', '1/3', '1/4', '1/5'],
      ['1/2', '1/3', '1/4', '1/5', '1/6'],
      ['1/3', '1/4', '1/5', '1/6', '1/7'],
      ['1/4', '1/5', '1/6', '1/7', '1/8'],
      ['1/5', '1/6', '1/7', '1/8', '1/9'],
    ]);
  });

  it('should build the Lehmer matrix (SageMath doctest)', () => {
    // sage: matrix.lehmer(3)
    // [  1 1/2 1/3]
    // [1/2   1 2/3]
    // [1/3 2/3   1]
    const L = lehmer(QQ, 3);
    expect(L.rows().map((r) => r.map((x) => x.toString()))).toEqual([
      ['1', '1/2', '1/3'],
      ['1/2', '1', '2/3'],
      ['1/3', '2/3', '1'],
    ]);
  });

  it('should compute prod_of_row_sums (SageMath doctests)', () => {
    // sage: a = matrix(QQ, 2, 2, [1,2,3,2]); a.prod_of_row_sums([0,1]) == 15
    const a = new Matrix(
      QQ,
      2,
      2,
      [1, 2, 3, 2].map((x) => QQ.__call__(BigInt(x)))
    );
    expect(prod_of_row_sums(a, [0, 1]).toString()).toBe('15');

    // sage: a = matrix(QQ, 2, 3, [1,2,3,2,5,6]); a.prod_of_row_sums([1,2]) == 55
    const b = new Matrix(
      QQ,
      2,
      3,
      [1, 2, 3, 2, 5, 6].map((x) => QQ.__call__(BigInt(x)))
    );
    expect(prod_of_row_sums(b, [1, 2]).toString()).toBe('55');
  });

  it('should compute the denominator of a rational matrix', () => {
    // sage: matrix(QQ, 2, [1/2, 1/3, 1/5, 1]).denominator() == 30
    const M = new Matrix(QQ, 2, 2, [
      QQ.__call__(1n, 2n),
      QQ.__call__(1n, 3n),
      QQ.__call__(1n, 5n),
      QQ.__call__(1n),
    ]);
    expect(denominator(M).toString()).toBe('30');
  });

  it('should apply a ring morphism entrywise', () => {
    const F7 = GF(7n);
    const M = new Matrix(
      F7,
      2,
      2,
      [0, 1, 2, 3].map((x) => F7.__call__(BigInt(x)))
    );
    const squared = apply_morphism(M, (x) => x.mul(x));
    expect(squared.rows().map((r) => r.map((x) => Number(x.value)))).toEqual([
      [0, 1],
      [4, 2],
    ]);
  });

  it('should substitute and differentiate polynomial entries', () => {
    const F7 = GF(7n);
    const R = new PolynomialRing(F7, 'x');
    const x = R.gen();
    const one = R.one();
    // [[x^2, x + 1], [1, x]]
    const M = new Matrix(R, 2, 2, [x.mul(x), x.add(one), one, x]);

    const D = derivative(M);
    expect(D.rows().map((r) => r.map((p) => p.toString()))).toEqual([
      ['2*x', '1'],
      ['0', '1'],
    ]);

    const S = subs(M, F7.__call__(3n));
    expect(S.rows().map((r) => r.map((p) => p.toString()))).toEqual([
      ['2', '4'],
      ['1', '3'],
    ]);
  });

  it('should reject numerical_approx and hadamard_bound over rings without them', () => {
    const F7 = GF(7n);
    const M = new Matrix(
      F7,
      2,
      2,
      [1, 2, 3, 4].map((x) => F7.__call__(BigInt(x)))
    );
    expect(() => numerical_approx(M)).toThrow('numerical_approx');
    expect(() => hadamard_bound(M)).toThrow('sqrt');
  });
});
