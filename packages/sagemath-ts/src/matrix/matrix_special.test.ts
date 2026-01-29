/**
 * Tests for matrix_special functions
 */

import { describe, expect, it } from 'vitest';
import { GF } from '../rings/finite_rings/finite_field_constructor.js';
import { IntegerRing } from '../rings/integer_ring.js';
import { Matrix } from './matrix_generic.js';
import {
  subdivide,
  subdivision,
  subdivisions,
  as_bipartite_graph,
  automorphisms_of_rows_and_columns,
  permutation_normal_form,
  is_permutation_of,
  block_matrix,
  block_diagonal_matrix,
  jordan_block,
  companion_matrix,
  circulant,
  toeplitz,
  hankel,
  vandermonde,
  ones_matrix,
  column_matrix,
  elementary_matrix,
  tensor_product,
  elementwise_product,
  set_block,
  find,
  apply_map,
  get_bandwidth,
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

      subdivide(M, [[1, 2], [2, 3]] as [number[], number[]]);

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

      const normal = permutation_normal_form(M) as Matrix<typeof M extends Matrix<infer R> ? R : never>;

      // The maximal permutation should have largest elements first
      // 4, 3 in first row, 2, 1 in second (or some permutation giving max lex order)
      expect(normal.get(0, 0).value).toBeGreaterThanOrEqual(normal.get(1, 0).value);
    });

    it('should return permutation with check=true', () => {
      const M = new Matrix(F7, 2, 2, [
        [F7.__call__(1n), F7.__call__(2n)],
        [F7.__call__(3n), F7.__call__(4n)],
      ]);

      const result = permutation_normal_form(M, true) as [Matrix<typeof M extends Matrix<infer R> ? R : never>, [number[], number[]]];

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

      const normal = permutation_normal_form(M) as Matrix<typeof M extends Matrix<infer R> ? R : never>;

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
    it('should create companion matrix for cubic polynomial', () => {
      // Polynomial x^3 - c2*x^2 - c1*x - c0
      // Coefficients [c0, c1, c2]
      const coeffs = [F7.__call__(1n), F7.__call__(2n), F7.__call__(3n)];

      const C = companion_matrix(F7, coeffs);

      expect(C.nrows).toBe(3);
      expect(C.ncols).toBe(3);

      // Right companion matrix form:
      // [0 0 c0]
      // [1 0 c1]
      // [0 1 c2]
      expect(C.get(0, 2).value).toBe(1n); // c0
      expect(C.get(1, 2).value).toBe(2n); // c1
      expect(C.get(2, 2).value).toBe(3n); // c2

      expect(C.get(1, 0).value).toBe(1n);
      expect(C.get(2, 1).value).toBe(1n);
    });

    it('should create left companion matrix', () => {
      const coeffs = [F7.__call__(1n), F7.__call__(2n), F7.__call__(3n)];

      const C = companion_matrix(F7, coeffs, 'left');

      expect(C.nrows).toBe(3);
      expect(C.ncols).toBe(3);

      // Left companion matrix form:
      // [c2 c1 c0]
      // [1  0  0 ]
      // [0  1  0 ]
      expect(C.get(0, 0).value).toBe(3n); // c2
      expect(C.get(0, 1).value).toBe(2n); // c1
      expect(C.get(0, 2).value).toBe(1n); // c0

      expect(C.get(1, 0).value).toBe(1n);
      expect(C.get(2, 1).value).toBe(1n);
    });

    it('should handle empty coefficients', () => {
      const C = companion_matrix(F7, []);

      expect(C.nrows).toBe(0);
      expect(C.ncols).toBe(0);
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
    it('should create Toeplitz matrix', () => {
      const c = [F7.__call__(1n), F7.__call__(2n), F7.__call__(3n)]; // First column
      const r = [F7.__call__(1n), F7.__call__(4n), F7.__call__(5n)]; // First row

      const T = toeplitz(F7, c, r);

      expect(T.nrows).toBe(3);
      expect(T.ncols).toBe(3);

      // Entry (i,j) = c[i-j] if i >= j, else r[j-i]
      expect(T.get(0, 0).value).toBe(1n); // c[0]
      expect(T.get(0, 1).value).toBe(4n); // r[1]
      expect(T.get(0, 2).value).toBe(5n); // r[2]
      expect(T.get(1, 0).value).toBe(2n); // c[1]
      expect(T.get(1, 1).value).toBe(1n); // c[0]
      expect(T.get(1, 2).value).toBe(4n); // r[1]
      expect(T.get(2, 0).value).toBe(3n); // c[2]
      expect(T.get(2, 1).value).toBe(2n); // c[1]
      expect(T.get(2, 2).value).toBe(1n); // c[0]
    });
  });

  describe('hankel', () => {
    it('should create Hankel matrix', () => {
      const c = [F7.__call__(1n), F7.__call__(2n), F7.__call__(3n)];
      const r = [F7.__call__(3n), F7.__call__(4n), F7.__call__(5n)]; // Last row, shares c[2]

      const H = hankel(F7, c, r);

      expect(H.nrows).toBe(3);
      expect(H.ncols).toBe(3);

      // Hankel: constant on anti-diagonals
      // Entry (i,j) depends on i+j
      expect(H.get(0, 0).value).toBe(1n); // c[0]
      expect(H.get(0, 1).value).toBe(2n); // c[1]
      expect(H.get(0, 2).value).toBe(3n); // c[2]
      expect(H.get(1, 0).value).toBe(2n); // c[1]
      expect(H.get(1, 1).value).toBe(3n); // c[2]
      expect(H.get(1, 2).value).toBe(4n); // r[1]
      expect(H.get(2, 0).value).toBe(3n); // c[2]
      expect(H.get(2, 1).value).toBe(4n); // r[1]
      expect(H.get(2, 2).value).toBe(5n); // r[2]
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
      const matching = find(M, (x) => x.value > 4n) as Array<typeof M extends Matrix<infer R> ? R : never>;

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
