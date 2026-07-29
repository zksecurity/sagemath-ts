/**
 * Tests for Z/nZ matrix operations.
 *
 * @see Reference: sage/matrix/matrix_modn_dense_template.pxi
 */

import { describe, expect, it } from 'vitest';
import {
  Matrix_modn_dense,
  identity_matrix_modn,
  matrix_modn_from_entries,
  zero_matrix_modn,
} from './matrix_modn.js';

/** Evaluate p(A) for a polynomial given constant-term-first. */
function polyEval(A: Matrix_modn_dense, coeffs: bigint[]): Matrix_modn_dense {
  const n = A.nrows;
  let result = new Matrix_modn_dense(n, n, A.modulus);
  let power = identity_matrix_modn(A.modulus, n);
  for (const c of coeffs) {
    result = result.add(power.scalar_mul(c));
    power = power.mul(A);
  }
  return result;
}

/** Deterministic LCG so the sweeps are reproducible. */
function makeRandom(seed: number): (m: number) => number {
  let state = seed;
  return (m: number) => {
    state = (state * 1103515245 + 12345) & 0x7fffffff;
    return state % m;
  };
}

describe('Matrix_modn_dense.charpoly', () => {
  it('computes det(x*I - A) for a diagonal matrix', () => {
    // charpoly(diag(2,3)) = (x-2)(x-3) = x^2 - 5x + 6
    const A = matrix_modn_from_entries(101n, [
      [2, 0],
      [0, 3],
    ]);
    expect(A.charpoly()).toEqual([6n, 96n, 1n]);
  });

  it('computes the charpoly of a Jordan block', () => {
    // charpoly([[1,1],[0,1]]) = (x-1)^2 = x^2 - 2x + 1
    const A = matrix_modn_from_entries(101n, [
      [1, 1],
      [0, 1],
    ]);
    expect(A.charpoly()).toEqual([1n, 99n, 1n]);
  });

  it('works over GF(2), where Faddeev-LeVerrier would divide by zero', () => {
    const A = matrix_modn_from_entries(2n, [
      [1, 1],
      [0, 1],
    ]);
    // (x-1)^2 = x^2 + 1 over GF(2)
    expect(A.charpoly()).toEqual([1n, 0n, 1n]);
  });

  it('works over GF(3) for a 4x4 matrix (i+1 hits the modulus)', () => {
    const A = matrix_modn_from_entries(3n, [
      [1, 2, 0, 1],
      [0, 1, 1, 0],
      [2, 0, 1, 1],
      [1, 1, 0, 2],
    ]);
    const cp = A.charpoly();
    expect(cp.length).toBe(5);
    expect(cp[4]).toBe(1n);
    // constant term is (-1)^n * det
    expect(cp[0]).toBe(A.determinant());
  });

  it('works over a composite modulus (Sage falls back to a generic path)', () => {
    const A = matrix_modn_from_entries(4n, [
      [1, 2, 3, 0],
      [2, 1, 0, 3],
      [3, 0, 1, 2],
      [0, 3, 2, 1],
    ]);
    const cp = A.charpoly();
    expect(cp.length).toBe(5);
    expect(cp[0]).toBe(A.determinant());
  });

  it('satisfies Cayley-Hamilton for random matrices over several moduli', () => {
    const rnd = makeRandom(12345);
    for (const p of [2n, 3n, 5n, 7n, 101n, 4n, 6n, 9n]) {
      for (let trial = 0; trial < 20; trial++) {
        const n = 1 + rnd(5);
        const entries: number[][] = [];
        for (let i = 0; i < n; i++) {
          entries.push([]);
          for (let j = 0; j < n; j++) entries[i]!.push(rnd(Number(p)));
        }
        const A = matrix_modn_from_entries(p, entries);
        const cp = A.charpoly();
        expect(cp.length).toBe(n + 1);
        expect(cp[n]).toBe(1n);
        // p(A) == 0
        expect(polyEval(A, cp).eq(new Matrix_modn_dense(n, n, p))).toBe(true);
        // constant coefficient is (-1)^n det(A)
        const det = A.determinant();
        const expected = n % 2 === 0 ? det : (p - det) % p;
        expect(cp[0]).toBe(expected);
      }
    }
  });

  it('rejects non-square matrices', () => {
    const A = matrix_modn_from_entries(7n, [
      [1, 2, 3],
      [4, 5, 6],
    ]);
    expect(() => A.charpoly()).toThrow();
  });
});

describe('Matrix_modn_dense.determinant', () => {
  it('handles n >= 4 without requiring an invertible pivot', () => {
    // Over Z/4 the leading pivots are not invertible; Sage lifts to ZZ.
    const A = matrix_modn_from_entries(4n, [
      [2, 2, 0, 0],
      [2, 2, 0, 0],
      [0, 0, 2, 2],
      [0, 0, 2, 2],
    ]);
    expect(A.determinant()).toBe(0n);
  });

  it('agrees with cofactor expansion over Z/6', () => {
    const entries = [
      [1, 2, 3, 4],
      [5, 0, 1, 2],
      [3, 3, 3, 3],
      [1, 0, 0, 5],
    ];
    const A = matrix_modn_from_entries(6n, entries);
    // exact determinant over ZZ by cofactor expansion
    const det = (M: number[][]): bigint => {
      const n = M.length;
      if (n === 1) return BigInt(M[0]![0]!);
      let s = 0n;
      for (let j = 0; j < n; j++) {
        const minor = M.slice(1).map((row) => row.filter((_, c) => c !== j));
        s += (j % 2 === 0 ? 1n : -1n) * BigInt(M[0]![j]!) * det(minor);
      }
      return s;
    };
    const exact = ((det(entries) % 6n) + 6n) % 6n;
    expect(A.determinant()).toBe(exact);
  });

  it('is multiplicative over GF(7)', () => {
    const rnd = makeRandom(99);
    for (let trial = 0; trial < 20; trial++) {
      const n = 4;
      const mk = (): Matrix_modn_dense => {
        const e: number[][] = [];
        for (let i = 0; i < n; i++) {
          e.push([]);
          for (let j = 0; j < n; j++) e[i]!.push(rnd(7));
        }
        return matrix_modn_from_entries(7n, e);
      };
      const A = mk();
      const B = mk();
      expect(A.mul(B).determinant()).toBe((A.determinant() * B.determinant()) % 7n);
      expect(A.transpose().determinant()).toBe(A.determinant());
    }
  });
});

describe('Matrix_modn_dense.echelonize', () => {
  it('raises NotImplementedError over a composite modulus', () => {
    // sage: matrix(Integers(4), 2, [2,2,2,2]).rank()
    // NotImplementedError: Echelon form not implemented over 'Ring of integers modulo 4'.
    const A = matrix_modn_from_entries(4n, [
      [2, 2],
      [2, 2],
    ]);
    expect(() => A.echelon_form()).toThrow(
      "Echelon form not implemented over 'Ring of integers modulo 4'."
    );
    expect(() => A.rank()).toThrow(
      "Echelon form not implemented over 'Ring of integers modulo 4'."
    );
  });

  it('reproduces Sage GF(97) doctest', () => {
    // sage: A = matrix(GF(97),3,4,range(12)); A.echelonize(); A
    // [ 1  0 96 95]
    // [ 0  1  2  3]
    // [ 0  0  0  0]
    const A = matrix_modn_from_entries(97n, [
      [0, 1, 2, 3],
      [4, 5, 6, 7],
      [8, 9, 10, 11],
    ]);
    const E = A.echelon_form();
    expect([E.get(0, 0), E.get(0, 1), E.get(0, 2), E.get(0, 3)]).toEqual([1n, 0n, 96n, 95n]);
    expect([E.get(1, 0), E.get(1, 1), E.get(1, 2), E.get(1, 3)]).toEqual([0n, 1n, 2n, 3n]);
    expect([E.get(2, 0), E.get(2, 1), E.get(2, 2), E.get(2, 3)]).toEqual([0n, 0n, 0n, 0n]);
    expect(A.pivots()).toEqual([0, 1]);
  });
});

describe('Matrix_modn_dense.right_kernel_matrix', () => {
  it('reproduces the Sage doctest over GF(5)', () => {
    // sage: M = matrix(GF(5),6,6,range(36))
    const entries: number[][] = [];
    for (let i = 0; i < 6; i++) {
      entries.push([]);
      for (let j = 0; j < 6; j++) entries[i]!.push(i * 6 + j);
    }
    const M = matrix_modn_from_entries(5n, entries);

    const rows = (A: Matrix_modn_dense): bigint[][] => {
      const r: bigint[][] = [];
      for (let i = 0; i < A.nrows; i++) {
        const row: bigint[] = [];
        for (let j = 0; j < A.ncols; j++) row.push(A.get(i, j));
        r.push(row);
      }
      return r;
    };

    expect(rows(M.right_kernel_matrix({ basis: 'computed' }))).toEqual([
      [4n, 2n, 4n, 0n, 0n, 0n],
      [3n, 3n, 0n, 4n, 0n, 0n],
      [2n, 4n, 0n, 0n, 4n, 0n],
      [1n, 0n, 0n, 0n, 0n, 4n],
    ]);
    expect(rows(M.right_kernel_matrix({ basis: 'pivot' }))).toEqual([
      [1n, 3n, 1n, 0n, 0n, 0n],
      [2n, 2n, 0n, 1n, 0n, 0n],
      [3n, 1n, 0n, 0n, 1n, 0n],
      [4n, 0n, 0n, 0n, 0n, 1n],
    ]);
    expect(rows(M.right_kernel_matrix())).toEqual([
      [1n, 0n, 0n, 0n, 0n, 4n],
      [0n, 1n, 0n, 0n, 1n, 3n],
      [0n, 0n, 1n, 0n, 2n, 2n],
      [0n, 0n, 0n, 1n, 3n, 1n],
    ]);

    // sage: M * M.right_kernel_matrix().transpose() == 0
    const K = M.right_kernel_matrix();
    expect(M.mul(K.transpose()).eq(new Matrix_modn_dense(6, 4, 5n))).toBe(true);
  });

  it('annihilates a rank-deficient 2x2 over GF(7)', () => {
    const A = matrix_modn_from_entries(7n, [
      [1, 2],
      [3, 6],
    ]);
    const K = A.right_kernel_matrix();
    expect(K.nrows).toBe(1);
    expect(A.mul(K.transpose()).eq(new Matrix_modn_dense(2, 1, 7n))).toBe(true);
  });

  it('handles a 1x2 matrix without throwing', () => {
    const A = matrix_modn_from_entries(7n, [[1, 1]]);
    const K = A.right_kernel_matrix();
    expect(K.nrows).toBe(1);
    expect(A.mul(K.transpose()).eq(new Matrix_modn_dense(1, 1, 7n))).toBe(true);
  });

  it('satisfies A*K^T = 0 and dim K = ncols - rank for random matrices', () => {
    const rnd = makeRandom(2468);
    for (let trial = 0; trial < 200; trial++) {
      const p = [2n, 3n, 5n, 7n, 11n, 101n][rnd(6)]!;
      const m = 1 + rnd(5);
      const n = 1 + rnd(5);
      const entries: number[][] = [];
      for (let i = 0; i < m; i++) {
        entries.push([]);
        for (let j = 0; j < n; j++) entries[i]!.push(rnd(Number(p)));
      }
      const A = matrix_modn_from_entries(p, entries);
      const K = A.right_kernel_matrix();
      expect(K.nrows).toBe(n - A.rank());
      expect(A.mul(K.transpose()).eq(new Matrix_modn_dense(m, K.nrows, p))).toBe(true);
    }
  });

  it('rejects an unknown basis format', () => {
    const A = matrix_modn_from_entries(7n, [[1, 1]]);
    // biome-ignore lint/suspicious/noExplicitAny: deliberately invalid input
    expect(() => A.right_kernel_matrix({ basis: 'bogus' as any })).toThrow(
      'matrix kernel basis format not recognized'
    );
  });
});

describe('Matrix_modn_dense.right_kernel_matrix over Z/nZ with composite n', () => {
  const rowsOf = (M: Matrix_modn_dense): bigint[][] =>
    Array.from({ length: M.nrows }, (_, i) =>
      Array.from({ length: M.ncols }, (_, j) => M.get(i, j))
    );

  it('reproduces the Sage doctest for Zmod(24480)', () => {
    // sage: A = matrix(Zmod(24480), [[1,2,3,4,5], [7,7,7,7,7]])
    // sage: A._right_kernel_matrix_over_integer_mod_ring()[1]
    // [    1 24478     1     0     0]
    // [    2 24477     0     1     0]
    // [    3 24476     0     0     1]
    const A = matrix_modn_from_entries(24480n, [
      [1, 2, 3, 4, 5],
      [7, 7, 7, 7, 7],
    ]);
    const K = A.right_kernel_matrix();
    expect(rowsOf(K)).toEqual([
      [1n, 24478n, 1n, 0n, 0n],
      [2n, 24477n, 0n, 1n, 0n],
      [3n, 24476n, 0n, 0n, 1n],
    ]);
    // sage: A * P.transpose() == 0
    expect(A.mul(K.transpose()).eq(new Matrix_modn_dense(2, 3, 24480n))).toBe(true);
  });

  it('matches SageMath basis for basis on a range of composite moduli', () => {
    // Golden values produced by SageMath 10.3:
    //   matrix(Zmod(d), A).right_kernel_matrix()
    const golden: [bigint, number[][], bigint[][]][] = [
      [
        4n,
        [
          [3, 1, 0, 0],
          [2, 0, 5, 6],
        ],
        [
          [2n, 2n, 0n, 0n],
          [1n, 1n, 2n, 0n],
          [1n, 1n, 0n, 1n],
        ],
      ],
      [
        6n,
        [
          [5, 5, 0, 0, 3],
          [0, 1, 4, 0, 0],
          [4, 0, 3, 1, 1],
        ],
        [
          [2n, 4n, 5n, 1n, 0n],
          [5n, 4n, 5n, 0n, 1n],
        ],
      ],
      [8n, [[6], [4]], [[4n]]],
      [
        9n,
        [[4, 1, 4, 5]],
        [
          [2n, 1n, 0n, 0n],
          [8n, 0n, 1n, 0n],
          [1n, 0n, 0n, 1n],
        ],
      ],
      [
        12n,
        [
          [4, 1],
          [2, 0],
          [2, 5],
        ],
        [[6n, 0n]],
      ],
      [
        24n,
        [[6, 4, 3, 6, 4]],
        [
          [4n, 0n, 0n, 0n, 0n],
          [2n, 3n, 0n, 0n, 0n],
          [3n, 0n, 2n, 0n, 0n],
          [3n, 0n, 0n, 1n, 0n],
          [2n, 2n, 0n, 0n, 1n],
        ],
      ],
      [36n, [[4]], [[9n]]],
      [
        100n,
        [
          [6, 0, 4, 0, 2],
          [1, 0, 4, 0, 5],
          [0, 5, 3, 3, 1],
          [6, 3, 6, 3, 6],
        ],
        [
          [20n, 15n, 70n, 5n, 0n],
          [63n, 56n, 3n, 2n, 5n],
        ],
      ],
      [
        360n,
        [
          [5, 3],
          [3, 5],
        ],
        [[45n, 45n]],
      ],
      [24480n, [[13529, 21667]], [[8197n, 1n]]],
      [
        65536n,
        [[59247, 46700, 56490, 37079]],
        [
          [26028n, 1n, 0n, 0n],
          [63242n, 0n, 1n, 0n],
          [1767n, 0n, 0n, 1n],
        ],
      ],
      [
        1000000n,
        [[1, 5, 4, 1]],
        [
          [999995n, 1n, 0n, 0n],
          [999996n, 0n, 1n, 0n],
          [999999n, 0n, 0n, 1n],
        ],
      ],
    ];

    for (const [d, entries, expected] of golden) {
      const A = matrix_modn_from_entries(d, entries);
      const K = A.right_kernel_matrix();
      expect(rowsOf(K)).toEqual(expected);
      expect(A.mul(K.transpose()).eq(new Matrix_modn_dense(A.nrows, K.nrows, d))).toBe(true);
    }
  });

  it('spans exactly the brute-force kernel for small composite moduli', () => {
    /** All v in (Z/d)^n with A*v = 0, as "a,b,c" keys. */
    const bruteKernel = (A: bigint[][], m: number, n: number, d: bigint): Set<string> => {
      const out = new Set<string>();
      const total = d ** BigInt(n);
      for (let t = 0n; t < total; t++) {
        const v: bigint[] = [];
        let x = t;
        for (let j = 0; j < n; j++) {
          v.push(x % d);
          x /= d;
        }
        let ok = true;
        for (let i = 0; i < m && ok; i++) {
          let s = 0n;
          for (let j = 0; j < n; j++) s += A[i]![j]! * v[j]!;
          if (s % d !== 0n) ok = false;
        }
        if (ok) out.add(v.join(','));
      }
      return out;
    };

    /** The Z/d-submodule generated by `gens` (closure under adding generators). */
    const spanOf = (gens: bigint[][], n: number, d: bigint): Set<string> => {
      const seen = new Set<string>([new Array(n).fill(0n).join(',')]);
      let frontier: bigint[][] = [new Array<bigint>(n).fill(0n)];
      while (frontier.length > 0) {
        const next: bigint[][] = [];
        for (const v of frontier) {
          for (const g of gens) {
            const w = v.map((x, j) => (x + g[j]!) % d);
            const key = w.join(',');
            if (!seen.has(key)) {
              seen.add(key);
              next.push(w);
            }
          }
        }
        frontier = next;
      }
      return seen;
    };

    const rnd = makeRandom(13579);
    let checked = 0;
    for (const d of [4n, 6n, 8n, 9n, 12n]) {
      for (let m = 1; m <= 3; m++) {
        for (let n = 1; n <= 3; n++) {
          for (let trial = 0; trial < 6; trial++) {
            const entries: number[][] = [];
            for (let i = 0; i < m; i++) {
              entries.push([]);
              for (let j = 0; j < n; j++) entries[i]!.push(rnd(Number(d)));
            }
            const A = matrix_modn_from_entries(d, entries);
            const K = A.right_kernel_matrix();
            expect(K.ncols).toBe(n);
            const basis = rowsOf(K);
            const raw = entries.map((r) => r.map((x) => BigInt(x)));
            // Every basis vector is in the kernel.
            expect(A.mul(K.transpose()).eq(new Matrix_modn_dense(m, K.nrows, d))).toBe(true);
            // ... and the module they generate is the whole kernel, not a submodule.
            const brute = bruteKernel(raw, m, n, d);
            const span = spanOf(basis, n, d);
            expect(span.size).toBe(brute.size);
            for (const key of brute) expect(span.has(key)).toBe(true);
            checked++;
          }
        }
      }
    }
    expect(checked).toBe(270);
  });

  it('handles degenerate shapes', () => {
    // 0 columns: the only kernel vector is the empty one.
    const A = new Matrix_modn_dense(2, 0, 6n);
    const K = A.right_kernel_matrix();
    expect(K.nrows).toBe(0);
    expect(K.ncols).toBe(0);

    // 0 rows: everything is in the kernel, so the basis is the identity.
    const B = new Matrix_modn_dense(0, 3, 6n);
    const KB = B.right_kernel_matrix();
    expect(rowsOf(KB)).toEqual([
      [1n, 0n, 0n],
      [0n, 1n, 0n],
      [0n, 0n, 1n],
    ]);

    // Zero matrix: the kernel is everything.
    const Z = new Matrix_modn_dense(2, 2, 4n);
    expect(rowsOf(Z.right_kernel_matrix())).toEqual([
      [1n, 0n],
      [0n, 1n],
    ]);

    // Invertible matrix over Z/nZ: trivial kernel, but ncols must be preserved
    // (matrix2.pyx:4924 patches up PARI's zero-column empty matrix).
    const U = matrix_modn_from_entries(9n, [
      [1, 0],
      [0, 1],
    ]);
    const KU = U.right_kernel_matrix();
    expect(KU.nrows).toBe(0);
    expect(KU.ncols).toBe(2);
  });

  it('returns the full kernel of a tall zero matrix (PARI matkermod bug)', () => {
    // PARI's matkermod shortcuts m > 2n through matimagemod of the transpose
    // (bb_hnf.c:1049).  For the zero matrix the image is empty, and a t_MAT with
    // zero columns has no row count, so PARI (and hence SageMath 10.3) returns
    // an EMPTY kernel basis:
    //   sage: matrix(Zmod(4), 3, 1, [0,0,0]).right_kernel_matrix()
    //   []
    // which does not span the kernel (the whole of Z/4).  We request the image
    // from matkermod, which disables that branch, and return the true kernel.
    const A = new Matrix_modn_dense(3, 1, 4n);
    expect(rowsOf(A.right_kernel_matrix())).toEqual([[1n]]);

    const B = new Matrix_modn_dense(5, 2, 6n);
    expect(rowsOf(B.right_kernel_matrix())).toEqual([
      [1n, 0n],
      [0n, 1n],
    ]);

    // Sanity: below the m > 2n threshold SageMath itself returns the identity,
    // so the two agree there.
    //   sage: matrix(Zmod(4), 2, 2, [0,0,0,0]).right_kernel_matrix()
    //   [1 0]
    //   [0 1]
    const C = new Matrix_modn_dense(2, 2, 4n);
    expect(rowsOf(C.right_kernel_matrix())).toEqual([
      [1n, 0n],
      [0n, 1n],
    ]);
  });
});

describe('Matrix_modn_dense.minpoly', () => {
  it('annihilates the matrix and divides the charpoly', () => {
    const rnd = makeRandom(777);
    // 2 is excluded: `matrix_modn_dense_template.pxi:1570` forces
    // `algorithm='generic'` when `self.p == 2`, and :1589 raises
    // `NotImplementedError: Minimal polynomials are not implemented for Z/nZ.`
    // Verified against SageMath 10.3 (`matrix(GF(2),3,3,...).minpoly()` raises).
    // The refusal itself is asserted in the dedicated test below.
    for (const p of [3n, 5n, 7n, 101n]) {
      for (let trial = 0; trial < 10; trial++) {
        const n = 1 + rnd(4);
        const entries: number[][] = [];
        for (let i = 0; i < n; i++) {
          entries.push([]);
          for (let j = 0; j < n; j++) entries[i]!.push(rnd(Number(p)));
        }
        const A = matrix_modn_from_entries(p, entries);
        const mp = A.minpoly();
        expect(polyEval(A, mp).eq(new Matrix_modn_dense(n, n, p))).toBe(true);
      }
    }
  });

  it('refuses characteristic 2 and composite moduli, as SageMath does', () => {
    expect(() =>
      matrix_modn_from_entries(2n, [
        [1, 1, 0],
        [0, 1, 1],
        [1, 0, 1],
      ]).minpoly()
    ).toThrow('Minimal polynomials are not implemented for Z/nZ.');
    expect(() =>
      matrix_modn_from_entries(6n, [
        [1, 2],
        [3, 4],
      ]).minpoly()
    ).toThrow('Minimal polynomials are not implemented for Z/nZ.');
  });

  it('computes the minimal polynomial of a scalar matrix', () => {
    const S = matrix_modn_from_entries(7n, [
      [3, 0],
      [0, 3],
    ]);
    expect(S.minpoly()).toEqual([4n, 1n]); // x - 3
  });
});

describe('factory functions', () => {
  it('creates zero and identity matrices', () => {
    const Z = zero_matrix_modn(7n, 2, 3);
    expect(Z.nrows).toBe(2);
    expect(Z.ncols).toBe(3);
    expect(Z.get(1, 2)).toBe(0n);

    const I = identity_matrix_modn(7n, 3);
    expect(I.get(0, 0)).toBe(1n);
    expect(I.get(0, 1)).toBe(0n);
    expect(I.determinant()).toBe(1n);
  });
});
