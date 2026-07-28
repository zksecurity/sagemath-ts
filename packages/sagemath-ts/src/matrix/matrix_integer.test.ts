/**
 * Tests for integer (ZZ) matrix normal forms and lattice reduction.
 *
 * @see Reference: sage/matrix/matrix_integer_dense.pyx
 * @see Reference: sage/matrix/symplectic_basis.py
 */

import { describe, expect, it } from 'vitest';
import { Rational } from '../rings/rational.js';
import {
  type IntegerMatrix,
  IntegerMatrixFromEntries,
  LLL,
  elementary_divisors_integer,
  frobenius_form_integer,
  hermite_normal_form,
  identity_integer_matrix,
  is_LLL_reduced,
  rank_integer,
  rational_reconstruction,
  smith_form_integer,
  symplectic_form_integer,
} from './matrix_integer.js';

function rowsOf(A: IntegerMatrix): bigint[][] {
  const rows: bigint[][] = [];
  for (let i = 0; i < A.nrows; i++) {
    const row: bigint[] = [];
    for (let j = 0; j < A.ncols; j++) row.push(A.get(i, j).value);
    rows.push(row);
  }
  return rows;
}

function eq(A: IntegerMatrix, B: IntegerMatrix): boolean {
  if (A.nrows !== B.nrows || A.ncols !== B.ncols) return false;
  for (let i = 0; i < A.nrows; i++) {
    for (let j = 0; j < A.ncols; j++) {
      if (A.get(i, j).value !== B.get(i, j).value) return false;
    }
  }
  return true;
}

/** Exact determinant over ZZ (fraction-free Bareiss). */
function det(M: bigint[][]): bigint {
  const n = M.length;
  if (n === 0) return 1n;
  const A = M.map((r) => [...r]);
  let sign = 1n;
  let prev = 1n;
  for (let k = 0; k < n - 1; k++) {
    if (A[k]![k] === 0n) {
      let sw = -1;
      for (let i = k + 1; i < n; i++) {
        if (A[i]![k] !== 0n) {
          sw = i;
          break;
        }
      }
      if (sw === -1) return 0n;
      [A[k], A[sw]] = [A[sw]!, A[k]!];
      sign = -sign;
    }
    for (let i = k + 1; i < n; i++) {
      for (let j = k + 1; j < n; j++) {
        A[i]![j] = (A[i]![j]! * A[k]![k]! - A[i]![k]! * A[k]![j]!) / prev;
      }
      A[i]![k] = 0n;
    }
    prev = A[k]![k]!;
  }
  return sign * A[n - 1]![n - 1]!;
}

function makeRandom(seed: number): (m: number) => number {
  let state = seed;
  return (m: number) => {
    state = (state * 1103515245 + 12345) & 0x7fffffff;
    return state % m;
  };
}

describe('smith_form_integer', () => {
  it('reproduces the Sage doctest matrix(ZZ,3,[1..9])', () => {
    const M = IntegerMatrixFromEntries([
      [1, 2, 3],
      [4, 5, 6],
      [7, 8, 9],
    ]);
    const [D, U, V] = smith_form_integer(M, true) as [IntegerMatrix, IntegerMatrix, IntegerMatrix];
    expect(rowsOf(D)).toEqual([
      [1n, 0n, 0n],
      [0n, 3n, 0n],
      [0n, 0n, 0n],
    ]);
    expect(eq(D, U.mul(M).mul(V))).toBe(true);
  });

  it('terminates on matrices that used to make the elimination loop cycle', () => {
    // Both of these hung forever before the pivot-reduction fix.
    const A = IntegerMatrixFromEntries([
      [-1, -1, -1],
      [-1, -1, -1],
      [0, -1, -1],
    ]);
    expect(rowsOf(smith_form_integer(A) as IntegerMatrix)).toEqual([
      [1n, 0n, 0n],
      [0n, 1n, 0n],
      [0n, 0n, 0n],
    ]);
    expect(rank_integer(A)).toBe(2);

    const B = IntegerMatrixFromEntries([
      [4, 6, -5, -4],
      [3, -1, 1, 0],
      [-4, 3, 5, -1],
    ]);
    expect(rowsOf(smith_form_integer(B) as IntegerMatrix)).toEqual([
      [1n, 0n, 0n, 0n],
      [0n, 1n, 0n, 0n],
      [0n, 0n, 1n, 0n],
    ]);
  });

  it('returns a diagonal D = U*M*V with a divisibility chain (random sweep)', () => {
    const rnd = makeRandom(20260727);
    for (let trial = 0; trial < 150; trial++) {
      const m = 1 + rnd(5);
      const n = 1 + rnd(5);
      const entries: number[][] = [];
      for (let i = 0; i < m; i++) {
        entries.push([]);
        for (let j = 0; j < n; j++) entries[i]!.push(rnd(13) - 6);
      }
      const M = IntegerMatrixFromEntries(entries);
      const [D, U, V] = smith_form_integer(M, true) as [
        IntegerMatrix,
        IntegerMatrix,
        IntegerMatrix,
      ];

      expect(eq(D, U.mul(M).mul(V))).toBe(true);
      expect(det(rowsOf(U)) ** 2n).toBe(1n);
      expect(det(rowsOf(V)) ** 2n).toBe(1n);

      const d = rowsOf(D);
      for (let i = 0; i < m; i++) {
        for (let j = 0; j < n; j++) {
          if (i !== j) expect(d[i]![j]).toBe(0n);
        }
      }
      for (let i = 0; i + 1 < Math.min(m, n); i++) {
        expect(d[i]![i]! >= 0n).toBe(true);
        if (d[i]![i] === 0n) {
          expect(d[i + 1]![i + 1]).toBe(0n);
        } else {
          expect(d[i + 1]![i + 1]! % d[i]![i]!).toBe(0n);
        }
      }
    }
  });
});

describe('elementary_divisors_integer', () => {
  it('returns one divisor per row, as PARI matsnf does', () => {
    // sage: M = Matrix([[3,0,1],[0,1,0]])
    // sage: M.elementary_divisors()          -> [1, 1]
    // sage: M.transpose().elementary_divisors() -> [1, 1, 0]
    const M = IntegerMatrixFromEntries([
      [3, 0, 1],
      [0, 1, 0],
    ]);
    expect(elementary_divisors_integer(M).map((d) => d.value)).toEqual([1n, 1n]);
    expect(elementary_divisors_integer(M.transpose()).map((d) => d.value)).toEqual([1n, 1n, 0n]);
  });

  it('reproduces Sage doctests for square matrices', () => {
    // sage: matrix(3, range(9)).elementary_divisors() -> [1, 3, 0]
    const A = IntegerMatrixFromEntries([
      [0, 1, 2],
      [3, 4, 5],
      [6, 7, 8],
    ]);
    expect(elementary_divisors_integer(A).map((d) => d.value)).toEqual([1n, 3n, 0n]);

    // sage: matrix(ZZ, 3, [1,5,7, 3,6,9, 0,1,2]).elementary_divisors() -> [1, 1, 6]
    const B = IntegerMatrixFromEntries([
      [1, 5, 7],
      [3, 6, 9],
      [0, 1, 2],
    ]);
    expect(elementary_divisors_integer(B).map((d) => d.value)).toEqual([1n, 1n, 6n]);

    // sage: MatrixSpace(ZZ,4)([3,4,5,6,7,3,8,10,14,5,6,7,2,2,10,9]).elementary_divisors()
    //   -> [1, 1, 1, 687]
    const C = IntegerMatrixFromEntries([
      [3, 4, 5, 6],
      [7, 3, 8, 10],
      [14, 5, 6, 7],
      [2, 2, 10, 9],
    ]);
    expect(elementary_divisors_integer(C).map((d) => d.value)).toEqual([1n, 1n, 1n, 687n]);
  });
});

describe('hermite_normal_form', () => {
  it('truncates the transformation matrix with include_zero_rows=false', () => {
    // sage: A = matrix(ZZ,5,3,[1..15])
    // sage: H, U = A.hermite_form(transformation=True, include_zero_rows=False)
    // sage: H
    // [1 2 3]
    // [0 3 6]
    // sage: U*A == H
    // True
    const A = IntegerMatrixFromEntries([
      [1, 2, 3],
      [4, 5, 6],
      [7, 8, 9],
      [10, 11, 12],
      [13, 14, 15],
    ]);
    const [H, U] = hermite_normal_form(A, undefined, undefined, false, true) as [
      IntegerMatrix,
      IntegerMatrix,
    ];
    expect(rowsOf(H)).toEqual([
      [1n, 2n, 3n],
      [0n, 3n, 6n],
    ]);
    expect(U.nrows).toBe(2);
    expect(U.ncols).toBe(5);
    expect(eq(U.mul(A), H)).toBe(true);
  });

  it('keeps U square when the zero rows are kept', () => {
    const A = IntegerMatrixFromEntries([
      [2, 0],
      [0, 2],
    ]);
    const [H, U] = hermite_normal_form(A, undefined, undefined, true, true) as [
      IntegerMatrix,
      IntegerMatrix,
    ];
    expect(rowsOf(H)).toEqual([
      [2n, 0n],
      [0n, 2n],
    ]);
    expect(eq(U.mul(A), H)).toBe(true);
  });
});

describe('LLL', () => {
  it('preserves the lattice for entries far beyond 2^53', () => {
    const big = 2n ** 60n;
    const M = IntegerMatrixFromEntries([
      [big + 1n, 1n, 0n],
      [0n, big + 3n, 1n],
      [1n, 0n, big + 7n],
    ]);
    const R = LLL(M) as IntegerMatrix;
    expect(det(rowsOf(R)) ** 2n).toBe(det(rowsOf(M)) ** 2n);
    expect(is_LLL_reduced(R, 0.99, 0.501)).toBe(true);
  });

  it('reproduces the Magma-handbook extended-gcd doctest', () => {
    // sage: Q = [67015143, 248934363018, ...]; X[i,0] = 100*Q[i]; X[i,i+1] = 1
    // sage: L = X.LLL(); L.row(n-1).list()[1:]
    // [-3, -1, 13, -1, -4, 2, 3, 4, 5, -1]
    const Q = [
      67015143n,
      248934363018n,
      109210n,
      25590011055n,
      74631449n,
      10230248n,
      709487n,
      68965012139n,
      972065n,
      864972271n,
    ];
    const n = Q.length;
    const X: bigint[][] = [];
    for (let i = 0; i < n; i++) {
      const row: bigint[] = new Array(n + 1).fill(0n);
      row[0] = 100n * Q[i]!;
      row[i + 1] = 1n;
      X.push(row);
    }
    const L = LLL(IntegerMatrixFromEntries(X)) as IntegerMatrix;
    const M = rowsOf(L)[n - 1]!.slice(1);
    expect(M).toEqual([-3n, -1n, 13n, -1n, -4n, 2n, 3n, 4n, 5n, -1n]);
    let s = 0n;
    for (let i = 0; i < n; i++) s += Q[i]! * M[i]!;
    expect(s).toBe(-1n);
  });

  it('puts the zero rows first for rank-deficient input', () => {
    // sage: matrix(4,3,[1,2,3,2,4,6,7,0,1,-1,-2,-3]).LLL()[0:2] is the zero matrix
    const M = IntegerMatrixFromEntries([
      [1, 2, 3],
      [2, 4, 6],
      [7, 0, 1],
      [-1, -2, -3],
    ]);
    const R = rowsOf(LLL(M) as IntegerMatrix);
    expect(R[0]).toEqual([0n, 0n, 0n]);
    expect(R[1]).toEqual([0n, 0n, 0n]);
    expect(R[2]!.some((x) => x !== 0n)).toBe(true);
  });

  it('rejects out-of-range delta the way Sage does', () => {
    const M = identity_integer_matrix(2);
    expect(() => LLL(M, 2)).toThrow('delta must be <= 1');
    expect(() => LLL(M, 0.25)).toThrow('delta must be > 0.25');
  });

  it('returns a unimodular transformation with U*A == R (random sweep)', () => {
    const rnd = makeRandom(424242);
    for (let trial = 0; trial < 40; trial++) {
      const n = 2 + rnd(3);
      const scale = [1n, 1n << 20n, 1n << 40n, 1n << 60n][rnd(4)]!;
      const entries: bigint[][] = [];
      for (let i = 0; i < n; i++) {
        entries.push([]);
        for (let j = 0; j < n; j++) {
          entries[i]!.push(BigInt(rnd(1000) - 500) * scale + BigInt(rnd(7)));
        }
      }
      const M = IntegerMatrixFromEntries(entries);
      const [R, U] = LLL(
        M,
        undefined,
        undefined,
        undefined,
        undefined,
        undefined,
        undefined,
        undefined,
        undefined,
        true
      ) as [IntegerMatrix, IntegerMatrix];

      expect(eq(U.mul(M), R)).toBe(true);
      expect(det(rowsOf(U)) ** 2n).toBe(1n);
      const d0 = det(entries);
      if (d0 !== 0n) {
        expect(det(rowsOf(R)) ** 2n).toBe(d0 ** 2n);
        expect(is_LLL_reduced(R, 0.99, 0.501)).toBe(true);
      }
    }
  });
});

describe('frobenius_form_integer', () => {
  /** Parse a Sage-style rational literal ("−23/15", "2") into a Rational. */
  const rat = (s: string): Rational => {
    const [n, d] = s.split('/');
    return new Rational(BigInt(n!), d === undefined ? 1n : BigInt(d));
  };
  const ratRows = (M: Rational[][]): string[][] =>
    M.map((row) => row.map((x) => (x.denominator === 1n ? `${x.numerator}` : `${x}`)));

  const ratMul = (A: Rational[][], B: Rational[][]): Rational[][] => {
    const n = A.length;
    return Array.from({ length: n }, (_, i) =>
      Array.from({ length: n }, (_, j) => {
        let s = new Rational(0n);
        for (let k = 0; k < n; k++) s = s.add(A[i]![k]!.mul(B[k]![j]!));
        return s;
      })
    );
  };
  /** Exact Gauss-Jordan inverse over QQ; throws if the matrix is singular. */
  const ratInv = (A: Rational[][]): Rational[][] => {
    const n = A.length;
    const M = A.map((row, i) => [
      ...row,
      ...Array.from({ length: n }, (_, j) => new Rational(i === j ? 1n : 0n)),
    ]);
    for (let c = 0; c < n; c++) {
      let p = -1;
      for (let i = c; i < n; i++) {
        if (!M[i]![c]!.isZero()) {
          p = i;
          break;
        }
      }
      if (p < 0) throw new Error('change-of-basis matrix is singular');
      if (p !== c) {
        const t = M[p]!;
        M[p] = M[c]!;
        M[c] = t;
      }
      const pinv = M[c]![c]!.inv();
      for (let j = 0; j < 2 * n; j++) M[c]![j] = M[c]![j]!.mul(pinv);
      for (let i = 0; i < n; i++) {
        if (i === c || M[i]![c]!.isZero()) continue;
        const f = M[i]![c]!;
        for (let j = 0; j < 2 * n; j++) M[i]![j] = M[i]![j]!.sub(f.mul(M[c]![j]!));
      }
    }
    return M.map((row) => row.slice(n));
  };

  /**
   * Golden values produced by SageMath 10.3 (i.e. by PARI's `matfrobenius`):
   * `[name, A, F, elementary divisors, B]`.
   */
  const golden: [string, number[][], string[][], string[][], string[][]][] = [
    [
      'MatrixSpace(ZZ,3)(range(9))',
      [
        [0, 1, 2],
        [3, 4, 5],
        [6, 7, 8],
      ],
      [
        ['0', '0', '0'],
        ['1', '0', '18'],
        ['0', '1', '12'],
      ],
      [['0', '-18', '-12', '1']],
      [
        ['1', '-2', '1'],
        ['0', '-23/15', '14/15'],
        ['0', '2/15', '-1/15'],
      ],
    ],
    [
      'diag(2,2,3)',
      [
        [2, 0, 0],
        [0, 2, 0],
        [0, 0, 3],
      ],
      [
        ['0', '-6', '0'],
        ['1', '5', '0'],
        ['0', '0', '2'],
      ],
      [
        ['6', '-5', '1'],
        ['-2', '1'],
      ],
      [
        ['-3/2', '-9/2', '-2'],
        ['1/2', '3/2', '1'],
        ['1/2', '1/2', '0'],
      ],
    ],
    [
      'diag(1,1,2)',
      [
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 2],
      ],
      [
        ['0', '-2', '0'],
        ['1', '3', '0'],
        ['0', '0', '1'],
      ],
      [
        ['2', '-3', '1'],
        ['-1', '1'],
      ],
      [
        ['-2', '-4', '-1'],
        ['1', '2', '1'],
        ['1', '1', '0'],
      ],
    ],
    [
      'identity 2x2',
      [
        [1, 0],
        [0, 1],
      ],
      [
        ['1', '0'],
        ['0', '1'],
      ],
      [
        ['-1', '1'],
        ['-1', '1'],
      ],
      [
        ['1', '0'],
        ['0', '1'],
      ],
    ],
    [
      '[[0,1],[2,3]]',
      [
        [0, 1],
        [2, 3],
      ],
      [
        ['0', '2'],
        ['1', '3'],
      ],
      [['-2', '-3', '1']],
      [
        ['1', '0'],
        ['0', '1/2'],
      ],
    ],
    [
      'diag(2,2)',
      [
        [2, 0],
        [0, 2],
      ],
      [
        ['2', '0'],
        ['0', '2'],
      ],
      [
        ['-2', '1'],
        ['-2', '1'],
      ],
      [
        ['1', '0'],
        ['0', '1'],
      ],
    ],
    [
      'a Jordan block plus a scalar',
      [
        [1, 1, 0],
        [0, 1, 0],
        [0, 0, 1],
      ],
      [
        ['0', '-1', '0'],
        ['1', '2', '0'],
        ['0', '0', '1'],
      ],
      [
        ['1', '-2', '1'],
        ['-1', '1'],
      ],
      [
        ['-1', '1', '0'],
        ['1', '0', '0'],
        ['0', '0', '1'],
      ],
    ],
    [
      'MatrixSpace(ZZ,4)(range(16))',
      [
        [0, 1, 2, 3],
        [4, 5, 6, 7],
        [8, 9, 10, 11],
        [12, 13, 14, 15],
      ],
      [
        ['0', '0', '0', '0'],
        ['1', '0', '80', '0'],
        ['0', '1', '30', '0'],
        ['0', '0', '0', '0'],
      ],
      [
        ['0', '-80', '-30', '1'],
        ['0', '1'],
      ],
      [
        ['1', '-16/7', '11/7', '-2/7'],
        ['0', '-59/56', '4/7', '3/56'],
        ['0', '1/28', '-1/56', '0'],
        ['0', '1', '-2', '1'],
      ],
    ],
    [
      'a nilpotent Jordan block',
      [
        [0, 1, 0],
        [0, 0, 1],
        [0, 0, 0],
      ],
      [
        ['0', '0', '0'],
        ['1', '0', '0'],
        ['0', '1', '0'],
      ],
      [['0', '0', '0', '1']],
      [
        ['0', '0', '1'],
        ['0', '1', '0'],
        ['1', '0', '0'],
      ],
    ],
    [
      'three blocks',
      [
        [2, 0, 0, 0],
        [0, 2, 0, 0],
        [0, 0, 2, 1],
        [0, 0, 0, 2],
      ],
      [
        ['0', '-4', '0', '0'],
        ['1', '4', '0', '0'],
        ['0', '0', '2', '0'],
        ['0', '0', '0', '2'],
      ],
      [
        ['4', '-4', '1'],
        ['-2', '1'],
        ['-2', '1'],
      ],
      [
        ['-1', '1', '-2', '4'],
        ['1/2', '-1/2', '1', '-3/2'],
        ['1/2', '0', '0', '-1/2'],
        ['0', '1/2', '0', '1'],
      ],
    ],
    [
      'a 5x5 with a degree-4 minimal polynomial',
      [
        [1, 0, 0, 0, 0],
        [0, 1, 0, 0, 0],
        [0, 0, 2, 3, 0],
        [0, 0, 0, 2, 0],
        [0, 0, 0, 0, 5],
      ],
      [
        ['0', '0', '0', '-20', '0'],
        ['1', '0', '0', '44', '0'],
        ['0', '1', '0', '-33', '0'],
        ['0', '0', '1', '10', '0'],
        ['0', '0', '0', '0', '1'],
      ],
      [
        ['20', '-44', '33', '-10', '1'],
        ['-1', '1'],
      ],
      [
        ['-5/9', '2/9', '10/27', '-35/27', '-1/9'],
        ['2/3', '-4/15', '-17/27', '52/27', '2/9'],
        ['-1/4', '1/10', '8/27', '-19/27', '-5/36'],
        ['1/36', '-1/90', '-1/27', '2/27', '1/36'],
        ['2', '-17/10', '0', '0', '0'],
      ],
    ],
  ];

  it.each(golden)('matches SageMath for %s', (_name, entries, F, polys, B) => {
    const A = IntegerMatrixFromEntries(entries);

    // flag = 0: the Frobenius form itself
    expect(rowsOf(frobenius_form_integer(A, 0) as IntegerMatrix)).toEqual(
      F.map((row) => row.map((x) => BigInt(x)))
    );

    // flag = 1: the elementary divisor polynomials, constant term first.
    // PARI orders them so that P_{i+1} divides P_i (minimal polynomial first).
    expect(frobenius_form_integer(A, 1)).toEqual(polys.map((p) => p.map((c) => BigInt(c))));

    // flag = 2: [F, B] over QQ
    const [F2, B2] = frobenius_form_integer(A, 2) as [Rational[][], Rational[][]];
    expect(ratRows(F2)).toEqual(F);
    expect(ratRows(B2)).toEqual(B);

    // sage: A == B^(-1)*F*B  -- the defining identity, checked exactly over QQ.
    const reconstructed = ratMul(ratMul(ratInv(B2), F2), B2);
    for (let i = 0; i < entries.length; i++) {
      for (let j = 0; j < entries.length; j++) {
        expect(reconstructed[i]![j]!.eq(BigInt(entries[i]![j]!))).toBe(true);
      }
    }
  });

  it('produces an invertible B with B^-1*F*B == A on random matrices', () => {
    let state = 20260728;
    const rnd = (m: number): number => {
      state = (state * 1103515245 + 12345) & 0x7fffffff;
      return state % m;
    };
    for (let trial = 0; trial < 150; trial++) {
      const n = 1 + rnd(5);
      const entries: number[][] = [];
      for (let i = 0; i < n; i++) {
        entries.push([]);
        for (let j = 0; j < n; j++) entries[i]!.push(trial % 3 === 0 ? rnd(3) - 1 : rnd(21) - 10);
      }
      const A = IntegerMatrixFromEntries(entries);
      const [F, B] = frobenius_form_integer(A, 2) as [Rational[][], Rational[][]];
      // ratInv throws on a singular matrix, so this also asserts B is invertible.
      const reconstructed = ratMul(ratMul(ratInv(B), F), B);
      for (let i = 0; i < n; i++) {
        for (let j = 0; j < n; j++) {
          expect(reconstructed[i]![j]!.eq(BigInt(entries[i]![j]!))).toBe(true);
        }
      }
      // F is the same matrix flag=0 returns, and it is integral.
      expect(rowsOf(frobenius_form_integer(A, 0) as IntegerMatrix)).toEqual(
        F.map((row) => row.map((x) => x.numerator))
      );
      for (const row of F) for (const x of row) expect(x.denominator).toBe(1n);
      // The elementary divisors have total degree n and form a divisibility
      // chain P_{i+1} | P_i (PARI's ordering).
      const polys = frobenius_form_integer(A, 1) as bigint[][];
      expect(polys.reduce((s, p) => s + p.length - 1, 0)).toBe(n);
    }
  });

  it('handles the empty matrix (sage: matrix([]).frobenius_form(2) == ([], []))', () => {
    const A = IntegerMatrixFromEntries([]);
    expect(rowsOf(frobenius_form_integer(A, 0) as IntegerMatrix)).toEqual([]);
    expect(frobenius_form_integer(A, 1)).toEqual([]);
    expect(frobenius_form_integer(A, 2)).toEqual([[], []]);
  });

  it('rejects non-square matrices', () => {
    const A = IntegerMatrixFromEntries([
      [1, 2, 3],
      [4, 5, 6],
    ]);
    expect(() => frobenius_form_integer(A)).toThrow(
      'frobenius matrix of non-square matrix not defined.'
    );
  });

  it('rejects an out-of-range flag (PARI pari_err_FLAG)', () => {
    const A = identity_integer_matrix(2);
    expect(() => frobenius_form_integer(A, 3)).toThrow('incorrect flag in matfrobenius');
  });
});

describe('symplectic_form_integer', () => {
  it('reproduces the Sage 5x5 doctest', () => {
    // sage: E = matrix(ZZ, 5, 5, [0,14,0,-8,-2, -14,0,-3,-11,4, 0,3,0,0,0,
    //                             8,11,0,0,8, 2,-4,0,-8,0])
    // sage: F, C = E.symplectic_form(); F
    // [ 0  0  1  0  0]
    // [ 0  0  0  2  0]
    // [-1  0  0  0  0]
    // [ 0 -2  0  0  0]
    // [ 0  0  0  0  0]
    const E = IntegerMatrixFromEntries([
      [0, 14, 0, -8, -2],
      [-14, 0, -3, -11, 4],
      [0, 3, 0, 0, 0],
      [8, 11, 0, 0, 8],
      [2, -4, 0, -8, 0],
    ]);
    const [F, C] = symplectic_form_integer(E);
    expect(rowsOf(F)).toEqual([
      [0n, 0n, 1n, 0n, 0n],
      [0n, 0n, 0n, 2n, 0n],
      [-1n, 0n, 0n, 0n, 0n],
      [0n, -2n, 0n, 0n, 0n],
      [0n, 0n, 0n, 0n, 0n],
    ]);
    expect(eq(F, C.mul(E).mul(C.transpose()))).toBe(true);
  });

  it('reproduces the 8x8 doctest from symplectic_basis.py', () => {
    const E = IntegerMatrixFromEntries([
      [0, 0, 0, 0, 1, 0, 0, 0],
      [0, 0, 0, 0, 0, 1, 0, 0],
      [0, 0, 0, 0, 0, 0, 1, 0],
      [0, 0, 0, 0, 0, 0, 0, 20191],
      [-1, 0, 0, 0, 0, 0, 0, 0],
      [0, -1, 0, 0, 0, 0, 0, 0],
      [0, 0, -1, 0, 0, 0, 0, 0],
      [0, 0, 0, -20191, 0, 0, 0, 0],
    ]);
    const [F, C] = symplectic_form_integer(E);
    expect(rowsOf(F)).toEqual(rowsOf(E));
    expect(eq(F, C.mul(E).mul(C.transpose()))).toBe(true);
  });

  it('satisfies F = C*M*C^T with a [[0,D],[-D,0]] layout (random sweep)', () => {
    const rnd = makeRandom(31337);
    for (let trial = 0; trial < 100; trial++) {
      const n = 1 + rnd(6);
      const e: bigint[][] = [];
      for (let i = 0; i < n; i++) e.push(new Array(n).fill(0n));
      for (let i = 0; i < n; i++) {
        for (let j = i + 1; j < n; j++) {
          const v = BigInt(rnd(21) - 10);
          e[i]![j] = v;
          e[j]![i] = -v;
        }
      }
      const M = IntegerMatrixFromEntries(e);
      const [F, C] = symplectic_form_integer(M);
      expect(eq(F, C.mul(M).mul(C.transpose()))).toBe(true);

      const a = rowsOf(F);
      let r = 0;
      for (const row of a) if (row.some((x) => x !== 0n)) r++;
      expect(r % 2).toBe(0);
      const half = r / 2;
      for (let i = 0; i < n; i++) {
        for (let j = 0; j < n; j++) {
          const inBlock = (i < half && j === i + half) || (j < half && i === j + half);
          if (!inBlock) expect(a[i]![j]).toBe(0n);
        }
      }
      for (let i = 0; i < half; i++) {
        expect(a[i]![i + half]! > 0n).toBe(true);
        if (i + 1 < half) {
          expect(a[i + 1]![i + 1 + half]! % a[i]![i + half]!).toBe(0n);
        }
      }
    }
  });

  it('rejects matrices that are not anti-symmetric or not alternating', () => {
    expect(() =>
      symplectic_form_integer(
        IntegerMatrixFromEntries([
          [0, 1],
          [1, 0],
        ])
      )
    ).toThrow('Can only find symplectic bases for anti-symmetric matrices');

    expect(() =>
      symplectic_form_integer(
        IntegerMatrixFromEntries([
          [1, 2, 3],
          [4, 5, 6],
        ])
      )
    ).toThrow('Can only find symplectic bases for square matrices');
  });
});

describe('rational_reconstruction', () => {
  it('raises ZeroDivisionError for a zero modulus (Sage issue #9345)', () => {
    const A = IntegerMatrixFromEntries([[1]]);
    expect(() => rational_reconstruction(A, 0n)).toThrow('The modulus cannot be zero');
  });

  it('raises ValueError when no reconstruction exists', () => {
    // sage: matrix(ZZ, 4, [4,-4,7,1, -1,1,-1,-12, -1,-1,1,-1, -3,1,5,-1])
    //         .rational_reconstruction(11)
    // ValueError: rational reconstruction does not exist
    const A = IntegerMatrixFromEntries([
      [4, -4, 7, 1],
      [-1, 1, -1, -12],
      [-1, -1, 1, -1],
      [-3, 1, 5, -1],
    ]);
    expect(() => rational_reconstruction(A, 11n)).toThrow('rational reconstruction does not exist');
  });

  it('reconstructs p/q with p/q = a mod N', () => {
    // 1/2 mod 101 == 51
    const A = IntegerMatrixFromEntries([[51]]);
    const { numerators, denominators } = rational_reconstruction(A, 101n);
    const p = numerators.get(0, 0).value;
    const q = denominators.get(0, 0).value;
    expect(q).not.toBe(0n);
    // p/q == 51 (mod 101), i.e. p == 51*q (mod 101)
    expect((((p - 51n * q) % 101n) + 101n) % 101n).toBe(0n);
  });
});
