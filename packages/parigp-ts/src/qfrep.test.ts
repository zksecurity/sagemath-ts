/**
 * Tests for {@link qfrep0} (PARI `qfrep`).
 *
 * Oracles, in decreasing order of authority:
 *
 *  1. PARI's own regression test suite, copied verbatim from
 *     `reference/pari/src/test/in/qf` + `reference/pari/src/test/32/qf`
 *     (the 12-dimensional form).
 *  2. cypari2's doctests, `reference/cypari2/cypari2/gen.pyx:4222-4251`.
 *  3. Exhaustive enumeration of every lattice vector in a provably sufficient
 *     box (see {@link bruteforce}).
 *  4. Closed-form theta series (E8, D4, A2, sums of squares).
 */

import { describe, expect, test } from 'bun:test';
import {
  PariDomainError,
  PariPrecError,
  PariTypeError,
  type ZM,
  ZM_det,
  isqrt,
  lllgramint,
  qf_ZM_apply,
  qfrep,
  qfrep0,
} from './qfrep.js';

/* ------------------------------------------------------------------ */
/* helpers                                                             */
/* ------------------------------------------------------------------ */

/** Build a PARI column-major matrix from rows. */
function fromRows(rows: (number | bigint)[][]): ZM {
  const m = rows.length;
  const n = m === 0 ? 0 : rows[0]!.length;
  const A: ZM = [];
  for (let j = 0; j < n; j++) {
    const col: bigint[] = [];
    for (let i = 0; i < m; i++) col.push(BigInt(rows[i]![j]!));
    A.push(col);
  }
  return A;
}

function nums(v: bigint[]): number[] {
  return v.map(Number);
}

/**
 * Exhaustive count of `q(v) = i`, `1 <= i <= B`, over ALL integer vectors.
 *
 * The search box is provably sufficient: for any `v`, Cauchy-Schwarz on
 * `v_i = e_i^T v` gives `v_i^2 <= (e_i^T a^{-1} e_i) * (v^T a v)`, i.e.
 * `v_i^2 * det(a) <= B * det(a with row i and column i removed)`.
 *
 * Returns half-counts, like `qfrep`.
 */
function bruteforce(rows: number[][], B: number): number[] {
  const n = rows.length;
  const A = fromRows(rows);
  const det = ZM_det(A);
  const lim: number[] = [];
  for (let i = 0; i < n; i++) {
    const minor: ZM = [];
    for (let j = 0; j < n; j++) {
      if (j === i) continue;
      const col: bigint[] = [];
      for (let r = 0; r < n; r++) {
        if (r === i) continue;
        col.push(A[j]![r]!);
      }
      minor.push(col);
    }
    const num = BigInt(B) * ZM_det(minor);
    lim.push(Number(isqrt(num / det)));
  }
  const out = new Array<number>(B).fill(0);
  const x = new Array<number>(n).fill(0);
  const rec = (i: number) => {
    if (i === n) {
      let s = 0;
      for (let r = 0; r < n; r++) {
        for (let c = 0; c < n; c++) s += rows[r]![c]! * x[r]! * x[c]!;
      }
      if (s >= 1 && s <= B) out[s - 1]! += 1;
      return;
    }
    for (let t = -lim[i]!; t <= lim[i]!; t++) {
      x[i] = t;
      rec(i + 1);
    }
  };
  rec(0);
  /* qfrep counts half of the +-pairs */
  return out.map((v) => {
    expect(v % 2).toBe(0);
    return v / 2;
  });
}

/* Standard lattices. */
const I = (n: number): number[][] =>
  Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => (i === j ? 1 : 0)));

/** A2: hexagonal lattice, minimal norm 2. */
const A2 = [
  [2, -1],
  [-1, 2],
];

/** D4: Gram matrix of the checkerboard lattice, minimal norm 2. */
const D4 = [
  [2, -1, 0, 0],
  [-1, 2, -1, -1],
  [0, -1, 2, 0],
  [0, -1, 0, 2],
];

/** E8: Gram matrix (Cartan matrix of E8), unimodular, even, minimal norm 2. */
const E8 = [
  [2, -1, 0, 0, 0, 0, 0, 0],
  [-1, 2, -1, 0, 0, 0, 0, 0],
  [0, -1, 2, -1, 0, 0, 0, 0],
  [0, 0, -1, 2, -1, 0, 0, 0],
  [0, 0, 0, -1, 2, -1, 0, -1],
  [0, 0, 0, 0, -1, 2, -1, 0],
  [0, 0, 0, 0, 0, -1, 2, 0],
  [0, 0, 0, 0, -1, 0, 0, 2],
];

function sigma(k: number, m: number): number {
  let s = 0;
  for (let d = 1; d <= m; d++) if (m % d === 0) s += d ** k;
  return s;
}
/** sum of the odd divisors of m */
function sigmaOdd(m: number): number {
  let s = 0;
  for (let d = 1; d <= m; d++) if (m % d === 0 && d % 2 === 1) s += d;
  return s;
}

/* ------------------------------------------------------------------ */

describe('isqrt', () => {
  test('exact integer square root', () => {
    for (let n = 0n; n < 200n; n++) {
      const r = isqrt(n);
      expect(r * r <= n).toBe(true);
      expect((r + 1n) * (r + 1n) > n).toBe(true);
    }
    expect(isqrt(10n ** 40n)).toBe(10n ** 20n);
    expect(isqrt(10n ** 40n - 1n)).toBe(10n ** 20n - 1n);
  });
});

describe('qfrep: PARI documentation example (functions/linear_algebra/qfrep)', () => {
  const q = fromRows([
    [2, 1],
    [1, 3],
  ]);

  test('qfrep(q, 5) = Vecsmall([0, 1, 2, 0, 0])', () => {
    expect(nums(qfrep0(q, 5))).toEqual([0, 1, 2, 0, 0]);
  });

  test('qfrep(q, 5, 1): the GP doc example is stale in its last entry', () => {
    /*
     * `reference/pari/src/functions/linear_algebra/qfrep` prints
     *     ? qfrep(q, 5, 1)
     *     %3 = Vecsmall([1, 0, 0, 1, 0])
     * but that last 0 is wrong for the code that is actually vendored here:
     * q(-1, 2) = 2 - 4 + 12 = 10 = 2*5, so entry 5 must be 1.  Verified three
     * ways: (a) brute force below, (b) `minim0_dolll` (bibli1.c:1324-1330)
     * allocates L of length B and *then* doubles the search bound to 2B, so
     * norm 2B is in range, (c) PARI's own regression test
     * `test/32/qf` has qfrep(Q,8,1)[8] = 4268 = qfrep(Q,16)[16] != 0,
     * i.e. the top entry is genuinely counted.
     */
    expect(nums(qfrep0(q, 5, 1))).toEqual([1, 0, 0, 1, 1]);

    const bf = bruteforce(
      [
        [2, 1],
        [1, 3],
      ],
      10
    );
    expect(nums(qfrep0(q, 5, 1))).toEqual([bf[1]!, bf[3]!, bf[5]!, bf[7]!, bf[9]!]);
  });
});

describe('qfrep: cypari2 doctests (cypari2/gen.pyx:4236-4241)', () => {
  const M = fromRows([
    [5, 1, 1],
    [1, 3, 1],
    [1, 1, 1],
  ]);

  test('M.qfrep(20)', () => {
    expect(nums(qfrep0(M, 20))).toEqual([
      1, 1, 2, 2, 2, 4, 4, 3, 3, 4, 2, 4, 6, 0, 4, 6, 4, 5, 6, 4,
    ]);
  });

  test('M.qfrep(20, flag=1)', () => {
    expect(nums(qfrep0(M, 20, 1))).toEqual([
      1, 2, 4, 3, 4, 4, 0, 6, 5, 4, 12, 4, 4, 8, 0, 3, 8, 6, 12, 12,
    ]);
  });

  test('flag=2 (cypari2: t_VEC vs t_VECSMALL) is the same data here', () => {
    expect(qfrep(M, 20, 2)).toEqual(qfrep0(M, 20, 0));
    expect(qfrep(M, 20, 3)).toEqual(qfrep0(M, 20, 1));
  });
});

describe('qfrep: PARI regression test (test/in/qf + test/32/qf)', () => {
  /* pari/src/test/in/qf:9 -- a 12-dimensional form */
  const Qrows = [
    [8, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4],
    [4, 8, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4],
    [4, 4, 8, 0, 0, 0, 3, 0, 0, 0, 0, 0],
    [4, 4, 0, 8, 4, 4, 1, 4, 4, 4, 4, 4],
    [4, 4, 0, 4, 8, 4, 4, 4, 4, 4, 4, 4],
    [4, 4, 0, 4, 4, 8, 4, 4, 4, 4, 4, 4],
    [4, 4, 3, 1, 4, 4, 8, 4, 1, 1, 1, 1],
    [4, 4, 0, 4, 4, 4, 4, 8, 4, 4, 4, 4],
    [4, 4, 0, 4, 4, 4, 1, 4, 8, 4, 4, 4],
    [4, 4, 0, 4, 4, 4, 1, 4, 4, 8, 4, 4],
    [4, 4, 0, 4, 4, 4, 1, 4, 4, 4, 8, 4],
    [4, 4, 0, 4, 4, 4, 1, 4, 4, 4, 4, 8],
  ];
  const Q = fromRows(Qrows);

  test('qfrep(Q, 16)', () => {
    expect(nums(qfrep0(Q, 16))).toEqual([
      0, 0, 0, 0, 0, 0, 0, 133, 0, 165, 0, 638, 0, 396, 0, 4268,
    ]);
  });

  test('qfrep(Q, 8, 1)', () => {
    expect(nums(qfrep0(Q, 8, 1))).toEqual([0, 0, 0, 133, 165, 638, 396, 4268]);
  });
});

describe('qfrep: exhaustive enumeration (identity forms, dimensions 2-5)', () => {
  for (const n of [1, 2, 3, 4, 5]) {
    test(`identity form of dimension ${n}`, () => {
      const B = n <= 3 ? 30 : 20;
      expect(nums(qfrep0(fromRows(I(n)), B))).toEqual(bruteforce(I(n), B));
    });
  }

  test('r_4(n) = 8 * (sum of divisors not divisible by 4)', () => {
    const B = 24;
    const got = nums(qfrep0(fromRows(I(4)), B));
    for (let n = 1; n <= B; n++) {
      let s = 0;
      for (let d = 1; d <= n; d++) if (n % d === 0 && d % 4 !== 0) s += d;
      expect(got[n - 1]).toBe(4 * s); /* half of 8*s */
    }
  });

  test('r_2(n) = 4 * (d_1(n) - d_3(n))', () => {
    const B = 40;
    const got = nums(qfrep0(fromRows(I(2)), B));
    for (let n = 1; n <= B; n++) {
      let s = 0;
      for (let d = 1; d <= n; d++) {
        if (n % d !== 0) continue;
        if (d % 4 === 1) s += 1;
        else if (d % 4 === 3) s -= 1;
      }
      expect(got[n - 1]).toBe(2 * s);
    }
  });
});

describe('qfrep: root lattices', () => {
  test('A2 against exhaustive enumeration and its theta series', () => {
    const B = 30;
    const got = nums(qfrep0(fromRows(A2), B));
    expect(got).toEqual(bruteforce(A2, B));
    /* This Gram matrix gives q(x,y) = 2(x^2 - xy + y^2), so all norms are even
     * and the number of vectors of norm 2m is the hexagonal theta coefficient
     * r(m) = 6 * sum_{d|m} chi_{-3}(d) (6 roots of norm 2, etc.). */
    for (let n = 1; n <= B; n++) {
      if (n % 2 === 1) {
        expect(got[n - 1]).toBe(0);
        continue;
      }
      const m = n / 2;
      let s = 0;
      for (let d = 1; d <= m; d++) {
        if (m % d !== 0) continue;
        if (d % 3 === 1) s += 1;
        else if (d % 3 === 2) s -= 1;
      }
      expect(got[n - 1]).toBe(3 * s);
    }
  });

  test('D4 against exhaustive enumeration and its theta series', () => {
    const B = 20;
    const got = nums(qfrep0(fromRows(D4), B));
    expect(got).toEqual(bruteforce(D4, B));
    /* r_{D4}(2m) = 24 * sigma_odd(m); odd norms do not occur */
    for (let n = 1; n <= B; n++) {
      if (n % 2 === 1) expect(got[n - 1]).toBe(0);
      else expect(got[n - 1]).toBe(12 * sigmaOdd(n / 2));
    }
  });

  test('E8 theta series: r(2m) = 240 * sigma_3(m)', () => {
    const B = 16;
    const got = nums(qfrep0(fromRows(E8), B));
    for (let n = 1; n <= B; n++) {
      if (n % 2 === 1) expect(got[n - 1]).toBe(0);
      else expect(got[n - 1]).toBe(120 * sigma(3, n / 2));
    }
    /* 240 vectors of norm 2 -- the E8 root system */
    expect(got[1]).toBe(120);
    /* flag = 1 view of the same data */
    expect(nums(qfrep0(fromRows(E8), B / 2, 1))).toEqual([
      got[1]!,
      got[3]!,
      got[5]!,
      got[7]!,
      got[9]!,
      got[11]!,
      got[13]!,
      got[15]!,
    ]);
  });

  test('E8 against exhaustive enumeration in the standard D8+ coordinates', () => {
    /*
     * Fully independent oracle: enumerate
     *     E8 = { x in Z^8 u (Z + 1/2)^8 : sum x_i in 2Z },  q(x) = sum x_i^2
     * directly, with no Gram matrix, no Cholesky and no LLL involved.  Work
     * with y = 2x so everything stays integral: q = (sum y_i^2) / 4.
     */
    const B = 16;
    const cnt = new Array<number>(4 * B + 1).fill(0);
    const rec = (i: number, s: number, sy: number, half: boolean) => {
      if (i === 8) {
        /* sum x_i = (sum y_i)/2 must be even, i.e. sum y_i = 0 mod 4 */
        if (sy % 4 === 0) cnt[s]! += 1;
        return;
      }
      const lim = Math.floor(Math.sqrt(4 * B - s));
      /* y_i all even (x integral) or all odd (x half-integral) */
      const first = half ? -(lim % 2 === 1 ? lim : lim - 1) : -(lim - (lim % 2));
      for (let v = first; v <= lim; v += 2) {
        if (s + v * v > 4 * B) continue;
        rec(i + 1, s + v * v, (((sy + v) % 4) + 4) % 4, half);
      }
    };
    rec(0, 0, 0, false);
    rec(0, 0, 0, true);
    const expected: number[] = [];
    for (let n = 1; n <= B; n++) expected.push(cnt[4 * n]! / 2);
    expect(nums(qfrep0(fromRows(E8), B))).toEqual(expected);
    /* sanity: the enumeration really found the 240 roots */
    expect(expected[1]).toBe(120);
  });
});

describe('qfrep: invariance under unimodular change of basis (exercises the LLL step)', () => {
  /* Deliberately skewed bases: without the lllgramint reduction of
   * bibli1.c:1345 the enumeration tree is astronomically larger. */
  const U3 = fromRows([
    [1, 17, -253],
    [0, 1, -14],
    [0, 0, 1],
  ]);
  const U4 = fromRows([
    [1, 0, 0, 0],
    [23, 1, 0, 0],
    [-7, 41, 1, 0],
    [113, -5, 29, 1],
  ]);

  test('identity form in dimension 3, skewed basis', () => {
    const G = qf_ZM_apply(fromRows(I(3)), U3);
    expect(nums(qfrep0(G, 25))).toEqual(nums(qfrep0(fromRows(I(3)), 25)));
  });

  test('D4 with a skewed basis', () => {
    const G = qf_ZM_apply(fromRows(D4), U4);
    expect(nums(qfrep0(G, 20))).toEqual(nums(qfrep0(fromRows(D4), 20)));
  });

  test('E8 with a skewed basis', () => {
    const U8: ZM = fromRows([
      [1, 3, 0, 0, 0, 0, 0, 0],
      [0, 1, 0, 0, 0, 0, 0, 0],
      [5, 0, 1, 0, 0, 0, 0, 0],
      [0, 0, 7, 1, 0, 0, 0, 0],
      [0, 0, 0, -2, 1, 0, 0, 0],
      [0, 11, 0, 0, 4, 1, 0, 0],
      [0, 0, 0, 0, 0, -6, 1, 0],
      [0, 0, 0, 0, 0, 0, 9, 1],
    ]);
    const G = qf_ZM_apply(fromRows(E8), U8);
    expect(nums(qfrep0(G, 12))).toEqual(nums(qfrep0(fromRows(E8), 12)));
  });
});

describe('lllgramint', () => {
  const forms: number[][][] = [
    I(2),
    I(5),
    A2,
    D4,
    E8,
    [
      [5, 1, 1],
      [1, 3, 1],
      [1, 1, 1],
    ],
  ];

  test('returns a unimodular transform reducing the Gram matrix', () => {
    for (const rows of forms) {
      const G = fromRows(rows);
      const u = lllgramint(G);
      expect(u).not.toBeNull();
      const det = ZM_det(u!);
      expect(det === 1n || det === -1n).toBe(true);
      const R = qf_ZM_apply(G, u!);
      /* reduced Gram matrix has the same determinant ... */
      expect(ZM_det(R)).toBe(ZM_det(G));
      /* ... and is symmetric with a positive diagonal */
      const n = R.length;
      for (let i = 0; i < n; i++) {
        expect(R[i]![i]! > 0n).toBe(true);
        for (let j = 0; j < n; j++) expect(R[j]![i]).toBe(R[i]![j]!);
      }
    }
  });

  test('skewed bases really are reduced', () => {
    const U = fromRows([
      [1, 17, -253],
      [0, 1, -14],
      [0, 0, 1],
    ]);
    const G = qf_ZM_apply(fromRows(I(3)), U);
    const u = lllgramint(G)!;
    const R = qf_ZM_apply(G, u);
    /* the identity lattice has minimum 1: a reduced basis has unit diagonal */
    for (let i = 0; i < 3; i++) expect(R[i]![i]).toBe(1n);
  });

  test('rejects forms that are not positive definite', () => {
    expect(
      lllgramint(
        fromRows([
          [1, 2],
          [2, 1],
        ])
      )
    ).toBeNull();
    expect(lllgramint(fromRows([[0]]))).toBeNull();
    expect(
      lllgramint(
        fromRows([
          [1, 1],
          [1, 1],
        ])
      )
    ).toBeNull();
  });
});

describe('qfrep: edge cases and errors', () => {
  test('B <= 0 gives the empty vector', () => {
    const q = fromRows([
      [2, 1],
      [1, 3],
    ]);
    expect(qfrep0(q, 0)).toEqual([]);
    expect(qfrep0(q, -5)).toEqual([]);
  });

  test('non-integral B is floored (gfloor)', () => {
    const q = fromRows([
      [2, 1],
      [1, 3],
    ]);
    expect(nums(qfrep0(q, 5.9))).toEqual([0, 1, 2, 0, 0]);
  });

  test('0-dimensional form gives all-zero counts', () => {
    expect(nums(qfrep0([], 3))).toEqual([0, 0, 0]);
  });

  test('1-dimensional form', () => {
    /* q(x) = 3x^2: norms 3, 12, 27, ... one +- pair each */
    expect(nums(qfrep0(fromRows([[3]]), 30))).toEqual(
      Array.from({ length: 30 }, (_, i) => {
        const n = i + 1;
        return n % 3 === 0 && Number.isInteger(Math.sqrt(n / 3)) ? 1 : 0;
      })
    );
  });

  test('forms that are not positive definite raise pari_err_DOMAIN', () => {
    const bad = fromRows([
      [1, 2],
      [2, 1],
    ]);
    expect(() => qfrep0(bad, 10)).toThrow(PariDomainError);
    expect(() => qfrep0(bad, 10)).toThrow('domain error in minim0: form is not positive definite');
    expect(() => qfrep0(fromRows([[0]]), 10)).toThrow(PariDomainError);
    expect(() => qfrep0(fromRows([[-1]]), 10)).toThrow(PariDomainError);
  });

  test('non-square / non-integral input raises pari_err_TYPE', () => {
    expect(() =>
      qfrep0(
        [
          [1n, 0n],
          [0n, 1n],
          [0n, 0n],
        ],
        5
      )
    ).toThrow(PariTypeError);
    // biome-ignore lint/suspicious/noExplicitAny: deliberately bad input
    expect(() => qfrep0([[1 as any]], 5)).toThrow(PariTypeError);
  });

  test('huge B raises pari_err_PREC, as in minim0_dolll (bibli1.c:1315,1371)', () => {
    const q = fromRows([
      [2, 1],
      [1, 3],
    ]);
    expect(() => qfrep0(q, 10n ** 16n)).toThrow(PariPrecError);
    expect(() => qfrep0(q, 10n ** 16n)).toThrow('precision too low in qfminim');
    expect(() => qfrep0(q, 1n << 70n)).toThrow(PariPrecError);
  });
});

describe('qfrep: flag 1 is the even-norm half of flag 0', () => {
  const forms: number[][][] = [
    A2,
    D4,
    [
      [5, 1, 1],
      [1, 3, 1],
      [1, 1, 1],
    ],
    [
      [4, 1, 0, 2],
      [1, 6, -1, 0],
      [0, -1, 3, 1],
      [2, 0, 1, 7],
    ],
  ];
  test('qfrep(q, B, 1)[i] = qfrep(q, 2B, 0)[2i]', () => {
    for (const rows of forms) {
      const q = fromRows(rows);
      const B = 12;
      const even = qfrep0(q, B, 1);
      const all = qfrep0(q, 2 * B, 0);
      for (let i = 1; i <= B; i++) expect(even[i - 1]).toBe(all[2 * i - 1]!);
    }
  });
});

describe('qfrep: random positive definite forms against exhaustive enumeration', () => {
  /* Deterministic LCG so the test is reproducible. */
  let seed = 0x2f6e2b1;
  const rnd = (m: number) => {
    seed = (seed * 1103515245 + 12345) & 0x7fffffff;
    return seed % m;
  };
  /** G = M^T M + k I, guaranteed positive definite. */
  const randomForm = (n: number, k: number): number[][] => {
    const M: number[][] = Array.from({ length: n }, () =>
      Array.from({ length: n }, () => rnd(7) - 3)
    );
    const G: number[][] = Array.from({ length: n }, () => new Array<number>(n).fill(0));
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        let s = 0;
        for (let r = 0; r < n; r++) s += M[r]![i]! * M[r]![j]!;
        G[i]![j] = s + (i === j ? k : 0);
      }
    }
    return G;
  };

  for (const n of [2, 3, 4]) {
    test(`dimension ${n}`, () => {
      for (let trial = 0; trial < 4; trial++) {
        const G = randomForm(n, 1);
        const B = n === 4 ? 12 : 20;
        expect(nums(qfrep0(fromRows(G), B))).toEqual(bruteforce(G, B));
      }
    });
  }

  test('dimension 5', () => {
    const G = randomForm(5, 2);
    const B = 12;
    expect(nums(qfrep0(fromRows(G), B))).toEqual(bruteforce(G, B));
  });

  test('600 random forms in dimensions 2-5, with and without a skewed basis', () => {
    let checked = 0;
    for (let iter = 0; iter < 600; iter++) {
      const n = 2 + rnd(4);
      const G = randomForm(n, rnd(3));
      if (ZM_det(fromRows(G)) <= 0n) continue; /* degenerate: skip */
      const B = 4 + rnd(n <= 3 ? 30 : 12);
      const bf = bruteforce(G, B);
      expect(nums(qfrep0(fromRows(G), B))).toEqual(bf);
      /* the same lattice in a deliberately bad basis must give the same counts */
      const U = fromRows(
        Array.from({ length: n }, (_, i) =>
          Array.from({ length: n }, (_, j) => (i === j ? 1 : i < j ? rnd(21) - 10 : 0))
        )
      );
      expect(nums(qfrep0(qf_ZM_apply(fromRows(G), U), B))).toEqual(bf);
      checked++;
    }
    expect(checked).toBeGreaterThan(400);
  });
});
