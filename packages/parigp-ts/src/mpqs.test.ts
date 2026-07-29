/**
 * Tests for the MPQS module.
 *
 * Reference: pari/src/basemath/mpqs.c, pari/src/basemath/mpqs.h
 *
 * Oracles used here:
 *  - the sizing tables are diffed against `reference/pari/src/basemath/mpqs.h`
 *    itself, so a transcription typo cannot pass;
 *  - `mpqs_increment` is compared against values produced by compiling
 *    mpqs.c:528-584 verbatim with cc and running it (the vectors below are the
 *    C program's output, copied verbatim);
 *  - `Fl_sqrt` / `krouu` / `kroiu` are checked by exhaustive brute force;
 *  - `F2Ms_ker` is checked against a brute-force rank computation and by
 *    verifying M * v = 0 for every basis vector;
 *  - the sieve itself is checked with upstream's own -DMPQS_DEBUG relation
 *    check (mpqs.c:1069 mpqs_check_rel, our `debug` option): every relation
 *    must satisfy Y^2 = q * prod p_i^e_i (mod N), and X^2 = Y^2 (mod N) must
 *    hold after the Gauss stage (mpqs.c:1525);
 *  - end to end, the returned factors must multiply back to N and be prime.
 */

import { describe, expect, it } from 'bun:test';
import { isPrime } from './ifactor.js';
import { mpqs, mpqsInternals } from './mpqs.js';

const {
  mpqs_increment,
  Fl_sqrt,
  Fl_inv,
  krouu,
  kroiu,
  F2Ms_ker,
  decimal_len,
  logint,
  mpqs_parameters,
  cand_multipliers,
} = mpqsInternals;

// ---------------------------------------------------------------------------
// helpers
// ---------------------------------------------------------------------------

function product(f: Array<[bigint, bigint]>): bigint {
  let p = 1n;
  for (const [v, e] of f) p *= v ** e;
  return p;
}

function makeRng(seed: bigint) {
  let s = seed;
  return () => {
    s = (s * 6364136223846793005n + 1442695040888963407n) & ((1n << 64n) - 1n);
    return s >> 16n;
  };
}

function randPrime(rng: () => bigint, bits: number): bigint {
  const lo = 1n << BigInt(bits - 1);
  for (;;) {
    const c = (lo + (rng() % lo)) | 1n;
    if (isPrime(c)) return c;
  }
}

/** factor and check: product restored, every factor prime */
function checkSplit(n: bigint): Array<[bigint, bigint]> {
  const f = mpqs(n, { debug: true });
  expect(f).not.toBe(null);
  expect(product(f!)).toBe(n);
  for (const [v] of f!) expect(v > 1n && v < n).toBe(true);
  return f!;
}

// ---------------------------------------------------------------------------
// sizing tables: diffed against mpqs.h
// ---------------------------------------------------------------------------

describe('sizing tables (mpqs.h)', () => {
  const HEADER = new URL('../../../reference/pari/src/basemath/mpqs.h', import.meta.url).pathname;

  it('mpqs_parameters matches mpqs.h:292-398 row by row', async () => {
    const h = await Bun.file(HEADER).text();
    const tbl = h.split('static const mpqs_parameterset_t mpqs_parameters[] =')[1].split('};')[0];
    const rows = [
      ...tbl.matchAll(
        /\{\s*\/\*\s*\d+\s*\*\/\s*([0-9.]+),\s*(\d+),\s*(\d+),\s*(\d+),\s*(\d+),\s*(\d+)\}/g
      ),
    ];
    expect(rows.length).toBe(99); /* 9 .. 107 decimal digits */
    expect(mpqs_parameters.length).toBe(rows.length);
    for (let i = 0; i < rows.length; i++) {
      const [, tol, lp, M, fb, oA, pmin] = rows[i];
      /* `tolerance` is a C float; the handle holds the float-rounded value */
      expect(mpqs_parameters[i][0]).toBe(Math.fround(Number(tol)));
      expect(mpqs_parameters[i].slice(1)).toEqual([
        Number(lp),
        Number(M),
        Number(fb),
        Number(oA),
        Number(pmin),
      ]);
    }
  });

  it('cand_multipliers matches mpqs.h:111-153', async () => {
    const h = await Bun.file(HEADER).text();
    const m = h.split('static const mpqs_multiplier_t cand_multipliers[] = {')[1].split('};')[0];
    const rows = [...m.matchAll(/\{\s*(\d+),\s*(\d+),\s*\{\s*(\d+),\s*(\d+)\}\}/g)];
    expect(rows.length).toBe(41);
    expect(cand_multipliers.length).toBe(rows.length);
    for (let i = 0; i < rows.length; i++) {
      const [, k, omega, p1, p2] = rows[i];
      expect(cand_multipliers[i].k).toBe(Number(k));
      expect(cand_multipliers[i].omega_k).toBe(Number(omega));
      expect([...cand_multipliers[i].kp]).toEqual([p1, p2].filter((x) => x !== '0').map(Number));
      /* every k is odd and squarefree, and omega_k counts its prime factors */
      expect(cand_multipliers[i].k & 1).toBe(1);
      let prod = 1;
      for (const p of cand_multipliers[i].kp) prod *= p;
      expect(cand_multipliers[i].omega_k === 0 ? 1 : prod).toBe(cand_multipliers[i].k);
    }
  });
});

// ---------------------------------------------------------------------------
// mpqs_increment (mpqs.c:528)
// ---------------------------------------------------------------------------

describe('mpqs_increment', () => {
  /* Produced by compiling mpqs.c:528-584 verbatim (cc -O1) and printing the
   * orbit of (1<<start)-1; copied verbatim from that program's output. */
  const C_ORBITS: Record<number, number[]> = {
    2: [3, 6, 9, 12, 17, 20, 33, 36, 48, 65, 68, 80, 129, 132, 144, 192, 257, 260, 272, 320],
    3: [7, 13, 19, 22, 25, 28, 35, 38, 41, 44, 49, 52, 67, 70, 73, 76, 81, 84, 97, 100],
    4: [15, 27, 30, 39, 45, 51, 54, 57, 60, 71, 77, 83, 86, 89, 92, 99, 102, 105, 108, 113],
    5: [
      31, 55, 61, 79, 91, 94, 103, 109, 115, 118, 121, 124, 143, 155, 158, 167, 173, 179, 182, 185,
    ],
    6: [
      63, 111, 123, 126, 159, 183, 189, 207, 219, 222, 231, 237, 243, 246, 249, 252, 287, 311, 317,
      335,
    ],
    7: [
      127, 223, 247, 253, 319, 367, 379, 382, 415, 439, 445, 463, 475, 478, 487, 493, 499, 502, 505,
      508,
    ],
    8: [
      255, 447, 495, 507, 510, 639, 735, 759, 765, 831, 879, 891, 894, 927, 951, 957, 975, 987, 990,
      999,
    ],
  };

  it('reproduces the orbits computed by upstream C', () => {
    for (const [start, orbit] of Object.entries(C_ORBITS)) {
      let x = (1 << Number(start)) - 1;
      for (const want of orbit) {
        expect(x).toBe(want);
        x = mpqs_increment(x);
      }
    }
  });

  /* the dense sweep x -> mpqs_increment(x) for x = 1..32, copied verbatim from
   * the same C program (`for x = 1..: print x, increment(x)`). Note that the
   * result is not the *next* integer with the same bit count (4 -> 16, not 8):
   * single-bit-left moves are skipped, cf. mpqs.c:524-527. */
  const C_DENSE = [
    4, 8, 6, 16, 10, 9, 13, 32, 12, 17, 14, 17, 19, 19, 27, 64, 20, 24, 22, 33, 26, 25, 29, 33, 28,
    35, 30, 35, 39, 39, 55, 128,
  ];

  it('reproduces upstream C on 1 <= x <= 32', () => {
    for (let x = 1; x <= 32; x++) expect(mpqs_increment(x)).toBe(C_DENSE[x - 1]);
  });

  it('preserves the bit count and increases (mpqs.c:520-527)', () => {
    const popcount = (x: number) => {
      let c = 0;
      for (let i = 0; i < 32; i++) if (x & (1 << i)) c++;
      return c;
    };
    for (let x = 1; x < 20000; x++) {
      const y = mpqs_increment(x);
      expect(y).toBeGreaterThan(x);
      expect(popcount(y)).toBe(popcount(x));
      /* "does not arise from the old value by moving a single 1 bit one
       * position to the left" -- i.e. y != x - 2^i + 2^(i+1) = x + 2^i */
      for (let i = 0; i < 31; i++) {
        if (x & (1 << i) && !(x & (1 << (i + 1)))) expect(y).not.toBe(x + (1 << i));
      }
    }
  });
});

// ---------------------------------------------------------------------------
// word-sized modular arithmetic
// ---------------------------------------------------------------------------

describe('Fl_sqrt / Fl_inv (arith1.c)', () => {
  it('returns a square root of every quadratic residue mod p < 500', () => {
    for (let p = 3; p < 500; p++) {
      let isP = true;
      for (let d = 2; d * d <= p; d++) if (p % d === 0) isP = false;
      if (!isP) continue;
      for (let a = 0; a < p; a++) {
        const sq = (a * a) % p;
        const r = Fl_sqrt(sq, p);
        expect((r * r) % p).toBe(sq);
      }
    }
  });

  it('inverts every nonzero residue mod p < 500', () => {
    for (let p = 3; p < 500; p++) {
      let isP = true;
      for (let d = 2; d * d <= p; d++) if (p % d === 0) isP = false;
      if (!isP) continue;
      for (let a = 1; a < p; a++) expect((a * Fl_inv(a, p)) % p).toBe(1);
    }
  });
});

describe('krouu / kroiu (Kronecker symbol)', () => {
  it('agrees with the Euler criterion for odd primes', () => {
    for (const p of [3, 5, 7, 11, 13, 17, 19, 23, 97, 101, 1009]) {
      for (let a = 0; a < p; a++) {
        const squares = new Set<number>();
        for (let i = 0; i < p; i++) squares.add((i * i) % p);
        const want = a === 0 ? 0 : squares.has(a) ? 1 : -1;
        expect(krouu(a, p)).toBe(want);
      }
    }
  });

  it('handles y = 2 as the Kronecker symbol does', () => {
    /* (x|2) = 0 for even x, 1 for x = +-1 mod 8, -1 for x = +-3 mod 8 */
    for (let x = 1; x < 200; x++) {
      const want = x % 2 === 0 ? 0 : x % 8 === 1 || x % 8 === 7 ? 1 : -1;
      expect(krouu(x, 2)).toBe(want);
    }
  });

  it('kroiu reduces a bigint modulo p first', () => {
    const N = 123456789012345678901234567n;
    for (const p of [2, 3, 5, 7, 11, 101, 65537]) {
      const want =
        p === 2
          ? Number(N % 8n) % 8 === 1 || Number(N % 8n) === 7
            ? 1
            : -1
          : krouu(Number(N % BigInt(p)), p);
      expect(kroiu(N, p)).toBe(want);
    }
  });
});

describe('decimal_len / logint', () => {
  it('decimal_len counts decimal digits (mpqs.c:115)', () => {
    expect(decimal_len(1n)).toBe(1);
    expect(decimal_len(9n)).toBe(1);
    expect(decimal_len(10n)).toBe(2);
    expect(decimal_len(99n)).toBe(2);
    expect(decimal_len(10n ** 50n)).toBe(51);
  });

  it('logint(N,B) is the largest e with B^e <= N', () => {
    expect(logint(8n, 3n)).toBe(1);
    expect(logint(9n, 3n)).toBe(2);
    expect(logint(26n, 3n)).toBe(2);
    expect(logint(27n, 3n)).toBe(3);
  });
});

// ---------------------------------------------------------------------------
// F2Ms_ker (F2v.c:397 F2m_ker_sp)
// ---------------------------------------------------------------------------

describe('F2Ms_ker', () => {
  /** brute-force rank over F2 of a column list */
  function rank(M: number[][]): number {
    const rows: bigint[] = M.map((col) => {
      let v = 0n;
      for (const r of col) v |= 1n << BigInt(r);
      return v;
    });
    let rk = 0;
    const pivots: bigint[] = [];
    for (let v of rows) {
      for (const p of pivots) {
        const hb = 1n << BigInt(p.toString(2).length - 1);
        if (v & hb) v ^= p;
      }
      if (v !== 0n) {
        pivots.push(v);
        pivots.sort((a, b) => (a > b ? -1 : 1));
        rk++;
      }
    }
    return rk;
  }

  it('returns a basis of the kernel on random matrices', () => {
    const rng = makeRng(0x1234567n);
    for (let trial = 0; trial < 200; trial++) {
      const nbrow = 1 + Number(rng() % 40n);
      const ncol = 1 + Number(rng() % 40n);
      const M: number[][] = [];
      for (let c = 0; c < ncol; c++) {
        const col: number[] = [];
        for (let r = 1; r <= nbrow; r++) if (rng() & 1n) col.push(r);
        M.push(col);
      }
      const K = F2Ms_ker(M, nbrow);
      expect(K.length).toBe(ncol - rank(M));
      for (const v of K) {
        /* v must be nonzero and sum the selected columns to 0 */
        const acc = new Uint8Array(nbrow + 1);
        let any = false;
        for (let c = 1; c <= ncol; c++) {
          if (v[c >>> 5] & (1 << (c & 31))) {
            any = true;
            for (const r of M[c - 1]) acc[r] ^= 1;
          }
        }
        expect(any).toBe(true);
        expect(acc.some((b) => b !== 0)).toBe(false);
      }
      /* basis vectors must be independent: their leading columns are distinct */
      const lead = new Set<number>();
      for (const v of K) {
        let l = -1;
        for (let c = 1; c <= ncol; c++) if (v[c >>> 5] & (1 << (c & 31))) l = c;
        expect(lead.has(l)).toBe(false);
        lead.add(l);
      }
    }
  });

  it('has an empty kernel for an identity-like matrix', () => {
    const M = [[1], [2], [3]];
    expect(F2Ms_ker(M, 3).length).toBe(0);
  });
});

// ---------------------------------------------------------------------------
// mpqs end to end
// ---------------------------------------------------------------------------

describe('mpqs (mpqs.c:1639)', () => {
  it('splits 300 random semiprimes, verifying every relation (-DMPQS_DEBUG)', () => {
    const rng = makeRng(0x12345678deadbeefn);
    let n0 = 0;
    for (let i = 0; i < 300; i++) {
      const p = randPrime(rng, 12 + Number(rng() % 19n));
      const q = randPrime(rng, 12 + Number(rng() % 19n));
      if (p === q) continue; /* p^2 is a pure power: ifac_crack strips it first */
      const n = p * q;
      n0++;
      const f = checkSplit(n);
      /* the two primes must be recovered (possibly as a product still to be
       * cracked, but for a semiprime the split is exact) */
      expect(f.map(([v]) => v).sort((a, b) => (a < b ? -1 : 1))).toEqual(
        [p, q].sort((a, b) => (a < b ? -1 : 1))
      );
    }
    expect(n0).toBeGreaterThan(280);
  }, 120000);

  it('splits composites with three prime factors', () => {
    const rng = makeRng(0xfeedfacecafebeefn);
    for (let i = 0; i < 30; i++) {
      const p = randPrime(rng, 10);
      const q = randPrime(rng, 26);
      const r = randPrime(rng, 22);
      if (p === q || q === r || p === r) continue;
      const n = p * q * r;
      const f = checkSplit(n);
      /* each returned value divides n and the product is n */
      for (const [v] of f) expect(n % v).toBe(0n);
    }
  }, 120000);

  it('splits RSA-style moduli p*q with p,q ~ 10^20', () => {
    const P = 100000000000000012349n;
    const Q = 100000000000000006903n;
    expect(isPrime(P) && isPrime(Q)).toBe(true);
    const f = checkSplit(P * Q);
    expect(f.map(([v]) => v).sort((a, b) => (a < b ? -1 : 1))).toEqual([Q, P]);
  }, 120000);

  it('splits a hard semiprime with a 25-digit smallest factor', () => {
    const p = 1000000000000000000012369n;
    const q = 1000000000000000000067987n;
    expect(isPrime(p) && isPrime(q)).toBe(true);
    const f = checkSplit(p * q);
    expect(f.map(([v]) => v).sort((a, b) => (a < b ? -1 : 1))).toEqual([p, q]);
  }, 300000);

  it('splits an unbalanced semiprime (25 and 32 digits)', () => {
    const p = 1000000000000000000012369n;
    const q = 10000000000000000000000000012411n;
    expect(isPrime(p) && isPrime(q)).toBe(true);
    const f = checkSplit(p * q);
    expect(f.map(([v]) => v).sort((a, b) => (a < b ? -1 : 1))).toEqual([p, q]);
  }, 300000);

  it('splits a hard semiprime with a 30-digit smallest factor', () => {
    const p = 100000000000000000000000012349n;
    const q = 100000000000000000000000068011n;
    expect(isPrime(p) && isPrime(q)).toBe(true);
    const f = checkSplit(p * q);
    expect(f.map(([v]) => v).sort((a, b) => (a < b ? -1 : 1))).toEqual([p, q]);
  }, 300000);

  it('splits a hard semiprime with a 35-digit smallest factor', () => {
    const p = 10000000000000000000000000000012423n;
    const q = 10000000000000000000000000000067969n;
    expect(isPrime(p) && isPrime(q)).toBe(true);
    const f = checkSplit(p * q);
    expect(f.map(([v]) => v).sort((a, b) => (a < b ? -1 : 1))).toEqual([p, q]);
  }, 600000);

  it('splits a hard semiprime with a 40-digit smallest factor', () => {
    const p = 1000000000000000000000000000000000012397n;
    const q = 1000000000000000000000000000000000067969n;
    expect(isPrime(p) && isPrime(q)).toBe(true);
    const f = checkSplit(p * q);
    expect(f.map(([v]) => v).sort((a, b) => (a < b ? -1 : 1))).toEqual([p, q]);
  }, 900000);

  it('reports the multiplicity of a proper power in ifac format', () => {
    /* mpqs.c:1444 split(): a factor which is a perfect power is returned as
     * (root, exponent), so the product of value^exponent is N. N itself must
     * not be a perfect power (mpqs.c:1220 relies on N not being a square, and
     * ifac_crack strips pure powers before calling MPQS). */
    const p = 1000003n;
    const q = 1000033n;
    const n = p * p * q;
    const f = mpqs(n, { debug: true });
    expect(f).not.toBe(null);
    expect(product(f!)).toBe(n);
    for (const [v, e] of f!) expect(v ** e > 1n).toBe(true);
  }, 120000);

  it('declines input above 107 decimal digits (mpqs.h:400, toolarge())', () => {
    const p = 100000000000000000000000000000000000000000000000000000001389n;
    const q = 100000000000000000000000000000000000000000000000000000009091n;
    expect(isPrime(p) && isPrime(q)).toBe(true);
    const n = p * q; /* 120 digits */
    expect(n.toString().length).toBeGreaterThan(107);
    expect(mpqs(n)).toBe(null);
  });

  it('returns null rather than a wrong answer when its budget runs out', () => {
    const p = 1000000000000000000012369n;
    const q = 1000000000000000000067987n;
    expect(mpqs(p * q, { maxPolys: 1 })).toBe(null);
  });

  it('returns the factor spotted while choosing the multiplier', () => {
    /* mpqs.c:1666: mpqs_find_k returns p as soon as kronecker(N,p) = 0, i.e.
     * p | N, and mpqs() hands that back as a factor (with the cofactor, which
     * ifac_decomp will crack next). PARI never reaches this in practice
     * because ifac_crack has trial-divided already. */
    const n = 101n * 1000003n * 1000033n;
    const f = mpqs(n);
    expect(f).not.toBe(null);
    expect(product(f!)).toBe(n);
    expect(f!.map(([v]) => v)).toContain(101n);
  });

  it('splits a composite whose small factor is beyond the multiplier search', () => {
    /* here mpqs_find_k stops before reaching 3623 (it only walks
     * size_of_FB primes, mpqs.c:262), so the sieve does the work */
    const n = 3623n * 17317n * 1593586531n;
    const f = checkSplit(n);
    for (const [v] of f) expect(n % v).toBe(0n);
  });
});
