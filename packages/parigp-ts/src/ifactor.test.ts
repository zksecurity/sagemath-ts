/**
 * Tests for ifactor module
 *
 * Reference: pari/src/basemath/ifactor1.c, pari/src/basemath/ispower.c
 *
 * Oracles used here:
 *  - exhaustive brute-force trial division for small inputs,
 *  - the identity `prod(p^e) === n` plus a BPSW test on every returned factor,
 *  - published factorizations of Mersenne/Fermat numbers.
 */

import { describe, expect, it } from 'bun:test';
import {
  NotImplementedError,
  Z_factor,
  Z_iroot,
  Z_isanypower,
  Z_issquareall,
  Z_pollardbrent,
  ellfacteur,
  factoru,
  formatFactorization,
  isPrime,
  is_357_power,
  is_kth_power,
  isprimepower,
  isqrt,
  pollardbrent,
  squfof,
  tridiv_bound,
} from './ifactor.js';

// ---------------------------------------------------------------------------
// helpers
// ---------------------------------------------------------------------------

/** Brute-force factorization oracle (trial division up to sqrt). */
function naiveFactor(n: bigint): Array<[bigint, bigint]> {
  const out: Array<[bigint, bigint]> = [];
  if (n < 0n) {
    out.push([-1n, 1n]);
    n = -n;
  }
  for (let p = 2n; p * p <= n; p++) {
    if (n % p === 0n) {
      let e = 0n;
      while (n % p === 0n) {
        n /= p;
        e++;
      }
      out.push([p, e]);
    }
  }
  if (n > 1n) out.push([n, 1n]);
  return out;
}

function product(f: Array<[bigint, bigint]>): bigint {
  let prod = 1n;
  for (const [p, e] of f) prod *= p ** e;
  return prod;
}

/** Every entry is a (BPSW) prime, apart from the leading sign. */
function allPrime(f: Array<[bigint, bigint]>): boolean {
  return f.every(([p], i) => (i === 0 && p === -1n) || isPrime(p));
}

/** Deterministic xorshift64 PRNG, so the random sweeps are reproducible. */
function makeRng(seed: bigint) {
  const MASK = (1n << 64n) - 1n;
  let state = seed;
  return () => {
    state ^= (state << 13n) & MASK;
    state ^= state >> 7n;
    state ^= (state << 17n) & MASK;
    return state;
  };
}

function randPrime(rng: () => bigint, bits: number): bigint {
  for (;;) {
    const c = (rng() % (1n << BigInt(bits))) | (1n << BigInt(bits - 1)) | 1n;
    if (isPrime(c)) return c;
  }
}

// ---------------------------------------------------------------------------

describe('isPrime', () => {
  it('correctly identifies small primes', () => {
    const primes = [2n, 3n, 5n, 7n, 11n, 13n, 17n, 19n, 23n, 29n, 31n];
    for (const p of primes) {
      expect(isPrime(p)).toBe(true);
    }
  });

  it('correctly identifies small non-primes', () => {
    const nonPrimes = [0n, 1n, 4n, 6n, 8n, 9n, 10n, 12n, 14n, 15n, 16n];
    for (const n of nonPrimes) {
      expect(isPrime(n)).toBe(false);
    }
  });

  it('correctly identifies medium primes', () => {
    const primes = [101n, 103n, 107n, 109n, 113n, 127n, 131n, 137n, 139n, 149n];
    for (const p of primes) {
      expect(isPrime(p)).toBe(true);
    }
  });

  it('correctly identifies primes ≡ 5 (mod 6)', () => {
    const primes = [5n, 11n, 17n, 23n, 29n, 41n, 47n, 53n, 59n, 71n, 83n, 89n];
    for (const p of primes) {
      expect(isPrime(p)).toBe(true);
    }
  });

  it('correctly identifies 761 and 1877 as prime', () => {
    expect(isPrime(761n)).toBe(true);
    expect(isPrime(1877n)).toBe(true);
  });

  it('correctly identifies large primes', () => {
    const primes = [
      104729n, // 10000th prime
      1000003n, // Just over 10^6
      15485863n, // 1 millionth prime
      2147483647n, // 2^31 - 1 (Mersenne prime)
      170141183460469231731687303715884105727n, // 2^127 - 1
    ];
    for (const p of primes) {
      expect(isPrime(p)).toBe(true);
    }
  });

  it('correctly identifies Carmichael numbers as composite', () => {
    const carmichaels = [561n, 1105n, 1729n, 2465n, 2821n, 6601n, 8911n];
    for (const n of carmichaels) {
      expect(isPrime(n)).toBe(false);
    }
  });

  it('agrees with brute force on 0..20000', () => {
    for (let n = 0n; n <= 20000n; n++) {
      const want = n >= 2n && naiveFactor(n).length === 1 && naiveFactor(n)[0][1] === 1n;
      expect(isPrime(n)).toBe(want);
    }
  });
});

describe('exact roots', () => {
  it('isqrt is the exact integer square root', () => {
    for (let n = 0n; n < 1000n; n++) {
      const r = isqrt(n);
      expect(r * r <= n && (r + 1n) * (r + 1n) > n).toBe(true);
    }
    const big = 10n ** 100n + 12345n;
    const r = isqrt(big);
    expect(r * r <= big && (r + 1n) * (r + 1n) > big).toBe(true);
  });

  it('Z_iroot is the exact integer k-th root', () => {
    for (const k of [2, 3, 5, 7, 11]) {
      for (const base of [2n, 3n, 10n, 1000003n, 10n ** 12n + 39n]) {
        expect(Z_iroot(base ** BigInt(k), k)).toEqual([base, true]);
        expect(Z_iroot(base ** BigInt(k) + 1n, k)[1]).toBe(false);
        expect(Z_iroot(base ** BigInt(k) - 1n, k)).toEqual([base - 1n, false]);
      }
    }
  });

  it('Z_issquareall detects squares', () => {
    for (let n = 0n; n < 5000n; n++) {
      const r = Z_issquareall(n);
      const s = isqrt(n);
      expect(r).toBe(s * s === n ? s : null);
    }
    expect(Z_issquareall(1000003n ** 2n)).toBe(1000003n);
  });

  it('is_357_power detects cubes, 5th and 7th powers', () => {
    expect(is_357_power(27n, 7)[0]).toBe(3);
    expect(is_357_power(3n ** 5n, 7)[0]).toBe(5);
    expect(is_357_power(3n ** 7n, 7)[0]).toBe(7);
    // priority to the higher power: 3^21 is reported as a 7th power
    expect(is_357_power(3n ** 21n, 7)).toEqual([7, 3n ** 3n, 7]);
    expect(is_357_power(3n ** 4n, 7)[0]).toBe(0);
    // the mask of failed exponents is cleared
    expect(is_357_power(11n, 7)[2]).toBe(0);
  });

  it('is_kth_power agrees with exact roots', () => {
    expect(is_kth_power(1000003n ** 11n, 11)).toBe(1000003n);
    expect(is_kth_power(1000003n ** 11n + 1n, 11)).toBe(null);
    expect(is_kth_power(2n ** 100n, 13)).toBe(null);
    expect(is_kth_power(7n ** 13n, 13)).toBe(7n);
  });
});

describe('Z_isanypower', () => {
  it('agrees with brute force on 2..20000', () => {
    for (let n = 2n; n <= 20000n; n++) {
      const f = naiveFactor(n);
      let g = 0n;
      for (const [, e] of f) {
        let [a, b] = [g, e];
        while (b) [a, b] = [b, a % b];
        g = a;
      }
      const [k, y] = Z_isanypower(n);
      expect(k).toBe(g > 1n ? Number(g) : 0);
      if (k) expect(y ** BigInt(k)).toBe(n);
    }
  });

  it('handles large perfect powers', () => {
    const p = 100000000000000000000117n; // 24-digit prime
    expect(Z_isanypower(p ** 6n)).toEqual([6, p]);
    expect(Z_isanypower((p * 1000003n) ** 3n)).toEqual([3, p * 1000003n]);
    expect(Z_isanypower(p ** 6n + 1n)[0]).toBe(0);
    expect(Z_isanypower(2n ** 64n)).toEqual([64, 2n]);
  });
});

describe('isprimepower', () => {
  it('agrees with brute force on 2..20000', () => {
    for (let n = 2n; n <= 20000n; n++) {
      const f = naiveFactor(n);
      const want: [bigint, number] | null = f.length === 1 ? [f[0][0], Number(f[0][1])] : null;
      expect(isprimepower(n)).toEqual(want);
    }
  });

  it('rejects 0, 1 and negatives', () => {
    expect(isprimepower(0n)).toBe(null);
    expect(isprimepower(1n)).toBe(null);
    expect(isprimepower(-4n)).toBe(null);
  });

  it('never factors: p*q with 24-digit p, q is not a prime power', () => {
    const p = 100000000000000000000117n;
    const q = 100000000000000000000213n;
    expect(isPrime(p) && isPrime(q)).toBe(true);
    expect(isprimepower(p * q)).toBe(null);
    expect(isprimepower(p * p)).toEqual([p, 2]);
    expect(isprimepower(p ** 3n)).toEqual([p, 3]);
    expect(isprimepower(p ** 6n)).toEqual([p, 6]);
    expect(isprimepower(1000003n ** 2n)).toEqual([1000003n, 2]);
    expect(isprimepower(1000003n ** 5n)).toEqual([1000003n, 5]);
  });
});

describe('tridiv_bound', () => {
  // Reference: ifactor1.c:3137 (tridiv_boundu) and ifactor1.c:3158
  it('matches PARI(64-bit) for word-sized input', () => {
    expect(tridiv_bound(1n << 20n)).toBe(1 << 12); // expu = 20 < 30
    expect(tridiv_bound(1n << 30n)).toBe(1 << 13); // 30 <= e < 34
    expect(tridiv_bound(1n << 34n)).toBe(1 << 14);
    expect(tridiv_bound(1n << 37n)).toBe(1 << 15);
    expect(tridiv_bound(1n << 42n)).toBe(1 << 16);
    expect(tridiv_bound(1n << 47n)).toBe(1 << 17);
    expect(tridiv_bound(1n << 56n)).toBe(1 << 19);
    expect(tridiv_bound(1n << 62n)).toBe(1 << 18);
  });

  it('matches PARI for multiword input', () => {
    expect(tridiv_bound(1n << 100n)).toBe((101 - 16) << 10);
    expect(tridiv_bound(1n << 600n)).toBe(1 << 19);
  });
});

describe('squfof', () => {
  // Reference: ifactor1.c:1474
  it('splits random semiprimes below 2^59', () => {
    const rng = makeRng(0x9e3779b97f4a7c15n);
    let split = 0;
    for (let i = 0; i < 200; i++) {
      const p = randPrime(rng, 28);
      const q = randPrime(rng, 28);
      if (p === q) continue;
      const n = p * q;
      const r = squfof(n);
      expect(r).not.toBe(null);
      let prod = 1n;
      for (const f of r!) prod *= f;
      expect(prod).toBe(n);
      expect(r!.includes(p) || r!.includes(q)).toBe(true);
      split++;
    }
    expect(split).toBeGreaterThan(190);
  });

  it('splits unbalanced semiprimes', () => {
    const rng = makeRng(0x2545f4914f6cdd1dn);
    for (const [pb, qb] of [
      [10, 40],
      [14, 30],
      [20, 30],
    ] as const) {
      for (let i = 0; i < 25; i++) {
        const p = randPrime(rng, pb);
        const q = randPrime(rng, qb);
        const n = p * q;
        const r = squfof(n);
        if (r === null) continue;
        let prod = 1n;
        for (const f of r) prod *= f;
        expect(prod).toBe(n);
      }
    }
  });

  it('declines input at or above 2^59', () => {
    expect(squfof((1n << 59n) + 1n)).toBe(null);
  });
});

describe('pollardbrent', () => {
  // Reference: ifactor1.c:1184 (pollardbrent_i), 1379 (Z_pollardbrent)
  it('splits semiprimes and returns a correct product', () => {
    const rng = makeRng(0xdeadbeefcafebaben);
    for (let i = 0; i < 30; i++) {
      const p = randPrime(rng, 24);
      const q = randPrime(rng, 40);
      const n = p * q;
      const r = Z_pollardbrent(n, 1 << 12, 0);
      expect(r).not.toBe(null);
      let prod = 1n;
      for (const f of r!) prod *= f;
      expect(prod).toBe(n);
      expect(r!.every((f) => f > 1n && f < n)).toBe(true);
    }
  }, 60000);

  it("every split found with PARI's own round budget is correct", () => {
    // PARI's budget (ifactor1.c:1370) is deliberately small: it is tuned to
    // hand the hard cases over to MPQS, so `pollardbrent` is allowed to
    // decline. What it must never do is return a wrong split.
    const rng = makeRng(0x5deece66dn);
    let found = 0;
    for (let i = 0; i < 20; i++) {
      const p = randPrime(rng, 16);
      const q = randPrime(rng, 70);
      const n = p * q;
      const r = pollardbrent(n);
      if (r === null) continue;
      found++;
      let prod = 1n;
      for (const f of r) prod *= f;
      expect(prod).toBe(n);
      expect(r.every((f) => f > 1n && f < n)).toBe(true);
    }
    expect(found).toBeGreaterThan(0);
  }, 60000);

  it('gives up (rather than lying) when its budget is exhausted', () => {
    const p = 100000000000000000000117n;
    const q = 100000000000000000000213n;
    expect(Z_pollardbrent(p * q, 4, 0)).toBe(null);
  });
});

describe('ellfacteur (ECM)', () => {
  // Reference: ifactor1.c:1038
  it('finds a 9-digit factor of a 192-bit composite', () => {
    const p = 982451653n;
    const r = 170141183460469231731687303715884105727n; // 2^127-1
    const N = p * 32416190071n * r;
    const g = ellfacteur(N, false);
    expect(g).not.toBe(null);
    expect(N % g!).toBe(0n);
    expect(g! > 1n && g! < N).toBe(true);
  }, 60000);

  it('declines inputs below 140 bits when not insisting', () => {
    // ifactor1.c:1091: "number too small to justify this stage"
    expect(ellfacteur(1000003n * 1000033n, false)).toBe(null);
  });
});

describe('Z_factor', () => {
  it('factors small primes correctly', () => {
    expect(Z_factor(2n)).toEqual([[2n, 1n]]);
    expect(Z_factor(3n)).toEqual([[3n, 1n]]);
    expect(Z_factor(5n)).toEqual([[5n, 1n]]);
    expect(Z_factor(7n)).toEqual([[7n, 1n]]);
    expect(Z_factor(11n)).toEqual([[11n, 1n]]);
  });

  it('factors 1 correctly', () => {
    expect(Z_factor(1n)).toEqual([]);
  });

  it('factors negative numbers correctly', () => {
    expect(Z_factor(-1n)).toEqual([[-1n, 1n]]);
    expect(Z_factor(-2n)).toEqual([
      [-1n, 1n],
      [2n, 1n],
    ]);
    expect(Z_factor(-12n)).toEqual([
      [-1n, 1n],
      [2n, 2n],
      [3n, 1n],
    ]);
  });

  it('factors prime powers correctly', () => {
    expect(Z_factor(4n)).toEqual([[2n, 2n]]);
    expect(Z_factor(8n)).toEqual([[2n, 3n]]);
    expect(Z_factor(9n)).toEqual([[3n, 2n]]);
    expect(Z_factor(27n)).toEqual([[3n, 3n]]);
    expect(Z_factor(32n)).toEqual([[2n, 5n]]);
    expect(Z_factor(1000003n ** 2n)).toEqual([[1000003n, 2n]]);
    expect(Z_factor(1000003n ** 7n)).toEqual([[1000003n, 7n]]);
  });

  it('factors composite numbers correctly', () => {
    expect(Z_factor(6n)).toEqual([
      [2n, 1n],
      [3n, 1n],
    ]);
    expect(Z_factor(12n)).toEqual([
      [2n, 2n],
      [3n, 1n],
    ]);
    expect(Z_factor(30n)).toEqual([
      [2n, 1n],
      [3n, 1n],
      [5n, 1n],
    ]);
    expect(Z_factor(100n)).toEqual([
      [2n, 2n],
      [5n, 2n],
    ]);
    expect(Z_factor(360n)).toEqual([
      [2n, 3n],
      [3n, 2n],
      [5n, 1n],
    ]);
  });

  it('factors 49993895 correctly (the original failing case)', () => {
    // 49993895 = 5 * 7 * 761 * 1877
    const result = Z_factor(49993895n);
    expect(result).toEqual([
      [5n, 1n],
      [7n, 1n],
      [761n, 1n],
      [1877n, 1n],
    ]);
    expect(product(result)).toBe(49993895n);
  });

  it('correctly factors numbers with primes ≡ 5 (mod 6)', () => {
    expect(Z_factor(55n)).toEqual([
      [5n, 1n],
      [11n, 1n],
    ]);
    expect(Z_factor(85n)).toEqual([
      [5n, 1n],
      [17n, 1n],
    ]);
    expect(Z_factor(115n)).toEqual([
      [5n, 1n],
      [23n, 1n],
    ]);
    expect(Z_factor(187n)).toEqual([
      [11n, 1n],
      [17n, 1n],
    ]);
    expect(Z_factor(1309n)).toEqual([
      [7n, 1n],
      [11n, 1n],
      [17n, 1n],
    ]);
  });

  it('throws on zero', () => {
    expect(() => Z_factor(0n)).toThrow();
  });

  it('agrees with brute force exhaustively on 1..5000', () => {
    for (let n = 1n; n <= 5000n; n++) {
      expect(Z_factor(n)).toEqual(naiveFactor(n));
      expect(Z_factor(-n)).toEqual(naiveFactor(-n));
    }
  }, 30000);

  it('agrees with brute force on 1000 random n < 10^9', () => {
    const rng = makeRng(0x123456789abcdefn);
    for (let i = 0; i < 1000; i++) {
      const n = (rng() % 10n ** 9n) + 1n;
      expect(Z_factor(n)).toEqual(naiveFactor(n));
    }
  }, 30000);

  it('factors 2^62-1 completely (each factor prime)', () => {
    // 2^62-1 = 3 * 715827883 * 2147483647; beyond the trial division bound,
    // this used to be returned as "3 * <composite declared prime>".
    const n = 2n ** 62n - 1n;
    const f = Z_factor(n);
    expect(f).toEqual([
      [3n, 1n],
      [715827883n, 1n],
      [2147483647n, 1n],
    ]);
    expect(product(f)).toBe(n);
    expect(allPrime(f)).toBe(true);
  }, 30000);

  it('matches published factorizations of Mersenne/Fermat numbers', () => {
    const cases: Array<[bigint, Array<[bigint, bigint]>]> = [
      // F6 = 2^64+1 = 274177 * 67280421310721 (Landry, 1880)
      [
        2n ** 64n + 1n,
        [
          [274177n, 1n],
          [67280421310721n, 1n],
        ],
      ],
      // M67 = 193707721 * 761838257287 (Cole, 1903)
      [
        2n ** 67n - 1n,
        [
          [193707721n, 1n],
          [761838257287n, 1n],
        ],
      ],
      // M71 = 228479 * 48544121 * 212885833
      [
        2n ** 71n - 1n,
        [
          [228479n, 1n],
          [48544121n, 1n],
          [212885833n, 1n],
        ],
      ],
      // M101 = 7432339208719 * 341117531003194129
      [
        2n ** 101n - 1n,
        [
          [7432339208719n, 1n],
          [341117531003194129n, 1n],
        ],
      ],
    ];
    for (const [n, want] of cases) {
      const f = Z_factor(n);
      expect(f).toEqual(want);
      expect(product(f)).toBe(n);
      expect(allPrime(f)).toBe(true);
    }
  }, 120000);

  it('factors composite prime powers whose base is composite', () => {
    const m = 1000003n * 1000033n;
    expect(Z_factor(m ** 4n)).toEqual([
      [1000003n, 4n],
      [1000033n, 4n],
    ]);
  }, 30000);

  it('2000 random semiprimes: product is restored and factors are prime', () => {
    const rng = makeRng(0xa5a5a5a5deadbeefn);
    for (let i = 0; i < 2000; i++) {
      const bitsP = 12 + Number(rng() % 19n); // 12..30 bits
      const bitsQ = 12 + Number(rng() % 19n);
      const p = randPrime(rng, bitsP);
      const q = randPrime(rng, bitsQ);
      const n = p * q;
      const f = Z_factor(n);
      expect(product(f)).toBe(n);
      expect(allPrime(f)).toBe(true);
      // the factorization must be exactly {p, q}
      const want =
        p === q
          ? [[p, 2n]]
          : p < q
            ? [
                [p, 1n],
                [q, 1n],
              ]
            : [
                [q, 1n],
                [p, 1n],
              ];
      expect(f).toEqual(want as Array<[bigint, bigint]>);
    }
  }, 300000);

  it('all factors are prime', () => {
    const testCases = [100n, 1000n, 10000n, 123456n, 654321n, 999999n, 49993895n];
    for (const n of testCases) {
      expect(allPrime(Z_factor(n))).toBe(true);
    }
  });

  it('factorization multiplies back to original', () => {
    const testCases = [1n, 2n, 6n, 12n, 100n, 1000n, 12345n, 49993895n, 123456789n];
    for (const n of testCases) {
      expect(product(Z_factor(n))).toBe(n);
    }
  });

  it('reports failure instead of declaring a composite prime', () => {
    // 48-digit semiprime: out of reach without MPQS. PARI would sieve it;
    // we must throw rather than return [[n, 1]].
    const p = 100000000000000000000117n;
    const q = 100000000000000000000213n;
    let threw: unknown = null;
    try {
      Z_factor(p * q, { ecmRounds: 0 });
    } catch (e) {
      threw = e;
    }
    expect(threw).toBeInstanceOf(NotImplementedError);
    expect((threw as Error).message).toContain('mpqs');
  }, 120000);
});

describe('formatFactorization', () => {
  it('formats empty factorization as 1', () => {
    expect(formatFactorization([])).toBe('1');
  });

  it('formats single prime correctly', () => {
    expect(formatFactorization([[2n, 1n]])).toBe('2');
    expect(formatFactorization([[3n, 2n]])).toBe('3^2');
  });

  it('formats multiple primes correctly', () => {
    expect(
      formatFactorization([
        [2n, 2n],
        [3n, 1n],
      ])
    ).toBe('2^2 * 3');
    expect(
      formatFactorization([
        [2n, 1n],
        [3n, 1n],
        [5n, 1n],
      ])
    ).toBe('2 * 3 * 5');
  });
});

describe('factoru', () => {
  it('is an alias for Z_factor', () => {
    expect(factoru(12n)).toEqual(Z_factor(12n));
    expect(factoru(100n)).toEqual(Z_factor(100n));
    expect(factoru(49993895n)).toEqual(Z_factor(49993895n));
  });
});
