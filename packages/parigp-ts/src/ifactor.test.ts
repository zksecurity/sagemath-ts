/**
 * Tests for ifactor module
 *
 * Reference: pari/src/basemath/ifactor1.c
 */

import { describe, expect, it } from 'bun:test';
import { Z_factor, factoru, formatFactorization, isPrime } from './ifactor.js';

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
    // These are the primes that were missed by the buggy wheel factorization
    const primes = [5n, 11n, 17n, 23n, 29n, 41n, 47n, 53n, 59n, 71n, 83n, 89n];
    for (const p of primes) {
      expect(isPrime(p)).toBe(true);
    }
  });

  it('correctly identifies 761 and 1877 as prime', () => {
    // These specific primes were missed in the original bug
    expect(isPrime(761n)).toBe(true);
    expect(isPrime(1877n)).toBe(true);
  });

  it('correctly identifies large primes', () => {
    const primes = [
      104729n, // 10000th prime
      1000003n, // Just over 10^6
      15485863n, // 1 millionth prime
      2147483647n, // 2^31 - 1 (Mersenne prime)
    ];
    for (const p of primes) {
      expect(isPrime(p)).toBe(true);
    }
  });

  it('correctly identifies Carmichael numbers as composite', () => {
    // Carmichael numbers pass Fermat test but are composite
    const carmichaels = [561n, 1105n, 1729n, 2465n, 2821n, 6601n, 8911n];
    for (const n of carmichaels) {
      expect(isPrime(n)).toBe(false);
    }
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
    // This was the counterexample that exposed the wheel factorization bug
    // 49993895 = 5 * 7 * 761 * 1877
    const result = Z_factor(49993895n);
    expect(result).toEqual([
      [5n, 1n],
      [7n, 1n],
      [761n, 1n],
      [1877n, 1n],
    ]);

    // Verify by multiplying back
    let product = 1n;
    for (const [p, e] of result) {
      product *= p ** e;
    }
    expect(product).toBe(49993895n);
  });

  it('correctly factors numbers with primes ≡ 5 (mod 6)', () => {
    // Test products of primes that are ≡ 5 (mod 6)
    expect(Z_factor(55n)).toEqual([
      [5n, 1n],
      [11n, 1n],
    ]); // 5 * 11
    expect(Z_factor(85n)).toEqual([
      [5n, 1n],
      [17n, 1n],
    ]); // 5 * 17
    expect(Z_factor(115n)).toEqual([
      [5n, 1n],
      [23n, 1n],
    ]); // 5 * 23
    expect(Z_factor(187n)).toEqual([
      [11n, 1n],
      [17n, 1n],
    ]); // 11 * 17
    expect(Z_factor(1309n)).toEqual([
      [7n, 1n],
      [11n, 1n],
      [17n, 1n],
    ]); // 7 * 11 * 17
  });

  it('throws on zero', () => {
    expect(() => Z_factor(0n)).toThrow();
  });

  it('all factors are prime', () => {
    const testCases = [
      100n,
      1000n,
      10000n,
      123456n,
      654321n,
      999999n,
      49993895n, // The original failing case
    ];

    for (const n of testCases) {
      const factors = Z_factor(n);
      for (const [p, _e] of factors) {
        if (p !== -1n) {
          expect(isPrime(p)).toBe(true);
        }
      }
    }
  });

  it('factorization multiplies back to original', () => {
    const testCases = [1n, 2n, 6n, 12n, 100n, 1000n, 12345n, 49993895n, 123456789n];

    for (const n of testCases) {
      const factors = Z_factor(n);
      let product = 1n;
      for (const [p, e] of factors) {
        product *= p ** e;
      }
      expect(product).toBe(n);
    }
  });
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
