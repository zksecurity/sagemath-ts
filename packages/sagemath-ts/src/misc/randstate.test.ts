/**
 * Unit tests for sage/misc/randstate
 *
 * Reference: reference/sage/src/sage/misc/randstate.pyx
 */
import { describe, expect, test } from 'bun:test';
import { ValueError } from '../errors.js';
import {
  RandState,
  SAGE_RAND_MAX,
  current_randstate,
  initial_seed,
  random,
  set_random_seed,
} from './randstate.js';

/**
 * Largest run of a repeating period-`p` pattern at the start of `seq`.
 */
function hasPeriod(seq: string[], p: number): boolean {
  if (seq.length <= p) return false;
  for (let i = p; i < seq.length; i++) {
    if (seq[i] !== seq[i - p]) return false;
  }
  return true;
}

describe('RandState', () => {
  describe('seed handling', () => {
    test('seed() returns the initial seed, not the evolving state', () => {
      // randstate.pyx:562-580
      const r = new RandState(314159);
      expect(r.seed()).toBe(314159n);
      r.c_rand_double();
      r.random_bits(64);
      r.randint(0n, 1000n);
      expect(r.seed()).toBe(314159n);
    });

    test('set_seed updates the reported seed and restarts the stream', () => {
      const r = new RandState(0);
      r.set_seed(12345);
      expect(r.seed()).toBe(12345n);
      const first = [r.random_bits(32), r.random_bits(32), r.random_bits(32)];
      r.set_seed(12345);
      const second = [r.random_bits(32), r.random_bits(32), r.random_bits(32)];
      expect(second).toEqual(first);
    });

    test('the sign of the seed is ignored, as in GMP mpz_export', () => {
      // randstate.pyx:955-969 shows set_random_seed(12345) and
      // set_random_seed(-12345) producing the same value.
      const a = new RandState(12345);
      const b = new RandState(-12345);
      expect(b.random_bits(64)).toBe(a.random_bits(64));
    });

    test('accepts arbitrarily large bigint seeds', () => {
      const big = 305866218880103397618377824640007711767n;
      const r = new RandState(big);
      expect(r.seed()).toBe(big);
      const s = new RandState(big);
      expect(s.random_bits(64)).toBe(r.random_bits(64));
    });

    test('different seeds give different streams', () => {
      const a = new RandState(1).random_bits(64);
      const b = new RandState(2).random_bits(64);
      expect(a).not.toBe(b);
    });
  });

  describe('no short cycles (C1 regression)', () => {
    test('random_below(2) does not alternate', () => {
      const r = new RandState(7);
      const bits: string[] = [];
      for (let i = 0; i < 200; i++) bits.push(r.random_below(2n).toString());
      expect(hasPeriod(bits, 2)).toBe(false);
      expect(hasPeriod(bits, 1)).toBe(false);
      const ones = bits.filter((b) => b === '1').length;
      expect(ones).toBeGreaterThan(70);
      expect(ones).toBeLessThan(130);
    });

    test('random_below(4) has no period-4 cycle', () => {
      const r = new RandState(7);
      const vals: string[] = [];
      for (let i = 0; i < 200; i++) vals.push(r.random_below(4n).toString());
      for (let p = 1; p <= 8; p++) {
        expect(hasPeriod(vals, p)).toBe(false);
      }
    });

    test('randint(0, 9) has no period-10 cycle', () => {
      const r = new RandState(7);
      const vals: string[] = [];
      for (let i = 0; i < 500; i++) vals.push(r.randint(0n, 9n).toString());
      expect(hasPeriod(vals, 10)).toBe(false);
    });

    test('every low bit is balanced', () => {
      // For an LCG mod 2^64, bit k has period 2^(k+1); the low bits of the
      // Mersenne Twister are not structured this way.
      const r = new RandState(99);
      const n = 20000;
      const ones = new Array(8).fill(0);
      for (let i = 0; i < n; i++) {
        const v = r.random_bits(64);
        for (let k = 0; k < 8; k++) {
          if ((v >> BigInt(k)) & 1n) ones[k]!++;
        }
      }
      for (let k = 0; k < 8; k++) {
        expect(Math.abs(ones[k]! / n - 0.5)).toBeLessThan(0.02);
      }
    });

    test('random_below(10) is uniform (chi-squared)', () => {
      const r = new RandState(12345);
      const n = 200000;
      const counts = new Array(10).fill(0);
      for (let i = 0; i < n; i++) counts[Number(r.random_below(10n))]!++;
      const expected = n / 10;
      let chi2 = 0;
      for (const c of counts) chi2 += ((c - expected) * (c - expected)) / expected;
      // 99.9% quantile of chi^2 with 9 degrees of freedom is 27.88.
      expect(chi2).toBeLessThan(27.88);
    });
  });

  describe('random_bits', () => {
    test('returns values in [0, 2^bits)', () => {
      const r = new RandState(3);
      for (const bits of [1, 7, 31, 32, 33, 64, 65, 200]) {
        for (let i = 0; i < 50; i++) {
          const v = r.random_bits(bits);
          expect(v >= 0n).toBe(true);
          expect(v < 1n << BigInt(bits)).toBe(true);
        }
      }
    });

    test('random_bits(0) is 0', () => {
      expect(new RandState(3).random_bits(0)).toBe(0n);
    });

    test('rejects a negative bit count', () => {
      expect(() => new RandState(3).random_bits(-1)).toThrow(ValueError);
    });

    test('large requests use the full width', () => {
      const r = new RandState(5);
      let maxSeen = 0n;
      for (let i = 0; i < 200; i++) {
        const v = r.random_bits(128);
        if (v > maxSeen) maxSeen = v;
      }
      expect(maxSeen > 1n << 120n).toBe(true);
    });
  });

  describe('c_random and c_rand_double', () => {
    test('c_random returns a 31-bit integer', () => {
      const r = new RandState(1207);
      for (let i = 0; i < 1000; i++) {
        const v = r.c_random();
        expect(Number.isInteger(v)).toBe(true);
        expect(v).toBeGreaterThanOrEqual(0);
        expect(v).toBeLessThanOrEqual(SAGE_RAND_MAX);
      }
    });

    test('c_rand_double lies in [0, 1) and consumes two draws', () => {
      // randstate.pyx:883-895
      const r = new RandState(2718281828);
      for (let i = 0; i < 1000; i++) {
        const d = r.c_rand_double();
        expect(d).toBeGreaterThanOrEqual(0);
        expect(d).toBeLessThan(1);
      }

      const a = new RandState(11);
      const b = new RandState(11);
      a.c_rand_double();
      b.random_bits(25);
      b.random_bits(28);
      // Both states consumed the same amount of randomness.
      expect(a.random_bits(32)).toBe(b.random_bits(32));
    });

    test('c_rand_double is uniform on [0, 1)', () => {
      const r = new RandState(4242);
      const n = 100000;
      const buckets = new Array(10).fill(0);
      let sum = 0;
      for (let i = 0; i < n; i++) {
        const d = r.c_rand_double();
        sum += d;
        buckets[Math.min(9, Math.floor(d * 10))]!++;
      }
      expect(Math.abs(sum / n - 0.5)).toBeLessThan(0.01);
      for (const b of buckets) {
        expect(Math.abs(b / n - 0.1)).toBeLessThan(0.01);
      }
    });

    test('random() is an alias for c_rand_double', () => {
      const a = new RandState(77);
      const b = new RandState(77);
      expect(a.random()).toBe(b.c_rand_double());
    });
  });

  describe('random_below and randint', () => {
    test('random_below rejects a non-positive bound', () => {
      expect(() => new RandState(1).random_below(0n)).toThrow(ValueError);
      expect(() => new RandState(1).random_below(-5n)).toThrow(ValueError);
    });

    test('randint covers the whole inclusive range', () => {
      const r = new RandState(8);
      const seen = new Set<string>();
      for (let i = 0; i < 2000; i++) {
        const v = r.randint(-3n, 4n);
        expect(v >= -3n).toBe(true);
        expect(v <= 4n).toBe(true);
        seen.add(v.toString());
      }
      expect(seen.size).toBe(8);
    });

    test('randint rejects min > max', () => {
      expect(() => new RandState(1).randint(5n, 4n)).toThrow(ValueError);
    });
  });
});

describe('module-level functions', () => {
  test('set_random_seed makes the global stream reproducible', () => {
    set_random_seed(0);
    const first = [random(), random(), random()];
    set_random_seed(0);
    const second = [random(), random(), random()];
    expect(second).toEqual(first);
  });

  test('random() returns a 31-bit integer', () => {
    // randstate.pyx:1008-1020
    set_random_seed(31);
    for (let i = 0; i < 100; i++) {
      const v = random();
      expect(typeof v).toBe('bigint');
      expect(v >= 0n).toBe(true);
      expect(v < 1n << 31n).toBe(true);
    }
  });

  test('initial_seed returns the seed of the current randstate', () => {
    // randstate.pyx:1024-1042
    set_random_seed(42);
    expect(initial_seed()).toBe(42n);
    random();
    current_randstate().random_bits(64);
    expect(initial_seed()).toBe(42n);

    set_random_seed(5);
    expect(initial_seed()).toBe(5n);
  });
});
