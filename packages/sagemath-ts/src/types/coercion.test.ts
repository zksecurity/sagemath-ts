/**
 * Unit tests for IntegerLike and RationalLike type coercion
 *
 * @see Deviation: Language and Type-System Adaptations
 */
import { describe, expect, test } from 'bun:test';
import {
  binomial,
  continued_fraction,
  crt,
  divisors,
  euler_phi,
  factor,
  factorial,
  gcd,
  inverse_mod,
  is_prime,
  isqrt,
  lcm,
  moebius,
  power_mod,
  rational_reconstruction,
  xgcd,
} from '../arith/misc.js';
import { Integer } from '../rings/integer_ring.js';
import { Rational } from '../rings/rational.js';
import { type IntegerLike, type RationalLike, toBigInt, toRational, toSafeNumber } from './coercion.js';

describe('toBigInt', () => {
  test('converts bigint to bigint', () => {
    expect(toBigInt(42n)).toBe(42n);
    expect(toBigInt(-100n)).toBe(-100n);
    expect(toBigInt(0n)).toBe(0n);
  });

  test('rejects JavaScript numbers (precision loss risk)', () => {
    // Numbers are rejected to prevent silent precision loss
    // @ts-expect-error - intentionally testing runtime behavior
    expect(() => toBigInt(42)).toThrow(TypeError);
    // @ts-expect-error - intentionally testing runtime behavior
    expect(() => toBigInt(42)).toThrow('JavaScript numbers are not accepted');
  });

  test('converts Integer to bigint', () => {
    expect(toBigInt(new Integer(42n))).toBe(42n);
    expect(toBigInt(new Integer(-100n))).toBe(-100n);
    expect(toBigInt(new Integer(0n))).toBe(0n);
  });
});

describe('toSafeNumber', () => {
  test('converts small bigints to number', () => {
    expect(toSafeNumber(42n)).toBe(42);
    expect(toSafeNumber(-100n)).toBe(-100);
    expect(toSafeNumber(0n)).toBe(0);
  });

  test('converts bigints at safe integer boundary', () => {
    expect(toSafeNumber(BigInt(Number.MAX_SAFE_INTEGER))).toBe(Number.MAX_SAFE_INTEGER);
    expect(toSafeNumber(BigInt(Number.MIN_SAFE_INTEGER))).toBe(Number.MIN_SAFE_INTEGER);
  });

  test('throws RangeError for bigints exceeding safe integer range', () => {
    const unsafeBigint = BigInt(Number.MAX_SAFE_INTEGER) + 1n;
    expect(() => toSafeNumber(unsafeBigint)).toThrow(RangeError);
    expect(() => toSafeNumber(unsafeBigint)).toThrow('exceeds safe integer range');

    // Negative side
    const negativeUnsafe = BigInt(Number.MIN_SAFE_INTEGER) - 1n;
    expect(() => toSafeNumber(negativeUnsafe)).toThrow(RangeError);
  });

  test('throws for very large bigints', () => {
    expect(() => toSafeNumber(10n ** 100n)).toThrow(RangeError);
    expect(() => toSafeNumber(-(10n ** 100n))).toThrow(RangeError);
  });
});

describe('IntegerLike acceptance in arith functions', () => {
  describe('gcd accepts IntegerLike', () => {
    test('accepts bigint', () => {
      expect(gcd(12n, 8n)).toBe(4n);
    });

    test('accepts Integer', () => {
      expect(gcd(new Integer(12n), new Integer(8n))).toBe(4n);
    });

    test('accepts mixed bigint and Integer', () => {
      expect(gcd(12n, new Integer(8n))).toBe(4n);
      expect(gcd(new Integer(12n), 8n)).toBe(4n);
    });

    test('accepts array of IntegerLike', () => {
      expect(gcd([12n, 8n, new Integer(4n)])).toBe(4n);
    });
  });

  describe('lcm accepts IntegerLike', () => {
    test('accepts bigint', () => {
      expect(lcm(4n, 6n)).toBe(12n);
    });

    test('accepts Integer', () => {
      expect(lcm(new Integer(4n), new Integer(6n))).toBe(12n);
    });

    test('accepts mixed bigint and Integer', () => {
      expect(lcm(4n, new Integer(6n))).toBe(12n);
    });
  });

  describe('xgcd accepts IntegerLike', () => {
    test('accepts bigint', () => {
      const [g, _s, _t] = xgcd(6n, 4n);
      expect(g).toBe(2n);
    });

    test('accepts Integer', () => {
      const [g, _s, _t] = xgcd(new Integer(6n), new Integer(4n));
      expect(g).toBe(2n);
    });
  });

  describe('is_prime accepts IntegerLike', () => {
    test('accepts bigint', () => {
      expect(is_prime(17n)).toBe(true);
      expect(is_prime(18n)).toBe(false);
    });

    test('accepts Integer', () => {
      expect(is_prime(new Integer(17n))).toBe(true);
      expect(is_prime(new Integer(18n))).toBe(false);
    });
  });

  describe('factor accepts IntegerLike', () => {
    test('accepts bigint', () => {
      expect(factor(12n)).toEqual([
        [2n, 2n],
        [3n, 1n],
      ]);
    });

    test('accepts Integer', () => {
      expect(factor(new Integer(12n))).toEqual([
        [2n, 2n],
        [3n, 1n],
      ]);
    });
  });

  describe('euler_phi accepts IntegerLike', () => {
    test('accepts bigint', () => {
      expect(euler_phi(12n)).toBe(4n);
    });

    test('accepts Integer', () => {
      expect(euler_phi(new Integer(12n))).toBe(4n);
    });
  });

  describe('inverse_mod accepts IntegerLike', () => {
    test('accepts bigint', () => {
      expect(inverse_mod(3n, 7n)).toBe(5n);
    });

    test('accepts Integer', () => {
      expect(inverse_mod(new Integer(3n), new Integer(7n))).toBe(5n);
    });

    test('accepts mixed bigint and Integer', () => {
      expect(inverse_mod(3n, new Integer(7n))).toBe(5n);
    });
  });

  describe('power_mod accepts IntegerLike', () => {
    test('accepts bigint', () => {
      expect(power_mod(2n, 10n, 1000n)).toBe(24n);
    });

    test('accepts Integer', () => {
      expect(power_mod(new Integer(2n), new Integer(10n), new Integer(1000n))).toBe(24n);
    });

    test('accepts mixed bigint and Integer', () => {
      expect(power_mod(2n, new Integer(10n), 1000n)).toBe(24n);
    });
  });

  describe('crt accepts IntegerLike', () => {
    test('accepts bigint', () => {
      expect(crt(2n, 3n, 3n, 5n)).toBe(8n);
    });

    test('accepts Integer', () => {
      expect(crt(new Integer(2n), new Integer(3n), new Integer(3n), new Integer(5n))).toBe(8n);
    });

    test('accepts mixed bigint and Integer', () => {
      expect(crt(2n, new Integer(3n), 3n, new Integer(5n))).toBe(8n);
    });
  });

  describe('isqrt accepts IntegerLike', () => {
    test('accepts bigint', () => {
      expect(isqrt(16n)).toBe(4n);
      expect(isqrt(17n)).toBe(4n);
    });

    test('accepts Integer', () => {
      expect(isqrt(new Integer(16n))).toBe(4n);
      expect(isqrt(new Integer(17n))).toBe(4n);
    });
  });

  describe('divisors accepts IntegerLike', () => {
    test('accepts bigint', () => {
      expect(divisors(12n)).toEqual([1n, 2n, 3n, 4n, 6n, 12n]);
    });

    test('accepts Integer', () => {
      expect(divisors(new Integer(12n))).toEqual([1n, 2n, 3n, 4n, 6n, 12n]);
    });
  });

  describe('moebius accepts IntegerLike', () => {
    test('accepts bigint', () => {
      expect(moebius(6n)).toBe(1n);
      expect(moebius(4n)).toBe(0n);
    });

    test('accepts Integer', () => {
      expect(moebius(new Integer(6n))).toBe(1n);
      expect(moebius(new Integer(4n))).toBe(0n);
    });
  });

  describe('binomial accepts IntegerLike', () => {
    test('accepts bigint', () => {
      expect(binomial(5n, 2n)).toBe(10n);
    });

    test('accepts Integer', () => {
      expect(binomial(new Integer(5n), new Integer(2n))).toBe(10n);
    });

    test('accepts mixed bigint and Integer', () => {
      expect(binomial(5n, new Integer(2n))).toBe(10n);
    });
  });

  describe('factorial accepts IntegerLike', () => {
    test('accepts bigint', () => {
      expect(factorial(5n)).toBe(120n);
    });

    test('accepts Integer', () => {
      expect(factorial(new Integer(5n))).toBe(120n);
    });
  });
});

describe('toRational', () => {
  test('returns Rational unchanged', () => {
    const r = new Rational(3n, 4n);
    expect(toRational(r)).toBe(r);
  });

  test('converts bigint to Rational with denominator 1', () => {
    const r = toRational(42n);
    expect(r.numerator).toBe(42n);
    expect(r.denominator).toBe(1n);
  });

  test('converts negative bigint to Rational', () => {
    const r = toRational(-100n);
    expect(r.numerator).toBe(-100n);
    expect(r.denominator).toBe(1n);
  });

  test('converts Integer to Rational', () => {
    const r = toRational(new Integer(42n));
    expect(r.numerator).toBe(42n);
    expect(r.denominator).toBe(1n);
  });

  test('converts zero to Rational', () => {
    const r = toRational(0n);
    expect(r.numerator).toBe(0n);
    expect(r.denominator).toBe(1n);
  });

  test('rejects JavaScript numbers', () => {
    // @ts-expect-error - intentionally testing runtime behavior
    expect(() => toRational(42)).toThrow(TypeError);
    // @ts-expect-error - intentionally testing runtime behavior
    expect(() => toRational(3.14)).toThrow('JavaScript numbers are not accepted');
  });
});

describe('RationalLike acceptance in arith functions', () => {
  describe('continued_fraction accepts RationalLike', () => {
    test('accepts bigint numerator and denominator', () => {
      expect(continued_fraction(13n, 9n)).toEqual([1n, 2n, 4n]);
    });

    test('accepts Integer numerator and denominator', () => {
      expect(continued_fraction(new Integer(13n), new Integer(9n))).toEqual([1n, 2n, 4n]);
    });

    test('accepts Rational directly', () => {
      const r = new Rational(13n, 9n);
      expect(continued_fraction(r)).toEqual([1n, 2n, 4n]);
    });

    test('accepts mixed bigint and Integer', () => {
      expect(continued_fraction(13n, new Integer(9n))).toEqual([1n, 2n, 4n]);
    });

    test('works with negative rationals', () => {
      expect(continued_fraction(-1n, 3n)).toEqual([-1n, 1n, 2n]);
      expect(continued_fraction(new Rational(-1n, 3n))).toEqual([-1n, 1n, 2n]);
    });

    test('works with integer values (denominator 1)', () => {
      expect(continued_fraction(5n)).toEqual([5n]);
      expect(continued_fraction(new Rational(5n, 1n))).toEqual([5n]);
    });
  });

  describe('rational_reconstruction accepts IntegerLike', () => {
    test('accepts bigint', () => {
      const [p, q] = rational_reconstruction(11323n, 100000n);
      expect(p).toBe(119n);
      expect(q).toBe(53n);
    });

    test('accepts Integer', () => {
      const [p, q] = rational_reconstruction(new Integer(11323n), new Integer(100000n));
      expect(p).toBe(119n);
      expect(q).toBe(53n);
    });

    test('accepts mixed bigint and Integer', () => {
      const [p, q] = rational_reconstruction(11323n, new Integer(100000n));
      expect(p).toBe(119n);
      expect(q).toBe(53n);
    });

    test('handles small values', () => {
      const [p, q] = rational_reconstruction(3n, 292393n);
      expect(p).toBe(3n);
      expect(q).toBe(1n);
    });
  });
});
