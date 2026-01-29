/**
 * Unit tests for IntegerLike and RationalLike type coercion
 */
import { describe, expect, test } from 'bun:test';
import { toBigInt, toRational, type IntegerLike, type RationalLike } from './coercion.js';
import { Integer } from '../rings/integer_ring.js';
import { Rational } from '../rings/rational.js';
import {
  gcd,
  lcm,
  xgcd,
  is_prime,
  factor,
  euler_phi,
  inverse_mod,
  power_mod,
  crt,
  isqrt,
  divisors,
  moebius,
  binomial,
  factorial,
  continued_fraction,
  rational_reconstruction,
} from '../arith/misc.js';

describe('toBigInt', () => {
  test('converts bigint to bigint', () => {
    expect(toBigInt(42n)).toBe(42n);
    expect(toBigInt(-100n)).toBe(-100n);
    expect(toBigInt(0n)).toBe(0n);
  });

  test('converts number to bigint', () => {
    expect(toBigInt(42)).toBe(42n);
    expect(toBigInt(-100)).toBe(-100n);
    expect(toBigInt(0)).toBe(0n);
  });

  test('converts Integer to bigint', () => {
    expect(toBigInt(new Integer(42n))).toBe(42n);
    expect(toBigInt(new Integer(-100))).toBe(-100n);
    expect(toBigInt(new Integer(0))).toBe(0n);
  });

  test('throws for non-integer numbers', () => {
    expect(() => toBigInt(3.14)).toThrow('cannot convert non-integer to Integer');
    expect(() => toBigInt(0.5)).toThrow('cannot convert non-integer to Integer');
  });

  test('throws RangeError for numbers exceeding safe integer range', () => {
    // This number loses precision when represented as a JavaScript number
    const unsafeNumber = Number.MAX_SAFE_INTEGER + 1;
    expect(() => toBigInt(unsafeNumber)).toThrow(RangeError);
    expect(() => toBigInt(unsafeNumber)).toThrow('exceeds safe integer range');

    // Negative side
    expect(() => toBigInt(Number.MIN_SAFE_INTEGER - 1)).toThrow(RangeError);
  });

  test('accepts numbers at the safe integer boundary', () => {
    expect(toBigInt(Number.MAX_SAFE_INTEGER)).toBe(BigInt(Number.MAX_SAFE_INTEGER));
    expect(toBigInt(Number.MIN_SAFE_INTEGER)).toBe(BigInt(Number.MIN_SAFE_INTEGER));
  });
});

describe('IntegerLike acceptance in arith functions', () => {
  describe('gcd accepts IntegerLike', () => {
    test('accepts bigint', () => {
      expect(gcd(12n, 8n)).toBe(4n);
    });

    test('accepts number', () => {
      expect(gcd(12, 8)).toBe(4n);
    });

    test('accepts Integer', () => {
      expect(gcd(new Integer(12), new Integer(8))).toBe(4n);
    });

    test('accepts mixed types', () => {
      expect(gcd(12n, 8)).toBe(4n);
      expect(gcd(12, new Integer(8))).toBe(4n);
      expect(gcd(new Integer(12), 8n)).toBe(4n);
    });

    test('accepts array of IntegerLike', () => {
      expect(gcd([12n, 8, new Integer(4)])).toBe(4n);
    });
  });

  describe('lcm accepts IntegerLike', () => {
    test('accepts bigint', () => {
      expect(lcm(4n, 6n)).toBe(12n);
    });

    test('accepts number', () => {
      expect(lcm(4, 6)).toBe(12n);
    });

    test('accepts Integer', () => {
      expect(lcm(new Integer(4), new Integer(6))).toBe(12n);
    });

    test('accepts mixed types', () => {
      expect(lcm(4n, 6)).toBe(12n);
      expect(lcm(4, new Integer(6))).toBe(12n);
    });
  });

  describe('xgcd accepts IntegerLike', () => {
    test('accepts bigint', () => {
      const [g, s, t] = xgcd(6n, 4n);
      expect(g).toBe(2n);
    });

    test('accepts number', () => {
      const [g, s, t] = xgcd(6, 4);
      expect(g).toBe(2n);
    });

    test('accepts Integer', () => {
      const [g, s, t] = xgcd(new Integer(6), new Integer(4));
      expect(g).toBe(2n);
    });
  });

  describe('is_prime accepts IntegerLike', () => {
    test('accepts bigint', () => {
      expect(is_prime(17n)).toBe(true);
      expect(is_prime(18n)).toBe(false);
    });

    test('accepts number', () => {
      expect(is_prime(17)).toBe(true);
      expect(is_prime(18)).toBe(false);
    });

    test('accepts Integer', () => {
      expect(is_prime(new Integer(17))).toBe(true);
      expect(is_prime(new Integer(18))).toBe(false);
    });
  });

  describe('factor accepts IntegerLike', () => {
    test('accepts bigint', () => {
      expect(factor(12n)).toEqual([[2n, 2n], [3n, 1n]]);
    });

    test('accepts number', () => {
      expect(factor(12)).toEqual([[2n, 2n], [3n, 1n]]);
    });

    test('accepts Integer', () => {
      expect(factor(new Integer(12))).toEqual([[2n, 2n], [3n, 1n]]);
    });
  });

  describe('euler_phi accepts IntegerLike', () => {
    test('accepts bigint', () => {
      expect(euler_phi(12n)).toBe(4n);
    });

    test('accepts number', () => {
      expect(euler_phi(12)).toBe(4n);
    });

    test('accepts Integer', () => {
      expect(euler_phi(new Integer(12))).toBe(4n);
    });
  });

  describe('inverse_mod accepts IntegerLike', () => {
    test('accepts bigint', () => {
      expect(inverse_mod(3n, 7n)).toBe(5n);
    });

    test('accepts number', () => {
      expect(inverse_mod(3, 7)).toBe(5n);
    });

    test('accepts Integer', () => {
      expect(inverse_mod(new Integer(3), new Integer(7))).toBe(5n);
    });

    test('accepts mixed types', () => {
      expect(inverse_mod(3n, 7)).toBe(5n);
      expect(inverse_mod(3, new Integer(7))).toBe(5n);
    });
  });

  describe('power_mod accepts IntegerLike', () => {
    test('accepts bigint', () => {
      expect(power_mod(2n, 10n, 1000n)).toBe(24n);
    });

    test('accepts number', () => {
      expect(power_mod(2, 10, 1000)).toBe(24n);
    });

    test('accepts Integer', () => {
      expect(power_mod(new Integer(2), new Integer(10), new Integer(1000))).toBe(24n);
    });

    test('accepts mixed types', () => {
      expect(power_mod(2n, 10, new Integer(1000))).toBe(24n);
    });
  });

  describe('crt accepts IntegerLike', () => {
    test('accepts bigint', () => {
      expect(crt(2n, 3n, 3n, 5n)).toBe(8n);
    });

    test('accepts number', () => {
      expect(crt(2, 3, 3, 5)).toBe(8n);
    });

    test('accepts Integer', () => {
      expect(crt(new Integer(2), new Integer(3), new Integer(3), new Integer(5))).toBe(8n);
    });

    test('accepts mixed types', () => {
      expect(crt(2n, 3, new Integer(3), 5n)).toBe(8n);
    });
  });

  describe('isqrt accepts IntegerLike', () => {
    test('accepts bigint', () => {
      expect(isqrt(16n)).toBe(4n);
      expect(isqrt(17n)).toBe(4n);
    });

    test('accepts number', () => {
      expect(isqrt(16)).toBe(4n);
      expect(isqrt(17)).toBe(4n);
    });

    test('accepts Integer', () => {
      expect(isqrt(new Integer(16))).toBe(4n);
      expect(isqrt(new Integer(17))).toBe(4n);
    });
  });

  describe('divisors accepts IntegerLike', () => {
    test('accepts bigint', () => {
      expect(divisors(12n)).toEqual([1n, 2n, 3n, 4n, 6n, 12n]);
    });

    test('accepts number', () => {
      expect(divisors(12)).toEqual([1n, 2n, 3n, 4n, 6n, 12n]);
    });

    test('accepts Integer', () => {
      expect(divisors(new Integer(12))).toEqual([1n, 2n, 3n, 4n, 6n, 12n]);
    });
  });

  describe('moebius accepts IntegerLike', () => {
    test('accepts bigint', () => {
      expect(moebius(6n)).toBe(1n);
      expect(moebius(4n)).toBe(0n);
    });

    test('accepts number', () => {
      expect(moebius(6)).toBe(1n);
      expect(moebius(4)).toBe(0n);
    });

    test('accepts Integer', () => {
      expect(moebius(new Integer(6))).toBe(1n);
      expect(moebius(new Integer(4))).toBe(0n);
    });
  });

  describe('binomial accepts IntegerLike', () => {
    test('accepts bigint', () => {
      expect(binomial(5n, 2n)).toBe(10n);
    });

    test('accepts number', () => {
      expect(binomial(5, 2)).toBe(10n);
    });

    test('accepts Integer', () => {
      expect(binomial(new Integer(5), new Integer(2))).toBe(10n);
    });

    test('accepts mixed types', () => {
      expect(binomial(5n, 2)).toBe(10n);
      expect(binomial(5, new Integer(2))).toBe(10n);
    });
  });

  describe('factorial accepts IntegerLike', () => {
    test('accepts bigint', () => {
      expect(factorial(5n)).toBe(120n);
    });

    test('accepts number', () => {
      expect(factorial(5)).toBe(120n);
    });

    test('accepts Integer', () => {
      expect(factorial(new Integer(5))).toBe(120n);
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

  test('converts number to Rational', () => {
    const r = toRational(42);
    expect(r.numerator).toBe(42n);
    expect(r.denominator).toBe(1n);
  });

  test('converts negative number to Rational', () => {
    const r = toRational(-100);
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

  test('throws for non-integer numbers', () => {
    expect(() => toRational(3.14)).toThrow('cannot convert non-integer to Integer');
    expect(() => toRational(0.5)).toThrow('cannot convert non-integer to Integer');
  });
});

describe('RationalLike acceptance in arith functions', () => {
  describe('continued_fraction accepts RationalLike', () => {
    test('accepts bigint numerator and denominator', () => {
      expect(continued_fraction(13n, 9n)).toEqual([1n, 2n, 4n]);
    });

    test('accepts number numerator and denominator', () => {
      expect(continued_fraction(13, 9)).toEqual([1n, 2n, 4n]);
    });

    test('accepts Integer numerator and denominator', () => {
      expect(continued_fraction(new Integer(13), new Integer(9))).toEqual([1n, 2n, 4n]);
    });

    test('accepts Rational directly', () => {
      const r = new Rational(13n, 9n);
      expect(continued_fraction(r)).toEqual([1n, 2n, 4n]);
    });

    test('accepts mixed types', () => {
      expect(continued_fraction(13n, 9)).toEqual([1n, 2n, 4n]);
      expect(continued_fraction(13, new Integer(9))).toEqual([1n, 2n, 4n]);
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

    test('accepts number', () => {
      const [p, q] = rational_reconstruction(11323, 100000);
      expect(p).toBe(119n);
      expect(q).toBe(53n);
    });

    test('accepts Integer', () => {
      const [p, q] = rational_reconstruction(new Integer(11323), new Integer(100000));
      expect(p).toBe(119n);
      expect(q).toBe(53n);
    });

    test('accepts mixed types', () => {
      const [p, q] = rational_reconstruction(11323n, 100000);
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
