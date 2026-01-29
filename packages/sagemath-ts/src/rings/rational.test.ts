/**
 * Tests for rational numbers (QQ)
 */

import { describe, expect, it } from 'vitest';
import { Rational } from './rational.js';
import { QQ, RationalField } from './rational_field.js';

describe('Rational', () => {
  describe('construction', () => {
    it('should create from two bigints', () => {
      const r = new Rational(3n, 4n);
      expect(r.numerator).toBe(3n);
      expect(r.denominator).toBe(4n);
    });

    it('should create from single bigint', () => {
      const r = new Rational(7n);
      expect(r.numerator).toBe(7n);
      expect(r.denominator).toBe(1n);
    });

    it('should throw on zero denominator', () => {
      expect(() => new Rational(1n, 0n)).toThrow('denominator must not be 0');
    });
  });

  describe('reduction to lowest terms', () => {
    it('should reduce 6/4 to 3/2', () => {
      const r = new Rational(6n, 4n);
      expect(r.numerator).toBe(3n);
      expect(r.denominator).toBe(2n);
    });

    it('should reduce -6/4 to -3/2', () => {
      const r = new Rational(-6n, 4n);
      expect(r.numerator).toBe(-3n);
      expect(r.denominator).toBe(2n);
    });

    it('should reduce 12/8 to 3/2', () => {
      const r = new Rational(12n, 8n);
      expect(r.numerator).toBe(3n);
      expect(r.denominator).toBe(2n);
    });

    it('should reduce 100/25 to 4/1', () => {
      const r = new Rational(100n, 25n);
      expect(r.numerator).toBe(4n);
      expect(r.denominator).toBe(1n);
    });

    it('should reduce gcd(12, 18) case', () => {
      const r = new Rational(12n, 18n);
      expect(r.numerator).toBe(2n);
      expect(r.denominator).toBe(3n);
    });
  });

  describe('denominator normalization', () => {
    it('should keep denominator positive when both negative', () => {
      const r = new Rational(-3n, -4n);
      expect(r.numerator).toBe(3n);
      expect(r.denominator).toBe(4n);
    });

    it('should move sign to numerator', () => {
      const r = new Rational(3n, -4n);
      expect(r.numerator).toBe(-3n);
      expect(r.denominator).toBe(4n);
    });
  });

  describe('zero handling', () => {
    it('should normalize 0/n to 0/1', () => {
      const r = new Rational(0n, 5n);
      expect(r.numerator).toBe(0n);
      expect(r.denominator).toBe(1n);
    });

    it('should normalize 0/-5 to 0/1', () => {
      const r = new Rational(0n, -5n);
      expect(r.numerator).toBe(0n);
      expect(r.denominator).toBe(1n);
    });
  });

  describe('properties', () => {
    it('should return correct sign for positive', () => {
      expect(new Rational(3n, 4n).sign).toBe(1n);
    });

    it('should return correct sign for negative', () => {
      expect(new Rational(-3n, 4n).sign).toBe(-1n);
    });

    it('should return correct sign for zero', () => {
      expect(new Rational(0n, 1n).sign).toBe(0n);
    });

    it('should have numer alias', () => {
      const r = new Rational(5n, 7n);
      expect(r.numer).toBe(5n);
    });

    it('should have denom alias', () => {
      const r = new Rational(5n, 7n);
      expect(r.denom).toBe(7n);
    });
  });

  describe('arithmetic: add', () => {
    it('should add two rationals', () => {
      const a = new Rational(1n, 3n);
      const b = new Rational(1n, 6n);
      const c = a.add(b);
      expect(c.numerator).toBe(1n);
      expect(c.denominator).toBe(2n);
    });

    it('should add with different denominators', () => {
      const a = new Rational(2n, 3n);
      const b = new Rational(1n, 4n);
      const c = a.add(b);
      expect(c.numerator).toBe(11n);
      expect(c.denominator).toBe(12n);
    });

    it('should add with bigint', () => {
      const a = new Rational(1n, 2n);
      const c = a.add(3n);
      expect(c.numerator).toBe(7n);
      expect(c.denominator).toBe(2n);
    });

    it('should add negative numbers', () => {
      const a = new Rational(-1n, 2n);
      const b = new Rational(1n, 3n);
      const c = a.add(b);
      expect(c.numerator).toBe(-1n);
      expect(c.denominator).toBe(6n);
    });
  });

  describe('arithmetic: sub', () => {
    it('should subtract two rationals', () => {
      const a = new Rational(1n, 2n);
      const b = new Rational(1n, 3n);
      const c = a.sub(b);
      expect(c.numerator).toBe(1n);
      expect(c.denominator).toBe(6n);
    });

    it('should subtract with bigint', () => {
      const a = new Rational(7n, 2n);
      const c = a.sub(3n);
      expect(c.numerator).toBe(1n);
      expect(c.denominator).toBe(2n);
    });
  });

  describe('arithmetic: mul', () => {
    it('should multiply two rationals', () => {
      const a = new Rational(2n, 3n);
      const b = new Rational(3n, 4n);
      const c = a.mul(b);
      expect(c.numerator).toBe(1n);
      expect(c.denominator).toBe(2n);
    });

    it('should multiply with bigint', () => {
      const a = new Rational(2n, 3n);
      const c = a.mul(6n);
      expect(c.numerator).toBe(4n);
      expect(c.denominator).toBe(1n);
    });

    it('should multiply by zero', () => {
      const a = new Rational(2n, 3n);
      const b = Rational.zero();
      const c = a.mul(b);
      expect(c.isZero()).toBe(true);
    });
  });

  describe('arithmetic: div', () => {
    it('should divide two rationals', () => {
      const a = new Rational(2n, 3n);
      const b = new Rational(4n, 5n);
      const c = a.div(b);
      expect(c.numerator).toBe(5n);
      expect(c.denominator).toBe(6n);
    });

    it('should divide with bigint', () => {
      const a = new Rational(2n, 3n);
      const c = a.div(2n);
      expect(c.numerator).toBe(1n);
      expect(c.denominator).toBe(3n);
    });

    it('should throw on division by zero', () => {
      const a = new Rational(2n, 3n);
      expect(() => a.div(Rational.zero())).toThrow('rational division by zero');
      expect(() => a.div(0n)).toThrow('rational division by zero');
    });
  });

  describe('arithmetic: neg', () => {
    it('should negate positive', () => {
      const r = new Rational(3n, 4n).neg();
      expect(r.numerator).toBe(-3n);
      expect(r.denominator).toBe(4n);
    });

    it('should negate negative', () => {
      const r = new Rational(-3n, 4n).neg();
      expect(r.numerator).toBe(3n);
      expect(r.denominator).toBe(4n);
    });

    it('should negate zero', () => {
      const r = Rational.zero().neg();
      expect(r.isZero()).toBe(true);
    });
  });

  describe('arithmetic: inv', () => {
    it('should invert positive rational', () => {
      const r = new Rational(3n, 4n).inv();
      expect(r.numerator).toBe(4n);
      expect(r.denominator).toBe(3n);
    });

    it('should invert negative rational', () => {
      const r = new Rational(-3n, 4n).inv();
      expect(r.numerator).toBe(-4n);
      expect(r.denominator).toBe(3n);
    });

    it('should throw on invert zero', () => {
      expect(() => Rational.zero().inv()).toThrow('rational division by zero');
    });
  });

  describe('arithmetic: pow', () => {
    it('should raise to positive power', () => {
      const r = new Rational(2n, 3n).pow(3n);
      expect(r.numerator).toBe(8n);
      expect(r.denominator).toBe(27n);
    });

    it('should raise to negative power', () => {
      const r = new Rational(2n, 3n).pow(-2n);
      expect(r.numerator).toBe(9n);
      expect(r.denominator).toBe(4n);
    });

    it('should raise to zero power', () => {
      const r = new Rational(2n, 3n).pow(0n);
      expect(r.numerator).toBe(1n);
      expect(r.denominator).toBe(1n);
    });

    it('should handle negative base', () => {
      const r = new Rational(-2n, 3n).pow(3n);
      expect(r.numerator).toBe(-8n);
      expect(r.denominator).toBe(27n);
    });

    it('should throw on 0^negative', () => {
      expect(() => Rational.zero().pow(-1n)).toThrow('rational division by zero');
    });
  });

  describe('arithmetic: abs', () => {
    it('should return absolute value of positive', () => {
      const r = new Rational(3n, 4n).abs();
      expect(r.numerator).toBe(3n);
    });

    it('should return absolute value of negative', () => {
      const r = new Rational(-3n, 4n).abs();
      expect(r.numerator).toBe(3n);
    });
  });

  describe('comparisons', () => {
    it('should test equality', () => {
      expect(new Rational(1n, 2n).eq(new Rational(2n, 4n))).toBe(true);
      expect(new Rational(1n, 2n).eq(new Rational(1n, 3n))).toBe(false);
    });

    it('should compare with bigint', () => {
      expect(new Rational(6n, 2n).eq(3n)).toBe(true);
      expect(new Rational(6n, 2n).eq(4n)).toBe(false);
    });

    it('should test less than', () => {
      expect(new Rational(1n, 3n).lt(new Rational(1n, 2n))).toBe(true);
      expect(new Rational(1n, 2n).lt(new Rational(1n, 3n))).toBe(false);
    });

    it('should test less than or equal', () => {
      expect(new Rational(1n, 2n).le(new Rational(1n, 2n))).toBe(true);
      expect(new Rational(1n, 3n).le(new Rational(1n, 2n))).toBe(true);
      expect(new Rational(2n, 3n).le(new Rational(1n, 2n))).toBe(false);
    });

    it('should test greater than', () => {
      expect(new Rational(2n, 3n).gt(new Rational(1n, 2n))).toBe(true);
      expect(new Rational(1n, 3n).gt(new Rational(1n, 2n))).toBe(false);
    });

    it('should test greater than or equal', () => {
      expect(new Rational(1n, 2n).ge(new Rational(1n, 2n))).toBe(true);
      expect(new Rational(2n, 3n).ge(new Rational(1n, 2n))).toBe(true);
      expect(new Rational(1n, 3n).ge(new Rational(1n, 2n))).toBe(false);
    });

    it('should compare with cmp', () => {
      expect(new Rational(1n, 3n).cmp(new Rational(1n, 2n))).toBe(-1);
      expect(new Rational(1n, 2n).cmp(new Rational(1n, 2n))).toBe(0);
      expect(new Rational(2n, 3n).cmp(new Rational(1n, 2n))).toBe(1);
    });

    it('should handle negative comparisons', () => {
      expect(new Rational(-1n, 2n).lt(new Rational(1n, 3n))).toBe(true);
      expect(new Rational(-1n, 2n).lt(new Rational(-1n, 3n))).toBe(true);
    });
  });

  describe('conversions', () => {
    it('should convert to string with denominator 1', () => {
      expect(new Rational(5n, 1n).toString()).toBe('5');
    });

    it('should convert to string with denominator > 1', () => {
      expect(new Rational(3n, 4n).toString()).toBe('3/4');
    });

    it('should convert negative to string', () => {
      expect(new Rational(-3n, 4n).toString()).toBe('-3/4');
    });

    it('should convert to number', () => {
      expect(new Rational(1n, 2n).toNumber()).toBe(0.5);
      expect(new Rational(1n, 3n).toNumber()).toBeCloseTo(0.333333, 5);
    });
  });

  describe('floor, ceil, round', () => {
    it('should compute floor of positive', () => {
      expect(new Rational(7n, 3n).floor()).toBe(2n); // 2.333... -> 2
      expect(new Rational(7n, 2n).floor()).toBe(3n); // 3.5 -> 3
    });

    it('should compute floor of negative', () => {
      expect(new Rational(-7n, 3n).floor()).toBe(-3n); // -2.333... -> -3
      expect(new Rational(-7n, 2n).floor()).toBe(-4n); // -3.5 -> -4
    });

    it('should compute floor of integer', () => {
      expect(new Rational(5n, 1n).floor()).toBe(5n);
      expect(new Rational(-5n, 1n).floor()).toBe(-5n);
    });

    it('should compute ceil of positive', () => {
      expect(new Rational(7n, 3n).ceil()).toBe(3n); // 2.333... -> 3
      expect(new Rational(7n, 2n).ceil()).toBe(4n); // 3.5 -> 4
    });

    it('should compute ceil of negative', () => {
      expect(new Rational(-7n, 3n).ceil()).toBe(-2n); // -2.333... -> -2
      expect(new Rational(-7n, 2n).ceil()).toBe(-3n); // -3.5 -> -3
    });

    it('should compute ceil of integer', () => {
      expect(new Rational(5n, 1n).ceil()).toBe(5n);
      expect(new Rational(-5n, 1n).ceil()).toBe(-5n);
    });

    it('should round with default even mode', () => {
      expect(new Rational(9n, 2n).round()).toBe(4n); // 4.5 -> 4 (even)
      expect(new Rational(7n, 2n).round()).toBe(4n); // 3.5 -> 4 (even)
      expect(new Rational(-5n, 2n).round()).toBe(-2n); // -2.5 -> -2 (even)
    });

    it('should round away from zero', () => {
      expect(new Rational(5n, 2n).round('away')).toBe(3n);
      expect(new Rational(-5n, 2n).round('away')).toBe(-3n);
    });

    it('should round toward zero', () => {
      expect(new Rational(5n, 2n).round('toward')).toBe(2n);
      expect(new Rational(-5n, 2n).round('toward')).toBe(-2n);
    });

    it('should round up', () => {
      expect(new Rational(5n, 2n).round('up')).toBe(3n);
      expect(new Rational(-5n, 2n).round('up')).toBe(-2n);
    });

    it('should round down', () => {
      expect(new Rational(5n, 2n).round('down')).toBe(2n);
      expect(new Rational(-5n, 2n).round('down')).toBe(-3n);
    });

    it('should round not-half cases correctly', () => {
      expect(new Rational(4n, 3n).round()).toBe(1n); // 1.333... -> 1
      expect(new Rational(5n, 3n).round()).toBe(2n); // 1.666... -> 2
    });

    it('should truncate toward zero', () => {
      expect(new Rational(7n, 3n).trunc()).toBe(2n);
      expect(new Rational(-7n, 3n).trunc()).toBe(-2n);
    });
  });

  describe('predicates', () => {
    it('should detect zero', () => {
      expect(Rational.zero().isZero()).toBe(true);
      expect(new Rational(1n, 2n).isZero()).toBe(false);
    });

    it('should detect one', () => {
      expect(Rational.one().isOne()).toBe(true);
      expect(new Rational(2n, 2n).isOne()).toBe(true);
      expect(new Rational(1n, 2n).isOne()).toBe(false);
    });

    it('should detect integer', () => {
      expect(new Rational(6n, 2n).isInteger()).toBe(true);
      expect(new Rational(1n, 2n).isInteger()).toBe(false);
    });

    it('should detect positive/negative', () => {
      expect(new Rational(1n, 2n).isPositive()).toBe(true);
      expect(new Rational(-1n, 2n).isPositive()).toBe(false);
      expect(Rational.zero().isPositive()).toBe(false);

      expect(new Rational(-1n, 2n).isNegative()).toBe(true);
      expect(new Rational(1n, 2n).isNegative()).toBe(false);
      expect(Rational.zero().isNegative()).toBe(false);
    });
  });

  describe('static methods', () => {
    it('should create zero and one', () => {
      expect(Rational.zero().numerator).toBe(0n);
      expect(Rational.one().numerator).toBe(1n);
    });

    it('should create from bigint', () => {
      const r = Rational.from(7n);
      expect(r.numerator).toBe(7n);
      expect(r.denominator).toBe(1n);
    });

    it('should create from integer number', () => {
      const r = Rational.from(5);
      expect(r.numerator).toBe(5n);
      expect(r.denominator).toBe(1n);
    });

    it('should create from fraction string', () => {
      const r = Rational.fromString('3/4');
      expect(r.numerator).toBe(3n);
      expect(r.denominator).toBe(4n);
    });

    it('should create from negative fraction string', () => {
      const r = Rational.fromString('-3/4');
      expect(r.numerator).toBe(-3n);
      expect(r.denominator).toBe(4n);
    });

    it('should create from decimal string', () => {
      const r = Rational.fromString('1.5');
      expect(r.numerator).toBe(3n);
      expect(r.denominator).toBe(2n);
    });

    it('should create from tuple', () => {
      const r = Rational.fromTuple([5n, 7n]);
      expect(r.numerator).toBe(5n);
      expect(r.denominator).toBe(7n);
    });
  });

  describe('asIntegerRatio', () => {
    it('should return tuple', () => {
      const [n, d] = new Rational(3n, 4n).asIntegerRatio();
      expect(n).toBe(3n);
      expect(d).toBe(4n);
    });
  });

  describe('height', () => {
    it('should return max of |numerator| and denominator', () => {
      expect(new Rational(2n, 3n).height()).toBe(3n);
      expect(new Rational(34n, 3n).height()).toBe(34n);
      expect(new Rational(-97n, 4n).height()).toBe(97n);
    });
  });

  describe('large numbers', () => {
    it('should handle large numerator', () => {
      const big = 2n ** 100n;
      const r = new Rational(big, 1n);
      expect(r.numerator).toBe(big);
    });

    it('should handle large denominator', () => {
      const big = 2n ** 100n;
      const r = new Rational(1n, big);
      expect(r.denominator).toBe(big);
    });

    it('should reduce large numbers', () => {
      const big = 2n ** 100n;
      const r = new Rational(big, big * 2n);
      expect(r.numerator).toBe(1n);
      expect(r.denominator).toBe(2n);
    });

    it('should multiply large numbers', () => {
      const a = new Rational(2n ** 50n, 3n ** 30n);
      const b = new Rational(3n ** 30n, 2n ** 50n);
      const c = a.mul(b);
      expect(c.isOne()).toBe(true);
    });
  });
});

describe('RationalField (QQ)', () => {
  describe('singleton', () => {
    it('should return same instance', () => {
      expect(RationalField.getInstance()).toBe(QQ);
    });
  });

  describe('__call__', () => {
    it('should create from bigint', () => {
      const r = QQ.__call__(7n);
      expect(r.numerator).toBe(7n);
      expect(r.denominator).toBe(1n);
    });

    it('should create from string', () => {
      const r = QQ.__call__('3/4');
      expect(r.numerator).toBe(3n);
      expect(r.denominator).toBe(4n);
    });

    it('should create from two bigints', () => {
      const r = QQ.__call__(3n, 4n);
      expect(r.numerator).toBe(3n);
      expect(r.denominator).toBe(4n);
    });

    it('should create from tuple', () => {
      const r = QQ.__call__([5n, 7n]);
      expect(r.numerator).toBe(5n);
      expect(r.denominator).toBe(7n);
    });

    it('should pass through Rational', () => {
      const r = new Rational(3n, 4n);
      expect(QQ.__call__(r)).toBe(r);
    });

    it('should create from number', () => {
      const r = QQ.__call__(1.5);
      expect(r.numerator).toBe(3n);
      expect(r.denominator).toBe(2n);
    });
  });

  describe('field properties', () => {
    it('should return correct zero', () => {
      expect(QQ.zero().isZero()).toBe(true);
    });

    it('should return correct one', () => {
      expect(QQ.one().isOne()).toBe(true);
    });

    it('should be a field', () => {
      expect(QQ.is_field()).toBe(true);
    });

    it('should have characteristic 0', () => {
      expect(QQ.characteristic()).toBe(0n);
    });

    it('should be a ring', () => {
      expect(QQ.is_ring()).toBe(true);
    });

    it('should be an integral domain', () => {
      expect(QQ.is_integral_domain()).toBe(true);
    });

    it('should be a prime field', () => {
      expect(QQ.is_prime_field()).toBe(true);
    });

    it('should be absolute', () => {
      expect(QQ.is_absolute()).toBe(true);
    });

    it('should have degree 1', () => {
      expect(QQ.degree()).toBe(1n);
      expect(QQ.absolute_degree()).toBe(1n);
    });
  });

  describe('generators', () => {
    it('should have 1 generator', () => {
      expect(QQ.ngens()).toBe(1n);
    });

    it('should return gens tuple', () => {
      const gens = QQ.gens();
      expect(gens.length).toBe(1);
      expect(gens[0].isOne()).toBe(true);
    });

    it('should return gen(0)', () => {
      expect(QQ.gen(0).isOne()).toBe(true);
    });

    it('should throw for gen(n) where n != 0', () => {
      expect(() => QQ.gen(1)).toThrow('n must be 0');
    });
  });

  describe('zeta (roots of unity)', () => {
    it('should return 1 for zeta(1)', () => {
      expect(QQ.zeta(1).isOne()).toBe(true);
    });

    it('should return -1 for zeta(2)', () => {
      const z = QQ.zeta(2);
      expect(z.numerator).toBe(-1n);
      expect(z.denominator).toBe(1n);
    });

    it('should throw for zeta(n) where n > 2', () => {
      expect(() => QQ.zeta(3)).toThrow('no n-th root of unity in rational field');
    });
  });

  describe('other properties', () => {
    it('should have discriminant 1', () => {
      expect(QQ.discriminant()).toBe(1n);
    });

    it('should have class number 1', () => {
      expect(QQ.class_number()).toBe(1n);
    });

    it('should have signature (1, 0)', () => {
      const [r, c] = QQ.signature();
      expect(r).toBe(1n);
      expect(c).toBe(0n);
    });
  });

  describe('string representation', () => {
    it('should have correct toString', () => {
      expect(QQ.toString()).toBe('Rational Field');
    });

    it('should have correct LaTeX', () => {
      expect(QQ._latex_()).toBe('\\Bold{Q}');
    });
  });

  describe('an_element and some_elements', () => {
    it('should return 1/2 for an_element', () => {
      const e = QQ.an_element();
      expect(e.numerator).toBe(1n);
      expect(e.denominator).toBe(2n);
    });

    it('should generate some_elements', () => {
      const elements = [...QQ.some_elements()];
      expect(elements.length).toBeGreaterThan(5);
      // Check that 0, 1, -1 are included
      expect(elements.some((e) => e.isZero())).toBe(true);
      expect(elements.some((e) => e.isOne())).toBe(true);
      expect(elements.some((e) => e.numerator === -1n && e.denominator === 1n)).toBe(true);
    });
  });
});

describe('continued_fraction_list', () => {
  it('should compute standard continued fraction of 13/9', () => {
    const cf = new Rational(13n, 9n).continued_fraction_list();
    expect(cf).toEqual([1n, 2n, 4n]);
    // Verify: 1 + 1/(2 + 1/4) = 1 + 1/(9/4) = 1 + 4/9 = 13/9
  });

  it('should compute standard continued fraction of 225/157', () => {
    const cf = new Rational(225n, 157n).continued_fraction_list();
    expect(cf).toEqual([1n, 2n, 3n, 4n, 5n]);
  });

  it('should compute continued fraction of -1/3', () => {
    const cf = new Rational(-1n, 3n).continued_fraction_list();
    expect(cf).toEqual([-1n, 1n, 2n]);
  });

  it('should compute continued fraction of integers', () => {
    expect(new Rational(1n).continued_fraction_list()).toEqual([1n]);
    expect(new Rational(0n).continued_fraction_list()).toEqual([0n]);
    expect(new Rational(-1n).continued_fraction_list()).toEqual([-1n]);
  });

  it('should compute Hirzebruch-Jung continued fraction', () => {
    const cf = new Rational(11n, 19n).continued_fraction_list('hj');
    // 11/19 = 1 - 1/(3 - 1/(2 - 1/(3 - 1/2)))
    expect(cf).toEqual([1n, 3n, 2n, 3n, 2n]);
  });
});

describe('valuation', () => {
  it('should compute valuation at prime dividing numerator', () => {
    expect(new Rational(-5n, 9n).valuation(5n)).toBe(1n);
  });

  it('should compute valuation at prime dividing denominator', () => {
    expect(new Rational(-5n, 9n).valuation(3n)).toBe(-2n);
  });

  it('should return 0 for prime not in factorization', () => {
    expect(new Rational(-5n, 9n).valuation(2n)).toBe(0n);
  });

  it('should return Infinity for zero', () => {
    expect(new Rational(0n).valuation(4n)).toBe('Infinity');
  });

  it('should handle non-prime bases', () => {
    expect(new Rational(7n, 16n).valuation(4n)).toBe(-2n);
  });
});

describe('is_square and sqrt', () => {
  it('should detect perfect squares', () => {
    expect(new Rational(9n, 4n).is_square()).toBe(true);
    expect(new Rational(49n, 100n).is_square()).toBe(true);
  });

  it('should reject non-squares', () => {
    expect(new Rational(4n, 3n).is_square()).toBe(false);
    expect(new Rational(-1n, 4n).is_square()).toBe(false);
  });

  it('should compute sqrt of perfect square', () => {
    const r = new Rational(25n, 9n).sqrt() as Rational;
    expect(r.numerator).toBe(5n);
    expect(r.denominator).toBe(3n);
  });

  it('should return both roots with all=true', () => {
    const roots = new Rational(25n, 9n).sqrt({ all: true }) as Rational[];
    expect(roots.length).toBe(2);
    expect(roots[0].numerator).toBe(5n);
    expect(roots[1].numerator).toBe(-5n);
  });

  it('should throw for non-rational sqrt when extend=false', () => {
    expect(() => new Rational(2n, 3n).sqrt({ extend: false })).toThrow(
      'square root of 2/3 not a rational number'
    );
  });
});

describe('is_nth_power and nth_root', () => {
  it('should detect perfect squares', () => {
    expect(new Rational(25n, 4n).is_nth_power(2n)).toBe(true);
    expect(new Rational(9n, 2n).is_nth_power(2n)).toBe(false);
  });

  it('should detect perfect cubes', () => {
    expect(new Rational(125n, 8n).is_nth_power(3n)).toBe(true);
    expect(new Rational(-125n, 8n).is_nth_power(3n)).toBe(true);
  });

  it('should detect even roots of negative numbers', () => {
    expect(new Rational(-25n).is_nth_power(2n)).toBe(false);
  });

  it('should compute nth root', () => {
    const r = new Rational(125n, 8n).nth_root(3n);
    expect(r.numerator).toBe(5n);
    expect(r.denominator).toBe(2n);
  });

  it('should compute nth root of negative', () => {
    const r = new Rational(-125n, 8n).nth_root(3n);
    expect(r.numerator).toBe(-5n);
    expect(r.denominator).toBe(2n);
  });

  it('should compute negative power nth root', () => {
    const r = new Rational(25n, 4n).nth_root(-2n);
    expect(r.numerator).toBe(2n);
    expect(r.denominator).toBe(5n);
  });

  it('should throw for non-perfect power', () => {
    expect(() => new Rational(9n, 2n).nth_root(2n)).toThrow('not a perfect 2nd power');
  });
});

describe('support', () => {
  it('should return primes in factorization', () => {
    expect(new Rational(-4n, 17n).support()).toEqual([2n, 17n]);
  });

  it('should return empty for 1/-1', () => {
    expect(new Rational(1n).support()).toEqual([]);
    expect(new Rational(-1n).support()).toEqual([]);
  });

  it('should throw for zero', () => {
    expect(() => new Rational(0n).support()).toThrow('Support of 0 not defined');
  });
});

describe('prime_to_S_part', () => {
  it('should remove powers of primes in S', () => {
    const r = new Rational(3n, 4n).prime_to_S_part([2n]);
    expect(r.numerator).toBe(3n);
    expect(r.denominator).toBe(1n);
  });

  it('should handle multiple primes', () => {
    const r = new Rational(700n, 99n).prime_to_S_part([2n, 3n, 5n]);
    expect(r.numerator).toBe(7n);
    expect(r.denominator).toBe(11n);
  });

  it('should handle zero', () => {
    const r = new Rational(0n).prime_to_S_part([2n, 3n]);
    expect(r.isZero()).toBe(true);
  });
});

describe('is_S_unit and is_S_integral', () => {
  it('should detect S-units', () => {
    expect(new Rational(1n, 2n).is_S_unit()).toBe(false);
    expect(new Rational(1n, 2n).is_S_unit([2n])).toBe(true);
  });

  it('should detect S-integrals', () => {
    expect(new Rational(1n, 2n).is_S_integral([])).toBe(false);
    expect(new Rational(1n, 2n).is_S_integral([2n])).toBe(true);
    expect(new Rational(3n).is_S_integral([])).toBe(true);
  });
});

describe('period', () => {
  it('should compute period of 1/7', () => {
    expect(new Rational(1n, 7n).period()).toBe(6n);
  });

  it('should return 1 for terminating decimals', () => {
    expect(new Rational(1n, 8n).period()).toBe(1n);
    expect(new Rational(1n, 5n).period()).toBe(1n);
    expect(new Rational(3n, 4n).period()).toBe(1n);
  });
});

describe('str', () => {
  it('should convert to base 10', () => {
    expect(new Rational(-4n, 17n).str()).toBe('-4/17');
    expect(new Rational(5n).str()).toBe('5');
  });

  it('should convert to base 2', () => {
    expect(new Rational(-4n, 17n).str(2)).toBe('-100/10001');
  });

  it('should throw for invalid base', () => {
    expect(() => new Rational(1n).str(1)).toThrow('base (=1) must be between 2 and 36');
    expect(() => new Rational(1n).str(40)).toThrow('base (=40) must be between 2 and 36');
  });
});

describe('ndigits and nbits', () => {
  it('should count digits', () => {
    // 123/45 reduces to 41/15, which has 2 + 2 = 4 digits
    expect(new Rational(123n, 45n).ndigits()).toBe(4n);
    expect(new Rational(0n, 1n).ndigits()).toBe(2n); // 1 + 1 digits (0/1)
  });

  it('should count bits', () => {
    const r = new Rational(8n, 4n); // reduces to 2/1
    // 2 in binary is "10" (2 bits), 1 in binary is "1" (1 bit)
    expect(r.nbits()).toBe(3n);
  });
});

describe('local_height and global_height', () => {
  it('should compute local height', () => {
    const r = new Rational(25n, 6n);
    expect(r.local_height(2n)).toBeCloseTo(Math.log(2), 10);
    expect(r.local_height(3n)).toBeCloseTo(Math.log(3), 10);
    expect(r.local_height(5n)).toBe(0); // 5 only in numerator
  });

  it('should compute global height', () => {
    const r = new Rational(6n, 25n);
    expect(r.global_height()).toBeCloseTo(Math.log(25), 10);
  });
});

describe('content', () => {
  it('should compute content of two rationals', () => {
    const a = new Rational(2n, 3n);
    const content = a.content(new Rational(2n, 3n));
    expect(content.numerator).toBe(2n);
    expect(content.denominator).toBe(3n);
  });

  it('should compute content with different rationals', () => {
    const a = new Rational(2n, 3n);
    const content = a.content(new Rational(1n, 5n));
    expect(content.numerator).toBe(1n);
    expect(content.denominator).toBe(15n);
  });
});

describe('factorial', () => {
  it('should compute factorial of integer', () => {
    const r = new Rational(5n).factorial();
    expect(r.numerator).toBe(120n);
    expect(r.denominator).toBe(1n);
  });

  it('should throw for non-integer', () => {
    expect(() => new Rational(1n, 2n).factorial()).toThrow('factorial only defined for integers');
  });

  it('should throw for negative', () => {
    expect(() => new Rational(-1n).factorial()).toThrow(
      'factorial not defined for negative integers'
    );
  });
});

describe('val_unit', () => {
  it('should return valuation and unit', () => {
    const [v, u] = new Rational(-4n, 17n).val_unit(2n);
    expect(v).toBe(2n);
    expect(u.numerator).toBe(-1n);
    expect(u.denominator).toBe(17n);
  });

  it('should handle zero', () => {
    const [v, u] = new Rational(0n).val_unit(2n);
    expect(v).toBe('Infinity');
    expect(u.isOne()).toBe(true);
  });
});

describe('integration tests', () => {
  it('should perform complex arithmetic', () => {
    // (2/3 + 1/4) * (5/6 - 1/2) = (11/12) * (1/3) = 11/36
    const a = QQ.__call__(2n, 3n);
    const b = QQ.__call__(1n, 4n);
    const c = QQ.__call__(5n, 6n);
    const d = QQ.__call__(1n, 2n);

    const result = a.add(b).mul(c.sub(d));
    expect(result.numerator).toBe(11n);
    expect(result.denominator).toBe(36n);
  });

  it('should handle SageMath-style example: 19/374 + 17/371', () => {
    const a = QQ.__call__('19/374');
    const b = QQ.__call__('17/371');
    const sum = a.add(b);
    expect(sum.numerator).toBe(13407n);
    expect(sum.denominator).toBe(138754n);
  });

  it('should compute 374 * (19/374) = 19', () => {
    const a = QQ.__call__('19/374');
    const result = a.mul(374n);
    expect(result.numerator).toBe(19n);
    expect(result.denominator).toBe(1n);
  });

  it('should handle inverse correctly: ~(-4/17) = -17/4', () => {
    const r = QQ.__call__(-4n, 17n);
    const inv = r.inv();
    expect(inv.numerator).toBe(-17n);
    expect(inv.denominator).toBe(4n);
  });
});

describe('Rational - minpoly and charpoly', () => {
  it('minpoly for 1/3 should be x - 1/3', () => {
    const r = new Rational(1n, 3n);
    const poly = r.minpoly();

    // The polynomial is [-self, 1] = [-1/3, 1]
    // Coefficient of x^0 is -1/3
    // Coefficient of x^1 is 1
    expect(poly.degree()).toBe(1);
    expect(poly.getCoeff(0).numerator).toBe(-1n);
    expect(poly.getCoeff(0).denominator).toBe(3n);
    expect(poly.getCoeff(1).numerator).toBe(1n);
    expect(poly.getCoeff(1).denominator).toBe(1n);
  });

  it('minpoly for 2 should be x - 2', () => {
    const r = new Rational(2n, 1n);
    const poly = r.minpoly();

    expect(poly.degree()).toBe(1);
    expect(poly.getCoeff(0).numerator).toBe(-2n);
    expect(poly.getCoeff(0).denominator).toBe(1n);
    expect(poly.getCoeff(1).numerator).toBe(1n);
  });

  it('charpoly equals minpoly for rationals', () => {
    const r = new Rational(5n, 7n);
    const minp = r.minpoly();
    const charp = r.charpoly();

    expect(charp.degree()).toBe(minp.degree());
    for (let i = 0; i <= charp.degree(); i++) {
      expect(charp.getCoeff(i).numerator).toBe(minp.getCoeff(i).numerator);
      expect(charp.getCoeff(i).denominator).toBe(minp.getCoeff(i).denominator);
    }
  });

  it('minpoly can use custom variable name', () => {
    const r = new Rational(3n, 4n);
    const poly = r.minpoly('y');
    // The polynomial is over QQ['y']
    expect(poly.parent.variable_name).toBe('y');
  });

  it('minpoly for 0 should be x', () => {
    const r = new Rational(0n, 1n);
    const poly = r.minpoly();

    expect(poly.degree()).toBe(1);
    expect(poly.getCoeff(0).numerator).toBe(0n);
    expect(poly.getCoeff(1).numerator).toBe(1n);
  });

  it('minpoly for negative rational', () => {
    const r = new Rational(-2n, 3n);
    const poly = r.minpoly();

    // x - (-2/3) = x + 2/3
    expect(poly.degree()).toBe(1);
    expect(poly.getCoeff(0).numerator).toBe(2n);
    expect(poly.getCoeff(0).denominator).toBe(3n);
  });
});

// ============================================
// Tests for new RationalField methods
// ============================================

import { is_RationalField } from './rational_field.js';
import { ZZ } from './integer_ring.js';

describe('RationalField additional methods', () => {
  describe('is_finite', () => {
    it('should return false since QQ is infinite', () => {
      expect(QQ.is_finite()).toBe(false);
    });
  });

  describe('discriminants', () => {
    it('absolute_discriminant should return 1', () => {
      expect(QQ.absolute_discriminant()).toBe(1n);
    });

    it('relative_discriminant should return 1', () => {
      expect(QQ.relative_discriminant()).toBe(1n);
    });

    it('discriminant should return 1', () => {
      expect(QQ.discriminant()).toBe(1n);
    });
  });

  describe('order', () => {
    it('should return Infinity', () => {
      expect(QQ.order()).toBe('Infinity');
    });
  });

  describe('maximal_order and ring_of_integers', () => {
    it('maximal_order should return ZZ', () => {
      expect(QQ.maximal_order()).toBe(ZZ);
    });

    it('ring_of_integers should return ZZ', () => {
      expect(QQ.ring_of_integers()).toBe(ZZ);
    });
  });

  describe('number_field', () => {
    it('should return self (QQ)', () => {
      expect(QQ.number_field()).toBe(QQ);
    });
  });

  describe('power_basis', () => {
    it('should return [1]', () => {
      const pb = QQ.power_basis();
      expect(pb.length).toBe(1);
      expect(pb[0].isOne()).toBe(true);
    });
  });

  describe('iteration by height', () => {
    it('should iterate in height order', () => {
      const rationals: string[] = [];
      let count = 0;
      for (const r of QQ) {
        rationals.push(r.toString());
        count++;
        if (count >= 17) break;
      }
      expect(rationals).toEqual([
        '0',
        '1',
        '-1',
        '1/2',
        '-1/2',
        '2',
        '-2',
        '1/3',
        '-1/3',
        '3',
        '-3',
        '2/3',
        '-2/3',
        '3/2',
        '-3/2',
        '1/4',
        '-1/4',
      ]);
    });

    it('should produce heights in order', () => {
      const heights: bigint[] = [];
      let count = 0;
      for (const r of QQ) {
        heights.push(r.height());
        count++;
        if (count >= 17) break;
      }
      expect(heights).toEqual([1n, 1n, 1n, 2n, 2n, 2n, 2n, 3n, 3n, 3n, 3n, 3n, 3n, 3n, 3n, 4n, 4n]);
    });
  });

  describe('range_by_height', () => {
    it('should return rationals with height < 4', () => {
      const list = [...QQ.range_by_height(4)].map((r) => r.toString());
      expect(list).toEqual([
        '0',
        '1',
        '-1',
        '1/2',
        '-1/2',
        '2',
        '-2',
        '1/3',
        '-1/3',
        '3',
        '-3',
        '2/3',
        '-2/3',
        '3/2',
        '-3/2',
      ]);
    });

    it('should return rationals with height exactly 2', () => {
      const list = [...QQ.range_by_height(2, 3)].map((r) => r.toString());
      expect(list).toEqual(['1/2', '-1/2', '2', '-2']);
    });

    it('should return empty list for invalid range', () => {
      expect([...QQ.range_by_height(3, 3)]).toEqual([]);
      expect([...QQ.range_by_height(10, 1)]).toEqual([]);
    });

    it('should return empty list for height <= 0', () => {
      expect([...QQ.range_by_height(-10, 1)]).toEqual([]);
    });
  });

  describe('primes_of_bounded_norm_iter', () => {
    it('should yield primes up to 10', () => {
      expect([...QQ.primes_of_bounded_norm_iter(10)]).toEqual([2n, 3n, 5n, 7n]);
    });

    it('should yield primes up to 20', () => {
      expect([...QQ.primes_of_bounded_norm_iter(20)]).toEqual([2n, 3n, 5n, 7n, 11n, 13n, 17n, 19n]);
    });

    it('should return empty for B < 2', () => {
      expect([...QQ.primes_of_bounded_norm_iter(1)]).toEqual([]);
      expect([...QQ.primes_of_bounded_norm_iter(0)]).toEqual([]);
    });

    it('should include boundary prime', () => {
      expect([...QQ.primes_of_bounded_norm_iter(7)]).toEqual([2n, 3n, 5n, 7n]);
    });
  });

  describe('selmer_generators', () => {
    it('should return [-1] for empty S and m=2', () => {
      expect(QQ.selmer_generators([], 2)).toEqual([-1n]);
    });

    it('should return [-1, 3] for S=[3] and m=2', () => {
      expect(QQ.selmer_generators([3n], 2)).toEqual([-1n, 3n]);
    });

    it('should return [-1, 5] for S=[5] and m=2', () => {
      expect(QQ.selmer_generators([5n], 2)).toEqual([-1n, 5n]);
    });

    it('should return with orders when requested', () => {
      const [gens, ords] = QQ.selmer_generators([2n, 3n, 5n, 7n], 2, { orders: true }) as [
        bigint[],
        bigint[],
      ];
      expect(gens).toEqual([-1n, 2n, 3n, 5n, 7n]);
      expect(ords).toEqual([2n, 2n, 2n, 2n, 2n]);
    });

    it('should not include -1 for odd m', () => {
      expect(QQ.selmer_generators([2n, 3n, 5n, 7n], 3)).toEqual([2n, 3n, 5n, 7n]);
    });

    it('should have correct orders for m=3', () => {
      const [gens, ords] = QQ.selmer_generators([2n, 3n, 5n, 7n], 3, { orders: true }) as [
        bigint[],
        bigint[],
      ];
      expect(gens).toEqual([2n, 3n, 5n, 7n]);
      expect(ords).toEqual([3n, 3n, 3n, 3n]);
    });
  });

  describe('selmer_group_iterator', () => {
    it('should iterate QQ(empty, 2)', () => {
      const list = [...QQ.selmer_group_iterator([], 2)].map((r) => r.toString());
      expect(list).toEqual(['1', '-1']);
    });

    it('should iterate QQ({2}, 2)', () => {
      const list = [...QQ.selmer_group_iterator([2n], 2)].map((r) => r.toString());
      expect(list).toEqual(['1', '2', '-1', '-2']);
    });

    it('should iterate QQ({2,3}, 2)', () => {
      const list = [...QQ.selmer_group_iterator([2n, 3n], 2)].map((r) => r.toString());
      expect(list).toEqual(['1', '3', '2', '6', '-1', '-3', '-2', '-6']);
    });

    it('should iterate QQ({5}, 2)', () => {
      const list = [...QQ.selmer_group_iterator([5n], 2)].map((r) => r.toString());
      expect(list).toEqual(['1', '5', '-1', '-5']);
    });
  });

  describe('quadratic_defect', () => {
    it('should return Infinity for 0', () => {
      expect(QQ.quadratic_defect(0, 7)).toBe('Infinity');
    });

    it('should return 0 for quadratic residue mod p', () => {
      // 5 is not a quadratic residue mod 7, so defect is 0
      expect(QQ.quadratic_defect(5, 7)).toBe(0n);
    });

    it('should return Infinity for quadratic residue', () => {
      // 2 is a quadratic residue mod 7 (since 3^2 = 9 = 2 mod 7)
      expect(QQ.quadratic_defect(2, 7)).toBe('Infinity');
    });

    it('should handle p=2 with u=1 mod 8', () => {
      // 1 mod 8: Infinity
      expect(QQ.quadratic_defect(1, 2)).toBe('Infinity');
    });

    it('should handle p=2 with u=5 mod 8', () => {
      // 5 mod 8: v + 2
      expect(QQ.quadratic_defect(5, 2)).toBe(2n);
    });

    it('should handle p=2 with u=3 mod 8', () => {
      // 3 mod 8: v + 1
      expect(QQ.quadratic_defect(3, 2)).toBe(1n);
    });

    it('should handle p=2 with u=7 mod 8', () => {
      // 7 mod 8: v + 1
      expect(QQ.quadratic_defect(7, 2)).toBe(1n);
    });

    it('should return valuation for odd valuation', () => {
      // valuation(5, 5) = 1 (odd), so return 1
      expect(QQ.quadratic_defect(5, 5)).toBe(1n);
    });

    it('should throw for non-prime p when check=true', () => {
      expect(() => QQ.quadratic_defect(5, 4)).toThrow('4 must be prime');
    });

    it('should accept Rational input', () => {
      const r = new Rational(5n, 1n);
      expect(QQ.quadratic_defect(r, 7)).toBe(0n);
    });
  });
});

describe('is_RationalField', () => {
  it('should return true for QQ', () => {
    expect(is_RationalField(QQ)).toBe(true);
  });

  it('should return false for ZZ', () => {
    expect(is_RationalField(ZZ)).toBe(false);
  });

  it('should return false for numbers', () => {
    expect(is_RationalField(42)).toBe(false);
  });

  it('should return false for null', () => {
    expect(is_RationalField(null)).toBe(false);
  });

  it('should return false for Rational instances', () => {
    expect(is_RationalField(new Rational(1n, 2n))).toBe(false);
  });
});
