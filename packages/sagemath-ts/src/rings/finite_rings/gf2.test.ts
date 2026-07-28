/**
 * Tests for sage/rings/finite_rings/gf2
 *
 * `GF2Element` had no test file at all, which is how the `Number(value)`
 * reduction (M22) survived: every input above 2^53 was rounded to a float
 * before the mod-2 reduction, so `GF(2)(2^64 + 1)` came out 0.
 */

import { describe, expect, test } from 'bun:test';
import { ZeroDivisionError } from '../../errors.js';
import { GF2, GF2Element } from './gf2.js';

describe('GF2 coercion', () => {
  test('small integers reduce mod 2', () => {
    expect(GF2.__call__(0).value).toBe(0);
    expect(GF2.__call__(1).value).toBe(1);
    expect(GF2.__call__(2).value).toBe(0);
    expect(GF2.__call__(7).value).toBe(1);
    expect(GF2.__call__(-1).value).toBe(1);
    expect(GF2.__call__(-2).value).toBe(0);
  });

  test('large bigints reduce exactly (M22)', () => {
    // sage: GF(2)(2^64 + 1)  ->  1
    expect(GF2.__call__(2n ** 64n + 1n).value).toBe(1);
    expect(GF2.__call__(2n ** 64n).value).toBe(0);
    // 2^53 + 1 is the first odd integer that is not exactly representable as
    // a double; Number() rounds it down to 2^53, an even number.
    expect(GF2.__call__(2n ** 53n + 1n).value).toBe(1);
    expect(GF2.__call__(2n ** 200n + 1n).value).toBe(1);
    expect(GF2.__call__(-(2n ** 64n + 1n)).value).toBe(1);
    const p = 2n ** 255n - 19n;
    expect(GF2.__call__(p).value).toBe(1);
  });

  test('booleans and elements round-trip', () => {
    expect(GF2.__call__(true).value).toBe(1);
    expect(GF2.__call__(false).value).toBe(0);
    expect(GF2.__call__(GF2.one()).value).toBe(1);
  });

  test('non-integers are rejected rather than silently kept', () => {
    expect(() => GF2.__call__(1.5)).toThrow('unable to convert 1.5 to an integer');
  });
});

describe('GF2 arithmetic', () => {
  test('addition table', () => {
    const [z, o] = [GF2.zero(), GF2.one()];
    expect(z.add(z).value).toBe(0);
    expect(z.add(o).value).toBe(1);
    expect(o.add(o).value).toBe(0);
    // subtraction is addition
    expect(o.sub(o).value).toBe(0);
    expect(o.sub(z).value).toBe(1);
  });

  test('multiplication table', () => {
    const [z, o] = [GF2.zero(), GF2.one()];
    expect(z.mul(z).value).toBe(0);
    expect(z.mul(o).value).toBe(0);
    expect(o.mul(o).value).toBe(1);
  });

  test('arithmetic with large bigint operands', () => {
    const o = GF2.one();
    expect(o.add(2n ** 64n).value).toBe(1);
    expect(o.add(2n ** 64n + 1n).value).toBe(0);
    expect(o.mul(2n ** 64n + 1n).value).toBe(1);
    expect(o.eq(2n ** 64n + 1n)).toBe(true);
    expect(o.eq(2n ** 64n)).toBe(false);
  });

  test('inverse and division', () => {
    expect(GF2.one().inv().value).toBe(1);
    expect(() => GF2.zero().inv()).toThrow(ZeroDivisionError);
    expect(GF2.one().div(GF2.one()).value).toBe(1);
    expect(() => GF2.one().div(GF2.zero())).toThrow(ZeroDivisionError);
    expect(() => GF2.one().div(2n ** 64n)).toThrow(ZeroDivisionError);
  });

  test('powers', () => {
    // sage: GF(2)(0)^0 -> 1 ; GF(2)(0)^5 -> 0 ; GF(2)(0)^-1 -> ZeroDivisionError
    expect(GF2.zero().pow(0).value).toBe(1);
    expect(GF2.zero().pow(5).value).toBe(0);
    expect(() => GF2.zero().pow(-1)).toThrow(ZeroDivisionError);
    expect(GF2.one().pow(0).value).toBe(1);
    expect(GF2.one().pow(10n ** 20n).value).toBe(1);
    expect(GF2.one().pow(-3).value).toBe(1);
  });

  test('negation is the identity', () => {
    expect(GF2.zero().neg().value).toBe(0);
    expect(GF2.one().neg().value).toBe(1);
  });
});

describe('GF2Field', () => {
  test('field data', () => {
    expect(GF2.characteristic).toBe(2n);
    expect(GF2.order).toBe(2n);
    expect(GF2.degree).toBe(1);
    expect(GF2.cardinality()).toBe(2n);
    expect(GF2.is_field()).toBe(true);
    expect(GF2.toString()).toBe('Finite Field of size 2');
    // sage: GF(2).gen() -> 1
    expect(GF2.gen().value).toBe(1);
  });

  test('iteration yields 0 then 1', () => {
    expect([...GF2].map((x) => x.value)).toEqual([0, 1]);
  });

  test('elements carry the parent', () => {
    const a = new GF2Element(1, GF2);
    expect(a.parent).toBe(GF2);
    expect(a.toBigInt()).toBe(1n);
    expect(a.toString()).toBe('1');
  });
});
