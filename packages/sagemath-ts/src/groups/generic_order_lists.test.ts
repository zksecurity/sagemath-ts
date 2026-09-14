import { describe, expect, test } from 'bun:test';
import { ArithmeticError } from '../errors.js';
import { Mod } from '../rings/finite_rings/integer_mod.js';
import { order_from_multiple } from './generic.js';

// Values, errors and callback traces are compared to bundled Sage in
// gg_order_list_trace, including former nonterminating valuation loops.
describe('Order-from-multiple factor-list regressions', () => {
  test('empty factor and prime lists request automatic factorization', () => {
    const a = Mod(1n, 12n);
    expect(order_from_multiple(a, 12n, [], '+')).toBe(12n);
    expect(
      order_from_multiple(a, 12n, undefined, '+', undefined, undefined, undefined, { plist: [] })
    ).toBe(12n);
    expect(
      order_from_multiple(a, 12n, [], '+', undefined, undefined, undefined, { plist: [2n, 3n] })
    ).toBe(12n);
  });

  test('invalid prime entries preserve Integer.valuation errors', () => {
    for (const p of [-2n, -1n, 0n, 1n]) {
      expect(() =>
        order_from_multiple(Mod(1n, 12n), 12n, undefined, '+', undefined, undefined, undefined, {
          plist: [p],
        })
      ).toThrow('You can only compute the valuation with respect to a integer larger than 1.');
    }
  });

  test('zero-multiple single-base valuations follow the native helper', () => {
    expect(
      order_from_multiple(Mod(3n, 8n), 0n, undefined, '+', undefined, undefined, undefined, {
        plist: [2n],
      })
    ).toBe(8n);
    expect(
      order_from_multiple(Mod(1n, 9n), 0n, undefined, '+', undefined, undefined, undefined, {
        plist: [3n],
      })
    ).toBe(9n);
  });

  test('infinite multi-factor costs preserve native diagnostics', () => {
    const call = (plist: bigint[]) =>
      order_from_multiple(Mod(1n, 12n), 0n, undefined, '+', undefined, undefined, undefined, {
        plist,
      });
    expect(() => call([2n, 3n])).toThrow(
      "unsupported operand parent(s) for ^: 'The Infinity Ring' and 'The Infinity Ring'"
    );
    expect(() => call([-2n, 2n])).toThrow('Python infinity cannot have complex phase.');
    try {
      call([0n, 2n]);
      throw new Error('expected SignError');
    } catch (error) {
      expect(error).toBeInstanceOf(ArithmeticError);
      expect((error as Error).name).toBe('SignError');
      expect((error as Error).message).toBe('cannot add infinity to minus infinity');
    }
    expect(() => call([1n, 2n])).toThrow('cannot multiply infinity by zero');
  });

  test('large finite valuations retain the exact order', () => {
    expect(
      order_from_multiple(
        Mod(3n, 16n),
        2n ** 1000n,
        undefined,
        '+',
        undefined,
        undefined,
        undefined,
        { plist: [2n] }
      )
    ).toBe(16n);
  });
});
