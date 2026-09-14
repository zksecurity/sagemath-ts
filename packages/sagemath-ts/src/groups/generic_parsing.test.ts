import { describe, expect, test } from 'bun:test';
import { Mod } from '../rings/finite_rings/integer_mod.js';
import {
  bsgs,
  discrete_log_lambda,
  discrete_log_rho,
  multiple,
  order_from_bounds,
} from './generic.js';

describe('Bundled group parsing and operation schedules', () => {
  test('missing custom operations precede invalid numeric arguments', () => {
    const a = Mod(2n, 11n),
      b = Mod(4n, 11n);
    const message =
      'identity, inverse and operation must all be specified when operation is neither addition nor multiplication';
    expect(() => bsgs(a, b, [-1n, 8n], 'other')).toThrow(message);
    expect(() => discrete_log_lambda(b, a, [-1n, 8n], 'other')).toThrow(message);
    expect(() => discrete_log_rho(b, a, 0n, 'other')).toThrow(message);
  });

  test('dividing a negative bound preserves its exact ceiling', () => {
    const d = 2n ** 54n;
    expect(() => order_from_bounds(Mod(1n, 11n), [-d, 0n], d, '+')).toThrow(
      'bsgs() requires 0<=lb<=ub'
    );
  });

  test('an idempotent custom operation is evaluated only once', () => {
    let calls = 0;
    const result = multiple(
      Mod(2n, 11n),
      2n ** 1000n,
      'other',
      Mod(10n, 11n),
      (x) => x.neg(),
      (x, y) => {
        calls++;
        return x.value < y.value ? x : y;
      }
    );
    expect(result.value).toBe(2n);
    expect(calls).toBe(1);
  });

  test('custom multiples follow the original binary schedule', () => {
    const trace: bigint[][] = [];
    const result = multiple(
      Mod(2n, 11n),
      5n,
      'other',
      Mod(0n, 11n),
      (x) => x.neg(),
      (x, y) => {
        trace.push([x.value, y.value]);
        return x.add(y);
      }
    );
    expect(result.value).toBe(10n);
    expect(trace).toEqual([
      [2n, 2n],
      [4n, 4n],
      [8n, 2n],
    ]);
  });

  test('order bounds do not repeat the completed order-multiple check', () => {
    const trace: bigint[][] = [];
    const result = order_from_bounds(
      Mod(2n, 11n),
      [1n, 22n],
      1n,
      'other',
      Mod(0n, 11n),
      (x) => x.neg(),
      (x, y) => {
        trace.push([x.value, y.value]);
        return x.add(y);
      }
    );
    expect(result).toBe(11n);
    expect(trace).toEqual([
      [0n, 2n],
      [2n, 2n],
      [2n, 4n],
      [2n, 6n],
      [2n, 8n],
      [2n, 10n],
      [2n, 1n],
      [2n, 3n],
      [2n, 5n],
      [2n, 7n],
      [2n, 9n],
    ]);
  });
});
