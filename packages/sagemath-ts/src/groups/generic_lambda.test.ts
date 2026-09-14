import { describe, expect, test } from 'bun:test';
import { set_random_seed } from '../misc/randstate.js';
import { Mod } from '../rings/finite_rings/integer_mod.js';
import { discrete_log_lambda } from './generic.js';

describe('Pollard lambda native walk regressions', () => {
  test('uses the same CPython step-size stream as Sage prandom', () => {
    // Live comparison: gg_lambda_trace(0,0,1009,1,0,4,2,0).
    set_random_seed(0);
    const trace: bigint[] = [];
    const result = discrete_log_lambda(
      Mod(2n, 1009n),
      Mod(1n, 1009n),
      [0n, 4n],
      '+',
      undefined,
      undefined,
      undefined,
      (v) => {
        trace.push(v.value);
        return v.value ** 2n + 1n;
      }
    );
    expect(result).toBe(2n);
    expect(trace).toEqual([4n, 5n, 6n, 2n, 3n, 4n, 5n, 6n]);
  });

  test('normalizes signed hash indices and preserves the selected logarithm', () => {
    // Live comparison: gg_lambda_trace(0,0,17,1,0,256,300,1).
    set_random_seed(0);
    expect(
      discrete_log_lambda(
        Mod(300n, 17n),
        Mod(1n, 17n),
        [0n, 256n],
        '+',
        undefined,
        undefined,
        undefined,
        (v) => -(v.value ** 2n + 1n)
      )
    ).toBe(232n);
  });

  test('preserves the native invalid-bounds error text', () => {
    expect(() => discrete_log_lambda(Mod(2n, 17n), Mod(1n, 17n), [-1n, 16n], '+')).toThrow(
      'discrete_log_lambda() requires 0<=lb<=ub'
    );
    expect(() => discrete_log_lambda(Mod(2n, 17n), Mod(1n, 17n), [16n, 8n], '+')).toThrow(
      'discrete_log_lambda() requires 0<=lb<=ub'
    );
  });
});
