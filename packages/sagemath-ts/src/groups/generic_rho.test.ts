import { describe, expect, test } from 'bun:test';
import { current_randstate, set_random_seed } from '../misc/randstate.js';
import { Mod } from '../rings/finite_rings/integer_mod.js';
import { discrete_log_rho } from './generic.js';

// Each state/error value is compared to bundled Sage by gg_rho_state, with
// the documented DJB2 string hash supplied explicitly on the Sage side.
describe('Pollard rho native draw and validation regressions', () => {
  test('draws m then n through the modular ring CPython stream', () => {
    set_random_seed(0);
    expect(discrete_log_rho(Mod(17n, 401n), Mod(1n, 401n), 401n, '+')).toBe(17n);
    expect(current_randstate().python_random().getrandbits(64)).toBe(731121675548719989n);
  });

  test('an identity target follows the native walk and consumes its draws', () => {
    set_random_seed(0);
    expect(discrete_log_rho(Mod(0n, 401n), Mod(1n, 401n), 401n, '+')).toBe(0n);
    expect(current_randstate().python_random().getrandbits(64)).toBe(731121675548719989n);
  });

  test('nonprime orders use the source prime-order error', () => {
    for (const order of [-7n, 0n, 1n, 4n]) {
      expect(() => discrete_log_rho(Mod(1n, 17n), Mod(2n, 17n), order, '+')).toThrow(
        'for Pollard rho algorithm the order of the group must be prime'
      );
    }
  });

  test('retry exhaustion retains its error and random-state effects', () => {
    set_random_seed(0);
    expect(() => discrete_log_rho(Mod(2n, 401n), Mod(0n, 401n), 401n, '+')).toThrow(
      'Pollard rho algorithm failed to find a logarithm'
    );
    expect(current_randstate().python_random().getrandbits(64)).toBe(1691914492695921706n);
  });

  test('the documented multiplicative example uses a prime-order subgroup', () => {
    set_random_seed(0);
    const base = Mod(4n, 1019n);
    expect(discrete_log_rho(base.pow(250n), base, 509n, '*')).toBe(250n);
  });
});
