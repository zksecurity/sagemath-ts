import { expect, test } from 'bun:test';
import {
  flint_rand_init,
  flint_rand_set_seed,
  flint_rand_get_seed,
  flint_rand_clear,
} from '../flint.js';
import { n_randlimb, n_randint, n_randbits, n_randtest_bits } from './randomisation.js';

test('FLINT default word stream and copied seed state match the bundled C routines', () => {
  const state = flint_rand_init();
  const saved = flint_rand_get_seed(state);
  expect(n_randlimb(state)).toBe(15737102703946861599n);
  expect(flint_rand_get_seed(state)).toEqual([14886521159149723162n, 15737102702505497715n]);
  expect(saved).toEqual([13845646450878251009n, 13142370077570254774n]);
  flint_rand_clear(state);
  expect(flint_rand_get_seed(state)).toEqual([14886521159149723162n, 15737102702505497715n]);
});

test('FLINT zero-bit samples leave state unchanged and word seeds wrap', () => {
  const state = flint_rand_init();
  flint_rand_set_seed(state, 1n << 80n, -1n);
  expect(flint_rand_get_seed(state)).toEqual([0n, (1n << 64n) - 1n]);
  expect(n_randbits(state, 0)).toBe(0n);
  expect(flint_rand_get_seed(state)).toEqual([0n, (1n << 64n) - 1n]);
  expect(() => n_randtest_bits(state, 65)).toThrow('bits must be between 0 and 64');
  expect(n_randint(flint_rand_init(), 17n)).toBe(14n);
});
