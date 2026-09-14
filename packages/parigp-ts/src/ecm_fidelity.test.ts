import { expect, test } from 'bun:test';
import { ECM, ECM_loop } from './_ecm.js';
import { ellfacteur } from './ifactor.js';

test('ECM batched inversions preserve which native factor is returned', () => {
  expect(ECM_loop(4295229443n, 8, 1, 142)).toBe(65537n);
});
test('ECM keeps the native whole-modulus failure behavior', () => {
  expect(ECM_loop(4295229443n, 8, 65537, 142)).toBeNull();
});
test('ECM continuation follows the native helix and product schedule', () => {
  expect(ECM_loop(100000054360729711807701683n, 16, 1, 2625)).toBe(1000000000217297n);
  expect(ECM_loop(10000000220546470070824045789n, 8, 1, 2625)).toBeNull();
});
test('ECM seeds advance exactly beyond the safe-number boundary', () => {
  expect(ECM_loop(4295229443n, 8, 2 ** 53, 142)).toBe(65539n);
  expect(ECM_loop(4295229443n, 4, 2 ** 60, 142)).toBe(65539n);
});
test('ECM preserves its native state across repeated rounds', () => {
  const state = new ECM(10000271611173029019820873n, 8, 1);
  for (let i = 0; i < 3; i++) expect(state.round(200)).toBeNull();
});
test('the insisting ECM driver retains its explicit round limit', () => {
  expect(ellfacteur(1000003n * 1000033n, true, 0)).toBeNull();
});
