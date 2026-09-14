import { expect, test } from 'bun:test';
import { addrr } from './add.js';
import { dbltor } from './mp_indep.js';
import { negr, real_0_bit } from '../../qfb.js';
import { divrr } from '../gmp/mp.js';

// Native exact-rational/exponent comparisons are in pari_real_arithmetic.
test('native one-word addition truncates alignment and keeps storage precision after cancellation', () => {
  expect(addrr(dbltor(1), dbltor(2 ** -64))).toEqual(dbltor(1));
  expect(addrr(dbltor(1), negr(dbltor(1 - 2 ** -52)))).toEqual(dbltor(2 ** -52));
  expect(addrr(dbltor(1), dbltor(-1))).toEqual(real_0_bit(-63));
  expect(addrr(dbltor(1), real_0_bit(0))).toEqual(real_0_bit(0));
});

test('native one-word division preserves the zero accuracy and inverse error', () => {
  expect(divrr(dbltor(3), dbltor(2))).toEqual(dbltor(1.5));
  expect(divrr(real_0_bit(-53), dbltor(2))).toEqual(real_0_bit(-54));
  expect(() => divrr(dbltor(1), real_0_bit(-53))).toThrow('impossible inverse in divrr');
});
