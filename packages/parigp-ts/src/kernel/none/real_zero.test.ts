import { expect, test } from 'bun:test';
import { divir, itor, mulir, real_0_bit, rtor, shiftr } from '../../qfb.js';
import { invr } from './mp_indep.js';
import { rdivii } from './level1.js';

test('native zero conversions retain allocation and clamp absolute accuracy', () => {
  expect(itor(0n, 128)).toEqual({ s: 0, e: -128, m: 0n, p: 128 });
  expect(rtor(real_0_bit(-179), 64)).toEqual({ s: 0, e: -179, m: 0n, p: 64 });
  expect(rtor(shiftr(itor(0n, 128), 77), 128)).toEqual({ s: 0, e: -128, m: 0n, p: 128 });
});

test('allocated zero accuracy affects integer multiplication and division', () => {
  const x = shiftr(itor(0n, 128), 77);
  expect(mulir(0n, x)).toEqual({ s: 0, e: -179, m: 0n, p: 0 });
  expect(divir(0n, x)).toEqual({ s: 0, e: -77, m: 0n, p: 0 });
});

test('allocated zero inverse failures originate in the division kernel', () => {
  expect(() => divir(1n, itor(0n, 128))).toThrow('impossible inverse in divrr');
  expect(() => invr(itor(0n, 128))).toThrow('impossible inverse in divrr');
  expect(() => invr(real_0_bit(-128))).toThrow('impossible inverse in invr');
});

test('native integer quotient assigns zero before checking the denominator', () => {
  expect(rdivii(0n, 0n, 128)).toEqual({ s: 0, e: -128, m: 0n, p: 128 });
});
