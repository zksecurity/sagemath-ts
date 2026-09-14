import { expect, test } from 'bun:test';
import {
  pari_init_rand,
  pari_rand,
  setrand,
  getrand,
  random_bits,
  random_Fl,
  randomi,
  random_F2x,
  random_zv,
} from './random.js';
import { random_Flx } from './Flx.js';
import { random_FpX } from './FpX.js';

test('native PARI seed-one stream and saved-state replay', () => {
  pari_init_rand();
  const initial = getrand();
  expect([pari_rand(), pari_rand(), pari_rand()]).toEqual([
    13282407956253574712n,
    7557322358563246340n,
    14991082624209354397n,
  ]);
  setrand(initial);
  expect(pari_rand()).toBe(13282407956253574712n);
});

test('native PARI signed word boundary and empty samplers', () => {
  setrand(1n);
  expect(random_bits(64)).toBe(-5164336117455976904n);
  setrand(1n);
  expect(random_zv(3)).toEqual([
    -5164336117455976904n,
    7557322358563246340n,
    -3455661449500197219n,
  ]);
  const before = getrand();
  expect(random_Fl(1n)).toBe(0n);
  expect(randomi(1n)).toBe(0n);
  expect(random_F2x(0)).toBe(0n);
  expect(random_zv(0)).toEqual([]);
  expect(random_Flx(0, 0n)).toEqual([]);
  expect(random_FpX(0, 0n)).toEqual([]);
  expect(getrand()).toBe(before);
});

test('native seed errors and undefined-shift guards preserve the stream', () => {
  setrand(42n);
  const before = getrand();
  expect(() => setrand(0n)).toThrow('domain error in setrand: n <= 0');
  expect(() => setrand(1n << 64n)).toThrow('domain error in setrand: n != getrand()');
  expect(() => random_bits(0)).toThrow('bits must be between 1 and 64');
  expect(() => random_Fl(0n)).toThrow('limit must be a positive unsigned word');
  expect(() => randomi(0n)).toThrow('limit must be positive');
  expect(getrand()).toBe(before);
});
