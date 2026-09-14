import { expect, test } from 'bun:test';
import { power } from './ZZ_pEX.js';
import { power as primePower } from './ZZ_pX.js';

test('NTL powers retain coefficient context and reduction', () => {
  expect(power([[0n, 1n]], 2n, [1n, 1n, 1n], 7n)).toEqual([[6n, 6n]]);
  expect(power([], 0n, [1n, 1n, 1n], 7n)).toEqual([[1n]]);
  expect(primePower([1n, 1n], 7n, 7n)).toEqual([1n, 0n, 0n, 0n, 0n, 0n, 0n, 1n]);
});
