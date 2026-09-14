import { expect, test } from 'bun:test';
import { _nmod_poly_resultant } from './resultant.js';

test('native modular resultant retains the same-array constant shortcut', () => {
  const a = [1n];
  expect(_nmod_poly_resultant(a, a, 7n)).toBe(0n);
  expect(_nmod_poly_resultant(a, [1n], 7n)).toBe(1n);
});
test('half-GCD resultant handles a constant remainder and odd-degree sign', () => {
  const b = [2n, ...Array<bigint>(626).fill(0n), 1n];
  const a = [...b];
  a[0] = 3n;
  expect(_nmod_poly_resultant(a, b, 7n)).toBe(6n);
});
