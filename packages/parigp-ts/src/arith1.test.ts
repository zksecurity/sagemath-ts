import { expect, test } from 'bun:test';
import { hilbert } from './arith1.js';

test('native Hilbert valuation parity, real place and modulus guard', () => {
  expect(hilbert(2n, 3n, 2n)).toBe(-1);
  expect(hilbert(-1n, -1n, 0n)).toBe(-1);
  expect(hilbert(0n, -1n, 0n)).toBe(0);
  expect(hilbert(2n ** 64n, 3n, 2n)).toBe(1);
  expect(hilbert(2n ** 65n, 3n, 2n)).toBe(-1);
  expect(() => hilbert(0n, 1n, -1n)).toThrow('not a prime number in hilbertii: -1');
});
