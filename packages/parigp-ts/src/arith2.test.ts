import { expect, test } from 'bun:test';
import { eulerphi, numdiv } from './arith2.js';

test('native arith2 signed values and zero conventions', () => {
  expect(eulerphi(-12n)).toBe(4n);
  expect(eulerphi(0n)).toBe(2n);
  expect(numdiv(-12n)).toBe(6n);
  expect(() => numdiv(0n)).toThrow('domain error in numdiv: argument = 0');
});

test('large prime-power totient and divisor count', () => {
  const n = (1n << 256n) * 3n ** 127n;
  expect(eulerphi(n)).toBe((1n << 256n) * 3n ** 126n);
  expect(numdiv(n)).toBe(257n * 128n);
});
