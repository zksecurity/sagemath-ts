import { expect, test } from 'bun:test';
import { Integer } from './integer_ring.js';
import {
  two_squares_pyx,
  three_squares_pyx,
  four_squares_pyx,
  is_sum_of_two_squares_pyx,
} from './sum_of_squares.js';

test('native sum-of-squares wrappers preserve source tuples and failures', () => {
  expect(two_squares_pyx(new Integer(106n))).toEqual([5n, 9n]);
  expect(three_squares_pyx(107n)).toEqual([1n, 5n, 9n]);
  expect(four_squares_pyx(15447n)).toEqual([2n, 5n, 17n, 123n]);
  expect(is_sum_of_two_squares_pyx(21n)).toBe(false);
  expect(is_sum_of_two_squares_pyx(106n)).toBe(true);
  expect(() => two_squares_pyx(21n)).toThrow('21 is not a sum of 2 squares');
  expect(() => three_squares_pyx(7n)).toThrow('7 is not a sum of 3 squares');
  for (const f of [
    two_squares_pyx,
    three_squares_pyx,
    four_squares_pyx,
    is_sum_of_two_squares_pyx,
  ]) {
    expect(() => f(-1n)).toThrow("can't convert negative value to uint32_t");
    expect(() => f(1n << 32n)).toThrow('value too large to convert to uint32_t');
    expect(() => f(1n << 64n)).toThrow('Python int too large to convert to C unsigned long');
  }
});
