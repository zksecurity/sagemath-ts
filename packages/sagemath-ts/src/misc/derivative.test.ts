import { expect, test } from 'bun:test';
import { Integer } from '../rings/integer_ring.js';
import { derivative_parse } from './derivative.js';

test('derivative parser keeps expanded singleton lists and separates consecutive counts', () => {
  const variables = [null, 'x', 'y'];
  expect(derivative_parse([variables])).toBe(variables);
  expect(derivative_parse(['x', 2n, 2n, 'y'])).toEqual(['x', 'x', null, null, 'y']);
  expect(derivative_parse([0n])).toEqual([]);
  expect(derivative_parse([])).toEqual([null]);
});

test('derivative parser follows Cython count conversion before checking negativity', () => {
  expect(() => derivative_parse([new Integer(-1n)])).toThrow(
    'derivative counts must be nonnegative'
  );
  expect(() => derivative_parse([-(1n << 31n) - 1n])).toThrow('value too large to convert to int');
  expect(() => derivative_parse([1n << 63n])).toThrow('Python int too large to convert to C long');
});
