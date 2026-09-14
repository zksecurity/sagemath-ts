import { expect, test } from 'bun:test';
import { GF } from '../finite_rings/index.js';
import { QQ } from '../rational_field.js';
import { FunctionField } from './index.js';

test('variable maps convert scalars, dictionaries and nested singleton lists', () => {
  const K = FunctionField(GF(3n) as never, 'x');
  const [L, fromL, toL] = K.change_variable_name('y');
  expect(toL(1n).parent).toBe(L);
  expect(String(toL([[1n, 2n]]))).toBe('2*y + 1');
  expect(String(fromL(new Map([[2n, 1n], [0n, 1n]])))).toBe('x^2 + 1');
  expect(() => toL(new Map([['bad', 1n]]))).toThrow("{'bad': 1} fails to convert into the map's domain");
  expect(() => toL(new Map([["a'b", 1n]]))).toThrow("{\"a'b\": 1} fails to convert into the map's domain");
});

test('identity variable maps coerce into their domain and retain existing elements', () => {
  const K = FunctionField(QQ as never, 'x');
  const [, from, to] = K.change_variable_name('x');
  expect(from).toBe(to);
  expect(from(K.gen())).toBe(K.gen());
  expect(to(1n).parent).toBe(K);
  expect(String(to(null))).toBe('0');
});
