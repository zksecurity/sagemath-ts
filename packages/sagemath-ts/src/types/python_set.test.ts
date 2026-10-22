import { expect, test } from 'bun:test';
import { python_int_hash, python_set } from './python_set.js';

const order = (xs: bigint[]): bigint[] => python_set(xs, python_int_hash, (a, b) => a === b);

// Expected orders are CPython 3.11.8's list(set(xs)).
test('python_set follows CPython slot order, probing and resizing', () => {
  expect(order([1n, 10000n])).toEqual([10000n, 1n]);
  expect(order([5n, 13n, 21n, 29n])).toEqual([29n, 13n, 21n, 5n]);
  expect(order([8n, 0n, 16n, 1n, 9n])).toEqual([0n, 1n, 8n, 9n, 16n]);
  expect(order([-1n, -2n, (1n << 61n) - 1n, (1n << 64n) + 8n, 3n])).toEqual([
    (1n << 61n) - 1n,
    3n,
    (1n << 64n) + 8n,
    -2n,
    -1n,
  ]);
  expect(order([7n, 7n, 7n])).toEqual([7n]);
});

test('python_int_hash reduces modulo 2^61 - 1 and avoids -1', () => {
  expect(python_int_hash((1n << 61n) - 1n)).toBe(0n);
  expect(python_int_hash(1n << 61n)).toBe(1n);
  expect(python_int_hash(-1n)).toBe(-2n);
  expect(python_int_hash(-5n)).toBe(-5n);
});
