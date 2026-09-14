import { expect, test } from 'bun:test';
import { sorted } from './python_sort.js';

test('Python sorting is stable and leaves its input alone', () => {
  const input = [
    { key: 2, id: 0 },
    { key: 1, id: 1 },
    { key: 1, id: 2 },
  ];
  expect(sorted(input, (a, b) => a.key < b.key).map((x) => x.id)).toEqual([1, 2, 0]);
  expect(input.map((x) => x.id)).toEqual([0, 1, 2]);
  const empty: number[] = [];
  expect(sorted(empty, (a, b) => a < b)).not.toBe(empty);
});

test('comparison errors propagate without changing the caller array', () => {
  const input = [3, 1, 2];
  expect(() =>
    sorted(input, () => {
      throw new Error('comparison failed');
    })
  ).toThrow('comparison failed');
  expect(input).toEqual([3, 1, 2]);
});
