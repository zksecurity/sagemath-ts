import { expect, test } from 'bun:test';
import { fixedfieldsurmer } from './galconj.js';

test('seventeen weights search beyond the first coefficient vector', () => {
  const w = [0, ...Array.from({ length: 17 }, (_, i) => i + 1)];
  const sums = [[], ...Array.from({ length: 17 }, () => [0n, 0n, 0n])];
  sums[1] = [0n, 0n, 1n];
  expect(fixedfieldsurmer(257n, sums, w)).toEqual({
    v: [0, 1, ...Array(15).fill(0), 1],
    w,
  });
});

test('search carries through seven base-four coefficient positions', () => {
  const w = [0, ...Array.from({ length: 23 }, (_, i) => i + 1)];
  const sums = [[], ...Array.from({ length: 23 }, () => [0n, 0n, 0n])];
  sums[7] = [0n, 0n, 1n];
  expect(fixedfieldsurmer(257n, sums, w)).toEqual({
    v: [0, ...Array(6).fill(0), 1, ...Array(15).fill(0), 1],
    w,
  });
});

test('a fully exhausted search returns native null', () => {
  expect(fixedfieldsurmer(3n, [[], [0n, 1n, 1n], [0n, 2n, 2n]], [0, 1, 2])).toBeNull();
});
