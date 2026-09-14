import { expect, test } from 'bun:test';
import {
  mzd_init,
  mzd_add,
  mzd_init_window,
  mzd_mul_naive,
  mzd_mul_m4rm,
  mzd_mul,
  mzd_make_table,
} from './index.js';

test('native packed products, alias squares, Gray tables and windows', () => {
  const a = mzd_init(2, 2, [3n, 2n]);
  for (const multiply of [mzd_mul_naive, mzd_mul_m4rm, mzd_mul])
    expect(multiply(a, a).rows).toEqual([1n, 2n]);
  expect(a.rows).toEqual([3n, 2n]);
  expect(mzd_add(a, a).rows).toEqual([0n, 0n]);
  expect(mzd_make_table(a, 0, 2)).toEqual({ rows: [0n, 3n, 1n, 2n], lookup: [0, 1, 3, 2] });
  expect(mzd_init_window(mzd_init(1, 65, [(1n << 64n) | 1n]), 0, 64, 1, 65).rows).toEqual([1n]);
});

test('pinned native zero-window cutoff faults become recoverable errors', () => {
  const a = mzd_init(86, 86);
  expect(() => mzd_mul(a, a, 64)).toThrow('Aborted');
  const left = mzd_init(128, 86),
    right = mzd_init(86, 128);
  expect(() => mzd_mul(left, right, 64)).toThrow('Segmentation fault');
  expect(mzd_mul(a, a, 128).rows).toEqual(Array(86).fill(0n));
});
