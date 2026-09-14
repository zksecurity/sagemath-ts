import { expect, test } from 'bun:test';
import { F2x_mul, F2x_sqr, F2x_sqrt, F2x_divrem, F2x_valrem, F2xq_powers } from './F2x.js';
import { F2m_ker, F2m_ker_sp } from './F2v.js';
import { F2x_factor, F2x_factor_squarefree } from './FpX_factor.js';
import { setrand, getrand } from './random.js';

test('PARI packed polynomial coefficients and native zero boundaries', () => {
  expect(F2x_mul(7n, 3n)).toBe(9n);
  expect(F2x_sqrt(F2x_sqr((1n << 128n) + 3n))).toBe((1n << 128n) + 3n);
  expect(F2x_divrem((1n << 8192n) + 7n, 1n)).toEqual([(1n << 8192n) + 7n, 0n]);
  expect(() => F2x_divrem(1n, 0n)).toThrow('impossible inverse in F2x_divrem: Vecsmall([0]).');
  expect(F2x_valrem(0n)).toEqual([(1n << 63n) - 1n, 0n]);
  expect(F2xq_powers(8n, 1, 3n)).toEqual([1n, 8n]);
});

test('PARI kernel basis follows native column order and mutation', () => {
  const columns = [3n, 3n];
  expect(F2m_ker_sp(columns, 2)).toEqual([3n]);
  expect(columns).toEqual([3n, 1n]);
  const fresh = [3n, 3n];
  expect(F2m_ker(fresh, 2)).toEqual([3n]);
  expect(fresh).toEqual([3n, 3n]);
  expect(F2m_ker_sp([1n, 2n], 2, 1)).toBeNull();
  expect(F2m_ker_sp([3n, 3n], 2, 1)).toBe(3n);
});

test('PARI binary factor unit, multiplicities and deterministic seed replay', () => {
  setrand(1n);
  const before = getrand();
  expect(F2x_factor(63n)).toEqual([
    [3n, 1],
    [7n, 2],
  ]);
  const after = getrand();
  setrand(before);
  expect(F2x_factor(63n)).toEqual([
    [3n, 1],
    [7n, 2],
  ]);
  expect(getrand()).toBe(after);
  expect(F2x_factor(0n)).toEqual([[0n, 1]]);
  expect(F2x_factor(1n)).toEqual([]);
  expect(F2x_factor_squarefree(63n)).toEqual([3n, 7n]);
});
