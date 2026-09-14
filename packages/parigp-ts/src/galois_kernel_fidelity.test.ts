import { expect, test } from 'bun:test';
import { FpM_ker } from './galconj.js';

const matrix = (n: number) =>
  Array.from({ length: n }, (_, i) =>
    Array.from({ length: n }, (_, j) => BigInt(i === 0 && j < 2))
  );
const basis = (n: number, c: bigint) =>
  Array.from({ length: n - 1 }, (_, i) =>
    Array.from({ length: n }, (_, j) => (j === i + 1 ? 1n : i === 0 && j === 0 ? c : 0n))
  );

test('generic recursive kernels preserve native signed coefficients', () => {
  const p = (1n << 127n) - 1n;
  expect(FpM_ker(matrix(5), 5, 5, p)).toEqual(basis(5, -1n));
  expect(FpM_ker(matrix(8), 8, 8, p)).toEqual(basis(8, -1n));
});

test('kernel thresholds preserve small and word-field representatives', () => {
  const p = (1n << 127n) - 1n;
  expect(FpM_ker(matrix(4), 4, 4, p)).toEqual(basis(4, p - 1n));
  expect(FpM_ker(matrix(5), 5, 5, 17n)).toEqual(basis(5, 16n));
  expect(FpM_ker(matrix(8), 8, 8, 17n)).toEqual(basis(8, 16n));
});

test('packed binary and ternary kernels retain inputs and basis order', () => {
  const input = matrix(65),
    saved = input.map((r) => r.slice());
  expect(FpM_ker(input, 65, 65, 2n)).toEqual(basis(65, 1n));
  expect(FpM_ker(input, 65, 65, 3n)).toEqual(basis(65, 2n));
  expect(input).toEqual(saved);
});
