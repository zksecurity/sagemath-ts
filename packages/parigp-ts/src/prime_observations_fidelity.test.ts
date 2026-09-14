import { expect, spyOn, test } from 'bun:test';
import { FpX_deriv, FpX_eval, FpX_center, FpX_div_by_X_x } from './galconj.js';
import * as group from './bb_group.js';
import { PariError } from './errors.js';

test('derivatives reduce the integer derivative after zero normalization', () => {
  expect(FpX_deriv([1n, 0n, 0n], 0n)).toEqual([]);
  expect(() => FpX_deriv([0n, 1n], 0n)).toThrow(new PariError('impossible inverse in dvmdii: 0.'));
});

test('evaluation shortcuts and linear cancellation precede modulus use', () => {
  expect(FpX_eval([], 123n, 0n)).toBe(0n);
  expect(FpX_eval([0n, 0n], 123n, 0n)).toBe(0n);
  expect(FpX_eval([-17n, 1n], 17n, 0n)).toBe(0n);
  expect(FpX_eval([0n, 1n], 0n, 0n)).toBe(0n);
  expect(FpX_eval([1n, 2n, 3n], -1n, -17n)).toBe(2n);
});

test('sparse evaluation powers zero gaps with the native logarithmic schedule', () => {
  const powers = spyOn(group, 'gen_powu_i');
  try {
    expect(FpX_eval([1n, ...Array<bigint>(127).fill(0n), 1n], 2n, 17n)).toBe(2n);
    expect(powers).toHaveBeenCalledTimes(1);
    expect(powers.mock.calls[0]![1]).toBe(128n);
    expect(FpX_eval([...Array<bigint>(128).fill(0n), 1n], 2n, 17n)).toBe(1n);
    expect(powers).toHaveBeenCalledTimes(2);
    expect(FpX_eval([0n, 0n, 1n], -9n, (1n << 64n) + 13n)).toBe(81n);
    expect(powers).toHaveBeenCalledTimes(2);
  } finally {
    powers.mockRestore();
  }
});

test('coefficient centering preserves native magnitude comparison and result storage', () => {
  expect(FpX_center([17n], 17n, 8n)).toEqual([0n]);
  expect(FpX_center([35n], 17n, 8n)).toEqual([18n]);
  expect(FpX_center([-10n], 17n, -8n)).toEqual([-27n]);
  const f = [0n, 0n];
  expect(FpX_center(f, 0n, -8n)).toEqual([]);
  expect(f).toEqual([0n, 0n]);
});

test('linear division does not evaluate an unrequested remainder', () => {
  const f = [1n, 1n];
  expect(FpX_div_by_X_x(f, 17n, 0n)).toEqual([1n]);
  expect(FpX_div_by_X_x([0n, 0n], 17n, 0n)).toEqual([]);
  expect(f).toEqual([1n, 1n]);
});
