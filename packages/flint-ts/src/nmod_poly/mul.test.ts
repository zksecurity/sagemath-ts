import { expect, test } from 'bun:test';
import {
  _nmod_poly_mul,
  _nmod_poly_mul_classical,
  _nmod_poly_mul_KS,
  _nmod_poly_mul_KS2,
  _nmod_poly_mul_KS4,
} from './mul.js';
import { _nmod_poly_add } from './add.js';
import { _nmod_poly_sub } from './sub.js';
import { _nmod_poly_pow } from './pow.js';

test('native product kernels preserve aliasing, inputs, and distributivity at word boundaries', () => {
  for (const p of [2n, 14n, 65537n, (1n << 64n) - 59n]) {
    const a = Array.from({ length: 97 }, (_, i) => p - BigInt(i % 7) - 1n),
      b = Array.from({ length: 103 }, (_, i) => BigInt(i * i) - p),
      c = [3n, -7n, 1n];
    const original = a.slice();
    const expected = _nmod_poly_mul_classical(a, b, p),
      square = _nmod_poly_mul_classical(a, a, p);
    for (const mul of [_nmod_poly_mul, _nmod_poly_mul_KS, _nmod_poly_mul_KS2, _nmod_poly_mul_KS4]) {
      expect(mul(a, b, p)).toEqual(expected);
      expect(mul(a, a, p)).toEqual(square);
      expect(mul(a, _nmod_poly_add(b, c, p), p)).toEqual(
        _nmod_poly_add(mul(a, b, p), mul(a, c, p), p)
      );
    }
    expect(a).toEqual(original);
    expect(_nmod_poly_sub(_nmod_poly_add(a, b, p), b, p)).toEqual(_nmod_poly_mul(a, [1n], p));
    expect(_nmod_poly_pow(a, 3n, p)).toEqual(_nmod_poly_mul(square, a, p));
  }
});
test('native power handles empty polynomials and unsigned-word exponents', () => {
  expect(_nmod_poly_pow([], 0n, 7n)).toEqual([1n]);
  expect(_nmod_poly_pow([], 3n, 7n)).toEqual([]);
  expect(_nmod_poly_pow([6n], (1n << 64n) - 1n, 7n)).toEqual([6n]);
  expect(() => _nmod_poly_pow([1n], -1n, 7n)).toThrow(RangeError);
  expect(() => _nmod_poly_pow([1n], 1n << 64n, 7n)).toThrow(RangeError);
  for (const p of [0n, 1n << 64n])
    for (const fn of [_nmod_poly_mul, _nmod_poly_add, _nmod_poly_sub])
      expect(() => fn([1n], [1n], p)).toThrow(RangeError);
});

test('native zero-ring kernels include the zero-to-zero power', () => {
  for (const f of [
    _nmod_poly_mul,
    _nmod_poly_mul_classical,
    _nmod_poly_mul_KS,
    _nmod_poly_mul_KS2,
    _nmod_poly_mul_KS4,
    _nmod_poly_add,
    _nmod_poly_sub,
  ])
    expect(f([1n, 2n], [3n, 4n], 1n)).toEqual([]);
  expect(_nmod_poly_pow([], 0n, 1n)).toEqual([]);
  expect(_nmod_poly_pow([1n, 1n], (1n << 64n) - 1n, 1n)).toEqual([]);
});
