import { expect, test } from 'bun:test';
import {
  _nmod_poly_sqrt, _nmod_poly_sqrt_series, _nmod_poly_mulhigh, _nmod_poly_mulhigh_classical,
  _gr_poly_sqrt_series_basecase, _gr_poly_sqrt_series_newton,
  _gr_poly_rsqrt_series_basecase, _gr_poly_inv_series_basecase,
} from '../index.js';

test('native polynomial roots retain the FLINT representative and valuation', () => {
  expect(_nmod_poly_sqrt([4n,8n,4n],17n)).toEqual([2n,2n]);
  expect(_nmod_poly_sqrt([0n,0n,4n,8n,4n],17n)).toEqual([0n,2n,2n]);
  expect(_nmod_poly_sqrt([1n,0n,1n],17n)).toBeNull();
  expect(_nmod_poly_sqrt([1n,0n,1n],2n)).toEqual([1n,1n]);
  expect(_nmod_poly_sqrt([],17n)).toEqual([]);
});
test('native series algorithms agree at a forced Newton crossover', () => {
  const f=[1n,2n,1n],inverse=[1n,16n,1n,16n,1n,16n,1n,16n];
  expect(_nmod_poly_sqrt_series(f,8,17n)).toEqual([1n,1n]);
  expect(_gr_poly_sqrt_series_basecase(f,8,17n)).toEqual([1n,1n]);
  expect(_gr_poly_sqrt_series_newton(f,8,2,17n)).toEqual([1n,1n]);
  expect(_gr_poly_rsqrt_series_basecase(f,8,17n)).toEqual(inverse);
  expect(_gr_poly_inv_series_basecase([1n,1n],8,17n)).toEqual(inverse);
});
test('native classical high products zero the unused low coefficients', () => {
  expect(_nmod_poly_mulhigh([1n,2n,3n],[4n,5n],3,17n)).toEqual([0n,0n,0n,15n]);
  expect(_nmod_poly_mulhigh_classical([1n,2n,3n],[4n,5n],3,17n)).toEqual([0n,0n,0n,15n]);
});
