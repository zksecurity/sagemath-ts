import { expect, test, spyOn } from 'bun:test';
import * as flint from '@sagemath-ts/flint-ts';
import { GF } from '../finite_rings/index.js';
import { PolynomialRing } from './polynomial_ring.js';

test('word-prime factors delegate to native FLINT and retain Sage ordering', () => {
  const factor = spyOn(flint, 'nmod_poly_factor');
  try {
    const R = new PolynomialRing<any>(GF(5n), 'x'),
      f = R.__call__([4n, 0n, 0n, 0n, 1n]);
    expect(f.factor().map(([q, e]) => [String(q), e])).toEqual([
      ['x + 1', 1],
      ['x + 2', 1],
      ['x + 3', 1],
      ['x + 4', 1],
    ]);
    expect(factor).toHaveBeenCalledWith([4n, 0n, 0n, 0n, 1n], 5n);
  } finally {
    factor.mockRestore();
  }
});
