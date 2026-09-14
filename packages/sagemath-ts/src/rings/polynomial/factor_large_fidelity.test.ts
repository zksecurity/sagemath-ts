import { expect, spyOn, test } from 'bun:test';
import * as pari from '@sagemath-ts/parigp-ts';
import { GF } from '../finite_rings/index.js';
import { PolynomialRing } from './polynomial_ring.js';

test('large-prime factors delegate to PARI and retain the unit and Sage order', () => {
  const call = spyOn(pari, 'FpX_factor');
  try {
    const p = 18446744073709551629n;
    const R = new PolynomialRing<any>(GF(p), 'x');
    const f = R.__call__([0n, 0n, 2n, 2n]);
    expect(f.factor().map(([g, e]) => [String(g), e])).toEqual([
      ['2', 1],
      ['x + 1', 1],
      ['x', 2],
    ]);
    expect(call).toHaveBeenCalledWith([0n, 0n, 2n, 2n], p);
  } finally {
    call.mockRestore();
  }
});
