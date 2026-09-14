import { expect, test, spyOn } from 'bun:test';
import * as pari from '@sagemath-ts/parigp-ts';
import { GF2 } from '../finite_rings/gf2.js';
import { PolynomialRing } from './polynomial_ring.js';

test('binary polynomial factors delegate to PARI and sort by Sage multiplicity', () => {
  const call = spyOn(pari, 'F2x_factor');
  try {
    const R = new PolynomialRing(GF2, 'x'),
      f = R.__call__([0n, 0n, 1n, 1n]);
    expect(f.factor().map(([g, e]) => [String(g), e])).toEqual([
      ['x + 1', 1],
      ['x', 2],
    ]);
    expect(call).toHaveBeenCalledWith(12n);
  } finally {
    call.mockRestore();
  }
});
