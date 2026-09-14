import { expect, spyOn, test } from 'bun:test';
import * as field from './FpX.js';
import * as word from './Flx.js';
import * as composition from './_polynomial_composition.js';
import { FpXQ_minpoly } from './galconj.js';
import { setrand, getrand } from './random.js';

test('Shoup minimal polynomials delegate word sampling and share one Barrett inverse', () => {
  const saved = getrand();
  const wr = spyOn(word, 'random_Flx'),
    fr = spyOn(field, 'random_FpX');
  const wi = spyOn(word, 'Flx_invBarrett'),
    fi = spyOn(field, 'FpX_invBarrett');
  try {
    for (const [p, degree, isWord] of [
      [7n, 91, true],
      [(1n << 64n) - 59n, 30, true],
      [(1n << 64n) + 13n, 36, false],
    ] as const) {
      for (const n of [degree - 1, degree]) {
        wr.mockClear();
        fr.mockClear();
        wi.mockClear();
        fi.mockClear();
        setrand(1n);
        const T = [1n, ...new Array<bigint>(n - 1).fill(0n), 1n];
        expect(FpXQ_minpoly([2n], T, p)).toEqual([p - 2n, 1n]);
        expect(isWord ? wr : fr).toHaveBeenCalled();
        expect(isWord ? fr : wr).not.toHaveBeenCalled();
        expect(isWord ? wi : fi).toHaveBeenCalledTimes(n === degree ? 1 : 0);
      }
    }
  } finally {
    wr.mockRestore();
    fr.mockRestore();
    wi.mockRestore();
    fi.mockRestore();
    setrand(saved);
  }
});

test('large-prime minimal polynomials reduce before the power table and restore reproducible streams', () => {
  const saved = getrand(),
    powers = spyOn(composition, 'quotientPowers');
  try {
    const p = (1n << 64n) + 13n;
    setrand(1n);
    expect(FpXQ_minpoly([-1n, 0n, 1n], [1n, 0n, 1n], p)).toEqual([2n, 1n]);
    expect(powers.mock.calls[0]![0]).toEqual([p - 2n]);
    const state = getrand();
    setrand(1n);
    expect(FpXQ_minpoly([p - 2n], [1n, 0n, 1n], p)).toEqual([2n, 1n]);
    expect(getrand()).toBe(state);
    expect(() => FpXQ_minpoly([-1n, 0n, 1n], [1n, 0n, 1n], 7n)).toThrow(
      'power polynomial exceeds modulus degree'
    );
  } finally {
    powers.mockRestore();
    setrand(saved);
  }
});
