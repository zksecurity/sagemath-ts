import { expect, spyOn, test } from 'bun:test';
import * as word from './Flx.js';
import * as field from './FpX.js';
import { Flx_ddf, Flx_nbfact_by_degree, FpX_nbfact } from './FpX_factor.js';
import * as powers from './bb_group.js';
import { FpX_mul } from './ffinit.js';
import { FpX_ddf, _galconj_factor_squarefree_irreducibles } from './galconj.js';
const f = [5n, 2n, 6n, 0n, 0n, 2n, 6n, 3n, 1n];
test('PARI word distinct-degree factors retain their native scale while factor callers normalize', () => {
  const expected: Array<[number, bigint[]]> = [
    [1, [2n, 4n, 5n]],
    [2, [4n, 4n, 3n]],
    [4, [5n, 1n, 3n, 6n, 1n]],
  ];
  expect([...FpX_ddf(f, 7n)]).toEqual(expected);
  expect(Flx_ddf(f, 7n).map(([g, d]) => [d, g])).toEqual(expected);
  expect(FpX_nbfact(f, 7n)).toBe(4);
  expect(Flx_nbfact_by_degree(f, 7n)).toEqual({ D: [0, 2, 1, 0, 1, 0, 0, 0, 0], nb: 4 });
  const factors = _galconj_factor_squarefree_irreducibles(f, 7n);
  expect(factors.every((g) => g.at(-1) === 1n)).toBe(true);
  expect(factors.reduce((a, b) => FpX_mul(a, b, 7n), [1n])).toEqual(f);
});
test('word DDF selects Frobenius powering or composition using the native cost threshold', () => {
  const schedule = spyOn(powers, 'gen_pow_i');
  try {
    const T = [1n, ...new Array<bigint>(49).fill(0n), 1n];
    Flx_ddf(T, 3n);
    expect(schedule).toHaveBeenCalledTimes(5);
    schedule.mockClear();
    Flx_ddf(T, (1n << 64n) - 59n);
    expect(schedule).toHaveBeenCalledTimes(1);
  } finally {
    schedule.mockRestore();
  }
});
test('FpX factor counts preserve and convert the native preinverse at word boundaries', () => {
  const wi = spyOn(word, 'Flx_invBarrett');
  const fi = spyOn(field, 'FpX_invBarrett');
  try {
    // FpX_get_red already supplies an inverse. Word conversion must retain it.
    const T = [1n, ...new Array<bigint>(90).fill(0n), 1n];
    FpX_nbfact(T, 7n);
    expect(fi).toHaveBeenCalledTimes(1);
    expect(wi).not.toHaveBeenCalled();
  } finally {
    wi.mockRestore();
    fi.mockRestore();
  }
});

test('dense DDF products split large packed intermediates before the slow engine regime', async () => {
  const packing = await import('./_polynomial_packing.js');
  const split = spyOn(packing, 'splitProduct');
  try {
    const p = (1n << 64n) - 59n;
    const n = 242;
    const a = new Array<bigint>(n).fill(p - 1n);
    const product = word.Flx_mul(a, a, p);
    expect(split).toHaveBeenCalled();
    expect(product).toEqual(
      Array.from({ length: 2 * n - 1 }, (_, i) => BigInt(i < n ? i + 1 : 2 * n - 1 - i))
    );
  } finally {
    split.mockRestore();
  }
});
