import { expect, spyOn, test } from 'bun:test';
import { ZM_mul, ZV_prod } from './ZV.js';
import * as modular from './FpV.js';
import { F2m_mul } from './F2v.js';

test('PARI matrix column adapters, empty dimensions and modulus one', () => {
  const a = [[], [0n, 1n, 3n], [0n, 2n, 4n]];
  const b = [[], [0n, 5n, 7n], [0n, 6n, 8n]];
  const copy = structuredClone([a, b]);
  expect(ZM_mul(a, b)).toEqual([[], [0n, 19n, 43n], [0n, 22n, 50n]]);
  expect(modular.FpM_mul(a, b, 1n)).toEqual([[], [0n, 0n, 0n], [0n, 0n, 0n]]);
  expect(
    modular.Flm_mul(
      a.map((c) => c.map(() => 0n)),
      b.map((c) => c.map(() => 0n)),
      1n
    )
  ).toEqual([[], [0n, 0n, 0n], [0n, 0n, 0n]]);
  expect(ZM_mul([[]], [[], [0n], [0n]])).toEqual([[], [0n], [0n]]);
  expect(F2m_mul([3n, 2n], [3n, 1n], 2)).toEqual([1n, 3n]);
  expect([a, b]).toEqual(copy);
  expect(ZV_prod([2n, 3n, 5n])).toBe(30n);
});

test('PARI large balanced integer matrices use modular reconstruction', () => {
  const a = [[], ...Array.from({ length: 70 }, () => [0n, ...new Array<bigint>(70).fill(1n)])];
  const spy = spyOn(modular, 'Flm_mul');
  try {
    const result = ZM_mul(a, a);
    expect(result.slice(1).every((c) => c.slice(1).every((v) => v === 70n))).toBe(true);
    expect(spy).toHaveBeenCalled();
    for (const [, , p] of spy.mock.calls) expect(p > 1n << 63n).toBe(true);
  } finally {
    spy.mockRestore();
  }
});
