import { describe, expect, test } from 'bun:test';
import { indexpartial, ZX_neg, ZX_Z_mul, ZX_is_squarefree } from './galconj.js';
import { absZ_factor_limit_strict_default } from './ifactor.js';
import { ZX_gcd } from './QX_factor.js';
import { ZpX_reduced_resultant, ZpX_reduced_resultant_fast } from './base2.js';
import { ZpM_echelon, zlm_echelon } from './hnf_snf.js';

describe('integer-polynomial regressions compared with bundled PARI', () => {
  test('polynomial observations normalize trailing zero coefficients', () => {
    expect(ZX_neg([1n, 0n, 0n])).toEqual([-1n]);
    expect(ZX_neg([0n])).toEqual([]);
    expect(ZX_Z_mul([2n, 0n], -3n)).toEqual([-6n]);
  });
  test('nonzero constants are squarefree, while repeated deflated roots are not', () => {
    for (const c of [-2n, -1n, 1n, 2n]) expect(ZX_is_squarefree([c, 0n])).toBe(true);
    expect(ZX_is_squarefree([])).toBe(false);
    expect(ZX_is_squarefree([0n, 0n, 1n])).toBe(false);
    expect(ZX_is_squarefree([1n, 0n, -2n, 0n, 1n])).toBe(false);
    expect(ZX_is_squarefree([1n, 0n, 0n, 0n, 1n])).toBe(true);
  });
  test('indexpartial uses zero factorization and refines prime powers', () => {
    expect(indexpartial([])).toBe(0n);
    expect(indexpartial([1n, -2n, 1n])).toBe(0n);
    expect(indexpartial([-2n, 0n, 0n, 0n, 1n])).toBe(8n);
    expect(indexpartial([0n, 1n], 16n)).toBe(1n);
    expect(() => indexpartial([0n, 0n], 16n)).toThrow('impossible inverse in Fl_inv: Mod(0, 4).');
  });
  test('strict partial factorization retains a large composite power', () => {
    const u = 500009n * 500029n;
    expect(absZ_factor_limit_strict_default(8n * u ** 3n)).toEqual([[[2n, 3n]], [u, 3n]]);
    expect(indexpartial([-2n, 0n, 1n], 8n * u ** 3n)).toBe(2n * u ** 2n);
    expect(absZ_factor_limit_strict_default(499979n * 500009n)).toEqual([
      [
        [499979n, 1n],
        [500009n, 1n],
      ],
      null,
    ]);
  });
  test('modular GCD reconstructs large primitive coefficients', () => {
    const c = (1n << 180n) + 3n;
    expect(ZX_gcd([2n * c, 3n * c + 2n, 3n], [5n * c, 7n * c + 5n, 7n])).toEqual([c, 1n]);
    expect(ZX_gcd([0n, 0n, 2n], [0n, 6n])).toEqual([0n, 2n]);
  });
  test('p-adic resultants agree across word and generic precisions', () => {
    const f = [-2n, 0n, 0n, 0n, 1n],
      df = [0n, 0n, 0n, 4n];
    expect(ZpX_reduced_resultant(f, df, 2n, 4n)).toBe(0n);
    expect(ZpX_reduced_resultant_fast(f, df, 2n, 2)).toBe(4n);
    expect(ZpX_reduced_resultant_fast(f, df, 2n, 65)).toBe(8n);
    expect(ZpX_reduced_resultant(f, df, 2n, 1n << 65n)).toBe(8n);
  });
  test('column echelon keeps native pivot order and aborts at a zero pivot', () => {
    const x = [
      [2n, 0n],
      [0n, 4n],
    ];
    expect(ZpM_echelon(x, false, 2n, 16n)).toEqual(x);
    expect(zlm_echelon(x, false, 2n, 16n)).toEqual(x);
    expect(
      ZpM_echelon(
        [
          [0n, 0n],
          [0n, 0n],
        ],
        true,
        2n,
        16n
      )
    ).toBeNull();
    expect(
      zlm_echelon(
        [
          [0n, 0n],
          [0n, 0n],
        ],
        true,
        2n,
        16n
      )
    ).toBeNull();
  });
});
