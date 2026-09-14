import { expect, test } from 'bun:test';
import { FpX_factor_squarefree } from './index.js';
import { FpX_factor_squarefree as legacyPath } from './galconj.js';

test('native squarefree components are indexed by multiplicity, not irreducibility', () => {
  expect(FpX_factor_squarefree([2n, 3n, 1n], 101n)).toEqual([[2n, 3n, 1n]]);
  expect(FpX_factor_squarefree([1n, 2n, 1n], 101n)).toEqual([[1n], [1n, 1n]]);
  expect(FpX_factor_squarefree([1n, 0n, 0n, 0n, 1n], 2n)).toEqual([
    [1n], [1n], [1n], [1n, 1n],
  ]);
  expect(legacyPath).toBe(FpX_factor_squarefree);
});

test('native dispatch preserves nonmonic squarefree input and distinct constant behavior', () => {
  const large = (1n << 127n) - 1n;
  for (const p of [101n, large]) {
    expect(FpX_factor_squarefree([2n, 2n], p)).toEqual([[2n, 2n]]);
  }
  expect(FpX_factor_squarefree([2n], 101n)).toEqual([]);
  expect(FpX_factor_squarefree([2n], large)).toEqual([[2n]]);
  expect(() => FpX_factor_squarefree([], 101n)).toThrow(
    'impossible inverse in Flx_divrem: Vecsmall([0])'
  );
  // Native large-prime zero input violates its low-level precondition and segfaults.
  expect(() => FpX_factor_squarefree([], large)).toThrow(RangeError);
});
