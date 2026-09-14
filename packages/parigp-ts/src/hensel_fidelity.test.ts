import { expect, test } from 'bun:test';
import {
  ZpX_liftroot,
  ZpX_liftroots,
  ZpX_roots,
  ZpX_liftfact,
  bezout_lift_fact,
  ZpX_ZpXQ_liftroot,
} from './galconj.js';
import { PariError } from './errors.js';

test('native single-factor normalization preserves an already monic polynomial', () => {
  expect(ZpX_liftfact([4n, 1n], [[], [0n, 1n]], 2n, 1)).toEqual([[], [4n, 1n]]);
  expect(ZpX_liftfact([-4n, -1n], [[], [0n, 1n]], 2n, 2)).toEqual([[], [0n, 1n]]);
});
test('native nonmonic factor and Bezout lifting', () => {
  const f = [-9n, -2n, -1n],
    Q = [[], [0n, 1n], [2n, 1n]];
  expect(ZpX_liftfact(f, Q, 3n, 3)).toEqual([[], [18n, 1n], [11n, 1n]]);
  expect(bezout_lift_fact(f, Q, 3n, 3)).toEqual([[], [10n, 23n], [18n, 4n]]);
});
test('root lift uses the original-prime inverse error', () => {
  const error = new PariError('impossible inverse in Fp_inv: Mod(3, 3).');
  expect(() => ZpX_liftroot([0n, 0n, 1n], 0n, 3n, 5)).toThrow(error);
  expect(() => ZpX_liftroots([0n, 0n, 1n], [0n, 0n], 3n, 5)).toThrow(error);
});
test('factor-tree validation follows native ordering', () => {
  expect(() => ZpX_liftfact([1n], [[]], 3n, 0)).toThrow(
    new PariError('domain error in MultiLift: #(modular factors) < 2')
  );
  expect(() => ZpX_liftfact([0n, 0n, 1n], [[], [0n, 1n], [0n, 1n]], 3n, 2)).toThrow(
    new PariError('elements not coprime in BuildTree:\n    x\n    x')
  );
  expect(ZpX_liftfact([0n, 0n, 1n], [[], [0n, 1n], [0n, 1n]], 3n, 1)).toEqual([
    [],
    [0n, 1n],
    [0n, 1n],
  ]);
});
test('quotient-root precision one retains the original representative', () => {
  expect(ZpX_ZpXQ_liftroot([-1n, 0n, 1n], [4n], [1n, 0n, 1n], 3n, 1)).toEqual([4n]);
  expect(ZpX_ZpXQ_liftroot([-1n, 0n, 1n], [4n], [1n, 0n, 1n], 3n, 3)).toEqual([1n]);
  expect(ZpX_ZpXQ_liftroot([-2n, 0n, 1n], [3n, 1n], [-2n, 0n, 1n], 3n, 5)).toEqual([0n, 1n]);
});
test('invalid quotient roots stop before native unchecked exact division', () => {
  expect(() => ZpX_ZpXQ_liftroot([-1n, 0n, 1n], [0n, 1n], [1n, 0n, 1n], 3n, 2)).toThrow(
    new RangeError('Hensel lifting requires exact polynomial division')
  );
});
test('root-list lifting dispatches through the full factor tree', () => {
  expect(ZpX_liftroots([4n, 1n], [0n, 0n], 2n, 1)).toEqual([0n, 0n]);
  expect(ZpX_roots([4n, 1n], 2n, 3)).toEqual([0n, 4n]);
});

test('multiple-factor Bezout precision one stops before the native uninitialized tree', () => {
  expect(() => bezout_lift_fact([0n, -1n, 1n], [[], [0n, 1n], [2n, 1n]], 3n, 1)).toThrow(
    new RangeError('bezout_lift_fact requires precision at least 2 for multiple factors')
  );
});
