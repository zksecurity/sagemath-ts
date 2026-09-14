import { expect, test } from 'bun:test';
import { FpV_invVandermonde } from './galconj.js';
import { FpV_inv, FpV_invVandermonde as dependencyVandermonde } from './FpX.js';
import { producttree_scheme, gen_product } from './bb_group.js';
import { PariError } from './errors.js';

test('Vandermonde facade retains one-indexed rows and columns', () => {
  expect(FpV_invVandermonde([0n, 0n, 1n], 1n, 5n)).toEqual([
    [0n, 0n, 0n],
    [0n, 1n, 0n],
    [0n, 4n, 1n],
  ]);
  expect(FpV_invVandermonde([0n, 0n, 1n], 2n, 5n)).toEqual([
    [0n, 0n, 0n],
    [0n, 2n, 0n],
    [0n, 3n, 2n],
  ]);
  expect(FpV_invVandermonde([0n], 0n, 0n)).toEqual([[0n]]);
  expect(FpV_invVandermonde([], 0n, 0n)).toEqual([]);
  expect(() => dependencyVandermonde([], 0n, 0n)).toThrow(
    new RangeError('FpV_invVandermonde requires a nonempty vector')
  );
});
test('singular Vandermonde vectors invert before multiplying a zero denominator', () => {
  for (const den of [0n, 1n, 5n])
    expect(() => FpV_invVandermonde([0n, 1n, 1n], den, 5n)).toThrow(
      new PariError('impossible inverse in Fp_inv: Mod(5, 5).')
    );
});
test('invalid modulus errors originate in native product/evaluation phases', () => {
  expect(() => FpV_invVandermonde([0n, 0n], 0n, 0n)).toThrow(
    new PariError('impossible inverse in dvmdii: 0.')
  );
});
test('batch inverse returns native residues and reports the product gcd', () => {
  expect(FpV_inv([2n, 3n, 4n], 7n)).toEqual([4n, 5n, 2n]);
  expect(() => FpV_inv([2n, 3n], 6n)).toThrow(
    new PariError('impossible inverse in Fp_inv: Mod(6, 6).')
  );
  expect(() => FpV_inv([], 7n)).toThrow(new RangeError('FpV_inv requires a nonempty vector'));
});
test('shared product scheme preserves odd splits and callback order', () => {
  expect(producttree_scheme(5)).toEqual([2, 1, 1, 1]);
  expect(producttree_scheme(0)).toEqual([0]);
  const calls: string[] = [];
  expect(
    gen_product(['a', 'b', 'c', 'd', 'e'], (a, b) => {
      const result = '(' + a + b + ')';
      calls.push(result);
      return result;
    })
  ).toBe('(((ab)c)(de))');
  expect(calls).toEqual(['(ab)', '((ab)c)', '(de)', '(((ab)c)(de))']);
});
