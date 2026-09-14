import { expect, test } from 'bun:test';
import { galoisinit, galoisconj4 } from './galconj.js';

test('Galois constant errors preserve native polynomial text and exception class', () => {
  for (const run of [galoisinit, galoisconj4]) {
    for (const [T, value] of [
      [[], '0'],
      [[-2n], '-2'],
    ] as const) {
      try {
        run([...T]);
        throw new Error('expected native rejection');
      } catch (e) {
        expect((e as Error).name).toBe('PariError');
        expect((e as Error).message).toBe(`not an irreducible polynomial in galoisinit: ${value}.`);
      }
    }
  }
});

test('squarefreeness precedes monicity in Galois validation', () => {
  for (const run of [galoisinit, galoisconj4]) {
    expect(() => run([2n, -4n, 2n])).toThrow('domain error in galoisinit: issquarefree(pol) = 0');
    expect(() => run([1n, 2n])).toThrow('sorry, galoisinit(nonmonic) is not yet implemented.');
  }
});

test('linear conjugates preserve native generic and cyclotomic shortcuts', () => {
  expect(galoisconj4([-2n, 1n], 0n)).toEqual([{ num: [0n, 1n], den: 1n }]);
  expect(galoisconj4([0n, 1n])).toEqual([{ num: [0n, 1n], den: 1n }]);
  expect(galoisconj4([-1n, 1n])).toEqual([{ num: [1n], den: 1n }]);
  expect(galoisconj4([1n, 1n])).toEqual([{ num: [-1n], den: 1n }]);
});

test('the non-Galois conjugate fallback returns the identity polynomial', () => {
  expect(galoisinit([-2n, 0n, 0n, 1n])).toBeNull();
  expect(galoisconj4([-2n, 0n, 0n, 1n])).toEqual([{ num: [0n, 1n], den: 1n }]);
});
