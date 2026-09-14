import { expect, test } from 'bun:test';
import { eulerphiu as wordTotient } from './arith2.js';
import { PariError } from './errors.js';
import { eulerphiu } from './galconj.js';

test('Galois word totient returns native two at zero', () => {
  expect(eulerphiu(0)).toBe(2);
  expect(wordTotient(0n)).toBe(2n);
});

test('Galois word totient uses native word magnitudes', () => {
  expect(eulerphiu(-1)).toBe(1);
  expect(eulerphiu(-12)).toBe(4);
  expect(eulerphiu(-32)).toBe(16);
});

test('word totients retain all 64 bits before facade rounding', () => {
  expect(wordTotient((1n << 64n) - 1n)).toBe(9208981628670443520n);
  expect(wordTotient(1n << 63n)).toBe(1n << 62n);
  expect(eulerphiu(2 ** 63)).toBe(2 ** 62);
});

test('word conversion rejects overflow before totient evaluation', () => {
  let caught: unknown;
  try {
    eulerphiu(2 ** 64);
  } catch (e) {
    caught = e;
  }
  expect(caught?.constructor).toBe(PariError);
  expect((caught as Error).message).toBe('overflow in t_INT-->ulong assignment.');
});
