import { expect, test } from 'bun:test';
import { mpqsInternals as P } from './mpqs.js';
import { PariError } from './errors.js';

test('MPQS word inverse failures preserve native PariError diagnostics', () => {
  for (const [a, p] of [
    [0, 2],
    [2, 4],
    [3, 9],
  ]) {
    let caught: unknown;
    try {
      P.Fl_inv(a!, p!);
    } catch (error) {
      caught = error;
    }
    expect(caught).toBeInstanceOf(PariError);
    expect((caught as Error).message).toBe(`impossible inverse in Fl_inv: Mod(${a}, ${p}).`);
  }
});

test('MPQS inverse retains native unit and modulus-one controls', () => {
  expect(P.Fl_inv(3, 97)).toBe(65);
  expect(P.Fl_inv(3, 100)).toBe(67);
  expect(P.Fl_inv(0, 1)).toBe(0);
});

test('MPQS factor-base construction returns factors and skips square divisors', () => {
  const h = P.newHandle();
  h.N = h.kN = 77n;
  h.size_of_FB = 6;
  h.index0_FB = 3;
  h.pmin_index1 = 3;
  expect(P.mpqs_create_FB(h, true)).toBe(7);
  h.N = h.kN = -36n;
  expect(P.mpqs_create_FB(h, false)).toBe(0);
  expect([...h.FB.p]).not.toContain(3);
});

test('MPQS candidate cap finishes the current native eight-byte block', () => {
  const h = P.newHandle();
  h.M = 2000;
  h.sieve_threshold = 128;
  P.mpqs_sieve_array_ctor(h);
  for (let i = 0; i < 2 * h.M; i++) h.sieve_array[i] = i % 8 === 7 ? 127 : 255;
  const n = P.mpqs_eval_sieve(h);
  expect(n).toBe(2002);
  expect(h.candidates[n - 1]).toBe(2286);
  expect(h.candidates[n]).toBe(0);
});
