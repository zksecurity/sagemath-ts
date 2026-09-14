import { expect, test } from 'bun:test';
import { mpqsInternals as P } from './mpqs.js';
import { PariError } from './errors.js';

function handle(n: bigint) {
  const h = P.newHandle();
  h.N = n;
  h.size_of_FB = 6;
  P.mpqs_FB_ctor(h).p.set([0, 0, 2, 3, 5, 7, 11, 13]);
  return h;
}
function nativeError(run: () => unknown, message: string) {
  let caught: unknown;
  try {
    run();
  } catch (error) {
    caught = error;
  }
  expect(caught).toBeInstanceOf(PariError);
  expect((caught as Error).message).toBe(message);
}

test('MPQS Kronecker symbols reduce signed inputs with native residues', () => {
  expect(P.kroiu(-100n, 3)).toBe(-1);
  expect(P.kroiu(-100n, 9)).toBe(1);
  expect(P.kroiu(-1n, 16)).toBe(1);
});

test('MPQS word symbols preserve bits beyond signed 32-bit shifts', () => {
  expect(P.krouu(2147483648, 2147483649)).toBe(1);
  expect(P.krouu(4294967294, 4294967295)).toBe(-1);
  expect(P.krouu(1099511627782, 1099511627783)).toBe(-1);
  expect(P.krouu(9007199254740990, 9007199254740991)).toBe(-1);
});

test('MPQS symbols terminate at word denominators with zero low limbs', () => {
  for (const k of [31, 32, 33, 40, 52]) {
    expect(P.krouu(1, 2 ** k)).toBe(1);
    expect(P.kroiu(1n, 2 ** k)).toBe(1);
  }
});

test('MPQS relation reconstruction supports signed moduli and inverse powers', () => {
  expect(P.mpqs_factorback(handle(-101n), [1048577])).toBe(100n);
  expect(P.mpqs_factorback(handle(-101n), [-1048574])).toBe(51n);
  expect(P.mpqs_factorback(handle(23n), [-1048574])).toBe(12n);
});

test('MPQS inverse relation powers retain the native inverse error', () => {
  nativeError(
    () => P.mpqs_factorback(handle(2n), [-1048574]),
    'impossible inverse in Fl_inv: Mod(0, 2).'
  );
});

test('MPQS debug relation failures retain native error type and text', () => {
  for (const q of [0, 1, 3]) {
    nativeError(
      () => P.mpqs_check_rel(handle(23n), { Y: 2n, relp: [1048578] }, q, 0),
      `bug in MPQS: wrong ${q ? 'large prime' : 'full'} relation found, please report.`
    );
  }
});

test('MPQS relation checker keeps valid and class-group early-return paths', () => {
  expect(P.mpqs_check_rel(handle(23n), { Y: 5n, relp: [1048578] }, 1, 0)).toBeUndefined();
  expect(P.mpqs_check_rel(handle(2n), { Y: 0n, relp: [-1048574] }, 1, 1)).toBeUndefined();
  expect(P.mpqs_check_rel(handle(2n), { Y: 1n, relp: [-1048574] }, 3, 1)).toBeUndefined();
});
