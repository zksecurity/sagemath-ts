import { expect, spyOn, test } from 'bun:test';
import * as P from './ff.js';
import * as factorization from './ifactor.js';
import { PariError } from './errors.js';

// All boundaries have permanent native counterparts in ff_pari_fp_predicate.
test('Fp_issquare uses the native Kronecker predicate at composite and zero moduli', () => {
  expect(P.Fp_issquare(2n, 6n)).toBe(true);
  expect(P.Fp_issquare(3n, 0n)).toBe(true);
  expect(P.Fp_issquare(-31n, -32n)).toBe(false);
  expect(P.Fp_issquare(-3n, -2n)).toBe(true);
});

test('Bezout chooses both coefficients zero for zero arguments', () => {
  expect(P.xgcd(0n, 0n)).toEqual([0n, 0n, 0n]);
  expect(P.xgcd(0n, -17n)).toEqual([17n, 0n, -1n]);
});

test('word order bounds normalize signs and default at zero or overflow', () => {
  const factor = spyOn(factorization, 'factoru');
  try {
    expect(P.Fp_order(2n, 0n, 17n)).toBe(8n);
    expect(P.Fp_order(2n, -16n, -17n)).toBe(8n);
    expect(P.Fp_order(2n, 1n << 64n, 17n)).toBe(8n);
    expect(factor.mock.calls.map((c) => c[0])).toEqual([16n, 16n, 16n]);
  } finally {
    factor.mockRestore();
  }
});

test('generic order validates before its identity shortcut and retains raw bases', () => {
  const p = (1n << 64n) + 13n;
  expect(() => P.Fp_order(1n, 0n, p)).toThrow(
    new PariError('incorrect type in generic discrete logarithm (order factorization) (t_INT).')
  );
  expect(P.Fp_order(p + 1n, 2n, p)).toBe(2n);
  expect(P.Fp_order(1n, 2n, 0n)).toBe(1n);
  expect(P.Fp_order(0n, 1n, 0n)).toBe(1n);
});

test('znorder preserves native modular coercion and full coprimality errors', () => {
  expect(P.znorder(2n, -17n, 0n)).toBe(8n);
  expect(() => P.znorder(1n, 0n)).toThrow(new PariError('impossible inverse in %: 0.'));
  expect(() => P.znorder(-2n, 6n)).toThrow(
    new PariError('elements not coprime in znorder:\n    4\n    6')
  );
  expect(() => P.znorder(2n, 10n ** 2100n)).toThrow(
    new PariError(
      'elements not coprime in znorder:\n    2\n    \n  ***  (...) Huge t_INT omitted; you can access it via dbg_err()'
    )
  );
});
