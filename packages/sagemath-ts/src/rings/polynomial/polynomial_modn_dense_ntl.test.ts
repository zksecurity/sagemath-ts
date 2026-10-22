import { expect, test } from 'bun:test';
import { next_prime } from '../../arith/misc.js';
import { ArithmeticError, AttributeError, ValueError } from '../../errors.js';
import { GFExtended as GF } from '../finite_rings/finite_field_extension.js';
import type { IntegerMod } from '../finite_rings/integer_mod.js';
import { IntegerModRing } from '../finite_rings/integer_mod_ring.js';
import type { CoefficientRing, Polynomial } from './polynomial_element.js';
import { type ModnCoefficient, small_roots } from './polynomial_modn_dense_ntl.js';
import { PolynomialRing } from './polynomial_ring.js';

const Zmod = (n: bigint) => new IntegerModRing(n);
const ring = (base: unknown, name = 'x') =>
  new PolynomialRing(base as CoefficientRing<ModnCoefficient>, name);

// Expected values are Sage 10.3's, from the small_roots docstring where possible.
test('small_roots finds the docstring root and returns roots in Python set order', () => {
  const R = new PolynomialRing(Zmod(10001n), 'x');
  const f = R.__call__([-222n, 5000n, 10n, 1n]);
  expect(f.small_roots().map(String)).toEqual(['4']);
  expect(small_roots(f).map(String)).toEqual(['4']);
  expect(f.small_roots({ X: 3n })).toEqual([]);
  expect(R.__call__([-1n, 0n, 1n]).small_roots().map(String)).toEqual(['10000', '1']);
  const S = ring(GF(10007n), 'y');
  expect(S.__call__([-4n, 0n, 1n]).small_roots().map(String)).toEqual(['2', '10005']);
});

test('small_roots recovers a stereotyped RSA message (e = 3, 56 unknown bits)', () => {
  const p = 2n ** 256n + 2n ** 8n + 2n ** 5n + 2n ** 3n + 1n;
  const q = 2n ** 256n + 2n ** 8n + 2n ** 5n + 2n ** 3n + 2n ** 2n + 1n;
  const N = p * q;
  const K = 0xc0ffee123456n;
  const padding = 2n ** 512n - 2n ** 56n;
  const R = new PolynomialRing(Zmod(N), 'x');
  const C = (padding + K) ** 3n % N;
  const f = (R.__call__([padding, 1n]).pow(3n) as Polynomial<IntegerMod>).sub(R.__call__([C]));
  expect(f.small_roots().map((r) => r.value)).toEqual([K]);
});

test('small_roots factors N with known high bits of q (beta = 0.5)', () => {
  const p = next_prime(2n ** 256n);
  const q = next_prime(3n * p + p / 7n);
  const qbar = q + 123456789n;
  const R = new PolynomialRing(Zmod(p * q), 'x');
  const d = R.__call__([-qbar, 1n]).small_roots({ X: 2n ** 110n - 1n, beta: 0.5 });
  expect(d.map((r) => qbar - r.value)).toEqual([q]);
});

test('small_roots keeps Sage errors', () => {
  const R = new PolynomialRing(Zmod(10001n), 'x');
  const f = R.__call__([-222n, 5000n, 10n, 1n]);
  expect(() => R.__call__([1n, 2n]).small_roots()).toThrow(
    new ArithmeticError('Polynomial must be monic.')
  );
  expect(() => f.small_roots({ beta: 0 })).toThrow(
    new ValueError('0.0 < beta <= 1.0 not satisfied.')
  );
  expect(() => R.one().small_roots()).toThrow(new ValueError('Calling ceil() on infinity or NaN'));
  expect(() => f.small_roots({ X: 0n })).toThrow('Integer division by zero');
  expect(() => ring(GF(2n)).gen().small_roots()).toThrow(
    new AttributeError(
      "'sage.rings.polynomial.polynomial_gf2x.Polynomial_GF2X' object has no attribute 'small_roots'"
    )
  );
});
