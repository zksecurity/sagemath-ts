import { expect, test } from 'bun:test';
import { idealmul } from '@sagemath-ts/parigp-ts/src/base4.js';
import { Integer } from '../integer_ring.js';
import { Rational } from '../rational.js';
import { NumberField, RationalPolynomial } from './number_field.js';

const field = (c = 2n, name = 'a') =>
  new NumberField(RationalPolynomial.fromBigInts([c, 0n, 1n]), name);

test('membership coerces scalars and returns false for incompatible foreign elements', () => {
  const K = field(),
    L = field(3n, 'b'),
    I = K.ideal(1n);
  expect(I.contains(2n)).toBe(true);
  expect(I.contains(new Rational(1n, 2n))).toBe(false);
  expect(I.contains([1n, 2n])).toBe(true);
  expect(I.contains(L.gen())).toBe(false);
  expect(I.contains(K.ideal(1n))).toBe(false);
  expect(() => I.contains([])).toThrow('Length must be equal to the degree of this number field');
  expect(() => I.contains([Infinity, 0n])).toThrow(
    'cannot convert NaN or infinity to rational number'
  );
});

test('division preserves float-zero and native-int fallback distinctions', () => {
  const K = field(),
    I = K.ideal(2n),
    zero = K.ideal(0n);
  expect(I.div(new Rational(1n, 2n)).norm().toString()).toBe('16');
  expect(() => I.div(0)).toThrow('number field element division by zero');
  expect(() => I.div(2)).toThrow(
    "unsupported operand type(s) for /: 'NumberFieldFractionalIdeal' and 'float'"
  );
  expect(() => zero.div(2n)).toThrow(
    "unsupported operand type(s) for /: 'NumberFieldIdeal' and 'int'"
  );
  expect(() => zero.div(new Integer(2n))).toThrow('[;] has unsupported PARI type t_MAT');
});

test('foreign ideal arithmetic keeps operation-specific native coercion', () => {
  const K = field(),
    L = field(3n, 'b'),
    I = K.ideal(3n, K.gen().add(1n));
  expect(I.add(L.ideal(2n)).norm().toString()).toBe('1');
  expect(I.is_coprime(L.ideal(2n))).toBe(true);
  expect(() => I.divides(L.ideal(2n))).toThrow('unsupported operand parent(s) for /');
  // Sage's multi-generator route interprets the RHS HNF in the LHS basis.
  expect(
    I.mul(L.ideal(L.gen().add(new Rational(3n, 2n))))
      .norm()
      .toString()
  ).toBe('9/4');
});

test('rational idealmul removes primitive content and preserves input matrices', () => {
  const nf = field(-2n)._pari_ideal_data();
  const I = [
      [12n, 0n],
      [0n, 12n],
    ],
    J = [
      [35n, 0n],
      [0n, 35n],
    ],
    original = structuredClone([I, J]);
  expect(idealmul(nf, I, 6n, J, 5n)).toEqual([
    [
      [14n, 0n],
      [0n, 14n],
    ],
    1n,
  ]);
  expect([I, J]).toEqual(original);
  expect(idealmul(nf, [], 1n, J, 5n)).toEqual([[], 1n]);
  expect(() => idealmul(nf, I, 0n, J, 1n)).toThrow(RangeError);
  expect(() => idealmul(nf, I, 1n, [[1n]], 1n)).toThrow('inconsistent dimensions in idealmul');
});
