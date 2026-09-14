import { expect, test } from 'bun:test';
import { NumberField, NumberFieldElement, RationalPolynomial } from './number_field.js';
import { Rational } from '../rational.js';
import { Integer } from '../integer_ring.js';
const field = (cs: bigint[]) => new NumberField(RationalPolynomial.fromBigInts(cs), 'a');

test('element construction reduces coefficients modulo the defining polynomial', () => {
  const K = field([-2n, 0n, 1n]);
  const a = new NumberFieldElement(
    K,
    [1n, 2n, 3n].map((c) => new Rational(c))
  );
  expect(a.list().map(String)).toEqual(['7', '2']);
  expect(new NumberFieldElement(K, []).list().map(String)).toEqual(['0', '0']);
  expect(new NumberFieldElement(K, [Rational.one()]).list().map(String)).toEqual(['1', '0']);
});

test('nonmonic degree-one generators and reduction retain their rational root', () => {
  const K = field([-1n, 2n]);
  expect(K.gen().list().map(String)).toEqual(['1/2']);
  expect(
    new NumberFieldElement(
      K,
      [1n, 2n, 3n].map((c) => new Rational(c))
    )
      .list()
      .map(String)
  ).toEqual(['11/4']);
  expect(field([-3n, 1n]).gen().list().map(String)).toEqual(['3']);
});

test('higher-degree coefficient reduction preserves exact large integers', () => {
  const K = field([-2n, 0n, 0n, 1n]);
  const a = new NumberFieldElement(
    K,
    [2n ** 80n + 1n, 0n, 0n, 1n].map((c) => new Rational(c, 3n))
  );
  expect(a.list().map(String)).toEqual([new Rational(2n ** 80n + 3n, 3n).toString(), '0', '0']);
});

test('field polynomial getters and display preserve original rational scaling', () => {
  const f = new RationalPolynomial([-2n, 0n, 1n].map((c) => new Rational(c, 3n)));
  const K = new NumberField(f, 'a');
  expect(K.polynomial()).toBe(f);
  expect(K.defining_polynomial()).toBe(f);
  expect(String(K)).toBe('Number Field in a with defining polynomial 1/3*x^2 - 2/3');
  expect(K.gen().pow(2n).eq(2n)).toBe(true);
});

test('field integer powers preserve Sage self aliasing and Integer coercion', () => {
  const K = field([-2n, 0n, 1n]),
    a = K.gen();
  expect(a.pow(1n)).toBe(a);
  expect(a.pow(new Integer(1n))).toBe(a);
  expect(a.pow(new Integer(3n)).list().map(String)).toEqual(['0', '2']);
  expect(a.pow(0n).is_one()).toBe(true);
  expect(a.pow(-1n).list().map(String)).toEqual(['0', '1/2']);
  expect(
    K.one()
      .pow(2n ** 100n + 7n)
      .is_one()
  ).toBe(true);
  expect(() => K.zero().pow(-1n)).toThrow('number field element division by zero');
});

test('scaled quadratic fields retain fundamental units, logs and regulators', () => {
  const K = new NumberField(
    new RationalPolynomial([new Rational(-1n, 6n), Rational.zero(), new Rational(1n, 2n)]),
    'a'
  );
  const U = K.unit_group();
  expect(U.fundamental_units()[0]!.list().map(String)).toEqual(['2', '3']);
  expect(U.log(U.exp([1n, -8n]))).toEqual([1n, -8n]);
  expect(Math.round(K.regulator() * 10 ** 8)).toBe(131695790);
  expect(field([2n, 0n, 2n]).unit_group().torsion_order()).toBe(4n);
});

test('scaled monogenic maximal orders retain different and codifferent norms', () => {
  const O = field([-4n, 0n, 2n]).maximal_order() as import('./order.js').AbsoluteOrder;
  expect(String(O.different().norm())).toBe('8');
  expect(String(O.codifferent().norm())).toBe('1/8');
});

test('Frobenius filtering accepts scaled presentations of the same integral generator', async () => {
  const { Frobenius_filter } = await import('../../schemes/elliptic_curves/isogeny_class.js');
  const K = field([2n, 0n, 2n]),
    a = K.gen();
  const E = {
    ainvs: () => [a.add(1n), a.neg(), a, a.mul(-240n).sub(399n), a.mul(2869n).add(2627n)],
    base_field: () => K,
  };
  expect(Frobenius_filter(E as never, [2n, 3n, 5n, 7n, 11n, 13n, 17n, 19n])).toEqual([2n, 3n]);
});
