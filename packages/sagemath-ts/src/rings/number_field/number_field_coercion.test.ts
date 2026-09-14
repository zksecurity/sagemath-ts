import { expect, test } from 'bun:test';
import { Integer } from '../integer_ring.js';
import { Rational } from '../rational.js';
import { NumberField, QuadraticField, RationalPolynomial } from './number_field.js';
const quadratic = () => QuadraticField.create(2n, 'a');
const cubic = () => new NumberField(RationalPolynomial.fromBigInts([-2n, 0n, 0n, 1n]), 'a');

test('field addition and subtraction coerce exact scalars', () => {
  const a = quadratic().gen();
  expect(a.add(3n).list().map(String)).toEqual(['3', '1']);
  expect(a.sub(new Integer(3n)).list().map(String)).toEqual(['-3', '1']);
  expect(a.add(new Rational(2n, 3n)).list().map(String)).toEqual(['2/3', '1']);
});

test('field division accepts scalar denominators', () => {
  expect(quadratic().gen().div(3n).list().map(String)).toEqual(['0', '1/3']);
  expect(cubic().gen().div(new Rational(2n, 3n)).list().map(String)).toEqual(['0', '3/2', '0']);
});

test('field equality coerces exact scalars', () => {
  const K = quadratic();
  expect(K.one().eq(1n)).toBe(true);
  expect(K.__call__(new Rational(2n, 3n)).eq(new Rational(2n, 3n))).toBe(true);
  expect(K.gen().eq(new Integer(1n))).toBe(false);
});

test('field constructors accept Integer wrappers', () => {
  expect(quadratic().__call__(new Integer(-3n)).list().map(String)).toEqual(['-3', '0']);
});

test('field constructors preserve rational values of real inputs', () => {
  const K = quadratic();
  expect(K.__call__(1.5).list().map(String)).toEqual(['3/2', '0']);
  expect(K.__call__(-0.1).list().map(String)).toEqual(['-1/10', '0']);
  expect(
    K.__call__(2 ** 100)
      .list()
      .map(String)
  ).toEqual(['1267650600228229331127959027712', '0']);
});

test('field construction rejects nonfinite reals with native diagnostics', () => {
  expect(() => cubic().__call__(Infinity)).toThrow(
    'unable to convert inf to Number Field in a with defining polynomial x^3 - 2'
  );
  expect(() => cubic().__call__(NaN)).toThrow(
    'unable to convert nan to Number Field in a with defining polynomial x^3 - 2'
  );
});

test('zero field inverses and negative powers preserve native errors', () => {
  for (const K of [quadratic(), cubic()]) {
    expect(() => K.zero().inv()).toThrow('number field element division by zero');
    expect(() => K.zero().pow(-1n)).toThrow('number field element division by zero');
  }
});

test('scalar and field zero division follow the native coercion route', () => {
  const K = quadratic(),
    L = cubic();
  expect(() => K.gen().div(0n)).toThrow('rational division by zero');
  expect(() => K.gen().div(K.zero())).toThrow('number field element division by zero');
  expect(() => L.gen().div(0n)).toThrow('number field element division by zero');
});
