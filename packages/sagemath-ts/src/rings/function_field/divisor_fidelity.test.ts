import { expect, test } from 'bun:test';
import { QQ } from '../rational_field.js';
import { GF } from '../finite_rings/index.js';
import { FunctionField, DivisorGroup, divisor } from './index.js';
import type { ConstantField, ConstantFieldElement } from './constant_field.js';

const K = FunctionField(QQ as unknown as ConstantField<ConstantFieldElement>, 'x');

test('divisor parents, place conversion and zero support follow native coercion', () => {
  const G = K.divisor_group(), P = K.maximal_order().ideal(K.gen()).place();
  const D = P.divisor(0n);
  expect(new DivisorGroup(K)).toBe(G);
  expect(G.__call__(D)).toBe(D);
  expect(G.zero()).toBe(G.zero());
  expect(String(G.__call__(P))).toBe('Place (x)');
  expect(D.support()).toHaveLength(1);
  expect(D.add(G.zero()).support()).toHaveLength(0);
  expect(D._format((v) => `<${v}>`, ' @ ', ' / ')).toBe('<0> @ <Place (x)>');
});

test('public function-space maps are cached and retain native vector conversion errors', () => {
  const D = divisor(K, []), A = D.function_space();
  expect(D.function_space()).toBe(A);
  expect(() => A[1]([])).toThrow("[] fails to convert into the map's domain Vector space of dimension 1 over Rational Field, but a `pushforward` method is not properly implemented");
  expect(() => A[2](K.gen())).toThrow('(0, 1) is not in list');
  const P = K.maximal_order_infinite().ideal(K.gen().inv()).place();
  expect(() => P.divisor(-1n).function_space()[1]([])).toThrow('Cannot convert int to sage.structure.element.Element');
});

test('dependent and zero echelon rows retain the native pivot protocol', () => {
  const D = K.divisor_group().zero();
  const [basis, coordinates] = D._echelon_basis([K.one(), K.one()]);
  expect(basis.map(String)).toEqual(['1', '0']);
  expect(coordinates(K.one()).map(String)).toEqual(['1', '0']);
  expect(() => D._echelon_basis([K.zero(), K.one()])).toThrow("'NoneType' object is not subscriptable");
});

test('generic dense vector empty lists expand at the exact native modulus cutoff', () => {
  const F = FunctionField(GF(2147483647n) as unknown as ConstantField<ConstantFieldElement>, 'x');
  expect(String(F.divisor_group().zero().function_space()[1]([]))).toBe('0');
});
