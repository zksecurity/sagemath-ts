import { expect, test } from 'bun:test';
import { QQ } from '../rational_field.js';
import type { ConstantField, ConstantFieldElement } from './constant_field.js';
import { FunctionFieldElement_rational } from './element_rational.js';
import { FunctionField, RationalFunctionField } from './index.js';

test('matrix, trace and norm run native operations on raw fractions', () => {
  const K = FunctionField(QQ as unknown as ConstantField<ConstantFieldElement>, 'x');
  const R = K._ring;
  const f = new FunctionFieldElement_rational(K, R.zero(), R.__call__([1n, 1n]), false);
  const M = f.matrix();
  expect(f.denominator().degree()).toBe(1);
  expect(M[0]![0]!.denominator().degree()).toBe(0);
  expect(f.trace().denominator().degree()).toBe(0);
  expect(f.norm().denominator().degree()).toBe(0);
  expect(f.denominator().degree()).toBe(1);
  expect(M[0]![0]).not.toBe(f);
  expect(f.trace()).not.toBe(f);
  expect(f.norm()).not.toBe(f);
  expect(f.matrix()).not.toBe(f.matrix());
});

test('function-field matrices reject writes and allow a mutable array copy', () => {
  const K = FunctionField(QQ as unknown as ConstantField<ConstantFieldElement>, 'x');
  const M = K.gen().matrix();
  expect(Object.isFrozen(M) && Object.isFrozen(M[0])).toBe(true);
  expect(() => { M[0]![0] = K.one(); }).toThrow('matrix is immutable; please change a copy instead');
  expect(() => { M[0] = []; }).toThrow('matrix is immutable; please change a copy instead');
  const C = M.map(row => row.slice());
  C[0]![0] = K.one();
  expect(String(C[0]![0])).toBe('1');
  expect(String(M[0]![0])).toBe('x');
});


test('matrix bases retain cold validation, shared warming and unhashable errors', () => {
  const K = new RationalFunctionField(QQ as never, 'matrixbase');
  const x = K.gen();
  expect(() => x.matrix(QQ)).toThrow('base must be the rational function field itself');
  expect(() => x.matrix(QQ)).toThrow('base must be the rational function field itself');
  K.one().trace();
  expect(String(x.matrix(QQ)[0]![0])).toBe('matrixbase');
  for (const [base, name] of [[[], 'list'], [{}, 'dict'], [new Map(), 'dict'],
    [Object.create(null), 'dict'], [new Set(), 'set']] as const) {
    expect(() => x.matrix(base)).toThrow(`unhashable type: '${name}'`);
  }
  const L = new RationalFunctionField(QQ as never, 'othermatrixbase');
  expect(() => L.gen().matrix(QQ)).toThrow('base must be the rational function field itself');
  expect(String(L.gen().matrix(L)[0]![0])).toBe('othermatrixbase');
});
