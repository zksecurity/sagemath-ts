import { expect, test } from 'bun:test';
import { Integer } from '../integer_ring.js';
import { NumberField, RationalPolynomial } from './number_field.js';
const field = (coeffs: bigint[]) => new NumberField(RationalPolynomial.fromBigInts(coeffs), 'a');

test('square-root coefficient tuples accept negative indices', () => {
  const a = field([-2n, 0n, 1n]).gen().add(2n);
  expect(String(a.__getitem__(-2))).toBe('2');
  expect(String(a.__getitem__(-1))).toBe('1');
  expect(() => a.__getitem__(2)).toThrow('index must be either 0 or 1');
});

test('general field coefficient indices enforce degree bounds', () => {
  const a = field([-2n, 0n, 0n, 1n]).gen();
  expect(() => a.__getitem__(-1)).toThrow('index must be between 0 and degree minus 1');
  expect(() => a.__getitem__(3)).toThrow('index must be between 0 and degree minus 1');
  expect(String(a.__getitem__(2))).toBe('0');
});

test('shifted quadratic generators use general coefficient indexing', () => {
  const a = field([-1n, -1n, 1n]).gen();
  expect(() => a.__getitem__(-1)).toThrow('index must be between 0 and degree minus 1');
});

test('a square factor changes the native quadratic representation', () => {
  const a = field([-8n, 0n, 1n]).gen();
  expect(() => a.__getitem__(-1)).toThrow('index must be between 0 and degree minus 1');
});

test('quadratic representation uses the native bounded square removal', () => {
  const a = field([-2n * 10007n ** 2n, 0n, 1n]).gen();
  const b = field([-2n * 9973n ** 2n, 0n, 1n]).gen();
  expect(String(a.__getitem__(-1))).toBe('1');
  expect(() => b.__getitem__(-1)).toThrow('index must be between 0 and degree minus 1');
});

test('original leading coefficient affects quadratic indexing', () => {
  const a = field([-4n, 0n, 2n]).gen();
  const b = field([4n, 0n, -2n]).gen();
  expect(String(a.__getitem__(-1))).toBe('1');
  expect(() => b.__getitem__(-1)).toThrow('index must be between 0 and degree minus 1');
});

test('fractional indices follow native trimmed-list and tuple branches', () => {
  const K = field([-2n, 0n, 0n, 1n]);
  expect(String(K.zero().__getitem__(0.5))).toBe('0');
  expect(() => K.gen().__getitem__(0.5)).toThrow(
    'list indices must be integers or slices, not float'
  );
  expect(() => field([-2n, 0n, 1n]).zero().__getitem__(0.5)).toThrow(
    'tuple indices must be integers or slices, not float'
  );
});

test('coefficient indexing accepts exact integer index types', () => {
  const a = field([-2n, 0n, 1n]).gen();
  expect(String(a.__getitem__(1n))).toBe('1');
  expect(String(a.__getitem__(new Integer(-1n)))).toBe('1');
  expect(() => a.__getitem__(2n ** 80n)).toThrow('index must be either 0 or 1');
});
