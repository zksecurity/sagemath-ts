import { expect, test } from 'bun:test';
import { QQ } from '../rings/rational_field.js';
import { RDF } from '../rings/real_double.js';
import { PolynomialRing } from '../rings/polynomial/polynomial_ring.js';
import { integer_to_real_double_dense } from './change_ring.js';
import { IntegerMatrixFromEntries } from './matrix_integer.js';

test('matrix evaluation uses scalar diagonal addition and preserves original input', () => {
  const point = IntegerMatrixFromEntries([
    [1n, 2n],
    [3n, 4n],
  ]);
  const f = new PolynomialRing(QQ, 'x').__call__([1n, 0n, 1n]);
  const value = f.evaluate(point);
  // Sage: (x^2+1)(matrix(ZZ,[[1,2],[3,4]])) == matrix(QQ,[[8,10],[15,23]])
  expect([value.get(0, 0), value.get(0, 1), value.get(1, 0), value.get(1, 1)].map(String)).toEqual([
    '8',
    '10',
    '15',
    '23',
  ]);
  expect(point.get(0, 0).value).toBe(1n);
});

test('native matrix conversion and scalar RDF construction use different rounding rules', () => {
  const value = 9007199254740995n;
  const h = integer_to_real_double_dense(IntegerMatrixFromEntries([[value]]));
  expect(h.get(0, 0).value).toBe(9007199254740994);
  expect(RDF.__call__(value).value).toBe(9007199254740996);
});
