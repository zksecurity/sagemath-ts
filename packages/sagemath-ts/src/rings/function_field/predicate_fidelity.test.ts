import { expect, test } from 'bun:test';
import { AttributeError, ValueError } from '../../errors.js';
import { Matrix } from '../../matrix/matrix_generic.js';
import { QQ } from '../rational_field.js';
import { FunctionField, is_FunctionField, is_FunctionFieldElement } from './index.js';

test('element predicates preserve both parent lookups and category membership', () => {
  const K = FunctionField(QQ as never, 'x');
  let calls = 0;
  const element = { parent: () => (++calls === 1 ? QQ : K) };
  expect(is_FunctionFieldElement(element)).toBe(true);
  expect(calls).toBe(2);
  expect(is_FunctionFieldElement(1n)).toBe(false);
  expect(is_FunctionFieldElement(new Matrix(QQ, 1, 1, [[QQ.one()]]))).toBe(false);
  expect(() => is_FunctionFieldElement({})).toThrow("'dict' object has no attribute 'parent'");
});

test('category predicates read dynamic attributes once and preserve error boundaries', () => {
  let reads = 0;
  const field = new Proxy(
    {},
    {
      get: (_target, key) => {
        if (key !== 'category') return undefined;
        reads++;
        return () => ({
          is_subcategory: (category: unknown) => String(category) === 'Category of function fields',
        });
      },
    }
  );
  expect(is_FunctionField(field)).toBe(true);
  expect(reads).toBe(1);
  expect(
    is_FunctionField({
      category() {
        throw new AttributeError('missing');
      },
    })
  ).toBe(false);
  expect(() =>
    is_FunctionField({
      category() {
        throw new ValueError('broken');
      },
    })
  ).toThrow('broken');
  expect(() => is_FunctionField({ category: () => ({}) })).toThrow(
    "'dict' object has no attribute 'is_subcategory'"
  );
});
