import { expect, test } from 'bun:test';
import {
  perm_powu,
  perm_cycles,
  cyc_pow,
  perm_orderu,
  vecsmall_lexcmp,
  vecsmall_uniq,
  group_order,
  group_domain,
  trivialgroup,
} from './galconj.js';
import { zv_prod } from './ZV.js';
import { PariError } from './errors.js';

test('unsigned permutation exponents use checked native magnitudes', () => {
  const p = [0, 2, 3, 1];
  expect(perm_powu(p, -1)).toEqual(p);
  expect(perm_powu(p, -2)).toEqual([0, 3, 1, 2]);
  expect(() => perm_powu([0], 2 ** 64)).toThrow(
    new PariError('overflow in t_INT-->ulong assignment.')
  );
});
test('signed cycle exponents validate before the empty-cycle shortcut', () => {
  expect(cyc_pow(perm_cycles([0, 2, 3, 1]), -1)).toEqual([[], [0, 1, 3, 2]]);
  for (const e of [2 ** 63, -(2 ** 63)])
    expect(() => cyc_pow([[]], e)).toThrow(new PariError('overflow in t_INT-->long assignment.'));
});
test('permutation order keeps exact unsigned intermediates before Number projection', () => {
  const p = [0];
  for (const n of [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53]) {
    const start = p.length;
    p.push(...Array.from({ length: n - 1 }, (_, i) => start + i + 1), start);
  }
  expect(perm_orderu(p)).toBe(Number(14142414403480493114n));
});
test('lexicographic prefix comparison returns only the native sign', () => {
  expect(vecsmall_lexcmp([0], [0, 1, 2, 3])).toBe(-1);
  expect(vecsmall_lexcmp([0, 1, 2, 3], [0])).toBe(1);
  expect(vecsmall_lexcmp([0, 1, 2], [0, 1, 2])).toBe(0);
});
test('signed-word group products round only at the public return', () => {
  const v = [-9, -14, 3, 15, -15, 10, -17, 11, -14, 11, -5, -13, 12, 9, 17, 15, -2, 2, -9, -3];
  expect(zv_prod(v)).toBe(Number(4735193091829200000n));
  expect(group_order({ gen: [[]], ord: [0, ...v] })).toBe(Number(4735193091829200000n));
  expect(zv_prod([-(2 ** 62), 2])).toBe(-(2 ** 63));
  expect(() => group_order({ gen: [[]], ord: [0, 2 ** 32, 2 ** 32, 0] })).toThrow(
    new RangeError('zv_prod requires a product within the signed-word range')
  );
});
test('trivial-group domain retains the native error class', () => {
  expect(() => group_domain(trivialgroup())).toThrow(
    new PariError('domain error in group_domain: #G = 1')
  );
});
test('unique integer vectors cover counting and comparison-sort paths', () => {
  expect(vecsmall_uniq([0, 3, 1, 3, 2, 0])).toEqual([0, 0, 1, 2, 3]);
  expect(vecsmall_uniq([0, 100, -1, 100, 3])).toEqual([0, -1, 3, 100]);
  expect(vecsmall_uniq([0, 0, 0])).toEqual([0, 0]);
});
