import { expect, test } from 'bun:test';
import { RelationTable, relationHash } from './_mpqs_hash.js';
import { mpqs } from './mpqs.js';
import { getrand, setrand } from './random.js';
import { PariError } from './errors.js';

test('MPQS preserves native ordered factors from hashed relations', () => {
  const saved = getrand();
  try {
    setrand(1n);
    expect(mpqs(1000509270188293n)).toEqual([
      [10005089n, 1n],
      [100000037n, 1n],
    ]);
  } finally {
    setrand(saved);
  }
});

test('native relation hash collisions retain distinct keys and suppress duplicates', () => {
  const a = { Y: 1n, relp: [1048578] },
    b = { Y: 2n, relp: [-403887955] };
  expect(relationHash(a)).toBe(6171988346495030103n);
  expect(relationHash(b)).toBe(relationHash(a));
  const table = new RelationTable(0);
  table.add(a);
  table.add(b);
  table.add({ Y: 1n, relp: [1048578] });
  expect(table.size).toBe(2);
  expect([...table.values()].map((r) => r.Y)).toEqual([2n, 1n]);
});

test('native relation table grows before insertion 36 and reorders its buckets', () => {
  const table = new RelationTable(0);
  for (let i = 0; i < 35; i++) table.add({ Y: BigInt(i), relp: [1048578] });
  expect([...table.values()].map((r) => Number(r.Y))).toEqual([
    20, 2, 19, 1, 0, 18, 17, 34, 16, 33, 15, 32, 14, 31, 13, 30, 12, 29, 11, 28, 10, 27, 9, 26, 8,
    25, 7, 24, 6, 23, 5, 22, 4, 21, 3,
  ]);
  table.add({ Y: 35n, relp: [1048578] });
  expect([...table.values()].map((r) => Number(r.Y))).toEqual([
    16, 23, 30, 7, 14, 21, 28, 35, 5, 12, 19, 26, 33, 3, 10, 17, 24, 31, 1, 8, 15, 22, 29, 0, 6, 13,
    20, 27, 34, 4, 11, 18, 25, 32, 2, 9,
  ]);
});

test('native relation table size overflow retains the PARI diagnostic', () => {
  expect(() => new RelationTable(1610612741)).toThrow(PariError);
  expect(() => new RelationTable(1610612741)).toThrow('overflow in hash table [too large].');
});
