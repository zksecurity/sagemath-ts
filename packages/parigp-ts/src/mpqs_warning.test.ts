import { expect, test } from 'bun:test';
import { RelationTable } from './_mpqs_hash.js';
import { mpqs, mpqsInternals as P } from './mpqs.js';

function solve(y: bigint) {
  const h = P.newHandle();
  h.N = 143n;
  h.size_of_FB = 6;
  h.debug = true;
  P.mpqs_FB_ctor(h).p.set([0, 0, 2, 3, 5, 7, 11, 13]);
  const relations = new RelationTable(20);
  relations.add({ Y: y, relp: [] });
  const warnings: unknown[] = [],
    saved = console.warn;
  console.warn = (message) => {
    warnings.push(message);
  };
  try {
    return { result: P.mpqs_solve_linear_system(h, relations), warnings };
  } finally {
    console.warn = saved;
  }
}

test('MPQS reports the native nonfatal post-Gauss warning and retains its factors', () => {
  expect(solve(10n)).toEqual({
    result: [
      [13n, 1n],
      [11n, 1n],
    ],
    warnings: ['MPQS: wrong relation found after Gauss'],
  });
});

test('MPQS consistent relation control emits no diagnostic', () => {
  expect(solve(1n)).toEqual({ result: null, warnings: [] });
});

test('MPQS returns a factor discovered after multiplier selection', () => {
  expect(mpqs(9999999999970063n)).toEqual([
    [541n, 1n],
    [18484288354843n, 1n],
  ]);
});
