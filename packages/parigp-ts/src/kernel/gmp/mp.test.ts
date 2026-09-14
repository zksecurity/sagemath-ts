import { expect, test } from 'bun:test';
import { divir, divri, divru } from '../../qfb.js';
import { divrr } from './mp.js';

const low = 1n << 127n;

test('native short real division uses its original remainder test', () => {
  expect(divrr(
    { s: 1, e: -65, p: 64, m: (1n << 64n) - 2n },
    { s: 1, e: 65, p: 128, m: low + (low - 1n) / 3n },
  )).toEqual({ s: 1, e: -130, p: 64, m: 13835058055282163710n });
});

test('integer over real preserves the conversion and quotient rounding order', () => {
  expect(divir(-3n, { s: -1, e: -65, p: 128, m: low + (low - 1n) / 3n }))
    .toEqual({ s: 1, e: 66, p: 128, m: 191408831393027885698148216680369618944n });
});

test('real over integer selects the GMP path at the signed-word boundary', () => {
  expect(divri({ s: -1, e: 65, p: 128, m: low + (low - 1n) / 3n * 2n }, -((1n << 63n) + 1n)))
    .toEqual({ s: 1, e: 2, p: 128, m: 283568639100782052855400932736957590188n });
});

test('unsigned-word division accepts the full word without Number coercion', () => {
  expect(divru({ s: 1, e: 0, p: 64, m: 1n << 63n }, (1n << 64n) - 1n))
    .toEqual({ s: 1, e: -64, p: 64, m: (1n << 63n) + 1n });
});
