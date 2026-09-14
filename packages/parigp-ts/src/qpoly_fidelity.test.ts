import { expect, test } from 'bun:test';
import { QPoly_normalize, QPoly_to_fractions, QPoly_to_FpX, permtopol } from './galconj.js';
import { PariError } from './errors.js';

// Permanent native counterparts: ff_pari_qpoly (including original static permtopol).
test('rational polynomials normalize denominator signs, zeros and content', () => {
  const q = { num: [2n, 4n, 0n], den: -6n };
  expect(QPoly_normalize(q)).toEqual({ num: [-1n, -2n], den: 3n });
  expect(QPoly_to_fractions(q)).toEqual([
    [-1n, 3n],
    [-2n, 3n],
  ]);
  expect(q).toEqual({ num: [2n, 4n, 0n], den: -6n });
  expect(QPoly_normalize({ num: [], den: -2n })).toEqual({ num: [], den: 1n });
  expect(QPoly_to_fractions({ num: [0n, 0n], den: 2n })).toEqual([]);
});

test('prime conversion cancels each rational coefficient before inversion', () => {
  expect(QPoly_to_FpX({ num: [2n, 4n, 6n], den: 2n }, 2n)).toEqual([1n, 0n, 1n]);
  expect(QPoly_to_FpX({ num: [17n], den: 34n }, 17n)).toEqual([9n]);
  expect(QPoly_to_FpX({ num: [0n], den: 17n }, 17n)).toEqual([]);
  expect(() => QPoly_to_FpX({ num: [1n], den: 17n }, 17n)).toThrow(
    new PariError('impossible inverse in Fl_inv: Mod(0, 17).')
  );
});

test('polynomial denominator zero is rejected before coefficient conversion', () => {
  expect(() => QPoly_normalize({ num: [], den: 0n })).toThrow(
    new PariError('impossible inverse in gdiv: 0.')
  );
  expect(() => QPoly_to_fractions({ num: [0n], den: 0n })).toThrow(
    new PariError('impossible inverse in gdiv: 0.')
  );
  expect(() => QPoly_to_FpX({ num: [1n], den: 0n }, 0n)).toThrow(
    new PariError('impossible inverse in gdiv: 0.')
  );
});

test('Galois centering preserves negative half ties and magnitude thresholds', () => {
  const M = [
    [0n, 0n],
    [0n, 1n],
  ];
  expect(permtopol([0, 1], [0n, -3n], M, 1n, 6n, 3n)).toEqual({ num: [-3n], den: 1n });
  expect(permtopol([0, 1], [0n, -3n], M, 1n, -6n, 3n)).toEqual({ num: [-3n], den: 1n });
  expect(permtopol([0, 1], [0n, 1n], M, 1n, 17n, -8n)).toEqual({ num: [1n], den: 1n });
});

test('permutation conversion preserves empty and per-coefficient error phases', () => {
  const M = [
    [0n, 0n],
    [0n, 1n],
  ];
  expect(permtopol([0], [0n], [[0n]], 0n, 0n, 0n)).toEqual({ num: [], den: 1n });
  expect(() => permtopol([0, 1], [0n, 0n], M, 1n, 0n, 0n)).toThrow(
    new PariError('impossible inverse in dvmdii: 0.')
  );
  expect(() => permtopol([0, 1], [0n, 1n], M, 0n, 17n, 8n)).toThrow(
    new PariError('impossible inverse in gdiv: 0.')
  );
  expect(() => permtopol([0, 1], [0n, 2n], M, 0n, 17n, 8n)).toThrow(
    new PariError('impossible inverse in dvmdii: 0.')
  );
  expect(() => permtopol([0, 1], [0n], [[0n]], 1n, 17n, 8n)).toThrow(
    new PariError('incorrect type in permtopol [permutation] (t_VECSMALL).')
  );
});
