import { test, expect } from 'bun:test';
import { legacyNativeFixtures } from '../../../tests/property/native-live.mjs';
const fixtures = await legacyNativeFixtures(import.meta.url, './lll_wrapper.fixtures.json');
import resources from './lll_wrapper.resources.json' with { type: 'json' };
import {
  flat,
  ZM_flatter,
  ZM_flatter_rank,
  ZM_flattergram,
  ZM_lll,
  lllfp,
  fplll_fast,
  fplll_dpe,
  fplll_heuristic,
  fplll,
} from './lll.js';
import { itor, shiftr, mkqfb, redimagsl2, type MpReal } from './qfb.js';
import { gram_matrix } from './RgV.js';
import { integerMatrixPivots } from './_matrix_inverse.js';

type Fixture = {
  group: string;
  probe: number;
  function?: string;
  args: unknown[];
  result?: string;
  error?: string;
  errorType?: string;
};
const matrix = (m: number, n: number, values: unknown): bigint[][] =>
  Array.from({ length: n }, (_, j) =>
    Array.from({ length: m }, (_, i) => BigInt((values as string[])[i * n + j]!))
  );
const encode = (x: unknown): unknown =>
  typeof x === 'bigint'
    ? String(x)
    : Array.isArray(x)
      ? x.map(encode)
      : x && typeof x === 'object'
        ? [(x as MpReal).s, String((x as MpReal).e), String((x as MpReal).m), (x as MpReal).p]
        : x;
function invoke(f: Fixture): unknown {
  const a = f.args,
    numbers = a.slice(0, -1).map(Number);
  if (f.group === 'stage') {
    const fast = f.function === 'pari_lll_fast',
      heuristic = f.function === 'pari_lll_heuristic';
    const [mode, m, n, dn, dd, en, ed, keep, track, want, p] = heuristic
      ? [0, ...numbers]
      : numbers;
    const B = matrix(m!, n!, a.at(-1)),
      saved = structuredClone(B);
    let result: unknown;
    if (fast) result = fplll_fast(B, dn! / dd!, en! / ed!, !!keep, !!track);
    else if (heuristic) result = fplll_heuristic(B, dn! / dd!, en! / ed!, !!keep, !!track, want, p);
    else {
      const G = mode ? (gram_matrix(B) as bigint[][]) : null,
        input = mode === 2 ? null : B;
      result =
        f.function === 'pari_lll_dpe'
          ? fplll_dpe(input, G, dn! / dd!, en! / ed!, !!keep, !!track, !!want)
          : fplll(input, G, dn! / dd!, en! / ed!, !!keep, !!track, !!want, p);
    }
    expect(B).toEqual(saved);
    return result;
  }
  const [op, flag, mode, m, n, p, shift] = numbers;
  const values = (a.at(-1) as string[]).map(BigInt);
  const B = matrix(m!, n!, a.at(-1)).map((c, j) =>
    c.map((v, i) => (mode === 1 || (mode === 2 && (i + j) % 2) ? shiftr(itor(v, p), shift!) : v))
  );
  const saved = structuredClone(B);
  try {
    if (op === 0) return ZM_lll(B as bigint[][], 0.99, flag);
    if (op === 1) return lllfp(B, 0.99, flag);
    if (op === 2) {
      const z = flat(B as bigint[][], flag);
      return [z[0], z[1] ?? 0n, z[2], z[3]];
    }
    if (op === 3) return ZM_flatter(B as bigint[][], flag);
    if (op === 4) {
      const rows = Array.from({ length: m! }, (_, i) => B.map((c) => c[i] as bigint));
      return ZM_flatter_rank(B as bigint[][], n! - integerMatrixPivots(rows, m!, n!)[1], flag);
    }
    if (op === 5) return ZM_flattergram(B as bigint[][], flag);
    const [aa, bb, cc] = values,
      z = redimagsl2(mkqfb(aa!, bb!, cc!, bb! * bb! - 4n * aa! * cc!));
    return [[z.Q.a, z.Q.b, z.Q.c], z.U[0]!.map((_, i) => z.U.map((c) => c[i]!))];
  } finally {
    expect(B).toEqual(saved);
  }
}
for (const fixture of fixtures as Fixture[])
  test(`native adaptive LLL ${fixture.group} ${fixture.probe}`, () => {
    const start = performance.now();
    let actual: unknown;
    try {
      actual = { result: JSON.stringify(encode(invoke(fixture))) };
    } catch (e) {
      actual = { error: (e as Error).message, errorType: (e as Error).name };
    }
    expect(actual).toEqual(
      fixture.error === undefined
        ? { result: fixture.result }
        : { error: fixture.error, errorType: fixture.errorType }
    );
    if (fixture.group === 'extra' && fixture.probe === 162) {
      // Native word-mantissa/shift ordering: the old expanded products took 15.2s,
      // while the corrected port took 2.1s on the audit machine. Keep a wide margin.
      expect(performance.now() - start).toBeLessThan(10000);
    }
  }, 15000);
for (const fixture of resources.observations as Fixture[])
  test(`native Gram FLATTER timeout ${fixture.probe} has an explicit resource boundary`, () => {
    expect(fixture.errorType).toBe('RuntimeError');
    expect(() => invoke(fixture)).toThrow(
      'rank-deficient Gram FLATTER exceeded 64 augmentation attempts'
    );
  }, 15000);
test('adaptive LLL adapter domains', () => {
  expect(() => ZM_lll([[1n], [2n, 3n]])).toThrow(RangeError);
  expect(() => ZM_lll([[1n], [2n]], 0.99, 258)).toThrow('square Gram matrix');
  expect(() => flat([])).toThrow('at least two columns');
  expect(() => ZM_flattergram([[1n]])).toThrow('at least two columns');
  expect(() =>
    lllfp(
      [
        [1n, itor(0n, 64)],
        [itor(0n, 64), itor(1n, 64)],
      ],
      0.99,
      258
    )
  ).toThrow('inexact Gram LLL requires real coefficients or exact zero');
});
