import { test, expect } from 'bun:test';
import { legacyNativeFixtures } from '../../../tests/property/native-live.mjs';
const fixtures = await legacyNativeFixtures(import.meta.url, './modular_linear.fixtures.json');
import { Flm_pivots, Flm_gauss } from './Flv.js';
import { ZM_pivots, ZM_rank, ZM_gauss } from './alglin1.js';
import {
  wordMatrixPivots,
  wordMatrixSolve,
  integerMatrixPivots,
  integerMatrixSolve,
} from './_matrix_inverse.js';
import { PariError } from './errors.js';
type Fixture = {
  kind: string;
  probe: number;
  args: (number | string | string[])[];
  expected?: unknown;
  error?: string;
};
for (const fixture of fixtures as Fixture[])
  test(`native modular ${fixture.kind} linear probe ${fixture.probe}`, () => {
    const [op, m, n, k] = fixture.args.slice(0, 4) as number[],
      word = fixture.kind === 'word';
    const p = word ? BigInt(fixture.args[4] as string) : 0n;
    const matrix = (m: number, n: number, flat: string[]) => [
      [],
      ...Array.from({ length: n }, (_, j) => [
        0n,
        ...Array.from({ length: m }, (_, i) => BigInt(flat[i * n + j]!)),
      ]),
    ];
    const A = matrix(m!, n!, fixture.args[word ? 5 : 4] as string[]),
      B = matrix(m!, k!, fixture.args[word ? 6 : 5] as string[]),
      saved = structuredClone([A, B]);
    let actual: unknown,
      error: unknown,
      X: bigint[][] | null = null,
      den = 1n;
    try {
      if (op === 0) {
        const [d, r] = word ? Flm_pivots(A, p) : ZM_pivots(A);
        actual = [d?.slice(1) ?? null, r];
      } else if (word) {
        X = Flm_gauss(A, B, p);
        actual = X?.slice(1).map((c) => c.slice(1)) ?? null;
      } else if (op === 1) actual = ZM_rank(A);
      else {
        const z = ZM_gauss(A, B);
        X = z?.[0] ?? null;
        den = z?.[1] ?? 1n;
        actual = z === null ? null : [z[0].slice(1).map((c) => c.slice(1)), z[1]];
      }
    } catch (e) {
      error = e;
    }
    if (fixture.error !== undefined) {
      expect(error).toBeInstanceOf(PariError);
      expect((error as Error).message).toBe(fixture.error);
    } else {
      expect(error).toBeUndefined();
      expect(
        JSON.parse(JSON.stringify(actual, (_, v) => (typeof v === 'bigint' ? String(v) : v)))
      ).toEqual(fixture.expected);
    }
    expect([A, B]).toEqual(saved);
    if (X && m === n && n) {
      for (let j = 1; j <= k!; j++)
        for (let i = 1; i <= n; i++) {
          let v = 0n;
          for (let c = 1; c <= n; c++) v += A[c]![i]! * X[j]![c]!;
          expect(word ? (((v - B[j]![i]!) % p) + p) % p : v).toBe(word ? 0n : B[j]![i]! * den);
        }
    }
  });
test('modular matrix adapter shape and word-prime guards', () => {
  for (const p of [0n, 1n, 1n << 64n]) {
    expect(() => Flm_pivots([[]], p)).toThrow(RangeError);
    expect(() => Flm_gauss([[]], [[]], p)).toThrow(RangeError);
  }
  for (const A of [
    [[], []],
    [[], [0n, 1n], [0n]],
  ]) {
    expect(() => Flm_pivots(A, 101n)).toThrow(RangeError);
    expect(() => ZM_pivots(A)).toThrow(RangeError);
  }
  expect(() => Flm_gauss([[], [0n, 1n]], [[], [0n]], 101n)).toThrow(RangeError);
  expect(() => ZM_gauss([[], [0n, 1n]], [[], [0n]])).toThrow(RangeError);
  expect(() => ZM_gauss([[], [0n, 1n], [0n, 1n]], [[]])).toThrow(RangeError);
  for (const f of [wordMatrixPivots, integerMatrixPivots])
    expect(() => f([], 1, 1, 101n)).toThrow(RangeError);
  expect(() => wordMatrixSolve([[1n]], [], 1, 1, 1, 101n)).toThrow(RangeError);
  expect(() => integerMatrixSolve([[1n]], [], 1, 1, 1)).toThrow(RangeError);
});
