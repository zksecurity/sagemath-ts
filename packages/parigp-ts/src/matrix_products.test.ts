import { test, expect } from 'bun:test';
import { legacyNativeFixtures } from '../../../tests/property/native-live.mjs';
const fixtures = await legacyNativeFixtures(import.meta.url, './matrix_products.fixtures.json');
import { RgM_rescale_to_int } from './polarit2.js';
import { gram_matrix, RgM_mul, RgV_dotsquare, RgV_dotproduct, type RgScalar } from './RgV.js';
import { itor, shiftr } from './qfb.js';
import { PariError } from './errors.js';
function encode(x: unknown): unknown {
  if (typeof x === 'bigint') return String(x);
  if (Array.isArray(x)) return x.map(encode);
  if (x !== null && typeof x === 'object') {
    const r = x as { s: number; e: number; m: bigint; p: number };
    return [r.s, String(r.e), String(r.m), r.p];
  }
  return x;
}
for (const fixture of fixtures)
  test(`native matrix product and rescaling ${fixture.probe}`, () => {
    const [op, mode, mode2, m, k, n, p, q, shift] = fixture.args.slice(0, 9) as number[];
    const matrix = (
      m: number,
      k: number,
      mode: number,
      p: number,
      shift: number,
      flat: string[]
    ): RgScalar[][] =>
      Array.from({ length: k }, (_, j) =>
        Array.from({ length: m }, (_, i) => {
          const x = BigInt(flat[i * k + j]!);
          return mode === 1 || (mode === 2 && (i + j) % 2) || (mode === 3 && x !== 0n)
            ? shiftr(itor(x, p), shift)
            : x;
        })
      );
    const A = matrix(m!, k!, mode!, p!, shift!, fixture.args[9] as string[]),
      B = matrix(k!, n!, mode2!, q!, -shift!, fixture.args[10] as string[]);
    const saved = structuredClone([A, B]);
    let value: unknown, error: unknown;
    try {
      value =
        op === 0
          ? RgM_rescale_to_int(A)
          : op === 1
            ? gram_matrix(A)
            : op === 2
              ? RgM_mul(A, B)
              : op === 3
                ? RgV_dotsquare(A[0]!)
                : op === 4
                  ? RgV_dotproduct(A[0]!, A[1]!)
                  : RgV_dotproduct(A[0]!, A[0]!);
    } catch (e) {
      error = e;
    }
    if (typeof fixture.error === 'string') {
      expect(error).toBeInstanceOf(PariError);
      expect((error as Error).message).toBe(fixture.error);
    } else {
      expect(error).toBeUndefined();
      expect(encode(value)).toEqual(fixture.expected);
    }
    expect([A, B]).toEqual(saved);
  });
test('integer/real matrix shape guards and empty products', () => {
  expect(() => RgM_rescale_to_int([[1n], [2n, 3n]])).toThrow(RangeError);
  expect(() => gram_matrix([[1n], [2n, 3n]])).toThrow(RangeError);
  expect(() => RgV_dotproduct([1n], [])).toThrow(RangeError);
  expect(() => RgM_mul([[1n]], [[1n, 2n]])).toThrow(RangeError);
  expect(() => RgM_mul([[1n], [2n, 3n]], [[1n, 2n]])).toThrow(RangeError);
  expect(RgM_mul([], [[], []])).toEqual([[], []]);
  expect(RgM_mul([[1n]], [])).toEqual([]);
  expect(RgV_dotsquare([])).toBe(0n);
  expect(RgV_dotproduct([], [])).toBe(0n);
});
