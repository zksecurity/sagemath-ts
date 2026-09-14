import { test, expect } from 'bun:test';
import { legacyNativeFixtures } from '../../../tests/property/native-live.mjs';
const fixtures = await legacyNativeFixtures(import.meta.url, './real_matrix.fixtures.json');
import { qfgaussred_positive, RgM_Cholesky } from './alglin2.js';
import { RgM_inv_upper } from './alglin1.js';
import {
  type MatrixReal,
  matrixRealMul,
  matrixRealSub,
  matrixRealNeg,
  matrixRealDiv,
  matrixRealInv,
} from './_real_matrix.js';
import { itor, shiftr, real_0_bit } from './qfb.js';
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
  test(`native real matrix and numeric branch ${fixture.probe}`, () => {
    const [op, mode, n, p, q, shift] = fixture.args.slice(0, 6) as number[];
    const flat = fixture.args[6] as string[];
    const A: MatrixReal[][] = Array.from({ length: n! }, (_, j) =>
      Array.from({ length: n! }, (_, i) => {
        const x = BigInt(flat[i * n! + j]!),
          bits = mode! & 2 && j === 0 ? p! : q!;
        return mode! & 1 && x === 0n
          ? 0n
          : mode! & 4 && x === 0n
            ? real_0_bit(shift! - bits)
            : shiftr(itor(x, bits), shift!);
      })
    );
    const saved = structuredClone(A);
    let value: unknown, error: unknown;
    try {
      const x = A[0]?.[0]!,
        y = A[1]?.[0]!;
      value =
        op === 0
          ? qfgaussred_positive(A)
          : op === 1
            ? RgM_Cholesky(A, p)
            : op === 2
              ? RgM_inv_upper(A)
              : op === 3
                ? matrixRealMul(x, y)
                : op === 4
                  ? matrixRealSub(x, y)
                  : op === 5
                    ? matrixRealNeg(x)
                    : op === 6
                      ? matrixRealDiv(x, y)
                      : matrixRealInv(x);
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
    expect(A).toEqual(saved);
  });
test('real matrix dimensions, precision and empty-matrix boundaries', () => {
  expect(qfgaussred_positive([])).toEqual([]);
  expect(RgM_Cholesky([])).toEqual([]);
  expect(RgM_inv_upper([])).toEqual([]);
  expect(() => qfgaussred_positive([[0n, 0n]])).toThrow(
    'inconsistent dimensions in qfgaussred_positive'
  );
  expect(() => RgM_inv_upper([[0n, 0n]])).toThrow(RangeError);
  for (const p of [0, -64, 65, NaN, Infinity])
    expect(() => RgM_Cholesky([], p)).toThrow(RangeError);
});
test('large allocated real-zero inversion keeps the documented deterministic rejection', () => {
  for (const e of [-9192, -8192, -7192]) {
    const z = { ...itor(0n, 8192), e };
    expect(() => matrixRealInv(z)).toThrow('impossible inverse in invr: 0');
    expect(() => RgM_inv_upper([[z]])).toThrow('impossible inverse in invr: 0');
  }
});
