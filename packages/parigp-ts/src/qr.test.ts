import { expect, test } from 'bun:test';
import { legacyNativeFixtures } from '../../../tests/property/native-live.mjs';
const fixtures = await legacyNativeFixtures(import.meta.url, './qr.fixtures.json');
import { QR_init, R_from_QR, gaussred_from_QR, type QrScalar } from './bibli1.js';
import { itor, shiftr } from './qfb.js';
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
  test(`native Householder QR probe ${fixture.probe}`, () => {
    const [op, mode, m, n, p, q, shift] = fixture.args.slice(0, 7) as number[];
    const flat = fixture.args[7] as string[];
    const A: QrScalar[][] = Array.from({ length: n! }, (_, j) =>
      Array.from({ length: m! }, (_, i) => {
        const x = BigInt(flat[i * n! + j]!);
        return mode === 1 || (mode === 2 && (i + j) % 2 !== 0) ? shiftr(itor(x, q!), shift!) : x;
      })
    );
    const saved = structuredClone(A);
    const actual = op === 0 ? QR_init(A, p) : op === 1 ? R_from_QR(A, p) : gaussred_from_QR(A, p);
    expect(encode(actual)).toEqual(fixture.expected);
    expect(A).toEqual(saved);
    if (op === 0 && actual?.[0] === 1) {
      const [, B, Q, L] = actual as ReturnType<typeof QR_init>;
      expect(B!.length).toBe(n!);
      expect(Q!.length).toBe(n! - 1);
      for (let j = 0; j < n!; j++) for (let i = 0; i < j; i++) expect(L![j]![i]).toBe(0n);
    }
  });
test('Householder typed adapter rejects unsupported shape and precision', () => {
  for (const A of [[], [[]], [[1n], [1n]], [[1n, 2n], [3n]]])
    expect(() => QR_init(A)).toThrow(RangeError);
  for (const p of [0, -64, 65, Infinity, NaN]) expect(() => QR_init([[1n]], p)).toThrow(RangeError);
  expect(() => gaussred_from_QR([[1n, 0n]])).toThrow(RangeError);
});
