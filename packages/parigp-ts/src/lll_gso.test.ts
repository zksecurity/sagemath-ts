import { test, expect } from 'bun:test';
import { legacyNativeFixtures } from '../../../tests/property/native-live.mjs';
const fixtures = await legacyNativeFixtures(import.meta.url, './lll_gso.fixtures.json');
import {
  drop,
  potential,
  spread,
  condition_bound,
  GS_extraprec,
  gramschmidt_upper,
  gramschmidt_dynprec,
  RgM_Cholesky_dynprec,
} from './lll.js';
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
  test(`native adaptive Gram-Schmidt ${fixture.probe}`, () => {
    const [op, mode, m, n, p, shift] = fixture.args.slice(0, 6) as number[];
    const flat = fixture.args[6] as string[];
    const A = Array.from({ length: n! }, (_, j) =>
      Array.from({ length: m! }, (_, i) => BigInt(flat[i * n! + j]!))
    );
    const saved = structuredClone(A);
    let z: unknown;
    if (op === 0) z = gramschmidt_upper(A);
    else if (op === 1) z = gramschmidt_dynprec(A);
    else if (op === 2)
      z = RgM_Cholesky_dynprec(
        A.map((a) => A.map((b) => a.reduce((s, v, i) => s + v * b[i]!, 0n)))
      );
    else {
      const M = mode ? A.map((c) => c.map((v) => shiftr(itor(v, p!), shift!))) : A;
      z = [
        drop(M),
        potential(M),
        spread(M),
        condition_bound(M),
        condition_bound(M, true),
        GS_extraprec(M),
        GS_extraprec(M, true),
      ];
    }
    expect(encode(z)).toEqual(fixture.expected);
    expect(A).toEqual(saved);
  });
test('adaptive Gram-Schmidt shape, type and unrepresentable precision boundaries', () => {
  for (const f of [
    drop,
    potential,
    spread,
    condition_bound,
    GS_extraprec,
    gramschmidt_upper,
    gramschmidt_dynprec,
    RgM_Cholesky_dynprec,
  ])
    expect(() => f([])).toThrow(RangeError);
  expect(() => gramschmidt_dynprec([[1n], [2n]])).toThrow(RangeError);
  expect(() => gramschmidt_dynprec([[1n, 2n], [3n]])).toThrow(RangeError);
  expect(() => RgM_Cholesky_dynprec([[1n, 0n]])).toThrow(RangeError);
  expect(() =>
    drop([
      [1n, 0n],
      [0n, itor(1n, 64)],
    ])
  ).toThrow('drop requires matching integer or real diagonal types');
  // Singular upper input violates the full-rank precondition: its exponent
  // sentinel would request an impossible native allocation.
  expect(() =>
    gramschmidt_upper([
      [0n, 0n],
      [0n, 1n],
    ])
  ).toThrow('LLL precision exceeds exact JavaScript integer range');
});
