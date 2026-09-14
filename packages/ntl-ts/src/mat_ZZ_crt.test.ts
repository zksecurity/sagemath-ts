import { nativeFixtures as loadLiveNative } from "../../../tests/property/native-live.mjs";
import { test, expect } from 'bun:test';


const fixtures = (await loadLiveNative(import.meta.url, "./mat_ZZ_crt.native.json.gz")) as { args: string[]; result: string | null; error: string | null; errorType: string | null }[];
import { CRT as scalarCRT, CRTInRange } from './ZZ.js';
import { CRT, DetBound } from './mat_ZZ.js';
import { determinant } from './mat_ZZ_p.js';
export function ntl_integer_crt(
  op0: bigint,
  a: bigint,
  p: bigint,
  dims: bigint[],
  g: bigint[],
  G: bigint[]
): string {
  const op = Number(op0),
    [n, m, r, c, badG, badR] = dims.map(Number);
  const make = (v: bigint[], rows: number, cols: number, bad: number): bigint[][] =>
    bad
      ? [[0n], [0n, 0n]]
      : Array.from({ length: Math.max(0, rows) }, (_, i) => v.slice(i * cols, (i + 1) * cols));
  let result: unknown;
  if (op <= 2) {
    result = g.map((gg, i) => {
      try {
        return [
          null,
          null,
          op === 0 ? CRTInRange(gg, a) : scalarCRT(gg, a, G[i]!, p, { word: op === 2 }),
        ];
      } catch (e) {
        return [(e as Error).name, (e as Error).message, null];
      }
    });
  } else {
    if (op === 3) {
      if (p <= 1n) throw new Error('zz_pContext: p must be > 1');
      if (p >= 1n << 60n) throw new Error('zz_pContext: modulus too big');
    }
    if (op === 5 && p <= 1n) throw new Error('ZZ_pContext: p must be > 1');
    if (n! < 0 || m! < 0 || (op === 3 && (r! < 0 || c! < 0))) throw new Error('SetDims: bad args');
    const M = make(g, n!, m!, badG!);
    if (op === 3) result = CRT(M, a, make(G, r!, c!, badR!), p, { columns: m, residueColumns: c });
    else if (op === 4) result = DetBound(M, { columns: m });
    else if (op === 5) result = determinant(M, p, { columns: m });
    else throw new Error('unknown integer CRT operation');
  }
  return JSON.stringify(result, (_, v) => (typeof v === 'bigint' ? String(v) : v));
}
const parse = (s: string): bigint[] =>
  s
    .slice(1, -1)
    .split(',')
    .filter((x) => x.trim())
    .map((x) => BigInt(x.trim()));
for (const [i, row] of fixtures.entries())
  test('native NTL integer CRT and modular determinant ' + i, () => {
    let result: string | null = null,
      error: string | null = null,
      errorType: string | null = null;
    try {
      result = ntl_integer_crt(
        BigInt(row.args[0]!),
        BigInt(row.args[1]!),
        BigInt(row.args[2]!),
        parse(row.args[3]!),
        parse(row.args[4]!),
        parse(row.args[5]!)
      );
    } catch (e) {
      error = (e as Error).message;
      errorType = (e as Error).name;
    }
    expect({ result, error, errorType }).toEqual({
      result: row.result,
      error: row.error,
      errorType: row.errorType,
    });
  }, 30000);

test('NTL matrix CRT owns output rows and preserves both input matrices', () => {
  const g = [[1n, -1n]],
    G = [[0n, 1n]],
    before = structuredClone([g, G]);
  const [changed, M, a] = CRT(g, 3n, G, 2n);
  expect([changed, M, a]).toEqual([1, [[-2n, -1n]], 6n]);
  M[0]![0] = 999n;
  expect([g, G]).toEqual(before);
});
test('NTL bounds and modular determinants preserve signed input rows on success and failure', () => {
  const A = [
      [3n, 4n],
      [0n, 1n],
    ],
    before = structuredClone(A);
  expect(DetBound(A)).toBe(3);
  expect(determinant(A, 17n)).toBe(3n);
  expect(() => determinant(A, 9n)).toThrow('InvMod: inverse undefined');
  expect(A).toEqual(before);
});
