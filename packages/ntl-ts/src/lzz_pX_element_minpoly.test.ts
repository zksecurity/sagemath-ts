import { nativeFixtures as loadLiveNative } from "../../../tests/property/native-live.mjs";


import { test, expect } from 'bun:test';
const fixtures = (await loadLiveNative(import.meta.url, "./lzz_pX_element_minpoly.native.json.gz")) as { args: (string)[]; error: null | string; errorType: null | string; function: string; result: null | string; seed: number }[];
import { zz_pXModulus } from './lzz_pX.js';
import { DoMinPolyMod, IrredPolyMod, ProbMinPolyMod, MinPolyMod } from './lzz_pX1.js';
import { RandomStream } from './ZZ.js';
export function ntl_element_minpoly(
  op: bigint,
  p: bigint,
  m: bigint,
  ready: bigint,
  f: bigint[],
  g: bigint[],
  r: bigint[],
  key: bigint[]
): string {
  const F = new zz_pXModulus(ready ? f : null, p),
    G = F.arithmetic.norm(g),
    R = r.map((x) => F.arithmetic.mod(x)),
    stream = new RandomStream(Uint8Array.from(key, Number));
  let errorType: string | null = null,
    error: string | null = null,
    result: bigint[] | null = null;
  try {
    if (op === 0n) result = DoMinPolyMod(G, F, Number(m), R);
    else if (op === 1n) result = IrredPolyMod(G, F, Number(m));
    else if (op === 2n) result = ProbMinPolyMod(G, F, Number(m), stream);
    else if (op === 3n) result = MinPolyMod(G, F, Number(m), stream);
    else if (op === 4n) result = IrredPolyMod(G, F);
    else if (op === 5n) result = ProbMinPolyMod(G, F, stream);
    else if (op === 6n) result = MinPolyMod(G, F, stream);
    else throw new Error('unknown element minimum polynomial operation');
  } catch (e) {
    error = (e as Error).message;
    errorType = (e as Error).name;
  }
  const tail = Array.from(stream.get(64), (x) => x.toString(16).padStart(2, '0')).join('');
  return JSON.stringify([errorType, error, result, tail], (_, v) =>
    typeof v === 'bigint' ? String(v) : v
  );
}

const parse = (s: string): bigint[] =>
  s
    .slice(1, -1)
    .split(',')
    .filter((x) => x.trim())
    .map((x) => BigInt(x.trim()));
for (const [i, row] of fixtures.entries())
  test('native NTL element minimum polynomial ' + i, () => {
    let result: string | null = null,
      error: string | null = null,
      errorType: string | null = null;
    try {
      result = ntl_element_minpoly(
        BigInt(row.args[0]!),
        BigInt(row.args[1]!),
        BigInt(row.args[2]!),
        BigInt(row.args[3]!),
        parse(row.args[4]!),
        parse(row.args[5]!),
        parse(row.args[6]!),
        parse(row.args[7]!)
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
test('NTL element minimum polynomials own outputs and preserve polynomial inputs', () => {
  const f = [1n, 1n, 1n],
    g = [0n, 1n],
    R = [1n],
    F = new zz_pXModulus(f, 2n);
  const h = DoMinPolyMod(g, F, 2, R);
  h.fill(0n);
  expect(IrredPolyMod(g, F)).toEqual([1n, 1n, 1n]);
  expect(f).toEqual([1n, 1n, 1n]);
  expect(g).toEqual([0n, 1n]);
  expect(R).toEqual([1n]);
  expect(F.f).toEqual(f);
});
test('NTL certified element minimum polynomial preserves explicit stream identity', () => {
  const F = new zz_pXModulus([0n, 0n, 1n, 0n, 1n], 2n),
    stream = new RandomStream(new Uint8Array(32)),
    alias = stream;
  expect(MinPolyMod([0n, 0n, 1n], F, stream)).toEqual([0n, 1n, 1n]);
  expect(alias).toBe(stream);
  expect(Array.from(alias.get(4))).toEqual([64, 93, 106, 229]);
});
