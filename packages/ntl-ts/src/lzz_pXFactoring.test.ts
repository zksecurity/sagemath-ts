import { nativeFixtures as loadLiveNative } from "../../../tests/property/native-live.mjs";


import { test, expect } from 'bun:test';
const fixtures = (await loadLiveNative(import.meta.url, "./lzz_pXFactoring.native.json.gz")) as { args: (string)[]; error: null | string; errorType: null | string; function: string; result: null | string; seed: number }[];
import { zz_pXModulus, PowerXMod } from './lzz_pX.js';
import { NewDDF, SFCanZass1 } from './lzz_pXFactoring.js';
export function ntl_word_ddf(op: bigint, p: bigint, f: bigint[], h: bigint[]): string {
  let x: unknown;
  if (op === 0n) x = NewDDF(f, h, p);
  else if (op === 1n) x = SFCanZass1(f, p);
  else if (op === 2n) x = NewDDF(f, PowerXMod(p, new zz_pXModulus(f, p)), p);
  else throw new Error('unknown NTL distinct-degree operation');
  return JSON.stringify(x, (_, v) =>
    typeof v === 'bigint' || typeof v === 'number' ? String(v) : v
  );
}

const arg = (s: string): bigint | bigint[] =>
  s.startsWith('[')
    ? s.slice(1, -1).trim()
      ? s
          .slice(1, -1)
          .split(',')
          .map((v) => BigInt(v.trim()))
      : []
    : BigInt(s);
for (const [i, row] of fixtures.entries())
  test('native NTL word distinct-degree ' + i, () => {
    let result: string | null = null,
      error: string | null = null,
      errorType: string | null = null;
    try {
      result = (ntl_word_ddf as Function)(...row.args.map(arg));
    } catch (e) {
      error = (e as Error).message;
      errorType = (e as Error).name;
    }
    expect({ result, error, errorType }).toEqual({
      result: row.result,
      error: row.error,
      errorType: row.errorType,
    });
  }, 10000);
test('NTL distinct-degree factors preserve inputs and return independent factors', () => {
  const f = [-2n, 2n, -1n, 1n],
    h = [3n, 4n, 4n],
    factors = NewDDF(f, h, 5n);
  expect(factors).toEqual([
    [[4n, 1n], 1],
    [[2n, 0n, 1n], 2],
  ]);
  factors[0]![0][0] = 99n;
  expect(f).toEqual([-2n, 2n, -1n, 1n]);
  expect(h).toEqual([3n, 4n, 4n]);
  expect(factors[1]).toEqual([[2n, 0n, 1n], 2]);
});
test('NTL squarefree first stage returns independent Frobenius data', () => {
  const f = [-2n, 2n, -1n, 1n],
    [factors, h] = SFCanZass1(f, 5n);
  h[0] = 99n;
  expect(factors).toEqual([
    [[4n, 1n], 1],
    [[2n, 0n, 1n], 2],
  ]);
  expect(f).toEqual([-2n, 2n, -1n, 1n]);
});
