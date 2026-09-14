import { nativeFixtures as loadLiveNative } from "../../../tests/property/native-live.mjs";


import { test, expect } from 'bun:test';
const fixtures = (await loadLiveNative(import.meta.url, "./lzz_pX_word_factor_recovery.native.json.gz")) as { args: (string)[]; error: null | string; errorType: null | string; function: string; result: null | string; seed: number }[];
import { RandomStream } from './ZZ.js';
import { zz_pXModulus } from './lzz_pX.js';
import {
  FindRoots,
  FindRoot,
  RootEDF,
  FindFactors,
  EDFSplit,
  EDF,
  SFCanZass2,
  SFCanZass,
} from './lzz_pXFactoring.js';
export function ntl_word_factor_recovery(
  op: bigint,
  p: bigint,
  d: bigint,
  f: bigint[],
  b: bigint[],
  roots: bigint[],
  packed: bigint[],
  key: bigint[]
): string {
  const k = new zz_pXModulus(null, p).arithmetic,
    F = k.norm(f),
    B = k.norm(b),
    R = roots.map(k.mod),
    groups: [bigint[], number][] = [];
  for (let i = 0; i < packed.length; ) {
    const n = Number(packed[i++]),
      degree = Number(packed[i++]);
    groups.push([k.norm(packed.slice(i, i + n)), degree]);
    i += n;
  }
  const stream = new RandomStream(Uint8Array.from(key, Number));
  let errorType: string | null = null,
    error: string | null = null,
    result: unknown = null;
  try {
    if (op === 0n) result = FindRoots(F, p, stream);
    else if (op === 1n) result = FindRoot(F, p, stream);
    else if (op === 2n) result = RootEDF(F, p, stream);
    else if (op === 3n) result = FindFactors(F, B, R, p);
    else if (op === 4n) result = EDFSplit(F, B, Number(d), p, stream);
    else if (op === 5n) result = EDF(F, B, Number(d), p, stream);
    else if (op === 6n) result = SFCanZass2(groups, B, p, stream);
    else if (op === 7n) result = SFCanZass(F, p, stream);
    else if (op === 8n) result = k.divrem(F, B);
    else throw new Error('unknown factor recovery operation');
  } catch (e) {
    errorType = (e as Error).name;
    error = (e as Error).message;
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
  test('native NTL word factor recovery ' + i, () => {
    let result: string | null = null,
      error: string | null = null,
      errorType: string | null = null;
    try {
      result = ntl_word_factor_recovery(
        BigInt(row.args[0]!),
        BigInt(row.args[1]!),
        BigInt(row.args[2]!),
        parse(row.args[3]!),
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
test('NTL word factor outputs preserve polynomial inputs and own factor storage', () => {
  const f = [2n, 1n, 0n, 1n, 1n],
    b = [0n, 0n, 0n, 1n],
    stream = new RandomStream(new Uint8Array(32));
  const out = EDF(f, b, 2, 3n, stream);
  out[0]!.fill(0n);
  expect(f).toEqual([2n, 1n, 0n, 1n, 1n]);
  expect(b).toEqual([0n, 0n, 0n, 1n]);
  expect(out[1]).toEqual([2n, 1n, 1n]);
});
test('NTL factor recovery preserves ordered roots and group inputs', () => {
  const roots = [0n, 1n],
    f = [0n, 1n, 1n],
    g = [0n, 1n];
  expect(FindFactors(f, g, roots, 2n)).toEqual([
    [0n, 1n],
    [1n, 1n],
  ]);
  expect(roots).toEqual([0n, 1n]);
  const groups: [bigint[], number][] = [[f, 1]];
  const out = SFCanZass2(groups, g, 2n, new RandomStream(new Uint8Array(32)));
  out[0]!.fill(0n);
  expect(groups).toEqual([[[0n, 1n, 1n], 1]]);
  expect(g).toEqual([0n, 1n]);
});
