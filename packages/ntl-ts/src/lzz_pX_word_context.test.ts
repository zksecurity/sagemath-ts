import { nativeFixtures as loadLiveNative } from "../../../tests/property/native-live.mjs";


import { test, expect } from 'bun:test';
const fixtures = (await loadLiveNative(import.meta.url, "./lzz_pX_word_context.native.json.gz")) as { args: (string)[]; error: null | string; errorType: null | string; function: string; result: null | string; seed: number }[];
import { RandomStream } from './ZZ.js';
import { zz_pXModulus, build } from './lzz_pX.js';
import { GCD, MinPolySeq, BuildFromRoots, MinPolyMod } from './lzz_pX1.js';
import { SFCanZass, EDF, FindFactors, SFCanZass1, NewDDF } from './lzz_pXFactoring.js';
export function ntl_word_context(
  op: bigint,
  p: bigint,
  maxroot: bigint,
  d: bigint,
  f: bigint[],
  g: bigint[],
  values: bigint[],
  key: bigint[]
): string {
  const options = { maxroot: Number(maxroot) },
    F = new zz_pXModulus(null, p, options),
    k = F.arithmetic,
    ff = k.norm(f),
    gg = k.norm(g),
    v = values.map(k.mod);
  const pc = F.PrimeCnt,
    profile = [
      pc,
      pc,
      F.MaxRoot,
      F.modCrossover,
      pc === 1 ? 150 : pc === 2 ? 300 : 500,
      pc === 1 ? 90 : pc === 2 ? 180 : 350,
      pc === 1 ? 400 : pc === 2 ? 800 : 1400,
      pc === 1 ? 480 : pc === 2 ? 900 : 1600,
    ];
  const stream = new RandomStream(Uint8Array.from(key, Number));
  let errorType: string | null = null,
    error: string | null = null,
    result: unknown = null;
  try {
    if (op === 0n) {
      build(F, ff);
      result = [F.n, F.n > F.modCrossover + 1 ? 1 : 0];
    } else if (op === 1n) result = GCD(ff, gg, p, options);
    else if (op === 2n) result = MinPolySeq(v, Number(d), p, options);
    else if (op === 3n) result = BuildFromRoots(v, p, options);
    else if (op === 4n) result = SFCanZass(ff, p, stream, options);
    else if (op === 5n) result = EDF(ff, gg, Number(d), p, stream, options);
    else if (op === 6n) {
      build(F, ff);
      result = MinPolyMod(gg, F, Number(d), stream);
    } else if (op === 7n) result = FindFactors(ff, gg, v, p, options);
    else if (op === 8n) result = SFCanZass1(ff, p, options);
    else if (op === 9n) result = [NewDDF(ff, gg, p, options), gg];
    else throw new Error('unknown NTL word context operation');
  } catch (e) {
    errorType = (e as Error).name;
    error = (e as Error).message;
  }
  const tail = Array.from(stream.get(64), (x) => x.toString(16).padStart(2, '0')).join('');
  return JSON.stringify([profile, errorType, error, result, tail], (_, v) =>
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
  test('native NTL word context ' + i, () => {
    let result: string | null = null,
      error: string | null = null,
      errorType: string | null = null;
    try {
      result = ntl_word_context(
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
test('NTL polynomial contexts capture constructor settings by value', () => {
  const settings = { maxroot: 6 },
    F = new zz_pXModulus(null, 101n, settings);
  settings.maxroot = 8;
  expect(F.options.maxroot).toBe(6);
  expect(Object.isFrozen(F.options)).toBe(true);
  expect(() => build(F, [1n, ...Array<bigint>(64).fill(0n), 1n])).toThrow(
    'Polynomial too big for FFT'
  );
});
test('NTL explicit context factors preserve input arrays and options', () => {
  const f = [2n, 1n, 0n, 1n, 1n],
    options = { maxroot: 60 },
    stream = new RandomStream(new Uint8Array(32));
  const out = SFCanZass(f, 3n, stream, options);
  out[0]!.fill(0n);
  expect(f).toEqual([2n, 1n, 0n, 1n, 1n]);
  expect(options).toEqual({ maxroot: 60 });
  expect(out[1]).toEqual([2n, 1n, 1n]);
});
