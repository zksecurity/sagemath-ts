import { nativeFixtures as loadLiveNative } from "../../../tests/property/native-live.mjs";


import { expect, test } from 'bun:test';
const fixtures = (await loadLiveNative(import.meta.url, "./gen2_ffelt.native.json.gz")) as { args: (string)[]; error: null; errorType: null; function: string; result: string; seed: number }[];
import { cmp_universal } from './gen2.js';
import { PariType } from './types.js';

for (const row of fixtures) {
  test(`bundled PARI finite-field comparison ${row.seed}`, () => {
    const parse = (s: string) =>
      s.slice(1, -1).trim()
        ? s
            .slice(1, -1)
            .split(',')
            .map((c) => BigInt(c.trim()))
        : [];
    const p = BigInt(row.args[0]!),
      T = parse(row.args[1]!),
      x = parse(row.args[2]!);
    const q = BigInt(row.args[3]!),
      U = parse(row.args[4]!),
      y = parse(row.args[5]!);
    const original = [[...T], [...x], [...U], [...y]];
    const scalar = row.function === 'pari_ffelt_compare_scalar';
    const result = cmp_universal(
      {
        type: PariType.t_FFELT,
        p,
        degree: T.length - 1,
        definingPoly: scalar ? undefined : T,
        value: scalar ? x[0]! : x,
      },
      {
        type: PariType.t_FFELT,
        p: q,
        degree: U.length - 1,
        definingPoly: scalar ? undefined : U,
        value: scalar ? y[0]! : y,
      }
    );
    expect({ result: String(result), error: null, errorType: null }).toEqual({
      result: row.result,
      error: row.error,
      errorType: row.errorType,
    });
    expect([T, x, U, y]).toEqual(original);
  });
}
