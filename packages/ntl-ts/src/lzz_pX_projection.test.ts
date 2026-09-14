import { nativeFixtures as loadLiveNative } from "../../../tests/property/native-live.mjs";


import { test, expect } from 'bun:test';
const fixtures = (await loadLiveNative(import.meta.url, "./lzz_pX_projection.native.json.gz")) as { args: (string)[]; error: null | string; errorType: null | string; function: string; result: null | string; seed: number }[];
import boundary from './lzz_pX_projection.baseline.json' with { type: 'json' };
import { zz_pXModulus, zz_pXMultiplier, build, MulMod } from './lzz_pX.js';
import { UpdateMap, ProjectPowers, zz_pXNewArgument, build as buildArgument } from './lzz_pX1.js';
export function ntl_word_projection(
  op: bigint,
  p: bigint,
  param: bigint[],
  packed: bigint[]
): string {
  const w: bigint[][] = [];
  for (let i = 0, k = 1; i < Number(packed[0]); i++) {
    const n = Number(packed[k++]!);
    w.push(packed.slice(k, k + n));
    k += n;
  }
  const F = new zz_pXModulus(param[2] ? w[0]! : null, p),
    B = new zz_pXMultiplier(),
    H = new zz_pXNewArgument();
  const m = Number(param[0]),
    count = Number(param[1]);
  let x: unknown;
  if (op <= 2n || op === 5n) {
    if (param[3]) build(B, w[1]!, F);
  }
  const state = () => [B.b, B.UseFFT, B.val()];
  if (op === 0n) x = state();
  else if (op === 1n) x = MulMod(w[2]!, B, F);
  else if (op === 2n) x = UpdateMap(w[2]!, B, F);
  else if (op === 3n) x = ProjectPowers(w[2]!, count, w[3]!, F);
  else if (op === 4n) {
    if (param[3]) buildArgument(H, w[3]!, F, m);
    x = ProjectPowers(w[2]!, count, H, F);
  } else if (op === 5n) {
    let error: string | null = null;
    try {
      build(B, w[4]!, F);
    } catch (e) {
      error = (e as Error).message;
    }
    const result = (run: () => bigint[]) => {
      try {
        return [null, run()];
      } catch (e) {
        return [(e as Error).message, null];
      }
    };
    x = [error, state(), result(() => MulMod(w[2]!, B, F)), result(() => UpdateMap(w[2]!, B, F))];
  } else throw new Error('unknown NTL projection operation');
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
  test('native NTL word projection ' + i, () => {
    let result: string | null = null,
      error: string | null = null,
      errorType: string | null = null;
    try {
      result = (ntl_word_projection as Function)(...row.args.map(arg));
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
test('NTL cached multiplier owns its input and returns independent val arrays', () => {
  const f = [1n, 0n, 1n],
    b = [1n, 1n],
    F = new zz_pXModulus(f, 5n),
    B = new zz_pXMultiplier(b, F),
    v = B.val();
  b[0] = 99n;
  v[1] = 99n;
  expect(B.b).toEqual([1n, 1n]);
  expect(MulMod([1n, 2n], B, F)).toEqual([4n, 3n]);
  expect(f).toEqual([1n, 0n, 1n]);
});
test('NTL projection outputs preserve input vectors and cache state', () => {
  const F = new zz_pXModulus([1n, 0n, 1n], 5n),
    B = new zz_pXMultiplier([1n, 1n], F),
    a = [1n, 2n],
    h = [0n, 1n];
  const x = UpdateMap(a, B, F),
    y = ProjectPowers(a, 5, h, F);
  x[0] = 99n;
  y[0] = 99n;
  expect(a).toEqual([1n, 2n]);
  expect(h).toEqual([0n, 1n]);
  expect(B.val()).toEqual([1n, 1n]);
});
test('native prepared-zero crash boundary has an explicit guard with native validation precedence', () => {
  expect(boundary.records.map((r) => r.returncode)).toEqual([0, -11]);
  const F = new zz_pXModulus([1n, 0n, 1n], 5n),
    H = new zz_pXNewArgument();
  buildArgument(H, [0n, 1n], F, 2);
  expect(ProjectPowers([1n, 2n], 0, [0n, 1n], F)).toEqual([]);
  expect(() => ProjectPowers([1n, 2n], 0, H, F)).toThrow(
    'ProjectPowers: prepared argument requires a positive count'
  );
  expect(() => ProjectPowers([1n, 2n, 0n], 0, H, F)).toThrow('ProjectPowers: bad args');
  expect(() => ProjectPowers([1n, 2n], 0, new zz_pXNewArgument(), F)).toThrow(
    'CompMod: uninitialized argument'
  );
});
