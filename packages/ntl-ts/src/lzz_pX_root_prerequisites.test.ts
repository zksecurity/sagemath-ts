import { nativeFixtures as loadLiveNative } from "../../../tests/property/native-live.mjs";


import { test, expect } from 'bun:test';
const fixtures = (await loadLiveNative(import.meta.url, "./lzz_pX_root_prerequisites.native.json.gz")) as { args: (string)[]; error: null | string; errorType: null | string; function: string; result: null | string; seed: number }[];
import { RandomStream, RandomBndGenerator, VectorRandomBnd } from './ZZ.js';
import { random as randomPolynomial, zz_pXModulus } from './lzz_pX.js';
import { BuildFromRoots } from './lzz_pX1.js';
export function ntl_root_prerequisites(
  op: bigint,
  p: bigint,
  n: bigint,
  data: bigint[],
  key: bigint[]
): string {
  if (op <= 1n) new zz_pXModulus(null, p);
  const stream = new RandomStream(Uint8Array.from(key, Number)),
    hex = (b: Uint8Array) => Array.from(b, (x) => x.toString(16).padStart(2, '0')).join('');
  const secondary = new RandomStream(Uint8Array.from(key, Number));
  secondary.set_nonce(1n);
  let errorType: string | null = null,
    error: string | null = null,
    value: unknown = null;
  try {
    if (op === 0n) value = BuildFromRoots(data, p);
    else if (op === 1n) value = randomPolynomial(Number(n), p, stream);
    else if (op === 2n) value = VectorRandomBnd(Number(n), p, stream);
    else if (op === 3n) {
      const g = [
          new RandomBndGenerator(n ? p : null, stream),
          new RandomBndGenerator(n ? p : null, n === 2n ? secondary : stream),
        ],
        trace: unknown[] = [];
      const state = () => g.map((x) => [x.p, x.p ? x.nb : null, x.p ? x.mask : null]);
      for (let i = 0; i < data.length; i += 3) {
        const action = data[i],
          slot = Number(data[i + 1]),
          arg = data[i + 2]!;
        let result: unknown = null,
          kind: string | null = null,
          message: string | null = null;
        try {
          if (action === 0n) g[slot]!.build(arg);
          else if (action === 1n) {
            const r: bigint[] = [];
            for (let j = 0; j < Number(arg); j++) r.push(g[slot]!.next());
            result = r;
          } else if (action === 2n) g[slot] = new RandomBndGenerator(g[Number(arg)]!);
          else if (action === 3n) g[slot] = g[Number(arg)]!;
          else if (action === 4n)
            g[slot] = new RandomBndGenerator(null, n === 2n && slot === 1 ? secondary : stream);
          else if (action === 5n) result = hex(stream.get(Number(arg)));
          else if (action === 6n) stream.set_nonce(arg);
          else if (action === 7n) g[slot]!.assign(g[Number(arg)]!);
          else throw new Error('unknown cached sampler command');
        } catch (e) {
          kind = (e as Error).name;
          message = (e as Error).message;
        }
        trace.push([kind, message, result, state()]);
      }
      value = trace;
    } else throw new Error('unknown root prerequisite operation');
  } catch (e) {
    errorType = (e as Error).name;
    error = (e as Error).message;
  }
  return JSON.stringify(
    [
      errorType,
      error,
      value,
      op === 3n && n === 2n ? [hex(stream.get(64)), hex(secondary.get(64))] : hex(stream.get(64)),
    ],
    (_, v) => (typeof v === 'bigint' ? String(v) : v)
  );
}

const parse = (s: string): bigint[] =>
  s
    .slice(1, -1)
    .split(',')
    .filter((x) => x.trim())
    .map((x) => BigInt(x.trim()));
for (const [i, row] of fixtures.entries())
  test('native NTL root prerequisites ' + i, () => {
    let result: string | null = null,
      error: string | null = null,
      errorType: string | null = null;
    try {
      result = ntl_root_prerequisites(
        BigInt(row.args[0]!),
        BigInt(row.args[1]!),
        BigInt(row.args[2]!),
        parse(row.args[3]!),
        parse(row.args[4]!)
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
test('NTL root products preserve inputs and own output storage', () => {
  const roots = [1n, 2n, 0n],
    f = BuildFromRoots(roots, 5n);
  f.fill(0n);
  expect(roots).toEqual([1n, 2n, 0n]);
  expect(BuildFromRoots(roots, 5n)).toEqual([0n, 2n, 2n, 1n]);
});
test('NTL cached sampler assignment preserves identity and copies stream reference', () => {
  const a = new RandomStream(new Uint8Array(32)),
    b = new RandomStream(new Uint8Array(32));
  b.set_nonce(1n);
  const x = new RandomBndGenerator(5n, a),
    y = new RandomBndGenerator(257n, b),
    alias = x,
    control = new RandomStream(b);
  expect(x.assign(y)).toBe(x);
  expect(alias.next()).toBe(VectorRandomBnd(1, 257n, control)[0]!);
  expect(y.next()).toBe(VectorRandomBnd(1, 257n, control)[0]!);
  expect(x.assign(x)).toBe(x);
});
test('NTL cached sampler copies share a stream and retain independent bounds', () => {
  const s = new RandomStream(new Uint8Array(32)),
    control = new RandomStream(s),
    a = new RandomBndGenerator(5n, s),
    b = new RandomBndGenerator(a);
  a.build(257n);
  expect(a.next()).toBe(VectorRandomBnd(1, 257n, control)[0]!);
  expect(b.next()).toBe(VectorRandomBnd(1, 5n, control)[0]!);
  expect(b.p).toBe(5n);
});
