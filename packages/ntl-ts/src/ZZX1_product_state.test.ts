import { nativeFixtures as loadLiveNative } from "../../../tests/property/native-live.mjs";
import { test, expect } from 'bun:test';


const fixtures = (await loadLiveNative(import.meta.url, "./ZZX1_product_state.native.json.gz")) as { args: string[]; result: string | null; error: string | null; errorType: string | null }[];
import { mul, sqr, ChooseSS } from './ZZX1.js';
import { RandomStream } from './ZZ.js';
import { FFTPrimeContext, UseFFTPrime, type FFTPrimeInfo } from './FFT.js';
export function ntl_integer_product_state(
  key: bigint[],
  a: bigint[],
  b: bigint[],
  commands: bigint[]
): string {
  const stream = new RandomStream(Uint8Array.from(key, Number)),
    context = new FFTPrimeContext(),
    state = { context, stream },
    out: unknown[] = [];
  const hex = (v: Uint8Array) => Array.from(v, (x) => x.toString(16).padStart(2, '0')).join('');
  const bits = (v: number) => {
    const d = new DataView(new ArrayBuffer(8));
    d.setFloat64(0, v, false);
    return d.getBigUint64(0, false).toString(16).padStart(16, '0');
  };
  const info = (v: FFTPrimeInfo) => [
    v.q,
    bits(v.qrecip),
    v.RootTable[0],
    v.RootTable[1],
    v.TwoInvTable,
  ];
  for (let i = 0; i < commands.length; i += 2) {
    const op = Number(commands[i]),
      p = commands[i + 1]!;
    let result: unknown = null,
      errorType: string | null = null,
      error: string | null = null;
    try {
      if (op === 0) result = mul(a, b, state);
      else if (op === 1) result = mul(a, a, state);
      else if (op === 2) result = sqr(a, state);
      else if (op === 3) result = mul(b, a, state);
      else if (op === 4 || op === 7) UseFFTPrime(Number(p), context, stream);
      else if (op === 5) stream.set_nonce(p);
      else if (op === 6) result = hex(stream.get(Number(p)));
      else if (op === 8) {
        if (p <= 1n) throw new Error('ZZ_pContext: p must be > 1');
      } else if (op === 9)
        result = ChooseSS(
          Number(a[0] ?? 0n),
          Number(a[1] ?? 0n),
          Number(a[2] ?? 0n),
          Number(a[3] ?? 0n)
        );
      else throw new Error('unknown integer polynomial state operation');
    } catch (e) {
      errorType = (e as Error).name;
      error = (e as Error).message;
    }
    out.push([
      errorType,
      error,
      result,
      Array.from({ length: context.length() }, (_, j) => info(context.get(j))),
      hex(new RandomStream(stream).get(64)),
    ]);
  }
  return JSON.stringify(out, (_, v) => (typeof v === 'bigint' ? String(v) : v));
}

const parse = (s: string): bigint[] =>
  s.slice(1, -1).trim()
    ? s
        .slice(1, -1)
        .split(',')
        .map((x) => BigInt(x.trim()))
    : [];
for (const [i, row] of fixtures.entries())
  test('native NTL integer product state and transform choice ' +
    i +
    ((i >= 1485 && i < 1490) || i >= 1495 ? ' — explicit native arithmetic guard' : ''), () => {
    let result: string | null = null,
      error: string | null = null,
      errorType: string | null = null;
    try {
      result = ntl_integer_product_state(
        parse(row.args[0]!),
        parse(row.args[1]!),
        parse(row.args[2]!),
        parse(row.args[3]!)
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
test('NTL integer products own coefficients and isolate explicit caches and streams', () => {
  const state = { context: new FFTPrimeContext(), stream: new RandomStream(new Uint8Array(32)) },
    other = { context: new FFTPrimeContext(), stream: new RandomStream(new Uint8Array(32)) },
    a = Array<bigint>(150).fill(1n),
    b = a.slice(),
    snapshot = structuredClone([a, b]),
    next = new RandomStream(other.stream).get(64);
  const product = mul(a, b, state),
    square = sqr(a, state);
  product[0] = 99n;
  square[1] = 99n;
  expect([a, b]).toEqual(snapshot);
  expect(state.context.length()).toBe(1);
  expect(other.context.length()).toBe(0);
  expect(new RandomStream(other.stream).get(64)).toEqual(next);
});
test('NTL ChooseSS rejects parameters not representable as native integer words', () => {
  for (const value of [NaN, Infinity, -Infinity, 1.5])
    expect(() => ChooseSS(79, value, 79, 1920)).toThrow(
      'SSRatio arguments must be signed native integers'
    );
});
