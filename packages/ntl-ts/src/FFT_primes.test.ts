import { nativeFixtures as loadLiveNative } from "../../../tests/property/native-live.mjs";


import { test, expect } from 'bun:test';
const fixtures = (await loadLiveNative(import.meta.url, "./FFT_primes.native.json.gz")) as { args: (string)[]; error: null; errorType: null; function: string; result: string; seed: number }[];
import { RandomStream } from './ZZ.js';
import {
  FFTPrimeContext,
  IsFFTPrime,
  NextFFTPrime,
  CalcMaxRoot,
  InitFFTPrimeInfo,
  UseFFTPrime,
  GetFFTPrime,
  GetFFTPrimeRecip,
  type FFTPrimeInfo,
} from './FFT.js';
export function ntl_fft_primes(key: bigint[], commands: bigint[]): string {
  const stream = new RandomStream(Uint8Array.from(key, Number)),
    context = new FFTPrimeContext(),
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
  for (let i = 0; i < commands.length; i += 4) {
    const op = Number(commands[i]),
      a = commands[i + 1]!,
      b = commands[i + 2]!;
    let result: unknown = null,
      error: string | null = null,
      errorType: string | null = null;
    try {
      if (op === 0) result = IsFFTPrime(a, stream, b);
      else if (op === 11) result = IsFFTPrime(a, stream);
      else if (op === 1) result = CalcMaxRoot(a);
      else if (op === 2) result = NextFFTPrime(Number(a), context, stream);
      else if (op === 3) {
        UseFFTPrime(Number(a), context, stream);
        result = null;
      } else if (op === 4) result = GetFFTPrime(Number(a), context);
      else if (op === 5) result = bits(GetFFTPrimeRecip(Number(a), context));
      else if (op === 6) result = info(context.get(Number(a)));
      else if (op === 7) result = info(InitFFTPrimeInfo(a, b));
      else if (op === 8) result = context.length();
      else if (op === 9) result = hex(stream.get(Number(a)));
      else if (op === 10) stream.set_nonce(a);
      else throw new Error('unknown FFT-prime operation');
    } catch (e) {
      error = (e as Error).message;
      errorType = (e as Error).name;
    }
    out.push([errorType, error, result, context.length(), hex(new RandomStream(stream).get(64))]);
  }
  return JSON.stringify(out, (_, v) => (typeof v === 'bigint' ? String(v) : v));
}

const parse = (s: string): bigint[] =>
  s
    .slice(1, -1)
    .split(',')
    .filter((x) => x.trim())
    .map((x) => BigInt(x.trim()));
for (const [i, row] of fixtures.entries())
  test('native NTL FFT prime cache ' +
    i +
    (i === 3254 ? ' — special-factor rejection after probabilistic tests' : ''), () => {
    let result: string | null = null,
      error: string | null = null,
      errorType: string | null = null;
    try {
      result = ntl_fft_primes(parse(row.args[0]!), parse(row.args[1]!));
    } catch (e) {
      error = (e as Error).message;
      errorType = (e as Error).name;
    }
    expect({ result, error, errorType }).toEqual({
      result: row.result,
      error: row.error,
      errorType: row.errorType,
    });
  });

test('NTL FFT cache entries and root tables are immutable and stable across growth', () => {
  const context = new FFTPrimeContext(),
    stream = new RandomStream(new Uint8Array(32));
  UseFFTPrime(0, context, stream);
  const first = context.get(0);
  expect(Object.isFrozen(first)).toBe(true);
  expect(Object.isFrozen(first.RootTable)).toBe(true);
  expect(Object.isFrozen(first.RootTable[0])).toBe(true);
  expect(Object.isFrozen(first.RootTable[1])).toBe(true);
  expect(Object.isFrozen(first.TwoInvTable)).toBe(true);
  expect(() => {
    (first.RootTable[0] as bigint[])[0] = 99n;
  }).toThrow(TypeError);
  UseFFTPrime(1, context, stream);
  expect(context.get(0)).toBe(first);
  expect(first.RootTable[0][0]).toBe(1n);
});
test('NTL explicit FFT caches initialize independently with caller-controlled streams', () => {
  const a = new FFTPrimeContext(),
    b = new FFTPrimeContext(),
    x = new RandomStream(new Uint8Array(32)),
    y = new RandomStream(new Uint8Array(32));
  UseFFTPrime(0, a, x);
  expect(b.length()).toBe(0);
  UseFFTPrime(0, b, y);
  expect(a.get(0)).toEqual(b.get(0));
  expect(x.get(64)).toEqual(y.get(64));
  UseFFTPrime(1, a, x);
  expect(b.length()).toBe(1);
});
