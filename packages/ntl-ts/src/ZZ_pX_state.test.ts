import { nativeFixtures as loadLiveNative } from "../../../tests/property/native-live.mjs";
import { test, expect } from 'bun:test';


const fixtures = (await loadLiveNative(import.meta.url, "./ZZ_pX_state.native.json.gz")) as { args: string[]; result: string | null; error: string | null; errorType: string | null }[];
import { mul } from './ZZ_pX.js';
import { mul as treeMul } from './ZZXFactoring.js';
import { SSRatio } from './ZZX1.js';
import { RandomStream } from './ZZ.js';
import { FFTPrimeContext, UseFFTPrime, type FFTPrimeInfo } from './FFT.js';
export function ntl_polynomial_state(
  key: bigint[],
  p: bigint,
  lengths: bigint[],
  coefficients: bigint[],
  indices: bigint[],
  commands: bigint[]
): string {
  let offset = 0;
  const W = lengths.map((n) => {
    const w = coefficients.slice(offset, offset + Number(n));
    offset += Number(n);
    return w;
  });
  if (p <= 1n) throw new Error('ZZ_pContext: p must be > 1');
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
  for (let i = 0; i < commands.length; i += 2) {
    const op = Number(commands[i]),
      a = commands[i + 1]!;
    let result: unknown = null,
      error: string | null = null,
      errorType: string | null = null;
    try {
      if (op === 0) result = mul(W[0]!, W[1]!, p, { context, stream });
      else if (op === 1) result = mul(W[0]!, W[0]!, p, { context, stream });
      else if (op === 2) result = treeMul(W, p, undefined, { context, stream });
      else if (op === 3) result = treeMul(W, p, indices.map(Number), { context, stream });
      else if (op === 4) UseFFTPrime(Number(a), context, stream);
      else if (op === 5) {
        if (a <= 1n) throw new Error('ZZ_pContext: p must be > 1');
        p = a;
      } else if (op === 6) result = hex(stream.get(Number(a)));
      else if (op === 7) stream.set_nonce(a);
      else if (op === 8) {
        const v = W[Number(a)]!;
        result = bits(
          SSRatio(Number(v[0] ?? 0n), Number(v[1] ?? 0n), Number(v[2] ?? 0n), Number(v[3] ?? 0n))
        );
      } else throw new Error('unknown polynomial state operation');
    } catch (e) {
      error = (e as Error).message;
      errorType = (e as Error).name;
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
  test('native NTL polynomial state ' +
    i +
    (i >= 667 && i <= 673 ? ' — explicit adapter guard' : ''), () => {
    let result: string | null = null,
      error: string | null = null,
      errorType: string | null = null;
    try {
      result = ntl_polynomial_state(
        parse(row.args[0]!),
        BigInt(row.args[1]!),
        parse(row.args[2]!),
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
test('NTL stateful products own coefficients and leave factor lists unchanged', () => {
  const p = (1n << 127n) - 1n,
    A = Array<bigint>(80).fill(1n),
    B = A.slice(),
    context = new FFTPrimeContext(),
    stream = new RandomStream(new Uint8Array(32));
  const expected = mul(A, B, p);
  const product = mul(A, A, p, { context, stream });
  expect(product).toEqual(expected);
  product[0] = 0n;
  expect(A[0]).toBe(1n);
  expect(B[0]).toBe(1n);
  const W = [A, B],
    snapshot = structuredClone(W),
    I = [1, 0, 1];
  const factors = treeMul(W, p, I, { context, stream });
  factors[0] = 0n;
  expect(W).toEqual(snapshot);
  expect(I).toEqual([1, 0, 1]);
});
test('NTL SSRatio rejects non-native numeric inputs', () => {
  for (const value of [NaN, Infinity, -Infinity, 1.5])
    expect(() => SSRatio(value, 1, 1, 1)).toThrow(
      'SSRatio arguments must be signed native integers'
    );
});
