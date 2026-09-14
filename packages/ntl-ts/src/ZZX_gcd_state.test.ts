import { nativeFixtures as loadLiveNative } from "../../../tests/property/native-live.mjs";
import { test, expect } from 'bun:test';


const fixtures = (await loadLiveNative(import.meta.url, "./ZZX_gcd_state.native.json.gz")) as { args: string[]; result: string | null; error: string | null; errorType: string | null }[];
import { SquareFreeDecomp } from './ZZXFactoring.js';
import { GCD, mul } from './ZZX1.js';
import { zz_pXModulus } from './lzz_pX.js';
import { FFTPrimeContext, type FFTPrimeInfo } from './FFT.js';
import { RandomStream } from './ZZ.js';
export function ntl_integer_gcd_state(
  key: bigint[],
  p: bigint,
  maxroot: bigint,
  lengths: bigint[],
  coefficients: bigint[],
  commands: bigint[]
): string {
  const stream = new RandomStream(Uint8Array.from(key, Number)),
    context = new FFTPrimeContext(),
    state = { context, stream };
  if (p !== 0n) new zz_pXModulus(null, p, { maxroot: Number(maxroot), state });
  let offset = 0;
  const polynomials = lengths.map((n) => {
    const a = coefficients.slice(offset, offset + Number(n));
    offset += Number(n);
    return a;
  });
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
  const out: unknown[] = [];
  for (let i = 0; i < commands.length; i += 2) {
    const op = Number(commands[i]),
      a = polynomials[Number(commands[i + 1])]!;
    let result: unknown = null,
      errorType: string | null = null,
      error: string | null = null;
    try {
      if (op === 0) result = GCD(polynomials[Number(a[0])]!, polynomials[Number(a[1])]!, state);
      else if (op === 1) result = SquareFreeDecomp(polynomials[Number(a[0])]!, state);
      else if (op === 2)
        result = mul(polynomials[Number(a[0])]!, polynomials[Number(a[1])]!, state);
      else throw new Error('unknown integer GCD state operation');
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
  test('native NTL integer GCD state ' + i, () => {
    let result: string | null = null,
      error: string | null = null,
      errorType: string | null = null;
    try {
      result = ntl_integer_gcd_state(
        parse(row.args[0]!),
        BigInt(row.args[1]!),
        BigInt(row.args[2]!),
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

test('integer GCD and squarefree outputs own their coefficients', () => {
  const a = [-2n, 0n, 2n],
    b = [2n, -4n, 2n],
    f = [1n, 4n, 4n];
  const state = { context: new FFTPrimeContext(), stream: new RandomStream(new Uint8Array(32)) };
  const g = GCD(a, b, state),
    factors = SquareFreeDecomp(f, state);
  expect(g).toEqual([-2n, 2n]);
  expect(factors).toEqual([[[1n, 2n], 2]]);
  g[0] = 0n;
  factors[0]![0][0] = 0n;
  expect([a, b, f]).toEqual([
    [-2n, 0n, 2n],
    [2n, -4n, 2n],
    [1n, 4n, 4n],
  ]);
});
