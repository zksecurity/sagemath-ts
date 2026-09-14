import { nativeFixtures as loadLiveNative } from "../../../tests/property/native-live.mjs";
import { test, expect } from 'bun:test';


const fixtures = (await loadLiveNative(import.meta.url, "./ZZXFactoring_hensel_state.native.json.gz")) as { args: string[]; result: string | null; error: string | null; errorType: string | null }[];
import { MultiLift, AdditionalLifting } from './ZZXFactoring.js';
import { mul } from './ZZX1.js';
import { zz_pXModulus } from './lzz_pX.js';
import { FFTPrimeContext, type FFTPrimeInfo } from './FFT.js';
import { RandomStream } from './ZZ.js';
export function ntl_hensel_state(
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
      if (op === 0) {
        if (p === 0n) throw new RangeError('MultiLift oracle requires an initial word context');
        const count = Number(a[1]);
        const factors = a.slice(2, 2 + count).map((j) => polynomials[Number(j)]!);
        result = MultiLift(factors, polynomials[Number(a[2 + count])]!, Number(a[0]), p, {
          maxroot: Number(maxroot),
          state,
        });
      } else if (op === 1) {
        const count = Number(a[5]),
          factors = a.slice(6, 6 + count).map((j) => polynomials[Number(j)]!);
        result = AdditionalLifting(
          a[0]!,
          Number(a[1]),
          factors,
          a[2]!,
          Number(a[3]),
          polynomials[Number(a[6 + count])]!,
          Boolean(a[4]),
          state
        );
      } else if (op === 2)
        result = mul(polynomials[Number(a[0])]!, polynomials[Number(a[1])]!, state);
      else throw new Error('unknown Hensel state operation');
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
  test('native NTL Hensel state ' + i, () => {
    let result: string | null = null,
      error: string | null = null,
      errorType: string | null = null;
    try {
      result = ntl_hensel_state(
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
  });

test('Hensel results and failed additional lifting retain independent inputs', () => {
  const factors = [
      [1n, 1n],
      [2n, 1n],
    ],
    target = [7n, 3n, 1n];
  const lifted = MultiLift(factors, target, 4, 5n);
  lifted[0]!.push(9n);
  expect(factors).toEqual([
    [1n, 1n],
    [2n, 1n],
  ]);
  expect(target).toEqual([7n, 3n, 1n]);
  const constants = [[1n], [1n]];
  const state = { context: new FFTPrimeContext(), stream: new RandomStream(new Uint8Array(32)) };
  expect(() => AdditionalLifting(5n, 1, constants, 5n, 3, [1n], false, state)).toThrow(
    'build: deg(f) must be at least 1'
  );
  expect(constants).toEqual([[1n], [1n]]);
  expect(state.context.length()).toBe(1);
});
