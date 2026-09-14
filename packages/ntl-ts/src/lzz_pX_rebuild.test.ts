import { nativeFixtures as loadLiveNative } from "../../../tests/property/native-live.mjs";
import { test, expect } from 'bun:test';


const fixtures = (await loadLiveNative(import.meta.url, "./lzz_pX_rebuild.native.json.gz")) as { args: string[]; result: string | null; error: string | null; errorType: string | null }[];
import { FFTRoundUp } from './FFT_impl.js';
import { RandomStream } from './ZZ.js';
import {
  FFTFwd_trunc,
  FFTRev1_trunc,
  FFTPrimeContext,
  UseFFTPrime,
  type FFTPrimeInfo,
} from './FFT.js';
import {
  FromModularRep,
  zz_pXMultiplier,
  zz_pXModulus,
  build,
  rem,
  MulMod,
  SqrMod,
} from './lzz_pX.js';
export function ntl_word_rebuild(
  key: bigint[],
  p: bigint,
  maxroot: bigint,
  lengths: bigint[],
  coefficients: bigint[],
  commands: bigint[]
): string {
  const stream = new RandomStream(Uint8Array.from(key, Number)),
    context = new FFTPrimeContext(),
    state = { context, stream },
    out: unknown[] = [];
  let F = new zz_pXModulus(null, p, { maxroot: Number(maxroot), state });
  const B = new zz_pXMultiplier();
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
  for (let i = 0; i < commands.length; i += 2) {
    const op = Number(commands[i]),
      a = polynomials[Number(commands[i + 1])]!;
    let result: unknown = null,
      errorType: string | null = null,
      error: string | null = null;
    try {
      if (op === 0) build(F, a);
      else if (op === 1) result = rem(a, F);
      else if (op === 2) result = MulMod(a, a, F);
      else if (op === 3) result = SqrMod(a, F);
      else if (op === 4) build(B, a, F);
      else if (op === 5) result = MulMod(a, B, F);
      else if (op === 6) result = [B.val(), BigInt(B.UseFFT)];
      else if (op === 7 || op === 8) {
        const [prime, k, yn, xn] = a.map(Number);
        UseFFTPrime(prime!, context, stream);
        result =
          op === 7
            ? FFTFwd_trunc(a.slice(4), k!, context.get(prime!), yn!, xn!)
            : FFTRev1_trunc(a.slice(4), k!, context.get(prime!), yn!);
      } else if (op === 9) result = BigInt(FFTRoundUp(Number(a[0]), Number(a[1])));
      else if (op === 10) F = new zz_pXModulus(null, p, { maxroot: Number(maxroot) });
      else if (op === 12) {
        const rows = Number(a[0]),
          cols = Number(a[1]);
        result = FromModularRep(
          Array.from({ length: rows }, (_, i) => a.slice(2 + i * cols, 2 + (i + 1) * cols)),
          F
        );
      } else throw new Error('unknown quotient rebuild operation');
    } catch (e) {
      errorType = (e as Error).name;
      error = (e as Error).message;
    }
    out.push([
      errorType,
      error,
      result,
      BigInt(F.n),
      F.f.slice(),
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
  test('native NTL word rebuild and truncated transform ' + i, () => {
    let result: string | null = null,
      error: string | null = null,
      errorType: string | null = null;
    try {
      result = ntl_word_rebuild(
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
test('NTL transform arrays are owned and inverse transform restores coefficients', () => {
  const stream = new RandomStream(new Uint8Array(32)),
    context = new FFTPrimeContext();
  UseFFTPrime(0, context, stream);
  const info = context.get(0),
    a = [1n, 2n, 3n, 4n],
    saved = a.slice(),
    freq = FFTFwd_trunc(a, 2, info, 4, 4),
    original = freq.slice();
  const back = FFTRev1_trunc(freq, 2, info, 4);
  expect(back).toEqual(a);
  expect(a).toEqual(saved);
  expect(freq).toEqual(original);
  back[0] = 99n;
  expect(a).toEqual(saved);
});
test('NTL transform numeric adapters reject fractional and non-finite dimensions', () => {
  const context = new FFTPrimeContext();
  UseFFTPrime(0, context, new RandomStream(new Uint8Array(32)));
  const info = context.get(0);
  for (const x of [NaN, Infinity, -Infinity, 1.5]) {
    expect(() => FFTRoundUp(x, 4)).toThrow('FFTRoundUp: length must be a signed native integer');
    expect(() => FFTRoundUp(16, x)).toThrow('FFTRoundUp: exponent must be between 0 and 62');
    expect(() => FFTFwd_trunc([1n], x, info, 1, 1)).toThrow(
      'FFT transform: exponent exceeds the prime root table'
    );
    expect(() => FFTRev1_trunc([1n], 0, info, x)).toThrow(
      'FFT transform: lengths must be admissible'
    );
  }
});

test('NTL CRT preserves coefficient columns and rejects ragged residue rows', () => {
  const F = new zz_pXModulus(null, 65535n),
    rows = [
      [0n, 1n],
      [0n, 1n],
    ],
    saved = rows.map((row) => row.slice());
  const result = FromModularRep(rows, F);
  expect(result).toEqual([0n, 1n]);
  result[0] = 99n;
  expect(rows).toEqual(saved);
  expect(() => FromModularRep([[1n], [2n, 3n]], F)).toThrow(
    'FromModularRep: inconsistent coefficient counts'
  );
});
