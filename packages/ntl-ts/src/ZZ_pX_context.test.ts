import { nativeFixtures as loadLiveNative } from "../../../tests/property/native-live.mjs";
import { test, expect } from 'bun:test';


const fixtures = (await loadLiveNative(import.meta.url, "./ZZ_pX_context.native.json.gz")) as { args: string[]; result: string | null; error: string | null; errorType: string | null }[];
import { RandomStream } from './ZZ.js';
import { FFTPrimeContext, type FFTPrimeInfo } from './FFT.js';
import {
  ZZ_pXModulus,
  ZZ_pXMultiplier,
  build,
  rem,
  MulMod,
  SqrMod,
  FromModularRep,
} from './ZZ_pX.js';
export function ntl_big_quotient(
  key: bigint[],
  p: bigint,
  lengths: bigint[],
  coefficients: bigint[],
  commands: bigint[]
): string {
  const context = new FFTPrimeContext(),
    stream = new RandomStream(Uint8Array.from(key, Number)),
    state = { context, stream };
  let F = new ZZ_pXModulus(null, p, state);
  let B = new ZZ_pXMultiplier();
  const trace: unknown[] = [];
  let offset = 0;
  const polynomials = lengths.map((n) => {
    const a = coefficients.slice(offset, offset + Number(n));
    offset += Number(n);
    return a;
  });
  const bits = (v: number) => {
    const view = new DataView(new ArrayBuffer(8));
    view.setFloat64(0, v, false);
    return view.getBigUint64(0, false).toString(16).padStart(16, '0');
  };
  const info = (v: FFTPrimeInfo) => [
    v.q,
    bits(v.qrecip),
    v.RootTable[0],
    v.RootTable[1],
    v.TwoInvTable,
  ];
  const hex = (v: Uint8Array) => Array.from(v, (x) => x.toString(16).padStart(2, '0')).join('');
  for (let i = 0; i < commands.length; i += 2) {
    const op = Number(commands[i]),
      a = polynomials[Number(commands[i + 1])]!;
    let result: unknown = null,
      error: string | null = null,
      errorType: string | null = null;
    try {
      if (op === 0) build(F, a);
      else if (op === 1) result = rem(a, F);
      else if (op === 2) result = MulMod(a, a, F);
      else if (op === 3) result = SqrMod(a, F);
      else if (op === 4) build(B, a, F);
      else if (op === 5) result = MulMod(a, B, F);
      else if (op === 6) result = [B.val(), BigInt(B.UseFFT)];
      else if (op === 7) {
        const rows = Number(a[0]),
          cols = Number(a[1]);
        result = FromModularRep(
          Array.from({ length: rows }, (_, j) => a.slice(2 + j * cols, 2 + (j + 1) * cols)),
          F
        );
      } else if (op === 8) F = new ZZ_pXModulus(null, p, state);
      else if (op === 9) result = F.val();
      else if (op === 10) new ZZ_pXModulus(a.slice(1), a[0]!, state);
      else if (op === 11) B = new ZZ_pXMultiplier(a, F);
      else if (op === 12) F = new ZZ_pXModulus(a, p, state);
      else if (op === 13) {
        const length = Number(a[1]);
        const other = new ZZ_pXModulus(a.slice(2, 2 + length), a[0]!, state);
        build(B, a.slice(2 + length), other);
      } else if (op === 15) {
        const length = Number(a[0]);
        result = MulMod(a.slice(1, length + 1), a.slice(length + 1), F);
      } else if (op === 16) result = MulMod(a, a.slice(), F);
      else if (op === 17) {
        const rows = Number(a[0]);
        let offset = 1 + rows;
        const residues = Array.from({ length: rows }, (_, j) => {
          const length = Number(a[j + 1]),
            row = a.slice(offset, offset + length);
          offset += length;
          return row;
        });
        result = FromModularRep(residues, F);
      } else throw new Error('unknown quotient rebuild operation');
    } catch (e) {
      error = (e as Error).message;
      errorType = (e as Error).name;
    }
    trace.push([
      errorType,
      error,
      result,
      BigInt(F.n),
      F.f.slice(),
      BigInt(F.UseFFT),
      Array.from({ length: context.length() }, (_, j) => info(context.get(j))),
      hex(new RandomStream(stream).get(64)),
    ]);
  }
  return JSON.stringify(trace, (_, v) => (typeof v === 'bigint' ? String(v) : v));
}

const parse = (s: string): bigint[] =>
  s.slice(1, -1).trim()
    ? s
        .slice(1, -1)
        .split(',')
        .map((x) => BigInt(x.trim()))
    : [];
for (const [i, row] of fixtures.entries())
  test('native NTL arbitrary-modulus multiplier context ' + i, () => {
    let result: string | null = null,
      error: string | null = null,
      errorType: string | null = null;
    try {
      result = ntl_big_quotient(
        parse(row.args[0]!),
        BigInt(row.args[1]!),
        parse(row.args[2]!),
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
  });
