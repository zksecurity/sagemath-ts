import { nativeFixtures as loadLiveNative } from "../../../tests/property/native-live.mjs";
import { test, expect } from 'bun:test';


const fixtures = (await loadLiveNative(import.meta.url, "./ZZX_factor_driver.native.json.gz")) as { args: string[]; result: string | null; error: string | null; errorType: string | null }[];
import { _ZZ_pX_euclidean_kernels } from './ZZ_pX1.js';
import {
  factor,
  SFFactor,
  ll_SFFactor,
  SmallPrimeFactorization,
  LocalInfoT,
  MultiLift,
  FindTrueFactors_vH,
} from './ZZXFactoring.js';
import { zz_pXModulus } from './lzz_pX.js';
import { FFTPrimeContext, type FFTPrimeInfo } from './FFT.js';
import { RandomStream } from './ZZ.js';
export function ntl_integer_factor_driver(
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
      const input = polynomials[Number(a[0])]!,
        bound = Number(a[1]);
      const options = {
        state,
        PowerHack: Number(a[2]),
        van_Hoeij: Number(a[3]),
        InitNumPrimes: Number(a[4]),
        MaxNumPrimes: Number(a[5]),
        MaxPrune: Number(a[6]),
        ok_to_abandon: Number(a[7] ?? 0n),
      };
      if (op === 0) result = factor(input, bound, options);
      else if (op === 1) result = SFFactor(input, bound, options);
      else if (op === 2) result = ll_SFFactor(input, bound, options);
      else if (op === 3) {
        const info = new LocalInfoT();
        if (a[9]) info.s.reset(a[9]);
        const small = SmallPrimeFactorization(info, input, stream, {
          InitNumPrimes: options.InitNumPrimes,
          MaxNumPrimes: options.MaxNumPrimes,
          context,
        });
        if (small === null) result = null;
        else {
          const p = info.context!.p;
          let exponent = Number(a[8]),
            modulus: bigint,
            factorBound = bound;
          if (exponent === 0) {
            const bits = (v: bigint) => (v === 0n ? 0 : (v < 0n ? -v : v).toString(2).length);
            const n = input.length - 1,
              bnd1 =
                input.reduce((best, c) => Math.max(best, bits(c)), 0) +
                Math.trunc((bits(BigInt(n + 1)) + 1) / 2);
            if (!factorBound || bnd1 < factorBound) factorBound = bnd1;
            let i = Math.trunc(n / 2);
            while (!((info.PossibleDegrees >> BigInt(i)) & 1n)) i--;
            const lc = bits(input.at(-1)!),
              coeffBound = factorBound + lc + i;
            const liftBound =
              Math.max(
                coeffBound + 15,
                factorBound + lc + 2 * bits(BigInt(n)) + 64,
                lc + bits(input.reduce((a, b) => a + b, 0n))
              ) + 2;
            modulus = 1n;
            while (bits(modulus) <= liftBound) {
              modulus *= p;
              exponent++;
            }
            factorBound = coeffBound;
          } else modulus = p ** BigInt(exponent);
          const ar = _ZZ_pX_euclidean_kernels(modulus),
            inverse = ar.inverse(ar.mod(input.at(-1)!));
          const target = input.map((c) => ar.mod(c * inverse));
          const lifted = MultiLift(small, target, exponent, p, {
            maxroot: info.context!.maxroot,
            state,
          });
          result = FindTrueFactors_vH(
            input,
            lifted,
            modulus,
            p,
            exponent,
            info,
            factorBound,
            options
          );
        }
      } else throw new Error('unknown integer factor driver operation');
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
  test('native NTL integer factor driver ' + i, () => {
    let result: string | null = null,
      error: string | null = null,
      errorType: string | null = null;
    try {
      result = ntl_integer_factor_driver(
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

test('integer factor drivers own their output coefficients', () => {
  const f = [-2n, -4n, -2n],
    primitive = [1n, 3n, 2n];
  const [content, factors] = factor(f),
    squarefree = SFFactor(primitive),
    plain = ll_SFFactor(primitive);
  expect([content, factors]).toEqual([-2n, [[[1n, 1n], 2]]]);
  expect(squarefree).toEqual([
    [1n, 1n],
    [1n, 2n],
  ]);
  expect(plain).toEqual(squarefree);
  factors[0]![0][0] = 99n;
  squarefree[0]![0] = 99n;
  plain[0]![0] = 99n;
  expect(f).toEqual([-2n, -4n, -2n]);
  expect(primitive).toEqual([1n, 3n, 2n]);
});
