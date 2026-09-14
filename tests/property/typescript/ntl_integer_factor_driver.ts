import { _ZZ_pX_euclidean_kernels } from '../../../packages/ntl-ts/src/ZZ_pX1.js';
import {
  factor,
  SFFactor,
  ll_SFFactor,
  SmallPrimeFactorization,
  LocalInfoT,
  MultiLift,
  FindTrueFactors_vH,
} from '../../../packages/ntl-ts/src/ZZXFactoring.js';
import { zz_pXModulus } from '../../../packages/ntl-ts/src/lzz_pX.js';
import { FFTPrimeContext, type FFTPrimeInfo } from '../../../packages/ntl-ts/src/FFT.js';
import { RandomStream } from '../../../packages/ntl-ts/src/ZZ.js';
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
