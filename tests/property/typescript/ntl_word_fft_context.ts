import { UpdateMap, XGCD, GCD, MinPolySeq } from '../../../packages/ntl-ts/src/lzz_pX1.js';
import { FFTRoundUp } from '../../../packages/ntl-ts/src/FFT_impl.js';
import { RandomStream } from '../../../packages/ntl-ts/src/ZZ.js';
import {
  FFTFwd_trunc,
  FFTFwd_trans,
  FFTRev1_trans,
  FFTRev1_trunc,
  FFTPrimeContext,
  UseFFTPrime,
  GetFFTPrime,
  type FFTPrimeInfo,
} from '../../../packages/ntl-ts/src/FFT.js';
import {
  FromModularRep,
  zz_pXMultiplier,
  zz_pXModulus,
  build,
  rem,
  MulMod,
  SqrMod,
  mul,
  sqr,
  InvMod,
  InvModStatus,
  PowerXMod,
  PowerXPlusAMod,
  PowerMod,
} from '../../../packages/ntl-ts/src/lzz_pX.js';
export function ntl_word_fft_context(
  key: bigint[],
  fftindex: bigint,
  lengths: bigint[],
  coefficients: bigint[],
  commands: bigint[]
): string {
  const stream = new RandomStream(Uint8Array.from(key, Number)),
    context = new FFTPrimeContext(),
    state = { context, stream },
    out: unknown[] = [];
  if (fftindex < 0n) throw new Error('bad FFT prime index');
  UseFFTPrime(Number(fftindex), context, stream);
  const p = GetFFTPrime(Number(fftindex), context);
  const options = { fftPrime: Number(fftindex), state };
  let F = new zz_pXModulus(null, p, options);
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
      else if (op === 10) F = new zz_pXModulus(null, p, options);
      else if (op === 11) result = UpdateMap(a, B, F);
      else if (op === 13 || op === 14) {
        const prime = Number(a[0]),
          k = Number(a[1]);
        UseFFTPrime(prime, context, stream);
        result = (op === 13 ? FFTFwd_trans : FFTRev1_trans)(a.slice(4), k, context.get(prime));
      } else if (op === 15) {
        const length = Number(a[2]);
        const other = new zz_pXModulus(a.slice(3, 3 + length), a[0]!, {
          maxroot: Number(a[1]),
          state,
        });
        build(B, a.slice(3 + length), other);
      } else if (op >= 16 && op <= 25) {
        const length = Number(a[0]),
          aa = a.slice(1, 1 + length),
          bb = a.slice(1 + length);
        if (op === 16) result = XGCD(aa, bb, p, options);
        else if (op === 17) result = mul(aa, bb, p, options);
        else if (op === 18) result = sqr(aa, p, options);
        else if (op === 19) result = mul(aa, aa, p, options);
        else if (op === 20) result = InvMod(aa, bb, p, options);
        else if (op === 21) result = InvModStatus(aa, bb, p, options);
        else if (op === 22) result = PowerXMod(a[1]!, F);
        else if (op === 23) result = PowerXPlusAMod(a[1]!, a[2]!, F);
        else if (op === 24) result = PowerMod(aa, a[1 + length]!, F);
        else result = GCD(aa, bb, p, options);
      } else if (op === 26)
        result = [F.p, BigInt(F.PrimeCnt), BigInt(F.NumPrimes), BigInt(F.MaxRoot)];
      else if (op === 27) result = MinPolySeq(a.slice(1), Number(a[0]), p, options);
      else throw new Error('unknown quotient rebuild operation');
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
