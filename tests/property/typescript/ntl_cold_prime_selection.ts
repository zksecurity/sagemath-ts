import {
  FFTPrimeContext,
  UseFFTPrime,
  type FFTPrimeInfo,
} from '../../../packages/ntl-ts/src/FFT.js';
import { LocalInfoT, SmallPrimeFactorization } from '../../../packages/ntl-ts/src/ZZXFactoring.js';
import { RandomStream } from '../../../packages/ntl-ts/src/ZZ.js';
import { zz_pXModulus } from '../../../packages/ntl-ts/src/lzz_pX.js';
export function ntl_cold_prime_selection(
  commands: bigint[],
  key: bigint[],
  capacity: bigint,
  warm: bigint
): string {
  const fresh = () => {
    const I = new LocalInfoT();
    I.p = Array<bigint>(Number(capacity)).fill(-17n);
    I.p = [];
    return I;
  };
  let I = fresh(),
    pos = 0;
  const stream = new RandomStream(Uint8Array.from(key, Number)),
    out: unknown[] = [];
  const cache = new FFTPrimeContext();
  if (warm > 0n) {
    const fallback = new RandomStream(new Uint8Array(32));
    for (let i = 0; i < Number(warm); i++) UseFFTPrime(i, cache, fallback);
  }
  const hex = (v: Uint8Array) => Array.from(v, (x) => x.toString(16).padStart(2, '0')).join('');
  const bits = (v: number) => {
    const d = new DataView(new ArrayBuffer(8));
    d.setFloat64(0, v, false);
    return d.getBigUint64(0, false).toString(16).padStart(16, '0');
  };
  const fftInfo = (v: FFTPrimeInfo) => [
    v.q,
    bits(v.qrecip),
    v.RootTable[0],
    v.RootTable[1],
    v.TwoInvTable,
  ];
  const scalar = () => commands[pos++]!,
    word = () => Number(scalar()),
    vector = () => {
      const n = word(),
        v = commands.slice(pos, pos + n);
      pos += n;
      return v;
    };
  while (pos < commands.length) {
    const op = word();
    let result: unknown = null,
      errorType: string | null = null,
      error: string | null = null;
    try {
      if (op === 0) {
        const InitNumPrimes = word(),
          MaxNumPrimes = word(),
          f = vector();
        result = SmallPrimeFactorization(I, f, stream, {
          InitNumPrimes,
          MaxNumPrimes,
          context: cache,
        });
      } else if (op === 9)
        result = SmallPrimeFactorization(I, vector(), stream, { context: cache });
      else if (op === 10) {
        const InitNumPrimes = word();
        result = SmallPrimeFactorization(I, vector(), stream, { InitNumPrimes, context: cache });
      } else if (op === 11) {
        const MaxNumPrimes = word();
        result = SmallPrimeFactorization(I, vector(), stream, { MaxNumPrimes, context: cache });
      } else if (op === 1) I.s.reset(scalar());
      else if (op === 2) I.PossibleDegrees = scalar();
      else if (op === 3) I.p = vector();
      else if (op === 4) {
        const n = word();
        I.pattern = Array.from({ length: n }, () => vector().map(Number));
      } else if (op === 5) I = fresh();
      else if (op === 6) result = Array.from({ length: word() }, () => I.s.next());
      else if (op === 7) {
        const p = scalar(),
          maxroot = word();
        new zz_pXModulus(null, p, { maxroot, state: { context: cache, stream } });
        I.context = Object.freeze({ p, maxroot });
      } else if (op === 8) {
        I.n = word();
        I.NumPrimes = word();
        I.NumFactors = word();
      } else if (op === 12) UseFFTPrime(word(), cache, stream);
      else if (op === 13) stream.set_nonce(scalar());
      else if (op === 14) result = hex(stream.get(word()));
      else if (op === 16) {
        const p = scalar(),
          maxroot = word(),
          f = vector();
        new zz_pXModulus(null, p, { maxroot, state: { context: cache, stream } });
        I.context = Object.freeze({ p, maxroot });
        const F = new zz_pXModulus(f, p, { maxroot, state: { context: cache, stream } });
        result = [F.n];
      } else throw new Error('unknown small-prime operation');
    } catch (e) {
      errorType = (e as Error).name;
      error = (e as Error).message;
    }
    let context: unknown = null;
    if (I.context) {
      const F = new zz_pXModulus(null, I.context.p, { maxroot: I.context.maxroot });
      context = [I.context.p, F.PrimeCnt, F.MaxRoot];
    }
    const state = [I.n, I.NumPrimes, I.NumFactors, I.PossibleDegrees, I.p, I.pattern, context];
    const tail = Array.from(new RandomStream(stream).get(64), (x) =>
      x.toString(16).padStart(2, '0')
    ).join('');
    out.push([
      errorType,
      error,
      result,
      state,
      tail,
      Array.from({ length: cache.length() }, (_, j) => fftInfo(cache.get(j))),
    ]);
  }
  return JSON.stringify(out, (_, v) => (typeof v === 'bigint' ? String(v) : v));
}
