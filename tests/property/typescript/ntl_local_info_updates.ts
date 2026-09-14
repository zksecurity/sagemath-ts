import {
  LocalInfoT,
  SmallPrimeFactorization,
  UpdateLocalInfo,
} from '../../../packages/ntl-ts/src/ZZXFactoring.js';
import { RandomStream } from '../../../packages/ntl-ts/src/ZZ.js';
import { zz_pXModulus } from '../../../packages/ntl-ts/src/lzz_pX.js';
export function ntl_local_info_updates(
  commands: bigint[],
  key: bigint[],
  capacity: bigint
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
        result = SmallPrimeFactorization(I, f, stream, { InitNumPrimes, MaxNumPrimes });
      } else if (op === 9) result = SmallPrimeFactorization(I, vector(), stream);
      else if (op === 10) {
        const InitNumPrimes = word();
        result = SmallPrimeFactorization(I, vector(), stream, { InitNumPrimes });
      } else if (op === 11) {
        const MaxNumPrimes = word();
        result = SmallPrimeFactorization(I, vector(), stream, { MaxNumPrimes });
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
        new zz_pXModulus(null, p, { maxroot });
        I.context = Object.freeze({ p, maxroot });
      } else if (op === 8) {
        I.n = word();
        I.NumPrimes = word();
        I.NumFactors = word();
      } else if (op >= 12 && op <= 15) {
        const options =
          op === 13
            ? { van_Hoeij: word(), MaxNumPrimes: word() }
            : op === 14
              ? { van_Hoeij: word() }
              : op === 15
                ? { MaxNumPrimes: word() }
                : undefined;
        const modulus = scalar(),
          count = word();
        const W = Array.from({ length: word() }, () => vector()),
          factors = Array.from({ length: word() }, () => vector()),
          f = vector();
        result = UpdateLocalInfo(I, W, factors, f, count, modulus, options);
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
    out.push([errorType, error, result, state, tail]);
  }
  return JSON.stringify(out, (_, v) => (typeof v === 'bigint' ? String(v) : v));
}
