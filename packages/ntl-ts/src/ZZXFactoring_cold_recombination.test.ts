import { nativeFixtures as loadLiveNative } from "../../../tests/property/native-live.mjs";
import { test, expect } from 'bun:test';


const fixtures = (await loadLiveNative(import.meta.url, "./ZZXFactoring_cold_recombination.native.json.gz")) as { args: string[]; result: string | null; error: string | null; errorType: string | null }[];
import { FFTPrimeContext, UseFFTPrime, type FFTPrimeInfo } from './FFT.js';
import {
  LocalInfoT,
  InvMul,
  SmallPrimeFactorization,
  UpdateLocalInfo,
  CardinalitySearch,
  CardinalitySearch1,
  FindTrueFactors,
} from './ZZXFactoring.js';
import { RandomStream } from './ZZ.js';
import { zz_pXModulus } from './lzz_pX.js';
export function ntl_cold_recombination(
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
  const shared = { context: cache, stream };
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
        new zz_pXModulus(null, p, { maxroot, state: shared });
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
        result = UpdateLocalInfo(I, W, factors, f, count, modulus, { ...options, state: shared });
      } else if (op === 16 || op === 17 || op === 18) {
        const modulus = scalar(),
          count = word(),
          bound = word(),
          van_Hoeij = word(),
          MaxNumPrimes = word(),
          MaxPrune = word();
        const factors = Array.from({ length: word() }, () => vector()),
          f = vector(),
          W = Array.from({ length: word() }, () => vector());
        const options = { van_Hoeij, MaxNumPrimes, MaxPrune, state: shared };
        result =
          op === 16
            ? CardinalitySearch(factors, f, W, I, count, bound, modulus, options)
            : op === 17
              ? CardinalitySearch1(factors, f, W, I, count, bound, modulus, options)
              : FindTrueFactors(f, W, modulus, I, bound, options);
      } else if (op === 20) UseFFTPrime(word(), cache, stream);
      else if (op === 21) stream.set_nonce(scalar());
      else if (op === 22) result = hex(stream.get(word()));
      else if (op === 23) {
        const p = scalar(),
          I = vector().map(Number),
          W = Array.from({ length: word() }, () => vector());
        result = InvMul(W, I, p, shared);
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

const parse = (s: string): bigint[] =>
  s.slice(1, -1).trim()
    ? s
        .slice(1, -1)
        .split(',')
        .map((x) => BigInt(x.trim()))
    : [];
for (const [i, row] of fixtures.entries())
  test('native NTL cold local updates and recombination ' +
    i +
    ([3798, 3850, 3851, 3852, 3853].includes(i) ? ' — explicit native indexing guard' : ''), () => {
    let result: string | null = null,
      error: string | null = null,
      errorType: string | null = null;
    try {
      result = ntl_cold_recombination(
        parse(row.args[0]!),
        parse(row.args[1]!),
        BigInt(row.args[2]!),
        BigInt(row.args[3]!)
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
test('NTL stateful recombination owns returned factors and complementary products', () => {
  const state = { context: new FFTPrimeContext(), stream: new RandomStream(new Uint8Array(32)) },
    info = new LocalInfoT();
  info.n = 20;
  info.PossibleDegrees = (1n << 21n) - 1n;
  const g = [1n, ...Array<bigint>(9).fill(0n), 1n],
    h = [2n, ...Array<bigint>(9).fill(0n), 1n],
    f = [2n, ...Array<bigint>(9).fill(0n), 3n, ...Array<bigint>(9).fill(0n), 1n],
    W = [g, h],
    previous = [[7n]],
    snapshot = structuredClone([W, f, previous]);
  const [factors, remaining, lifted] = CardinalitySearch(previous, f, W, info, 1, 100, 65537n, {
    state,
  });
  expect(factors).toEqual([[7n], g]);
  expect(state.context.length()).toBe(2);
  factors[0]![0] = 99n;
  factors[1]![0] = 99n;
  remaining[0] = 99n;
  lifted[0]![0] = 99n;
  const complement = InvMul(W, [0], 65537n, state);
  complement[0] = 99n;
  expect([W, f, previous]).toEqual(snapshot);
});
