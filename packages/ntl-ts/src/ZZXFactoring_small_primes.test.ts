import { nativeFixtures as loadLiveNative } from "../../../tests/property/native-live.mjs";


import { test, expect } from 'bun:test';
const fixtures = (await loadLiveNative(import.meta.url, "./ZZXFactoring_small_primes.native.json.gz")) as { args: (string)[]; error: null; errorType: null; function: string; result: string; seed: number }[];
import { LocalInfoT, SmallPrimeFactorization } from './ZZXFactoring.js';
import { RandomStream } from './ZZ.js';
import { zz_pXModulus } from './lzz_pX.js';
export function ntl_small_prime_selection(
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

const parse = (s: string): bigint[] =>
  s
    .slice(1, -1)
    .split(',')
    .filter((x) => x.trim())
    .map((x) => BigInt(x.trim()));
for (const [i, row] of fixtures.entries())
  test('native NTL stateful small-prime selection ' + i, () => {
    let result: string | null = null,
      error: string | null = null,
      errorType: string | null = null;
    try {
      result = ntl_small_prime_selection(
        parse(row.args[0]!),
        parse(row.args[1]!),
        BigInt(row.args[2]!)
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

test('NTL local information vector assignments and snapshots own their storage', () => {
  const I = new LocalInfoT(),
    p = [101n, 103n],
    pattern = [[2, 3], [4]];
  I.p = p;
  I.pattern = pattern;
  p[0] = 0n;
  pattern[0]![0] = 0;
  const primes = I.p,
    patterns = I.pattern;
  primes[1] = 0n;
  patterns[1]![0] = 0;
  expect(I.p).toEqual([101n, 103n]);
  expect(I.pattern).toEqual([[2, 3], [4]]);
});
test('NTL local information instances have independent sequence and vector state', () => {
  const a = new LocalInfoT(),
    b = new LocalInfoT();
  a.p = [17n];
  a.pattern = [[9]];
  expect([a.s.next(), a.s.next(), b.s.next()]).toEqual([2n, 3n, 2n]);
  expect(b.p).toEqual([]);
  expect(b.pattern).toEqual([]);
});
