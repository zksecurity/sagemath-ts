import { nativeFixtures as loadLiveNative } from "../../../tests/property/native-live.mjs";


import { test, expect } from 'bun:test';
const fixtures = (await loadLiveNative(import.meta.url, "./ZZ_prime_generation.native.json.gz")) as { args: (string)[]; error: null; errorType: null; function: string; result: string; seed: number }[];
import {
  RandomStream,
  ComputePrimeBound,
  ErrBoundTest,
  ProbPrime,
  RandomPrime,
  OldRandomPrime,
  RandomPrime_long,
  GenPrime,
  GenPrime_long,
} from './ZZ.js';
export function ntl_prime_generation(key: bigint[], cmd: bigint[]): string {
  const bytes = Uint8Array.from(key, Number),
    s = [new RandomStream(bytes), new RandomStream(bytes), new RandomStream(bytes)],
    out: unknown[] = [];
  const hex = (a: Uint8Array) => Array.from(a, (x) => x.toString(16).padStart(2, '0')).join('');
  for (let i = 0; i < cmd.length; i += 5) {
    const op = Number(cmd[i]),
      slot = Number(cmd[i + 1]),
      a = cmd[i + 2]!,
      b = Number(cmd[i + 3]),
      c = Number(cmd[i + 4]);
    let result: unknown = null,
      error: string | null = null,
      errorType: string | null = null;
    try {
      if (op === 0) result = ProbPrime(a, s[slot]!, { NumTrials: b, word: !!c });
      else if (op === 1) result = ProbPrime(a, s[slot]!);
      else if (op === 2) result = ComputePrimeBound(Number(a));
      else if (op === 3) result = ErrBoundTest(Number(a), b, c);
      else if (op === 4)
        result = RandomPrime(Number(a), s[slot]!, c ? undefined : { NumTrials: b });
      else if (op === 5)
        result = OldRandomPrime(Number(a), s[slot]!, c ? undefined : { NumTrials: b });
      else if (op === 6)
        result = RandomPrime_long(Number(a), s[slot]!, c ? undefined : { NumTrials: b });
      else if (op === 7) result = GenPrime(Number(a), s[slot]!, c ? undefined : { err: b });
      else if (op === 8) result = GenPrime_long(Number(a), s[slot]!, c ? undefined : { err: b });
      else if (op === 9) result = hex(s[slot]!.get(Number(a)));
      else if (op === 10) s[slot]!.set_nonce(a);
      else if (op === 11) s[slot] = new RandomStream(s[Number(a)]!);
      else if (op === 12) s[slot]!.assign(s[Number(a)]!);
      else if (op === 13) s[slot] = s[Number(a)]!;
      else throw new Error('unknown prime-generation operation');
    } catch (e) {
      error = (e as Error).message;
      errorType = (e as Error).name;
    }
    out.push([errorType, error, result]);
  }
  for (const stream of s) out.push([null, null, hex(stream.get(64))]);
  return JSON.stringify(out, (_, v) => (typeof v === 'bigint' ? String(v) : v));
}

const parse = (s: string): bigint[] =>
  s
    .slice(1, -1)
    .split(',')
    .filter((x) => x.trim())
    .map((x) => BigInt(x.trim()));
for (const [i, row] of fixtures.entries())
  test('native NTL prime generation ' + i, () => {
    let result: string | null = null,
      error: string | null = null,
      errorType: string | null = null;
    try {
      result = ntl_prime_generation(parse(row.args[0]!), parse(row.args[1]!));
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

test('NTL existing class-wrapper overloads retain their unimplemented errors', () => {
  // ZZ itself is still unconstructible; the old wrapper unconditionally delegates.
  expect(() => ProbPrime(Object.create(null), 2)).toThrow('NTL_NOT_IMPLEMENTED: ZZ.ProbPrime');
  expect(() => RandomPrime(32, 2)).toThrow('NTL_NOT_IMPLEMENTED: ZZ.RandomPrime');
  expect(() => RandomPrime(32)).toThrow('NTL_NOT_IMPLEMENTED: ZZ.RandomPrime');
});
