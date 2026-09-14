import { nativeFixtures as loadLiveNative } from "../../../tests/property/native-live.mjs";


import { test, expect } from 'bun:test';
const fixtures = (await loadLiveNative(import.meta.url, "./ZZ_prime_sequence.native.json.gz")) as { args: (string)[]; error: null; errorType: null; function: string; result: string; seed: number }[];
import { PrimeSeq } from './ZZ.js';
export function ntl_prime_sequence(commands: bigint[]): string {
  const seq = [new PrimeSeq(), new PrimeSeq()],
    out: unknown[] = [];
  for (let i = 0; i < commands.length; i += 3) {
    const op = commands[i],
      slot = Number(commands[i + 1]),
      arg = commands[i + 2]!;
    if (op === 0n) out.push(Array.from({ length: Number(arg) }, () => seq[slot]!.next()));
    else if (op === 1n) {
      seq[slot]!.reset(arg);
      out.push(null);
    } else if (op === 2n) {
      seq[slot] = new PrimeSeq();
      out.push(null);
    } else throw new Error('unknown NTL prime sequence operation');
  }
  return JSON.stringify(out, (_, v) => (typeof v === 'bigint' ? String(v) : v));
}

for (const [i, row] of fixtures.entries())
  test('native NTL prime sequence ' + i, () => {
    const commands = row.args[0]!.slice(1, -1)
      .split(',')
      .filter((x) => x.trim())
      .map((x) => BigInt(x.trim()));
    let result: string | null = null,
      error: string | null = null,
      errorType: string | null = null;
    try {
      result = ntl_prime_sequence(commands);
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
