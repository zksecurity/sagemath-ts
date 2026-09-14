import { nativeFixtures as loadLiveNative } from "../../../tests/property/native-live.mjs";


import { test, expect } from 'bun:test';
const fixtures = (await loadLiveNative(import.meta.url, "./ZZ_random_stream.native.json.gz")) as { args: (string)[]; error: null | string; errorType: null | string; function: string; result: null | string; seed: number }[];
import {
  sha256,
  hmac_sha256,
  DeriveKey,
  RandomStream,
  salsa20_core,
  salsa20_apply,
  salsa20_init,
} from './ZZ.js';
export function ntl_random_stream(
  op: bigint,
  n: bigint,
  av: bigint[],
  bv: bigint[],
  commands: bigint[]
): string {
  const a = Uint8Array.from(av, Number),
    b = Uint8Array.from(bv, Number);
  const hex = (data: Uint8Array) =>
    Array.from(data, (x) => x.toString(16).padStart(2, '0')).join('');
  const words = (state: Uint32Array) => Array.from(state, String);
  let out: unknown;
  if (op === 0n) out = hex(sha256(a, Number(n)));
  else if (op === 1n) out = hex(hmac_sha256(b, a, Number(n)));
  else if (op === 2n) out = hex(DeriveKey(a, Number(n)));
  else if (op === 7n) {
    // A view with a reported native pointer length; only its low 32 bits are read.
    Object.defineProperty(a, 'length', { value: Number(commands[0]) });
    out = hex(hmac_sha256(b, a, Number(n)));
  } else if (op === 3n) {
    const streams = [new RandomStream(a), new RandomStream(a), new RandomStream(a)],
      trace: unknown[] = [];
    for (let i = 0; i < commands.length; i += 3) {
      const action = commands[i],
        slot = Number(commands[i + 1]),
        arg = commands[i + 2]!;
      try {
        if (action === 0n) trace.push([null, hex(streams[slot]!.get(Number(arg)))]);
        else {
          if (action === 1n) streams[slot]!.set_nonce(arg);
          else if (action === 2n) streams[slot] = new RandomStream(streams[Number(arg)]!);
          else if (action === 3n) streams[slot] = streams[Number(arg)]!;
          else if (action === 4n) streams[slot]!.assign(streams[Number(arg)]!);
          else throw new Error('unknown stream command');
          trace.push([null, null]);
        }
      } catch (e) {
        trace.push([(e as Error).message, null]);
      }
    }
    out = trace;
  } else if (op === 4n || op === 5n) {
    const state = Uint32Array.from(av, Number);
    if (op === 4n) {
      salsa20_core(state);
      out = words(state);
    } else {
      const block = salsa20_apply(state);
      out = [words(state), words(block)];
    }
  } else if (op === 6n) out = words(salsa20_init(a));
  else throw new Error('unknown NTL random stream operation');
  return JSON.stringify(out);
}

const parse = (s: string): bigint[] =>
  s
    .slice(1, -1)
    .split(',')
    .filter((x) => x.trim())
    .map((x) => BigInt(x.trim()));
for (const [i, row] of fixtures.entries())
  test('native NTL byte stream ' + i, () => {
    let result: string | null = null,
      error: string | null = null,
      errorType: string | null = null;
    try {
      result = ntl_random_stream(
        BigInt(row.args[0]!),
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
test('NTL byte stream key view offset and returned buffer ownership', () => {
  const raw = Uint8Array.from({ length: 80 }, (_, i) => i),
    key = raw.subarray(11, 43),
    original = key.slice();
  const stream = new RandomStream(key),
    control = new RandomStream(original);
  raw.fill(255);
  const first = stream.get(17);
  expect(first).toEqual(control.get(17));
  first.fill(0);
  expect(stream.get(97)).toEqual(control.get(97));
});
test('NTL byte stream assignment returns its receiver and preserves aliases', () => {
  const key = new Uint8Array(32),
    a = new RandomStream(key),
    b = new RandomStream(key),
    alias = a;
  b.get(17);
  expect(a.assign(b)).toBe(a);
  expect(alias.get(63)).toEqual(b.get(63));
  expect(a.assign(a)).toBe(a);
  expect(a.get(65)).toEqual(b.get(65));
});
test('NTL byte helpers retain typed-array views and own their results', () => {
  const data = Uint8Array.from({ length: 100 }, (_, i) => i),
    view = data.subarray(13, 78),
    copy = view.slice();
  expect(sha256(view)).toEqual(sha256(copy));
  expect(hmac_sha256(view, view)).toEqual(hmac_sha256(copy, copy));
  expect(DeriveKey(view, 65)).toEqual(DeriveKey(copy, 65));
  const state = salsa20_init(view),
    before = state.slice(),
    block = salsa20_apply(state);
  block.fill(0);
  expect(Array.from(state.slice(0, 12))).toEqual(Array.from(before.slice(0, 12)));
  expect(state[12]).toBe(1);
});
