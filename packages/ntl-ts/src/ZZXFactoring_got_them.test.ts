import { nativeFixtures as loadLiveNative } from "../../../tests/property/native-live.mjs";
import { test, expect } from 'bun:test';


const fixtures = (await loadLiveNative(import.meta.url, "./ZZXFactoring_got_them.native.json.gz")) as { regressionId?: string; args: string[]; result: string | null; error: string | null; errorType: string | null }[];
import { GotThem } from './ZZXFactoring.js';
import { RandomStream } from './ZZ.js';
import { FFTPrimeContext, UseFFTPrime, type FFTPrimeInfo } from './FFT.js';
export function ntl_got_them(
  key: bigint[],
  dims: bigint[],
  flat: bigint[],
  lengths: bigint[],
  coefficients: bigint[],
  f: bigint[],
  previousLengths: bigint[],
  previousCoefficients: bigint[],
  bound: bigint,
  p: bigint,
  commands: bigint[]
): string {
  const [n, m] = dims.map(Number);
  if (n! < 0 || m! < 0) throw new Error('SetDims: bad args');
  const B = Array.from({ length: n! }, (_, i) => flat.slice(i * m!, (i + 1) * m!));
  const polys = (lengths: bigint[], coefficients: bigint[]) => {
    let offset = 0;
    return lengths.map((v) => {
      const out = coefficients.slice(offset, offset + Number(v));
      offset += Number(v);
      return out;
    });
  };
  const W = polys(lengths, coefficients),
    previous = polys(previousLengths, previousCoefficients);
  if (p <= 1n) throw new Error('ZZ_pContext: p must be > 1');
  const stream = new RandomStream(Uint8Array.from(key, Number)),
    context = new FFTPrimeContext(),
    out: unknown[] = [];
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
      a = commands[i + 1]!;
    let result: unknown = null,
      error: string | null = null,
      errorType: string | null = null,
      diagnostics = '';
    const oldError = console.error;
    console.error = (...args: unknown[]) => {
      diagnostics += args.map(String).join(' ') + '\n';
    };
    try {
      if (op === 0) {
        const [status, factors] = GotThem(previous, B, W, f, Number(bound), p, context, stream, {
          columns: m,
        });
        result = [BigInt(status), factors];
      } else if (op === 1) UseFFTPrime(Number(a), context, stream);
      else if (op === 2) stream.set_nonce(a);
      else if (op === 3) result = hex(stream.get(Number(a)));
      else if (op === 4) {
        if (a < 0n) throw new Error('bad FFT prime index');
        UseFFTPrime(Number(a), context, stream);
      } else throw new Error('unknown GotThem operation');
    } catch (e) {
      error = (e as Error).message;
      errorType = (e as Error).name;
    } finally {
      console.error = oldError;
    }
    out.push([
      errorType,
      error,
      result,
      Array.from({ length: context.length() }, (_, j) => info(context.get(j))),
      hex(new RandomStream(stream).get(64)),
      hex(new TextEncoder().encode(diagnostics)),
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
  test('native NTL GotThem acceptance ' +
    i +
    ([1661, 1663].includes(i) ? ' — explicit bound adapter guard' : ''), () => {
    let result: string | null = null,
      error: string | null = null,
      errorType: string | null = null;
    try {
      result = ntl_got_them(
        parse(row.args[0]!),
        parse(row.args[1]!),
        parse(row.args[2]!),
        parse(row.args[3]!),
        parse(row.args[4]!),
        parse(row.args[5]!),
        parse(row.args[6]!),
        parse(row.args[7]!),
        BigInt(row.args[8]!),
        BigInt(row.args[9]!),
        parse(row.args[10]!)
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
test('NTL GotThem owns appended and retained factors and preserves its input arrays', () => {
  const previous = [[7n]],
    B = [
      [1n, 1n, 1n, 1n, 0n, 0n, 0n, 0n],
      [0n, 0n, 0n, 0n, 1n, 1n, 1n, 1n],
    ],
    W = [
      ...Array.from({ length: 4 }, () => [1n, 1n]),
      ...Array.from({ length: 4 }, () => [-1n, 1n]),
    ],
    f = [1n, 0n, -4n, 0n, 6n, 0n, -4n, 0n, 1n],
    snapshot = structuredClone([previous, B, W, f]);
  const context = new FFTPrimeContext(),
    stream = new RandomStream(new Uint8Array(32));
  const [status, factors] = GotThem(previous, B, W, f, 16, 65537n, context, stream);
  expect(status).toBe(1);
  expect(factors[0]).toEqual([7n]);
  factors[0]![0] = 99n;
  factors[1]![0] = 99n;
  const [rejected, retained] = GotThem(previous, B, W, f, -1, 65537n, context, stream);
  expect([rejected, retained]).toEqual([0, [[7n]]]);
  retained[0]![0] = 99n;
  expect([previous, B, W, f]).toEqual(snapshot);
});
test('NTL GotThem rejects non-native numeric coefficient bounds', () => {
  for (const bound of [NaN, Infinity, -Infinity, 1.5])
    expect(() =>
      GotThem(
        [],
        [[1n]],
        [[1n, 1n]],
        [1n],
        bound,
        17n,
        new FFTPrimeContext(),
        new RandomStream(new Uint8Array(32))
      )
    ).toThrow('GotThem bound must be a signed native integer');
});

test('NTL GotThem rejects invalid moduli with the native context error', () => {
  for (const i of [1657, 1658, 1659]) {
    const row = fixtures.find(row => row.regressionId === String(i))!;
    expect(() =>
      GotThem(
        [],
        [[1n, 1n, 1n, 1n]],
        [
          [1n, 1n],
          [1n, 1n],
          [1n, 1n],
          [1n, 1n],
        ],
        [1n],
        20,
        BigInt(row.args[9]!),
        new FFTPrimeContext(),
        new RandomStream(new Uint8Array(32))
      )
    ).toThrow(row.error!);
  }
});
