import { nativeFixtures as loadLiveNative } from "../../../tests/property/native-live.mjs";
import { test, expect } from 'bun:test';


const fixtures = (await loadLiveNative(import.meta.url, "./ZZX1_division.native.json.gz")) as { regressionId?: string; args: string[]; result: string | null; error: string | null; errorType: string | null }[];
import * as polynomials from './ZZX1.js';
import { RandomStream } from './ZZ.js';
import { FFTPrimeContext, UseFFTPrime, type FFTPrimeInfo } from './FFT.js';
export function ntl_integer_polynomial_division(
  key: bigint[],
  A: bigint[],
  B: bigint[],
  previous: bigint[],
  commands: bigint[]
): string {
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
      errorType: string | null = null;
    try {
      if (op <= 3 || op >= 9) {
        const options = { previous, state: { context, stream } };
        if (op === 12) result = polynomials._ZZX_kernels.divide(A, B, options.state);
        else {
          const [status, q] =
            op === 1 || op === 10
              ? polynomials.PlainDivide(A, B, options)
              : op === 2 || op === 11
                ? polynomials.HomDivide(A, B, options)
                : polynomials.divide(A, op === 3 ? (B[0] ?? 0n) : B, options);
          result = op >= 9 ? BigInt(status) : [BigInt(status), q];
        }
      } else if (op === 4) UseFFTPrime(Number(a), context, stream);
      else if (op === 5) {
        if (a < 0n) throw new Error('bad FFT prime index');
        UseFFTPrime(Number(a), context, stream);
      } else if (op === 6) result = hex(stream.get(Number(a)));
      else if (op === 7) stream.set_nonce(a);
      else if (op === 8) {
        if (a <= 1n) throw new Error('ZZ_pContext: p must be > 1');
      } else throw new Error('unknown integer polynomial division operation');
    } catch (e) {
      error = (e as Error).message;
      errorType = (e as Error).name;
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

const parse = (s: string): bigint[] =>
  s.slice(1, -1).trim()
    ? s
        .slice(1, -1)
        .split(',')
        .map((x) => BigInt(x.trim()))
    : [];
for (const [i, row] of fixtures.entries())
  test('native NTL integer polynomial exact division ' + i, () => {
    let result: string | null = null,
      error: string | null = null,
      errorType: string | null = null;
    try {
      result = ntl_integer_polynomial_division(
        parse(row.args[0]!),
        parse(row.args[1]!),
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
  }, 30000);
test('NTL exact division owns successful and retained outputs and keeps argument aliases intact', () => {
  const A = [1n, 2n, 1n],
    B = [1n, 1n],
    previous = [7n, -3n],
    snapshot = structuredClone([A, B, previous]);
  for (const fn of [polynomials.divide, polynomials.PlainDivide, polynomials.HomDivide]) {
    const context = new FFTPrimeContext(),
      stream = new RandomStream(new Uint8Array(32)),
      options = { previous, state: { context, stream } };
    const [status, q] = fn(A, B, options);
    expect([status, q]).toEqual([1, [1n, 1n]]);
    q[0] = 99n;
    const [failed, retained] = fn([1n, 1n], [1n, 0n, 1n], options);
    expect([failed, retained]).toEqual([0, previous]);
    retained[0] = 99n;
    expect(fn(A, A, options)).toEqual([1, [1n]]);
    expect([A, B, previous]).toEqual(snapshot);
  }
});

// The coefficient-only API intentionally has no externally observable stream.
// Compare its complete status/error/quotient with the same native division rows.
for (const index of [1797, 1896, 1929, 1962, 2094]) {
  test('native coefficient-only homomorphic division ' + index, () => {
    const row = fixtures.find(row => row.regressionId === String(index))!,
      commands = parse(row.args[4]!);
    const operation = commands.findIndex((x, i) => i % 2 === 0 && x === 2n) / 2;
    const expected = JSON.parse(row.result!)[operation].slice(0, 3);
    let result: unknown = null,
      error: string | null = null,
      errorType: string | null = null;
    try {
      const [status, q] = polynomials.HomDivide(parse(row.args[1]!), parse(row.args[2]!), {
        previous: parse(row.args[3]!),
      });
      result = [String(status), q.map(String)];
    } catch (e) {
      error = (e as Error).message;
      errorType = (e as Error).name;
    }
    expect([errorType, error, result]).toEqual(expected);
  }, 30000);
}
