import { GotThem } from '../../../packages/ntl-ts/src/ZZXFactoring.js';
import { RandomStream } from '../../../packages/ntl-ts/src/ZZ.js';
import {
  FFTPrimeContext,
  UseFFTPrime,
  type FFTPrimeInfo,
} from '../../../packages/ntl-ts/src/FFT.js';
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
