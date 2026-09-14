import { determinant, inv } from '../../../packages/ntl-ts/src/mat_ZZ.js';
import { RandomStream, GenPrime } from '../../../packages/ntl-ts/src/ZZ.js';
import {
  FFTPrimeContext,
  NextFFTPrime,
  UseFFTPrime,
  type FFTPrimeInfo,
} from '../../../packages/ntl-ts/src/FFT.js';
export function ntl_integer_reconstruction(
  key: bigint[],
  dims: bigint[],
  av: bigint[],
  prev: bigint[],
  commands: bigint[]
): string {
  const [n, m, r, c, badA, badPrev] = dims.map(Number);
  const make = (v: bigint[], rows: number, cols: number, bad: number): bigint[][] => {
    if (bad) throw new Error('nonrectangular matrix');
    if (rows < 0 || cols < 0) throw new Error('SetDims: bad args');
    return Array.from({ length: rows }, (_, i) => v.slice(i * cols, (i + 1) * cols));
  };
  const A = make(av, n!, m!, badA!),
    previous = make(prev, r!, c!, badPrev!);
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
      if (op === 0) result = determinant(A, context, stream, { columns: m });
      else if (op === 1)
        result = determinant(A, context, stream, { columns: m, deterministic: !!a });
      else if (op === 2)
        result = inv(A, context, stream, {
          columns: m,
          status: true,
          deterministic: !!a,
          previous,
        });
      else if (op === 3) result = inv(A, context, stream, { columns: m });
      else if (op === 4) UseFFTPrime(Number(a), context, stream);
      else if (op === 5) result = NextFFTPrime(Number(a), context, stream);
      else if (op === 6) result = hex(stream.get(Number(a)));
      else if (op === 7) stream.set_nonce(a);
      else if (op === 8) {
        if (a < 0n) throw new Error('bad FFT prime index');
        UseFFTPrime(Number(a), context, stream);
      } else if (op === 9) {
        if (a <= 1n) throw new Error('ZZ_pContext: p must be > 1');
      } else if (op === 10) result = GenPrime(Number(a), stream, { err: 90 });
      else throw new Error('unknown integer reconstruction operation');
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
