import { mul } from '../../../packages/ntl-ts/src/ZZ_pX.js';
import { mul as treeMul } from '../../../packages/ntl-ts/src/ZZXFactoring.js';
import { SSRatio } from '../../../packages/ntl-ts/src/ZZX1.js';
import { RandomStream } from '../../../packages/ntl-ts/src/ZZ.js';
import {
  FFTPrimeContext,
  UseFFTPrime,
  type FFTPrimeInfo,
} from '../../../packages/ntl-ts/src/FFT.js';
export function ntl_polynomial_state(
  key: bigint[],
  p: bigint,
  lengths: bigint[],
  coefficients: bigint[],
  indices: bigint[],
  commands: bigint[]
): string {
  let offset = 0;
  const W = lengths.map((n) => {
    const w = coefficients.slice(offset, offset + Number(n));
    offset += Number(n);
    return w;
  });
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
      errorType: string | null = null;
    try {
      if (op === 0) result = mul(W[0]!, W[1]!, p, { context, stream });
      else if (op === 1) result = mul(W[0]!, W[0]!, p, { context, stream });
      else if (op === 2) result = treeMul(W, p, undefined, { context, stream });
      else if (op === 3) result = treeMul(W, p, indices.map(Number), { context, stream });
      else if (op === 4) UseFFTPrime(Number(a), context, stream);
      else if (op === 5) {
        if (a <= 1n) throw new Error('ZZ_pContext: p must be > 1');
        p = a;
      } else if (op === 6) result = hex(stream.get(Number(a)));
      else if (op === 7) stream.set_nonce(a);
      else if (op === 8) {
        const v = W[Number(a)]!;
        result = bits(
          SSRatio(Number(v[0] ?? 0n), Number(v[1] ?? 0n), Number(v[2] ?? 0n), Number(v[3] ?? 0n))
        );
      } else throw new Error('unknown polynomial state operation');
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
