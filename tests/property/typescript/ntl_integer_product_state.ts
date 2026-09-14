import { mul, sqr, ChooseSS } from '../../../packages/ntl-ts/src/ZZX1.js';
import { RandomStream } from '../../../packages/ntl-ts/src/ZZ.js';
import {
  FFTPrimeContext,
  UseFFTPrime,
  type FFTPrimeInfo,
} from '../../../packages/ntl-ts/src/FFT.js';
export function ntl_integer_product_state(
  key: bigint[],
  a: bigint[],
  b: bigint[],
  commands: bigint[]
): string {
  const stream = new RandomStream(Uint8Array.from(key, Number)),
    context = new FFTPrimeContext(),
    state = { context, stream },
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
      p = commands[i + 1]!;
    let result: unknown = null,
      errorType: string | null = null,
      error: string | null = null;
    try {
      if (op === 0) result = mul(a, b, state);
      else if (op === 1) result = mul(a, a, state);
      else if (op === 2) result = sqr(a, state);
      else if (op === 3) result = mul(b, a, state);
      else if (op === 4 || op === 7) UseFFTPrime(Number(p), context, stream);
      else if (op === 5) stream.set_nonce(p);
      else if (op === 6) result = hex(stream.get(Number(p)));
      else if (op === 8) {
        if (p <= 1n) throw new Error('ZZ_pContext: p must be > 1');
      } else if (op === 9)
        result = ChooseSS(
          Number(a[0] ?? 0n),
          Number(a[1] ?? 0n),
          Number(a[2] ?? 0n),
          Number(a[3] ?? 0n)
        );
      else throw new Error('unknown integer polynomial state operation');
    } catch (e) {
      errorType = (e as Error).name;
      error = (e as Error).message;
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
