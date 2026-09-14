import { SquareFreeDecomp } from '../../../packages/ntl-ts/src/ZZXFactoring.js';
import { GCD, mul } from '../../../packages/ntl-ts/src/ZZX1.js';
import { zz_pXModulus } from '../../../packages/ntl-ts/src/lzz_pX.js';
import { FFTPrimeContext, type FFTPrimeInfo } from '../../../packages/ntl-ts/src/FFT.js';
import { RandomStream } from '../../../packages/ntl-ts/src/ZZ.js';
export function ntl_integer_gcd_state(
  key: bigint[],
  p: bigint,
  maxroot: bigint,
  lengths: bigint[],
  coefficients: bigint[],
  commands: bigint[]
): string {
  const stream = new RandomStream(Uint8Array.from(key, Number)),
    context = new FFTPrimeContext(),
    state = { context, stream };
  if (p !== 0n) new zz_pXModulus(null, p, { maxroot: Number(maxroot), state });
  let offset = 0;
  const polynomials = lengths.map((n) => {
    const a = coefficients.slice(offset, offset + Number(n));
    offset += Number(n);
    return a;
  });
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
  const out: unknown[] = [];
  for (let i = 0; i < commands.length; i += 2) {
    const op = Number(commands[i]),
      a = polynomials[Number(commands[i + 1])]!;
    let result: unknown = null,
      errorType: string | null = null,
      error: string | null = null;
    try {
      if (op === 0) result = GCD(polynomials[Number(a[0])]!, polynomials[Number(a[1])]!, state);
      else if (op === 1) result = SquareFreeDecomp(polynomials[Number(a[0])]!, state);
      else if (op === 2)
        result = mul(polynomials[Number(a[0])]!, polynomials[Number(a[1])]!, state);
      else throw new Error('unknown integer GCD state operation');
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
