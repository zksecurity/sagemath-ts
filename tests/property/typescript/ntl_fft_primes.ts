import { RandomStream } from '../../../packages/ntl-ts/src/ZZ.js';
import {
  FFTPrimeContext,
  IsFFTPrime,
  NextFFTPrime,
  CalcMaxRoot,
  InitFFTPrimeInfo,
  UseFFTPrime,
  GetFFTPrime,
  GetFFTPrimeRecip,
  type FFTPrimeInfo,
} from '../../../packages/ntl-ts/src/FFT.js';
export function ntl_fft_primes(key: bigint[], commands: bigint[]): string {
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
  for (let i = 0; i < commands.length; i += 4) {
    const op = Number(commands[i]),
      a = commands[i + 1]!,
      b = commands[i + 2]!;
    let result: unknown = null,
      error: string | null = null,
      errorType: string | null = null;
    try {
      if (op === 0) result = IsFFTPrime(a, stream, b);
      else if (op === 11) result = IsFFTPrime(a, stream);
      else if (op === 1) result = CalcMaxRoot(a);
      else if (op === 2) result = NextFFTPrime(Number(a), context, stream);
      else if (op === 3) {
        UseFFTPrime(Number(a), context, stream);
        result = null;
      } else if (op === 4) result = GetFFTPrime(Number(a), context);
      else if (op === 5) result = bits(GetFFTPrimeRecip(Number(a), context));
      else if (op === 6) result = info(context.get(Number(a)));
      else if (op === 7) result = info(InitFFTPrimeInfo(a, b));
      else if (op === 8) result = context.length();
      else if (op === 9) result = hex(stream.get(Number(a)));
      else if (op === 10) stream.set_nonce(a);
      else throw new Error('unknown FFT-prime operation');
    } catch (e) {
      error = (e as Error).message;
      errorType = (e as Error).name;
    }
    out.push([errorType, error, result, context.length(), hex(new RandomStream(stream).get(64))]);
  }
  return JSON.stringify(out, (_, v) => (typeof v === 'bigint' ? String(v) : v));
}
