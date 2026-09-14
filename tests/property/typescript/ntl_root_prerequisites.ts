import {
  RandomStream,
  RandomBndGenerator,
  VectorRandomBnd,
} from '../../../packages/ntl-ts/src/ZZ.js';
import { random as randomPolynomial, zz_pXModulus } from '../../../packages/ntl-ts/src/lzz_pX.js';
import { BuildFromRoots } from '../../../packages/ntl-ts/src/lzz_pX1.js';
export function ntl_root_prerequisites(
  op: bigint,
  p: bigint,
  n: bigint,
  data: bigint[],
  key: bigint[]
): string {
  if (op <= 1n) new zz_pXModulus(null, p);
  const stream = new RandomStream(Uint8Array.from(key, Number)),
    hex = (b: Uint8Array) => Array.from(b, (x) => x.toString(16).padStart(2, '0')).join('');
  const secondary = new RandomStream(Uint8Array.from(key, Number));
  secondary.set_nonce(1n);
  let errorType: string | null = null,
    error: string | null = null,
    value: unknown = null;
  try {
    if (op === 0n) value = BuildFromRoots(data, p);
    else if (op === 1n) value = randomPolynomial(Number(n), p, stream);
    else if (op === 2n) value = VectorRandomBnd(Number(n), p, stream);
    else if (op === 3n) {
      const g = [
          new RandomBndGenerator(n ? p : null, stream),
          new RandomBndGenerator(n ? p : null, n === 2n ? secondary : stream),
        ],
        trace: unknown[] = [];
      const state = () => g.map((x) => [x.p, x.p ? x.nb : null, x.p ? x.mask : null]);
      for (let i = 0; i < data.length; i += 3) {
        const action = data[i],
          slot = Number(data[i + 1]),
          arg = data[i + 2]!;
        let result: unknown = null,
          kind: string | null = null,
          message: string | null = null;
        try {
          if (action === 0n) g[slot]!.build(arg);
          else if (action === 1n) {
            const r: bigint[] = [];
            for (let j = 0; j < Number(arg); j++) r.push(g[slot]!.next());
            result = r;
          } else if (action === 2n) g[slot] = new RandomBndGenerator(g[Number(arg)]!);
          else if (action === 3n) g[slot] = g[Number(arg)]!;
          else if (action === 4n)
            g[slot] = new RandomBndGenerator(null, n === 2n && slot === 1 ? secondary : stream);
          else if (action === 5n) result = hex(stream.get(Number(arg)));
          else if (action === 6n) stream.set_nonce(arg);
          else if (action === 7n) g[slot]!.assign(g[Number(arg)]!);
          else throw new Error('unknown cached sampler command');
        } catch (e) {
          kind = (e as Error).name;
          message = (e as Error).message;
        }
        trace.push([kind, message, result, state()]);
      }
      value = trace;
    } else throw new Error('unknown root prerequisite operation');
  } catch (e) {
    errorType = (e as Error).name;
    error = (e as Error).message;
  }
  return JSON.stringify(
    [
      errorType,
      error,
      value,
      op === 3n && n === 2n ? [hex(stream.get(64)), hex(secondary.get(64))] : hex(stream.get(64)),
    ],
    (_, v) => (typeof v === 'bigint' ? String(v) : v)
  );
}
