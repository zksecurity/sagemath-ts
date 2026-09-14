import { zz_pXModulus } from '../../../packages/ntl-ts/src/lzz_pX.js';
import {
  DoMinPolyMod,
  IrredPolyMod,
  ProbMinPolyMod,
  MinPolyMod,
} from '../../../packages/ntl-ts/src/lzz_pX1.js';
import { RandomStream } from '../../../packages/ntl-ts/src/ZZ.js';
export function ntl_element_minpoly(
  op: bigint,
  p: bigint,
  m: bigint,
  ready: bigint,
  f: bigint[],
  g: bigint[],
  r: bigint[],
  key: bigint[]
): string {
  const F = new zz_pXModulus(ready ? f : null, p),
    G = F.arithmetic.norm(g),
    R = r.map((x) => F.arithmetic.mod(x)),
    stream = new RandomStream(Uint8Array.from(key, Number));
  let errorType: string | null = null,
    error: string | null = null,
    result: bigint[] | null = null;
  try {
    if (op === 0n) result = DoMinPolyMod(G, F, Number(m), R);
    else if (op === 1n) result = IrredPolyMod(G, F, Number(m));
    else if (op === 2n) result = ProbMinPolyMod(G, F, Number(m), stream);
    else if (op === 3n) result = MinPolyMod(G, F, Number(m), stream);
    else if (op === 4n) result = IrredPolyMod(G, F);
    else if (op === 5n) result = ProbMinPolyMod(G, F, stream);
    else if (op === 6n) result = MinPolyMod(G, F, stream);
    else throw new Error('unknown element minimum polynomial operation');
  } catch (e) {
    error = (e as Error).message;
    errorType = (e as Error).name;
  }
  const tail = Array.from(stream.get(64), (x) => x.toString(16).padStart(2, '0')).join('');
  return JSON.stringify([errorType, error, result, tail], (_, v) =>
    typeof v === 'bigint' ? String(v) : v
  );
}
