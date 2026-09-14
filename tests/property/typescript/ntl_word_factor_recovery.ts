import { RandomStream } from '../../../packages/ntl-ts/src/ZZ.js';
import { zz_pXModulus } from '../../../packages/ntl-ts/src/lzz_pX.js';
import {
  FindRoots,
  FindRoot,
  RootEDF,
  FindFactors,
  EDFSplit,
  EDF,
  SFCanZass2,
  SFCanZass,
} from '../../../packages/ntl-ts/src/lzz_pXFactoring.js';
export function ntl_word_factor_recovery(
  op: bigint,
  p: bigint,
  d: bigint,
  f: bigint[],
  b: bigint[],
  roots: bigint[],
  packed: bigint[],
  key: bigint[]
): string {
  const k = new zz_pXModulus(null, p).arithmetic,
    F = k.norm(f),
    B = k.norm(b),
    R = roots.map(k.mod),
    groups: [bigint[], number][] = [];
  for (let i = 0; i < packed.length; ) {
    const n = Number(packed[i++]),
      degree = Number(packed[i++]);
    groups.push([k.norm(packed.slice(i, i + n)), degree]);
    i += n;
  }
  const stream = new RandomStream(Uint8Array.from(key, Number));
  let errorType: string | null = null,
    error: string | null = null,
    result: unknown = null;
  try {
    if (op === 0n) result = FindRoots(F, p, stream);
    else if (op === 1n) result = FindRoot(F, p, stream);
    else if (op === 2n) result = RootEDF(F, p, stream);
    else if (op === 3n) result = FindFactors(F, B, R, p);
    else if (op === 4n) result = EDFSplit(F, B, Number(d), p, stream);
    else if (op === 5n) result = EDF(F, B, Number(d), p, stream);
    else if (op === 6n) result = SFCanZass2(groups, B, p, stream);
    else if (op === 7n) result = SFCanZass(F, p, stream);
    else if (op === 8n) result = k.divrem(F, B);
    else throw new Error('unknown factor recovery operation');
  } catch (e) {
    errorType = (e as Error).name;
    error = (e as Error).message;
  }
  const tail = Array.from(stream.get(64), (x) => x.toString(16).padStart(2, '0')).join('');
  return JSON.stringify([errorType, error, result, tail], (_, v) =>
    typeof v === 'bigint' ? String(v) : v
  );
}
