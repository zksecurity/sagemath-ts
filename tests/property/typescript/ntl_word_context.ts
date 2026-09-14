import { RandomStream } from '../../../packages/ntl-ts/src/ZZ.js';
import { zz_pXModulus, build } from '../../../packages/ntl-ts/src/lzz_pX.js';
import {
  GCD,
  MinPolySeq,
  BuildFromRoots,
  MinPolyMod,
} from '../../../packages/ntl-ts/src/lzz_pX1.js';
import {
  SFCanZass,
  EDF,
  FindFactors,
  SFCanZass1,
  NewDDF,
} from '../../../packages/ntl-ts/src/lzz_pXFactoring.js';
export function ntl_word_context(
  op: bigint,
  p: bigint,
  maxroot: bigint,
  d: bigint,
  f: bigint[],
  g: bigint[],
  values: bigint[],
  key: bigint[]
): string {
  const options = { maxroot: Number(maxroot) },
    F = new zz_pXModulus(null, p, options),
    k = F.arithmetic,
    ff = k.norm(f),
    gg = k.norm(g),
    v = values.map(k.mod);
  const pc = F.PrimeCnt,
    profile = [
      pc,
      pc,
      F.MaxRoot,
      F.modCrossover,
      pc === 1 ? 150 : pc === 2 ? 300 : 500,
      pc === 1 ? 90 : pc === 2 ? 180 : 350,
      pc === 1 ? 400 : pc === 2 ? 800 : 1400,
      pc === 1 ? 480 : pc === 2 ? 900 : 1600,
    ];
  const stream = new RandomStream(Uint8Array.from(key, Number));
  let errorType: string | null = null,
    error: string | null = null,
    result: unknown = null;
  try {
    if (op === 0n) {
      build(F, ff);
      result = [F.n, F.n > F.modCrossover + 1 ? 1 : 0];
    } else if (op === 1n) result = GCD(ff, gg, p, options);
    else if (op === 2n) result = MinPolySeq(v, Number(d), p, options);
    else if (op === 3n) result = BuildFromRoots(v, p, options);
    else if (op === 4n) result = SFCanZass(ff, p, stream, options);
    else if (op === 5n) result = EDF(ff, gg, Number(d), p, stream, options);
    else if (op === 6n) {
      build(F, ff);
      result = MinPolyMod(gg, F, Number(d), stream);
    } else if (op === 7n) result = FindFactors(ff, gg, v, p, options);
    else if (op === 8n) result = SFCanZass1(ff, p, options);
    else if (op === 9n) result = [NewDDF(ff, gg, p, options), gg];
    else throw new Error('unknown NTL word context operation');
  } catch (e) {
    errorType = (e as Error).name;
    error = (e as Error).message;
  }
  const tail = Array.from(stream.get(64), (x) => x.toString(16).padStart(2, '0')).join('');
  return JSON.stringify([profile, errorType, error, result, tail], (_, v) =>
    typeof v === 'bigint' ? String(v) : v
  );
}
