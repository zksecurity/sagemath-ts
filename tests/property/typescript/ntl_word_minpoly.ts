import { MinPolySeq } from '../../../packages/ntl-ts/src/lzz_pX1.js';
export function ntl_word_minpoly(p: bigint, m: bigint, a: bigint[]): string {
  return JSON.stringify(MinPolySeq(a, Number(m), p), (_, v) =>
    typeof v === 'bigint' ? String(v) : v
  );
}
