import { LLL_check_progress } from '../../../packages/parigp-ts/src/QX_factor.js';
import { itor, shiftr } from '../../../packages/parigp-ts/src/qfb.js';
export function pari_lll_progress(
  n0: bigint,
  final: bigint,
  m: bigint,
  n: bigint,
  p: bigint,
  shift: bigint,
  bound: bigint,
  flat: bigint[]
): string {
  const A = Array.from({ length: Number(n) }, (_, j) =>
    Array.from({ length: Number(m) }, (_, i) => flat[i * Number(n) + j]!)
  );
  return JSON.stringify(
    LLL_check_progress(shiftr(itor(bound, Number(p)), Number(shift)), Number(n0), A, !!final),
    (_, v) => (typeof v === 'bigint' ? String(v) : v)
  );
}
