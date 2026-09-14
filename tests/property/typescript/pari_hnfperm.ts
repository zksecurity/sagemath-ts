import { ZM_hnfperm, ZM_hnf_knapsack, hnfperm } from '../../../packages/parigp-ts/src/hnf_snf.js';
export function pari_hnfperm(op: bigint, m: bigint, n: bigint, flat: bigint[]): string {
  const A = Array.from({ length: Number(n) }, (_, j) =>
    Array.from({ length: Number(m) }, (_, i) => flat[i * Number(n) + j]!)
  );
  const out =
    op === 4n
      ? ZM_hnf_knapsack(A)
      : op === 5n
        ? hnfperm(A)
        : ZM_hnfperm(A, !!(op & 1n), !!(op & 2n));
  return JSON.stringify(out, (_, v) =>
    typeof v === 'bigint' || typeof v === 'number' ? String(v) : v
  );
}
