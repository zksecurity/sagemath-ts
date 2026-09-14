/** Direct comparisons with the bundled PARI C matrix and half-GCD routines. */
import {
  integerMatrixInverse,
  rationalLift,
  wordMatrixAdjoint,
} from '../../../packages/parigp-ts/src/_matrix_inverse.js';
import { hnf_divscale } from '../../../packages/parigp-ts/src/hnf_snf.js';
import { halfgcdii } from '../../../packages/parigp-ts/src/kernel/none/halfgcd.js';

export function pari_ideal_matrix(
  op: bigint,
  p: bigint,
  m: bigint,
  n: bigint,
  flat: bigint[]
): string {
  const nr = Number(m),
    nc = Number(n);
  let result: unknown;
  if (op === 2n) result = halfgcdii(flat[0]!, flat[1]!);
  else if (op === 3n) result = rationalLift(flat[0]!, flat[1]!, flat[2]!);
  else if (op === 4n) {
    const columns = (start: number) =>
      Array.from({ length: nc }, (_, j) =>
        Array.from({ length: nr }, (_, i) => flat[start + i * nc + j]!)
      );
    const C = hnf_divscale(columns(0), columns(nr * nc), flat[2 * nr * nc]!);
    result = [Array.from({ length: nr }, (_, i) => C.map((c) => c[i]!)), 1n];
  } else {
    const A = Array.from({ length: nr }, (_, i) => flat.slice(i * nc, (i + 1) * nc));
    result = op === 0n ? [wordMatrixAdjoint(A, p), 1n] : integerMatrixInverse(A);
  }
  return JSON.stringify(result, (_, v) => (typeof v === 'bigint' ? String(v) : v));
}
