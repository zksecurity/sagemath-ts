import { RgM_rescale_to_int } from '../../../packages/parigp-ts/src/polarit2.js';
import {
  gram_matrix,
  RgM_mul,
  RgV_dotsquare,
  RgV_dotproduct,
  type RgScalar,
} from '../../../packages/parigp-ts/src/RgV.js';
import { itor, shiftr } from '../../../packages/parigp-ts/src/qfb.js';
export function pari_matrix_products(
  op: bigint,
  mode: bigint,
  mode2: bigint,
  m: bigint,
  k: bigint,
  n: bigint,
  p: bigint,
  q: bigint,
  shift: bigint,
  flat: bigint[],
  other: bigint[]
): string {
  const matrix = (
    m: number,
    k: number,
    mode: bigint,
    p: bigint,
    shift: bigint,
    flat: bigint[]
  ): RgScalar[][] =>
    Array.from({ length: k }, (_, j) =>
      Array.from({ length: m }, (_, i) => {
        const x = flat[i * k + j]!;
        return mode === 1n || (mode === 2n && (i + j) % 2) || (mode === 3n && x !== 0n)
          ? shiftr(itor(x, Number(p)), Number(shift))
          : x;
      })
    );
  const A = matrix(Number(m), Number(k), mode, p, shift, flat),
    B = matrix(Number(k), Number(n), mode2, q, -shift, other);
  const encode = (x: any): any =>
    typeof x === 'bigint'
      ? String(x)
      : Array.isArray(x)
        ? x.map(encode)
        : x === null
          ? null
          : typeof x === 'object'
            ? [x.s, String(x.e), String(x.m), x.p]
            : x;
  const z =
    op === 0n
      ? RgM_rescale_to_int(A)
      : op === 1n
        ? gram_matrix(A)
        : op === 2n
          ? RgM_mul(A, B)
          : op === 3n
            ? RgV_dotsquare(A[0]!)
            : op === 4n
              ? RgV_dotproduct(A[0]!, A[1]!)
              : RgV_dotproduct(A[0]!, A[0]!);
  return JSON.stringify(encode(z));
}
