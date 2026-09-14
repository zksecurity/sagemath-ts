import { qfgaussred_positive, RgM_Cholesky } from '../../../packages/parigp-ts/src/alglin2.js';
import { RgM_inv_upper } from '../../../packages/parigp-ts/src/alglin1.js';
import {
  matrixRealMul,
  matrixRealSub,
  matrixRealNeg,
  matrixRealDiv,
  matrixRealInv,
  type MatrixReal,
} from '../../../packages/parigp-ts/src/_real_matrix.js';
import { itor, shiftr, real_0_bit } from '../../../packages/parigp-ts/src/qfb.js';
export function pari_real_linear(
  op: bigint,
  mode: bigint,
  n: bigint,
  p: bigint,
  q: bigint,
  shift: bigint,
  flat: bigint[]
): string {
  const A: MatrixReal[][] = Array.from({ length: Number(n) }, (_, j) =>
    Array.from({ length: Number(n) }, (_, i) => {
      const x = flat[i * Number(n) + j]!;
      const bits = mode & 2n && j === 0 ? Number(p) : Number(q);
      return mode & 1n && x === 0n
        ? 0n
        : mode & 4n && x === 0n
          ? real_0_bit(Number(shift) - bits)
          : shiftr(itor(x, bits), Number(shift));
    })
  );
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
  return JSON.stringify(
    encode(
      op === 0n
        ? qfgaussred_positive(A)
        : op === 1n
          ? RgM_Cholesky(A, Number(p))
          : op === 2n
            ? RgM_inv_upper(A)
            : op === 3n
              ? matrixRealMul(A[0]![0]!, A[1]![0]!)
              : op === 4n
                ? matrixRealSub(A[0]![0]!, A[1]![0]!)
                : op === 5n
                  ? matrixRealNeg(A[0]![0]!)
                  : op === 6n
                    ? matrixRealDiv(A[0]![0]!, A[1]![0]!)
                    : matrixRealInv(A[0]![0]!)
    )
  );
}
