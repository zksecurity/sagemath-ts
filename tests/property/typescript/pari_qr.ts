import {
  QR_init,
  R_from_QR,
  gaussred_from_QR,
  type QrScalar,
} from '../../../packages/parigp-ts/src/bibli1.js';
import { itor, shiftr } from '../../../packages/parigp-ts/src/qfb.js';
export function pari_qr(
  op: bigint,
  mode: bigint,
  m: bigint,
  n: bigint,
  p: bigint,
  q: bigint,
  shift: bigint,
  flat: bigint[]
): string {
  const A: QrScalar[][] = Array.from({ length: Number(n) }, (_, j) =>
    Array.from({ length: Number(m) }, (_, i) => {
      const x = flat[i * Number(n) + j]!;
      return mode === 1n || (mode === 2n && (i + j) % 2 !== 0)
        ? shiftr(itor(x, Number(q)), Number(shift))
        : x;
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
        ? QR_init(A, Number(p))
        : op === 1n
          ? R_from_QR(A, Number(p))
          : gaussred_from_QR(A, Number(p))
    )
  );
}
