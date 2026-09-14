import {
  drop,
  potential,
  spread,
  condition_bound,
  GS_extraprec,
  gramschmidt_upper,
  gramschmidt_dynprec,
  RgM_Cholesky_dynprec,
} from '../../../packages/parigp-ts/src/lll.js';
import { itor, shiftr } from '../../../packages/parigp-ts/src/qfb.js';
export function pari_lll_gso(
  op: bigint,
  mode: bigint,
  m: bigint,
  n: bigint,
  p: bigint,
  shift: bigint,
  flat: bigint[]
): string {
  const A = Array.from({ length: Number(n) }, (_, j) =>
    Array.from({ length: Number(m) }, (_, i) => flat[i * Number(n) + j]!)
  );
  let z: unknown;
  if (op === 0n) z = gramschmidt_upper(A);
  else if (op === 1n) z = gramschmidt_dynprec(A);
  else if (op === 2n)
    z = RgM_Cholesky_dynprec(A.map((a) => A.map((b) => a.reduce((s, v, i) => s + v * b[i]!, 0n))));
  else {
    const M = mode ? A.map((c) => c.map((v) => shiftr(itor(v, Number(p)), Number(shift)))) : A;
    z = [
      drop(M),
      potential(M),
      spread(M),
      condition_bound(M),
      condition_bound(M, true),
      GS_extraprec(M),
      GS_extraprec(M, true),
    ];
  }
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
  return JSON.stringify(encode(z));
}
