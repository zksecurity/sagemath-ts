import {
  flat,
  ZM_flatter,
  ZM_flatter_rank,
  ZM_flattergram,
  ZM_lll,
  lllfp,
} from '../../../packages/parigp-ts/src/lll.js';
import { integerMatrixPivots } from '../../../packages/parigp-ts/src/_matrix_inverse.js';
import {
  itor,
  shiftr,
  mkqfb,
  redimagsl2,
  type MpReal,
} from '../../../packages/parigp-ts/src/qfb.js';
export function native_pari_lll_wrapper(
  op: bigint,
  flag: bigint,
  mode: bigint,
  m: bigint,
  n: bigint,
  p: bigint,
  shift: bigint,
  values: bigint[]
): string {
  const rows = Number(m),
    cols = Number(n);
  const A = Array.from({ length: cols }, (_, j) =>
    Array.from({ length: rows }, (_, i) => {
      const x = values[i * cols + j]!;
      return mode === 1n || (mode === 2n && (i + j) % 2)
        ? shiftr(itor(x, Number(p)), Number(shift))
        : x;
    })
  );
  const f = Number(flag);
  let r: unknown;
  if (op === 0n) r = ZM_lll(A as bigint[][], 0.99, f);
  else if (op === 1n) r = lllfp(A, 0.99, f);
  else if (op === 2n) {
    const z = flat(A as bigint[][], f);
    r = [z[0], z[1] ?? 0n, z[2], z[3]];
  } else if (op === 3n) r = ZM_flatter(A as bigint[][], f);
  else if (op === 4n) {
    const R = Array.from({ length: rows }, (_, i) => A.map((c) => c[i] as bigint));
    r = ZM_flatter_rank(A as bigint[][], cols - integerMatrixPivots(R, rows, cols)[1], f);
  } else if (op === 5n) r = ZM_flattergram(A as bigint[][], f);
  else {
    const [a, b, c] = values;
    const z = redimagsl2(mkqfb(a!, b!, c!, b! * b! - 4n * a! * c!));
    r = [[z.Q.a, z.Q.b, z.Q.c], z.U[0]!.map((_, i) => z.U.map((c) => c[i]!))];
  }
  const out = (x: any): unknown =>
    typeof x === 'bigint'
      ? String(x)
      : Array.isArray(x)
        ? x.map(out)
        : x && typeof x === 'object'
          ? [x.s, String(x.e), String(x.m), x.p]
          : x;
  return JSON.stringify(out(r));
}
