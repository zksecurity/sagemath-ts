import { lllgramint } from '../../../packages/parigp-ts/src/qfrep.ts';
import { ZM_hnflll, ZM_snf_group } from '../../../packages/parigp-ts/src/buch.ts';

export function pari_lll_dependents(
  op: bigint,
  degree: bigint,
  flag: bigint,
  flat: bigint[]
): string {
  const n = Number(degree),
    rows = n ? (flat.length - (op === 3n ? 2 : 0)) / n : 0,
    M = Array.from({ length: n }, (_, j) =>
      Array.from({ length: rows }, (_, i) => flat[i * n + j]!)
    );
  const B = [[], ...M.map((c) => [0n, ...c])];
  const strip = (A: bigint[][] | null) => A?.slice(1).map((c) => c.slice(1)) ?? null;
  let v: unknown;
  if (op === 0n) v = lllgramint(M);
  else if (op === 3n) v = lllgramint(M, { n: flat[flat.length - 2]!, d: flat[flat.length - 1]! });
  else if (op === 4n) v = strip(ZM_hnflll(B, false, flag !== 0n).H);
  else if (op === 1n) {
    const { H, B: U } = ZM_hnflll(B, true, flag !== 0n);
    v = [strip(H), strip(U)];
  } else {
    const { D, Ui } = ZM_snf_group(B);
    v = [D.slice(1), strip(Ui)];
  }
  return JSON.stringify(v, (_, x) => (typeof x === 'bigint' ? String(x) : x));
}
