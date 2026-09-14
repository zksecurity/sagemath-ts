import { polsym, polsym_gen } from '../../../packages/parigp-ts/src/polarit2.js';
export function pari_newton_sums(
  op: bigint,
  n: bigint,
  N: bigint,
  use: bigint,
  P: bigint[],
  c: bigint[]
): string {
  let result: unknown;
  if (op) result = polsym(P, Number(n));
  else if (N) result = polsym_gen(P, use ? c : null, Number(n), null, N);
  else {
    const y = use
      ? Array.from(
          { length: c.length / 2 },
          (_, i) => [c[2 * i]!, c[2 * i + 1]!] as [bigint, bigint]
        )
      : null;
    result = polsym_gen(P, y, Number(n), null, null);
  }
  return JSON.stringify(result, (_, v) => (typeof v === 'bigint' ? String(v) : v));
}
