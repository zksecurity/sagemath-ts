import { Flm_pivots, Flm_gauss } from '../../../packages/parigp-ts/src/Flv.js';
export function pari_word_linear(
  op: bigint,
  m: bigint,
  n: bigint,
  k: bigint,
  p: bigint,
  flat: bigint[],
  other: bigint[]
): string {
  const matrix = (m: number, n: number, flat: bigint[]) => [
    [],
    ...Array.from({ length: n }, (_, j) => [
      0n,
      ...Array.from({ length: m }, (_, i) => flat[i * n + j]!),
    ]),
  ];
  const A = matrix(Number(m), Number(n), flat),
    B = matrix(Number(m), Number(k), other);
  let z;
  if (op === 0n) {
    const [d, r] = Flm_pivots(A, p);
    z = [d?.slice(1) ?? null, r];
  } else {
    const X = Flm_gauss(A, B, p);
    z = X?.slice(1).map((c) => c.slice(1)) ?? null;
  }
  return JSON.stringify(z, (_, v) => (typeof v === 'bigint' ? String(v) : v));
}
