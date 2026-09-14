import { ZM_pivots, ZM_rank, ZM_gauss } from '../../../packages/parigp-ts/src/alglin1.js';
export function pari_integer_linear(
  op: bigint,
  m: bigint,
  n: bigint,
  k: bigint,
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
    const [d, r] = ZM_pivots(A);
    z = [d?.slice(1) ?? null, r];
  } else if (op === 1n) z = ZM_rank(A);
  else {
    const X = ZM_gauss(A, B);
    z = X === null ? null : [X[0].slice(1).map((c) => c.slice(1)), X[1]];
  }
  return JSON.stringify(z, (_, v) => (typeof v === 'bigint' ? String(v) : v));
}
