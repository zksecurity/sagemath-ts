/** Exact native DPE stage state; resource failures have a separate dispatcher. */
import { fplll_dpe } from '../../../packages/parigp-ts/src/lll.js';

export function pari_lll_dpe(
  mode: bigint,
  m: bigint,
  n: bigint,
  dn: bigint,
  dd: bigint,
  en: bigint,
  ed: bigint,
  keep: bigint,
  track: bigint,
  want: bigint,
  flat: bigint[]
): string {
  const rows = Number(m),
    columns = Number(n);
  const B = Array.from({ length: columns }, (_, j) =>
    Array.from({ length: rows }, (_, i) => flat[i * columns + j]!)
  );
  const G = mode ? B.map((a) => B.map((b) => a.reduce((sum, x, i) => sum + x * b[i]!, 0n))) : null;
  const [status, gram, basis, U, norms] = fplll_dpe(
    mode === 2n ? null : B,
    G,
    Number(dn) / Number(dd),
    Number(en) / Number(ed),
    !!keep,
    !!track,
    !!want
  );
  const result = [
    status,
    gram,
    basis,
    U,
    norms?.map((x) => [x.s, String(x.e), String(x.m), x.p]) ?? null,
  ];
  return JSON.stringify(result, (_, value) => (typeof value === 'bigint' ? String(value) : value));
}

export function pari_lll_dpe_resource(...args: Parameters<typeof pari_lll_dpe>): string {
  try {
    pari_lll_dpe(...args);
  } catch (error) {
    if (
      error instanceof RangeError &&
      error.message === 'fplll_dpe: Gram-Schmidt coefficient requires an unrepresentable shift'
    )
      return 'resource_failure';
    throw error;
  }
  return 'unexpected_success';
}
