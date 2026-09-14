/** Native real LLL precision, state, rounding and absolute comparison. */
import { fplll, fplll_heuristic } from '../../../packages/parigp-ts/src/lll.js';
import { roundr_safe } from '../../../packages/parigp-ts/src/gen3.js';
import { abscmprr } from '../../../packages/parigp-ts/src/kernel/none/cmp.js';

const columns = (m: bigint, n: bigint, flat: bigint[]) =>
  Array.from({ length: Number(n) }, (_, j) =>
    Array.from({ length: Number(m) }, (_, i) => flat[i * Number(n) + j]!)
  );
const encode = (value: unknown): string =>
  JSON.stringify(value, (_, x) => (typeof x === 'bigint' ? String(x) : x));

export function pari_lll_heuristic(
  m: bigint,
  n: bigint,
  dn: bigint,
  dd: bigint,
  en: bigint,
  ed: bigint,
  keep: bigint,
  track: bigint,
  precision: bigint,
  gramPrecision: bigint,
  flat: bigint[]
): string {
  return encode(
    fplll_heuristic(
      columns(m, n, flat),
      Number(dn) / Number(dd),
      Number(en) / Number(ed),
      !!keep,
      !!track,
      Number(precision),
      Number(gramPrecision)
    )
  );
}

export function pari_lll_proved(
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
  precision: bigint,
  flat: bigint[]
): string {
  const B = columns(m, n, flat);
  const G = mode ? B.map((a) => B.map((b) => a.reduce((sum, x, i) => sum + x * b[i]!, 0n))) : null;
  const [status, gram, basis, U, norms] = fplll(
    mode === 2n ? null : B,
    G,
    Number(dn) / Number(dd),
    Number(en) / Number(ed),
    !!keep,
    !!track,
    !!want,
    Number(precision)
  );
  return encode([
    status,
    gram,
    basis,
    U,
    norms?.map((x) => [x.s, String(x.e), String(x.m), x.p]) ?? null,
  ]);
}

export function pari_lll_real_scalar(
  op: bigint,
  p: bigint,
  e: bigint,
  m: bigint,
  s: bigint,
  q: bigint,
  f: bigint,
  n: bigint,
  t: bigint
): string {
  const x = { p: Number(p), e: Number(e), m, s: Number(s) as -1 | 0 | 1 };
  const y = { p: Number(q), e: Number(f), m: n, s: Number(t) as -1 | 0 | 1 };
  return String(op ? abscmprr(x, y) : roundr_safe(x));
}
