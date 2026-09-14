import { fma, frexp, ldexp } from '../../../packages/parigp-ts/src/_binary64.js';
/** Columns, partial failure state and optional U from the original PARI fast stage. */
import { fplll_fast } from '../../../packages/parigp-ts/src/lll.js';

export function pari_lll_fast(
  _op: bigint,
  m: bigint,
  n: bigint,
  dn: bigint,
  dd: bigint,
  en: bigint,
  ed: bigint,
  keep: bigint,
  track: bigint,
  flat: bigint[]
): string {
  if (_op) {
    const view = new DataView(new ArrayBuffer(8));
    const read = (bits: bigint) => {
      view.setBigUint64(0, bits, false);
      return view.getFloat64(0, false);
    };
    const write = (x: number) => {
      if (Number.isNaN(x)) return 'NaN';
      view.setFloat64(0, x, false);
      return String(view.getBigUint64(0, false));
    };
    const a = read(flat[0]!);
    if (_op === 1n) return JSON.stringify([write(fma(a, read(flat[1]!), read(flat[2]!)))]);
    if (_op === 2n) return JSON.stringify([write(ldexp(a, Number(flat[1]!)))]);
    const [fraction, exponent] = frexp(a);
    return JSON.stringify([write(fraction), exponent]);
  }
  const rows = Number(m),
    columns = Number(n);
  const B = Array.from({ length: columns }, (_, j) =>
    Array.from({ length: rows }, (_, i) => flat[i * columns + j]!)
  );
  const result = fplll_fast(B, Number(dn) / Number(dd), Number(en) / Number(ed), !!keep, !!track);
  return JSON.stringify(result, (_, value) => (typeof value === 'bigint' ? String(value) : value));
}
