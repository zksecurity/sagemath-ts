import { CRT as scalarCRT, CRTInRange } from '../../../packages/ntl-ts/src/ZZ.js';
import { CRT, DetBound } from '../../../packages/ntl-ts/src/mat_ZZ.js';
import { determinant } from '../../../packages/ntl-ts/src/mat_ZZ_p.js';
export function ntl_integer_crt(
  op0: bigint,
  a: bigint,
  p: bigint,
  dims: bigint[],
  g: bigint[],
  G: bigint[]
): string {
  const op = Number(op0),
    [n, m, r, c, badG, badR] = dims.map(Number);
  const make = (v: bigint[], rows: number, cols: number, bad: number): bigint[][] =>
    bad
      ? [[0n], [0n, 0n]]
      : Array.from({ length: Math.max(0, rows) }, (_, i) => v.slice(i * cols, (i + 1) * cols));
  let result: unknown;
  if (op <= 2) {
    result = g.map((gg, i) => {
      try {
        return [
          null,
          null,
          op === 0 ? CRTInRange(gg, a) : scalarCRT(gg, a, G[i]!, p, { word: op === 2 }),
        ];
      } catch (e) {
        return [(e as Error).name, (e as Error).message, null];
      }
    });
  } else {
    if (op === 3) {
      if (p <= 1n) throw new Error('zz_pContext: p must be > 1');
      if (p >= 1n << 60n) throw new Error('zz_pContext: modulus too big');
    }
    if (op === 5 && p <= 1n) throw new Error('ZZ_pContext: p must be > 1');
    if (n! < 0 || m! < 0 || (op === 3 && (r! < 0 || c! < 0))) throw new Error('SetDims: bad args');
    const M = make(g, n!, m!, badG!);
    if (op === 3) result = CRT(M, a, make(G, r!, c!, badR!), p, { columns: m, residueColumns: c });
    else if (op === 4) result = DetBound(M, { columns: m });
    else if (op === 5) result = determinant(M, p, { columns: m });
    else throw new Error('unknown integer CRT operation');
  }
  return JSON.stringify(result, (_, v) => (typeof v === 'bigint' ? String(v) : v));
}
