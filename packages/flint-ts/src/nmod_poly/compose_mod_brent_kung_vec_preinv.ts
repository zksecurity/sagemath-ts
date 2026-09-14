/** Shared matrix product from nmod_poly/compose_mod_brent_kung_vec_preinv.c.
 * @see Deviation: Native FLINT factor array kernels
 */
import { _nmod_poly_xgcd_kernels as k } from './gcd.js';
import { _nmod_poly_mul } from './mul.js';
import { _nmod_poly_divrem } from './divrem.js';
import { _nmod_poly_preinv_remainder } from './powmod_binexp_preinv.js';
import { nmod_mat_mul } from '../nmod_mat/mul.js';
import { _nmod_poly_mod_matrix_rows_evaluate } from './mod_matrix_rows_evaluate.js';
export function nmod_poly_compose_mod_brent_kung_vec_preinv(
  polys: readonly (readonly bigint[])[],
  count: number,
  inner: readonly bigint[],
  f: readonly bigint[],
  finv: readonly bigint[],
  p: bigint
): bigint[][] {
  const F = k.normalized(f, p),
    a = polys.map((v) => k.normalized(v, p));
  if (a.some((v) => v.length >= F.length))
    throw new RangeError('outer degree must be smaller than modulus degree');
  if (!Number.isSafeInteger(count) || count < 0 || count > a.length)
    throw new RangeError('invalid polynomial count');
  if (!count) return [];
  if (F.length === 1) return Array.from({ length: count }, () => []);
  if (F.length === 2) return a.slice(0, count);
  const n = F.length - 1,
    m = Math.floor(Math.sqrt(n * count)) + 1,
    blocksPerPoly = Math.floor(F.length / m) + 1;
  const reduce = (v: bigint[]) => _nmod_poly_preinv_remainder(v, F, finv, p);
  const powers: bigint[][] = [[1n], _nmod_poly_divrem(k.normalized(inner, p), F, p)[1]];
  for (let i = 2; i < m; i++)
    powers.push(reduce(_nmod_poly_mul(powers[Math.floor(i / 2)]!, powers[Math.ceil(i / 2)]!, p)));
  const A = powers.map((v) => Array.from({ length: n }, (_, i) => v[i] ?? 0n));
  const B = Array.from({ length: blocksPerPoly * count }, (_, i) =>
    Array.from(
      { length: m },
      (_, j) => a[Math.floor(i / blocksPerPoly)]![(i % blocksPerPoly) * m + j] ?? 0n
    )
  );
  const C = nmod_mat_mul(B, A, p);
  const h = reduce(_nmod_poly_mul(A[m - 1]!, A[1]!, p));
  return Array.from({ length: count }, (_, i) =>
    _nmod_poly_mod_matrix_rows_evaluate(
      C.slice(i * blocksPerPoly, (i + 1) * blocksPerPoly),
      h,
      F,
      finv,
      p
    )
  );
}
