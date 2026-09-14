/** nmod_poly/compose_mod_brent_kung_precomp_preinv.c, dense matrix adapter.
 * @see Deviation: Native FLINT factor array kernels
 */
import { _nmod_poly_xgcd_kernels as k } from './gcd.js';
import { _nmod_poly_divrem } from './divrem.js';
import { _nmod_poly_mul } from './mul.js';
import { _nmod_poly_preinv_remainder } from './powmod_binexp_preinv.js';
import { nmod_mat_mul } from '../nmod_mat/mul.js';
import { _nmod_poly_mod_matrix_rows_evaluate } from './mod_matrix_rows_evaluate.js';
export function nmod_poly_precompute_matrix(
  inner: readonly bigint[],
  f: readonly bigint[],
  finv: readonly bigint[],
  p: bigint
): bigint[][] {
  const F = k.normalized(f, p);
  if (!F.length) throw new RangeError('polynomial division by zero');
  const n = F.length - 1,
    m = Math.floor(Math.sqrt(n)) + 1;
  if (!n) return [[]];
  const g = _nmod_poly_divrem(k.normalized(inner, p), F, p)[1];
  const powers: bigint[][] = [[1n], g];
  for (let i = 2; i < m; i++)
    powers.push(
      _nmod_poly_preinv_remainder(
        _nmod_poly_mul(powers[Math.floor(i / 2)]!, powers[Math.ceil(i / 2)]!, p),
        F,
        finv,
        p
      )
    );
  return powers.slice(0, m).map((row) => Array.from({ length: n }, (_, i) => row[i] ?? 0n));
}
export function _nmod_poly_reduce_matrix_mod_poly(
  rows: readonly (readonly bigint[])[],
  f: readonly bigint[],
  p: bigint
): bigint[][] {
  const F = k.normalized(f, p),
    n = F.length - 1,
    m = Math.floor(Math.sqrt(n)) + 1;
  if (n <= 0) throw new RangeError('modulus must have positive degree');
  return Array.from({ length: m }, (_, i) => {
    const row = i === 0 ? [1n] : _nmod_poly_divrem(k.normalized(rows[i]!, p), F, p)[1];
    return Array.from({ length: n }, (_, j) => row[j] ?? 0n);
  });
}
export function nmod_poly_compose_mod_brent_kung_precomp_preinv(
  outer: readonly bigint[],
  rows: readonly (readonly bigint[])[],
  f: readonly bigint[],
  finv: readonly bigint[],
  p: bigint
): bigint[] {
  const F = k.normalized(f, p),
    a = k.normalized(outer, p),
    n = F.length - 1;
  if (!F.length) throw new RangeError('polynomial division by zero');
  if (a.length >= F.length)
    throw new RangeError('outer degree must be smaller than modulus degree');
  if (!a.length || !n) return [];
  if (a.length === 1) return a;
  const m = Math.floor(Math.sqrt(n)) + 1;
  if (rows.length !== m || rows.some((r) => r.length !== n))
    throw new RangeError('incompatible matrix dimensions');
  const blocks = Array.from({ length: m }, (_, i) =>
    Array.from({ length: m }, (_, j) => a[i * m + j] ?? 0n)
  );
  const products = nmod_mat_mul(blocks, rows, p);
  const h = _nmod_poly_preinv_remainder(
    _nmod_poly_mul(rows[Math.floor(m / 2)]!, rows[m - Math.floor(m / 2)]!, p),
    F,
    finv,
    p
  );
  return _nmod_poly_mod_matrix_rows_evaluate(products, h, F, finv, p);
}
