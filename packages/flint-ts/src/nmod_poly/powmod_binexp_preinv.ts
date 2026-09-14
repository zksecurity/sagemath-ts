/** nmod_poly/powmod_binexp_preinv.c: reuse the reversed-modulus reciprocal.
 * @see Deviation: Polynomial Modular Powers
 */
import { _nmod_poly_mul } from './mul.js';
import { _nmod_poly_mullow } from './mullow.js';
import { _nmod_poly_divrem } from './divrem.js';
import { _nmod_poly_xgcd_kernels as k } from './gcd.js';
export function _nmod_poly_powmod_fmpz_binexp_preinv(
  a: readonly bigint[],
  e: bigint,
  f: readonly bigint[],
  finv: readonly bigint[],
  p: bigint
): bigint[] {
  if (p < 2n) throw new RangeError('modulus must be at least 2');
  if (e < 0n) throw new RangeError('exponent must be nonnegative');
  const F = k.normalized(f, p);
  let A = k.normalized(a, p);
  if (!F.length)
    throw new Error('Exception (nmod_poly_powmod_fmpz_binexp_preinv). Divide by zero.');
  if (F.length === 1) return [];
  if (A.length >= F.length) A = _nmod_poly_divrem(A, F, p)[1];
  if (e === 0n) return [1n];
  if (e === 1n) return A;
  if (!A.length) return [];
  const reduce = (v: bigint[]) => _nmod_poly_preinv_remainder(v, F, finv, p);
  let out = A;
  for (const bit of e.toString(2).slice(1)) {
    out = reduce(_nmod_poly_mul(out, out, p));
    if (bit === '1') out = reduce(_nmod_poly_mul(out, A, p));
  }
  return out;
}
/** Shared Newton quotient for products of reduced operands (degree < 2*deg(f)). */
export function _nmod_poly_preinv_remainder(
  a: readonly bigint[],
  f: readonly bigint[],
  finv: readonly bigint[],
  p: bigint
): bigint[] {
  if (a.length < f.length) return a.slice();
  const n = a.length - f.length + 1,
    q = _nmod_poly_mullow(a.slice().reverse().slice(0, n), finv, n, p);
  while (q.length < n) q.push(0n);
  q.reverse();
  const product = _nmod_poly_mullow(q, f, f.length - 1, p);
  return k.normalized(
    Array.from({ length: f.length - 1 }, (_, i) => (a[i] ?? 0n) - (product[i] ?? 0n)),
    p
  );
}
