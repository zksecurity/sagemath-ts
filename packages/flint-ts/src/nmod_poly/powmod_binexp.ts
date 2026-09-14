/** nmod_poly/powmod_binexp.c: left-to-right binary modular power.
 * @see Deviation: Polynomial Modular Powers
 */
import { _nmod_poly_mul } from './mul.js';
import { _nmod_poly_divrem } from './divrem.js';
import { _nmod_poly_xgcd_kernels as k } from './gcd.js';
export function _nmod_poly_powmod_ui_binexp(
  a: readonly bigint[],
  e: bigint,
  f: readonly bigint[],
  p: bigint
): bigint[] {
  if (p < 2n) throw new RangeError('modulus must be at least 2');
  if (e < 0n || e >= 1n << 64n) throw new RangeError('exponent must fit an unsigned word');
  const F = k.normalized(f, p);
  let A = k.normalized(a, p);
  if (!F.length) throw new Error('Exception (nmod_poly_powmod_ui_binexp). Divide by zero.');
  if (F.length === 1) return [];
  if (A.length >= F.length) A = _nmod_poly_divrem(A, F, p)[1];
  if (e === 0n) return [1n];
  if (e === 1n) return A;
  if (!A.length) return [];
  let out = A;
  for (const bit of e.toString(2).slice(1)) {
    out = _nmod_poly_divrem(_nmod_poly_mul(out, out, p), F, p)[1];
    if (bit === '1') out = _nmod_poly_divrem(_nmod_poly_mul(out, A, p), F, p)[1];
  }
  return out;
}
