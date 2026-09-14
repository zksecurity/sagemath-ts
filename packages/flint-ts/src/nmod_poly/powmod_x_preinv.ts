/** nmod_poly/powmod_x_preinv.c: native shift windows for powers of x.
 * @see Deviation: Polynomial Modular Powers
 */
import { _nmod_poly_mul } from './mul.js';
import {
  _nmod_poly_powmod_fmpz_binexp_preinv,
  _nmod_poly_preinv_remainder,
} from './powmod_binexp_preinv.js';
import { _nmod_poly_xgcd_kernels as k } from './gcd.js';
export function _nmod_poly_powmod_x_fmpz_preinv(
  e: bigint,
  f: readonly bigint[],
  finv: readonly bigint[],
  p: bigint
): bigint[] {
  if (p < 2n) throw new RangeError('modulus must be at least 2');
  if (e < 0n) throw new RangeError('exponent must be nonnegative');
  const F = k.normalized(f, p);
  if (!F.length) throw new Error('Exception (nmod_poly_powmod_x_fmpz_preinv). Divide by zero.');
  if (F.length <= 2 || e <= 2n)
    return _nmod_poly_powmod_fmpz_binexp_preinv([0n, 1n], e, F, finv, p);
  const reduce = (v: bigint[]) => _nmod_poly_preinv_remainder(v, F, finv, p);
  let out = [1n],
    l = BigInt(F.length - 1).toString(2).length - 2,
    window = 2 ** l,
    c = l;
  const bits = e.toString(2);
  let i = bits.length - 2;
  if (i <= l) {
    window = 2 ** i;
    c = i;
    l = i;
  }
  if (c === 0) {
    out = reduce(Array<bigint>(window).fill(0n).concat(out));
    c = l + 1;
    window = 0;
  }
  for (; i >= 0; i--) {
    out = reduce(_nmod_poly_mul(out, out, p));
    c--;
    if (bits[bits.length - 1 - i] === '1') {
      if (window === 0 && i <= l - 1) c = i;
      if (c >= 0) window += 2 ** c;
    } else if (window === 0) c = l + 1;
    if (c === 0) {
      out = reduce(Array<bigint>(window).fill(0n).concat(out));
      c = l + 1;
      window = 0;
    }
  }
  return out;
}
