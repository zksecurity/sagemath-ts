/** @see Deviation: Polynomial Roots and Truncated Series */
/** nmod_poly/pow_trunc.c, including the native zero-polynomial precedence at e=0. */
import { _nmod_poly_mullow } from './mullow.js';
import { _nmod_poly_pow } from './pow.js';
import { _fmpz_poly_resultant_kernels as k } from '../fmpz_poly/gcd.js';
export function _nmod_poly_pow_trunc(
  a: readonly bigint[],
  e: bigint,
  n: number,
  p: bigint
): bigint[] {
  if (e < 0n || e >= 1n << 64n) throw new RangeError('exponent must fit an unsigned word');
  if (!Number.isSafeInteger(n) || n < 0) throw new RangeError('invalid truncation length');
  if (p < 1n || p >= 1n << 64n) throw new RangeError('modulus must fit a positive unsigned word');
  const A = k.normalized(a.map((c) => ((c % p) + p) % p));
  if (!n || !A.length) return [];
  if (A.length === 1) return _nmod_poly_pow(A, e, p);
  if (!e) return [1n];
  const base = A.slice(0, n);
  let result = base;
  for (const bit of e.toString(2).slice(1)) {
    result = _nmod_poly_mullow(result, result, n, p);
    if (bit === '1') result = _nmod_poly_mullow(result, base, n, p);
  }
  return k.normalized(result);
}
