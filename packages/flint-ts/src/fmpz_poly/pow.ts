/** fmpz_poly/pow.c: small, binomial, multinomial and binary dispatch. */
import { _fmpz_poly_resultant_kernels as k } from './gcd.js';
import { _fmpz_poly_pow_small } from './pow_small.js';
import { _fmpz_poly_pow_binomial } from './pow_binomial.js';
import { _fmpz_poly_pow_multinomial } from './pow_multinomial.js';
import { _fmpz_poly_pow_binexp } from './pow_binexp.js';
export function _fmpz_poly_pow(a: readonly bigint[], e: bigint): bigint[] {
  if (e < 0n || e >= 1n << 64n) throw new RangeError('exponent must fit an unsigned word');
  const A = k.normalized([...a]);
  if (e < 5n) return _fmpz_poly_pow_small(A, e);
  if (!A.length) return [];
  if (A.length === 1) return [A[0]! ** e];
  if (A.length === 2) return _fmpz_poly_pow_binomial(A, e);
  const limbs = BigInt(Math.ceil(k.maxBits(A) / 64));
  return limbs < ((3n * e) / 2n + 150n) / BigInt(A.length)
    ? _fmpz_poly_pow_multinomial(A, e)
    : _fmpz_poly_pow_binexp(A, e);
}
