/** fmpz_poly/pow_binomial.c: symmetric binomial coefficient and power recurrence. */
import { _fmpz_poly_pow_small } from './pow_small.js';
export function _fmpz_poly_pow_binomial(a: readonly bigint[], e: bigint): bigint[] {
  if (a.length !== 2 || a[1] === 0n) throw new RangeError('binomial input must have length 2');
  if (e < 0n || e >= 1n << 64n) throw new RangeError('exponent must fit an unsigned word');
  if (e < 3n) return _fmpz_poly_pow_small(a, e);
  const n = Number(e),
    result = Array<bigint>(n + 1).fill(0n);
  result[0] = result[n] = 1n;
  let A = 1n,
    B = 1n,
    C = 1n,
    i = 1,
    f = n - 1;
  for (; i <= Math.floor((n - 1) / 2); i++, f--) {
    A *= a[0]!;
    B *= a[1]!;
    C = (C * BigInt(f + 1)) / BigInt(i);
    result[i] = B * C;
    result[f] = A * C;
  }
  if ((e & 1n) === 0n) {
    A *= a[0]!;
    B *= a[1]!;
    C = (C * BigInt(f + 1)) / BigInt(i);
    result[i] = B * C * A;
    i++;
    f--;
  }
  for (; i <= n; i++, f--) {
    A *= a[0]!;
    B *= a[1]!;
    result[i] = result[i]! * B;
    result[f] = result[f]! * A;
  }
  return result;
}
