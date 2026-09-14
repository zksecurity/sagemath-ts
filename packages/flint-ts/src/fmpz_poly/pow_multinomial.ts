/** fmpz_poly/pow_multinomial.c: exact coefficient recurrence with valuation removal. */
import { _fmpz_poly_resultant_kernels as k } from './gcd.js';
import { _fmpz_poly_pow_small } from './pow_small.js';
export function _fmpz_poly_pow_multinomial(a: readonly bigint[], e: bigint): bigint[] {
  if (e < 0n || e >= 1n << 64n) throw new RangeError('exponent must fit an unsigned word');
  const A = k.normalized([...a]);
  if (e < 3n) return _fmpz_poly_pow_small(A, e);
  if (!A.length) return [];
  if (A.length === 1) return [A[0]! ** e];
  const low = A.findIndex((c) => c !== 0n),
    P = A.slice(low);
  const offset = Number(e) * low,
    rlen = Number(e) * (P.length - 1) + 1;
  const result = Array<bigint>(rlen).fill(0n);
  result[0] = P[0]! ** e;
  let denominator = 0n;
  for (let j = 1; j < rlen; j++) {
    let u = -BigInt(j),
      sum = 0n;
    for (let i = 1; i <= Math.min(j, P.length - 1); i++) {
      u += e + 1n;
      sum += P[i]! * result[j - i]! * u;
    }
    denominator += P[0]!;
    result[j] = sum / denominator;
  }
  return Array<bigint>(offset).fill(0n).concat(result);
}
