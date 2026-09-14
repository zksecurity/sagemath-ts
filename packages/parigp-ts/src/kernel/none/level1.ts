import { divri, divru, itor, rtor, shiftr, type MpReal } from '../../qfb.js';

/** PARI rdiviiz/rdivii: integer quotient with native conversion and guard order. */
export function rdivii(x: bigint, y: bigint, p: number): MpReal {
  const bits = (n: bigint) => (n === 0n ? 0 : (n < 0n ? -n : n).toString(2).length);
  if (x === 0n) return itor(0n, p);
  if (y === 0n) throw new RangeError('rdivii requires a nonzero denominator');
  const bx = bits(x),
    by = bits(y);
  if (by <= 64) return rtor(divri(itor(x, p), y), p);
  if (Math.ceil(bx / 64) * 64 > p + 64 || Math.ceil(by / 64) * 64 > p + 64)
    return rtor(divri(itor(x, p), y), p);
  const shift = p + by - bx + 1;
  const q = (shift > 0 ? x << BigInt(shift) : x) / y;
  const result = itor(q, p);
  return shift > 0 ? shiftr(result, -shift) : result;
}

import { cmprr } from './cmp.js';
/** Native integer/real comparison rounds the integer to the real's precision.
 * @see Deviation: PARI factorization scalar and bound adapters
 */
export function cmpir(x: bigint, y: MpReal): number {
  if (!x) return y.s === 0 ? 0 : -y.s;
  if (!y.s) {
    const e = (x < 0n ? -x : x).toString(2).length - 1;
    return y.e >= e ? 0 : x < 0n ? -1 : 1;
  }
  return cmprr(itor(x, y.p), y);
}
export function cmpri(x: MpReal, y: bigint): number {
  const c = cmpir(y, x);
  return c === 0 ? 0 : -c;
}
