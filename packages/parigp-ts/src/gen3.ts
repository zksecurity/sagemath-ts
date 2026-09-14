import { rationalPair } from './_rational_polynomial.js';
import type { MpReal } from './qfb.js';
import { PariError } from './errors.js';

/** Native integer quotient rounded to nearest, with ties toward positive infinity. */
export function diviiround(x: bigint, y: bigint): bigint {
  if (y === 0n) throw new PariError('impossible inverse in dvmdii: 0');
  let q = x / y;
  const r = x % y,
    twice = 2n * (r < 0n ? -r : r),
    denominator = y < 0n ? -y : y;
  const positive = x < 0n === y < 0n;
  if (twice > denominator || (twice === denominator && positive)) q += positive ? 1n : -1n;
  return q;
}

/** PARI roundr_safe, allowing integer conversion beyond mantissa precision. */
export function roundr_safe(x: MpReal): bigint {
  // gen3.c:2436-2493 rounds toward +infinity on a half-integer tie.
  if (!x.s || x.e < -1) return 0n;
  const signed = x.s < 0 ? -x.m : x.m;
  const shift = x.e + 1 - x.p;
  if (shift >= 0) return signed << BigInt(shift);
  const bits = BigInt(-shift);
  return (signed + (1n << (bits - 1n))) >> bits;
}

/** Native ceil_safe for integer, rational-pair and real scalars, including accuracy bounds.
 * @see Deviation: PARI factorization scalar and bound adapters
 */
export function ceil_safe(x: bigint | [bigint, bigint] | MpReal): bigint {
  if (typeof x === 'bigint') return x;
  if (Array.isArray(x)) {
    const [a, b] = rationalPair(x[0], x[1]),
      q = a / b,
      r = a % b;
    return q + (r !== 0n && a < 0n === b < 0n ? 1n : 0n);
  }
  const shift = x.e + 1 - x.p;
  const magnitude =
    !x.s || x.e < 0 ? 0n : shift >= 0 ? x.m << BigInt(shift) : x.m >> BigInt(-shift);
  const integral = x.s < 0 ? -magnitude : magnitude;
  return x.s < 0
    ? integral
    : integral + (1n << BigInt(Math.max(0, x.e < 0 ? x.e : shift > 0 ? shift : !x.s ? x.e : -1)));
}
