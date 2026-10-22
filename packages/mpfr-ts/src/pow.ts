import { roundedDyadic, singular } from './arithmetic.js';
import { mpfr_cmp } from './cmp.js';
import { rounding, unsupported } from './round_raw.js';
import type { mpfr_rnd_t, mpfr_t } from './types.js';

/** Largest exact integer-power result materialized as a BigInt, in bits. */
const EXACT_POWER_BITS = 1 << 22;
/** Working-precision ceiling for the Ziv loop. */
const MAX_WORKING_BITS = 1 << 16;

/** A finite nonzero value as `sign * odd * 2^exp`. */
interface Dyadic {
  sign: 1 | -1;
  odd: bigint;
  exp: bigint;
}

function dyadic(x: mpfr_t): Dyadic {
  let odd = x.mantissa;
  let exp = BigInt(x.exponent - x.precision);
  while ((odd & 1n) === 0n) {
    odd >>= 1n;
    exp++;
  }
  return { sign: x.sign, odd, exp };
}

function bitLength(n: bigint): number {
  return n === 0n ? 0 : (n < 0n ? -n : n).toString(2).length;
}

function safeExponent(e: bigint): number {
  if (e > BigInt(Number.MAX_SAFE_INTEGER) || e < BigInt(Number.MIN_SAFE_INTEGER))
    unsupported('exponent outside safe integer range');
  return Number(e);
}

/** Integer r-th root of n >= 0, rounded down. */
function iroot(n: bigint, r: bigint): bigint {
  if (n < 2n) return n;
  let x = 1n << BigInt(Math.ceil(bitLength(n) / Number(r)));
  for (;;) {
    const y = ((r - 1n) * x + n / x ** (r - 1n)) / r;
    if (y >= x) return x;
    x = y;
  }
}

/** Round `num / den * 2^exp` (num, den > 0) into `z` with a sticky bit. */
function roundedQuotient(z: mpfr_t, num: bigint, den: bigint, exp: bigint, sign: 1 | -1): number {
  const shift = BigInt(Math.max(0, z.precision + 2 + bitLength(den) - bitLength(num)));
  const scaled = num << shift;
  const q = scaled / den;
  const sticky = scaled % den === 0n ? 0n : 1n;
  return roundedDyadic(z, 2n * q + sticky, safeExponent(exp - shift - 1n), sign);
}

/** `S * atanh(p / q)` for 0 <= p/q <= 1/3 as [value, error bound], both in ulps of 1/S. */
function atanhScaled(p: bigint, q: bigint, W: bigint): [bigint, bigint] {
  if (p === 0n) return [0n, 0n];
  const p2 = p * p;
  const q2 = q * q;
  let power = (p << W) / q;
  let powerError = 1n;
  let sum = power;
  let error = 1n;
  for (let i = 1n; ; i++) {
    power = (power * p2) / q2;
    powerError++;
    if (power === 0n) break;
    sum += power / (2n * i + 1n);
    error += powerError + 1n;
  }
  // The truncated tail is below powerError * sum(s^2k) <= 9/8 * powerError.
  return [sum, error + 2n * powerError];
}

/** `x^(Y / 2^j)` for finite x > 0, odd Y and j >= 1 that is not a dyadic number. */
function powZiv(z: mpfr_t, x: mpfr_t, Y: bigint, j: bigint): number {
  const half = 1n << BigInt(x.precision - 1);
  const k = BigInt(x.exponent - 1);
  for (let bits = z.precision + 64; bits <= MAX_WORKING_BITS; bits *= 2) {
    const W = BigInt(bits);
    const S = 1n << W;
    // x = f * 2^k with f in [1, 2), so log2(x) = k + ln(f) / ln(2).
    const [l3, e3] = atanhScaled(1n, 3n, W);
    const vL = 2n * l3;
    const eL = 2n * e3;
    const [lf, ef] = atanhScaled(x.mantissa - half, x.mantissa + half, W);
    const va = 2n * lf;
    const ea = 2n * ef;
    const vLog = (va << W) / vL;
    const eLog = ((ea * vL + va * eL) * S) / (vL * (vL - eL)) + 2n;
    // t = y * log2(x), split as n + phi with phi in [0, 1).
    const T = ((k << W) + vLog) * Y;
    const vT = T >> j;
    const eT = ((eLog * (Y < 0n ? -Y : Y)) >> j) + 2n;
    let n = vT >> W;
    const phi = vT - (n << W);
    let vU = (phi * vL) >> W;
    let eU = ((eT * vL + phi * eL + eT * eL) >> W) + 2n;
    if (vU < eU) {
      n--;
      vU += vL;
      eU += eL;
    }
    // exp(u) by its Taylor series; each truncated term adds at most one ulp.
    let term = S;
    let sum = S;
    let terms = 0n;
    for (let i = 1n; term > 0n; i++) {
      term = (term * vU) / (i << W);
      sum += term;
      terms++;
    }
    const eE = (terms * (terms + 1n)) / 2n + 2n * terms + 2n + 3n * eU;
    const lo = { ...z };
    const hi = { ...z };
    const e2 = safeExponent(n - W);
    roundedDyadic(lo, sum - eE, e2, 1);
    roundedDyadic(hi, sum + eE, e2, 1);
    if (mpfr_cmp(lo, hi) !== 0) continue;
    // The exact value is not dyadic, so a tight enough interval excludes it.
    const scaled = (v: bigint, e: number): mpfr_t => ({
      kind: 'finite',
      sign: 1,
      precision: Math.max(1, bitLength(v)),
      mantissa: v,
      exponent: e + bitLength(v),
    });
    const ternary =
      mpfr_cmp(lo, scaled(sum - eE, e2)) < 0 ? -1 : mpfr_cmp(lo, scaled(sum + eE, e2)) > 0 ? 1 : 0;
    if (ternary === 0) continue;
    Object.assign(z, lo);
    return ternary;
  }
  return unsupported('mpfr_pow working precision');
}

/** MPFR 4.2.1 pow.c: correctly rounded `x^y`, including the C99 special values.
 * Exact and integer powers are materialized; other powers use a Ziv loop over
 * interval bounds of exp(y * log(x)) in fixed-point BigInt arithmetic.
 * Supported: nearest-even rounding and exact intermediate powers up to 2^22 bits.
 * @see Deviation: Native MPFR power
 */
export function mpfr_pow(z: mpfr_t, x: mpfr_t, y: mpfr_t, mode: mpfr_rnd_t = 'RNDN'): number {
  rounding(mode);
  const one = (): number => roundedDyadic(z, 1n, 0, 1);
  if (y.kind === 'zero') return one();
  if (
    x.kind === 'finite' &&
    x.sign === 1 &&
    x.mantissa === 1n << BigInt(x.precision - 1) &&
    x.exponent === 1
  )
    return one();
  if (x.kind === 'nan' || y.kind === 'nan') return singular(z, 'nan');

  const yd = y.kind === 'finite' ? dyadic(y) : null;
  const yInteger = yd !== null && yd.exp >= 0n;
  const yOddInteger = yd !== null && yd.exp === 0n;
  const yPositive = y.sign === 1;

  if (y.kind === 'inf') {
    if (x.kind === 'inf') return singular(z, yPositive ? 'inf' : 'zero', 1);
    if (x.kind === 'zero') return singular(z, yPositive ? 'zero' : 'inf', 1);
    const xd = dyadic(x);
    // |x| = 1 here only for x = -1.
    if (xd.odd === 1n && xd.exp === 0n) return one();
    return singular(z, x.exponent >= 1 === yPositive ? 'inf' : 'zero', 1);
  }
  if (x.kind === 'zero' || x.kind === 'inf') {
    const sign = x.sign === -1 && yOddInteger ? -1 : 1;
    const large = (x.kind === 'inf') === yPositive;
    return singular(z, large ? 'inf' : 'zero', sign);
  }

  const xd = dyadic(x);
  if (x.sign === -1 && !yInteger) return singular(z, 'nan');
  const Y = yd!.sign === 1 ? yd!.odd : -yd!.odd;

  if (yInteger) {
    const n = Y << yd!.exp;
    const sign = (x.sign === -1 && (n & 1n) === 1n ? -1 : 1) as 1 | -1;
    const magnitude = n < 0n ? -n : n;
    if (xd.odd !== 1n && BigInt(bitLength(xd.odd)) * magnitude > BigInt(EXACT_POWER_BITS))
      unsupported('mpfr_pow exact power too large');
    const power = xd.odd === 1n ? 1n : xd.odd ** magnitude;
    if (n >= 0n) return roundedDyadic(z, power, safeExponent(xd.exp * n), sign);
    return roundedQuotient(z, 1n, power, xd.exp * n, sign);
  }

  // y = Y / 2^j with Y odd: x^y is dyadic iff the odd part of x is a 2^j-th
  // power (with Y > 0, or that part is 1) and 2^j divides the scaled exponent.
  const j = -yd!.exp;
  const root = 1n << j;
  if ((xd.exp * Y) % root === 0n) {
    const exp = (xd.exp * Y) / root;
    if (xd.odd === 1n) return roundedDyadic(z, 1n, safeExponent(exp), 1);
    if (Y > 0n && j <= 5n) {
      const w = iroot(xd.odd, root);
      if (w ** root === xd.odd) {
        if (BigInt(bitLength(w)) * Y > BigInt(EXACT_POWER_BITS))
          unsupported('mpfr_pow exact power too large');
        return roundedDyadic(z, w ** Y, safeExponent(exp), 1);
      }
    }
  }
  return powZiv(z, x, Y, j);
}
