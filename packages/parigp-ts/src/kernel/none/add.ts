import { type MpReal, real_0_bit } from '../../qfb.js';
/**
 * PARI addrr_sign restricted to one 64-bit mantissa word.
 * Native alignment truncates the smaller mantissa before addition; cancellation
 * normalizes the result without reducing its whole-word storage precision.
 * @see Deviation: Native Single-Word Real Elimination
 */
export function addrr(x: MpReal, y: MpReal): MpReal {
  if ((x.s && x.p !== 64) || (y.s && y.p !== 64))
    throw new Error('PARI_NOT_IMPLEMENTED: multiword real arithmetic');
  if (!x.s) return !y.s ? real_0_bit(Math.max(x.e, y.e)) : x.e >= y.e ? real_0_bit(x.e) : { ...y };
  if (!y.s) return y.e >= x.e ? real_0_bit(y.e) : { ...x };
  if (x.e > y.e) [x, y] = [y, x];
  const shift = y.e - x.e;
  if (shift >= 64) return { ...y };
  const a = x.m >> BigInt(shift);
  let z = BigInt(y.s) * y.m + BigInt(x.s) * a;
  if (!z) return real_0_bit(y.e - 63);
  const s = z < 0n ? -1 : 1;
  if (z < 0n) z = -z;
  const bits = z.toString(2).length;
  return {
    s,
    e: y.e + bits - 64,
    m: bits > 64 ? z >> BigInt(bits - 64) : z << BigInt(64 - bits),
    p: 64,
  };
}
