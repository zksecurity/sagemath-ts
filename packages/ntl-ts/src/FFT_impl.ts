/** Bundled generic FFT admissible-size rounding.
 * @see Deviation: NTL portable truncated transforms
 */
export function FFTRoundUp(xn: number, k: number): number {
  if (!Number.isInteger(k) || k < 0 || k > 62)
    throw new RangeError('FFTRoundUp: exponent must be between 0 and 62');
  if (!Number.isInteger(xn) || xn < -(2 ** 63) || xn >= 2 ** 63)
    throw new RangeError('FFTRoundUp: length must be a signed native integer');
  const value = xn <= 0 ? 1n : BigInt(xn);
  if (value >= (1n << 63n) - 15n) throw new RangeError('FFTRoundUp: length rounding overflows');
  const n = 1n << BigInt(k),
    rounded = ((value + 15n) >> 4n) << 4n;
  return Number(rounded > n - (n >> BigInt(k >= 10 ? 4 : 3)) ? n : rounded);
}
