/** GMP 6.3.0 mpz/remove.c and mpn/generic/remove.c using native BigInt.
 * Copyright 1998-2002, 2009, 2012-2014, 2017 Free Software Foundation, Inc.
 * Adapted under GPL-2.0-or-later, matching this package.
 * Source: https://ftp.gnu.org/gnu/gmp/gmp-6.3.0.tar.xz
 * @see Deviation: Valuation dispatch and native GMP factor removal
 * Internal positive-factor domain: factor >= 2. Return [multiplicity, unit].
 */
import { ValueError } from '../errors.js';
export function mpz_remove(value: bigint, factor: bigint): [bigint, bigint] {
  if (factor < 2n)
    throw new ValueError(
      'You can only compute the valuation with respect to a integer larger than 1.'
    );
  if (value === 0n) return [0n, 0n];
  if (factor === 2n) {
    const exponent = BigInt((value & -value).toString(2).length - 1);
    return [exponent, value >> exponent];
  }
  let unit = value,
    divisor = factor;
  const powers: bigint[] = [];
  // Successively remove factor^(1), factor^(2), factor^(4), ... .
  while ((unit < 0n ? -unit : unit) >= divisor) {
    if (unit % divisor !== 0n) break;
    unit /= divisor;
    powers.push(divisor);
    divisor *= divisor;
  }
  let exponent = (1n << BigInt(powers.length)) - 1n;
  // Resolve the remaining multiplicity from high powers down to factor.
  for (let i = powers.length - 1; i >= 0; i--) {
    const power = powers[i]!;
    if (unit % power === 0n) {
      unit /= power;
      exponent += 1n << BigInt(i);
    }
  }
  return [exponent, unit];
}
