/** Integer logarithms from PARI basemath/ispower.c.
 * @see Deviation: PARI Galois integer helper boundaries
 */
import { PariError } from './errors.js';

/** Checked integer entry point, preserving base validation before the operand. */
export function logint0(B: bigint, y: bigint): number {
  if (y < 2n) throw new PariError('domain error in logint: b <= 1');
  if (B <= 0n) throw new PariError('domain error in logint: x <= 0');
  return logintall(B, y);
}

/** Internal native contract: B > 0, y > 1; omit the optional power output. */
export function logintall(B: bigint, y: bigint): number {
  const word = 1n << 64n;
  const eB = B.toString(2).length - 1;
  if (B < word) {
    if (y >= word) return 0;
    if (y === 2n) return eB;
    let r = y;
    for (let e = 1; ; e++) {
      if (r >= B) return r === B ? e : e - 1;
      r *= y;
      if (r >= word) return e;
    }
  }
  if (y === 2n) return eB;
  const ey = y.toString(2).length - 1;
  const emax = Math.floor(eB / ey);
  if (emax <= 13) {
    let r = y;
    for (let e = 1; ; e++) {
      if (r >= B) return r === B ? e : e - 1;
      r *= y;
    }
  }
  // Binary splitting: build y^(2^i), then determine the exponent bit by bit.
  const powers = [y];
  let i = 0,
    q = y;
  for (;;) {
    const r = powers[i]!;
    if (r === B) return 2 ** i;
    if (r > B) {
      i--;
      break;
    }
    q = r;
    if (2 ** (i + 1) > emax) break;
    powers[++i] = q * q;
  }
  let e = 2 ** i;
  while (--i >= 0) {
    const r = q * powers[i]!;
    if (r <= B) {
      e += 2 ** i;
      q = r;
      if (r === B) break;
    }
  }
  return e;
}
