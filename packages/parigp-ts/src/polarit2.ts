import { dbltor, rtodbl } from './kernel/none/mp_indep.js';
import { real_0_bit, type MpReal } from './qfb.js';
import { det } from './alglin1.js';

/**
 * Inexact PARI resultant (constant-first binary64 coefficient arrays).
 * @see Deviation: Polynomial Resultant Delegation and Real Double Coefficients
 */
export function resultant(a: readonly number[], b: readonly number[]): number {
  // Conversion happens before resultant_fast/init_resultant can inspect zero inputs.
  // cypari2 double_to_REAL deliberately uses 53-bit accuracy for a zero.
  const convert = (x: number): MpReal => (x === 0 ? real_0_bit(-53) : dbltor(x));
  const A = a.map(convert),
    B = b.map(convert);
  // Polrev/normalizepol_lg retain inexact leading zeros but clear the sign
  // when every coefficient is zero. Sage normalizes its own arrays earlier.
  if (!A.some((c) => c.s) || !B.some((c) => c.s)) return 0;
  const m = A.length - 1,
    n = B.length - 1,
    size = m + n;
  const matrix = Array.from({ length: size }, () => Array<MpReal | null>(size).fill(null));
  // RgX_sylvestermatrix: shifted descending coefficients are columns.
  for (let col = 0; col < n; col++) for (let j = 0; j <= m; j++) matrix[col + j]![col] = A[m - j]!;
  for (let col = 0; col < m; col++)
    for (let j = 0; j <= n; j++) matrix[col + j]![n + col] = B[n - j]!;
  const r = det(matrix);
  return r === null ? 0 : rtodbl(r);
}


import { roundr_safe } from './gen3.js';
import { PariError } from './errors.js';

/** Native integer/real rescaling to an integer column matrix.
 * @see Deviation: PARI integer and real matrix products and rescaling
 */
export function RgM_rescale_to_int(x: (bigint | MpReal)[][]): bigint[][] {
  if (x.some((c) => c.length !== x[0]!.length))
    throw new RangeError('RgM_rescale_to_int requires rectangular columns');
  let exact = true,
    emin = Infinity;
  for (const c of x)
    for (const v of c) {
      if (typeof v === 'bigint') {
        if (v) emin = Math.min(emin, (v < 0n ? -v : v).toString(2).length - 1);
      } else {
        exact = false;
        if (!v.s) continue;
        // Native strips zero low mantissa words then adds vals(lowest word).
        let m = v.m,
          e = v.e + 1 - v.p;
        while ((m & ((1n << 64n) - 1n)) === 0n) {
          m >>= 64n;
          e += 64;
        }
        while (!(m & 1n)) {
          m >>= 1n;
          e++;
        }
        emin = Math.min(emin, e);
      }
    }
  if (exact) return x.map((c) => c.slice()) as bigint[][];
  if (emin === Infinity) {
    // Native shifts every inexact zero by -HIGHEXPOBIT. Negative accuracy
    // then crosses the representable exponent boundary before grndtoi.
    if (x.some((c) => c.some((v) => typeof v !== 'bigint' && v.e < 0)))
      throw new PariError('overflow in expo()');
    return x.map((c) => c.map(() => 0n));
  }
  return x.map((c) =>
    c.map((v) => {
      if (typeof v !== 'bigint') return !v.s ? 0n : roundr_safe({ ...v, e: v.e - emin });
      if (!v) return 0n;
      if (emin <= 0) return v << BigInt(-emin);
      const shift = BigInt(emin);
      return (v + (1n << (shift - 1n))) >> shift;
    })
  );
}

/** Native signed remainder and optional centered half-modulus comparison.
 * @see Deviation: PARI bounded factor recombination adapters
 */
export function centermodii(x: bigint, p: bigint, half: bigint | null): bigint {
  if (!p) throw new PariError('impossible inverse in dvmdii: 0');
  let y = x % p;
  const beyond = half !== null && (y < 0n ? -y : y) > (half < 0n ? -half : half);
  if (y > 0n && beyond) y -= p;
  else if (y < 0n && (half === null || beyond)) y += p;
  return y;
}

import {
  rationalPair,
  rationalProduct,
  rationalDifference,
  rationalQuotient,
  type RationalPair,
} from './_rational_polynomial.js';
import { trimPolynomial } from './_polynomial_packing.js';
import { inverseCoefficient, residue } from './_polynomial_division.js';
/** Native integer-polynomial Newton sums and cache continuation.
 * A null modulus returns normalized rational pairs; modular output keeps native
 * negative residues. The extension-modulus parameter is null in this adapter.
 * @see Deviation: PARI integer polynomial Newton sum adapters
 */
export function polsym_gen(
  P: bigint[],
  y0: bigint[] | null,
  n: number,
  T: null,
  N: bigint
): bigint[];
export function polsym_gen(
  P: bigint[],
  y0: RationalPair[] | null,
  n: number,
  T: null,
  N: null
): RationalPair[];
export function polsym_gen(
  P: bigint[],
  y0: bigint[] | RationalPair[] | null,
  n: number,
  T: null,
  N: bigint | null
): bigint[] | RationalPair[] {
  if (n < 0) throw new PariError('sorry, polsym of a negative n is not yet implemented');
  P = trimPolynomial(P);
  if (!P.length) throw new PariError('zero polynomial in polsym');
  if (y0 && y0.length > n + 1)
    throw new RangeError('polsym_gen prefix exceeds the requested output length');
  const degree = P.length - 1,
    lead = P[degree]!,
    m = y0 === null ? 1 : y0.length;
  if (N !== null) {
    const y = y0 === null ? [BigInt(degree)] : (y0 as bigint[]).slice();
    let inv: bigint | null = null;
    if (lead !== 1n) {
      try {
        inv = inverseCoefficient(lead, N, false);
      } catch (e) {
        if (e instanceof PariError) throw new PariError(e.message.replace(/\.$/, ''));
        throw e;
      }
    }
    for (let k = m; k <= n; k++) {
      let s = k <= degree ? BigInt(k) * P[degree - k]! : 0n;
      for (let i = 1; i < k && i <= degree; i++) s += y[k - i]! * P[degree - i]!;
      s = residue(s, N);
      if (inv !== null) s = residue(s * inv, N);
      y[k] = -s;
    }
    return y;
  }
  const y: RationalPair[] =
    y0 === null
      ? [[BigInt(degree), 1n]]
      : (y0 as RationalPair[]).map((c) => rationalPair(c[0], c[1]));
  for (let k = m; k <= n; k++) {
    let s: RationalPair = [k <= degree ? BigInt(k) * P[degree - k]! : 0n, 1n];
    for (let i = 1; i < k && i <= degree; i++)
      s = rationalDifference(s, rationalProduct(y[k - i]!, [-P[degree - i]!, 1n]));
    if (lead !== 1n) s = rationalQuotient(s, [lead, 1n]);
    y[k] = [-s[0], s[1]];
  }
  return y;
}
/** Native uncached characteristic-zero Newton sums, as rational coefficient pairs.
 * @see Deviation: PARI integer polynomial Newton sum adapters
 */
export function polsym(P: bigint[], n: number): RationalPair[] {
  return polsym_gen(P, null, n, null, null);
}
