/** Integer valuation and finite-field universal comparison, PARI basemath/gen2.c. */

/** Native divide-and-conquer valuation for positive x and q>1. */
function Z_pvalrem_DC(x: bigint, q: bigint): [number, bigint] {
  if (x % q !== 0n) return [0, x];
  const z = x / q;
  const limbs = (n: bigint) => Math.ceil(n.toString(2).length / 64) + 2;
  let v = 0,
    unit = z;
  if (2 * limbs(q) <= limbs(z) + 3) {
    const reduced = Z_pvalrem_DC(z, q * q);
    v = 2 * reduced[0];
    unit = reduced[1];
  }
  return unit % q !== 0n ? [v + 1, unit] : [v + 2, unit / q];
}

/**
 * Return [valuation, unit] for nonzero n and p>1.
 * Follows Z_lvalrem's 16-division threshold and divide-and-conquer path
 * for word-size p; larger p follows Z_pvalrem's original division loop.
 * @see Reference: reference/pari/src/basemath/gen2.c:Z_lvalrem/Z_pvalrem
 * @see Deviation: Hilbert symbol dependency domain
 */
export function Z_pvalrem(n: bigint, p: bigint): [number, bigint] {
  if (n === 0n || p <= 1n) throw new RangeError('Z_pvalrem requires nonzero n and p > 1');
  if (p === 2n) {
    const exponent = (n & -n).toString(2).length - 1;
    return [exponent, n >> BigInt(exponent)];
  }
  const negative = n < 0n;
  let unit = negative ? -n : n,
    v = 0;
  const large = unit >= 1n << 64n,
    word = p < 1n << 64n;
  while (unit % p === 0n) {
    unit /= p;
    v++;
    if (large && word && v === 16) {
      const reduced = Z_pvalrem_DC(unit, p * p);
      v += 2 * reduced[0];
      unit = reduced[1];
      if (unit % p === 0n) {
        v++;
        unit /= p;
      }
      break;
    }
  }
  return [v, negative ? -unit : unit];
}

import type { PariFfelt } from './types.js';

/**
 * PARI cmp_universal for reduced finite-field elements at one variable index.
 * Inputs use PariFfelt coefficient arrays and defining polynomials; a missing
 * defining polynomial denotes the prime-field modulus X. Other GEN types are
 * outside this adapter. Values must already be reduced in their own field.
 * @see Reference: reference/pari/src/basemath/gen2.c:cmp_universal
 * @see Deviation: PARI finite-field universal comparison
 */
export function cmp_universal(x: PariFfelt, y: PariFfelt): number {
  // FF.c:ffgen selects FpXQ=0, Flxq=1, F2xq=2 on the 64-bit native profile.
  const tx = x.p === 2n ? 2 : x.p < 1n << 64n ? 1 : 0;
  const ty = y.p === 2n ? 2 : y.p < 1n << 64n ? 1 : 0;
  if (tx !== ty) return tx < ty ? -1 : 1;
  const value = cmp_FF_raw(
    typeof x.value === 'bigint' ? [x.value] : x.value,
    typeof y.value === 'bigint' ? [y.value] : y.value,
    tx
  );
  if (value) return value;
  const modulus = cmp_FF_raw(x.definingPoly ?? [0n, 1n], y.definingPoly ?? [0n, 1n], tx);
  return modulus || (x.p < y.p ? -1 : x.p > y.p ? 1 : 0);
}

/** Same-backend raw comparison: length, then native words/coefficients low first. */
function cmp_FF_raw(x: readonly bigint[], y: readonly bigint[], backend: number): number {
  let nx = x.length,
    ny = y.length;
  while (nx && x[nx - 1] === 0n) nx--;
  while (ny && y[ny - 1] === 0n) ny--;
  // F2x stores packed 64-bit words, not one word per coefficient.
  const lx = backend === 2 ? Math.ceil(nx / 64) : nx;
  const ly = backend === 2 ? Math.ceil(ny / 64) : ny;
  if (lx !== ly) return lx < ly ? -1 : 1;
  for (let i = 0; i < lx; i++) {
    let a = x[i]!,
      b = y[i]!;
    if (backend === 2) {
      a = 0n;
      b = 0n;
      for (let j = 0; j < 64; j++) {
        a |= (x[64 * i + j] ?? 0n) << BigInt(j);
        b |= (y[64 * i + j] ?? 0n) << BigInt(j);
      }
    }
    // t_VECSMALL compares signed long slots, including Flx coefficients above
    // 2^63 and F2x words whose highest bit is set. FpX slots are signed bigints.
    if (backend !== 0) {
      a = BigInt.asIntN(64, a);
      b = BigInt.asIntN(64, b);
    }
    if (a !== b) return a < b ? -1 : 1;
  }
  return 0;
}
