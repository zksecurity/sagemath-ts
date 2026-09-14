/**
 * FLINT's continued-fraction Dedekind sum.
 * Reference: reference/flint/src/fmpq/dedekind_sum.c:fmpq_dedekind_sum.
 * @see Deviation: Dedekind sum backend arithmetic
 */
export function fmpq_dedekind_sum(h: bigint, k: bigint): [bigint, bigint] {
  if (k <= 2n || h === 0n) return [0n, 1n];
  let a = k;
  let b = ((h % k) + k) % k;
  let m11 = 1n;
  let m12 = 0n;
  let m21 = 0n;
  let m22 = 1n;
  let altSum = 0n;
  let determinant = 1n;
  // Both native branches accumulate this matrix and alternating quotient sum.
  // BigInt performs their arithmetic without machine-word overflow/limb packing.
  while (b !== 0n) {
    const q = a / b;
    [a, b] = [b, a % b];
    altSum += determinant * q;
    [m11, m12, m21, m22] = [m12 + q * m11, m11, m22 + q * m21, m21];
    determinant = -determinant;
  }
  if (determinant === -1n) altSum -= 3n;
  let numerator = altSum * m11 + m21 - determinant * m12;
  let denominator = 12n * m11;
  let g = numerator < 0n ? -numerator : numerator;
  let d = denominator;
  while (d !== 0n) [g, d] = [d, g % d];
  numerator /= g;
  denominator /= g;
  return [numerator, denominator];
}
