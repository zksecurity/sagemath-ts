/** PARI odd-extension elliptic point kernels.
 * Copyright (C) The PARI group; adapted under GPL-2.0-or-later.
 * Source: FpE.c:1519–1686.
 * Inputs are reduced ascending coefficient arrays over a valid odd finite field.
 * @see Deviation: PARI odd-extension elliptic kernels
 */
import { oddElliptic, type FqEllipticPoint, type FqEllipticChange } from './_odd_elliptic.js';
export function FpXQE_add(
  P: FqEllipticPoint,
  Q: FqEllipticPoint,
  a: bigint[],
  T: bigint[],
  p: bigint
): FqEllipticPoint {
  return oddElliptic(0, T, p).add(P, Q, a);
}

export function FpXQE_dbl(
  P: FqEllipticPoint,
  a: bigint[],
  T: bigint[],
  p: bigint
): FqEllipticPoint {
  return oddElliptic(0, T, p).dbl(P, a);
}

export function FpXQE_neg(P: FqEllipticPoint, T: bigint[], p: bigint): FqEllipticPoint {
  return oddElliptic(0, T, p).neg(P);
}

export function FpXQE_sub(
  P: FqEllipticPoint,
  Q: FqEllipticPoint,
  a: bigint[],
  T: bigint[],
  p: bigint
): FqEllipticPoint {
  return oddElliptic(0, T, p).sub(P, Q, a);
}

export function FpXQE_mul(
  P: FqEllipticPoint,
  n: bigint,
  a: bigint[],
  T: bigint[],
  p: bigint
): FqEllipticPoint {
  return oddElliptic(0, T, p).mul(P, n, a);
}

export function FpXQE_changepoint(
  P: FqEllipticPoint,
  ch: FqEllipticChange,
  T: bigint[],
  p: bigint
): FqEllipticPoint {
  return oddElliptic(0, T, p).change(P, ch);
}

export function FpXQE_changepointinv(
  P: FqEllipticPoint,
  ch: FqEllipticChange,
  T: bigint[],
  p: bigint
): FqEllipticPoint {
  return oddElliptic(0, T, p).changeinv(P, ch);
}
