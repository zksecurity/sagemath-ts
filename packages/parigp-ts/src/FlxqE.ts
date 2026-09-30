/** PARI odd-extension elliptic point kernels.
 * Copyright (C) The PARI group; adapted under GPL-2.0-or-later.
 * Source: FlxqE.c:40–246.
 * Inputs are reduced ascending coefficient arrays over a valid odd finite field.
 * @see Deviation: PARI odd-extension elliptic kernels
 */
import {
  oddElliptic,
  type FqEllipticPoint,
  type FqEllipticChange,
  type FlxqECoefficient,
} from './_odd_elliptic.js';
export function FlxqE_add(
  P: FqEllipticPoint,
  Q: FqEllipticPoint,
  a: FlxqECoefficient,
  T: bigint[],
  p: bigint
): FqEllipticPoint {
  return oddElliptic(1, T, p).add(P, Q, a);
}

export function FlxqE_dbl(
  P: FqEllipticPoint,
  a: FlxqECoefficient,
  T: bigint[],
  p: bigint
): FqEllipticPoint {
  return oddElliptic(1, T, p).dbl(P, a);
}

export function FlxqE_neg(P: FqEllipticPoint, T: bigint[], p: bigint): FqEllipticPoint {
  return oddElliptic(1, T, p).neg(P);
}

export function FlxqE_sub(
  P: FqEllipticPoint,
  Q: FqEllipticPoint,
  a: FlxqECoefficient,
  T: bigint[],
  p: bigint
): FqEllipticPoint {
  return oddElliptic(1, T, p).sub(P, Q, a);
}

export function FlxqE_mul(
  P: FqEllipticPoint,
  n: bigint,
  a: FlxqECoefficient,
  T: bigint[],
  p: bigint
): FqEllipticPoint {
  return oddElliptic(1, T, p).mul(P, n, a);
}

export function FlxqE_changepoint(
  P: FqEllipticPoint,
  ch: FqEllipticChange,
  T: bigint[],
  p: bigint
): FqEllipticPoint {
  return oddElliptic(1, T, p).change(P, ch);
}

export function FlxqE_changepointinv(
  P: FqEllipticPoint,
  ch: FqEllipticChange,
  T: bigint[],
  p: bigint
): FqEllipticPoint {
  return oddElliptic(1, T, p).changeinv(P, ch);
}
