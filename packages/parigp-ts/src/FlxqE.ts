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

import { gen_order, type GroupOrder } from './bb_group.js';
/** FlxqE.c:313: exact order from a supplied annihilating multiple.
 * @see Deviation: PARI generic and extension-curve order adapters
 */
export function FlxqE_order(P: FqEllipticPoint, order: GroupOrder, a: FlxqECoefficient, T: bigint[], p: bigint): bigint {
  return gen_order(P, order, (Q, n) => FlxqE_mul(Q, n, a, T, p), Q => Q.isInfinity);
}

import { extensionField } from './_extension_field.js';
import { Flxq_issquare, Flxq_sqrt, Flxq_trace } from './Flx.js';
import { gen_pow_i } from './bb_group.js';
/** Supersingular characteristic-three count; FlxqE.c:1348–1375.
 * q=3^n, n=degree(T), and a4 is nonzero. Corrects the native odd-degree
 * twist-sign omission for a randomly selected nonsquare sqrt(-a4).
 * @see Deviation: PARI characteristic-three supersingular count correction
 */
export function F3xq_ellcardj(a4: bigint[], a6: bigint[], T: bigint[], q: bigint, n: number): bigint {
  const F = extensionField(1, T, 3n), q1 = q + 1n;
  const minusA = F.neg(a4) as bigint[];
  if (!Flxq_issquare(minusA, T, 3n)) return q1;
  const root = Flxq_sqrt(minusA, T, 3n)!;
  const denominator = F.mul(minusA, root);
  const c = F.mul(a6, F.inv(denominator)) as bigint[];
  const t = Flxq_trace(c, T, 3n);
  if (n % 2 === 1) {
    if (t === 0n) return q1;
    const q3 = 3n ** BigInt((n + 1) / 2);
    let trace = (t === 1n) !== (n % 4 === 1) ? q3 : -q3;
    // x = root*X introduces the quadratic twist root^3. Its character is
    // chi(root); replacing root by -root must not change the curve's count.
    if (!Flxq_issquare(root, T, 3n)) trace = -trace;
    return q1 - trace;
  }
  const q2 = 3n ** BigInt(n / 2);
  const W = gen_pow_i(a4, q >> 2n, a => F.sqr(a) as bigint[],
    (a, b) => F.mul(a, b) as bigint[]);
  const sign = (W[0] === 1n) !== (n % 4 === 2);
  if (t !== 0n) return sign ? q1 + q2 : q1 - q2;
  return sign ? q1 - 2n * q2 : q1 + 2n * q2;
}

import { constantJCard } from './_extension_elliptic_cardinality.js';
/** Native constant-j extension count for p>3; q=p^n, n=degree(T).
 * @see Deviation: PARI constant-j extension counting adapters
 */
export function Flxq_ellcardj(a4: bigint[], a6: bigint[], j: bigint, T: bigint[], q: bigint, p: bigint, n: number): bigint {
  return constantJCard(1,a4,a6,j,T,q,p,n);
}
