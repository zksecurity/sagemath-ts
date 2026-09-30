/** PARI basemath/F2xqE.c: binary elliptic point arithmetic.
 * Copyright (C) 2012 The PARI group; adapted under GPL-2.0-or-later.
 * Packed nonnegative bigint coefficient bits use the existing F2x convention.
 * @see Deviation: PARI binary elliptic point kernels
 */
import { F2xq_mul, F2xq_sqr, F2xq_inv, F2xq_div } from './F2x.js';
import { gen_pow_i } from './bb_group.js';
import { type EllipticPoint, ellinf } from './elliptic/points.js';

/** Ordinary a2, or supersingular [a3,a4,inverse(a3)]. */
export type F2xqECoefficient = bigint | readonly [bigint, bigint, bigint];
export type F2xqEChange = readonly [bigint, bigint, bigint, bigint];

export function F2xqE_changepoint(P: EllipticPoint, ch: F2xqEChange, T: bigint): EllipticPoint {
  if (P.isInfinity) return P;
  const [u, r, s, t] = ch;
  const v = F2xq_inv(u, T),
    v2 = F2xq_sqr(v, T),
    v3 = F2xq_mul(v, v2, T);
  const c = P.x ^ r;
  return {
    isInfinity: false,
    x: F2xq_mul(v2, c, T),
    y: F2xq_mul(v3, P.y ^ F2xq_mul(s, c, T) ^ t, T),
  };
}
export function F2xqE_changepointinv(P: EllipticPoint, ch: F2xqEChange, T: bigint): EllipticPoint {
  if (P.isInfinity) return P;
  const [u, r, s, t] = ch;
  const u2 = F2xq_sqr(u, T),
    u3 = F2xq_mul(u, u2, T),
    c = F2xq_mul(u2, P.x, T);
  return { isInfinity: false, x: c ^ r, y: F2xq_mul(u3, P.y, T) ^ F2xq_mul(s, c, T) ^ t };
}
export function F2xqE_dbl(P: EllipticPoint, a: F2xqECoefficient, T: bigint): EllipticPoint {
  if (P.isInfinity) return ellinf();
  if (typeof a === 'bigint') {
    if (P.x === 0n) return ellinf();
    const slope = P.x ^ F2xq_div(P.y, P.x, T);
    const x = F2xq_sqr(slope, T) ^ slope ^ a;
    return { isInfinity: false, x, y: F2xq_mul(slope, P.x ^ x, T) ^ P.y ^ x };
  }
  const [a3, a4, a3i] = a;
  const slope = F2xq_mul(a4 ^ F2xq_sqr(P.x, T), a3i, T),
    x = F2xq_sqr(slope, T);
  return { isInfinity: false, x, y: F2xq_mul(slope, P.x ^ x, T) ^ P.y ^ a3 };
}
export function F2xqE_add(
  P: EllipticPoint,
  Q: EllipticPoint,
  a: F2xqECoefficient,
  T: bigint
): EllipticPoint {
  if (P.isInfinity) return Q;
  if (Q.isInfinity) return P;
  if (P.x === Q.x) return P.y === Q.y ? F2xqE_dbl(P, a, T) : ellinf();
  const slope = F2xq_div(P.y ^ Q.y, P.x ^ Q.x, T);
  const x = F2xq_sqr(slope, T) ^ P.x ^ Q.x ^ (typeof a === 'bigint' ? slope ^ a : 0n);
  return {
    isInfinity: false,
    x,
    y: F2xq_mul(slope, P.x ^ x, T) ^ P.y ^ (typeof a === 'bigint' ? x : a[0]),
  };
}
export function F2xqE_neg(P: EllipticPoint, a: F2xqECoefficient, T: bigint): EllipticPoint {
  if (P.isInfinity) return ellinf();
  return { isInfinity: false, x: P.x, y: P.y ^ (typeof a === 'bigint' ? P.x : a[0]) };
}
export function F2xqE_sub(
  P: EllipticPoint,
  Q: EllipticPoint,
  a: F2xqECoefficient,
  T: bigint
): EllipticPoint {
  return F2xqE_add(P, F2xqE_neg(Q, a, T), a, T);
}
export function F2xqE_mul(
  P: EllipticPoint,
  n: bigint,
  a: F2xqECoefficient,
  T: bigint
): EllipticPoint {
  if (n === 0n || P.isInfinity) return ellinf();
  if (n < 0n) P = F2xqE_neg(P, a, T);
  if (n === 1n || n === -1n) return { ...P };
  return gen_pow_i(
    P,
    n,
    (Q) => F2xqE_dbl(Q, a, T),
    (Q, R) => F2xqE_add(Q, R, a, T)
  );
}

import { gen_order, type GroupOrder } from './bb_group.js';
/** F2xqE.c:298: exact order from a supplied annihilating multiple.
 * @see Deviation: PARI generic and extension-curve order adapters
 */
export function F2xqE_order(P: EllipticPoint, order: GroupOrder, a: F2xqECoefficient, T: bigint): bigint {
  return gen_order(P, order, (Q, n) => F2xqE_mul(Q, n, a, T), Q => Q.isInfinity);
}
