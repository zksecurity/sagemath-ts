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

import { gen_order, gen_powu_i, type GroupOrder } from './bb_group.js';
import { ellcard } from './elliptic/group.js';
/** FpE.c:1727: exact order from a supplied annihilating multiple.
 * @see Deviation: PARI generic and extension-curve order adapters
 */
export function FpXQE_order(P: FqEllipticPoint, order: GroupOrder, a: bigint[], T: bigint[], p: bigint): bigint {
  return gen_order(P, order, (Q, n) => FpXQE_mul(Q, n, a, T, p), Q => Q.isInfinity);
}


/** FpE.c:2046: Frobenius trace over an extension of nonnegative degree n.
 * @see Deviation: PARI base-field extension cardinality adapters
 */
export function elltrace_extension(t: bigint, n: number, q: bigint): bigint {
  // RgXQ_powu(X, n, X^2-tX+q), represented by its two exact coefficients.
  type Pair = readonly [bigint, bigint];
  const multiply = ([a, b]: Pair, [c, d]: Pair): Pair => {
    const bd = b*d;
    return [a*c-q*bd, a*d+b*c+t*bd];
  };
  const v = n===0 ? [1n, 0n] as const : gen_powu_i<Pair>([0n, 1n], BigInt(n),
    a=>multiply(a,a), multiply);
  return 2n*v[0]+t*v[1];
}
/** FpE.c:2055: count a base-field curve over the specified extension.
 * @see Deviation: PARI base-field extension cardinality adapters
 */
export function Fp_ffellcard(a4: bigint, a6: bigint, q: bigint, n: number, p: bigint): bigint {
  const trace = p+1n-ellcard({a4,a6,p});
  return q+1n-elltrace_extension(trace,n,p);
}

import { constantJCard } from './_extension_elliptic_cardinality.js';
/** Native constant-j extension count for p>3; q=p^n, n=degree(T).
 * @see Deviation: PARI constant-j extension counting adapters
 */
export function FpXQ_ellcardj(a4: bigint[], a6: bigint[], j: bigint, T: bigint[], q: bigint, p: bigint, n: number): bigint {
  return constantJCard(0,a4,a6,j,T,q,p,n);
}
