/** Finite-field elliptic adapters from PARI basemath/ff.c.
 * Copyright (C) The PARI group; adapted under GPL-2.0-or-later.
 */
import { type PariFfelt } from './types.js';
import { EllCurveType } from './elliptic/init.js';
import { oddFFInit, oddFFMul, oddInitFq, oddFFOrder } from './_odd_elliptic_model.js';
import { type FlxqECoefficient, type FqEllipticChange } from './_odd_elliptic.js';
import { F2x_rem, F2xq_mul, F2xq_sqr, F2xq_inv, F2xq_div } from './F2x.js';
import {
  type F2xqECoefficient,
  type F2xqEChange,
  F2xqE_changepoint,
  F2xqE_changepointinv,
  F2xqE_mul,
  F2xqE_order,
} from './F2xqE.js';

export type FFEllipticScalar = bigint | PariFfelt;
export type FFEllipticInvariants = Readonly<
  Record<
    'a1' | 'a2' | 'a3' | 'a4' | 'a6' | 'b2' | 'b4' | 'b6' | 'b8' | 'c4' | 'c6' | 'disc',
    FFEllipticScalar
  >
>;
/** FF elements stay distinct from the packed F2x coefficients of the model. */
type FFInvariants = {
  readonly [K in keyof FFEllipticInvariants]: PariFfelt;
} & {
  readonly j: PariFfelt;
  readonly type: EllCurveType.t_ELL_Fq;
  readonly field: PariFfelt;
};
export type BinaryFFEllipticCurve = FFInvariants & {
  readonly binaryModel: readonly [F2xqECoefficient, bigint, F2xqEChange];
};
export type OddFFEllipticCurve = FFInvariants & {
  readonly oddModel: readonly [FlxqECoefficient, bigint[], FqEllipticChange];
};
export type FFEllipticCurve = BinaryFFEllipticCurve | OddFFEllipticCurve;
export type FFEllipticPoint =
  | { readonly isInfinity: true }
  | {
      readonly isInfinity: false;
      readonly x: PariFfelt;
      readonly y: PariFfelt;
    };
export type FFEllipticInputPoint =
  | { readonly isInfinity: true }
  | {
      readonly isInfinity: false;
      readonly x: FFEllipticScalar;
      readonly y: FFEllipticScalar;
    };

const pack = (a: readonly bigint[]) => a.reduce((v, c, i) => v | ((c & 1n) << BigInt(i)), 0n);
const invariantKeys = [
  'a1',
  'a2',
  'a3',
  'a4',
  'a6',
  'b2',
  'b4',
  'b6',
  'b8',
  'c4',
  'c6',
  'disc',
] as const;
function binaryModulus(fg: PariFfelt): bigint {
  return pack(fg.definingPoly ?? [0n, 1n]);
}
function coefficient(x: FFEllipticScalar, T: bigint): bigint {
  // A native integer is a constant, never a packed polynomial.
  if (typeof x === 'bigint') return x & 1n;
  return F2x_rem(typeof x.value === 'bigint' ? x.value & 1n : pack(x.value), T);
}
function element(x: bigint, fg: PariFfelt): PariFfelt {
  const value: bigint[] = [];
  while (x) {
    value.push(x & 1n);
    x >>= 1n;
  }
  return { ...fg, value };
}

/** ff.c:1200. Coefficients are already Rg_to_F2xq-converted here. */
function F2xq_ell_to_a4a6(
  [a1, a2, a3, a4, a6]: readonly [bigint, bigint, bigint, bigint, bigint],
  T: bigint
): BinaryFFEllipticCurve['binaryModel'] {
  const mul = (x: bigint, y: bigint) => F2xq_mul(x, y, T);
  const sqr = (x: bigint) => F2xq_sqr(x, T);
  if (a1 !== 0n) {
    const a1i = F2xq_inv(a1, T),
      a1i2 = sqr(a1i),
      a1i3 = mul(a1i, a1i2),
      a1i6 = sqr(a1i3);
    const d = mul(a3, a1i),
      dd = mul(d, a1i2);
    const e = mul(a4 ^ sqr(d), a1i),
      ee = mul(e, a1i3),
      da2 = a2 ^ d;
    return [mul(da2, a1i2), mul(mul(mul(da2, d) ^ a4, d) ^ a6 ^ sqr(e), a1i6), [a1i, dd, 0n, ee]];
  }
  // Native FF_ellinit also accepts singular curves (a1=a3=0).
  return [[a3, sqr(a2) ^ a4, a3 ? F2xq_inv(a3, T) : a3], mul(a2, a4) ^ a6, [1n, a2, 0n, 0n]];
}

/** ff.c:1363. Initialize an initsmall record; singular records have j=0.
 * Dispatches binary, word and arbitrary-prime polynomial backends.
 * @see Deviation: PARI binary elliptic model adapters
 * @see Deviation: PARI odd-extension elliptic model adapters
 */
export function FF_ellinit(E: FFEllipticInvariants, fg: PariFfelt): FFEllipticCurve {
  if (fg.p !== 2n) return oddFFInit(E, fg);
  const T = binaryModulus(fg);
  const binaryModel = F2xq_ell_to_a4a6(
    [
      coefficient(E.a1, T),
      coefficient(E.a2, T),
      coefficient(E.a3, T),
      coefficient(E.a4, T),
      coefficient(E.a6, T),
    ],
    T
  );
  const values = Object.fromEntries(
    invariantKeys.map((k) => [k, element(coefficient(E[k], T), fg)])
  ) as { [K in keyof FFEllipticInvariants]: PariFfelt };
  const c4 = coefficient(E.c4, T),
    D = coefficient(E.disc, T);
  return {
    ...values,
    type: EllCurveType.t_ELL_Fq,
    field: fg,
    j: element(D === 0n ? 0n : F2xq_div(F2xq_mul(F2xq_sqr(c4, T), c4, T), D, T), fg),
    binaryModel,
  };
}
/** ff.c:1466: convert, multiply on the native model, restore coordinates.
 * @see Deviation: PARI binary elliptic model adapters
 * @see Deviation: PARI odd-extension elliptic model adapters
 */
export function FF_ellmul(E: FFEllipticCurve, P: FFEllipticInputPoint, n: bigint): FFEllipticPoint {
  if ('oddModel' in E) return oddFFMul(E, P, n);
  const T = binaryModulus(E.field),
    [a, , ch] = E.binaryModel;
  const point = P.isInfinity
    ? P
    : {
        isInfinity: false as const,
        x: coefficient(P.x, T),
        y: coefficient(P.y, T),
      };
  const Q = F2xqE_changepoint(F2xqE_mul(F2xqE_changepointinv(point, ch, T), n, a, T), ch, T);
  return Q.isInfinity
    ? Q
    : { isInfinity: false, x: element(Q.x, E.field), y: element(Q.y, E.field) };
}
/** elliptic.c:798/521/6572: initialize from j, [a4,a6] or five coefficients.
 * Coefficient arithmetic dispatches by characteristic and native word size.
 * Integer inputs denote prime-field constants.
 * @see Deviation: PARI binary elliptic model adapters
 * @see Deviation: PARI odd-extension elliptic model adapters
 */
export function ellinit_Fq(
  x:
    | readonly [FFEllipticScalar]
    | readonly [FFEllipticScalar, FFEllipticScalar]
    | readonly [
        FFEllipticScalar,
        FFEllipticScalar,
        FFEllipticScalar,
        FFEllipticScalar,
        FFEllipticScalar,
      ],
  fg: PariFfelt
): FFEllipticCurve | null {
  if (x.length === 2) x = [0n, 0n, 0n, x[0], x[1]];
  if (fg.p !== 2n) return oddInitFq(x, fg);
  const T = binaryModulus(fg);
  if (x.length === 1) {
    const j = coefficient(x[0], T);
    x = j === 0n ? [0n, 0n, 1n, 0n, 0n] : [1n, 0n, 0n, 0n, element(F2xq_inv(j, T), fg)];
  }
  const [a1, a2, a3, a4, a6] = x.map((c) => coefficient(c, T));
  const mul = (a: bigint, b: bigint) => F2xq_mul(a, b, T),
    sqr = (a: bigint) => F2xq_sqr(a, T);
  const b2 = sqr(a1!),
    b4 = mul(a1!, a3!),
    b6 = sqr(a3!);
  const b8 = mul(b2, a6!) ^ mul(b6, a2!) ^ mul(a4!, a4! ^ b4);
  const c4 = sqr(b2),
    c6 = mul(b2, c4),
    disc = mul(b4, mul(b2, b6)) ^ mul(c4, b8) ^ sqr(b6);
  if (disc === 0n) return null;
  const values = { a1: a1!, a2: a2!, a3: a3!, a4: a4!, a6: a6!, b2, b4, b6, b8, c4, c6, disc };
  const E = Object.fromEntries(
    Object.entries(values).map(([k, v]) => [k, element(v, fg)])
  ) as FFEllipticInvariants;
  return FF_ellinit(E, fg);
}

import { type GroupOrder } from './bb_group.js';
/** ff.c:1495: convert to the native model before computing the supplied-bound order.
 * @see Deviation: PARI generic and extension-curve order adapters
 */
export function FF_ellorder(E: FFEllipticCurve, P: FFEllipticInputPoint, order: GroupOrder): bigint {
  if ('oddModel' in E) return oddFFOrder(E, P, order);
  const T = binaryModulus(E.field), [a, , ch] = E.binaryModel;
  const point = P.isInfinity ? P : {isInfinity: false as const, x: coefficient(P.x, T), y: coefficient(P.y, T)};
  return F2xqE_order(F2xqE_changepointinv(point, ch, T), order, a, T);
}
