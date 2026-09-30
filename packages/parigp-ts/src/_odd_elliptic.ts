/** Shared formulas from PARI FlxqE.c and FpE.c (FpXQE).
 * Copyright (C) The PARI group; adapted under GPL-2.0-or-later.
 */
import { extensionField } from './_extension_field.js';
import { gen_pow_i } from './bb_group.js';

export type FqEllipticPoint =
  | { readonly isInfinity: true }
  | {
      readonly isInfinity: false;
      readonly x: bigint[];
      readonly y: bigint[];
    };
export type FqEllipticChange = readonly [bigint[], bigint[], bigint[], bigint[]];
/** Polynomial a4, or native characteristic-three ordinary vector [a2]. */
export type FlxqECoefficient = bigint[] | readonly [bigint[]];

export function oddElliptic(mode: 0 | 1, T: bigint[], p: bigint) {
  const F = extensionField(mode, T, p);
  const infinity = (): FqEllipticPoint => ({ isInfinity: true });
  const sum = (a: bigint[], b: bigint[]) => F.add(a, b) as bigint[];
  const difference = (a: bigint[], b: bigint[]) => F.sub(a, b) as bigint[];
  const product = (a: bigint[], b: bigint[]) => F.mul(a, b) as bigint[];
  const square = (a: bigint[]) => F.sqr(a) as bigint[];
  const inverse = (a: bigint[]) => F.inv(a) as bigint[];
  const quotient = (a: bigint[], b: bigint[]) => product(a, inverse(b));
  const scale = (a: bigint[], n: bigint) => F.mul(a, n) as bigint[];
  const minus = (a: bigint[]) => F.neg(a) as bigint[];
  const equal = (a: bigint[], b: bigint[]) =>
    a.length === b.length && a.every((v, i) => v === b[i]);
  const ordinary = (a: FlxqECoefficient): a is readonly [bigint[]] => Array.isArray(a[0]);
  const point = (x: bigint[], y: bigint[]): FqEllipticPoint => ({ isInfinity: false, x, y });

  // FlxqE.c:115 / FpE.c:1580. The word backend retains the F3 ordinary x² term.
  const dbl = (P: FqEllipticPoint, a: FlxqECoefficient): FqEllipticPoint => {
    if (P.isInfinity || P.y.length === 0) return infinity();
    const s =
      mode === 1 && p === 3n
        ? ordinary(a)
          ? quotient(product(P.x, a[0]), P.y)
          : quotient(a, minus(P.y))
        : quotient(sum(scale(square(P.x), 3n), a as bigint[]), scale(P.y, 2n));
    let x = difference(square(s), scale(P.x, 2n));
    if (mode === 1 && ordinary(a)) x = difference(x, a[0]);
    return point(x, difference(product(s, difference(P.x, x)), P.y));
  };
  const add = (P: FqEllipticPoint, Q: FqEllipticPoint, a: FlxqECoefficient): FqEllipticPoint => {
    if (P.isInfinity) return Q;
    if (Q.isInfinity) return P;
    if (equal(P.x, Q.x)) return equal(P.y, Q.y) ? dbl(P, a) : infinity();
    const s = quotient(difference(P.y, Q.y), difference(P.x, Q.x));
    let x = difference(difference(square(s), P.x), Q.x);
    if (mode === 1 && ordinary(a)) x = difference(x, a[0]);
    return point(x, difference(product(s, difference(P.x, x)), P.y));
  };
  const neg = (P: FqEllipticPoint): FqEllipticPoint =>
    P.isInfinity ? infinity() : point([...P.x], minus(P.y));
  const sub = (P: FqEllipticPoint, Q: FqEllipticPoint, a: FlxqECoefficient): FqEllipticPoint =>
    add(P, neg(Q), a);
  const mul = (P: FqEllipticPoint, n: bigint, a: FlxqECoefficient): FqEllipticPoint => {
    if (n === 0n || P.isInfinity) return infinity();
    if (n < 0n) P = neg(P);
    if (n === 1n || n === -1n) return P.isInfinity ? infinity() : point([...P.x], [...P.y]);
    return gen_pow_i(
      P,
      n,
      (Q) => dbl(Q, a),
      (Q, R) => add(Q, R, a)
    );
  };
  const change = (P: FqEllipticPoint, [u, r, s, t]: FqEllipticChange): FqEllipticPoint => {
    if (P.isInfinity) return P;
    const v = inverse(u),
      v2 = square(v),
      v3 = product(v, v2),
      c = difference(P.x, r);
    return point(product(v2, c), product(v3, difference(P.y, sum(product(s, c), t))));
  };
  const changeinv = (P: FqEllipticPoint, [u, r, s, t]: FqEllipticChange): FqEllipticPoint => {
    if (P.isInfinity) return P;
    const u2 = square(u),
      u3 = product(u, u2),
      c = product(u2, P.x);
    return point(sum(c, r), sum(product(u3, P.y), sum(product(s, c), t)));
  };
  return { add, dbl, neg, sub, mul, change, changeinv };
}
