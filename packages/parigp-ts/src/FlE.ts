/** Word-prime elliptic kernels, PARI FlE.c (GPL-2.0-or-later).
 * Word values use exact bigint residues. The precomputed reduction word pi is
 * retained in signatures; BigInt modular reduction does not need it.
 * @see Deviation: PARI word-prime elliptic adapters
 */
import { Fp_inv, Fp_neg } from './ff.js';
import { gen_order, gen_pow_i, type GroupOrder } from './bb_group.js';
import { FpJ_add, FpJ_dbl } from './elliptic/point.js';
import { ellinf, mkpoint, type EllipticPoint, type JacobianPoint } from './elliptic/points.js';
const red = (x: bigint, p: bigint) => ((x % p) + p) % p;

/** FlE.c:35–71; finite formulas are identical to the FpJ formulas. */
export function Flj_dbl_pre(P: JacobianPoint, a4: bigint, p: bigint, _pi: bigint): JacobianPoint {
  return P.Z === 0n ? { ...P } : FpJ_dbl(P, a4, p);
}
/** FlE.c:77–129. */
export function Flj_add_pre(P: JacobianPoint, Q: JacobianPoint, a4: bigint, p: bigint, _pi: bigint): JacobianPoint {
  return FpJ_add(P, Q, a4, p);
}
export function Flj_neg(P: JacobianPoint, p: bigint): JacobianPoint {
  return { X: P.X, Y: Fp_neg(P.Y, p), Z: P.Z };
}
/** FlE.c:136–189; left-to-right non-adjacent signed binary powering. */
export function Flj_mulu_pre(P: JacobianPoint, n: bigint, a4: bigint, p: bigint, pi: bigint): JacobianPoint {
  if (n === 0n) return { X: 1n, Y: 1n, Z: 0n };
  if (n === 1n) return { ...P };
  const digits: number[] = [];
  for (let k = n; k > 0n; k >>= 1n) {
    const digit = (k & 1n) === 0n ? 0 : Number(2n - (k & 3n));
    digits.push(digit); k -= BigInt(digit);
  }
  let R = Flj_dbl_pre(P, a4, p, pi);
  const minus = Flj_neg(P, p);
  for (let i = digits.length - 3; i >= 0; i--) {
    R = Flj_dbl_pre(R, a4, p, pi);
    if (digits[i]) R = Flj_add_pre(R, digits[i]! > 0 ? P : minus, a4, p, pi);
  }
  return R;
}
export function Fle_to_Flj(P: EllipticPoint): JacobianPoint {
  return P.isInfinity ? { X: 1n, Y: 1n, Z: 0n } : { X: P.x, Y: P.y, Z: 1n };
}
export function Flj_to_Fle_pre(P: JacobianPoint, p: bigint, _pi: bigint): EllipticPoint {
  if (P.Z === 0n) return ellinf();
  const z = Fp_inv(P.Z, p), z2 = z * z % p;
  return mkpoint(P.X * z2 % p, P.Y * z2 % p * z % p);
}
/** FlE.c:346–365. */
export function Fle_dbl(P: EllipticPoint, a4: bigint, p: bigint): EllipticPoint {
  if (P.isInfinity || P.y === 0n) return ellinf();
  const slope = red((3n * P.x * P.x + a4) * Fp_inv(2n * P.y % p, p), p);
  const x = red(slope * slope - 2n * P.x, p);
  return mkpoint(x, red(slope * (P.x - x) - P.y, p));
}
/** FlE.c:369–392. */
export function Fle_add(P: EllipticPoint, Q: EllipticPoint, a4: bigint, p: bigint): EllipticPoint {
  if (P.isInfinity) return Q;
  if (Q.isInfinity) return P;
  if (P.x === Q.x) return P.y === Q.y ? Fle_dbl(P, a4, p) : ellinf();
  const slope = red((P.y - Q.y) * Fp_inv(red(P.x - Q.x, p), p), p);
  const x = red(slope * slope - P.x - Q.x, p);
  return mkpoint(x, red(slope * (P.x - x) - P.y, p));
}
/** FlE.c:428–435. n is an unsigned 64-bit word. */
export function Fle_mulu(P: EllipticPoint, n: bigint, a4: bigint, p: bigint): EllipticPoint {
  if (n === 0n || P.isInfinity) return ellinf();
  if (n === 1n) return { ...P };
  if (n === 2n) return Fle_dbl(P, a4, p);
  return Flj_to_Fle_pre(Flj_mulu_pre(Fle_to_Flj(P), n, a4, p, 0n), p, 0n);
}
/** FlE.c:439–452,484–491; native generic-order recursion. */
export function Fle_order(P: EllipticPoint, order: GroupOrder, a4: bigint, p: bigint): bigint {
  const multiply = (Q: EllipticPoint, n: bigint): EllipticPoint => {
    if (n === 0n || Q.isInfinity) return ellinf();
    if (n < 0n) { Q = mkpoint(Q.x, Fp_neg(Q.y, p)); n = -n; }
    if (n < (1n << 64n)) return Fle_mulu(Q, n, a4, p);
    return gen_pow_i<EllipticPoint>(Q, n, R => Fle_dbl(R, a4, p), (R, S) => Fle_add(R, S, a4, p));
  };
  return gen_order(P, order, multiply, Q => Q.isInfinity);
}
