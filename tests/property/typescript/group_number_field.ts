/** Number-field actions exercised through the generic group API. */
import * as g from '../../../packages/sagemath-ts/src/groups/generic.js';
import {
  NumberField,
  QuadraticField,
  RationalPolynomial,
} from '../../../packages/sagemath-ts/src/rings/number_field/number_field.js';
import { Rational } from '../../../packages/sagemath-ts/src/rings/rational.js';
import { Integer } from '../../../packages/sagemath-ts/src/rings/integer_ring.js';
const field = (d: bigint) =>
  d === 0n
    ? new NumberField(RationalPolynomial.fromBigInts([-2n, 0n, 0n, 1n]), 'a')
    : QuadraticField.create(d, 'a');
export function nf_scalar(d: bigint, c: bigint, n: bigint, den: bigint, kind: bigint): string {
  const K = field(d),
    a = K.gen().add(K.__call__(c)).scalarMul(new Rational(1n, 3n));
  const scalar = kind === 4n ? new Rational(n, den) : kind === 3n ? new Integer(n) : n;
  return JSON.stringify(a.mul(scalar).list().map(String));
}
export function nf_group(d: bigint, fn: bigint, mode: bigint, n: bigint, variant: bigint): string {
  const K = field(d),
    op = mode === 0n ? '+' : '*';
  const a =
    variant === 0n ? K.gen() : variant === 1n ? K.__call__(-1n) : op === '+' ? K.zero() : K.one();
  const target = op === '+' ? a.scalarMul(new Rational(n)) : a.pow(n);
  let r: unknown;
  if (fn === 0n) r = g.multiple(a, n, op).list().map(String);
  else if (fn === 1n) r = String(g.bsgs(a, target, [0n, 64n], op));
  else if (fn === 2n) r = g.has_order(a, n, op);
  else if (fn === 3n) r = String(g.order_from_multiple(a, n, undefined, op));
  else if (fn === 4n) r = String(g.discrete_log(target, a, 4n, op));
  else if (fn === 5n) r = String(g.pohlig_hellman(target, a, 4n, undefined, op));
  else r = String(g.order_from_bounds(a, [1n, 16n], undefined, op));
  return JSON.stringify(r);
}
