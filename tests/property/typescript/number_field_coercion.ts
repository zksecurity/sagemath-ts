/** Shared scalar coercion and zero-division probes. */
import {
  NumberField,
  RationalPolynomial,
} from '../../../packages/sagemath-ts/src/rings/number_field/number_field.js';
import { Rational } from '../../../packages/sagemath-ts/src/rings/rational.js';
import { Integer } from '../../../packages/sagemath-ts/src/rings/integer_ring.js';
const field = (d: bigint) =>
  new NumberField(RationalPolynomial.fromBigInts(d === 0n ? [-2n, 0n, 0n, 1n] : [-d, 0n, 1n]), 'a');
export function nf_coercion(
  d: bigint,
  fn: bigint,
  c: bigint,
  n: bigint,
  den: bigint,
  kind: bigint
): string {
  const K = field(d),
    a = K.gen()
      .scalarMul(new Rational(c, 3n))
      .add(K.__call__(new Rational(2n, 3n)));
  const q = new Rational(n, den),
    scalar = kind === 0n ? n : kind === 1n ? new Integer(n) : kind === 2n ? q : K.__call__(q);
  let r: unknown;
  if (fn === 0n) r = a.add(scalar).list().map(String);
  else if (fn === 1n) r = a.sub(scalar).list().map(String);
  else if (fn === 2n) r = a.div(scalar).list().map(String);
  else if (fn === 3n) r = a.eq(scalar);
  else r = K.__call__(scalar).list().map(String);
  return JSON.stringify(r);
}
export function nf_zero(d: bigint, fn: bigint, kind: bigint): string {
  const K = field(d),
    z = K.zero(),
    a = K.gen();
  const scalar =
    kind === 0n ? 0n : kind === 1n ? new Integer(0n) : kind === 2n ? Rational.zero() : z;
  const r =
    fn === 0n ? z.inv() : fn === 1n ? a.div(scalar) : fn === 2n ? z.div(scalar) : z.pow(-1n);
  return JSON.stringify(r.list().map(String));
}
export function nf_float(d: bigint, n: bigint, den: bigint, scale: bigint): string {
  const v = (Number(n) / Number(den)) * 2 ** Number(scale);
  return JSON.stringify(field(d).__call__(v).list().map(String));
}
