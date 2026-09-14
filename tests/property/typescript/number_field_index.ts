/** Native coefficient indexing, including quadratic representation selection. */
import {
  NumberField,
  RationalPolynomial,
} from '../../../packages/sagemath-ts/src/rings/number_field/number_field.js';
import { Rational } from '../../../packages/sagemath-ts/src/rings/rational.js';
import { Integer } from '../../../packages/sagemath-ts/src/rings/integer_ring.js';
export function nf_index(
  coeffs: bigint[],
  denominator: bigint,
  value: bigint,
  n: bigint,
  den: bigint,
  kind: bigint
): string {
  const K = new NumberField(
    new RationalPolynomial(coeffs.map((c) => new Rational(c, denominator))),
    'a'
  );
  const a =
    value === 0n
      ? K.zero()
      : value === 1n
        ? K.__call__(new Rational(2n, 3n))
        : K.gen().add(new Rational(2n, 3n));
  const i = kind === 0n ? Number(n) / Number(den) : kind === 1n ? n : new Integer(n);
  return String(a.__getitem__(i));
}
