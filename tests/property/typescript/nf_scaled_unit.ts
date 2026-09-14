import {
  NumberField,
  RationalPolynomial,
} from '../../../packages/sagemath-ts/src/rings/number_field/number_field.js';
import { Rational } from '../../../packages/sagemath-ts/src/rings/rational.js';
export function nf_scaled_unit(
  op: bigint,
  cs: bigint[],
  num: bigint,
  den: bigint,
  exponent: bigint
): string {
  const K = new NumberField(new RationalPolynomial(cs.map((c) => new Rational(c * num, den))), 'a');
  let r: unknown,
    error = null;
  try {
    const U = K.unit_group();
    if (op === 0n)
      r = K.signature()[0]
        ? U.fundamental_units()[0]!.list().map(String)
        : U.roots_of_unity()
            .map((z) => z.list().map(String).join(','))
            .sort();
    else if (op === 1n) {
      const v = [1n, ...(K.signature()[0] ? [exponent] : [])];
      r = U.log(U.exp(v)).map(String);
    } else if (op === 2n) r = String(Math.round(K.regulator() * 10 ** 8));
    else r = String(U.torsion_order());
  } catch (e) {
    r = null;
    error = e.name + ': ' + e.message;
  }
  return JSON.stringify([r, error]);
}
