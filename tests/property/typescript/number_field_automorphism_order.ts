import {
  NumberField,
  RationalPolynomial,
} from '../../../packages/sagemath-ts/src/rings/number_field/number_field.js';
import { Rational } from '../../../packages/sagemath-ts/src/rings/rational.js';

export function nf_automorphism_order(coefficients: bigint[], denominator: bigint): string {
  const K = new NumberField(
    new RationalPolynomial(coefficients.map((c) => new Rational(c, denominator))),
    'a'
  );
  return JSON.stringify(
    K.automorphisms().map((sigma) => sigma.__call__(K.gen()).list().map(String))
  );
}
