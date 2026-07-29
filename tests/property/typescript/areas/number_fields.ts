/**
 * sagemath-ts side of the `number_fields` property-test area.
 *
 * Cases: tests/property/cases/number_fields.cases.json
 * SageMath counterpart: tests/property/python/areas/number_fields.py
 */

import { fundamental_discriminant } from '../../../../packages/sagemath-ts/src/arith/misc.js';
import {
  NumberField,
  QuadraticField,
  RationalPolynomial,
} from '../../../../packages/sagemath-ts/src/rings/number_field/number_field.js';
import { Rational } from '../../../../packages/sagemath-ts/src/rings/rational.js';

export const functions = {
  quadratic_discriminant: (d: bigint) => {
    const K = QuadraticField.create(d);
    return K.discriminant();
  },
  nf_degree: (coeffs: bigint[]) => {
    const poly = new RationalPolynomial(coeffs.map((c) => new Rational(c)));
    const K = new NumberField(poly, 'a');
    return BigInt(K.degree());
  },
  quadratic_signature: (d: bigint) => {
    const K = QuadraticField.create(d);
    const sig = K.signature();
    return sig;
  },
  is_totally_real: (d: bigint) => {
    const K = QuadraticField.create(d);
    return K.is_totally_real();
  },
  polynomial_discriminant: (coeffs: bigint[]) => {
    const poly = new RationalPolynomial(coeffs.map((c) => new Rational(c)));
    const disc = poly.discriminant();
    // Return as integer (numerator/denominator)
    return disc.numerator / disc.denominator;
  },
  fundamental_discriminant_test: (d: bigint) => fundamental_discriminant(d),
};
