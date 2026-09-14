import {
  NumberField,
  NumberFieldElement,
  RationalPolynomial,
} from '../../../packages/sagemath-ts/src/rings/number_field/number_field.js';
import { Rational } from '../../../packages/sagemath-ts/src/rings/rational.js';
import { Integer } from '../../../packages/sagemath-ts/src/rings/integer_ring.js';
export function nf_element(
  fn: bigint,
  fcoeffs: bigint[],
  fden: bigint,
  coeffs: bigint[],
  den: bigint,
  exponent: bigint,
  kind: bigint
): string {
  const K = new NumberField(new RationalPolynomial(fcoeffs.map((c) => new Rational(c, fden))), 'a'),
    a = new NumberFieldElement(
      K,
      coeffs.map((c) => new Rational(c, den))
    );
  const poly = (p: RationalPolynomial) => p.coeffs.map(String),
    elem = (x: NumberFieldElement) => x.list().map(String);
  let r: unknown,
    error = null;
  try {
    if (fn === 0n) r = elem(a);
    else if (fn === 1n) r = elem(K.gen());
    else if (fn === 2n) r = poly(a.polynomial());
    else if (fn === 3n) r = poly(K.defining_polynomial());
    else if (fn === 4n) r = poly(K.polynomial());
    else if (fn === 5n) r = String(K);
    else if (fn === 6n) r = String(a);
    else if (fn === 7n) r = String(a.norm());
    else if (fn === 8n) r = String(a.trace());
    else if (fn === 9n) r = poly(a.charpoly());
    else if (fn === 10n) r = poly(a.minpoly());
    else if (fn === 11n) r = a.is_integral();
    else if (fn === 12n) r = a.is_unit();
    else if (fn === 13n) r = a.is_integral_unit();
    else if (fn === 14n) r = String(a.denominator());
    else if (fn === 15n) r = elem(a.numerator());
    else if (fn === 16n) {
      const result = a.pow(kind === 0n ? exponent : new Integer(exponent));
      r = [elem(result), result === a];
    } else if (fn === 17n) r = String(a.absolute_norm());
    else if (fn === 18n) r = String(a.relative_norm());
    else if (fn === 19n) r = String(a.absolute_trace());
    else if (fn === 20n) r = String(a.relative_trace());
    else if (fn === 21n) r = a.parent() === K;
    else if (fn === 22n) r = a.is_zero();
    else if (fn === 23n) r = a.is_one();
    else r = poly(a.minimal_polynomial());
  } catch (e) {
    r = null;
    error = e.name + ': ' + e.message;
  }
  return JSON.stringify([r, error]);
}
