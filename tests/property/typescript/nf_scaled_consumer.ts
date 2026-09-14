import {
  NumberField,
  RationalPolynomial,
} from '../../../packages/sagemath-ts/src/rings/number_field/number_field.js';
import { Rational } from '../../../packages/sagemath-ts/src/rings/rational.js';
import { Frobenius_filter } from '../../../packages/sagemath-ts/src/schemes/elliptic_curves/isogeny_class.js';
export function nf_scaled_consumer(op: bigint, cs: bigint[], num: bigint, den: bigint): string {
  const K = new NumberField(new RationalPolynomial(cs.map((c) => new Rational(c * num, den))), 'a'),
    a = K.gen();
  let r: unknown,
    error = null;
  try {
    if (op === 0n) {
      const O = K.maximal_order() as import(
        '../../../packages/sagemath-ts/src/rings/number_field/order.js'
      ).AbsoluteOrder;
      r = [String(O.different().norm()), String(O.codifferent().norm())];
    } else {
      const E = { ainvs: () => [K.zero(), K.zero(), K.zero(), a.add(1n), a], base_field: () => K };
      r = Frobenius_filter(E as never, [2n, 3n, 5n, 7n, 11n, 13n]).map(String);
    }
  } catch (e) {
    r = null;
    error = e.name + ': ' + e.message;
  }
  return JSON.stringify([r, error]);
}
