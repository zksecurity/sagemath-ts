import { Integer } from '../../../packages/sagemath-ts/src/rings/integer_ring.js';
import {
  NumberField,
  NumberFieldElement,
  RationalPolynomial,
} from '../../../packages/sagemath-ts/src/rings/number_field/number_field.js';
import { Rational } from '../../../packages/sagemath-ts/src/rings/rational.js';
import { hnf } from '../../../packages/sagemath-ts/src/rings/number_field/pari_nf.js';
import { lcm } from '../../../packages/sagemath-ts/src/arith/misc.js';
export function nf_ideal_audit(
  op: bigint,
  cs: bigint[],
  fn: bigint,
  fd: bigint,
  flat: bigint[],
  gd: bigint,
  count: bigint,
  exponent: bigint
): string {
  const K = new NumberField(new RationalPolynomial(cs.map((c) => new Rational(c * fn, fd))), 'a'),
    n = K.degree(),
    gens = Array.from(
      { length: Number(count) },
      (_, j) =>
        new NumberFieldElement(
          K,
          flat.slice(j * n, (j + 1) * n).map((c) => new Rational(c, gd))
        )
    ),
    I = K.ideal(...gens);
  const ideal = (J: typeof I) => {
    if (J.is_zero()) return ['0', []];
    const v = J.zk_basis().map((z) => z.list());
    let d = 1n;
    for (const row of v) for (const c of row) d = lcm(d, c.denominator);
    const H = hnf(
      v.map((row) => row.map((c) => c.numerator * (d / c.denominator))),
      n
    );
    return [String(J.norm()), H.map((row) => row.map((c) => String(new Rational(c, d))))];
  };
  let r: unknown,
    error = null;
  try {
    const O = () =>
      K.maximal_order() as import(
        '../../../packages/sagemath-ts/src/rings/number_field/order.js'
      ).AbsoluteOrder;
    if (op === 0n) r = ideal(I);
    else if (op === 1n) r = ideal(I.inverse());
    else if (op === 2n) r = ideal(I.mul(I.inverse()));
    else if (op === 3n) r = ideal(I.pow(exponent));
    else if (op === 4n) r = ideal(O().different());
    else if (op === 5n) r = ideal(O().codifferent());
    else if (op === 6n) r = [ideal(K.different()), K.different() === K.different()];
    else {
      const result = I.pow(op === 7n ? exponent : new Integer(exponent));
      r = [ideal(result), result === I];
    }
  } catch (caught) {
    const e = caught as Error;
    r = null;
    error = e.name + ': ' + e.message;
  }
  return JSON.stringify([r, error]);
}
