import {
  NumberField,
  NumberFieldElement,
  RationalPolynomial,
} from '../../../packages/sagemath-ts/src/rings/number_field/number_field.js';
import { Rational } from '../../../packages/sagemath-ts/src/rings/rational.js';
export function nf_ideal_construction(
  cs: bigint[],
  fn: bigint,
  fd: bigint,
  flat: bigint[],
  den: bigint,
  count: bigint
): string {
  const K = new NumberField(new RationalPolynomial(cs.map((c) => new Rational(c * fn, fd))), 'a'),
    n = K.degree(),
    gens = Array.from(
      { length: Number(count) },
      (_, j) =>
        new NumberFieldElement(
          K,
          flat.slice(j * n, (j + 1) * n).map((c) => new Rational(c, den))
        )
    );
  let r: unknown,
    error = null;
  try {
    const I = K.ideal(...gens);
    r = [
      I.gens().map((a) => a.list().map(String)),
      I.ngens(),
      I === K.ideal(...gens),
      I === K.ideal(0n),
    ];
  } catch (caught) {
    const e = caught as Error;
    r = null;
    error = e.name + ': ' + e.message;
  }
  return JSON.stringify([r, error]);
}
