/** Compare ideal lattices canonically despite differing integral-basis choices. */
import {
  NumberField,
  NumberFieldElement,
  RationalPolynomial,
} from '../../../packages/sagemath-ts/src/rings/number_field/number_field.js';
import { Rational } from '../../../packages/sagemath-ts/src/rings/rational.js';
import { hnf } from '../../../packages/sagemath-ts/src/rings/number_field/pari_nf.js';
import { lcm } from '../../../packages/sagemath-ts/src/arith/misc.js';

export function nf_ideal_basis(
  op: bigint,
  cs: bigint[],
  fn: bigint,
  fd: bigint,
  flat: bigint[],
  gd: bigint,
  count: bigint,
  _exponent: bigint
): string {
  const K = new NumberField(new RationalPolynomial(cs.map((c) => new Rational(c * fn, fd))), 'a'),
    n = K.degree();
  const I = K.ideal(
    ...Array.from(
      { length: Number(count) },
      (_, j) =>
        new NumberFieldElement(
          K,
          flat.slice(j * n, (j + 1) * n).map((c) => new Rational(c, gd))
        )
    )
  );
  const lattice = (basis: NumberFieldElement[]) => {
    const rows = basis.map((z) => z.list()).filter((row) => row.some((c) => !c.isZero()));
    let d = 1n;
    for (const row of rows) for (const c of row) d = lcm(d, c.denominator);
    return rows.length
      ? hnf(
          rows.map((row) => row.map((c) => c.numerator * (d / c.denominator))),
          n
        ).map((row) => row.map((c) => String(new Rational(c, d))))
      : [];
  };
  let result: unknown,
    error: string | null = null;
  try {
    if (op === 12n) result = I.free_module() === I.free_module();
    else if (op === 9n || op === 10n) {
      const basis = op === 9n ? I.integral_basis() : I.free_module().basis;
      result = [basis.length, lattice(basis)];
    } else {
      const pair = I.gens_two(),
        [a, b] = pair;
      result = [
        String(a),
        JSON.stringify(lattice(K.ideal(a, b).zk_basis())) === JSON.stringify(lattice(I.zk_basis())),
        b.is_zero(),
        a.parent() === K,
        b.parent() === K,
        I.gens_two() === pair,
      ];
    }
  } catch (caught) {
    const e = caught as Error;
    result = null;
    error = e.name + ': ' + e.message;
  }
  return JSON.stringify([result, error]);
}
