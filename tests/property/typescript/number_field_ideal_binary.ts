import { idealdiv } from '../../../packages/parigp-ts/src/base4.js';
import {
  NumberField,
  NumberFieldElement,
  RationalPolynomial,
} from '../../../packages/sagemath-ts/src/rings/number_field/number_field.js';
import { Rational } from '../../../packages/sagemath-ts/src/rings/rational.js';
import { hnf } from '../../../packages/sagemath-ts/src/rings/number_field/pari_nf.js';
import { lcm } from '../../../packages/sagemath-ts/src/arith/misc.js';
export function nf_ideal_binary(
  op: bigint,
  cs: bigint[],
  fn: bigint,
  fd: bigint,
  aa: bigint[],
  ad: bigint,
  ac: bigint,
  bb: bigint[],
  bd: bigint,
  bc: bigint
): string {
  const K = new NumberField(new RationalPolynomial(cs.map((c) => new Rational(c * fn, fd))), 'a'),
    n = K.degree();
  const ideal = (flat: bigint[], d: bigint, count: bigint) =>
      K.ideal(
        ...Array.from(
          { length: Number(count) },
          (_, j) =>
            new NumberFieldElement(
              K,
              flat.slice(j * n, (j + 1) * n).map((c) => new Rational(c, d))
            )
        )
      ),
    I = ideal(aa, ad, ac),
    J = ideal(bb, bd, bc);
  const fmt = (L: typeof I) => {
    if (L.is_zero()) return ['0', []];
    const rows = L.zk_basis().map((a) => a.list());
    let d = 1n;
    for (const row of rows) for (const c of row) d = lcm(d, c.denominator);
    const H = hnf(
      rows.map((row) => row.map((c) => c.numerator * (d / c.denominator))),
      n
    );
    return [String(L.norm()), H.map((row) => row.map((c) => String(new Rational(c, d))))];
  };
  let r: unknown,
    error = null;
  try {
    if (op === 9n) {
      const hnfOf = (L: typeof I) =>
        L.is_zero()
          ? { entries: [], denominator: 1n }
          : (
              L as unknown as { _computeHNF(): { entries: bigint[][]; denominator: bigint } }
            )._computeHNF();
      const A = hnfOf(I),
        B = hnfOf(J),
        [C, d] = idealdiv(K._pari_ideal_data(), A.entries, A.denominator, B.entries, B.denominator);
      const basis = K._pari_integral_basis();
      r = C.length
        ? fmt(
            K.ideal(
              ...C.map((column) =>
                column.reduce((sum, c, i) => sum.add(basis[i]!.mul(new Rational(c, d))), K.zero())
              )
            )
          )
        : fmt(K.ideal(0n));
    } else if (op === 8n) r = fmt(I.mul(J));
    else if (op === 0n) r = fmt(I.add(J));
    else if (op === 1n) r = fmt(I.intersection(J));
    else if (op === 2n) r = I.is_coprime(J);
    else if (op === 3n) r = I.divides(J);
    else if (op === 4n) r = fmt(I.div(J));
    else if (op === 5n) {
      const D = I.denominator();
      r = [fmt(D), D === I.denominator(), D.number_field() === K];
    } else if (op === 6n) {
      const N = I.numerator();
      r = [fmt(N), N === I.numerator(), N === I];
    } else {
      const N = I.numerator(),
        D = I.denominator();
      r = [fmt(N.div(D)), N.is_coprime(D)];
    }
  } catch (caught) {
    const e = caught as Error;
    r = null;
    error = e.name + ': ' + e.message;
  }
  return JSON.stringify([r, error]);
}
