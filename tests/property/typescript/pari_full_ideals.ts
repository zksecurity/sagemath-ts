import {
  NumberField,
  RationalPolynomial,
} from '../../../packages/sagemath-ts/src/rings/number_field/number_field.js';
import { idealprimedec_limit_f } from '../../../packages/parigp-ts/src/base2.js';
import { pr_hnf } from '../../../packages/parigp-ts/src/base4.js';
export function pari_full_primedec(cs: bigint[], p: bigint, limit: bigint): string {
  const K = new NumberField(RationalPolynomial.fromBigInts(cs), 'a'),
    nf = K._pari_ideal_data();
  return JSON.stringify(
    idealprimedec_limit_f(nf, p, limit).map((P) => [
      [P.p, P.generator, P.e, P.f, P.tau],
      pr_hnf(nf, P),
    ]),
    (_, x) => (typeof x === 'bigint' ? String(x) : x)
  );
}

import { Rational } from '../../../packages/sagemath-ts/src/rings/rational.js';
import { idealfactor, idealnumden, idealismaximal } from '../../../packages/parigp-ts/src/base4.js';
export function pari_full_ideal(
  op: bigint,
  cs: bigint[],
  flat: bigint[],
  d: bigint,
  count: bigint
): string {
  const K = new NumberField(RationalPolynomial.fromBigInts(cs), 'a'),
    n = K.degree(),
    nf = K._pari_ideal_data();
  const I = K.ideal(
    Array.from({ length: Number(count) }, (_, i) =>
      K.__call__(flat.slice(i * n, (i + 1) * n).map((c) => new Rational(c, d)))
    )
  );
  if (op === 5n) {
    const inverse = ratInverse(K._pari_integral_basis().map((b) => b.list())),
      power = flat.slice(0, n).map((c) => new Rational(c, d));
    const coordinates = Array.from({ length: n }, (_, i) =>
        power.reduce((s, c, j) => s.add(c.mul(inverse[j]![i]!)), Rational.zero())
      ),
      den = coordinates.reduce((d, c) => lcm(d, c.denominator), 1n),
      M = zk_multable(
        nf.multiplication,
        coordinates.map((c) => c.numerator * (den / c.denominator))
      ),
      [a, b] = zkmultable_inv(M);
    return JSON.stringify([a.map((c) => String(new Rational(c, b))), String(zkmultable_capZ(M))]);
  }
  const H = I.is_zero() ? { entries: [], denominator: 1n } : (I as any)._computeHNF();
  if (op === 4n)
    return JSON.stringify([H.entries, H.denominator], (_, x) =>
      typeof x === 'bigint' ? String(x) : x
    );
  const prime = (P: any) => [P.p, P.generator, P.e, P.f, P.tau];
  const P = op === 3n ? idealismaximal(nf, H.entries, H.denominator) : null;
  const result =
    op === 1n
      ? idealfactor(nf, H.entries, H.denominator).map(([P, e]) => [prime(P), pr_hnf(nf, P), e])
      : op === 2n
        ? idealnumden(nf, H.entries, H.denominator)
        : P
          ? prime(P)
          : null;
  return JSON.stringify(result, (_, x) => (typeof x === 'bigint' ? String(x) : x));
}

import { ratInverse } from '../../../packages/sagemath-ts/src/rings/number_field/pari_nf.js';
import { zk_multable } from '../../../packages/parigp-ts/src/base4.js';
import { zkmultable_inv, zkmultable_capZ } from '../../../packages/parigp-ts/src/base3.js';
import { lcm } from '../../../packages/sagemath-ts/src/arith/misc.js';

export function pari_full_primedec_no_index(cs:bigint[],p:bigint,limit:bigint):string {
 const K=new NumberField(RationalPolynomial.fromBigInts(cs),'a'),nf={...K._pari_ideal_data()};
 delete nf.index;
 return JSON.stringify(idealprimedec_limit_f(nf,p,limit).map(P=>[[P.p,P.generator,P.e,P.f,P.tau],pr_hnf(nf,P)]),(_,x)=>typeof x==='bigint'?String(x):x);
}
