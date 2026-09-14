import {
  NumberField,
  RationalPolynomial,
} from '../../../packages/sagemath-ts/src/rings/number_field/number_field.js';
import { Rational } from '../../../packages/sagemath-ts/src/rings/rational.js';
import { primedec_end } from '../../../packages/parigp-ts/src/base2.js';
import { idealval } from '../../../packages/parigp-ts/src/base4.js';
const encode = (v: unknown) => JSON.stringify(v, (_, x) => (typeof x === 'bigint' ? String(x) : x));
const compare = (A: bigint[][], B: bigint[][]) => {
  for (let j = 0; j < A.length; j++)
    for (let i = 0; i < A[j]!.length; i++) {
      if (A[j]![i]! < B[j]![i]!) return -1;
      if (A[j]![i]! > B[j]![i]!) return 1;
    }
  return 0;
};
function data(cs: bigint[], p: bigint) {
  const K = new NumberField(RationalPolynomial.fromBigInts(cs), 'a'),
    nf = K._pari_ideal_data(),
    H = K.decomposition(p)
      .map(([P]) => (P as any)._computeHNF().entries as bigint[][])
      .sort(compare);
  return { K, nf, H, primes: primedec_end(nf, H, p) };
}
export function pari_prime_valuation_data(cs: bigint[], p: bigint): string {
  const { H, primes } = data(cs, p);
  return encode(primes.map((P, i) => [H[i], [P.p, P.generator, P.e, P.f, P.tau]]));
}
export function pari_ideal_valuation(
  cs: bigint[],
  p: bigint,
  flat: bigint[],
  d: bigint,
  count: bigint
): string {
  const { K, nf, H, primes } = data(cs, p),
    n = K.degree(),
    gens = Array.from({ length: Number(count) }, (_, i) =>
      K.__call__(flat.slice(i * n, (i + 1) * n).map((c) => new Rational(c, d)))
    ),
    I = K.ideal(gens),
    h = I.is_zero() ? { entries: [], denominator: 1n } : (I as any)._computeHNF();
  return encode(
    primes.map((P, i) => {
      const v = idealval(nf, h.entries, h.denominator, P);
      return [H[i], v === 'Infinity' ? '+oo' : v];
    })
  );
}

import { ZC_nfval } from '../../../packages/parigp-ts/src/base3.js';
export function pari_prime_element_valuation(cs: bigint[], p: bigint, flat: bigint[]): string {
  const { H, primes } = data(cs, p);
  return encode(primes.map((P, i) => [H[i], ZC_nfval(flat, P)]));
}
