import {
  NumberField,
  RationalPolynomial,
} from '../../../packages/sagemath-ts/src/rings/number_field/number_field.js';
import { NumberFieldIdeal } from '../../../packages/sagemath-ts/src/rings/number_field/number_field_ideal.js';
import { Rational } from '../../../packages/sagemath-ts/src/rings/rational.js';
import { ValueError } from '../../../packages/sagemath-ts/src/errors.js';
import { Integer } from '../../../packages/sagemath-ts/src/rings/integer_ring.js';
import { hnf } from '../../../packages/sagemath-ts/src/rings/number_field/pari_nf.js';
import { lcm } from '../../../packages/sagemath-ts/src/arith/misc.js';
function lattice(I: NumberFieldIdeal): string[][] {
  if (I.is_zero()) return [];
  const rows = I.zk_basis().map((x) => x.list());
  let d = 1n;
  for (const row of rows) for (const c of row) d = lcm(d, c.denominator);
  return hnf(
    rows.map((row) => row.map((c) => c.numerator * (d / c.denominator))),
    I.number_field().degree()
  ).map((row) => row.map((c) => String(new Rational(c, d))));
}
export function nf_ideal_valuation(
  cs: bigint[],
  p: bigint,
  flat: bigint[],
  d: bigint,
  count: bigint,
  base: bigint
): string {
  const K = new NumberField(RationalPolynomial.fromBigInts(cs), 'a'),
    n = K.degree(),
    gens = Array.from({ length: Number(count) }, (_, i) =>
      K.__call__(flat.slice(i * n, (i + 1) * n).map((c) => new Rational(c, d)))
    ),
    I = base ? new NumberFieldIdeal(K, gens) : K.ideal(gens);
  const out = K.decomposition(p).map(([P]) => {
    const v = I.valuation(P as any);
    return [lattice(P), v === 'Infinity' ? '+Infinity' : String(v)];
  });
  out.sort((a, b) =>
    JSON.stringify(a[0]) < JSON.stringify(b[0])
      ? -1
      : JSON.stringify(a[0]) > JSON.stringify(b[0])
        ? 1
        : 0
  );
  return JSON.stringify(out);
}
export function nf_ideal_valuation_input(
  cs: bigint[],
  receiver: bigint,
  kind: bigint,
  n: bigint,
  d: bigint
): string {
  const K = new NumberField(RationalPolynomial.fromBigInts(cs), 'a'),
    L = new NumberField(RationalPolynomial.fromBigInts(cs), 'b'),
    q = new Rational(n, d);
  const I = [
    K.ideal(0n),
    K.ideal(1n),
    K.ideal(6n),
    K.ideal(K.gen().add(1n)),
    new NumberFieldIdeal(K, [K.zero()]),
    new NumberFieldIdeal(K, [K.__call__(6n)]),
  ][Number(receiver)]!;
  const values = [
    n,
    new Integer(n),
    q,
    Number(n) / Number(d),
    K.__call__(q),
    K.gen().add(q),
    [],
    [q],
    Array(K.degree()).fill(q),
    K.ideal(q),
    K.ideal(0n),
    new NumberFieldIdeal(K, [K.__call__(q)]),
    L.__call__(q),
    L.gen().add(q),
    L.ideal(q),
    L.decomposition(2n)[0]![0],
    '2',
    [K.gen(), q],
  ];
  try {
    const v = I.valuation(values[Number(kind)] as any);
    return v === 'Infinity' ? '+Infinity' : String(v);
  } catch (error) {
    if (
      error instanceof ValueError &&
      (kind === 8n || kind === 17n) &&
      error.message.startsWith('p (= ')
    ) {
      const P = K.ideal(values[Number(kind)] as any);
      throw new ValueError(
        'p (= ' +
          JSON.stringify(lattice(P)) +
          ')' +
          error.message.slice(error.message.lastIndexOf(')') + 1)
      );
    }
    throw error;
  }
}

export function nf_prime_below(cs: bigint[], p: bigint): string {
  const K = new NumberField(RationalPolynomial.fromBigInts(cs), 'a'),
    out = K.decomposition(p).map(([P]) => [lattice(P), String(P.prime_below())]);
  out.sort((a, b) =>
    JSON.stringify(a[0]) < JSON.stringify(b[0])
      ? -1
      : JSON.stringify(a[0]) > JSON.stringify(b[0])
        ? 1
        : 0
  );
  return JSON.stringify(out);
}
