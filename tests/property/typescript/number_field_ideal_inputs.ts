import {
  NumberField,
  NumberFieldElement,
  RationalPolynomial,
} from '../../../packages/sagemath-ts/src/rings/number_field/number_field.js';
import {
  NumberFieldIdeal,
  NumberFieldFractionalIdeal,
} from '../../../packages/sagemath-ts/src/rings/number_field/number_field_ideal.js';
import { Rational } from '../../../packages/sagemath-ts/src/rings/rational.js';
import { hnf } from '../../../packages/sagemath-ts/src/rings/number_field/pari_nf.js';
import { lcm } from '../../../packages/sagemath-ts/src/arith/misc.js';
import { setrand, randomi } from '../../../packages/parigp-ts/src/random.js';

export function nf_ideal_intersection_rng(
  cs: bigint[],
  fn: bigint,
  fd: bigint,
  aa: bigint[],
  ad: bigint,
  ac: bigint,
  bb: bigint[],
  bd: bigint,
  bc: bigint,
  seed: bigint
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
    );
  const I = ideal(aa, ad, ac),
    J = ideal(bb, bd, bc);
  K._pari_ideal_data();
  I.zk_basis();
  J.zk_basis();
  setrand(seed);
  const result = I.intersection(J),
    next = randomi(1n << 128n);
  return JSON.stringify([lattice(result), String(next)]);
}

const field = (cs: bigint[], name = 'a') =>
  new NumberField(new RationalPolynomial(cs.map((c) => new Rational(c))), name);
function lattice(I: NumberFieldIdeal): string[][] {
  if (I.is_zero()) return [];
  const rows = I.zk_basis().map((v) => v.list());
  let d = 1n;
  for (const row of rows) for (const c of row) d = lcm(d, c.denominator);
  return hnf(
    rows.map((row) => row.map((c) => c.numerator * (d / c.denominator))),
    I.number_field().degree()
  ).map((row) => row.map((v) => String(new Rational(v, d))));
}
export function nf_ideal_inputs(
  op: bigint,
  cs: bigint[],
  aa: bigint[],
  bb: bigint[],
  d: bigint
): string {
  const K = field(cs),
    a = new NumberFieldElement(
      K,
      aa.map((v) => new Rational(v, d))
    ),
    b = new NumberFieldElement(
      K,
      bb.map((v) => new Rational(v, d))
    ),
    I = K.ideal(a),
    J = K.ideal(b);
  let r: unknown,
    error: string | null = null;
  try {
    let L: NumberFieldIdeal;
    if (op === 0n) L = I.intersection(J);
    else if (op === 1n) L = I.intersection(new Rational(3n, d));
    else if (op === 2n) L = I.intersection(b);
    else if (op === 3n) L = I.intersection([b, new Rational(3n, d)]);
    else if (op === 4n) L = I.intersection([]);
    else if (op === 5n) L = I.intersection([0n]);
    else if (op === 6n) L = K.fractional_ideal(0n);
    else if (op === 7n) L = new NumberFieldFractionalIdeal(K, [K.zero()]);
    else if (op === 8n) L = K.ideal([a, b]);
    else if (op === 9n) L = K.ideal(I);
    else if (op === 10n) L = K.fractional_ideal(I);
    else if (op === 11n) L = K.fractional_ideal([a, b]);
    else if (op === 12n) L = new NumberFieldFractionalIdeal(K, []);
    else if (op === 13n) L = K.ideal(new NumberFieldIdeal(K, [K.__call__(3n)]));
    else if (op === 14n) L = I.intersection(new NumberFieldIdeal(K, [K.__call__(3n)]));
    else if (op === 15n) {
      r = K.ideal(0n).toString();
      return JSON.stringify([r, null]);
    } else if (op === 16n) {
      r = new NumberFieldIdeal(K, [K.__call__(3n)]).toString();
      return JSON.stringify([r, null]);
    } else if (op === 17n) L = I.add(J);
    else if (op === 18n) L = K.ideal(I.add(J));
    else L = I.intersection(I.add(J));
    r = [lattice(L), L === I, L instanceof NumberFieldFractionalIdeal];
    if (op === 17n) (r as unknown[]).push(L.gens().map((v) => v.list().map(String)));
  } catch (e) {
    r = null;
    error = (e as Error).name + ': ' + (e as Error).message;
  }
  return JSON.stringify([r, error]);
}
export function nf_cross_ideal(cs: bigint[], mode: bigint, op: bigint): string {
  const K = field([2n, 0n, 1n]),
    L = field(cs, 'b');
  const v = [L.zero(), L.__call__(2n), L.gen(), L.ideal(0n), L.ideal(2n), L.ideal(L.gen())][
    Number(mode)
  ]!;
  let r: unknown,
    error: string | null = null;
  try {
    r =
      op === 0n
        ? K.__call__(v as NumberFieldElement)
            .list()
            .map(String)
        : lattice(op === 1n ? K.ideal(v) : K.fractional_ideal(v));
  } catch (e) {
    r = null;
    error = (e as Error).name + ': ' + (e as Error).message;
  }
  return JSON.stringify([r, error]);
}

import { Integer } from '../../../packages/sagemath-ts/src/rings/integer_ring.js';
export function nf_ideal_factory(
  cs: bigint[],
  op: bigint,
  kind: bigint,
  layout: bigint,
  n: bigint,
  d: bigint
): string {
  const K = field(cs);
  let r: unknown,
    error: string | null = null;
  try {
    const value =
      kind === 0n
        ? n
        : kind === 1n
          ? new Integer(n)
          : kind === 2n
            ? new Rational(n, d)
            : kind === 3n
              ? K.__call__(new Rational(n, d))
              : Number(n) / Number(d);
    const args: any[] =
      layout === 0n
        ? [value]
        : layout === 1n
          ? [[value]]
          : layout === 2n
            ? [value, value]
            : layout === 3n
              ? []
              : [[]];
    const I =
      op === 0n
        ? K.ideal(...args)
        : op === 1n
          ? K.fractional_ideal(...args)
          : op === 2n
            ? new NumberFieldIdeal(K, args)
            : new NumberFieldFractionalIdeal(K, args);
    r = [
      lattice(I),
      I instanceof NumberFieldFractionalIdeal,
      I.gens().map((v) => v.list().map(String)),
    ];
  } catch (e) {
    r = null;
    error = (e as Error).name + ': ' + (e as Error).message;
  }
  return JSON.stringify([r, error]);
}
