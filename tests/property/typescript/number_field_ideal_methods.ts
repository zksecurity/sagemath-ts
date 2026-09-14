import {
  NumberField,
  RationalPolynomial,
} from '../../../packages/sagemath-ts/src/rings/number_field/number_field.js';
import { Integer } from '../../../packages/sagemath-ts/src/rings/integer_ring.js';
import { Rational } from '../../../packages/sagemath-ts/src/rings/rational.js';
import { hnf } from '../../../packages/sagemath-ts/src/rings/number_field/pari_nf.js';
import { lcm } from '../../../packages/sagemath-ts/src/arith/misc.js';
import { NumberFieldIdeal } from '../../../packages/sagemath-ts/src/rings/number_field/number_field_ideal.js';
const field = (cs: bigint[], name = 'a') =>
  new NumberField(RationalPolynomial.fromBigInts(cs), name);
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
export function nf_ideal_method_input(
  cs: bigint[],
  idealKind: bigint,
  op: bigint,
  kind: bigint,
  n: bigint,
  d: bigint,
  otherCs: bigint[] = cs
): string {
  const K = field(cs),
    L = field(otherCs, 'b'),
    q = new Rational(n, d);
  const I = [
    K.ideal(0n),
    K.ideal(1n),
    K.ideal(2n),
    K.ideal(new Rational(1n, 2n)),
    K.ideal(K.gen()),
    K.ideal(3n, K.gen().add(1n)),
  ][Number(idealKind)]!;
  const values = [
    n,
    new Integer(n),
    q,
    Number(n) / Number(d),
    K.__call__(q),
    K.gen().add(q),
    [q],
    [],
    Array(K.degree()).fill(q),
    [q, K.gen()],
    K.ideal(q),
    K.ideal(0n),
    L.__call__(q),
    L.gen().add(q),
    L.ideal(q),
    L.ideal(L.gen().add(q)),
    [0n],
    K.ideal(3n, K.gen().add(1n)),
    new NumberFieldIdeal(K, [K.__call__(q)]),
  ];
  const x = values[Number(kind)] as any;
  let r: unknown,
    error: string | null = null;
  try {
    r =
      op === 0n
        ? I.contains(x)
        : op === 1n
          ? I.divides(x)
          : op === 2n
            ? I.is_coprime(x)
            : op === 3n
              ? lattice(I.add(x))
              : op === 4n
                ? lattice(I.mul(x))
                : lattice(I.div(x));
  } catch (e) {
    r = null;
    error = (e as Error).name + ': ' + (e as Error).message;
  }
  return JSON.stringify([r, error]);
}
export function nf_field_vector(
  cs: bigint[],
  kind: bigint,
  length: bigint,
  n: bigint,
  d: bigint
): string {
  const K = field(cs),
    L = field(cs, 'b');
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
              ? Number(n) / Number(d)
              : kind === 4n
                ? K.__call__(new Rational(n, d))
                : kind === 5n
                  ? L.__call__(new Rational(n, d))
                  : kind === 6n
                    ? K.gen()
                    : L.gen();
    r = K.__call__(Array(Number(length)).fill(value))
      .list()
      .map(String);
  } catch (e) {
    r = null;
    error = (e as Error).name + ': ' + (e as Error).message;
  }
  return JSON.stringify([r, error]);
}

export function nf_basis_class_number(cs: bigint[], fn: bigint, fd: bigint, op: bigint): string {
  const K = new NumberField(new RationalPolynomial(cs.map((c) => new Rational(c * fn, fd))), 'a');
  let r: unknown,
    error: string | null = null;
  try {
    r = String(op === 0n ? K.class_number() : K.class_group().order());
  } catch (e) {
    r = null;
    error = (e as Error).name + ': ' + (e as Error).message;
  }
  return JSON.stringify([r, error]);
}
