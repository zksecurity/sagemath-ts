import { znorder } from '../../../../packages/parigp-ts/src/ff.js';
import { set_random_seed } from '../../../../packages/sagemath-ts/src/misc/randstate.js';
import { GF } from '../../../../packages/sagemath-ts/src/rings/finite_rings/finite_field_constructor.js';
import { FiniteFieldPrime } from '../../../../packages/sagemath-ts/src/rings/finite_rings/finite_field_prime.js';
/** Direct residue-ring comparisons, including result parents and exact errors. */
import {
  IntegerMod,
  Mod,
} from '../../../../packages/sagemath-ts/src/rings/finite_rings/integer_mod.js';
import {
  IntegerModRing,
  IntegerModRingFactory,
  Integers,
  Zmod,
} from '../../../../packages/sagemath-ts/src/rings/finite_rings/integer_mod_ring.js';
import { Integer, ZZ } from '../../../../packages/sagemath-ts/src/rings/integer_ring.js';
import { Rational } from '../../../../packages/sagemath-ts/src/rings/rational.js';

function comparison(run: () => unknown): string {
  const normalize = (x: unknown): unknown =>
    Array.isArray(x) ? x.map(normalize) : typeof x === 'boolean' || x === null ? x : String(x);
  try {
    return JSON.stringify({ value: normalize(run()) });
  } catch (e) {
    return JSON.stringify({ error: (e as Error).name, message: (e as Error).message });
  }
}
function scalar(mode: bigint, a: bigint, b: bigint): unknown {
  if (mode === 0n) return new Integer(a);
  if (mode === 1n) return new Rational(a, b);
  if (mode === 2n) return Boolean(a);
  if (mode === 3n) return null;
  if (mode === 4n) return String(a);
  if (mode === 5n) return Number(a) / Number(b);
  if (mode === 6n) return Zmod(b).__call__(a);
  if (mode === 7n) return GF(b).__call__(a);
  return [];
}
function invoke(x: object, name: string, ...args: unknown[]): unknown {
  return (x as Record<string, (...args: unknown[]) => unknown>)[name]!.apply(x, args);
}
export const functions: Record<string, (...args: never[]) => unknown> = {
  mi_constructor: (n: bigint, mode: bigint, a: bigint, b: bigint, direct: bigint) =>
    comparison(() =>
      direct
        ? new IntegerMod(scalar(mode, a, b) as bigint, Zmod(n))
        : Zmod(n).__call__(scalar(mode, a, b))
    ),
  mi_ring: (n: bigint) =>
    comparison(() => {
      const R = Zmod(n);
      return [
        String(R),
        R.characteristic,
        R.order,
        R.zero(),
        R.one(),
        R.gen(),
        R.is_field(),
        R.cardinality(),
        [...R],
        R.list(),
        R.units(),
        R.multiplicative_group_is_cyclic(),
        R.unit_gens(),
      ];
    }),
  mi_random: (n: bigint, seed: bigint, bound: bigint, mode: bigint) =>
    comparison(() => {
      set_random_seed(seed);
      const R = Zmod(n);
      return Array.from({ length: 20 }, () =>
        mode === 0n ? R.random_element() : invoke(R, 'random_element', bound)
      );
    }),
  mi_ring_direct: (n: bigint) => comparison(() => String(new IntegerModRing(n))),
};
for (const method of ['add', 'sub', 'mul', 'div', 'eq', 'pow']) {
  functions['mi_' + method] = (n: bigint, a: bigint, mode: bigint, b: bigint, d: bigint) =>
    comparison(() => {
      const z = invoke(Zmod(n).__call__(a), method, scalar(mode, b, d));
      return typeof z === 'boolean'
        ? z
        : [
            z,
            typeof z === 'string'
              ? "<class 'str'>"
              : Array.isArray(z)
                ? "<class 'list'>"
                : String((z as IntegerMod).parent),
          ];
    });
}
for (const method of [
  'neg',
  'inv',
  'isZero',
  'isOne',
  'isUnit',
  'modulus',
  'lift',
  'toBigInt',
  'toString',
  'repr',
  'multiplicative_order',
]) {
  functions['mi_' + method] = (n: bigint, a: bigint) =>
    comparison(() => {
      const x = Zmod(n).__call__(a);
      return method === 'modulus' ? x.modulus : invoke(x, method);
    });
}

functions.mi_pari_order = (a: bigint, n: bigint, order: bigint) =>
  comparison(() => znorder(a, n, order === 0n ? undefined : order));

functions.mi_string = (n: bigint, codes: bigint[]) =>
  comparison(() => Zmod(n).__call__(String.fromCodePoint(...codes.map(Number))));

functions.mi_factory = (mode: bigint, kind: bigint, a: bigint, b: bigint) =>
  comparison(() => {
    const factory = [Zmod, Integers, IntegerModRingFactory][Number(mode)]!;
    const R =
      kind === 9n
        ? (factory as () => unknown)()
        : (factory as (n: unknown) => unknown)(scalar(kind, a, b));
    const ring = R as IntegerModRing;
    return [
      String(R),
      R === ZZ,
      typeof ring.characteristic === 'function'
        ? (ring.characteristic as () => bigint)()
        : ring.characteristic,
      ring.__call__(12n),
    ];
  });
functions.mi_mod_factory = (
  n: bigint,
  kind: bigint,
  a: bigint,
  b: bigint,
  parentOrder: bigint,
  withParent: bigint
) =>
  comparison(() => {
    const value = scalar(kind, a, b);
    const make = Mod as (x: unknown, n: bigint, parent?: unknown) => unknown;
    const result = withParent ? make(value, n, Zmod(parentOrder)) : make(value, n);
    return n === 0n ? [result, result === value] : [result, String((result as IntegerMod).parent)];
  });

functions.mi_factory_identity = (n: bigint) =>
  comparison(() => [
    Zmod(n) === Zmod(new Integer(n)),
    Zmod(n) === Zmod(-n),
    Zmod(n) === Zmod(new Rational(n)),
    Zmod(n) === Integers(n),
    (Mod(1n, n) as IntegerMod).parent === Zmod(n),
  ]);

functions.mi_mod_prime_parent = (
  n: bigint,
  p: bigint,
  kind: bigint,
  a: bigint,
  b: bigint,
  legacy: bigint
) =>
  comparison(() => {
    const parent = legacy ? new FiniteFieldPrime(p) : GF(p);
    const result = (Mod as (value: unknown, n: bigint, parent: unknown) => unknown)(
      scalar(kind, a, b),
      n,
      parent
    );
    if (n === 0n) return result;
    return [
      String(result),
      String((result as IntegerMod).parent),
      invoke(result as object, 'is_square'),
    ];
  });

import { PolynomialRing } from '../../../../packages/sagemath-ts/src/rings/polynomial/polynomial_ring.js';
import type {
  CoefficientRing,
  RingElement,
} from '../../../../packages/sagemath-ts/src/rings/polynomial/polynomial_element.js';
functions.mi_polynomial_roots = (n: bigint, coefficients: bigint[], multiplicities: bigint) =>
  comparison(() => {
    const K = Zmod(n),
      R = new PolynomialRing(K as unknown as CoefficientRing<IntegerMod & RingElement>, 'x');
    const f = R.__call__(coefficients);
    if (multiplicities)
      return f
        .roots()
        .map(([r, m]) => [
          String(r),
          String(m),
          String(r.parent),
          typeof (K as unknown as { field?: () => unknown }).field === 'function' &&
            r.parent === (K as unknown as { field: () => unknown }).field(),
        ]);
    return f
      .roots({ multiplicities: false })
      .map((r) => [String(r), String(r.parent), r.parent === K]);
  });
functions.mi_field = (n: bigint, direct: bigint, values: bigint[]) =>
  comparison(() => {
    const K = direct ? new IntegerModRing(n) : (Zmod(n) as IntegerModRing),
      F = K.field();
    return [
      String(F),
      F.order,
      F.characteristic,
      F.degree,
      F === K.field(),
      values.map((v) => [String(F.__call__(v)), F.__call__(v).parent === F]),
    ];
  });
functions.mi_factored_order = (n: bigint, direct: bigint) =>
  comparison(() => {
    const K = direct ? new IntegerModRing(n) : (Zmod(n) as IntegerModRing),
      result = K.factored_order();
    return [result, result === K.factored_order()];
  });
functions.mi_residue_root_lift = (p: bigint, e: bigint, coefficients: bigint[], root: bigint) =>
  comparison(() => {
    const K = Zmod(p ** e) as IntegerModRing,
      R = new PolynomialRing(K as unknown as CoefficientRing<IntegerMod & RingElement>, 'x');
    const f = R.__call__(coefficients),
      r = (Zmod(p) as IntegerModRing).__call__(root);
    return IntegerModRing._lift_residue_field_root(p, e, f, f.derivative(), r).map((v) => [
      String(v),
      String(v.parent),
      v.parent === K,
    ]);
  });

functions.mi_modular_roots_hook = (
  n: bigint,
  coefficients: bigint[],
  multiplicities: bigint,
  hook: bigint
) =>
  comparison(() => {
    const K = Zmod(n) as IntegerModRing,
      R = new PolynomialRing(K as unknown as CoefficientRing<IntegerMod & RingElement>, 'x');
    const f = R.__call__(coefficients);
    const result = K._roots_univariate_polynomial(f, {
      multiplicities: !!multiplicities,
      ring: hook === 0n ? null : hook === 1n ? K : (Zmod(n + 1n) as IntegerModRing),
      algorithm: 'ignored',
    });
    return multiplicities
      ? (
          result as Array<
            [
              import(
                '../../../../packages/sagemath-ts/src/rings/finite_rings/finite_field_prime.js'
              ).FiniteFieldElement,
              number,
            ]
          >
        ).map(([r, m]) => [String(r), String(m), String(r.parent), r.parent === K.field()])
      : (result as IntegerMod[]).map((r) => [String(r), String(r.parent), r.parent === K]);
  });
