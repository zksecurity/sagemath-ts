import {flint_remainder} from '../flint_remainder.js';
import { rdf_scalar, rdf_alias } from '../real_double_group.js';
import { AttributeError } from '../../../../packages/sagemath-ts/src/errors.js';
import { GF2 } from '../../../../packages/sagemath-ts/src/rings/finite_rings/gf2.js';
/**
 * sagemath-ts side of the `polynomials` property-test area (GF(p)[x]).
 *
 * Cases: tests/property/cases/polynomials.cases.json
 * SageMath counterpart: tests/property/python/areas/polynomials.py
 */

import {
  polyAdd,
  polyDerivative,
  polyEval,
  polyFactor,
  polyGcd,
  polyIsIrreducible,
  polyMod,
  polyMul,
  polyPow,
  polyQuoRem,
  polyRoots,
} from './_helpers.js';

export const functions: Record<string, (...args: never[]) => unknown> = {
  poly_add: (p: bigint, coeffs1: bigint[], coeffs2: bigint[]) => polyAdd(p, coeffs1, coeffs2),
  poly_mul: (p: bigint, coeffs1: bigint[], coeffs2: bigint[]) => polyMul(p, coeffs1, coeffs2),
  poly_quo_rem: (p: bigint, coeffs1: bigint[], coeffs2: bigint[]) =>
    polyQuoRem(p, coeffs1, coeffs2),
  poly_mod: (p: bigint, coeffs1: bigint[], coeffs2: bigint[]) => polyMod(p, coeffs1, coeffs2),
  poly_gcd: (p: bigint, coeffs1: bigint[], coeffs2: bigint[]) => polyGcd(p, coeffs1, coeffs2),
  poly_eval: (p: bigint, coeffs: bigint[], x: bigint) => polyEval(p, coeffs, x),
  poly_factor: (p: bigint, coeffs: bigint[]) => polyFactor(p, coeffs),
  poly_derivative: (p: bigint, coeffs: bigint[]) => polyDerivative(p, coeffs),
  poly_is_irreducible: (p: bigint, coeffs: bigint[]) => polyIsIrreducible(p, coeffs),
  poly_roots: (p: bigint, coeffs: bigint[]) => polyRoots(p, coeffs),
  poly_pow: (p: bigint, coeffs: bigint[], n: bigint) => polyPow(p, coeffs, n),
};

// Constructor comparisons deliberately pass uncoerced inputs to the public API.
import {
  GFpn,
  PrimeField,
} from '../../../../packages/sagemath-ts/src/rings/finite_rings/finite_field_extension.js';
import { FiniteFieldPrime } from '../../../../packages/sagemath-ts/src/rings/finite_rings/finite_field_prime.js';
import { Zmod } from '../../../../packages/sagemath-ts/src/rings/finite_rings/integer_mod_ring.js';
import { Integer, ZZ } from '../../../../packages/sagemath-ts/src/rings/integer_ring.js';
import { Polynomial } from '../../../../packages/sagemath-ts/src/rings/polynomial/polynomial_element.js';
import type {
  CoefficientRing,
  RingElement,
} from '../../../../packages/sagemath-ts/src/rings/polynomial/polynomial_element.js';
import { PolynomialRing } from '../../../../packages/sagemath-ts/src/rings/polynomial/polynomial_ring.js';
import { Rational } from '../../../../packages/sagemath-ts/src/rings/rational.js';
import { QQ } from '../../../../packages/sagemath-ts/src/rings/rational_field.js';

function comparison(run: () => unknown): string {
  const normalized = (x: unknown): unknown =>
    Array.isArray(x) ? x.map(normalized) : typeof x === 'boolean' || x === null ? x : String(x);
  try {
    return JSON.stringify({ value: normalized(run()) });
  } catch (e) {
    return JSON.stringify({ error: (e as Error).name, message: (e as Error).message });
  }
}
function polynomialBase(kind: bigint): CoefficientRing<RingElement> {
  const base =
    kind === 0n
      ? QQ
      : kind === 1n
        ? new PrimeField(7n)
        : kind === 2n
          ? new FiniteFieldPrime(7n)
          : kind === 3n
            ? Zmod(14n)
            : kind === 4n
              ? GFpn(7n, 2, undefined, 'a')
              : new PrimeField(2n);
  return base as unknown as CoefficientRing<RingElement>;
}
function polynomialInput(kind: bigint, n: bigint): unknown {
  if (kind === 26n) return Zmod(1n).__call__(n);
  if (kind === 24n) return 1 / 3;
  if (kind === 25n) return Number.MIN_VALUE;
  if (kind === 21n) return GFpn(7n, 2, undefined, 'a').__call__(n);
  if (kind === 22n) return GF2.__call__(n);
  if (kind === 23n) return new FiniteFieldPrime(7n).__call__(n);
  return kind === 0n
    ? n
    : kind === 1n
      ? new Integer(n)
      : kind === 2n
        ? new Rational(n, 2n)
        : kind === 3n
          ? Boolean(n)
          : kind === 4n
            ? String(n)
            : kind === 5n
              ? null
              : kind === 6n
                ? Number(n) + 0.5
                : kind === 7n
                  ? new PrimeField(7n).__call__(n)
                  : kind === 8n
                    ? Zmod(14n).__call__(n)
                    : kind === 9n
                      ? new PrimeField(5n).__call__(n)
                      : kind === 10n
                        ? []
                        : kind === 11n
                          ? [n]
                          : kind === 12n
                            ? GFpn(7n, 2, undefined, 'a').gen().add(n)
                            : kind === 13n
                              ? [n, 1n]
                              : kind === 14n
                                ? [[n, 1n]]
                                : kind === 15n
                                  ? [[n]]
                                  : kind === 16n
                                    ? [null]
                                    : kind === 17n
                                      ? [[]]
                                      : Number.NaN;
}
functions.polynomial_constructor_input = (base: bigint, kind: bigint, n: bigint, shape: bigint) =>
  comparison(() => {
    const R = new PolynomialRing(polynomialBase(base), 'x');
    const c = polynomialInput(kind, n);
    const input =
      shape === 0n
        ? [c]
        : shape === 1n
          ? [c, 0n, 1n, 0n]
          : shape === 2n
            ? []
            : shape === 3n
              ? null
              : undefined;
    const f = R.__call__(input);
    return [
      f.coeffs.map(String),
      f.degree(),
      String(f.parent),
      f.coeffs.every((a) => a.eq(R.base_ring.__call__(a))),
    ];
  });
functions.polynomial_gen_input = (base: bigint, kind: bigint, n: bigint) =>
  comparison(() => {
    const R = new PolynomialRing(polynomialBase(base), 'x');
    return String((R.gen as (...args: unknown[]) => unknown)(polynomialInput(kind, n)));
  });

functions.ff_prime__rational_ = (p: bigint, n: bigint, legacy: bigint) =>
  comparison(() => (legacy ? new FiniteFieldPrime(p) : new PrimeField(p)).__call__(n)._rational_());
functions.modular__rational_ = (p: bigint, n: bigint) =>
  comparison(() => Zmod(p).__call__(n)._rational_());
functions.polynomial_coefficient_conversion = (base: bigint, source: bigint, n: bigint) =>
  comparison(() => {
    const R = polynomialBase(base);
    const x =
      source === 0n
        ? Zmod(14n).__call__(n)
        : source === 1n
          ? new FiniteFieldPrime(7n).__call__(n)
          : source === 2n
            ? GFpn(7n, 2).fromInteger(n)
            : Zmod(5n).__call__(n);
    const y = R.__call__(x);
    return [String(y), String(R)];
  });

functions.polynomial_parent_constructor = (base: bigint, nested: bigint) =>
  comparison(() => {
    const R = new PolynomialRing(polynomialBase(base), 'y');
    const f = R.__call__([2n, 1n]);
    if (!nested) return R.__call__(f) === f;
    const S = new PolynomialRing(R, 'x');
    const g = S.__call__(f);
    return [g.degree(), g.coeffs.map(String), String(g.parent)];
  });

const polynomialIntegerBase = {
  zero: () => new Integer(0n),
  one: () => new Integer(1n),
  __call__: (x: unknown) => new Integer(x as bigint),
  is_field: () => false,
  toString: () => 'Integer Ring',
};
function unaryPolynomialBase(kind: bigint): CoefficientRing<RingElement> {
  if (kind === 13n || kind === 14n)
    return Zmod((1n << 127n) - (kind === 13n ? 2n : 1n)) as CoefficientRing<RingElement>;
  if (kind === 10n) return GF2 as CoefficientRing<RingElement>;
  return (
    kind === 6n
      ? polynomialIntegerBase
      : kind === 7n
        ? Zmod(1n)
        : kind === 8n
          ? new PrimeField((1n << 127n) - 1n)
          : kind === 9n
            ? GFpn(2n, 3, undefined, 'a')
            : polynomialBase(kind)
  ) as CoefficientRing<RingElement>;
}
functions.polynomial_monic = (base: bigint, coeffs: bigint[], denominator: bigint) =>
  comparison(() => {
    const R = new PolynomialRing(unaryPolynomialBase(base), 'x');
    const f = R.__call__(
      coeffs.map((c) => (denominator === 1n ? c : new Rational(c, denominator)))
    );
    const g = f.monic();
    return [g.coeffs.map(String), String(g.parent), g === f];
  });
functions.polynomial_is_monic = (base: bigint, coeffs: bigint[]) =>
  comparison(() => {
    return new PolynomialRing(unaryPolynomialBase(base), 'x').__call__(coeffs).is_monic();
  });
functions.polynomial_monic_extension = (p: bigint, degree: bigint, coeffs: bigint[]) =>
  comparison(() => {
    const F = GFpn(p, Number(degree), undefined, 'a');
    const R = new PolynomialRing(F, 'x');
    const f = R.__call__(coeffs.map((c) => F.fromInteger(c)));
    const g = f.monic();
    return [g.coeffs.map(String), String(g.parent), g === f];
  });

functions.gf2_zero_protocol = (operation: bigint, a: bigint, b: bigint) =>
  comparison(() => {
    const x = GF2.__call__(a);
    return String(
      operation === 0n
        ? x.div(GF2.__call__(b))
        : operation === 1n
          ? x.inv()
          : x.pow(operation === 3n ? Number(b) + 0.5 : b)
    );
  });

functions.polynomial_index_operation = (
  base: bigint,
  operation: bigint,
  coeffs: bigint[],
  kind: bigint,
  n: bigint
) =>
  comparison(() => {
    const B =
      base === 11n
        ? (new PolynomialRing(QQ, 'y') as CoefficientRing<RingElement>)
        : unaryPolynomialBase(base);
    const R = new PolynomialRing(B, 'x');
    const f = R.__call__(coeffs);
    const index =
      kind === 19n
        ? Number.POSITIVE_INFINITY
        : kind === 20n
          ? Number.NEGATIVE_INFINITY
          : polynomialInput(kind, n);
    if (operation === 3n) return String(f.getCoeff(index as number));
    const g =
      operation === 0n
        ? f.shift(index as number)
        : operation === 1n
          ? f.truncate(index as number)
          : f.reverse(index as number);
    return g == null ? null : [g.coeffs.map(String), String(g.parent), g === f];
  });

functions.polynomial_mixed_arithmetic = (
  left: bigint,
  right: bigint,
  operation: bigint,
  a: bigint[],
  b: bigint[],
  differentVariable: bigint,
  shape: bigint = 0n,
  variableKind: bigint = -1n,
  variableN: bigint = 0n
) =>
  comparison(() => {
    const base = (k: bigint) =>
      k === 17n
        ? RDF
        : k >= 12n && k <= 14n
          ? (Zmod([6n, 10n, 15n][Number(k - 12n)]!) as CoefficientRing<RingElement>)
          : k === 15n
            ? new PrimeField(5n)
            : k === 16n
              ? new PrimeField(11n)
              : k === 11n
                ? (new PolynomialRing(QQ, 't') as unknown as CoefficientRing<RingElement>)
                : unaryPolynomialBase(k);
    const R = new PolynomialRing(base(left), differentVariable === 3n ? 't' : 'x');
    const S = new PolynomialRing(
      base(right),
      differentVariable === 2n ? 't' : differentVariable === 1n ? 'y' : 'x'
    );
    const coeffs = (B: CoefficientRing<RingElement>, k: bigint, cs: bigint[]) =>
      cs.map((c) =>
        shape === 3n && (k === 4n || k === 9n)
          ? (B as FiniteFieldExtension).gen().add(c)
          : shape === 1n && k === 11n
            ? (B as PolynomialRing<RingElement>).gen().add(B.__call__(c) as never)
            : shape === 2n && (k === 0n || k === 11n)
              ? new Rational(c, 2n)
              : c
      );
    const f = R.__call__(coeffs(R.base_ring, left, a));
    const h = differentVariable === 4n ? f : S.__call__(coeffs(S.base_ring, right, b));
    if (operation === 11n) {
      const g = f.multiplication_trunc(h, variableN);
      return [g.coeffs.map(String), String(g.parent), g === f, g === h];
    }
    if (operation === 3n) return f.eq(h);
    if (operation === 10n) {
      const result = f.pseudo_quo_rem(h);
      return result.map((g) => [
        String(g),
        String(g.parent),
        g === f,
        g === h,
        result.map((k) => g === k),
      ]);
    }
    if (operation === 9n) {
      const g =
        variableKind === -1n
          ? f.resultant(h)
          : (f.resultant as (x: typeof h, options: { proof: unknown }) => RingElement)(h, {
              proof: polynomialInput(variableKind, variableN),
            });
      return polynomialScalarResult(g);
    }
    if (operation === 8n) {
      const variable =
        variableKind === 26n
          ? R.gen()
          : variableKind === 27n
            ? S.gen()
            : variableKind === 28n
              ? new PolynomialRing(QQ, 'z').gen()
              : polynomialInput(variableKind, variableN);
      const matrix =
        variableKind === -1n
          ? f.sylvester_matrix(h)
          : (f.sylvester_matrix as (h: typeof f, variable: unknown) => RingElement[][])(
              h,
              variable
            );
      return matrix.map((row) => row.map(String));
    }
    if (operation === 7n) {
      const result = f.xgcd(h);
      return result.map((g) => [
        g instanceof Polynomial ? g.coeffs.map(String) : String(g),
        g instanceof Polynomial ? String(g.parent) : 'Integer Ring',
        g === f,
        g === h,
        result.map((x) => g === x),
      ]);
    }
    if (operation === 4n)
      return f.quo_rem(h).map((g) => [g.coeffs.map(String), String(g.parent), g === f, g === h]);
    const g =
      operation === 6n
        ? f.gcd(h)
        : operation === 5n
          ? f.mod(h)
          : operation === 0n
            ? f.add(h)
            : operation === 1n
              ? f.sub(h)
              : f.mul(h);
    return [g.coeffs.map(String), String(g.parent), g === f, g === h];
  });

functions.polynomial_repr = (base: bigint, coeffs: bigint[], shape: bigint) =>
  comparison(() => {
    const B =
      base === 11n
        ? (new PolynomialRing(QQ, 't') as unknown as CoefficientRing<RingElement>)
        : unaryPolynomialBase(base);
    const R = new PolynomialRing(B, 'x');
    const values = coeffs.map((c, i) => {
      if (shape === 1n) return new Rational(c, 2n);
      if (shape === 2n) {
        const P = B as PolynomialRing<RingElement>;
        return P.gen()
          .pow(i + 1)
          .add(P.one())
          .scalar_mul(P.base_ring.__call__(c));
      }
      return c;
    });
    return String(R.__call__(values));
  });

functions.polynomial_scalar_equal = (
  base: bigint,
  coeffs: bigint[],
  kind: bigint,
  n: bigint,
  denominator: bigint = 1n
) =>
  comparison(() => {
    const B =
      base === 11n
        ? (new PolynomialRing(QQ, 't') as unknown as CoefficientRing<RingElement>)
        : unaryPolynomialBase(base);
    const scalar = kind === 19n ? Infinity : kind === 20n ? -Infinity : polynomialInput(kind, n);
    return new PolynomialRing(B, 'x')
      .__call__(coeffs.map((c) => (denominator === 1n ? c : new Rational(c, denominator))))
      .eq(scalar);
  });

functions.polynomial_sylvester_operand = (
  base: bigint,
  kind: bigint,
  n: bigint,
  coeffs: bigint[]
) =>
  comparison(() => {
    const R = new PolynomialRing(
      base === 11n
        ? (new PolynomialRing(QQ, 't') as unknown as CoefficientRing<RingElement>)
        : unaryPolynomialBase(base),
      'x'
    );
    const f = R.__call__(coeffs);
    return (f.sylvester_matrix as (x: unknown) => RingElement[][])(polynomialInput(kind, n)).map(
      (row) => row.map(String)
    );
  });

functions.polynomial_constructor_scalar = (base: bigint, kind: bigint, n: bigint) =>
  comparison(() => {
    const R = new PolynomialRing(
      base === 11n
        ? (new PolynomialRing(QQ, 't') as unknown as CoefficientRing<RingElement>)
        : unaryPolynomialBase(base),
      'x'
    );
    const f = R.__call__(polynomialInput(kind, n));
    return [f.degree(), f.coeffs.map(String), String(f.parent)];
  });

functions.polynomial_sylvester_degree_one = (p: bigint, a: bigint[], b: bigint[]) =>
  comparison(() => {
    const R = new PolynomialRing(GFpn(p, 1), 'x');
    return R.__call__(a)
      .sylvester_matrix(R.__call__(b))
      .map((row) => row.map(String));
  });

function polynomialScalarResult(g: unknown): unknown[] {
  return [
    String(g),
    g instanceof Rational
      ? 'Rational Field'
      : g instanceof Integer || typeof g === 'bigint'
        ? 'Integer Ring'
        : String((g as { parent: unknown }).parent),
  ];
}
functions.polynomial_resultant_operand = (
  base: bigint,
  kind: bigint,
  n: bigint,
  coeffs: bigint[]
) =>
  comparison(() => {
    const R = new PolynomialRing(
      base === 11n
        ? (new PolynomialRing(QQ, 't') as unknown as CoefficientRing<RingElement>)
        : unaryPolynomialBase(base),
      'x'
    );
    const f = R.__call__(coeffs);
    return polynomialScalarResult(
      (f.resultant as (x: unknown) => RingElement)(polynomialInput(kind, n))
    );
  });

import { RDF, RealDoubleElement } from '../../../../packages/sagemath-ts/src/rings/real_double.js';
function doubleFromBits(bits: bigint): number {
  const v = new DataView(new ArrayBuffer(8));
  v.setBigUint64(0, bits, false);
  return v.getFloat64(0, false);
}
function doubleBits(x: number): string {
  if (Number.isNaN(x)) return 'NaN';
  const v = new DataView(new ArrayBuffer(8));
  v.setFloat64(0, x, false);
  return String(v.getBigUint64(0, false));
}
functions.rdf_coefficient = (operation: bigint, a: bigint, b: bigint) =>
  comparison(() => {
    const x = RDF.__call__(doubleFromBits(a)),
      y = RDF.__call__(doubleFromBits(b));
    if (operation === 6n) return x.eq(y);
    if (operation === 9n)
      return [
        String(RDF),
        RDF.characteristic(),
        RDF.is_field(),
        RDF.zero() === RDF.zero(),
        RDF.one() === RDF.one(),
        RDF.__call__(x) === x,
      ];
    const r =
      operation === 0n
        ? x
        : operation === 1n
          ? x.add(y)
          : operation === 2n
            ? x.sub(y)
            : operation === 3n
              ? x.mul(y)
              : operation === 4n
                ? x.div(y)
                : operation === 5n
                  ? x.neg()
                  : operation === 7n
                    ? RDF.zero()
                    : operation === 8n
                      ? RDF.one()
                      : operation === 10n
                        ? RDF.__call__(new Rational(a, b))
                        : operation === 11n
                          ? RDF.__call__(a)
                          : RDF.__call__(Boolean(a));
    return [doubleBits(r.value), String(r), r.isZero(), String(r.parent)];
  });

functions.rdf_polynomial_resultant = (a: bigint[], b: bigint[]) =>
  comparison(() => {
    const R = new PolynomialRing(RDF, 'x');
    const r = R.__call__(a.map((x) => RDF.__call__(doubleFromBits(x)))).resultant(
      R.__call__(b.map((x) => RDF.__call__(doubleFromBits(x))))
    );
    return [doubleBits(r.value), String(r), String(r.parent)];
  });

function derivativeArguments(R: PolynomialRing<RingElement>, form: bigint, n: bigint): unknown[] {
  const x = R.gen(),
    B = R.base_ring as CoefficientRing<RingElement> & { gen?: () => unknown };
  const z = B.gen ? B.gen() : B.__call__(1n);
  return form === 0n
    ? []
    : form === 1n
      ? [n]
      : form === 2n
        ? [new Integer(n)]
        : form === 3n
          ? [new Rational(n, 2n)]
          : form === 4n
            ? [Boolean(n)]
            : form === 5n
              ? [String(n)]
              : form === 6n
                ? [null]
                : form === 7n
                  ? [x]
                  : form === 8n
                    ? [R.__call__(n).mul(x)]
                    : form === 9n
                      ? [new PolynomialRing(B, 'y').gen()]
                      : form === 10n
                        ? [new PolynomialRing(QQ, 'x').gen()]
                        : form === 11n
                          ? [z]
                          : form === 12n
                            ? [[]]
                            : form === 13n
                              ? [[x]]
                              : form === 14n
                                ? [[x, n]]
                                : form === 15n
                                  ? [x, n]
                                  : form === 16n
                                    ? [x, n, n]
                                    : form === 17n
                                      ? [n, x, n]
                                      : form === 18n
                                        ? [n, null, n]
                                        : form === 19n
                                          ? [[null, x]]
                                          : form === 20n
                                            ? [[null, []]]
                                            : form === 21n
                                              ? [Number.NaN]
                                              : [Number(n) + 0.5];
}
functions.polynomial_derivative_protocol = (
  base: bigint,
  method: bigint,
  form: bigint,
  n: bigint,
  coeffs: bigint[],
  shape: bigint
) =>
  comparison(() => {
    const B = (
      base === 17n
        ? RDF
        : base === 18n
          ? new PolynomialRing(
              polynomialIntegerBase as unknown as CoefficientRing<RingElement>,
              't'
            )
          : base === 19n
            ? new PolynomialRing(RDF, 't')
            : base === 11n
              ? new PolynomialRing(QQ, 't')
              : unaryPolynomialBase(base)
    ) as CoefficientRing<RingElement> & { gen: () => RingElement };
    const R = new PolynomialRing(B, 'x');
    const values = coeffs.map((c) =>
      shape === 3n
        ? RDF.__call__(doubleFromBits(c))
        : shape === 2n
          ? B.gen().add(B.__call__(c))
          : shape === 1n
            ? new Rational(c, 2n)
            : c
    );
    const f = R.__call__(values);
    const methods = f as unknown as Record<
      string,
      (...args: unknown[]) => Polynomial<RingElement> | Polynomial<RingElement>[]
    >;
    if (method === 5n)
      return [methods.diff === methods.derivative, methods.differentiate === methods.derivative];
    const args = derivativeArguments(R, form, n);
    const result =
      method === 6n
        ? multi_derivative(f, args)
        : methods[
            ['derivative', '_derivative', 'diff', 'differentiate', 'gradient'][Number(method)]!
          ]!(...args);
    const frame = (g: Polynomial<RingElement>) => [g.coeffs.map(String), String(g.parent), g === f];
    return Array.isArray(result) ? result.map(frame) : frame(result);
  });

import {
  derivative_parse,
  multi_derivative,
} from '../../../../packages/sagemath-ts/src/misc/derivative.js';
functions.derivative_parse_protocol = (form: bigint, n: bigint) =>
  comparison(() => {
    const R = new PolynomialRing(QQ, 'x') as unknown as PolynomialRing<RingElement>;
    const args = derivativeArguments(R, form, n),
      out = derivative_parse(args);
    const frame = (x: unknown): unknown =>
      Array.isArray(x)
        ? x.map(frame)
        : x == null
          ? null
          : typeof x === 'number' && Number.isNaN(x)
            ? 'nan'
            : String(x);
    return [frame(out), args.length === 1 && out === args[0]];
  });

functions.polynomial_pseudo_operand = (base: bigint, kind: bigint, n: bigint, coeffs: bigint[]) =>
  comparison(() => {
    const B = (
      base === 11n ? new PolynomialRing(QQ, 't') : unaryPolynomialBase(base)
    ) as CoefficientRing<RingElement>;
    const R = new PolynomialRing(B, 'x'),
      f = R.__call__(coeffs),
      h = polynomialInput(kind, n);
    const result = f.pseudo_quo_rem(h as Polynomial<RingElement>);
    return result.map((g) => [
      String(g),
      String(g.parent),
      g === f,
      g === h,
      result.map((k) => g === k),
    ]);
  });

import {
  FractionField,
  FractionField_generic,
} from '../../../../packages/sagemath-ts/src/rings/fraction_field.js';
functions.polynomial_constant_fraction = (
  base: bigint,
  a: bigint[],
  b: bigint[],
  d: bigint,
  op: bigint
) =>
  comparison(() => {
    const R = new PolynomialRing(unaryPolynomialBase(base), 'x'),
      K = new FractionField_generic(R);
    const f = K.__call__(a, d),
      g = K.__call__(b, d);
    if (op === 6n) return f.eq(g);
    const h =
      op === 1n
        ? f.add(g)
        : op === 2n
          ? f.sub(g)
          : op === 3n
            ? f.mul(g)
            : op === 4n
              ? f.div(g)
              : op === 5n
                ? f.neg()
                : f;
    return [
      String(h),
      String(h.numerator()),
      String(h.denominator()),
      String(h.parent),
      h.isZero(),
      K.is_field(),
      K.ring() === R,
      String(K.zero()),
      String(K.one()),
      K.__call__(f) === f,
    ];
  });

import {
  FractionFieldElement,
  type FractionElement,
} from '../../../../packages/sagemath-ts/src/rings/fraction_field_element.js';
functions.polynomial_fraction_protocol = (
  base: bigint,
  op: bigint,
  a: bigint[],
  ad: bigint[],
  b: bigint[],
  bd: bigint[],
  flavor: bigint,
  mode: bigint,
  shape: bigint = 0n
) =>
  comparison(() => {
    const B = unaryPolynomialBase(base),
      R = new PolynomialRing(B, 'x'),
      K = flavor ? new FractionField_generic(R) : R.fraction_field();
    const coefficients = (values: bigint[]) =>
      values.map((c) =>
        shape && (base === 4n || base === 9n)
          ? (B as unknown as { gen(): { add(c: bigint): RingElement } }).gen().add(c)
          : shape && base === 0n
            ? new Rational(c, 2n)
            : c
      );
    const [A, D, E, F] = [a, ad, b, bd].map((v) => R.__call__(coefficients(v)));
    const f =
      mode === 0n ? K.__call__(A, D) : new K._element_class(K, A, D, { reduce: mode === 1n });
    const g =
      mode === 0n ? K.__call__(E, F) : new K._element_class(K, E, F, { reduce: mode === 1n });
    if (op === 6n) return f.eq(g);
    const protocol = f as unknown as {
      inv(): FractionFieldElement<RingElement>;
      pow(n: bigint): FractionFieldElement<RingElement>;
      reduce(): void;
    };
    let h: FractionElement<RingElement>;
    if (op === 11n) {
      if (typeof protocol.reduce !== 'function')
        throw new AttributeError(
          "'sage.rings.fraction_field_FpT.FpTElement' object has no attribute 'reduce'"
        );
      protocol.reduce();
      h = f;
    } else
      h =
        op === 1n
          ? f.add(g)
          : op === 2n
            ? f.sub(g)
            : op === 3n
              ? f.mul(g)
              : op === 4n
                ? f.div(g)
                : op === 5n
                  ? f.neg()
                  : op === 7n
                    ? protocol.inv()
                    : op === 8n
                      ? protocol.pow(2n)
                      : op === 9n
                        ? protocol.pow(0n)
                        : op === 10n
                          ? protocol.pow(-2n)
                          : f;
    const n = h.numerator(),
      d = h.denominator();
    return [
      String(h),
      n.coeffs.map(String),
      d.coeffs.map(String),
      String(h.parent),
      h.constructor.name,
      h.parent === K,
      h === f,
      h === g,
      n === h.numerator(),
      d === h.denominator(),
      K.zero() === K.zero(),
      K.one() === K.one(),
      h.isZero(),
      String(K.base_ring()),
      String(K.characteristic()),
      K.is_exact(),
      String(K.gen()),
      K.ngens(),
      R.fraction_field() === R.fraction_field(),
      FractionField(R) === R.fraction_field(),
      h.constructor.name !== 'FpTElement' ||
        ((h as FpTElement<RingElement>).numer().eq(n) &&
          (h as FpTElement<RingElement>).denom().eq(d)),
    ];
  });

import { FpT, FpTElement } from '../../../../packages/sagemath-ts/src/rings/fraction_field_FpT.js';
functions.polynomial_fraction_parent = (base: bigint, flavor: bigint, index: bigint) =>
  comparison(() => {
    const R = new PolynomialRing(unaryPolynomialBase(base), 'x'),
      K =
        flavor === 2n
          ? new FpT(R)
          : flavor === 1n
            ? new FractionField_generic(R)
            : FractionField(R);
    return [
      String(K),
      K.constructor.name,
      String(K.gen(index)),
      K.ngens(),
      String(K.base_ring()),
      String(K.characteristic()),
      K.is_exact(),
    ];
  });

import {
  FreeModule,
  type FreeModulePID,
} from '../../../../packages/sagemath-ts/src/modules/free_module.js';
functions.polynomial_fraction_span = (base: bigint, n: bigint, c: bigint) =>
  comparison(() => {
    const R = new PolynomialRing(unaryPolynomialBase(base), 'x'),
      x = R.gen(),
      A = FreeModule(R as never, 3) as FreeModulePID;
    const L = A.span([
      A.createElement([x, x.pow(n).add(R.__call__(c)), R.zero()]),
      A.createElement([R.zero(), R.zero(), x]),
    ]);
    const V = L.vectorSpaceSpan(L.basis());
    return [
      String(A.baseField()),
      V.dimension(),
      (V.basisMatrix() as unknown[][]).map((row) => row.map(String)),
    ];
  });

functions.polynomial_fraction_conversion = (
  source: bigint,
  target: bigint,
  mode: bigint,
  operation: bigint,
  a: bigint[],
  d: bigint[],
  variable: bigint = 0n
) =>
  comparison(() => {
    const B = unaryPolynomialBase(source),
      R = new PolynomialRing(B, 'x'),
      T =
        source === target && !variable
          ? R
          : new PolynomialRing(unaryPolynomialBase(target), variable ? 'y' : 'x');
    const K = mode === 1n || mode === 2n ? new FractionField_generic(R) : R.fraction_field();
    const f = new K._element_class(K, R.__call__(a), R.__call__(d), {
      reduce: mode === 0n || mode === 1n,
    });
    const n = f.numerator(),
      den = f.denominator();
    let out: unknown[];
    try {
      const h =
        operation === 0n
          ? T.__call__(f)
          : operation === 1n
            ? T.base_ring.__call__(f)
            : operation === 3n
              ? (f as FractionFieldElement<RingElement>)._conversion(
                  target === 6n ? ZZ : T.base_ring
                )
              : T.fraction_field().__call__(f);
      out = [
        'value',
        String(h),
        String(
          operation === 1n || operation === 3n
            ? T.base_ring
            : (h as unknown as { parent: unknown }).parent
        ),
        h === n,
      ];
    } catch (e) {
      out = ['error', (e as Error).name, (e as Error).message];
    }
    return [
      out,
      String(f),
      f.numerator().coeffs.map(String),
      f.denominator().coeffs.map(String),
      n === f.numerator(),
      den === f.denominator(),
    ];
  });

functions.polynomial_fraction_import = (mode: bigint) =>
  comparison(() => {
    const spec = [
      ['fraction_field', ['FractionField', 'FractionField_generic', 'FractionField_1poly_field']],
      ['fraction_field_element', ['FractionFieldElement', 'FractionFieldElement_1poly_field']],
      ['fraction_field_FpT', ['FpT', 'FpTElement']],
    ][Number(mode)]!;
    const [file, names] = spec as [string, string[]];
    const path = new URL(`../../../../packages/sagemath-ts/src/rings/${file}.ts`, import.meta.url)
      .pathname;
    const code = `import * as m from ${JSON.stringify(path)};console.log(JSON.stringify(${JSON.stringify(names)}.map(k=>typeof m[k]==='function')));`;
    const child = Bun.spawnSync([process.execPath, '-e', code]);
    if (child.exitCode !== 0) throw new Error(`isolated ${file} import failed`);
    return JSON.parse(child.stdout.toString());
  });

functions.polynomial_ring_coercion = (source: bigint, target: bigint, variable: bigint) =>
  comparison(() =>
    new PolynomialRing(unaryPolynomialBase(target), variable ? 'y' : 'x').has_coerce_map_from(
      new PolynomialRing(unaryPolynomialBase(source), 'x')
    )
  );

functions.polynomial_fraction_tower = (
  base: bigint,
  mode: bigint,
  operation: bigint,
  a: bigint[],
  d: bigint[]
) =>
  comparison(() => {
    const U = new PolynomialRing(unaryPolynomialBase(base), 't'),
      F = mode === 1n || mode === 2n ? new FractionField_generic(U) : U.fraction_field();
    const f = new F._element_class(F, U.__call__(a), U.__call__(d), {
      reduce: mode === 0n || mode === 1n,
    });
    const n = f.numerator(),
      den = f.denominator(),
      R = new PolynomialRing(
        (operation === 1n || operation === 3n || operation === 5n || operation === 7n
          ? U
          : F) as CoefficientRing<RingElement>,
        'x'
      );
    let out: unknown[];
    try {
      let result: unknown;
      if (operation < 2n) {
        const h = R.__call__(f);
        result = [String(h), String(h.parent), h.coeffs.map(String)];
      } else if (operation === 10n) {
        result = [
          R.__call__(f).eq(f),
          R.gen().eq(f),
          U.gen().eq(f),
          R.__call__(f).eq(R.__call__(f)),
          f.eq(U.gen()),
        ];
      } else if (operation >= 6n) {
        const z = operation === 6n || operation === 7n ? R.gen() : R.__call__(f);
        const h = operation === 6n || operation === 8n ? F.__call__(z) : U.__call__(z);
        result = [String(h), String(h.parent), h === f];
      } else {
        const g = R.__call__(operation < 4n ? [U.gen()] : [U.one(), U.zero(), U.gen()]);
        const values = (operation < 4n ? R.zero() : R.one()).pseudo_quo_rem(g);
        result = values.map((v) => [String(v), String(v.parent), values.map((w) => v === w)]);
      }
      out = ['value', result];
    } catch (e) {
      out = ['error', (e as Error).name, (e as Error).message];
    }
    return [
      out,
      String(f),
      f.numerator().coeffs.map(String),
      f.denominator().coeffs.map(String),
      n === f.numerator(),
      den === f.denominator(),
    ];
  });

functions.polynomial_scalar_conversion = (
  source: bigint,
  target: bigint,
  operation: bigint,
  coeffs: bigint[],
  denominator: bigint
) =>
  comparison(() => {
    const R = new PolynomialRing(unaryPolynomialBase(source), 'x');
    const f = R.__call__(
      coeffs.map((c) => (denominator === 1n ? c : new Rational(c, denominator)))
    );
    const T = target === 6n ? ZZ : unaryPolynomialBase(target);
    const h =
      operation === 0n
        ? T.__call__(f)
        : operation === 1n
          ? f._scalar_conversion(T)
          : operation === 2n
            ? f._integer_(ZZ)
            : f._rational_();
    return [String(h), String(operation === 2n ? ZZ : operation === 3n ? QQ : T)];
  });
functions.gf2_scalar_hook = (value: bigint, operation: bigint) =>
  comparison(() =>
    String(operation === 0n ? GF2.__call__(value)._integer_() : GF2.__call__(value)._rational_())
  );

functions.polynomial_power_integer = (
  base: bigint,
  kind: bigint,
  exponent: bigint,
  coeffs: bigint[]
) =>
  comparison(() => {
    const R = new PolynomialRing(unaryPolynomialBase(base), 'x');
    const f = R.__call__(coeffs);
    const power = kind === 0n ? exponent : kind === 1n ? Number(exponent) : new Integer(exponent);
    const h = f.pow(power as bigint);
    return [
      String(h),
      String(h.parent),
      h === f,
      h instanceof Polynomial
        ? ['polynomial', h.coeffs.map(String)]
        : [
            'fraction',
            (h as FractionElement<RingElement>).numerator().coeffs.map(String),
            (h as FractionElement<RingElement>).denominator().coeffs.map(String),
          ],
    ];
  });

functions.polynomial_power_strict_rational = (
  coeffs: bigint[],
  denominator: bigint,
  exponent: bigint
) =>
  comparison(() => {
    // The coefficient-ring contract allows integer-only construction plus division.
    const base: CoefficientRing<Rational> = {
      zero: () => new Rational(0n),
      one: () => new Rational(1n),
      __call__: (x: unknown) => {
        if (typeof x !== 'bigint') throw new Error('expected an integer');
        return new Rational(x);
      },
      is_field: () => true,
      toString: () => 'Rational Field',
    };
    const R = new PolynomialRing(base, 'x');
    const f = new Polynomial(
      coeffs.map((c) => new Rational(c, denominator)),
      R
    );
    const h = f.pow(exponent) as Polynomial<Rational>;
    return [String(h), h.coeffs.map(String), h.coeffs.every((c) => c instanceof Rational)];
  });

functions.polynomial_power_exponent = (
  base: bigint,
  kind: bigint,
  numerator: bigint,
  denominator: bigint,
  coeffs: bigint[]
) =>
  comparison(() => {
    const R = new PolynomialRing(unaryPolynomialBase(base), 'x'),
      f = R.__call__(coeffs);
    const exponent: unknown =
      kind === 0n
        ? numerator
        : kind === 1n
          ? Number(numerator) / Number(denominator)
          : kind === 2n
            ? new Rational(numerator, denominator)
            : kind === 3n
              ? Boolean(numerator)
              : kind === 4n
                ? String(numerator)
                : kind === 5n
                  ? null
                  : kind === 6n
                    ? [numerator]
                    : kind === 7n
                      ? R.__call__(numerator)
                      : kind === 8n
                        ? new FiniteFieldPrime(7n).__call__(numerator)
                        : kind === 9n
                          ? GF2.__call__(numerator)
                          : kind === 10n
                            ? RDF.__call__(Number(numerator) / Number(denominator))
                            : kind === 11n
                              ? Number.NaN
                              : kind === 12n
                                ? Number.POSITIVE_INFINITY
                                : kind === 13n
                                  ? Number.NEGATIVE_INFINITY
                                  : R.__call__([numerator, denominator]);
    const h = (f.pow as (n: unknown) => Polynomial<RingElement> | FractionElement<RingElement>)(
      exponent
    );
    return [
      String(h),
      String(h.parent),
      h === f,
      h instanceof Polynomial
        ? ['polynomial', h.coeffs.map(String)]
        : ['fraction', h.numerator().coeffs.map(String), h.denominator().coeffs.map(String)],
    ];
  });

functions.polynomial_root_series = (
  base: bigint,
  operation: bigint,
  coeffs: bigint[],
  denominator: bigint,
  n: bigint,
  precision: bigint,
  other: bigint[],
  startKind: bigint
) =>
  comparison(() => {
    const R = new PolynomialRing(unaryPolynomialBase(base), 'x');
    const f = R.__call__(
        coeffs.map((c) => (denominator === 1n ? c : new Rational(c, denominator)))
      ),
      g = R.__call__(other);
    const h =
      operation === 7n
        ? generic_power_trunc(f, n, Number(precision))
        : operation === 0n
          ? f.nth_root(n)
          : operation === 1n
            ? f._nth_root_series(
                n,
                precision,
                startKind === 0n ? undefined : startKind === 1n ? R.one() : g
              )
            : operation === 2n
              ? f.inverse_series_trunc(precision)
              : operation === 3n
                ? f.power_trunc(n, precision)
                : operation === 4n
                  ? f._power_trunc(n, precision)
                  : operation === 5n
                    ? f._mul_trunc_(g, precision)
                    : f.multiplication_trunc(g, precision);
    return [String(h), String(h.parent), h === f, h === g, h.coeffs.map(String)];
  });

import { _coefficient_nth_root } from '../../../../packages/sagemath-ts/src/rings/finite_rings/element_base.js';
functions.finite_coefficient_nth_root = (
  kind: bigint,
  p: bigint,
  degree: bigint,
  n: bigint,
  value: bigint
) =>
  comparison(() => {
    const B =
      kind === 0n
        ? Zmod(p)
        : kind === 1n
          ? new PrimeField(p)
          : kind === 2n
            ? new FiniteFieldPrime(p)
            : kind === 3n
              ? GF2
              : GFpn(p, Number(degree), undefined, 'a');
    const c = kind === 4n ? (B as ReturnType<typeof GFpn>).fromInteger(value) : B.__call__(value);
    const h = _coefficient_nth_root(c as RingElement, n);
    return [String(h), String(B)];
  });

import { generic_power_trunc } from '../../../../packages/sagemath-ts/src/rings/polynomial/polynomial_element.js';

functions.polynomial_factor_order = (base: bigint, factors: bigint[], multiplicities: bigint[]) =>
  comparison(() => {
    const R = new PolynomialRing(unaryPolynomialBase(base), 'x');
    let f = R.one();
    multiplicities.forEach((m, j) => {
      f = f.mul(R.__call__(factors.slice(3 * j, 3 * j + 3)).pow(m) as Polynomial<RingElement>);
    });
    return f.factor().map(([g, e]) => [String(g), String(e)]);
  });

functions.polynomial_series_input = (
  base: bigint,
  operation: bigint,
  slot: bigint,
  kind: bigint,
  value: bigint,
  coeffs: bigint[]
) =>
  comparison(() => {
    const R = new PolynomialRing(unaryPolynomialBase(base), 'x'),
      f = R.__call__(coeffs);
    const v = polynomialInput(kind, value),
      n = slot === 0n ? v : 2n,
      precision = slot === 1n ? v : 3n;
    const h =
      operation === 0n
        ? f.nth_root(n)
        : operation === 1n
          ? f._nth_root_series(n, precision)
          : operation === 2n
            ? f.inverse_series_trunc(precision)
            : operation === 3n
              ? f.power_trunc(n, precision)
              : operation === 4n
                ? f._power_trunc(n, precision)
                : operation === 5n
                  ? f._mul_trunc_(f, precision)
                  : f.multiplication_trunc(f, precision);
    return [String(h), String(h.parent), h === f];
  });

functions.polynomial_truncated_operand = (
  base: bigint,
  kind: bigint,
  value: bigint,
  coeffs: bigint[],
  precision: bigint,
  internal: bigint,
  precisionKind: bigint = -1n
) =>
  comparison(() => {
    const R = new PolynomialRing(base === 11n ? RDF : unaryPolynomialBase(base), 'x'),
      f = R.__call__(coeffs);
    const other =
      kind === 27n
        ? RDF.__call__(Number(value) / 2)
        : kind === 28n
          ? R.__call__(value)
          : kind === 29n
            ? R.fraction_field().__call__(1n, R.gen().add(R.one()))
            : kind === 30n
              ? R.__call__([value, 1n])
              : polynomialInput(kind, value);
    const n = precisionKind === -1n ? precision : polynomialInput(precisionKind, precision);
    const g = internal
      ? f._mul_trunc_(other as Polynomial<RingElement>, n)
      : f.multiplication_trunc(other, n);
    return [String(g), String(g.parent), g === f, g === other];
  });

functions.polynomial_real_truncated_product = (
  a: bigint[],
  b: bigint[],
  n: bigint,
  publicMethod: bigint
) =>
  comparison(() => {
    const view = new DataView(new ArrayBuffer(8));
    const unpack = (x: bigint) => {
      view.setBigUint64(0, x);
      return RDF.__call__(view.getFloat64(0));
    };
    const pack = (x: RingElement) => {
      const value = (x as RealDoubleElement).value;
      if (Number.isNaN(value)) return 'nan';
      view.setFloat64(0, value);
      return view.getBigUint64(0).toString(16).padStart(16, '0');
    };
    const R = new PolynomialRing(RDF, 'x'),
      f = R.__call__(a.map(unpack)),
      g = [3n, 4n, 5n].includes(publicMethod) ? f : R.__call__(b.map(unpack));
    const h = [2n, 3n].includes(publicMethod)
      ? f.mul(g)
      : [1n, 5n].includes(publicMethod)
        ? f.multiplication_trunc(g, n)
        : f._mul_trunc_(g, n);
    return [h.coeffs.map(pack), String(h)];
  });

import { createHash } from 'node:crypto';
functions.polynomial_exact_product = (
  base: bigint,
  a: bigint[],
  b: bigint[],
  same: bigint,
  shape: bigint,
  denominator: bigint
) =>
  comparison(() => {
    const inner = new PolynomialRing(QQ, 't');
    const B =
      base === 11n ? inner : base === 12n ? inner.fraction_field() : unaryPolynomialBase(base);
    const R = new PolynomialRing(B, 'x');
    const coefficient = (c: bigint): unknown => {
      if ((base === 4n || base === 9n) && shape) {
        const p = (B as FiniteFieldExtension).characteristic,
          v: bigint[] = [];
        while (c) {
          v.push(c % p);
          c /= p;
        }
        return B.__call__(v);
      }
      if (base === 11n || base === 12n)
        return shape
          ? (B.__call__(inner.gen()) as RingElement).add(
              B.__call__(new Rational(c, denominator)) as never
            )
          : B.__call__(new Rational(c, denominator));
      return base === 0n ? new Rational(c, denominator) : c;
    };
    const f = R.__call__(a.map(coefficient)),
      g = same ? f : R.__call__(b.map(coefficient)),
      h = f.mul(g);
    return [
      h.coeffs.length,
      createHash('sha256').update(h.coeffs.map(String).join(',')).digest('hex'),
      String(h.parent),
      h === f,
      h === g,
    ];
  });

import { IntegerMatrix } from '../../../../packages/sagemath-ts/src/matrix/matrix_integer.js';
functions.polynomial_matrix_product = (a: bigint[], b: bigint[], same: bigint) =>
  comparison(() => {
    // Protocol adapter only: all matrix arithmetic uses production IntegerMatrix.
    class Coefficient implements RingElement {
      constructor(readonly matrix: IntegerMatrix) {}
      add(b: this): this {
        return new Coefficient(this.matrix.add(b.matrix)) as this;
      }
      sub(b: this): this {
        return new Coefficient(this.matrix.sub(b.matrix)) as this;
      }
      mul(b: this): this {
        return new Coefficient(this.matrix.mul(b.matrix)) as this;
      }
      neg(): this {
        return new Coefficient(this.matrix.neg()) as this;
      }
      eq(b: this | number): boolean {
        return typeof b === 'number'
          ? this.matrix.eq(
              new IntegerMatrix(2, 2, [
                [BigInt(b), 0n],
                [0n, BigInt(b)],
              ])
            )
          : this.matrix.eq(b.matrix);
      }
      isZero(): boolean {
        return this.eq(0);
      }
      toString(): string {
        return this.matrix.toString();
      }
    }
    const base = {
      zero: () => new Coefficient(new IntegerMatrix(2, 2)),
      one: () => base.__call__(1n),
      __call__: (x: unknown): Coefficient =>
        x instanceof Coefficient
          ? x
          : new Coefficient(
              new IntegerMatrix(2, 2, [
                [BigInt(x as bigint), 0n],
                [0n, BigInt(x as bigint)],
              ])
            ),
      is_field: () => false,
    };
    const R = new PolynomialRing(base, 'x');
    const poly = (v: bigint[]) =>
      R.__call__(
        Array.from(
          { length: v.length / 4 },
          (_, i) =>
            new Coefficient(
              new IntegerMatrix(2, 2, [v.slice(i * 4, i * 4 + 2), v.slice(i * 4 + 2, i * 4 + 4)])
            )
        )
      );
    const f = poly(a),
      g = same ? f : poly(b);
    return f
      .mul(g)
      .coeffs.map((c) =>
        [c.matrix.get(0, 0), c.matrix.get(0, 1), c.matrix.get(1, 0), c.matrix.get(1, 1)].map(String)
      );
  });

functions.polynomial_modular_power = (
  base: bigint,
  a: bigint[],
  kind: bigint,
  value: bigint,
  m: bigint[],
  exponentKind: bigint,
  exponent: bigint
) =>
  comparison(() => {
    const B =
      base === 11n ? new PolynomialRing(QQ, 't') : base === 12n ? RDF : unaryPolynomialBase(base);
    const R = new PolynomialRing(B, 'x'),
      f = R.__call__(a);
    const modulus =
      kind === 29n
        ? new PolynomialRing(QQ, 'x').__call__(m)
        : kind === 30n
          ? new PolynomialRing(new PrimeField(5n), 'x').__call__(m)
          : kind === 31n
            ? new PolynomialRing(B, 'y').__call__(m)
            : kind === 32n
              ? new PolynomialRing(new PolynomialRing(QQ, 't'), 'x').__call__(m)
              : kind === 0n
                ? R.__call__(m)
                : kind === 1n
                  ? f
                  : kind >= 2n
                    ? polynomialInput(kind - 2n, value)
                    : null;
    const e = polynomialInput(exponentKind, exponent);
    const h = (
      f.pow as (n: unknown, m: unknown) => Polynomial<RingElement> | FractionElement<RingElement>
    )(e, modulus);
    return [String(h), String(h.parent), h === f, h === modulus];
  });

import { gsl_pow_int, gsl_sf_log, gsl_sf_exp } from '../../../../packages/gsl-ts/src/index.js';
functions.polynomial_real_integer_power = (bits: bigint, n: bigint, method: bigint) =>
  comparison(() => {
    const v = new DataView(new ArrayBuffer(8));
    v.setBigUint64(0, bits);
    const x = RDF.__call__(v.getFloat64(0));
    const h =
      method === 0n
        ? x.pow(n).value
        : method === 1n
          ? x.inv().value
          : method === 2n
            ? (
                new PolynomialRing(RDF, 'x').__call__([x]).pow(n) as Polynomial<RealDoubleElement>
              ).getCoeff(0).value
            : method === 3n
              ? gsl_pow_int(x.value, Number(n))
              : method === 4n
                ? gsl_sf_log(x.value)
                : gsl_sf_exp(x.value);
    if (Number.isNaN(h)) return 'nan';
    v.setFloat64(0, h);
    return v.getBigUint64(0).toString(16).padStart(16, '0');
  });

functions.polynomial_real_quotient = (a: bigint[], b: bigint[], same: bigint) =>
  comparison(() => {
    const v = new DataView(new ArrayBuffer(8));
    const unpack = (x: bigint) => {
      v.setBigUint64(0, x);
      return RDF.__call__(v.getFloat64(0));
    };
    const pack = (x: RealDoubleElement) => {
      if (Number.isNaN(x.value)) return 'nan';
      v.setFloat64(0, x.value);
      return v.getBigUint64(0).toString(16).padStart(16, '0');
    };
    const R = new PolynomialRing(RDF, 'x'),
      f = R.__call__(a.map(unpack)),
      g = same ? f : R.__call__(b.map(unpack));
    return f.quo_rem(g).map((h) => [h.coeffs.map(pack), String(h)]);
  });

Object.assign(functions, {
  polynomial_real_evaluate: (a: bigint[], point: bigint) =>
    comparison(() => {
      const unpack = (v: bigint) => {
        const b = Buffer.alloc(8);
        b.writeBigUInt64BE(v);
        return RDF.__call__(b.readDoubleBE());
      };
      const value = new PolynomialRing(RDF, 'x')
        .__call__(a.map(unpack))
        .evaluate(unpack(point)).value;
      if (Number.isNaN(value)) return 'nan';
      const b = Buffer.alloc(8);
      b.writeDoubleBE(value);
      return b.toString('hex');
    }),
});

import { CompiledPolynomialFunction } from '../../../../packages/sagemath-ts/src/rings/polynomial/polynomial_compiled.js';
Object.assign(functions, {
  polynomial_compiled_evaluation: (a: bigint[], point: bigint, algorithm: bigint, method: bigint) =>
    comparison(() => {
      const unpack = (v: bigint) => {
        const b = Buffer.alloc(8);
        b.writeBigUInt64BE(v);
        return RDF.__call__(b.readDoubleBE());
      };
      const pack = (v: ReturnType<typeof unpack>) => {
        if (Number.isNaN(v.value)) return 'nan';
        const b = Buffer.alloc(8);
        b.writeDoubleBE(v.value);
        return b.toString('hex');
      };
      const f = new CompiledPolynomialFunction(
          a.map(unpack),
          ['binary', 'pippenger', 'other', ''][Number(algorithm)]!
        ),
        x = unpack(point);
      return method === 1n
        ? f.toString()
        : method === 2n
          ? [pack(f.__call__(x)), pack(f.__call__(x.neg())), pack(f.__call__(x))]
          : pack(f.eval(x));
    }),
});

Object.assign(functions, {
  polynomial_evaluation_input: (base: bigint, a: bigint[], kind: bigint, n: bigint) =>
    comparison(() => {
      const B =
        base === 11n ? new PolynomialRing(QQ, 't') : base === 12n ? RDF : unaryPolynomialBase(base);
      const R = new PolynomialRing(B as CoefficientRing<RingElement>, 'x'),
        f = R.__call__(a);
      const x =
        kind === -1n
          ? B.__call__(n)
          : kind === 27n
            ? RDF.__call__(n)
            : kind === 28n
              ? RDF.__call__(Number.NaN)
              : kind === 29n
                ? RDF.__call__(n < 0n ? Number.NEGATIVE_INFINITY : Number.POSITIVE_INFINITY)
                : kind === 30n
                  ? f
                  : kind === 31n
                    ? R.gen()
                    : kind === 32n
                      ? R.zero()
                      : kind === 33n
                        ? new PolynomialRing(QQ, 'y').__call__([n, 1n])
                        : kind === 34n
                          ? new PolynomialRing(B as CoefficientRing<RingElement>, 'y').__call__([
                              n,
                              1n,
                            ])
                          : kind === 35n
                            ? R.__call__([n, 1n])
                            : kind === 36n
                              ? R.__call__([n, 0n, 1n])
                              : polynomialInput(kind, n);
      const h = (f.evaluate as (x: unknown) => RingElement)(x);
      const P =
        h instanceof Integer || typeof h === 'bigint'
          ? 'Integer Ring'
          : h instanceof Rational
            ? 'Rational Field'
            : typeof h === 'number'
              ? "<class 'float'>"
              : String((h as unknown as { parent: unknown }).parent);
      let value = String(h);
      if (typeof h === 'number') {
        const b = Buffer.alloc(8);
        b.writeDoubleBE(h);
        value = 'float:' + (Number.isNaN(h) ? 'nan' : b.toString('hex'));
      }
      return [value, P, h === f, h instanceof Polynomial && h === x];
    }),
});

Object.assign(functions, {
  polynomial_generator_identity: (base: bigint, kind: bigint) =>
    comparison(() => {
      const B =
        base === 11n ? new PolynomialRing(QQ, 't') : base === 12n ? RDF : unaryPolynomialBase(base);
      const R = new PolynomialRing(B as CoefficientRing<RingElement>, 'x'),
        x = R.gen();
      const values = [
        x,
        R.__call__([0n, 1n]),
        R.__call__(x),
        x.add(R.zero()),
        x.mul(R.one()),
        x.neg().neg(),
        x.pow(1n),
        x.pow(0n),
        R.zero(),
        R.one(),
        x.compose(x),
        R.__call__([0n, 1n]).compose(x),
        x.compose(R.__call__([0n, 1n])),
      ];
      const f = values[Number(kind)]!;
      return [f.is_gen(), f === x, x === R.gen(), String(f)];
    }),
  polynomial_composition_input: (
    base: bigint,
    a: bigint[],
    b: bigint[],
    target: bigint,
    method: bigint
  ) =>
    comparison(() => {
      const B =
        base === 11n ? new PolynomialRing(QQ, 't') : base === 12n ? RDF : unaryPolynomialBase(base);
      const R = new PolynomialRing(B as CoefficientRing<RingElement>, 'x');
      const S =
        target === 0n
          ? R
          : new PolynomialRing(
              (target === 1n ? B : target === 2n ? QQ : RDF) as CoefficientRing<RingElement>,
              target === 1n ? 'y' : 'x'
            );
      const f = R.__call__(a),
        g = S.__call__(b);
      const h = method === 0n ? f.compose(g) : f.evaluate(g);
      const pack = (p: Polynomial<RingElement>) =>
        p.coeffs.map((c) => {
          if (c instanceof RealDoubleElement) {
            const bytes = Buffer.alloc(8);
            bytes.writeDoubleBE(c.value);
            return Number.isNaN(c.value) ? 'nan' : bytes.toString('hex');
          }
          return String(c);
        });
      return [
        String(h),
        String((h as Polynomial<RingElement>).parent),
        h === f,
        h === g,
        pack(h as Polynomial<RingElement>),
      ];
    }),
});

Object.assign(functions, {
  polynomial_evaluation_stress: (
    base: bigint,
    length: bigint,
    pattern: bigint,
    den: bigint,
    num: bigint,
    xden: bigint,
    kind: bigint
  ) =>
    comparison(() => {
      const B = base === 12n ? RDF : unaryPolynomialBase(base),
        N = Number(length);
      const coeffs = Array.from({ length: N }, (_, i) => {
        const n =
          (pattern === 1n && i !== 0 && i !== N - 1) || (pattern === 2n && i % 13 !== 0)
            ? 0n
            : BigInt((i % 7) - 3) * (pattern === 3n ? 1n << 300n : 1n);
        return base === 6n ? B.__call__(n) : B.__call__(new Rational(n, den));
      });
      const f = new PolynomialRing(B as CoefficientRing<RingElement>, 'x').__call__(coeffs),
        q = new Rational(num, xden);
      const x =
        kind === 0n
          ? q
          : kind === 1n
            ? RDF.__call__(q)
            : kind === 2n
              ? RDF.__call__(num < 0n ? -Infinity : Infinity)
              : kind === 3n
                ? q.toNumber()
                : Number.NaN;
      const pack = (h: RingElement | number) => {
        if (typeof h === 'number' || h instanceof RealDoubleElement) {
          const v = typeof h === 'number' ? h : h.value,
            b = Buffer.alloc(8);
          b.writeDoubleBE(v);
          return [
            Number.isNaN(v) ? 'nan' : b.toString('hex'),
            typeof h === 'number' ? "<class 'float'>" : 'Real Double Field',
          ];
        }
        return [
          String(h),
          h instanceof Rational
            ? 'Rational Field'
            : h instanceof Integer
              ? 'Integer Ring'
              : String((h as { parent: unknown }).parent),
        ];
      };
      return [pack(f.evaluate(x)), pack(f.evaluate(x))];
    }),
});

Object.assign(functions, {
  polynomial_generator_constructor: (base: bigint, a: bigint[], flag: bigint) =>
    comparison(() => {
      const B =
        base === 11n ? new PolynomialRing(QQ, 't') : base === 12n ? RDF : unaryPolynomialBase(base);
      const R = new PolynomialRing(B as CoefficientRing<RingElement>, 'x');
      const f = new Polynomial(a.map((c) => B.__call__(c)) as RingElement[], R, flag !== 0n);
      return [String(f), f.degree(), f.is_gen(), f.coeffs.map(String), String(f.evaluate(2n))];
    }),
});

import { Matrix as GenericMatrix } from '../../../../packages/sagemath-ts/src/matrix/matrix_generic.js';
import { Matrix_modn_dense } from '../../../../packages/sagemath-ts/src/matrix/matrix_modn.js';
Object.assign(functions, {
  polynomial_matrix_evaluation: (
    base: bigint,
    target: bigint,
    rows: bigint,
    cols: bigint,
    pattern: bigint,
    a: bigint[]
  ) =>
    comparison(() => {
      const B =
        base === 11n ? new PolynomialRing(QQ, 't') : base === 12n ? RDF : unaryPolynomialBase(base);
      const T =
        target === 0n
          ? polynomialIntegerBase
          : target === 1n
            ? QQ
            : target === 2n
              ? RDF
              : target === 3n
                ? Zmod(7n)
                : target === 4n
                  ? Zmod(14n)
                  : target === 5n
                    ? unaryPolynomialBase(1n)
                    : target === 6n
                      ? unaryPolynomialBase(4n)
                      : target === 7n || target === 9n
                        ? unaryPolynomialBase(5n)
                        : Zmod(1n);
      const m = Number(rows),
        n = Number(cols),
        entries = Array.from({ length: m }, (_, i) =>
          Array.from({ length: n }, (_, j) =>
            pattern === 0n
              ? 0n
              : pattern === 1n
                ? i === j
                  ? 1n
                  : 0n
                : BigInt(((i * n + j) % 7) - 2)
          )
        );
      const point =
        target === 9n
          ? new Matrix_mod2_dense(
              m,
              n,
              entries.map((row) => row.map((v) => Number(v & 1n)))
            )
          : target === 0n
            ? new IntegerMatrix(m, n, entries)
            : target === 3n || target === 4n
              ? new Matrix_modn_dense(m, n, target === 3n ? 7n : 14n, entries)
              : new GenericMatrix(
                  T as CoefficientRing<RingElement>,
                  m,
                  n,
                  entries.map((row) => row.map((c) => T.__call__(c)) as RingElement[])
                );
      const f = new PolynomialRing(B as CoefficientRing<RingElement>, 'x').__call__(a),
        h = f.evaluate(point) as unknown as
          | IntegerMatrix
          | Matrix_modn_dense
          | GenericMatrix<RingElement>
          | Matrix_mod2_dense;
      const resultBase =
        h instanceof Matrix_mod2_dense
          ? 'Finite Field of size 2'
          : h instanceof IntegerMatrix
            ? 'Integer Ring'
            : h instanceof Matrix_modn_dense
              ? String(Zmod(h.modulus))
              : String(h.base_ring);
      return [
        Array.from({ length: h.nrows }, (_, i) =>
          Array.from({ length: h.ncols }, (_, j) => String(h.get(i, j)))
        ),
        resultBase,
        h.nrows,
        h.ncols,
        h === point,
      ];
    }),
});
functions.polynomial_matrix_actions = (
  base: bigint,
  target: bigint,
  n: bigint,
  wrap: bigint,
  degree: bigint,
  den: bigint,
  a: bigint[],
  b: bigint[]
) =>
  comparison(() => {
    const B = (
      base === 11n ? new PolynomialRing(QQ, 't') : base === 12n ? RDF : unaryPolynomialBase(base)
    ) as CoefficientRing<RingElement>;
    const T = (
        target === 0n ? polynomialIntegerBase : target === 1n ? QQ : RDF
      ) as CoefficientRing<RingElement>,
      m = Number(n);
    const entries = Array.from({ length: m }, (_, i) =>
      Array.from({ length: m }, (_, j) => {
        const v = b[(i * m + j) % b.length]!;
        return target === 2n
          ? RDF.__call__(doubleFromBits(v))
          : target === 1n
            ? QQ.__call__(v).div(QQ.__call__(den))
            : T.__call__(v);
      })
    );
    const point =
      target === 0n
        ? new IntegerMatrix(
            m,
            m,
            entries.map((row) => row.map((c) => (c as Integer).value))
          )
        : new GenericMatrix(T, m, m, entries);
    const cs: RingElement[] = Array.from({ length: Number(degree) + 1 }, () => B.zero());
    a.forEach((v, i) => {
      cs[Math.floor((i * Number(degree)) / Math.max(1, a.length - 1))] =
        base === 12n
          ? RDF.__call__(doubleFromBits(v))
          : base === 11n
            ? B.__call__([
                QQ.__call__(v).div(QQ.__call__(den)),
                QQ.__call__(BigInt(i + 1)).div(QQ.__call__(den)),
              ])
            : B.__call__(QQ.__call__(v).div(QQ.__call__(den)));
    });
    const f = new PolynomialRing(B, 'x').__call__(cs);
    const frame = (h: EvaluationActionMatrix) => {
      const hb =
        h instanceof IntegerMatrix
          ? polynomialIntegerBase
          : h instanceof Matrix_modn_dense
            ? Zmod(h.modulus)
            : h.base_ring;
      return [
        Array.from({ length: h.nrows }, (_, i) =>
          Array.from({ length: h.ncols }, (_, j) =>
            hb === RDF ? doubleBits((h.get(i, j) as RealDoubleElement).value) : String(h.get(i, j))
          )
        ),
        String(hb),
        h.nrows,
        h.ncols,
        h === point,
      ];
    };
    const first = frame(f.evaluate(wrap ? [point] : point) as EvaluationActionMatrix);
    f.evaluate(B.__call__(2n));
    const second = frame(f.evaluate(wrap ? [point] : point) as EvaluationActionMatrix);
    return [first, second];
  });
type EvaluationActionMatrix = IntegerMatrix | Matrix_modn_dense | GenericMatrix<RingElement>;

import { integer_to_real_double_dense } from '../../../../packages/sagemath-ts/src/matrix/change_ring.js';
functions.polynomial_matrix_double_conversion = (rows: bigint, cols: bigint, a: bigint[]) =>
  comparison(() => {
    const m = Number(rows),
      n = Number(cols),
      h = integer_to_real_double_dense(
        new IntegerMatrix(
          m,
          n,
          Array.from({ length: m }, (_, i) =>
            Array.from({ length: n }, (_, j) => a[(i * n + j) % a.length]!)
          )
        )
      );
    return [
      Array.from({ length: m }, (_, i) =>
        Array.from({ length: n }, (_, j) => doubleBits(h.get(i, j).value))
      ),
      String(h.base_ring),
      h.nrows,
      h.ncols,
    ];
  });

import { Matrix_mod2_dense } from '../../../../packages/sagemath-ts/src/matrix/matrix_mod2.js';

functions.rdf_scalar = rdf_scalar;

functions.rdf_alias = rdf_alias;

functions.flint_remainder = flint_remainder;

functions.polynomial_roots_fidelity = (kind: bigint, coefficients: bigint[], denominator: bigint) =>
  comparison(() => {
    const R = new PolynomialRing(
      (kind === 0n ? polynomialIntegerBase : QQ) as unknown as CoefficientRing<RingElement>,
      'x'
    );
    const f = R.__call__(coefficients.map((c) => (kind === 0n ? c : new Rational(c, denominator))));
    return f
      .roots()
      .map(([r, m]) => [
        String(r),
        String(m),
        r instanceof Integer ? 'Integer Ring' : r instanceof Rational ? 'Rational Field' : typeof r,
      ]);
  });

functions.polynomial_factor_fidelity = (
  op: bigint,
  kind: bigint,
  coefficients: bigint[],
  denominator: bigint
) =>
  comparison(() => {
    const B =
      kind === 0n
        ? polynomialIntegerBase
        : kind === 1n
          ? QQ
          : kind === 2n
            ? new PrimeField(7n)
            : kind === 3n
              ? GF2
              : new PrimeField((1n << 89n) - 1n);
    const R = new PolynomialRing(B as unknown as CoefficientRing<RingElement>, 'x');
    const f = R.__call__(
      coefficients.map((c) => (denominator === 1n ? c : new Rational(c, denominator)))
    );
    if (op === 1n) return f.is_irreducible();
    return f.factor().map(([g, e]) => [g.coeffs.map(String), String(e)]);
  });
import {
  getrand as polynomialGetrand,
  setrand as polynomialSetrand,
} from '../../../../packages/parigp-ts/src/random.js';
functions.polynomial_factor_state = (
  op: bigint,
  kind: bigint,
  coefficients: bigint[],
  denominator: bigint
) => {
  const saved = polynomialGetrand();
  try {
    polynomialSetrand(1n);
    const R = new PolynomialRing(
      (kind === 0n ? polynomialIntegerBase : QQ) as unknown as CoefficientRing<RingElement>,
      'x'
    );
    const f = R.__call__(
      coefficients.map((c) => (denominator === 1n ? c : new Rational(c, denominator)))
    );
    let value: unknown;
    if (op === 0n)
      value = [f.factor().map(([g, e]) => [g.coeffs.map(String), String(e)]), polynomialGetrand()];
    else {
      value = [
        [f.is_irreducible(), polynomialGetrand()],
        [f.is_irreducible(), polynomialGetrand()],
        [new Polynomial([...f.coeffs], f.parent).is_irreducible(), polynomialGetrand()],
      ];
    }
    return JSON.stringify(value, (_, v) => (typeof v === 'bigint' ? String(v) : v));
  } finally {
    polynomialSetrand(saved);
  }
};

functions.finite_polynomial_roots = (
  kind: bigint,
  p: bigint,
  degree: bigint,
  modulus: bigint[],
  coefficients: bigint[]
) =>
  comparison(() => {
    const B = kind === 0n ? new PrimeField(p)
      : kind === 1n ? new FiniteFieldPrime(p)
        : kind === 2n ? GF2 : GFpn(p, Number(degree), modulus.slice(0, -1).map(Number), 'a');
    const R = new PolynomialRing(B as unknown as CoefficientRing<RingElement>, 'x');
    const f = R.__call__(coefficients.map(c => kind === 3n
      ? (B as ReturnType<typeof GFpn>).fromInteger(c) : B.__call__(c)));
    return f.roots().map(([r, m]) => [
      kind === 3n ? (r as ReturnType<ReturnType<typeof GFpn>['fromInteger']>).integer_representation() : String(r),
      String(m), (r as unknown as { parent: unknown }).parent === B,
    ]);
  });

import { cmp_universal as pariFiniteCompare } from '../../../../packages/parigp-ts/src/gen2.js';
import { PariType } from '../../../../packages/parigp-ts/src/types.js';
functions.pari_ffelt_compare = (p: bigint, T: bigint[], x: bigint[], q: bigint, U: bigint[], y: bigint[]) =>
  pariFiniteCompare({ type: PariType.t_FFELT, p, degree: T.length - 1, definingPoly: T, value: x },
    { type: PariType.t_FFELT, p: q, degree: U.length - 1, definingPoly: U, value: y });
functions.finite_polynomial_factors = (
  kind: bigint,
  p: bigint,
  degree: bigint,
  modulus: bigint[],
  coefficients: bigint[]
) =>
  comparison(() => {
    const B = kind === 0n ? new PrimeField(p)
      : kind === 1n ? new FiniteFieldPrime(p)
        : kind === 2n ? GF2 : GFpn(p, Number(degree), modulus.slice(0, -1).map(Number), 'a');
    const R = new PolynomialRing(B as unknown as CoefficientRing<RingElement>, 'x');
    const f = R.__call__(coefficients.map(c => kind === 3n
      ? (B as ReturnType<typeof GFpn>).fromInteger(c) : B.__call__(c)));
    return f.factor().map(([g, e]) => [g.coeffs.map(c => kind === 3n
      ? String((c as ReturnType<ReturnType<typeof GFpn>['fromInteger']>).integer_representation()) : String(c)), String(e)]);
  });

functions.pari_ffelt_compare_scalar = (p: bigint, _T: bigint[], x: bigint[], q: bigint, _U: bigint[], y: bigint[]) =>
  pariFiniteCompare({ type: PariType.t_FFELT, p, degree: 1, value: x[0]! },
    { type: PariType.t_FFELT, p: q, degree: 1, value: y[0]! });

functions.finite_polynomial_roots_state = (kind: bigint, p: bigint, coefficients: bigint[]) => {
  const saved = polynomialGetrand();
  try {
    polynomialSetrand(1n);
    const B = kind === 0n ? new PrimeField(p) : kind === 1n ? new FiniteFieldPrime(p) : GF2;
    const R = new PolynomialRing(B as unknown as CoefficientRing<RingElement>, 'x');
    return comparison(() => [R.__call__(coefficients).roots().map(([r, m]) => [String(r), String(m),
      (r as unknown as { parent: unknown }).parent === B]), polynomialGetrand()]);
  } finally { polynomialSetrand(saved); }
};
functions.finite_polynomial_distinct_roots = (
  kind: bigint,
  p: bigint,
  degree: bigint,
  modulus: bigint[],
  coefficients: bigint[]
) =>
  comparison(() => {
    const B = kind === 0n ? new PrimeField(p)
      : kind === 1n ? new FiniteFieldPrime(p)
        : kind === 2n ? GF2 : GFpn(p, Number(degree), modulus.slice(0, -1).map(Number), 'a');
    const R = new PolynomialRing(B as unknown as CoefficientRing<RingElement>, 'x');
    const f = R.__call__(coefficients.map(c => kind === 3n
      ? (B as ReturnType<typeof GFpn>).fromInteger(c) : B.__call__(c)));
    return f.roots({ multiplicities: false }).map((r) => [
      kind === 3n ? (r as ReturnType<ReturnType<typeof GFpn>['fromInteger']>).integer_representation() : String(r),
      (r as unknown as { parent: unknown }).parent === B,
    ]);
  });

functions.polynomial_distinct_roots_fidelity = (kind: bigint, coefficients: bigint[], denominator: bigint) =>
  comparison(() => {
    const R = new PolynomialRing((kind === 0n ? polynomialIntegerBase : QQ) as unknown as CoefficientRing<RingElement>, 'x');
    const f = R.__call__(coefficients.map(c => kind === 0n ? c : new Rational(c, denominator)));
    return f.roots({ multiplicities: false }).map(r => [String(r),
      r instanceof Integer ? 'Integer Ring' : r instanceof Rational ? 'Rational Field' : typeof r]);
  });

functions.finite_polynomial_distinct_roots_state = (kind: bigint, p: bigint, coefficients: bigint[]) => {
  const saved = polynomialGetrand();
  try {
    polynomialSetrand(1n);
    const B = kind === 0n ? new PrimeField(p) : kind === 1n ? new FiniteFieldPrime(p) : GF2;
    const R = new PolynomialRing(B as unknown as CoefficientRing<RingElement>, 'x');
    return comparison(() => [R.__call__(coefficients).roots({ multiplicities: false }).map(r => [String(r),
      (r as unknown as { parent: unknown }).parent === B]), polynomialGetrand()]);
  } finally { polynomialSetrand(saved); }
};
