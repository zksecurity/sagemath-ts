/**
 * sagemath-ts side of the `padics_series` property-test area.
 *
 * Covers `src/rings/padics/*`, `src/rings/power_series_ring.ts` and
 * `src/rings/laurent_series_ring.ts`.
 *
 * Cases: tests/property/cases/padics_series.cases.json
 * SageMath counterpart: tests/property/python/areas/padics_series.py
 *
 * Every function mirrors the SageMath module argument-for-argument and
 * string-for-string; see that file for the encoding and output conventions.
 *
 * The `QQ` adapter below is the only thing in this file that is not the port:
 * `PowerSeriesRing` takes its coefficient ring as an interface, and SageMath's
 * `QQ` is what the ported algorithms expect (`sqrt`, `is_square` and `nth_root`
 * on the coefficients are all used by `power_series_ring.ts`, exactly as
 * `power_series_ring_element.pyx` uses them). It models `QQ` and nothing more --
 * no arithmetic under test is implemented here.
 */

import { GF } from '../../../../packages/sagemath-ts/src/rings/finite_rings/finite_field_constructor.js';
import {
  type LaurentSeriesElement,
  LaurentSeriesRing,
} from '../../../../packages/sagemath-ts/src/rings/laurent_series_ring.js';
import {
  Qp,
  Zp,
  type pAdicGeneric,
  type pAdicGenericElement,
} from '../../../../packages/sagemath-ts/src/rings/padics/index.js';
import {
  type CoefficientRing,
  type PowerSeriesElement,
  PowerSeriesRing,
  type RingElement,
} from '../../../../packages/sagemath-ts/src/rings/power_series_ring.js';
import { Rational } from '../../../../packages/sagemath-ts/src/rings/rational.js';

/** Sentinel meaning "argument omitted" (SageMath `None`). */
const NONE = -999999n;

// ---------------------------------------------------------------------------
// formatting helpers (must produce byte-identical strings to the Python side)
// ---------------------------------------------------------------------------

/** Run `fn`, turning any exception into the same string Python's `_guard` builds. */
function guard(fn: () => string): string {
  try {
    return fn();
  } catch (e) {
    const name =
      e instanceof Error ? (e.constructor?.name ?? e.name) : (typeof e as unknown as string);
    const message = e instanceof Error ? e.message : String(e);
    return `ERROR:${name}: ${message}`;
  }
}

/** SageMath prints `Infinity` as `+Infinity`. */
function num(v: number | bigint): string {
  if (typeof v === 'bigint') return v.toString();
  if (v === Number.POSITIVE_INFINITY) return '+Infinity';
  if (v === Number.NEGATIVE_INFINITY) return '-Infinity';
  return String(v);
}

function bool(b: boolean): string {
  return b ? 'True' : 'False';
}

/** Decode an optional integer argument. */
function opt(v: bigint): number | null {
  return v === NONE ? null : Number(v);
}

function fmtPadic(x: pAdicGenericElement): string {
  return `${x.toString()} [v=${num(x.valuation())}, rp=${num(x.precision_relative())}, ap=${num(
    x.precision_absolute()
  )}]`;
}

function fmtSeries(f: {
  toString(): string;
  valuation(): number;
  prec(): number;
  degree(): number;
}): string {
  return `${f.toString()} [v=${num(f.valuation())}, prec=${num(f.prec())}, deg=${num(f.degree())}]`;
}

function fmtList(items: readonly { toString(): string }[]): string {
  return `[${items.map((x) => x.toString()).join(', ')}]`;
}

// ---------------------------------------------------------------------------
// p-adic helpers
// ---------------------------------------------------------------------------

function padicRing(p: bigint, prec: bigint, isField: bigint): pAdicGeneric {
  return isField !== 0n ? Qp(p, Number(prec)) : Zp(p, Number(prec));
}

/** Decode `[p, prec, is_field, num, den]` into `[R, R(num)/R(den)]`. */
function padicElt(spec: bigint[]): [pAdicGeneric, pAdicGenericElement] {
  const R = padicRing(spec[0]!, spec[1]!, spec[2]!);
  const a = R.__call__(spec[3]!);
  if (spec[4] === 1n) {
    return [R, a];
  }
  return [R, a.div(R.__call__(spec[4]!))];
}

// ---------------------------------------------------------------------------
// p-adic rings
// ---------------------------------------------------------------------------

function padic_ring_info(spec: bigint[]): string {
  return guard(() => {
    const R = padicRing(spec[0]!, spec[1]!, spec[2]!);
    return [
      R.toString(),
      `p=${R.prime()}`,
      `cap=${R.precision_cap()}`,
      `field=${bool(R.is_field())}`,
      `char=${R.characteristic()}`,
      `res_char=${R.residue_characteristic()}`,
      `e=${R.absolute_e()}`,
      `f=${R.absolute_f()}`,
      `unif=${R.uniformizer()}`,
      `zero=${fmtPadic(R.zero())}`,
      `one=${fmtPadic(R.one())}`,
    ].join(' | ');
  });
}

function padic_uniformizer_pow(spec: bigint[]): string {
  return guard(() => {
    const R = padicRing(spec[0]!, spec[1]!, spec[2]!);
    return fmtPadic(R.uniformizer_pow(spec[3]!));
  });
}

function padic_teichmuller(spec: bigint[]): string {
  return guard(() => {
    const R = padicRing(spec[0]!, spec[1]!, spec[2]!);
    return fmtPadic(R.teichmuller(spec[3]!));
  });
}

function padic_teichmuller_system(spec: bigint[]): string {
  return guard(() => {
    const R = padicRing(spec[0]!, spec[1]!, spec[2]!);
    return fmtList(R.teichmuller_system());
  });
}

function padic_roots_of_unity(spec: bigint[]): string {
  return guard(() => {
    const R = padicRing(spec[0]!, spec[1]!, spec[2]!);
    const n = spec[3]! === NONE ? undefined : spec[3]!;
    return fmtList(R.roots_of_unity(n));
  });
}

// ---------------------------------------------------------------------------
// p-adic elements
// ---------------------------------------------------------------------------

function padic_repr(spec: bigint[]): string {
  return guard(() => fmtPadic(padicElt(spec)[1]));
}

function padic_parent(spec: bigint[]): string {
  return guard(() => padicElt(spec)[1].parent().toString());
}

function padic_from_int(p: bigint, prec: bigint, value: bigint): string {
  return guard(() => fmtPadic(Zp(p, Number(prec)).__call__(value)));
}

function padic_from_int_field(p: bigint, prec: bigint, value: bigint): string {
  return guard(() => fmtPadic(Qp(p, Number(prec)).__call__(value)));
}

function padic_expansion(spec: bigint[]): string {
  return guard(() => fmtList(padicElt(spec)[1].expansion()));
}

function padic_flags(spec: bigint[]): string {
  return guard(() => {
    const x = padicElt(spec)[1];
    return [
      `zero=${bool(x.is_zero())}`,
      `one=${bool(x.is_one())}`,
      `unit=${bool(x.is_unit())}`,
      `integral=${bool(x.is_integral())}`,
      `val=${num(x.valuation())}`,
      `lift=${guard(() => x.lift().toString())}`,
    ].join(' | ');
  });
}

function padic_unit_part(spec: bigint[]): string {
  return guard(() => fmtPadic(padicElt(spec)[1].unit_part()));
}

function binaryPadic(
  spec: bigint[],
  op: (a: pAdicGenericElement, b: pAdicGenericElement) => pAdicGenericElement
): string {
  return guard(() => {
    const a = padicElt(spec)[1];
    const b = padicElt([spec[0]!, spec[1]!, spec[2]!, spec[5]!, spec[6]!])[1];
    return fmtPadic(op(a, b));
  });
}

function padic_add(spec: bigint[]): string {
  return binaryPadic(spec, (a, b) => a.add(b));
}

function padic_sub(spec: bigint[]): string {
  return binaryPadic(spec, (a, b) => a.sub(b));
}

function padic_mul(spec: bigint[]): string {
  return binaryPadic(spec, (a, b) => a.mul(b));
}

function padic_div(spec: bigint[]): string {
  return binaryPadic(spec, (a, b) => a.div(b));
}

function padic_pow(spec: bigint[]): string {
  return guard(() => fmtPadic(padicElt(spec)[1].pow(spec[5]!)));
}

function padic_inv(spec: bigint[]): string {
  return guard(() => fmtPadic(padicElt(spec)[1].pow(-1n)));
}

function padic_residue(spec: bigint[]): string {
  return guard(() => padicElt(spec)[1].residue(Number(spec[5]!)).toString());
}

function padic_getitem(spec: bigint[]): string {
  return guard(() => padicElt(spec)[1].__getitem__(Number(spec[5]!)).toString());
}

function padic_slice(spec: bigint[]): string {
  return guard(() =>
    fmtPadic(padicElt(spec)[1].slice(opt(spec[5]!), opt(spec[6]!), opt(spec[7]!)))
  );
}

function padic_add_bigoh(spec: bigint[]): string {
  return guard(() => fmtPadic(padicElt(spec)[1].add_bigoh(Number(spec[5]!))));
}

function padic_lift_to_precision(spec: bigint[]): string {
  return guard(() =>
    fmtPadic(padicElt(spec)[1].add_bigoh(Number(spec[5]!)).lift_to_precision(Number(spec[6]!)))
  );
}

function padic_inexact_zero(spec: bigint[]): string {
  return guard(() => {
    const R = padicRing(spec[0]!, spec[1]!, spec[2]!);
    const z = R.__call__(0n, Number(spec[3]!));
    return [
      fmtPadic(z),
      `is_zero=${bool(z.is_zero())}`,
      `exp=${fmtList(z.expansion())}`,
      `add_order=${guard(() => num(z.additive_order()))}`,
      `plus_one=${fmtPadic(z.add(R.__call__(1n)))}`,
    ].join(' | ');
  });
}

function padic_sqrt(spec: bigint[]): string {
  return guard(() => fmtPadic(padicElt(spec)[1].sqrt()));
}

function padic_sqrt_all(spec: bigint[]): string {
  return guard(() => fmtList(padicElt(spec)[1].square_root_all()));
}

function padic_is_square(spec: bigint[]): string {
  return guard(() => bool(padicElt(spec)[1].is_square()));
}

function padic_nth_root(spec: bigint[]): string {
  return guard(() => fmtPadic(padicElt(spec)[1].nth_root(spec[5]!)));
}

function padic_nth_root_all(spec: bigint[]): string {
  return guard(() => fmtList(padicElt(spec)[1].nth_root_all(spec[5]!)));
}

function padic_log(spec: bigint[]): string {
  return guard(() => fmtPadic(padicElt(spec)[1].log()));
}

function padic_exp(spec: bigint[]): string {
  return guard(() => fmtPadic(padicElt(spec)[1].exp()));
}

function padic_artin_hasse_exp(spec: bigint[]): string {
  return guard(() => fmtPadic(padicElt(spec)[1].artin_hasse_exp()));
}

function padic_multiplicative_order(spec: bigint[]): string {
  return guard(() => num(padicElt(spec)[1].multiplicative_order()));
}

function padic_additive_order(spec: bigint[]): string {
  return guard(() => num(padicElt(spec)[1].additive_order()));
}

// ---------------------------------------------------------------------------
// coefficient rings for the series tests
// ---------------------------------------------------------------------------

/** Exact integer `n`-th root of a non-negative bigint, or `null`. */
function exactRoot(a: bigint, n: number): bigint | null {
  if (a < 0n) return null;
  if (a < 2n) return a;
  let lo = 1n;
  let hi = 1n;
  while (hi ** BigInt(n) <= a) hi *= 2n;
  while (lo < hi) {
    const mid = (lo + hi + 1n) / 2n;
    if (mid ** BigInt(n) <= a) lo = mid;
    else hi = mid - 1n;
  }
  return lo ** BigInt(n) === a ? lo : null;
}

/** `Rational` wrapped in the `RingElement` interface (SageMath's `QQ`). */
class RationalElement implements RingElement {
  readonly value: Rational;

  constructor(value: Rational) {
    this.value = value;
  }

  /**
   * SageMath's `Rational.is_square()`.
   * `QQ(4).is_square()` is `True`, `QQ(2).is_square()` is `False`.
   */
  is_square(): boolean {
    if (this.value.numerator < 0n) return false;
    return (
      exactRoot(this.value.numerator, 2) !== null && exactRoot(this.value.denominator, 2) !== null
    );
  }

  /**
   * SageMath's `Rational.sqrt(extend=False)`: an exact rational root, or
   * `ValueError: square root of x not a rational number`.
   */
  sqrt(): RingElement {
    const n = exactRoot(this.value.numerator, 2);
    const d = exactRoot(this.value.denominator, 2);
    if (n === null || d === null) {
      throw new Error(`square root of ${this.value} not a rational number`);
    }
    return new RationalElement(new Rational(n, d));
  }

  /**
   * SageMath's `Rational.nth_root(k)`: `QQ(8).nth_root(3)` is `2`,
   * `QQ(2).nth_root(3)` raises `ValueError: not a perfect 3rd power`.
   */
  nth_root(k: number): RingElement {
    const neg = this.value.numerator < 0n;
    if (neg && k % 2 === 0) {
      throw new Error('cannot take even root of negative number');
    }
    const n = exactRoot(neg ? -this.value.numerator : this.value.numerator, k);
    const d = exactRoot(this.value.denominator, k);
    if (n === null || d === null) {
      throw new Error(`not a perfect ${k}th power`);
    }
    return new RationalElement(new Rational(neg ? -n : n, d));
  }

  add(other: RingElement): RingElement {
    return new RationalElement(this.value.add((other as RationalElement).value));
  }
  sub(other: RingElement): RingElement {
    return new RationalElement(this.value.sub((other as RationalElement).value));
  }
  mul(other: RingElement): RingElement {
    return new RationalElement(this.value.mul((other as RationalElement).value));
  }
  div(other: RingElement): RingElement {
    return new RationalElement(this.value.div((other as RationalElement).value));
  }
  neg(): RingElement {
    return new RationalElement(this.value.neg());
  }
  inv(): RingElement {
    return new RationalElement(this.value.inv());
  }
  eq(other: RingElement | number | bigint): boolean {
    if (typeof other === 'number') return this.value.eq(BigInt(other));
    if (typeof other === 'bigint') return this.value.eq(other);
    return this.value.eq((other as RationalElement).value);
  }
  isZero(): boolean {
    return this.value.eq(0n);
  }
  isOne(): boolean {
    return this.value.eq(1n);
  }
  isUnit(): boolean {
    return !this.isZero();
  }
  toString(): string {
    return this.value.toString();
  }
}

const QQ: CoefficientRing<RingElement> = {
  zero: () => new RationalElement(Rational.zero()),
  one: () => new RationalElement(Rational.one()),
  __call__: (x: unknown): RingElement => {
    if (x instanceof RationalElement) return x;
    if (x instanceof Rational) return new RationalElement(x);
    if (typeof x === 'bigint') return new RationalElement(new Rational(x, 1n));
    if (typeof x === 'number') return new RationalElement(new Rational(BigInt(x), 1n));
    if (typeof x === 'string') return new RationalElement(Rational.fromString(x));
    throw new TypeError(`cannot convert ${typeof x} to a rational`);
  },
  is_field: () => true,
  characteristic: () => 0n,
  element_is_atomic: () => true,
  toString: () => 'Rational Field',
};

/** `0` -> QQ, otherwise the prime field `GF(code)`. */
function baseRing(code: bigint): CoefficientRing<RingElement> {
  if (code === 0n) return QQ;
  return GF(code) as unknown as CoefficientRing<RingElement>;
}

// ---------------------------------------------------------------------------
// power series
// ---------------------------------------------------------------------------

function psRing(params: bigint[]): PowerSeriesRing<RingElement> {
  return new PowerSeriesRing<RingElement>(baseRing(params[0]!), 'x', Number(params[1]!));
}

function ps(
  params: bigint[],
  coeffs: bigint[],
  precIndex: number
): [PowerSeriesRing<RingElement>, PowerSeriesElement<RingElement>] {
  const R = psRing(params);
  const prec = opt(params[precIndex]!);
  return [R, prec === null ? R.__call__(coeffs) : R.__call__(coeffs, prec)];
}

function ps_repr(coeffs: bigint[], params: bigint[]): string {
  return guard(() => fmtSeries(ps(params, coeffs, 2)[1]));
}

function ps_info(coeffs: bigint[], params: bigint[]): string {
  return guard(() => {
    const f = ps(params, coeffs, 2)[1];
    return [
      fmtSeries(f),
      `ap=${num(f.precision_absolute())}`,
      `rp=${num(f.precision_relative())}`,
      `zero=${bool(f.is_zero())}`,
      `one=${bool(f.is_one())}`,
      `unit=${bool(f.is_unit())}`,
      `monomial=${bool(f.is_monomial())}`,
      `list=${fmtList(f.list())}`,
    ].join(' | ');
  });
}

function binaryPs(
  c1: bigint[],
  c2: bigint[],
  params: bigint[],
  op: (
    f: PowerSeriesElement<RingElement>,
    g: PowerSeriesElement<RingElement>
  ) => PowerSeriesElement<RingElement>
): string {
  return guard(() => {
    const f = ps(params, c1, 2)[1];
    const g = ps(params, c2, 3)[1];
    return fmtSeries(op(f, g));
  });
}

function ps_add(c1: bigint[], c2: bigint[], params: bigint[]): string {
  return binaryPs(c1, c2, params, (f, g) => f.add(g));
}

function ps_sub(c1: bigint[], c2: bigint[], params: bigint[]): string {
  return binaryPs(c1, c2, params, (f, g) => f.sub(g));
}

function ps_mul(c1: bigint[], c2: bigint[], params: bigint[]): string {
  return binaryPs(c1, c2, params, (f, g) => f.mul(g));
}

function ps_div(c1: bigint[], c2: bigint[], params: bigint[]): string {
  return binaryPs(c1, c2, params, (f, g) => f.div(g));
}

function ps_compose(c1: bigint[], c2: bigint[], params: bigint[]): string {
  return binaryPs(c1, c2, params, (f, g) => f.__call__(g));
}

function ps_inv(coeffs: bigint[], params: bigint[]): string {
  return guard(() => fmtSeries(ps(params, coeffs, 2)[1].inv()));
}

function ps_pow(coeffs: bigint[], params: bigint[]): string {
  return guard(() => fmtSeries(ps(params, coeffs, 2)[1].pow(params[3]!)));
}

function ps_derivative(coeffs: bigint[], params: bigint[]): string {
  return guard(() => fmtSeries(ps(params, coeffs, 2)[1].derivative()));
}

function ps_integral(coeffs: bigint[], params: bigint[]): string {
  return guard(() => fmtSeries(ps(params, coeffs, 2)[1].integral()));
}

function ps_exp(coeffs: bigint[], params: bigint[]): string {
  return guard(() => {
    const f = ps(params, coeffs, 2)[1];
    const target = opt(params[3]!);
    return fmtSeries(target === null ? f.exp() : f.exp(target));
  });
}

function ps_log(coeffs: bigint[], params: bigint[]): string {
  return guard(() => {
    const f = ps(params, coeffs, 2)[1];
    const target = opt(params[3]!);
    return fmtSeries(target === null ? f.log() : f.log(target));
  });
}

function ps_sqrt(coeffs: bigint[], params: bigint[]): string {
  return guard(() => {
    const f = ps(params, coeffs, 2)[1];
    const target = opt(params[3]!);
    return fmtSeries(target === null ? f.sqrt() : f.sqrt(target));
  });
}

function ps_nth_root(coeffs: bigint[], params: bigint[]): string {
  return guard(() => {
    const f = ps(params, coeffs, 2)[1];
    const n = Number(params[3]!);
    const target = opt(params[4]!);
    return fmtSeries(target === null ? f.nth_root(n) : f.nth_root(n, target));
  });
}

function ps_pade(coeffs: bigint[], params: bigint[]): string {
  return guard(() =>
    ps(params, coeffs, 2)[1].pade(Number(params[3]!), Number(params[4]!)).toString()
  );
}

function ps_reverse(coeffs: bigint[], params: bigint[]): string {
  return guard(() => {
    const f = ps(params, coeffs, 2)[1];
    const precision = opt(params[3]!);
    return fmtSeries(precision === null ? f.reversion() : f.reversion(precision));
  });
}

function ps_V(coeffs: bigint[], params: bigint[]): string {
  return guard(() => fmtSeries(ps(params, coeffs, 2)[1].V(Number(params[3]!))));
}

function ps_shift(coeffs: bigint[], params: bigint[]): string {
  return guard(() => fmtSeries(ps(params, coeffs, 2)[1].shift(Number(params[3]!))));
}

function ps_truncate(coeffs: bigint[], params: bigint[]): string {
  return guard(() => {
    const f = ps(params, coeffs, 2)[1];
    const n = opt(params[3]!);
    return (n === null ? f.truncate() : f.truncate(n)).toString();
  });
}

function ps_truncate_powerseries(coeffs: bigint[], params: bigint[]): string {
  return guard(() => fmtSeries(ps(params, coeffs, 2)[1].truncate_powerseries(Number(params[3]!))));
}

function ps_add_bigoh(coeffs: bigint[], params: bigint[]): string {
  return guard(() => fmtSeries(ps(params, coeffs, 2)[1].add_bigoh(Number(params[3]!))));
}

function ps_O(coeffs: bigint[], params: bigint[]): string {
  return guard(() => fmtSeries(ps(params, coeffs, 2)[1].O(Number(params[3]!))));
}

function ps_valuation_zero_part(coeffs: bigint[], params: bigint[]): string {
  return guard(() => fmtSeries(ps(params, coeffs, 2)[1].valuation_zero_part()));
}

function ps_getitem(coeffs: bigint[], params: bigint[]): string {
  return guard(() => ps(params, coeffs, 2)[1].__getitem__(Number(params[3]!)).toString());
}

function ps_is_square(coeffs: bigint[], params: bigint[]): string {
  return guard(() => bool(ps(params, coeffs, 2)[1].is_square()));
}

// ---------------------------------------------------------------------------
// Laurent series
// ---------------------------------------------------------------------------

function lsRing(params: bigint[]): LaurentSeriesRing<RingElement> {
  return new LaurentSeriesRing<RingElement>(baseRing(params[0]!), 'x', Number(params[1]!));
}

function ls(
  params: bigint[],
  coeffs: bigint[],
  shiftIndex: number,
  precIndex: number
): [LaurentSeriesRing<RingElement>, LaurentSeriesElement<RingElement>] {
  const L = lsRing(params);
  const R = L.power_series_ring();
  const n = Number(params[shiftIndex]!);
  const prec = opt(params[precIndex]!);
  const f = R.__call__(coeffs);
  const g = L.__call__(f, n);
  return [L, prec === null ? g : g.add_bigoh(prec)];
}

function ls_repr(coeffs: bigint[], params: bigint[]): string {
  return guard(() => fmtSeries(ls(params, coeffs, 2, 3)[1]));
}

function ls_info(coeffs: bigint[], params: bigint[]): string {
  return guard(() => {
    const f = ls(params, coeffs, 2, 3)[1];
    return [
      fmtSeries(f),
      `ap=${num(f.precision_absolute())}`,
      `rp=${num(f.precision_relative())}`,
      `zero=${bool(f.is_zero())}`,
      `unit=${bool(f.is_unit())}`,
      `monomial=${bool(f.is_monomial())}`,
      `exponents=${fmtList(f.exponents())}`,
      `coefficients=${fmtList(f.coefficients())}`,
      `residue=${f.residue()}`,
      `vzp=${f.valuation_zero_part()}`,
    ].join(' | ');
  });
}

function binaryLs(
  c1: bigint[],
  c2: bigint[],
  params: bigint[],
  op: (
    f: LaurentSeriesElement<RingElement>,
    g: LaurentSeriesElement<RingElement>
  ) => LaurentSeriesElement<RingElement>
): string {
  return guard(() => {
    const f = ls(params, c1, 2, 3)[1];
    const g = ls(params, c2, 4, 5)[1];
    return fmtSeries(op(f, g));
  });
}

function ls_add(c1: bigint[], c2: bigint[], params: bigint[]): string {
  return binaryLs(c1, c2, params, (f, g) => f.add(g));
}

function ls_sub(c1: bigint[], c2: bigint[], params: bigint[]): string {
  return binaryLs(c1, c2, params, (f, g) => f.sub(g));
}

function ls_mul(c1: bigint[], c2: bigint[], params: bigint[]): string {
  return binaryLs(c1, c2, params, (f, g) => f.mul(g));
}

function ls_div(c1: bigint[], c2: bigint[], params: bigint[]): string {
  return binaryLs(c1, c2, params, (f, g) => f.div(g));
}

function ls_compose(c1: bigint[], c2: bigint[], params: bigint[]): string {
  return binaryLs(c1, c2, params, (f, g) => f.__call__(g));
}

function ls_inv(coeffs: bigint[], params: bigint[]): string {
  return guard(() => fmtSeries(ls(params, coeffs, 2, 3)[1].inverse()));
}

function ls_pow(coeffs: bigint[], params: bigint[]): string {
  return guard(() => fmtSeries(ls(params, coeffs, 2, 3)[1].pow(params[4]!)));
}

function ls_derivative(coeffs: bigint[], params: bigint[]): string {
  return guard(() => fmtSeries(ls(params, coeffs, 2, 3)[1].derivative()));
}

function ls_integral(coeffs: bigint[], params: bigint[]): string {
  return guard(() => fmtSeries(ls(params, coeffs, 2, 3)[1].integral()));
}

function ls_shift(coeffs: bigint[], params: bigint[]): string {
  return guard(() => fmtSeries(ls(params, coeffs, 2, 3)[1].shift(Number(params[4]!))));
}

function ls_truncate(coeffs: bigint[], params: bigint[]): string {
  return guard(() => ls(params, coeffs, 2, 3)[1].truncate(Number(params[4]!)).toString());
}

function ls_truncate_laurentseries(coeffs: bigint[], params: bigint[]): string {
  return guard(() =>
    fmtSeries(ls(params, coeffs, 2, 3)[1].truncate_laurentseries(Number(params[4]!)))
  );
}

function ls_truncate_neg(coeffs: bigint[], params: bigint[]): string {
  return guard(() => fmtSeries(ls(params, coeffs, 2, 3)[1].truncate_neg(Number(params[4]!))));
}

function ls_add_bigoh(coeffs: bigint[], params: bigint[]): string {
  return guard(() => fmtSeries(ls(params, coeffs, 2, 3)[1].add_bigoh(Number(params[4]!))));
}

function ls_verschiebung(coeffs: bigint[], params: bigint[]): string {
  return guard(() => fmtSeries(ls(params, coeffs, 2, 3)[1].verschiebung(Number(params[4]!))));
}

function ls_reverse(coeffs: bigint[], params: bigint[]): string {
  return guard(() => {
    const f = ls(params, coeffs, 2, 3)[1];
    const precision = opt(params[4]!);
    return fmtSeries(precision === null ? f.reverse() : f.reverse(precision));
  });
}

function ls_nth_root(coeffs: bigint[], params: bigint[]): string {
  return guard(() => {
    const f = ls(params, coeffs, 2, 3)[1];
    const k = Number(params[4]!);
    const target = opt(params[5]!);
    return fmtSeries(target === null ? f.nth_root(k) : f.nth_root(k, target));
  });
}

function ls_is_square(coeffs: bigint[], params: bigint[]): string {
  return guard(() => bool(ls(params, coeffs, 2, 3)[1].is_square()));
}

function ls_power_series(coeffs: bigint[], params: bigint[]): string {
  return guard(() => fmtSeries(ls(params, coeffs, 2, 3)[1].power_series()));
}

function ls_lift_to_precision(coeffs: bigint[], params: bigint[]): string {
  return guard(() => {
    const f = ls(params, coeffs, 2, 3)[1];
    const absprec = opt(params[4]!);
    return fmtSeries(absprec === null ? f.lift_to_precision() : f.lift_to_precision(absprec));
  });
}

function ls_getitem(coeffs: bigint[], params: bigint[]): string {
  return guard(() => ls(params, coeffs, 2, 3)[1].__getitem__(Number(params[4]!)).toString());
}

function ls_ring_info(params: bigint[]): string {
  return guard(() => {
    const L = lsRing(params);
    return [
      L.toString(),
      `field=${bool(L.is_field())}`,
      `exact=${bool(L.is_exact())}`,
      `char=${L.characteristic()}`,
      `default_prec=${L.default_prec()}`,
      `gen=${L.gen()}`,
      `unif=${guard(() => L.uniformizer().toString())}`,
      `zero=${fmtSeries(L.zero())}`,
      `one=${fmtSeries(L.one())}`,
    ].join(' | ');
  });
}

export const functions = {
  // p-adic rings
  padic_ring_info,
  padic_uniformizer_pow,
  padic_teichmuller,
  padic_teichmuller_system,
  padic_roots_of_unity,
  // p-adic elements
  padic_repr,
  padic_parent,
  padic_from_int,
  padic_from_int_field,
  padic_expansion,
  padic_flags,
  padic_unit_part,
  padic_add,
  padic_sub,
  padic_mul,
  padic_div,
  padic_pow,
  padic_inv,
  padic_residue,
  padic_getitem,
  padic_slice,
  padic_add_bigoh,
  padic_lift_to_precision,
  padic_inexact_zero,
  padic_sqrt,
  padic_sqrt_all,
  padic_is_square,
  padic_nth_root,
  padic_nth_root_all,
  padic_log,
  padic_exp,
  padic_artin_hasse_exp,
  padic_multiplicative_order,
  padic_additive_order,
  // power series
  ps_repr,
  ps_info,
  ps_add,
  ps_sub,
  ps_mul,
  ps_div,
  ps_compose,
  ps_inv,
  ps_pow,
  ps_derivative,
  ps_integral,
  ps_exp,
  ps_log,
  ps_sqrt,
  ps_nth_root,
  ps_pade,
  ps_reverse,
  ps_V,
  ps_shift,
  ps_truncate,
  ps_truncate_powerseries,
  ps_add_bigoh,
  ps_O,
  ps_valuation_zero_part,
  ps_getitem,
  ps_is_square,
  // Laurent series
  ls_ring_info,
  ls_repr,
  ls_info,
  ls_add,
  ls_sub,
  ls_mul,
  ls_div,
  ls_inv,
  ls_pow,
  ls_derivative,
  ls_integral,
  ls_shift,
  ls_truncate,
  ls_truncate_laurentseries,
  ls_truncate_neg,
  ls_add_bigoh,
  ls_verschiebung,
  ls_reverse,
  ls_nth_root,
  ls_is_square,
  ls_power_series,
  ls_lift_to_precision,
  ls_getitem,
  ls_compose,
};
