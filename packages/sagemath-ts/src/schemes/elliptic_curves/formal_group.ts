/**
 * @module sage/schemes/elliptic_curves/formal_group
 * @description Formal groups of elliptic curves
 *
 * Port of: sage/schemes/elliptic_curves/formal_group.py
 *
 * This module implements formal groups associated to elliptic curves.
 * The formal group is a power series that encodes the group law
 * near the identity element.
 */

import { NotImplementedError, ValueError, ZeroDivisionError } from '../../errors.js';
import {
  type CoefficientRing,
  type LaurentSeriesElement,
  LaurentSeriesRing,
  MPowerSeries,
  MPowerSeriesRing,
  type PowerSeriesElement,
  type PowerSeriesRing,
  type RingElement,
} from '../../rings/power_series_ring.js';
import {
  type EllipticCurveInterface,
  type EllipticCurvePoint,
  type FieldElement,
  type FieldRing,
  affinePoint,
  pointAtInfinity,
} from './ell_point.js';

/**
 * Interface for an elliptic curve that provides what we need
 */
interface EllipticCurveForFormalGroup {
  base_ring: CoefficientRing<RingElement>;
  ainvs(): [RingElement, RingElement, RingElement, RingElement, RingElement];
  a1(): RingElement;
  a2(): RingElement;
  a3(): RingElement;
  a4(): RingElement;
  a6(): RingElement;
  b_invariants(): [RingElement, RingElement, RingElement, RingElement];
  toString(): string;
}

/**
 * Helper function to compute Newton method iteration sizes.
 * Given a target precision, returns a list of intermediate precisions
 * to use in Newton iteration (each roughly doubling).
 *
 * @see Reference: sage/misc/misc.py:newton_method_sizes
 */
function newton_method_sizes(n: number): number[] {
  if (n <= 1) return [1];
  const sizes: number[] = [];
  let current = n;
  while (current > 1) {
    sizes.unshift(current);
    current = Math.ceil(current / 2);
  }
  sizes.unshift(1);
  return sizes;
}

/** Convert an integer to an element of the coefficient ring. */
function ringInt(k: CoefficientRing<RingElement>, n: number | bigint): RingElement {
  return k.__call__(typeof n === 'bigint' ? n : BigInt(n));
}

// ---------------------------------------------------------------------------
// Two-variable power series: SageMath's ``PowerSeriesRing(R, 2, 't1,t2')``.
// ---------------------------------------------------------------------------

/** ``PowerSeriesRing(k, 2, 't1,t2')``, one per coefficient ring. */
const _t1t2Rings = new WeakMap<object, MPowerSeriesRing<RingElement>>();

function bivariateRing(k: CoefficientRing<RingElement>): MPowerSeriesRing<RingElement> {
  let R = _t1t2Rings.get(k as unknown as object);
  if (R === undefined) {
    R = new MPowerSeriesRing<RingElement>(k, 't1,t2');
    _t1t2Rings.set(k as unknown as object, R);
  }
  return R;
}

/**
 * A power series in the two variables ``t1``, ``t2`` over the base ring, i.e.
 * an element of SageMath's ``PowerSeriesRing(R, 2, 't1,t2')``.
 *
 * This is a thin naming layer over {@link MPowerSeries} (the port of
 * `sage/rings/multi_power_series_ring_element.py`): all the arithmetic,
 * precision bookkeeping and printing come from there.  The only things added
 * here are the two-variable spellings `terms` (coefficients keyed `"i,j"`) and
 * `coefficient(i, j)` that the elliptic-curve code and its tests use, plus
 * covariant overrides so that arithmetic on two-variable series keeps this
 * type.
 *
 * @see Reference: sage/rings/multi_power_series_ring_element.py:202 (MPowerSeries)
 * @see Reference: sage/schemes/elliptic_curves/formal_group.py:529 (group_law)
 */
export class BivariatePowerSeries extends MPowerSeries<RingElement> {
  /** Re-type an `MPowerSeries` in ``k[[t1,t2]]`` as a `BivariatePowerSeries`. */
  static of(f: MPowerSeries<RingElement>): BivariatePowerSeries {
    return new BivariatePowerSeries(f.parent(), f.monomial_coefficients(), f.prec());
  }

  /** The zero series of ``k[[t1,t2]]``, to precision `prec`. */
  static zero(k: CoefficientRing<RingElement>, prec: number): BivariatePowerSeries {
    return new BivariatePowerSeries(bivariateRing(k), [], prec);
  }

  /** The one of ``k[[t1,t2]]``, to precision `prec`. */
  static one(k: CoefficientRing<RingElement>, prec: number): BivariatePowerSeries {
    return BivariatePowerSeries.monomial(k, prec, 0, 0, k.one());
  }

  /** The monomial ``c * t1^i * t2^j``, to precision `prec`. */
  static monomial(
    k: CoefficientRing<RingElement>,
    prec: number,
    i: number,
    j: number,
    c: RingElement
  ): BivariatePowerSeries {
    return new BivariatePowerSeries(bivariateRing(k), [[[i, j], c]], prec);
  }

  /** The coefficients, keyed `` `${i},${j}` `` for the monomial ``t1^i t2^j``. */
  get terms(): Map<string, RingElement> {
    const m = new Map<string, RingElement>();
    for (const [e, c] of this.monomial_coefficients()) {
      m.set(e.join(','), c);
    }
    return m;
  }

  /** The coefficient of ``t1^i * t2^j`` (zero beyond the precision). */
  coefficient(i: number, j: number): RingElement {
    if (i + j >= this.prec()) {
      return this.base_ring().zero();
    }
    return this.__getitem__([i, j]);
  }

  override add(right: MPowerSeries<RingElement>): BivariatePowerSeries {
    return BivariatePowerSeries.of(super.add(right));
  }

  override sub(right: MPowerSeries<RingElement>): BivariatePowerSeries {
    return BivariatePowerSeries.of(super.sub(right));
  }

  override neg(): BivariatePowerSeries {
    return BivariatePowerSeries.of(super.neg());
  }

  override mul(right: MPowerSeries<RingElement>): BivariatePowerSeries {
    return BivariatePowerSeries.of(super.mul(right));
  }

  override pow(n: number | bigint): BivariatePowerSeries {
    return BivariatePowerSeries.of(super.pow(n));
  }

  /** ``~self``; see {@link bivariateInverse} for why it is not `super.inv()`. */
  override inv(): BivariatePowerSeries {
    return BivariatePowerSeries.of(bivariateInverse(this));
  }

  override add_bigoh(prec: number): BivariatePowerSeries {
    return BivariatePowerSeries.of(super.add_bigoh(prec));
  }

  /** Multiply by an element of the base ring (SageMath's ``_lmul_``). */
  scalarMul(c: RingElement): BivariatePowerSeries {
    return BivariatePowerSeries.of(this.scalar_mul(c));
  }

  /**
   * Substitute ``t1 -> A``, ``t2 -> B``.
   *
   * This is SageMath's ``F(A, B)``; when `A` and `B` do not live in the parent
   * of `self` (e.g. they are three-variable series) SageMath's `__call__`
   * falls through to `_subs_formal`, which is what is used here.
   *
   * @see Reference: sage/rings/multi_power_series_ring_element.py:512 (_subs_formal)
   */
  subs<S extends MPowerSeries<RingElement>>(A: S, B: S): S {
    return this._subs_formal(A, B) as S;
  }
}

/**
 * ``~f`` for a two-variable series `f` whose constant term is a unit.
 *
 * SageMath inverts the *background univariate* series
 * (`MPowerSeries.__invert__` is ``self.parent(~self._bg_value)``,
 * multi_power_series_ring_element.py:725), so the answer has precision exactly
 * `f.prec()`.  {@link MPowerSeries.inv} in this port sums the geometric series
 * `sum (-z)^m` instead, and `MPowerSeries.mul`'s precision rule
 * `min(p1+v2, p2+v1)` makes the precision of the running power grow by
 * `val(z)` at every step -- so it never becomes zero and the loop runs the full
 * `prec` iterations on ever larger series.  Truncating each step at `f.prec()`
 * restores both upstream's precision and upstream's cost (`group_law(50)` over
 * `GF(17)`: 12 s -> 0.2 s), and drops no coefficient the answer keeps.
 */
function bivariateInverse(f: MPowerSeries<RingElement>): MPowerSeries<RingElement> {
  const R = f.parent();
  const prec = f.prec() === Number.POSITIVE_INFINITY ? R.default_prec() : f.prec();
  const zeroExp = new Array(R.ngens()).fill(0) as number[];
  const c = f.__getitem__(zeroExp);
  if (c.isZero()) {
    throw new ZeroDivisionError('cannot invert a series with zero constant term');
  }
  const cinv = c.inv ? c.inv() : R.base_ring().one().div(c);
  const one = R.one().add_bigoh(prec);
  const z = f.add_bigoh(prec).scalar_mul(cinv).sub(one).neg();
  let result = one;
  let zp = one;
  for (let m = 1; m < prec; m++) {
    zp = zp.mul(z).add_bigoh(prec);
    if (zp.is_zero()) break;
    result = result.add(zp);
  }
  return result.scalar_mul(cinv).add_bigoh(prec);
}

/**
 * Evaluate the univariate power series `f` at the multivariate series `a`.
 *
 * SageMath writes this as plain function application (``w(t1)`` and
 * ``inv(t3)`` in ``group_law``); the coercion framework sends it to
 * `PowerSeries_poly.__call__`.  That routine — ported verbatim as
 * {@link PowerSeriesElement.__call__} — only touches the part of the interface
 * that {@link MPowerSeries} also implements (`parent`, `valuation`, `prec`,
 * `add_bigoh`, `mul`, `add`), so it evaluates a multivariate argument
 * unchanged, including its truncation of the argument to `(s - r + 1) t`.
 * The cast exists only because TypeScript types the parameter univariately.
 *
 * @see Reference: sage/rings/power_series_poly.pyx:176 (__call__)
 */
function evalAt(
  f: PowerSeriesElement<RingElement>,
  a: MPowerSeries<RingElement>
): MPowerSeries<RingElement> {
  return (f as unknown as { __call__(x: unknown): MPowerSeries<RingElement> }).__call__(a);
}

// ---------------------------------------------------------------------------
// The base change E -> E x_R R((t)) used by ``mult_by_n`` in characteristic 0.
// ---------------------------------------------------------------------------

/**
 * A {@link LaurentSeriesElement} presented as an elliptic-curve `FieldElement`.
 *
 * SageMath's characteristic-zero `mult_by_n` base-changes `E` to the Laurent
 * series ring `R = x.parent()` and runs the curve's own group law there.  The
 * port's `LaurentSeriesElement` (from `laurent_series_ring_element.pyx`) is
 * already the field arithmetic that needs; the elliptic-curve layer merely
 * spells two things differently (`parent` as a property rather than a method,
 * `isZero` rather than `is_zero`).  This wrapper renames them and does no
 * arithmetic of its own.
 *
 * @see Reference: sage/schemes/elliptic_curves/formal_group.py:655 (mult_by_n)
 */
class LaurentFieldElement implements FieldElement {
  readonly parent: LaurentFieldRing;
  readonly series: LaurentSeriesElement<RingElement>;

  constructor(parent: LaurentFieldRing, series: LaurentSeriesElement<RingElement>) {
    this.parent = parent;
    this.series = series;
  }

  private _s(other: FieldElement | number | bigint): LaurentSeriesElement<RingElement> {
    if (other instanceof LaurentFieldElement) {
      return other.series;
    }
    return this.parent.__call__(other as number | bigint).series;
  }

  add(other: FieldElement | number | bigint): LaurentFieldElement {
    return new LaurentFieldElement(this.parent, this.series.add(this._s(other)));
  }

  sub(other: FieldElement | number | bigint): LaurentFieldElement {
    return new LaurentFieldElement(this.parent, this.series.sub(this._s(other)));
  }

  mul(other: FieldElement | number | bigint): LaurentFieldElement {
    return new LaurentFieldElement(this.parent, this.series.mul(this._s(other)));
  }

  div(other: FieldElement | number | bigint): LaurentFieldElement {
    return new LaurentFieldElement(this.parent, this.series.div(this._s(other)));
  }

  neg(): LaurentFieldElement {
    return new LaurentFieldElement(this.parent, this.series.neg());
  }

  inv(): LaurentFieldElement {
    return new LaurentFieldElement(this.parent, this.series.inv());
  }

  pow(n: bigint | number): LaurentFieldElement {
    return new LaurentFieldElement(this.parent, this.series.pow(n));
  }

  isZero(): boolean {
    return this.series.is_zero();
  }

  eq(other: FieldElement): boolean {
    return this.series.eq(this._s(other));
  }

  toString(): string {
    return this.series.toString();
  }
}

/** A {@link LaurentSeriesRing} presented as an elliptic-curve `FieldRing`. */
class LaurentFieldRing implements FieldRing {
  readonly series_ring: LaurentSeriesRing<RingElement>;
  readonly characteristic: bigint;

  constructor(series_ring: LaurentSeriesRing<RingElement>) {
    this.series_ring = series_ring;
    this.characteristic = series_ring.characteristic();
  }

  /** Wrap a Laurent series as a field element of this ring. */
  fromSeries(s: LaurentSeriesElement<RingElement>): LaurentFieldElement {
    return new LaurentFieldElement(this, s);
  }

  zero(): LaurentFieldElement {
    return new LaurentFieldElement(this, this.series_ring.zero());
  }

  one(): LaurentFieldElement {
    return new LaurentFieldElement(this, this.series_ring.one());
  }

  __call__(value: bigint | number | FieldElement): LaurentFieldElement {
    if (value instanceof LaurentFieldElement) {
      return value;
    }
    return new LaurentFieldElement(this, this.series_ring.__call__(value));
  }

  toString(): string {
    return this.series_ring.toString();
  }
}

/**
 * ``E.change_ring(R)`` for `R` a Laurent series ring: the same Weierstrass
 * equation, with the `a`-invariants read as constant Laurent series.
 *
 * Its points are ordinary {@link EllipticCurvePoint}s, so `n*P` below is the
 * port of SageMath's own group law (`ell_point.py:_add_`, `_neg_` and the
 * double-and-add of `IntegerMulAction`), which is exactly what upstream's
 * characteristic-zero `mult_by_n` uses.
 *
 * @see Reference: sage/schemes/elliptic_curves/ell_generic.py:change_ring
 */
class EllipticCurveOverLaurentSeries implements EllipticCurveInterface<LaurentFieldElement> {
  readonly base_ring: LaurentFieldRing;
  private readonly _ainvs: [
    LaurentFieldElement,
    LaurentFieldElement,
    LaurentFieldElement,
    LaurentFieldElement,
    LaurentFieldElement,
  ];
  private _infinity: EllipticCurvePoint<LaurentFieldElement> | null = null;

  constructor(
    K: LaurentFieldRing,
    ainvs: [
      LaurentFieldElement,
      LaurentFieldElement,
      LaurentFieldElement,
      LaurentFieldElement,
      LaurentFieldElement,
    ]
  ) {
    this.base_ring = K;
    this._ainvs = ainvs;
  }

  a1(): LaurentFieldElement {
    return this._ainvs[0];
  }

  a2(): LaurentFieldElement {
    return this._ainvs[1];
  }

  a3(): LaurentFieldElement {
    return this._ainvs[2];
  }

  a4(): LaurentFieldElement {
    return this._ainvs[3];
  }

  a6(): LaurentFieldElement {
    return this._ainvs[4];
  }

  a_invariants(): [
    LaurentFieldElement,
    LaurentFieldElement,
    LaurentFieldElement,
    LaurentFieldElement,
    LaurentFieldElement,
  ] {
    return [...this._ainvs];
  }

  ainvs(): [
    LaurentFieldElement,
    LaurentFieldElement,
    LaurentFieldElement,
    LaurentFieldElement,
    LaurentFieldElement,
  ] {
    return this.a_invariants();
  }

  zero(): EllipticCurvePoint<LaurentFieldElement> {
    if (this._infinity === null) {
      this._infinity = pointAtInfinity<LaurentFieldElement>(this);
    }
    return this._infinity;
  }

  /** ``y^2 + a1 x y + a3 y - (x^3 + a2 x^2 + a4 x + a6) == 0``. */
  is_on_curve(x: LaurentFieldElement, y: LaurentFieldElement): boolean {
    const [a1, a2, a3, a4, a6] = this._ainvs;
    const lhs = y.mul(y).add(a1.mul(x).mul(y)).add(a3.mul(y));
    const rhs = x.mul(x).mul(x).add(a2.mul(x).mul(x)).add(a4.mul(x)).add(a6);
    return lhs.sub(rhs).isZero();
  }

  toString(): string {
    return `Elliptic Curve over ${this.base_ring}`;
  }
}

/**
 * The formal group associated to an elliptic curve.
 *
 * The formal group F of an elliptic curve E is a one-dimensional
 * formal group law F(X, Y) = X + Y + higher order terms that encodes
 * the addition law on E near the origin.
 *
 * @see Reference: sage/schemes/elliptic_curves/formal_group.py:EllipticCurveFormalGroup
 */
export class EllipticCurveFormalGroup {
  private readonly _E: EllipticCurveForFormalGroup;

  // Cached computed values
  private _cachedW: { prec: number; value: PowerSeriesElement<RingElement> } | null = null;
  private _cachedY: { prec: number; value: LaurentSeriesElement<RingElement> } | null = null;
  private _cachedOmega: { prec: number; value: PowerSeriesElement<RingElement> } | null = null;
  private _cachedInverse: { prec: number; value: PowerSeriesElement<RingElement> } | null = null;
  private _cachedGroupLaw: { prec: number; value: BivariatePowerSeries } | null = null;
  private _lsRingCache: LaurentSeriesRing<RingElement> | null = null;

  /**
   * Create the formal group for an elliptic curve.
   *
   * INPUT:
   * - E: an elliptic curve
   *
   * @see Reference: sage/schemes/elliptic_curves/formal_group.py:__init__
   */
  constructor(E: EllipticCurveForFormalGroup) {
    this._E = E;
  }

  /**
   * Return the elliptic curve this formal group is associated to.
   *
   * @see Reference: sage/schemes/elliptic_curves/formal_group.py:curve
   */
  curve(): EllipticCurveForFormalGroup {
    return this._E;
  }

  /**
   * Return the formal group power series w(t).
   *
   * This is the formal expansion of w = -1/y about the formal parameter
   * t = -x/y at infinity:
   *
   *   w(t) = t^3 + a_1*t^4 + (a_2 + a_1^2)*t^5 + ...
   *
   * INPUT:
   * - prec: precision (default 20)
   *
   * OUTPUT: a power series w(t) to precision O(t^prec)
   *
   * ALGORITHM: Uses Newton's method to solve the elliptic curve
   * equation at the origin.
   *
   * @see Reference: sage/schemes/elliptic_curves/formal_group.py:w
   */
  w(prec: number = 20): PowerSeriesElement<RingElement> {
    prec = Math.max(prec, 0);
    const k = this._E.base_ring;
    const R = this._psRing();

    // Check cache
    if (this._cachedW !== null && prec <= this._cachedW.prec) {
      return this._cachedW.value.add_bigoh(prec);
    }

    // Initialize: w = t^3 (to precision 4)
    let w = R.__call__([k.zero(), k.zero(), k.zero(), k.one()], 4);
    let currentPrec = 4;

    if (prec < currentPrec) {
      return w.add_bigoh(prec);
    }

    // Newton iteration to compute w(t) satisfying:
    // w = t^3 * (1 + a_1*t + (a_2 + a_1^2)*t^2 + ...)
    // From the Weierstrass equation, we have:
    // w^(-1) + a_1*t*w^(-1) + a_3*w = t^(-2)*w^(-1) + a_2*t^(-1)*w^(-1) + a_4*w + a_6*w^2
    // Rearranging: w = t^3 / (1 - a_1*t - a_2*t^2 - 2*a_3*w - 2*a_4*t*w - 3*a_6*w^2)
    // We iterate using Newton's method

    const [a1, a2, a3, a4, a6] = this._E.ainvs();

    // Constants for the iteration
    const t3_poly = R.__call__([k.zero(), k.zero(), k.zero(), k.one()], Number.POSITIVE_INFINITY);
    const const_poly_1 = R.__call__([k.one()], Number.POSITIVE_INFINITY);
    const neg_a1_t = R.__call__([k.zero(), a1.neg()], Number.POSITIVE_INFINITY);
    const neg_a2_t2 = R.__call__([k.zero(), k.zero(), a2.neg()], Number.POSITIVE_INFINITY);

    const sizes = newton_method_sizes(prec);

    for (const nextPrec of sizes) {
      if (nextPrec <= currentPrec) continue;

      // Current w truncated appropriately
      const wTrunc = w.truncate(nextPrec);

      // Compute w^2 and w^3
      const w2 = wTrunc.mul(wTrunc).add_bigoh(nextPrec);

      // Numerator: t^3 - a_3*w^2 - a_4*t*w^2 - 2*a_6*w^3
      const w3 = w2.mul(wTrunc).add_bigoh(nextPrec);
      const a3_w2 = w2
        ._shiftLeft(0)
        .mul(R.__call__([a3], Number.POSITIVE_INFINITY))
        .add_bigoh(nextPrec);
      const a4_t_w2 = w2
        ._shiftLeft(1)
        .mul(R.__call__([a4], Number.POSITIVE_INFINITY))
        .add_bigoh(nextPrec);
      const two_a6_w3 = w3
        .mul(R.__call__([a6.add(a6)], Number.POSITIVE_INFINITY))
        .add_bigoh(nextPrec);
      const numerator = t3_poly.add_bigoh(nextPrec).sub(a3_w2).sub(a4_t_w2).sub(two_a6_w3);

      // Denominator: 1 - a_1*t - a_2*t^2 - 2*a_3*w - 2*a_4*t*w - 3*a_6*w^2
      const two_a3_w = wTrunc
        .mul(R.__call__([a3.add(a3)], Number.POSITIVE_INFINITY))
        .add_bigoh(nextPrec);
      const two_a4_t_w = wTrunc
        ._shiftLeft(1)
        .mul(R.__call__([a4.add(a4)], Number.POSITIVE_INFINITY))
        .add_bigoh(nextPrec);
      const three_a6 = a6.add(a6).add(a6);
      const three_a6_w2 = w2
        .mul(R.__call__([three_a6], Number.POSITIVE_INFINITY))
        .add_bigoh(nextPrec);
      const denominator = const_poly_1
        .add_bigoh(nextPrec)
        .add(neg_a1_t.add_bigoh(nextPrec))
        .add(neg_a2_t2.add_bigoh(nextPrec))
        .sub(two_a3_w)
        .sub(two_a4_t_w)
        .sub(three_a6_w2);

      // w = numerator / denominator
      const inv = denominator.inv();
      w = numerator.mul(inv).add_bigoh(nextPrec);
      currentPrec = nextPrec;
    }

    // Cache and return
    this._cachedW = { prec: currentPrec, value: w };
    return w.add_bigoh(prec);
  }

  /**
   * Return the formal series x(t) = t/w(t) in terms of the local parameter
   * t = -x/y at infinity:
   *
   *   x(t) = t^(-2) - a_1*t^(-1) - a_2 - a_3*t - ...
   *
   * INPUT:
   * - prec: precision (default 20)
   *
   * OUTPUT: a Laurent series x(t) to precision O(t^prec)
   *
   * @see Reference: sage/schemes/elliptic_curves/formal_group.py:233 (x)
   */
  x(prec: number = 20): LaurentSeriesElement<RingElement> {
    prec = Math.max(prec, 0);
    // Sage: y = self.y(prec); t = y.parent().gen(); return -t*y + O(t**prec)
    const y = this.y(prec);
    const t = this._lsRing().gen();
    return t.neg().mul(y).add_bigoh(prec);
  }

  /**
   * Return the formal series y(t) = -1/w(t) in terms of the local parameter
   * t = -x/y at infinity:
   *
   *   y(t) = -t^(-3) + a_1*t^(-2) + a_2*t^(-1) + a_3 + ...
   *
   * @see Reference: sage/schemes/elliptic_curves/formal_group.py:268 (y)
   */
  y(prec: number = 20): LaurentSeriesElement<RingElement> {
    prec = Math.max(prec, 0);
    if (this._cachedY !== null && prec <= this._cachedY.prec) {
      return this._cachedY.value.add_bigoh(prec);
    }
    // Sage: w = self.w(prec+6); t = w.parent().gen(); y = -(w**(-1)) + O(t**prec)
    const w = this.w(prec + 6);
    const L = this._lsRing();
    const y = L.__call__(w).inv().neg().add_bigoh(prec);
    this._cachedY = { prec, value: y };
    return y;
  }

  /**
   * Return the coefficients of ``x(t)`` from its valuation ``-2`` up to
   * ``t^(prec-1)``, i.e. ``[c_{-2}, c_{-1}, c_0, ..., c_{prec-1}]``.
   *
   * In Sage this is ``E.formal_group().x(prec).list()``; the Laurent series
   * `list()` starts at the valuation, so this is the same list.  It is kept as
   * a named accessor because it predates the Laurent series port.
   *
   * @see Reference: sage/rings/laurent_series_ring_element.pyx:571 (list)
   */
  x_list(prec: number = 20): RingElement[] {
    prec = Math.max(prec, 0);
    const x = this.x(prec);
    const out: RingElement[] = [];
    for (let n = x.valuation(); n < prec; n++) {
      out.push(x.__getitem__(n));
    }
    return out;
  }

  /**
   * Return the coefficients of ``y(t)`` from its valuation ``-3`` up to
   * ``t^(prec-1)``, i.e. ``[c_{-3}, c_{-2}, c_{-1}, c_0, ..., c_{prec-1}]``.
   *
   * @see Reference: sage/rings/laurent_series_ring_element.pyx:571 (list)
   */
  y_list(prec: number = 20): RingElement[] {
    prec = Math.max(prec, 0);
    const y = this.y(prec);
    const out: RingElement[] = [];
    for (let n = y.valuation(); n < prec; n++) {
      out.push(y.__getitem__(n));
    }
    return out;
  }

  /**
   * Return the power series f(t) = 1 + ... such that f(t) dt is the usual
   * invariant differential dx/(2y + a_1 x + a_3).
   *
   * @see Reference: sage/schemes/elliptic_curves/formal_group.py:315 (differential)
   */
  differential(prec: number = 20): PowerSeriesElement<RingElement> {
    prec = Math.max(prec, 0);

    if (this._cachedOmega !== null && prec <= this._cachedOmega.prec) {
      return this._cachedOmega.value.add_bigoh(prec);
    }

    const k = this._E.base_ring;
    const a = this._E.ainvs();
    const L = this._lsRing();

    // Sage: x = self.x(prec+1); y = self.y(prec+1); xprime = x.derivative()
    //       g = xprime / (2*y + a[0]*x + a[2])
    const x = this.x(prec + 1);
    const y = this.y(prec + 1);
    const xprime = x.derivative();
    const denom = y.scalar_mul(ringInt(k, 2)).add(x.scalar_mul(a[0])).add(L.__call__(a[2]));

    const g = xprime.div(denom).power_series().add_bigoh(prec);
    this._cachedOmega = { prec, value: g };
    return g;
  }

  /**
   * Return the power series f(t) = t + ... which is an isomorphism to the
   * additive formal group (the formal logarithm).
   *
   * Generally this only makes sense in characteristic zero, although the
   * terms before t^p may work in characteristic p.
   *
   * @see Reference: sage/schemes/elliptic_curves/formal_group.py:372 (log)
   */
  log(prec: number = 20): PowerSeriesElement<RingElement> {
    return this.differential(prec - 1)
      .integral()
      .add_bigoh(prec);
  }

  /**
   * Return the formal group inverse law i(t), which satisfies F(t, i(t)) = 0.
   *
   * i(t) = -t + a1*t^2 + ... to precision O(t^prec) (page 114 of [Sil2009]).
   *
   * @see Reference: sage/schemes/elliptic_curves/formal_group.py:397 (inverse)
   */
  inverse(prec: number = 20): PowerSeriesElement<RingElement> {
    prec = Math.max(prec, 0);

    if (this._cachedInverse !== null && prec <= this._cachedInverse.prec) {
      return this._cachedInverse.value.add_bigoh(prec);
    }

    const [a1, , a3] = this._E.ainvs();
    const L = this._lsRing();

    // Sage: inv = x / (y + a1*x + a3)   (page 114 of Silverman, AEC I)
    const x = this.x(prec);
    const y = this.y(prec);
    const denom = y.add(x.scalar_mul(a1)).add(L.__call__(a3));
    const inv = x.div(denom).power_series().add_bigoh(prec);
    this._cachedInverse = { prec, value: inv };
    return inv;
  }

  /**
   * Return the formal group law F(t1, t2) = t1 + t2 - a1*t1*t2 - ... to
   * precision O(t1, t2)^prec (page 115 of [Sil2009]).
   *
   * @see Reference: sage/schemes/elliptic_curves/formal_group.py:450 (group_law)
   */
  group_law(prec: number = 10): BivariatePowerSeries {
    prec = Math.max(prec, 0);
    if (prec <= 0) {
      throw new ValueError('The precision must be positive.');
    }

    const k = this._E.base_ring;
    const [a1, a2, a3, a4, a6] = this._E.ainvs();

    // Sage: R = PowerSeriesRing(self.curve().base_ring(), 2, 't1,t2')
    const R = bivariateRing(k);
    const [t1, t2] = R.gens() as [MPowerSeries<RingElement>, MPowerSeries<RingElement>];

    if (prec === 1) {
      // Sage: return R(0)
      return BivariatePowerSeries.of(R.zero());
    }
    if (prec === 2) {
      // Sage: return t1 + t2 - self.curve().a1()*t1*t2
      return BivariatePowerSeries.of(t1.add(t2).sub(t1.mul(t2).scalar_mul(a1)));
    }

    if (this._cachedGroupLaw !== null && prec <= this._cachedGroupLaw.prec) {
      return this._cachedGroupLaw.value.add_bigoh(prec);
    }

    const w = this.w(prec + 1);

    // Sage: lam = sum([w[n]*sum(t2**m * t1**(n-m-1) for m in range(n))
    //                  for n in range(3, prec+1)]); lam = lam.add_bigoh(prec)
    let lam: MPowerSeries<RingElement> = R.zero();
    for (let n = 3; n <= prec; n++) {
      const c = w.__getitem__(n);
      if (c.isZero()) continue;
      for (let m = 0; m < n; m++) {
        lam = lam.add(new MPowerSeries<RingElement>(R, [[[n - m - 1, m], c]]));
      }
    }
    lam = lam.add_bigoh(prec);

    // Sage: nu = w(t1) - lam*t1
    const nu = evalAt(w, t1).sub(lam.mul(t1));

    const lam2 = lam.mul(lam);
    const lam3 = lam2.mul(lam);

    // Note that the following formula differs from the one in Silverman
    // page 119.  See github issue 9646 for the explanation and justification.
    // Sage: t3 = -t1 - t2
    //            - (a1*lam + a3*lam2 + a2*nu + 2*a4*lam*nu + 3*a6*lam2*nu)
    //              / (1 + a2*lam + a4*lam2 + a6*lam3)
    const numer = lam
      .scalar_mul(a1)
      .add(lam2.scalar_mul(a3))
      .add(nu.scalar_mul(a2))
      .add(lam.mul(nu).scalar_mul(a4.mul(ringInt(k, 2))))
      .add(lam2.mul(nu).scalar_mul(a6.mul(ringInt(k, 3))));
    const denom = R.one().add(lam.scalar_mul(a2)).add(lam2.scalar_mul(a4)).add(lam3.scalar_mul(a6));
    const t3 = t1
      .neg()
      .sub(t2)
      .sub(numer.mul(bivariateInverse(denom)));

    // Sage: inv = self.inverse(prec); F = inv(t3).add_bigoh(prec)
    const inv = this.inverse(prec);
    const F = BivariatePowerSeries.of(evalAt(inv, t3).add_bigoh(prec));

    this._cachedGroupLaw = { prec, value: F };
    return F;
  }

  /**
   * Return the formal "multiplication by n" endomorphism [n](t) = n*t + ...
   * to precision O(t^prec) (Proposition 2.3 of [Sil2009]).
   *
   * ALGORITHM: exactly upstream's two branches.  Over a field of
   * characteristic zero the formal point ``(x(t), y(t))`` is formed on ``E``
   * base-changed to the Laurent series ring and multiplied by ``n`` with the
   * curve's own group law, giving ``[n](t) = -Q[0]/Q[1]``; otherwise the
   * general double-and-add with the formal group law is used.
   *
   * @see Deviation: upstream's characteristic-zero branch returns an element of
   *      the Laurent series ring (of valuation 1); this port converts it to the
   *      power series it is, so that both branches share one return type.  Same
   *      coefficients, same precision.
   * @see Reference: sage/schemes/elliptic_curves/formal_group.py:562 (mult_by_n)
   */
  mult_by_n(n: bigint | number, prec: number = 10): PowerSeriesElement<RingElement> {
    let nVal = typeof n === 'number' ? BigInt(n) : n;
    const R = this._psRing();
    const t = R.gen();

    // Sage: if self.curve().base_ring().is_field()
    //          and self.curve().base_ring().characteristic() == 0 and n != 0
    if (nVal !== 0n && this._baseIsCharacteristicZeroField()) {
      // The following algorithm only works over a field of characteristic
      // zero.  It is much faster than using the formal group law. -- dmharvey

      // Create a "formal point" on the original curve E.  Our answer only
      // needs prec-1 coefficients (since the lowest term is t^1), and
      // x(t) = t^(-2) + ... and y(t) = t^(-3) + ..., so we only need
      // x(t) mod t^(prec-3) and y(t) mod t^(prec-4).
      const x = this.x(prec - 3);
      const y = this.y(prec - 4);
      const K = new LaurentFieldRing(x.parent()); // the Laurent series ring over the base ring
      const X = this._change_ring_to_laurent(K); // self.curve().change_ring(R)
      const P = affinePoint<LaurentFieldElement>(X, K.fromSeries(x), K.fromSeries(y));

      // and multiply it by n, using the group law on E
      const Q = P.mul(nVal);

      // express it in terms of the formal parameter
      if (Q.is_zero()) {
        // Unreachable over a field of characteristic zero: the formal point is
        // of infinite order.  Guarded so a wrong answer can never be returned.
        throw new NotImplementedError(
          'SAGE_NOT_IMPLEMENTED: mult_by_n: n*P is the point at infinity of the formal point'
        );
      }
      const [qx, qy] = [Q.x(), Q.y()];
      return qx.neg().div(qy).series.power_series();
    }

    // Now the general case, not necessarily over a field.

    if (nVal === 1n) {
      return t.add_bigoh(prec);
    }
    if (nVal === 0n) {
      return R.zero().add_bigoh(prec);
    }
    if (nVal === -1n) {
      return this.inverse(prec);
    }
    if (nVal < 0n) {
      return this.inverse(prec).__call__(this.mult_by_n(-nVal, prec));
    }

    const F = this.group_law(prec);
    /**
     * ``F(A, B)`` for univariate power series `A`, `B`: in SageMath the
     * arguments do not coerce into `F`'s parent, so `MPowerSeries.__call__`
     * falls through to `_subs_formal`, which accumulates
     * ``sum_m c_m A^{m_1} B^{m_2}`` in the parent of its first argument and
     * finishes with ``add_bigoh(F.prec())``.
     *
     * This is that sum, with each intermediate product truncated at `prec` as
     * well.  That is harmless — the answer is truncated at `prec` anyway — and
     * it keeps the powers `A^i`, `B^j` from growing past `prec` coefficients,
     * which is what makes ``mult_by_n(10, 50)`` finish in a couple of seconds
     * rather than a couple of dozen.
     *
     * @see Reference: sage/rings/multi_power_series_ring_element.py:512 (_subs_formal)
     */
    const applyF = (
      A: PowerSeriesElement<RingElement>,
      B: PowerSeriesElement<RingElement>
    ): PowerSeriesElement<RingElement> => {
      let result = R.zero().add_bigoh(prec);
      const Apows: PowerSeriesElement<RingElement>[] = [R.one().add_bigoh(prec)];
      const Bpows: PowerSeriesElement<RingElement>[] = [R.one().add_bigoh(prec)];
      for (const [e, c] of F.monomial_coefficients()) {
        const [i, j] = e as [number, number];
        while (Apows.length <= i) Apows.push(Apows[Apows.length - 1]!.mul(A).add_bigoh(prec));
        while (Bpows.length <= j) Bpows.push(Bpows[Bpows.length - 1]!.mul(B).add_bigoh(prec));
        result = result.add(Apows[i]!.mul(Bpows[j]!)._scalarMul(c)).add_bigoh(prec);
      }
      return result;
    };

    let result: PowerSeriesElement<RingElement>;
    if (nVal < 4n) {
      result = t.add_bigoh(prec);
      for (let m = 1n; m < nVal; m++) {
        result = applyF(result, t.add_bigoh(prec));
      }
      return result;
    }

    // Double and add is faster than the naive method when n >= 4.
    let g = t.add_bigoh(prec);
    result = (nVal & 1n) === 1n ? g : R.zero().add_bigoh(prec);
    nVal >>= 1n;
    while (nVal > 0n) {
      g = applyF(g, g);
      if ((nVal & 1n) === 1n) {
        result = applyF(result, g);
      }
      nVal >>= 1n;
    }
    return result;
  }

  /**
   * Return the Weierstrass sigma function as a formal power series solution of
   * d^2 log(sigma)/dz^2 = -wp(z), expressed in the formal-group parameter t.
   *
   * @see Reference: sage/schemes/elliptic_curves/formal_group.py:708 (sigma)
   */
  sigma(prec: number = 10): PowerSeriesElement<RingElement> {
    const k = this._E.base_ring;
    const [a1, a2] = this._E.ainvs();
    const T = this._psRing();

    const fl = this.log(prec);
    const F0 = fl.reversion(prec);

    // Sage: S = LaurentSeriesRing(k,'z'); z = S.gen(); F = F(z + O(z**prec))
    const S = new LaurentSeriesRing<RingElement>(k, 'z');
    const z = S.gen();
    const F = S.__call__(F0).__call__(z.add_bigoh(prec));

    // Sage: wp = self.x()(F) + (a1**2 + 4*a2)/12
    const twelve = ringInt(k, 12);
    if (twelve.isZero()) {
      throw new ZeroDivisionError('sigma requires 12 to be invertible in the base ring');
    }
    const shift = a1
      .mul(a1)
      .add(a2.mul(ringInt(k, 4)))
      .div(twelve);
    const wp = S.__call__(this.x().__call__(F)).add(S.__call__(shift));

    // Sage: g = (1/z**2 - wp).power_series(); h = g.integral().integral()
    const g = S.one().div(z.mul(z)).sub(wp).power_series();
    const h = g.integral().integral();

    // Sage: sigma_of_z = z.power_series() * h.exp()
    const sigmaOfZ = h.exp().mul(S.__call__(z).power_series());

    // Sage: fl = fl(T.gen()+O(T.gen()**prec)); sigma_of_t = sigma_of_z(fl)
    const flT = fl.__call__(T.gen().add_bigoh(prec));
    return sigmaOfZ.__call__(flT);
  }

  /**
   * Whether the base ring is a field of characteristic zero, i.e. whether
   * upstream's fast `mult_by_n` branch applies.
   *
   * The port spells `characteristic` as a method on some rings (`QQ`) and as a
   * property on others (`GF(p)`), so both are accepted here.
   */
  private _baseIsCharacteristicZeroField(): boolean {
    const k = this._E.base_ring as CoefficientRing<RingElement> & {
      is_field?: () => boolean;
      characteristic?: bigint | (() => bigint);
    };
    if (typeof k.is_field !== 'function' || !k.is_field()) {
      return false;
    }
    const ch = typeof k.characteristic === 'function' ? k.characteristic() : k.characteristic;
    return ch === 0n;
  }

  /**
   * ``self.curve().change_ring(R)`` for `R` the Laurent series ring, i.e. the
   * same Weierstrass equation with constant Laurent series coefficients.
   */
  private _change_ring_to_laurent(K: LaurentFieldRing): EllipticCurveOverLaurentSeries {
    const R = K.series_ring;
    const a = this._E.ainvs().map((ai) => K.fromSeries(R.__call__(ai))) as [
      LaurentFieldElement,
      LaurentFieldElement,
      LaurentFieldElement,
      LaurentFieldElement,
      LaurentFieldElement,
    ];
    return new EllipticCurveOverLaurentSeries(K, a);
  }

  /** The Laurent series ring R((t)) over the base ring. */
  private _lsRing(): LaurentSeriesRing<RingElement> {
    if (this._lsRingCache === null) {
      this._lsRingCache = new LaurentSeriesRing<RingElement>(this._E.base_ring, 't');
    }
    return this._lsRingCache;
  }

  /** The power series ring R[[t]] over the base ring. */
  private _psRing(): PowerSeriesRing<RingElement> {
    return this._lsRing().power_series_ring();
  }

  /**
   * Check equality with another formal group.
   */
  equals(other: EllipticCurveFormalGroup): boolean {
    if (!(other instanceof EllipticCurveFormalGroup)) {
      return false;
    }
    // Two formal groups are equal if their curves are equal
    return this._E === other._E;
  }

  /**
   * String representation.
   */
  toString(): string {
    return `Formal Group associated to the ${this._E}`;
  }
}
