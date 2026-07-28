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
  LaurentSeriesElement,
  LaurentSeriesRing,
  type PowerSeriesElement,
  PowerSeriesRing,
  type RingElement,
} from '../../rings/power_series_ring.js';

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

/**
 * A Laurent series represented as ``t^v * s(t)`` with ``s`` a power series.
 *
 * The port's LaurentSeriesElement has no arithmetic, so the formal-group code
 * carries the shift explicitly; this mirrors what Sage does inside its Laurent
 * series ring.
 */
interface Lau {
  s: PowerSeriesElement<RingElement>;
  v: number;
}

/** Convert an integer to an element of the coefficient ring. */
function ringInt(k: CoefficientRing<RingElement>, n: number | bigint): RingElement {
  return k.__call__(typeof n === 'bigint' ? n : BigInt(n));
}

function lauAdd(A: Lau, B: Lau): Lau {
  const v = Math.min(A.v, B.v);
  const a = A.s._shiftLeft(A.v - v);
  const b = B.s._shiftLeft(B.v - v);
  return { s: a.add(b), v };
}

function lauSub(A: Lau, B: Lau): Lau {
  return lauAdd(A, { s: B.s.neg(), v: B.v });
}

function lauMul(A: Lau, B: Lau): Lau {
  return { s: A.s.mul(B.s), v: A.v + B.v };
}

function lauScalar(A: Lau, c: RingElement): Lau {
  return { s: A.s._scalarMul(c), v: A.v };
}

function lauDiv(A: Lau, B: Lau): Lau {
  const bv = B.s.valuation();
  const unit = B.s._shiftRight(bv);
  return { s: A.s.mul(unit.inv()), v: A.v - B.v - bv };
}

/** d/dt of t^v * s(t) = t^(v-1) * (v*s + t*s'). */
function lauDeriv(A: Lau, k: CoefficientRing<RingElement>): Lau {
  const vEl = ringInt(k, A.v);
  return { s: A.s._scalarMul(vEl).add(A.s.derivative()._shiftLeft(1)), v: A.v - 1 };
}

/** Convert a Laurent series with nonnegative valuation to a power series. */
function lauToPowerSeries(A: Lau): PowerSeriesElement<RingElement> {
  if (A.v >= 0) {
    return A.s._shiftLeft(A.v);
  }
  const drop = -A.v;
  for (let i = 0; i < drop; i++) {
    if (!A.s.__getitem__(i).isZero()) {
      throw new ValueError('Laurent series has a pole; cannot convert to a power series');
    }
  }
  return A.s._shiftRight(drop);
}

/**
 * A truncated power series in two variables t1, t2 over the coefficient ring,
 * with total-degree truncation: everything of total degree >= prec is dropped.
 *
 * This stands in for Sage's ``PowerSeriesRing(R, 2, 't1,t2')``.
 *
 * @see Reference: sage/rings/multi_power_series_ring_element.py
 */
export class BivariatePowerSeries {
  readonly k: CoefficientRing<RingElement>;
  readonly prec: number;
  /** key `${i},${j}` -> coefficient of t1^i * t2^j */
  readonly terms: Map<string, RingElement>;

  constructor(
    k: CoefficientRing<RingElement>,
    prec: number,
    terms: Map<string, RingElement> = new Map()
  ) {
    this.k = k;
    this.prec = prec;
    this.terms = new Map();
    for (const [key, c] of terms) {
      if (c.isZero()) continue;
      const [i, j] = key.split(',').map(Number) as [number, number];
      if (i + j >= prec) continue;
      this.terms.set(key, c);
    }
  }

  static zero(k: CoefficientRing<RingElement>, prec: number): BivariatePowerSeries {
    return new BivariatePowerSeries(k, prec);
  }

  static one(k: CoefficientRing<RingElement>, prec: number): BivariatePowerSeries {
    return BivariatePowerSeries.monomial(k, prec, 0, 0, k.one());
  }

  static monomial(
    k: CoefficientRing<RingElement>,
    prec: number,
    i: number,
    j: number,
    c: RingElement
  ): BivariatePowerSeries {
    const m = new Map<string, RingElement>();
    m.set(`${i},${j}`, c);
    return new BivariatePowerSeries(k, prec, m);
  }

  /** The coefficient of t1^i * t2^j. */
  coefficient(i: number, j: number): RingElement {
    return this.terms.get(`${i},${j}`) ?? this.k.zero();
  }

  add(other: BivariatePowerSeries): BivariatePowerSeries {
    const prec = Math.min(this.prec, other.prec);
    const m = new Map(this.terms);
    for (const [key, c] of other.terms) {
      const cur = m.get(key);
      m.set(key, cur === undefined ? c : cur.add(c));
    }
    return new BivariatePowerSeries(this.k, prec, m);
  }

  sub(other: BivariatePowerSeries): BivariatePowerSeries {
    return this.add(other.neg());
  }

  neg(): BivariatePowerSeries {
    const m = new Map<string, RingElement>();
    for (const [key, c] of this.terms) {
      m.set(key, c.neg());
    }
    return new BivariatePowerSeries(this.k, this.prec, m);
  }

  scalarMul(c: RingElement): BivariatePowerSeries {
    const m = new Map<string, RingElement>();
    for (const [key, a] of this.terms) {
      m.set(key, a.mul(c));
    }
    return new BivariatePowerSeries(this.k, this.prec, m);
  }

  mul(other: BivariatePowerSeries): BivariatePowerSeries {
    const prec = Math.min(this.prec, other.prec);
    const m = new Map<string, RingElement>();
    for (const [k1, c1] of this.terms) {
      const [i1, j1] = k1.split(',').map(Number) as [number, number];
      for (const [k2, c2] of other.terms) {
        const [i2, j2] = k2.split(',').map(Number) as [number, number];
        const i = i1 + i2;
        const j = j1 + j2;
        if (i + j >= prec) continue;
        const key = `${i},${j}`;
        const cur = m.get(key);
        const prod = c1.mul(c2);
        m.set(key, cur === undefined ? prod : cur.add(prod));
      }
    }
    return new BivariatePowerSeries(this.k, prec, m);
  }

  pow(n: number): BivariatePowerSeries {
    let result = BivariatePowerSeries.one(this.k, this.prec);
    let base: BivariatePowerSeries = this;
    let e = n;
    while (e > 0) {
      if (e & 1) result = result.mul(base);
      base = base.mul(base);
      e >>= 1;
    }
    return result;
  }

  /** The smallest total degree occurring, or ``prec`` if the series is zero. */
  valuation(): number {
    let v = this.prec;
    for (const key of this.terms.keys()) {
      const [i, j] = key.split(',').map(Number) as [number, number];
      v = Math.min(v, i + j);
    }
    return v;
  }

  /** Multiplicative inverse; the constant term must be a unit. */
  inv(): BivariatePowerSeries {
    const c = this.coefficient(0, 0);
    if (c.isZero()) {
      throw new ZeroDivisionError('cannot invert a series with zero constant term');
    }
    const cinv = this.k.one().div(c);
    // self = c*(1 + z) with z of positive valuation, so 1/self = c^-1 * sum (-z)^m.
    const z = this.scalarMul(cinv).sub(BivariatePowerSeries.one(this.k, this.prec));
    let result = BivariatePowerSeries.one(this.k, this.prec);
    let zp = BivariatePowerSeries.one(this.k, this.prec);
    for (let m = 1; m < this.prec; m++) {
      zp = zp.mul(z.neg());
      if (zp.terms.size === 0) break;
      result = result.add(zp);
    }
    return result.scalarMul(cinv);
  }

  /** Truncate to O(t1,t2)^n. */
  add_bigoh(n: number): BivariatePowerSeries {
    return new BivariatePowerSeries(this.k, Math.min(this.prec, n), this.terms);
  }

  /** Substitute t1 -> A, t2 -> B. */
  subs(A: BivariatePowerSeries, B: BivariatePowerSeries): BivariatePowerSeries {
    const prec = Math.min(this.prec, A.prec, B.prec);
    let result = BivariatePowerSeries.zero(this.k, prec);
    const Apows: BivariatePowerSeries[] = [BivariatePowerSeries.one(this.k, prec)];
    const Bpows: BivariatePowerSeries[] = [BivariatePowerSeries.one(this.k, prec)];
    for (const [key, c] of this.terms) {
      const [i, j] = key.split(',').map(Number) as [number, number];
      while (Apows.length <= i) Apows.push(Apows[Apows.length - 1]!.mul(A));
      while (Bpows.length <= j) Bpows.push(Bpows[Bpows.length - 1]!.mul(B));
      result = result.add(Apows[i]!.mul(Bpows[j]!).scalarMul(c));
    }
    return result;
  }

  toString(): string {
    const keys = [...this.terms.keys()].sort((x, y) => {
      const [i1, j1] = x.split(',').map(Number) as [number, number];
      const [i2, j2] = y.split(',').map(Number) as [number, number];
      if (i1 + j1 !== i2 + j2) return i1 + j1 - (i2 + j2);
      return i2 - i1;
    });
    const parts: string[] = [];
    for (const key of keys) {
      const [i, j] = key.split(',').map(Number) as [number, number];
      const c = this.terms.get(key)!;
      let mon = '';
      if (i > 0) mon += i === 1 ? 't1' : `t1^${i}`;
      if (j > 0) mon += (mon ? '*' : '') + (j === 1 ? 't2' : `t2^${j}`);
      const cs = c.toString();
      if (mon === '') parts.push(cs);
      else if (cs === '1') parts.push(mon);
      else parts.push(`${cs}*${mon}`);
    }
    if (parts.length === 0) return `0 + O(t1, t2)^${this.prec}`;
    return `${parts.join(' + ')} + O(t1, t2)^${this.prec}`;
  }
}

/**
 * Substitute the bivariate series ``B`` (of positive valuation) into the
 * univariate power series ``f``.
 */
function substUnivariate(
  f: PowerSeriesElement<RingElement>,
  B: BivariatePowerSeries,
  k: CoefficientRing<RingElement>,
  prec: number
): BivariatePowerSeries {
  let result = BivariatePowerSeries.zero(k, prec);
  let Bp = BivariatePowerSeries.one(k, prec);
  const coeffs = f.list();
  for (let n = 0; n < prec; n++) {
    const c = coeffs[n];
    if (c !== undefined && !c.isZero()) {
      result = result.add(Bp.scalarMul(c));
    }
    Bp = Bp.mul(B);
    if (Bp.terms.size === 0) break;
  }
  return result;
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
  private _psRingCache: PowerSeriesRing<RingElement> | null = null;

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
    const R = new PowerSeriesRing<RingElement>(k, 't');

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
   * Return the formal group power series x(t).
   *
   * This is the formal expansion of x = t/w about the origin:
   *
   *   x(t) = t^(-2) - a_1*t^(-1) - a_2 - a_3*t - ...
   *
   * INPUT:
   * - prec: precision (default 20)
   *
   * OUTPUT: a Laurent series x(t)
   *
   * @see Reference: sage/schemes/elliptic_curves/formal_group.py:x
   */
  x(prec: number = 20): LaurentSeriesElement<RingElement> {
    prec = Math.max(prec, 0);
    const lsRing = new LaurentSeriesRing<RingElement>(this._E.base_ring, 't');
    const A = this._xLaurent(prec);
    return new LaurentSeriesElement<RingElement>(lsRing, A.s.add_bigoh(prec - A.v), A.v);
  }

  /**
   * Return the formal series y(t) = -1/w(t) in terms of the local parameter
   * t = -x/y at infinity:
   *
   *   y(t) = -t^(-3) + a_1*t^(-2) + a_2*t^(-1) + a_3 + ...
   *
   * @see Reference: sage/schemes/elliptic_curves/formal_group.py:y
   */
  y(prec: number = 20): LaurentSeriesElement<RingElement> {
    prec = Math.max(prec, 0);
    const lsRing = new LaurentSeriesRing<RingElement>(this._E.base_ring, 't');
    const A = this._yLaurent(prec);
    return new LaurentSeriesElement<RingElement>(lsRing, A.s.add_bigoh(prec - A.v), A.v);
  }

  /**
   * Return the power series f(t) = 1 + ... such that f(t) dt is the usual
   * invariant differential dx/(2y + a_1 x + a_3).
   *
   * @see Reference: sage/schemes/elliptic_curves/formal_group.py:differential
   */
  differential(prec: number = 20): PowerSeriesElement<RingElement> {
    prec = Math.max(prec, 0);

    if (this._cachedOmega !== null && prec <= this._cachedOmega.prec) {
      return this._cachedOmega.value.add_bigoh(prec);
    }

    const k = this._E.base_ring;
    const a = this._E.ainvs();

    // Sage: x = self.x(prec+1); y = self.y(prec+1)
    //       g = x.derivative() / (2*y + a[0]*x + a[2])
    const x = this._xLaurent(prec + 1);
    const y = this._yLaurent(prec + 1);
    const xprime = lauDeriv(x, k);

    const two = ringInt(k, 2);
    const denom = lauAdd(lauAdd(lauScalar(y, two), lauScalar(x, a[0])), {
      s: this._psRing().__call__([a[2]], Number.POSITIVE_INFINITY),
      v: 0,
    });

    const g = lauToPowerSeries(lauDiv(xprime, denom)).add_bigoh(prec);
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
   * @see Reference: sage/schemes/elliptic_curves/formal_group.py:log
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
   * @see Reference: sage/schemes/elliptic_curves/formal_group.py:inverse
   */
  inverse(prec: number = 20): PowerSeriesElement<RingElement> {
    prec = Math.max(prec, 0);

    if (this._cachedInverse !== null && prec <= this._cachedInverse.prec) {
      return this._cachedInverse.value.add_bigoh(prec);
    }

    const k = this._E.base_ring;
    const [a1, , a3] = this._E.ainvs();

    // Sage: inv = x / (y + a1*x + a3)   (page 114 of Silverman, AEC I)
    const x = this._xLaurent(prec);
    const y = this._yLaurent(prec);
    const denom = lauAdd(lauAdd(y, lauScalar(x, a1)), {
      s: this._psRing().__call__([a3], Number.POSITIVE_INFINITY),
      v: 0,
    });
    void k;
    const inv = lauToPowerSeries(lauDiv(x, denom)).add_bigoh(prec);
    this._cachedInverse = { prec, value: inv };
    return inv;
  }

  /**
   * Return the formal group law F(t1, t2) = t1 + t2 - a1*t1*t2 - ... to
   * precision O(t1, t2)^prec (page 115 of [Sil2009]).
   *
   * @see Reference: sage/schemes/elliptic_curves/formal_group.py:group_law
   */
  group_law(prec: number = 10): BivariatePowerSeries {
    prec = Math.max(prec, 0);
    if (prec <= 0) {
      throw new ValueError('The precision must be positive.');
    }

    const k = this._E.base_ring;
    const [a1, a2, a3, a4, a6] = this._E.ainvs();

    if (prec === 1) {
      return BivariatePowerSeries.zero(k, 1);
    }
    if (prec === 2) {
      const t1 = BivariatePowerSeries.monomial(k, 2, 1, 0, k.one());
      const t2 = BivariatePowerSeries.monomial(k, 2, 0, 1, k.one());
      const t1t2 = BivariatePowerSeries.monomial(k, 2, 1, 1, a1.neg());
      return t1.add(t2).add(t1t2);
    }

    if (this._cachedGroupLaw !== null && prec <= this._cachedGroupLaw.prec) {
      return this._cachedGroupLaw.value.add_bigoh(prec);
    }

    const T1 = BivariatePowerSeries.monomial(k, prec, 1, 0, k.one());
    const T2 = BivariatePowerSeries.monomial(k, prec, 0, 1, k.one());

    const w = this.w(prec + 1);

    // lam = sum_{n=3}^{prec} w[n] * sum_{m=0}^{n-1} t2^m * t1^(n-m-1)
    let lam = BivariatePowerSeries.zero(k, prec);
    for (let n = 3; n <= prec; n++) {
      const c = w.__getitem__(n);
      if (c.isZero()) continue;
      for (let m = 0; m < n; m++) {
        lam = lam.add(BivariatePowerSeries.monomial(k, prec, n - m - 1, m, c));
      }
    }

    // nu = w(t1) - lam*t1
    const nu = substUnivariate(w, T1, k, prec).sub(lam.mul(T1));

    const lam2 = lam.mul(lam);
    const lam3 = lam2.mul(lam);

    // Note: this formula differs from Silverman p.119; see Sage issue 9646.
    const numer = lam
      .scalarMul(a1)
      .add(lam2.scalarMul(a3))
      .add(nu.scalarMul(a2))
      .add(lam.mul(nu).scalarMul(a4.mul(ringInt(k, 2))))
      .add(lam2.mul(nu).scalarMul(a6.mul(ringInt(k, 3))));
    const denom = BivariatePowerSeries.one(k, prec)
      .add(lam.scalarMul(a2))
      .add(lam2.scalarMul(a4))
      .add(lam3.scalarMul(a6));

    const t3 = T1.neg().sub(T2).sub(numer.mul(denom.inv()));

    const inv = this.inverse(prec);
    const F = substUnivariate(inv, t3, k, prec);

    this._cachedGroupLaw = { prec, value: F };
    return F;
  }

  /**
   * Return the formal "multiplication by n" endomorphism [n](t) = n*t + ...
   * to precision O(t^prec) (Proposition 2.3 of [Sil2009]).
   *
   * @see Reference: sage/schemes/elliptic_curves/formal_group.py:mult_by_n
   */
  mult_by_n(n: bigint | number, prec: number = 10): PowerSeriesElement<RingElement> {
    let nVal = typeof n === 'number' ? BigInt(n) : n;
    const k = this._E.base_ring;
    const R = this._psRing();
    const t = R.gen();

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
    // Composition with the group law: [m+1](t) = F([m](t), t).
    const applyF = (
      A: PowerSeriesElement<RingElement>,
      B: PowerSeriesElement<RingElement>
    ): PowerSeriesElement<RingElement> => {
      let result = R.zero().add_bigoh(prec);
      const Apows: PowerSeriesElement<RingElement>[] = [R.one().add_bigoh(prec)];
      const Bpows: PowerSeriesElement<RingElement>[] = [R.one().add_bigoh(prec)];
      for (const [key, c] of F.terms) {
        const [i, j] = key.split(',').map(Number) as [number, number];
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

    // Double-and-add is faster than the naive method when n >= 4.
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
    void k;
    return result;
  }

  /**
   * Return the Weierstrass sigma function as a formal power series solution of
   * d^2 log(sigma)/dz^2 = -wp(z), expressed in the formal-group parameter t.
   *
   * @see Reference: sage/schemes/elliptic_curves/formal_group.py:sigma
   */
  sigma(prec: number = 10): PowerSeriesElement<RingElement> {
    const k = this._E.base_ring;
    const [a1, a2] = this._E.ainvs();
    const R = this._psRing();

    const fl = this.log(prec);
    const F = fl.reversion(prec);

    // wp = x(F) + (a1^2 + 4*a2)/12, as a Laurent series in z.
    const x = this._xLaurent(prec + 4);
    const twelve = ringInt(k, 12);
    if (twelve.isZero()) {
      throw new ZeroDivisionError('sigma requires 12 to be invertible in the base ring');
    }
    const shift = a1
      .mul(a1)
      .add(a2.mul(ringInt(k, 4)))
      .div(twelve);

    // x = t^-2 * u(t)  =>  x(F) = F^-2 * u(F)
    const uF = x.s.__call__(F);
    const Fv = F.valuation();
    const Funit = F._shiftRight(Fv);
    const FinvPow = Funit.inv().pow(-x.v); // (1/unit)^2
    const wp: Lau = { s: uF.mul(FinvPow), v: x.v * Fv };

    // g = 1/z^2 - wp  (a power series)
    const oneOverZ2: Lau = { s: R.one(), v: -2 };
    const g = lauToPowerSeries(
      lauSub(oneOverZ2, lauAdd(wp, { s: R.__call__([shift], Number.POSITIVE_INFINITY), v: 0 }))
    );

    const h = g.integral().integral();
    // sigma(z) = z * exp(h)
    const sigmaOfZ = h.exp().add_bigoh(prec)._shiftLeft(1);

    // sigma(t) = sigma(log(t))
    return sigmaOfZ.__call__(fl).add_bigoh(prec);
  }

  /** The power series ring R[[t]] over the base ring. */
  private _psRing(): PowerSeriesRing<RingElement> {
    if (this._psRingCache === null) {
      this._psRingCache = new PowerSeriesRing<RingElement>(this._E.base_ring, 't');
    }
    return this._psRingCache;
  }

  /** x(t) = t^-2 * u(t) with u = (w/t^3)^-1. */
  private _xLaurent(prec: number): Lau {
    return { s: this._unitInverse(prec + 3), v: -2 };
  }

  /** y(t) = -t^-3 * u(t). */
  private _yLaurent(prec: number): Lau {
    return { s: this._unitInverse(prec + 3).neg(), v: -3 };
  }

  /** The power series u(t) = (w(t)/t^3)^-1, which has constant term 1. */
  private _unitInverse(prec: number): PowerSeriesElement<RingElement> {
    const w = this.w(prec + 6);
    return w._shiftRight(3).add_bigoh(prec).inv();
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
