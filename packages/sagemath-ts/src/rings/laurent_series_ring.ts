/**
 * @module sage/rings/laurent_series_ring
 * @description Laurent series rings and their elements
 *
 * Port of: sage/rings/laurent_series_ring.py and
 *          sage/rings/laurent_series_ring_element.pyx
 *
 * Reference: reference/sage/src/sage/rings/laurent_series_ring.py
 * Reference: reference/sage/src/sage/rings/laurent_series_ring_element.pyx
 *
 * Laurent series are represented exactly as in SageMath: internally a power of
 * the variable times a power series part.  If a Laurent series `f` is written
 * `f = t^n * u` where `t` is the variable and `u` has nonzero constant term,
 * then `u` is {@link LaurentSeriesElement.valuation_zero_part} and `n` is
 * {@link LaurentSeriesElement.valuation}.
 *
 * @see Reference: sage/rings/laurent_series_ring_element.pyx:1-10
 */

import {
  ArithmeticError,
  IndexError,
  NotImplementedError,
  ValueError,
  ZeroDivisionError,
} from '../errors.js';
import {
  type CoefficientRing,
  PowerSeriesElement,
  PowerSeriesRing,
  type RingElement,
} from './power_series_ring.js';
import { ringCharacteristic } from './power_series_ring.js';
import { Rational } from './rational.js';

/**
 * Univariate Laurent series ring.
 *
 * SageMath's `LaurentSeriesRing` is a `UniqueRepresentation` parent built from
 * the corresponding power series ring; this port stores the power series ring
 * the same way (`self._power_series_ring`).
 *
 * @see Reference: sage/rings/laurent_series_ring.py:88 (LaurentSeriesRing)
 */
export class LaurentSeriesRing<T extends RingElement = RingElement> {
  private readonly _base_ring: CoefficientRing<T>;
  private readonly _name: string;
  private readonly _power_series_ring: PowerSeriesRing<T>;

  constructor(base_ring: CoefficientRing<T>, name: string = 'x', default_prec: number = 20) {
    this._base_ring = base_ring;
    this._name = name;
    this._power_series_ring = new PowerSeriesRing<T>(base_ring, name, default_prec);
  }

  /**
   * Return the base ring.
   * @see Reference: sage/rings/laurent_series_ring.py:233 (__init__)
   */
  base_ring(): CoefficientRing<T> {
    return this._base_ring;
  }

  /**
   * Return the variable name.
   * @see Reference: sage/rings/laurent_series_ring.py:397 (_repr_)
   */
  variable_name(): string {
    return this._name;
  }

  /**
   * If this is the Laurent series ring `R((t))`, return the power series ring
   * `R[[t]]`.
   * @see Reference: sage/rings/laurent_series_ring.py:944 (power_series_ring)
   */
  power_series_ring(): PowerSeriesRing<T> {
    return this._power_series_ring;
  }

  /**
   * Get the precision to which exact elements are truncated when necessary
   * (most frequently when inverting).
   * @see Reference: sage/rings/laurent_series_ring.py:832 (default_prec)
   */
  default_prec(): number {
    return this._power_series_ring.default_prec();
  }

  /**
   * A Laurent series ring is a field if and only if the base ring is a field.
   * @see Reference: sage/rings/laurent_series_ring.py:373 (is_field)
   */
  is_field(): boolean {
    return this._base_ring.is_field?.() ?? false;
  }

  /**
   * Laurent series rings are inexact.
   * @see Reference: sage/rings/laurent_series_ring.py:845 (is_exact)
   */
  is_exact(): boolean {
    return false;
  }

  /**
   * @see Reference: sage/rings/laurent_series_ring.py:799 (characteristic)
   */
  characteristic(): bigint {
    return ringCharacteristic(this._base_ring);
  }

  /**
   * Return the residue field of this Laurent series field if it is a complete
   * discrete valuation field (i.e. if the base ring is a field, in which case
   * it is also the residue field).
   * @see Reference: sage/rings/laurent_series_ring.py:809 (residue_field)
   */
  residue_field(): CoefficientRing<T> {
    if (!this.is_field()) {
      throw new TypeError('the base ring is not a field');
    }
    return this._base_ring;
  }

  /**
   * Return a uniformizer of this Laurent series field if it is a discrete
   * valuation field (i.e. if the base ring is actually a field).  Otherwise an
   * error is raised.
   * @see Reference: sage/rings/laurent_series_ring.py:880 (uniformizer)
   */
  uniformizer(): LaurentSeriesElement<T> {
    if (!this.is_field()) {
      throw new TypeError('the base ring is not a field');
    }
    return this.gen();
  }

  /**
   * @see Reference: sage/rings/laurent_series_ring.py:858 (gen)
   */
  gen(n: number = 0): LaurentSeriesElement<T> {
    if (n !== 0) {
      throw new RangeError(`generator ${n} not defined`);
    }
    return new LaurentSeriesElement<T>(this, [this._base_ring.zero(), this._base_ring.one()]);
  }

  /**
   * @see Reference: sage/rings/laurent_series_ring.py:870 (gens)
   */
  gens(): [LaurentSeriesElement<T>] {
    return [this.gen()];
  }

  /**
   * Laurent series rings are univariate.
   * @see Reference: sage/rings/laurent_series_ring.py:902 (ngens)
   */
  ngens(): number {
    return 1;
  }

  /**
   * Return the zero element.
   */
  zero(): LaurentSeriesElement<T> {
    return new LaurentSeriesElement<T>(this, this._power_series_ring.zero());
  }

  /**
   * Return the one element.
   * @see Reference: sage/rings/laurent_series_ring.py:233 (__init__, _one_element)
   */
  one(): LaurentSeriesElement<T> {
    return new LaurentSeriesElement<T>(this, this._power_series_ring.one());
  }

  /**
   * Return the Laurent series ring over `R` in the same variable.
   * @see Reference: sage/rings/laurent_series_ring.py:347 (change_ring)
   */
  change_ring<S extends RingElement>(R: CoefficientRing<S>): LaurentSeriesRing<S> {
    return new LaurentSeriesRing<S>(R, this._name, this.default_prec());
  }

  /**
   * Construct a Laurent series from `x`.
   *
   * INPUT:
   * - `x` -- object that can be converted into a Laurent series
   * - `n` -- (default: 0) multiply the result by `t^n`
   * - `prec` -- (default: `Infinity`) the precision of the series
   *
   * @see Reference: sage/rings/laurent_series_ring.py:438 (_element_constructor_)
   */
  __call__(
    x: unknown,
    n: number = 0,
    prec: number = Number.POSITIVE_INFINITY
  ): LaurentSeriesElement<T> {
    if (x instanceof LaurentSeriesElement && n === 0 && x.parent() === this) {
      // ok, since Laurent series are immutable (no need to make a copy)
      return (x as LaurentSeriesElement<T>).add_bigoh(prec);
    }
    return new LaurentSeriesElement<T>(this, x, n).add_bigoh(prec);
  }

  /**
   * @see Reference: sage/rings/laurent_series_ring.py:397 (_repr_)
   */
  toString(): string {
    return `Laurent Series Ring in ${this._name} over ${this._base_ring}`;
  }
}

/**
 * A Laurent series.
 *
 * We consider a Laurent series of the form `f = t^n * u` where `u` is a power
 * series with nonzero constant term.
 *
 * INPUT:
 * - `parent` -- a Laurent series ring
 * - `f` -- a power series (or something that can be coerced to one); note that
 *   `f` does *not* have to be a unit
 * - `n` -- (default: 0) integer
 *
 * @see Reference: sage/rings/laurent_series_ring_element.pyx:89 (LaurentSeries)
 */
export class LaurentSeriesElement<T extends RingElement = RingElement> {
  private readonly _parent: LaurentSeriesRing<T>;
  /** The unit part `u` (SageMath's `__u`). */
  private readonly _u: PowerSeriesElement<T>;
  /** The power of the variable `n` (SageMath's `__n`). */
  private readonly _n: number;

  /**
   * @see Reference: sage/rings/laurent_series_ring_element.pyx:104 (__init__)
   */
  constructor(parent: LaurentSeriesRing<T>, f: unknown, n: number = 0) {
    this._parent = parent;
    const R = parent.power_series_ring();

    let g: PowerSeriesElement<T>;
    if (f instanceof LaurentSeriesElement) {
      const lf = f as LaurentSeriesElement<T>;
      n += lf._n;
      // SageMath compares parents with ``is``; power series rings are unique
      // parents there, so the faithful analogue here is structural equality.
      g = R.is_identical_to(lf._u.parent()) ? lf._u : R.__call__(lf._u);
    } else if (f instanceof Map) {
      // Sanitize input to make sure all exponents are nonnegative, adjusting n
      // to match (SageMath's ``dict`` branch).
      const d = f as Map<number, unknown>;
      const keys = [...d.keys()];
      if (keys.length === 0) {
        g = R.zero();
      } else {
        let n1 = Math.min(...keys);
        if (n1 >= 0) {
          n1 = 0;
        } else {
          n += n1;
        }
        const maxKey = Math.max(...keys);
        const coeffs: T[] = [];
        for (let e = 0; e <= maxKey - n1; e++) {
          const c = d.get(e + n1);
          coeffs.push(c === undefined ? parent.base_ring().zero() : parent.base_ring().__call__(c));
        }
        g = R.__call__(coeffs);
      }
    } else if (!(f instanceof PowerSeriesElement)) {
      g = R.__call__(f);
    } else if (!R.is_identical_to((f as PowerSeriesElement<T>).parent())) {
      g = R.__call__(f);
    } else {
      g = f as PowerSeriesElement<T>;
    }

    // self is that t^n * u:
    if (g.is_zero()) {
      if (n === Number.POSITIVE_INFINITY) {
        this._n = 0;
        this._u = R.zero();
      } else {
        this._n = n;
        this._u = g;
      }
    } else {
      const val = g.valuation();
      if (val === Number.POSITIVE_INFINITY) {
        this._n = 0;
        this._u = g;
      } else if (val === 0) {
        this._n = n;
        this._u = g;
      } else {
        this._n = n + val;
        this._u = g._shiftRight(val);
      }
    }
  }

  /** Return the parent ring. */
  parent(): LaurentSeriesRing<T> {
    return this._parent;
  }

  /** Return the base ring of the parent. */
  base_ring(): CoefficientRing<T> {
    return this._parent.base_ring();
  }

  /**
   * Return `True` if this Laurent series is a unit in this ring.
   *
   * ALGORITHM: A Laurent series is a unit if and only if its "unit part" is a
   * unit.
   *
   * @see Reference: sage/rings/laurent_series_ring_element.pyx:187 (is_unit)
   */
  is_unit(): boolean {
    return this._u.is_unit();
  }

  /**
   * @see Reference: sage/rings/laurent_series_ring_element.pyx:222 (is_zero)
   */
  is_zero(): boolean {
    return this._u.is_zero();
  }

  /**
   * Return `True` if this element is a monomial, that is if it is `x^n` for
   * some integer `n`.
   * @see Reference: sage/rings/laurent_series_ring_element.pyx:237 (is_monomial)
   */
  is_monomial(): boolean {
    return this._u.is_monomial();
  }

  /**
   * `bool(self)`.
   * @see Reference: sage/rings/laurent_series_ring_element.pyx:257 (__bool__)
   */
  bool(): boolean {
    return !this._u.is_zero();
  }

  /**
   * The unit part `u` in `self = t^n u`.
   * @see Reference: sage/rings/laurent_series_ring_element.pyx:1293 (valuation_zero_part)
   */
  valuation_zero_part(): PowerSeriesElement<T> {
    return this._u;
  }

  /**
   * The valuation of `self`; `+Infinity` for an element indistinguishable
   * from zero.
   * @see Reference: sage/rings/laurent_series_ring_element.pyx:1310 (valuation)
   */
  valuation(): number {
    if (this.is_zero()) {
      return Number.POSITIVE_INFINITY;
    }
    return this._n;
  }

  /**
   * Return the name of the variable.
   * @see Reference: sage/rings/laurent_series_ring_element.pyx:1342 (variable)
   */
  variable(): string {
    return this._parent.variable_name();
  }

  /**
   * Return the `n` so that the Laurent series is of the form
   * `(stuff) + O(t^n)`.  It does not matter how many negative powers appear
   * in the expansion; in particular `prec` can be negative.
   * @see Reference: sage/rings/laurent_series_ring_element.pyx:1354 (prec)
   */
  prec(): number {
    return this._u.prec() + this._n;
  }

  /**
   * Return the absolute precision of this series (by definition the `r` in
   * `... + O(x^r)`).
   * @see Reference: sage/rings/laurent_series_ring_element.pyx:1373 (precision_absolute)
   */
  precision_absolute(): number {
    return this.prec();
  }

  /**
   * Return the relative precision of this series, that is the difference
   * between its absolute precision and its valuation.  By convention the
   * relative precision of `0` (or `O(x^r)`) is `0`.
   * @see Reference: sage/rings/laurent_series_ring_element.pyx:1390 (precision_relative)
   */
  precision_relative(): number {
    if (this.is_zero()) {
      return 0;
    }
    return this.prec() - this.valuation();
  }

  /**
   * Return the degree of a polynomial equivalent to this power series modulo
   * big oh of the precision.
   * @see Reference: sage/rings/laurent_series_ring_element.pyx:884 (degree)
   */
  degree(): number {
    return this._u.degree() + this._n;
  }

  /**
   * The coefficient of `t^i`.
   * @see Reference: sage/rings/laurent_series_ring_element.pyx:501 (__getitem__)
   */
  __getitem__(i: number): T {
    return this._u.__getitem__(i - this._n);
  }

  /**
   * The list of coefficients of the unit part, i.e. of `t^valuation`,
   * `t^(valuation+1)`, ...
   * @see Reference: sage/rings/laurent_series_ring_element.pyx:571 (list)
   */
  list(): T[] {
    return this._u.list();
  }

  /**
   * Return the nonzero coefficients of `self`.
   * @see Reference: sage/rings/laurent_series_ring_element.pyx:583 (coefficients)
   */
  coefficients(): T[] {
    return this.list().filter((c) => !c.isZero());
  }

  /**
   * Return the exponents appearing in `self` with nonzero coefficients.
   * @see Reference: sage/rings/laurent_series_ring_element.pyx:610 (exponents)
   */
  exponents(): number[] {
    const v = this.valuation();
    const out: number[] = [];
    const l = this.list();
    for (let i = 0; i < l.length; i++) {
      if (!l[i]!.isZero()) {
        out.push(i + v);
      }
    }
    return out;
  }

  /**
   * Return the residue of `self`, i.e. the coefficient of `t^-1`.
   * @see Reference: sage/rings/laurent_series_ring_element.pyx:594 (residue)
   */
  residue(): T {
    return this.__getitem__(-1);
  }

  /**
   * Return a congruent Laurent series with absolute precision at least
   * `absprec`.  If `absprec` is omitted, lift to an exact element.
   * @see Reference: sage/rings/laurent_series_ring_element.pyx:657 (lift_to_precision)
   */
  lift_to_precision(absprec?: number): LaurentSeriesElement<T> {
    if (absprec !== undefined && absprec <= this.precision_absolute()) {
      return this;
    }
    const P = this._parent;
    const exact = this.is_zero()
      ? P.__call__(0)
      : new LaurentSeriesElement<T>(P, P.power_series_ring().__call__(this.list()), this._n);
    if (absprec === undefined) {
      return exact;
    }
    return exact.add_bigoh(absprec);
  }

  /**
   * Add two Laurent series.
   *
   * ALGORITHM: Shift the unit parts to align them, then add.
   *
   * @see Reference: sage/rings/laurent_series_ring_element.pyx:735 (_add_)
   */
  add(right: LaurentSeriesElement<T>): LaurentSeriesElement<T> {
    // 1. Special case when one or the other is 0.
    if (!right.bool()) {
      return this.add_bigoh(right.prec());
    }
    if (!this.bool()) {
      return right.add_bigoh(this.prec());
    }

    // 2. Align the unit parts.
    let m: number;
    let f1: PowerSeriesElement<T>;
    let f2: PowerSeriesElement<T>;
    if (this._n < right._n) {
      m = this._n;
      f1 = this._u;
      f2 = right._u._shiftLeft(right._n - m);
    } else if (this._n > right._n) {
      m = right._n;
      f1 = this._u._shiftLeft(this._n - m);
      f2 = right._u;
    } else {
      m = this._n;
      f1 = this._u;
      f2 = right._u;
    }
    // 3. Add
    return new LaurentSeriesElement<T>(this._parent, f1.add(f2), m);
  }

  /**
   * Subtract two Laurent series.
   *
   * ALGORITHM: Shift the unit parts to align them, then subtract.
   *
   * @see Reference: sage/rings/laurent_series_ring_element.pyx:790 (_sub_)
   */
  sub(right: LaurentSeriesElement<T>): LaurentSeriesElement<T> {
    // 1. Special case when one or the other is 0.
    if (!right.bool()) {
      return this.add_bigoh(right.prec());
    }
    if (!this.bool()) {
      return right.neg().add_bigoh(this.prec());
    }

    // 2. Align the unit parts.
    let m: number;
    let f1: PowerSeriesElement<T>;
    let f2: PowerSeriesElement<T>;
    if (this._n < right._n) {
      m = this._n;
      f1 = this._u;
      f2 = right._u._shiftLeft(right._n - m);
    } else {
      m = right._n;
      f1 = this._u._shiftLeft(this._n - m);
      f2 = right._u;
    }
    // 3. Subtract
    return new LaurentSeriesElement<T>(this._parent, f1.sub(f2), m);
  }

  /**
   * @see Reference: sage/rings/laurent_series_ring_element.pyx:894 (__neg__)
   */
  neg(): LaurentSeriesElement<T> {
    return new LaurentSeriesElement<T>(this._parent, this._u.neg(), this._n);
  }

  /**
   * @see Reference: sage/rings/laurent_series_ring_element.pyx:906 (_mul_)
   */
  mul(right: LaurentSeriesElement<T>): LaurentSeriesElement<T> {
    return new LaurentSeriesElement<T>(this._parent, this._u.mul(right._u), this._n + right._n);
  }

  /**
   * Multiply by an element of the base ring.
   * @see Reference: sage/rings/laurent_series_ring_element.pyx:920 (_rmul_/_lmul_)
   */
  scalar_mul(c: T): LaurentSeriesElement<T> {
    return new LaurentSeriesElement<T>(this._parent, this._u.scalar_mul(c), this._n);
  }

  /**
   * @see Reference: sage/rings/laurent_series_ring_element.pyx:1152 (_div_)
   */
  div(right: LaurentSeriesElement<T>): LaurentSeriesElement<T> {
    if (right._u.is_zero()) {
      throw new ZeroDivisionError('');
    }
    return new LaurentSeriesElement<T>(this._parent, this._u.div(right._u), this._n - right._n);
  }

  /**
   * Return the inverse of `self`, i.e. `self^(-1)`.
   *
   * SageMath computes `~self`, which for a ring element is
   * `self.parent().one() / self`.
   *
   * @see Reference: sage/rings/laurent_series_ring_element.pyx:1858 (inverse)
   */
  inverse(): LaurentSeriesElement<T> {
    return this._parent.one().div(this);
  }

  /** Alias of {@link LaurentSeriesElement.inverse}. */
  inv(): LaurentSeriesElement<T> {
    return this.inverse();
  }

  /**
   * `self^r` for a rational exponent `r`.
   * @see Reference: sage/rings/laurent_series_ring_element.pyx:931 (__pow__)
   */
  pow(r: bigint | number | Rational): LaurentSeriesElement<T> {
    const right = r instanceof Rational ? r : new Rational(typeof r === 'number' ? BigInt(r) : r);

    if (right.denominator === 1n) {
      const e = right.numerator;
      return new LaurentSeriesElement<T>(this._parent, this._u.pow(e), this._n * Number(e));
    }

    if (this.is_zero()) {
      // SageMath: ``self._parent(0).O(self.prec()*right)``
      const p = this.prec();
      const scaled =
        p === Number.POSITIVE_INFINITY
          ? p
          : (p * Number(right.numerator)) / Number(right.denominator);
      return this._parent.__call__(0).O(scaled);
    }

    const d = Number(right.denominator);
    const n = right.numerator;

    const val = this.valuation();

    if (val % d !== 0) {
      throw new ValueError('power series valuation would be fractional');
    }

    const u = this.valuation_zero_part().nth_root(d);

    const s = new LaurentSeriesElement<T>(this._parent, u, val / d);

    return s.pow(n);
  }

  /**
   * Return this Laurent series multiplied by `t^k`.
   * @see Reference: sage/rings/laurent_series_ring_element.pyx:1006 (shift)
   */
  shift(k: number): LaurentSeriesElement<T> {
    return new LaurentSeriesElement<T>(this._parent, this._u, this._n + k);
  }

  /**
   * `self << k`, i.e. multiplication by `t^k`.
   * @see Reference: sage/rings/laurent_series_ring_element.pyx:1039 (__lshift__)
   */
  lshift(k: number): LaurentSeriesElement<T> {
    return new LaurentSeriesElement<T>(this._parent, this._u, this._n + k);
  }

  /**
   * `self >> k`, i.e. multiplication by `t^-k`.
   * @see Reference: sage/rings/laurent_series_ring_element.pyx:1042 (__rshift__)
   */
  rshift(k: number): LaurentSeriesElement<T> {
    return new LaurentSeriesElement<T>(this._parent, this._u, this._n - k);
  }

  /**
   * Return the truncated series at chosen precision `prec`.
   * @see Reference: sage/rings/laurent_series_ring_element.pyx:822 (add_bigoh)
   */
  add_bigoh(prec: number): LaurentSeriesElement<T> {
    if (prec === Number.POSITIVE_INFINITY || prec >= this.prec()) {
      return this;
    }
    const P = this._parent;
    if (!this.bool() || prec < this._n) {
      return new LaurentSeriesElement<T>(P, P.power_series_ring().__call__(0, 0), prec);
    }
    const u = this._u.add_bigoh(prec - this._n);
    return new LaurentSeriesElement<T>(P, u, this._n);
  }

  /**
   * Return the Laurent series of precision at most `prec` obtained by adding
   * `O(q^prec)`.
   * @see Reference: sage/rings/laurent_series_ring_element.pyx:865 (O)
   */
  O(prec: number): LaurentSeriesElement<T> {
    return this.add_bigoh(prec);
  }

  /**
   * Return the Laurent series of degree `< n` which is equivalent to `self`
   * modulo `x^n`.
   * @see Reference: sage/rings/laurent_series_ring_element.pyx:1078 (truncate)
   */
  truncate(n: number): LaurentSeriesElement<T> {
    if (n <= this._n) {
      return this._parent.zero();
    }
    return new LaurentSeriesElement<T>(this._parent, this._u._truncateSeries(n - this._n), this._n);
  }

  /**
   * Replace any terms of degree >= `n` by big oh.
   * @see Reference: sage/rings/laurent_series_ring_element.pyx:1096 (truncate_laurentseries)
   */
  truncate_laurentseries(n: number): LaurentSeriesElement<T> {
    if (n <= this._n) {
      return this._parent.zero();
    }
    return new LaurentSeriesElement<T>(
      this._parent,
      this._u.truncate_powerseries(n - this._n),
      this._n
    );
  }

  /**
   * Return the Laurent series equivalent to `self` except without any terms of
   * degree `< n`.  Equivalent to `self - self.truncate(n)`.
   * @see Reference: sage/rings/laurent_series_ring_element.pyx:1116 (truncate_neg)
   */
  truncate_neg(n: number): LaurentSeriesElement<T> {
    return new LaurentSeriesElement<T>(this._parent, this._u._shiftRight(n - this._n), n);
  }

  /**
   * Return the minimum precision of `self` and `other`.
   * @see Reference: sage/rings/laurent_series_ring_element.pyx:1180 (common_prec)
   */
  common_prec(other: LaurentSeriesElement<T>): number {
    return Math.min(this.prec(), other.prec());
  }

  /**
   * Return the minimum valuation of `self` and `other`.
   * @see Reference: sage/rings/laurent_series_ring_element.pyx:1231 (common_valuation)
   */
  common_valuation(other: LaurentSeriesElement<T>): number {
    return Math.min(this.valuation(), other.valuation());
  }

  /**
   * Comparison of `self` and `right`.
   *
   * Two approximate Laurent series are equal if they agree for all coefficients
   * up to the *minimum* of the precisions of each.
   *
   * SageMath's `_richcmp_` also orders elements (dictionary order from lowest
   * to highest degree); ordering is not ported because the `RingElement`
   * interface of this port has no order on coefficients.
   *
   * @see Reference: sage/rings/laurent_series_ring_element.pyx:1284 (_richcmp_)
   */
  eq(right: LaurentSeriesElement<T>): boolean {
    const val = this.common_valuation(right);
    if (val === Number.POSITIVE_INFINITY) {
      return true; // Both arguments are zero
    }

    let deg = Math.max(this.degree(), right.degree());
    const prec = this.common_prec(right);
    if (deg >= prec) {
      deg = prec - 1;
    }

    for (let i = val; i <= deg; i++) {
      const li = this.__getitem__(i);
      const ri = right.__getitem__(i);
      if (!li.eq(ri)) {
        return false;
      }
    }
    return true;
  }

  /**
   * Return the `n`-th Verschiebung of `self`: if `f = sum a_m x^m` then this
   * returns `sum a_m x^{mn}`.
   * @see Reference: sage/rings/laurent_series_ring_element.pyx:366 (verschiebung)
   */
  verschiebung(n: number): LaurentSeriesElement<T> {
    if (n === 0) {
      throw new ValueError('n must be nonzero');
    }

    if (n < 0) {
      if (this.prec() !== Number.POSITIVE_INFINITY) {
        throw new ValueError('For finite precision only positive arguments allowed');
      }
      const exponents = this.exponents().map((e) => e * n);
      const u = Math.min(...exponents);
      const shifted = exponents.map((e) => e - u);
      const coefficients = this.coefficients();
      const zero = this.base_ring().zero();
      const w: T[] = new Array(Math.max(...shifted) + 1).fill(zero);
      for (let i = 0; i < shifted.length; i++) {
        w[shifted[i]!] = coefficients[i]!;
      }
      return new LaurentSeriesElement<T>(this._parent, w, u);
    }
    return new LaurentSeriesElement<T>(this._parent, this._u.V(n), this._n * n);
  }

  /** Alias of {@link LaurentSeriesElement.verschiebung}. */
  V(n: number): LaurentSeriesElement<T> {
    return this.verschiebung(n);
  }

  /**
   * Return the reverse of `f`, i.e. the series `g` such that `g(f(x)) = x`.
   *
   * This is only possible if the valuation of `self` is exactly 1.
   *
   * @see Reference: sage/rings/laurent_series_ring_element.pyx:1538 (reverse)
   */
  reverse(precision?: number): LaurentSeriesElement<T> {
    const val = this.valuation();
    if (val !== 1) {
      throw new ValueError('Series must have valuation one for reversion.');
    }
    let u = this.valuation_zero_part();
    u = u.parent().gen(0).mul(u);

    const rev = u.reversion(precision);

    return this._parent.__call__(rev);
  }

  /**
   * Return whether this Laurent series is a square.
   *
   * @param root -- if `true`, return a pair `[true, sqrt]` if this element is a
   *   square and `[false, null]` otherwise.
   *
   * @see Reference: sage/rings/laurent_series_ring_element.pyx:1637 (is_square)
   */
  is_square(root?: false): boolean;
  is_square(root: true): [boolean, LaurentSeriesElement<T> | null];
  is_square(root: boolean = false): boolean | [boolean, LaurentSeriesElement<T> | null] {
    // Case 1: Handle Zero
    if (this.is_zero()) {
      if (root) {
        return [true, this];
      }
      return true;
    }

    // Case 2: Valuation must be even
    const v = this.valuation();
    if (((v % 2) + 2) % 2 !== 0) {
      if (root) {
        return [false, null];
      }
      return false;
    }

    // Case 3: The unit part must be a square
    const unit_part = this.rshift(v).power_series();

    // We use a try-except block to handle inconsistent API in base rings
    let is_sq: boolean;
    try {
      is_sq = unit_part.is_square();
    } catch (e) {
      if (
        e instanceof TypeError ||
        e instanceof ValueError ||
        e instanceof ArithmeticError ||
        e instanceof NotImplementedError
      ) {
        if (root) {
          return [false, null];
        }
        return false;
      }
      throw e;
    }

    if (!root) {
      return is_sq;
    }

    if (is_sq) {
      let sqrt_unit: PowerSeriesElement<T>;
      try {
        sqrt_unit = unit_part.sqrt();
      } catch (e) {
        if (e instanceof ValueError || e instanceof ArithmeticError) {
          return [false, null];
        }
        throw e;
      }
      // Reconstruct: t^(v/2) * sqrt(unit)
      return [true, this._parent.__call__(sqrt_unit).lshift(v / 2)];
    }
    return [false, null];
  }

  /**
   * The formal derivative of this Laurent series with respect to the variable.
   * @see Reference: sage/rings/laurent_series_ring_element.pyx:1790 (_derivative)
   */
  derivative(): LaurentSeriesElement<T> {
    const R = this._parent.power_series_ring();
    const k = this._parent.base_ring();
    const n = this._n;
    if (this.is_zero()) {
      const p = this._u.prec();
      return new LaurentSeriesElement<T>(
        this._parent,
        R.zero(),
        p === Number.POSITIVE_INFINITY ? p : p - 1
      );
    }
    const a = this._u.list();
    const v: T[] = [];
    for (let m = 0; m < a.length; m++) {
      v.push(k.__call__(BigInt(n + m)).mul(a[m]!) as T);
    }
    const u = R.__call__(v, this._u.prec());
    return new LaurentSeriesElement<T>(this._parent, u, n - 1);
  }

  /**
   * The formal integral of this Laurent series with 0 constant term.
   * @see Reference: sage/rings/laurent_series_ring_element.pyx:1857 (integral)
   */
  integral(): LaurentSeriesElement<T> {
    const R = this._parent.power_series_ring();
    const k = this._parent.base_ring();
    const n = this._n;
    const a = this._u.list();
    if (!this.__getitem__(-1).isZero()) {
      throw new ArithmeticError(
        'The integral of is not a Laurent series, since t^-1 has nonzero coefficient.'
      );
    }

    const v: T[] = [];
    let u: PowerSeriesElement<T>;
    try {
      if (n < 0) {
        for (let i = 0; i < Math.min(-1 - n, a.length); i++) {
          v.push(a[i]!.div(k.__call__(BigInt(n + i + 1))) as T);
        }
        v.push(k.zero());
      }
      for (let i = Math.max(-n, 0); i < a.length; i++) {
        v.push(a[i]!.div(k.__call__(BigInt(n + i + 1))) as T);
      }
      u = R.__call__(v, this._u.prec());
    } catch {
      throw new ArithmeticError('Coefficients of integral cannot be coerced into the base ring');
    }
    return new LaurentSeriesElement<T>(this._parent, u, n + 1);
  }

  /**
   * Return the `n`-th root of this Laurent power series.
   * @see Reference: sage/rings/laurent_series_ring_element.pyx:1913 (nth_root)
   */
  nth_root(n: number, prec?: number): LaurentSeriesElement<T> {
    let p: number;
    if (prec === undefined) {
      p = this.prec();
      if (p === Number.POSITIVE_INFINITY) {
        p = this.parent().default_prec();
      }
    } else {
      p = Math.min(this.prec(), prec);
    }

    if (n <= 0) {
      throw new ValueError('n must be positive');
    }

    const i = this.valuation();
    if (i % n !== 0) {
      throw new ValueError('valuation must be divisible by n');
    }

    const q = this._u.nth_root(n, p);
    return new LaurentSeriesElement<T>(this._parent, q.add_bigoh(p), i / n);
  }

  /**
   * Convert this Laurent series to a power series.
   *
   * An error is raised if the Laurent series has a term (or an error term
   * `O(x^k)`) whose exponent is negative.
   *
   * @see Reference: sage/rings/laurent_series_ring_element.pyx:1954 (power_series)
   */
  power_series(): PowerSeriesElement<T> {
    if (this._n < 0) {
      if (this._u.is_zero() && this._u.prec() >= -this._n) {
        return this._u._shiftRight(-this._n);
      }
      throw new TypeError('self is not a power series');
    }
    return this._u._shiftLeft(this._n);
  }

  /**
   * Compute the value of this Laurent series at `x`.
   *
   * It is only possible to substitute elements of positive valuation.
   *
   * @see Reference: sage/rings/laurent_series_ring_element.pyx:2020 (__call__)
   */
  __call__(x: LaurentSeriesElement<T>): LaurentSeriesElement<T> {
    // Upstream (`laurent_series_ring_element.pyx:1955`) is
    // `self.__u(*x) * (x[0]**self.__n)`: the argument is handed STRAIGHT to the
    // power series `__call__`, which permits any argument when `self` has
    // infinite precision.  Converting via `x.power_series()` first (as this used
    // to) rejected arguments of negative valuation even when the unit part is an
    // exact polynomial, so `(x^-1 + 1)(x^-2)` failed instead of giving `1 + x^2`.
    const P = x.parent();
    let composed: LaurentSeriesElement<T>;
    if (this._u.prec() === Number.POSITIVE_INFINITY) {
      // Exact polynomial unit part: Horner directly in the Laurent ring.
      const coeffs = this._u.list();
      composed = P.zero();
      for (let i = coeffs.length - 1; i >= 0; i--) {
        composed = composed.mul(x).add(P.__call__(coeffs[i]!));
      }
    } else {
      const val = x.valuation();
      if (val <= 0) {
        throw new ValueError('Can only substitute elements of positive valuation');
      }
      composed = P.__call__(this._u.__call__(x.power_series()));
    }
    return composed.mul(x.pow(this._n));
  }

  /**
   * @see Reference: sage/rings/laurent_series_ring_element.pyx:325 (_repr_)
   */
  toString(): string {
    if (this.is_zero()) {
      if (this.prec() === Number.POSITIVE_INFINITY) {
        return '0';
      }
      return `O(${this._parent.variable_name()}^${this.prec()})`;
    }
    let s = ' ';
    const v = this._u.list();
    const valuation = this._n;
    const m = v.length;
    const X = this._parent.variable_name();
    const atomic_repr = this._parent.base_ring().element_is_atomic?.() ?? true;
    let first = true;
    for (let n = 0; n < m; n++) {
      let x = v[n]!.toString();
      const e = n + valuation;
      if (x !== '0') {
        if (!first) {
          s += ' + ';
        }
        if (!atomic_repr && (x.slice(1).includes('+') || x.slice(1).includes('-'))) {
          x = `(${x})`;
        }
        let vr: string;
        if (e === 1) {
          vr = `*${X}`;
        } else if (e === 0) {
          vr = '';
        } else {
          vr = `*${X}^${e}`;
        }
        s += `${x}${vr}`;
        first = false;
      }
    }
    s = s.replaceAll(' + -', ' - ');
    s = s.replaceAll(' 1*', ' ');
    s = s.replaceAll(' -1*', ' -');
    const pr = this.prec();
    let bigoh: string;
    if (pr === 0) {
      bigoh = 'O(1)';
    } else if (pr === 1) {
      bigoh = `O(${X})`;
    } else {
      bigoh = `O(${X}^${pr})`;
    }
    if (pr !== Number.POSITIVE_INFINITY) {
      if (s === ' ') {
        return bigoh;
      }
      s += ` + ${bigoh}`;
    }
    return s.slice(1);
  }
}
