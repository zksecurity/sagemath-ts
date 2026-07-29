/**
 * @module sage/rings/power_series_ring
 * @description Power series rings
 *
 * Port of: sage/rings/power_series_ring.py
 * Reference: reference/sage/src/sage/rings/power_series_ring.py
 */

import { ArithmeticError, NotImplementedError, ValueError, ZeroDivisionError } from '../errors.js';
import { LaurentSeriesRing } from './laurent_series_ring.js';

// Laurent series live in their own module (mirroring
// ``sage/rings/laurent_series_ring.py``); they are re-exported here because
// they used to be defined in this file.
export { LaurentSeriesElement, LaurentSeriesRing } from './laurent_series_ring.js';

/**
 * Return a sequence of integers `1 = a_1 <= a_2 <= ... <= a_n = N` such that
 * `a_{i+1} <= 2 a_i`, suitable for driving a Newton iteration.
 *
 * @see Reference: sage/misc/misc.py:newton_method_sizes
 */
function newton_method_sizes(N: number): number[] {
  if (N < 1) {
    throw new ValueError(`N (=${N}) must be a positive integer`);
  }
  const output: number[] = [];
  let n = Math.floor(N);
  while (n > 1) {
    output.push(n);
    n = (n + 1) >> 1;
  }
  output.push(1);
  output.reverse();
  return output;
}

/**
 * Interface for ring elements that can be used as coefficients.
 */
export interface RingElement {
  add(other: RingElement): RingElement;
  sub(other: RingElement): RingElement;
  mul(other: RingElement): RingElement;
  div(other: RingElement): RingElement;
  neg(): RingElement;
  eq(other: RingElement | number | bigint): boolean;
  isZero(): boolean;
  isOne?(): boolean;
  isUnit?(): boolean;
  inv?(): RingElement;
  toString(): string;
}

/**
 * Interface for a coefficient ring.
 */
export interface CoefficientRing<T extends RingElement = RingElement> {
  zero(): T;
  one(): T;
  __call__(x: unknown): T;
  is_field?(): boolean;
  characteristic?(): bigint;
  /**
   * SageMath's ``R._repr_option('element_is_atomic')``: whether elements print
   * without needing parentheses inside a sum.  Defaults to ``true`` (as for
   * ZZ, QQ and finite fields).
   *
   * @see Reference: sage/rings/laurent_series_ring_element.pyx:333 (_repr_)
   */
  element_is_atomic?(): boolean;
  /** How the ring prints itself inside a ring's `_repr_`. */
  toString?(): string;
}

/**
 * A ring of power series over a base ring.
 *
 * Power series are formal infinite sums sum_{n>=0} a_n * x^n.
 *
 * @see Reference: sage/rings/power_series_ring.py:PowerSeriesRing_generic
 */
export class PowerSeriesRing<T extends RingElement = RingElement> {
  private readonly _base_ring: CoefficientRing<T>;
  private readonly _name: string;
  private readonly _default_prec: number;
  private readonly _generator: PowerSeriesElement<T>;

  constructor(base_ring: CoefficientRing<T>, name: string = 'x', default_prec: number = 20) {
    if (default_prec < 0) {
      throw new ValueError(`default_prec (= ${default_prec}) must be nonnegative`);
    }
    this._base_ring = base_ring;
    this._name = name;
    this._default_prec = default_prec;
    // Create the generator x with infinite precision
    this._generator = new PowerSeriesElement<T>(
      this,
      [base_ring.zero(), base_ring.one()],
      Number.POSITIVE_INFINITY
    );
  }

  /**
   * Return the base ring.
   * @see Reference: sage/rings/power_series_ring.py:base_ring
   */
  base_ring(): CoefficientRing<T> {
    return this._base_ring;
  }

  /**
   * Return the variable name.
   * @see Reference: sage/rings/power_series_ring.py:variable_name
   */
  variable_name(): string {
    return this._name;
  }

  /**
   * Return the default precision.
   * @see Reference: sage/rings/power_series_ring.py:default_prec
   */
  default_prec(): number {
    return this._default_prec;
  }

  /**
   * Return the generator (the variable x).
   * @see Reference: sage/rings/power_series_ring.py:gen
   */
  gen(n: number = 0): PowerSeriesElement<T> {
    if (n !== 0) {
      throw new Error('generator n>0 not defined');
    }
    return this._generator;
  }

  /**
   * Return the number of generators (always 1).
   * @see Reference: sage/rings/power_series_ring.py:ngens
   */
  ngens(): number {
    return 1;
  }

  /**
   * Whether `other` is the same parent as `self`.
   *
   * SageMath's power series rings are `UniqueRepresentation` parents, so it
   * compares them with `is`; this port has no parent cache, so the faithful
   * analogue is equality of the defining data (base ring, variable name and
   * default precision -- exactly SageMath's `UniqueRepresentation` key).
   *
   * @see Reference: sage/rings/power_series_ring.py:PowerSeriesRing (unique parents)
   */
  is_identical_to(other: PowerSeriesRing<T>): boolean {
    return (
      this === other ||
      (this._base_ring === other._base_ring &&
        this._name === other._name &&
        this._default_prec === other._default_prec)
    );
  }

  /**
   * Return the characteristic.
   * @see Reference: sage/rings/power_series_ring.py:characteristic
   */
  characteristic(): bigint {
    if (this._base_ring.characteristic) {
      return this._base_ring.characteristic();
    }
    return 0n;
  }

  /**
   * Coerce an element to this ring.
   * @see Reference: sage/rings/power_series_ring.py:__call__
   */
  __call__(f: unknown, prec?: number): PowerSeriesElement<T> {
    if (prec !== undefined && prec < 0) {
      throw new ValueError(`prec (= ${prec}) must be nonnegative`);
    }

    const actualPrec = prec ?? Number.POSITIVE_INFINITY;

    if (f instanceof PowerSeriesElement) {
      // If already a power series from this ring, possibly truncate
      if (f.parent() === this) {
        if (actualPrec >= f.prec()) {
          return f as PowerSeriesElement<T>;
        }
        return f.truncate(actualPrec) as unknown as PowerSeriesElement<T>;
      }
      // From another power series ring - convert coefficients
      const coeffs: T[] = [];
      const fList = f.list() as RingElement[];
      for (const c of fList) {
        coeffs.push(this._base_ring.__call__(c));
      }
      return new PowerSeriesElement<T>(this, coeffs, Math.min(actualPrec, f.prec()));
    }

    if (Array.isArray(f)) {
      // From a list of coefficients
      const coeffs: T[] = f.map((c) => this._base_ring.__call__(c));
      return new PowerSeriesElement<T>(this, coeffs, actualPrec);
    }

    // From a scalar
    const coeff = this._base_ring.__call__(f);
    return new PowerSeriesElement<T>(this, [coeff], actualPrec);
  }

  /**
   * Return zero.
   * @see Reference: sage/rings/power_series_ring.py:zero
   */
  zero(): PowerSeriesElement<T> {
    return new PowerSeriesElement<T>(this, [], Number.POSITIVE_INFINITY);
  }

  /**
   * Return one.
   * @see Reference: sage/rings/power_series_ring.py:one
   */
  one(): PowerSeriesElement<T> {
    return new PowerSeriesElement<T>(this, [this._base_ring.one()], Number.POSITIVE_INFINITY);
  }

  /**
   * Return a random element.
   * @see Reference: sage/rings/power_series_ring.py:random_element
   */
  random_element(degree?: number): PowerSeriesElement<T> {
    // Reference: sage/rings/power_series_ring.py:random_element
    // Generate a random power series with the specified precision

    const prec = degree ?? this._default_prec;
    const coeffs: T[] = [];

    // Check if base ring has random_element method
    if (
      'random_element' in this._base_ring &&
      typeof (this._base_ring as unknown as { random_element: () => T }).random_element ===
        'function'
    ) {
      for (let i = 0; i < prec; i++) {
        coeffs.push((this._base_ring as unknown as { random_element: () => T }).random_element());
      }
    } else {
      // Fallback: use 0 and 1 randomly
      for (let i = 0; i < prec; i++) {
        coeffs.push(Math.random() < 0.5 ? this._base_ring.zero() : this._base_ring.one());
      }
    }

    return new PowerSeriesElement<T>(this, coeffs, prec);
  }

  /**
   * Return the corresponding Laurent series ring.
   * @see Reference: sage/rings/power_series_ring.py:laurent_series_ring
   */
  laurent_series_ring(): LaurentSeriesRing<T> {
    return new LaurentSeriesRing<T>(this._base_ring, this._name, this._default_prec);
  }

  toString(): string {
    return `Power Series Ring in ${this._name} over ${this._base_ring}`;
  }
}

/**
 * An element of a power series ring.
 * @see Reference: sage/rings/power_series_ring_element.pyx:PowerSeries
 */
export class PowerSeriesElement<T extends RingElement = RingElement> {
  private readonly _parent: PowerSeriesRing<T>;
  private readonly _coefficients: T[];
  private readonly _prec: number;

  constructor(parent: PowerSeriesRing<T>, coefficients: T[], prec: number) {
    this._parent = parent;
    // SageMath's ``PowerSeries_poly.__init__`` stores ``f.truncate(prec)``:
    // coefficients of degree >= prec are not part of the element.
    // Reference: sage/rings/power_series_poly.pyx:PowerSeries_poly.__init__
    const known =
      prec === Number.POSITIVE_INFINITY ? coefficients : coefficients.slice(0, Math.max(0, prec));
    // Strip trailing zeros
    let lastNonZero = -1;
    for (let i = known.length - 1; i >= 0; i--) {
      if (!known[i]!.isZero()) {
        lastNonZero = i;
        break;
      }
    }
    this._coefficients = lastNonZero < 0 ? [] : known.slice(0, lastNonZero + 1);
    this._prec = prec;
  }

  /**
   * Return the parent ring.
   * @see Reference: sage/rings/power_series_ring_element.pyx:parent
   */
  parent(): PowerSeriesRing<T> {
    return this._parent;
  }

  /**
   * Return the precision.
   * @see Reference: sage/rings/power_series_ring_element.pyx:prec
   */
  prec(): number {
    return this._prec;
  }

  /**
   * Return the coefficients as a list.
   * @see Reference: sage/rings/power_series_ring_element.pyx:list
   */
  list(): T[] {
    return [...this._coefficients];
  }

  /**
   * Return the i-th coefficient.
   * @see Reference: sage/rings/power_series_ring_element.pyx:__getitem__
   */
  __getitem__(n: number): T {
    if (n < 0) {
      return this._parent.base_ring().zero();
    }
    if (n >= this._coefficients.length) {
      if (this._prec > n) {
        return this._parent.base_ring().zero();
      } else {
        throw new Error('coefficient not known');
      }
    }
    return this._coefficients[n]!;
  }

  /**
   * Return the valuation (smallest n with non-zero coefficient).
   * @see Reference: sage/rings/power_series_ring_element.pyx:valuation
   */
  valuation(): number {
    // SageMath: ``if self.__f == 0: return self._prec`` -- only the exact zero
    // (whose precision is +Infinity) has infinite valuation; O(x^r) has
    // valuation r.
    // Reference: sage/rings/power_series_poly.pyx:valuation
    if (this._coefficients.length === 0) {
      return this._prec;
    }
    for (let i = 0; i < this._coefficients.length; i++) {
      if (!this._coefficients[i]!.isZero()) {
        return i;
      }
    }
    return this._prec;
  }

  /**
   * Return the absolute precision of this series (by definition the `r` in
   * `... + O(x^r)`).
   * @see Reference: sage/rings/power_series_ring_element.pyx:precision_absolute
   */
  precision_absolute(): number {
    return this._prec;
  }

  /**
   * Return the relative precision, i.e. the difference between the absolute
   * precision and the valuation. By convention this is 0 for `O(x^r)`.
   * @see Reference: sage/rings/power_series_ring_element.pyx:precision_relative
   */
  precision_relative(): number {
    if (this.is_zero()) {
      return 0;
    }
    return this._prec - this.valuation();
  }

  /**
   * Return the degree (largest n with non-zero coefficient before prec).
   * @see Reference: sage/rings/power_series_ring_element.pyx:degree
   */
  degree(): number {
    if (this._coefficients.length === 0) {
      return -1;
    }
    return this._coefficients.length - 1;
  }

  /**
   * Check if this is zero.
   * @see Reference: sage/rings/power_series_ring_element.pyx:is_zero
   */
  is_zero(): boolean {
    return this._coefficients.length === 0;
  }

  /**
   * Check if this is one.
   * @see Reference: sage/rings/power_series_ring_element.pyx:is_one
   */
  is_one(): boolean {
    if (this._coefficients.length !== 1) {
      return false;
    }
    const c = this._coefficients[0]!;
    if (c.isOne) {
      return c.isOne();
    }
    return c.eq(1);
  }

  /**
   * Check if this is a unit.
   * A power series is invertible iff its constant term is invertible.
   * @see Reference: sage/rings/power_series_ring_element.pyx:is_unit
   */
  is_unit(): boolean {
    if (this._prec === 0) {
      return false;
    }
    if (this._coefficients.length === 0) {
      return false;
    }
    const c0 = this._coefficients[0]!;
    if (c0.isUnit) {
      return c0.isUnit();
    }
    // If the base ring is a field, any non-zero element is a unit
    if (this._parent.base_ring().is_field?.()) {
      return !c0.isZero();
    }
    // For integers, only +/-1 are units
    return c0.eq(1) || c0.eq(-1);
  }

  /**
   * Return the truncation to precision n.
   * Returns a polynomial (i.e., power series with infinite precision truncated to degree < n).
   * @see Reference: sage/rings/power_series_ring_element.pyx:truncate
   */
  truncate(n?: number): PowerSeriesElement<T> {
    const prec = n ?? this._prec;
    if (prec === Number.POSITIVE_INFINITY) {
      return this;
    }
    const coeffs = this._coefficients.slice(0, prec);
    return new PowerSeriesElement<T>(this._parent, coeffs, Number.POSITIVE_INFINITY);
  }

  /**
   * Given input `prec` = `n`, return the power series of degree `< n` which is
   * equivalent to `self` modulo `x^n` (keeping the big-oh term).
   * @see Reference: sage/rings/power_series_poly.pyx:765 (truncate_powerseries)
   */
  truncate_powerseries(prec: number): PowerSeriesElement<T> {
    return new PowerSeriesElement<T>(
      this._parent,
      this._coefficients.slice(0, prec === Number.POSITIVE_INFINITY ? undefined : prec),
      Math.min(this._prec, prec)
    );
  }

  /**
   * Factor `self` as `q^n (a_0 + a_1 q + ...)` with `a_0` nonzero and return
   * `a_0 + a_1 q + ...`.
   * @see Reference: sage/rings/power_series_ring_element.pyx:1012 (valuation_zero_part)
   */
  valuation_zero_part(): PowerSeriesElement<T> {
    if (this.is_zero()) {
      throw new ValueError('power series has no valuation 0 part');
    }
    const n = this.valuation();
    if (n === 0) {
      return this;
    }
    return this._shiftRight(n);
  }

  /**
   * Return `true` if this element is a monomial, that is `c*x^n`.
   * @see Reference: sage/rings/power_series_ring_element.pyx:1257 (is_monomial)
   */
  is_monomial(): boolean {
    // SageMath: ``self.polynomial().is_monomial()`` -- a single term with
    // coefficient one.
    if (this._coefficients.length === 0) {
      return false;
    }
    const d = this._coefficients.length - 1;
    for (let i = 0; i < d; i++) {
      if (!this._coefficients[i]!.isZero()) {
        return false;
      }
    }
    const c = this._coefficients[d]!;
    return c.isOne ? c.isOne() : c.eq(1);
  }

  /**
   * If `f = sum a_m x^m` then return `sum a_m x^{nm}`.
   * @see Reference: sage/rings/power_series_ring_element.pyx:2695 (V)
   */
  V(n: number): PowerSeriesElement<T> {
    const v = this.list();
    const zero = this._parent.base_ring().zero();
    const w: T[] = [];
    let m = 0;
    for (let i = 0; i < v.length * n; i++) {
      if (i % n !== 0) {
        w.push(zero);
      } else {
        w.push(v[m]!);
        m += 1;
      }
    }
    return new PowerSeriesElement<T>(
      this._parent,
      w,
      this._prec === Number.POSITIVE_INFINITY ? Number.POSITIVE_INFINITY : this._prec * n
    );
  }

  /**
   * Return this power series multiplied by `x^n`.
   * @see Reference: sage/rings/power_series_ring_element.pyx:1180 (shift)
   */
  shift(n: number): PowerSeriesElement<T> {
    return n >= 0 ? this._shiftLeft(n) : this._shiftRight(-n);
  }

  /**
   * Return `true` if this power series has a square root in this ring.
   *
   * ALGORITHM: if the base ring is a field, this is true whenever the power
   * series has even valuation and the leading coefficient is a perfect square.
   *
   * @see Reference: sage/rings/power_series_ring_element.pyx:1564 (is_square)
   */
  is_square(): boolean {
    const val = this.valuation();
    if (val !== Number.POSITIVE_INFINITY && ((val % 2) + 2) % 2 === 1) {
      return false;
    }
    const lead = this.__getitem__(val === Number.POSITIVE_INFINITY ? 0 : val);
    const leadAny = lead as unknown as { is_square?: () => boolean };
    if (typeof leadAny.is_square !== 'function') {
      // SageMath calls ``self[val].is_square()``; without that method it raises
      // AttributeError, which callers such as ``LaurentSeries.is_square`` treat
      // as "cannot decide".
      throw new NotImplementedError(
        'SAGE_NOT_IMPLEMENTED: is_square of a power series whose base ring elements have no is_square method'
      );
    }
    if (!leadAny.is_square()) {
      return false;
    }
    if (this._parent.base_ring().is_field?.()) {
      return true;
    }
    try {
      this._parent.__call__(this.sqrt());
      return true;
    } catch {
      return false;
    }
  }

  /**
   * Return the series with precision truncated to n.
   * @see Reference: sage/rings/power_series_ring_element.pyx:add_bigoh
   */
  add_bigoh(n: number): PowerSeriesElement<T> {
    if (n === Number.POSITIVE_INFINITY || n > this._prec) {
      return this;
    }
    const coeffs = this._coefficients.slice(0, n);
    return new PowerSeriesElement<T>(this._parent, coeffs, n);
  }

  /**
   * Return the derivative.
   * @see Reference: sage/rings/power_series_ring_element.pyx:derivative
   */
  derivative(): PowerSeriesElement<T> {
    if (this._coefficients.length <= 1) {
      return new PowerSeriesElement<T>(
        this._parent,
        [],
        this._prec === Number.POSITIVE_INFINITY ? Number.POSITIVE_INFINITY : this._prec - 1
      );
    }

    const baseRing = this._parent.base_ring();
    const newCoeffs: T[] = [];
    for (let i = 1; i < this._coefficients.length; i++) {
      // Multiply coefficient by i
      let coeff = this._coefficients[i]!;
      for (let j = 1; j < i; j++) {
        coeff = coeff.add(this._coefficients[i]!) as T;
      }
      newCoeffs.push(coeff);
    }

    return new PowerSeriesElement<T>(
      this._parent,
      newCoeffs,
      this._prec === Number.POSITIVE_INFINITY ? Number.POSITIVE_INFINITY : this._prec - 1
    );
  }

  /**
   * Return the integral (with constant term 0).
   * @see Reference: sage/rings/power_series_ring_element.pyx:integral
   */
  integral(): PowerSeriesElement<T> {
    if (this._coefficients.length === 0) {
      return new PowerSeriesElement<T>(
        this._parent,
        [],
        this._prec === Number.POSITIVE_INFINITY ? Number.POSITIVE_INFINITY : this._prec + 1
      );
    }

    const baseRing = this._parent.base_ring();
    const newCoeffs: T[] = [baseRing.zero()]; // constant term is 0

    for (let i = 0; i < this._coefficients.length; i++) {
      // Divide coefficient by (i+1)
      // For rational coefficients, we need division
      const coeff = this._coefficients[i]!;
      // Create (i+1) as a ring element and divide
      const divisor = baseRing.one();
      let sum = divisor;
      for (let j = 0; j < i; j++) {
        sum = sum.add(divisor) as T;
      }
      newCoeffs.push(coeff.div(sum) as T);
    }

    return new PowerSeriesElement<T>(
      this._parent,
      newCoeffs,
      this._prec === Number.POSITIVE_INFINITY ? Number.POSITIVE_INFINITY : this._prec + 1
    );
  }

  /**
   * Return the exponential exp(self).
   * Requires the constant term to be 0 (or requires special handling).
   * @see Reference: sage/rings/power_series_ring_element.pyx:exp
   */
  exp(prec?: number): PowerSeriesElement<T> {
    const targetPrec = prec ?? this._parent.default_prec();
    const computePrec = Math.min(targetPrec, this._prec);

    // Check for non-zero constant term
    const c0 =
      this._coefficients.length > 0 ? this._coefficients[0]! : this._parent.base_ring().zero();
    if (!c0.isZero()) {
      throw new ArithmeticError(
        'can only compute exp of power series with zero constant term (or use a ring that supports exp of the constant term)'
      );
    }

    // Use the differential equation: if f = exp(self), then f' = self' * f
    // Start with f = 1
    const baseRing = this._parent.base_ring();
    const one = baseRing.one();
    const zero = baseRing.zero();

    // Initialize result with 1
    const result: T[] = [one];

    // Compute coefficients iteratively
    // If g = this and f = exp(g), then f' = g' * f
    // So f[n] = (1/n) * sum_{k=0}^{n-1} (n-k) * g'[n-1-k] * f[k]
    // = (1/n) * sum_{k=0}^{n-1} (n-k) * (n-k) * g[n-k] * f[k]
    // Wait, let's use a simpler approach:
    // f' = g' * f means:
    // (n+1) * f[n+1] = sum_{k=0}^{n} g'[k] * f[n-k]
    // where g'[k] = (k+1) * g[k+1]
    // So: (n+1) * f[n+1] = sum_{k=0}^{n} (k+1) * g[k+1] * f[n-k]

    const selfDerivCoeffs: T[] = [];
    for (let i = 1; i < Math.min(this._coefficients.length, computePrec); i++) {
      // i * this._coefficients[i]
      let coeff = this._coefficients[i]!;
      for (let j = 1; j < i; j++) {
        coeff = coeff.add(this._coefficients[i]!) as T;
      }
      selfDerivCoeffs.push(coeff);
    }

    for (let n = 0; n < computePrec - 1; n++) {
      // Compute f[n+1]
      let sum = zero;
      for (let k = 0; k <= Math.min(n, selfDerivCoeffs.length - 1); k++) {
        const gDeriv = selfDerivCoeffs[k] ?? zero;
        const fCoeff = result[n - k] ?? zero;
        sum = sum.add(gDeriv.mul(fCoeff)) as T;
      }
      // Divide by (n+1)
      let divisor = one;
      for (let j = 0; j < n; j++) {
        divisor = divisor.add(one) as T;
      }
      result.push(sum.div(divisor) as T);
    }

    return new PowerSeriesElement<T>(this._parent, result, computePrec);
  }

  /**
   * Return the logarithm log(self).
   * Requires the constant term to be 1.
   * @see Reference: sage/rings/power_series_ring_element.pyx:log
   */
  log(prec?: number): PowerSeriesElement<T> {
    const targetPrec = prec ?? this._parent.default_prec();

    // Check that constant term is 1
    const c0 =
      this._coefficients.length > 0 ? this._coefficients[0]! : this._parent.base_ring().zero();
    if (c0.isZero() || !c0.eq(1)) {
      throw new ArithmeticError('constant term of power series is not 1');
    }

    // SageMath computes ``zero.solve_linear_de(prec, b=self.derivative()/self, f0=0)``,
    // i.e. the solution of t' = f'/f with t(0) = 0, which is the integral of
    // the logarithmic derivative.
    // Reference: sage/rings/power_series_ring_element.pyx:log
    const t = this.derivative().div(this).integral();
    return t.add_bigoh(targetPrec);
  }

  /**
   * Return the square root.
   * @see Reference: sage/rings/power_series_ring_element.pyx:sqrt
   */
  sqrt(prec?: number): PowerSeriesElement<T> {
    const targetPrec = prec ?? this._parent.default_prec();
    const computePrec = Math.min(targetPrec, this._prec);

    if (this.is_zero()) {
      const newPrec =
        this._prec === Number.POSITIVE_INFINITY
          ? Number.POSITIVE_INFINITY
          : Math.floor(this._prec / 2);
      return new PowerSeriesElement<T>(this._parent, [], newPrec);
    }

    const val = this.valuation();
    if (val === Number.POSITIVE_INFINITY) {
      return new PowerSeriesElement<T>(this._parent, [], Number.POSITIVE_INFINITY);
    }

    if (val % 2 !== 0) {
      throw new ValueError('power series does not have a square root since it has odd valuation');
    }

    const baseRing = this._parent.base_ring();
    const one = baseRing.one();
    // SageMath uses ``half = ~R(2)``, which for R = ZZ lands in QQ.  This port
    // cannot change the coefficient ring on the fly, so ``1/2`` is only formed
    // when a coefficient actually needs it (never for an exact square such as
    // ``(t^-4).valuation_zero_part() == 1``).
    // Reference: sage/rings/power_series_ring_element.pyx:1770
    let _half: T | null = null;
    const half = (): T => {
      if (_half === null) {
        _half = one.div(one.add(one)) as T;
      }
      return _half;
    };

    // Get the valuation zero part
    const valuationZeroPart = val > 0 ? this._shiftRight(val) : this;
    const c0 = valuationZeroPart.__getitem__(0);

    // Try to compute square root of constant term
    // If the base ring supports it, we can handle any square constant term
    let c0Sqrt: T;
    try {
      // Try to call sqrt on the base ring element
      if ('sqrt' in c0 && typeof (c0 as unknown as { sqrt: () => T }).sqrt === 'function') {
        c0Sqrt = (c0 as unknown as { sqrt: () => T }).sqrt();
      } else if (c0.eq(1)) {
        c0Sqrt = one;
      } else {
        throw new ValueError(`unable to take the square root of ${c0}`);
      }
    } catch {
      throw new ValueError(`unable to take the square root of ${c0}`);
    }

    // Use Newton's method: x_{n+1} = (x_n + f/x_n) / 2
    // Or equivalently: x_{n+1} = x_n * (3 - f * x_n^2) / 2 for sqrt
    // Or: compute using series expansion directly
    // sqrt(1 + g) = 1 + g/2 - g^2/8 + g^3/16 - 5g^4/128 + ...
    // Using the formula: sqrt(1+g) = sum_{n>=0} C(1/2, n) * g^n
    // where C(1/2, n) = (1/2)(1/2-1)(1/2-2)...(1/2-n+1) / n!

    // Normalize to sqrt(c0) * sqrt(1 + g) where g = (f/c0 - 1)
    const normalizedPart = valuationZeroPart._scalarDiv(c0);
    const g = normalizedPart.sub(this._parent.one()).add_bigoh(computePrec);
    const result: T[] = [c0Sqrt];

    // Compute coefficients using the recurrence for sqrt(c0) * sqrt(1+g)
    // If sqrt(1+g) = sum b_n x^n and g = sum a_n x^n (with a_0 = 0)
    // Then 2*b_0 * b_n = a_n - sum_{k=1}^{n-1} b_k * b_{n-k}
    // Since b_0 = 1 for the normalized sqrt(1+g), and we scale by sqrt(c0)
    const aCoeffs = g.list();
    const c0SqrtInv = one.div(c0Sqrt) as T;

    // Compute intermediate coefficients for sqrt(1+g)
    // g = a_0 + a_1*x + a_2*x^2 + ... where aCoeffs = [a_0, a_1, a_2, ...]
    // Since g = f/c0 - 1 and f starts with c0, we have a_0 = 0
    // y = sqrt(1+g) = b_0 + b_1*x + b_2*x^2 + ... with b_0 = 1
    // y^2 = 1 + g gives: 2*b_0*b_n + sum_{k=1}^{n-1} b_k*b_{n-k} = a_n
    // So: b_n = (a_n - sum_{k=1}^{n-1} b_k*b_{n-k}) / 2
    const bCoeffs: T[] = [one];
    for (let n = 1; n < computePrec; n++) {
      // a_n is the coefficient of x^n in g, which is at index n in aCoeffs
      const an = n < aCoeffs.length ? (aCoeffs[n] ?? baseRing.zero()) : baseRing.zero();
      let sum = an;
      for (let k = 1; k < n; k++) {
        const bk = bCoeffs[k]!;
        const bn_k = bCoeffs[n - k]!;
        sum = sum.sub(bk.mul(bn_k)) as T;
      }
      bCoeffs.push(sum.isZero() ? sum : (sum.mul(half()) as T));
    }

    // Scale by sqrt(c0)
    for (let n = 1; n < bCoeffs.length; n++) {
      result.push(c0Sqrt.mul(bCoeffs[n]!) as T);
    }

    let finalResult = new PowerSeriesElement<T>(this._parent, result, computePrec);

    // Shift back if original had valuation > 0
    if (val > 0) {
      finalResult = finalResult._shiftLeft(val / 2);
    }

    // SageMath's ``test_exact``: when the input was exact and the answer is
    // short enough to be checked, an exact square root is returned exactly.
    // Reference: sage/rings/power_series_ring_element.pyx:1782
    if (this._prec === Number.POSITIVE_INFINITY && finalResult.degree() < computePrec / 2) {
      const sq = finalResult.mul(finalResult);
      if (sq.sub(this).is_zero() && sq.prec() >= this.degree() + 1) {
        finalResult = finalResult.truncate();
      }
    }

    return finalResult;
  }

  /**
   * Return the series with precision truncated to `n` (alias of `add_bigoh`).
   * @see Reference: sage/rings/power_series_ring_element.pyx:O
   */
  O(n: number): PowerSeriesElement<T> {
    return this.add_bigoh(n);
  }

  /**
   * Return the n-th root.
   * @see Reference: sage/rings/power_series_ring_element.pyx:nth_root
   */
  nth_root(n: number, prec?: number): PowerSeriesElement<T> {
    // Reference: sage/rings/power_series_ring_element.pyx:nth_root
    const val = this.valuation();

    if (this.is_zero()) {
      if (val === Number.POSITIVE_INFINITY) {
        return this;
      }
      if (n <= 0) {
        throw new ValueError(`n (=${n}) must be positive`);
      }
      return this._parent.zero().add_bigoh(Math.floor(val / n));
    }

    if (n <= 0) {
      // SageMath reaches Polynomial._nth_root_series, which rejects n <= 0.
      throw new ValueError(`n (=${n}) must be positive`);
    }

    if (val !== Number.POSITIVE_INFINITY && val % n !== 0) {
      throw new ValueError(`power series valuation is not a multiple of ${n}`);
    }

    const maxprec = Math.floor(val / n) + this.precision_relative();
    let targetPrec: number;
    if (prec === undefined) {
      targetPrec = maxprec === Number.POSITIVE_INFINITY ? this._parent.default_prec() : maxprec;
    } else {
      targetPrec = Math.min(maxprec, prec);
    }

    const p = this.truncate();
    const q = p._nth_root_series(n, targetPrec);
    let ans = q;
    if (
      !(
        this._prec === Number.POSITIVE_INFINITY &&
        q.degree() * n <= targetPrec &&
        q.pow(n).sub(p).is_zero()
      )
    ) {
      ans = ans.add_bigoh(targetPrec);
    }
    return ans;
  }

  /**
   * Return the first ``prec`` coefficients of the ``n``-th root series of this
   * polynomial (i.e. of this series regarded as an exact polynomial).
   *
   * ALGORITHM: Newton's method for the fixed point of `F(x) = x^{-n} - a^{-1}`,
   * which requires only one series inversion at the very end.
   *
   * @see Reference: sage/rings/polynomial/polynomial_element.pyx:_nth_root_series
   */
  private _nth_root_series(n: number, prec: number, start?: T): PowerSeriesElement<T> {
    const baseRing = this._parent.base_ring();
    const m = n;
    if (m <= 0) {
      throw new ValueError(`n (=${m}) must be positive`);
    }
    if (m === 1 || this.is_zero() || this.is_one()) {
      return this;
    }

    const c0 = this.__getitem__(0);
    if (c0.isZero()) {
      // p = x^i q, so p^(1/m) = x^(i/m) q^(1/m)
      const i = this.valuation();
      if (i % m !== 0) {
        throw new ValueError(`not a ${m}th power`);
      }
      return this._shiftRight(i)
        ._nth_root_series(m, prec - i / m)
        ._shiftLeft(i / m);
    }

    // SageMath additionally handles the case where the characteristic of the
    // base ring divides n; that requires n-th roots of the coefficients.
    const characteristic = baseRing.characteristic ? baseRing.characteristic() : 0n;
    if (characteristic !== 0n && BigInt(m) % characteristic === 0n) {
      throw new NotImplementedError(
        'SAGE_NOT_IMPLEMENTED: nth_root_series when the characteristic of the base ring divides n'
      );
    }

    const one = baseRing.one();

    // The constant term of the root
    let a: T;
    if (start !== undefined) {
      a = start;
    } else if (c0.isOne ? c0.isOne() : c0.eq(1)) {
      a = one;
    } else {
      const c0AsAny = c0 as unknown as { nth_root?: (k: number) => T };
      if (typeof c0AsAny.nth_root === 'function') {
        a = c0AsAny.nth_root(m);
      } else {
        throw new NotImplementedError(
          `SAGE_NOT_IMPLEMENTED: nth root of the constant coefficient ${c0}`
        );
      }
    }

    let mi: T;
    try {
      mi = one.div(baseRing.__call__(BigInt(m))) as T;
    } catch {
      throw new ArithmeticError('exponent not invertible in base ring');
    }

    let q = this._parent.__call__([one.div(a) as T]);
    const mp1 = m + 1;
    const mp1El = baseRing.__call__(BigInt(mp1));
    for (const i of newton_method_sizes(prec)) {
      // q = mi * ((m+1)*q - p * q^(m+1))   truncated at x^i
      const qPow = q._power_trunc(mp1, i);
      const rhs = q._scalarMul(mp1El).sub(this._mul_trunc(qPow, i));
      q = rhs.truncate(i)._scalarMul(mi);
    }
    return q.add_bigoh(prec).inv().truncate(prec);
  }

  /**
   * Return `self * other` truncated at `x^n` (as an exact polynomial).
   * @see Reference: sage/rings/polynomial/polynomial_element.pyx:_mul_trunc_
   */
  private _mul_trunc(other: PowerSeriesElement<T>, n: number): PowerSeriesElement<T> {
    return this.add_bigoh(n).mul(other.add_bigoh(n)).truncate(n);
  }

  /**
   * Return `self^e` truncated at `x^n` (as an exact polynomial).
   * @see Reference: sage/rings/polynomial/polynomial_element.pyx:_power_trunc
   */
  private _power_trunc(e: number, n: number): PowerSeriesElement<T> {
    let result = this._parent.one().add_bigoh(n);
    let base = this.add_bigoh(n);
    let k = e;
    while (k > 0) {
      if (k & 1) {
        result = result.mul(base).add_bigoh(n);
      }
      base = base.mul(base).add_bigoh(n);
      k >>= 1;
    }
    return result.truncate(n);
  }

  /**
   * Evaluate the series at `x = a`, i.e. return the composition `f(a)` where
   * `f = this`.
   *
   * To substitute a value it must have valuation at least 1, unless `self` has
   * infinite precision (i.e. is a polynomial).
   *
   * ALGORITHM: SageMath's `PowerSeries_poly.__call__`: truncate the argument
   * to the precision `(s - r + 1) t` that the answer can see (`s` the precision
   * of `self`, `r` the valuation of `self - self[0]`, `t` the valuation of the
   * argument), then evaluate the *polynomial* part at the argument -- power
   * series multiplication then produces the correct precision by itself.
   *
   * NOTE: SageMath returns a base ring element when `a` has infinite valuation;
   * this port always returns an element of `a`'s parent (the constant series),
   * because TypeScript needs a single return type.
   *
   * @see Reference: sage/rings/power_series_poly.pyx:176 (__call__)
   */
  __call__(a: PowerSeriesElement<T>): PowerSeriesElement<T> {
    const s = this._prec;
    if (s === Number.POSITIVE_INFINITY) {
      return this._polynomialEval(a);
    }

    const t = a.valuation();

    if (t === Number.POSITIVE_INFINITY) {
      return a.parent().__call__(this.__getitem__(0));
    }

    if (t <= 0) {
      throw new ValueError('Can only substitute elements of positive valuation');
    }

    const r = this.sub(this._parent.__call__(this.__getitem__(0))).valuation();
    if (r === s) {
      // self is constant + O(x^s)
      return a
        .parent()
        .__call__(this.__getitem__(0))
        .add_bigoh(s * t);
    }

    const u = a.prec();
    const n = (s - r + 1) * t;
    let arg = a;
    if (n < u) {
      arg = a.add_bigoh(n);
    }
    return this._polynomialEval(arg);
  }

  /**
   * Evaluate the (exact) polynomial part of this series at `a` by Horner's
   * rule, in the parent of `a`.
   *
   * @see Reference: sage/rings/power_series_poly.pyx:176 (__call__, ``self.__f(x)``)
   */
  private _polynomialEval(a: PowerSeriesElement<T>): PowerSeriesElement<T> {
    const P = a.parent();
    let result = P.zero();
    for (let i = this._coefficients.length - 1; i >= 0; i--) {
      result = result.mul(a).add(P.__call__(this._coefficients[i]!));
    }
    return result;
  }

  /**
   * Return the reversion (compositional inverse).
   * If this = f, return g such that f(g(x)) = x.
   * Requires f to have valuation 1.
   * @see Reference: sage/rings/power_series_poly.pyx:reverse
   */
  reversion(prec?: number): PowerSeriesElement<T> {
    if (this.valuation() !== 1) {
      throw new ValueError('Series must have valuation one for reversion.');
    }

    const targetPrec =
      prec ?? (this._prec === Number.POSITIVE_INFINITY ? this._parent.default_prec() : this._prec);
    const computePrec = Math.min(targetPrec, this._prec);

    // Use Lagrange inversion formula
    // If f(x) = a_1*x + a_2*x^2 + ..., then g(x) = sum_{n>=1} b_n * x^n
    // where b_n = (1/n) * [t^{n-1}] (t/f(t))^n
    const baseRing = this._parent.base_ring();
    const one = baseRing.one();
    const x = this._parent.gen();

    // First, extract the leading coefficient a_1 and check it's a unit
    const a1 = this.__getitem__(1);
    if (a1.isZero()) {
      throw new ValueError('Series must have non-zero linear coefficient for reversion.');
    }

    // h = x / f = 1/(a_1 + a_2*x + a_3*x^2 + ...)
    // To compute x/f, we first shift f right by 1, then invert, then the x cancels
    // f = x * (a_1 + a_2*x + ...) = x * f_shifted
    // x/f = 1/f_shifted
    const fShifted = this._shiftRight(1);
    const h = fShifted.inv().add_bigoh(computePrec);

    const resultCoeffs: T[] = [baseRing.zero()]; // g has no constant term

    // h^n gives coefficients we need
    let hPower = h.add_bigoh(computePrec);

    for (let n = 1; n < computePrec; n++) {
      // Get coefficient of x^{n-1} in h^n
      let coeff: T;
      if (n === 1) {
        // h^1, get coefficient of x^0
        coeff = hPower.__getitem__(0);
      } else {
        coeff = hPower.__getitem__(n - 1);
      }
      // Divide by n
      let divisor = one;
      for (let j = 1; j < n; j++) {
        divisor = divisor.add(one) as T;
      }
      resultCoeffs.push(coeff.div(divisor) as T);

      // Compute next power of h
      if (n < computePrec - 1) {
        hPower = hPower.mul(h).add_bigoh(computePrec);
      }
    }

    return new PowerSeriesElement<T>(this._parent, resultCoeffs, computePrec);
  }

  /**
   * Return the Pade approximant [m/n].
   * Returns [Q, P] such that deg(Q) <= m, deg(P) <= n, and f - Q/P = O(x^{m+n+1})
   * @see Reference: sage/rings/power_series_poly.pyx:pade
   */
  pade(m: number, n: number): PadeApproximant<T> {
    // Reference: sage/rings/power_series_poly.pyx:pade
    // The Pade approximant uses the extended Euclidean algorithm (rational reconstruction)

    if (this.precision_absolute() < n + m + 1) {
      throw new ValueError('the precision of the series is not large enough');
    }

    const baseRing = this._parent.base_ring();
    const one = baseRing.one();
    const zero = baseRing.zero();

    // We need to find Q, P such that Q/P ≡ f (mod x^{m+n+1})
    // This is equivalent to: Q ≡ f * P (mod x^{m+n+1})
    // Using extended GCD on f and x^{m+n+1}

    // Represent polynomials as coefficient arrays
    // f = sum f[i] x^i
    const fCoeffs: T[] = [];
    for (let i = 0; i <= m + n; i++) {
      fCoeffs.push(this.__getitem__(i));
    }

    // Extended Euclidean algorithm to find (Q, P) with deg(Q) <= m, deg(P) <= n
    // such that Q ≡ f * P (mod x^{m+n+1})

    // Start with: r0 = x^{m+n+1}, r1 = f (truncated)
    // t0 = 0, t1 = 1
    // We iterate until deg(r_k) <= m

    // Build x^{m+n+1}
    const modCoeffs: T[] = new Array(m + n + 1).fill(zero);
    modCoeffs.push(one);

    // Use polynomial division algorithm
    // r = [r0, r1], t = [t0, t1]
    let r0 = modCoeffs;
    let r1 = [...fCoeffs];
    let t0: T[] = [zero];
    let t1: T[] = [one];

    // Helper: polynomial degree
    const polyDeg = (p: T[]): number => {
      for (let i = p.length - 1; i >= 0; i--) {
        if (!p[i]!.isZero()) return i;
      }
      return -1;
    };

    // Helper: polynomial leading coefficient
    const polyLC = (p: T[]): T => {
      const d = polyDeg(p);
      return d < 0 ? zero : p[d]!;
    };

    // Helper: polynomial subtraction
    const polySub = (a: T[], b: T[]): T[] => {
      const result: T[] = [];
      const len = Math.max(a.length, b.length);
      for (let i = 0; i < len; i++) {
        const ai = i < a.length ? a[i]! : zero;
        const bi = i < b.length ? b[i]! : zero;
        result.push(ai.sub(bi) as T);
      }
      return result;
    };

    // Helper: polynomial scalar multiplication
    const polyScalarMul = (a: T[], c: T): T[] => a.map((x) => x.mul(c) as T);

    // Helper: polynomial shift (multiply by x^k)
    const polyShift = (a: T[], k: number): T[] => {
      const result: T[] = new Array(k).fill(zero);
      return result.concat(a);
    };

    // Extended Euclidean algorithm
    while (polyDeg(r1) > m) {
      const degR0 = polyDeg(r0);
      const degR1 = polyDeg(r1);
      if (degR1 < 0) break;

      const shift = degR0 - degR1;
      const factor = polyLC(r0).div(polyLC(r1)) as T;

      const r1Shifted = polyShift(r1, shift);
      const t1Shifted = polyShift(t1, shift);

      const newR = polySub(r0, polyScalarMul(r1Shifted, factor));
      const newT = polySub(t0, polyScalarMul(t1Shifted, factor));

      r0 = r1;
      r1 = newR;
      t0 = t1;
      t1 = newT;
    }

    // Q = r1, P = t1 (up to normalization)
    // Trim trailing zeros
    while (r1.length > 0 && r1[r1.length - 1]!.isZero()) r1.pop();
    while (t1.length > 0 && t1[t1.length - 1]!.isZero()) t1.pop();

    if (t1.length === 0) {
      t1 = [one];
    }

    // SageMath returns ``u / v`` in the fraction field of the polynomial ring,
    // which normalizes the denominator to be monic.
    const lc = polyLC(t1);
    const numer = r1.length === 0 ? [zero] : r1.map((c) => c.div(lc) as T);
    const denom = t1.map((c) => c.div(lc) as T);

    const Q = new PowerSeriesElement<T>(this._parent, numer, Number.POSITIVE_INFINITY);
    const P = new PowerSeriesElement<T>(this._parent, denom, Number.POSITIVE_INFINITY);

    return new PadeApproximant<T>(Q, P);
  }

  // Arithmetic operations

  /**
   * Add two power series.
   * @see Reference: sage/rings/power_series_ring_element.pyx:__add__
   */
  add(other: PowerSeriesElement<T>): PowerSeriesElement<T> {
    const newPrec = Math.min(this._prec, other._prec);
    const maxLen = Math.max(this._coefficients.length, other._coefficients.length);
    const baseRing = this._parent.base_ring();
    const newCoeffs: T[] = [];

    for (let i = 0; i < maxLen; i++) {
      const a = i < this._coefficients.length ? this._coefficients[i]! : baseRing.zero();
      const b = i < other._coefficients.length ? other._coefficients[i]! : baseRing.zero();
      newCoeffs.push(a.add(b) as T);
    }

    return new PowerSeriesElement<T>(this._parent, newCoeffs, newPrec);
  }

  /**
   * Subtract two power series.
   * @see Reference: sage/rings/power_series_ring_element.pyx:__sub__
   */
  sub(other: PowerSeriesElement<T>): PowerSeriesElement<T> {
    const newPrec = Math.min(this._prec, other._prec);
    const maxLen = Math.max(this._coefficients.length, other._coefficients.length);
    const baseRing = this._parent.base_ring();
    const newCoeffs: T[] = [];

    for (let i = 0; i < maxLen; i++) {
      const a = i < this._coefficients.length ? this._coefficients[i]! : baseRing.zero();
      const b = i < other._coefficients.length ? other._coefficients[i]! : baseRing.zero();
      newCoeffs.push(a.sub(b) as T);
    }

    return new PowerSeriesElement<T>(this._parent, newCoeffs, newPrec);
  }

  /**
   * Multiply two power series.
   * @see Reference: sage/rings/power_series_ring_element.pyx:__mul__
   */
  mul(other: PowerSeriesElement<T>): PowerSeriesElement<T> {
    if (this.is_zero() || other.is_zero()) {
      return this._parent.zero().add_bigoh(this._computeMulPrec(other));
    }

    const selfVal = this.valuation();
    const otherVal = other.valuation();
    let newPrec: number;
    if (this._prec === Number.POSITIVE_INFINITY) {
      if (other._prec === Number.POSITIVE_INFINITY) {
        newPrec = Number.POSITIVE_INFINITY;
      } else {
        newPrec = other._prec + selfVal;
      }
    } else {
      if (other._prec === Number.POSITIVE_INFINITY) {
        newPrec = this._prec + otherVal;
      } else {
        newPrec = Math.min(other._prec + selfVal, this._prec + otherVal);
      }
    }

    const maxDeg =
      newPrec === Number.POSITIVE_INFINITY
        ? this._coefficients.length + other._coefficients.length - 1
        : newPrec - 1;

    const baseRing = this._parent.base_ring();
    const newCoeffs: T[] = [];

    for (let k = 0; k <= maxDeg; k++) {
      let sum = baseRing.zero();
      const jMin = Math.max(0, k - other._coefficients.length + 1);
      const jMax = Math.min(k, this._coefficients.length - 1);
      for (let j = jMin; j <= jMax; j++) {
        const a = this._coefficients[j]!;
        const b = other._coefficients[k - j]!;
        sum = sum.add(a.mul(b)) as T;
      }
      newCoeffs.push(sum);
    }

    return new PowerSeriesElement<T>(this._parent, newCoeffs, newPrec);
  }

  private _computeMulPrec(other: PowerSeriesElement<T>): number {
    const selfVal = this.valuation();
    const otherVal = other.valuation();
    if (this._prec === Number.POSITIVE_INFINITY) {
      if (other._prec === Number.POSITIVE_INFINITY) {
        return Number.POSITIVE_INFINITY;
      }
      return other._prec + selfVal;
    }
    if (other._prec === Number.POSITIVE_INFINITY) {
      return this._prec + otherVal;
    }
    return Math.min(other._prec + selfVal, this._prec + otherVal);
  }

  /**
   * Divide two power series.
   * @see Reference: sage/rings/power_series_ring_element.pyx:__truediv__
   */
  div(other: PowerSeriesElement<T>): PowerSeriesElement<T> {
    if (other.is_zero()) {
      throw new ZeroDivisionError("Can't divide by something indistinguishable from 0");
    }

    const otherVal = other.valuation();
    const selfVal = this.valuation();

    if (otherVal > selfVal) {
      // Would need Laurent series
      throw new NotImplementedError('Division would produce Laurent series (negative powers)');
    }

    // Cancel common factors and invert the denominator
    const shiftedSelf = otherVal > 0 ? this._shiftRight(otherVal) : this;
    const valuationZeroPart = other._shiftRight(otherVal);
    const inv = valuationZeroPart.inv();

    return shiftedSelf.mul(inv);
  }

  /**
   * Return the negation.
   * @see Reference: sage/rings/power_series_ring_element.pyx:__neg__
   */
  neg(): PowerSeriesElement<T> {
    const newCoeffs: T[] = this._coefficients.map((c) => c.neg() as T);
    return new PowerSeriesElement<T>(this._parent, newCoeffs, this._prec);
  }

  /**
   * Return the multiplicative inverse.
   * @see Reference: sage/rings/power_series_ring_element.pyx:__invert__
   */
  inv(): PowerSeriesElement<T> {
    // SageMath: ``if self.is_one(): return self`` -- the inverse of an exact 1
    // is exact.
    // Reference: sage/rings/power_series_poly.pyx:705
    if (this.is_one()) {
      return this;
    }
    if (this.is_zero()) {
      throw new ZeroDivisionError('Power series is not invertible (constant term is not a unit)');
    }
    if (!this.is_unit()) {
      // SageMath ends up in ``Polynomial.inverse_series_trunc``, which reports
      // the offending constant term.
      // Reference: sage/rings/polynomial/polynomial_element.pyx:1773
      throw new ValueError(`constant term ${this._coefficients[0]!} is not a unit`);
    }

    const computePrec =
      this._prec === Number.POSITIVE_INFINITY ? this._parent.default_prec() : this._prec;
    const baseRing = this._parent.base_ring();
    const c0 = this._coefficients[0]!;

    // Compute inverse using the formula:
    // If f = c0 * (1 + g) where g has positive valuation,
    // then 1/f = (1/c0) * sum_{k>=0} (-g)^k
    // More directly, use the recurrence:
    // If f = sum a_n x^n and 1/f = sum b_n x^n,
    // then b_0 = 1/a_0 and b_n = -(1/a_0) * sum_{k=1}^{n} a_k * b_{n-k}

    const c0Inv = c0.inv ? (c0.inv() as T) : (baseRing.one().div(c0) as T);
    const result: T[] = [c0Inv];

    for (let n = 1; n < computePrec; n++) {
      let sum = baseRing.zero();
      for (let k = 1; k <= Math.min(n, this._coefficients.length - 1); k++) {
        const ak = this._coefficients[k]!;
        const bn_k = result[n - k]!;
        sum = sum.add(ak.mul(bn_k)) as T;
      }
      result.push(sum.neg().mul(c0Inv) as T);
    }

    return new PowerSeriesElement<T>(this._parent, result, computePrec);
  }

  /**
   * Return self raised to power n.
   * @see Reference: sage/rings/power_series_ring_element.pyx:__pow__
   */
  pow(n: bigint | number): PowerSeriesElement<T> {
    const exp = typeof n === 'number' ? BigInt(n) : n;

    if (exp === 0n) {
      return this._parent.one().add_bigoh(this._prec);
    }

    if (exp < 0n) {
      return this.inv().pow(-exp);
    }

    // Binary exponentiation
    let result = this._parent.one().add_bigoh(this._prec);
    let base: PowerSeriesElement<T> = this;
    let e = exp;

    while (e > 0n) {
      if ((e & 1n) === 1n) {
        result = result.mul(base);
      }
      base = base.mul(base);
      e >>= 1n;
    }

    return result;
  }

  // Helper methods

  /**
   * Multiply by a scalar.
   * @see Reference: sage/rings/power_series_ring_element.pyx
   */
  _scalarMul(scalar: T): PowerSeriesElement<T> {
    const newCoeffs = this._coefficients.map((c) => c.mul(scalar) as T);
    return new PowerSeriesElement<T>(this._parent, newCoeffs, this._prec);
  }

  /**
   * Multiply by a scalar (public alias).
   * @see Reference: sage/rings/power_series_ring_element.pyx
   */
  scalar_mul(scalar: T): PowerSeriesElement<T> {
    return this._scalarMul(scalar);
  }

  /**
   * Divide by a scalar.
   * @see Reference: sage/rings/power_series_ring_element.pyx
   */
  _scalarDiv(scalar: T): PowerSeriesElement<T> {
    const newCoeffs = this._coefficients.map((c) => c.div(scalar) as T);
    return new PowerSeriesElement<T>(this._parent, newCoeffs, this._prec);
  }

  /**
   * Divide by a scalar (public alias).
   * @see Reference: sage/rings/power_series_ring_element.pyx
   */
  scalar_div(scalar: T): PowerSeriesElement<T> {
    return this._scalarDiv(scalar);
  }

  /**
   * Shift right (divide by x^n, discarding terms below x^n).
   *
   * A negative `n` shifts left, exactly as SageMath's `f >> n`.
   *
   * @see Reference: sage/rings/power_series_poly.pyx:598 (__rshift__)
   */
  _shiftRight(n: number): PowerSeriesElement<T> {
    if (n === 0) return this;
    if (n < 0) return this._shiftLeft(-n);
    const newCoeffs = this._coefficients.slice(n);
    const newPrec =
      this._prec === Number.POSITIVE_INFINITY
        ? Number.POSITIVE_INFINITY
        : Math.max(0, this._prec - n);
    return new PowerSeriesElement<T>(this._parent, newCoeffs, newPrec);
  }

  /**
   * Shift left (multiply by x^n).
   *
   * A negative `n` shifts right.  SageMath's `f << n` would give the resulting
   * series the precision `prec + n`, which can be negative; this port routes
   * negative shifts through `f >> -n`, whose precision SageMath clamps at 0.
   *
   * @see Reference: sage/rings/power_series_poly.pyx:582 (__lshift__)
   */
  _shiftLeft(n: number): PowerSeriesElement<T> {
    if (n === 0) return this;
    if (n < 0) return this._shiftRight(-n);
    const baseRing = this._parent.base_ring();
    const newCoeffs: T[] = [];
    for (let i = 0; i < n; i++) {
      newCoeffs.push(baseRing.zero());
    }
    newCoeffs.push(...this._coefficients);
    return new PowerSeriesElement<T>(
      this._parent,
      newCoeffs,
      this._prec === Number.POSITIVE_INFINITY ? Number.POSITIVE_INFINITY : this._prec + n
    );
  }

  toString(): string {
    if (this._coefficients.length === 0) {
      if (this._prec === Number.POSITIVE_INFINITY) {
        return '0';
      }
      return `O(${this._parent.variable_name()}^${this._prec})`;
    }

    const varName = this._parent.variable_name();
    const terms: string[] = [];

    for (let i = 0; i < this._coefficients.length; i++) {
      const c = this._coefficients[i]!;
      if (c.isZero()) continue;

      let termStr: string;
      const coeffStr = c.toString();

      if (i === 0) {
        termStr = coeffStr;
      } else if (i === 1) {
        if (c.eq(1)) {
          termStr = varName;
        } else if (c.eq(-1)) {
          termStr = `-${varName}`;
        } else {
          termStr = `${coeffStr}*${varName}`;
        }
      } else {
        if (c.eq(1)) {
          termStr = `${varName}^${i}`;
        } else if (c.eq(-1)) {
          termStr = `-${varName}^${i}`;
        } else {
          termStr = `${coeffStr}*${varName}^${i}`;
        }
      }
      terms.push(termStr);
    }

    let result = terms.join(' + ').replace(/\+ -/g, '- ');

    if (this._prec !== Number.POSITIVE_INFINITY) {
      if (result === '') {
        result = `O(${varName}^${this._prec})`;
      } else {
        result += ` + O(${varName}^${this._prec})`;
      }
    }

    return result || '0';
  }
}

/**
 * The Pade approximant `Q/P` of a power series.
 *
 * SageMath's `pade` returns an element of `Frac(R[z])`; this port has no
 * fraction field of a polynomial ring, so the quotient is represented by this
 * pair-with-accessors, normalised exactly as SageMath's fraction field
 * normalises it (monic denominator).
 *
 * @see Reference: sage/rings/power_series_poly.pyx:pade
 */
export class PadeApproximant<T extends RingElement = RingElement> {
  private readonly _numerator: PowerSeriesElement<T>;
  private readonly _denominator: PowerSeriesElement<T>;

  constructor(numerator: PowerSeriesElement<T>, denominator: PowerSeriesElement<T>) {
    this._numerator = numerator;
    this._denominator = denominator;
  }

  /** The numerator polynomial `Q`. */
  numerator(): PowerSeriesElement<T> {
    return this._numerator;
  }

  /** The denominator polynomial `P`. */
  denominator(): PowerSeriesElement<T> {
    return this._denominator;
  }

  /**
   * Expand `Q/P` as a power series to the given precision.
   */
  power_series(prec?: number): PowerSeriesElement<T> {
    const target = prec ?? this._numerator.parent().default_prec();
    return this._numerator.add_bigoh(target).div(this._denominator.add_bigoh(target));
  }

  /**
   * Render as SageMath renders an element of `Frac(R[z])`: polynomials in
   * descending degree order, denominator omitted when it is 1.
   */
  toString(): string {
    const num = this._polyStr(this._numerator.list());
    if (this._denominator.is_one()) {
      return num;
    }
    const den = this._polyStr(this._denominator.list());
    const wrap = (s: string): string => (s.includes(' ') ? `(${s})` : s);
    return `${wrap(num)}/${wrap(den)}`;
  }

  private _polyStr(coeffs: T[]): string {
    const varName = this._numerator.parent().variable_name();
    const terms: string[] = [];
    for (let i = coeffs.length - 1; i >= 0; i--) {
      const c = coeffs[i]!;
      if (c.isZero()) continue;
      let s: string;
      if (i === 0) {
        s = c.toString();
      } else if (c.eq(1)) {
        s = i === 1 ? varName : `${varName}^${i}`;
      } else if (c.eq(-1)) {
        s = i === 1 ? `-${varName}` : `-${varName}^${i}`;
      } else {
        s = i === 1 ? `${c}*${varName}` : `${c}*${varName}^${i}`;
      }
      terms.push(s);
    }
    if (terms.length === 0) {
      return '0';
    }
    return terms.join(' + ').replace(/\+ -/g, '- ');
  }
}

// ===========================================================================
// Multivariate power series
//
// Port of: sage/rings/multi_power_series_ring.py and
//          sage/rings/multi_power_series_ring_element.py
//
// SageMath represents a multivariate power series by a "background" univariate
// power series in an auxiliary variable T over the multivariate polynomial
// ring, where the coefficient of T^d is the degree-d homogeneous part
// (``MPowerSeries._bg_value``).  Precision is therefore a bound on the *total*
// degree, and the arithmetic precision rules are the univariate ones applied to
// that grading.  This port stores the same data directly as a map from
// exponent vectors to coefficients plus the total-degree precision; all
// precision rules below are the ones the background ring produces.
//
// Reference: reference/sage/src/sage/rings/multi_power_series_ring_element.py
// ===========================================================================

/** Total degree of an exponent vector. */
function _totalDegree(e: readonly number[]): number {
  let d = 0;
  for (const x of e) d += x;
  return d;
}

function _expKey(e: readonly number[]): string {
  return e.join(',');
}

function _keyExp(k: string): number[] {
  return k.split(',').map(Number);
}

/**
 * Multivariate power series ring `R[[x_1, ..., x_n]]`.
 *
 * @see Reference: sage/rings/multi_power_series_ring.py:MPowerSeriesRing_generic
 */
export class MPowerSeriesRing<T extends RingElement = RingElement> {
  private readonly _base_ring: CoefficientRing<T>;
  private readonly _names: string[];
  private readonly _default_prec: number;

  /**
   * @param base_ring - the coefficient ring
   * @param names - variable names, either a list or a comma separated string
   * @param default_prec - SageMath's default for multivariate rings is 10
   *
   * @see Reference: sage/rings/multi_power_series_ring.py:311 (__classcall__)
   */
  constructor(base_ring: CoefficientRing<T>, names: string | string[], default_prec: number = 10) {
    this._base_ring = base_ring;
    this._names = typeof names === 'string' ? names.split(',').map((s) => s.trim()) : [...names];
    if (this._names.length === 0) {
      throw new ValueError('multivariate power series rings must have at least one variable');
    }
    this._default_prec = default_prec;
  }

  base_ring(): CoefficientRing<T> {
    return this._base_ring;
  }

  /** @see Reference: sage/rings/multi_power_series_ring.py:variable_names */
  variable_names(): string[] {
    return [...this._names];
  }

  /** @see Reference: sage/rings/multi_power_series_ring.py:ngens */
  ngens(): number {
    return this._names.length;
  }

  /** @see Reference: sage/rings/multi_power_series_ring.py:default_prec */
  default_prec(): number {
    return this._default_prec;
  }

  characteristic(): bigint {
    return this._base_ring.characteristic?.() ?? 0n;
  }

  /** Whether `other` is the same parent (see {@link PowerSeriesRing.is_identical_to}). */
  is_identical_to(other: MPowerSeriesRing<T>): boolean {
    return (
      this === other ||
      (this._base_ring === other._base_ring &&
        this._names.length === other._names.length &&
        this._names.every((n, i) => n === other._names[i]) &&
        this._default_prec === other._default_prec)
    );
  }

  /** @see Reference: sage/rings/multi_power_series_ring.py:gen */
  gen(i: number = 0): MPowerSeries<T> {
    if (i < 0 || i >= this._names.length) {
      throw new ValueError('generator not defined');
    }
    const e = new Array(this._names.length).fill(0);
    e[i] = 1;
    return new MPowerSeries<T>(this, [[e, this._base_ring.one()]], Number.POSITIVE_INFINITY);
  }

  /** @see Reference: sage/rings/multi_power_series_ring.py:gens */
  gens(): MPowerSeries<T>[] {
    return this._names.map((_, i) => this.gen(i));
  }

  zero(): MPowerSeries<T> {
    return new MPowerSeries<T>(this, [], Number.POSITIVE_INFINITY);
  }

  one(): MPowerSeries<T> {
    return new MPowerSeries<T>(
      this,
      [[new Array(this._names.length).fill(0), this._base_ring.one()]],
      Number.POSITIVE_INFINITY
    );
  }

  /**
   * `R.O(prec)`: the zero series of precision `prec`.
   * @see Reference: sage/rings/multi_power_series_ring.py:O
   */
  O(prec: number): MPowerSeries<T> {
    if (prec === Number.POSITIVE_INFINITY) {
      return this.zero();
    }
    if (prec < 0) {
      throw new ValueError('prec (= %s) must be nonnegative'.replace('%s', String(prec)));
    }
    return new MPowerSeries<T>(this, [], prec);
  }

  /**
   * Convert `x` into this ring.  Accepts an element of this ring, an element of
   * the base ring, or a map from exponent vectors (as arrays or comma-joined
   * strings) to coefficients.
   *
   * @see Reference: sage/rings/multi_power_series_ring_element.py:308 (__init__)
   */
  __call__(x: unknown, prec: number = Number.POSITIVE_INFINITY): MPowerSeries<T> {
    if (x instanceof MPowerSeries) {
      const f = x as MPowerSeries<T>;
      if (f.parent().ngens() !== this.ngens()) {
        throw new TypeError('cannot coerce input to polynomial ring');
      }
      return new MPowerSeries<T>(this, f.monomial_coefficients(), Math.min(prec, f.prec()));
    }
    if (x instanceof Map) {
      const terms: [number[], T][] = [];
      for (const [k, c] of x as Map<string | number[], unknown>) {
        const e = typeof k === 'string' ? _keyExp(k) : [...k];
        if (e.length !== this.ngens()) {
          throw new ValueError('exponent vector has the wrong length');
        }
        terms.push([e, this._base_ring.__call__(c)]);
      }
      return new MPowerSeries<T>(this, terms, prec);
    }
    // a scalar
    const c = this._base_ring.__call__(x);
    return new MPowerSeries<T>(this, [[new Array(this.ngens()).fill(0), c]], prec);
  }

  /**
   * @see Reference: sage/rings/multi_power_series_ring.py:_repr_
   */
  toString(): string {
    return `Multivariate Power Series Ring in ${this._names.join(', ')} over ${this._base_ring}`;
  }
}

/**
 * An element of a multivariate power series ring.
 *
 * @see Reference: sage/rings/multi_power_series_ring_element.py:202 (MPowerSeries)
 */
export class MPowerSeries<T extends RingElement = RingElement> {
  private readonly _parent: MPowerSeriesRing<T>;
  /** exponent-vector key -> coefficient; only terms of total degree < prec. */
  private readonly _terms: Map<string, T>;
  private readonly _prec: number;

  constructor(
    parent: MPowerSeriesRing<T>,
    terms: Iterable<[number[] | string, T]>,
    prec: number = Number.POSITIVE_INFINITY
  ) {
    this._parent = parent;
    this._prec = prec;
    this._terms = new Map<string, T>();
    for (const [k, c] of terms) {
      if (c.isZero()) continue;
      const e = typeof k === 'string' ? _keyExp(k) : k;
      if (_totalDegree(e) >= prec) continue; // background truncation
      const key = _expKey(e);
      const cur = this._terms.get(key);
      const val = cur === undefined ? c : (cur.add(c) as T);
      if (val.isZero()) {
        this._terms.delete(key);
      } else {
        this._terms.set(key, val);
      }
    }
  }

  parent(): MPowerSeriesRing<T> {
    return this._parent;
  }

  base_ring(): CoefficientRing<T> {
    return this._parent.base_ring();
  }

  /**
   * Return the precision of `self` (a bound on the total degree).
   * @see Reference: sage/rings/multi_power_series_ring_element.py:1337 (prec)
   */
  prec(): number {
    return this._prec;
  }

  /** @see Reference: sage/rings/multi_power_series_ring_element.py:1337 (prec) */
  precision_absolute(): number {
    return this._prec;
  }

  /**
   * Return the dictionary with keys the exponents and values the coefficients.
   * @see Reference: sage/rings/multi_power_series_ring_element.py:1129 (monomial_coefficients)
   */
  monomial_coefficients(): [number[], T][] {
    return [...this._terms].map(([k, c]) => [_keyExp(k), c] as [number[], T]);
  }

  /**
   * Return the coefficient of the monomial `x1^e1 ... xk^ek`.
   * @see Reference: sage/rings/multi_power_series_ring_element.py:676 (__getitem__)
   */
  __getitem__(e: number[]): T {
    if (_totalDegree(e) >= this._prec) {
      throw new RangeError(
        'Cannot return the coefficients of terms of total degree greater than or equal to precision of self.'
      );
    }
    return this._terms.get(_expKey(e)) ?? this._parent.base_ring().zero();
  }

  /**
   * Return the valuation of `self`, i.e. the smallest total degree of a nonzero
   * term (the precision if there is none).
   * @see Reference: sage/rings/multi_power_series_ring_element.py:1432 (valuation)
   */
  valuation(): number {
    let v = Number.POSITIVE_INFINITY;
    for (const k of this._terms.keys()) {
      v = Math.min(v, _totalDegree(_keyExp(k)));
    }
    return this._terms.size === 0 ? this._prec : v;
  }

  /**
   * Return the (total) degree of the underlying polynomial.
   * @see Reference: sage/rings/multi_power_series_ring_element.py:1520 (degree)
   */
  degree(): number {
    let d = -1;
    for (const k of this._terms.keys()) {
      d = Math.max(d, _totalDegree(_keyExp(k)));
    }
    return d;
  }

  is_zero(): boolean {
    return this._terms.size === 0;
  }

  /**
   * A multivariate power series is a unit if and only if its constant
   * coefficient is a unit.
   * @see Reference: sage/rings/multi_power_series_ring_element.py:1537 (is_unit)
   */
  is_unit(): boolean {
    if (this.precision_absolute() === 0) {
      return false;
    }
    const c = this._terms.get(_expKey(new Array(this._parent.ngens()).fill(0)));
    if (c === undefined) {
      return false;
    }
    if (c.isUnit) {
      return c.isUnit();
    }
    if (this._parent.base_ring().is_field?.()) {
      return !c.isZero();
    }
    return c.eq(1) || c.eq(-1);
  }

  private _checkParent(other: MPowerSeries<T>): void {
    if (!this._parent.is_identical_to(other.parent())) {
      throw new TypeError('the two power series have different parents');
    }
  }

  /**
   * @see Reference: sage/rings/multi_power_series_ring_element.py:781 (_add_)
   */
  add(right: MPowerSeries<T>): MPowerSeries<T> {
    this._checkParent(right);
    return new MPowerSeries<T>(
      this._parent,
      [...this.monomial_coefficients(), ...right.monomial_coefficients()],
      Math.min(this._prec, right._prec)
    );
  }

  /**
   * @see Reference: sage/rings/multi_power_series_ring_element.py:801 (_sub_)
   */
  sub(right: MPowerSeries<T>): MPowerSeries<T> {
    return this.add(right.neg());
  }

  neg(): MPowerSeries<T> {
    return new MPowerSeries<T>(
      this._parent,
      this.monomial_coefficients().map(([e, c]) => [e, c.neg() as T]),
      this._prec
    );
  }

  /**
   * @see Reference: sage/rings/multi_power_series_ring_element.py:821 (_mul_)
   */
  mul(right: MPowerSeries<T>): MPowerSeries<T> {
    this._checkParent(right);
    // Precision of a product in the background univariate ring:
    // min(prec1 + val2, prec2 + val1).
    const v1 = this.valuation();
    const v2 = right.valuation();
    let prec: number;
    if (this._prec === Number.POSITIVE_INFINITY && right._prec === Number.POSITIVE_INFINITY) {
      prec = Number.POSITIVE_INFINITY;
    } else if (this._prec === Number.POSITIVE_INFINITY) {
      prec = right._prec + v1;
    } else if (right._prec === Number.POSITIVE_INFINITY) {
      prec = this._prec + v2;
    } else {
      prec = Math.min(this._prec + v2, right._prec + v1);
    }

    const terms: [number[], T][] = [];
    for (const [e1, c1] of this.monomial_coefficients()) {
      for (const [e2, c2] of right.monomial_coefficients()) {
        const e = e1.map((x, i) => x + e2[i]!);
        if (_totalDegree(e) >= prec) continue;
        terms.push([e, c1.mul(c2) as T]);
      }
    }
    return new MPowerSeries<T>(this._parent, terms, prec);
  }

  /**
   * Multiply by an element of the base ring.
   * @see Reference: sage/rings/multi_power_series_ring_element.py:846 (_lmul_)
   */
  scalar_mul(c: T): MPowerSeries<T> {
    return new MPowerSeries<T>(
      this._parent,
      this.monomial_coefficients().map(([e, a]) => [e, a.mul(c) as T]),
      this._prec
    );
  }

  /** `self^n` for a nonnegative integer `n` (binary powering). */
  pow(n: number | bigint): MPowerSeries<T> {
    let e = typeof n === 'bigint' ? n : BigInt(n);
    if (e < 0n) {
      return this.inv().pow(-e);
    }
    let result = this._parent.one();
    let base: MPowerSeries<T> = this;
    while (e > 0n) {
      if (e & 1n) result = result.mul(base);
      base = base.mul(base);
      e >>= 1n;
    }
    return result;
  }

  /**
   * Return the multiplicative inverse of this multivariate power series.
   *
   * Currently implemented only if the constant coefficient is a unit, exactly
   * as in SageMath.
   *
   * @see Reference: sage/rings/multi_power_series_ring_element.py:725 (__invert__)
   */
  inv(): MPowerSeries<T> {
    if (this.valuation() !== 0) {
      throw new NotImplementedError(
        'Multiplicative inverse of multivariate power series currently implemented only if constant coefficient is a unit.'
      );
    }
    const R = this._parent.base_ring();
    const zeroExp = new Array(this._parent.ngens()).fill(0);
    const c = this._terms.get(_expKey(zeroExp))!;
    const cinv = (c.inv ? c.inv() : R.one().div(c)) as T;

    // Precision of the inverse in the background univariate ring: the series
    // precision, or the ring's default precision for an exact input.
    const prec = this._prec === Number.POSITIVE_INFINITY ? this._parent.default_prec() : this._prec;

    // self = c*(1 + z) with z of positive valuation, so 1/self = c^-1 sum (-z)^m.
    const one = this._parent.one().add_bigoh(prec);
    const z = this.add_bigoh(prec).scalar_mul(cinv).sub(one);
    let result = one;
    let zp = one;
    for (let m = 1; m < prec; m++) {
      zp = zp.mul(z.neg());
      if (zp.is_zero()) break;
      result = result.add(zp);
    }
    return result.scalar_mul(cinv).add_bigoh(prec);
  }

  /**
   * Division in the ring of power series.
   *
   * SageMath falls back on `quo_rem` (multivariate polynomial division, which
   * needs Singular) when the denominator is not a unit; that fallback is not
   * ported.
   *
   * @see Reference: sage/rings/multi_power_series_ring_element.py:1059 (_div_)
   */
  div(denom: MPowerSeries<T>): MPowerSeries<T> {
    this._checkParent(denom);
    if (denom.is_zero()) {
      throw new ZeroDivisionError('');
    }
    if (denom.is_unit()) {
      return this.mul(denom.inv());
    }
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: division of multivariate power series by a non-unit (SageMath uses MPowerSeries.quo_rem)'
    );
  }

  /**
   * Return a multivariate power series of precision `prec` obtained by
   * truncating `self` at precision `prec`.
   * @see Reference: sage/rings/multi_power_series_ring_element.py:1353 (add_bigoh)
   */
  add_bigoh(prec: number): MPowerSeries<T> {
    return new MPowerSeries<T>(
      this._parent,
      this.monomial_coefficients(),
      Math.min(this._prec, prec)
    );
  }

  /** Alias of {@link MPowerSeries.add_bigoh}. */
  O(prec: number): MPowerSeries<T> {
    return this.add_bigoh(prec);
  }

  /**
   * Return the infinite precision multivariate power series formed by
   * truncating `self` at precision `prec`.
   * @see Reference: sage/rings/multi_power_series_ring_element.py:1401 (truncate)
   */
  truncate(prec: number = Number.POSITIVE_INFINITY): MPowerSeries<T> {
    return new MPowerSeries<T>(
      this._parent,
      this.add_bigoh(prec).monomial_coefficients(),
      Number.POSITIVE_INFINITY
    );
  }

  /**
   * Compare `self` to `other`: the two series are equal when they agree in
   * every total degree below the smaller of the two precisions.
   * @see Reference: sage/rings/multi_power_series_ring_element.py:746 (_richcmp_)
   */
  eq(other: MPowerSeries<T>): boolean {
    const prec = Math.min(this._prec, other._prec);
    const a = this.add_bigoh(prec);
    const b = other.add_bigoh(prec);
    if (a._terms.size !== b._terms.size) {
      return false;
    }
    for (const [k, c] of a._terms) {
      const d = b._terms.get(k);
      if (d === undefined || !c.eq(d)) {
        return false;
      }
    }
    return true;
  }

  /**
   * Evaluate `self`.
   *
   * If every argument lives in the parent of `self`, SageMath requires the
   * arguments to have positive valuation (unless `self` has infinite
   * precision) and gives the answer the precision
   * `self.prec() * min(valuations)`; otherwise it falls back on the formal
   * substitution {@link MPowerSeries._subs_formal}.
   *
   * @see Reference: sage/rings/multi_power_series_ring_element.py:442 (__call__)
   */
  __call__(...x: MPowerSeries<T>[]): MPowerSeries<T> {
    if (x.length === 1 && Array.isArray(x[0])) {
      x = x[0] as unknown as MPowerSeries<T>[];
    }
    if (x.length !== this._parent.ngens()) {
      throw new ValueError('Number of arguments does not match number of variables in parent.');
    }
    if (!x.every((xi) => this._parent.is_identical_to(xi.parent()))) {
      // Input does not coerce to parent ring of self: attempt formal substitution
      return this._subs_formal(...x);
    }

    const args: MPowerSeries<T>[] = [];
    const valn_list: number[] = [];
    for (const xi of x) {
      const v = xi.valuation();
      if (v === 0 && this._prec !== Number.POSITIVE_INFINITY) {
        throw new TypeError(
          'Substitution defined only for elements of positive valuation, unless self has infinite precision.'
        );
      }
      if (v > 0) {
        args.push(xi.add_bigoh(v * this._prec));
        valn_list.push(v);
      } else {
        args.push(xi);
      }
    }
    const newprec =
      this._prec === Number.POSITIVE_INFINITY
        ? Number.POSITIVE_INFINITY
        : this._prec * Math.min(...valn_list);
    // SageMath substitutes into ``self._value()``, the exact polynomial part.
    return this.truncate()
      ._subs_formal(...args)
      .add_bigoh(newprec);
  }

  /**
   * Substitution of the inputs as variables of `self`.
   *
   * The inputs need not be elements of the same ring as `self`; the result is
   * `sum c * prod x_i^{m_i}` with `O(...)^{self.prec()}` added at the end.
   *
   * @see Reference: sage/rings/multi_power_series_ring_element.py:512 (_subs_formal)
   */
  _subs_formal(...x: MPowerSeries<T>[]): MPowerSeries<T> {
    if (x.length === 1 && Array.isArray(x[0])) {
      x = x[0] as unknown as MPowerSeries<T>[];
    }
    const n = this._parent.ngens();
    if (x.length !== n) {
      throw new ValueError('Input must be of correct length.');
    }
    if (n === 0) {
      return this;
    }
    const P = x[0]!.parent();
    let y = P.zero();
    // Cache the powers of each argument.
    const powers: MPowerSeries<T>[][] = x.map((xi) => [P.one(), xi]);
    for (const [m, c] of this.monomial_coefficients()) {
      let term = P.__call__(c);
      for (let i = 0; i < n; i++) {
        if (m[i] === 0) continue;
        const pi = powers[i]!;
        while (pi.length <= m[i]!) {
          pi.push(pi[pi.length - 1]!.mul(x[i]!));
        }
        term = term.mul(pi[m[i]!]!);
      }
      y = y.add(term);
    }
    if (this._prec === Number.POSITIVE_INFINITY) {
      return y;
    }
    return y.add_bigoh(this._prec);
  }

  /**
   * String representation.
   *
   * The monomials are printed in SageMath's `negdeglex` order (the term order
   * of the foreground polynomial ring of a multivariate power series ring):
   * increasing total degree, then lexicographically decreasing.
   *
   * @see Reference: sage/rings/multi_power_series_ring_element.py:597 (_repr_)
   * @see Reference: sage/rings/multi_power_series_ring.py:311 (order='negdeglex')
   */
  toString(): string {
    const names = this._parent.variable_names();
    const keys = [...this._terms.keys()].sort((k1, k2) => {
      const e1 = _keyExp(k1);
      const e2 = _keyExp(k2);
      const d1 = _totalDegree(e1);
      const d2 = _totalDegree(e2);
      if (d1 !== d2) return d1 - d2;
      for (let i = 0; i < e1.length; i++) {
        if (e1[i] !== e2[i]) return e2[i]! - e1[i]!;
      }
      return 0;
    });
    const parts: string[] = [];
    for (const key of keys) {
      const e = _keyExp(key);
      const c = this._terms.get(key)!;
      const mon = e
        .map((k, i) => (k === 0 ? '' : k === 1 ? names[i]! : `${names[i]}^${k}`))
        .filter((s) => s !== '')
        .join('*');
      const cs = c.toString();
      if (mon === '') {
        parts.push(cs);
      } else if (cs === '1') {
        parts.push(mon);
      } else if (cs === '-1') {
        parts.push(`-${mon}`);
      } else {
        parts.push(`${cs}*${mon}`);
      }
    }
    const value = parts.length === 0 ? '0' : parts.join(' + ').replaceAll('+ -', '- ');
    if (this._prec === Number.POSITIVE_INFINITY) {
      return value;
    }
    return `${value} + O(${names.join(', ')})^${this._prec}`;
  }
}
