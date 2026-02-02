/**
 * @module sage/rings/power_series_ring
 * @description Power series rings
 *
 * Port of: sage/rings/power_series_ring.py
 * Reference: reference/sage/src/sage/rings/power_series_ring.py
 */

import { ArithmeticError, NotImplementedError, ValueError, ZeroDivisionError } from '../errors.js';

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
    // Strip trailing zeros
    let lastNonZero = -1;
    for (let i = coefficients.length - 1; i >= 0; i--) {
      if (!coefficients[i]!.isZero()) {
        lastNonZero = i;
        break;
      }
    }
    this._coefficients = lastNonZero < 0 ? [] : coefficients.slice(0, lastNonZero + 1);
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
    if (this._coefficients.length === 0) {
      return Number.POSITIVE_INFINITY;
    }
    for (let i = 0; i < this._coefficients.length; i++) {
      if (!this._coefficients[i]!.isZero()) {
        return i;
      }
    }
    return Number.POSITIVE_INFINITY;
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
    const computePrec = Math.min(targetPrec, this._prec);

    // Check that constant term is 1
    const c0 =
      this._coefficients.length > 0 ? this._coefficients[0]! : this._parent.base_ring().zero();
    if (c0.isZero() || !c0.eq(1)) {
      throw new ArithmeticError('constant term of power series is not 1');
    }

    // Use: log(1 + g) = g - g^2/2 + g^3/3 - g^4/4 + ...
    // where g = f - 1
    const baseRing = this._parent.base_ring();
    const one = baseRing.one();
    const g = this.sub(this._parent.one()).add_bigoh(computePrec);

    const result: T[] = [baseRing.zero()]; // constant term is 0

    // We use the recurrence for log coefficients
    // If f = 1 + sum_{n>=1} a_n x^n and log(f) = sum_{n>=1} b_n x^n
    // Then n*b_n = a_n - sum_{k=1}^{n-1} k*b_k * a_{n-k}
    // Note: g.list() returns coefficients starting from x^0, so a_n = aCoeffs[n]
    const aCoeffs = g.list();
    const bCoeffs: T[] = [];

    for (let n = 1; n < computePrec; n++) {
      // a_n is the coefficient of x^n in g, which is aCoeffs[n]
      const an = n < aCoeffs.length ? (aCoeffs[n] ?? baseRing.zero()) : baseRing.zero();
      let sum = an;
      for (let k = 1; k < n; k++) {
        const bk = bCoeffs[k - 1] ?? baseRing.zero();
        // a_{n-k} is the coefficient of x^{n-k} in g, which is aCoeffs[n-k]
        const an_k =
          n - k < aCoeffs.length ? (aCoeffs[n - k] ?? baseRing.zero()) : baseRing.zero();
        // k * b_k * a_{n-k}
        let kBk = bk;
        for (let j = 1; j < k; j++) {
          kBk = kBk.add(bk) as T;
        }
        sum = sum.sub(kBk.mul(an_k)) as T;
      }
      // Divide by n
      let divisor = one;
      for (let j = 1; j < n; j++) {
        divisor = divisor.add(one) as T;
      }
      bCoeffs.push(sum.div(divisor) as T);
      result.push(bCoeffs[bCoeffs.length - 1]!);
    }

    return new PowerSeriesElement<T>(this._parent, result, computePrec);
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
    const two = one.add(one) as T;
    const half = one.div(two) as T;

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
      bCoeffs.push(sum.mul(half) as T);
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

    return finalResult;
  }

  /**
   * Return the n-th root.
   * @see Reference: sage/rings/power_series_ring_element.pyx:nth_root
   */
  nth_root(n: number, prec?: number): PowerSeriesElement<T> {
    // Reference: sage/rings/power_series_ring_element.pyx:nth_root
    // Compute f^(1/n) using the formula f^(1/n) = c0^(1/n) * (1 + g)^(1/n)
    // where g = f/c0 - 1

    if (n === 0) {
      throw new ValueError('n must be nonzero');
    }

    if (n < 0) {
      // f^(1/n) = (f^(-1))^(-1/n)
      return this.inverse().nth_root(-n, prec);
    }

    if (n === 1) {
      return this;
    }

    if (n === 2) {
      return this.sqrt(prec);
    }

    const targetPrec = prec ?? this._parent.default_prec();
    const computePrec = Math.min(targetPrec, this._prec);

    if (this.is_zero()) {
      const newPrec =
        this._prec === Number.POSITIVE_INFINITY
          ? Number.POSITIVE_INFINITY
          : Math.floor(this._prec / n);
      return new PowerSeriesElement<T>(this._parent, [], newPrec);
    }

    const val = this.valuation();
    if (val === Number.POSITIVE_INFINITY) {
      return new PowerSeriesElement<T>(this._parent, [], Number.POSITIVE_INFINITY);
    }

    if (val % n !== 0) {
      throw new ValueError(
        `power series does not have an ${n}th root since valuation is not divisible by ${n}`
      );
    }

    const baseRing = this._parent.base_ring();
    const one = baseRing.one();

    // Get the valuation zero part
    const valuationZeroPart = val > 0 ? this._shiftRight(val) : this;
    const c0 = valuationZeroPart.__getitem__(0);

    // For n-th root of c0, we need base ring support
    // For now, only handle c0 = 1
    if (!c0.eq(1)) {
      throw new NotImplementedError(`nth_root for constant term != 1 (got ${c0})`);
    }

    // Compute (1 + g)^(1/n) using the binomial series
    // (1 + g)^(1/n) = sum_{k>=0} C(1/n, k) * g^k
    // where C(1/n, k) = (1/n)(1/n-1)...(1/n-k+1) / k!
    const g = valuationZeroPart.sub(this._parent.one()).add_bigoh(computePrec);
    const gCoeffs = g.list();

    // Build coefficients iteratively
    const result: T[] = [one];

    // For each coefficient, we compute using the generalized binomial theorem
    // (1+g)^(1/n) = 1 + (1/n)g + (1/n)(1/n-1)/2! g^2 + ...
    // Coefficient of x^m in (1+g)^(1/n) where g = a_1 x + a_2 x^2 + ...
    // is computed by tracking powers of g

    // Use the recurrence: if y = (1+g)^(1/n), then ny' = y'g + (1+g)^((1-n)/n) * g'
    // More simply: n * y^(n-1) * y' = g' * (1+g)^((1-n)/n) => n y' = g'/(1+g) * y
    // Or use: n * d/dx[(1+g)^(1/n)] = (1+g)^((1-n)/n) * g'

    // Simpler: directly compute using convolution
    // If y = sum b_k x^k and g = sum a_k x^k (a_0 = 0)
    // Then y^n = 1 + g, which gives a recurrence
    // But this is complex. Instead use exp/log approach if available.

    // For now, use Newton iteration: y_{k+1} = y_k - (y_k^n - f) / (n * y_k^(n-1))
    // = y_k * (1 + (1 - f/y_k^n)/(n)) (approximately)
    // = y_k * ((n-1) + f/y_k^n) / n

    let y = this._parent.one().add_bigoh(computePrec);
    const nVal = baseRing.__call__(BigInt(n));

    // Newton iterations
    for (let iter = 0; iter < Math.ceil(Math.log2(computePrec)) + 2; iter++) {
      // y_new = y * ((n-1) + f * y^(-n)) / n
      // = y * ((n-1)*y^n + f) / (n * y^n)
      // = ((n-1)*y + f*y^(1-n)) / n
      const yPowN = y._pow(n);
      const numerator = y._scalarMul(nVal.sub(one) as T).add(valuationZeroPart.div(yPowN));
      y = numerator._scalarDiv(nVal).add_bigoh(computePrec);
    }

    // Shift back if original had valuation > 0
    if (val > 0) {
      y = y._shiftLeft(val / n);
    }

    return y;
  }

  /**
   * Return the composition f(g(x)) where f = this.
   * @see Reference: sage/rings/power_series_ring_element.pyx:__call__
   */
  __call__(g: PowerSeriesElement<T>): PowerSeriesElement<T> {
    // g must have positive valuation for composition to converge
    if (g.valuation() <= 0) {
      throw new ValueError('can only compose with series of positive valuation');
    }

    const computePrec = Math.min(this._prec, g.prec(), this._parent.default_prec());

    // Compute f(g) = sum_{n>=0} f[n] * g^n
    let result = this._parent.zero().add_bigoh(computePrec);
    let gPower = this._parent.one().add_bigoh(computePrec);

    for (let n = 0; n < this._coefficients.length && n < computePrec; n++) {
      const coeff = this._coefficients[n]!;
      if (!coeff.isZero()) {
        // result += coeff * g^n
        result = result.add(gPower._scalarMul(coeff)).add_bigoh(computePrec);
      }
      gPower = gPower.mul(g).add_bigoh(computePrec);
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
  pade(m: number, n: number): [PowerSeriesElement<T>, PowerSeriesElement<T>] {
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

    // Create power series for Q and P
    const Q = new PowerSeriesElement<T>(
      this._parent,
      r1.length === 0 ? [zero] : r1,
      Number.POSITIVE_INFINITY
    );
    const P = new PowerSeriesElement<T>(
      this._parent,
      t1.length === 0 ? [one] : t1,
      Number.POSITIVE_INFINITY
    );

    return [Q, P];
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
    if (!this.is_unit()) {
      throw new ZeroDivisionError('Power series is not invertible (constant term is not a unit)');
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
   */
  _shiftRight(n: number): PowerSeriesElement<T> {
    if (n <= 0) return this;
    const newCoeffs = this._coefficients.slice(n);
    const newPrec =
      this._prec === Number.POSITIVE_INFINITY
        ? Number.POSITIVE_INFINITY
        : Math.max(0, this._prec - n);
    return new PowerSeriesElement<T>(this._parent, newCoeffs, newPrec);
  }

  /**
   * Shift left (multiply by x^n).
   */
  _shiftLeft(n: number): PowerSeriesElement<T> {
    if (n <= 0) return this;
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
 * A ring of Laurent series over a base ring.
 * @see Reference: sage/rings/laurent_series_ring.py:LaurentSeriesRing
 */
export class LaurentSeriesRing<T extends RingElement = RingElement> {
  private readonly _base_ring: CoefficientRing<T>;
  private readonly _name: string;
  private readonly _default_prec: number;
  private readonly _power_series_ring: PowerSeriesRing<T>;

  constructor(base_ring: CoefficientRing<T>, name: string = 'x', default_prec: number = 20) {
    this._base_ring = base_ring;
    this._name = name;
    this._default_prec = default_prec;
    this._power_series_ring = new PowerSeriesRing<T>(base_ring, name, default_prec);
  }

  /**
   * Return the base ring.
   * @see Reference: sage/rings/laurent_series_ring.py:base_ring
   */
  base_ring(): CoefficientRing<T> {
    return this._base_ring;
  }

  /**
   * Return the variable name.
   * @see Reference: sage/rings/laurent_series_ring.py:variable_name
   */
  variable_name(): string {
    return this._name;
  }

  /**
   * Return the power series ring.
   * @see Reference: sage/rings/laurent_series_ring.py:power_series_ring
   */
  power_series_ring(): PowerSeriesRing<T> {
    return this._power_series_ring;
  }

  /**
   * Return the generator.
   * @see Reference: sage/rings/laurent_series_ring.py:gen
   */
  gen(): LaurentSeriesElement<T> {
    return new LaurentSeriesElement<T>(this, this._power_series_ring.gen(), 0);
  }

  /**
   * Coerce an element to this ring.
   * @see Reference: sage/rings/laurent_series_ring.py:__call__
   */
  __call__(f: unknown, prec?: number): LaurentSeriesElement<T> {
    if (f instanceof LaurentSeriesElement) {
      return f as LaurentSeriesElement<T>;
    }
    if (f instanceof PowerSeriesElement) {
      return new LaurentSeriesElement<T>(this, f as PowerSeriesElement<T>, 0);
    }
    const ps = this._power_series_ring.__call__(f, prec);
    return new LaurentSeriesElement<T>(this, ps, 0);
  }

  toString(): string {
    return `Laurent Series Ring in ${this._name} over ${this._base_ring}`;
  }
}

/**
 * An element of a Laurent series ring.
 * @see Reference: sage/rings/laurent_series_ring_element.pyx:LaurentSeries
 */
export class LaurentSeriesElement<T extends RingElement = RingElement> {
  private readonly _parent: LaurentSeriesRing<T>;
  private readonly _power_series: PowerSeriesElement<T>;
  private readonly _valuation_shift: number;

  constructor(
    parent: LaurentSeriesRing<T>,
    power_series: PowerSeriesElement<T>,
    valuation_shift: number
  ) {
    this._parent = parent;
    this._power_series = power_series;
    this._valuation_shift = valuation_shift;
  }

  /**
   * Return the parent ring.
   * @see Reference: sage/rings/laurent_series_ring_element.pyx:parent
   */
  parent(): LaurentSeriesRing<T> {
    return this._parent;
  }

  /**
   * Return the valuation (can be negative for Laurent series).
   * @see Reference: sage/rings/laurent_series_ring_element.pyx:valuation
   */
  valuation(): number {
    const psVal = this._power_series.valuation();
    if (psVal === Number.POSITIVE_INFINITY) {
      return Number.POSITIVE_INFINITY;
    }
    return psVal + this._valuation_shift;
  }

  /**
   * Return the power series part (positive powers).
   * @see Reference: sage/rings/laurent_series_ring_element.pyx:power_series
   */
  power_series(): PowerSeriesElement<T> {
    if (this._valuation_shift >= 0) {
      return this._power_series;
    }
    throw new Error('Laurent series has negative valuation; cannot convert to power series');
  }

  /**
   * Return the principal part (negative powers).
   * The principal part is the sum of terms with negative exponents.
   * @see Reference: sage/rings/laurent_series_ring_element.pyx:principal_part
   */
  principal_part(): LaurentSeriesElement<T> {
    // The principal part consists of terms x^k where k < 0
    // In our representation, these are coefficients from index 0 to |valuation_shift| - 1
    // in the underlying power series, when valuation_shift < 0

    if (this._valuation_shift >= 0) {
      // No negative powers
      return new LaurentSeriesElement<T>(this._parent, this._parent.power_series_ring().zero(), 0);
    }

    // Extract coefficients for negative powers
    // valuation_shift is negative, so we have terms x^{valuation_shift}, ..., x^{-1}
    const numNegTerms = -this._valuation_shift;
    const baseRing = this._parent.base_ring();
    const coeffs: T[] = [];

    for (let i = 0; i < numNegTerms; i++) {
      coeffs.push(this._power_series.__getitem__(i));
    }

    const principalPS = new PowerSeriesElement<T>(
      this._parent.power_series_ring(),
      coeffs,
      Number.POSITIVE_INFINITY
    );

    return new LaurentSeriesElement<T>(this._parent, principalPS, this._valuation_shift);
  }

  /**
   * Return the residue (coefficient of x^-1).
   * @see Reference: sage/rings/laurent_series_ring_element.pyx:residue
   */
  residue(): T {
    // Coefficient of x^{-1} is the coefficient at index -1 - valuation_shift in the power series
    const idx = -1 - this._valuation_shift;
    if (idx < 0) {
      return this._parent.base_ring().zero();
    }
    return this._power_series.__getitem__(idx);
  }

  toString(): string {
    return `[LaurentSeriesElement with shift ${this._valuation_shift}]`;
  }
}
