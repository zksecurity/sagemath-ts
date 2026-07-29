/**
 * @module sage/rings/polynomial/polynomial_element
 * @description Polynomial elements over arbitrary coefficient rings
 *
 * Port of: sage/rings/polynomial/polynomial_element.pyx
 */

import {
  factor as factorInteger,
  gcd as gcdBigInt,
  is_prime,
  next_prime,
} from '../../arith/misc.js';
import {
  ArithmeticError,
  AssertionError,
  NotImplementedError,
  ValueError,
  ZeroDivisionError,
} from '../../errors.js';
import { current_randstate } from '../../misc/randstate.js';

/**
 * Interface for coefficient rings/fields.
 */
export interface CoefficientRing<T> {
  zero(): T;
  one(): T;
  __call__(x: unknown): T;
  is_field?(): boolean;
}

/**
 * Interface for ring elements that can be used as polynomial coefficients.
 */
export interface RingElement {
  add(other: this): this;
  sub(other: this): this;
  mul(other: this): this;
  neg(): this;
  eq(other: this | number): boolean;
  isZero(): boolean;
  toString(): string;
}

/**
 * A polynomial with coefficients in a ring R.
 *
 * Internally stored as an array of coefficients where coeffs[i] is the
 * coefficient of x^i. Trailing zeros are removed.
 */
export class Polynomial<C extends RingElement> {
  readonly coeffs: readonly C[];
  readonly parent: PolynomialRingBase<C>;

  constructor(coeffs: C[], parent: PolynomialRingBase<C>) {
    this.parent = parent;

    // Remove trailing zeros
    let len = coeffs.length;
    while (len > 0 && coeffs[len - 1]!.isZero()) {
      len--;
    }
    this.coeffs = coeffs.slice(0, len);
  }

  /**
   * Return the degree of this polynomial.
   * The zero polynomial has degree -1.
   */
  degree(): number {
    return this.coeffs.length - 1;
  }

  /**
   * Return the leading coefficient.
   */
  leading_coefficient(): C {
    if (this.coeffs.length === 0) {
      return this.parent.base_ring.zero() as C;
    }
    return this.coeffs[this.coeffs.length - 1]!;
  }

  /**
   * Return the coefficient of x^n.
   */
  getCoeff(n: number): C {
    if (n < 0 || n >= this.coeffs.length) {
      return this.parent.base_ring.zero() as C;
    }
    return this.coeffs[n]!;
  }

  /**
   * Check if this is the zero polynomial.
   */
  isZero(): boolean {
    return this.coeffs.length === 0;
  }

  /**
   * Check if this is a constant polynomial.
   */
  isConstant(): boolean {
    return this.coeffs.length <= 1;
  }

  /**
   * Check if this is monic (leading coefficient is 1).
   */
  is_monic(): boolean {
    if (this.coeffs.length === 0) {
      return false;
    }
    return this.leading_coefficient().eq(1);
  }

  /**
   * Add two polynomials.
   */
  add(other: Polynomial<C>): Polynomial<C> {
    const maxLen = Math.max(this.coeffs.length, other.coeffs.length);
    const result: C[] = [];

    for (let i = 0; i < maxLen; i++) {
      const a = this.getCoeff(i);
      const b = other.getCoeff(i);
      result.push(a.add(b) as C);
    }

    return new Polynomial(result, this.parent);
  }

  /**
   * Subtract two polynomials.
   */
  sub(other: Polynomial<C>): Polynomial<C> {
    const maxLen = Math.max(this.coeffs.length, other.coeffs.length);
    const result: C[] = [];

    for (let i = 0; i < maxLen; i++) {
      const a = this.getCoeff(i);
      const b = other.getCoeff(i);
      result.push(a.sub(b) as C);
    }

    return new Polynomial(result, this.parent);
  }

  /**
   * Negate the polynomial.
   */
  neg(): Polynomial<C> {
    return new Polynomial(
      this.coeffs.map((c) => c.neg() as C),
      this.parent
    );
  }

  /**
   * Multiply two polynomials.
   */
  mul(other: Polynomial<C>): Polynomial<C> {
    if (this.isZero() || other.isZero()) {
      return this.parent.zero();
    }

    const resultLen = this.coeffs.length + other.coeffs.length - 1;
    const result: C[] = [];

    // Initialize with zeros
    for (let i = 0; i < resultLen; i++) {
      result.push(this.parent.base_ring.zero() as C);
    }

    // Schoolbook multiplication
    for (let i = 0; i < this.coeffs.length; i++) {
      for (let j = 0; j < other.coeffs.length; j++) {
        const prod = this.coeffs[i]!.mul(other.coeffs[j]!) as C;
        result[i + j] = result[i + j]!.add(prod) as C;
      }
    }

    return new Polynomial(result, this.parent);
  }

  /**
   * Multiply by a scalar.
   */
  scalar_mul(c: C): Polynomial<C> {
    if (c.isZero()) {
      return this.parent.zero();
    }
    return new Polynomial(
      this.coeffs.map((coeff) => coeff.mul(c) as C),
      this.parent
    );
  }

  /**
   * Compute the remainder of this polynomial divided by other.
   * Only works over fields.
   */
  mod(other: Polynomial<C>): Polynomial<C> {
    const [_q, r] = this.quo_rem(other);
    return r;
  }

  /**
   * Compute quotient and remainder of the Euclidean division.
   *
   * Raises a {@link ZeroDivisionError} if `other` is zero, and an
   * {@link ArithmeticError} if the division is not exact (i.e. a quotient
   * coefficient does not lie in the base ring).
   *
   * @see Reference: sage/rings/polynomial/polynomial_element.pyx:12548 (quo_rem)
   */
  quo_rem(other: Polynomial<C>): [Polynomial<C>, Polynomial<C>] {
    if (other.isZero()) {
      throw new ZeroDivisionError('division by zero polynomial');
    }

    if (this.degree() < other.degree()) {
      return [this.parent.zero(), this];
    }

    // Make a mutable copy of coefficients
    const remainder = [...this.coeffs] as C[];
    const divisorLC = other.leading_coefficient();
    const divisorDeg = other.degree();
    const quotientCoeffs: C[] = [];

    // Sage first tries ``inverse_of_unit()`` on the leading coefficient; when
    // that succeeds every quotient coefficient is automatically in the base
    // ring and no further check is needed.  Only in the fallback branch
    // ("convert") does it verify that the quotient coefficient lies in R.
    const lcInverse = inverseOfUnit(divisorLC, this.parent.base_ring);

    // Initialize quotient with zeros
    for (let i = 0; i <= this.degree() - divisorDeg; i++) {
      quotientCoeffs.push(this.parent.base_ring.zero() as C);
    }

    for (let i = this.degree(); i >= divisorDeg; i--) {
      if (remainder[i]?.isZero()) {
        continue;
      }

      // Compute quotient coefficient
      // This requires the coefficient ring to support division
      let qCoeff: C;
      if (lcInverse !== null) {
        qCoeff = remainder[i]!.mul(lcInverse) as C;
      } else {
        qCoeff = divideCoeffs(remainder[i]!, divisorLC);
        // Sage raises here when the quotient does not lie in the base ring
        // (`polynomial_element.pyx:12634-12640`); a coefficient ring whose
        // division truncates (e.g. ZZ) would otherwise silently return garbage.
        if (!qCoeff.mul(divisorLC).eq(remainder[i]!)) {
          throw new ArithmeticError(
            'division non exact (consider coercing to polynomials over the fraction field)'
          );
        }
      }
      quotientCoeffs[i - divisorDeg] = qCoeff;

      // Subtract qCoeff * other * x^(i - divisorDeg) from remainder
      for (let j = 0; j <= divisorDeg; j++) {
        const prod = qCoeff.mul(other.coeffs[j]!) as C;
        remainder[i - divisorDeg + j] = remainder[i - divisorDeg + j]!.sub(prod) as C;
      }
    }

    return [new Polynomial(quotientCoeffs, this.parent), new Polynomial(remainder, this.parent)];
  }

  /**
   * Compute the pseudo-division of two polynomials.
   *
   * Returns `[Q, R]` such that `l^(m-n+1) * self = Q*other + R` with
   * `deg(R) < deg(other)`, where `m = deg(self)`, `n = deg(other)` and `l` is
   * the leading coefficient of `other`.  Unlike {@link quo_rem} this needs no
   * division in the base ring.
   *
   * Algorithm 3.1.2 in [Coh1993].
   *
   * @see Reference: sage/rings/polynomial/polynomial_element.pyx:5375 (pseudo_quo_rem)
   */
  pseudo_quo_rem(other: Polynomial<C>): [Polynomial<C>, Polynomial<C>] {
    if (other.isZero()) {
      throw new ZeroDivisionError('Pseudo-division by zero is not possible');
    }

    // If other is a constant then R = 0 and Q = self * other^deg(self)
    if (other.degree() === 0) {
      const c = other.getCoeff(0);
      let scale = this.parent.base_ring.one() as C;
      for (let i = 0; i < this.degree(); i++) {
        scale = scale.mul(c) as C;
      }
      return [this.scalar_mul(scale), this.parent.zero()];
    }

    let R: Polynomial<C> = this;
    const B = other;
    let Q = this.parent.zero();
    let e = this.degree() - other.degree() + 1;
    const d = B.leading_coefficient();

    while (R.degree() >= B.degree() && !R.isZero()) {
      const c = R.leading_coefficient();
      const diffdeg = R.degree() - B.degree();
      Q = Q.scalar_mul(d).add(new Polynomial([c], this.parent).shift(diffdeg));
      R = R.scalar_mul(d).sub(B.scalar_mul(c).shift(diffdeg));
      e -= 1;
    }

    let q = this.parent.base_ring.one() as C;
    for (let i = 0; i < e; i++) {
      q = q.mul(d) as C;
    }

    return [Q.scalar_mul(q), R.scalar_mul(q)];
  }

  /**
   * Compute this^n.
   */
  pow(n: number | bigint): Polynomial<C> {
    let exp = typeof n === 'bigint' ? n : BigInt(n);

    if (exp < 0n) {
      throw new ValueError('negative exponent not supported for polynomials');
    }

    if (exp === 0n) {
      return this.parent.one();
    }

    // Binary exponentiation
    let result = this.parent.one();
    let base: Polynomial<C> = this;

    while (exp > 0n) {
      if ((exp & 1n) === 1n) {
        result = result.mul(base);
      }
      base = base.mul(base);
      exp >>= 1n;
    }

    return result;
  }

  /**
   * Evaluate the polynomial at a point.
   */
  evaluate(x: C): C {
    if (this.coeffs.length === 0) {
      return this.parent.base_ring.zero() as C;
    }

    // Horner's method
    let result = this.coeffs[this.coeffs.length - 1]!;
    for (let i = this.coeffs.length - 2; i >= 0; i--) {
      result = result.mul(x).add(this.coeffs[i]!) as C;
    }

    return result;
  }

  /**
   * Check equality.
   */
  eq(other: Polynomial<C>): boolean {
    if (this.coeffs.length !== other.coeffs.length) {
      return false;
    }
    for (let i = 0; i < this.coeffs.length; i++) {
      if (!this.coeffs[i]!.eq(other.coeffs[i]!)) {
        return false;
      }
    }
    return true;
  }

  /**
   * Return the derivative of this polynomial.
   *
   * @returns The formal derivative d/dx of this polynomial
   *
   * @example
   * ```typescript
   * // If f = x^3 + 2x + 1, then f.derivative() = 3x^2 + 2
   * ```
   *
   * @see Reference: sage/rings/polynomial/polynomial_element.pyx:derivative
   */
  derivative(): Polynomial<C> {
    if (this.coeffs.length <= 1) {
      return this.parent.zero();
    }

    const result: C[] = [];
    for (let i = 1; i < this.coeffs.length; i++) {
      // Multiply coefficient by i (the power).  Sage computes ``n * self[n]``
      // in the base ring (`polynomial_element.pyx:_derivative`); we use
      // double-and-add so that no coercion of the integer ``n`` is required,
      // which keeps this O(d log d) instead of O(d^2).
      result.push(mulByInteger(this.coeffs[i]!, i, this.parent.base_ring));
    }

    return new Polynomial(result, this.parent);
  }

  /**
   * Return the GCD of this polynomial and other.
   *
   * Uses the Euclidean algorithm. Only works over fields.
   *
   * @param other - Another polynomial in the same ring
   * @returns The monic GCD of this and other
   *
   * @see Reference: sage/rings/polynomial/polynomial_element.pyx:gcd
   */
  gcd(other: Polynomial<C>): Polynomial<C> {
    const baseRing = this.parent.base_ring;

    // Over ZZ, Sage delegates to FLINT's fmpz_poly_gcd (a subresultant PRS on
    // the primitive parts, times the gcd of the contents), which is *not* the
    // Euclidean algorithm: coefficient division is not exact in ZZ.
    if (isIntegerRing(baseRing)) {
      if (this.isZero() && other.isZero()) {
        return this.parent.zero();
      }
      const a = extractIntegerCoeffs(this);
      const b = extractIntegerCoeffs(other);
      const g = intPolyGcdWithContent(a, b);
      return new Polynomial(
        g.map((c) => baseRing.__call__(c) as C),
        this.parent
      );
    }

    if (!ringIsField(baseRing)) {
      throw new NotImplementedError(
        `${baseRing} does not provide a gcd implementation for univariate polynomials`
      );
    }

    // Fields: Euclidean algorithm (sage/categories/fields.py:_gcd_univariate_polynomial)
    let a: Polynomial<C> = this;
    let b: Polynomial<C> = other;

    while (!b.isZero()) {
      const [_q, r] = a.quo_rem(b);
      a = b;
      b = r;
    }

    // Return monic GCD (zero stays zero)
    return a.isZero() ? a : a._monic();
  }

  /**
   * Return the extended GCD of this polynomial and other.
   *
   * Returns (g, s, t) such that g = gcd(this, other) = s*this + t*other.
   *
   * @param other - Another polynomial in the same ring
   * @returns Tuple [g, s, t] where g = s*this + t*other
   *
   * @see Reference: sage/rings/polynomial/polynomial_element.pyx:xgcd
   */
  xgcd(other: Polynomial<C>): [Polynomial<C>, Polynomial<C>, Polynomial<C>] {
    const R = this.parent;
    const baseRing = R.base_ring;

    if (!ringIsField(baseRing)) {
      throw new NotImplementedError(
        `${baseRing} does not provide an xgcd implementation for univariate polynomials`
      );
    }

    const zero = R.zero();
    const one = R.one();

    // sage/categories/fields.py:526-543 (_xgcd_univariate_polynomial)
    if (other.isZero()) {
      if (this.isZero()) {
        return [zero, zero, zero];
      }
      const c = divideCoeffs(baseRing.one() as C, this.leading_coefficient());
      return [this.scalar_mul(c), R.__call__(c), zero];
    }
    if (this.isZero()) {
      const c = divideCoeffs(baseRing.one() as C, other.leading_coefficient());
      return [other.scalar_mul(c), zero, R.__call__(c)];
    }

    let u = one;
    let d: Polynomial<C> = this;
    let v1 = zero;
    let v3 = other;

    while (!v3.isZero()) {
      const [q, r] = d.quo_rem(v3);
      const newU = v1;
      const newD = v3;
      v1 = u.sub(v1.mul(q));
      v3 = r;
      u = newU;
      d = newD;
    }

    // v = (d - a*u) // b
    let v = d.sub(this.mul(u)).quo_rem(other)[0];

    if (!d.isZero()) {
      const c = divideCoeffs(baseRing.one() as C, d.leading_coefficient());
      d = d.scalar_mul(c);
      u = u.scalar_mul(c);
      v = v.scalar_mul(c);
    }

    return [d, u, v];
  }

  /**
   * Return the composition f(g) where this polynomial is f.
   *
   * @param other - The polynomial g to substitute for x
   * @returns f(g(x))
   *
   * @example
   * ```typescript
   * // If f = x^2 + 1 and g = x + 1, then f.compose(g) = (x+1)^2 + 1 = x^2 + 2x + 2
   * ```
   *
   * @see Reference: sage/rings/polynomial/polynomial_element.pyx:compose_trunc
   */
  compose(other: Polynomial<C>): Polynomial<C> {
    if (this.isZero()) {
      return this.parent.zero();
    }

    // Use Horner's method for polynomial composition
    let result = this.parent.__call__(this.coeffs[this.coeffs.length - 1]!);
    for (let i = this.coeffs.length - 2; i >= 0; i--) {
      result = result.mul(other).add(this.parent.__call__(this.coeffs[i]!));
    }

    return result;
  }

  /**
   * Return a monic version of this polynomial (leading coefficient = 1).
   *
   * A monic polynomial has leading coefficient 1. This method divides
   * the polynomial by its leading coefficient.
   *
   * @returns Monic polynomial equal to this / leading_coefficient
   *
   * @example
   * ```typescript
   * // 2x^2 + 4x + 2 becomes x^2 + 2x + 1
   * ```
   *
   * @see Reference: sage/rings/polynomial/polynomial_element.pyx:monic
   */
  monic(): Polynomial<C> {
    if (this.isZero()) {
      return this;
    }
    const lc = this.leading_coefficient();
    if (lc.eq(1)) {
      return this;
    }
    const lcInv = divideCoeffs(this.parent.base_ring.one() as C, lc);
    return this.scalar_mul(lcInv);
  }

  /**
   * Internal alias for monic() for backward compatibility.
   */
  _monic(): Polynomial<C> {
    return this.monic();
  }

  /**
   * Return the content of this polynomial (GCD of all coefficients).
   *
   * The content is the GCD of all coefficients. For the zero polynomial,
   * returns zero. Requires that the coefficient ring supports a gcd method.
   *
   * @returns The content (GCD of coefficients)
   *
   * @example
   * ```typescript
   * // If f = 6x^2 + 4x + 2 over ZZ
   * // f.content() = 2
   * ```
   *
   * @see Reference: sage/rings/polynomial/polynomial_integer_dense_flint.pyx:474 (content)
   */
  content(): C {
    if (this.isZero()) {
      return this.parent.base_ring.zero() as C;
    }

    // Start with the first coefficient
    let g = this.coeffs[0]!;

    // Compute GCD with all other coefficients
    for (let i = 1; i < this.coeffs.length; i++) {
      g = gcdCoeffs(g, this.coeffs[i]!);
      // If GCD is 1 (or a unit), we can stop early
      if (g.eq(1)) {
        break;
      }
    }

    // The sign of the content is the sign of the leading coefficient
    // (`polynomial_integer_dense_flint.pyx:477`, issue #13053):
    //     R(-1).content() == -1,  (-2*x^2-4).content() == -2
    if (isNegative(g) !== isNegative(this.leading_coefficient())) {
      g = g.neg() as C;
    }

    return g;
  }

  /**
   * Return the primitive part of this polynomial (this / content).
   *
   * The primitive part is the polynomial divided by its content, so that
   * the GCD of the resulting coefficients is 1.
   *
   * @returns The primitive part
   *
   * @example
   * ```typescript
   * // If f = 6x^2 + 4x + 2 over ZZ
   * // f.primitive_part() = 3x^2 + 2x + 1
   * ```
   *
   * The leading coefficient of the primitive part is always positive, since
   * {@link content} carries the sign of the leading coefficient (this matches
   * FLINT's `fmpz_poly_primitive_part`, see
   * `polynomial_integer_dense_flint.pyx:1535`).
   *
   * @see Reference: sage/libs/flint/fmpz_poly.pxd (fmpz_poly_primitive_part)
   */
  primitive_part(): Polynomial<C> {
    if (this.isZero()) {
      return this;
    }

    const c = this.content();

    // If content is 1, return self
    if (c.eq(1)) {
      return this;
    }

    // Divide each coefficient by the content
    const newCoeffs = this.coeffs.map((coeff) => divideCoeffs(coeff, c));
    return new Polynomial(newCoeffs, this.parent);
  }

  /**
   * Return this polynomial shifted by n (multiplied by x^n).
   *
   * If n is positive, this is equivalent to multiplying by x^n.
   * If n is negative, terms below x^(-n) are discarded (integer division by x^(-n)).
   *
   * @param n - The shift amount (can be negative for division by x^n)
   * @returns x^n * this (or floor division if n < 0)
   *
   * @example
   * ```typescript
   * // If f = x^2 + 2x + 4
   * // f.shift(2) = x^4 + 2x^3 + 4x^2
   * // f.shift(-1) = x + 2
   * ```
   *
   * @see Reference: sage/rings/polynomial/polynomial_element.pyx:shift
   */
  shift(n: number): Polynomial<C> {
    // Zero polynomial or no shift returns self (immutable)
    if (n === 0 || this.degree() < 0) {
      return this;
    }

    if (n > 0) {
      // Multiply by x^n: prepend n zeros
      const output: C[] = [];
      for (let i = 0; i < n; i++) {
        output.push(this.parent.base_ring.zero() as C);
      }
      output.push(...this.coeffs);
      return new Polynomial(output, this.parent);
    }

    // n < 0: divide by x^(-n), i.e., drop lowest (-n) coefficients
    const dropCount = -n;
    if (dropCount > this.coeffs.length - 1) {
      // All coefficients are dropped
      return this.parent.zero();
    }
    return new Polynomial([...this.coeffs.slice(dropCount)] as C[], this.parent);
  }

  /**
   * Return the truncation of this polynomial to degree < n.
   *
   * Returns the polynomial with all terms of degree >= n removed.
   *
   * @param n - The degree bound (must be non-negative)
   * @returns Polynomial with terms of degree >= n removed
   *
   * @example
   * ```typescript
   * // If f = x^3 + 2x^2 + 3x + 4
   * // f.truncate(2) = 3x + 4 (terms of degree < 2)
   * // f.truncate(0) = 0 (no terms of degree < 0)
   * ```
   *
   * @see Reference: sage/rings/polynomial/polynomial_element.pyx:truncate
   */
  truncate(n: number): Polynomial<C> {
    if (n <= 0) {
      return this.parent.zero();
    }

    if (n >= this.coeffs.length) {
      // All coefficients have degree < n, return self
      return this;
    }

    // Keep only coefficients with degree < n
    return new Polynomial([...this.coeffs.slice(0, n)] as C[], this.parent);
  }

  /**
   * Return the reverse of this polynomial.
   *
   * If f(x) = a_0 + a_1*x + ... + a_n*x^n, then
   * reverse(f)(x) = a_n + a_{n-1}*x + ... + a_0*x^n = x^n * f(1/x).
   *
   * If an optional degree argument is given, the coefficient list will be
   * truncated or zero-padded as necessary before reversing.
   *
   * @param degree - Optional degree to use (pads with zeros or truncates)
   * @returns Polynomial with coefficients in reverse order
   *
   * @example
   * ```typescript
   * // If f = x^3 + 2x + 3 (coeffs: [3, 2, 0, 1])
   * // f.reverse() = 1 + 0*x + 2*x^2 + 3*x^3 = 1 + 2x^2 + 3x^3
   * ```
   *
   * @see Reference: sage/rings/polynomial/polynomial_element.pyx:reverse
   */
  reverse(degree?: number): Polynomial<C> {
    if (this.isZero()) {
      return this;
    }

    let v = [...this.coeffs] as C[];

    if (degree !== undefined) {
      if (degree < 0) {
        throw new ValueError(`degree argument must be a nonnegative integer, got ${degree}`);
      }

      const targetLen = degree + 1;
      if (v.length < targetLen) {
        // Reverse first, then prepend zeros
        v.reverse();
        const padding: C[] = [];
        for (let i = 0; i < targetLen - v.length; i++) {
          padding.push(this.parent.base_ring.zero() as C);
        }
        v = [...padding, ...v];
      } else if (v.length > targetLen) {
        // Truncate to first (degree+1) coefficients, then reverse
        v = v.slice(0, targetLen);
        v.reverse();
      } else {
        // v.length === targetLen
        v.reverse();
      }
    } else {
      v.reverse();
    }

    return new Polynomial(v, this.parent);
  }

  /**
   * Return the resultant of this polynomial and other.
   *
   * The resultant of two polynomials f and g is the determinant of their
   * Sylvester matrix. It is zero if and only if f and g have a common root.
   *
   * @param other - Another polynomial in the same ring
   * @returns The resultant (an element of the base ring)
   *
   * @example
   * ```typescript
   * // If f = x^3 + x + 1 and g = x^3 - x - 1
   * // f.resultant(g) = -8
   * ```
   *
   * @see Reference: sage/rings/polynomial/polynomial_element.pyx:resultant
   */
  resultant(other: Polynomial<C>): C {
    // Handle zero polynomials
    if (this.isZero() || other.isZero()) {
      return this.parent.base_ring.zero() as C;
    }

    const m = this.degree();
    const n = other.degree();

    // If both are constants, resultant is 1 (empty matrix determinant)
    if (m === 0 && n === 0) {
      return this.parent.base_ring.one() as C;
    }

    // Handle constant polynomial cases
    if (m === 0) {
      // Res(c, g) = c^deg(g)
      let result = this.coeffs[0]!;
      for (let i = 1; i < n; i++) {
        result = result.mul(this.coeffs[0]!) as C;
      }
      return result;
    }

    if (n === 0) {
      // Res(f, c) = c^deg(f)
      let result = other.coeffs[0]!;
      for (let i = 1; i < m; i++) {
        result = result.mul(other.coeffs[0]!) as C;
      }
      return result;
    }

    // Build and compute Sylvester matrix determinant
    return matrixDeterminant(this.sylvester_matrix(other), this.parent.base_ring);
  }

  /**
   * Return the Sylvester matrix of this polynomial and `other`.
   *
   * For `deg(self) = m` and `deg(other) = n` this is the `(m+n) x (m+n)`
   * matrix whose first `n` rows hold the coefficients of `x^i * self` and
   * whose last `m` rows hold the coefficients of `x^i * other`.
   *
   * @see Reference: sage/rings/polynomial/polynomial_element.pyx:sylvester_matrix
   */
  sylvester_matrix(other: Polynomial<C>): C[][] {
    const m = this.degree();
    const n = other.degree();
    const size = m + n;
    const matrix: C[][] = [];

    // Initialize matrix with zeros
    for (let i = 0; i < size; i++) {
      const row: C[] = [];
      for (let j = 0; j < size; j++) {
        row.push(this.parent.base_ring.zero() as C);
      }
      matrix.push(row);
    }

    // Fill in rows for f (n rows)
    // Row i contains coefficients of x^i * f
    for (let i = 0; i < n; i++) {
      for (let j = 0; j <= m; j++) {
        matrix[i]![i + (m - j)] = this.getCoeff(j);
      }
    }

    // Fill in rows for g (m rows)
    // Row n+i contains coefficients of x^i * g
    for (let i = 0; i < m; i++) {
      for (let j = 0; j <= n; j++) {
        matrix[n + i]![i + (n - j)] = other.getCoeff(j);
      }
    }

    return matrix;
  }

  /**
   * Return the discriminant of this polynomial.
   *
   * The discriminant is defined as:
   *   disc(f) = (-1)^(n(n-1)/2) * Res(f, f') / a_n
   *
   * where n is the degree, a_n is the leading coefficient, and f' is the derivative.
   *
   * The discriminant is zero if and only if the polynomial has a repeated root.
   *
   * @returns The discriminant (an element of the base ring)
   *
   * @example
   * ```typescript
   * // If f = x^3 + x + 1
   * // f.discriminant() = -31
   * ```
   *
   * @see Reference: sage/rings/polynomial/polynomial_element.pyx:discriminant
   */
  discriminant(): C {
    if (this.isZero()) {
      return this.parent.base_ring.zero() as C;
    }

    const n = this.degree();
    if (n <= 0) {
      // Constant polynomial has discriminant 1 (or the leading coefficient for degree 0)
      return this.parent.base_ring.one() as C;
    }

    const d = this.derivative();
    const k = d.degree();
    const an = this.leading_coefficient();

    // Compute sign: (-1)^(n*(n-1)/2)
    // n*(n-1)/2 mod 2:
    //   n=0: 0 -> +1
    //   n=1: 0 -> +1
    //   n=2: 1 -> -1
    //   n=3: 3 -> -1
    //   n=4: 6 -> +1
    //   n=5: 10 -> +1
    // Pattern: sign is -1 when n mod 4 is 2 or 3
    const r = n % 4;
    const signIsNegative = r === 2 || r === 3;

    // Compute a_n^(n - k - 2) where k = deg(f')
    // Normally k = n - 1, so n - k - 2 = n - (n-1) - 2 = -1
    // This means we need to divide by a_n
    const exponent = n - k - 2;

    // Compute resultant
    const res = this.resultant(d);

    let result: C;
    if (exponent >= 0) {
      // Multiply by a_n^exponent
      let anPower = this.parent.base_ring.one() as C;
      for (let i = 0; i < exponent; i++) {
        anPower = anPower.mul(an) as C;
      }
      result = res.mul(anPower) as C;
    } else {
      // Divide by a_n^(-exponent)
      let anPower = an;
      for (let i = 1; i < -exponent; i++) {
        anPower = anPower.mul(an) as C;
      }
      const quotient = divideCoeffs(res, anPower);
      if (quotient.mul(anPower).eq(res)) {
        result = quotient;
      } else {
        // Division by the leading coefficient is not exact in the base ring.
        // Rather than dividing the resultant, alter the Sylvester matrix
        // (Sage issue #11782, `polynomial_element.pyx:8094-8099`).
        if (exponent !== -1) {
          throw new ArithmeticError('discriminant: division by the leading coefficient failed');
        }
        const mat = this.sylvester_matrix(d);
        mat[0]![0] = this.parent.base_ring.one() as C;
        mat[n - 1]![0] = mulByInteger(this.parent.base_ring.one() as C, n, this.parent.base_ring);
        result = matrixDeterminant(mat, this.parent.base_ring);
      }
    }

    // Apply sign
    if (signIsNegative) {
      result = result.neg() as C;
    }

    return result;
  }

  /**
   * Return the roots of this polynomial in the base ring.
   *
   * For finite fields, this finds all roots by either:
   * - Trying all elements (for small fields)
   * - Using factorization for larger fields
   *
   * @returns List of roots with multiplicities as Array<[root, multiplicity]>
   *
   * @example
   * ```typescript
   * // Over GF(7): x^2 - 1 = (x-1)(x+1) has roots 1 and 6
   * const p = x.pow(2).sub(R.one());
   * const roots = p.roots(); // [[1, 1], [6, 1]]
   * ```
   *
   * @see Reference: sage/rings/polynomial/polynomial_element.pyx:roots
   * @see Deviation: Polynomial Roots and Factorization
   */
  roots(): Array<[C, number]> {
    if (this.isZero()) {
      throw new ValueError('roots of zero polynomial are not defined');
    }

    if (this.degree() === 0) {
      return []; // Constant non-zero polynomial has no roots
    }

    // Check if base ring is a finite field (has cardinality method and is iterable)
    const baseRing = this.parent.base_ring;

    // Handle integer polynomials (ZZ[x])
    if (isIntegerRing(baseRing)) {
      const coeffs = extractIntegerCoeffs(this);
      const intRoots = findIntegerRoots(coeffs);
      // Convert back to ring elements
      return sortRootsSageOrder(
        intRoots.map(([root, mult]) => [baseRing.__call__(root) as C, mult])
      );
    }

    // Handle rational polynomials (QQ[x]) - find rational roots
    if (isRationalField(baseRing)) {
      // For QQ[x], we can find rational roots using the rational root theorem
      // First, clear denominators to get an integer polynomial
      const rationalRoots = findRationalRoots(this);
      return sortRootsSageOrder(rationalRoots);
    }

    // Check if it's a finite field we can iterate over
    if (!isFiniteField(baseRing)) {
      throw new NotImplementedError('roots only implemented for finite fields, ZZ, and QQ');
    }

    const roots: Array<[C, number]> = [];

    // For small fields, try all elements
    const order = getFieldOrder(baseRing);

    if (order <= 10000n) {
      // Small field: try all elements directly
      let f = this as Polynomial<C>;

      for (const elem of iterateField(baseRing)) {
        if (f.evaluate(elem as C).isZero()) {
          // Found a root, count its multiplicity
          let mult = 0;
          const linearFactor = this._linearFactor(elem as C);

          while (f.degree() >= 1) {
            const [q, r] = f.quo_rem(linearFactor);
            if (!r.isZero()) {
              break;
            }
            f = q;
            mult++;
          }

          if (mult > 0) {
            roots.push([elem as C, mult]);
          }
        }
      }
    } else {
      // Large field: use factorization
      const factors = this.factor();

      for (const [fac, mult] of factors) {
        if (fac.degree() === 1) {
          // Linear factor (x - a), extract root a = -c_0 / c_1
          const c0 = fac.getCoeff(0);
          const c1 = fac.getCoeff(1);
          const root = divideCoeffs(c0.neg() as C, c1);
          roots.push([root, mult]);
        }
      }
    }

    return sortRootsSageOrder(roots);
  }

  /**
   * Return (x - a) where a is the given root.
   * Helper for roots computation.
   */
  private _linearFactor(root: C): Polynomial<C> {
    return new Polynomial([root.neg() as C, this.parent.base_ring.one() as C], this.parent);
  }

  /**
   * Factor this polynomial over the integers ZZ.
   * Returns pairs [irreducible_factor, multiplicity].
   *
   * The content is returned as its own prime factors, as Sage does
   * (`(12*(x^2+1)^3*(x+2)).factor()` is `2^2 * 3 * (x + 2) * (x^2 + 1)^3`),
   * and a negative content contributes the unit `-1`.
   *
   * @see Deviation: Polynomial Roots and Factorization
   */
  private _factorOverIntegers(): Array<[Polynomial<C>, number]> {
    const coeffs = extractIntegerCoeffs(this);
    const [content, factors] = factorIntegerPolynomial(coeffs);

    const result: Array<[Polynomial<C>, number]> = [];

    // The sign of the content is the unit of the factorization; Sage keeps it
    // in ``Factorization.unit()`` (e.g. ``(-x^2+4).factor() == (-1)*(x-2)*(x+2)``).
    if (content < 0n) {
      result.push([new Polynomial([this.parent.base_ring.__call__(-1n) as C], this.parent), 1]);
    }

    // Add content as a factor if it's not 1 or -1
    if (content !== 1n && content !== -1n) {
      // Factor the integer content
      const intFactors = factorInteger(content < 0n ? -content : content);
      for (const [p, e] of intFactors) {
        if (p > 1n) {
          const constPoly = new Polynomial([this.parent.base_ring.__call__(p) as C], this.parent);
          result.push([constPoly, Number(e)]);
        }
      }
    }

    // Convert polynomial factors back to Polynomial<C>
    for (const [facCoeffs, mult] of factors) {
      const polyCoeffs = facCoeffs.map((c) => this.parent.base_ring.__call__(c) as C);
      const poly = new Polynomial(polyCoeffs, this.parent);
      result.push([poly, mult]);
    }

    // Sort factors by degree, then lexicographically
    result.sort((a, b) => {
      if (a[0].degree() !== b[0].degree()) {
        return a[0].degree() - b[0].degree();
      }
      return a[0].toString().localeCompare(b[0].toString());
    });

    return result;
  }

  /**
   * Factor this polynomial over the rationals QQ.
   * Returns pairs [monic_irreducible_factor, multiplicity].
   *
   * Sage keeps the leading coefficient in `Factorization.unit()`; we return it
   * as a degree-0 factor so that the product of the returned factors is again
   * `self`.
   *
   * @see Deviation: Polynomial Roots and Factorization
   */
  private _factorOverRationals(): Array<[Polynomial<C>, number]> {
    // Clear denominators to get an integer polynomial, factor that over ZZ and
    // convert back to monic factors over QQ.  Neither the common denominator
    // nor the integer content matters here: both are units of QQ[x] and are
    // accounted for by the leading coefficient added below.
    const [intCoeffs] = clearDenominators(this);

    const [, factors] = factorIntegerPolynomial(intCoeffs);

    const result: Array<[Polynomial<C>, number]> = [];

    // Convert polynomial factors back to monic Polynomial<C> over QQ.
    // The quotient is formed with the base ring's own division so that any
    // rational-like coefficient ring works, not only ones whose `__call__`
    // understands a numerator/denominator pair.
    for (const [facCoeffs, mult] of factors) {
      // Make monic by dividing by leading coefficient
      const lc = this.parent.base_ring.__call__(facCoeffs[facCoeffs.length - 1]!) as C;
      const monicCoeffs = facCoeffs.map((c) =>
        divideCoeffs(this.parent.base_ring.__call__(c) as C, lc)
      );
      const poly = new Polynomial(monicCoeffs, this.parent);
      result.push([poly, mult]);
    }

    // All factors are monic, so the leading coefficient of ``self`` is the
    // unit of the factorization (Sage keeps it in ``Factorization.unit()``).
    const unit = this.leading_coefficient();
    if (!unit.eq(1)) {
      result.push([new Polynomial([unit], this.parent), 1]);
    }

    // Sort factors by degree, then lexicographically
    result.sort((a, b) => {
      if (a[0].degree() !== b[0].degree()) {
        return a[0].degree() - b[0].degree();
      }
      return a[0].toString().localeCompare(b[0].toString());
    });

    return result;
  }

  /**
   * Return the squarefree decomposition of this polynomial.
   *
   * Returns a list of pairs (f_i, i) where this polynomial equals
   * prod(f_i^i) and each f_i is squarefree and coprime to the others.
   *
   * Uses the standard algorithm based on gcd with the derivative.
   *
   * @returns Array of [squarefree_factor, multiplicity] pairs
   *
   * @example
   * ```typescript
   * // (x-1)^2 * (x-2) has squarefree decomposition [(x-2, 1), ((x-1), 2)]
   * ```
   *
   * @see Reference: sage/rings/polynomial/polynomial_element.pyx:squarefree_decomposition
   */
  squarefree_decomposition(): Array<[Polynomial<C>, number]> {
    if (this.isZero()) {
      throw new ValueError('squarefree decomposition of zero polynomial is not defined');
    }

    if (this.degree() === 0) {
      // Constant polynomial is squarefree
      return [[this, 1]];
    }

    const result: Array<[Polynomial<C>, number]> = [];

    // Get the characteristic of the base field (if finite)
    const baseRing = this.parent.base_ring;
    const p = getCharacteristic(baseRing);

    // Make polynomial monic for easier computation
    const f = this._monic();

    // Standard squarefree decomposition algorithm
    const d = f.derivative();

    if (d.isZero()) {
      // Derivative is zero, meaning all exponents are divisible by p
      // This can only happen in characteristic p > 0
      if (p === 0n) {
        // Should not happen for non-zero derivative
        return [[f, 1]];
      }

      // f = g^p for some g, find g by taking p-th roots of coefficients
      const g = this._pthRoot(p);
      const subDecomp = g.squarefree_decomposition();

      // Multiply all multiplicities by p
      for (const [factor, mult] of subDecomp) {
        result.push([factor, mult * Number(p)]);
      }

      return result;
    }

    // gcd(f, f') gives us the product of repeated factors
    let g = f.gcd(d);
    let h = f.quo_rem(g)[0]; // f / gcd(f, f')

    let i = 1;

    while (!h.eq(this.parent.one())) {
      // g_i = gcd(g, h)
      const gi = g.gcd(h);
      // h_i = h / g_i (the squarefree part with multiplicity i)
      const hi = h.quo_rem(gi)[0];

      if (!hi.eq(this.parent.one())) {
        result.push([hi, i]);
      }

      // Update for next iteration
      g = g.quo_rem(gi)[0];
      h = gi;
      i++;
    }

    // If g is not 1, then g = (product of factors)^p in characteristic p
    if (!g.eq(this.parent.one()) && p > 0n) {
      const gRoot = g._pthRoot(p);
      const subDecomp = gRoot.squarefree_decomposition();

      for (const [factor, mult] of subDecomp) {
        result.push([factor, mult * Number(p)]);
      }
    }

    return result;
  }

  /**
   * Compute the p-th root of a polynomial (in characteristic p).
   * Used for squarefree decomposition when derivative is zero.
   */
  private _pthRoot(p: bigint): Polynomial<C> {
    // In characteristic p, if f = sum(a_i x^{ip}), then f^{1/p} = sum(a_i^{1/p} x^i)
    const newCoeffs: C[] = [];

    for (let i = 0; i < this.coeffs.length; i++) {
      if (BigInt(i) % p === 0n) {
        // Take p-th root of coefficient
        // For finite fields, this is a^{q/p} where q = |F|
        const coeff = this.coeffs[i]!;
        const rootCoeff = pthRootCoeff(coeff, p);
        newCoeffs.push(rootCoeff);
      }
    }

    return new Polynomial(newCoeffs, this.parent);
  }

  /**
   * Return the distinct-degree factorization of this polynomial.
   *
   * Assumes the polynomial is squarefree. Returns a list of pairs (g_d, d)
   * where g_d is the product of all irreducible factors of degree d.
   *
   * This is a key step in the Berlekamp and Cantor-Zassenhaus algorithms.
   *
   * @returns Array of [product_of_degree_d_factors, d] pairs
   *
   * @example
   * ```typescript
   * // Over GF(163), factoring (x+162)(x^3+7x+161)(x^7+9x+161)
   * // gives [(x+162, 1), (x^3+7x+161, 3), (x^7+9x+161, 7)]
   * ```
   *
   * @see Reference: sage/rings/polynomial/polynomial_element.pyx:_distinct_degree_factorisation_squarefree
   */
  distinct_degree_factorization(): Array<[Polynomial<C>, number]> {
    if (this.isZero()) {
      throw new ValueError('distinct-degree factorization of zero polynomial is not defined');
    }

    const baseRing = this.parent.base_ring;

    if (!isFiniteField(baseRing)) {
      throw new NotImplementedError(
        'distinct-degree factorization only implemented for finite fields'
      );
    }

    const q = getFieldOrder(baseRing);
    const result: Array<[Polynomial<C>, number]> = [];

    // x
    const x = this.parent.gen();
    // Work with monic polynomial
    let v = this._monic();
    // w = x mod v
    let w = x.mod(v);
    let d = 0;
    let e = v.degree();

    // Iterate over all possible degrees
    while (2 * (d + 1) <= e) {
      d = d + 1;

      // w = w^q mod v (computing x^{q^d} mod v)
      w = powerMod(w, q, v);

      // a_d = gcd(v, w - x)
      const wMinusX = w.sub(x);
      const ad = v.gcd(wMinusX);

      if (!ad.eq(this.parent.one())) {
        result.push([ad._monic(), d]);
        // v = v / a_d
        v = v.quo_rem(ad)[0];
      }

      e = v.degree();
    }

    // If v still has positive degree, it's irreducible of degree e
    if (e > 0) {
      result.push([v._monic(), e]);
    }

    return result;
  }

  /**
   * Return the factorization of this polynomial.
   *
   * For finite fields, uses:
   * 1. Squarefree decomposition
   * 2. Distinct-degree factorization
   * 3. Cantor-Zassenhaus algorithm for equal-degree factorization
   *
   * @returns Factorization as list of [factor, multiplicity] pairs
   *
   * @example
   * ```typescript
   * // x^4 - 1 over GF(5) factors as (x-1)(x+1)(x-2)(x+2)
   * const factors = (x.pow(4).sub(R.one())).factor();
   * ```
   *
   * @see Reference: sage/rings/polynomial/polynomial_element.pyx:factor
   * @see Deviation: Polynomial Roots and Factorization
   */
  factor(): Array<[Polynomial<C>, number]> {
    if (this.isZero()) {
      throw new ValueError('factorization of zero polynomial is not defined');
    }

    if (this.degree() === 0) {
      // Constant polynomial
      return [[this, 1]];
    }

    const baseRing = this.parent.base_ring;

    // Handle integer polynomials (ZZ[x])
    if (isIntegerRing(baseRing)) {
      return this._factorOverIntegers();
    }

    // Handle rational polynomials (QQ[x])
    if (isRationalField(baseRing)) {
      return this._factorOverRationals();
    }

    if (!isFiniteField(baseRing)) {
      throw new NotImplementedError('factorization only implemented for finite fields, ZZ, and QQ');
    }

    const result: Array<[Polynomial<C>, number]> = [];

    // Step 1: Squarefree decomposition
    const sqfree = this.squarefree_decomposition();

    // Step 2: For each squarefree factor, do distinct-degree and equal-degree factorization
    for (const [sqfFactor, mult] of sqfree) {
      if (sqfFactor.degree() === 0) {
        // Constant factor
        if (!sqfFactor.leading_coefficient().eq(1)) {
          result.push([sqfFactor, mult]);
        }
        continue;
      }

      // Distinct-degree factorization
      const ddf = sqfFactor.distinct_degree_factorization();

      // Equal-degree factorization for each degree
      for (const [ddFactor, degree] of ddf) {
        if (ddFactor.degree() === degree) {
          // ddFactor is already irreducible
          result.push([ddFactor, mult]);
        } else {
          // Need to split ddFactor into irreducible factors of the given degree
          const irreducibles = cantorZassenhausFactorization(ddFactor, degree, this.parent);

          for (const irr of irreducibles) {
            result.push([irr, mult]);
          }
        }
      }
    }

    // The factors above are all monic, so the leading coefficient of ``self``
    // is the unit of the factorization.  Sage keeps it in
    // ``Factorization.unit()``; we return it as a degree-0 factor so that the
    // product of the returned factors is again ``self``.
    const unit = this.leading_coefficient();
    if (!unit.eq(1)) {
      result.push([new Polynomial([unit], this.parent), 1]);
    }

    // Sort factors by degree, then lexicographically
    result.sort((a, b) => {
      if (a[0].degree() !== b[0].degree()) {
        return a[0].degree() - b[0].degree();
      }
      return a[0].toString().localeCompare(b[0].toString());
    });

    return result;
  }

  /**
   * Test if this polynomial is irreducible.
   *
   * Follows Sage: the zero polynomial and units are reducible, a constant is
   * irreducible iff it is irreducible in the base ring, and otherwise the
   * polynomial is factored (over finite fields we use Rabin's test, which is
   * what FLINT's `nmod_poly_is_irreducible_rabin` does).
   *
   * @returns true if irreducible
   *
   * @example
   * ```typescript
   * // x^2 + x + 1 is irreducible over GF(2)
   * const p = x.pow(2).add(x).add(R.one());
   * p.is_irreducible(); // true
   * ```
   *
   * @see Reference: sage/rings/polynomial/polynomial_element.pyx:10182 (is_irreducible)
   */
  is_irreducible(): boolean {
    if (this.isZero()) {
      return false;
    }

    const baseRing = this.parent.base_ring;
    const n = this.degree();

    if (n === 0) {
      // Sage: ``if self.is_unit(): return False`` then defers to the base
      // ring, so ZZ(5) is irreducible while ZZ(4), ZZ(1) and any nonzero
      // element of a field are not.
      const c = this.coeffs[0]!;
      if (isIntegerRing(baseRing)) {
        const v = extractIntegerCoeffs(this)[0]!;
        const a = v < 0n ? -v : v;
        return a > 1n && is_prime(a);
      }
      if (
        'is_irreducible' in c &&
        typeof (c as unknown as { is_irreducible: () => boolean }).is_irreducible === 'function'
      ) {
        return (c as unknown as { is_irreducible: () => boolean }).is_irreducible();
      }
      // Every nonzero constant is a unit over a field.
      return false;
    }

    // Handle integer polynomials (ZZ[x])
    if (isIntegerRing(baseRing)) {
      // A polynomial is irreducible over ZZ iff it is primitive and
      // irreducible over QQ (Gauss's lemma)
      const coeffs = extractIntegerCoeffs(this);
      const content = intPolyContent(coeffs);
      if (content !== 1n && content !== -1n) {
        return false; // Not primitive
      }
      if (n === 1) {
        return true;
      }
      const [_, factors] = factorIntegerPolynomial(coeffs);
      return factors.length === 1 && factors[0]![1] === 1;
    }

    // Handle rational polynomials (QQ[x])
    if (isRationalField(baseRing)) {
      if (n === 1) {
        return true;
      }
      // Sage tests ``len(F) > 1 or F[0][1] > 1`` on a ``Factorization`` whose
      // unit is kept apart (polynomial_element.pyx:is_irreducible); we return
      // that unit as a degree-0 factor instead, so it must not be counted --
      // over a field it is a unit.  This is why ``2*x^2 + 2`` is irreducible
      // in ``QQ[x]`` (and not in ``ZZ[x]``).
      const factors = this.factor().filter(([g]) => g.degree() > 0);
      return factors.length === 1 && factors[0]![1] === 1;
    }

    if (!isFiniteField(baseRing)) {
      throw new NotImplementedError(
        'is_irreducible only implemented for finite fields, ZZ, and QQ'
      );
    }

    if (n === 1) {
      return true; // Linear polynomials are irreducible over a field
    }

    // Rabin's irreducibility test over GF(q) (FLINT
    // `nmod_poly_factor/is_irreducible.c:nmod_poly_is_irreducible_rabin`):
    // f of degree n is irreducible iff x^(q^n) = x mod f and
    // gcd(x^(q^(n/l)) - x, f) = 1 for every prime l | n.
    const q = getFieldOrder(baseRing);
    const x = this.parent.gen();
    const monic = this._monic();

    // x^(q^n) mod f
    const xqn = powerModIterated(x, q, n, monic);
    if (!xqn.eq(x)) {
      return false;
    }

    for (const [l] of factorInteger(BigInt(n))) {
      if (l <= 1n) continue;
      const a = powerModIterated(x, q, n / Number(l), monic).sub(x);
      if (a.isZero()) {
        return false;
      }
      const g = monic.gcd(a);
      if (g.degree() > 0) {
        return false;
      }
    }

    return true;
  }

  /**
   * String representation.
   */
  toString(): string {
    if (this.coeffs.length === 0) {
      return '0';
    }

    const varName = this.parent.variable_name;
    const terms: string[] = [];

    for (let i = this.coeffs.length - 1; i >= 0; i--) {
      const c = this.coeffs[i]!;
      if (c.isZero()) {
        continue;
      }

      let term: string;
      const cStr = c.toString();

      if (i === 0) {
        term = cStr;
      } else if (i === 1) {
        if (c.eq(1)) {
          term = varName;
        } else {
          term = needsParens(cStr) ? `(${cStr})*${varName}` : `${cStr}*${varName}`;
        }
      } else {
        if (c.eq(1)) {
          term = `${varName}^${i}`;
        } else {
          term = needsParens(cStr) ? `(${cStr})*${varName}^${i}` : `${cStr}*${varName}^${i}`;
        }
      }

      terms.push(term);
    }

    if (terms.length === 0) {
      return '0';
    }

    return terms.join(' + ');
  }
}

/**
 * Check if a coefficient string needs parentheses when multiplied.
 */
function needsParens(s: string): boolean {
  return s.includes('+') || s.includes('-') || s.includes('*');
}

/**
 * Sort a root list into SageMath's order.
 *
 * `Polynomial.roots()` reads the linear factors off `factor()`, and
 * `Factorization.sort` (`sage/structure/factorization.py:671-741`) orders the
 * factors by `(prime.degree(), exponent, prime)`.  Every linear factor has
 * degree 1, so what remains is `(multiplicity, x - r)`, and Sage compares
 * polynomials by their coefficient list in ASCENDING order
 * (`polynomial_element.pyx` `_richcmp_` -> `richcmp(self.list(), ...)`), whose
 * first entry is `-r`.
 *
 * The net effect is: multiplicity ascending, then the root DESCENDING in the
 * ring's own representation.  Concretely `(y^2-1).roots()` over `GF(11)` is
 * `[10, 1]`, `(y^2+7y+3).roots()` is `[3, 1]`, and over `QQ`
 * `(x^5+x^3-x^2-x).roots()` is `[1, 0]`.
 *
 * This is load-bearing wherever a caller takes `roots()[0]`: Jordan block
 * ordering (`matrix2.pyx:12251-12254`), `cantor_reduction`'s choice of point at
 * infinity, `odd_degree_model`'s choice of root.
 *
 * Elements whose value cannot be read as an exact rational keep their
 * insertion order (the sort is stable).
 */
function sortRootsSageOrder<C extends RingElement>(roots: Array<[C, number]>): Array<[C, number]> {
  const keyed = roots.map((entry, idx) => ({ entry, idx, key: exactRootKey(entry[0]) }));
  keyed.sort((a, b) => {
    if (a.entry[1] !== b.entry[1]) {
      return a.entry[1] - b.entry[1];
    }
    if (a.key === null || b.key === null) {
      return a.idx - b.idx;
    }
    // Compare `-r` ascending, i.e. `r` descending.
    const lhs = -a.key.n * b.key.d;
    const rhs = -b.key.n * a.key.d;
    if (lhs < rhs) {
      return -1;
    }
    if (lhs > rhs) {
      return 1;
    }
    return a.idx - b.idx;
  });
  return keyed.map((k) => k.entry);
}

/** The exact rational value of a coefficient, or `null` if it has none. */
function exactRootKey(x: unknown): { n: bigint; d: bigint } | null {
  const v = (x as { value?: unknown }).value;
  if (typeof v === 'bigint') {
    return { n: v, d: 1n };
  }
  const num = (x as { numerator?: unknown }).numerator;
  const den = (x as { denominator?: unknown }).denominator;
  const numV = typeof num === 'function' ? (num as () => unknown).call(x) : num;
  const denV = typeof den === 'function' ? (den as () => unknown).call(x) : den;
  if (typeof numV === 'bigint' && typeof denV === 'bigint' && denV > 0n) {
    return { n: numV, d: denV };
  }
  return null;
}

/**
 * Divide two coefficients. Assumes the ring supports division.
 */
function divideCoeffs<C extends RingElement>(a: C, b: C): C {
  // Try to call div method if it exists
  if ('div' in a && typeof (a as unknown as { div: (b: C) => C }).div === 'function') {
    return (a as unknown as { div: (b: C) => C }).div(b);
  }

  // Try inv method for field elements
  if ('inv' in b && typeof (b as unknown as { inv: () => C }).inv === 'function') {
    const bInv = (b as unknown as { inv: () => C }).inv();
    return a.mul(bInv) as C;
  }

  throw new ValueError('coefficient ring does not support division');
}

/**
 * Return the inverse of `c` when `c` is a unit of the coefficient ring, and
 * `null` otherwise.
 *
 * This is Sage's `inverse_of_unit()`: it must not succeed for a non-unit (over
 * ZZ, `2` has no inverse), so the candidate inverse is verified.  Inexact
 * rings (where `c * c^-1` is only approximately one) keep working because the
 * verification is done with the ring's own equality.
 */
function inverseOfUnit<C extends RingElement>(c: C, ring: CoefficientRing<C>): C | null {
  const withInv = c as unknown as { inv?: () => C };
  if (typeof withInv.inv !== 'function') {
    return null;
  }
  let inv: C;
  try {
    inv = withInv.inv();
  } catch {
    // Sage catches ArithmeticError/ValueError from inverse_of_unit here.
    return null;
  }
  if (!c.mul(inv).eq(ring.one())) {
    return null;
  }
  return inv;
}

/**
 * Multiply a ring element by a non-negative integer using double-and-add.
 *
 * This is `n * c` in the base ring, computed with O(log n) additions instead
 * of n-1 of them.
 */
function mulByInteger<C extends RingElement>(coeff: C, n: number, ring: CoefficientRing<C>): C {
  if (n === 0) {
    return ring.zero() as C;
  }
  let k = n;
  let acc: C | null = null;
  let addend = coeff;
  while (k > 0) {
    if (k & 1) {
      acc = acc === null ? addend : (acc.add(addend) as C);
    }
    k >>= 1;
    if (k > 0) {
      addend = addend.add(addend) as C;
    }
  }
  return acc ?? (ring.zero() as C);
}

/**
 * Return whether a coefficient is negative (meaningful only in ordered rings
 * such as ZZ and QQ; always false elsewhere).
 */
function isNegative<C extends RingElement>(c: C): boolean {
  if ('value' in c) {
    const v = (c as unknown as { value: unknown }).value;
    if (typeof v === 'bigint') return v < 0n;
    if (typeof v === 'number') return v < 0;
  }
  if ('numerator' in c) {
    const num = (c as unknown as { numerator: unknown }).numerator;
    if (typeof num === 'bigint') return num < 0n;
  }
  return c.toString().startsWith('-');
}

/**
 * Return whether a coefficient ring is a field.
 */
function ringIsField<C extends RingElement>(ring: CoefficientRing<C>): boolean {
  if (typeof ring.is_field === 'function') {
    return ring.is_field();
  }
  if (isRationalField(ring)) {
    return true;
  }
  if (isIntegerRing(ring)) {
    return false;
  }
  // Fall back on whether elements can be inverted.
  const one = ring.one();
  return 'inv' in one && typeof (one as unknown as { inv: unknown }).inv === 'function';
}

/**
 * Compute GCD of two coefficients. Assumes the ring supports a gcd method.
 */
function gcdCoeffs<C extends RingElement>(a: C, b: C): C {
  // Try to call gcd method if it exists
  if ('gcd' in a && typeof (a as unknown as { gcd: (b: C) => C }).gcd === 'function') {
    return (a as unknown as { gcd: (b: C) => C }).gcd(b);
  }

  // For fields, GCD is always 1 (or a unit)
  // Check if the ring is a field by looking for inv method
  if ('inv' in a && typeof (a as unknown as { inv: () => C }).inv === 'function') {
    // In a field, gcd(a, b) = 1 for any nonzero a, b
    if (!a.isZero()) {
      // Return the multiplicative identity
      return a.mul((a as unknown as { inv: () => C }).inv()) as C; // This gives 1
    }
    if (!b.isZero()) {
      return b.mul((b as unknown as { inv: () => C }).inv()) as C;
    }
    return a; // Both zero, return zero
  }

  throw new ValueError('coefficient ring does not support gcd');
}

/**
 * Compute the determinant of a square matrix over an integral domain using
 * fraction-free (Bareiss) Gaussian elimination.
 *
 * Every division performed here is exact, so this is valid over any integral
 * domain -- in particular over ZZ, where the previous division-based
 * elimination silently truncated and returned wrong resultants/discriminants.
 *
 * @see Reference: sage/matrix/matrix2.pyx (determinant, "df" / Bareiss)
 */
function matrixDeterminant<C extends RingElement>(matrix: C[][], ring: CoefficientRing<C>): C {
  const n = matrix.length;
  if (n === 0) {
    return ring.one() as C;
  }

  // Make a copy of the matrix
  const M: C[][] = matrix.map((row) => [...row]);

  let sign = 1;
  let prevPivot = ring.one() as C;

  for (let col = 0; col < n - 1; col++) {
    // Find pivot
    let pivotRow = -1;
    for (let row = col; row < n; row++) {
      if (!M[row]![col]!.isZero()) {
        pivotRow = row;
        break;
      }
    }

    if (pivotRow === -1) {
      // Column is all zeros on and below the diagonal: determinant is 0
      return ring.zero() as C;
    }

    // Swap rows if needed
    if (pivotRow !== col) {
      [M[col], M[pivotRow]] = [M[pivotRow]!, M[col]!];
      sign = -sign;
    }

    const pivot = M[col]![col]!;

    for (let row = col + 1; row < n; row++) {
      for (let j = col + 1; j < n; j++) {
        // M[row][j] = (M[row][j]*pivot - M[row][col]*M[col][j]) / prevPivot
        const numer = M[row]![j]!.mul(pivot).sub(M[row]![col]!.mul(M[col]![j]!) as C) as C;
        M[row]![j] = prevPivot.eq(1) ? numer : divideCoeffs(numer, prevPivot);
      }
      M[row]![col] = ring.zero() as C;
    }

    prevPivot = pivot;
  }

  const det = M[n - 1]![n - 1]!;
  return sign === -1 ? (det.neg() as C) : det;
}

/**
 * Get characteristic from a ring (handles both property and method).
 */
function getRingCharacteristic<C extends RingElement>(ring: CoefficientRing<C>): bigint | null {
  if (!('characteristic' in ring)) return null;

  const char = (ring as { characteristic: unknown }).characteristic;

  // If it's a function, call it
  if (typeof char === 'function') {
    const result = (char as () => bigint | number)();
    return typeof result === 'number' ? BigInt(result) : result;
  }

  // If it's already a value
  if (typeof char === 'bigint') return char;
  if (typeof char === 'number') return BigInt(char);

  return null;
}

/**
 * Check if a ring is the integer ring ZZ.
 */
function isIntegerRing<C extends RingElement>(ring: CoefficientRing<C>): boolean {
  // Check for IntegerRing signature: is_field() returns false, characteristic() returns 0
  // and toString() returns 'Integer Ring'
  if (ring.toString && ring.toString() === 'Integer Ring') {
    return true;
  }
  // Check for duck typing: has is_field, is_integral_domain, characteristic
  if ('is_field' in ring && 'is_integral_domain' in ring && 'characteristic' in ring) {
    const r = ring as { is_field: () => boolean; is_integral_domain: () => boolean };
    const char = getRingCharacteristic(ring);
    if (!r.is_field() && r.is_integral_domain() && char === 0n) {
      // Check it's not QQ by verifying there's no fraction field marker
      if (!('is_absolute' in ring)) {
        return true;
      }
    }
  }
  return false;
}

/**
 * Check if a ring is the rational field QQ.
 */
function isRationalField<C extends RingElement>(ring: CoefficientRing<C>): boolean {
  // Check for RationalField signature
  if (ring.toString && ring.toString() === 'Rational Field') {
    return true;
  }
  // Duck typing: is_field returns true, characteristic returns 0
  if ('is_field' in ring && 'characteristic' in ring) {
    const r = ring as { is_field: () => boolean };
    const char = getRingCharacteristic(ring);
    if (r.is_field() && char === 0n) {
      return true;
    }
  }
  return false;
}

/**
 * Check if a ring is a finite field.
 */
function isFiniteField<C extends RingElement>(ring: CoefficientRing<C>): boolean {
  // Check for is_field method
  if (ring.is_field && !ring.is_field()) {
    return false;
  }

  // Check for cardinality or order method
  if (
    !(
      'cardinality' in ring ||
      'order' in ring ||
      'characteristic' in ring ||
      Symbol.iterator in ring
    )
  ) {
    return false;
  }

  return true;
}

/**
 * Get the order (cardinality) of a finite field.
 */
function getFieldOrder<C extends RingElement>(ring: CoefficientRing<C>): bigint {
  if (
    'cardinality' in ring &&
    typeof (ring as { cardinality: () => bigint }).cardinality === 'function'
  ) {
    return (ring as { cardinality: () => bigint }).cardinality();
  }

  if ('order' in ring) {
    const order = (ring as { order: bigint | number }).order;
    return typeof order === 'number' ? BigInt(order) : order;
  }

  if ('characteristic' in ring && 'degree' in ring) {
    const p = (ring as { characteristic: bigint }).characteristic;
    const n = (ring as { degree: number }).degree;
    return p ** BigInt(n);
  }

  throw new ValueError('cannot determine field order');
}

/**
 * Get the characteristic of the base ring.
 */
function getCharacteristic<C extends RingElement>(ring: CoefficientRing<C>): bigint {
  if ('characteristic' in ring) {
    const p = (ring as { characteristic: bigint | number }).characteristic;
    return typeof p === 'number' ? BigInt(p) : p;
  }

  // Try to determine characteristic by computing 1 + 1 + ... until we get 0
  let sum = ring.one();
  const one = ring.one();

  for (let i = 1; i < 1000; i++) {
    sum = sum.add(one) as C;
    if (sum.isZero()) {
      return BigInt(i + 1);
    }
  }

  return 0n; // Assume characteristic 0 if not found
}

/**
 * Iterate over all elements of a finite field.
 */
function* iterateField<C extends RingElement>(ring: CoefficientRing<C>): Generator<C> {
  if (Symbol.iterator in ring) {
    yield* ring as Iterable<C>;
    return;
  }

  throw new ValueError('cannot iterate over field elements');
}

/**
 * Compute the p-th root of a coefficient in a finite field.
 * In GF(q) where q = p^k, the p-th root of a is a^{q/p}.
 */
function pthRootCoeff<C extends RingElement>(coeff: C, p: bigint): C {
  if (coeff.isZero()) {
    return coeff;
  }

  // Get the field order
  if ('parent' in coeff) {
    const parent = (coeff as { parent: CoefficientRing<C> }).parent;
    const q = getFieldOrder(parent);
    const exp = q / p;

    // Use pow method if available
    if ('pow' in coeff && typeof (coeff as { pow: (n: bigint) => C }).pow === 'function') {
      return (coeff as { pow: (n: bigint) => C }).pow(exp);
    }
  }

  // Fallback: in GF(p), p-th root is identity (Frobenius inverse)
  return coeff;
}

/**
 * Compute base^exp mod modulus for polynomials.
 */
function powerMod<C extends RingElement>(
  base: Polynomial<C>,
  exp: bigint,
  modulus: Polynomial<C>
): Polynomial<C> {
  if (exp === 0n) {
    return base.parent.one();
  }

  let result = base.parent.one();
  let b = base.mod(modulus);

  while (exp > 0n) {
    if ((exp & 1n) === 1n) {
      result = result.mul(b).mod(modulus);
    }
    b = b.mul(b).mod(modulus);
    exp >>= 1n;
  }

  return result;
}

/**
 * Compute base^(q^k) mod modulus by iterating k q-th powers.
 *
 * This is FLINT's `nmod_poly_powpowmod` (`nmod_poly_factor/is_irreducible.c`).
 */
function powerModIterated<C extends RingElement>(
  base: Polynomial<C>,
  q: bigint,
  k: number,
  modulus: Polynomial<C>
): Polynomial<C> {
  let result = base.mod(modulus);
  for (let i = 0; i < k; i++) {
    result = powerMod(result, q, modulus);
  }
  return result;
}

/**
 * Cantor-Zassenhaus algorithm for equal-degree factorization.
 *
 * Given a polynomial f that is a product of distinct irreducible polynomials
 * all of the same degree d, find all irreducible factors.
 */
function cantorZassenhausFactorization<C extends RingElement>(
  f: Polynomial<C>,
  degree: number,
  ring: PolynomialRingBase<C>
): Polynomial<C>[] {
  const n = f.degree();

  if (n === 0) {
    return [];
  }

  if (n === degree) {
    // f is already irreducible
    return [f._monic()];
  }

  const baseRing = ring.base_ring;
  const q = getFieldOrder(baseRing);
  const p = getCharacteristic(baseRing);

  // We expect to succeed with probability > 1/2 per attempt, so 100 failures
  // means there is a bug (`polynomial_element.pyx:2205`).
  const maxAttempts = 100;

  for (let attempt = 0; attempt < maxAttempts; attempt++) {
    // Sample T uniformly from R of degree exactly 2*degree + 1, then make it
    // monic (`polynomial_element.pyx:2209`).
    const t = randomPolynomial(ring, 2 * degree + 1).monic();

    let h: Polynomial<C>;

    if (p === 2n) {
      // Characteristic 2: use trace
      // Compute T + T^2 + T^4 + ... + T^{2^{dk-1}} mod f
      // where k is the degree of the base field over GF(p)
      const fieldDegree = getFieldDegree(baseRing);
      const numTerms = degree * fieldDegree;

      let c = t.mod(f);
      let tt = t.mod(f);

      for (let i = 1; i < numTerms; i++) {
        tt = tt.mul(tt).mod(f); // T^{2^i}
        c = c.add(tt);
      }

      h = f.gcd(c);
    } else {
      // Odd characteristic: use (q^d - 1)/2 power
      const exponent = (q ** BigInt(degree) - 1n) / 2n;
      const tPow = powerMod(t, exponent, f);
      const tPowMinus1 = tPow.sub(ring.one());
      h = f.gcd(tPowMinus1);
    }

    const hd = h.degree();

    // Check if we found a non-trivial factor
    if (hd > 0 && hd < n) {
      // Recursively factor both parts
      const factors1 = cantorZassenhausFactorization(h._monic(), degree, ring);
      const quotient = f.quo_rem(h)[0];
      const factors2 = cantorZassenhausFactorization(quotient._monic(), degree, ring);

      return [...factors1, ...factors2];
    }
  }

  // Sage raises an AssertionError here rather than returning an unsplit
  // factor (`polynomial_element.pyx:2236`): reaching this point means the
  // input was not a product of distinct irreducibles of the given degree,
  // or that the sampler is broken.
  throw new AssertionError(`no splitting of degree ${degree} found for ${f}`);
}

/**
 * Return a uniformly random element of a (finite) coefficient ring.
 *
 * Sampling via `baseRing.__call__(someNumber)` is wrong for extension
 * fields: the number is routed through the prime subfield, so the sampled
 * element never leaves GF(p) and Cantor-Zassenhaus can never split a
 * polynomial over GF(p^k).
 */
function randomRingElement<C extends RingElement>(ring: CoefficientRing<C>): C {
  if (
    'random_element' in ring &&
    typeof (ring as { random_element: () => C }).random_element === 'function'
  ) {
    return (ring as { random_element: () => C }).random_element();
  }

  if (Symbol.iterator in ring) {
    const elements = [...(ring as unknown as Iterable<C>)];
    if (elements.length === 0) {
      throw new ValueError('cannot sample from an empty ring');
    }
    const index = Number(current_randstate().random_below(BigInt(elements.length)));
    return elements[index]!;
  }

  throw new NotImplementedError(`cannot sample a random element of ${ring}`);
}

/**
 * Generate a random polynomial of degree exactly `degree`, sampling every
 * coefficient uniformly from the base ring.
 *
 * @see Reference: sage/rings/polynomial/polynomial_ring.py:1344 (random_element)
 */
function randomPolynomial<C extends RingElement>(
  ring: PolynomialRingBase<C>,
  degree: number
): Polynomial<C> {
  const baseRing = ring.base_ring;

  if (degree < 0) {
    return ring.zero();
  }

  const coeffs: C[] = [];
  for (let i = 0; i < degree; i++) {
    coeffs.push(randomRingElement(baseRing));
  }

  // The leading coefficient must be nonzero so that the degree is exactly
  // `degree` (Sage's `random_element(d)` samples until this holds).
  let lead = randomRingElement(baseRing);
  while (lead.isZero()) {
    lead = randomRingElement(baseRing);
  }
  coeffs.push(lead);

  return new Polynomial(coeffs, ring);
}

/**
 * Get the degree of a finite field over its prime field.
 */
function getFieldDegree<C extends RingElement>(ring: CoefficientRing<C>): number {
  if ('degree' in ring) {
    return (ring as { degree: number }).degree;
  }

  // Compute from characteristic and order
  const q = getFieldOrder(ring);
  const p = getCharacteristic(ring);

  if (p === 0n) {
    return 1;
  }

  let degree = 0;
  let power = 1n;

  while (power < q) {
    power *= p;
    degree++;
  }

  return degree;
}

// ============================================
// Integer/Rational Polynomial Factorization
// ============================================

/**
 * Extract bigint coefficients from a polynomial over ZZ.
 */
function extractIntegerCoeffs<C extends RingElement>(poly: Polynomial<C>): bigint[] {
  return poly.coeffs.map((c) => {
    // Handle Integer wrapper class
    if ('value' in c && typeof (c as { value: bigint }).value === 'bigint') {
      return (c as { value: bigint }).value;
    }
    // Handle raw bigint (shouldn't happen but for safety)
    if (typeof c === 'bigint') {
      return c;
    }
    // Try toString and parse
    return BigInt(c.toString());
  });
}

/**
 * Compute the content of an integer polynomial (GCD of coefficients).
 */
function intPolyContent(coeffs: bigint[]): bigint {
  if (coeffs.length === 0) return 0n;
  let g = coeffs[0]!;
  for (let i = 1; i < coeffs.length; i++) {
    g = gcdBigInt(g, coeffs[i]!);
    if (g === 1n || g === -1n) return 1n;
  }
  return g < 0n ? -g : g;
}

/**
 * Divide all coefficients by a constant.
 */
function intPolyDivideByConstant(coeffs: bigint[], c: bigint): bigint[] {
  return coeffs.map((coeff) => coeff / c);
}

/**
 * Make polynomial primitive (divide by content).
 */
function intPolyPrimitive(coeffs: bigint[]): [bigint, bigint[]] {
  if (coeffs.length === 0) return [1n, []];
  const content = intPolyContent(coeffs);
  if (content === 0n) return [1n, coeffs];
  // Make leading coefficient positive
  const lc = coeffs[coeffs.length - 1]!;
  const sign = lc < 0n ? -1n : 1n;
  const adjustedContent = content * sign;
  return [adjustedContent, intPolyDivideByConstant(coeffs, adjustedContent)];
}

/**
 * Multiply two integer polynomials.
 */
function intPolyMul(a: bigint[], b: bigint[]): bigint[] {
  if (a.length === 0 || b.length === 0) return [];
  const result = new Array(a.length + b.length - 1).fill(0n);
  for (let i = 0; i < a.length; i++) {
    for (let j = 0; j < b.length; j++) {
      result[i + j] += a[i]! * b[j]!;
    }
  }
  return result;
}

/**
 * Compute quotient and remainder of integer polynomial division.
 * Returns [quotient, remainder] or null if division is not exact in ZZ.
 */
function intPolyQuoRem(a: bigint[], b: bigint[]): [bigint[], bigint[]] | null {
  if (b.length === 0) throw new ZeroDivisionError('polynomial division by zero');

  // Remove trailing zeros
  while (a.length > 0 && a[a.length - 1] === 0n) a = a.slice(0, -1);
  while (b.length > 0 && b[b.length - 1] === 0n) b = b.slice(0, -1);

  if (a.length < b.length) return [[0n], a];

  const degA = a.length - 1;
  const degB = b.length - 1;
  const lcB = b[degB]!;

  const quotient = new Array(degA - degB + 1).fill(0n);
  const remainder = [...a];

  for (let i = degA; i >= degB; i--) {
    if (remainder[i] === 0n) continue;

    // Check if division is exact
    if (remainder[i]! % lcB !== 0n) {
      return null; // Not exact division in ZZ
    }

    const qCoeff = remainder[i]! / lcB;
    quotient[i - degB] = qCoeff;

    for (let j = 0; j <= degB; j++) {
      remainder[i - degB + j] -= qCoeff * b[j]!;
    }
  }

  // Remove trailing zeros from remainder
  while (remainder.length > 0 && remainder[remainder.length - 1] === 0n) {
    remainder.pop();
  }

  return [quotient, remainder];
}

/**
 * Evaluate integer polynomial at a point.
 */
function intPolyEval(coeffs: bigint[], x: bigint): bigint {
  if (coeffs.length === 0) return 0n;
  let result = coeffs[coeffs.length - 1]!;
  for (let i = coeffs.length - 2; i >= 0; i--) {
    result = result * x + coeffs[i]!;
  }
  return result;
}

/**
 * Compute polynomial modulo a prime (reduce coefficients mod p).
 */
function intPolyModP(coeffs: bigint[], p: bigint): bigint[] {
  // Ensure p is bigint
  const pBig = typeof p === 'bigint' ? p : BigInt(p);
  const result = coeffs.map((c) => {
    // Ensure c is bigint
    const cBig = typeof c === 'bigint' ? c : BigInt(c);
    let r = cBig % pBig;
    if (r < 0n) r += pBig;
    return r;
  });
  // Remove trailing zeros
  while (result.length > 0 && result[result.length - 1] === 0n) {
    result.pop();
  }
  return result;
}

/**
 * Modular polynomial multiplication.
 */
function modPolyMul(a: bigint[], b: bigint[], p: bigint): bigint[] {
  if (a.length === 0 || b.length === 0) return [];
  const pBig = toBigIntSafe(p);
  const result: bigint[] = new Array(a.length + b.length - 1).fill(0n);
  for (let i = 0; i < a.length; i++) {
    for (let j = 0; j < b.length; j++) {
      const ai = toBigIntSafe(a[i]);
      const bj = toBigIntSafe(b[j]);
      result[i + j] = (result[i + j]! + ai * bj) % pBig;
    }
  }
  // Remove trailing zeros
  while (result.length > 0 && result[result.length - 1] === 0n) {
    result.pop();
  }
  return result;
}

/**
 * Modular polynomial quotient and remainder.
 */
function modPolyQuoRem(a: bigint[], b: bigint[], p: bigint): [bigint[], bigint[]] {
  if (b.length === 0) throw new ZeroDivisionError('polynomial division by zero');

  const pBig = toBigIntSafe(p);

  // Ensure all values are bigints and remove trailing zeros
  let aCopy = a.map((c) => toBigIntSafe(c));
  let bCopy = b.map((c) => toBigIntSafe(c));

  while (aCopy.length > 0 && aCopy[aCopy.length - 1] === 0n) aCopy = aCopy.slice(0, -1);
  while (bCopy.length > 0 && bCopy[bCopy.length - 1] === 0n) bCopy = bCopy.slice(0, -1);

  if (aCopy.length === 0 || aCopy.length < bCopy.length) return [[0n], aCopy];

  const degA = aCopy.length - 1;
  const degB = bCopy.length - 1;
  const lcB = bCopy[degB]!;
  const lcBInv = modInverse(lcB, pBig);

  const quotient: bigint[] = new Array(degA - degB + 1).fill(0n);
  const remainder = aCopy.map((c) => ((c % pBig) + pBig) % pBig);

  for (let i = degA; i >= degB; i--) {
    if (remainder[i] === 0n) continue;

    const qCoeff = (((remainder[i]! * lcBInv) % pBig) + pBig) % pBig;
    quotient[i - degB] = qCoeff;

    for (let j = 0; j <= degB; j++) {
      remainder[i - degB + j] =
        (((remainder[i - degB + j]! - qCoeff * bCopy[j]!) % pBig) + pBig) % pBig;
    }
  }

  // Remove trailing zeros
  while (remainder.length > 0 && remainder[remainder.length - 1] === 0n) {
    remainder.pop();
  }
  while (quotient.length > 0 && quotient[quotient.length - 1] === 0n) {
    quotient.pop();
  }

  return [quotient.length > 0 ? quotient : [0n], remainder];
}

/**
 * Convert a value to bigint safely.
 */
function toBigIntSafe(x: unknown): bigint {
  if (typeof x === 'bigint') return x;
  if (typeof x === 'number') return BigInt(Math.floor(x));
  if (typeof x === 'string') return BigInt(x);
  if (x && typeof x === 'object' && 'value' in x) {
    const v = (x as { value: unknown }).value;
    if (typeof v === 'bigint') return v;
    if (typeof v === 'number') return BigInt(Math.floor(v));
  }
  throw new Error(`Cannot convert ${typeof x} to bigint: ${x}`);
}

/**
 * Modular inverse using extended GCD.
 */
function modInverse(a: bigint, p: bigint): bigint {
  // Ensure inputs are bigints
  const aBig = toBigIntSafe(a);
  const pBig = toBigIntSafe(p);

  let [oldR, r] = [((aBig % pBig) + pBig) % pBig, pBig];
  let [oldS, s] = [1n, 0n];

  while (r !== 0n) {
    const q = oldR / r;
    [oldR, r] = [r, oldR - q * r];
    [oldS, s] = [s, oldS - q * s];
  }

  return ((oldS % pBig) + pBig) % pBig;
}

/**
 * Modular GCD of two polynomials.
 */
function modPolyGcd(a: bigint[], b: bigint[], p: bigint): bigint[] {
  const pBig = toBigIntSafe(p);
  while (b.length > 0) {
    const [_, rem] = modPolyQuoRem(a, b, pBig);
    a = b;
    b = rem;
  }
  // Make monic
  if (a.length > 0 && a[a.length - 1] !== 0n) {
    const lcInv = modInverse(a[a.length - 1]!, pBig);
    a = a.map((c) => (((toBigIntSafe(c) * lcInv) % pBig) + pBig) % pBig);
  }
  return a;
}

/**
 * Integer square root rounded down: the largest `s >= 0` with `s*s <= n`.
 *
 * This is FLINT's `fmpz_sqrt`, computed exactly (a `Math.sqrt` of a bigint both
 * overflows and rounds the wrong way).
 */
function isqrtFloor(n: bigint): bigint {
  if (n <= 0n) return 0n;
  if (n < 4n) return 1n;
  // Newton's iteration for floor(sqrt(n)), started above the root.
  let x = 1n << BigInt((n.toString(2).length >> 1) + 1);
  let y = (x + n / x) >> 1n;
  while (y < x) {
    x = y;
    y = (x + n / x) >> 1n;
  }
  return x;
}

/**
 * Mignotte's bound on the coefficients of a factor of an integer polynomial.
 *
 * Verbatim transcription of FLINT's `_fmpz_poly_factor_mignotte`
 * (`reference/flint/src/fmpz_poly_factor/factor_zassenhaus.c:47-88`), including
 * its initialisation `b = m - 1` (the comment there describes
 * `b = binomial(m-1, j-1)`, so the value FLINT actually computes is `m - 1`
 * times that; the result is therefore a valid -- merely more generous -- upper
 * bound, and we reproduce it exactly rather than "fixing" it).
 *
 * For `g` dividing `f` of degree `m`, every coefficient `b_j` of `g` satisfies
 * `|b_j| <= B`.
 */
function fmpz_poly_factor_mignotte(f: bigint[]): bigint {
  const m = f.length - 1;
  if (m < 0) return 0n;
  const abs = (x: bigint) => (x < 0n ? -x : x);

  let f2 = 0n;
  for (let j = 0; j <= m; j++) f2 += f[j]! * f[j]!;
  f2 = isqrtFloor(f2) + 1n;

  const lc = abs(f[m]!);
  let B = abs(f[0]!);

  let b = BigInt(m - 1);
  for (let j = 1; j < m; j++) {
    const t = b * lc;
    b = (b * BigInt(m - j)) / BigInt(j);
    const s = b * f2 + t;
    if (B < s) B = s;
  }

  if (B < lc) B = lc;
  return B;
}

/**
 * Reduce a polynomial's coefficients into `[0, p)` and strip trailing zeros.
 */
function modPolyNormalize(a: bigint[], p: bigint): bigint[] {
  const result = a.map((c) => ((c % p) + p) % p);
  while (result.length > 0 && result[result.length - 1] === 0n) result.pop();
  return result;
}

/**
 * Add two polynomials over Z/pZ.
 */
function modPolyAdd(a: bigint[], b: bigint[], p: bigint): bigint[] {
  const len = Math.max(a.length, b.length);
  const result: bigint[] = new Array(len).fill(0n);
  for (let i = 0; i < len; i++) {
    result[i] = ((((a[i] ?? 0n) + (b[i] ?? 0n)) % p) + p) % p;
  }
  while (result.length > 0 && result[result.length - 1] === 0n) result.pop();
  return result;
}

/**
 * Subtract two polynomials over Z/pZ.
 */
function modPolySub(a: bigint[], b: bigint[], p: bigint): bigint[] {
  const len = Math.max(a.length, b.length);
  const result: bigint[] = new Array(len).fill(0n);
  for (let i = 0; i < len; i++) {
    result[i] = ((((a[i] ?? 0n) - (b[i] ?? 0n)) % p) + p) % p;
  }
  while (result.length > 0 && result[result.length - 1] === 0n) result.pop();
  return result;
}

/**
 * Make a nonzero polynomial monic over Z/pZ.
 */
function modPolyMonic(a: bigint[], p: bigint): bigint[] {
  const f = modPolyNormalize(a, p);
  if (f.length === 0) return f;
  const lc = f[f.length - 1]!;
  if (lc === 1n) return f;
  const inv = modInverse(lc, p);
  return f.map((c) => (c * inv) % p);
}

/**
 * `base^e mod m` over Z/pZ.
 */
function modPolyPowMod(base: bigint[], e: bigint, m: bigint[], p: bigint): bigint[] {
  let result: bigint[] = [1n];
  let b = modPolyQuoRem(modPolyNormalize(base, p), m, p)[1];
  let k = e;
  while (k > 0n) {
    if ((k & 1n) === 1n) {
      result = modPolyQuoRem(modPolyMul(result, b, p), m, p)[1];
    }
    k >>= 1n;
    if (k > 0n) {
      b = modPolyQuoRem(modPolyMul(b, b, p), m, p)[1];
    }
  }
  return result;
}

/**
 * Inverse of `a` modulo `m` over Z/pZ, via the extended Euclidean algorithm.
 *
 * @throws {ArithmeticError} if `a` is not invertible modulo `m`
 */
function modPolyInvMod(a: bigint[], m: bigint[], p: bigint): bigint[] {
  let r0 = modPolyNormalize(m, p);
  let r1 = modPolyQuoRem(modPolyNormalize(a, p), m, p)[1];
  let s0: bigint[] = [];
  let s1: bigint[] = [1n];

  while (r1.length > 0) {
    const [q, rem] = modPolyQuoRem(r0, r1, p);
    const s = modPolySub(s0, modPolyMul(modPolyNormalize(q, p), s1, p), p);
    r0 = r1;
    r1 = rem;
    s0 = s1;
    s1 = s;
  }

  if (r0.length !== 1) {
    throw new ArithmeticError('polynomial is not invertible modulo the given modulus');
  }
  const inv = modInverse(r0[0]!, p);
  return modPolyNormalize(
    s0.map((c) => c * inv),
    p
  );
}

/**
 * Distinct-degree factorization of a monic squarefree polynomial over Z/pZ.
 *
 * Returns pairs `[g_d, d]` where `g_d` is the product of all monic irreducible
 * factors of degree `d`.
 *
 * @see Reference: sage/rings/polynomial/polynomial_element.pyx:_distinct_degree_factorisation_squarefree
 */
function modpDistinctDegree(f: bigint[], p: bigint): Array<[bigint[], number]> {
  const result: Array<[bigint[], number]> = [];
  let v = modPolyMonic(f, p);
  const x: bigint[] = [0n, 1n];
  // w = x^(p^d) mod v, maintained across iterations.
  let w = modPolyQuoRem(x, v, p)[1];
  let d = 0;

  while (2 * (d + 1) <= v.length - 1) {
    d += 1;
    w = modPolyPowMod(w, p, v, p);
    const g = modPolyGcd(v, modPolySub(w, x, p), p);
    if (g.length > 1) {
      result.push([g, d]);
      v = modPolyQuoRem(v, g, p)[0];
      w = modPolyQuoRem(w, v, p)[1];
    }
  }

  if (v.length > 1) {
    result.push([v, v.length - 1]);
  }

  return result;
}

/**
 * Equal-degree (Cantor-Zassenhaus) splitting over Z/pZ.
 *
 * `f` is monic, squarefree and a product of monic irreducible factors that all
 * have degree `d`.  For odd `p` the splitting element is `a^((p^d-1)/2) - 1`;
 * for `p = 2` it is the trace `a + a^2 + ... + a^(2^(d-1))` (von zur Gathen &
 * Gerhard, Algorithms 14.8/14.10 -- the same case distinction FLINT makes in
 * `nmod_poly_factor_equal_deg`).
 */
function modpEqualDegree(f: bigint[], d: number, p: bigint): bigint[][] {
  const monic = modPolyMonic(f, p);
  const n = monic.length - 1;
  if (n === d) return [monic];

  const rstate = current_randstate();
  // One trial splits with probability >= 1/2, so 512 consecutive failures do
  // not happen; we raise rather than silently drop a factor.
  const maxAttempts = 512;

  for (let attempt = 0; attempt < maxAttempts; attempt++) {
    const a: bigint[] = [];
    for (let i = 0; i < n; i++) {
      a.push(rstate.random_below(p));
    }
    const aNorm = modPolyNormalize(a, p);
    if (aNorm.length <= 1) continue;

    // A common factor with a random polynomial already splits f.
    let g = modPolyGcd(monic, aNorm, p);
    if (g.length <= 1 || g.length >= monic.length) {
      let b: bigint[];
      if (p === 2n) {
        // Trace map: a + a^2 + a^4 + ... + a^(2^(d-1)) mod f
        b = [];
        let term = modPolyQuoRem(aNorm, monic, p)[1];
        for (let i = 0; i < d; i++) {
          b = modPolyAdd(b, term, p);
          term = modPolyQuoRem(modPolyMul(term, term, p), monic, p)[1];
        }
      } else {
        const e = (p ** BigInt(d) - 1n) / 2n;
        b = modPolySub(modPolyPowMod(aNorm, e, monic, p), [1n], p);
      }
      g = modPolyGcd(monic, b, p);
    }

    if (g.length > 1 && g.length < monic.length) {
      const h = modPolyQuoRem(monic, g, p)[0];
      return [...modpEqualDegree(g, d, p), ...modpEqualDegree(h, d, p)];
    }
  }

  throw new ArithmeticError(
    `Cantor-Zassenhaus failed to split a degree-${n} polynomial modulo ${p}`
  );
}

/**
 * Factor a squarefree polynomial over Z/pZ into monic irreducible factors.
 */
function modpFactorSquarefree(f: bigint[], p: bigint): bigint[][] {
  const monic = modPolyMonic(f, p);
  if (monic.length <= 1) return [];
  if (monic.length === 2) return [monic];

  const factors: bigint[][] = [];
  for (const [g, d] of modpDistinctDegree(monic, p)) {
    if (g.length - 1 === d) {
      factors.push(g);
    } else {
      factors.push(...modpEqualDegree(g, d, p));
    }
  }
  return factors;
}

/**
 * Multifactor Hensel lifting.
 *
 * Given `f = lc(f) * h_1 * ... * h_r (mod p)` with the `h_i` monic, pairwise
 * coprime and `p` not dividing `lc(f)`, return monic `H_i = h_i (mod p)` with
 * `f = lc(f) * H_1 * ... * H_r (mod p^k)`.
 *
 * Linear lifting: writing `H_i' = H_i + p^j d_i`, the corrections satisfy
 * `sum_i d_i * prod_{l != i} h_l = (f - lc(f) * prod_i H_i) / (p^j * lc(f))`
 * modulo `p`, which is solved by `d_i = (c * s_i) mod h_i` where `s_i` is the
 * inverse of `prod_{l != i} h_l` modulo `h_i` and `c` is the right-hand side.
 *
 * @see Reference: von zur Gathen & Gerhard, "Modern Computer Algebra",
 *   Algorithm 15.17 (multifactor Hensel lifting); FLINT
 *   `fmpz_poly_factor/hensel_lift.c`.
 */
function henselLiftFactors(f: bigint[], hs: bigint[][], p: bigint, k: number): bigint[][] {
  const r = hs.length;
  const b = f[f.length - 1]!;
  const bInv = modInverse(((b % p) + p) % p, p);

  // s_i = (prod_{l != i} h_l)^-1 mod h_i, over Z/pZ
  const s: bigint[][] = [];
  for (let i = 0; i < r; i++) {
    let u: bigint[] = [1n];
    for (let l = 0; l < r; l++) {
      if (l !== i) u = modPolyQuoRem(modPolyMul(u, hs[l]!, p), hs[i]!, p)[1];
    }
    s.push(modPolyInvMod(u, hs[i]!, p));
  }

  const H = hs.map((h) => modPolyNormalize(h, p));
  let pj = p;

  for (let j = 1; j < k; j++) {
    const pj1 = pj * p;

    // e = (f - b * prod H_i) / p^j  (mod p)
    let prod: bigint[] = [b];
    for (const h of H) prod = intPolyMul(prod, h);
    const diff: bigint[] = [];
    for (let i = 0; i < Math.max(f.length, prod.length); i++) {
      let c = ((f[i] ?? 0n) - (prod[i] ?? 0n)) % pj1;
      if (c < 0n) c += pj1;
      diff.push(c);
    }

    if (diff.some((c) => c !== 0n)) {
      const e = modPolyNormalize(
        diff.map((c) => {
          if (c % pj !== 0n) {
            throw new ArithmeticError('Hensel lifting lost its invariant');
          }
          return c / pj;
        }),
        p
      );
      const c = modPolyNormalize(
        e.map((coeff) => coeff * bInv),
        p
      );

      for (let i = 0; i < r; i++) {
        const di = modPolyQuoRem(modPolyMul(c, s[i]!, p), hs[i]!, p)[1];
        const lifted = H[i]!.slice();
        for (let idx = 0; idx < di.length; idx++) {
          lifted[idx] = ((lifted[idx] ?? 0n) + pj * di[idx]!) % pj1;
        }
        H[i] = lifted;
      }
    }

    pj = pj1;
  }

  return H;
}

/* ===================================================================== *
 * Factoring in ZZ[x]: Zassenhaus recombination and van Hoeij / LLL
 *
 * Sage's `ZZ[x].factor()` delegates to FLINT's `fmpz_poly_factor`, so the code
 * from here to `factorSquarefreeIntPoly` is a transcription of FLINT, routine
 * by routine:
 *
 *   reference/flint/src/fmpz_poly_factor/factor_zassenhaus.c
 *   reference/flint/src/fmpz_poly_factor/factor_zassenhaus_recombination.c
 *   reference/flint/src/fmpz_poly_factor/zassenhaus_subset.c
 *   reference/flint/src/fmpz_poly_factor/zassenhaus_prune.c
 *   reference/flint/src/fmpz_poly_factor/factor_van_hoeij.c
 *   reference/flint/src/fmpz_poly_factor/CLD_mat.c
 *   reference/flint/src/fmpz_poly_factor/van_hoeij_check_if_solved.c
 *   reference/flint/src/fmpz_poly/CLD_bound.c
 *   reference/flint/src/fmpz_poly/divlow_smodp.c
 *   reference/flint/src/fmpz_poly/divhigh_smodp.c
 *   reference/flint/src/fmpz_mat/next_col_van_hoeij.c
 *   reference/flint/src/fmpz_mat/col_partition.c
 *   reference/flint/src/fmpz_lll/  (LLL with removal, `knapsack` variant)
 * ===================================================================== */

/** Number of bits of `|x|` (FLINT `fmpz_bits`; `0` for `x = 0`). */
function fmpzBits(x: bigint): number {
  const a = x < 0n ? -x : x;
  return a === 0n ? 0 : a.toString(2).length;
}

/** `FLINT_BIT_COUNT`: one more than the index of the highest set bit. */
function bitCountU(x: number): number {
  let n = 0;
  let v = Math.floor(x);
  while (v > 0) {
    n++;
    v = Math.floor(v / 2);
  }
  return n;
}

/** `|x|` for bigints. */
function absBig(x: bigint): bigint {
  return x < 0n ? -x : x;
}

/**
 * Symmetric remainder: the `r` with `r == a (mod m)` and `-m/2 < r <= m/2`
 * (FLINT `fmpz_smod`, `reference/flint/src/fmpz/smod.c`).
 */
function fmpzSmod(a: bigint, m: bigint): bigint {
  let r = a % m;
  if (r < 0n) r += m;
  if (r > m / 2n) r -= m;
  return r;
}

/** `fmpz_tdiv_q_2exp`: division by `2^k`, truncating towards zero. */
function tdivQ2exp(a: bigint, k: number): bigint {
  if (k <= 0) return a << BigInt(-k);
  return a < 0n ? -(-a >> BigInt(k)) : a >> BigInt(k);
}

/** Strip trailing zero coefficients (`_fmpz_poly_normalise`). */
function zpNormalise(a: bigint[]): bigint[] {
  let n = a.length;
  while (n > 0 && a[n - 1] === 0n) n--;
  return n === a.length ? a : a.slice(0, n);
}

/** Coefficientwise symmetric remainder (`fmpz_poly_scalar_smod_fmpz`). */
function zpSmod(a: bigint[], m: bigint): bigint[] {
  return zpNormalise(a.map((c) => fmpzSmod(c, m)));
}

/** `gcd(a, m)` together with `a^-1 mod m` when the gcd is 1 (`fmpz_gcdinv`). */
function gcdInvMod(a: bigint, m: bigint): [bigint, bigint] {
  let [oldR, r] = [((a % m) + m) % m, m];
  let [oldS, s] = [1n, 0n];
  while (r !== 0n) {
    const q = oldR / r;
    [oldR, r] = [r, oldR - q * r];
    [oldS, s] = [s, oldS - q * s];
  }
  return [oldR, ((oldS % m) + m) % m];
}

/* ------------------------------------------------------------------ *
 * zassenhaus_prune  (fmpz_poly_factor/zassenhaus_prune.c and
 *                    fmpz_poly_factor.h:97-145)
 * ------------------------------------------------------------------ */

interface ZassenhausPrune {
  deg: number;
  /** bit 0 = possible degree so far, bit 1 = possible for the current prime */
  posDegs: Uint8Array;
  newLength: number;
  newTotal: number;
  newDegs: number[];
}

/** `zassenhaus_prune_set_degree` (zassenhaus_prune.c:22-52). */
function zassenhausPruneSetDegree(d: number): ZassenhausPrune {
  if (d < 1) throw new ValueError('zassenhaus_prune_set_degree');
  const posDegs = new Uint8Array(d + 1);
  posDegs.fill(1);
  return { deg: d, posDegs, newLength: 0, newTotal: 0, newDegs: [] };
}

/** `zassenhaus_prune_start_add_factors` (fmpz_poly_factor.h:123-128). */
function zassenhausPruneStartAddFactors(Z: ZassenhausPrune): void {
  Z.newLength = 0;
  Z.newTotal = 0;
}

/** `zassenhaus_prune_add_factor` (zassenhaus_prune.c:54-72). */
function zassenhausPruneAddFactor(Z: ZassenhausPrune, deg: number, exp: number): void {
  if (exp < 1 || deg < 1) return;
  for (let i = 0; i < exp; i++) {
    if (Z.newLength >= Z.deg) throw new ArithmeticError('zassenhaus_prune_add_factor');
    Z.newTotal += deg;
    Z.newDegs[Z.newLength] = deg;
    Z.newLength++;
  }
}

/** `zassenhaus_prune_end_add_factors` (zassenhaus_prune.c:74-114). */
function zassenhausPruneEndAddFactors(Z: ZassenhausPrune): void {
  const a = Z.posDegs;
  const posMask = 1;
  const newMask = 2;

  if (Z.newTotal !== Z.deg) throw new ArithmeticError('zassenhaus_prune_add_factor');

  a[0] = a[0]! | newMask;
  for (let j = 1; j <= Z.deg; j++) a[j] = a[j]! & ~newMask;

  for (let i = 0; i < Z.newLength; i++) {
    const d = Z.newDegs[i]!;
    for (let j = Z.deg; j >= 0; j--) {
      if ((a[j]! & newMask) !== 0) {
        if (j + d > Z.deg) throw new ArithmeticError('zassenhaus_prune_add_factor');
        a[j + d] = a[j + d]! | newMask;
      }
    }
  }

  // merge new possibilities with old
  for (let j = 0; j <= Z.deg; j++) a[j] = a[j]! & (a[j]! >> 1);

  // 0 and deg should always be possible
  if (a[0] !== posMask || a[Z.deg] !== posMask) {
    throw new ArithmeticError('zassenhaus_prune_add_factor');
  }
}

/** `zassenhaus_prune_degree_is_possible` (fmpz_poly_factor.h:137-146). */
function zassenhausPruneDegreeIsPossible(Z: ZassenhausPrune, d: number): boolean {
  if (d <= 0) return d === 0;
  if (d >= Z.deg) return d === Z.deg;
  return Z.posDegs[d] !== 0;
}

/* ------------------------------------------------------------------ *
 * zassenhaus_subset  (fmpz_poly_factor/zassenhaus_subset.c)
 *
 * A subset of {0,...,r-1} is encoded in `s[0..r-1]`: `s[i] >= 0` means the
 * index `s[i]` is *in* the subset, `s[i] < 0` means the index `-s[i]-1` is out.
 * ------------------------------------------------------------------ */

/** `zassenhaus_subset_first` (zassenhaus_subset.c:15-27). */
function zassenhausSubsetFirst(s: number[], r: number, m: number): void {
  for (let i = 0; i < r; i++) {
    if (i >= m) s[i] = s[i]! < 0 ? s[i]! : -s[i]! - 1;
    else s[i] = s[i]! >= 0 ? s[i]! : -s[i]! - 1;
  }
}

/** `zassenhaus_subset_next` (zassenhaus_subset.c:29-59). */
function zassenhausSubsetNext(s: number[], r: number): boolean {
  let i = 0;
  while (i < r && s[i]! < 0) i++;
  const j = i;
  while (i < r && s[i]! >= 0) i++;
  const k = i;

  if (k === 0 || k >= r) return false;

  s[k] = -s[k]! - 1;
  s[k - 1] = -s[k - 1]! - 1;

  if (j > 0) {
    for (let t = 0; t < k - j - 1; t++) if (s[t]! < 0) s[t] = -s[t]! - 1;
    for (let t = k - j - 1; t < k - 1; t++) if (s[t]! >= 0) s[t] = -s[t]! - 1;
  }
  return true;
}

/** `zassenhaus_subset_next_disjoint` (zassenhaus_subset.c:61-95). */
function zassenhausSubsetNextDisjoint(s: number[], r: number): boolean {
  let total = 0;
  let last = r - 1;
  for (let i = 0; i < r; i++) {
    if (s[i]! >= 0) {
      total++;
      last = i;
    }
  }

  let j = 0;
  for (let i = 0; i < r; i++) if (s[i]! < 0) s[j++] = s[i]!;

  if (r - total < total || total < 1 || last === r - 1) return false;

  const min = Math.min(total - 1, last - total + 1);

  for (let i = 0; i < min; i++) s[i] = -s[i]! - 1;
  for (let i = last - total + 1; i < last - min + 1; i++) s[i] = -s[i]! - 1;

  return true;
}

/**
 * `_fmpz_poly_product` (factor_zassenhaus_recombination.c:18-77): the product
 * of the selected lifted factors, scaled by `leadf` and reduced into the
 * symmetric range modulo `P`.
 */
function zassenhausSubsetProduct(
  liftedFac: bigint[][],
  s: number[],
  len: number,
  P: bigint,
  leadf: bigint
): bigint[] {
  let res: bigint[] = [1n];
  for (let i = 0; i < len; i++) {
    if (s[i]! < 0) continue;
    res = zpSmod(intPolyMul(res, liftedFac[s[i]!]!), P);
  }
  return zpSmod(
    res.map((c) => c * leadf),
    P
  );
}

/**
 * `fmpz_poly_factor_zassenhaus_recombination_with_prune`
 * (factor_zassenhaus_recombination.c:156-233).
 *
 * @param liftedFac - monic factors of `F/lc(F)` modulo `P`
 * @param F - the primitive squarefree polynomial being factored
 * @param P - the lifting modulus `p^a`
 * @param Z - degree pruning data, or `null` for the unpruned variant
 */
function zassenhausRecombination(
  liftedFac: bigint[][],
  F: bigint[],
  P: bigint,
  Z: ZassenhausPrune | null
): bigint[][] {
  const r = liftedFac.length;
  const subset: number[] = [];
  for (let k = 0; k < r; k++) subset.push(k);

  const out: bigint[][] = [];
  let f = F;
  let len = r;

  for (let k = 1; k <= Math.floor(len / 2); k++) {
    zassenhausSubsetFirst(subset, len, k);
    for (;;) {
      if (Z !== null) {
        let total = 0;
        for (let i = 0; i < len; i++) {
          if (subset[i]! >= 0) total += liftedFac[subset[i]!]!.length - 1;
        }
        if (!zassenhausPruneDegreeIsPossible(Z, total)) {
          if (!zassenhausSubsetNext(subset, len)) break;
          continue;
        }
      }

      const tryme = intPolyPrimitive(
        zassenhausSubsetProduct(liftedFac, subset, len, P, f[f.length - 1]!)
      )[1];

      const div = tryme.length > 0 ? intPolyQuoRem(f, tryme) : null;
      if (div !== null && zpNormalise(div[1]).length === 0) {
        out.push(tryme);
        f = div[0];
        len -= k;
        if (!zassenhausSubsetNextDisjoint(subset, len + k)) break;
      } else {
        if (!zassenhausSubsetNext(subset, len)) break;
      }
    }
  }

  if (f.length - 1 > 0) out.push(f);
  return out;
}

/* ------------------------------------------------------------------ *
 * Double + exponent arithmetic used by fmpz_poly_CLD_bound
 * (fmpz_poly/CLD_bound.c and fmpz_poly/evaluate_horner_d_2exp.c).
 *
 * Upstream is explicitly inexact here -- the CLD bound is computed in C
 * `double`s -- so we mirror it with JavaScript numbers.  Only the SPEED of
 * van Hoeij depends on it: `max{B_1(r), B_2(r)} * N` is a valid bound for
 * every `r > 0`, so an over-estimate merely discards more columns, and the
 * factorisation itself is certified by exact trial division.
 * ------------------------------------------------------------------ */

/** `ldexp` without the overflow of a single `Math.pow(2, e)`. */
function ldexpD(m: number, e: number): number {
  let r = m;
  let k = e;
  while (k > 1000) {
    r *= 2 ** 1000;
    k -= 1000;
  }
  while (k < -1000) {
    r *= 2 ** -1000;
    k += 1000;
  }
  return r * 2 ** k;
}

/** `frexp`: `x = m * 2^e` with `|m|` in `[0.5, 1)` (and `m = 0` for `x = 0`). */
function frexpD(x: number): [number, number] {
  if (x === 0 || !Number.isFinite(x)) return [x, 0];
  const a = Math.abs(x);
  let e = Math.floor(Math.log2(a)) + 1;
  let m = ldexpD(a, -e);
  while (m >= 1) {
    m /= 2;
    e++;
  }
  while (m < 0.5) {
    m *= 2;
    e--;
  }
  return [x < 0 ? -m : m, e];
}

/** `fmpz_get_d`. */
function fmpzGetD(x: bigint): number {
  return Number(x);
}

/** `fmpz_get_d_2exp`: `x ~ m * 2^e` with `|m|` in `[0.5, 1)`. */
function fmpzGetD2exp(x: bigint): [number, number] {
  if (x === 0n) return [0, 0];
  const neg = x < 0n;
  const a = neg ? -x : x;
  const e = a.toString(2).length;
  const shift = e - 53;
  const top = shift > 0 ? a >> BigInt(shift) : a << BigInt(-shift);
  const m = Number(top) / 9007199254740992; // 2^53
  return [neg ? -m : m, e];
}

/** `fmpz_set_d`: truncation towards zero. */
function fmpzSetD(d: number): bigint {
  if (!Number.isFinite(d)) throw new ArithmeticError('CLD bound overflowed a double');
  if (Math.abs(d) < 1) return 0n;
  const neg = d < 0;
  const [m, e] = frexpD(Math.abs(d));
  const mant = BigInt(Math.floor(m * 9007199254740992)); // m * 2^53, an integer
  const sh = e - 53;
  const r = sh >= 0 ? mant << BigInt(sh) : mant >> BigInt(-sh);
  return neg ? -r : r;
}

/** `fmpz_set_d_2exp` (fmpz/set.c). */
function fmpzSetD2exp(m: number, exp: number): bigint {
  const [mm, e2] = frexpD(m);
  const e = exp + e2;
  if (e >= 53) return fmpzSetD(mm * 9007199254740992) << BigInt(e - 53);
  if (e < 0) return 0n;
  return fmpzSetD(ldexpD(mm, e));
}

/** `d_polyval` (double_extras.h:40-50). */
function dPolyval(poly: number[], x: number): number {
  const len = poly.length;
  if (len === 0) return 0;
  let t = poly[len - 1]!;
  for (let i = len - 2; i >= 0; i--) t = poly[i]! + x * t;
  return t;
}

/**
 * `_fmpz_poly_evaluate_horner_d_2exp2_precomp`
 * (fmpz_poly/evaluate_horner_d_2exp.c:114-149), including the
 * `ADJUSTMENT_DELAY = 16` delayed normalisation.
 */
function hornerD2exp(polyM: number[], polyE: number[], d: number, dexp: number): [number, number] {
  const n = polyM.length;
  if (n === 0) return [0, 0];
  if (d === 0) return [polyM[0]!, polyE[0]!];

  const ADJ = 16;
  let [xm, xe] = frexpD(d);
  xe += dexp;

  let sm = polyM[n - 1]!;
  let se = polyE[n - 1]!;

  for (let i = n - 2; i >= 0; i--) {
    // dpe_mul
    sm = sm * xm;
    se = se + xe;

    if (polyM[i] !== 0) {
      // dpe_add
      const tm = polyM[i]!;
      const te = polyE[i]!;
      const diff = se - te;
      if (sm === 0) {
        sm = tm;
        se = te;
      } else if (diff >= 0) {
        if (diff <= 53 + ADJ) sm = sm + ldexpD(tm, -diff);
      } else {
        const nd = -diff;
        if (nd > 53 + ADJ) {
          sm = tm;
          se = te;
        } else {
          sm = tm + ldexpD(sm, -nd);
          se = te;
        }
      }
    }

    if (i % ADJ === 0) {
      const [am, ae] = frexpD(sm);
      sm = am;
      se += ae;
    }
  }

  const [am, ae] = frexpD(sm);
  return [am, se + ae];
}

/** `_d_cmp_2exp` (fmpz_poly/CLD_bound.c:31-68). */
function dCmp2exp(a: number, aExp: number, b: number, bExp: number): number {
  const log2 = (n: number) => Math.log(n) / Math.log(2);

  if (aExp === 0) {
    if (bExp === 0) {
      if (a > 1.5 * b) return 2;
      if (b > 1.5 * a) return -2;
      return a >= b ? 1 : -1;
    }
    const t = 1 + Math.trunc(log2(a));
    if (t >= bExp + 2) return 2;
    if (bExp >= t + 2) return -2;
    return dCmp2exp(a / 4, 0, b * 2 ** (bExp - 2), 0);
  } else if (bExp === 0) {
    return -dCmp2exp(b, bExp, a, aExp);
  } else {
    if (aExp >= bExp + 2) return 2;
    if (bExp >= aExp + 2) return -2;
    if (aExp >= bExp) return dCmp2exp(a, aExp - bExp, b, 0);
    return -dCmp2exp(b, bExp - aExp, a, 0);
  }
}

/**
 * `fmpz_poly_CLD_bound` (fmpz_poly/CLD_bound.c:70-263).
 *
 * Given `f = a_0 + ... + a_N x^N` and `0 <= n < N`, minimise
 * `max{B_1(r), B_2(r)}` over `r > 0`, where
 * `B_1(r) = (|a_0| + ... + |a_n| r^n) / r^{n+1}` and
 * `B_2(r) = (|a_{n+1}| r^{n+1} + ... + |a_N| r^N) / r^{n+1}`,
 * and return `N` times that value: a bound for the `n`-th coefficient of the
 * "coefficient of the logarithmic derivative" `f g' / g` of any factor `g`.
 */
function fmpzPolyCLDBound(f: bigint[], n: number): bigint {
  const CLD_EPS = 0.00000001;
  const flen = f.length;

  // lo(x) = |a_n| x + |a_{n-1}| x^2 + ... + |a_0| x^{n+1}
  const loArr: bigint[] = new Array(n + 2).fill(0n);
  for (let i = 1; i <= n + 1; i++) loArr[i] = absBig(f[n + 1 - i] ?? 0n);
  const lo = zpNormalise(loArr);
  // hi(x) = |a_{n+1}| + ... + |a_N| x^{N-n-1}
  const hi = zpNormalise(f.slice(n + 1).map(absBig));

  let sizeF = 0;
  for (const c of f) sizeF = Math.max(sizeF, fmpzBits(c));

  const loD = lo.map(fmpzGetD);
  const hiD = hi.map(fmpzGetD);
  const loM: number[] = [];
  const loE: number[] = [];
  for (const c of lo) {
    const [m, e] = fmpzGetD2exp(c);
    loM.push(m);
    loE.push(e);
  }
  const hiM: number[] = [];
  const hiE: number[] = [];
  for (const c of hi) {
    const [m, e] = fmpzGetD2exp(c);
    hiM.push(m);
    hiE.push(e);
  }

  const fudge = (flen - 1) * (1.0 + flen * 2 ** -50);

  let rpow = 0.0;
  let step = 1.0;
  let rexp = 0;
  let tooMuch = false;

  // The bound is valid for every r > 0, so a bailout after many refinement
  // steps is safe; upstream's loop is unbounded.
  for (let iter = 0; iter < 5000; iter++) {
    const rbits = Math.abs(rpow);
    const maxExp = rbits * hi.length + sizeF + 1;

    let hiEval: number;
    let loEval: number;
    let hiExp: number;
    let loExp: number;

    if (rbits > 200 || rexp !== 0) {
      const r = 2 ** rpow;
      [hiEval, hiExp] = hornerD2exp(hiM, hiE, r, rexp);
      [loEval, loExp] = hornerD2exp(loM, loE, 1 / r, -rexp);
    } else if (maxExp > 950 || tooMuch) {
      const r = 2 ** rpow;
      [hiEval, hiExp] = hornerD2exp(hiM, hiE, r, 0);
      [loEval, loExp] = hornerD2exp(loM, loE, 1 / r, 0);
    } else {
      const r = 2 ** rpow;
      hiEval = dPolyval(hiD, r);
      loEval = dPolyval(loD, 1 / r);
      hiExp = 0;
      loExp = 0;
    }

    if (hiExp === 0 && loExp === 0) {
      if (1.5 * loEval < hiEval) {
        if (step >= 0.0) step = -step / 2.0;
        rpow += step;
      } else if (loEval > 1.5 * hiEval) {
        if (step < 0.0) step = -step / 2.0;
        rpow += step;
      } else if (Number.isNaN(hiEval) || Number.isNaN(loEval)) {
        tooMuch = true;
      } else {
        const maxEval = hiEval > loEval ? hiEval : loEval;
        return fmpzSetD(maxEval * fudge);
      }
    } else {
      let cmp: number;
      if (rbits > 200 || rexp !== 0) {
        const l2 = (x: number) => Math.log(x) / Math.log(2);
        if (hiExp + l2(hiEval) > 1.01 * (loExp + l2(loEval))) cmp = 2;
        else if (loExp + l2(loEval) > 1.01 * (hiExp + l2(hiEval))) cmp = -2;
        else if (hiExp + l2(hiEval) >= loExp + l2(loEval)) cmp = 1;
        else cmp = -1;
      } else {
        cmp = dCmp2exp(hiEval, hiExp, loEval, loExp);
      }

      if (Math.abs(step) < CLD_EPS) cmp = cmp === 2 ? 1 : -1;

      if (cmp === 2) {
        if (step >= 0.0) step = -step / 2.0;
        rpow += step;
      } else if (cmp === -2) {
        if (step < 0.0) step = -step / 2.0;
        rpow += step;
      } else if (cmp === 1) {
        return fmpzSetD2exp(hiEval * fudge, hiExp);
      } else {
        return fmpzSetD2exp(loEval * fudge, loExp);
      }
    }

    if (rpow > 1000.0) {
      rpow -= 1000.0;
      rexp += 1000;
    } else if (rpow < -1000.0) {
      rpow += 1000.0;
      rexp -= 1000;
    }
  }

  // Bailout: `sum |a_i|` times N is a (crude but valid) bound for r = 1.
  let s = 0n;
  for (const c of f) s += absBig(c);
  return s * BigInt(flen - 1) + 1n;
}

/**
 * `fmpz_poly_divlow_smodp` (fmpz_poly/divlow_smodp.c:15-70): the low `n`
 * coefficients of `f / g` modulo `p`, in the symmetric range.
 */
function fmpzPolyDivlowSmodp(f: bigint[], g: bigint[], p: bigint, n: number): bigint[] {
  let i = 0;
  while (g[i] === 0n) i++;
  const zeroes = i;

  const tf: bigint[] = new Array(n + zeroes).fill(0n);
  for (let j = 0; j < Math.min(f.length, n + zeroes); j++) tf[j] = f[j]!;

  const c0 = g[zeroes]! >= 0n ? g[zeroes]! : g[zeroes]! + p;
  const [d, cinv] = gcdInvMod(c0, p);
  if (d !== 1n) {
    throw new ArithmeticError('Exception (fmpz_poly_divlow_smodp). Impossible inverse.');
  }

  const res: bigint[] = new Array(n).fill(0n);
  for (let k = 0; k < n; i++, k++) {
    res[k] = fmpzSmod(tf[i]! * cinv, p);
    const m = Math.min(g.length - zeroes, n - k);
    for (let j = 0; j < m; j++) {
      tf[i + j] = fmpzSmod(tf[i + j]! - g[zeroes + j]! * res[k]!, p);
    }
  }
  return res;
}

/**
 * `fmpz_poly_divhigh_smodp` (fmpz_poly/divhigh_smodp.c:15-59): the high `n`
 * coefficients of `f / g` modulo `p`, in the symmetric range.
 */
function fmpzPolyDivhighSmodp(
  f: bigint[],
  fLen: number,
  g: bigint[],
  p: bigint,
  n: number
): bigint[] {
  const lenG = g.length;
  const tf: bigint[] = new Array(fLen).fill(0n);
  for (let j = 0; j < Math.min(f.length, fLen); j++) tf[j] = f[j]!;

  const [d, cinv] = gcdInvMod(g[lenG - 1]!, p);
  if (d !== 1n) {
    throw new ArithmeticError('Exception (fmpz_poly_divhigh_smodp). Impossible inverse.');
  }

  const res: bigint[] = new Array(n).fill(0n);
  let start = 0;
  for (let k = n - 1, i = fLen - lenG; k >= 0; i--, k--) {
    if (i < fLen - n) start++;
    res[k] = fmpzSmod(tf[i + lenG - 1]! * cinv, p);
    for (let j = start; j < lenG; j++) {
      tf[i + j] = fmpzSmod(tf[i + j]! - g[j]! * res[k]!, p);
    }
  }
  return res;
}

/**
 * `_fmpz_poly_factor_CLD_mat` (fmpz_poly_factor/CLD_mat.c:19-136).
 *
 * Returns an `(r+1) x 2k` matrix whose column `j` holds the `j`-th coefficient
 * of the logarithmic derivative `f g_i' / g_i` (mod `P`) for each lifted factor
 * `g_i`, plus the CLD bound for that column in the last row; and the number of
 * usable columns (those whose bound is small enough compared to `P`).
 */
function fmpzPolyFactorCLDMat(
  f: bigint[],
  liftedFac: bigint[][],
  P: bigint,
  k: number
): { data: bigint[][]; numDataCols: number } {
  const r = liftedFac.length;
  const bitR = Math.max(r, 20);
  const flen = f.length;

  const res: bigint[][] = [];
  for (let i = 0; i <= r; i++) res.push(new Array<bigint>(2 * k).fill(0n));

  for (let i = 0; i < k; i++) {
    res[r]![i] = fmpzPolyCLDBound(f, i);
    res[r]![2 * k - i - 1] = fmpzPolyCLDBound(f, flen - i - 2);
  }

  const bound = fmpzBits(P) - bitR - Math.floor(bitR / 2);
  const sqLen = BigInt(Math.trunc(Math.sqrt(flen)));

  let loN = 0;
  for (; loN < k; loN++) {
    if (fmpzBits(res[r]![loN]! * sqLen) > bound) break;
  }
  let hiN = 0;
  for (; hiN < k; hiN++) {
    if (fmpzBits(res[r]![2 * k - hiN - 1]! * sqLen) > bound) break;
  }

  if (loN > 0) {
    for (let i = 0; i < r; i++) {
      const fac = liftedFac[i]!;
      let zeroes = 0;
      while (fac[zeroes] === 0n) zeroes++;

      const truncLen = Math.min(fac.length, loN + zeroes + 1);
      const truncFac = fac.slice(0, truncLen);
      const gd = derivativeZ(truncFac);
      const gcld = mullowZ(f, gd, loN + zeroes);
      const col = fmpzPolyDivlowSmodp(gcld, truncFac, P, loN);
      for (let j = 0; j < loN; j++) res[i]![j] = col[j]!;
    }
  }

  if (hiN > 0) {
    const truncF = f.slice(flen - hiN);
    for (let i = 0; i < r; i++) {
      const fac = liftedFac[i]!;
      const len = fac.length - hiN - 1;
      let g: bigint[];
      if (len < 0) {
        g = new Array<bigint>(-len).fill(0n).concat(fac);
      } else {
        g = fac.slice(len);
      }
      const gd = derivativeZ(g);
      const gcld = intPolyMul(truncF, gd);
      const gcldLen = truncF.length + gd.length - 1;
      const col = fmpzPolyDivhighSmodp(gcld, gcldLen, g, P, hiN);
      for (let j = 0; j < hiN; j++) res[i]![loN + j] = col[j]!;
    }
  }

  if (hiN > 0) {
    for (let i = 0; i < hiN; i++) res[r]![loN + i] = res[r]![2 * k - hiN + i]!;
  }

  return { data: res, numDataCols: loN + hiN };
}

/** `fmpz_poly_derivative`, without stripping trailing zeros. */
function derivativeZ(a: bigint[]): bigint[] {
  const out: bigint[] = [];
  for (let i = 1; i < a.length; i++) out.push(a[i]! * BigInt(i));
  return out;
}

/** `fmpz_poly_mullow`: the product truncated to `n` coefficients. */
function mullowZ(a: bigint[], b: bigint[], n: number): bigint[] {
  const out: bigint[] = new Array(n).fill(0n);
  for (let i = 0; i < a.length && i < n; i++) {
    for (let j = 0; j < b.length && i + j < n; j++) {
      out[i + j] = out[i + j]! + a[i]! * b[j]!;
    }
  }
  return out;
}

/**
 * `fmpz_mat_next_col_van_hoeij` (fmpz_mat/next_col_van_hoeij.c:35-101).
 *
 * Appends the scaled knapsack column `col` (with the new relation row) to `M`.
 * Returns `false` when running LLL would not be justified yet.
 */
function fmpzMatNextColVanHoeij(
  M: bigint[][],
  P: bigint,
  col: bigint[],
  exp: number,
  uExp: number
): boolean {
  const r = col.length;
  const bitR = Math.max(r, 20);
  const s = M.length;

  let k = fmpzBits(P) - bitR - Math.floor(bitR / 2);

  // check if LLL justified
  if (k < exp + bitCountU(r + 1)) return false;

  k -= uExp; // we want this many bits beyond the radix point

  let x: bigint[];
  let pTrunc: bigint;
  if (k >= 0) {
    x = col.map((c) => tdivQ2exp(c, k));
    pTrunc = tdivQ2exp(P, k);
  } else {
    x = col.map((c) => c << BigInt(-k));
    pTrunc = P << BigInt(-k);
  }

  // y = U x, where U is the combinatorial (first r columns) part of M
  const y: bigint[] = new Array(s).fill(0n);
  for (let j = 0; j < s; j++) {
    let acc = 0n;
    const row = M[j]!;
    for (let i = 0; i < r; i++) acc += row[i]! * x[i]!;
    y[j] = fmpzSmod(tdivQ2exp(acc, uExp), pTrunc);
  }

  // resize M: a new zero row on top, and a new column
  const c = M[0]!.length;
  for (let j = s - 1; j >= 0; j--) M[j + 1] = M[j]!.concat([y[j]!]);
  M[0] = new Array<bigint>(c).fill(0n).concat([pTrunc]);

  return true;
}

/**
 * LLL reduction of the rows of `B` with removal, `delta = 0.99`.
 *
 * This is `fmpz_lll_wrapper_with_removal_knapsack`
 * (fmpz_lll/wrapper_with_removal_knapsack.c:17-40).  FLINT reduces with a
 * chain of floating-point implementations (`fmpz_lll_d_with_removal_knapsack`,
 * `fmpz_lll_d_heuristic_with_removal`, `fmpz_lll_mpf_with_removal`) whose
 * output it then *verifies* against the exact rational predicate
 * `fmpz_mat_is_reduced_with_removal` / `gr_mat_is_row_lll_reduced_with_removal_naive`
 * (fmpz_mat/is_reduced_with_removal.c:18-70, gr_mat/is_lll_reduced.c:18-100).
 * We implement that exact contract directly with the integral Gram-Schmidt LLL
 * (Cohen, *A Course in Computational Algebraic Number Theory*, Algorithm 2.6.7)
 * at `delta = fl->delta = 0.99`, so the result satisfies the predicate FLINT
 * checks by construction.
 *
 * @see Deviation: Polynomial Roots and Factorization
 *
 * The removal rule is upstream's (`fmpz_lll/lll_d.c:392-406`): scanning from
 * the last row backwards, a row is dropped while `||b_i^*||^2 / 2 > gs_B`.
 * Since that implies `||b_i^*||^2 > gs_B`, no lattice vector of squared norm
 * at most `gs_B` can involve the dropped rows, which is exactly the property
 * van Hoeij needs.
 *
 * @param B - the basis, one lattice vector per row; reduced in place
 * @param gsB - the removal bound
 * @returns the number of rows to keep
 */
function lllWithRemovalKnapsack(B: bigint[][], gsB: bigint): number {
  const d = B.length;
  if (d === 0) return 0;

  const dot = (u: bigint[], v: bigint[]): bigint => {
    let s = 0n;
    for (let i = 0; i < u.length; i++) s += u[i]! * v[i]!;
    return s;
  };

  const dd: bigint[] = new Array(d + 1).fill(0n);
  dd[0] = 1n;
  const lam: bigint[][] = [];
  for (let i = 0; i < d; i++) lam.push(new Array<bigint>(d).fill(0n));

  dd[1] = dot(B[0]!, B[0]!);
  if (dd[1] === 0n) {
    throw new ArithmeticError('LLL: the van Hoeij lattice basis is linearly dependent');
  }

  // round to nearest integer (ties away from -infinity), exactly
  const roundDiv = (a: bigint, b: bigint): bigint => {
    const num = 2n * a + b;
    const den = 2n * b;
    let q = num / den;
    if (num % den !== 0n && num < 0n !== den < 0n) q -= 1n;
    return q;
  };

  const RED = (kk: number, l: number): void => {
    const lkl = lam[kk]![l]!;
    const dl1 = dd[l + 1]!;
    if (2n * absBig(lkl) <= dl1) return;
    const q = roundDiv(lkl, dl1);
    const bk = B[kk]!;
    const bl = B[l]!;
    for (let i = 0; i < bk.length; i++) bk[i] = bk[i]! - q * bl[i]!;
    lam[kk]![l] = lkl - q * dl1;
    for (let i = 0; i < l; i++) lam[kk]![i] = lam[kk]![i]! - q * lam[l]![i]!;
  };

  let k = 1;
  let kmax = 0;

  while (k < d) {
    if (k > kmax) {
      kmax = k;
      for (let j = 0; j <= k; j++) {
        let u = dot(B[k]!, B[j]!);
        for (let i = 0; i < j; i++) {
          u = (dd[i + 1]! * u - lam[k]![i]! * lam[j]![i]!) / dd[i]!;
        }
        if (j < k) lam[k]![j] = u;
        else dd[k + 1] = u;
      }
      if (dd[k + 1] === 0n) {
        throw new ArithmeticError('LLL: the van Hoeij lattice basis is linearly dependent');
      }
    }

    RED(k, k - 1);

    // Lovasz condition with delta = 99/100
    const lamk = lam[k]![k - 1]!;
    if (100n * dd[k + 1]! * dd[k - 1]! < 99n * dd[k]! * dd[k]! - 100n * lamk * lamk) {
      // SWAPI(k)
      const tmpRow = B[k]!;
      B[k] = B[k - 1]!;
      B[k - 1] = tmpRow;
      for (let j = 0; j <= k - 2; j++) {
        const t = lam[k]![j]!;
        lam[k]![j] = lam[k - 1]![j]!;
        lam[k - 1]![j] = t;
      }
      const BB = (dd[k - 1]! * dd[k + 1]! + lamk * lamk) / dd[k]!;
      for (let i = k + 1; i <= kmax; i++) {
        const t = lam[i]![k]!;
        lam[i]![k] = (dd[k + 1]! * lam[i]![k - 1]! - lamk * t) / dd[k]!;
        lam[i]![k - 1] = (BB * t + lamk * lam[i]![k]!) / dd[k + 1]!;
      }
      dd[k] = BB;
      k = Math.max(1, k - 1);
    } else {
      for (let l = k - 2; l >= 0; l--) RED(k, l);
      k++;
    }
  }

  // removal: ||b_i^*||^2 = dd[i+1]/dd[i]
  let newd = d;
  for (let i = d - 1; i >= 0; i--) {
    if (dd[i + 1]! > 2n * gsB * dd[i]!) newd--;
    else break;
  }
  return newd;
}

/**
 * `fmpz_mat_col_partition` (fmpz_mat/col_partition.c:62-133) with
 * `short_circuit = 1`: partition the columns of `M` into classes of equal
 * columns, numbered from 1; return the number of classes, or `0` if there are
 * more classes than `M` has rows.
 *
 * Upstream numbers the classes in the order produced by sorting a cheap hash of
 * each column; we number them by first occurrence.  Only the numbering differs
 * -- the partition itself, and hence every trial factor built from it, is the
 * same (the trial factors are sorted by degree straight afterwards).
 */
function fmpzMatColPartition(part: number[], M: bigint[][], nrows: number): number {
  const c = part.length;
  const keys: string[] = [];
  for (let j = 0; j < c; j++) {
    const col: string[] = [];
    for (let i = 0; i < M.length; i++) col.push(M[i]![j]!.toString(36));
    keys.push(col.join(','));
  }

  const seen = new Map<string, number>();
  let p = 0;
  for (let j = 0; j < c; j++) {
    let id = seen.get(keys[j]!);
    if (id === undefined) {
      p++;
      if (p > nrows) return 0;
      id = p;
      seen.set(keys[j]!, id);
    }
    part[j] = id;
  }
  return p;
}

/**
 * `fmpz_poly_factor_van_hoeij_check_if_solved`
 * (fmpz_poly_factor/van_hoeij_check_if_solved.c:29-142).
 *
 * If the combinatorial part of `M` already describes a 0-1 basis of the true
 * factor combinations, build the factors and verify them by exact division;
 * otherwise return `null`.
 */
function vanHoeijCheckIfSolved(
  M: bigint[][],
  liftedFac: bigint[][],
  f: bigint[],
  P: bigint,
  lc: bigint
): bigint[][] | null {
  const r = liftedFac.length;
  const U = M.map((row) => row.slice(0, r));
  const part: number[] = new Array(r).fill(0);

  const numFacs = fmpzMatColPartition(part, U, M.length);
  if (numFacs === 0 || numFacs > r) return null;

  if (numFacs === 1) {
    // f is irreducible
    return [f];
  }

  // there is a potential 0-1 basis, so make the potential factors
  const trial: bigint[][] = [];
  let tempLc = lc;
  for (let i = 1; i <= numFacs; i++) {
    let prod: bigint[] = [tempLc];
    for (let j = 0; j < r; j++) {
      if (part[j] === i) prod = zpSmod(intPolyMul(prod, liftedFac[j]!), P);
    }
    if (prod.length === 0) return null;
    tempLc = absBig(intPolyContent(prod));
    if (tempLc === 0n) return null;
    trial.push(prod.map((cc) => cc / tempLc));
  }

  // sort factors by length
  trial.sort((a, b) => a.length - b.length);

  // trial divide potential factors
  let fCopy = f;
  let remaining = numFacs;
  let i = 0;
  for (; i < trial.length && remaining > 1; i++) {
    const div = intPolyQuoRem(fCopy, trial[i]!);
    if (div !== null && zpNormalise(div[1]).length === 0) {
      fCopy = div[0];
      remaining--;
    } else {
      return null;
    }
  }

  if (remaining === 1) {
    const out = trial.slice(0, i);
    out.push(fCopy);
    return out;
  }
  return null;
}

/**
 * `_heuristic_van_hoeij_starting_precision` (factor_van_hoeij.c:24-42).
 */
function heuristicVanHoeijStartingPrecision(f: bigint[], r: number, p: bigint): number {
  const leadB = fmpzPolyCLDBound(f, f.length - 2);
  const trailB = fmpzPolyCLDBound(f, 0);
  const minB = Math.min(fmpzBits(leadB), fmpzBits(trailB));
  // C truncates the inner expression to slong before dividing by log(p)
  const inner = Math.trunc((2.5 * r + minB) * Math.LN2 + Math.log(f.length) / 2.0);
  return Math.trunc(inner / Math.log(Number(p)));
}

/** Smallest `a` with `p^a >= B` (`fmpz_clog_ui`). */
function clogUi(B: bigint, p: bigint): number {
  if (B <= 1n) return 0;
  let a = 0;
  let q = 1n;
  while (q < B) {
    q *= p;
    a++;
  }
  return a;
}

/**
 * `fmpz_poly_factor_van_hoeij` (fmpz_poly_factor/factor_van_hoeij.c:62-232).
 *
 * Recombination by LLL on the knapsack lattice of the lifted modular factors,
 * which is polynomial time where plain Zassenhaus is exponential.
 *
 * @param fac - the monic irreducible factors of `f` modulo `p`
 * @param f - the primitive squarefree polynomial to factor
 * @param p - the factorisation prime
 */
function fmpzPolyFactorVanHoeij(fac: bigint[][], f: bigint[], p: bigint): bigint[][] {
  const r = fac.length;
  const bitR = Math.max(r, 20);

  // set to identity, prescaled by 2^U_exp
  const uExp = bitCountU(bitR);
  let M: bigint[][] = [];
  for (let i = 0; i < r; i++) {
    const row = new Array<bigint>(r).fill(0n);
    row[i] = 1n << BigInt(uExp);
    M.push(row);
  }

  // compute Mignotte bound
  let B = fmpz_poly_factor_mignotte(f);
  B = absBig(B * f[f.length - 1]!) * 2n + 1n;
  let a = clogUi(B, p);

  // compute heuristic starting precision
  a = Math.min(a, heuristicVanHoeijStartingPrecision(f, r, p));
  if (a < 1) a = 1;

  let lifted = liftModularFactors(f, fac, p, a);

  // compute bound
  const gsB = BigInt(r + 1) << BigInt(2 * uExp);

  const N = f.length - 1;
  const sqN = BigInt(Math.trunc(Math.sqrt(N)));
  const lc = f[N]!;

  let henselLoops = 0;
  let P = p ** BigInt(a);

  for (;;) {
    const solved = vanHoeijCheckIfSolved(M, lifted, f, P, lc);
    if (solved !== null) return solved;

    let numCoeffs: number;
    if (henselLoops < 3 && 3 * r > N + 1) numCoeffs = r > 200 ? 50 : 30;
    else numCoeffs = 10;

    numCoeffs = Math.min(numCoeffs, Math.floor((N + 1) / 2));
    let prevNumCoeffs = 0;

    do {
      const { data, numDataCols } = fmpzPolyFactorCLDMat(f, lifted, P, numCoeffs);

      for (let nextCol = prevNumCoeffs; nextCol < numDataCols - prevNumCoeffs; nextCol++) {
        // we alternate taking columns from the right and left
        const diff = nextCol - prevNumCoeffs;
        const altCol =
          diff % 2 === 0
            ? prevNumCoeffs + diff / 2
            : numDataCols - prevNumCoeffs - Math.floor((diff + 1) / 2);

        const boundSum = data[r]![altCol]! * sqN;
        const worstExp = fmpzBits(boundSum);

        const col: bigint[] = [];
        for (let i = 0; i < r; i++) col.push(data[i]![altCol]!);

        if (fmpzMatNextColVanHoeij(M, P, col, worstExp, uExp)) {
          const numRows = lllWithRemovalKnapsack(M, gsB);
          M = M.slice(0, numRows);

          const done = vanHoeijCheckIfSolved(M, lifted, f, P, lc);
          if (done !== null) return done;
        }
      }

      prevNumCoeffs = numCoeffs;
      numCoeffs = Math.min(2 * numCoeffs, Math.floor((N + 1) / 2));
    } while (numCoeffs !== prevNumCoeffs);

    henselLoops++;
    if (henselLoops > 32) {
      throw new ArithmeticError(
        'van Hoeij recombination failed to converge (fmpz_poly_factor_van_hoeij)'
      );
    }

    a = 2 * a;
    P = p ** BigInt(a);
    lifted = liftModularFactors(f, fac, p, a);
  }
}

/**
 * Hensel-lift the modular factorisation of `f` to `p^a` and put the factors in
 * the symmetric range, matching what `_fmpz_poly_hensel_start_lift`
 * (fmpz_poly/hensel_start_lift.c:19-104) leaves in `lifted_fac`: monic `H_i`
 * with `H_1 ... H_r = f / lc(f) (mod p^a)`.
 */
function liftModularFactors(f: bigint[], fac: bigint[][], p: bigint, a: number): bigint[][] {
  const P = p ** BigInt(a);
  return henselLiftFactors(f, fac, p, a).map((h) => {
    const s = h.map((c) => fmpzSmod(c, P));
    s[s.length - 1] = 1n; // the lifted factors are monic
    return s;
  });
}

/**
 * Choose the factorisation prime and factor `f` modulo it.
 *
 * Transcription of the prime search in `_fmpz_poly_factor_zassenhaus`
 * (fmpz_poly_factor/factor_zassenhaus.c:113-166): three rounds over the primes,
 * each starting where the previous stopped, keeping the factorisation with the
 * fewest factors.  `p` must not divide the leading *or* the constant
 * coefficient, and `f` must stay squarefree modulo `p`.
 *
 * The search is unbounded, exactly as upstream's `for ( ; ; p = n_nextprime(p, 0))`
 * is: for a squarefree `f` with `f(0) != 0` all but finitely many primes work
 * (those not dividing `lc(f) f(0) disc(f)`), so it always terminates.
 */
function chooseFactorizationPrime(coeffs: bigint[]): {
  p: bigint;
  fac: bigint[][];
  Z: ZassenhausPrune;
} {
  const lenF = coeffs.length;
  const Z = zassenhausPruneSetDegree(lenF - 1);

  let r = lenF;
  let p = 2n;
  let bestP = 0n;
  let bestFac: bigint[][] = [];

  // A prime is rejected only when it divides `lc(f) * f(0) * disc(f)`, so a
  // squarefree `f` with `f(0) != 0` has at most `bits(lc f(0) disc(f))` bad
  // primes.  This is a (very generous) bound on that count, used only so that a
  // caller who violates the precondition gets an error instead of a hang;
  // upstream's loop (factor_zassenhaus.c:120) has no bound at all.
  let maxBits = 0;
  for (const c of coeffs) maxBits = Math.max(maxBits, fmpzBits(c));
  const badLimit = 1000 + 4 * lenF * (maxBits + 10);
  let bad = 0;

  for (let i = 0; i < 3; i++) {
    for (; ; p = next_prime(p)) {
      if (bad++ > badLimit) {
        throw new ValueError(
          'no factorization prime found: the polynomial is not squarefree, or f(0) = 0'
        );
      }
      const t = intPolyModP(coeffs, p);
      if (t.length === lenF && t[0] !== 0n) {
        const d = modPolyNormalize(
          t.slice(1).map((c, j) => c * BigInt(j + 1)),
          p
        );
        const g = d.length === 0 ? [0n] : modPolyGcd(t, d, p);

        if (g.length === 1 && g[0] === 1n) {
          const tempFac = modpFactorSquarefree(t, p);

          zassenhausPruneStartAddFactors(Z);
          for (const h of tempFac) zassenhausPruneAddFactor(Z, h.length - 1, 1);
          zassenhausPruneEndAddFactors(Z);

          if (tempFac.length <= r) {
            r = tempFac.length;
            bestFac = tempFac;
            bestP = p;
          }
          break;
        }
      }
    }
    p = next_prime(p);
  }

  return { p: bestP, fac: bestFac, Z };
}

/**
 * Factor a primitive squarefree integer polynomial into irreducible factors.
 *
 * ALGORITHM: `_fmpz_poly_factor_zassenhaus` with `cutoff = 8` and
 * `use_van_hoeij = 1` -- which is exactly how `fmpz_poly_factor` (the routine
 * behind Sage's `ZZ[x].factor()`) calls it, see
 * `reference/flint/src/fmpz_poly_factor/factor.c:98-104` and
 * `reference/flint/src/fmpz_poly_factor/factor_zassenhaus.c:90-210`.
 *
 * Factor `f` modulo a well-chosen prime `p`, then:
 *   - `r = 1`: `f` is irreducible over `Z`;
 *   - `r <= 8`: Hensel-lift to `p^a` with `p^a > 2|lc B_Mignotte| + 1` and run
 *     subset recombination with degree pruning;
 *   - `r > 8`: van Hoeij, i.e. LLL on the knapsack lattice of the lifted
 *     factors, which is polynomial rather than exponential in `r`.
 *
 * @param coeffs - Primitive squarefree polynomial coefficients, constant first
 * @returns Irreducible factors, each primitive with positive leading
 *   coefficient, whose product is the primitive part of `coeffs`
 */
function factorSquarefreeIntPoly(coeffs: bigint[]): bigint[][] {
  const n = coeffs.length - 1; // degree
  if (n <= 0) return coeffs.length > 0 && coeffs[0] !== 0n ? [[coeffs[0]!]] : [];
  if (n === 1) return [coeffs];

  // `fmpz_poly_factor` strips the x^k part before the squarefree
  // decomposition (factor.c:52-63); a squarefree f can only have x^1, and the
  // prime search below needs f(0) != 0.
  if (coeffs[0] === 0n) {
    const rest = factorSquarefreeIntPoly(coeffs.slice(1));
    return [[0n, 1n], ...rest];
  }

  const { p, fac, Z } = chooseFactorizationPrime(coeffs);
  const r = fac.length;
  const cutoff = 8;

  let factors: bigint[][];

  if (r === 1 && r <= cutoff) {
    // irreducible modulo p, hence irreducible over Z
    factors = [coeffs];
  } else if (r > cutoff) {
    factors = fmpzPolyFactorVanHoeij(fac, coeffs, p);
  } else {
    // bound adjustment: we multiply true factors (which might be monic) by the
    // leading coefficient of f (factor_zassenhaus.c:186-196)
    let T = fmpz_poly_factor_mignotte(coeffs);
    T = absBig(T * coeffs[n]!) * 2n + 1n;
    const a = clogUi(T, p);

    const lifted = liftModularFactors(coeffs, fac, p, a);
    const P = p ** BigInt(a);

    factors = zassenhausRecombination(lifted, coeffs, P, Z);
  }

  // Normalise: primitive, positive leading coefficient.  (FLINT keeps whatever
  // sign the recombination produced; the number of sign flips is even, so the
  // product is unchanged.)
  factors = factors.map((g) => intPolyPrimitive(g)[1]);

  // The factors must reproduce the input exactly; silently returning a wrong
  // factorization would be worse than raising.
  let check: bigint[] = [1n];
  for (const g of factors) check = intPolyMul(check, g);
  const checkPrimitive = intPolyPrimitive(check)[1];
  const inputPrimitive = intPolyPrimitive(coeffs)[1];
  if (
    checkPrimitive.length !== inputPrimitive.length ||
    checkPrimitive.some((c, i) => c !== inputPrimitive[i])
  ) {
    throw new ArithmeticError('integer polynomial factorization failed to reproduce its input');
  }

  return factors.length > 0 ? factors : [coeffs];
}

/**
 * Compute squarefree factorization of an integer polynomial.
 * Returns pairs [squarefree_factor, multiplicity].
 */
function squarefreeFactorIntPoly(coeffs: bigint[]): Array<[bigint[], number]> {
  if (coeffs.length === 0) return [];
  if (coeffs.length === 1) return [[coeffs, 1]];

  // Make primitive
  const [content, primitive] = intPolyPrimitive(coeffs);

  // For simple cases, just return the polynomial
  if (primitive.length <= 2) {
    return [[primitive, 1]];
  }

  // Compute derivative
  const deriv: bigint[] = [];
  for (let i = 1; i < primitive.length; i++) {
    deriv.push(primitive[i]! * BigInt(i));
  }

  // If derivative is zero or constant, polynomial is already squarefree (or p-th power in char p)
  if (deriv.length === 0 || (deriv.length === 1 && deriv[0] === 0n)) {
    return [[primitive, 1]];
  }

  // gcd(f, f')
  let g = intPolyGcd(primitive, deriv);

  // If gcd is constant, f is squarefree
  if (g.length <= 1) {
    return [[primitive, 1]];
  }

  // f / gcd(f, f')
  const divResult = intPolyQuoRem(primitive, g);
  if (divResult === null) {
    // Shouldn't happen for well-formed input
    return [[primitive, 1]];
  }
  let h = divResult[0];

  const result: Array<[bigint[], number]> = [];
  let i = 1;
  // Every multiplicity is at most the degree, so the loop cannot run longer
  // than that; an earlier fixed cap of 20 silently truncated the
  // decomposition of things like (x-1)^25, which then reached the Zassenhaus
  // code with a non-squarefree input.
  const maxIter = primitive.length;

  while (h.length > 1 && i <= maxIter) {
    // gcd(g, h)
    const gi = intPolyGcd(g, h);

    // h / gi
    const hDivGi = intPolyQuoRem(h, gi);
    // Both are primitive and gi divides h in Q[x], so by Gauss's lemma the
    // division is exact over Z; a failure means a broken invariant, and
    // continuing would silently return a wrong decomposition.
    if (hDivGi === null) {
      throw new ArithmeticError('squarefree decomposition: inexact division over ZZ');
    }
    const hi = hDivGi[0];

    if (hi.length > 1) {
      result.push([hi, i]);
    }

    // g = g / gi
    const gDivGi = intPolyQuoRem(g, gi);
    if (gDivGi === null) {
      throw new ArithmeticError('squarefree decomposition: inexact division over ZZ');
    }
    g = gDivGi[0];
    h = gi;
    i++;
  }

  // g might still have content
  if (g.length > 1) {
    result.push([g, i]);
  }

  // If no factors found, return original as squarefree
  if (result.length === 0) {
    return [[primitive, 1]];
  }

  return result;
}

/**
 * GCD of the *primitive parts* of two integer polynomials, using the
 * primitive PRS (pseudo-remainder sequence).
 *
 * The result is primitive with a positive leading coefficient.
 */
function intPolyGcd(a: bigint[], b: bigint[]): bigint[] {
  // Remove trailing zeros
  while (a.length > 0 && a[a.length - 1] === 0n) a = a.slice(0, -1);
  while (b.length > 0 && b[b.length - 1] === 0n) b = b.slice(0, -1);

  if (b.length === 0) return a.length > 0 ? intPolyPrimitive(a)[1] : [1n];
  if (a.length === 0) return b.length > 0 ? intPolyPrimitive(b)[1] : [1n];
  if (a.length < b.length) [a, b] = [b, a];

  // Primitive PRS: deg(b) strictly decreases at every step, so this
  // terminates after at most deg(a) iterations.
  while (b.length > 0) {
    const [_, rem] = pseudoDivide(a, b);
    if (rem.length === 0) {
      // b divides a exactly: b is the gcd (returning `a` here dropped one
      // Euclid step and produced a *multiple* of the gcd).
      a = b;
      break;
    }
    // Make primitive to avoid coefficient explosion
    const [__, primRem] = intPolyPrimitive(rem);
    if (primRem.length === 0) {
      a = b;
      break;
    }
    a = b;
    b = primRem;
  }

  // Make primitive and positive leading coefficient
  const [_, primA] = intPolyPrimitive(a);
  return primA;
}

/**
 * GCD of two integer polynomials in ZZ[x], i.e. including the content:
 * `gcd(f, g) = gcd(cont(f), cont(g)) * gcd(pp(f), pp(g))`.
 *
 * The result has a positive leading coefficient, matching FLINT's
 * `fmpz_poly_gcd` (which is what Sage's `ZZ[x].gcd` delegates to).
 */
function intPolyGcdWithContent(a: bigint[], b: bigint[]): bigint[] {
  while (a.length > 0 && a[a.length - 1] === 0n) a = a.slice(0, -1);
  while (b.length > 0 && b[b.length - 1] === 0n) b = b.slice(0, -1);

  if (a.length === 0 && b.length === 0) return [];
  if (a.length === 0) return intPolyPrimitive(b)[1].map((c) => c * intPolyContent(b));
  if (b.length === 0) return intPolyPrimitive(a)[1].map((c) => c * intPolyContent(a));

  const contentGcd = gcdBigInt(intPolyContent(a), intPolyContent(b));
  const primitiveGcd = intPolyGcd(a, b);
  return primitiveGcd.map((c) => c * contentGcd);
}

/**
 * Pseudo-division: compute q, r such that b_n^{m-n+1} * a = q * b + r
 * where b_n is the leading coefficient of b.
 */
function pseudoDivide(a: bigint[], b: bigint[]): [bigint[], bigint[]] {
  if (b.length === 0) throw new ZeroDivisionError('division by zero');

  const m = a.length - 1;
  const n = b.length - 1;

  if (m < n) return [[0n], a];

  const bn = b[n]!;
  const d = m - n;

  let r = [...a];
  let q = new Array(d + 1).fill(0n);

  for (let i = m; i >= n; i--) {
    // Invariant: bn^k * a = q*b + r after k completed iterations.  Each
    // iteration replaces r by bn*r - qCoeff*x^(i-n)*b, so the quotient
    // accumulated so far must be scaled by bn as well.
    q = q.map((c) => c * bn);

    if (r[i] === undefined || r[i] === 0n) {
      // Multiply r by bn
      r = r.map((c) => c * bn);
      continue;
    }

    const qCoeff = r[i]!;
    q[i - n] = qCoeff;

    // r = bn * r - qCoeff * x^{i-n} * b
    for (let j = 0; j <= i; j++) {
      if (j >= i - n && j <= i) {
        r[j] = bn * r[j]! - qCoeff * (b[j - (i - n)] || 0n);
      } else {
        r[j] = bn * r[j]!;
      }
    }
  }

  // Remove trailing zeros
  while (r.length > 0 && r[r.length - 1] === 0n) r.pop();

  return [q, r];
}

/**
 * Factor an integer polynomial completely: squarefree decomposition followed by
 * Zassenhaus on each squarefree part.
 *
 * Returns `[content, factors]` where `factors` is an array of
 * `[irreducible_factor, multiplicity]` and `content * prod(factors^mult)` is
 * the input.
 *
 * @see Deviation: Polynomial Roots and Factorization
 */
function factorIntegerPolynomial(coeffs: bigint[]): [bigint, Array<[bigint[], number]>] {
  if (coeffs.length === 0) return [0n, []];

  // Extract content and make primitive
  const [content, primitive] = intPolyPrimitive(coeffs);

  if (primitive.length <= 1) {
    return [content, []];
  }

  // Squarefree factorization
  const sqfree = squarefreeFactorIntPoly(primitive);

  // Factor each squarefree part
  const result: Array<[bigint[], number]> = [];

  for (const [sqfFactor, mult] of sqfree) {
    if (sqfFactor.length <= 1) continue;

    const irredFactors = factorSquarefreeIntPoly(sqfFactor);
    for (const irredFactor of irredFactors) {
      result.push([irredFactor, mult]);
    }
  }

  return [content, result];
}

/**
 * Find integer roots of a polynomial using rational root theorem.
 * A rational root p/q must have p dividing the constant term and q dividing the leading coefficient.
 * For integer roots, q = 1, so we only need p dividing the constant term.
 */
function findIntegerRoots(coeffs: bigint[]): Array<[bigint, number]> {
  if (coeffs.length === 0) return [];
  if (coeffs.length === 1) return []; // Constant polynomial has no roots

  const constant = coeffs[0]!;
  if (constant === 0n) {
    // 0 is a root, find its multiplicity
    let mult = 0;
    let f = coeffs;
    while (f.length > 0 && f[0] === 0n) {
      mult++;
      f = f.slice(1);
    }
    const roots: Array<[bigint, number]> = [[0n, mult]];

    // Recursively find other roots
    if (f.length > 1) {
      roots.push(...findIntegerRoots(f));
    }
    return roots;
  }

  // Get divisors of constant term
  const divisors = getDivisorsBigInt(constant < 0n ? -constant : constant);
  const roots: Array<[bigint, number]> = [];

  let f = coeffs;

  for (const d of divisors) {
    // Try both d and -d
    for (const candidate of [d, -d]) {
      if (intPolyEval(f, candidate) === 0n) {
        // Found a root, find its multiplicity
        let mult = 0;
        const linearFactor = [-candidate, 1n];
        let divResult = intPolyQuoRem(f, linearFactor);

        while (
          divResult !== null &&
          (divResult[1].length === 0 || divResult[1].every((c) => c === 0n))
        ) {
          mult++;
          f = divResult[0];
          divResult = intPolyQuoRem(f, linearFactor);
        }

        if (mult > 0) {
          roots.push([candidate, mult]);
        }
      }
    }
  }

  return roots;
}

/**
 * Get all positive divisors of a positive integer.
 *
 * Uses the prime factorization (which delegates to PARI) instead of
 * trial dividing up to sqrt(n): the latter made `roots()` over ZZ/QQ take
 * Theta(sqrt(|a_0|)) time, e.g. 10 s for `roots(x - 10^17)`.
 *
 * @see Reference: sage/arith/misc.py:divisors
 */
function getDivisorsBigInt(n: bigint): bigint[] {
  if (n <= 0n) return [];
  if (n === 1n) return [1n];

  let divisors: bigint[] = [1n];
  for (const [p, e] of factorInteger(n)) {
    if (p <= 1n) continue;
    const next: bigint[] = [];
    let pk = 1n;
    for (let i = 0n; i <= e; i++) {
      for (const d of divisors) {
        next.push(d * pk);
      }
      pk *= p;
    }
    divisors = next;
  }
  return divisors.sort((a, b) => (a < b ? -1 : a > b ? 1 : 0));
}

/**
 * Clear denominators of a polynomial over QQ, returning integer coefficients
 * and the LCM of denominators.
 */
function clearDenominators<C extends RingElement>(poly: Polynomial<C>): [bigint[], bigint] {
  // Extract rational coefficients as [numerator, denominator] pairs
  const rats: Array<[bigint, bigint]> = poly.coeffs.map((c) => {
    // Handle Rational class
    if ('numerator' in c && 'denominator' in c) {
      const r = c as unknown as { numerator: bigint; denominator: bigint };
      return [r.numerator, r.denominator];
    }
    if ('numer' in c && 'denom' in c) {
      const r = c as unknown as { numer: bigint; denom: bigint };
      return [r.numer, r.denom];
    }
    // Handle _numerator and _denominator (private fields)
    if ('_numerator' in c && '_denominator' in c) {
      const r = c as unknown as { _numerator: bigint; _denominator: bigint };
      return [r._numerator, r._denominator];
    }
    // Assume integer
    const val = 'value' in c ? (c as { value: bigint }).value : BigInt(c.toString());
    return [val, 1n];
  });

  // Compute LCM of denominators
  let lcmDenom = 1n;
  for (const [_, d] of rats) {
    lcmDenom = lcm(lcmDenom, d);
  }

  // Multiply each coefficient by lcm / denom
  const intCoeffs: bigint[] = rats.map(([n, d]) => n * (lcmDenom / d));

  return [intCoeffs, lcmDenom];
}

/**
 * LCM of two bigints.
 */
function lcm(a: bigint, b: bigint): bigint {
  if (a === 0n || b === 0n) return 0n;
  const absA = a < 0n ? -a : a;
  const absB = b < 0n ? -b : b;
  return (absA / gcdBigInt(absA, absB)) * absB;
}

/**
 * Find rational roots of a polynomial over QQ.
 * Uses the rational root theorem: if p/q is a root in lowest terms,
 * then p divides the constant term and q divides the leading coefficient.
 */
function findRationalRoots<C extends RingElement>(poly: Polynomial<C>): Array<[C, number]> {
  // Clear denominators to get integer coefficients
  const [intCoeffs, _] = clearDenominators(poly);

  if (intCoeffs.length === 0) return [];
  if (intCoeffs.length === 1) return [];

  // Get divisors of constant term (for numerators)
  const constant = intCoeffs[0]!;
  const leading = intCoeffs[intCoeffs.length - 1]!;

  if (constant === 0n) {
    // 0 is a root
    let mult = 0;
    let f = intCoeffs;
    while (f.length > 0 && f[0] === 0n) {
      mult++;
      f = f.slice(1);
    }

    const roots: Array<[C, number]> = [[poly.parent.base_ring.__call__(0) as C, mult]];

    // Recursively find other roots in the deflated polynomial
    if (f.length > 1) {
      // Create deflated polynomial
      const deflatedCoeffs = f.map((c) => poly.parent.base_ring.__call__(c) as C);
      const deflated = new Polynomial(deflatedCoeffs, poly.parent);
      roots.push(...findRationalRoots(deflated));
    }

    return roots;
  }

  const numerDivisors = getDivisorsBigInt(constant < 0n ? -constant : constant);
  const denomDivisors = getDivisorsBigInt(leading < 0n ? -leading : leading);

  const roots: Array<[C, number]> = [];
  let f = intCoeffs;

  // Try all possible rational roots p/q
  for (const p of numerDivisors) {
    for (const q of denomDivisors) {
      // Try both p/q and -p/q
      for (const sign of [1n, -1n]) {
        const numer = sign * p;
        const denom = q;

        // Evaluate f at numer/denom
        // f(p/q) = sum(a_i * p^i * q^{n-i}) / q^n
        // We only need to check if the numerator is zero
        const n = f.length - 1;
        let numeratorSum = 0n;
        let pPow = 1n;
        let qPow = 1n;
        for (let i = 0; i < n; i++) qPow *= denom;

        for (let i = 0; i <= n; i++) {
          numeratorSum += f[i]! * pPow * qPow;
          pPow *= numer;
          if (i < n) qPow /= denom;
        }

        if (numeratorSum === 0n) {
          // Found a root, find its multiplicity
          let mult = 0;
          // Linear factor is (qx - p) = q(x - p/q)
          const linearFactor = [-numer, denom];

          let divResult = intPolyQuoRem(f, linearFactor);
          while (
            divResult !== null &&
            (divResult[1].length === 0 || divResult[1].every((c) => c === 0n))
          ) {
            mult++;
            f = divResult[0];
            divResult = intPolyQuoRem(f, linearFactor);
          }

          if (mult > 0) {
            // Create the rational root as a ring element
            // Try to call the ring with a rational-like object
            let rootElem: C;
            try {
              rootElem = poly.parent.base_ring.__call__({ numer, denom }) as C;
            } catch {
              // Fall back to string representation
              const g = gcdBigInt(numer < 0n ? -numer : numer, denom);
              const reducedNumer = numer / g;
              const reducedDenom = denom / g;
              if (reducedDenom === 1n) {
                rootElem = poly.parent.base_ring.__call__(reducedNumer) as C;
              } else {
                rootElem = poly.parent.base_ring.__call__(`${reducedNumer}/${reducedDenom}`) as C;
              }
            }
            roots.push([rootElem, mult]);
          }
        }
      }
    }
  }

  return roots;
}

/**
 * Base class for polynomial rings.
 */
export interface PolynomialRingBase<C extends RingElement> {
  readonly base_ring: CoefficientRing<C>;
  readonly variable_name: string;
  zero(): Polynomial<C>;
  one(): Polynomial<C>;
  gen(): Polynomial<C>;
  __call__(x: C | C[] | Polynomial<C> | number): Polynomial<C>;
}

/**
 * Test-only surface for the FLINT transcriptions used by `ZZ[x].factor()`.
 *
 * Not part of the public API: exported purely so that
 * `polynomial_factorization.test.ts` can exercise the individual routines
 * (subset enumeration, degree pruning, the CLD bound, LLL, van Hoeij) against
 * brute force and against upstream's own invariants.  Do NOT re-export this
 * from `rings/index.ts`.
 */
export const _zz_factor_internal = {
  fmpz_poly_factor_mignotte,
  fmpzPolyCLDBound,
  fmpzPolyDivlowSmodp,
  fmpzPolyDivhighSmodp,
  fmpzPolyFactorCLDMat,
  fmpzMatNextColVanHoeij,
  fmpzMatColPartition,
  lllWithRemovalKnapsack,
  vanHoeijCheckIfSolved,
  fmpzPolyFactorVanHoeij,
  liftModularFactors,
  chooseFactorizationPrime,
  zassenhausSubsetFirst,
  zassenhausSubsetNext,
  zassenhausSubsetNextDisjoint,
  zassenhausRecombination,
  zassenhausPruneSetDegree,
  zassenhausPruneStartAddFactors,
  zassenhausPruneAddFactor,
  zassenhausPruneEndAddFactors,
  zassenhausPruneDegreeIsPossible,
  factorSquarefreeIntPoly,
  factorIntegerPolynomial,
  isqrtFloor,
  fmpzSmod,
  clogUi,
};
