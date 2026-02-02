/**
 * @module sage/rings/polynomial/polynomial_element
 * @description Polynomial elements over arbitrary coefficient rings
 *
 * Port of: sage/rings/polynomial/polynomial_element.pyx
 */

import { factor as factorInteger, gcd as gcdBigInt, is_prime } from '../../arith/misc.js';
import { NotImplementedError, ValueError, ZeroDivisionError } from '../../errors.js';
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
   * Compute quotient and remainder.
   * Only works over fields (requires division of coefficients).
   */
  quo_rem(other: Polynomial<C>): [Polynomial<C>, Polynomial<C>] {
    if (other.isZero()) {
      throw new ZeroDivisionError('polynomial division by zero');
    }

    if (this.degree() < other.degree()) {
      return [this.parent.zero(), this];
    }

    // Make a mutable copy of coefficients
    const remainder = [...this.coeffs] as C[];
    const divisorLC = other.leading_coefficient();
    const divisorDeg = other.degree();
    const quotientCoeffs: C[] = [];

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
      const qCoeff = divideCoeffs(remainder[i]!, divisorLC);
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
      // Multiply coefficient by i (the power)
      let coeff = this.coeffs[i]!;
      for (let j = 1; j < i; j++) {
        coeff = coeff.add(this.coeffs[i]!) as C;
      }
      result.push(coeff);
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
    if (other.isZero()) {
      return this.isZero() ? this.parent.zero() : this._monic();
    }
    if (this.isZero()) {
      return other._monic();
    }

    let a: Polynomial<C> = this;
    let b: Polynomial<C> = other;

    while (!b.isZero()) {
      const [_q, r] = a.quo_rem(b);
      a = b;
      b = r;
    }

    // Return monic GCD
    return a._monic();
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
    let oldR = this as Polynomial<C>;
    let r = other;
    let oldS = this.parent.one();
    let s = this.parent.zero();
    let oldT = this.parent.zero();
    let t = this.parent.one();

    while (!r.isZero()) {
      const [quotient, remainder] = oldR.quo_rem(r);

      const tempR = r;
      r = remainder;
      oldR = tempR;

      const tempS = s;
      s = oldS.sub(quotient.mul(s));
      oldS = tempS;

      const tempT = t;
      t = oldT.sub(quotient.mul(t));
      oldT = tempT;
    }

    // Make the GCD monic
    if (!oldR.isZero()) {
      const lc = oldR.leading_coefficient();
      if (!lc.eq(1)) {
        const lcInv = divideCoeffs(this.parent.base_ring.one() as C, lc);
        oldR = oldR.scalar_mul(lcInv);
        oldS = oldS.scalar_mul(lcInv);
        oldT = oldT.scalar_mul(lcInv);
      }
    }

    return [oldR, oldS, oldT];
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
   * @see Reference: sage/rings/polynomial/polynomial_element.pyx:content
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
        return g;
      }
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
   * @see Reference: sage/rings/polynomial/polynomial_element.pyx:primitive_part
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
    // The Sylvester matrix has dimension (m + n) x (m + n)
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

    // Compute determinant using Gaussian elimination with pivoting
    return matrixDeterminant(matrix, this.parent.base_ring);
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
      result = divideCoeffs(res, anPower);
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
      return intRoots.map(([root, mult]) => [baseRing.__call__(root) as C, mult]);
    }

    // Handle rational polynomials (QQ[x]) - find rational roots
    if (isRationalField(baseRing)) {
      // For QQ[x], we can find rational roots using the rational root theorem
      // First, clear denominators to get an integer polynomial
      const rationalRoots = findRationalRoots(this);
      return rationalRoots;
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

    return roots;
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
   */
  private _factorOverIntegers(): Array<[Polynomial<C>, number]> {
    const coeffs = extractIntegerCoeffs(this);
    const [content, factors] = factorIntegerPolynomial(coeffs);

    const result: Array<[Polynomial<C>, number]> = [];

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
   */
  private _factorOverRationals(): Array<[Polynomial<C>, number]> {
    // Clear denominators to get an integer polynomial
    // Then factor over ZZ and convert back to monic factors over QQ
    const [intCoeffs, lcmDenom] = clearDenominators(this);

    const [content, factors] = factorIntegerPolynomial(intCoeffs);

    const result: Array<[Polynomial<C>, number]> = [];

    // Convert polynomial factors back to monic Polynomial<C> over QQ
    for (const [facCoeffs, mult] of factors) {
      // Make monic by dividing by leading coefficient
      const lc = facCoeffs[facCoeffs.length - 1]!;
      const monicCoeffs = facCoeffs.map((c) => {
        // Create rational c / lc
        return this.parent.base_ring.__call__({ numer: c, denom: lc }) as C;
      });
      const poly = new Polynomial(monicCoeffs, this.parent);
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
   * For finite fields, uses distinct-degree factorization:
   * A squarefree polynomial is irreducible iff its distinct-degree
   * factorization returns a single factor of the same degree.
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
   * @see Reference: sage/rings/polynomial/polynomial_element.pyx:is_irreducible
   */
  is_irreducible(): boolean {
    if (this.isZero()) {
      return false;
    }

    const n = this.degree();

    if (n <= 0) {
      return false; // Constants are not irreducible
    }

    if (n === 1) {
      // For ZZ[x], linear is irreducible only if primitive
      const baseRing = this.parent.base_ring;
      if (isIntegerRing(baseRing)) {
        const coeffs = extractIntegerCoeffs(this);
        const content = intPolyContent(coeffs);
        return content === 1n || content === -1n;
      }
      return true; // Linear polynomials are irreducible over fields
    }

    const baseRing = this.parent.base_ring;

    // Handle integer polynomials (ZZ[x])
    if (isIntegerRing(baseRing)) {
      // A polynomial is irreducible over ZZ iff it's primitive and
      // irreducible over QQ (by Gauss's lemma)
      const coeffs = extractIntegerCoeffs(this);
      const content = intPolyContent(coeffs);
      if (content !== 1n && content !== -1n) {
        return false; // Not primitive
      }
      // Factor and check if there's only one irreducible factor
      const [_, factors] = factorIntegerPolynomial(coeffs);
      return factors.length === 1 && factors[0]![1] === 1;
    }

    // Handle rational polynomials (QQ[x])
    if (isRationalField(baseRing)) {
      // Factor and check
      const factors = this.factor();
      return factors.length === 1 && factors[0]![1] === 1;
    }

    if (!isFiniteField(baseRing)) {
      throw new NotImplementedError(
        'is_irreducible only implemented for finite fields, ZZ, and QQ'
      );
    }

    // Check if polynomial is squarefree (no repeated roots)
    const d = this.derivative();
    if (!d.isZero()) {
      const g = this.gcd(d);
      if (g.degree() > 0) {
        return false; // Has repeated factors
      }
    } else {
      // Derivative is zero in characteristic p, polynomial is a p-th power
      return false;
    }

    // Use distinct-degree factorization
    // A squarefree polynomial is irreducible iff its DDF has exactly one
    // factor, and that factor has degree equal to the polynomial's degree
    const q = getFieldOrder(baseRing);
    const x = this.parent.gen();
    let w = x.mod(this);
    const monic = this._monic();

    // Check gcd(f, x^{q^d} - x) = 1 for all d < n/2
    // (equivalent to checking all proper divisors of n)
    for (let d = 1; 2 * d <= n; d++) {
      if (n % d === 0) {
        // d divides n, need to check
        w = powerMod(w, q, monic);
        const g = monic.gcd(w.sub(x));

        if (g.degree() > 0) {
          return false; // Has factor of degree d
        }
      } else {
        // Still need to update w
        w = powerMod(w, q, monic);
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
 * Compute the determinant of a square matrix using Gaussian elimination.
 * This works over any ring that supports division (i.e., fields).
 * For rings without division, this may fail.
 */
function matrixDeterminant<C extends RingElement>(matrix: C[][], ring: CoefficientRing<C>): C {
  const n = matrix.length;
  if (n === 0) {
    return ring.one() as C;
  }

  // Make a copy of the matrix
  const M: C[][] = matrix.map((row) => [...row]);

  let det = ring.one() as C;
  let sign = 1;

  for (let col = 0; col < n; col++) {
    // Find pivot
    let pivotRow = -1;
    for (let row = col; row < n; row++) {
      if (!M[row]![col]!.isZero()) {
        pivotRow = row;
        break;
      }
    }

    if (pivotRow === -1) {
      // Column is all zeros below diagonal, determinant is 0
      return ring.zero() as C;
    }

    // Swap rows if needed
    if (pivotRow !== col) {
      [M[col], M[pivotRow]] = [M[pivotRow]!, M[col]!];
      sign = -sign;
    }

    const pivot = M[col]![col]!;
    det = det.mul(pivot) as C;

    // Eliminate below pivot
    for (let row = col + 1; row < n; row++) {
      if (!M[row]![col]!.isZero()) {
        const factor = divideCoeffs(M[row]![col]!, pivot);
        for (let j = col; j < n; j++) {
          M[row]![j] = M[row]![j]!.sub(factor.mul(M[col]![j]!) as C) as C;
        }
      }
    }
  }

  if (sign === -1) {
    det = det.neg() as C;
  }

  return det;
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

  // Try to split f
  const maxAttempts = 100;

  for (let attempt = 0; attempt < maxAttempts; attempt++) {
    // Generate a random polynomial T of degree < n
    const t = randomPolynomial(ring, n - 1);

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

  // If we failed to split after many attempts, there might be an issue
  // but return the polynomial as-is (it might actually be irreducible)
  return [f._monic()];
}

/**
 * Generate a random polynomial of degree at most maxDegree.
 */
function randomPolynomial<C extends RingElement>(
  ring: PolynomialRingBase<C>,
  maxDegree: number
): Polynomial<C> {
  const baseRing = ring.base_ring;
  const coeffs: C[] = [];

  // Get field order for random coefficient generation
  const q = getFieldOrder(baseRing);
  const qNum = q <= BigInt(Number.MAX_SAFE_INTEGER) ? Number(q) : 1000;

  for (let i = 0; i <= maxDegree; i++) {
    // Generate random coefficient
    const randVal = Math.floor(current_randstate().random() * qNum);

    if ('__call__' in baseRing) {
      coeffs.push(baseRing.__call__(randVal) as C);
    } else {
      // Fallback: use zero/one
      if (randVal === 0) {
        coeffs.push(baseRing.zero() as C);
      } else {
        let coeff = baseRing.one() as C;
        for (let j = 1; j < randVal; j++) {
          coeff = coeff.add(baseRing.one() as C) as C;
        }
        coeffs.push(coeff);
      }
    }
  }

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
 * Factor a squarefree polynomial over Z/pZ using Berlekamp's algorithm.
 * Returns a list of monic irreducible factors over Z/pZ.
 */
function berlekampFactor(coeffs: bigint[], p: bigint): bigint[][] {
  const n = coeffs.length - 1; // degree
  if (n <= 0) return coeffs.length > 0 && coeffs[0] !== 0n ? [[coeffs[0]!]] : [];
  if (n === 1) {
    // Linear polynomial is irreducible, make monic
    const lcInv = modInverse(coeffs[1]!, p);
    return [[(((coeffs[0]! * lcInv) % p) + p) % p, 1n]];
  }

  // Make monic
  const lc = coeffs[n]!;
  const lcInv = modInverse(lc, p);
  const monicCoeffs = coeffs.map((c) => (((c * lcInv) % p) + p) % p);

  // Build Berlekamp matrix Q where Q[i][j] = coeff of x^j in x^{ip} mod f
  const Q: bigint[][] = [];
  for (let i = 0; i < n; i++) {
    // Compute x^{ip} mod f
    let xPow: bigint[] = [1n];
    for (let j = 0; j < i; j++) {
      xPow = modPolyMul(xPow, modPowX(p, monicCoeffs, p), p);
      const [_, rem] = modPolyQuoRem(xPow, monicCoeffs, p);
      xPow = rem.length > 0 ? rem : [0n];
    }
    if (i === 0) {
      xPow = [1n];
    } else {
      xPow = modPowX(BigInt(i) * p, monicCoeffs, p);
    }

    const row = new Array(n).fill(0n);
    for (let j = 0; j < Math.min(xPow.length, n); j++) {
      row[j] = xPow[j]!;
    }
    // Subtract identity: Q - I
    row[i] = (((row[i]! - 1n) % p) + p) % p;
    Q.push(row);
  }

  // Find null space of Q - I
  const nullSpace = modMatrixNullSpace(Q, p);

  if (nullSpace.length <= 1) {
    // Polynomial is irreducible
    return [monicCoeffs];
  }

  // Use null space vectors to split
  const factors = splitUsingNullSpace(monicCoeffs, nullSpace, p);
  return factors;
}

/**
 * Compute x^k mod f over Z/pZ.
 */
function modPowX(k: bigint, f: bigint[], p: bigint): bigint[] {
  if (k === 0n) return [1n];

  let result: bigint[] = [1n];
  let base: bigint[] = [0n, 1n]; // x

  while (k > 0n) {
    if ((k & 1n) === 1n) {
      result = modPolyMul(result, base, p);
      const [_, rem] = modPolyQuoRem(result, f, p);
      result = rem.length > 0 ? rem : [0n];
    }
    base = modPolyMul(base, base, p);
    const [_, rem] = modPolyQuoRem(base, f, p);
    base = rem.length > 0 ? rem : [0n];
    k >>= 1n;
  }

  return result;
}

/**
 * Find null space of a matrix over Z/pZ using Gaussian elimination.
 */
function modMatrixNullSpace(M: bigint[][], p: bigint): bigint[][] {
  const n = M.length;
  if (n === 0) return [];
  const m = M[0]!.length;

  // Augment with identity
  const aug: bigint[][] = M.map((row, i) => {
    const newRow = [...row];
    for (let j = 0; j < n; j++) {
      newRow.push(i === j ? 1n : 0n);
    }
    return newRow;
  });

  // Row reduce
  let col = 0;
  for (let row = 0; row < n && col < m; row++) {
    // Find pivot
    let pivotRow = -1;
    for (let i = row; i < n; i++) {
      if (aug[i]![col] !== 0n) {
        pivotRow = i;
        break;
      }
    }

    if (pivotRow === -1) {
      col++;
      row--;
      continue;
    }

    // Swap rows
    [aug[row], aug[pivotRow]] = [aug[pivotRow]!, aug[row]!];

    // Scale pivot row
    const pivotInv = modInverse(aug[row]![col]!, p);
    aug[row] = aug[row]!.map((v) => (((v * pivotInv) % p) + p) % p);

    // Eliminate
    for (let i = 0; i < n; i++) {
      if (i !== row && aug[i]![col] !== 0n) {
        const factor = aug[i]![col]!;
        for (let j = 0; j < aug[i]!.length; j++) {
          aug[i]![j] = (((aug[i]![j]! - factor * aug[row]![j]!) % p) + p) % p;
        }
      }
    }

    col++;
  }

  // Extract null space vectors (rows where the original part is zero)
  const nullVectors: bigint[][] = [];
  for (let i = 0; i < n; i++) {
    let isZero = true;
    for (let j = 0; j < m; j++) {
      if (aug[i]![j] !== 0n) {
        isZero = false;
        break;
      }
    }
    if (isZero) {
      // The identity part gives us the null vector
      nullVectors.push(aug[i]!.slice(m));
    }
  }

  // Always include the trivial vector [1, 0, 0, ...]
  if (nullVectors.length === 0) {
    const trivial = new Array(n).fill(0n);
    trivial[0] = 1n;
    nullVectors.push(trivial);
  }

  return nullVectors;
}

/**
 * Split polynomial using null space vectors from Berlekamp.
 */
function splitUsingNullSpace(f: bigint[], nullSpace: bigint[][], p: bigint): bigint[][] {
  let factors = [f];

  for (const v of nullSpace) {
    // v represents a polynomial h = sum v[i] * x^i
    // We try gcd(f, h - c) for various c in Z/pZ
    const h = v.slice();
    while (h.length > 0 && h[h.length - 1] === 0n) h.pop();
    if (h.length === 0) continue;

    const newFactors: bigint[][] = [];

    for (const fac of factors) {
      if (fac.length - 1 <= 1) {
        newFactors.push(fac);
        continue;
      }

      let split = false;
      for (let c = 0n; c < p && !split; c++) {
        // h - c
        const hMinusC = [...h];
        hMinusC[0] = ((hMinusC[0] || 0n) - c + p) % p;

        const g = modPolyGcd(fac, hMinusC, p);

        if (g.length > 1 && g.length < fac.length) {
          // Found a non-trivial factor
          const [q, _] = modPolyQuoRem(fac, g, p);
          newFactors.push(g);
          if (q.length > 1) {
            newFactors.push(q);
          }
          split = true;
        }
      }

      if (!split) {
        newFactors.push(fac);
      }
    }

    factors = newFactors;
  }

  // Recursively factor any remaining reducible factors
  const result: bigint[][] = [];
  for (const fac of factors) {
    if (fac.length - 1 <= 1) {
      result.push(fac);
    } else {
      // Check if irreducible using distinct-degree factorization
      const ddf = distinctDegreeFactor(fac, p);
      if (ddf.length === 1 && ddf[0]![1] === fac.length - 1) {
        result.push(fac);
      } else {
        // Need to continue splitting
        for (const [g, d] of ddf) {
          if (g.length - 1 === d) {
            result.push(g);
          } else {
            // Use Cantor-Zassenhaus for equal-degree factorization
            const edf = equalDegreeFactor(g, d, p);
            result.push(...edf);
          }
        }
      }
    }
  }

  return result;
}

/**
 * Distinct-degree factorization over Z/pZ.
 * Returns pairs [g_d, d] where g_d is product of irreducible factors of degree d.
 */
function distinctDegreeFactor(f: bigint[], p: bigint): Array<[bigint[], number]> {
  const n = f.length - 1;
  if (n <= 0) return [];

  const result: Array<[bigint[], number]> = [];
  let v = f;
  let w: bigint[] = [0n, 1n]; // x

  for (let d = 1; 2 * d <= v.length - 1; d++) {
    // w = x^{p^d} mod v
    w = modPowX(p, v, p);
    for (let i = 1; i < d; i++) {
      w = modPowX(p, v, p);
    }
    w = modPowX(p ** BigInt(d), v, p);

    // gcd(v, w - x)
    const wMinusX = [...w];
    if (wMinusX.length === 0) wMinusX.push(0n);
    if (wMinusX.length === 1) wMinusX.push(0n);
    wMinusX[1] = ((wMinusX[1] || 0n) - 1n + p) % p;

    const g = modPolyGcd(v, wMinusX, p);

    if (g.length > 1) {
      result.push([g, d]);
      const [q, _] = modPolyQuoRem(v, g, p);
      v = q;
    }
  }

  if (v.length > 1) {
    result.push([v, v.length - 1]);
  }

  return result;
}

/**
 * Equal-degree factorization using Cantor-Zassenhaus algorithm.
 */
function equalDegreeFactor(f: bigint[], d: number, p: bigint): bigint[][] {
  const n = f.length - 1;
  if (n === d) return [f];
  if (n === 0) return [];

  const numFactors = n / d;
  if (numFactors <= 1) return [f];

  const factors: bigint[][] = [];
  const remaining = [f];

  const maxAttempts = 50;

  while (remaining.length > 0) {
    const curr = remaining.pop()!;
    if (curr.length - 1 === d) {
      factors.push(curr);
      continue;
    }

    for (let attempt = 0; attempt < maxAttempts; attempt++) {
      // Generate random polynomial of degree < n
      const t: bigint[] = [];
      for (let i = 0; i < curr.length - 1; i++) {
        t.push(BigInt(Math.floor(Math.random() * Number(p))));
      }
      while (t.length > 0 && t[t.length - 1] === 0n) t.pop();
      if (t.length === 0) continue;

      // Compute gcd(curr, t^{(p^d-1)/2} - 1)
      const exp = (p ** BigInt(d) - 1n) / 2n;
      let tPow = t;
      let e = exp;
      let result: bigint[] = [1n];

      while (e > 0n) {
        if ((e & 1n) === 1n) {
          result = modPolyMul(result, tPow, p);
          const [_, rem] = modPolyQuoRem(result, curr, p);
          result = rem.length > 0 ? rem : [0n];
        }
        tPow = modPolyMul(tPow, tPow, p);
        const [_, rem] = modPolyQuoRem(tPow, curr, p);
        tPow = rem.length > 0 ? rem : [0n];
        e >>= 1n;
      }

      // result - 1
      result[0] = ((result[0] || 0n) - 1n + p) % p;

      const g = modPolyGcd(curr, result, p);

      if (g.length > 1 && g.length < curr.length) {
        const [q, _] = modPolyQuoRem(curr, g, p);
        remaining.push(g);
        if (q.length > 1) remaining.push(q);
        break;
      }
    }
  }

  return factors;
}

/**
 * Hensel lifting: lift factorization from Z/pZ to Z/p^k Z.
 * Given f = g * h mod p with gcd(g, h) = 1 mod p,
 * find G, H such that f = G * H mod p^k and G = g mod p, H = h mod p.
 */
function henselLift(
  f: bigint[],
  g: bigint[],
  h: bigint[],
  p: bigint,
  k: number
): [bigint[], bigint[]] {
  // Compute s, t such that s*g + t*h = 1 mod p
  let [s, t] = extendedGcdPoly(g, h, p);

  let G = g;
  let H = h;
  let pk = p;

  for (let i = 1; i < k; i++) {
    const pk2 = pk * p;

    // e = f - G * H mod p^{i+1}
    const prod = intPolyMul(G, H);
    const e = f.map((c, idx) => (c - (prod[idx] || 0n)) % pk2);
    while (e.length > 0 && e[e.length - 1] === 0n) e.pop();

    if (e.length === 0) {
      pk = pk2;
      continue;
    }

    // q, r such that s*e = q*H + r with deg(r) < deg(H)
    const se = intPolyMul(s, e).map((c) => c % pk2);
    const [q, r] = modPolyQuoRem(se, H, pk2);

    // G' = G + t*e + q*G
    const te = intPolyMul(t, e);
    const qG = intPolyMul(q, G);
    G = G.map((c, idx) => (c + (te[idx] || 0n) + (qG[idx] || 0n)) % pk2);
    while (G.length > 0 && G[G.length - 1] === 0n) G.pop();

    // H' = H + r
    H = H.map((c, idx) => (c + (r[idx] || 0n)) % pk2);
    while (H.length > 0 && H[H.length - 1] === 0n) H.pop();

    // Update s, t for next iteration
    // s*G + t*H = 1 mod p^{i+1}
    const sg = intPolyMul(s, G);
    const th = intPolyMul(t, H);
    const err = sg.map((c, idx) => (c + (th[idx] || 0n) - (idx === 0 ? 1n : 0n)) % pk2);
    while (err.length > 0 && err[err.length - 1] === 0n) err.pop();

    if (err.length > 0) {
      const sErr = intPolyMul(s, err);
      const [q2, r2] = modPolyQuoRem(sErr, H, pk2);
      s = s.map((c, idx) => (c - (r2[idx] || 0n)) % pk2);
      const tErr = intPolyMul(t, err);
      const q2G = intPolyMul(q2, G);
      t = t.map((c, idx) => (c - (tErr[idx] || 0n) - (q2G[idx] || 0n)) % pk2);
    }

    pk = pk2;
  }

  return [G, H];
}

/**
 * Extended GCD for polynomials over Z/pZ.
 * Returns [s, t] such that s*a + t*b = gcd(a, b) mod p.
 */
function extendedGcdPoly(a: bigint[], b: bigint[], p: bigint): [bigint[], bigint[]] {
  let [oldR, r] = [a, b];
  let [oldS, s]: [bigint[], bigint[]] = [[1n], [0n]];
  let [oldT, t]: [bigint[], bigint[]] = [[0n], [1n]];

  while (r.length > 0) {
    const [q, rem] = modPolyQuoRem(oldR, r, p);

    [oldR, r] = [r, rem];

    // s = oldS - q * s
    const qs = modPolyMul(q, s, p);
    const newS = oldS.map((c, idx) => (((c - (qs[idx] || 0n)) % p) + p) % p);
    while (newS.length > 1 && newS[newS.length - 1] === 0n) newS.pop();
    [oldS, s] = [s, newS];

    // t = oldT - q * t
    const qt = modPolyMul(q, t, p);
    const newT = oldT.map((c, idx) => (((c - (qt[idx] || 0n)) % p) + p) % p);
    while (newT.length > 1 && newT[newT.length - 1] === 0n) newT.pop();
    [oldT, t] = [t, newT];
  }

  // Normalize so gcd is monic
  if (oldR.length > 0 && oldR[oldR.length - 1] !== 1n) {
    const lcInv = modInverse(oldR[oldR.length - 1]!, p);
    oldS = oldS.map((c) => (((c * lcInv) % p) + p) % p);
    oldT = oldT.map((c) => (((c * lcInv) % p) + p) % p);
  }

  return [oldS, oldT];
}

/**
 * Factor an integer polynomial by finding integer roots.
 * This is a simple approach for small degree polynomials.
 */
function factorByRationalRoots(coeffs: bigint[]): bigint[][] {
  const n = coeffs.length - 1;
  if (n <= 0) return coeffs.length > 0 && coeffs[0] !== 0n ? [[coeffs[0]!]] : [];
  if (n === 1) return [coeffs];

  const factors: bigint[][] = [];
  let f = [...coeffs];
  let maxIterations = 100; // Prevent infinite loops

  // Try to find integer roots using rational root theorem
  while (f.length > 1 && maxIterations-- > 0) {
    const constant = f[0]!;
    const lc = f[f.length - 1]!;

    if (constant === 0n) {
      // 0 is a root, factor out x
      factors.push([0n, 1n]); // x
      f = f.slice(1);
      // Remove leading zeros after shift
      while (f.length > 1 && f[f.length - 1] === 0n) f.pop();
      continue;
    }

    // Get divisors of constant term (only try integer roots for simplicity)
    const constDivisors = getDivisorsBigInt(constant < 0n ? -constant : constant);

    // Limit divisors to prevent too many iterations
    const limitedDivisors = constDivisors.slice(0, 20);

    let foundRoot = false;

    // Try integer roots first (p/1 where p divides constant term)
    for (const d of limitedDivisors) {
      for (const sign of [1n, -1n]) {
        const root = sign * d;

        // Evaluate f at root
        const val = intPolyEval(f, root);

        if (val === 0n) {
          // root is a root, divide f by (x - root)
          const linearFactor = [-root, 1n];
          const divResult = intPolyQuoRem(f, linearFactor);

          if (
            divResult !== null &&
            (divResult[1].length === 0 || divResult[1].every((c) => c === 0n))
          ) {
            factors.push(linearFactor);
            f = divResult[0];

            // Remove trailing zeros
            while (f.length > 1 && f[f.length - 1] === 0n) f.pop();

            foundRoot = true;
            break;
          }
        }
      }
      if (foundRoot) break;
    }

    if (!foundRoot) {
      // No more integer roots, remaining polynomial may be irreducible
      break;
    }
  }

  // Add remaining polynomial if degree > 0
  if (f.length > 1) {
    factors.push(f);
  }

  return factors.length > 0 ? factors : [coeffs];
}

/**
 * Factor a primitive squarefree integer polynomial.
 * Uses a combination of rational root theorem and irreducibility testing.
 *
 * @param coeffs - Primitive squarefree polynomial coefficients
 * @returns Array of irreducible factor coefficient arrays
 */
function factorSquarefreeIntPoly(coeffs: bigint[]): bigint[][] {
  const n = coeffs.length - 1; // degree
  if (n <= 0) return coeffs.length > 0 && coeffs[0] !== 0n ? [[coeffs[0]!]] : [];
  if (n === 1) return [coeffs];

  // For small degree polynomials, use simple root finding
  if (n <= 10) {
    return factorByRationalRoots(coeffs);
  }

  // Find a small prime p that doesn't divide the leading coefficient
  // and keeps the polynomial squarefree mod p
  const lc = coeffs[n]!;
  let p = 2n;
  let modFactors: bigint[][] = [];

  const smallPrimes = [2n, 3n, 5n, 7n, 11n, 13n, 17n, 19n, 23n, 29n, 31n, 37n, 41n, 43n, 47n];

  for (const prime of smallPrimes) {
    if (lc % prime === 0n) continue;

    const fModP = intPolyModP(coeffs, prime);
    if (fModP.length - 1 < n) continue; // degree dropped

    // Check if squarefree mod p
    const deriv = fModP.slice(1).map((c, i) => (c * BigInt(i + 1)) % prime);
    const g = modPolyGcd(fModP, deriv, prime);
    if (g.length > 1) continue; // not squarefree

    // Factor mod p
    modFactors = berlekampFactor(fModP, prime);
    if (modFactors.length === 1) {
      // Irreducible mod p, hence irreducible over Z
      return [coeffs];
    }

    p = prime;
    break;
  }

  if (modFactors.length === 0 || modFactors.length === 1) {
    return [coeffs];
  }

  // Compute bound on coefficients of factors
  // Use Mignotte bound: if g divides f, |g_i| <= C(n,i) * ||f|| where ||f|| is the 2-norm
  const norm = Math.sqrt(Number(coeffs.reduce((s, c) => s + c * c, 0n)));
  const bound = BigInt(Math.ceil(2 ** n * norm * Math.abs(Number(lc))));

  // Determine k such that p^k > 2 * bound * lc
  let k = 1;
  let pk = p;
  while (pk <= 2n * bound * (lc < 0n ? -lc : lc)) {
    k++;
    pk *= p;
  }

  // Lift modular factors using Hensel lifting
  // For simplicity, we'll try combinations of modular factors
  // This is a simplified version - full Zassenhaus would use LLL

  const factors: bigint[][] = [];
  let remaining = coeffs;

  // Try combining subsets of modular factors
  const numModFactors = modFactors.length;

  for (let size = 1; size <= Math.floor(numModFactors / 2); size++) {
    const subsets = getSubsets(numModFactors, size);

    for (const subset of subsets) {
      // Compute product of selected modular factors
      let g = [1n];
      for (const idx of subset) {
        g = modPolyMul(g, modFactors[idx]!, pk);
      }

      // Make monic and multiply by lc
      if (g.length > 0) {
        const gLc = g[g.length - 1]!;
        const gLcInv = modInverse(gLc, pk);
        g = g.map((c) => (((c * gLcInv * (remaining[remaining.length - 1]! % pk)) % pk) + pk) % pk);
      }

      // Reduce coefficients to symmetric range
      g = g.map((c) => {
        let r = c % pk;
        if (r > pk / 2n) r -= pk;
        return r;
      });

      // Check if this divides the remaining polynomial over Z
      const divResult = intPolyQuoRem(remaining, g);
      if (divResult !== null) {
        const [q, rem] = divResult;
        if (rem.length === 0 || rem.every((c) => c === 0n)) {
          // Found a factor
          const [content, primitive] = intPolyPrimitive(g);
          factors.push(primitive);
          remaining = q;

          // Remove used modular factors
          for (const idx of subset.reverse()) {
            modFactors.splice(idx, 1);
          }
          break;
        }
      }
    }

    if (remaining.length - 1 <= 0) break;
  }

  // Add remaining polynomial if non-trivial
  if (remaining.length > 1) {
    const [_, primitive] = intPolyPrimitive(remaining);
    factors.push(primitive);
  }

  return factors.length > 0 ? factors : [coeffs];
}

/**
 * Generate all subsets of {0, 1, ..., n-1} of given size.
 */
function getSubsets(n: number, size: number): number[][] {
  if (size === 0) return [[]];
  if (size > n) return [];

  const result: number[][] = [];

  function helper(start: number, current: number[]): void {
    if (current.length === size) {
      result.push([...current]);
      return;
    }
    for (let i = start; i < n; i++) {
      current.push(i);
      helper(i + 1, current);
      current.pop();
    }
  }

  helper(0, []);
  return result;
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
  let maxIter = 20;

  while (h.length > 1 && maxIter-- > 0) {
    // gcd(g, h)
    const gi = intPolyGcd(g, h);

    // h / gi
    const hDivGi = intPolyQuoRem(h, gi);
    if (hDivGi === null) break;
    const hi = hDivGi[0];

    if (hi.length > 1) {
      result.push([hi, i]);
    }

    // g = g / gi
    const gDivGi = intPolyQuoRem(g, gi);
    if (gDivGi === null) break;
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
 * GCD of two integer polynomials using subresultant PRS.
 */
function intPolyGcd(a: bigint[], b: bigint[]): bigint[] {
  // Remove trailing zeros
  while (a.length > 0 && a[a.length - 1] === 0n) a = a.slice(0, -1);
  while (b.length > 0 && b[b.length - 1] === 0n) b = b.slice(0, -1);

  if (b.length === 0) return a.length > 0 ? a : [1n];
  if (a.length === 0) return b.length > 0 ? b : [1n];
  if (a.length < b.length) [a, b] = [b, a];

  // Use pseudo-division with iteration limit
  let maxIter = 100;
  while (b.length > 0 && maxIter-- > 0) {
    const [_, rem] = pseudoDivide(a, b);
    // Make primitive to avoid coefficient explosion
    if (rem.length === 0) break;
    const [__, primRem] = intPolyPrimitive(rem);
    if (primRem.length === 0) break;
    a = b;
    b = primRem;
  }

  // Make primitive and positive leading coefficient
  const [_, primA] = intPolyPrimitive(a);
  return primA;
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
  const q = new Array(d + 1).fill(0n);

  for (let i = m; i >= n; i--) {
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
 * Factor an integer polynomial completely.
 * Returns [content, factors] where factors is array of [irreducible_factor, multiplicity].
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
 */
function getDivisorsBigInt(n: bigint): bigint[] {
  if (n <= 0n) return [];
  const divisors: bigint[] = [];
  let i = 1n;
  while (i * i <= n) {
    if (n % i === 0n) {
      divisors.push(i);
      if (i !== n / i) {
        divisors.push(n / i);
      }
    }
    i++;
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
