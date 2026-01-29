/**
 * @module sage/rings/rational
 * @description Rational numbers
 *
 * Port of: sage/rings/rational.pyx
 * Reference: reference/sage/src/sage/rings/rational.pyx
 */

import {
  gcd,
  lcm,
  valuation as integer_valuation,
  is_square as integer_is_square,
  isqrt,
  prime_factors,
  factorial as integer_factorial,
} from '../arith/misc.js';
import { type IntegerLike, toBigInt } from '../types/coercion.js';
import {
  ArithmeticError,
  NotImplementedError,
  TypeError,
  ValueError,
  ZeroDivisionError,
} from '../errors.js';
import { Polynomial, type PolynomialRingBase } from './polynomial/polynomial_element.js';
import { PolynomialRing } from './polynomial/polynomial_ring.js';

/**
 * Compute the integer n-th root of a non-negative integer.
 *
 * @param x - Non-negative integer
 * @param n - Positive integer root
 * @returns [exact, root] where exact is true if root^n === x
 */
function integerNthRoot(x: bigint, n: bigint): [boolean, bigint] {
  if (x === 0n) {
    return [true, 0n];
  }
  if (x === 1n) {
    return [true, 1n];
  }
  if (n === 1n) {
    return [true, x];
  }
  if (n === 2n) {
    const root = isqrt(x);
    return [root * root === x, root];
  }

  // Newton's method for n-th root
  // Initial estimate using bit length
  const bitLen = x.toString(2).length;
  let root = 1n << (BigInt(Math.ceil(bitLen / Number(n))) + 1n);

  // Newton iteration: root = ((n-1)*root + x/root^(n-1)) / n
  const nMinus1 = n - 1n;
  let prevRoot = 0n;

  while (root !== prevRoot) {
    prevRoot = root;
    const rootPowNMinus1 = bigintPow(root, nMinus1);
    const newRoot = (nMinus1 * root + x / rootPowNMinus1) / n;
    root = newRoot;
  }

  // The actual root might be off by 1 due to integer division
  // Check root and root+1
  let rootPowN = bigintPow(root, n);
  if (rootPowN === x) {
    return [true, root];
  }
  if (rootPowN > x) {
    root--;
    rootPowN = bigintPow(root, n);
    if (rootPowN === x) {
      return [true, root];
    }
  } else {
    const nextRoot = root + 1n;
    const nextPow = bigintPow(nextRoot, n);
    if (nextPow === x) {
      return [true, nextRoot];
    }
    if (nextPow < x) {
      root = nextRoot;
    }
  }

  return [false, root];
}

/**
 * Compute base^exp for bigints.
 */
function bigintPow(base: bigint, exp: bigint): bigint {
  if (exp === 0n) return 1n;
  if (exp === 1n) return base;

  let result = 1n;
  let b = base;
  let e = exp;

  while (e > 0n) {
    if ((e & 1n) === 1n) {
      result *= b;
    }
    b *= b;
    e >>= 1n;
  }

  return result;
}

/**
 * Return the ordinal string for a number (1st, 2nd, 3rd, etc.)
 */
function ordinalStr(n: bigint): string {
  const num = Number(n);
  const lastDigit = num % 10;
  const lastTwoDigits = num % 100;

  if (lastTwoDigits >= 11 && lastTwoDigits <= 13) {
    return `${num}th`;
  }

  switch (lastDigit) {
    case 1:
      return `${num}st`;
    case 2:
      return `${num}nd`;
    case 3:
      return `${num}rd`;
    default:
      return `${num}th`;
  }
}

/**
 * A rational number.
 *
 * Rational numbers are stored as a numerator and denominator,
 * always in lowest terms with the denominator positive.
 *
 * @example
 * ```typescript
 * const r = new Rational(3n, 4n);    // 3/4
 * const s = new Rational(-2n, 6n);   // reduced to -1/3
 * const t = Rational.fromString("5/7");
 * ```
 */
export class Rational {
  private readonly _numerator: bigint;
  private readonly _denominator: bigint;

  /**
   * Create a new rational number.
   *
   * The result is always stored in lowest terms with a positive denominator.
   *
   * @param numerator - The numerator
   * @param denominator - The denominator (must be non-zero)
   * @throws {ZeroDivisionError} If denominator is zero
   */
  constructor(numerator: IntegerLike, denominator: IntegerLike = 1n) {
    let num = toBigInt(numerator);
    let den = toBigInt(denominator);

    if (den === 0n) {
      throw new ZeroDivisionError('denominator must not be 0');
    }

    // Normalize: always keep denominator positive
    if (den < 0n) {
      num = -num;
      den = -den;
    }

    // Reduce to lowest terms
    if (num === 0n) {
      this._numerator = 0n;
      this._denominator = 1n;
    } else {
      const g = gcd(num < 0n ? -num : num, den);
      this._numerator = num / g;
      this._denominator = den / g;
    }
  }

  /**
   * Return the numerator of this rational number.
   */
  get numerator(): bigint {
    return this._numerator;
  }

  /**
   * Alias for numerator.
   */
  get numer(): bigint {
    return this._numerator;
  }

  /**
   * Return the denominator of this rational number.
   */
  get denominator(): bigint {
    return this._denominator;
  }

  /**
   * Alias for denominator.
   */
  get denom(): bigint {
    return this._denominator;
  }

  /**
   * Return the sign of this rational number (-1, 0, or 1).
   */
  get sign(): bigint {
    if (this._numerator < 0n) return -1n;
    if (this._numerator > 0n) return 1n;
    return 0n;
  }

  /**
   * Create a Rational from various input types.
   *
   * @param value - Can be a number, bigint, string ("n/d" format), or Rational
   * @returns A new Rational
   */
  static from(value: number | bigint | string | Rational): Rational {
    if (value instanceof Rational) {
      return value;
    }

    if (typeof value === 'bigint') {
      return new Rational(value, 1n);
    }

    if (typeof value === 'number') {
      if (!Number.isFinite(value)) {
        throw new TypeError('cannot convert infinity or NaN to rational');
      }
      if (Number.isInteger(value)) {
        return new Rational(BigInt(value), 1n);
      }
      // Convert float to rational using continued fractions approach
      // For simplicity, we use a string-based approach
      const str = value.toString();
      return Rational.fromString(str);
    }

    if (typeof value === 'string') {
      return Rational.fromString(value);
    }

    throw new TypeError(`cannot convert ${typeof value} to Rational`);
  }

  /**
   * Create a Rational from a string.
   *
   * Accepts formats: "n", "n/d", or decimal like "1.5"
   *
   * @param str - The string representation
   * @returns A new Rational
   */
  static fromString(str: string): Rational {
    str = str.trim();

    if (str.includes('/')) {
      const parts = str.split('/');
      if (parts.length !== 2) {
        throw new TypeError(`unable to convert '${str}' to a rational`);
      }
      const num = BigInt(parts[0]!.trim());
      const den = BigInt(parts[1]!.trim());
      return new Rational(num, den);
    }

    if (str.includes('.')) {
      // Handle decimal notation
      const parts = str.split('.');
      if (parts.length !== 2) {
        throw new TypeError(`unable to convert '${str}' to a rational`);
      }
      const intPart = parts[0]!;
      const fracPart = parts[1]!;
      const sign = str.startsWith('-') ? -1n : 1n;
      const absInt = intPart.startsWith('-') ? intPart.slice(1) : intPart;

      const denominator = 10n ** BigInt(fracPart.length);
      const numerator = sign * (BigInt(absInt || '0') * denominator + BigInt(fracPart));

      return new Rational(numerator, denominator);
    }

    // Plain integer
    return new Rational(BigInt(str), 1n);
  }

  /**
   * Create a Rational from a tuple [numerator, denominator].
   */
  static fromTuple(tuple: [bigint, bigint]): Rational {
    return new Rational(tuple[0], tuple[1]);
  }

  /**
   * Return zero.
   */
  static zero(): Rational {
    return new Rational(0n, 1n);
  }

  /**
   * Return one.
   */
  static one(): Rational {
    return new Rational(1n, 1n);
  }

  // ============================================
  // Arithmetic operations
  // ============================================

  /**
   * Add two rational numbers.
   */
  add(other: Rational | IntegerLike): Rational {
    if (other instanceof Rational) {
      // a/b + c/d = (a*d + c*b) / (b*d)
      const num = this._numerator * other._denominator + other._numerator * this._denominator;
      const den = this._denominator * other._denominator;
      return new Rational(num, den);
    }
    const otherBig = toBigInt(other);
    return this.add(new Rational(otherBig, 1n));
  }

  /**
   * Subtract a rational number from this one.
   */
  sub(other: Rational | IntegerLike): Rational {
    if (other instanceof Rational) {
      // a/b - c/d = (a*d - c*b) / (b*d)
      const num = this._numerator * other._denominator - other._numerator * this._denominator;
      const den = this._denominator * other._denominator;
      return new Rational(num, den);
    }
    const otherBig = toBigInt(other);
    return this.sub(new Rational(otherBig, 1n));
  }

  /**
   * Multiply two rational numbers.
   */
  mul(other: Rational | IntegerLike): Rational {
    if (other instanceof Rational) {
      return new Rational(this._numerator * other._numerator, this._denominator * other._denominator);
    }
    const otherBig = toBigInt(other);
    return this.mul(new Rational(otherBig, 1n));
  }

  /**
   * Divide this rational by another.
   */
  div(other: Rational | IntegerLike): Rational {
    if (other instanceof Rational) {
      if (other._numerator === 0n) {
        throw new ZeroDivisionError('rational division by zero');
      }
      return new Rational(this._numerator * other._denominator, this._denominator * other._numerator);
    }
    const otherBig = toBigInt(other);
    return this.div(new Rational(otherBig, 1n));
  }

  /**
   * Return the negation of this rational.
   */
  neg(): Rational {
    return new Rational(-this._numerator, this._denominator);
  }

  /**
   * Return the multiplicative inverse (reciprocal) of this rational.
   */
  inv(): Rational {
    if (this._numerator === 0n) {
      throw new ZeroDivisionError('rational division by zero');
    }
    return new Rational(this._denominator, this._numerator);
  }

  /**
   * Return this rational raised to an integer power.
   *
   * @param n - The exponent (can be negative)
   */
  pow(n: IntegerLike): Rational {
    let nBig = toBigInt(n);
    if (nBig === 0n) {
      return Rational.one();
    }

    if (nBig < 0n) {
      if (this._numerator === 0n) {
        throw new ZeroDivisionError('rational division by zero');
      }
      // (a/b)^(-n) = (b/a)^n
      return this.inv().pow(-nBig);
    }

    // Use binary exponentiation
    let result = Rational.one();
    let base: Rational = this;

    while (nBig > 0n) {
      if ((nBig & 1n) === 1n) {
        result = result.mul(base);
      }
      base = base.mul(base);
      nBig >>= 1n;
    }

    return result;
  }

  /**
   * Return the absolute value of this rational.
   */
  abs(): Rational {
    if (this._numerator < 0n) {
      return new Rational(-this._numerator, this._denominator);
    }
    return this;
  }

  // ============================================
  // Comparison operations
  // ============================================

  /**
   * Test equality with another rational or integer.
   */
  eq(other: Rational | IntegerLike): boolean {
    if (other instanceof Rational) {
      return this._numerator === other._numerator && this._denominator === other._denominator;
    }
    const otherBig = toBigInt(other);
    return this._numerator === otherBig && this._denominator === 1n;
  }

  /**
   * Test if this rational is less than another.
   */
  lt(other: Rational | IntegerLike): boolean {
    if (other instanceof Rational) {
      // a/b < c/d iff a*d < c*b (when denominators are positive)
      return this._numerator * other._denominator < other._numerator * this._denominator;
    }
    const otherBig = toBigInt(other);
    return this._numerator < otherBig * this._denominator;
  }

  /**
   * Test if this rational is less than or equal to another.
   */
  le(other: Rational | IntegerLike): boolean {
    if (other instanceof Rational) {
      return this._numerator * other._denominator <= other._numerator * this._denominator;
    }
    const otherBig = toBigInt(other);
    return this._numerator <= otherBig * this._denominator;
  }

  /**
   * Test if this rational is greater than another.
   */
  gt(other: Rational | IntegerLike): boolean {
    if (other instanceof Rational) {
      return this._numerator * other._denominator > other._numerator * this._denominator;
    }
    const otherBig = toBigInt(other);
    return this._numerator > otherBig * this._denominator;
  }

  /**
   * Test if this rational is greater than or equal to another.
   */
  ge(other: Rational | IntegerLike): boolean {
    if (other instanceof Rational) {
      return this._numerator * other._denominator >= other._numerator * this._denominator;
    }
    const otherBig = toBigInt(other);
    return this._numerator >= otherBig * this._denominator;
  }

  /**
   * Compare this rational with another.
   *
   * @returns -1 if this < other, 0 if equal, 1 if this > other
   */
  cmp(other: Rational | IntegerLike): -1 | 0 | 1 {
    let rhs: bigint;
    let lhs: bigint;
    if (other instanceof Rational) {
      lhs = this._numerator * other._denominator;
      rhs = other._numerator * this._denominator;
    } else {
      const otherBig = toBigInt(other);
      lhs = this._numerator;
      rhs = otherBig * this._denominator;
    }
    if (lhs < rhs) return -1;
    if (lhs > rhs) return 1;
    return 0;
  }

  // ============================================
  // Conversion methods
  // ============================================

  /**
   * Return a string representation of this rational.
   */
  toString(): string {
    if (this._denominator === 1n) {
      return this._numerator.toString();
    }
    return `${this._numerator}/${this._denominator}`;
  }

  /**
   * Return a floating point approximation.
   */
  toNumber(): number {
    return Number(this._numerator) / Number(this._denominator);
  }

  /**
   * Return the floor of this rational (greatest integer <= this).
   */
  floor(): bigint {
    if (this._denominator === 1n) {
      return this._numerator;
    }
    // For positive numbers: floor(a/b) = a div b
    // For negative numbers: floor(a/b) = (a - (b-1)) div b = a div b - 1 if remainder != 0
    const q = this._numerator / this._denominator;
    if (this._numerator >= 0n) {
      return q;
    }
    // Negative case: if there's a remainder, subtract 1
    const rem = this._numerator % this._denominator;
    if (rem !== 0n) {
      return q - 1n;
    }
    return q;
  }

  /**
   * Return the ceiling of this rational (least integer >= this).
   */
  ceil(): bigint {
    if (this._denominator === 1n) {
      return this._numerator;
    }
    // For negative numbers: ceil(a/b) = a div b
    // For positive numbers: ceil(a/b) = a div b + 1 if remainder != 0
    const q = this._numerator / this._denominator;
    if (this._numerator <= 0n) {
      return q;
    }
    // Positive case: if there's a remainder, add 1
    const rem = this._numerator % this._denominator;
    if (rem !== 0n) {
      return q + 1n;
    }
    return q;
  }

  /**
   * Return the nearest integer to this rational.
   *
   * @param mode - Rounding mode for half integers:
   *   - 'even': round toward even (default, banker's rounding)
   *   - 'away': round away from zero
   *   - 'toward': round toward zero
   *   - 'up': round up (toward +infinity)
   *   - 'down': round down (toward -infinity)
   *   - 'odd': round toward odd
   */
  round(mode: 'even' | 'away' | 'toward' | 'up' | 'down' | 'odd' = 'even'): bigint {
    if (this._denominator === 1n) {
      return this._numerator;
    }

    // Check if exactly at a half (denominator is 2)
    if (this._denominator === 2n) {
      // For n/2 where n is odd:
      // - floor is (n-1)/2 for positive, (n-1)/2 for negative (but n is negative so it's smaller)
      // - ceil is (n+1)/2
      // Actually we need floor and ceil of the rational value
      const floorVal = this.floor();
      const ceilVal = this.ceil();

      switch (mode) {
        case 'down':
          return floorVal;
        case 'up':
          return ceilVal;
        case 'toward':
          // Toward zero: use ceil for negative, floor for positive
          return this._numerator > 0n ? floorVal : ceilVal;
        case 'away':
          // Away from zero: use floor for negative, ceil for positive
          return this._numerator > 0n ? ceilVal : floorVal;
        case 'even':
          // Round toward even
          return (floorVal & 1n) === 0n ? floorVal : ceilVal;
        case 'odd':
          // Round toward odd
          return (floorVal & 1n) === 1n ? floorVal : ceilVal;
      }
    }

    // Not at a half - round to nearest
    const floorVal = this.floor();
    const ceilVal = this.ceil();

    // Distance to floor and ceil
    const distToFloor = this.sub(new Rational(floorVal, 1n)).abs();
    const distToCeil = new Rational(ceilVal, 1n).sub(this).abs();

    // Return the closer one
    if (distToFloor.lt(distToCeil)) {
      return floorVal;
    } else {
      return ceilVal;
    }
  }

  /**
   * Truncate toward zero (return integer part).
   */
  trunc(): bigint {
    return this._numerator / this._denominator;
  }

  /**
   * Check if this rational is zero.
   */
  isZero(): boolean {
    return this._numerator === 0n;
  }

  /**
   * Check if this rational is one.
   */
  isOne(): boolean {
    return this._numerator === 1n && this._denominator === 1n;
  }

  /**
   * Check if this rational is an integer (denominator is 1).
   */
  isInteger(): boolean {
    return this._denominator === 1n;
  }

  /**
   * Check if this rational is positive.
   */
  isPositive(): boolean {
    return this._numerator > 0n;
  }

  /**
   * Check if this rational is negative.
   */
  isNegative(): boolean {
    return this._numerator < 0n;
  }

  /**
   * Return the pair (numerator, denominator) as a tuple.
   */
  asIntegerRatio(): [bigint, bigint] {
    return [this._numerator, this._denominator];
  }

  /**
   * Return the height of this rational (max of |numerator|, denominator).
   */
  height(): bigint {
    const absNum = this._numerator < 0n ? -this._numerator : this._numerator;
    return absNum > this._denominator ? absNum : this._denominator;
  }

  // ============================================
  // Additional SageMath Rational methods (stubs)
  // ============================================

  /**
   * Return the continued fraction list of partial quotients.
   *
   * @param type - Either 'std' (standard) or 'hj' (Hirzebruch-Jung)
   * @returns List of partial quotients
   *
   * @example
   * ```typescript
   * new Rational(13n, 9n).continued_fraction_list()  // [1n, 2n, 4n]
   * // 13/9 = 1 + 1/(2 + 1/4)
   * ```
   *
   * @see Reference: sage/rings/rational.pyx:continued_fraction_list
   */
  continued_fraction_list(type: 'std' | 'hj' = 'std'): bigint[] {
    const result: bigint[] = [];
    let p = this._numerator;
    let q = this._denominator;

    // Helper: floor division (rounds toward -infinity)
    const fdiv = (a: bigint, b: bigint): [bigint, bigint] => {
      // For positive b: floor(a/b)
      if (b > 0n) {
        if (a >= 0n) {
          const quotient = a / b;
          return [quotient, a - quotient * b];
        } else {
          // For negative a, floor division: (a - (b-1)) / b
          const quotient = (a - b + 1n) / b;
          return [quotient, a - quotient * b];
        }
      } else {
        // b < 0: flip signs
        const quotient = a >= 0n ? (a - b - 1n) / b : a / b;
        return [quotient, a - quotient * b];
      }
    };

    // Helper: ceiling division (rounds toward +infinity)
    const cdiv = (a: bigint, b: bigint): [bigint, bigint] => {
      // For positive b: ceil(a/b)
      if (b > 0n) {
        if (a >= 0n) {
          const quotient = (a + b - 1n) / b;
          return [quotient, a - quotient * b];
        } else {
          const quotient = a / b;
          return [quotient, a - quotient * b];
        }
      } else {
        // b < 0
        const quotient = a >= 0n ? a / b : (a + b + 1n) / b;
        return [quotient, a - quotient * b];
      }
    };

    if (type === 'std') {
      // Standard continued fraction using floor division
      while (q !== 0n) {
        const [quotient, tmp] = fdiv(p, q);
        result.push(quotient);
        p = q;
        q = tmp;
      }
    } else if (type === 'hj') {
      // Hirzebruch-Jung continued fraction
      // Alternates: ceiling division, then floor division with negation
      while (q !== 0n) {
        // Ceiling division
        const [cq, ctmp] = cdiv(p, q);
        result.push(cq);
        p = q;
        q = ctmp;

        if (q === 0n) break;

        // Floor division, then negate the quotient
        const [fq, ftmp] = fdiv(p, q);
        result.push(-fq);
        p = q;
        q = ftmp;
      }
    } else {
      throw new ValueError("the type must be one of 'std', 'hj'");
    }

    return result;
  }

  /**
   * Return the continued fraction of this rational.
   *
   * This is an alias for continued_fraction_list with standard type.
   *
   * @see Reference: sage/rings/rational.pyx:continued_fraction
   */
  continued_fraction(): bigint[] {
    return this.continued_fraction_list('std');
  }

  /**
   * Return the valuation at prime p.
   *
   * The valuation is the power of p in the factorization of this rational.
   * For a/b, it's valuation(a, p) - valuation(b, p).
   *
   * @param p - A prime number (must be >= 2)
   * @returns The p-adic valuation, or 'Infinity' if self is zero
   *
   * @example
   * ```typescript
   * new Rational(-5n, 9n).valuation(5n)  // 1n (numerator has 5^1)
   * new Rational(-5n, 9n).valuation(3n)  // -2n (denominator has 3^2)
   * ```
   *
   * @see Reference: sage/rings/rational.pyx:valuation
   */
  valuation(p: IntegerLike): bigint | 'Infinity' {
    const pBig = toBigInt(p);
    if (pBig < 2n) {
      throw new ValueError('p must be at least 2');
    }

    if (this._numerator === 0n) {
      return 'Infinity';
    }

    const numVal = integer_valuation(this._numerator, pBig);
    const denVal = integer_valuation(this._denominator, pBig);
    return numVal - denVal;
  }

  /**
   * Alias for valuation.
   * @see Reference: sage/rings/rational.pyx:ord
   */
  ord(p: IntegerLike): bigint {
    const result = this.valuation(p);
    if (result === 'Infinity') {
      throw new ValueError('valuation of 0 is not defined');
    }
    return result;
  }

  /**
   * Return the local height at prime p.
   *
   * The local height at p is max(-valuation(self, p), 0) * log(p).
   *
   * @param p - A prime number
   * @param prec - Precision (not used in this implementation)
   * @returns The local height as a floating point number
   *
   * @example
   * ```typescript
   * new Rational(25n, 6n).local_height(2n)  // log(2) ≈ 0.693
   * new Rational(25n, 6n).local_height(3n)  // log(3) ≈ 1.099
   * ```
   *
   * @see Reference: sage/rings/rational.pyx:local_height
   */
  local_height(p: IntegerLike, prec?: number): number {
    const pBig = toBigInt(p);
    if (this._numerator === 0n) {
      return 0;
    }

    const val = this.valuation(pBig);
    if (val === 'Infinity' || val >= 0n) {
      return 0;
    }

    return Number(-val) * Math.log(Number(pBig));
  }

  /**
   * Return the global (absolute logarithmic) height of this rational.
   *
   * The height is max(log|numerator|, log|denominator|).
   *
   * @param prec - Precision (not used in this implementation)
   * @returns The global height as a floating point number
   *
   * @example
   * ```typescript
   * new Rational(6n, 25n).global_height()  // log(25) ≈ 3.219
   * new Rational(0n, 1n).global_height()  // 0
   * ```
   *
   * @see Reference: sage/rings/rational.pyx:global_height
   */
  global_height(prec?: number): number {
    const absNum = this._numerator < 0n ? -this._numerator : this._numerator;
    const maxVal = absNum > this._denominator ? absNum : this._denominator;
    return Math.log(Number(maxVal));
  }

  /**
   * Check if this rational is a perfect square.
   *
   * A rational a/b is a perfect square if both a and b are perfect squares
   * (with a >= 0).
   *
   * @returns true if this rational is a perfect square
   *
   * @example
   * ```typescript
   * new Rational(9n, 4n).is_square()   // true (= (3/2)^2)
   * new Rational(4n, 3n).is_square()   // false
   * new Rational(-1n, 4n).is_square()  // false (negative)
   * ```
   *
   * @see Reference: sage/rings/rational.pyx:is_square
   */
  is_square(): boolean {
    if (this._numerator < 0n) {
      return false;
    }
    if (this._numerator === 0n) {
      return true;
    }
    return integer_is_square(this._numerator) && integer_is_square(this._denominator);
  }

  /**
   * Return the square root if this rational is a perfect square.
   *
   * @param options - Options for square root computation
   * @param options.extend - If true, allow non-rational results (throws if false and not a square)
   * @param options.all - If true, return all square roots (both positive and negative)
   * @returns The square root(s)
   * @throws {ValueError} If not a perfect square and extend is false
   *
   * @example
   * ```typescript
   * new Rational(25n, 9n).sqrt()  // 5/3
   * new Rational(25n, 9n).sqrt({ all: true })  // [5/3, -5/3]
   * ```
   *
   * @see Reference: sage/rings/rational.pyx:sqrt
   */
  sqrt(options?: { extend?: boolean; all?: boolean }): Rational | Rational[] {
    const extend = options?.extend ?? true;
    const all = options?.all ?? false;

    if (this._numerator === 0n) {
      return all ? [this] : this;
    }

    if (this._numerator < 0n) {
      if (extend) {
        throw new NotImplementedError('Complex square roots not supported');
      }
      if (all) {
        return [];
      }
      throw new ValueError('square root of negative number not rational');
    }

    const numSqrt = isqrt(this._numerator);
    if (numSqrt * numSqrt !== this._numerator) {
      if (extend) {
        throw new NotImplementedError('Non-rational square roots not supported');
      }
      if (all) {
        return [];
      }
      throw new ValueError(`square root of ${this.toString()} not a rational number`);
    }

    const denSqrt = isqrt(this._denominator);
    if (denSqrt * denSqrt !== this._denominator) {
      if (extend) {
        throw new NotImplementedError('Non-rational square roots not supported');
      }
      if (all) {
        return [];
      }
      throw new ValueError(`square root of ${this.toString()} not a rational number`);
    }

    const result = new Rational(numSqrt, denSqrt);
    if (all) {
      return [result, result.neg()];
    }
    return result;
  }

  /**
   * Check if this rational is an n-th power.
   *
   * @param n - The power to check (can be negative)
   * @returns true if this rational is an n-th power
   *
   * @example
   * ```typescript
   * new Rational(25n, 4n).is_nth_power(2n)   // true
   * new Rational(125n, 8n).is_nth_power(3n)  // true
   * new Rational(-125n, 8n).is_nth_power(3n) // true
   * new Rational(9n, 2n).is_nth_power(2n)    // false
   * ```
   *
   * @see Reference: sage/rings/rational.pyx:is_nth_power
   */
  is_nth_power(n: IntegerLike): boolean {
    let nBig = toBigInt(n);
    if (nBig === 0n) {
      throw new ValueError('n cannot be zero');
    }

    if (nBig < 0n) {
      nBig = -nBig;
    }

    // Even roots of negative numbers don't exist in rationals
    if (nBig % 2n === 0n && this._numerator < 0n) {
      return false;
    }

    const [numIsNthPower] = integerNthRoot(
      this._numerator < 0n ? -this._numerator : this._numerator,
      nBig
    );
    if (!numIsNthPower) {
      return false;
    }

    const [denIsNthPower] = integerNthRoot(this._denominator, nBig);
    return denIsNthPower;
  }

  /**
   * Compute the n-th root of this rational, or raise an error if not a perfect power.
   *
   * @param n - The root to compute (can be negative)
   * @returns The n-th root
   * @throws {ValueError} If n is zero, or if not a perfect n-th power
   *
   * @example
   * ```typescript
   * new Rational(25n, 4n).nth_root(2n)   // 5/2
   * new Rational(125n, 8n).nth_root(3n)  // 5/2
   * new Rational(-125n, 8n).nth_root(3n) // -5/2
   * new Rational(25n, 4n).nth_root(-2n)  // 2/5
   * ```
   *
   * @see Reference: sage/rings/rational.pyx:nth_root
   */
  nth_root(n: IntegerLike): Rational {
    let nBig = toBigInt(n);
    if (nBig === 0n) {
      throw new ValueError('n cannot be zero');
    }

    let negative = false;
    if (nBig < 0n) {
      nBig = -nBig;
      negative = true;
    }

    // Even roots of negative numbers
    if (nBig % 2n === 0n && this._numerator < 0n) {
      throw new ValueError('cannot take even root of negative number');
    }

    const absNum = this._numerator < 0n ? -this._numerator : this._numerator;
    const numSign = this._numerator < 0n ? -1n : 1n;

    const [numExact, numRoot] = integerNthRoot(absNum, nBig);
    if (!numExact) {
      throw new ValueError(`not a perfect ${ordinalStr(nBig)} power`);
    }

    const [denExact, denRoot] = integerNthRoot(this._denominator, nBig);
    if (!denExact) {
      throw new ValueError(`not a perfect ${ordinalStr(nBig)} power`);
    }

    if (negative) {
      return new Rational(numSign * denRoot, numRoot);
    } else {
      return new Rational(numSign * numRoot, denRoot);
    }
  }

  /**
   * Return the period of the repeating part of the decimal expansion.
   *
   * For a rational n/d with (n, d) = 1, the period is the multiplicative
   * order of 10 modulo d (after removing factors of 2 and 5).
   *
   * @returns The period length
   *
   * @example
   * ```typescript
   * new Rational(1n, 7n).period()  // 6n (1/7 = 0.142857...)
   * new Rational(1n, 8n).period()  // 1n (1/8 = 0.125, terminates)
   * ```
   *
   * @see Reference: sage/rings/rational.pyx:period
   */
  period(): bigint {
    if (this._numerator === 0n) {
      return 1n;
    }

    // Remove factors of 2 and 5 from denominator
    let d = this._denominator;
    while (d % 2n === 0n) {
      d /= 2n;
    }
    while (d % 5n === 0n) {
      d /= 5n;
    }

    if (d === 1n) {
      return 1n; // Terminating decimal
    }

    // Find the multiplicative order of 10 mod d
    let t = 1n;
    let power = 10n % d;
    while (power !== 1n) {
      power = (power * 10n) % d;
      t++;
    }

    return t;
  }

  /**
   * Return the real part (self for rationals).
   * @see Reference: sage/rings/rational.pyx:real
   */
  real(): Rational {
    return this;
  }

  /**
   * Return the imaginary part (0 for rationals).
   * @see Reference: sage/rings/rational.pyx:imag
   */
  imag(): Rational {
    return Rational.zero();
  }

  /**
   * Return the conjugate (self for rationals).
   * @see Reference: sage/rings/rational.pyx:conjugate
   */
  conjugate(): Rational {
    return this;
  }

  /**
   * Return the norm (self^2 for rationals).
   * @see Reference: sage/rings/rational.pyx:norm
   */
  norm(): Rational {
    return this.mul(this);
  }

  /**
   * Return the trace (2*self for rationals).
   * @see Reference: sage/rings/rational.pyx:trace
   */
  trace(): Rational {
    return this.add(this);
  }

  /**
   * Return the minimal polynomial of this rational.
   *
   * For a rational r, the minimal polynomial is always (x - r).
   *
   * @param varName - Variable name (default: 'x')
   * @returns The minimal polynomial (x - r) over QQ
   *
   * @example
   * ```typescript
   * new Rational(1n, 3n).minpoly()  // x - 1/3
   * ```
   *
   * @see Reference: sage/rings/rational.pyx:minpoly
   */
  minpoly(varName: string = 'x'): Polynomial<Rational> {
    // QQ[varName]([-self, 1]) = (x - self)
    const QQRing = this._qqRing();
    const R = new PolynomialRing(QQRing, varName);
    // Polynomial is [-self, 1] which represents -self + x
    return R.__call__([this.neg(), QQRing.one()]);
  }

  /**
   * Helper to create a coefficient ring for QQ.
   */
  private _qqRing(): {
    zero: () => Rational;
    one: () => Rational;
    __call__: (x: unknown) => Rational;
    is_field: () => boolean;
  } {
    return {
      zero: () => new Rational(0n, 1n),
      one: () => new Rational(1n, 1n),
      __call__: (x: unknown): Rational => {
        if (x instanceof Rational) return x;
        if (typeof x === 'bigint') return new Rational(x, 1n);
        if (typeof x === 'number') return Rational.from(x);
        throw new ValueError(`cannot coerce ${x} to Rational`);
      },
      is_field: () => true,
    };
  }

  /**
   * Return the characteristic polynomial of this rational.
   *
   * For a rational r, the characteristic polynomial is always (x - r).
   *
   * @param varName - Variable name (default: 'x')
   * @returns The characteristic polynomial (x - r) over QQ
   *
   * @example
   * ```typescript
   * new Rational(2n, 1n).charpoly()  // x - 2
   * ```
   *
   * @see Reference: sage/rings/rational.pyx:charpoly
   */
  charpoly(varName: string = 'x'): Polynomial<Rational> {
    // For rationals, charpoly and minpoly are the same
    return this.minpoly(varName);
  }

  /**
   * Return the content of self and other.
   *
   * The content is the unique positive rational c such that self/c and other/c
   * are coprime integers.
   *
   * @param other - Another rational or list of rationals
   * @returns The content as a Rational
   *
   * @example
   * ```typescript
   * new Rational(2n, 3n).content(new Rational(2n, 3n))  // 2/3
   * new Rational(2n, 3n).content(new Rational(1n, 5n))  // 1/15
   * ```
   *
   * @see Reference: sage/rings/rational.pyx:content
   */
  content(other: Rational | Rational[]): Rational {
    const others = Array.isArray(other) ? other : [other];
    const all = [this, ...others];

    const nums = all.map((r) => r._numerator);
    const denoms = all.map((r) => r._denominator);

    const numGcd = gcd(nums.map((n) => (n < 0n ? -n : n)));
    const denomLcm = lcm(denoms);

    return new Rational(numGcd, denomLcm);
  }

  /**
   * Return self with all powers of all primes in S removed.
   *
   * @param S - List or tuple of primes
   * @returns The rational with all powers of primes in S removed
   *
   * @example
   * ```typescript
   * new Rational(3n, 4n).prime_to_S_part([2n])  // 3
   * new Rational(-3n, 4n).prime_to_S_part([3n]) // -1/4
   * new Rational(700n, 99n).prime_to_S_part([2n, 3n, 5n]) // 7/11
   * ```
   *
   * @see Reference: sage/rings/rational.pyx:prime_to_S_part
   */
  prime_to_S_part(S: bigint[] = []): Rational {
    if (this._numerator === 0n) {
      return this;
    }

    let result: Rational = this;
    for (const p of S) {
      const [, unit] = result.val_unit(p);
      result = unit;
    }
    return result;
  }

  /**
   * Return the valuation and p-adic unit part.
   *
   * @param p - A prime (must be >= 2)
   * @returns [valuation, unit] where self = unit * p^valuation
   *
   * @example
   * ```typescript
   * new Rational(-4n, 17n).val_unit(2n)  // [2n, -1/17]
   * new Rational(-4n, 17n).val_unit(17n) // [-1n, -4]
   * ```
   *
   * @see Reference: sage/rings/rational.pyx:val_unit
   */
  val_unit(p: bigint): [bigint | 'Infinity', Rational] {
    if (p < 2n) {
      throw new ValueError('p must be at least 2');
    }

    if (this._numerator === 0n) {
      return ['Infinity', Rational.one()];
    }

    const v = this.valuation(p);
    if (v === 'Infinity') {
      return ['Infinity', Rational.one()];
    }

    // unit = self / p^v
    const pPowV = bigintPow(p, v < 0n ? -v : v);
    let unit: Rational;
    if (v >= 0n) {
      unit = this.div(new Rational(pPowV, 1n));
    } else {
      unit = this.mul(new Rational(pPowV, 1n));
    }

    return [v, unit];
  }

  /**
   * Return the support (list of primes dividing numerator or denominator).
   *
   * @returns Sorted list of primes appearing in the factorization
   * @throws {ArithmeticError} If self is zero
   *
   * @example
   * ```typescript
   * new Rational(-4n, 17n).support()  // [2n, 17n]
   * ```
   *
   * @see Reference: sage/rings/rational.pyx:support
   */
  support(): bigint[] {
    if (this._numerator === 0n) {
      throw new ArithmeticError('Support of 0 not defined.');
    }

    const absNum = this._numerator < 0n ? -this._numerator : this._numerator;
    const numFactors = absNum === 1n ? [] : prime_factors(absNum);
    const denFactors = this._denominator === 1n ? [] : prime_factors(this._denominator);

    // Combine and deduplicate
    const combined = new Set([...numFactors, ...denFactors]);
    return [...combined].sort((a, b) => (a < b ? -1 : a > b ? 1 : 0));
  }

  /**
   * Return the GCD with another rational.
   *
   * For rationals a/b and c/d, gcd = gcd(a, c) / lcm(b, d).
   *
   * @param other - Another rational
   * @returns The GCD as a positive rational
   *
   * @see Reference: sage/rings/rational.pyx
   */
  rational_gcd(other: Rational): Rational {
    const numGcd = gcd(
      this._numerator < 0n ? -this._numerator : this._numerator,
      other._numerator < 0n ? -other._numerator : other._numerator
    );
    const denomLcm = lcm(this._denominator, other._denominator);
    return new Rational(numGcd, denomLcm);
  }

  /**
   * Return the LCM with another rational.
   *
   * For rationals a/b and c/d, lcm = lcm(a, c) / gcd(b, d).
   *
   * @param other - Another rational
   * @returns The LCM as a positive rational
   *
   * @see Reference: sage/rings/rational.pyx
   */
  rational_lcm(other: Rational): Rational {
    const numLcm = lcm(
      this._numerator < 0n ? -this._numerator : this._numerator,
      other._numerator < 0n ? -other._numerator : other._numerator
    );
    const denomGcd = gcd(this._denominator, other._denominator);
    return new Rational(numLcm, denomGcd);
  }

  /**
   * Return the factorial if this rational is a non-negative integer.
   *
   * @returns n! as a Rational
   * @throws {ValueError} If this is not a non-negative integer
   *
   * @example
   * ```typescript
   * new Rational(5n).factorial()  // 120
   * ```
   *
   * @see Reference: sage/rings/rational.pyx:factorial
   */
  factorial(): Rational {
    if (this._denominator !== 1n) {
      throw new ValueError('factorial only defined for integers');
    }
    if (this._numerator < 0n) {
      throw new ValueError('factorial not defined for negative integers');
    }
    return new Rational(integer_factorial(this._numerator), 1n);
  }

  /**
   * Return the Gamma function value.
   *
   * For positive integers n, gamma(n) = (n-1)!
   *
   * @see Reference: sage/rings/rational.pyx:gamma
   */
  gamma(): Rational {
    if (this._denominator !== 1n) {
      throw new NotImplementedError('gamma only implemented for integers');
    }
    if (this._numerator <= 0n) {
      throw new ValueError('gamma not defined for non-positive integers');
    }
    return new Rational(integer_factorial(this._numerator - 1n), 1n);
  }

  /**
   * Round to a given number of decimal digits.
   *
   * @param ndigits - Number of decimal places (default: 0)
   * @returns Rounded value as a Rational
   *
   * @example
   * ```typescript
   * new Rational(7n, 3n).roundToRational(2)  // 233/100 (≈ 2.33)
   * new Rational(7n, 3n).roundToRational(0)  // 2
   * ```
   *
   * @see Reference: sage/rings/rational.pyx
   */
  roundToRational(ndigits: number = 0): Rational {
    if (ndigits === 0) {
      return new Rational(this.round(), 1n);
    }

    const factor = bigintPow(10n, BigInt(ndigits));
    const scaled = this.mul(new Rational(factor, 1n));
    const rounded = scaled.round();
    return new Rational(rounded, factor);
  }

  /**
   * Determine if this rational is an S-unit.
   *
   * x is an S-unit if x.valuation(p) == 0 for all p not in S.
   *
   * @param S - List of primes (or undefined for empty S)
   * @returns true if this is an S-unit
   *
   * @example
   * ```typescript
   * new Rational(1n, 2n).is_S_unit([])     // false
   * new Rational(1n, 2n).is_S_unit([2n])   // true
   * ```
   *
   * @see Reference: sage/rings/rational.pyx:is_S_unit
   */
  is_S_unit(S?: bigint[]): boolean {
    const a = this.abs();
    if (a.eq(1n)) {
      return true;
    }
    if (!S || S.length === 0) {
      return false;
    }
    return a.prime_to_S_part(S).eq(1n);
  }

  /**
   * Determine if this rational is an S-integer.
   *
   * x is S-integral if x.valuation(p) >= 0 for all p not in S.
   *
   * @param S - List of primes
   * @returns true if this is S-integral
   *
   * @example
   * ```typescript
   * new Rational(1n, 2n).is_S_integral([])    // false
   * new Rational(1n, 2n).is_S_integral([2n])  // true
   * ```
   *
   * @see Reference: sage/rings/rational.pyx:is_S_integral
   */
  is_S_integral(S: bigint[] = []): boolean {
    if (this.is_integral()) {
      return true;
    }
    return this.prime_to_S_part(S).is_integral();
  }

  /**
   * Return the p-adic valuation.
   * @see Reference: sage/rings/rational.pyx:padic_valuation
   */
  padic_valuation(p: bigint): bigint | 'Infinity' {
    return this.valuation(p);
  }

  /**
   * Return the denominator valuation at p.
   *
   * @param p - A prime (must be >= 2)
   * @returns The valuation of the denominator at p
   *
   * @see Reference: sage/rings/rational.pyx:denominator_valuation
   */
  denominator_valuation(p: bigint): bigint {
    if (p < 2n) {
      throw new ValueError('p must be at least 2');
    }
    return integer_valuation(this._denominator, p);
  }

  /**
   * Return the numerator valuation at p.
   *
   * @param p - A prime (must be >= 2)
   * @returns The valuation of the numerator at p (or Infinity if numerator is 0)
   *
   * @see Reference: sage/rings/rational.pyx:numerator_valuation
   */
  numerator_valuation(p: bigint): bigint | 'Infinity' {
    if (p < 2n) {
      throw new ValueError('p must be at least 2');
    }
    if (this._numerator === 0n) {
      return 'Infinity';
    }
    return integer_valuation(this._numerator < 0n ? -this._numerator : this._numerator, p);
  }

  /**
   * Return the sign as an integer.
   * @see Reference: sage/rings/rational.pyx:sign
   */
  signAsInt(): bigint {
    return this.sign;
  }

  /**
   * Convert to a string in given base.
   *
   * @param base - The base (default: 10, must be between 2 and 36)
   * @returns String representation in the given base
   *
   * @example
   * ```typescript
   * new Rational(-4n, 17n).str()    // "-4/17"
   * new Rational(-4n, 17n).str(2)   // "-100/10001"
   * ```
   *
   * @see Reference: sage/rings/rational.pyx:str
   */
  str(base: number = 10): string {
    if (base < 2 || base > 36) {
      throw new ValueError(`base (=${base}) must be between 2 and 36`);
    }

    const numStr = this._numerator.toString(base);
    if (this._denominator === 1n) {
      return numStr;
    }
    const denStr = this._denominator.toString(base);
    return `${numStr}/${denStr}`;
  }

  /**
   * Return the number of digits in the numerator plus denominator.
   *
   * @param base - The base (default: 10)
   * @returns Sum of digit counts for numerator and denominator
   *
   * @example
   * ```typescript
   * new Rational(123n, 45n).ndigits()  // 5n (3 + 2 digits)
   * ```
   *
   * @see Reference: sage/rings/rational.pyx:ndigits
   */
  ndigits(base: IntegerLike = 10n): bigint {
    const baseBig = toBigInt(base);
    if (baseBig < 2n) {
      throw new ValueError('base must be at least 2');
    }

    const absNum = this._numerator < 0n ? -this._numerator : this._numerator;

    // Count digits in numerator
    let numDigits = 0n;
    if (absNum === 0n) {
      numDigits = 1n;
    } else {
      let temp = absNum;
      while (temp > 0n) {
        temp /= baseBig;
        numDigits++;
      }
    }

    // Count digits in denominator
    let denDigits = 0n;
    let temp = this._denominator;
    while (temp > 0n) {
      temp /= baseBig;
      denDigits++;
    }

    return numDigits + denDigits;
  }

  /**
   * Return the number of bits needed to represent this rational.
   *
   * This is the sum of the bit lengths of numerator and denominator.
   *
   * @returns Total number of bits
   *
   * @example
   * ```typescript
   * new Rational(8n, 4n).nbits()  // reduced to 2/1, so nbits of 2 + nbits of 1
   * ```
   *
   * @see Reference: sage/rings/rational.pyx:nbits
   */
  nbits(): bigint {
    const absNum = this._numerator < 0n ? -this._numerator : this._numerator;

    // Bit length of numerator
    let numBits: bigint;
    if (absNum === 0n) {
      numBits = 0n;
    } else {
      numBits = BigInt(absNum.toString(2).length);
    }

    // Bit length of denominator
    const denBits = BigInt(this._denominator.toString(2).length);

    return numBits + denBits;
  }

  /**
   * Return the list [numerator, denominator].
   * @see Reference: sage/rings/rational.pyx:list
   */
  list(): bigint[] {
    return [this._numerator, this._denominator];
  }

  /**
   * Test if this rational is a unit (always true for non-zero rationals).
   * @see Reference: sage/rings/rational.pyx:is_unit
   */
  is_unit(): boolean {
    return this._numerator !== 0n;
  }

  /**
   * Test if this is an integral element (denominator is 1).
   * @see Reference: sage/rings/rational.pyx:is_integral
   */
  is_integral(): boolean {
    return this._denominator === 1n;
  }

  /**
   * Return the additive order (infinity for non-zero, 1 for zero).
   * @see Reference: sage/rings/rational.pyx:additive_order
   */
  additive_order(): bigint | string {
    if (this._numerator === 0n) {
      return 1n;
    }
    return 'Infinity';
  }

  /**
   * Return the multiplicative order (infinity for |self| != 1).
   * @see Reference: sage/rings/rational.pyx:multiplicative_order
   */
  multiplicative_order(): bigint | string {
    if (this._numerator === 1n && this._denominator === 1n) {
      return 1n;
    }
    if (this._numerator === -1n && this._denominator === 1n) {
      return 2n;
    }
    return 'Infinity';
  }

  /**
   * Simplify and return self (rationals are always simplified).
   * @see Reference: sage/rings/rational.pyx:simplify
   */
  simplify(): Rational {
    return this;
  }

  /**
   * Return a floating point approximation with given precision.
   * @see Reference: sage/rings/rational.pyx:n
   */
  n(prec?: number): number {
    return this.toNumber();
  }

  /**
   * Alias for n().
   * @see Reference: sage/rings/rational.pyx:numerical_approx
   */
  numerical_approx(prec?: number): number {
    return this.n(prec);
  }
}
