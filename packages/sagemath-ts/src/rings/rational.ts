import { mpz_remove } from '../types/gmp.js';
/**
 * @module sage/rings/rational
 * @description Rational numbers
 *
 * Port of: sage/rings/rational.pyx
 * Reference: reference/sage/src/sage/rings/rational.pyx
 */

import {
  gcd,
  factorial as integer_factorial,
  is_square as integer_is_square,
  isqrt,
  lcm,
  prime_factors,
} from '../arith/misc.js';
import {
  ArithmeticError,
  NotImplementedError,
  OverflowError,
  TypeError,
  ValueError,
  ZeroDivisionError,
} from '../errors.js';
import { type IntegerLike, toBigInt } from '../types/coercion.js';
import { UnsignedInfinity, type UnsignedInfinityElement } from './complex_mpfr.js';
import { IntegerModRing } from './finite_rings/integer_mod_ring.js';
import { Integer } from './integer_ring.js';
import type { Polynomial, PolynomialRingBase } from './polynomial/polynomial_element.js';
import { PolynomialRing } from './polynomial/polynomial_ring.js';
import { RR } from './real_mpfr.js';
import { RDF } from './real_double.js';

/** The rational.pyx root methods take a signed C int exponent. */
function rootExponent(n: IntegerLike): bigint {
  const value = toBigInt(n);
  if (value < -2147483648n || value > 2147483647n) {
    throw new OverflowError('value too large to convert to int');
  }
  return value;
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
  // Use bigint modulo to avoid precision loss for large numbers
  const lastDigit = Number(n % 10n);
  const lastTwoDigits = Number(n % 100n);

  if (lastTwoDigits >= 11 && lastTwoDigits <= 13) {
    return `${n}th`;
  }

  switch (lastDigit) {
    case 1:
      return `${n}st`;
    case 2:
      return `${n}nd`;
    case 3:
      return `${n}rd`;
    default:
      return `${n}th`;
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
   * @param numerator - The numerator (default: zero)
   * @param denominator - The denominator (must be non-zero)
   * @throws {ValueError} If denominator is zero
   */
  constructor(numerator: IntegerLike = 0n, denominator: IntegerLike = 1n) {
    let num = toBigInt(numerator);
    let den = toBigInt(denominator);

    if (den === 0n) {
      throw new ValueError('denominator must not be 0');
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
   * @param value - IntegerLike, number, string, Rational, boolean, or null (default: zero)
   * @returns A new Rational
   */
  static from(value: number | IntegerLike | string | Rational | boolean | null = 0n): Rational {
    if (value === null) return Rational.zero();
    if (typeof value === 'boolean') return new Rational(value ? 1n : 0n);
    if (value instanceof Integer) return new Rational(value.value);
    if (value instanceof Rational) {
      return value;
    }

    if (typeof value === 'bigint') {
      return new Rational(value, 1n);
    }

    if (typeof value === 'number') {
      // rational.pyx:681 -> RealNumber -> simplest_rational (or exact integer).
      if (Number.isInteger(value)) return new Rational(BigInt(value), 1n);
      const [numerator, denominator] = RR().__call__(value).simplest_rational();
      return new Rational(numerator, denominator);
    }

    if (typeof value === 'string') {
      return Rational.fromString(value);
    }

    throw new TypeError(`cannot convert ${typeof value} to Rational`);
  }

  /**
   * Create a Rational from a string.
   *
   * Accepts integer and quotient strings using Sage/GMP base-prefix rules.
   * Decimal strings such as "1.5" are rejected; Rational.from(1.5) accepts a float.
   *
   * @param str - The string representation
   * @returns A new Rational
   */
  static fromString(str: string): Rational {
    // rational.pyx:630-639 uses mpq_set_str with base=0. GMP accepts integer
    // numerators/denominators, ASCII whitespace between digits and base
    // prefixes; it rejects decimal floats, plus signs and digit separators.
    const invalid = (): never => {
      throw new TypeError(`unable to convert '${str}' to a rational`);
    };
    const parseInteger = (part: string): bigint => {
      let text = part.replace(/^[ \t\n\r\v\f]+|[ \t\n\r\v\f]+$/g, '');
      const negative = text.startsWith('-');
      if (negative) text = text.slice(1);
      if (!/^[0-9]/.test(text)) return invalid();
      let prefix = '';
      let pattern = /^[0-9]+$/;
      if (/^0[xX]/.test(text)) {
        prefix = '0x';
        pattern = /^[0-9a-fA-F]+$/;
        text = text.slice(2);
      } else if (/^0[bB]/.test(text)) {
        prefix = '0b';
        pattern = /^[01]+$/;
        text = text.slice(2);
      } else if (text.startsWith('0')) {
        prefix = '0o';
        pattern = /^[0-7]+$/;
      }
      const digits = text.replace(/[ \t\n\r\v\f]/g, '');
      if (!pattern.test(digits)) return invalid();
      const value = BigInt(prefix + digits);
      return negative ? -value : value;
    };
    const parts = str.split('/');
    if (parts.length > 2) return invalid();
    const numerator = parseInteger(parts[0]!);
    const denominator = parts.length === 2 ? parseInteger(parts[1]!) : 1n;
    if (denominator === 0n) return invalid();
    return new Rational(numerator, denominator);
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
      return new Rational(
        this._numerator * other._numerator,
        this._denominator * other._denominator
      );
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
      return new Rational(
        this._numerator * other._denominator,
        this._denominator * other._numerator
      );
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
   * Host numbers compare after Sage's coercion to the Real Double Field.
   * @see Reference: real_double.pyx:RealDoubleElement._richcmp_
   */
  eq(other: Rational | IntegerLike | number): boolean {
    if (other instanceof Rational) {
      return this._numerator === other._numerator && this._denominator === other._denominator;
    }
    // Sage coerces QQ and Python float into RDF before comparing, including
    // rounded integral values, subnormals, infinities, and NaN.
    if (typeof other === 'number') return RDF.__call__(this).eq(other);
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
   * Return the nearest binary64 value, rounding ties to even.
   * @see Reference: sage/rings/rational.pyx:3898 (mpq_get_d_nearest)
   */
  toNumber(): number {
    if (this._numerator === 0n) return 0;
    const negative = this._numerator < 0n;
    const a = negative ? -this._numerator : this._numerator;
    const b = this._denominator;
    const sa = a.toString(2).length;
    const sb = b.toString(2).length;
    if (sa <= 53 && sb <= 53) return Number(this._numerator) / Number(b);

    // Sage keeps 54 or 55 quotient bits and a sticky remainder bit. Both
    // integer divisions must contribute to that remainder before rounding.
    let shift = sa - sb - 54;
    if (shift <= -1130) return negative ? -0 : 0;
    if (shift >= 971) return negative ? Number.NEGATIVE_INFINITY : Number.POSITIVE_INFINITY;
    let q: bigint;
    let remainderIsZero = true;
    if (shift > 0) {
      const divisor = 1n << BigInt(shift);
      remainderIsZero = a % divisor === 0n;
      q = a / divisor;
    } else {
      q = a << BigInt(-shift);
    }
    remainderIsZero = remainderIsZero && q % b === 0n;
    q /= b;

    let addShift = q < 1n << 54n ? 0 : 1;
    // Keep one rounding bit at the subnormal boundary; avoid double rounding.
    if (shift + addShift < -1075) addShift = -1075 - shift;
    if (addShift !== 0) {
      shift += addShift;
      const mask = (1n << BigInt(addShift)) - 1n;
      remainderIsZero = remainderIsZero && (q & mask) === 0n;
      q >>= BigInt(addShift);
    }
    if ((q & 1n) !== 0n) q += remainderIsZero ? (q & 2n) - 1n : 1n;

    // q is now even and exactly representable. Fold its low zero bit into
    // the exponent so JS never evaluates the unrepresentable 2**-1075.
    const significand = Number(q >> 1n);
    return (negative ? -significand : significand) * 2 ** (shift + 1);
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
      throw new ValueError(
        'You can only compute the valuation with respect to a integer larger than 1.'
      );
    }

    if (this._numerator === 0n) {
      return 'Infinity';
    }

    // Both components are nonzero here, so their valuations are finite.
    const numVal = new Integer(this._numerator).valuation(pBig) as bigint;
    const denVal = new Integer(this._denominator).valuation(pBig) as bigint;
    return numVal - denVal;
  }

  /**
   * Alias for valuation.
   * @see Deviation: Infinity Representation
   * @see Reference: sage/rings/rational.pyx:ord
   */
  ord(p: IntegerLike): bigint | 'Infinity' {
    return this.valuation(p);
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
   * @see Deviation: Rational Comparative Adapters
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

    return Number(-val) * new Integer(pBig).real_log();
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
   * @see Deviation: Rational Comparative Adapters
   * @see Reference: sage/rings/rational.pyx:global_height
   */
  global_height(prec?: number): number {
    const absNum = this._numerator < 0n ? -this._numerator : this._numerator;
    const maxVal = absNum > this._denominator ? absNum : this._denominator;
    return new Integer(maxVal).real_log();
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
    let exponent = rootExponent(n);
    if (exponent === 0n) throw new ValueError('n cannot be zero');
    // Preserve the C int minimum's wraparound in the vendored implementation.
    if (exponent < 0n) exponent = BigInt.asIntN(32, -exponent);
    if (exponent % 2n === 0n && this._numerator < 0n) return false;
    return (
      new Integer(this._numerator).nth_root(exponent, true)[1] &&
      new Integer(this._denominator).nth_root(exponent, true)[1]
    );
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
    const signedExponent = rootExponent(n);
    if (signedExponent === 0n) throw new ValueError('n cannot be zero');
    const inverse = signedExponent < 0n;
    const exponent = inverse ? BigInt.asIntN(32, -signedExponent) : signedExponent;
    const [num, numExact] = new Integer(this._numerator).nth_root(exponent, true);
    if (!numExact) throw new ValueError(`not a perfect ${ordinalStr(exponent)} power`);
    const [den, denExact] = new Integer(this._denominator).nth_root(exponent, true);
    if (!denExact) throw new ValueError(`not a perfect ${ordinalStr(exponent)} power`);
    if (inverse) {
      if (num.value === 0n) throw new ZeroDivisionError('rational division by zero');
      return new Rational(den, num);
    }
    return new Rational(num, den);
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
   * @see Reference: sage/rings/rational.pyx:2034 (period)
   */
  period(): bigint {
    // Sage: d = self.denominator(); d = d.val_unit(2)[1]; d = d.val_unit(5)[1]
    //       return Mod(10, d).multiplicative_order()
    let d = this._denominator;
    while (d % 2n === 0n) {
      d /= 2n;
    }
    while (d % 5n === 0n) {
      d /= 5n;
    }

    // Mod(10, 1) is the zero ring: its unique element has order 1.
    if (d === 1n) {
      return 1n;
    }

    // d is positive; the modular-element order delegates to PARI znorder.
    return new IntegerModRing(d).__call__(10n).multiplicative_order();
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
   * Return the norm from QQ to QQ of self, which is just self.
   *
   * Added for compatibility with NumberField.
   *
   * @see Reference: sage/rings/rational.pyx:2880 (norm)
   */
  norm(): Rational {
    return this;
  }

  /**
   * Return the norm from QQ to QQ of self, which is just self.
   * @see Reference: sage/rings/rational.pyx:2898 (relative_norm)
   */
  relative_norm(): Rational {
    return this;
  }

  /**
   * Return the norm from QQ to QQ of self, which is just self.
   * @see Reference: sage/rings/rational.pyx:2913 (absolute_norm)
   */
  absolute_norm(): Rational {
    return this;
  }

  /**
   * Return the trace from QQ to QQ of self, which is just self.
   *
   * Added for compatibility with NumberField.
   *
   * @see Reference: sage/rings/rational.pyx:2928 (trace)
   */
  trace(): Rational {
    return this;
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
   * @see Deviation: Valuation dispatch and native GMP factor removal
   */
  val_unit(p: IntegerLike): [bigint | 'Infinity', Rational] {
    const prime = toBigInt(p);
    if (prime < 2n) throw new ValueError('p must be at least 2.');
    if (this._numerator === 0n) return ['Infinity', Rational.one()];
    const [numeratorVal, numeratorUnit] = mpz_remove(this._numerator, prime);
    if (numeratorVal !== 0n) {
      return [numeratorVal, new Rational(numeratorUnit, this._denominator)];
    }
    const [denominatorVal, denominatorUnit] = mpz_remove(this._denominator, prime);
    return [-denominatorVal, new Rational(numeratorUnit, denominatorUnit)];
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
   * @see Deviation: Exact Return Types and Numeric Backends
   * @see Reference: sage/rings/rational.pyx:gamma
   */
  gamma(): Rational | UnsignedInfinityElement {
    if (this._denominator !== 1n) {
      throw new NotImplementedError('gamma only implemented for integers');
    }
    if (this._numerator <= 0n) {
      return UnsignedInfinity;
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

    const factor = new Rational(10n).pow(BigInt(ndigits));
    const rounded = this.mul(factor).round();
    return new Rational(rounded).div(factor);
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
   * @see Reference: sage/rings/rational.pyx:valuation
   */
  padic_valuation(p: IntegerLike): bigint | 'Infinity' {
    return this.valuation(p);
  }

  /**
   * Return the denominator valuation at p.
   *
   * @param p - A prime (must be >= 2)
   * @returns The valuation of the denominator at p
   *
   * @see Reference: sage/rings/integer.pyx:valuation (on the denominator)
   */
  denominator_valuation(p: IntegerLike): bigint {
    // A rational denominator is nonzero, hence this valuation is finite.
    return new Integer(this._denominator).valuation(p) as bigint;
  }

  /**
   * Return the numerator valuation at p.
   *
   * @param p - A prime (must be >= 2)
   * @returns The valuation of the numerator at p (or Infinity if numerator is 0)
   *
   * @see Reference: sage/rings/integer.pyx:valuation (on the numerator)
   */
  numerator_valuation(p: IntegerLike): bigint | 'Infinity' {
    return new Integer(this._numerator).valuation(p);
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
   * Return a list with this rational in it, for compatibility with the method
   * for number fields: the list ``[self]``.
   *
   * @see Reference: sage/rings/rational.pyx:709 (list)
   */
  list(): Rational[] {
    return [this];
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
