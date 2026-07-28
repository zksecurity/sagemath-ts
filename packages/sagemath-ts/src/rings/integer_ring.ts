/**
 * @module sage/rings/integer_ring
 * @description The ring ZZ of integers
 *
 * Port of: sage/rings/integer_ring.pyx
 * Reference: reference/sage/src/sage/rings/integer_ring.pyx
 */

import {
  type Factorization,
  binomial as _binomial,
  carmichael_lambda as _carmichael_lambda,
  divisors as _divisors,
  euler_phi as _euler_phi,
  factor as _factor,
  factorial as _factorial,
  fibonacci as _fibonacci,
  gcd as _gcd,
  inverse_mod as _inverse_mod,
  is_prime as _is_prime,
  is_prime_power as _is_prime_power,
  is_pseudoprime as _is_pseudoprime,
  is_squarefree as _is_squarefree,
  is_strong_probable_prime as _is_strong_probable_prime,
  jacobi_symbol as _jacobi_symbol,
  kronecker_symbol as _kronecker_symbol,
  lcm as _lcm,
  legendre_symbol as _legendre_symbol,
  lucas_number as _lucas_number,
  moebius as _moebius,
  next_prime as _next_prime,
  next_prime_power as _next_prime_power,
  nth_prime as _nth_prime,
  number_of_divisors as _number_of_divisors,
  power_mod as _power_mod,
  previous_prime as _previous_prime,
  prime_factors as _prime_factors,
  prime_to_m_part as _prime_to_m_part,
  primitive_root as _primitive_root,
  radical as _radical,
  sigma as _sigma,
  sqrt_mod as _sqrt_mod,
  squarefree_part as _squarefree_part,
  trial_division as _trial_division,
  xgcd as _xgcd,
  isqrt,
} from '../arith/misc.js';
import {
  ArithmeticError,
  NotImplementedError,
  TypeError as SageTypeError,
  ValueError,
  ZeroDivisionError,
} from '../errors.js';
import { SAGE_RAND_MAX, current_randstate } from '../misc/randstate.js';
import { DiscreteGaussianDistributionIntegerSampler } from '../stats/distributions/discrete_gaussian_integer.js';
import { type IntegerLike, toBigInt } from '../types/coercion.js';

let prevGaussianSampler: {
  sigma: number;
  sampler: DiscreteGaussianDistributionIntegerSampler;
} | null = null;

/**
 * Floor division of bigints, mirroring GMP's ``mpz_fdiv_q``.
 *
 * BigInt's ``/`` truncates towards zero; GMP (and therefore Sage's
 * ``Integer.__floordiv__`` / ``quo_rem``) rounds towards `-oo`.
 *
 * @internal
 */
function fdiv_q(a: bigint, b: bigint): bigint {
  let q = a / b;
  if (a % b !== 0n && a < 0n !== b < 0n) {
    q -= 1n;
  }
  return q;
}

/**
 * Remainder of floor division, mirroring GMP's ``mpz_fdiv_r``: the result is
 * zero or has the same sign as ``b``.
 *
 * @internal
 */
function fdiv_r(a: bigint, b: bigint): bigint {
  const r = a % b;
  if (r !== 0n && r < 0n !== b < 0n) {
    return r + b;
  }
  return r;
}

/**
 * The ring of integers ZZ.
 *
 * This is a singleton class representing the mathematical ring of integers.
 */
export class IntegerRing {
  private static instance: IntegerRing;

  private constructor() {}

  static getInstance(): IntegerRing {
    if (!IntegerRing.instance) {
      IntegerRing.instance = new IntegerRing();
    }
    return IntegerRing.instance;
  }

  /**
   * Coerce a value to an integer.
   */
  __call__(x: bigint | number | string): bigint {
    if (typeof x === 'bigint') {
      return x;
    }
    if (typeof x === 'number') {
      if (!Number.isInteger(x)) {
        throw new SageTypeError('cannot convert non-integer to Integer');
      }
      return BigInt(x);
    }
    if (typeof x === 'string') {
      return BigInt(x);
    }
    throw new SageTypeError(`cannot convert ${typeof x} to Integer`);
  }

  /**
   * Return zero element.
   */
  zero(): bigint {
    return 0n;
  }

  /**
   * Return one element.
   */
  one(): bigint {
    return 1n;
  }

  /**
   * Return the characteristic of ZZ (which is 0).
   */
  characteristic(): bigint {
    return 0n;
  }

  /**
   * Check if ZZ is a field (it's not).
   */
  is_field(): boolean {
    return false;
  }

  /**
   * Check if ZZ is a ring (it is).
   */
  is_ring(): boolean {
    return true;
  }

  /**
   * Check if ZZ is a domain (it is).
   */
  is_integral_domain(): boolean {
    return true;
  }

  /**
   * Return a random integer (Sage-compatible semantics).
   *
   * If only x is given, returns an integer in [0, x).
   * If x and y are given, returns an integer in [x, y).
   * If distribution is '1/n' or no bounds are given, uses Sage's default distribution.
   */
  random_element(
    x?: bigint | number,
    y?: bigint | number,
    distribution?: 'uniform' | 'mpz_rrandomb' | '1/n' | 'gaussian'
  ): bigint {
    // Mirror Sage's argument normalization rules.
    if (distribution === '1/n') {
      x = undefined;
      y = undefined;
    } else if (distribution === 'mpz_rrandomb' || distribution === 'gaussian') {
      y = undefined;
    }

    const rstate = current_randstate();

    if (distribution === 'gaussian') {
      if (x === undefined) {
        throw new ValueError("must specify x to use 'distribution=gaussian'");
      }
      const sigma = typeof x === 'number' ? x : Number(x);
      if (!Number.isFinite(sigma) || sigma <= 0) {
        throw new SageTypeError('x must be > 0');
      }
      if (prevGaussianSampler?.sigma === sigma) {
        return prevGaussianSampler.sampler.sample();
      }
      const sampler = new DiscreteGaussianDistributionIntegerSampler({ sigma });
      prevGaussianSampler = { sigma, sampler };
      return sampler.sample();
    }

    const xVal = x !== undefined ? this.__call__(x) : undefined;
    const yVal = y !== undefined ? this.__call__(y) : undefined;

    if (xVal !== undefined && yVal === undefined && xVal <= 0n) {
      throw new SageTypeError('x must be > 0');
    }
    if (xVal !== undefined && yVal !== undefined && xVal >= yVal) {
      throw new SageTypeError('x must be < y');
    }

    if ((distribution === undefined && xVal === undefined) || distribution === '1/n') {
      const half = Math.floor(SAGE_RAND_MAX / 2);
      let den = BigInt(rstate.c_random() - half);
      if (den === 0n) {
        den = 1n;
      }
      const numerator = BigInt(Math.floor((SAGE_RAND_MAX / 5) * 2));
      return numerator / den;
    }

    if (distribution === undefined || distribution === 'uniform') {
      if (yVal === undefined) {
        if (xVal === undefined) {
          return BigInt((rstate.c_random() % 5) - 2);
        }
        return rstate.random_below(xVal);
      }
      let nMin = xVal as bigint;
      let nWidth = (yVal as bigint) - nMin;
      if (nWidth <= 0n) {
        nMin = -2n;
        nWidth = 5n;
      }
      return nMin + rstate.random_below(nWidth);
    }

    if (distribution === 'mpz_rrandomb') {
      if (xVal === undefined) {
        throw new ValueError("must specify x to use 'distribution=mpz_rrandomb'");
      }
      const bits = Number(xVal);
      if (!Number.isFinite(bits) || bits < 0) {
        throw new SageTypeError('x must be >= 0');
      }
      return rstate.random_bits(bits);
    }

    throw new ValueError(`Unknown distribution for the integers: ${distribution}`);
  }

  toString(): string {
    return 'Integer Ring';
  }
}

/**
 * The ring of integers ZZ.
 */
export const ZZ = IntegerRing.getInstance();

/**
 * Integer class wrapping bigint with SageMath-compatible methods.
 *
 * This provides method-based access to arithmetic functions,
 * matching SageMath's Integer class interface.
 */
export class Integer {
  readonly value: bigint;

  constructor(value: bigint | number | string) {
    if (typeof value === 'bigint') {
      this.value = value;
    } else if (typeof value === 'number') {
      if (!Number.isInteger(value)) {
        throw new SageTypeError('cannot convert non-integer to Integer');
      }
      this.value = BigInt(value);
    } else {
      this.value = BigInt(value);
    }
  }

  /**
   * Return the absolute value.
   */
  abs(): Integer {
    return new Integer(this.value < 0n ? -this.value : this.value);
  }

  /**
   * Return the sign of this integer (-1, 0, or 1).
   */
  sign(): bigint {
    if (this.value < 0n) return -1n;
    if (this.value > 0n) return 1n;
    return 0n;
  }

  /**
   * Return the GCD of this integer and n.
   */
  gcd(n: Integer | bigint): Integer {
    const other = n instanceof Integer ? n.value : n;
    return new Integer(_gcd(this.value, other));
  }

  /**
   * Return the LCM of this integer and n.
   */
  lcm(n: Integer | bigint): Integer {
    const other = n instanceof Integer ? n.value : n;
    return new Integer(_lcm(this.value, other));
  }

  /**
   * Return the extended GCD of this integer and n.
   */
  xgcd(n: Integer | bigint): [Integer, Integer, Integer] {
    const other = n instanceof Integer ? n.value : n;
    const [g, s, t] = _xgcd(this.value, other);
    return [new Integer(g), new Integer(s), new Integer(t)];
  }

  /**
   * Return the prime factorization.
   *
   * @see Deviation: PARI Factorization Algorithms Limited (parigp-ts)
   */
  factor(): Factorization {
    return _factor(this.value);
  }

  /**
   * Test if this integer is prime.
   */
  is_prime(): boolean {
    return _is_prime(this.value);
  }

  /**
   * Test if this integer is a unit (±1).
   */
  is_unit(): boolean {
    return this.value === 1n || this.value === -1n;
  }

  /**
   * Return the integer square root.
   */
  isqrt(): Integer {
    if (this.value < 0n) {
      throw new ValueError('isqrt() argument must be nonnegative');
    }
    return new Integer(isqrt(this.value));
  }

  /**
   * Test if this integer is a perfect square.
   */
  is_square(): boolean {
    if (this.value < 0n) {
      return false;
    }
    const s = isqrt(this.value);
    return s * s === this.value;
  }

  /**
   * Return this mod n.
   *
   * Mirrors GMP's ``mpz_fdiv_r`` (as used by ``Integer.__mod__``): the
   * remainder is zero or has the same sign as the divisor.
   *
   * @see Reference: sage/rings/integer.pyx:__mod__
   */
  mod(n: Integer | bigint): Integer {
    const other = n instanceof Integer ? n.value : n;
    if (other === 0n) {
      throw new ZeroDivisionError('Integer modulo by zero');
    }
    return new Integer(fdiv_r(this.value, other));
  }

  /**
   * Return this divided by n (floor division).
   *
   * @see Reference: sage/rings/integer.pyx:_floordiv_
   */
  div(n: Integer | bigint): Integer {
    const other = n instanceof Integer ? n.value : n;
    if (other === 0n) {
      throw new ZeroDivisionError('Integer division by zero');
    }
    return new Integer(fdiv_q(this.value, other));
  }

  /**
   * Return quotient and remainder of this divided by n.
   *
   * As in Sage (``mpz_fdiv_qr``) the remainder returned is always either zero
   * or of the same sign as ``other``.
   *
   * @see Reference: sage/rings/integer.pyx:quo_rem
   */
  quo_rem(n: Integer | bigint): [Integer, Integer] {
    const other = n instanceof Integer ? n.value : n;
    if (other === 0n) {
      throw new ZeroDivisionError('Integer division by zero');
    }
    const q = fdiv_q(this.value, other);
    return [new Integer(q), new Integer(this.value - q * other)];
  }

  /**
   * Return the number of digits in base b.
   */
  ndigits(b: IntegerLike = 10n): bigint {
    const bBig = toBigInt(b);
    if (bBig <= 1n) {
      throw new ValueError('base must be at least 2');
    }

    let n = this.value < 0n ? -this.value : this.value;
    if (n === 0n) {
      // Sage returns ``self`` (i.e. 0) for zero -- integer.pyx:1833
      return 0n;
    }

    let count = 0n;
    while (n > 0n) {
      n /= bBig;
      count++;
    }
    return count;
  }

  /**
   * Return the number of bits needed to represent this integer.
   */
  nbits(): bigint {
    let n = this.value < 0n ? -this.value : this.value;
    if (n === 0n) {
      return 0n;
    }

    let count = 0n;
    while (n > 0n) {
      n >>= 1n;
      count++;
    }
    return count;
  }

  /**
   * Return the valuation of this integer at prime p.
   * This is the largest k such that p^k divides this integer.
   */
  valuation(p: Integer | bigint): bigint {
    const prime = p instanceof Integer ? p.value : p;

    if (prime <= 1n) {
      throw new ValueError('p must be at least 2');
    }

    let n = this.value < 0n ? -this.value : this.value;
    if (n === 0n) {
      throw new ValueError('valuation of 0 is not defined');
    }

    let k = 0n;
    while (n % prime === 0n) {
      n /= prime;
      k++;
    }
    return k;
  }

  // Arithmetic operations returning Integer
  add(n: Integer | bigint): Integer {
    const other = n instanceof Integer ? n.value : n;
    return new Integer(this.value + other);
  }

  sub(n: Integer | bigint): Integer {
    const other = n instanceof Integer ? n.value : n;
    return new Integer(this.value - other);
  }

  mul(n: Integer | bigint): Integer {
    const other = n instanceof Integer ? n.value : n;
    return new Integer(this.value * other);
  }

  neg(): Integer {
    return new Integer(-this.value);
  }

  pow(n: IntegerLike): Integer {
    const nBig = toBigInt(n);
    if (nBig < 0n) {
      throw new ValueError('negative exponent not allowed for Integer');
    }
    return new Integer(this.value ** nBig);
  }

  // Comparison
  eq(n: Integer | bigint): boolean {
    const other = n instanceof Integer ? n.value : n;
    return this.value === other;
  }

  lt(n: Integer | bigint): boolean {
    const other = n instanceof Integer ? n.value : n;
    return this.value < other;
  }

  le(n: Integer | bigint): boolean {
    const other = n instanceof Integer ? n.value : n;
    return this.value <= other;
  }

  gt(n: Integer | bigint): boolean {
    const other = n instanceof Integer ? n.value : n;
    return this.value > other;
  }

  ge(n: Integer | bigint): boolean {
    const other = n instanceof Integer ? n.value : n;
    return this.value >= other;
  }

  /**
   * Check if this integer is zero.
   */
  isZero(): boolean {
    return this.value === 0n;
  }

  toString(): string {
    return this.value.toString();
  }

  valueOf(): bigint {
    return this.value;
  }

  // ============================================
  // Additional SageMath Integer methods (stubs)
  // ============================================

  /**
   * Return the n-th root of this integer.
   *
   * If truncate_mode is false (default), returns the exact n-th root
   * if self is an n-th power, or throws ValueError if it is not.
   *
   * If truncate_mode is true, returns [root, exact_flag] where root is the
   * truncated n-th root (rounded towards zero) and exact_flag indicates
   * whether the root extraction was exact.
   *
   * @param n - Positive integer n >= 1
   * @param truncate_mode - Whether to return truncated root with exactness flag
   * @returns The n-th root or [root, exact] tuple
   * @throws {ValueError} If n < 1 or if taking even root of negative number
   * @see Reference: sage/rings/integer.pyx:nth_root
   */
  nth_root(n: IntegerLike, truncate_mode?: false): Integer;
  nth_root(n: IntegerLike, truncate_mode: true): [Integer, boolean];
  nth_root(n: IntegerLike, truncate_mode: boolean = false): Integer | [Integer, boolean] {
    const nBig = toBigInt(n);
    if (nBig < 1n) {
      throw new ValueError(`n (=${nBig}) must be positive`);
    }

    const isNegative = this.value < 0n;

    // Cannot take even root of negative number
    if (isNegative && (nBig & 1n) === 0n) {
      throw new ValueError('cannot take even root of negative number');
    }

    // Work with absolute value
    const absVal = isNegative ? -this.value : this.value;

    // Special cases
    if (absVal === 0n) {
      if (truncate_mode) {
        return [new Integer(0n), true];
      }
      return new Integer(0n);
    }

    if (absVal === 1n) {
      const result = isNegative ? -1n : 1n;
      if (truncate_mode) {
        return [new Integer(result), true];
      }
      return new Integer(result);
    }

    if (nBig === 1n) {
      if (truncate_mode) {
        return [new Integer(this.value), true];
      }
      return new Integer(this.value);
    }

    // Binary search for the n-th root
    // We find the largest x such that x^n <= absVal
    let low = 1n;
    let high = absVal;

    // Start with a better initial guess using bit length
    const bitLen = absVal.toString(2).length;
    const approxRootBits = BigInt(Math.ceil(bitLen / Number(nBig)));
    high = 1n << (approxRootBits + 1n);
    if (high > absVal) high = absVal;

    while (low < high) {
      const mid = (low + high + 1n) / 2n;
      const midPow = mid ** nBig;
      if (midPow <= absVal) {
        low = mid;
      } else {
        high = mid - 1n;
      }
    }

    const root = isNegative ? -low : low;
    const rootPow = low ** nBig;
    const isExact = rootPow === absVal;

    if (truncate_mode) {
      return [new Integer(root), isExact];
    }

    if (isExact) {
      return new Integer(root);
    }

    // Generate ordinal suffix for error message
    const nNum = Number(nBig);
    let suffix = 'th';
    if (nNum % 100 < 11 || nNum % 100 > 13) {
      if (nNum % 10 === 1) suffix = 'st';
      else if (nNum % 10 === 2) suffix = 'nd';
      else if (nNum % 10 === 3) suffix = 'rd';
    }

    throw new ValueError(`${this.value} is not a ${nBig}${suffix} power`);
  }

  /**
   * Return the largest integer k such that b^k <= self.
   *
   * This is the floor of log_b(self), computed exactly using integer arithmetic.
   *
   * @param b - Base (must be >= 2)
   * @returns Floor of log_b(self)
   * @throws {ValueError} If self <= 0 or b < 2
   * @see Reference: sage/rings/integer.pyx:exact_log
   */
  exact_log(b: Integer | bigint): bigint {
    const base = b instanceof Integer ? b.value : b;

    if (base < 2n) {
      throw new ValueError('base must be >= 2');
    }

    if (this.value <= 0n) {
      throw new ValueError('self must be positive');
    }

    if (this.value < base) {
      return 0n;
    }

    if (this.value === base) {
      return 1n;
    }

    // Binary search for the exact log
    // Find k such that base^k <= self < base^(k+1)

    // Get initial bounds using bit lengths
    const selfBits = BigInt(this.value.toString(2).length);
    const baseBits = BigInt(base.toString(2).length);

    // Lower bound: selfBits / (baseBits + 1) since base^k has at most k*(baseBits) bits
    let low = selfBits / (baseBits + 1n);
    // Upper bound: selfBits / (baseBits - 1) + 1
    let high = baseBits > 1n ? selfBits / (baseBits - 1n) + 1n : selfBits;

    // Ensure low >= 0
    if (low < 0n) low = 0n;

    // Binary search
    while (low < high) {
      const mid = (low + high + 1n) / 2n;
      const midPow = base ** mid;
      if (midPow <= this.value) {
        low = mid;
      } else {
        high = mid - 1n;
      }
    }

    return low;
  }

  /**
   * Return the largest divisor of this integer coprime to m.
   *
   * This is n / gcd(n, m^infinity), i.e., n with all prime factors
   * that also divide m removed.
   *
   * @param m - Integer
   * @returns The prime-to-m part of self
   * @see Reference: sage/rings/integer.pyx:prime_to_m_part
   */
  prime_to_m_part(m: Integer | bigint): Integer {
    const mVal = m instanceof Integer ? m.value : m;
    return new Integer(_prime_to_m_part(this.value, mVal));
  }

  /**
   * Return the list of prime divisors of this integer.
   *
   * @returns Sorted list of distinct prime divisors
   * @throws {ValueError} If self is 0
   * @see Reference: sage/rings/integer.pyx:prime_divisors
   */
  prime_divisors(): Integer[] {
    return _prime_factors(this.value).map((p) => new Integer(p));
  }

  /**
   * Return the list of all positive divisors of this integer.
   *
   * @returns Sorted list of positive divisors
   * @throws {ValueError} If self is 0
   * @see Reference: sage/rings/integer.pyx:divisors
   */
  divisors(): Integer[] {
    return _divisors(this.value).map((d) => new Integer(d));
  }

  /**
   * Return the Jacobi symbol (self/n).
   *
   * @param n - Positive odd integer
   * @returns -1, 0, or 1
   * @throws {ValueError} If n is not a positive odd integer
   * @see Reference: sage/rings/integer.pyx:jacobi
   */
  jacobi(n: Integer | bigint): bigint {
    const nVal = n instanceof Integer ? n.value : n;
    return _jacobi_symbol(this.value, nVal);
  }

  /**
   * Return the Kronecker symbol (self/n).
   *
   * The Kronecker symbol is an extension of the Jacobi symbol to all integers.
   *
   * @param n - Integer
   * @returns -1, 0, or 1
   * @see Reference: sage/rings/integer.pyx:kronecker
   */
  kronecker(n: Integer | bigint): bigint {
    const nVal = n instanceof Integer ? n.value : n;
    return _kronecker_symbol(this.value, nVal);
  }

  /**
   * Return the class number of the quadratic order with this discriminant.
   *
   * INPUT:
   * - self: an integer congruent to 0 or 1 mod 4 which is not a perfect square
   *
   * OUTPUT: the class number of the quadratic order with this discriminant
   *
   * NOTE: For positive D, this is the narrow class number. For negative D,
   * this equals the number of classes of primitive binary quadratic forms.
   *
   * @returns Class number
   * @throws {ValueError} If self is a perfect square or not congruent to 0 or 1 mod 4
   * @see Reference: sage/rings/integer.pyx:class_number
   * @see Deviation: Unimplemented Number-Theoretic Functions
   */
  class_number(): bigint {
    const D = this.value;

    // Check that D is not a perfect square
    if (D >= 0n) {
      const s = isqrt(D);
      if (s * s === D) {
        throw new ValueError('class_number not defined for square integers');
      }
    }

    // Check that D ≡ 0 or 1 (mod 4)
    const mod4 = ((D % 4n) + 4n) % 4n;
    if (mod4 !== 0n && mod4 !== 1n) {
      throw new ValueError('class_number only defined for integers congruent to 0 or 1 modulo 4');
    }

    // Known class numbers for small imaginary quadratic discriminants
    const knownImaginary: Record<string, bigint> = {
      '-3': 1n,
      '-4': 1n,
      '-7': 1n,
      '-8': 1n,
      '-11': 1n,
      '-12': 1n,
      '-15': 2n,
      '-16': 1n,
      '-19': 1n,
      '-20': 2n,
      '-23': 3n,
      '-24': 2n,
      '-27': 1n,
      '-28': 1n,
      '-31': 3n,
      '-35': 2n,
      '-36': 2n,
      '-39': 4n,
      '-40': 2n,
      '-43': 1n,
      '-44': 3n,
      '-47': 5n,
      '-48': 2n,
      '-51': 2n,
      '-52': 2n,
      '-55': 4n,
      '-56': 4n,
      '-59': 3n,
      '-60': 2n,
      '-63': 4n,
      '-67': 1n,
      '-68': 4n,
      '-71': 7n,
      '-72': 2n,
      '-75': 2n,
      '-76': 3n,
      '-79': 5n,
      '-80': 4n,
      '-83': 3n,
      '-84': 4n,
      '-87': 6n,
      '-88': 2n,
      '-91': 2n,
      '-92': 3n,
      '-95': 8n,
      '-96': 4n,
      '-99': 2n,
      '-100': 2n,
      '-103': 5n,
      '-104': 6n,
      '-107': 3n,
      '-108': 3n,
      '-111': 8n,
      '-112': 2n,
      '-115': 2n,
      '-116': 6n,
      '-119': 10n,
      '-120': 4n,
      '-123': 2n,
      '-124': 3n,
      '-127': 5n,
      '-128': 4n,
      '-131': 5n,
      '-132': 4n,
      '-135': 6n,
      '-136': 4n,
      '-139': 3n,
      '-140': 6n,
      '-143': 10n,
      '-144': 4n,
      '-147': 2n,
      '-148': 2n,
      '-151': 7n,
      '-152': 6n,
      '-155': 4n,
      '-156': 4n,
      '-159': 10n,
      '-160': 4n,
      '-163': 1n,
      '-164': 8n,
      '-167': 11n,
    };

    // Known class numbers for small real quadratic discriminants
    const knownReal: Record<string, bigint> = {
      '5': 1n,
      '8': 1n,
      '12': 1n,
      '13': 1n,
      '17': 1n,
      '21': 1n,
      '24': 1n,
      '28': 1n,
      '29': 1n,
      '33': 1n,
      '37': 1n,
      '40': 2n,
      '41': 1n,
      '44': 1n,
      '53': 1n,
      '56': 1n,
      '57': 1n,
      '60': 2n,
      '61': 1n,
      '65': 2n,
      '69': 1n,
      '73': 1n,
      '76': 1n,
      '77': 1n,
      '85': 2n,
      '88': 1n,
      '89': 1n,
      '92': 1n,
      '93': 1n,
      '97': 1n,
    };

    const key = D.toString();
    if (D < 0n) {
      if (key in knownImaginary) {
        return knownImaginary[key]!;
      }
    } else {
      if (key in knownReal) {
        return knownReal[key]!;
      }
    }

    // For larger discriminants, we would need PARI's qfbclassno
    throw new NotImplementedError(
      `class_number: computation for discriminant ${D} requires PARI integration`
    );
  }

  /**
   * Return the squarefree part of this integer.
   *
   * The squarefree part is the unique integer z such that n = z * y^2
   * where y^2 is a perfect square and z is squarefree.
   *
   * @returns The squarefree part of self
   * @see Reference: sage/rings/integer.pyx:squarefree_part
   */
  squarefree_part(): Integer {
    return new Integer(_squarefree_part(this.value));
  }

  /**
   * Return the next prime greater than this integer.
   * @see Reference: sage/rings/integer.pyx:next_prime
   */
  next_prime(): Integer {
    return new Integer(_next_prime(this.value));
  }

  /**
   * Return the next prime power greater than this integer.
   *
   * A prime power is a prime raised to a positive power (p, p^2, p^3, ...).
   *
   * @returns The smallest prime power > self
   * @see Reference: sage/rings/integer.pyx:next_prime_power
   */
  next_prime_power(): Integer {
    return new Integer(_next_prime_power(this.value));
  }

  /**
   * Test if this integer is a prime power (p^k for prime p and k >= 1).
   *
   * Note: 1 is not a prime power.
   *
   * @returns true if self is a prime power
   * @see Reference: sage/rings/integer.pyx:is_prime_power
   */
  is_prime_power(): boolean {
    return _is_prime_power(this.value);
  }

  /**
   * Test if this integer is a perfect power (n = a^k for some k >= 2).
   *
   * Returns true if self = a^k for some integers a and k with k >= 2.
   * Note: 0, 1, and -1 are perfect powers.
   *
   * @returns true if self is a perfect power
   * @see Reference: sage/rings/integer.pyx:is_perfect_power
   */
  is_perfect_power(): boolean {
    const n = this.value;

    // 0, 1, -1 are perfect powers
    if (n === 0n || n === 1n || n === -1n) {
      return true;
    }

    const isNegative = n < 0n;
    let absN = isNegative ? -n : n;

    // For negative numbers, we need to check if absN is a perfect odd power
    // First remove all square factors to get to a non-square
    if (isNegative) {
      // Remove all square factors
      let s = isqrt(absN);
      while (s * s === absN) {
        absN = s;
        s = isqrt(absN);
      }
      // Now absN is not a perfect square; check if it's a perfect (odd) power
    }

    // Check if absN is a perfect power
    // Try small exponents: 2, 3, 5, 7, 11, ... up to log2(absN)
    const bitLen = absN.toString(2).length;

    // Check for perfect square (only if positive)
    if (!isNegative) {
      const s = isqrt(absN);
      if (s * s === absN) {
        return true;
      }
    }

    // Check odd exponents starting from 3
    for (let k = 3; k <= bitLen; k += 2) {
      const [root, exact] = new Integer(absN).nth_root(BigInt(k), true);
      if (exact) {
        return true;
      }
    }

    // Also check even exponents if positive
    if (!isNegative) {
      for (let k = 4; k <= bitLen; k += 2) {
        const [root, exact] = new Integer(absN).nth_root(BigInt(k), true);
        if (exact) {
          return true;
        }
      }
    }

    return false;
  }

  /**
   * Test if this integer is irreducible.
   *
   * In the integers, an element is irreducible iff it is +-prime.
   * This differs from is_prime() which requires positivity.
   *
   * @returns true if self is irreducible (|self| is prime)
   * @see Reference: sage/rings/integer.pyx:is_irreducible
   */
  is_irreducible(): boolean {
    const absVal = this.value < 0n ? -this.value : this.value;
    return _is_prime(absVal);
  }

  /**
   * Test if this integer is a pseudoprime.
   *
   * This uses probabilistic primality testing (Miller-Rabin).
   * The result is NOT proven correct.
   *
   * @returns true if self is a pseudoprime
   * @see Reference: sage/rings/integer.pyx:is_pseudoprime
   */
  is_pseudoprime(): boolean {
    return _is_pseudoprime(this.value);
  }

  /**
   * Test if this integer is squarefree (not divisible by any perfect square > 1).
   *
   * @returns true if self is squarefree
   * @see Reference: sage/rings/integer.pyx:is_squarefree
   */
  is_squarefree(): boolean {
    return _is_squarefree(this.value);
  }

  /**
   * Test if this integer is a discriminant.
   *
   * A discriminant is an integer congruent to 0 or 1 modulo 4. Note that this
   * includes 0, 1 and perfect squares such as 100 (see the Sage doctests).
   *
   * @returns true if self is a discriminant
   * @see Reference: sage/rings/integer.pyx:6295 (is_discriminant)
   */
  is_discriminant(): boolean {
    // Sage is literally ``self % 4 in [0, 1]`` (Python's non-negative mod).
    const mod4 = ((this.value % 4n) + 4n) % 4n;
    return mod4 === 0n || mod4 === 1n;
  }

  /**
   * Test if this integer is a fundamental discriminant.
   *
   * A fundamental discriminant is a discriminant, not 0 or 1, and not a square
   * multiple of a smaller discriminant.
   *
   * @returns true if self is a fundamental discriminant
   * @see Reference: sage/rings/integer.pyx:6325 (is_fundamental_discriminant)
   */
  is_fundamental_discriminant(): boolean {
    const D = this.value;

    if (D === 0n || D === 1n) {
      return false;
    }

    const mod4 = ((D % 4n) + 4n) % 4n;

    if (mod4 === 2n || mod4 === 3n) {
      return false;
    }

    if (mod4 === 1n) {
      // D ≡ 1 (mod 4): must be squarefree
      return _is_squarefree(D);
    }

    // D ≡ 0 (mod 4): d = D // 4 (floor division) must satisfy
    // d % 4 in [2, 3] and d squarefree.
    const d = fdiv_q(D, 4n);
    const dmod4 = ((d % 4n) + 4n) % 4n;
    return (dmod4 === 2n || dmod4 === 3n) && _is_squarefree(d);
  }

  /**
   * Return the binomial coefficient C(n, k) where n = self.
   *
   * C(n, k) = n! / (k! * (n-k)!)
   *
   * @param k - Non-negative integer
   * @returns The binomial coefficient
   * @see Reference: sage/rings/integer.pyx:binomial
   */
  binomial(k: Integer | bigint): Integer {
    const kVal = k instanceof Integer ? k.value : k;
    return new Integer(_binomial(this.value, kVal));
  }

  /**
   * Return the factorial n! where n = self.
   *
   * @returns n!
   * @throws {ValueError} If self is negative
   * @see Reference: sage/rings/integer.pyx:factorial
   */
  factorial(): Integer {
    return new Integer(_factorial(this.value));
  }

  /**
   * Return Euler's totient function phi(n).
   *
   * phi(n) is the number of integers k with 1 <= k <= n and gcd(k, n) = 1.
   *
   * @returns phi(self)
   * @see Reference: sage/rings/integer.pyx:euler_phi
   */
  euler_phi(): Integer {
    return new Integer(_euler_phi(this.value));
  }

  /**
   * Return the sum of divisors function sigma_k(n).
   *
   * sigma_k(n) = sum of d^k for all divisors d of n.
   * When k=1 (default), this is the sum of divisors.
   * When k=0, this is the number of divisors.
   *
   * @param k - Power (default: 1)
   * @returns sigma_k(self)
   * @see Reference: sage/rings/integer.pyx:sigma
   */
  sigma(k: IntegerLike = 1n): Integer {
    const kBig = toBigInt(k);
    return new Integer(_sigma(this.value, kBig));
  }

  /**
   * Return the Moebius mu function.
   *
   * mu(n) = 0 if n has a squared prime factor
   * mu(n) = (-1)^k if n is a product of k distinct primes
   * mu(1) = 1
   *
   * @returns -1, 0, or 1
   * @see Reference: sage/rings/integer.pyx:moebius
   */
  moebius(): bigint {
    return _moebius(this.value);
  }

  /**
   * Return the radical of this integer (product of distinct prime factors).
   *
   * rad(n) = product of p for all primes p dividing n.
   *
   * @returns The radical of self
   * @see Reference: sage/rings/integer.pyx:radical
   */
  radical(): Integer {
    return new Integer(_radical(this.value));
  }

  /**
   * Return the number of divisors of this integer.
   *
   * @returns Number of positive divisors
   * @see Reference: sage/rings/integer.pyx:number_of_divisors
   */
  number_of_divisors(): bigint {
    const absVal = this.value < 0n ? -this.value : this.value;
    return _number_of_divisors(absVal);
  }

  /**
   * Return the digit sum in base b.
   *
   * This is the sum of the absolute values of the digits in base b.
   *
   * @param b - Base (default: 10, must be >= 2)
   * @returns Sum of digits
   * @see Reference: sage/rings/integer.pyx:digit_sum
   */
  digit_sum(b: IntegerLike = 10n): bigint {
    const bBig = toBigInt(b);
    if (bBig < 2n) {
      throw new ValueError('base must be >= 2');
    }

    if (this.value === 0n) {
      return 0n;
    }

    let n = this.value < 0n ? -this.value : this.value;
    let sum = 0n;

    while (n > 0n) {
      sum += n % bBig;
      n = n / bBig;
    }

    return sum;
  }

  /**
   * Return the digits in base b in little-endian order.
   *
   * For negative numbers, the digits are negated (matching SageMath behavior).
   *
   * @param b - Base (default: 10, must be >= 2)
   * @returns Array of digits in little-endian order
   * @see Reference: sage/rings/integer.pyx:digits
   */
  digits(b: IntegerLike = 10n): bigint[] {
    const bBig = toBigInt(b);
    if (bBig < 2n) {
      throw new ValueError('base must be >= 2');
    }

    if (this.value === 0n) {
      return [];
    }

    const isNegative = this.value < 0n;
    let n = isNegative ? -this.value : this.value;
    const result: bigint[] = [];

    while (n > 0n) {
      const digit = n % bBig;
      result.push(isNegative ? -digit : digit);
      n = n / bBig;
    }

    return result;
  }

  /**
   * Return the popcount (number of 1 bits in binary representation).
   *
   * For negative numbers, SageMath returns Infinity. Since we can't represent
   * Infinity as bigint, we throw an error for negative inputs.
   *
   * @returns Number of 1 bits
   * @throws {ValueError} If self is negative
   * @see Reference: sage/rings/integer.pyx:popcount
   */
  popcount(): bigint {
    if (this.value < 0n) {
      throw new ValueError('popcount of negative integer is infinite');
    }

    let n = this.value;
    let count = 0n;

    while (n > 0n) {
      count += n & 1n;
      n >>= 1n;
    }

    return count;
  }

  /**
   * Return the Hamming weight (same as popcount for non-negative integers).
   *
   * For non-negative integers, this is the number of 1 bits.
   * For negative integers, SageMath returns Infinity; we throw an error.
   *
   * @returns Hamming weight
   * @throws {ValueError} If self is negative
   * @see Reference: sage/rings/integer.pyx:hamming_weight
   */
  hamming_weight(): bigint {
    return this.popcount();
  }

  /**
   * Return the square of this integer.
   * @see Reference: sage/rings/integer.pyx:square
   */
  square(): Integer {
    return new Integer(this.value * this.value);
  }

  /**
   * Return the cube of this integer.
   * @see Reference: sage/rings/integer.pyx:cube
   */
  cube(): Integer {
    return new Integer(this.value * this.value * this.value);
  }

  /**
   * Return self modulo m.
   *
   * As in Python/GMP (``mpz_fdiv_r``), the result is zero or has the same sign
   * as ``m``: ``5 % -7 == -2`` and ``(-5) % -7 == -5``.
   *
   * @see Reference: sage/rings/integer.pyx:__mod__
   */
  __mod__(m: Integer | bigint): Integer {
    const other = m instanceof Integer ? m.value : m;
    if (other === 0n) {
      throw new ZeroDivisionError('Integer modulo by zero');
    }
    return new Integer(fdiv_r(this.value, other));
  }

  /**
   * Return the floor division of self by other (``mpz_fdiv_q``).
   * @see Reference: sage/rings/integer.pyx:_floordiv_
   */
  __floordiv__(other: Integer | bigint): Integer {
    const o = other instanceof Integer ? other.value : other;
    if (o === 0n) {
      throw new ZeroDivisionError('Integer division by zero');
    }
    return new Integer(fdiv_q(this.value, o));
  }

  /**
   * Test divisibility: return True if m divides self.
   * @see Reference: sage/rings/integer.pyx:divides
   */
  divides(m: Integer | bigint): boolean {
    const other = m instanceof Integer ? m.value : m;
    if (this.value === 0n) {
      return other === 0n;
    }
    return other % this.value === 0n;
  }

  /**
   * Return the content of this integer (absolute value).
   * @see Reference: sage/rings/integer.pyx:content
   */
  content(): Integer {
    return this.abs();
  }

  /**
   * Return the primitive part of this integer (sign).
   * @see Reference: sage/rings/integer.pyx:primitive_part
   */
  primitive_part(): Integer {
    if (this.value === 0n) return new Integer(0n);
    return new Integer(this.value < 0n ? -1n : 1n);
  }

  /**
   * Return a square root of self modulo n if it exists.
   *
   * @param n - The modulus (must be a prime)
   * @returns A square root of self mod n, or null if none exists
   * @see Reference: sage/rings/integer.pyx:sqrt_mod
   */
  sqrt_mod(n: Integer | bigint): Integer | null {
    const nVal = n instanceof Integer ? n.value : n;
    const result = _sqrt_mod(this.value, nVal);
    return result !== null ? new Integer(result) : null;
  }

  /**
   * Return an n-th root of self modulo p if it exists.
   *
   * Uses the Adleman-Manders-Miller / Johnston algorithm for computing
   * n-th roots in finite fields.
   *
   * @param n - The root index (must be positive)
   * @param p - The modulus (must be a prime)
   * @returns An n-th root of self mod p
   * @throws {ValueError} If no n-th root exists or p is not prime
   * @see Reference: sage/rings/finite_rings/integer_mod.pyx:nth_root
   * @see Reference: sage/rings/finite_rings/element_base.pyx:_nth_root_common
   */
  nth_root_mod(n: Integer | bigint, p: Integer | bigint): Integer {
    const nVal = n instanceof Integer ? n.value : n;
    const pVal = p instanceof Integer ? p.value : p;

    if (nVal <= 0n) {
      throw new ValueError('n must be positive');
    }

    if (pVal < 2n || !_is_prime(pVal)) {
      throw new ValueError('p must be a prime');
    }

    // Normalize a to [0, p)
    let a = ((this.value % pVal) + pVal) % pVal;

    // Special case: a = 0
    if (a === 0n) {
      return new Integer(0n);
    }

    // Special case: n = 1
    if (nVal === 1n) {
      return new Integer(a);
    }

    // For n=2, use optimized sqrt_mod
    if (nVal === 2n) {
      const result = _sqrt_mod(a, pVal);
      if (result === null) {
        throw new ValueError('no n-th root');
      }
      return new Integer(result);
    }

    // General case: transcription of Sage's
    // FiniteRingElement._nth_root_common (Johnston / Adleman-Manders-Miller),
    // reference/sage/src/sage/rings/finite_rings/element_base.pyx:29-114.
    const q = pVal - 1n; // order of the multiplicative group, i.e. Sage's q-1
    const gcd0 = _gcd(nVal, q);

    if (a === 1n) {
      // Sage returns ``K.zeta(gcd)``; 1 is always an n-th root of 1 and this
      // method returns a single root, so keep the canonical one.
      return new Integer(1n);
    }

    if (gcd0 === q) {
      throw new ValueError('no n-th root');
    }

    // gcd = alpha*n + beta*(q-1), so 1/n = alpha/gcd (mod q-1)
    const [, alpha] = _xgcd(nVal, q);

    if (gcd0 === 1n) {
      return new Integer(_power_mod(a, alpha, pVal));
    }

    const nRed = gcd0;
    const q1overn = q / nRed;
    if (_power_mod(a, q1overn, pVal) !== 1n) {
      throw new ValueError('no n-th root');
    }

    a = _power_mod(a, alpha, pVal);

    for (const [r, v] of _factor(nRed)) {
      if (r === -1n) continue;

      // (q-1).val_unit(r): q-1 = r^k * h with gcd(r, h) = 1.  0 < v <= k.
      let k = 0n;
      let h = q;
      while (h % r === 0n) {
        h /= r;
        k++;
      }

      const rv = r ** v;
      // hinv = (-h)^(-1) mod r^v
      const hinv = _inverse_mod(((-h % rv) + rv) % rv, rv);
      const z = h * hinv;
      const x = (1n + z) / rv;

      if (k === v) {
        a = _power_mod(a, x, pVal);
      } else {
        // We need an element of order exactly r^k (Sage's ``K.zeta(r**k)``;
        // ``g^h`` in Johnston's article).  Checking ``c^h != 1`` alone is not
        // enough -- that only bounds the order below by r.
        const gh = zetaPrimePower(pVal, r, k, h);
        const t = discreteLog(_power_mod(a, h, pVal), _power_mod(gh, rv, pVal), r ** (k - v), pVal);
        a = (_power_mod(a, x, pVal) * _power_mod(gh, -hinv * t, pVal)) % pVal;
      }
    }

    return new Integer(a);
  }

  /**
   * Return the multiplicative order of self **in the ring of integers**.
   *
   * As in Sage this is 1 for 1, 2 for -1, and ``+Infinity`` for every other
   * integer (no other integer is a unit of infinite order in ZZ).
   *
   * For the order of a residue class modulo `n`, use
   * ``Mod(a, n).multiplicative_order()`` from
   * ``rings/finite_rings/integer_mod`` -- exactly as in Sage.
   *
   * @returns 1n, 2n or the string 'Infinity'
   * @see Reference: sage/rings/integer.pyx:6257 (multiplicative_order)
   */
  multiplicative_order(): bigint | 'Infinity' {
    if (this.value === 1n) {
      return 1n;
    }
    if (this.value === -1n) {
      return 2n;
    }
    return 'Infinity';
  }

  /**
   * Test if self is a primitive root modulo n.
   *
   * self is a primitive root mod n iff the multiplicative order of self
   * equals phi(n).
   *
   * @param n - The modulus
   * @returns true if self is a primitive root mod n
   * @see Reference: sage/rings/integer.pyx:is_primitive_root
   */
  is_primitive_root(n: Integer | bigint): boolean {
    const nVal = n instanceof Integer ? n.value : n;

    if (nVal <= 1n) {
      return false;
    }

    const a = ((this.value % nVal) + nVal) % nVal;

    // Check gcd(a, n) = 1
    if (_gcd(a, nVal) !== 1n) {
      return false;
    }

    // Check if order equals phi(n)
    const phi = _euler_phi(nVal);
    const order = multiplicativeOrderMod(a, nVal);

    return order === phi;
  }

  /**
   * Return the inverse of self modulo n.
   *
   * @param n - The modulus
   * @returns The inverse of self mod n
   * @throws {ZeroDivisionError} If gcd(self, n) != 1
   * @see Reference: sage/rings/integer.pyx:inverse_mod
   */
  inverse_mod(n: Integer | bigint): Integer {
    const nVal = n instanceof Integer ? n.value : n;
    return new Integer(_inverse_mod(this.value, nVal));
  }

  /**
   * Return self raised to power e modulo n.
   *
   * @param e - The exponent
   * @param n - The modulus
   * @returns self^e mod n
   * @see Reference: sage/rings/integer.pyx:powermod
   */
  powermod(e: Integer | bigint, n: Integer | bigint): Integer {
    const eVal = e instanceof Integer ? e.value : e;
    const nVal = n instanceof Integer ? n.value : n;
    return new Integer(_power_mod(this.value, eVal, nVal));
  }

  /**
   * Return the integer part of the log base b of self.
   *
   * This is the same as exact_log().
   *
   * @param b - Base (default: e, which returns floor(ln(self)))
   * @returns Floor of log_b(self)
   * @see Reference: sage/rings/integer.pyx:log
   * @see Deviation: `Integer.log` always floors. Sage returns an *exact*
   *   Integer only when `b^k == self`, and otherwise a real/symbolic
   *   logarithm (e.g. `Integer(8).log(2) == 3` but `Integer(9).log(2)` is
   *   `log(9)/log(2)`). We have no symbolic ring, so we return
   *   `exact_log(b) = floor(log_b(self))` unconditionally.
   */
  log(b?: Integer | bigint): Integer {
    if (b === undefined) {
      // Natural log - use floating point approximation
      if (this.value <= 0n) {
        throw new ValueError('log of non-positive number');
      }
      const ln = Math.log(Number(this.value));
      return new Integer(BigInt(Math.floor(ln)));
    }

    return new Integer(this.exact_log(b));
  }

  /**
   * Return the real-valued natural log of this integer.
   *
   * @returns ln(self) as a floating point number
   * @throws {ValueError} If self <= 0
   * @see Reference: sage/rings/integer.pyx:real_log
   */
  real_log(): number {
    if (this.value <= 0n) {
      throw new ValueError('log of non-positive number');
    }

    // For values that do not fit exactly in a double, split off the top 53
    // bits: ln(n) = ln(n >> s) + s*ln(2) with s = bitlen(n) - 53.
    const bits = BigInt(this.value.toString(2).length);
    if (bits > 53n) {
      const shift = bits - 53n;
      const mantissa = Number(this.value >> shift);
      return Math.log(mantissa) + Number(shift) * Math.LN2;
    }

    return Math.log(Number(this.value));
  }

  /**
   * Return the continued fraction expansion of self.
   *
   * For an integer, this is just [self].
   *
   * @returns Array containing self
   * @see Reference: sage/rings/integer.pyx:continued_fraction
   */
  continued_fraction(): bigint[] {
    return [this.value];
  }

  /**
   * Return the p-adic valuation of self.
   * @see Reference: sage/rings/integer.pyx:ord
   */
  ord(p: Integer | bigint): bigint {
    return this.valuation(p);
  }

  /**
   * Return the integer square root and remainder.
   *
   * Returns [s, r] where self = s^2 + r and 0 <= r <= 2*s.
   *
   * @returns [sqrt, remainder]
   * @throws {ValueError} If self is negative
   * @see Reference: sage/rings/integer.pyx:sqrtrem
   */
  sqrtrem(): [Integer, Integer] {
    if (this.value < 0n) {
      throw new ValueError('sqrtrem requires non-negative input');
    }

    const s = isqrt(this.value);
    const r = this.value - s * s;
    return [new Integer(s), new Integer(r)];
  }

  /**
   * Return whether self is a quadratic residue modulo p.
   *
   * @param p - A prime modulus
   * @returns true if self is a quadratic residue mod p
   * @see Reference: sage/rings/integer.pyx:is_quadratic_residue
   */
  is_quadratic_residue(p: Integer | bigint): boolean {
    const pVal = p instanceof Integer ? p.value : p;
    const ls = _legendre_symbol(this.value, pVal);
    return ls !== -1n;
  }

  /**
   * Return Legendre symbol (self/p).
   *
   * @param p - An odd prime
   * @returns -1, 0, or 1
   * @see Reference: sage/rings/integer.pyx:legendre_symbol
   */
  legendre_symbol(p: Integer | bigint): bigint {
    const pVal = p instanceof Integer ? p.value : p;
    return _legendre_symbol(this.value, pVal);
  }

  /**
   * Return the bit at position n (0-indexed from least significant bit).
   *
   * @param n - Bit position (must be non-negative)
   * @returns 0n or 1n
   * @see Reference: sage/rings/integer.pyx:bit
   */
  bit(n: bigint): bigint {
    if (n < 0n) {
      return 0n;
    }
    return (this.value >> n) & 1n;
  }

  /**
   * Return the bits as a list (little-endian order).
   *
   * This is equivalent to digits(2).
   *
   * @returns Array of bits (0n or 1n, or -1n for negative numbers)
   * @see Reference: sage/rings/integer.pyx:bits
   */
  bits(): bigint[] {
    return this.digits(2n);
  }

  /**
   * Return the bit length (same as nbits).
   * @see Reference: sage/rings/integer.pyx:bit_length
   */
  bit_length(): bigint {
    return this.nbits();
  }

  /**
   * Return the numerator (self for integers).
   * @see Reference: sage/rings/integer.pyx:numerator
   */
  numerator(): Integer {
    return this;
  }

  /**
   * Return the denominator (1 for integers).
   * @see Reference: sage/rings/integer.pyx:denominator
   */
  denominator(): Integer {
    return new Integer(1n);
  }

  /**
   * Return the floor (self for integers).
   * @see Reference: sage/rings/integer.pyx:floor
   */
  floor(): Integer {
    return this;
  }

  /**
   * Return the ceiling (self for integers).
   * @see Reference: sage/rings/integer.pyx:ceil
   */
  ceil(): Integer {
    return this;
  }

  /**
   * Round to nearest integer (self for integers).
   * @see Reference: sage/rings/integer.pyx:round
   */
  round(): Integer {
    return this;
  }

  /**
   * Return self truncated toward zero (self for integers).
   * @see Reference: sage/rings/integer.pyx:trunc
   */
  trunc(): Integer {
    return this;
  }

  /**
   * Return the fractional part (0 for integers).
   * @see Reference: sage/rings/integer.pyx:frac
   */
  frac(): Integer {
    return new Integer(0n);
  }

  /**
   * Test if self is coprime to other.
   * @see Reference: sage/rings/integer.pyx:is_coprime
   */
  is_coprime(other: Integer | bigint): boolean {
    return this.gcd(other).value === 1n;
  }

  /**
   * Return the smallest prime factor of self via trial division.
   *
   * If bound is given, only check primes up to that bound.
   * If no factor is found, returns self.
   *
   * @param bound - Maximum prime to check
   * @returns Smallest prime factor (or self if none found)
   * @see Reference: sage/rings/integer.pyx:trial_division
   */
  trial_division(bound?: bigint): Integer {
    return new Integer(_trial_division(this.value, bound));
  }

  /**
   * Return whether self is a power of 2.
   * @see Reference: sage/rings/integer.pyx:is_power_of_two
   */
  is_power_of_two(): boolean {
    if (this.value <= 0n) return false;
    return (this.value & (this.value - 1n)) === 0n;
  }

  /**
   * Return whether self is even.
   * @see Reference: sage/rings/integer.pyx:is_even
   */
  is_even(): boolean {
    return (this.value & 1n) === 0n;
  }

  /**
   * Return whether self is odd.
   * @see Reference: sage/rings/integer.pyx:is_odd
   */
  is_odd(): boolean {
    return (this.value & 1n) === 1n;
  }

  /**
   * Return the multiplicative inverse of this integer.
   *
   * In the ring of integers ZZ, only the units +1 and -1 have multiplicative
   * inverses (themselves). For other integers, inversion would require
   * working in the rationals QQ.
   *
   * NOTE: In SageMath, Integer.__invert__ returns a Rational 1/self.
   * Since we don't have Rational type here, we return Integer for units
   * and throw for non-units.
   *
   * @returns The inverse (only for +/-1)
   * @throws {ArithmeticError} If self is not a unit
   * @see Reference: sage/rings/integer.pyx:__invert__
   */
  __invert__(): Integer {
    if (this.value === 1n || this.value === -1n) {
      return this;
    }
    // In SageMath, this would return Rational(1, self.value)
    // Since we're in ZZ and don't have Rational here, we throw
    throw new ArithmeticError(
      `Integer ${this.value} is not invertible in ZZ (only +/-1 are units)`
    );
  }

  /**
   * Return the n-th Bell number if self = n.
   *
   * Bell numbers count the number of partitions of a set.
   * B_0 = 1, B_1 = 1, B_2 = 2, B_3 = 5, B_4 = 15, ...
   *
   * @returns The n-th Bell number
   * @throws {ValueError} If self is negative
   * @see Reference: sage/rings/integer.pyx:bell_number
   */
  bell_number(): Integer {
    const n = this.value;

    if (n < 0n) {
      throw new ValueError('Bell number only defined for non-negative integers');
    }

    if (n === 0n || n === 1n) {
      return new Integer(1n);
    }

    // Use Bell triangle (similar to Pascal's triangle)
    // B_n is the first element of the n-th row
    let row: bigint[] = [1n];

    for (let i = 1n; i <= n; i++) {
      const newRow: bigint[] = [row[row.length - 1]!];
      for (let j = 0; j < row.length; j++) {
        newRow.push(newRow[j]! + row[j]!);
      }
      row = newRow;
    }

    return new Integer(row[0]!);
  }

  /**
   * Return the n-th Catalan number if self = n.
   *
   * C_n = (2n)! / ((n+1)! * n!) = binomial(2n, n) / (n+1)
   *
   * @returns The n-th Catalan number
   * @throws {ValueError} If self is negative
   * @see Reference: sage/rings/integer.pyx:catalan_number
   */
  catalan_number(): Integer {
    const n = this.value;

    if (n < 0n) {
      throw new ValueError('Catalan number only defined for non-negative integers');
    }

    // C_n = binomial(2n, n) / (n + 1)
    const binom = _binomial(2n * n, n);
    return new Integer(binom / (n + 1n));
  }

  /**
   * Return the n-th Fibonacci number if self = n.
   *
   * F_0 = 0, F_1 = 1, F_n = F_{n-1} + F_{n-2}
   *
   * @returns The n-th Fibonacci number
   * @see Reference: sage/rings/integer.pyx:fibonacci
   */
  fibonacci(): Integer {
    return new Integer(_fibonacci(this.value));
  }

  /**
   * Return the n-th Lucas number if self = n.
   *
   * L_0 = 2, L_1 = 1, L_n = L_{n-1} + L_{n-2}
   *
   * @returns The n-th Lucas number
   * @see Reference: sage/rings/integer.pyx:lucas_number
   */
  lucas_number(): Integer {
    return new Integer(_lucas_number(this.value));
  }

  /**
   * Return the partition number p(n) if self = n.
   *
   * p(n) is the number of ways to write n as a sum of positive integers.
   *
   * @returns The partition number
   * @throws {ValueError} If self is negative
   * @see Reference: sage/rings/integer.pyx:number_of_partitions
   * @see Deviation: Combinatorial Function Limits
   */
  number_of_partitions(): Integer {
    const n = this.value;

    if (n < 0n) {
      return new Integer(0n);
    }

    if (n === 0n) {
      return new Integer(1n);
    }

    // Use dynamic programming for small n
    const nNum = Number(n);
    if (nNum > 10000) {
      throw new NotImplementedError('number_of_partitions: value too large');
    }

    // p[k] = number of partitions of k
    const p: bigint[] = new Array(nNum + 1).fill(0n);
    p[0] = 1n;

    for (let k = 1; k <= nNum; k++) {
      for (let i = k; i <= nNum; i++) {
        p[i] += p[i - k]!;
      }
    }

    return new Integer(p[nNum]!);
  }

  /**
   * Return a primitive root modulo self.
   *
   * A primitive root exists iff self = 1, 2, 4, p^k, or 2*p^k for odd prime p.
   *
   * @returns A primitive root mod self
   * @throws {ValueError} If no primitive root exists
   * @see Reference: sage/rings/integer.pyx:primitive_root
   */
  primitive_root(): Integer {
    return new Integer(_primitive_root(this.value));
  }

  /**
   * Return the previous prime less than self.
   *
   * @throws {ValueError} If self <= 2
   * @see Reference: sage/rings/integer.pyx:previous_prime
   */
  previous_prime(): Integer {
    return new Integer(_previous_prime(this.value));
  }

  /**
   * Return the n-th prime if self = n.
   *
   * The primes are 1-indexed: the 1st prime is 2.
   *
   * @returns The n-th prime
   * @throws {ValueError} If self <= 0
   * @see Reference: sage/rings/integer.pyx:nth_prime
   */
  nth_prime(): Integer {
    return new Integer(_nth_prime(this.value));
  }

  /**
   * Return the prime counting function pi(n) if self = n.
   *
   * pi(n) = number of primes <= n.
   *
   * @returns pi(self)
   * @see Reference: sage/rings/integer.pyx:prime_pi
   * @see Deviation: Combinatorial Function Limits
   */
  prime_pi(): Integer {
    const n = this.value;

    if (n < 2n) {
      return new Integer(0n);
    }

    // Simple counting for reasonable values
    // For very large values, more sophisticated algorithms exist (e.g., Meissel-Lehmer)
    const nNum = Number(n);
    if (nNum > 10000000) {
      throw new NotImplementedError('prime_pi: value too large for naive counting');
    }

    let count = 0n;
    for (let i = 2n; i <= n; i++) {
      if (_is_prime(i)) {
        count++;
      }
    }

    return new Integer(count);
  }

  /**
   * Check if self passes Miller-Rabin primality test with given base.
   *
   * @param base - The base for Miller-Rabin test
   * @returns true if self is a strong pseudoprime to the given base
   * @see Reference: sage/rings/integer.pyx:is_strong_pseudoprime
   */
  is_strong_pseudoprime(base: Integer | bigint): boolean {
    const baseVal = base instanceof Integer ? base.value : base;
    return _is_strong_probable_prime(this.value, baseVal);
  }

  /**
   * Return the t-th core of this integer.
   *
   * The t-th core of n is n / (largest t-th power dividing n).
   * For t=2 (default), this is the squarefree part.
   *
   * @param t - Power (default: 2)
   * @returns The t-th core of self
   * @see Reference: sage/rings/integer.pyx:core
   */
  core(t: IntegerLike = 2n): Integer {
    const tBig = toBigInt(t);
    if (tBig < 1n) {
      throw new ValueError('t must be positive');
    }

    const n = this.value;
    if (n === 0n) {
      return new Integer(0n);
    }

    const sign = n < 0n ? -1n : 1n;
    const absN = n < 0n ? -n : n;

    if (absN === 1n) {
      return new Integer(sign);
    }

    const factors = _factor(absN);
    let result = sign;

    for (const [p, e] of factors) {
      if (p === -1n) continue;
      // Include p^(e mod t) in the core
      const remainder = e % tBig;
      result *= p ** remainder;
    }

    return new Integer(result);
  }

  /**
   * Return Carmichael's lambda function.
   *
   * lambda(n) is the smallest positive integer k such that a^k ≡ 1 (mod n)
   * for all a coprime to n.
   *
   * @returns lambda(self)
   * @see Reference: sage/rings/integer.pyx:carmichael_lambda
   */
  carmichael_lambda(): Integer {
    return new Integer(_carmichael_lambda(this.value));
  }

  /**
   * Return the global height (logarithmic).
   *
   * For an integer n, the global height is log(max(1, |n|)).
   *
   * @returns log(max(1, |self|))
   * @see Reference: sage/rings/integer.pyx:global_height
   */
  global_height(): number {
    const absVal = this.value < 0n ? -this.value : this.value;
    if (absVal <= 1n) {
      return 0;
    }
    return new Integer(absVal).real_log();
  }
}

/**
 * Return the multiplicative order of ``a`` modulo ``n``.
 *
 * Factorization-based (never an O(order) loop): start from phi(n) and divide
 * out prime factors as long as the reduced exponent still gives 1.
 *
 * @internal
 */
function multiplicativeOrderMod(a: bigint, n: bigint): bigint {
  if (n === 1n) {
    return 1n;
  }
  const x = ((a % n) + n) % n;
  if (_gcd(x, n) !== 1n) {
    throw new ValueError(`${a} is not a unit modulo ${n}`);
  }
  if (x === 1n) {
    return 1n;
  }

  let order = _euler_phi(n);
  for (const [p, e] of _factor(order)) {
    if (p === -1n) continue;
    for (let i = 0n; i < e; i++) {
      if (_power_mod(x, order / p, n) !== 1n) break;
      order /= p;
    }
  }
  return order;
}

/**
 * Return an element of exact multiplicative order ``r^k`` in ``(Z/pZ)*``.
 *
 * This is the analogue of Sage's ``K.zeta(r**k)`` used by
 * ``FiniteRingElement._nth_root_common``. Here ``p - 1 = r^k * h`` with
 * ``gcd(r, h) = 1``; a candidate ``c^h`` has order dividing ``r^k`` and has
 * order exactly ``r^k`` iff ``(c^h)^(r^(k-1)) != 1``.
 *
 * @internal
 */
function zetaPrimePower(p: bigint, r: bigint, k: bigint, h: bigint): bigint {
  const rk1 = r ** (k - 1n);
  for (let c = 2n; c < p; c++) {
    const e = _power_mod(c, h, p);
    if (e !== 1n && _power_mod(e, rk1, p) !== 1n) {
      return e;
    }
  }
  throw new ValueError(`no element of order ${r ** k} modulo ${p}`);
}

/**
 * Compute discrete logarithm of a in base g of known order using baby-step giant-step.
 *
 * Finds x such that g^x = a (mod p).
 *
 * @param a - Target element
 * @param g - Generator
 * @param order - Known order of g
 * @param p - Prime modulus
 * @returns x such that g^x = a (mod p)
 * @throws {ValueError} If no discrete log exists
 * @internal
 */
function discreteLog(a: bigint, g: bigint, order: bigint, p: bigint): bigint {
  // Normalize
  a = ((a % p) + p) % p;
  g = ((g % p) + p) % p;

  if (a === 1n) {
    return 0n;
  }

  if (g === 1n) {
    if (a === 1n) return 0n;
    throw new ValueError('no discrete log');
  }

  // Baby-step giant-step algorithm
  const m = isqrt(order) + 1n;

  // Baby step: compute g^j for j = 0, 1, ..., m-1
  const table = new Map<string, bigint>();
  let gj = 1n;
  for (let j = 0n; j < m; j++) {
    table.set(gj.toString(), j);
    gj = (gj * g) % p;
  }

  // Giant step: compute a * (g^(-m))^i and look for match
  // g^(-m) = g^(order - m) mod p if order divides p-1
  const gInvM = _power_mod(g, order - (m % order), p);
  let gamma = a;
  for (let i = 0n; i < m; i++) {
    const key = gamma.toString();
    if (table.has(key)) {
      const j = table.get(key)!;
      const x = (i * m + j) % order;
      // Verify
      if (_power_mod(g, x, p) === a) {
        return x;
      }
    }
    gamma = (gamma * gInvM) % p;
  }

  throw new ValueError('no discrete log');
}
