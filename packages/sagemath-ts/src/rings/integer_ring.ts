import { mpz_remove } from '../types/gmp.js';
import { mpz_fac_ui } from '../types/gmp_factorial.js';
/**
 * @module sage/rings/integer_ring
 * @description The ring ZZ of integers
 *
 * Port of: sage/rings/integer_ring.pyx
 * Reference: reference/sage/src/sage/rings/integer_ring.pyx
 */

import { _fraction_native_integer } from './fraction_field_element.js';
import { factor_trial_division } from './factorint.js';
import { quadclassno, precprime } from '@sagemath-ts/parigp-ts';
import {
  type Factorization,
  binomial as _binomial,
  carmichael_lambda as _carmichael_lambda,
  divisors as _divisors,
  euler_phi as _euler_phi,
  factor as _factor,
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
  OverflowError,
  TypeError as SageTypeError,
  ValueError,
  ZeroDivisionError,
} from '../errors.js';
import { SAGE_RAND_MAX, current_randstate } from '../misc/randstate.js';
import { DiscreteGaussianDistributionIntegerSampler } from '../stats/distributions/discrete_gaussian_integer.js';
import { type IntegerLike, toBigInt } from '../types/coercion.js';
import { Rational } from './rational.js';

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
 * GMP's precision-doubling Newton root extraction using native BigInt limbs.
 *
 * Reference: GMP 6.3.0 mpn/generic/rootrem.c, mpn_rootrem_internal.
 * https://ftp.gnu.org/gnu/gmp/gmp-6.3.0.tar.xz
 * The size schedule, prefix/remainder update, derivative quotient, clipping
 * and downward correction follow that routine. An exact one-bit seed replaces
 * logbased_root's approximate nine-bit lookup seed; no floating point is used.
 *
 * @internal
 * @see Deviation: Integer Root Backend
 */
function mpz_root_positive(value: bigint, k: bigint, approximate = false): [bigint, boolean] {
  if (value <= 1n || k === 1n) return [value, true];
  if (k === 2n) {
    const root = isqrt(value);
    return [root, root * root === value];
  }
  const log = BigInt(value.toString(2).length - 1);
  // GMP mpn_rootrem's remp=NULL optimization: pad with k zero limbs,
  // compute one extra root limb, and usually avoid the final full-size power.
  if (!approximate && (log + 64n) / 64n > 3n * k) {
    const [root, exact] = mpz_root_positive(value << (64n * k), k, true);
    return [root >> 64n, exact];
  }
  // GMP checks this before constructing any power. It is essential when the
  // exponent dwarfs the input: the result is 1, regardless of k's magnitude.
  if (log < k) return [1n, false];
  const rootBits = log / k;
  const logk = BigInt((k - 1n).toString(2).length + 1);
  const sizes = [rootBits];
  while (sizes[sizes.length - 1]! > 0n) {
    const last = sizes[sizes.length - 1]!;
    sizes.push(last > logk ? (last + logk) / 2n : last - 1n);
  }
  let root = 1n;
  let truncated = k * rootBits;
  let bits = 0n;
  for (let i = sizes.length - 1; i > 0; i--) {
    truncated -= (k - 1n) * bits;
    const prefix = value >> truncated;
    let power = root ** (k - 1n);
    let fullPower = power * root;
    // Correct the previous approximation before determining the next bits.
    while (fullPower > prefix) {
      root--;
      power = root ** (k - 1n);
      fullPower = power * root;
    }
    bits = sizes[i - 1]! - sizes[i]!;
    truncated -= bits;
    const mask = (1n << bits) - 1n;
    const remainder = ((prefix - fullPower) << bits) | ((value >> truncated) & mask);
    const delta = remainder / (k * power);
    root = (root << bits) + (delta > mask ? mask : delta);
  }
  // The approximation is the floor root or one too large. Unless its low
  // limb is 0 or 1, truncating that limb proves a non-perfect-power result.
  if (approximate && (root & ((1n << 64n) - 1n)) > 1n) return [root, false];
  let power = root ** k;
  while (power > value) {
    root--;
    power = root ** k;
  }
  return [root, power === value];
}

/** Inputs handled by the port's Integer constructor and ZZ coercion. */
type IntegerInput =
  | IntegerLike
  | Rational
  | number
  | string
  | boolean
  | null
  | readonly unknown[]
  | { _integer_: (parent: IntegerRing) => IntegerLike }
  | { lift: () => IntegerLike };

/** Python repr for string-conversion diagnostics (after underscore removal). */
function integerStringRepr(value: string): string {
  const quote = value.includes("'") && !value.includes('"') ? '"' : "'";
  let result = quote;
  for (const character of value) {
    const code = character.codePointAt(0)!;
    if (character === quote || character === '\\') result += '\\' + character;
    else if (character === '\n') result += '\\n';
    else if (character === '\r') result += '\\r';
    else if (character === '\t') result += '\\t';
    else if (character !== ' ' && /[\p{C}\p{Z}]/u.test(character)) {
      result +=
        code <= 255
          ? `\\x${code.toString(16).padStart(2, '0')}`
          : code <= 65535
            ? `\\u${code.toString(16).padStart(4, '0')}`
            : `\\U${code.toString(16).padStart(8, '0')}`;
    } else result += character;
  }
  return result + quote;
}

/** integer.pyx:7402-7537, mpz_set_str_python and GMP string parsing. */
function mpz_set_str_python(input: string, base: bigint): bigint {
  if (base < -2147483648n || base > 2147483647n) {
    throw new OverflowError('value too large to convert to int');
  }
  if (base !== 0n && (base < 2n || base > 36n)) {
    throw new ValueError(`base (=${base}) must be 0 or between 2 and 36`);
  }
  // Integer.__init__ strips underscores before passing the C string to GMP.
  const source = input.replaceAll('_', '').split('\0', 1)[0]!;
  const invalid = (): never => {
    throw new SageTypeError(`unable to convert ${integerStringRepr(source)} to an integer`);
  };
  let text = source.replace(/^ +/, '');
  let sign = 1n;
  if (text[0] === '-' || text[0] === '+') {
    if (text[0] === '-') sign = -1n;
    text = text.slice(1);
  }
  text = text.replace(/^ +/, '');
  if (base === 0n) {
    if (/^0[bBoOxX]/.test(text)) {
      base = text[1]!.toLowerCase() === 'b' ? 2n : text[1]!.toLowerCase() === 'o' ? 8n : 16n;
      text = text.slice(2);
    } else base = 10n;
  }
  text = text.replace(/^ +/, '');
  if (text[0] === '-' || text[0] === '+') return invalid();
  // GMP accepts ASCII whitespace and an optional minus sign; the wrapper
  // above deliberately strips only ordinary spaces before checking prefixes.
  text = text.replace(/^[ \t\n\r\v\f]+|[ \t\n\r\v\f]+$/g, '');
  if (text[0] === '-') {
    sign = -sign;
    text = text.slice(1);
    if (!text || /^[ \t\n\r\v\f]/.test(text)) return invalid();
  }
  text = text.replace(/[ \t\n\r\v\f]/g, '').toLowerCase();
  if (!text || !/^[0-9a-z]+$/.test(text)) return invalid();
  const digits = [...text].map((c) => c.charCodeAt(0) - (c <= '9' ? 48 : 87));
  if (digits.some((d) => BigInt(d) >= base)) return invalid();
  if (base === 10n) return sign * BigInt(text);
  if (base === 2n || base === 8n || base === 16n) {
    return sign * BigInt((base === 2n ? '0b' : base === 8n ? '0o' : '0x') + text);
  }
  // GMP uses divide-and-conquer conversion for large non-power-of-two bases.
  const powers = new Map<number, bigint>();
  function convert(start: number, end: number): bigint {
    if (end - start <= 32) {
      let value = 0n;
      for (let i = start; i < end; i++) value = value * base + BigInt(digits[i]!);
      return value;
    }
    const mid = (start + end) >>> 1;
    const length = end - mid;
    let power = powers.get(length);
    if (power === undefined) {
      power = base ** BigInt(length);
      powers.set(length, power);
    }
    return convert(start, mid) * power + convert(mid, end);
  }
  return sign * convert(0, digits.length);
}

/** Shared Integer.__init__ coercion (integer.pyx:698-804). */
function integerValue(value: unknown, base: IntegerLike = 0n): bigint {
  if (value === undefined || value === null) return 0n;
  if (typeof value === 'bigint') return value;
  if (typeof value === 'boolean') return value ? 1n : 0n;
  if (typeof value === 'number') {
    if (Number.isNaN(value)) throw new ValueError('cannot convert float NaN to integer');
    if (!Number.isFinite(value))
      throw new OverflowError('cannot convert float infinity to integer');
    if (!Number.isInteger(value))
      throw new SageTypeError('cannot convert non-integral float to integer');
    return BigInt(value);
  }
  if (value instanceof Integer) return value.value;
  if (value instanceof Rational) {
    if (!value.isInteger()) throw new SageTypeError('no conversion of this rational to integer');
    return value.numerator;
  }
  if (typeof value === 'string') return mpz_set_str_python(value, toBigInt(base));
  if (Array.isArray(value)) {
    const radix = toBigInt(base);
    if (radix <= 1n) throw new SageTypeError("unable to coerce <class 'list'> to an integer");
    const digits = value.map((digit) => integerValue(digit));
    if (radix === 2n && digits.every((d) => d === 0n || d === 1n)) {
      return digits.length ? BigInt('0b' + digits.reverse().join('')) : 0n;
    }
    let result = 0n;
    for (let i = 0; i < digits.length; i++) result += digits[i]! * radix ** BigInt(i);
    return result;
  }
  if (typeof value === 'object') {
    const poly = value as {
      coeffs?: unknown[];
      degree?: () => number;
      getCoeff?: (i: number) => unknown;
      parent?: { base_ring: { zero(): unknown } };
    };
    if (Array.isArray(poly.coeffs) && poly.degree && poly.getCoeff && poly.parent) {
      const zero = poly.parent.base_ring.zero() as { value?: unknown };
      // Sections of ZZ -> QQ/prime coefficients -> polynomial injection.
      if (
        zero instanceof Rational ||
        zero instanceof Integer ||
        typeof zero.value === 'bigint' ||
        typeof zero.value === 'number'
      ) {
        if (poly.degree() > 0) throw new SageTypeError(`${value} is not a constant polynomial`);
        return integerValue(poly.getCoeff(0));
      }
    }
    const nativeFraction = _fraction_native_integer(value);
    if (nativeFraction !== undefined) return nativeFraction;
    const object = value as { _integer_?: (parent: IntegerRing) => unknown; lift?: () => unknown };
    if (typeof object._integer_ === 'function') {
      const result = object._integer_(ZZ);
      if (typeof result === 'bigint') return result;
      if (result instanceof Integer) return result.value;
      throw new SageTypeError('integer conversion hook must return an Integer');
    }
    if (typeof object.lift === 'function') {
      const result = object.lift();
      if (typeof result === 'bigint') return result;
      if (result instanceof Integer) return result.value;
    }
  }
  // An unrecognized plain object models Python's bare object() in scalar APIs.
  const typeName = typeof value === 'object' &&
    (Object.getPrototypeOf(value) === Object.prototype || Object.getPrototypeOf(value) === null)
    ? "<class 'object'>" : typeof value;
  throw new SageTypeError(`unable to coerce ${typeName} to an integer`);
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
   * Coerce a value to an integer; a missing value is zero.
   * Strings accept a base (0 auto-detects prefixes); arrays are digit lists.
   * @see Deviation: Scalar Constructor Adaptations
   */
  __call__(x?: IntegerInput, base?: IntegerLike): bigint {
    // Q_to_Z inherits Map._call_with_args, which rejects every extra argument.
    // Integer(QQ(...), base) goes straight through __init__ and differs here.
    if (base !== undefined && x instanceof Rational) {
      throw new NotImplementedError(
        "_call_with_args not overridden to accept arguments for <class 'sage.rings.rational.Q_to_Z'>"
      );
    }
    return integerValue(x, base);
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

    // `integer_ring.pyx:801` computes `den` UNCONDITIONALLY, at the top of
    // `_randomize_mpz`, before any distribution branch:
    //
    //     cdef int den = rstate.c_random()-SAGE_RAND_MAX/2
    //     if den == 0: den = 1
    //
    // Every call therefore burns exactly one 31-bit draw, whether or not the
    // `1/n` branch goes on to use it.  Computing it lazily inside that branch
    // put the uniform / mpz_rrandomb / gaussian streams one draw ahead of Sage
    // forever.
    let den = BigInt(rstate.c_random() - Math.floor(SAGE_RAND_MAX / 2));
    if (den === 0n) {
      den = 1n;
    }

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
      // `integer_ring.pyx:829` pins `algorithm='uniform+logtable'`; the
      // constructor's own default picks a table/online algorithm from
      // sigma*tau, which draws a different number of words.
      const sampler = new DiscreteGaussianDistributionIntegerSampler({
        sigma,
        algorithm: 'uniform+logtable',
      });
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
      // `integer_ring.pyx:822` calls GMP's `mpz_rrandomb`, the *runs*
      // generator, not `mpz_urandomb`.
      return rstate.random_bits_rrandomb(bits);
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
  readonly #bitLength: bigint;

  /**
   * Construct an Integer with Sage's scalar, string and digit-list coercion.
   * @see Deviation: Scalar Constructor Adaptations
   */
  constructor(value?: IntegerInput, base: IntegerLike = 0n) {
    this.value = integerValue(value, base);
    // GMP stores magnitude-size metadata beside the limbs. BigInt exposes no
    // equivalent query, so capture it once when constructing the immutable wrapper.
    if (value instanceof Integer) {
      this.#bitLength = value.#bitLength;
    } else {
      const magnitude = this.value < 0n ? -this.value : this.value;
      this.#bitLength = magnitude === 0n ? 0n : BigInt(magnitude.toString(2).length);
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
   * Return the prime factorization, or a partial trial factorization with limit.
   *
   * @see Deviation: PARI Integer Factorization (parigp-ts)
   * @see Deviation: Bounded factorization for elliptic hybrid orders
   */
  factor(options?: { limit?: IntegerLike }): Factorization {
    if (this.value === 0n) throw new ArithmeticError('factorization of 0 is not defined');
    if (options?.limit !== undefined) return factor_trial_division(this.value, options.limit);
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
   * Reduce modulo the ideal generated by n, whose generator in ZZ is |n|.
   * The zero ideal leaves self unchanged. Unlike __mod__, negative n does
   * not make the representative negative.
   * @see Reference: sage/structure/element.pyx:mod; sage/rings/ideal.py:reduce
   */
  mod(n: Integer | bigint): Integer {
    const other = n instanceof Integer ? n.value : n;
    if (other === 0n) return this;
    return new Integer(fdiv_r(this.value, other < 0n ? -other : other));
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
    if (this.value === 0n) return 0n;
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
   * @see Reference: sage/rings/integer.pyx:nbits
   * @see Deviation: Valuation dispatch and native GMP factor removal
   */
  nbits(): bigint {
    return this.bit_length();
  }

  /**
   * Return the valuation of this integer at prime p.
   * This is the largest k such that p^k divides this integer.
   * @see Reference: sage/rings/integer.pyx:valuation
   * @see Deviation: Valuation dispatch and native GMP factor removal
   */
  valuation(p: IntegerLike): bigint | 'Infinity' {
    const prime = toBigInt(p);
    // Integer(p) conversion precedes _valuation; zero precedes base validation.
    if (this.value === 0n) return 'Infinity';
    if (prime < 2n) {
      throw new ValueError(
        'You can only compute the valuation with respect to a integer larger than 1.'
      );
    }
    return mpz_remove(this.value, prime)[0];
  }

  // Arithmetic operations returning Integer
  add(n: Integer | bigint): Integer;
  add(n: Rational): Rational;
  add(n: Integer | bigint | Rational): Integer | Rational;
  add(n: Integer | bigint | Rational): Integer | Rational {
    if (n instanceof Rational) return new Rational(this.value).add(n);
    const other = n instanceof Integer ? n.value : n;
    return new Integer(this.value + other);
  }

  sub(n: Integer | bigint): Integer;
  sub(n: Rational): Rational;
  sub(n: Integer | bigint | Rational): Integer | Rational;
  sub(n: Integer | bigint | Rational): Integer | Rational {
    if (n instanceof Rational) return new Rational(this.value).sub(n);
    const other = n instanceof Integer ? n.value : n;
    return new Integer(this.value - other);
  }

  mul(n: Integer | bigint): Integer;
  mul(n: Rational): Rational;
  mul(n: Integer | bigint | Rational): Integer | Rational;
  mul(n: Integer | bigint | Rational): Integer | Rational {
    if (n instanceof Rational) return new Rational(this.value).mul(n);
    const other = n instanceof Integer ? n.value : n;
    return new Integer(this.value * other);
  }

  neg(): Integer {
    return new Integer(-this.value);
  }

  /** Exact integer powers; negative exponents belong to QQ, as in integer.pyx:_pow_long. */
  pow(n: IntegerLike): Integer | Rational {
    const nBig = toBigInt(n);
    if (nBig < 0n) {
      if (this.value === 0n) throw new ZeroDivisionError('rational division by zero');
      return new Rational(1n, this.value ** -nBig);
    }
    return new Integer(this.value ** nBig);
  }

  // Comparison
  eq(n: Integer | bigint | Rational): boolean;
  eq(n: Integer | bigint): boolean;
  eq(n: Integer | bigint | Rational): boolean {
    if (n instanceof Rational) return new Rational(this.value).eq(n);
    const other = n instanceof Integer ? n.value : n;
    return this.value === other;
  }

  lt(n: Integer | bigint | Rational): boolean;
  lt(n: Integer | bigint): boolean;
  lt(n: Integer | bigint | Rational): boolean {
    if (n instanceof Rational) return new Rational(this.value).lt(n);
    const other = n instanceof Integer ? n.value : n;
    return this.value < other;
  }

  le(n: Integer | bigint | Rational): boolean;
  le(n: Integer | bigint): boolean;
  le(n: Integer | bigint | Rational): boolean {
    if (n instanceof Rational) return new Rational(this.value).le(n);
    const other = n instanceof Integer ? n.value : n;
    return this.value <= other;
  }

  gt(n: Integer | bigint | Rational): boolean;
  gt(n: Integer | bigint): boolean;
  gt(n: Integer | bigint | Rational): boolean {
    if (n instanceof Rational) return new Rational(this.value).gt(n);
    const other = n instanceof Integer ? n.value : n;
    return this.value > other;
  }

  ge(n: Integer | bigint | Rational): boolean;
  ge(n: Integer | bigint): boolean;
  ge(n: Integer | bigint | Rational): boolean {
    if (n instanceof Rational) return new Rational(this.value).ge(n);
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
    if (nBig < -2147483648n || nBig > 2147483647n) {
      throw new OverflowError('value too large to convert to int');
    }
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

    const [magnitude, isExact] = mpz_root_positive(absVal, nBig);
    const root = isNegative ? -magnitude : magnitude;

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
   * @returns Floor of log_b(self), or negative infinity for zero
   * @throws {ValueError} If self < 0 or b < 2
   * @see Deviation: Infinity Representation
   * @see Reference: sage/rings/integer.pyx:exact_log
   */
  exact_log(b: Integer | bigint): bigint | '-Infinity' {
    const base = b instanceof Integer ? b.value : b;

    // Sage handles zero before validating the signs and base (integer.pyx:2803).
    if (this.value === 0n) return '-Infinity';
    if (this.value < 0n || base <= 0n) {
      throw new ValueError('self must be nonnegative and m must be positive');
    }
    if (base < 2n) {
      throw new ValueError('m must be at least 2');
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
   * NOTE: For positive D, this is the ordinary class number, which may be
   * half the narrow class number. For negative D, the two agree.
   *
   * @returns Class number
   * @throws {ValueError} If self is a perfect square or not congruent to 0 or 1 mod 4
   * @see Reference: sage/rings/integer.pyx:class_number
   * @see Deviation: Integer Quadratic Class Number Backend
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

    // Use the ported PARI quadratic class-group backend instead of a finite
    // table. It computes the ordinary class number for either sign of D.
    return quadclassno(D);
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
   * @see Deviation: GMP factorial kernels and native resource bounds
   */
  factorial(): Integer {
    if (this.value < 0n) throw new ValueError('factorial only defined for nonnegative integers');
    if (this.value >= 1n << 64n) throw new OverflowError('argument too large for factorial');
    return new Integer(mpz_fac_ui(this.value));
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
   * Negative inputs have infinitely many set bits in two's complement.
   * @see Deviation: Infinity Representation
   *
   * @returns Number of 1 bits
   * @see Reference: sage/rings/integer.pyx:popcount
   */
  popcount(): bigint | 'Infinity' {
    if (this.value < 0n) return 'Infinity';

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
   * For negative integers, returns the same Infinity sentinel as popcount.
   *
   * @returns Hamming weight
   * @see Deviation: Infinity Representation
   * @see Reference: sage/rings/integer.pyx:hamming_weight
   */
  hamming_weight(): bigint | 'Infinity' {
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
  log(b?: Integer | bigint): Integer | '-Infinity' {
    if (b === undefined) {
      if (this.value === 0n) return '-Infinity';
      // Natural log - use the existing overflow-safe approximation
      if (this.value <= 0n) {
        throw new ValueError('log of non-positive number');
      }
      const ln = this.real_log();
      return new Integer(BigInt(Math.floor(ln)));
    }

    const result = this.exact_log(b);
    return result === '-Infinity' ? result : new Integer(result);
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
  ord(p: Integer | bigint): bigint | 'Infinity' {
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
    if (pVal === 2n) return true;
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
    return this.#bitLength;
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
   * Return the exact reciprocal in QQ (including rational 1 and -1 for units).
   * @see Reference: sage/rings/integer.pyx:7031 (__invert__)
   */
  __invert__(): Rational {
    if (this.value === 0n) throw new ZeroDivisionError('rational division by zero');
    return new Rational(1n, this.value);
  }

  /**
   * Return the n-th Bell number if self = n.
   *
   * Bell numbers count the number of partitions of a set.
   * B_0 = 1, B_1 = 1, B_2 = 2, B_3 = 5, B_4 = 15, ...
   *
   * @returns The n-th Bell number
   * @throws {ArithmeticError} If self is negative
   * @see Reference: sage/combinat/combinat.py:150 (bell_number)
   */
  bell_number(): Integer {
    const n = this.value;

    if (n < 0n) {
      throw new ArithmeticError('Bell numbers not defined for negative indices');
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
   * @returns The n-th Catalan number, 0 for n < -1, and -1/2 for n = -1
   * @see Reference: sage/combinat/combinat.py:398 (catalan_number)
   */
  catalan_number(): Integer | Rational {
    const n = this.value;

    if (n < -1n) return new Integer(0n);
    if (n === -1n) return new Rational(-1n, 2n);

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
   * @see Deviation: Arithmetic Functions Not Delegated to PARI/FLINT
   */
  number_of_partitions(): Integer {
    const n = this.value;

    if (n < 0n) {
      throw new ValueError(`n (=${n}) must be a nonnegative integer`);
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
    if (this.value <= 2n) throw new ValueError('no prime less than 2');
    return new Integer(precprime(this.value - 1n));
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
   * @see Deviation: Arithmetic Functions Not Delegated to PARI/FLINT
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
