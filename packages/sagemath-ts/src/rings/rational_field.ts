/**
 * @module sage/rings/rational_field
 * @description Field QQ of rational numbers
 *
 * Port of: sage/rings/rational_field.py
 * Reference: reference/sage/src/sage/rings/rational_field.py
 */

import { gcd, is_prime, legendre_symbol, next_prime } from '../arith/misc.js';
import { TypeError, ValueError } from '../errors.js';
import { current_randstate } from '../misc/randstate.js';
import { type IntegerLike, toBigInt } from '../types/coercion.js';
import { type IntegerRing, ZZ } from './integer_ring.js';
import { Rational } from './rational.js';

/**
 * The class RationalField represents the field QQ of rational numbers.
 *
 * This is a singleton class - use the exported `QQ` constant.
 *
 * @example
 * ```typescript
 * import { QQ } from './rational_field';
 *
 * const r = QQ.__call__('3/4');        // Create from string
 * const s = QQ.__call__(7n);            // Create from bigint
 * const t = QQ.__call__(3, 4);          // Create from numerator/denominator
 *
 * console.log(QQ.is_field());           // true
 * console.log(QQ.characteristic());     // 0n
 * ```
 */
export class RationalField {
  private static instance: RationalField;

  private constructor() {}

  /**
   * Get the singleton instance.
   */
  static getInstance(): RationalField {
    if (!RationalField.instance) {
      RationalField.instance = new RationalField();
    }
    return RationalField.instance;
  }

  /**
   * Coerce a value to a rational number.
   *
   * Accepts:
   * - Rational objects (returned as-is)
   * - bigint values
   * - number values (must be finite)
   * - strings in format "n", "n/d", or decimal like "1.5"
   * - Two arguments (numerator, denominator)
   * - Tuple [numerator, denominator]
   *
   * @example
   * ```typescript
   * QQ.__call__(3n)           // Rational 3/1
   * QQ.__call__('3/4')        // Rational 3/4
   * QQ.__call__(3n, 4n)       // Rational 3/4
   * QQ.__call__([3n, 4n])     // Rational 3/4
   * QQ.__call__(1.5)          // Rational 3/2
   * ```
   */
  __call__(x: Rational): Rational;
  __call__(x: bigint | number | string): Rational;
  __call__(x: [bigint, bigint]): Rational;
  __call__(numerator: bigint, denominator: bigint): Rational;
  __call__(
    x: Rational | bigint | number | string | [bigint, bigint],
    denominator?: bigint
  ): Rational {
    // Two-argument form: numerator and denominator
    if (denominator !== undefined) {
      if (typeof x !== 'bigint') {
        throw new TypeError('numerator must be a bigint when denominator is provided');
      }
      return new Rational(x, denominator);
    }

    // Already a Rational
    if (x instanceof Rational) {
      return x;
    }

    // Tuple form [numerator, denominator]
    if (Array.isArray(x)) {
      if (x.length !== 2) {
        throw new TypeError('tuple must have exactly 2 elements');
      }
      return new Rational(x[0], x[1]);
    }

    // Use Rational.from for other types
    return Rational.from(x);
  }

  /**
   * Return the zero element of QQ.
   */
  zero(): Rational {
    return Rational.zero();
  }

  /**
   * Return the one element of QQ.
   */
  one(): Rational {
    return Rational.one();
  }

  /**
   * Return True since QQ is a field.
   */
  is_field(): boolean {
    return true;
  }

  /**
   * Return the characteristic of QQ, which is 0.
   */
  characteristic(): bigint {
    return 0n;
  }

  /**
   * Return True since QQ is a ring.
   */
  is_ring(): boolean {
    return true;
  }

  /**
   * Return True since QQ is an integral domain.
   */
  is_integral_domain(): boolean {
    return true;
  }

  /**
   * Return True since QQ is a prime field.
   */
  is_prime_field(): boolean {
    return true;
  }

  /**
   * Return True since QQ is the base field for QQ.
   */
  is_absolute(): boolean {
    return true;
  }

  /**
   * Return False since QQ is infinite.
   *
   * @see Reference: sage/rings/rational_field.py:is_finite
   */
  is_finite(): boolean {
    return false;
  }

  /**
   * Return the degree of QQ, which is 1.
   */
  degree(): bigint {
    return 1n;
  }

  /**
   * Return the absolute degree of QQ, which is 1.
   */
  absolute_degree(): bigint {
    return 1n;
  }

  /**
   * Return the number of generators of QQ, which is 1.
   */
  ngens(): bigint {
    return 1n;
  }

  /**
   * Return a tuple of generators of QQ, which is just (1,).
   */
  gens(): [Rational] {
    return [this.one()];
  }

  /**
   * Return the n-th generator of QQ.
   *
   * @param n - Must be 0
   * @throws {ValueError} If n is not 0
   */
  gen(n: number = 0): Rational {
    if (n === 0) {
      return this.one();
    }
    throw new ValueError('n must be 0');
  }

  /**
   * Return a random element of QQ.
   *
   * @param numBound - Bound on the absolute value of the numerator
   * @param denBound - Bound on the value of the denominator (default: numBound)
   */
  random_element(numBound: bigint = 100n, denBound?: bigint): Rational {
    if (denBound === undefined) {
      denBound = numBound;
    }

    if (numBound <= 0n || denBound <= 0n) {
      throw new ValueError('bounds must be positive');
    }

    const rstate = current_randstate();

    // Generate random numerator in [-numBound, numBound]
    const numRange = 2n * numBound + 1n;
    const num = rstate.random_below(numRange) - numBound;

    // Generate random denominator in [1, denBound]
    const den = rstate.randint(1n, denBound);

    return new Rational(num, den);
  }

  /**
   * Return a typical element of QQ (used for testing).
   */
  an_element(): Rational {
    return new Rational(1n, 2n);
  }

  /**
   * Generate some elements of QQ (used for testing).
   */
  *some_elements(): Generator<Rational> {
    yield this.an_element();
    yield this.an_element().neg();
    yield this.an_element().inv();
    yield this.an_element().neg().inv();
    yield this.zero();
    yield this.one();
    yield this.one().neg();
    yield new Rational(42n, 1n);

    // Additional elements
    for (let n = 1n; n <= 10n; n++) {
      const a = 2n * n;
      const b = 2n * n + 1n;
      yield new Rational(a, b);
      yield new Rational(-a, b);
      yield new Rational(b, a);
      yield new Rational(-b, a);
    }
  }

  /**
   * Return a primitive n-th root of unity in QQ.
   *
   * Only n=1 (returns 1) and n=2 (returns -1) are valid.
   *
   * @param n - The order of the root of unity
   * @throws {ValueError} If n is not 1 or 2
   */
  zeta(n: number = 2): Rational {
    if (n === 1) {
      return this.one();
    }
    if (n === 2) {
      return this.one().neg();
    }
    throw new ValueError('no n-th root of unity in rational field');
  }

  /**
   * Return the discriminant of QQ, which is 1.
   */
  discriminant(): bigint {
    return 1n;
  }

  /**
   * Return the absolute discriminant of QQ, which is 1.
   *
   * @see Reference: sage/rings/rational_field.py:absolute_discriminant
   */
  absolute_discriminant(): bigint {
    return this.discriminant();
  }

  /**
   * Return the relative discriminant of QQ, which is 1.
   *
   * @see Reference: sage/rings/rational_field.py:relative_discriminant
   */
  relative_discriminant(): bigint {
    return this.discriminant();
  }

  /**
   * Return the class number of QQ, which is 1.
   */
  class_number(): bigint {
    return 1n;
  }

  /**
   * Return the signature of QQ, which is (1, 0).
   */
  signature(): [bigint, bigint] {
    return [1n, 0n];
  }

  /**
   * Return the order of QQ, which is Infinity.
   *
   * @see Reference: sage/rings/rational_field.py:order
   */
  order(): 'Infinity' {
    return 'Infinity';
  }

  /**
   * Return the maximal order of QQ, which is ZZ.
   *
   * @see Reference: sage/rings/rational_field.py:maximal_order
   */
  maximal_order(): IntegerRing {
    return ZZ;
  }

  /**
   * Return the ring of integers of QQ, which is ZZ.
   *
   * @see Reference: sage/rings/rational_field.py:maximal_order
   */
  ring_of_integers(): IntegerRing {
    return this.maximal_order();
  }

  /**
   * Return the number field associated to QQ (which is QQ itself).
   *
   * @see Reference: sage/rings/rational_field.py:number_field
   */
  number_field(): RationalField {
    return this;
  }

  /**
   * Return the power basis for QQ over its base field.
   *
   * The power basis is always [1] for the rational field.
   *
   * @see Reference: sage/rings/rational_field.py:power_basis
   */
  power_basis(): [Rational] {
    return [this.gen()];
  }

  /**
   * Iterate over rational numbers ordered by height.
   *
   * The first elements are: 0, 1, -1, 1/2, -1/2, 2, -2, 1/3, -1/3, 3, -3, ...
   *
   * @see Reference: sage/rings/rational_field.py:__iter__
   *
   * @example
   * ```typescript
   * // Get first 10 rationals by height
   * const rationals = [];
   * for (const r of QQ) {
   *   rationals.push(r);
   *   if (rationals.length >= 10) break;
   * }
   * // [0, 1, -1, 1/2, -1/2, 2, -2, 1/3, -1/3, 3]
   * ```
   */
  *[Symbol.iterator](): Generator<Rational> {
    yield this.__call__(0n);
    yield this.__call__(1n);
    yield this.__call__(-1n);

    let height = 1n;
    while (true) {
      height = height + 1n;
      for (let other = 1n; other < height; other++) {
        if (gcd(height, other) === 1n) {
          yield new Rational(other, height);
          yield new Rational(-other, height);
          yield new Rational(height, other);
          yield new Rational(-height, other);
        }
      }
    }
  }

  /**
   * Range function for rational numbers, ordered by height.
   *
   * Returns a generator for rational numbers with heights in range(start, end).
   *
   * @param start - Start height (or end if end is undefined)
   * @param end - End height (exclusive)
   *
   * @see Reference: sage/rings/rational_field.py:range_by_height
   *
   * @example
   * ```typescript
   * // All rationals with height < 4
   * [...QQ.range_by_height(4)]
   * // [0, 1, -1, 1/2, -1/2, 2, -2, 1/3, -1/3, 3, -3, 2/3, -2/3, 3/2, -3/2]
   *
   * // Rationals with height exactly 2
   * [...QQ.range_by_height(2, 3)]
   * // [1/2, -1/2, 2, -2]
   * ```
   */
  *range_by_height(start: IntegerLike, end?: IntegerLike): Generator<Rational> {
    let startVal = toBigInt(start);
    let endVal: bigint;

    if (end === undefined) {
      endVal = startVal;
      startVal = 1n;
    } else {
      endVal = toBigInt(end);
    }

    // Ensure start is at least 1
    if (startVal < 1n) {
      startVal = 1n;
    }

    for (let height = startVal; height < endVal; height++) {
      if (height === 1n) {
        yield this.__call__(0n);
        yield this.__call__(1n);
        yield this.__call__(-1n);
      } else {
        for (let other = 1n; other < height; other++) {
          if (gcd(height, other) === 1n) {
            yield new Rational(other, height);
            yield new Rational(-other, height);
            yield new Rational(height, other);
            yield new Rational(-height, other);
          }
        }
      }
    }
  }

  /**
   * Iterator yielding all primes less than or equal to B.
   *
   * This function exists for compatibility with the related number
   * field method, though it returns prime integers, not ideals.
   *
   * @param B - Positive integer upper bound on primes
   *
   * @see Reference: sage/rings/rational_field.py:primes_of_bounded_norm_iter
   *
   * @example
   * ```typescript
   * [...QQ.primes_of_bounded_norm_iter(10)]
   * // [2n, 3n, 5n, 7n]
   * ```
   */
  *primes_of_bounded_norm_iter(B: IntegerLike): Generator<bigint> {
    const bound = toBigInt(B);

    if (bound < 2n) {
      return;
    }

    // Generate primes up to bound
    let p = 2n;
    while (p <= bound) {
      yield p;
      p = next_prime(p);
    }
  }

  /**
   * Return generators of the group QQ(S,m).
   *
   * QQ(S,m) is the subgroup of QQ^* modulo (QQ^*)^m consisting of elements a
   * such that the valuation of a is divisible by m at all primes not in S.
   * It is equal to the group of S-units modulo m-th powers.
   *
   * @param S - Set of primes
   * @param m - Positive integer
   * @param options - Options object
   * @param options.proof - Ignored (for compatibility)
   * @param options.orders - If true, return [generators, orders]
   *
   * @see Reference: sage/rings/rational_field.py:selmer_generators
   *
   * @example
   * ```typescript
   * QQ.selmer_generators([], 2)
   * // [-1n]
   *
   * QQ.selmer_generators([3n], 2)
   * // [-1n, 3n]
   *
   * QQ.selmer_generators([2n, 3n, 5n, 7n], 2, { orders: true })
   * // [[-1n, 2n, 3n, 5n, 7n], [2n, 2n, 2n, 2n, 2n]]
   * ```
   */
  selmer_generators(
    S: IntegerLike[],
    m: IntegerLike,
    options?: { proof?: boolean; orders?: boolean }
  ): bigint[] | [bigint[], bigint[]] {
    const mVal = toBigInt(m);
    const primes = S.map((p) => toBigInt(p));
    const orders = options?.orders ?? false;

    const gens: bigint[] = [...primes];
    const ords: bigint[] = new Array(primes.length).fill(mVal);

    // When m is even, -1 is a generator of order 2
    if (mVal % 2n === 0n) {
      gens.unshift(-1n);
      ords.unshift(2n);
    }

    return orders ? [gens, ords] : gens;
  }

  /**
   * Return an iterator through elements of the finite group QQ(S,m).
   *
   * @param S - Set of primes
   * @param m - Positive integer
   * @param _proof - Ignored (for compatibility)
   *
   * @see Reference: sage/rings/rational_field.py:selmer_group_iterator
   *
   * @example
   * ```typescript
   * [...QQ.selmer_group_iterator([], 2)]
   * // [Rational(1), Rational(-1)]
   *
   * [...QQ.selmer_group_iterator([2n], 2)]
   * // [Rational(1), Rational(2), Rational(-1), Rational(-2)]
   * ```
   */
  *selmer_group_iterator(S: IntegerLike[], m: IntegerLike, _proof?: boolean): Generator<Rational> {
    const [gens, ords] = this.selmer_generators(S, m, { orders: true }) as [bigint[], bigint[]];

    // Generate all products of generators^exponent for exponents in [0, order)
    const indices = ords.map((o) => Number(o));

    // Recursive generator for Cartesian product
    function* cartesianProduct(arrays: number[], current: number[] = []): Generator<number[]> {
      if (current.length === arrays.length) {
        yield current;
        return;
      }
      const idx = current.length;
      for (let i = 0; i < arrays[idx]!; i++) {
        yield* cartesianProduct(arrays, [...current, i]);
      }
    }

    for (const exponents of cartesianProduct(indices)) {
      let result = 1n;
      for (let i = 0; i < gens.length; i++) {
        result *= gens[i]! ** BigInt(exponents[i]!);
      }
      yield new Rational(result, 1n);
    }
  }

  /**
   * Return the valuation of the quadratic defect of a at p.
   *
   * @param a - An element of QQ
   * @param p - A prime number
   * @param check - If true, check that p is prime
   *
   * @see Reference: sage/rings/rational_field.py:quadratic_defect
   *
   * @example
   * ```typescript
   * QQ.quadratic_defect(0, 7)   // 'Infinity'
   * QQ.quadratic_defect(5, 7)   // 0n
   * QQ.quadratic_defect(5, 2)   // 2n
   * QQ.quadratic_defect(5, 5)   // 1n
   * ```
   */
  quadratic_defect(
    a: Rational | IntegerLike,
    p: IntegerLike,
    check: boolean = true
  ): bigint | 'Infinity' {
    const pVal = toBigInt(p);

    // Convert a to Rational
    const aRat = a instanceof Rational ? a : new Rational(toBigInt(a), 1n);

    if (check && !is_prime(pVal)) {
      throw new ValueError(`${pVal} must be prime`);
    }

    if (aRat.isZero()) {
      return 'Infinity';
    }

    const [v, u] = aRat.val_unit(pVal);

    if (v === 'Infinity') {
      return 'Infinity';
    }

    // If valuation is odd, return it
    if (v % 2n === 1n) {
      return v;
    }

    // For p != 2
    if (pVal !== 2n) {
      // u is a Rational, get its integer value for Legendre symbol
      const uInt = u.numerator; // Since unit has denominator coprime to p
      if (legendre_symbol(uInt, pVal) === 1n) {
        return 'Infinity';
      } else {
        return v;
      }
    }

    // For p = 2
    // Get u mod 8
    const uInt = u.numerator;
    const uMod8 = ((uInt % 8n) + 8n) % 8n;

    if (uMod8 === 1n) {
      return 'Infinity';
    }
    if (uMod8 === 5n) {
      return v + 2n;
    }
    if (uMod8 === 3n || uMod8 === 7n) {
      return v + 1n;
    }

    // Should not reach here for valid inputs
    return v;
  }

  /**
   * String representation.
   */
  toString(): string {
    return 'Rational Field';
  }

  /**
   * LaTeX representation.
   */
  _latex_(): string {
    return '\\Bold{Q}';
  }
}

/**
 * The field of rational numbers QQ.
 *
 * This is the singleton instance of RationalField.
 *
 * @example
 * ```typescript
 * import { QQ } from './rational_field';
 *
 * const r = QQ.__call__('3/4');
 * const s = QQ.__call__(7n, 3n);
 * console.log(r.add(s).toString()); // '37/12'
 * ```
 */
export const QQ = RationalField.getInstance();

/**
 * Alias for QQ.
 */
export const Q = QQ;

/**
 * Check if x is the rational field QQ.
 *
 * @param x - Value to check
 * @returns true if x is the RationalField singleton
 *
 * @see Reference: sage/rings/rational_field.py:is_RationalField
 *
 * @example
 * ```typescript
 * is_RationalField(QQ)  // true
 * is_RationalField(ZZ)  // false
 * ```
 */
export function is_RationalField(x: unknown): x is RationalField {
  return x instanceof RationalField;
}
