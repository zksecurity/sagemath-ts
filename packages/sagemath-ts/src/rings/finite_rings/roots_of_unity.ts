/**
 * @module sage/rings/finite_rings/roots_of_unity
 * @description Roots of unity, multiplicative subgroups, and FFT domains for ZK applications
 *
 * This module provides functionality for working with roots of unity in finite fields,
 * which is essential for NTT/FFT-based polynomial multiplication in ZK-SNARKs and STARKs.
 *
 * Port of concepts from:
 * - sage/rings/polynomial/cyclotomic.pyx
 * - sage/rings/finite_rings/finite_field_base.pyx (multiplicative_generator, zeta)
 *
 * Key concepts:
 * - In GF(p), n-th roots of unity exist iff n | (p-1)
 * - In GF(p^k), n-th roots of unity exist iff n | (p^k - 1)
 * - A primitive n-th root omega satisfies: omega^n = 1 and omega^k != 1 for 0 < k < n
 * - For FFT, we need n to be a power of 2, and the domain consists of [1, omega, omega^2, ..., omega^(n-1)]
 */

import { divisors, euler_phi, factor, gcd, moebius, power_mod } from '../../arith/misc.js';
import { ValueError } from '../../errors.js';
import type { CoefficientRing, RingElement } from '../polynomial/polynomial_element.js';
import { Polynomial } from '../polynomial/polynomial_element.js';
import type { PolynomialRing } from '../polynomial/polynomial_ring.js';

/**
 * Type for finite field elements that support the required operations.
 */
export interface FiniteFieldElementLike extends RingElement {
  parent: FiniteFieldLike;
  value?: bigint;
  pow(n: bigint | number): FiniteFieldElementLike;
  mul(other: FiniteFieldElementLike): FiniteFieldElementLike;
  isOne(): boolean;
  isZero(): boolean;
  eq(other: FiniteFieldElementLike | number | bigint): boolean;
}

/**
 * Type for finite fields that support the required operations.
 */
export interface FiniteFieldLike extends CoefficientRing<FiniteFieldElementLike> {
  characteristic: bigint;
  order: bigint;
  degree: number;
  one(): FiniteFieldElementLike;
  zero(): FiniteFieldElementLike;
  __call__(x: unknown): FiniteFieldElementLike;
  multiplicative_generator?(): FiniteFieldElementLike;
  gen?(): FiniteFieldElementLike;
}

/**
 * Check if n-th primitive roots of unity exist in the given finite field.
 *
 * For GF(p), n-th roots exist iff n | (p - 1).
 * For GF(p^k), n-th roots exist iff n | (p^k - 1).
 *
 * @param field - The finite field
 * @param n - The order of roots to check for
 * @returns True if primitive n-th roots of unity exist
 *
 * @example
 * ```typescript
 * const F17 = GF(17n);
 * has_primitive_root(F17, 16n)  // true, since 16 | 16
 * has_primitive_root(F17, 8n)   // true, since 8 | 16
 * has_primitive_root(F17, 5n)   // false, since 5 does not divide 16
 * ```
 */
export function has_primitive_root(field: FiniteFieldLike, n: bigint | number): boolean {
  const order = typeof n === 'number' ? BigInt(n) : n;

  if (order <= 0n) {
    return false;
  }

  if (order === 1n) {
    return true; // 1 is always a 1st root of unity
  }

  // The multiplicative group has order (field.order - 1)
  const groupOrder = field.order - 1n;

  // n-th roots exist iff n divides the group order
  return groupOrder % order === 0n;
}

/**
 * Find a primitive n-th root of unity in the given finite field.
 *
 * A primitive n-th root omega satisfies:
 * - omega^n = 1
 * - omega^k != 1 for all 0 < k < n
 *
 * Uses the formula: omega = g^((|F*|)/n) where g is a multiplicative generator.
 *
 * @param field - The finite field
 * @param n - The order of the root to find
 * @returns A primitive n-th root of unity
 * @throws {ValueError} If n-th roots do not exist in the field
 *
 * @example
 * ```typescript
 * const F17 = GF(17n);
 * const omega = primitive_nth_root(F17, 8n);
 * // omega^8 = 1, omega^4 != 1, omega^2 != 1, omega != 1
 * ```
 */
export function primitive_nth_root(
  field: FiniteFieldLike,
  n: bigint | number
): FiniteFieldElementLike {
  const order = typeof n === 'number' ? BigInt(n) : n;

  if (order <= 0n) {
    throw new ValueError('n must be positive');
  }

  if (order === 1n) {
    return field.one();
  }

  const groupOrder = field.order - 1n;

  // Check that n-th roots exist
  if (groupOrder % order !== 0n) {
    throw new ValueError(
      `No primitive ${order}-th root of unity exists in the field. ` +
        `${order} does not divide ${groupOrder} = |F*|.`
    );
  }

  // Find a multiplicative generator
  const g = findMultiplicativeGenerator(field);

  // omega = g^((|F*|)/n) is a primitive n-th root
  const exp = groupOrder / order;
  return g.pow(exp);
}

/**
 * Return all n-th roots of unity in the given finite field.
 *
 * Returns [1, omega, omega^2, ..., omega^(n-1)] where omega is a primitive n-th root.
 *
 * @param field - The finite field
 * @param n - The order of the roots
 * @returns Array of all n-th roots of unity
 * @throws {ValueError} If n-th roots do not exist in the field
 *
 * @example
 * ```typescript
 * const F7 = GF(7n);
 * const roots = roots_of_unity(F7, 3n);  // [1, 2, 4] (since 2^3 = 8 ≡ 1 mod 7)
 * ```
 */
export function roots_of_unity(
  field: FiniteFieldLike,
  n: bigint | number
): FiniteFieldElementLike[] {
  const order = typeof n === 'number' ? BigInt(n) : n;

  if (order <= 0n) {
    throw new ValueError('n must be positive');
  }

  if (order === 1n) {
    return [field.one()];
  }

  const omega = primitive_nth_root(field, order);
  const roots: FiniteFieldElementLike[] = [];

  let current = field.one();
  for (let i = 0n; i < order; i++) {
    roots.push(current);
    current = current.mul(omega);
  }

  return roots;
}

/**
 * Compute the multiplicative order of an element in its parent field.
 *
 * The multiplicative order is the smallest positive integer k such that element^k = 1.
 *
 * @param element - A non-zero field element
 * @returns The multiplicative order of the element
 * @throws {ValueError} If the element is zero
 *
 * @example
 * ```typescript
 * const F7 = GF(7n);
 * multiplicative_order(F7(3n))  // 6n (3 is a generator)
 * multiplicative_order(F7(2n))  // 3n (2^3 = 8 ≡ 1)
 * ```
 */
export function multiplicative_order(element: FiniteFieldElementLike): bigint {
  if (element.isZero()) {
    throw new ValueError('multiplicative order of zero is not defined');
  }

  if (element.isOne()) {
    return 1n;
  }

  const field = element.parent;
  const groupOrder = field.order - 1n;

  // The order divides |F*|
  // Find the smallest divisor that works
  const divs = divisors(groupOrder);

  for (const d of divs) {
    if (element.pow(d).isOne()) {
      return d;
    }
  }

  // Should never reach here for valid field elements
  return groupOrder;
}

/**
 * Find a multiplicative generator (primitive element) of the finite field.
 *
 * A multiplicative generator g satisfies: g^(|F*|) = 1 and g^k != 1 for 0 < k < |F*|.
 * In other words, g generates the entire multiplicative group.
 *
 * @param field - The finite field
 * @returns A multiplicative generator of the field
 */
export function findMultiplicativeGenerator(field: FiniteFieldLike): FiniteFieldElementLike {
  // If the field provides a multiplicative generator, use it
  if (field.multiplicative_generator) {
    return field.multiplicative_generator();
  }

  const groupOrder = field.order - 1n;

  // Factor the group order to check generator property efficiently
  const factors = factor(groupOrder);
  const primeFactors = factors.map(([p]) => p);

  // Try candidates starting from 2 (or the field generator for extensions)
  const candidates: FiniteFieldElementLike[] = [];

  // For extension fields, try the generator first
  if (field.gen && field.degree > 1) {
    candidates.push(field.gen());
  }

  // Try small integers
  for (let i = 2n; i < field.order && i < 1000n; i++) {
    candidates.push(field.__call__(i));
  }

  for (const g of candidates) {
    if (g.isZero()) continue;

    let isGenerator = true;
    for (const p of primeFactors) {
      const exp = groupOrder / p;
      if (g.pow(exp).isOne()) {
        isGenerator = false;
        break;
      }
    }

    if (isGenerator) {
      return g;
    }
  }

  throw new ValueError('Could not find multiplicative generator');
}

/**
 * Compute the n-th cyclotomic polynomial over a given ring.
 *
 * The n-th cyclotomic polynomial Phi_n(x) is the minimal polynomial of
 * primitive n-th roots of unity over the rationals. It has degree phi(n)
 * (Euler's totient function).
 *
 * Uses the formula: Phi_n(x) = prod_{d|n} (x^d - 1)^{mu(n/d)}
 * where mu is the Mobius function.
 *
 * @param n - The index of the cyclotomic polynomial
 * @param ring - The polynomial ring to use (optional, defaults to integers)
 * @returns The n-th cyclotomic polynomial
 *
 * @example
 * ```typescript
 * cyclotomic_polynomial(1)   // x - 1
 * cyclotomic_polynomial(2)   // x + 1
 * cyclotomic_polynomial(3)   // x^2 + x + 1
 * cyclotomic_polynomial(4)   // x^2 + 1
 * cyclotomic_polynomial(6)   // x^2 - x + 1
 * ```
 */
export function cyclotomic_polynomial<C extends RingElement>(
  n: bigint | number,
  ring?: PolynomialRing<C>
): number[] | Polynomial<C> {
  const order = typeof n === 'number' ? BigInt(n) : n;

  if (order < 1n) {
    throw new ValueError('n must be positive');
  }

  // Compute the coefficients
  const coeffs = cyclotomicCoefficients(order);

  // If no ring is specified, return the coefficients as numbers
  if (!ring) {
    return coeffs;
  }

  // Convert to a polynomial over the given ring
  const polyCoeffs = coeffs.map((c) => ring.base_ring.__call__(c) as C);
  return new Polynomial(polyCoeffs, ring);
}

/**
 * Compute the coefficients of the n-th cyclotomic polynomial.
 *
 * Uses the formula from Sage:
 * Phi_n(x) = prod_{d|n} (x^d - 1)^{mu(n/d)}
 *
 * For efficient computation, we use:
 * - For n with repeated prime factors: Phi_n(x) = Phi_rad(x^(n/rad)) where rad = radical(n)
 * - Multiply by (x^d - 1) when mu(n/d) = 1
 * - Divide by (x^d - 1) when mu(n/d) = -1
 *
 * @param n - The index of the cyclotomic polynomial
 * @returns Array of coefficients [a_0, a_1, ..., a_deg]
 */
function cyclotomicCoefficients(n: bigint): number[] {
  if (n === 1n) {
    return [-1, 1]; // x - 1
  }

  if (n === 2n) {
    return [1, 1]; // x + 1
  }

  // Factor n
  const factors = factor(n);

  // Check if n has repeated prime factors
  const hasRepeatedFactors = factors.some(([_p, e]) => e > 1n);

  if (hasRepeatedFactors) {
    // Phi_n(x) = Phi_rad(x^(n/rad)) where rad is the radical (squarefree part)
    let rad = 1n;
    for (const [p, _e] of factors) {
      rad *= p;
    }
    const power = n / rad;

    // Get coefficients of Phi_rad
    const radCoeffs = cyclotomicCoefficients(rad);

    // Substitute x -> x^power
    const deg = (radCoeffs.length - 1) * Number(power);
    const result: number[] = Array(deg + 1).fill(0);

    for (let i = 0; i < radCoeffs.length; i++) {
      result[i * Number(power)] = radCoeffs[i]!;
    }

    return result;
  }

  // n is squarefree - use the product formula
  // For prime n: Phi_n(x) = 1 + x + x^2 + ... + x^(n-1) = (x^n - 1)/(x - 1)
  if (factors.length === 1 && factors[0]![1] === 1n) {
    const deg = Number(n - 1n);
    return Array(deg + 1).fill(1);
  }

  // General squarefree case
  // We compute prod_{d|n} (x^d - 1)^{mu(n/d)}
  // Group divisors by mobius value and process accordingly

  const divs = divisors(n);
  const primes = factors.map(([p]) => p);

  // Subsets of primes determine divisors and their mobius values
  // mu(n/d) = mu(product of primes not in d) = (-1)^(number of primes not in d)

  // Start with f = 1
  let coeffs: number[] = [1];

  // Process divisors with mu(n/d) = 1 (multiply by x^d - 1)
  // These are d where n/d has even number of prime factors
  for (const d of divs) {
    const nd = n / d;
    const mu = moebius(nd);
    if (mu === 1n) {
      // Multiply by (x^d - 1)
      coeffs = multiplyByXdMinus1(coeffs, Number(d));
    }
  }

  // Process divisors with mu(n/d) = -1 (divide by x^d - 1)
  // These are d where n/d has odd number of prime factors
  // Process in reverse order (largest to smallest) for numerical stability
  const divsReversed = [...divs].reverse();
  for (const d of divsReversed) {
    const nd = n / d;
    const mu = moebius(nd);
    if (mu === -1n) {
      // Divide by (x^d - 1)
      coeffs = divideByXdMinus1(coeffs, Number(d));
    }
  }

  return coeffs;
}

/**
 * Multiply a polynomial by (x^d - 1).
 * If f(x) = sum a_i x^i, then f(x) * (x^d - 1) = sum a_i x^(i+d) - sum a_i x^i
 */
function multiplyByXdMinus1(coeffs: number[], d: number): number[] {
  const deg = coeffs.length - 1;
  const newDeg = deg + d;
  const result: number[] = Array(newDeg + 1).fill(0);

  for (let i = 0; i <= deg; i++) {
    // -a_i * x^i
    result[i] -= coeffs[i]!;
    // a_i * x^(i+d)
    result[i + d] += coeffs[i]!;
  }

  // Remove leading zeros
  while (result.length > 1 && result[result.length - 1] === 0) {
    result.pop();
  }

  return result;
}

/**
 * Divide a polynomial by (x^d - 1).
 * This is exact division - we assume (x^d - 1) divides the polynomial.
 *
 * If f(x) = (x^d - 1) * q(x), then we recover q(x).
 * Working from low to high degree:
 * f_i = q_{i-d} - q_i (where q_j = 0 for j < 0)
 * So q_i = q_{i-d} - f_i
 */
function divideByXdMinus1(coeffs: number[], d: number): number[] {
  const deg = coeffs.length - 1;
  const newDeg = deg - d;

  if (newDeg < 0) {
    // Result is 0 or input is smaller than divisor
    return [0];
  }

  const result: number[] = Array(newDeg + 1).fill(0);

  // Working from lowest degree to highest
  for (let i = 0; i <= newDeg; i++) {
    // q_i = q_{i-d} - f_i
    const qPrev = i >= d ? result[i - d]! : 0;
    result[i] = qPrev - coeffs[i]!;
  }

  // Remove leading zeros
  while (result.length > 1 && result[result.length - 1] === 0) {
    result.pop();
  }

  return result;
}

/**
 * An FFT domain for ZK applications.
 *
 * An FFT domain consists of the powers of a primitive root of unity:
 * {1, omega, omega^2, ..., omega^(n-1)}
 *
 * This is used for:
 * - Number Theoretic Transform (NTT)
 * - Polynomial commitment schemes
 * - STARK/SNARK proof systems
 *
 * @example
 * ```typescript
 * const F = GF(65537n);  // 65537 = 2^16 + 1, Fermat prime
 * const domain = new FFTDomain(F, 16);  // Domain of size 16
 *
 * // Get all domain elements
 * const elems = domain.elements();
 *
 * // Check membership
 * domain.contains(elems[0])  // true
 *
 * // Get a coset
 * const coset = domain.coset(F(7n));
 * ```
 */
export class FFTDomain {
  readonly field: FiniteFieldLike;
  readonly size: bigint;
  readonly generator: FiniteFieldElementLike;
  private _elements: FiniteFieldElementLike[] | null = null;

  /**
   * Create an FFT domain of the given size.
   *
   * @param field - The finite field
   * @param size - The size of the domain (must be a power of 2 that divides |F*|)
   * @throws {ValueError} If the size is not valid
   */
  constructor(field: FiniteFieldLike, size: bigint | number) {
    const domainSize = typeof size === 'number' ? BigInt(size) : size;

    if (domainSize <= 0n) {
      throw new ValueError('domain size must be positive');
    }

    // Check that size is a power of 2
    if ((domainSize & (domainSize - 1n)) !== 0n) {
      throw new ValueError(`domain size ${domainSize} must be a power of 2`);
    }

    // Check that roots of unity exist
    if (!has_primitive_root(field, domainSize)) {
      throw new ValueError(
        `No ${domainSize}-th roots of unity in the field. ` +
          `${domainSize} does not divide ${field.order - 1n}.`
      );
    }

    this.field = field;
    this.size = domainSize;
    this.generator = primitive_nth_root(field, domainSize);
  }

  /**
   * Return all elements of the domain: [1, omega, omega^2, ..., omega^(n-1)].
   */
  elements(): FiniteFieldElementLike[] {
    if (this._elements !== null) {
      return this._elements;
    }

    this._elements = [];
    let current = this.field.one();

    for (let i = 0n; i < this.size; i++) {
      this._elements.push(current);
      current = current.mul(this.generator);
    }

    return this._elements;
  }

  /**
   * Check if an element is in the domain.
   *
   * An element x is in the domain iff x^n = 1.
   */
  contains(element: FiniteFieldElementLike): boolean {
    return element.pow(this.size).isOne();
  }

  /**
   * Return a coset of the domain: {offset * omega^i | i = 0, 1, ..., n-1}.
   *
   * Cosets are used in FRI (Fast Reed-Solomon IOPP) for domain extension.
   */
  coset(offset: FiniteFieldElementLike): CosetDomain {
    return new CosetDomain(this.field, this.size, this.generator, offset);
  }

  /**
   * Return the vanishing polynomial of this domain: Z(x) = x^n - 1.
   *
   * This polynomial has all domain elements as roots.
   */
  vanishingPolynomial<C extends RingElement>(ring: PolynomialRing<C>): Polynomial<C> {
    // x^n - 1
    const coeffs: C[] = [];
    const n = Number(this.size);

    // Constant term: -1
    coeffs.push(ring.base_ring.__call__(-1) as C);

    // Zeros for x^1 through x^(n-1)
    for (let i = 1; i < n; i++) {
      coeffs.push(ring.base_ring.zero() as C);
    }

    // Leading coefficient: 1
    coeffs.push(ring.base_ring.one() as C);

    return new Polynomial(coeffs, ring);
  }

  /**
   * Get the size of the domain (for iteration).
   */
  get length(): number {
    return Number(this.size);
  }

  /**
   * Iterate over domain elements.
   */
  *[Symbol.iterator](): Iterator<FiniteFieldElementLike> {
    yield* this.elements();
  }

  /**
   * Return the log base 2 of the domain size.
   */
  get log2Size(): number {
    let size = this.size;
    let log = 0;
    while (size > 1n) {
      size >>= 1n;
      log++;
    }
    return log;
  }

  toString(): string {
    return `FFTDomain of size ${this.size} over ${this.field}`;
  }
}

/**
 * A coset domain for FRI and other ZK protocols.
 *
 * A coset is of the form {offset * omega^i | i = 0, 1, ..., n-1}
 * where omega is a primitive n-th root of unity.
 *
 * Cosets are disjoint from the original domain (if offset is not in the domain)
 * and are used in:
 * - FRI protocol for domain extension
 * - Low-degree testing
 * - STARK proof systems
 */
export class CosetDomain {
  readonly field: FiniteFieldLike;
  readonly size: bigint;
  readonly generator: FiniteFieldElementLike;
  readonly offset: FiniteFieldElementLike;
  private _elements: FiniteFieldElementLike[] | null = null;

  /**
   * Create a coset domain.
   *
   * @param field - The finite field
   * @param size - The size of the domain
   * @param generator - The primitive root of unity
   * @param offset - The coset offset
   */
  constructor(
    field: FiniteFieldLike,
    size: bigint | number,
    generator: FiniteFieldElementLike,
    offset: FiniteFieldElementLike
  ) {
    this.field = field;
    this.size = typeof size === 'number' ? BigInt(size) : size;
    this.generator = generator;
    this.offset = offset;
  }

  /**
   * Return all elements of the coset: [offset, offset*omega, ..., offset*omega^(n-1)].
   */
  elements(): FiniteFieldElementLike[] {
    if (this._elements !== null) {
      return this._elements;
    }

    this._elements = [];
    let current = this.offset;

    for (let i = 0n; i < this.size; i++) {
      this._elements.push(current);
      current = current.mul(this.generator);
    }

    return this._elements;
  }

  /**
   * Check if an element is in this coset.
   *
   * An element x is in the coset iff (x/offset)^n = 1,
   * i.e., x = offset * omega^k for some k.
   */
  contains(element: FiniteFieldElementLike): boolean {
    // We need to check if element / offset is a power of generator
    // This is expensive, so we check if element is in the computed list
    const elems = this.elements();
    for (const e of elems) {
      if (element.eq(e)) {
        return true;
      }
    }
    return false;
  }

  /**
   * Fold the coset domain by half.
   *
   * This is a key operation in the FRI protocol. Given a domain of size n
   * with generator omega, folding creates a domain of size n/2 with
   * generator omega^2 and offset adjusted by the folding challenge.
   *
   * @param challenge - The folding randomness (alpha in FRI)
   * @returns A new coset domain of half the size
   * @throws {ValueError} If the domain cannot be folded (size = 1)
   */
  fold(challenge: FiniteFieldElementLike): CosetDomain {
    if (this.size <= 1n) {
      throw new ValueError('cannot fold domain of size 1');
    }

    const newSize = this.size / 2n;
    const newGenerator = this.generator.mul(this.generator); // omega^2
    const newOffset = this.offset.mul(this.offset); // offset^2

    return new CosetDomain(this.field, newSize, newGenerator, newOffset);
  }

  /**
   * Return the vanishing polynomial of this coset.
   *
   * For a coset {offset * omega^i}, the vanishing polynomial is:
   * Z(x) = x^n - offset^n
   */
  vanishingPolynomial<C extends RingElement>(ring: PolynomialRing<C>): Polynomial<C> {
    const n = Number(this.size);

    // Compute offset^n
    const offsetN = this.offset.pow(this.size);

    // x^n - offset^n
    const coeffs: C[] = [];

    // Get the value as appropriate for the ring
    let negOffsetN: C;
    if ('value' in offsetN && typeof (offsetN as { value: bigint }).value === 'bigint') {
      negOffsetN = ring.base_ring.__call__(-(offsetN as { value: bigint }).value) as C;
    } else {
      // For extension fields, we need to handle this differently
      // Assume we can negate and convert
      negOffsetN = ring.base_ring.__call__(-1) as C;
      const offsetNCoerced = ring.base_ring.__call__(offsetN) as C;
      negOffsetN = negOffsetN.mul(offsetNCoerced) as C;
    }

    // Constant term: -offset^n
    coeffs.push(negOffsetN);

    // Zeros for x^1 through x^(n-1)
    for (let i = 1; i < n; i++) {
      coeffs.push(ring.base_ring.zero() as C);
    }

    // Leading coefficient: 1
    coeffs.push(ring.base_ring.one() as C);

    return new Polynomial(coeffs, ring);
  }

  /**
   * Get the size of the coset (for iteration).
   */
  get length(): number {
    return Number(this.size);
  }

  /**
   * Iterate over coset elements.
   */
  *[Symbol.iterator](): Iterator<FiniteFieldElementLike> {
    yield* this.elements();
  }

  toString(): string {
    return `CosetDomain of size ${this.size} with offset ${this.offset}`;
  }
}

/**
 * Find the largest power-of-2 FFT domain size supported by the field.
 *
 * This is useful for determining the maximum FFT size for a given field.
 *
 * @param field - The finite field
 * @returns The largest power of 2 that divides |F*|
 *
 * @example
 * ```typescript
 * maxFFTSize(GF(65537n))  // 65536n = 2^16
 * maxFFTSize(GF(17n))     // 16n = 2^4
 * ```
 */
export function maxFFTSize(field: FiniteFieldLike): bigint {
  const groupOrder = field.order - 1n;

  // Find the largest power of 2 dividing the group order
  let power = 1n;
  let remaining = groupOrder;

  while ((remaining & 1n) === 0n) {
    power <<= 1n;
    remaining >>= 1n;
  }

  return power;
}

/**
 * Find all power-of-2 sizes that are valid FFT domain sizes for the field.
 *
 * @param field - The finite field
 * @returns Array of valid FFT sizes [1, 2, 4, 8, ..., maxSize]
 */
export function validFFTSizes(field: FiniteFieldLike): bigint[] {
  const max = maxFFTSize(field);
  const sizes: bigint[] = [];

  let size = 1n;
  while (size <= max) {
    sizes.push(size);
    size <<= 1n;
  }

  return sizes;
}

/**
 * Compute the two-adicity of a finite field.
 *
 * The two-adicity is the largest k such that 2^k divides |F*|.
 * This determines the maximum FFT size supported by the field.
 *
 * @param field - The finite field
 * @returns The two-adicity (log base 2 of maxFFTSize)
 *
 * @example
 * ```typescript
 * twoAdicity(GF(65537n))  // 16 (since 65537 - 1 = 2^16)
 * twoAdicity(GF(17n))     // 4 (since 17 - 1 = 16 = 2^4)
 * ```
 */
export function twoAdicity(field: FiniteFieldLike): number {
  const groupOrder = field.order - 1n;

  let count = 0;
  let remaining = groupOrder;

  while ((remaining & 1n) === 0n) {
    count++;
    remaining >>= 1n;
  }

  return count;
}

/**
 * Create an FFT domain of the specified size, or the largest valid size <= requested.
 *
 * @param field - The finite field
 * @param size - The desired domain size (will be rounded down to power of 2 if needed)
 * @param exact - If true, throw error if exact size not possible (default: false)
 * @returns An FFT domain
 *
 * @example
 * ```typescript
 * const F = GF(17n);
 * const domain = createFFTDomain(F, 16n);  // Creates domain of size 16
 * const domain2 = createFFTDomain(F, 10n); // Creates domain of size 8 (largest power of 2 <= 10)
 * ```
 */
export function createFFTDomain(
  field: FiniteFieldLike,
  size: bigint | number,
  exact: boolean = false
): FFTDomain {
  let targetSize = typeof size === 'number' ? BigInt(size) : size;

  if (targetSize <= 0n) {
    throw new ValueError('domain size must be positive');
  }

  // Round down to nearest power of 2 if not exact
  if ((targetSize & (targetSize - 1n)) !== 0n) {
    if (exact) {
      throw new ValueError(`${targetSize} is not a power of 2`);
    }
    // Find largest power of 2 <= targetSize
    let pow2 = 1n;
    while (pow2 * 2n <= targetSize) {
      pow2 *= 2n;
    }
    targetSize = pow2;
  }

  // Check if this size is valid
  const maxSize = maxFFTSize(field);
  if (targetSize > maxSize) {
    if (exact) {
      throw new ValueError(
        `FFT domain of size ${targetSize} not supported. ` +
          `Maximum size for this field is ${maxSize}.`
      );
    }
    targetSize = maxSize;
  }

  return new FFTDomain(field, targetSize);
}

/**
 * Find an element with specific multiplicative order in the field.
 *
 * @param field - The finite field
 * @param order - The desired multiplicative order
 * @returns An element with the specified order
 * @throws {ValueError} If no such element exists
 */
export function elementOfOrder(
  field: FiniteFieldLike,
  order: bigint | number
): FiniteFieldElementLike {
  const targetOrder = typeof order === 'number' ? BigInt(order) : order;

  if (targetOrder <= 0n) {
    throw new ValueError('order must be positive');
  }

  const groupOrder = field.order - 1n;

  // Check that such an element can exist
  if (groupOrder % targetOrder !== 0n) {
    throw new ValueError(
      `No element of order ${targetOrder} exists. ` +
        `${targetOrder} does not divide ${groupOrder}.`
    );
  }

  // g^(|F*|/order) has the desired order
  const g = findMultiplicativeGenerator(field);
  const exp = groupOrder / targetOrder;
  const result = g.pow(exp);

  // Verify (for debugging)
  if (!result.pow(targetOrder).isOne()) {
    throw new ValueError('internal error: computed element does not have correct order');
  }

  return result;
}
