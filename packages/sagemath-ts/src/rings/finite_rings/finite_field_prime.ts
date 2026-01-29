/**
 * @module sage/rings/finite_rings/finite_field_prime
 * @description Prime finite fields GF(p)
 *
 * Port of: sage/rings/finite_rings/finite_field_prime_modn.py
 */

import { factor, gcd, is_prime, power_mod, xgcd } from '../../arith/misc.js';
import { ArithmeticError, ValueError, ZeroDivisionError } from '../../errors.js';
import { current_randstate } from '../../misc/randstate.js';
import type { CoefficientRing, RingElement } from '../polynomial/polynomial_element.js';

// Import PARI functions for finite field operations
import { Fp_sqrt } from '@sagemath-ts/parigp-ts';

/**
 * An element of a prime field GF(p).
 *
 * Elements are represented as integers in the range [0, p).
 * All non-zero elements are invertible since p is prime.
 *
 * @example
 * ```typescript
 * const F7 = GF(7n);
 * const a = F7(3n);
 * const b = F7(5n);
 * console.log(a.mul(b));  // 1 (since 3*5 = 15 ≡ 1 mod 7)
 * console.log(a.inv());   // 5 (since 3*5 ≡ 1 mod 7)
 * ```
 */
export class FiniteFieldElement implements RingElement {
  readonly value: bigint;
  readonly parent: FiniteFieldPrime;

  /**
   * Create an element of GF(p).
   *
   * @param value - The integer value (will be reduced modulo p)
   * @param parent - The parent field GF(p)
   */
  constructor(value: bigint | number | FiniteFieldElement, parent: FiniteFieldPrime) {
    this.parent = parent;

    if (value instanceof FiniteFieldElement) {
      this.value = mod(value.value, parent.characteristic);
    } else {
      const v = typeof value === 'number' ? BigInt(value) : value;
      this.value = mod(v, parent.characteristic);
    }
  }

  /**
   * Return the characteristic (prime modulus) of the field.
   */
  get p(): bigint {
    return this.parent.characteristic;
  }

  /**
   * Add two elements.
   */
  add(other: FiniteFieldElement | number | bigint): FiniteFieldElement {
    const otherVal = this.coerceValue(other);
    return new FiniteFieldElement(mod(this.value + otherVal, this.p), this.parent);
  }

  /**
   * Subtract two elements.
   */
  sub(other: FiniteFieldElement | number | bigint): FiniteFieldElement {
    const otherVal = this.coerceValue(other);
    return new FiniteFieldElement(mod(this.value - otherVal, this.p), this.parent);
  }

  /**
   * Multiply two elements.
   */
  mul(other: FiniteFieldElement | number | bigint): FiniteFieldElement {
    const otherVal = this.coerceValue(other);
    return new FiniteFieldElement(mod(this.value * otherVal, this.p), this.parent);
  }

  /**
   * Divide two elements.
   *
   * @throws {ZeroDivisionError} If dividing by zero
   */
  div(other: FiniteFieldElement | number | bigint): FiniteFieldElement {
    const otherVal = this.coerceValue(other);
    if (otherVal === 0n) {
      throw new ZeroDivisionError('division by zero in GF(p)');
    }

    // In a prime field, we can always compute inverse via extended Euclidean algorithm
    const [_g, s] = xgcd(otherVal, this.p);
    // g is always 1 since p is prime and otherVal != 0

    return new FiniteFieldElement(mod(this.value * s, this.p), this.parent);
  }

  /**
   * Return the additive inverse (-self).
   */
  neg(): FiniteFieldElement {
    if (this.value === 0n) {
      return this;
    }
    return new FiniteFieldElement(this.p - this.value, this.parent);
  }

  /**
   * Return the multiplicative inverse (1/self).
   *
   * @throws {ZeroDivisionError} If self is zero
   */
  inv(): FiniteFieldElement {
    if (this.value === 0n) {
      throw new ZeroDivisionError('division by zero in GF(p)');
    }

    // Extended Euclidean algorithm
    const [_g, s] = xgcd(this.value, this.p);
    // g is always 1 since p is prime and value != 0

    return new FiniteFieldElement(mod(s, this.p), this.parent);
  }

  /**
   * Return self^n.
   *
   * @param n - The exponent (can be negative)
   */
  pow(n: number | bigint): FiniteFieldElement {
    const exp = typeof n === 'number' ? BigInt(n) : n;

    if (exp === 0n) {
      return this.parent.one();
    }

    if (exp < 0n) {
      // For negative exponents, compute inverse first
      return this.inv().pow(-exp);
    }

    // Use fast modular exponentiation
    const result = power_mod(this.value, exp, this.p);
    return new FiniteFieldElement(result, this.parent);
  }

  /**
   * Check equality with another element.
   */
  eq(other: FiniteFieldElement | number | bigint): boolean {
    const otherVal = this.coerceValue(other);
    return this.value === otherVal;
  }

  /**
   * Check if this element is zero.
   */
  isZero(): boolean {
    return this.value === 0n;
  }

  /**
   * Check if this element is one.
   */
  isOne(): boolean {
    return this.value === 1n;
  }

  /**
   * Check if this element is a unit (always true for non-zero elements in a field).
   */
  isUnit(): boolean {
    return this.value !== 0n;
  }

  /**
   * Lift this element to an integer.
   */
  lift(): bigint {
    return this.value;
  }

  /**
   * Return the integer value as a bigint.
   */
  toBigInt(): bigint {
    return this.value;
  }

  /**
   * Return the multiplicative order of this element.
   *
   * For a non-zero element in GF(p), the order divides p-1.
   *
   * @throws {ValueError} If the element is zero
   */
  multiplicative_order(): bigint {
    if (this.value === 0n) {
      throw new ValueError('multiplicative order of zero is not defined');
    }

    if (this.value === 1n) {
      return 1n;
    }

    // The order divides p-1
    const pMinus1 = this.p - 1n;

    // Try divisors of p-1
    const factors = primeFactorsSimple(pMinus1);
    let order = pMinus1;

    for (const q of factors) {
      while (order % q === 0n) {
        const testOrder = order / q;
        if (this.pow(testOrder).isOne()) {
          order = testOrder;
        } else {
          break;
        }
      }
    }

    return order;
  }

  /**
   * Check if this element is a quadratic residue (has a square root).
   *
   * For p > 2, this uses the Legendre symbol.
   */
  is_square(): boolean {
    if (this.value === 0n) {
      return true;
    }

    if (this.p === 2n) {
      return true; // Every element of GF(2) is a square
    }

    // Euler's criterion: a is a square iff a^((p-1)/2) = 1
    return this.pow((this.p - 1n) / 2n).isOne();
  }

  /**
   * Compute a square root of this element.
   *
   * @throws {ValueError} If the element is not a quadratic residue
   */
  sqrt(): FiniteFieldElement {
    if (this.value === 0n) {
      return this;
    }

    if (!this.is_square()) {
      throw new ValueError(`${this.value} is not a square in GF(${this.p})`);
    }

    if (this.p === 2n) {
      return this;
    }

    // Use PARI's Fp_sqrt (Tonelli-Shanks algorithm)
    return sqrtMod(this.value, this.p, this.parent);
  }

  /**
   * String representation.
   */
  toString(): string {
    return this.value.toString();
  }

  /**
   * Repr for debugging.
   */
  repr(): string {
    return this.value.toString();
  }

  /**
   * Coerce a value to a bigint in [0, p).
   */
  private coerceValue(other: FiniteFieldElement | number | bigint): bigint {
    if (other instanceof FiniteFieldElement) {
      return mod(other.value, this.p);
    }
    const v = typeof other === 'number' ? BigInt(other) : other;
    return mod(v, this.p);
  }
}

/**
 * A prime finite field GF(p).
 *
 * This is the field of integers modulo a prime p.
 *
 * @example
 * ```typescript
 * const F7 = new FiniteFieldPrime(7n);
 * const a = F7(3n);
 * const b = F7(5n);
 *
 * console.log(a.mul(b));        // 1
 * console.log(a.inv());         // 5
 * console.log(F7.cardinality()); // 7n
 *
 * // Iterate over all elements
 * for (const x of F7) {
 *   console.log(x.value);
 * }
 * ```
 */
export class FiniteFieldPrime implements CoefficientRing<FiniteFieldElement> {
  readonly characteristic: bigint;
  readonly order: bigint;
  readonly degree: number = 1;

  private _zero: FiniteFieldElement | null = null;
  private _one: FiniteFieldElement | null = null;

  /**
   * Create a prime finite field GF(p).
   *
   * @param p - The prime modulus
   * @param check - Whether to verify that p is prime (default: true)
   * @throws {ArithmeticError} If p is not prime and check is true
   */
  constructor(p: bigint | number, check: boolean = true) {
    const prime = typeof p === 'number' ? BigInt(p) : p;

    if (prime < 2n) {
      throw new ArithmeticError('p must be at least 2');
    }

    if (check && !is_prime(prime)) {
      throw new ArithmeticError('p must be prime');
    }

    this.characteristic = prime;
    this.order = prime;
  }

  /**
   * Create an element of this field.
   *
   * @param x - The value to convert to an element
   */
  __call__(x: unknown): FiniteFieldElement {
    if (typeof x === 'number') {
      return new FiniteFieldElement(BigInt(x), this);
    }
    if (typeof x === 'bigint') {
      return new FiniteFieldElement(x, this);
    }
    if (x instanceof FiniteFieldElement) {
      return new FiniteFieldElement(x.value, this);
    }
    if (typeof x === 'boolean') {
      return new FiniteFieldElement(x ? 1n : 0n, this);
    }
    throw new ValueError(`cannot coerce ${x} to GF(${this.characteristic})`);
  }

  /**
   * Return the zero element.
   */
  zero(): FiniteFieldElement {
    if (this._zero === null) {
      this._zero = new FiniteFieldElement(0n, this);
    }
    return this._zero;
  }

  /**
   * Return the one element.
   */
  one(): FiniteFieldElement {
    if (this._one === null) {
      this._one = new FiniteFieldElement(1n, this);
    }
    return this._one;
  }

  /**
   * Return a generator of the field.
   *
   * For a prime field, this returns 1 (the additive generator).
   * For a multiplicative generator, use multiplicative_generator().
   */
  gen(): FiniteFieldElement {
    return this.one();
  }

  /**
   * Return a multiplicative generator (primitive root) of GF(p)*.
   *
   * This finds the smallest positive integer g such that g generates
   * the multiplicative group of the field.
   */
  multiplicative_generator(): FiniteFieldElement {
    if (this.characteristic === 2n) {
      return this.one();
    }

    const pMinus1 = this.characteristic - 1n;
    const factors = primeFactorsSimple(pMinus1);

    // Try candidates starting from 2
    for (let g = 2n; g < this.characteristic; g++) {
      let isGenerator = true;

      for (const q of factors) {
        const exp = pMinus1 / q;
        if (power_mod(g, exp, this.characteristic) === 1n) {
          isGenerator = false;
          break;
        }
      }

      if (isGenerator) {
        return new FiniteFieldElement(g, this);
      }
    }

    // Should never reach here for valid prime
    throw new ArithmeticError('no multiplicative generator found');
  }

  /**
   * This is a field.
   */
  is_field(): boolean {
    return true;
  }

  /**
   * Return the cardinality of this field.
   */
  cardinality(): bigint {
    return this.characteristic;
  }

  /**
   * Iterate over all elements of this field.
   */
  *[Symbol.iterator](): Iterator<FiniteFieldElement> {
    for (let i = 0n; i < this.characteristic; i++) {
      yield new FiniteFieldElement(i, this);
    }
  }

  /**
   * Return a list of all elements.
   */
  list(): FiniteFieldElement[] {
    return Array.from(this);
  }

  /**
   * Iterate over all elements of this field.
   * Alias for [Symbol.iterator] for compatibility.
   */
  elements(): IterableIterator<FiniteFieldElement> {
    return this[Symbol.iterator]();
  }

  /**
   * Return a random element of this field.
   */
  random_element(): FiniteFieldElement {
    const rstate = current_randstate();
    const randomInt = rstate.random_below(this.characteristic);
    return new FiniteFieldElement(randomInt, this);
  }

  /**
   * Return a non-residue (an element that is not a quadratic residue).
   *
   * Useful for constructing quadratic extensions.
   */
  quadratic_non_residue(): FiniteFieldElement {
    if (this.characteristic === 2n) {
      throw new ValueError('no quadratic non-residue in GF(2)');
    }

    // Find the first non-residue
    for (let a = 2n; a < this.characteristic; a++) {
      const elem = new FiniteFieldElement(a, this);
      if (!elem.is_square()) {
        return elem;
      }
    }

    throw new ArithmeticError('no quadratic non-residue found');
  }

  /**
   * String representation.
   */
  toString(): string {
    return `Finite Field of size ${this.characteristic}`;
  }
}

/**
 * Compute a mod n, ensuring the result is in [0, n).
 */
function mod(a: bigint, n: bigint): bigint {
  const result = a % n;
  return result < 0n ? result + n : result;
}

/**
 * Get the unique prime factors of n (without multiplicity).
 *
 * Uses PARI's factorization via the factor() function.
 */
function primeFactorsSimple(n: bigint): bigint[] {
  if (n <= 1n) {
    return [];
  }

  // Use PARI's factorization and extract just the primes
  const factorization = factor(n);
  return factorization
    .filter(([p, _]) => p > 0n) // Exclude -1 sign factor
    .map(([p, _]) => p);
}

/**
 * Compute square root modulo a prime using PARI's Fp_sqrt.
 *
 * Uses PARI's Tonelli-Shanks implementation via @sagemath-ts/parigp-ts.
 *
 * @see Implementation: parigp-ts/src/ff.ts (port of pari/src/basemath/arith1.c)
 */
function sqrtMod(a: bigint, p: bigint, parent: FiniteFieldPrime): FiniteFieldElement {
  if (a === 0n) {
    return new FiniteFieldElement(0n, parent);
  }

  // Use PARI's Fp_sqrt (Tonelli-Shanks algorithm)
  const result = Fp_sqrt(a, p);

  if (result === null) {
    // a is not a quadratic residue - this shouldn't happen if called correctly
    throw new ArithmeticError(`${a} is not a quadratic residue mod ${p}`);
  }

  return new FiniteFieldElement(result, parent);
}
