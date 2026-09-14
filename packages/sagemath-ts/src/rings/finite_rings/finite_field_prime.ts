import { isFractionElement } from '../fraction_field_element.js';
/**
 * @module sage/rings/finite_rings/finite_field_prime
 * @description Prime finite fields GF(p)
 * @see Deviation: Finite Field Coercion and Backend Boundaries
 *
 * Port of: sage/rings/finite_rings/finite_field_prime_modn.py
 */

import { factor, gcd, inverse_mod, is_prime, power_mod, xgcd } from '../../arith/misc.js';
import { ArithmeticError, ValueError, ZeroDivisionError } from '../../errors.js';
import { current_randstate } from '../../misc/randstate.js';
import type { IntegerLike } from '../../types/coercion.js';
import { type Integer, ZZ } from '../integer_ring.js';
import type { CoefficientRing, RingElement } from '../polynomial/polynomial_element.js';
import { Rational } from '../rational.js';
import {
  FiniteFieldElement as ExtensionElement,
  PrimeField,
  PrimeFieldElement,
} from './finite_field_extension.js';
import {
  type IntegerMod,
  canonicalFiniteOperands,
  checkFiniteGeneratorIndex,
  finiteArithmeticEquals,
  repeatFiniteSequence,
} from './integer_mod.js';
import { IntegerModRing } from './integer_mod_ring.js';

// Import PARI functions for finite field operations
import { Fp_order } from '@sagemath-ts/parigp-ts';

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
  /** IntegerMod._rational_: lift the canonical residue into QQ. */
  _rational_(): Rational {
    return new Rational(this.value);
  }

  readonly value: bigint;
  readonly parent: FiniteFieldPrime;

  /**
   * Create an element of GF(p).
   *
   * @param value - The integer value (will be reduced modulo p)
   * @param parent - The parent field GF(p)
   */
  constructor(value: unknown, parent: FiniteFieldPrime) {
    this.parent = parent;
    const v = typeof value === 'bigint' ? value : parent.__call__(value).value;
    this.value = mod(v, parent.characteristic);
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
  add(other: ExtensionElement): ExtensionElement;
  add(
    other: FiniteFieldElement | PrimeFieldElement | IntegerMod | Integer | boolean | number | bigint
  ): FiniteFieldElement;
  add(
    other:
      | FiniteFieldElement
      | PrimeFieldElement
      | ExtensionElement
      | IntegerMod
      | Integer
      | boolean
      | number
      | bigint
  ): FiniteFieldElement | ExtensionElement;
  add(other: FiniteFieldElement): FiniteFieldElement;
  add(other: unknown): FiniteFieldElement | ExtensionElement {
    const [left, right] = canonicalFiniteOperands(this, other, '+');
    if (left instanceof ExtensionElement) return left.add(right as ExtensionElement);
    const operand = right as FiniteFieldElement;
    return new FiniteFieldElement(this.value + operand.value, this.parent);
  }

  /**
   * Subtract two elements.
   */
  sub(other: ExtensionElement): ExtensionElement;
  sub(
    other: FiniteFieldElement | PrimeFieldElement | IntegerMod | Integer | boolean | number | bigint
  ): FiniteFieldElement;
  sub(
    other:
      | FiniteFieldElement
      | PrimeFieldElement
      | ExtensionElement
      | IntegerMod
      | Integer
      | boolean
      | number
      | bigint
  ): FiniteFieldElement | ExtensionElement;
  sub(other: FiniteFieldElement): FiniteFieldElement;
  sub(other: unknown): FiniteFieldElement | ExtensionElement {
    const [left, right] = canonicalFiniteOperands(this, other, '-');
    if (left instanceof ExtensionElement) return left.sub(right as ExtensionElement);
    const operand = right as FiniteFieldElement;
    return new FiniteFieldElement(this.value - operand.value, this.parent);
  }

  /**
   * Multiply two elements.
   */
  mul(other: string): string;
  mul<T>(other: readonly T[]): T[];
  mul(other: ExtensionElement): ExtensionElement;
  mul(
    other: FiniteFieldElement | PrimeFieldElement | IntegerMod | Integer | boolean | number | bigint
  ): FiniteFieldElement;
  mul(
    other:
      | FiniteFieldElement
      | PrimeFieldElement
      | ExtensionElement
      | IntegerMod
      | Integer
      | boolean
      | number
      | bigint
  ): FiniteFieldElement | ExtensionElement;
  mul(other: FiniteFieldElement): FiniteFieldElement;
  mul(other: unknown): FiniteFieldElement | ExtensionElement | string | unknown[] {
    if (typeof other === 'string' || Array.isArray(other))
      return repeatFiniteSequence(this.value, other);
    const [left, right] = canonicalFiniteOperands(this, other, '*');
    if (left instanceof ExtensionElement) return left.mul(right as ExtensionElement);
    const operand = right as FiniteFieldElement;
    return new FiniteFieldElement(this.value * operand.value, this.parent);
  }

  /**
   * Divide two elements.
   *
   * @throws {ZeroDivisionError} If dividing by zero
   */
  div(other: ExtensionElement): ExtensionElement;
  div(
    other: FiniteFieldElement | PrimeFieldElement | IntegerMod | Integer | boolean | number | bigint
  ): FiniteFieldElement;
  div(
    other:
      | FiniteFieldElement
      | PrimeFieldElement
      | ExtensionElement
      | IntegerMod
      | Integer
      | boolean
      | number
      | bigint
  ): FiniteFieldElement | ExtensionElement;
  div(other: FiniteFieldElement): FiniteFieldElement;
  div(other: unknown): FiniteFieldElement | ExtensionElement {
    const [left, right] = canonicalFiniteOperands(this, other, '/');
    if (left instanceof ExtensionElement) return left.div(right as ExtensionElement);
    const operand = right as FiniteFieldElement;
    return this.mul(operand.inv());
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
      throw new ZeroDivisionError(`inverse of Mod(0, ${this.p}) does not exist`);
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
  pow(n: IntegerLike | number | Rational | boolean | string | null): FiniteFieldElement {
    // Sage prime-field elements inherit IntegerMod's native/GMP power dispatch.
    const result = new IntegerModRing(this.p).__call__(this.value).pow(n);
    return new FiniteFieldElement(result.value, this.parent);
  }

  /**
   * Check equality with another element.
   */
  eq(other: unknown): boolean {
    return finiteArithmeticEquals(this, other);
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
      throw new ArithmeticError(
        `multiplicative order of 0 not defined since it is not a unit modulo ${this.p}`
      );
    }
    // Sage IntegerMod.multiplicative_order -> PARI znorder -> Fp_order
    // (integer_mod.pyx:1896, arith1.c:2652).
    return Fp_order(this.value, this.p - 1n, this.p);
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
   * Compute one or all square roots using Sage's prime-field options/defaults.
   * Preserve this implementation's parent for roots in the base field.
   * @see Deviation: Finite Field Coercion and Backend Boundaries
   */
  sqrt(options: { all: true; extend?: boolean }): FiniteFieldElement[];
  sqrt(options: { extend: false; all?: false }): FiniteFieldElement;
  sqrt(options?: { extend?: boolean; all?: false }): FiniteFieldElement | ExtensionElement;
  sqrt(options: { extend?: boolean; all?: boolean }):
    | FiniteFieldElement
    | ExtensionElement
    | FiniteFieldElement[];
  sqrt(
    options: { extend?: boolean; all?: boolean } = {}
  ): FiniteFieldElement | ExtensionElement | FiniteFieldElement[] {
    const result = new PrimeField(this.p).__call__(this.value).sqrt(options);
    if (Array.isArray(result)) return result.map((x) => this.parent.__call__(x.value));
    return result instanceof PrimeFieldElement ? this.parent.__call__(result.value) : result;
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
  constructor(p: unknown, check: boolean = true) {
    const prime = ZZ.__call__(p as Parameters<typeof ZZ.__call__>[0]);
    if (check && !is_prime(prime)) throw new ArithmeticError('p must be prime');
    // Sage checks primality before constructing the generic positive-order ring.
    if (prime <= 0n) throw new ZeroDivisionError('order must be positive');

    this.characteristic = prime;
    this.order = prime;
  }

  /**
   * Create an element of this field.
   *
   * @param x - The value to convert to an element
   */
  /** IntegerMod.__init__ (integer_mod.pyx:377-404). */
  __call__(x?: unknown): FiniteFieldElement {
    if (
      isFractionElement(x) ||
      (typeof x === 'object' && x !== null && 'coeffs' in x && 'getCoeff' in x)
    )
      return new FiniteFieldElement(
        new IntegerModRing(this.characteristic).__call__(x).value,
        this
      );
    if (x instanceof FiniteFieldElement) return new FiniteFieldElement(x.value, this);
    if (x instanceof ExtensionElement)
      return new FiniteFieldElement(new PrimeField(this.characteristic).__call__(x).value, this);
    if (x instanceof Rational) {
      const d = mod(x.denominator, this.characteristic);
      if (d === 0n)
        throw new ZeroDivisionError(`inverse of Mod(0, ${this.characteristic}) does not exist`);
      return new FiniteFieldElement(x.numerator * inverse_mod(d, this.characteristic), this);
    }
    return new FiniteFieldElement(ZZ.__call__(x as Parameters<typeof ZZ.__call__>[0]), this);
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
  gen(n: unknown = 0): FiniteFieldElement {
    checkFiniteGeneratorIndex(n, true);
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
    const randomInt = rstate.python_random().randrange(this.characteristic);
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
