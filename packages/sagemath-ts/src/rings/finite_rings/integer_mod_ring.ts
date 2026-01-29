/**
 * @module sage/rings/finite_rings/integer_mod_ring
 * @description The ring Z/nZ of integers modulo n
 *
 * Port of: sage/rings/finite_rings/integer_mod_ring.py
 */

import { gcd, is_prime } from '../../arith/misc.js';
import { ValueError } from '../../errors.js';
import { current_randstate } from '../../misc/randstate.js';
import type { CoefficientRing, RingElement } from '../polynomial/polynomial_element.js';
import { IntegerMod, type IntegerModRingBase } from './integer_mod.js';

/**
 * The ring Z/nZ of integers modulo n.
 *
 * This class implements the ring of integers modulo n, where n is a positive integer.
 * When n is prime, this ring is a field (but for fields, use GF(p) for additional
 * field-specific functionality).
 *
 * @example
 * ```typescript
 * const Z7 = Zmod(7n);
 * const a = Z7(3n);
 * const b = Z7(4n);
 * console.log(a.mul(b));  // 5 (since 3*4 = 12 ≡ 5 mod 7)
 *
 * // Iterate over elements
 * for (const x of Z7) {
 *   console.log(x.value);  // 0, 1, 2, 3, 4, 5, 6
 * }
 * ```
 */
export class IntegerModRing implements IntegerModRingBase, CoefficientRing<IntegerMod> {
  readonly modulus: bigint;
  readonly characteristic: bigint;
  readonly order: bigint;

  private _zero: IntegerMod | null = null;
  private _one: IntegerMod | null = null;

  /**
   * Create the ring Z/nZ.
   *
   * @param n - The modulus (must be a positive integer)
   */
  constructor(n: bigint | number) {
    const modulus = typeof n === 'number' ? BigInt(n) : n;

    if (modulus <= 0n) {
      throw new ValueError('modulus must be positive');
    }

    this.modulus = modulus;
    this.characteristic = modulus;
    this.order = modulus;
  }

  /**
   * Create an element of this ring.
   *
   * @param x - The value to convert to an element
   */
  __call__(x: unknown): IntegerMod {
    if (typeof x === 'number') {
      return new IntegerMod(BigInt(x), this);
    }
    if (typeof x === 'bigint') {
      return new IntegerMod(x, this);
    }
    if (x instanceof IntegerMod) {
      return new IntegerMod(x.value, this);
    }
    if (typeof x === 'boolean') {
      return new IntegerMod(x ? 1n : 0n, this);
    }
    throw new ValueError(`cannot coerce ${x} to Z/${this.modulus}Z`);
  }

  /**
   * Return the zero element.
   */
  zero(): IntegerMod {
    if (this._zero === null) {
      this._zero = new IntegerMod(0n, this);
    }
    return this._zero;
  }

  /**
   * Return the one element.
   */
  one(): IntegerMod {
    if (this._one === null) {
      this._one = new IntegerMod(1n, this);
    }
    return this._one;
  }

  /**
   * Return a generator of the ring (which is just 1).
   */
  gen(): IntegerMod {
    return this.one();
  }

  /**
   * Check if this ring is a field.
   * Z/nZ is a field if and only if n is prime.
   */
  is_field(): boolean {
    return is_prime(this.modulus);
  }

  /**
   * Return the cardinality of this ring.
   */
  cardinality(): bigint {
    return this.modulus;
  }

  /**
   * Iterate over all elements of this ring.
   *
   * @example
   * ```typescript
   * const Z5 = Zmod(5n);
   * for (const x of Z5) {
   *   console.log(x.value);  // 0, 1, 2, 3, 4
   * }
   * ```
   */
  *[Symbol.iterator](): Iterator<IntegerMod> {
    for (let i = 0n; i < this.modulus; i++) {
      yield new IntegerMod(i, this);
    }
  }

  /**
   * Return a list of all elements.
   */
  list(): IntegerMod[] {
    return Array.from(this);
  }

  /**
   * Return a list of all units (invertible elements).
   */
  units(): IntegerMod[] {
    const result: IntegerMod[] = [];
    for (let i = 1n; i < this.modulus; i++) {
      if (gcd(i, this.modulus) === 1n) {
        result.push(new IntegerMod(i, this));
      }
    }
    return result;
  }

  /**
   * Return a random element of this ring.
   */
  random_element(): IntegerMod {
    const rstate = current_randstate();
    const randomInt = rstate.random_below(this.modulus);
    return new IntegerMod(randomInt, this);
  }

  /**
   * String representation.
   */
  toString(): string {
    return `Ring of integers modulo ${this.modulus}`;
  }
}

/**
 * Factory function to create Z/nZ.
 *
 * @param n - The modulus (must be a positive integer)
 * @returns The ring Z/nZ
 *
 * @example
 * ```typescript
 * const Z12 = Zmod(12n);
 * const a = Z12(7n);
 * const b = Z12(8n);
 * console.log(a.add(b));  // 3 (since 7+8 = 15 ≡ 3 mod 12)
 * ```
 */
export function Zmod(n: bigint | number): IntegerModRing {
  return new IntegerModRing(n);
}

/**
 * Alias for Zmod - IntegerModRing factory.
 */
export const IntegerModRingFactory = Zmod;

/**
 * Alias for Zmod - Integers function.
 */
export const Integers = Zmod;
