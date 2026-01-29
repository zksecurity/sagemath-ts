/**
 * @module sage/rings/finite_rings/gf2
 * @description The finite field GF(2) = {0, 1}
 *
 * Port of: sage/rings/finite_rings/finite_field_prime_modn.py (for p=2)
 */

import { ValueError, ZeroDivisionError } from '../../errors.js';

/**
 * Element of GF(2).
 */
export class GF2Element {
  readonly value: 0 | 1;
  readonly parent: GF2Field;

  constructor(value: number | bigint | boolean | GF2Element, parent: GF2Field) {
    this.parent = parent;

    if (value instanceof GF2Element) {
      this.value = value.value;
    } else if (typeof value === 'boolean') {
      this.value = value ? 1 : 0;
    } else {
      const v = typeof value === 'bigint' ? Number(value) : value;
      this.value = (((v % 2) + 2) % 2) as 0 | 1;
    }
  }

  add(other: GF2Element | number): GF2Element {
    const otherVal = other instanceof GF2Element ? other.value : ((other % 2) + 2) % 2;
    return new GF2Element((this.value + otherVal) % 2, this.parent);
  }

  sub(other: GF2Element | number): GF2Element {
    // In GF(2), subtraction is the same as addition
    return this.add(other);
  }

  mul(other: GF2Element | number): GF2Element {
    const otherVal = other instanceof GF2Element ? other.value : ((other % 2) + 2) % 2;
    return new GF2Element(((this.value * otherVal) % 2) as 0 | 1, this.parent);
  }

  div(other: GF2Element | number): GF2Element {
    const otherVal = other instanceof GF2Element ? other.value : ((other % 2) + 2) % 2;
    if (otherVal === 0) {
      throw new ZeroDivisionError('division by zero in GF(2)');
    }
    return new GF2Element(this.value, this.parent);
  }

  neg(): GF2Element {
    // In GF(2), -x = x
    return new GF2Element(this.value, this.parent);
  }

  inv(): GF2Element {
    if (this.value === 0) {
      throw new ZeroDivisionError('division by zero in GF(2)');
    }
    return new GF2Element(1, this.parent);
  }

  pow(n: number | bigint): GF2Element {
    if (this.value === 0) {
      if (n === 0 || n === 0n) {
        return new GF2Element(1, this.parent);
      }
      return new GF2Element(0, this.parent);
    }
    return new GF2Element(1, this.parent);
  }

  eq(other: GF2Element | number): boolean {
    const otherVal = other instanceof GF2Element ? other.value : ((other % 2) + 2) % 2;
    return this.value === otherVal;
  }

  isZero(): boolean {
    return this.value === 0;
  }

  isOne(): boolean {
    return this.value === 1;
  }

  toString(): string {
    return this.value.toString();
  }

  repr(): string {
    return this.value.toString();
  }

  // For use as polynomial coefficient
  toBigInt(): bigint {
    return BigInt(this.value);
  }
}

/**
 * The finite field GF(2).
 */
export class GF2Field {
  private static instance: GF2Field;

  readonly characteristic = 2n;
  readonly order = 2n;
  readonly degree = 1;

  private constructor() {}

  static getInstance(): GF2Field {
    if (!GF2Field.instance) {
      GF2Field.instance = new GF2Field();
    }
    return GF2Field.instance;
  }

  /**
   * Create an element of GF(2).
   */
  __call__(x: number | bigint | boolean | GF2Element): GF2Element {
    return new GF2Element(x, this);
  }

  /**
   * Return the zero element.
   */
  zero(): GF2Element {
    return new GF2Element(0, this);
  }

  /**
   * Return the one element.
   */
  one(): GF2Element {
    return new GF2Element(1, this);
  }

  /**
   * Return the generator (which is 1 for GF(2)).
   */
  gen(): GF2Element {
    return new GF2Element(1, this);
  }

  /**
   * Iterate over all elements.
   */
  *[Symbol.iterator](): Iterator<GF2Element> {
    yield this.zero();
    yield this.one();
  }

  /**
   * Return the number of elements.
   */
  cardinality(): bigint {
    return 2n;
  }

  is_field(): boolean {
    return true;
  }

  toString(): string {
    return 'Finite Field of size 2';
  }
}

/**
 * The finite field GF(2).
 */
export const GF2 = GF2Field.getInstance();
