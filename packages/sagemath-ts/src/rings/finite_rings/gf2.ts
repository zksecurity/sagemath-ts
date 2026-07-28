/**
 * @module sage/rings/finite_rings/gf2
 * @description The finite field GF(2) = {0, 1}
 *
 * Port of: sage/rings/finite_rings/finite_field_prime_modn.py (for p=2)
 */

import { ValueError, ZeroDivisionError } from '../../errors.js';

/**
 * Reduce an integer input modulo 2.
 *
 * The reduction is done in `bigint` arithmetic: converting to `number` first
 * silently rounds every input above 2^53, so e.g. `GF(2)(2^64 + 1)` came out
 * as 0 instead of 1.
 */
function toBit(x: number | bigint | boolean | GF2Element): 0 | 1 {
  if (x instanceof GF2Element) {
    return x.value;
  }
  if (typeof x === 'boolean') {
    return x ? 1 : 0;
  }
  let v: bigint;
  if (typeof x === 'bigint') {
    v = x;
  } else {
    if (!Number.isInteger(x)) {
      throw new ValueError(`unable to convert ${x} to an integer`);
    }
    v = BigInt(x);
  }
  return Number(((v % 2n) + 2n) % 2n) as 0 | 1;
}

/**
 * Element of GF(2).
 */
export class GF2Element {
  readonly value: 0 | 1;
  readonly parent: GF2Field;

  constructor(value: number | bigint | boolean | GF2Element, parent: GF2Field) {
    this.parent = parent;
    this.value = toBit(value);
  }

  add(other: GF2Element | number | bigint | boolean): GF2Element {
    return new GF2Element(((this.value + toBit(other)) % 2) as 0 | 1, this.parent);
  }

  sub(other: GF2Element | number | bigint | boolean): GF2Element {
    // In GF(2), subtraction is the same as addition
    return this.add(other);
  }

  mul(other: GF2Element | number | bigint | boolean): GF2Element {
    return new GF2Element(((this.value * toBit(other)) % 2) as 0 | 1, this.parent);
  }

  div(other: GF2Element | number | bigint | boolean): GF2Element {
    if (toBit(other) === 0) {
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
    const e = typeof n === 'bigint' ? n : BigInt(n);
    if (this.value === 0) {
      if (e === 0n) {
        return new GF2Element(1, this.parent);
      }
      if (e < 0n) {
        // Sage: GF(2)(0)^-1 raises ZeroDivisionError
        throw new ZeroDivisionError('division by zero in GF(2)');
      }
      return new GF2Element(0, this.parent);
    }
    return new GF2Element(1, this.parent);
  }

  eq(other: GF2Element | number | bigint | boolean): boolean {
    return this.value === toBit(other);
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
