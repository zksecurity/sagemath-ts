/**
 * @module sage/rings/finite_rings/gf2
 * @description The finite field GF(2) = {0, 1}
 *
 * Port of: sage/rings/finite_rings/finite_field_prime_modn.py (for p=2)
 */

import { ValueError, ZeroDivisionError } from '../../errors.js';
import { Mod } from './integer_mod.js';
import { Rational } from '../rational.js';
import { PrimeField } from './finite_field_extension.js';

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

  _integer_(_ZZ?: unknown): bigint {
    return BigInt(this.value);
  }
  _rational_(): Rational {
    return new Rational(BigInt(this.value));
  }

  add(other: GF2Element | number | bigint | boolean): GF2Element {
    return this.parent.__call__((this.value + toBit(other)) % 2);
  }

  sub(other: GF2Element | number | bigint | boolean): GF2Element {
    // In GF(2), subtraction is the same as addition
    return this.add(other);
  }

  mul(other: GF2Element | number | bigint | boolean): GF2Element {
    return this.parent.__call__((this.value * toBit(other)) % 2);
  }

  div(other: GF2Element | number | bigint | boolean): GF2Element {
    if (toBit(other) === 0) {
      throw new ZeroDivisionError('inverse of Mod(0, 2) does not exist');
    }
    return this.parent.__call__(this.value);
  }

  neg(): GF2Element {
    // IntegerMod_int preserves an uncached zero; nonzero results use its table.
    return this.value === 0 ? this : this.parent.__call__(this.value);
  }

  inv(): GF2Element {
    if (this.value === 0) {
      throw new ZeroDivisionError('inverse of Mod(0, 2) does not exist');
    }
    return this.parent.one();
  }

  /** Sage GF(2) uses IntegerMod exponent conversion and native/GMP error paths. */
  pow(n: number | bigint): GF2Element {
    return this.parent.__call__(Mod(BigInt(this.value), 2n).pow(n).value);
  }

  /** Both residues are squares, as in IntegerMod_int.is_square modulo two. */
  is_square(): boolean { return true; }
  /** Native tiny-modulus roots return the cached coefficient element.
   * @see Deviation: Binary field coefficient roots and caching
   */
  sqrt(options: {all: true; extend?: boolean}): GF2Element[];
  sqrt(options?: {all?: false; extend?: boolean}): GF2Element;
  sqrt(options: {all?: boolean; extend?: boolean}): GF2Element | GF2Element[];
  sqrt(options: {all?: boolean; extend?: boolean} = {}): GF2Element | GF2Element[] {
    for (const key of Object.keys(options))
      if (key !== 'all' && key !== 'extend') throw new TypeError(`sqrt() got an unexpected keyword argument '${key}'`);
    const root = this.parent.__call__(this.value);
    return options.all ? [root] : root;
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

  // Native IntegerMod_int uses a two-entry modulus table for GF(2).
  private readonly elements: [GF2Element, GF2Element] = [new GF2Element(0, this), new GF2Element(1, this)];
  private constructor() {}

  static getInstance(): GF2Field {
    if (!GF2Field.instance) {
      GF2Field.instance = new GF2Field();
    }
    return GF2Field.instance;
  }

  /**
   * Create an element of GF(2).
   * Nontrivial scalar conversions use the shared prime-field constructor,
   * matching IntegerMod.__init__ (including strings and None/undefined).
   */
  __call__(x?: unknown): GF2Element {
    if (x instanceof GF2Element) return x;
    if (
      typeof x === 'bigint' || typeof x === 'boolean' ||
      (typeof x === 'number' && Number.isInteger(x))
    ) return this.elements[toBit(x)];
    return this.elements[Number(new PrimeField(2n).__call__(x).value)]!;
  }

  /**
   * Return the zero element.
   */
  zero(): GF2Element {
    return this.elements[0];
  }

  /**
   * Return the one element.
   */
  one(): GF2Element {
    return this.elements[1];
  }

  /**
   * Return the generator (which is 1 for GF(2)).
   */
  gen(): GF2Element {
    return this.elements[1];
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
