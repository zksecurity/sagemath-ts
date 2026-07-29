/**
 * @module sage/rings/function_field/place
 * @description Places of function fields
 *
 * The places of a function field correspond, one-to-one, to valuation rings of
 * the function field, each of which defines a discrete valuation for the
 * elements of the function field.  "Finite" places are in one-to-one
 * correspondence with the prime ideals of the finite maximal order while places
 * "at infinity" are in one-to-one correspondence with the prime ideals of the
 * infinite maximal order.
 *
 * Port of: sage/rings/function_field/place.py
 */

import { ValueError } from '../../errors.js';
import type { ConstantField, ConstantFieldElement } from './constant_field.js';
import type { FunctionFieldDivisor } from './divisor.js';
import { prime_divisor } from './divisor.js';
import type { FunctionFieldElement } from './element.js';
import type { FunctionField } from './function_field.js';
import { FunctionFieldIdeal } from './ideal.js';
import { FunctionFieldOrderInfinite } from './order.js';
import type { FunctionFieldValuationRing } from './valuation_ring.js';

/**
 * Places of function fields.
 *
 * @see Reference: sage/rings/function_field/place.py:67 (FunctionFieldPlace)
 */
export abstract class FunctionFieldPlace<C extends ConstantFieldElement> {
  readonly _parent: PlaceSet<C>;
  readonly _prime: FunctionFieldIdeal<C>;

  constructor(parent: PlaceSet<C>, prime: FunctionFieldIdeal<C>) {
    this._parent = parent;
    this._prime = prime;
  }

  parent(): PlaceSet<C> {
    return this._parent;
  }

  /**
   * @see Reference: sage/rings/function_field/place.py:113 (_repr_)
   */
  _repr_(): string {
    const gens = this._prime.gens();
    return `Place (${gens.map((g) => g.toString()).join(', ')})`;
  }

  toString(): string {
    return this._repr_();
  }

  /**
   * A key identifying this place inside its function field, used where
   * SageMath uses ``__hash__``.
   */
  _key(): string {
    return `${this._isInfiniteOrder() ? 0 : 1}|${this._repr_()}`;
  }

  private _isInfiniteOrder(): boolean {
    return this._prime.ring() instanceof FunctionFieldOrderInfinite;
  }

  /**
   * Compare this place with ``other``: returns -1, 0 or 1.
   *
   * Places at infinity are ordered first, exactly as in SageMath.
   *
   * @see Reference: sage/rings/function_field/place.py:146 (_richcmp_)
   */
  cmp(other: FunctionFieldPlace<C>): number {
    const s = this._isInfiniteOrder() ? 0 : 1;
    const o = other._isInfiniteOrder() ? 0 : 1;
    if (s !== o) {
      return s < o ? -1 : 1;
    }
    return this._prime.cmp(other._prime);
  }

  eq(other: FunctionFieldPlace<C>): boolean {
    return this.cmp(other) === 0;
  }

  /**
   * Return the function field to which the place belongs.
   *
   * @see Reference: sage/rings/function_field/place.py:273 (function_field)
   */
  function_field(): FunctionField<C> {
    return this._parent._field;
  }

  /**
   * Return the prime ideal associated with the place.
   *
   * @see Reference: sage/rings/function_field/place.py:287 (prime_ideal)
   */
  prime_ideal(): FunctionFieldIdeal<C> {
    return this._prime;
  }

  /**
   * Return the prime divisor corresponding to the place.
   *
   * @see Reference: sage/rings/function_field/place.py:302 (divisor)
   */
  divisor(multiplicity: bigint | number = 1): FunctionFieldDivisor<C> {
    return prime_divisor(this.function_field(), this, BigInt(multiplicity));
  }

  /**
   * Return the divisor ``m * self``.
   *
   * @see Reference: sage/rings/function_field/place.py:170 (_acted_upon_)
   */
  scalar_mul(m: bigint | number): FunctionFieldDivisor<C> {
    return this.divisor().scalar_mul(BigInt(m));
  }

  /**
   * Return the negative of the prime divisor of this place.
   *
   * @see Reference: sage/rings/function_field/place.py:192 (_neg_)
   */
  neg(): FunctionFieldDivisor<C> {
    return this.divisor(-1);
  }

  /**
   * Return the divisor that is the sum of the place and ``other``.
   *
   * @see Reference: sage/rings/function_field/place.py:209 (_add_)
   */
  add(other: FunctionFieldPlace<C> | FunctionFieldDivisor<C>): FunctionFieldDivisor<C> {
    return this.divisor().add(other);
  }

  /**
   * Return the divisor that is this place minus ``other``.
   *
   * @see Reference: sage/rings/function_field/place.py:227 (_sub_)
   */
  sub(other: FunctionFieldPlace<C> | FunctionFieldDivisor<C>): FunctionFieldDivisor<C> {
    return this.divisor().sub(other);
  }

  /** Return the degree of the place. */
  abstract degree(): bigint;

  /** Return ``true`` if the place is at infinity. */
  abstract is_infinite_place(): boolean;

  /** Return a local uniformizer of the place. */
  abstract local_uniformizer(): FunctionFieldElement<C>;

  /**
   * Return the residue field of the place together with the maps from and to
   * it.
   */
  abstract _residue_field(
    name?: string
  ): [ConstantField<C>, (e: C) => FunctionFieldElement<C>, (f: FunctionFieldElement<C>) => C];

  /** Return the residue field of the place. */
  abstract residue_field(
    name?: string
  ): [ConstantField<C>, (e: C) => FunctionFieldElement<C>, (f: FunctionFieldElement<C>) => C];

  /** Return the valuation ring at the place. */
  abstract valuation_ring(): FunctionFieldValuationRing<C>;
}

/**
 * Sets of places of function fields.
 *
 * @see Reference: sage/rings/function_field/place.py:321 (PlaceSet)
 */
export class PlaceSet<C extends ConstantFieldElement> {
  readonly _field: FunctionField<C>;

  constructor(field: FunctionField<C>) {
    this._field = field;
  }

  /**
   * @see Reference: sage/rings/function_field/place.py:354 (_repr_)
   */
  _repr_(): string {
    return `Set of places of ${this._field}`;
  }

  toString(): string {
    return this._repr_();
  }

  /**
   * Create a place from ``x`` if ``x`` is a prime ideal.
   *
   * @see Reference: sage/rings/function_field/place.py:367 (_element_constructor_)
   */
  __call__(x: unknown): FunctionFieldPlace<C> {
    if (x instanceof FunctionFieldIdeal && x.is_prime()) {
      return this._field._place_class(this, x as FunctionFieldIdeal<C>);
    }
    throw new ValueError('not a prime ideal');
  }

  /**
   * Return the function field to which this place set belongs.
   *
   * @see Reference: sage/rings/function_field/place.py:412 (function_field)
   */
  function_field(): FunctionField<C> {
    return this._field;
  }
}
