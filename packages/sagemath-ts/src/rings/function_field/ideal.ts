/**
 * @module sage/rings/function_field/ideal
 * @description Ideals of function fields
 *
 * Port of: sage/rings/function_field/ideal.py
 */

import { NotImplementedError, ValueError } from '../../errors.js';
import type { ConstantFieldElement } from './constant_field.js';
import { divisor } from './divisor.js';
import type { FunctionFieldDivisor } from './divisor.js';
import type { FunctionFieldElement } from './element.js';
import type { FunctionFieldOrder_base } from './order.js';
import type { FunctionFieldPlace } from './place.js';

/**
 * Base class of fractional ideals of function fields.
 *
 * @see Reference: sage/rings/function_field/ideal.py:107 (FunctionFieldIdeal)
 */
export abstract class FunctionFieldIdeal<C extends ConstantFieldElement> {
  readonly _ring: FunctionFieldOrder_base<C>;

  constructor(ring: FunctionFieldOrder_base<C>) {
    this._ring = ring;
  }

  /** Return the generators of this ideal. */
  abstract gens(): Array<FunctionFieldElement<C>>;

  /**
   * Return the generators of this ideal as a module over the maximal order of
   * the base rational function field.
   */
  abstract gens_over_base(): Array<FunctionFieldElement<C>>;

  /**
   * Return ``true`` if this is the zero ideal.
   *
   * SageMath inherits `is_zero()` from `Element`, whose `__bool__` falls back
   * to ``True`` when the parent has no meaningful zero — and the parent here is
   * the *multiplicative* `IdealMonoid`.  So upstream `O.ideal(0).is_zero()` is
   * `False`, and the zero ideal prints as ``Ideal (0) of ...`` rather than
   * ``Zero ideal of ...``.  Verified by executing SageMath 10.3.
   *
   * @see Reference: sage/structure/element.pyx:998 (Element.__bool__)
   */
  is_zero(): boolean {
    return false;
  }

  /** Return ``true`` if this is a prime ideal. */
  abstract is_prime(): boolean;

  /**
   * Return the list of prime ideal / multiplicity pairs of the factorization
   * of this ideal.
   */
  abstract _factor(): Array<[FunctionFieldIdeal<C>, bigint]>;

  /** Multiply this ideal with ``other``. */
  abstract mul(other: FunctionFieldIdeal<C>): FunctionFieldIdeal<C>;

  /** Add this ideal to ``other``. */
  abstract add(other: FunctionFieldIdeal<C>): FunctionFieldIdeal<C>;

  /** Return the inverse of this fractional ideal. */
  abstract inv(): FunctionFieldIdeal<C>;

  /** Compare this ideal with ``other``; returns -1, 0 or 1. */
  abstract cmp(other: FunctionFieldIdeal<C>): number;

  /**
   * Return the valuation of ``ideal`` at this prime ideal.
   *
   * The value is `Number.POSITIVE_INFINITY` for the zero ideal, mirroring
   * SageMath's `+Infinity`.
   */
  abstract valuation(ideal: FunctionFieldIdeal<C>): bigint | number;

  /**
   * Return this ideal divided by the ``other`` ideal.
   *
   * @see Reference: sage/rings/function_field/ideal.py:233 (_div_)
   */
  div(other: FunctionFieldIdeal<C>): FunctionFieldIdeal<C> {
    return this.mul(other.inv());
  }

  /**
   * Return this ideal raised to the ``n``-th power (``n`` may be negative).
   *
   * SageMath gets this from the generic monoid element power; a negative
   * exponent goes through ``__invert__`` first.
   */
  pow(n: bigint | number): FunctionFieldIdeal<C> {
    let e = BigInt(n);
    // eslint-disable-next-line @typescript-eslint/no-this-alias
    let base: FunctionFieldIdeal<C> = this;
    if (e < 0n) {
      base = base.inv();
      e = -e;
    }
    let result = this._ring.ideal(this._ring.function_field().one());
    while (e > 0n) {
      if (e & 1n) {
        result = result.mul(base);
      }
      base = base.mul(base);
      e >>= 1n;
    }
    return result;
  }

  /** Return ``true`` if the two ideals are equal. */
  eq(other: FunctionFieldIdeal<C>): boolean {
    return this.cmp(other) === 0;
  }

  /**
   * Return reduced generators.
   *
   * @see Reference: sage/rings/function_field/ideal.py:252 (gens_reduced)
   */
  gens_reduced(): Array<FunctionFieldElement<C>> {
    const gens = this.gens();
    if (gens.length === 1) {
      return gens;
    }
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: FunctionFieldIdeal.gens_reduced for ideals with more than one generator'
    );
  }

  /**
   * Return the ring to which this ideal belongs.
   *
   * @see Reference: sage/rings/function_field/ideal.py:283 (ring)
   */
  ring(): FunctionFieldOrder_base<C> {
    return this._ring;
  }

  /**
   * Return the base ring of this ideal.
   *
   * @see Reference: sage/rings/function_field/ideal.py:297 (base_ring)
   */
  base_ring(): FunctionFieldOrder_base<C> {
    return this.ring();
  }

  /**
   * Return a string representation of this ideal that does not include the
   * name of the ambient ring.
   *
   * @see Reference: sage/rings/function_field/ideal.py:137 (_repr_short)
   */
  _repr_short(): string {
    if (this.is_zero()) {
      return '(0)';
    }
    return `(${this.gens_reduced()
      .map((g) => g.toString())
      .join(', ')})`;
  }

  /**
   * @see Reference: sage/rings/function_field/ideal.py:157 (_repr_)
   */
  _repr_(): string {
    if (this.is_zero()) {
      return `Zero ideal of ${this._ring}`;
    }
    return `Ideal ${this._repr_short()} of ${this.ring()}`;
  }

  toString(): string {
    return this._repr_();
  }

  /**
   * Return the place associated with this prime ideal.
   *
   * @see Reference: sage/rings/function_field/ideal.py:313 (place)
   */
  place(): FunctionFieldPlace<C> {
    if (!this.is_prime()) {
      throw new TypeError('not a prime ideal');
    }
    const place_set = this.ring().fraction_field().place_set();
    return place_set.__call__(this);
  }

  /**
   * Return the factorization of this ideal.
   *
   * @see Reference: sage/rings/function_field/ideal.py:389 (factor)
   */
  factor(): Array<[FunctionFieldIdeal<C>, bigint]> {
    // SageMath wraps `_factor()` in a `Factorization`, whose constructor sorts
    // the factors.  Ideals have neither `dimension` nor `degree`, so the sort
    // key falls through to the prime itself.
    // @see Reference: sage/structure/factorization.py:671 (sort)
    return [...this._factor()].sort((a, b) => a[0].cmp(b[0]));
  }

  /**
   * Return the divisor corresponding to the ideal.
   *
   * @see Reference: sage/rings/function_field/ideal.py:451 (divisor)
   */
  divisor(): FunctionFieldDivisor<C> {
    if (this.is_zero()) {
      throw new ValueError('not defined for zero ideal');
    }
    const F = this.ring().fraction_field();
    const data: Array<[FunctionFieldPlace<C>, bigint]> = this._factor().map(
      ([prime, multiplicity]) => [prime.place(), multiplicity]
    );
    return divisor(F, data);
  }

  /**
   * Return the divisor of zeros corresponding to the ideal.
   *
   * @see Reference: sage/rings/function_field/ideal.py:510 (divisor_of_zeros)
   */
  divisor_of_zeros(): FunctionFieldDivisor<C> {
    if (this.is_zero()) {
      throw new ValueError('not defined for zero ideal');
    }
    const F = this.ring().fraction_field();
    const data: Array<[FunctionFieldPlace<C>, bigint]> = this._factor()
      .filter(([, m]) => m > 0n)
      .map(([prime, m]) => [prime.place(), m]);
    return divisor(F, data);
  }

  /**
   * Return the divisor of poles corresponding to the ideal.
   *
   * @see Reference: sage/rings/function_field/ideal.py:547 (divisor_of_poles)
   */
  divisor_of_poles(): FunctionFieldDivisor<C> {
    if (this.is_zero()) {
      throw new ValueError('not defined for zero ideal');
    }
    const F = this.ring().fraction_field();
    const data: Array<[FunctionFieldPlace<C>, bigint]> = this._factor()
      .filter(([, m]) => m < 0n)
      .map(([prime, m]) => [prime.place(), -m]);
    return divisor(F, data);
  }
}

/**
 * Base class of ideals of maximal infinite orders.
 *
 * @see Reference: sage/rings/function_field/ideal.py:884 (FunctionFieldIdealInfinite)
 */
export abstract class FunctionFieldIdealInfinite<
  C extends ConstantFieldElement,
> extends FunctionFieldIdeal<C> {}

/**
 * The monoid of ideals in orders of function fields.
 *
 * @see Reference: sage/rings/function_field/ideal.py:1017 (IdealMonoid)
 */
export class IdealMonoid<C extends ConstantFieldElement> {
  private readonly __R: FunctionFieldOrder_base<C>;

  constructor(R: FunctionFieldOrder_base<C>) {
    this.__R = R;
  }

  /**
   * @see Reference: sage/rings/function_field/ideal.py:1035 (_repr_)
   */
  _repr_(): string {
    return `Monoid of ideals of ${this.__R}`;
  }

  toString(): string {
    return this._repr_();
  }

  /**
   * Return the ring of which this is the ideal monoid.
   *
   * @see Reference: sage/rings/function_field/ideal.py:1048 (ring)
   */
  ring(): FunctionFieldOrder_base<C> {
    return this.__R;
  }

  /**
   * Create an ideal in the monoid from ``x``.
   *
   * @see Reference: sage/rings/function_field/ideal.py:1062 (_element_constructor_)
   */
  __call__(x: unknown): FunctionFieldIdeal<C> {
    if (x instanceof FunctionFieldIdeal) {
      return this.__R.ideal(x.gens());
    }
    return this.__R.ideal(x);
  }
}
