/**
 * @module sage/rings/function_field/ideal_rational
 * @description Ideals of function fields: rational
 *
 * Port of: sage/rings/function_field/ideal_rational.py
 */

import type { ConstantFieldElement } from './constant_field.js';
import type { FunctionFieldElement_rational } from './element_rational.js';
import { FunctionFieldIdeal, FunctionFieldIdealInfinite } from './ideal.js';
import type { FunctionFieldOrder_base } from './order.js';

/**
 * Fractional ideals of the maximal order of a rational function field.
 *
 * @see Reference: sage/rings/function_field/ideal_rational.py:26 (FunctionFieldIdeal_rational)
 */
export class FunctionFieldIdeal_rational<
  C extends ConstantFieldElement,
> extends FunctionFieldIdeal<C> {
  readonly _gen: FunctionFieldElement_rational<C>;

  constructor(ring: FunctionFieldOrder_base<C>, gen: FunctionFieldElement_rational<C>) {
    super(ring);
    this._gen = gen;
  }

  /**
   * Test if ``element`` is in this ideal.
   *
   * @see Reference: sage/rings/function_field/ideal_rational.py:70 (__contains__)
   */
  contains(element: FunctionFieldElement_rational<C>): boolean {
    if (this._gen.is_zero()) {
      return element.is_zero();
    }
    return this._ring.contains(element.div(this._gen));
  }

  /**
   * @see Reference: sage/rings/function_field/ideal_rational.py:88 (_richcmp_)
   */
  override cmp(other: FunctionFieldIdeal<C>): number {
    return this._gen.cmp((other as FunctionFieldIdeal_rational<C>)._gen);
  }

  /**
   * @see Reference: sage/rings/function_field/ideal_rational.py:114 (_add_)
   */
  override add(other: FunctionFieldIdeal<C>): FunctionFieldIdeal<C> {
    return this._ring.ideal([this._gen, (other as FunctionFieldIdeal_rational<C>)._gen]);
  }

  /**
   * @see Reference: sage/rings/function_field/ideal_rational.py:133 (_mul_)
   */
  override mul(other: FunctionFieldIdeal<C>): FunctionFieldIdeal<C> {
    return this._ring.ideal([this._gen.mul((other as FunctionFieldIdeal_rational<C>)._gen)]);
  }

  /**
   * Multiply ``other`` (a function field element) with this ideal.
   *
   * @see Reference: sage/rings/function_field/ideal_rational.py:152 (_acted_upon_)
   */
  acted_upon(other: FunctionFieldElement_rational<C>): FunctionFieldIdeal<C> {
    return this._ring.ideal([other.mul(this._gen)]);
  }

  /**
   * @see Reference: sage/rings/function_field/ideal_rational.py:172 (__invert__)
   */
  override inv(): FunctionFieldIdeal<C> {
    return this._ring.ideal([this._gen.inv()]);
  }

  /**
   * Return the denominator of this fractional ideal.
   *
   * @see Reference: sage/rings/function_field/ideal_rational.py:187 (denominator)
   */
  denominator(): ReturnType<FunctionFieldElement_rational<C>['denominator']> {
    return this._gen.denominator();
  }

  /**
   * @see Reference: sage/rings/function_field/ideal_rational.py:201 (is_prime)
   */
  override is_prime(): boolean {
    return this._gen.denominator().degree() === 0 && this._gen.numerator().is_irreducible();
  }

  /**
   * @see Reference: sage/rings/function_field/ideal_rational.py:240 (gen)
   */
  gen(): FunctionFieldElement_rational<C> {
    return this._gen;
  }

  /**
   * @see Reference: sage/rings/function_field/ideal_rational.py:255 (gens)
   */
  override gens(): Array<FunctionFieldElement_rational<C>> {
    return [this._gen];
  }

  /**
   * @see Reference: sage/rings/function_field/ideal_rational.py:270 (gens_over_base)
   */
  override gens_over_base(): Array<FunctionFieldElement_rational<C>> {
    return [this._gen];
  }

  /**
   * Return the valuation of ``ideal`` at this prime ideal.
   *
   * @see Reference: sage/rings/function_field/ideal_rational.py:286 (valuation)
   */
  valuation(ideal: FunctionFieldIdeal<C>): bigint | number {
    if (!this.is_prime()) {
      throw new TypeError('not a prime ideal');
    }
    const O = this.ring();
    const K = O.function_field();
    const d = (ideal as FunctionFieldIdeal_rational<C>).denominator();
    const dElt = K.__call__(d) as FunctionFieldElement_rational<C>;
    const scaled = (ideal as FunctionFieldIdeal_rational<C>).acted_upon(dElt);
    const a = this._valuation(scaled);
    const b = this._valuation(O.ideal(dElt));
    if (a === Number.POSITIVE_INFINITY) {
      return Number.POSITIVE_INFINITY;
    }
    return (a as bigint) - (b as bigint);
  }

  /**
   * Return the valuation of the integral ideal at this prime ideal.
   *
   * @see Reference: sage/rings/function_field/ideal_rational.py:309 (_valuation)
   */
  _valuation(ideal: FunctionFieldIdeal<C>): bigint | number {
    return (ideal as FunctionFieldIdeal_rational<C>)._gen.valuation(this._gen);
  }

  /**
   * @see Reference: sage/rings/function_field/ideal_rational.py:333 (_factor)
   */
  override _factor(): Array<[FunctionFieldIdeal<C>, bigint]> {
    return this._gen
      .factor()
      .factors.map(([f, m]) => [this.ring().ideal(f), m] as [FunctionFieldIdeal<C>, bigint]);
  }
}

/**
 * Fractional ideals of the maximal infinite order of a rational function field.
 *
 * Note that the infinite maximal order is a principal ideal domain.
 *
 * @see Reference: sage/rings/function_field/ideal_rational.py:353 (FunctionFieldIdealInfinite_rational)
 */
export class FunctionFieldIdealInfinite_rational<
  C extends ConstantFieldElement,
> extends FunctionFieldIdealInfinite<C> {
  readonly _gen: FunctionFieldElement_rational<C>;

  constructor(ring: FunctionFieldOrder_base<C>, gen: FunctionFieldElement_rational<C>) {
    super(ring);
    this._gen = gen;
  }

  /**
   * @see Reference: sage/rings/function_field/ideal_rational.py:402 (__contains__)
   */
  contains(element: FunctionFieldElement_rational<C>): boolean {
    if (this._gen.is_zero()) {
      return element.is_zero();
    }
    return this._ring.contains(element.div(this._gen));
  }

  /**
   * @see Reference: sage/rings/function_field/ideal_rational.py:426 (_richcmp_)
   */
  override cmp(other: FunctionFieldIdeal<C>): number {
    return this._gen.cmp((other as FunctionFieldIdealInfinite_rational<C>)._gen);
  }

  /**
   * @see Reference: sage/rings/function_field/ideal_rational.py:442 (_add_)
   */
  override add(other: FunctionFieldIdeal<C>): FunctionFieldIdeal<C> {
    return this._ring.ideal([this._gen, (other as FunctionFieldIdealInfinite_rational<C>)._gen]);
  }

  /**
   * @see Reference: sage/rings/function_field/ideal_rational.py:463 (_mul_)
   */
  override mul(other: FunctionFieldIdeal<C>): FunctionFieldIdeal<C> {
    return this._ring.ideal([
      this._gen.mul((other as FunctionFieldIdealInfinite_rational<C>)._gen),
    ]);
  }

  /**
   * @see Reference: sage/rings/function_field/ideal_rational.py:484 (_acted_upon_)
   */
  acted_upon(other: FunctionFieldElement_rational<C>): FunctionFieldIdeal<C> {
    return this._ring.ideal([other.mul(this._gen)]);
  }

  /**
   * @see Reference: sage/rings/function_field/ideal_rational.py:504 (__invert__)
   */
  override inv(): FunctionFieldIdeal<C> {
    return this._ring.ideal([this._gen.inv()]);
  }

  /**
   * @see Reference: sage/rings/function_field/ideal_rational.py:520 (is_prime)
   */
  override is_prime(): boolean {
    const x = this._ring.fraction_field().gen() as FunctionFieldElement_rational<C>;
    return this._gen.eq(x.inv());
  }

  /**
   * @see Reference: sage/rings/function_field/ideal_rational.py:536 (gen)
   */
  gen(): FunctionFieldElement_rational<C> {
    return this._gen;
  }

  /**
   * @see Reference: sage/rings/function_field/ideal_rational.py:551 (gens)
   */
  override gens(): Array<FunctionFieldElement_rational<C>> {
    return [this._gen];
  }

  /**
   * @see Reference: sage/rings/function_field/ideal_rational.py:566 (gens_over_base)
   */
  override gens_over_base(): Array<FunctionFieldElement_rational<C>> {
    return [this._gen];
  }

  /**
   * Return the valuation of ``ideal`` at this prime ideal.
   *
   * @see Reference: sage/rings/function_field/ideal_rational.py:582 (valuation)
   */
  valuation(ideal: FunctionFieldIdeal<C>): bigint | number {
    if (!this.is_prime()) {
      throw new TypeError('not a prime ideal');
    }
    const f = (ideal as FunctionFieldIdealInfinite_rational<C>).gen();
    if (f.is_zero()) {
      return Number.POSITIVE_INFINITY;
    }
    return BigInt(f.denominator().degree() - f.numerator().degree());
  }

  /**
   * @see Reference: sage/rings/function_field/ideal_rational.py:609 (_factor)
   */
  override _factor(): Array<[FunctionFieldIdeal<C>, bigint]> {
    const g = (this.ring().fraction_field().gen() as FunctionFieldElement_rational<C>).inv();
    const m = BigInt(this._gen.denominator().degree() - this._gen.numerator().degree());
    if (m === 0n) {
      return [];
    }
    return [[this.ring().ideal(g), m]];
  }
}
