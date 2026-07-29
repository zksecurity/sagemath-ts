/**
 * @module sage/rings/function_field/place_rational
 * @description Places of function fields: rational
 *
 * Port of: sage/rings/function_field/place_rational.py
 */

import type { ConstantField, ConstantFieldElement } from './constant_field.js';
import { divide_constants } from './constant_field.js';
import type { FunctionFieldElement } from './element.js';
import type { FunctionFieldElement_rational } from './element_rational.js';
import type { RationalFunctionField } from './function_field_rational.js';
import type { FunctionFieldMaximalOrder_rational } from './order_rational.js';
import { FunctionFieldPlace } from './place.js';
import { FunctionFieldValuationRing } from './valuation_ring.js';

/**
 * Places of rational function fields.
 *
 * @see Reference: sage/rings/function_field/place_rational.py:20 (FunctionFieldPlace_rational)
 */
export class FunctionFieldPlace_rational<
  C extends ConstantFieldElement,
> extends FunctionFieldPlace<C> {
  /**
   * Return the degree of the place.
   *
   * @see Reference: sage/rings/function_field/place_rational.py:24 (degree)
   */
  override degree(): bigint {
    if (this.is_infinite_place()) {
      return 1n;
    }
    const gen = this._prime.gens()[0] as FunctionFieldElement_rational<C>;
    return BigInt(gen.numerator().degree());
  }

  /**
   * Return ``true`` if the place is at infinity.
   *
   * @see Reference: sage/rings/function_field/place_rational.py:42 (is_infinite_place)
   */
  override is_infinite_place(): boolean {
    const F = this.function_field();
    return this.prime_ideal().ring() === F.maximal_order_infinite();
  }

  /**
   * Return a local uniformizer of the place.
   *
   * @see Reference: sage/rings/function_field/place_rational.py:57 (local_uniformizer)
   */
  override local_uniformizer(): FunctionFieldElement_rational<C> {
    return this.prime_ideal().gens()[0] as FunctionFieldElement_rational<C>;
  }

  /**
   * Return the residue field of the place.
   *
   * @see Reference: sage/rings/function_field/place_rational.py:71 (residue_field)
   */
  override residue_field(
    name?: string
  ): [ConstantField<C>, (e: C) => FunctionFieldElement<C>, (f: FunctionFieldElement<C>) => C] {
    return this.valuation_ring().residue_field(name);
  }

  /**
   * Return the residue field of the place along with the maps from and to it.
   *
   * @see Reference: sage/rings/function_field/place_rational.py:94 (_residue_field)
   */
  override _residue_field(
    name?: string
  ): [ConstantField<C>, (e: C) => FunctionFieldElement<C>, (f: FunctionFieldElement<C>) => C] {
    const F = this.function_field() as RationalFunctionField<C>;
    const prime = this.prime_ideal();

    if (this.is_infinite_place()) {
      const K = F.constant_base_field();
      const from_K = (e: C): FunctionFieldElement_rational<C> => F.__call__(e);
      const to_K = (fe: FunctionFieldElement<C>): C => {
        const f = fe as FunctionFieldElement_rational<C>;
        const n = f.numerator();
        const d = f.denominator();
        const n_deg = n.degree();
        const d_deg = d.degree();
        if (n_deg < d_deg) {
          return K.zero();
        }
        if (n_deg === d_deg) {
          return divide_constants(n.leading_coefficient(), d.leading_coefficient());
        }
        throw new TypeError('not in the valuation ring');
      };
      return [K, from_K, to_K];
    }

    const O = F.maximal_order() as FunctionFieldMaximalOrder_rational<C>;
    const [K, from_K, _to_K] = O._residue_field(prime, name);

    const to_K = (fe: FunctionFieldElement<C>): C => {
      const f = fe as FunctionFieldElement_rational<C>;
      if (O.contains(f)) {
        // f.denominator() is 1
        return _to_K(f.numerator());
      }
      const d = F.__call__(f.denominator());
      const n = d.mul(f);

      const nv = (
        prime as unknown as {
          valuation(i: unknown): bigint | number;
        }
      ).valuation(O.ideal(n));
      const dv = (
        prime as unknown as {
          valuation(i: unknown): bigint | number;
        }
      ).valuation(O.ideal(d));

      if (nv > dv) {
        return K.zero();
      }
      if (dv > nv) {
        throw new TypeError('not in the valuation ring');
      }

      const s = (prime.gens()[0] as FunctionFieldElement_rational<C>).inv();
      const rd = d.mul(s.pow(dv as bigint));
      const rn = n.mul(s.pow(nv as bigint));
      return divide_constants(to_K(rn), to_K(rd));
    };

    return [K, from_K, to_K];
  }

  /**
   * Return the valuation ring at the place.
   *
   * @see Reference: sage/rings/function_field/place_rational.py:166 (valuation_ring)
   */
  override valuation_ring(): FunctionFieldValuationRing<C> {
    return new FunctionFieldValuationRing(this.function_field(), this);
  }
}
