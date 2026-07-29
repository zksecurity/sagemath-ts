/**
 * @module sage/rings/function_field/valuation_ring
 * @description Valuation rings of function fields
 *
 * Port of: sage/rings/function_field/valuation_ring.py
 */

import type { ConstantField, ConstantFieldElement } from './constant_field.js';
import type { FunctionFieldElement } from './element.js';
import type { FunctionField } from './function_field.js';
import type { FunctionFieldPlace } from './place.js';

/**
 * Valuation ring of a function field at a place.
 *
 * @see Reference: sage/rings/function_field/valuation_ring.py:76 (FunctionFieldValuationRing)
 */
export class FunctionFieldValuationRing<C extends ConstantFieldElement> {
  readonly _field: FunctionField<C>;
  readonly _place: FunctionFieldPlace<C>;

  constructor(field: FunctionField<C>, place: FunctionFieldPlace<C>) {
    this._field = field;
    this._place = place;
  }

  /**
   * Construct an element of the function field belonging to the valuation ring.
   *
   * @see Reference: sage/rings/function_field/valuation_ring.py:111 (_element_constructor_)
   */
  __call__(x: unknown): FunctionFieldElement<C> {
    const e = this._field.__call__(x);
    if (e.valuation(this._place) >= 0) {
      return e;
    }
    throw new TypeError();
  }

  /**
   * @see Reference: sage/rings/function_field/valuation_ring.py:137 (_repr_)
   */
  _repr_(): string {
    return `Valuation ring at ${this._place}`;
  }

  toString(): string {
    return this._repr_();
  }

  /**
   * Return the place associated with the valuation ring.
   *
   * @see Reference: sage/rings/function_field/valuation_ring.py:151 (place)
   */
  place(): FunctionFieldPlace<C> {
    return this._place;
  }

  /**
   * Return the residue field of the valuation ring together with the maps from
   * and to it.
   *
   * SageMath wraps the two maps in `FunctionFieldRingMorphism` objects; this
   * port returns plain functions.
   *
   * @see Reference: sage/rings/function_field/valuation_ring.py:166 (residue_field)
   * @see Deviation: function-field residue field maps returned as plain functions
   */
  residue_field(
    name?: string
  ): [ConstantField<C>, (e: C) => FunctionFieldElement<C>, (f: FunctionFieldElement<C>) => C] {
    return this._place._residue_field(name);
  }
}
