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
import { NotImplementedError } from '../../errors.js';

/**
 * Valuation ring of a function field at a place.
 *
 * @see Reference: sage/rings/function_field/valuation_ring.py:76 (FunctionFieldValuationRing)
 */
export class FunctionFieldValuationRing<C extends ConstantFieldElement> {
  private static readonly _cache = new WeakMap<object, Map<string, unknown>>();
  private readonly _residue_fields = new Map<string | undefined,
    [ConstantField<C>, (e: C) => FunctionFieldElement<C>, (f: FunctionFieldElement<C>) => C]>();
  readonly _field: FunctionField<C>;
  readonly _place: FunctionFieldPlace<C>;

  constructor(field: FunctionField<C>, place: FunctionFieldPlace<C>) {
    this._field = field;
    this._place = place;
    let cache = FunctionFieldValuationRing._cache.get(field);
    if (!cache) {
      cache = new Map();
      FunctionFieldValuationRing._cache.set(field, cache);
    }
    const key = place._key();
    const existing = cache.get(key);
    if (existing) return existing as FunctionFieldValuationRing<C>;
    cache.set(key, this);
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
    const existing = this._residue_fields.get(name);
    if (existing) return existing;
    const [k, from_k, to_k] = this._place._residue_field(name);
    // categories/map.pyx: Map.__call__ converts to the map's domain first.
    const convert = <T>(x: unknown, domain: { __call__(x: unknown): T; toString(): string }): T => {
      try {
        return domain.__call__(x);
      } catch (error) {
        if (!(error instanceof TypeError || error instanceof NotImplementedError)) throw error;
        throw new TypeError(`${x} fails to convert into the map's domain ${domain}, but a \`pushforward\` method is not properly implemented`);
      }
    };
    const result: [ConstantField<C>, (e: C) => FunctionFieldElement<C>, (f: FunctionFieldElement<C>) => C] =
      [k, (e) => from_k(convert(e, k)), (f) => to_k(convert(f, this))];
    this._residue_fields.set(name, result);
    return result;
  }
}
