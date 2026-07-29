/**
 * @module sage/rings/function_field/function_field
 * @description Function fields
 *
 * A function field (of one variable) is a finitely generated field extension of
 * transcendence degree one.  In SageMath a function field can be a rational
 * function field or a finite extension of a function field; this port currently
 * implements the rational case (see `function_field_rational.ts`).
 *
 * Port of: sage/rings/function_field/function_field.py
 */

import { NotImplementedError } from '../../errors.js';
import { constant_field_characteristic, constant_field_is_finite } from './constant_field.js';
import type { ConstantField, ConstantFieldElement } from './constant_field.js';
import { DivisorGroup } from './divisor.js';
import type { FunctionFieldDivisor } from './divisor.js';
import type { FunctionFieldElement } from './element.js';
import type { FunctionFieldIdeal } from './ideal.js';
import type { FunctionFieldOrder_base } from './order.js';
import { PlaceSet } from './place.js';
import type { FunctionFieldPlace } from './place.js';

/**
 * Abstract base class for all function fields.
 *
 * @see Reference: sage/rings/function_field/function_field.py:288 (FunctionField)
 */
export abstract class FunctionField<C extends ConstantFieldElement> {
  protected _divisor_group: DivisorGroup<C> | null = null;
  protected _place_set: PlaceSet<C> | null = null;

  /** Return the constant field of which this field is a transcendental extension. */
  abstract constant_base_field(): ConstantField<C>;

  /** Return the name of the generator. */
  abstract variable_name(): string;

  /** Return the base field. */
  abstract base_field(): FunctionField<C>;

  /** Return the rational function field at the bottom of the tower. */
  abstract rational_function_field(): FunctionField<C>;

  /** Return the `n`-th generator. */
  abstract gen(n?: number): FunctionFieldElement<C>;

  /** Return the number of generators. */
  abstract ngens(): number;

  /** Return the degree over ``base``. */
  abstract degree(base?: FunctionField<C>): bigint;

  /** Return the genus of the function field. */
  abstract genus(): bigint;

  /** Return the maximal (finite) order. */
  abstract maximal_order(): FunctionFieldOrder_base<C>;

  /** Return the maximal infinite order. */
  abstract maximal_order_infinite(): FunctionFieldOrder_base<C>;

  /** Return the different of the function field, as a divisor. */
  abstract different(): FunctionFieldDivisor<C>;

  /** Coerce ``x`` into an element of the function field. */
  abstract __call__(x: unknown): FunctionFieldElement<C>;

  abstract zero(): FunctionFieldElement<C>;
  abstract one(): FunctionFieldElement<C>;

  abstract _repr_(): string;

  /**
   * Build the place class of this function field.
   *
   * SageMath stores the class itself in ``self._place_class``
   * (`function_field_rational.py:162`); TypeScript has no comparable
   * "class attribute plus `element_class`" mechanism, so subclasses implement
   * this factory instead.
   */
  abstract _place_class(parent: PlaceSet<C>, prime: FunctionFieldIdeal<C>): FunctionFieldPlace<C>;

  toString(): string {
    return this._repr_();
  }

  /** Alias of {@link constant_base_field}, as in SageMath. */
  constant_field(): ConstantField<C> {
    return this.constant_base_field();
  }

  /**
   * Return the characteristic of the function field.
   *
   * @see Reference: sage/rings/function_field/function_field.py:389 (characteristic)
   */
  characteristic(): bigint {
    return constant_field_characteristic(this.constant_base_field());
  }

  /**
   * Return whether the function field is finite, which is false.
   *
   * @see Reference: sage/rings/function_field/function_field.py:411 (is_finite)
   */
  is_finite(): false {
    return false;
  }

  /**
   * Return whether the function field is global, that is, whether the constant
   * field is finite.
   *
   * @see Reference: sage/rings/function_field/function_field.py:426 (is_global)
   */
  is_global(): boolean {
    return constant_field_is_finite(this.constant_base_field());
  }

  /**
   * Return whether the field is perfect, i.e. its characteristic is zero.
   *
   * @see Reference: sage/rings/function_field/function_field.py:326 (is_perfect)
   */
  is_perfect(): boolean {
    return this.characteristic() === 0n;
  }

  /** Return ``true``; function fields are fields. */
  is_field(_proof: boolean = true): boolean {
    return true;
  }

  /**
   * Return the group of divisors attached to the function field.
   *
   * @see Reference: sage/rings/function_field/function_field.py:1102 (divisor_group)
   */
  divisor_group(): DivisorGroup<C> {
    if (this._divisor_group === null) {
      this._divisor_group = new DivisorGroup(this);
    }
    return this._divisor_group;
  }

  /**
   * Return the set of all places of the function field.
   *
   * @see Reference: sage/rings/function_field/function_field.py:1125 (place_set)
   */
  place_set(): PlaceSet<C> {
    if (this._place_set === null) {
      this._place_set = new PlaceSet(this);
    }
    return this._place_set;
  }

  /**
   * Create an extension field by adjoining a root of ``f``.
   *
   * @see Reference: sage/rings/function_field/function_field.py:445 (extension)
   */
  extension(_f: unknown, _names?: string): never {
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: FunctionField.extension ' +
        '(sage/rings/function_field/function_field_polymod.py)'
    );
  }

  /**
   * Return the space of differentials of the function field.
   *
   * @see Reference: sage/rings/function_field/function_field.py:1029 (space_of_differentials)
   */
  space_of_differentials(): never {
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: FunctionField.space_of_differentials ' +
        '(sage/rings/function_field/differential.py)'
    );
  }

  /**
   * Return the space of holomorphic differentials of the function field.
   *
   * @see Reference: sage/rings/function_field/function_field.py:1047 (space_of_holomorphic_differentials)
   */
  space_of_holomorphic_differentials(): never {
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: FunctionField.space_of_holomorphic_differentials ' +
        '(sage/rings/function_field/differential.py)'
    );
  }

  /**
   * Return the completion of the function field at a place.
   *
   * @see Reference: sage/rings/function_field/function_field.py:1148 (completion)
   */
  completion(_place: unknown, _name?: string, _prec?: number, _gen_name?: string): never {
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: FunctionField.completion (sage/rings/function_field/maps.py)'
    );
  }

  /**
   * Return the valuation of the function field at ``prime``.
   *
   * @see Reference: sage/rings/function_field/function_field.py:907 (valuation)
   */
  valuation(_prime: unknown): never {
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: FunctionField.valuation (sage/rings/function_field/valuation.py)'
    );
  }

  /**
   * Return the Jacobian of the function field.
   *
   * @see Reference: sage/rings/function_field/function_field.py:1391 (jacobian)
   */
  jacobian(_model?: string, _base_div?: unknown): never {
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: FunctionField.jacobian ' +
        '(sage/rings/function_field/jacobian_hess.py, jacobian_khuri_makdisi.py)'
    );
  }

  /**
   * Return the Hilbert symbol `(a, b)_P`.
   *
   * @see Reference: sage/rings/function_field/function_field.py:1251 (hilbert_symbol)
   */
  hilbert_symbol(_a: unknown, _b: unknown, _P: unknown): never {
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: FunctionField.hilbert_symbol (needs completions)'
    );
  }

  /**
   * Return the higher derivation (Hasse-Schmidt derivation) of the field.
   *
   * @see Reference: sage/rings/function_field/function_field_rational.py:828 (higher_derivation)
   */
  higher_derivation(): never {
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: FunctionField.higher_derivation ' +
        '(sage/rings/function_field/derivations_polymod.py)'
    );
  }
}

/**
 * Return ``true`` if ``x`` is a function field.
 *
 * @see Reference: sage/rings/function_field/function_field.py:264 (is_FunctionField)
 */
export function is_FunctionField(x: unknown): boolean {
  return x instanceof FunctionField;
}
