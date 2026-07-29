/**
 * @module sage/rings/function_field/element
 * @description Elements of function fields
 *
 * Port of: sage/rings/function_field/element.pyx
 */

import { NotImplementedError, ValueError } from '../../errors.js';
import type { ConstantFieldElement } from './constant_field.js';
import type { FunctionFieldDivisor } from './divisor.js';
import type { FunctionField } from './function_field.js';
import type { FunctionFieldPlace } from './place.js';

/**
 * Abstract base class of function field elements.
 *
 * @see Reference: sage/rings/function_field/element.pyx:117 (FunctionFieldElement)
 */
export abstract class FunctionFieldElement<C extends ConstantFieldElement> {
  abstract get parent(): FunctionField<C>;

  abstract is_zero(): boolean;
  abstract is_one(): boolean;
  abstract toString(): string;

  // Field arithmetic; SageMath inherits these from ``FieldElement``.
  abstract add(other: FunctionFieldElement<C>): FunctionFieldElement<C>;
  abstract sub(other: FunctionFieldElement<C>): FunctionFieldElement<C>;
  abstract mul(other: FunctionFieldElement<C>): FunctionFieldElement<C>;
  abstract div(other: FunctionFieldElement<C>): FunctionFieldElement<C>;
  abstract neg(): FunctionFieldElement<C>;
  abstract inv(): FunctionFieldElement<C>;
  abstract pow(n: bigint | number): FunctionFieldElement<C>;
  abstract eq(other: FunctionFieldElement<C>): boolean;
  /** Multiply by an element of the constant field. */
  abstract scalar_mul(c: C): FunctionFieldElement<C>;

  /**
   * Return the divisor of the element.
   *
   * @see Reference: sage/rings/function_field/element.pyx:756 (divisor)
   */
  divisor(): FunctionFieldDivisor<C> {
    if (this.is_zero()) {
      throw new ValueError('divisor not defined for zero');
    }
    const F = this.parent;
    const I = F.maximal_order().ideal(this);
    const J = F.maximal_order_infinite().ideal(this);
    return I.divisor().add(J.divisor());
  }

  /**
   * Return the divisor of zeros for the element.
   *
   * @see Reference: sage/rings/function_field/element.pyx:786 (divisor_of_zeros)
   */
  divisor_of_zeros(): FunctionFieldDivisor<C> {
    if (this.is_zero()) {
      throw new ValueError('divisor of zeros not defined for zero');
    }
    const F = this.parent;
    const I = F.maximal_order().ideal(this);
    const J = F.maximal_order_infinite().ideal(this);
    return I.divisor_of_zeros().add(J.divisor_of_zeros());
  }

  /**
   * Return the divisor of poles for the element.
   *
   * @see Reference: sage/rings/function_field/element.pyx:812 (divisor_of_poles)
   */
  divisor_of_poles(): FunctionFieldDivisor<C> {
    if (this.is_zero()) {
      throw new ValueError('divisor of poles not defined for zero');
    }
    const F = this.parent;
    const I = F.maximal_order().ideal(this);
    const J = F.maximal_order_infinite().ideal(this);
    return I.divisor_of_poles().add(J.divisor_of_poles());
  }

  /**
   * Return the list of the zeros of the element.
   *
   * @see Reference: sage/rings/function_field/element.pyx:839 (zeros)
   */
  zeros(): Array<FunctionFieldPlace<C>> {
    return this.divisor_of_zeros().support();
  }

  /**
   * Return the list of the poles of the element.
   *
   * @see Reference: sage/rings/function_field/element.pyx:859 (poles)
   */
  poles(): Array<FunctionFieldPlace<C>> {
    return this.divisor_of_poles().support();
  }

  /**
   * Return the value of the element at the place.
   *
   * If the element is in the valuation ring at the place an element of the
   * residue field is returned; otherwise a {@link ValueError} is raised.
   *
   * @see Reference: sage/rings/function_field/element.pyx:909 (evaluate)
   */
  evaluate(place: FunctionFieldPlace<C>): C {
    const [R, , to_R] = place._residue_field();

    const v = this.valuation(place);
    if (v > 0) {
      return R.zero();
    }
    if (v < 0) {
      throw new ValueError('has a pole at the place');
    }
    return to_R(this);
  }

  /**
   * Return the valuation of the element at the place.
   *
   * The value is `Number.POSITIVE_INFINITY` for the zero element, mirroring
   * SageMath's `+Infinity`.
   *
   * @see Reference: sage/rings/function_field/element.pyx:879 (valuation)
   */
  abstract valuation(place: FunctionFieldPlace<C>): bigint | number;

  /**
   * Return whether this element is an ``n``-th power in the function field.
   *
   * @see Reference: sage/rings/function_field/element.pyx:955 (is_nth_power)
   */
  is_nth_power(_n: bigint | number): boolean {
    throw new NotImplementedError('is_nth_power() not implemented for generic elements');
  }

  /**
   * Return an ``n``-th root of this element in the function field.
   *
   * @see Reference: sage/rings/function_field/element.pyx:982 (nth_root)
   */
  nth_root(_n: bigint | number): FunctionFieldElement<C> {
    throw new NotImplementedError('nth_root() not implemented for generic elements');
  }

  /**
   * Return the differential ``d(self)``.
   *
   * @see Reference: sage/rings/function_field/element.pyx:666 (differential)
   */
  differential(): never {
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: FunctionFieldElement.differential (sage/rings/function_field/differential.py)'
    );
  }

  /**
   * Return the higher derivative ``D^{(i)}(self)``.
   *
   * @see Reference: sage/rings/function_field/element.pyx:726 (higher_derivative)
   */
  higher_derivative(_i: number, _separating_element?: unknown): never {
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: FunctionFieldElement.higher_derivative (sage/rings/function_field/derivations_polymod.py)'
    );
  }
}

/**
 * Return ``true`` if ``x`` is a function field element.
 *
 * @see Reference: sage/rings/function_field/element.pyx:75 (is_FunctionFieldElement)
 */
export function is_FunctionFieldElement(x: unknown): boolean {
  return x instanceof FunctionFieldElement;
}
