/**
 * @module sage/rings/number_field/number_field_element
 * @description Elements of number fields
 *
 * Port of: sage/rings/number_field/number_field_element.pyx
 * Reference: reference/sage/src/sage/rings/number_field/number_field_element.pyx
 *
 * This module re-exports the NumberFieldElement from number_field.ts
 * and adds additional element-specific functionality.
 */

// Re-export the main NumberFieldElement from number_field.ts
export { NumberFieldElement, RationalPolynomial } from './number_field.js';

import { NotImplementedError } from '../../errors.js';
import { Rational } from '../rational.js';
import {
  NumberField,
  NumberFieldElement,
  type QuadraticField,
  RationalPolynomial,
} from './number_field.js';

/**
 * An element of a quadratic number field.
 *
 * For Q(sqrt(d)), elements have the form a + b*sqrt(d) where a, b are rational.
 *
 * @see Reference: sage/rings/number_field/number_field_element_quadratic.pyx
 */
export class NumberFieldElement_quadratic extends NumberFieldElement {
  /**
   * Return the "real" part (coefficient of 1).
   */
  real(): Rational {
    return this.__getitem__(0);
  }

  /**
   * Return the "imaginary" part (coefficient of sqrt(d)).
   */
  imag(): Rational {
    return this.__getitem__(1);
  }

  /**
   * Return the complex conjugate (for imaginary quadratic fields).
   * For a + b*sqrt(d), returns a - b*sqrt(d).
   */
  conjugate(): NumberFieldElement_quadratic {
    const coeffs = this.list();
    if (coeffs.length >= 2) {
      coeffs[1] = coeffs[1]!.neg();
    }
    return new NumberFieldElement_quadratic(this.parent() as QuadraticField, coeffs);
  }
}

// Additional standalone functions for number field elements

/**
 * Check if an element is a root of unity.
 *
 * An element is a root of unity if it has finite multiplicative order.
 */
export function is_root_of_unity(element: NumberFieldElement): boolean {
  // An element is a root of unity iff it's a unit with all conjugates on the unit circle
  // For simplicity, we check if the absolute norm is 1 and it's integral
  if (!element.is_integral()) return false;

  const norm = element.norm();
  if (!norm.eq(Rational.one()) && !norm.eq(new Rational(-1n))) {
    return false;
  }

  // More sophisticated check would verify all conjugates have absolute value 1
  // For now, just check if some power equals 1
  const maxOrder = 10000n; // Reasonable bound for roots of unity in number fields
  let current = element;

  for (let i = 1n; i <= maxOrder; i++) {
    if (current.is_one()) return true;
    current = current.mul(element);
  }

  return false;
}

/**
 * Compute the multiplicative order of a root of unity.
 *
 * Returns null if the element is not a root of unity.
 */
export function multiplicative_order(element: NumberFieldElement): bigint | null {
  if (!element.is_integral()) return null;

  const norm = element.norm();
  if (!norm.eq(Rational.one()) && !norm.eq(new Rational(-1n))) {
    return null;
  }

  const maxOrder = 10000n;
  let current = element;

  for (let i = 1n; i <= maxOrder; i++) {
    if (current.is_one()) return i;
    current = current.mul(element);
  }

  return null;
}

/**
 * Check if an element is totally positive (all real embeddings are positive).
 */
export function is_totally_positive(element: NumberFieldElement): boolean {
  // Full implementation requires computing real embeddings
  // For quadratic fields, we can check directly
  const K = element.parent();

  if (K.degree() === 2) {
    const [r1] = K.signature();
    if (r1 === 0) {
      // Imaginary quadratic field - no real embeddings
      return true;
    }

    // Real quadratic field - need both embeddings positive
    // For a + b*sqrt(d), embeddings are a + b*sqrt(d) and a - b*sqrt(d)
    const a = element.__getitem__(0);
    const b = element.__getitem__(1);

    // Check norm > 0 and trace > 0
    const norm = element.norm();
    const trace = element.trace();

    return norm.numerator > 0n && trace.numerator > 0n;
  }

  // General case: not implemented
  throw new NotImplementedError('is_totally_positive not implemented for degree > 2');
}

/**
 * Compute the GCD of two number field elements.
 *
 * In a number field, gcd(a, b) = 1 for any non-zero elements
 * (since the field has no non-unit ideals).
 * For the ring of integers, this is more interesting.
 */
export function element_gcd(a: NumberFieldElement, b: NumberFieldElement): NumberFieldElement {
  // In a field, all non-zero elements are units, so gcd is 1
  if (a.is_zero() && b.is_zero()) {
    return a.parent().zero();
  }
  return a.parent().one();
}

/**
 * Compute the LCM of two number field elements.
 */
export function element_lcm(a: NumberFieldElement, b: NumberFieldElement): NumberFieldElement {
  // In a field, lcm(a, b) = a * b / gcd(a, b) = a * b
  return a.mul(b);
}
