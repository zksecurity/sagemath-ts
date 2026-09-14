/**
 * @module sage/schemes/elliptic_curves/types
 * @description Type definitions for elliptic curves
 */

import {
  FiniteFieldExtension,
  PrimeField,
} from '../../rings/finite_rings/finite_field_extension.js';
import { FiniteFieldPrime } from '../../rings/finite_rings/finite_field_prime.js';
import { GF2Field } from '../../rings/finite_rings/gf2.js';

/** Parent equality for the field implementations used by the curve port.
 * @see Deviation: Generic Curve Isomorphism Ordering
 */
export function _same_base_ring(left: FieldRing, right: FieldRing): boolean {
  if (left === right) return true;
  if (left instanceof FiniteFieldExtension && right instanceof FiniteFieldExtension) {
    return (
      left.characteristic === right.characteristic &&
      left.degree === right.degree &&
      left.variableName === right.variableName &&
      left.modulus.eq(right.modulus)
    );
  }
  const isPrimeParent = (K: FieldRing) =>
    K instanceof PrimeField || K instanceof FiniteFieldPrime || K instanceof GF2Field;
  return (
    isPrimeParent(left) &&
    isPrimeParent(right) &&
    left.characteristic === right.characteristic &&
    (left instanceof PrimeField ? left.gen().value : 1n) ===
      (right instanceof PrimeField ? right.gen().value : 1n)
  );
}

/**
 * Interface for field elements that can be used in elliptic curve computations.
 *
 * This interface abstracts over different field implementations (prime fields,
 * extension fields, etc.) to allow elliptic curves to work with any field.
 */
export interface FieldElement {
  /** The parent ring/field */
  readonly parent: FieldRing;

  /** Add another element */
  add(other: FieldElement | number | bigint): FieldElement;

  /** Subtract another element */
  sub(other: FieldElement | number | bigint): FieldElement;

  /** Multiply by another element */
  mul(other: FieldElement | number | bigint): FieldElement;

  /** Divide by another element (must be non-zero) */
  div(other: FieldElement | number | bigint): FieldElement;

  /** Return the additive inverse */
  neg(): FieldElement;

  /** Return the multiplicative inverse (must be non-zero) */
  inv(): FieldElement;

  /** Raise to a power */
  pow(n: bigint | number): FieldElement;

  /** Check if this element is zero */
  isZero(): boolean;

  /** Check if this element equals another */
  eq(other: FieldElement): boolean;

  /** String representation */
  toString(): string;
}

/**
 * Interface for a ring/field that can be used as the base for elliptic curves.
 */
export interface FieldRing {
  /** Return the zero element */
  zero(): FieldElement;

  /** Return the one element */
  one(): FieldElement;

  /** Create an element from a value */
  __call__(value: bigint | number | FieldElement): FieldElement;

  /** The characteristic of the field */
  readonly characteristic: bigint;

  /** String representation */
  toString(): string;
}
