/**
 * @module sage/rings/function_field/constant_field
 * @description Structural interface for the constant field of a function field.
 *
 * SageMath has no such file: there, any object in ``Fields()`` can serve as the
 * constant field of a rational function field and the coercion framework sorts
 * out the rest.  TypeScript has no coercion framework, so this module pins down
 * (structurally) exactly what a constant field has to provide, and supplies the
 * few accessors that have to tolerate the fact that our own field classes spell
 * ``characteristic`` sometimes as a property (``FiniteFieldPrime``,
 * ``FiniteFieldExtension``) and sometimes as a method (``RationalField``).
 *
 * @see DESIGN.md - dependency architecture
 */

import type { CoefficientRing, RingElement } from '../polynomial/polynomial_element.js';

/**
 * An element of a constant field.  On top of {@link RingElement} we need
 * inversion (the constant field is a field) and a total order (SageMath
 * compares polynomials coefficient-wise, and places are sorted).
 */
export interface ConstantFieldElement extends RingElement {
  inv(): this;
  isOne?(): boolean;
}

/**
 * A field usable as the constant field of a rational function field.
 *
 * Satisfied structurally by `FiniteFieldPrime`, `PrimeField`,
 * `FiniteFieldExtension` and `RationalField`.
 */
export interface ConstantField<C extends ConstantFieldElement> extends CoefficientRing<C> {
  zero(): C;
  one(): C;
  __call__(x: unknown): C;
  is_field?(): boolean;
  /** Present as a property on our finite fields, as a method on `RationalField`. */
  characteristic?: bigint | (() => bigint);
  cardinality?(): bigint;
  order?: bigint | (() => unknown);
  is_finite?(): boolean;
  elements?(): IterableIterator<C>;
  list?(): C[];
  degree?: number | (() => bigint);
  toString(): string;
}

/**
 * Return the characteristic of ``k`` as a bigint.
 *
 * Port of ``FunctionField.characteristic``
 * (`reference/sage/src/sage/rings/function_field/function_field.py:389`), which
 * is just ``self.constant_base_field().characteristic()``.
 */
export function constant_field_characteristic<C extends ConstantFieldElement>(
  k: ConstantField<C>
): bigint {
  const c = (k as { characteristic?: unknown }).characteristic;
  if (typeof c === 'bigint') {
    return c;
  }
  if (typeof c === 'function') {
    return BigInt((c as () => bigint).call(k));
  }
  if (typeof c === 'number') {
    return BigInt(c);
  }
  throw new TypeError(`cannot determine the characteristic of ${k}`);
}

/**
 * Return whether ``k`` is a finite field.
 *
 * Used by ``FunctionField.is_global``
 * (`reference/sage/src/sage/rings/function_field/function_field.py:426`).
 */
export function constant_field_is_finite<C extends ConstantFieldElement>(
  k: ConstantField<C>
): boolean {
  if (typeof k.is_finite === 'function') {
    return k.is_finite();
  }
  // Our finite fields expose ``cardinality()``; ``RationalField`` exposes
  // ``is_finite()`` returning false, so it never reaches here.
  return typeof k.cardinality === 'function';
}

/**
 * Return the cardinality of a finite constant field.
 */
export function constant_field_cardinality<C extends ConstantFieldElement>(
  k: ConstantField<C>
): bigint {
  if (typeof k.cardinality === 'function') {
    return k.cardinality();
  }
  const o = (k as { order?: unknown }).order;
  if (typeof o === 'bigint') {
    return o;
  }
  throw new TypeError(`${k} is not a finite field`);
}

/**
 * Return the elements of a finite constant field, in the field's own iteration
 * order.
 *
 * SageMath's ``PolynomialRing._polys_max`` iterates ``for c in base``
 * (`reference/sage/src/sage/rings/polynomial/polynomial_ring.py:1548`), so the
 * enumeration order of places of a given degree is inherited from the constant
 * field's iteration order.  For prime fields both SageMath and we iterate
 * `0, 1, ..., p-1`.
 */
export function constant_field_element_list<C extends ConstantFieldElement>(
  k: ConstantField<C>
): C[] {
  if (typeof k.list === 'function') {
    return k.list();
  }
  if (typeof k.elements === 'function') {
    return [...k.elements()];
  }
  throw new TypeError(`cannot enumerate the elements of ${k}`);
}

/**
 * Total order on constant field elements, mirroring SageMath's ``richcmp``.
 *
 * For prime fields SageMath compares the canonical lifts `0 <= a < p`; this is
 * what ``lift()``/``toBigInt()``/``integer_representation()`` give us.  For any
 * other field we fall back to the string representation, which orders
 * deterministically but need not agree with SageMath.
 *
 * @see Deviation: constant-field element ordering outside prime fields
 */
export function compare_constants<C extends ConstantFieldElement>(a: C, b: C): number {
  const av = constantSortKey(a);
  const bv = constantSortKey(b);
  if (av !== null && bv !== null) {
    return av < bv ? -1 : av > bv ? 1 : 0;
  }
  const as = a.toString();
  const bs = b.toString();
  return as < bs ? -1 : as > bs ? 1 : 0;
}

function constantSortKey(a: unknown): bigint | null {
  const o = a as {
    lift?: unknown;
    toBigInt?: () => bigint;
    integer_representation?: () => bigint;
    numerator?: () => bigint;
  };
  if (typeof o.toBigInt === 'function') {
    try {
      return o.toBigInt();
    } catch {
      /* fall through */
    }
  }
  if (typeof o.lift === 'function') {
    const v = (o.lift as () => unknown)();
    if (typeof v === 'bigint') {
      return v;
    }
  } else if (typeof o.lift === 'bigint') {
    return o.lift;
  }
  if (typeof o.integer_representation === 'function') {
    return o.integer_representation();
  }
  return null;
}

/**
 * Divide two constant field elements.
 */
export function divide_constants<C extends ConstantFieldElement>(a: C, b: C): C {
  const withDiv = a as unknown as { div?: (o: C) => C };
  if (typeof withDiv.div === 'function') {
    return withDiv.div(b);
  }
  return a.mul(b.inv()) as C;
}
