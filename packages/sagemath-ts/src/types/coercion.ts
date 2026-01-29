/**
 * @module sage/types/coercion
 * @description Type coercion utilities for flexible type acceptance
 *
 * This module provides the IntegerLike type and toBigInt function to enable
 * SageMath-style flexible type acceptance. In SageMath, functions accept both
 * Python int and Integer types; this module enables the same pattern in TypeScript
 * with bigint, number, and Integer types.
 *
 * Similarly, RationalLike and toRational provide flexible acceptance of values
 * that can be coerced to Rational numbers.
 */

import { Integer } from '../rings/integer_ring.js';
import { Rational } from '../rings/rational.js';

/**
 * Types that can be coerced to a bigint integer.
 * Mirrors SageMath's acceptance of both Python int and Integer.
 */
export type IntegerLike = bigint | number | Integer;

/**
 * Types that can be coerced to a Rational number.
 * Includes IntegerLike types (which become rationals with denominator 1)
 * and Rational numbers themselves.
 */
export type RationalLike = IntegerLike | Rational;

/**
 * Coerce an IntegerLike value to bigint.
 * @param x - Value to coerce
 * @returns The bigint value
 * @throws {TypeError} If x is a non-integer number or cannot be coerced
 */
export function toBigInt(x: IntegerLike): bigint {
  if (typeof x === 'bigint') return x;
  if (typeof x === 'number') {
    if (!Number.isInteger(x)) {
      throw new TypeError('cannot convert non-integer to Integer');
    }
    // Check for potential precision loss
    if (x > Number.MAX_SAFE_INTEGER || x < Number.MIN_SAFE_INTEGER) {
      throw new RangeError(
        `number ${x} exceeds safe integer range; use bigint literal (e.g., ${x}n) for large values`
      );
    }
    return BigInt(x);
  }
  if (x instanceof Integer) return x.value;
  throw new TypeError(`cannot coerce ${typeof x} to Integer`);
}

/**
 * Coerce a RationalLike value to Rational.
 * @param x - Value to coerce
 * @returns The Rational value
 * @throws {TypeError} If x is a non-integer number or cannot be coerced
 */
export function toRational(x: RationalLike): Rational {
  if (x instanceof Rational) return x;
  // For IntegerLike, convert to Rational with denominator 1
  const n = toBigInt(x);
  return new Rational(n, 1n);
}
