/**
 * @module sage/types/coercion
 * @description Type coercion utilities for flexible type acceptance
 *
 * This module provides the IntegerLike type and toBigInt function to enable
 * SageMath-style flexible type acceptance. In SageMath, functions accept both
 * Python int and Integer types; this module enables the same pattern in TypeScript
 * with bigint and Integer types.
 *
 * IMPORTANT: Unlike SageMath, we do NOT accept JavaScript numbers because they
 * are IEEE 754 floats that silently lose precision beyond 2^53-1. This is
 * critical for cryptographic correctness. Use bigint literals (e.g., 123n).
 *
 * @see Deviation: no-number-coercion
 */

import { Integer } from '../rings/integer_ring.js';
import { Rational } from '../rings/rational.js';

/**
 * Types that can be coerced to a bigint integer.
 *
 * NOTE: JavaScript `number` is intentionally excluded because it silently
 * loses precision for values > 2^53-1. Use bigint literals (e.g., 123n).
 *
 * @see Deviation: no-number-coercion
 */
export type IntegerLike = bigint | Integer;

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
 * @throws {TypeError} If x is a number (use bigint literals instead) or cannot be coerced
 */
export function toBigInt(x: IntegerLike): bigint {
  if (typeof x === 'bigint') return x;
  if (typeof x === 'number') {
    throw new TypeError(
      'JavaScript numbers are not accepted due to precision loss risk; use bigint literals (e.g., 123n)'
    );
  }
  if (x instanceof Integer) return x.value;
  throw new TypeError(`cannot coerce ${typeof x} to Integer`);
}

/**
 * Safely convert a bigint to a JavaScript number.
 *
 * Use this instead of Number(bigint) to prevent silent precision loss.
 * For cases where precision loss is intentional (e.g., floating-point math),
 * use Number() directly but document the reason.
 *
 * @param x - The bigint to convert
 * @returns The number value
 * @throws {RangeError} If x exceeds the safe integer range (±2^53-1)
 */
export function toSafeNumber(x: bigint): number {
  if (x > BigInt(Number.MAX_SAFE_INTEGER) || x < BigInt(Number.MIN_SAFE_INTEGER)) {
    throw new RangeError(
      `bigint ${x} exceeds safe integer range (±${Number.MAX_SAFE_INTEGER}); cannot safely convert to number`
    );
  }
  return Number(x);
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
