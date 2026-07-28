/**
 * @module sage/rings/finite_rings/finite_field_constructor
 * @description GF() constructor for finite fields
 *
 * Port of: sage/rings/finite_rings/finite_field_constructor.py
 */

import { is_prime, is_prime_power } from '../../arith/misc.js';
import { ArithmeticError, ValueError } from '../../errors.js';
import { GFExtended } from './finite_field_extension.js';
import { FiniteFieldElement, FiniteFieldPrime } from './finite_field_prime.js';

/**
 * Decompose q as p^n with p prime, or return `[q, 0n]` when q is not a prime
 * power.
 *
 * `arith.is_prime_power` is trusted for the decomposition but its answer is
 * re-verified here (`is_prime(p)` and `p^n === q`): PARI's `Z_factor` gives up
 * on composites past its trial-division bound and then reports them as prime,
 * which made e.g. a 13-digit semiprime look like a valid prime field order.
 */
function primePowerData(q: bigint): [bigint, bigint] {
  const [p, n] = is_prime_power(q, true);
  if (n === 0n || !is_prime(p) || p ** n !== q) {
    return [q, 0n];
  }
  return [p, n];
}

/**
 * Create a finite field.
 *
 * Currently supports only prime fields GF(p). Extension fields GF(p^n) for n > 1
 * will be supported in the future.
 *
 * @param order - The order of the field (must be a prime or prime power)
 * @param options - Optional configuration
 * @returns The finite field
 *
 * @example
 * ```typescript
 * // Prime field
 * const F7 = GF(7n);
 * const a = F7(3n);
 * const b = F7(5n);
 * console.log(a.mul(b)); // 1
 *
 * // Large prime field (Curve25519 base field)
 * const p = 2n ** 255n - 19n;
 * const Fp = GF(p);
 * ```
 *
 * @see Deviation: Finite Field Constructors and Display
 *
 * @throws {ValueError} If order is not a prime power
 * @throws {ValueError} If order is a prime power p^n with n > 1 (not yet supported)
 */
export function GF(order: bigint | number, options?: GFOptions): FiniteFieldPrime {
  const q = typeof order === 'number' ? BigInt(order) : order;
  const check = options?.check ?? true;

  if (q < 2n) {
    throw new ValueError('the order of a finite field must be at least 2');
  }

  // Check if order is a prime power
  if (check) {
    // primePowerData returns [p, k] where q = p^k, or [q, 0] if not a prime power
    const [p, n] = primePowerData(q);

    if (n === 0n) {
      throw new ValueError('the order of a finite field must be a prime power');
    }

    if (n > 1n) {
      throw new ValueError(
        `extension fields GF(${p}^${n}) are not yet supported; ` +
          'only prime fields GF(p) are currently available'
      );
    }
  }

  // For prime fields, create directly
  return new FiniteFieldPrime(q, check);
}

/**
 * Options for the GF constructor.
 */
export interface GFOptions {
  /**
   * Whether to check that the order is a prime power (default: true).
   * Set to false if you are certain the order is prime to skip the check.
   */
  check?: boolean;

  /**
   * Variable name for the generator (for extension fields).
   * Currently unused for prime fields.
   */
  name?: string;
}

/**
 * `FiniteField` is Sage's alias for `GF`
 * (`finite_field_constructor.py`: `FiniteField = FiniteFieldFactory("FiniteField")`,
 * and `GF = FiniteField`), so it must accept extension orders as well.
 *
 * The local `GF` above is the prime-field-only constructor returning the
 * narrowly typed `FiniteFieldPrime`; many call sites depend on that type, so
 * the general constructor is `GFExtended` from `finite_field_extension.ts`,
 * which is what the package exports as `GF`.
 */
export const FiniteField = GFExtended;

/**
 * Check if a value is a valid finite field order.
 *
 * @param q - The value to check
 * @returns True if q is a prime power >= 2
 */
export function isValidFiniteFieldOrder(q: bigint | number): boolean {
  const order = typeof q === 'number' ? BigInt(q) : q;

  if (order < 2n) {
    return false;
  }

  // primePowerData returns [p, k] where k=0 means not a prime power
  const [_p, k] = primePowerData(order);
  return k > 0n;
}

/**
 * Get information about a finite field order.
 *
 * @param q - The order to analyze
 * @returns Information about the field, or null if not a valid order
 */
export function analyzeFiniteFieldOrder(q: bigint | number): FiniteFieldInfo | null {
  const order = typeof q === 'number' ? BigInt(q) : q;

  if (order < 2n) {
    return null;
  }

  // primePowerData returns [p, k] where q = p^k, or [q, 0] if not a prime power
  const [p, n] = primePowerData(order);

  if (n === 0n) {
    return null;
  }

  return {
    order,
    characteristic: p,
    degree: n,
    isPrimeField: n === 1n,
  };
}

/**
 * Information about a finite field order.
 */
export interface FiniteFieldInfo {
  /** The total order q = p^n */
  order: bigint;
  /** The characteristic p */
  characteristic: bigint;
  /** The extension degree n */
  degree: bigint;
  /** Whether this is a prime field (n = 1) */
  isPrimeField: boolean;
}

// Re-export types for convenience
export { FiniteFieldPrime, FiniteFieldElement } from './finite_field_prime.js';
