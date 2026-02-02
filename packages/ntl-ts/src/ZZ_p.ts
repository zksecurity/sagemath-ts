/**
 * NTL integers mod p (ZZ_p).
 * @see Reference: ntl/src/ZZ_p.cpp
 *
 * ZZ_p represents integers modulo a prime p. The modulus is set
 * globally using ZZ_p.init(p) before performing operations.
 */

import type { ZZ } from './ZZ.js';

/**
 * Context for ZZ_p modulus.
 * Stores the modulus and related precomputed values.
 */
export class ZZ_pContext {
  private _modulus: ZZ;

  /**
   * Creates a context with the given modulus.
   * @param p - The modulus (must be > 1)
   */
  constructor(p: ZZ) {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_pContext');
  }

  /**
   * Saves the current context.
   */
  save(): void {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_pContext.save');
  }

  /**
   * Restores this context as the current modulus.
   */
  restore(): void {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_pContext.restore');
  }
}

/**
 * Backup class for saving/restoring ZZ_p modulus.
 */
export class ZZ_pBak {
  private _context: ZZ_pContext | null = null;
  private _mustRestore: boolean = false;

  /**
   * Saves the current modulus.
   */
  save(): void {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_pBak.save');
  }

  /**
   * Restores the saved modulus.
   */
  restore(): void {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_pBak.restore');
  }
}

/**
 * Integers mod p.
 * TypeScript port of NTL's ZZ_p class.
 *
 * IMPORTANT: Call ZZ_p.init(p) to set the modulus before operations.
 */
export class ZZ_p {
  private _rep: ZZ;

  /**
   * Creates a new ZZ_p element.
   * @param value - Initial value (will be reduced mod p)
   */
  constructor(value?: ZZ | bigint | number) {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_p');
  }

  // ============================================
  // Modulus Management
  // ============================================

  /**
   * Sets the global modulus for ZZ_p operations.
   * @param p - The modulus (must be > 1)
   */
  static init(p: ZZ): void {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_p.init');
  }

  /**
   * Returns the current modulus.
   * @returns The modulus
   */
  static modulus(): ZZ {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_p.modulus');
  }

  /**
   * Returns the size of the modulus in words.
   * @returns Modulus size
   */
  static ModulusSize(): number {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_p.ModulusSize');
  }

  // ============================================
  // Element Access
  // ============================================

  /**
   * Returns the zero element.
   */
  static zero(): ZZ_p {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_p.zero');
  }

  /**
   * Returns the one element.
   */
  static one(): ZZ_p {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_p.one');
  }

  /**
   * Returns the representative (reduced mod p).
   * @returns The ZZ representative
   */
  rep(): ZZ {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_p.rep');
  }

  // ============================================
  // Arithmetic Operations
  // ============================================

  /**
   * Adds two ZZ_p values.
   * @param other - Value to add
   * @returns The sum
   */
  add(other: ZZ_p): ZZ_p {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_p.add');
  }

  /**
   * Subtracts a ZZ_p value.
   * @param other - Value to subtract
   * @returns The difference
   */
  sub(other: ZZ_p): ZZ_p {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_p.sub');
  }

  /**
   * Multiplies two ZZ_p values.
   * @param other - Value to multiply
   * @returns The product
   */
  mul(other: ZZ_p): ZZ_p {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_p.mul');
  }

  /**
   * Divides by a ZZ_p value.
   * @param other - Divisor
   * @returns The quotient
   */
  div(other: ZZ_p): ZZ_p {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_p.div');
  }

  /**
   * Negates the value.
   * @returns The negation
   */
  negate(): ZZ_p {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_p.negate');
  }

  /**
   * Squares the value.
   * @returns The square
   */
  sqr(): ZZ_p {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_p.sqr');
  }

  /**
   * Computes multiplicative inverse.
   * @returns The inverse
   * @throws If element is not invertible
   */
  inv(): ZZ_p {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_p.inv');
  }

  /**
   * Computes a^e (exponentiation).
   * @param e - Exponent
   * @returns a^e
   */
  power(e: ZZ | bigint | number): ZZ_p {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_p.power');
  }

  // ============================================
  // Comparison
  // ============================================

  /**
   * Checks if value is zero.
   * @returns True if zero
   */
  IsZero(): boolean {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_p.IsZero');
  }

  /**
   * Checks if value is one.
   * @returns True if one
   */
  IsOne(): boolean {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_p.IsOne');
  }

  /**
   * Checks equality.
   * @param other - Value to compare
   * @returns True if equal
   */
  equals(other: ZZ_p): boolean {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_p.equals');
  }

  // ============================================
  // Conversion
  // ============================================

  /**
   * Converts a ZZ to ZZ_p.
   * @param a - Value to convert
   * @returns The ZZ_p element
   */
  static conv(a: ZZ | bigint | number): ZZ_p {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_p.conv');
  }

  /**
   * Converts to string.
   * @returns String representation
   */
  toString(): string {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_p.toString');
  }

  // ============================================
  // Random
  // ============================================

  /**
   * Generates a random ZZ_p element.
   * @returns Random element
   */
  static random(): ZZ_p {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_p.random');
  }
}

// ============================================
// Standalone Functions (NTL-style API)
// ============================================

/**
 * Adds two ZZ_p values.
 */
export function add(a: ZZ_p, b: ZZ_p): ZZ_p {
  return a.add(b);
}

/**
 * Subtracts two ZZ_p values.
 */
export function sub(a: ZZ_p, b: ZZ_p): ZZ_p {
  return a.sub(b);
}

/**
 * Multiplies two ZZ_p values.
 */
export function mul(a: ZZ_p, b: ZZ_p): ZZ_p {
  return a.mul(b);
}

/**
 * Divides two ZZ_p values.
 */
export function div(a: ZZ_p, b: ZZ_p): ZZ_p {
  return a.div(b);
}

/**
 * Computes inverse.
 */
export function inv(a: ZZ_p): ZZ_p {
  return a.inv();
}

/**
 * Negates a ZZ_p value.
 */
export function negate(a: ZZ_p): ZZ_p {
  return a.negate();
}

/**
 * Squares a ZZ_p value.
 */
export function sqr(a: ZZ_p): ZZ_p {
  return a.sqr();
}

/**
 * Computes power.
 */
export function power(a: ZZ_p, e: ZZ | bigint | number): ZZ_p {
  return a.power(e);
}

/**
 * Checks if zero.
 */
export function IsZero(a: ZZ_p): boolean {
  return a.IsZero();
}

/**
 * Checks if one.
 */
export function IsOne(a: ZZ_p): boolean {
  return a.IsOne();
}

/**
 * Converts to ZZ_p.
 */
export function conv(a: ZZ | bigint | number): ZZ_p {
  return ZZ_p.conv(a);
}
