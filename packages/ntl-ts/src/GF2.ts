/**
 * NTL binary field GF(2).
 * @see Reference: ntl/src/GF2.cpp
 *
 * GF2 represents elements of the binary field {0, 1} with
 * arithmetic modulo 2.
 */

/**
 * Binary field element.
 * TypeScript port of NTL's GF2 class.
 *
 * Elements are 0 or 1, with XOR as addition and AND as multiplication.
 */
export class GF2 {
  private _value: 0 | 1;

  /**
   * Creates a new GF2 element.
   * @param value - Initial value (0 or 1)
   */
  constructor(value?: number | boolean) {
    if (value === undefined) {
      this._value = 0;
    } else if (typeof value === 'boolean') {
      this._value = value ? 1 : 0;
    } else {
      // NTL: GF2(INIT_VAL, a) stores a & 1.
      this._value = (Number(value) & 1) as 0 | 1;
    }
  }

  // ============================================
  // Constants
  // ============================================

  /**
   * Returns the zero element.
   */
  static zero(): GF2 {
    return new GF2(0);
  }

  /**
   * Returns the one element.
   */
  static one(): GF2 {
    return new GF2(1);
  }

  // ============================================
  // Arithmetic Operations
  // ============================================

  /**
   * Adds two GF2 elements (XOR).
   * @param other - Element to add
   * @returns The sum
   */
  add(other: GF2): GF2 {
    return new GF2(this._value ^ other._value);
  }

  /**
   * Subtracts a GF2 element (same as add in characteristic 2).
   * @param other - Element to subtract
   * @returns The difference
   */
  sub(other: GF2): GF2 {
    return new GF2(this._value ^ other._value);
  }

  /**
   * Multiplies two GF2 elements (AND).
   * @param other - Element to multiply
   * @returns The product
   */
  mul(other: GF2): GF2 {
    return new GF2(this._value & other._value);
  }

  /**
   * Divides by a GF2 element.
   * @param other - Divisor (must be 1)
   * @returns The quotient
   * @throws ArithmeticError "GF2: division by zero" if `other` is 0
   */
  div(other: GF2): GF2 {
    // NTL GF2.h:220 -- operator/(GF2 a, GF2 b)
    if (other.IsZero()) throw new Error('GF2: division by zero');
    return new GF2(this._value);
  }

  /**
   * Negates the element (identity in characteristic 2).
   * @returns The negation
   */
  negate(): GF2 {
    return new GF2(this._value);
  }

  /**
   * Computes multiplicative inverse.
   * @returns The inverse (must be 1)
   * @throws ArithmeticError "GF2: division by zero" if this is 0
   */
  inv(): GF2 {
    // NTL GF2.h:234 -- inv(a) == 1/a
    return GF2.one().div(this);
  }

  /**
   * Computes a^e.
   * @param e - Exponent
   * @returns a^e
   */
  power(e: number | bigint): GF2 {
    // NTL GF2.cpp:8 -- power(a, 0) == 1, otherwise a (0^negative is an error)
    const _e = BigInt(e);
    if (_e === 0n) return GF2.one();
    if (_e < 0n && this.IsZero()) throw new Error('GF2: division by zero');
    return new GF2(this._value);
  }

  /**
   * Squares the element.
   * @returns The square (identity in GF2)
   */
  sqr(): GF2 {
    return new GF2(this._value);
  }

  // ============================================
  // Comparison
  // ============================================

  /**
   * Checks if element is zero.
   * @returns True if zero
   */
  IsZero(): boolean {
    return this._value === 0;
  }

  /**
   * Checks if element is one.
   * @returns True if one
   */
  IsOne(): boolean {
    return this._value === 1;
  }

  /**
   * Checks equality.
   * @param other - Element to compare
   * @returns True if equal
   */
  equals(other: GF2): boolean {
    return this._value === other._value;
  }

  // ============================================
  // Conversion
  // ============================================

  /**
   * Returns the representative (0 or 1).
   * @returns 0 or 1
   */
  rep(): number {
    return this._value;
  }

  /**
   * Converts to number.
   * @returns 0 or 1
   */
  toNumber(): number {
    return this._value;
  }

  /**
   * Converts to boolean.
   * @returns False if 0, true if 1
   */
  toBoolean(): boolean {
    return this._value === 1;
  }

  /**
   * Converts to string.
   * @returns "0" or "1"
   */
  toString(): string {
    return this._value === 0 ? '0' : '1';
  }

  /**
   * Converts from various types.
   * @param value - Value to convert
   * @returns GF2 element
   */
  static conv(value: number | bigint | boolean): GF2 {
    // NTL GF2.h:172/176 -- to_GF2(long a) keeps a & 1, to_GF2(ZZ a) keeps IsOdd(a)
    if (typeof value === 'bigint') {
      return new GF2(value % 2n === 0n ? 0 : 1);
    }
    return new GF2(value as number | boolean);
  }

  // ============================================
  // Random
  // ============================================

  /**
   * Generates random GF2 element.
   * @returns Random element (0 or 1)
   */
  static random(): GF2 {
    throw new Error('NTL_NOT_IMPLEMENTED: GF2.random');
  }
}

// ============================================
// Standalone Functions (NTL-style API)
// ============================================

/**
 * Converts to GF2.
 */
export function to_GF2(value: number | bigint | boolean): GF2 {
  return GF2.conv(value);
}

/**
 * Returns representative.
 */
export function rep(a: GF2): number {
  return a.rep();
}

/**
 * Adds two GF2 elements.
 */
export function add(a: GF2, b: GF2): GF2 {
  return a.add(b);
}

/**
 * Subtracts two GF2 elements.
 */
export function sub(a: GF2, b: GF2): GF2 {
  return a.sub(b);
}

/**
 * Multiplies two GF2 elements.
 */
export function mul(a: GF2, b: GF2): GF2 {
  return a.mul(b);
}

/**
 * Divides two GF2 elements.
 */
export function div(a: GF2, b: GF2): GF2 {
  return a.div(b);
}

/**
 * Computes power.
 */
export function power(a: GF2, e: number | bigint): GF2 {
  return a.power(e);
}

/**
 * Checks if zero.
 */
export function IsZero(a: GF2): boolean {
  return a.IsZero();
}

/**
 * Checks if one.
 */
export function IsOne(a: GF2): boolean {
  return a.IsOne();
}

/**
 * Computes inverse.
 */
export function inv(a: GF2): GF2 {
  return a.inv();
}
