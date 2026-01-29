/**
 * @module sage/rings/finite_rings/integer_mod
 * @description Elements of Z/nZ (integers modulo n)
 *
 * Port of: sage/rings/finite_rings/integer_mod.pyx
 */

import { crt, euler_phi, factor, gcd, power_mod, xgcd } from '../../arith/misc.js';
import { ValueError, ZeroDivisionError } from '../../errors.js';
import { discrete_log, order_from_multiple } from '../../groups/generic.js';
import type { RingElement } from '../polynomial/polynomial_element.js';

/**
 * Forward declaration for parent ring type.
 */
export interface IntegerModRingBase {
  readonly modulus: bigint;
  zero(): IntegerMod;
  one(): IntegerMod;
  __call__(x: unknown): IntegerMod;
  is_field?(): boolean;
}

/**
 * An element of Z/nZ.
 *
 * Elements are represented as integers in the range [0, n).
 *
 * @example
 * ```typescript
 * const Zmod5 = Zmod(5n);
 * const a = Zmod5(3n);
 * const b = Zmod5(4n);
 * console.log(a.add(b).value); // 2n (since 3+4=7 ≡ 2 mod 5)
 * ```
 */
export class IntegerMod implements RingElement {
  readonly value: bigint;
  readonly parent: IntegerModRingBase;

  /**
   * Create an element of Z/nZ.
   *
   * @param value - The integer value (will be reduced modulo n)
   * @param parent - The parent ring Z/nZ
   */
  constructor(value: bigint | number | IntegerMod, parent: IntegerModRingBase) {
    this.parent = parent;

    if (value instanceof IntegerMod) {
      // Coerce from another IntegerMod
      this.value = mod(value.value, parent.modulus);
    } else {
      const v = typeof value === 'number' ? BigInt(value) : value;
      this.value = mod(v, parent.modulus);
    }
  }

  /**
   * Return the modulus of the parent ring.
   */
  get modulus(): bigint {
    return this.parent.modulus;
  }

  /**
   * Add two elements.
   */
  add(other: IntegerMod | number | bigint): IntegerMod {
    const otherVal = this.coerceValue(other);
    return new IntegerMod(mod(this.value + otherVal, this.modulus), this.parent);
  }

  /**
   * Subtract two elements.
   */
  sub(other: IntegerMod | number | bigint): IntegerMod {
    const otherVal = this.coerceValue(other);
    return new IntegerMod(mod(this.value - otherVal, this.modulus), this.parent);
  }

  /**
   * Multiply two elements.
   */
  mul(other: IntegerMod | number | bigint): IntegerMod {
    const otherVal = this.coerceValue(other);
    return new IntegerMod(mod(this.value * otherVal, this.modulus), this.parent);
  }

  /**
   * Divide two elements.
   * This requires the divisor to be invertible (gcd(divisor, n) = 1).
   *
   * @throws {ZeroDivisionError} If the divisor is not invertible
   */
  div(other: IntegerMod | number | bigint): IntegerMod {
    const otherVal = this.coerceValue(other);
    if (otherVal === 0n) {
      throw new ZeroDivisionError('division by zero in Z/nZ');
    }

    const [g, s] = xgcd(otherVal, this.modulus);
    if (g !== 1n) {
      throw new ZeroDivisionError(`inverse of Mod(${otherVal}, ${this.modulus}) does not exist`);
    }

    return new IntegerMod(mod(this.value * s, this.modulus), this.parent);
  }

  /**
   * Return the additive inverse (-self).
   */
  neg(): IntegerMod {
    if (this.value === 0n) {
      return this;
    }
    return new IntegerMod(this.modulus - this.value, this.parent);
  }

  /**
   * Return the multiplicative inverse (1/self).
   *
   * @throws {ZeroDivisionError} If self is not invertible
   */
  inv(): IntegerMod {
    if (this.value === 0n) {
      throw new ZeroDivisionError('division by zero in Z/nZ');
    }

    const [g, s] = xgcd(this.value, this.modulus);
    if (g !== 1n) {
      throw new ZeroDivisionError(`inverse of Mod(${this.value}, ${this.modulus}) does not exist`);
    }

    return new IntegerMod(mod(s, this.modulus), this.parent);
  }

  /**
   * Return self^n.
   *
   * @param n - The exponent (can be negative if self is invertible)
   */
  pow(n: number | bigint): IntegerMod {
    const exp = typeof n === 'number' ? BigInt(n) : n;

    if (exp === 0n) {
      return this.parent.one();
    }

    if (exp < 0n) {
      // For negative exponents, compute inverse first
      return this.inv().pow(-exp);
    }

    // Use the optimized modular exponentiation
    const result = power_mod(this.value, exp, this.modulus);
    return new IntegerMod(result, this.parent);
  }

  /**
   * Check equality with another element.
   */
  eq(other: IntegerMod | number | bigint): boolean {
    const otherVal = this.coerceValue(other);
    return this.value === otherVal;
  }

  /**
   * Check if this element is zero.
   */
  isZero(): boolean {
    return this.value === 0n;
  }

  /**
   * Check if this element is one.
   */
  isOne(): boolean {
    return this.value === 1n;
  }

  /**
   * Check if this element is a unit (invertible).
   */
  isUnit(): boolean {
    return gcd(this.value, this.modulus) === 1n;
  }

  /**
   * Lift this element to an integer.
   */
  lift(): bigint {
    return this.value;
  }

  /**
   * Return the integer value as a bigint.
   */
  toBigInt(): bigint {
    return this.value;
  }

  /**
   * String representation.
   */
  toString(): string {
    return this.value.toString();
  }

  /**
   * Repr for debugging.
   */
  repr(): string {
    return `Mod(${this.value}, ${this.modulus})`;
  }

  /**
   * Return the multiplicative order of this element.
   *
   * The multiplicative order is the smallest positive integer k such that
   * self^k = 1 (mod n).
   *
   * @returns The multiplicative order
   * @throws {ValueError} If the element is not a unit
   *
   * @example
   * ```typescript
   * Mod(2n, 7n).multiplicative_order(); // 3n (since 2^3 = 8 = 1 mod 7)
   * Mod(3n, 7n).multiplicative_order(); // 6n (primitive root)
   * ```
   */
  multiplicative_order(): bigint {
    if (!this.isUnit()) {
      throw new ValueError(
        `multiplicative order of ${this.value} not defined since it is not a unit modulo ${this.modulus}`
      );
    }

    if (this.value === 1n) {
      return 1n;
    }

    // The order divides phi(n), use order_from_multiple for efficiency
    const phi = euler_phi(this.modulus);
    return order_from_multiple(this, phi, undefined, '*');
  }

  /**
   * Return the discrete logarithm of self with respect to base b.
   *
   * Find x such that b^x = self (mod n).
   *
   * @param b - The base (default: a generator of the multiplicative group if available)
   * @param order - The order of b (computed if not provided)
   * @returns x such that b^x = self
   * @throws {ValueError} If no such x exists or if self/b is not a unit
   *
   * @example
   * ```typescript
   * const a = Mod(5n, 37n);
   * const b = Mod(2n, 37n);
   * // If 5 = 2^x mod 37, find x
   * const x = a.log(b); // Find x such that 2^x = 5 mod 37
   * ```
   */
  log(b?: IntegerMod | bigint | number, order?: bigint): bigint {
    if (!this.isUnit()) {
      throw new ValueError(
        `logarithm of ${this.value} is not defined since it is not a unit modulo ${this.modulus}`
      );
    }

    // Convert base to IntegerMod if needed
    let base: IntegerMod;
    if (b === undefined) {
      // Use a primitive root if we can find one (for prime moduli)
      base = this._findPrimitiveRoot();
    } else if (b instanceof IntegerMod) {
      base = new IntegerMod(b.value, this.parent);
    } else {
      base = new IntegerMod(b, this.parent);
    }

    if (!base.isUnit()) {
      throw new ValueError(
        `logarithm with base ${base.value} is not defined since it is not a unit modulo ${this.modulus}`
      );
    }

    // Use the generic discrete_log algorithm
    // For composite moduli, we use CRT to combine results
    const factors = factor(this.modulus);

    // Filter out -1 factor
    const primeFactors = factors.filter(([p]) => p > 0n);

    if (primeFactors.length === 1) {
      // Prime power modulus - use direct discrete log
      const baseOrder = order ?? base.multiplicative_order();
      return discrete_log(this, base, baseOrder, '*');
    }

    // Composite modulus - use CRT
    const residues: bigint[] = [];
    const moduli: bigint[] = [];

    for (const [p, e] of primeFactors) {
      const q = p ** e;
      const selfReduced = new IntegerMod(this.value, createParent(q));
      const baseReduced = new IntegerMod(base.value, createParent(q));

      // Check if solution exists locally
      const selfOrder = selfReduced.multiplicative_order();
      const baseOrder = baseReduced.multiplicative_order();

      if (selfOrder % gcd(selfOrder, baseOrder) !== 0n) {
        throw new ValueError(
          `no logarithm of ${this.value} found to base ${base.value} modulo ${this.modulus} (no solution modulo ${q})`
        );
      }

      // Solve locally
      try {
        const localLog = discrete_log(selfReduced, baseReduced, baseOrder, '*');
        residues.push(localLog);
        moduli.push(baseOrder);
      } catch {
        throw new ValueError(
          `no logarithm of ${this.value} found to base ${base.value} modulo ${this.modulus} (no solution modulo ${q})`
        );
      }
    }

    // Combine using CRT
    try {
      return crt(residues[0]!, residues[1] ?? 0n, moduli[0]!, moduli[1] ?? 1n);
    } catch {
      throw new ValueError(
        `no logarithm of ${this.value} found to base ${base.value} modulo ${this.modulus} (incompatible local solutions)`
      );
    }
  }

  /**
   * Find a primitive root modulo n (if one exists).
   * @private
   */
  private _findPrimitiveRoot(): IntegerMod {
    const phi = euler_phi(this.modulus);

    // Try small values first
    for (let g = 2n; g < this.modulus; g++) {
      const elem = new IntegerMod(g, this.parent);
      if (elem.isUnit()) {
        const ord = elem.multiplicative_order();
        if (ord === phi) {
          return elem;
        }
      }
    }

    throw new ValueError(`no primitive root found modulo ${this.modulus}`);
  }

  /**
   * Coerce a value to a bigint in [0, modulus).
   */
  private coerceValue(other: IntegerMod | number | bigint): bigint {
    if (other instanceof IntegerMod) {
      return mod(other.value, this.modulus);
    }
    const v = typeof other === 'number' ? BigInt(other) : other;
    return mod(v, this.modulus);
  }
}

/**
 * Compute a mod n, ensuring the result is in [0, n).
 */
function mod(a: bigint, n: bigint): bigint {
  const result = a % n;
  return result < 0n ? result + n : result;
}

/**
 * Create a minimal parent ring for a given modulus.
 * @private
 */
function createParent(modulus: bigint): IntegerModRingBase {
  return {
    modulus,
    zero(): IntegerMod {
      return new IntegerMod(0n, this);
    },
    one(): IntegerMod {
      return new IntegerMod(1n, this);
    },
    __call__(x: unknown): IntegerMod {
      if (typeof x === 'number' || typeof x === 'bigint') {
        return new IntegerMod(x, this);
      }
      if (x instanceof IntegerMod) {
        return new IntegerMod(x.value, this);
      }
      throw new ValueError(`cannot coerce ${x} to IntegerMod`);
    },
  };
}

/**
 * Create an IntegerMod element. Shorthand for the SageMath Mod() function.
 *
 * @param value - The integer value
 * @param modulus - The modulus
 * @returns An IntegerMod element
 *
 * @example
 * ```typescript
 * const a = Mod(3n, 7n);  // 3 mod 7
 * console.log(a.inv());   // 5 (since 3*5 = 15 ≡ 1 mod 7)
 * ```
 */
export function Mod(value: bigint | number, modulus: bigint | number): IntegerMod {
  const m = typeof modulus === 'number' ? BigInt(modulus) : modulus;
  if (m <= 0n) {
    throw new ValueError('modulus must be positive');
  }

  // Create a temporary parent ring
  const parent: IntegerModRingBase = {
    modulus: m,
    zero(): IntegerMod {
      return new IntegerMod(0n, this);
    },
    one(): IntegerMod {
      return new IntegerMod(1n, this);
    },
    __call__(x: unknown): IntegerMod {
      if (typeof x === 'number' || typeof x === 'bigint') {
        return new IntegerMod(x, this);
      }
      if (x instanceof IntegerMod) {
        return new IntegerMod(x.value, this);
      }
      throw new ValueError(`cannot coerce ${x} to IntegerMod`);
    },
  };

  return new IntegerMod(value, parent);
}
