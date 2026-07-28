/**
 * @module sage/rings/finite_rings/integer_mod
 * @description Elements of Z/nZ (integers modulo n)
 *
 * Port of: sage/rings/finite_rings/integer_mod.pyx
 */

import {
  crt,
  euler_phi,
  factor,
  gcd,
  is_prime,
  is_prime_power,
  lcm,
  power_mod,
  primitive_root,
  xgcd,
} from '../../arith/misc.js';
import { ArithmeticError, ValueError, ZeroDivisionError } from '../../errors.js';
import { discrete_log, has_order, order_from_multiple } from '../../groups/generic.js';
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
   * @param b - The base (default: `parent.multiplicative_generator()`)
   * @param order - The claimed order of `b`; only consulted when `check` is set
   *   (as in SageMath, where `order` is passed straight to `has_order`)
   * @param options.check - Verify that `b` really has order `order`
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
   *
   * @see Reference: sage/rings/finite_rings/integer_mod.pyx:log (lines 786-833)
   */
  log(b?: IntegerMod | bigint | number, order?: bigint, options?: { check?: boolean }): bigint {
    if (!this.isUnit()) {
      throw new ValueError(
        `logarithm of ${this.value} is not defined since it is not a unit modulo ${this.modulus}`
      );
    }

    // Convert base to IntegerMod if needed
    let base: IntegerMod;
    if (b === undefined) {
      base = new IntegerMod(multiplicative_generator(this.modulus), this.parent);
    } else {
      base = b instanceof IntegerMod ? new IntegerMod(b.value, this.parent) : new IntegerMod(b, this.parent);
      if (!base.isUnit()) {
        throw new ValueError(
          `logarithm with base ${base.value} is not defined since it is not a unit modulo ${this.modulus}`
        );
      }
    }

    if (options?.check) {
      if (order === undefined || !has_order(base, order, '*')) {
        throw new ValueError('base does not have the provided order');
      }
    }

    // Solve the DLP modulo every prime power dividing the modulus and combine
    // the answers with a *running* CRT, exactly as integer_mod.pyx:806-831:
    //
    //     n = crt(n, v, m, nb); m = lcm(m, nb)
    //
    // (the previous code kept only the first two components, so any modulus
    // with three or more prime factors returned a wrong exponent).
    let n = 0n;
    let m = 1n;

    for (const [p, e] of factor(this.modulus).filter(([q]) => q > 0n)) {
      const q = p ** e;
      const suffix = q !== this.modulus ? ` (no solution modulo ${q})` : '';
      const noLog = `no logarithm of ${this.value} found to base ${base.value} modulo ${this.modulus}`;

      const aRed = new IntegerMod(this.value, createParent(q));
      const bRed = new IntegerMod(base.value, createParent(q));

      const na = aRed.multiplicative_order();
      const nb = bRed.multiplicative_order();
      // Sage: `if not na.divides(nb)` -- self cannot be a power of b unless
      // ord(self) | ord(b).
      if (nb % na !== 0n) {
        throw new ValueError(noLog + suffix);
      }

      let v: bigint;
      try {
        v = discrete_log(aRed, bRed, nb, '*');
      } catch {
        throw new ValueError(noLog + suffix);
      }

      try {
        n = crt(n, v, m, nb);
      } catch {
        throw new ValueError(
          `no logarithm of ${this.value} found to base ${base.value} modulo ${this.modulus} (incompatible local solutions)`
        );
      }
      m = lcm(m, nb);
    }

    return n;
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
 * Return whether `(Z/nZ)*` is cyclic.
 *
 * Port of `sage/rings/finite_rings/integer_mod_ring.py:837-846`
 * (`IntegerModRing_generic.multiplicative_group_is_cyclic`): true exactly when
 * n < 8, or n is a power of an odd prime, or twice such a power.
 *
 * Lives here rather than in `integer_mod_ring.ts` so that `IntegerMod.log`
 * (whose parent may be the minimal `IntegerModRingBase`) can use it without an
 * import cycle; `IntegerModRing` re-exposes it as a method.
 */
export function multiplicative_group_is_cyclic(modulus: bigint): boolean {
  let n = modulus;
  if (n < 8n) {
    return true;
  }
  if (n % 4n === 0n) {
    return false; // n > 7, so the n = 4 case is not a problem
  }
  if (n % 4n === 2n) {
    n = n / 2n;
  }
  return is_prime_power(n);
}

/**
 * Generators of `(Z/nZ)*` together with their orders.
 *
 * Port of `integer_mod_ring.py:259-284` (`_unit_gens_primepowercase`) combined
 * with the CRT loop of `unit_gens` (`integer_mod_ring.py:1500-1510`).
 *
 * `sage: Integers(75).unit_gens()` -> `(26, 52)`;
 * `sage: Integers(162).unit_gens()` -> `(83,)`.
 */
export function unit_gens(modulus: bigint): Array<[bigint, bigint]> {
  if (modulus <= 1n) {
    return [];
  }
  const gens: Array<[bigint, bigint]> = [];
  for (const [p, e] of factor(modulus).filter(([q]) => q > 0n)) {
    const pr = p ** e;
    const m = modulus / pr;
    const local: Array<[bigint, bigint]> = [];
    if (p === 2n) {
      if (e === 1n) {
        // no generators
      } else if (e === 2n) {
        local.push([3n, 2n]);
      } else {
        local.push([pr - 1n, 2n]);
        local.push([5n, 2n ** (e - 2n)]);
      }
    } else {
      local.push([primitive_root(pr, false), p ** (e - 1n) * (p - 1n)]);
    }
    for (const [g, o] of local) {
      // Sage: `g.crt(Mod(1, m))` -- lift g to Z/nZ fixing 1 modulo the rest
      gens.push([m === 1n ? g : crt(g, 1n, pr, m), o]);
    }
  }
  return gens;
}

/**
 * A generator of `(Z/nZ)*`, assuming that group is cyclic.
 *
 * Port of `integer_mod_ring.py:848-895`
 * (`IntegerModRing_generic.multiplicative_generator`).  Replaces a scan over
 * all residues that computed a full multiplicative order for each one and
 * reported "no primitive root found modulo n" for non-cyclic groups.
 *
 * `sage: Integers(8).multiplicative_generator()` ->
 * `ValueError: multiplicative group of this ring is not cyclic`.
 */
export function multiplicative_generator(modulus: bigint): bigint {
  if (is_prime(modulus)) {
    // Sage: `self.field().multiplicative_generator()`, i.e. primitive_root(p)
    return primitive_root(modulus, false);
  }
  if (multiplicative_group_is_cyclic(modulus)) {
    const v = unit_gens(modulus);
    if (v.length !== 1) {
      // (Z/1Z)* and (Z/2Z)* are trivial: their generator is 1.
      if (v.length === 0) {
        return 1n;
      }
      throw new ArithmeticError(`expected one generator modulo ${modulus}, got ${v.length}`);
    }
    return v[0]![0];
  }
  throw new ValueError('multiplicative group of this ring is not cyclic');
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
