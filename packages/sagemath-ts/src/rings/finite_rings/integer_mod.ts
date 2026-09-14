/**
 * @module sage/rings/finite_rings/integer_mod
 * @description Elements of Z/nZ (integers modulo n)
 *
 * @see Deviation: Modular Integer Coercion and Factories
 *
 * Port of: sage/rings/finite_rings/integer_mod.pyx
 */

import { znorder } from '@sagemath-ts/parigp-ts';
import {
  crt,
  factor,
  gcd,
  is_prime,
  is_prime_power,
  lcm,
  power_mod,
  primitive_root,
  xgcd,
} from '../../arith/misc.js';
import {
  ArithmeticError,
  AttributeError,
  IndexError,
  OverflowError,
  ValueError,
  ZeroDivisionError,
} from '../../errors.js';
import { discrete_log, has_order } from '../../groups/generic.js';
import type { IntegerLike } from '../../types/coercion.js';
import { Integer, type IntegerRing, ZZ } from '../integer_ring.js';
import type { RingElement } from '../polynomial/polynomial_element.js';
import { Rational } from '../rational.js';
import { QQ } from '../rational_field.js';
import { isFractionElement, _fraction_native_integer } from '../fraction_field_element.js';
import {
  FiniteFieldElement as ExtensionElement,
  PrimeField,
  PrimeFieldElement,
} from './finite_field_extension.js';
import {
  FiniteFieldPrime,
  FiniteFieldElement as LegacyPrimeElement,
} from './finite_field_prime.js';
import { IntegerModRing, Zmod } from './integer_mod_ring.js';

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
  /** IntegerMod._rational_: lift the canonical residue into QQ. */
  _rational_(): Rational {
    return new Rational(this.value);
  }

  readonly value: bigint;
  readonly parent: IntegerModRingBase;

  /**
   * Create an element of Z/nZ.
   *
   * @param value - The integer value (will be reduced modulo n)
   * @param parent - The parent ring Z/nZ
   */
  constructor(value: unknown, parent: IntegerModRingBase) {
    this.parent = parent;

    if (value instanceof ExtensionElement && value.lift.degree() > 0)
      throw new TypeError(`unable to convert ${value} to a rational`);
    const polynomial = value as {
      coeffs?: unknown[];
      getCoeff?: (i: number) => unknown;
      degree?: () => number;
      parent?: { base_ring: { zero(): unknown; toString(): string } };
    } | null;
    const isPolynomial = !!polynomial && Array.isArray(polynomial.coeffs) && !!polynomial.getCoeff;
    if (isPolynomial) {
      const base = polynomial!.parent!.base_ring;
      const zero = base.zero() as { parent?: { characteristic?: bigint }; modulus?: bigint };
      if (zero.parent?.characteristic === parent.modulus || zero.modulus === parent.modulus) {
        if (polynomial!.degree!() > 0) throw new TypeError(`${value} is not a constant polynomial`);
        value = polynomial!.getCoeff!(0);
      }
    }
    if (isFractionElement(value) && value.parent.characteristic() === parent.modulus) {
      // Canonical FpT prime-field section precedes IntegerMod's fallback constructor.
      const native = _fraction_native_integer(value);
      if (native !== undefined) value = native;
    }
    if (isFractionElement(value) || isPolynomial) {
      try {
        value = ZZ.__call__(value as unknown as Parameters<typeof ZZ.__call__>[0]);
      } catch (e) {
        if (!(e instanceof TypeError || e instanceof ValueError)) throw e;
        value = QQ.__call__(value as Parameters<typeof QQ.__call__>[0]);
      }
    }
    if (value instanceof Rational) {
      // Rational.__mod__: reduce the denominator before attempting inversion.
      const denominator = mod(value.denominator, parent.modulus);
      const [g, inverse] = xgcd(denominator, parent.modulus);
      if (g !== 1n)
        throw new ZeroDivisionError(
          `inverse of Mod(${denominator}, ${parent.modulus}) does not exist`
        );
      this.value = mod(value.numerator * inverse, parent.modulus);
    } else {
      this.value = mod(ZZ.__call__(value as Parameters<typeof ZZ.__call__>[0]), parent.modulus);
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
  add(other: ExtensionElement): ExtensionElement;
  add(other: PrimeFieldElement): PrimeFieldElement;
  add(other: LegacyPrimeElement): LegacyPrimeElement;
  add(other: IntegerMod | IntegerLike | number | boolean): IntegerMod;
  add(other: FiniteArithmeticElement | IntegerLike | number | boolean): FiniteArithmeticElement;
  add(other: IntegerMod): IntegerMod;
  add(other: unknown): FiniteArithmeticElement {
    const [left, operand] = canonicalFiniteOperands(this, other, '+');
    if (left instanceof ExtensionElement) return left.add(operand as ExtensionElement);
    const right = (operand as IntegerMod | PrimeFieldElement | LegacyPrimeElement).value;
    if (left instanceof PrimeFieldElement) return left.add(right);
    if (left instanceof LegacyPrimeElement) return left.add(right);
    return new IntegerMod(left.value + right, left.parent);
  }

  /**
   * Subtract two elements.
   */
  sub(other: ExtensionElement): ExtensionElement;
  sub(other: PrimeFieldElement): PrimeFieldElement;
  sub(other: LegacyPrimeElement): LegacyPrimeElement;
  sub(other: IntegerMod | IntegerLike | number | boolean): IntegerMod;
  sub(other: FiniteArithmeticElement | IntegerLike | number | boolean): FiniteArithmeticElement;
  sub(other: IntegerMod): IntegerMod;
  sub(other: unknown): FiniteArithmeticElement {
    const [left, operand] = canonicalFiniteOperands(this, other, '-');
    if (left instanceof ExtensionElement) return left.sub(operand as ExtensionElement);
    const right = (operand as IntegerMod | PrimeFieldElement | LegacyPrimeElement).value;
    if (left instanceof PrimeFieldElement) return left.sub(right);
    if (left instanceof LegacyPrimeElement) return left.sub(right);
    return new IntegerMod(left.value - right, left.parent);
  }

  /**
   * Multiply two elements.
   */
  mul(other: string): string;
  mul<T>(other: readonly T[]): T[];
  mul(other: ExtensionElement): ExtensionElement;
  mul(other: PrimeFieldElement): PrimeFieldElement;
  mul(other: LegacyPrimeElement): LegacyPrimeElement;
  mul(other: IntegerMod | IntegerLike | number | boolean): IntegerMod;
  mul(other: FiniteArithmeticElement | IntegerLike | number | boolean): FiniteArithmeticElement;
  mul(other: IntegerMod): IntegerMod;
  mul(other: unknown): FiniteArithmeticElement | string | unknown[] {
    if (typeof other === 'string' || Array.isArray(other))
      return repeatFiniteSequence(this.value, other);
    const [left, operand] = canonicalFiniteOperands(this, other, '*');
    if (left instanceof ExtensionElement) return left.mul(operand as ExtensionElement);
    const right = (operand as IntegerMod | PrimeFieldElement | LegacyPrimeElement).value;
    if (left instanceof PrimeFieldElement) return left.mul(right);
    if (left instanceof LegacyPrimeElement) return left.mul(right);
    return new IntegerMod(left.value * right, left.parent);
  }

  /**
   * Divide two elements.
   * This requires the divisor to be invertible (gcd(divisor, n) = 1).
   *
   * @throws {ZeroDivisionError} If the divisor is not invertible
   */
  div(other: ExtensionElement): ExtensionElement;
  div(other: PrimeFieldElement): PrimeFieldElement;
  div(other: LegacyPrimeElement): LegacyPrimeElement;
  div(other: IntegerMod | IntegerLike | number | boolean): IntegerMod;
  div(other: FiniteArithmeticElement | IntegerLike | number | boolean): FiniteArithmeticElement;
  div(other: IntegerMod): IntegerMod;
  div(other: unknown): FiniteArithmeticElement {
    const [left, operand] = canonicalFiniteOperands(this, other, '/');
    if (left instanceof ExtensionElement) return left.div(operand as ExtensionElement);
    const right = (operand as IntegerMod | PrimeFieldElement | LegacyPrimeElement).value;
    if (left instanceof PrimeFieldElement) return left.div(right);
    if (left instanceof LegacyPrimeElement) return left.div(right);
    const [g, inverse] = xgcd(right, left.modulus);
    if (g !== 1n)
      throw new ZeroDivisionError(`inverse of Mod(${right}, ${left.modulus}) does not exist`);
    return new IntegerMod(left.value * inverse, left.parent);
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
    // See `div` above: no zero shortcut, matching `integer_mod.pyx:2375`.
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
  pow(n: IntegerLike | number | Rational | boolean | string | null): IntegerMod {
    const exp = ZZ.__call__(n);
    const result = new IntegerMod(
      power_mod(this.value, exp < 0n ? -exp : exp, this.modulus),
      this.parent
    );
    if (exp >= 0n) return result;
    // Native backends invert the powered residue; GMP uses mpz_pow_helper.
    // The exponent cutover is strict, even for a small modulus.
    const nativeExponent = typeof n === 'bigint' || typeof n === 'number' || n instanceof Integer;
    if ((!nativeExponent || this.modulus > 2147483647n || exp <= -100000n) && !result.isUnit()) {
      throw new ZeroDivisionError('Inverse does not exist.');
    }
    return result.inv();
  }

  /**
   * Check equality with another element.
   */
  eq(other: unknown): boolean {
    return finiteArithmeticEquals(this, other);
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
    // Modulo 1 the one element is 0, so compare against `1 mod n` rather than 1.
    return this.value === mod(1n, this.modulus);
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
      // `integer_mod.pyx:1897` raises ArithmeticError, not ValueError.
      throw new ArithmeticError(
        `multiplicative order of ${this.value} not defined since it is not a unit modulo ${this.modulus}`
      );
    }

    return znorder(this.value, this.modulus);
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
      // An explicit base conversion uses the constructor, not arithmetic coercion.
      base = new IntegerMod(b, this.parent);
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
      // `integer_mod_ring.py:891-892` raises a bare ArithmeticError.
      throw new ArithmeticError('');
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
function createParent(modulus: bigint): IntegerModRing {
  return new IntegerModRing(modulus);
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
/** A dynamic zero modulus can return the caller's original value unchanged. */
type ModResult<T, N, E = IntegerMod> = N extends 0 | 0n | false
  ? T
  : N extends bigint
    ? bigint extends N
      ? T | E
      : E
    : N extends number
      ? number extends N
        ? T | E
        : E
      : T | E;

type ModParent =
  | IntegerModRingBase
  | IntegerRing
  | PrimeField
  | FiniteFieldPrime
  | null
  | undefined;
type ModParentElement<P> = P extends PrimeField
  ? PrimeFieldElement
  : P extends FiniteFieldPrime
    ? LegacyPrimeElement
    : IntegerMod;

export function Mod<
  T,
  N extends IntegerLike | number | boolean | Rational,
  P extends ModParent = undefined,
>(value: T, modulus: N, parent?: P): ModResult<T, N, ModParentElement<P>> {
  // Sage checks equality with zero before looking at the optional parent.
  const zero =
    typeof modulus === 'number' || typeof modulus === 'bigint' || typeof modulus === 'boolean'
      ? !modulus
      : modulus instanceof Rational
        ? modulus.numerator === 0n
        : modulus instanceof Integer
          ? modulus.value === 0n
          : false;
  if (zero) return value as ModResult<T, N, ModParentElement<P>>;
  const ring = parent ?? Zmod(modulus);
  if (ring instanceof PrimeField || ring instanceof FiniteFieldPrime) {
    return ring.__call__(value) as ModResult<T, N, ModParentElement<P>>;
  }
  if (!('modulus' in ring)) {
    throw new AttributeError(
      "'sage.rings.integer_ring.IntegerRing_class' object has no attribute '_pyx_order'"
    );
  }
  return new IntegerMod(value, ring) as ModResult<T, N, ModParentElement<P>>;
}

/** @internal Finite-ring subset of Sage's canonical coercion model. */
export type FiniteArithmeticElement =
  | IntegerMod
  | PrimeFieldElement
  | LegacyPrimeElement
  | ExtensionElement;

/**
 * @internal Canonical maps from finite_field_prime_modn.py:_coerce_map_from_,
 * finite_field_base.pyx:_coerce_map_from_, and QuotientFunctor.merge.
 * Named extension fields have no implicit embeddings between different degrees.
 * @see Deviation: Finite Field Coercion and Backend Boundaries
 */
export function canonicalFiniteOperands(
  left: FiniteArithmeticElement,
  right: unknown,
  operation: string
): [FiniteArithmeticElement, FiniteArithmeticElement] {
  if (
    typeof right === 'bigint' ||
    typeof right === 'boolean' ||
    right instanceof Integer ||
    (typeof right === 'number' && Number.isInteger(right))
  ) {
    return [left, left.parent.__call__(ZZ.__call__(right))];
  }
  const parentError = (parent: unknown): never => {
    throw new TypeError(
      `unsupported operand parent(s) for ${operation}: '${left.parent}' and '${parent}'`
    );
  };
  if (right instanceof Rational) return parentError('Rational Field');
  if (
    !(
      right instanceof IntegerMod ||
      right instanceof PrimeFieldElement ||
      right instanceof LegacyPrimeElement ||
      right instanceof ExtensionElement
    )
  ) {
    const type =
      right == null
        ? 'NoneType'
        : typeof right === 'string'
          ? 'str'
          : Array.isArray(right)
            ? 'list'
            : typeof right === 'number'
              ? 'float'
              : 'object';
    const modulus = left instanceof IntegerMod ? left.modulus : left.parent.characteristic;
    const backend = modulus <= 46341n ? 'int' : modulus <= 2147483647n ? 'int64' : 'gmp';
    const elementType =
      left instanceof ExtensionElement
        ? 'sage.rings.finite_rings.element_pari_ffelt.FiniteFieldElement_pari_ffelt'
        : `sage.rings.finite_rings.integer_mod.IntegerMod_${backend}`;
    if (operation === '*' && (type === 'str' || type === 'list'))
      throw new TypeError(`can't multiply sequence by non-int of type '${elementType}'`);
    throw new TypeError(
      `unsupported operand type(s) for ${operation}: '${elementType}' and '${type}'`
    );
  }
  if (left.parent === right.parent) return [left, right];
  if (left instanceof ExtensionElement || right instanceof ExtensionElement) {
    if (left instanceof ExtensionElement && right instanceof ExtensionElement) {
      const a = left.parent;
      const b = right.parent;
      if (
        a.characteristic !== b.characteristic ||
        a.degree !== b.degree ||
        a.variableName !== b.variableName ||
        !a.modulus.eq(b.modulus)
      )
        return parentError(b);
      return [left, a.__call__(right)];
    }
    const extension = left instanceof ExtensionElement ? left : (right as ExtensionElement);
    const scalar = (left instanceof ExtensionElement ? right : left) as
      | IntegerMod
      | PrimeFieldElement
      | LegacyPrimeElement;
    const modulus = scalar instanceof IntegerMod ? scalar.modulus : scalar.parent.characteristic;
    if (modulus % extension.parent.characteristic !== 0n) return parentError(right.parent);
    const converted = extension.parent.__call__(scalar.value);
    return left instanceof ExtensionElement ? [left, converted] : [converted, extension];
  }
  if (left instanceof IntegerMod && right instanceof IntegerMod) {
    if (right.modulus % left.modulus === 0n) return [left, left.parent.__call__(right.value)];
    if (left.modulus % right.modulus === 0n) return [right.parent.__call__(left.value), right];
    const common = gcd(left.modulus, right.modulus);
    if (common === 1n) return parentError(right.parent);
    const parent = new IntegerModRing(common);
    return [parent.__call__(left.value), parent.__call__(right.value)];
  }
  const prime =
    left instanceof IntegerMod ? (right as PrimeFieldElement | LegacyPrimeElement) : left;
  const other = left instanceof IntegerMod ? left : right;
  const modulus = other instanceof IntegerMod ? other.modulus : other.parent.characteristic;
  if (modulus % prime.parent.characteristic !== 0n) return parentError(right.parent);
  return [prime.parent.__call__(left.value), prime.parent.__call__(right.value)];
}

/** @internal Sage richcmp returns false when no canonical common parent exists. */
export function finiteArithmeticEquals(left: FiniteArithmeticElement, right: unknown): boolean {
  try {
    const [a, b] = canonicalFiniteOperands(left, right, '==');
    if (a instanceof ExtensionElement) return a.lift.eq((b as ExtensionElement).lift);
    return a.value === (b as IntegerMod | PrimeFieldElement | LegacyPrimeElement).value;
  } catch (error) {
    if (error instanceof TypeError) return false;
    throw error;
  }
}

/** @internal Python's reflected sequence multiplication uses IntegerMod.__index__. */
export function repeatFiniteSequence<T>(
  count: bigint,
  sequence: string | readonly T[]
): string | T[] {
  if (count > 9223372036854775807n)
    throw new OverflowError(
      "cannot fit 'sage.rings.finite_rings.integer_mod.IntegerMod_gmp' into an index-sized integer"
    );
  if (sequence.length === 0) return typeof sequence === 'string' ? '' : [];
  if (typeof sequence === 'string') return sequence.repeat(Number(count));
  const result: T[] = [];
  for (let i = 0n; i < count; i++) for (const value of sequence) result.push(value);
  return result;
}

/**
 * @internal Validate a single generator index. Fields use Python truthiness
 * (finite_field_prime_modn.py:gen, finite_field_pari_ffelt.py:gen); quotient
 * rings use ZZ.gen's equality-to-zero test through QuotientRing.gen.
 */
export function checkFiniteGeneratorIndex(index: unknown, field: boolean): void {
  const scalarZero =
    index === 0 ||
    index === 0n ||
    index === false ||
    (index instanceof Integer && index.value === 0n) ||
    (index instanceof Rational && index.numerator === 0n) ||
    ((index instanceof IntegerMod ||
      index instanceof PrimeFieldElement ||
      index instanceof LegacyPrimeElement ||
      index instanceof ExtensionElement) &&
      index.isZero());
  if (
    scalarZero ||
    (field && (index == null || index === '' || (Array.isArray(index) && index.length === 0)))
  )
    return;
  throw new IndexError(field ? 'only one generator' : 'n must be 0');
}
