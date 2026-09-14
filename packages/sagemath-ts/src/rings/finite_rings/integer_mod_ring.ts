/**
 * @module sage/rings/finite_rings/integer_mod_ring
 * @description The ring Z/nZ of integers modulo n
 *
 * @see Deviation: Modular Integer Coercion and Factories
 *
 * Port of: sage/rings/finite_rings/integer_mod_ring.py
 */

import { gcd, is_prime, factor, CRT_basis, type Factorization } from '../../arith/misc.js';
import { ZeroDivisionError, ValueError, NotImplementedError } from '../../errors.js';
import { current_randstate } from '../../misc/randstate.js';
import { toBigInt, type IntegerLike } from '../../types/coercion.js';
import { type IntegerRing, ZZ } from '../integer_ring.js';
import {
  type Polynomial,
  type CoefficientRing,
  type RingElement,
} from '../polynomial/polynomial_element.js';
import { PolynomialRing } from '../polynomial/polynomial_ring.js';
import { GF } from './finite_field_constructor.js';
import type { FiniteFieldPrime, FiniteFieldElement } from './finite_field_prime.js';
import type { Rational } from '../rational.js';
import {
  IntegerMod,
  type IntegerModRingBase,
  checkFiniteGeneratorIndex,
  multiplicative_generator,
  multiplicative_group_is_cyclic,
  unit_gens,
} from './integer_mod.js';

/**
 * The ring Z/nZ of integers modulo n.
 *
 * This class implements the ring of integers modulo n, where n is a positive integer.
 * When n is prime, this ring is a field (but for fields, use GF(p) for additional
 * field-specific functionality).
 *
 * @example
 * ```typescript
 * const Z7 = Zmod(7n);
 * const a = Z7(3n);
 * const b = Z7(4n);
 * console.log(a.mul(b));  // 5 (since 3*4 = 12 ≡ 5 mod 7)
 *
 * // Iterate over elements
 * for (const x of Z7) {
 *   console.log(x.value);  // 0, 1, 2, 3, 4, 5, 6
 * }
 * ```
 */
export class IntegerModRing implements IntegerModRingBase, CoefficientRing<IntegerMod> {
  readonly modulus: bigint;
  readonly characteristic: bigint;
  readonly order: bigint;

  private _zero: IntegerMod | null = null;
  private _one: IntegerMod | null = null;

  /**
   * Create the ring Z/nZ.
   *
   * @param n - The modulus (must be a positive integer)
   */
  constructor(n: IntegerLike | number) {
    const modulus = ZZ.__call__(n);

    if (modulus <= 0n) {
      throw new ZeroDivisionError('order must be positive');
    }

    this.modulus = modulus;
    this.characteristic = modulus;
    this.order = modulus;
  }

  /**
   * Create an element of this ring.
   *
   * @param x - The value to convert to an element
   */
  __call__(x?: unknown): IntegerMod {
    return new IntegerMod(x, this);
  }

  /**
   * Return the zero element.
   */
  zero(): IntegerMod {
    if (this._zero === null) {
      this._zero = new IntegerMod(0n, this);
    }
    return this._zero;
  }

  /**
   * Return the one element.
   */
  one(): IntegerMod {
    if (this._one === null) {
      this._one = new IntegerMod(1n, this);
    }
    return this._one;
  }

  /**
   * Return a generator of the ring (which is just 1).
   */
  gen(n: unknown = 0): IntegerMod {
    checkFiniteGeneratorIndex(n, false);
    return this.one();
  }

  /**
   * Check if this ring is a field.
   * Z/nZ is a field if and only if n is prime.
   */
  is_field(): boolean {
    return is_prime(this.modulus);
  }

  /**
   * Corresponding finite field, cached on this ring.
   * @see Reference: sage/rings/finite_rings/integer_mod_ring.py:field
   * @see Deviation: Modular Polynomial Roots and Hensel Lifting
   */
  field(): FiniteFieldPrime {
    const cached = correspondingFields.get(this);
    if (cached) return cached;
    if (!this.is_field()) throw new ValueError('self must be a field');
    const field = GF(this.order);
    correspondingFields.set(this, field);
    return field;
  }

  /**
   * Cached factorization of the ring order, using the PARI-backed factor port.
   * @see Reference: sage/rings/finite_rings/integer_mod_ring.py:factored_order
   * @see Deviation: Modular Polynomial Roots and Hensel Lifting
   */
  factored_order(): Factorization {
    const cached = factoredOrders.get(this);
    if (cached) return cached;
    const result = factor(this.order);
    factoredOrders.set(this, result);
    return result;
  }

  /**
   * Lift a residue-field root from p to p^e, where p is prime and e >= 1.
   * Nonzero derivative doubles precision; singular roots lift one digit at a time.
   * The input root must lie in Zmod(p) and vanish under f modulo p.
   * @see Reference: sage/rings/finite_rings/integer_mod_ring.py:_lift_residue_field_root
   * @see Deviation: Modular Polynomial Roots and Hensel Lifting
   */
  static _lift_residue_field_root(
    p: IntegerLike,
    e: IntegerLike,
    f: Polynomial<IntegerMod & RingElement>,
    fprime: Polynomial<IntegerMod & RingElement>,
    root: IntegerMod
  ): IntegerMod[] {
    const prime = toBigInt(p),
      exponent = toBigInt(e);
    if (exponent === 1n) return [root];
    const deriv = fprime.evaluate(root as IntegerMod & RingElement);
    if (!deriv.isZero()) {
      let precision = 1n;
      while (true) {
        precision = 2n * precision < exponent ? 2n * precision : exponent;
        const K = Zmod(prime ** precision) as IntegerModRing;
        root = K.__call__(root.value);
        const step = f
          .evaluate(root as IntegerMod & RingElement)
          .div(fprime.evaluate(root as IntegerMod & RingElement));
        root = root.sub(step);
        if (precision >= exponent) return [root];
      }
    }
    let modulus = prime,
      roots = [root];
    for (let precision = 1n; precision < exponent; precision++) {
      const increment = modulus;
      modulus *= prime;
      const K = Zmod(modulus) as IntegerModRing,
        next: IntegerMod[] = [];
      for (const previous of roots) {
        let candidate = K.__call__(previous.value);
        if (!f.evaluate(candidate as IntegerMod & RingElement).isZero()) continue;
        next.push(candidate);
        for (let digit = 1n; digit < prime; digit++) {
          candidate = candidate.add(K.__call__(increment));
          next.push(candidate);
        }
      }
      roots = next;
    }
    return roots;
  }

  /**
   * Base-ring roots using finite-field factorization or CRT and Hensel lifting.
   * Default roots belong to field(); distinct roots belong to this ring.
   * Retains the bundled nonunit-linear recursion, including its upstream bug.
   * @see Reference: sage/rings/finite_rings/integer_mod_ring.py:_roots_univariate_polynomial
   * @see Deviation: Modular Polynomial Roots and Hensel Lifting
   */
  _roots_univariate_polynomial(
    f: Polynomial<IntegerMod & RingElement>,
    options: { ring?: IntegerModRing | null; multiplicities?: boolean; algorithm?: unknown } = {}
  ): Array<[FiniteFieldElement, number]> | IntegerMod[] {
    if (options.ring != null && options.ring !== this) throw new NotImplementedError('');
    const multiplicities = options.multiplicities ?? true;
    const degree = f.degree();
    if (multiplicities) {
      if (degree < 0 || !this.is_field())
        throw new NotImplementedError(
          'root finding with multiplicities for this polynomial not implemented (try the multiplicities=False option)'
        );
      return new PolynomialRing(this.field(), f.parent.variable_name)
        .__call__(f.coeffs.map((c) => c.value))
        .roots();
    }
    if (degree < 0) return [...this];
    if (degree === 0) return [];
    if (degree === 1) {
      const b = f.getCoeff(0),
        a = f.getCoeff(1);
      if (a.isUnit()) return [b.neg().mul(a.inv())];
      const g = gcd(this.order, a.value);
      if (b.value % g !== 0n) return [];
      const quotient = this.order / g,
        K = Zmod(quotient) as IntegerModRing;
      // Preserve the bundled source's recursion on f itself. It computes a/g
      // but never divides f by g, so some returned values are not roots.
      // See the explicit upstream-behavior note in DEVIATIONS.md.
      const reduced = new PolynomialRing(
        K as unknown as CoefficientRing<IntegerMod & RingElement>,
        f.parent.variable_name
      ).__call__(f.coeffs.map((c) => c.value));
      const first = this.__call__(reduced.roots({ multiplicities: false })[0]!.value);
      const increment = this.__call__(quotient),
        result: IntegerMod[] = [];
      for (let k = 0n; k < g; k++) result.push(first.add(increment.mul(k)));
      return result;
    }
    if (this.is_field()) {
      return new PolynomialRing(this.field(), f.parent.variable_name)
        .__call__(f.coeffs.map((c) => c.value))
        .roots({ multiplicities: false })
        .map((root) => this.__call__(root.value));
    }
    const factors = this.factored_order(),
      primePowerRoots: IntegerMod[][] = [];
    for (const [p, e] of factors) {
      const K = Zmod(p ** e) as IntegerModRing,
        Fp = Zmod(p) as IntegerModRing;
      const lifted = new PolynomialRing(
        K as unknown as CoefficientRing<IntegerMod & RingElement>,
        f.parent.variable_name
      ).__call__(f.coeffs.map((c) => c.value));
      const derivative = lifted.derivative();
      const reduced = new PolynomialRing(
        Fp as unknown as CoefficientRing<IntegerMod & RingElement>,
        f.parent.variable_name
      ).__call__(f.coeffs.map((c) => c.value));
      const values: IntegerMod[] = [];
      for (const root of reduced.roots({ multiplicities: false })) {
        values.push(...IntegerModRing._lift_residue_field_root(p, e, lifted, derivative, root));
      }
      primePowerRoots.push(values);
    }
    const basis = CRT_basis(factors.map(([p, e]) => p ** e)) as bigint[];
    let result = [this.zero()];
    // Cartesian-product ordering: the last prime-power component varies fastest.
    for (let i = 0; i < primePowerRoots.length; i++) {
      const next: IntegerMod[] = [];
      for (const previous of result)
        for (const root of primePowerRoots[i]!)
          next.push(previous.add(this.__call__(basis[i]! * root.value)));
      result = next;
    }
    return result;
  }

  /**
   * Return the cardinality of this ring.
   */
  cardinality(): bigint {
    return this.modulus;
  }

  /**
   * Iterate over all elements of this ring.
   *
   * @example
   * ```typescript
   * const Z5 = Zmod(5n);
   * for (const x of Z5) {
   *   console.log(x.value);  // 0, 1, 2, 3, 4
   * }
   * ```
   */
  *[Symbol.iterator](): Iterator<IntegerMod> {
    for (let i = 0n; i < this.modulus; i++) {
      yield new IntegerMod(i, this);
    }
  }

  /**
   * Return a list of all elements.
   */
  list(): IntegerMod[] {
    return Array.from(this);
  }

  /**
   * Return a list of all units (invertible elements).
   */
  units(): IntegerMod[] {
    const result: IntegerMod[] = [];
    // Start at 0: modulo 1 the element 0 is a unit (gcd(0, 1) = 1), and Sage's
    // `Integers(1).list_of_elements_of_multiplicative_group()` returns `[0]`.
    for (let i = 0n; i < this.modulus; i++) {
      if (gcd(i, this.modulus) === 1n) {
        result.push(new IntegerMod(i, this));
      }
    }
    return result;
  }

  /**
   * Return ``true`` if the multiplicative group of this ring is cyclic.
   *
   * @see Reference: sage/rings/finite_rings/integer_mod_ring.py:810
   */
  multiplicative_group_is_cyclic(): boolean {
    return multiplicative_group_is_cyclic(this.modulus);
  }

  /**
   * Return generators for the unit group `(Z/nZ)*`.
   *
   * `sage: Integers(75).unit_gens()` -> `(26, 52)`.
   *
   * @see Reference: sage/rings/finite_rings/integer_mod_ring.py:1442
   */
  unit_gens(): IntegerMod[] {
    return unit_gens(this.modulus).map(([g]) => new IntegerMod(g, this));
  }

  /**
   * Return a generator for the multiplicative group of this ring, assuming
   * the multiplicative group is cyclic.
   *
   * @throws {ValueError} `multiplicative group of this ring is not cyclic`
   * @see Reference: sage/rings/finite_rings/integer_mod_ring.py:849
   */
  multiplicative_generator(): IntegerMod {
    return new IntegerMod(multiplicative_generator(this.modulus), this);
  }

  /**
   * Return a random element of this ring.
   */
  random_element(bound?: IntegerLike | number | null): IntegerMod {
    const random = current_randstate().python_random();
    const limit = bound == null ? undefined : ZZ.__call__(bound);
    return this.__call__(
      limit === undefined ? random.randint(0n, this.modulus - 1n) : random.randint(-limit, limit)
    );
  }

  /**
   * String representation.
   */
  toString(): string {
    return `Ring of integers modulo ${this.modulus}`;
  }
}

/**
 * Factory function to create Z/nZ.
 *
 * @param n - The modulus (must be a positive integer)
 * @returns The ring Z/nZ
 *
 * @example
 * ```typescript
 * const Z12 = Zmod(12n);
 * const a = Z12(7n);
 * const b = Z12(8n);
 * console.log(a.add(b));  // 3 (since 7+8 = 15 ≡ 3 mod 12)
 * ```
 */
/** Zero literals narrow to ZZ; a dynamic order may produce either parent. */
type ModularRingFor<N> = N extends 0 | 0n | false | undefined
  ? IntegerRing
  : N extends bigint
    ? bigint extends N
      ? IntegerModRing | IntegerRing
      : IntegerModRing
    : N extends number
      ? number extends N
        ? IntegerModRing | IntegerRing
        : IntegerModRing
      : N extends true
        ? IntegerModRing
        : IntegerModRing | IntegerRing;

// UniqueFactory uses a WeakValueDictionary: unused parent objects may be collected.
const modularRingCache = new Map<bigint, WeakRef<IntegerModRing>>();
const modularRingFinalizer = new FinalizationRegistry<{
  order: bigint;
  reference: WeakRef<IntegerModRing>;
}>(({ order, reference }) => {
  if (modularRingCache.get(order) === reference) modularRingCache.delete(order);
});

/** IntegerModFactory.create_object, integer_mod_ring.py:231-245. */
export function Zmod<
  N extends IntegerLike | number | boolean | Rational | null | undefined = undefined,
>(n?: N): ModularRingFor<N> {
  // The factory compares the order before applying the generic ring's ZZ coercion.
  if (n === null || typeof n === 'string') {
    throw new TypeError(
      `'<' not supported between instances of '${n === null ? 'NoneType' : 'str'}' and 'int'`
    );
  }
  const order = ZZ.__call__(n);
  if (order === 0n) return ZZ as ModularRingFor<N>;
  // Sage's factory key is the supplied signed order; normalization is later.
  let ring = modularRingCache.get(order)?.deref();
  if (ring === undefined) {
    ring = new IntegerModRing(order < 0n ? -order : order);
    const reference = new WeakRef(ring);
    modularRingCache.set(order, reference);
    modularRingFinalizer.register(ring, { order, reference });
  }
  return ring as ModularRingFor<N>;
}

/**
 * Alias for Zmod - IntegerModRing factory.
 */
export const IntegerModRingFactory = Zmod;

/**
 * Alias for Zmod - Integers function.
 */
export const Integers = Zmod;

const correspondingFields = new WeakMap<IntegerModRing, FiniteFieldPrime>();
const factoredOrders = new WeakMap<IntegerModRing, Factorization>();
