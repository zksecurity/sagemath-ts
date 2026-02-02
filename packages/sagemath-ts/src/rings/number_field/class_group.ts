/**
 * @module sage/rings/number_field/class_group
 * @description Class groups of number fields
 *
 * Port of: sage/rings/number_field/class_group.py
 * Reference: reference/sage/src/sage/rings/number_field/class_group.py
 *
 * NOTE: Full class group computation requires PARI's bnfinit.
 * This implementation provides the structure and basic operations.
 * For full PARI compatibility, see DEVIATIONS.md.
 */

import { NotImplementedError, ValueError } from '../../errors.js';
import type { NumberField } from './number_field.js';
import type { NumberFieldIdeal } from './number_field_ideal.js';

/**
 * The class group of a number field.
 *
 * The class group is the quotient of the group of fractional ideals by
 * the subgroup of principal ideals.
 *
 * Class group computation is one of the fundamental problems in
 * algebraic number theory. The structure is:
 *   Cl(K) ≅ Z/d_1 × Z/d_2 × ... × Z/d_r
 * where d_i | d_{i+1} (elementary divisors).
 *
 * The class number h(K) = |Cl(K)| = d_1 * d_2 * ... * d_r.
 *
 * @see Reference: sage/rings/number_field/class_group.py:ClassGroup
 */
export class ClassGroup {
  private readonly _number_field: NumberField;
  private readonly _invariants: readonly bigint[];
  private readonly _gens: readonly ClassGroupElement[];

  constructor(
    number_field: NumberField,
    invariants: bigint[] = [],
    gens: ClassGroupElement[] = []
  ) {
    this._number_field = number_field;
    this._invariants = Object.freeze([...invariants]);
    this._gens = Object.freeze([...gens]);
  }

  /**
   * Return the number field this class group belongs to.
   */
  number_field(): NumberField {
    return this._number_field;
  }

  /**
   * Return the class number (order of the class group).
   *
   * The class number h(K) = |Cl(K)| is a fundamental invariant of the number field.
   * h(K) = 1 means the ring of integers is a unique factorization domain.
   *
   * @see Reference: sage/rings/number_field/class_group.py:order
   */
  order(): bigint {
    if (this._invariants.length === 0) {
      return 1n;
    }
    return this._invariants.reduce((a, b) => a * b, 1n);
  }

  /**
   * Return the class number.
   *
   * Alias for order().
   *
   * @see Reference: sage/rings/number_field/number_field.py:class_number
   */
  class_number(): bigint {
    return this.order();
  }

  /**
   * Return the exponent of the class group.
   *
   * The exponent is the smallest positive integer e such that
   * x^e = 1 for all x in the class group. It equals the largest invariant.
   */
  exponent(): bigint {
    if (this._invariants.length === 0) {
      return 1n;
    }
    return this._invariants[this._invariants.length - 1]!;
  }

  /**
   * Return the invariants (elementary divisors) of the class group.
   *
   * These are the d_i such that Cl(K) ≅ Z/d_1 × Z/d_2 × ... × Z/d_r
   * with d_i | d_{i+1}.
   */
  invariants(): bigint[] {
    return [...this._invariants];
  }

  /**
   * Return the generators of the class group.
   *
   * These are ideal classes that generate the class group.
   * Each generator has order equal to the corresponding invariant.
   *
   * @see Reference: sage/rings/number_field/class_group.py:gens
   */
  gens(): ClassGroupElement[] {
    return [...this._gens];
  }

  /**
   * Return the generators of the class group.
   *
   * Alias for gens(). Returns ideal class generators that generate the class group.
   *
   * @see Reference: sage/rings/number_field/class_group.py:gens
   */
  class_group_generators(): ClassGroupElement[] {
    return this.gens();
  }

  /**
   * Return the number of generators.
   */
  ngens(): number {
    return this._gens.length;
  }

  /**
   * Return the i-th generator.
   */
  gen(i: number): ClassGroupElement {
    if (i < 0 || i >= this._gens.length) {
      throw new ValueError(`generator index ${i} out of range`);
    }
    return this._gens[i]!;
  }

  /**
   * Return the identity element.
   */
  identity(): ClassGroupElement {
    // The identity is represented by zero exponents for each invariant
    return new ClassGroupElement(
      this,
      Array(this._invariants.length)
        .fill(null)
        .map(() => 0n)
    );
  }

  /**
   * Return the class of an ideal.
   *
   * Maps an ideal to its equivalence class in the class group.
   * Two ideals are equivalent if their quotient is a principal ideal.
   *
   * @see Reference: sage/rings/number_field/class_group.py:__call__
   */
  __call__(ideal: NumberFieldIdeal): ClassGroupElement {
    // For trivial class group (class number 1), everything is principal
    if (this.order() === 1n) {
      return this.identity();
    }

    // For principal ideals (1 generator), return identity
    if (ideal.ngens() === 1) {
      return this.identity();
    }

    // For quadratic fields with known class groups, we can compute the class
    const K = this._number_field;
    if (K.degree() === 2) {
      const norm = ideal.norm();
      const invariants = this.invariants();

      // For cyclic class groups, find the exponent
      if (invariants.length === 1) {
        const h = invariants[0]!;

        // The class of an ideal is determined by the prime factorization
        // For small class numbers, we can search
        if (h <= 10n) {
          // Check each class
          for (let e = 0n; e < h; e++) {
            const classElem = new ClassGroupElement(this, [e]);
            // Would need to compare representative ideals
            // For now, return based on norm congruence
          }
        }
      }
    }

    throw new NotImplementedError('__call__ requires PARI bnfisprincipal for general ideals');
  }

  /**
   * Return a random element of the class group.
   *
   * Generates random exponents for each generator.
   */
  random_element(): ClassGroupElement {
    const invariants = this._invariants;
    const exponents: bigint[] = [];

    for (const d of invariants) {
      // Random exponent in [0, d)
      const randomExp = BigInt(Math.floor(Math.random() * Number(d)));
      exponents.push(randomExp);
    }

    return new ClassGroupElement(this, exponents);
  }

  /**
   * Check if the class group is trivial (class number is 1).
   */
  is_trivial(): boolean {
    return this.order() === 1n;
  }

  /**
   * Check if the class group is cyclic.
   */
  is_cyclic(): boolean {
    return this._invariants.length <= 1;
  }

  /**
   * Return a list of all elements.
   *
   * WARNING: This can be very expensive for large class groups.
   */
  list(): ClassGroupElement[] {
    // Generate all combinations of exponents
    const elements: ClassGroupElement[] = [];
    const numInvariants = this._invariants.length;

    if (numInvariants === 0) {
      return [this.identity()];
    }

    const generateElements = (idx: number, exponents: bigint[]): void => {
      if (idx === numInvariants) {
        elements.push(new ClassGroupElement(this, [...exponents]));
        return;
      }

      const order = this._invariants[idx]!;
      for (let e = 0n; e < order; e++) {
        exponents[idx] = e;
        generateElements(idx + 1, exponents);
      }
    };

    generateElements(0, Array(numInvariants).fill(0n));
    return elements;
  }

  /**
   * Iterate over all elements.
   */
  *[Symbol.iterator](): Iterator<ClassGroupElement> {
    for (const element of this.list()) {
      yield element;
    }
  }

  toString(): string {
    const order = this.order();
    if (order === 1n) {
      return `Trivial class group of ${this._number_field}`;
    }

    const structure = this._invariants.map((d) => `C${d}`).join(' x ');
    return `Class group of order ${order} with structure ${structure} of ${this._number_field}`;
  }
}

/**
 * An element of a class group.
 *
 * Represented by exponents with respect to the generators.
 *
 * @see Reference: sage/rings/number_field/class_group.py:FractionalIdealClass
 */
export class ClassGroupElement {
  private readonly _parent: ClassGroup;
  private readonly _exponents: readonly bigint[];
  private _ideal?: NumberFieldIdeal;

  constructor(parent: ClassGroup, exponents: bigint[], ideal?: NumberFieldIdeal) {
    this._parent = parent;
    // Reduce exponents modulo the invariants
    const invariants = parent.invariants();
    const normalizedExps = exponents.map((e, i) => {
      const d = invariants[i];
      if (d === undefined || d === 0n) return e;
      return ((e % d) + d) % d;
    });
    this._exponents = Object.freeze(normalizedExps);
    this._ideal = ideal;
  }

  /**
   * Return the parent class group.
   */
  parent(): ClassGroup {
    return this._parent;
  }

  /**
   * Return a representative ideal for this class.
   */
  ideal(): NumberFieldIdeal {
    if (!this._ideal) {
      throw new NotImplementedError('ideal computation requires PARI');
    }
    return this._ideal;
  }

  /**
   * Return the exponent vector with respect to generators.
   */
  exponents(): bigint[] {
    return [...this._exponents];
  }

  /**
   * Return the order of this element in the class group.
   *
   * The order is the smallest positive integer n such that this^n = identity.
   */
  order(): bigint {
    // The order divides the exponent of the group
    const groupOrder = this._parent.order();
    const invariants = this._parent.invariants();

    if (this.is_trivial()) {
      return 1n;
    }

    // Compute order as LCM of orders in each cyclic factor
    let order = 1n;
    for (let i = 0; i < this._exponents.length; i++) {
      const e = this._exponents[i]!;
      const d = invariants[i]!;
      if (e !== 0n) {
        const g = gcdBigInt(e, d);
        const localOrder = d / g;
        order = lcmBigInt(order, localOrder);
      }
    }

    return order;
  }

  /**
   * Check if this is the trivial class.
   */
  is_trivial(): boolean {
    return this._exponents.every((e) => e === 0n);
  }

  /**
   * Check if this is the trivial class (alias for is_trivial).
   */
  is_principal(): boolean {
    return this.is_trivial();
  }

  /**
   * Multiply two class group elements.
   */
  mul(other: ClassGroupElement): ClassGroupElement {
    if (this._parent !== other._parent) {
      throw new ValueError('elements must be in the same class group');
    }
    const newExponents = this._exponents.map((e, i) => e + other._exponents[i]!);
    return new ClassGroupElement(this._parent, newExponents);
  }

  /**
   * Return the inverse.
   */
  inv(): ClassGroupElement {
    const newExponents = this._exponents.map((e) => -e);
    return new ClassGroupElement(this._parent, newExponents);
  }

  /**
   * Return self raised to power n.
   */
  pow(n: bigint): ClassGroupElement {
    const newExponents = this._exponents.map((e) => e * n);
    return new ClassGroupElement(this._parent, newExponents);
  }

  /**
   * Check equality.
   */
  eq(other: ClassGroupElement): boolean {
    if (this._parent !== other._parent) {
      return false;
    }
    for (let i = 0; i < this._exponents.length; i++) {
      if (this._exponents[i] !== other._exponents[i]) {
        return false;
      }
    }
    return true;
  }

  toString(): string {
    if (this.is_trivial()) {
      return 'Trivial principal fractional ideal class';
    }
    const exps = this._exponents.map((e, i) => `g${i}^${e}`).join(' * ');
    return `Fractional ideal class (${exps})`;
  }
}

/**
 * The narrow (or strict) class group of a number field.
 *
 * The narrow class group uses totally positive principal ideals
 * instead of all principal ideals.
 *
 * @see Reference: sage/rings/number_field/class_group.py:NarrowClassGroup
 */
export class NarrowClassGroup extends ClassGroup {
  // The narrow class group considers only totally positive elements
  // as generators of principal ideals
}

/**
 * The ray class group modulo a modulus m.
 *
 * @see Reference: sage/rings/number_field/class_group.py:RayClassGroup
 */
export class RayClassGroup extends ClassGroup {
  private readonly _modulus: unknown;

  constructor(
    number_field: NumberField,
    modulus: unknown,
    invariants: bigint[] = [],
    gens: ClassGroupElement[] = []
  ) {
    super(number_field, invariants, gens);
    this._modulus = modulus;
  }

  /**
   * Return the modulus.
   */
  modulus(): unknown {
    return this._modulus;
  }
}

// Helper functions

function gcdBigInt(a: bigint, b: bigint): bigint {
  a = a < 0n ? -a : a;
  b = b < 0n ? -b : b;
  while (b !== 0n) {
    [a, b] = [b, a % b];
  }
  return a;
}

function lcmBigInt(a: bigint, b: bigint): bigint {
  if (a === 0n || b === 0n) return 0n;
  return (a * b) / gcdBigInt(a, b);
}
