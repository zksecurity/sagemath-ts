/**
 * @module sage/rings/polynomial/multi_polynomial_ring
 * @description Multivariate polynomial rings over arbitrary coefficient rings
 *
 * Port of: sage/rings/polynomial/multi_polynomial_ring.py
 *
 * Provides multivariate polynomial ring construction suitable for
 * ZK constraint systems where polynomials in multiple variables are common.
 *
 * @example
 * ```typescript
 * // Create a polynomial ring in three variables over a prime field
 * const F = GF(101n);
 * const R = new MPolynomialRing(F, ['x', 'y', 'z']);
 * const [x, y, z] = R.gens();
 *
 * // Build a constraint: x^2 + y*z - 1 = 0
 * const constraint = x.pow(2).add(y.mul(z)).sub(R.one());
 *
 * // Evaluate at a witness
 * const result = constraint.evaluate({ x: F(3n), y: F(5n), z: F(2n) });
 * ```
 */

import { ValueError } from '../../errors.js';
import type { CoefficientRing, RingElement } from './polynomial_element.js';
import {
  MPolynomial,
  type MPolynomialRingBase,
  type TermOrder,
  exponentToKey,
  keyToExponent,
} from './multi_polynomial_element.js';

/**
 * A multivariate polynomial ring R[x1, x2, ..., xn] over a base ring R.
 *
 * Supports customizable term orderings for monomial comparisons,
 * which is important for leading term computations and Groebner basis algorithms.
 *
 * @example
 * ```typescript
 * // Create GF(p)[x, y, z] with degrevlex order (default)
 * const R = new MPolynomialRing(GF(p), ['x', 'y', 'z']);
 *
 * // Create with lexicographic order
 * const S = new MPolynomialRing(GF(p), ['x', 'y'], 'lex');
 *
 * // Get generators
 * const [x, y, z] = R.gens();
 *
 * // Build polynomials
 * const f = x.mul(y).add(z);  // x*y + z
 * ```
 */
export class MPolynomialRing<C extends RingElement> implements MPolynomialRingBase<C> {
  /**
   * The coefficient ring.
   */
  readonly base_ring: CoefficientRing<C>;

  /**
   * The variable names.
   */
  readonly names: readonly string[];

  /**
   * The number of generators (variables).
   */
  readonly ngens_value: number;

  /**
   * The term ordering used for monomial comparison.
   */
  readonly term_order: TermOrder;

  /**
   * Cached generators.
   */
  private _gens: MPolynomial<C>[] | null = null;

  /**
   * Cached zero polynomial.
   */
  private _zero: MPolynomial<C> | null = null;

  /**
   * Cached one polynomial.
   */
  private _one: MPolynomial<C> | null = null;

  /**
   * Create a multivariate polynomial ring.
   *
   * @param base_ring - The coefficient ring
   * @param names - Array of variable names, e.g., ['x', 'y', 'z']
   * @param order - Term ordering: 'lex', 'deglex', or 'degrevlex' (default)
   *
   * @example
   * ```typescript
   * const R = new MPolynomialRing(GF(p), ['x', 'y', 'z']);
   * const S = new MPolynomialRing(QQ, ['a', 'b'], 'lex');
   * ```
   */
  constructor(base_ring: CoefficientRing<C>, names: string[], order: TermOrder = 'degrevlex') {
    if (names.length === 0) {
      throw new ValueError('must have at least one variable');
    }

    // Check for duplicate names
    const nameSet = new Set(names);
    if (nameSet.size !== names.length) {
      throw new ValueError('variable names must be distinct');
    }

    this.base_ring = base_ring;
    this.names = Object.freeze([...names]);
    this.ngens_value = names.length;
    this.term_order = order;
  }

  /**
   * Return the number of generators (variables).
   */
  ngens(): number {
    return this.ngens_value;
  }

  /**
   * Return the generators (variables) as polynomials.
   *
   * @returns Array of generator polynomials [x1, x2, ..., xn]
   *
   * @example
   * ```typescript
   * const R = new MPolynomialRing(F, ['x', 'y', 'z']);
   * const [x, y, z] = R.gens();
   * ```
   */
  gens(): MPolynomial<C>[] {
    if (this._gens === null) {
      this._gens = [];
      const one = this.base_ring.one() as C;

      for (let i = 0; i < this.ngens_value; i++) {
        const exp = new Array(this.ngens_value).fill(0);
        exp[i] = 1;
        const terms = new Map<string, C>();
        terms.set(exponentToKey(exp), one);
        this._gens.push(new MPolynomial(terms, this));
      }
    }
    return [...this._gens];
  }

  /**
   * Return the i-th generator.
   *
   * @param i - Index of the generator (0-based)
   * @returns The i-th generator polynomial
   */
  gen(i: number): MPolynomial<C> {
    if (i < 0 || i >= this.ngens_value) {
      throw new ValueError(`generator index ${i} out of range [0, ${this.ngens_value})`);
    }
    return this.gens()[i]!;
  }

  /**
   * Return the zero polynomial.
   */
  zero(): MPolynomial<C> {
    if (this._zero === null) {
      this._zero = new MPolynomial(new Map(), this);
    }
    return this._zero;
  }

  /**
   * Return the one polynomial (constant 1).
   */
  one(): MPolynomial<C> {
    if (this._one === null) {
      const one = this.base_ring.one() as C;
      const zeroExp = new Array(this.ngens_value).fill(0);
      const terms = new Map<string, C>();
      terms.set(exponentToKey(zeroExp), one);
      this._one = new MPolynomial(terms, this);
    }
    return this._one;
  }

  /**
   * Create a polynomial from various inputs.
   *
   * @param x - Can be:
   *   - A number/bigint (creates constant polynomial)
   *   - An element of the base ring (creates constant polynomial)
   *   - A Map from exponent keys to coefficients
   *   - An object mapping exponent tuples to coefficients
   *   - Another MPolynomial (converts/coerces)
   *
   * @example
   * ```typescript
   * // Constant polynomial
   * const c = R.__call__(5);
   *
   * // From coefficient element
   * const c2 = R.__call__(F(7n));
   *
   * // From dictionary
   * const f = R.__call__({ '2,0,0': F(1n), '0,1,1': F(2n) });
   * // Creates: x^2 + 2*y*z
   * ```
   */
  __call__(x: unknown): MPolynomial<C> {
    // Handle numbers and bigints as constants
    if (typeof x === 'number' || typeof x === 'bigint') {
      const coeff = this.base_ring.__call__(x) as C;
      if (coeff.isZero()) {
        return this.zero();
      }
      const zeroExp = new Array(this.ngens_value).fill(0);
      const terms = new Map<string, C>();
      terms.set(exponentToKey(zeroExp), coeff);
      return new MPolynomial(terms, this);
    }

    // Handle base ring elements as constants
    if (this.isBaseRingElement(x)) {
      const coeff = x as C;
      if (coeff.isZero()) {
        return this.zero();
      }
      const zeroExp = new Array(this.ngens_value).fill(0);
      const terms = new Map<string, C>();
      terms.set(exponentToKey(zeroExp), coeff);
      return new MPolynomial(terms, this);
    }

    // Handle MPolynomial from this ring
    if (x instanceof MPolynomial) {
      if (x.parent === this) {
        return x;
      }
      // Convert from another ring (simple coefficient conversion)
      const terms = new Map<string, C>();
      for (const [key, coeff] of x.monomial_coefficients()) {
        const newCoeff = this.base_ring.__call__(coeff) as C;
        if (!newCoeff.isZero()) {
          terms.set(key, newCoeff);
        }
      }
      return new MPolynomial(terms, this);
    }

    // Handle Map input
    if (x instanceof Map) {
      const terms = new Map<string, C>();
      for (const [key, coeff] of x) {
        const c = this.base_ring.__call__(coeff) as C;
        if (!c.isZero()) {
          // Ensure key is properly formatted
          const exp = keyToExponent(key as string);
          while (exp.length < this.ngens_value) {
            exp.push(0);
          }
          terms.set(exponentToKey(exp), c);
        }
      }
      return new MPolynomial(terms, this);
    }

    // Handle plain object as { 'exp': coeff } dictionary
    if (typeof x === 'object' && x !== null && !Array.isArray(x)) {
      const terms = new Map<string, C>();
      for (const [key, coeff] of Object.entries(x)) {
        const c = this.base_ring.__call__(coeff) as C;
        if (!c.isZero()) {
          // Ensure key is properly formatted
          const exp = keyToExponent(key);
          while (exp.length < this.ngens_value) {
            exp.push(0);
          }
          terms.set(exponentToKey(exp), c);
        }
      }
      return new MPolynomial(terms, this);
    }

    throw new ValueError(`cannot convert ${typeof x} to multivariate polynomial`);
  }

  /**
   * Create a monomial from an exponent tuple and optional coefficient.
   *
   * @param exponent - Array of exponents [e1, e2, ..., en]
   * @param coeff - Optional coefficient (default: 1)
   * @returns Monomial polynomial c * x1^e1 * x2^e2 * ... * xn^en
   *
   * @example
   * ```typescript
   * const R = new MPolynomialRing(F, ['x', 'y', 'z']);
   * const m = R.monomial([2, 1, 0], F(3n));  // 3*x^2*y
   * ```
   */
  monomial(exponent: number[], coeff?: C): MPolynomial<C> {
    const c = coeff ?? (this.base_ring.one() as C);
    if (c.isZero()) {
      return this.zero();
    }

    // Pad exponent to correct length
    const exp = [...exponent];
    while (exp.length < this.ngens_value) {
      exp.push(0);
    }

    const terms = new Map<string, C>();
    terms.set(exponentToKey(exp), c);
    return new MPolynomial(terms, this);
  }

  /**
   * Create a polynomial from terms (coefficient, exponent pairs).
   *
   * @param terms - Array of [coefficient, exponent] pairs
   * @returns Polynomial as sum of the terms
   *
   * @example
   * ```typescript
   * const R = new MPolynomialRing(F, ['x', 'y']);
   * const f = R.fromTerms([
   *   [F(1n), [2, 0]],  // x^2
   *   [F(2n), [1, 1]],  // 2*x*y
   *   [F(3n), [0, 0]],  // 3
   * ]);
   * // Creates: x^2 + 2*x*y + 3
   * ```
   */
  fromTerms(terms: Array<[C, number[]]>): MPolynomial<C> {
    const termMap = new Map<string, C>();

    for (const [coeff, exp] of terms) {
      if (coeff.isZero()) continue;

      // Pad exponent to correct length
      const paddedExp = [...exp];
      while (paddedExp.length < this.ngens_value) {
        paddedExp.push(0);
      }

      const key = exponentToKey(paddedExp);
      const existing = termMap.get(key);
      if (existing) {
        const sum = existing.add(coeff) as C;
        if (sum.isZero()) {
          termMap.delete(key);
        } else {
          termMap.set(key, sum);
        }
      } else {
        termMap.set(key, coeff);
      }
    }

    return new MPolynomial(termMap, this);
  }

  /**
   * Return the variable names as an array.
   */
  variable_names(): readonly string[] {
    return this.names;
  }

  /**
   * Return a dictionary mapping variable names to generators.
   *
   * @example
   * ```typescript
   * const { x, y, z } = R.gens_dict();
   * const f = x.mul(y).add(z);
   * ```
   */
  gens_dict(): Record<string, MPolynomial<C>> {
    const result: Record<string, MPolynomial<C>> = {};
    const generators = this.gens();
    for (let i = 0; i < this.names.length; i++) {
      result[this.names[i]!] = generators[i]!;
    }
    return result;
  }

  /**
   * Check if this ring is a field (always false for polynomial rings with variables).
   */
  is_field(): boolean {
    return false;
  }

  /**
   * Check if x is an element of the base ring.
   */
  private isBaseRingElement(x: unknown): boolean {
    if (x === null || x === undefined) return false;
    if (typeof x !== 'object') return false;

    // Check if it has the RingElement methods
    const obj = x as Record<string, unknown>;
    return (
      typeof obj.add === 'function' &&
      typeof obj.sub === 'function' &&
      typeof obj.mul === 'function' &&
      typeof obj.neg === 'function' &&
      typeof obj.eq === 'function' &&
      typeof obj.isZero === 'function'
    );
  }

  /**
   * String representation.
   */
  toString(): string {
    const varStr = this.names.join(', ');
    return `Multivariate Polynomial Ring in ${varStr} over ${this.base_ring}`;
  }
}

/**
 * Create a multivariate polynomial ring and return it with its generators.
 *
 * This is a convenience function similar to Sage's `PolynomialRing(R, names)`
 * that returns both the ring and its generators in a tuple.
 *
 * @param base_ring - The coefficient ring
 * @param names - Array of variable names
 * @param order - Term ordering (default: 'degrevlex')
 * @returns Tuple [ring, ...generators]
 *
 * @example
 * ```typescript
 * const [R, x, y, z] = MPolynomialRingConstructor(GF(p), ['x', 'y', 'z']);
 *
 * const f = x.pow(2).add(y.mul(z)).sub(R.one());  // x^2 + yz - 1
 * ```
 */
export function MPolynomialRingConstructor<C extends RingElement>(
  base_ring: CoefficientRing<C>,
  names: string[],
  order: TermOrder = 'degrevlex'
): [MPolynomialRing<C>, ...MPolynomial<C>[]] {
  const ring = new MPolynomialRing(base_ring, names, order);
  return [ring, ...ring.gens()];
}
