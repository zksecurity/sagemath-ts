/**
 * @module sage/rings/polynomial/multi_polynomial_element
 * @description Multivariate polynomial elements over arbitrary coefficient rings
 *
 * Port of: sage/rings/polynomial/multi_polynomial_element.py
 *          sage/rings/polynomial/multi_polynomial.pyx
 *
 * Provides sparse representation of multivariate polynomials suitable for
 * ZK constraint systems.
 */

import { ValueError, ZeroDivisionError } from '../../errors.js';
import type { CoefficientRing, RingElement } from './polynomial_element.js';

/**
 * Term ordering types for multivariate polynomials.
 * - 'lex': Lexicographic order (x > y > z)
 * - 'deglex': Total degree first, then lexicographic
 * - 'degrevlex': Total degree first, then reverse lexicographic (default)
 */
export type TermOrder = 'lex' | 'deglex' | 'degrevlex';

/**
 * Exponent tuple represented as an array of non-negative integers.
 * For a polynomial in variables x, y, z:
 * - [2, 1, 0] represents x^2*y
 * - [0, 0, 3] represents z^3
 */
export type Exponent = readonly number[];

/**
 * Convert exponent tuple to a string key for Map storage.
 */
export function exponentToKey(exp: Exponent): string {
  return exp.join(',');
}

/**
 * Convert string key back to exponent tuple.
 */
export function keyToExponent(key: string): number[] {
  if (key === '') return [];
  return key.split(',').map(Number);
}

/**
 * Compute the total degree of an exponent tuple.
 */
export function totalDegree(exp: Exponent): number {
  return exp.reduce((sum, e) => sum + e, 0);
}

/**
 * Compare two exponents using lexicographic order.
 * Returns negative if a < b, positive if a > b, 0 if equal.
 */
function compareLex(a: Exponent, b: Exponent): number {
  const n = Math.max(a.length, b.length);
  for (let i = 0; i < n; i++) {
    const ai = a[i] ?? 0;
    const bi = b[i] ?? 0;
    if (ai !== bi) {
      return ai - bi;
    }
  }
  return 0;
}

/**
 * Compare two exponents using degree lexicographic order.
 * First compares total degree, then uses lexicographic order.
 */
function compareDegLex(a: Exponent, b: Exponent): number {
  const degA = totalDegree(a);
  const degB = totalDegree(b);
  if (degA !== degB) {
    return degA - degB;
  }
  return compareLex(a, b);
}

/**
 * Compare two exponents using degree reverse lexicographic order.
 * First compares total degree, then uses reverse lexicographic order
 * (comparing from the last variable first, with smaller being "greater").
 */
function compareDegRevLex(a: Exponent, b: Exponent): number {
  const degA = totalDegree(a);
  const degB = totalDegree(b);
  if (degA !== degB) {
    return degA - degB;
  }
  // Reverse lexicographic: compare from end, smaller exponent is "greater"
  const n = Math.max(a.length, b.length);
  for (let i = n - 1; i >= 0; i--) {
    const ai = a[i] ?? 0;
    const bi = b[i] ?? 0;
    if (ai !== bi) {
      // Note: reversed comparison - smaller is "greater"
      return bi - ai;
    }
  }
  return 0;
}

/**
 * Get the comparison function for a given term order.
 */
export function getTermOrderComparator(order: TermOrder): (a: Exponent, b: Exponent) => number {
  switch (order) {
    case 'lex':
      return compareLex;
    case 'deglex':
      return compareDegLex;
    case 'degrevlex':
      return compareDegRevLex;
    default:
      return compareDegRevLex;
  }
}

/**
 * Add two exponent tuples (for multiplication of monomials).
 */
export function addExponents(a: Exponent, b: Exponent): number[] {
  const n = Math.max(a.length, b.length);
  const result: number[] = [];
  for (let i = 0; i < n; i++) {
    result.push((a[i] ?? 0) + (b[i] ?? 0));
  }
  return result;
}

/**
 * Subtract two exponent tuples (for division of monomials).
 * Returns null if result would have negative exponents.
 */
export function subtractExponents(a: Exponent, b: Exponent): number[] | null {
  const n = Math.max(a.length, b.length);
  const result: number[] = [];
  for (let i = 0; i < n; i++) {
    const diff = (a[i] ?? 0) - (b[i] ?? 0);
    if (diff < 0) {
      return null;
    }
    result.push(diff);
  }
  return result;
}

/**
 * Check if exponent a is divisible by exponent b.
 */
export function exponentDivides(a: Exponent, b: Exponent): boolean {
  const n = Math.max(a.length, b.length);
  for (let i = 0; i < n; i++) {
    if ((a[i] ?? 0) < (b[i] ?? 0)) {
      return false;
    }
  }
  return true;
}

/**
 * Check if an exponent is all zeros.
 */
export function isZeroExponent(exp: Exponent): boolean {
  return exp.every((e) => e === 0);
}

/**
 * Interface for multivariate polynomial rings.
 */
export interface MPolynomialRingBase<C extends RingElement> {
  readonly base_ring: CoefficientRing<C>;
  readonly names: readonly string[];
  readonly ngens_value: number;
  readonly term_order: TermOrder;

  zero(): MPolynomial<C>;
  one(): MPolynomial<C>;
  gens(): MPolynomial<C>[];
  ngens(): number;
  __call__(x: unknown): MPolynomial<C>;
  gen(i: number): MPolynomial<C>;
}

/**
 * A multivariate polynomial with coefficients in a ring R.
 *
 * Uses sparse representation with a Map from exponent keys to coefficients.
 * This is efficient for sparse polynomials common in ZK constraint systems.
 *
 * @example
 * ```typescript
 * const R = new MPolynomialRing(GF(p), ['x', 'y', 'z']);
 * const [x, y, z] = R.gens();
 * const f = x.mul(x).add(y.mul(z)).add(R.one());  // x^2 + yz + 1
 * ```
 */
export class MPolynomial<C extends RingElement> {
  /**
   * Sparse representation: maps exponent keys to coefficients.
   * Keys are comma-separated exponent values, e.g., "2,1,0" for x^2*y.
   */
  private readonly terms: Map<string, C>;

  /**
   * The parent polynomial ring.
   */
  readonly parent: MPolynomialRingBase<C>;

  /**
   * Create a multivariate polynomial.
   *
   * @param terms - Map from exponent keys to coefficients
   * @param parent - The parent polynomial ring
   */
  constructor(terms: Map<string, C>, parent: MPolynomialRingBase<C>) {
    this.parent = parent;

    // Remove zero coefficients
    this.terms = new Map();
    for (const [key, coeff] of terms) {
      if (!coeff.isZero()) {
        this.terms.set(key, coeff);
      }
    }
  }

  /**
   * Return the number of terms (monomials with non-zero coefficients).
   */
  numberOfTerms(): number {
    return this.terms.size;
  }

  /**
   * Check if this is the zero polynomial.
   */
  isZero(): boolean {
    return this.terms.size === 0;
  }

  /**
   * Check if this is a constant polynomial.
   */
  isConstant(): boolean {
    if (this.terms.size === 0) return true;
    if (this.terms.size > 1) return false;

    const key = this.terms.keys().next().value;
    if (key === undefined) return true;
    return isZeroExponent(keyToExponent(key));
  }

  /**
   * Check if this polynomial is a single term (monomial times coefficient).
   */
  isTerm(): boolean {
    return this.terms.size === 1;
  }

  /**
   * Check if this polynomial is a monomial (single term with coefficient 1).
   */
  isMonomial(): boolean {
    if (this.terms.size !== 1) return false;
    const coeff = this.terms.values().next().value;
    return coeff?.eq(1);
  }

  /**
   * Return the total degree of this polynomial.
   * The zero polynomial has degree -1.
   */
  degree(): number {
    if (this.terms.size === 0) {
      return -1;
    }

    let maxDeg = -1;
    for (const key of this.terms.keys()) {
      const exp = keyToExponent(key);
      const deg = totalDegree(exp);
      if (deg > maxDeg) {
        maxDeg = deg;
      }
    }
    return maxDeg;
  }

  /**
   * Return the degree in a specific variable.
   *
   * @param varIndex - The index of the variable (0-based)
   * @returns The maximum power of that variable, or -1 for zero polynomial
   */
  degreeIn(varIndex: number): number {
    if (this.terms.size === 0) {
      return -1;
    }

    let maxDeg = 0;
    for (const key of this.terms.keys()) {
      const exp = keyToExponent(key);
      const deg = exp[varIndex] ?? 0;
      if (deg > maxDeg) {
        maxDeg = deg;
      }
    }
    return maxDeg;
  }

  /**
   * Return the maximum degree in each variable as an array.
   *
   * For a polynomial f = x^2*y + y^3*z, returns [2, 3, 1].
   * For the zero polynomial, returns an array of -1s.
   *
   * SageMath equivalent: `poly.degrees()`
   *
   * @returns Array of maximum degrees for each variable (0 to ngens-1)
   */
  degrees(): number[] {
    const result: number[] = [];
    for (let i = 0; i < this.parent.ngens_value; i++) {
      result.push(this.degreeIn(i));
    }
    return result;
  }

  /**
   * Return the variables that actually appear in this polynomial.
   *
   * For a polynomial f = x^2 + z, returns [x, z] (not y).
   *
   * SageMath equivalent: `poly.variables()`
   *
   * @returns Array of polynomial ring generators that appear with non-zero degree
   */
  variables(): MPolynomial<C>[] {
    if (this.isZero()) {
      return [];
    }

    const gens = this.parent.gens();
    const result: MPolynomial<C>[] = [];

    for (let i = 0; i < this.parent.ngens_value; i++) {
      const deg = this.degreeIn(i);
      if (deg > 0) {
        result.push(gens[i]!);
      }
    }

    return result;
  }

  /**
   * Return the number of variables that actually appear in this polynomial.
   *
   * SageMath equivalent: `len(poly.variables())`
   *
   * @returns Count of variables with non-zero degree
   */
  nvariables(): number {
    if (this.isZero()) {
      return 0;
    }

    let count = 0;
    for (let i = 0; i < this.parent.ngens_value; i++) {
      if (this.degreeIn(i) > 0) {
        count++;
      }
    }
    return count;
  }

  /**
   * Return all generators of the parent polynomial ring.
   *
   * This is an alias for `this.parent.gens()`, provided for SageMath
   * compatibility where `poly.args()` is commonly used to get the
   * variable list in sumcheck protocols.
   *
   * SageMath equivalent: `poly.args()`
   *
   * @returns Array of all ring generators (variables)
   */
  args(): MPolynomial<C>[] {
    return this.parent.gens();
  }

  /**
   * Substitute values for some variables.
   *
   * This is an alias for `partialEvaluate()`, provided for SageMath
   * compatibility where `.subs({x: 5, y: 3})` syntax is common.
   *
   * SageMath equivalent: `poly.subs({x: 5, y: 3})`
   *
   * @param values - An object mapping variable names to values
   * @returns A new polynomial with those variables substituted
   */
  subs(values: Record<string, C>): MPolynomial<C> {
    return this.partialEvaluate(values);
  }

  /**
   * Return the constant coefficient (the term with all zero exponents).
   */
  constantCoefficient(): C {
    const zeroKey = exponentToKey(new Array(this.parent.ngens_value).fill(0));
    const coeff = this.terms.get(zeroKey);
    return coeff ?? (this.parent.base_ring.zero() as C);
  }

  /**
   * Return the coefficient of a given monomial.
   *
   * @param exponent - The exponent tuple
   * @returns The coefficient, or zero if not present
   */
  coefficient(exponent: Exponent): C {
    // Pad exponent to match number of variables
    const paddedExp = [...exponent];
    while (paddedExp.length < this.parent.ngens_value) {
      paddedExp.push(0);
    }

    const key = exponentToKey(paddedExp);
    const coeff = this.terms.get(key);
    return coeff ?? (this.parent.base_ring.zero() as C);
  }

  /**
   * Return all monomials in this polynomial, sorted by term order.
   *
   * @returns Array of polynomials, each consisting of a single monomial with coefficient 1
   */
  monomials(): MPolynomial<C>[] {
    const exps = this.exponents();
    const one = this.parent.base_ring.one() as C;

    return exps.map((exp) => {
      const terms = new Map<string, C>();
      terms.set(exponentToKey(exp), one);
      return new MPolynomial(terms, this.parent);
    });
  }

  /**
   * Return the coefficients of this polynomial, sorted by term order.
   */
  coefficients(): C[] {
    const exps = this.exponents();
    return exps.map((exp) => this.coefficient(exp));
  }

  /**
   * Return the exponents of all monomials, sorted by term order (descending).
   */
  exponents(): number[][] {
    const comparator = getTermOrderComparator(this.parent.term_order);
    const exps = Array.from(this.terms.keys()).map(keyToExponent);
    // Sort in descending order (leading term first)
    exps.sort((a, b) => -comparator(a, b));
    return exps;
  }

  /**
   * Return the internal terms map (for advanced use).
   */
  monomial_coefficients(): Map<string, C> {
    return new Map(this.terms);
  }

  /**
   * Return the leading monomial (with respect to term order).
   *
   * @returns A polynomial consisting of the leading monomial with coefficient 1
   */
  lm(): MPolynomial<C> {
    if (this.isZero()) {
      return this.parent.zero();
    }

    const leadingExp = this.leadingExponent();
    const one = this.parent.base_ring.one() as C;
    const terms = new Map<string, C>();
    terms.set(exponentToKey(leadingExp), one);
    return new MPolynomial(terms, this.parent);
  }

  /**
   * Return the leading coefficient (with respect to term order).
   */
  lc(): C {
    if (this.isZero()) {
      return this.parent.base_ring.zero() as C;
    }

    const leadingExp = this.leadingExponent();
    return this.coefficient(leadingExp);
  }

  /**
   * Return the leading term (leading coefficient times leading monomial).
   */
  lt(): MPolynomial<C> {
    if (this.isZero()) {
      return this.parent.zero();
    }

    const leadingExp = this.leadingExponent();
    const coeff = this.coefficient(leadingExp);
    const terms = new Map<string, C>();
    terms.set(exponentToKey(leadingExp), coeff);
    return new MPolynomial(terms, this.parent);
  }

  /**
   * Get the leading exponent according to term order.
   */
  private leadingExponent(): number[] {
    const comparator = getTermOrderComparator(this.parent.term_order);
    let leading: number[] | null = null;

    for (const key of this.terms.keys()) {
      const exp = keyToExponent(key);
      if (leading === null || comparator(exp, leading) > 0) {
        leading = exp;
      }
    }

    return leading ?? new Array(this.parent.ngens_value).fill(0);
  }

  /**
   * Add two multivariate polynomials.
   */
  add(other: MPolynomial<C>): MPolynomial<C> {
    const result = new Map<string, C>(this.terms);

    for (const [key, coeff] of other.terms) {
      const existing = result.get(key);
      if (existing) {
        const sum = existing.add(coeff) as C;
        if (sum.isZero()) {
          result.delete(key);
        } else {
          result.set(key, sum);
        }
      } else {
        result.set(key, coeff);
      }
    }

    return new MPolynomial(result, this.parent);
  }

  /**
   * Subtract two multivariate polynomials.
   */
  sub(other: MPolynomial<C>): MPolynomial<C> {
    return this.add(other.neg());
  }

  /**
   * Negate this polynomial.
   */
  neg(): MPolynomial<C> {
    const result = new Map<string, C>();
    for (const [key, coeff] of this.terms) {
      result.set(key, coeff.neg() as C);
    }
    return new MPolynomial(result, this.parent);
  }

  /**
   * Multiply two multivariate polynomials.
   */
  mul(other: MPolynomial<C>): MPolynomial<C> {
    if (this.isZero() || other.isZero()) {
      return this.parent.zero();
    }

    const result = new Map<string, C>();

    for (const [keyA, coeffA] of this.terms) {
      const expA = keyToExponent(keyA);

      for (const [keyB, coeffB] of other.terms) {
        const expB = keyToExponent(keyB);
        const newExp = addExponents(expA, expB);
        const newKey = exponentToKey(newExp);
        const newCoeff = coeffA.mul(coeffB) as C;

        const existing = result.get(newKey);
        if (existing) {
          const sum = existing.add(newCoeff) as C;
          if (sum.isZero()) {
            result.delete(newKey);
          } else {
            result.set(newKey, sum);
          }
        } else if (!newCoeff.isZero()) {
          result.set(newKey, newCoeff);
        }
      }
    }

    return new MPolynomial(result, this.parent);
  }

  /**
   * Multiply by a scalar (element of the base ring).
   */
  scalarMul(c: C): MPolynomial<C> {
    if (c.isZero()) {
      return this.parent.zero();
    }

    const result = new Map<string, C>();
    for (const [key, coeff] of this.terms) {
      const prod = coeff.mul(c) as C;
      if (!prod.isZero()) {
        result.set(key, prod);
      }
    }
    return new MPolynomial(result, this.parent);
  }

  /**
   * Compute this^n.
   */
  pow(n: number | bigint): MPolynomial<C> {
    let exp = typeof n === 'bigint' ? n : BigInt(n);

    if (exp < 0n) {
      throw new ValueError('negative exponent not supported for polynomials');
    }

    if (exp === 0n) {
      return this.parent.one();
    }

    // Binary exponentiation
    let result = this.parent.one();
    let base: MPolynomial<C> = this;

    while (exp > 0n) {
      if ((exp & 1n) === 1n) {
        result = result.mul(base);
      }
      base = base.mul(base);
      exp >>= 1n;
    }

    return result;
  }

  /**
   * Evaluate this polynomial at a given point.
   *
   * @param values - An object mapping variable names to values, or an array of values
   * @returns The result of substituting all variables
   */
  evaluate(values: Record<string, C> | C[]): C {
    let result = this.parent.base_ring.zero() as C;

    // Convert to array form if given as object
    let valueArray: C[];
    if (Array.isArray(values)) {
      valueArray = values;
    } else {
      valueArray = this.parent.names.map(
        (name) => values[name] ?? (this.parent.base_ring.zero() as C)
      );
    }

    for (const [key, coeff] of this.terms) {
      const exp = keyToExponent(key);

      // Compute the monomial value: x1^e1 * x2^e2 * ...
      let monomialValue = this.parent.base_ring.one() as C;
      for (let i = 0; i < exp.length; i++) {
        const e = exp[i]!;
        if (e > 0) {
          const xi = valueArray[i] ?? (this.parent.base_ring.zero() as C);
          // Compute xi^e
          let power = this.parent.base_ring.one() as C;
          for (let j = 0; j < e; j++) {
            power = power.mul(xi) as C;
          }
          monomialValue = monomialValue.mul(power) as C;
        }
      }

      result = result.add(coeff.mul(monomialValue) as C) as C;
    }

    return result;
  }

  /**
   * Partially evaluate this polynomial by substituting some variables.
   *
   * @param values - An object mapping variable names to values (only some variables)
   * @returns A new polynomial with those variables substituted
   */
  partialEvaluate(values: Record<string, C>): MPolynomial<C> {
    // Find which variable indices are being substituted
    const substituteIndices = new Set<number>();
    for (let i = 0; i < this.parent.names.length; i++) {
      if (this.parent.names[i]! in values) {
        substituteIndices.add(i);
      }
    }

    const result = new Map<string, C>();

    for (const [key, coeff] of this.terms) {
      const exp = keyToExponent(key);

      // Compute the coefficient from substituted variables
      let substCoeff = coeff;
      for (const i of substituteIndices) {
        const e = exp[i] ?? 0;
        if (e > 0) {
          const xi = values[this.parent.names[i]!]!;
          // Compute xi^e and multiply
          let power = this.parent.base_ring.one() as C;
          for (let j = 0; j < e; j++) {
            power = power.mul(xi) as C;
          }
          substCoeff = substCoeff.mul(power) as C;
        }
      }

      // Create new exponent with substituted variables set to 0
      const newExp: number[] = [];
      for (let i = 0; i < this.parent.ngens_value; i++) {
        if (substituteIndices.has(i)) {
          newExp.push(0);
        } else {
          newExp.push(exp[i] ?? 0);
        }
      }

      const newKey = exponentToKey(newExp);

      const existing = result.get(newKey);
      if (existing) {
        const sum = existing.add(substCoeff) as C;
        if (sum.isZero()) {
          result.delete(newKey);
        } else {
          result.set(newKey, sum);
        }
      } else if (!substCoeff.isZero()) {
        result.set(newKey, substCoeff);
      }
    }

    return new MPolynomial(result, this.parent);
  }

  /**
   * Check equality with another polynomial.
   */
  eq(other: MPolynomial<C> | number | bigint): boolean {
    if (typeof other === 'number' || typeof other === 'bigint') {
      if (!this.isConstant()) return false;
      const constCoeff = this.constantCoefficient();
      return constCoeff.eq(other);
    }

    if (this.terms.size !== other.terms.size) {
      return false;
    }

    for (const [key, coeff] of this.terms) {
      const otherCoeff = other.terms.get(key);
      if (!otherCoeff || !coeff.eq(otherCoeff)) {
        return false;
      }
    }

    return true;
  }

  /**
   * Return the derivative with respect to a variable.
   *
   * @param varIndex - The index of the variable to differentiate by
   * @returns The partial derivative
   */
  derivative(varIndex: number): MPolynomial<C> {
    const result = new Map<string, C>();

    for (const [key, coeff] of this.terms) {
      const exp = keyToExponent(key);
      const e = exp[varIndex] ?? 0;

      if (e > 0) {
        // New exponent with one less power in this variable
        const newExp = [...exp];
        newExp[varIndex] = e - 1;

        // Multiply coefficient by e
        let newCoeff = coeff;
        for (let i = 1; i < e; i++) {
          newCoeff = newCoeff.add(coeff) as C;
        }

        if (!newCoeff.isZero()) {
          const newKey = exponentToKey(newExp);
          const existing = result.get(newKey);
          if (existing) {
            const sum = existing.add(newCoeff) as C;
            if (!sum.isZero()) {
              result.set(newKey, sum);
            } else {
              result.delete(newKey);
            }
          } else {
            result.set(newKey, newCoeff);
          }
        }
      }
    }

    return new MPolynomial(result, this.parent);
  }

  /**
   * Check if this polynomial is homogeneous (all terms have the same total degree).
   */
  isHomogeneous(): boolean {
    if (this.terms.size <= 1) return true;

    let degree: number | null = null;
    for (const key of this.terms.keys()) {
      const exp = keyToExponent(key);
      const deg = totalDegree(exp);
      if (degree === null) {
        degree = deg;
      } else if (deg !== degree) {
        return false;
      }
    }
    return true;
  }

  /**
   * Return a string representation of this polynomial.
   */
  toString(): string {
    if (this.terms.size === 0) {
      return '0';
    }

    const exps = this.exponents();
    const terms: string[] = [];

    for (const exp of exps) {
      const coeff = this.coefficient(exp);
      const coeffStr = coeff.toString();
      const isOne = coeff.eq(1);
      const isMinusOne = coeff.eq(-1);
      const isNegative = coeffStr.startsWith('-');

      // Build monomial string
      const monomialParts: string[] = [];
      for (let i = 0; i < exp.length; i++) {
        const e = exp[i]!;
        if (e > 0) {
          const varName = this.parent.names[i] ?? `x${i}`;
          if (e === 1) {
            monomialParts.push(varName);
          } else {
            monomialParts.push(`${varName}^${e}`);
          }
        }
      }
      const monomialStr = monomialParts.join('*');

      let termStr: string;
      if (monomialStr === '') {
        // Constant term
        termStr = coeffStr;
      } else if (isOne) {
        termStr = monomialStr;
      } else if (isMinusOne) {
        termStr = `-${monomialStr}`;
      } else {
        // Check if coefficient needs parentheses
        const needsParens =
          coeffStr.includes('+') || (coeffStr.includes('-') && !coeffStr.startsWith('-'));
        if (needsParens) {
          termStr = `(${coeffStr})*${monomialStr}`;
        } else {
          termStr = `${coeffStr}*${monomialStr}`;
        }
      }

      terms.push(termStr);
    }

    // Join with + and handle negative terms
    let result = terms[0]!;
    for (let i = 1; i < terms.length; i++) {
      const term = terms[i]!;
      if (term.startsWith('-')) {
        result += ` - ${term.slice(1)}`;
      } else {
        result += ` + ${term}`;
      }
    }

    return result;
  }
}
