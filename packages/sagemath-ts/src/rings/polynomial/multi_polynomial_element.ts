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

import { NotImplementedError, ValueError, ZeroDivisionError } from '../../errors.js';
import type { CoefficientRing, RingElement } from './polynomial_element.js';

/**
 * Term ordering types for multivariate polynomials.
 * - 'lex': Lexicographic order (x > y > z)
 * - 'deglex': Total degree first, then lexicographic
 * - 'degrevlex': Total degree first, then reverse lexicographic (default)
 *
 * SageMath's `sage.rings.polynomial.term_order.TermOrder` supports twelve
 * orders plus block and weighted (`wdeglex`, `wdegrevlex`, ...) gradings.
 * Only the three global orders above are implemented here; every other name
 * is rejected at runtime with SageMath's message
 * (``unknown term order 'name'``) rather than silently falling back.
 *
 * @see Deviation: Polynomials — Printing, Factor Shape, Term Orders and Base Rings
 */
export type TermOrder = 'lex' | 'deglex' | 'degrevlex';

/**
 * The term orders this port understands, in SageMath's spelling.
 */
export const TERM_ORDERS: readonly TermOrder[] = Object.freeze([
  'lex',
  'deglex',
  'degrevlex',
] as TermOrder[]);

/**
 * Validate a term order name, raising SageMath's error for unknown names.
 *
 * SageMath: `term_order.py:796` -- ``raise ValueError("unknown term order {!r}".format(name))``.
 *
 * @param order - The candidate term order name
 * @returns The validated term order
 */
export function validateTermOrder(order: unknown): TermOrder {
  if (typeof order === 'string' && (TERM_ORDERS as readonly string[]).includes(order)) {
    return order as TermOrder;
  }
  throw new ValueError(`unknown term order '${String(order)}'`);
}

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
      // Do not silently fall back: an unrecognised order would otherwise
      // produce leading terms for a different ordering than requested.
      validateTermOrder(order);
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
   * For the zero polynomial, returns an array of zeros, matching SageMath
   * (`multi_polynomial_element.py:571`: `R(0).degrees()` is `(0, 0, 0, 0)`).
   * Note that this differs from `degreeIn`/`degree(x)`, which SageMath
   * deliberately defines as -1 for the zero polynomial.
   *
   * SageMath equivalent: `poly.degrees()`
   *
   * @returns Array of maximum degrees for each variable (0 to ngens-1)
   */
  degrees(): number[] {
    if (this.isZero()) {
      return new Array<number>(this.parent.ngens_value).fill(0);
    }
    const result: number[] = [];
    for (let i = 0; i < this.parent.ngens_value; i++) {
      result.push(this.degreeIn(i));
    }
    return result;
  }

  /**
   * Return the total degree of this polynomial.
   *
   * SageMath equivalent: `poly.total_degree()` (`multi_polynomial_element.py:708`).
   */
  total_degree(): number {
    return this.degree();
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
   * Return the coefficient in the base ring of the monomial `mon`.
   *
   * SageMath: `MPolynomial_polydict.monomial_coefficient`
   * (`multi_polynomial_element.py:748`) -- the argument is a monomial with the
   * same parent as `self` and the result lies in the base ring. This contrasts
   * with {@link coefficient}, which returns an element of the polynomial ring.
   *
   * As an extension we also accept the raw exponent tuple of the monomial,
   * which is the internal representation used throughout this port.
   *
   * @param mon - A monomial of the parent ring, or its exponent tuple
   * @returns The coefficient, or zero if not present
   */
  monomial_coefficient(mon: MPolynomial<C> | Exponent): C {
    let exponent: Exponent;
    if (mon instanceof MPolynomial) {
      if (mon.isZero()) {
        return this.parent.base_ring.zero() as C;
      }
      if (!mon.isTerm()) {
        throw new ValueError('input must be a monomial');
      }
      exponent = keyToExponent(mon.terms.keys().next().value as string);
    } else {
      exponent = mon;
    }

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
   * Return the coefficient of the variables with the degrees specified in
   * `degrees`, as an element of the parent polynomial ring.
   *
   * SageMath: `MPolynomial_polydict.coefficient`
   * (`multi_polynomial_element.py:946`), which delegates to
   * `PolyDict.polynomial_coefficient` (`polydict.pyx:531`): the terms whose
   * exponent matches every *restricted* variable are kept, the restricted
   * positions are zeroed, and the result is rebuilt in the same ring.
   *
   * ```
   * sage: R.<x, y> = QQ[]
   * sage: f = y^2 - x^9 - 7*x + 5*x*y
   * sage: f.coefficient({y: 1})
   * 5*x
   * sage: f.coefficient({y: 0})
   * -x^9 + (-7)*x
   * ```
   *
   * @param degrees - one of
   *   - a record mapping variable names to the required degree
   *     (SageMath keys this dictionary by the generators themselves; this port
   *     keys it by variable name, `@see Deviation: Multivariate coefficient() keys`),
   *   - a list of degree restrictions with `null` in the unrestricted positions,
   *   - a monomial of the parent ring (its positive exponents are the restrictions)
   * @returns Element of the parent ring
   */
  coefficient(
    degrees: MPolynomial<C> | (number | null)[] | Record<string, number>
  ): MPolynomial<C> {
    const n = this.parent.ngens_value;
    let looking_for: (number | null)[] | null = null;

    if (degrees instanceof MPolynomial) {
      if (degrees.parent === this.parent && degrees.isMonomial()) {
        const exp = keyToExponent(degrees.terms.keys().next().value as string);
        looking_for = [];
        for (let i = 0; i < n; i++) {
          const e = exp[i] ?? 0;
          looking_for.push(e > 0 ? e : null);
        }
      }
    } else if (Array.isArray(degrees)) {
      looking_for = [...degrees];
    } else if (typeof degrees === 'object' && degrees !== null) {
      looking_for = new Array<number | null>(n).fill(null);
      for (const [name, exp] of Object.entries(degrees)) {
        const i = this.parent.names.indexOf(name);
        if (i >= 0) {
          looking_for[i] = exp;
        }
      }
    }

    if (looking_for === null || looking_for.length === 0) {
      throw new ValueError('You must pass a dictionary list or monomial.');
    }

    // polydict.pyx:551-569
    const nz: number[] = [];
    for (let i = 0; i < looking_for.length; i++) {
      const d = looking_for[i];
      if (d !== null && d !== undefined) {
        nz.push(i);
      }
    }

    const ans = new Map<string, C>();
    for (const [key, coeff] of this.terms) {
      const exp = keyToExponent(key);
      let exactlyDivides = true;
      for (const j of nz) {
        if ((exp[j] ?? 0) !== looking_for[j]) {
          exactlyDivides = false;
          break;
        }
      }
      if (exactlyDivides) {
        const t = [...exp];
        while (t.length < n) t.push(0);
        for (const m of nz) {
          t[m] = 0;
        }
        const newKey = exponentToKey(t);
        const existing = ans.get(newKey);
        ans.set(newKey, existing ? (existing.add(coeff) as C) : coeff);
      }
    }

    return new MPolynomial(ans, this.parent);
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
    return exps.map((exp) => this.monomial_coefficient(exp));
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
   * Return the dictionary of this polynomial: exponent tuple -> coefficient.
   *
   * SageMath (`multi_polynomial_element.py:811`) returns a Python `dict`
   * keyed by `ETuple` exponent tuples. JavaScript `Map` keys are compared by
   * identity, so array keys would never collide correctly; this port keys the
   * map by the canonical comma-joined exponent string produced by
   * {@link exponentToKey} (use {@link keyToExponent} to recover the tuple).
   *
   * @see Deviation: Multivariate monomial_coefficients() key encoding
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
    return this.monomial_coefficient(leadingExp);
  }

  /**
   * Return the leading term (leading coefficient times leading monomial).
   */
  lt(): MPolynomial<C> {
    if (this.isZero()) {
      return this.parent.zero();
    }

    const leadingExp = this.leadingExponent();
    const coeff = this.monomial_coefficient(leadingExp);
    const terms = new Map<string, C>();
    terms.set(exponentToKey(leadingExp), coeff);
    return new MPolynomial(terms, this.parent);
  }

  /**
   * Return the leading term as a tuple [exponent, coefficient].
   *
   * This is a convenience method combining lm() and lc() information.
   * Returns null for the zero polynomial.
   *
   * @returns Tuple of [exponent array, coefficient] or null
   */
  leadingTerm(): [number[], C] | null {
    if (this.isZero()) {
      return null;
    }
    const exp = this.leadingExponent();
    const coeff = this.monomial_coefficient(exp);
    return [exp, coeff];
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
   * Return ``true`` if this polynomial is a generator of its parent.
   *
   * SageMath: `multi_polynomial_element.py:1358` -- a single term whose
   * exponent has exactly one nonzero entry, equal to 1, and whose coefficient
   * is one.
   */
  is_gen(): boolean {
    if (this.terms.size !== 1) {
      return false;
    }
    const entry = this.terms.entries().next().value as [string, C];
    const nonzero = keyToExponent(entry[0]).filter((e) => e !== 0);
    return nonzero.length === 1 && nonzero[0] === 1 && entry[1].eq(1);
  }

  /**
   * Return ``true`` if this multivariate polynomial is univariate.
   *
   * SageMath: `multi_polynomial_element.py:1543`. Constants (including zero)
   * are univariate.
   */
  is_univariate(): boolean {
    let found = -1;
    for (const key of this.terms.keys()) {
      const exp = keyToExponent(key);
      for (let i = 0; i < exp.length; i++) {
        if ((exp[i] ?? 0) === 0) continue;
        if (found !== i) {
          if (found !== -1) {
            return false;
          }
          found = i;
        }
      }
    }
    return true;
  }

  /**
   * Return the i-th variable occurring in this polynomial.
   *
   * SageMath: `multi_polynomial_element.py:1678` -- `self.variables()[i]`.
   */
  variable(i: number): MPolynomial<C> {
    const vars = this.variables();
    const v = vars[i];
    if (v === undefined) {
      throw new ValueError(`index ${i} out of range: this polynomial has ${vars.length} variables`);
    }
    return v;
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

  // ---------------------------------------------------------------------
  // Not-yet-implemented members of SageMath's MPolynomial interface.
  //
  // CLAUDE.md rule 7 requires every function of the mirrored module to exist,
  // raising NotImplementedError rather than being silently absent (callers
  // would otherwise get "x.factor is not a function"). Each stub names the
  // upstream location it must be ported from.
  // ---------------------------------------------------------------------

  /** SageMath: `multi_polynomial_element.py:2057` (delegates to Singular). */
  factor(_options?: { proof?: boolean }): never {
    throw new NotImplementedError('SAGE_NOT_IMPLEMENTED: MPolynomial.factor');
  }

  /** SageMath: `multi_polynomial_element.py:2262` (delegates to Singular). */
  quo_rem(_right: MPolynomial<C>): never {
    throw new NotImplementedError('SAGE_NOT_IMPLEMENTED: MPolynomial.quo_rem');
  }

  /** SageMath: `multi_polynomial_element.py:2221` (delegates to Singular). */
  lift(_I: unknown): never {
    throw new NotImplementedError('SAGE_NOT_IMPLEMENTED: MPolynomial.lift');
  }

  /** SageMath: `multi_polynomial_element.py:2423`; see `multi_polynomial_ideal.reduce`. */
  reduce(_I: unknown): never {
    throw new NotImplementedError('SAGE_NOT_IMPLEMENTED: MPolynomial.reduce');
  }

  /** SageMath: `multi_polynomial_element.py:2312` (delegates to Singular/Macaulay2). */
  resultant(_other: MPolynomial<C>, _variable?: MPolynomial<C>): never {
    throw new NotImplementedError('SAGE_NOT_IMPLEMENTED: MPolynomial.resultant');
  }

  /** SageMath: `multi_polynomial_element.py:2391` (delegates to Singular). */
  subresultants(_other: MPolynomial<C>, _variable?: MPolynomial<C>): never {
    throw new NotImplementedError('SAGE_NOT_IMPLEMENTED: MPolynomial.subresultants');
  }

  /** SageMath: `multi_polynomial_element.py:1962`. */
  integral(_var?: MPolynomial<C>): never {
    throw new NotImplementedError('SAGE_NOT_IMPLEMENTED: MPolynomial.integral');
  }

  /** SageMath: `multi_polynomial_element.py:1575`. */
  univariate_polynomial(_R?: unknown): never {
    throw new NotImplementedError('SAGE_NOT_IMPLEMENTED: MPolynomial.univariate_polynomial');
  }

  /** SageMath: `multi_polynomial_element.py:1276`. */
  inverse_of_unit(): never {
    throw new NotImplementedError('SAGE_NOT_IMPLEMENTED: MPolynomial.inverse_of_unit');
  }

  /** SageMath: `multi_polynomial_element.py:1036`. */
  global_height(_prec?: number): never {
    throw new NotImplementedError('SAGE_NOT_IMPLEMENTED: MPolynomial.global_height');
  }

  /** SageMath: `multi_polynomial_element.py:1116`. */
  local_height(_v: unknown, _prec?: number): never {
    throw new NotImplementedError('SAGE_NOT_IMPLEMENTED: MPolynomial.local_height');
  }

  /** SageMath: `multi_polynomial_element.py:1166`. */
  local_height_arch(_i: number, _prec?: number): never {
    throw new NotImplementedError('SAGE_NOT_IMPLEMENTED: MPolynomial.local_height_arch');
  }

  /** SageMath: `multi_polynomial.pyx` -- `MPolynomial.gcd` (delegates to Singular). */
  gcd(_other: MPolynomial<C>): never {
    throw new NotImplementedError('SAGE_NOT_IMPLEMENTED: MPolynomial.gcd');
  }

  /** SageMath: `multi_polynomial.pyx` -- `MPolynomial.lcm` (delegates to Singular). */
  lcm(_other: MPolynomial<C>): never {
    throw new NotImplementedError('SAGE_NOT_IMPLEMENTED: MPolynomial.lcm');
  }

  /** SageMath: `multi_polynomial.pyx` -- `MPolynomial.homogenize`. */
  homogenize(_var?: unknown): never {
    throw new NotImplementedError('SAGE_NOT_IMPLEMENTED: MPolynomial.homogenize');
  }

  /** SageMath: `multi_polynomial.pyx` -- `MPolynomial.numerator`. */
  numerator(): never {
    throw new NotImplementedError('SAGE_NOT_IMPLEMENTED: MPolynomial.numerator');
  }

  /** SageMath: `multi_polynomial.pyx` -- `MPolynomial.denominator`. */
  denominator(): never {
    throw new NotImplementedError('SAGE_NOT_IMPLEMENTED: MPolynomial.denominator');
  }

  /** SageMath: `multi_polynomial.pyx` -- `MPolynomial.is_squarefree`. */
  is_squarefree(): never {
    throw new NotImplementedError('SAGE_NOT_IMPLEMENTED: MPolynomial.is_squarefree');
  }

  /** SageMath: `multi_polynomial.pyx` -- `MPolynomial.is_unit`. */
  is_unit(): never {
    throw new NotImplementedError('SAGE_NOT_IMPLEMENTED: MPolynomial.is_unit');
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
      const coeff = this.monomial_coefficient(exp);
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

/**
 * Type alias for MPolynomial, used in multi_polynomial_ideal.
 * This provides compatibility with code expecting MPolynomialElement.
 */
export type MPolynomialElement<R extends CoefficientRing, E extends RingElement> = MPolynomial<E>;
