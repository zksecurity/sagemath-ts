/**
 * @module sage/rings/polynomial/polynomial_ring
 * @description Polynomial rings over arbitrary coefficient rings
 *
 * Port of: sage/rings/polynomial/polynomial_ring.py
 */

import { factor } from '../../arith/misc.js';
import { ArithmeticError, TypeError, ValueError } from '../../errors.js';
import {
  type CoefficientRing,
  Polynomial,
  type PolynomialRingBase,
  type RingElement,
} from './polynomial_element.js';

/**
 * Interface for field elements that support division/inverse.
 */
interface FieldElement extends RingElement {
  inv(): this;
  div?(other: this): this;
}

/**
 * A polynomial ring R[x] over a base ring R.
 */
export class PolynomialRing<C extends RingElement> implements PolynomialRingBase<C> {
  readonly base_ring: CoefficientRing<C>;
  readonly variable_name: string;

  constructor(base_ring: CoefficientRing<C>, variable_name: string = 'x') {
    this.base_ring = base_ring;
    this.variable_name = variable_name;
  }

  /**
   * Create a polynomial from coefficients or another polynomial.
   */
  __call__(x: C | C[] | Polynomial<C> | number | bigint): Polynomial<C> {
    if (x instanceof Polynomial) {
      // Convert polynomial from potentially different ring
      return new Polynomial(
        x.coeffs.map((c) => this.base_ring.__call__(c) as C),
        this
      );
    }

    if (Array.isArray(x)) {
      return new Polynomial(x, this);
    }

    // Single coefficient (constant polynomial)
    const coeff = this.base_ring.__call__(x) as C;
    return new Polynomial([coeff], this);
  }

  /**
   * Return the zero polynomial.
   */
  zero(): Polynomial<C> {
    return new Polynomial([], this);
  }

  /**
   * Return the one polynomial (constant 1).
   */
  one(): Polynomial<C> {
    return new Polynomial([this.base_ring.one() as C], this);
  }

  /**
   * Return the generator (the variable x).
   */
  gen(): Polynomial<C> {
    return new Polynomial([this.base_ring.zero() as C, this.base_ring.one() as C], this);
  }

  /**
   * Create a polynomial from a list of (coefficient, degree) pairs.
   */
  fromTerms(terms: Array<[C, number]>): Polynomial<C> {
    if (terms.length === 0) {
      return this.zero();
    }

    const maxDeg = Math.max(...terms.map(([_, d]) => d));
    const coeffs: C[] = [];

    for (let i = 0; i <= maxDeg; i++) {
      coeffs.push(this.base_ring.zero() as C);
    }

    for (const [c, d] of terms) {
      coeffs[d] = coeffs[d]!.add(c) as C;
    }

    return new Polynomial(coeffs, this);
  }

  /**
   * Return the monomial x^n.
   */
  monomial(n: number): Polynomial<C> {
    const coeffs: C[] = [];
    for (let i = 0; i < n; i++) {
      coeffs.push(this.base_ring.zero() as C);
    }
    coeffs.push(this.base_ring.one() as C);
    return new Polynomial(coeffs, this);
  }

  /**
   * Check if this is a polynomial ring (yes).
   */
  is_ring(): boolean {
    return true;
  }

  /**
   * Check if this is a field (no, unless base is trivial).
   */
  is_field(): boolean {
    return false;
  }

  /**
   * Return the Lagrange interpolation polynomial through the given points.
   *
   * This computes the unique polynomial P of degree at most n-1 such that
   * P(x_i) = y_i for all given points (x_i, y_i).
   *
   * With `algorithm='neville'` the **whole last row of the Neville table** is
   * returned (a list of polynomials), exactly as Sage does; the interpolating
   * polynomial is its last entry.
   *
   * @param points - Array of [x, y] pairs where the x values must be distinct
   *   (`x_i - x_j` must be invertible); as in Sage there is no precheck, a
   *   repeated x surfaces as a {@link ZeroDivisionError} from the base ring.
   * @param algorithm - `'divided_difference'` (default), `'neville'` or `'pari'`
   * @param previous_row - only used with `'neville'`: the last row of a
   *   previous Neville computation, so that it can be extended incrementally
   * @returns The interpolating polynomial, or the last Neville row
   *
   * @example
   * ```typescript
   * const [R, x] = PolynomialRingConstructor(F7, 'x');
   * const p = R.lagrange_polynomial([[F7(0), F7(1)], [F7(2), F7(2)], [F7(3), F7(6)]]);
   * // p(0) = 1, p(2) = 2, p(3) = 6
   * ```
   *
   * @see Reference: sage/rings/polynomial/polynomial_ring.py:2295 (lagrange_polynomial)
   */
  lagrange_polynomial(
    points: Array<[C, C]>,
    algorithm?: 'divided_difference' | 'pari'
  ): Polynomial<C>;
  lagrange_polynomial(
    points: Array<[C, C]>,
    algorithm: 'neville',
    previous_row?: Polynomial<C>[]
  ): Polynomial<C>[];
  lagrange_polynomial(
    points: Array<[C, C]>,
    algorithm: 'divided_difference' | 'neville' | 'pari' = 'divided_difference',
    previous_row?: Polynomial<C>[]
  ): Polynomial<C> | Polynomial<C>[] {
    const n = points.length;

    if (algorithm === 'divided_difference') {
      if (n === 0) {
        return this.zero();
      }
      return this._lagrange_divided_difference(points);
    }

    if (algorithm === 'neville') {
      return this._lagrange_neville(points, previous_row);
    }

    if (algorithm === 'pari') {
      // PARI's polinterpolate is Newton's divided-difference scheme, which is
      // what `_lagrange_divided_difference` implements.
      if (n === 0) {
        return this.zero();
      }
      return this._lagrange_divided_difference(points);
    }

    throw new ValueError("algorithm can be 'divided_difference', 'neville' or 'pari'");
  }

  /**
   * Compute Lagrange interpolation using divided differences.
   *
   * Uses Newton's form with divided differences, evaluated using
   * Horner's method for efficiency.
   */
  private _lagrange_divided_difference(points: Array<[C, C]>): Polynomial<C> {
    const n = points.length;
    if (n === 0) {
      return this.zero();
    }

    // Compute divided differences
    const F = this.divided_difference(points);

    // Evaluate in nested form (Horner's method)
    // P(x) = F[n-1] + (x - x_{n-2}) * (F[n-2] + (x - x_{n-3}) * (...))
    let P = this.__call__(F[n - 1]!);
    const x = this.gen();

    for (let i = n - 2; i >= 0; i--) {
      // P = P * (x - x_i) + F[i]
      const xi = this.__call__(points[i]![0]);
      P = P.mul(x.sub(xi));
      P = P.add(this.__call__(F[i]!));
    }

    return P;
  }

  /**
   * Return the Newton divided-difference coefficients.
   *
   * These are the coefficients F[0,0], F[1,1], ..., F[n-1,n-1] such that
   * P(x) = sum_{i=0}^{n-1} F[i,i] * prod_{j=0}^{i-1} (x - x_j)
   *
   * @param points - Array of [x, y] pairs
   * @param full_table - if true, return the full divided-difference table
   *   instead of only its main diagonal
   * @returns Array of divided difference coefficients (or the full table)
   *
   * @see Reference: sage/rings/polynomial/polynomial_ring.py:2206 (divided_difference)
   */
  divided_difference(points: Array<[C, C]>, full_table?: false): C[];
  divided_difference(points: Array<[C, C]>, full_table: true): C[][];
  divided_difference(points: Array<[C, C]>, full_table: boolean = false): C[] | C[][] {
    const n = points.length;
    if (n === 0) {
      return [];
    }

    // F[i][j] stores the divided difference f[x_{i-j}, ..., x_i]
    const F: C[][] = [];
    for (let i = 0; i < n; i++) {
      F.push([points[i]![1]]);
    }

    for (let i = 1; i < n; i++) {
      for (let j = 1; j <= i; j++) {
        const numer = F[i]![j - 1]!.sub(F[i - 1]![j - 1]!) as C;
        const denom = points[i]![0].sub(points[i - j]![0]) as C;
        const quotient = divideElements(numer, denom);
        F[i]!.push(quotient);
      }
    }

    if (full_table) {
      return F;
    }

    // Return diagonal elements F[i][i]
    return F.map((row, i) => row[i]!);
  }

  /**
   * Compute Lagrange interpolation using Neville's method.
   *
   * Returns the last row of the Neville table, where the last element
   * is the interpolating polynomial.
   *
   * @param points - Array of [x, y] pairs
   * @param previousRow - Optional previous row for incremental computation
   * @returns Array of polynomials (last row of Neville table)
   */
  private _lagrange_neville(points: Array<[C, C]>, previousRow?: Polynomial<C>[]): Polynomial<C>[] {
    const N = points.length;
    const M = previousRow?.length ?? 0;

    // P keeps track of the previous row, Q keeps track of the current row
    const P: Polynomial<C>[] = previousRow
      ? [...previousRow, ...new Array(N - M).fill(null)]
      : new Array(N).fill(null);
    const Q: Polynomial<C>[] = new Array(N).fill(null);

    const x = this.gen();

    for (let i = M; i < N; i++) {
      // Start populating the current row
      Q[0] = this.__call__(points[i]![1]);

      for (let j = 1; j <= i; j++) {
        // Q[j] = ((x - x_{i-j}) * Q[j-1] - (x - x_i) * P[j-1]) / (x_i - x_{i-j})
        const xiMinusJ = this.__call__(points[i - j]![0]);
        const xi = this.__call__(points[i]![0]);

        const numer = x
          .sub(xiMinusJ)
          .mul(Q[j - 1]!)
          .sub(x.sub(xi).mul(P[j - 1]!));

        const denom = points[i]![0].sub(points[i - j]![0]) as C;

        // Divide polynomial by scalar
        Q[j] = scalarDividePolynomial(numer, denom);
      }

      // Swap P and Q for next iteration
      for (let k = 0; k <= i; k++) {
        P[k] = Q[k]!;
      }
    }

    return P.filter((p) => p !== null) as Polynomial<C>[];
  }

  /**
   * Newton interpolation polynomial through the given points.
   *
   * This is more efficient than Lagrange when adding new points incrementally.
   * Uses divided differences internally.
   *
   * @param points - Array of [x, y] pairs where x values must be distinct
   * @returns The interpolating polynomial in Newton form
   *
   * @example
   * ```typescript
   * const [R, x] = PolynomialRingConstructor(F7, 'x');
   * const p = R.newton_interpolation([[F7(0), F7(1)], [F7(1), F7(3)], [F7(2), F7(7)]]);
   * ```
   *
   * @see Reference: sage/rings/polynomial/polynomial_ring.py:divided_difference
   */
  newton_interpolation(points: Array<[C, C]>): Polynomial<C> {
    // Newton interpolation is the same as Lagrange using divided differences
    // The difference is mainly conceptual - Newton form is easier to extend
    return this.lagrange_polynomial(points, 'divided_difference');
  }

  /**
   * Return the vanishing polynomial for a given domain.
   *
   * The vanishing polynomial is Z_H(X) = prod_{h in H} (X - h),
   * which evaluates to 0 for all elements in the domain H.
   *
   * @param domain - Array of field elements
   * @returns The vanishing polynomial
   *
   * @example
   * ```typescript
   * const [R, x] = PolynomialRingConstructor(F7, 'x');
   * const Z = R.vanishing_polynomial([F7(1), F7(2), F7(3)]);
   * // Z = (x - 1)(x - 2)(x - 3)
   * // Z(1) = Z(2) = Z(3) = 0
   * ```
   */
  vanishing_polynomial(domain: C[]): Polynomial<C> {
    if (domain.length === 0) {
      return this.one();
    }

    // Compute product (X - h) for all h in domain
    const x = this.gen();
    let result = this.one();

    for (const h of domain) {
      const factor = x.sub(this.__call__(h));
      result = result.mul(factor);
    }

    return result;
  }

  /**
   * Evaluate the interpolating polynomial at a point using barycentric interpolation.
   *
   * This is more efficient than constructing the full polynomial when you only
   * need the value at a single point. Complexity is O(n) for evaluation after
   * O(n) preprocessing, vs O(n^2) for full polynomial construction.
   *
   * Uses the second form of the barycentric interpolation formula:
   * L(x) = (sum_{j=0}^{n-1} w_j * y_j / (x - x_j)) / (sum_{j=0}^{n-1} w_j / (x - x_j))
   *
   * where w_j = 1 / prod_{k != j} (x_j - x_k) are the barycentric weights.
   *
   * @param points - Array of [x, y] pairs where x values must be distinct
   * @param evalPoint - The point at which to evaluate
   * @returns The value of the interpolating polynomial at evalPoint
   *
   * @throws {ValueError} If evalPoint coincides with one of the x values
   *
   * @example
   * ```typescript
   * const [R, x] = PolynomialRingConstructor(F7, 'x');
   * const val = R.barycentric_interpolation([[F7(0), F7(1)], [F7(1), F7(2)]], F7(3));
   * ```
   */
  barycentric_interpolation(points: Array<[C, C]>, evalPoint: C): C {
    const n = points.length;

    if (n === 0) {
      return this.base_ring.zero() as C;
    }

    // Check if evalPoint coincides with any x value
    for (let j = 0; j < n; j++) {
      if (evalPoint.eq(points[j]![0])) {
        // Return the corresponding y value directly
        return points[j]![1];
      }
    }

    // Compute barycentric weights: w_j = 1 / prod_{k != j} (x_j - x_k)
    const weights: C[] = [];
    for (let j = 0; j < n; j++) {
      let w = this.base_ring.one() as C;
      for (let k = 0; k < n; k++) {
        if (k !== j) {
          const diff = points[j]![0].sub(points[k]![0]) as C;
          w = divideElements(w, diff);
        }
      }
      weights.push(w);
    }

    // Compute L(x) using barycentric formula
    let numerator = this.base_ring.zero() as C;
    let denominator = this.base_ring.zero() as C;

    for (let j = 0; j < n; j++) {
      const diff = evalPoint.sub(points[j]![0]) as C;
      const term = divideElements(weights[j]!, diff);

      numerator = numerator.add(term.mul(points[j]![1]) as C) as C;
      denominator = denominator.add(term) as C;
    }

    return divideElements(numerator, denominator);
  }

  /**
   * Create a polynomial with the given roots.
   *
   * Returns the monic polynomial prod_{r in roots} (X - r).
   *
   * @param roots - Array of roots
   * @returns The polynomial with the given roots
   *
   * @example
   * ```typescript
   * const [R, x] = PolynomialRingConstructor(F7, 'x');
   * const p = R.from_roots([F7(1), F7(2), F7(3)]);
   * // p = (x - 1)(x - 2)(x - 3) = x^3 - 6x^2 + 11x - 6
   * // p(1) = p(2) = p(3) = 0
   * ```
   */
  from_roots(roots: C[]): Polynomial<C> {
    // This is the same as the vanishing polynomial
    return this.vanishing_polynomial(roots);
  }

  /**
   * Return the n-th cyclotomic polynomial.
   *
   * The n-th cyclotomic polynomial Phi_n(x) is the minimal polynomial over Q
   * of a primitive n-th root of unity. Its roots are exactly the primitive
   * n-th roots of unity.
   *
   * The polynomial is computed using the formula:
   * Phi_n(x) = prod_{d|n} (x^d - 1)^{mu(n/d)}
   *
   * where mu is the Moebius function.
   *
   * @param n - A positive integer
   * @returns The n-th cyclotomic polynomial
   *
   * @throws {ValueError} If n is not a positive integer
   *
   * @example
   * ```typescript
   * const [R, x] = PolynomialRingConstructor(F7, 'x');
   * const phi6 = R.cyclotomic_polynomial(6);
   * // Phi_6(x) = x^2 - x + 1
   * ```
   *
   * @see Reference: sage/rings/polynomial/cyclotomic.pyx:cyclotomic_coeffs
   */
  cyclotomic_polynomial(n: number): Polynomial<C> {
    if (n <= 0) {
      // sage/rings/polynomial/polynomial_ring.py:1188
      throw new ArithmeticError(`n=${n} must be positive`);
    }
    if (!Number.isInteger(n)) {
      throw new TypeError(`n=${n} must be an integer`);
    }

    // Compute the cyclotomic polynomial coefficients
    const coeffs = cyclotomicCoeffs(n);

    // Convert integer coefficients to ring elements
    const ringCoeffs: C[] = coeffs.map((c) => this.base_ring.__call__(c) as C);

    return new Polynomial(ringCoeffs, this);
  }

  toString(): string {
    return `Univariate Polynomial Ring in ${this.variable_name} over ${this.base_ring}`;
  }
}

/**
 * Divide two field elements.
 */
function divideElements<C extends RingElement>(a: C, b: C): C {
  // Try to call div method if it exists
  if ('div' in a && typeof (a as unknown as { div: (b: C) => C }).div === 'function') {
    return (a as unknown as { div: (b: C) => C }).div(b);
  }

  // Try inv method for field elements
  if ('inv' in b && typeof (b as unknown as { inv: () => C }).inv === 'function') {
    const bInv = (b as unknown as { inv: () => C }).inv();
    return a.mul(bInv) as C;
  }

  throw new ValueError('coefficient ring does not support division');
}

/**
 * Divide a polynomial by a scalar.
 */
function scalarDividePolynomial<C extends RingElement>(p: Polynomial<C>, c: C): Polynomial<C> {
  if (c.isZero()) {
    throw new ValueError('division by zero');
  }

  const cInv = divideElements(p.parent.base_ring.one() as C, c);
  return p.scalar_mul(cInv);
}

/**
 * Compute the coefficients of the n-th cyclotomic polynomial.
 *
 * Uses the relation: x^n - 1 = prod_{d|n} Phi_d(x)
 * So: Phi_n(x) = (x^n - 1) / prod_{d|n, d<n} Phi_d(x)
 *
 * For prime powers p^k: Phi_{p^k}(x) = Phi_p(x^{p^{k-1}})
 * For non-squarefree n = m * rad where rad is the radical:
 *   Phi_n(x) = Phi_rad(x^{n/rad})
 *
 * @param n - A positive integer
 * @returns Array of integer coefficients [c_0, c_1, ..., c_deg]
 */
function cyclotomicCoeffs(n: number): number[] {
  if (n === 1) {
    return [-1, 1]; // x - 1
  }

  // Factor n
  const factorization = primeFactorization(n);

  // Check for non-squarefree case first
  const isSquarefree = factorization.every(([_, exp]) => exp === 1);

  if (!isSquarefree) {
    // For non-squarefree n: Phi_n(x) = Phi_rad(x^{n/rad})
    // where rad is the radical (product of distinct prime factors)
    const rad = factorization.map(([p, _]) => p).reduce((a, b) => a * b, 1);
    const pow = n / rad;
    const phiRad = cyclotomicCoeffs(rad);

    // Substitute x^pow for x
    const degree = (phiRad.length - 1) * pow;
    const result = new Array(degree + 1).fill(0);
    for (let i = 0; i < phiRad.length; i++) {
      result[i * pow] = phiRad[i]!;
    }
    return result;
  }

  // Now n is squarefree
  // For prime p: Phi_p(x) = 1 + x + x^2 + ... + x^{p-1}
  if (factorization.length === 1) {
    const p = factorization[0]![0];
    return new Array(p).fill(1);
  }

  // For squarefree n with multiple prime factors, use recursion:
  // Phi_n(x) = (x^n - 1) / prod_{d|n, d<n} Phi_d(x)

  // Get all proper divisors of n (divisors < n)
  const divisors = getDivisors(n).filter((d) => d < n);

  // Start with x^n - 1
  let result = new Array(n + 1).fill(0);
  result[0] = -1;
  result[n] = 1;

  // Divide by Phi_d for each proper divisor d
  for (const d of divisors) {
    const phiD = cyclotomicCoeffs(d);
    result = polynomialDivide(result, phiD);
  }

  return result;
}

/**
 * Get all divisors of n.
 */
function getDivisors(n: number): number[] {
  const divisors: number[] = [];
  for (let i = 1; i * i <= n; i++) {
    if (n % i === 0) {
      divisors.push(i);
      if (i !== n / i) {
        divisors.push(n / i);
      }
    }
  }
  return divisors.sort((a, b) => a - b);
}

/**
 * Divide polynomial a by polynomial b (exact division).
 * Assumes b divides a exactly with no remainder.
 */
function polynomialDivide(a: number[], b: number[]): number[] {
  // Remove trailing zeros
  while (a.length > 1 && a[a.length - 1] === 0) {
    a.pop();
  }
  while (b.length > 1 && b[b.length - 1] === 0) {
    b.pop();
  }

  if (b.length === 0 || (b.length === 1 && b[0] === 0)) {
    throw new Error('Division by zero polynomial');
  }

  const degA = a.length - 1;
  const degB = b.length - 1;

  if (degA < degB) {
    return [0];
  }

  const degQ = degA - degB;
  const quotient = new Array(degQ + 1).fill(0);
  const remainder = [...a];

  const lcB = b[degB]!;

  for (let i = degQ; i >= 0; i--) {
    const coeff = remainder[i + degB]! / lcB;
    quotient[i] = coeff;

    for (let j = 0; j <= degB; j++) {
      remainder[i + j] -= coeff * b[j]!;
    }
  }

  // Remove trailing zeros from quotient
  while (quotient.length > 1 && quotient[quotient.length - 1] === 0) {
    quotient.pop();
  }

  return quotient;
}

/**
 * Get the prime factorization of n as [prime, exponent] pairs.
 *
 * Uses PARI's factorization via the factor() function.
 */
function primeFactorization(n: number): Array<[number, number]> {
  if (n <= 1) {
    return [];
  }

  // Use PARI's factorization and convert to number pairs
  const factorization = factor(BigInt(n));
  return factorization
    .filter(([p, _]) => p > 0n) // Exclude -1 sign factor
    .map(([p, e]) => [Number(p), Number(e)]);
}

/**
 * Create a polynomial ring over a base ring.
 *
 * @param base_ring - The coefficient ring
 * @param names - Variable name(s)
 * @returns A polynomial ring and its generator
 *
 * @example
 * ```typescript
 * const [R, x] = PolynomialRingConstructor(GF2, 'x');
 * const p = x.pow(2).add(x).add(R.one());  // x^2 + x + 1
 * ```
 */
export function PolynomialRingConstructor<C extends RingElement>(
  base_ring: CoefficientRing<C>,
  names: string = 'x'
): [PolynomialRing<C>, Polynomial<C>] {
  const ring = new PolynomialRing(base_ring, names);
  return [ring, ring.gen()];
}
