/**
 * @module sage/rings/polynomial/polynomial_commitment
 * @description Polynomial operations for polynomial commitment schemes (KZG, FRI)
 *
 * ADDITION, NOT A PORT. SageMath has no `sage/rings/polynomial/polynomial_commitment.py`
 * and none of the functions below correspond to an upstream SageMath function.
 * The module lives in the mirrored `rings/polynomial` tree only because every
 * operation is plain arithmetic in a univariate polynomial ring; the natural
 * home would be `src/zk/`, alongside `sumcheck` and `multilinear`.
 * Nothing here may be taken as evidence of SageMath behaviour.
 *
 * @see Deviation: polynomial_commitment.ts is an addition with no SageMath counterpart
 *
 * This module provides specialized polynomial operations used in zero-knowledge
 * proof systems, particularly for polynomial commitment schemes like KZG and FRI.
 *
 * The operations include:
 * - Quotient polynomial computation for opening proofs
 * - Batch quotient polynomials with vanishing polynomial division
 * - Random linear combinations (linearization)
 * - Lagrange basis evaluation with barycentric weights
 * - Polynomial splitting for multi-part commitments
 * - Polynomial composition with linear functions (FRI folding)
 */

import { ValueError, ZeroDivisionError } from '../../errors.js';
import type {
  CoefficientRing,
  Polynomial,
  PolynomialRingBase,
  RingElement,
} from './polynomial_element.js';

/**
 * Interface for field elements that support division and inversion.
 * Required for most polynomial commitment operations.
 */
export interface FieldElement extends RingElement {
  div(other: this): this;
  inv(): this;
}

/**
 * A point in the domain with its corresponding evaluation.
 */
export interface EvaluationPoint<C extends RingElement> {
  x: C;
  y: C;
}

/**
 * Compute the quotient polynomial q(X) = (f(X) - y) / (X - z)
 *
 * This is the core operation in KZG opening proofs. Given a polynomial f,
 * an evaluation point z, and the claimed value y = f(z), this computes
 * the quotient polynomial that proves the evaluation is correct.
 *
 * The division must be exact (no remainder), which is guaranteed when f(z) = y.
 *
 * @param f - The polynomial f(X)
 * @param z - The evaluation point
 * @param y - The claimed evaluation value (must equal f(z))
 * @returns The quotient polynomial q(X) such that f(X) - y = q(X) * (X - z)
 *
 * @throws {ValueError} If f(z) !== y (division has remainder)
 *
 * @example
 * ```typescript
 * // f(X) = X^2 - 1, z = 1, y = f(1) = 0
 * // q(X) = (X^2 - 1 - 0) / (X - 1) = X + 1
 * const q = compute_quotient(f, one, zero);
 * ```
 */
export function compute_quotient<C extends FieldElement>(
  f: Polynomial<C>,
  z: C,
  y: C
): Polynomial<C> {
  const ring = f.parent;

  // Compute f(X) - y
  const yPoly = ring.__call__(y);
  const numerator = f.sub(yPoly);

  // Create (X - z)
  const divisor = ring.gen().sub(ring.__call__(z));

  // Perform the division
  const [quotient, remainder] = numerator.quo_rem(divisor);

  // Check that the division is exact
  if (!remainder.isZero()) {
    throw new ValueError(
      'compute_quotient: f(z) != y, division has non-zero remainder. ' +
        `Expected f(${z}) = ${y}, but got non-zero remainder.`
    );
  }

  return quotient;
}

/**
 * Compute the batch quotient polynomial q(X) = (f(X) - I(X)) / Z_H(X)
 *
 * This is used for batch opening proofs where multiple points are opened at once.
 * I(X) is the interpolation polynomial through the given points, and Z_H(X)
 * is the vanishing polynomial on the domain.
 *
 * @param f - The polynomial f(X)
 * @param points - Array of (x, y) evaluation points where y = f(x)
 * @returns The quotient polynomial q(X)
 *
 * @throws {ValueError} If division has remainder (evaluations don't match f)
 *
 * @example
 * ```typescript
 * // f(X) that evaluates to y0 at x0, y1 at x1, y2 at x2
 * const points = [{ x: x0, y: y0 }, { x: x1, y: y1 }, { x: x2, y: y2 }];
 * const q = batch_quotient(f, points);
 * // q(X) = (f(X) - I(X)) / ((X - x0)(X - x1)(X - x2))
 * ```
 */
export function batch_quotient<C extends FieldElement>(
  f: Polynomial<C>,
  points: EvaluationPoint<C>[]
): Polynomial<C> {
  if (points.length === 0) {
    throw new ValueError('batch_quotient: points array cannot be empty');
  }

  const ring = f.parent;

  // Compute the vanishing polynomial Z_H(X) = product of (X - x_i)
  const vanishing = compute_vanishing_polynomial(
    ring,
    points.map((p) => p.x)
  );

  // Compute the interpolation polynomial I(X) through the points
  const interpolation = lagrange_interpolation(ring, points);

  // Compute f(X) - I(X)
  const numerator = f.sub(interpolation);

  // Divide by Z_H(X)
  const [quotient, remainder] = numerator.quo_rem(vanishing);

  if (!remainder.isZero()) {
    throw new ValueError('batch_quotient: f does not evaluate to the claimed values at all points');
  }

  return quotient;
}

/**
 * Compute the vanishing polynomial Z_H(X) = product of (X - x_i) for all x_i in domain.
 *
 * The vanishing polynomial is zero at all points in the domain.
 *
 * @param ring - The polynomial ring
 * @param domain - The domain points
 * @returns Z_H(X) = (X - x_0)(X - x_1)...(X - x_n)
 */
export function compute_vanishing_polynomial<C extends FieldElement>(
  ring: PolynomialRingBase<C>,
  domain: C[]
): Polynomial<C> {
  if (domain.length === 0) {
    return ring.one();
  }

  let result = ring.one();
  const x = ring.gen();

  for (const point of domain) {
    // (X - point)
    const factor = x.sub(ring.__call__(point));
    result = result.mul(factor);
  }

  return result;
}

/**
 * Compute Lagrange interpolation polynomial through the given points.
 *
 * Finds the unique polynomial I(X) of degree < n that passes through
 * all n given points.
 *
 * @param ring - The polynomial ring
 * @param points - Array of (x, y) points
 * @returns The interpolation polynomial I(X)
 *
 * @throws {ValueError} If domain points are not distinct
 */
export function lagrange_interpolation<C extends FieldElement>(
  ring: PolynomialRingBase<C>,
  points: EvaluationPoint<C>[]
): Polynomial<C> {
  if (points.length === 0) {
    return ring.zero();
  }

  if (points.length === 1) {
    return ring.__call__(points[0]!.y);
  }

  const n = points.length;
  let result = ring.zero();

  for (let i = 0; i < n; i++) {
    // Compute the i-th Lagrange basis polynomial L_i(X)
    const basis = lagrange_basis_polynomial(ring, points, i);
    // Add y_i * L_i(X) to the result
    result = result.add(basis.scalar_mul(points[i]!.y));
  }

  return result;
}

/**
 * Compute the i-th Lagrange basis polynomial L_i(X).
 *
 * L_i(X) = product_{j != i} (X - x_j) / (x_i - x_j)
 *
 * This polynomial equals 1 at x_i and 0 at all other x_j.
 *
 * @param ring - The polynomial ring
 * @param points - Array of evaluation points
 * @param i - Index of the basis polynomial
 * @returns L_i(X)
 */
function lagrange_basis_polynomial<C extends FieldElement>(
  ring: PolynomialRingBase<C>,
  points: EvaluationPoint<C>[],
  i: number
): Polynomial<C> {
  const x = ring.gen();
  const xi = points[i]!.x;
  let numerator = ring.one();
  let denominator = ring.base_ring.one() as C;

  for (let j = 0; j < points.length; j++) {
    if (j !== i) {
      const xj = points[j]!.x;
      // numerator *= (X - x_j)
      numerator = numerator.mul(x.sub(ring.__call__(xj)));
      // denominator *= (x_i - x_j)
      const diff = xi.sub(xj) as C;
      if (diff.isZero()) {
        throw new ValueError('lagrange_basis_polynomial: domain points must be distinct');
      }
      denominator = denominator.mul(diff) as C;
    }
  }

  // Divide the polynomial by the scalar denominator
  const denomInv = denominator.inv();
  return numerator.scalar_mul(denomInv);
}

/**
 * Compute a random linear combination of polynomials.
 *
 * L(X) = sum_{i} challenges[i] * polys[i](X)
 *
 * This is used in various ZK protocols to compress multiple polynomial
 * checks into a single check using random challenges.
 *
 * @param polys - Array of polynomials
 * @param challenges - Array of random challenge scalars
 * @returns The linear combination polynomial
 *
 * @throws {ValueError} If arrays have different lengths
 *
 * @example
 * ```typescript
 * // Combine polynomials f0, f1, f2 with challenges alpha, alpha^2, alpha^3
 * const L = linearization([f0, f1, f2], [alpha, alpha2, alpha3]);
 * // L(X) = alpha * f0(X) + alpha^2 * f1(X) + alpha^3 * f2(X)
 * ```
 */
export function linearization<C extends FieldElement>(
  polys: Polynomial<C>[],
  challenges: C[]
): Polynomial<C> {
  if (polys.length === 0) {
    throw new ValueError('linearization: polys array cannot be empty');
  }

  if (polys.length !== challenges.length) {
    throw new ValueError(
      'linearization: polys and challenges must have same length ' +
        `(got ${polys.length} and ${challenges.length})`
    );
  }

  const ring = polys[0]!.parent;
  let result = ring.zero();

  for (let i = 0; i < polys.length; i++) {
    result = result.add(polys[i]!.scalar_mul(challenges[i]!));
  }

  return result;
}

/**
 * Evaluate all Lagrange basis polynomials at a point.
 *
 * Returns [L_0(point), L_1(point), ..., L_{n-1}(point)] where L_i is the
 * i-th Lagrange basis polynomial for the given domain.
 *
 * This is useful for computing polynomial evaluations from coefficient
 * representations using the formula: f(point) = sum_i f_i * L_i(point)
 *
 * @param domain - The domain points (roots of unity or arbitrary points)
 * @param point - The evaluation point
 * @returns Array of Lagrange basis evaluations at the point
 *
 * @example
 * ```typescript
 * // For a domain of 4th roots of unity [1, omega, omega^2, omega^3]
 * const basis = evaluate_basis(domain, z);
 * // basis[i] = L_i(z) where L_i(omega^j) = (i == j) ? 1 : 0
 * ```
 */
export function evaluate_basis<C extends FieldElement>(domain: C[], point: C): C[] {
  const n = domain.length;
  if (n === 0) {
    return [];
  }

  // Use barycentric formula for efficiency
  const weights = barycentric_weights(domain);
  return evaluate_basis_with_weights(domain, weights, point);
}

/**
 * Evaluate Lagrange basis using precomputed barycentric weights.
 *
 * Uses the barycentric interpolation formula:
 * L_i(x) = (w_i / (x - x_i)) / sum_j (w_j / (x - x_j))
 *
 * @param domain - The domain points
 * @param weights - Precomputed barycentric weights
 * @param point - The evaluation point
 * @returns Array of Lagrange basis evaluations
 */
export function evaluate_basis_with_weights<C extends FieldElement>(
  domain: C[],
  weights: C[],
  point: C
): C[] {
  const n = domain.length;

  // Check if point is in the domain
  for (let i = 0; i < n; i++) {
    if (point.eq(domain[i]!)) {
      // Point is exactly x_i, so L_i = 1 and L_j = 0 for j != i
      const result: C[] = [];
      // Get 1 and 0 from the field using weights (which are always non-zero)
      const one = weights[0]!.mul(weights[0]!.inv()) as C;
      const zero = one.sub(one) as C;
      for (let j = 0; j < n; j++) {
        result.push(j === i ? one : zero);
      }
      return result;
    }
  }

  // Compute w_i / (point - x_i) for each i
  const terms: C[] = [];
  for (let i = 0; i < n; i++) {
    const diff = point.sub(domain[i]!) as C;
    terms.push(weights[i]!.div(diff));
  }

  // Compute the sum of all terms
  let sum = terms[0]!;
  for (let i = 1; i < n; i++) {
    sum = sum.add(terms[i]!) as C;
  }

  // L_i(point) = terms[i] / sum
  const sumInv = sum.inv();
  const result: C[] = [];
  for (let i = 0; i < n; i++) {
    result.push(terms[i]!.mul(sumInv) as C);
  }

  return result;
}

/**
 * Compute barycentric weights for a domain.
 *
 * The barycentric weight for x_i is:
 * w_i = 1 / product_{j != i} (x_i - x_j)
 *
 * These weights enable O(n) polynomial evaluation using the barycentric formula.
 *
 * @param domain - The domain points
 * @returns Array of barycentric weights
 *
 * @example
 * ```typescript
 * // For roots of unity, weights have a nice closed form
 * const weights = barycentric_weights(domain);
 * // weights[i] = 1 / product_{j != i} (omega^i - omega^j)
 * ```
 */
export function barycentric_weights<C extends FieldElement>(domain: C[]): C[] {
  const n = domain.length;
  if (n === 0) {
    return [];
  }

  // The multiplicative identity comes from the parent ring, so a one-point
  // domain such as [0] is legitimate: its weight is 1/(empty product) = 1.
  const one = getOne(domain[0]!);

  const weights: C[] = [];

  for (let i = 0; i < n; i++) {
    // Compute product_{j != i} (x_i - x_j)
    let product = one;

    for (let j = 0; j < n; j++) {
      if (j !== i) {
        const diff = domain[i]!.sub(domain[j]!) as C;
        if (diff.isZero()) {
          throw new ValueError('barycentric_weights: domain points must be distinct');
        }
        product = product.mul(diff) as C;
      }
    }

    // w_i = 1 / product
    weights.push(product.inv());
  }

  return weights;
}

/**
 * Helper to get the multiplicative identity of the field an element lives in.
 *
 * Prefer the parent ring's `one()` (SageMath's `x.parent().one()`); `x * x^-1`
 * is only a fallback for elements that do not expose a parent, and it fails
 * for zero — which made `generate_powers(0, n)` and `barycentric_weights([0])`
 * throw on perfectly legitimate inputs (audit L30).
 */
function getOne<C extends FieldElement>(element: C): C {
  const parent = (element as unknown as { parent?: { one?: () => unknown } }).parent;
  if (parent && typeof parent.one === 'function') {
    return parent.one() as C;
  }
  // a * a^(-1) = 1
  if (element.isZero()) {
    throw new ValueError('cannot get one from zero element');
  }
  return element.mul(element.inv()) as C;
}

/**
 * Get the verified degree bound of a polynomial.
 *
 * This returns the degree of the polynomial, which can be used to verify
 * that a polynomial commitment satisfies a claimed degree bound.
 *
 * @param poly - The polynomial
 * @returns The degree of the polynomial (-1 for zero polynomial)
 *
 * @example
 * ```typescript
 * // Verify that committed polynomial has degree < n
 * const deg = degree_bound(f);
 * assert(deg < n, 'Polynomial exceeds degree bound');
 * ```
 */
export function degree_bound<C extends RingElement>(poly: Polynomial<C>): number {
  return poly.degree();
}

/**
 * Check if a polynomial satisfies a degree bound.
 *
 * @param poly - The polynomial to check
 * @param bound - The degree bound (exclusive)
 * @returns true if deg(poly) < bound
 */
export function check_degree_bound<C extends RingElement>(
  poly: Polynomial<C>,
  bound: number
): boolean {
  return poly.degree() < bound;
}

/**
 * Split a polynomial into chunks for multi-part commitments.
 *
 * Given f(X) of degree d, split into ceil((d+1)/chunk_size) polynomials
 * f_0, f_1, ..., f_k such that:
 *   f(X) = f_0(X) + X^chunk_size * f_1(X) + X^(2*chunk_size) * f_2(X) + ...
 *
 * This is used when the polynomial is too large to commit to in a single
 * group element, requiring multiple commitments.
 *
 * @param poly - The polynomial to split
 * @param chunk_size - The size of each chunk (number of coefficients)
 * @returns Array of polynomial chunks
 *
 * @example
 * ```typescript
 * // f(X) = a0 + a1*X + a2*X^2 + a3*X^3 + a4*X^4 + a5*X^5
 * const chunks = split_poly(f, 2);
 * // chunks[0] = a0 + a1*X
 * // chunks[1] = a2 + a3*X
 * // chunks[2] = a4 + a5*X
 * // f(X) = chunks[0] + X^2 * chunks[1] + X^4 * chunks[2]
 * ```
 */
export function split_poly<C extends RingElement>(
  poly: Polynomial<C>,
  chunk_size: number
): Polynomial<C>[] {
  if (chunk_size <= 0) {
    throw new ValueError('split_poly: chunk_size must be positive');
  }

  if (poly.isZero()) {
    return [poly];
  }

  const ring = poly.parent;
  const chunks: Polynomial<C>[] = [];
  const coeffs = poly.coeffs;

  for (let start = 0; start < coeffs.length; start += chunk_size) {
    const end = Math.min(start + chunk_size, coeffs.length);
    const chunkCoeffs: C[] = [];

    for (let i = start; i < end; i++) {
      chunkCoeffs.push(coeffs[i] as C);
    }

    chunks.push(ring.__call__(chunkCoeffs));
  }

  return chunks;
}

/**
 * Recombine polynomial chunks back into the original polynomial.
 *
 * Given chunks f_0, f_1, ..., f_k, compute:
 *   f(X) = f_0(X) + X^chunk_size * f_1(X) + X^(2*chunk_size) * f_2(X) + ...
 *
 * @param chunks - Array of polynomial chunks
 * @param chunk_size - The size of each chunk
 * @returns The recombined polynomial
 */
export function recombine_chunks<C extends RingElement>(
  chunks: Polynomial<C>[],
  chunk_size: number
): Polynomial<C> {
  if (chunks.length === 0) {
    throw new ValueError('recombine_chunks: chunks array cannot be empty');
  }

  if (chunk_size <= 0) {
    throw new ValueError('recombine_chunks: chunk_size must be positive');
  }

  const ring = chunks[0]!.parent;
  let result = ring.zero();

  for (let i = chunks.length - 1; i >= 0; i--) {
    if (i < chunks.length - 1) {
      // Multiply current result by X^chunk_size
      result = result.shift(chunk_size);
    }
    result = result.add(chunks[i]!);
  }

  return result;
}

/**
 * Compute f(aX + b) efficiently.
 *
 * This operation is used in FRI (Fast Reed-Solomon Interactive Oracle Proof)
 * for the folding step, where a polynomial f(X) is transformed to f(alpha * X + beta)
 * for random challenges alpha and beta.
 *
 * Uses the Taylor expansion approach for efficiency:
 * If f(X) = sum_i c_i * X^i, then
 * f(aX + b) = sum_i c_i * (aX + b)^i
 *
 * @param f - The input polynomial f(X)
 * @param a - The linear coefficient
 * @param b - The constant offset
 * @returns The composed polynomial f(aX + b)
 *
 * @example
 * ```typescript
 * // FRI folding: combine even and odd coefficients
 * // f(X) = f_even(X^2) + X * f_odd(X^2)
 * // f(X) + f(-X) = 2 * f_even(X^2)
 * // At round r with challenge alpha:
 * // f_next(X) = f_even(X) + alpha * f_odd(X)
 *
 * const folded = compose_with_linear(f, alpha, beta);
 * ```
 */
export function compose_with_linear<C extends FieldElement>(
  f: Polynomial<C>,
  a: C,
  b: C
): Polynomial<C> {
  if (f.isZero()) {
    return f;
  }

  const ring = f.parent;
  const coeffs = f.coeffs;
  const n = coeffs.length;

  // Build the linear polynomial (aX + b)
  const linear = ring.__call__([b, a] as C[]);

  // Use Horner's method: f(aX + b) = c_n * (aX + b)^n + ... + c_1 * (aX + b) + c_0
  // = ((...((c_n * (aX + b) + c_{n-1}) * (aX + b) + c_{n-2}) * ...) * (aX + b) + c_0
  let result = ring.__call__(coeffs[n - 1]!);

  for (let i = n - 2; i >= 0; i--) {
    result = result.mul(linear).add(ring.__call__(coeffs[i]!));
  }

  return result;
}

/**
 * FRI-specific folding operation.
 *
 * Given a polynomial f(X) and a random challenge alpha, compute the
 * FRI folding:
 *   f_folded(X) = f_even(X) + alpha * f_odd(X)
 *
 * where f(X) = f_even(X^2) + X * f_odd(X^2)
 *
 * @param f - The polynomial to fold
 * @param alpha - The random folding challenge
 * @returns The folded polynomial of half the degree
 *
 * @example
 * ```typescript
 * // f(X) = a0 + a1*X + a2*X^2 + a3*X^3
 * // f_even(Y) = a0 + a2*Y
 * // f_odd(Y) = a1 + a3*Y
 * // fri_fold(f, alpha) = (a0 + alpha*a1) + (a2 + alpha*a3)*Y
 * ```
 */
export function fri_fold<C extends FieldElement>(f: Polynomial<C>, alpha: C): Polynomial<C> {
  if (f.isZero()) {
    return f;
  }

  const ring = f.parent;
  const coeffs = f.coeffs;

  // Split into even and odd coefficients
  const evenCoeffs: C[] = [];
  const oddCoeffs: C[] = [];

  for (let i = 0; i < coeffs.length; i++) {
    if (i % 2 === 0) {
      evenCoeffs.push(coeffs[i] as C);
    } else {
      oddCoeffs.push(coeffs[i] as C);
    }
  }

  // Compute f_even + alpha * f_odd
  const fEven = ring.__call__(evenCoeffs);
  const fOdd = ring.__call__(oddCoeffs);

  return fEven.add(fOdd.scalar_mul(alpha));
}

/**
 * Compute the derivative polynomial f'(X).
 *
 * Used in some polynomial commitment schemes for checking polynomial
 * identities.
 *
 * @param f - The polynomial
 * @returns The derivative f'(X)
 */
export function polynomial_derivative<C extends RingElement>(f: Polynomial<C>): Polynomial<C> {
  return f.derivative();
}

/**
 * Evaluate a polynomial at multiple points efficiently.
 *
 * For a polynomial f and points [z_0, z_1, ..., z_k], returns
 * [f(z_0), f(z_1), ..., f(z_k)].
 *
 * @param f - The polynomial
 * @param points - Array of evaluation points
 * @returns Array of evaluations
 */
export function multi_evaluate<C extends FieldElement>(f: Polynomial<C>, points: C[]): C[] {
  return points.map((z) => f.evaluate(z));
}

/**
 * Compute the polynomial that vanishes on roots of unity.
 *
 * For n-th roots of unity, the vanishing polynomial is X^n - 1.
 *
 * @param ring - The polynomial ring
 * @param n - The order of the roots of unity
 * @returns The vanishing polynomial X^n - 1
 */
export function roots_of_unity_vanishing<C extends FieldElement>(
  ring: PolynomialRingBase<C>,
  n: number
): Polynomial<C> {
  if (n <= 0) {
    throw new ValueError('roots_of_unity_vanishing: n must be positive');
  }

  // X^n - 1
  const one = ring.one();
  const xn = ring.monomial(n);
  return xn.sub(one);
}

/**
 * Compute a coset vanishing polynomial.
 *
 * For a coset k * H where H is the group of n-th roots of unity,
 * the vanishing polynomial is X^n - k^n.
 *
 * @param ring - The polynomial ring
 * @param n - The order of the roots of unity
 * @param coset_offset - The coset offset k
 * @returns The vanishing polynomial X^n - k^n
 */
export function coset_vanishing<C extends FieldElement>(
  ring: PolynomialRingBase<C>,
  n: number,
  coset_offset: C
): Polynomial<C> {
  if (n <= 0) {
    throw new ValueError('coset_vanishing: n must be positive');
  }

  // Compute k^n
  let kn = coset_offset;
  for (let i = 1; i < n; i++) {
    kn = kn.mul(coset_offset) as C;
  }

  // X^n - k^n
  const xn = ring.monomial(n);
  const knPoly = ring.__call__(kn);
  return xn.sub(knPoly);
}

/**
 * Verify that a quotient polynomial proof is valid.
 *
 * Checks that f(X) - y = q(X) * (X - z) by verifying:
 * 1. The polynomial equation holds
 * 2. deg(q) = deg(f) - 1
 *
 * @param f - The committed polynomial
 * @param z - The evaluation point
 * @param y - The claimed evaluation
 * @param q - The quotient polynomial (proof)
 * @returns true if the proof is valid
 */
export function verify_quotient_proof<C extends FieldElement>(
  f: Polynomial<C>,
  z: C,
  y: C,
  q: Polynomial<C>
): boolean {
  const ring = f.parent;

  // Create (X - z)
  const divisor = ring.gen().sub(ring.__call__(z));

  // Compute q(X) * (X - z) + y
  const reconstructed = q.mul(divisor).add(ring.__call__(y));

  // Check equality with f
  if (!reconstructed.eq(f)) {
    return false;
  }

  // Check degree bound: deg(q) should be deg(f) - 1
  if (f.degree() >= 0 && q.degree() !== f.degree() - 1) {
    // Special case: if f is constant and y = f, then q should be zero
    if (f.degree() === 0 && q.isZero()) {
      return true;
    }
    return false;
  }

  return true;
}

/**
 * Generate powers of a field element.
 *
 * Returns [1, alpha, alpha^2, ..., alpha^(n-1)]
 *
 * Useful for generating challenges in linearization.
 *
 * @param alpha - The base element
 * @param n - Number of powers to generate
 * @returns Array of powers
 */
export function generate_powers<C extends FieldElement>(alpha: C, n: number): C[] {
  if (n <= 0) {
    return [];
  }

  const one = getOne(alpha);
  const powers: C[] = [one];

  let current = alpha;
  for (let i = 1; i < n; i++) {
    powers.push(current);
    current = current.mul(alpha) as C;
  }

  return powers;
}

// Export all functions
export {
  compute_quotient as computeQuotient,
  batch_quotient as batchQuotient,
  compute_vanishing_polynomial as computeVanishingPolynomial,
  lagrange_interpolation as lagrangeInterpolation,
  evaluate_basis as evaluateBasis,
  barycentric_weights as barycentricWeights,
  degree_bound as degreeBound,
  check_degree_bound as checkDegreeBound,
  split_poly as splitPoly,
  recombine_chunks as recombineChunks,
  compose_with_linear as composeWithLinear,
  fri_fold as friFold,
  polynomial_derivative as polynomialDerivative,
  multi_evaluate as multiEvaluate,
  roots_of_unity_vanishing as rootsOfUnityVanishing,
  coset_vanishing as cosetVanishing,
  verify_quotient_proof as verifyQuotientProof,
  generate_powers as generatePowers,
};
