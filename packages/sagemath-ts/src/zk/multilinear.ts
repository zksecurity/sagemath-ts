/**
 * @module sage/zk/multilinear
 * @description Multi-Linear Extensions (MLEs) and Boolean Hypercube Utilities
 *
 * MLEs are the foundation of many modern proof systems (GKR, sumcheck, etc.)
 * An MLE uniquely extends a function f: {0,1}^n -> F to F^n -> F
 * such that the polynomial is multilinear (degree at most 1 in each variable).
 *
 * @example
 * ```typescript
 * import { GF, multilinearExtension, booleanHypercube, eqPolynomial } from 'sagemath-ts';
 *
 * const F = GF(101n);
 * const values = [F(9n), F(2n), F(5n), F(4n)];
 *
 * // Interpolate the MLE from values on the boolean hypercube
 * const poly = multilinearExtension(values, F);
 *
 * // Verify: poly evaluates back to the original values
 * for (let i = 0; i < values.length; i++) {
 *   const bits = intToBinary(i, 2);
 *   const evalPt = { x0: F(bits[0]), x1: F(bits[1]) };
 *   console.log(poly.evaluate(evalPt)); // matches values[i]
 * }
 * ```
 */

import { ValueError } from '../errors.js';
import type { FiniteFieldElement } from '../rings/finite_rings/finite_field_prime.js';
import type { MPolynomial, MPolynomialRingBase } from '../rings/polynomial/multi_polynomial_element.js';
import type { MPolynomialRing } from '../rings/polynomial/multi_polynomial_ring.js';
import { MPolynomialRingConstructor } from '../rings/polynomial/multi_polynomial_ring.js';
import type { CoefficientRing, RingElement } from '../rings/polynomial/polynomial_element.js';

// ============================================
// 1. Helper functions for binary/integer conversion
// ============================================

/**
 * Compute ceil(log2(n)).
 *
 * Returns the smallest integer k such that 2^k >= n.
 * This determines the number of variables needed for an MLE
 * interpolating n values.
 *
 * @param n - A positive integer
 * @returns ceil(log2(n)), or 0 if n <= 1
 *
 * @example
 * ```typescript
 * closestPowerOfTwo(1);   // 0
 * closestPowerOfTwo(2);   // 1
 * closestPowerOfTwo(3);   // 2
 * closestPowerOfTwo(4);   // 2
 * closestPowerOfTwo(5);   // 3
 * closestPowerOfTwo(100); // 7
 * ```
 */
export function closestPowerOfTwo(n: number): number {
  if (n <= 1) return 0;
  return Math.ceil(Math.log2(n));
}

/**
 * Convert an integer to its binary representation as an array of bits.
 *
 * Returns bits in big-endian order: [most significant, ..., least significant].
 *
 * @param val - The integer value to convert (non-negative)
 * @param numBits - The number of bits in the output (zero-padded if needed)
 * @returns Array of 0s and 1s representing the binary value
 *
 * @example
 * ```typescript
 * intToBinary(5, 4);  // [0, 1, 0, 1]
 * intToBinary(5, 8);  // [0, 0, 0, 0, 0, 1, 0, 1]
 * intToBinary(0, 3);  // [0, 0, 0]
 * intToBinary(7, 3);  // [1, 1, 1]
 * ```
 */
export function intToBinary(val: number, numBits: number): number[] {
  if (val < 0) {
    throw new ValueError('val must be non-negative');
  }
  if (numBits < 0) {
    throw new ValueError('numBits must be non-negative');
  }
  const binary = val.toString(2).padStart(numBits, '0');
  return binary.split('').map(Number);
}

/**
 * Convert a binary array back to an integer.
 *
 * Expects bits in big-endian order: [most significant, ..., least significant].
 *
 * @param bits - Array of 0s and 1s
 * @returns The integer value
 *
 * The result is a `bigint`: the bit array can be longer than 53 bits (the MLE
 * of a large table has one index per hypercube point), and a `number` would
 * silently lose precision there.
 *
 * @throws {ValueError} If any entry is not 0 or 1
 *
 * @example
 * ```typescript
 * binaryToInt([0, 1, 0, 1]); // 5n
 * binaryToInt([1, 1, 1]);    // 7n
 * binaryToInt([0, 0, 0]);    // 0n
 * binaryToInt([]);           // 0n
 * ```
 */
export function binaryToInt(bits: number[]): bigint {
  let result = 0n;
  for (const bit of bits) {
    if (bit !== 0 && bit !== 1) {
      throw new ValueError(`bits must contain only 0 and 1 (got ${bit})`);
    }
    result = (result << 1n) | BigInt(bit);
  }
  return result;
}

/**
 * Largest dimension accepted by {@link booleanHypercube}.
 *
 * 2^25 points already allocate roughly a gigabyte of small arrays; anything
 * beyond that is a programming error rather than a computation.
 */
export const MAX_HYPERCUBE_DIM = 25;

/**
 * Generate all points in the n-dimensional boolean hypercube {0,1}^n.
 *
 * Returns points in lexicographic order by their binary representation,
 * i.e., [0,...,0], [0,...,1], ..., [1,...,1].
 *
 * This is the port of SageMath's `Tuples([0, 1], n)` as used by the sumcheck
 * blueprint (`reference/sage_blueprints/sumcheck.sage:116`). Materialising the
 * whole hypercube costs 2^n arrays, so dimensions beyond
 * {@link MAX_HYPERCUBE_DIM} are rejected instead of exhausting memory.
 *
 * @param n - The dimension (number of coordinates per point)
 * @returns Array of 2^n points, each point being an array of n bits
 * @throws {ValueError} If n is negative or larger than {@link MAX_HYPERCUBE_DIM}
 *
 * @example
 * ```typescript
 * booleanHypercube(2);
 * // [[0, 0], [0, 1], [1, 0], [1, 1]]
 *
 * booleanHypercube(3);
 * // [[0,0,0], [0,0,1], [0,1,0], [0,1,1], [1,0,0], [1,0,1], [1,1,0], [1,1,1]]
 * ```
 */
export function booleanHypercube(n: number): number[][] {
  if (n < 0) {
    throw new ValueError('n must be non-negative');
  }
  if (!Number.isInteger(n)) {
    throw new ValueError('n must be an integer');
  }
  if (n > MAX_HYPERCUBE_DIM) {
    throw new ValueError(
      `n must be at most ${MAX_HYPERCUBE_DIM} (2^${n} points cannot be materialized)`
    );
  }
  if (n === 0) {
    return [[]]; // Single point: the empty tuple
  }

  // 2 ** n, not 1 << n: the shift operates on 32-bit signed integers, so
  // n = 31 would yield a negative count and n = 32 would wrap around to 1.
  const numPoints = 2 ** n;
  const result: number[][] = [];

  for (let i = 0; i < numPoints; i++) {
    result.push(intToBinary(i, n));
  }

  return result;
}

// ============================================
// 2. Equality polynomial (selector polynomial)
// ============================================

/**
 * Compute the equality polynomial (selector polynomial):
 *
 *   eq(b, x) = prod_{i=1}^{n} x_i^{b_i} * (1 - x_i)^{1 - b_i}
 *
 * This polynomial evaluates to 1 when x = b (as binary vectors) and 0
 * otherwise on the boolean hypercube. It is the Lagrange basis polynomial
 * for the point b in {0,1}^n.
 *
 * @param base - The binary point, either as an integer (to be bit-decomposed)
 *               or as an array of bits [b_0, b_1, ..., b_{n-1}]
 * @param variables - The polynomial variables [x_0, x_1, ..., x_{n-1}]
 * @param ring - The polynomial ring
 * @returns The equality polynomial eq(base, variables)
 *
 * @example
 * ```typescript
 * const F = GF(101n);
 * const [R, x0, x1] = MPolynomialRingConstructor(F, ['x0', 'x1']);
 *
 * // eq([0, 1], [x0, x1]) = (1 - x0) * x1
 * const eq01 = eqPolynomial([0, 1], [x0, x1], R);
 *
 * // Evaluate at (0, 1) -> 1
 * eq01.evaluate({ x0: F(0), x1: F(1) }); // 1
 *
 * // Evaluate at (1, 0) -> 0
 * eq01.evaluate({ x0: F(1), x1: F(0) }); // 0
 * ```
 */
export function eqPolynomial<C extends RingElement>(
  base: number | number[],
  variables: MPolynomial<C>[],
  ring: MPolynomialRingBase<C>
): MPolynomial<C> {
  // Bit decompose base if necessary
  const bits = typeof base === 'number' ? intToBinary(base, variables.length) : base;

  if (bits.length !== variables.length) {
    throw new ValueError(
      `base must have the same length as variables (got ${bits.length} vs ${variables.length})`
    );
  }

  let result = ring.one();
  for (let i = 0; i < bits.length; i++) {
    const xi = variables[i]!;
    const bi = bits[i]!;

    // x_i if b_i = 1, else (1 - x_i)
    const term = bi === 1 ? xi : ring.one().sub(xi);
    result = result.mul(term);
  }

  return result;
}

// ============================================
// 3. MLE interpolation
// ============================================

/**
 * Interpolate a multilinear extension (MLE) from function values.
 *
 * Given values f(0), f(1), ..., f(2^n - 1) representing a function on the
 * boolean hypercube, computes the unique multilinear polynomial p such that
 * p(b) = f(b) for all b in {0,1}^n.
 *
 * The MLE is computed as:
 *
 *   p(x) = sum_{b in {0,1}^n} f(b) * eq(b, x)
 *
 * where eq(b, x) is the equality polynomial for point b.
 *
 * @param values - Array of function values. The length determines the number
 *                 of variables: ceil(log2(values.length)). If the length is
 *                 not a power of 2, missing values are treated as 0.
 * @param field - The finite field (coefficient ring)
 * @param ring - Optional polynomial ring. If not provided, one is created
 *               with variables named x0, x1, ..., x_{n-1}.
 * @param variables - Optional polynomial variables to use
 * @returns The multilinear extension polynomial
 *
 * @example
 * ```typescript
 * const F = GF(101n);
 *
 * // Interpolate from 4 values (2 variables)
 * const values = [F(1n), F(2n), F(3n), F(4n)];
 * const poly = multilinearExtension(values, F);
 *
 * // Verify: poly(0,0) = 1, poly(0,1) = 2, poly(1,0) = 3, poly(1,1) = 4
 * ```
 */
export function multilinearExtension(
  values: FiniteFieldElement[],
  field: CoefficientRing<FiniteFieldElement>,
  ring?: MPolynomialRing<FiniteFieldElement>,
  variables?: MPolynomial<FiniteFieldElement>[]
): MPolynomial<FiniteFieldElement> {
  if (values.length === 0) {
    throw new ValueError('values array cannot be empty');
  }

  // Calculate number of variables required
  const numVars = closestPowerOfTwo(values.length);

  // Create the polynomial ring if not provided
  if (!ring) {
    // Handle edge case: if numVars is 0, we need at least one variable
    const varNames =
      numVars > 0 ? Array.from({ length: numVars }, (_, i) => `x${i}`) : ['x0'];
    [ring] = MPolynomialRingConstructor(field, varNames);
    variables = ring.gens();
  }

  if (!variables) {
    variables = ring.gens();
  }

  // Edge case: constant polynomial (single value)
  if (values.length === 1) {
    return ring.__call__(values[0]!);
  }

  // Handle case where numVars is 0 but we have exactly 1 value
  if (numVars === 0) {
    return ring.__call__(values[0]!);
  }

  // Interpolate: sum_{b} f(b) * eq(b, x)
  let result = ring.zero();
  for (let idx = 0; idx < values.length; idx++) {
    const val = values[idx]!;
    if (!val.isZero()) {
      const eq = eqPolynomial(idx, variables, ring);
      // Convert scalar to constant polynomial and multiply
      const valPoly = ring.__call__(val);
      result = result.add(eq.mul(valPoly));
    }
  }

  // Sanity check (blueprint `mle.sage:59-62`): repeated or otherwise
  // degenerate `variables` produce a polynomial of degree > 1 in some
  // variable, which is not a multilinear extension.
  if (!isMultilinear(result)) {
    throw new ValueError('Polynomial is not linear');
  }

  return result;
}

// ============================================
// 4. Sparse MLE interpolation
// ============================================

/**
 * Interpolate a sparse multilinear extension.
 *
 * For polynomials where only some basis elements are non-zero, this is more
 * efficient than the dense interpolation. Supports two modes:
 *
 * 1. **Selector polynomial**: Pass an array of indices to get a polynomial
 *    that is 1 at those indices and 0 elsewhere.
 *
 * 2. **Sparse values**: Pass a Map from indices to values.
 *
 * @param entries - Either:
 *   - An array of indices [i1, i2, ...] for selector polynomial (value = 1 at each index)
 *   - A Map<number, C> from indices to values
 * @param variables - The polynomial variables [x_0, x_1, ..., x_{n-1}]
 * @param ring - The polynomial ring
 * @returns The sparse multilinear extension polynomial
 *
 * @remarks
 * The blueprint (`reference/sage_blueprints/mle.sage:108-112`) special-cases a
 * one-element non-list `entries` and returns the *constant* polynomial
 * `R(entries[0])`, so `interpolate_sparse([5], vars)` is the constant 5 rather
 * than the selector for index 5. That branch is unreachable for a genuine
 * selector of a single index, is inconsistent with the multi-index branch
 * right below it, and is left over from debugging (it still contains `print`
 * statements). We treat a one-element array like any other selector list, so
 * `sparseMultilinearExtension([5], vars, R)` is `eq(5, x)`.
 *
 * @see Deviation: Sparse MLE single-index selector
 *
 * @example
 * ```typescript
 * const F = GF(101n);
 * const [R, x0, x1, x2] = MPolynomialRingConstructor(F, ['x0', 'x1', 'x2']);
 * const vars = [x0, x1, x2];
 *
 * // Selector polynomial: 1 at indices 0, 2, 5; 0 elsewhere
 * const selector = sparseMultilinearExtension([0, 2, 5], vars, R);
 *
 * // Sparse interpolation with values
 * const sparseVals = new Map([
 *   [0, F(9n)],
 *   [2, F(5n)],
 *   [5, F(2n)],
 * ]);
 * const sparse = sparseMultilinearExtension(sparseVals, vars, R);
 * ```
 */
export function sparseMultilinearExtension<C extends RingElement>(
  entries: number[] | Map<number, C>,
  variables: MPolynomial<C>[],
  ring: MPolynomialRingBase<C>
): MPolynomial<C> {
  if (variables.length === 0) {
    throw new ValueError('variables array cannot be empty');
  }

  // Handle selector polynomial (array of indices)
  if (Array.isArray(entries)) {
    if (entries.length === 0) {
      return ring.zero();
    }

    let result = ring.zero();
    for (const entry of entries) {
      result = result.add(eqPolynomial(entry, variables, ring));
    }
    return checkMultilinear(result);
  }

  // Handle sparse interpolation (Map of index -> value)
  if (entries.size === 0) {
    return ring.zero();
  }

  let result = ring.zero();
  for (const [entry, val] of entries) {
    const eq = eqPolynomial(entry, variables, ring);
    const valPoly = ring.__call__(val);
    result = result.add(eq.mul(valPoly));
  }
  return checkMultilinear(result);
}

/**
 * Blueprint sanity check (`mle.sage:127-129`): the interpolant must be
 * multilinear.
 */
function checkMultilinear<C extends RingElement>(poly: MPolynomial<C>): MPolynomial<C> {
  if (!isMultilinear(poly)) {
    throw new ValueError('Polynomial is not linear');
  }
  return poly;
}

// ============================================
// 5. Verification utilities
// ============================================

/**
 * Check if a polynomial is multilinear.
 *
 * A polynomial is multilinear if its degree in each variable is at most 1.
 *
 * @param poly - The polynomial to check
 * @returns true if the polynomial is multilinear, false otherwise
 *
 * @example
 * ```typescript
 * const F = GF(101n);
 * const [R, x, y] = MPolynomialRingConstructor(F, ['x', 'y']);
 *
 * isMultilinear(x.mul(y).add(x).add(y));  // true: x*y + x + y
 * isMultilinear(x.pow(2));                 // false: x^2 has degree 2 in x
 * ```
 */
export function isMultilinear<C extends RingElement>(poly: MPolynomial<C>): boolean {
  const degrees = poly.degrees();
  return degrees.every((d) => d <= 1);
}

/**
 * Evaluate an MLE at a point and verify it matches the expected value.
 *
 * @param poly - The MLE polynomial
 * @param point - The evaluation point as a record mapping variable names to field elements
 * @param expected - The expected value
 * @returns true if poly(point) equals expected, false otherwise
 *
 * @example
 * ```typescript
 * const F = GF(101n);
 * const values = [F(1n), F(2n), F(3n), F(4n)];
 * const poly = multilinearExtension(values, F);
 *
 * // Verify evaluation at (0, 0) equals 1
 * verifyMLEEvaluation(poly, { x0: F(0n), x1: F(0n) }, F(1n)); // true
 * ```
 */
export function verifyMLEEvaluation<C extends RingElement>(
  poly: MPolynomial<C>,
  point: Record<string, C>,
  expected: C
): boolean {
  const result = poly.evaluate(point) as C;
  return result.eq(expected);
}
