/**
 * @module sage/rings/polynomial/convolution
 * @description Fast Fourier Transform (FFT) and Number Theoretic Transform (NTT)
 * for polynomial multiplication and evaluation.
 *
 * Port of: sage/rings/polynomial/convolution.py
 *
 * This module provides O(n log n) polynomial multiplication using FFT-based
 * algorithms. For finite fields GF(p), NTT is used which operates entirely
 * within the field.
 *
 * Key functions:
 * - fft/ifft: Classic FFT using roots of unity
 * - ntt/intt: Number Theoretic Transform for finite fields
 * - fft_multiply: Fast polynomial multiplication
 * - evaluate_on_domain: Multi-point evaluation on FFT domain
 * - interpolate_from_domain: Interpolation from FFT evaluations
 *
 * @example
 * ```typescript
 * // NTT-based polynomial multiplication over GF(p)
 * const F = GF(17n);  // 17 = 1 + 16 = 1 + 2^4, so has 16th roots of unity
 * const [R, x] = PolynomialRingConstructor(F, 'x');
 *
 * const f = x.pow(2).add(x).add(R.one());     // x^2 + x + 1
 * const g = x.add(R.__call__(2));              // x + 2
 *
 * const product = ntt_multiply(f, g, F);
 * // Result: x^3 + 3x^2 + 3x + 2
 * ```
 */

import { gcd, inverse_mod, power_mod } from '../../arith/misc.js';
import { NotImplementedError, ValueError, ZeroDivisionError } from '../../errors.js';
import type { FiniteFieldElement, FiniteFieldPrime } from '../finite_rings/finite_field_prime.js';
import type { CoefficientRing, PolynomialRingBase, RingElement } from './polynomial_element.js';
import { Polynomial } from './polynomial_element.js';

// ============================================================================
// Helper Functions
// ============================================================================

/**
 * Check if n is a power of 2.
 */
function isPowerOfTwo(n: number): boolean {
  return n > 0 && (n & (n - 1)) === 0;
}

/**
 * Return the next power of 2 >= n.
 */
function nextPowerOfTwo(n: number): number {
  if (n <= 1) return 1;
  let p = 1;
  while (p < n) {
    p <<= 1;
  }
  return p;
}

/**
 * Compute the number of trailing zeros in the binary representation of n.
 * This gives the largest k such that 2^k divides n.
 */
function trailingZeros(n: bigint): number {
  if (n === 0n) return 0;
  let count = 0;
  while ((n & 1n) === 0n) {
    n >>= 1n;
    count++;
  }
  return count;
}

/**
 * Compute bit-reversal of index i with respect to log2(n) bits.
 */
function bitReverse(i: number, log2n: number): number {
  let result = 0;
  for (let j = 0; j < log2n; j++) {
    result = (result << 1) | ((i >> j) & 1);
  }
  return result;
}

/**
 * Perform bit-reversal permutation of array in-place.
 */
function bitReversalPermutation<T>(arr: T[]): void {
  const n = arr.length;
  const log2n = Math.log2(n);

  for (let i = 0; i < n; i++) {
    const j = bitReverse(i, log2n);
    if (i < j) {
      const temp = arr[i]!;
      arr[i] = arr[j]!;
      arr[j] = temp;
    }
  }
}

// ============================================================================
// Root of Unity Computation
// ============================================================================

/**
 * Find a primitive n-th root of unity in GF(p).
 *
 * A primitive n-th root of unity omega satisfies:
 * - omega^n = 1
 * - omega^k != 1 for 0 < k < n
 *
 * Such an omega exists in GF(p) if and only if n divides p-1.
 *
 * @param n - The order of the root of unity (must divide p-1)
 * @param field - The prime field GF(p)
 * @returns A primitive n-th root of unity
 * @throws {ValueError} If n does not divide p-1 or no primitive root exists
 *
 * @example
 * ```typescript
 * const F17 = GF(17n);
 * const omega = find_primitive_root(8, F17);  // 8 divides 16 = 17 - 1
 * // omega^8 = 1, but omega^4 != 1
 * ```
 */
export function find_primitive_root(n: number, field: FiniteFieldPrime): FiniteFieldElement {
  const p = field.characteristic;
  const pMinus1 = p - 1n;
  const bigN = BigInt(n);

  // Check that n divides p - 1
  if (pMinus1 % bigN !== 0n) {
    throw new ValueError(
      `No primitive ${n}-th root of unity exists in GF(${p}): ${n} does not divide ${pMinus1}`
    );
  }

  // Find a generator of the multiplicative group
  const g = field.multiplicative_generator();

  // omega = g^((p-1)/n) is a primitive n-th root of unity
  const exp = pMinus1 / bigN;
  const omega = g.pow(exp);

  return omega;
}

/**
 * Find the maximum power of 2 that divides p-1 for a prime p.
 * This determines the largest FFT size supported by GF(p).
 *
 * @param field - The prime field GF(p)
 * @returns The largest k such that 2^k divides p-1
 */
export function max_fft_size(field: FiniteFieldPrime): number {
  const pMinus1 = field.characteristic - 1n;
  return trailingZeros(pMinus1);
}

// ============================================================================
// Number Theoretic Transform (NTT)
// ============================================================================

/**
 * Compute the Number Theoretic Transform of a sequence.
 *
 * The NTT is the FFT over a finite field. Given a polynomial represented
 * by its coefficients [a_0, a_1, ..., a_{n-1}], compute its evaluations
 * at omega^0, omega^1, ..., omega^{n-1} where omega is a primitive n-th
 * root of unity.
 *
 * The result is: y_k = sum_{j=0}^{n-1} a_j * omega^{jk} for k = 0, ..., n-1
 *
 * @param coeffs - Array of field elements (length must be power of 2)
 * @param omega - A primitive n-th root of unity in the field
 * @param field - The prime field GF(p)
 * @returns The NTT of the input sequence
 *
 * @example
 * ```typescript
 * const F17 = GF(17n);
 * const omega = find_primitive_root(4, F17);
 * const coeffs = [F17(1n), F17(2n), F17(3n), F17(4n)];
 * const result = ntt(coeffs, omega, F17);
 * // result[k] = sum of coeffs[j] * omega^(j*k)
 * ```
 */
export function ntt(
  coeffs: FiniteFieldElement[],
  omega: FiniteFieldElement,
  field: FiniteFieldPrime
): FiniteFieldElement[] {
  const n = coeffs.length;

  if (n === 0) {
    return [];
  }

  if (n === 1) {
    return [coeffs[0]!];
  }

  if (!isPowerOfTwo(n)) {
    throw new ValueError(`NTT length must be a power of 2, got ${n}`);
  }

  // Verify omega is an n-th root of unity
  const omegaN = omega.pow(n);
  if (!omegaN.isOne()) {
    throw new ValueError('omega must be an n-th root of unity');
  }

  // Make a working copy
  const result = coeffs.map((c) => field.__call__(c.value));

  // Cooley-Tukey iterative FFT algorithm (decimation in time)
  const log2n = Math.log2(n);

  // Bit-reversal permutation
  bitReversalPermutation(result);

  // Butterfly operations
  for (let s = 1; s <= log2n; s++) {
    const m = 1 << s; // m = 2^s
    const halfM = m >> 1;

    // omega_m = omega^(n/m) is a primitive m-th root of unity
    const omegaM = omega.pow(n / m);

    for (let k = 0; k < n; k += m) {
      let omegaPower = field.one();

      for (let j = 0; j < halfM; j++) {
        const t = omegaPower.mul(result[k + j + halfM]!);
        const u = result[k + j]!;

        result[k + j] = u.add(t);
        result[k + j + halfM] = u.sub(t);

        omegaPower = omegaPower.mul(omegaM);
      }
    }
  }

  return result;
}

/**
 * Compute the Inverse Number Theoretic Transform.
 *
 * Given evaluations [y_0, y_1, ..., y_{n-1}] at omega^0, ..., omega^{n-1},
 * recover the polynomial coefficients [a_0, a_1, ..., a_{n-1}].
 *
 * The inverse NTT is: a_j = (1/n) * sum_{k=0}^{n-1} y_k * omega^{-jk}
 *
 * @param values - Array of field elements (length must be power of 2)
 * @param omega - The same primitive n-th root of unity used in forward NTT
 * @param field - The prime field GF(p)
 * @returns The inverse NTT (original polynomial coefficients)
 *
 * @example
 * ```typescript
 * const F17 = GF(17n);
 * const omega = find_primitive_root(4, F17);
 * const coeffs = [F17(1n), F17(2n), F17(3n), F17(4n)];
 *
 * const transformed = ntt(coeffs, omega, F17);
 * const recovered = intt(transformed, omega, F17);
 * // recovered equals coeffs
 * ```
 */
export function intt(
  values: FiniteFieldElement[],
  omega: FiniteFieldElement,
  field: FiniteFieldPrime
): FiniteFieldElement[] {
  const n = values.length;

  if (n === 0) {
    return [];
  }

  if (n === 1) {
    return [values[0]!];
  }

  // Compute INTT using NTT with omega^(-1)
  const omegaInv = omega.inv();
  const result = ntt(values, omegaInv, field);

  // Normalize by dividing by n
  const nInv = field.__call__(inverse_mod(BigInt(n), field.characteristic));

  return result.map((v) => v.mul(nInv));
}

// ============================================================================
// Generic FFT (for complex-like fields)
// ============================================================================

/**
 * Interface for elements that support the basic operations needed for FFT.
 * This allows FFT to work with any field-like structure.
 */
export interface FFTElement<T> {
  add(other: T): T;
  sub(other: T): T;
  mul(other: T): T;
  inv(): T;
  isZero(): boolean;
}

/**
 * Compute the Forward FFT of a sequence using a root of unity.
 *
 * This is a generic FFT that works with any field elements that support
 * the required arithmetic operations.
 *
 * @param poly - Polynomial coefficients (length must be power of 2)
 * @param omega - A primitive n-th root of unity
 * @returns FFT of the input (evaluations at powers of omega)
 */
export function fft<T extends FFTElement<T>>(poly: T[], omega: T): T[] {
  const n = poly.length;

  if (n === 0) {
    return [];
  }

  if (n === 1) {
    return [poly[0]!];
  }

  if (!isPowerOfTwo(n)) {
    throw new ValueError(`FFT length must be a power of 2, got ${n}`);
  }

  // Make a working copy
  const result: T[] = [...poly];

  // Cooley-Tukey iterative FFT
  const log2n = Math.log2(n);

  // Bit-reversal permutation
  bitReversalPermutation(result);

  // Butterfly operations
  for (let s = 1; s <= log2n; s++) {
    const m = 1 << s;
    const halfM = m >> 1;

    // Compute omega_m = omega^(n/m)
    let omegaM = omega;
    for (let i = 0; i < log2n - s; i++) {
      omegaM = omegaM.mul(omegaM);
    }

    for (let k = 0; k < n; k += m) {
      const w = result[0]!; // Start with identity (we'll compute properly)
      // Actually need a "one" element - compute omega^0 as omega/omega
      // This is a limitation of the generic interface
      // For now, track power manually

      let omegaPower: T | null = null;

      for (let j = 0; j < halfM; j++) {
        if (omegaPower === null) {
          // First iteration: omega^0 = 1
          // We need to extract "one" from the elements
          // A workaround: omega * omega^(-1) = 1
          const one = omega.mul(omega.inv());
          omegaPower = one;
        }

        const t = omegaPower.mul(result[k + j + halfM]!);
        const u = result[k + j]!;

        result[k + j] = u.add(t);
        result[k + j + halfM] = u.sub(t);

        omegaPower = omegaPower.mul(omegaM);
      }
    }
  }

  return result;
}

/**
 * Compute the Inverse FFT.
 *
 * @param values - FFT values (length must be power of 2)
 * @param omega - The same primitive n-th root of unity used in forward FFT
 * @returns Inverse FFT (original polynomial coefficients, scaled by n)
 *
 * Note: The result is NOT normalized. You need to divide each element by n
 * to get the original coefficients.
 */
export function ifft<T extends FFTElement<T>>(values: T[], omega: T): T[] {
  const n = values.length;

  if (n === 0) {
    return [];
  }

  if (n === 1) {
    return [values[0]!];
  }

  // IFFT = FFT with omega^(-1), then divide by n
  const omegaInv = omega.inv();
  return fft(values, omegaInv);
}

// ============================================================================
// FFT-based Polynomial Multiplication
// ============================================================================

/**
 * Multiply two polynomials using NTT over a finite field.
 *
 * This achieves O(n log n) multiplication complexity for polynomials
 * over GF(p) where p has suitable roots of unity.
 *
 * Algorithm:
 * 1. Pad both polynomials to length n (next power of 2 >= deg(f) + deg(g) + 1)
 * 2. Compute NTT of both polynomials
 * 3. Multiply pointwise in frequency domain
 * 4. Compute INTT to get product coefficients
 *
 * @param f - First polynomial
 * @param g - Second polynomial
 * @param field - The prime field (must have suitable roots of unity)
 * @returns Product polynomial f * g
 *
 * @example
 * ```typescript
 * const F17 = GF(17n);
 * const [R, x] = PolynomialRingConstructor(F17, 'x');
 *
 * const f = x.pow(2).add(x);           // x^2 + x
 * const g = x.add(R.__call__(F17(1))); // x + 1
 *
 * const product = ntt_multiply(f, g, F17);
 * // Result: x^3 + 2x^2 + x
 * ```
 */
export function ntt_multiply(
  f: Polynomial<FiniteFieldElement>,
  g: Polynomial<FiniteFieldElement>,
  field: FiniteFieldPrime
): Polynomial<FiniteFieldElement> {
  // Handle zero polynomials
  if (f.isZero() || g.isZero()) {
    return f.parent.zero();
  }

  // Compute required size
  const resultDegree = f.degree() + g.degree();
  const n = nextPowerOfTwo(resultDegree + 1);

  // Check that the field supports this FFT size.  ``1 << k`` overflows int32
  // for k >= 31 (Goldilocks has k = 32), so compare in bigint.
  const maxSize = 1n << BigInt(max_fft_size(field));
  if (BigInt(n) > maxSize) {
    throw new ValueError(
      `Polynomial product requires NTT size ${n}, but GF(${field.characteristic}) only supports up to ${maxSize}`
    );
  }

  // Find primitive n-th root of unity
  const omega = find_primitive_root(n, field);

  // Pad polynomials to length n
  const fCoeffs: FiniteFieldElement[] = [];
  const gCoeffs: FiniteFieldElement[] = [];

  for (let i = 0; i < n; i++) {
    fCoeffs.push(i <= f.degree() ? f.getCoeff(i) : field.zero());
    gCoeffs.push(i <= g.degree() ? g.getCoeff(i) : field.zero());
  }

  // Forward NTT
  const fNTT = ntt(fCoeffs, omega, field);
  const gNTT = ntt(gCoeffs, omega, field);

  // Pointwise multiplication
  const productNTT: FiniteFieldElement[] = [];
  for (let i = 0; i < n; i++) {
    productNTT.push(fNTT[i]!.mul(gNTT[i]!));
  }

  // Inverse NTT
  const productCoeffs = intt(productNTT, omega, field);

  // Create result polynomial (truncate to actual degree)
  const trimmedCoeffs = productCoeffs.slice(0, resultDegree + 1);

  return new Polynomial(trimmedCoeffs, f.parent);
}

/**
 * Fast polynomial multiplication using FFT/NTT.
 *
 * This is a convenience wrapper that automatically selects the appropriate
 * algorithm based on the coefficient ring.
 *
 * @param f - First polynomial
 * @param g - Second polynomial
 * @returns Product polynomial f * g
 */
export function fft_multiply<C extends RingElement>(
  f: Polynomial<C>,
  g: Polynomial<C>
): Polynomial<C> {
  // Check if we're working over a finite field
  const baseRing = f.parent.base_ring;

  // Try to detect if this is a finite field with suitable structure
  if ('characteristic' in baseRing && 'multiplicative_generator' in baseRing) {
    const field = baseRing as unknown as FiniteFieldPrime;
    const fCast = f as unknown as Polynomial<FiniteFieldElement>;
    const gCast = g as unknown as Polynomial<FiniteFieldElement>;

    try {
      const result = ntt_multiply(fCast, gCast, field);
      return result as unknown as Polynomial<C>;
    } catch {
      // Fall back to schoolbook multiplication if NTT fails
      // (e.g., field doesn't have suitable roots of unity)
    }
  }

  // Fall back to standard multiplication
  return f.mul(g);
}

// ============================================================================
// Multi-point Evaluation and Interpolation
// ============================================================================

/**
 * Evaluate a polynomial on an FFT domain.
 *
 * Given a polynomial f and a primitive n-th root of unity omega,
 * compute [f(1), f(omega), f(omega^2), ..., f(omega^{n-1})].
 *
 * This is exactly what the NTT computes when the polynomial
 * has degree < n.
 *
 * @param poly - The polynomial to evaluate
 * @param omega - A primitive n-th root of unity
 * @param field - The prime field
 * @param domainSize - Size of the domain (must be power of 2 and >= degree + 1)
 * @returns Array of evaluations at powers of omega
 *
 * @example
 * ```typescript
 * const F17 = GF(17n);
 * const [R, x] = PolynomialRingConstructor(F17, 'x');
 * const f = x.pow(2).add(x).add(R.one()); // x^2 + x + 1
 *
 * const omega = find_primitive_root(4, F17);
 * const evals = evaluate_on_domain(f, omega, F17, 4);
 * // evals[i] = f(omega^i) for i = 0, 1, 2, 3
 * ```
 */
export function evaluate_on_domain(
  poly: Polynomial<FiniteFieldElement>,
  omega: FiniteFieldElement,
  field: FiniteFieldPrime,
  domainSize?: number
): FiniteFieldElement[] {
  const n = domainSize ?? nextPowerOfTwo(poly.degree() + 1);

  if (!isPowerOfTwo(n)) {
    throw new ValueError(`Domain size must be a power of 2, got ${n}`);
  }

  if (n < poly.degree() + 1) {
    throw new ValueError(`Domain size ${n} is too small for polynomial of degree ${poly.degree()}`);
  }

  // Pad polynomial coefficients to length n
  const coeffs: FiniteFieldElement[] = [];
  for (let i = 0; i < n; i++) {
    coeffs.push(i <= poly.degree() ? poly.getCoeff(i) : field.zero());
  }

  // NTT gives evaluations at powers of omega
  return ntt(coeffs, omega, field);
}

/**
 * Interpolate a polynomial from its evaluations on an FFT domain.
 *
 * Given evaluations [y_0, y_1, ..., y_{n-1}] where y_i = f(omega^i),
 * recover the unique polynomial f of degree < n.
 *
 * This is exactly what the inverse NTT computes.
 *
 * @param values - Evaluations at powers of omega
 * @param omega - The primitive n-th root of unity used for evaluation
 * @param field - The prime field
 * @param parent - The polynomial ring
 * @returns The interpolated polynomial
 *
 * @example
 * ```typescript
 * const F17 = GF(17n);
 * const [R, x] = PolynomialRingConstructor(F17, 'x');
 * const f = x.pow(2).add(x).add(R.one());
 *
 * const omega = find_primitive_root(4, F17);
 * const evals = evaluate_on_domain(f, omega, F17, 4);
 * const recovered = interpolate_from_domain(evals, omega, F17, R);
 * // recovered equals f
 * ```
 */
export function interpolate_from_domain(
  values: FiniteFieldElement[],
  omega: FiniteFieldElement,
  field: FiniteFieldPrime,
  parent: PolynomialRingBase<FiniteFieldElement>
): Polynomial<FiniteFieldElement> {
  const n = values.length;

  if (n === 0) {
    return parent.zero();
  }

  if (!isPowerOfTwo(n)) {
    throw new ValueError(`Number of values must be a power of 2, got ${n}`);
  }

  // INTT recovers the polynomial coefficients
  const coeffs = intt(values, omega, field);

  return new Polynomial(coeffs, parent);
}

// ============================================================================
// Naive DFT for Comparison (Testing)
// ============================================================================

/**
 * Compute the DFT naively (O(n^2)) for testing purposes.
 *
 * This implements the same transform as NTT but with O(n^2) complexity.
 * Used to verify correctness of the FFT implementation.
 *
 * @param coeffs - Polynomial coefficients
 * @param omega - Primitive n-th root of unity
 * @param field - The prime field
 * @returns DFT values
 */
export function dft_naive(
  coeffs: FiniteFieldElement[],
  omega: FiniteFieldElement,
  field: FiniteFieldPrime
): FiniteFieldElement[] {
  const n = coeffs.length;
  const result: FiniteFieldElement[] = [];

  for (let k = 0; k < n; k++) {
    // y_k = sum_{j=0}^{n-1} a_j * omega^{jk}
    let sum = field.zero();
    const omegaK = omega.pow(k); // omega^k

    let omegaJK = field.one(); // omega^{0*k} = 1

    for (let j = 0; j < n; j++) {
      sum = sum.add(coeffs[j]!.mul(omegaJK));
      omegaJK = omegaJK.mul(omegaK); // omega^{(j+1)*k}
    }

    result.push(sum);
  }

  return result;
}

/**
 * Compute the inverse DFT naively for testing purposes.
 *
 * @param values - DFT values
 * @param omega - Primitive n-th root of unity
 * @param field - The prime field
 * @returns Original coefficients
 */
export function idft_naive(
  values: FiniteFieldElement[],
  omega: FiniteFieldElement,
  field: FiniteFieldPrime
): FiniteFieldElement[] {
  const n = values.length;
  const omegaInv = omega.inv();

  // Compute DFT with omega^{-1}
  const unnormalized = dft_naive(values, omegaInv, field);

  // Normalize by 1/n
  const nInv = field.__call__(inverse_mod(BigInt(n), field.characteristic));

  return unnormalized.map((v) => v.mul(nInv));
}

/**
 * Naive polynomial multiplication using convolution formula.
 *
 * This has O(n^2) complexity and is used for testing the FFT-based version.
 *
 * @param fCoeffs - Coefficients of first polynomial
 * @param gCoeffs - Coefficients of second polynomial
 * @param field - The prime field
 * @returns Coefficients of the product
 */
export function convolve_naive(
  fCoeffs: FiniteFieldElement[],
  gCoeffs: FiniteFieldElement[],
  field: FiniteFieldPrime
): FiniteFieldElement[] {
  const m = fCoeffs.length;
  const n = gCoeffs.length;

  if (m === 0 || n === 0) {
    return [];
  }

  const resultLen = m + n - 1;
  const result: FiniteFieldElement[] = [];

  for (let k = 0; k < resultLen; k++) {
    let sum = field.zero();

    const jStart = Math.max(0, k - n + 1);
    const jEnd = Math.min(k + 1, m);

    for (let j = jStart; j < jEnd; j++) {
      sum = sum.add(fCoeffs[j]!.mul(gCoeffs[k - j]!));
    }

    result.push(sum);
  }

  return result;
}

// ============================================================================
// Utility Functions for Working with FFT-friendly Primes
// ============================================================================

/**
 * Check if a prime p is suitable for NTT of a given size.
 *
 * A prime p supports NTT of size n if n divides p-1.
 *
 * @param p - The prime
 * @param n - The desired NTT size (must be power of 2)
 * @returns true if p supports NTT of size n
 */
export function supports_ntt_size(p: bigint, n: number): boolean {
  if (!isPowerOfTwo(n)) {
    return false;
  }
  return (p - 1n) % BigInt(n) === 0n;
}

/**
 * Find a prime suitable for NTT of a given size.
 *
 * Searches for a prime of the form p = k * n + 1 where p is prime
 * and has a certain number of bits.
 *
 * @param n - Required NTT size (power of 2)
 * @param minBits - Minimum number of bits for the prime (default: 32)
 * @returns A suitable prime p such that n divides p-1
 */
export function find_ntt_prime(n: number, minBits: number = 32): bigint {
  if (!isPowerOfTwo(n)) {
    throw new ValueError(`NTT size must be a power of 2, got ${n}`);
  }

  const bigN = BigInt(n);
  const minValue = 1n << BigInt(minBits - 1);

  // Start searching from k such that k*n + 1 >= 2^(minBits-1)
  let k = (minValue - 1n) / bigN;
  if (k * bigN + 1n < minValue) {
    k++;
  }

  // Simple primality test using trial division for small primes
  // and Miller-Rabin for larger ones
  const isPrime = (p: bigint): boolean => {
    if (p < 2n) return false;
    if (p === 2n) return true;
    if (p % 2n === 0n) return false;

    const smallPrimes = [3n, 5n, 7n, 11n, 13n, 17n, 19n, 23n, 29n, 31n];
    for (const sp of smallPrimes) {
      if (p === sp) return true;
      if (p % sp === 0n) return false;
    }

    // Miller-Rabin test with a few witnesses
    let d = p - 1n;
    let r = 0;
    while (d % 2n === 0n) {
      d >>= 1n;
      r++;
    }

    const witnesses = [2n, 3n, 5n, 7n, 11n];
    for (const a of witnesses) {
      if (a >= p) continue;

      let x = power_mod(a, d, p);
      if (x === 1n || x === p - 1n) continue;

      let composite = true;
      for (let i = 0; i < r - 1; i++) {
        x = power_mod(x, 2n, p);
        if (x === p - 1n) {
          composite = false;
          break;
        }
      }

      if (composite) return false;
    }

    return true;
  };

  // Search for a prime
  while (true) {
    const candidate = k * bigN + 1n;
    if (isPrime(candidate)) {
      return candidate;
    }
    k++;
  }
}

/**
 * Some well-known NTT-friendly primes.
 *
 * These primes have the form p = k * 2^n + 1 and are commonly used
 * in cryptographic applications requiring NTT.
 */
export const NTT_FRIENDLY_PRIMES = {
  // p = 3 * 2^30 + 1 = 3221225473, supports NTT up to size 2^30
  P_30: 3221225473n,

  // p = 5 * 2^25 + 1 = 167772161, supports NTT up to size 2^25
  P_25: 167772161n,

  // p = 7 * 2^26 + 1 = 469762049, supports NTT up to size 2^26
  P_26: 469762049n,

  // p = 17 (small, supports NTT up to size 16)
  P_SMALL: 17n,

  // p = 97 (small, supports NTT up to size 32)
  P_SMALL_2: 97n,

  // p = 257 (small, supports NTT up to size 256)
  P_SMALL_3: 257n,

  // BN254 scalar field prime (commonly used in zkSNARKs)
  // p = 21888242871839275222246405745257275088548364400416034343698204186575808495617
  // This is 2^254 - 2^128 * p' + 1 for some small p'
  // Supports NTT up to size 2^28
  BN254_SCALAR: 21888242871839275222246405745257275088548364400416034343698204186575808495617n,
};

// ============================================================================
// Generic convolution (port of sage/rings/polynomial/convolution.py)
// ============================================================================

/**
 * Arithmetic used by the generic convolution routines.
 *
 * Sage works over "any commutative ring in which the multiply-by-two map is
 * injective"; the only non-ring operation it needs is `R(x / 2^k)`.
 */
interface ConvolutionOps<T> {
  zero: T;
  add(a: T, b: T): T;
  sub(a: T, b: T): T;
  mul(a: T, b: T): T;
  /** Exact division by 2^k. */
  divByPowerOfTwo(a: T, k: number): T;
}

/**
 * Derive the arithmetic operations from a sample element.
 */
function convolutionOps<T>(sample: T): ConvolutionOps<T> {
  if (typeof sample === 'bigint') {
    return {
      zero: 0n as unknown as T,
      add: (a, b) => ((a as bigint) + (b as bigint)) as unknown as T,
      sub: (a, b) => ((a as bigint) - (b as bigint)) as unknown as T,
      mul: (a, b) => ((a as bigint) * (b as bigint)) as unknown as T,
      divByPowerOfTwo: (a, k) => ((a as bigint) / (1n << BigInt(k))) as unknown as T,
    };
  }

  const elem = sample as unknown as {
    add(o: unknown): unknown;
    sub(o: unknown): unknown;
    mul(o: unknown): unknown;
    div?(o: unknown): unknown;
    parent?: { zero(): unknown; one(): unknown; __call__(x: unknown): unknown };
  };

  const parent = elem.parent;
  if (
    parent === undefined ||
    typeof parent.zero !== 'function' ||
    typeof parent.__call__ !== 'function'
  ) {
    throw new NotImplementedError(
      'convolution: elements must be bigints or ring elements exposing a parent'
    );
  }

  const zero = parent.zero() as T;
  const divide = (a: T, d: unknown): T => {
    const x = a as unknown as { div?(o: unknown): unknown };
    if (typeof x.div === 'function') {
      return x.div(d) as T;
    }
    const inv = (d as { inv?(): unknown }).inv;
    if (typeof inv === 'function') {
      return (a as unknown as { mul(o: unknown): unknown }).mul(inv.call(d)) as T;
    }
    throw new NotImplementedError('convolution: base ring does not support division by 2^k');
  };

  return {
    zero,
    add: (a, b) => (a as unknown as { add(o: unknown): unknown }).add(b) as T,
    sub: (a, b) => (a as unknown as { sub(o: unknown): unknown }).sub(b) as T,
    mul: (a, b) => (a as unknown as { mul(o: unknown): unknown }).mul(b) as T,
    divByPowerOfTwo: (a, k) => divide(a, parent.__call__(1n << BigInt(k))),
  };
}

/**
 * Return the convolution of two non-empty lists.
 *
 * ```
 * convolution([1, 2, 3, 4, 5], [6, 7]) == [6, 19, 32, 45, 58, 35]
 * ```
 *
 * @see Reference: sage/rings/polynomial/convolution.py:57 (convolution)
 */
export function convolution<T>(L1: T[], L2: T[]): T[] {
  if (L1.length === 0 || L2.length === 0) {
    throw new ValueError('cannot compute convolution of empty lists');
  }
  if (L1.length <= 100 && L2.length <= 100) {
    // very arbitrary cutoff (as in Sage)
    return _convolution_naive(L1, L2);
  }
  return _convolution_fft(L1, L2);
}

/**
 * Convolution of two non-empty lists, using the naive algorithm.
 *
 * @see Reference: sage/rings/polynomial/convolution.py:89 (_convolution_naive)
 */
export function _convolution_naive<T>(L1: T[], L2: T[]): T[] {
  const ops = convolutionOps<T>(L1[0] as T);
  const m1 = L1.length;
  const m2 = L2.length;
  const out: T[] = [];
  for (let k = 0; k < m1 + m2 - 1; k++) {
    let acc = ops.zero;
    for (let i = Math.max(0, k - m2 + 1); i < Math.min(k + 1, m1); i++) {
      acc = ops.add(acc, ops.mul(L1[i]!, L2[k - i]!));
    }
    out.push(acc);
  }
  return out;
}

/**
 * Negacyclic convolution of two equal-length lists, naive algorithm.
 *
 * @see Reference: sage/rings/polynomial/convolution.py:120 (_negaconvolution_naive)
 */
export function _negaconvolution_naive<T>(L1: T[], L2: T[]): T[] {
  if (L1.length === 0 || L1.length !== L2.length) {
    throw new ValueError('_negaconvolution_naive requires two non-empty lists of equal length');
  }
  const ops = convolutionOps<T>(L1[0] as T);
  const N = L1.length;
  const out: T[] = [];
  for (let j = 0; j < N; j++) {
    let acc = ops.zero;
    for (let i = 0; i <= j; i++) {
      acc = ops.add(acc, ops.mul(L1[i]!, L2[j - i]!));
    }
    for (let i = j + 1; i < N; i++) {
      acc = ops.sub(acc, ops.mul(L1[i]!, L2[N + j - i]!));
    }
    out.push(acc);
  }
  return out;
}

/**
 * `(L1 + y^r L2, L1 - y^r L2)` in `S = R[y]/(y^K + 1)`.
 *
 * @see Reference: sage/rings/polynomial/convolution.py:150 (_forward_butterfly)
 */
function _forward_butterfly<T>(L1: T[], L2: T[], r: number, ops: ConvolutionOps<T>): [T[], T[]] {
  const K = L1.length;
  const a: T[] = [];
  const b: T[] = [];
  for (let i = 0; i < r; i++) {
    a.push(ops.sub(L1[i]!, L2[i + K - r]!));
    b.push(ops.add(L1[i]!, L2[i + K - r]!));
  }
  for (let i = r; i < K; i++) {
    a.push(ops.add(L1[i]!, L2[i - r]!));
    b.push(ops.sub(L1[i]!, L2[i - r]!));
  }
  return [a, b];
}

/**
 * `(L1 + L2, y^{-r} (L1 - L2))` in `S = R[y]/(y^K + 1)`.
 *
 * @see Reference: sage/rings/polynomial/convolution.py:167 (_inverse_butterfly)
 */
function _inverse_butterfly<T>(L1: T[], L2: T[], r: number, ops: ConvolutionOps<T>): [T[], T[]] {
  const K = L1.length;
  const a: T[] = [];
  for (let i = 0; i < K; i++) {
    a.push(ops.add(L1[i]!, L2[i]!));
  }
  const b: T[] = [];
  for (let i = r; i < K; i++) {
    b.push(ops.sub(L1[i]!, L2[i]!));
  }
  for (let i = 0; i < r; i++) {
    b.push(ops.sub(L2[i]!, L1[i]!));
  }
  return [a, b];
}

/**
 * In-place FFT over `S[x]/(x^D - y^(2 root))`, results in bit-reversed order.
 *
 * @see Reference: sage/rings/polynomial/convolution.py:183 (_fft)
 */
function _fft<T>(
  L: T[][],
  K: number,
  start: number,
  depth: number,
  root: number,
  ops: ConvolutionOps<T>
): void {
  const half = 1 << (depth - 1);
  const start2 = start + half;

  for (let i = 0; i < half; i++) {
    const [a, b] = _forward_butterfly(L[start + i]!, L[start2 + i]!, root, ops);
    L[start + i] = a;
    L[start2 + i] = b;
  }

  if (depth >= 2) {
    _fft(L, K, start, depth - 1, root >> 1, ops);
    _fft(L, K, start2, depth - 1, (root + K) >> 1, ops);
  }
}

/**
 * Inverse of {@link _fft} (up to a factor of `2^depth`).
 *
 * @see Reference: sage/rings/polynomial/convolution.py:209 (_ifft)
 */
function _ifft<T>(
  L: T[][],
  K: number,
  start: number,
  depth: number,
  root: number,
  ops: ConvolutionOps<T>
): void {
  const half = 1 << (depth - 1);
  const start2 = start + half;

  if (depth >= 2) {
    _ifft(L, K, start, depth - 1, root >> 1, ops);
    _ifft(L, K, start2, depth - 1, (root + K) >> 1, ops);
  }

  for (let i = 0; i < half; i++) {
    const [a, b] = _inverse_butterfly(L[start + i]!, L[start2 + i]!, root, ops);
    L[start + i] = a;
    L[start2 + i] = b;
  }
}

/**
 * Split a list of length `2^(m+k-1)` into `2^m` zero-padded lists of length `2^k`.
 *
 * @see Reference: sage/rings/polynomial/convolution.py:228 (_split)
 */
export function _split<T>(L: T[], m: number, k: number, ops?: ConvolutionOps<T>): T[][] {
  const o = ops ?? convolutionOps<T>(L[0] as T);
  const K = 1 << (k - 1);
  const out: T[][] = [];
  for (let i = 0; i < K << m; i += K) {
    const piece: T[] = [];
    for (let j = 0; j < K; j++) {
      piece.push(L[i + j]!);
    }
    for (let j = 0; j < K; j++) {
      piece.push(o.zero);
    }
    out.push(piece);
  }
  return out;
}

/**
 * Inverse of {@link _split}, overlaying the pieces (second half of the last
 * piece is discarded).
 *
 * @see Reference: sage/rings/polynomial/convolution.py:248 (_combine)
 */
function _combine<T>(L: T[][], m: number, k: number, ops: ConvolutionOps<T>): T[] {
  const M = 1 << m;
  const halfK = 1 << (k - 1);
  const out: T[] = [];
  for (let j = 0; j < halfK; j++) {
    out.push(L[0]![j]!);
  }
  for (let i = 0; i < M - 1; i++) {
    for (let j = 0; j < halfK; j++) {
      out.push(ops.add(L[i + 1]![j]!, L[i]![j + halfK]!));
    }
  }
  return out;
}

/**
 * Like {@link _combine} but wrapping the last piece around negacyclically.
 *
 * @see Reference: sage/rings/polynomial/convolution.py:265 (_nega_combine)
 */
function _nega_combine<T>(L: T[][], m: number, k: number, ops: ConvolutionOps<T>): T[] {
  const M = 1 << m;
  const halfK = 1 << (k - 1);
  const out: T[] = [];
  for (let j = 0; j < halfK; j++) {
    out.push(ops.sub(L[0]![j]!, L[M - 1]![j + halfK]!));
  }
  for (let i = 0; i < M - 1; i++) {
    for (let j = 0; j < halfK; j++) {
      out.push(ops.add(L[i + 1]![j]!, L[i]![j + halfK]!));
    }
  }
  return out;
}

/**
 * Negacyclic convolution of two lists of length `2^n`.
 *
 * @see Reference: sage/rings/polynomial/convolution.py:285 (_negaconvolution)
 */
export function _negaconvolution<T>(L1: T[], L2: T[], n: number): T[] {
  if (n <= 3) {
    // arbitrary cutoff (as in Sage)
    return _negaconvolution_naive(L1, L2);
  }
  return _negaconvolution_fft(L1, L2, n);
}

/**
 * Negacyclic convolution of two lists of length `2^n`, `n >= 3`, via FFT.
 *
 * @see Reference: sage/rings/polynomial/convolution.py:294 (_negaconvolution_fft)
 */
export function _negaconvolution_fft<T>(L1: T[], L2: T[], n: number): T[] {
  if (n < 3) {
    throw new ValueError('_negaconvolution_fft requires n >= 3');
  }
  const ops = convolutionOps<T>(L1[0] as T);

  // split into 2^m pieces of 2^(k-1) coefficients each, with k as small as
  // possible, subject to m <= k
  const m = (n + 1) >> 1;
  const k = n + 1 - m;
  const M = 1 << m;
  const K = 1 << k;

  const P1 = _split(L1, m, k, ops);
  const P2 = _split(L2, m, k, ops);

  _fft(P1, K, 0, m, K >> 1, ops);
  _fft(P2, K, 0, m, K >> 1, ops);

  const P3: T[][] = [];
  for (let i = 0; i < M; i++) {
    P3.push(_negaconvolution(P1[i]!, P2[i]!, k));
  }

  _ifft(P3, K, 0, m, K >> 1, ops);

  const combined = _nega_combine(P3, m, k, ops);
  return combined.map((x) => ops.divByPowerOfTwo(x, m));
}

/**
 * Convolution of two non-empty lists, using the Schoenhage-style FFT.
 *
 * @see Reference: sage/rings/polynomial/convolution.py:346 (_convolution_fft)
 */
export function _convolution_fft<T>(L1: T[], L2: T[]): T[] {
  const ops = convolutionOps<T>(L1[0] as T);

  const len1 = L1.length;
  const len2 = L2.length;
  const outlen = len1 + len2 - 1;

  // choose n so that the output convolution length is at most 2^n
  let n = 0;
  while (1 << n < outlen) {
    n++;
  }

  // split into 2^m pieces of 2^(k-1) coefficients each, with k as small as
  // possible, subject to m <= k + 1
  const m = (n >> 1) + 1;
  const k = n + 1 - m;

  const N = 1 << n;
  const K = 1 << k;

  const padded1 = [...L1];
  const padded2 = [...L2];
  while (padded1.length < N) padded1.push(ops.zero);
  while (padded2.length < N) padded2.push(ops.zero);

  const P1 = _split(padded1, m, k, ops);
  const P2 = _split(padded2, m, k, ops);

  _fft(P1, K, 0, m, K, ops);
  _fft(P2, K, 0, m, K, ops);

  const M = 1 << m;
  const P3: T[][] = [];
  for (let i = 0; i < M; i++) {
    P3.push(_negaconvolution(P1[i]!, P2[i]!, k));
  }

  _ifft(P3, K, 0, m, K, ops);

  const combined = _combine(P3, m, k, ops);
  const out: T[] = [];
  for (let i = 0; i < outlen; i++) {
    out.push(ops.divByPowerOfTwo(combined[i]!, m));
  }
  return out;
}
