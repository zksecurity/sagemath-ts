/**
 * FLINT polynomials over the integers (fmpz_poly)
 *
 * Polynomials with fmpz (arbitrary precision integer) coefficients.
 * Uses dense representation with automatic memory management.
 *
 * @see Reference: flint/src/fmpz_poly.h, flint/src/fmpz_types.h
 * @module
 */

import { fmpz } from './fmpz.js';

/**
 * Result of polynomial factorization.
 *
 * @see Reference: flint/src/fmpz_types.h (fmpz_poly_factor_struct)
 */
export interface fmpz_poly_factor {
  /** Content (leading constant factor) */
  readonly content: fmpz;
  /** Irreducible polynomial factors */
  readonly factors: readonly fmpz_poly[];
  /** Corresponding exponents */
  readonly exponents: readonly bigint[];
  /** Number of factors */
  readonly num: number;
}

/**
 * Polynomial over the integers with fmpz coefficients.
 *
 * Stored as a dense array of coefficients, from degree 0 up.
 * The polynomial c_0 + c_1*x + c_2*x^2 + ... + c_n*x^n is stored
 * as [c_0, c_1, c_2, ..., c_n].
 *
 * @see Reference: flint/src/fmpz_poly.h
 */
export class fmpz_poly {
  /** Coefficients array */
  private _coeffs: fmpz[];
  /** Allocated length */
  private _alloc: number;
  /** Actual length (degree + 1) */
  private _length: number;

  /**
   * Initialize a new polynomial.
   * Corresponds to fmpz_poly_init.
   */
  constructor() {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly constructor');
  }

  // ============================================================
  // Memory management
  // ============================================================

  /**
   * Initialize polynomial to zero.
   * Corresponds to fmpz_poly_init.
   */
  static init(): fmpz_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_init');
  }

  /**
   * Initialize with pre-allocated space.
   * Corresponds to fmpz_poly_init2.
   *
   * @param alloc - Number of coefficients to allocate space for
   */
  static init2(alloc: number): fmpz_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_init2');
  }

  /**
   * Clear/free polynomial resources.
   * Corresponds to fmpz_poly_clear.
   */
  clear(): void {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_clear');
  }

  /**
   * Reallocate to hold at least alloc coefficients.
   * Corresponds to fmpz_poly_realloc.
   */
  realloc(alloc: number): void {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_realloc');
  }

  /**
   * Ensure space for at least len coefficients.
   * Corresponds to fmpz_poly_fit_length.
   */
  fitLength(len: number): void {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_fit_length');
  }

  /**
   * Create a copy of this polynomial.
   */
  clone(): fmpz_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly clone');
  }

  // ============================================================
  // Polynomial parameters
  // ============================================================

  /**
   * Get the length (degree + 1, or 0 for zero polynomial).
   * Corresponds to fmpz_poly_length.
   */
  length(): number {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_length');
  }

  /**
   * Get the degree (length - 1, or -1 for zero polynomial).
   * Corresponds to fmpz_poly_degree.
   */
  degree(): number {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_degree');
  }

  // ============================================================
  // Assignment and basic manipulation
  // ============================================================

  /**
   * Set from another polynomial.
   * Corresponds to fmpz_poly_set.
   */
  set(other: fmpz_poly): void {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_set');
  }

  /**
   * Set to zero polynomial.
   * Corresponds to fmpz_poly_zero.
   */
  zero(): void {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_zero');
  }

  /**
   * Set to one (constant 1).
   * Corresponds to fmpz_poly_one.
   */
  one(): void {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_one');
  }

  /**
   * Set from a constant fmpz.
   * Corresponds to fmpz_poly_set_fmpz.
   */
  setFmpz(c: fmpz): void {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_set_fmpz');
  }

  /**
   * Set from string representation.
   * Corresponds to fmpz_poly_set_str.
   *
   * @returns 0 on success, -1 on failure
   */
  setStr(str: string): number {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_set_str');
  }

  /**
   * Get string representation.
   * Corresponds to fmpz_poly_get_str.
   */
  toString(): string {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_get_str');
  }

  /**
   * Get pretty string representation.
   * Corresponds to fmpz_poly_get_str_pretty.
   *
   * @param varname - Variable name (default "x")
   */
  toPrettyString(varname?: string): string {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_get_str_pretty');
  }

  /**
   * Swap two polynomials.
   * Corresponds to fmpz_poly_swap.
   */
  swap(other: fmpz_poly): void {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_swap');
  }

  /**
   * Truncate to length n.
   * Corresponds to fmpz_poly_truncate.
   */
  truncate(n: number): void {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_truncate');
  }

  /**
   * Reverse coefficients.
   * Corresponds to fmpz_poly_reverse.
   *
   * @param n - Length to reverse (pads with zeros if needed)
   */
  static reverse(poly: fmpz_poly, n: number): fmpz_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_reverse');
  }

  // ============================================================
  // Getting and setting coefficients
  // ============================================================

  /**
   * Get coefficient at position n.
   * Corresponds to fmpz_poly_get_coeff_fmpz.
   */
  getCoeff(n: number): fmpz {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_get_coeff_fmpz');
  }

  /**
   * Set coefficient at position n.
   * Corresponds to fmpz_poly_set_coeff_fmpz.
   */
  setCoeff(n: number, x: fmpz): void {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_set_coeff_fmpz');
  }

  /**
   * Set coefficient at position n from bigint.
   * Corresponds to fmpz_poly_set_coeff_si / fmpz_poly_set_coeff_ui.
   */
  setCoeffInt(n: number, x: bigint | number): void {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_set_coeff_si');
  }

  /**
   * Get leading coefficient.
   * Returns null for zero polynomial.
   */
  leadCoeff(): fmpz | null {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_lead');
  }

  // ============================================================
  // Comparison
  // ============================================================

  /**
   * Check if zero polynomial.
   * Corresponds to fmpz_poly_is_zero.
   */
  isZero(): boolean {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_is_zero');
  }

  /**
   * Check if constant 1.
   * Corresponds to fmpz_poly_is_one.
   */
  isOne(): boolean {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_is_one');
  }

  /**
   * Check if unit (constant +/-1).
   * Corresponds to fmpz_poly_is_unit.
   */
  isUnit(): boolean {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_is_unit');
  }

  /**
   * Check if generator (x).
   * Corresponds to fmpz_poly_is_gen.
   */
  isGen(): boolean {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_is_gen');
  }

  /**
   * Check equality.
   * Corresponds to fmpz_poly_equal.
   */
  equals(other: fmpz_poly): boolean {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_equal');
  }

  // ============================================================
  // Addition and subtraction
  // ============================================================

  /**
   * Add two polynomials.
   * Corresponds to fmpz_poly_add.
   */
  static add(poly1: fmpz_poly, poly2: fmpz_poly): fmpz_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_add');
  }

  /**
   * Subtract two polynomials.
   * Corresponds to fmpz_poly_sub.
   */
  static sub(poly1: fmpz_poly, poly2: fmpz_poly): fmpz_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_sub');
  }

  /**
   * Negate polynomial.
   * Corresponds to fmpz_poly_neg.
   */
  static neg(poly: fmpz_poly): fmpz_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_neg');
  }

  /**
   * Add constant fmpz.
   * Corresponds to fmpz_poly_add_fmpz.
   */
  static addFmpz(poly: fmpz_poly, c: fmpz): fmpz_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_add_fmpz');
  }

  /**
   * Subtract constant fmpz.
   * Corresponds to fmpz_poly_sub_fmpz.
   */
  static subFmpz(poly: fmpz_poly, c: fmpz): fmpz_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_sub_fmpz');
  }

  // ============================================================
  // Scalar multiplication
  // ============================================================

  /**
   * Multiply by fmpz scalar.
   * Corresponds to fmpz_poly_scalar_mul_fmpz.
   */
  static scalarMul(poly: fmpz_poly, c: fmpz): fmpz_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_scalar_mul_fmpz');
  }

  /**
   * Exact division by fmpz scalar.
   * Corresponds to fmpz_poly_scalar_divexact_fmpz.
   */
  static scalarDivexact(poly: fmpz_poly, c: fmpz): fmpz_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_scalar_divexact_fmpz');
  }

  /**
   * Floor division by fmpz scalar.
   * Corresponds to fmpz_poly_scalar_fdiv_fmpz.
   */
  static scalarFdiv(poly: fmpz_poly, c: fmpz): fmpz_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_scalar_fdiv_fmpz');
  }

  /**
   * Reduce coefficients modulo fmpz.
   * Corresponds to fmpz_poly_scalar_mod_fmpz.
   */
  static scalarMod(poly: fmpz_poly, c: fmpz): fmpz_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_scalar_mod_fmpz');
  }

  // ============================================================
  // Multiplication
  // ============================================================

  /**
   * Multiply two polynomials.
   * Corresponds to fmpz_poly_mul.
   */
  static mul(poly1: fmpz_poly, poly2: fmpz_poly): fmpz_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_mul');
  }

  /**
   * Multiply and truncate to n terms.
   * Corresponds to fmpz_poly_mullow.
   */
  static mullow(poly1: fmpz_poly, poly2: fmpz_poly, n: number): fmpz_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_mullow');
  }

  /**
   * Multiply and keep only high terms.
   * Corresponds to fmpz_poly_mulhigh_n.
   */
  static mulhigh(poly1: fmpz_poly, poly2: fmpz_poly, n: number): fmpz_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_mulhigh_n');
  }

  /**
   * Square a polynomial.
   * Corresponds to fmpz_poly_sqr.
   */
  static sqr(poly: fmpz_poly): fmpz_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_sqr');
  }

  // ============================================================
  // Powering
  // ============================================================

  /**
   * Compute polynomial power.
   * Corresponds to fmpz_poly_pow.
   */
  static pow(poly: fmpz_poly, e: bigint | number): fmpz_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_pow');
  }

  /**
   * Compute power and truncate.
   * Corresponds to fmpz_poly_pow_trunc.
   */
  static powTrunc(poly: fmpz_poly, e: bigint | number, n: number): fmpz_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_pow_trunc');
  }

  // ============================================================
  // Shifting
  // ============================================================

  /**
   * Shift left (multiply by x^n).
   * Corresponds to fmpz_poly_shift_left.
   */
  static shiftLeft(poly: fmpz_poly, n: number): fmpz_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_shift_left');
  }

  /**
   * Shift right (divide by x^n, discarding remainder).
   * Corresponds to fmpz_poly_shift_right.
   */
  static shiftRight(poly: fmpz_poly, n: number): fmpz_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_shift_right');
  }

  // ============================================================
  // Division
  // ============================================================

  /**
   * Division with remainder.
   * Corresponds to fmpz_poly_divrem.
   *
   * @returns [quotient, remainder] where A = Q*B + R
   */
  static divrem(A: fmpz_poly, B: fmpz_poly): [fmpz_poly, fmpz_poly] {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_divrem');
  }

  /**
   * Division quotient.
   * Corresponds to fmpz_poly_div.
   */
  static div(A: fmpz_poly, B: fmpz_poly): fmpz_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_div');
  }

  /**
   * Division remainder.
   * Corresponds to fmpz_poly_rem.
   */
  static rem(A: fmpz_poly, B: fmpz_poly): fmpz_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_rem');
  }

  /**
   * Exact division (assumes B divides A exactly).
   * Corresponds to fmpz_poly_divexact.
   */
  static divexact(A: fmpz_poly, B: fmpz_poly): fmpz_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_divexact');
  }

  /**
   * Pseudo-division with remainder.
   * Corresponds to fmpz_poly_pseudo_divrem.
   *
   * Computes Q, R, d such that l^d * A = Q*B + R where l is lead(B).
   *
   * @returns [quotient, remainder, d]
   */
  static pseudoDivrem(A: fmpz_poly, B: fmpz_poly): [fmpz_poly, fmpz_poly, bigint] {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_pseudo_divrem');
  }

  /**
   * Division by root: computes Q such that A = Q*(x - c).
   * Corresponds to fmpz_poly_div_root_fmpz.
   */
  static divRoot(A: fmpz_poly, c: fmpz): fmpz_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_div_root_fmpz');
  }

  // ============================================================
  // GCD and related
  // ============================================================

  /**
   * Greatest common divisor.
   * Corresponds to fmpz_poly_gcd.
   */
  static gcd(poly1: fmpz_poly, poly2: fmpz_poly): fmpz_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_gcd');
  }

  /**
   * Least common multiple.
   * Corresponds to fmpz_poly_lcm.
   */
  static lcm(poly1: fmpz_poly, poly2: fmpz_poly): fmpz_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_lcm');
  }

  /**
   * Extended GCD.
   * Corresponds to fmpz_poly_xgcd.
   *
   * @returns [r, s, t] where r = gcd(poly1, poly2) = s*poly1 + t*poly2
   */
  static xgcd(poly1: fmpz_poly, poly2: fmpz_poly): [fmpz, fmpz_poly, fmpz_poly] {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_xgcd');
  }

  /**
   * Resultant of two polynomials.
   * Corresponds to fmpz_poly_resultant.
   */
  static resultant(poly1: fmpz_poly, poly2: fmpz_poly): fmpz {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_resultant');
  }

  /**
   * Discriminant of a polynomial.
   * Corresponds to fmpz_poly_discriminant.
   */
  static discriminant(poly: fmpz_poly): fmpz {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_discriminant');
  }

  // ============================================================
  // Content and primitive part
  // ============================================================

  /**
   * Content (GCD of coefficients).
   * Corresponds to fmpz_poly_content.
   */
  static content(poly: fmpz_poly): fmpz {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_content');
  }

  /**
   * Primitive part (poly / content).
   * Corresponds to fmpz_poly_primitive_part.
   */
  static primitivePart(poly: fmpz_poly): fmpz_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_primitive_part');
  }

  /**
   * Check if squarefree.
   * Corresponds to fmpz_poly_is_squarefree.
   */
  static isSquarefree(poly: fmpz_poly): boolean {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_is_squarefree');
  }

  // ============================================================
  // Derivative and integral
  // ============================================================

  /**
   * Derivative.
   * Corresponds to fmpz_poly_derivative.
   */
  static derivative(poly: fmpz_poly): fmpz_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_derivative');
  }

  /**
   * N-th derivative.
   * Corresponds to fmpz_poly_nth_derivative.
   */
  static nthDerivative(poly: fmpz_poly, n: bigint | number): fmpz_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_nth_derivative');
  }

  // ============================================================
  // Evaluation
  // ============================================================

  /**
   * Evaluate at fmpz point.
   * Corresponds to fmpz_poly_evaluate_fmpz.
   */
  static evaluateFmpz(poly: fmpz_poly, a: fmpz): fmpz {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_evaluate_fmpz');
  }

  /**
   * Evaluate at multiple points.
   * Corresponds to fmpz_poly_evaluate_fmpz_vec.
   */
  static evaluateFmpzVec(poly: fmpz_poly, points: fmpz[]): fmpz[] {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_evaluate_fmpz_vec');
  }

  // ============================================================
  // Composition
  // ============================================================

  /**
   * Compose: compute poly1(poly2(x)).
   * Corresponds to fmpz_poly_compose.
   */
  static compose(poly1: fmpz_poly, poly2: fmpz_poly): fmpz_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_compose');
  }

  /**
   * Taylor shift: compute poly(x + c).
   * Corresponds to fmpz_poly_taylor_shift.
   */
  static taylorShift(poly: fmpz_poly, c: fmpz): fmpz_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_taylor_shift');
  }

  // ============================================================
  // Power series operations
  // ============================================================

  /**
   * Power series inverse mod x^n.
   * Corresponds to fmpz_poly_inv_series.
   */
  static invSeries(poly: fmpz_poly, n: number): fmpz_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_inv_series');
  }

  /**
   * Power series division mod x^n.
   * Corresponds to fmpz_poly_div_series.
   */
  static divSeries(A: fmpz_poly, B: fmpz_poly, n: number): fmpz_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_div_series');
  }

  /**
   * Power series composition mod x^n.
   * Corresponds to fmpz_poly_compose_series.
   */
  static composeSeries(poly1: fmpz_poly, poly2: fmpz_poly, n: number): fmpz_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_compose_series');
  }

  /**
   * Compositional inverse mod x^n.
   * Corresponds to fmpz_poly_revert_series.
   */
  static revertSeries(poly: fmpz_poly, n: number): fmpz_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_revert_series');
  }

  // ============================================================
  // Square root
  // ============================================================

  /**
   * Polynomial square root (if it exists).
   * Corresponds to fmpz_poly_sqrt.
   *
   * @returns square root if exists, null otherwise
   */
  static sqrt(poly: fmpz_poly): fmpz_poly | null {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_sqrt');
  }

  /**
   * Power series square root mod x^n.
   * Corresponds to fmpz_poly_sqrt_series.
   */
  static sqrtSeries(poly: fmpz_poly, n: number): fmpz_poly | null {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_sqrt_series');
  }

  // ============================================================
  // Norms
  // ============================================================

  /**
   * 2-norm (Euclidean norm of coefficient vector).
   * Corresponds to fmpz_poly_2norm.
   */
  static norm2(poly: fmpz_poly): fmpz {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_2norm');
  }

  /**
   * Height (max absolute coefficient).
   * Corresponds to fmpz_poly_height.
   */
  static height(poly: fmpz_poly): fmpz {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_height');
  }

  /**
   * Maximum coefficient bits.
   * Corresponds to fmpz_poly_max_bits.
   */
  static maxBits(poly: fmpz_poly): number {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_max_bits');
  }

  // ============================================================
  // Interpolation
  // ============================================================

  /**
   * Interpolate polynomial through points.
   * Corresponds to fmpz_poly_interpolate_fmpz_vec.
   */
  static interpolate(xs: fmpz[], ys: fmpz[]): fmpz_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_interpolate_fmpz_vec');
  }

  // ============================================================
  // Roots
  // ============================================================

  /**
   * Bound on absolute value of roots.
   * Corresponds to fmpz_poly_bound_roots.
   */
  static boundRoots(poly: fmpz_poly): fmpz {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_bound_roots');
  }

  /**
   * Number of real roots.
   * Corresponds to fmpz_poly_num_real_roots.
   */
  static numRealRoots(poly: fmpz_poly): number {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_num_real_roots');
  }

  /**
   * Signature (number of real, complex root pairs).
   * Corresponds to fmpz_poly_signature.
   *
   * @returns [r1, r2] where r1 = real roots, r2 = complex conjugate pairs
   */
  static signature(poly: fmpz_poly): [number, number] {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_signature');
  }

  // ============================================================
  // Special polynomials
  // ============================================================

  /**
   * N-th cyclotomic polynomial.
   * Corresponds to fmpz_poly_cyclotomic.
   */
  static cyclotomic(n: bigint | number): fmpz_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_cyclotomic');
  }

  /**
   * Check if cyclotomic polynomial.
   * Corresponds to fmpz_poly_is_cyclotomic.
   *
   * @returns n if poly is Phi_n, 0 otherwise
   */
  static isCyclotomic(poly: fmpz_poly): bigint {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_is_cyclotomic');
  }

  /**
   * Chebyshev polynomial T_n.
   * Corresponds to fmpz_poly_chebyshev_t.
   */
  static chebyshevT(n: bigint | number): fmpz_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_chebyshev_t');
  }

  /**
   * Chebyshev polynomial U_n.
   * Corresponds to fmpz_poly_chebyshev_u.
   */
  static chebyshevU(n: bigint | number): fmpz_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_chebyshev_u');
  }

  /**
   * Hermite polynomial H_n.
   * Corresponds to fmpz_poly_hermite_h.
   */
  static hermiteH(n: bigint | number): fmpz_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_hermite_h');
  }

  /**
   * Fibonacci polynomial F_n.
   * Corresponds to fmpz_poly_fibonacci.
   */
  static fibonacci(n: bigint | number): fmpz_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_fibonacci');
  }

  // ============================================================
  // Products from roots
  // ============================================================

  /**
   * Product of (x - roots[i]).
   * Corresponds to fmpz_poly_product_roots_fmpz_vec.
   */
  static productRoots(roots: fmpz[]): fmpz_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_poly_product_roots_fmpz_vec');
  }
}

// ============================================================
// Functional API (FLINT-style)
// ============================================================

/**
 * Initialize polynomial.
 * @see fmpz_poly.init
 */
export function fmpz_poly_init(): fmpz_poly {
  return fmpz_poly.init();
}

/**
 * Clear polynomial.
 * @see fmpz_poly.clear
 */
export function fmpz_poly_clear(poly: fmpz_poly): void {
  poly.clear();
}

/**
 * Set coefficient at position n.
 * @see fmpz_poly.setCoeff
 */
export function fmpz_poly_set_coeff_fmpz(poly: fmpz_poly, n: number, x: fmpz): void {
  poly.setCoeff(n, x);
}

/**
 * Add two polynomials.
 * @see fmpz_poly.add
 */
export function fmpz_poly_add(poly1: fmpz_poly, poly2: fmpz_poly): fmpz_poly {
  return fmpz_poly.add(poly1, poly2);
}

/**
 * Multiply two polynomials.
 * @see fmpz_poly.mul
 */
export function fmpz_poly_mul(poly1: fmpz_poly, poly2: fmpz_poly): fmpz_poly {
  return fmpz_poly.mul(poly1, poly2);
}

/**
 * GCD of two polynomials.
 * @see fmpz_poly.gcd
 */
export function fmpz_poly_gcd(poly1: fmpz_poly, poly2: fmpz_poly): fmpz_poly {
  return fmpz_poly.gcd(poly1, poly2);
}
