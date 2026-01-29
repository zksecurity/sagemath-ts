/**
 * FLINT polynomials modulo a word-size prime (nmod_poly)
 *
 * Polynomials with coefficients in Z/nZ where n fits in a single word.
 * Optimized for fast arithmetic using precomputed inverses.
 *
 * @see Reference: flint/src/nmod_poly.h, flint/src/nmod_types.h
 * @module
 */

import { fmpz } from './fmpz.js';

/**
 * Modular arithmetic context for single-word modulus.
 *
 * Stores the modulus n along with precomputed inverse for Barrett reduction.
 *
 * @see Reference: flint/src/nmod.h (nmod_t)
 */
export interface nmod_t {
  /** The modulus */
  readonly n: bigint;
  /** Precomputed inverse for Barrett reduction */
  readonly ninv: bigint;
  /** Normalization count (leading zeros) */
  readonly norm: number;
}

/**
 * Result of polynomial factorization.
 *
 * @see Reference: flint/src/nmod_types.h (nmod_poly_factor_struct)
 */
export interface nmod_poly_factor {
  /** Irreducible polynomial factors */
  readonly factors: readonly nmod_poly[];
  /** Corresponding exponents */
  readonly exponents: readonly bigint[];
  /** Number of factors */
  readonly num: number;
}

/**
 * Polynomial over Z/nZ where n is a word-size modulus.
 *
 * Stored as a dense array of word-size coefficients.
 * The polynomial c_0 + c_1*x + ... + c_d*x^d is stored as [c_0, c_1, ..., c_d].
 *
 * @see Reference: flint/src/nmod_poly.h
 */
export class nmod_poly {
  /** Coefficients array */
  private _coeffs: bigint[];
  /** Allocated length */
  private _alloc: number;
  /** Actual length (degree + 1) */
  private _length: number;
  /** Modular arithmetic context */
  private _mod: nmod_t;

  /**
   * Initialize polynomial with modulus.
   * Corresponds to nmod_poly_init.
   *
   * @param n - The modulus (must be > 0)
   */
  constructor(n: bigint | number) {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly constructor');
  }

  // ============================================================
  // Memory management
  // ============================================================

  /**
   * Initialize polynomial to zero with modulus n.
   * Corresponds to nmod_poly_init.
   */
  static init(n: bigint | number): nmod_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_init');
  }

  /**
   * Initialize with precomputed inverse.
   * Corresponds to nmod_poly_init_preinv.
   */
  static initPreinv(n: bigint | number, ninv: bigint): nmod_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_init_preinv');
  }

  /**
   * Initialize with pre-allocated space.
   * Corresponds to nmod_poly_init2.
   */
  static init2(n: bigint | number, alloc: number): nmod_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_init2');
  }

  /**
   * Initialize from nmod_t context.
   * Corresponds to nmod_poly_init_mod.
   */
  static initMod(mod: nmod_t): nmod_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_init_mod');
  }

  /**
   * Clear polynomial resources.
   * Corresponds to nmod_poly_clear.
   */
  clear(): void {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_clear');
  }

  /**
   * Reallocate to hold at least alloc coefficients.
   * Corresponds to nmod_poly_realloc.
   */
  realloc(alloc: number): void {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_realloc');
  }

  /**
   * Ensure space for at least alloc coefficients.
   * Corresponds to nmod_poly_fit_length.
   */
  fitLength(alloc: number): void {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_fit_length');
  }

  /**
   * Create a copy.
   */
  clone(): nmod_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly clone');
  }

  // ============================================================
  // Polynomial parameters
  // ============================================================

  /**
   * Get the length (degree + 1, or 0 for zero polynomial).
   * Corresponds to nmod_poly_length.
   */
  length(): number {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_length');
  }

  /**
   * Get the degree (length - 1, or -1 for zero polynomial).
   * Corresponds to nmod_poly_degree.
   */
  degree(): number {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_degree');
  }

  /**
   * Get the modulus.
   * Corresponds to nmod_poly_modulus.
   */
  modulus(): bigint {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_modulus');
  }

  /**
   * Get maximum coefficient bit size.
   * Corresponds to nmod_poly_max_bits.
   */
  maxBits(): number {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_max_bits');
  }

  /**
   * Get leading coefficient (null for zero polynomial).
   * Corresponds to nmod_poly_lead.
   */
  lead(): bigint | null {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_lead');
  }

  // ============================================================
  // Assignment and basic manipulation
  // ============================================================

  /**
   * Set from another polynomial.
   * Corresponds to nmod_poly_set.
   */
  set(other: nmod_poly): void {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_set');
  }

  /**
   * Swap two polynomials.
   * Corresponds to nmod_poly_swap.
   */
  swap(other: nmod_poly): void {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_swap');
  }

  /**
   * Set to zero polynomial.
   * Corresponds to nmod_poly_zero.
   */
  zero(): void {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_zero');
  }

  /**
   * Set to one (constant 1).
   * Corresponds to nmod_poly_one.
   */
  one(): void {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_one');
  }

  /**
   * Truncate to length n.
   * Corresponds to nmod_poly_truncate.
   */
  truncate(n: number): void {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_truncate');
  }

  /**
   * Set truncated copy.
   * Corresponds to nmod_poly_set_trunc.
   */
  static setTrunc(poly: nmod_poly, len: number): nmod_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_set_trunc');
  }

  /**
   * Reverse coefficients.
   * Corresponds to nmod_poly_reverse.
   */
  static reverse(poly: nmod_poly, m: number): nmod_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_reverse');
  }

  // ============================================================
  // Comparison
  // ============================================================

  /**
   * Check equality.
   * Corresponds to nmod_poly_equal.
   */
  equals(other: nmod_poly): boolean {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_equal');
  }

  /**
   * Check if zero.
   * Corresponds to nmod_poly_is_zero.
   */
  isZero(): boolean {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_is_zero');
  }

  /**
   * Check if one.
   * Corresponds to nmod_poly_is_one.
   */
  isOne(): boolean {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_is_one');
  }

  /**
   * Check if unit (non-zero constant).
   * Corresponds to nmod_poly_is_unit.
   */
  isUnit(): boolean {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_is_unit');
  }

  /**
   * Check if generator (x).
   * Corresponds to nmod_poly_is_gen.
   */
  isGen(): boolean {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_is_gen');
  }

  /**
   * Check if monic (leading coefficient is 1).
   * Corresponds to nmod_poly_is_monic.
   */
  isMonic(): boolean {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_is_monic');
  }

  // ============================================================
  // Getting and setting coefficients
  // ============================================================

  /**
   * Get coefficient at position j.
   * Corresponds to nmod_poly_get_coeff_ui.
   */
  getCoeff(j: number): bigint {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_get_coeff_ui');
  }

  /**
   * Set coefficient at position j.
   * Corresponds to nmod_poly_set_coeff_ui.
   */
  setCoeff(j: number, c: bigint | number): void {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_set_coeff_ui');
  }

  // ============================================================
  // Input and output
  // ============================================================

  /**
   * Get string representation.
   * Corresponds to nmod_poly_get_str.
   */
  toString(): string {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_get_str');
  }

  /**
   * Get pretty string representation.
   * Corresponds to nmod_poly_get_str_pretty.
   */
  toPrettyString(varname?: string): string {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_get_str_pretty');
  }

  /**
   * Set from string.
   * Corresponds to nmod_poly_set_str.
   */
  setStr(str: string): number {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_set_str');
  }

  // ============================================================
  // Shifting
  // ============================================================

  /**
   * Shift left (multiply by x^k).
   * Corresponds to nmod_poly_shift_left.
   */
  static shiftLeft(poly: nmod_poly, k: number): nmod_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_shift_left');
  }

  /**
   * Shift right (divide by x^k).
   * Corresponds to nmod_poly_shift_right.
   */
  static shiftRight(poly: nmod_poly, k: number): nmod_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_shift_right');
  }

  // ============================================================
  // Addition and subtraction
  // ============================================================

  /**
   * Add two polynomials.
   * Corresponds to nmod_poly_add.
   */
  static add(poly1: nmod_poly, poly2: nmod_poly): nmod_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_add');
  }

  /**
   * Subtract two polynomials.
   * Corresponds to nmod_poly_sub.
   */
  static sub(poly1: nmod_poly, poly2: nmod_poly): nmod_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_sub');
  }

  /**
   * Negate polynomial.
   * Corresponds to nmod_poly_neg.
   */
  static neg(poly: nmod_poly): nmod_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_neg');
  }

  /**
   * Add constant.
   * Corresponds to nmod_poly_add_ui.
   */
  static addUi(poly: nmod_poly, c: bigint | number): nmod_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_add_ui');
  }

  /**
   * Subtract constant.
   * Corresponds to nmod_poly_sub_ui.
   */
  static subUi(poly: nmod_poly, c: bigint | number): nmod_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_sub_ui');
  }

  // ============================================================
  // Scalar multiplication
  // ============================================================

  /**
   * Multiply by scalar.
   * Corresponds to nmod_poly_scalar_mul_nmod.
   */
  static scalarMul(poly: nmod_poly, c: bigint | number): nmod_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_scalar_mul_nmod');
  }

  /**
   * Make monic (divide by leading coefficient).
   * Corresponds to nmod_poly_make_monic.
   */
  static makeMonic(poly: nmod_poly): nmod_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_make_monic');
  }

  // ============================================================
  // Multiplication
  // ============================================================

  /**
   * Multiply two polynomials.
   * Corresponds to nmod_poly_mul.
   */
  static mul(poly1: nmod_poly, poly2: nmod_poly): nmod_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_mul');
  }

  /**
   * Multiply and truncate to n terms.
   * Corresponds to nmod_poly_mullow.
   */
  static mullow(poly1: nmod_poly, poly2: nmod_poly, n: number): nmod_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_mullow');
  }

  /**
   * Multiply and keep only high terms.
   * Corresponds to nmod_poly_mulhigh.
   */
  static mulhigh(poly1: nmod_poly, poly2: nmod_poly, n: number): nmod_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_mulhigh');
  }

  /**
   * Multiply modulo f: result = poly1 * poly2 mod f.
   * Corresponds to nmod_poly_mulmod.
   */
  static mulmod(poly1: nmod_poly, poly2: nmod_poly, f: nmod_poly): nmod_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_mulmod');
  }

  /**
   * Square a polynomial.
   * Corresponds to nmod_poly_sqr (via pow).
   */
  static sqr(poly: nmod_poly): nmod_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_sqr');
  }

  // ============================================================
  // Powering
  // ============================================================

  /**
   * Compute polynomial power.
   * Corresponds to nmod_poly_pow.
   */
  static pow(poly: nmod_poly, e: bigint | number): nmod_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_pow');
  }

  /**
   * Compute power and truncate.
   * Corresponds to nmod_poly_pow_trunc.
   */
  static powTrunc(poly: nmod_poly, e: bigint | number, n: number): nmod_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_pow_trunc');
  }

  /**
   * Compute power modulo f: result = poly^e mod f.
   * Corresponds to nmod_poly_powmod_ui_binexp.
   */
  static powmod(poly: nmod_poly, e: bigint | number, f: nmod_poly): nmod_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_powmod_ui_binexp');
  }

  /**
   * Compute x^e mod f efficiently.
   * Corresponds to nmod_poly_powmod_x_ui_preinv.
   */
  static powmodX(e: bigint | number, f: nmod_poly): nmod_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_powmod_x_ui_preinv');
  }

  // ============================================================
  // Division
  // ============================================================

  /**
   * Division with remainder.
   * Corresponds to nmod_poly_divrem.
   *
   * @returns [quotient, remainder]
   */
  static divrem(A: nmod_poly, B: nmod_poly): [nmod_poly, nmod_poly] {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_divrem');
  }

  /**
   * Division quotient.
   * Corresponds to nmod_poly_div.
   */
  static div(A: nmod_poly, B: nmod_poly): nmod_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_div');
  }

  /**
   * Division remainder.
   * Corresponds to nmod_poly_rem.
   */
  static rem(A: nmod_poly, B: nmod_poly): nmod_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_rem');
  }

  /**
   * Exact division.
   * Corresponds to nmod_poly_divexact.
   */
  static divexact(A: nmod_poly, B: nmod_poly): nmod_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_divexact');
  }

  /**
   * Division by (x - c).
   * Corresponds to nmod_poly_div_root.
   *
   * @returns [quotient, evaluation_at_c]
   */
  static divRoot(A: nmod_poly, c: bigint | number): [nmod_poly, bigint] {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_div_root');
  }

  // ============================================================
  // Power series operations
  // ============================================================

  /**
   * Power series inverse mod x^n.
   * Corresponds to nmod_poly_inv_series.
   */
  static invSeries(Q: nmod_poly, n: number): nmod_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_inv_series');
  }

  /**
   * Power series division mod x^n.
   * Corresponds to nmod_poly_div_series.
   */
  static divSeries(A: nmod_poly, B: nmod_poly, n: number): nmod_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_div_series');
  }

  // ============================================================
  // Divisibility
  // ============================================================

  /**
   * Check if B divides A.
   * Corresponds to nmod_poly_divides.
   *
   * @returns quotient if B divides A, null otherwise
   */
  static divides(A: nmod_poly, B: nmod_poly): nmod_poly | null {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_divides');
  }

  /**
   * Remove highest power of p dividing f.
   * Corresponds to nmod_poly_remove.
   *
   * @returns multiplicity of p in f
   */
  static remove(f: nmod_poly, p: nmod_poly): [nmod_poly, bigint] {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_remove');
  }

  // ============================================================
  // Derivative and integral
  // ============================================================

  /**
   * Derivative.
   * Corresponds to nmod_poly_derivative.
   */
  static derivative(poly: nmod_poly): nmod_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_derivative');
  }

  /**
   * Integral.
   * Corresponds to nmod_poly_integral.
   */
  static integral(poly: nmod_poly): nmod_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_integral');
  }

  // ============================================================
  // Evaluation
  // ============================================================

  /**
   * Evaluate at point.
   * Corresponds to nmod_poly_evaluate_nmod.
   */
  static evaluate(poly: nmod_poly, c: bigint | number): bigint {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_evaluate_nmod');
  }

  /**
   * Evaluate at multiple points.
   * Corresponds to nmod_poly_evaluate_nmod_vec.
   */
  static evaluateVec(poly: nmod_poly, xs: bigint[]): bigint[] {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_evaluate_nmod_vec');
  }

  // ============================================================
  // Composition
  // ============================================================

  /**
   * Compose: compute poly1(poly2(x)).
   * Corresponds to nmod_poly_compose.
   */
  static compose(poly1: nmod_poly, poly2: nmod_poly): nmod_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_compose');
  }

  /**
   * Compose modulo h: compute poly1(poly2(x)) mod h.
   * Corresponds to nmod_poly_compose_mod.
   */
  static composeMod(poly1: nmod_poly, poly2: nmod_poly, h: nmod_poly): nmod_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_compose_mod');
  }

  /**
   * Taylor shift: compute poly(x + c).
   * Corresponds to nmod_poly_taylor_shift.
   */
  static taylorShift(poly: nmod_poly, c: bigint | number): nmod_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_taylor_shift');
  }

  // ============================================================
  // Power series composition
  // ============================================================

  /**
   * Power series composition mod x^n.
   * Corresponds to nmod_poly_compose_series.
   */
  static composeSeries(poly1: nmod_poly, poly2: nmod_poly, n: number): nmod_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_compose_series');
  }

  /**
   * Compositional inverse mod x^n.
   * Corresponds to nmod_poly_revert_series.
   */
  static revertSeries(poly: nmod_poly, n: number): nmod_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_revert_series');
  }

  // ============================================================
  // GCD and related
  // ============================================================

  /**
   * Greatest common divisor.
   * Corresponds to nmod_poly_gcd.
   */
  static gcd(poly1: nmod_poly, poly2: nmod_poly): nmod_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_gcd');
  }

  /**
   * Extended GCD.
   * Corresponds to nmod_poly_xgcd.
   *
   * @returns [gcd, s, t] where gcd = s*poly1 + t*poly2
   */
  static xgcd(poly1: nmod_poly, poly2: nmod_poly): [nmod_poly, nmod_poly, nmod_poly] {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_xgcd');
  }

  /**
   * GCD with inverse.
   * Corresponds to nmod_poly_gcdinv.
   *
   * @returns [gcd, s] where gcd = s*A mod B
   */
  static gcdinv(A: nmod_poly, B: nmod_poly): [nmod_poly, nmod_poly] {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_gcdinv');
  }

  /**
   * Resultant.
   * Corresponds to nmod_poly_resultant.
   */
  static resultant(poly1: nmod_poly, poly2: nmod_poly): bigint {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_resultant');
  }

  /**
   * Discriminant.
   * Corresponds to nmod_poly_discriminant.
   */
  static discriminant(poly: nmod_poly): bigint {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_discriminant');
  }

  // ============================================================
  // Modular inverse
  // ============================================================

  /**
   * Modular inverse: compute A^(-1) mod P.
   * Corresponds to nmod_poly_invmod.
   *
   * @returns inverse if it exists, null otherwise
   */
  static invmod(A: nmod_poly, P: nmod_poly): nmod_poly | null {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_invmod');
  }

  // ============================================================
  // Square roots
  // ============================================================

  /**
   * Polynomial square root.
   * Corresponds to nmod_poly_sqrt.
   *
   * @returns square root if it exists, null otherwise
   */
  static sqrt(poly: nmod_poly): nmod_poly | null {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_sqrt');
  }

  /**
   * Power series square root mod x^n.
   * Corresponds to nmod_poly_sqrt_series.
   */
  static sqrtSeries(poly: nmod_poly, n: number): nmod_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_sqrt_series');
  }

  /**
   * Power series inverse square root mod x^n.
   * Corresponds to nmod_poly_invsqrt_series.
   */
  static invsqrtSeries(poly: nmod_poly, n: number): nmod_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_invsqrt_series');
  }

  // ============================================================
  // Interpolation
  // ============================================================

  /**
   * Interpolate polynomial through points.
   * Corresponds to nmod_poly_interpolate_nmod_vec.
   */
  static interpolate(n: bigint | number, xs: bigint[], ys: bigint[]): nmod_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_interpolate_nmod_vec');
  }

  // ============================================================
  // Products from roots
  // ============================================================

  /**
   * Product of (x - roots[i]).
   * Corresponds to nmod_poly_product_roots_nmod_vec.
   */
  static productRoots(n: bigint | number, roots: bigint[]): nmod_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_product_roots_nmod_vec');
  }

  // ============================================================
  // Inflation and deflation
  // ============================================================

  /**
   * Get deflation (gcd of exponents with non-zero coefficients).
   * Corresponds to nmod_poly_deflation.
   */
  static deflation(poly: nmod_poly): number {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_deflation');
  }

  /**
   * Deflate: poly(x) -> poly(x^(1/d)).
   * Corresponds to nmod_poly_deflate.
   */
  static deflate(poly: nmod_poly, d: number): nmod_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_deflate');
  }

  /**
   * Inflate: poly(x) -> poly(x^d).
   * Corresponds to nmod_poly_inflate.
   */
  static inflate(poly: nmod_poly, d: number): nmod_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_inflate');
  }

  // ============================================================
  // Irreducibility and special polynomials
  // ============================================================

  /**
   * Generate random monic irreducible polynomial.
   * Corresponds to nmod_poly_randtest_monic_irreducible.
   */
  static randomIrreducible(n: bigint | number, degree: number): nmod_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_randtest_monic_irreducible');
  }

  /**
   * Generate random monic primitive polynomial.
   * Corresponds to nmod_poly_randtest_monic_primitive.
   */
  static randomPrimitive(n: bigint | number, degree: number): nmod_poly {
    throw new Error('FLINT_NOT_IMPLEMENTED: nmod_poly_randtest_monic_primitive');
  }

  /**
   * Get Conway polynomial (standard irreducible).
   * Corresponds to _nmod_poly_conway.
   *
   * @returns Conway polynomial Phi_{p,d} if it exists, null otherwise
   */
  static conway(p: bigint | number, d: number): nmod_poly | null {
    throw new Error('FLINT_NOT_IMPLEMENTED: _nmod_poly_conway');
  }
}

// ============================================================
// Functional API (FLINT-style)
// ============================================================

/**
 * Initialize nmod_poly.
 * @see nmod_poly.init
 */
export function nmod_poly_init(n: bigint | number): nmod_poly {
  return nmod_poly.init(n);
}

/**
 * Clear nmod_poly.
 * @see nmod_poly.clear
 */
export function nmod_poly_clear(poly: nmod_poly): void {
  poly.clear();
}

/**
 * Add two nmod_poly.
 * @see nmod_poly.add
 */
export function nmod_poly_add(poly1: nmod_poly, poly2: nmod_poly): nmod_poly {
  return nmod_poly.add(poly1, poly2);
}

/**
 * Subtract two nmod_poly.
 * @see nmod_poly.sub
 */
export function nmod_poly_sub(poly1: nmod_poly, poly2: nmod_poly): nmod_poly {
  return nmod_poly.sub(poly1, poly2);
}

/**
 * Multiply two nmod_poly.
 * @see nmod_poly.mul
 */
export function nmod_poly_mul(poly1: nmod_poly, poly2: nmod_poly): nmod_poly {
  return nmod_poly.mul(poly1, poly2);
}

/**
 * GCD of two nmod_poly.
 * @see nmod_poly.gcd
 */
export function nmod_poly_gcd(poly1: nmod_poly, poly2: nmod_poly): nmod_poly {
  return nmod_poly.gcd(poly1, poly2);
}

/**
 * Power of nmod_poly.
 * @see nmod_poly.pow
 */
export function nmod_poly_pow(poly: nmod_poly, e: bigint | number): nmod_poly {
  return nmod_poly.pow(poly, e);
}

/**
 * Evaluate nmod_poly at point.
 * @see nmod_poly.evaluate
 */
export function nmod_poly_evaluate_nmod(poly: nmod_poly, c: bigint | number): bigint {
  return nmod_poly.evaluate(poly, c);
}
