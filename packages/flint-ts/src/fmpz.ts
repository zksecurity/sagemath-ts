/**
 * FLINT arbitrary precision integers (fmpz)
 *
 * The fmpz type is FLINT's main arbitrary precision integer type.
 * It uses a small/large representation: small integers (fitting in a signed word)
 * are stored inline, while larger integers use GMP's mpz_t.
 *
 * @see Reference: flint/src/fmpz.h, flint/src/fmpz_types.h
 * @module
 */

/**
 * Result of integer factorization.
 * Contains prime factors and their exponents.
 *
 * @see Reference: flint/src/fmpz_types.h (fmpz_factor_struct)
 */
export interface fmpz_factor {
  /** Sign of the factored integer (-1, 0, or 1) */
  readonly sign: -1 | 0 | 1;
  /** Prime factors */
  readonly primes: readonly fmpz[];
  /** Corresponding exponents */
  readonly exponents: readonly bigint[];
  /** Number of prime factors */
  readonly num: number;
}

/**
 * FLINT arbitrary precision integer.
 *
 * In FLINT, fmpz uses a tagged representation:
 * - Small integers (|n| <= COEFF_MAX) are stored inline
 * - Large integers are stored as pointers to mpz_t structures
 *
 * This TypeScript implementation uses bigint internally.
 *
 * @see Reference: flint/src/fmpz.h
 */
export class fmpz {
  /** Internal value storage */
  private _value: bigint;

  /**
   * Initialize a new fmpz.
   * Corresponds to fmpz_init / fmpz_init_set.
   *
   * @param value - Initial value (default: 0)
   * @throws Error - FLINT_NOT_IMPLEMENTED
   */
  constructor(value?: bigint | number | string) {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz constructor');
    // Implementation would be:
    // this._value = value !== undefined ? BigInt(value) : 0n;
  }

  /**
   * Get the value as bigint.
   */
  get value(): bigint {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz.value getter');
  }

  // ============================================================
  // Memory management
  // ============================================================

  /**
   * Initialize fmpz to zero.
   * Corresponds to fmpz_init.
   *
   * @returns New fmpz initialized to zero
   */
  static init(): fmpz {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_init');
  }

  /**
   * Clear/free fmpz resources.
   * Corresponds to fmpz_clear.
   *
   * In TypeScript this is a no-op due to garbage collection,
   * but included for API compatibility.
   */
  clear(): void {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_clear');
  }

  /**
   * Create a copy of this fmpz.
   * Corresponds to fmpz_init_set.
   */
  clone(): fmpz {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_init_set');
  }

  // ============================================================
  // Assignment and conversion
  // ============================================================

  /**
   * Set fmpz to zero.
   * Corresponds to fmpz_zero.
   */
  zero(): void {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_zero');
  }

  /**
   * Set fmpz to one.
   * Corresponds to fmpz_one.
   */
  one(): void {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_one');
  }

  /**
   * Set value from another fmpz.
   * Corresponds to fmpz_set.
   */
  set(other: fmpz): void {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_set');
  }

  /**
   * Set value from bigint/number.
   * Corresponds to fmpz_set_si / fmpz_set_ui.
   */
  setInt(value: bigint | number): void {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_set_si');
  }

  /**
   * Set value from string.
   * Corresponds to fmpz_set_str.
   *
   * @param str - String representation
   * @param base - Base (2-36, default 10)
   * @returns 0 on success, -1 on failure
   */
  setStr(str: string, base?: number): number {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_set_str');
  }

  /**
   * Get string representation.
   * Corresponds to fmpz_get_str.
   *
   * @param base - Base (2-36, default 10)
   */
  toString(base?: number): string {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_get_str');
  }

  /**
   * Get as JavaScript number (may lose precision).
   * Corresponds to fmpz_get_d.
   */
  toNumber(): number {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_get_d');
  }

  // ============================================================
  // Comparison
  // ============================================================

  /**
   * Check if zero.
   * Corresponds to fmpz_is_zero.
   */
  isZero(): boolean {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_is_zero');
  }

  /**
   * Check if one.
   * Corresponds to fmpz_is_one.
   */
  isOne(): boolean {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_is_one');
  }

  /**
   * Check if plus or minus one.
   * Corresponds to fmpz_is_pm1.
   */
  isPm1(): boolean {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_is_pm1');
  }

  /**
   * Check if even.
   * Corresponds to fmpz_is_even.
   */
  isEven(): boolean {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_is_even');
  }

  /**
   * Check if odd.
   * Corresponds to fmpz_is_odd.
   */
  isOdd(): boolean {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_is_odd');
  }

  /**
   * Compare two fmpz values.
   * Corresponds to fmpz_cmp.
   *
   * @returns negative if this < other, 0 if equal, positive if this > other
   */
  cmp(other: fmpz): number {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_cmp');
  }

  /**
   * Check equality.
   * Corresponds to fmpz_equal.
   */
  equals(other: fmpz): boolean {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_equal');
  }

  /**
   * Get sign.
   * Corresponds to fmpz_sgn.
   *
   * @returns -1, 0, or 1
   */
  sgn(): -1 | 0 | 1 {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_sgn');
  }

  // ============================================================
  // Basic properties
  // ============================================================

  /**
   * Get number of bits.
   * Corresponds to fmpz_bits.
   */
  bits(): bigint {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_bits');
  }

  /**
   * Get size in base b.
   * Corresponds to fmpz_sizeinbase.
   */
  sizeInBase(b: number): number {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_sizeinbase');
  }

  /**
   * Get valuation at 2 (largest power of 2 dividing this).
   * Corresponds to fmpz_val2.
   */
  val2(): bigint {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_val2');
  }

  // ============================================================
  // Basic arithmetic
  // ============================================================

  /**
   * Add two fmpz values: result = a + b.
   * Corresponds to fmpz_add.
   */
  static add(a: fmpz, b: fmpz): fmpz {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_add');
  }

  /**
   * Subtract two fmpz values: result = a - b.
   * Corresponds to fmpz_sub.
   */
  static sub(a: fmpz, b: fmpz): fmpz {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_sub');
  }

  /**
   * Multiply two fmpz values: result = a * b.
   * Corresponds to fmpz_mul.
   */
  static mul(a: fmpz, b: fmpz): fmpz {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_mul');
  }

  /**
   * Add multiply: result = f + g * h.
   * Corresponds to fmpz_addmul.
   */
  static addmul(f: fmpz, g: fmpz, h: fmpz): fmpz {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_addmul');
  }

  /**
   * Subtract multiply: result = f - g * h.
   * Corresponds to fmpz_submul.
   */
  static submul(f: fmpz, g: fmpz, h: fmpz): fmpz {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_submul');
  }

  /**
   * Negate: result = -a.
   * Corresponds to fmpz_neg.
   */
  static neg(a: fmpz): fmpz {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_neg');
  }

  /**
   * Absolute value: result = |a|.
   * Corresponds to fmpz_abs.
   */
  static abs(a: fmpz): fmpz {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_abs');
  }

  // ============================================================
  // Division
  // ============================================================

  /**
   * Truncated division quotient: q = trunc(a / b).
   * Corresponds to fmpz_tdiv_q.
   */
  static tdiv_q(a: fmpz, b: fmpz): fmpz {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_tdiv_q');
  }

  /**
   * Truncated division with remainder.
   * Corresponds to fmpz_tdiv_qr.
   *
   * @returns [quotient, remainder] where a = q*b + r
   */
  static tdiv_qr(a: fmpz, b: fmpz): [fmpz, fmpz] {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_tdiv_qr');
  }

  /**
   * Floor division quotient: q = floor(a / b).
   * Corresponds to fmpz_fdiv_q.
   */
  static fdiv_q(a: fmpz, b: fmpz): fmpz {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_fdiv_q');
  }

  /**
   * Floor division with remainder.
   * Corresponds to fmpz_fdiv_qr.
   *
   * @returns [quotient, remainder]
   */
  static fdiv_qr(a: fmpz, b: fmpz): [fmpz, fmpz] {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_fdiv_qr');
  }

  /**
   * Ceiling division quotient: q = ceil(a / b).
   * Corresponds to fmpz_cdiv_q.
   */
  static cdiv_q(a: fmpz, b: fmpz): fmpz {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_cdiv_q');
  }

  /**
   * Exact division (assumes b divides a exactly).
   * Corresponds to fmpz_divexact.
   */
  static divexact(a: fmpz, b: fmpz): fmpz {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_divexact');
  }

  /**
   * Check if b divides a.
   * Corresponds to fmpz_divisible.
   */
  static divisible(a: fmpz, b: fmpz): boolean {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_divisible');
  }

  /**
   * Compute quotient if division is exact.
   * Corresponds to fmpz_divides.
   *
   * @returns quotient if b divides a, null otherwise
   */
  static divides(a: fmpz, b: fmpz): fmpz | null {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_divides');
  }

  /**
   * Modulo operation: result = a mod b (always non-negative if b > 0).
   * Corresponds to fmpz_mod.
   */
  static mod(a: fmpz, b: fmpz): fmpz {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_mod');
  }

  /**
   * Symmetric modulo: result in (-b/2, b/2].
   * Corresponds to fmpz_smod.
   */
  static smod(a: fmpz, b: fmpz): fmpz {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_smod');
  }

  // ============================================================
  // Powers and roots
  // ============================================================

  /**
   * Power: result = base^exp.
   * Corresponds to fmpz_pow_ui.
   */
  static pow(base: fmpz, exp: bigint | number): fmpz {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_pow_ui');
  }

  /**
   * Integer square root: result = floor(sqrt(a)).
   * Corresponds to fmpz_sqrt.
   */
  static sqrt(a: fmpz): fmpz {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_sqrt');
  }

  /**
   * Integer square root with remainder.
   * Corresponds to fmpz_sqrtrem.
   *
   * @returns [sqrt, remainder] where a = sqrt^2 + remainder
   */
  static sqrtrem(a: fmpz): [fmpz, fmpz] {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_sqrtrem');
  }

  /**
   * Check if perfect square.
   * Corresponds to fmpz_is_square.
   */
  static isSquare(a: fmpz): boolean {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_is_square');
  }

  /**
   * Integer n-th root: result = floor(a^(1/n)).
   * Corresponds to fmpz_root.
   *
   * @returns [root, exact] where exact is true if a = root^n
   */
  static root(a: fmpz, n: number): [fmpz, boolean] {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_root');
  }

  /**
   * Check if perfect power.
   * Corresponds to fmpz_is_perfect_power.
   *
   * @returns [root, true] if a = root^k for some k > 1, [a, false] otherwise
   */
  static isPerfectPower(a: fmpz): [fmpz, boolean] {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_is_perfect_power');
  }

  // ============================================================
  // GCD and LCM
  // ============================================================

  /**
   * Greatest common divisor.
   * Corresponds to fmpz_gcd.
   */
  static gcd(a: fmpz, b: fmpz): fmpz {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_gcd');
  }

  /**
   * Least common multiple.
   * Corresponds to fmpz_lcm.
   */
  static lcm(a: fmpz, b: fmpz): fmpz {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_lcm');
  }

  /**
   * Extended GCD: returns [d, a, b] where d = gcd(f, g) = a*f + b*g.
   * Corresponds to fmpz_xgcd.
   */
  static xgcd(f: fmpz, g: fmpz): [fmpz, fmpz, fmpz] {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_xgcd');
  }

  /**
   * GCD with inverse: returns [d, a] where d = gcd(f, g) and a*f = d mod g.
   * Corresponds to fmpz_gcdinv.
   */
  static gcdinv(f: fmpz, g: fmpz): [fmpz, fmpz] {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_gcdinv');
  }

  // ============================================================
  // Modular arithmetic
  // ============================================================

  /**
   * Modular inverse: result = a^(-1) mod m.
   * Corresponds to fmpz_invmod.
   *
   * @returns inverse if it exists, null otherwise
   */
  static invmod(a: fmpz, m: fmpz): fmpz | null {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_invmod');
  }

  /**
   * Modular exponentiation: result = base^exp mod m.
   * Corresponds to fmpz_powm.
   */
  static powm(base: fmpz, exp: fmpz, m: fmpz): fmpz {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_powm');
  }

  /**
   * Modular square root: find x such that x^2 = a mod p.
   * Corresponds to fmpz_sqrtmod.
   *
   * @returns square root if it exists, null otherwise
   */
  static sqrtmod(a: fmpz, p: fmpz): fmpz | null {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_sqrtmod');
  }

  // ============================================================
  // Primality testing
  // ============================================================

  /**
   * Probabilistic primality test (BPSW).
   * Corresponds to fmpz_is_probabprime.
   */
  static isProbabPrime(n: fmpz): boolean {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_is_probabprime');
  }

  /**
   * Proven primality test.
   * Corresponds to fmpz_is_prime.
   */
  static isPrime(n: fmpz): boolean {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_is_prime');
  }

  /**
   * Find next prime >= n.
   * Corresponds to fmpz_nextprime.
   *
   * @param proved - If true, use proven primality test
   */
  static nextprime(n: fmpz, proved?: boolean): fmpz {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_nextprime');
  }

  // ============================================================
  // Factorization
  // ============================================================

  /**
   * Factor an integer.
   * Corresponds to fmpz_factor.
   *
   * @see Reference: flint/src/fmpz_factor.h
   */
  static factor(n: fmpz): fmpz_factor {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_factor');
  }

  /**
   * Expand a factorization back to the integer.
   * Corresponds to fmpz_factor_expand.
   */
  static factorExpand(factor: fmpz_factor): fmpz {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_factor_expand');
  }

  // ============================================================
  // Number theoretic functions
  // ============================================================

  /**
   * Jacobi symbol (a/p).
   * Corresponds to fmpz_jacobi.
   */
  static jacobi(a: fmpz, p: fmpz): -1 | 0 | 1 {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_jacobi');
  }

  /**
   * Kronecker symbol (a/n).
   * Corresponds to fmpz_kronecker.
   */
  static kronecker(a: fmpz, n: fmpz): -1 | 0 | 1 {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_kronecker');
  }

  /**
   * Euler's totient function.
   * Corresponds to fmpz_euler_phi.
   */
  static eulerPhi(n: fmpz): fmpz {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_euler_phi');
  }

  /**
   * Moebius function.
   * Corresponds to fmpz_moebius_mu.
   */
  static moebiusMu(n: fmpz): -1 | 0 | 1 {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_moebius_mu');
  }

  /**
   * Divisor sigma function sigma_k(n).
   * Corresponds to fmpz_divisor_sigma.
   */
  static divisorSigma(k: bigint | number, n: fmpz): fmpz {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_divisor_sigma');
  }

  // ============================================================
  // Combinatorics
  // ============================================================

  /**
   * Factorial: result = n!.
   * Corresponds to fmpz_fac_ui.
   */
  static factorial(n: bigint | number): fmpz {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_fac_ui');
  }

  /**
   * Binomial coefficient: result = C(n, k).
   * Corresponds to fmpz_bin_uiui.
   */
  static binomial(n: bigint | number, k: bigint | number): fmpz {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_bin_uiui');
  }

  /**
   * Fibonacci number: result = F_n.
   * Corresponds to fmpz_fib_ui.
   */
  static fibonacci(n: bigint | number): fmpz {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_fib_ui');
  }

  /**
   * Primorial: result = n# (product of primes <= n).
   * Corresponds to fmpz_primorial.
   */
  static primorial(n: bigint | number): fmpz {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_primorial');
  }

  // ============================================================
  // Bitwise operations
  // ============================================================

  /**
   * Bitwise AND.
   * Corresponds to fmpz_and.
   */
  static and(a: fmpz, b: fmpz): fmpz {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_and');
  }

  /**
   * Bitwise OR.
   * Corresponds to fmpz_or.
   */
  static or(a: fmpz, b: fmpz): fmpz {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_or');
  }

  /**
   * Bitwise XOR.
   * Corresponds to fmpz_xor.
   */
  static xor(a: fmpz, b: fmpz): fmpz {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_xor');
  }

  /**
   * Bitwise complement.
   * Corresponds to fmpz_complement.
   */
  static complement(a: fmpz): fmpz {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_complement');
  }

  /**
   * Left shift: result = a * 2^exp.
   * Corresponds to fmpz_mul_2exp.
   */
  static shiftLeft(a: fmpz, exp: bigint | number): fmpz {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_mul_2exp');
  }

  /**
   * Right shift (floor division by 2^exp).
   * Corresponds to fmpz_fdiv_q_2exp.
   */
  static shiftRight(a: fmpz, exp: bigint | number): fmpz {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_fdiv_q_2exp');
  }

  /**
   * Test bit at position i.
   * Corresponds to fmpz_tstbit.
   */
  static testBit(a: fmpz, i: bigint | number): boolean {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_tstbit');
  }

  /**
   * Set bit at position i.
   * Corresponds to fmpz_setbit.
   */
  setBit(i: bigint | number): void {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_setbit');
  }

  /**
   * Clear bit at position i.
   * Corresponds to fmpz_clrbit.
   */
  clearBit(i: bigint | number): void {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_clrbit');
  }

  /**
   * Population count (number of set bits).
   * Corresponds to fmpz_popcnt.
   */
  static popcnt(a: fmpz): bigint {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_popcnt');
  }

  // ============================================================
  // CRT (Chinese Remainder Theorem)
  // ============================================================

  /**
   * Chinese Remainder Theorem combining.
   * Corresponds to fmpz_CRT.
   *
   * Given r1 mod m1 and r2 mod m2, compute r mod (m1*m2).
   *
   * @param sign - If true, result is in symmetric range
   */
  static CRT(r1: fmpz, m1: fmpz, r2: fmpz, m2: fmpz, sign?: boolean): fmpz {
    throw new Error('FLINT_NOT_IMPLEMENTED: fmpz_CRT');
  }
}

// ============================================================
// Functional API (FLINT-style)
// ============================================================

/**
 * Initialize fmpz to zero.
 * @see fmpz.init
 */
export function fmpz_init(): fmpz {
  return fmpz.init();
}

/**
 * Clear fmpz resources.
 * @see fmpz.clear
 */
export function fmpz_clear(f: fmpz): void {
  f.clear();
}

/**
 * Add two fmpz: f = g + h.
 * @see fmpz.add
 */
export function fmpz_add(g: fmpz, h: fmpz): fmpz {
  return fmpz.add(g, h);
}

/**
 * Subtract two fmpz: f = g - h.
 * @see fmpz.sub
 */
export function fmpz_sub(g: fmpz, h: fmpz): fmpz {
  return fmpz.sub(g, h);
}

/**
 * Multiply two fmpz: f = g * h.
 * @see fmpz.mul
 */
export function fmpz_mul(g: fmpz, h: fmpz): fmpz {
  return fmpz.mul(g, h);
}

/**
 * Truncated division: f = trunc(g / h).
 * @see fmpz.tdiv_q
 */
export function fmpz_tdiv_q(g: fmpz, h: fmpz): fmpz {
  return fmpz.tdiv_q(g, h);
}

/**
 * GCD of two fmpz.
 * @see fmpz.gcd
 */
export function fmpz_gcd(g: fmpz, h: fmpz): fmpz {
  return fmpz.gcd(g, h);
}

/**
 * LCM of two fmpz.
 * @see fmpz.lcm
 */
export function fmpz_lcm(g: fmpz, h: fmpz): fmpz {
  return fmpz.lcm(g, h);
}

/**
 * Check if n is prime.
 * @see fmpz.isPrime
 */
export function fmpz_is_prime(n: fmpz): boolean {
  return fmpz.isPrime(n);
}

/**
 * Factor an integer.
 * @see fmpz.factor
 */
export function fmpz_factor_func(n: fmpz): fmpz_factor {
  return fmpz.factor(n);
}
