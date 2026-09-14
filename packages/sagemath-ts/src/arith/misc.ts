/**
 * @module sage/arith/misc
 * @description Miscellaneous arithmetic functions
 *
 * Port of: sage/arith/misc.py
 * Reference: reference/sage/src/sage/arith/misc.py
 *
 * NOTE: SageMath delegates many functions to PARI/GP. This module follows
 * the same pattern, using @sagemath-ts/parigp-ts for the underlying algorithms.
 */

import {
  ArithmeticError,
  AssertionError,
  AttributeError,
  IndexError,
  RuntimeError,
  NotImplementedError,
  OverflowError,
  ValueError,
  ZeroDivisionError,
} from '../errors.js';
import { IntegerMatrix, LLL } from '../matrix/index.js';
import { current_randstate } from '../misc/randstate.js';
import { IntegerMod, Mod } from '../rings/finite_rings/integer_mod.js';
import { PrimeFieldElement } from '../rings/finite_rings/finite_field_extension.js';
import { FiniteFieldElement as LegacyPrimeElement } from '../rings/finite_rings/finite_field_prime.js';
import { two_squares_pyx, three_squares_pyx, four_squares_pyx } from '../rings/sum_of_squares.js';
import { ProductTree } from '../rings/generic.js';
import { arith_int, arith_llong } from '../rings/fast_arith.js';
import { RealNumber } from '../rings/real_mpfr.js';
import {
  mpfr_init2,
  mpfr_set_d,
  mpfr_set_z,
  mpfr_mul,
  mpfr_add,
  mpfr_round,
  mpfr_get_z,
  mpfr_cmp,
  type mpfr_t,
} from '@sagemath-ts/mpfr-ts';
import { Rational } from '../rings/rational.js';
import { Integer, ZZ } from '../rings/integer_ring.js';
import type { CoefficientRing, RingElement } from '../rings/polynomial/polynomial_element.js';
import { PolynomialRing } from '../rings/polynomial/polynomial_ring.js';
import { sorted as python_sorted } from '../types/python_sort.js';
import { type FloatInput, float as python_float } from '../types/python_float.js';
import {
  type IntegerLike,
  type RationalLike,
  toBigInt,
  toRational,
  toSafeNumber,
} from '../types/coercion.js';

// Import PARI functions for factorization and primality testing
// Reference: sage/arith/misc.py uses PARI via cypari2
import {
  type Factorization as PariFactorization,
  type MpReal,
  mpfactr as pari_mpfactr,
  algdep as pari_algdep,
  eulerphi as pari_eulerphi,
  numdiv as pari_numdiv,
  hilbert as pari_hilbert,
  Z_factor as pari_Z_factor,
  isPrime as pari_isPrime,
  isprimepower as pari_isprimepower,
  nextprime as pari_nextprime,
  prime as pari_prime,
  sumdedekind as pari_sumdedekind,
} from '@sagemath-ts/parigp-ts';
import { fmpq_dedekind_sum } from '@sagemath-ts/flint-ts';

/**
 * Return the non-negative remainder of a mod m.
 * Unlike JavaScript's %, this always returns a non-negative result.
 * @param a - Dividend
 * @param m - Modulus (must be positive)
 * @returns a mod m (always in range [0, m))
 */
function mod(a: bigint, m: bigint): bigint {
  const r = a % m;
  return r < 0n ? r + m : r;
}

/**
 * Return the greatest common divisor of a and b.
 *
 * The result is always non-negative.
 *
 * @param a - First integer
 * @param b - Second integer (or array if computing GCD of multiple values)
 * @returns The GCD of a and b
 *
 * @example
 * ```typescript
 * gcd(12n, 8n)  // 4n
 * gcd(-4n, 6n)  // 2n
 * gcd([2n, 4n, 6n, 8n])  // 2n
 * ```
 *
 * @see Reference: sage/arith/misc.py:gcd
 */
export function gcd(a: IntegerLike, b?: IntegerLike): bigint;
export function gcd(values: IntegerLike[]): bigint;
export function gcd(a: IntegerLike | IntegerLike[], b?: IntegerLike): bigint {
  if (Array.isArray(a)) {
    // GCD of a list
    if (a.length === 0) {
      return 0n;
    }
    let result = toBigInt(a[0]!);
    for (let i = 1; i < a.length; i++) {
      result = gcd(result, toBigInt(a[i]!));
      if (result === 1n) {
        return 1n; // Early exit optimization
      }
    }
    return result;
  }

  const _a = toBigInt(a);

  if (b === undefined) {
    throw new TypeError("'bigint' object is not iterable");
  }

  const _b = toBigInt(b);

  // Binary GCD algorithm (Stein's algorithm)
  // This matches the behavior of GMP's mpz_gcd
  let x = _a < 0n ? -_a : _a;
  let y = _b < 0n ? -_b : _b;

  if (x === 0n) return y;
  if (y === 0n) return x;

  // Find common factors of 2
  let shift = 0n;
  while (((x | y) & 1n) === 0n) {
    x >>= 1n;
    y >>= 1n;
    shift++;
  }

  // Remove remaining factors of 2 from x
  while ((x & 1n) === 0n) {
    x >>= 1n;
  }

  while (y !== 0n) {
    // Remove factors of 2 from y
    while ((y & 1n) === 0n) {
      y >>= 1n;
    }

    // Ensure x <= y
    if (x > y) {
      const temp = x;
      x = y;
      y = temp;
    }

    y -= x;
  }

  return x << shift;
}

// Alias for compatibility with SageMath
export const GCD = gcd;

/**
 * Return the least common multiple of a and b.
 *
 * @param a - First integer
 * @param b - Second integer
 * @returns The LCM of a and b
 *
 * @example
 * ```typescript
 * lcm(4n, 6n)  // 12n
 * lcm(0n, 5n)  // 0n
 * ```
 */
export function lcm(a: IntegerLike, b?: IntegerLike): bigint;
export function lcm(values: IntegerLike[]): bigint;
export function lcm(a: IntegerLike | IntegerLike[], b?: IntegerLike): bigint {
  if (Array.isArray(a)) {
    if (a.length === 0) {
      return 1n;
    }
    let result = toBigInt(a[0]!);
    for (let i = 1; i < a.length; i++) {
      result = lcm(result, toBigInt(a[i]!));
    }
    return result;
  }

  const _a = toBigInt(a);

  if (b === undefined) {
    throw new TypeError("'bigint' object is not iterable");
  }

  const _b = toBigInt(b);

  if (_a === 0n || _b === 0n) {
    return 0n;
  }

  const absA = _a < 0n ? -_a : _a;
  const absB = _b < 0n ? -_b : _b;

  return (absA / gcd(absA, absB)) * absB;
}

export const LCM = lcm;

/**
 * Return the extended gcd of a and b.
 *
 * Returns a triple (g, s, t) such that g = gcd(a, b) and g = s*a + t*b.
 *
 * @param a - First integer
 * @param b - Second integer
 * @returns A tuple [g, s, t] where g = gcd(a,b) = s*a + t*b
 *
 * @example
 * ```typescript
 * xgcd(6n, 4n)  // [2n, 1n, -1n]
 * // Verify: 2 = 1*6 + (-1)*4
 * ```
 *
 * @see Reference: sage/rings/integer.pyx:xgcd
 */
export function xgcd(a: IntegerLike, b: IntegerLike): [bigint, bigint, bigint] {
  const _a = toBigInt(a);
  const _b = toBigInt(b);
  // GMP's mpz_gcdext returns three zeros for (0,0).
  if (_a === 0n && _b === 0n) return [0n, 0n, 0n];
  // Extended Euclidean algorithm
  let oldR = _a;
  let r = _b;
  let oldS = 1n;
  let s = 0n;
  let oldT = 0n;
  let t = 1n;

  while (r !== 0n) {
    const quotient = oldR / r;

    const tempR = r;
    r = oldR - quotient * r;
    oldR = tempR;

    const tempS = s;
    s = oldS - quotient * s;
    oldS = tempS;

    const tempT = t;
    t = oldT - quotient * t;
    oldT = tempT;
  }

  // Ensure gcd is non-negative
  if (oldR < 0n) {
    return [-oldR, -oldS, -oldT];
  }

  return [oldR, oldS, oldT];
}

/**
 * Return the modular inverse of a modulo m.
 *
 * @param a - The integer to invert
 * @param m - The modulus
 * @returns The inverse of a modulo m
 * @throws {ZeroDivisionError} If gcd(a, m) != 1
 *
 * @example
 * ```typescript
 * inverse_mod(3n, 7n)  // 5n (since 3*5 = 15 ≡ 1 mod 7)
 * ```
 *
 * @see Reference: sage/arith/misc.py:inverse_mod
 */
export function inverse_mod(a: IntegerLike, m: IntegerLike): bigint {
  const _a = toBigInt(a);
  const _m = toBigInt(m);
  if (_m === 0n) {
    // Modulo the zero ideal only the integer units have inverses (mpz_invert).
    if (_a === 1n || _a === -1n) return _a;
    throw new ZeroDivisionError(`inverse of Mod(${_a}, ${_m}) does not exist`);
  }

  const absM = _m < 0n ? -_m : _m;
  const [g, s] = xgcd(_a, absM);

  if (g !== 1n) {
    throw new ZeroDivisionError(`inverse of Mod(${_a}, ${_m}) does not exist`);
  }

  // Ensure result is in [0, |m|)
  let result = s % absM;
  if (result < 0n) {
    result += absM;
  }

  return result;
}

/**
 * Return a^n mod m.
 *
 * @param a - Base
 * @param n - Exponent (can be negative if inverse exists)
 * @param m - Modulus
 * @returns a^n mod m
 *
 * @example
 * ```typescript
 * power_mod(2n, 10n, 1000n)  // 24n
 * power_mod(3n, -1n, 7n)     // 5n (modular inverse)
 * ```
 *
 * @see Reference: sage/arith/misc.py:power_mod
 */
export function power_mod(a: IntegerLike, n: IntegerLike, m: IntegerLike): bigint {
  let _a = toBigInt(a);
  let _n = toBigInt(n);
  const _m = toBigInt(m);
  if (_m === 0n) {
    throw new ZeroDivisionError('modulus must be nonzero');
  }

  if (_m === 1n || _m === -1n) {
    return 0n;
  }

  const absM = _m < 0n ? -_m : _m;

  // Handle negative exponent
  if (_n < 0n) {
    _a = inverse_mod(_a, absM);
    _n = -_n;
  }

  // Reduce a modulo m first
  _a = ((_a % absM) + absM) % absM;

  if (_n === 0n) {
    return 1n;
  }

  // Binary exponentiation (square-and-multiply)
  let result = 1n;
  let base = _a;

  while (_n > 0n) {
    if ((_n & 1n) === 1n) {
      result = (result * base) % absM;
    }
    base = (base * base) % absM;
    _n >>= 1n;
  }

  return result;
}

/** Largest value of a C ``signed long``; SageMath's default ``bound``. */
const LONG_MAX = 9223372036854775807n;

/**
 * Return the smallest prime divisor of n up to `bound`, beginning checking at
 * `start`, or |n| if no such divisor is found.
 *
 * This is a direct transcription of SageMath's
 * ``Integer.trial_division`` (`sage/rings/integer.pyx:3828`), which trial
 * divides by 2, 3, 5 and then by every integer congruent to
 * 1, 7, 11, 13, 17, 19, 23, 29 mod 30 up to min(bound, isqrt(|n|)).
 * It does **not** factor n.
 *
 * @param n - Nonzero integer
 * @param bound - Positive upper bound for the trial divisors (default: LONG_MAX)
 * @param start - Positive integer to start checking at (default: 2)
 * @returns The smallest prime factor ≤ bound, or |n| if no such factor exists
 * @throws {ValueError} if bound <= 0 or n == 0
 *
 * @example
 * ```typescript
 * trial_division(15n)        // 3n
 * trial_division(387833n, 300n)  // 387833n
 * trial_division(387833n, 400n)  // 389n
 * trial_division(3n * 5n * 101n * 103n, undefined, 50n)  // 101n
 * ```
 *
 * @see Reference: sage/arith/misc.py:trial_division -> sage/rings/integer.pyx:trial_division
 */
export function trial_division(
  n: IntegerLike,
  bound?: IntegerLike,
  start: IntegerLike = 2n
): bigint {
  const _n = toBigInt(n);
  const _bound = bound !== undefined ? toBigInt(bound) : LONG_MAX;
  const _start = toBigInt(start);

  if (_bound <= 0n) {
    throw new ValueError('bound must be positive');
  }
  if (_n === 0n) {
    throw new ValueError('self must be nonzero');
  }

  const N = _n < 0n ? -_n : _n;
  if (N === 1n) {
    return 1n;
  }

  if (_start <= 2n && N % 2n === 0n) {
    return 2n;
  }
  if (_start <= 3n && N % 3n === 0n) {
    return 3n;
  }
  if (_start <= 5n && N % 5n === 0n) {
    return 5n;
  }

  // Only trial divide by numbers congruent to 1,7,11,13,17,19,23,29 mod 30.
  const dif = [6n, 4n, 2n, 4n, 2n, 4n, 6n, 2n];
  let m = 7n;
  let i = 0;
  if (_start > 7n) {
    // Find the wheel position corresponding to ``start``.
    const r = _start % 30n;
    if (r <= 1n) {
      i = 0;
      m = _start + (1n - r);
    } else if (r <= 7n) {
      i = 1;
      m = _start + (7n - r);
    } else if (r <= 11n) {
      i = 2;
      m = _start + (11n - r);
    } else if (r <= 13n) {
      i = 3;
      m = _start + (13n - r);
    } else if (r <= 17n) {
      i = 4;
      m = _start + (17n - r);
    } else if (r <= 19n) {
      i = 5;
      m = _start + (19n - r);
    } else if (r <= 23n) {
      i = 6;
      m = _start + (23n - r);
    } else {
      i = 7;
      m = _start + (29n - r);
    }
  } else {
    // SageMath starts at m = 7 with i = 1 (so the first step is dif[1] = 4).
    i = 1;
  }

  const sqrtN = isqrt(N);
  const limit = sqrtN < _bound ? sqrtN : _bound;

  while (m <= limit) {
    if (N % m === 0n) {
      return m;
    }
    m += dif[i % 8]!;
    i++;
  }

  return N;
}

/**
 * Integer square root (floor of square root).
 *
 * @param n - Non-negative integer
 * @returns Floor of the square root of n
 */
export function isqrt(n: IntegerLike): bigint {
  const _n = toBigInt(n);
  if (_n < 0n) {
    throw new ValueError('isqrt() argument must be nonnegative');
  }

  if (_n < 2n) {
    return _n;
  }

  // Newton's method
  let x = _n;
  let y = (x + 1n) >> 1n;

  while (y < x) {
    x = y;
    y = (x + _n / x) >> 1n;
  }

  return x;
}

/**
 * Miller-Rabin primality test for a single witness.
 *
 * @param n - Number to test
 * @param a - Witness
 * @returns false if n is definitely composite, true if probably prime
 */
function millerRabinWitness(n: bigint, a: bigint): boolean {
  // Write n-1 = 2^r * d
  let d = n - 1n;
  let r = 0n;

  while ((d & 1n) === 0n) {
    d >>= 1n;
    r++;
  }

  // Compute a^d mod n
  let x = power_mod(a, d, n);

  if (x === 1n || x === n - 1n) {
    return true;
  }

  for (let i = 1n; i < r; i++) {
    x = (x * x) % n;
    if (x === n - 1n) {
      return true;
    }
    if (x === 1n) {
      return false;
    }
  }

  return false;
}

/**
 * Test if n is a strong probable prime to base a (Miller-Rabin test).
 *
 * @param n - Integer to test
 * @param a - Base for the test
 * @returns true if n is a strong probable prime to base a
 */
export function is_strong_probable_prime(n: IntegerLike, a: IntegerLike): boolean {
  const _n = toBigInt(n);
  let _a = toBigInt(a);
  if (_n < 2n) {
    return false;
  }

  if (_n === 2n) {
    return true;
  }

  if ((_n & 1n) === 0n) {
    return false;
  }

  _a = ((_a % _n) + _n) % _n;
  if (_a === 0n) {
    return true; // a ≡ 0 mod n, not a valid witness
  }

  return millerRabinWitness(_n, _a);
}

// Small primes for trial division and Miller-Rabin
const SMALL_PRIMES = [
  2n,
  3n,
  5n,
  7n,
  11n,
  13n,
  17n,
  19n,
  23n,
  29n,
  31n,
  37n,
  41n,
  43n,
  47n,
  53n,
  59n,
  61n,
  67n,
  71n,
  73n,
  79n,
  83n,
  89n,
  97n,
  101n,
  103n,
  107n,
  109n,
  113n,
  127n,
  131n,
  137n,
  139n,
  149n,
  151n,
  157n,
  163n,
  167n,
  173n,
  179n,
  181n,
  191n,
  193n,
  197n,
  199n,
  211n,
  223n,
  227n,
  229n,
  233n,
  239n,
  241n,
  251n,
  257n,
  263n,
  269n,
  271n,
  277n,
  281n,
  283n,
  293n,
];

/**
 * Deterministic Miller-Rabin witnesses for numbers up to given bounds.
 * These witness sets are proven to give correct results.
 */
const MILLER_RABIN_WITNESSES: Array<[bigint, bigint[]]> = [
  [2047n, [2n]],
  [1373653n, [2n, 3n]],
  [9080191n, [31n, 73n]],
  [25326001n, [2n, 3n, 5n]],
  [3215031751n, [2n, 3n, 5n, 7n]],
  [4759123141n, [2n, 7n, 61n]],
  [1122004669633n, [2n, 13n, 23n, 1662803n]],
  [2152302898747n, [2n, 3n, 5n, 7n, 11n]],
  [3474749660383n, [2n, 3n, 5n, 7n, 11n, 13n]],
  [341550071728321n, [2n, 3n, 5n, 7n, 11n, 13n, 17n]],
  [3825123056546413051n, [2n, 3n, 5n, 7n, 11n, 13n, 17n, 19n, 23n]],
];

/**
 * Test whether n is a prime integer.
 *
 * Uses PARI's BPSW (Baillie-PSW) primality test via @sagemath-ts/parigp-ts.
 * BPSW combines Miller-Rabin base 2 with a Strong Lucas test.
 * No known composite passes BPSW as of 2024.
 *
 * @param n - Integer to test
 * @returns true if n is prime, false otherwise
 *
 * @example
 * ```typescript
 * is_prime(2n)    // true
 * is_prime(4n)    // false
 * is_prime(17n)   // true
 * is_prime(-5n)   // false (negative numbers are not prime)
 * ```
 *
 * @see Reference: sage/arith/misc.py:is_prime
 * @see Implementation: parigp-ts/src/ifactor.ts (port of pari/src/basemath/arith2.c)
 */
export function is_prime(n: IntegerLike): boolean {
  const _n = toBigInt(n);
  // Negative numbers are not prime (matches SageMath)
  if (_n <= 1n) {
    return false;
  }

  // Delegate to PARI's BPSW implementation
  return pari_isPrime(_n);
}

/**
 * Test whether n is a prime power `p^k` (k >= 1).
 *
 * SageMath (`sage/rings/integer.pyx:5399`) delegates to PARI's
 * ``isprimepower``, which extracts perfect powers and then runs a BPSW
 * primality test on the base -- it never factors n. The complete extraction
 * delegates to the existing `parigp-ts` implementation of ``isprimepower_i``.
 *
 * @param n - Integer to test
 * @param get_data - If true, return `[p, k]` with `n = p^k`, or `[n, 0]` when n
 *   is not a prime power
 * @returns Whether n is a prime power
 *
 * @example
 * ```typescript
 * is_prime_power(4n, true)    // [2n, 2n]
 * is_prime_power(512n, true)  // [2n, 9n]
 * is_prime_power(6n)          // false
 * ```
 *
 * @see Reference: sage/arith/misc.py:is_prime_power -> sage/rings/integer.pyx:is_prime_power
 */
export function is_prime_power(n: IntegerLike, get_data?: false): boolean;
export function is_prime_power(n: IntegerLike, get_data: true): [bigint, bigint];
export function is_prime_power(
  n: IntegerLike,
  get_data: boolean = false
): boolean | [bigint, bigint] {
  const _n = toBigInt(n);
  const data = pari_isprimepower(_n);
  if (data === null) {
    // SageMath returns ``(self, 0)`` when n is not a prime power.
    return get_data ? [_n, 0n] : false;
  }
  return get_data ? [data[0], BigInt(data[1])] : true;
}

/**
 * Return the next prime greater than n.
 *
 * @param n - Starting integer
 * @returns The smallest prime > n
 *
 * @example
 * ```typescript
 * next_prime(10n)  // 11n
 * next_prime(11n)  // 13n
 * ```
 *
 * @see Reference: sage/arith/misc.py:next_prime
 */
export function next_prime(n: IntegerLike): bigint {
  return pari_nextprime(toBigInt(n) + 1n);
}

/**
 * Return the largest prime less than n.
 *
 * @param n - Starting integer
 * @returns The largest prime < n
 * @throws {ValueError} If n <= 2
 *
 * @example
 * ```typescript
 * previous_prime(10n)  // 7n
 * previous_prime(11n)  // 7n
 * ```
 *
 * @see Reference: sage/arith/misc.py:previous_prime
 */
export function previous_prime(n: IntegerLike): bigint {
  let candidate = toBigInt(n) - 1n;
  if (candidate <= 1n) throw new ValueError('no previous prime');
  if (candidate <= 3n) return candidate;
  if ((candidate & 1n) === 0n) candidate--;
  while (!is_prime(candidate)) candidate -= 2n;
  return candidate;
}

/**
 * Type representing a prime factorization.
 * Array of [prime, exponent] pairs, sorted by prime.
 */
export type Factorization = Array<[bigint, bigint]>;

/**
 * Return the prime factorization of n.
 *
 * Uses PARI's factorization algorithm via @sagemath-ts/parigp-ts.
 *
 * @param n - Integer to factor (must be nonzero)
 * @returns Array of [prime, exponent] pairs
 *
 * @example
 * ```typescript
 * factor(12n)   // [[2n, 2n], [3n, 1n]]
 * factor(-12n)  // [[-1n, 1n], [2n, 2n], [3n, 1n]]
 * ```
 *
 * @see Reference: sage/arith/misc.py:factor
 * @see Implementation: parigp-ts/src/ifactor.ts (port of pari/src/basemath/ifactor1.c)
 * @see Deviation: PARI Integer Factorization (parigp-ts)
 */
export function factor(n: IntegerLike): Factorization {
  const _n = toBigInt(n);
  if (_n === 0n) {
    // SageMath raises ArithmeticError here (`sage/arith/misc.py` factor).
    throw new ArithmeticError('factorization of 0 is not defined');
  }

  // Delegate to PARI's Z_factor
  // This matches SageMath's behavior which also delegates to PARI
  return pari_Z_factor(_n);
}

/**
 * Format a factorization as a string matching SageMath's output.
 *
 * @param f - Factorization to format
 * @returns String representation like "2^2 * 3"
 */
export function formatFactorization(f: Factorization): string {
  if (f.length === 0) {
    return '1';
  }

  return f.map(([prime, exp]) => (exp === 1n ? `${prime}` : `${prime}^${exp}`)).join(' * ');
}

/**
 * Return Euler's totient function phi(n).
 *
 * @param n - Positive integer
 * @returns The number of integers k with 1 <= k <= n and gcd(k,n) = 1
 *
 * @example
 * ```typescript
 * euler_phi(12n)  // 4n
 * euler_phi(1n)   // 1n
 * ```
 *
 * @see Reference: sage/arith/misc.py:Euler_Phi
 */
export function euler_phi(n: IntegerLike): bigint {
  const value = toBigInt(n);
  if (value <= 0n) return 0n;
  if (value <= 2n) return 1n;
  return pari_eulerphi(value);
}

/**
 * Return the radical of n (product of distinct prime factors).
 *
 * @param n - Integer
 * @returns Product of distinct prime divisors of n
 *
 * @example
 * ```typescript
 * radical(12n)  // 6n (= 2 * 3)
 * ```
 */
export function radical(n: IntegerLike): bigint {
  const _n = toBigInt(n);
  // SageMath returns 0 for radical(0)
  if (_n === 0n) {
    return 0n;
  }

  const factors = factor(_n < 0n ? -_n : _n);
  let result = 1n;

  for (const [p] of factors) {
    if (p > 0n) {
      result *= p;
    }
  }

  return result;
}

/**
 * Return the Kronecker symbol (a/n).
 *
 * @param a - Integer
 * @param n - Integer
 * @returns The Kronecker symbol value (-1, 0, or 1)
 *
 * @see Reference: sage/arith/misc.py:kronecker_symbol
 */
export function kronecker_symbol(a: IntegerLike, n: IntegerLike): bigint {
  let _a = toBigInt(a);
  let _n = toBigInt(n);
  // Handle special cases
  if (_n === 0n) {
    return _a === 1n || _a === -1n ? 1n : 0n;
  }

  if (_n === 1n) {
    return 1n;
  }

  // Handle negative n
  let result = 1n;
  if (_n < 0n) {
    _n = -_n;
    if (_a < 0n) {
      result = -1n;
    }
  }

  // Remove factors of 2 from n
  let v = 0n;
  while ((_n & 1n) === 0n) {
    _n >>= 1n;
    v++;
  }

  if (v > 0n) {
    // (a|2) = 0 for even a, so (a|n) = 0 as soon as 2 | n and 2 | a.
    // GMP's mpz_kronecker does the same; without this the symbol would be
    // nonzero although gcd(a, n) > 1.
    if ((_a & 1n) === 0n) {
      return 0n;
    }
    // (a|2) = +1 if a = +-1 mod 8, -1 if a = +-3 mod 8.
    const aMod8 = ((_a % 8n) + 8n) % 8n;
    if ((v & 1n) === 1n) {
      if (aMod8 === 3n || aMod8 === 5n) {
        result = -result;
      }
    }
  }

  // Now n is odd
  _a = ((_a % _n) + _n) % _n;

  // Use quadratic reciprocity
  while (_a !== 0n) {
    // Remove factors of 2 from a
    let u = 0n;
    while ((_a & 1n) === 0n) {
      _a >>= 1n;
      u++;
    }

    if ((u & 1n) === 1n) {
      const nMod8 = _n & 7n;
      if (nMod8 === 3n || nMod8 === 5n) {
        result = -result;
      }
    }

    // Quadratic reciprocity
    if ((_a & 3n) === 3n && (_n & 3n) === 3n) {
      result = -result;
    }

    const temp = _a;
    _a = _n % temp;
    _n = temp;
  }

  return _n === 1n ? result : 0n;
}

/**
 * Return the Legendre symbol (a/p), for p an odd prime.
 *
 * Unlike {@link kronecker_symbol}, this rejects composite and even moduli
 * exactly as SageMath does.
 *
 * @param a - Integer
 * @param p - Odd prime
 * @returns The Legendre symbol value (-1, 0, or 1)
 * @throws {ValueError} `p must be a prime` if p is not prime,
 *   `p must be odd` if p = 2
 *
 * @example
 * ```typescript
 * legendre_symbol(2n, 3n)   // -1n
 * legendre_symbol(2n, 15n)  // ValueError: p must be a prime
 * legendre_symbol(1n, 2n)   // ValueError: p must be odd
 * ```
 *
 * @see Reference: sage/arith/misc.py:legendre_symbol
 */
export function legendre_symbol(a: IntegerLike, p: IntegerLike): bigint {
  const _p = toBigInt(p);
  if (!is_prime(_p)) {
    throw new ValueError('p must be a prime');
  }
  if (_p === 2n) {
    throw new ValueError('p must be odd');
  }
  return kronecker_symbol(a, _p);
}

/**
 * Return the Jacobi symbol (a/b), for b odd.
 *
 * SageMath only requires b to be odd; negative odd b is accepted (the value is
 * then the Kronecker symbol, which carries the (a|-1) factor).
 *
 * @param a - Integer
 * @param b - Odd integer
 * @returns The Jacobi symbol value (-1, 0, or 1)
 * @throws {ValueError} `second input must be odd, <b> is not odd`
 *
 * @example
 * ```typescript
 * jacobi_symbol(10n, 777n)  // -1n
 * jacobi_symbol(10n, -3n)   // 1n
 * jacobi_symbol(10n, 2n)    // ValueError
 * ```
 *
 * @see Reference: sage/arith/misc.py:jacobi_symbol
 */
export function jacobi_symbol(a: IntegerLike, b: IntegerLike): bigint {
  const _b = toBigInt(b);
  if (_b % 2n === 0n) {
    throw new ValueError(`second input must be odd, ${_b} is not odd`);
  }
  return kronecker_symbol(a, _b);
}

/**
 * Return the Chinese Remainder Theorem solution.
 *
 * @param a - First residue
 * @param b - Second residue
 * @param m - First modulus
 * @param n - Second modulus
 * @returns x such that x ≡ a (mod m) and x ≡ b (mod n)
 *
 * @example
 * ```typescript
 * crt(2n, 3n, 3n, 5n)  // 8n
 * ```
 *
 * @see Reference: sage/arith/misc.py:crt
 */
export function crt(a: IntegerLike[], b: IntegerLike[]): bigint;
export function crt(a: IntegerLike, b: IntegerLike, m: IntegerLike, n: IntegerLike): bigint;
export function crt(
  a: IntegerLike | IntegerLike[],
  b: IntegerLike | IntegerLike[],
  m?: IntegerLike,
  n?: IntegerLike
): bigint {
  if (Array.isArray(a)) return CRT_list(a, b as IntegerLike[]);
  const _a = toBigInt(a);
  const _b = toBigInt(b as IntegerLike);
  const _m = toBigInt(m!);
  const _n = toBigInt(n!);
  const [g, s, _t] = xgcd(_m, _n);
  if (g === 0n) throw new ZeroDivisionError('Integer division by zero');

  if ((_a - _b) % g !== 0n) {
    // SageMath's message, verbatim (sage/arith/misc.py:3493).
    throw new ValueError(
      `no solution to crt problem since gcd(${_m},${_n}) does not divide ${_a}-${_b}`
    );
  }

  const lcmMN = lcm(_m, _n);
  if (lcmMN === 0n) throw new ZeroDivisionError('Integer modulo by zero');
  let result = (_a + _m * (((_b - _a) / g) * s)) % lcmMN;

  if (result < 0n) {
    result += lcmMN;
  }

  return result;
}

/** Integer.__mod__ (including signed divisors) at the CRT boundary. */
function crtRemainder(value: bigint, modulus: bigint): bigint {
  if (modulus === 0n) throw new ZeroDivisionError('Integer modulo by zero');
  const r = value % modulus;
  return r !== 0n && r < 0n !== modulus < 0n ? r + modulus : r;
}

type CRTModularValue = IntegerMod | PrimeFieldElement | LegacyPrimeElement;

/**
 * Combine residues through Sage's balanced binary tree of CRT operations.
 * A single modular-element argument preserves its parent and singleton identity.
 * @see Reference: sage/arith/misc.py:CRT_list
 */
export function CRT_list(residues: IntegerLike[], moduli: IntegerLike[]): bigint;
export function CRT_list<T extends CRTModularValue>(residues: T[]): T | IntegerMod;
export function CRT_list(
  residues: IntegerLike[] | CRTModularValue[],
  moduli?: IntegerLike[] | null
): bigint | CRTModularValue {
  if (!Array.isArray(residues) || (moduli != null && !Array.isArray(moduli))) {
    throw new ValueError('arguments to CRT_list should be lists');
  }
  const returnMod = moduli == null;
  let values: bigint[], mods: bigint[];
  if (returnMod) {
    if (residues.length === 0) return Mod(0n, 1n);
    if (
      !residues.every(
        (v) =>
          v instanceof IntegerMod ||
          v instanceof PrimeFieldElement ||
          v instanceof LegacyPrimeElement
      )
    ) {
      throw new TypeError('if one argument is given, it should be a list of IntegerMod');
    }
    const elements = residues as CRTModularValue[];
    if (elements.length === 1) return elements[0]!;
    mods = elements.map((v) => (v instanceof IntegerMod ? v.modulus : v.parent.characteristic));
    values = elements.map((v) => v.value);
  } else {
    if (residues.length !== moduli!.length) {
      throw new ValueError('arguments to CRT_list should be lists of the same length');
    }
    if (residues.length === 0) return 0n;
    if (residues.length === 1) return toBigInt(residues[0] as IntegerLike);
    values = (residues as IntegerLike[]).map(toBigInt);
    mods = moduli!.map(toBigInt);
  }
  while (values.length > 1) {
    const nextValues: bigint[] = [],
      nextMods: bigint[] = [];
    for (let i = 0; i < values.length; i += 2) {
      if (i + 1 === values.length) {
        nextValues.push(values[i]!);
        nextMods.push(mods[i]!);
      } else {
        nextValues.push(crt(values[i]!, values[i + 1]!, mods[i]!, mods[i + 1]!));
        nextMods.push(lcm(mods[i]!, mods[i + 1]!));
      }
    }
    values = nextValues;
    mods = nextMods;
  }
  return returnMod ? Mod(values[0]!, mods[0]!) : crtRemainder(values[0]!, mods[0]!);
}

/**
 * Test whether n is a perfect square.
 *
 * With `root = true` this mirrors SageMath's `(True, sqrt)` / `(False, None)`:
 * the second component is `null` (Python's `None`) when n is not a square, so
 * it can never be confused with the genuine root of 0.
 *
 * @param n - Integer to test
 * @param root - If true, return [isSquare, root] instead of just boolean
 * @returns Whether n is a perfect square (and optionally its root, or null)
 *
 * @example
 * ```typescript
 * is_square(16n, true)  // [true, 4n]
 * is_square(15n, true)  // [false, null]
 * is_square(0n, true)   // [true, 0n]
 * ```
 *
 * @see Reference: sage/arith/misc.py:is_square -> sage/rings/integer.pyx:is_square
 */
export function is_square(n: IntegerLike, root?: false): boolean;
export function is_square(n: IntegerLike, root: true): [boolean, bigint | null];
export function is_square(
  n: IntegerLike,
  root: boolean = false
): boolean | [boolean, bigint | null] {
  const _n = toBigInt(n);
  if (_n < 0n) {
    if (root) {
      return [false, null];
    }
    return false;
  }

  if (_n === 0n) {
    if (root) {
      return [true, 0n];
    }
    return true;
  }

  const s = isqrt(_n);
  const isSquareVal = s * s === _n;

  if (root) {
    return [isSquareVal, isSquareVal ? s : null];
  }

  return isSquareVal;
}

/**
 * Test whether n is squarefree (not divisible by any perfect square > 1).
 *
 * @param n - Integer to test
 * @returns true if n is squarefree
 */
export function is_squarefree(n: IntegerLike): boolean {
  let _n = toBigInt(n);
  if (_n === 0n) {
    return false;
  }

  _n = _n < 0n ? -_n : _n;

  if (_n === 1n) {
    return true;
  }

  const factors = factor(_n);

  for (const [_p, e] of factors) {
    if (e > 1n) {
      return false;
    }
  }

  return true;
}

/**
 * Return the list of all positive divisors of n.
 *
 * @param n - Positive integer
 * @returns Sorted list of positive divisors
 */
export function divisors(n: IntegerLike): bigint[] {
  let _n = toBigInt(n);
  if (_n === 0n) {
    throw new ValueError('divisors of 0 is not defined');
  }

  // SageMath returns divisors of |n| for negative numbers
  if (_n < 0n) {
    _n = -_n;
  }

  if (_n === 1n) {
    return [1n];
  }

  const factors = factor(_n);
  let divs = [1n];

  for (const [p, e] of factors) {
    if (p === -1n) continue;

    const newDivs: bigint[] = [];
    let pk = 1n;

    for (let k = 0n; k <= e; k++) {
      for (const d of divs) {
        newDivs.push(d * pk);
      }
      pk *= p;
    }

    divs = newDivs;
  }

  return divs.sort((a, b) => (a < b ? -1 : a > b ? 1 : 0));
}

/**
 * Return the number of divisors of the (nonzero) integer n.
 *
 * SageMath delegates to PARI's ``numdiv``, which ignores the sign of n; only
 * n = 0 is rejected.
 *
 * @param n - Nonzero integer
 * @returns The number of positive divisors of |n|
 * @throws {ValueError} `input must be nonzero` when n = 0
 *
 * @example
 * ```typescript
 * number_of_divisors(100n)   // 9n
 * number_of_divisors(-720n)  // 30n
 * ```
 *
 * @see Reference: sage/arith/misc.py:number_of_divisors
 */
export function number_of_divisors(n: IntegerLike): bigint {
  const value = toBigInt(n);
  if (value === 0n) throw new ValueError('input must be nonzero');
  return pari_numdiv(value);
}

/**
 * Return the sum of the k-th powers of the divisors of n.
 *
 * SageMath's `Sigma.__call__` (`sage/arith/misc.py:1653`) evaluates the
 * multiplicative closed form over `factor(n)` -- it never enumerates the
 * divisors, so `sigma(factorial(100), 0)` (about 3.9e16 divisors) is instant.
 * The unit -1 of a negative factorization is not part of the underlying list,
 * hence `sigma(-4) = sigma(4) = 7`.
 *
 * @param n - Nonzero integer
 * @param k - Power (default 1)
 * @returns Sum of d^k over the positive divisors d of |n|
 *
 * @example
 * ```typescript
 * sigma(6n)       // 12n
 * sigma(6n, 2n)   // 50n
 * sigma(100n, 4n) // 106811523n
 * sigma(-4n)      // 7n
 * ```
 *
 * @see Reference: sage/arith/misc.py:Sigma
 */
export function sigma(n: IntegerLike, k: IntegerLike = 1n): bigint {
  const _n = toBigInt(n);
  const _k = toBigInt(k);

  // factor() raises for 0, exactly as SageMath's factor(0) does.
  const factors = factor(_n);

  let result = 1n;
  for (const [p, e] of factors) {
    // The unit (-1) is not part of a Factorization's underlying list.
    if (p === -1n) continue;
    // Sage Sigma.__call__ invokes an Integer-only method on a Rational for
    // negative k. Preserve that observable failure, including the unit cases
    // (where this loop is empty). See Deviation: Integer Audit Oracle Boundaries.
    if (_k < 0n) {
      throw new AttributeError(
        "'sage.rings.rational.Rational' object has no attribute 'divide_knowing_divisible_by'"
      );
    }
    if (_k === 0n) {
      result *= e + 1n;
    } else if (_k === 1n) {
      result *= (p ** (e + 1n) - 1n) / (p - 1n);
    } else {
      result *= (p ** ((e + 1n) * _k) - 1n) / (p ** _k - 1n);
    }
  }

  return result;
}

/**
 * Find a quadratic non-residue modulo p.
 *
 * @param p - An odd prime
 * @returns A quadratic non-residue modulo p
 */
function find_quadratic_nonresidue(p: bigint): bigint {
  // For p ≡ 3 (mod 4), -1 is a non-residue
  if (p % 4n === 3n) {
    return p - 1n;
  }
  // For p ≡ 5 (mod 8), 2 is a non-residue
  if (p % 8n === 5n) {
    return 2n;
  }
  // Otherwise, search for a non-residue
  for (let n = 2n; n < p; n++) {
    if (legendre_symbol(n, p) === -1n) {
      return n;
    }
  }
  // Should never reach here for a prime p
  throw new ValueError('Could not find quadratic non-residue');
}

/**
 * Compute a square root of a modulo p using the Tonelli-Shanks algorithm.
 *
 * Returns a square root of a modulo p, or null if a is not a quadratic residue.
 * If all_roots is true, returns an array of all square roots (0, 1, or 2 values).
 *
 * @param a - The value to find the square root of
 * @param p - A prime modulus
 * @param all_roots - If true, return all square roots; otherwise return one or null
 * @returns A square root, array of square roots, or null
 *
 * @example
 * ```typescript
 * sqrt_mod(2n, 7n)           // 3n or 4n (since 3^2 = 9 ≡ 2 mod 7)
 * sqrt_mod(2n, 7n, true)     // [3n, 4n]
 * sqrt_mod(3n, 7n)           // null (3 is not a QR mod 7)
 * sqrt_mod(3n, 7n, true)     // []
 * sqrt_mod(0n, 7n)           // 0n
 * sqrt_mod(0n, 7n, true)     // [0n]
 * ```
 *
 * @see Reference: sage/rings/finite_rings/integer_mod.pyx:square_root_mod_prime
 */
export function sqrt_mod(a: IntegerLike, p: IntegerLike, all_roots?: false): bigint | null;
export function sqrt_mod(a: IntegerLike, p: IntegerLike, all_roots: true): bigint[];
export function sqrt_mod(
  a: IntegerLike,
  p: IntegerLike,
  all_roots: boolean = false
): bigint | bigint[] | null {
  let _a = toBigInt(a);
  const _p = toBigInt(p);
  // Normalize a to be in [0, p)
  _a = ((_a % _p) + _p) % _p;

  // Handle edge case: a = 0
  if (_a === 0n) {
    return all_roots ? [0n] : 0n;
  }

  // Handle edge case: p = 2
  if (_p === 2n) {
    // In Z/2Z, every element is its own square root
    return all_roots ? [_a] : _a;
  }

  // Check if a is a quadratic residue using Legendre symbol
  const ls = legendre_symbol(_a, _p);
  if (ls === 0n) {
    // a ≡ 0 (mod p), but we already handled a = 0 above
    return all_roots ? [0n] : 0n;
  }
  if (ls === -1n) {
    // a is not a quadratic residue
    return all_roots ? [] : null;
  }

  // Now we know a is a quadratic residue (ls === 1)
  let root: bigint;

  // Case: p ≡ 3 (mod 4)
  // sqrt(a) = a^((p+1)/4)
  if (_p % 4n === 3n) {
    root = power_mod(_a, (_p + 1n) / 4n, _p);
  }
  // Case: p ≡ 5 (mod 8)
  // Use the formula with i = sqrt(-1)
  else if (_p % 8n === 5n) {
    const two_a = (2n * _a) % _p;
    const zeta = power_mod(two_a, (_p - 5n) / 8n, _p);
    const i = (((zeta * zeta) % _p) * two_a) % _p; // = (2a)^((p-1)/4)
    root = (((zeta * _a) % _p) * ((i - 1n + _p) % _p)) % _p;
  }
  // General case: Tonelli-Shanks algorithm
  else {
    // Write p - 1 = 2^r * q where q is odd
    let q = _p - 1n;
    let r = 0n;
    while ((q & 1n) === 0n) {
      q >>= 1n;
      r++;
    }

    // Find a quadratic non-residue n
    const nqr = find_quadratic_nonresidue(_p);

    // v = n^q mod p (a primitive 2^r-th root of unity)
    let v = power_mod(nqr, q, _p);

    // x = a^((q-1)/2) mod p
    const x = power_mod(_a, (q - 1n) / 2n, _p);

    // b = a * x^2 = a^q mod p
    let b = (((_a * x) % _p) * x) % _p;

    // res = a * x = a^((q+1)/2) mod p
    let res = (_a * x) % _p;

    // Main loop: adjust res until b = 1
    while (b !== 1n) {
      // Find the smallest m such that b^(2^m) = 1
      let m = 1n;
      let bpow = (b * b) % _p;
      while (bpow !== 1n) {
        bpow = (bpow * bpow) % _p;
        m++;
      }

      // g = v^(2^(r-m-1))
      const exp = 1n << (r - m - 1n);
      const g = power_mod(v, exp, _p);

      // Update values
      res = (res * g) % _p;
      v = (g * g) % _p;
      b = (b * v) % _p;
      r = m;
    }

    root = res;
  }

  // Ensure we return the canonical square root (smaller of the two)
  const negRoot = (_p - root) % _p;
  if (negRoot < root) {
    root = negRoot;
  }

  if (all_roots) {
    if (root === 0n) {
      return [0n];
    }
    const roots = [root, (_p - root) % _p];
    roots.sort((x, y) => (x < y ? -1 : x > y ? 1 : 0));
    return roots;
  }

  return root;
}

/**
 * Generate all primes in a range.
 *
 * If only one argument is given, returns all primes less than that value.
 * If two arguments are given, returns all primes in [start, stop).
 *
 * @param start - If stop is undefined, upper bound (exclusive). Otherwise, lower bound (inclusive).
 * @param stop - Upper bound (exclusive)
 * @returns Array of prime numbers as bigints
 *
 * @example
 * ```typescript
 * prime_range(10n)        // [2n, 3n, 5n, 7n]
 * prime_range(10n, 20n)   // [11n, 13n, 17n, 19n]
 * prime_range(2n, 3n)     // [2n]
 * ```
 *
 * @see Reference: sage/rings/fast_arith.pyx:prime_range
 */
export function prime_range(start: IntegerLike, stop?: IntegerLike): bigint[] {
  let _start = toBigInt(start);
  let _stop = stop !== undefined ? toBigInt(stop) : undefined;
  // Handle single argument case: primes less than start
  if (_stop === undefined) {
    _stop = _start;
    _start = 2n;
  }

  if (_stop <= _start || _stop <= 2n) {
    return [];
  }

  if (_start < 2n) {
    _start = 2n;
  }

  // For small ranges, use Sieve of Eratosthenes
  const SIEVE_LIMIT = 1000000n; // Use sieve for ranges up to 1 million

  if (_stop <= SIEVE_LIMIT) {
    return sieveOfEratosthenes(_start, _stop);
  }

  // For larger ranges, iterate with is_prime checks
  const primes: bigint[] = [];

  // Handle 2 separately
  if (_start <= 2n && _stop > 2n) {
    primes.push(2n);
  }

  // Start from next odd number >= start
  let candidate = _start;
  if (candidate <= 2n) {
    candidate = 3n;
  } else if ((candidate & 1n) === 0n) {
    candidate++;
  }

  while (candidate < _stop) {
    if (is_prime(candidate)) {
      primes.push(candidate);
    }
    candidate += 2n;
  }

  return primes;
}

/**
 * Sieve of Eratosthenes for generating primes in a range.
 */
function sieveOfEratosthenes(start: bigint, stop: bigint): bigint[] {
  const n = toSafeNumber(stop);
  const startNum = toSafeNumber(start);

  // Create sieve array (index i represents number i)
  const sieve = new Uint8Array(n);

  // 0 and 1 are not prime
  sieve[0] = 1;
  sieve[1] = 1;

  // Mark composites
  for (let i = 2; i * i < n; i++) {
    if (sieve[i] === 0) {
      for (let j = i * i; j < n; j += i) {
        sieve[j] = 1;
      }
    }
  }

  // Collect primes in range
  const primes: bigint[] = [];
  for (let i = Math.max(2, startNum); i < n; i++) {
    if (sieve[i] === 0) {
      primes.push(BigInt(i));
    }
  }

  return primes;
}

/**
 * Return the value of the Moebius function mu(n).
 *
 * mu(n) is defined as:
 * - mu(1) = 1
 * - mu(n) = 0 if n has a squared prime factor
 * - mu(n) = (-1)^k if n is a product of k distinct primes
 *
 * For simplicity, mu(0) = 0.
 *
 * @param n - Integer (absolute value is used)
 * @returns 0, 1, or -1
 *
 * @example
 * ```typescript
 * moebius(1n)   // 1n
 * moebius(6n)   // 1n (6 = 2 * 3, two distinct primes)
 * moebius(5n)   // -1n (5 is prime)
 * moebius(4n)   // 0n (4 = 2^2, has squared factor)
 * moebius(30n)  // -1n (30 = 2 * 3 * 5, three distinct primes)
 * ```
 *
 * @see Reference: sage/arith/misc.py:Moebius
 */
export function moebius(n: IntegerLike): bigint {
  let _n = toBigInt(n);
  // Handle special cases
  if (_n === 0n) {
    return 0n;
  }

  // Work with absolute value
  if (_n < 0n) {
    _n = -_n;
  }

  if (_n === 1n) {
    return 1n;
  }

  // Factor n and check for squared factors
  const factors = factor(_n);
  let numPrimes = 0n;

  for (const [p, e] of factors) {
    if (p === -1n) continue;
    if (e >= 2n) {
      return 0n; // n has a squared prime factor
    }
    numPrimes++;
  }

  // mu(n) = (-1)^k where k is the number of distinct prime factors
  return numPrimes % 2n === 0n ? 1n : -1n;
}

/**
 * Return the squarefree part of n.
 *
 * The squarefree part of n is the unique integer z such that n = z * y^2
 * where y^2 is a perfect square and z is squarefree.
 *
 * @param n - Integer
 * @returns The squarefree part of n
 *
 * @example
 * ```typescript
 * squarefree_part(12n)   // 3n (12 = 3 * 4 = 3 * 2^2)
 * squarefree_part(72n)   // 2n (72 = 2 * 36 = 2 * 6^2)
 * squarefree_part(100n)  // 1n (100 = 1 * 10^2)
 * squarefree_part(-12n)  // -3n
 * ```
 *
 * @see Reference: sage/rings/integer.pyx:squarefree_part
 */
export function squarefree_part(n: IntegerLike): bigint {
  let _n = toBigInt(n);
  if (_n === 0n) {
    return 0n;
  }

  const sign = _n < 0n ? -1n : 1n;
  _n = _n < 0n ? -_n : _n;

  if (_n === 1n) {
    return sign;
  }

  const factors = factor(_n);
  let result = 1n;

  for (const [p, e] of factors) {
    if (p === -1n) continue;
    // Include p if its exponent is odd
    if (e % 2n === 1n) {
      result *= p;
    }
  }

  return sign * result;
}

/**
 * Return the list of prime factors of n (without multiplicities).
 *
 * @param n - Integer to factor
 * @returns Sorted list of distinct prime divisors
 *
 * @example
 * ```typescript
 * prime_factors(12n)   // [2n, 3n]
 * prime_factors(100n)  // [2n, 5n]
 * prime_factors(1n)    // []
 * prime_factors(-12n)  // [2n, 3n]
 * ```
 *
 * @see Reference: sage/arith/misc.py:prime_divisors
 */
export function prime_factors(n: IntegerLike): bigint[] {
  let _n = toBigInt(n);
  if (_n === 0n) {
    throw new ArithmeticError('factorization of 0 is not defined');
  }

  if (_n < 0n) {
    _n = -_n;
  }

  if (_n === 1n) {
    return [];
  }

  const factors = factor(_n);
  const primes: bigint[] = [];

  for (const [p] of factors) {
    if (p > 0n) {
      primes.push(p);
    }
  }

  return primes;
}

/**
 * Alias for prime_factors.
 *
 * @see prime_factors
 */
export const prime_divisors = prime_factors;

/**
 * Return the largest k such that p^k divides n.
 *
 * @param n - Integer
 * @param p - Prime (or any integer > 1)
 * @returns The p-adic valuation of n
 *
 * @example
 * ```typescript
 * valuation(24n, 2n)   // 3n (24 = 2^3 * 3)
 * valuation(100n, 5n)  // 2n (100 = 4 * 5^2)
 * valuation(7n, 2n)    // 0n (7 is odd)
 * valuation(0n, 2n)    // 'Infinity'
 * ```
 *
 * @see Reference: sage/arith/misc.py:valuation
 * @see Deviation: Valuation dispatch and native GMP factor removal
 */
export function valuation<T>(n: { valuation(p: IntegerLike): T }, p: IntegerLike): T;
export function valuation(
  n: IntegerLike | { _integer_: (parent: typeof ZZ) => IntegerLike },
  p: IntegerLike
): bigint | 'Infinity';
export function valuation(
  n: IntegerLike | { _integer_: (parent: typeof ZZ) => IntegerLike } | { valuation(p: IntegerLike): unknown },
  p: IntegerLike
): unknown {
  if (typeof n === 'bigint') return new Integer(n).valuation(p);
  try {
    // Missing JavaScript properties do not raise Python's AttributeError.
    if ('valuation' in n) return n.valuation(p);
  } catch (error) {
    // The source also falls back if the method raises AttributeError internally.
    if (!(error instanceof AttributeError)) throw error;
  }
  return new Integer(n as IntegerLike).valuation(p);
}

// ============================================================================
// STUB FUNCTIONS - Not yet implemented
// ============================================================================

/** Precision and proof options for Sage's real-field relation lattice. */
type AlgebraicDependencyOptions = {
  known_bits?: bigint;
  use_bits?: bigint;
  known_digits?: bigint;
  use_digits?: bigint;
  height_bound?: bigint;
  proof?: boolean;
};

/**
 * Return ascending coefficients of the best-fitting irreducible relation.
 * Python-float/JavaScript-number inputs delegate to PARI. Real-field inputs
 * use Sage's precision-controlled LLL lattice; exact inputs return a linear
 * relation. A height bound can exclude the relation, returning null.
 * @see Reference: sage/arith/misc.py:algebraic_dependency
 * @see Deviation: PARI algebraic dependencies
 */
export function algebraic_dependency(
  z: number | IntegerLike | Rational | RealNumber,
  degree: IntegerLike
): bigint[];
export function algebraic_dependency(
  z: number,
  degree: IntegerLike,
  options: AlgebraicDependencyOptions
): bigint[];
export function algebraic_dependency(
  z: number | IntegerLike | Rational | RealNumber,
  degree: IntegerLike,
  options?: AlgebraicDependencyOptions
): bigint[] | null;
export function algebraic_dependency(
  z: number | IntegerLike | Rational | RealNumber,
  degree: IntegerLike,
  options?: AlgebraicDependencyOptions
): bigint[] | null {
  const opts = options ?? {};
  const height = opts.height_bound;
  if (opts.proof && !height) throw new ValueError('height_bound must be given for proof=True');
  const abs = (x: bigint) => (x < 0n ? -x : x);
  // Sage returns an integer relation before coercing the requested degree.
  if (typeof z === 'bigint' || z instanceof Integer) {
    const value = toBigInt(z);
    return height && abs(value) >= height ? null : [-value, 1n];
  }
  const n = toBigInt(degree);
  if (z instanceof Rational) {
    return height && (abs(z.numerator) >= height || z.denominator >= height)
      ? null
      : [-z.numerator, z.denominator];
  }
  let coefficients: bigint[];
  let real: mpfr_t | undefined;
  if (z instanceof RealNumber) {
    let bits = z.precision() - 6;
    let known = opts.known_bits === undefined ? undefined : Number(opts.known_bits);
    if (opts.known_digits !== undefined) known = Number(opts.known_digits) * Math.log2(10);
    let use = opts.use_bits === undefined ? undefined : Number(opts.use_bits);
    if (known !== undefined) use = known * 0.8;
    if (opts.use_digits !== undefined) use = Number(opts.use_digits) * Math.log2(10);
    if (use !== undefined) bits = Math.trunc(use);
    if (n < -1n) throw new ArithmeticError('number of rows must be non-negative');
    if (n === -1n) throw new IndexError('index out of range');
    const size = toSafeNumber(n + 1n);
    const rows: bigint[][] = Array.from({ length: size }, (_, i) =>
      Array.from({ length: size + 1 }, (_, j) => (i === j ? 1n : 0n))
    );
    const scale = 1n << BigInt(bits);
    rows[0]![size] = scale;
    const [sign, mantissa, exponent] = z.sign_mantissa_exponent();
    real = {
      kind: z.is_NaN() ? 'nan' : z.is_infinity() ? 'inf' : mantissa ? 'finite' : 'zero',
      precision: z.precision(),
      sign: sign < 0 ? -1 : 1,
      mantissa,
      exponent: mantissa ? Number(exponent) + z.precision() : 0,
    };
    const power = mpfr_init2(z.precision()),
      rounded = mpfr_init2(z.precision());
    mpfr_set_z(power, scale);
    for (let k = 1; k < size; k++) {
      mpfr_mul(power, power, real);
      if (power.kind === 'nan' || power.kind === 'inf')
        throw new ValueError('cannot convert infinity or NaN to Sage Integer');
      mpfr_round(rounded, power);
      rows[k]![size] = mpfr_get_z(rounded, 'RNDZ')[0];
    }
    const reduced = LLL(new IntegerMatrix(size, size + 1, rows), 0.75) as IntegerMatrix;
    const row = (i: number) => Array.from({ length: size + 1 }, (_, j) => reduced.get(i, j).value);
    coefficients = row(0).slice(0, size);
    if (coefficients.slice(1).every((c) => c === 0n)) {
      if (size === 1) throw new IndexError('matrix index out of range');
      coefficients = row(1).slice(0, size);
    }
    if (height) {
      const norm2 = (i: number) => row(i).reduce((sum, c) => sum + c * c, 0n);
      const factor = 1n << n;
      const bound2 = BigInt(size) * height * height;
      if (coefficients.some((c) => abs(c) > height)) {
        if (opts.proof && height > 0n && norm2(0) <= factor * bound2)
          throw new ValueError('insufficient precision for non-existence proof');
        return null;
      }
      // max(RIF(sqrt(norm2)), sqrt(n)*height) compares an interval with
      // a symbolic irrational. At equality Sage's symbolic/Maxima conversion
      // of that interval raises TypeError rather than certifying the relation.
      if (opts.proof && norm2(0) === bound2 && isqrt(BigInt(size)) ** 2n !== BigInt(size))
        throw new TypeError('');
      if (opts.proof && norm2(1) < factor * (norm2(0) > bound2 ? norm2(0) : bound2))
        throw new ValueError('insufficient precision for uniqueness proof');
    }
    if (coefficients[size - 1]! < 0n) coefficients = coefficients.map((c) => -c);
  } else {
    if (opts.proof || height)
      throw new NotImplementedError(
        'proof and height bound only implemented for real and complex numbers'
      );
    coefficients = pari_algdep(z, n);
  }
  // Polynomial coefficients use Integer objects; the public ZZ factory returns bigint.
  const integerRing = {
    zero: () => new Integer(0n),
    one: () => new Integer(1n),
    __call__: (value: unknown) =>
      new Integer(ZZ.__call__(value as ConstructorParameters<typeof Integer>[0])),
    is_field: () => false,
    is_integral_domain: () => true,
    characteristic: () => 0n,
    toString: () => 'Integer Ring',
  };
  // Integer intentionally restricts eq to IntegerLike, while the generic
  // polynomial interface also permits number. This ring supplies only Integers.
  const ring = new PolynomialRing<Integer & RingElement>(
    integerRing as CoefficientRing<Integer & RingElement>,
    'x'
  );
  const polynomial = ring.__call__(
    coefficients.map((c) => integerRing.__call__(c) as Integer & RingElement)
  );
  const factors = polynomial.factor().filter(([factor]) => factor.degree() > 0);
  if (factors.length === 0) throw new ValueError('min() arg is an empty sequence');
  let best: bigint[] | undefined;
  let bestValue = Infinity;
  let bestReal: mpfr_t | undefined;
  for (const [factor] of factors) {
    const coeffs = factor.coeffs.map((c) => c.value);
    if (real) {
      const value = mpfr_init2(real.precision),
        coefficient = mpfr_init2(real.precision);
      mpfr_set_z(value, 0n);
      for (let i = coeffs.length - 1; i >= 0; i--) {
        mpfr_mul(value, value, real);
        mpfr_set_z(coefficient, coeffs[i]!);
        mpfr_add(value, value, coefficient);
      }
      value.sign = 1;
      if (bestReal === undefined || mpfr_cmp(value, bestReal) < 0) {
        best = coeffs;
        bestReal = value;
      }
    } else {
      let value = 0;
      for (let i = coeffs.length - 1; i >= 0; i--)
        value = value * (z as number) + Number(coeffs[i]!);
      const residual = Math.abs(value);
      if (best === undefined || residual < bestValue) {
        best = coeffs;
        bestValue = residual;
      }
    }
  }
  return best!;
}

/**
 * Alias for algebraic_dependency.
 * @see algebraic_dependency
 */
export const algdep = algebraic_dependency;

/**
 * Return the n-th Bernoulli number, as a rational number.
 *
 * The Bernoulli numbers B_n are defined by the generating function:
 * x / (e^x - 1) = sum_{n=0}^{infty} B_n x^n / n!
 *
 * @param n - Integer
 * @param algorithm - Algorithm to use (default: 'default')
 * @param num_threads - Number of threads for parallel algorithms
 * @returns The n-th Bernoulli number
 *
 * @example
 * ```typescript
 * bernoulli(12n)  // -691n/2730n as a rational
 * ```
 *
 * @see Reference: sage/arith/misc.py:bernoulli
 * @see Deviation: Bernoulli Numbers (single algorithm, size limits)
 */
export function bernoulli(
  n: bigint,
  algorithm?: 'default' | 'arb' | 'flint' | 'pari' | 'gap' | 'gp' | 'bernmm',
  num_threads?: number
): { numerator: bigint; denominator: bigint } {
  // Reference: sage/arith/misc.py:bernoulli
  // SageMath uses PARI's bernfrac() for the 'pari' algorithm
  // For TypeScript, we implement the algorithm directly using the formula:
  // B_n = sum_{k=0}^{n} sum_{v=0}^{k} (-1)^v * C(k,v) * v^n / (k+1)
  // with the Von Staudt-Clausen theorem for the denominator

  if (n < 0n) {
    throw new ValueError('Bernoulli number index must be non-negative');
  }

  // B_0 = 1
  if (n === 0n) {
    return { numerator: 1n, denominator: 1n };
  }

  // B_1 = -1/2
  if (n === 1n) {
    return { numerator: -1n, denominator: 2n };
  }

  // B_n = 0 for odd n > 1
  if (n % 2n === 1n) {
    return { numerator: 0n, denominator: 1n };
  }

  // For even n >= 2, compute using the Akiyama-Tanigawa algorithm
  // This is more efficient for smaller values
  const nNum = toSafeNumber(n);

  // Use a table-based approach for Bernoulli numbers
  // B_n = sum_{k=0}^{n} 1/(k+1) * sum_{j=0}^{k} (-1)^j * C(k,j) * j^n
  // We use the recurrence relation: B_n = -1/(n+1) * sum_{k=0}^{n-1} C(n+1,k) * B_k

  // Store Bernoulli numbers as [numerator, denominator] pairs
  const B: Array<[bigint, bigint]> = [
    [1n, 1n],
    [-1n, 2n],
  ];

  for (let m = 2n; m <= n; m++) {
    if (m % 2n === 1n) {
      B.push([0n, 1n]);
      continue;
    }

    // B_m = -1/(m+1) * sum_{k=0}^{m-1} C(m+1,k) * B_k
    let sumNum = 0n;
    let sumDen = 1n;

    for (let k = 0n; k < m; k++) {
      const bk = B[toSafeNumber(k)]!;
      const coeff = binomial(m + 1n, k);

      // Add coeff * B_k to sum
      const termNum = coeff * bk[0];
      const termDen = bk[1];

      // sumNum/sumDen + termNum/termDen
      sumNum = sumNum * termDen + termNum * sumDen;
      sumDen = sumDen * termDen;

      // Simplify
      const g = gcd(sumNum < 0n ? -sumNum : sumNum, sumDen);
      sumNum = sumNum / g;
      sumDen = sumDen / g;
    }

    // B_m = -sumNum / ((m+1) * sumDen)
    let bmNum = -sumNum;
    let bmDen = (m + 1n) * sumDen;

    // Simplify
    const g = gcd(bmNum < 0n ? -bmNum : bmNum, bmDen);
    bmNum = bmNum / g;
    bmDen = bmDen / g;

    B.push([bmNum, bmDen]);
  }

  const result = B[nNum]!;
  return { numerator: result[0], denominator: result[1] };
}

/**
 * Compute the factorial of n, which is the product 1 * 2 * 3 * ... * (n-1) * n.
 *
 * @param n - Non-negative integer
 * @param algorithm - Algorithm to use ('gmp' or 'pari')
 * @returns n!
 *
 * @example
 * ```typescript
 * factorial(5n)  // 120n
 * factorial(0n)  // 1n
 * ```
 *
 * @see Reference: sage/arith/misc.py:factorial
 * @see Deviation: PARI factorial real representation and transcendental dependencies
 */
export function factorial(n: IntegerLike, algorithm?: 'gmp'): bigint;
export function factorial(n: IntegerLike, algorithm: 'pari'): MpReal<bigint>;
export function factorial(n: IntegerLike, algorithm?: 'gmp' | 'pari'): bigint | MpReal<bigint>;
export function factorial(n: IntegerLike, algorithm: 'gmp' | 'pari' = 'gmp'): bigint | MpReal<bigint> {
  const _n = toBigInt(n);
  if (_n < 0n) throw new ValueError('factorial -- must be nonnegative');
  if (algorithm === 'gmp') return new Integer(_n).factorial().value;
  if (algorithm !== 'pari') throw new ValueError('unknown algorithm');
  // cypari2 converts the Python argument to a signed C long before mpfactr.
  if (_n >= 1n << 63n) throw new OverflowError('Python int too large to convert to C long');
  return pari_mpfactr(_n, 64);
}

/**
 * Test whether n is a pseudo-prime.
 *
 * The result is NOT proven correct - this is a pseudo-primality test!
 *
 * SageMath (`sage/arith/misc.py:607` -> `sage/rings/integer.pyx:5607`) delegates
 * to PARI's ``ispseudoprime``, i.e. the Baillie-PSW test.  We delegate to the
 * very same BPSW entry point that `is_prime` uses, so the two agree everywhere.
 * A fixed set of Miller-Rabin bases would not: 318665857834031151167461 is
 * composite yet a strong pseudoprime to all of the first 12 prime bases.
 *
 * @param n - Integer to test
 * @returns true if n is a pseudo-prime
 *
 * @example
 * ```typescript
 * is_pseudoprime(389n)   // true
 * is_pseudoprime(2000n)  // false
 * is_pseudoprime(2n)     // true
 * is_pseudoprime(-1n)    // false
 * ```
 *
 * @see Reference: sage/arith/misc.py:is_pseudoprime
 * @see Implementation: parigp-ts/src/ifactor.ts (BPSW, port of pari/src/basemath/ifactor1.c)
 */
export function is_pseudoprime(n: IntegerLike): boolean {
  const _n = toBigInt(n);
  // Negative numbers and values <= 1 are not prime
  if (_n <= 1n) {
    return false;
  }

  // Delegate to PARI's BPSW implementation (same entry point as is_prime).
  return pari_isPrime(_n);
}

/**
 * Test if n is a power of a pseudoprime.
 *
 * The result is NOT proven correct - this IS a pseudo-primality test!
 *
 * SageMath calls ``ZZ(n).is_prime_power(proof=False, get_data=get_data)``, i.e.
 * PARI's ``ispseudoprimepower``; that differs from ``isprimepower`` only in
 * using BPSW instead of a primality proof, so it shares our implementation.
 *
 * @param n - Integer to test
 * @param get_data - If true, return (p, k) such that n = p^k
 * @returns Whether n is a pseudoprime power
 *
 * @see Reference: sage/arith/misc.py:is_pseudoprime_power
 */
export function is_pseudoprime_power(n: IntegerLike, get_data?: false): boolean;
export function is_pseudoprime_power(n: IntegerLike, get_data: true): [bigint, bigint];
export function is_pseudoprime_power(
  n: IntegerLike,
  get_data: boolean = false
): boolean | [bigint, bigint] {
  const value = toBigInt(n);
  const data = pari_isprimepower(value);
  if (data === null) {
    return get_data ? [value, 0n] : false;
  }
  return get_data ? [data[0], BigInt(data[1])] : true;
}

/**
 * List of all positive prime powers between start and stop-1, inclusive.
 *
 * @param start - Lower bound (or upper bound if stop is undefined)
 * @param stop - Upper bound (exclusive)
 * @returns Array of prime powers
 *
 * @example
 * ```typescript
 * prime_powers(20n)  // [2n, 3n, 4n, 5n, 7n, 8n, 9n, 11n, 13n, 16n, 17n, 19n]
 * ```
 *
 * @see Reference: sage/arith/misc.py:prime_powers
 */
export function prime_powers(start: IntegerLike, stop?: IntegerLike): bigint[] {
  let lower = toBigInt(start);
  const upper = stop === undefined ? lower : toBigInt(stop);
  if (stop === undefined) lower = 2n;
  if (upper <= 2n || lower >= upper) return [];

  // SageMath walks the primes below ``stop`` and emits their powers, rather
  // than testing every integer for prime-power-ness.
  const output: bigint[] = [];
  for (const p of prime_range(upper)) {
    let q = p;
    while (q < lower) {
      q *= p;
    }
    while (q < upper) {
      output.push(q);
      q *= p;
    }
  }

  output.sort((a, b) => (a < b ? -1 : a > b ? 1 : 0));
  return output;
}

/**
 * Return the first n primes.
 *
 * @param n - Non-negative integer
 * @returns List of the first n prime numbers
 *
 * @example
 * ```typescript
 * primes_first_n(10)  // [2n, 3n, 5n, 7n, 11n, 13n, 17n, 19n, 23n, 29n]
 * ```
 *
 * @see Reference: sage/arith/misc.py:primes_first_n
 */
export function primes_first_n(n: IntegerLike | number): bigint[] {
  // Keep the existing numeric-count API. Sage checks these bounds before PARI
  // converts a positive fractional count to its truncated integer value.
  const value = typeof n === 'number' ? n : toBigInt(n);
  if (value < 0) throw new ValueError('n must be nonnegative');
  if (value < 1) return [];
  if (typeof value === 'number' && !Number.isFinite(value)) {
    if (Number.isNaN(value)) throw new ValueError('cannot convert float NaN to integer');
    throw new OverflowError('cannot convert float infinity to integer');
  }
  const count = typeof value === 'number' ? BigInt(Math.trunc(value)) : value;
  return prime_range(nth_prime(count) + 1n);
}

/**
 * Return a list of the primes <= n.
 *
 * This is extremely slow and is for educational purposes only.
 * Uses the Sieve of Eratosthenes algorithm.
 *
 * @param n - Positive integer
 * @returns List of primes <= n
 *
 * @see Reference: sage/arith/misc.py:eratosthenes
 */
export function eratosthenes(n: IntegerLike): bigint[] {
  const value = toBigInt(n);
  if (value < 2n) return [];
  if (value === 2n) return [2n];
  const limit = toSafeNumber(value);
  const half = Math.floor((limit + 1) / 2);
  const sieve = Array.from({ length: half }, (_, i) => 2 * i + 3);
  const root = toSafeNumber(isqrt(value));
  for (let i = 0, m = 3; m <= root; i++, m = 2 * i + 3) {
    if (sieve[i]) for (let j = (m * m - 3) / 2; j < half; j += m) sieve[j] = 0;
  }
  return [2n, ...sieve.filter((x) => x !== 0 && x <= limit).map(BigInt)];
}

/**
 * Return an iterator over all primes between start and stop-1, inclusive.
 *
 * @param start - Lower bound, or exclusive upper bound when stop is omitted (default: 2)
 * @param stop - Exclusive upper bound; Infinity selects an unbounded iterator
 * @param proof - Whether to use the proven or probable-prime path
 * @returns Iterator over primes
 *
 * @see Reference: sage/arith/misc.py:primes
 */
export function* primes(
  start: IntegerLike = 2n,
  stop?: IntegerLike | number,
  proof?: boolean
): Generator<bigint, void, unknown> {
  let lower = toBigInt(start);
  const upper =
    stop === undefined ? lower : stop === Infinity ? Infinity : toBigInt(stop as IntegerLike);
  if (stop === undefined) lower = 2n;
  let n = lower - 1n;
  // The original advances next_prime once per draw, even for finite bounds.
  while (true) {
    n = proof === false ? next_probable_prime(n) : next_prime(n);
    if (n < upper) yield n;
    else return;
  }
}

/**
 * Return the smallest prime power greater than n.
 *
 * @param n - Integer
 * @returns Smallest prime power > n
 *
 * @example
 * ```typescript
 * next_prime_power(7n)  // 8n
 * next_prime_power(10n) // 11n
 * next_prime_power(1n)  // 2n
 * next_prime_power(99n) // 101n
 * ```
 *
 * @see Reference: sage/arith/misc.py:next_prime_power
 */
export function next_prime_power(n: IntegerLike): bigint {
  const value = toBigInt(n);
  if (value < 2n) return 2n;
  const bound = 1n << BigInt(value.toString(2).length);
  let candidate = value + ((value & 1n) === 0n ? 1n : 2n);
  for (; candidate < bound; candidate += 2n) if (is_prime_power(candidate)) return candidate;
  return bound;
}

/**
 * Return the next probable prime after n, as determined by pseudoprimality test.
 *
 * @param n - Integer
 * @returns Next probable prime
 *
 * @see Reference: sage/arith/misc.py:next_probable_prime
 */
export function next_probable_prime(n: IntegerLike): bigint {
  return pari_nextprime(toBigInt(n) + 1n);
}

/**
 * Return the largest prime power smaller than n.
 *
 * @param n - Integer > 2
 * @returns Largest prime power < n
 * @throws {ValueError} If n <= 2
 *
 * @example
 * ```typescript
 * previous_prime_power(3n)   // 2n
 * previous_prime_power(10n)  // 9n
 * previous_prime_power(7n)   // 5n
 * previous_prime_power(127n) // 125n
 * ```
 *
 * @see Reference: sage/arith/misc.py:previous_prime_power
 */
export function previous_prime_power(n: IntegerLike): bigint {
  const value = toBigInt(n);
  if (value <= 2n) throw new ValueError('no prime power less than 2');
  const previous = value - 1n;
  const bound = 1n << BigInt(previous.toString(2).length - 1);
  let candidate = (previous & 1n) === 0n ? previous - 1n : previous;
  for (; candidate > bound; candidate -= 2n) if (is_prime_power(candidate)) return candidate;
  return bound;
}

/**
 * Return a random prime p between lbound and n (inclusive).
 *
 * Uses rejection sampling: repeatedly picks random numbers in the range
 * and tests for primality until a prime is found.
 *
 * @param n - Integer >= 2, upper bound
 * @param proof - Whether to use proven primality test (default: true)
 * @param lbound - Lower bound (default: 2)
 * @returns Random prime in [lbound, n]
 * @throws {ValueError} If n < 2 or if no primes exist in the range
 *
 * @example
 * ```typescript
 * random_prime(100n)            // some prime <= 100
 * random_prime(200n, true, 100n) // some prime in [100, 200]
 * ```
 *
 * @see Reference: sage/arith/misc.py:random_prime
 */
export function random_prime(n: bigint, proof: boolean = true, lbound: bigint = 2n): bigint {
  if (n < 2n) {
    throw new ValueError('n must be greater than or equal to 2');
  }

  if (n < lbound) {
    throw new ValueError('n must be at least lbound: ' + lbound.toString());
  }

  if (lbound < 2n) {
    lbound = 2n;
  }

  // Special case: if n == 2 and lbound <= 2
  if (n === 2n) {
    return 2n;
  }

  // Check that there exists a prime in the range using Bertrand's postulate
  // For small ranges, verify by finding first prime
  const primeTest = proof ? is_prime : is_pseudoprime;

  // Find the smallest prime >= lbound to verify primes exist in range
  let firstPrime = lbound;
  if (firstPrime <= 2n) {
    firstPrime = 2n;
  } else if ((firstPrime & 1n) === 0n) {
    firstPrime++;
  }
  while (firstPrime <= n && !primeTest(firstPrime)) {
    if (firstPrime === 2n) {
      firstPrime = 3n;
    } else {
      firstPrime += 2n;
    }
  }
  if (firstPrime > n) {
    throw new ValueError(
      'there are no primes between ' + lbound.toString() + ' and ' + n.toString() + ' (inclusive)'
    );
  }

  // Rejection sampling: pick random numbers until we find a prime.
  // Mirror Sage's use of ZZ.random_element(lbound, n), which samples in [lbound, n).
  const range = n - lbound;
  if (range <= 0n) {
    throw new TypeError('x must be < y');
  }
  const rstate = current_randstate();

  const getRandomInRange = (): bigint => {
    return lbound + rstate.random_below(range);
  };

  // Keep trying until we find a prime.
  while (true) {
    const p = getRandomInRange();
    if (primeTest(p)) {
      return p;
    }
  }
}

/**
 * Extended lcm function: given two positive integers m, n, returns
 * a triple (l, m1, n1) such that l = lcm(m, n) = m1 * n1 where
 * m1 | m, n1 | n and gcd(m1, n1) = 1, all with no factorization.
 *
 * @param m - Positive integer
 * @param n - Positive integer
 * @returns Triple (lcm, m1, n1)
 *
 * @example
 * ```typescript
 * xlcm(12n, 18n)  // [36n, 4n, 9n] since 36 = 4 * 9, 4 | 12, 9 | 18, gcd(4,9) = 1
 * ```
 *
 * @see Reference: sage/arith/misc.py:xlcm
 */
export function xlcm(m: IntegerLike, n: IntegerLike): [bigint, bigint, bigint] {
  let left = toBigInt(m);
  const right = toBigInt(n);
  let g = gcd(left, right);
  if (g === 0n) throw new ZeroDivisionError('Integer division by zero');
  const multiple = (left * right) / g;
  g = gcd(left, right / g);
  while (g !== 1n) {
    left /= g;
    g = gcd(left, g);
  }
  if (left === 0n) throw new ZeroDivisionError('Integer division by zero');
  return [multiple, left, multiple / left];
}

/**
 * Return a CRT basis for the given moduli.
 *
 * Given moduli [m1, m2, ..., mn], returns basis elements [e1, e2, ..., en]
 * such that e_i ≡ 1 (mod m_i) and e_i ≡ 0 (mod m_j) for j ≠ i.
 *
 * This is useful for precomputation in repeated CRT applications.
 *
 * When `require_coprime_moduli` is false and the moduli are *not* pairwise
 * coprime, SageMath falls back to the construction
 * `e_i = CRT(0, 1, M_i/d_i, m_i/d_i)` (with `M_i` the running lcm and
 * `d_i = gcd(M_i, m_i)`) combined with a table of partial products of
 * `1 - e_i`; the resulting `a_i` are *not* reduced modulo the lcm and only
 * solve the system when a solution exists.
 *
 * @param moduli - List of moduli (must be pairwise coprime if require_coprime_moduli is true)
 * @param require_coprime_moduli - Whether moduli must be pairwise coprime (default: true)
 * @returns List of CRT basis elements, or [basis, coprime_flag] if require_coprime_moduli is false
 *
 * @example
 * ```typescript
 * // For coprime moduli [5, 13]:
 * const [c1, c2] = CRT_basis([5n, 13n]);
 * // c1 ≡ 1 (mod 5), c1 ≡ 0 (mod 13)
 * // c2 ≡ 0 (mod 5), c2 ≡ 1 (mod 13)
 * // To combine residues a1 mod 5 and a2 mod 13:
 * // result = (a1 * c1 + a2 * c2) % 65
 *
 * CRT_basis([60n, 90n, 150n], false);  // [[15n, -20n, 6n], false]
 * ```
 *
 * @see Reference: sage/arith/misc.py:CRT_basis
 */
export function CRT_basis(
  moduli: IntegerLike[],
  require_coprime_moduli: boolean = true
): bigint[] | [bigint[], boolean] {
  const n = moduli.length;
  const mods = moduli.map(toBigInt);
  if (n === 0) {
    return [];
  }

  const cs: bigint[] = [];
  let coprime = true;

  // Compute M = product of all moduli
  let M = 1n;
  for (const m of mods) {
    M *= m;
  }

  for (const m of mods) {
    if (m === 0n) throw new ZeroDivisionError('Integer division by zero');
    const Mm = M / m;
    const [d, , v] = xgcd(m, Mm);
    if (d !== 1n) {
      if (require_coprime_moduli) {
        throw new ValueError('moduli must be coprime');
      }
      coprime = false;
      break;
    }
    // e_i = v * M_i mod M, where M_i = M / m_i
    cs.push(crtRemainder(v * Mm, M));
  }

  if (coprime) {
    return require_coprime_moduli ? cs : [cs, true];
  }

  // Preserve the bundled source's partial prefix when entering the fallback.
  // It is observable through both CRT_basis and CRT_vectors.
  const e: bigint[] = [1n];
  let M_i = mods[0]!;
  for (let i = 1; i < n; i++) {
    const m_i = mods[i]!;
    const d_i = gcd(M_i, m_i);
    e.push(crt(0n, 1n, M_i / d_i, m_i / d_i));
    M_i = lcm(M_i, m_i);
  }

  const partial_prod_table: bigint[] = [1n];
  for (let i = 1; i < n; i++) {
    partial_prod_table.push((1n - e[n - i]!) * partial_prod_table[i - 1]!);
  }

  for (let i = 0; i < n; i++) {
    cs.push(e[i]! * partial_prod_table[n - i - 1]!);
  }

  return [cs, false];
}

/**
 * Vector form of the Chinese Remainder Theorem.
 *
 * Given a list of integer vectors X = [v1, v2, ..., vn] and moduli [m1, m2, ..., mn],
 * returns a vector w such that w ≡ vi (mod mi) for all i (component-wise).
 *
 * This is more efficient than applying CRT to each component separately when
 * the same moduli are used repeatedly.
 *
 * @param X - List of vectors of integers (all must have the same length)
 * @param moduli - List of moduli (must have same length as X)
 * @returns Vector satisfying CRT conditions
 * @throws {ValueError} If no solution exists or if lengths don't match
 *
 * @example
 * ```typescript
 * CRT_vectors([[3n, 5n, 7n], [3n, 5n, 11n]], [2n, 3n]);
 * // Returns [3n, 5n, 5n] - each component satisfies CRT
 * ```
 *
 * @see Reference: sage/arith/misc.py:CRT_vectors
 */
export function CRT_vectors(X: IntegerLike[][], moduli: IntegerLike[]): bigint[] {
  if (X.length === 0 || X[0]!.length === 0) return [];
  const n = X.length;
  if (n !== moduli.length) throw new ValueError('number of moduli must equal length of X');
  const [basis, coprime] = CRT_basis(moduli, false) as [bigint[], boolean];
  const mods = moduli.map(toBigInt);
  const modulus = lcm(mods);
  const candidate: bigint[] = [];
  for (let j = 0; j < X[0]!.length; j++) {
    let sum = 0n;
    for (let i = 0; i < n; i++) {
      if (j >= X[i]!.length) throw new IndexError('list index out of range');
      sum += basis[i]! * toBigInt(X[i]![j]!);
    }
    candidate.push(crtRemainder(sum, modulus));
  }
  if (!coprime) {
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < X[i]!.length; j++) {
        if (j >= candidate.length) throw new IndexError('list index out of range');
        if (crtRemainder(toBigInt(X[i]![j]!) - candidate[j]!, mods[i]!) !== 0n) {
          throw new ValueError('solution does not exist');
        }
      }
    }
  }
  return candidate;
}

/**
 * Rational reconstruction: recover p/q from a mod m.
 *
 * Given a and m, attempts to find integers p and q such that:
 * - a ≡ p/q (mod m)
 * - |p|, q < sqrt(m/2)
 * - gcd(p, q) = 1
 *
 * This uses the extended Euclidean algorithm with early termination,
 * as described in Knuth Vol 2, 3rd ed, pages 656-657.
 *
 * @param a - The residue (IntegerLike)
 * @param m - The modulus (must be positive, IntegerLike)
 * @returns [p, q] where a ≡ p/q (mod m), or throws if no solution exists
 * @throws {ZeroDivisionError} If m is zero
 * @throws {ArithmeticError} If no valid rational reconstruction exists
 *
 * @example
 * ```typescript
 * // 11323 ≡ 119/53 (mod 100000) since 119 * inverse_mod(53, 100000) = 11323
 * rational_reconstruction(11323n, 100000n);  // [119n, 53n]
 * rational_reconstruction(11323, 100000);    // [119n, 53n] (with numbers)
 *
 * // Small values are found directly
 * rational_reconstruction(3n, 292393n);  // [3n, 1n]
 * ```
 *
 * @see Reference: sage/arith/misc.py:rational_reconstruction
 */
export function rational_reconstruction(a: IntegerLike, m: IntegerLike): [bigint, bigint] {
  let _a = toBigInt(a);
  let _m = toBigInt(m);

  if (_m === 0n) {
    throw new ZeroDivisionError('rational reconstruction with zero modulus');
  }

  // Work with positive modulus
  _m = _m < 0n ? -_m : _m;

  // Reduce a mod m to be in [0, m)
  _a = ((_a % _m) + _m) % _m;

  // Special case: a = 0 means p/q = 0/1
  if (_a === 0n) {
    return [0n, 1n];
  }

  // Compute bound = floor(sqrt(m/2))
  // We want |p| <= bound and q <= bound
  const bound = isqrt(_m / 2n);

  // Extended Euclidean algorithm with early termination
  // We maintain: v3 = v2 * a + (some multiple of m)
  // When v3 <= bound, we check if we have a valid solution
  let u2 = 0n;
  let v2 = 1n;
  let u3 = _m;
  let v3 = _a;

  while (v3 > bound) {
    const quot = u3 / v3;
    const t2 = u2 - quot * v2;
    const t3 = u3 - quot * v3;
    u2 = v2;
    u3 = v3;
    v2 = t2;
    v3 = t3;
  }

  // Now v3 <= bound
  // The candidate is (p, q) = (v3 * sign(v2), |v2|)
  // We need to verify:
  // 1. |v2| <= bound
  // 2. gcd(v3, |v2|) = 1

  const q = v2 < 0n ? -v2 : v2;
  const p = v2 < 0n ? -v3 : v3;

  // Check the bound on q
  if (q > bound) {
    throw new ArithmeticError(
      `rational reconstruction of ${((_a % _m) + _m) % _m} (mod ${_m}) does not exist`
    );
  }

  // Check coprimality
  if (gcd(p < 0n ? -p : p, q) !== 1n) {
    throw new ArithmeticError(
      `rational reconstruction of ${((_a % _m) + _m) % _m} (mod ${_m}) does not exist`
    );
  }

  return [p, q];
}

/**
 * Half-GCD algorithm: compute a transformation matrix for GCD computation.
 *
 * Given integers a >= b >= 0, returns a 2x2 matrix [[r, s], [t, u]] such that:
 * - [a', b'] = [[r, s], [t, u]] * [a, b] (matrix-vector multiplication)
 * - a' and b' are the values at approximately the halfway point of the
 *   extended Euclidean algorithm
 *
 * This is a key subroutine in fast (quasi-linear) GCD algorithms and
 * is useful for fast rational reconstruction.
 *
 * @param a - First integer (a >= b >= 0)
 * @param b - Second integer
 * @returns Transformation matrix [[r, s], [t, u]] as a flat array [r, s, t, u]
 *
 * @example
 * ```typescript
 * const [r, s, t, u] = half_gcd(1000n, 37n);
 * // The matrix [[r, s], [t, u]] transforms [1000, 37] toward their GCD
 * ```
 *
 * @see Reference: Modern Computer Algebra, von zur Gathen & Gerhard
 */
export function half_gcd(a: bigint, b: bigint): [bigint, bigint, bigint, bigint] {
  // Ensure a >= b >= 0
  if (a < 0n) a = -a;
  if (b < 0n) b = -b;
  if (a < b) {
    [a, b] = [b, a];
  }

  // Base case: if b is small enough, return identity matrix
  if (b === 0n) {
    return [1n, 0n, 0n, 1n];
  }

  // Compute the "half" point - we want to reduce until the result
  // has about half the bits of the original
  const targetBits = BigInt(a.toString(2).length) / 2n;
  const targetBound = 1n << targetBits;

  // Run extended Euclidean algorithm until we reach the halfway point
  let r0 = 1n;
  let s0 = 0n;
  let r1 = 0n;
  let s1 = 1n;
  let a0 = a;
  let a1 = b;

  while (a1 >= targetBound && a1 !== 0n) {
    const q = a0 / a1;
    const temp_a = a0 - q * a1;
    const temp_r = r0 - q * r1;
    const temp_s = s0 - q * s1;

    a0 = a1;
    a1 = temp_a;
    r0 = r1;
    r1 = temp_r;
    s0 = s1;
    s1 = temp_s;
  }

  return [r0, s0, r1, s1];
}

/**
 * Compute the continued fraction expansion of a rational number.
 *
 * For a rational number p/q, returns the continued fraction [a0; a1, a2, ..., an]
 * such that p/q = a0 + 1/(a1 + 1/(a2 + ... + 1/an)).
 *
 * @param x - A RationalLike value (Rational, bigint, number, or Integer), or numerator if q is provided
 * @param q - Denominator (default: 1n, ignored if x is a Rational)
 * @param bound - Optional maximum number of terms to compute
 * @returns Array of partial quotients [a0, a1, a2, ...]
 *
 * @example
 * ```typescript
 * continued_fraction(13n, 9n);  // [1n, 2n, 4n] since 13/9 = 1 + 1/(2 + 1/4)
 * continued_fraction(225n, 157n);  // [1n, 2n, 3n, 4n, 5n]
 * continued_fraction(-1n, 3n);  // [-1n, 1n, 2n] since -1/3 = -1 + 1/(1 + 1/2)
 * continued_fraction(new Rational(13n, 9n));  // [1n, 2n, 4n]
 * ```
 *
 * @see Reference: sage/rings/rational.pyx:continued_fraction_list
 */
export function continued_fraction(x: RationalLike, q?: IntegerLike, bound?: number): bigint[] {
  let p: bigint;
  let _q: bigint;

  if (x instanceof Rational) {
    // If x is a Rational, use its numerator and denominator directly
    p = x.numerator;
    _q = x.denominator;
  } else {
    // x is IntegerLike
    p = toBigInt(x);
    _q = q !== undefined ? toBigInt(q) : 1n;
  }

  if (_q === 0n) {
    throw new ZeroDivisionError('denominator cannot be zero');
  }

  // Normalize so q > 0
  if (_q < 0n) {
    p = -p;
    _q = -_q;
  }

  const result: bigint[] = [];
  const maxTerms = bound ?? Number.POSITIVE_INFINITY;
  let termCount = 0;

  while (_q !== 0n && termCount < maxTerms) {
    // Use floor division (fdiv)
    // For negative numbers, we need proper floor division
    let quotient: bigint;
    if (p >= 0n) {
      quotient = p / _q;
    } else {
      // Floor division for negative: -((-p-1)/q + 1) when p < 0
      quotient = -((-p - 1n) / _q + 1n);
    }

    result.push(quotient);

    // Update: (p, q) <- (q, p - quotient * q)
    const newP = _q;
    const newQ = p - quotient * _q;
    p = newP;
    _q = newQ;
    termCount++;
  }

  return result;
}

/**
 * Compute the convergents of a continued fraction.
 *
 * Given a continued fraction [a0; a1, a2, ..., an], returns the sequence
 * of convergents p_i/q_i as pairs [p_i, q_i].
 *
 * The convergents satisfy:
 * - p_{-1} = 1, p_0 = a_0, p_n = a_n * p_{n-1} + p_{n-2}
 * - q_{-1} = 0, q_0 = 1, q_n = a_n * q_{n-1} + q_{n-2}
 *
 * @param cf - Array of partial quotients [a0, a1, a2, ...]
 * @returns Array of convergent pairs [[p0, q0], [p1, q1], ...]
 *
 * @example
 * ```typescript
 * convergents([1n, 2n, 4n]);
 * // Returns [[1n, 1n], [3n, 2n], [13n, 9n]]
 * // Representing 1/1, 3/2, 13/9
 *
 * convergents([0n, 6n, 1n, 4n, 1n, 3n]);
 * // Returns [[0n, 1n], [1n, 6n], [1n, 7n], [5n, 34n], [6n, 41n], [23n, 157n]]
 * // For the continued fraction of 23/157
 * ```
 *
 * @see Reference: sage/rings/continued_fraction.py:convergents
 */
export function convergents(cf: bigint[]): Array<[bigint, bigint]> {
  if (cf.length === 0) {
    return [];
  }

  const result: Array<[bigint, bigint]> = [];

  // Initialize: p_{-1} = 1, p_0 = a_0
  //             q_{-1} = 0, q_0 = 1
  let pPrev = 1n;
  let pCurr = cf[0]!;
  let qPrev = 0n;
  let qCurr = 1n;

  result.push([pCurr, qCurr]);

  for (let i = 1; i < cf.length; i++) {
    const a = cf[i]!;

    // p_n = a_n * p_{n-1} + p_{n-2}
    // q_n = a_n * q_{n-1} + q_{n-2}
    const pNext = a * pCurr + pPrev;
    const qNext = a * qCurr + qPrev;

    pPrev = pCurr;
    pCurr = pNext;
    qPrev = qCurr;
    qCurr = qNext;

    result.push([pCurr, qCurr]);
  }

  return result;
}

/**
 * Evaluate a continued fraction to get its rational value.
 *
 * @param cf - Array of partial quotients [a0, a1, a2, ...]
 * @returns [numerator, denominator] of the rational value
 *
 * @example
 * ```typescript
 * continued_fraction_value([1n, 2n, 4n]);  // [13n, 9n]
 * continued_fraction_value([3n, 7n, 16n]); // [355n, 113n] (approximation to pi)
 * ```
 */
export function continued_fraction_value(cf: bigint[]): [bigint, bigint] {
  if (cf.length === 0) {
    throw new ValueError('continued fraction cannot be empty');
  }

  const convs = convergents(cf);
  return convs[convs.length - 1]!;
}

/**
 * Return the binomial coefficient C(x, m) = x(x-1)...(x-m+1) / m!
 *
 * @param x - Integer or number
 * @param m - Non-negative integer
 * @returns Binomial coefficient
 *
 * @example
 * ```typescript
 * binomial(5n, 2n)  // 10n
 * binomial(20n, 10n) // 184756n
 * ```
 *
 * @see Reference: sage/arith/misc.py:binomial
 */
export function binomial(x: IntegerLike, m: IntegerLike): bigint {
  const _x = toBigInt(x);
  let _m = toBigInt(m);
  // If m < 0, return 0 (SageMath behavior)
  if (_m < 0n) {
    return 0n;
  }

  // If m == 0, return 1
  if (_m === 0n) {
    return 1n;
  }

  // For negative x, use the identity: binomial(-n, k) = (-1)^k * binomial(n+k-1, k)
  if (_x < 0n) {
    const sign = _m % 2n === 0n ? 1n : -1n;
    return sign * binomial(-_x + _m - 1n, _m);
  }

  // If m > x (for non-negative x), return 0
  if (_m > _x) {
    return 0n;
  }

  // Optimization: use symmetry C(n, k) = C(n, n-k) to minimize multiplications
  if (_m > _x - _m) {
    _m = _x - _m;
  }

  // Compute the binomial coefficient using the product formula
  // C(x, m) = x * (x-1) * ... * (x-m+1) / (m!)
  // We compute this incrementally to avoid overflow in intermediate results
  let result = 1n;
  for (let i = 0n; i < _m; i++) {
    result = (result * (_x - i)) / (i + 1n);
  }

  return result;
}

/**
 * Return the multinomial coefficient.
 *
 * multinomial(k1, k2, ..., kn) = (k1 + k2 + ... + kn)! / (k1! * k2! * ... * kn!)
 *
 * @param ks - Non-negative integers
 * @returns Multinomial coefficient
 *
 * @example
 * ```typescript
 * multinomial(2n, 3n, 4n)  // 1260n (ways to arrange 2 A's, 3 B's, 4 C's)
 * ```
 *
 * @see Reference: sage/arith/misc.py:multinomial
 */
export function multinomial(...ks: IntegerLike[]): bigint {
  if (ks.length === 0) {
    return 1n;
  }

  const _ks = ks.map(toBigInt);
  for (const k of _ks) {
    if (k < 0n) {
      return 0n;
    }
  }

  // Sum of all k's
  let total = 0n;
  for (const k of _ks) {
    total += k;
  }

  // Compute as product of binomial coefficients
  // C(k1+k2+...+kn, k1) * C(k2+...+kn, k2) * ... * C(kn, kn)
  let result = 1n;
  let remaining = total;

  for (const k of _ks) {
    result = result * binomial(remaining, k);
    remaining -= k;
  }

  return result;
}

/**
 * Return a dictionary containing pairs {(k1, k2): C(n, k1)} where C(n, k1)
 * are binomial coefficients and k1 + k2 = n.
 *
 * @param n - Non-negative integer
 * @returns Map from (k1, k2) pairs to binomial coefficients
 *
 * @example
 * ```typescript
 * binomial_coefficients(3n)
 * // Map { [0, 3] => 1, [1, 2] => 3, [2, 1] => 3, [3, 0] => 1 }
 * ```
 *
 * @see Reference: sage/arith/misc.py:binomial_coefficients
 */
export function binomial_coefficients(n: bigint): Map<string, bigint> {
  const result = new Map<string, bigint>();

  if (n < 0n) {
    return result;
  }

  for (let k = 0n; k <= n; k++) {
    const key = `${k},${n - k}`;
    result.set(key, binomial(n, k));
  }

  return result;
}

/**
 * Return a dictionary containing multinomial coefficients.
 *
 * Returns all multinomial(k1, k2, ..., km) where k1 + k2 + ... + km = n.
 *
 * @param m - Number of parts
 * @param n - Sum of parts
 * @returns Map from tuples (k1, ..., km) to multinomial coefficients
 *
 * @see Reference: sage/arith/misc.py:multinomial_coefficients
 */
export function multinomial_coefficients(m: bigint, n: bigint): Map<string, bigint> {
  const result = new Map<string, bigint>();

  if (m <= 0n || n < 0n) {
    return result;
  }

  // Generate all partitions of n into m non-negative parts
  function* partitions(sum: bigint, parts: bigint, prefix: bigint[]): Generator<bigint[]> {
    if (parts === 1n) {
      yield [...prefix, sum];
      return;
    }

    for (let k = 0n; k <= sum; k++) {
      yield* partitions(sum - k, parts - 1n, [...prefix, k]);
    }
  }

  for (const tuple of partitions(n, m, [])) {
    const key = tuple.join(',');
    result.set(key, multinomial(...tuple));
  }

  return result;
}

/**
 * Return a positive integer that generates the multiplicative group
 * of integers modulo n, if one exists; otherwise, raise a ValueError.
 *
 * A primitive root exists if n=4 or n=p^k or n=2p^k, where p is an odd prime.
 *
 * @param n - Nonzero integer
 * @param check - Whether to check if n has a primitive root
 * @returns A primitive root of n
 *
 * @see Reference: sage/arith/misc.py:primitive_root
 */
export function primitive_root(n: bigint, check: boolean = true): bigint {
  // Use absolute value
  n = n < 0n ? -n : n;

  // n = 0 has no primitive root
  if (n === 0n) {
    throw new ValueError('no primitive root');
  }

  // Special cases for small n: n in {1, 2, 3, 4}
  // n-1 is a primitive root for these
  if (n <= 4n) {
    return n - 1n;
  }

  // Check if n has a primitive root (only if check is true)
  // A primitive root exists iff n = 1, 2, 4, p^k, or 2*p^k for odd prime p
  if (check) {
    let hasPrimitiveRoot = false;

    if (n % 2n === 1n) {
      // n is odd: check if n = p^k for some prime p
      hasPrimitiveRoot = is_prime_power(n);
    } else {
      // n is even
      const m = n / 2n;
      if (m % 2n === 1n) {
        // n = 2 * m where m is odd: check if m = p^k for some prime p
        hasPrimitiveRoot = m === 1n || is_prime_power(m);
      }
      // Otherwise n = 2^k * m where k >= 2, which has no primitive root
    }

    if (!hasPrimitiveRoot) {
      throw new ValueError('no primitive root');
    }
  }

  // Compute phi(n) to find the order of the multiplicative group
  const phi = euler_phi(n);

  // Get the prime factorization of phi(n) for testing orders
  const phiFactors = factor(phi);
  const primeFactorsOfPhi: bigint[] = [];
  for (const [p] of phiFactors) {
    if (p > 0n) {
      primeFactorsOfPhi.push(p);
    }
  }

  // Find the smallest primitive root by testing candidates
  // g is a primitive root mod n iff g^(phi(n)/p) != 1 mod n for all prime p | phi(n)
  for (let g = 2n; g < n; g++) {
    // Check if gcd(g, n) = 1
    if (gcd(g, n) !== 1n) {
      continue;
    }

    // Check if g is a primitive root
    let isPrimitiveRoot = true;
    for (const p of primeFactorsOfPhi) {
      const exp = phi / p;
      if (power_mod(g, exp, n) === 1n) {
        isPrimitiveRoot = false;
        break;
      }
    }

    if (isPrimitiveRoot) {
      return g;
    }
  }

  throw new ValueError('no primitive root');
}

/**
 * Return the n-th prime number (1-indexed, so 2 is the 1st prime).
 *
 * @param n - Positive integer
 * @returns The n-th prime
 *
 * @example
 * ```typescript
 * nth_prime(1n)  // 2n
 * nth_prime(10n) // 29n
 * ```
 *
 * @see Reference: sage/arith/misc.py:nth_prime
 */
export function nth_prime(n: IntegerLike): bigint {
  const value = toBigInt(n);
  if (value <= 0n) {
    throw new ValueError('nth prime meaningless for nonpositive n (=' + value.toString() + ')');
  }
  return pari_prime(value);
}

/**
 * Return a sorted list of all squares modulo the integer n in the range 0 <= x < |n|.
 *
 * @param n - Integer
 * @returns List of quadratic residues
 *
 * @example
 * ```typescript
 * quadratic_residues(11n)  // [0n, 1n, 3n, 4n, 5n, 9n]
 * quadratic_residues(1n)   // [0n]
 * quadratic_residues(2n)   // [0n, 1n]
 * quadratic_residues(8n)   // [0n, 1n, 4n]
 * ```
 *
 * @see Reference: sage/arith/misc.py:quadratic_residues
 */
export function quadratic_residues(n: IntegerLike): bigint[] {
  let value = toBigInt(n);
  if (value < 0n) value = -value;
  if (value === 0n) throw new ZeroDivisionError('integer modulo by zero');
  const residues = new Set<bigint>();
  for (let a = 0n; a <= value / 2n; a++) residues.add((a * a) % value);
  return [...residues].sort((a, b) => (a < b ? -1 : a > b ? 1 : 0));
}

/**
 * Return the continuant of the sequence v.
 *
 * The continuant is defined by:
 * - K_0() = 1
 * - K_1(x_1) = x_1
 * - K_n(x_1, ..., x_n) = K_{n-1}(x_1, ..., x_{n-1}) * x_n + K_{n-2}(x_1, ..., x_{n-2})
 *
 * The continuant is the numerator of the continued fraction [x_1; x_2, ..., x_n].
 *
 * @param v - List or tuple of elements
 * @param n - Optional length (uses first n elements of v)
 * @returns The continuant
 *
 * @example
 * ```typescript
 * continuant([1n, 2n, 3n])  // 10n (since [1; 2, 3] = 10/7)
 * ```
 *
 * @see Reference: sage/arith/misc.py:continuant
 */
export function continuant(v: IntegerLike[], n?: IntegerLike): bigint {
  const requested = n === undefined ? BigInt(v.length) : toBigInt(n);
  const count = requested > BigInt(v.length) ? BigInt(v.length) : requested;
  if (count === 0n) return 1n;
  // Sage also returns the first entry for a negative order; an empty input
  // reaches the original list-index error instead of the zero-order case.
  if (v.length === 0) throw new IndexError('list index out of range');
  let previous = 1n,
    current = toBigInt(v[0]!);
  if (count <= 1n) return current;
  for (let i = 1; i < Number(count); i++) {
    const next = previous + current * toBigInt(v[i]!);
    previous = current;
    current = next;
  }
  return current;
}

/**
 * Return 1 if ax^2 + by^2 p-adically represents a nonzero square,
 * otherwise returns -1. If either a or b is 0, returns 0.
 *
 * @param a - Integer or rational
 * @param b - Integer or rational
 * @param p - Prime or -1 (representing the archimedean place)
 * @param algorithm - Algorithm to use ('pari', 'direct', or 'all')
 * @returns 0, -1, or 1
 *
 * @see Reference: sage/arith/misc.py:hilbert_symbol
 * @see Deviation: Hilbert symbol dependency domain
 */
export function hilbert_symbol(
  a: RationalLike,
  b: RationalLike,
  p: IntegerLike,
  algorithm: 'pari' | 'direct' | 'all' = 'pari'
): bigint {
  const prime = toBigInt(p);
  if (prime !== -1n && !is_prime(prime)) throw new ValueError('p must be prime or -1');
  const aq = toRational(a),
    bq = toRational(b);
  let x = aq.numerator * aq.denominator,
    y = bq.numerator * bq.denominator;
  if (algorithm === 'pari') return BigInt(pari_hilbert(x, y, prime === -1n ? 0n : prime));
  if (algorithm === 'all') {
    const pari = hilbert_symbol(x, y, prime, 'pari'),
      direct = hilbert_symbol(x, y, prime, 'direct');
    if (pari !== direct)
      throw new RuntimeError(
        `there is a bug in hilbert_symbol; two ways of computing the Hilbert symbol (${x},${y})_${prime} disagree`
      );
    return pari;
  }
  if (algorithm !== 'direct')
    throw new ValueError(`algorithm ${algorithm === null ? 'None' : algorithm} not defined`);
  if (x === 0n || y === 0n) return 0n;
  if (prime !== -1n) {
    const square = prime * prime;
    while (x % square === 0n) x /= square;
    while (y % square === 0n) y /= square;
  }
  if (prime !== 2n && [x, y, x + y].some((v) => kronecker_symbol(v, prime) === 1n)) return 1n;
  if (x % prime === 0n) {
    if (y % prime === 0n)
      return hilbert_symbol(prime, -(y / prime), prime) * hilbert_symbol(x / prime, y, prime);
    if (prime === 2n && mod(y, 4n) === 3n) {
      if (kronecker_symbol(x + y, prime) === -1n) return -1n;
    } else if (kronecker_symbol(y, prime) === -1n) return -1n;
  } else if (y % prime === 0n) {
    if (prime === 2n && mod(x, 4n) === 3n) {
      if (kronecker_symbol(x + y, prime) === -1n) return -1n;
    } else if (kronecker_symbol(x, prime) === -1n) return -1n;
  } else if (prime === 2n && mod(x, 4n) === 3n && mod(y, 4n) === 3n) return -1n;
  return 1n;
}

/**
 * Return the product of all (finite) primes where the Hilbert symbol is -1.
 *
 * @param a - Integer
 * @param b - Integer
 * @returns Squarefree positive integer
 *
 * @see Reference: sage/arith/misc.py:hilbert_conductor
 */
export function hilbert_conductor(a: IntegerLike, b: IntegerLike): bigint {
  a = toBigInt(a);
  b = toBigInt(b);
  // Reference: sage/arith/misc.py:hilbert_conductor
  // Return the product of all finite primes where the Hilbert symbol is -1

  // Get prime divisors of a and b, plus 2
  const primesSet = new Set<bigint>([2n]);

  for (const p of prime_factors(a)) {
    primesSet.add(p);
  }
  for (const p of prime_factors(b)) {
    primesSet.add(p);
  }

  let product = 1n;
  for (const p of primesSet) {
    if (hilbert_symbol(a, b, p) === -1n) {
      product *= p;
    }
  }

  return product;
}

/**
 * Find a pair of integers (a, b) such that hilbert_conductor(a, b) == d.
 *
 * @param d - Square-free positive integer
 * @returns Pair of integers
 *
 * @see Reference: sage/arith/misc.py:hilbert_conductor_inverse
 */
export function hilbert_conductor_inverse(d: IntegerLike): [bigint, bigint] {
  d = toBigInt(d);
  // Reference: sage/arith/misc.py:hilbert_conductor_inverse
  // Find (a, b) such that hilbert_conductor(a, b) == d

  if (d <= 0n) {
    throw new ValueError('d needs to be positive');
  }

  if (d === 1n) {
    return [-1n, 1n];
  }

  if (d === 2n) {
    return [-1n, -1n];
  }

  if (is_prime(d)) {
    if (mod(d, 4n) === 3n) {
      return [-1n, -d];
    }
    if (mod(d, 8n) === 5n) {
      return [-2n, -d];
    }
    // Find q such that q ≡ 3 (mod 4) and kronecker(d, q) = -1
    let q = 3n;
    while (mod(q, 4n) !== 3n || kronecker_symbol(d, q) !== -1n) {
      q = next_prime(q);
    }
    return [-q, -d];
  }

  // Check if d is squarefree
  const mo = moebius(d);
  if (mo === 0n) {
    throw new ValueError('d needs to be squarefree');
  }

  // d is composite and squarefree
  let dd: bigint;
  if (d % 2n === 0n && mod(mo * d, 16n) !== 2n) {
    dd = (mo * d) / 2n;
  } else {
    dd = mo * d;
  }

  let q = 1n;
  while (hilbert_conductor(-q, dd) !== d) {
    q += 1n;
  }

  if (dd % q === 0n) {
    dd = dd / q;
  }

  return [-q, dd];
}

/**
 * Return the falling factorial (x)_a = x(x-1)...(x-a+1).
 *
 * @param x - Element of a ring
 * @param a - Non-negative integer
 * @returns The falling factorial
 *
 * @example
 * ```typescript
 * falling_factorial(10n, 3n)  // 720n
 * ```
 *
 * @see Reference: sage/arith/misc.py:falling_factorial
 */
export function falling_factorial(x: bigint, a: bigint): bigint {
  // For integer a >= 0: x(x-1)(x-2)...(x-a+1)
  if (a < 0n) {
    throw new ValueError('a must be a non-negative integer');
  }

  if (a === 0n) {
    return 1n;
  }

  let result = 1n;
  for (let i = 0n; i < a; i++) {
    result *= x - i;
  }
  return result;
}

/**
 * Return the rising factorial (x)^a = x(x+1)...(x+a-1).
 *
 * Also known as the Pochhammer symbol.
 *
 * @param x - Element of a ring
 * @param a - Non-negative integer
 * @returns The rising factorial
 *
 * @example
 * ```typescript
 * rising_factorial(10n, 3n)  // 1320n
 * ```
 *
 * @see Reference: sage/arith/misc.py:rising_factorial
 */
export function rising_factorial(x: bigint, a: bigint): bigint {
  // For integer a >= 0: x(x+1)(x+2)...(x+a-1)
  // Also known as the Pochhammer symbol
  if (a < 0n) {
    throw new ValueError('a must be a non-negative integer');
  }

  if (a === 0n) {
    return 1n;
  }

  let result = 1n;
  for (let i = 0n; i < a; i++) {
    result *= x + i;
  }
  return result;
}

/**
 * Return the ceiling of x (smallest integer >= x).
 *
 * @param x - A number
 * @returns The ceiling of x
 *
 * @see Reference: sage/arith/misc.py:integer_ceil
 */
export function integer_ceil(x: FloatInput): bigint {
  return integer_rounding(x, 'ceil');
}

/** Method dispatch followed by Sage's math.floor/ceil(float(x)) fallback. */
function integer_rounding(x: FloatInput, operation: 'floor' | 'ceil'): bigint {
  if (typeof x === 'bigint') return x;
  try {
    if (x !== null && typeof x === 'object' && operation in x) {
      const method = (x as { floor?: () => IntegerLike; ceil?: () => IntegerLike })[operation];
      if (typeof method === 'function') return toBigInt(method.call(x));
    }
  } catch (error) {
    // Sage falls back when the method itself raises AttributeError as well.
    if (!(error instanceof AttributeError)) throw error;
  }
  let value: number;
  try {
    value = python_float(x);
  } catch (error) {
    if (!(error instanceof TypeError)) throw error;
    throw new NotImplementedError(
      `computation of ${operation} of ${x === null ? 'None' : String(x)} not implemented`
    );
  }
  if (Number.isNaN(value)) throw new ValueError('cannot convert float NaN to integer');
  if (!Number.isFinite(value)) throw new OverflowError('cannot convert float infinity to integer');
  return BigInt(operation === 'floor' ? Math.floor(value) : Math.ceil(value));
}

/**
 * Return the largest integer <= x.
 *
 * @param x - A number
 * @returns The floor of x
 *
 * @see Reference: sage/arith/misc.py:integer_floor
 */
export function integer_floor(x: FloatInput): bigint {
  return integer_rounding(x, 'floor');
}

/**
 * Truncate to the integer closer to zero.
 *
 * @param x - A number
 * @returns The truncation of x
 *
 * @see Reference: sage/arith/misc.py:integer_trunc
 */
export function integer_trunc(x: FloatInput): bigint {
  // The source compares with zero before invoking floor/ceil; strings and bytes
  // therefore differ from the coercion accepted by integer_floor/integer_ceil.
  if (x === null || typeof x === 'string' || x instanceof Uint8Array) {
    const name = x === null ? 'NoneType' : typeof x === 'string' ? 'str' : 'bytes';
    throw new TypeError(`'>=' not supported between instances of '${name}' and 'int'`);
  }
  let nonnegative: boolean;
  if (typeof x === 'number' || typeof x === 'bigint' || typeof x === 'boolean') {
    nonnegative = typeof x === 'boolean' || x >= 0;
  } else if (x instanceof Rational) {
    nonnegative = x.sign >= 0n;
  } else if ('sign' in x) {
    nonnegative = !('is_NaN' in x && x.is_NaN()) && x.sign() >= 0;
  } else {
    nonnegative = x.value >= 0n;
  }
  return nonnegative ? integer_floor(x) : integer_ceil(x);
}

/** 2^32 — SageMath switches to `sage.rings.sum_of_squares` below this bound. */
const SUM_OF_SQUARES_CUTOFF = 4294967296n;

/**
 * Write the integer n as a sum of two integer squares if possible;
 * otherwise raise ValueError.
 *
 * A number n can be written as a sum of two squares if and only if
 * all prime powers p^e in its factorization with p ≡ 3 (mod 4) have e even.
 *
 * For n < 2^32 SageMath dispatches to `sum_of_squares.two_squares_pyx`, whose
 * answer generally differs from the Cornacchia-based one, so we dispatch too.
 *
 * @param n - Integer
 * @returns Tuple [a, b] with a <= b such that n = a^2 + b^2
 *
 * @example
 * ```typescript
 * two_squares(389n)  // [10n, 17n]
 * two_squares(21n)   // ValueError: 21 is not a sum of 2 squares
 * two_squares(0n)    // [0n, 0n]
 * ```
 *
 * @see Reference: sage/arith/misc.py:two_squares
 * @see Deviation: Square decomposition input mapping and native bounds
 */
export function two_squares(input: IntegerLike): [bigint, bigint] {
  const n = toBigInt(input);
  if (n < 0n) {
    throw new ValueError(`${n} is not a sum of 2 squares`);
  }

  if (n === 0n) {
    return [0n, 0n];
  }

  if (n < SUM_OF_SQUARES_CUTOFF) {
    return two_squares_pyx(n);
  }

  // Factor n
  const F = factor(n);

  // Check whether it is possible to write n as a sum of two squares:
  // All prime powers p^e must have p = 2 or p ≡ 1 (mod 4) or e even.
  for (const [p, e] of F) {
    if (p === -1n) continue;
    if (e % 2n === 1n && p % 4n === 3n) {
      throw new ValueError(`${n} is not a sum of 2 squares`);
    }
  }

  // We run over all factors of n, write each factor p^e as a sum of 2 squares
  // and accumulate the product (using multiplication in Z[i]) in a^2 + b^2.
  let a = 1n;
  let b = 0n;

  for (const [p, e] of F) {
    if (p === -1n) continue;

    // Handle even powers: just multiply by p^(e/2)
    if (e >= 2n) {
      const m = p ** (e / 2n);
      a *= m;
      b *= m;
    }

    // Handle odd power
    if (e % 2n === 1n) {
      if (p === 2n) {
        // 2 = 1^2 + 1^2, so (a + bi) *= (1 + i)
        const newA = a - b;
        const newB = a + b;
        a = newA;
        b = newB;
      } else {
        // p ≡ 1 (mod 4): use Cornacchia's algorithm
        // Find a square root of -1 mod p
        let s = findSqrtMinusOne(p);

        // Apply Cornacchia's algorithm to write p as r^2 + s^2
        let r = p;
        while (s * s > p) {
          const temp = s;
          s = r % s;
          r = temp;
        }
        r = r % s;

        // Multiply (a + bi) by (r + si)
        const newA = a * r - b * s;
        const newB = b * r + a * s;
        a = newA;
        b = newB;
      }
    }
  }

  // Take absolute values
  if (a < 0n) a = -a;
  if (b < 0n) b = -b;

  // Return sorted
  return a <= b ? [a, b] : [b, a];
}

/**
 * Find a square root of -1 modulo p, where p ≡ 1 (mod 4).
 * Uses the fact that if y is a quadratic non-residue, then y^((p-1)/4) is a sqrt of -1.
 */
function findSqrtMinusOne(p: bigint): bigint {
  // Try small values starting from 2
  for (let y = 2n; y < p; y++) {
    const s = power_mod(y, (p - 1n) / 4n, p);
    // Check if s^2 ≡ -1 (mod p)
    if ((s * s + 1n) % p === 0n) {
      return s;
    }
  }
  // Should never reach here for p ≡ 1 (mod 4)
  throw new ValueError('Could not find sqrt(-1) mod p');
}

/**
 * Write the integer n as a sum of three integer squares if possible;
 * otherwise raise ValueError.
 *
 * By Legendre's three-square theorem, a positive integer n can be expressed
 * as a sum of three squares if and only if n is NOT of the form 4^a(8b+7).
 *
 * @param n - Integer
 * @returns Tuple [a, b, c] with a <= b <= c such that n = a^2 + b^2 + c^2
 *
 * @example
 * ```typescript
 * three_squares(389n)  // [1n, 8n, 18n]
 * three_squares(7n)    // ValueError: 7 is not a sum of 3 squares
 * three_squares(0n)    // [0n, 0n, 0n]
 * ```
 *
 * @see Reference: sage/arith/misc.py:three_squares
 * @see Deviation: Square decomposition input mapping and native bounds
 */
export function three_squares(input: IntegerLike): [bigint, bigint, bigint] {
  const n = toBigInt(input);
  if (n < 0n) {
    throw new ValueError(`${n} is not a sum of 3 squares`);
  }

  if (n === 0n) {
    return [0n, 0n, 0n];
  }

  if (n < SUM_OF_SQUARES_CUTOFF) {
    return three_squares_pyx(n);
  }

  // First, remove all factors of 4 from n
  // e = valuation(n, 2) // 2
  // Integer.valuation(2) uses GMP's low-bit scan, followed by one shift.
  const e = BigInt((n & -n).toString(2).length - 1) / 2n;
  const N = n >> (2n * e);
  const m = 1n << e; // 2^e

  // Let x be the largest integer at most sqrt(N)
  const sqrtN = isqrt(N);

  // Check if N is a perfect square - special case to avoid factoring
  if (sqrtN * sqrtN === N) {
    return [0n, 0n, sqrtN * m];
  }

  // Check if N ≡ 7 (mod 8) - by Legendre's theorem, n is not a sum of 3 squares
  if (N % 8n === 7n) {
    throw new ValueError(`${n} is not a sum of 3 squares`);
  }

  // Find x such that N - x^2 can be written as sum of 2 squares.
  // We use the Rabin-Shallit approach: find x such that N - x^2 is either
  // p or 2p, where p is a prime ≡ 1 (mod 4). This makes two_squares fast.
  // Algorithm: https://schorn.ch/lagrange.html
  let x = sqrtN;

  if (N % 4n === 1n) {
    // Write N = x^2 + p with x even, p = 1 mod 4 prime
    // Note: when N ≡ 1 (mod 4) and x is even, x^2 ≡ 0 (mod 4),
    // so p = N - x^2 ≡ 1 (mod 4) automatically
    if (x % 2n === 1n) {
      x -= 1n;
    }
    while (x >= 0n) {
      const p = N - x * x;
      if (is_pseudoprime(p)) {
        break;
      }
      x -= 2n;
    }
  } else if (N % 4n === 2n) {
    // Write N = x^2 + p with x odd, p = 1 mod 4 prime
    // Note: when N ≡ 2 (mod 4) and x is odd, x^2 ≡ 1 (mod 4),
    // so p = N - x^2 ≡ 1 (mod 4) automatically
    if (x % 2n === 0n) {
      x -= 1n;
    }
    while (x >= 0n) {
      const p = N - x * x;
      if (is_pseudoprime(p)) {
        break;
      }
      x -= 2n;
    }
  } else if (N % 8n === 3n) {
    // Write N = x^2 + 2p with x odd, p = 1 mod 4 prime
    // Note: when N ≡ 3 (mod 8) and x is odd, x^2 ≡ 1 (mod 8),
    // so N - x^2 ≡ 2 (mod 8), and (N - x^2)/2 ≡ 1 (mod 4) automatically
    if (x % 2n === 0n) {
      x -= 1n;
    }
    while (x >= 0n) {
      const p = (N - x * x) >> 1n; // (N - x^2) / 2
      if (is_pseudoprime(p)) {
        break;
      }
      x -= 2n;
    }
  }

  // If we found no good x using the prime-finding approach, fall back to brute force.
  // This should only happen for small values of N (numerical experiments suggest
  // 9634 is the largest integer that may need brute force).
  if (x < 0n) {
    x = sqrtN;
  }

  // Try to write N - x^2 as sum of 2 squares.
  // In the usual case, this loop executes only once since we already found
  // the right x above. This only really loops in the brute force fallback case.
  while (true) {
    try {
      const [a, b] = two_squares(N - x * x);
      if (x >= b) return [a * m, b * m, x * m];
      if (x >= a) return [a * m, x * m, b * m];
      return [x * m, a * m, b * m];
    } catch (error) {
      if (!(error instanceof ValueError)) throw error;
      x -= 1n;
      if (x < 0n) throw new AssertionError('');
    }
  }
}

/**
 * Write the integer n as a sum of four integer squares.
 *
 * By Lagrange's four square theorem, every non-negative integer can be
 * expressed as a sum of four integer squares.
 *
 * @param n - Non-negative integer
 * @returns Tuple [a, b, c, d] with a <= b <= c <= d such that n = a^2 + b^2 + c^2 + d^2; raises ValueError if n is negative
 *
 * @example
 * ```typescript
 * four_squares(3n)    // [0n, 1n, 1n, 1n]
 * four_squares(13n)   // [0n, 0n, 2n, 3n]
 * four_squares(130n)  // [0n, 0n, 3n, 11n]
 * ```
 *
 * @see Reference: sage/arith/misc.py:four_squares
 * @see Deviation: Square decomposition input mapping and native bounds
 */
export function four_squares(input: IntegerLike): [bigint, bigint, bigint, bigint] {
  const n = toBigInt(input);
  if (n < 0n) {
    throw new ValueError(`${n} is not a sum of 4 squares`);
  }

  if (n > 0n && n < SUM_OF_SQUARES_CUTOFF) {
    return four_squares_pyx(n);
  }

  if (n === 0n) {
    return [0n, 0n, 0n, 0n];
  }

  // First, remove all factors of 4 from n
  // Integer.valuation(2) uses GMP's low-bit scan, followed by one shift.
  const e = BigInt((n & -n).toString(2).length - 1) / 2n;
  const N = n >> (2n * e);
  const m = 1n << e; // 2^e

  // Find x such that N - x^2 can be written as sum of 3 squares
  // N - x^2 must be 1, 2, 3, 5, or 6 mod 8 (not 0, 4, or 7)
  let x = isqrt(N);
  let y = N - x * x;

  // If y >= 7 and (y ≡ 0 (mod 4) or y ≡ 7 (mod 8)), adjust x
  if (y >= 7n && (y % 4n === 0n || y % 8n === 7n)) {
    x -= 1n;
    y = N - x * x;
  }

  const [a, b, c] = three_squares(y);
  return [a * m, b * m, c * m, x * m];
}

/**
 * Write the integer n as a sum of k integer squares if possible;
 * otherwise raise ValueError.
 *
 * For k >= 4, this always succeeds for non-negative n (Lagrange's theorem).
 * For k = 3, succeeds iff n is not of form 4^a(8b+7) (Legendre's theorem).
 * For k = 2, succeeds iff all primes p ≡ 3 (mod 4) in factorization have even exponent.
 * For k = 1, succeeds iff n is a perfect square.
 * For k = 0, succeeds iff n = 0.
 *
 * @param k - Non-negative integer
 * @param n - Integer
 * @returns Array [x_1, ..., x_k] of non-negative integers whose squares sum to n
 *
 * @example
 * ```typescript
 * sum_of_k_squares(2, 9634n)  // [15n, 97n]
 * sum_of_k_squares(4, 9634n)  // [1n, 2n, 5n, 98n]
 * sum_of_k_squares(1, 9n)     // [3n]
 * sum_of_k_squares(1, 10n)    // ValueError: 10 is not a sum of 1 square
 * sum_of_k_squares(0, 0n)     // []
 * ```
 *
 * @see Reference: sage/arith/misc.py:sum_of_k_squares
 * @see Deviation: Square decomposition input mapping and native bounds
 */
export function sum_of_k_squares(count: IntegerLike | number, input: IntegerLike): bigint[] {
  let n = toBigInt(input);
  let k = typeof count === 'number' ? integer_trunc(count) : toBigInt(count);
  if (k <= 4n) {
    if (k === 4n) return four_squares(n);
    if (k === 3n) return three_squares(n);
    if (k === 2n) return two_squares(n);
    if (k === 1n) {
      if (n >= 0n) {
        const x = isqrt(n);
        if (x * x === n) return [x];
      }
      throw new ValueError(`${n} is not a sum of 1 square`);
    }
    if (k === 0n) {
      if (n === 0n) return [];
      throw new ValueError(`${n} is not a sum of 0 squares`);
    }
    throw new ValueError(`k = ${k} must be nonnegative`);
  }
  if (n < 0n) throw new ValueError(`${n} is not a sum of ${k} squares`);
  const extras: bigint[] = [];
  while (k > 4n) {
    const x = isqrt(n);
    extras.push(x);
    n -= x * x;
    k--;
  }
  return [...four_squares(n), ...extras.reverse()];
}

/**
 * Return the subfactorial (derangement number) of n.
 *
 * The subfactorial !n is the number of permutations of n elements with no fixed points.
 * !n = n! * sum_{k=0}^n (-1)^k / k! = floor(n! / e + 1/2) for n >= 1
 *
 * @param n - Non-negative integer
 * @returns The subfactorial of n
 *
 * @example
 * ```typescript
 * subfactorial(0n)  // 1n
 * subfactorial(1n)  // 0n
 * subfactorial(2n)  // 1n
 * subfactorial(3n)  // 2n
 * subfactorial(4n)  // 9n
 * subfactorial(5n)  // 44n
 * ```
 *
 * @see Reference: sage/arith/misc.py:subfactorial
 */
export function subfactorial(n: IntegerLike): bigint {
  const value = toBigInt(n);
  if (value < 0n) {
    throw new ValueError('factorial -- must be nonnegative');
  }

  if (value === 0n) {
    return 1n;
  }

  if (value === 1n) {
    return 0n;
  }

  // Use recurrence: !n = (n-1) * (!(n-1) + !(n-2))
  let prev2 = 1n; // !0
  let prev1 = 0n; // !1

  for (let k = 2n; k <= value; k++) {
    const curr = (k - 1n) * (prev1 + prev2);
    prev2 = prev1;
    prev1 = curr;
  }

  return prev1;
}

/**
 * Return whether n is a power of 2.
 *
 * @param n - Integer
 * @returns true if n is a power of 2
 *
 * @see Reference: sage/arith/misc.py:is_power_of_two
 */
export function is_power_of_two(n: IntegerLike): boolean {
  const value = toBigInt(n);
  // A number is a power of 2 if it has exactly one bit set (popcount == 1)
  // This is equivalent to: n > 0 && (n & (n - 1)) === 0
  return value > 0n && (value & (value - 1n)) === 0n;
}

/**
 * Return the n successive differences of the elements in lis.
 *
 * @param lis - List of numbers
 * @param n - Number of differences (default: 1)
 * @returns List of differences
 *
 * @example
 * ```typescript
 * differences([1n, 4n, 9n, 16n])     // [3n, 5n, 7n]
 * differences([1n, 4n, 9n, 16n], 2n) // [2n, 2n]
 * ```
 *
 * @see Deviation: Successive differences recursion limit
 * @see Reference: sage/arith/misc.py:differences
 */
export function differences(lis: IntegerLike[], n: IntegerLike = 1n): bigint[] {
  const order = toBigInt(n);
  if (order < 1n) {
    throw new ValueError('n must be greater than 0');
  }
  let result = lis.map(toBigInt);
  // Iteration avoids Python's recursion limit; once empty all further differences are empty.
  for (let i = 0n; i < order && result.length > 0; i++) {
    const next: bigint[] = [];
    for (let j = 0; j < result.length - 1; j++) {
      next.push(result[j + 1]! - result[j]!);
    }
    result = next;
  }
  return result;
}

type DisplayComplex = { re: number; im: number };
/** Binary64 input mapping of sage/arith/misc.py:_key_complex_for_display.
 * Nine display digits correspond to 34 binary bits. Keep MPFR keys so rounding
 * a finite binary64 maximum cannot collapse it into infinity.
 */
function _key_complex_for_display(a: DisplayComplex): [number, mpfr_t, number] {
  const ar = a.re,
    ai = a.im;
  const real = mpfr_init2(ai === 0 ? 53 : 34);
  mpfr_set_d(real, ai !== 0 && Math.abs(ar) < 1e-10 ? 0 : ar);
  return [ai === 0 ? 0 : 1, real, ai];
}
/** Stable display order, with original tuple and empty-list identity behavior.
 * @see Reference: sage/arith/misc.py:sort_complex_numbers_for_display
 * @see Deviation: Complex display ordering and Python sorting
 */
export function sort_complex_numbers_for_display<T extends DisplayComplex>(nums: T[]): T[];
export function sort_complex_numbers_for_display<T extends readonly [DisplayComplex, ...unknown[]]>(
  nums: T[]
): T[];
export function sort_complex_numbers_for_display<
  T extends DisplayComplex | readonly [DisplayComplex, ...unknown[]],
>(nums: T[]): T[] {
  if (nums.length === 0) return nums;
  const tuples = Array.isArray(nums[0]);
  const entries = nums.map((value) => ({
    value,
    key: _key_complex_for_display(
      tuples ? (value as readonly [DisplayComplex])[0] : (value as DisplayComplex)
    ),
  }));
  const ordered = python_sorted(entries, (a, b) => {
    if (a.key[0] !== b.key[0]) return a.key[0] < b.key[0];
    if (a.key[1].kind === 'nan' || b.key[1].kind === 'nan') return false;
    const comparison = mpfr_cmp(a.key[1], b.key[1]);
    return comparison !== 0 ? comparison < 0 : a.key[2] < b.key[2];
  });
  return ordered.map((entry) => entry.value);
}

/**
 * Return the fundamental discriminant of Q(sqrt(D)).
 *
 * The fundamental discriminant of Q(sqrt(D)) is the discriminant of its ring
 * of integers. If D is squarefree, this is D if D ≡ 1 (mod 4), and 4D otherwise.
 *
 * @param D - Integer
 * @returns The fundamental discriminant
 *
 * @example
 * ```typescript
 * fundamental_discriminant(5n)   // 5n (5 ≡ 1 mod 4)
 * fundamental_discriminant(2n)   // 8n (2 ≡ 2 mod 4, so 4*2)
 * fundamental_discriminant(-3n)  // -3n (-3 ≡ 1 mod 4)
 * fundamental_discriminant(-7n)  // -7n (-7 ≡ 1 mod 4)
 * ```
 *
 * @see Reference: sage/arith/misc.py:fundamental_discriminant
 */
export function fundamental_discriminant(D: IntegerLike): bigint {
  D = toBigInt(D);

  // First get the squarefree part
  const sf = squarefree_part(D);

  // The fundamental discriminant is:
  // - sf if sf ≡ 1 (mod 4)
  // - 4*sf otherwise
  const mod4 = ((sf % 4n) + 4n) % 4n;
  if (mod4 === 1n) {
    return sf;
  }
  return 4n * sf;
}

/**
 * Return an iterator over the squarefree divisors of x.
 *
 * A squarefree divisor is a divisor that is not divisible by any perfect square > 1.
 *
 * @param x - An integer (nonzero)
 * @returns Iterator over squarefree divisors
 *
 * @example
 * ```typescript
 * [...squarefree_divisors(12n)]  // [1n, 2n, 3n, 6n]
 * ```
 *
 * @see Reference: sage/arith/misc.py:squarefree_divisors
 */
export function* squarefree_divisors(x: IntegerLike): Generator<bigint, void, unknown> {
  // Sage's powerset grows its mask width as primes are consumed. BigInt
  // preserves that order without the JavaScript 32-bit bitwise cutoff.
  const primes = prime_factors(toBigInt(x));
  yield 1n;
  const pairs: [bigint, bigint][] = [];
  let power2 = 1n;
  for (const p of primes) {
    pairs.push([power2, p]);
    const nextPower2 = power2 << 1n;
    for (let mask = power2; mask < nextPower2; mask++) {
      let divisor = 1n;
      for (const [bit, prime] of pairs) if (mask & bit) divisor *= prime;
      yield divisor;
    }
    power2 = nextPower2;
  }
}

/**
 * Return the Dedekind sum s(p, q).
 *
 * @param p - Integer
 * @param q - Integer
 * @param algorithm - Algorithm to use ('default', 'flint', or 'pari')
 * @returns The Dedekind sum as a rational
 *
 * @see Reference: sage/arith/misc.py:dedekind_sum
 * @see Deviation: Dedekind sum backend arithmetic
 */
export function dedekind_sum(
  p: IntegerLike,
  q: IntegerLike,
  algorithm: 'default' | 'flint' | 'pari' = 'default'
): { numerator: bigint; denominator: bigint } {
  // Select before coercion, as Sage does, including unknown-algorithm errors.
  if (algorithm !== 'default' && algorithm !== 'flint' && algorithm !== 'pari') {
    throw new ValueError('unknown algorithm');
  }
  const h = toBigInt(p);
  const k = toBigInt(q);
  const [numerator, denominator] =
    algorithm === 'pari' ? pari_sumdedekind(h, k) : fmpq_dedekind_sum(h, k);
  return { numerator, denominator };
}

/**
 * Interface representing a finite field element.
 */
interface FiniteFieldElement {
  trace(): { lift(): number };
  mul?(other: this): this;
}

/**
 * Interface representing a finite field.
 */
interface FiniteField {
  cardinality(): bigint;
  characteristic(): bigint;
  multiplicative_generator(): FiniteFieldElement;
  one(): FiniteFieldElement;
}

/**
 * Interface for character values with ring operations.
 */
interface CharacterValue {
  parent(): { zero(): CharacterValue; zeta(n: bigint): { powers(m: number): CharacterValue[] } };
  add?(other: CharacterValue): CharacterValue;
  mul?(other: CharacterValue): CharacterValue;
}

/**
 * Return the Gauss sum for a general finite field.
 *
 * For a finite field F of characteristic p, the Gauss sum associated
 * to a multiplicative character chi (with values in a ring K) is defined as:
 *
 *   sum_{x in F^*} chi(x) * zeta_p^{Tr(x)}
 *
 * where zeta_p in K is a primitive p-th root of unity and Tr is the
 * trace map from F to its prime field GF(p).
 *
 * @param char_value - Value of multiplicative character on the generator
 * @param finite_field - A finite field
 * @returns The Gauss sum (in the same ring as char_value)
 *
 * @example
 * ```typescript
 * // For GF(5), computing g(chi) where chi(g) = zeta_4
 * const F = GF(5n);
 * const zq = ComplexField.zeta(4);  // primitive 4th root of unity
 * const g = gauss_sum(zq, F);
 * // g * conjugate(g) should equal 5
 * ```
 *
 * @see Reference: sage/arith/misc.py:gauss_sum
 * @see Deviation: Arithmetic Functions Not Delegated to PARI/FLINT
 */
export function gauss_sum(char_value: CharacterValue, finite_field: FiniteField): CharacterValue {
  // Validate that finite_field is actually a finite field
  if (
    typeof finite_field.cardinality !== 'function' ||
    typeof finite_field.characteristic !== 'function'
  ) {
    throw new ValueError('second input must be a finite field');
  }

  const ring = char_value.parent();
  const q = finite_field.cardinality();
  const p = finite_field.characteristic();
  const gen = finite_field.multiplicative_generator();

  // Get powers of zeta_p for the additive character
  const zeta_p_powers = ring.zeta(p).powers(Number(p));
  const zeta_q = char_value;

  // Compute the Gauss sum:
  // sum_{k=0}^{q-2} chi(g^k) * zeta_p^{Tr(g^k)}
  // where chi(g^k) = zeta_q^k

  let resu = ring.zero();
  let gen_power = finite_field.one();
  let zq_power: CharacterValue = ring.zeta(1n).powers(1)[0]!; // Start with 1

  // We need to track powers more carefully
  // In the reference implementation:
  // gen_power = finite_field.one()
  // zq_power = ring.one()
  // for k in range(q - 1):
  //     resu += zq_power * zeta_p_powers[gen_power.trace().lift()]
  //     gen_power *= gen
  //     zq_power *= zeta_q

  // Since we can't easily do this without a proper ring implementation,
  // we'll provide a simplified version that works with numeric types
  // For a full implementation, you'd need proper cyclotomic field support

  for (let k = 0n; k < q - 1n; k++) {
    // Get the trace of gen^k
    const traceVal = gen_power.trace().lift();
    const zeta_p_power = zeta_p_powers[traceVal % Number(p)]!;

    // Add zq_power * zeta_p_power to result
    if (resu.add && zq_power.mul) {
      const term = zq_power.mul(zeta_p_power as unknown as CharacterValue);
      resu = resu.add(term);
    }

    // Update powers
    if (gen_power.mul) {
      gen_power = gen_power.mul(gen);
    }
    if (zq_power.mul) {
      zq_power = zq_power.mul(zeta_q);
    }
  }

  return resu;
}

/**
 * Return the value of the Dedekind psi function at N.
 *
 * psi(n) = n * product_{p|n, p prime}(1 + 1/p)
 *
 * @param N - Positive integer
 * @returns The Dedekind psi value
 *
 * @example
 * ```typescript
 * dedekind_psi(6n)   // 12n (6 * (1 + 1/2) * (1 + 1/3) = 6 * 3/2 * 4/3 = 12)
 * dedekind_psi(12n)  // 24n
 * ```
 *
 * @see Reference: sage/arith/misc.py:dedekind_psi
 */
export function dedekind_psi(N: IntegerLike): bigint {
  const value = toBigInt(N);
  const primes = prime_factors(value);
  let numerator = value;
  let denominator = 1n;
  for (const p of primes) {
    numerator *= p + 1n;
    denominator *= p;
  }
  return numerator / denominator;
}

/**
 * Return the largest divisor of x that is smooth over the factor base,
 * as a factorization over `base`.
 *
 * The entries of `base` need not be prime: SageMath divides out the full
 * valuation at each element of the factor base in turn, exactly as done here,
 * so `smooth_part(240, [6])` is `6^1` and not `2^4 * 3`.
 *
 * @param x - Element of a Euclidean domain
 * @param base - Factor base (sequence of elements)
 * @returns Factorization of the smooth part (unit part is always 1)
 *
 * @example
 * ```typescript
 * smooth_part(240n, [2n, 3n])  // [[2n, 4n], [3n, 1n]]
 * smooth_part(240n, [6n])      // [[6n, 1n]]
 * ```
 *
 * @see Deviation: Product trees and factor-base nontermination
 * @see Reference: sage/arith/misc.py:smooth_part
 */
export function smooth_part(
  x: IntegerLike,
  base: Iterable<IntegerLike> | ProductTree
): Factorization {
  const tree = base instanceof ProductTree ? base : new ProductTree(base);
  let remaining = toBigInt(x);
  const remainders = tree.remainders(remaining);
  const factors: Factorization = [];
  let index = 0;
  for (const p of tree) {
    if (remainders[index++] !== 0n) continue;
    remaining = floorQuotient(remaining, p);
    // The original loops forever for units or a zero running quotient.
    if (p === 1n || p === -1n || remaining === 0n) {
      throw new NotImplementedError(
        'SAGE_NOT_IMPLEMENTED: smooth_part: original does not terminate for this factor base'
      );
    }
    let exponent = 1n;
    while (remaining % p === 0n) {
      remaining /= p;
      exponent++;
    }
    factors.push([p, exponent]);
  }
  // Factorization(fs) sorts by the integer factor, then combines equal factors.
  factors.sort(([a], [b]) => (a < b ? -1 : a > b ? 1 : 0));
  const result: Factorization = [];
  for (const [p, exponent] of factors) {
    const previous = result[result.length - 1];
    if (previous?.[0] === p) previous[1] += exponent;
    else result.push([p, exponent]);
  }
  return result;
}

/** Integer floor division with a nonzero divisor; factors can be signed or overlap. */
function floorQuotient(value: bigint, divisor: bigint): bigint {
  const quotient = value / divisor;
  return value % divisor !== 0n && value < 0n !== divisor < 0n ? quotient - 1n : quotient;
}

/**
 * Return the largest divisor of x that is not divisible by any element of base.
 *
 * ALGORITHM: divide x by {@link smooth_part}, exactly as SageMath does. Note
 * that when the elements of `base` are composite this is *not* the same as
 * removing every prime factor of every base element: `coprime_part(240, [6])`
 * is 40 (only one factor 6 is stripped), not 5.
 *
 * @param x - Element of a Euclidean domain
 * @param base - Factor base
 * @returns The floor quotient by the smooth part (overlapping bases retain Sage's quirks)
 *
 * @example
 * ```typescript
 * coprime_part(240n, [2n, 3n])  // 5n
 * coprime_part(240n, [6n])      // 40n
 * ```
 *
 * @see Deviation: Product trees and factor-base nontermination
 * @see Reference: sage/arith/misc.py:coprime_part
 */
export function coprime_part(x: IntegerLike, base: Iterable<IntegerLike> | ProductTree): bigint {
  const value = toBigInt(x);
  let smooth = 1n;
  for (const [p, e] of smooth_part(value, base)) smooth *= p ** e;
  return floorQuotient(value, smooth);
}

/**
 * Return the Carmichael function of a positive integer n.
 *
 * The Carmichael function lambda(n) is the smallest positive integer k
 * such that a^k ≡ 1 (mod n) for all a coprime to n.
 *
 * @param n - Positive integer
 * @returns lambda(n)
 *
 * @example
 * ```typescript
 * carmichael_lambda(12n)  // 2n
 * ```
 *
 * @see Reference: sage/arith/misc.py:carmichael_lambda
 */
export function carmichael_lambda(n: IntegerLike): bigint {
  const value = toBigInt(n);
  // The Carmichael function lambda(n) is the smallest positive integer k
  // such that a^k ≡ 1 (mod n) for all a coprime to n.
  //
  // Algorithm:
  // - If n = 2 or 4: lambda(n) = phi(n)
  // - If n = 2^k with k >= 3: lambda(2^k) = 2^(k-2)
  // - If n = p^k for odd prime p: lambda(p^k) = phi(p^k) = p^(k-1) * (p-1)
  // - For composite n = p1^k1 * p2^k2 * ...: lambda(n) = lcm(lambda(p1^k1), lambda(p2^k2), ...)

  if (value < 1n) {
    throw new ValueError('Input n must be a positive integer.');
  }

  if (value === 1n) {
    return 1n;
  }

  const factors = factor(value);
  const lambdaValues: bigint[] = [];

  for (const [p, k] of factors) {
    if (p === -1n) {
      continue; // Skip the sign factor
    }

    if (p === 2n) {
      // Special case for powers of 2
      // lambda(2) = 1, lambda(4) = 2, lambda(2^k) = 2^(k-2) for k >= 3
      if (k === 1n) {
        lambdaValues.push(1n);
      } else if (k === 2n) {
        lambdaValues.push(2n);
      } else {
        // k >= 3
        lambdaValues.push(1n << (k - 2n));
      }
    } else {
      // Odd prime p: lambda(p^k) = phi(p^k) = p^(k-1) * (p - 1)
      lambdaValues.push(p ** (k - 1n) * (p - 1n));
    }
  }

  if (lambdaValues.length === 0) {
    return 1n;
  }

  // Return lcm of all lambda values
  return lcm(lambdaValues);
}

/**
 * Return the odd part of the integer n.
 *
 * This is n / 2^v, where v = valuation(n, 2).
 *
 * @param n - Integer
 * @returns The odd part of n
 *
 * @example
 * ```typescript
 * odd_part(24n)  // 3n
 * odd_part(5n)   // 5n
 * ```
 *
 * @see Reference: sage/arith/misc.py:odd_part
 */
export function odd_part(n: IntegerLike): bigint {
  const value = toBigInt(n);
  if (value === 0n) return 0n;
  // GMP mpz_scan1 followed by mpz_tdiv_q_2exp. The isolated low bit
  // identifies the exact shift, including negative integers.
  const bits = (value & -value).toString(2).length - 1;
  return value >> BigInt(bits);
}

/**
 * Return the prime-to-m part of n, i.e. the largest divisor of n coprime to m.
 *
 * Transcribed from SageMath's `Integer.prime_to_m_part`
 * (`sage/rings/integer.pyx:3033`): n = 0 is an error, m = 0 gives 1, and a unit
 * m (+-1) gives n back unchanged. The sign of n is preserved.
 *
 * @param n - Nonzero integer
 * @param m - Integer
 * @returns The prime-to-m part of n
 * @throws {ArithmeticError} `self must be nonzero` when n = 0
 *
 * @example
 * ```typescript
 * prime_to_m_part(43434n, 20n)  // 21717n
 * prime_to_m_part(2048n, 2n)    // 1n
 * prime_to_m_part(2048n, 3n)    // 2048n
 * prime_to_m_part(240n, 0n)     // 1n
 * ```
 *
 * @see Reference: sage/rings/integer.pyx:prime_to_m_part
 */
export function prime_to_m_part(n: IntegerLike, m: IntegerLike): bigint {
  n = toBigInt(n);
  if (n === 0n) {
    throw new ArithmeticError('self must be nonzero');
  }
  m = toBigInt(m);

  if (m === 0n) {
    return 1n;
  }

  if (m === 1n || m === -1n) {
    return n;
  }

  // mpz_gcd is nonnegative, and mpz_divexact preserves the sign of n.
  let g = gcd(n, m);
  if (g === 1n) {
    return n;
  }

  let result = n / g;
  while (g !== 1n) {
    g = gcd(result, g);
    result /= g;
  }

  return result;
}

/**
 * Return the native-width gcd implementation selected by the size bound.
 * @see Reference: sage/arith/misc.py:get_gcd
 * @see Deviation: Bounded native arithmetic and Cython error sentinels
 */
export function get_gcd(order: IntegerLike): (a: IntegerLike, b: IntegerLike) => bigint {
  const size = toBigInt(order);
  if (size <= 46340n) {
    const arith = new arith_int();
    return arith.gcd_int.bind(arith);
  }
  if (size <= 2147483647n) {
    const arith = new arith_llong();
    return arith.gcd_longlong.bind(arith);
  }
  return gcd;
}

/**
 * Return the native-width inverse implementation selected by the size bound.
 * @see Reference: sage/arith/misc.py:get_inverse_mod
 * @see Deviation: Bounded native arithmetic and Cython error sentinels
 */
export function get_inverse_mod(order: IntegerLike): (a: IntegerLike, m: IntegerLike) => bigint {
  const size = toBigInt(order);
  if (size <= 46340n) {
    const arith = new arith_int();
    return arith.inverse_mod_int.bind(arith);
  }
  if (size <= 2147483647n) {
    const arith = new arith_llong();
    return arith.inverse_mod_longlong.bind(arith);
  }
  return inverse_mod;
}

/**
 * Maximal Quotient Rational Reconstruction.
 *
 * @param u - Integer with m > u >= 0
 * @param m - Modulus
 * @param T - Positive integer bound
 * @returns (n, d) or null
 *
 * @see Reference: sage/arith/misc.py:mqrr_rational_reconstruction
 * @see Deviation: Maximal-quotient reconstruction and exact source division
 */
export function mqrr_rational_reconstruction(
  input: IntegerLike,
  modulus: IntegerLike,
  threshold: IntegerLike
): [bigint, bigint] | null {
  const u = toBigInt(input);
  const m = toBigInt(modulus);
  const T = toBigInt(threshold);
  if (u === 0n) return m > T ? [0n, 1n] : null;

  let n = new Rational(0n);
  let d = new Rational(0n);
  let t0 = new Rational(0n);
  let r0 = new Rational(m);
  let t1 = new Rational(1n);
  let r1 = new Rational(u);
  let bound = new Rational(T);
  while (!r1.eq(0n) && r0.gt(bound)) {
    // The original Python body uses / on Sage integers: this is exact division,
    // despite its stale "C division implicit floor" comment. Preserve that behavior.
    const q = r0.div(r1);
    if (q.gt(bound)) {
      n = r1;
      d = t1;
      bound = q;
    }
    [r0, r1] = [r1, r0.sub(q.mul(r1))];
    [t0, t1] = [t1, t0.sub(q.mul(t1))];
  }
  // Integer inputs leave r1 zero after the first exact division, so any chosen
  // n and d are integers. Do not normalize the returned pair's signs.
  if (!d.eq(0n) && gcd(n.numerator, d.numerator) === 1n) {
    return [n.numerator, d.numerator];
  }
  return null;
}

/**
 * Alias for kronecker_symbol.
 * @see kronecker_symbol
 */
export const kronecker = kronecker_symbol;

/**
 * Alias for CRT_list.
 * @see crt
 */
export const CRT = crt;

// Fibonacci-related functions (from other SageMath modules but commonly used)

/**
 * Return the n-th Fibonacci number.
 *
 * @param n - Non-negative integer
 * @returns F_n
 *
 * @see Reference: sage/rings/integer.pyx:fibonacci
 */
export function fibonacci(n: bigint): bigint {
  // Handle negative n: F_{-n} = (-1)^{n+1} * F_n
  if (n < 0n) {
    const absN = -n;
    const fib = fibonacci(absN);
    // (-1)^{n+1} means: if absN is even, return -fib; if absN is odd, return fib
    return absN % 2n === 0n ? -fib : fib;
  }

  if (n === 0n) {
    return 0n;
  }

  if (n === 1n || n === 2n) {
    return 1n;
  }

  // Use fast doubling method: O(log n) time
  // F(2k) = F(k) * (2*F(k+1) - F(k))
  // F(2k+1) = F(k)^2 + F(k+1)^2
  function fibPair(k: bigint): [bigint, bigint] {
    if (k === 0n) {
      return [0n, 1n]; // [F(0), F(1)]
    }

    const [a, b] = fibPair(k >> 1n); // [F(k/2), F(k/2 + 1)]
    const c = a * (2n * b - a); // F(2k/2) = F(k)
    const d = a * a + b * b; // F(2k/2 + 1) = F(k+1)

    if ((k & 1n) === 0n) {
      return [c, d]; // [F(k), F(k+1)]
    } else {
      return [d, c + d]; // [F(k+1), F(k+2)]
    }
  }

  return fibPair(n)[0];
}

/**
 * Return the n-th Lucas number.
 *
 * The Lucas numbers are defined by L_0 = 2, L_1 = 1, L_n = L_{n-1} + L_{n-2}.
 *
 * @param n - Non-negative integer
 * @returns L_n
 *
 * @see Reference: sage/rings/integer.pyx:lucas_number1
 */
export function lucas_number(n: bigint): bigint {
  // Handle negative n: L_{-n} = (-1)^n * L_n
  if (n < 0n) {
    const absN = -n;
    const luc = lucas_number(absN);
    // (-1)^n means: if absN is odd, negate; if even, keep same
    return absN % 2n === 1n ? -luc : luc;
  }

  if (n === 0n) {
    return 2n;
  }

  if (n === 1n) {
    return 1n;
  }

  // Use identity: L_n = F_{n-1} + F_{n+1}
  return fibonacci(n - 1n) + fibonacci(n + 1n);
}
