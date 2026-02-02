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

import { ArithmeticError, NotImplementedError, ValueError, ZeroDivisionError } from '../errors.js';
import { current_randstate } from '../misc/randstate.js';
import { Rational } from '../rings/rational.js';
import { type IntegerLike, type RationalLike, toBigInt, toRational, toSafeNumber } from '../types/coercion.js';

// Import PARI functions for factorization and primality testing
// Reference: sage/arith/misc.py uses PARI via cypari2
import {
  type Factorization as PariFactorization,
  Z_factor as pari_Z_factor,
  isPrime as pari_isPrime,
} from '@sagemath-ts/parigp-ts';

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
    throw new ZeroDivisionError('inverse_mod(a, 0) is not defined');
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

/**
 * Trial division factorization.
 *
 * Returns the smallest prime factor of n, or n if n is prime or 1.
 * If bound is specified, only considers primes up to that bound.
 *
 * Uses PARI's factorization via @sagemath-ts/parigp-ts.
 *
 * @param n - Integer to factor
 * @param bound - Maximum prime to try (default: all primes up to sqrt(n))
 * @returns The smallest prime factor ≤ bound, or n if no such factor exists
 *
 * @see Reference: sage/arith/misc.py:trial_division
 * @see Implementation: Uses factor() which delegates to PARI's Z_factor
 */
export function trial_division(n: IntegerLike, bound?: IntegerLike): bigint {
  let _n = toBigInt(n);
  const _bound = bound !== undefined ? toBigInt(bound) : undefined;
  if (_n < 0n) {
    _n = -_n;
  }

  if (_n <= 1n) {
    return _n;
  }

  // Get full factorization using PARI
  const factors = pari_Z_factor(_n);

  // Filter out the -1 sign factor and find the smallest prime
  for (const [p, _] of factors) {
    if (p > 0n) {
      // If bound is specified and the smallest prime exceeds it, return _n
      if (_bound !== undefined && p > _bound) {
        return _n;
      }
      return p;
    }
  }

  // n is 1 or -1
  return _n;
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
 * Test whether n is a perfect power (n = m^k for some k > 1).
 *
 * @param n - Integer to test
 * @param get_data - If true, return [true, base, exponent] or [false, n, 1]
 * @returns Whether n is a perfect power
 */
export function is_prime_power(n: IntegerLike, get_data?: false): boolean;
export function is_prime_power(n: IntegerLike, get_data: true): [bigint, bigint];
export function is_prime_power(
  n: IntegerLike,
  get_data: boolean = false
): boolean | [bigint, bigint] {
  const _n = toBigInt(n);
  // SageMath behavior for get_data=True:
  // - If prime power: returns (p, k) where n = p^k
  // - If not prime power: returns (n, 0)
  if (_n <= 1n) {
    if (get_data) {
      return [_n, 0n];
    }
    return false;
  }

  // Check if n is prime
  if (is_prime(_n)) {
    if (get_data) {
      return [_n, 1n];
    }
    return true;
  }

  // Find the smallest prime factor
  const p = trial_division(_n);
  if (p === _n) {
    // n is prime (shouldn't happen if is_prime is correct)
    if (get_data) {
      return [_n, 1n];
    }
    return true;
  }

  // Check if n = p^k
  let k = 0n;
  let m = _n;
  while (m % p === 0n) {
    m /= p;
    k++;
  }

  if (m === 1n) {
    // n = p^k
    if (get_data) {
      return [p, k];
    }
    return true;
  }

  // n has multiple distinct prime factors
  if (get_data) {
    return [_n, 0n];
  }
  return false;
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
  const _n = toBigInt(n);
  if (_n < 2n) {
    return 2n;
  }

  // Start at n+1, or n+2 if n+1 is even
  let candidate = _n + 1n;
  if ((candidate & 1n) === 0n) {
    candidate++;
  }

  while (!is_prime(candidate)) {
    candidate += 2n;
  }

  return candidate;
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
  const _n = toBigInt(n);
  if (_n <= 2n) {
    throw new ValueError('no prime less than 2');
  }

  if (_n === 3n) {
    return 2n;
  }

  let candidate = _n - 1n;
  if ((candidate & 1n) === 0n) {
    candidate--;
  }

  while (!is_prime(candidate)) {
    candidate -= 2n;
  }

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
 */
export function factor(n: IntegerLike): Factorization {
  const _n = toBigInt(n);
  if (_n === 0n) {
    throw new ValueError('factorization of 0 is not defined');
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
  const _n = toBigInt(n);
  // SageMath returns 0 for n <= 0
  if (_n <= 0n) {
    return 0n;
  }

  if (_n === 1n) {
    return 1n;
  }

  const factors = factor(_n);
  let result = 1n;

  for (const [p, e] of factors) {
    if (p === -1n) continue;
    // phi(p^e) = p^(e-1) * (p - 1)
    result *= (p - 1n) * p ** (e - 1n);
  }

  return result;
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
    // Handle (a/2)
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
 * Return the Legendre symbol (a/p).
 *
 * @param a - Integer
 * @param p - Odd prime
 * @returns The Legendre symbol value (-1, 0, or 1)
 *
 * @see Reference: sage/arith/misc.py:legendre_symbol
 */
export function legendre_symbol(a: IntegerLike, p: IntegerLike): bigint {
  const _p = toBigInt(p);
  if (_p === 2n) {
    throw new ValueError('p must be an odd prime');
  }
  return kronecker_symbol(a, _p);
}

/**
 * Return the Jacobi symbol (a/n).
 *
 * @param a - Integer
 * @param n - Positive odd integer
 * @returns The Jacobi symbol value (-1, 0, or 1)
 *
 * @see Reference: sage/arith/misc.py:jacobi_symbol
 */
export function jacobi_symbol(a: IntegerLike, n: IntegerLike): bigint {
  const _n = toBigInt(n);
  if (_n <= 0n || (_n & 1n) === 0n) {
    throw new ValueError('n must be a positive odd integer');
  }
  return kronecker_symbol(a, _n);
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
export function crt(a: IntegerLike, b: IntegerLike, m: IntegerLike, n: IntegerLike): bigint {
  const _a = toBigInt(a);
  const _b = toBigInt(b);
  const _m = toBigInt(m);
  const _n = toBigInt(n);
  const [g, s, _t] = xgcd(_m, _n);

  if ((_a - _b) % g !== 0n) {
    throw new ValueError('no solution exists');
  }

  const lcmMN = (_m / g) * _n;
  let result = (_a + _m * (((_b - _a) / g) * s)) % lcmMN;

  if (result < 0n) {
    result += lcmMN;
  }

  return result;
}

/**
 * Chinese Remainder Theorem for a list of residues and moduli.
 *
 * @param residues - List of residues
 * @param moduli - List of moduli (must be pairwise coprime)
 * @returns x such that x ≡ residues[i] (mod moduli[i]) for all i
 */
export function CRT_list(residues: IntegerLike[], moduli: IntegerLike[]): bigint {
  if (residues.length !== moduli.length) {
    throw new ValueError('residues and moduli must have the same length');
  }

  if (residues.length === 0) {
    return 0n;
  }

  let result = toBigInt(residues[0]!);
  let modulus = toBigInt(moduli[0]!);

  for (let i = 1; i < residues.length; i++) {
    result = crt(result, toBigInt(residues[i]!), modulus, toBigInt(moduli[i]!));
    modulus = lcm(modulus, toBigInt(moduli[i]!));
  }

  return result;
}

/**
 * Test whether n is a perfect square.
 *
 * @param n - Integer to test
 * @param root - If true, return [isSquare, root] instead of just boolean
 * @returns Whether n is a perfect square (and optionally its root)
 */
export function is_square(n: IntegerLike, root?: false): boolean;
export function is_square(n: IntegerLike, root: true): [boolean, bigint];
export function is_square(n: IntegerLike, root: boolean = false): boolean | [boolean, bigint] {
  const _n = toBigInt(n);
  if (_n < 0n) {
    if (root) {
      return [false, 0n];
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
    return [isSquareVal, isSquareVal ? s : 0n];
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
 * Return the number of positive divisors of n.
 *
 * @param n - Positive integer
 * @returns The number of positive divisors
 */
export function number_of_divisors(n: IntegerLike): bigint {
  const _n = toBigInt(n);
  if (_n <= 0n) {
    throw new ValueError('n must be positive');
  }

  if (_n === 1n) {
    return 1n;
  }

  const factors = factor(_n);
  let result = 1n;

  for (const [_p, e] of factors) {
    result *= e + 1n;
  }

  return result;
}

/**
 * Return the sum of the k-th powers of divisors of n.
 *
 * @param n - Positive integer
 * @param k - Power (default 1)
 * @returns Sum of d^k for all divisors d of n
 */
export function sigma(n: IntegerLike, k: IntegerLike = 1n): bigint {
  const _n = toBigInt(n);
  const _k = toBigInt(k);
  if (_n <= 0n) {
    throw new ValueError('n must be positive');
  }

  const divs = divisors(_n);
  let result = 0n;

  for (const d of divs) {
    result += d ** _k;
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
    throw new ValueError('prime_factors of 0 is not defined');
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
 * valuation(0n, 2n)    // Infinity is not representable, throws error
 * ```
 *
 * @see Reference: sage/arith/misc.py:valuation
 */
export function valuation(n: IntegerLike, p: IntegerLike): bigint {
  let _n = toBigInt(n);
  const _p = toBigInt(p);
  if (_p <= 1n) {
    throw new ValueError(
      'You can only compute the valuation with respect to an integer larger than 1.'
    );
  }

  if (_n === 0n) {
    // In SageMath this returns +Infinity, but we can't represent that
    // Following the convention of returning a very large value or throwing
    throw new ValueError('valuation of 0 is infinite');
  }

  // Work with absolute value
  if (_n < 0n) {
    _n = -_n;
  }

  let k = 0n;
  while (_n % _p === 0n) {
    _n /= _p;
    k++;
  }

  return k;
}

// ============================================================================
// STUB FUNCTIONS - Not yet implemented
// ============================================================================

/**
 * Return an irreducible polynomial of degree at most `degree` which
 * is approximately satisfied by the number `z`.
 *
 * You can specify the number of known bits or digits of `z` with
 * `known_bits=k` or `known_digits=k`. PARI is then told to
 * compute the result using `0.8k` of these bits/digits. Or, you can
 * specify the precision to use directly with `use_bits=k` or
 * `use_digits=k`. If none of these are specified, then the precision
 * is taken from the input value.
 *
 * ALGORITHM: Uses LLL for real/complex inputs.
 *
 * @param z - Real or complex number (as JavaScript number)
 * @param degree - Maximum degree of polynomial
 * @param options - Optional parameters for precision control
 * @returns Coefficients of an irreducible polynomial approximately satisfied by z
 *
 * @example
 * ```typescript
 * algebraic_dependency(1.888888888888888, 1n)  // coefficients for 9*x - 17
 * algebraic_dependency(Math.sqrt(2), 2n)       // coefficients for x^2 - 2
 * ```
 *
 * @see Reference: sage/arith/misc.py:algebraic_dependency
 */
export function algebraic_dependency(
  z: number,
  degree: bigint,
  options?: {
    known_bits?: bigint;
    use_bits?: bigint;
    known_digits?: bigint;
    use_digits?: bigint;
    height_bound?: bigint;
    proof?: boolean;
  }
): bigint[] {
  const opts = options ?? {};
  const log2_10 = Math.log(10) / Math.log(2);

  // Handle special cases
  if (!Number.isFinite(z)) {
    throw new ValueError('z must be a finite number');
  }

  const degreeNum = toSafeNumber(degree);
  if (degreeNum < 1) {
    throw new ValueError('degree must be at least 1');
  }

  // For integers, return x - value
  if (Number.isInteger(z)) {
    const intZ = BigInt(Math.round(z));
    if (opts.height_bound && (intZ < 0n ? -intZ : intZ) >= opts.height_bound) {
      throw new ValueError('no polynomial found within height bound');
    }
    return [-intZ, 1n];
  }

  // Determine precision to use
  let prec = 53 - 6; // Default IEEE 754 double precision minus safety margin
  if (opts.known_digits !== undefined) {
    prec = Math.floor(Number(opts.known_digits) * log2_10 * 0.8);
  } else if (opts.known_bits !== undefined) {
    prec = Math.floor(Number(opts.known_bits) * 0.8);
  } else if (opts.use_digits !== undefined) {
    prec = Math.floor(Number(opts.use_digits) * log2_10);
  } else if (opts.use_bits !== undefined) {
    prec = Number(opts.use_bits);
  }

  // Build the LLL matrix
  // M is (degree+1) x (degree+2) for real numbers
  const n = degreeNum + 1;
  const scale = 2 ** prec;
  const M: number[][] = [];

  for (let i = 0; i < n; i++) {
    const row: number[] = [];
    for (let j = 0; j < n; j++) {
      row.push(i === j ? 1 : 0);
    }
    // Last column: scaled power of z
    row.push(Math.round(scale * z ** i));
    M.push(row);
  }

  // Perform LLL reduction
  const lllReduced = _lllReduceSimple(M, 0.75);

  // Get coefficients from the first row (the shortest vector)
  let coeffs = lllReduced[0]!.slice(0, n);

  // If constant polynomial, try the second row
  let allButFirstZero = true;
  for (let i = 1; i < coeffs.length; i++) {
    if (Math.round(coeffs[i]!) !== 0) {
      allButFirstZero = false;
      break;
    }
  }
  if (allButFirstZero && lllReduced.length > 1) {
    coeffs = lllReduced[1]!.slice(0, n);
  }

  // Check height bound
  if (opts.height_bound) {
    const maxCoeff = Math.max(...coeffs.map((c) => Math.abs(Math.round(c))));
    if (maxCoeff > Number(opts.height_bound)) {
      throw new ValueError('no polynomial found within height bound');
    }
  }

  // Convert to bigint
  const result = coeffs.map((c) => BigInt(Math.round(c)));

  // Make leading coefficient positive
  let lastNonzero = -1;
  for (let i = result.length - 1; i >= 0; i--) {
    if (result[i] !== 0n) {
      lastNonzero = i;
      break;
    }
  }
  if (lastNonzero >= 0 && result[lastNonzero]! < 0n) {
    return result.map((c) => -c);
  }

  return result;
}

/**
 * Simple LLL reduction for algebraic_dependency.
 * Uses Lenstra-Lenstra-Lovasz algorithm with Gram-Schmidt orthogonalization.
 */
function _lllReduceSimple(M: number[][], delta: number): number[][] {
  const n = M.length;
  const m = M[0]!.length;
  const B = M.map((row) => [...row]);

  // Gram-Schmidt orthogonalization
  const GStar: number[][] = [];
  const mu: number[][] = [];
  const Bnorms: number[] = [];

  const computeGS = (): void => {
    GStar.length = 0;
    mu.length = 0;
    Bnorms.length = 0;

    for (let i = 0; i < n; i++) {
      GStar.push([...B[i]!]);
      mu.push(new Array(n).fill(0));

      for (let j = 0; j < i; j++) {
        let dot1 = 0;
        let dot2 = 0;
        for (let k = 0; k < m; k++) {
          dot1 += B[i]![k]! * GStar[j]![k]!;
          dot2 += GStar[j]![k]! * GStar[j]![k]!;
        }
        mu[i]![j] = dot2 !== 0 ? dot1 / dot2 : 0;

        for (let k = 0; k < m; k++) {
          GStar[i]![k] -= mu[i]![j]! * GStar[j]![k]!;
        }
      }

      let norm = 0;
      for (let k = 0; k < m; k++) {
        norm += GStar[i]![k]! * GStar[i]![k]!;
      }
      Bnorms.push(norm);
    }
  };

  const sizeReduce = (i: number, j: number): void => {
    if (Math.abs(mu[i]![j]!) > 0.5) {
      const q = Math.round(mu[i]![j]!);
      for (let k = 0; k < m; k++) {
        B[i]![k] -= q * B[j]![k]!;
      }
      for (let k = 0; k <= j; k++) {
        mu[i]![k] -= q * mu[j]![k]!;
      }
    }
  };

  computeGS();

  let k = 1;
  while (k < n) {
    // Size reduce B[k] against B[k-1]
    sizeReduce(k, k - 1);

    // Check Lovasz condition
    const lovaszCond = Bnorms[k]! >= (delta - mu[k]![k - 1]! * mu[k]![k - 1]!) * Bnorms[k - 1]!;

    if (lovaszCond) {
      // Size reduce B[k] against B[0], ..., B[k-2]
      for (let j = k - 2; j >= 0; j--) {
        sizeReduce(k, j);
      }
      k++;
    } else {
      // Swap B[k] and B[k-1]
      [B[k], B[k - 1]] = [B[k - 1]!, B[k]!];
      computeGS();
      k = Math.max(k - 1, 1);
    }
  }

  return B;
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
 */
export function factorial(n: IntegerLike, algorithm?: 'gmp' | 'pari'): bigint {
  const _n = toBigInt(n);
  if (_n < 0n) {
    throw new ValueError('factorial -- must be nonnegative');
  }

  if (_n === 0n || _n === 1n) {
    return 1n;
  }

  let result = 1n;
  for (let i = 2n; i <= _n; i++) {
    result *= i;
  }
  return result;
}

/**
 * Test whether n is a pseudo-prime.
 *
 * The result is NOT proven correct - this is a pseudo-primality test!
 * Uses Miller-Rabin with several bases.
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
 */
export function is_pseudoprime(n: IntegerLike): boolean {
  const _n = toBigInt(n);
  // Negative numbers and values <= 1 are not prime
  if (_n <= 1n) {
    return false;
  }

  // Check small primes
  if (_n === 2n || _n === 3n) {
    return true;
  }

  if ((_n & 1n) === 0n) {
    return false;
  }

  // Miller-Rabin with probabilistic witnesses
  const witnesses = [2n, 3n, 5n, 7n, 11n, 13n, 17n, 19n, 23n, 29n, 31n, 37n];

  for (const a of witnesses) {
    if (a >= _n) {
      continue;
    }
    if (!is_strong_probable_prime(_n, a)) {
      return false;
    }
  }

  return true;
}

/**
 * Test if n is a power of a pseudoprime.
 *
 * The result is NOT proven correct - this IS a pseudo-primality test!
 *
 * @param n - Integer to test
 * @param get_data - If true, return (p, k) such that n = p^k
 * @returns Whether n is a pseudoprime power
 *
 * @see Reference: sage/arith/misc.py:is_pseudoprime_power
 */
export function is_pseudoprime_power(n: bigint, get_data?: false): boolean;
export function is_pseudoprime_power(n: bigint, get_data: true): [bigint, bigint];
export function is_pseudoprime_power(
  n: bigint,
  get_data: boolean = false
): boolean | [bigint, bigint] {
  if (n <= 1n) {
    if (get_data) {
      return [n, 0n];
    }
    return false;
  }

  // Check if n is a pseudoprime
  if (is_pseudoprime(n)) {
    if (get_data) {
      return [n, 1n];
    }
    return true;
  }

  // Check if n is a perfect power of a pseudoprime
  // Find the smallest prime factor using trial division
  const p = trial_division(n);
  if (p === n) {
    // n is prime or pseudoprime (trial division didn't find a factor)
    if (is_pseudoprime(n)) {
      if (get_data) {
        return [n, 1n];
      }
      return true;
    }
    if (get_data) {
      return [n, 0n];
    }
    return false;
  }

  // Check if n = p^k
  let k = 0n;
  let m = n;
  while (m % p === 0n) {
    m /= p;
    k++;
  }

  if (m === 1n && is_pseudoprime(p)) {
    if (get_data) {
      return [p, k];
    }
    return true;
  }

  if (get_data) {
    return [n, 0n];
  }
  return false;
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
export function prime_powers(start: bigint, stop?: bigint): bigint[] {
  // Handle single argument case
  if (stop === undefined) {
    stop = start;
    start = 1n;
  }

  if (stop <= start || stop <= 2n) {
    return [];
  }

  if (start < 2n) {
    start = 2n;
  }

  const result: bigint[] = [];

  for (let n = start; n < stop; n++) {
    if (is_prime_power(n)) {
      result.push(n);
    }
  }

  return result;
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
export function primes_first_n(n: number): bigint[] {
  if (n < 0) {
    throw new ValueError('n must be nonnegative');
  }

  if (n === 0) {
    return [];
  }

  const primes: bigint[] = [];
  let candidate = 2n;

  while (primes.length < n) {
    if (is_prime(candidate)) {
      primes.push(candidate);
    }
    candidate++;
  }

  return primes;
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
export function eratosthenes(n: bigint): bigint[] {
  if (n < 2n) {
    return [];
  }

  return prime_range(2n, n + 1n);
}

/**
 * Return an iterator over all primes between start and stop-1, inclusive.
 *
 * @param start - Lower bound (default: 2)
 * @param stop - Upper bound (exclusive)
 * @param proof - Whether to use proven primality test (ignored, always uses BPSW)
 * @returns Iterator over primes
 *
 * @see Reference: sage/arith/misc.py:primes
 */
export function* primes(
  start?: bigint,
  stop?: bigint,
  proof?: boolean
): Generator<bigint, void, unknown> {
  // Handle argument patterns
  if (start === undefined) {
    start = 2n;
  }

  if (stop === undefined) {
    // Infinite sequence starting from start
    let candidate = start < 2n ? 2n : start;

    // Handle 2 specially
    if (candidate === 2n) {
      yield 2n;
      candidate = 3n;
    } else if ((candidate & 1n) === 0n) {
      candidate++;
    }

    while (true) {
      if (is_prime(candidate)) {
        yield candidate;
      }
      candidate += 2n;
    }
  } else {
    // Bounded sequence
    for (const p of prime_range(start, stop)) {
      yield p;
    }
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
export function next_prime_power(n: bigint): bigint {
  // For n < 2, the smallest prime power is 2
  if (n < 2n) {
    return 2n;
  }

  // Start searching from n+1
  let candidate = n + 1n;

  while (true) {
    if (is_prime_power(candidate)) {
      return candidate;
    }
    candidate++;
  }
}

/**
 * Return the next probable prime after n, as determined by pseudoprimality test.
 *
 * @param n - Integer
 * @returns Next probable prime
 *
 * @see Reference: sage/arith/misc.py:next_probable_prime
 */
export function next_probable_prime(n: bigint): bigint {
  if (n < 2n) {
    return 2n;
  }

  // Start at n+1, or n+2 if n+1 is even
  let candidate = n + 1n;
  if ((candidate & 1n) === 0n) {
    candidate++;
  }

  while (!is_pseudoprime(candidate)) {
    candidate += 2n;
  }

  return candidate;
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
export function previous_prime_power(n: bigint): bigint {
  if (n <= 2n) {
    throw new ValueError('no prime power less than 2');
  }

  // Start searching from n-1
  let candidate = n - 1n;

  while (candidate >= 2n) {
    if (is_prime_power(candidate)) {
      return candidate;
    }
    candidate--;
  }

  // Should never reach here since 2 is a prime power
  throw new ValueError('no prime power less than ' + n.toString());
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
export function xlcm(m: bigint, n: bigint): [bigint, bigint, bigint] {
  if (m <= 0n || n <= 0n) {
    throw new ValueError('xlcm requires positive integers');
  }

  const l = lcm(m, n);

  // We need to find m1, n1 such that:
  // - l = m1 * n1
  // - m1 | m
  // - n1 | n
  // - gcd(m1, n1) = 1

  // Start with m1 = m, n1 = l/m = n/gcd(m,n)
  // Then iteratively move common factors from m1 to n1 or vice versa

  let m1 = m;
  let n1 = l / m;

  // While gcd(m1, n1) > 1, move the common factor
  let g = gcd(m1, n1);
  while (g > 1n) {
    // Move g from m1 to n1, but we need to ensure n1 still divides n
    // Actually, we need a different approach: extract coprime parts

    // Remove common factors from m1
    while (gcd(m1, n1) > 1n) {
      const common = gcd(m1, n1);
      m1 /= common;
    }

    // Recalculate n1
    n1 = l / m1;

    g = gcd(m1, n1);
  }

  return [l, m1, n1];
}

/**
 * Return a CRT basis for the given moduli.
 *
 * Given moduli [m1, m2, ..., mn], returns basis elements [e1, e2, ..., en]
 * such that e_i ≡ 1 (mod m_i) and e_i ≡ 0 (mod m_j) for j ≠ i.
 *
 * This is useful for precomputation in repeated CRT applications.
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
 * ```
 *
 * @see Reference: sage/arith/misc.py:CRT_basis
 */
export function CRT_basis(
  moduli: bigint[],
  require_coprime_moduli: boolean = true
): bigint[] | [bigint[], boolean] {
  const n = moduli.length;
  if (n === 0) {
    return require_coprime_moduli ? [] : [[], true];
  }

  const cs: bigint[] = [];

  try {
    // Compute M = product of all moduli
    let M = 1n;
    for (const m of moduli) {
      M *= m;
    }

    for (const m of moduli) {
      const Mm = M / m;
      const [d, , v] = xgcd(m, Mm);
      if (d !== 1n) {
        throw new ValueError('moduli must be coprime');
      }
      // e_i = v * M_i mod M, where M_i = M / m_i
      let basis = (v * Mm) % M;
      if (basis < 0n) {
        basis += M;
      }
      cs.push(basis);
    }

    if (require_coprime_moduli) {
      return cs;
    }
    return [cs, true];
  } catch (e) {
    if (require_coprime_moduli) {
      throw e;
    }
    // For non-coprime moduli, compute a weaker basis using iterative approach
    // This follows SageMath's approach for non-coprime case
    const basis: bigint[] = [];
    let prevLcm = 1n;
    let prevBasis: bigint[] = [];

    for (let i = 0; i < n; i++) {
      const m = moduli[i]!;
      const newLcm = lcm(prevLcm, m);

      // Build new basis
      const newBasis: bigint[] = [];

      // Update previous basis elements
      for (let j = 0; j < i; j++) {
        const [g, s] = xgcd(prevLcm, m);
        // Scale previous basis element
        newBasis.push((((prevBasis[j]! * (1n - s * (prevLcm / g))) % newLcm) + newLcm) % newLcm);
      }

      // Add new basis element
      const [g, s] = xgcd(prevLcm, m);
      const newElem = (((s * prevLcm) % newLcm) + newLcm) % newLcm;
      newBasis.push(newElem);

      prevLcm = newLcm;
      prevBasis = newBasis;
    }

    return [prevBasis, false];
  }
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
export function CRT_vectors(X: bigint[][], moduli: bigint[]): bigint[] {
  if (X.length === 0 || X[0]!.length === 0) {
    return [];
  }

  const n = X.length;
  if (n !== moduli.length) {
    throw new ValueError('number of moduli must equal length of X');
  }

  // Get the CRT basis (allowing non-coprime moduli)
  const res = CRT_basis(moduli, false) as [bigint[], boolean];
  const a = res[0];
  const coprime = res[1];

  // Compute LCM of all moduli
  let modulus = 1n;
  for (const m of moduli) {
    modulus = lcm(modulus, m);
  }

  const vectorLen = X[0]!.length;

  // Compute candidate solution for each component
  const candidate: bigint[] = [];
  for (let j = 0; j < vectorLen; j++) {
    let sum = 0n;
    for (let i = 0; i < n; i++) {
      sum += a[i]! * X[i]![j]!;
    }
    let val = sum % modulus;
    if (val < 0n) {
      val += modulus;
    }
    candidate.push(val);
  }

  // If moduli are not coprime, verify the solution
  if (!coprime) {
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < vectorLen; j++) {
        if ((X[i]![j]! - candidate[j]!) % moduli[i]! !== 0n) {
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
export function nth_prime(n: bigint): bigint {
  if (n <= 0n) {
    throw new ValueError('nth prime meaningless for nonpositive n (=' + n.toString() + ')');
  }

  let count = 0n;
  let candidate = 2n;

  while (count < n) {
    if (is_prime(candidate)) {
      count++;
      if (count === n) {
        return candidate;
      }
    }
    candidate++;
  }

  // Should never reach here
  throw new ValueError('nth prime not found');
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
export function quadratic_residues(n: bigint): bigint[] {
  // Use absolute value of n
  n = n < 0n ? -n : n;

  if (n === 0n) {
    return [];
  }

  // Compute all squares modulo n
  // We only need to check a from 0 to n/2 since (n-a)^2 = a^2 mod n
  const residues = new Set<bigint>();
  const limit = n / 2n + 1n;

  for (let a = 0n; a <= limit; a++) {
    const square = (a * a) % n;
    residues.add(square);
  }

  // Convert to sorted array
  const result = Array.from(residues);
  result.sort((a, b) => (a < b ? -1 : a > b ? 1 : 0));
  return result;
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
export function continuant(v: bigint[], n?: bigint): bigint {
  const len = n !== undefined ? Number(n) : v.length;

  if (len <= 0 || v.length === 0) {
    return 1n;
  }

  if (len === 1) {
    return v[0]!;
  }

  // K_n = x_n * K_{n-1} + K_{n-2}
  let kPrev2 = 1n; // K_0
  let kPrev1 = v[0]!; // K_1

  for (let i = 1; i < len && i < v.length; i++) {
    const kCurr = v[i]! * kPrev1 + kPrev2;
    kPrev2 = kPrev1;
    kPrev1 = kCurr;
  }

  return kPrev1;
}

/**
 * Return 1 if ax^2 + by^2 p-adically represents a nonzero square,
 * otherwise returns -1. If either a or b is 0, returns 0.
 *
 * @param a - Integer
 * @param b - Integer
 * @param p - Prime or -1 (representing the archimedean place)
 * @param algorithm - Algorithm to use ('pari', 'direct', or 'all')
 * @returns 0, -1, or 1
 *
 * @see Reference: sage/arith/misc.py:hilbert_symbol
 */
export function hilbert_symbol(
  a: bigint,
  b: bigint,
  p: bigint,
  algorithm: 'pari' | 'direct' | 'all' = 'direct'
): bigint {
  // Reference: sage/arith/misc.py:hilbert_symbol
  // Reference: reference/pari/src/basemath/arith1.c:hilbert

  // Validate p
  if (p !== -1n && !is_prime(p)) {
    throw new ValueError('p must be prime or -1');
  }

  // Handle zero cases
  if (a === 0n || b === 0n) {
    return 0n;
  }

  // For the archimedean place (p = -1)
  if (p === -1n) {
    // hilbert(a, b, infinity) = -1 iff a < 0 and b < 0
    return a < 0n && b < 0n ? -1n : 1n;
  }

  const alg = algorithm || 'direct';

  if (alg === 'all') {
    // Both algorithms should agree
    const pariResult = hilbert_symbol(a, b, p, 'direct');
    return pariResult;
  }

  // 'direct' algorithm implementation (matching SageMath's direct algorithm)
  // First, remove p^2 factors
  const pSqr = p * p;
  while (a % pSqr === 0n) {
    a = a / pSqr;
  }
  while (b % pSqr === 0n) {
    b = b / pSqr;
  }

  // Check easy cases using Kronecker symbol
  if (p !== 2n) {
    // If any of a, b, or a+b is a quadratic residue mod p, return 1
    if (
      kronecker_symbol(a, p) === 1n ||
      kronecker_symbol(b, p) === 1n ||
      kronecker_symbol(a + b, p) === 1n
    ) {
      return 1n;
    }
  }

  // Check divisibility conditions
  const aDivP = a % p === 0n;
  const bDivP = b % p === 0n;

  if (aDivP) {
    if (bDivP) {
      // Both a and b divisible by p
      return hilbert_symbol(p, -(b / p), p, 'direct') * hilbert_symbol(a / p, b, p, 'direct');
    } else {
      // Only a divisible by p
      if (p === 2n && mod(b, 4n) === 3n) {
        if (kronecker_symbol(a + b, p) === -1n) {
          return -1n;
        }
      } else if (kronecker_symbol(b, p) === -1n) {
        return -1n;
      }
    }
  } else if (bDivP) {
    // Only b divisible by p
    if (p === 2n && mod(a, 4n) === 3n) {
      if (kronecker_symbol(a + b, p) === -1n) {
        return -1n;
      }
    } else if (kronecker_symbol(a, p) === -1n) {
      return -1n;
    }
  } else if (p === 2n && mod(a, 4n) === 3n && mod(b, 4n) === 3n) {
    // Neither divisible by p, but both congruent to 3 mod 4 at p=2
    return -1n;
  }

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
export function hilbert_conductor(a: bigint, b: bigint): bigint {
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
export function hilbert_conductor_inverse(d: bigint): [bigint, bigint] {
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
export function integer_ceil(x: number | bigint): bigint {
  if (typeof x === 'bigint') {
    return x;
  }

  if (!Number.isFinite(x)) {
    throw new ValueError('integer_ceil requires a finite input');
  }

  return BigInt(Math.ceil(x));
}

/**
 * Return the largest integer <= x.
 *
 * @param x - A number
 * @returns The floor of x
 *
 * @see Reference: sage/arith/misc.py:integer_floor
 */
export function integer_floor(x: number | bigint): bigint {
  if (typeof x === 'bigint') {
    return x;
  }

  if (!Number.isFinite(x)) {
    throw new ValueError('integer_floor requires a finite input');
  }

  return BigInt(Math.floor(x));
}

/**
 * Truncate to the integer closer to zero.
 *
 * @param x - A number
 * @returns The truncation of x
 *
 * @see Reference: sage/arith/misc.py:integer_trunc
 */
export function integer_trunc(x: number | bigint): bigint {
  if (typeof x === 'bigint') {
    return x;
  }

  if (!Number.isFinite(x)) {
    throw new ValueError('integer_trunc requires a finite input');
  }

  return BigInt(Math.trunc(x));
}

/**
 * Write the integer n as a sum of two integer squares if possible;
 * otherwise return null.
 *
 * A number n can be written as a sum of two squares if and only if
 * all prime powers p^e in its factorization with p ≡ 3 (mod 4) have e even.
 *
 * @param n - Integer
 * @returns Tuple [a, b] with a <= b such that n = a^2 + b^2, or null if impossible
 *
 * @example
 * ```typescript
 * two_squares(389n)  // [10n, 17n]
 * two_squares(21n)   // null (21 is not a sum of 2 squares)
 * two_squares(0n)    // [0n, 0n]
 * ```
 *
 * @see Reference: sage/arith/misc.py:two_squares
 */
export function two_squares(n: bigint): [bigint, bigint] | null {
  if (n < 0n) {
    return null;
  }

  if (n === 0n) {
    return [0n, 0n];
  }

  // Factor n
  const F = factor(n);

  // Check whether it is possible to write n as a sum of two squares:
  // All prime powers p^e must have p = 2 or p ≡ 1 (mod 4) or e even.
  for (const [p, e] of F) {
    if (p === -1n) continue;
    if (e % 2n === 1n && p % 4n === 3n) {
      return null;
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
 * otherwise return null.
 *
 * By Legendre's three-square theorem, a positive integer n can be expressed
 * as a sum of three squares if and only if n is NOT of the form 4^a(8b+7).
 *
 * @param n - Integer
 * @returns Tuple [a, b, c] with a <= b <= c such that n = a^2 + b^2 + c^2, or null if impossible
 *
 * @example
 * ```typescript
 * three_squares(389n)  // [1n, 8n, 18n]
 * three_squares(7n)    // null (7 = 4^0 * (8*0 + 7))
 * three_squares(0n)    // [0n, 0n, 0n]
 * ```
 *
 * @see Reference: sage/arith/misc.py:three_squares
 */
export function three_squares(n: bigint): [bigint, bigint, bigint] | null {
  if (n < 0n) {
    return null;
  }

  if (n === 0n) {
    return [0n, 0n, 0n];
  }

  // First, remove all factors of 4 from n
  // e = valuation(n, 2) // 2
  let e = 0n;
  let N = n;
  while (N % 4n === 0n) {
    N /= 4n;
    e++;
  }
  const m = 1n << e; // 2^e

  // Check if N is a perfect square
  const sqrtN = isqrt(N);
  if (sqrtN * sqrtN === N) {
    return [0n, 0n, sqrtN * m];
  }

  // Check if N ≡ 7 (mod 8) - then n is not a sum of 3 squares
  if (N % 8n === 7n) {
    return null;
  }

  // Find x such that N - x^2 can be written as sum of 2 squares
  // Strategy depends on N mod 8
  let x = sqrtN;

  if (N % 4n === 1n) {
    // Write N = x^2 + p with x even, p ≡ 1 (mod 4) prime
    if (x % 2n === 1n) {
      x -= 1n;
    }
    while (x >= 0n) {
      const p = N - x * x;
      if (is_prime(p) && p % 4n === 1n) {
        break;
      }
      x -= 2n;
    }
  } else if (N % 4n === 2n) {
    // Write N = x^2 + p with x odd, p ≡ 1 (mod 4) prime
    if (x % 2n === 0n) {
      x -= 1n;
    }
    while (x >= 0n) {
      const p = N - x * x;
      if (is_prime(p) && p % 4n === 1n) {
        break;
      }
      x -= 2n;
    }
  } else if (N % 8n === 3n) {
    // Write N = x^2 + 2p with x odd, p ≡ 1 (mod 4) prime
    if (x % 2n === 0n) {
      x -= 1n;
    }
    while (x >= 0n) {
      const remainder = N - x * x;
      if (remainder % 2n === 0n) {
        const p = remainder / 2n;
        if (is_prime(p) && p % 4n === 1n) {
          break;
        }
      }
      x -= 2n;
    }
  }

  // If we didn't find a good x, brute force
  if (x < 0n) {
    x = sqrtN;
  }

  // Try to write N - x^2 as sum of 2 squares
  while (x >= 0n) {
    const remainder = N - x * x;
    const twoSq = two_squares(remainder);
    if (twoSq !== null) {
      const [a, b] = twoSq;
      // Sort and return
      const result = [a, b, x].sort((p, q) => (p < q ? -1 : p > q ? 1 : 0));
      return [result[0]! * m, result[1]! * m, result[2]! * m];
    }
    x -= 1n;
  }

  // Should not reach here if the algorithm is correct
  return null;
}

/**
 * Write the integer n as a sum of four integer squares.
 *
 * By Lagrange's four square theorem, every non-negative integer can be
 * expressed as a sum of four integer squares.
 *
 * @param n - Non-negative integer
 * @returns Tuple [a, b, c, d] with a <= b <= c <= d such that n = a^2 + b^2 + c^2 + d^2,
 *          or null if n is negative
 *
 * @example
 * ```typescript
 * four_squares(3n)    // [0n, 1n, 1n, 1n]
 * four_squares(13n)   // [0n, 0n, 2n, 3n]
 * four_squares(130n)  // [0n, 0n, 3n, 11n]
 * ```
 *
 * @see Reference: sage/arith/misc.py:four_squares
 */
export function four_squares(n: bigint): [bigint, bigint, bigint, bigint] | null {
  if (n < 0n) {
    return null;
  }

  if (n === 0n) {
    return [0n, 0n, 0n, 0n];
  }

  // First, remove all factors of 4 from n
  let e = 0n;
  let N = n;
  while (N % 4n === 0n) {
    N /= 4n;
    e++;
  }
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

  const threeSq = three_squares(y);
  if (threeSq === null) {
    // This shouldn't happen by Lagrange's theorem, but handle it
    return null;
  }

  const [a, b, c] = threeSq;

  // Result is [a, b, c, x], all scaled by m
  const result = [a * m, b * m, c * m, x * m].sort((p, q) => (p < q ? -1 : p > q ? 1 : 0));
  return [result[0]!, result[1]!, result[2]!, result[3]!];
}

/**
 * Write the integer n as a sum of k integer squares if possible;
 * otherwise return null.
 *
 * For k >= 4, this always succeeds for non-negative n (Lagrange's theorem).
 * For k = 3, succeeds iff n is not of form 4^a(8b+7) (Legendre's theorem).
 * For k = 2, succeeds iff all primes p ≡ 3 (mod 4) in factorization have even exponent.
 * For k = 1, succeeds iff n is a perfect square.
 * For k = 0, succeeds iff n = 0.
 *
 * @param k - Non-negative integer
 * @param n - Integer
 * @returns Array [x_1, ..., x_k] of non-negative integers summing to n, or null if impossible
 *
 * @example
 * ```typescript
 * sum_of_k_squares(2, 9634n)  // [15n, 97n]
 * sum_of_k_squares(4, 9634n)  // [1n, 2n, 5n, 98n]
 * sum_of_k_squares(1, 9n)     // [3n]
 * sum_of_k_squares(1, 10n)    // null
 * sum_of_k_squares(0, 0n)     // []
 * ```
 *
 * @see Reference: sage/arith/misc.py:sum_of_k_squares
 */
export function sum_of_k_squares(k: number, n: bigint): bigint[] | null {
  if (k < 0) {
    return null;
  }

  // Handle base cases
  if (k === 0) {
    return n === 0n ? [] : null;
  }

  if (k === 1) {
    if (n < 0n) {
      return null;
    }
    const sqrtN = isqrt(n);
    if (sqrtN * sqrtN === n) {
      return [sqrtN];
    }
    return null;
  }

  if (k === 2) {
    const result = two_squares(n);
    return result ? [result[0], result[1]] : null;
  }

  if (k === 3) {
    const result = three_squares(n);
    return result ? [result[0], result[1], result[2]] : null;
  }

  if (k === 4) {
    const result = four_squares(n);
    return result ? [result[0], result[1], result[2], result[3]] : null;
  }

  // For k > 4, recursively subtract the largest square until k = 4
  if (n < 0n) {
    return null;
  }

  const extras: bigint[] = [];
  let remaining = n;
  let kRemaining = k;

  while (kRemaining > 4) {
    const x = isqrt(remaining);
    extras.push(x);
    remaining -= x * x;
    kRemaining--;
  }

  // Now write remaining as sum of 4 squares
  const fourSq = four_squares(remaining);
  if (fourSq === null) {
    return null;
  }

  // Combine and return
  const result = [...fourSq, ...extras.reverse()];
  return result;
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
export function subfactorial(n: bigint): bigint {
  if (n < 0n) {
    throw new ValueError('subfactorial only defined for non-negative integers');
  }

  if (n === 0n) {
    return 1n;
  }

  if (n === 1n) {
    return 0n;
  }

  // Use recurrence: !n = (n-1) * (!(n-1) + !(n-2))
  let prev2 = 1n; // !0
  let prev1 = 0n; // !1

  for (let k = 2n; k <= n; k++) {
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
export function is_power_of_two(n: bigint): boolean {
  // A number is a power of 2 if it has exactly one bit set (popcount == 1)
  // This is equivalent to: n > 0 && (n & (n - 1)) === 0
  if (n <= 0n) {
    return false;
  }
  return (n & (n - 1n)) === 0n;
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
 * @see Reference: sage/arith/misc.py:differences
 */
export function differences(lis: bigint[], n: bigint = 1n): bigint[] {
  if (n < 0n) {
    throw new ValueError('n must be non-negative');
  }

  if (n === 0n || lis.length === 0) {
    return [...lis];
  }

  let result = [...lis];

  for (let i = 0n; i < n && result.length > 1; i++) {
    const newResult: bigint[] = [];
    for (let j = 0; j < result.length - 1; j++) {
      newResult.push(result[j + 1]! - result[j]!);
    }
    result = newResult;
  }

  return result;
}

/**
 * Sort complex numbers in a "pretty" order for display.
 *
 * @param nums - List of complex numbers
 * @returns Sorted list
 *
 * @see Reference: sage/arith/misc.py:sort_complex_numbers_for_display
 */
export function sort_complex_numbers_for_display(nums: Array<{ re: number; im: number }>): Array<{
  re: number;
  im: number;
}> {
  // Reference: sage/arith/misc.py:sort_complex_numbers_for_display
  // Real numbers come before complex numbers, sorted by value
  // Complex numbers are sorted by real part (truncated for near-zero), then by imaginary part

  if (nums.length === 0) {
    return [];
  }

  const epsilon = 1e-10;

  // Key function for sorting
  const keyFunc = (a: { re: number; im: number }): [number, number, number] => {
    const ar = a.re;
    const ai = a.im;

    if (Math.abs(ai) < epsilon) {
      // Real number
      return [0, ar, 0];
    }

    // Complex number
    // Truncate real part if close to zero or if we have enough precision
    let arTruncated: number;
    if (Math.abs(ar) < epsilon) {
      arTruncated = 0;
    } else {
      // Truncate to ~9 significant figures
      arTruncated = Number.parseFloat(ar.toPrecision(9));
    }

    return [1, arTruncated, ai];
  };

  return [...nums].sort((a, b) => {
    const keyA = keyFunc(a);
    const keyB = keyFunc(b);

    // Compare element by element
    for (let i = 0; i < 3; i++) {
      if (keyA[i]! < keyB[i]!) return -1;
      if (keyA[i]! > keyB[i]!) return 1;
    }
    return 0;
  });
}

/**
 * Return the fundamental discriminant of Q(sqrt(D)).
 *
 * The fundamental discriminant of Q(sqrt(D)) is the discriminant of its ring
 * of integers. If D is squarefree, this is D if D ≡ 1 (mod 4), and 4D otherwise.
 *
 * @param D - Nonzero integer
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
export function fundamental_discriminant(D: bigint): bigint {
  if (D === 0n) {
    throw new ValueError('D must be nonzero');
  }

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
export function* squarefree_divisors(x: bigint): Generator<bigint, void, unknown> {
  if (x === 0n) {
    throw new ValueError('squarefree_divisors of 0 is not defined');
  }

  x = x < 0n ? -x : x;

  // Get prime factors
  const primeFactorsList = prime_factors(x);

  // Generate all subsets of prime factors (2^k combinations)
  const n = primeFactorsList.length;
  const total = 1 << n;

  for (let mask = 0; mask < total; mask++) {
    let divisor = 1n;
    for (let i = 0; i < n; i++) {
      if (mask & (1 << i)) {
        divisor *= primeFactorsList[i]!;
      }
    }
    yield divisor;
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
 */
export function dedekind_sum(
  p: bigint,
  q: bigint,
  algorithm?: 'default' | 'flint' | 'pari'
): { numerator: bigint; denominator: bigint } {
  // Reference: sage/arith/misc.py:dedekind_sum
  // Reference: reference/pari/src/basemath/elltrans.c:sumdedekind
  // s(h, k) = sum(n=1 to k-1, (n/k) * (frac(h*n/k) - 1/2))
  // Using Knuth's algorithm for coprime h, k

  // First reduce to coprime case
  const d = gcd(p < 0n ? -p : p, q < 0n ? -q : q);
  let h = p / d;
  let k = q / d;

  // Handle sign: s(-h, k) = -s(h, k), s(h, -k) = s(h, k)
  if (k < 0n) {
    k = -k;
  }
  const negateResult = h < 0n;
  if (h < 0n) {
    h = -h;
  }

  // Handle trivial cases
  if (k === 0n) {
    return { numerator: 0n, denominator: 1n };
  }
  if (k === 1n) {
    return { numerator: 0n, denominator: 1n };
  }

  h = mod(h, k);

  // Knuth's algorithm for computing s(h, k)
  // Reference: PARI's u_sumdedekind_coprime
  // Returns s(h, k) = (s2 + k * s1) / (12 * k)
  let sign = 1n;
  let s1 = 0n;
  let s2 = h;
  let pVal = 1n;
  let pp = 0n;

  while (h !== 0n) {
    const nexth = k % h;
    const a = k / h; // a >= 1, a >= 2 if h == 1

    // When h == 1, this is the last iteration
    if (h === 1n) {
      s2 = s2 + pVal * sign;
    }
    s1 = s1 + a * sign;
    sign = -sign;
    k = h;
    h = nexth;
    const r = a * pVal + pp;
    pp = pVal;
    pVal = r;
  }

  // At this point pVal equals the original k
  if (sign < 0n) {
    s1 = s1 - 3n;
  }

  // s(h, k) = (s2 + pVal * s1) / (12 * pVal)
  let num = s2 + pVal * s1;
  let den = 12n * pVal;

  // Simplify
  const g = gcd(num < 0n ? -num : num, den);
  num = num / g;
  den = den / g;

  if (negateResult) {
    num = -num;
  }

  return { numerator: num, denominator: den };
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
export function dedekind_psi(N: bigint): bigint {
  if (N <= 0n) {
    throw new ValueError('Dedekind psi function requires a positive integer');
  }

  if (N === 1n) {
    return 1n;
  }

  // psi(n) = n * product_{p|n}(1 + 1/p) = n * product_{p|n}((p+1)/p)
  const primeFactorsList = prime_factors(N);

  let numerator = N;
  let denominator = 1n;

  for (const p of primeFactorsList) {
    numerator *= p + 1n;
    denominator *= p;
  }

  return numerator / denominator;
}

/**
 * Return the largest divisor of x that is smooth over the factor base.
 *
 * The smooth part is the largest divisor of x whose prime factors all appear in base.
 *
 * @param x - Element of a Euclidean domain
 * @param base - Factor base (sequence of primes)
 * @returns Factorization of smooth part
 *
 * @example
 * ```typescript
 * smooth_part(240n, [2n, 3n])  // [[2n, 4n], [3n, 1n]] (240 = 2^4 * 3 * 5, smooth part is 2^4 * 3)
 * ```
 *
 * @see Reference: sage/arith/misc.py:smooth_part
 */
export function smooth_part(x: bigint, base: bigint[]): Factorization {
  if (x === 0n) {
    return [];
  }

  const sign = x < 0n ? -1n : 1n;
  let remaining = x < 0n ? -x : x;
  const result: Factorization = [];

  if (sign === -1n) {
    result.push([-1n, 1n]);
  }

  // Create a set from base for efficient lookup
  const baseSet = new Set(base.map((b) => (b < 0n ? -b : b)));

  // Extract factors that are in the base
  for (const b of base) {
    const p = b < 0n ? -b : b;
    if (p <= 1n) continue;

    let e = 0n;
    while (remaining % p === 0n) {
      remaining /= p;
      e++;
    }

    if (e > 0n) {
      result.push([p, e]);
    }
  }

  // Sort by prime
  result.sort((a, b) => {
    if (a[0] < b[0]) return -1;
    if (a[0] > b[0]) return 1;
    return 0;
  });

  return result;
}

/**
 * Return the largest divisor of x coprime to all elements of base.
 *
 * This is x with all prime factors that appear in any element of base removed.
 *
 * @param x - Element of a Euclidean domain
 * @param base - Factor base
 * @returns The coprime part
 *
 * @example
 * ```typescript
 * coprime_part(240n, [2n, 3n])  // 5n
 * coprime_part(240n, [6n])      // 5n (same, since 6 = 2*3)
 * ```
 *
 * @see Reference: sage/arith/misc.py:coprime_part
 */
export function coprime_part(x: bigint, base: bigint[]): bigint {
  if (x === 0n) {
    return 0n;
  }

  const sign = x < 0n ? -1n : 1n;
  let result = x < 0n ? -x : x;

  for (const b of base) {
    result = prime_to_m_part(result, b);
  }

  return sign * result;
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
export function carmichael_lambda(n: bigint): bigint {
  // The Carmichael function lambda(n) is the smallest positive integer k
  // such that a^k ≡ 1 (mod n) for all a coprime to n.
  //
  // Algorithm:
  // - If n = 2 or 4: lambda(n) = phi(n)
  // - If n = 2^k with k >= 3: lambda(2^k) = 2^(k-2)
  // - If n = p^k for odd prime p: lambda(p^k) = phi(p^k) = p^(k-1) * (p-1)
  // - For composite n = p1^k1 * p2^k2 * ...: lambda(n) = lcm(lambda(p1^k1), lambda(p2^k2), ...)

  if (n < 1n) {
    throw new ValueError('Input n must be a positive integer.');
  }

  if (n === 1n) {
    return 1n;
  }

  const factors = factor(n);
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
export function odd_part(n: bigint): bigint {
  // Return n / 2^v where v = valuation(n, 2)
  // This is n with all factors of 2 removed

  if (n === 0n) {
    return 0n;
  }

  // Handle negative numbers - return the odd part of |n| with the same sign
  const sign = n < 0n ? -1n : 1n;
  let result = n < 0n ? -n : n;

  // Divide out all factors of 2
  while ((result & 1n) === 0n) {
    result >>= 1n;
  }

  return sign * result;
}

/**
 * Return the prime-to-m part of n.
 *
 * This is the largest divisor of n that is coprime to m.
 * It equals n / gcd(n, m^infinity), i.e., n with all prime factors
 * that also divide m removed.
 *
 * @param n - Integer (nonzero)
 * @param m - Integer
 * @returns The prime-to-m part of n
 *
 * @example
 * ```typescript
 * prime_to_m_part(240n, 2n)  // 15n (240 = 2^4 * 15)
 * prime_to_m_part(240n, 6n)  // 5n (240 = 2^4 * 3 * 5, remove 2 and 3)
 * ```
 *
 * @see Reference: sage/arith/misc.py:prime_to_m_part
 */
export function prime_to_m_part(n: bigint, m: bigint): bigint {
  if (n === 0n) {
    return 0n;
  }

  if (m === 0n || m === 1n || m === -1n) {
    return n;
  }

  // Work with absolute values
  const sign = n < 0n ? -1n : 1n;
  let result = n < 0n ? -n : n;
  m = m < 0n ? -m : m;

  // Repeatedly divide out gcd(result, m) until result is coprime to m
  let g = gcd(result, m);
  while (g > 1n) {
    // Divide out all powers of g from result
    while (result % g === 0n) {
      result /= g;
    }
    g = gcd(result, m);
  }

  return sign * result;
}

/**
 * Return the fastest gcd function for integers of size no larger than order.
 *
 * In this implementation, we simply return the standard gcd function since
 * JavaScript's BigInt handles all sizes uniformly.
 *
 * @param order - Maximum size of integers (ignored)
 * @returns A gcd function
 *
 * @see Reference: sage/arith/misc.py:get_gcd
 */
export function get_gcd(order: bigint): (a: bigint, b: bigint) => bigint {
  return gcd;
}

/**
 * Return the fastest inverse_mod function for integers of size no larger than order.
 *
 * In this implementation, we simply return the standard inverse_mod function.
 *
 * @param order - Maximum size of integers (ignored)
 * @returns An inverse_mod function
 *
 * @see Reference: sage/arith/misc.py:get_inverse_mod
 */
export function get_inverse_mod(order: bigint): (a: bigint, m: bigint) => bigint {
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
 */
export function mqrr_rational_reconstruction(
  u: bigint,
  m: bigint,
  T: bigint
): [bigint, bigint] | null {
  // Reference: sage/arith/misc.py:mqrr_rational_reconstruction
  // Maximal Quotient Rational Reconstruction
  // Reference: Monagan, "Maximal Quotient Rational Reconstruction"

  // Handle u = 0 case
  if (u === 0n) {
    if (m > T) {
      return [0n, 1n];
    } else {
      return null;
    }
  }

  let n = 0n;
  let d = 0n;
  let t0 = 0n;
  let r0 = m;
  let t1 = 1n;
  let r1 = u;

  while (r1 !== 0n && r0 > T) {
    const q = r0 / r1; // Integer division (floor)
    if (q > T) {
      n = r1;
      d = t1;
      T = q;
    }
    const tempR = r1;
    r1 = r0 - q * r1;
    r0 = tempR;

    const tempT = t1;
    t1 = t0 - q * t1;
    t0 = tempT;
  }

  if (d !== 0n && gcd(n < 0n ? -n : n, d < 0n ? -d : d) === 1n) {
    // Ensure d is positive
    if (d < 0n) {
      return [-n, -d];
    }
    return [n, d];
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
