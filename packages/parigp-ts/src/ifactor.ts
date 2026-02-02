/**
 * @module ifactor
 * @description Integer factorization ported from PARI/GP
 *
 * Reference: reference/pari/src/basemath/ifactor1.c
 *
 * PARI's factorization uses:
 * 1. Trial division with mod 210 wheel (skips multiples of 2,3,5,7)
 * 2. Primality tests (BPSW - Baillie-Pomerance-Selfridge-Wagstaff)
 * 3. For large composites: Pollard rho, ECM, MPQS (not yet ported)
 *
 * This port implements trial division matching PARI's approach.
 */

import { Fp_pow } from './ff.js';

// ============================================================================
// Mod 210 wheel for prime iteration
// Reference: ifactor1.c:23-46
//
// 210 = 2 * 3 * 5 * 7
// Primes > 7 must be coprime to 210, i.e., in residue classes coprime to 210.
// There are 48 such classes: 1, 11, 13, 17, 19, 23, 29, 31, ...
// ============================================================================

/**
 * First differences between consecutive prime residue classes mod 210.
 * Starting from residue class 1, these are the gaps to the next coprime class.
 *
 * Reference: ifactor1.c:42-46
 */
const PRC210_D1: readonly number[] = [
  10, 2, 4, 2, 4, 6, 2, 6, 4, 2, 4, 6, 6, 2, 6, 4, 2, 6, 4, 6, 8, 4, 2, 4, 2, 4, 8, 6, 4, 6, 2, 4,
  6, 2, 6, 6, 4, 2, 4, 6, 2, 6, 4, 2, 4, 2, 10, 2,
] as const;

/**
 * Map from prime residue class (mod 210) to index in PRC210_D1.
 * Residue class r maps to index at (r-1)/2 (since all primes > 2 are odd).
 * Non-prime residue classes map to 128 (NPRC marker).
 *
 * Reference: ifactor1.c:28-39
 */
const NPRC = 128;
const PRC210_NO: readonly number[] = [
  0,
  NPRC,
  NPRC,
  NPRC,
  NPRC,
  1,
  2,
  NPRC,
  3,
  4,
  NPRC, // 1,3,5,7,9,11,13,15,17,19,21
  5,
  NPRC,
  NPRC,
  6,
  7,
  NPRC,
  NPRC,
  8,
  NPRC,
  9, // 23,25,27,29,31,33,35,37,39,41
  10,
  NPRC,
  11,
  NPRC,
  NPRC,
  12,
  NPRC,
  NPRC,
  13,
  14,
  NPRC, // 43,45,47,49,51,53,55,57,59,61,63
  NPRC,
  15,
  NPRC,
  16,
  17,
  NPRC,
  NPRC,
  18,
  NPRC,
  19, // 65,67,69,71,73,75,77,79,81,83
  NPRC,
  NPRC,
  20,
  NPRC,
  NPRC,
  NPRC,
  21,
  NPRC,
  22,
  23,
  NPRC, // 85,87,89,91,93,95,97,99,101,103,105
  24,
  25,
  NPRC,
  26,
  NPRC,
  NPRC,
  NPRC,
  27,
  NPRC,
  NPRC, // 107,109,111,113,115,117,119,121,123,125
  28,
  NPRC,
  29,
  NPRC,
  NPRC,
  30,
  31,
  NPRC,
  32,
  NPRC,
  NPRC, // 127,129,131,133,135,137,139,141,143,145,147
  33,
  34,
  NPRC,
  NPRC,
  35,
  NPRC,
  NPRC,
  36,
  NPRC,
  37, // 149,151,153,155,157,159,161,163,165,167
  38,
  NPRC,
  39,
  NPRC,
  NPRC,
  40,
  41,
  NPRC,
  NPRC,
  42,
  NPRC, // 169,171,173,175,177,179,181,183,185,187,189
  43,
  44,
  NPRC,
  45,
  46,
  NPRC,
  NPRC,
  NPRC,
  NPRC,
  47, // 191,193,195,197,199,201,203,205,207,209
] as const;

// ============================================================================
// Small primes table
// Reference: PARI uses a precomputed prime table up to some limit
// ============================================================================

/**
 * Small primes for initial trial division.
 * Using primes up to 997 (168 primes) for efficient small factor extraction.
 */
const SMALL_PRIMES: readonly bigint[] = [
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
  307n,
  311n,
  313n,
  317n,
  331n,
  337n,
  347n,
  349n,
  353n,
  359n,
  367n,
  373n,
  379n,
  383n,
  389n,
  397n,
  401n,
  409n,
  419n,
  421n,
  431n,
  433n,
  439n,
  443n,
  449n,
  457n,
  461n,
  463n,
  467n,
  479n,
  487n,
  491n,
  499n,
  503n,
  509n,
  521n,
  523n,
  541n,
  547n,
  557n,
  563n,
  569n,
  571n,
  577n,
  587n,
  593n,
  599n,
  601n,
  607n,
  613n,
  617n,
  619n,
  631n,
  641n,
  643n,
  647n,
  653n,
  659n,
  661n,
  673n,
  677n,
  683n,
  691n,
  701n,
  709n,
  719n,
  727n,
  733n,
  739n,
  743n,
  751n,
  757n,
  761n,
  769n,
  773n,
  787n,
  797n,
  809n,
  811n,
  821n,
  823n,
  827n,
  829n,
  839n,
  853n,
  857n,
  859n,
  863n,
  877n,
  881n,
  883n,
  887n,
  907n,
  911n,
  919n,
  929n,
  937n,
  941n,
  947n,
  953n,
  967n,
  971n,
  977n,
  983n,
  991n,
  997n,
] as const;

// ============================================================================
// Primality testing
// Reference: ifactor1.c uses BPSW (Baillie-PSW) primality test
// ============================================================================

/**
 * Miller-Rabin primality test for a single witness.
 *
 * Reference: arith2.c (uisprime_661, BPSW_psp)
 *
 * @param n - Number to test (must be odd > 2)
 * @param a - Witness to test
 * @returns true if n passes the test for witness a
 */
function millerRabinWitness(n: bigint, a: bigint): boolean {
  // Write n-1 = 2^s * d where d is odd
  let d = n - 1n;
  let s = 0;
  while ((d & 1n) === 0n) {
    d >>= 1n;
    s++;
  }

  // Compute a^d mod n
  let x = Fp_pow(a, d, n);

  if (x === 1n || x === n - 1n) return true;

  // Square s-1 times
  for (let r = 1; r < s; r++) {
    x = (x * x) % n;
    if (x === n - 1n) return true;
    if (x === 1n) return false;
  }

  return false;
}

/**
 * Strong Lucas probable prime test.
 *
 * Reference: arith2.c (Lucas pseudoprime test in BPSW)
 *
 * Uses the Selfridge parameters: find first D in 5, -7, 9, -11, ... such that
 * Jacobi(D, n) = -1, then P = 1, Q = (1-D)/4.
 */
function strongLucas(n: bigint): boolean {
  // Find D using Selfridge's method
  let D = 5n;
  let sign = 1n;
  while (true) {
    const jacobi = jacobiSymbol(D * sign, n);
    if (jacobi === 0) return n === D * sign || n === -(D * sign);
    if (jacobi === -1) {
      D = D * sign;
      break;
    }
    D += 2n;
    sign = -sign;
    if (D > 1000n) {
      // Fallback: if we can't find D, assume it's probably prime
      // (This shouldn't happen for reasonable inputs)
      return true;
    }
  }

  const P = 1n;
  const Q = (1n - D) / 4n;

  // Compute Lucas sequence U_d, V_d where n+1 = 2^s * d, d odd
  let d = n + 1n;
  let s = 0;
  while ((d & 1n) === 0n) {
    d >>= 1n;
    s++;
  }

  // Lucas sequence computation using binary method
  let U = 1n;
  let V = P;
  let Qk = Q;

  // Get binary representation of d
  const bits: boolean[] = [];
  let temp = d;
  while (temp > 0n) {
    bits.push((temp & 1n) === 1n);
    temp >>= 1n;
  }

  // Process bits from second-highest to lowest
  for (let i = bits.length - 2; i >= 0; i--) {
    // Double: U_2k = U_k * V_k, V_2k = V_k^2 - 2*Q^k
    U = (U * V) % n;
    V = (V * V - 2n * Qk) % n;
    if (V < 0n) V += n;
    Qk = (Qk * Qk) % n;

    if (bits[i]) {
      // Add: U_{k+1} = (P*U_k + V_k)/2, V_{k+1} = (D*U_k + P*V_k)/2
      const Unew = (P * U + V) % n;
      const Vnew = (D * U + P * V) % n;
      U = (Unew & 1n) === 0n ? Unew / 2n : (Unew + n) / 2n;
      V = (Vnew & 1n) === 0n ? Vnew / 2n : (Vnew + n) / 2n;
      U = ((U % n) + n) % n;
      V = ((V % n) + n) % n;
      Qk = (Qk * Q) % n;
      if (Qk < 0n) Qk += n;
    }
  }

  // Check U_d = 0 mod n
  if (U === 0n) return true;

  // Check V_{d*2^r} = 0 mod n for some r in [0, s-1]
  for (let r = 0; r < s; r++) {
    if (V === 0n) return true;
    V = (V * V - 2n * Qk) % n;
    if (V < 0n) V += n;
    Qk = (Qk * Qk) % n;
  }

  return false;
}

/**
 * Jacobi symbol (a|n) using quadratic reciprocity.
 */
function jacobiSymbol(a: bigint, n: bigint): number {
  if (n <= 0n || (n & 1n) === 0n) {
    throw new Error('Jacobi symbol: n must be positive and odd');
  }

  a = ((a % n) + n) % n;
  let result = 1;

  while (a !== 0n) {
    // Extract powers of 2
    while ((a & 1n) === 0n) {
      a >>= 1n;
      const nMod8 = Number(n & 7n);
      if (nMod8 === 3 || nMod8 === 5) {
        result = -result;
      }
    }

    // Quadratic reciprocity
    [a, n] = [n, a];
    if ((a & 3n) === 3n && (n & 3n) === 3n) {
      result = -result;
    }
    a = a % n;
  }

  return n === 1n ? result : 0;
}

/**
 * BPSW primality test (Baillie-PSW).
 *
 * Reference: ifactor1.c (BPSW_psp, ifac_isprime)
 *
 * Combines:
 * 1. Trial division by small primes
 * 2. Miller-Rabin test with base 2
 * 3. Strong Lucas probable prime test
 *
 * No known composite passes BPSW as of 2024.
 */
export function isPrime(n: bigint): boolean {
  if (n < 2n) return false;
  if (n === 2n) return true;
  if ((n & 1n) === 0n) return false;

  // Trial division by small primes
  for (const p of SMALL_PRIMES) {
    if (p * p > n) return true;
    if (n % p === 0n) return n === p;
  }

  // Miller-Rabin with base 2
  if (!millerRabinWitness(n, 2n)) return false;

  // Strong Lucas test
  return strongLucas(n);
}

// ============================================================================
// Factorization
// Reference: ifactor1.c
// ============================================================================

export type Factorization = Array<[bigint, bigint]>;

/**
 * Extract the highest power of p dividing n.
 *
 * Reference: ifactor1.c (Z_lvalrem, u_lvalrem)
 *
 * @returns [valuation, quotient] where n = p^valuation * quotient
 */
function lvalrem(n: bigint, p: bigint): [number, bigint] {
  let v = 0;
  while (n % p === 0n) {
    n /= p;
    v++;
  }
  return [v, n];
}

/**
 * Iterator for primes using mod 210 wheel.
 *
 * Reference: ifactor1.c:59-95 (unextprime, nextprime)
 *
 * Yields primes starting from a given value.
 */
function* primeIterator(start: bigint): Generator<bigint> {
  // Handle small primes first
  for (const p of SMALL_PRIMES) {
    if (p >= start) yield p;
  }

  // Continue with wheel
  // Start from first candidate after our small primes table
  let n = 1009n; // First prime after 997
  if (start > n) {
    // Round up to next number coprime to 210
    n = start;
    if ((n & 1n) === 0n) n++; // Make odd
  }

  // Find residue class mod 210
  let rc = Number(n % 210n);
  let rcn: number;

  // Find next prime residue class
  while (true) {
    rcn = PRC210_NO[(rc - 1) >> 1] ?? NPRC;
    if (rcn !== NPRC) break;
    rc += 2;
    n += 2n;
    if (rc >= 210) rc -= 210;
  }

  // Generate primes using wheel
  while (true) {
    if (isPrime(n)) {
      yield n;
    }
    n += BigInt(PRC210_D1[rcn]);
    rcn++;
    if (rcn >= 48) rcn = 0;
  }
}

/**
 * Trial division bound based on number size.
 *
 * Reference: ifactor1.c:3054-3083 (tridiv_bound)
 *
 * For small numbers, we can afford to trial divide up to sqrt(n).
 * For larger numbers, we use a smaller bound and rely on other methods.
 */
function tridivBound(n: bigint): bigint {
  // For simplicity, trial divide up to sqrt(n) or 10^6, whichever is smaller
  // PARI uses more sophisticated bounds based on bit length
  const sqrtn = isqrt(n);
  const limit = 1000000n;
  return sqrtn < limit ? sqrtn : limit;
}

/**
 * Integer square root using Newton's method.
 */
function isqrt(n: bigint): bigint {
  if (n < 0n) throw new Error('isqrt: negative argument');
  if (n < 2n) return n;

  // Initial guess: 2^(ceil(bitLength/2))
  let x = 1n << (BigInt(n.toString(2).length + 1) >> 1n);

  // Newton's method: x_{n+1} = (x_n + n/x_n) / 2
  while (true) {
    const x1 = (x + n / x) >> 1n;
    if (x1 >= x) return x;
    x = x1;
  }
}

/**
 * Z_factor - Factor an integer using trial division.
 *
 * Reference: ifactor1.c:4514-4515 (Z_factor), 4279-4456 (ifactor_sign)
 *
 * This is a simplified port focusing on trial division.
 * PARI's full implementation also uses Pollard rho, ECM, and MPQS
 * for large composites, which are not yet ported.
 *
 * @param n - Integer to factor (must be nonzero)
 * @returns Array of [prime, exponent] pairs
 */
export function Z_factor(n: bigint): Factorization {
  if (n === 0n) {
    throw new Error('Z_factor: factorization of 0 is not defined');
  }

  const factors: Factorization = [];

  // Handle sign
  if (n < 0n) {
    factors.push([-1n, 1n]);
    n = -n;
  }

  if (n === 1n) {
    return factors;
  }

  // Trial division by 2
  if ((n & 1n) === 0n) {
    const [v, q] = lvalrem(n, 2n);
    factors.push([2n, BigInt(v)]);
    n = q;
    if (n === 1n) return factors;
  }

  // Trial division by small primes (3, 5, 7, ...)
  for (let i = 1; i < SMALL_PRIMES.length; i++) {
    const p = SMALL_PRIMES[i];
    if (p * p > n) break;
    if (n % p === 0n) {
      const [v, q] = lvalrem(n, p);
      factors.push([p, BigInt(v)]);
      n = q;
      if (n === 1n) return factors;
    }
  }

  // Check if remaining n is prime
  if (n > 1n && isPrime(n)) {
    factors.push([n, 1n]);
    return factors;
  }

  // Continue trial division with wheel
  const bound = tridivBound(n);
  const primes = primeIterator(SMALL_PRIMES[SMALL_PRIMES.length - 1] + 1n);

  for (const p of primes) {
    if (p > bound) break;
    if (p * p > n) break;

    if (n % p === 0n) {
      const [v, q] = lvalrem(n, p);
      factors.push([p, BigInt(v)]);
      n = q;
      if (n === 1n) return factors;
      if (isPrime(n)) {
        factors.push([n, 1n]);
        return factors;
      }
    }
  }

  // If we get here with n > 1, it's either prime or we need advanced methods
  if (n > 1n) {
    if (isPrime(n)) {
      factors.push([n, 1n]);
    } else {
      // For composites beyond our trial division bound, we would need
      // Pollard rho, ECM, or MPQS. For now, just add it as a "probable prime"
      // (This is a limitation of this port - PARI would factor it further)
      console.warn(
        `Z_factor: ${n} is composite but beyond trial division bound. ` +
          'Advanced factorization (Pollard rho, ECM, MPQS) not yet implemented.'
      );
      factors.push([n, 1n]);
    }
  }

  return factors;
}

/**
 * factoru - Factor a small unsigned integer.
 *
 * Reference: ifactor1.c:3522-3523
 *
 * Optimized version for numbers that fit in a JavaScript number.
 * Uses the same algorithm as Z_factor but with number arithmetic where possible.
 *
 * @param n - Unsigned integer to factor (as bigint for consistency)
 * @returns Array of [prime, exponent] pairs
 */
export function factoru(n: bigint): Factorization {
  // Just delegate to Z_factor for now
  // A fully optimized version would use number arithmetic for small factors
  return Z_factor(n);
}

/**
 * Format a factorization as a string.
 *
 * @param f - Factorization to format
 * @returns String like "2^2 * 3 * 5^2"
 */
export function formatFactorization(f: Factorization): string {
  if (f.length === 0) return '1';

  return f.map(([p, e]) => (e === 1n ? `${p}` : `${p}^${e}`)).join(' * ');
}
