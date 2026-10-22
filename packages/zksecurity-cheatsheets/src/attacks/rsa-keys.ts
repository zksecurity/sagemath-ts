/**
 * @module rsa-keys
 * @description Checks RSA public keys for known key-generation flaws, alone or as a batch.
 *
 * A check that factors n or recovers d reports `broken`. A check that recognizes a known
 * flaw without exploiting it reports `vulnerable`; its false-positive rate is negligible.
 * Policy findings report `warning`.
 *
 * References:
 * - Fermat factorization; B. de Weger, "Cryptanalysis of RSA with small prime difference" (2002)
 * - R. S. Lehman, "Factoring large integers" (1974): Fermat with a multiplier
 * - M. Wiener, "Cryptanalysis of short RSA secret exponents" (1990)
 * - J. M. Pollard, "Theorems on factorization and primality testing" (1974): p − 1
 * - M. Nemec et al., "The Return of Coppersmith's Attack" (CCS 2017): ROCA
 * - N. Heninger et al., "Mining Your Ps and Qs" (USENIX Security 2012): batch GCD
 *
 * @example
 * ```typescript
 * import { checkRsaKey } from '@zksecurity/cheatsheets';
 *
 * for (const f of checkRsaKey({ n, e: 65537n })) console.log(f.check, f.severity, f.message);
 * ```
 */

import {
  continued_fraction,
  convergents,
  gcd,
  isqrt,
  lcm,
  power_mod,
} from '../../../sagemath-ts/src/arith/misc.js';
import { ValueError } from '../../../sagemath-ts/src/errors.js';
import { discrete_log, order_from_multiple } from '../../../sagemath-ts/src/groups/generic.js';
import { IntegerModRing } from '../../../sagemath-ts/src/rings/finite_rings/integer_mod_ring.js';
import { type RsaBatchFinding, findSharedNMinus1Factors, findSharedPrimes } from './rsa-batch.js';
import {
  bitLength,
  exactSqrt,
  fermatSearch,
  gcdReduced,
  primeProduct,
  primesUpTo,
} from './rsa-math.js';

// =============================================================================
// Types
// =============================================================================

export interface RsaPublicKey {
  n: bigint;
  e: bigint;
}

/** `broken`: the private key was recovered. `vulnerable`: a known flaw, not exploited here. */
export type RsaSeverity = 'broken' | 'vulnerable' | 'warning';

export interface RsaFinding {
  check: string;
  severity: RsaSeverity;
  message: string;
  /** Set when the check factored n = p·q, with p ≤ q. */
  p?: bigint;
  q?: bigint;
  /** Set when the check recovered the private exponent. */
  d?: bigint;
}

export interface RsaCheckOptions {
  /** Plain Fermat steps. Default 2^20. */
  fermatIterations?: number;
  /** Fermat steps for each guessed count of shared low bits. Default 2^12. */
  sharedBitsIterations?: number;
  /** Largest prime tried by the small-factor check. Default 2^20. */
  smallFactorBound?: number;
  /** Pollard p − 1 stage-1 bound. Default 100 000. */
  pm1Bound?: number;
  /** Pollard p − 1 stage-2 bound. Default 1 000 000. */
  pm1Bound2?: number;
  /** Largest denominator, in bits, of the rational-approximation check. Default 128. */
  maxDenominatorBits?: number;
  /** Checks to run. Default `RSA_KEY_CHECKS`. */
  checks?: readonly RsaKeyCheck[];
}

export interface RsaKeyCheck {
  name: string;
  run(key: RsaPublicKey, options: RsaCheckOptions): Omit<RsaFinding, 'check'> | null;
}

// =============================================================================
// Attacks
// =============================================================================

function ordered(a: bigint, b: bigint): [bigint, bigint] {
  return a <= b ? [a, b] : [b, a];
}

function positiveMod(x: bigint, m: bigint): bigint {
  return ((x % m) + m) % m;
}

/** x with x² ≡ n (mod 2^bits), for n ≡ 1 (mod 8). */
function twoAdicSqrt(n: bigint, bits: number): bigint {
  const low = n & ((1n << BigInt(bits)) - 1n);
  let x = 1n;
  for (let k = 3; k < bits; k++) {
    if (((x * x - low) >> BigInt(k)) & 1n) x += 1n << BigInt(k - 1);
  }
  return x;
}

/** Residues of (p + q)/2 modulo 2^(l−1) when p ≡ q (mod 2^l), given x² ≡ n (mod 2^l). */
function halfSumClasses(n: bigint, lowBits: number, x: bigint): bigint[] {
  if (lowBits < 2) return [0n];
  if (lowBits === 2) return n % 4n === 1n ? [1n] : [];
  if (n % 8n !== 1n) return [];
  const m = 1n << BigInt(lowBits - 1);
  const r = x % m;
  return r === m - r ? [r] : [r, m - r];
}

/**
 * Fermat factorization: finds n = p·q when |p − q| is small, in about
 * (p − q)² / (8√n) steps.
 *
 * With `lowBits = l`, assumes p ≡ q (mod 2^l) and steps (p + q)/2 by 2^(l−1), so primes
 * sharing h high bits and l low bits of b-bit halves fall in about 2^(b − 2h − l) steps.
 * A list of guesses is tried in order, `maxIterations` steps each.
 */
export function fermatFactor(
  n: bigint,
  options: { maxIterations?: number; lowBits?: number | readonly number[] } = {}
): [bigint, bigint] | null {
  const { maxIterations = 1 << 20, lowBits = 0 } = options;
  if (n < 9n || (n & 1n) === 0n) return null;
  const root = isqrt(n);
  if (root * root === n) return [root, root];
  const guesses = typeof lowBits === 'number' ? [lowBits] : lowBits;
  const maxLowBits = Math.max(0, ...guesses);
  const x = maxLowBits >= 3 && n % 8n === 1n ? twoAdicSqrt(n, maxLowBits) : 1n;
  for (const l of guesses) {
    const step = l >= 2 ? 1n << BigInt(l - 1) : 1n;
    for (const c of halfSumClasses(n, l, x)) {
      const a0 = root + 1n + positiveMod(c - root - 1n, step);
      const found = fermatSearch(n, a0, step, maxIterations, (a, b) =>
        a - b > 1n ? ordered(a - b, a + b) : null
      );
      if (found) return found;
    }
  }
  return null;
}

/** Distinct primes up to `bound` dividing n, smallest first. */
export function smallFactors(n: bigint, bound = 1 << 20): bigint[] {
  let g = gcdReduced(n, primeProduct(bound));
  const found: bigint[] = [];
  for (const r of primesUpTo(bound)) {
    if (g === 1n) break;
    if (g % r === 0n) {
      found.push(r);
      g /= r;
    }
  }
  return found;
}

/** Wiener's attack: recovers d < n^(1/4)/3 from the continued fraction of e/n. */
export function wienerAttack(n: bigint, e: bigint): { d: bigint; p: bigint; q: bigint } | null {
  for (const [k, d] of convergents(continued_fraction(e, n))) {
    if (k === 0n) continue;
    const ed1 = e * d - 1n;
    if (ed1 % k !== 0n) continue;
    const s = n - ed1 / k + 1n;
    const r = exactSqrt(s * s - 4n * n);
    if (r === null || ((s + r) & 1n) === 1n) continue;
    const [p, q] = ordered((s - r) / 2n, (s + r) / 2n);
    if (p > 1n && p * q === n) return { d, p, q };
  }
  return null;
}

export interface RationalApproximation {
  /** n ≈ (numerator / denominator) · 2^exponent. */
  numerator: bigint;
  denominator: bigint;
  exponent: number;
  /** Bits of agreement beyond what a random n gives at this denominator. */
  precisionBits: number;
}

/**
 * Finds n ≈ (u/v)·2^k with small v far closer than chance allows. Primes with repeating
 * bit patterns are close to a small rational times a power of two, and so is their product.
 * A random n reaches `minPrecisionBits` = 64 with probability below 2^−56.
 */
export function rationalApproximation(
  n: bigint,
  options: { maxDenominatorBits?: number; minPrecisionBits?: number } = {}
): RationalApproximation | null {
  const { maxDenominatorBits = 128, minPrecisionBits = 64 } = options;
  const exponent = bitLength(n);
  const terms = continued_fraction(
    n,
    1n << BigInt(exponent),
    Math.ceil(1.5 * maxDenominatorBits) + 3
  );
  const convs = convergents(terms);
  const limit = 1n << BigInt(maxDenominatorBits);
  let best: RationalApproximation | null = null;
  for (let i = 1; i + 1 < terms.length; i++) {
    const [numerator, denominator] = convs[i]!;
    if (denominator > limit) break;
    const precisionBits = bitLength(terms[i + 1]!);
    if (
      precisionBits >= minPrecisionBits &&
      (best === null || precisionBits > best.precisionBits)
    ) {
      best = { numerator, denominator, exponent, precisionBits };
    }
  }
  return best;
}

/**
 * Factors n from a rational approximation (Lehman's method with a guessed multiplier). If
 * p ≈ x·2^i and q ≈ y·2^j for small rationals x, y, then p/q ≈ s/t, and uv·w²·z = st·d²
 * for some small w, z ∈ {1, 2} and integer d. Fermat on 4·uv·w²·z·n then finds (2dtp)(2dsq).
 */
export function rationalApproximationFactor(
  n: bigint,
  approximation: RationalApproximation,
  options: { maxCofactor?: number; iterations?: number } = {}
): [bigint, bigint] | null {
  const { maxCofactor = 256, iterations = 256 } = options;
  const uv = approximation.numerator * approximation.denominator;
  for (let w = 1n; w <= BigInt(maxCofactor); w++) {
    for (const z of [1n, 2n]) {
      const N = 4n * uv * w * w * z * n;
      const found = fermatSearch(N, isqrt(N), 1n, iterations, (a, b) => {
        const g = gcd(a + b, n);
        return g > 1n && g < n ? ordered(g, n / g) : null;
      });
      if (found) return found;
    }
  }
  return null;
}

function primePower(r: bigint, bound: bigint): bigint {
  let power = r;
  while (power * r <= bound) power *= r;
  return power;
}

/**
 * Pollard's p − 1: splits n when p − 1 is `bound`-smooth apart from at most one prime
 * up to `bound2`.
 */
export function pollardPm1(
  n: bigint,
  options: { bound?: number; bound2?: number } = {}
): [bigint, bigint] | null {
  const { bound = 100_000, bound2 = 1_000_000 } = options;
  if ((n & 1n) === 0n) return n > 2n ? [2n, n / 2n] : null;
  const split = (g: bigint) => (g > 1n && g < n ? ordered(g, n / g) : null);
  const primes = primesUpTo(Math.max(bound, bound2));
  const B1 = BigInt(bound);
  const B2 = BigInt(bound2);

  let a = 2n;
  let i = 0;
  while (i < primes.length && primes[i]! <= B1) {
    const start = i;
    const saved = a;
    let exponent = 1n;
    for (; i < primes.length && primes[i]! <= B1 && i - start < 256; i++) {
      exponent *= primePower(primes[i]!, B1);
    }
    a = power_mod(a, exponent, n);
    const g = gcd(a - 1n, n);
    if (g === 1n) continue;
    if (g !== n) return split(g);
    a = saved;
    for (let j = start; j < i; j++) {
      a = power_mod(a, primePower(primes[j]!, B1), n);
      const h = gcd(a - 1n, n);
      if (h !== 1n) return split(h);
    }
  }

  const gapPowers = new Map<bigint, bigint>();
  const gapPower = (gap: bigint) => {
    let v = gapPowers.get(gap);
    if (v === undefined) {
      v = power_mod(a, gap, n);
      gapPowers.set(gap, v);
    }
    return v;
  };
  const next = (x: bigint, j: number) =>
    j + 1 < primes.length ? (x * gapPower(primes[j + 1]! - primes[j]!)) % n : x;
  let x = i < primes.length ? power_mod(a, primes[i]!, n) : 1n;
  while (i < primes.length && primes[i]! <= B2) {
    const start = i;
    const saved = x;
    let acc = 1n;
    for (; i < primes.length && primes[i]! <= B2 && i - start < 1024; i++) {
      acc = (acc * (x - 1n)) % n;
      x = next(x, i);
    }
    const g = gcd(acc, n);
    if (g === 1n) continue;
    if (g !== n) return split(g);
    x = saved;
    for (let j = start; j < i; j++) {
      const h = gcd(x - 1n, n);
      if (h !== 1n) return split(h);
      x = next(x, j);
    }
  }
  return null;
}

const ROCA_GENERATOR = 65537n;
let rocaGroup: { ring: IntegerModRing; order: bigint } | undefined;

function getRocaGroup(): { ring: IntegerModRing; order: bigint } {
  if (rocaGroup === undefined) {
    // Infineon's M is a product of the first 39, 71, 126 or 225 primes by key size; the
    // first 39 (up to 167) divide every one of them.
    let modulus = 1n;
    let order = 1n;
    for (const r of primesUpTo(167)) {
      modulus *= r;
      const g = new IntegerModRing(r).__call__(ROCA_GENERATOR);
      order = lcm(order, order_from_multiple(g, r - 1n, undefined, '*'));
    }
    rocaGroup = { ring: new IntegerModRing(modulus), order };
  }
  return rocaGroup;
}

/**
 * ROCA fingerprint (CVE-2017-15361): Infineon RSALib primes are k·M + (65537^a mod M), so
 * n mod M lies in the subgroup generated by 65537. Returns c with n ≡ 65537^c (mod M),
 * or null when no such c exists. A random n passes with probability about 2^−154.
 */
export function rocaDiscreteLog(n: bigint): bigint | null {
  const { ring, order } = getRocaGroup();
  const modulus = ring.modulus;
  if (gcd(n % modulus, modulus) !== 1n || power_mod(n, order, modulus) !== 1n) return null;
  try {
    return discrete_log(ring.__call__(n), ring.__call__(ROCA_GENERATOR), order, '*');
  } catch (e) {
    if (e instanceof ValueError) return null;
    throw e;
  }
}

// =============================================================================
// Checks
// =============================================================================

function broken(message: string, [p, q]: [bigint, bigint]): Omit<RsaFinding, 'check'> {
  return { severity: 'broken', message, p, q };
}

export const RSA_KEY_CHECKS: readonly RsaKeyCheck[] = [
  {
    name: 'size',
    run: ({ n }) => {
      const bits = bitLength(n);
      return bits < 2048
        ? { severity: 'warning', message: `${bits}-bit modulus is below 2048 bits` }
        : null;
    },
  },
  {
    name: 'exponent',
    run: ({ n, e }) => {
      if (e < 3n || (e & 1n) === 0n || e >= n) {
        return { severity: 'warning', message: `e = ${e} is not a valid RSA public exponent` };
      }
      return e === 65537n ? null : { severity: 'warning', message: `e = ${e} is not 65537` };
    },
  },
  {
    name: 'small-factors',
    run: ({ n }, options) => {
      const factors = smallFactors(n, options.smallFactorBound).filter((r) => r < n);
      return factors.length === 0
        ? null
        : broken(
            `n has small prime factors ${factors.join(', ')}`,
            ordered(factors[0]!, n / factors[0]!)
          );
    },
  },
  {
    name: 'fermat',
    run: ({ n }, options) => {
      const factors = fermatFactor(n, { maxIterations: options.fermatIterations });
      return factors === null
        ? null
        : broken(
            `p and q are close (|p − q| has ${bitLength(factors[1] - factors[0])} bits): Fermat factors n`,
            factors
          );
    },
  },
  {
    name: 'fermat-shared-bits',
    run: ({ n }, options) => {
      const lowBits: number[] = [];
      for (let l = Math.floor(bitLength(n) / 16) * 8; l >= 8; l -= 8) lowBits.push(l);
      const factors = fermatFactor(n, {
        maxIterations: options.sharedBitsIterations ?? 1 << 12,
        lowBits,
      });
      if (factors === null) return null;
      const diff = factors[1] - factors[0];
      const low = diff === 0n ? bitLength(factors[0]) : bitLength(diff & -diff) - 1;
      const high = bitLength(factors[0]) - bitLength(diff);
      return broken(
        `p and q share their top ${high} and bottom ${low} bits: Fermat factors n`,
        factors
      );
    },
  },
  {
    name: 'wiener',
    run: ({ n, e }) => {
      const result = wienerAttack(n, e);
      return result === null
        ? null
        : {
            severity: 'broken',
            message: `the private exponent d has ${bitLength(result.d)} bits: Wiener's attack recovers it`,
            p: result.p,
            q: result.q,
            d: result.d,
          };
    },
  },
  {
    name: 'rational-approximation',
    run: ({ n }, options) => {
      const approx = rationalApproximation(n, { maxDenominatorBits: options.maxDenominatorBits });
      if (approx === null) return null;
      const fit = `n ≈ (${approx.numerator}/${approx.denominator})·2^${approx.exponent} to ${approx.precisionBits} bits beyond chance`;
      const factors = rationalApproximationFactor(n, approx);
      return factors === null
        ? {
            severity: 'vulnerable',
            message: `${fit}: the primes likely have repeating bit patterns`,
          }
        : broken(
            `${fit}: the primes have repeating bit patterns and Fermat with a multiplier factors n`,
            factors
          );
    },
  },
  {
    name: 'pollard-p-1',
    run: ({ n }, options) => {
      const factors = pollardPm1(n, { bound: options.pm1Bound, bound2: options.pm1Bound2 });
      return factors === null ? null : broken('p − 1 is smooth: Pollard p − 1 factors n', factors);
    },
  },
  {
    name: 'roca',
    run: ({ n }) => {
      const c = rocaDiscreteLog(n);
      return c === null
        ? null
        : {
            severity: 'vulnerable',
            message: `ROCA (CVE-2017-15361): n ≡ 65537^${c} modulo the Infineon primorial, so the primes come from Infineon RSALib and n falls to Coppersmith's method`,
          };
    },
  },
  // Coppersmith-based checks (patterned primes, partially known primes) belong here.
];

/** Runs every check on one key. Takes well under a second for a 2048-bit key. */
export function checkRsaKey(key: RsaPublicKey, options: RsaCheckOptions = {}): RsaFinding[] {
  if (key.n < 2n) throw new Error('n must be at least 2');
  const findings: RsaFinding[] = [];
  for (const check of options.checks ?? RSA_KEY_CHECKS) {
    const finding = check.run(key, options);
    if (finding !== null) findings.push({ check: check.name, ...finding });
  }
  return findings;
}

export interface RsaBatchReport {
  /** Findings for `keys[i]`, from `checkRsaKey`. */
  keys: RsaFinding[][];
  /** Findings that involve several keys at once. */
  batch: RsaBatchFinding[];
}

/** Runs the batch checks on all keys and, unless `perKey` is false, `checkRsaKey` on each. */
export function checkRsaKeys(
  keys: readonly RsaPublicKey[],
  options: RsaCheckOptions & { perKey?: boolean } = {}
): RsaBatchReport {
  const moduli = keys.map((k) => k.n);
  return {
    keys: keys.map((k) => (options.perKey === false ? [] : checkRsaKey(k, options))),
    batch: [...findSharedPrimes(moduli), ...findSharedNMinus1Factors(moduli)],
  };
}
