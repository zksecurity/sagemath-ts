import { pariErrorPayload } from './_error_display.js';
/**
 * @module ff
 * @description Finite field arithmetic (Fp operations)
 *
 * This file ports the Fp_ functions from PARI/GP.
 * Reference: reference/pari/src/headers/pariinl.h (lines 1668-1850)
 *            reference/pari/src/basemath/arith1.c (Tonelli-Shanks, kronecker)
 *            reference/pari/src/basemath/FpV.c (vector operations)
 */

import { scalarPower } from './_scalar_power.js';
import { modularSquareRoot, wordSquareRoot } from './_modular_sqrt.js';
import { Z_factor, factoru } from './ifactor.js';
import { mod4, mod8, vali } from './types.js';
import { inverseCoefficient, residue } from './_polynomial_division.js';
import { PariError } from './errors.js';

// ============================================================================
// Core Fp operations
// Reference: pariinl.h lines 1668-1850
// ============================================================================

/**
 * Fp_red - Reduce a mod p to [0, p-1]
 *
 * Reference: pariinl.h:1668
 * GEN Fp_red(GEN a, GEN m) { return modii(a, m); }
 * @see Deviation: PARI exported scalar arithmetic boundaries
 */
export function Fp_red(a: bigint, p: bigint): bigint {
  return residue(a, p);
}

/**
 * Fp_add - Addition mod p
 *
 * Reference: pariinl.h:1670-1690
 *
 * PARI takes the fast path when 0 <= a, b < p, but always falls back to a
 * full `remii`/`modii` reduction otherwise, so the result is in [0, p) for
 * arbitrary integer inputs. We mirror that exactly.
 * @see Deviation: PARI exported scalar arithmetic boundaries
 */
export function Fp_add(a: bigint, b: bigint, p: bigint): bigint {
  const sum = a + b;
  if (sum === 0n) return 0n;
  if (sum > 0n) {
    const t = sum - p;
    if (t === 0n) return 0n;
    if (t < 0n) return sum;
    if (t < p) return t; /* general case ! */
    return residue(t, p); /* t > 0: native remii has a nonnegative remainder. */
  }
  return Fp_red(sum, p);
}

/**
 * Fp_sub - Subtraction mod p
 *
 * Reference: pariinl.h:1696-1714
 *
 * As for Fp_add: PARI reduces fully when the fast path does not apply.
 * @see Deviation: PARI exported scalar arithmetic boundaries
 */
export function Fp_sub(a: bigint, b: bigint, p: bigint): bigint {
  const diff = a - b;
  if (diff === 0n) return 0n;
  if (diff > 0n) {
    if (diff < p) return diff; /* general case ! */
    return residue(diff, p);
  }
  return Fp_red(diff + p, p);
}

/**
 * Fp_neg - Negation mod p
 *
 * Reference: pariinl.h:1717-1730
 *
 * As for Fp_add: PARI reduces fully when the fast path does not apply.
 * @see Deviation: PARI exported scalar arithmetic boundaries
 */
export function Fp_neg(a: bigint, p: bigint): bigint {
  if (a === 0n) return 0n;
  if (a > 0n) {
    const t = p - a;
    if (t >= 0n) return t; /* general case ! */
    return Fp_red(t, p);
  }
  return residue(-a, p); /* -a > 0: native remii has a nonnegative remainder. */
}

/**
 * Fp_mul - Multiplication mod p
 *
 * Reference: pariinl.h:1761-1768
 * Computes (a * b) mod p
 * @see Deviation: PARI exported scalar arithmetic boundaries
 */
export function Fp_mul(a: bigint, b: bigint, p: bigint): bigint {
  return residue(a * b, p);
}

/**
 * Fp_sqr - Square mod p
 *
 * Reference: pariinl.h:1770-1777
 * Computes a^2 mod p
 * @see Deviation: PARI exported scalar arithmetic boundaries
 */
export function Fp_sqr(a: bigint, p: bigint): bigint {
  if (p === 0n) throw new PariError('impossible inverse in dvmdii: 0.');
  return residue(a * a, p);
}

/**
 * Fp_inv - Modular inverse using extended Euclidean algorithm
 *
 * Reference: pariinl.h:1818-1823
 * Computes a^(-1) mod p
 * Reports the native Fp_inv error when invmod fails, including modulus zero.
 * @see Deviation: PARI exported scalar arithmetic boundaries
 */
export function Fp_inv(a: bigint, p: bigint): bigint {
  return inverseCoefficient(a, p, false);
}

/** Native unsigned-divisor branch; a is a nonzero word magnitude.
 * @see Deviation: PARI exported scalar arithmetic boundaries
 */
function Fp_divu(x: bigint, a: bigint, p: bigint): bigint {
  const pp = p < 0n ? -p : p;
  if (pp > 0n && pp < 1n << 64n) {
    x = residue(x, pp);
    if (x === 0n) return 0n;
    return residue(x * inverseCoefficient(a % pp, pp, true), pp);
  }
  x = residue(x, p);
  const t = residue(residue(x, a) * inverseCoefficient(residue(p, a), a, true), a);
  const b = t === 0n ? 0n : a - t;
  return (x + p * b) / a;
}

/**
 * Fp_div - Division mod p
 *
 * Reference: pariinl.h:1832-1846
 * Native word-divisor and generic inverse paths preserve distinct shortcuts.
 * @see Deviation: PARI exported scalar arithmetic boundaries
 */
export function Fp_div(a: bigint, b: bigint, p: bigint): bigint {
  const divisor = b < 0n ? -b : b;
  if (divisor > 0n && divisor < 1n << 64n) {
    const result = Fp_divu(a, divisor, p);
    return b < 0n ? Fp_neg(result, p) : result;
  }
  return Fp_mul(a, Fp_inv(b, p), p);
}

/**
 * Fp_pow - Integer modular powering with PARI's native dispatch and reductions.
 * Reference: basemath/arith1.c:Fp_pow, Fp_powu, Fp_select_red
 * @see Deviation: PARI scalar exponentiation boundaries
 */
export function Fp_pow(a: bigint, n: bigint, p: bigint): bigint {
  return scalarPower(a, n, p);
}

// ============================================================================
// Kronecker/Legendre symbol
// Reference: arith1.c:361-574
// ============================================================================

/**
 * Helper: t = 3,5 mod 8? (= 2 not a square mod t)
 * Reference: arith1.c:362-371
 */
function ome(t: number): boolean {
  const r = t & 7;
  return r === 3 || r === 5;
}

/**
 * Helper for bigint version of ome
 */
function gome(t: bigint): boolean {
  if (t === 0n) return false;
  return ome(Number(((t % 8n) + 8n) % 8n));
}

/**
 * krouu_s - Kronecker symbol for unsigned integers with sign accumulator
 *
 * Reference: arith1.c:378-394
 * Assume y odd, return kronecker(x,y) * s
 */
function krouu_s(x: bigint, y: bigint, s: number): number {
  if (x < 0n) x = ((x % y) + y) % y;
  x = x % y;

  while (x !== 0n) {
    // Extract powers of 2 from x
    let r = 0;
    while ((x & 1n) === 0n) {
      r++;
      x >>= 1n;
    }
    if (r & 1) {
      // odd power of 2
      if (gome(y)) s = -s;
    }
    // Quadratic reciprocity: (x|y)(y|x) = (-1)^((x-1)/2 * (y-1)/2)
    // Both x and y are odd here
    if ((x & 2n) !== 0n && (y & 2n) !== 0n) s = -s;
    const z = y % x;
    y = x;
    x = z;
  }

  return y === 1n ? s : 0;
}

/**
 * kronecker - Kronecker symbol (x|y)
 *
 * Reference: arith1.c:396-447
 *
 * The Kronecker symbol generalizes the Jacobi symbol.
 * For odd prime p: (a|p) = Legendre symbol
 * Returns: 1 if a is a quadratic residue mod p
 *         -1 if a is a non-residue
 *          0 if gcd(a,p) > 1
 */
export function kronecker(x: bigint, y: bigint): number {
  let s = 1;

  // Handle sign of y
  if (y < 0n) {
    y = -y;
    if (x < 0n) s = -1;
  }

  // y = 0 case
  if (y === 0n) {
    return x === 1n || x === -1n ? 1 : 0;
  }

  // Extract powers of 2 from y
  let r = 0;
  while ((y & 1n) === 0n) {
    r++;
    y >>= 1n;
  }

  if (r > 0) {
    // x must be odd for (x|2^r) to be nonzero
    if ((x & 1n) === 0n) return 0;
    // (x|2) = (-1)^((x^2-1)/8) = 1 if x = +/-1 mod 8, -1 if x = +/-3 mod 8
    if ((r & 1) !== 0 && gome(x)) s = -s;
  }

  // Now y is odd and positive
  // Reduce x mod y
  x = ((x % y) + y) % y;

  // Main loop - quadratic reciprocity
  return krouu_s(x, y, s);
}

/**
 * Fp_issquare - PARI's characteristic-two or Kronecker predicate.
 * Reference: basemath/FpX.c:Fp_issquare. The symbol may be zero at composites.
 * @see Deviation: PARI scalar predicates and order boundaries
 */
export function Fp_issquare(a: bigint, p: bigint): boolean {
  return p === 2n || p === -2n || kronecker(a, p) !== -1;
}

// ============================================================================
// Square root in Fp - Tonelli-Shanks algorithm
// Reference: arith1.c:762-858, 1189-1258
// ============================================================================

/** Square root modulo a prime, choosing the native smallest root.
 * Delegates to the word, signed-small, Atkin, Cipolla or Tonelli-Shanks kernel.
 * Reference: arith1.c:Fp_sqrt / Fp_sqrt_i.
 * @see Deviation: PARI modular square-root adapters
 */
export function Fp_sqrt(a: bigint, p: bigint): bigint | null {
  return modularSquareRoot(a, null, p);
}

/** As Fp_sqrt, with an optional generator of the 2-Sylow subgroup.
 * The native word route ignores y; null requests deterministic generator search.
 */
export function Fp_sqrt_i(a: bigint, y: bigint | null, p: bigint): bigint | null {
  return modularSquareRoot(a, y, p);
}

/** Reduced word-prime square root; null represents the native ULONG_MAX sentinel.
 * @see Deviation: PARI modular square-root adapters
 */
export function Fl_sqrt(a: bigint, p: bigint): bigint | null {
  return wordSquareRoot(a, p);
}

// ============================================================================
// Additional utility functions
// ============================================================================

/**
 * Fp_center - Center an already reduced coefficient by magnitude comparison.
 *
 * Reference: pariinl.h (Fp_center)
 * @see Deviation: PARI exported scalar arithmetic boundaries
 */
export function Fp_center(a: bigint, p: bigint): bigint {
  const half = (p < 0n ? -p : p) / 2n;
  return (a < 0n ? -a : a) <= half ? a : a - p;
}

/**
 * Fp_mulu - Multiply by a small unsigned integer
 * @see Deviation: PARI exported scalar arithmetic boundaries
 */
export function Fp_mulu(a: bigint, b: number, p: bigint): bigint {
  const pp = p < 0n ? -p : p;
  if (pp > 0n && pp < 1n << 64n) return residue(residue(a, pp) * BigInt(b), pp);
  return Fp_mul(a, BigInt(b), p);
}

/**
 * Fp_addmul - Compute x + y*z mod p
 *
 * Reference: pariinl.h:1751-1758
 * @see Deviation: PARI exported scalar arithmetic boundaries
 */
export function Fp_addmul(x: bigint, y: bigint, z: bigint, p: bigint): bigint {
  if (y === 0n || z === 0n) return Fp_red(x, p);
  if (x === 0n) return Fp_mul(z, y, p);
  return residue(x + y * z, p);
}

/**
 * Fp_double - Double an already reduced coefficient with one subtraction.
 * @see Deviation: PARI exported scalar arithmetic boundaries
 */
export function Fp_double(a: bigint, p: bigint): bigint {
  const doubled = a * 2n;
  return doubled >= p ? doubled - p : doubled;
}

/**
 * Fp_halve - Compute a/2 mod p
 *
 * If a is even, just divide. If a is odd, compute (a+p)/2.
 */
export function Fp_halve(a: bigint, p: bigint): bigint {
  if ((a & 1n) === 0n) {
    return a / 2n;
  }
  return (a + p) / 2n;
}

/**
 * Fp_eq - Check equality of two Fp elements
 */
export function Fp_eq(a: bigint, b: bigint): boolean {
  return a === b;
}

/**
 * GCD using Euclidean algorithm
 */
export function gcd(a: bigint, b: bigint): bigint {
  if (a < 0n) a = -a;
  if (b < 0n) b = -b;
  while (b !== 0n) {
    const t = b;
    b = a % b;
    a = t;
  }
  return a;
}

/**
 * Extended GCD: returns [g, x, y] such that a*x + b*y = g = gcd(a,b).
 * PARI bezout(0,0) chooses both coefficients zero.
 * @see Deviation: PARI scalar predicates and order boundaries
 */
export function xgcd(a: bigint, b: bigint): [bigint, bigint, bigint] {
  if (a === 0n && b === 0n) return [0n, 0n, 0n];
  let [old_r, r] = [a, b];
  let [old_s, s] = [1n, 0n];
  let [old_t, t] = [0n, 1n];

  while (r !== 0n) {
    const quotient = old_r / r;
    [old_r, r] = [r, old_r - quotient * r];
    [old_s, s] = [s, old_s - quotient * s];
    [old_t, t] = [t, old_t - quotient * t];
  }

  // Make gcd positive
  if (old_r < 0n) {
    return [-old_r, -old_s, -old_t];
  }
  return [old_r, old_s, old_t];
}

/**
 * Order of a nonzero element, given a positive multiple N of its order.
 * arith1.c:2603-2635 (Fl_order/Fp_order), bb_group.c:670-713 (rec_order).
 * Uses PARI's machine-word path below 2^64 and its recursive factor splitting
 * for larger moduli. Inputs must satisfy a^N = 1 modulo p.
 * Word moduli use |N| when it fits a nonzero word, otherwise default to |p|-1.
 * The generic branch requires N > 0 and retains the original unreduced base.
 * @see Deviation: PARI scalar predicates and order boundaries
 */
export function Fp_order(a: bigint, N: bigint, p: bigint): bigint {
  const pp = p < 0n ? -p : p;
  if (pp > 0n && pp < 1n << 64n) {
    a = Fp_red(a, pp);
    if (a === 1n) return 1n;
    N = N < 0n ? -N : N;
    if (N === 0n || N >= 1n << 64n) N = pp - 1n;
    p = pp;
    const factors = factoru(N).filter(([q]) => q > 1n);
    let order = N;
    for (let i = factors.length - 1; i >= 0; i--) {
      const [q, e] = factors[i]!;
      const t = order / q ** e;
      let y = Fp_pow(a, t, p);
      if (y === 1n) order = t;
      else
        for (let j = 1n; j < e; j++) {
          y = Fp_pow(y, q, p);
          if (y === 1n) {
            order = t * q ** j;
            break;
          }
        }
    }
    return order;
  }
  if (N <= 0n)
    throw new PariError('incorrect type in generic discrete logarithm (order factorization) (t_INT).');
  const factors = Z_factor(N).filter(([q]) => q > 1n);
  if (N === 1n) return 1n;
  function recOrder(x: bigint, lo: number, hi: number): bigint {
    if (x === 1n) return 1n;
    if (lo === hi) {
      const [q, e] = factors[lo]!;
      let result = 1n;
      for (let i = 0n; i < e; i++) {
        if (x === 1n) return result;
        x = Fp_pow(x, q, p);
        result *= q;
      }
      return result;
    }
    const mid = Math.floor((lo + hi) / 2);
    let cofactor = 1n;
    for (let i = lo; i <= mid; i++) {
      const [q, e] = factors[i]!;
      cofactor *= q ** e;
    }
    const order1 = recOrder(Fp_pow(x, cofactor, p), mid + 1, hi);
    return order1 * recOrder(Fp_pow(x, order1, p), lo, mid);
  }
  return recOrder(a, 0, factors.length - 1);
}

/** PARI arith1.c:2639-2692: order of Mod(a,n), optionally given a multiple.
 * @see Deviation: PARI scalar predicates and order boundaries
 */
export function znorder(a: bigint, n: bigint, order?: bigint): bigint {
  n = n < 0n ? -n : n;
  if (n === 0n) throw new PariError('impossible inverse in %: 0.');
  a = Fp_red(a, n);
  if (gcd(a, n) !== 1n)
    throw new PariError(
      'elements not coprime in znorder:\n    ' + pariErrorPayload(String(a), 't_INT') +
      '\n    ' + pariErrorPayload(String(n), 't_INT')
    );
  if (order !== undefined) return Fp_order(a, order, n);
  let result = 1n;
  for (const [p, e] of Z_factor(n)) {
    if (p <= 1n) continue;
    const pe = p ** e;
    let residue = Fp_red(a, pe);
    let local: bigint;
    // Zp_order: strip the residue-field order, then use the p-adic
    // valuation to determine the exact remaining prime-power factor.
    if (p === 2n) {
      if (e === 1n) local = 1n;
      else if (e === 2n) local = residue % 4n === 1n ? 1n : 2n;
      else {
        local = residue % 4n === 1n ? 1n : 2n;
        if (local === 2n) residue = Fp_sqr(residue, pe);
        if (residue !== 1n) {
          let difference = residue - 1n;
          let valuation = 0n;
          while (difference % p === 0n) {
            difference /= p;
            valuation++;
          }
          local *= p ** (e - valuation);
        }
      }
    } else {
      local = Fp_order(Fp_red(residue, p), p - 1n, p);
      if (e > 1n) {
        residue = Fp_pow(residue, local, pe);
        if (residue !== 1n) {
          let difference = residue - 1n;
          let valuation = 0n;
          while (difference % p === 0n) {
            difference /= p;
            valuation++;
          }
          local *= p ** (e - valuation);
        }
      }
    }
    result = (result / gcd(result, local)) * local;
  }
  return result;
}

import type { PariFfelt } from './types.js';
import { finiteFieldSquareRoot } from './_finite_field_square_root.js';
/** FF_issquareall with the root output pointer present; null means no root.
 * The element uses a reduced polynomial representation over an actual field.
 * @see Deviation: PARI finite-field scalar square-root adapter
 */
export function FF_issquareall(x: PariFfelt): PariFfelt | null {
  return finiteFieldSquareRoot(x);
}

/** Native finite-field square predicate; this call does not sample randomness.
 * @see Deviation: PARI finite-field norm and square-predicate adapters
 */
export function FF_issquare(x: PariFfelt): boolean {
  if (x.p===2n) return true;
  const a=typeof x.value==='bigint'?[x.value]:[...x.value], T=[...(x.definingPoly??[0n,1n])];
  return x.p<1n<<64n ? Flxq_issquare(a,T,x.p):FpXQ_issquare(a,T,x.p);
}
/** Native field norm down to its prime field, returned as a bigint residue.
 * @see Deviation: PARI finite-field norm and square-predicate adapters
 */
export function FF_norm(x: PariFfelt): bigint {
  const a=typeof x.value==='bigint'?[x.value]:[...x.value], T=[...(x.definingPoly??[0n,1n])];
  if (x.p===2n) return a.some(c=>c!==0n)?1n:0n;
  return x.p<1n<<64n ? Flxq_norm(a,T,x.p):FpXQ_norm(a,T,x.p);
}

import { Flxq_issquare, Flxq_norm } from './Flx.js';
import { FpXQ_issquare, FpXQ_norm } from './FpX.js';
