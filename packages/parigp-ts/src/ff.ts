/**
 * @module ff
 * @description Finite field arithmetic (Fp operations)
 *
 * This file ports the Fp_ functions from PARI/GP.
 * Reference: reference/pari/src/headers/pariinl.h (lines 1668-1850)
 *            reference/pari/src/basemath/arith1.c (Tonelli-Shanks, kronecker)
 *            reference/pari/src/basemath/FpV.c (vector operations)
 */

import { mod4, mod8, vali } from './types.js';

// ============================================================================
// Core Fp operations
// Reference: pariinl.h lines 1668-1850
// ============================================================================

/**
 * Fp_red - Reduce a mod p to [0, p-1]
 *
 * Reference: pariinl.h:1668
 * GEN Fp_red(GEN a, GEN m) { return modii(a, m); }
 */
export function Fp_red(a: bigint, p: bigint): bigint {
  const r = a % p;
  return r < 0n ? r + p : r;
}

/**
 * Fp_add - Addition mod p
 *
 * Reference: pariinl.h:1670-1693
 * Computes (a + b) mod p, assuming 0 <= a, b < p
 */
export function Fp_add(a: bigint, b: bigint, p: bigint): bigint {
  const sum = a + b;
  // If sum >= p, subtract p
  return sum >= p ? sum - p : sum;
}

/**
 * Fp_sub - Subtraction mod p
 *
 * Reference: pariinl.h:1696-1714
 * Computes (a - b) mod p, assuming 0 <= a, b < p
 */
export function Fp_sub(a: bigint, b: bigint, p: bigint): bigint {
  const diff = a - b;
  // If diff < 0, add p
  return diff < 0n ? diff + p : diff;
}

/**
 * Fp_neg - Negation mod p
 *
 * Reference: pariinl.h:1717-1730
 * Computes (-a) mod p, assuming 0 <= a < p
 */
export function Fp_neg(a: bigint, p: bigint): bigint {
  return a === 0n ? 0n : p - a;
}

/**
 * Fp_mul - Multiplication mod p
 *
 * Reference: pariinl.h:1761-1768
 * Computes (a * b) mod p
 */
export function Fp_mul(a: bigint, b: bigint, p: bigint): bigint {
  return (a * b) % p;
}

/**
 * Fp_sqr - Square mod p
 *
 * Reference: pariinl.h:1770-1777
 * Computes a^2 mod p
 */
export function Fp_sqr(a: bigint, p: bigint): bigint {
  return (a * a) % p;
}

/**
 * Fp_inv - Modular inverse using extended Euclidean algorithm
 *
 * Reference: pariinl.h:1818-1823
 * Computes a^(-1) mod p
 * Throws if gcd(a, p) != 1
 */
export function Fp_inv(a: bigint, p: bigint): bigint {
  // Extended Euclidean algorithm
  // Find x such that a*x = gcd(a,p) (mod p)
  // If gcd(a,p) = 1, then a*x = 1 (mod p), so x = a^(-1)

  let [old_r, r] = [a, p];
  let [old_s, s] = [1n, 0n];

  while (r !== 0n) {
    const quotient = old_r / r;
    [old_r, r] = [r, old_r - quotient * r];
    [old_s, s] = [s, old_s - quotient * s];
  }

  // old_r is now gcd(a, p)
  if (old_r !== 1n && old_r !== -1n) {
    throw new Error(`Fp_inv: ${a} is not invertible mod ${p} (gcd = ${old_r})`);
  }

  // old_s is the inverse (may be negative)
  let result = old_s % p;
  if (result < 0n) result += p;
  return result;
}

/**
 * Fp_div - Division mod p
 *
 * Reference: pariinl.h:1832-1846
 * Computes a / b mod p = a * b^(-1) mod p
 */
export function Fp_div(a: bigint, b: bigint, p: bigint): bigint {
  return Fp_mul(a, Fp_inv(b, p), p);
}

/**
 * Fp_pow - Modular exponentiation using binary method (square-and-multiply)
 *
 * Reference: gen1.c:Fp_pow
 * Computes a^n mod p
 */
export function Fp_pow(a: bigint, n: bigint, p: bigint): bigint {
  if (n < 0n) {
    a = Fp_inv(a, p);
    n = -n;
  }

  if (n === 0n) return 1n;
  if (n === 1n) return Fp_red(a, p);

  let result = 1n;
  a = Fp_red(a, p);

  while (n > 0n) {
    if ((n & 1n) === 1n) {
      result = Fp_mul(result, a, p);
    }
    n >>= 1n;
    if (n > 0n) {
      a = Fp_sqr(a, p);
    }
  }

  return result;
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
 * Fp_issquare - Check if a is a quadratic residue mod p
 *
 * Reference: Uses Euler's criterion: a^((p-1)/2) = 1 if a is a QR
 * Or equivalently, kronecker(a, p) = 1
 *
 * Returns true if a is a square mod p
 */
export function Fp_issquare(a: bigint, p: bigint): boolean {
  if (a === 0n) return true;
  a = Fp_red(a, p);
  if (a === 0n) return true;

  // Use Euler's criterion: a is a square iff a^((p-1)/2) = 1 (mod p)
  const exp = (p - 1n) / 2n;
  return Fp_pow(a, exp, p) === 1n;
}

// ============================================================================
// Square root in Fp - Tonelli-Shanks algorithm
// Reference: arith1.c:762-858, 1189-1258
// ============================================================================

/**
 * Find a non-square mod p (a quadratic non-residue)
 *
 * Reference: arith1.c:1085-1096 (nonsquare_Fp)
 */
function nonsquare_Fp(p: bigint): bigint {
  // Special cases for small p mod 8
  if (mod4(p) === 3) return p - 1n; // -1 is a non-square if p = 3 mod 4
  if (mod8(p) === 5) return 2n; // 2 is a non-square if p = 5 mod 8

  // General case: search for a non-square
  for (let a = 2n; ; a++) {
    if (kronecker(a, p) < 0) return a;
  }
}

/**
 * Fp_sqrt - Square root mod p using Tonelli-Shanks algorithm
 *
 * Reference: arith1.c:801-858 (Fl_sqrt_i), 1189-1258 (Fp_sqrt_ii)
 *
 * Returns sqrt(a) mod p, or null if a is not a quadratic residue
 *
 * Algorithm (Tonelli-Shanks):
 * Write p-1 = 2^e * q where q is odd
 * Find a generator y of the 2-Sylow subgroup of (Z/pZ)*
 * Then compute the square root using the structure of the 2-Sylow
 */
export function Fp_sqrt(a: bigint, p: bigint): bigint | null {
  a = Fp_red(a, p);
  if (a === 0n) return 0n;

  // Check if a is a square
  if (!Fp_issquare(a, p)) return null;

  // Find e and q such that p-1 = 2^e * q, q odd
  const p1 = p - 1n;
  const e = vali(p1);

  // Handle p = 2 case
  if (e === 0) {
    if (p !== 2n) throw new Error('Fp_sqrt: p is not prime');
    return (a & 1n) === 0n ? 0n : 1n;
  }

  // Case e = 1: p = 3 mod 4, direct formula
  if (e === 1) {
    // sqrt(a) = a^((p+1)/4) mod p
    const exp = (p + 1n) / 4n;
    const v = Fp_pow(a, exp, p);
    // Verify
    if (Fp_sqr(v, p) !== a) return null;
    // Return the smaller of v and p-v
    const pv = p - v;
    return v > pv ? pv : v;
  }

  // Case e = 2: p = 5 mod 8, Atkin's formula
  if (e === 2) {
    // a2 = 2a
    const a2 = Fp_add(a, a, p);
    // v = (2a)^((p-5)/8) mod p
    const exp = (p - 5n) / 8n;
    let v = Fp_pow(a2, exp, p);
    // I = 2a * v^2 (should satisfy I^2 = -1)
    const I = Fp_mul(a2, Fp_sqr(v, p), p);
    // sqrt(a) = a * v * (I - 1)
    v = Fp_mul(a, Fp_mul(v, Fp_sub(I, 1n, p), p), p);
    // Verify
    if (Fp_sqr(v, p) !== a) return null;
    const pv = p - v;
    return v > pv ? pv : v;
  }

  // General case: Tonelli-Shanks
  // q = (p-1) / 2^e, q odd
  const q = p1 >> BigInt(e);

  // Find a generator y of the 2-Sylow subgroup
  // y = (non-square)^q is a primitive 2^e-th root of unity
  const ns = nonsquare_Fp(p);
  let y = Fp_pow(ns, q, p);

  // p1 = a^((q-1)/2)
  const p1_val = Fp_pow(a, (q - 1n) / 2n, p);

  // v = a * p1 = a^((q+1)/2)
  let v = Fp_mul(a, p1_val, p);

  // w = v * p1 = a^q
  let w = Fp_mul(v, p1_val, p);

  let currentE = e;

  // Main loop: reduce the order of w until w = 1
  while (w !== 1n) {
    // Find k such that w^(2^k) = 1
    let temp = Fp_sqr(w, p);
    let k = 1;
    while (temp !== 1n && k < currentE) {
      temp = Fp_sqr(temp, p);
      k++;
    }

    if (k === currentE) {
      // This means a is not a square (shouldn't happen if we checked earlier)
      return null;
    }

    // y1 = y^(2^(e-k-1))
    let y1 = y;
    for (let i = 1; i < currentE - k; i++) {
      y1 = Fp_sqr(y1, p);
    }

    // Update: y = y1^2, e = k, w = y * w, v = y1 * v
    y = Fp_sqr(y1, p);
    currentE = k;
    w = Fp_mul(y, w, p);
    v = Fp_mul(v, y1, p);
  }

  // Return the smaller of v and p-v
  const pv = p - v;
  return v > pv ? pv : v;
}

// ============================================================================
// Additional utility functions
// ============================================================================

/**
 * Fp_center - Return the representative in (-p/2, p/2]
 *
 * Reference: pariinl.h (Fp_center)
 */
export function Fp_center(a: bigint, p: bigint): bigint {
  const half = p / 2n;
  if (a > half) return a - p;
  return a;
}

/**
 * Fp_mulu - Multiply by a small unsigned integer
 */
export function Fp_mulu(a: bigint, b: number, p: bigint): bigint {
  return Fp_mul(a, BigInt(b), p);
}

/**
 * Fp_addmul - Compute x + y*z mod p
 *
 * Reference: pariinl.h:1751-1758
 */
export function Fp_addmul(x: bigint, y: bigint, z: bigint, p: bigint): bigint {
  if (y === 0n || z === 0n) return Fp_red(x, p);
  if (x === 0n) return Fp_mul(z, y, p);
  return (((x + y * z) % p) + p) % p;
}

/**
 * Fp_double - Compute 2*a mod p
 */
export function Fp_double(a: bigint, p: bigint): bigint {
  return Fp_add(a, a, p);
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
 * Extended GCD: returns [g, x, y] such that a*x + b*y = g = gcd(a,b)
 */
export function xgcd(a: bigint, b: bigint): [bigint, bigint, bigint] {
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
