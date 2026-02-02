/**
 * @module types
 * @description Core PARI types ported to TypeScript
 *
 * This file defines the fundamental types used in PARI/GP.
 * Reference: reference/pari/src/headers/parigen.h
 *
 * In PARI, GEN is a pointer to a long, representing the generic type.
 * We use TypeScript discriminated unions to represent this.
 */

/**
 * PARI type codes from parigen.h:175-198
 * Note: We only implement the subset needed for finite field arithmetic
 */
export enum PariType {
  t_INT = 1, // Integer
  t_REAL = 2, // Real number
  t_INTMOD = 3, // Integer mod n
  t_FRAC = 4, // Fraction
  t_FFELT = 5, // Finite field element
  t_COMPLEX = 6, // Complex number
  t_PADIC = 7, // p-adic number
  t_QUAD = 8, // Quadratic number
  t_POLMOD = 9, // Polynomial mod
  t_POL = 10, // Polynomial
  t_SER = 11, // Power series
  t_RFRAC = 13, // Rational function
  t_QFB = 15, // Quadratic binary form
  t_VEC = 17, // Row vector
  t_COL = 18, // Column vector
  t_MAT = 19, // Matrix
  t_LIST = 20, // List
  t_STR = 21, // String
  t_VECSMALL = 22, // Small integer vector
  t_CLOSURE = 23, // Closure
  t_ERROR = 24, // Error
  t_INFINITY = 25, // Infinity
}

/**
 * t_INT - PARI integer type
 *
 * In PARI, integers are stored as arrays of machine words.
 * In TypeScript, we use bigint for arbitrary precision.
 */
export interface PariInt {
  readonly type: PariType.t_INT;
  readonly value: bigint;
}

/**
 * t_REAL - PARI real number type
 *
 * In PARI, reals are multi-precision floating point.
 * For now, we use number as a placeholder.
 */
export interface PariReal {
  readonly type: PariType.t_REAL;
  readonly value: number;
  readonly precision: number;
}

/**
 * t_FFELT - Finite field element
 *
 * In PARI, this represents an element of a finite field Fq.
 * For prime fields Fp, the element is just an integer mod p.
 * For extension fields Fq where q = p^n, the element is a polynomial.
 *
 * Reference: reference/pari/src/basemath/FF.c
 */
export interface PariFfelt {
  readonly type: PariType.t_FFELT;
  /** The element value (integer for Fp, polynomial coefficients for Fq) */
  readonly value: bigint | readonly bigint[];
  /** The characteristic p */
  readonly p: bigint;
  /** The degree n (1 for prime fields) */
  readonly degree: number;
  /** For extension fields: the defining polynomial T */
  readonly definingPoly?: readonly bigint[];
}

/**
 * t_VEC - Row vector
 */
export interface PariVec {
  readonly type: PariType.t_VEC;
  readonly elements: readonly GEN[];
}

/**
 * t_COL - Column vector
 */
export interface PariCol {
  readonly type: PariType.t_COL;
  readonly elements: readonly GEN[];
}

/**
 * t_MAT - Matrix
 */
export interface PariMat {
  readonly type: PariType.t_MAT;
  readonly columns: readonly PariCol[];
}

/**
 * GEN - Generic PARI type
 *
 * This is the fundamental type in PARI. In C, GEN is `long*`.
 * In TypeScript, we use a discriminated union.
 */
export type GEN = PariInt | PariReal | PariFfelt | PariVec | PariCol | PariMat;

// ============================================================================
// Helper functions for creating GEN values
// ============================================================================

/**
 * Create a t_INT from a bigint
 */
export function mkInt(value: bigint): PariInt {
  return { type: PariType.t_INT, value };
}

/**
 * Create a t_INT from a number (JavaScript integer)
 */
export function stoi(value: number): PariInt {
  return { type: PariType.t_INT, value: BigInt(value) };
}

/**
 * Convert a t_INT to bigint
 */
export function itos(x: PariInt): bigint {
  return x.value;
}

/**
 * Create a t_FFELT for a prime field element
 */
export function mkFfeltFp(value: bigint, p: bigint): PariFfelt {
  // Reduce the value mod p
  let v = value % p;
  if (v < 0n) v += p;
  return {
    type: PariType.t_FFELT,
    value: v,
    p,
    degree: 1,
  };
}

/**
 * Create a t_VEC (row vector)
 */
export function mkvec(...elements: readonly GEN[]): PariVec {
  return { type: PariType.t_VEC, elements };
}

/**
 * Create a t_COL (column vector)
 */
export function mkcol(...elements: readonly GEN[]): PariCol {
  return { type: PariType.t_COL, elements };
}

/**
 * Create a t_MAT (matrix from columns)
 */
export function mkmat(...columns: readonly PariCol[]): PariMat {
  return { type: PariType.t_MAT, columns };
}

// ============================================================================
// Type guards
// ============================================================================

export function isInt(x: GEN): x is PariInt {
  return x.type === PariType.t_INT;
}

export function isReal(x: GEN): x is PariReal {
  return x.type === PariType.t_REAL;
}

export function isFfelt(x: GEN): x is PariFfelt {
  return x.type === PariType.t_FFELT;
}

export function isVec(x: GEN): x is PariVec {
  return x.type === PariType.t_VEC;
}

export function isCol(x: GEN): x is PariCol {
  return x.type === PariType.t_COL;
}

export function isMat(x: GEN): x is PariMat {
  return x.type === PariType.t_MAT;
}

// ============================================================================
// Integer utilities
// ============================================================================

/**
 * Sign of an integer: -1, 0, or 1
 */
export function signe(x: bigint): -1 | 0 | 1 {
  if (x < 0n) return -1;
  if (x > 0n) return 1;
  return 0;
}

/**
 * Absolute value
 */
export function absi(x: bigint): bigint {
  return x < 0n ? -x : x;
}

/**
 * Compare absolute values: returns -1, 0, or 1
 */
export function abscmpii(x: bigint, y: bigint): number {
  const ax = absi(x);
  const ay = absi(y);
  if (ax < ay) return -1;
  if (ax > ay) return 1;
  return 0;
}

/**
 * Check if integer equals 1
 */
export function equali1(x: bigint): boolean {
  return x === 1n;
}

/**
 * Check if integer is +1 or -1
 */
export function is_pm1(x: bigint): boolean {
  return x === 1n || x === -1n;
}

/**
 * gen_0, gen_1, gen_m1 - Common constants
 */
export const gen_0: PariInt = mkInt(0n);
export const gen_1: PariInt = mkInt(1n);
export const gen_m1: PariInt = mkInt(-1n);
export const gen_2: PariInt = mkInt(2n);

/**
 * Check if x is zero
 */
export function isZero(x: bigint): boolean {
  return x === 0n;
}

/**
 * mod2 - Return x mod 2 (as 0 or 1)
 */
export function mod2(x: bigint): 0 | 1 {
  return (x & 1n) === 0n ? 0 : 1;
}

/**
 * mod4 - Return x mod 4 (as 0, 1, 2, or 3)
 */
export function mod4(x: bigint): 0 | 1 | 2 | 3 {
  const r = Number(((x % 4n) + 4n) % 4n);
  return r as 0 | 1 | 2 | 3;
}

/**
 * mod8 - Return x mod 8
 */
export function mod8(x: bigint): number {
  return Number(((x % 8n) + 8n) % 8n);
}

/**
 * vali - 2-adic valuation: largest k such that 2^k divides x
 * Returns -1 for x = 0 (representing infinity)
 */
export function vali(x: bigint): number {
  if (x === 0n) return -1; // Represents infinity
  let k = 0;
  while ((x & 1n) === 0n) {
    x >>= 1n;
    k++;
  }
  return k;
}

/**
 * vals - 2-adic valuation for JavaScript numbers
 */
export function vals(x: number): number {
  if (x === 0) return -1;
  let k = 0;
  while ((x & 1) === 0) {
    x >>= 1;
    k++;
  }
  return k;
}
