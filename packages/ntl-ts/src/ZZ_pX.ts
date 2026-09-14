/** @see Deviation: Polynomial Integer Powers and Portable Native Products */
/**
 * NTL polynomials over Z_p (ZZ_pX).
 * @see Reference: ntl/src/ZZ_pX.cpp
 *
 * ZZ_pX represents polynomials with coefficients in Z_p.
 * The modulus p is determined by the current ZZ_p modulus.
 */

import type { ZZ } from './ZZ.js';
import type { ZZ_p } from './ZZ_p.js';
import { RandomStream } from './ZZ.js';
import { _ZZ_pX_euclidean_kernels } from './ZZ_pX1.js';
import { FFTRoundUp } from './FFT_impl.js';
import {
  FFTPrimeContext,
  UseFFTPrime,
  CalcMaxRoot,
  FFTFwd_trunc,
  FFTRev1_trunc,
  type FFTPrimeInfo,
} from './FFT.js';
import { SSRatio } from './ZZX1.js';

/**
 * Polynomials over Z_p.
 * TypeScript port of NTL's ZZ_pX class.
 *
 * IMPORTANT: ZZ_p.init(p) must be called before using ZZ_pX.
 */
export class ZZ_pX {
  private _rep: ZZ_p[];

  /**
   * Creates a new polynomial.
   * @param coeffs - Coefficients (constant term first)
   */
  constructor(coeffs?: ZZ_p[]) {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_pX');
  }

  // ============================================
  // Basic Properties
  // ============================================

  /**
   * Returns the zero polynomial.
   */
  static zero(): ZZ_pX {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_pX.zero');
  }

  /**
   * Returns the one polynomial.
   */
  static one(): ZZ_pX {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_pX.one');
  }

  /**
   * Returns the polynomial X.
   */
  static X(): ZZ_pX {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_pX.X');
  }

  /**
   * Returns the degree (-1 for zero polynomial).
   * @returns The degree
   */
  deg(): number {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_pX.deg');
  }

  /**
   * Normalizes the polynomial (removes leading zeros).
   */
  normalize(): void {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_pX.normalize');
  }

  /**
   * Gets coefficient at index i.
   * @param i - Index
   * @returns The coefficient
   */
  coeff(i: number): ZZ_p {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_pX.coeff');
  }

  /**
   * Sets coefficient at index i.
   * @param i - Index
   * @param a - Value
   */
  SetCoeff(i: number, a: ZZ_p): void {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_pX.SetCoeff');
  }

  /**
   * Returns the leading coefficient.
   * @returns Leading coefficient
   */
  LeadCoeff(): ZZ_p {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_pX.LeadCoeff');
  }

  /**
   * Returns the constant term.
   * @returns Constant term
   */
  ConstTerm(): ZZ_p {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_pX.ConstTerm');
  }

  // ============================================
  // Arithmetic Operations
  // ============================================

  /**
   * Adds two polynomials.
   * @param other - Polynomial to add
   * @returns The sum
   */
  add(other: ZZ_pX): ZZ_pX {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_pX.add');
  }

  /**
   * Subtracts a polynomial.
   * @param other - Polynomial to subtract
   * @returns The difference
   */
  sub(other: ZZ_pX): ZZ_pX {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_pX.sub');
  }

  /**
   * Multiplies two polynomials.
   * @param other - Polynomial to multiply
   * @returns The product
   */
  mul(other: ZZ_pX): ZZ_pX {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_pX.mul');
  }

  /**
   * Multiplies by a scalar.
   * @param c - Scalar
   * @returns The product
   */
  mulScalar(c: ZZ_p): ZZ_pX {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_pX.mulScalar');
  }

  /**
   * Squares the polynomial.
   * @returns The square
   */
  sqr(): ZZ_pX {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_pX.sqr');
  }

  /**
   * Negates the polynomial.
   * @returns The negation
   */
  negate(): ZZ_pX {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_pX.negate');
  }

  // ============================================
  // Division Operations
  // ============================================

  /**
   * Polynomial division with remainder.
   * @param b - Divisor
   * @returns [quotient, remainder]
   */
  DivRem(b: ZZ_pX): [ZZ_pX, ZZ_pX] {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_pX.DivRem');
  }

  /**
   * Polynomial division (quotient only).
   * @param b - Divisor
   * @returns The quotient
   */
  div(b: ZZ_pX): ZZ_pX {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_pX.div');
  }

  /**
   * Polynomial remainder.
   * @param b - Divisor
   * @returns The remainder
   */
  rem(b: ZZ_pX): ZZ_pX {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_pX.rem');
  }

  // ============================================
  // GCD and Related
  // ============================================

  /**
   * Computes GCD of two polynomials.
   * @param a - First polynomial
   * @param b - Second polynomial
   * @returns The GCD (monic)
   */
  static GCD(a: ZZ_pX, b: ZZ_pX): ZZ_pX {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_pX.GCD');
  }

  /**
   * Extended GCD.
   * Computes d, s, t such that d = GCD(a, b) = s*a + t*b
   * @param a - First polynomial
   * @param b - Second polynomial
   * @returns [d, s, t]
   */
  static XGCD(a: ZZ_pX, b: ZZ_pX): [ZZ_pX, ZZ_pX, ZZ_pX] {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_pX.XGCD');
  }

  // ============================================
  // Modular Arithmetic
  // ============================================

  /**
   * Computes (a * b) mod f.
   * @param a - First polynomial
   * @param b - Second polynomial
   * @param f - Modulus polynomial
   * @returns (a * b) mod f
   */
  static MulMod(a: ZZ_pX, b: ZZ_pX, f: ZZ_pX): ZZ_pX {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_pX.MulMod');
  }

  /**
   * Computes a^2 mod f.
   * @param a - Polynomial
   * @param f - Modulus polynomial
   * @returns a^2 mod f
   */
  static SqrMod(a: ZZ_pX, f: ZZ_pX): ZZ_pX {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_pX.SqrMod');
  }

  /**
   * Computes modular inverse.
   * @param a - Polynomial to invert
   * @param f - Modulus polynomial
   * @returns a^(-1) mod f
   */
  static InvMod(a: ZZ_pX, f: ZZ_pX): ZZ_pX {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_pX.InvMod');
  }

  /**
   * Computes a^e mod f.
   * @param a - Base polynomial
   * @param e - Exponent
   * @param f - Modulus polynomial
   * @returns a^e mod f
   */
  static PowerMod(a: ZZ_pX, e: ZZ, f: ZZ_pX): ZZ_pX {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_pX.PowerMod');
  }

  // ============================================
  // Evaluation and Interpolation
  // ============================================

  /**
   * Evaluates polynomial at a point.
   * @param a - Point to evaluate at
   * @returns f(a)
   */
  eval(a: ZZ_p): ZZ_p {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_pX.eval');
  }

  /**
   * Interpolates polynomial through points.
   * @param x - X coordinates
   * @param y - Y coordinates
   * @returns Interpolating polynomial
   */
  static interpolate(x: ZZ_p[], y: ZZ_p[]): ZZ_pX {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_pX.interpolate');
  }

  // ============================================
  // Factoring
  // ============================================

  /**
   * Factors a polynomial.
   * @returns Array of [factor, multiplicity] pairs
   */
  factor(): Array<[ZZ_pX, number]> {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_pX.factor');
  }

  /**
   * Computes squarefree factorization.
   * @returns Array of [factor, multiplicity] pairs
   */
  SquareFreeDecomp(): Array<[ZZ_pX, number]> {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_pX.SquareFreeDecomp');
  }

  /**
   * Distinct degree factorization.
   * @returns Array of [product of degree-d irreducibles, d] pairs
   */
  DistinctDegFactor(): Array<[ZZ_pX, number]> {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_pX.DistinctDegFactor');
  }

  /**
   * Equal degree factorization.
   * Factors polynomial into irreducibles of specified degree.
   * @param d - Degree of factors
   * @returns Array of irreducible factors
   */
  EqualDegFactor(d: number): ZZ_pX[] {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_pX.EqualDegFactor');
  }

  /**
   * Berlekamp factorization.
   * @returns Array of irreducible factors
   */
  BerlekampFactor(): ZZ_pX[] {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_pX.BerlekampFactor');
  }

  /**
   * Cantor-Zassenhaus factorization.
   * @returns Array of irreducible factors
   */
  CanZassFactor(): ZZ_pX[] {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_pX.CanZassFactor');
  }

  // ============================================
  // Roots
  // ============================================

  /**
   * Finds roots of polynomial.
   * @returns Array of roots in Z_p
   */
  FindRoots(): ZZ_p[] {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_pX.FindRoots');
  }

  // ============================================
  // Special Polynomials
  // ============================================

  /**
   * Computes derivative.
   * @returns The derivative
   */
  diff(): ZZ_pX {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_pX.diff');
  }

  /**
   * Makes polynomial monic (leading coeff = 1).
   * @returns Monic polynomial
   */
  MakeMonic(): ZZ_pX {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_pX.MakeMonic');
  }

  /**
   * Reverses coefficients.
   * @returns Reversed polynomial
   */
  reverse(): ZZ_pX {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_pX.reverse');
  }

  /**
   * Left shift (multiply by X^n).
   * @param n - Shift amount
   * @returns Shifted polynomial
   */
  LeftShift(n: number): ZZ_pX {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_pX.LeftShift');
  }

  /**
   * Right shift (divide by X^n).
   * @param n - Shift amount
   * @returns Shifted polynomial
   */
  RightShift(n: number): ZZ_pX {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_pX.RightShift');
  }

  // ============================================
  // Comparison and Utility
  // ============================================

  /**
   * Checks if polynomial is zero.
   * @returns True if zero
   */
  IsZero(): boolean {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_pX.IsZero');
  }

  /**
   * Checks if polynomial is one.
   * @returns True if one
   */
  IsOne(): boolean {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_pX.IsOne');
  }

  /**
   * Checks if polynomial is X.
   * @returns True if X
   */
  IsX(): boolean {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_pX.IsX');
  }

  /**
   * Checks equality.
   * @param other - Polynomial to compare
   * @returns True if equal
   */
  equals(other: ZZ_pX): boolean {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_pX.equals');
  }

  /**
   * Converts to string.
   * @returns String representation
   */
  toString(): string {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_pX.toString');
  }

  // ============================================
  // Random
  // ============================================

  /**
   * Generates a random polynomial of degree < n, i.e. with n random
   * coefficients (NTL: `random(ZZ_pX& x, long n)`, ntl/doc/ZZ_pX.txt:391).
   * @param n - Number of random coefficients; the result has degree < n
   * @returns Random polynomial of degree < n
   */
  static random(n: number): ZZ_pX {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_pX.random');
  }

  /**
   * Builds from coefficient array.
   * @param coeffs - Coefficients
   * @returns Polynomial
   */
  static BuildFromVec(coeffs: ZZ_p[]): ZZ_pX {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ_pX.BuildFromVec');
  }
}

// ============================================
// Standalone Functions (NTL-style API)
// ============================================

/**
 * Returns degree of polynomial.
 */
export function deg(f: ZZ_pX): number {
  return f.deg();
}

/**
 * Returns coefficient at index i.
 */
export function coeff(f: ZZ_pX, i: number): ZZ_p {
  return f.coeff(i);
}

/**
 * Returns leading coefficient.
 */
export function LeadCoeff(f: ZZ_pX): ZZ_p {
  return f.LeadCoeff();
}

/**
 * Returns constant term.
 */
export function ConstTerm(f: ZZ_pX): ZZ_p {
  return f.ConstTerm();
}

/**
 * Computes GCD.
 */
export function GCD(a: ZZ_pX, b: ZZ_pX): ZZ_pX {
  return ZZ_pX.GCD(a, b);
}

/**
 * Extended GCD.
 */
export function XGCD(a: ZZ_pX, b: ZZ_pX): [ZZ_pX, ZZ_pX, ZZ_pX] {
  return ZZ_pX.XGCD(a, b);
}

/**
 * Checks if zero.
 */
export function IsZero(f: ZZ_pX): boolean {
  return f.IsZero();
}

/**
 * Checks if one.
 */
export function IsOne(f: ZZ_pX): boolean {
  return f.IsOne();
}

/**
 * Evaluates polynomial.
 */
export function eval_(f: ZZ_pX, a: ZZ_p): ZZ_p {
  return f.eval(a);
}

/**
 * Computes derivative.
 */
export function diff(f: ZZ_pX): ZZ_pX {
  return f.diff();
}

/**
 * Computes modular product.
 * @see Deviation: NTL arbitrary-modulus polynomial quotient adapters
 */
export function MulMod(a: ZZ_pX, b: ZZ_pX, f: ZZ_pX): ZZ_pX;
export function MulMod(
  a: readonly bigint[],
  b: readonly bigint[] | ZZ_pXMultiplier,
  f: ZZ_pXModulus
): bigint[];
export function MulMod(
  a: ZZ_pX | readonly bigint[],
  b: ZZ_pX | readonly bigint[] | ZZ_pXMultiplier,
  f: ZZ_pX | ZZ_pXModulus
): ZZ_pX | bigint[] {
  if (f instanceof ZZ_pXModulus)
    return quotientMulMod(a as readonly bigint[], b as readonly bigint[] | ZZ_pXMultiplier, f);
  return ZZ_pX.MulMod(a as ZZ_pX, b as ZZ_pX, f);
}

/**
 * Computes modular power.
 */
export function PowerMod(a: ZZ_pX, e: ZZ, f: ZZ_pX): ZZ_pX {
  return ZZ_pX.PowerMod(a, e, f);
}

/**
 * Computes modular inverse.
 */
export function InvMod(a: ZZ_pX, f: ZZ_pX): ZZ_pX {
  return ZZ_pX.InvMod(a, f);
}

/** Dense-array adapter for NTL ZZ_pX.cpp:power. */
export function power(a: readonly bigint[], e: bigint, p: bigint): bigint[] {
  if (p <= 1n) throw new RangeError('modulus must exceed one');
  if (e < 0n) throw new Error('power: negative exponent');
  if (e >= 1n << 63n) throw new RangeError('exponent must fit a signed word');
  const A = a.map((c) => ((c % p) + p) % p);
  while (A.length && A[A.length - 1] === 0n) A.pop();
  if (e === 0n) return [1n];
  if (!A.length || (A.length === 1 && A[0] === 1n)) return A;
  if (A.length === 1) {
    let n = e,
      result = 1n,
      base = A[0]!;
    while (n) {
      if (n & 1n) result = (result * base) % p;
      n >>= 1n;
      if (n) base = (base * base) % p;
    }
    return result ? [result] : [];
  }
  if (BigInt(A.length - 1) > ((1n << 63n) - 2n) / e) throw new Error('overflow in power');
  let result = [1n];
  for (const bit of e.toString(2)) {
    result = denseProduct(result, result, p);
    if (bit === '1') result = denseProduct(result, A, p);
  }
  return result;
}

/** Portable dense product boundary used by ZZ_pX power. */
function denseProduct(a: readonly bigint[], b: readonly bigint[], p: bigint): bigint[] {
  if (!a.length || !b.length) return [];
  const length = a.length + b.length - 1;
  let result: bigint[];
  if (Math.min(a.length, b.length) <= 20) {
    result = Array<bigint>(length).fill(0n);
    if (a === b) {
      for (let i = 0; i < a.length; i++) {
        result[2 * i] = result[2 * i]! + a[i]! * a[i]!;
        for (let j = i + 1; j < a.length; j++) result[i + j] = result[i + j]! + 2n * a[i]! * a[j]!;
      }
    } else {
      for (let i = 0; i < a.length; i++)
        for (let j = 0; j < b.length; j++) result[i + j] = result[i + j]! + a[i]! * b[j]!;
    }
  } else {
    // Exact packed BigInt convolution replaces the native FFT integer buffers.
    const width =
      2 * p.toString(2).length + BigInt(Math.min(a.length, b.length)).toString(2).length;
    if (length * width > 1 << 19) {
      // ZZX.cpp KarMul: split the longer operand in half, with a two-product
      // branch for an unbalanced pair. Bound each packed temporary, not the
      // output degree: valid individual coefficients can fit when packing fails.
      const x = a.length >= b.length ? a : b,
        y = a.length >= b.length ? b : a;
      const half = Math.ceil(x.length / 2),
        lo = x.slice(0, half),
        hi = x.slice(half);
      result = Array<bigint>(length).fill(0n);
      const accumulate = (v: readonly bigint[], offset: number, sign: bigint = 1n) => {
        for (let i = 0; i < v.length; i++) result[i + offset] = result[i + offset]! + sign * v[i]!;
      };
      if (y.length <= half) {
        accumulate(denseProduct(lo, y, p), 0);
        accumulate(denseProduct(hi, y, p), half);
      } else {
        const yl = y.slice(0, half),
          yh = y.slice(half);
        const fold = (l: readonly bigint[], h: readonly bigint[]) =>
          l.map((c, i) => (c + (h[i] ?? 0n)) % p);
        const low = denseProduct(lo, yl, p),
          high = denseProduct(hi, yh, p);
        const middle = denseProduct(fold(lo, hi), fold(yl, yh), p);
        accumulate(low, 0);
        accumulate(high, 2 * half);
        accumulate(middle, half);
        accumulate(low, half, -1n);
        accumulate(high, half, -1n);
      }
      result = result.map((c) => ((c % p) + p) % p);
      while (result.length && result[result.length - 1] === 0n) result.pop();
      return result;
    }
    const pack = (v: readonly bigint[]) =>
      BigInt(
        '0b' +
          v
            .map((c) => c.toString(2).padStart(width, '0'))
            .reverse()
            .join('')
      );
    const A = pack(a),
      B = a === b ? A : pack(b),
      bits = (A * B).toString(2).padStart(length * width, '0');
    result = [];
    for (let i = 0; i < length; i++) {
      const end = bits.length - i * width;
      result.push(BigInt('0b' + bits.slice(end - width, end)));
    }
  }
  result = result.map((c) => c % p);
  while (result.length && result[result.length - 1] === 0n) result.pop();
  return result;
}
/** @internal Shared native product boundary for extension polynomial Kronecker substitution. */
export const _ZZ_pX_power_kernels = { multiply: denseProduct };

/** Explicit native FFT cache and random stream shared by modular products.
 * @see Deviation: NTL stateful modular polynomial products
 */
export interface PolynomialProductState {
  context: FFTPrimeContext;
  stream: RandomStream;
}

/** ZZ_pX.cpp multiplication/squaring dispatch and ZZ_p.cpp lazy FFT installation. */
function installProductFFT(
  a: bigint[],
  b: bigint[],
  p: bigint,
  state: PolynomialProductState
): void {
  if (!a.length || !b.length) return;
  const square = a === b,
    size = Math.min(a.length, b.length);
  // Every native PlainMul/PlainSqr crossover lies below the Karatsuba limit.
  if (size < (square ? 80 : 200)) return;
  const mbits = p.toString(2).length,
    limbs = Math.ceil(mbits / 64);
  const ratio = SSRatio(a.length - 1, mbits, b.length - 1, mbits);
  const ss = square
    ? (limbs >= 53 && ratio < 1.2) ||
      (limbs >= 106 && ratio < 1.3) ||
      (limbs >= 212 && ratio < 1.75)
    : (limbs >= 106 && ratio < 1.5) || (limbs >= 212 && ratio < 1.75);
  if (ss) return;
  const info = quotientFFTInfo(p, state);
  if (Math.ceil(Math.log2(a.length + b.length - 1)) > info.MaxRoot)
    throw new Error('Polynomial too big for FFT');
}

/** Exact modular product, optionally preserving the native FFT initialization stream.
 * @see Deviation: NTL stateful modular polynomial products
 */
export function mul(
  a: readonly bigint[],
  b: readonly bigint[],
  p: bigint,
  state?: PolynomialProductState
): bigint[] {
  if (p <= 1n) throw new RangeError('modulus must exceed one');
  const normalize = (v: readonly bigint[]) => {
    const out = v.map((c) => ((c % p) + p) % p);
    while (out.length && out[out.length - 1] === 0n) out.pop();
    return out;
  };
  const A = normalize(a),
    B = a === b ? A : normalize(b);
  if (state) installProductFFT(A, B, p, state);
  return denseProduct(A, B, p);
}

/** Array adapter for NTL ZZ_pX1.cpp eval. @see Deviation: Polynomial Evaluation and Composition */
export function evaluate(a: readonly bigint[], x: bigint, p: bigint): bigint {
  if (p <= 1n) throw new RangeError('modulus must exceed one');
  const point = ((x % p) + p) % p;
  let acc = 0n;
  for (let i = a.length - 1; i >= 0; i--) acc = (((acc * point + a[i]!) % p) + p) % p;
  return acc;
}

/** Explicit arbitrary-modulus polynomial quotient and retained native caches.
 * @see Deviation: NTL arbitrary-modulus polynomial quotient adapters
 */
export class ZZ_pXModulus {
  f: bigint[] = [];
  n = -1;
  UseFFT = false;
  readonly p: bigint;
  readonly state?: PolynomialProductState;
  /** @internal Exact coefficient arithmetic with native scalar inverse errors. */
  readonly arithmetic: ReturnType<typeof _ZZ_pX_euclidean_kernels>;
  constructor(f: readonly bigint[] | null, p: bigint, state?: PolynomialProductState) {
    this.arithmetic = quotientArithmetic(p);
    this.p = p;
    this.state = state;
    if (f !== null) build(this, f);
  }
  /** Return an independent copy of the current public modulus polynomial. */
  val(): bigint[] {
    return this.f.slice();
  }
}
/** Cached arbitrary-modulus multiplier.
 * @see Deviation: NTL arbitrary-modulus polynomial quotient adapters
 */
export class ZZ_pXMultiplier {
  b: bigint[] = [];
  UseFFT = false;
  constructor();
  constructor(b: readonly bigint[], F: ZZ_pXModulus);
  constructor(b?: readonly bigint[], F?: ZZ_pXModulus) {
    if (b !== undefined) build(this, b, F!);
  }
  /** Return an independent copy of the current public multiplier polynomial. */
  val(): bigint[] {
    return this.b.slice();
  }
}
function quotientArithmetic(p: bigint) {
  if (p <= 1n) throw new Error('ZZ_pContext: p must be > 1');
  const base = _ZZ_pX_euclidean_kernels(p);
  const call = <T>(f: () => T): T => {
    try {
      return f();
    } catch (e) {
      if (e instanceof Error && e.message === "ZZ_pX InvMod: can't compute multiplicative inverse")
        throw new Error('ZZ_p: division by non-invertible element');
      throw e;
    }
  };
  return {
    ...base,
    inverse: (a: bigint) => call(() => base.inverse(a)),
    divrem: (a: bigint[], b: bigint[]) => call(() => base.divrem(a, b)),
  };
}
interface QuotientFFTInfo {
  p: bigint;
  primes: FFTPrimeInfo[];
  MaxRoot: number;
  reconstruct: (rows: readonly (readonly bigint[])[], column: number) => bigint;
}
// Keep only the most recent modulus per global prime context. Live quotient
// objects retain their own information, like restored native ZZ_pContexts.
const quotientInstalled = new WeakMap<FFTPrimeContext, QuotientFFTInfo>();
const quotientDoubleView = new DataView(new ArrayBuffer(8));
function quotientFusedEstimate(a: number, b: number, c: number): number {
  if (c === 0) return a * b;
  if (a === 0) return c;
  const parts = (x: number): [bigint, number] => {
    quotientDoubleView.setFloat64(0, x, false);
    const bits = quotientDoubleView.getBigUint64(0, false);
    return [(bits & ((1n << 52n) - 1n)) | (1n << 52n), Number((bits >> 52n) & 2047n) - 1075];
  };
  const [am, ae] = parts(a),
    [bm, be] = parts(b),
    [cm, ce] = parts(c),
    exponent = Math.min(ae + be, ce);
  const exact = ((am * bm) << BigInt(ae + be - exponent)) + (cm << BigInt(ce - exponent));
  // CRT operands are positive normals or zero, and the sum is at most 800.
  return Number(exact) * 2 ** exponent;
}
interface QuotientCRTNode {
  start: number;
  end: number;
  product: bigint;
  left?: QuotientCRTNode;
  right?: QuotientCRTNode;
  coefficients?: bigint[];
  inverses?: bigint[];
}
function quotientCRTTree(primes: FFTPrimeInfo[], p: bigint) {
  let levels = 0;
  while (Math.floor(primes.length / 2 ** levels) >= 16) levels++;
  const make = (start: number, end: number, depth: number): QuotientCRTNode => {
    if (depth === levels - 1) {
      let product = 1n;
      for (let i = start; i < end; i++) product *= primes[i]!.q;
      return {
        start,
        end,
        product,
        coefficients: primes.slice(start, end).map((info) => product / info.q),
      };
    }
    const mid = start + Math.floor((end - start) / 2),
      left = make(start, mid, depth + 1),
      right = make(mid, end, depth + 1);
    return { start, end, product: left.product * right.product, left, right };
  };
  const root = make(0, primes.length, 0);
  const prepare = (node: QuotientCRTNode, outside: bigint): void => {
    if (node.left && node.right) {
      prepare(node.left, ((outside % node.left.product) * node.right.product) % node.left.product);
      prepare(
        node.right,
        ((outside % node.right.product) * node.left.product) % node.right.product
      );
    } else {
      node.inverses = node.coefficients!.map((cofactor, j) => {
        const q = primes[node.start + j]!.q;
        return _ZZ_pX_euclidean_kernels(q).inverse(((cofactor % q) * (outside % q)) % q);
      });
    }
  };
  prepare(root, 1n);
  return (rows: readonly (readonly bigint[])[], column: number): bigint => {
    const evaluate = (node: QuotientCRTNode): bigint => {
      if (node.left && node.right)
        return node.left.product * evaluate(node.right) + node.right.product * evaluate(node.left);
      let sum = 0n;
      for (let i = node.start; i < node.end; i++) {
        const q = primes[i]!.q,
          value = (((rows[i]![column]! % q) + q) * node.inverses![i - node.start]!) % q;
        sum += value * node.coefficients![i - node.start]!;
      }
      return sum;
    };
    let value = evaluate(root) % root.product;
    if (value > root.product - value) value -= root.product;
    return ((value % p) + p) % p;
  };
}
// ZZ_p.cpp DoInstall, including the native small/fast CRT construction threshold.
function quotientFFTInfo(p: bigint, state: PolynomialProductState): QuotientFFTInfo {
  const previous = quotientInstalled.get(state.context);
  if (previous?.p === p) return previous;
  const bound = (p * p) << 29n,
    primes: FFTPrimeInfo[] = [];
  let product = 1n;
  while (product <= bound) {
    UseFFTPrime(primes.length, state.context, state.stream);
    const info = state.context.get(primes.length);
    primes.push(info);
    product *= info.q;
  }
  if (8 * primes.length * (primes.length + 48) > 2 ** 52) throw new Error('modulus too big');
  let reconstruct: QuotientFFTInfo['reconstruct'];
  if (primes.length > 800) reconstruct = quotientCRTTree(primes, p);
  else {
    const cofactors = primes.map((info) => product / info.q),
      inverses = primes.map((info, i) =>
        _ZZ_pX_euclidean_kernels(info.q).inverse(cofactors[i]! % info.q)
      ),
      weights = cofactors.map((x) => x % p),
      productModP = product % p;
    reconstruct = (rows, column) => {
      let sum = 0n,
        estimate = 0;
      for (let i = 0; i < primes.length; i++) {
        const info = primes[i]!,
          value = (((rows[i]![column]! % info.q) + info.q) * inverses[i]!) % info.q;
        estimate = quotientFusedEstimate(Number(value), info.qrecip, estimate);
        sum += value * weights[i]!;
      }
      const value = (sum - BigInt(Math.floor(estimate + 0.5)) * productModP) % p;
      return value < 0n ? value + p : value;
    };
  }
  const info = { p, primes, MaxRoot: CalcMaxRoot(primes[primes.length - 1]!.q), reconstruct };
  quotientInstalled.set(state.context, info);
  return info;
}
interface QuotientCache {
  fallback?: PolynomialProductState;
  info?: QuotientFFTInfo;
  k: number;
  l: number;
  fk: number;
  hk: number;
  f: bigint[];
  h: bigint[];
}
const quotientCaches = new WeakMap<ZZ_pXModulus, QuotientCache>();
function quotientCache(F: ZZ_pXModulus): QuotientCache {
  let cache = quotientCaches.get(F);
  if (!cache) {
    cache = { k: 0, l: 0, fk: -1, hk: -1, f: [], h: [] };
    quotientCaches.set(F, cache);
  }
  return cache;
}
interface QuotientMultiplierCache {
  primes: number;
  b1k: number;
  b2k: number;
  b1: bigint[];
  b2: bigint[];
}
const quotientMultiplierCaches = new WeakMap<ZZ_pXMultiplier, QuotientMultiplierCache>();
function quotientMultiplierCache(B: ZZ_pXMultiplier): QuotientMultiplierCache {
  let cache = quotientMultiplierCaches.get(B);
  if (!cache) {
    cache = { primes: 0, b1k: -1, b2k: -1, b1: [], b2: [] };
    quotientMultiplierCaches.set(B, cache);
  }
  return cache;
}
function quotientContext(F: ZZ_pXModulus): QuotientFFTInfo {
  const cache = quotientCache(F);
  return (cache.info ??= quotientFFTInfo(
    F.p,
    F.state ??
      (cache.fallback ??= {
        context: new FFTPrimeContext(),
        stream: new RandomStream(new Uint8Array(32)),
      })
  ));
}
function quotientFFTSize(length: number): number {
  return length <= 1 ? 0 : BigInt(length - 1).toString(2).length;
}
function quotientCheckSize(size: number, F: ZZ_pXModulus): void {
  if (size > quotientContext(F).MaxRoot) throw new Error('Polynomial too big for FFT');
}
function quotientFold(a: readonly bigint[], size: number, F: ZZ_pXModulus): bigint[] {
  const n = 2 ** size;
  if (a.length <= n) return F.arithmetic.norm(a);
  const out = Array<bigint>(n).fill(0n);
  for (let i = 0; i < a.length; i++) out[i % n] += a[i]!;
  return F.arithmetic.norm(out);
}
/** Native arbitrary-modulus CRT, including the 800-prime algorithm boundary.
 * @see Deviation: NTL arbitrary-modulus CRT reconstruction
 */
export function FromModularRep(
  residues: readonly (readonly bigint[])[],
  F: ZZ_pXModulus
): bigint[] {
  const info = quotientContext(F);
  if (residues.length !== info.primes.length)
    throw new RangeError('FromModularRep: incorrect number of prime rows');
  const columns = residues[0]!.length;
  if (residues.some((row) => row.length !== columns))
    throw new RangeError('FromModularRep: inconsistent coefficient counts');
  return Array.from({ length: columns }, (_, j) => info.reconstruct(residues, j));
}
function quotientTruncatedProduct(
  a: bigint[],
  b: bigint[],
  size: number,
  requested: number,
  F: ZZ_pXModulus
): bigint[] {
  const len = FFTRoundUp(requested, size);
  if (len === 2 ** size || a.length + b.length - 1 <= len)
    return quotientFold(F.arithmetic.mul(a, F.arithmetic.norm(b)), size, F);
  const residues = quotientContext(F).primes.map((info) => {
    const forward = (input: bigint[], retained: boolean) => {
      // Cached HRep/B1 coefficients retain their original modulus. Converting
      // them under the current p changes a partially inverted FFT prefix.
      const values = retained ? input.slice() : quotientFold(input, size, F),
        xn = FFTRoundUp(values.length, size);
      while (values.length < xn) values.push(0n);
      return FFTFwd_trunc(values, size, info, len, xn);
    };
    const A = forward(a, false),
      B = forward(b, true);
    return FFTRev1_trunc(
      A.map((x, i) => (x * B[i]!) % info.q),
      size,
      info,
      len
    );
  });
  return F.arithmetic.norm(FromModularRep(residues, F));
}
function quotientInverse(a: bigint[], m: number, F: ZZ_pXModulus): bigint[] {
  const k = F.arithmetic,
    scale = k.inverse(a[0]!);
  const plain = (length: number) => {
    const out = Array<bigint>(length).fill(0n);
    out[0] = scale;
    for (let j = 1; j < length; j++) {
      let sum = 0n;
      for (let i = Math.max(0, j - a.length + 1); i < j; i++) sum += out[i]! * a[j - i]!;
      out[j] = k.mod(-sum * scale);
    }
    return k.norm(out);
  };
  if (a.length === 1) return [scale];
  if (m <= 45) return plain(m);
  let out = plain(32),
    length = 32;
  while (length < m) {
    const next = Math.min(2 * length, m);
    quotientCheckSize(quotientFFTSize(2 * length), F);
    const high = k.mul(a.slice(0, next), out).slice(length, next),
      correction = k.mul(high, out).slice(0, next - length);
    out = k.norm(
      Array.from({ length: next }, (_, i) =>
        i < length ? (out[i] ?? 0n) : -(correction[i - length] ?? 0n)
      )
    );
    length = next;
  }
  return out;
}
/** Native quotient/multiplier construction with independent retained FFT caches.
 * @see Deviation: NTL arbitrary-modulus polynomial quotient adapters
 */
export function build(F: ZZ_pXModulus, f: readonly bigint[]): void;
export function build(B: ZZ_pXMultiplier, b: readonly bigint[], F: ZZ_pXModulus): void;
export function build(
  out: ZZ_pXModulus | ZZ_pXMultiplier,
  a: readonly bigint[],
  context?: ZZ_pXModulus
): void {
  if (out instanceof ZZ_pXModulus) {
    const F = out,
      k = F.arithmetic,
      cache = quotientCache(F);
    F.f = k.norm(a);
    F.n = F.f.length - 1;
    if (F.n <= 0) throw new Error('build: deg(f) must be at least 1');
    if (F.n <= 21) {
      F.UseFFT = false;
      return;
    }
    F.UseFFT = true;
    cache.k = quotientFFTSize(F.n);
    cache.l = quotientFFTSize(2 * F.n - 3);
    quotientCheckSize(cache.k, F);
    cache.f = quotientFold(F.f, cache.k, F);
    cache.fk = cache.k;
    const reciprocal = quotientInverse(k.norm(F.f.slice().reverse()), F.n - 1, F);
    quotientCheckSize(cache.l, F);
    cache.h = Array.from({ length: F.n - 1 }, (_, i) => reciprocal[F.n - 2 - i] ?? 0n);
    cache.hk = cache.l;
    return;
  }
  const B = out,
    F = context!;
  if (F.n < 0) throw new Error('build ZZ_pXMultiplier: uninitialized modulus');
  B.b = F.arithmetic.norm(a);
  if (B.b.length - 1 >= F.n) throw new Error('build ZZ_pXMultiplier: deg(b) >= deg(f)');
  B.UseFFT = F.UseFFT && B.b.length - 1 > 20;
  if (!B.UseFFT) return;
  const cache = quotientCache(F),
    multiplier = quotientMultiplierCache(B);
  quotientCheckSize(cache.l, F);
  // FFTRep::DoSetSize retains NumPrimes across plain and failed rebuilds.
  // B2 is resized before the inverse transform is read or B1 is overwritten.
  const primes = quotientContext(F).primes.length;
  if (multiplier.primes !== 0 && multiplier.primes !== primes)
    throw new Error('FFTRep: inconsistent use');
  multiplier.primes = primes;
  multiplier.b2 = quotientFold(B.b, cache.k, F);
  multiplier.b2k = cache.k;
  if (cache.l !== cache.hk) throw new Error('FFT rep mismatch');
  const product = quotientTruncatedProduct(B.b, cache.h, cache.l, 2 * F.n - 2, F);
  multiplier.b1 = product.slice(F.n - 1, 2 * F.n - 2);
  multiplier.b1k = cache.l;
}
function quotientFFTRemainder(a: bigint[], F: ZZ_pXModulus): bigint[] {
  const cache = quotientCache(F),
    k = F.arithmetic,
    n = F.n;
  quotientCheckSize(cache.l, F);
  if (cache.l !== cache.hk) throw new Error('FFT rep mismatch');
  const high = quotientTruncatedProduct(a.slice(n), cache.h, cache.l, 2 * n - 3, F),
    q = high.slice(n - 2, 2 * n - 3);
  quotientCheckSize(cache.k, F);
  if (cache.k !== cache.fk) throw new Error('FFT rep mismatch');
  return k.norm(
    k.sub(quotientFold(a, cache.k, F), quotientFold(k.mul(q, cache.f), cache.k, F)).slice(0, n)
  );
}
/** Native block remainder using an arbitrary-modulus quotient.
 * @see Deviation: NTL arbitrary-modulus polynomial quotient adapters
 */
export function rem(a: readonly bigint[], F: ZZ_pXModulus): bigint[] {
  if (F.n < 0) throw new Error('rem: unitialized modulus');
  const k = F.arithmetic,
    A = k.norm(a),
    n = F.n;
  const reduce = (v: bigint[]): bigint[] => {
    if (v.length <= n) return v.slice();
    if (!F.UseFFT || v.length - 1 - n <= 20) return k.divrem(v, F.f)[1];
    return quotientFFTRemainder(v, F);
  };
  if (A.length <= 2 * n - 1 || !F.UseFFT || A.length - 1 - n <= 20) return reduce(A);
  if (n === 0) throw new Error('negative length in vector::SetLength');
  let buf: bigint[] = [],
    left = A.length;
  while (left > 0) {
    const amount = Math.min(2 * n - 1 - buf.length, left);
    buf = reduce(k.norm([...A.slice(left - amount, left), ...buf]));
    left -= amount;
  }
  return buf;
}
function quotientMulMod(
  a: readonly bigint[],
  b: readonly bigint[] | ZZ_pXMultiplier,
  F: ZZ_pXModulus
): bigint[] {
  const k = F.arithmetic,
    A = k.norm(a),
    cache = quotientCache(F);
  if (b instanceof ZZ_pXMultiplier) {
    if (A.length - 1 >= F.n)
      throw new Error(' bad args to MulMod(ZZ_pX,ZZ_pX,ZZ_pXMultiplier,ZZ_pXModulus)');
    if (!A.length) return [];
    if (!b.UseFFT || !F.UseFFT || A.length - 1 <= 20) return rem(mul(A, b.b, F.p, F.state), F);
    const multiplier = quotientMultiplierCache(b);
    quotientCheckSize(cache.l, F);
    if (cache.l !== multiplier.b1k) throw new Error('FFT rep mismatch');
    const high = quotientTruncatedProduct(
        A,
        multiplier.b1,
        cache.l,
        Math.max(2 ** cache.k, 2 * F.n - 2),
        F
      ),
      q = high.slice(F.n - 1, 2 * F.n - 2);
    if (cache.k !== multiplier.b2k) throw new Error('FFT rep mismatch');
    const product = quotientFold(k.mul(A, k.norm(multiplier.b2)), cache.k, F);
    quotientCheckSize(cache.k, F);
    if (cache.k !== cache.fk) throw new Error('FFT rep mismatch');
    return k.norm(k.sub(product, quotientFold(k.mul(q, cache.f), cache.k, F)).slice(0, F.n));
  }
  if (F.n < 0) throw new Error('MulMod: uninitialized modulus');
  const B = a === b ? A : k.norm(b);
  if (A.length > F.n || B.length > F.n)
    throw new Error('bad args to MulMod(ZZ_pX,ZZ_pX,ZZ_pX,ZZ_pXModulus)');
  if (!A.length || !B.length) return [];
  if (!F.UseFFT || A.length - 1 <= 20 || B.length - 1 <= 20) return rem(mul(A, B, F.p, F.state), F);
  quotientCheckSize(Math.max(quotientFFTSize(A.length + B.length - 1), cache.k), F);
  return quotientFFTRemainder(k.mul(A, B), F);
}
/** Native modular square with a prebuilt arbitrary-modulus quotient.
 * @see Deviation: NTL arbitrary-modulus polynomial quotient adapters
 */
export function SqrMod(a: readonly bigint[], F: ZZ_pXModulus): bigint[] {
  if (F.n < 0) throw new Error('SqrMod: uninitailized modulus');
  const k = F.arithmetic,
    A = k.norm(a),
    cache = quotientCache(F);
  if (A.length > F.n) throw new Error('bad args to SqrMod(ZZ_pX,ZZ_pX,ZZ_pXModulus)');
  if (!F.UseFFT || A.length - 1 <= 20) return rem(mul(A, A, F.p, F.state), F);
  quotientCheckSize(Math.max(quotientFFTSize(2 * A.length - 1), cache.k), F);
  return quotientFFTRemainder(k.mul(A, A), F);
}
