/**
 * @module sage/stats/distributions/discrete_gaussian_lattice
 * @description Discrete Gaussian Distribution over Lattices
 *
 * This module implements the GPV (Gentry-Peikert-Vaikuntanathan) algorithm
 * for sampling from a discrete Gaussian distribution over a lattice, plus
 * Peikert's non-spherical sampler for a covariance matrix `Sigma`.
 *
 * The discrete Gaussian distribution D_{L,sigma,c} over a lattice L is defined by:
 *   P(x) = rho_{sigma,c}(x) / rho_{sigma,c}(L)
 *
 * where rho_{sigma,c}(x) = exp(-||x-c||^2 / (2*sigma^2)) is the Gaussian function.
 *
 * @see Reference: sage/stats/distributions/discrete_gaussian_lattice.py
 * @see [GPV2008] Gentry, Peikert, Vaikuntanathan, "Trapdoors for Hard Lattices...", 2008
 * @see [Pei2010] Peikert, "An Efficient and Parallel Gaussian Sampler for Lattices", 2010
 */

// PARI's `qfrep`, i.e. `minim0(a, borne, gen_0, min_VECSMALL)`
// (`reference/pari/src/basemath/bibli1.c:1649`).  SageMath's
// `_normalisation_factor_zz` delegates the theta series to PARI
// (`Q.__pari__().qfrep(B, 0)`, `discrete_gaussian_lattice.py:340,351`), so we
// delegate to our port of PARI, as required by CLAUDE.md's "Architecture
// Fidelity" rule.
import { qfrep0 as pari_qfrep0 } from '@sagemath-ts/parigp-ts';
import {
  NotImplementedError,
  RuntimeError,
  TypeError as SageTypeError,
  ValueError,
} from '../../errors.js';
import { current_randstate } from '../../misc/randstate.js';
import { Rational } from '../../rings/rational.js';
import { type IntegerLike, toBigInt, toSafeNumber } from '../../types/coercion.js';
import {
  DiscreteGaussianDistributionIntegerSampler,
  type DiscreteGaussianOptions,
} from './discrete_gaussian_integer.js';

/* ====================================================================== */
/* RealField(prec): a minimal binary floating point layer                  */
/*                                                                         */
/* `_normalisation_factor_zz` builds its sum in `RealField(prec)`           */
/* (`discrete_gaussian_lattice.py:325-355`) and its doctest                 */
/*                                                                         */
/*     sage: D = DGL(ZZ^8, 1000)                                            */
/*     sage: round(D._normalisation_factor_zz(prec=100))                    */
/*     1558545456544038969634991553                                         */
/*                                                                         */
/* needs 91 significant bits, so a JavaScript `number` cannot carry the     */
/* answer.  SageMath's `RealField` is MPFR, which is *not* vendored under   */
/* `reference/`; what follows is therefore not a transcription but a        */
/* re-implementation of MPFR's *semantics*: a sign/mantissa/exponent triple */
/* with a mantissa of exactly `p` bits and round-to-nearest, ties-to-even   */
/* on every operation.  Every operation rounds the *exact* result, so the   */
/* four basic operations are correctly rounded by construction.             */
/*                                                                         */
/* @see Deviation: `exp`, `sqrt`, `pi` and `log 2` are computed with 96     */
/* guard bits and then rounded, so their last bit may differ from MPFR's    */
/* (which is correctly rounded).  The printing rules and `round()` follow   */
/* `reference/sage/src/sage/rings/real_mpfr.pyx` exactly.                   */
/* ====================================================================== */

/** Number of bits in `|n|` (0 for `n = 0`). */
function bitLength(n: bigint): number {
  const x = n < 0n ? -n : n;
  if (x === 0n) return 0;
  const s = x.toString(16);
  return (s.length - 1) * 4 + (32 - Math.clz32(Number.parseInt(s[0]!, 16)));
}

/** Integer square root (floor). */
function isqrtBig(n: bigint): bigint {
  if (n < 0n) throw new ValueError('isqrt of a negative number');
  if (n < 2n) return n;
  let x = 1n << BigInt((bitLength(n) + 1) >> 1);
  for (;;) {
    const y = (x + n / x) >> 1n;
    if (y >= x) return x;
    x = y;
  }
}

/**
 * An element of `RealField(p)`: the value is `s * m * 2^e` with `m` having
 * exactly `p` bits (or `s = 0`, `m = 0`).
 */
export class RealNumberMP {
  constructor(
    /** Sign: `-1`, `0` or `1`. */
    readonly s: number,
    /** Mantissa, exactly `p` bits when `s != 0`. */
    readonly m: bigint,
    /** Binary exponent. */
    readonly e: number,
    /** Precision of the parent `RealField`. */
    readonly p: number
  ) {}

  /** `RealNumber.precision()` (`real_mpfr.pyx`). */
  precision(): number {
    return this.p;
  }

  is_zero(): boolean {
    return this.s === 0;
  }

  /** Nearest `double`. */
  toNumber(): number {
    if (this.s === 0) return 0;
    const x = rnSetPrec(this, 53);
    let v = Number(x.m);
    let ee = x.e;
    while (ee < -1000) {
      v *= 2 ** -1000;
      ee += 1000;
    }
    while (ee > 1000) {
      v *= 2 ** 1000;
      ee -= 1000;
    }
    return this.s * v * 2 ** ee;
  }

  /** So that `x / y`, `x < y`, `Math.abs(x)` … keep working. */
  valueOf(): number {
    return this.toNumber();
  }

  /**
   * `RealNumber.round()` (`real_mpfr.pyx:3034-3056`): round to the nearest
   * integer, halfway cases away from zero.
   */
  round(): bigint {
    if (this.s === 0) return 0n;
    const s = BigInt(this.s);
    if (this.e >= 0) return s * (this.m << BigInt(this.e));
    const sh = BigInt(-this.e);
    const q = this.m >> sh;
    const r = this.m & ((1n << sh) - 1n);
    const half = 1n << (sh - 1n);
    return s * (r >= half ? q + 1n : q);
  }

  /** `RealNumber.floor()`. */
  floor(): bigint {
    if (this.s === 0) return 0n;
    if (this.e >= 0) return BigInt(this.s) * (this.m << BigInt(this.e));
    const sh = BigInt(-this.e);
    const q = this.m >> sh;
    if (this.s > 0) return q;
    return this.m === q << sh ? -q : -q - 1n;
  }

  /**
   * `RealNumber.str(truncate=True)`, which is `_repr_`
   * (`real_mpfr.pyx:1897-2149`): print
   * `digits = max(2, floor((prec - 1) * log10(2)))` decimal digits, in
   * scientific notation iff `|exponent - 1| >= 6`.
   */
  str(): string {
    let digits = Math.floor((this.p - 1) * Math.LN2 * Math.LOG10E);
    if (digits < 2) digits = 2;
    if (this.s === 0) {
      // real_mpfr.pyx:2101: "For backwards compatibility, add one extra digit
      // for 0.0"; real_mpfr.pyx:2126 treats the exponent as 1.
      const t = '0'.repeat(digits + 1);
      return `${t.slice(0, 1)}.${t.slice(1)}`;
    }
    const sgn = this.s < 0 ? '-' : '';
    // First guess for the decimal exponent E with |value| in [10^(E-1), 10^E).
    const bl = bitLength(this.m);
    const lead =
      bl > 53 ? Number(this.m >> BigInt(bl - 53)) / 2 ** 52 : Number(this.m) / 2 ** (bl - 1);
    const log10v = (this.e + bl - 1) * Math.LN2 * Math.LOG10E + Math.log10(lead);
    let E = Math.floor(log10v) + 1;
    let t = '';
    for (let iter = 0; iter < 8; iter++) {
      const k = digits - E;
      if (Math.abs(k) > 100000) {
        throw new NotImplementedError(
          'SAGE_NOT_IMPLEMENTED: decimal printing of a RealField element whose ' +
            'decimal exponent exceeds 100000 in absolute value'
        );
      }
      const N = decimalRoundHalfEven(this.m, this.e, k);
      if (N === 0n) {
        // The guess for E was far too large; take a big step down.
        E -= digits;
        continue;
      }
      const str = N.toString();
      if (str.length === digits) {
        t = str;
        break;
      }
      E += str.length - digits;
    }
    if (t === '') {
      throw new RuntimeError('unable to convert an mpfr number to a string');
    }
    // real_mpfr.pyx:2129-2149
    if (Math.abs(E - 1) >= 6) {
      const mant = t.length > 1 ? `${t[0]}.${t.slice(1)}` : t;
      return `${sgn}${mant}e${E - 1}`;
    }
    if (E <= 0) return `${sgn}0.${'0'.repeat(-E)}${t}`;
    if (E >= t.length) return `${sgn}${t}${'0'.repeat(E - t.length)}.`;
    return `${sgn}${t.slice(0, E)}.${t.slice(E)}`;
  }

  toString(): string {
    return this.str();
  }
}

/** `round(m * 2^e * 10^k)` with ties to even, `m >= 0`. */
function decimalRoundHalfEven(m: bigint, e: number, k: number): bigint {
  let num = m;
  let den = 1n;
  if (k >= 0) num *= 10n ** BigInt(k);
  else den *= 10n ** BigInt(-k);
  if (e >= 0) num <<= BigInt(e);
  else den <<= BigInt(-e);
  let q = num / den;
  const r = num - q * den;
  const twice = 2n * r;
  if (twice > den || (twice === den && (q & 1n) === 1n)) q += 1n;
  return q;
}

/** Compare two non-negative bigints. */
function cmpBig(a: bigint, b: bigint): number {
  return a < b ? -1 : a > b ? 1 : 0;
}

/**
 * Round `sign * (num / den) * 2^extraExp` to `p` significant bits
 * (round-to-nearest, ties to even).  `num > 0`, `den > 0`.
 */
function rnMake(sign: number, num: bigint, den: bigint, extraExp: number, p: number): RealNumberMP {
  if (num === 0n) return rnZero(p);
  // e0 with 2^e0 <= num/den < 2^(e0+1)
  let e0 = bitLength(num) - bitLength(den);
  const cmpShifted = (k: number): number =>
    k >= 0 ? cmpBig(num, den << BigInt(k)) : cmpBig(num << BigInt(-k), den);
  if (cmpShifted(e0) < 0) e0 -= 1;
  const sh = p - 1 - e0;
  let N = num;
  let D = den;
  if (sh >= 0) N = num << BigInt(sh);
  else D = den << BigInt(-sh);
  let q = N / D;
  const r = N - q * D;
  const twice = 2n * r;
  if (twice > D || (twice === D && (q & 1n) === 1n)) q += 1n;
  let e = e0 - (p - 1) + extraExp;
  if (bitLength(q) > p) {
    // q === 2^p, so the shift is exact.
    q >>= 1n;
    e += 1;
  }
  return new RealNumberMP(sign < 0 ? -1 : 1, q, e, p);
}

function rnZero(p: number): RealNumberMP {
  return new RealNumberMP(0, 0n, 0, p);
}

function rnOne(p: number): RealNumberMP {
  return new RealNumberMP(1, 1n << BigInt(p - 1), -(p - 1), p);
}

function rnFromBigInt(n: bigint, p: number): RealNumberMP {
  if (n === 0n) return rnZero(p);
  return rnMake(n < 0n ? -1 : 1, n < 0n ? -n : n, 1n, 0, p);
}

function rnFromRational(q: Rational, p: number): RealNumberMP {
  const num = q.numerator;
  if (num === 0n) return rnZero(p);
  return rnMake(num < 0n ? -1 : 1, num < 0n ? -num : num, q.denominator, 0, p);
}

/** Change the working precision (exact when `p >= x.p`). */
function rnSetPrec(x: RealNumberMP, p: number): RealNumberMP {
  if (x.p === p) return x;
  if (x.s === 0) return rnZero(p);
  return rnMake(x.s, x.m, 1n, x.e, p);
}

/** Multiply by `2^k`, exactly. */
function rnScalePow2(x: RealNumberMP, k: number): RealNumberMP {
  if (x.s === 0) return x;
  return new RealNumberMP(x.s, x.m, x.e + k, x.p);
}

function rnNeg(x: RealNumberMP): RealNumberMP {
  return new RealNumberMP(-x.s, x.m, x.e, x.p);
}

function rnAdd(x: RealNumberMP, y: RealNumberMP, p: number): RealNumberMP {
  if (x.s === 0) return rnSetPrec(y, p);
  if (y.s === 0) return rnSetPrec(x, p);
  // If one operand is below half the smallest gap around the other, the sum
  // rounds back to that operand.  This keeps `1 + 2^-28500000` cheap.
  if (y.e + y.p <= x.e + x.p - p - 2) return rnSetPrec(x, p);
  if (x.e + x.p <= y.e + y.p - p - 2) return rnSetPrec(y, p);
  const emin = Math.min(x.e, y.e);
  const nx = BigInt(x.s) * (x.m << BigInt(x.e - emin));
  const ny = BigInt(y.s) * (y.m << BigInt(y.e - emin));
  const n = nx + ny;
  if (n === 0n) return rnZero(p);
  return rnMake(n < 0n ? -1 : 1, n < 0n ? -n : n, 1n, emin, p);
}

function rnSub(x: RealNumberMP, y: RealNumberMP, p: number): RealNumberMP {
  return rnAdd(x, rnNeg(y), p);
}

function rnMul(x: RealNumberMP, y: RealNumberMP, p: number): RealNumberMP {
  if (x.s === 0 || y.s === 0) return rnZero(p);
  return rnMake(x.s * y.s, x.m * y.m, 1n, x.e + y.e, p);
}

function rnDiv(x: RealNumberMP, y: RealNumberMP, p: number): RealNumberMP {
  if (y.s === 0) throw new ValueError('division by zero');
  if (x.s === 0) return rnZero(p);
  return rnMake(x.s * y.s, x.m, y.m, x.e - y.e, p);
}

function rnCmp(x: RealNumberMP, y: RealNumberMP): number {
  if (x.s !== y.s) return x.s < y.s ? -1 : 1;
  if (x.s === 0) return 0;
  const magnitude = (a: RealNumberMP): number => a.e + a.p;
  const mx = magnitude(x);
  const my = magnitude(y);
  if (mx !== my) return (mx < my ? -1 : 1) * x.s;
  // Same binade: align (the shift is bounded by |p_x - p_y|).
  const emin = Math.min(x.e, y.e);
  const a = x.m << BigInt(x.e - emin);
  const b = y.m << BigInt(y.e - emin);
  return cmpBig(a, b) * x.s;
}

/** Structural equality, i.e. `RealField` equality of two normalised elements. */
function rnEq(x: RealNumberMP, y: RealNumberMP): boolean {
  return x.s === y.s && x.m === y.m && (x.s === 0 || x.e === y.e);
}

function rnPowInt(x: RealNumberMP, k: number, p: number): RealNumberMP {
  if (k === 0) return rnOne(p);
  if (k < 0) return rnDiv(rnOne(p), rnPowInt(x, -k, p), p);
  let base = x;
  let result: RealNumberMP | null = null;
  let n = k;
  while (n > 0) {
    if (n & 1) result = result === null ? base : rnMul(result, base, p);
    n >>= 1;
    if (n > 0) base = rnMul(base, base, p);
  }
  return rnSetPrec(result!, p);
}

/** `atan(1/x) * scale` in fixed point. */
function atanInvFixed(x: bigint, scale: bigint): bigint {
  let power = scale / x;
  const x2 = x * x;
  let total = power;
  let n = 1n;
  let sign = 1n;
  while (power !== 0n) {
    power /= x2;
    n += 2n;
    sign = -sign;
    total += sign * (power / n);
  }
  return total;
}

/** `atanh(1/x) * scale` in fixed point. */
function atanhInvFixed(x: bigint, scale: bigint): bigint {
  let power = scale / x;
  const x2 = x * x;
  let total = power;
  let n = 1n;
  while (power !== 0n) {
    power /= x2;
    n += 2n;
    total += power / n;
  }
  return total;
}

const _piCache = new Map<number, RealNumberMP>();

/** `RealField(p).pi()`, by Machin's formula. */
function rnPi(p: number): RealNumberMP {
  const cached = _piCache.get(p);
  if (cached !== undefined) return cached;
  const g = p + 96;
  const scale = 1n << BigInt(g);
  const value = 16n * atanInvFixed(5n, scale) - 4n * atanInvFixed(239n, scale);
  const out = rnMake(1, value, 1n, -g, p);
  _piCache.set(p, out);
  return out;
}

const _log2Cache = new Map<number, RealNumberMP>();

/** `RealField(p).log2()`, from `log 2 = 2 atanh(1/3)`. */
function rnLog2(p: number): RealNumberMP {
  const cached = _log2Cache.get(p);
  if (cached !== undefined) return cached;
  const g = p + 96;
  const scale = 1n << BigInt(g);
  const out = rnMake(1, 2n * atanhInvFixed(3n, scale), 1n, -g, p);
  _log2Cache.set(p, out);
  return out;
}

/** `RealNumber.sqrt()` for a non-negative argument. */
function rnSqrt(x: RealNumberMP, p: number): RealNumberMP {
  if (x.s === 0) return rnZero(p);
  if (x.s < 0) throw new ValueError('sqrt of a negative number');
  let m = x.m;
  let e = x.e;
  if (((e % 2) + 2) % 2 !== 0) {
    m <<= 1n;
    e -= 1;
  }
  const target = 2 * p + 128;
  const t = Math.max(0, Math.ceil((target - bitLength(m)) / 2));
  const q = isqrtBig(m << BigInt(2 * t));
  return rnMake(1, q, 1n, e / 2 - t, p);
}

/**
 * `RealNumber.exp()`: `exp(x) = 2^k exp(r)` with `k = round(x / log 2)` and
 * `|r| <= (log 2)/2`, then `exp(r) = (exp(r/2^j))^(2^j)` evaluated by its
 * Taylor series.
 */
function rnExp(x: RealNumberMP, p: number): RealNumberMP {
  if (x.s === 0) return rnOne(p);
  const k = rnDiv(rnSetPrec(x, p + 64), rnLog2(p + 64), p + 64).round();
  const kb = bitLength(k);
  const wp = p + 96 + kb;
  const r =
    k === 0n
      ? rnSetPrec(x, wp)
      : rnSub(rnSetPrec(x, wp), rnMul(rnFromBigInt(k, wp), rnLog2(wp), wp), wp);
  const j = r.s === 0 ? 0 : Math.max(0, r.e + r.p + 16);
  const rr = rnScalePow2(r, -j);
  let term = rnOne(wp);
  let sum = rnOne(wp);
  for (let i = 1; i < 1000000; i++) {
    term = rnDiv(rnMul(term, rr, wp), rnFromBigInt(BigInt(i), wp), wp);
    if (term.s === 0) break;
    sum = rnAdd(sum, term, wp);
    if (term.e + term.p < sum.e + sum.p - wp - 4) break;
  }
  for (let t = 0; t < j; t++) sum = rnMul(sum, sum, wp);
  return rnSetPrec(rnScalePow2(sum, Number(k)), p);
}

/**
 * Test-only surface for the `RealField(prec)` layer above.  Not part of the
 * SageMath API; do not re-export from `index.ts`.
 */
export const _mp = {
  bitLength,
  isqrtBig,
  rnZero,
  rnOne,
  rnFromBigInt,
  rnFromRational,
  rnSetPrec,
  rnScalePow2,
  rnNeg,
  rnAdd,
  rnSub,
  rnMul,
  rnDiv,
  rnCmp,
  rnEq,
  rnPowInt,
  rnPi,
  rnLog2,
  rnSqrt,
  rnExp,
  exactFromDouble: (x: number): Rational => exactFromDouble(x),
};

/**
 * The exact value of a JavaScript `number`, i.e. of the `RealField(53)`
 * element SageMath would build from the same literal.
 */
function exactFromDouble(x: number): Rational {
  if (!Number.isFinite(x)) {
    throw new SageTypeError(`expected a finite number, got ${x}`);
  }
  if (x === 0) return Rational.zero();
  const view = new DataView(new ArrayBuffer(8));
  view.setFloat64(0, x);
  const hi = view.getUint32(0);
  const lo = view.getUint32(4);
  const sign = hi >>> 31 ? -1n : 1n;
  const expBits = (hi >>> 20) & 0x7ff;
  let mant = (BigInt(hi & 0xfffff) << 32n) | BigInt(lo);
  let e: number;
  if (expBits === 0) {
    e = -1074;
  } else {
    mant |= 1n << 52n;
    e = expBits - 1075;
  }
  const num = sign * mant;
  return e >= 0 ? new Rational(num << BigInt(e), 1n) : new Rational(num, 1n << BigInt(-e));
}

/**
 * Convert a basis/center entry to an exact rational.
 *
 * Sage's sampler takes a matrix over an exact ring (typically `ZZ` or `QQ`)
 * and does all lattice arithmetic exactly; only `sigma` lives in `RealField`.
 * We mirror that by keeping every basis and center entry as a {@link Rational}.
 */
function toRationalEntry(x: number | bigint | Rational): Rational {
  if (x instanceof Rational) {
    return x;
  }
  if (typeof x === 'bigint') {
    return new Rational(x, 1n);
  }
  if (typeof x !== 'number' || !Number.isFinite(x)) {
    throw new SageTypeError(`entries must be finite numbers, got ${x}`);
  }
  if (Number.isInteger(x)) {
    return new Rational(BigInt(x), 1n);
  }
  const str = x.toString();
  const eIndex = str.search(/[eE]/);
  if (eIndex < 0) {
    return Rational.from(str);
  }
  // Exponential notation, e.g. "1.5e-7": expand it exactly in base 10.
  const mantissa = Rational.from(str.slice(0, eIndex));
  const exponent = Number.parseInt(str.slice(eIndex + 1), 10);
  const power = new Rational(10n ** BigInt(Math.abs(exponent)), 1n);
  return exponent >= 0 ? mantissa.mul(power) : mantissa.div(power);
}

/**
 * Render a `RealField(53)` element the way Sage does: 15 significant digits.
 */
function rrRepr(x: number): string {
  if (x === 0) return '0.000000000000000';
  if (!Number.isFinite(x)) return String(x);
  return x.toPrecision(15);
}

/**
 * Render a matrix the way Sage does: one bracketed row per line, each column
 * padded to its widest entry and right-justified.
 */
function matrixRepr(rows: string[][]): string {
  if (rows.length === 0) return '[]';
  const ncols = rows[0]!.length;
  const widths: number[] = [];
  for (let j = 0; j < ncols; j++) {
    widths.push(Math.max(...rows.map((r) => r[j]!.length)));
  }
  return rows.map((r) => `[${r.map((e, j) => e.padStart(widths[j]!)).join(' ')}]`).join('\n');
}

/**
 * Result of Gram-Schmidt orthogonalization (exact).
 */
interface GramSchmidtResult {
  /** The orthogonalized basis B* (rows are b_1*, ..., b_n*) */
  bStar: Rational[][];
  /** The mu coefficients matrix where mu[i][j] = <b_i, b_j*> / <b_j*, b_j*> */
  mu: Rational[][];
  /** The squared norms of the orthogonal vectors |b_i*|^2 */
  bStarNormsSq: Rational[];
}

/**
 * Compute the exact Gram-Schmidt orthogonalization of a basis.
 *
 * This mirrors Sage's `B.gram_schmidt()` on an exact matrix
 * (`discrete_gaussian_lattice.py:587`), which returns rationals — not the
 * floating-point `gram_schmidt(..., orthonormal=True)` variant.
 *
 * @param basis - A matrix whose rows form the basis
 * @returns GramSchmidtResult containing orthogonal basis, mu matrix, and squared norms
 */
function gramSchmidt(basis: Rational[][]): GramSchmidtResult {
  const n = basis.length;
  if (n === 0) {
    return { bStar: [], mu: [], bStarNormsSq: [] };
  }
  const m = basis[0]!.length;

  const bStar: Rational[][] = [];
  const mu: Rational[][] = [];
  const bStarNormsSq: Rational[] = [];

  for (let i = 0; i < n; i++) {
    mu.push(new Array(n).fill(Rational.zero()));
    mu[i]![i] = Rational.one();

    // Start with b_i
    const bi = basis[i]!;
    const biStar = [...bi];

    // Subtract projections onto previous orthogonal vectors
    for (let j = 0; j < i; j++) {
      // mu[i][j] = <b_i, b_j*> / <b_j*, b_j*>
      let dotProduct = Rational.zero();
      for (let k = 0; k < m; k++) {
        dotProduct = dotProduct.add(bi[k]!.mul(bStar[j]![k]!));
      }
      if (bStarNormsSq[j]!.isZero()) {
        throw new ValueError('basis vectors must be linearly independent');
      }
      mu[i]![j] = dotProduct.div(bStarNormsSq[j]!);

      // b_i* = b_i* - mu[i][j] * b_j*
      for (let k = 0; k < m; k++) {
        biStar[k] = biStar[k]!.sub(mu[i]![j]!.mul(bStar[j]![k]!));
      }
    }

    bStar.push(biStar);

    // Compute |b_i*|^2
    let normSq = Rational.zero();
    for (let k = 0; k < m; k++) {
      normSq = normSq.add(biStar[k]!.mul(biStar[k]!));
    }
    bStarNormsSq.push(normSq);
  }

  return { bStar, mu, bStarNormsSq };
}

/** Exact non-negative integer power of a rational. */
function ratPowInt(x: Rational, k: number): Rational {
  let result = Rational.one();
  let base = x;
  let n = k;
  while (n > 0) {
    if (n & 1) result = result.mul(base);
    n >>= 1;
    if (n > 0) base = base.mul(base);
  }
  return result;
}

/** Exact dot product of two rational vectors. */
function ratDot(a: Rational[], b: Rational[]): Rational {
  let sum = Rational.zero();
  for (let i = 0; i < a.length; i++) {
    sum = sum.add(a[i]!.mul(b[i]!));
  }
  return sum;
}

/** Exact matrix product. */
function ratMatMul(A: Rational[][], B: Rational[][]): Rational[][] {
  const n = A.length;
  const k = B.length;
  const m = B[0]!.length;
  const out: Rational[][] = [];
  for (let i = 0; i < n; i++) {
    const row: Rational[] = [];
    for (let j = 0; j < m; j++) {
      let s = Rational.zero();
      for (let t = 0; t < k; t++) {
        s = s.add(A[i]![t]!.mul(B[t]![j]!));
      }
      row.push(s);
    }
    out.push(row);
  }
  return out;
}

/** Exact transpose. */
function ratTranspose(A: Rational[][]): Rational[][] {
  const n = A.length;
  const m = A[0]!.length;
  const out: Rational[][] = [];
  for (let j = 0; j < m; j++) {
    const row: Rational[] = [];
    for (let i = 0; i < n; i++) row.push(A[i]![j]!);
    out.push(row);
  }
  return out;
}

/** Exact inverse by Gauss-Jordan elimination over QQ. */
function ratInverse(A: Rational[][]): Rational[][] {
  const n = A.length;
  if (n === 0 || A[0]!.length !== n) {
    throw new ValueError('matrix must be square to be inverted');
  }
  const M = A.map((r, i) => [
    ...r,
    ...Array.from({ length: n }, (_, j) => (i === j ? Rational.one() : Rational.zero())),
  ]);
  for (let col = 0; col < n; col++) {
    let pivot = -1;
    for (let row = col; row < n; row++) {
      if (!M[row]![col]!.isZero()) {
        pivot = row;
        break;
      }
    }
    if (pivot < 0) {
      throw new ValueError('matrix is not invertible');
    }
    if (pivot !== col) {
      const tmp = M[col]!;
      M[col] = M[pivot]!;
      M[pivot] = tmp;
    }
    const inv = M[col]![col]!.inv();
    M[col] = M[col]!.map((x) => x.mul(inv));
    for (let row = 0; row < n; row++) {
      if (row === col) continue;
      const factor = M[row]![col]!;
      if (factor.isZero()) continue;
      M[row] = M[row]!.map((x, j) => x.sub(factor.mul(M[col]![j]!)));
    }
  }
  return M.map((r) => r.slice(n));
}

/** Exact determinant by fraction-free-ish Gaussian elimination over QQ. */
function ratDet(A: Rational[][]): Rational {
  const n = A.length;
  if (n === 0) return Rational.one();
  const M = A.map((r) => [...r]);
  let det = Rational.one();
  for (let col = 0; col < n; col++) {
    let pivot = -1;
    for (let row = col; row < n; row++) {
      if (!M[row]![col]!.isZero()) {
        pivot = row;
        break;
      }
    }
    if (pivot < 0) return Rational.zero();
    if (pivot !== col) {
      const tmp = M[col]!;
      M[col] = M[pivot]!;
      M[pivot] = tmp;
      det = det.neg();
    }
    det = det.mul(M[col]![col]!);
    const inv = M[col]![col]!.inv();
    for (let row = col + 1; row < n; row++) {
      const factor = M[row]![col]!.mul(inv);
      if (factor.isZero()) continue;
      M[row] = M[row]!.map((x, j) => x.sub(factor.mul(M[col]![j]!)));
    }
  }
  return det;
}

/**
 * Solve `x * A = b` exactly (Sage's `A.solve_left(b)`), for a square
 * invertible `A`.
 */
function ratSolveLeft(A: Rational[][], b: Rational[]): Rational[] {
  const Ainv = ratInverse(A);
  return Ainv[0]!.map((_, j) =>
    b.reduce((acc, bi, i) => acc.add(bi.mul(Ainv[i]![j]!)), Rational.zero())
  );
}

/** Double-precision matrix product. */
function rrMatMul(A: number[][], B: number[][]): number[][] {
  const n = A.length;
  const k = B.length;
  const m = B[0]!.length;
  const out: number[][] = [];
  for (let i = 0; i < n; i++) {
    const row: number[] = [];
    for (let j = 0; j < m; j++) {
      let s = 0;
      for (let t = 0; t < k; t++) s += A[i]![t]! * B[t]![j]!;
      row.push(s);
    }
    out.push(row);
  }
  return out;
}

/** Double-precision matrix inverse via Gauss-Jordan with partial pivoting. */
function rrInverse(A: number[][]): number[][] {
  const n = A.length;
  const M = A.map((r, i) => [...r, ...Array.from({ length: n }, (_, j) => (i === j ? 1 : 0))]);
  for (let col = 0; col < n; col++) {
    let pivot = col;
    for (let row = col + 1; row < n; row++) {
      if (Math.abs(M[row]![col]!) > Math.abs(M[pivot]![col]!)) pivot = row;
    }
    if (M[pivot]![col] === 0) {
      throw new ValueError('matrix is not invertible');
    }
    if (pivot !== col) {
      const tmp = M[col]!;
      M[col] = M[pivot]!;
      M[pivot] = tmp;
    }
    const d = M[col]![col]!;
    M[col] = M[col]!.map((x) => x / d);
    for (let row = 0; row < n; row++) {
      if (row === col) continue;
      const f = M[row]![col]!;
      if (f === 0) continue;
      M[row] = M[row]!.map((x, j) => x - f * M[col]![j]!);
    }
  }
  return M.map((r) => r.slice(n));
}

/**
 * Cholesky decomposition `A = L L^T` with `L` lower triangular.
 *
 * Mirrors Sage's `Matrix.cholesky()` over `RealField`, which raises
 * `ValueError` when the matrix is not positive definite.
 */
function rrCholesky(A: number[][]): number[][] {
  const n = A.length;
  const L: number[][] = Array.from({ length: n }, () => new Array<number>(n).fill(0));
  for (let i = 0; i < n; i++) {
    for (let j = 0; j <= i; j++) {
      let s = A[i]![j]!;
      for (let k = 0; k < j; k++) s -= L[i]![k]! * L[j]![k]!;
      if (i === j) {
        if (s <= 0) {
          throw new ValueError('matrix is not positive definite');
        }
        L[i]![j] = Math.sqrt(s);
      } else {
        L[i]![j] = s / L[j]![j]!;
      }
    }
  }
  return L;
}

/**
 * Sage's `Matrix.is_positive_definite()` for a real symmetric matrix.
 */
function rrIsPositiveDefinite(A: number[][]): boolean {
  const n = A.length;
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      if (Math.abs(A[i]![j]! - A[j]![i]!) > 0) return false;
    }
  }
  try {
    rrCholesky(A);
    return true;
  } catch {
    return false;
  }
}

/**
 * LLL reduction of an integral row basis, mirroring `B.LLL()`.
 *
 * @see Deviation: Sage delegates to fpLLL with `delta = 0.99, eta = 0.501`;
 * this is the classical exact-rational LLL with the same `delta`.  The two
 * agree on the reduced *lattice* but may return different reduced bases for
 * inputs with several near-optimal reductions.  It is used only by
 * {@link DiscreteGaussianDistributionLatticeSampler._normalisation_factor_zz}'s
 * non-spherical branch, which needs it to pick an enumeration centre.
 */
function lllReduce(basis: Rational[][]): Rational[][] {
  const delta = new Rational(99n, 100n);
  const B = basis.map((r) => [...r]);
  const n = B.length;
  if (n <= 1) return B;
  let gs = gramSchmidt(B);
  let k = 1;
  let guard = 0;
  const maxIter = 1000 * n * n + 1000;
  while (k < n && guard++ < maxIter) {
    for (let j = k - 1; j >= 0; j--) {
      const mu = gs.mu[k]![j]!;
      const q = mu.round();
      if (q !== 0n) {
        const qr = new Rational(q, 1n);
        B[k] = B[k]!.map((x, t) => x.sub(qr.mul(B[j]![t]!)));
        gs = gramSchmidt(B);
      }
    }
    const mukk1 = gs.mu[k]![k - 1]!;
    const lhs = gs.bStarNormsSq[k]!;
    const rhs = delta.sub(mukk1.mul(mukk1)).mul(gs.bStarNormsSq[k - 1]!);
    if (lhs.cmp(rhs) >= 0) {
      k++;
    } else {
      const tmp = B[k]!;
      B[k] = B[k - 1]!;
      B[k - 1] = tmp;
      gs = gramSchmidt(B);
      k = Math.max(k - 1, 1);
    }
  }
  return B;
}

/**
 * Theta series of a positive definite integral quadratic form.
 *
 * `qfrep(Q, B)` returns, for `m = 1..floor(B)`, *half* the number of integer
 * vectors `x != 0` with `x Q x^T = m` — that is one representative of each
 * `{x, -x}` pair.  This is PARI's `qfrep(Q, B, 0)`, which SageMath calls at
 * `discrete_gaussian_lattice.py:340` and `:351` as `Q.__pari__().qfrep(B, 0)`.
 *
 * This is a thin adapter over `parigp-ts`'s port of PARI's
 * `qfrep0 -> minim0 -> minim0_dolll` (`reference/pari/src/basemath/bibli1.c`
 * `:1649` and `:1299-1462`); the enumeration itself lives there, where it
 * belongs.
 *
 * @param Q     the (symmetric, positive definite, **integral**) Gram matrix,
 *              given by rows.  Entries must be integers: PARI's `qfrep` takes
 *              a `t_MAT` of `t_INT`s, and SageMath only reaches this call
 *              after checking `self.B.base_ring() == ZZ`.
 * @param bound the bound `B`; non-integral bounds are floored (PARI `gfloor`).
 */
export function qfrep(Q: (number | bigint | Rational)[][], bound: number | bigint): bigint[] {
  const n = Q.length;
  const a: bigint[][] = [];
  for (let j = 0; j < n; j++) {
    const col: bigint[] = [];
    for (let i = 0; i < n; i++) {
      const entry = toRationalEntry(Q[i]![j]!);
      if (!entry.isInteger()) {
        throw new SageTypeError('incorrect type in qfminim');
      }
      col.push(entry.numerator);
    }
    a.push(col);
  }
  // PARI's `t_MAT` is column-major; the form is symmetric, so `a` is `Q`.
  return pari_qfrep0(a, bound, 0);
}

/**
 * Options for constructing a discrete Gaussian lattice sampler.
 */
export interface DiscreteGaussianLatticeOptions {
  /**
   * Gaussian parameter, one of
   * - a real number `sigma > 0` (spherical),
   * - a positive definite matrix `Sigma` (non-spherical), or
   * - any matrix-like `S`, equivalent to `Sigma = S S^T`, when
   *   {@link sigma_basis} is set.
   */
  sigma: number | number[][] | Rational[][];

  /**
   * Center of the distribution (default: origin).
   * Should be a vector of the same dimension as the lattice.
   */
  c?: number[] | bigint[] | Rational[];

  /**
   * Tail cutoff parameter tau >= 1 (default: 6).
   *
   * @see Deviation: extra knob, forwarded to the per-coordinate integer
   * samplers.  Sage's lattice sampler always leaves `tau` at its default.
   */
  tau?: IntegerLike;

  /**
   * Rounding parameter `r` as defined in [Pei2010]_; ignored for a spherical
   * Gaussian parameter.  If not provided, set to the maximal possible value
   * such that `Sigma - r^2 B B^T` is positive definite.
   */
  r?: number;

  /** Bit precision `>= 53`. */
  precision?: number;

  /**
   * When set, `sigma` is treated as a (row) basis, i.e. the covariance matrix
   * is `Sigma = S S^T`.
   */
  sigma_basis?: boolean;
}

/**
 * Discrete Gaussian Distribution Sampler over a Lattice.
 *
 * Implements the GPV algorithm for a scalar `sigma` and Peikert's sampler for
 * a covariance matrix `Sigma`.
 *
 * @example
 * ```typescript
 * const D = new DiscreteGaussianDistributionLatticeSampler([[1, 0], [0, 1]], { sigma: 3 });
 * const v = D.sample(); // e.g., [2n, -1n]
 * ```
 *
 * @see Reference: sage/stats/distributions/discrete_gaussian_lattice.py
 */
export class DiscreteGaussianDistributionLatticeSampler {
  /**
   * The lattice basis (rows are basis vectors), as floating point.
   *
   * This is a lossy view kept for convenience; all sampling arithmetic uses
   * {@link basisExact}.
   */
  public readonly basis: number[][];

  /**
   * The lattice basis (rows are basis vectors), exactly.  Sage calls this `B`.
   */
  public readonly basisExact: Rational[][];

  /** Sage's `self.B`: the exact (row) basis. */
  public readonly B: Rational[][];

  /** Sage's `self.Q = B * B.T` (`discrete_gaussian_lattice.py:586`). */
  public readonly Q: Rational[][];

  /** Sage's `self.n = B.ncols()`. */
  public readonly n: number;

  /**
   * Whether every basis entry is an integer, i.e. the lattice sits inside Z^d.
   *
   * When false, {@link sample} cannot return integer coordinates and
   * {@link sampleExact} must be used instead.
   */
  public readonly isIntegral: boolean;

  /**
   * Whether the Gaussian parameter is a scalar (`true`) or a covariance
   * matrix (`false`).  Sage: `self.is_spherical`.
   */
  public readonly is_spherical: boolean;

  /**
   * Tail cutoff parameter.
   */
  public readonly tau: number;

  /**
   * Dimension of the lattice (number of basis vectors).
   */
  public readonly rank: number;

  /**
   * Dimension of the ambient space.
   */
  public readonly degree: number;

  /** Working precision, Sage's `compute_precision(precision, sigma)`. */
  public readonly precision: number;

  /**
   * Offline samples of `B^{-1} D_1`, used by the non-spherical sampler
   * (`discrete_gaussian_lattice.py:874-892`).
   */
  public offline_samples: number[][] = [];

  /** Peikert's rounding parameter (non-spherical only). */
  public r: number | null;

  /** Gaussian parameter: a scalar, or a covariance matrix in `RealField`. */
  private _sigma: number | number[][];

  /** Inverse of the covariance matrix (non-spherical only). */
  private sigma_inv: number[][] | null = null;

  /**
   * The covariance matrix, exactly (non-spherical only).
   *
   * @see Deviation: SageMath stores `matrix(RealField(prec), sigma)` and
   * inverts it in floating point (`discrete_gaussian_lattice.py:641`); we keep
   * the exact rational matrix as well so that `f_R` — and hence
   * `_normalisation_factor_zz`'s non-spherical branch — can be evaluated at
   * any requested precision instead of only at 53 bits.
   */
  private _sigmaMatrixExact: Rational[][] | null = null;

  /** Exact inverse of {@link _sigmaMatrixExact}. */
  private _sigmaInvExact: Rational[][] | null = null;

  /** Center. */
  private _c: Rational[];

  /** `c * B^{-1}` (non-spherical only), kept exact. */
  private _c_mul_B_inv: Rational[] | null = null;

  /** Exact inverse of the basis (non-spherical only). */
  private B_inv: Rational[][] | null = null;

  /** `Sigma_2 = Sigma - r^2 Q` Cholesky factor, transposed. */
  private B2: number[][] | null = null;

  /** `B2 * B^{-1}` (non-spherical only). */
  private B2_B_inv: number[][] | null = null;

  /** Whether `c` is in the lattice and `B^* == 1`. */
  private _c_in_lattice_and_lattice_trivial = false;

  /** Per-coordinate samplers used by `_call_simple`. */
  private D: DiscreteGaussianDistributionIntegerSampler[] | null = null;

  /**
   * Gram-Schmidt orthogonalization data (Sage's `self._G`, plus norms).
   */
  private readonly gs: GramSchmidtResult;

  /** Cached `_maximal_r()` (Sage caches it with `@cached_method`). */
  private _maximal_r_cache: number | null = null;

  /**
   * Construct a new discrete Gaussian lattice sampler.
   *
   * @param basis - The lattice basis (rows are basis vectors)
   * @param options - Configuration options
   */
  constructor(
    basis: number[][] | bigint[][] | Rational[][],
    options: DiscreteGaussianLatticeOptions
  ) {
    // Validate basis
    if (!Array.isArray(basis) || basis.length === 0) {
      throw new ValueError('basis must be a non-empty array of vectors');
    }

    // Keep the basis exactly; the number[][] view is only for reporting.
    this.basisExact = (basis as (number | bigint | Rational)[][]).map((row) =>
      row.map((x) => toRationalEntry(x))
    );
    this.B = this.basisExact;
    this.basis = this.basisExact.map((row) => row.map((x) => x.toNumber()));
    this.isIntegral = this.basisExact.every((row) => row.every((x) => x.isInteger()));
    this.rank = this.basisExact.length;
    this.degree = this.basisExact[0]!.length;
    this.n = this.degree;

    // Validate all rows have the same length
    for (const row of this.basisExact) {
      if (row.length !== this.degree) {
        throw new ValueError('all basis vectors must have the same dimension');
      }
    }

    // --- sigma: scalar or covariance matrix (discrete_gaussian_lattice.py:553-571)
    if (options.sigma === undefined || options.sigma === null) {
      throw new SageTypeError('sigma is required');
    }
    this.precision = DiscreteGaussianDistributionLatticeSampler.compute_precision(
      options.precision,
      options.sigma
    );
    if (typeof options.sigma === 'number') {
      if (!Number.isFinite(options.sigma)) {
        throw new SageTypeError(`sigma must be a finite number, got ${options.sigma}`);
      }
      if (options.sigma <= 0) {
        throw new ValueError(`sigma must be > 0, got ${options.sigma}`);
      }
      this._sigma = options.sigma;
      this.is_spherical = true;
    } else {
      let Sx = (options.sigma as (number | bigint | Rational)[][]).map((row) =>
        row.map((x) => toRationalEntry(x))
      );
      // A scaled identity is handled as a scalar (py:564-565).
      const s00 = Sx[0]![0]!;
      const isScalarMatrix = Sx.every((row, i) =>
        row.every((v, j) => v.eq(i === j ? s00 : Rational.zero()))
      );
      if (isScalarMatrix) {
        this._sigma = s00.toNumber();
        this.is_spherical = true;
      } else {
        if (options.sigma_basis) {
          Sx = ratMatMul(Sx, ratTranspose(Sx));
        }
        const S = Sx.map((r) => r.map((x) => x.toNumber()));
        if (!rrIsPositiveDefinite(S)) {
          // Sage raises RuntimeError here (discrete_gaussian_lattice.py:570).
          throw new RuntimeError(
            `Sigma(=${matrixRepr(S.map((r) => r.map(rrRepr)))}) is not positive definite`
          );
        }
        this._sigma = S;
        this._sigmaMatrixExact = Sx;
        this._sigmaInvExact = ratInverse(Sx);
        this.is_spherical = false;
      }
    }

    // Validate and store tau (default: 6)
    const tauValue = options.tau !== undefined ? toSafeNumber(toBigInt(options.tau)) : 6;
    if (tauValue < 1) {
      throw new ValueError(`tau must be >= 1, got ${tauValue}`);
    }
    this.tau = tauValue;

    // Q = B * B^T and the exact Gram-Schmidt orthogonalization.
    this.Q = ratMatMul(this.basisExact, ratTranspose(this.basisExact));
    this.gs = gramSchmidt(this.basisExact);

    this.r = options.r ?? null;
    this._c = new Array(this.degree).fill(Rational.zero());
    this.set_c(options.c ?? null);
  }

  /**
   * Compute precision to use.
   *
   * Reference: `discrete_gaussian_lattice.py:155-188`.
   */
  static compute_precision(precision: number | undefined | null, sigma: unknown): number {
    if (precision === undefined || precision === null) {
      // `try: precision = ZZ(sigma.precision()) / except AttributeError: return 53`
      if (sigma instanceof RealNumberMP) {
        precision = sigma.precision();
      } else {
        // A JavaScript `number` is a `RealField(53)` element and has no
        // `.precision()` of its own, so Sage's `except AttributeError` fires.
        return 53;
      }
    }
    return Math.max(53, precision);
  }

  /**
   * Gaussian parameter `sigma`.
   *
   * Reference: `discrete_gaussian_lattice.py:723-736`.
   */
  sigma(): number | number[][] {
    return this._sigma;
  }

  /**
   * Center `c`.
   *
   * Reference: `discrete_gaussian_lattice.py:738-756`.
   */
  c(): Rational[] {
    return this._c;
  }

  /**
   * Center `c` as a floating-point vector (convenience view).
   */
  cNumeric(): number[] {
    return this._c.map((x) => x.toNumber());
  }

  /**
   * Modify the center `c`, recomputing the cached data.
   *
   * Reference: `discrete_gaussian_lattice.py:758-792`.
   */
  set_c(c: number[] | bigint[] | Rational[] | null): void {
    if (c === null || c === undefined) {
      this._c = new Array(this.degree).fill(Rational.zero());
    } else {
      if (!Array.isArray(c) || c.length !== this.degree) {
        throw new ValueError(`c must be a vector of dimension ${this.degree}`);
      }
      this._c = (c as (number | bigint | Rational)[]).map((x) => toRationalEntry(x));
    }
    this._precompute_data();
  }

  /**
   * Precompute basis data.  Do not call directly.
   *
   * Reference: `discrete_gaussian_lattice.py:597-660`.
   */
  private _precompute_data(): void {
    this.D = null;
    this._c_in_lattice_and_lattice_trivial = false;
    if (this.is_spherical) {
      const sigma = this._sigma as number;
      const cIsZero = this._c.every((x) => x.isZero());
      const gIsIdentity = this._gramSchmidtIsIdentity();
      if (cIsZero && gIsIdentity) {
        this._c_in_lattice_and_lattice_trivial = true;
        const D = new DiscreteGaussianDistributionIntegerSampler({ sigma, tau: this.tau });
        this.D = new Array(this.rank).fill(D);
      } else if (gIsIdentity) {
        // w = B.solve_left(c); trivial iff w is integral.
        let w: Rational[] | null = null;
        try {
          w = ratSolveLeft(this.basisExact, this._c);
        } catch {
          w = null;
        }
        if (w?.every((x) => x.isInteger())) {
          this._c_in_lattice_and_lattice_trivial = true;
          const D: DiscreteGaussianDistributionIntegerSampler[] = [];
          for (let i = 0; i < this.rank; i++) {
            const sigma_ = sigma / Math.sqrt(this.gs.bStarNormsSq[i]!.toNumber());
            D.push(
              new DiscreteGaussianDistributionIntegerSampler({ sigma: sigma_, tau: this.tau })
            );
          }
          this.D = D;
        }
      }
      return;
    }

    // Non-spherical: [Pei2010] (discrete_gaussian_lattice.py:636-660).
    const Sigma = this._sigma as number[][];
    this.offline_samples = [];
    this.B_inv = ratInverse(this.basisExact);
    this.sigma_inv = this._sigmaInvExact!.map((r) => r.map((x) => x.toNumber()));
    this._c_mul_B_inv = this._c.map((_, j) =>
      this._c.reduce((acc, ci, i) => acc.add(ci.mul(this.B_inv![i]![j]!)), Rational.zero())
    );

    if (this.r === null) {
      this.r = this._maximal_r() * 0.9999;
    }
    const r = this.r;
    const Qd = this.Q.map((row) => row.map((x) => x.toNumber()));
    const Sigma2 = Sigma.map((row, i) => row.map((v, j) => v - r * r * Qd[i]![j]!));
    let L: number[][];
    try {
      L = rrCholesky(Sigma2);
    } catch {
      throw new ValueError(
        `Σ₂ is not positive definite. Is your r(=${rrRepr(r)}) too large? It should be at most ${rrRepr(this._maximal_r())}`
      );
    }
    // Sage: self.B2 = Sigma2.cholesky().T
    this.B2 = L[0]!.map((_, j) => L.map((row) => row[j]!));
    this.B2_B_inv = rrMatMul(
      this.B2,
      this.B_inv.map((row) => row.map((x) => x.toNumber()))
    );
  }

  /** Whether the exact Gram-Schmidt matrix `_G` equals the identity. */
  private _gramSchmidtIsIdentity(): boolean {
    if (this.rank !== this.degree) return false;
    for (let i = 0; i < this.rank; i++) {
      for (let j = 0; j < this.degree; j++) {
        const want = i === j ? 1n : 0n;
        const v = this.gs.bStar[i]![j]!;
        if (!v.isInteger() || v.numerator !== want) return false;
      }
    }
    return true;
  }

  /**
   * The largest `r > 0` such that `Sigma - r^2 B B^T` is positive definite,
   * found by power iteration.
   *
   * Reference: `discrete_gaussian_lattice.py:357-389`.
   */
  _maximal_r(): number {
    if (!this.is_spherical) {
      if (this._maximal_r_cache !== null) return this._maximal_r_cache;
      const Sigma = this._sigma as number[][];
      const Qd = this.Q.map((row) => row.map((x) => x.toNumber()));
      // Q = self.Q / self._sigma  (i.e. Q * Sigma^{-1})
      const Qh = rrMatMul(Qd, rrInverse(Sigma));
      let v = [...Qh[0]!];
      const apply = (m: number[][], x: number[]): number[] =>
        m.map((row) => row.reduce((acc, e, j) => acc + e * x[j]!, 0));
      const norm = (x: number[]) => Math.sqrt(x.reduce((a, e) => a + e * e, 0));
      for (let cnt = 0; cnt < 10000; cnt++) {
        const w = apply(Qh, v);
        const nw = norm(w);
        const nv = w.map((e) => e / nw);
        const diff = Math.sqrt(nv.reduce((a, e, i) => a + (e - v[i]!) ** 2, 0));
        if (diff < 1e-12) {
          v = nv;
          break;
        }
        v = nv;
      }
      const res = Math.sqrt(v[0]! / apply(Qh, v)[0]!);
      this._maximal_r_cache = res;
      return res;
    }
    throw new ValueError('_maximal_r is only defined for non-spherical samplers');
  }

  /**
   * Compute the Gaussian `rho_{Lambda, c, Sigma}`.
   *
   * Reference: `discrete_gaussian_lattice.py:698-721`.
   */
  f(x: number[] | bigint[] | Rational[]): number {
    const xr = (x as (number | bigint | Rational)[]).map((v) => toRationalEntry(v));
    if (xr.length !== this.n) {
      throw new ValueError(`expected a vector of dimension ${this.n}`);
    }
    const d = xr.map((v, i) => v.sub(this._c[i]!));
    if (this.is_spherical) {
      const sigma = this._sigma as number;
      const normSq = ratDot(d, d).toNumber();
      return Math.exp(-normSq / (2 * sigma * sigma));
    }
    const dn = d.map((v) => v.toNumber());
    const S = this.sigma_inv!;
    let q = 0;
    for (let i = 0; i < dn.length; i++) {
      for (let j = 0; j < dn.length; j++) q += dn[i]! * S[i]![j]! * dn[j]!;
    }
    return Math.exp(-q / 2);
  }

  /**
   * The Gaussian `rho_{Lambda, c, Sigma}` as an element of `RealField(p)`.
   *
   * This is Sage's {@link f} (`discrete_gaussian_lattice.py:698-721`).  Both of
   * its branches produce a `RealField(self._RR.prec())` element — `x.norm()**2`
   * is an exact Integer and `self.sigma_inv` a `RealField` matrix, so the whole
   * expression is evaluated at *this sampler's* precision and there is no
   * precision argument, exactly as upstream has none.
   *
   * @see Deviation: upstream forms the quadratic form `x * sigma_inv * x` with
   * a floating point `Sigma^-1` and floating point matrix arithmetic; we build
   * it exactly over `Q` and round once.  The last bit may therefore differ.
   */
  f_R(x: number[] | bigint[] | Rational[]): RealNumberMP {
    const xr = (x as (number | bigint | Rational)[]).map((v) => toRationalEntry(v));
    if (xr.length !== this.n) {
      throw new ValueError(`expected a vector of dimension ${this.n}`);
    }
    const p = this.precision;
    const d = xr.map((v, i) => v.sub(this._c[i]!));
    if (this.is_spherical) {
      // exp(-x.norm()^2 / (2 sigma^2)); `x.norm()^2` is exact in Sage.
      const normSq = ratDot(d, d);
      return rnExp(rnNeg(rnDiv(rnFromRational(normSq, p), this._twoSigmaSq(), p)), p);
    }
    // exp(-x * sigma_inv * x / 2)
    const S = this._sigmaInvExact!;
    let q = Rational.zero();
    for (let i = 0; i < d.length; i++) {
      for (let j = 0; j < d.length; j++) {
        q = q.add(d[i]!.mul(S[i]![j]!).mul(d[j]!));
      }
    }
    return rnExp(rnNeg(rnScalePow2(rnFromRational(q, p), -1)), p);
  }

  /**
   * `2 * sigma**2` as an element of `RealField(self.precision)`, exactly as
   * Sage forms it from `self._sigma` (`discrete_gaussian_lattice.py:288-292`).
   */
  private _twoSigmaSq(): RealNumberMP {
    const sig = rnFromRational(exactFromDouble(this._sigma as number), this.precision);
    return rnScalePow2(rnMul(sig, sig, this.precision), 1);
  }

  /**
   * An approximation of `sum_{x in B} exp(-|x|_2^2 / (2 sigma^2))`, i.e. the
   * normalisation factor such that the sum over all probabilities is 1 for
   * `B`, via Poisson summation.
   *
   * Reference: `discrete_gaussian_lattice.py:190-355`.  The theta series is
   * delegated to PARI's `qfrep`, exactly as Sage does
   * (`Q.__pari__().qfrep(tau * sigma, 0)` at `py:340`, `Q.__pari__().qfrep(bound, 0)`
   * at `py:351`), and the sum is accumulated in `RealField(prec)`.
   *
   * @param tau - all vectors `v` with `|v|_2^2 <= tau sigma` are enumerated;
   *   if none is provided, enumerate vectors with increasing norm until the
   *   sum converges to the given precision
   * @param prec - passed to `compute_precision`; when omitted the precision of
   *   `sigma` is used, i.e. this sampler's precision
   *
   * @see Deviation: this used to return a `number`.  It now returns a
   * {@link RealNumberMP}, i.e. the `RealField(prec)` element SageMath returns —
   * a `number` cannot hold the 91 significant bits of the `prec = 100` doctest.
   * The class defines `valueOf()`, so `nf / x`, `1 / nf` and `nf < y` are
   * unchanged; use `.toNumber()` where a `number` is required, `.toString()`
   * for SageMath's printed form and `.round()` for SageMath's `round()`.
   */
  _normalisation_factor_zz(tau?: number, prec?: number): RealNumberMP {
    if (!this.is_spherical) {
      // discrete_gaussian_lattice.py:295-313.  Upstream flags this branch as
      // "only a poor approximation placeholder" and emits the warning below;
      // we reproduce the placeholder (and the warning) verbatim.  Note that
      // upstream returns *before* `prec` is looked at, so the non-spherical
      // result always lives in `self._RR`, i.e. at `this.precision` bits.
      console.warn(
        'Note: `_normalisation_factor_zz` has not been properly implemented for non-spherical distributions.'
      );
      const p = this.precision;
      const basis = lllReduce(this.basisExact);
      const solved = ratSolveLeft(basis, this._c);
      const base = solved.map((v) => new Rational(v.round(), 1n));
      // BOUND is the largest integer such that |coords| <= 10^4 (py:305-309).
      let BOUND = Math.max(1, Math.floor((Math.ceil(10 ** (4 / this.n)) - 1) / 2));
      BOUND = Math.min(BOUND, 10);
      let total = rnZero(p);
      const u = new Array<number>(this.n).fill(-BOUND);
      const rec = (i: number): void => {
        if (i === this.n) {
          const coords = u.map((v, j) => new Rational(BigInt(v), 1n).add(base[j]!));
          // py:313: `(vector(u) + base) * self.B` — the *original* basis, not
          // the LLL-reduced one used to pick `base`.
          const pt = coords.map((_, j) =>
            coords.reduce((acc, ck, k) => acc.add(ck.mul(this.basisExact[k]![j]!)), Rational.zero())
          );
          total = rnAdd(total, this.f_R(pt), p);
          return;
        }
        for (let v = -BOUND; v <= BOUND; v++) {
          u[i] = v;
          rec(i + 1);
        }
      };
      rec(0);
      return total;
    }

    if (this.rank !== this.degree) {
      throw new NotImplementedError('basis must be a square matrix');
    }
    if (!this.isIntegral) {
      throw new NotImplementedError('lattice must be integral');
    }
    if (!this._c_in_lattice_and_lattice_trivial) {
      throw new NotImplementedError('center must be at zero and basis must be trivial');
    }

    // py:322-324: prec = compute_precision(prec, self._sigma).  `self._sigma`
    // lives in RealField(this.precision), so `sigma.precision()` is
    // `this.precision`.
    const P =
      prec === undefined || prec === null
        ? this.precision
        : DiscreteGaussianDistributionLatticeSampler.compute_precision(prec, this._sigma);
    // Working precision for the symbolic sub-expressions Sage only evaluates
    // at the very end (`norm_factor`, and the argument of `exp`).
    const wp = P + 64;

    const sigmaExact = exactFromDouble(this._sigma as number);
    const sigmaGt1 = sigmaExact.cmp(Rational.one()) > 0;
    const sp = this.precision;
    const twoSigmaSq = this._twoSigmaSq();
    const piW = rnPi(wp);

    // f_or_hat (py:285-293)
    //
    // sigma <= 1: `R(exp(-x / (2 * sigma**2)))`.  `x` is a Rational and
    // `2*sigma**2` a `RealField(self._RR.prec())` element, so the quotient AND
    // the `exp` are evaluated at *sigma's* precision and only then widened to
    // `R` — a `prec` above 53 buys no accuracy here.  Verified against a live
    // SageMath 10.3: `RealField(100)(exp(-QQ(1)/(2*RealField(53)(0.5)**2)))`
    // is `0.13533528323661270231781372786`, i.e. `exp(-2)` at 53 bits.
    //
    // sigma > 1: `R(exp(-pi**2 * (2 * sigma**2) * x))` is a *symbolic*
    // expression (`pi` is `sage.symbolic.constants.pi`).
    //
    // @see Deviation: pynac rewrites `exp(-y)` as `cosh(y) - sinh(y)`, so
    // SageMath evaluates that symbolic exponential with a catastrophic
    // cancellation of about `2y/log 2` bits: for `sigma = 1.1, x = 1` a live
    // SageMath 10.3 returns `RealField(53)(exp(-2.42*pi^2)) == 0` and
    // `RealField(100)(...) == 4.2375843245100...e-11` where the true value is
    // `4.2375843253587...e-11`.  Upstream's own source comment ("Fun fact: ...
    // RR(1 + 100 * exp(-5.0 * pi^2)) == 0", py:280-284) records the same wart.
    // We evaluate the exponential correctly instead.  This changes nothing for
    // any doctest value — every `sigma > 1` doctest has increments far below
    // the ulp of the sum either way — but for `1 < sigma < ~1.3` our sum has
    // the terms SageMath silently drops.
    const f_or_hat = (x: Rational): RealNumberMP => {
      if (!sigmaGt1) {
        const arg = rnNeg(rnDiv(rnFromRational(x, sp), twoSigmaSq, sp));
        return rnSetPrec(rnExp(arg, sp), P);
      }
      const arg = rnMul(
        rnMul(piW, piW, wp),
        rnMul(rnSetPrec(twoSigmaSq, wp), rnFromRational(x, wp), wp),
        wp
      );
      return rnSetPrec(rnExp(rnNeg(arg), wp), P);
    };

    // py:326-332: `norm_factor = (sigma * sqrt(2 * pi))**self.n / det`.
    //
    // `sigma * sqrt(2*pi)` is symbolic, and pynac pulls the numeric factor out
    // of the power: `(1000.0 * sqrt(2*pi))^8` normalises to
    // `1.60000000000000e25 * pi^4`, i.e. `sigma^n` is rounded to *sigma's*
    // precision before the (exactly evaluated) `(2 pi)^(n/2)` is applied.
    // 10^24 does not fit in 53 bits, which is exactly why the doctest reads
    // 1558545456544038969634991553 and not the mathematically correct
    // 1558545456544038995783045323.  Verified against a live SageMath for
    // (n, sigma) = (2, 3.0), (3, 2.5), (5, 1.5) and (8, 1000).
    let det = 1n;
    let norm_factor = rnOne(wp);
    if (sigmaGt1) {
      const d = ratDet(this.basisExact);
      if (!d.isInteger()) {
        throw new NotImplementedError('lattice must be integral');
      }
      det = d.numerator;
      // `sigma**n` at sigma's own precision, rounded once from the exact value
      // (which is what MPFR's `pow` does).
      const sigmaPowN = rnFromRational(ratPowInt(sigmaExact, this.n), sp);
      const twoPi = rnScalePow2(piW, 1);
      const twoPiHalfN =
        this.n % 2 === 0
          ? rnPowInt(twoPi, this.n / 2, wp)
          : rnMul(rnPowInt(twoPi, (this.n - 1) / 2, wp), rnSqrt(twoPi, wp), wp);
      norm_factor = rnDiv(
        rnMul(rnSetPrec(sigmaPowN, wp), twoPiHalfN, wp),
        rnFromBigInt(det, wp),
        wp
      );
    }
    const detN = det ** BigInt(this.n);
    const finish = (res: RealNumberMP): RealNumberMP =>
      rnSetPrec(rnMul(norm_factor, rnSetPrec(res, wp), wp), P);

    // py:334-336: qfrep computes the theta series of a quadratic form, which
    // is *half* the generating function of number of vectors with given norm.
    const Q = this.Q;
    if (tau !== undefined && tau !== null) {
      // py:338: freq = Q.__pari__().qfrep(tau * sigma, 0)
      const bound = rnMul(
        rnFromRational(exactFromDouble(tau), this.precision),
        rnFromRational(sigmaExact, this.precision),
        this.precision
      ).floor();
      const freq = qfrep(Q, bound);
      let res = rnOne(P);
      for (let x = 0; x < freq.length; x++) {
        const inc = rnMul(
          rnFromBigInt(2n * freq[x]!, P),
          f_or_hat(new Rational(BigInt(x + 1), detN)),
          P
        );
        res = rnAdd(res, inc, P);
      }
      return finish(res);
    }

    // py:345-355
    let res = rnOne(P);
    let bound = 0;
    for (;;) {
      bound += 1;
      const cnt = qfrep(Q, bound)[bound - 1]!;
      const inc = rnMul(rnFromBigInt(2n * cnt, P), f_or_hat(new Rational(BigInt(bound), detN)), P);
      const next = rnAdd(res, inc, P);
      if (cnt > 0n && rnEq(res, next)) {
        return finish(res);
      }
      res = next;
      if (bound > 100000) {
        // Sage would loop forever here; refuse rather than hang.
        throw new NotImplementedError(
          'SAGE_NOT_IMPLEMENTED: _normalisation_factor_zz did not converge; pass tau'
        );
      }
    }
  }

  /**
   * Randomly round to the lattice coset `Z + v` with Gaussian parameter `r`.
   *
   * Reference: `discrete_gaussian_lattice.py:391-407`.
   */
  _randomise(v: number[]): bigint[] {
    return v.map((vi) =>
      new DiscreteGaussianDistributionIntegerSampler({ sigma: this.r!, c: vi }).sample()
    );
  }

  /**
   * Precompute samples from `B^{-1} D_1` for {@link _call_non_spherical}.
   *
   * Reference: `discrete_gaussian_lattice.py:874-892`.
   */
  add_offline_samples(cnt = 1): void {
    if (this.is_spherical) {
      throw new ValueError('offline samples are only used by the non-spherical sampler');
    }
    const p = current_randstate().python_random();
    for (let i = 0; i < cnt; i++) {
      const coord = Array.from({ length: this.n }, () => p.normalvariate(0, 1));
      const M = this.B2_B_inv!;
      // vector(coord) * B2_B_inv
      this.offline_samples.push(
        M[0]!.map((_, j) => coord.reduce((acc, ci, i2) => acc + ci * M[i2]![j]!, 0))
      );
    }
  }

  /**
   * Sample from the discrete Gaussian distribution over the lattice, exactly.
   *
   * Dispatches exactly as Sage's `__call__`
   * (`discrete_gaussian_lattice.py:662-696`).
   */
  sampleExact(): Rational[] {
    if (!this.is_spherical) {
      return this._call_non_spherical();
    }
    if (this._c_in_lattice_and_lattice_trivial) {
      return this._call_simple();
    }
    return this._call();
  }

  /**
   * `_call_simple`: `c` is in the lattice and `B^* == 1`.
   *
   * Reference: `discrete_gaussian_lattice.py:822-840`.
   */
  _call_simple(): Rational[] {
    const w = this.D!.map((d) => new Rational(d.sample(), 1n));
    return this._c.map((ci, j) =>
      w.reduce((acc, wi, i) => acc.add(wi.mul(this.basisExact[i]![j]!)), ci)
    );
  }

  /**
   * `_call`: the general GPV loop.
   *
   * Reference: `discrete_gaussian_lattice.py:842-872`.
   */
  _call(): Rational[] {
    let c = [...this._c];
    let v: Rational[] = new Array(this.degree).fill(Rational.zero());
    const sigma = this._sigma as number;

    for (let i = this.rank - 1; i >= 0; i--) {
      const bStarI = this.gs.bStar[i]!;
      const bStarNormSq = this.gs.bStarNormsSq[i]!;
      const ci = ratDot(c, bStarI).div(bStarNormSq);
      const sigmaI = sigma / Math.sqrt(bStarNormSq.toNumber());
      const Di = new DiscreteGaussianDistributionIntegerSampler({
        sigma: sigmaI,
        c: ci.toNumber(),
        tau: this.tau,
        // Sage passes this explicitly (py:869).
        algorithm: 'uniform+online',
      });
      const z = Di.sample();

      const zRat = new Rational(z, 1n);
      const bi = this.basisExact[i]!;
      c = c.map((x, j) => x.sub(zRat.mul(bi[j]!)));
      v = v.map((x, j) => x.add(zRat.mul(bi[j]!)));
    }

    return v;
  }

  /**
   * `_call_non_spherical`: Peikert's sampler.
   *
   * Reference: `discrete_gaussian_lattice.py:894-915`.
   */
  _call_non_spherical(): Rational[] {
    if (this.offline_samples.length === 0) {
      this.add_offline_samples();
    }
    const off = this.offline_samples.pop()!;
    const vec = this._c_mul_B_inv!.map((x, i) => x.toNumber() - off[i]!);
    const z = this._randomise(vec);
    return this.basisExact[0]!.map((_, j) =>
      z.reduce(
        (acc, zi, i) => acc.add(new Rational(zi, 1n).mul(this.basisExact[i]![j]!)),
        Rational.zero()
      )
    );
  }

  /**
   * Sample from the discrete Gaussian distribution over the lattice.
   *
   * @returns A sample as a bigint vector
   * @throws ValueError if the sample is not integral; use {@link sampleExact}
   */
  sample(): bigint[] {
    if (!this.isIntegral) {
      throw new ValueError(
        'lattice basis is not integral; use sampleExact() for exact rational samples'
      );
    }
    const v = this.sampleExact();
    if (!v.every((x) => x.isInteger())) {
      // Possible with a non-integral centre on an integral lattice.
      throw new ValueError('sample is not integral; use sampleExact() for exact rational samples');
    }
    return v.map((x) => x.numerator);
  }

  /**
   * Alias for sample() - makes the sampler callable like a function.
   *
   * @returns A sample from the discrete Gaussian distribution
   */
  call(): bigint[] {
    return this.sample();
  }

  /**
   * Generate multiple samples.
   *
   * @param n - Number of samples to generate
   * @returns Array of samples
   */
  samples(n: number): bigint[][] {
    const result: bigint[][] = [];
    for (let i = 0; i < n; i++) {
      result.push(this.sample());
    }
    return result;
  }

  /**
   * Generate multiple exact samples.
   *
   * @param n - Number of samples to generate
   * @returns Array of exact samples
   */
  samplesExact(n: number): Rational[][] {
    const result: Rational[][] = [];
    for (let i = 0; i < n; i++) {
      result.push(this.sampleExact());
    }
    return result;
  }

  /**
   * Compute the smoothing parameter eta_epsilon(L) for this lattice.
   *
   * For practical purposes, we use the approximation:
   *   eta_epsilon(L) >= ||b_n*|| * sqrt(ln(2n(1 + 1/epsilon)) / pi)
   *
   * @see Deviation: not present in Sage; a convenience estimate.
   *
   * @param epsilon - The statistical parameter (default: 2^{-n})
   * @returns An estimate of the smoothing parameter
   */
  smoothingParameter(epsilon?: number): number {
    const n = this.rank;
    const eps = epsilon ?? 2 ** -n;
    const maxBStarNorm = Math.sqrt(Math.max(...this.gs.bStarNormsSq.map((x) => x.toNumber())));
    return maxBStarNorm * Math.sqrt(Math.log(2 * n * (1 + 1 / eps)) / Math.PI);
  }

  /**
   * Check if sigma is above the smoothing parameter.
   *
   * @param epsilon - The statistical parameter (default: 2^{-n})
   * @returns True if sigma >= eta_epsilon(L); always false for a covariance
   *   matrix, for which the scalar comparison is meaningless.
   */
  isAboveSmoothingParameter(epsilon?: number): boolean {
    if (!this.is_spherical) return false;
    return (this._sigma as number) >= this.smoothingParameter(epsilon);
  }

  /**
   * String representation.
   *
   * Reference: `discrete_gaussian_lattice.py:794-820`.
   */
  repr(): string {
    const sigma_str = this.is_spherical
      ? `σ = ${rrRepr(this._sigma as number)}`
      : `Σ =\n${matrixRepr((this._sigma as number[][]).map((r) => r.map(rrRepr)))}`;
    const cStr = `(${this._c.map((x) => x.toString()).join(', ')})`;
    const bStr = matrixRepr(this.basisExact.map((r) => r.map((x) => x.toString())));
    return `Discrete Gaussian sampler with Gaussian parameter ${sigma_str}, c=${cStr} over lattice with basis\n\n${bStr}`;
  }

  /**
   * String representation for console output.
   */
  toString(): string {
    return this.repr();
  }
}

/**
 * Factory function to create a discrete Gaussian lattice sampler.
 *
 * @param basis - The lattice basis
 * @param sigma - Gaussian parameter (scalar or covariance matrix)
 * @param c - Center (default: origin)
 * @param tau - Tail cutoff (default: 6)
 * @returns A discrete Gaussian lattice sampler
 */
export function DiscreteGaussianLattice(
  basis: number[][] | bigint[][] | Rational[][],
  sigma: number | number[][] | Rational[][],
  c?: number[] | bigint[] | Rational[],
  tau: IntegerLike = 6n
): DiscreteGaussianDistributionLatticeSampler {
  return new DiscreteGaussianDistributionLatticeSampler(basis, { sigma, c, tau });
}

/**
 * Sample a short vector from a lattice using discrete Gaussian sampling.
 *
 * This is a convenience function that creates a sampler and returns one sample.
 *
 * @param basis - The lattice basis
 * @param sigma - Standard deviation
 * @returns A short lattice vector
 */
export function sampleShortVector(
  basis: number[][] | bigint[][] | Rational[][],
  sigma: number
): bigint[] {
  const D = new DiscreteGaussianDistributionLatticeSampler(basis, { sigma });
  return D.sample();
}

/**
 * Preimage sampling: given a target t, sample a lattice vector close to t.
 *
 * This samples from D_{L, sigma, c} where c is chosen such that the sample
 * is statistically close to the nearest lattice vector to t.
 *
 * @param basis - The lattice basis
 * @param sigma - Standard deviation
 * @param target - The target vector
 * @returns A lattice vector close to target
 */
export function samplePreimage(
  basis: number[][] | bigint[][] | Rational[][],
  sigma: number,
  target: number[] | bigint[] | Rational[]
): bigint[] {
  const D = new DiscreteGaussianDistributionLatticeSampler(basis, { sigma, c: target });
  return D.sample();
}

/**
 * Discrete Gaussian sampler for polynomials (coefficient-wise).
 *
 * Samples a polynomial where each coefficient is drawn from a 1D discrete Gaussian.
 *
 * @see Deviation: SageMath's `DiscreteGaussianDistributionPolynomialSampler`
 * lives in `sage.crypto.lwe` (ported at `crypto/lwe.ts`), takes
 * `(P, n, sigma)` and returns an element of the polynomial ring `P`. This
 * class is an extra convenience wrapper living in the lattice module with
 * signature `(n, options)` that returns a coefficient array. It has no
 * SageMath counterpart; prefer `crypto/lwe.ts`'s class for parity.
 */
export class DiscreteGaussianDistributionPolynomialSampler {
  /**
   * The integer sampler used for each coefficient.
   */
  public readonly integerSampler: DiscreteGaussianDistributionIntegerSampler;

  /**
   * Number of coefficients.
   */
  public readonly n: number;

  /**
   * Construct a polynomial sampler.
   *
   * @param n - Number of coefficients
   * @param options - Options for the underlying integer sampler
   */
  constructor(n: number, options: DiscreteGaussianOptions) {
    if (n <= 0) {
      throw new ValueError(`n must be > 0, got ${n}`);
    }
    this.n = n;
    this.integerSampler = new DiscreteGaussianDistributionIntegerSampler(options);
  }

  /**
   * Sample a polynomial (as an array of coefficients).
   *
   * @returns Array of coefficients [a_0, a_1, ..., a_{n-1}]
   */
  sample(): bigint[] {
    return this.integerSampler.samples(this.n);
  }

  /**
   * Alias for sample().
   */
  call(): bigint[] {
    return this.sample();
  }

  /**
   * Generate multiple polynomial samples.
   *
   * @param count - Number of polynomials to sample
   * @returns Array of polynomial coefficient arrays
   */
  samples(count: number): bigint[][] {
    const result: bigint[][] = [];
    for (let i = 0; i < count; i++) {
      result.push(this.sample());
    }
    return result;
  }

  /**
   * Get the sigma parameter.
   */
  get sigma(): number {
    return this.integerSampler.sigma;
  }

  /**
   * Get the c parameter.
   */
  get c(): number {
    return this.integerSampler.c;
  }

  /**
   * String representation.
   */
  repr(): string {
    return `DiscreteGaussianDistributionPolynomialSampler(n=${this.n}, sigma=${this.sigma}, c=${this.c})`;
  }

  toString(): string {
    return this.repr();
  }
}

/**
 * Factory function to create a discrete Gaussian polynomial sampler.
 *
 * @param n - Number of coefficients
 * @param sigma - Standard deviation
 * @param c - Center (default: 0)
 * @param tau - Tail cutoff (default: 6)
 * @returns A discrete Gaussian polynomial sampler
 */
export function DiscreteGaussianPolynomial(
  n: number,
  sigma: number,
  c: IntegerLike = 0n,
  tau: IntegerLike = 6n
): DiscreteGaussianDistributionPolynomialSampler {
  return new DiscreteGaussianDistributionPolynomialSampler(n, { sigma, c, tau });
}
