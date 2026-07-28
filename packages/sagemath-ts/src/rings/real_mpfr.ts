/**
 * @module sage/rings/real_mpfr
 * @description Arbitrary precision floating point real numbers using GNU MPFR
 *
 * Port of: sage/rings/real_mpfr.pyx
 * Reference: reference/sage/src/sage/rings/real_mpfr.pyx
 *
 * DEVIATION from SageMath:
 * - Uses JavaScript's native number type (IEEE 754 double precision) instead of MPFR
 * - Precision parameter is accepted but has no effect beyond 53 bits (JavaScript limit)
 * - Some special functions (Bessel, erf) use approximations or throw for complex cases
 * - See DEVIATIONS.md for full details
 */

import { algebraic_dependency as arith_algebraic_dependency } from '../arith/misc.js';
import { ValueError } from '../errors.js';
import type { ComplexField } from './complex_mpfr.js';

// Mathematical constants with high precision (as much as JavaScript can handle)
const PI = Math.PI;
const E = Math.E;
const LN2 = Math.LN2;
const LN10 = Math.LN10;
const EULER_CONSTANT = 0.5772156649015329; // Euler-Mascheroni constant
const CATALAN_CONSTANT = 0.915965594177219; // Catalan's constant

/**
 * MPFR's exponent range on a 64-bit platform (`mpfr_get_exp_min()` /
 * `mpfr_get_exp_max()`).  Used by {@link RealNumber.fp_rank}.
 */
const MPFR_EXP_MIN = 1n - (1n << 62n);
const MPFR_EXP_MAX = (1n << 62n) - 1n;

/** log(Gamma(1/2)) = log(sqrt(pi)) */
const LN_GAMMA_HALF = 0.5723649429247001;

/** 2/sqrt(pi), the leading coefficient of the Taylor series of erf at 0. */
const TWO_OVER_SQRT_PI = 1.1283791670955126;

/**
 * Regularised lower incomplete gamma function `P(1/2, y)` for `y >= 0`,
 * evaluated with the ascending series
 * `P(a,y) = e^{-y} y^a / Gamma(a) * sum_{n>=0} y^n / (a (a+1) ... (a+n))`.
 *
 * Only accurate (and only used) for `y < a + 1 = 3/2`.
 */
function gammp_half(y: number): number {
  if (y <= 0) return 0;
  const a = 0.5;
  let ap = a;
  let del = 1 / a;
  let sum = del;
  for (let n = 0; n < 1000; n++) {
    ap += 1;
    del *= y / ap;
    sum += del;
    if (Math.abs(del) < Math.abs(sum) * 1e-18) break;
  }
  return sum * Math.exp(-y + a * Math.log(y) - LN_GAMMA_HALF);
}

/**
 * Regularised upper incomplete gamma function `Q(1/2, y)` for `y > 0`,
 * evaluated with the continued fraction
 * `Q(a,y) = e^{-y} y^a / Gamma(a) * (1/(y+1-a - 1*(1-a)/(y+3-a - 2*(2-a)/(...))))`
 * using the modified Lentz algorithm.
 *
 * Only accurate (and only used) for `y >= a + 1 = 3/2`.
 */
function gammq_half(y: number): number {
  const a = 0.5;
  const TINY = 1e-300;
  let b = y + 1 - a;
  let c = 1 / TINY;
  let d = 1 / b;
  let h = d;
  for (let i = 1; i < 1000; i++) {
    const an = -i * (i - a);
    b += 2;
    d = an * d + b;
    if (Math.abs(d) < TINY) d = TINY;
    c = b + an / c;
    if (Math.abs(c) < TINY) c = TINY;
    d = 1 / d;
    const del = d * c;
    h *= del;
    if (Math.abs(del - 1) < 1e-18) break;
  }
  return Math.exp(-y + a * Math.log(y) - LN_GAMMA_HALF) * h;
}

/**
 * Riemann zeta on `s >= 0`, `s != 1`, via the alternating (eta) series
 * accelerated by Borwein's algorithm 2 (P. Borwein, *An efficient algorithm
 * for the Riemann zeta function*, 1991):
 *
 * ```
 * d_k   = n * sum_{i=0}^{k} (n+i-1)! 4^i / ((n-i)! (2i)!)
 * eta(s) ~ (1 / d_n) * sum_{k=0}^{n-1} (-1)^k (d_k - d_n) / (k+1)^s
 * zeta(s) = eta(s) / (1 - 2^{1-s})
 * ```
 *
 * The error is bounded by `(3/(3+sqrt(8))^n) / |1 - 2^{1-s}|` for `Re(s) >= 1/2`,
 * so `n = 32` is far more than enough for double precision.  Unlike the naive
 * alternating series it also converges throughout the critical strip, which is
 * what `mpfr_zeta` covers.
 */
function borwein_zeta(s: number): number {
  const n = 32;
  // d_k, computed by the standard recurrence on the summand
  //   t_i = (n+i-1)! 4^i / ((n-i)! (2i)!),  t_0 = 1,
  //   t_{i+1} = t_i * 4 * (n+i) * (n-i) / ((2i+1)(2i+2)).
  const d: number[] = new Array(n + 1);
  let t = 1;
  let acc = 1;
  d[0] = n * acc;
  for (let i = 0; i < n; i++) {
    t = (t * 4 * (n + i) * (n - i)) / ((2 * i + 1) * (2 * i + 2));
    acc += t;
    d[i + 1] = n * acc;
  }

  const dn = d[n]!;
  let sum = 0;
  for (let k = 0; k < n; k++) {
    const sign = k % 2 === 0 ? 1 : -1;
    sum += (sign * (d[k]! - dn)) / (k + 1) ** s;
  }
  const eta = -sum / dn;
  return eta / (1 - 2 ** (1 - s));
}

// ---------------------------------------------------------------------------
// Exact rational helpers, used by simplest_rational()/nearby_rational().
// A rational is a normalised [numerator, denominator] pair with denominator > 0
// and gcd(numerator, denominator) = 1.
// ---------------------------------------------------------------------------

type Rat = [bigint, bigint];

function bigAbs(a: bigint): bigint {
  return a < 0n ? -a : a;
}

function bigGcd(a: bigint, b: bigint): bigint {
  let x = bigAbs(a);
  let y = bigAbs(b);
  while (y !== 0n) {
    const t = x % y;
    x = y;
    y = t;
  }
  return x;
}

function ratMake(num: bigint, den: bigint): Rat {
  if (den === 0n) {
    throw new ValueError('rational with zero denominator');
  }
  if (den < 0n) {
    num = -num;
    den = -den;
  }
  const g = bigGcd(num, den);
  if (g > 1n) {
    num /= g;
    den /= g;
  }
  return [num, den];
}

/** The rational `m * 2^e` for an integer exponent `e` of either sign. */
function ratPow2(m: bigint, e: bigint): Rat {
  return e >= 0n ? ratMake(m << e, 1n) : ratMake(m, 1n << -e);
}

function ratAdd(a: Rat, b: Rat): Rat {
  return ratMake(a[0] * b[1] + b[0] * a[1], a[1] * b[1]);
}

function ratSub(a: Rat, b: Rat): Rat {
  return ratMake(a[0] * b[1] - b[0] * a[1], a[1] * b[1]);
}

function ratNeg(a: Rat): Rat {
  return [-a[0], a[1]];
}

function ratInv(a: Rat): Rat {
  return ratMake(a[1], a[0]);
}

/** Compare a and b; returns -1, 0 or 1. */
function ratCmp(a: Rat, b: Rat): number {
  const l = a[0] * b[1];
  const r = b[0] * a[1];
  return l < r ? -1 : l > r ? 1 : 0;
}

/** Compare a with the integer n. */
function ratCmpInt(a: Rat, n: bigint): number {
  const r = n * a[1];
  return a[0] < r ? -1 : a[0] > r ? 1 : 0;
}

function bigFloorDiv(a: bigint, b: bigint): bigint {
  // b > 0 assumed
  const q = a / b;
  return a % b !== 0n && a < 0n ? q - 1n : q;
}

function ratFloor(a: Rat): bigint {
  return bigFloorDiv(a[0], a[1]);
}

function ratCeil(a: Rat): bigint {
  return -bigFloorDiv(-a[0], a[1]);
}

/**
 * Return the simplest rational between `low` and `high`.  May return `low` or
 * `high` unless `lowOpen`/`highOpen` are set.  Both endpoints must be
 * nonnegative and `high > low`.
 *
 * Direct port of `sage/rings/real_mpfi.pyx:_simplest_rational_exact`.
 */
function simplestRationalExact(low: Rat, high: Rat, lowOpen: boolean, highOpen: boolean): Rat {
  if (ratCmpInt(low, 1n) < 0) {
    if (low[0] === 0n) {
      if (lowOpen) {
        if (ratCmpInt(high, 1n) > 0) {
          return [1n, 1n];
        }
        const invHigh = ratInv(high);
        if (highOpen) {
          return ratInv([ratFloor(invHigh) + 1n, 1n]);
        }
        return ratInv([ratCeil(invHigh), 1n]);
      }
      return [0n, 1n];
    }

    if (ratCmpInt(high, 1n) > 0) {
      return [1n, 1n];
    }

    const r = simplestRationalExact(ratInv(high), ratInv(low), highOpen, lowOpen);
    return ratInv(r);
  }

  const fl = ratFloor(low);
  const shifted = simplestRationalExact(
    ratSub(low, [fl, 1n]),
    ratSub(high, [fl, 1n]),
    lowOpen,
    highOpen
  );
  return ratAdd([fl, 1n], shifted);
}

/**
 * The simplest rational in the interval `[low, high]` (open at either end as
 * requested), without the sign/zero preconditions of
 * {@link simplestRationalExact}.
 *
 * Port of `sage/rings/real_mpfi.pyx:RealIntervalFieldElement.simplest_rational`
 * (its exact branch: the floating point fast path returns the same value).
 */
function simplestRationalInInterval(low: Rat, high: Rat, lowOpen: boolean, highOpen: boolean): Rat {
  if (ratCmp(low, high) === 0) {
    if (lowOpen || highOpen) {
      throw new ValueError('simplest_rational() on open, empty interval');
    }
    return low;
  }
  // The interval contains zero.
  if (ratCmpInt(low, 0n) <= 0 && ratCmpInt(high, 0n) >= 0) {
    return [0n, 1n];
  }
  if (ratCmpInt(high, 0n) < 0) {
    return ratNeg(simplestRationalInInterval(ratNeg(high), ratNeg(low), highOpen, lowOpen));
  }
  return simplestRationalExact(low, high, lowOpen, highOpen);
}

/**
 * Rounding modes for MPFR operations
 */
export enum RoundingMode {
  /** Round to nearest (ties go to even) */
  RNDN = 0,
  /** Round towards zero */
  RNDZ = 1,
  /** Round towards minus infinity */
  RNDD = 2,
  /** Round towards plus infinity */
  RNDU = 3,
  /** Round away from zero */
  RNDA = 4,
  /** Faithful rounding */
  RNDF = 5,
}

/**
 * Return the minimum MPFR precision.
 * @see Reference: sage/rings/real_mpfr.pyx:mpfr_prec_min
 */
export function mpfr_prec_min(): number {
  return 1;
}

/**
 * Return the maximum MPFR precision.
 * @see Reference: sage/rings/real_mpfr.pyx:mpfr_prec_max
 */
export function mpfr_prec_max(): bigint {
  // 64-bit value
  return 9223372036854775551n;
}

/**
 * An approximation to the field of real numbers using floating point
 * numbers with any specified precision.
 * @see Reference: sage/rings/real_mpfr.pyx:RealField_class
 */
export class RealField {
  private readonly _prec: number;
  private readonly _sci_not: boolean;
  private readonly _rnd: RoundingMode;

  constructor(prec: number = 53, sci_not: boolean = false, rnd: RoundingMode = RoundingMode.RNDN) {
    if (prec < 1 || prec > Number(mpfr_prec_max())) {
      throw new Error(`prec (=${prec}) must be >= 1 and <= ${mpfr_prec_max()}`);
    }
    this._prec = prec;
    this._sci_not = sci_not;
    this._rnd = rnd;
  }

  /**
   * Return the precision of this field.
   * @see Reference: sage/rings/real_mpfr.pyx:precision
   */
  precision(): number {
    return this._prec;
  }

  /**
   * Return the characteristic of the field (which is 0).
   * @see Reference: sage/rings/real_mpfr.pyx:characteristic
   */
  characteristic(): bigint {
    return 0n;
  }

  /**
   * Return the rounding mode.
   * @see Reference: sage/rings/real_mpfr.pyx:rounding_mode
   */
  rounding_mode(): RoundingMode {
    return this._rnd;
  }

  /**
   * Create an element in this field.
   * @see Reference: sage/rings/real_mpfr.pyx:_element_constructor_
   */
  __call__(x: number | bigint | string): RealNumber {
    return new RealNumber(this, x);
  }

  /**
   * Return the complex field with the same precision.
   * @see Reference: sage/rings/real_mpfr.pyx:complex_field
   */
  complex_field(): ComplexField {
    // Dynamic import to avoid circular dependency
    // eslint-disable-next-line @typescript-eslint/no-require-imports
    const { ComplexField: CF } = require('./complex_mpfr.js');
    return new CF(this._prec);
  }

  /**
   * Return pi to the precision of this field.
   * @see Reference: sage/rings/real_mpfr.pyx:pi
   */
  pi(): RealNumber {
    return new RealNumber(this, PI);
  }

  /**
   * Return Euler's constant to the precision of this field.
   * @see Reference: sage/rings/real_mpfr.pyx:euler_constant
   */
  euler_constant(): RealNumber {
    return new RealNumber(this, EULER_CONSTANT);
  }

  /**
   * Return Catalan's constant to the precision of this field.
   * @see Reference: sage/rings/real_mpfr.pyx:catalan_constant
   */
  catalan_constant(): RealNumber {
    return new RealNumber(this, CATALAN_CONSTANT);
  }

  /**
   * Return log(2) to the precision of this field.
   * @see Reference: sage/rings/real_mpfr.pyx:log2
   */
  log2(): RealNumber {
    return new RealNumber(this, LN2);
  }

  /**
   * Return a random element in [min, max].
   * @see Reference: sage/rings/real_mpfr.pyx:random_element
   */
  random_element(min: number = -1, max: number = 1): RealNumber {
    const value = Math.random() * (max - min) + min;
    return new RealNumber(this, value);
  }

  /**
   * Return n! as an element of this field.
   * @see Reference: sage/rings/real_mpfr.pyx:factorial
   */
  factorial(n: number): RealNumber {
    if (n < 0) {
      throw new Error('n must be nonnegative');
    }
    if (n === 0 || n === 1) {
      return new RealNumber(this, 1);
    }
    // Use Stirling's approximation for large n, otherwise compute directly
    if (n <= 170) {
      // Direct computation (170! is about the max for IEEE 754)
      let result = 1;
      for (let i = 2; i <= n; i++) {
        result *= i;
      }
      return new RealNumber(this, result);
    }
    // Stirling's approximation for large n: n! ~ sqrt(2*pi*n) * (n/e)^n
    const approx = Math.sqrt(2 * PI * n) * (n / E) ** n * (1 + 1 / (12 * n));
    return new RealNumber(this, approx);
  }

  /**
   * Return a primitive n-th root of unity.
   * In the real field, only n=1 (returns 1) and n=2 (returns -1) are valid.
   * @see Reference: sage/rings/real_mpfr.pyx:zeta
   */
  zeta(n: number = 2): RealNumber {
    if (n === 1) {
      return new RealNumber(this, 1);
    }
    if (n === 2) {
      return new RealNumber(this, -1);
    }
    throw new Error(`No ${n}th root of unity in self`);
  }

  toString(): string {
    if (this._rnd !== RoundingMode.RNDN) {
      return `Real Field with ${this._prec} bits of precision and rounding ${RoundingMode[this._rnd]}`;
    }
    return `Real Field with ${this._prec} bits of precision`;
  }
}

/**
 * An arbitrary precision floating point real number.
 * @see Reference: sage/rings/real_mpfr.pyx:RealNumber
 */
export class RealNumber {
  private readonly _parent: RealField;
  private readonly _value: number; // Placeholder - would use MPFR in real implementation

  constructor(parent: RealField, value: number | bigint | string = 0) {
    this._parent = parent;
    if (typeof value === 'bigint') {
      this._value = Number(value);
    } else if (typeof value === 'string') {
      this._value = Number.parseFloat(value);
    } else {
      this._value = value;
    }
  }

  /**
   * Return the precision of this number.
   * @see Reference: sage/rings/real_mpfr.pyx:precision
   */
  precision(): number {
    return this._parent.precision();
  }

  /**
   * Return the real part (which is self).
   * @see Reference: sage/rings/real_mpfr.pyx:real
   */
  real(): RealNumber {
    return this;
  }

  /**
   * Return the imaginary part (which is zero).
   * @see Reference: sage/rings/real_mpfr.pyx:imag
   */
  imag(): RealNumber {
    return new RealNumber(this._parent, 0);
  }

  /**
   * Return the sign of this number.
   * @see Reference: sage/rings/real_mpfr.pyx:sign
   */
  sign(): number {
    if (this._value < 0) return -1;
    if (this._value > 0) return 1;
    return 0;
  }

  /**
   * Return the absolute value.
   * @see Reference: sage/rings/real_mpfr.pyx:__abs__
   */
  abs(): RealNumber {
    return new RealNumber(this._parent, Math.abs(this._value));
  }

  /**
   * Return the floor of this number.
   * @see Reference: sage/rings/real_mpfr.pyx:floor
   */
  floor(): bigint {
    return BigInt(Math.floor(this._value));
  }

  /**
   * Return the ceiling of this number.
   * @see Reference: sage/rings/real_mpfr.pyx:ceil
   */
  ceil(): bigint {
    return BigInt(Math.ceil(this._value));
  }

  /**
   * Return ``self`` rounded to the nearest representable integer, rounding
   * halfway cases **away from zero** (this is ``mpfr_round``, not IEEE
   * round-half-to-even and not JavaScript's ``Math.round`` which rounds
   * halfway cases towards `+\infty`).
   *
   * The rounding mode of the parent field does not affect the result.
   *
   * ```
   * sage: RR(0.49).round()   ->  0
   * sage: RR(0.5).round()    ->  1
   * sage: RR(-0.49).round()  ->  0
   * sage: RR(-0.5).round()   -> -1
   * ```
   *
   * @see Reference: sage/rings/real_mpfr.pyx:round
   */
  round(): bigint {
    const x = this._value;
    if (!Number.isFinite(x)) {
      throw new ValueError('Cannot convert NaN or infinity to Integer');
    }
    // mpfr_round: ties away from zero.
    const rounded = Math.round(Math.abs(x));
    return x < 0 ? -BigInt(rounded) : BigInt(rounded);
  }

  /**
   * Return the truncation towards zero.
   * @see Reference: sage/rings/real_mpfr.pyx:trunc
   */
  trunc(): bigint {
    return BigInt(Math.trunc(this._value));
  }

  /**
   * Return the fractional part.
   * @see Reference: sage/rings/real_mpfr.pyx:frac
   */
  frac(): RealNumber {
    return new RealNumber(this._parent, this._value - Math.trunc(this._value));
  }

  /**
   * Return the square root.
   * For negative numbers, returns NaN (in SageMath this returns a complex number).
   * @see Reference: sage/rings/real_mpfr.pyx:sqrt
   */
  sqrt(): RealNumber {
    if (this._value < 0) {
      // In SageMath, this would return a complex number
      // We return NaN to indicate the result is not a real number
      return new RealNumber(this._parent, Number.NaN);
    }
    return new RealNumber(this._parent, Math.sqrt(this._value));
  }

  /**
   * Return the cube root.
   * @see Reference: sage/rings/real_mpfr.pyx:cube_root
   */
  cube_root(): RealNumber {
    return new RealNumber(this._parent, Math.cbrt(this._value));
  }

  /**
   * Return the n-th root.
   * @see Reference: sage/rings/real_mpfr.pyx:nth_root
   */
  nth_root(n: number): RealNumber {
    if (n <= 0) {
      throw new Error('n must be positive');
    }
    if (this._value === 0) {
      return new RealNumber(this._parent, 0);
    }
    if (this._value < 0) {
      if (n % 2 === 0) {
        throw new Error('taking an even root of a negative number');
      }
      // For odd roots of negative numbers
      return new RealNumber(this._parent, -((-this._value) ** (1 / n)));
    }
    return new RealNumber(this._parent, this._value ** (1 / n));
  }

  /**
   * Return the natural logarithm.
   * For negative numbers, returns NaN (in SageMath this returns a complex number).
   * @see Reference: sage/rings/real_mpfr.pyx:log
   */
  log(base?: number): RealNumber {
    if (Number.isNaN(this._value)) {
      return this;
    }
    if (this._value < 0) {
      // In SageMath, this would return a complex number
      return new RealNumber(this._parent, Number.NaN);
    }
    if (base === undefined || base === null) {
      return new RealNumber(this._parent, Math.log(this._value));
    }
    if (base === 10) {
      return this.log10();
    }
    if (base === 2) {
      return this.log2();
    }
    // log_base(x) = ln(x) / ln(base)
    return new RealNumber(this._parent, Math.log(this._value) / Math.log(base));
  }

  /**
   * Return log base 2.
   * @see Reference: sage/rings/real_mpfr.pyx:log2
   */
  log2(): RealNumber {
    if (this._value < 0) {
      return new RealNumber(this._parent, Number.NaN);
    }
    return new RealNumber(this._parent, Math.log2(this._value));
  }

  /**
   * Return log base 10.
   * @see Reference: sage/rings/real_mpfr.pyx:log10
   */
  log10(): RealNumber {
    if (this._value < 0) {
      return new RealNumber(this._parent, Number.NaN);
    }
    return new RealNumber(this._parent, Math.log10(this._value));
  }

  /**
   * Return log(1 + self).
   * More accurate than log(1 + x) for small x.
   * @see Reference: sage/rings/real_mpfr.pyx:log1p
   */
  log1p(): RealNumber {
    if (this._value < -1) {
      return new RealNumber(this._parent, Number.NaN);
    }
    return new RealNumber(this._parent, Math.log1p(this._value));
  }

  /**
   * Return the exponential of this number (e^self).
   * @see Reference: sage/rings/real_mpfr.pyx:exp
   */
  exp(): RealNumber {
    return new RealNumber(this._parent, Math.exp(this._value));
  }

  /**
   * Return 2^self.
   * @see Reference: sage/rings/real_mpfr.pyx:exp2
   */
  exp2(): RealNumber {
    return new RealNumber(this._parent, 2 ** this._value);
  }

  /**
   * Return 10^self.
   * @see Reference: sage/rings/real_mpfr.pyx:exp10
   */
  exp10(): RealNumber {
    return new RealNumber(this._parent, 10 ** this._value);
  }

  /**
   * Return exp(self) - 1.
   * More accurate than exp(x) - 1 for small x.
   * @see Reference: sage/rings/real_mpfr.pyx:expm1
   */
  expm1(): RealNumber {
    return new RealNumber(this._parent, Math.expm1(this._value));
  }

  /**
   * Return the cosine.
   * @see Reference: sage/rings/real_mpfr.pyx:cos
   */
  cos(): RealNumber {
    return new RealNumber(this._parent, Math.cos(this._value));
  }

  /**
   * Return the sine.
   * @see Reference: sage/rings/real_mpfr.pyx:sin
   */
  sin(): RealNumber {
    return new RealNumber(this._parent, Math.sin(this._value));
  }

  /**
   * Return the tangent.
   * @see Reference: sage/rings/real_mpfr.pyx:tan
   */
  tan(): RealNumber {
    return new RealNumber(this._parent, Math.tan(this._value));
  }

  /**
   * Return (sin(self), cos(self)).
   * @see Reference: sage/rings/real_mpfr.pyx:sincos
   */
  sincos(): [RealNumber, RealNumber] {
    return [
      new RealNumber(this._parent, Math.sin(this._value)),
      new RealNumber(this._parent, Math.cos(this._value)),
    ];
  }

  /**
   * Return the arccosine.
   * @see Reference: sage/rings/real_mpfr.pyx:arccos
   */
  arccos(): RealNumber {
    return new RealNumber(this._parent, Math.acos(this._value));
  }

  /**
   * Return the arcsine.
   * @see Reference: sage/rings/real_mpfr.pyx:arcsin
   */
  arcsin(): RealNumber {
    return new RealNumber(this._parent, Math.asin(this._value));
  }

  /**
   * Return the arctangent.
   * @see Reference: sage/rings/real_mpfr.pyx:arctan
   */
  arctan(): RealNumber {
    return new RealNumber(this._parent, Math.atan(this._value));
  }

  /**
   * Return the hyperbolic cosine.
   * @see Reference: sage/rings/real_mpfr.pyx:cosh
   */
  cosh(): RealNumber {
    return new RealNumber(this._parent, Math.cosh(this._value));
  }

  /**
   * Return the hyperbolic sine.
   * @see Reference: sage/rings/real_mpfr.pyx:sinh
   */
  sinh(): RealNumber {
    return new RealNumber(this._parent, Math.sinh(this._value));
  }

  /**
   * Return the hyperbolic tangent.
   * @see Reference: sage/rings/real_mpfr.pyx:tanh
   */
  tanh(): RealNumber {
    return new RealNumber(this._parent, Math.tanh(this._value));
  }

  /**
   * Return the inverse hyperbolic cosine.
   * @see Reference: sage/rings/real_mpfr.pyx:arccosh
   */
  arccosh(): RealNumber {
    return new RealNumber(this._parent, Math.acosh(this._value));
  }

  /**
   * Return the inverse hyperbolic sine.
   * @see Reference: sage/rings/real_mpfr.pyx:arcsinh
   */
  arcsinh(): RealNumber {
    return new RealNumber(this._parent, Math.asinh(this._value));
  }

  /**
   * Return the inverse hyperbolic tangent.
   * @see Reference: sage/rings/real_mpfr.pyx:arctanh
   */
  arctanh(): RealNumber {
    return new RealNumber(this._parent, Math.atanh(this._value));
  }

  /**
   * Return the hyperbolic cotangent.
   * coth(x) = 1/tanh(x) = cosh(x)/sinh(x)
   * @see Reference: sage/rings/real_mpfr.pyx:coth
   */
  coth(): RealNumber {
    return new RealNumber(this._parent, 1 / Math.tanh(this._value));
  }

  /**
   * Return the cotangent.
   * cot(x) = 1/tan(x) = cos(x)/sin(x)
   * @see Reference: sage/rings/real_mpfr.pyx:cot
   */
  cot(): RealNumber {
    return new RealNumber(this._parent, 1 / Math.tan(this._value));
  }

  /**
   * Return the secant.
   * sec(x) = 1/cos(x)
   * @see Reference: sage/rings/real_mpfr.pyx:sec
   */
  sec(): RealNumber {
    return new RealNumber(this._parent, 1 / Math.cos(this._value));
  }

  /**
   * Return the cosecant.
   * csc(x) = 1/sin(x)
   * @see Reference: sage/rings/real_mpfr.pyx:csc
   */
  csc(): RealNumber {
    return new RealNumber(this._parent, 1 / Math.sin(this._value));
  }

  /**
   * Return the hyperbolic secant.
   * sech(x) = 1/cosh(x)
   * @see Reference: sage/rings/real_mpfr.pyx:sech
   */
  sech(): RealNumber {
    return new RealNumber(this._parent, 1 / Math.cosh(this._value));
  }

  /**
   * Return the hyperbolic cosecant.
   * csch(x) = 1/sinh(x)
   * @see Reference: sage/rings/real_mpfr.pyx:csch
   */
  csch(): RealNumber {
    return new RealNumber(this._parent, 1 / Math.sinh(this._value));
  }

  /**
   * Return the arithmetic-geometric mean of self and other.
   * The AGM is the common limit of the sequences u_n and v_n where:
   * u_0 = self, v_0 = other, u_{n+1} = (u_n + v_n)/2, v_{n+1} = sqrt(u_n * v_n)
   * @see Reference: sage/rings/real_mpfr.pyx:agm
   */
  agm(other: RealNumber | number): RealNumber {
    let a = this._value;
    let b = typeof other === 'number' ? other : other.toNumber();

    if (a < 0 || b < 0) {
      return new RealNumber(this._parent, Number.NaN);
    }

    // Iterate until convergence
    const eps = Number.EPSILON;
    while (Math.abs(a - b) > eps * Math.max(Math.abs(a), Math.abs(b))) {
      const an = (a + b) / 2;
      const bn = Math.sqrt(a * b);
      a = an;
      b = bn;
    }
    return new RealNumber(this._parent, (a + b) / 2);
  }

  /**
   * Return the error function of this number.
   * erf(x) = (2/sqrt(pi)) * integral from 0 to x of exp(-t^2) dt
   *
   * SageMath delegates to `mpfr_erf`, which is correctly rounded.  We obtain
   * full double precision through the regularised incomplete gamma function:
   * `erf(x) = P(1/2, x^2)` and `erfc(x) = Q(1/2, x^2)` for `x >= 0`, using the
   * series for `P` and the continued fraction for `Q`.
   *
   * @see Reference: sage/rings/real_mpfr.pyx:erf
   */
  erf(): RealNumber {
    const x = this._value;

    // Handle special cases
    if (x === 0) {
      return new RealNumber(this._parent, 0);
    }
    if (Number.isNaN(x)) {
      return new RealNumber(this._parent, Number.NaN);
    }
    if (!Number.isFinite(x)) {
      return new RealNumber(this._parent, x > 0 ? 1 : -1);
    }

    const ax = Math.abs(x);
    const x2 = ax * ax;
    if (x2 === 0) {
      // x^2 underflowed; erf(x) = 2x/sqrt(pi) + O(x^3) there.
      return new RealNumber(this._parent, TWO_OVER_SQRT_PI * x);
    }
    // P(a, y) is best computed by its series for y < a + 1 = 3/2.
    const value = x2 < 1.5 ? gammp_half(x2) : 1 - gammq_half(x2);
    return new RealNumber(this._parent, x > 0 ? value : -value);
  }

  /**
   * Return the complementary error function `erfc(x) = 1 - erf(x)`.
   *
   * Computing this as `1 - erf(x)` loses all significance once `erf(x)` is
   * within a rounding error of 1; SageMath's `mpfr_erfc` does not.  We use
   * the continued fraction for `Q(1/2, x^2)` for `x^2 >= 3/2` instead, which
   * keeps full relative accuracy: `R(6).erfc()` = 2.15e-17,
   * `R(10).erfc()` = 2.09e-45.
   *
   * @see Reference: sage/rings/real_mpfr.pyx:erfc
   */
  erfc(): RealNumber {
    const x = this._value;
    if (x === 0) {
      return new RealNumber(this._parent, 1);
    }
    if (Number.isNaN(x)) {
      return new RealNumber(this._parent, Number.NaN);
    }
    if (!Number.isFinite(x)) {
      return new RealNumber(this._parent, x > 0 ? 0 : 2);
    }

    const ax = Math.abs(x);
    const x2 = ax * ax;
    // erfc(|x|) = Q(1/2, x^2); erfc(-|x|) = 2 - erfc(|x|).
    const positiveTail = x2 < 1.5 ? 1 - gammp_half(x2) : gammq_half(x2);
    return new RealNumber(this._parent, x > 0 ? positiveTail : 2 - positiveTail);
  }

  /**
   * Return the Bessel function J_0(self).
   * @see Reference: sage/rings/real_mpfr.pyx:j0
   */
  j0(): RealNumber {
    const x = this._value;

    // Special case: J_0(0) = 1
    if (x === 0) {
      return new RealNumber(this._parent, 1);
    }

    // Use polynomial approximations for different ranges
    if (Math.abs(x) < 8) {
      // Use series expansion for small x
      const x2 = x * x;
      const num =
        57568490574.0 +
        x2 *
          (-13362590354.0 +
            x2 * (651619640.7 + x2 * (-11214424.18 + x2 * (77392.33017 + x2 * -184.9052456))));
      const den =
        57568490411.0 +
        x2 * (1029532985.0 + x2 * (9494680.718 + x2 * (59272.64853 + x2 * (267.8532712 + x2))));
      return new RealNumber(this._parent, num / den);
    } else {
      // Asymptotic expansion for large x
      const ax = Math.abs(x);
      const z = 8 / ax;
      const z2 = z * z;
      const xx = ax - 0.785398163397448;
      const p =
        1 +
        z2 *
          (-0.001098628627 +
            z2 * (0.00002734510407 + z2 * (-0.000002073370639 + z2 * 2.093887211e-7)));
      const q =
        -0.01562499995 +
        z2 *
          (0.0001430488765 +
            z2 * (-0.000006911147651 + z2 * (7.621095161e-7 + z2 * -9.34945152e-8)));
      return new RealNumber(
        this._parent,
        Math.sqrt(0.636619772367581 / ax) * (Math.cos(xx) * p - z * Math.sin(xx) * q)
      );
    }
  }

  /**
   * Return the Bessel function J_1(self).
   * @see Reference: sage/rings/real_mpfr.pyx:j1
   */
  j1(): RealNumber {
    const x = this._value;
    if (Math.abs(x) < 8) {
      // Series expansion for small x
      const x2 = x * x;
      const num =
        x *
        (72362614232.0 +
          x2 *
            (-7895059235.0 +
              x2 * (242396853.1 + x2 * (-2972611.439 + x2 * (15704.4826 + x2 * -30.16036606)))));
      const den =
        144725228442.0 +
        x2 * (2300535178.0 + x2 * (18583304.74 + x2 * (99447.43394 + x2 * (376.9991397 + x2))));
      return new RealNumber(this._parent, num / den);
    } else {
      // Asymptotic expansion for large x
      const ax = Math.abs(x);
      const z = 8 / ax;
      const z2 = z * z;
      const xx = ax - 2.356194490192345;
      const p =
        1 +
        z2 *
          (0.00183105 + z2 * (-0.00003516396496 + z2 * (0.000002457520174 + z2 * -2.40337019e-7)));
      const q =
        0.04687499995 +
        z2 *
          (-0.0002002690873 + z2 * (0.000008449199096 + z2 * (-8.8228987e-7 + z2 * 1.05787412e-7)));
      const result = Math.sqrt(0.636619772367581 / ax) * (Math.cos(xx) * p - z * Math.sin(xx) * q);
      return new RealNumber(this._parent, x < 0 ? -result : result);
    }
  }

  /**
   * Return the Bessel function J_n(self).
   * Uses the recurrence relation for integer orders.
   * @see Reference: sage/rings/real_mpfr.pyx:jn
   */
  jn(n: number): RealNumber {
    const x = this._value;
    if (n === 0) return this.j0();
    if (n === 1) return this.j1();

    const sign = n < 0 && n % 2 !== 0 ? -1 : 1;
    n = Math.abs(n);

    if (x === 0) {
      return new RealNumber(this._parent, n === 0 ? 1 : 0);
    }

    // Miller's algorithm using downward recurrence.
    // The recurrence is seeded with J_{m+1} = 0, J_m = 1 (arbitrary scale) and
    // renormalised at the end via 1 = J_0 + 2*(J_2 + J_4 + ...).
    const ax = Math.abs(x);
    if (n > Math.floor(ax)) {
      // Numerical Recipes uses ACC = 40, which only targets single precision.
      // mpfr_jn is correctly rounded, so we start the downward recurrence
      // further out (the extra terms cost nothing) to reach double precision.
      const ACC = 160;
      const BIGNO = 1e10;
      const BIGNI = 1e-10;
      const tox = 2 / ax;
      let bjp = 0; // J_{m+1} = 0 seeds the downward recurrence
      let bj = 1; // J_m = 1 (unnormalised)
      let sum = 0;
      let ans = 0;
      const m = 2 * Math.floor((n + Math.floor(Math.sqrt(ACC * n))) / 2);

      for (let j = m; j > 0; j--) {
        const bjm = j * tox * bj - bjp;
        bjp = bj;
        bj = bjm;
        if (Math.abs(bj) > BIGNO) {
          bj *= BIGNI;
          bjp *= BIGNI;
          ans *= BIGNI;
          sum *= BIGNI;
        }
        // After the update `bj` holds J_{j-1}; accumulate the even-order terms
        // of the normalisation 1 = J_0 + 2*(J_2 + J_4 + ...).
        if (j % 2 !== 0) sum += bj;
        if (j === n) ans = bjp;
      }
      sum = 2 * sum - bj;
      ans /= sum;
      return new RealNumber(this._parent, x < 0 && n % 2 !== 0 ? -sign * ans : sign * ans);
    } else {
      // Upward recurrence, seeded at |x| (J_1 is odd, so the sign of x must be
      // applied once, at the end, and not through the seed).
      const tox = 2 / ax;
      const abs = new RealNumber(this._parent, ax);
      let bjm = abs.j0().toNumber();
      let bj = abs.j1().toNumber();
      for (let j = 1; j < n; j++) {
        const bjp = j * tox * bj - bjm;
        bjm = bj;
        bj = bjp;
      }
      return new RealNumber(this._parent, x < 0 && n % 2 !== 0 ? -sign * bj : sign * bj);
    }
  }

  /**
   * Return the Bessel function Y_0(self).
   * @see Reference: sage/rings/real_mpfr.pyx:y0
   */
  y0(): RealNumber {
    const x = this._value;
    if (x <= 0) {
      return new RealNumber(this._parent, x === 0 ? Number.NEGATIVE_INFINITY : Number.NaN);
    }
    if (x < 8) {
      const x2 = x * x;
      const num =
        -2957821389.0 +
        x2 *
          (7062834065.0 +
            x2 * (-512359803.6 + x2 * (10879881.29 + x2 * (-86327.92757 + x2 * 228.4622733))));
      const den =
        40076544269.0 +
        x2 * (745249964.8 + x2 * (7189466.438 + x2 * (47447.2647 + x2 * (226.1030244 + x2))));
      return new RealNumber(
        this._parent,
        num / den + 0.636619772367581 * this.j0().toNumber() * Math.log(x)
      );
    } else {
      const z = 8 / x;
      const z2 = z * z;
      const xx = x - 0.785398163397448;
      const p =
        1 +
        z2 *
          (-0.001098628627 +
            z2 * (0.00002734510407 + z2 * (-0.000002073370639 + z2 * 2.093887211e-7)));
      const q =
        -0.01562499995 +
        z2 *
          (0.0001430488765 +
            z2 * (-0.000006911147651 + z2 * (7.621095161e-7 + z2 * -9.34945152e-8)));
      return new RealNumber(
        this._parent,
        Math.sqrt(0.636619772367581 / x) * (Math.sin(xx) * p + z * Math.cos(xx) * q)
      );
    }
  }

  /**
   * Return the Bessel function Y_1(self).
   * @see Reference: sage/rings/real_mpfr.pyx:y1
   */
  y1(): RealNumber {
    const x = this._value;
    if (x <= 0) {
      return new RealNumber(this._parent, x === 0 ? Number.NEGATIVE_INFINITY : Number.NaN);
    }
    if (x < 8) {
      const x2 = x * x;
      const num =
        x *
        (-0.4900604943e13 +
          x2 *
            (0.127527439e13 +
              x2 *
                (-0.5153438139e11 +
                  x2 * (0.7349264551e9 + x2 * (-0.4237922726e7 + x2 * 0.8511937935e4)))));
      const den =
        0.249958057e14 +
        x2 *
          (0.4244419664e12 +
            x2 *
              (0.3733650367e10 +
                x2 * (0.2245904002e8 + x2 * (0.102042605e6 + x2 * (0.3549632885e3 + x2)))));
      return new RealNumber(
        this._parent,
        num / den + 0.636619772367581 * (this.j1().toNumber() * Math.log(x) - 1 / x)
      );
    } else {
      const z = 8 / x;
      const z2 = z * z;
      const xx = x - 2.356194490192345;
      const p =
        1 +
        z2 *
          (0.00183105 + z2 * (-0.00003516396496 + z2 * (0.000002457520174 + z2 * -2.40337019e-7)));
      const q =
        0.04687499995 +
        z2 *
          (-0.0002002690873 + z2 * (0.000008449199096 + z2 * (-8.8228987e-7 + z2 * 1.05787412e-7)));
      return new RealNumber(
        this._parent,
        Math.sqrt(0.636619772367581 / x) * (Math.sin(xx) * p + z * Math.cos(xx) * q)
      );
    }
  }

  /**
   * Return the Bessel function Y_n(self).
   * Uses the recurrence relation for integer orders.
   * @see Reference: sage/rings/real_mpfr.pyx:yn
   */
  yn(n: number): RealNumber {
    const x = this._value;
    if (x <= 0) {
      return new RealNumber(this._parent, x === 0 ? Number.NEGATIVE_INFINITY : Number.NaN);
    }
    if (n === 0) return this.y0();
    if (n === 1) return this.y1();

    const sign = n < 0 && n % 2 !== 0 ? -1 : 1;
    n = Math.abs(n);

    // Upward recurrence
    const tox = 2 / x;
    let bym = this.y0().toNumber();
    let by = this.y1().toNumber();
    for (let j = 1; j < n; j++) {
      const byp = j * tox * by - bym;
      bym = by;
      by = byp;
    }
    return new RealNumber(this._parent, sign * by);
  }

  /**
   * Return the gamma function of this number.
   * Gamma(n) = (n-1)! for positive integers.
   * @see Reference: sage/rings/real_mpfr.pyx:gamma
   */
  gamma(): RealNumber {
    const x = this._value;

    // Use Lanczos approximation
    if (x <= 0 && Number.isInteger(x)) {
      // Gamma has poles at non-positive integers
      return new RealNumber(this._parent, x === 0 ? Number.POSITIVE_INFINITY : Number.NaN);
    }

    // Lanczos coefficients
    const g = 7;
    const c = [
      0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313,
      -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6,
      1.5056327351493116e-7,
    ];

    if (x < 0.5) {
      // Reflection formula: Gamma(x) * Gamma(1-x) = pi / sin(pi*x)
      return new RealNumber(
        this._parent,
        PI / (Math.sin(PI * x) * new RealNumber(this._parent, 1 - x).gamma().toNumber())
      );
    }

    const xm1 = x - 1;
    let sum = c[0]!;
    for (let i = 1; i < g + 2; i++) {
      sum += c[i]! / (xm1 + i);
    }

    const t = xm1 + g + 0.5;
    const result = Math.sqrt(2 * PI) * t ** (xm1 + 0.5) * Math.exp(-t) * sum;
    return new RealNumber(this._parent, result);
  }

  /**
   * Return log(gamma(self)).
   * More stable than log(gamma(x)) for large x.
   * @see Reference: sage/rings/real_mpfr.pyx:log_gamma
   */
  log_gamma(): RealNumber {
    const x = this._value;

    if (x <= 0) {
      // For negative values, the principal branch involves complex numbers
      return new RealNumber(this._parent, Number.NaN);
    }

    // log(gamma(1)) = log(1) = 0
    if (x === 1) {
      return new RealNumber(this._parent, 0);
    }

    // Use Stirling's approximation for large x
    if (x > 10) {
      // log(Gamma(x)) ≈ (x-0.5)*ln(x) - x + 0.5*ln(2*pi) + 1/(12*x) - ...
      const lnx = Math.log(x);
      return new RealNumber(
        this._parent,
        (x - 0.5) * lnx -
          x +
          0.5 * Math.log(2 * PI) +
          1 / (12 * x) -
          1 / (360 * x * x * x) +
          1 / (1260 * x * x * x * x * x)
      );
    }

    // For smaller values, use log(gamma(x))
    return new RealNumber(this._parent, Math.log(this.gamma().toNumber()));
  }

  /**
   * Return the Riemann zeta function at this number.
   * @see Reference: sage/rings/real_mpfr.pyx:zeta
   */
  zeta(): RealNumber {
    const s = this._value;

    if (Number.isNaN(s)) {
      return new RealNumber(this._parent, Number.NaN);
    }
    if (s === 1) {
      return new RealNumber(this._parent, Number.POSITIVE_INFINITY);
    }

    // For negative even integers, zeta is 0 (trivial zeros).
    if (s < 0 && Number.isInteger(s) && s % 2 === 0) {
      return new RealNumber(this._parent, 0);
    }

    // mpfr_zeta is defined on the whole real line.  Only reflect for s < 0;
    // reflecting on the critical strip 0 < s < 1 maps it onto itself and
    // never terminates.
    if (s < 0) {
      // zeta(s) = 2^s * pi^(s-1) * sin(pi*s/2) * gamma(1-s) * zeta(1-s)
      const oneMinusS = 1 - s;
      const zetaOneMinusS = new RealNumber(this._parent, oneMinusS).zeta().toNumber();
      const gammaOneMinusS = new RealNumber(this._parent, oneMinusS).gamma().toNumber();
      const result =
        2 ** s * PI ** (s - 1) * Math.sin((PI * s) / 2) * gammaOneMinusS * zetaOneMinusS;
      return new RealNumber(this._parent, result);
    }

    return new RealNumber(this._parent, borwein_zeta(s));
  }

  /**
   * Return the exponential integral Ei(self).
   * Ei(x) = -integral from -x to infinity of exp(-t)/t dt
   * @see Reference: sage/rings/real_mpfr.pyx:eint
   */
  eint(): RealNumber {
    const x = this._value;

    if (x === 0) {
      return new RealNumber(this._parent, Number.NEGATIVE_INFINITY);
    }

    // For small |x|, use series expansion
    if (Math.abs(x) < 40) {
      // Series expansion: Ei(x) = gamma + ln|x| + sum_{k=1}^inf x^k / (k * k!)
      let sum = EULER_CONSTANT + Math.log(Math.abs(x));
      let term = x;
      let factorial = 1;
      for (let k = 1; k < 100; k++) {
        sum += term / (k * factorial);
        factorial *= k + 1;
        term *= x;
        if (Math.abs(term / (k * factorial)) < Number.EPSILON * Math.abs(sum)) break;
      }
      return new RealNumber(this._parent, sum);
    }

    // Asymptotic expansion for large |x|
    // Ei(x) ~ exp(x)/x * (1 + 1!/x + 2!/x^2 + 3!/x^3 + ...)
    let sum = 0;
    let term = 1;
    for (let k = 0; k < 50; k++) {
      const oldSum = sum;
      sum += term;
      if (Math.abs(term) < Number.EPSILON * Math.abs(sum) || sum === oldSum) break;
      term *= (k + 1) / x;
    }
    return new RealNumber(this._parent, (Math.exp(x) / x) * sum);
  }

  /**
   * Return the exact rational representation as [numerator, denominator].
   * For IEEE 754 doubles, returns the exact fraction represented.
   * @see Reference: sage/rings/real_mpfr.pyx:exact_rational
   */
  exact_rational(): [bigint, bigint] {
    const x = this._value;

    if (!Number.isFinite(x)) {
      throw new Error('Cannot convert NaN or infinity to rational');
    }

    if (x === 0) {
      return [0n, 1n];
    }

    // Get the exact representation as sign * mantissa * 2^exponent
    const [sign, mantissa, exponent] = this.sign_mantissa_exponent();

    if (exponent >= 0n) {
      // Result is an integer: mantissa * 2^exponent
      return [BigInt(sign) * mantissa * (1n << exponent), 1n];
    } else {
      // Result is a fraction: mantissa / 2^(-exponent)
      const denom = 1n << -exponent;
      return [BigInt(sign) * mantissa, denom];
    }
  }

  /**
   * Return the simplest rational which is equal to `self` (in the Sage sense),
   * i.e. the simplest rational which rounds back to `self` in this field.
   * Given `a/b` and `c/d` in lowest terms, the former is simpler if `b < d`, or
   * if `b = d` and `|a| < |c|`.
   *
   * In round-to-nearest mode, this is the simplest rational in the interval
   * bounded by the two `(prec+1)`-bit neighbours of `self`; the interval is
   * closed when the mantissa of `self` is even and open when it is odd (because
   * MPFR breaks ties towards an even mantissa).
   *
   * ```
   * sage: RR(pi).simplest_rational()       -> 245850922/78256779
   * sage: RR(2).sqrt().simplest_rational() -> 131836323/93222358
   * sage: RR(1234).simplest_rational()     -> 1234
   * ```
   *
   * The result always rounds back to `self`, unlike a bare continued-fraction
   * truncation.
   *
   * @see Reference: sage/rings/real_mpfr.pyx:simplest_rational
   */
  simplest_rational(): [bigint, bigint] {
    const x = this._value;

    if (!Number.isFinite(x)) {
      throw new ValueError('cannot convert NaN or infinity to rational number');
    }
    if (x === 0) {
      return [0n, 1n];
    }
    if (x < 0) {
      const [num, den] = new RealNumber(this._parent, -x).simplest_rational();
      return [-num, den];
    }

    // Decompose |self| = m * 2^e with m having exactly 53 bits (the storage
    // precision of a JavaScript double).
    let [, m, e] = this.sign_mantissa_exponent();
    if (m < 0n) m = -m;
    const bits = BigInt(m.toString(2).length);
    if (bits < 53n) {
      const shift = 53n - bits;
      m <<= shift;
      e -= shift;
    }

    // The two neighbours at precision 54 straddle self by half an ulp; at a
    // power of two the gap below is half as wide.
    const powerOfTwo = m === 1n << 52n;
    const low: Rat = powerOfTwo ? ratPow2((1n << 54n) - 1n, e - 2n) : ratPow2(2n * m - 1n, e - 1n);
    const high: Rat = ratPow2(2n * m + 1n, e - 1n);

    // Ties round to an even mantissa, so the neighbours round back to self
    // exactly when the mantissa of self is even; otherwise the interval is open.
    const odd = (m & 1n) === 1n;
    return simplestRationalInInterval(low, high, odd, odd);
  }

  /**
   * Find a rational near to `self`.  Exactly one of `max_error` or
   * `max_denominator` must be specified.
   *
   * If `max_error` is given, return the simplest rational in the range
   * `[self - max_error, self + max_error]`.  If `max_denominator` is given,
   * return the rational closest to `self` with denominator at most
   * `max_denominator` (ties broken towards the simpler rational).
   *
   * ```
   * sage: (0.333).nearby_rational(max_error=0.001)             -> 1/3
   * sage: (0.333).nearby_rational(max_error=1)                 -> 0
   * sage: (-0.333).nearby_rational(max_error=0.0001)           -> -257/772
   * sage: RR(1/3 + 1/1000000).nearby_rational(max_denominator=2999999)
   *                                                            -> 777780/2333333
   * sage: RR(pi).nearby_rational(max_denominator=100000)       -> 312689/99532
   * sage: RR(-3.5).nearby_rational(max_denominator=1)          -> -3
   * ```
   *
   * @see Reference: sage/rings/real_mpfr.pyx:nearby_rational
   */
  nearby_rational(max_error?: number, max_denominator?: bigint): [bigint, bigint] {
    const x = this._value;

    if (!Number.isFinite(x)) {
      throw new ValueError('cannot convert NaN or infinity to rational number');
    }

    const hasError = max_error !== undefined && max_error !== null;
    const hasDenom = max_denominator !== undefined && max_denominator !== null;
    if (hasError === hasDenom) {
      throw new ValueError(
        'Must specify exactly one of max_error or max_denominator in nearby_rational()'
      );
    }

    if (hasError) {
      const low = new RealNumber(this._parent, x - max_error).exact_rational();
      const high = new RealNumber(this._parent, x + max_error).exact_rational();
      return simplestRationalInInterval(low, high, false, false);
    }

    const maxDenominator = max_denominator!;
    if (maxDenominator < 1n) {
      throw new ValueError('max_denominator must be at least 1');
    }
    if (x === 0) {
      return [0n, 1n];
    }

    const negative = x < 0;
    let selfR: Rat = new RealNumber(this._parent, x).exact_rational();
    if (selfR[1] <= maxDenominator) {
      return selfR;
    }
    if (negative) {
      selfR = ratNeg(selfR);
    }

    const fl = ratFloor(selfR);
    const target = ratSub(selfR, [fl, 1n]);

    // Stern-Brocot search for the two neighbours of `target` with denominator
    // at most max_denominator, taking many tree steps at a time.
    let a = 0n;
    let b = 1n;
    const c = target[0];
    const d = target[1];
    let e = 1n;
    let f = 1n;

    let lowDone = false;
    let highDone = false;

    while (!lowDone || !highDone) {
      // Move the low side.
      let k = (c * b - d * a) / (d * e - c * f);
      if (b + k * f > maxDenominator) {
        k = (maxDenominator - b) / f;
        lowDone = true;
      }
      if (k === 0n) {
        lowDone = true;
      }
      a = a + k * e;
      b = b + k * f;

      // Move the high side.
      k = (d * e - c * f) / (c * b - d * a);
      if (k * b + f >= maxDenominator) {
        k = (maxDenominator - f) / b;
        highDone = true;
      }
      if (k === 0n) {
        highDone = true;
      }
      e = k * a + e;
      f = k * b + f;
    }

    const low: Rat = ratMake(a, b);
    const high: Rat = ratMake(e, f);
    const d0 = ratSub(target, low);
    const d1 = ratSub(high, target);

    let result: Rat;
    const cmp = ratCmp(d1, d0);
    if (cmp < 0) {
      result = high;
    } else if (cmp > 0) {
      result = low;
    } else if (f < b) {
      result = high;
    } else if (b < f) {
      result = low;
    } else if (e < a) {
      result = high;
    } else {
      result = low;
    }

    result = ratAdd(result, [fl, 1n]);
    return negative ? ratNeg(result) : result;
  }

  /**
   * Check if this is NaN.
   * @see Reference: sage/rings/real_mpfr.pyx:is_NaN
   */
  is_NaN(): boolean {
    return Number.isNaN(this._value);
  }

  /**
   * Check if this is positive infinity.
   * @see Reference: sage/rings/real_mpfr.pyx:is_positive_infinity
   */
  is_positive_infinity(): boolean {
    return this._value === Number.POSITIVE_INFINITY;
  }

  /**
   * Check if this is negative infinity.
   * @see Reference: sage/rings/real_mpfr.pyx:is_negative_infinity
   */
  is_negative_infinity(): boolean {
    return this._value === Number.NEGATIVE_INFINITY;
  }

  /**
   * Check if this is infinite (positive or negative).
   * @see Reference: sage/rings/real_mpfr.pyx:is_infinity
   */
  is_infinity(): boolean {
    return !Number.isFinite(this._value) && !Number.isNaN(this._value);
  }

  /**
   * Check if this is an integer.
   * @see Reference: sage/rings/real_mpfr.pyx:is_integer
   */
  is_integer(): boolean {
    return Number.isInteger(this._value);
  }

  /**
   * Check if this is a square.
   * @see Reference: sage/rings/real_mpfr.pyx:is_square
   */
  is_square(): boolean {
    return this._value >= 0;
  }

  /**
   * Return the conjugate (which is self for real numbers).
   * @see Reference: sage/rings/real_mpfr.pyx:conjugate
   */
  conjugate(): RealNumber {
    return this;
  }

  /**
   * Return the multiplicative order.
   * In the real field, only 1 and -1 have finite multiplicative order;
   * everything else has order `+Infinity` (SageMath returns
   * `sage.rings.infinity.infinity`, it does not raise).
   *
   * ```
   * sage: RR(1).multiplicative_order()   ->  1
   * sage: RR(-1).multiplicative_order()  ->  2
   * sage: RR(3).multiplicative_order()   ->  +Infinity
   * ```
   *
   * @see Reference: sage/rings/real_mpfr.pyx:multiplicative_order
   */
  multiplicative_order(): number {
    if (this._value === 1) {
      return 1;
    }
    if (this._value === -1) {
      return 2;
    }
    return Number.POSITIVE_INFINITY;
  }

  /**
   * Return the next floating point number towards other.
   * @see Reference: sage/rings/real_mpfr.pyx:nexttoward
   */
  nexttoward(other: RealNumber | number): RealNumber {
    const direction = typeof other === 'number' ? other : other.toNumber();

    if (this._value === direction) {
      return this;
    }

    if (direction > this._value) {
      return this.nextabove();
    } else {
      return this.nextbelow();
    }
  }

  /**
   * Return the next floating point number above self.
   * @see Reference: sage/rings/real_mpfr.pyx:nextabove
   */
  nextabove(): RealNumber {
    const x = this._value;

    if (Number.isNaN(x)) {
      return this;
    }
    if (x === Number.POSITIVE_INFINITY) {
      return this;
    }
    if (x === Number.NEGATIVE_INFINITY) {
      return new RealNumber(this._parent, -Number.MAX_VALUE);
    }
    if (x === 0) {
      return new RealNumber(this._parent, Number.MIN_VALUE);
    }

    // Get the bit representation and add/subtract 1
    const buffer = new ArrayBuffer(8);
    const f64 = new Float64Array(buffer);
    const i64 = new BigInt64Array(buffer);

    f64[0] = x;
    const bits = i64[0]!;

    if (x > 0) {
      i64[0] = bits + 1n;
    } else {
      i64[0] = bits - 1n;
    }

    return new RealNumber(this._parent, f64[0]!);
  }

  /**
   * Return the next floating point number below self.
   * @see Reference: sage/rings/real_mpfr.pyx:nextbelow
   */
  nextbelow(): RealNumber {
    const x = this._value;

    if (Number.isNaN(x)) {
      return this;
    }
    if (x === Number.NEGATIVE_INFINITY) {
      return this;
    }
    if (x === Number.POSITIVE_INFINITY) {
      return new RealNumber(this._parent, Number.MAX_VALUE);
    }
    if (x === 0) {
      return new RealNumber(this._parent, -Number.MIN_VALUE);
    }

    // Get the bit representation and add/subtract 1
    const buffer = new ArrayBuffer(8);
    const f64 = new Float64Array(buffer);
    const i64 = new BigInt64Array(buffer);

    f64[0] = x;
    const bits = i64[0]!;

    if (x > 0) {
      i64[0] = bits - 1n;
    } else {
      i64[0] = bits + 1n;
    }

    return new RealNumber(this._parent, f64[0]!);
  }

  /**
   * Return the unit in the last place (ulp).
   * The ulp is the distance to the next representable floating point number.
   * @see Reference: sage/rings/real_mpfr.pyx:ulp
   */
  ulp(): RealNumber {
    const x = this._value;

    if (Number.isNaN(x)) {
      return new RealNumber(this._parent, Number.NaN);
    }
    if (!Number.isFinite(x)) {
      // mpfr_set_inf(x, 1): +infinity for both +oo and -oo
      return new RealNumber(this._parent, Number.POSITIVE_INFINITY);
    }
    if (x === 0) {
      // The ulp of zero is the smallest representable positive number.
      return new RealNumber(this._parent, Number.MIN_VALUE);
    }

    // e = mpfr_get_exp(self) - prec, result = 2^e rounded up.
    // mpfr's exponent convention is self = m * 2^e with 1/2 <= |m| < 1, i.e.
    // e = floor(log2|self|) + 1.  Derived exactly from the bit pattern rather
    // than from a logarithm.
    const e = this._mpfr_exponent() - this._parent.precision();
    const value = 2 ** e;
    // "Round up in case of underflow" (MPFR_RNDU).
    return new RealNumber(this._parent, value === 0 ? Number.MIN_VALUE : value);
  }

  /**
   * Return `abs(self)` divided by `2^b`, where `b` is the precision in bits of
   * `self`.  Equivalently, `abs(self)` multiplied by the ulp of 1.  This is a
   * scale-invariant version of {@link ulp}: it lies in `[u/2, u)` where `u` is
   * `self.ulp()`.
   *
   * ```
   * sage: RR(2^53).epsilon()  ->  1.00000000000000
   * sage: RR(0).epsilon()     ->  0.000000000000000
   * sage: RR.pi().epsilon()   ->  3.48786849800863e-16
   * ```
   *
   * Note that this is **not** the constant `2^(1-prec)`.
   *
   * @see Reference: sage/rings/real_mpfr.pyx:epsilon
   */
  epsilon(): RealNumber {
    const x = this._value;
    if (Number.isNaN(x)) {
      return new RealNumber(this._parent, Number.NaN);
    }
    if (!Number.isFinite(x)) {
      return new RealNumber(this._parent, Number.POSITIVE_INFINITY);
    }
    const prec = this._parent.precision();
    return new RealNumber(this._parent, Math.abs(x) / 2 ** prec);
  }

  /**
   * Return the MPFR exponent of this (finite, nonzero) number: the unique `e`
   * with `self = m * 2^e` and `1/2 <= |m| < 1`.
   */
  private _mpfr_exponent(): number {
    const [, mantissa, exp] = this.sign_mantissa_exponent();
    // self = mantissa * 2^exp, so self = (mantissa / 2^bits) * 2^(exp + bits)
    // where bits = bitlength(mantissa) puts the fraction in [1/2, 1).
    return Number(exp) + mantissa.toString(2).length;
  }

  /**
   * Return (sign, mantissa, exponent) such that self = sign * mantissa * 2^exponent.
   * For IEEE 754 double precision: self = sign * (1.mantissa) * 2^(exponent - 1023)
   * @see Reference: sage/rings/real_mpfr.pyx:sign_mantissa_exponent
   */
  sign_mantissa_exponent(): [number, bigint, bigint] {
    const x = this._value;

    if (Number.isNaN(x)) {
      throw new Error('Cannot decompose NaN');
    }
    if (!Number.isFinite(x)) {
      throw new Error('Cannot decompose infinity');
    }

    if (x === 0) {
      return [1, 0n, 0n];
    }

    const sign = x < 0 ? -1 : 1;
    const absX = Math.abs(x);

    // Get the bit representation
    const buffer = new ArrayBuffer(8);
    const f64 = new Float64Array(buffer);
    const u64 = new BigUint64Array(buffer);

    f64[0] = absX;
    const bits = u64[0]!;

    // IEEE 754 double: 1 sign bit, 11 exponent bits, 52 mantissa bits
    const expBits = Number((bits >> 52n) & 0x7ffn);
    const mantissaBits = bits & 0xfffffffffffffn;

    if (expBits === 0) {
      // Subnormal number: mantissa has no implicit leading 1
      const exp = -1022n - 52n;
      return [sign, mantissaBits, exp];
    } else {
      // Normal number: mantissa has implicit leading 1
      const exp = BigInt(expBits - 1023 - 52);
      const mantissa = (1n << 52n) | mantissaBits;
      return [sign, mantissa, exp];
    }
  }

  /**
   * Return the floating-point rank of this number: if you list the
   * floating-point numbers of this precision in order, and number them
   * starting with `0.0 -> 0`, extending to positive and negative infinity,
   * this is the number corresponding to this floating-point number.
   *
   * SageMath's formula (`real_mpfr.pyx:fp_rank`), with `p` the precision and
   * `EXP_MIN`/`EXP_MAX` the MPFR exponent range:
   *
   * ```
   * fp_rank(0)          = 0
   * fp_rank(m * 2^e)    = m + (e + p - EXP_MIN - 1) * 2^(p-1) + 1
   * fp_rank(infinity)   = (EXP_MAX + 1 - EXP_MIN) * 2^(p-1) + 1
   * fp_rank(-x)         = -fp_rank(x)
   * ```
   *
   * where `m` is the `p`-bit mantissa of `mpfr_get_z_exp`.  This is *not* the
   * raw IEEE bit pattern: `RR(1).fp_rank()` is
   * `20769187434139310514121985316880385`, not `4607182418800017408`.
   *
   * @see Reference: sage/rings/real_mpfr.pyx:fp_rank
   */
  fp_rank(): bigint {
    const x = this._value;

    if (Number.isNaN(x)) {
      throw new ValueError('Cannot compute fp_rank of NaN');
    }

    if (x === 0) {
      return 0n;
    }

    const prec = BigInt(this._parent.precision());
    const half = 1n << (prec - 1n);
    const sgn = x < 0 ? -1n : 1n;

    if (!Number.isFinite(x)) {
      const rank = (MPFR_EXP_MAX + 1n - MPFR_EXP_MIN) * half + 1n;
      return sgn * rank;
    }

    // Normalise to the mpfr_get_z_exp form: |self| = mantissa * 2^exponent
    // with mantissa having exactly `prec` bits.
    let [, mantissa, exponent] = this.sign_mantissa_exponent();
    if (mantissa < 0n) mantissa = -mantissa;
    let bits = BigInt(mantissa.toString(2).length);
    if (bits < prec) {
      const shift = prec - bits;
      mantissa <<= shift;
      exponent -= shift;
    } else if (bits > prec) {
      // The value is stored as a double; report the rank of its rounding to
      // `prec` bits (round to nearest, ties to even), as MPFR would store it.
      const shift = bits - prec;
      const dropped = mantissa & ((1n << shift) - 1n);
      const halfway = 1n << (shift - 1n);
      mantissa >>= shift;
      exponent += shift;
      if (dropped > halfway || (dropped === halfway && (mantissa & 1n) === 1n)) {
        mantissa += 1n;
        if (mantissa === 1n << prec) {
          mantissa >>= 1n;
          exponent += 1n;
        }
      }
      bits = prec;
    }

    let rank = (exponent + prec - MPFR_EXP_MIN - 1n) * half + 1n;
    rank += mantissa;
    return sgn * rank;
  }

  /**
   * Return an algebraic dependency.
   * Uses the LLL algorithm to find a polynomial of degree at most n
   * that this number approximately satisfies.
   *
   * @param n - Maximum degree of the polynomial
   * @returns Coefficients of the polynomial [a_0, a_1, ..., a_n] such that
   *          a_0 + a_1*x + ... + a_n*x^n is approximately satisfied by this number
   *
   * @example
   * ```typescript
   * const R = new RealField();
   * const sqrt2 = R.__call__(Math.sqrt(2));
   * const poly = sqrt2.algebraic_dependency(5);
   * // Returns approximately [-2, 0, 1] representing x^2 - 2
   * ```
   *
   * @see Reference: sage/rings/real_mpfr.pyx:algebraic_dependency
   */
  algebraic_dependency(n: number): bigint[] {
    // Delegate to the fixed arith/misc implementation which uses proper IntegerMatrix LLL
    // Reference: sage/arith/misc.py:algebraic_dependency
    return arith_algebraic_dependency(this._value, BigInt(n));
  }

  /**
   * Alias for algebraic_dependency.
   * @see algebraic_dependency
   */
  algdep(n: number): bigint[] {
    return this.algebraic_dependency(n);
  }

  /**
   * Convert to JavaScript number.
   */
  toNumber(): number {
    return this._value;
  }

  toString(): string {
    return this._value.toString();
  }

  /**
   * Return self raised to the power of exponent.
   * @see Reference: sage/rings/real_mpfr.pyx:__pow__
   */
  pow(exponent: RealNumber | number): RealNumber {
    const exp = typeof exponent === 'number' ? exponent : exponent.toNumber();
    const result = this._value ** exp;

    // If the result is NaN and we had a negative base, the actual result is complex
    if (Number.isNaN(result) && this._value < 0) {
      // In SageMath this would return a complex number
      return new RealNumber(this._parent, Number.NaN);
    }

    return new RealNumber(this._parent, result);
  }

  /**
   * Return self + other.
   */
  add(other: RealNumber | number): RealNumber {
    const otherVal = typeof other === 'number' ? other : other.toNumber();
    return new RealNumber(this._parent, this._value + otherVal);
  }

  /**
   * Return self - other.
   */
  sub(other: RealNumber | number): RealNumber {
    const otherVal = typeof other === 'number' ? other : other.toNumber();
    return new RealNumber(this._parent, this._value - otherVal);
  }

  /**
   * Return self * other.
   */
  mul(other: RealNumber | number): RealNumber {
    const otherVal = typeof other === 'number' ? other : other.toNumber();
    return new RealNumber(this._parent, this._value * otherVal);
  }

  /**
   * Return self / other.
   */
  div(other: RealNumber | number): RealNumber {
    const otherVal = typeof other === 'number' ? other : other.toNumber();
    return new RealNumber(this._parent, this._value / otherVal);
  }

  /**
   * Return -self.
   */
  neg(): RealNumber {
    return new RealNumber(this._parent, -this._value);
  }

  /**
   * Return the parent field.
   */
  parent(): RealField {
    return this._parent;
  }

  /**
   * Compare self to other.
   * Returns -1 if self < other, 0 if self == other, 1 if self > other.
   */
  cmp(other: RealNumber | number): number {
    const otherVal = typeof other === 'number' ? other : other.toNumber();
    if (this._value < otherVal) return -1;
    if (this._value > otherVal) return 1;
    return 0;
  }

  /**
   * Check equality with other.
   */
  equals(other: RealNumber | number): boolean {
    const otherVal = typeof other === 'number' ? other : other.toNumber();
    return this._value === otherVal;
  }
}

/**
 * Create a real field with given precision.
 * @see Reference: sage/rings/real_mpfr.pyx:RealField
 */
export function RR(prec: number = 53): RealField {
  return new RealField(prec);
}

// Default real field with 53 bits of precision
export const RealFieldDefault = new RealField(53);
