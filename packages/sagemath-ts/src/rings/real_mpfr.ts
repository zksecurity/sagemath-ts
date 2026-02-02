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

import { NotImplementedError } from '../errors.js';
import type { ComplexField } from './complex_mpfr.js';

// Mathematical constants with high precision (as much as JavaScript can handle)
const PI = Math.PI;
const E = Math.E;
const LN2 = Math.LN2;
const LN10 = Math.LN10;
const EULER_CONSTANT = 0.5772156649015329; // Euler-Mascheroni constant
const CATALAN_CONSTANT = 0.915965594177219; // Catalan's constant

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
   * Return this number rounded to the nearest integer.
   * @see Reference: sage/rings/real_mpfr.pyx:round
   */
  round(): bigint {
    return BigInt(Math.round(this._value));
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
   * @see Reference: sage/rings/real_mpfr.pyx:erf
   */
  erf(): RealNumber {
    const x = this._value;

    // Handle special cases
    if (x === 0) {
      return new RealNumber(this._parent, 0);
    }

    // Use a Horner form polynomial approximation (Abramowitz and Stegun)
    const t = 1 / (1 + 0.5 * Math.abs(x));

    // Coefficients for the approximation
    const tau =
      t *
      Math.exp(
        -x * x -
          1.26551223 +
          t *
            (1.00002368 +
              t *
                (0.37409196 +
                  t *
                    (0.09678418 +
                      t *
                        (-0.18628806 +
                          t *
                            (0.27886807 +
                              t *
                                (-1.13520398 +
                                  t * (1.48851587 + t * (-0.82215223 + t * 0.17087277))))))))
      );

    const result = x >= 0 ? 1 - tau : tau - 1;
    return new RealNumber(this._parent, result);
  }

  /**
   * Return the complementary error function.
   * erfc(x) = 1 - erf(x)
   * @see Reference: sage/rings/real_mpfr.pyx:erfc
   */
  erfc(): RealNumber {
    const x = this._value;
    if (x === 0) {
      return new RealNumber(this._parent, 1);
    }
    return new RealNumber(this._parent, 1 - this.erf().toNumber());
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

    // Miller's algorithm using downward recurrence
    const ax = Math.abs(x);
    if (n > Math.floor(ax)) {
      const ACC = 40;
      const BIGNO = 1e10;
      const BIGNI = 1e-10;
      const tox = 2 / ax;
      let bjp: number;
      let bj = 1;
      let bjm: number;
      let sum = 0;
      let ans = 0;
      const m = 2 * Math.floor((n + Math.floor(Math.sqrt(ACC * n))) / 2);

      for (let j = m; j > 0; j--) {
        bjm = j * tox * bj - bjp!;
        bjp = bj;
        bj = bjm;
        if (Math.abs(bj) > BIGNO) {
          bj *= BIGNI;
          bjp! *= BIGNI;
          ans *= BIGNI;
          sum *= BIGNI;
        }
        if (j % 2 !== 0) sum += bj;
        if (j === n) ans = bjp;
      }
      sum = 2 * sum - bj;
      ans /= sum;
      return new RealNumber(this._parent, x < 0 && n % 2 !== 0 ? -sign * ans : sign * ans);
    } else {
      // Upward recurrence
      const tox = 2 / ax;
      let bjm = this.j0().toNumber();
      let bj = this.j1().toNumber();
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

    if (s === 1) {
      return new RealNumber(this._parent, Number.POSITIVE_INFINITY);
    }

    // For negative even integers, zeta is 0
    if (s < 0 && Number.isInteger(s) && s % 2 === 0) {
      return new RealNumber(this._parent, 0);
    }

    // For s > 1, use direct summation with Euler-Maclaurin correction
    if (s > 1) {
      // Use many more terms for better accuracy
      const n = 10000;
      let sum = 0;
      for (let k = 1; k <= n; k++) {
        sum += k ** -s;
      }
      // Euler-Maclaurin tail correction with more terms
      // zeta(s) = sum_{k=1}^n k^-s + n^(1-s)/(s-1) + n^-s/2 + sum_k B_{2k}/(2k)! * s(s+1)...(s+2k-2) * n^-(s+2k-1)
      const nms = n ** -s;
      sum += n ** (1 - s) / (s - 1);
      sum += nms / 2;
      // Bernoulli correction terms
      sum += ((s / 12) * nms) / n;
      sum -= (((s * (s + 1) * (s + 2)) / 720) * nms) / n ** 3;
      sum += (((s * (s + 1) * (s + 2) * (s + 3) * (s + 4)) / 30240) * nms) / n ** 5;
      return new RealNumber(this._parent, sum);
    }

    // Use functional equation for s < 1: zeta(s) = 2^s * pi^(s-1) * sin(pi*s/2) * gamma(1-s) * zeta(1-s)
    const oneMinusS = 1 - s;
    const zetaOneMinusS = new RealNumber(this._parent, oneMinusS).zeta().toNumber();
    const gammaOneMinusS = new RealNumber(this._parent, oneMinusS).gamma().toNumber();
    const result = 2 ** s * PI ** (s - 1) * Math.sin((PI * s) / 2) * gammaOneMinusS * zetaOneMinusS;
    return new RealNumber(this._parent, result);
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
   * Return the simplest rational within the interval defined by this float's uncertainty.
   * @see Reference: sage/rings/real_mpfr.pyx:simplest_rational
   */
  simplest_rational(): [bigint, bigint] {
    // For standard float precision, find a simple rational approximation
    // Use continued fraction to find the simplest rational
    return this.nearby_rational(undefined, 10000000n);
  }

  /**
   * Return the nearby rational with bounded denominator.
   * Uses continued fraction expansion to find best rational approximation.
   * @see Reference: sage/rings/real_mpfr.pyx:nearby_rational
   */
  nearby_rational(max_error?: number, max_denominator?: bigint): [bigint, bigint] {
    const x = this._value;

    if (!Number.isFinite(x)) {
      throw new Error('Cannot convert NaN or infinity to rational');
    }

    const maxDenom = max_denominator ?? 1000000n;
    const maxErr = max_error ?? Number.EPSILON;

    // Simple continued fraction algorithm
    let a = Math.floor(x);
    let h1 = BigInt(a);
    let h2 = 1n;
    let k1 = 1n;
    let k2 = 0n;

    let val = x - a;
    while (k1 <= maxDenom && val !== 0) {
      val = 1 / val;
      a = Math.floor(val);
      const h0 = h1;
      const k0 = k1;
      h1 = BigInt(a) * h1 + h2;
      k1 = BigInt(a) * k1 + k2;
      h2 = h0;
      k2 = k0;

      if (k1 > maxDenom) {
        h1 = h2;
        k1 = k2;
        break;
      }

      val = val - a;
      if (Math.abs(Number(h1) / Number(k1) - x) < maxErr) {
        break;
      }
    }

    return [h1, k1];
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
   * In the real field, only 1 and -1 have finite multiplicative order.
   * @see Reference: sage/rings/real_mpfr.pyx:multiplicative_order
   */
  multiplicative_order(): number {
    if (this._value === 1) {
      return 1;
    }
    if (this._value === -1) {
      return 2;
    }
    throw new Error('Element does not have finite multiplicative order');
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

    if (!Number.isFinite(x)) {
      return new RealNumber(this._parent, Number.NaN);
    }

    if (x === 0) {
      return new RealNumber(this._parent, Number.MIN_VALUE);
    }

    const next = this.nextabove().toNumber();
    const prev = this.nextbelow().toNumber();

    // ulp is the minimum of the two distances
    const ulpValue = Math.min(Math.abs(next - x), Math.abs(x - prev));
    return new RealNumber(this._parent, ulpValue);
  }

  /**
   * Return the machine epsilon for this precision.
   * This is 2^(1-p) where p is the precision in bits.
   * @see Reference: sage/rings/real_mpfr.pyx:epsilon
   */
  epsilon(): RealNumber {
    // For standard 53-bit precision (IEEE 754 double)
    // epsilon = 2^-52 ≈ 2.220446049250313e-16
    const prec = this._parent.precision();
    const eps = 2 ** (1 - prec);
    return new RealNumber(this._parent, eps);
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
   * Return the floating point rank.
   * This is an integer that uniquely identifies the float's position in the
   * ordered set of all IEEE 754 doubles.
   * @see Reference: sage/rings/real_mpfr.pyx:fp_rank
   */
  fp_rank(): bigint {
    const x = this._value;

    if (Number.isNaN(x)) {
      throw new Error('NaN does not have a floating point rank');
    }

    // Get the bit representation
    const buffer = new ArrayBuffer(8);
    const f64 = new Float64Array(buffer);
    const i64 = new BigInt64Array(buffer);

    f64[0] = x;
    const bits = i64[0]!;

    // For positive numbers and positive zero, the rank is the bit pattern
    // For negative numbers, we need to flip the pattern
    if (x >= 0) {
      return bits;
    } else {
      // Negative numbers: flip sign bit and negate
      return -(bits & 0x7fffffffffffffffn);
    }
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
    if (n < 1) {
      throw new Error('degree must be at least 1');
    }

    const z = this._value;

    // Handle special cases
    if (!Number.isFinite(z)) {
      throw new Error('cannot compute algebraic dependency for non-finite number');
    }

    // For integers, return x - value
    if (Number.isInteger(z)) {
      return [-BigInt(Math.round(z)), 1n];
    }

    // Use LLL-based algorithm (matching SageMath's implementation)
    // Reference: sage/arith/misc.py:algebraic_dependency
    //
    // Build matrix M where:
    // - M is (n+1) x (n+2)
    // - First n+1 columns are identity
    // - Last column contains scaled powers: [r, r*z, r*z^2, ..., r*z^n]
    //   where r = 2^prec
    const prec = this._parent.precision() - 6;
    const dim = n + 1;

    // Build the matrix
    const M: number[][] = [];
    let r = 2 ** prec;

    for (let k = 0; k <= n; k++) {
      const row: number[] = [];
      for (let j = 0; j < dim; j++) {
        row.push(k === j ? 1 : 0);
      }
      // Last column: scaled power r * z^k
      row.push(Math.round(r));
      M.push(row);
      r *= z; // r *= z for next iteration
    }

    // Perform LLL reduction
    const lllReduced = this._lllReduce(M, 0.75);

    // The first row of the LLL-reduced basis gives the coefficients
    let coeffs = lllReduced[0]!.slice(0, dim);

    // If all coefficients except the first are zero, try the second row
    let allButFirstZero = true;
    for (let i = 1; i < coeffs.length; i++) {
      if (Math.abs(Math.round(coeffs[i]!)) !== 0) {
        allButFirstZero = false;
        break;
      }
    }
    if (allButFirstZero && lllReduced.length > 1) {
      coeffs = lllReduced[1]!.slice(0, dim);
    }

    // Convert to bigint
    const result: bigint[] = coeffs.map((c) => BigInt(Math.round(c)));

    // Make sure leading coefficient is positive (standard form)
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
   * Simple LLL reduction for the algebraic_dependency method.
   * Uses Lenstra-Lenstra-Lovasz algorithm with Gram-Schmidt orthogonalization.
   */
  private _lllReduce(M: number[][], delta: number): number[][] {
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
