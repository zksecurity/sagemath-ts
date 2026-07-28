/**
 * @module sage/rings/complex_mpfr
 * @description Arbitrary precision floating point complex numbers using GNU MPFR
 *
 * Port of: sage/rings/complex_mpfr.pyx
 * Reference: reference/sage/src/sage/rings/complex_mpfr.pyx
 */

import { NotImplementedError, ValueError } from '../errors.js';
import { IntegerMatrix, LLL } from '../matrix/matrix_integer.js';
import { Integer } from './integer_ring.js';
import { PolynomialRing } from './polynomial/polynomial_ring.js';
import { RealField, RealNumber } from './real_mpfr.js';

/**
 * A minimal `ZZ` coefficient ring for the `ZZ[x]` used to factor the output of
 * {@link ComplexNumber.algebraic_dependency}.
 */
const ZZ_FOR_ALGDEP = {
  zero: () => new Integer(0n),
  one: () => new Integer(1n),
  __call__: (x: unknown): Integer => {
    if (x instanceof Integer) return x;
    if (typeof x === 'bigint') return new Integer(x);
    if (typeof x === 'number') return new Integer(BigInt(x));
    if (typeof x === 'string') return new Integer(BigInt(x));
    throw new ValueError(`cannot coerce ${String(x)} to Integer`);
  },
  is_field: () => false,
  is_integral_domain: () => true,
  characteristic: () => 0n,
  toString: () => 'Integer Ring',
};

/**
 * `B_{2n} / (2n+1)!` for `n = 1, 2, ...`, the coefficients of the Bernoulli
 * series for the dilogarithm
 * `Li_2(z) = u - u^2/4 + sum_{n>=1} B_{2n} u^{2n+1}/(2n+1)!`, `u = -log(1-z)`.
 */
const DILOG_BERNOULLI = [
  0.027777777777777776, -0.0002777777777777778, 4.72411186696901e-6, -9.185773074661964e-8,
  1.8978869988971e-9, -4.0647616451442256e-11, 8.921691020456452e-13, -1.9939295860721074e-14,
  4.518980029619918e-16, -1.0356517612181247e-17, 2.395218621026187e-19, -5.581785874325009e-21,
  1.3091507554183213e-22, -3.0874198024267403e-24, 7.315975652702203e-26, -1.740845657234001e-27,
  4.1576356446139e-29, -9.962148488284622e-31, 2.3940344248961652e-32, -5.76834735536739e-34,
  1.393179479647008e-35, -3.3721219654850894e-37,
];

/**
 * An approximation to the field of complex numbers using floating
 * point numbers with any specified precision.
 * @see Reference: sage/rings/complex_mpfr.pyx:ComplexField_class
 */
export class ComplexField {
  private readonly _prec: number;
  private _realField?: RealField;

  constructor(prec: number = 53) {
    this._prec = prec;
  }

  /**
   * Return the precision of this field.
   * @see Reference: sage/rings/complex_mpfr.pyx:prec
   */
  prec(): number {
    return this._prec;
  }

  /**
   * Alias for prec().
   * @see Reference: sage/rings/complex_mpfr.pyx:precision
   */
  precision(): number {
    return this._prec;
  }

  /**
   * Return whether this field is exact (always false).
   * @see Reference: sage/rings/complex_mpfr.pyx:is_exact
   */
  is_exact(): boolean {
    return false;
  }

  /**
   * Return the characteristic (which is 0).
   * @see Reference: sage/rings/complex_mpfr.pyx:characteristic
   */
  characteristic(): bigint {
    return 0n;
  }

  /**
   * Return the underlying real field.
   * @see Reference: sage/rings/complex_mpfr.pyx:_real_field
   */
  real_field(): RealField {
    if (!this._realField) {
      this._realField = new RealField(this._prec);
    }
    return this._realField;
  }

  /**
   * Return the generator I.
   * @see Reference: sage/rings/complex_mpfr.pyx:gen
   */
  gen(): ComplexNumber {
    return new ComplexNumber(this, 0, 1);
  }

  /**
   * Return the number of generators (which is 1).
   * @see Reference: sage/rings/complex_mpfr.pyx:ngens
   */
  ngens(): number {
    return 1;
  }

  /**
   * Return complex field with specified precision.
   * @see Reference: sage/rings/complex_mpfr.pyx:to_prec
   */
  to_prec(prec: number): ComplexField {
    return new ComplexField(prec);
  }

  /**
   * Return pi as a complex number.
   * @see Reference: sage/rings/complex_mpfr.pyx:pi
   */
  pi(): ComplexNumber {
    return new ComplexNumber(this, Math.PI, 0);
  }

  /**
   * Return a primitive n-th root of unity.
   * @see Reference: sage/rings/complex_mpfr.pyx:zeta
   */
  zeta(n: number = 2): ComplexNumber {
    if (n === 1) {
      return new ComplexNumber(this, 1, 0);
    } else if (n === 2) {
      return new ComplexNumber(this, -1, 0);
    } else {
      // e^(2*pi*i/n) = cos(2*pi/n) + i*sin(2*pi/n)
      const theta = (2 * Math.PI) / n;
      return new ComplexNumber(this, Math.cos(theta), Math.sin(theta));
    }
  }

  /**
   * Return a random element.
   * @see Reference: sage/rings/complex_mpfr.pyx:random_element
   */
  random_element(component_max: number = 1): ComplexNumber {
    const real = (Math.random() * 2 - 1) * component_max;
    const imag = (Math.random() * 2 - 1) * component_max;
    return new ComplexNumber(this, real, imag);
  }

  /**
   * Return the algebraic closure.
   * @see Reference: sage/rings/complex_mpfr.pyx:algebraic_closure
   */
  algebraic_closure(): ComplexField {
    return this;
  }

  /**
   * Create a complex number.
   * @see Reference: sage/rings/complex_mpfr.pyx:__call__
   */
  __call__(real: number | bigint | string = 0, imag?: number | bigint | string): ComplexNumber {
    const r =
      typeof real === 'bigint'
        ? Number(real)
        : typeof real === 'string'
          ? Number.parseFloat(real)
          : real;
    const i =
      imag === undefined
        ? 0
        : typeof imag === 'bigint'
          ? Number(imag)
          : typeof imag === 'string'
            ? Number.parseFloat(imag)
            : imag;
    return new ComplexNumber(this, r, i);
  }

  toString(): string {
    return `Complex Field with ${this._prec} bits of precision`;
  }
}

/**
 * A complex number with arbitrary precision real and imaginary parts.
 * @see Reference: sage/rings/complex_mpfr.pyx:ComplexNumber
 */
export class ComplexNumber {
  private readonly _parent: ComplexField;
  private readonly _real: number;
  private readonly _imag: number;

  constructor(parent: ComplexField, real: number = 0, imag: number = 0) {
    this._parent = parent;
    this._real = real;
    this._imag = imag;
  }

  /**
   * Return the precision.
   * @see Reference: sage/rings/complex_mpfr.pyx:prec
   */
  prec(): number {
    return this._parent.prec();
  }

  /**
   * Return the real part.
   * @see Reference: sage/rings/complex_mpfr.pyx:real
   */
  real(): number {
    return this._real;
  }

  /**
   * Return the imaginary part.
   * @see Reference: sage/rings/complex_mpfr.pyx:imag
   */
  imag(): number {
    return this._imag;
  }

  /**
   * Return the norm (square of absolute value).
   * @see Reference: sage/rings/complex_mpfr.pyx:norm
   */
  norm(): number {
    return this._real * this._real + this._imag * this._imag;
  }

  /**
   * Return the absolute value.
   * @see Reference: sage/rings/complex_mpfr.pyx:__abs__
   */
  abs(): number {
    return Math.sqrt(this.norm());
  }

  /**
   * Return the argument (phase angle).
   * @see Reference: sage/rings/complex_mpfr.pyx:argument
   */
  argument(): number {
    return Math.atan2(this._imag, this._real);
  }

  /**
   * Alias for argument().
   * @see Reference: sage/rings/complex_mpfr.pyx:arg
   */
  arg(): number {
    return this.argument();
  }

  /**
   * Return the complex conjugate.
   * @see Reference: sage/rings/complex_mpfr.pyx:conjugate
   */
  conjugate(): ComplexNumber {
    return new ComplexNumber(this._parent, this._real, -this._imag);
  }

  /**
   * Return the negation.
   * @see Reference: sage/rings/complex_mpfr.pyx:__neg__
   */
  neg(): ComplexNumber {
    return new ComplexNumber(this._parent, -this._real, -this._imag);
  }

  /**
   * Return the multiplicative inverse.
   * @see Reference: sage/rings/complex_mpfr.pyx:__invert__
   */
  inv(): ComplexNumber {
    const n = this.norm();
    return new ComplexNumber(this._parent, this._real / n, -this._imag / n);
  }

  /**
   * Add two complex numbers.
   */
  add(other: ComplexNumber): ComplexNumber {
    return new ComplexNumber(this._parent, this._real + other._real, this._imag + other._imag);
  }

  /**
   * Subtract two complex numbers.
   */
  sub(other: ComplexNumber): ComplexNumber {
    return new ComplexNumber(this._parent, this._real - other._real, this._imag - other._imag);
  }

  /**
   * Multiply two complex numbers.
   */
  mul(other: ComplexNumber): ComplexNumber {
    const r = this._real * other._real - this._imag * other._imag;
    const i = this._real * other._imag + this._imag * other._real;
    return new ComplexNumber(this._parent, r, i);
  }

  /**
   * Divide two complex numbers.
   */
  div(other: ComplexNumber): ComplexNumber {
    const n = other.norm();
    const r = (this._real * other._real + this._imag * other._imag) / n;
    const i = (this._imag * other._real - this._real * other._imag) / n;
    return new ComplexNumber(this._parent, r, i);
  }

  /**
   * Return the exponential.
   * exp(a + bi) = exp(a) * (cos(b) + i*sin(b))
   * @see Reference: sage/rings/complex_mpfr.pyx:exp
   */
  exp(): ComplexNumber {
    const r = Math.exp(this._real);
    return new ComplexNumber(this._parent, r * Math.cos(this._imag), r * Math.sin(this._imag));
  }

  /**
   * Return the logarithm.
   * log(z) = log(|z|) + i*arg(z)
   * @see Reference: sage/rings/complex_mpfr.pyx:log
   */
  log(base?: number): ComplexNumber {
    if (Number.isNaN(this._real) || Number.isNaN(this._imag)) {
      return new ComplexNumber(this._parent, Number.NaN, Number.NaN);
    }
    const theta = this.argument();
    const rho = this.abs();
    const logRho = Math.log(rho);
    if (base === undefined) {
      return new ComplexNumber(this._parent, logRho, theta);
    } else {
      const logBase = Math.log(base);
      return new ComplexNumber(this._parent, logRho / logBase, theta / logBase);
    }
  }

  /**
   * Return the square root.
   * sqrt(z) = sqrt(|z|) * e^(i*arg(z)/2)
   * @see Reference: sage/rings/complex_mpfr.pyx:sqrt
   */
  sqrt(all: boolean = false): ComplexNumber | ComplexNumber[] {
    // Handle purely real case
    if (this._imag === 0) {
      if (this._real >= 0) {
        const z = new ComplexNumber(this._parent, Math.sqrt(this._real), 0);
        if (all) {
          return this._real === 0 ? [z] : [z, z.neg()];
        }
        return z;
      } else {
        // sqrt of negative real
        const z = new ComplexNumber(this._parent, 0, Math.sqrt(-this._real));
        if (all) {
          return [z, z.neg()];
        }
        return z;
      }
    }

    // General complex case using formula:
    // a^2 = (x + sqrt(x^2+y^2))/2
    // b = y/(2a)
    const x = this._real;
    const y = this._imag;
    const r = Math.sqrt(x * x + y * y);

    // Use stable computation avoiding cancellation near negative real axis
    const avoidBranch = x < 0 && Math.abs(y) < Math.abs(x);

    let a: number;
    let b: number;
    if (avoidBranch) {
      // x + sqrt(x^2+y^2) is numerically unstable for x near the negative real
      // axis, so we compute sqrt(-z) and shift by i at the end.
      const a2 = (r - x) / 2;
      a = Math.sqrt(a2);
      b = y / (2 * a);
      // mpfr_swap(re, im): note that y (hence b) was never negated, so we have
      // a + b*i = i*sqrt(self); swapping the parts (WITHOUT negating either)
      // divides by i.
      const tempA = a;
      a = b;
      b = tempA;
      // If we were below the branch cut, we want the other branch.
      if (y < 0) {
        a = -a;
        b = -b;
      }
    } else {
      const a2 = (r + x) / 2;
      a = Math.sqrt(a2);
      b = y / (2 * a);
    }

    const z = new ComplexNumber(this._parent, a, b);
    if (all) {
      return [z, z.neg()];
    }
    return z;
  }

  /**
   * Return the n-th root.
   * nth_root(z) = |z|^(1/n) * e^(i*arg(z)/n)
   * @see Reference: sage/rings/complex_mpfr.pyx:nth_root
   */
  nth_root(n: number, all: boolean = false): ComplexNumber | ComplexNumber[] {
    if (this._real === 0 && this._imag === 0) {
      const z = new ComplexNumber(this._parent, 0, 0);
      return all ? [z] : z;
    }

    const rho = this.abs();
    const arg = this.argument() / n;
    const r = rho ** (1 / n);

    const z = new ComplexNumber(this._parent, r * Math.cos(arg), r * Math.sin(arg));

    if (!all) {
      return z;
    }

    // Return all n-th roots
    const roots: ComplexNumber[] = [z];
    const theta = (2 * Math.PI) / n;
    let currentArg = arg;
    for (let k = 1; k < n; k++) {
      currentArg += theta;
      roots.push(
        new ComplexNumber(this._parent, r * Math.cos(currentArg), r * Math.sin(currentArg))
      );
    }
    return roots;
  }

  /**
   * Return the cosine.
   * cos(a + bi) = cosh(b)*cos(a) - i*sinh(b)*sin(a)
   * @see Reference: sage/rings/complex_mpfr.pyx:cos
   */
  cos(): ComplexNumber {
    const sh = Math.sinh(this._imag);
    const ch = Math.cosh(this._imag);
    const sa = Math.sin(this._real);
    const ca = Math.cos(this._real);
    return new ComplexNumber(this._parent, ch * ca, -sh * sa);
  }

  /**
   * Return the sine.
   * sin(a + bi) = cosh(b)*sin(a) + i*sinh(b)*cos(a)
   * @see Reference: sage/rings/complex_mpfr.pyx:sin
   */
  sin(): ComplexNumber {
    const sh = Math.sinh(this._imag);
    const ch = Math.cosh(this._imag);
    const sa = Math.sin(this._real);
    const ca = Math.cos(this._real);
    return new ComplexNumber(this._parent, ch * sa, sh * ca);
  }

  /**
   * Return the tangent.
   * tan(a + bi) = [cos(a)*sin(a) + i*cosh(b)*sinh(b)] / [sinh^2(b) + cos^2(a)]
   * @see Reference: sage/rings/complex_mpfr.pyx:tan
   */
  tan(): ComplexNumber {
    const sh = Math.sinh(this._imag);
    const ch = Math.cosh(this._imag);
    const sa = Math.sin(this._real);
    const ca = Math.cos(this._real);
    const denom = sh * sh + ca * ca;
    return new ComplexNumber(this._parent, (ca * sa) / denom, (ch * sh) / denom);
  }

  /**
   * Return the hyperbolic cosine.
   * cosh(a + bi) = cosh(a)*cos(b) + i*sinh(a)*sin(b)
   * @see Reference: sage/rings/complex_mpfr.pyx:cosh
   */
  cosh(): ComplexNumber {
    const sha = Math.sinh(this._real);
    const cha = Math.cosh(this._real);
    const sb = Math.sin(this._imag);
    const cb = Math.cos(this._imag);
    return new ComplexNumber(this._parent, cha * cb, sha * sb);
  }

  /**
   * Return the hyperbolic sine.
   * sinh(a + bi) = sinh(a)*cos(b) + i*cosh(a)*sin(b)
   * @see Reference: sage/rings/complex_mpfr.pyx:sinh
   */
  sinh(): ComplexNumber {
    const sha = Math.sinh(this._real);
    const cha = Math.cosh(this._real);
    const sb = Math.sin(this._imag);
    const cb = Math.cos(this._imag);
    return new ComplexNumber(this._parent, sha * cb, cha * sb);
  }

  /**
   * Return the hyperbolic tangent.
   * tanh(a + bi) = [cosh(a)*sinh(a) + i*cos(b)*sin(b)] / [sinh^2(a) + cos^2(b)]
   * @see Reference: sage/rings/complex_mpfr.pyx:tanh
   */
  tanh(): ComplexNumber {
    const sha = Math.sinh(this._real);
    const cha = Math.cosh(this._real);
    const sb = Math.sin(this._imag);
    const cb = Math.cos(this._imag);
    const denom = sha * sha + cb * cb;
    return new ComplexNumber(this._parent, (cha * sha) / denom, (cb * sb) / denom);
  }

  /**
   * Return the arccosine.
   * arccos(z) = -i * log(z + i*sqrt(1 - z^2))
   * @see Reference: sage/rings/complex_mpfr.pyx:arccos
   */
  arccos(): ComplexNumber {
    // arccos(z) = -i * log(z + i*sqrt(1 - z^2))
    const one = new ComplexNumber(this._parent, 1, 0);
    const i = new ComplexNumber(this._parent, 0, 1);
    const z2 = this.mul(this);
    const oneMinusZ2 = one.sub(z2);
    const sqrtPart = oneMinusZ2.sqrt() as ComplexNumber;
    const inside = this.add(i.mul(sqrtPart));
    const logPart = inside.log();
    // -i * log = multiply by -i: (a + bi) * (-i) = b - ai
    return new ComplexNumber(this._parent, logPart.imag(), -logPart.real());
  }

  /**
   * Return the arcsine.
   * arcsin(z) = -i * log(i*z + sqrt(1 - z^2))
   * @see Reference: sage/rings/complex_mpfr.pyx:arcsin
   */
  arcsin(): ComplexNumber {
    // arcsin(z) = -i * log(i*z + sqrt(1 - z^2))
    const one = new ComplexNumber(this._parent, 1, 0);
    const i = new ComplexNumber(this._parent, 0, 1);
    const z2 = this.mul(this);
    const oneMinusZ2 = one.sub(z2);
    const sqrtPart = oneMinusZ2.sqrt() as ComplexNumber;
    const iz = i.mul(this);
    const inside = iz.add(sqrtPart);
    const logPart = inside.log();
    // -i * log = (a + bi) * (-i) = b - ai
    return new ComplexNumber(this._parent, logPart.imag(), -logPart.real());
  }

  /**
   * Return the arctangent.
   * arctan(z) = (i/2) * log((1-iz)/(1+iz))
   * @see Reference: sage/rings/complex_mpfr.pyx:arctan
   */
  arctan(): ComplexNumber {
    // arctan(z) = (i/2) * log((1-iz)/(1+iz))
    const one = new ComplexNumber(this._parent, 1, 0);
    const i = new ComplexNumber(this._parent, 0, 1);
    const iz = i.mul(this);
    const num = one.sub(iz);
    const den = one.add(iz);
    const logPart = num.div(den).log();
    // (i/2) * log = (0 + 0.5i) * (a + bi) = -0.5*b + 0.5*a*i
    return new ComplexNumber(this._parent, -0.5 * logPart.imag(), 0.5 * logPart.real());
  }

  /**
   * Return the inverse hyperbolic cosine.
   * arccosh(z) = log(z + sqrt(z^2 - 1))
   * @see Reference: sage/rings/complex_mpfr.pyx:arccosh
   */
  arccosh(): ComplexNumber {
    // arccosh(z) = log(z + sqrt(z^2 - 1))
    const one = new ComplexNumber(this._parent, 1, 0);
    const z2 = this.mul(this);
    const z2minus1 = z2.sub(one);
    const sqrtPart = z2minus1.sqrt() as ComplexNumber;
    return this.add(sqrtPart).log();
  }

  /**
   * Return the inverse hyperbolic sine.
   * arcsinh(z) = log(z + sqrt(z^2 + 1))
   * @see Reference: sage/rings/complex_mpfr.pyx:arcsinh
   */
  arcsinh(): ComplexNumber {
    // arcsinh(z) = log(z + sqrt(z^2 + 1))
    const one = new ComplexNumber(this._parent, 1, 0);
    const z2 = this.mul(this);
    const z2plus1 = z2.add(one);
    const sqrtPart = z2plus1.sqrt() as ComplexNumber;
    return this.add(sqrtPart).log();
  }

  /**
   * Return the inverse hyperbolic tangent.
   * arctanh(z) = (1/2) * log((1+z)/(1-z))
   * @see Reference: sage/rings/complex_mpfr.pyx:arctanh
   */
  arctanh(): ComplexNumber {
    // arctanh(z) = (1/2) * log((1+z)/(1-z))
    const one = new ComplexNumber(this._parent, 1, 0);
    const num = one.add(this);
    const den = one.sub(this);
    const logPart = num.div(den).log();
    return new ComplexNumber(this._parent, 0.5 * logPart.real(), 0.5 * logPart.imag());
  }

  /**
   * Return the cotangent.
   * cot(z) = 1/tan(z)
   * @see Reference: sage/rings/complex_mpfr.pyx:cot
   */
  cot(): ComplexNumber {
    return this.tan().inv();
  }

  /**
   * Return the hyperbolic cotangent.
   * coth(z) = 1/tanh(z)
   * @see Reference: sage/rings/complex_mpfr.pyx:coth
   */
  coth(): ComplexNumber {
    return this.tanh().inv();
  }

  /**
   * Return the secant.
   * sec(z) = 1/cos(z)
   * @see Reference: sage/rings/complex_mpfr.pyx:sec
   */
  sec(): ComplexNumber {
    return this.cos().inv();
  }

  /**
   * Return the hyperbolic secant.
   * sech(z) = 1/cosh(z)
   * @see Reference: sage/rings/complex_mpfr.pyx:sech
   */
  sech(): ComplexNumber {
    return this.cosh().inv();
  }

  /**
   * Return the cosecant.
   * csc(z) = 1/sin(z)
   * @see Reference: sage/rings/complex_mpfr.pyx:csc
   */
  csc(): ComplexNumber {
    return this.sin().inv();
  }

  /**
   * Return the hyperbolic cosecant.
   * csch(z) = 1/sinh(z)
   * @see Reference: sage/rings/complex_mpfr.pyx:csch
   */
  csch(): ComplexNumber {
    return this.sinh().inv();
  }

  /**
   * Return the gamma function.
   * Uses Stirling's approximation for complex numbers.
   * @see Reference: sage/rings/complex_mpfr.pyx:gamma
   */
  gamma(): ComplexNumber {
    // Handle special cases
    if (this._imag === 0) {
      if (this._real > 0 && Number.isInteger(this._real)) {
        // Gamma(n) = (n-1)! for positive integers
        let result = 1;
        for (let i = 2; i < this._real; i++) {
          result *= i;
        }
        return new ComplexNumber(this._parent, result, 0);
      }
      if (this._real <= 0 && Number.isInteger(this._real)) {
        // Gamma has poles at non-positive integers
        return new ComplexNumber(this._parent, Number.POSITIVE_INFINITY, 0);
      }
    }

    // Use Lanczos approximation for complex gamma
    // This is a standard approximation for the gamma function
    const g = 7;
    const c = [
      0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313,
      -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6,
      1.5056327351493116e-7,
    ];

    let z = this;
    if (this._real < 0.5) {
      // Use reflection formula: Gamma(z)*Gamma(1-z) = pi/sin(pi*z)
      const piZ = new ComplexNumber(this._parent, Math.PI * z._real, Math.PI * z._imag);
      const sinPiZ = piZ.sin();
      const oneMinusZ = new ComplexNumber(this._parent, 1 - z._real, -z._imag);
      const gammaOneMinusZ = oneMinusZ.gamma();
      const pi = new ComplexNumber(this._parent, Math.PI, 0);
      return pi.div(sinPiZ.mul(gammaOneMinusZ));
    }

    z = new ComplexNumber(this._parent, z._real - 1, z._imag);
    let x = new ComplexNumber(this._parent, c[0], 0);
    for (let i = 1; i < g + 2; i++) {
      const denom = new ComplexNumber(this._parent, z._real + i, z._imag);
      x = x.add(new ComplexNumber(this._parent, c[i], 0).div(denom));
    }

    const t = new ComplexNumber(this._parent, z._real + g + 0.5, z._imag);
    const sqrtTwoPi = new ComplexNumber(this._parent, Math.sqrt(2 * Math.PI), 0);
    const tPowZPlusHalf = t
      .log()
      .mul(new ComplexNumber(this._parent, z._real + 0.5, z._imag))
      .exp();
    const expNegT = t.neg().exp();

    return sqrtTwoPi.mul(tPowZPlusHalf).mul(expNegT).mul(x);
  }

  /**
   * Return the incomplete gamma function.
   * Uses series expansion for small t, continued fraction for large t.
   * @see Reference: sage/rings/complex_mpfr.pyx:gamma_inc
   */
  gamma_inc(t: ComplexNumber): ComplexNumber {
    const x = t;

    if (x._real === 0 && x._imag === 0) {
      return this.gamma();
    }

    // x^a * e^(-x), the common prefactor of both expansions.
    const prefactor = x.log().mul(this).exp().mul(x.neg().exp());

    // For |x| > |a| + 1 the ascending series for the *lower* incomplete gamma
    // cancels catastrophically against Gamma(a) (Gamma(2,50) came out as
    // 7.6e-11 instead of 9.8e-21).  Use Legendre's continued fraction for the
    // upper incomplete gamma directly:
    //   Gamma(a,x) = x^a e^{-x} / (x+1-a - 1*(1-a)/(x+3-a - 2*(2-a)/(...)))
    // evaluated by the modified Lentz algorithm.
    if (x.abs() > this.abs() + 1) {
      const TINY = 1e-300;
      const tiny = new ComplexNumber(this._parent, TINY, 0);
      const two = new ComplexNumber(this._parent, 2, 0);
      const one = new ComplexNumber(this._parent, 1, 0);

      let b = x.add(one).sub(this);
      let c = new ComplexNumber(this._parent, 1 / TINY, 0);
      let d = b.inv();
      let h = d;

      for (let i = 1; i < 1000; i++) {
        const iC = new ComplexNumber(this._parent, i, 0);
        const an = iC.mul(iC.sub(this)).neg();
        b = b.add(two);
        d = an.mul(d).add(b);
        if (d.abs() < TINY) d = tiny;
        c = b.add(an.div(c));
        if (c.abs() < TINY) c = tiny;
        d = d.inv();
        const del = d.mul(c);
        h = h.mul(del);
        if (del.sub(one).abs() < 1e-17) break;
      }

      return prefactor.mul(h);
    }

    // Ascending series for the lower incomplete gamma:
    //   gamma(a,x) = x^a e^{-x} / a * sum_{n>=0} x^n / ((a+1)...(a+n))
    // then Gamma(a,x) = Gamma(a) - gamma(a,x).
    const maxIter = 1000;
    const eps = 1e-17;

    let term = new ComplexNumber(this._parent, 1, 0);
    let sum = new ComplexNumber(this._parent, 1, 0);
    let aPlusN = new ComplexNumber(this._parent, this._real, this._imag);

    for (let n = 1; n < maxIter; n++) {
      aPlusN = new ComplexNumber(this._parent, aPlusN._real + 1, aPlusN._imag);
      term = term.mul(x).div(aPlusN);
      sum = sum.add(term);
      if (term.abs() < eps * sum.abs()) break;
    }

    const gammaLower = prefactor.mul(sum).div(this);
    return this.gamma().sub(gammaLower);
  }

  /**
   * Return the Riemann zeta function.
   * Uses the Borwein algorithm for fast convergence.
   * @see Reference: sage/rings/complex_mpfr.pyx:zeta
   */
  zeta(): ComplexNumber {
    // Special case: zeta(1) is infinity
    if (this._imag === 0 && this._real === 1) {
      return new ComplexNumber(this._parent, Number.POSITIVE_INFINITY, 0);
    }

    // Use reflection formula for Re(s) < 0:
    // zeta(s) = 2^s * pi^(s-1) * sin(pi*s/2) * Gamma(1-s) * zeta(1-s)
    if (this._real < 0) {
      const one = new ComplexNumber(this._parent, 1, 0);
      const oneMinusS = one.sub(this);
      const zetaOneMinusS = oneMinusS.zeta();
      const two = new ComplexNumber(this._parent, 2, 0);
      const pi = new ComplexNumber(this._parent, Math.PI, 0);
      const twoPowS = two.log().mul(this).exp();
      const piPowSMinus1 = pi.log().mul(one.sub(this).neg()).exp();
      const piSOver2 = pi.mul(this).mul(new ComplexNumber(this._parent, 0.5, 0));
      const sinPart = piSOver2.sin();
      const gammaOneMinusS = oneMinusS.gamma();
      return twoPowS.mul(piPowSMinus1).mul(sinPart).mul(gammaOneMinusS).mul(zetaOneMinusS);
    }

    // Use the Borwein algorithm with Chebyshev-like coefficients
    // This gives much faster convergence than the naive alternating series
    const n = 30; // Number of terms

    // Compute the d_k coefficients
    const d: number[] = new Array(n + 1);
    d[0] = 1;
    for (let k = 1; k <= n; k++) {
      d[k] = d[k - 1] + n ** k / this._factorial(k);
    }
    // Actually, use the simpler formula for eta acceleration
    // d_k = n * sum_{i=0}^{k} (n+i-1)! * 4^i / ((n-i)! * (2i)!)

    // Simpler approach: use Euler's transformation for the eta function
    // eta(s) = sum_{k=0}^inf (1/2^{k+1}) * sum_{j=0}^k C(k,j) * (-1)^j / (j+1)^s

    const maxK = 50;
    let eta = new ComplexNumber(this._parent, 0, 0);
    let halfPow = 0.5; // 1/2^{k+1}

    for (let k = 0; k < maxK; k++) {
      let innerSum = new ComplexNumber(this._parent, 0, 0);
      let binom = 1; // C(k, j)

      for (let j = 0; j <= k; j++) {
        const jPlus1 = new ComplexNumber(this._parent, j + 1, 0);
        const term = jPlus1.log().mul(this).neg().exp();
        const sign = j % 2 === 0 ? 1 : -1;
        innerSum = innerSum.add(term.mul(new ComplexNumber(this._parent, sign * binom, 0)));
        binom = (binom * (k - j)) / (j + 1);
      }

      eta = eta.add(innerSum.mul(new ComplexNumber(this._parent, halfPow, 0)));
      halfPow /= 2;

      if (halfPow < 1e-16) break;
    }

    // zeta(s) = eta(s) / (1 - 2^(1-s))
    const one = new ComplexNumber(this._parent, 1, 0);
    const two = new ComplexNumber(this._parent, 2, 0);
    const oneMinusS = one.sub(this);
    const twoPow1MinusS = two.log().mul(oneMinusS).exp();
    const denom = one.sub(twoPow1MinusS);

    return eta.div(denom);
  }

  /**
   * Helper function for factorial
   */
  private _factorial(n: number): number {
    let result = 1;
    for (let i = 2; i <= n; i++) {
      result *= i;
    }
    return result;
  }

  /**
   * Return the dilogarithm.
   * Li_2(z) = -integral from 0 to z of log(1-t)/t dt = sum_{k=1}^inf z^k/k^2
   * @see Reference: sage/rings/complex_mpfr.pyx:dilog
   */
  dilog(): ComplexNumber {
    // Special case: Li_2(1) = pi^2/6
    if (this._real === 1 && this._imag === 0) {
      return new ComplexNumber(this._parent, (Math.PI * Math.PI) / 6, 0);
    }

    // Special case: Li_2(0) = 0
    if (this._real === 0 && this._imag === 0) {
      return new ComplexNumber(this._parent, 0, 0);
    }

    const pi2Over6 = new ComplexNumber(this._parent, (Math.PI * Math.PI) / 6, 0);
    const half = new ComplexNumber(this._parent, 0.5, 0);
    const one = new ComplexNumber(this._parent, 1, 0);

    // Inversion, for EVERY |z| > 1 (not only |z| > 2; the series diverges as
    // soon as |z| > 1):
    //   Li_2(z) = -Li_2(1/z) - pi^2/6 - (1/2) log(-z)^2
    // The principal branch of log(-z) puts the cut of Li_2 on [1, +oo), which
    // is PARI's (and hence SageMath's) convention.
    if (this.abs() > 1) {
      // Build -z with a *positive* zero imaginary part when Im(z) == 0, so that
      // arg(-z) = +pi for real z > 1.  That selects the value of Li_2 on the
      // cut which PARI (and hence SageMath) reports:
      // CC(2).dilog() = 2.46740110027234 - 2.17758609030360*I.
      const negZ = new ComplexNumber(this._parent, -this._real, this._imag === 0 ? 0 : -this._imag);
      const logNegZ = negZ.log();
      return this.inv()
        .dilog()
        .neg()
        .sub(pi2Over6)
        .sub(half.mul(logNegZ.mul(logNegZ)));
    }

    // Reflection, to push Re(z) down to at most 1/2:
    //   Li_2(z) = pi^2/6 - log(z) log(1-z) - Li_2(1-z)
    // If |z| <= 1 and Re(z) > 1/2 then |1-z|^2 = 1 - 2 Re(z) + |z|^2 <= 1,
    // so the recursion lands in the series branch below.
    if (this._real > 0.5) {
      const oneMinusZ = one.sub(this);
      return pi2Over6.sub(this.log().mul(oneMinusZ.log())).sub(oneMinusZ.dilog());
    }

    // Here |z| <= 1 and Re(z) <= 1/2, so u = -log(1-z) satisfies |u| < 1.8 and
    // the Bernoulli series converges geometrically:
    //   Li_2(z) = sum_{k>=0} B_k u^{k+1}/(k+1)!  =  u - u^2/4
    //             + sum_{n>=1} B_{2n} u^{2n+1}/(2n+1)!
    return this._dilogBernoulli();
  }

  /**
   * Bernoulli-series evaluation of `Li_2`, valid for `|log(1-z)| < 2*pi`
   * (in practice: `|z| <= 1` and `Re(z) <= 1/2`).
   */
  private _dilogBernoulli(): ComplexNumber {
    const one = new ComplexNumber(this._parent, 1, 0);
    // u = -log(1 - z)
    const u = one.sub(this).log().neg();
    const u2 = u.mul(u);

    // u - u^2/4
    let sum = u.sub(u2.mul(new ComplexNumber(this._parent, 0.25, 0)));

    // sum_{n>=1} B_{2n}/(2n+1)! * u^{2n+1}
    let uPow = u.mul(u2); // u^3
    for (let n = 0; n < DILOG_BERNOULLI.length; n++) {
      const c = DILOG_BERNOULLI[n]!;
      const term = uPow.mul(new ComplexNumber(this._parent, c, 0));
      sum = sum.add(term);
      if (term.abs() < 1e-18 * sum.abs()) break;
      uPow = uPow.mul(u2);
    }
    return sum;
  }

  /**
   * Return the Dedekind eta function.
   * eta(z) = e^(pi*i*z/12) * prod_{n=1}^inf (1 - e^(2*pi*i*n*z))
   * @see Reference: sage/rings/complex_mpfr.pyx:eta
   */
  eta(omit_frac: boolean = false): ComplexNumber {
    // eta is only defined for z in upper half plane
    if (this._imag <= 0) {
      throw new Error('value must be in the upper half plane');
    }
    const i = new ComplexNumber(this._parent, 0, 1);
    const pi = Math.PI;

    // Compute the product: prod_{n=1}^inf (1 - e^(2*pi*i*n*z))
    const maxIter = 100;
    const eps = 1e-15;
    let prod = new ComplexNumber(this._parent, 1, 0);

    for (let n = 1; n <= maxIter; n++) {
      // q = e^(2*pi*i*z)
      // term = 1 - q^n = 1 - e^(2*pi*i*n*z)
      const exp_arg = new ComplexNumber(
        this._parent,
        -2 * pi * n * this._imag,
        2 * pi * n * this._real
      );
      const qPowN = exp_arg.exp();
      const one = new ComplexNumber(this._parent, 1, 0);
      const term = one.sub(qPowN);
      prod = prod.mul(term);

      if (Math.abs(qPowN.abs() - 0) < eps) break;
    }

    if (omit_frac) {
      return prod;
    }

    // Multiply by e^(pi*i*z/12)
    const fracExp = new ComplexNumber(
      this._parent,
      (-pi * this._imag) / 12,
      (pi * this._real) / 12
    );
    const frac = fracExp.exp();
    return frac.mul(prod);
  }

  /**
   * Return the arithmetic-geometric mean.
   * @see Reference: sage/rings/complex_mpfr.pyx:agm
   */
  agm(right: ComplexNumber, algorithm: string = 'optimal'): ComplexNumber {
    // AGM(a, 0) = AGM(0, a) = AGM(a, -a) = 0
    if (right._real === 0 && right._imag === 0) {
      return new ComplexNumber(this._parent, 0, 0);
    }
    if (this._real === 0 && this._imag === 0) {
      return new ComplexNumber(this._parent, 0, 0);
    }
    if (this._real === -right._real && this._imag === -right._imag) {
      return new ComplexNumber(this._parent, 0, 0);
    }

    let a = this;
    let b = right;
    const maxIter = 100;
    const eps = 1e-15;

    for (let i = 0; i < maxIter; i++) {
      const aNew = a.add(b).mul(new ComplexNumber(this._parent, 0.5, 0));
      const sqrtAB = a.mul(b).sqrt() as ComplexNumber;

      // Choose the sign of sqrt based on algorithm
      let bNew: ComplexNumber;
      if (algorithm === 'optimal') {
        // Choose sign so that |a_new - b_new| <= |a_new + b_new|
        const bPos = sqrtAB;
        const bNeg = sqrtAB.neg();
        const diffPos = aNew.sub(bPos).abs();
        const diffNeg = aNew.sub(bNeg).abs();
        bNew = diffPos <= diffNeg ? bPos : bNeg;
      } else if (algorithm === 'principal') {
        // Principal branch: Re(sqrt) >= 0
        bNew = sqrtAB._real >= 0 ? sqrtAB : sqrtAB.neg();
      } else {
        // Default to principal
        bNew = sqrtAB._real >= 0 ? sqrtAB : sqrtAB.neg();
      }

      // Check convergence
      const diff = aNew.sub(bNew).abs();
      if (diff < eps * aNew.abs()) {
        return aNew;
      }

      a = aNew;
      b = bNew;
    }

    return a;
  }

  /**
   * Check if this is a real number.
   * @see Reference: sage/rings/complex_mpfr.pyx:is_real
   */
  is_real(): boolean {
    return this._imag === 0;
  }

  /**
   * Check if this is purely imaginary.
   * @see Reference: sage/rings/complex_mpfr.pyx:is_imaginary
   */
  is_imaginary(): boolean {
    return this._real === 0 && this._imag !== 0;
  }

  /**
   * Check if this is an integer.
   * @see Reference: sage/rings/complex_mpfr.pyx:is_integer
   */
  is_integer(): boolean {
    return this._imag === 0 && Number.isInteger(this._real);
  }

  /**
   * Check if this is a square.
   * @see Reference: sage/rings/complex_mpfr.pyx:is_square
   */
  is_square(): boolean {
    return true; // All complex numbers are squares
  }

  /**
   * Check if this is positive infinity.
   * @see Reference: sage/rings/complex_mpfr.pyx:is_positive_infinity
   */
  is_positive_infinity(): boolean {
    return this._real === Number.POSITIVE_INFINITY && this._imag === 0;
  }

  /**
   * Check if this is negative infinity.
   * @see Reference: sage/rings/complex_mpfr.pyx:is_negative_infinity
   */
  is_negative_infinity(): boolean {
    return this._real === Number.NEGATIVE_INFINITY && this._imag === 0;
  }

  /**
   * Check if this is infinity.
   * @see Reference: sage/rings/complex_mpfr.pyx:is_infinity
   */
  is_infinity(): boolean {
    return !Number.isFinite(this._real) || !Number.isFinite(this._imag);
  }

  /**
   * Check if this is NaN.
   * @see Reference: sage/rings/complex_mpfr.pyx:is_NaN
   */
  is_NaN(): boolean {
    return Number.isNaN(this._real) || Number.isNaN(this._imag);
  }

  /**
   * Return the multiplicative order.
   * For complex numbers, only roots of unity have finite order.
   * @see Reference: sage/rings/complex_mpfr.pyx:multiplicative_order
   */
  multiplicative_order(): number {
    if (this._real === 1 && this._imag === 0) {
      return 1;
    }
    if (this._real === -1 && this._imag === 0) {
      return 2;
    }
    // self == C.gen() (= I) or self == -C.gen()
    if (this._real === 0 && this._imag === 1) {
      return 4;
    }
    if (this._real === 0 && this._imag === -1) {
      return 4;
    }
    // Clearly not a root of unity.
    if (Math.abs(this.abs() - 1) > 0.1) {
      return Number.POSITIVE_INFINITY;
    }
    // SageMath does NOT search for the order: a floating point number on the
    // unit circle carries no proof that it is a root of unity, so
    // e.g. ((1 + sqrt(-3))/2).multiplicative_order() raises.
    throw new NotImplementedError('order of element not known');
  }

  /**
   * Return the additive order: 1 for zero, `+Infinity` otherwise.
   *
   * ```
   * sage: CC(0).additive_order()      -> 1
   * sage: CC.gen().additive_order()   -> +Infinity
   * ```
   *
   * @see Reference: sage/rings/complex_mpfr.pyx:additive_order
   */
  additive_order(): number {
    if (this._real === 0 && this._imag === 0) {
      return 1;
    }
    return Number.POSITIVE_INFINITY;
  }

  /**
   * Return an algebraic dependency.
   * Uses LLL algorithm to find integer relations.
   * @see Reference: sage/rings/complex_mpfr.pyx:algebraic_dependency
   */
  algebraic_dependency(degree: number): bigint[] {
    // Integers: x - z.
    if (this._imag === 0 && Number.isInteger(this._real)) {
      const r = BigInt(Math.round(this._real));
      return [-r, 1n];
    }

    if (degree < 1) {
      throw new ValueError('degree must be at least 1');
    }

    // prec = z.prec() - 6, exactly as sage/arith/misc.py.
    // prec = z.prec() - 6, capped at the 53 bits a JavaScript double actually
    // carries: scaling by 2^prec with prec > 53 would only amplify noise.
    const prec = Math.min(this._parent.prec(), 53) - 6;
    const n = degree + 1;

    // M is n x (n + 2): the identity on the left, then the rounded imaginary
    // and real parts of r = 2^prec * z^k in the last two columns.
    const data: bigint[][] = [];
    for (let i = 0; i < n; i++) {
      data.push(new Array<bigint>(n + 2).fill(0n));
    }
    let rRe = 2 ** prec;
    let rIm = 0;
    data[0]![0] = 1n;
    data[0]![n + 1] = BigInt(Math.round(rRe));
    for (let k = 1; k <= degree; k++) {
      data[k]![k] = 1n;
      const newRe = rRe * this._real - rIm * this._imag;
      const newIm = rRe * this._imag + rIm * this._real;
      rRe = newRe;
      rIm = newIm;
      data[k]![n + 1] = BigInt(Math.round(rRe));
      data[k]![n] = BigInt(Math.round(rIm));
    }

    const reduced = LLL(new IntegerMatrix(n, n + 2, data), 0.75) as IntegerMatrix;

    let coeffs: bigint[] = [];
    for (let j = 0; j < n; j++) {
      coeffs.push(reduced.get(0, j).value);
    }
    // We're supposed to find an irreducible polynomial, so we cannot return a
    // constant one.  If the first LLL basis vector gives a constant polynomial,
    // use the next one.
    if (coeffs.slice(1).every((c) => c === 0n)) {
      coeffs = [];
      for (let j = 0; j < n; j++) {
        coeffs.push(reduced.get(1, j).value);
      }
    }

    if (coeffs[degree]! < 0n) {
      coeffs = coeffs.map((c) => -c);
    }

    // f might be reducible; return the best fitting irreducible factor.
    return this._bestIrreducibleFactor(coeffs);
  }

  /**
   * Given the integer coefficients of a polynomial `f`, return the irreducible
   * factor of `f` over `ZZ` minimising `|g(self)|`.
   *
   * Port of the last two lines of `sage/arith/misc.py:algebraic_dependency`:
   * `min((p for p, _ in R(f).factor()), key=lambda f: abs(f(z)))`.
   */
  private _bestIrreducibleFactor(coeffs: bigint[]): bigint[] {
    while (coeffs.length > 1 && coeffs[coeffs.length - 1] === 0n) {
      coeffs = coeffs.slice(0, -1);
    }
    if (coeffs.length <= 2) {
      return coeffs;
    }

    const ZZx = new PolynomialRing(ZZ_FOR_ALGDEP, 'x');
    const f = ZZx.__call__(coeffs.map((c) => new Integer(c)));
    const factors = f.factor().filter(([g]) => g.degree() > 0);
    if (factors.length === 0) {
      return coeffs;
    }

    let best: bigint[] | null = null;
    let bestValue = Number.POSITIVE_INFINITY;
    for (const [g] of factors) {
      const gCoeffs: bigint[] = g.coeffs.map((c) => c.value);
      const value = this._evaluateIntegerPolynomial(gCoeffs).abs();
      if (value < bestValue) {
        bestValue = value;
        best = gCoeffs;
      }
    }
    return best ?? coeffs;
  }

  /** Evaluate an integer polynomial (ascending coefficients) at ``self``. */
  private _evaluateIntegerPolynomial(coeffs: bigint[]): ComplexNumber {
    let result = new ComplexNumber(this._parent, 0, 0);
    for (let i = coeffs.length - 1; i >= 0; i--) {
      result = result.mul(this).add(new ComplexNumber(this._parent, Number(coeffs[i]!), 0));
    }
    return result;
  }

  toString(): string {
    if (this._imag === 0) {
      return this._real.toString();
    }
    if (this._real === 0) {
      return `${this._imag}*I`;
    }
    if (this._imag < 0) {
      return `${this._real} - ${-this._imag}*I`;
    }
    return `${this._real} + ${this._imag}*I`;
  }
}

/**
 * Create a complex field with given precision.
 * @see Reference: sage/rings/complex_mpfr.pyx:ComplexField
 */
export function CC(prec: number = 53): ComplexField {
  return new ComplexField(prec);
}

// Default complex field with 53 bits of precision
export const ComplexFieldDefault = new ComplexField(53);
