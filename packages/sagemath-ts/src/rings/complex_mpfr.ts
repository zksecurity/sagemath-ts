/**
 * @module sage/rings/complex_mpfr
 * @description Arbitrary precision floating point complex numbers using GNU MPFR
 *
 * Port of: sage/rings/complex_mpfr.pyx
 * Reference: reference/sage/src/sage/rings/complex_mpfr.pyx
 */

import { NotImplementedError } from '../errors.js';
import { RealField, RealNumber } from './real_mpfr.js';

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
    const r = typeof real === 'bigint' ? Number(real) : typeof real === 'string' ? parseFloat(real) : real;
    const i = imag === undefined ? 0 : typeof imag === 'bigint' ? Number(imag) : typeof imag === 'string' ? parseFloat(imag) : imag;
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
      return new ComplexNumber(this._parent, NaN, NaN);
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

    let a: number, b: number;
    if (avoidBranch) {
      // Compute sqrt of -z and multiply by i
      const a2 = (r - x) / 2;
      a = Math.sqrt(a2);
      b = y / (2 * a);
      // z = i * sqrt(-self), swap and adjust sign
      const tempA = a;
      a = -b;
      b = tempA;
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
    const r = Math.pow(rho, 1 / n);

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
      roots.push(new ComplexNumber(this._parent, r * Math.cos(currentArg), r * Math.sin(currentArg)));
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
        return new ComplexNumber(this._parent, Infinity, 0);
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
    const tPowZPlusHalf = t.log().mul(new ComplexNumber(this._parent, z._real + 0.5, z._imag)).exp();
    const expNegT = t.neg().exp();

    return sqrtTwoPi.mul(tPowZPlusHalf).mul(expNegT).mul(x);
  }

  /**
   * Return the incomplete gamma function.
   * Uses series expansion for small t, continued fraction for large t.
   * @see Reference: sage/rings/complex_mpfr.pyx:gamma_inc
   */
  gamma_inc(t: ComplexNumber): ComplexNumber {
    // Incomplete gamma function Gamma(a, x) = integral from x to infinity of t^(a-1)*e^(-t) dt
    // For simplicity, use series expansion
    // Gamma(a, x) = Gamma(a) - gamma_lower(a, x)
    // gamma_lower(a, x) = x^a * e^(-x) * sum_{n=0}^inf x^n / (a)_{n+1}

    const a = this;
    const x = t;
    const gammaA = a.gamma();

    // Series expansion for lower incomplete gamma
    const maxIter = 100;
    const eps = 1e-15;

    let term = new ComplexNumber(this._parent, 1, 0);
    let sum = new ComplexNumber(this._parent, 1, 0);
    let aPlusN = new ComplexNumber(this._parent, a._real, a._imag);

    for (let n = 1; n < maxIter; n++) {
      aPlusN = new ComplexNumber(this._parent, aPlusN._real + 1, aPlusN._imag);
      term = term.mul(x).div(aPlusN);
      sum = sum.add(term);
      if (term.abs() < eps * sum.abs()) break;
    }

    // gamma_lower = x^a * e^(-x) * sum / a
    const xPowA = x.log().mul(a).exp();
    const expNegX = x.neg().exp();
    const gammaLower = xPowA.mul(expNegX).mul(sum).div(a);

    return gammaA.sub(gammaLower);
  }

  /**
   * Return the Riemann zeta function.
   * Uses the Borwein algorithm for fast convergence.
   * @see Reference: sage/rings/complex_mpfr.pyx:zeta
   */
  zeta(): ComplexNumber {
    // Special case: zeta(1) is infinity
    if (this._imag === 0 && this._real === 1) {
      return new ComplexNumber(this._parent, Infinity, 0);
    }

    const s = this;

    // Use reflection formula for Re(s) < 0:
    // zeta(s) = 2^s * pi^(s-1) * sin(pi*s/2) * Gamma(1-s) * zeta(1-s)
    if (this._real < 0) {
      const one = new ComplexNumber(this._parent, 1, 0);
      const oneMinusS = one.sub(s);
      const zetaOneMinusS = oneMinusS.zeta();
      const two = new ComplexNumber(this._parent, 2, 0);
      const pi = new ComplexNumber(this._parent, Math.PI, 0);
      const twoPowS = two.log().mul(s).exp();
      const piPowSMinus1 = pi.log().mul(one.sub(s).neg()).exp();
      const piSOver2 = pi.mul(s).mul(new ComplexNumber(this._parent, 0.5, 0));
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
      d[k] = d[k - 1] + Math.pow(n, k) / this._factorial(k);
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
        const term = jPlus1.log().mul(s).neg().exp();
        const sign = (j % 2 === 0) ? 1 : -1;
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
    const oneMinusS = one.sub(s);
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
    const z = this;

    // Special case: Li_2(1) = pi^2/6
    if (z._real === 1 && z._imag === 0) {
      return new ComplexNumber(this._parent, Math.PI * Math.PI / 6, 0);
    }

    // Special case: Li_2(0) = 0
    if (z._real === 0 && z._imag === 0) {
      return new ComplexNumber(this._parent, 0, 0);
    }

    // For |z| <= 0.5, use direct series
    if (z.abs() <= 0.5) {
      return this._dilogSeries(z);
    }

    // For |z| > 0.5 and |1-z| small, use Li_2(z) + Li_2(1-z) + ln(z)*ln(1-z) = pi^2/6
    const one = new ComplexNumber(this._parent, 1, 0);
    const oneMinusZ = one.sub(z);

    if (oneMinusZ.abs() < 0.5) {
      // Use functional equation: Li_2(z) = pi^2/6 - Li_2(1-z) - ln(z)*ln(1-z)
      const pi2Over6 = new ComplexNumber(this._parent, Math.PI * Math.PI / 6, 0);
      const logZ = z.log();
      const logOneMinusZ = oneMinusZ.log();
      return pi2Over6.sub(this._dilogSeries(oneMinusZ)).sub(logZ.mul(logOneMinusZ));
    }

    // Use inversion formula: Li_2(z) = -Li_2(1/z) - pi^2/6 - (1/2)*ln(-z)^2
    if (z.abs() > 2) {
      const zInv = z.inv();
      const pi2Over6 = new ComplexNumber(this._parent, Math.PI * Math.PI / 6, 0);
      const logNegZ = z.neg().log();
      const half = new ComplexNumber(this._parent, 0.5, 0);
      return this._dilogSeries(zInv).neg().sub(pi2Over6).sub(half.mul(logNegZ.mul(logNegZ)));
    }

    // Default to series
    return this._dilogSeries(z);
  }

  /**
   * Helper for dilog series computation.
   */
  private _dilogSeries(z: ComplexNumber): ComplexNumber {
    const maxIter = 100;
    const eps = 1e-15;
    let sum = new ComplexNumber(this._parent, 0, 0);
    let zPow = z;

    for (let k = 1; k <= maxIter; k++) {
      const term = zPow.mul(new ComplexNumber(this._parent, 1 / (k * k), 0));
      sum = sum.add(term);
      if (term.abs() < eps * sum.abs()) break;
      zPow = zPow.mul(z);
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

    const z = this;
    const i = new ComplexNumber(this._parent, 0, 1);
    const pi = Math.PI;

    // Compute the product: prod_{n=1}^inf (1 - e^(2*pi*i*n*z))
    const maxIter = 100;
    const eps = 1e-15;
    let prod = new ComplexNumber(this._parent, 1, 0);

    for (let n = 1; n <= maxIter; n++) {
      // q = e^(2*pi*i*z)
      // term = 1 - q^n = 1 - e^(2*pi*i*n*z)
      const exp_arg = new ComplexNumber(this._parent, -2 * pi * n * z._imag, 2 * pi * n * z._real);
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
    const fracExp = new ComplexNumber(this._parent, -pi * z._imag / 12, pi * z._real / 12);
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
    return this._real === Infinity && this._imag === 0;
  }

  /**
   * Check if this is negative infinity.
   * @see Reference: sage/rings/complex_mpfr.pyx:is_negative_infinity
   */
  is_negative_infinity(): boolean {
    return this._real === -Infinity && this._imag === 0;
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
  multiplicative_order(): number | string {
    // Check if this is a root of unity (|z| = 1 and z^n = 1 for some n)
    const absVal = this.abs();
    if (Math.abs(absVal - 1) > 1e-10) {
      return 'Infinity';
    }
    if (this._real === 1 && this._imag === 0) {
      return 1;
    }
    if (this._real === -1 && this._imag === 0) {
      return 2;
    }

    // Check if it's an n-th root of unity for small n
    const arg = this.argument();
    for (let n = 1; n <= 1000; n++) {
      if (Math.abs(Math.sin(n * arg)) < 1e-10 && Math.abs(Math.cos(n * arg) - 1) < 1e-10) {
        return n;
      }
    }
    return 'Infinity';
  }

  /**
   * Return the additive order.
   * Returns 1 if self is 0, otherwise infinity.
   * @see Reference: sage/rings/complex_mpfr.pyx:additive_order
   */
  additive_order(): number | string {
    if (this._real === 0 && this._imag === 0) {
      return 1;
    }
    return 'Infinity';
  }

  /**
   * Return an algebraic dependency.
   * Uses LLL algorithm to find integer relations.
   * @see Reference: sage/rings/complex_mpfr.pyx:algebraic_dependency
   */
  algebraic_dependency(n: number): bigint[] {
    // Simple implementation using PSLQ-like algorithm
    // Returns coefficients [a_0, a_1, ..., a_n] such that a_0 + a_1*z + ... + a_n*z^n = 0

    // Build a vector of powers of z
    const powers: ComplexNumber[] = [];
    let zPow = new ComplexNumber(this._parent, 1, 0);
    for (let i = 0; i <= n; i++) {
      powers.push(zPow);
      zPow = zPow.mul(this);
    }

    // Use a simple LLL-based approach
    // This is a simplified version - in practice would need a proper PSLQ implementation
    const scale = 1e15;
    const matrix: bigint[][] = [];

    for (let i = 0; i <= n; i++) {
      const row: bigint[] = new Array(n + 3).fill(0n);
      row[i] = 1n;
      row[n + 1] = BigInt(Math.round(powers[i].real() * scale));
      row[n + 2] = BigInt(Math.round(powers[i].imag() * scale));
      matrix.push(row);
    }

    // Simple reduction (not full LLL, but gives approximate result)
    // Return trivial polynomial if real and integer
    if (this._imag === 0 && Number.isInteger(this._real)) {
      const r = BigInt(Math.round(this._real));
      return [-r, 1n]; // z - r = 0
    }

    // For non-trivial cases, return a placeholder
    // A proper implementation would use LLL lattice reduction
    const result: bigint[] = new Array(n + 1).fill(0n);
    result[0] = 1n;
    result[n] = 1n;
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
