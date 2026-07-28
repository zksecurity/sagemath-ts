/**
 * @module sage/rings/number_field/number_field
 * @description Number fields (algebraic number fields)
 *
 * Port of: sage/rings/number_field/number_field.py
 * Reference: reference/sage/src/sage/rings/number_field/number_field.py
 *
 * NOTE: SageMath delegates most number field operations to PARI/GP (nfinit, bnfinit).
 * This implementation provides a pure TypeScript version with core functionality.
 * For full PARI compatibility, see DEVIATIONS.md.
 *
 * @see Deviation: Number Field Implementation Without PARI
 */

import { gcd as intGcd, lcm as intLcm, is_prime, isqrt } from '../../arith/misc.js';
import { NotImplementedError, ValueError, ZeroDivisionError } from '../../errors.js';
import { Rational } from '../rational.js';
import type { ClassGroup } from './class_group.js';
import type { GaloisGroup } from './galois_group.js';
import type { NumberFieldIdeal } from './number_field_ideal.js';
import type { AbsoluteOrder, Order } from './order.js';
import {
  type MulTable,
  type NfBasisResult,
  type PrimeDecEntry,
  type ZPoly,
  fpFactor,
  hnf,
  integralDefiningPolynomial,
  nfbasis,
  nfgaloisconj,
  primedec,
  ratInverse,
  zpIsIrreducibleOverQ,
} from './pari_nf.js';
import type { UnitGroup } from './unit_group.js';

/**
 * A polynomial over the rationals, represented as an array of Rational coefficients.
 * coeffs[i] is the coefficient of x^i.
 */
export class RationalPolynomial {
  readonly coeffs: readonly Rational[];

  constructor(coeffs: Rational[]) {
    // Remove trailing zeros
    let len = coeffs.length;
    while (len > 0 && coeffs[len - 1]!.isZero()) {
      len--;
    }
    this.coeffs = Object.freeze(coeffs.slice(0, len));
  }

  /**
   * Create a polynomial from bigint coefficients.
   */
  static fromBigInts(coeffs: bigint[]): RationalPolynomial {
    return new RationalPolynomial(coeffs.map((c) => new Rational(c)));
  }

  /**
   * Create the zero polynomial.
   */
  static zero(): RationalPolynomial {
    return new RationalPolynomial([]);
  }

  /**
   * Create a constant polynomial.
   */
  static constant(c: Rational): RationalPolynomial {
    if (c.isZero()) return RationalPolynomial.zero();
    return new RationalPolynomial([c]);
  }

  /**
   * Create the polynomial x.
   */
  static x(): RationalPolynomial {
    return new RationalPolynomial([Rational.zero(), Rational.one()]);
  }

  /**
   * Return the degree. Zero polynomial has degree -1.
   */
  degree(): number {
    return this.coeffs.length - 1;
  }

  /**
   * Check if this is the zero polynomial.
   */
  isZero(): boolean {
    return this.coeffs.length === 0;
  }

  /**
   * Check if this is a constant polynomial.
   */
  isConstant(): boolean {
    return this.coeffs.length <= 1;
  }

  /**
   * Return the leading coefficient.
   */
  leadingCoefficient(): Rational {
    if (this.coeffs.length === 0) return Rational.zero();
    return this.coeffs[this.coeffs.length - 1]!;
  }

  /**
   * Return coefficient of x^n.
   */
  getCoeff(n: number): Rational {
    if (n < 0 || n >= this.coeffs.length) return Rational.zero();
    return this.coeffs[n]!;
  }

  /**
   * Check if monic (leading coeff = 1).
   */
  isMonic(): boolean {
    if (this.coeffs.length === 0) return false;
    return this.leadingCoefficient().eq(Rational.one());
  }

  /**
   * Return the monic version of this polynomial.
   */
  monic(): RationalPolynomial {
    if (this.isZero()) {
      throw new ZeroDivisionError('cannot make zero polynomial monic');
    }
    const lc = this.leadingCoefficient();
    if (lc.eq(Rational.one())) return this;
    return new RationalPolynomial(this.coeffs.map((c) => c.div(lc)));
  }

  /**
   * Add two polynomials.
   */
  add(other: RationalPolynomial): RationalPolynomial {
    const maxLen = Math.max(this.coeffs.length, other.coeffs.length);
    const result: Rational[] = [];
    for (let i = 0; i < maxLen; i++) {
      const a = this.getCoeff(i);
      const b = other.getCoeff(i);
      result.push(a.add(b));
    }
    return new RationalPolynomial(result);
  }

  /**
   * Subtract two polynomials.
   */
  sub(other: RationalPolynomial): RationalPolynomial {
    const maxLen = Math.max(this.coeffs.length, other.coeffs.length);
    const result: Rational[] = [];
    for (let i = 0; i < maxLen; i++) {
      const a = this.getCoeff(i);
      const b = other.getCoeff(i);
      result.push(a.sub(b));
    }
    return new RationalPolynomial(result);
  }

  /**
   * Negate the polynomial.
   */
  neg(): RationalPolynomial {
    return new RationalPolynomial(this.coeffs.map((c) => c.neg()));
  }

  /**
   * Multiply two polynomials.
   */
  mul(other: RationalPolynomial): RationalPolynomial {
    if (this.isZero() || other.isZero()) {
      return RationalPolynomial.zero();
    }
    const n = this.coeffs.length;
    const m = other.coeffs.length;
    const result: Rational[] = Array(n + m - 1)
      .fill(null)
      .map(() => Rational.zero());

    for (let i = 0; i < n; i++) {
      for (let j = 0; j < m; j++) {
        result[i + j] = result[i + j]!.add(this.coeffs[i]!.mul(other.coeffs[j]!));
      }
    }
    return new RationalPolynomial(result);
  }

  /**
   * Multiply by a scalar.
   */
  scale(c: Rational): RationalPolynomial {
    if (c.isZero()) return RationalPolynomial.zero();
    return new RationalPolynomial(this.coeffs.map((x) => x.mul(c)));
  }

  /**
   * Divide with remainder: returns [quotient, remainder] such that
   * this = quotient * other + remainder with deg(remainder) < deg(other).
   */
  divmod(other: RationalPolynomial): [RationalPolynomial, RationalPolynomial] {
    if (other.isZero()) {
      throw new ZeroDivisionError('polynomial division by zero');
    }

    if (this.degree() < other.degree()) {
      return [RationalPolynomial.zero(), this];
    }

    let remainder = new RationalPolynomial([...this.coeffs]);
    const quotientCoeffs: Rational[] = Array(this.degree() - other.degree() + 1)
      .fill(null)
      .map(() => Rational.zero());

    const divisorLC = other.leadingCoefficient();
    const divisorDeg = other.degree();

    while (!remainder.isZero() && remainder.degree() >= divisorDeg) {
      const shift = remainder.degree() - divisorDeg;
      const coeff = remainder.leadingCoefficient().div(divisorLC);
      quotientCoeffs[shift] = coeff;

      // Subtract coeff * x^shift * other from remainder
      const subtractCoeffs: Rational[] = Array(shift)
        .fill(null)
        .map(() => Rational.zero());
      for (let i = 0; i < other.coeffs.length; i++) {
        subtractCoeffs.push(other.coeffs[i]!.mul(coeff));
      }
      remainder = remainder.sub(new RationalPolynomial(subtractCoeffs));
    }

    return [new RationalPolynomial(quotientCoeffs), remainder];
  }

  /**
   * Return this polynomial modulo other.
   */
  mod(other: RationalPolynomial): RationalPolynomial {
    return this.divmod(other)[1];
  }

  /**
   * Evaluate at a Rational value.
   */
  evaluate(x: Rational): Rational {
    if (this.isZero()) return Rational.zero();
    // Horner's method
    let result = this.coeffs[this.coeffs.length - 1]!;
    for (let i = this.coeffs.length - 2; i >= 0; i--) {
      result = result.mul(x).add(this.coeffs[i]!);
    }
    return result;
  }

  /**
   * Return the derivative.
   */
  derivative(): RationalPolynomial {
    if (this.coeffs.length <= 1) return RationalPolynomial.zero();
    const result: Rational[] = [];
    for (let i = 1; i < this.coeffs.length; i++) {
      result.push(this.coeffs[i]!.mul(new Rational(BigInt(i))));
    }
    return new RationalPolynomial(result);
  }

  /**
   * Return GCD of two polynomials using Euclidean algorithm.
   */
  gcd(other: RationalPolynomial): RationalPolynomial {
    let a: RationalPolynomial = this;
    let b: RationalPolynomial = other;
    while (!b.isZero()) {
      const r = a.mod(b);
      a = b;
      b = r;
    }
    if (a.isZero()) return a;
    return a.monic();
  }

  /**
   * Return the integer polynomial obtained by clearing denominators.
   * The result is primitive, with the same roots as this polynomial.
   */
  integerCoefficients(): bigint[] {
    let denomLcm = 1n;
    for (const c of this.coeffs) {
      denomLcm = intLcm(denomLcm, c.denominator);
    }
    const out = this.coeffs.map((c) => c.numerator * (denomLcm / c.denominator));
    let g = 0n;
    for (const c of out) g = intGcd(g, c < 0n ? -c : c);
    if (g > 1n) return out.map((c) => c / g);
    return out;
  }

  /**
   * Check if this polynomial is irreducible over Q.
   *
   * SageMath delegates to PARI's `polisirreducible`; this runs the same
   * algorithm (squarefree test followed by Zassenhaus factorisation).
   *
   * @see Reference: sage/rings/polynomial/polynomial_element.pyx:is_irreducible
   */
  isIrreducible(): boolean {
    if (this.degree() <= 0) return false;
    if (this.degree() === 1) return true;
    return zpIsIrreducibleOverQ(this.integerCoefficients());
  }

  /**
   * Return the content (GCD of all coefficients).
   */
  content(): Rational {
    if (this.isZero()) return Rational.zero();

    // Compute GCD of all numerators and LCM of all denominators
    let numerGcd = 0n;
    let denomLcm = 1n;

    for (const c of this.coeffs) {
      const n = c.numerator < 0n ? -c.numerator : c.numerator;
      const d = c.denominator;
      numerGcd = numerGcd === 0n ? n : intGcd(numerGcd, n);
      denomLcm = (denomLcm * d) / intGcd(denomLcm, d);
    }

    return new Rational(numerGcd, denomLcm);
  }

  /**
   * Return the primitive part (this / content).
   */
  primitivePart(): RationalPolynomial {
    const c = this.content();
    if (c.isZero()) return this;
    return this.scale(c.inv());
  }

  /**
   * Check equality.
   */
  eq(other: RationalPolynomial): boolean {
    if (this.coeffs.length !== other.coeffs.length) return false;
    for (let i = 0; i < this.coeffs.length; i++) {
      if (!this.coeffs[i]!.eq(other.coeffs[i]!)) return false;
    }
    return true;
  }

  /**
   * Return the discriminant of this polynomial.
   * disc(f) = (-1)^(n(n-1)/2) * Res(f, f') / a_n
   * where n = deg(f) and a_n is the leading coefficient.
   */
  discriminant(): Rational {
    const n = this.degree();
    if (n < 1) {
      throw new ValueError('discriminant undefined for constant polynomials');
    }
    if (n === 1) {
      return Rational.one();
    }

    const fp = this.derivative();
    const res = this.resultant(fp);
    const an = this.leadingCoefficient();
    const sign = ((n * (n - 1)) / 2) % 2 === 0 ? 1n : -1n;

    return res.mul(new Rational(sign)).div(an);
  }

  /**
   * Compute the resultant of this polynomial and other.
   * Uses the Sylvester matrix determinant.
   */
  resultant(other: RationalPolynomial): Rational {
    const n = this.degree();
    const m = other.degree();

    if (n < 0 || m < 0) {
      throw new ValueError('resultant undefined for zero polynomials');
    }

    // Build Sylvester matrix of size (n + m) x (n + m)
    const size = n + m;
    const matrix: Rational[][] = Array(size)
      .fill(null)
      .map(() =>
        Array(size)
          .fill(null)
          .map(() => Rational.zero())
      );

    // First m rows: coefficients of this polynomial
    for (let i = 0; i < m; i++) {
      for (let j = 0; j <= n; j++) {
        matrix[i]![i + j] = this.getCoeff(n - j);
      }
    }

    // Next n rows: coefficients of other polynomial
    for (let i = 0; i < n; i++) {
      for (let j = 0; j <= m; j++) {
        matrix[m + i]![i + j] = other.getCoeff(m - j);
      }
    }

    // Compute determinant using Gaussian elimination
    return matrixDeterminant(matrix);
  }

  toString(): string {
    if (this.isZero()) return '0';

    const terms: string[] = [];
    for (let i = this.coeffs.length - 1; i >= 0; i--) {
      const c = this.coeffs[i]!;
      if (c.isZero()) continue;

      let term = '';
      if (i === 0) {
        term = c.toString();
      } else if (i === 1) {
        if (c.eq(Rational.one())) {
          term = 'x';
        } else if (c.eq(new Rational(-1n))) {
          term = '-x';
        } else {
          term = `${c}*x`;
        }
      } else {
        if (c.eq(Rational.one())) {
          term = `x^${i}`;
        } else if (c.eq(new Rational(-1n))) {
          term = `-x^${i}`;
        } else {
          term = `${c}*x^${i}`;
        }
      }

      terms.push(term);
    }

    return terms.length > 0 ? terms.join(' + ').replace(/\+ -/g, '- ') : '0';
  }
}

/**
 * Helper: compute divisors of n.
 */
function divisors(n: bigint): bigint[] {
  if (n <= 0n) return [1n];
  const result: bigint[] = [];
  let i = 1n;
  while (i * i <= n) {
    if (n % i === 0n) {
      result.push(i);
      if (i !== n / i) {
        result.push(n / i);
      }
    }
    i++;
  }
  return result.sort((a, b) => (a < b ? -1 : a > b ? 1 : 0));
}

/**
 * Helper: compute determinant of a rational matrix using Gaussian elimination.
 */
function matrixDeterminant(m: Rational[][]): Rational {
  const n = m.length;
  if (n === 0) return Rational.one();

  // Make a copy
  const matrix = m.map((row) => [...row]);

  let det = Rational.one();
  let sign = 1;

  for (let col = 0; col < n; col++) {
    // Find pivot
    let pivotRow = -1;
    for (let row = col; row < n; row++) {
      if (!matrix[row]![col]!.isZero()) {
        pivotRow = row;
        break;
      }
    }

    if (pivotRow === -1) {
      return Rational.zero(); // Singular matrix
    }

    // Swap rows if needed
    if (pivotRow !== col) {
      [matrix[col], matrix[pivotRow]] = [matrix[pivotRow]!, matrix[col]!];
      sign = -sign;
    }

    const pivot = matrix[col]![col]!;
    det = det.mul(pivot);

    // Eliminate column below pivot
    for (let row = col + 1; row < n; row++) {
      const factor = matrix[row]![col]!.div(pivot);
      for (let j = col; j < n; j++) {
        matrix[row]![j] = matrix[row]![j]!.sub(factor.mul(matrix[col]![j]!));
      }
    }
  }

  return sign === 1 ? det : det.neg();
}

/**
 * An element of a number field.
 *
 * Represented as a polynomial in the generator (alpha) modulo the defining polynomial.
 */
export class NumberFieldElement {
  private readonly _parent: NumberField;
  private readonly _coeffs: readonly Rational[]; // Coefficients in power basis [a_0, a_1, ..., a_{n-1}]

  constructor(parent: NumberField, coeffs: Rational[]) {
    this._parent = parent;
    const n = parent.degree();
    // Ensure coeffs has exactly n elements (pad with zeros or truncate)
    const normalizedCoeffs = Array(n)
      .fill(null)
      .map((_, i) => (i < coeffs.length ? coeffs[i]! : Rational.zero()));
    this._coeffs = Object.freeze(normalizedCoeffs);
  }

  /**
   * Create from a single rational (constant element).
   */
  static fromRational(parent: NumberField, r: Rational): NumberFieldElement {
    return new NumberFieldElement(parent, [r]);
  }

  /**
   * Create from a bigint.
   */
  static fromBigInt(parent: NumberField, n: bigint): NumberFieldElement {
    return NumberFieldElement.fromRational(parent, new Rational(n));
  }

  /**
   * Return the polynomial representation as coefficients.
   */
  polynomial(): RationalPolynomial {
    return new RationalPolynomial([...this._coeffs]);
  }

  /**
   * Return the parent number field.
   */
  parent(): NumberField {
    return this._parent;
  }

  /**
   * Return the coefficients in the power basis.
   */
  list(): Rational[] {
    return [...this._coeffs];
  }

  /**
   * Return the i-th coefficient.
   */
  __getitem__(i: number): Rational {
    if (i < 0 || i >= this._coeffs.length) return Rational.zero();
    return this._coeffs[i]!;
  }

  /**
   * Add two elements.
   */
  add(other: NumberFieldElement): NumberFieldElement {
    this._checkParent(other);
    const result = this._coeffs.map((c, i) => c.add(other._coeffs[i]!));
    return new NumberFieldElement(this._parent, result);
  }

  /**
   * Subtract two elements.
   */
  sub(other: NumberFieldElement): NumberFieldElement {
    this._checkParent(other);
    const result = this._coeffs.map((c, i) => c.sub(other._coeffs[i]!));
    return new NumberFieldElement(this._parent, result);
  }

  /**
   * Negate the element.
   */
  neg(): NumberFieldElement {
    return new NumberFieldElement(
      this._parent,
      this._coeffs.map((c) => c.neg())
    );
  }

  /**
   * Multiply two elements.
   */
  mul(other: NumberFieldElement): NumberFieldElement {
    this._checkParent(other);
    const n = this._parent.degree();
    const poly = this._parent.polynomial();

    // Multiply as polynomials
    const productPoly = this.polynomial().mul(other.polynomial());

    // Reduce modulo the defining polynomial
    const reduced = productPoly.mod(poly);

    // Convert back to coefficients
    const coeffs = Array(n)
      .fill(null)
      .map((_, i) => reduced.getCoeff(i));
    return new NumberFieldElement(this._parent, coeffs);
  }

  /**
   * Multiply by a rational.
   */
  scalarMul(r: Rational): NumberFieldElement {
    return new NumberFieldElement(
      this._parent,
      this._coeffs.map((c) => c.mul(r))
    );
  }

  /**
   * Divide by another element.
   */
  div(other: NumberFieldElement): NumberFieldElement {
    return this.mul(other.inv());
  }

  /**
   * Return the multiplicative inverse using extended Euclidean algorithm.
   */
  inv(): NumberFieldElement {
    if (this.is_zero()) {
      throw new ZeroDivisionError('cannot invert zero');
    }

    const poly = this._parent.polynomial();
    const thisPoly = this.polynomial();

    // Extended GCD: find s, t such that s * thisPoly + t * poly = gcd
    // Since poly is irreducible and thisPoly is non-zero, gcd = 1
    const [_, s] = extendedGcd(thisPoly, poly);

    // s is the inverse modulo poly
    const reduced = s.mod(poly);
    const n = this._parent.degree();
    const coeffs = Array(n)
      .fill(null)
      .map((_, i) => reduced.getCoeff(i));
    return new NumberFieldElement(this._parent, coeffs);
  }

  /**
   * Return self raised to power n.
   */
  pow(n: bigint): NumberFieldElement {
    if (n === 0n) {
      return this._parent.one();
    }

    if (n < 0n) {
      return this.inv().pow(-n);
    }

    // Binary exponentiation
    let result = this._parent.one();
    let base: NumberFieldElement = this;

    while (n > 0n) {
      if (n % 2n === 1n) {
        result = result.mul(base);
      }
      base = base.mul(base);
      n = n / 2n;
    }

    return result;
  }

  /**
   * Check equality.
   */
  eq(other: NumberFieldElement): boolean {
    if (this._parent !== other._parent) return false;
    for (let i = 0; i < this._coeffs.length; i++) {
      if (!this._coeffs[i]!.eq(other._coeffs[i]!)) return false;
    }
    return true;
  }

  /**
   * Check if this is zero.
   */
  is_zero(): boolean {
    return this._coeffs.every((c) => c.isZero());
  }

  /**
   * Check if this is one.
   */
  is_one(): boolean {
    if (this._coeffs.length === 0) return false;
    if (!this._coeffs[0]!.eq(Rational.one())) return false;
    for (let i = 1; i < this._coeffs.length; i++) {
      if (!this._coeffs[i]!.isZero()) return false;
    }
    return true;
  }

  /**
   * Return the trace of this element.
   * trace(a) = sum of all conjugates = -coefficient of x^{n-1} in charpoly.
   */
  trace(): Rational {
    const cp = this.charpoly();
    const n = this._parent.degree();
    // For monic polynomial x^n + a_{n-1}x^{n-1} + ..., trace = -a_{n-1}
    return cp.getCoeff(n - 1).neg();
  }

  /**
   * Alias for trace.
   */
  absolute_trace(): Rational {
    return this.trace();
  }

  /**
   * Return the relative trace (for relative extensions - same as trace for absolute).
   */
  relative_trace(): Rational {
    return this.trace();
  }

  /**
   * Return the norm of this element.
   * norm(a) = product of all conjugates = (-1)^n * constant term of charpoly.
   */
  norm(): Rational {
    const cp = this.charpoly();
    const n = this._parent.degree();
    const sign = n % 2 === 0 ? 1n : -1n;
    return cp.getCoeff(0).mul(new Rational(sign));
  }

  /**
   * Return the absolute norm.
   */
  absolute_norm(): Rational {
    return this.norm();
  }

  /**
   * Return the relative norm (for relative extensions - same as norm for absolute).
   */
  relative_norm(): Rational {
    return this.norm();
  }

  /**
   * Return the characteristic polynomial of this element.
   *
   * The characteristic polynomial is the determinant of (xI - M) where M is the
   * matrix representing multiplication by this element.
   */
  charpoly(): RationalPolynomial {
    const n = this._parent.degree();

    // Build the multiplication matrix
    const matrix = this._multiplicationMatrix();

    // Compute det(xI - M) using Faddeev-LeVerrier algorithm or direct expansion
    // For simplicity, use direct computation for small n
    return characteristicPolynomial(matrix);
  }

  /**
   * Return the minimal polynomial of this element.
   */
  minpoly(): RationalPolynomial {
    const cp = this.charpoly();

    // The minimal polynomial divides the characteristic polynomial
    // For a field element, it's the monic polynomial of smallest degree
    // that annihilates the element

    // Start with the characteristic polynomial and try to find a proper divisor
    // that still annihilates the element
    let mp = cp;

    // Factor out repeated factors by computing GCD with derivative
    let prev = mp;
    while (true) {
      const g = mp.gcd(mp.derivative());
      if (g.degree() <= 0) break;
      const [q, r] = mp.divmod(g);
      if (!r.isZero()) break;
      mp = q;
      if (mp.eq(prev)) break;
      prev = mp;
    }

    return mp.monic();
  }

  /**
   * Alias for minpoly.
   */
  minimal_polynomial(): RationalPolynomial {
    return this.minpoly();
  }

  /**
   * Return the matrix representing multiplication by this element.
   */
  private _multiplicationMatrix(): Rational[][] {
    const n = this._parent.degree();
    const matrix: Rational[][] = Array(n)
      .fill(null)
      .map(() =>
        Array(n)
          .fill(null)
          .map(() => Rational.zero())
      );

    // For each basis element alpha^i, compute this * alpha^i
    const alpha = this._parent.gen();
    let alphaPower = this._parent.one();

    for (let i = 0; i < n; i++) {
      const product = this.mul(alphaPower);
      for (let j = 0; j < n; j++) {
        matrix[j]![i] = product._coeffs[j]!;
      }
      alphaPower = alphaPower.mul(alpha);
    }

    return matrix;
  }

  /**
   * Check if this element is integral (in the ring of integers).
   * An element is integral if its minimal polynomial has integer coefficients.
   */
  is_integral(): boolean {
    const mp = this.minpoly();
    // Check if all coefficients are integers
    for (const c of mp.coeffs) {
      if (c.denominator !== 1n) return false;
    }
    return true;
  }

  /**
   * Check whether this element is a unit *in the ring where it is defined*.
   *
   * The parent of a `NumberFieldElement` is a field, so every nonzero element
   * is a unit; `is_integral_unit()` is the test for the ring of integers.
   *
   * @see Reference: sage/rings/number_field/number_field_element.pyx:1592
   *   (`if self.parent().is_field(): return bool(self)`)
   */
  is_unit(): boolean {
    return !this.is_zero();
  }

  /**
   * Check whether this element is a unit of the ring of integers, i.e. is an
   * algebraic integer of norm +/-1.  This is what
   * `O_K(x).is_unit()` returns in SageMath.
   */
  is_integral_unit(): boolean {
    if (!this.is_integral()) return false;
    const n = this.norm();
    return n.eq(Rational.one()) || n.eq(new Rational(-1n));
  }

  /**
   * Return the denominator (LCM of coefficient denominators).
   */
  denominator(): bigint {
    let lcm = 1n;
    for (const c of this._coeffs) {
      lcm = (lcm * c.denominator) / intGcd(lcm, c.denominator);
    }
    return lcm;
  }

  /**
   * Return the numerator (this * denominator).
   */
  numerator(): NumberFieldElement {
    const d = this.denominator();
    return this.scalarMul(new Rational(d));
  }

  private _checkParent(other: NumberFieldElement): void {
    if (this._parent !== other._parent) {
      throw new ValueError('elements must be in the same number field');
    }
  }

  toString(): string {
    const alpha = this._parent._name;
    const terms: string[] = [];

    for (let i = this._coeffs.length - 1; i >= 0; i--) {
      const c = this._coeffs[i]!;
      if (c.isZero()) continue;

      let term = '';
      if (i === 0) {
        term = c.toString();
      } else if (i === 1) {
        if (c.eq(Rational.one())) {
          term = alpha;
        } else if (c.eq(new Rational(-1n))) {
          term = `-${alpha}`;
        } else {
          term = `${c}*${alpha}`;
        }
      } else {
        if (c.eq(Rational.one())) {
          term = `${alpha}^${i}`;
        } else if (c.eq(new Rational(-1n))) {
          term = `-${alpha}^${i}`;
        } else {
          term = `${c}*${alpha}^${i}`;
        }
      }

      terms.push(term);
    }

    return terms.length > 0 ? terms.join(' + ').replace(/\+ -/g, '- ') : '0';
  }
}

/**
 * Extended GCD for polynomials: returns [gcd, s, t] such that s*a + t*b = gcd.
 */
function extendedGcd(
  a: RationalPolynomial,
  b: RationalPolynomial
): [RationalPolynomial, RationalPolynomial, RationalPolynomial] {
  let old_r = a;
  let r = b;
  let old_s = RationalPolynomial.constant(Rational.one());
  let s = RationalPolynomial.zero();
  let old_t = RationalPolynomial.zero();
  let t = RationalPolynomial.constant(Rational.one());

  while (!r.isZero()) {
    const [q] = old_r.divmod(r);
    [old_r, r] = [r, old_r.sub(q.mul(r))];
    [old_s, s] = [s, old_s.sub(q.mul(s))];
    [old_t, t] = [t, old_t.sub(q.mul(t))];
  }

  // Make monic
  if (!old_r.isZero()) {
    const lc = old_r.leadingCoefficient();
    old_r = old_r.scale(lc.inv());
    old_s = old_s.scale(lc.inv());
    old_t = old_t.scale(lc.inv());
  }

  return [old_r, old_s, old_t];
}

/**
 * Compute characteristic polynomial of a matrix using Faddeev-LeVerrier algorithm.
 */
function characteristicPolynomial(m: Rational[][]): RationalPolynomial {
  const n = m.length;
  if (n === 0) return RationalPolynomial.constant(Rational.one());

  // For small matrices, use direct formulas
  if (n === 1) {
    // det(xI - M) = x - m[0][0]
    return new RationalPolynomial([m[0]![0]!.neg(), Rational.one()]);
  }

  if (n === 2) {
    // det(xI - M) = x^2 - tr(M)*x + det(M)
    const tr = m[0]![0]!.add(m[1]![1]!);
    const det = m[0]![0]!.mul(m[1]![1]!).sub(m[0]![1]!.mul(m[1]![0]!));
    return new RationalPolynomial([det, tr.neg(), Rational.one()]);
  }

  if (n === 3) {
    // det(xI - M) = x^3 - tr(M)*x^2 + (sum of 2x2 minors)*x - det(M)
    const a = m[0]![0]!;
    const b = m[0]![1]!;
    const c = m[0]![2]!;
    const d = m[1]![0]!;
    const e = m[1]![1]!;
    const f = m[1]![2]!;
    const g = m[2]![0]!;
    const h = m[2]![1]!;
    const i = m[2]![2]!;

    const tr = a.add(e).add(i);

    // Sum of principal 2x2 minors
    const m11 = e.mul(i).sub(f.mul(h)); // minor removing row/col 0
    const m22 = a.mul(i).sub(c.mul(g)); // minor removing row/col 1
    const m33 = a.mul(e).sub(b.mul(d)); // minor removing row/col 2
    const sumMinors = m11.add(m22).add(m33);

    // 3x3 determinant
    const det = a
      .mul(e.mul(i).sub(f.mul(h)))
      .sub(b.mul(d.mul(i).sub(f.mul(g))))
      .add(c.mul(d.mul(h).sub(e.mul(g))));

    return new RationalPolynomial([det.neg(), sumMinors, tr.neg(), Rational.one()]);
  }

  // For larger matrices, use Faddeev-LeVerrier algorithm
  // B_k = M^k - c_{n-1}*M^{k-1} - ... - c_{n-k+1}*I
  // c_{n-k} = tr(M * B_{k-1}) / k

  const coeffs: Rational[] = Array(n + 1)
    .fill(null)
    .map(() => Rational.zero());
  coeffs[n] = Rational.one();

  // B_0 = I (identity)
  let B: Rational[][] = Array(n)
    .fill(null)
    .map((_, i) =>
      Array(n)
        .fill(null)
        .map((_, j) => (i === j ? Rational.one() : Rational.zero()))
    );

  for (let k = 1; k <= n; k++) {
    // Compute M * B
    const MB: Rational[][] = Array(n)
      .fill(null)
      .map(() =>
        Array(n)
          .fill(null)
          .map(() => Rational.zero())
      );

    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        for (let l = 0; l < n; l++) {
          MB[i]![j] = MB[i]![j]!.add(m[i]![l]!.mul(B[l]![j]!));
        }
      }
    }

    // c_{n-k} = -tr(M*B) / k
    let tr = Rational.zero();
    for (let i = 0; i < n; i++) {
      tr = tr.add(MB[i]![i]!);
    }
    coeffs[n - k] = tr.neg().div(new Rational(BigInt(k)));

    // B = M*B + c_{n-k} * I
    if (k < n) {
      for (let i = 0; i < n; i++) {
        MB[i]![i] = MB[i]![i]!.add(coeffs[n - k]!);
      }
      B = MB;
    }
  }

  return new RationalPolynomial(coeffs);
}

/**
 * A number field (a finite extension of the rational numbers).
 *
 * A number field can be created by specifying a polynomial over the rationals
 * and a name for the generator of the field.
 *
 * @example
 * ```typescript
 * // Create Q(sqrt(2))
 * const poly = RationalPolynomial.fromBigInts([-2n, 0n, 1n]); // x^2 - 2
 * const K = new NumberField(poly, 'a');
 * const a = K.gen();
 * const two = a.mul(a); // a^2 = 2
 * ```
 *
 * @see Reference: sage/rings/number_field/number_field.py:NumberField_generic
 */
export class NumberField {
  private readonly _polynomial: RationalPolynomial;
  readonly _name: string;
  private readonly _embedding?: unknown;
  private _cachedGen?: NumberFieldElement;
  private _cachedSignature?: [number, number];
  private _cachedIntegralData?: { nf: NfBasisResult; scale: bigint };
  private _cachedIntegralBasis?: NumberFieldElement[];
  private _cachedPariBasis?: NumberFieldElement[];
  private _cachedMulTable?: MulTable;

  constructor(polynomial: RationalPolynomial, name: string, embedding?: unknown, check = true) {
    // Validate polynomial
    if (polynomial.degree() < 1) {
      throw new ValueError('defining polynomial must have degree at least 1');
    }

    if (check && !polynomial.isIrreducible()) {
      throw new ValueError(`defining polynomial (${polynomial}) must be irreducible`);
    }

    // Store the monic version
    this._polynomial = polynomial.monic();
    this._name = name;
    this._embedding = embedding;
  }

  /**
   * The integral data of this field: a monic *integral* polynomial `g` with
   * root `theta = scale * alpha`, an integral basis of `O_K` in the power
   * basis of `theta`, the index `[O_K : Z[theta]]` and `disc(O_K)`.
   *
   * This is the port's stand-in for Sage's `K.pari_nf()`.
   *
   * @see Reference: sage/rings/number_field/number_field.py:pari_nf
   */
  protected _integralData(): { nf: NfBasisResult; scale: bigint } {
    if (this._cachedIntegralData) return this._cachedIntegralData;
    const { g, scale } = integralDefiningPolynomial([...this._polynomial.coeffs]);
    this._cachedIntegralData = { nf: nfbasis(g), scale };
    return this._cachedIntegralData;
  }

  /**
   * The monic integral polynomial defining this field, as an array of
   * coefficients (low degree first).  Its root is `scale * alpha`.
   *
   * @see Reference: sage/rings/number_field/number_field.py:pari_polynomial
   */
  pari_polynomial(): ZPoly {
    return [...this._integralData().nf.g];
  }

  /**
   * Return the defining polynomial of this number field.
   */
  defining_polynomial(): RationalPolynomial {
    return this._polynomial;
  }

  /**
   * Return the minimal polynomial of this number field (same as defining).
   */
  polynomial(): RationalPolynomial {
    return this._polynomial;
  }

  /**
   * Return the degree of this number field over Q.
   */
  degree(): number {
    return this._polynomial.degree();
  }

  /**
   * Return the absolute degree of this number field.
   */
  absolute_degree(): number {
    return this.degree();
  }

  /**
   * Return the relative degree of this number field.
   */
  relative_degree(): number {
    return this.degree();
  }

  /**
   * Return the generator of this number field.
   */
  gen(): NumberFieldElement {
    if (!this._cachedGen) {
      // The generator is alpha, represented as [0, 1, 0, ..., 0]
      const coeffs: Rational[] = Array(this.degree())
        .fill(null)
        .map(() => Rational.zero());
      if (coeffs.length > 1) {
        coeffs[1] = Rational.one();
      }
      this._cachedGen = new NumberFieldElement(this, coeffs);
    }
    return this._cachedGen;
  }

  /**
   * Return the number of generators.
   */
  ngens(): number {
    return 1;
  }

  /**
   * Return the zero element.
   */
  zero(): NumberFieldElement {
    return new NumberFieldElement(this, []);
  }

  /**
   * Return the one element.
   */
  one(): NumberFieldElement {
    return new NumberFieldElement(this, [Rational.one()]);
  }

  /**
   * Return the discriminant of the ring of integers of this number field.
   *
   * This is `disc(O_K)`, *not* the discriminant of the defining polynomial:
   * `NumberField(x^3 + x^2 - 2*x + 8).disc()` is `-503`, while the polynomial
   * discriminant is `-2012`.
   *
   * SageMath computes `ZZ(self.pari_polynomial().nfdisc())`; the port runs the
   * Round 2 algorithm behind PARI's `nfdisc`.
   *
   * @see Reference: sage/rings/number_field/number_field.py:5760 (discriminant)
   */
  discriminant(): bigint {
    return this._integralData().nf.disc;
  }

  /**
   * Return the absolute discriminant.
   */
  absolute_discriminant(): bigint {
    return this.discriminant();
  }

  /**
   * Return the signature (r1, r2) of this number field.
   *
   * r1 = number of real embeddings
   * r2 = number of pairs of complex conjugate embeddings
   *
   * Note: This requires finding roots of the defining polynomial, which
   * is expensive in general. For common cases, we use shortcuts.
   */
  signature(): [number, number] {
    if (this._cachedSignature) {
      return this._cachedSignature;
    }

    const n = this.degree();

    // For quadratic fields, check the discriminant
    if (n === 2) {
      // x^2 - bx + c = 0 has discriminant b^2 - 4c
      const b = this._polynomial.getCoeff(1);
      const c = this._polynomial.getCoeff(0);
      const disc = b.mul(b).sub(c.mul(new Rational(4n)));
      if (disc.numerator > 0n) {
        this._cachedSignature = [2, 0];
      } else {
        this._cachedSignature = [0, 1];
      }
      return this._cachedSignature;
    }

    // For higher degrees, we'd need to count real roots
    // Use Sturm's theorem or numerical methods
    // For now, use a heuristic based on the constant term sign
    // This is NOT correct in general!
    const realRoots = this._countRealRoots();
    const r1 = realRoots;
    const r2 = (n - r1) / 2;
    this._cachedSignature = [r1, r2];
    return this._cachedSignature;
  }

  /**
   * Count real roots using Sturm's theorem.
   */
  private _countRealRoots(): number {
    const f = this._polynomial;
    const fp = f.derivative();

    // Build Sturm sequence
    const sturmSeq: RationalPolynomial[] = [f, fp];
    while (!sturmSeq[sturmSeq.length - 1]!.isZero()) {
      const n = sturmSeq.length;
      const remainder = sturmSeq[n - 2]!.mod(sturmSeq[n - 1]!);
      sturmSeq.push(remainder.neg());
    }
    sturmSeq.pop(); // Remove the zero polynomial

    // Count sign changes at -infinity and +infinity
    const atNegInf = sturmSeq.map((p) => {
      const lc = p.leadingCoefficient();
      const deg = p.degree();
      return deg % 2 === 0 ? lc.sign : -lc.sign;
    });

    const atPosInf = sturmSeq.map((p) => p.leadingCoefficient().sign);

    const countChanges = (signs: bigint[]): number => {
      let changes = 0;
      let lastNonZero = 0n;
      for (const s of signs) {
        if (s !== 0n) {
          if (lastNonZero !== 0n && s !== lastNonZero) {
            changes++;
          }
          lastNonZero = s;
        }
      }
      return changes;
    };

    return countChanges(atNegInf) - countChanges(atPosInf);
  }

  /**
   * Check if this field is totally real (all embeddings are real).
   */
  is_totally_real(): boolean {
    const [r1, r2] = this.signature();
    return r2 === 0;
  }

  /**
   * Check if this field is totally imaginary (no real embeddings).
   */
  is_totally_imaginary(): boolean {
    const [r1] = this.signature();
    return r1 === 0;
  }

  /**
   * Return an integral basis for this number field.
   *
   * Returns power basis [1, alpha, alpha^2, ..., alpha^{n-1}].
   * Note: This is the power basis, not necessarily an integral basis for O_K.
   */
  power_basis(): NumberFieldElement[] {
    const n = this.degree();
    const basis: NumberFieldElement[] = [];
    let current = this.one();

    for (let i = 0; i < n; i++) {
      basis.push(current);
      if (i < n - 1) {
        current = current.mul(this.gen());
      }
    }

    return basis;
  }

  /**
   * Coerce x into this number field.
   */
  __call__(x: bigint | Rational | number | NumberFieldElement): NumberFieldElement {
    if (x instanceof NumberFieldElement) {
      if (x.parent() === this) {
        return x;
      }
      throw new ValueError('element is from a different number field');
    }

    if (typeof x === 'number') {
      x = BigInt(Math.floor(x));
    }

    if (typeof x === 'bigint') {
      return NumberFieldElement.fromBigInt(this, x);
    }

    if (x instanceof Rational) {
      return NumberFieldElement.fromRational(this, x);
    }

    throw new ValueError(`cannot coerce ${x} into number field`);
  }

  // The following methods provide number field functionality.
  // NOTE: Full PARI-based implementation requires nfinit/bnfinit which are not yet in parigp-ts.
  // See DEVIATIONS.md for details on the pure TypeScript implementations.

  /**
   * Return the ring of integers of this number field.
   *
   * For the equation order Z[alpha], this returns the order generated by
   * the power basis. For the actual maximal order, PARI's nfinit is needed.
   *
   * @see Reference: sage/rings/number_field/number_field.py:ring_of_integers
   */
  ring_of_integers(): Order {
    // Dynamically import to avoid circular dependency
    const { AbsoluteOrder } = require('./order.js');
    return new AbsoluteOrder(this, this.integral_basis());
  }

  /**
   * Alias for ring_of_integers.
   * @see Reference: sage/rings/number_field/number_field.py:maximal_order
   */
  maximal_order(): Order {
    return this.ring_of_integers();
  }

  /**
   * Return the class group of this number field.
   *
   * NOTE: Full implementation requires PARI bnfinit. Currently returns
   * trivial class group for quadratic fields where class number computation
   * is available, throws NotImplementedError otherwise.
   *
   * @see Reference: sage/rings/number_field/number_field.py:class_group
   */
  class_group(): ClassGroup {
    const { ClassGroup, quadraticClassGroupInvariants } = require('./class_group.js');

    if (this.degree() === 2) {
      const disc = this.discriminant();
      const invariants: bigint[] = quadraticClassGroupInvariants(disc);
      return new ClassGroup(this, invariants, []);
    }

    if (this._classGroupIsProvablyTrivial()) {
      return new ClassGroup(this, [], []);
    }

    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: class_group of a field of degree > 2 needs PARI bnfinit ' +
        '(Buchmann subexponential relation search over a factor base, HNF/SNF of the ' +
        'relation matrix and the regulator); only the case where the Minkowski bound ' +
        'admits no prime ideal, hence a provably trivial class group, is implemented'
    );
  }

  /**
   * Is the class group provably trivial by the Minkowski bound alone?
   *
   * Every ideal class contains an integral ideal of norm at most
   * `M_K = sqrt(|d_K|) * (4/pi)^r2 * n!/n^n`, and every such ideal is a product
   * of prime ideals of norm at most `M_K`.  So if *no* prime ideal has norm
   * `<= M_K`, the only integral ideal of norm `<= M_K` is `O_K` itself and the
   * class group is trivial.  That is a proof, not an estimate: `B` below is a
   * rational **upper** bound for `M_K` computed with `4/pi <= 1.2733` and
   * `sqrt(|d_K|) <= isqrt(|d_K|) + 1`, so the test is conservative.
   *
   * Returns `false` whenever the criterion is inconclusive -- it never claims a
   * class group is trivial without the proof above.
   */
  private _classGroupIsProvablyTrivial(): boolean {
    const n = this.degree();
    const [, r2] = this.signature();
    const d = this.discriminant();
    const absd = d < 0n ? -d : d;
    // n! / n^n  (exact)
    let fact = 1n;
    for (let i = 2n; i <= BigInt(n); i++) fact *= i;
    let B = new Rational(fact, BigInt(n) ** BigInt(n));
    B = B.mul(new Rational(isqrt(absd) + 1n));
    for (let i = 0; i < r2; i++) B = B.mul(new Rational(12733n, 10000n));
    // ceil(B)
    let bound = B.numerator / B.denominator;
    if (bound * B.denominator < B.numerator) bound += 1n;
    if (bound > 1000000n) return false; // hopeless, and certainly not conclusive
    for (let p = 2n; p <= bound; p++) {
      if (!is_prime(p)) continue;
      let dec: Array<[NumberFieldIdeal, bigint]>;
      try {
        dec = this.decomposition(p);
      } catch {
        return false; // cannot decompose p: inconclusive, never claim triviality
      }
      for (const [P] of dec) {
        const nrm = P.norm();
        if (nrm.denominator === 1n && nrm.numerator <= bound) return false;
      }
    }
    return true;
  }

  /**
   * Return the class number of this number field.
   *
   * For quadratic fields, uses analytic class number formula.
   * For other fields, requires PARI.
   *
   * @see Reference: sage/rings/number_field/number_field.py:class_number
   */
  class_number(): bigint {
    if (this.degree() === 2) {
      const disc = this.discriminant();
      return this._computeQuadraticClassNumber(disc);
    }
    if (this._classGroupIsProvablyTrivial()) return 1n;
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: class_number of a field of degree > 2 needs PARI bnfinit ' +
        '(Buchmann subexponential relation search over a factor base, HNF/SNF of the ' +
        'relation matrix and the regulator); only the case where the Minkowski bound ' +
        'admits no prime ideal, hence a provable class number of 1, is implemented'
    );
  }

  /**
   * Compute the class number of a quadratic field from its (fundamental)
   * discriminant.  SageMath calls PARI's `qfbclassno`/`bnfinit`; this counts
   * classes of binary quadratic forms of the same discriminant.
   */
  private _computeQuadraticClassNumber(D: bigint): bigint {
    const { quadraticClassNumber } = require('./class_group.js');
    return quadraticClassNumber(D) as bigint;
  }

  /**
   * Return the unit group of this number field.
   *
   * NOTE: Full implementation requires PARI bnfinit. Currently provides
   * basic structure for quadratic fields.
   *
   * @see Reference: sage/rings/number_field/number_field.py:unit_group
   */
  unit_group(): UnitGroup {
    const { quadraticUnitGroup, UnitGroup } = require('./unit_group.js');

    if (this.degree() === 2) {
      return quadraticUnitGroup(this);
    }

    // For higher degree fields, create unit group with minimal info
    const [r1, r2] = this.signature();
    const rank = r1 + r2 - 1;

    // Torsion is at least {+/-1}
    const torsionOrder = 2n;
    const torsionGen = this.__call__(-1n);

    return new UnitGroup(this, torsionOrder, torsionGen, []);
  }

  /**
   * Return the regulator of this number field.
   *
   * The regulator is the determinant of the logarithmic embedding matrix
   * of the fundamental units.
   *
   * NOTE: Requires numerical computation and fundamental units.
   *
   * @see Reference: sage/rings/number_field/number_field.py:regulator
   */
  regulator(): number {
    // Sage: RealField(53)(self.pari_bnf(proof).bnf_get_reg()).  For a quadratic
    // field PARI's bnfinit gets the regulator from `quadregulator`, which is
    // log(epsilon) for the fundamental unit that `quadunit` produces; the unit
    // group knows how to embed it.  Higher degrees still need bnfinit.
    if (this.degree() === 2) {
      return this.unit_group().regulator();
    }
    throw new NotImplementedError('regulator requires numerical computation with PARI');
  }

  /**
   * Return the Galois group of this number field.
   *
   * @see Reference: sage/rings/number_field/number_field.py:galois_group
   */
  galois_group(): GaloisGroup {
    const { GaloisGroup } = require('./galois_group.js');
    return new GaloisGroup(this);
  }

  /**
   * Return the automorphisms of this number field.
   *
   * For a Galois extension, returns all automorphisms.
   * For a non-Galois extension, returns only the identity.
   *
   * @see Reference: sage/rings/number_field/number_field.py:automorphisms
   */
  automorphisms(): NumberFieldAutomorphism[] {
    const n = this.degree();
    const alpha = this.gen();

    if (n === 1) {
      return [new NumberFieldAutomorphism(this, alpha)];
    }

    if (n === 2) {
      // For x^2 + bx + c = 0, if alpha is one root, -b - alpha is the other.
      const b = this._polynomial.getCoeff(1);
      const otherRoot = alpha.neg().sub(NumberFieldElement.fromRational(this, b));
      const auts = [new NumberFieldAutomorphism(this, alpha)];
      if (!otherRoot.eq(alpha)) {
        auts.push(new NumberFieldAutomorphism(this, otherRoot));
      }
      return auts;
    }

    // General case: PARI's nfgaloisconj on the integral model, translated back
    // to the power basis of alpha via theta = scale * alpha.
    const { scale } = this._integralData();
    const conjugates = nfgaloisconj(this._integralData().nf.g);
    const scaleR = new Rational(scale);
    const auts: NumberFieldAutomorphism[] = [];
    for (const c of conjugates) {
      // beta_theta = sum c_i theta^i = sum c_i scale^i alpha^i, and the image of
      // alpha is beta_theta / scale.
      const coeffs: Rational[] = [];
      let scalePow = Rational.one();
      for (let i = 0; i < n; i++) {
        coeffs.push((c[i] ?? Rational.zero()).mul(scalePow).div(scaleR));
        scalePow = scalePow.mul(scaleR);
      }
      auts.push(new NumberFieldAutomorphism(this, new NumberFieldElement(this, coeffs)));
    }
    return auts;
  }

  /**
   * Return embeddings into a target field.
   *
   * NOTE: Full implementation requires numerical root finding.
   *
   * @see Reference: sage/rings/number_field/number_field.py:embeddings
   */
  embeddings(target?: unknown): unknown[] {
    throw new NotImplementedError('embeddings requires numerical root finding');
  }

  /**
   * Return the real embeddings.
   * @see Reference: sage/rings/number_field/number_field.py:real_embeddings
   */
  real_embeddings(): unknown[] {
    const [r1] = this.signature();
    if (r1 === 0) return [];
    throw new NotImplementedError('real_embeddings requires numerical computation');
  }

  /**
   * Return the complex embeddings.
   * @see Reference: sage/rings/number_field/number_field.py:complex_embeddings
   */
  complex_embeddings(): unknown[] {
    const [, r2] = this.signature();
    if (r2 === 0) return [];
    throw new NotImplementedError('complex_embeddings requires numerical computation');
  }

  /**
   * Return the different ideal.
   *
   * The different is the inverse of the codifferent and measures ramification.
   *
   * @see Reference: sage/rings/number_field/number_field.py:different
   */
  different(): NumberFieldIdeal {
    throw new NotImplementedError('different requires PARI idealdiff');
  }

  /**
   * Return the prime ideals above a rational prime p.
   *
   * This factors pO_K into prime ideals.
   *
   * @see Reference: sage/rings/number_field/number_field.py:primes_above
   */
  primes_above(p: bigint): NumberFieldIdeal[] {
    return this.decomposition(p).map(([P]) => P);
  }

  /**
   * Return *one* prime ideal above `p` (the first one found), matching
   * SageMath's `K.prime_above(p)`.
   *
   * @see Reference: sage/rings/number_field/number_field.py:prime_above
   */
  prime_above(p: bigint, options?: { degree?: bigint }): NumberFieldIdeal {
    const primes = this.primes_above(p);
    if (options?.degree !== undefined) {
      const match = primes.filter((P) => P.residue_class_degree() === options.degree);
      if (match.length === 0) {
        throw new ValueError(`No prime of degree ${options.degree} above ${p}`);
      }
      return match[match.length - 1]!;
    }
    return primes[primes.length - 1]!;
  }

  /**
   * The multiplication table of `O_K` in the basis `_pari_integral_basis()`:
   * `mul[i][j]` are the (integer) coordinates of `w_i * w_j`.  Cached.
   */
  private _orderMulTable(): MulTable {
    if (this._cachedMulTable) return this._cachedMulTable;
    const n = this.degree();
    const basis = this._pari_integral_basis();
    if (!basis[0]!.eq(this.one())) {
      throw new ValueError('the integral basis does not start with 1');
    }
    const W: Rational[][] = basis.map((b) => b.list());
    const Winv = ratInverse(W);
    const coords = (x: NumberFieldElement): bigint[] => {
      const l = x.list();
      const row: bigint[] = [];
      for (let k = 0; k < n; k++) {
        let acc = Rational.zero();
        for (let t = 0; t < n; t++) acc = acc.add(l[t]!.mul(Winv[t]![k]!));
        if (acc.denominator !== 1n) {
          throw new ValueError('the integral basis is not closed under multiplication');
        }
        row.push(acc.numerator);
      }
      return row;
    };
    const table: MulTable = [];
    for (let i = 0; i < n; i++) {
      const row: bigint[][] = [];
      for (let j = 0; j < n; j++) row.push(coords(basis[i]!.mul(basis[j]!)));
      table.push(row);
    }
    this._cachedMulTable = table;
    return table;
  }

  /**
   * Turn one `primedec` entry into a `NumberFieldIdeal`.
   *
   * SageMath/PARI return a two-element representation `(p, alpha)`; we look for
   * one by testing candidate `alpha` in `P` until `N(p, alpha) = p^f`, which
   * proves `(p, alpha) = P` because `(p, alpha) subseteq P` already.  If no
   * such `alpha` turns up, the ideal is returned on the full generating set
   * `p, w'_1, ..., w'_r` -- still exactly `P`, only less pretty.
   *
   * @see Reference: reference/pari/src/basemath/base2.c:2085 (idealprimedec_kummer)
   */
  private _idealFromPrimeDec(p: bigint, entry: PrimeDecEntry): NumberFieldIdeal {
    const { NumberFieldIdeal } = require('./number_field_ideal.js');
    const n = this.degree();
    const basis = this._pari_integral_basis();
    const pElem = this.__call__(p);
    const target = new Rational(p ** entry.f);
    const build = (coeffs: bigint[]): NumberFieldElement => {
      let acc = this.zero();
      for (let i = 0; i < n; i++) {
        if (coeffs[i] === 0n) continue;
        acc = acc.add(basis[i]!.scalarMul(new Rational(coeffs[i]!)));
      }
      return acc;
    };
    const candidates: bigint[][] = entry.gens.map((g) => [...g]);
    // small combinations of the generators, deterministically enumerated
    for (const g of entry.gens) {
      for (const h of entry.gens) {
        if (g === h) continue;
        candidates.push(g.map((x, i) => x + h[i]!));
        candidates.push(g.map((x, i) => x - h[i]!));
      }
    }
    let seed = 1n;
    for (let t = 0; t < 40; t++) {
      const c = new Array<bigint>(n).fill(0n);
      for (const g of entry.gens) {
        seed = (seed * 1103515245n + 12345n) % 2147483648n;
        const m = (seed % (2n * p + 1n)) - p;
        for (let i = 0; i < n; i++) c[i] = c[i]! + m * g[i]!;
      }
      candidates.push(c);
    }
    for (const c of candidates) {
      const alpha = build(c);
      if (alpha.is_zero()) continue;
      const I = new NumberFieldIdeal(this, [pElem, alpha]);
      if (I.norm().eq(target)) return I;
    }
    return new NumberFieldIdeal(this, [pElem, ...entry.gens.map(build)]);
  }

  /**
   * Factor `p * O_K` into prime ideals: an array of `[P, e]` pairs.
   *
   * SageMath delegates to PARI's `idealprimedec`, which has two branches and so
   * do we:
   *
   * - `p` prime to `[O_K : Z[theta]]`: the Dedekind--Kummer theorem, i.e.
   *   factoring `g mod p = prod gbar_i^{e_i}` gives
   *   `p O_K = prod (p, g_i(theta))^{e_i}` with residue degree `deg g_i`;
   * - otherwise: the Buchmann--Lenstra "round 4" algorithm (`primedec`), which
   *   needs no monogenic generator and therefore also handles *inessential
   *   discriminant divisors*, where no `gamma` with `p` prime to
   *   `[O_K : Z[gamma]]` exists at all (the smallest example being `p = 2` in
   *   `Q[x]/(x^3 - x^2 - 2x - 8)`, of discriminant -503).
   *
   * @see Reference: reference/pari/src/basemath/base2.c:2248 (primedec_aux)
   */
  decomposition(p: bigint): Array<[NumberFieldIdeal, bigint]> {
    const { NumberFieldIdeal } = require('./number_field_ideal.js');
    if (this._integralData().nf.index % p === 0n) {
      const dec = primedec(this._orderMulTable(), p);
      const out: Array<[NumberFieldIdeal, bigint]> = [];
      for (const entry of dec) out.push([this._idealFromPrimeDec(p, entry), entry.e]);
      return out;
    }
    const { gamma, minpoly } = this._decompositionGenerator(p);
    const pElem = this.__call__(p);
    const out: Array<[NumberFieldIdeal, bigint]> = [];
    for (const [gi, e] of fpFactor(minpoly, p)) {
      // g_i(gamma), reduced in K
      let acc = this.zero();
      let pw = this.one();
      for (let i = 0; i < gi.length; i++) {
        acc = acc.add(pw.scalarMul(new Rational(gi[i]!)));
        pw = pw.mul(gamma);
      }
      out.push([new NumberFieldIdeal(this, [pElem, acc]), BigInt(e)]);
    }
    return out;
  }

  /**
   * The generator of `O_K` over `Z` used by the Dedekind--Kummer branch of
   * `decomposition`: `theta = scale * alpha`, which works exactly when `p` is
   * prime to `[O_K : Z[theta]]`.  This is PARI's `p_2` branch; the caller has
   * already routed the remaining primes to `primedec` (round 4), so the
   * fallback here is only a guard.
   *
   * @see Reference: reference/pari/src/basemath/base2.c:2248 (primedec_aux)
   */
  private _decompositionGenerator(p: bigint): { gamma: NumberFieldElement; minpoly: ZPoly } {
    const { nf, scale } = this._integralData();
    if (nf.index % p !== 0n) {
      const theta = this.gen().scalarMul(new Rational(scale));
      return { gamma: theta, minpoly: [...nf.g] };
    }
    throw new ValueError(
      `${p} divides the index [O_K : Z[theta]]; use the Buchmann-Lenstra branch`
    );
  }

  /**
   * Factor an ideal.
   * @see Reference: sage/rings/number_field/number_field.py:factor
   */
  factor(ideal: NumberFieldIdeal): Array<[NumberFieldIdeal, bigint]> {
    throw new NotImplementedError('factor requires PARI idealfactor');
  }

  /**
   * Return an ideal of this number field.
   * @see Reference: sage/rings/number_field/number_field.py:ideal
   */
  ideal(...gens: (bigint | Rational | NumberFieldElement)[]): NumberFieldIdeal {
    const { NumberFieldIdeal } = require('./number_field_ideal.js');

    if (gens.length === 0) {
      throw new ValueError('must specify at least one generator');
    }

    const nfGens: NumberFieldElement[] = gens.map((g) => {
      if (g instanceof NumberFieldElement) {
        return g;
      }
      return this.__call__(g);
    });

    return new NumberFieldIdeal(this, nfGens);
  }

  /**
   * Return a fractional ideal.
   * @see Reference: sage/rings/number_field/number_field.py:fractional_ideal
   */
  fractional_ideal(...gens: (bigint | Rational | NumberFieldElement)[]): NumberFieldIdeal {
    return this.ideal(...gens);
  }

  /**
   * Return an integral basis for this number field.
   *
   * For the equation order, this is the power basis.
   * For the maximal order, requires PARI.
   *
   * @see Reference: sage/rings/number_field/number_field.py:integral_basis
   */
  integral_basis(): NumberFieldElement[] {
    if (this._cachedIntegralBasis) return [...this._cachedIntegralBasis];
    const n = this.degree();
    const { nf } = this._integralData();
    // Sage returns `maximal_order().basis()`, which is the *row echelon* basis
    // of the lattice (pivots read left to right), e.g.
    // NumberField(x^3 + x^2 - 2*x + 8).integral_basis() == [1, 1/2*a^2 + 1/2*a, a^2].
    // PARI's nfbasis returns the transposed (increasing-degree) shape, so
    // re-normalise here.
    void nf;
    const pari = this._pari_integral_basis();
    let common = 1n;
    for (const b of pari) {
      for (const c of b.list()) common = intLcm(common, c.denominator);
    }
    const rows = pari.map((b) => b.list().map((c) => c.numerator * (common / c.denominator)));
    const H = hnf(rows, n);
    const basis = H.map(
      (row) =>
        new NumberFieldElement(
          this,
          row.map((x) => new Rational(x, common))
        )
    );
    this._cachedIntegralBasis = basis;
    return [...basis];
  }

  /**
   * The integral basis in PARI's `nfbasis` shape: `w_0 = 1` and the `i`-th
   * element involves only `alpha^0, ..., alpha^i`.  This is the basis used for
   * ideal Hermite normal forms, so that the first HNF row generates `I cap Q`.
   *
   * @see Reference: sage/rings/number_field/number_field.py:_pari_integral_basis
   */
  _pari_integral_basis(): NumberFieldElement[] {
    if (this._cachedPariBasis) return [...this._cachedPariBasis];
    const { nf, scale } = this._integralData();
    const n = this.degree();
    const scaleR = new Rational(scale);
    const basis: NumberFieldElement[] = [];
    for (let i = 0; i < n; i++) {
      // (1/den) * sum_j M[i][j] theta^j  with theta = scale * alpha
      const coeffs: Rational[] = [];
      let scalePow = Rational.one();
      for (let j = 0; j < n; j++) {
        coeffs.push(new Rational(nf.basis[i]![j]!, nf.den).mul(scalePow));
        scalePow = scalePow.mul(scaleR);
      }
      basis.push(new NumberFieldElement(this, coeffs));
    }
    this._cachedPariBasis = basis;
    return [...basis];
  }

  /**
   * Check if this field is a CM field.
   *
   * A CM field is a totally imaginary quadratic extension of a totally real field.
   *
   * @see Reference: sage/rings/number_field/number_field.py:is_CM
   */
  is_CM(): boolean {
    const n = this.degree();

    // Degree must be even for CM
    if (n % 2 !== 0) return false;

    // For quadratic fields
    if (n === 2) {
      // Imaginary quadratic fields are CM (over Q which is totally real)
      return this.is_totally_imaginary();
    }

    // For higher degree, need to check for totally real subfield
    // This requires computing subfields and checking their signatures
    throw new NotImplementedError('is_CM for degree > 2 requires subfield computation');
  }

  /**
   * Check if this field is Galois over Q.
   *
   * A number field K/Q is Galois iff it has [K:Q] automorphisms.
   *
   * @see Reference: sage/rings/number_field/number_field.py:is_galois
   */
  is_galois(): boolean {
    const n = this.degree();

    // Degree 1 is always Galois
    if (n === 1) return true;

    // Degree 2 is always Galois
    if (n === 2) return true;

    // For higher degrees, count automorphisms
    try {
      const auts = this.automorphisms();
      return auts.length === n;
    } catch {
      // Can't compute automorphisms, use polynomial criterion
      // K is Galois iff the defining polynomial splits completely in K
      throw new NotImplementedError('is_galois requires polynomial factorization for degree > 2');
    }
  }

  /**
   * Check if this field is abelian over Q.
   *
   * An extension is abelian if its Galois group is abelian.
   *
   * @see Reference: sage/rings/number_field/number_field.py:is_abelian
   */
  is_abelian(): boolean {
    // Must be Galois first
    if (!this.is_galois()) return false;

    // Degree 1 or 2 fields are always abelian
    if (this.degree() <= 2) return true;

    // Cyclotomic fields are abelian
    if (this instanceof CyclotomicField) return true;

    // Check if Galois group is abelian
    try {
      const G = this.galois_group();
      return G.is_abelian();
    } catch {
      throw new NotImplementedError('is_abelian requires Galois group computation');
    }
  }

  /**
   * Return the Galois closure.
   *
   * The Galois closure is the smallest Galois extension containing this field.
   *
   * @see Reference: sage/rings/number_field/number_field.py:galois_closure
   */
  galois_closure(): NumberField {
    // For Galois fields, return self
    if (this.is_galois()) {
      return this;
    }

    throw new NotImplementedError('galois_closure requires polynomial resultant computation');
  }

  /**
   * Return the splitting field.
   *
   * The splitting field is the field over which the defining polynomial
   * factors completely into linear factors. Same as Galois closure.
   *
   * @see Reference: sage/rings/number_field/number_field.py:splitting_field
   */
  splitting_field(): NumberField {
    return this.galois_closure();
  }

  /**
   * Return the compositum with another number field.
   *
   * @see Reference: sage/rings/number_field/number_field.py:composite_fields
   */
  composite_fields(other: NumberField): NumberField[] {
    throw new NotImplementedError('composite_fields requires resultant computation');
  }

  /**
   * Return subfields of this number field.
   *
   * @see Reference: sage/rings/number_field/number_field.py:subfields
   */
  subfields(degree?: number): Array<[NumberField, NumberFieldElement]> {
    const n = this.degree();

    // For quadratic fields, only subfield is Q
    if (n === 2) {
      return [];
    }

    throw new NotImplementedError('subfields requires polynomial factorization');
  }

  /**
   * Return the Dedekind zeta function.
   *
   * @see Reference: sage/rings/number_field/number_field.py:zeta_function
   */
  zeta_function(): unknown {
    throw new NotImplementedError('zeta_function requires special function computation');
  }

  toString(): string {
    return `Number Field in ${this._name} with defining polynomial ${this._polynomial}`;
  }
}

/**
 * A quadratic number field Q(sqrt(d)).
 */
export class QuadraticField extends NumberField {
  private readonly _D: bigint;

  private constructor(D: bigint, name: string) {
    // Sage uses x^2 - D verbatim: QuadraticField(8).gen()^2 == 8.
    const poly = new RationalPolynomial([new Rational(-D), Rational.zero(), Rational.one()]);
    super(poly, name, undefined, false);
    this._D = D;
  }

  /**
   * Create a quadratic field Q(sqrt(D)) as `NumberField(x^2 - D)`.
   *
   * @see Reference: sage/rings/number_field/number_field.py:970 (QuadraticField)
   */
  static create(D: bigint, name: string = 'a'): QuadraticField {
    const abs = D < 0n ? -D : D;
    const r = isqrtBigInt(abs);
    if (D >= 0n && r * r === D) {
      throw new ValueError('D must not be a perfect square.');
    }
    return new QuadraticField(D, name);
  }

  /**
   * Return D where this is Q(sqrt(D)), i.e. the constant term of `x^2 - D`.
   */
  get D(): bigint {
    return this._D;
  }

  /**
   * The squarefree integer `d` with `Q(sqrt(D)) = Q(sqrt(d))`.
   */
  get d(): bigint {
    return squarefreePart(this._D);
  }

  /**
   * Check if this is a real quadratic field.
   */
  override is_totally_real(): boolean {
    return this._D > 0n;
  }
}

/** Integer square root (floor). */
function isqrtBigInt(n: bigint): bigint {
  if (n < 0n) throw new ValueError('isqrt of a negative number');
  if (n < 2n) return n;
  let x = n;
  let y = (x + 1n) / 2n;
  while (y < x) {
    x = y;
    y = (x + n / x) / 2n;
  }
  return x;
}

/**
 * A cyclotomic field Q(zeta_n).
 */
export class CyclotomicField extends NumberField {
  private readonly _order: bigint;

  constructor(n: bigint, name: string = 'zeta') {
    if (n < 1n) {
      throw new ValueError('n must be at least 1');
    }
    // Compute the n-th cyclotomic polynomial (irreducible by construction, so
    // Sage passes check=False here as well).
    const poly = cyclotomicPolynomial(n);
    super(poly, name, undefined, false);
    this._order = n;
  }

  /**
   * Create a cyclotomic field Q(zeta_n).
   */
  static create(n: bigint, name: string = 'zeta'): CyclotomicField {
    return new CyclotomicField(n, name);
  }

  /**
   * Return n where this is Q(zeta_n).
   */
  get n(): bigint {
    return this._order;
  }

  /**
   * Return the degree phi(n).
   */
  override degree(): number {
    return Number(eulerPhi(this._order));
  }

  /**
   * The automorphisms of `Q(zeta_n)` are `zeta -> zeta^k` for `k` coprime to
   * `n`; the Galois group is `(Z/nZ)^*`.
   *
   * @see Reference: sage/rings/number_field/number_field.py:automorphisms
   */
  override automorphisms(): NumberFieldAutomorphism[] {
    const zeta = this.gen();
    const out: NumberFieldAutomorphism[] = [];
    for (const k of this.galois_exponents()) {
      out.push(new NumberFieldAutomorphism(this, zeta.pow(k)));
    }
    return out;
  }

  /**
   * The residues `k` coprime to `n`, in increasing order; `k = 1` first.
   * These index the automorphisms `zeta -> zeta^k`.
   */
  galois_exponents(): bigint[] {
    const n = this._order;
    const out: bigint[] = [];
    for (let k = 1n; k <= n; k++) {
      if (intGcd(k, n) === 1n) out.push(k);
    }
    if (n === 1n) return [1n];
    return out;
  }

  /**
   * Return the primitive root of unity.
   */
  zeta(k?: bigint): NumberFieldElement {
    if (k === undefined || k === 1n) {
      return this.gen();
    }
    return this.gen().pow(k);
  }
}

// Helper functions

/**
 * Compute the squarefree part of n.
 */
function squarefreePart(n: bigint): bigint {
  if (n === 0n) return 0n;

  const sign = n < 0n ? -1n : 1n;
  n = n < 0n ? -n : n;

  let result = 1n;
  let d = 2n;

  while (d * d <= n) {
    let exp = 0n;
    while (n % d === 0n) {
      n = n / d;
      exp++;
    }
    if (exp % 2n === 1n) {
      result *= d;
    }
    d++;
  }

  if (n > 1n) {
    result *= n;
  }

  return sign * result;
}

/**
 * Compute Euler's totient function phi(n).
 */
function eulerPhi(n: bigint): bigint {
  if (n <= 0n) {
    throw new ValueError('n must be positive');
  }
  if (n === 1n) return 1n;

  let result = n;
  let temp = n;
  let p = 2n;

  while (p * p <= temp) {
    if (temp % p === 0n) {
      while (temp % p === 0n) {
        temp = temp / p;
      }
      result = result - result / p;
    }
    p++;
  }

  if (temp > 1n) {
    result = result - result / temp;
  }

  return result;
}

/**
 * Compute the n-th cyclotomic polynomial.
 *
 * Phi_n(x) = product over primitive n-th roots of unity (x - zeta)
 */
function cyclotomicPolynomial(n: bigint): RationalPolynomial {
  if (n === 1n) {
    // Phi_1(x) = x - 1
    return new RationalPolynomial([new Rational(-1n), Rational.one()]);
  }

  if (n === 2n) {
    // Phi_2(x) = x + 1
    return new RationalPolynomial([Rational.one(), Rational.one()]);
  }

  // Use the formula: Phi_n(x) = (x^n - 1) / product_{d|n, d<n} Phi_d(x)
  // Start with x^n - 1
  const coeffs: Rational[] = Array(Number(n) + 1)
    .fill(null)
    .map(() => Rational.zero());
  coeffs[0] = new Rational(-1n);
  coeffs[Number(n)] = Rational.one();
  let result = new RationalPolynomial(coeffs);

  // Divide by Phi_d for each proper divisor d of n
  const divs = divisors(n);
  for (const d of divs) {
    if (d < n) {
      const phiD = cyclotomicPolynomial(d);
      const [q, r] = result.divmod(phiD);
      if (!r.isZero()) {
        throw new Error(`cyclotomicPolynomial: internal error dividing by Phi_${d}`);
      }
      result = q;
    }
  }

  return result;
}

/**
 * An automorphism of a number field.
 *
 * An automorphism is determined by where it sends the generator.
 *
 * @see Reference: sage/rings/number_field/morphism.py
 */
export class NumberFieldAutomorphism {
  private readonly _domain: NumberField;
  private readonly _image_of_gen: NumberFieldElement;

  constructor(domain: NumberField, image_of_gen: NumberFieldElement) {
    this._domain = domain;
    this._image_of_gen = image_of_gen;
  }

  /**
   * Return the domain of this automorphism.
   */
  domain(): NumberField {
    return this._domain;
  }

  /**
   * Return the codomain (same as domain for automorphisms).
   */
  codomain(): NumberField {
    return this._domain;
  }

  /**
   * Return the image of the generator.
   */
  im_gens(): NumberFieldElement[] {
    return [this._image_of_gen];
  }

  /**
   * Apply this automorphism to an element.
   *
   * For an element a_0 + a_1*alpha + ... + a_{n-1}*alpha^{n-1},
   * returns a_0 + a_1*beta + ... + a_{n-1}*beta^{n-1}
   * where beta is the image of alpha.
   */
  __call__(x: NumberFieldElement): NumberFieldElement {
    const coeffs = x.list();
    const beta = this._image_of_gen;

    let result = this._domain.zero();
    let betaPower = this._domain.one();

    for (const coeff of coeffs) {
      if (!coeff.isZero()) {
        result = result.add(betaPower.scalarMul(coeff));
      }
      betaPower = betaPower.mul(beta);
    }

    return result;
  }

  /**
   * Check if this is the identity automorphism.
   */
  is_identity(): boolean {
    return this._image_of_gen.eq(this._domain.gen());
  }

  /**
   * Check equality.
   */
  eq(other: NumberFieldAutomorphism): boolean {
    return this._domain === other._domain && this._image_of_gen.eq(other._image_of_gen);
  }

  toString(): string {
    if (this.is_identity()) {
      return `Ring endomorphism of ${this._domain}\n  Defn: ${this._domain._name} |--> ${this._domain._name}`;
    }
    return `Ring endomorphism of ${this._domain}\n  Defn: ${this._domain._name} |--> ${this._image_of_gen}`;
  }
}

/**
 * Create a number field from a polynomial.
 */
export function NumberFieldConstructor(
  polynomial: RationalPolynomial | bigint[],
  name: string,
  embedding?: unknown
): NumberField {
  if (Array.isArray(polynomial)) {
    polynomial = RationalPolynomial.fromBigInts(polynomial);
  }
  return new NumberField(polynomial, name, embedding);
}
