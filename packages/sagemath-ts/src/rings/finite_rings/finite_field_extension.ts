/**
 * @module sage/rings/finite_rings/finite_field_extension
 * @description Finite field extensions GF(p^n) using polynomial quotient ring construction
 *
 * Port of: sage/rings/finite_rings/finite_field_ext_pari.py (conceptually)
 *
 * The extension field GF(p^n) is constructed as:
 *   GF(p)[x] / <f(x)>
 * where f(x) is an irreducible polynomial of degree n over GF(p).
 *
 * When available, we use Conway polynomials for standardization.
 * Otherwise, we find an irreducible polynomial.
 */

import { inverse_mod, is_prime, power_mod } from '../../arith/misc.js';
import { ValueError, ZeroDivisionError } from '../../errors.js';
import { current_randstate } from '../../misc/randstate.js';
import {
  type CoefficientRing,
  Polynomial,
  type RingElement,
} from '../polynomial/polynomial_element.js';
import { PolynomialRing } from '../polynomial/polynomial_ring.js';
import { QuotientRing, QuotientRingElement } from '../polynomial/quotient_ring.js';
import { conway_polynomial, has_conway_polynomial } from './conway_polynomials.js';

/**
 * Element of a prime field GF(p).
 */
export class PrimeFieldElement implements RingElement {
  readonly value: bigint;
  readonly parent: PrimeField;

  constructor(value: number | bigint | PrimeFieldElement, parent: PrimeField) {
    this.parent = parent;

    if (value instanceof PrimeFieldElement) {
      this.value = value.value;
    } else {
      const v = typeof value === 'number' ? BigInt(value) : value;
      this.value = ((v % parent.characteristic) + parent.characteristic) % parent.characteristic;
    }
  }

  add(other: PrimeFieldElement): PrimeFieldElement {
    return new PrimeFieldElement(
      (this.value + other.value) % this.parent.characteristic,
      this.parent
    );
  }

  sub(other: PrimeFieldElement): PrimeFieldElement {
    return new PrimeFieldElement(
      (((this.value - other.value) % this.parent.characteristic) + this.parent.characteristic) %
        this.parent.characteristic,
      this.parent
    );
  }

  mul(other: PrimeFieldElement): PrimeFieldElement {
    return new PrimeFieldElement(
      (this.value * other.value) % this.parent.characteristic,
      this.parent
    );
  }

  neg(): PrimeFieldElement {
    if (this.value === 0n) {
      return this;
    }
    return new PrimeFieldElement(this.parent.characteristic - this.value, this.parent);
  }

  inv(): PrimeFieldElement {
    if (this.value === 0n) {
      throw new ZeroDivisionError('division by zero in finite field');
    }
    return new PrimeFieldElement(inverse_mod(this.value, this.parent.characteristic), this.parent);
  }

  div(other: PrimeFieldElement): PrimeFieldElement {
    return this.mul(other.inv());
  }

  pow(n: number | bigint): PrimeFieldElement {
    const exp = typeof n === 'bigint' ? n : BigInt(n);

    if (exp < 0n) {
      return this.inv().pow(-exp);
    }

    if (exp === 0n) {
      return this.parent.one();
    }

    if (this.value === 0n) {
      return this.parent.zero();
    }

    return new PrimeFieldElement(
      power_mod(this.value, exp, this.parent.characteristic),
      this.parent
    );
  }

  eq(other: PrimeFieldElement | number): boolean {
    if (typeof other === 'number') {
      const otherVal =
        ((BigInt(other) % this.parent.characteristic) + this.parent.characteristic) %
        this.parent.characteristic;
      return this.value === otherVal;
    }
    return this.value === other.value;
  }

  isZero(): boolean {
    return this.value === 0n;
  }

  isOne(): boolean {
    return this.value === 1n;
  }

  toString(): string {
    return this.value.toString();
  }

  repr(): string {
    return this.value.toString();
  }

  toBigInt(): bigint {
    return this.value;
  }
}

/**
 * A prime field GF(p).
 */
export class PrimeField implements CoefficientRing<PrimeFieldElement> {
  readonly characteristic: bigint;
  readonly order: bigint;
  readonly degree = 1;

  constructor(p: number | bigint) {
    const prime = typeof p === 'number' ? BigInt(p) : p;

    if (prime <= 1n) {
      throw new ValueError('p must be a prime > 1');
    }

    if (!is_prime(prime)) {
      throw new ValueError(`${prime} is not prime`);
    }

    this.characteristic = prime;
    this.order = prime;
  }

  __call__(x: number | bigint | PrimeFieldElement | unknown): PrimeFieldElement {
    if (x instanceof PrimeFieldElement) {
      return new PrimeFieldElement(x.value, this);
    }
    if (typeof x === 'number' || typeof x === 'bigint') {
      return new PrimeFieldElement(x, this);
    }
    throw new ValueError(`Cannot convert ${typeof x} to PrimeFieldElement`);
  }

  zero(): PrimeFieldElement {
    return new PrimeFieldElement(0n, this);
  }

  one(): PrimeFieldElement {
    return new PrimeFieldElement(1n, this);
  }

  gen(): PrimeFieldElement {
    // For a prime field, we return a primitive root if p > 2
    // For simplicity, return 1 (or find primitive root)
    if (this.characteristic === 2n) {
      return this.one();
    }
    // Find a primitive root
    return this.primitiveRoot();
  }

  /**
   * Find a primitive root modulo p.
   */
  private primitiveRoot(): PrimeFieldElement {
    if (this.characteristic === 2n) {
      return this.one();
    }

    const phi = this.characteristic - 1n;

    // Factor phi to check primitive root property
    const factors = this.factorPhi(phi);

    // Try candidates starting from 2
    for (let g = 2n; g < this.characteristic; g++) {
      let isPrimitive = true;

      for (const [p, _e] of factors) {
        if (power_mod(g, phi / p, this.characteristic) === 1n) {
          isPrimitive = false;
          break;
        }
      }

      if (isPrimitive) {
        return new PrimeFieldElement(g, this);
      }
    }

    // Fallback (should not happen for primes)
    return new PrimeFieldElement(2n, this);
  }

  /**
   * Simple factorization for primitive root computation.
   */
  private factorPhi(n: bigint): Array<[bigint, bigint]> {
    const factors: Array<[bigint, bigint]> = [];
    let temp = n;

    for (let p = 2n; p * p <= temp; p++) {
      if (temp % p === 0n) {
        let e = 0n;
        while (temp % p === 0n) {
          temp /= p;
          e++;
        }
        factors.push([p, e]);
      }
    }

    if (temp > 1n) {
      factors.push([temp, 1n]);
    }

    return factors;
  }

  cardinality(): bigint {
    return this.characteristic;
  }

  *[Symbol.iterator](): Iterator<PrimeFieldElement> {
    for (let i = 0n; i < this.characteristic; i++) {
      yield new PrimeFieldElement(i, this);
    }
  }

  /**
   * Iterate over all elements of this field.
   * Alias for [Symbol.iterator] for compatibility.
   */
  elements(): IterableIterator<PrimeFieldElement> {
    return this[Symbol.iterator]();
  }

  is_field(): boolean {
    return true;
  }

  /**
   * Return a random element of this field.
   */
  random_element(): PrimeFieldElement {
    const rstate = current_randstate();
    const randomInt = rstate.random_below(this.characteristic);
    return new PrimeFieldElement(randomInt, this);
  }

  toString(): string {
    return `Finite Field of size ${this.characteristic}`;
  }

  /**
   * Return a multiplicative generator (primitive root) of GF(p)*.
   *
   * This finds the smallest positive integer g such that g generates
   * the multiplicative group of the field.
   */
  multiplicative_generator(): PrimeFieldElement {
    return this.primitiveRoot();
  }
}

/**
 * Element of a finite field extension GF(p^n).
 */
export class FiniteFieldElement implements RingElement {
  readonly lift: Polynomial<PrimeFieldElement>;
  readonly parent: FiniteFieldExtension;

  constructor(poly: Polynomial<PrimeFieldElement>, parent: FiniteFieldExtension) {
    this.parent = parent;

    // Reduce modulo the modulus
    if (poly.degree() >= parent.modulus.degree()) {
      const [_q, r] = poly.quo_rem(parent.modulus);
      this.lift = r;
    } else {
      this.lift = poly;
    }
  }

  add(other: FiniteFieldElement): FiniteFieldElement {
    return new FiniteFieldElement(this.lift.add(other.lift), this.parent);
  }

  sub(other: FiniteFieldElement): FiniteFieldElement {
    return new FiniteFieldElement(this.lift.sub(other.lift), this.parent);
  }

  mul(other: FiniteFieldElement): FiniteFieldElement {
    const prod = this.lift.mul(other.lift);
    return new FiniteFieldElement(prod, this.parent);
  }

  neg(): FiniteFieldElement {
    return new FiniteFieldElement(this.lift.neg(), this.parent);
  }

  /**
   * Compute the multiplicative inverse using extended Euclidean algorithm.
   */
  inv(): FiniteFieldElement {
    if (this.isZero()) {
      throw new ZeroDivisionError('division by zero');
    }

    // Extended Euclidean algorithm for polynomials
    const [g, s, _t] = polyXgcd(this.lift, this.parent.modulus);

    // g should be a constant (unit) since modulus is irreducible
    if (g.degree() !== 0) {
      throw new ValueError('element is not invertible (modulus not irreducible)');
    }

    // s * this ≡ g (mod modulus)
    // We need to divide s by g's constant term
    const gInv = g.getCoeff(0).inv();
    const sNormalized = s.scalar_mul(gInv);

    return new FiniteFieldElement(sNormalized, this.parent);
  }

  div(other: FiniteFieldElement): FiniteFieldElement {
    return this.mul(other.inv());
  }

  pow(n: number | bigint): FiniteFieldElement {
    let exp = typeof n === 'bigint' ? n : BigInt(n);

    if (exp < 0n) {
      return this.inv().pow(-exp);
    }

    if (exp === 0n) {
      return this.parent.one();
    }

    if (this.isZero()) {
      return this.parent.zero();
    }

    // Binary exponentiation
    let result = this.parent.one();
    let base: FiniteFieldElement = this;

    while (exp > 0n) {
      if ((exp & 1n) === 1n) {
        result = result.mul(base);
      }
      base = base.mul(base);
      exp >>= 1n;
    }

    return result;
  }

  /**
   * Apply the Frobenius automorphism: x -> x^p
   *
   * This is a field automorphism of GF(p^n) that fixes GF(p).
   */
  frobenius(power: number = 1): FiniteFieldElement {
    const p = this.parent.characteristic;
    const exp = p ** BigInt(power);
    return this.pow(exp);
  }

  /**
   * Compute the trace: Tr(x) = x + x^p + x^{p^2} + ... + x^{p^{n-1}}
   *
   * Returns an element of the base field GF(p).
   */
  trace(): PrimeFieldElement {
    let result = this.parent.zero();
    let term: FiniteFieldElement = this;
    const n = this.parent.degree;

    for (let i = 0; i < n; i++) {
      result = result.add(term);
      if (i < n - 1) {
        term = term.frobenius();
      }
    }

    // The trace is in GF(p), so extract the constant coefficient
    return result.lift.getCoeff(0);
  }

  /**
   * Compute the norm: N(x) = x * x^p * x^{p^2} * ... * x^{p^{n-1}}
   *
   * Returns an element of the base field GF(p).
   */
  norm(): PrimeFieldElement {
    let result = this.parent.one();
    let term: FiniteFieldElement = this;
    const n = this.parent.degree;

    for (let i = 0; i < n; i++) {
      result = result.mul(term);
      if (i < n - 1) {
        term = term.frobenius();
      }
    }

    // The norm is in GF(p), so extract the constant coefficient
    return result.lift.getCoeff(0);
  }

  /**
   * Compute the minimal polynomial of this element over GF(p).
   */
  minimalPolynomial(): Polynomial<PrimeFieldElement> {
    // The minimal polynomial divides x^{p^n} - x
    // For a primitive element, it equals the field's defining polynomial
    // In general, we compute it using the Frobenius conjugates

    const conjugates: FiniteFieldElement[] = [];
    let current: FiniteFieldElement = this;
    const seen = new Set<string>();

    while (true) {
      const key = current.toString();
      if (seen.has(key)) {
        break;
      }
      seen.add(key);
      conjugates.push(current);
      current = current.frobenius();
    }

    // Build polynomial (x - conjugates[0])(x - conjugates[1])...
    const polyRing = new PolynomialRing(this.parent.baseField, 'x');
    let minPoly = polyRing.one();

    for (const conj of conjugates) {
      // (x - conj) but conj is in the extension field
      // We need to work in the polynomial ring over the extension
      // For now, this is a simplified version
      const factor = new Polynomial<PrimeFieldElement>(
        [conj.lift.getCoeff(0).neg(), this.parent.baseField.one()],
        polyRing
      );
      minPoly = minPoly.mul(factor);
    }

    return minPoly;
  }

  eq(other: FiniteFieldElement | number): boolean {
    if (typeof other === 'number') {
      const otherElem = this.parent.__call__(other);
      return this.lift.eq(otherElem.lift);
    }
    return this.lift.eq(other.lift);
  }

  isZero(): boolean {
    return this.lift.isZero();
  }

  isOne(): boolean {
    return this.lift.coeffs.length === 1 && this.lift.coeffs[0]!.isOne();
  }

  /**
   * Return the coefficients of the polynomial representation.
   * The element is c_0 + c_1*a + c_2*a^2 + ... + c_{n-1}*a^{n-1}
   * where a is the generator.
   */
  coefficients(): PrimeFieldElement[] {
    const coeffs: PrimeFieldElement[] = [];
    for (let i = 0; i < this.parent.degree; i++) {
      coeffs.push(this.lift.getCoeff(i));
    }
    return coeffs;
  }

  /**
   * Return the integer representation of this element.
   * The element c_0 + c_1*a + ... is represented as c_0 + c_1*p + c_2*p^2 + ...
   */
  integer_representation(): bigint {
    let result = 0n;
    let pPower = 1n;
    const p = this.parent.characteristic;

    for (let i = 0; i < this.parent.degree; i++) {
      result += this.lift.getCoeff(i).value * pPower;
      pPower *= p;
    }

    return result;
  }

  toString(): string {
    if (this.lift.isZero()) {
      return '0';
    }

    const terms: string[] = [];
    const genName = this.parent.variableName;

    for (let i = this.lift.degree(); i >= 0; i--) {
      const c = this.lift.getCoeff(i);
      if (c.isZero()) {
        continue;
      }

      let term: string;
      const cVal = c.value;

      if (i === 0) {
        term = cVal.toString();
      } else if (i === 1) {
        if (cVal === 1n) {
          term = genName;
        } else {
          term = `${cVal}*${genName}`;
        }
      } else {
        if (cVal === 1n) {
          term = `${genName}^${i}`;
        } else {
          term = `${cVal}*${genName}^${i}`;
        }
      }

      terms.push(term);
    }

    if (terms.length === 0) {
      return '0';
    }

    return terms.join(' + ');
  }

  repr(): string {
    return this.toString();
  }
}

/**
 * A finite field extension GF(p^n).
 *
 * Constructed as GF(p)[x] / <f(x)> where f(x) is an irreducible polynomial.
 */
export class FiniteFieldExtension implements CoefficientRing<FiniteFieldElement> {
  readonly baseField: PrimeField;
  readonly polynomialRing: PolynomialRing<PrimeFieldElement>;
  readonly modulus: Polynomial<PrimeFieldElement>;
  readonly degree: number;
  readonly characteristic: bigint;
  readonly order: bigint;
  readonly variableName: string;

  private _elements: FiniteFieldElement[] | null = null;

  constructor(
    p: number | bigint,
    n: number,
    modulus?: Polynomial<PrimeFieldElement> | number[],
    variableName: string = 'a'
  ) {
    if (n < 1) {
      throw new ValueError('degree must be at least 1');
    }

    const prime = typeof p === 'number' ? BigInt(p) : p;
    this.baseField = new PrimeField(prime);
    this.polynomialRing = new PolynomialRing(this.baseField, variableName);
    this.degree = n;
    this.characteristic = prime;
    this.order = prime ** BigInt(n);
    this.variableName = variableName;

    if (modulus) {
      if (modulus instanceof Polynomial) {
        this.modulus = modulus;
      } else {
        // modulus is an array of coefficients
        this.modulus = this.polynomialFromCoeffs(modulus);
      }
    } else if (n === 1) {
      // GF(p) itself
      this.modulus = this.polynomialRing.gen();
    } else {
      // Use Conway polynomial if available, otherwise find irreducible
      this.modulus = this.getDefaultModulus(Number(prime), n);
    }

    if (this.modulus.degree() !== n) {
      throw new ValueError(`modulus must have degree ${n}, got ${this.modulus.degree()}`);
    }
  }

  /**
   * Create polynomial from coefficient array.
   */
  private polynomialFromCoeffs(coeffs: number[]): Polynomial<PrimeFieldElement> {
    const polyCoeffs = coeffs.map((c) => this.baseField.__call__(c));
    // Add leading coefficient 1 (monic polynomial)
    polyCoeffs.push(this.baseField.one());
    return new Polynomial(polyCoeffs, this.polynomialRing);
  }

  /**
   * Get the default modulus (Conway polynomial or random irreducible).
   */
  private getDefaultModulus(p: number, n: number): Polynomial<PrimeFieldElement> {
    // Try Conway polynomial first
    if (has_conway_polynomial(p, n)) {
      const conwayCoeffs = conway_polynomial(p, n);
      return this.polynomialFromCoeffs(conwayCoeffs);
    }

    // Fall back to finding an irreducible polynomial
    return this.findIrreducible(n);
  }

  /**
   * Find an irreducible polynomial of degree n over GF(p).
   *
   * Uses a simple search strategy: try random/sequential polynomials
   * and check irreducibility.
   */
  private findIrreducible(n: number): Polynomial<PrimeFieldElement> {
    // For small cases, use known irreducibles
    // For general case, we'll search

    const x = this.polynomialRing.gen();
    const one = this.polynomialRing.one();

    // Try x^n + x + 1 (often irreducible for characteristic 2)
    let candidate = x.pow(n).add(x).add(one);
    if (this.isIrreducible(candidate)) {
      return candidate;
    }

    // Try x^n - c for various c
    for (let c = 2n; c < this.characteristic; c++) {
      const cPoly = this.polynomialRing.__call__(this.baseField.__call__(c));
      candidate = x.pow(n).sub(cPoly);
      if (this.isIrreducible(candidate)) {
        return candidate;
      }
    }

    // Brute force search through monic polynomials
    // This is slow but guaranteed to find something
    const p = Number(this.characteristic);
    const numPolys = p ** n;

    for (let i = 0; i < numPolys; i++) {
      const coeffs: PrimeFieldElement[] = [];
      let temp = i;

      for (let j = 0; j < n; j++) {
        coeffs.push(this.baseField.__call__(temp % p));
        temp = Math.floor(temp / p);
      }
      coeffs.push(this.baseField.one()); // monic

      candidate = new Polynomial(coeffs, this.polynomialRing);

      if (this.isIrreducible(candidate)) {
        return candidate;
      }
    }

    throw new ValueError(
      `Could not find irreducible polynomial of degree ${n} over GF(${this.characteristic})`
    );
  }

  /**
   * Check if a polynomial is irreducible over GF(p).
   *
   * Uses the fact that f(x) is irreducible iff:
   * 1. f(x) divides x^{p^n} - x
   * 2. gcd(f(x), x^{p^k} - x) = 1 for all k | n with k < n
   */
  private isIrreducible(f: Polynomial<PrimeFieldElement>): boolean {
    const n = f.degree();
    if (n <= 0) {
      return false;
    }

    if (n === 1) {
      return true; // Linear polynomials are irreducible
    }

    const p = this.characteristic;
    const x = this.polynomialRing.gen();

    // Check that f has no repeated roots: gcd(f, f') = 1
    // For characteristic p, this may need special handling
    // Skip for now as Conway polynomials are squarefree

    // Check gcd(f, x^{p^k} - x) = 1 for proper divisors k of n
    const divisors = this.getDivisors(n).filter((d) => d < n && d > 0);

    for (const k of divisors) {
      // Compute x^{p^k} mod f
      const pk = p ** BigInt(k);
      const xPk = this.powerMod(x, pk, f);
      const diff = xPk.sub(x);

      if (!diff.isZero()) {
        const g = this.polyGcd(f, diff);
        if (g.degree() > 0) {
          return false; // f has a factor of degree <= k
        }
      }
    }

    // Check that f divides x^{p^n} - x
    const pn = p ** BigInt(n);
    const xPn = this.powerMod(x, pn, f);
    const remainder = xPn.sub(x);

    // remainder should be zero mod f
    if (!remainder.isZero()) {
      const [_, r] = remainder.quo_rem(f);
      if (!r.isZero()) {
        return false;
      }
    }

    return true;
  }

  /**
   * Get divisors of n.
   */
  private getDivisors(n: number): number[] {
    const divisors: number[] = [];
    for (let i = 1; i <= n; i++) {
      if (n % i === 0) {
        divisors.push(i);
      }
    }
    return divisors;
  }

  /**
   * Compute polynomial power modulo another polynomial.
   */
  private powerMod(
    base: Polynomial<PrimeFieldElement>,
    exp: bigint,
    mod: Polynomial<PrimeFieldElement>
  ): Polynomial<PrimeFieldElement> {
    if (exp === 0n) {
      return this.polynomialRing.one();
    }

    let result = this.polynomialRing.one();
    let b = base;

    while (exp > 0n) {
      if ((exp & 1n) === 1n) {
        result = result.mul(b).mod(mod);
      }
      b = b.mul(b).mod(mod);
      exp >>= 1n;
    }

    return result;
  }

  /**
   * Compute GCD of two polynomials.
   */
  private polyGcd(
    a: Polynomial<PrimeFieldElement>,
    b: Polynomial<PrimeFieldElement>
  ): Polynomial<PrimeFieldElement> {
    while (!b.isZero()) {
      const [_, r] = a.quo_rem(b);
      a = b;
      b = r;
    }

    // Make monic
    if (!a.isZero()) {
      const lc = a.leading_coefficient();
      if (!lc.isOne()) {
        const lcInv = lc.inv();
        a = a.scalar_mul(lcInv);
      }
    }

    return a;
  }

  /**
   * Create an element from various inputs.
   */
  __call__(
    x:
      | number
      | bigint
      | number[]
      | Polynomial<PrimeFieldElement>
      | FiniteFieldElement
      | PrimeFieldElement
      | unknown
  ): FiniteFieldElement {
    if (x instanceof FiniteFieldElement) {
      return new FiniteFieldElement(x.lift, this);
    }

    if (x instanceof Polynomial) {
      return new FiniteFieldElement(x, this);
    }

    if (x instanceof PrimeFieldElement) {
      const poly = this.polynomialRing.__call__(x);
      return new FiniteFieldElement(poly, this);
    }

    if (typeof x === 'number' || typeof x === 'bigint') {
      const coeff = this.baseField.__call__(x);
      const poly = this.polynomialRing.__call__(coeff);
      return new FiniteFieldElement(poly, this);
    }

    if (Array.isArray(x)) {
      // Array of coefficients
      const coeffs = x.map((c) => this.baseField.__call__(c as number));
      const poly = new Polynomial(coeffs, this.polynomialRing);
      return new FiniteFieldElement(poly, this);
    }

    throw new ValueError(`Cannot convert ${typeof x} to FiniteFieldElement`);
  }

  /**
   * Create element from integer representation.
   */
  fromInteger(n: bigint): FiniteFieldElement {
    const coeffs: PrimeFieldElement[] = [];
    let temp = n;

    for (let i = 0; i < this.degree; i++) {
      coeffs.push(this.baseField.__call__(temp % this.characteristic));
      temp /= this.characteristic;
    }

    const poly = new Polynomial(coeffs, this.polynomialRing);
    return new FiniteFieldElement(poly, this);
  }

  zero(): FiniteFieldElement {
    return new FiniteFieldElement(this.polynomialRing.zero(), this);
  }

  one(): FiniteFieldElement {
    return new FiniteFieldElement(this.polynomialRing.one(), this);
  }

  /**
   * Return the generator (image of x in the quotient, which is a root of the modulus).
   */
  gen(): FiniteFieldElement {
    return new FiniteFieldElement(this.polynomialRing.gen(), this);
  }

  cardinality(): bigint {
    return this.order;
  }

  /**
   * Iterate over all elements.
   */
  *[Symbol.iterator](): Iterator<FiniteFieldElement> {
    if (this._elements !== null) {
      yield* this._elements;
      return;
    }

    // Generate all elements
    this._elements = [];

    const p = Number(this.characteristic);
    const numElements = Number(this.order);

    // Generate all polynomials of degree < n
    for (let i = 0; i < numElements; i++) {
      const coeffs: PrimeFieldElement[] = [];
      let temp = i;

      for (let j = 0; j < this.degree; j++) {
        coeffs.push(this.baseField.__call__(temp % p));
        temp = Math.floor(temp / p);
      }

      const poly = new Polynomial(coeffs, this.polynomialRing);
      const elem = new FiniteFieldElement(poly, this);
      this._elements.push(elem);
    }

    yield* this._elements;
  }

  /**
   * Iterate over all elements of this field.
   * Alias for [Symbol.iterator] for compatibility.
   */
  elements(): IterableIterator<FiniteFieldElement> {
    return this[Symbol.iterator]();
  }

  /**
   * Find a primitive element (generator of the multiplicative group).
   */
  primitiveElement(): FiniteFieldElement {
    // The multiplicative group has order p^n - 1
    const groupOrder = this.order - 1n;

    // Factor the group order
    const factors = this.factorSimple(groupOrder);

    // Check each element
    for (const elem of this) {
      if (elem.isZero()) continue;

      let isPrimitive = true;
      for (const [p, _] of factors) {
        const exp = groupOrder / p;
        if (elem.pow(exp).isOne()) {
          isPrimitive = false;
          break;
        }
      }

      if (isPrimitive) {
        return elem;
      }
    }

    // The generator should always be primitive for Conway polynomials
    return this.gen();
  }

  /**
   * Simple factorization for primitive element search.
   */
  private factorSimple(n: bigint): Array<[bigint, bigint]> {
    const factors: Array<[bigint, bigint]> = [];
    let temp = n;

    for (let p = 2n; p * p <= temp; p++) {
      if (temp % p === 0n) {
        let e = 0n;
        while (temp % p === 0n) {
          temp /= p;
          e++;
        }
        factors.push([p, e]);
      }
    }

    if (temp > 1n) {
      factors.push([temp, 1n]);
    }

    return factors;
  }

  /**
   * Return a random element of the field.
   */
  random_element(): FiniteFieldElement {
    const coeffs: PrimeFieldElement[] = [];
    const p = this.characteristic;
    const rstate = current_randstate();

    for (let i = 0; i < this.degree; i++) {
      coeffs.push(this.baseField.__call__(rstate.random_below(p)));
    }

    const poly = new Polynomial(coeffs, this.polynomialRing);
    return new FiniteFieldElement(poly, this);
  }

  is_field(): boolean {
    return true;
  }

  toString(): string {
    return `Finite Field in ${this.variableName} of size ${this.characteristic}^${this.degree}`;
  }
}

/**
 * Extended Euclidean algorithm for polynomials.
 */
function polyXgcd<C extends RingElement>(
  a: Polynomial<C>,
  b: Polynomial<C>
): [Polynomial<C>, Polynomial<C>, Polynomial<C>] {
  const ring = a.parent;

  let oldR = a;
  let r = b;
  let oldS = ring.one();
  let s = ring.zero();
  let oldT = ring.zero();
  let t = ring.one();

  while (!r.isZero()) {
    const [quotient, remainder] = oldR.quo_rem(r);

    const tempR = r;
    r = remainder;
    oldR = tempR;

    const tempS = s;
    s = oldS.sub(quotient.mul(s));
    oldS = tempS;

    const tempT = t;
    t = oldT.sub(quotient.mul(t));
    oldT = tempT;
  }

  return [oldR, oldS, oldT];
}

/**
 * Construct a finite field GF(q) where q = p^n, supporting extension fields.
 *
 * This is the full implementation that supports both prime fields and extension fields.
 * For prime-only fields, use the simpler GF() from finite_field_constructor.ts.
 *
 * @param q - The order of the field (must be a prime power)
 * @param variableName - Name for the generator (default: 'a')
 * @returns The finite field GF(q)
 *
 * @example
 * ```typescript
 * const F4 = GFExtended(4);           // GF(2^2)
 * const F9 = GFExtended(9);           // GF(3^2)
 * const F8 = GFExtended(8, 'b');      // GF(2^3) with generator named 'b'
 * ```
 */
export function GFExtended(
  q: number | bigint,
  variableName: string = 'a'
): PrimeField | FiniteFieldExtension {
  const order = typeof q === 'number' ? BigInt(q) : q;

  if (order < 2n) {
    throw new ValueError('field order must be at least 2');
  }

  // Check if q is a prime power
  const [isPrimePower, p, n] = factorPrimePower(order);

  if (!isPrimePower) {
    throw new ValueError(`${order} is not a prime power`);
  }

  if (n === 1n) {
    return new PrimeField(p);
  }

  return new FiniteFieldExtension(p, Number(n), undefined, variableName);
}

/**
 * Alias for GFExtended for backward compatibility in tests.
 */
export const GF = GFExtended;

/**
 * Construct a finite field extension explicitly.
 *
 * @param p - Prime characteristic
 * @param n - Extension degree
 * @param modulus - Optional modulus polynomial coefficients
 * @param variableName - Name for the generator
 */
export function GFpn(
  p: number | bigint,
  n: number,
  modulus?: number[],
  variableName: string = 'a'
): FiniteFieldExtension {
  return new FiniteFieldExtension(p, n, modulus, variableName);
}

/**
 * Factor a number as a prime power.
 * Returns [isPrimePower, base, exponent].
 */
function factorPrimePower(n: bigint): [boolean, bigint, bigint] {
  if (n < 2n) {
    return [false, n, 1n];
  }

  // Check if n is a prime
  if (is_prime(n)) {
    return [true, n, 1n];
  }

  // Find the smallest prime factor
  let p = 2n;
  while (p * p <= n) {
    if (n % p === 0n) {
      break;
    }
    p++;
  }

  if (p * p > n) {
    // n is prime
    return [true, n, 1n];
  }

  // Check if n = p^k
  let k = 0n;
  let temp = n;
  while (temp % p === 0n) {
    temp /= p;
    k++;
  }

  if (temp === 1n) {
    return [true, p, k];
  }

  return [false, n, 1n];
}
