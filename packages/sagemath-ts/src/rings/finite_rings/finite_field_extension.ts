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

import { factor, inverse_mod, is_prime, power_mod, primitive_root } from '../../arith/misc.js';
import { ValueError, ZeroDivisionError } from '../../errors.js';
import { current_randstate } from '../../misc/randstate.js';
import {
  type CoefficientRing,
  Polynomial,
  type RingElement,
} from '../polynomial/polynomial_element.js';
import { PolynomialRing } from '../polynomial/polynomial_ring.js';
import { QuotientRing, QuotientRingElement } from '../polynomial/quotient_ring.js';
import { Integer } from '../integer_ring.js';
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

  /**
   * Reduce an operand to its canonical representative in [0, p).
   *
   * SageMath coerces plain integers into GF(p) automatically, so `3 * x` is
   * valid there, and the `FieldElement` contract that the elliptic-curve code
   * programs against (`schemes/elliptic_curves/types.ts`) declares
   * `add`/`sub`/`mul`/`div` as accepting `FieldElement | number | bigint`. The
   * sibling `FiniteFieldElement` (`finite_field_prime.ts`) already coerces.
   * Without this, callers reaching a `PrimeFieldElement` through the
   * `FieldElement` interface -- e.g. Velu's formulas in
   * `ell_curve_isogeny.ts`, which write `xQ.mul(xQ).mul(3)` -- threw
   * `TypeError: Invalid mix of BigInt and other type in multiplication`.
   *
   * The declared parameter types below stay `PrimeFieldElement` on purpose:
   * `RingElement` (`rings/polynomial/polynomial_element.ts`) specifies
   * `add(other: this): this`, and widening the declaration would stop
   * `PrimeFieldElement` satisfying the `C extends RingElement` constraint used
   * by `Polynomial<C>`. The coercion is therefore applied in the body only.
   */
  private _coerceValue(other: PrimeFieldElement | number | bigint): bigint {
    if (other instanceof PrimeFieldElement) {
      return other.value;
    }
    const v = typeof other === 'number' ? BigInt(other) : other;
    return (
      ((v % this.parent.characteristic) + this.parent.characteristic) % this.parent.characteristic
    );
  }

  add(other: PrimeFieldElement): PrimeFieldElement {
    return new PrimeFieldElement(
      (this.value + this._coerceValue(other)) % this.parent.characteristic,
      this.parent
    );
  }

  sub(other: PrimeFieldElement): PrimeFieldElement {
    return new PrimeFieldElement(
      (((this.value - this._coerceValue(other)) % this.parent.characteristic) +
        this.parent.characteristic) %
        this.parent.characteristic,
      this.parent
    );
  }

  mul(other: PrimeFieldElement): PrimeFieldElement {
    return new PrimeFieldElement(
      (this.value * this._coerceValue(other)) % this.parent.characteristic,
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
    const d =
      other instanceof PrimeFieldElement
        ? other
        : new PrimeFieldElement(this._coerceValue(other), this.parent);
    return this.mul(d.inv());
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

  /**
   * Return a generator of this field over its prime field, i.e. a root of the
   * modulus.  For GF(p) with the default modulus x - 1 this is `1`.
   *
   * This is **not** a generator of the multiplicative group; use
   * {@link multiplicative_generator} for that.
   *
   * Port of `sage/rings/finite_rings/finite_field_prime_modn.py:gen`
   * (`sage: GF(13).gen()` -> `1`).
   */
  gen(): PrimeFieldElement {
    return this.one();
  }

  /**
   * Find a primitive root modulo p.
   *
   * Sage's `multiplicative_generator` for a degree-1 field is
   * `self(primitive_root(self.order()))`
   * (`finite_field_base.pyx:723-725`), so delegate to `arith.primitive_root`.
   */
  private primitiveRoot(): PrimeFieldElement {
    return new PrimeFieldElement(primitive_root(this.characteristic), this);
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

  /**
   * Alias for {@link multiplicative_generator}
   * (`finite_field_base.pyx:729`: `primitive_element = multiplicative_generator`).
   */
  primitive_element(): PrimeFieldElement {
    return this.multiplicative_generator();
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
   *
   * @see Deviation: Finite Field Extension Minimal Polynomial Simplified
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
 *
 * @see Deviation: Conway Polynomial Database Limited
 */
export class FiniteFieldExtension implements CoefficientRing<FiniteFieldElement> {
  readonly baseField: PrimeField;
  readonly polynomialRing: PolynomialRing<PrimeFieldElement>;
  readonly modulus: Polynomial<PrimeFieldElement>;
  readonly degree: number;
  readonly characteristic: bigint;
  readonly order: bigint;
  readonly variableName: string;


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
   * This is SageMath's `algorithm='first_lexicographic'`
   * (`sage/rings/polynomial/polynomial_ring.py:2677-2681`):
   *
   *     for g in self.polynomials(max_degree=n-1):
   *         f = self.gen()**n + g
   *         if f.is_irreducible():
   *             return f
   *
   * `polynomials(max_degree=d)` enumerates by `_polys_max`
   * (`polynomial_ring.py:1548-1557`), i.e. the constant term varies fastest —
   * so `g` runs through `0, 1, ..., p-1, x, x+1, ...`, exactly the base-`p`
   * counter used below.  Sage's doctest
   * `GF(19)['x'].irreducible_element(21, algorithm='first_lexicographic')`
   * gives `x^21 + x + 5`, which this reproduces.
   *
   * Sage's *default* here would be `pari(p).ffinit(n)` (Adleman-Lenstra), or
   * NTL's `GF2X_BuildSparseIrred` for p = 2; neither is available in our
   * PARI/NTL ports yet, so this deterministic search is used instead of the
   * former random search (which could fail outright, and depended on the
   * random state).
   *
   * @see Deviation: Irreducible Modulus Search
   */
  private findIrreducible(n: number): Polynomial<PrimeFieldElement> {
    const x = this.polynomialRing.gen();
    const xPowN = x.pow(n);
    const p = this.characteristic;

    // Number of monic candidates of degree n is p^n; an irreducible one always
    // exists, so this loop terminates well before the bound.
    const bound = p ** BigInt(n);
    for (let k = 0n; k < bound; k++) {
      const coeffs: PrimeFieldElement[] = [];
      let temp = k;
      for (let j = 0; j < n; j++) {
        coeffs.push(this.baseField.__call__(temp % p));
        temp /= p;
      }
      const candidate = xPowN.add(new Polynomial(coeffs, this.polynomialRing));
      if (this.isIrreducible(candidate)) {
        return candidate;
      }
    }

    // Unreachable: monic irreducible polynomials of every degree exist over
    // every finite field.
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

      // gcd(f, 0) = f, so a zero difference must be rejected here: it means f
      // divides x^{p^k} - x outright, i.e. f splits into distinct factors of
      // degree dividing k < n. FLINT guards only its make_monic call on this
      // condition (nmod_poly_factor/is_irreducible.c:238), never the gcd test.
      const g = this.polyGcd(f, diff);
      if (g.degree() > 0) {
        return false; // f has a factor of degree <= k
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
    // Lazy: the constant elements come first (as they do from PARI's finite
    // field iterator, which Sage relies on in `_element_of_factored_order`),
    // and nothing is materialised, so consuming a prefix of a huge field is
    // cheap.
    const p = this.characteristic;

    for (let i = 0n; i < this.order; i++) {
      const coeffs: PrimeFieldElement[] = [];
      let temp = i;

      for (let j = 0; j < this.degree; j++) {
        coeffs.push(this.baseField.__call__(temp % p));
        temp /= p;
      }

      const poly = new Polynomial(coeffs, this.polynomialRing);
      yield new FiniteFieldElement(poly, this);
    }
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
   *
   * Port of `finite_field_base.pyx:731-778` (`_element_of_factored_order`,
   * called by `multiplicative_generator`) with `n = self.order() - 1`, hence
   * `c = 1`: we test `g + x` for `x` running through the field, starting with
   * `x = 0`, so a Conway modulus (whose root is primitive by construction)
   * returns immediately.
   */
  primitiveElement(): FiniteFieldElement {
    // The multiplicative group has order p^n - 1
    const groupOrder = this.order - 1n;
    const primes = factorSimple(groupOrder).map(([p]) => p);

    const g = this.gen();
    for (const x of this) {
      const a = g.add(x);
      if (a.isZero()) continue;
      if (primes.every((p) => !a.pow(groupOrder / p).isOne())) {
        return a;
      }
    }

    throw new ValueError('no element found');
  }

  /**
   * Alias for {@link primitiveElement}, matching Sage's spelling
   * (`finite_field_base.pyx:689`).
   */
  multiplicative_generator(): FiniteFieldElement {
    return this.primitiveElement();
  }

  /**
   * Alias for {@link primitiveElement} (`finite_field_base.pyx:729`).
   */
  primitive_element(): FiniteFieldElement {
    return this.primitiveElement();
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
 *
 * Delegates to `arith.is_prime_power` (which delegates to PARI) instead of
 * trial dividing up to sqrt(n), and re-verifies the answer with `is_prime`
 * and `p^k === n` so that a composite PARI fails to split cannot be mistaken
 * for a prime power.
 */
function factorPrimePower(n: bigint): [boolean, bigint, bigint] {
  if (n < 2n) {
    return [false, n, 1n];
  }

  if (is_prime(n)) {
    return [true, n, 1n];
  }

  // n is composite, so n = p^k needs k >= 2.  Try every exponent from the
  // largest possible one downwards: the first exact k-th root of a genuine
  // prime power p^k is reached at k, where the root is p itself.
  const maxExp = BigInt(n.toString(2).length - 1);
  for (let k = maxExp; k >= 2n; k--) {
    const [root, exact] = new Integer(n).nth_root(k, true);
    if (exact && is_prime(root.value)) {
      return [true, root.value, k];
    }
  }

  return [false, n, 1n];
}

/**
 * Prime factorization used by the primitive-element search.
 *
 * Delegates to `arith.factor` (PARI's `Z_factor`) rather than trial dividing
 * to sqrt(n), which was hopeless for p^n - 1 with p^n of cryptographic size.
 */
function factorSimple(n: bigint): Array<[bigint, bigint]> {
  return factor(n).filter(([p]) => p > 0n);
}
