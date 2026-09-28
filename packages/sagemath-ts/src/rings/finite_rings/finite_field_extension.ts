import { isFractionElement } from '../fraction_field_element.js';
import { GF2Element } from './gf2.js';
/**
 * @module sage/rings/finite_rings/finite_field_extension
 * @description Finite field extensions GF(p^n) using polynomial quotient ring construction
 * @see Deviation: Finite Field Coercion and Backend Boundaries
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

import { GF2X_BuildIrred, GF2X_BuildSparseIrred } from '@sagemath-ts/ntl-ts';
import {
  FF_issquareall,
  FF_trace,
  FF_charpoly,
  FF_issquare,
  PariType,
  FpXQ_inv,
  FpXQ_minpoly,
  FpXQ_mul,
  FpXQ_pow,
  FpX_add,
  FpX_neg,
  FpX_sub,
  Z_isanypower,
  ffinit,
} from '@sagemath-ts/parigp-ts';
import {
  factor,
  inverse_mod,
  is_prime,
  power_mod,
  primitive_root,
  sqrt_mod,
} from '../../arith/misc.js';
import {
  ArithmeticError,
  IndexError,
  NotImplementedError,
  ValueError,
  ZeroDivisionError,
} from '../../errors.js';
import { current_randstate } from '../../misc/randstate.js';
import type { IntegerLike } from '../../types/coercion.js';
import { Integer, ZZ } from '../integer_ring.js';
import {
  type CoefficientRing,
  Polynomial,
  type RingElement,
} from '../polynomial/polynomial_element.js';
import { PolynomialRing } from '../polynomial/polynomial_ring.js';
import { QuotientRing, QuotientRingElement } from '../polynomial/quotient_ring.js';
import { Rational } from '../rational.js';
import { conway_polynomial, has_conway_polynomial } from './conway_polynomials.js';
import { FiniteFieldElement as LegacyPrimeElement } from './finite_field_prime.js';
import {
  IntegerMod,
  canonicalFiniteOperands,
  checkFiniteGeneratorIndex,
  finiteArithmeticEquals,
  repeatFiniteSequence,
} from './integer_mod.js';
import { IntegerModRing } from './integer_mod_ring.js';

/**
 * Element of a prime field GF(p).
 */
export class PrimeFieldElement implements RingElement {
  /** IntegerMod._rational_: lift the canonical residue into QQ. */
  _rational_(): Rational {
    return new Rational(this.value);
  }

  readonly value: bigint;
  readonly parent: PrimeField;

  constructor(value: unknown, parent: PrimeField) {
    this.parent = parent;

    const v = typeof value === 'bigint' ? value : parent.__call__(value).value;
    this.value = ((v % parent.characteristic) + parent.characteristic) % parent.characteristic;
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
   * Arithmetic requires compatible parents; explicit ring construction can still
   * convert the integer representative from another characteristic.
   */
  add(other: FiniteFieldElement): FiniteFieldElement;
  add(
    other: PrimeFieldElement | LegacyPrimeElement | IntegerMod | Integer | boolean | number | bigint
  ): PrimeFieldElement;
  add(
    other:
      | PrimeFieldElement
      | LegacyPrimeElement
      | FiniteFieldElement
      | IntegerMod
      | Integer
      | boolean
      | number
      | bigint
  ): PrimeFieldElement | FiniteFieldElement;
  add(other: PrimeFieldElement): PrimeFieldElement;
  add(other: unknown): PrimeFieldElement | FiniteFieldElement {
    const [left, right] = canonicalFiniteOperands(this, other, '+');
    if (left instanceof FiniteFieldElement) return left.add(right as FiniteFieldElement);
    const operand = right as PrimeFieldElement;
    return new PrimeFieldElement(this.value + operand.value, this.parent);
  }

  sub(other: FiniteFieldElement): FiniteFieldElement;
  sub(
    other: PrimeFieldElement | LegacyPrimeElement | IntegerMod | Integer | boolean | number | bigint
  ): PrimeFieldElement;
  sub(
    other:
      | PrimeFieldElement
      | LegacyPrimeElement
      | FiniteFieldElement
      | IntegerMod
      | Integer
      | boolean
      | number
      | bigint
  ): PrimeFieldElement | FiniteFieldElement;
  sub(other: PrimeFieldElement): PrimeFieldElement;
  sub(other: unknown): PrimeFieldElement | FiniteFieldElement {
    const [left, right] = canonicalFiniteOperands(this, other, '-');
    if (left instanceof FiniteFieldElement) return left.sub(right as FiniteFieldElement);
    const operand = right as PrimeFieldElement;
    return new PrimeFieldElement(this.value - operand.value, this.parent);
  }

  mul(other: string): string;
  mul<T>(other: readonly T[]): T[];
  mul(other: FiniteFieldElement): FiniteFieldElement;
  mul(
    other: PrimeFieldElement | LegacyPrimeElement | IntegerMod | Integer | boolean | number | bigint
  ): PrimeFieldElement;
  mul(
    other:
      | PrimeFieldElement
      | LegacyPrimeElement
      | FiniteFieldElement
      | IntegerMod
      | Integer
      | boolean
      | number
      | bigint
  ): PrimeFieldElement | FiniteFieldElement;
  mul(other: PrimeFieldElement): PrimeFieldElement;
  mul(other: unknown): PrimeFieldElement | FiniteFieldElement | string | unknown[] {
    if (typeof other === 'string' || Array.isArray(other))
      return repeatFiniteSequence(this.value, other);
    const [left, right] = canonicalFiniteOperands(this, other, '*');
    if (left instanceof FiniteFieldElement) return left.mul(right as FiniteFieldElement);
    const operand = right as PrimeFieldElement;
    return new PrimeFieldElement(this.value * operand.value, this.parent);
  }

  neg(): PrimeFieldElement {
    if (this.value === 0n) {
      return this;
    }
    return new PrimeFieldElement(this.parent.characteristic - this.value, this.parent);
  }

  inv(): PrimeFieldElement {
    if (this.value === 0n) {
      throw new ZeroDivisionError(
        `inverse of Mod(0, ${this.parent.characteristic}) does not exist`
      );
    }
    return new PrimeFieldElement(inverse_mod(this.value, this.parent.characteristic), this.parent);
  }

  div(other: FiniteFieldElement): FiniteFieldElement;
  div(
    other: PrimeFieldElement | LegacyPrimeElement | IntegerMod | Integer | boolean | number | bigint
  ): PrimeFieldElement;
  div(
    other:
      | PrimeFieldElement
      | LegacyPrimeElement
      | FiniteFieldElement
      | IntegerMod
      | Integer
      | boolean
      | number
      | bigint
  ): PrimeFieldElement | FiniteFieldElement;
  div(other: PrimeFieldElement): PrimeFieldElement;
  div(other: unknown): PrimeFieldElement | FiniteFieldElement {
    const [left, right] = canonicalFiniteOperands(this, other, '/');
    if (left instanceof FiniteFieldElement) return left.div(right as FiniteFieldElement);
    const operand = right as PrimeFieldElement;
    return this.mul(operand.inv());
  }

  pow(n: IntegerLike | number | Rational | boolean | string | null): PrimeFieldElement {
    // Sage prime-field elements inherit IntegerMod's native/GMP power dispatch.
    const result = new IntegerModRing(this.parent.characteristic).__call__(this.value).pow(n);
    return new PrimeFieldElement(result.value, this.parent);
  }

  eq(other: unknown): boolean {
    return finiteArithmeticEquals(this, other);
  }

  isZero(): boolean {
    return this.value === 0n;
  }

  isOne(): boolean {
    return this.value === 1n;
  }

  /**
   * Return whether this element is a square.
   *
   * This is the degree-one branch of Sage's finite-field element
   * `is_square()`: zero and every element of GF(2) are squares; over an odd
   * prime field Euler's criterion decides quadratic residuosity.
   */
  is_square(): boolean {
    if (this.value === 0n || this.parent.characteristic === 2n) {
      return true;
    }
    return this.pow((this.parent.characteristic - 1n) / 2n).isOne();
  }

  /**
   * Return one square root, or all roots when requested.
   *
   * Sage's prime-field default extends to GF(p^2) for a nonsquare. The all-roots
   * extension branch preserves integer_mod.pyx's NotImplementedError and its
   * small-modulus search shortcut. Base-field roots use the PARI-backed sqrt_mod.
   * @see Deviation: Finite Field Coercion and Backend Boundaries
   */
  sqrt(options: { all: true; extend?: boolean }): PrimeFieldElement[];
  sqrt(options: { extend: false; all?: false }): PrimeFieldElement;
  sqrt(options?: { extend?: boolean; all?: false }): PrimeFieldElement | FiniteFieldElement;
  sqrt(options: { extend?: boolean; all?: boolean }):
    | PrimeFieldElement
    | FiniteFieldElement
    | PrimeFieldElement[];
  sqrt(
    options: { extend?: boolean; all?: boolean } = {}
  ): PrimeFieldElement | FiniteFieldElement | PrimeFieldElement[] {
    for (const key of Object.keys(options)) {
      if (key !== 'extend' && key !== 'all')
        throw new TypeError(`sqrt() got an unexpected keyword argument '${key}'`);
    }
    const p = this.parent.characteristic;
    const root = sqrt_mod(this.value, p);
    if (root !== null) {
      const r = root <= p - root ? root : p - root;
      const x = this.parent.__call__(r);
      return options.all ? (r === 0n || p === 2n ? [x] : [x, this.parent.__call__(p - r)]) : x;
    }
    const extend = options.extend ?? true;
    // integer_mod.pyx:3045-3076: the int32 all-roots search returns []
    // immediately for nonsquares, even if extend=True. For larger p=3 mod 4
    // the preceding optimized branch takes priority and falls through instead.
    const smallSearch = p <= 100n || (p < 10000n && p % 4n !== 3n);
    if (options.all && (!extend || smallSearch)) return [];
    if (!extend) throw new ValueError('self must be a square');
    if (options.all) {
      throw new NotImplementedError(
        'Finding all square roots in extensions is not implemented; try extend=False to find only roots in the base ring Zmod(n).'
      );
    }
    const ring = new PolynomialRing(this.parent, 'x');
    const modulus = new Polynomial<PrimeFieldElement>(
      [this.neg(), this.parent.zero(), this.parent.one()],
      ring
    );
    return new FiniteFieldExtension(p, 2, modulus, `sqrt${this.value}`).gen();
  }

  toString(): string {
    return this.value.toString();
  }

  repr(): string {
    return this.value.toString();
  }

  /** IntegerMod._integer_; return the canonical integer representative. */
  _integer_(_ZZ?: unknown): bigint {
    return this.value;
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
  private readonly _generator: bigint;

  constructor(
    p: unknown,
    options?: { modulus?: Polynomial<PrimeFieldElement> | FiniteFieldModulusAlgorithm }
  ) {
    const prime = ZZ.__call__(p as Parameters<typeof ZZ.__call__>[0]);

    if (!is_prime(prime)) {
      throw new ArithmeticError('p must be prime');
    }

    this.characteristic = prime;
    this.order = prime;
    // finite_field_prime_modn.py:gen returns -modulus[0] for a custom
    // degree-one modulus. Reuse the same modulus normalization/algorithm path.
    this._generator =
      options?.modulus === undefined
        ? 1n
        : new FiniteFieldExtension(prime, 1, options.modulus).gen().lift.getCoeff(0).value;
  }

  /** Coerce through IntegerMod.__init__ (integer_mod.pyx:377-404). */
  __call__(x?: unknown): PrimeFieldElement {
    if (
      isFractionElement(x) ||
      (typeof x === 'object' && x !== null && 'coeffs' in x && 'getCoeff' in x)
    )
      return new PrimeFieldElement(new IntegerModRing(this.characteristic).__call__(x).value, this);
    if (x instanceof PrimeFieldElement) {
      return new PrimeFieldElement(x.value, this);
    }
    if (x instanceof FiniteFieldElement) {
      if (x.lift.degree() > 0) {
        if (x.parent.characteristic === this.characteristic)
          throw new ValueError(
            `${x} is not in the image of (map internal to coercion system -- copy before use)\nRing morphism:\n  From: ${this}\n  To:   ${x.parent}`
          );
        throw new TypeError(`unable to convert ${x} to a rational`);
      }
      return new PrimeFieldElement(x.lift.getCoeff(0).value, this);
    }
    if (x instanceof Rational) {
      const d = ((x.denominator % this.characteristic) + this.characteristic) % this.characteristic;
      if (d === 0n) {
        throw new ZeroDivisionError(`inverse of Mod(0, ${this.characteristic}) does not exist`);
      }
      return new PrimeFieldElement(x.numerator * inverse_mod(d, this.characteristic), this);
    }
    return new PrimeFieldElement(ZZ.__call__(x as Parameters<typeof ZZ.__call__>[0]), this);
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
  gen(n: unknown = 0): PrimeFieldElement {
    checkFiniteGeneratorIndex(n, true);
    return this.__call__(this._generator);
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
    const randomInt = rstate.python_random().randrange(this.characteristic);
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

    // element_pari_ffelt.pyx:492-496 changes coefficient rings before
    // evaluating at the generator. A low-degree foreign polynomial must not
    // retain coefficients from its old characteristic.
    if (poly.parent !== parent.polynomialRing) {
      poly = new Polynomial<PrimeFieldElement>(
        poly.coeffs.map((c) => parent.baseField.__call__(c)),
        parent.polynomialRing
      );
    }
    // Reduce modulo the modulus
    if (poly.degree() >= parent.modulus.degree()) {
      const [_q, r] = poly.quo_rem(parent.modulus);
      this.lift = r;
    } else {
      this.lift = poly;
    }
  }

  /** Test squareness through the bundled PARI FF_issquare predicate.
   * @see Reference: sage/rings/finite_rings/element_pari_ffelt.pyx:is_square
   * @see Deviation: PARI finite-field norm and square-predicate adapters
   */
  is_square(): boolean {
    return FF_issquare({type:PariType.t_FFELT,p:this.parent.characteristic,degree:this.parent.degree,
      value:this._pariCoefficients(),definingPoly:this._pariModulus()});
  }

  /** Return a PARI-selected root, or all roots in native [r,-r] order.
   * The default extend=false and unsupported extend=true follow pari_ffelt.
   * @see Reference: sage/rings/finite_rings/element_pari_ffelt.pyx:sqrt
   * @see Deviation: PARI finite-field scalar square-root adapter
   */
  sqrt(options: { all: true; extend?: boolean }): FiniteFieldElement[];
  sqrt(options?: { all?: false; extend?: boolean }): FiniteFieldElement;
  sqrt(options: { all?: boolean; extend?: boolean }): FiniteFieldElement | FiniteFieldElement[];
  sqrt(options: { all?: boolean; extend?: boolean } = {}): FiniteFieldElement | FiniteFieldElement[] {
    for (const key of Object.keys(options)) {
      if (key !== 'extend' && key !== 'all')
        throw new TypeError(`sqrt() got an unexpected keyword argument '${key}'`);
    }
    if (options.extend) throw new NotImplementedError('');
    const root = FF_issquareall({
      type: PariType.t_FFELT, p: this.parent.characteristic,
      degree: this.parent.degree, value: this._pariCoefficients(),
      definingPoly: this._pariModulus(),
    });
    if (root === null) {
      if (options.all) return [];
      throw new ValueError('element is not a square');
    }
    const value = this._fromPari(typeof root.value === 'bigint' ? [root.value] : [...root.value]);
    return options.all
      ? value.isZero() || this.parent.characteristic === 2n ? [value] : [value, value.neg()]
      : value;
  }

  /** PARI GEN polynomial coefficients, in ascending degree order. */
  private _pariCoefficients(): bigint[] {
    return this.lift.coeffs.map((coefficient) => coefficient.value);
  }

  private _pariModulus(): bigint[] {
    return this.parent.modulus.coeffs.map((coefficient) => coefficient.value);
  }

  private _fromPari(coefficients: bigint[]): FiniteFieldElement {
    return new FiniteFieldElement(
      this.parent.polynomialRing.__call__(
        coefficients.map((c) => this.parent.baseField.__call__(c))
      ),
      this.parent
    );
  }

  /** Canonical arithmetic coercion; explicit __call__ conversions are separate. */
  private _coerceOperand(other: unknown, operation: string): FiniteFieldElement {
    return canonicalFiniteOperands(this, other, operation)[1] as FiniteFieldElement;
  }

  add(
    other:
      | FiniteFieldElement
      | PrimeFieldElement
      | LegacyPrimeElement
      | IntegerMod
      | Integer
      | boolean
      | number
      | bigint
  ): FiniteFieldElement;
  add(other: FiniteFieldElement): FiniteFieldElement;
  add(
    other:
      | FiniteFieldElement
      | PrimeFieldElement
      | LegacyPrimeElement
      | IntegerMod
      | Integer
      | boolean
      | number
      | bigint
  ): FiniteFieldElement {
    const operand = this._coerceOperand(other, '+');
    return this._fromPari(
      FpX_add(this._pariCoefficients(), operand._pariCoefficients(), this.parent.characteristic)
    );
  }

  sub(
    other:
      | FiniteFieldElement
      | PrimeFieldElement
      | LegacyPrimeElement
      | IntegerMod
      | Integer
      | boolean
      | number
      | bigint
  ): FiniteFieldElement;
  sub(other: FiniteFieldElement): FiniteFieldElement;
  sub(
    other:
      | FiniteFieldElement
      | PrimeFieldElement
      | LegacyPrimeElement
      | IntegerMod
      | Integer
      | boolean
      | number
      | bigint
  ): FiniteFieldElement {
    const operand = this._coerceOperand(other, '-');
    return this._fromPari(
      FpX_sub(this._pariCoefficients(), operand._pariCoefficients(), this.parent.characteristic)
    );
  }

  mul(
    other:
      | FiniteFieldElement
      | PrimeFieldElement
      | LegacyPrimeElement
      | IntegerMod
      | Integer
      | boolean
      | number
      | bigint
  ): FiniteFieldElement;
  mul(other: FiniteFieldElement): FiniteFieldElement;
  mul(
    other:
      | FiniteFieldElement
      | PrimeFieldElement
      | LegacyPrimeElement
      | IntegerMod
      | Integer
      | boolean
      | number
      | bigint
  ): FiniteFieldElement {
    const operand = this._coerceOperand(other, '*');
    return this._fromPari(
      FpXQ_mul(
        this._pariCoefficients(),
        operand._pariCoefficients(),
        this._pariModulus(),
        this.parent.characteristic
      )
    );
  }

  neg(): FiniteFieldElement {
    return this._fromPari(FpX_neg(this._pariCoefficients(), this.parent.characteristic));
  }

  /** Sage FF_inv's FpXQ backend; reject zero before entering PARI. */
  inv(): FiniteFieldElement {
    if (this.isZero()) throw new ZeroDivisionError('');
    return this._fromPari(
      FpXQ_inv(this._pariCoefficients(), this._pariModulus(), this.parent.characteristic)
    );
  }

  div(
    other:
      | FiniteFieldElement
      | PrimeFieldElement
      | LegacyPrimeElement
      | IntegerMod
      | Integer
      | boolean
      | number
      | bigint
  ): FiniteFieldElement;
  div(other: FiniteFieldElement): FiniteFieldElement;
  div(
    other:
      | FiniteFieldElement
      | PrimeFieldElement
      | LegacyPrimeElement
      | IntegerMod
      | Integer
      | boolean
      | number
      | bigint
  ): FiniteFieldElement {
    const operand = this._coerceOperand(other, '/');
    return this.mul(operand.inv());
  }

  /**
   * element_pari_ffelt.pyx:826-868 compares the exponent before Integer coercion.
   * The field's PARI quotient kernel then performs signed exponentiation.
   * @see Deviation: Extension Arithmetic and PARI Quotient Kernels
   */
  pow(n: unknown): FiniteFieldElement {
    const finite =
      n instanceof IntegerMod ||
      n instanceof PrimeFieldElement ||
      n instanceof LegacyPrimeElement ||
      n instanceof FiniteFieldElement;
    if (
      n === 0 ||
      n === 0n ||
      n === false ||
      (n instanceof Integer && n.value === 0n) ||
      (n instanceof Rational && n.numerator === 0n) ||
      (finite && n.isZero())
    )
      return this.parent.one();
    if (
      !(
        typeof n === 'number' ||
        typeof n === 'bigint' ||
        typeof n === 'boolean' ||
        n instanceof Integer ||
        n instanceof Rational ||
        finite
      )
    ) {
      const type =
        n == null
          ? 'NoneType'
          : typeof n === 'string'
            ? 'str'
            : Array.isArray(n)
              ? 'list'
              : 'object';
      throw new TypeError(`'<' not supported between instances of '${type}' and 'int'`);
    }
    const negative =
      n instanceof Rational
        ? n.numerator < 0n
        : n instanceof Integer
          ? n.value < 0n
          : (typeof n === 'bigint' || typeof n === 'number') && n < 0;
    if (negative && this.isZero()) throw new ZeroDivisionError('');
    const exponent = ZZ.__call__(n as Parameters<typeof ZZ.__call__>[0]);
    return this._fromPari(
      FpXQ_pow(this._pariCoefficients(), exponent, this._pariModulus(), this.parent.characteristic)
    );
  }

  /**
   * Apply the Frobenius automorphism: x -> x^p
   *
   * This is a field automorphism of GF(p^n) that fixes GF(p).
   */
  frobenius(power: number = 1): FiniteFieldElement {
    const p = this.parent.characteristic;
    const degree = BigInt(this.parent.degree);
    const reducedPower = ((BigInt(power) % degree) + degree) % degree;
    const exp = p ** reducedPower;
    return this.pow(exp);
  }

  /**
   * Compute the prime-subfield trace through PARI FF_trace, as in Sage.
   * @see Deviation: PARI finite-field trace adapters
   *
   * Returns an element of the base field GF(p).
   */
  trace(): PrimeFieldElement {
    return this.parent.baseField.__call__(FF_trace({
      type: PariType.t_FFELT,
      p: this.parent.characteristic,
      degree: this.parent.degree,
      value: this._pariCoefficients(),
      definingPoly: this._pariModulus(),
    }));
  }

  /** Characteristic polynomial over the prime subfield (element_pari_ffelt.pyx:982).
   * @see Deviation: PARI bivariate polynomial storage
   */
  charpoly(varName: string = 'x'): Polynomial<PrimeFieldElement> {
    const coefficients = FF_charpoly({
      type: PariType.t_FFELT,
      p: this.parent.characteristic,
      degree: this.parent.degree,
      value: this._pariCoefficients(),
      definingPoly: this._pariModulus(),
    });
    return new Polynomial(
      coefficients.map((c) => this.parent.baseField.__call__(c)),
      new PolynomialRing(this.parent.baseField, varName)
    );
  }

  /** Prime-subfield norm, from the signed constant term of charpoly (element_base.pyx). */
  norm(): PrimeFieldElement {
    const f = this.charpoly('x'),
      n = f.getCoeff(0);
    return f.degree() % 2 ? n.neg() : n;
  }

  /**
   * Compute the minimal polynomial of this element over GF(p).
   *
   * Port of `FinitePolyExtElement.minpoly`
   * (`sage/rings/finite_rings/element_base.pyx:204`) and
   * `FiniteFieldElement_pari_ffelt.minpoly`
   * (`element_pari_ffelt.pyx:963`).
   *
   * Sage's default implementation delegates to PARI's `FF_minpoly`, whose
   * prime-field branch calls `FpXQ_minpoly` (`pari/src/basemath/FF.c:1045`).
   * Delegate to the same routine in parigp-ts, passing the coefficient vectors
   * of this element and the defining modulus in constant-term-first order.
   *
   * The optional `algorithm` argument is retained for compatibility with
   * `FinitePolyExtElement`: both of Sage's accepted algorithms compute the
   * same canonical monic polynomial. The port's matrix layer does not expose
   * the left-multiplication matrix of a finite-field element, so both accepted
   * values use the already-ported PARI routine.
   */
  minpoly(varName: string = 'x', algorithm: string = 'pari'): Polynomial<PrimeFieldElement> {
    if (algorithm !== 'pari' && algorithm !== 'matrix') {
      throw new ValueError(`unknown algorithm '${algorithm}'`);
    }

    const element = this.lift.coeffs.map((coefficient) => coefficient.value);
    const modulus = this.parent.modulus.coeffs.map((coefficient) => coefficient.value);
    const coefficients = FpXQ_minpoly(element, modulus, this.parent.characteristic);
    const polynomialRing = new PolynomialRing(this.parent.baseField, varName);

    return new Polynomial(
      coefficients.map((coefficient) => this.parent.baseField.__call__(coefficient)),
      polynomialRing
    );
  }

  /**
   * Sage-compatible alias for {@link minpoly}.
   */
  minimal_polynomial(varName: string = 'x'): Polynomial<PrimeFieldElement> {
    return this.minpoly(varName);
  }

  /**
   * Backwards-compatible camelCase alias used by earlier sagemath-ts releases.
   */
  minimalPolynomial(
    varName: string = 'x',
    algorithm: string = 'pari'
  ): Polynomial<PrimeFieldElement> {
    return this.minpoly(varName, algorithm);
  }

  eq(other: unknown): boolean {
    return finiteArithmeticEquals(this, other);
  }

  isZero(): boolean {
    return this.lift.isZero();
  }

  isOne(): boolean {
    return this.lift.coeffs.length === 1 && this.lift.coeffs[0]!.isOne();
  }

  /** Lift a prime-subfield element to an integer; the optional parent is ignored. */
  _integer_(_ZZ?: unknown): bigint {
    // element_pari_ffelt.pyx:1259-1274 delegates to the prime-subfield lift.
    if (this.lift.degree() > 0) throw new ValueError('element is not in the prime field');
    return this.lift.getCoeff(0).value;
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
 * The `algorithm` strings SageMath's `irreducible_element` accepts
 * (`reference/sage/src/sage/rings/polynomial/polynomial_ring.py:3560`), which
 * `GF(q, modulus=<string>)` forwards
 * (`finite_field_constructor.py:729-734`).
 */
export type FiniteFieldModulusAlgorithm =
  | 'conway'
  | 'adleman-lenstra'
  | 'primitive'
  | 'first_lexicographic'
  | 'minimal_weight'
  | 'ffprimroot'
  | 'random';

/**
 * Coefficient list of an NTL `GF2X`, constant term first, padded to `n + 1`
 * entries -- the shape of Sage's `GF2X_Build*Irred_list`
 * (`polynomial_gf2x.pyx:296`: `[GF2(...) for i in range(n + 1)]`).
 */
function gf2xToCoeffs(f: { coeff(i: number): { rep(): number } }, n: number): number[] {
  const out: number[] = [];
  for (let i = 0; i <= n; i++) {
    out.push(f.coeff(i).rep());
  }
  return out;
}

/**
 * A finite field extension GF(p^n).
 *
 * Constructed as GF(p)[x] / <f(x)> where f(x) is an irreducible polynomial.
 *
 * @see Deviation: Finite Fields — Conway Table and Constructor Algorithms
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
    modulus?: Polynomial<PrimeFieldElement> | number[] | FiniteFieldModulusAlgorithm,
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

    if (modulus !== undefined) {
      if (modulus instanceof Polynomial) {
        const polynomial = this.polynomialRing.__call__(modulus);
        if (polynomial.isZero()) {
          // NTL's generic monic path in characteristic two inverts zero;
          // Executed Sage 10.3 FLINT reports a noninvertible lead; newer
          // FLINT's exception mapping remains a documented version boundary.
          if (prime === 2n) throw new ZeroDivisionError('inverse of Mod(0, 2) does not exist');
          throw new ValueError('leading coefficient must be invertible');
        }
        this.modulus = polynomial.monic();
      } else if (typeof modulus === 'string') {
        // Sage: `if isinstance(modulus, str): modulus = R.irreducible_element(n, algorithm=modulus)`
        // (`finite_field_constructor.py:729-734`).
        this.modulus = this.irreducible_element(n, modulus);
      } else {
        // modulus is an array of coefficients
        this.modulus = this.polynomialFromCoeffs(modulus);
      }
    } else {
      // Sage: `modulus = R.irreducible_element(n)` (`finite_field_constructor.py:728`).
      this.modulus = this.irreducible_element(n);
    }

    if (this.modulus.degree() !== n) {
      throw new ValueError('the degree of the modulus does not equal the degree of the field');
    }
    // User-supplied polynomials need the constructor's irreducibility check.
    // Generated moduli have already been certified by their construction algorithm.
    if (modulus !== undefined && typeof modulus !== 'string' && !this.modulus.is_irreducible()) {
      throw new ValueError('finite field modulus must be irreducible but it is not');
    }
  }

  /**
   * Create polynomial from an array of the coefficients *below* the leading
   * one, constant term first; the monic leading coefficient is appended.
   *
   * This is the layout of our Conway polynomial table (which, like Sage's
   * `conway_polynomials` package, stores the `n` low coefficients only).
   */
  private polynomialFromCoeffs(coeffs: number[]): Polynomial<PrimeFieldElement> {
    const polyCoeffs = coeffs.map((c) => this.baseField.__call__(c));
    // Add leading coefficient 1 (monic polynomial)
    polyCoeffs.push(this.baseField.one());
    return new Polynomial(polyCoeffs, this.polynomialRing);
  }

  /**
   * Create a polynomial from a *complete* coefficient list, constant term
   * first (the layout PARI's `ffinit` and NTL's `GF2X` use).
   */
  private polynomialFromFullCoeffs(
    coeffs: ReadonlyArray<number | bigint>
  ): Polynomial<PrimeFieldElement> {
    return new Polynomial(
      coeffs.map((c) => this.baseField.__call__(c)),
      this.polynomialRing
    );
  }

  /**
   * Construct a monic irreducible polynomial of degree `n` over GF(p).
   *
   * Port of `PolynomialRing_dense_mod_p.irreducible_element`
   * (`reference/sage/src/sage/rings/polynomial/polynomial_ring.py:3560-3626`),
   * falling through to `PolynomialRing_dense_finite_field.irreducible_element`
   * (`polynomial_ring.py:2628-2681`) exactly where Sage does.
   *
   * Sage's default (`algorithm=None`) is fully deterministic:
   *
   * 1. `n == 1` -> `x - 1`;
   * 2. a Conway polynomial when one exists (`algorithm='conway'`);
   * 3. `p == 2` -> NTL's `GF2X_BuildSparseIrred` (`algorithm='minimal_weight'`);
   * 4. otherwise PARI's `ffinit` (`algorithm='adleman-lenstra'`).
   *
   * Steps 3 and 4 delegate to our ports of the very libraries Sage delegates
   * to: `@sagemath-ts/ntl-ts`'s `GF2X_BuildSparseIrred`
   * (`ntl/src/GF2XFactoring.cpp:900`) and `@sagemath-ts/parigp-ts`'s `ffinit`
   * (`pari/src/basemath/polarit3.c:3520`).
   *
   * @param n - degree of the polynomial to construct
   * @param algorithm - Sage's `algorithm` keyword, or `undefined` for the default
   */
  private irreducible_element(
    n: number,
    algorithm?: FiniteFieldModulusAlgorithm
  ): Polynomial<PrimeFieldElement> {
    const p = this.characteristic;

    if (n < 1) {
      throw new ValueError('degree must be at least 1');
    }

    let algo: FiniteFieldModulusAlgorithm | undefined = algorithm;

    if (algo === undefined) {
      if (n === 1) {
        // Sage: `return self((-1,1))`  # Polynomial x - 1
        return this.polynomialFromFullCoeffs([p - 1n, 1n]);
      }
      if (this.hasConwayPolynomial(n)) {
        algo = 'conway';
      } else if (p === 2n) {
        algo = 'minimal_weight';
      } else {
        algo = 'adleman-lenstra';
      }
    } else if (algo === 'primitive') {
      algo = this.hasConwayPolynomial(n) ? 'conway' : 'ffprimroot';
    }

    if (algo === 'adleman-lenstra') {
      // Sage: `return self(pari(p).ffinit(n))`
      return this.polynomialFromFullCoeffs(ffinit(p, n));
    }
    if (algo === 'conway') {
      if (!this.hasConwayPolynomial(n)) {
        // Sage raises RuntimeError from the conway_polynomials database.
        throw new ValueError(`Conway polynomial for GF(${p}^${n}) is not in the database`);
      }
      return this.polynomialFromCoeffs(conway_polynomial(Number(p), n));
    }
    if (algo === 'minimal_weight') {
      if (p !== 2n) {
        throw new NotImplementedError("'minimal_weight' option only implemented for p = 2");
      }
      // Sage: `return self(GF2X_BuildSparseIrred_list(n))`
      return this.polynomialFromFullCoeffs(gf2xToCoeffs(GF2X_BuildSparseIrred(n), n));
    }
    if (algo === 'ffprimroot') {
      // Sage: `self(pari(p).ffinit(n).ffgen().ffprimroot().charpoly())`.
      throw new NotImplementedError(
        "SAGE_NOT_IMPLEMENTED: irreducible_element(algorithm='ffprimroot') needs PARI's " +
          'ffgen/ffprimroot/charpoly, which are not in @sagemath-ts/parigp-ts'
      );
    }
    if (algo === 'first_lexicographic' && p === 2n) {
      // Sage: `return self(GF2X_BuildIrred_list(n))`
      return this.polynomialFromFullCoeffs(gf2xToCoeffs(GF2X_BuildIrred(n), n));
    }
    if (algo === 'random' && p === 2n) {
      // Sage: `return self(GF2X_BuildRandomIrred_list(n))`.  ntl-ts has no
      // `BuildRandomIrred` (it needs NTL's `IrredPolyMod`/`GF2XModulus`), so
      // fall through to the generic random search below, which is what Sage
      // itself does when the NTL import fails (`polynomial_ring.py:3615-3620`).
    }

    // "No suitable algorithm found, try algorithms from the base class."
    // (`polynomial_ring.py:3624-3625`) ->
    // PolynomialRing_dense_finite_field.irreducible_element
    if (algo === 'random') {
      return this.randomIrreducible(n);
    }
    if (algo === 'first_lexicographic') {
      return this.findIrreducible(n);
    }
    throw new ValueError(`no such algorithm for finding an irreducible polynomial: ${algo}`);
  }

  /**
   * Whether our Conway polynomial table has an entry for GF(p^n).
   *
   * Sage calls `exists_conway_polynomial(p, n)` (`polynomial_ring.py:3577`).
   *
   * @see Deviation: Finite Fields — Conway Table and Constructor Algorithms
   */
  private hasConwayPolynomial(n: number): boolean {
    const p = this.characteristic;
    if (p > BigInt(Number.MAX_SAFE_INTEGER)) {
      return false;
    }
    return has_conway_polynomial(Number(p), n);
  }

  /**
   * `algorithm='random'` of `PolynomialRing_dense_finite_field.irreducible_element`
   * (`polynomial_ring.py:2672-2676`):
   *
   *     while True:
   *         f = self.gen()**n + self.random_element(degree=(0, n - 1))
   *         if f.is_irreducible():
   *             return f
   */
  private randomIrreducible(n: number): Polynomial<PrimeFieldElement> {
    const x = this.polynomialRing.gen();
    const xPowN = x.pow(n) as Polynomial<PrimeFieldElement>;
    const randstate = current_randstate();
    const p = this.characteristic;
    for (;;) {
      const coeffs: PrimeFieldElement[] = [];
      for (let j = 0; j < n; j++) {
        coeffs.push(this.baseField.__call__(randstate.random_below(p)));
      }
      const candidate = xPowN.add(new Polynomial(coeffs, this.polynomialRing));
      if (this.isIrreducible(candidate)) {
        return candidate;
      }
    }
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
   * This is *not* the default any more: `irreducible_element` now delegates to
   * NTL / PARI exactly as Sage does.
   */
  private findIrreducible(n: number): Polynomial<PrimeFieldElement> {
    const x = this.polynomialRing.gen();
    const xPowN = x.pow(n) as Polynomial<PrimeFieldElement>;
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
    x?:
      | number
      | bigint
      | number[]
      | Polynomial<PrimeFieldElement>
      | FiniteFieldElement
      | PrimeFieldElement
      | unknown
  ): FiniteFieldElement {
    if (x instanceof FiniteFieldElement) {
      if (
        x.parent.characteristic !== this.characteristic ||
        x.parent.degree !== this.degree ||
        x.parent.variableName !== this.variableName ||
        !x.parent.modulus.eq(this.modulus)
      ) {
        throw new TypeError('no coercion defined');
      }
      return new FiniteFieldElement(x.lift, this);
    }

    if (x instanceof Polynomial) {
      if (String(x.parent.base_ring) === String(this)) {
        if (x.degree() > 0) throw new TypeError(`${x} is not a constant polynomial`);
        return this.__call__(x.getCoeff(0));
      }
      return new FiniteFieldElement(x, this);
    }

    if (typeof x === 'string') {
      return new FiniteFieldElement(this.polynomialRing.__call__(x), this);
    }

    if (
      x instanceof PrimeFieldElement ||
      x instanceof LegacyPrimeElement ||
      x instanceof IntegerMod ||
      x instanceof GF2Element
    ) {
      const modulus = x instanceof IntegerMod ? x.modulus : x.parent.characteristic;
      if (modulus % this.characteristic !== 0n) throw new TypeError('no coercion defined');
      const poly = this.polynomialRing.__call__(
        this.baseField.__call__(x instanceof GF2Element ? BigInt(x.value) : x.value)
      );
      return new FiniteFieldElement(poly, this);
    }

    if (typeof x === 'number' && !Number.isInteger(x)) {
      throw new TypeError('no coercion defined');
    }
    if (
      x == null ||
      x instanceof Integer ||
      x instanceof Rational ||
      typeof x === 'boolean' ||
      typeof x === 'number' ||
      typeof x === 'bigint'
    ) {
      const coeff = this.baseField.__call__(x);
      const poly = this.polynomialRing.__call__(coeff);
      return new FiniteFieldElement(poly, this);
    }

    if (Array.isArray(x)) {
      // Array of coefficients
      const coeffs = x.map((c) => {
        // A full coordinate vector in characteristic two uses vector_mod2_dense,
        // whose rational-denominator error differs from IntegerMod.__init__.
        // element_pari_ffelt.pyx:507; vector_mod2_dense.pyx:242-245.
        if (
          x.length === this.degree &&
          this.characteristic === 2n &&
          c instanceof Rational &&
          c.denominator % 2n === 0n
        ) {
          throw new ZeroDivisionError('inverse does not exist');
        }
        return this.baseField.__call__(c);
      });
      const poly = new Polynomial(coeffs, this.polynomialRing);
      return new FiniteFieldElement(poly, this);
    }

    throw new TypeError('no coercion defined');
  }

  /**
   * Create an element from its base-p integer representation.
   * finite_field_base.pyx:477-482 rejects values outside [0, order).
   */
  fromInteger(n: bigint): FiniteFieldElement {
    if (n < 0n || n >= this.order) {
      throw new ValueError('n must be between 0 and self.order()');
    }
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
  gen(n: unknown = 0): FiniteFieldElement {
    checkFiniteGeneratorIndex(n, true);
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
    // finite_field_base.pyx:1058-1059: a single Python randrange draw.
    // Sage 10.3 still used its vector-space randomizer here; the bundled
    // reference has since changed the no-argument path.
    return this.fromInteger(current_randstate().python_random().randrange(this.order));
  }

  is_field(): boolean {
    return true;
  }

  toString(): string {
    return `Finite Field in ${this.variableName} of size ${this.characteristic}^${this.degree}`;
  }
}

/** Options corresponding to Sage's name and modulus keyword arguments. */
export interface FiniteFieldOptions {
  name?: string;
  /** Full coefficient list, constant term first, including the leading term. */
  modulus?:
    | readonly (IntegerLike | number)[]
    | Polynomial<PrimeFieldElement>
    | FiniteFieldModulusAlgorithm
    | null;
}

/** Adapt certify_names/normalize_names(1, name), category_object.pyx:1010-1093. */
function finiteFieldName(name: string, extension: boolean): string {
  // Python str.strip includes C0 separators and NEL, but does not strip BOM.
  const pythonSpace =
    // biome-ignore lint/suspicious/noControlCharactersInRegex: Python strips these exact control characters.
    /^[\u0009-\u000d\u001c-\u0020\u0085\u00a0\u1680\u2000-\u200a\u2028\u2029\u202f\u205f\u3000]+|[\u0009-\u000d\u001c-\u0020\u0085\u00a0\u1680\u2000-\u200a\u2028\u2029\u202f\u205f\u3000]+$/g;
  const names =
    extension && name.includes(',')
      ? name.split(',').map((s) => s.replace(pythonSpace, ''))
      : [name];
  const seen = new Set<string>();
  for (const n of names) {
    const quote = n.includes("'") && !n.includes('"') ? '"' : "'";
    let repr = quote;
    for (const ch of n) {
      const code = ch.codePointAt(0)!;
      if (ch === quote || ch === '\\') repr += '\\' + ch;
      else if (ch === '\n') repr += '\\n';
      else if (ch === '\r') repr += '\\r';
      else if (ch === '\t') repr += '\\t';
      else if (ch !== ' ' && /[\p{C}\p{Z}]/u.test(ch)) {
        repr +=
          code <= 255
            ? `\\x${code.toString(16).padStart(2, '0')}`
            : code <= 65535
              ? `\\u${code.toString(16).padStart(4, '0')}`
              : `\\U${code.toString(16).padStart(8, '0')}`;
      } else repr += ch;
    }
    repr += quote;
    if (!n) throw new ValueError('variable name must be nonempty');
    if (!/^[\p{L}\p{N}_]+$/u.test(n))
      throw new ValueError(`variable name ${repr} is not alphanumeric`);
    if (!/^\p{L}/u.test(n))
      throw new ValueError(`variable name ${repr} does not start with a letter`);
    if (seen.has(n)) throw new ValueError(`variable name ${repr} appears more than once`);
    seen.add(n);
  }
  if (names.length !== 1)
    throw new IndexError('the number of names must equal the number of generators');
  return names[0]!;
}

/**
 * Construct a finite field GF(q) where q = p^n, supporting extension fields.
 *
 * This is the full implementation that supports both prime fields and extension fields.
 * For prime-only fields, use the simpler GF() from finite_field_constructor.ts.
 *
 * @param q - The order of the field (must be a prime power)
 * @param nameOrOptions - Generator name or name/modulus options (default: 'a')
 * @param options - Options when a generator name is supplied separately
 * @returns The finite field GF(q)
 * @see Deviation: Finite Field Coercion and Backend Boundaries
 *
 * @example
 * ```typescript
 * const F4 = GFExtended(4);           // GF(2^2)
 * const F9 = GFExtended(9);           // GF(3^2)
 * const F8 = GFExtended(8, 'b');      // GF(2^3) with generator named 'b'
 * ```
 */
export function GFExtended(
  q: IntegerLike | number,
  nameOrOptions: string | FiniteFieldOptions = 'a',
  options?: FiniteFieldOptions
): PrimeField | FiniteFieldExtension {
  const opts = typeof nameOrOptions === 'string' ? (options ?? {}) : nameOrOptions;
  // Never silently turn keyword arguments into a generator name.
  for (const key of Object.keys(opts)) {
    if (key !== 'name' && key !== 'modulus') {
      throw new TypeError(
        `create_key_and_extra_args() got an unexpected keyword argument '${key}'`
      );
    }
  }
  const order = ZZ.__call__(q);
  if (order < 2n) throw new ValueError('the order of a finite field must be at least 2');

  // integer.pyx:5030-5054 delegates perfect-power decomposition to PARI.
  // Decompose before validating the name, then check primality as the factory
  // does; even an invalid field order can have an invalid name first.
  const [exponent, p] = Z_isanypower(order);
  const n = exponent || 1;
  const name = finiteFieldName(
    opts.name ?? (typeof nameOrOptions === 'string' ? nameOrOptions : 'a'),
    n > 1
  );
  if (!is_prime(p)) throw new ValueError('the order of a finite field must be a prime power');

  let modulus: Polynomial<PrimeFieldElement> | FiniteFieldModulusAlgorithm | undefined;
  if (opts.modulus != null) {
    if (typeof opts.modulus === 'string' || opts.modulus instanceof Polynomial)
      modulus = opts.modulus;
    else {
      const base = new PrimeField(p);
      const ring = new PolynomialRing(base, 'x');
      modulus = new Polynomial<PrimeFieldElement>(
        opts.modulus.map((c) => base.__call__(c)),
        ring
      );
    }
  }
  if (n === 1) return new PrimeField(p, { modulus });
  return new FiniteFieldExtension(p, n, modulus, name);
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
 * @param modulus - Optional modulus polynomial coefficients (constant term
 *   first, without the monic leading 1), or one of SageMath's `algorithm`
 *   strings (`'conway'`, `'minimal_weight'`, `'adleman-lenstra'`, ...)
 * @param variableName - Name for the generator
 */
export function GFpn(
  p: number | bigint,
  n: number,
  modulus?: number[] | FiniteFieldModulusAlgorithm,
  variableName: string = 'a'
): FiniteFieldExtension {
  return new FiniteFieldExtension(p, n, modulus, variableName);
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
