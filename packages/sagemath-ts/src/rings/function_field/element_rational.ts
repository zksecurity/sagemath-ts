/**
 * @module sage/rings/function_field/element_rational
 * @description Elements of function fields: rational
 *
 * Port of: sage/rings/function_field/element_rational.pyx
 *
 * Elements wrap the selected `Frac(k[x])` element, preserving its identity,
 * normalization state and arithmetic dispatch.
 */

import { ArithmeticError, AssertionError, NotImplementedError, TypeError, ValueError } from '../../errors.js';
import { isFractionElement, type FractionElement } from '../fraction_field_element.js';
import { Polynomial } from '../polynomial/polynomial_element.js';
import { compare_constants, divide_constants, constant_field_uses_fpt } from './constant_field.js';
import type { ConstantFieldElement } from './constant_field.js';
import { FunctionFieldElement } from './element.js';
import type { RationalFunctionField } from './function_field_rational.js';
import type { FunctionFieldIdeal } from './ideal.js';
import { FunctionFieldPlace } from './place.js';

// RationalFunctionField.free_module caches its map=True result by that flag
// alone. The one-dimensional matrix specialization retains this parent-level
// initialization state, independently of its fresh immutable matrix arrays.
const initializedMatrixSpaces = new WeakSet<object>();

/**
 * Compare two polynomials the way SageMath's ``Polynomial._richcmp_`` does:
 * first by degree (treating the zero polynomial as smaller than any nonzero
 * constant only through the constant comparison), then in dictionary order
 * starting with the coefficient of largest degree.
 *
 * @see Reference: sage/rings/polynomial/polynomial_element.pyx:960 (_richcmp_)
 */
export function compare_polynomials<C extends ConstantFieldElement>(
  a: Polynomial<C>,
  b: Polynomial<C>
): number {
  const d1 = a.degree();
  const d2 = b.degree();
  const zero = a.parent.base_ring.zero() as C;

  if (d1 === -1) {
    if (d2 === -1) {
      return 0;
    }
    if (d2 === 0) {
      return compare_constants(zero, b.getCoeff(0));
    }
    return -1;
  }
  if (d1 === 0) {
    if (d2 === -1) {
      return compare_constants(a.getCoeff(0), zero);
    }
    if (d2 === 0) {
      return compare_constants(a.getCoeff(0), b.getCoeff(0));
    }
    return -1;
  }
  if (d1 !== d2) {
    return d1 < d2 ? -1 : 1;
  }
  for (let i = d1; i >= 0; i--) {
    const c = compare_constants(a.getCoeff(i), b.getCoeff(i));
    if (c !== 0) {
      return c;
    }
  }
  return 0;
}

/** FpT compares degree then coefficients from the constant term upwards. */
function compare_fpt_polynomials<C extends ConstantFieldElement>(a: Polynomial<C>, b: Polynomial<C>): number {
  if (a.degree() !== b.degree()) return a.degree() < b.degree() ? -1 : 1;
  for (let i = 0; i <= a.degree(); i++) {
    const c = compare_constants(a.getCoeff(i), b.getCoeff(i));
    if (c) return c;
  }
  return 0;
}

/** Multiplicity of the irreducible ``p`` in ``f``. */
function polynomial_valuation<C extends ConstantFieldElement>(
  f: Polynomial<C>,
  p: Polynomial<C>
): bigint | number {
  if (f.isZero()) return Number.POSITIVE_INFINITY;
  if (p.degree() === 0) {
    throw new ArithmeticError('The polynomial, p, must be non-constant.');
  }
  let v = 0n;
  let cur = f;
  for (;;) {
    const [q, r] = cur.quo_rem(p);
    if (!r.isZero()) {
      return v;
    }
    cur = q;
    v += 1n;
  }
}

/**
 * Elements of a rational function field.
 *
 * @see Reference: sage/rings/function_field/element_rational.pyx:29 (FunctionFieldElement_rational)
 */
export class FunctionFieldElement_rational<
  C extends ConstantFieldElement,
> extends FunctionFieldElement<C> {
  readonly _parent: RationalFunctionField<C>;
  readonly _x: FractionElement<C>;
  get _num(): Polynomial<C> { return this._x.numerator(); }
  get _den(): Polynomial<C> { return this._x.denominator(); }

  constructor(
    parent: RationalFunctionField<C>,
    num: Polynomial<C> | FractionElement<C>,
    den?: Polynomial<C>,
    reduce: boolean = true
  ) {
    super();
    this._parent = parent;
    const F = parent.field();
    this._x = isFractionElement(num) && den === undefined ? num
      : new F._element_class(F, num, den ?? parent._ring.one(), { coerce: false, reduce });
  }

  override get parent(): RationalFunctionField<C> {
    return this._parent;
  }

  private _new(num: Polynomial<C>, den?: Polynomial<C>): FunctionFieldElement_rational<C> {
    return new FunctionFieldElement_rational(this._parent, num, den);
  }

  /**
   * Return the underlying representation of the element.
   *
   * @see Reference: sage/rings/function_field/element_rational.pyx:67 (element)
   */
  element(): FractionElement<C> {
    return this._x;
  }

  /**
   * Return a list with just the element.
   *
   * @see Reference: sage/rings/function_field/element_rational.pyx:88 (list)
   */
  list(): Array<FunctionFieldElement_rational<C>> {
    return [this];
  }

  /**
   * Return the numerator of the rational function.
   *
   * @see Reference: sage/rings/function_field/element_rational.pyx:250 (numerator)
   */
  numerator(): Polynomial<C> {
    return this._num;
  }

  /**
   * Return the denominator of the rational function.
   *
   * @see Reference: sage/rings/function_field/element_rational.pyx:264 (denominator)
   */
  denominator(): Polynomial<C> {
    return this._den;
  }

  override is_zero(): boolean {
    return this._x.isZero();
  }

  override is_one(): boolean {
    return this.eq(this._parent.one());
  }

  /**
   * Return the string representation of the element.
   *
   * Mirrors SageMath's `FractionFieldElement._repr_`
   * (`reference/sage/src/sage/rings/fraction_field_element.pyx:523`).
   *
   * @see Reference: sage/rings/function_field/element_rational.pyx:103 (_repr_)
   */
  _repr_(): string {
    return this._x.toString();
  }

  override toString(): string {
    return this._repr_();
  }

  /**
   * Compare the element with ``other``.
   *
   * @see Reference: sage/rings/function_field/element_rational.pyx:146 (_richcmp_)
   * @see Reference: sage/rings/fraction_field_element.pyx:994 (_richcmp_)
   */
  cmp(other: FunctionFieldElement_rational<C>): number {
    if (constant_field_uses_fpt(this._parent.constant_base_field())) {
      return compare_fpt_polynomials(this._num, other._num) || compare_fpt_polynomials(this._den, other._den);
    }
    return compare_polynomials(this._num.mul(other._den), this._den.mul(other._num));
  }

  override eq(other: FunctionFieldElement<C>): boolean {
    const o = other as FunctionFieldElement_rational<C>;
    return this.cmp(o) === 0;
  }

  /**
   * @see Reference: sage/rings/function_field/element_rational.pyx:177 (_add_)
   */
  override add(other: FunctionFieldElement<C>): FunctionFieldElement_rational<C> {
    const o = other as FunctionFieldElement_rational<C>;
    return new FunctionFieldElement_rational(this._parent, this._x.add(o._x));
  }

  /**
   * @see Reference: sage/rings/function_field/element_rational.pyx:195 (_sub_)
   */
  override sub(other: FunctionFieldElement<C>): FunctionFieldElement_rational<C> {
    const o = other as FunctionFieldElement_rational<C>;
    return new FunctionFieldElement_rational(this._parent, this._x.sub(o._x));
  }

  /**
   * @see Reference: sage/rings/function_field/element_rational.pyx:213 (_mul_)
   */
  override mul(other: FunctionFieldElement<C>): FunctionFieldElement_rational<C> {
    const o = other as FunctionFieldElement_rational<C>;
    return new FunctionFieldElement_rational(this._parent, this._x.mul(o._x));
  }

  /**
   * @see Reference: sage/rings/function_field/element_rational.pyx:231 (_div_)
   */
  override div(other: FunctionFieldElement<C>): FunctionFieldElement_rational<C> {
    const o = other as FunctionFieldElement_rational<C>;
    return new FunctionFieldElement_rational(this._parent, this._x.div(o._x));
  }

  override neg(): FunctionFieldElement_rational<C> {
    // structure/element.pyx:_neg_ uses scalar multiplication by -1.
    return this.mul(this._parent.__call__(-1n));
  }

  override inv(): FunctionFieldElement_rational<C> {
    return this._parent.one().div(this);
  }

  override pow(n: bigint | number): FunctionFieldElement_rational<C> {
    let e = BigInt(n);
    // eslint-disable-next-line @typescript-eslint/no-this-alias
    let base: FunctionFieldElement_rational<C> = this;
    if (e < 0n) {
      base = base.inv();
      e = -e;
    }
    // arith/power.pyx:generic_power_pos starts at the least set bit.
    // In particular f^1 returns f without normalizing an unreduced fraction.
    if (e === 0n) return this._parent.one();
    while (!(e & 1n)) {
      base = base.mul(base);
      e >>= 1n;
    }
    let result = base;
    e >>= 1n;
    while (e > 0n) {
      base = base.mul(base);
      if (e & 1n) result = base.mul(result);
      e >>= 1n;
    }
    return result;
  }

  override scalar_mul(c: C): FunctionFieldElement_rational<C> {
    return this.mul(this._parent.__call__(c));
  }

  /**
   * Return the max degree between the denominator and numerator.
   *
   * @see Reference: sage/rings/function_field/element.pyx:561 (degree)
   */
  degree(): bigint {
    return BigInt(Math.max(this._den.degree(), this._num.degree()));
  }

  /**
   * Return the matrix of multiplication by this element, over the base field.
   *
   * For a rational function field the base field is the field itself, so the
   * matrix is the `1 x 1` matrix `[self]`.
   *
   * @see Reference: sage/rings/function_field/element.pyx:444 (matrix)
   */
  matrix(base?: unknown): Array<Array<FunctionFieldElement_rational<C>>> {
    // cached_method hashes arguments before entering matrix or free_module.
    const unhashable = Array.isArray(base) ? 'list' : base instanceof Set ? 'set'
      : base instanceof Map || (base !== null && typeof base === 'object'
        && (Object.getPrototypeOf(base) === Object.prototype || Object.getPrototypeOf(base) === null)) ? 'dict' : null;
    if (unhashable !== null) throw new TypeError(`unhashable type: '${unhashable}'`);
    if (!initializedMatrixSpaces.has(this.parent)) {
      if (base !== undefined && base !== null && base !== this.parent)
        throw new ValueError('base must be the rational function field itself');
      initializedMatrixSpaces.add(this.parent);
    }
    // The native basis vector maps to one; multiplying by it still runs the
    // fraction backend, which can normalize an unreduced input.
    const entry = this.mul(this.parent.one());
    return immutable_matrix_array([immutable_matrix_array([entry])]);
  }

  /**
   * @see Reference: sage/rings/function_field/element.pyx:525 (trace)
   */
  trace(): FunctionFieldElement_rational<C> {
    // matrix2.pyx:trace starts its diagonal sum at the base-ring zero.
    return this.parent.zero().add(this.matrix()[0]![0]!);
  }

  /**
   * @see Reference: sage/rings/function_field/element.pyx:538 (norm)
   */
  norm(): FunctionFieldElement_rational<C> {
    // matrix2.pyx:determinant returns the sole entry for a 1-by-1 matrix.
    return this.matrix()[0]![0]!;
  }

  /**
   * Return the valuation of the rational function at the place.
   *
   * ``place`` may be a place of the function field or a nonconstant
   * polynomial (given as a polynomial or as an integral element of the field).
   *
   * @see Reference: sage/rings/function_field/element_rational.pyx:278 (valuation)
   */
  override valuation(
    place: FunctionFieldPlace<C> | Polynomial<C> | FunctionFieldElement_rational<C>
  ): bigint | number {
    if (!(place instanceof FunctionFieldPlace)) {
      // Sage first converts the entire argument to the polynomial ring. A
      // rational-function denominator must be a unit; numerator-only conversion
      // silently changes the requested valuation.
      const R = this._parent._ring;
      const p = place instanceof Polynomial ? R.__call__(place)
        : R.__call__(R.fraction_field().__call__(place.numerator(), place.denominator()));
      const numerator = polynomial_valuation(this._num, p);
      const denominator = polynomial_valuation(this._den, p);
      // Evaluate both valuations: even zero must validate its denominator at p.
      return numerator === Number.POSITIVE_INFINITY ? numerator
        : (numerator as bigint) - (denominator as bigint);
    }
    const prime = place.prime_ideal();
    const ideal = prime.ring().ideal(this);
    return prime.valuation(ideal);
  }

  /**
   * Return whether the element is a square.
   *
   * @see Reference: sage/rings/function_field/element_rational.pyx:316 (is_square)
   */
  is_square(): boolean {
    return this._x.is_square();
  }

  /**
   * Return the square root of the rational function.
   *
   * @see Reference: sage/rings/function_field/element_rational.pyx:338 (sqrt)
   */
  sqrt(
    all: boolean = false
  ): FunctionFieldElement_rational<C> | Array<FunctionFieldElement_rational<C>> {
    if (all) return (this._x.sqrt(true, true) as FractionElement<C>[]).map(r => this._parent.__call__(r));
    return this._parent.__call__(this._x.sqrt() as FractionElement<C>);
  }

  /**
   * Factor the rational function.
   *
   * SageMath returns a `Factorization`; this port returns the unit together
   * with the list of (monic irreducible, exponent) pairs, exponents being
   * negative for the factors of the denominator.
   *
   * @see Reference: sage/rings/function_field/element_rational.pyx:471 (factor)
   * @see Deviation: function-field factorizations returned as plain data
   */
  factor(): {
    unit: FunctionFieldElement_rational<C>;
    factors: Array<[FunctionFieldElement_rational<C>, bigint]>;
  } {
    // The bundled polynomial backends reject zero. Do not bypass their guard
    // when filtering the constant unit from the factorization.
    if (this.is_zero()) this._num.factor();
    const R = this._parent._ring;
    const raw: Array<[Polynomial<C>, bigint]> = [];
    let unit = this._num.leading_coefficient();
    for (const [f, e] of monic_irreducible_factors(this._num)) {
      raw.push([f, e]);
    }
    for (const [f, e] of monic_irreducible_factors(this._den)) {
      raw.push([f, -e]);
    }
    // `Factorization.sort()` orders by (degree, exponent, prime) whenever the
    // primes have a `degree` method.
    // @see Reference: sage/structure/factorization.py:671 (sort)
    raw.sort((a, b) => {
      if (a[0].degree() !== b[0].degree()) {
        return a[0].degree() - b[0].degree();
      }
      if (a[1] !== b[1]) {
        return a[1] < b[1] ? -1 : 1;
      }
      return constant_field_uses_fpt(this._parent.constant_base_field())
        ? compare_fpt_polynomials(a[0], b[0]) : compare_polynomials(a[0], b[0]);
    });
    const factors: Array<[FunctionFieldElement_rational<C>, bigint]> = raw.map(([f, e]) => [
      this._new(f),
      e,
    ]);
    unit = divide_constants(unit, this._den.leading_coefficient());
    return { unit: this._new(R.__call__(unit)), factors };
  }

  /**
   * Return an inverse of the element modulo the integral ideal `I`.
   *
   * @see Reference: sage/rings/function_field/element_rational.pyx:494 (inverse_mod)
   */
  inverse_mod(I: FunctionFieldIdeal<C>): FunctionFieldElement_rational<C> {
    const gens = I.gens();
    if (gens.length !== 1) {
      throw new AssertionError();
    }
    const f = gens[0] as FunctionFieldElement_rational<C>;
    if (!f._den.eq(1n)) {
      throw new AssertionError();
    }
    if (!this._den.eq(1n)) {
      throw new AssertionError();
    }
    const [g, s] = this._num.xgcd(f._num);
    if (g.degree() !== 0) {
      // SageMath's `Polynomial.inverse_mod`
      // (`reference/sage/src/sage/rings/polynomial/polynomial_element.pyx:1644`)
      throw new ValueError('Impossible inverse modulo');
    }
    // s * self = g mod f. A constant gcd can be nonmonic on the
    // native zero-modulus path, so divide the Bezout coefficient by that unit.
    // Polynomial.inverse_mod returns the Bezout coefficient directly, including
    // inversion of a constant modulo the zero polynomial.
    return this._new(s.scalar_mul(g.leading_coefficient().inv()));
  }

  /**
   * Return whether this element is an ``n``-th power.
   *
   * @see Reference: sage/rings/function_field/element_rational.pyx:362 (is_nth_power)
   */
  override is_nth_power(n: bigint | number): boolean {
    const e = BigInt(n);
    if (e === 1n) {
      return true;
    }
    if (e < 0n) {
      // The source cpdef is noexcept: inversion failure produces false.
      return this.is_zero() ? false : this.inv().is_nth_power(-e);
    }
    const p = this._parent.characteristic();
    if (e === 2n) {
      return this.is_square();
    }
    if (p !== 0n && e === p) {
      throw new NotImplementedError(
        'SAGE_NOT_IMPLEMENTED: FunctionFieldElement_rational.is_nth_power for n equal to the ' +
          'characteristic (needs sage/rings/function_field/derivations_rational.py)'
      );
    }
    throw new NotImplementedError('is_nth_power() not implemented for the given n');
  }

  /**
   * Return an ``n``-th root of this element.
   *
   * @see Reference: sage/rings/function_field/element_rational.pyx:415 (nth_root)
   */
  override nth_root(n: bigint | number): FunctionFieldElement_rational<C> {
    const e = BigInt(n);
    if (e === 0n) {
      if (!this.is_one()) {
        throw new ValueError('element is not a 0-th power');
      }
      return this;
    }
    if (e === 1n) {
      return this;
    }
    if (e < 0n) {
      return this.inv().nth_root(-e);
    }
    if (e === 2n) {
      if (this._parent.characteristic() === 2n && !this.is_square()) {
        throw new ValueError('element is not an n-th power');
      }
      return this.sqrt() as FunctionFieldElement_rational<C>;
    }
    throw new NotImplementedError(`nth_root() not implemented for ${e}`);
  }
}

/** Return the monic irreducible factors of ``f`` with their multiplicities. */
function monic_irreducible_factors<C extends ConstantFieldElement>(
  f: Polynomial<C>
): Array<[Polynomial<C>, bigint]> {
  if (f.degree() <= 0) {
    return [];
  }
  return f
    .factor()
    .filter(([g]) => g.degree() > 0)
    .map(([g, e]) => [g, BigInt(e)] as [Polynomial<C>, bigint]);
}

/** Nested-array adapter for Sage's immutable matrix entry storage. */
function immutable_matrix_array<T>(entries: T[]): T[] {
  const target = Object.freeze(entries);
  return new Proxy(target, {
    set() {
      throw new ValueError('matrix is immutable; please change a copy instead (i.e., use copy(M) to change a copy of M).');
    },
  }) as T[];
}
