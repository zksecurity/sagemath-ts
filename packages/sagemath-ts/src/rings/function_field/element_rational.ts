/**
 * @module sage/rings/function_field/element_rational
 * @description Elements of function fields: rational
 *
 * Port of: sage/rings/function_field/element_rational.pyx
 *
 * The underlying object is an element of `Frac(k[x])`.  This port has no
 * fraction-field type, so the numerator/denominator pair is stored directly on
 * the element and normalised (coprime, monic denominator) exactly the way
 * SageMath's `FractionFieldElement` normalises it.
 */

import { NotImplementedError, ValueError, ZeroDivisionError } from '../../errors.js';
import { Polynomial } from '../polynomial/polynomial_element.js';
import { compare_constants, divide_constants } from './constant_field.js';
import type { ConstantFieldElement } from './constant_field.js';
import { FunctionFieldElement } from './element.js';
import type { RationalFunctionField } from './function_field_rational.js';
import type { FunctionFieldIdeal } from './ideal.js';
import { FunctionFieldPlace } from './place.js';

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

/** Multiplicity of the irreducible ``p`` in ``f``. */
function polynomial_valuation<C extends ConstantFieldElement>(
  f: Polynomial<C>,
  p: Polynomial<C>
): bigint {
  if (f.isZero()) {
    throw new ValueError('valuation of the zero polynomial is not defined');
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
 * SageMath's `_is_atomic` for polynomials: a polynomial prints without needing
 * parentheses exactly when it is a single monomial with an atomic coefficient.
 *
 * @see Reference: sage/rings/polynomial/polynomial_element.pyx:3415 (_is_atomic)
 */
function polynomial_is_atomic<C extends ConstantFieldElement>(f: Polynomial<C>): boolean {
  if (f.isZero()) {
    return true;
  }
  // degree == valuation, i.e. exactly one nonzero coefficient
  let nonzero = 0;
  for (let i = 0; i <= f.degree(); i++) {
    if (!f.getCoeff(i).isZero()) {
      nonzero++;
    }
  }
  if (nonzero !== 1) {
    return false;
  }
  const s = f.leading_coefficient().toString();
  return !s.includes('+') && !s.includes('-') && !s.includes(' ');
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
  readonly _num: Polynomial<C>;
  readonly _den: Polynomial<C>;

  constructor(
    parent: RationalFunctionField<C>,
    num: Polynomial<C>,
    den?: Polynomial<C>,
    reduce: boolean = true
  ) {
    super();
    this._parent = parent;
    const R = parent._ring;
    let n = num;
    let d = den ?? R.one();

    if (d.isZero()) {
      throw new ZeroDivisionError('fraction field element division by zero');
    }

    if (reduce) {
      if (n.isZero()) {
        d = R.one();
      } else {
        const g = n.gcd(d);
        if (g.degree() > 0 || !g.leading_coefficient().eq(1)) {
          n = n.quo_rem(g)[0];
          d = d.quo_rem(g)[0];
        }
        // SageMath normalises the denominator of Frac(k[x]) to be monic.
        const lc = d.leading_coefficient();
        if (!lc.eq(1)) {
          const inv = lc.inv();
          n = n.scalar_mul(inv);
          d = d.scalar_mul(inv);
        }
      }
    }
    this._num = n;
    this._den = d;
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
   * SageMath returns the `Frac(k[x])` element; this port has no fraction-field
   * type, so we return the element itself.
   *
   * @see Reference: sage/rings/function_field/element_rational.pyx:67 (element)
   * @see Deviation: function-field rational elements carry their own fraction
   */
  element(): FunctionFieldElement_rational<C> {
    return this;
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
    return this._num.isZero();
  }

  override is_one(): boolean {
    return this._den.eq(this._num) && !this._num.isZero();
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
    if (this.is_zero()) {
      return '0';
    }
    let s = this._num.toString();
    if (!(this._den.degree() === 0 && this._den.getCoeff(0).eq(1))) {
      const denom_string = this._den.toString();
      const numStr = polynomial_is_atomic(this._num) ? s : `(${s})`;
      if (
        polynomial_is_atomic(this._den) &&
        !denom_string.includes('*') &&
        !denom_string.includes('/')
      ) {
        s = `${numStr}/${denom_string}`;
      } else {
        s = `${numStr}/(${denom_string})`;
      }
    }
    return s;
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
    return this._new(this._num.mul(o._den).add(o._num.mul(this._den)), this._den.mul(o._den));
  }

  /**
   * @see Reference: sage/rings/function_field/element_rational.pyx:195 (_sub_)
   */
  override sub(other: FunctionFieldElement<C>): FunctionFieldElement_rational<C> {
    const o = other as FunctionFieldElement_rational<C>;
    return this._new(this._num.mul(o._den).sub(o._num.mul(this._den)), this._den.mul(o._den));
  }

  /**
   * @see Reference: sage/rings/function_field/element_rational.pyx:213 (_mul_)
   */
  override mul(other: FunctionFieldElement<C>): FunctionFieldElement_rational<C> {
    const o = other as FunctionFieldElement_rational<C>;
    return this._new(this._num.mul(o._num), this._den.mul(o._den));
  }

  /**
   * @see Reference: sage/rings/function_field/element_rational.pyx:231 (_div_)
   */
  override div(other: FunctionFieldElement<C>): FunctionFieldElement_rational<C> {
    const o = other as FunctionFieldElement_rational<C>;
    if (o.is_zero()) {
      throw new ZeroDivisionError('division by zero');
    }
    return this._new(this._num.mul(o._den), this._den.mul(o._num));
  }

  override neg(): FunctionFieldElement_rational<C> {
    return this._new(this._num.neg(), this._den);
  }

  override inv(): FunctionFieldElement_rational<C> {
    if (this.is_zero()) {
      throw new ZeroDivisionError('cannot invert zero');
    }
    return this._new(this._den, this._num);
  }

  override pow(n: bigint | number): FunctionFieldElement_rational<C> {
    let e = BigInt(n);
    // eslint-disable-next-line @typescript-eslint/no-this-alias
    let base: FunctionFieldElement_rational<C> = this;
    if (e < 0n) {
      base = base.inv();
      e = -e;
    }
    let result = this._parent.one();
    while (e > 0n) {
      if (e & 1n) {
        result = result.mul(base);
      }
      base = base.mul(base);
      e >>= 1n;
    }
    return result;
  }

  override scalar_mul(c: C): FunctionFieldElement_rational<C> {
    return this._new(this._num.scalar_mul(c), this._den);
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
  matrix(): Array<Array<FunctionFieldElement_rational<C>>> {
    return [[this]];
  }

  /**
   * @see Reference: sage/rings/function_field/element.pyx:525 (trace)
   */
  trace(): FunctionFieldElement_rational<C> {
    return this;
  }

  /**
   * @see Reference: sage/rings/function_field/element.pyx:538 (norm)
   */
  norm(): FunctionFieldElement_rational<C> {
    return this;
  }

  /**
   * Return the valuation of the rational function at the place.
   *
   * ``place`` may be a place of the function field or an irreducible
   * polynomial (given as a polynomial or as an element of the field).
   *
   * @see Reference: sage/rings/function_field/element_rational.pyx:278 (valuation)
   */
  override valuation(
    place: FunctionFieldPlace<C> | Polynomial<C> | FunctionFieldElement_rational<C>
  ): bigint | number {
    if (!(place instanceof FunctionFieldPlace)) {
      // place is an irreducible polynomial
      const p =
        place instanceof Polynomial ? place : (place as FunctionFieldElement_rational<C>)._num;
      if (this.is_zero()) {
        return Number.POSITIVE_INFINITY;
      }
      return polynomial_valuation(this._num, p) - polynomial_valuation(this._den, p);
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
    return polynomial_is_square(this._num) && polynomial_is_square(this._den);
  }

  /**
   * Return the square root of the rational function.
   *
   * @see Reference: sage/rings/function_field/element_rational.pyx:338 (sqrt)
   */
  sqrt(
    all: boolean = false
  ): FunctionFieldElement_rational<C> | Array<FunctionFieldElement_rational<C>> {
    if (this.is_zero()) {
      return all ? [this] : this;
    }
    if (!this.is_square()) {
      throw new ValueError('element is not a square');
    }
    const r = this._new(polynomial_sqrt(this._num), polynomial_sqrt(this._den));
    if (all) {
      return r.is_zero() ? [r] : [r, r.neg()];
    }
    return r;
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
    // SageMath does not special-case zero here: `K(0).factor()` returns the
    // empty factorization with unit 0 (verified against SageMath 10.3).
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
      return compare_polynomials(a[0], b[0]);
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
      throw new ValueError('ideal must be principal');
    }
    const f = gens[0] as FunctionFieldElement_rational<C>;
    if (f._den.degree() !== 0) {
      throw new ValueError('ideal generator must be integral');
    }
    if (this._den.degree() !== 0) {
      throw new ValueError('element must be integral');
    }
    const [g, s] = this._num.xgcd(f._num);
    if (g.degree() !== 0) {
      // SageMath's `Polynomial.inverse_mod`
      // (`reference/sage/src/sage/rings/polynomial/polynomial_element.pyx:1644`)
      throw new ValueError('Impossible inverse modulo');
    }
    // s * self = g mod f, so s/g is the inverse; g is monic (degree 0 -> 1).
    return this._new(s.mod(f._num));
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
      return this.inv().is_nth_power(-e);
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

/** Whether ``f`` is a square in `k[x]`. */
function polynomial_is_square<C extends ConstantFieldElement>(f: Polynomial<C>): boolean {
  if (f.isZero()) {
    return true;
  }
  for (const [, e] of monic_irreducible_factors(f)) {
    if (e % 2n !== 0n) {
      return false;
    }
  }
  return constant_is_square(f.leading_coefficient());
}

/** A square root of ``f`` in `k[x]`; ``f`` must be a square. */
function polynomial_sqrt<C extends ConstantFieldElement>(f: Polynomial<C>): Polynomial<C> {
  const R = f.parent;
  if (f.isZero()) {
    return R.zero();
  }
  let r = R.one();
  for (const [g, e] of monic_irreducible_factors(f)) {
    r = r.mul(g.pow(Number(e / 2n)));
  }
  return r.scalar_mul(constant_sqrt(f.leading_coefficient()));
}

function constant_is_square<C extends ConstantFieldElement>(c: C): boolean {
  const o = c as unknown as { is_square?: () => boolean };
  if (typeof o.is_square === 'function') {
    return o.is_square();
  }
  throw new NotImplementedError(
    'SAGE_NOT_IMPLEMENTED: is_square for this constant field (no is_square() on its elements)'
  );
}

function constant_sqrt<C extends ConstantFieldElement>(c: C): C {
  const o = c as unknown as { sqrt?: () => C };
  if (typeof o.sqrt === 'function') {
    return o.sqrt();
  }
  throw new NotImplementedError(
    'SAGE_NOT_IMPLEMENTED: sqrt for this constant field (no sqrt() on its elements)'
  );
}
