/**
 * @module sage/rings/function_field/function_field_rational
 * @description Function Fields: rational
 *
 * Port of: sage/rings/function_field/function_field_rational.py
 */

import { NotImplementedError, ValueError } from '../../errors.js';
import { Polynomial } from '../polynomial/polynomial_element.js';
import { PolynomialRing } from '../polynomial/polynomial_ring.js';
import {
  constant_field_cardinality,
  constant_field_element_list,
  constant_field_is_finite,
} from './constant_field.js';
import type { ConstantField, ConstantFieldElement } from './constant_field.js';
import type { FunctionFieldDivisor } from './divisor.js';
import { FunctionFieldElement_rational } from './element_rational.js';
import { FunctionField } from './function_field.js';
import type { FunctionFieldIdeal } from './ideal.js';
import {
  FunctionFieldMaximalOrderInfinite_rational,
  FunctionFieldMaximalOrder_rational,
} from './order_rational.js';
import type { FunctionFieldPlace, PlaceSet } from './place.js';
import { FunctionFieldPlace_rational } from './place_rational.js';

/**
 * Rational function field in one variable, over an arbitrary base field.
 *
 * @see Reference: sage/rings/function_field/function_field_rational.py:42 (RationalFunctionField)
 */
export class RationalFunctionField<C extends ConstantFieldElement> extends FunctionField<C> {
  readonly _constant_field: ConstantField<C>;
  readonly _names: [string];
  readonly _ring: PolynomialRing<C>;

  private _gen_cache: FunctionFieldElement_rational<C> | null = null;
  private _maximal_order_cache: FunctionFieldMaximalOrder_rational<C> | null = null;
  private _maximal_order_infinite_cache: FunctionFieldMaximalOrderInfinite_rational<C> | null =
    null;

  /**
   * @see Reference: sage/rings/function_field/function_field_rational.py:130 (__init__)
   */
  constructor(constant_field: ConstantField<C>, names: string | [string]) {
    super();
    if (names === null || names === undefined) {
      throw new ValueError('variable name must be specified');
    }
    const nameTuple: [string] = Array.isArray(names) ? names : [names];
    if (typeof constant_field.is_field === 'function' && !constant_field.is_field()) {
      throw new TypeError('constant_field must be a field');
    }
    this._constant_field = constant_field;
    this._names = nameTuple;
    this._ring = new PolynomialRing(constant_field, nameTuple[0]);
  }

  /**
   * @see Reference: sage/rings/function_field/function_field_rational.py:210 (_repr_)
   */
  override _repr_(): string {
    return `Rational function field in ${this.variable_name()} over ${this._constant_field}`;
  }

  override variable_name(): string {
    return this._names[0];
  }

  variable_names(): [string] {
    return this._names;
  }

  /**
   * Coerce ``x`` into an element of the function field.
   *
   * @see Reference: sage/rings/function_field/function_field_rational.py:223 (_element_constructor_)
   */
  override __call__(x: unknown, den?: unknown): FunctionFieldElement_rational<C> {
    if (den !== undefined) {
      return this.__call__(x).div(this.__call__(den));
    }
    if (x instanceof FunctionFieldElement_rational) {
      const e = x as FunctionFieldElement_rational<C>;
      return new FunctionFieldElement_rational(this, e.numerator(), e.denominator(), false);
    }
    if (x instanceof Polynomial) {
      return new FunctionFieldElement_rational(this, this._ring.__call__(x));
    }
    if (Array.isArray(x)) {
      return new FunctionFieldElement_rational(this, this._ring.__call__(x as C[]));
    }
    // constant field element, number or bigint
    return new FunctionFieldElement_rational(this, this._ring.__call__(x as never));
  }

  override zero(): FunctionFieldElement_rational<C> {
    return new FunctionFieldElement_rational(this, this._ring.zero());
  }

  override one(): FunctionFieldElement_rational<C> {
    return new FunctionFieldElement_rational(this, this._ring.one());
  }

  /**
   * Return the ``n``-th generator of the function field.
   *
   * @see Reference: sage/rings/function_field/function_field_rational.py:556 (gen)
   */
  override gen(n: number = 0): FunctionFieldElement_rational<C> {
    if (n !== 0) {
      throw new RangeError('Only one generator.');
    }
    if (this._gen_cache === null) {
      this._gen_cache = new FunctionFieldElement_rational(this, this._ring.gen());
    }
    return this._gen_cache;
  }

  /**
   * @see Reference: sage/rings/function_field/function_field_rational.py:576 (ngens)
   */
  override ngens(): number {
    return 1;
  }

  /**
   * Return the base field of the rational function field, which is the field
   * itself.
   *
   * @see Reference: sage/rings/function_field/function_field_rational.py:588 (base_field)
   */
  override base_field(): RationalFunctionField<C> {
    return this;
  }

  /**
   * @see Reference: sage/rings/function_field/function_field.py:882 (rational_function_field)
   */
  override rational_function_field(): RationalFunctionField<C> {
    return this;
  }

  /**
   * Return the degree over the base field: 1.
   *
   * @see Reference: sage/rings/function_field/function_field_rational.py:532 (degree)
   */
  override degree(base?: FunctionField<C>): bigint {
    if (base !== undefined && base !== this) {
      throw new ValueError('base must be the rational function field itself');
    }
    return 1n;
  }

  /**
   * Return the genus of the function field, namely 0.
   *
   * @see Reference: sage/rings/function_field/function_field_rational.py:743 (genus)
   */
  override genus(): bigint {
    return 0n;
  }

  /**
   * Return the field of which the rational function field is a transcendental
   * extension.
   *
   * @see Reference: sage/rings/function_field/function_field_rational.py:713 (constant_base_field)
   */
  override constant_base_field(): ConstantField<C> {
    return this._constant_field;
  }

  /**
   * Return the underlying polynomial ring `k[x]`.
   *
   * SageMath returns `Frac(k[x])`; this port has no fraction-field type, so we
   * return the polynomial ring whose fraction field it is.
   *
   * @see Reference: sage/rings/function_field/function_field_rational.py:654 (field)
   * @see Deviation: function-field `field()` returns the polynomial ring
   */
  field(): PolynomialRing<C> {
    return this._ring;
  }

  /**
   * Return a polynomial ring in one variable over the rational function field.
   *
   * @see Reference: sage/rings/function_field/function_field_rational.py:434 (polynomial_ring)
   */
  polynomial_ring(_var: string = 'x'): never {
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: RationalFunctionField.polynomial_ring ' +
        '(needs a PolynomialRing over a function field)'
    );
  }

  /**
   * Return the maximal order of the function field, namely `k[x]`.
   *
   * @see Reference: sage/rings/function_field/function_field_rational.py:672 (maximal_order)
   */
  override maximal_order(): FunctionFieldMaximalOrder_rational<C> {
    if (this._maximal_order_cache === null) {
      this._maximal_order_cache = new FunctionFieldMaximalOrder_rational(this);
    }
    return this._maximal_order_cache;
  }

  /** Alias of {@link maximal_order}, as in SageMath. */
  equation_order(): FunctionFieldMaximalOrder_rational<C> {
    return this.maximal_order();
  }

  /**
   * Return the maximal infinite order of the function field.
   *
   * @see Reference: sage/rings/function_field/function_field_rational.py:693 (maximal_order_infinite)
   */
  override maximal_order_infinite(): FunctionFieldMaximalOrderInfinite_rational<C> {
    if (this._maximal_order_infinite_cache === null) {
      this._maximal_order_infinite_cache = new FunctionFieldMaximalOrderInfinite_rational(this);
    }
    return this._maximal_order_infinite_cache;
  }

  /** Alias of {@link maximal_order_infinite}, as in SageMath. */
  equation_order_infinite(): FunctionFieldMaximalOrderInfinite_rational<C> {
    return this.maximal_order_infinite();
  }

  /**
   * Return the different of the rational function field: the zero divisor.
   *
   * @see Reference: sage/rings/function_field/function_field_rational.py:728 (different)
   */
  override different(): FunctionFieldDivisor<C> {
    return this.divisor_group().zero();
  }

  /**
   * Return a field isomorphic to this field with variable ``name``.
   *
   * @see Reference: sage/rings/function_field/function_field_rational.py:755 (change_variable_name)
   */
  change_variable_name(name: string | [string]): RationalFunctionField<C> {
    const n = Array.isArray(name) ? name[0] : name;
    if (n === this.variable_name()) {
      return this;
    }
    return makeRationalFunctionField(this._constant_field, n);
  }

  /**
   * Return the residue field of the place.
   *
   * @see Reference: sage/rings/function_field/function_field_rational.py:799 (residue_field)
   */
  residue_field(
    place: FunctionFieldPlace<C>,
    name?: string
  ): ReturnType<FunctionFieldPlace<C>['residue_field']> {
    return place.residue_field(name);
  }

  override _place_class(parent: PlaceSet<C>, prime: FunctionFieldIdeal<C>): FunctionFieldPlace<C> {
    return new FunctionFieldPlace_rational(parent, prime);
  }
}

/**
 * Rational function fields of characteristic zero.
 *
 * @see Reference: sage/rings/function_field/function_field_rational.py:823 (RationalFunctionField_char_zero)
 */
export class RationalFunctionField_char_zero<
  C extends ConstantFieldElement,
> extends RationalFunctionField<C> {}

/**
 * Rational function field over finite fields.
 *
 * @see Reference: sage/rings/function_field/function_field_rational.py:847 (RationalFunctionField_global)
 */
export class RationalFunctionField_global<
  C extends ConstantFieldElement,
> extends RationalFunctionField<C> {
  /**
   * Return all places of the degree.
   *
   * @see Reference: sage/rings/function_field/function_field_rational.py:853 (places)
   */
  places(degree: number = 1): Array<FunctionFieldPlace<C>> {
    if (degree === 1) {
      return [this.place_infinite(), ...this.places_finite(degree)];
    }
    return this.places_finite(degree);
  }

  /**
   * Return the finite places of the degree.
   *
   * @see Reference: sage/rings/function_field/function_field_rational.py:877 (places_finite)
   */
  places_finite(degree: number = 1): Array<FunctionFieldPlace<C>> {
    return [...this._places_finite(degree)];
  }

  /**
   * Return a generator for the places attached to all monic irreducible
   * polynomials of the given degree.
   *
   * The enumeration order is SageMath's: `R.polynomials(max_degree=degree-1)`
   * runs a little-endian odometer over the constant field's own iteration
   * order, with the *constant* coefficient varying fastest
   * (`reference/sage/src/sage/rings/polynomial/polynomial_ring.py:1548`).
   *
   * @see Reference: sage/rings/function_field/function_field_rational.py:893 (_places_finite)
   */
  *_places_finite(degree: number = 1): IterableIterator<FunctionFieldPlace<C>> {
    if (degree < 1) {
      throw new ValueError('degree must be a positive integer');
    }
    const O = this.maximal_order();
    const R = O._ring;
    const k = this.constant_base_field();
    const els = constant_field_element_list(k);
    const q = BigInt(els.length);

    const lm = R.monomial(degree);
    const n = degree; // number of free coefficients: x^0 .. x^{degree-1}
    const total = q ** BigInt(n);
    for (let i = 0n; i < total; i++) {
      const coeffs: C[] = [];
      let t = i;
      for (let j = 0; j < n; j++) {
        coeffs.push(els[Number(t % q)]!);
        t /= q;
      }
      const g = new Polynomial<C>(coeffs, R);
      const h = lm.add(g);
      if (h.is_irreducible()) {
        yield O.ideal(h).place();
      }
    }
  }

  /**
   * Return the unique place at infinity.
   *
   * @see Reference: sage/rings/function_field/function_field_rational.py:916 (place_infinite)
   */
  place_infinite(): FunctionFieldPlace<C> {
    return this.maximal_order_infinite().prime_ideal().place();
  }

  /**
   * Return a place of ``degree``.
   *
   * @see Reference: sage/rings/function_field/function_field_rational.py:928 (get_place)
   */
  get_place(degree: number): FunctionFieldPlace<C> {
    for (const p of this._places_finite(degree)) {
      return p;
    }
    throw new ValueError('there is a bug around');
  }

  /** Number of elements of the constant field. */
  constant_field_order(): bigint {
    return constant_field_cardinality(this.constant_base_field());
  }
}

/**
 * Build the right `RationalFunctionField` subclass for ``constant_field``.
 *
 * @see Reference: sage/rings/function_field/constructor.py:104 (FunctionFieldFactory.create_object)
 */
export function makeRationalFunctionField<C extends ConstantFieldElement>(
  constant_field: ConstantField<C>,
  names: string | [string]
): RationalFunctionField<C> {
  if (constant_field_is_finite(constant_field)) {
    return new RationalFunctionField_global(constant_field, names);
  }
  const c = (constant_field as { characteristic?: unknown }).characteristic;
  const char = typeof c === 'function' ? (c as () => bigint).call(constant_field) : c;
  if (char === 0n || char === 0) {
    return new RationalFunctionField_char_zero(constant_field, names);
  }
  return new RationalFunctionField(constant_field, names);
}
