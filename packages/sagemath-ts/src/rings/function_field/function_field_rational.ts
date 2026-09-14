/**
 * @module sage/rings/function_field/function_field_rational
 * @description Function Fields: rational
 *
 * Port of: sage/rings/function_field/function_field_rational.py
 */

import { AssertionError, IndexError, NotImplementedError, ValueError } from '../../errors.js';
import type { FractionField_generic } from '../fraction_field.js';
import { Polynomial } from '../polynomial/polynomial_element.js';
import { PolynomialRing } from '../polynomial/polynomial_ring.js';
import { constant_field_cardinality, constant_field_is_finite } from './constant_field.js';
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
  readonly _field: FractionField_generic<C>;

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
    const nameTuple: [string] = Array.isArray(names) ? [names[0]] : [names];
    if (typeof constant_field.is_field === 'function' && !constant_field.is_field()) {
      throw new TypeError('constant_field must be a field');
    }
    this._constant_field = constant_field;
    this._names = nameTuple;
    this._ring = new PolynomialRing(constant_field, nameTuple[0]);
    this._field = this._ring.fraction_field();
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
    return [...this._names];
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
      if (x.parent === this) return x;
      x = x.element();
    }
    return new FunctionFieldElement_rational(this, this._field.__call__(x));
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
      throw new IndexError('Only one generator.');
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
   * Return the underlying fraction field `Frac(k[x])`.
   *
   * @see Reference: sage/rings/function_field/function_field_rational.py:654 (field)
   */
  field(): FractionField_generic<C> {
    return this._field;
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
   * Return [isomorphic field, map to this field, inverse map] for variable ``name``.
   *
   * @see Reference: sage/rings/function_field/function_field_rational.py:755 (change_variable_name)
   */
  change_variable_name(
    name: string | [string]
  ): [
    RationalFunctionField<C>,
    (x: unknown) => FunctionFieldElement_rational<C>,
    (x: unknown) => FunctionFieldElement_rational<C>,
  ] {
    if (Array.isArray(name) && name.length !== 1) {
      throw new ValueError('names must be a tuple with a single string');
    }
    const n = Array.isArray(name) ? name[0] : name;
    // Map.__call__ converts into the domain before invoking even an identity map.
    const convert = (source: RationalFunctionField<C>, x: unknown) => {
      try {
        return source.__call__(x);
      } catch (error) {
        if (!(error instanceof TypeError || error instanceof NotImplementedError)) throw error;
        const display = (value: unknown, nested = false): string => {
          if (value == null) return 'None';
          if (typeof value === 'boolean') return value ? 'True' : 'False';
          if (typeof value === 'string' && nested) {
            const quoted = JSON.stringify(value);
            if (value.includes("'") && !value.includes('"')) return quoted;
            return "'" + quoted.slice(1, -1).replace(/\\"/g, '"').replace(/'/g, "\\'") + "'";
          }
          if (Array.isArray(value)) return `[${value.map(v => display(v, true)).join(', ')}]`;
          if (value instanceof Map)
            return `{${[...value].map(([k, v]) => `${display(k, true)}: ${display(v, true)}`).join(', ')}}`;
          if (typeof value === 'object' && (Object.getPrototypeOf(value) === Object.prototype
            || Object.getPrototypeOf(value) === null))
            return `{${Object.entries(value).map(([k, v]) => `${/^-?\d+$/.test(k) ? k : display(k, true)}: ${display(v, true)}`).join(', ')}}`;
          return String(value);
        };
        const value = display(x);
        throw new TypeError(`${value} fails to convert into the map's domain ${source}, but a \`pushforward\` method is not properly implemented`);
      }
    };
    if (n === this.variable_name()) {
      const identity = (x: unknown) => convert(this, x);
      return [this, identity, identity];
    }
    const L = makeRationalFunctionField(this._constant_field, n);
    const map = (target: RationalFunctionField<C>, x: FunctionFieldElement_rational<C>) =>
      new FunctionFieldElement_rational(
        target,
        target._ring.__call__(x.numerator().coeffs),
        target._ring.__call__(x.denominator().coeffs)
      );
    return [L, (x) => map(this, convert(L, x)), (x) => map(L, convert(this, x))];
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
    const O = this.maximal_order();
    const R = O._ring;
    const G = R.polynomials({ max_degree: degree - 1 });
    const lm = R.monomial(degree);
    for (const g of G) {
      const h = lm.add(g);
      if (h.is_irreducible()) yield O.ideal(h).place();
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
    throw new AssertionError('there is a bug around');
  }

  /** Number of elements of the constant field. */
  constant_field_order(): bigint {
    return constant_field_cardinality(this.constant_base_field());
  }
}

// Sage's UniqueFactory keys by the constant-field parent and variable name.
const rationalFunctionFields = new WeakMap<object, Map<string, unknown>>();

/**
 * Build the right `RationalFunctionField` subclass for ``constant_field``.
 *
 * @see Reference: sage/rings/function_field/constructor.py:104 (FunctionFieldFactory.create_object)
 */
export function makeRationalFunctionField<C extends ConstantFieldElement>(
  constant_field: ConstantField<C>,
  names: string | [string]
): RationalFunctionField<C> {
  const name = Array.isArray(names) ? names[0] : names;
  let cache = rationalFunctionFields.get(constant_field);
  const cached = cache?.get(name);
  if (cached) return cached as RationalFunctionField<C>;
  let field: RationalFunctionField<C>;
  if (constant_field_is_finite(constant_field)) {
    field = new RationalFunctionField_global(constant_field, names);
  } else {
    const c = (constant_field as { characteristic?: unknown }).characteristic;
    const char = typeof c === 'function' ? (c as () => bigint).call(constant_field) : c;
    field =
      char === 0n || char === 0
        ? new RationalFunctionField_char_zero(constant_field, names)
        : new RationalFunctionField(constant_field, names);
  }
  if (!cache) {
    cache = new Map();
    rationalFunctionFields.set(constant_field, cache);
  }
  cache.set(name, field);
  return field;
}
