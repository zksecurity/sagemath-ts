/**
 * @see Deviation: Polynomial String Parsing
 * @module sage/rings/polynomial/polynomial_ring
 * @description Polynomial rings over arbitrary coefficient rings
 *
 * Port of: sage/rings/polynomial/polynomial_ring.py
 */

import { factor, factorial } from '../../arith/misc.js';
import {
  ArithmeticError,
  AttributeError,
  IndexError,
  NotImplementedError,
  OverflowError,
  TypeError,
  ValueError,
  ZeroDivisionError,
} from '../../errors.js';
import { _xmrange_iter } from '../../misc/mrange.js';
import { LookupNameMaker, Parser } from '../../misc/parser.js';
import type { ConstantField, ConstantFieldElement } from '../function_field/constant_field.js';
import { FunctionFieldElement_rational } from '../function_field/element_rational.js';
import { RationalFunctionField } from '../function_field/function_field_rational.js';
import {
  FractionField_generic,
  FractionField_1poly_field,
  _fraction_field_domain,
  _fraction_characteristic,
} from '../fraction_field.js';
import { FpT, FpTElement } from '../fraction_field_FpT.js';
import { isFractionElement, _fraction_polynomial } from '../fraction_field_element.js';
import { ZZ } from '../integer_ring.js';
import { Rational } from '../rational.js';
import {
  PrimeField,
  PrimeFieldElement,
  FiniteFieldElement as ExtensionElement,
} from '../finite_rings/finite_field_extension.js';
import { GF2Element } from '../finite_rings/gf2.js';
import { IntegerMod } from '../finite_rings/integer_mod.js';
import { FiniteFieldElement as LegacyPrimeElement } from '../finite_rings/finite_field_prime.js';
import {
  type CoefficientRing,
  Polynomial,
  _polynomial_coercion,
  type PolynomialRingBase,
  type RingElement,
} from './polynomial_element.js';

/**
 * Interface for field elements that support division/inverse.
 */
interface FieldElement extends RingElement {
  inv(): this;
  div?(other: this): this;
}

/**
 * A polynomial ring R[x] over a base ring R.
 */
export class PolynomialRing<C extends RingElement> implements PolynomialRingBase<C> {
  readonly base_ring: CoefficientRing<C>;
  readonly variable_name: string;
  private _generator?: Polynomial<C>;
  private _fraction_field?: FractionField_generic<C>;
  has_coerce_map_from(other: CoefficientRing<RingElement>): boolean {
    return _polynomial_coercion.canCoerce(this, other);
  }
  fraction_field(): FractionField_generic<C> {
    if (this._fraction_field) return this._fraction_field;
    if (!_fraction_field_domain(this)) {
      const z = this.base_ring.zero();
      const cls =
        z instanceof IntegerMod && z.modulus !== 1n
          ? 'PolynomialRing_dense_mod_n'
          : 'PolynomialRing_commutative';
      throw new AttributeError(`'${cls}_with_category' object has no attribute 'fraction_field'`);
    }
    const p = _fraction_characteristic(this.base_ring),
      c = this.base_ring.zero();
    return (this._fraction_field = this.base_ring.is_field?.()
      ? p > 2n &&
        p < 46341n &&
        (c instanceof PrimeFieldElement ||
          c instanceof LegacyPrimeElement ||
          c instanceof IntegerMod)
        ? new FpT(this)
        : new FractionField_1poly_field(this)
      : new FractionField_generic(this));
  }

  constructor(base_ring: CoefficientRing<C>, variable_name: string = 'x') {
    this.base_ring = base_ring;
    this.variable_name = variable_name;
  }

  /**
   * Create a polynomial from coefficients, an exponent dictionary, or an arithmetic string.
   * @see Deviation: Negative Polynomial Monomial Degrees
   * @see Deviation: Polynomial Sylvester Matrices and Explicit Scalar Construction
   */
  __call__(x?: unknown): Polynomial<C> {
    if (x === undefined || x === null) return this.zero();
    if (typeof x === 'string') return this._parse_string(x);
    if (isFractionElement(x)) {
      const origin = x.parent.ring();
      if (_polynomial_coercion.canCoerce(this.base_ring, x.parent))
        return new Polynomial([this.base_ring.__call__(x)], this);
      const same =
        this.variable_name === origin.variable_name &&
        _polynomial_coercion.canCoerce(this.base_ring, origin.base_ring) &&
        _polynomial_coercion.canCoerce(origin.base_ring, this.base_ring);
      if (x instanceof FpTElement) {
        if (this.base_ring instanceof PolynomialRing)
          return new Polynomial([this.base_ring.__call__(x) as C], this);
        if (same) return this.__call__(FpTElement._section(x, false));
        const zero = this.base_ring.zero();
        if (zero instanceof Rational) throw new TypeError(`unable to convert ${x} to a rational`);
        if (zero instanceof ExtensionElement) throw new TypeError('no coercion defined');
        let c: bigint;
        try {
          c = FpTElement._section(x, true) as bigint;
        } catch (e) {
          if (!(e instanceof ValueError)) throw e;
          if (
            String(this.base_ring) === 'Integer Ring' ||
            _fraction_characteristic(this.base_ring) === x.parent.p
          )
            throw e;
          throw new TypeError(`unable to convert ${x} to a rational`);
        }
        return this.__call__(c);
      }
      const polynomial = _fraction_polynomial(x, same);
      return same ? (polynomial as Polynomial<C>) : this.__call__(polynomial);
    }

    if (
      (x instanceof PrimeFieldElement ||
        x instanceof IntegerMod ||
        x instanceof LegacyPrimeElement) &&
      x.isZero()
    )
      return this.zero();
    // Parent.__call__ first uses canonical coefficient maps. Without one,
    // polynomial_ring.py:435-439 expands finite elements via polynomial().
    if (x instanceof ExtensionElement && !_polynomial_coercion.canCoerce(this.base_ring, x.parent))
      return this.__call__(x.lift);
    if (x instanceof GF2Element) {
      if (this.base_ring.zero() instanceof ExtensionElement)
        return new Polynomial(
          x.isZero() ? [] : [this.base_ring.__call__(new PrimeField(2n).__call__(x.toBigInt()))],
          this
        );
      x = x.toBigInt();
    }
    const coefficient = (c: unknown): C => {
      if (c instanceof GF2Element)
        c =
          this.base_ring.zero() instanceof ExtensionElement
            ? new PrimeField(2n).__call__(c.toBigInt())
            : c.toBigInt();
      if (this.base_ring.zero() instanceof GF2Element)
        return this.base_ring.__call__(new PrimeField(2n).__call__(c as bigint).value);
      return this.base_ring.__call__(c);
    };
    if (x instanceof Map || (typeof x === 'object' && x !== null
      && (Object.getPrototypeOf(x) === Object.prototype || Object.getPrototypeOf(x) === null))) {
      // Python coefficient dictionaries use explicit exponent keys. Records are
      // the numeric-key adapter; Maps also retain tuple and nonnumeric keys.
      const entries: Array<[unknown, unknown]> = x instanceof Map ? [...x.entries()]
        : Object.entries(x).map(([key, value]) => [/^-?\d+$/.test(key) ? BigInt(key) : key, value]);
      if (entries.length === 0) return this.zero();
      const zero = this.base_ring.zero();
      const rational = zero instanceof Rational;
      const integer = String(this.base_ring) === 'Integer Ring';
      const p = _fraction_characteristic(this.base_ring);
      const template = zero instanceof GF2Element || zero instanceof ExtensionElement
        || (p > 1n && p < 1n << 64n);
      const exponent = (key: unknown): number => {
        if (!rational && !template && Array.isArray(key)) {
          if (key.length === 0) throw new IndexError('tuple index out of range');
          key = key[0];
        }
        if (typeof key === 'number' && !Number.isInteger(key) && !rational)
          throw new TypeError(integer || template ? "'float' object cannot be interpreted as an integer"
            : "can't multiply sequence by non-int of type 'float'");
        if (typeof key !== 'bigint' && typeof key !== 'number') {
          if (rational) throw new TypeError('an integer is required');
          if (!integer && !template && typeof key === 'string')
            throw new TypeError('can only concatenate str (not "int") to str');
          throw new TypeError(`'${Array.isArray(key) ? 'tuple' : 'str'}' object cannot be interpreted as an integer`);
        }
        const exact = typeof key === 'bigint' ? key : BigInt(Math.trunc(key));
        if (exact < -(1n << 63n) || exact >= 1n << 63n)
          throw new OverflowError(rational ? 'Python int too large to convert to C long'
            : integer || template ? 'Python int too large to convert to C ssize_t'
            : "cannot fit 'int' into an index-sized integer");
        const n = Number(exact);
        if (integer && n < 0) throw new ValueError(`Negative monomial degrees not allowed: ${key}`);
        return n;
      };
      const dictCoefficient = (value: unknown): C => {
        if (rational && Array.isArray(value))
          throw new TypeError(`unable to convert [${value.map(String).join(', ')}] to a rational`);
        if (value instanceof Polynomial && value.degree() > 0)
          throw new TypeError(rational ? 'cannot convert nonconstant polynomial' : 'not a constant polynomial');
        return coefficient(value);
      };
      if (template) {
        // polynomial_template.pxi accumulates element_class(coef) * gen^degree.
        let result = this.zero();
        for (const [key, value] of entries) {
          const monomial = this.monomial(exponent(key));
          // Direct element construction treats strings as scalars, not expressions.
          result = result.add(monomial.mul(this.__call__(typeof value === 'string' ? coefficient(value) : value)));
        }
        return result;
      }
      if (rational) {
        const coefficients: C[] = [];
        for (const [key, value] of entries) {
          const n = exponent(key);
          if (n < 0) throw new RangeError('negative monomial degree is outside the native FLINT contract');
          const c = dictCoefficient(value);
          while (coefficients.length <= n) coefficients.push(zero as C);
          coefficients[n] = c;
        }
        return new Polynomial(coefficients, this);
      }
      if (!integer) {
        // _dict_to_list calls max(keys) before allocating or converting coefficients.
        let maxKey = entries[0]![0];
        for (const [key] of entries.slice(1)) {
          const numeric = (v: unknown) => typeof v === 'bigint' || typeof v === 'number';
          if (numeric(key) !== numeric(maxKey)) {
            const other = numeric(key) ? maxKey : key;
            const type = Array.isArray(other) ? 'tuple' : 'str';
            throw new TypeError(`unsupported operand parent(s) for ${numeric(key) ? '>' : '<'}: 'Integer Ring' and '<class '${type}'>'`);
          }
          if (numeric(key) ? Number(key) > Number(maxKey) : String(key) > String(maxKey)) maxKey = key;
        }
      }
      const exponents = entries.map(([key]) => exponent(key));
      const coefficients = Array.from({ length: Math.max(0, Math.max(...exponents) + 1) }, () => zero as C);
      for (let i = 0; i < entries.length; i++) {
        let n = exponents[i]!;
        if (n < 0) {
          if (rational) throw new RangeError('negative monomial degree is outside the native FLINT contract');
          n += coefficients.length;
          if (n < 0) throw new IndexError('list assignment index out of range');
        }
        // Generic dense constructors assign all raw values before coefficient conversion.
        coefficients[n] = (integer ? dictCoefficient(entries[i]![1]) : entries[i]![1]) as C;
      }
      return new Polynomial(integer ? coefficients : coefficients.map(dictCoefficient), this);
    }
    if (x instanceof Polynomial) {
      // polynomial_ring.py:388-403 checks the exact parent before coefficient
      // conversion, including constant embedding from a polynomial base ring.
      if (x.parent === this) return x;
      if ((x.parent as unknown) === this.base_ring)
        return new Polynomial([x as unknown as C], this);
      if (_polynomial_coercion.canCoerce(x.parent.base_ring, this)) {
        if (x.degree() > 0) throw new TypeError(`${x} is not a constant polynomial`);
        return this.__call__(x.getCoeff(0));
      }
      // Convert polynomial from potentially different ring
      return new Polynomial(x.coeffs.map(coefficient), this);
    }

    if (Array.isArray(x)) {
      const zero = this.base_ring.zero();
      // FLINT's QQ constructor unwraps singleton lists recursively.
      if (zero instanceof Rational && x.length === 1) {
        // Native QQ singleton recursion stays inside the element constructor;
        // it does not re-enter the parent's finite-ring polynomial hook.
        if (x[0] instanceof ExtensionElement)
          return new Polynomial([this.base_ring.__call__(x[0])], this);
        return this.__call__(x[0]);
      }
      // The GF2X template constructs each list entry as a polynomial and
      // accumulates it at its index, so nested lists can have positive degree.
      const base = this.base_ring as CoefficientRing<C> & { characteristic?: bigint };
      if (base.characteristic === 2n && 'value' in zero && x.some(Array.isArray)) {
        let result = this.zero();
        for (let i = 0; i < x.length; i++) result = result.add(this.__call__(x[i]).shift(i));
        return result;
      }
      // polynomial_generic_dense.__init__ and polynomial_zmod_flint.__init__
      // convert every coefficient before stripping trailing zeros.
      return new Polynomial(
        x.map((c) => {
          if (zero instanceof Rational && Array.isArray(c)) {
            // These are Python coefficient lists, not QQ's explicit tuple API.
            while (Array.isArray(c) && c.length === 1) c = c[0];
            if (Array.isArray(c)) {
              const repr = (v: unknown): string =>
                Array.isArray(v)
                  ? `[${v.map(repr).join(', ')}]`
                  : v == null
                    ? 'None'
                    : typeof v === 'string'
                      ? `'${v}'`
                      : typeof v === 'boolean'
                        ? v
                          ? 'True'
                          : 'False'
                        : String(v);
              throw new TypeError(`unable to convert ${repr(c)} to a rational`);
            }
            // Rational([None]) differs from Rational(None).
            if (c == null) throw new TypeError('unable to convert None to a rational');
          }
          return coefficient(c);
        }),
        this
      );
    }

    // Single coefficient (constant polynomial)
    const coeff = coefficient(x);
    return new Polynomial([coeff], this);
  }

  /** Polynomial-ring string constructor, polynomial_ring.py:417-425. */
  /** @internal Shared expression evaluator for fraction_field.py string fallback. */
  _parse_fraction_string(source: string): unknown {
    // sage_eval compiles the complete expression before resolving its names.
    // Parse once without arithmetic or name lookup to preserve that ordering.
    // Its preparser does not enable polynomial implicit multiplication.
    source = source.replace(/\\\n/g, '');
    for (let i = 0; i < source.length; i++) {
      const quote = source[i];
      if (quote === '\\') {
        const line = source.slice(0, i).split('\n').length;
        throw new SyntaxError(
          `unexpected character after line continuation character (<string>, line ${line})`
        );
      }
      if (quote !== "'" && quote !== '"') continue;
      const line = source.slice(0, i).split('\n').length;
      let j = i + 1;
      for (; j < source.length && source[j] !== quote && source[j] !== '\n'; j++)
        if (source[j] === '\\') j++;
      if (j >= source.length || source[j] === '\n')
        throw new SyntaxError(
          `unterminated string literal (detected at line ${line}) (<string>, line ${line})`
        );
      i = j;
    }
    const syntax = new Parser(
      () => 0,
      () => 0,
      () => 0,
      {},
      false,
      { binary: () => 0, unary: () => 0 }
    );
    try {
      syntax.parse(source);
    } catch (e) {
      if (e instanceof SyntaxError) {
        const position = (e as SyntaxError & { position?: number }).position ?? source.length;
        const line = Array.from(source).slice(0, position).join('').split('\n').length;
        throw new SyntaxError(`invalid syntax (<string>, line ${line})`);
      } else if (!(e instanceof NotImplementedError) && (e as Error).name !== 'NameError') throw e;
    }
    return this._parse_string(source, true);
  }

  private _parse_string(source: string): Polynomial<C>;
  private _parse_string(source: string, fractionResult: true): unknown;
  private _parse_string(source: string, fractionResult = false): unknown {
    type Fraction = FunctionFieldElement_rational<ConstantFieldElement>;
    type Value = bigint | boolean | Rational | C | Polynomial<C> | Fraction;
    const base = this.base_ring as CoefficientRing<C> & {
      characteristic?: bigint | (() => bigint);
    };
    const p =
      typeof base.characteristic === 'function'
        ? base.characteristic()
        : (base.characteristic ?? 0n);
    // polynomial_ring.py:3640-3650 selects FpT only for small odd prime fields.
    const zero = this.base_ring.zero();
    const primeBase = typeof zero === 'object' && zero !== null && 'value' in zero;
    const specialized = primeBase && p > 2n && p < 46341n;
    let fractionField: RationalFunctionField<ConstantFieldElement> | undefined;
    const rational = (x: Value): Rational =>
      x instanceof Rational
        ? x
        : new Rational(typeof x === 'boolean' ? (x ? 1n : 0n) : (x as bigint));
    const scalar = (x: Value): boolean =>
      typeof x === 'bigint' || typeof x === 'boolean' || x instanceof Rational;
    const polynomial = (x: Value): Polynomial<C> => this.__call__(x as C | Polynomial<C>);
    const fraction = (x: Value): Fraction => {
      if (x instanceof FunctionFieldElement_rational) return x;
      if (!fractionField)
        fractionField = new RationalFunctionField(
          this.base_ring as unknown as ConstantField<ConstantFieldElement>,
          this.variable_name
        );
      return fractionField.__call__(x);
    };
    const finish = (x: Value): unknown => {
      if (fractionResult) {
        if (!(x instanceof FunctionFieldElement_rational)) return x;
        // The temporary evaluator parent must not escape into the result.
        // Reparent its polynomials while preserving the fraction's reduction state.
        const F = this.fraction_field();
        return new F._element_class(F, x.numerator(), x.denominator(), { reduce: false });
      }
      if (x instanceof FunctionFieldElement_rational) {
        if (x.denominator().degree() > 0) {
          if (specialized) throw new ValueError('not integral');
          throw new TypeError('fraction must have unit denominator');
        }
        return this.__call__(x.numerator() as unknown as Polynomial<C>);
      }
      return polynomial(x);
    };
    const binary = (op: string, left: Value, right: Value): Value => {
      if (op === '^') {
        let exponent: bigint;
        try {
          exponent = ZZ.__call__(right as Parameters<typeof ZZ.__call__>[0]);
        } catch (error) {
          if (left instanceof Polynomial || left instanceof FunctionFieldElement_rational) {
            if (p < 2n ** 63n && right instanceof Polynomial && !right.isConstant())
              throw new TypeError('cannot convert nonconstant polynomial');
            throw new TypeError(
              primeBase && p < 2n ** 63n
                ? 'Only integral powers defined.'
                : 'non-integral exponents not supported'
            );
          }
          throw error;
        }
        if (scalar(left)) {
          const result = rational(left).pow(exponent);
          return typeof left === 'bigint' && exponent >= 0n ? result.numerator : result;
        }
        if (left instanceof Polynomial) {
          // polynomial_template.pxi:646 returns zero before reciprocation.
          if (left.isZero() && primeBase && p < 2n ** 63n)
            return exponent === 0n ? this.one() : this.zero();
          if (exponent >= 0n) return left.pow(exponent) as Polynomial<C>;
          if (left.isZero() && primeBase) throw new ZeroDivisionError('Inverse does not exist.');
        }
        return fraction(left).pow(exponent);
      }
      const comparisons = ['=', 'NOT_EQ', '<', 'LESS_EQ', '>', 'GREATER_EQ'];
      if (comparisons.includes(op)) {
        const cmp =
          scalar(left) && scalar(right)
            ? rational(left).cmp(rational(right))
            : fraction(left).cmp(fraction(right));
        return op === '='
          ? cmp === 0
          : op === 'NOT_EQ'
            ? cmp !== 0
            : op === '<'
              ? cmp < 0
              : op === 'LESS_EQ'
                ? cmp <= 0
                : op === '>'
                  ? cmp > 0
                  : cmp >= 0;
      }
      if (scalar(left) && scalar(right)) {
        const a = rational(left);
        const b = rational(right);
        const result =
          op === '+' ? a.add(b) : op === '-' ? a.sub(b) : op === '*' ? a.mul(b) : a.div(b);
        return op !== '/' && !(left instanceof Rational) && !(right instanceof Rational)
          ? result.numerator
          : result;
      }
      // QQ has no canonical map into positive-characteristic polynomial rings,
      // even when the particular rational happens to be integral.
      if (p > 0n && (left instanceof Rational || right instanceof Rational)) {
        const parent = (x: Value) =>
          x instanceof Rational
            ? 'Rational Field'
            : x instanceof FunctionFieldElement_rational
              ? 'Fraction Field of ' + this.toString()
              : this.toString();
        throw new TypeError(
          'unsupported operand parent(s) for ' +
            op +
            ": '" +
            parent(left) +
            "' and '" +
            parent(right) +
            "'"
        );
      }
      if (
        left instanceof FunctionFieldElement_rational ||
        right instanceof FunctionFieldElement_rational
      ) {
        const a = fraction(left);
        const b = fraction(right);
        if (op === '+') return a.add(b);
        if (op === '-') return a.sub(b);
        if (op === '*') return a.mul(b);
        if (b.is_zero())
          throw new ZeroDivisionError(specialized ? '' : 'fraction field element division by zero');
        return a.div(b);
      }
      const a = polynomial(left);
      const b = polynomial(right);
      if (op === '+') return a.add(b);
      if (op === '-') return a.sub(b);
      if (op === '*') return a.mul(b);
      // Polynomial scalar division stays in the polynomial ring, including
      // the coefficient ring's own inverse-of-zero diagnostic.
      if (b.isConstant() && !(right instanceof Polynomial)) {
        const c = b.getCoeff(0) as C & { inv?: () => C; div?: (x: C) => C };
        if (c.inv) return a.scalar_mul(c.inv());
        if (c.div) return a.scalar_mul((this.base_ring.one() as typeof c).div!(c));
      }
      if (b.isZero())
        throw new ZeroDivisionError(
          specialized ? 'fraction has denominator 0' : 'fraction field element division by zero'
        );
      return fraction(a).div(fraction(b));
    };
    const unary = (op: string, x: Value): Value => {
      if (op === '!') {
        if (scalar(x)) return factorial(ZZ.__call__(x as Parameters<typeof ZZ.__call__>[0]));
        throw new TypeError(
          'cannot coerce arguments: positive characteristic not allowed in symbolic computations'
        );
      }
      if (typeof x === 'bigint') return -x;
      if (typeof x === 'boolean') return x ? -1n : 0n;
      return x.neg() as Value;
    };
    const names = new LookupNameMaker<Value>(
      { [this.variable_name]: fractionResult ? fraction(this.gen()) : this.gen() },
      fractionResult ? undefined : (name) => this.base_ring.__call__(name)
    );
    const parser = new Parser<Value>(
      (s) => ZZ.__call__(s),
      (s) => this.base_ring.__call__(s),
      (name) => names.__call__(name),
      {},
      !fractionResult,
      { binary, unary }
    );
    try {
      return finish(parser.parse(source));
    } catch (error) {
      if (!fractionResult && (error as Error).name === 'NameError')
        throw new TypeError('Unable to coerce string');
      throw error;
    }
  }

  /**
   * Return the zero polynomial.
   */
  zero(): Polynomial<C> {
    return new Polynomial([], this);
  }

  /**
   * Return the one polynomial (constant 1).
   */
  one(): Polynomial<C> {
    return new Polynomial([this.base_ring.one() as C], this);
  }

  /**
   * Return the generator (the variable x).
   */
  gen(n: unknown = 0): Polynomial<C> {
    // The cached Sage method rejects unhashable indices before its n != 0 test.
    if (Array.isArray(n)) throw new TypeError("unhashable type: 'list'");
    const comparable =
      n !== null && typeof n === 'object' ? (n as { eq?: (other: bigint) => unknown }) : undefined;
    if (
      !(
        n === 0 ||
        n === 0n ||
        n === false ||
        (typeof comparable?.eq === 'function' && comparable.eq(0n) === true)
      )
    )
      throw new IndexError('generator n not defined');
    return (this._generator ??= new Polynomial(
      [this.base_ring.zero() as C, this.base_ring.one() as C],
      this,
      true
    ));
  }

  /**
   * Create a polynomial from a list of (coefficient, degree) pairs.
   */
  fromTerms(terms: Array<[C, number]>): Polynomial<C> {
    if (terms.length === 0) {
      return this.zero();
    }

    const maxDeg = Math.max(...terms.map(([_, d]) => d));
    const coeffs: C[] = [];

    for (let i = 0; i <= maxDeg; i++) {
      coeffs.push(this.base_ring.zero() as C);
    }

    for (const [c, d] of terms) {
      coeffs[d] = coeffs[d]!.add(c) as C;
    }

    return new Polynomial(coeffs, this);
  }

  /**
   * Return the monomial x^n.
   * @see Deviation: Negative Polynomial Monomial Degrees
   */
  monomial(n: number): Polynomial<C> {
    if (n < 0) {
      const one = this.base_ring.one();
      if (one instanceof GF2Element || one instanceof ExtensionElement) return this.zero();
      const p =
        one instanceof IntegerMod
          ? one.modulus
          : one instanceof PrimeFieldElement
            ? one.parent.characteristic
            : one instanceof LegacyPrimeElement
              ? one.p
              : undefined;
      // NTL GF2X/ZZ_pEX discard negative monomials. FLINT's dictionary
      // constructor passes the negative exponent to unchecked native kernels.
      if (p === 2n) return this.zero();
      if (one instanceof Rational || (p !== undefined && p > 2n && p < 1n << 64n))
        throw new RangeError('negative monomial degree is outside the native FLINT contract');
      throw new IndexError('list assignment index out of range');
    }
    const coeffs: C[] = [];
    for (let i = 0; i < n; i++) {
      coeffs.push(this.base_ring.zero() as C);
    }
    coeffs.push(this.base_ring.one() as C);
    return new Polynomial(coeffs, this);
  }

  /** Enumerate polynomials with exactly one degree bound, as in Sage. */
  polynomials(
    options: { of_degree?: number | null; max_degree?: number | null } = {}
  ): IterableIterator<Polynomial<C>> {
    this._check_enumeration_base();
    const { of_degree, max_degree } = options;
    if (of_degree != null && max_degree == null) return this._polys_degree(of_degree);
    if (max_degree != null && of_degree == null) return this._polys_max(max_degree);
    throw new ValueError('you should pass exactly one of of_degree and max_degree');
  }

  /** Enumerate monic polynomials with exactly one degree bound, as in Sage. */
  monics(
    options: { of_degree?: number | null; max_degree?: number | null } = {}
  ): IterableIterator<Polynomial<C>> {
    this._check_enumeration_base();
    const { of_degree, max_degree } = options;
    if (of_degree != null && max_degree == null) return this._monics_degree(of_degree);
    if (max_degree != null && of_degree == null) return this._monics_max(max_degree);
    throw new ValueError('you should pass exactly one of of_degree and max_degree');
  }

  private _check_enumeration_base(): void {
    const base = this.base_ring as CoefficientRing<C> & {
      order?: unknown;
      is_finite?: () => boolean;
    };
    const order = typeof base.order === 'function' ? base.order() : base.order;
    if (order === Infinity || order === 'Infinity' || base.is_finite?.() === false)
      throw new NotImplementedError('');
  }

  *_polys_degree(of_degree: number): IterableIterator<Polynomial<C>> {
    const base = this.base_ring as CoefficientRing<C> & Iterable<C>;
    const zero = base.zero();
    for (const leading of base) {
      if (!leading.eq(zero)) {
        for (const lower of _xmrange_iter(
          Array.from({ length: Math.max(0, of_degree) }, () => base)
        ))
          yield this.__call__([leading, ...lower].reverse());
      }
    }
  }

  *_polys_max(max_degree: number): IterableIterator<Polynomial<C>> {
    const base = this.base_ring as CoefficientRing<C> & Iterable<C>;
    for (const coeffs of _xmrange_iter(
      Array.from({ length: Math.max(0, max_degree + 1) }, () => base)
    ))
      yield this.__call__(coeffs.reverse());
  }

  *_monics_degree(of_degree: number): IterableIterator<Polynomial<C>> {
    const base = this.base_ring as CoefficientRing<C> & Iterable<C>;
    for (const coeffs of _xmrange_iter<C>([
      [base.one()],
      ...Array.from({ length: Math.max(0, of_degree) }, () => base),
    ]))
      yield this.__call__(coeffs.reverse());
  }

  *_monics_max(max_degree: number): IterableIterator<Polynomial<C>> {
    for (let degree = 0; degree <= max_degree; degree++) yield* this._monics_degree(degree);
  }

  /**
   * Check if this is a polynomial ring (yes).
   */
  is_ring(): boolean {
    return true;
  }

  /**
   * Check if this is a field (no, unless base is trivial).
   */
  is_field(): boolean {
    return false;
  }

  /**
   * Return the Lagrange interpolation polynomial through the given points.
   *
   * This computes the unique polynomial P of degree at most n-1 such that
   * P(x_i) = y_i for all given points (x_i, y_i).
   *
   * With `algorithm='neville'` the **whole last row of the Neville table** is
   * returned (a list of polynomials), exactly as Sage does; the interpolating
   * polynomial is its last entry.
   *
   * @param points - Array of [x, y] pairs where the x values must be distinct
   *   (`x_i - x_j` must be invertible); as in Sage there is no precheck, a
   *   repeated x surfaces as a {@link ZeroDivisionError} from the base ring.
   * @param algorithm - `'divided_difference'` (default), `'neville'` or `'pari'`
   * @param previous_row - only used with `'neville'`: the last row of a
   *   previous Neville computation, so that it can be extended incrementally
   * @returns The interpolating polynomial, or the last Neville row
   *
   * @example
   * ```typescript
   * const [R, x] = PolynomialRingConstructor(F7, 'x');
   * const p = R.lagrange_polynomial([[F7(0), F7(1)], [F7(2), F7(2)], [F7(3), F7(6)]]);
   * // p(0) = 1, p(2) = 2, p(3) = 6
   * ```
   *
   * @see Reference: sage/rings/polynomial/polynomial_ring.py:2295 (lagrange_polynomial)
   */
  lagrange_polynomial(
    points: Array<[C, C]>,
    algorithm?: 'divided_difference' | 'pari'
  ): Polynomial<C>;
  lagrange_polynomial(
    points: Array<[C, C]>,
    algorithm: 'neville',
    previous_row?: Polynomial<C>[]
  ): Polynomial<C>[];
  lagrange_polynomial(
    points: Array<[C, C]>,
    algorithm: 'divided_difference' | 'neville' | 'pari' = 'divided_difference',
    previous_row?: Polynomial<C>[]
  ): Polynomial<C> | Polynomial<C>[] {
    const n = points.length;

    if (algorithm === 'divided_difference') {
      if (n === 0) {
        return this.zero();
      }
      return this._lagrange_divided_difference(points);
    }

    if (algorithm === 'neville') {
      return this._lagrange_neville(points, previous_row);
    }

    if (algorithm === 'pari') {
      // PARI's polinterpolate is Newton's divided-difference scheme, which is
      // what `_lagrange_divided_difference` implements.
      if (n === 0) {
        return this.zero();
      }
      return this._lagrange_divided_difference(points);
    }

    throw new ValueError("algorithm can be 'divided_difference', 'neville' or 'pari'");
  }

  /**
   * Compute Lagrange interpolation using divided differences.
   *
   * Uses Newton's form with divided differences, evaluated using
   * Horner's method for efficiency.
   */
  private _lagrange_divided_difference(points: Array<[C, C]>): Polynomial<C> {
    const n = points.length;
    if (n === 0) {
      return this.zero();
    }

    // Compute divided differences
    const F = this.divided_difference(points);

    // Evaluate in nested form (Horner's method)
    // P(x) = F[n-1] + (x - x_{n-2}) * (F[n-2] + (x - x_{n-3}) * (...))
    let P = this.__call__(F[n - 1]!);
    const x = this.gen();

    for (let i = n - 2; i >= 0; i--) {
      // P = P * (x - x_i) + F[i]
      const xi = this.__call__(points[i]![0]);
      P = P.mul(x.sub(xi));
      P = P.add(this.__call__(F[i]!));
    }

    return P;
  }

  /**
   * Return the Newton divided-difference coefficients.
   *
   * These are the coefficients F[0,0], F[1,1], ..., F[n-1,n-1] such that
   * P(x) = sum_{i=0}^{n-1} F[i,i] * prod_{j=0}^{i-1} (x - x_j)
   *
   * @param points - Array of [x, y] pairs
   * @param full_table - if true, return the full divided-difference table
   *   instead of only its main diagonal
   * @returns Array of divided difference coefficients (or the full table)
   *
   * @see Reference: sage/rings/polynomial/polynomial_ring.py:2206 (divided_difference)
   */
  divided_difference(points: Array<[C, C]>, full_table?: false): C[];
  divided_difference(points: Array<[C, C]>, full_table: true): C[][];
  divided_difference(points: Array<[C, C]>, full_table: boolean = false): C[] | C[][] {
    const n = points.length;
    if (n === 0) {
      return [];
    }

    // F[i][j] stores the divided difference f[x_{i-j}, ..., x_i]
    const F: C[][] = [];
    for (let i = 0; i < n; i++) {
      F.push([points[i]![1]]);
    }

    for (let i = 1; i < n; i++) {
      for (let j = 1; j <= i; j++) {
        const numer = F[i]![j - 1]!.sub(F[i - 1]![j - 1]!) as C;
        const denom = points[i]![0].sub(points[i - j]![0]) as C;
        const quotient = divideElements(numer, denom);
        F[i]!.push(quotient);
      }
    }

    if (full_table) {
      return F;
    }

    // Return diagonal elements F[i][i]
    return F.map((row, i) => row[i]!);
  }

  /**
   * Compute Lagrange interpolation using Neville's method.
   *
   * Returns the last row of the Neville table, where the last element
   * is the interpolating polynomial.
   *
   * @param points - Array of [x, y] pairs
   * @param previousRow - Optional previous row for incremental computation
   * @returns Array of polynomials (last row of Neville table)
   */
  private _lagrange_neville(points: Array<[C, C]>, previousRow?: Polynomial<C>[]): Polynomial<C>[] {
    const N = points.length;
    const M = previousRow?.length ?? 0;

    // P keeps track of the previous row, Q keeps track of the current row
    const P: Polynomial<C>[] = previousRow
      ? [...previousRow, ...new Array(N - M).fill(null)]
      : new Array(N).fill(null);
    const Q: Polynomial<C>[] = new Array(N).fill(null);

    const x = this.gen();

    for (let i = M; i < N; i++) {
      // Start populating the current row
      Q[0] = this.__call__(points[i]![1]);

      for (let j = 1; j <= i; j++) {
        // Q[j] = ((x - x_{i-j}) * Q[j-1] - (x - x_i) * P[j-1]) / (x_i - x_{i-j})
        const xiMinusJ = this.__call__(points[i - j]![0]);
        const xi = this.__call__(points[i]![0]);

        const numer = x
          .sub(xiMinusJ)
          .mul(Q[j - 1]!)
          .sub(x.sub(xi).mul(P[j - 1]!));

        const denom = points[i]![0].sub(points[i - j]![0]) as C;

        // Divide polynomial by scalar
        Q[j] = scalarDividePolynomial(numer, denom);
      }

      // Swap P and Q for next iteration
      for (let k = 0; k <= i; k++) {
        P[k] = Q[k]!;
      }
    }

    return P.filter((p) => p !== null) as Polynomial<C>[];
  }

  /**
   * Newton interpolation polynomial through the given points.
   *
   * This is more efficient than Lagrange when adding new points incrementally.
   * Uses divided differences internally.
   *
   * @param points - Array of [x, y] pairs where x values must be distinct
   * @returns The interpolating polynomial in Newton form
   *
   * @example
   * ```typescript
   * const [R, x] = PolynomialRingConstructor(F7, 'x');
   * const p = R.newton_interpolation([[F7(0), F7(1)], [F7(1), F7(3)], [F7(2), F7(7)]]);
   * ```
   *
   * @see Reference: sage/rings/polynomial/polynomial_ring.py:divided_difference
   */
  newton_interpolation(points: Array<[C, C]>): Polynomial<C> {
    // Newton interpolation is the same as Lagrange using divided differences
    // The difference is mainly conceptual - Newton form is easier to extend
    return this.lagrange_polynomial(points, 'divided_difference');
  }

  /**
   * Return the vanishing polynomial for a given domain.
   *
   * The vanishing polynomial is Z_H(X) = prod_{h in H} (X - h),
   * which evaluates to 0 for all elements in the domain H.
   *
   * @param domain - Array of field elements
   * @returns The vanishing polynomial
   *
   * @example
   * ```typescript
   * const [R, x] = PolynomialRingConstructor(F7, 'x');
   * const Z = R.vanishing_polynomial([F7(1), F7(2), F7(3)]);
   * // Z = (x - 1)(x - 2)(x - 3)
   * // Z(1) = Z(2) = Z(3) = 0
   * ```
   */
  vanishing_polynomial(domain: C[]): Polynomial<C> {
    if (domain.length === 0) {
      return this.one();
    }

    // Compute product (X - h) for all h in domain
    const x = this.gen();
    let result = this.one();

    for (const h of domain) {
      const factor = x.sub(this.__call__(h));
      result = result.mul(factor);
    }

    return result;
  }

  /**
   * Evaluate the interpolating polynomial at a point using barycentric interpolation.
   *
   * This is more efficient than constructing the full polynomial when you only
   * need the value at a single point. Complexity is O(n) for evaluation after
   * O(n) preprocessing, vs O(n^2) for full polynomial construction.
   *
   * Uses the second form of the barycentric interpolation formula:
   * L(x) = (sum_{j=0}^{n-1} w_j * y_j / (x - x_j)) / (sum_{j=0}^{n-1} w_j / (x - x_j))
   *
   * where w_j = 1 / prod_{k != j} (x_j - x_k) are the barycentric weights.
   *
   * @param points - Array of [x, y] pairs where x values must be distinct
   * @param evalPoint - The point at which to evaluate
   * @returns The value of the interpolating polynomial at evalPoint
   *
   * @throws {ValueError} If evalPoint coincides with one of the x values
   *
   * @example
   * ```typescript
   * const [R, x] = PolynomialRingConstructor(F7, 'x');
   * const val = R.barycentric_interpolation([[F7(0), F7(1)], [F7(1), F7(2)]], F7(3));
   * ```
   */
  barycentric_interpolation(points: Array<[C, C]>, evalPoint: C): C {
    const n = points.length;

    if (n === 0) {
      return this.base_ring.zero() as C;
    }

    // Check if evalPoint coincides with any x value
    for (let j = 0; j < n; j++) {
      if (evalPoint.eq(points[j]![0])) {
        // Return the corresponding y value directly
        return points[j]![1];
      }
    }

    // Compute barycentric weights: w_j = 1 / prod_{k != j} (x_j - x_k)
    const weights: C[] = [];
    for (let j = 0; j < n; j++) {
      let w = this.base_ring.one() as C;
      for (let k = 0; k < n; k++) {
        if (k !== j) {
          const diff = points[j]![0].sub(points[k]![0]) as C;
          w = divideElements(w, diff);
        }
      }
      weights.push(w);
    }

    // Compute L(x) using barycentric formula
    let numerator = this.base_ring.zero() as C;
    let denominator = this.base_ring.zero() as C;

    for (let j = 0; j < n; j++) {
      const diff = evalPoint.sub(points[j]![0]) as C;
      const term = divideElements(weights[j]!, diff);

      numerator = numerator.add(term.mul(points[j]![1]) as C) as C;
      denominator = denominator.add(term) as C;
    }

    return divideElements(numerator, denominator);
  }

  /**
   * Create a polynomial with the given roots.
   *
   * Returns the monic polynomial prod_{r in roots} (X - r).
   *
   * @param roots - Array of roots
   * @returns The polynomial with the given roots
   *
   * @example
   * ```typescript
   * const [R, x] = PolynomialRingConstructor(F7, 'x');
   * const p = R.from_roots([F7(1), F7(2), F7(3)]);
   * // p = (x - 1)(x - 2)(x - 3) = x^3 - 6x^2 + 11x - 6
   * // p(1) = p(2) = p(3) = 0
   * ```
   */
  from_roots(roots: C[]): Polynomial<C> {
    // This is the same as the vanishing polynomial
    return this.vanishing_polynomial(roots);
  }

  /**
   * Return the n-th cyclotomic polynomial.
   *
   * The n-th cyclotomic polynomial Phi_n(x) is the minimal polynomial over Q
   * of a primitive n-th root of unity. Its roots are exactly the primitive
   * n-th roots of unity.
   *
   * The polynomial is computed using the formula:
   * Phi_n(x) = prod_{d|n} (x^d - 1)^{mu(n/d)}
   *
   * where mu is the Moebius function.
   *
   * @param n - A positive integer
   * @returns The n-th cyclotomic polynomial
   *
   * @throws {ValueError} If n is not a positive integer
   *
   * @example
   * ```typescript
   * const [R, x] = PolynomialRingConstructor(F7, 'x');
   * const phi6 = R.cyclotomic_polynomial(6);
   * // Phi_6(x) = x^2 - x + 1
   * ```
   *
   * @see Reference: sage/rings/polynomial/cyclotomic.pyx:cyclotomic_coeffs
   */
  cyclotomic_polynomial(n: number): Polynomial<C> {
    if (n <= 0) {
      // sage/rings/polynomial/polynomial_ring.py:1188
      throw new ArithmeticError(`n=${n} must be positive`);
    }
    if (!Number.isInteger(n)) {
      throw new TypeError(`n=${n} must be an integer`);
    }

    // Compute the cyclotomic polynomial coefficients
    const coeffs = cyclotomicCoeffs(n);

    // Convert integer coefficients to ring elements
    const ringCoeffs: C[] = coeffs.map((c) => this.base_ring.__call__(c) as C);

    return new Polynomial(ringCoeffs, this);
  }

  toString(): string {
    const base = this.base_ring as CoefficientRing<C> & {
      characteristic?: bigint | (() => bigint);
    };
    const p =
      typeof base.characteristic === 'function' ? base.characteristic() : base.characteristic;
    const zero = base.zero();
    const prime = typeof zero === 'object' && zero !== null && 'value' in zero;
    const suffix =
      prime && p === 2n
        ? ' (using GF2X)'
        : prime && p !== undefined && p >= 2n ** 63n
          ? ' (using NTL)'
          : '';
    return `Univariate Polynomial Ring in ${this.variable_name} over ${this.base_ring}${suffix}`;
  }
}

/**
 * Divide two field elements.
 */
function divideElements<C extends RingElement>(a: C, b: C): C {
  // Try to call div method if it exists
  if ('div' in a && typeof (a as unknown as { div: (b: C) => C }).div === 'function') {
    return (a as unknown as { div: (b: C) => C }).div(b);
  }

  // Try inv method for field elements
  if ('inv' in b && typeof (b as unknown as { inv: () => C }).inv === 'function') {
    const bInv = (b as unknown as { inv: () => C }).inv();
    return a.mul(bInv) as C;
  }

  throw new ValueError('coefficient ring does not support division');
}

/**
 * Divide a polynomial by a scalar.
 */
function scalarDividePolynomial<C extends RingElement>(p: Polynomial<C>, c: C): Polynomial<C> {
  if (c.isZero()) {
    throw new ValueError('division by zero');
  }

  const cInv = divideElements(p.parent.base_ring.one() as C, c);
  return p.scalar_mul(cInv);
}

/**
 * Compute the coefficients of the n-th cyclotomic polynomial.
 *
 * Uses the relation: x^n - 1 = prod_{d|n} Phi_d(x)
 * So: Phi_n(x) = (x^n - 1) / prod_{d|n, d<n} Phi_d(x)
 *
 * For prime powers p^k: Phi_{p^k}(x) = Phi_p(x^{p^{k-1}})
 * For non-squarefree n = m * rad where rad is the radical:
 *   Phi_n(x) = Phi_rad(x^{n/rad})
 *
 * @param n - A positive integer
 * @returns Array of integer coefficients [c_0, c_1, ..., c_deg]
 */
function cyclotomicCoeffs(n: number): number[] {
  if (n === 1) {
    return [-1, 1]; // x - 1
  }

  // Factor n
  const factorization = primeFactorization(n);

  // Check for non-squarefree case first
  const isSquarefree = factorization.every(([_, exp]) => exp === 1);

  if (!isSquarefree) {
    // For non-squarefree n: Phi_n(x) = Phi_rad(x^{n/rad})
    // where rad is the radical (product of distinct prime factors)
    const rad = factorization.map(([p, _]) => p).reduce((a, b) => a * b, 1);
    const pow = n / rad;
    const phiRad = cyclotomicCoeffs(rad);

    // Substitute x^pow for x
    const degree = (phiRad.length - 1) * pow;
    const result = new Array(degree + 1).fill(0);
    for (let i = 0; i < phiRad.length; i++) {
      result[i * pow] = phiRad[i]!;
    }
    return result;
  }

  // Now n is squarefree
  // For prime p: Phi_p(x) = 1 + x + x^2 + ... + x^{p-1}
  if (factorization.length === 1) {
    const p = factorization[0]![0];
    return new Array(p).fill(1);
  }

  // For squarefree n with multiple prime factors, use recursion:
  // Phi_n(x) = (x^n - 1) / prod_{d|n, d<n} Phi_d(x)

  // Get all proper divisors of n (divisors < n)
  const divisors = getDivisors(n).filter((d) => d < n);

  // Start with x^n - 1
  let result = new Array(n + 1).fill(0);
  result[0] = -1;
  result[n] = 1;

  // Divide by Phi_d for each proper divisor d
  for (const d of divisors) {
    const phiD = cyclotomicCoeffs(d);
    result = polynomialDivide(result, phiD);
  }

  return result;
}

/**
 * Get all divisors of n.
 */
function getDivisors(n: number): number[] {
  const divisors: number[] = [];
  for (let i = 1; i * i <= n; i++) {
    if (n % i === 0) {
      divisors.push(i);
      if (i !== n / i) {
        divisors.push(n / i);
      }
    }
  }
  return divisors.sort((a, b) => a - b);
}

/**
 * Divide polynomial a by polynomial b (exact division).
 * Assumes b divides a exactly with no remainder.
 */
function polynomialDivide(a: number[], b: number[]): number[] {
  // Remove trailing zeros
  while (a.length > 1 && a[a.length - 1] === 0) {
    a.pop();
  }
  while (b.length > 1 && b[b.length - 1] === 0) {
    b.pop();
  }

  if (b.length === 0 || (b.length === 1 && b[0] === 0)) {
    throw new Error('Division by zero polynomial');
  }

  const degA = a.length - 1;
  const degB = b.length - 1;

  if (degA < degB) {
    return [0];
  }

  const degQ = degA - degB;
  const quotient = new Array(degQ + 1).fill(0);
  const remainder = [...a];

  const lcB = b[degB]!;

  for (let i = degQ; i >= 0; i--) {
    const coeff = remainder[i + degB]! / lcB;
    quotient[i] = coeff;

    for (let j = 0; j <= degB; j++) {
      remainder[i + j] -= coeff * b[j]!;
    }
  }

  // Remove trailing zeros from quotient
  while (quotient.length > 1 && quotient[quotient.length - 1] === 0) {
    quotient.pop();
  }

  return quotient;
}

/**
 * Get the prime factorization of n as [prime, exponent] pairs.
 *
 * Uses PARI's factorization via the factor() function.
 */
function primeFactorization(n: number): Array<[number, number]> {
  if (n <= 1) {
    return [];
  }

  // Use PARI's factorization and convert to number pairs
  const factorization = factor(BigInt(n));
  return factorization
    .filter(([p, _]) => p > 0n) // Exclude -1 sign factor
    .map(([p, e]) => [Number(p), Number(e)]);
}

/**
 * Create a polynomial ring over a base ring.
 *
 * @param base_ring - The coefficient ring
 * @param names - Variable name(s)
 * @returns A polynomial ring and its generator
 *
 * @example
 * ```typescript
 * const [R, x] = PolynomialRingConstructor(GF2, 'x');
 * const p = x.pow(2).add(x).add(R.one());  // x^2 + x + 1
 * ```
 */
export function PolynomialRingConstructor<C extends RingElement>(
  base_ring: CoefficientRing<C>,
  names: string = 'x'
): [PolynomialRing<C>, Polynomial<C>] {
  const ring = new PolynomialRing(base_ring, names);
  return [ring, ring.gen()];
}
