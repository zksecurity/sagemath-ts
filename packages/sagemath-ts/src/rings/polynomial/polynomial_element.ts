import { cmp_universal as pariFiniteCompare } from '@sagemath-ts/parigp-ts/src/gen2.js';
import { PariType } from '@sagemath-ts/parigp-ts/src/types.js';
import { factor as ntlIntegerFactor } from '@sagemath-ts/ntl-ts/src/ZZXFactoring.js';
import {
  ZX_factor as pariIntegerFactor,
  QX_factor as pariRationalFactor,
} from '@sagemath-ts/parigp-ts/src/QX_factor.js';
import { _squarefree_decomposition_univariate_polynomial as fieldSquarefree } from '../../categories/fields.js';
import { _squarefree_decomposition_univariate_polynomial as finiteSquarefree } from '../finite_rings/finite_field_base.js';
import { CompiledPolynomialFunction } from './polynomial_compiled.js';
import { _coefficient_nth_root } from '../finite_rings/element_base.js';
/**
 * @module sage/rings/polynomial/polynomial_element
 * @description Polynomial elements over arbitrary coefficient rings
 *
 * Port of: sage/rings/polynomial/polynomial_element.pyx
 */

import {
  fmpq_poly_get_numerator,
  fmpq_poly_get_denominator,
  _fmpq_poly_lcm,
  nmod_poly_factor_squarefree,
  nmod_poly_factor,
  _fmpz_poly_evaluate_fmpz,
  _fmpq_poly_evaluate_fmpz,
  _fmpq_poly_evaluate_fmpq,
  _fmpq_poly_compose,
  _nmod_poly_evaluate_nmod,
  _nmod_poly_compose,
  _fmpz_poly_mullow,
  _fmpq_poly_mullow,
  _fmpq_poly_mul,
  _fmpz_poly_mul,
  _nmod_poly_mul,
  _nmod_poly_inv_series_newton,
  _nmod_poly_powmod_ui_binexp,
  _nmod_poly_powmod_fmpz_binexp_preinv,
  _nmod_poly_powmod_x_fmpz_preinv,
  _nmod_poly_mullow,
  _fmpz_poly_pow_trunc,
  _nmod_poly_pow_trunc,
  _fmpz_poly_inv_series,
  _fmpq_poly_inv_series_newton,
  _fmpz_poly_pow,
  _fmpq_poly_pow,
  _nmod_poly_pow,
  _fmpz_poly_gcd,
  _fmpq_poly_gcd,
  _nmod_poly_gcd,
  _nmod_poly_xgcd,
  _fmpz_poly_xgcd,
  _fmpq_poly_xgcd,
  _fmpz_poly_resultant,
  _fmpq_poly_resultant,
  _fmpq_poly_derivative,
  _nmod_poly_resultant,
  _nmod_poly_make_monic,
  _fmpz_poly_divrem,
  _nmod_poly_divrem,
} from '@sagemath-ts/flint-ts';
import { Integer, ZZ } from '../integer_ring.js';
type IntegerInput = Parameters<typeof ZZ.__call__>[0];
import {
  ZZX_SquareFreeDecomp,
  _ZZX_kernels,
  GF2X,
  ZZ_pX_mul,
  ZZ_pEX_mul,
  ZZ_pEX_eval,
  ZZ_pX_evaluate,
  ZZ_pX_power,
  ZZ_pEX_power,
  ZZ_pEX_InvTrunc,
  ZZ_pEX_PowerMod,
  ZZ_pEX_PowerXMod,
} from '@sagemath-ts/ntl-ts';
import { PariError, F2x_factor, FpX_factor as pariFpXFactor, resultant as pariResultant } from '@sagemath-ts/parigp-ts';
import { FractionField, FractionField_generic } from '../fraction_field.js';
import { isFractionElement, type FractionElement } from '../fraction_field_element.js';
import { RDF, RealDoubleElement } from '../real_double.js';
import { FiniteFieldElement, PrimeField } from '../finite_rings/finite_field_extension.js';
import { Rational } from '../rational.js';
import { GF2, GF2Element } from '../finite_rings/gf2.js';
import { canonicalFiniteOperands, IntegerMod } from '../finite_rings/integer_mod.js';
import { IntegerModRing } from '../finite_rings/integer_mod_ring.js';
import { FiniteFieldElement as LegacyPrimeElement } from '../finite_rings/finite_field_prime.js';
import { PrimeFieldElement } from '../finite_rings/finite_field_extension.js';
import { QQ } from '../rational_field.js';
import { PolynomialRing } from './polynomial_ring.js';
import { Matrix_modn_dense } from '../../matrix/matrix_modn.js';
import { Matrix_mod2_dense } from '../../matrix/matrix_mod2.js';
import { Matrix as EvaluationGenericMatrix } from '../../matrix/matrix_generic.js';
import { integer_to_real_double_dense } from '../../matrix/change_ring.js';
import { IntegerMatrix as EvaluationIntegerMatrix } from '../../matrix/matrix_integer.js';
import { Zmod } from '../finite_rings/integer_mod_ring.js';
type EvaluationMatrix =
  | EvaluationGenericMatrix<RingElement>
  | EvaluationIntegerMatrix
  | Matrix_modn_dense
  | Matrix_mod2_dense;
import { multi_derivative } from '../../misc/derivative.js';
import {
  factor as factorInteger,
  gcd as gcdBigInt,
  is_prime,
  next_prime,
} from '../../arith/misc.js';
import {
  ArithmeticError,
  AttributeError,
  AssertionError,
  IndexError,
  NotImplementedError,
  NTLError,
  RuntimeError,
  OverflowError,
  ValueError,
  ZeroDivisionError,
} from '../../errors.js';
import { current_randstate } from '../../misc/randstate.js';

/**
 * Interface for coefficient rings/fields.
 */
export interface CoefficientRing<T> {
  zero(): T;
  one(): T;
  __call__(x: unknown): T;
  is_field?(): boolean;
}

/**
 * Interface for ring elements that can be used as polynomial coefficients.
 */
export interface RingElement {
  add(other: this): this;
  sub(other: this): this;
  mul(other: this): this;
  neg(): this;
  eq(other: this | number): boolean;
  isZero(): boolean;
  toString(): string;
}

/**
 * A polynomial with coefficients in a ring R.
 *
 * Internally stored as an array of coefficients where coeffs[i] is the
 * coefficient of x^i. Trailing zeros are removed.
 */
type NonnegativePowerLiteral<N extends number | bigint> = number extends N
  ? never
  : bigint extends N
    ? never
    : Extract<`${N}`, `-${string}`> extends never
      ? N
      : never;

export class Polynomial<C extends RingElement> {
  /** Generic scalar hooks; parent constructors may choose an earlier partial section.
   * @see Deviation: Native Modular Polynomial Products and Fraction Fields
   */
  _scalar_conversion<T>(R: { __call__(x: unknown): T }): T {
    if (this.degree() > 0) throw new TypeError('cannot convert nonconstant polynomial');
    return R.__call__(this.getCoeff(0));
  }
  _integer_(R: typeof ZZ): bigint {
    return this._scalar_conversion(R);
  }
  _rational_(): Rational {
    return this._scalar_conversion(QQ);
  }

  readonly coeffs: readonly C[];
  readonly parent: PolynomialRingBase<C>;

  constructor(coeffs: C[], parent: PolynomialRingBase<C>, is_gen = false) {
    if (is_gen) {
      polynomialGenerators.add(this);
      const backend = polynomialBackend(parent.base_ring);
      if (backend === 'rational' || backend === 'binary')
        coeffs = [parent.base_ring.zero(), parent.base_ring.one()];
    }
    this.parent = parent;

    // Remove trailing zeros
    let len = coeffs.length;
    while (len > 0 && coeffs[len - 1]!.isZero()) {
      len--;
    }
    this.coeffs = coeffs.slice(0, len);
  }

  /** Original distinguished-generator predicate, including native equality overrides. */
  is_gen(): boolean | bigint {
    const backend = polynomialBackend(this.parent.base_ring);
    if (
      backend === 'word' ||
      backend === 'binary' ||
      backend === 'extension' ||
      (backend === 'large' && this.parent.base_ring.is_field?.() !== true)
    ) {
      const result =
        this.degree() === 1 &&
        this.getCoeff(0).isZero() &&
        this.getCoeff(1).eq(this.parent.base_ring.one());
      return backend === 'large' ? BigInt(result) : result;
    }
    return polynomialGenerators.has(this);
  }

  /** Generic _new_c deliberately preserves unchecked coefficient storage. */
  private _new_c(coeffs: C[], parent: PolynomialRingBase<C>): Polynomial<C> {
    return Object.create(Polynomial.prototype, {
      coeffs: { value: coeffs.slice(), enumerable: true },
      parent: { value: parent, enumerable: true },
    }) as Polynomial<C>;
  }

  /**
   * Return the degree of this polynomial.
   * The zero polynomial has degree -1.
   */
  degree(): number {
    return this.coeffs.length - 1;
  }

  /**
   * Return the leading coefficient.
   */
  leading_coefficient(): C {
    if (this.coeffs.length === 0) {
      return this.parent.base_ring.zero() as C;
    }
    return this.coeffs[this.coeffs.length - 1]!;
  }

  /**
   * Return the coefficient of x^n.
   */
  getCoeff(n: unknown): C {
    const index = polynomialInteger(n, 'index');
    if (index < 0n || index >= BigInt(this.coeffs.length)) return this.parent.base_ring.zero();
    return this.coeffs[Number(index)]!;
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
   * Check if this is monic (leading coefficient is 1).
   */
  is_monic(): boolean {
    return !this.isZero() && this.leading_coefficient().eq(this.parent.base_ring.one());
  }

  private _hasCompatibleParent(other: Polynomial<C>): boolean {
    if (!(other instanceof Polynomial)) return false;
    if (this.parent === other.parent) return true;
    return polynomialCommonBase(this.parent, other.parent) !== null;
  }

  private _checkVariable(other: Polynomial<C>, operation: string): void {
    if (!this._hasCompatibleParent(other)) {
      throw new TypeError(
        `unsupported operand parent(s) for ${operation}: '${this.parent}' and '${other.parent}'`
      );
    }
  }

  /**
   * Add two polynomials.
   * @see Deviation: Polynomial Common Coefficient Parents and Representation
   */
  add(other: Polynomial<C>): Polynomial<C>;
  add<D extends RingElement>(other: Polynomial<D>): Polynomial<C | D>;
  add(other: Polynomial<RingElement>): Polynomial<RingElement> {
    if (this.parent !== other.parent) {
      const [a, b] = polynomialCommonOperands(this, other, '+');
      return a.add(b);
    }
    const maxLen = Math.max(this.coeffs.length, other.coeffs.length);
    const result: C[] = [];

    for (let i = 0; i < maxLen; i++) {
      const a = this.getCoeff(i);
      const b = other.getCoeff(i);
      result.push(a.add(b as C) as C);
    }

    return new Polynomial(result, this.parent);
  }

  /**
   * Subtract two polynomials.
   * @see Deviation: Polynomial Common Coefficient Parents and Representation
   */
  sub(other: Polynomial<C>): Polynomial<C>;
  sub<D extends RingElement>(other: Polynomial<D>): Polynomial<C | D>;
  sub(other: Polynomial<RingElement>): Polynomial<RingElement> {
    if (this.parent !== other.parent) {
      const [a, b] = polynomialCommonOperands(this, other, '-');
      return a.sub(b);
    }
    const maxLen = Math.max(this.coeffs.length, other.coeffs.length);
    const result: C[] = [];

    for (let i = 0; i < maxLen; i++) {
      const a = this.getCoeff(i);
      const b = other.getCoeff(i);
      result.push(a.sub(b as C) as C);
    }

    return new Polynomial(result, this.parent);
  }

  /**
   * Negate the polynomial.
   */
  neg(): Polynomial<C> {
    return new Polynomial(
      this.coeffs.map((c) => c.neg() as C),
      this.parent
    );
  }

  /**
   * Multiply two polynomials. @see Deviation: Polynomial Integer Powers and Portable Native Products
   * @see Deviation: Truncated Multiplication Parents and Real Coefficients
   * @see Deviation: Polynomial Common Coefficient Parents and Representation
   */
  mul(other: Polynomial<C>): Polynomial<C>;
  mul<D extends RingElement>(other: Polynomial<D>): Polynomial<C | D>;
  mul(other: Polynomial<RingElement>): Polynomial<RingElement> {
    if (this.parent !== other.parent) {
      const [a, b] = polynomialCommonOperands(this, other, '*');
      // Generic-dense scalar actions preserve an already-zero outer polynomial.
      if (
        this.parent.variable_name !== other.parent.variable_name &&
        polynomialBackend(a.parent.base_ring) === 'generic'
      ) {
        const outer = a.parent.variable_name === this.parent.variable_name ? a : b;
        if (outer.isZero()) return outer;
      }
      return a.mul(b);
    }
    if (this.isZero() || other.isZero()) {
      return this.parent.zero();
    }

    const base = this.parent.base_ring,
      backend = polynomialBackend(base);
    if (backend === 'integer' || backend === 'word' || backend === 'large') {
      const a = extractIntegerCoeffs(this),
        b = this === other ? a : extractIntegerCoeffs(other);
      const out =
        backend === 'integer'
          ? _fmpz_poly_mul(a, b)
          : backend === 'word'
            ? _nmod_poly_mul(a, b, getRingCharacteristic(base)!)
            : ZZ_pX_mul(a, b, getRingCharacteristic(base)!);
      return new Polynomial(
        out.map((c) => base.__call__(c)),
        this.parent
      );
    }
    if (backend === 'rational') {
      const [a, da] = polynomialRationalData(this),
        [b, db] = this === other ? [a, da] : polynomialRationalData(other);
      const [out, den] = _fmpq_poly_mul(a, da, b, db),
        D = base.__call__(den);
      return new Polynomial(
        out.map((c) => divideCoeffs(base.__call__(c), D)),
        this.parent
      );
    }
    if (backend === 'binary') {
      const pack = (v: Polynomial<RingElement>) =>
        new GF2X(v.coeffs.map((c) => Number((c as unknown as { value: bigint | number }).value)));
      const a = pack(this),
        b = this === other ? a : pack(other),
        rep = a.mul(b).rep();
      const out =
        rep === 0n
          ? []
          : rep
              .toString(2)
              .split('')
              .reverse()
              .map((c) => base.__call__(BigInt(c)));
      return new Polynomial(out, this.parent);
    }
    if (backend === 'extension') {
      const B = base as unknown as FiniteFieldElement['parent'],
        f = B.modulus.coeffs.map((c) => c.value);
      const pack = (v: Polynomial<RingElement>) =>
        v.coeffs.map((c) => (c as unknown as FiniteFieldElement).lift.coeffs.map((v) => v.value));
      const a = pack(this),
        b = this === other ? a : pack(other);
      return new Polynomial(
        ZZ_pEX_mul(a, b, f, B.characteristic).map((c) => B.__call__(c) as unknown as C),
        this.parent
      );
    }
    {
      // Generic multiplication starts with ordered term actions; exact rings
      // use Karatsuba and inexact rings retain the specialized square order.
      const leftTerm = this.coeffs.filter((c) => !c.isZero()).length === 1;
      const rightTerm = other.coeffs.filter((c) => !c.isZero()).length === 1;
      if (leftTerm || rightTerm) {
        const term = leftTerm ? this : other,
          poly = leftTerm ? other : this,
          c = term.leading_coefficient();
        const out = poly.coeffs.map((a) => (leftTerm ? c.mul(a) : a.mul(c)));
        return new Polynomial(
          Array<RingElement>(term.degree()).fill(this.parent.base_ring.zero()).concat(out),
          this.parent
        );
      }
      if (!polynomialRdfBase(base)) {
        const threshold =
          base instanceof PolynomialRing
            ? 0
            : base instanceof FractionField_generic
              ? Number.MAX_SAFE_INTEGER
              : 8;
        return new Polynomial(
          do_karatsuba_product(this.coeffs, other.coeffs, threshold),
          this.parent
        );
      }
      if (this === other) {
        const out = Array<C>(2 * this.degree() + 1).fill(this.parent.base_ring.zero()),
          two = this.parent.base_ring.__call__(2n);
        for (let i = 0; i < this.coeffs.length; i++) {
          out[2 * i] = this.coeffs[i]!.mul(this.coeffs[i]!);
          for (let j = 0; j < i; j++)
            out[i + j] = out[i + j]!.add(two.mul(this.coeffs[i]!).mul(this.coeffs[j]!));
        }
        return new Polynomial(out, this.parent);
      }
      return new Polynomial(do_schoolbook_product(this.coeffs, other.coeffs, -1), this.parent);
    }
  }

  /**
   * Multiply by a scalar.
   */
  scalar_mul(c: C): Polynomial<C> {
    if (c.isZero()) {
      return this.parent.zero();
    }
    return new Polynomial(
      this.coeffs.map((coeff) => coeff.mul(c) as C),
      this.parent
    );
  }

  /** Native LCM: ZZ retains its signed product; fields return a monic result.
   * @see Reference: polynomial_element.pyx:5533; polynomial_integer_dense_flint.pyx:835;
   * polynomial_rational_flint.pyx:1003
   */
  lcm(other: Polynomial<C>): Polynomial<C> {
    if (this.isZero() || other.isZero()) return this.parent.zero();
    const backend = polynomialBackend(this.parent.base_ring);
    if (backend === 'rational') {
      const [a] = polynomialRationalData(this), [b] = polynomialRationalData(other);
      const [n, d] = _fmpq_poly_lcm(a, b), k = this.parent.base_ring;
      return this.parent.__call__(n.map((c) => divideCoeffs(k.__call__(c), k.__call__(d))));
    }
    const g = this.gcd(other);
    if (backend === 'integer') return this.quo_rem(g)[0].mul(other);
    return this.mul(other).quo_rem(g)[0].monic();
  }

  /** Native coefficient denominator (QQ uses FLINT's stored positive integer).
   * @see Reference: polynomial_element.pyx:4026; polynomial_rational_flint.pyx:1493
   */
  denominator(): bigint | RingElement {
    if (polynomialBackend(this.parent.base_ring) === 'rational') {
      const [a, den] = polynomialRationalData(this);
      return fmpq_poly_get_denominator(a, den);
    }
    if (this.isZero()) return this.parent.base_ring.one();
    let d: unknown;
    for (const c of this.coeffs.filter((c) => !c.isZero())) {
      const method = (c as unknown as { denominator?: unknown }).denominator;
      if (typeof method !== 'function') return this.parent.base_ring.one();
      const next = method.call(c);
      if (d === undefined) d = next;
      const raw = (x: unknown) => typeof x === 'bigint' ? x : x instanceof Integer ? x.value : null;
      const a = raw(d), b = raw(next);
      if (a !== null && b !== null) d = a / gcdBigInt(a, b) * b;
      else if (d && typeof (d as { lcm?: unknown }).lcm === 'function')
        d = (d as { lcm(x: unknown): unknown }).lcm(next);
      else return this.parent.base_ring.one();
    }
    return d as bigint | RingElement;
  }

  /** QQ's numerator belongs to ZZ[x]; other backends multiply by the coefficient denominator.
   * @see Reference: polynomial_element.pyx:4115; polynomial_rational_flint.pyx:1468
   */
  numerator(): Polynomial<C> | Polynomial<RingElement> {
    if (polynomialBackend(this.parent.base_ring) === 'rational') {
      const [a, den] = polynomialRationalData(this);
      let R = polynomialIntegerNumeratorParents.get(this.parent.variable_name);
      if (!R) {
        R = new PolynomialRing(polynomialIntegerCoefficientRing, this.parent.variable_name);
        polynomialIntegerNumeratorParents.set(this.parent.variable_name, R);
      }
      return R.__call__(fmpq_poly_get_numerator(a, den));
    }
    return this.scalar_mul(this.parent.base_ring.__call__(this.denominator()));
  }

  /**
   * Compute the remainder of this polynomial divided by other.
   * Only works over fields.
   */
  mod(other: Polynomial<C>): Polynomial<C>;
  mod<D extends RingElement>(other: Polynomial<D>): Polynomial<C | D>;
  mod(operand: Polynomial<RingElement>): Polynomial<RingElement> {
    const other = operand as Polynomial<C>;
    const [_q, r] = this.quo_rem(other);
    return r;
  }

  /**
   * Compute quotient and remainder of the Euclidean division.
   *
   * Uses the backend's zero-divisor error. Generic coefficient division raises
   * ArithmeticError if a quotient coefficient does not lie in the base ring;
   * native ZZ division instead keeps a possibly high-degree remainder.
   * @see Deviation: Polynomial Quotient and Remainder Backends
   *
   * @see Reference: sage/rings/polynomial/polynomial_element.pyx:12548 (quo_rem)
   */
  quo_rem(other: Polynomial<C>): [Polynomial<C>, Polynomial<C>];
  quo_rem<D extends RingElement>(other: Polynomial<D>): [Polynomial<C | D>, Polynomial<C | D>];
  quo_rem(operand: Polynomial<RingElement>): [Polynomial<RingElement>, Polynomial<RingElement>] {
    const other = operand as Polynomial<C>;
    if (this.parent !== other.parent) {
      const [a, b] = polynomialCommonOperands(this, other, 'quo_rem');
      return a.quo_rem(b);
    }
    const backend = polynomialBackend(this.parent.base_ring);
    if (other.isZero()) {
      if (backend === 'large') throw new NTLError('ZZ_pX: division by zero');
      throw new ZeroDivisionError(
        backend === 'word' || backend === 'binary' || backend === 'extension'
          ? ''
          : 'division by zero polynomial'
      );
    }
    if (this.isZero() && (backend === 'integer' || backend === 'rational' || backend === 'generic'))
      return [this, this];
    if (backend === 'integer') {
      const [q, r] = _fmpz_poly_divrem(extractIntegerCoeffs(this), extractIntegerCoeffs(other));
      return [
        new Polynomial(
          q.map((c) => this.parent.base_ring.__call__(c)),
          this.parent
        ),
        new Polynomial(
          r.map((c) => this.parent.base_ring.__call__(c)),
          this.parent
        ),
      ];
    }
    if (backend === 'word') {
      const n = getRingCharacteristic(this.parent.base_ring)!;
      const lift = (c: C) => (c as unknown as { value: bigint }).value;
      if (gcdBigInt(lift(other.leading_coefficient()), n) !== 1n)
        throw new ValueError('Leading coefficient of a must be invertible.');
      const [q, r] = _nmod_poly_divrem(this.coeffs.map(lift), other.coeffs.map(lift), n);
      return [
        new Polynomial(
          q.map((c) => this.parent.base_ring.__call__(c)),
          this.parent
        ),
        new Polynomial(
          r.map((c) => this.parent.base_ring.__call__(c)),
          this.parent
        ),
      ];
    }
    if (backend === 'binary') {
      const packed = (f: Polynomial<C>) =>
        new GF2X(
          f.isZero()
            ? 0n
            : BigInt(
                '0b' +
                  f.coeffs
                    .map((c) => String((c as unknown as { value: bigint | number }).value))
                    .reverse()
                    .join('')
              )
        );
      return packed(this)
        .DivRem(packed(other))
        .map((f) => {
          const bits = f.rep() === 0n ? [] : f.rep().toString(2).split('').reverse();
          return new Polynomial(
            bits.map((c) => this.parent.base_ring.__call__(BigInt(c))),
            this.parent
          );
        }) as [Polynomial<C>, Polynomial<C>];
    }
    if (this.degree() < other.degree()) {
      return [
        this.parent.zero(),
        backend === 'generic' ? this : new Polynomial([...this.coeffs], this.parent),
      ];
    }

    // Make a mutable copy of coefficients
    const remainder = [...this.coeffs] as C[];
    const divisorLC = other.leading_coefficient();
    const divisorDeg = other.degree();
    const quotientCoeffs: C[] = [];

    // Sage first tries ``inverse_of_unit()`` on the leading coefficient; when
    // that succeeds every quotient coefficient is automatically in the base
    // ring and no further check is needed.  Only in the fallback branch
    // ("convert") does it verify that the quotient coefficient lies in R.
    const lcInverse = inverseOfUnit(divisorLC, this.parent.base_ring);

    // Initialize quotient with zeros
    for (let i = 0; i <= this.degree() - divisorDeg; i++) {
      quotientCoeffs.push(this.parent.base_ring.zero() as C);
    }

    for (let i = this.degree(); i >= divisorDeg; i--) {
      if (backend !== 'generic' && remainder[i]?.isZero()) {
        continue;
      }

      // Compute quotient coefficient
      // This requires the coefficient ring to support division
      let qCoeff: C;
      if (lcInverse !== null) {
        qCoeff = remainder[i]!.mul(lcInverse) as C;
      } else {
        qCoeff = divideCoeffs(remainder[i]!, divisorLC);
        // Sage raises here when the quotient does not lie in the base ring
        // (`polynomial_element.pyx:12634-12640`); a coefficient ring whose
        // division truncates (e.g. ZZ) would otherwise silently return garbage.
        if (!qCoeff.mul(divisorLC).eq(remainder[i]!)) {
          throw new ArithmeticError(
            'division non exact (consider coercing to polynomials over the fraction field)'
          );
        }
      }
      quotientCoeffs[i - divisorDeg] = qCoeff;

      // Subtract qCoeff * other * x^(i - divisorDeg) from remainder
      for (let j = divisorDeg - 1; j >= 0; j--) {
        const prod = qCoeff.mul(other.coeffs[j]!) as C;
        remainder[i - divisorDeg + j] = remainder[i - divisorDeg + j]!.sub(prod) as C;
      }
    }

    return [
      backend === 'generic'
        ? this._new_c(quotientCoeffs, this.parent)
        : new Polynomial(quotientCoeffs, this.parent),
      new Polynomial(remainder.slice(0, divisorDeg), this.parent),
    ];
  }

  /**
   * Compute the pseudo-division of two polynomials.
   *
   * Returns `[Q, R]` such that `l^(m-n+1) * self = Q*other + R` with
   * `deg(R) < deg(other)`, where `m = deg(self)`, `n = deg(other)` and `l` is
   * the leading coefficient of `other`.  Unlike {@link quo_rem} this needs no
   * division during the cancellation loop. Degree gaps can require negative
   * coefficient powers, and constant operands can produce a fraction-field quotient.
   *
   * Algorithm 3.1.2 in [Coh1993].
   *
   * @see Reference: sage/rings/polynomial/polynomial_element.pyx:5375 (pseudo_quo_rem)
   * @see Deviation: Polynomial Pseudo-Division and Fraction Parents
   */
  pseudo_quo_rem(other: unknown): [Polynomial<C> | FractionElement<C>, Polynomial<C>] {
    let B = other;
    if (typeof B === 'bigint' || (typeof B === 'number' && Number.isInteger(B))) B = new Integer(B);
    if (B === null || B === undefined || typeof (B as RingElement).isZero !== 'function')
      throw new AttributeError(`'${polynomialPseudoType(B)}' object has no attribute 'is_zero'`);
    if ((B as RingElement).isZero())
      throw new ZeroDivisionError('Pseudo-division by zero is not possible');
    if (polynomialCoefficientContains(this.parent.base_ring, B)) {
      if (B instanceof Polynomial) {
        if (this.degree() < 0) {
          const power = polynomialCoefficientPower(B, this.degree());
          if (isFractionElement(power)) {
            const original = power.parent.ring();
            const common = polynomialCommonBase(this.parent, original);
            if (!(common instanceof PolynomialRing) || !polynomialDomain(common))
              throw new TypeError(
                `unsupported operand parent(s) for *: '${this.parent}' and '${power.parent}'`
              );
            const field = FractionField(common);
            return [
              field.__call__(this, common.__call__(power.denominator())),
              this.parent.zero(),
            ] as [FractionElement<C>, Polynomial<C>];
          }
          const Q = this.mul(power as Polynomial<C>);
          return [
            Q,
            polynomialBackend(B.parent.base_ring) === 'generic' && Q.parent === this.parent
              ? Q
              : this.parent.zero(),
          ];
        }
        return [this.mul(B.pow(this.degree()) as Polynomial<C>), this.parent.zero()];
      }
      const scale = polynomialCoefficientPower(B as RingElement, this.degree());
      return [polynomialScalarProduct(this, scale, false) as Polynomial<C>, this.parent.zero()];
    }
    if (!(B instanceof Polynomial))
      throw new AttributeError(`'${polynomialPseudoType(B)}' object has no attribute 'degree'`);
    let R: Polynomial<RingElement> = this;
    let Q: Polynomial<RingElement> = this.parent.zero();
    let e = this.degree() - B.degree() + 1;
    const d = B.leading_coefficient();
    while (R.degree() >= B.degree()) {
      const c = R.leading_coefficient();
      const diffdeg = R.degree() - B.degree();
      Q = polynomialScalarProduct(Q, d, true).add(this.parent.__call__(c).shift(diffdeg));
      R = polynomialScalarProduct(R, d, true).sub(
        polynomialScalarProduct(B.shift(diffdeg), c, true)
      );
      e -= 1;
    }
    const q = polynomialCoefficientPower(d, e, B.parent.base_ring);
    return [polynomialScalarProduct(Q, q, true), polynomialScalarProduct(R, q, true)] as [
      Polynomial<C>,
      Polynomial<C>,
    ];
  }

  /**
   * Compute this^n, optionally reduced modulo a polynomial.
   * @see Deviation: Polynomial Modular Powers
   * Negative integer powers may return an element of the fraction field.
   * @see Deviation: Polynomial Integer Powers and Portable Native Products
   * @see Deviation: Polynomial Roots and Truncated Series
   */
  pow<N extends number | bigint>(
    n: N & NonnegativePowerLiteral<N>,
    modulus?: Polynomial<C> | null
  ): Polynomial<C>;
  pow(n: unknown, modulus?: unknown): Polynomial<C> | FractionElement<C>;
  pow(n: unknown, modulus?: unknown): Polynomial<C> | FractionElement<C> {
    const backend = polynomialBackend(this.parent.base_ring);
    if (
      modulus !== undefined &&
      modulus !== null &&
      (backend === 'integer' || backend === 'rational')
    )
      throw new NotImplementedError('pow() with a modulus is not implemented for this ring');
    let exponent: bigint;
    if (backend === 'integer' || backend === 'rational') {
      if (n instanceof Rational && n.denominator !== 1n) {
        if (this.degree() === 0) {
          const c = this.getCoeff(0) as unknown as {
            nth_root(e: bigint): { pow(e: bigint): unknown };
          };
          return this.parent.__call__(c.nth_root(n.denominator).pow(n.numerator));
        }
        return this.nth_root(n.denominator).pow(n.numerator);
      }
      if (
        n === null ||
        n === undefined ||
        typeof n === 'string' ||
        Array.isArray(n) ||
        n instanceof Polynomial ||
        n instanceof RealDoubleElement ||
        (typeof n === 'number' && !Number.isInteger(n))
      ) {
        const source =
          n instanceof Polynomial || n instanceof RealDoubleElement
            ? String(n.parent)
            : `<class '${n === null || n === undefined ? 'NoneType' : Array.isArray(n) ? 'list' : typeof n === 'string' ? 'str' : 'float'}'>`;
        throw new TypeError(`no canonical coercion from ${source} to Rational Field`);
      }
      exponent = ZZ.__call__(n as IntegerInput);
    } else if (
      modulus !== undefined &&
      modulus !== null &&
      (backend === 'word' || backend === 'extension')
    ) {
      exponent = ZZ.__call__(n as IntegerInput);
    } else if (backend === 'generic' || backend === 'large') {
      try {
        exponent = ZZ.__call__((n instanceof RealDoubleElement ? n.value : n) as IntegerInput);
      } catch (e) {
        if (e instanceof TypeError) throw new TypeError('non-integral exponents not supported');
        throw e;
      }
    } else {
      if (n === null || n === undefined || typeof n === 'string' || Array.isArray(n))
        throw new TypeError('an integer is required');
      const value = n instanceof RealDoubleElement ? n.value : n;
      if (
        (typeof value === 'number' && Number.isFinite(value) && !Number.isInteger(value)) ||
        (value instanceof Rational && value.denominator !== 1n)
      )
        throw new TypeError('Only integral powers defined.');
      exponent =
        value instanceof Polynomial ? value._integer_(ZZ) : ZZ.__call__(value as IntegerInput);
    }
    if (modulus !== undefined && modulus !== null)
      return polynomialModularPower(this, exponent, modulus);
    const minimum = -(1n << 63n),
      maximum = (1n << 63n) - 1n;
    if (
      (backend === 'integer' || backend === 'rational') &&
      (exponent < minimum || exponent > maximum)
    )
      throw new OverflowError(
        typeof n === 'number' || n instanceof Rational
          ? 'Python int too large to convert to C long'
          : 'Sage Integer too large to convert to C long'
      );
    const scalarPower = (c: C, e: bigint): C => (c as unknown as { pow(e: bigint): C }).pow(e);
    const inverse = (f: Polynomial<C>): FractionElement<C> => {
      if (!polynomialDomain(f.parent))
        throw new TypeError(`unsupported operand parent(s) for /: '${f.parent}' and '${f.parent}'`);
      return FractionField(f.parent as PolynomialRing<C>).__call__(1n, f);
    };
    // Generic polynomial powers also provide the overflow fallback of the templates.
    const generic = (): Polynomial<C> | FractionElement<C> => {
      if (this.degree() <= 0) return this.parent.__call__(scalarPower(this.getCoeff(0), exponent));
      if (exponent < 0n) return inverse(this).pow(-exponent);
      if (
        this.degree() === 1 &&
        this.getCoeff(0).isZero() &&
        this.leading_coefficient().eq(this.parent.base_ring.one())
      ) {
        if (exponent > maximum)
          throw new OverflowError("cannot fit 'int' into an index-sized integer");
        const coeffs = Array<C>(Number(exponent) + 1).fill(this.parent.base_ring.zero());
        coeffs[Number(exponent)] = this.parent.base_ring.one();
        return new Polynomial(coeffs, this.parent);
      }
      const power = (a: Polynomial<C>, e: bigint): Polynomial<C> => {
        if (e === 0n) return this.parent.one();
        if (e === 1n) return a;
        let base = a,
          n = e;
        while (!(n & 1n)) {
          base = base.mul(base);
          n >>= 1n;
        }
        let result = base;
        n >>= 1n;
        while (n) {
          base = base.mul(base);
          if (n & 1n) result = result.mul(base);
          n >>= 1n;
        }
        return result;
      };
      const characteristic = getRingCharacteristic(this.parent.base_ring);
      if (
        exponent > 20n &&
        characteristic !== null &&
        characteristic > 0n &&
        characteristic <= exponent &&
        (polynomialDomain(this.parent.base_ring) || is_prime(characteristic))
      ) {
        let result = this.parent.one(),
          q = exponent,
          e = 1n;
        while (q) {
          const r = q % characteristic;
          q /= characteristic;
          if (r) {
            if (e * BigInt(this.degree()) >= maximum)
              throw new OverflowError("cannot fit 'int' into an index-sized integer");
            const coeffs = Array<C>(Number(e) * this.degree() + 1).fill(
              this.parent.base_ring.zero()
            );
            for (let i = 0; i < this.coeffs.length; i++)
              coeffs[Number(e) * i] = scalarPower(this.coeffs[i]!, e);
            result = result.mul(power(new Polynomial(coeffs, this.parent), r));
          }
          e *= characteristic;
        }
        return result;
      }
      return power(this, exponent);
    };
    if (backend === 'generic') return generic();
    if (backend === 'large') {
      if (this.degree() <= 0) {
        if (this.isZero() && exponent < 0n) throw new ZeroDivisionError('Inverse does not exist.');
        return this.parent.__call__(scalarPower(this.getCoeff(0), exponent));
      }
      if (exponent < 0n) {
        const positive = this.pow(-exponent) as Polynomial<C>;
        return inverse(positive);
      }
      if (exponent > maximum) throw new OverflowError('Python int too large to convert to C long');
      return new Polynomial(
        ZZ_pX_power(
          this.coeffs.map((c) => (c as unknown as { value: bigint }).value),
          exponent,
          getRingCharacteristic(this.parent.base_ring)!
        ).map((c) => this.parent.base_ring.__call__(c)),
        this.parent
      );
    }
    if (exponent < minimum || exponent > maximum) return generic();
    if (this.isZero()) {
      if (exponent < 0n && (backend === 'integer' || backend === 'rational'))
        throw new ZeroDivisionError('negative exponent in power of zero');
      return exponent === 0n ? this.parent.one() : this.parent.zero();
    }
    const e = exponent < 0n ? -exponent : exponent;
    let result: Polynomial<C>;
    if (backend === 'integer')
      result = new Polynomial(
        _fmpz_poly_pow(extractIntegerCoeffs(this), e).map((c) => this.parent.base_ring.__call__(c)),
        this.parent
      );
    else if (backend === 'rational') {
      const [a, d] = polynomialRationalData(this),
        [coeffs, den] = _fmpq_poly_pow(a, d, e);
      const base = this.parent.base_ring,
        denominator = base.__call__(den);
      result = new Polynomial(
        coeffs.map((c) => (base.__call__(c) as unknown as { div(d: C): C }).div(denominator)),
        this.parent
      );
    } else if (backend === 'word')
      result = new Polynomial(
        _nmod_poly_pow(
          this.coeffs.map((c) => (c as unknown as { value: bigint }).value),
          e,
          getRingCharacteristic(this.parent.base_ring)!
        ).map((c) => this.parent.base_ring.__call__(c)),
        this.parent
      );
    else {
      if (e > maximum) throw new NTLError('power: negative exponent');
      if (backend === 'binary') {
        const f = GF2X.power(
          new GF2X(
            this.coeffs.map((c) => Number((c as unknown as { value: bigint | number }).value))
          ),
          e
        );
        const coeffs =
          f.rep() === 0n
            ? []
            : f
                .rep()
                .toString(2)
                .split('')
                .reverse()
                .map((c) => this.parent.base_ring.__call__(BigInt(c)));
        result = new Polynomial(coeffs, this.parent);
      } else {
        const base = this.parent.base_ring as unknown as FiniteFieldElement['parent'];
        const f = base.modulus.coeffs.map((c) => c.value);
        const a = this.coeffs.map((c) =>
          (c as unknown as FiniteFieldElement).lift.coeffs.map((v) => v.value)
        );
        result = new Polynomial(
          ZZ_pEX_power(a, e, f, base.characteristic).map((c) => base.__call__(c) as unknown as C),
          this.parent
        );
      }
    }
    return exponent < 0n ? inverse(result) : result;
  }

  /** Exact polynomial root. @see Deviation: Polynomial Roots and Truncated Series */
  nth_root(n: unknown): Polynomial<C> {
    const base = this.parent.base_ring;
    if (!polynomialDomain(base))
      throw new ValueError(
        'n-th root of polynomials over rings with zero divisors not implemented'
      );
    if (n === null || n === undefined || typeof n === 'string' || Array.isArray(n))
      throw new TypeError(
        `'<=' not supported between instances of '${polynomialScalarType(n)}' and 'int'`
      );
    const raw =
      n instanceof FiniteFieldElement
        ? n.integer_representation()
        : n instanceof Rational
          ? n.numerator
          : typeof n === 'number' || typeof n === 'bigint'
            ? n
            : n instanceof GF2Element
              ? BigInt(n.value)
              : typeof n === 'boolean'
                ? BigInt(n)
                : ZZ.__call__(n as IntegerInput);
    if (raw <= 0)
      throw new ValueError(
        `n (=${typeof n === 'boolean' ? (n ? 'True' : 'False') : n}) must be positive`
      );
    const one = n instanceof Rational ? n.numerator === n.denominator : raw === 1 || raw === 1n;
    if (one || this.isZero() || this.eq(this.parent.one())) return this;
    if (typeof n === 'number' && !Number.isInteger(n))
      throw new TypeError(
        "unsupported operand type(s) for %: 'sage.rings.integer.Integer' and 'float'"
      );
    if (n instanceof FiniteFieldElement)
      throw new TypeError(`unsupported operand parent(s) for %: '${n.parent}' and '${n.parent}'`);
    const e = ZZ.__call__(n as IntegerInput),
      ordinal = `${e}${e % 100n !== 11n && e % 10n === 1n ? 'st' : e % 100n !== 12n && e % 10n === 2n ? 'nd' : e % 100n !== 13n && e % 10n === 3n ? 'rd' : 'th'}`;
    if (
      n instanceof IntegerMod ||
      n instanceof PrimeFieldElement ||
      n instanceof LegacyPrimeElement
    ) {
      const characteristic =
        n instanceof IntegerMod
          ? n.modulus
          : n instanceof PrimeFieldElement
            ? n.parent.characteristic
            : n.p;
      if (characteristic % e !== 0n) throw new ArithmeticError(`reduction modulo ${n} not defined`);
    }
    if (BigInt(this.degree()) % e) throw new ValueError(`not a ${ordinal} power`);
    if (this.getCoeff(0).isZero()) {
      const valuation = this.coeffs.findIndex((c) => !c.isZero());
      if (BigInt(valuation) % e) throw new ValueError(`not a ${ordinal} power`);
      return this.shift(-valuation)!
        .nth_root(e)
        .shift(BigInt(valuation) / e)!;
    }
    const c = this.getCoeff(0),
      start = c.eq(base.one())
        ? this.parent.one()
        : this.parent.__call__(_coefficient_nth_root(c, e) as C);
    let p: Polynomial<RingElement> = this;
    if (polynomialBackend(base) === 'integer')
      p = new PolynomialRing(QQ, this.parent.variable_name).__call__(
        this
      ) as unknown as Polynomial<RingElement>;
    else if (base instanceof PolynomialRing)
      p = new PolynomialRing(FractionField(base), this.parent.variable_name).__call__(this);
    const q = p._nth_root_series(e, BigInt(this.degree()) / e + 1n, start);
    if (q.pow(e).eq(p)) return this.parent.__call__(q);
    throw new ValueError(`not a ${ordinal} power`);
  }
  /** Newton root series. @see Deviation: Polynomial Roots and Truncated Series */
  _nth_root_series(n: unknown, prec: unknown, start?: unknown): Polynomial<C> {
    let m = polynomialInteger(n, 'long');
    const precision = Number(polynomialInteger(prec, 'long')),
      base = this.parent.base_ring;
    const ordinal = `${m}${m % 100n !== 11n && m % 10n === 1n ? 'st' : m % 100n !== 12n && m % 10n === 2n ? 'nd' : m % 100n !== 13n && m % 10n === 3n ? 'rd' : 'th'}`;
    if (m <= 0n) throw new ValueError(`n (=${m}) must be positive`);
    if (m === 1n || this.isZero() || this.eq(this.parent.one())) return this;
    if (this.getCoeff(0).isZero()) {
      const valuation = this.coeffs.findIndex((c) => !c.isZero());
      if (BigInt(valuation) % m) throw new ValueError(`not a ${ordinal} power`);
      return this.shift(-valuation)!
        ._nth_root_series(m, precision - Number(BigInt(valuation) / m))
        .shift(BigInt(valuation) / m)!;
    }
    let p: Polynomial<C> = this;
    const characteristic = getRingCharacteristic(base) ?? 0n;
    if (characteristic > 0n && m % characteristic === 0n) {
      let cc = 1n;
      while (m % characteristic === 0n) {
        cc *= characteristic;
        m /= characteristic;
      }
      const out = Array<C>(Math.floor(this.degree() / Number(cc)) + 1).fill(base.zero());
      for (let i = 0; i < this.coeffs.length; i++) {
        const c = this.coeffs[i]!;
        if (c.isZero()) continue;
        if (BigInt(i) % cc) throw new ValueError(`not a ${ordinal} power`);
        out[Number(BigInt(i) / cc)] = _coefficient_nth_root(c, cc) as C;
      }
      p = new Polynomial(out, this.parent);
      if (m === 1n) return p;
    }
    const a =
      start !== undefined && start !== null
        ? base.__call__(start)
        : p.getCoeff(0).eq(base.one())
          ? base.one()
          : (_coefficient_nth_root(p.getCoeff(0), m) as C);
    const inverse = (c: C): C | null => {
      if (polynomialBackend(base) === 'integer') {
        const value = ZZ.__call__(c as unknown as IntegerInput);
        return value === 1n || value === -1n ? base.__call__(value) : null;
      }
      return inverseOfUnit(c, base);
    };
    const ai = inverse(a);
    if (ai === null) throw new ArithmeticError('constant coefficient not invertible in base ring');
    const mi = inverse(base.__call__(m));
    if (mi === null) throw new ArithmeticError('exponent not invertible in base ring');
    if (precision < 1) throw new ValueError(`N (=${precision}) must be a positive integer`);
    const sizes = [precision];
    while (sizes[sizes.length - 1]! > 1) sizes.push(Math.ceil(sizes[sizes.length - 1]! / 2));
    let q = this.parent.__call__(ai);
    for (const size of sizes.reverse())
      q = q
        .scalar_mul(base.__call__(m + 1n))
        .sub(p._mul_trunc_(q._power_trunc(m + 1n, size), size))
        .scalar_mul(mi);
    return q.inverse_series_trunc(precision);
  }
  /** Reciprocal series. @see Deviation: Polynomial Roots and Truncated Series */
  inverse_series_trunc(prec: unknown): Polynomial<C> {
    const base = this.parent.base_ring,
      backend = polynomialBackend(base);
    if (backend === 'extension') {
      if (prec === null || prec === undefined || typeof prec === 'string' || Array.isArray(prec))
        throw new TypeError(
          `'<=' not supported between instances of '${polynomialScalarType(prec)}' and 'int'`
        );
      const raw =
        prec instanceof Rational
          ? Number(prec.numerator) / Number(prec.denominator)
          : prec instanceof FiniteFieldElement
            ? prec.integer_representation()
            : typeof prec === 'number' || typeof prec === 'bigint'
              ? prec
              : ZZ.__call__(prec as IntegerInput);
      if (raw <= 0)
        throw new ValueError(
          `the precision must be positive, got ${typeof prec === 'boolean' ? (prec ? 'True' : 'False') : prec}`
        );
      const constant = this.getCoeff(0);
      if (constant.isZero()) throw new ValueError(`constant term ${constant} is not a unit`);
      if (typeof raw === 'number' && (Number.isNaN(raw) || raw < 1)) return this.parent.zero();
      if (prec instanceof FiniteFieldElement) {
        // Sage's generated C-long conversion runs while the PARI stack is guarded.
        // Its temporary conversion fails with this deterministic wrapper error.
        const error = new Error('calling remove_from_pari_stack() inside sig_on()');
        error.name = 'SystemError';
        throw error;
      }
    }
    const n = Number(polynomialInteger(prec, 'long'));
    if (n <= 0) throw new ValueError(`the precision must be positive, got ${n}`);
    if (backend === 'integer') {
      if (this.isZero()) throw new ValueError('constant term is zero');
      const a = extractIntegerCoeffs(this);
      if (a[0] !== 1n && a[0] !== -1n)
        throw new ValueError(`constant term ${this.getCoeff(0)} is not a unit`);
      return new Polynomial(
        _fmpz_poly_inv_series(a, n).map((c) => base.__call__(c)),
        this.parent
      );
    }
    if (backend === 'rational') {
      if (this.getCoeff(0).isZero()) throw new ValueError('constant term is zero');
      const [a, d] = polynomialRationalData(this),
        [out, den] = _fmpq_poly_inv_series_newton(a, d, n),
        D = base.__call__(den);
      return new Polynomial(
        out.map((c) => divideCoeffs(base.__call__(c), D)),
        this.parent
      );
    }
    if (backend === 'extension') {
      const field = base as unknown as FiniteFieldElement['parent'];
      const out = ZZ_pEX_InvTrunc(
        this.coeffs.map((c) =>
          (c as unknown as FiniteFieldElement).lift.coeffs.map((x) => x.value)
        ),
        n,
        field.modulus.coeffs.map((c) => c.value),
        field.characteristic
      );
      return new Polynomial(
        out.map((c) => base.__call__(c)),
        this.parent
      );
    }
    const c = this.getCoeff(0),
      first = inverseOfUnit(c, base);
    if (first === null) throw new ValueError(`constant term ${c} is not a unit`);
    let current = this.parent.__call__(first);
    const sizes = [n];
    while (sizes[sizes.length - 1]! > 1) sizes.push(Math.ceil(sizes[sizes.length - 1]! / 2));
    for (const next of sizes.reverse().slice(1)) {
      const z = current._mul_trunc_(this, next)._mul_trunc_(current, next);
      current = current.add(current).sub(z);
    }
    return current;
  }
  power_trunc(n: unknown, prec: unknown): Polynomial<C> {
    const e = ZZ.__call__(n as IntegerInput);
    if (e >= 0n && e < 1n << 64n) return this._power_trunc(e, prec);
    return generic_power_trunc(this, e, Number(polynomialInteger(prec, 'index')));
  }
  _power_trunc(n: unknown, prec: unknown): Polynomial<C> {
    const e = polynomialInteger(n, 'unsigned'),
      precision = Number(polynomialInteger(prec, 'long'));
    const base = this.parent.base_ring,
      backend = polynomialBackend(base);
    if (backend === 'integer' || backend === 'word') {
      if (precision <= 0) return this.parent.zero();
      const a = extractIntegerCoeffs(this);
      const out =
        backend === 'integer'
          ? _fmpz_poly_pow_trunc(a, e, precision)
          : _nmod_poly_pow_trunc(a, e, precision, getRingCharacteristic(base)!);
      return new Polynomial(
        out.map((c) => base.__call__(c)),
        this.parent
      );
    }
    return generic_power_trunc(this, e, precision);
  }
  /** @see Deviation: Truncated Multiplication Parents and Real Coefficients */
  _mul_trunc_(right: Polynomial<C> | null, n: unknown): Polynomial<C>;
  _mul_trunc_(right: unknown, n: unknown): Polynomial<C> {
    const precision = Number(polynomialInteger(n, 'long')),
      base = this.parent.base_ring,
      backend = polynomialBackend(base);
    if (right === null || right === undefined) right = this.parent.zero();
    if (!(right instanceof Polynomial)) {
      const type =
        right instanceof RealDoubleElement
          ? 'sage.rings.real_double_element_gsl.RealDoubleElement_gsl'
          : isFractionElement(right)
            ? right.constructor.name === 'FpTElement'
              ? 'sage.rings.fraction_field_FpT.FpTElement'
              : `sage.rings.fraction_field_element.${right.constructor.name}`
            : polynomialPseudoType(right);
      throw new TypeError(
        `Argument 'right' has incorrect type (expected sage.rings.polynomial.polynomial_element.Polynomial, got ${type})`
      );
    }
    if (backend === 'integer' || backend === 'rational' || backend === 'word') {
      if (precision <= 0)
        throw new ValueError(backend === 'rational' ? 'n must be > 0' : 'length must be > 0');
      if (backend === 'rational') {
        const [a, da] = polynomialRationalData(this),
          [b, db] = polynomialRationalData(right),
          [out, d] = _fmpq_poly_mullow(a, da, b, db, precision),
          D = base.__call__(d);
        return new Polynomial(
          out.map((c) => divideCoeffs(base.__call__(c), D)),
          this.parent
        );
      }
      const a = extractIntegerCoeffs(this),
        b = right === this ? a : extractIntegerCoeffs(right);
      const out =
        backend === 'integer'
          ? _fmpz_poly_mullow(a, b, precision)
          : _nmod_poly_mullow(a, b, precision, getRingCharacteristic(base)!);
      return new Polynomial(
        out.map((c) => base.__call__(c)),
        this.parent
      );
    }
    if (this.isZero() || right.isZero()) return this.parent.zero();
    // Generic polynomial_element.pyx uses schoolbook below its ring threshold,
    // otherwise multiplying truncated operands has the native generic complexity.
    const threshold =
      base instanceof PolynomialRing
        ? 0
        : base instanceof FractionField_generic
          ? Number.MAX_SAFE_INTEGER
          : 8;
    if (precision < threshold) {
      return new Polynomial(
        do_schoolbook_product(this.coeffs, right.coeffs, precision),
        this.parent
      );
    }
    return this.truncate(precision).mul(right.truncate(precision)).truncate(precision);
  }
  /** @see Deviation: Truncated Multiplication Parents and Real Coefficients */
  multiplication_trunc(other: unknown, n: unknown): Polynomial<RingElement> {
    if (other instanceof Polynomial) {
      const [a, r] = polynomialCommonOperands(this, other, 'multiplication_trunc');
      return a._mul_trunc_(r, polynomialInteger(n, 'index'));
    }
    const primitiveFloat = typeof other === 'number' && !Number.isInteger(other);
    const primitiveInvalid =
      other === null || other === undefined || typeof other === 'string' || Array.isArray(other);
    const input =
      typeof other === 'boolean'
        ? BigInt(other)
        : primitiveFloat
          ? RDF.__call__(other as number)
          : other;
    const rightParent = primitiveInvalid ? null : polynomialElementParent(input);
    const zero = this.parent.base_ring.zero();
    const scalarZeroRing =
      other instanceof IntegerMod &&
      other.modulus === 1n &&
      zero instanceof IntegerMod &&
      zero.modulus !== 1n;
    const common =
      rightParent === null || scalarZeroRing
        ? null
        : polynomialCommonBase(this.parent, rightParent);
    if (common === null) {
      const label =
        primitiveInvalid || primitiveFloat
          ? `<class '${polynomialScalarType(other)}'>`
          : String(rightParent);
      throw new TypeError(
        `no common canonical parent for objects with parents: '${this.parent}' and '${label}'`
      );
    }
    if (common instanceof FractionField_generic) {
      const value = common.__call__(this);
      const type =
        value.constructor.name === 'FpTElement'
          ? 'sage.rings.fraction_field_FpT.FpTElement'
          : `sage.rings.fraction_field_element.${value.constructor.name}`;
      throw new TypeError(
        `Cannot convert ${type} to sage.rings.polynomial.polynomial_element.Polynomial`
      );
    }
    const parent = common as PolynomialRing<RingElement>;
    const a = parent.__call__(this),
      b = parent.__call__(input);
    return a._mul_trunc_(b, polynomialInteger(n, 'index'));
  }

  /** Evaluate using the native backend and Sage's canonical scalar coercions.
   * @see Deviation: Polynomial Evaluation and Composition
   * @see Deviation: Polynomial Matrix Evaluation Actions
   */
  evaluate(x: C): C;
  evaluate(x: bigint | Integer | boolean): C;
  evaluate(x: EvaluationMatrix): EvaluationMatrix;
  evaluate(x: readonly unknown[]): RingElement | number | EvaluationMatrix;
  evaluate<T>(
    x: T
  ): 0 extends 1 & T
    ? C
    : T extends C
      ? C
      : T extends EvaluationMatrix
        ? EvaluationMatrix
        : T extends object
          ? RingElement
          : number | RingElement | EvaluationMatrix;
  evaluate(x: { isZero(): boolean }): RingElement;
  evaluate(x?: unknown): RingElement | number | EvaluationMatrix;
  evaluate(x?: unknown): RingElement | number | EvaluationMatrix {
    return polynomialEvaluate(this, x);
  }

  /**
   * Check equality.
   * @see Deviation: Polynomial Common Coefficient Parents and Representation
   * @see Deviation: Polynomial Scalar Equality and Coefficient Embedding
   */
  eq(other: unknown): boolean {
    if (!(other instanceof Polynomial)) {
      if (isFractionElement(other)) {
        const common = polynomialCommonBase(this.parent, other.parent);
        if (common instanceof FractionField_generic)
          return common.__call__(this).eq(common.__call__(other));
        if (common instanceof PolynomialRing)
          return common.__call__(this).eq(common.__call__(other));
        return false;
      }
      // A zero-ring scalar has no canonical map into a nontrivial modular
      // polynomial ring. Wrapping it in a polynomial would introduce a
      // quotient-parent pushout which Sage's scalar comparison does not use.
      if (other instanceof IntegerMod && other.modulus === 1n) {
        const zero = this.parent.base_ring.zero();
        if (zero instanceof IntegerMod && zero.modulus !== 1n) return false;
      }
      if (typeof other === 'number' && !Number.isInteger(other)) {
        if (!polynomialRealBase(this.parent.base_ring)) return false;
        const constant = this.coeffs.length ? polynomialRealCoefficient(this.coeffs[0]!) : 0;
        return (
          constant === other &&
          this.coeffs.slice(1).every((c) => polynomialRealCoefficient(c) === 0)
        );
      }
      if (
        typeof other === 'bigint' ||
        typeof other === 'boolean' ||
        typeof other === 'number' ||
        other instanceof Integer
      ) {
        const scalar = other instanceof Integer ? other.value : BigInt(other);
        return this.eq(this.parent.__call__(scalar));
      }
      let scalarBase: CoefficientRing<RingElement>;
      if (other instanceof Rational) scalarBase = QQ as unknown as CoefficientRing<RingElement>;
      else if (
        other instanceof IntegerMod ||
        other instanceof PrimeFieldElement ||
        other instanceof LegacyPrimeElement ||
        other instanceof FiniteFieldElement ||
        other instanceof GF2Element
      )
        scalarBase = other.parent as unknown as CoefficientRing<RingElement>;
      else return false;
      return this.eq(
        new Polynomial(
          [other as RingElement],
          new PolynomialRing(scalarBase, this.parent.variable_name)
        )
      );
    }
    if (this.parent !== other.parent) {
      if (!this._hasCompatibleParent(other)) return false;
      const [a, b] = polynomialCommonOperands(this, other, '==');
      return a.eq(b);
    }
    if (this.coeffs.length !== other.coeffs.length) {
      return false;
    }
    for (let i = 0; i < this.coeffs.length; i++) {
      if (!this.coeffs[i]!.eq(other.coeffs[i]!)) {
        return false;
      }
    }
    return true;
  }

  /**
   * Differentiate with respect to a sequence of variables and repetition counts.
   *
   * @returns The formal derivative d/dx of this polynomial
   *
   * @example
   * ```typescript
   * // If f = x^3 + 2x + 1, then f.derivative() = 3x^2 + 2
   * ```
   *
   * @see Reference: sage/rings/polynomial/polynomial_element.pyx:derivative
   * @see Deviation: Polynomial Derivative Protocol and Native Coefficients
   */
  derivative(...args: unknown[]): Polynomial<C> {
    return multi_derivative<Polynomial<C>>(this, args);
  }

  get diff(): (...args: unknown[]) => Polynomial<C> {
    return this.derivative;
  }
  get differentiate(): (...args: unknown[]) => Polynomial<C> {
    return this.derivative;
  }

  /**
   * Differentiate once with respect to a generator, recursing into coefficients.
   * @see Deviation: Polynomial Derivative Protocol and Native Coefficients
   */
  _derivative(variable?: unknown): Polynomial<C>;
  _derivative(...args: unknown[]): Polynomial<C> {
    if (args.length > 1)
      throw new TypeError(
        `_derivative() takes at most 1 positional argument (${args.length} given)`
      );
    const variable = args[0],
      base = this.parent.base_ring;
    const rational = polynomialBackend(base) === 'rational';
    if (variable !== undefined && variable !== null && !this.parent.gen().eq(variable)) {
      if (rational)
        throw new ValueError(
          `cannot differentiate with respect to ${polynomialDerivativeVariable(variable)}`
        );
      try {
        return this.parent.__call__(
          this.coeffs.map((c) => {
            const f = (c as unknown as { _derivative?: (variable: unknown) => C })._derivative;
            if (typeof f !== 'function')
              throw new AttributeError('coefficient has no _derivative method');
            return f.call(c, variable);
          })
        );
      } catch (e) {
        if (e instanceof AttributeError)
          throw new ValueError(
            `cannot differentiate with respect to ${polynomialDerivativeVariable(variable)}`
          );
        throw e;
      }
    }
    if (rational) {
      const [a, den] = polynomialRationalData(this);
      const [coeffs, d] = _fmpq_poly_derivative(a, den);
      return new Polynomial(
        coeffs.map((c) => base.__call__(new Rational(c, d))),
        this.parent
      );
    }
    if (this.isZero()) return this;
    if (this.isConstant()) return this.parent.zero();
    // Sage multiplies the integer degree by the coefficient in its base ring.
    // A single scalar product also preserves RDF rounding, unlike double-and-add.
    return new Polynomial(
      this.coeffs.slice(1).map((c, i) => base.__call__(BigInt(i + 1)).mul(c)),
      this.parent
    );
  }

  /** The one partial derivative of a univariate polynomial. */
  gradient(...args: unknown[]): Polynomial<C>[] {
    if (args.length)
      throw new TypeError(`gradient() takes exactly 0 positional arguments (${args.length} given)`);
    return [this.diff()];
  }

  /**
   * Return the GCD of this polynomial and other.
   *
   * Delegates ZZ/QQ/word-modular coefficients to FLINT and binary coefficients to NTL.
   *
   * @param other - Another polynomial in the same ring
   * @returns The GCD; native finite zero shortcuts preserve the nonmonic operand.
   *
   * @see Deviation: Polynomial GCD Backend Coercion and Identity
   * @see Reference: sage/rings/polynomial/polynomial_element.pyx:gcd
   */
  gcd(other: Polynomial<C>): Polynomial<C>;
  gcd<D extends RingElement>(other: Polynomial<D>): Polynomial<C | D>;
  gcd(other: Polynomial<RingElement>): Polynomial<RingElement> {
    if (this.parent !== other.parent) {
      const [a, b] = polynomialCommonOperands(this, other, 'gcd');
      const g = a.gcd(b);
      // Sage interns these parents. Preserve operand-returning shortcuts when
      // separate TypeScript instances have canonical maps in both directions.
      if (
        g === a &&
        polynomialCanCoerce(this.parent, g.parent) &&
        polynomialCanCoerce(g.parent, this.parent)
      )
        return this;
      if (
        g === b &&
        polynomialCanCoerce(other.parent, g.parent) &&
        polynomialCanCoerce(g.parent, other.parent)
      )
        return other;
      return g;
    }
    const rhs = other as Polynomial<C>;
    const baseRing = this.parent.base_ring;

    // Sage's integer backend dispatches among FLINT's subresultant, heuristic
    // and modular kernels; preserve that dispatch in the dependency port.
    if (isIntegerRing(baseRing)) {
      // polynomial_integer_dense_flint.pyx:822-825 returns the operand
      // directly on zero/one shortcuts, before FLINT normalizes its sign.
      const isOne = (f: Polynomial<C>) => f.degree() === 0 && f.getCoeff(0).eq(baseRing.one());
      if (this.isZero() || isOne(rhs)) return rhs;
      if (isOne(this) || rhs.isZero()) return this;
      const a = extractIntegerCoeffs(this);
      const b = extractIntegerCoeffs(rhs);
      const g = _fmpz_poly_gcd(a, b);
      return new Polynomial(
        g.map((c) => baseRing.__call__(c) as C),
        this.parent
      );
    }

    const backend = polynomialBackend(baseRing);
    if (backend === 'rational') {
      const [g, den] = _fmpq_poly_gcd(
        polynomialRationalData(this)[0],
        polynomialRationalData(rhs)[0]
      );
      return new Polynomial(
        g.map((c) => baseRing.__call__(new Rational(c, den))),
        this.parent
      );
    }
    if (backend === 'word' || backend === 'binary' || backend === 'extension') {
      // polynomial_template.pxi:gcd handles these before native dispatch.
      if (this.isZero()) return rhs;
      if (rhs.isZero()) return this;
      if (this.eq(rhs)) return this.monic();
    }
    if (backend === 'word') {
      const n = getRingCharacteristic(baseRing)!;
      const lift = (c: C) => (c as unknown as { value: bigint }).value;
      try {
        const g = _nmod_poly_gcd(this.coeffs.map(lift), rhs.coeffs.map(lift), n);
        const monic = g.length === 1 ? [1n] : _nmod_poly_make_monic(g, n);
        return new Polynomial(
          monic.map((c) => baseRing.__call__(c)),
          this.parent
        );
      } catch (e) {
        if (e instanceof RangeError) throw new RuntimeError('FLINT gcd calculation failed');
        throw e;
      }
    }
    if (backend === 'binary') {
      const packed = (f: Polynomial<C>) =>
        new GF2X(
          BigInt(
            '0b' +
              f.coeffs
                .map((c) => String((c as unknown as { value: bigint | number }).value))
                .reverse()
                .join('')
          )
        );
      const g = GF2X.GCD(packed(this), packed(rhs)).rep();
      return new Polynomial(
        g === 0n
          ? []
          : g
              .toString(2)
              .split('')
              .reverse()
              .map((c) => baseRing.__call__(BigInt(c))),
        this.parent
      );
    }
    if (!ringIsField(baseRing)) {
      throw new NotImplementedError(
        `${baseRing} does not provide a gcd implementation for univariate polynomials`
      );
    }

    // Fields: Euclidean algorithm (sage/categories/fields.py:_gcd_univariate_polynomial)
    let a: Polynomial<C> = this;
    let b: Polynomial<C> = rhs;

    while (!b.isZero()) {
      const [_q, r] = a.quo_rem(b);
      a = b;
      b = r;
    }

    // Return monic GCD (zero stays zero)
    const g = a.isZero() ? a : a._monic();
    // These NTL wrappers allocate after a nontrivial native GCD call.
    return backend === 'large' || backend === 'extension'
      ? new Polynomial([...g.coeffs], this.parent)
      : g;
  }

  /**
   * Return the extended GCD of this polynomial and other.
   *
   * Returns (g, s, t) with g = s*this + t*other, preserving native zero conventions.
   *
   * @param other - Another polynomial in the same ring
   * @returns Tuple [g, s, t] where g = s*this + t*other
   *
   * Over ZZ, constant/zero branches can return Integer components. Nonconstant
   * inputs use FLINT's length-ordered resultant or a denominator-cleared QQ triple.
   *
   * @see Deviation: Polynomial Integer and Rational Extended GCD
   * @see Deviation: Polynomial Extended GCD Finite Backends
   * @see Reference: sage/rings/polynomial/polynomial_element.pyx:xgcd
   */
  xgcd(
    this: Polynomial<Integer & RingElement>,
    other: Polynomial<Integer & RingElement>
  ): [
    Polynomial<Integer & RingElement> | Integer,
    Polynomial<Integer & RingElement> | Integer,
    Polynomial<Integer & RingElement> | Integer,
  ];
  xgcd(other: Polynomial<C>): [Polynomial<C>, Polynomial<C>, Polynomial<C>];
  xgcd<D extends RingElement>(
    other: Polynomial<D>
  ): [Polynomial<C | D>, Polynomial<C | D>, Polynomial<C | D>];
  xgcd(
    operand: Polynomial<RingElement>
  ): [
    Polynomial<RingElement> | Integer,
    Polynomial<RingElement> | Integer,
    Polynomial<RingElement> | Integer,
  ] {
    if (this.parent !== operand.parent) {
      const [a, b] = polynomialCommonOperands(this, operand, 'xgcd');
      const result = a.xgcd(b);
      return result.map((g) => {
        if (
          g === a &&
          polynomialCanCoerce(this.parent, g.parent) &&
          polynomialCanCoerce(g.parent, this.parent)
        )
          return this;
        if (
          g === b &&
          polynomialCanCoerce(operand.parent, g.parent) &&
          polynomialCanCoerce(g.parent, operand.parent)
        )
          return operand;
        return g;
      }) as [Polynomial<RingElement>, Polynomial<RingElement>, Polynomial<RingElement>];
    }
    const other = operand as Polynomial<C>;
    const R = this.parent;
    const baseRing = R.base_ring;

    const backend = polynomialBackend(baseRing);
    const materialize = (coeffs: bigint[]) =>
      new Polynomial(
        coeffs.map((c) => baseRing.__call__(c)),
        R
      );
    if (backend === 'integer') {
      if (this.isZero()) return [other, new Integer(0n), new Integer(1n)];
      if (other.isZero()) return [this, new Integer(1n), new Integer(0n)];
      if (this.isConstant() && other.isConstant())
        return (this.getCoeff(0) as unknown as Integer).xgcd(
          other.getCoeff(0) as unknown as Integer
        );
      const A = extractIntegerCoeffs(this),
        B = extractIntegerCoeffs(other);
      const [r, S, T] = _fmpz_poly_xgcd(A, B);
      if (r) return [materialize([r]), materialize(S), materialize(T)];
      // Sage clears a common denominator from the QQ Bezout triple when
      // the native integer resultant vanishes (the input polynomials share a factor).
      const Q = new PolynomialRing(QQ, R.variable_name);
      const rational = Q.__call__(A).xgcd(Q.__call__(B)).map(polynomialRationalData);
      let d = 1n;
      for (const [, den] of rational) d = (d / gcdBigInt(d, den)) * den;
      return rational.map(([poly, den]) => materialize(poly.map((c) => c * (d / den)))) as [
        Polynomial<C>,
        Polynomial<C>,
        Polynomial<C>,
      ];
    }
    if (backend === 'rational') {
      const [A, denA] = polynomialRationalData(this),
        [B, denB] = polynomialRationalData(other);
      return _fmpq_poly_xgcd(A, denA, B, denB).map(
        ([poly, den]) =>
          new Polynomial(
            poly.map((c) => baseRing.__call__(new Rational(c, den))),
            R
          )
      ) as [Polynomial<C>, Polynomial<C>, Polynomial<C>];
    }
    if (backend === 'word' || backend === 'binary' || backend === 'extension') {
      if (this.isZero()) return [other, R.zero(), R.one()];
      if (other.isZero()) return [this, R.one(), R.zero()];
    }
    if (backend === 'word') {
      const n = getRingCharacteristic(baseRing)!;
      const lift = (c: C) => (c as unknown as { value: bigint }).value;
      try {
        return _nmod_poly_xgcd(this.coeffs.map(lift), other.coeffs.map(lift), n).map(
          materialize
        ) as [Polynomial<C>, Polynomial<C>, Polynomial<C>];
      } catch (e) {
        if (e instanceof RangeError)
          throw new ValueError('non-invertible elements encountered during XGCD');
        throw e;
      }
    }
    if (backend === 'binary') {
      const packed = (f: Polynomial<C>) =>
        new GF2X(
          BigInt(
            '0b' +
              f.coeffs
                .map((c) => String((c as unknown as { value: bigint | number }).value))
                .reverse()
                .join('')
          )
        );
      return GF2X.XGCD(packed(this), packed(other)).map((f) =>
        materialize(f.rep() === 0n ? [] : f.rep().toString(2).split('').reverse().map(BigInt))
      ) as [Polynomial<C>, Polynomial<C>, Polynomial<C>];
    }
    if (!ringIsField(baseRing)) {
      throw new NotImplementedError(
        `${baseRing} does not provide an xgcd implementation for univariate polynomials`
      );
    }

    const zero = R.zero();
    const one = R.one();

    // sage/categories/fields.py:526-543 (_xgcd_univariate_polynomial)
    if (other.isZero()) {
      if (this.isZero()) {
        return backend === 'large' ? [zero, one, R.zero()] : [zero, zero, zero];
      }
      const c = divideCoeffs(baseRing.one() as C, this.leading_coefficient());
      return [this.scalar_mul(c), R.__call__(c), zero];
    }
    if (this.isZero()) {
      const c = divideCoeffs(baseRing.one() as C, other.leading_coefficient());
      return [other.scalar_mul(c), zero, R.__call__(c)];
    }

    let u = one;
    let d: Polynomial<C> = this;
    let v1 = zero;
    let v3 = other;

    while (!v3.isZero()) {
      const [q, r] = d.quo_rem(v3);
      const newU = v1;
      const newD = v3;
      v1 = u.sub(v1.mul(q));
      v3 = r;
      u = newU;
      d = newD;
    }

    // v = (d - a*u) // b
    let v = d.sub(this.mul(u)).quo_rem(other)[0];

    if (!d.isZero()) {
      const c = divideCoeffs(baseRing.one() as C, d.leading_coefficient());
      d = d.scalar_mul(c);
      u = u.scalar_mul(c);
      v = v.scalar_mul(c);
    }

    return [d, u, v];
  }

  /**
   * Return the composition f(g) where this polynomial is f.
   *
   * @param other - The polynomial g to substitute for x
   * @returns f(g(x))
   *
   * @example
   * ```typescript
   * // If f = x^2 + 1 and g = x + 1, then f.compose(g) = (x+1)^2 + 1 = x^2 + 2x + 2
   * ```
   *
   * @see Reference: sage/rings/polynomial/polynomial_element.pyx:__call__
   * @see Deviation: Polynomial Evaluation and Composition
   */
  compose(other: Polynomial<C>): Polynomial<C>;
  compose<D extends RingElement>(other: Polynomial<D>): Polynomial<RingElement>;
  compose(other: Polynomial<RingElement>): Polynomial<RingElement> {
    return polynomialEvaluate(this, other) as Polynomial<RingElement>;
  }

  /**
   * Return a monic version of this polynomial (leading coefficient = 1).
   *
   * A monic polynomial has leading coefficient 1. This method divides
   * the polynomial by its leading coefficient.
   *
   * @returns Monic polynomial equal to this / leading_coefficient
   *
   * @example
   * ```typescript
   * // 2x^2 + 4x + 2 becomes x^2 + 2x + 1
   * ```
   *
   * @see Deviation: Polynomial Monic Normalization
   * @see Reference: sage/rings/polynomial/polynomial_element.pyx:monic
   */
  monic(
    this: Polynomial<Integer & RingElement>
  ): Polynomial<Integer & RingElement> | Polynomial<Rational & RingElement>;
  monic(): Polynomial<C>;
  monic(): Polynomial<C> | Polynomial<Integer & RingElement> | Polynomial<Rational & RingElement> {
    const base = this.parent.base_ring;
    const p = getRingCharacteristic(base);
    const zero = base.zero();
    // polynomial_zmod_flint.pyx:838-867 has a backend-specific path and
    // always allocates a result, including already-monic inputs.
    if (p !== null && p > 2n && p < 2n ** 63n && 'value' in zero) {
      const lc = this.leading_coefficient() as C & { value: bigint };
      if (this.isZero() || gcdBigInt(lc.value, p) !== 1n)
        throw new ValueError('leading coefficient must be invertible');
      const result = _nmod_poly_make_monic(
        this.coeffs.map((c) => (c as C & { value: bigint }).value),
        p
      );
      return new Polynomial(
        result.map((c) => base.__call__(c)),
        this.parent
      );
    }
    if (this.is_monic() || base.one().isZero()) return this;
    // Generic Sage monic() inverts the leading coefficient before choosing
    // its result parent. Integer inversion lands in QQ, even for -1.
    if (isIntegerRing(base)) {
      const coeffs = extractIntegerCoeffs(this);
      if (!coeffs.length) throw new ZeroDivisionError('rational division by zero');
      const lead = coeffs[coeffs.length - 1]!;
      const R = new PolynomialRing(QQ, this.parent.variable_name);
      // Rational implements these coefficient operations; the intersection bridges
      // the existing polymorphic-this RingElement constraint at this boundary.
      return R.__call__(coeffs.map((c) => new Rational(c, lead))) as unknown as Polynomial<
        Rational & RingElement
      >;
    }
    const lcInv = divideCoeffs(base.one(), this.leading_coefficient());
    return this.scalar_mul(lcInv);
  }

  /**
   * Internal alias for monic() for backward compatibility.
   */
  _monic(
    this: Polynomial<Integer & RingElement>
  ): Polynomial<Integer & RingElement> | Polynomial<Rational & RingElement>;
  _monic(): Polynomial<C>;
  _monic(): Polynomial<C> | Polynomial<Integer & RingElement> | Polynomial<Rational & RingElement> {
    return this.monic();
  }

  /**
   * Return the content of this polynomial (GCD of all coefficients).
   *
   * The content is the GCD of all coefficients. For the zero polynomial,
   * returns zero. Requires that the coefficient ring supports a gcd method.
   *
   * @returns The content (GCD of coefficients)
   *
   * @example
   * ```typescript
   * // If f = 6x^2 + 4x + 2 over ZZ
   * // f.content() = 2
   * ```
   *
   * @see Reference: sage/rings/polynomial/polynomial_integer_dense_flint.pyx:474 (content)
   */
  content(): C {
    if (this.isZero()) {
      return this.parent.base_ring.zero() as C;
    }

    // Start with the first coefficient
    let g = this.coeffs[0]!;

    // Compute GCD with all other coefficients
    for (let i = 1; i < this.coeffs.length; i++) {
      g = gcdCoeffs(g, this.coeffs[i]!);
      // If GCD is 1 (or a unit), we can stop early
      if (g.eq(1)) {
        break;
      }
    }

    // The sign of the content is the sign of the leading coefficient
    // (`polynomial_integer_dense_flint.pyx:477`, issue #13053):
    //     R(-1).content() == -1,  (-2*x^2-4).content() == -2
    if (isNegative(g) !== isNegative(this.leading_coefficient())) {
      g = g.neg() as C;
    }

    return g;
  }

  /**
   * Return the primitive part of this polynomial (this / content).
   *
   * The primitive part is the polynomial divided by its content, so that
   * the GCD of the resulting coefficients is 1.
   *
   * @returns The primitive part
   *
   * @example
   * ```typescript
   * // If f = 6x^2 + 4x + 2 over ZZ
   * // f.primitive_part() = 3x^2 + 2x + 1
   * ```
   *
   * The leading coefficient of the primitive part is always positive, since
   * {@link content} carries the sign of the leading coefficient (this matches
   * FLINT's `fmpz_poly_primitive_part`, see
   * `polynomial_integer_dense_flint.pyx:1535`).
   *
   * @see Reference: sage/libs/flint/fmpz_poly.pxd (fmpz_poly_primitive_part)
   */
  primitive_part(): Polynomial<C> {
    if (this.isZero()) {
      return this;
    }

    const c = this.content();

    // If content is 1, return self
    if (c.eq(1)) {
      return this;
    }

    // Divide each coefficient by the content
    const newCoeffs = this.coeffs.map((coeff) => divideCoeffs(coeff, c));
    return new Polynomial(newCoeffs, this.parent);
  }

  /**
   * Return this polynomial shifted by n (multiplied by x^n).
   *
   * If n is positive, this is equivalent to multiplying by x^n.
   * If n is negative, terms below x^(-n) are discarded (integer division by x^(-n)).
   *
   * @param n - The shift amount (can be negative for division by x^n)
   * @returns x^n * this (or floor division if n < 0)
   *
   * @example
   * ```typescript
   * // If f = x^2 + 2x + 4
   * // f.shift(2) = x^4 + 2x^3 + 4x^2
   * // f.shift(-1) = x + 2
   * ```
   *
   * @see Deviation: Polynomial Index Conversion and Storage
   * @see Reference: sage/rings/polynomial/polynomial_element.pyx:shift
   */
  shift(n: number | bigint): Polynomial<C>;
  shift(n: unknown): Polynomial<C> | null;
  shift(n: unknown): Polynomial<C> | null {
    const backend = polynomialBackend(this.parent.base_ring);
    let amount: bigint;
    if (backend === 'integer' || backend === 'rational') {
      // Generic Polynomial.shift checks equality/zero before ordering or
      // sequence repetition. NaN falls through both comparisons to None.
      if (this.isZero()) return this;
      const sign = polynomialScalarSign(n, '>');
      if (sign === 0) return this;
      if (sign === null) return null;
      if (n instanceof FiniteFieldElement)
        throw new TypeError(
          `can't multiply sequence by non-int of type '${polynomialScalarType(n)}'`
        );
      if (sign > 0 && (n instanceof Rational || (typeof n === 'number' && !Number.isInteger(n)))) {
        if (n instanceof Rational)
          throw new TypeError(
            "unsupported operand parent(s) for *: '<class 'list'>' and 'Rational Field'"
          );
        throw new TypeError("can't multiply sequence by non-int of type 'float'");
      }
      // Negative generic shifts use int(n), then Python slice clipping.
      if (n instanceof Rational) amount = n.numerator / n.denominator;
      else if (typeof n === 'number') {
        if (!Number.isFinite(n))
          throw new OverflowError('cannot convert float infinity to integer');
        amount = BigInt(Math.trunc(n));
      } else
        amount =
          n instanceof Integer ? n.value : typeof n === 'boolean' ? BigInt(n) : (n as bigint);
      if (amount > (1n << 63n) - 1n)
        throw new IndexError("cannot fit 'sage.rings.integer.Integer' into an index-sized integer");
    } else {
      if (backend === 'large') {
        // Polynomial_dense_mod_n.shift tests the original n before NTL's
        // long conversion. A nonzero fraction truncated to zero allocates.
        if (
          this.isZero() ||
          n === 0n ||
          n === 0 ||
          n === false ||
          (n instanceof Integer && n.value === 0n) ||
          (n instanceof Rational && n.numerator === 0n) ||
          (typeof n === 'object' &&
            n !== null &&
            'isZero' in n &&
            typeof n.isZero === 'function' &&
            n.isZero())
        )
          return this;
      }
      amount = polynomialInteger(
        n,
        backend === 'generic' ? 'ssize' : backend === 'large' ? 'long' : 'int'
      );
      if (
        backend !== 'extension' &&
        backend !== 'large' &&
        (amount === 0n || (this.isZero() && backend === 'generic'))
      )
        return this;
    }
    if (this.isZero()) return new Polynomial([], this.parent);
    if (amount <= 0n) {
      const count = -amount >= BigInt(this.coeffs.length) ? this.coeffs.length : Number(-amount);
      return new Polynomial(this.coeffs.slice(count), this.parent);
    }
    return new Polynomial(
      Array.from({ length: Number(amount) }, () => this.parent.base_ring.zero()).concat(
        this.coeffs
      ),
      this.parent
    );
  }

  /**
   * Return the truncation of this polynomial to degree < n.
   *
   * Returns the polynomial with all terms of degree >= n removed.
   *
   * @param n - The C-long degree bound (generic dense rings retain negative slice bounds)
   * @returns Polynomial with terms of degree >= n removed
   *
   * @example
   * ```typescript
   * // If f = x^3 + 2x^2 + 3x + 4
   * // f.truncate(2) = 3x + 4 (terms of degree < 2)
   * // f.truncate(0) = 0 (no terms of degree < 0)
   * ```
   *
   * @see Reference: sage/rings/polynomial/polynomial_element.pyx:truncate
   */
  truncate(n: unknown): Polynomial<C> {
    const length = polynomialInteger(n, 'long');
    const backend = polynomialBackend(this.parent.base_ring);
    if (
      backend !== 'generic' &&
      backend !== 'integer' &&
      backend !== 'large' &&
      length >= BigInt(this.coeffs.length)
    )
      return this;
    if (backend === 'generic') {
      // Polynomial_generic_dense uses Python slicing even for negative n.
      const clipped =
        length < -BigInt(this.coeffs.length)
          ? -this.coeffs.length
          : length > BigInt(this.coeffs.length)
            ? this.coeffs.length
            : Number(length);
      return new Polynomial(this.coeffs.slice(0, clipped), this.parent);
    }
    const clipped =
      length <= 0n ? 0 : length >= BigInt(this.coeffs.length) ? this.coeffs.length : Number(length);
    return new Polynomial(this.coeffs.slice(0, clipped), this.parent);
  }

  /**
   * Return the reverse of this polynomial.
   *
   * If f(x) = a_0 + a_1*x + ... + a_n*x^n, then
   * reverse(f)(x) = a_n + a_{n-1}*x + ... + a_0*x^n = x^n * f(1/x).
   *
   * If an optional degree argument is given, the coefficient list will be
   * truncated or zero-padded as necessary before reversing.
   *
   * @param degree - Optional degree to use (pads with zeros or truncates)
   * @returns Polynomial with coefficients in reverse order
   *
   * @example
   * ```typescript
   * // If f = x^3 + 2x + 3 (coeffs: [3, 2, 0, 1])
   * // f.reverse() = 1 + 0*x + 2*x^2 + 3*x^3 = 1 + 2x^2 + 3x^3
   * ```
   *
   * @see Deviation: Polynomial Index Conversion and Storage
   * @see Reference: sage/rings/polynomial/polynomial_element.pyx:reverse
   */
  reverse(degree?: unknown): Polynomial<C> {
    const backend = polynomialBackend(this.parent.base_ring);
    let length = BigInt(this.coeffs.length);
    if (degree !== undefined && degree !== null) {
      if (backend === 'rational') {
        // FLINT QQ reverse converts degree+1, so e.g. -1 becomes length zero
        // and fractional values are truncated after the addition.
        try {
          const plusOne =
            degree instanceof Rational
              ? degree.add(1n)
              : degree instanceof Integer
                ? degree.value + 1n
                : typeof degree === 'bigint'
                  ? degree + 1n
                  : typeof degree === 'boolean'
                    ? BigInt(degree) + 1n
                    : typeof degree === 'number'
                      ? degree + 1
                      : polynomialScalarAddOne(degree);
          if (plusOne === null) throw new TypeError('an integer is required');
          length = polynomialInteger(plusOne, 'unsigned');
        } catch (e) {
          if (e instanceof TypeError || e instanceof ValueError)
            throw new ValueError('degree must be convertible to long');
          throw e;
        }
      } else {
        const invalid = () =>
          new ValueError(
            `degree argument must be a nonnegative integer, got ${polynomialScalarRepr(degree)}`
          );
        if (polynomialScalarSign(degree, '<') === -1) throw invalid();
        const index = polynomialInteger(degree, 'unsigned');
        if (
          (degree instanceof Rational && degree.denominator !== 1n) ||
          (typeof degree === 'number' && !Number.isInteger(degree))
        )
          throw invalid();
        length =
          backend === 'integer' || backend === 'word' || backend === 'extension'
            ? BigInt.asUintN(64, index + 1n)
            : index + 1n;
        if (backend === 'generic' || backend === 'binary' || backend === 'large') {
          const incremented = polynomialScalarAddOne(degree);
          if (incremented !== null) {
            length = polynomialInteger(incremented, 'long');
            // Sage compares len(v) in the scalar's parent, and computes the
            // padding count there too. Reduction can change the branch.
            const scalar = degree as { parent: { __call__(x: bigint): unknown } };
            const reducedLength = polynomialInteger(
              scalar.parent.__call__(BigInt(this.coeffs.length)),
              'long'
            );
            if (degree instanceof FiniteFieldElement && reducedLength !== length) {
              if (reducedLength < length)
                throw new TypeError(
                  `can't multiply sequence by non-int of type '${polynomialScalarType(degree)}'`
                );
              // Slicing requires __index__, which PARI extension elements lack.
              polynomialInteger(incremented, 'index');
            }
            if (reducedLength <= length)
              length = BigInt(this.coeffs.length) + length - reducedLength;
          }
        }
        if (
          (backend === 'generic' || backend === 'binary' || backend === 'large') &&
          degree instanceof Rational &&
          length > BigInt(this.coeffs.length)
        )
          throw new TypeError(
            "unsupported operand parent(s) for *: '<class 'list'>' and 'Rational Field'"
          );
      }
    }
    if (length - BigInt(this.coeffs.length) >= 1n << 63n)
      throw new OverflowError("cannot fit 'int' into an index-sized integer");
    const count = Number(length);
    const result = Array.from(
      { length: count },
      (_, i) => this.coeffs[count - i - 1] ?? this.parent.base_ring.zero()
    );
    return new Polynomial(result, this.parent);
  }

  /**
   * Return the resultant of this polynomial and other.
   *
   * The resultant of two polynomials f and g is the determinant of their
   * Sylvester matrix. It is zero if and only if f and g have a common root.
   *
   * @param other - A polynomial or scalar with a common canonical parent
   * @returns The resultant (an element of the base ring)
   *
   * @example
   * ```typescript
   * // If f = x^3 + x + 1 and g = x^3 - x - 1
   * // f.resultant(g) = -8
   * ```
   *
   * @see Reference: sage/rings/polynomial/polynomial_element.pyx:resultant
   * @see Deviation: Polynomial Resultant Delegation and Real Double Coefficients
   */
  resultant(other: Polynomial<C>, options?: { proof?: unknown }): C;
  resultant<D extends RingElement>(other: Polynomial<D>, options?: { proof?: unknown }): C | D;
  resultant(other: unknown, options?: { proof?: unknown }): RingElement;
  resultant(operand: unknown, options?: { proof?: unknown }): RingElement {
    const backend = polynomialBackend(this.parent.base_ring);
    const proof = options !== undefined && Object.hasOwn(options, 'proof');
    if (backend === 'extension') {
      if (proof) throw new TypeError("resultant() got an unexpected keyword argument 'proof'");
      const parent = polynomialElementParent(operand);
      if (!polynomialCanCoerce(this.parent, parent))
        throw new TypeError(`no canonical coercion from ${parent} to ${this.parent}`);
    }
    if (!(operand instanceof Polynomial)) {
      if (
        typeof operand === 'bigint' ||
        operand instanceof Integer ||
        typeof operand === 'boolean' ||
        (typeof operand === 'number' && Number.isInteger(operand))
      )
        return this.resultant(
          this.parent.__call__(operand instanceof Integer ? operand.value : BigInt(operand)),
          options
        );
      let parent: CoefficientRing<RingElement>;
      if (typeof operand === 'number' && polynomialRealBase(this.parent.base_ring)) {
        const scalar = new Polynomial(
          [RDF.__call__(operand)],
          new PolynomialRing(RDF, this.parent.variable_name)
        );
        const [a, b] = polynomialCommonOperands<RingElement>(this, scalar, 'resultant');
        return a.resultant(b, options);
      }
      try {
        parent = polynomialElementParent(operand);
      } catch (e) {
        if (!(e instanceof AttributeError)) throw e;
        throw new TypeError(
          `no common canonical parent for objects with parents: '${this.parent}' and '<class '${polynomialScalarType(operand)}'>'`
        );
      }
      if (!polynomialCommonBase(this.parent, parent))
        throw new TypeError(
          `no common canonical parent for objects with parents: '${this.parent}' and '${parent}'`
        );
      const scalar = new Polynomial(
        [operand as RingElement],
        new PolynomialRing(parent, this.parent.variable_name)
      );
      const [a, b] = polynomialCommonOperands<RingElement>(this, scalar, 'resultant');
      return a.resultant(b, options);
    }
    if (this.parent !== operand.parent) {
      const [a, b] = polynomialCommonOperands(this, operand, 'resultant');
      return a.resultant(b, options);
    }
    const other = operand as Polynomial<C>;
    const base = this.parent.base_ring;
    if (proof && backend !== 'integer')
      throw new TypeError("resultant() got an unexpected keyword argument 'proof'");
    if (polynomialRdfBase(base)) {
      const realA = this.coeffs.map(polynomialRealCoefficient);
      const realB = other.coeffs.map(polynomialRealCoefficient);
      if (realA.some((c) => c === null) || realB.some((c) => c === null))
        throw new NotImplementedError(
          'SAGE_NOT_IMPLEMENTED: resultant with nonconstant real polynomial coefficients'
        );
      try {
        return base.__call__(RDF.__call__(pariResultant(realA as number[], realB as number[])));
      } catch (e) {
        if (!(e instanceof PariError)) throw e;
        if (this.isZero() || other.isZero())
          throw new ValueError('The Sylvester matrix is not defined for zero polynomials');
        if (this.isConstant() && other.isConstant()) return base.one();
        if (base === (RDF as unknown)) {
          if ([...realA, ...realB].some((c) => !Number.isFinite(c)))
            throw new ValueError('array must not contain infs or NaNs');
          throw new NotImplementedError('SAGE_NOT_IMPLEMENTED: RDF resultant overflow fallback');
        }
        // A polynomial coefficient ring uses generic determinant arithmetic,
        // which retains NaNs rather than SciPy's finite-entry validation.
      }
    }
    if (backend === 'integer')
      return base.__call__(
        _fmpz_poly_resultant(extractIntegerCoeffs(this), extractIntegerCoeffs(other))
      );
    if (backend === 'rational') {
      const [a, denA] = polynomialRationalData(this),
        [b, denB] = polynomialRationalData(other);
      const [num, den] = _fmpq_poly_resultant(a, denA, b, denB);
      return base.__call__(new Rational(num, den));
    }
    if (backend === 'word') {
      const p = getCharacteristic(base);
      const lift = (x: C) => (x as unknown as { value: bigint }).value;
      if (ringIsField(base)) {
        const a = this.coeffs.map(lift),
          b = this === other ? a : other.coeffs.map(lift);
        return base.__call__(_nmod_poly_resultant(a, b, p));
      }
      const rows = this.sylvester_matrix(other).map((row) => row.map(lift));
      return base.__call__(new Matrix_modn_dense(rows.length, rows.length, p, rows).determinant());
    }
    // Handle zero polynomials
    if (this.isZero() || other.isZero()) {
      return this.parent.base_ring.zero() as C;
    }

    const m = this.degree();
    const n = other.degree();

    // If both are constants, resultant is 1 (empty matrix determinant)
    if (m === 0 && n === 0) {
      return this.parent.base_ring.one() as C;
    }

    // Handle constant polynomial cases
    if (m === 0) {
      // Res(c, g) = c^deg(g)
      let result = this.coeffs[0]!;
      for (let i = 1; i < n; i++) {
        result = result.mul(this.coeffs[0]!) as C;
      }
      return result;
    }

    if (n === 0) {
      // Res(f, c) = c^deg(f)
      let result = other.coeffs[0]!;
      for (let i = 1; i < m; i++) {
        result = result.mul(other.coeffs[0]!) as C;
      }
      return result;
    }

    // Build and compute Sylvester matrix determinant
    return matrixDeterminant(polynomialSylvesterEntries(this, other), this.parent.base_ring);
  }

  /**
   * Return the Sylvester matrix of this polynomial and `other`.
   *
   * For `deg(self) = m` and `deg(other) = n` this is the `(m+n) x (m+n)`
   * matrix whose first `n` rows hold the coefficients of `x^i * self` and
   * whose last `m` rows hold the coefficients of `x^i * other`.
   *
   * Mixed/scalar parents are canonically coerced; zero polynomials raise.
   * The optional variable is converted as in Sage and otherwise does not alter
   * the univariate matrix. Entries are returned as a dense array.
   * @see Deviation: Polynomial Sylvester Matrices and Explicit Scalar Construction
   * @see Reference: sage/rings/polynomial/polynomial_element.pyx:sylvester_matrix
   */
  sylvester_matrix(other: Polynomial<C>, variable?: unknown): C[][];
  sylvester_matrix<D extends RingElement>(other: Polynomial<D>, variable?: unknown): (C | D)[][];
  sylvester_matrix(other: unknown, variable?: unknown): RingElement[][];
  sylvester_matrix(other: unknown, variable?: unknown): RingElement[][] {
    const rightParent = polynomialElementParent(other);
    const sameParent =
      other instanceof Polynomial &&
      polynomialCanCoerce(this.parent, rightParent) &&
      polynomialCanCoerce(rightParent, this.parent);
    if (!sameParent) {
      const common = polynomialCommonBase(this.parent, rightParent);
      if (!common)
        throw new TypeError(
          `no common canonical parent for objects with parents: '${this.parent}' and '${rightParent}'`
        );
      // Sage calls self.variables()[0] after canonical coercion and replaces
      // the supplied variable. Constant polynomials have an empty variable tuple.
      let a: Polynomial<RingElement>, b: Polynomial<RingElement>;
      if (other instanceof Polynomial)
        [a, b] = polynomialCommonOperands(this, other, 'sylvester_matrix');
      else if (rightParent === (ZZ as unknown)) {
        a = this;
        b = this.parent.__call__(other instanceof Integer ? other.value : BigInt(other as bigint));
      } else {
        const scalar = new Polynomial(
          [other as RingElement],
          new PolynomialRing(rightParent, this.parent.variable_name)
        );
        [a, b] = polynomialCommonOperands(this, scalar, 'sylvester_matrix');
      }
      if (this.isConstant()) throw new IndexError('tuple index out of range');
      return a.sylvester_matrix(b, a.parent.gen());
    }
    if (polynomialTruth(variable)) {
      const variableParent = polynomialElementParent(variable);
      if (variableParent !== this.parent) this.parent.__call__(variable);
    }
    const rhs = polynomialCommonOperands(this, other as Polynomial<C>, 'sylvester_matrix')[1];
    if (this.isZero() || rhs.isZero())
      throw new ValueError('The Sylvester matrix is not defined for zero polynomials');
    // The reference's default small extension matrix implementations require
    // Givaro's _cache, which its PARI-FFELT parent does not provide.
    const base = this.parent.base_ring;
    if (base.zero() instanceof FiniteFieldElement) {
      const q = getFieldOrder(base),
        p = getCharacteristic(base);
      if (q > p && (p === 2n ? q <= 65536n && this.degree() + rhs.degree() > 0 : q < 256n))
        throw new AttributeError(
          "'FiniteField_pari_ffelt_with_category' object has no attribute '_cache'"
        );
    }
    return polynomialSylvesterEntries(this, rhs);
  }

  /**
   * Return the discriminant of this polynomial.
   *
   * The discriminant is defined as:
   *   disc(f) = (-1)^(n(n-1)/2) * Res(f, f') / a_n
   *
   * where n is the degree, a_n is the leading coefficient, and f' is the derivative.
   *
   * The discriminant is zero if and only if the polynomial has a repeated root.
   *
   * @returns The discriminant (an element of the base ring)
   *
   * @example
   * ```typescript
   * // If f = x^3 + x + 1
   * // f.discriminant() = -31
   * ```
   *
   * @see Reference: sage/rings/polynomial/polynomial_element.pyx:discriminant
   */
  discriminant(): C {
    if (this.isZero()) {
      return this.parent.base_ring.zero() as C;
    }

    const n = this.degree();
    if (n <= 0) {
      // Constant polynomial has discriminant 1 (or the leading coefficient for degree 0)
      return this.parent.base_ring.one() as C;
    }

    const d = this.derivative();
    const k = d.degree();
    const an = this.leading_coefficient();

    // Compute sign: (-1)^(n*(n-1)/2)
    // n*(n-1)/2 mod 2:
    //   n=0: 0 -> +1
    //   n=1: 0 -> +1
    //   n=2: 1 -> -1
    //   n=3: 3 -> -1
    //   n=4: 6 -> +1
    //   n=5: 10 -> +1
    // Pattern: sign is -1 when n mod 4 is 2 or 3
    const r = n % 4;
    const signIsNegative = r === 2 || r === 3;

    // Compute a_n^(n - k - 2) where k = deg(f')
    // Normally k = n - 1, so n - k - 2 = n - (n-1) - 2 = -1
    // This means we need to divide by a_n
    const exponent = n - k - 2;

    // Compute resultant
    const res = this.resultant(d);

    let result: C;
    if (exponent >= 0) {
      // Multiply by a_n^exponent
      let anPower = this.parent.base_ring.one() as C;
      for (let i = 0; i < exponent; i++) {
        anPower = anPower.mul(an) as C;
      }
      result = res.mul(anPower) as C;
    } else {
      // Divide by a_n^(-exponent)
      let anPower = an;
      for (let i = 1; i < -exponent; i++) {
        anPower = anPower.mul(an) as C;
      }
      const quotient = divideCoeffs(res, anPower);
      if (quotient.mul(anPower).eq(res)) {
        result = quotient;
      } else {
        // Division by the leading coefficient is not exact in the base ring.
        // Rather than dividing the resultant, alter the Sylvester matrix
        // (Sage issue #11782, `polynomial_element.pyx:8094-8099`).
        if (exponent !== -1) {
          throw new ArithmeticError('discriminant: division by the leading coefficient failed');
        }
        const mat = polynomialSylvesterEntries(this, d);
        mat[0]![0] = this.parent.base_ring.one() as C;
        mat[n - 1]![0] = mulByInteger(this.parent.base_ring.one() as C, n, this.parent.base_ring);
        result = matrixDeterminant(mat, this.parent.base_ring);
      }
    }

    // Apply sign
    if (signIsNegative) {
      result = result.neg() as C;
    }

    return result;
  }

  /**
   * Return the roots of this polynomial in the base ring.
   *
   * ZZ uses Sage's degree-100 dense/sparse dispatch and NTL/PARI factorization;
   * QQ uses PARI factorization. Sparse ZZ root order is preserved.
   * Finite fields extract linear roots in the order returned by factorization.
   * Over IntegerModRing, default roots have the cached GF parent; distinct roots
   * retain the residue-ring parent and use its CRT/Hensel method.
   *
   * Set multiplicities:false to return distinct roots. Over finite fields this
   * first computes gcd(self, x^q-x), as in Sage, then factors that polynomial.
   *
   * @returns Root/multiplicity pairs by default, or a list of roots when false.
   *
   * @example
   * ```typescript
   * // Over GF(7): x^2 - 1 = (x-1)(x+1) has roots 1 and 6
   * const p = x.pow(2).sub(R.one());
   * const roots = p.roots(); // [[6, 1], [1, 1]]
   * ```
   *
   * @see Reference: sage/rings/polynomial/polynomial_element.pyx:roots
   * @see Deviation: Polynomial Roots and Factorization
   * @see Deviation: Modular Polynomial Roots and Hensel Lifting
   */
  roots(this: Polynomial<IntegerMod & RingElement>, options?: { multiplicities?: true }): Array<[LegacyPrimeElement, number]>;
  roots(this: Polynomial<IntegerMod & RingElement>, options: { multiplicities: false }): IntegerMod[];
  roots(this: Polynomial<IntegerMod & RingElement>, options: { multiplicities?: boolean }): Array<[LegacyPrimeElement, number]> | IntegerMod[];
  roots(options?: { multiplicities?: true }): Array<[C, number]>;
  roots(options: { multiplicities: false }): C[];
  roots(options: { multiplicities?: boolean }): Array<[C, number]> | C[];
  roots(options: { multiplicities?: boolean } = {}): unknown[] {
    if (this.parent.base_ring instanceof IntegerModRing) {
      return this.parent.base_ring._roots_univariate_polynomial(
        this as unknown as Polynomial<IntegerMod & RingElement>, options);
    }
    if (options.multiplicities === false) {
      const base = this.parent.base_ring;
      if (!isIntegerRing(base) && !isRationalField(base) && isFiniteField(base)) {
        if (this.isZero()) {
          // Preserve the bundled modular-power wrappers' distinct zero errors.
          const backend = polynomialBackend(base);
          if (backend === 'word') throw new ZeroDivisionError('');
          if (backend === 'large') throw new NTLError('ZZ_pX: division by zero');
          throw new ZeroDivisionError('modulus must be nonzero');
        }
        if (this.degree() === 0) return [];
        // finite_field_base: remove all nonlinear factors and multiplicities
        // before factoring. In particular, the factorization degree is <= q.
        const x = this.parent.gen();
        const g = this.gcd(powerMod(x, getFieldOrder(base), this).sub(x));
        return g.factor().filter(([f]) => f.degree() === 1).map(([f]) =>
          base.__call__(divideCoeffs(f.getCoeff(0).neg() as C, f.getCoeff(1))) as C);
      }
      return this.roots().map(([root]) => root);
    }
    if (this.isZero()) {
      const ring = this.parent.base_ring;
      if (isIntegerRing(ring)) throw new ValueError('roots of 0 are not defined');
      if (isRationalField(ring))
        throw new NotImplementedError('root finding for this polynomial not implemented');
      if (isFiniteField(ring)) this.factor(); // Preserve the field backend's zero error.
      throw new ValueError('roots of zero polynomial are not defined');
    }

    if (this.degree() === 0) {
      return []; // Constant non-zero polynomial has no roots
    }

    // Check if base ring is a finite field (has cardinality method and is iterable)
    const baseRing = this.parent.base_ring;

    // Handle integer polynomials (ZZ[x])
    if (isIntegerRing(baseRing)) {
      const coeffs = extractIntegerCoeffs(this);
      const intRoots = findIntegerRoots(coeffs);
      // Convert back to ring elements
      return intRoots.map(([root, mult]): [C, number] => [baseRing.__call__(root) as C, mult]);
    }

    // Handle rational polynomials (QQ[x]) - find rational roots
    if (isRationalField(baseRing)) {
      // QQ factors through PARI, then makes factors monic before sorting.
      const rationalRoots = findRationalRoots(this);
      return sortRootsSageOrder(rationalRoots);
    }

    // Check if it's a finite field we can iterate over
    if (!isFiniteField(baseRing)) {
      throw new NotImplementedError('roots only implemented for finite fields, ZZ, and QQ');
    }

    const roots: Array<[C, number]> = [];

    // finite_field_base._roots_univariate_polynomial uses factorization even
    // for small fields. Factor order compares -r in the coefficient field;
    // sorting numeric root representatives would move zero to the wrong end.
    for (const [fac, mult] of this.factor()) {
      if (fac.degree() === 1) {
        const root = divideCoeffs(fac.getCoeff(0).neg() as C, fac.getCoeff(1));
        roots.push([baseRing.__call__(root) as C, mult]);
      }
    }
    return roots;
  }

  /**
   * Factor this polynomial over the integers ZZ.
   * Returns pairs [irreducible_factor, multiplicity].
   *
   * The content is returned as its own prime factors, as Sage does
   * (`(12*(x^2+1)^3*(x+2)).factor()` is `2^2 * 3 * (x + 2) * (x^2 + 1)^3`),
   * and a negative content contributes the unit `-1`.
   *
   * @see Deviation: Polynomial Roots and Factorization
   */
  private _factorOverIntegers(): Array<[Polynomial<C>, number]> {
    const coeffs = extractIntegerCoeffs(this);
    const positiveContent = intPolyContent(coeffs);
    const content = coeffs[coeffs.length - 1]! < 0n ? -positiveContent : positiveContent;
    const result: Array<[Polynomial<C>, number]> = [];

    // The sign of the content is the unit of the factorization; Sage keeps it
    // in ``Factorization.unit()`` (e.g. ``(-x^2+4).factor() == (-1)*(x-2)*(x+2)``).
    if (content < 0n) {
      result.push([new Polynomial([this.parent.base_ring.__call__(-1n) as C], this.parent), 1]);
    }

    // Add content as a factor if it's not 1 or -1
    if (content !== 1n && content !== -1n) {
      // Factor the integer content
      const intFactors = factorInteger(content < 0n ? -content : content);
      for (const [p, e] of intFactors) {
        if (p > 1n) {
          const constPoly = new Polynomial([this.parent.base_ring.__call__(p) as C], this.parent);
          result.push([constPoly, Number(e)]);
        }
      }
    }

    // Sage factors content before calling its degree-selected polynomial backend.
    const primitive = coeffs.map((c) => c / positiveContent);
    const degree = this.degree();
    const factors =
      degree < 30 || degree > 300 ? ntlIntegerFactor(primitive)[1] : pariIntegerFactor(primitive);

    // Convert polynomial factors back to Polynomial<C>
    for (const [facCoeffs, mult] of factors) {
      const polyCoeffs = facCoeffs.map((c) => this.parent.base_ring.__call__(c) as C);
      const poly = new Polynomial(polyCoeffs, this.parent);
      result.push([poly, mult]);
    }

    // Factorization.sort: degree, multiplicity, then polynomial comparison.
    result.sort(comparePolynomialFactors);

    return result;
  }

  /**
   * Factor this polynomial over the rationals QQ.
   * Returns pairs [monic_irreducible_factor, multiplicity].
   *
   * Sage keeps the leading coefficient in `Factorization.unit()`; we return it
   * as a degree-0 factor so that the product of the returned factors is again
   * `self`.
   *
   * @see Deviation: Polynomial Roots and Factorization
   */
  private _factorOverRationals(): Array<[Polynomial<C>, number]> {
    // QQ._factor_univariate_polynomial always delegates to PARI, including
    // degrees for which integer-polynomial factorization selects NTL.
    const factors = pariRationalFactor(clearDenominators(this));

    const result: Array<[Polynomial<C>, number]> = [];

    // Convert polynomial factors back to monic Polynomial<C> over QQ.
    // The quotient is formed with the base ring's own division so that any
    // rational-like coefficient ring works, not only ones whose `__call__`
    // understands a numerator/denominator pair.
    for (const [facCoeffs, mult] of factors) {
      // Make monic by dividing by leading coefficient
      const lc = this.parent.base_ring.__call__(facCoeffs[facCoeffs.length - 1]!) as C;
      const monicCoeffs = facCoeffs.map((c) =>
        divideCoeffs(this.parent.base_ring.__call__(c) as C, lc)
      );
      const poly = new Polynomial(monicCoeffs, this.parent);
      result.push([poly, mult]);
    }

    // All factors are monic, so the leading coefficient of ``self`` is the
    // unit of the factorization (Sage keeps it in ``Factorization.unit()``).
    const unit = this.leading_coefficient();
    if (!unit.eq(1)) {
      result.push([new Polynomial([unit], this.parent), 1]);
    }

    // Factorization.sort: degree, multiplicity, then polynomial comparison.
    result.sort(comparePolynomialFactors);

    return result;
  }

  /** Square test and optional root via Sage's squarefree-decomposition algorithm.
   * @see Reference: sage/rings/polynomial/polynomial_element.pyx:1995
   */
  is_square(root?: false): boolean;
  is_square(root: true): [boolean, Polynomial<C> | null];
  is_square(root: boolean): boolean | [boolean, Polynomial<C> | null];
  is_square(root = false): boolean | [boolean, Polynomial<C> | null] {
    if (this.isZero()) return root ? [true, this] : true;
    let factors: Array<[Polynomial<C>, number]>;
    try { factors = this.squarefree_decomposition(); }
    catch (e) {
      if (!(e instanceof NotImplementedError)) throw e;
      factors = this.factor();
    }
    const unit = factors.filter(([f]) => f.degree() === 0).reduce((c, [f, e]) => c.mul(polynomialCoefficientPower(f.getCoeff(0), e, this.parent.base_ring) as C) as C, this.parent.base_ring.one());
    factors = factors.filter(([f]) => f.degree() > 0);
    const u = unit as C & {
      is_square?: () => boolean; sqrt?: () => C; isqrt?: () => C;
    };
    if (factors.some(([, e]) => e % 2 !== 0)) return root ? [false, null] : false;
    if (!u.is_square) throw new NotImplementedError('SAGE_NOT_IMPLEMENTED: coefficient is_square');
    if (!u.is_square()) return root ? [false, null] : false;
    const c = u.sqrt ? u.sqrt() : u.isqrt ? u.isqrt() : undefined;
    if (c === undefined) throw new NotImplementedError('SAGE_NOT_IMPLEMENTED: coefficient sqrt');
    let g = new Polynomial([c], this.parent);
    for (const [f, e] of factors) g = g.mul(f.pow(e / 2) as Polynomial<C>);
    return root ? [true, g] : true;
  }

  /**
   * Return the squarefree decomposition of this polynomial.
   *
   * Returns a list of pairs (f_i, i) where this polynomial equals
   * prod(f_i^i) and each f_i is squarefree and coprime to the others.
   *
   * Uses the standard algorithm based on gcd with the derivative.
   *
   * @returns Array of [squarefree_factor, multiplicity] pairs
   *
   * @example
   * ```typescript
   * // (x-1)^2 * (x-2) has squarefree decomposition [(x-2, 1), ((x-1), 2)]
   * ```
   *
   * @see Reference: sage/rings/polynomial/polynomial_element.pyx:squarefree_decomposition
   */
  squarefree_decomposition(): Array<[Polynomial<C>, number]> {
    const base = this.parent.base_ring, backend = polynomialBackend(base);
    if (backend === 'integer') {
      if (this.degree() <= 0) return [[this, 1]];
      const a = this.coeffs.map(c => BigInt(String(c))), content = _ZZX_kernels.content(a);
      const factors = ZZX_SquareFreeDecomp(a.map(c => c / content)).map(([f, e]) =>
        [new Polynomial(f.map(c => base.__call__(c)), this.parent), e] as [Polynomial<C>, number]);
      if (content !== 1n) factors.unshift([new Polynomial([base.__call__(content)], this.parent), 1]);
      return factors;
    }
    if (this.isZero()) {
      if (backend === 'word') throw new ArithmeticError('square-free decomposition of 0 is not defined');
      throw new ValueError('square-free decomposition not defined for zero polynomial');
    }
    if (backend === 'word') {
      const p = getRingCharacteristic(base)!;
      if (!ringIsField(base)) throw new NotImplementedError('square free factorization of polynomials over rings with composite characteristic is not implemented');
      const factors = nmod_poly_factor_squarefree(this.coeffs.map(c => BigInt(String(c))), p).map(([f, e]) =>
        [new Polynomial(f.map(c => base.__call__(c)), this.parent), e] as [Polynomial<C>, number]);
      const unit = this.leading_coefficient();
      if (!unit.eq(base.one())) factors.unshift([new Polynomial([unit], this.parent), 1]);
      return factors;
    }
    // IntegerModRing does not provide the finite-field base-class method.
    if (base.zero() instanceof IntegerMod) throw new NotImplementedError('square-free decomposition not implemented for this polynomial');
    if (getRingCharacteristic(base) === 0n && ringIsField(base)) return fieldSquarefree(this);
    if (isFiniteField(base)) return finiteSquarefree(this, getRingCharacteristic(base)!, getFieldOrder(base));
    throw new NotImplementedError('square-free decomposition not implemented for this polynomial');
  }

  /**
   * Return the distinct-degree factorization of this polynomial.
   *
   * Assumes the polynomial is squarefree. Returns a list of pairs (g_d, d)
   * where g_d is the product of all irreducible factors of degree d.
   *
   * This is a key step in the Berlekamp and Cantor-Zassenhaus algorithms.
   *
   * @returns Array of [product_of_degree_d_factors, d] pairs
   *
   * @example
   * ```typescript
   * // Over GF(163), factoring (x+162)(x^3+7x+161)(x^7+9x+161)
   * // gives [(x+162, 1), (x^3+7x+161, 3), (x^7+9x+161, 7)]
   * ```
   *
   * @see Reference: sage/rings/polynomial/polynomial_element.pyx:_distinct_degree_factorisation_squarefree
   */
  distinct_degree_factorization(): Array<[Polynomial<C>, number]> {
    if (this.isZero()) {
      throw new ValueError('distinct-degree factorization of zero polynomial is not defined');
    }

    const baseRing = this.parent.base_ring;

    if (!isFiniteField(baseRing)) {
      throw new NotImplementedError(
        'distinct-degree factorization only implemented for finite fields'
      );
    }

    const q = getFieldOrder(baseRing);
    const result: Array<[Polynomial<C>, number]> = [];

    // x
    const x = this.parent.gen();
    // Work with monic polynomial
    let v = this._monic();
    // w = x mod v
    let w = x.mod(v);
    let d = 0;
    let e = v.degree();

    // Iterate over all possible degrees
    while (2 * (d + 1) <= e) {
      d = d + 1;

      // w = w^q mod v (computing x^{q^d} mod v)
      w = powerMod(w, q, v);

      // a_d = gcd(v, w - x)
      const wMinusX = w.sub(x);
      const ad = v.gcd(wMinusX);

      if (!ad.eq(this.parent.one())) {
        result.push([ad._monic(), d]);
        // v = v / a_d
        v = v.quo_rem(ad)[0];
      }

      e = v.degree();
    }

    // If v still has positive degree, it's irreducible of degree e
    if (e > 0) {
      result.push([v._monic(), e]);
    }

    return result;
  }

  /**
   * Return the factorization of this polynomial.
   *
   * For finite fields, uses:
   * 1. Squarefree decomposition
   * 2. Distinct-degree factorization
   * 3. Cantor-Zassenhaus algorithm for equal-degree factorization
   *
   * @returns Factorization as list of [factor, multiplicity] pairs
   *
   * @example
   * ```typescript
   * // x^4 - 1 over GF(5) factors as (x-1)(x+1)(x-2)(x+2)
   * const factors = (x.pow(4).sub(R.one())).factor();
   * ```
   *
   * @see Reference: sage/rings/polynomial/polynomial_element.pyx:factor
   * @see Deviation: Polynomial Roots and Factorization
   */
  factor(): Array<[Polynomial<C>, number]> {
    const baseRing = this.parent.base_ring;
    if (this.isZero()) {
      // Dense ZZ factor() divides by its zero content before calling a backend.
      if (isIntegerRing(baseRing)) throw new ZeroDivisionError('division by zero');
      throw new ArithmeticError('factorization of 0 is not defined');
    }

    // Integer constants have prime factors too; dispatch before the field-unit case.
    if (isIntegerRing(baseRing)) return this._factorOverIntegers();
    if (this.degree() === 0) return this.coeffs[0]!.eq(baseRing.one() as C) ? [] : [[this, 1]];

    // Handle rational polynomials (QQ[x])
    if (isRationalField(baseRing)) {
      return this._factorOverRationals();
    }

    if (!isFiniteField(baseRing)) {
      const characteristic = getRingCharacteristic(baseRing);
      if (characteristic !== null && characteristic > 0n && !is_prime(characteristic))
        throw new NotImplementedError(
          'factorization of polynomials over rings with composite characteristic is not implemented'
        );
      throw new NotImplementedError('factorization only implemented for finite fields, ZZ, and QQ');
    }

    const result: Array<[Polynomial<C>, number]> = [];

    // polynomial_zmod_flint.factor -> nmod_poly_linkage.factor_helper.
    if (polynomialBackend(baseRing) === 'word') {
      const p = getRingCharacteristic(baseRing)!;
      const [unit, factors] = nmod_poly_factor(this.coeffs.map(c => BigInt(String(c))), p);
      for (const [f, e] of factors)
        result.push([new Polynomial(f.map(c => baseRing.__call__(c)), this.parent), e]);
      if (unit !== 1n) result.push([this.parent.__call__(unit), 1]);
      result.sort(comparePolynomialFactors);
      return result;
    }

    // GF(2) inherits polynomial_element.factor, which delegates to PARI.
    if (polynomialBackend(baseRing) === 'binary') {
      const bits = this.coeffs.reduce((v, c, i) => v | (BigInt(String(c)) << BigInt(i)), 0n);
      for (const [f, e] of F2x_factor(bits)) {
        const coefficients = Array.from({ length: f.toString(2).length }, (_, i) =>
          baseRing.__call__((f >> BigInt(i)) & 1n));
        result.push([new Polynomial(coefficients, this.parent), e]);
      }
      result.sort(comparePolynomialFactors);
      return result;
    }

    // polynomial_element.factor -> PARI for the generic large-prime backend.
    if (polynomialBackend(baseRing) === 'large') {
      const p = getRingCharacteristic(baseRing)!;
      for (const [f, e] of pariFpXFactor(this.coeffs.map((c) => BigInt(String(c))), p))
        result.push([new Polynomial(f.map((c) => baseRing.__call__(c)), this.parent), e]);
      const unit = this.leading_coefficient();
      if (!unit.eq(1)) result.push([this.parent.__call__(unit), 1]);
      result.sort(comparePolynomialFactors);
      return result;
    }

    // Step 1: Squarefree decomposition
    const sqfree = baseRing.zero() instanceof IntegerMod && polynomialBackend(baseRing) !== 'word'
      ? finiteSquarefree(this, getRingCharacteristic(baseRing)!, getFieldOrder(baseRing))
      : this.squarefree_decomposition();

    // Step 2: For each squarefree factor, do distinct-degree and equal-degree factorization
    for (const [sqfFactor, mult] of sqfree) {
      if (sqfFactor.degree() === 0) {
        // Constant factor
        if (!sqfFactor.leading_coefficient().eq(1)) {
          result.push([sqfFactor, mult]);
        }
        continue;
      }

      // Distinct-degree factorization
      const ddf = sqfFactor.distinct_degree_factorization();

      // Equal-degree factorization for each degree
      for (const [ddFactor, degree] of ddf) {
        if (ddFactor.degree() === degree) {
          // ddFactor is already irreducible
          result.push([ddFactor, mult]);
        } else {
          // Need to split ddFactor into irreducible factors of the given degree
          const irreducibles = cantorZassenhausFactorization(ddFactor, degree, this.parent);

          for (const irr of irreducibles) {
            result.push([irr, mult]);
          }
        }
      }
    }

    // Factorization.sort: degree, multiplicity, then polynomial comparison.
    result.sort(comparePolynomialFactors);

    return result;
  }

  /**
   * Test if this polynomial is irreducible.
   *
   * Follows Sage: the zero polynomial and units are reducible, a constant is
   * irreducible iff it is irreducible in the base ring, and otherwise the
   * polynomial is factored (over finite fields we use Rabin's test, which is
   * what FLINT's `nmod_poly_is_irreducible_rabin` does).
   *
   * @returns true if irreducible
   *
   * @example
   * ```typescript
   * // x^2 + x + 1 is irreducible over GF(2)
   * const p = x.pow(2).add(x).add(R.one());
   * p.is_irreducible(); // true
   * ```
   *
   * @see Reference: sage/rings/polynomial/polynomial_element.pyx:is_irreducible
   * @see Deviation: Polynomial Roots and Factorization
   */
  is_irreducible(): boolean {
    const ring = this.parent.base_ring;
    // This adapter implements the ZZ/QQ cached methods. Finite-field methods
    // have separate backend/algorithm contracts and remain on their existing path.
    if (!isIntegerRing(ring) && !isRationalField(ring)) return _is_irreducible_uncached(this);
    const known = polynomialIrreducibility.get(this);
    if (known !== undefined) return known;
    const result = _is_irreducible_uncached(this);
    polynomialIrreducibility.set(this, result);
    return result;
  }


  /**
   * String representation.
   * @see Deviation: Polynomial Common Coefficient Parents and Representation
   */
  toString(): string {
    if (polynomialGenerators.has(this) && polynomialBackend(this.parent.base_ring) !== 'integer')
      return this.parent.variable_name;
    if (this.coeffs.length === 0) {
      return '0';
    }

    const varName = this.parent.variable_name;
    const terms: string[] = [];

    for (let i = this.coeffs.length - 1; i >= 0; i--) {
      const c = this.coeffs[i]!;
      if (c.isZero()) {
        continue;
      }

      const coefficient = c.toString();
      const term =
        i === 0
          ? coefficient
          : `${needsParens(coefficient) ? `(${coefficient})` : coefficient}*${varName}${i > 1 ? `^${i}` : ''}`;

      terms.push(term);
    }

    if (terms.length === 0) {
      return '0';
    }

    // polynomial_element.pyx:_repr performs sign/unit simplification after
    // adding terms, preserving parentheses around compound coefficients.
    return (' ' + (this.coeffs[this.coeffs.length - 1]!.isZero() ? ' + ' : '') + terms.join(' + '))
      .replaceAll(' + -', ' - ')
      .replace(/ 1(\.0+)?\*/g, ' ')
      .replace(/ -1(\.0+)?\*/g, ' -')
      .slice(1);
  }
}

/** Underlying predicate; the ZZ/QQ wrapper caches only successful boolean results. */
function _is_irreducible_uncached<C extends RingElement>(poly: Polynomial<C>): boolean {
  if (poly.isZero()) {
    return false;
  }

  const baseRing = poly.parent.base_ring;
  const n = poly.degree();

  if (n === 0) {
    // Sage: ``if self.is_unit(): return False`` then defers to the base
    // ring, so ZZ(5) is irreducible while ZZ(4), ZZ(1) and any nonzero
    // element of a field are not.
    const c = poly.coeffs[0]!;
    if (isIntegerRing(baseRing)) {
      const v = extractIntegerCoeffs(poly)[0]!;
      const a = v < 0n ? -v : v;
      return a > 1n && is_prime(a);
    }
    if (
      'is_irreducible' in c &&
      typeof (c as unknown as { is_irreducible: () => boolean }).is_irreducible === 'function'
    ) {
      return (c as unknown as { is_irreducible: () => boolean }).is_irreducible();
    }
    // Every nonzero constant is a unit over a field.
    return false;
  }

  // Generic ZZ irreducibility uses the full factorization, including content.
  // Only the separately stored unit is ignored in our array representation.
  if (isIntegerRing(baseRing)) {
    const factors = poly.factor().filter(
      ([g]) =>
        g.degree() !== 0 ||
        (!g.getCoeff(0).eq(baseRing.one() as C) && !g.getCoeff(0).eq(baseRing.__call__(-1n) as C))
    );
    return factors.length === 1 && factors[0]![1] === 1;
  }

  // QQ's specialized method factors its positive-leading primitive numerator
  // over ZZ; it does not use QQ.factor()'s unconditional PARI route.
  if (isRationalField(baseRing)) {
    if (n === 1) return true;
    const primitive = intPolyPrimitive(clearDenominators(poly)[0])[1];
    const factors =
      n < 30 || n > 300 ? ntlIntegerFactor(primitive)[1] : pariIntegerFactor(primitive);
    return factors.length === 1 && factors[0]![1] === 1;
  }

  if (!isFiniteField(baseRing)) {
    throw new NotImplementedError(
      'is_irreducible only implemented for finite fields, ZZ, and QQ'
    );
  }

  if (n === 1) {
    return true; // Linear polynomials are irreducible over a field
  }

  // Rabin's irreducibility test over GF(q) (FLINT
  // `nmod_poly_factor/is_irreducible.c:nmod_poly_is_irreducible_rabin`):
  // f of degree n is irreducible iff x^(q^n) = x mod f and
  // gcd(x^(q^(n/l)) - x, f) = 1 for every prime l | n.
  const q = getFieldOrder(baseRing);
  const x = poly.parent.gen();
  const monic = poly._monic();

  // x^(q^n) mod f
  const xqn = powerModIterated(x, q, n, monic);
  if (!xqn.eq(x)) {
    return false;
  }

  for (const [l] of factorInteger(BigInt(n))) {
    if (l <= 1n) continue;
    const a = powerModIterated(x, q, n / Number(l), monic).sub(x);
    if (a.isZero()) {
      return false;
    }
    const g = monic.gcd(a);
    if (g.degree() > 0) {
      return false;
    }
  }

  return true;
}

/**
 * Check if a coefficient string needs parentheses when multiplied.
 */
function needsParens(s: string): boolean {
  const magnitude = s.startsWith('-') ? s.slice(1) : s;
  if (magnitude === '+infinity' || magnitude === 'infinity') return false;
  const expression = magnitude.replace(/(\d|\.)[eE][+-]\d+/g, '$1');
  return expression.includes('+') || expression.includes('-');
}

/**
 * Sort QQ roots by multiplicity and the exact rational constant of x-r.
 * This reverses rational root order, while finite-field roots retain factor
 * order directly because modular negation does not reverse integer order.
 */
function sortRootsSageOrder<C extends RingElement>(roots: Array<[C, number]>): Array<[C, number]> {
  const keyed = roots.map((entry, idx) => ({ entry, idx, key: exactRootKey(entry[0]) }));
  keyed.sort((a, b) => {
    if (a.entry[1] !== b.entry[1]) {
      return a.entry[1] - b.entry[1];
    }
    if (a.key === null || b.key === null) {
      return a.idx - b.idx;
    }
    // Compare `-r` ascending, i.e. `r` descending.
    const lhs = -a.key.n * b.key.d;
    const rhs = -b.key.n * a.key.d;
    if (lhs < rhs) {
      return -1;
    }
    if (lhs > rhs) {
      return 1;
    }
    return a.idx - b.idx;
  });
  return keyed.map((k) => k.entry);
}

/**
 * Sage polynomial comparison starts at the highest coefficient.
 * @see Deviation: PARI finite-field universal comparison
 */
function comparePolynomialFactors<C extends RingElement>(
  a: [Polynomial<C>, number],
  b: [Polynomial<C>, number]
): number {
  const degree = a[0].degree() - b[0].degree();
  if (degree) return degree;
  if (a[1] !== b[1]) return a[1] - b[1];
  for (let i = a[0].degree(); i >= 0; i--) {
    const x = a[0].getCoeff(i),
      y = b[0].getCoeff(i);
    const u = exactRootKey(x),
      v = exactRootKey(y);
    if (u && v) {
      const left = u.n * v.d,
        right = v.n * u.d;
      if (left !== right) return left < right ? -1 : 1;
    } else if (x instanceof FiniteFieldElement && y instanceof FiniteFieldElement) {
      // PARI-backed Sage elements use cmp_universal, whose raw polynomial
      // comparison reads low coefficients first (and packed words for GF(2^n)).
      const value = pariFiniteCompare(
        { type: PariType.t_FFELT, p: x.parent.characteristic, degree: x.parent.degree,
          value: x.coefficients().map(c => c.value), definingPoly: x.parent.modulus.coeffs.map(c => c.value) },
        { type: PariType.t_FFELT, p: y.parent.characteristic, degree: y.parent.degree,
          value: y.coefficients().map(c => c.value), definingPoly: y.parent.modulus.coeffs.map(c => c.value) }
      );
      if (value) return value;
    }
  }
  return 0;
}

/** polynomial_element.pyx:do_schoolbook_product, including ordered scalar leaves. */
function do_schoolbook_product<C extends RingElement>(
  a: readonly C[],
  b: readonly C[],
  precision: number
): C[] {
  if (!a.length || !b.length) return [];
  const fullLength = a.length + b.length - 1,
    n = precision < 0 ? fullLength : Math.min(precision, fullLength);
  if (a.length === 1) return b.slice(0, n).map((c) => a[0]!.mul(c));
  if (b.length === 1) return a.slice(0, n).map((c) => c.mul(b[0]!));
  const out: C[] = [];
  for (let k = 0; k < n; k++) {
    const start = Math.max(0, k - b.length + 1),
      end = Math.min(k, a.length - 1);
    let sum = a[start]!.mul(b[k - start]!);
    for (let i = start + 1; i <= end; i++) sum = sum.add(a[i]!.mul(b[k - i]!));
    out.push(sum);
  }
  return out;
}

/** The exact rational value of a coefficient, or `null` if it has none. */
function exactRootKey(x: unknown): { n: bigint; d: bigint } | null {
  if (x instanceof GF2Element) return { n: BigInt(x.value), d: 1n };
  const v = (x as { value?: unknown }).value;
  if (typeof v === 'bigint') {
    return { n: v, d: 1n };
  }
  const num = (x as { numerator?: unknown }).numerator;
  const den = (x as { denominator?: unknown }).denominator;
  const numV = typeof num === 'function' ? (num as () => unknown).call(x) : num;
  const denV = typeof den === 'function' ? (den as () => unknown).call(x) : den;
  if (typeof numV === 'bigint' && typeof denV === 'bigint' && denV > 0n) {
    return { n: numV, d: denV };
  }
  return null;
}

/**
 * Divide two coefficients. Assumes the ring supports division.
 */
function divideCoeffs<C extends RingElement>(a: C, b: C): C {
  if (a instanceof Polynomial && b instanceof Polynomial) {
    const [quotient, remainder] = a.quo_rem(b);
    if (!remainder.isZero())
      throw new ArithmeticError(
        'division non exact (consider coercing to polynomials over the fraction field)'
      );
    return quotient as unknown as C;
  }
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
 * Return the inverse of `c` when `c` is a unit of the coefficient ring, and
 * `null` otherwise.
 *
 * This is Sage's `inverse_of_unit()`: it must not succeed for a non-unit (over
 * ZZ, `2` has no inverse), so the candidate inverse is verified.  Inexact
 * rings (where `c * c^-1` is only approximately one) keep working because the
 * verification is done with the ring's own equality.
 */
function inverseOfUnit<C extends RingElement>(c: C, ring: CoefficientRing<C>): C | null {
  if (c instanceof RealDoubleElement) return c.inv() as unknown as C;
  if (c instanceof Polynomial) {
    if (c.degree() !== 0) return null;
    const scalar = c.getCoeff(0);
    const inverse = inverseOfUnit(scalar, c.parent.base_ring);
    if (inverse !== null) return c.parent.__call__(inverse) as unknown as C;
    // Rational.div implements inversion even where the wrapper has no inv method.
    try {
      const candidate = divideCoeffs(c.parent.base_ring.one(), scalar);
      return candidate.mul(scalar).eq(c.parent.base_ring.one())
        ? (c.parent.__call__(candidate) as unknown as C)
        : null;
    } catch (e) {
      if (e instanceof ArithmeticError || e instanceof ValueError || e instanceof ZeroDivisionError)
        return null;
      throw e;
    }
  }
  const withInv = c as unknown as { inv?: () => C };
  if (typeof withInv.inv !== 'function') {
    return null;
  }
  let inv: C;
  try {
    inv = withInv.inv();
  } catch {
    // Sage catches ArithmeticError/ValueError from inverse_of_unit here.
    return null;
  }
  if (!c.mul(inv).eq(ring.one())) {
    return null;
  }
  return inv;
}

/**
 * Multiply a ring element by a non-negative integer using double-and-add.
 *
 * This is `n * c` in the base ring, computed with O(log n) additions instead
 * of n-1 of them.
 */
function mulByInteger<C extends RingElement>(coeff: C, n: number, ring: CoefficientRing<C>): C {
  if (n === 0) {
    return ring.zero() as C;
  }
  let k = n;
  let acc: C | null = null;
  let addend = coeff;
  while (k > 0) {
    if (k & 1) {
      acc = acc === null ? addend : (acc.add(addend) as C);
    }
    k >>= 1;
    if (k > 0) {
      addend = addend.add(addend) as C;
    }
  }
  return acc ?? (ring.zero() as C);
}

/**
 * Return whether a coefficient is negative (meaningful only in ordered rings
 * such as ZZ and QQ; always false elsewhere).
 */
function isNegative<C extends RingElement>(c: C): boolean {
  if ('value' in c) {
    const v = (c as unknown as { value: unknown }).value;
    if (typeof v === 'bigint') return v < 0n;
    if (typeof v === 'number') return v < 0;
  }
  if ('numerator' in c) {
    const num = (c as unknown as { numerator: unknown }).numerator;
    if (typeof num === 'bigint') return num < 0n;
  }
  return c.toString().startsWith('-');
}

/**
 * Return whether a coefficient ring is a field.
 */
function ringIsField<C extends RingElement>(ring: CoefficientRing<C>): boolean {
  if (typeof ring.is_field === 'function') {
    return ring.is_field();
  }
  if (isRationalField(ring)) {
    return true;
  }
  if (isIntegerRing(ring)) {
    return false;
  }
  // Fall back on whether elements can be inverted.
  const one = ring.one();
  return 'inv' in one && typeof (one as unknown as { inv: unknown }).inv === 'function';
}

/**
 * Compute GCD of two coefficients. Assumes the ring supports a gcd method.
 */
function gcdCoeffs<C extends RingElement>(a: C, b: C): C {
  // Try to call gcd method if it exists
  if ('gcd' in a && typeof (a as unknown as { gcd: (b: C) => C }).gcd === 'function') {
    return (a as unknown as { gcd: (b: C) => C }).gcd(b);
  }

  // For fields, GCD is always 1 (or a unit)
  // Check if the ring is a field by looking for inv method
  if ('inv' in a && typeof (a as unknown as { inv: () => C }).inv === 'function') {
    // In a field, gcd(a, b) = 1 for any nonzero a, b
    if (!a.isZero()) {
      // Return the multiplicative identity
      return a.mul((a as unknown as { inv: () => C }).inv()) as C; // This gives 1
    }
    if (!b.isZero()) {
      return b.mul((b as unknown as { inv: () => C }).inv()) as C;
    }
    return a; // Both zero, return zero
  }

  throw new ValueError('coefficient ring does not support gcd');
}

/**
 * Compute the determinant of a square matrix over an integral domain using
 * fraction-free (Bareiss) Gaussian elimination.
 *
 * Every division performed here is exact, so this is valid over any integral
 * domain -- in particular over ZZ, where the previous division-based
 * elimination silently truncated and returned wrong resultants/discriminants.
 *
 * @see Reference: sage/matrix/matrix2.pyx (determinant, "df" / Bareiss)
 */
function matrixDeterminant<C extends RingElement>(matrix: C[][], ring: CoefficientRing<C>): C {
  const n = matrix.length;
  if (n === 0) {
    return ring.one() as C;
  }

  // Make a copy of the matrix
  const M: C[][] = matrix.map((row) => [...row]);

  let sign = 1;
  let prevPivot = ring.one() as C;

  for (let col = 0; col < n - 1; col++) {
    // Find pivot
    let pivotRow = -1;
    for (let row = col; row < n; row++) {
      if (!M[row]![col]!.isZero()) {
        pivotRow = row;
        break;
      }
    }

    if (pivotRow === -1) {
      // Column is all zeros on and below the diagonal: determinant is 0
      return ring.zero() as C;
    }

    // Swap rows if needed
    if (pivotRow !== col) {
      [M[col], M[pivotRow]] = [M[pivotRow]!, M[col]!];
      sign = -sign;
    }

    const pivot = M[col]![col]!;

    for (let row = col + 1; row < n; row++) {
      for (let j = col + 1; j < n; j++) {
        // M[row][j] = (M[row][j]*pivot - M[row][col]*M[col][j]) / prevPivot
        const numer = M[row]![j]!.mul(pivot).sub(M[row]![col]!.mul(M[col]![j]!) as C) as C;
        M[row]![j] = prevPivot.eq(1) ? numer : divideCoeffs(numer, prevPivot);
      }
      M[row]![col] = ring.zero() as C;
    }

    prevPivot = pivot;
  }

  const det = M[n - 1]![n - 1]!;
  return sign === -1 ? (det.neg() as C) : det;
}

/**
 * Get characteristic from a ring (handles both property and method).
 */
function getRingCharacteristic<C extends RingElement>(ring: CoefficientRing<C>): bigint | null {
  if (!('characteristic' in ring)) return null;

  const char = (ring as { characteristic: unknown }).characteristic;

  // If it's a function, call it
  if (typeof char === 'function') {
    const result = (char as () => bigint | number).call(ring);
    return typeof result === 'number' ? BigInt(result) : result;
  }

  // If it's already a value
  if (typeof char === 'bigint') return char;
  if (typeof char === 'number') return BigInt(char);

  return null;
}

/**
 * Check if a ring is the integer ring ZZ.
 */
function isIntegerRing<C extends RingElement>(ring: CoefficientRing<C>): boolean {
  // Check for IntegerRing signature: is_field() returns false, characteristic() returns 0
  // and toString() returns 'Integer Ring'
  if (ring.toString && ring.toString() === 'Integer Ring') {
    return true;
  }
  // Check for duck typing: has is_field, is_integral_domain, characteristic
  if ('is_field' in ring && 'is_integral_domain' in ring && 'characteristic' in ring) {
    const r = ring as { is_field: () => boolean; is_integral_domain: () => boolean };
    const char = getRingCharacteristic(ring);
    if (!r.is_field() && r.is_integral_domain() && char === 0n) {
      // Check it's not QQ by verifying there's no fraction field marker
      if (!('is_absolute' in ring)) {
        return true;
      }
    }
  }
  return false;
}

/**
 * Check if a ring is the rational field QQ.
 */
function isRationalField<C extends RingElement>(ring: CoefficientRing<C>): boolean {
  if (ring instanceof FractionField_generic) return false;
  // Check for RationalField signature
  if (ring.toString && ring.toString() === 'Rational Field') {
    return true;
  }
  // Duck typing: is_field returns true, characteristic returns 0
  if ('is_field' in ring && 'characteristic' in ring) {
    const r = ring as { is_field: () => boolean };
    const char = getRingCharacteristic(ring);
    if (r.is_field() && char === 0n) {
      return true;
    }
  }
  return false;
}

/**
 * Check if a ring is a finite field.
 */
function isFiniteField<C extends RingElement>(ring: CoefficientRing<C>): boolean {
  // Check for is_field method
  if (ring.is_field && !ring.is_field()) {
    return false;
  }

  // Check for cardinality or order method
  if (
    !(
      'cardinality' in ring ||
      'order' in ring ||
      'characteristic' in ring ||
      Symbol.iterator in ring
    )
  ) {
    return false;
  }

  return true;
}

/**
 * Get the order (cardinality) of a finite field.
 */
function getFieldOrder<C extends RingElement>(ring: CoefficientRing<C>): bigint {
  if (
    'cardinality' in ring &&
    typeof (ring as { cardinality: () => bigint }).cardinality === 'function'
  ) {
    return (ring as { cardinality: () => bigint }).cardinality();
  }

  if ('order' in ring) {
    const order = (ring as { order: bigint | number }).order;
    return typeof order === 'number' ? BigInt(order) : order;
  }

  if ('characteristic' in ring && 'degree' in ring) {
    const p = (ring as { characteristic: bigint }).characteristic;
    const n = (ring as { degree: number }).degree;
    return p ** BigInt(n);
  }

  throw new ValueError('cannot determine field order');
}

/**
 * Get the characteristic of the base ring.
 */
function getCharacteristic<C extends RingElement>(ring: CoefficientRing<C>): bigint {
  if ('characteristic' in ring) {
    const p = (ring as { characteristic: bigint | number }).characteristic;
    return typeof p === 'number' ? BigInt(p) : p;
  }

  // Try to determine characteristic by computing 1 + 1 + ... until we get 0
  let sum = ring.one();
  const one = ring.one();

  for (let i = 1; i < 1000; i++) {
    sum = sum.add(one) as C;
    if (sum.isZero()) {
      return BigInt(i + 1);
    }
  }

  return 0n; // Assume characteristic 0 if not found
}

/**
 * Compute base^exp mod modulus for polynomials.
 */
function powerMod<C extends RingElement>(
  base: Polynomial<C>,
  exp: bigint,
  modulus: Polynomial<C>
): Polynomial<C> {
  if (exp === 0n) {
    return base.parent.one();
  }

  let result = base.parent.one();
  let b = base.mod(modulus);

  while (exp > 0n) {
    if ((exp & 1n) === 1n) {
      result = result.mul(b).mod(modulus);
    }
    b = b.mul(b).mod(modulus);
    exp >>= 1n;
  }

  return result;
}

/**
 * Compute base^(q^k) mod modulus by iterating k q-th powers.
 *
 * This is FLINT's `nmod_poly_powpowmod` (`nmod_poly_factor/is_irreducible.c`).
 */
function powerModIterated<C extends RingElement>(
  base: Polynomial<C>,
  q: bigint,
  k: number,
  modulus: Polynomial<C>
): Polynomial<C> {
  let result = base.mod(modulus);
  for (let i = 0; i < k; i++) {
    result = powerMod(result, q, modulus);
  }
  return result;
}

/**
 * Cantor-Zassenhaus algorithm for equal-degree factorization.
 *
 * Given a polynomial f that is a product of distinct irreducible polynomials
 * all of the same degree d, find all irreducible factors.
 */
function cantorZassenhausFactorization<C extends RingElement>(
  f: Polynomial<C>,
  degree: number,
  ring: PolynomialRingBase<C>
): Polynomial<C>[] {
  const n = f.degree();

  if (n === 0) {
    return [];
  }

  if (n === degree) {
    // f is already irreducible
    return [f._monic()];
  }

  const baseRing = ring.base_ring;
  const q = getFieldOrder(baseRing);
  const p = getCharacteristic(baseRing);

  // We expect to succeed with probability > 1/2 per attempt, so 100 failures
  // means there is a bug (`polynomial_element.pyx:2205`).
  const maxAttempts = 100;

  for (let attempt = 0; attempt < maxAttempts; attempt++) {
    // Sample T uniformly from R of degree exactly 2*degree + 1, then make it
    // monic (`polynomial_element.pyx:2209`).
    const t = randomPolynomial(ring, 2 * degree + 1).monic();

    let h: Polynomial<C>;

    if (p === 2n) {
      // Characteristic 2: use trace
      // Compute T + T^2 + T^4 + ... + T^{2^{dk-1}} mod f
      // where k is the degree of the base field over GF(p)
      const fieldDegree = getFieldDegree(baseRing);
      const numTerms = degree * fieldDegree;

      let c = t.mod(f);
      let tt = t.mod(f);

      for (let i = 1; i < numTerms; i++) {
        tt = tt.mul(tt).mod(f); // T^{2^i}
        c = c.add(tt);
      }

      h = f.gcd(c);
    } else {
      // Odd characteristic: use (q^d - 1)/2 power
      const exponent = (q ** BigInt(degree) - 1n) / 2n;
      const tPow = powerMod(t, exponent, f);
      const tPowMinus1 = tPow.sub(ring.one());
      h = f.gcd(tPowMinus1);
    }

    const hd = h.degree();

    // Check if we found a non-trivial factor
    if (hd > 0 && hd < n) {
      // Recursively factor both parts
      const factors1 = cantorZassenhausFactorization(h._monic(), degree, ring);
      const quotient = f.quo_rem(h)[0];
      const factors2 = cantorZassenhausFactorization(quotient._monic(), degree, ring);

      return [...factors1, ...factors2];
    }
  }

  // Sage raises an AssertionError here rather than returning an unsplit
  // factor (`polynomial_element.pyx:2236`): reaching this point means the
  // input was not a product of distinct irreducibles of the given degree,
  // or that the sampler is broken.
  throw new AssertionError(`no splitting of degree ${degree} found for ${f}`);
}

/**
 * Return a uniformly random element of a (finite) coefficient ring.
 *
 * Sampling via `baseRing.__call__(someNumber)` is wrong for extension
 * fields: the number is routed through the prime subfield, so the sampled
 * element never leaves GF(p) and Cantor-Zassenhaus can never split a
 * polynomial over GF(p^k).
 */
function randomRingElement<C extends RingElement>(ring: CoefficientRing<C>): C {
  if (
    'random_element' in ring &&
    typeof (ring as { random_element: () => C }).random_element === 'function'
  ) {
    return (ring as { random_element: () => C }).random_element();
  }

  if (Symbol.iterator in ring) {
    const elements = [...(ring as unknown as Iterable<C>)];
    if (elements.length === 0) {
      throw new ValueError('cannot sample from an empty ring');
    }
    const index = Number(current_randstate().random_below(BigInt(elements.length)));
    return elements[index]!;
  }

  throw new NotImplementedError(`cannot sample a random element of ${ring}`);
}

/**
 * Generate a random polynomial of degree exactly `degree`, sampling every
 * coefficient uniformly from the base ring.
 *
 * @see Reference: sage/rings/polynomial/polynomial_ring.py:1344 (random_element)
 */
function randomPolynomial<C extends RingElement>(
  ring: PolynomialRingBase<C>,
  degree: number
): Polynomial<C> {
  const baseRing = ring.base_ring;

  if (degree < 0) {
    return ring.zero();
  }

  const coeffs: C[] = [];
  for (let i = 0; i < degree; i++) {
    coeffs.push(randomRingElement(baseRing));
  }

  // The leading coefficient must be nonzero so that the degree is exactly
  // `degree` (Sage's `random_element(d)` samples until this holds).
  let lead = randomRingElement(baseRing);
  while (lead.isZero()) {
    lead = randomRingElement(baseRing);
  }
  coeffs.push(lead);

  return new Polynomial(coeffs, ring);
}

/**
 * Get the degree of a finite field over its prime field.
 */
function getFieldDegree<C extends RingElement>(ring: CoefficientRing<C>): number {
  if ('degree' in ring) {
    return (ring as { degree: number }).degree;
  }

  // Compute from characteristic and order
  const q = getFieldOrder(ring);
  const p = getCharacteristic(ring);

  if (p === 0n) {
    return 1;
  }

  let degree = 0;
  let power = 1n;

  while (power < q) {
    power *= p;
    degree++;
  }

  return degree;
}

// ============================================
// Integer/Rational Polynomial Factorization
// ============================================

/**
 * Extract bigint coefficients from a polynomial over ZZ.
 */
function extractIntegerCoeffs<C extends RingElement>(poly: Polynomial<C>): bigint[] {
  return poly.coeffs.map((c) => {
    // Handle Integer wrapper class
    if ('value' in c && typeof (c as { value: bigint }).value === 'bigint') {
      return (c as { value: bigint }).value;
    }
    // Handle raw bigint (shouldn't happen but for safety)
    if (typeof c === 'bigint') {
      return c;
    }
    // Try toString and parse
    return BigInt(c.toString());
  });
}

/**
 * Compute the content of an integer polynomial (GCD of coefficients).
 */
function intPolyContent(coeffs: bigint[]): bigint {
  if (coeffs.length === 0) return 0n;
  let g = coeffs[0]!;
  for (let i = 1; i < coeffs.length; i++) {
    g = gcdBigInt(g, coeffs[i]!);
    if (g === 1n || g === -1n) return 1n;
  }
  return g < 0n ? -g : g;
}

/**
 * Divide all coefficients by a constant.
 */
function intPolyDivideByConstant(coeffs: bigint[], c: bigint): bigint[] {
  return coeffs.map((coeff) => coeff / c);
}

/**
 * Make polynomial primitive (divide by content).
 */
function intPolyPrimitive(coeffs: bigint[]): [bigint, bigint[]] {
  if (coeffs.length === 0) return [1n, []];
  const content = intPolyContent(coeffs);
  if (content === 0n) return [1n, coeffs];
  // Make leading coefficient positive
  const lc = coeffs[coeffs.length - 1]!;
  const sign = lc < 0n ? -1n : 1n;
  const adjustedContent = content * sign;
  return [adjustedContent, intPolyDivideByConstant(coeffs, adjustedContent)];
}

/**
 * Multiply two integer polynomials.
 */
function intPolyMul(a: bigint[], b: bigint[]): bigint[] {
  if (a.length === 0 || b.length === 0) return [];
  const result = new Array(a.length + b.length - 1).fill(0n);
  for (let i = 0; i < a.length; i++) {
    for (let j = 0; j < b.length; j++) {
      result[i + j] += a[i]! * b[j]!;
    }
  }
  return result;
}

/**
 * Compute quotient and remainder of integer polynomial division.
 * Returns [quotient, remainder] or null if division is not exact in ZZ.
 */
function intPolyQuoRem(a: bigint[], b: bigint[]): [bigint[], bigint[]] | null {
  if (b.length === 0) throw new ZeroDivisionError('polynomial division by zero');

  // Remove trailing zeros
  while (a.length > 0 && a[a.length - 1] === 0n) a = a.slice(0, -1);
  while (b.length > 0 && b[b.length - 1] === 0n) b = b.slice(0, -1);

  if (a.length < b.length) return [[0n], a];

  const degA = a.length - 1;
  const degB = b.length - 1;
  const lcB = b[degB]!;

  const quotient = new Array(degA - degB + 1).fill(0n);
  const remainder = [...a];

  for (let i = degA; i >= degB; i--) {
    if (remainder[i] === 0n) continue;

    // Check if division is exact
    if (remainder[i]! % lcB !== 0n) {
      return null; // Not exact division in ZZ
    }

    const qCoeff = remainder[i]! / lcB;
    quotient[i - degB] = qCoeff;

    for (let j = 0; j <= degB; j++) {
      remainder[i - degB + j] -= qCoeff * b[j]!;
    }
  }

  // Remove trailing zeros from remainder
  while (remainder.length > 0 && remainder[remainder.length - 1] === 0n) {
    remainder.pop();
  }

  return [quotient, remainder];
}

/**
 * Evaluate integer polynomial at a point.
 */
function intPolyEval(coeffs: bigint[], x: bigint): bigint {
  if (coeffs.length === 0) return 0n;
  let result = coeffs[coeffs.length - 1]!;
  for (let i = coeffs.length - 2; i >= 0; i--) {
    result = result * x + coeffs[i]!;
  }
  return result;
}

/**
 * Compute polynomial modulo a prime (reduce coefficients mod p).
 */
function intPolyModP(coeffs: bigint[], p: bigint): bigint[] {
  // Ensure p is bigint
  const pBig = typeof p === 'bigint' ? p : BigInt(p);
  const result = coeffs.map((c) => {
    // Ensure c is bigint
    const cBig = typeof c === 'bigint' ? c : BigInt(c);
    let r = cBig % pBig;
    if (r < 0n) r += pBig;
    return r;
  });
  // Remove trailing zeros
  while (result.length > 0 && result[result.length - 1] === 0n) {
    result.pop();
  }
  return result;
}

/**
 * Modular polynomial multiplication.
 */
function modPolyMul(a: bigint[], b: bigint[], p: bigint): bigint[] {
  if (a.length === 0 || b.length === 0) return [];
  const pBig = toBigIntSafe(p);
  const result: bigint[] = new Array(a.length + b.length - 1).fill(0n);
  for (let i = 0; i < a.length; i++) {
    for (let j = 0; j < b.length; j++) {
      const ai = toBigIntSafe(a[i]);
      const bj = toBigIntSafe(b[j]);
      result[i + j] = (result[i + j]! + ai * bj) % pBig;
    }
  }
  // Remove trailing zeros
  while (result.length > 0 && result[result.length - 1] === 0n) {
    result.pop();
  }
  return result;
}

/**
 * Modular polynomial quotient and remainder.
 */
function modPolyQuoRem(a: bigint[], b: bigint[], p: bigint): [bigint[], bigint[]] {
  if (b.length === 0) throw new ZeroDivisionError('polynomial division by zero');

  const pBig = toBigIntSafe(p);

  // Ensure all values are bigints and remove trailing zeros
  let aCopy = a.map((c) => toBigIntSafe(c));
  let bCopy = b.map((c) => toBigIntSafe(c));

  while (aCopy.length > 0 && aCopy[aCopy.length - 1] === 0n) aCopy = aCopy.slice(0, -1);
  while (bCopy.length > 0 && bCopy[bCopy.length - 1] === 0n) bCopy = bCopy.slice(0, -1);

  if (aCopy.length === 0 || aCopy.length < bCopy.length) return [[0n], aCopy];

  const degA = aCopy.length - 1;
  const degB = bCopy.length - 1;
  const lcB = bCopy[degB]!;
  const lcBInv = modInverse(lcB, pBig);

  const quotient: bigint[] = new Array(degA - degB + 1).fill(0n);
  const remainder = aCopy.map((c) => ((c % pBig) + pBig) % pBig);

  for (let i = degA; i >= degB; i--) {
    if (remainder[i] === 0n) continue;

    const qCoeff = (((remainder[i]! * lcBInv) % pBig) + pBig) % pBig;
    quotient[i - degB] = qCoeff;

    for (let j = 0; j <= degB; j++) {
      remainder[i - degB + j] =
        (((remainder[i - degB + j]! - qCoeff * bCopy[j]!) % pBig) + pBig) % pBig;
    }
  }

  // Remove trailing zeros
  while (remainder.length > 0 && remainder[remainder.length - 1] === 0n) {
    remainder.pop();
  }
  while (quotient.length > 0 && quotient[quotient.length - 1] === 0n) {
    quotient.pop();
  }

  return [quotient.length > 0 ? quotient : [0n], remainder];
}

/**
 * Convert a value to bigint safely.
 */
function toBigIntSafe(x: unknown): bigint {
  if (typeof x === 'bigint') return x;
  if (typeof x === 'number') return BigInt(Math.floor(x));
  if (typeof x === 'string') return BigInt(x);
  if (x && typeof x === 'object' && 'value' in x) {
    const v = (x as { value: unknown }).value;
    if (typeof v === 'bigint') return v;
    if (typeof v === 'number') return BigInt(Math.floor(v));
  }
  throw new Error(`Cannot convert ${typeof x} to bigint: ${x}`);
}

/**
 * Modular inverse using extended GCD.
 */
function modInverse(a: bigint, p: bigint): bigint {
  // Ensure inputs are bigints
  const aBig = toBigIntSafe(a);
  const pBig = toBigIntSafe(p);

  let [oldR, r] = [((aBig % pBig) + pBig) % pBig, pBig];
  let [oldS, s] = [1n, 0n];

  while (r !== 0n) {
    const q = oldR / r;
    [oldR, r] = [r, oldR - q * r];
    [oldS, s] = [s, oldS - q * s];
  }

  return ((oldS % pBig) + pBig) % pBig;
}

/**
 * Modular GCD of two polynomials.
 */
function modPolyGcd(a: bigint[], b: bigint[], p: bigint): bigint[] {
  const pBig = toBigIntSafe(p);
  while (b.length > 0) {
    const [_, rem] = modPolyQuoRem(a, b, pBig);
    a = b;
    b = rem;
  }
  // Make monic
  if (a.length > 0 && a[a.length - 1] !== 0n) {
    const lcInv = modInverse(a[a.length - 1]!, pBig);
    a = a.map((c) => (((toBigIntSafe(c) * lcInv) % pBig) + pBig) % pBig);
  }
  return a;
}

/**
 * Integer square root rounded down: the largest `s >= 0` with `s*s <= n`.
 *
 * This is FLINT's `fmpz_sqrt`, computed exactly (a `Math.sqrt` of a bigint both
 * overflows and rounds the wrong way).
 */
function isqrtFloor(n: bigint): bigint {
  if (n <= 0n) return 0n;
  if (n < 4n) return 1n;
  // Newton's iteration for floor(sqrt(n)), started above the root.
  let x = 1n << BigInt((n.toString(2).length >> 1) + 1);
  let y = (x + n / x) >> 1n;
  while (y < x) {
    x = y;
    y = (x + n / x) >> 1n;
  }
  return x;
}

/**
 * Mignotte's bound on the coefficients of a factor of an integer polynomial.
 *
 * Verbatim transcription of FLINT's `_fmpz_poly_factor_mignotte`
 * (`reference/flint/src/fmpz_poly_factor/factor_zassenhaus.c:47-88`), including
 * its initialisation `b = m - 1` (the comment there describes
 * `b = binomial(m-1, j-1)`, so the value FLINT actually computes is `m - 1`
 * times that; the result is therefore a valid -- merely more generous -- upper
 * bound, and we reproduce it exactly rather than "fixing" it).
 *
 * For `g` dividing `f` of degree `m`, every coefficient `b_j` of `g` satisfies
 * `|b_j| <= B`.
 */
function fmpz_poly_factor_mignotte(f: bigint[]): bigint {
  const m = f.length - 1;
  if (m < 0) return 0n;
  const abs = (x: bigint) => (x < 0n ? -x : x);

  let f2 = 0n;
  for (let j = 0; j <= m; j++) f2 += f[j]! * f[j]!;
  f2 = isqrtFloor(f2) + 1n;

  const lc = abs(f[m]!);
  let B = abs(f[0]!);

  let b = BigInt(m - 1);
  for (let j = 1; j < m; j++) {
    const t = b * lc;
    b = (b * BigInt(m - j)) / BigInt(j);
    const s = b * f2 + t;
    if (B < s) B = s;
  }

  if (B < lc) B = lc;
  return B;
}

/**
 * Reduce a polynomial's coefficients into `[0, p)` and strip trailing zeros.
 */
function modPolyNormalize(a: bigint[], p: bigint): bigint[] {
  const result = a.map((c) => ((c % p) + p) % p);
  while (result.length > 0 && result[result.length - 1] === 0n) result.pop();
  return result;
}

/**
 * Add two polynomials over Z/pZ.
 */
function modPolyAdd(a: bigint[], b: bigint[], p: bigint): bigint[] {
  const len = Math.max(a.length, b.length);
  const result: bigint[] = new Array(len).fill(0n);
  for (let i = 0; i < len; i++) {
    result[i] = ((((a[i] ?? 0n) + (b[i] ?? 0n)) % p) + p) % p;
  }
  while (result.length > 0 && result[result.length - 1] === 0n) result.pop();
  return result;
}

/**
 * Subtract two polynomials over Z/pZ.
 */
function modPolySub(a: bigint[], b: bigint[], p: bigint): bigint[] {
  const len = Math.max(a.length, b.length);
  const result: bigint[] = new Array(len).fill(0n);
  for (let i = 0; i < len; i++) {
    result[i] = ((((a[i] ?? 0n) - (b[i] ?? 0n)) % p) + p) % p;
  }
  while (result.length > 0 && result[result.length - 1] === 0n) result.pop();
  return result;
}

/**
 * Make a nonzero polynomial monic over Z/pZ.
 */
function modPolyMonic(a: bigint[], p: bigint): bigint[] {
  const f = modPolyNormalize(a, p);
  if (f.length === 0) return f;
  const lc = f[f.length - 1]!;
  if (lc === 1n) return f;
  const inv = modInverse(lc, p);
  return f.map((c) => (c * inv) % p);
}

/**
 * `base^e mod m` over Z/pZ.
 */
function modPolyPowMod(base: bigint[], e: bigint, m: bigint[], p: bigint): bigint[] {
  let result: bigint[] = [1n];
  let b = modPolyQuoRem(modPolyNormalize(base, p), m, p)[1];
  let k = e;
  while (k > 0n) {
    if ((k & 1n) === 1n) {
      result = modPolyQuoRem(modPolyMul(result, b, p), m, p)[1];
    }
    k >>= 1n;
    if (k > 0n) {
      b = modPolyQuoRem(modPolyMul(b, b, p), m, p)[1];
    }
  }
  return result;
}

/**
 * Inverse of `a` modulo `m` over Z/pZ, via the extended Euclidean algorithm.
 *
 * @throws {ArithmeticError} if `a` is not invertible modulo `m`
 */
function modPolyInvMod(a: bigint[], m: bigint[], p: bigint): bigint[] {
  let r0 = modPolyNormalize(m, p);
  let r1 = modPolyQuoRem(modPolyNormalize(a, p), m, p)[1];
  let s0: bigint[] = [];
  let s1: bigint[] = [1n];

  while (r1.length > 0) {
    const [q, rem] = modPolyQuoRem(r0, r1, p);
    const s = modPolySub(s0, modPolyMul(modPolyNormalize(q, p), s1, p), p);
    r0 = r1;
    r1 = rem;
    s0 = s1;
    s1 = s;
  }

  if (r0.length !== 1) {
    throw new ArithmeticError('polynomial is not invertible modulo the given modulus');
  }
  const inv = modInverse(r0[0]!, p);
  return modPolyNormalize(
    s0.map((c) => c * inv),
    p
  );
}

/**
 * Distinct-degree factorization of a monic squarefree polynomial over Z/pZ.
 *
 * Returns pairs `[g_d, d]` where `g_d` is the product of all monic irreducible
 * factors of degree `d`.
 *
 * @see Reference: sage/rings/polynomial/polynomial_element.pyx:_distinct_degree_factorisation_squarefree
 */
function modpDistinctDegree(f: bigint[], p: bigint): Array<[bigint[], number]> {
  const result: Array<[bigint[], number]> = [];
  let v = modPolyMonic(f, p);
  const x: bigint[] = [0n, 1n];
  // w = x^(p^d) mod v, maintained across iterations.
  let w = modPolyQuoRem(x, v, p)[1];
  let d = 0;

  while (2 * (d + 1) <= v.length - 1) {
    d += 1;
    w = modPolyPowMod(w, p, v, p);
    const g = modPolyGcd(v, modPolySub(w, x, p), p);
    if (g.length > 1) {
      result.push([g, d]);
      v = modPolyQuoRem(v, g, p)[0];
      w = modPolyQuoRem(w, v, p)[1];
    }
  }

  if (v.length > 1) {
    result.push([v, v.length - 1]);
  }

  return result;
}

/**
 * Equal-degree (Cantor-Zassenhaus) splitting over Z/pZ.
 *
 * `f` is monic, squarefree and a product of monic irreducible factors that all
 * have degree `d`.  For odd `p` the splitting element is `a^((p^d-1)/2) - 1`;
 * for `p = 2` it is the trace `a + a^2 + ... + a^(2^(d-1))` (von zur Gathen &
 * Gerhard, Algorithms 14.8/14.10 -- the same case distinction FLINT makes in
 * `nmod_poly_factor_equal_deg`).
 */
function modpEqualDegree(f: bigint[], d: number, p: bigint): bigint[][] {
  const monic = modPolyMonic(f, p);
  const n = monic.length - 1;
  if (n === d) return [monic];

  const rstate = current_randstate();
  // One trial splits with probability >= 1/2, so 512 consecutive failures do
  // not happen; we raise rather than silently drop a factor.
  const maxAttempts = 512;

  for (let attempt = 0; attempt < maxAttempts; attempt++) {
    const a: bigint[] = [];
    for (let i = 0; i < n; i++) {
      a.push(rstate.random_below(p));
    }
    const aNorm = modPolyNormalize(a, p);
    if (aNorm.length <= 1) continue;

    // A common factor with a random polynomial already splits f.
    let g = modPolyGcd(monic, aNorm, p);
    if (g.length <= 1 || g.length >= monic.length) {
      let b: bigint[];
      if (p === 2n) {
        // Trace map: a + a^2 + a^4 + ... + a^(2^(d-1)) mod f
        b = [];
        let term = modPolyQuoRem(aNorm, monic, p)[1];
        for (let i = 0; i < d; i++) {
          b = modPolyAdd(b, term, p);
          term = modPolyQuoRem(modPolyMul(term, term, p), monic, p)[1];
        }
      } else {
        const e = (p ** BigInt(d) - 1n) / 2n;
        b = modPolySub(modPolyPowMod(aNorm, e, monic, p), [1n], p);
      }
      g = modPolyGcd(monic, b, p);
    }

    if (g.length > 1 && g.length < monic.length) {
      const h = modPolyQuoRem(monic, g, p)[0];
      return [...modpEqualDegree(g, d, p), ...modpEqualDegree(h, d, p)];
    }
  }

  throw new ArithmeticError(
    `Cantor-Zassenhaus failed to split a degree-${n} polynomial modulo ${p}`
  );
}

/**
 * Factor a squarefree polynomial over Z/pZ into monic irreducible factors.
 */
function modpFactorSquarefree(f: bigint[], p: bigint): bigint[][] {
  const monic = modPolyMonic(f, p);
  if (monic.length <= 1) return [];
  if (monic.length === 2) return [monic];

  const factors: bigint[][] = [];
  for (const [g, d] of modpDistinctDegree(monic, p)) {
    if (g.length - 1 === d) {
      factors.push(g);
    } else {
      factors.push(...modpEqualDegree(g, d, p));
    }
  }
  return factors;
}

/**
 * Multifactor Hensel lifting.
 *
 * Given `f = lc(f) * h_1 * ... * h_r (mod p)` with the `h_i` monic, pairwise
 * coprime and `p` not dividing `lc(f)`, return monic `H_i = h_i (mod p)` with
 * `f = lc(f) * H_1 * ... * H_r (mod p^k)`.
 *
 * Linear lifting: writing `H_i' = H_i + p^j d_i`, the corrections satisfy
 * `sum_i d_i * prod_{l != i} h_l = (f - lc(f) * prod_i H_i) / (p^j * lc(f))`
 * modulo `p`, which is solved by `d_i = (c * s_i) mod h_i` where `s_i` is the
 * inverse of `prod_{l != i} h_l` modulo `h_i` and `c` is the right-hand side.
 *
 * @see Reference: von zur Gathen & Gerhard, "Modern Computer Algebra",
 *   Algorithm 15.17 (multifactor Hensel lifting); FLINT
 *   `fmpz_poly_factor/hensel_lift.c`.
 */
function henselLiftFactors(f: bigint[], hs: bigint[][], p: bigint, k: number): bigint[][] {
  const r = hs.length;
  const b = f[f.length - 1]!;
  const bInv = modInverse(((b % p) + p) % p, p);

  // s_i = (prod_{l != i} h_l)^-1 mod h_i, over Z/pZ
  const s: bigint[][] = [];
  for (let i = 0; i < r; i++) {
    let u: bigint[] = [1n];
    for (let l = 0; l < r; l++) {
      if (l !== i) u = modPolyQuoRem(modPolyMul(u, hs[l]!, p), hs[i]!, p)[1];
    }
    s.push(modPolyInvMod(u, hs[i]!, p));
  }

  const H = hs.map((h) => modPolyNormalize(h, p));
  let pj = p;

  for (let j = 1; j < k; j++) {
    const pj1 = pj * p;

    // e = (f - b * prod H_i) / p^j  (mod p)
    let prod: bigint[] = [b];
    for (const h of H) prod = intPolyMul(prod, h);
    const diff: bigint[] = [];
    for (let i = 0; i < Math.max(f.length, prod.length); i++) {
      let c = ((f[i] ?? 0n) - (prod[i] ?? 0n)) % pj1;
      if (c < 0n) c += pj1;
      diff.push(c);
    }

    if (diff.some((c) => c !== 0n)) {
      const e = modPolyNormalize(
        diff.map((c) => {
          if (c % pj !== 0n) {
            throw new ArithmeticError('Hensel lifting lost its invariant');
          }
          return c / pj;
        }),
        p
      );
      const c = modPolyNormalize(
        e.map((coeff) => coeff * bInv),
        p
      );

      for (let i = 0; i < r; i++) {
        const di = modPolyQuoRem(modPolyMul(c, s[i]!, p), hs[i]!, p)[1];
        const lifted = H[i]!.slice();
        for (let idx = 0; idx < di.length; idx++) {
          lifted[idx] = ((lifted[idx] ?? 0n) + pj * di[idx]!) % pj1;
        }
        H[i] = lifted;
      }
    }

    pj = pj1;
  }

  return H;
}

/* ===================================================================== *
 * Factoring in ZZ[x]: Zassenhaus recombination and van Hoeij / LLL
 *
 * Sage's `ZZ[x].factor()` delegates to FLINT's `fmpz_poly_factor`, so the code
 * from here to `factorSquarefreeIntPoly` is a transcription of FLINT, routine
 * by routine:
 *
 *   reference/flint/src/fmpz_poly_factor/factor_zassenhaus.c
 *   reference/flint/src/fmpz_poly_factor/factor_zassenhaus_recombination.c
 *   reference/flint/src/fmpz_poly_factor/zassenhaus_subset.c
 *   reference/flint/src/fmpz_poly_factor/zassenhaus_prune.c
 *   reference/flint/src/fmpz_poly_factor/factor_van_hoeij.c
 *   reference/flint/src/fmpz_poly_factor/CLD_mat.c
 *   reference/flint/src/fmpz_poly_factor/van_hoeij_check_if_solved.c
 *   reference/flint/src/fmpz_poly/CLD_bound.c
 *   reference/flint/src/fmpz_poly/divlow_smodp.c
 *   reference/flint/src/fmpz_poly/divhigh_smodp.c
 *   reference/flint/src/fmpz_mat/next_col_van_hoeij.c
 *   reference/flint/src/fmpz_mat/col_partition.c
 *   reference/flint/src/fmpz_lll/  (LLL with removal, `knapsack` variant)
 * ===================================================================== */

/** Number of bits of `|x|` (FLINT `fmpz_bits`; `0` for `x = 0`). */
function fmpzBits(x: bigint): number {
  const a = x < 0n ? -x : x;
  return a === 0n ? 0 : a.toString(2).length;
}

/** `FLINT_BIT_COUNT`: one more than the index of the highest set bit. */
function bitCountU(x: number): number {
  let n = 0;
  let v = Math.floor(x);
  while (v > 0) {
    n++;
    v = Math.floor(v / 2);
  }
  return n;
}

/** `|x|` for bigints. */
function absBig(x: bigint): bigint {
  return x < 0n ? -x : x;
}

/**
 * Symmetric remainder: the `r` with `r == a (mod m)` and `-m/2 < r <= m/2`
 * (FLINT `fmpz_smod`, `reference/flint/src/fmpz/smod.c`).
 */
function fmpzSmod(a: bigint, m: bigint): bigint {
  let r = a % m;
  if (r < 0n) r += m;
  if (r > m / 2n) r -= m;
  return r;
}

/** `fmpz_tdiv_q_2exp`: division by `2^k`, truncating towards zero. */
function tdivQ2exp(a: bigint, k: number): bigint {
  if (k <= 0) return a << BigInt(-k);
  return a < 0n ? -(-a >> BigInt(k)) : a >> BigInt(k);
}

/** Strip trailing zero coefficients (`_fmpz_poly_normalise`). */
function zpNormalise(a: bigint[]): bigint[] {
  let n = a.length;
  while (n > 0 && a[n - 1] === 0n) n--;
  return n === a.length ? a : a.slice(0, n);
}

/** Coefficientwise symmetric remainder (`fmpz_poly_scalar_smod_fmpz`). */
function zpSmod(a: bigint[], m: bigint): bigint[] {
  return zpNormalise(a.map((c) => fmpzSmod(c, m)));
}

/** `gcd(a, m)` together with `a^-1 mod m` when the gcd is 1 (`fmpz_gcdinv`). */
function gcdInvMod(a: bigint, m: bigint): [bigint, bigint] {
  let [oldR, r] = [((a % m) + m) % m, m];
  let [oldS, s] = [1n, 0n];
  while (r !== 0n) {
    const q = oldR / r;
    [oldR, r] = [r, oldR - q * r];
    [oldS, s] = [s, oldS - q * s];
  }
  return [oldR, ((oldS % m) + m) % m];
}

/* ------------------------------------------------------------------ *
 * zassenhaus_prune  (fmpz_poly_factor/zassenhaus_prune.c and
 *                    fmpz_poly_factor.h:97-145)
 * ------------------------------------------------------------------ */

interface ZassenhausPrune {
  deg: number;
  /** bit 0 = possible degree so far, bit 1 = possible for the current prime */
  posDegs: Uint8Array;
  newLength: number;
  newTotal: number;
  newDegs: number[];
}

/** `zassenhaus_prune_set_degree` (zassenhaus_prune.c:22-52). */
function zassenhausPruneSetDegree(d: number): ZassenhausPrune {
  if (d < 1) throw new ValueError('zassenhaus_prune_set_degree');
  const posDegs = new Uint8Array(d + 1);
  posDegs.fill(1);
  return { deg: d, posDegs, newLength: 0, newTotal: 0, newDegs: [] };
}

/** `zassenhaus_prune_start_add_factors` (fmpz_poly_factor.h:123-128). */
function zassenhausPruneStartAddFactors(Z: ZassenhausPrune): void {
  Z.newLength = 0;
  Z.newTotal = 0;
}

/** `zassenhaus_prune_add_factor` (zassenhaus_prune.c:54-72). */
function zassenhausPruneAddFactor(Z: ZassenhausPrune, deg: number, exp: number): void {
  if (exp < 1 || deg < 1) return;
  for (let i = 0; i < exp; i++) {
    if (Z.newLength >= Z.deg) throw new ArithmeticError('zassenhaus_prune_add_factor');
    Z.newTotal += deg;
    Z.newDegs[Z.newLength] = deg;
    Z.newLength++;
  }
}

/** `zassenhaus_prune_end_add_factors` (zassenhaus_prune.c:74-114). */
function zassenhausPruneEndAddFactors(Z: ZassenhausPrune): void {
  const a = Z.posDegs;
  const posMask = 1;
  const newMask = 2;

  if (Z.newTotal !== Z.deg) throw new ArithmeticError('zassenhaus_prune_add_factor');

  a[0] = a[0]! | newMask;
  for (let j = 1; j <= Z.deg; j++) a[j] = a[j]! & ~newMask;

  for (let i = 0; i < Z.newLength; i++) {
    const d = Z.newDegs[i]!;
    for (let j = Z.deg; j >= 0; j--) {
      if ((a[j]! & newMask) !== 0) {
        if (j + d > Z.deg) throw new ArithmeticError('zassenhaus_prune_add_factor');
        a[j + d] = a[j + d]! | newMask;
      }
    }
  }

  // merge new possibilities with old
  for (let j = 0; j <= Z.deg; j++) a[j] = a[j]! & (a[j]! >> 1);

  // 0 and deg should always be possible
  if (a[0] !== posMask || a[Z.deg] !== posMask) {
    throw new ArithmeticError('zassenhaus_prune_add_factor');
  }
}

/** `zassenhaus_prune_degree_is_possible` (fmpz_poly_factor.h:137-146). */
function zassenhausPruneDegreeIsPossible(Z: ZassenhausPrune, d: number): boolean {
  if (d <= 0) return d === 0;
  if (d >= Z.deg) return d === Z.deg;
  return Z.posDegs[d] !== 0;
}

/* ------------------------------------------------------------------ *
 * zassenhaus_subset  (fmpz_poly_factor/zassenhaus_subset.c)
 *
 * A subset of {0,...,r-1} is encoded in `s[0..r-1]`: `s[i] >= 0` means the
 * index `s[i]` is *in* the subset, `s[i] < 0` means the index `-s[i]-1` is out.
 * ------------------------------------------------------------------ */

/** `zassenhaus_subset_first` (zassenhaus_subset.c:15-27). */
function zassenhausSubsetFirst(s: number[], r: number, m: number): void {
  for (let i = 0; i < r; i++) {
    if (i >= m) s[i] = s[i]! < 0 ? s[i]! : -s[i]! - 1;
    else s[i] = s[i]! >= 0 ? s[i]! : -s[i]! - 1;
  }
}

/** `zassenhaus_subset_next` (zassenhaus_subset.c:29-59). */
function zassenhausSubsetNext(s: number[], r: number): boolean {
  let i = 0;
  while (i < r && s[i]! < 0) i++;
  const j = i;
  while (i < r && s[i]! >= 0) i++;
  const k = i;

  if (k === 0 || k >= r) return false;

  s[k] = -s[k]! - 1;
  s[k - 1] = -s[k - 1]! - 1;

  if (j > 0) {
    for (let t = 0; t < k - j - 1; t++) if (s[t]! < 0) s[t] = -s[t]! - 1;
    for (let t = k - j - 1; t < k - 1; t++) if (s[t]! >= 0) s[t] = -s[t]! - 1;
  }
  return true;
}

/** `zassenhaus_subset_next_disjoint` (zassenhaus_subset.c:61-95). */
function zassenhausSubsetNextDisjoint(s: number[], r: number): boolean {
  let total = 0;
  let last = r - 1;
  for (let i = 0; i < r; i++) {
    if (s[i]! >= 0) {
      total++;
      last = i;
    }
  }

  let j = 0;
  for (let i = 0; i < r; i++) if (s[i]! < 0) s[j++] = s[i]!;

  if (r - total < total || total < 1 || last === r - 1) return false;

  const min = Math.min(total - 1, last - total + 1);

  for (let i = 0; i < min; i++) s[i] = -s[i]! - 1;
  for (let i = last - total + 1; i < last - min + 1; i++) s[i] = -s[i]! - 1;

  return true;
}

/**
 * `_fmpz_poly_product` (factor_zassenhaus_recombination.c:18-77): the product
 * of the selected lifted factors, scaled by `leadf` and reduced into the
 * symmetric range modulo `P`.
 */
function zassenhausSubsetProduct(
  liftedFac: bigint[][],
  s: number[],
  len: number,
  P: bigint,
  leadf: bigint
): bigint[] {
  let res: bigint[] = [1n];
  for (let i = 0; i < len; i++) {
    if (s[i]! < 0) continue;
    res = zpSmod(intPolyMul(res, liftedFac[s[i]!]!), P);
  }
  return zpSmod(
    res.map((c) => c * leadf),
    P
  );
}

/**
 * `fmpz_poly_factor_zassenhaus_recombination_with_prune`
 * (factor_zassenhaus_recombination.c:156-233).
 *
 * @param liftedFac - monic factors of `F/lc(F)` modulo `P`
 * @param F - the primitive squarefree polynomial being factored
 * @param P - the lifting modulus `p^a`
 * @param Z - degree pruning data, or `null` for the unpruned variant
 */
function zassenhausRecombination(
  liftedFac: bigint[][],
  F: bigint[],
  P: bigint,
  Z: ZassenhausPrune | null
): bigint[][] {
  const r = liftedFac.length;
  const subset: number[] = [];
  for (let k = 0; k < r; k++) subset.push(k);

  const out: bigint[][] = [];
  let f = F;
  let len = r;

  for (let k = 1; k <= Math.floor(len / 2); k++) {
    zassenhausSubsetFirst(subset, len, k);
    for (;;) {
      if (Z !== null) {
        let total = 0;
        for (let i = 0; i < len; i++) {
          if (subset[i]! >= 0) total += liftedFac[subset[i]!]!.length - 1;
        }
        if (!zassenhausPruneDegreeIsPossible(Z, total)) {
          if (!zassenhausSubsetNext(subset, len)) break;
          continue;
        }
      }

      const tryme = intPolyPrimitive(
        zassenhausSubsetProduct(liftedFac, subset, len, P, f[f.length - 1]!)
      )[1];

      const div = tryme.length > 0 ? intPolyQuoRem(f, tryme) : null;
      if (div !== null && zpNormalise(div[1]).length === 0) {
        out.push(tryme);
        f = div[0];
        len -= k;
        if (!zassenhausSubsetNextDisjoint(subset, len + k)) break;
      } else {
        if (!zassenhausSubsetNext(subset, len)) break;
      }
    }
  }

  if (f.length - 1 > 0) out.push(f);
  return out;
}

/* ------------------------------------------------------------------ *
 * Double + exponent arithmetic used by fmpz_poly_CLD_bound
 * (fmpz_poly/CLD_bound.c and fmpz_poly/evaluate_horner_d_2exp.c).
 *
 * Upstream is explicitly inexact here -- the CLD bound is computed in C
 * `double`s -- so we mirror it with JavaScript numbers.  Only the SPEED of
 * van Hoeij depends on it: `max{B_1(r), B_2(r)} * N` is a valid bound for
 * every `r > 0`, so an over-estimate merely discards more columns, and the
 * factorisation itself is certified by exact trial division.
 * ------------------------------------------------------------------ */

/** `ldexp` without the overflow of a single `Math.pow(2, e)`. */
function ldexpD(m: number, e: number): number {
  let r = m;
  let k = e;
  while (k > 1000) {
    r *= 2 ** 1000;
    k -= 1000;
  }
  while (k < -1000) {
    r *= 2 ** -1000;
    k += 1000;
  }
  return r * 2 ** k;
}

/** `frexp`: `x = m * 2^e` with `|m|` in `[0.5, 1)` (and `m = 0` for `x = 0`). */
function frexpD(x: number): [number, number] {
  if (x === 0 || !Number.isFinite(x)) return [x, 0];
  const a = Math.abs(x);
  let e = Math.floor(Math.log2(a)) + 1;
  let m = ldexpD(a, -e);
  while (m >= 1) {
    m /= 2;
    e++;
  }
  while (m < 0.5) {
    m *= 2;
    e--;
  }
  return [x < 0 ? -m : m, e];
}

/** `fmpz_get_d`. */
function fmpzGetD(x: bigint): number {
  return Number(x);
}

/** `fmpz_get_d_2exp`: `x ~ m * 2^e` with `|m|` in `[0.5, 1)`. */
function fmpzGetD2exp(x: bigint): [number, number] {
  if (x === 0n) return [0, 0];
  const neg = x < 0n;
  const a = neg ? -x : x;
  const e = a.toString(2).length;
  const shift = e - 53;
  const top = shift > 0 ? a >> BigInt(shift) : a << BigInt(-shift);
  const m = Number(top) / 9007199254740992; // 2^53
  return [neg ? -m : m, e];
}

/** `fmpz_set_d`: truncation towards zero. */
function fmpzSetD(d: number): bigint {
  if (!Number.isFinite(d)) throw new ArithmeticError('CLD bound overflowed a double');
  if (Math.abs(d) < 1) return 0n;
  const neg = d < 0;
  const [m, e] = frexpD(Math.abs(d));
  const mant = BigInt(Math.floor(m * 9007199254740992)); // m * 2^53, an integer
  const sh = e - 53;
  const r = sh >= 0 ? mant << BigInt(sh) : mant >> BigInt(-sh);
  return neg ? -r : r;
}

/** `fmpz_set_d_2exp` (fmpz/set.c). */
function fmpzSetD2exp(m: number, exp: number): bigint {
  const [mm, e2] = frexpD(m);
  const e = exp + e2;
  if (e >= 53) return fmpzSetD(mm * 9007199254740992) << BigInt(e - 53);
  if (e < 0) return 0n;
  return fmpzSetD(ldexpD(mm, e));
}

/** `d_polyval` (double_extras.h:40-50). */
function dPolyval(poly: number[], x: number): number {
  const len = poly.length;
  if (len === 0) return 0;
  let t = poly[len - 1]!;
  for (let i = len - 2; i >= 0; i--) t = poly[i]! + x * t;
  return t;
}

/**
 * `_fmpz_poly_evaluate_horner_d_2exp2_precomp`
 * (fmpz_poly/evaluate_horner_d_2exp.c:114-149), including the
 * `ADJUSTMENT_DELAY = 16` delayed normalisation.
 */
function hornerD2exp(polyM: number[], polyE: number[], d: number, dexp: number): [number, number] {
  const n = polyM.length;
  if (n === 0) return [0, 0];
  if (d === 0) return [polyM[0]!, polyE[0]!];

  const ADJ = 16;
  let [xm, xe] = frexpD(d);
  xe += dexp;

  let sm = polyM[n - 1]!;
  let se = polyE[n - 1]!;

  for (let i = n - 2; i >= 0; i--) {
    // dpe_mul
    sm = sm * xm;
    se = se + xe;

    if (polyM[i] !== 0) {
      // dpe_add
      const tm = polyM[i]!;
      const te = polyE[i]!;
      const diff = se - te;
      if (sm === 0) {
        sm = tm;
        se = te;
      } else if (diff >= 0) {
        if (diff <= 53 + ADJ) sm = sm + ldexpD(tm, -diff);
      } else {
        const nd = -diff;
        if (nd > 53 + ADJ) {
          sm = tm;
          se = te;
        } else {
          sm = tm + ldexpD(sm, -nd);
          se = te;
        }
      }
    }

    if (i % ADJ === 0) {
      const [am, ae] = frexpD(sm);
      sm = am;
      se += ae;
    }
  }

  const [am, ae] = frexpD(sm);
  return [am, se + ae];
}

/** `_d_cmp_2exp` (fmpz_poly/CLD_bound.c:31-68). */
function dCmp2exp(a: number, aExp: number, b: number, bExp: number): number {
  const log2 = (n: number) => Math.log(n) / Math.log(2);

  if (aExp === 0) {
    if (bExp === 0) {
      if (a > 1.5 * b) return 2;
      if (b > 1.5 * a) return -2;
      return a >= b ? 1 : -1;
    }
    const t = 1 + Math.trunc(log2(a));
    if (t >= bExp + 2) return 2;
    if (bExp >= t + 2) return -2;
    return dCmp2exp(a / 4, 0, b * 2 ** (bExp - 2), 0);
  } else if (bExp === 0) {
    return -dCmp2exp(b, bExp, a, aExp);
  } else {
    if (aExp >= bExp + 2) return 2;
    if (bExp >= aExp + 2) return -2;
    if (aExp >= bExp) return dCmp2exp(a, aExp - bExp, b, 0);
    return -dCmp2exp(b, bExp - aExp, a, 0);
  }
}

/**
 * `fmpz_poly_CLD_bound` (fmpz_poly/CLD_bound.c:70-263).
 *
 * Given `f = a_0 + ... + a_N x^N` and `0 <= n < N`, minimise
 * `max{B_1(r), B_2(r)}` over `r > 0`, where
 * `B_1(r) = (|a_0| + ... + |a_n| r^n) / r^{n+1}` and
 * `B_2(r) = (|a_{n+1}| r^{n+1} + ... + |a_N| r^N) / r^{n+1}`,
 * and return `N` times that value: a bound for the `n`-th coefficient of the
 * "coefficient of the logarithmic derivative" `f g' / g` of any factor `g`.
 */
function fmpzPolyCLDBound(f: bigint[], n: number): bigint {
  const CLD_EPS = 0.00000001;
  const flen = f.length;

  // lo(x) = |a_n| x + |a_{n-1}| x^2 + ... + |a_0| x^{n+1}
  const loArr: bigint[] = new Array(n + 2).fill(0n);
  for (let i = 1; i <= n + 1; i++) loArr[i] = absBig(f[n + 1 - i] ?? 0n);
  const lo = zpNormalise(loArr);
  // hi(x) = |a_{n+1}| + ... + |a_N| x^{N-n-1}
  const hi = zpNormalise(f.slice(n + 1).map(absBig));

  let sizeF = 0;
  for (const c of f) sizeF = Math.max(sizeF, fmpzBits(c));

  const loD = lo.map(fmpzGetD);
  const hiD = hi.map(fmpzGetD);
  const loM: number[] = [];
  const loE: number[] = [];
  for (const c of lo) {
    const [m, e] = fmpzGetD2exp(c);
    loM.push(m);
    loE.push(e);
  }
  const hiM: number[] = [];
  const hiE: number[] = [];
  for (const c of hi) {
    const [m, e] = fmpzGetD2exp(c);
    hiM.push(m);
    hiE.push(e);
  }

  const fudge = (flen - 1) * (1.0 + flen * 2 ** -50);

  let rpow = 0.0;
  let step = 1.0;
  let rexp = 0;
  let tooMuch = false;

  // The bound is valid for every r > 0, so a bailout after many refinement
  // steps is safe; upstream's loop is unbounded.
  for (let iter = 0; iter < 5000; iter++) {
    const rbits = Math.abs(rpow);
    const maxExp = rbits * hi.length + sizeF + 1;

    let hiEval: number;
    let loEval: number;
    let hiExp: number;
    let loExp: number;

    if (rbits > 200 || rexp !== 0) {
      const r = 2 ** rpow;
      [hiEval, hiExp] = hornerD2exp(hiM, hiE, r, rexp);
      [loEval, loExp] = hornerD2exp(loM, loE, 1 / r, -rexp);
    } else if (maxExp > 950 || tooMuch) {
      const r = 2 ** rpow;
      [hiEval, hiExp] = hornerD2exp(hiM, hiE, r, 0);
      [loEval, loExp] = hornerD2exp(loM, loE, 1 / r, 0);
    } else {
      const r = 2 ** rpow;
      hiEval = dPolyval(hiD, r);
      loEval = dPolyval(loD, 1 / r);
      hiExp = 0;
      loExp = 0;
    }

    if (hiExp === 0 && loExp === 0) {
      if (1.5 * loEval < hiEval) {
        if (step >= 0.0) step = -step / 2.0;
        rpow += step;
      } else if (loEval > 1.5 * hiEval) {
        if (step < 0.0) step = -step / 2.0;
        rpow += step;
      } else if (Number.isNaN(hiEval) || Number.isNaN(loEval)) {
        tooMuch = true;
      } else {
        const maxEval = hiEval > loEval ? hiEval : loEval;
        return fmpzSetD(maxEval * fudge);
      }
    } else {
      let cmp: number;
      if (rbits > 200 || rexp !== 0) {
        const l2 = (x: number) => Math.log(x) / Math.log(2);
        if (hiExp + l2(hiEval) > 1.01 * (loExp + l2(loEval))) cmp = 2;
        else if (loExp + l2(loEval) > 1.01 * (hiExp + l2(hiEval))) cmp = -2;
        else if (hiExp + l2(hiEval) >= loExp + l2(loEval)) cmp = 1;
        else cmp = -1;
      } else {
        cmp = dCmp2exp(hiEval, hiExp, loEval, loExp);
      }

      if (Math.abs(step) < CLD_EPS) cmp = cmp === 2 ? 1 : -1;

      if (cmp === 2) {
        if (step >= 0.0) step = -step / 2.0;
        rpow += step;
      } else if (cmp === -2) {
        if (step < 0.0) step = -step / 2.0;
        rpow += step;
      } else if (cmp === 1) {
        return fmpzSetD2exp(hiEval * fudge, hiExp);
      } else {
        return fmpzSetD2exp(loEval * fudge, loExp);
      }
    }

    if (rpow > 1000.0) {
      rpow -= 1000.0;
      rexp += 1000;
    } else if (rpow < -1000.0) {
      rpow += 1000.0;
      rexp -= 1000;
    }
  }

  // Bailout: `sum |a_i|` times N is a (crude but valid) bound for r = 1.
  let s = 0n;
  for (const c of f) s += absBig(c);
  return s * BigInt(flen - 1) + 1n;
}

/**
 * `fmpz_poly_divlow_smodp` (fmpz_poly/divlow_smodp.c:15-70): the low `n`
 * coefficients of `f / g` modulo `p`, in the symmetric range.
 */
function fmpzPolyDivlowSmodp(f: bigint[], g: bigint[], p: bigint, n: number): bigint[] {
  let i = 0;
  while (g[i] === 0n) i++;
  const zeroes = i;

  const tf: bigint[] = new Array(n + zeroes).fill(0n);
  for (let j = 0; j < Math.min(f.length, n + zeroes); j++) tf[j] = f[j]!;

  const c0 = g[zeroes]! >= 0n ? g[zeroes]! : g[zeroes]! + p;
  const [d, cinv] = gcdInvMod(c0, p);
  if (d !== 1n) {
    throw new ArithmeticError('Exception (fmpz_poly_divlow_smodp). Impossible inverse.');
  }

  const res: bigint[] = new Array(n).fill(0n);
  for (let k = 0; k < n; i++, k++) {
    res[k] = fmpzSmod(tf[i]! * cinv, p);
    const m = Math.min(g.length - zeroes, n - k);
    for (let j = 0; j < m; j++) {
      tf[i + j] = fmpzSmod(tf[i + j]! - g[zeroes + j]! * res[k]!, p);
    }
  }
  return res;
}

/**
 * `fmpz_poly_divhigh_smodp` (fmpz_poly/divhigh_smodp.c:15-59): the high `n`
 * coefficients of `f / g` modulo `p`, in the symmetric range.
 */
function fmpzPolyDivhighSmodp(
  f: bigint[],
  fLen: number,
  g: bigint[],
  p: bigint,
  n: number
): bigint[] {
  const lenG = g.length;
  const tf: bigint[] = new Array(fLen).fill(0n);
  for (let j = 0; j < Math.min(f.length, fLen); j++) tf[j] = f[j]!;

  const [d, cinv] = gcdInvMod(g[lenG - 1]!, p);
  if (d !== 1n) {
    throw new ArithmeticError('Exception (fmpz_poly_divhigh_smodp). Impossible inverse.');
  }

  const res: bigint[] = new Array(n).fill(0n);
  let start = 0;
  for (let k = n - 1, i = fLen - lenG; k >= 0; i--, k--) {
    if (i < fLen - n) start++;
    res[k] = fmpzSmod(tf[i + lenG - 1]! * cinv, p);
    for (let j = start; j < lenG; j++) {
      tf[i + j] = fmpzSmod(tf[i + j]! - g[j]! * res[k]!, p);
    }
  }
  return res;
}

/**
 * `_fmpz_poly_factor_CLD_mat` (fmpz_poly_factor/CLD_mat.c:19-136).
 *
 * Returns an `(r+1) x 2k` matrix whose column `j` holds the `j`-th coefficient
 * of the logarithmic derivative `f g_i' / g_i` (mod `P`) for each lifted factor
 * `g_i`, plus the CLD bound for that column in the last row; and the number of
 * usable columns (those whose bound is small enough compared to `P`).
 */
function fmpzPolyFactorCLDMat(
  f: bigint[],
  liftedFac: bigint[][],
  P: bigint,
  k: number
): { data: bigint[][]; numDataCols: number } {
  const r = liftedFac.length;
  const bitR = Math.max(r, 20);
  const flen = f.length;

  const res: bigint[][] = [];
  for (let i = 0; i <= r; i++) res.push(new Array<bigint>(2 * k).fill(0n));

  for (let i = 0; i < k; i++) {
    res[r]![i] = fmpzPolyCLDBound(f, i);
    res[r]![2 * k - i - 1] = fmpzPolyCLDBound(f, flen - i - 2);
  }

  const bound = fmpzBits(P) - bitR - Math.floor(bitR / 2);
  const sqLen = BigInt(Math.trunc(Math.sqrt(flen)));

  let loN = 0;
  for (; loN < k; loN++) {
    if (fmpzBits(res[r]![loN]! * sqLen) > bound) break;
  }
  let hiN = 0;
  for (; hiN < k; hiN++) {
    if (fmpzBits(res[r]![2 * k - hiN - 1]! * sqLen) > bound) break;
  }

  if (loN > 0) {
    for (let i = 0; i < r; i++) {
      const fac = liftedFac[i]!;
      let zeroes = 0;
      while (fac[zeroes] === 0n) zeroes++;

      const truncLen = Math.min(fac.length, loN + zeroes + 1);
      const truncFac = fac.slice(0, truncLen);
      const gd = derivativeZ(truncFac);
      const gcld = mullowZ(f, gd, loN + zeroes);
      const col = fmpzPolyDivlowSmodp(gcld, truncFac, P, loN);
      for (let j = 0; j < loN; j++) res[i]![j] = col[j]!;
    }
  }

  if (hiN > 0) {
    const truncF = f.slice(flen - hiN);
    for (let i = 0; i < r; i++) {
      const fac = liftedFac[i]!;
      const len = fac.length - hiN - 1;
      let g: bigint[];
      if (len < 0) {
        g = new Array<bigint>(-len).fill(0n).concat(fac);
      } else {
        g = fac.slice(len);
      }
      const gd = derivativeZ(g);
      const gcld = intPolyMul(truncF, gd);
      const gcldLen = truncF.length + gd.length - 1;
      const col = fmpzPolyDivhighSmodp(gcld, gcldLen, g, P, hiN);
      for (let j = 0; j < hiN; j++) res[i]![loN + j] = col[j]!;
    }
  }

  if (hiN > 0) {
    for (let i = 0; i < hiN; i++) res[r]![loN + i] = res[r]![2 * k - hiN + i]!;
  }

  return { data: res, numDataCols: loN + hiN };
}

/** `fmpz_poly_derivative`, without stripping trailing zeros. */
function derivativeZ(a: bigint[]): bigint[] {
  const out: bigint[] = [];
  for (let i = 1; i < a.length; i++) out.push(a[i]! * BigInt(i));
  return out;
}

/** `fmpz_poly_mullow`: the product truncated to `n` coefficients. */
function mullowZ(a: bigint[], b: bigint[], n: number): bigint[] {
  const out: bigint[] = new Array(n).fill(0n);
  for (let i = 0; i < a.length && i < n; i++) {
    for (let j = 0; j < b.length && i + j < n; j++) {
      out[i + j] = out[i + j]! + a[i]! * b[j]!;
    }
  }
  return out;
}

/**
 * `fmpz_mat_next_col_van_hoeij` (fmpz_mat/next_col_van_hoeij.c:35-101).
 *
 * Appends the scaled knapsack column `col` (with the new relation row) to `M`.
 * Returns `false` when running LLL would not be justified yet.
 */
function fmpzMatNextColVanHoeij(
  M: bigint[][],
  P: bigint,
  col: bigint[],
  exp: number,
  uExp: number
): boolean {
  const r = col.length;
  const bitR = Math.max(r, 20);
  const s = M.length;

  let k = fmpzBits(P) - bitR - Math.floor(bitR / 2);

  // check if LLL justified
  if (k < exp + bitCountU(r + 1)) return false;

  k -= uExp; // we want this many bits beyond the radix point

  let x: bigint[];
  let pTrunc: bigint;
  if (k >= 0) {
    x = col.map((c) => tdivQ2exp(c, k));
    pTrunc = tdivQ2exp(P, k);
  } else {
    x = col.map((c) => c << BigInt(-k));
    pTrunc = P << BigInt(-k);
  }

  // y = U x, where U is the combinatorial (first r columns) part of M
  const y: bigint[] = new Array(s).fill(0n);
  for (let j = 0; j < s; j++) {
    let acc = 0n;
    const row = M[j]!;
    for (let i = 0; i < r; i++) acc += row[i]! * x[i]!;
    y[j] = fmpzSmod(tdivQ2exp(acc, uExp), pTrunc);
  }

  // resize M: a new zero row on top, and a new column
  const c = M[0]!.length;
  for (let j = s - 1; j >= 0; j--) M[j + 1] = M[j]!.concat([y[j]!]);
  M[0] = new Array<bigint>(c).fill(0n).concat([pTrunc]);

  return true;
}

/**
 * LLL reduction of the rows of `B` with removal, `delta = 0.99`.
 *
 * This is `fmpz_lll_wrapper_with_removal_knapsack`
 * (fmpz_lll/wrapper_with_removal_knapsack.c:17-40).  FLINT reduces with a
 * chain of floating-point implementations (`fmpz_lll_d_with_removal_knapsack`,
 * `fmpz_lll_d_heuristic_with_removal`, `fmpz_lll_mpf_with_removal`) whose
 * output it then *verifies* against the exact rational predicate
 * `fmpz_mat_is_reduced_with_removal` / `gr_mat_is_row_lll_reduced_with_removal_naive`
 * (fmpz_mat/is_reduced_with_removal.c:18-70, gr_mat/is_lll_reduced.c:18-100).
 * We implement that exact contract directly with the integral Gram-Schmidt LLL
 * (Cohen, *A Course in Computational Algebraic Number Theory*, Algorithm 2.6.7)
 * at `delta = fl->delta = 0.99`, so the result satisfies the predicate FLINT
 * checks by construction.
 *
 * @see Deviation: Polynomial Roots and Factorization
 *
 * The removal rule is upstream's (`fmpz_lll/lll_d.c:392-406`): scanning from
 * the last row backwards, a row is dropped while `||b_i^*||^2 / 2 > gs_B`.
 * Since that implies `||b_i^*||^2 > gs_B`, no lattice vector of squared norm
 * at most `gs_B` can involve the dropped rows, which is exactly the property
 * van Hoeij needs.
 *
 * @param B - the basis, one lattice vector per row; reduced in place
 * @param gsB - the removal bound
 * @returns the number of rows to keep
 */
function lllWithRemovalKnapsack(B: bigint[][], gsB: bigint): number {
  const d = B.length;
  if (d === 0) return 0;

  const dot = (u: bigint[], v: bigint[]): bigint => {
    let s = 0n;
    for (let i = 0; i < u.length; i++) s += u[i]! * v[i]!;
    return s;
  };

  const dd: bigint[] = new Array(d + 1).fill(0n);
  dd[0] = 1n;
  const lam: bigint[][] = [];
  for (let i = 0; i < d; i++) lam.push(new Array<bigint>(d).fill(0n));

  dd[1] = dot(B[0]!, B[0]!);
  if (dd[1] === 0n) {
    throw new ArithmeticError('LLL: the van Hoeij lattice basis is linearly dependent');
  }

  // round to nearest integer (ties away from -infinity), exactly
  const roundDiv = (a: bigint, b: bigint): bigint => {
    const num = 2n * a + b;
    const den = 2n * b;
    let q = num / den;
    if (num % den !== 0n && num < 0n !== den < 0n) q -= 1n;
    return q;
  };

  const RED = (kk: number, l: number): void => {
    const lkl = lam[kk]![l]!;
    const dl1 = dd[l + 1]!;
    if (2n * absBig(lkl) <= dl1) return;
    const q = roundDiv(lkl, dl1);
    const bk = B[kk]!;
    const bl = B[l]!;
    for (let i = 0; i < bk.length; i++) bk[i] = bk[i]! - q * bl[i]!;
    lam[kk]![l] = lkl - q * dl1;
    for (let i = 0; i < l; i++) lam[kk]![i] = lam[kk]![i]! - q * lam[l]![i]!;
  };

  let k = 1;
  let kmax = 0;

  while (k < d) {
    if (k > kmax) {
      kmax = k;
      for (let j = 0; j <= k; j++) {
        let u = dot(B[k]!, B[j]!);
        for (let i = 0; i < j; i++) {
          u = (dd[i + 1]! * u - lam[k]![i]! * lam[j]![i]!) / dd[i]!;
        }
        if (j < k) lam[k]![j] = u;
        else dd[k + 1] = u;
      }
      if (dd[k + 1] === 0n) {
        throw new ArithmeticError('LLL: the van Hoeij lattice basis is linearly dependent');
      }
    }

    RED(k, k - 1);

    // Lovasz condition with delta = 99/100
    const lamk = lam[k]![k - 1]!;
    if (100n * dd[k + 1]! * dd[k - 1]! < 99n * dd[k]! * dd[k]! - 100n * lamk * lamk) {
      // SWAPI(k)
      const tmpRow = B[k]!;
      B[k] = B[k - 1]!;
      B[k - 1] = tmpRow;
      for (let j = 0; j <= k - 2; j++) {
        const t = lam[k]![j]!;
        lam[k]![j] = lam[k - 1]![j]!;
        lam[k - 1]![j] = t;
      }
      const BB = (dd[k - 1]! * dd[k + 1]! + lamk * lamk) / dd[k]!;
      for (let i = k + 1; i <= kmax; i++) {
        const t = lam[i]![k]!;
        lam[i]![k] = (dd[k + 1]! * lam[i]![k - 1]! - lamk * t) / dd[k]!;
        lam[i]![k - 1] = (BB * t + lamk * lam[i]![k]!) / dd[k + 1]!;
      }
      dd[k] = BB;
      k = Math.max(1, k - 1);
    } else {
      for (let l = k - 2; l >= 0; l--) RED(k, l);
      k++;
    }
  }

  // removal: ||b_i^*||^2 = dd[i+1]/dd[i]
  let newd = d;
  for (let i = d - 1; i >= 0; i--) {
    if (dd[i + 1]! > 2n * gsB * dd[i]!) newd--;
    else break;
  }
  return newd;
}

/**
 * `fmpz_mat_col_partition` (fmpz_mat/col_partition.c:62-133) with
 * `short_circuit = 1`: partition the columns of `M` into classes of equal
 * columns, numbered from 1; return the number of classes, or `0` if there are
 * more classes than `M` has rows.
 *
 * Upstream numbers the classes in the order produced by sorting a cheap hash of
 * each column; we number them by first occurrence.  Only the numbering differs
 * -- the partition itself, and hence every trial factor built from it, is the
 * same (the trial factors are sorted by degree straight afterwards).
 */
function fmpzMatColPartition(part: number[], M: bigint[][], nrows: number): number {
  const c = part.length;
  const keys: string[] = [];
  for (let j = 0; j < c; j++) {
    const col: string[] = [];
    for (let i = 0; i < M.length; i++) col.push(M[i]![j]!.toString(36));
    keys.push(col.join(','));
  }

  const seen = new Map<string, number>();
  let p = 0;
  for (let j = 0; j < c; j++) {
    let id = seen.get(keys[j]!);
    if (id === undefined) {
      p++;
      if (p > nrows) return 0;
      id = p;
      seen.set(keys[j]!, id);
    }
    part[j] = id;
  }
  return p;
}

/**
 * `fmpz_poly_factor_van_hoeij_check_if_solved`
 * (fmpz_poly_factor/van_hoeij_check_if_solved.c:29-142).
 *
 * If the combinatorial part of `M` already describes a 0-1 basis of the true
 * factor combinations, build the factors and verify them by exact division;
 * otherwise return `null`.
 */
function vanHoeijCheckIfSolved(
  M: bigint[][],
  liftedFac: bigint[][],
  f: bigint[],
  P: bigint,
  lc: bigint
): bigint[][] | null {
  const r = liftedFac.length;
  const U = M.map((row) => row.slice(0, r));
  const part: number[] = new Array(r).fill(0);

  const numFacs = fmpzMatColPartition(part, U, M.length);
  if (numFacs === 0 || numFacs > r) return null;

  if (numFacs === 1) {
    // f is irreducible
    return [f];
  }

  // there is a potential 0-1 basis, so make the potential factors
  const trial: bigint[][] = [];
  let tempLc = lc;
  for (let i = 1; i <= numFacs; i++) {
    let prod: bigint[] = [tempLc];
    for (let j = 0; j < r; j++) {
      if (part[j] === i) prod = zpSmod(intPolyMul(prod, liftedFac[j]!), P);
    }
    if (prod.length === 0) return null;
    tempLc = absBig(intPolyContent(prod));
    if (tempLc === 0n) return null;
    trial.push(prod.map((cc) => cc / tempLc));
  }

  // sort factors by length
  trial.sort((a, b) => a.length - b.length);

  // trial divide potential factors
  let fCopy = f;
  let remaining = numFacs;
  let i = 0;
  for (; i < trial.length && remaining > 1; i++) {
    const div = intPolyQuoRem(fCopy, trial[i]!);
    if (div !== null && zpNormalise(div[1]).length === 0) {
      fCopy = div[0];
      remaining--;
    } else {
      return null;
    }
  }

  if (remaining === 1) {
    const out = trial.slice(0, i);
    out.push(fCopy);
    return out;
  }
  return null;
}

/**
 * `_heuristic_van_hoeij_starting_precision` (factor_van_hoeij.c:24-42).
 */
function heuristicVanHoeijStartingPrecision(f: bigint[], r: number, p: bigint): number {
  const leadB = fmpzPolyCLDBound(f, f.length - 2);
  const trailB = fmpzPolyCLDBound(f, 0);
  const minB = Math.min(fmpzBits(leadB), fmpzBits(trailB));
  // C truncates the inner expression to slong before dividing by log(p)
  const inner = Math.trunc((2.5 * r + minB) * Math.LN2 + Math.log(f.length) / 2.0);
  return Math.trunc(inner / Math.log(Number(p)));
}

/** Smallest `a` with `p^a >= B` (`fmpz_clog_ui`). */
function clogUi(B: bigint, p: bigint): number {
  if (B <= 1n) return 0;
  let a = 0;
  let q = 1n;
  while (q < B) {
    q *= p;
    a++;
  }
  return a;
}

/**
 * `fmpz_poly_factor_van_hoeij` (fmpz_poly_factor/factor_van_hoeij.c:62-232).
 *
 * Recombination by LLL on the knapsack lattice of the lifted modular factors,
 * which is polynomial time where plain Zassenhaus is exponential.
 *
 * @param fac - the monic irreducible factors of `f` modulo `p`
 * @param f - the primitive squarefree polynomial to factor
 * @param p - the factorisation prime
 */
function fmpzPolyFactorVanHoeij(fac: bigint[][], f: bigint[], p: bigint): bigint[][] {
  const r = fac.length;
  const bitR = Math.max(r, 20);

  // set to identity, prescaled by 2^U_exp
  const uExp = bitCountU(bitR);
  let M: bigint[][] = [];
  for (let i = 0; i < r; i++) {
    const row = new Array<bigint>(r).fill(0n);
    row[i] = 1n << BigInt(uExp);
    M.push(row);
  }

  // compute Mignotte bound
  let B = fmpz_poly_factor_mignotte(f);
  B = absBig(B * f[f.length - 1]!) * 2n + 1n;
  let a = clogUi(B, p);

  // compute heuristic starting precision
  a = Math.min(a, heuristicVanHoeijStartingPrecision(f, r, p));
  if (a < 1) a = 1;

  let lifted = liftModularFactors(f, fac, p, a);

  // compute bound
  const gsB = BigInt(r + 1) << BigInt(2 * uExp);

  const N = f.length - 1;
  const sqN = BigInt(Math.trunc(Math.sqrt(N)));
  const lc = f[N]!;

  let henselLoops = 0;
  let P = p ** BigInt(a);

  for (;;) {
    const solved = vanHoeijCheckIfSolved(M, lifted, f, P, lc);
    if (solved !== null) return solved;

    let numCoeffs: number;
    if (henselLoops < 3 && 3 * r > N + 1) numCoeffs = r > 200 ? 50 : 30;
    else numCoeffs = 10;

    numCoeffs = Math.min(numCoeffs, Math.floor((N + 1) / 2));
    let prevNumCoeffs = 0;

    do {
      const { data, numDataCols } = fmpzPolyFactorCLDMat(f, lifted, P, numCoeffs);

      for (let nextCol = prevNumCoeffs; nextCol < numDataCols - prevNumCoeffs; nextCol++) {
        // we alternate taking columns from the right and left
        const diff = nextCol - prevNumCoeffs;
        const altCol =
          diff % 2 === 0
            ? prevNumCoeffs + diff / 2
            : numDataCols - prevNumCoeffs - Math.floor((diff + 1) / 2);

        const boundSum = data[r]![altCol]! * sqN;
        const worstExp = fmpzBits(boundSum);

        const col: bigint[] = [];
        for (let i = 0; i < r; i++) col.push(data[i]![altCol]!);

        if (fmpzMatNextColVanHoeij(M, P, col, worstExp, uExp)) {
          const numRows = lllWithRemovalKnapsack(M, gsB);
          M = M.slice(0, numRows);

          const done = vanHoeijCheckIfSolved(M, lifted, f, P, lc);
          if (done !== null) return done;
        }
      }

      prevNumCoeffs = numCoeffs;
      numCoeffs = Math.min(2 * numCoeffs, Math.floor((N + 1) / 2));
    } while (numCoeffs !== prevNumCoeffs);

    henselLoops++;
    if (henselLoops > 32) {
      throw new ArithmeticError(
        'van Hoeij recombination failed to converge (fmpz_poly_factor_van_hoeij)'
      );
    }

    a = 2 * a;
    P = p ** BigInt(a);
    lifted = liftModularFactors(f, fac, p, a);
  }
}

/**
 * Hensel-lift the modular factorisation of `f` to `p^a` and put the factors in
 * the symmetric range, matching what `_fmpz_poly_hensel_start_lift`
 * (fmpz_poly/hensel_start_lift.c:19-104) leaves in `lifted_fac`: monic `H_i`
 * with `H_1 ... H_r = f / lc(f) (mod p^a)`.
 */
function liftModularFactors(f: bigint[], fac: bigint[][], p: bigint, a: number): bigint[][] {
  const P = p ** BigInt(a);
  return henselLiftFactors(f, fac, p, a).map((h) => {
    const s = h.map((c) => fmpzSmod(c, P));
    s[s.length - 1] = 1n; // the lifted factors are monic
    return s;
  });
}

/**
 * Choose the factorisation prime and factor `f` modulo it.
 *
 * Transcription of the prime search in `_fmpz_poly_factor_zassenhaus`
 * (fmpz_poly_factor/factor_zassenhaus.c:113-166): three rounds over the primes,
 * each starting where the previous stopped, keeping the factorisation with the
 * fewest factors.  `p` must not divide the leading *or* the constant
 * coefficient, and `f` must stay squarefree modulo `p`.
 *
 * The search is unbounded, exactly as upstream's `for ( ; ; p = n_nextprime(p, 0))`
 * is: for a squarefree `f` with `f(0) != 0` all but finitely many primes work
 * (those not dividing `lc(f) f(0) disc(f)`), so it always terminates.
 */
function chooseFactorizationPrime(coeffs: bigint[]): {
  p: bigint;
  fac: bigint[][];
  Z: ZassenhausPrune;
} {
  const lenF = coeffs.length;
  const Z = zassenhausPruneSetDegree(lenF - 1);

  let r = lenF;
  let p = 2n;
  let bestP = 0n;
  let bestFac: bigint[][] = [];

  // A prime is rejected only when it divides `lc(f) * f(0) * disc(f)`, so a
  // squarefree `f` with `f(0) != 0` has at most `bits(lc f(0) disc(f))` bad
  // primes.  This is a (very generous) bound on that count, used only so that a
  // caller who violates the precondition gets an error instead of a hang;
  // upstream's loop (factor_zassenhaus.c:120) has no bound at all.
  let maxBits = 0;
  for (const c of coeffs) maxBits = Math.max(maxBits, fmpzBits(c));
  const badLimit = 1000 + 4 * lenF * (maxBits + 10);
  let bad = 0;

  for (let i = 0; i < 3; i++) {
    for (; ; p = next_prime(p)) {
      if (bad++ > badLimit) {
        throw new ValueError(
          'no factorization prime found: the polynomial is not squarefree, or f(0) = 0'
        );
      }
      const t = intPolyModP(coeffs, p);
      if (t.length === lenF && t[0] !== 0n) {
        const d = modPolyNormalize(
          t.slice(1).map((c, j) => c * BigInt(j + 1)),
          p
        );
        const g = d.length === 0 ? [0n] : modPolyGcd(t, d, p);

        if (g.length === 1 && g[0] === 1n) {
          const tempFac = modpFactorSquarefree(t, p);

          zassenhausPruneStartAddFactors(Z);
          for (const h of tempFac) zassenhausPruneAddFactor(Z, h.length - 1, 1);
          zassenhausPruneEndAddFactors(Z);

          if (tempFac.length <= r) {
            r = tempFac.length;
            bestFac = tempFac;
            bestP = p;
          }
          break;
        }
      }
    }
    p = next_prime(p);
  }

  return { p: bestP, fac: bestFac, Z };
}

/**
 * Factor a primitive squarefree integer polynomial into irreducible factors.
 *
 * ALGORITHM: `_fmpz_poly_factor_zassenhaus` with `cutoff = 8` and
 * `use_van_hoeij = 1` -- which is exactly how `fmpz_poly_factor` (the routine
 * behind Sage's `ZZ[x].factor()`) calls it, see
 * `reference/flint/src/fmpz_poly_factor/factor.c:98-104` and
 * `reference/flint/src/fmpz_poly_factor/factor_zassenhaus.c:90-210`.
 *
 * Factor `f` modulo a well-chosen prime `p`, then:
 *   - `r = 1`: `f` is irreducible over `Z`;
 *   - `r <= 8`: Hensel-lift to `p^a` with `p^a > 2|lc B_Mignotte| + 1` and run
 *     subset recombination with degree pruning;
 *   - `r > 8`: van Hoeij, i.e. LLL on the knapsack lattice of the lifted
 *     factors, which is polynomial rather than exponential in `r`.
 *
 * @param coeffs - Primitive squarefree polynomial coefficients, constant first
 * @returns Irreducible factors, each primitive with positive leading
 *   coefficient, whose product is the primitive part of `coeffs`
 */
function factorSquarefreeIntPoly(coeffs: bigint[]): bigint[][] {
  const n = coeffs.length - 1; // degree
  if (n <= 0) return coeffs.length > 0 && coeffs[0] !== 0n ? [[coeffs[0]!]] : [];
  if (n === 1) return [coeffs];

  // `fmpz_poly_factor` strips the x^k part before the squarefree
  // decomposition (factor.c:52-63); a squarefree f can only have x^1, and the
  // prime search below needs f(0) != 0.
  if (coeffs[0] === 0n) {
    const rest = factorSquarefreeIntPoly(coeffs.slice(1));
    return [[0n, 1n], ...rest];
  }

  const { p, fac, Z } = chooseFactorizationPrime(coeffs);
  const r = fac.length;
  const cutoff = 8;

  let factors: bigint[][];

  if (r === 1 && r <= cutoff) {
    // irreducible modulo p, hence irreducible over Z
    factors = [coeffs];
  } else if (r > cutoff) {
    factors = fmpzPolyFactorVanHoeij(fac, coeffs, p);
  } else {
    // bound adjustment: we multiply true factors (which might be monic) by the
    // leading coefficient of f (factor_zassenhaus.c:186-196)
    let T = fmpz_poly_factor_mignotte(coeffs);
    T = absBig(T * coeffs[n]!) * 2n + 1n;
    const a = clogUi(T, p);

    const lifted = liftModularFactors(coeffs, fac, p, a);
    const P = p ** BigInt(a);

    factors = zassenhausRecombination(lifted, coeffs, P, Z);
  }

  // Normalise: primitive, positive leading coefficient.  (FLINT keeps whatever
  // sign the recombination produced; the number of sign flips is even, so the
  // product is unchanged.)
  factors = factors.map((g) => intPolyPrimitive(g)[1]);

  // The factors must reproduce the input exactly; silently returning a wrong
  // factorization would be worse than raising.
  let check: bigint[] = [1n];
  for (const g of factors) check = intPolyMul(check, g);
  const checkPrimitive = intPolyPrimitive(check)[1];
  const inputPrimitive = intPolyPrimitive(coeffs)[1];
  if (
    checkPrimitive.length !== inputPrimitive.length ||
    checkPrimitive.some((c, i) => c !== inputPrimitive[i])
  ) {
    throw new ArithmeticError('integer polynomial factorization failed to reproduce its input');
  }

  return factors.length > 0 ? factors : [coeffs];
}

/**
 * Compute squarefree factorization of an integer polynomial.
 * Returns pairs [squarefree_factor, multiplicity].
 */
function squarefreeFactorIntPoly(coeffs: bigint[]): Array<[bigint[], number]> {
  if (coeffs.length === 0) return [];
  if (coeffs.length === 1) return [[coeffs, 1]];

  // Make primitive
  const [content, primitive] = intPolyPrimitive(coeffs);

  // For simple cases, just return the polynomial
  if (primitive.length <= 2) {
    return [[primitive, 1]];
  }

  // Compute derivative
  const deriv: bigint[] = [];
  for (let i = 1; i < primitive.length; i++) {
    deriv.push(primitive[i]! * BigInt(i));
  }

  // If derivative is zero or constant, polynomial is already squarefree (or p-th power in char p)
  if (deriv.length === 0 || (deriv.length === 1 && deriv[0] === 0n)) {
    return [[primitive, 1]];
  }

  // gcd(f, f')
  let g = intPolyGcd(primitive, deriv);

  // If gcd is constant, f is squarefree
  if (g.length <= 1) {
    return [[primitive, 1]];
  }

  // f / gcd(f, f')
  const divResult = intPolyQuoRem(primitive, g);
  if (divResult === null) {
    // Shouldn't happen for well-formed input
    return [[primitive, 1]];
  }
  let h = divResult[0];

  const result: Array<[bigint[], number]> = [];
  let i = 1;
  // Every multiplicity is at most the degree, so the loop cannot run longer
  // than that; an earlier fixed cap of 20 silently truncated the
  // decomposition of things like (x-1)^25, which then reached the Zassenhaus
  // code with a non-squarefree input.
  const maxIter = primitive.length;

  while (h.length > 1 && i <= maxIter) {
    // gcd(g, h)
    const gi = intPolyGcd(g, h);

    // h / gi
    const hDivGi = intPolyQuoRem(h, gi);
    // Both are primitive and gi divides h in Q[x], so by Gauss's lemma the
    // division is exact over Z; a failure means a broken invariant, and
    // continuing would silently return a wrong decomposition.
    if (hDivGi === null) {
      throw new ArithmeticError('squarefree decomposition: inexact division over ZZ');
    }
    const hi = hDivGi[0];

    if (hi.length > 1) {
      result.push([hi, i]);
    }

    // g = g / gi
    const gDivGi = intPolyQuoRem(g, gi);
    if (gDivGi === null) {
      throw new ArithmeticError('squarefree decomposition: inexact division over ZZ');
    }
    g = gDivGi[0];
    h = gi;
    i++;
  }

  // g might still have content
  if (g.length > 1) {
    result.push([g, i]);
  }

  // If no factors found, return original as squarefree
  if (result.length === 0) {
    return [[primitive, 1]];
  }

  return result;
}

/** Primitive GCD for squarefree factorization, delegated to FLINT. */
function intPolyGcd(a: bigint[], b: bigint[]): bigint[] {
  const g = _fmpz_poly_gcd(a, b);
  return g.length ? intPolyPrimitive(g)[1] : [1n];
}

/**
 * Factor an integer polynomial completely: squarefree decomposition followed by
 * Zassenhaus on each squarefree part.
 *
 * Returns `[content, factors]` where `factors` is an array of
 * `[irreducible_factor, multiplicity]` and `content * prod(factors^mult)` is
 * the input.
 *
 * @see Deviation: Polynomial Roots and Factorization
 */
function factorIntegerPolynomial(coeffs: bigint[]): [bigint, Array<[bigint[], number]>] {
  if (coeffs.length === 0) return [0n, []];

  // Extract content and make primitive
  const [content, primitive] = intPolyPrimitive(coeffs);

  if (primitive.length <= 1) {
    return [content, []];
  }

  // Squarefree factorization
  const sqfree = squarefreeFactorIntPoly(primitive);

  // Factor each squarefree part
  const result: Array<[bigint[], number]> = [];

  for (const [sqfFactor, mult] of sqfree) {
    if (sqfFactor.length <= 1) continue;

    const irredFactors = factorSquarefreeIntPoly(sqfFactor);
    for (const irredFactor of irredFactors) {
      result.push([irredFactor, mult]);
    }
  }

  return [content, result];
}

/** Integer-ring factor dispatch, retaining only roots which coerce back to ZZ. */
function integerRootsFromFactorization(coeffs: bigint[]): Array<[bigint, number]> {
  if (coeffs.length === 0) return [];
  const content = coeffs.reduce((g, c) => gcdBigInt(g, c), 0n);
  // Sage factor() factors content even when roots ignores the constant factors.
  // The dense root path and sparse gap path remove it before reaching this helper.
  if (content > 1n) factorInteger(content);
  if (coeffs.length === 1) return [];
  const primitive = coeffs.map((c) => c / content);
  const degree = primitive.length - 1;
  const factors =
    degree < 30 || degree > 300 ? ntlIntegerFactor(primitive)[1] : pariIntegerFactor(primitive);
  const roots: Array<[bigint, number]> = [];
  for (const [f, multiplicity] of factors) {
    if (f.length === 2 && f[0]! % f[1]! === 0n) roots.push([-f[0]! / f[1]!, multiplicity]);
  }
  // Primitive integral linear factors of integral roots are monic.
  return roots.sort(([a, m], [b, n]) => m - n || (a > b ? -1 : a < b ? 1 : 0));
}

/**
 * IntegerRing._roots_univariate_polynomial: dense factorization up to degree 100,
 * then the Cucker–Koiran–Smale exponent-gap algorithm and sparse derivatives.
 * @see Reference: sage/rings/integer_ring.pyx:_roots_univariate_polynomial
 */
function findIntegerRoots(coeffs: bigint[]): Array<[bigint, number]> {
  if (coeffs.length <= 1) return [];
  if (coeffs.length <= 101) {
    const content = coeffs.reduce((g, c) => gcdBigInt(g, c), 0n);
    return integerRootsFromFactorization(coeffs.map((c) => c / content));
  }
  let valuation = 0;
  while (coeffs[valuation] === 0n) valuation++;
  const roots: Array<[bigint, number]> = valuation ? [[0n, valuation]] : [];
  const p = coeffs.slice(valuation);
  if (p.length === 1) return roots;
  const e = p.flatMap((c, i) => (c === 0n ? [] : [i]));
  if (e.length === p.length) return roots.concat(integerRootsFromFactorization(p));
  const content = p.reduce((g, c) => gcdBigInt(g, c), 0n);
  const c = e.map((i) => p[i]! / content);
  const k = e.length;
  const nbits = (n: bigint): number => (n < 0n ? -n : n).toString(2).length;
  const block = (start: number, end: number): bigint[] => {
    const f = Array<bigint>(e[end - 1]! - e[start]! + 1).fill(0n);
    for (let j = start; j < end; j++) f[e[j]! - e[start]!] = c[j]!;
    return f;
  };
  let maxBits = nbits(c[0]!);
  let first = 0;
  let g: bigint[] = [];
  for (let i = 1; i < k; i++) {
    if (e[i]! - e[i - 1]! > maxBits) {
      g = _fmpz_poly_gcd(g, block(first, i));
      if (g.length === 1 && g[0] === 1n) break;
      first = i;
      maxBits = nbits(c[i]!);
    } else maxBits = Math.max(maxBits, nbits(c[i]!));
  }
  if (!g.length) return roots.concat(integerRootsFromFactorization(p.map((x) => x / content)));
  g = _fmpz_poly_gcd(g, block(first, k));
  let cc = c;
  let ee = e;
  let m1 = 0,
    m2 = 0;
  let b1 = true,
    b2 = true;
  for (let i = 0; i < k; i++) {
    let s1 = 0n,
      s2 = 0n;
    for (let j = 0; j < k - i; j++) {
      if (b1) s1 += cc[j]!;
      if (b2) s2 += ee[j]! % 2 ? -cc[j]! : cc[j]!;
    }
    if (b1 && s1 !== 0n) {
      m1 = i;
      b1 = false;
    }
    if (b2 && s2 !== 0n) {
      m2 = i;
      b2 = false;
    }
    if (!b1 && !b2) break;
    ee = ee.slice(1).map((x) => x - ee[0]! - 1);
    cc = ee.map((x, j) => BigInt(x + 1) * cc[j + 1]!);
  }
  if (m1 > 0) roots.push([1n, m1]);
  if (m2 > 0) roots.push([-1n, m2]);
  roots.push(...integerRootsFromFactorization(g).filter(([r]) => r > 1n || r < -1n));
  return roots;
}

/**
 * Clear denominators of a polynomial over QQ, returning integer coefficients
 * and the LCM of denominators.
 */
function clearDenominators<C extends RingElement>(poly: Polynomial<C>): [bigint[], bigint] {
  // Extract rational coefficients as [numerator, denominator] pairs
  const rats: Array<[bigint, bigint]> = poly.coeffs.map((c) => {
    // Handle Rational class
    if ('numerator' in c && 'denominator' in c) {
      const r = c as unknown as { numerator: bigint; denominator: bigint };
      return [r.numerator, r.denominator];
    }
    if ('numer' in c && 'denom' in c) {
      const r = c as unknown as { numer: bigint; denom: bigint };
      return [r.numer, r.denom];
    }
    // Handle _numerator and _denominator (private fields)
    if ('_numerator' in c && '_denominator' in c) {
      const r = c as unknown as { _numerator: bigint; _denominator: bigint };
      return [r._numerator, r._denominator];
    }
    // Assume integer
    const val = 'value' in c ? (c as { value: bigint }).value : BigInt(c.toString());
    return [val, 1n];
  });

  // Compute LCM of denominators
  let lcmDenom = 1n;
  for (const [_, d] of rats) {
    lcmDenom = lcm(lcmDenom, d);
  }

  // Multiply each coefficient by lcm / denom
  const intCoeffs: bigint[] = rats.map(([n, d]) => n * (lcmDenom / d));

  return [intCoeffs, lcmDenom];
}

/**
 * LCM of two bigints.
 */
function lcm(a: bigint, b: bigint): bigint {
  if (a === 0n || b === 0n) return 0n;
  const absA = a < 0n ? -a : a;
  const absB = b < 0n ? -b : b;
  return (absA / gcdBigInt(absA, absB)) * absB;
}

/** QQ._factor_univariate_polynomial delegates rational factorization to PARI. */
function findRationalRoots<C extends RingElement>(poly: Polynomial<C>): Array<[C, number]> {
  const factors = pariRationalFactor(clearDenominators(poly));
  const roots: Array<[C, number]> = [];
  for (const [f, multiplicity] of factors) {
    if (f.length === 2) {
      const root = new Rational(-f[0]!, f[1]!);
      roots.push([poly.parent.base_ring.__call__(root) as C, multiplicity]);
    }
  }
  return roots;
}

/**
 * Base class for polynomial rings.
 */
export interface PolynomialRingBase<C extends RingElement> {
  readonly base_ring: CoefficientRing<C>;
  readonly variable_name: string;
  zero(): Polynomial<C>;
  one(): Polynomial<C>;
  gen(n?: unknown): Polynomial<C>;
  __call__(x?: unknown): Polynomial<C>;
}

/**
 * Test-only surface for the FLINT transcriptions used by `ZZ[x].factor()`.
 *
 * Not part of the public API: exported purely so that
 * `polynomial_factorization.test.ts` can exercise the individual routines
 * (subset enumeration, degree pruning, the CLD bound, LLL, van Hoeij) against
 * brute force and against upstream's own invariants.  Do NOT re-export this
 * from `rings/index.ts`.
 */
export const _zz_factor_internal = {
  fmpz_poly_factor_mignotte,
  fmpzPolyCLDBound,
  fmpzPolyDivlowSmodp,
  fmpzPolyDivhighSmodp,
  fmpzPolyFactorCLDMat,
  fmpzMatNextColVanHoeij,
  fmpzMatColPartition,
  lllWithRemovalKnapsack,
  vanHoeijCheckIfSolved,
  fmpzPolyFactorVanHoeij,
  liftModularFactors,
  chooseFactorizationPrime,
  zassenhausSubsetFirst,
  zassenhausSubsetNext,
  zassenhausSubsetNextDisjoint,
  zassenhausRecombination,
  zassenhausPruneSetDegree,
  zassenhausPruneStartAddFactors,
  zassenhausPruneAddFactor,
  zassenhausPruneEndAddFactors,
  zassenhausPruneDegreeIsPossible,
  factorSquarefreeIntPoly,
  factorIntegerPolynomial,
  isqrtFloor,
  fmpzSmod,
  clogUi,
};

/** Scalar index adapters for polynomial_element.pyx and Cython's integer arguments. */
function polynomialDerivativeVariable(x: unknown): string {
  if (x === undefined || x === null) return 'None';
  if (typeof x === 'boolean') return x ? 'True' : 'False';
  if (Array.isArray(x))
    return `[${x.map((c) => (typeof c === 'string' ? `'${c.replaceAll('\\', '\\\\').replaceAll("'", "\\'")}'` : polynomialDerivativeVariable(c))).join(', ')}]`;
  return polynomialScalarRepr(x);
}
function polynomialScalarType(x: unknown): string {
  if (x === null || x === undefined) return 'NoneType';
  if (Array.isArray(x)) return 'list';
  if (typeof x === 'string') return 'str';
  if (typeof x === 'number') return Number.isInteger(x) ? 'sage.rings.integer.Integer' : 'float';
  if (typeof x === 'bigint' || x instanceof Integer) return 'sage.rings.integer.Integer';
  if (typeof x === 'boolean') return 'bool';
  if (x instanceof Rational) return 'sage.rings.rational.Rational';
  if (x instanceof FiniteFieldElement)
    return 'sage.rings.finite_rings.element_pari_ffelt.FiniteFieldElement_pari_ffelt';
  return typeof x === 'object' ? x.constructor.name : typeof x;
}
function polynomialScalarRepr(x: unknown): string {
  if (typeof x === 'number') {
    if (Number.isNaN(x)) return 'nan';
    if (x === Infinity) return 'inf';
    if (x === -Infinity) return '-inf';
  }
  return String(x);
}
/** arith/long.pxd:pyobject_to_long, plus the distinct generated Cython converters. */
function polynomialInteger(
  x: unknown,
  kind: 'index' | 'long' | 'int' | 'ssize' | 'unsigned'
): bigint {
  const strict = kind === 'index' || kind === 'ssize';
  let n: bigint;
  if (typeof x === 'bigint') n = x;
  else if (x instanceof Integer) n = x.value;
  else if (typeof x === 'boolean') n = BigInt(x);
  else if (x instanceof Rational) {
    if (strict && x.denominator !== 1n)
      throw new TypeError(`unable to convert rational ${x} to an integer`);
    n = x.numerator / x.denominator;
  } else if (typeof x === 'number') {
    if (strict && !Number.isInteger(x))
      throw new TypeError("'float' object cannot be interpreted as an integer");
    if (Number.isNaN(x)) throw new ValueError('cannot convert float NaN to integer');
    if (!Number.isFinite(x)) throw new OverflowError('cannot convert float infinity to integer');
    n = BigInt(Math.trunc(x));
  } else if (!strict && x instanceof FiniteFieldElement) {
    n = x._integer_();
  } else if (
    x !== null &&
    typeof x === 'object' &&
    'toBigInt' in x &&
    typeof x.toBigInt === 'function'
  ) {
    n = x.toBigInt();
  } else {
    throw new TypeError(
      strict
        ? `'${polynomialScalarType(x)}' object cannot be interpreted as an integer`
        : 'an integer is required'
    );
  }
  if (kind === 'unsigned') {
    if (n < 0n) throw new OverflowError("can't convert negative value to unsigned long");
    if (n >= 1n << 64n)
      throw new OverflowError('Python int too large to convert to C unsigned long');
  } else {
    if (n < -(1n << 63n) || n >= 1n << 63n) {
      const label =
        kind === 'ssize'
          ? 'Python int too large to convert to C ssize_t'
          : kind === 'index' &&
              (typeof x === 'bigint' || x instanceof Integer || typeof x === 'number')
            ? 'Sage Integer too large to convert to C long'
            : 'Python int too large to convert to C long';
      throw new OverflowError(label);
    }
    if (kind === 'int' && (n < -(1n << 31n) || n >= 1n << 31n))
      throw new OverflowError('value too large to convert to int');
  }
  return n;
}
function polynomialBackend<C extends RingElement>(
  base: CoefficientRing<C>
): 'integer' | 'rational' | 'word' | 'binary' | 'extension' | 'large' | 'generic' {
  if (base === (RDF as unknown)) return 'generic';
  if (isIntegerRing(base)) return 'integer';
  if (isRationalField(base)) return 'rational';
  const p = getRingCharacteristic(base),
    zero = base.zero();
  if (p === null || p <= 1n) return 'generic';
  if ('value' in zero) return p === 2n ? 'binary' : p < 1n << 63n ? 'word' : 'large';
  if ('degree' in base && typeof base.degree === 'number' && base.degree > 1) return 'extension';
  return 'generic';
}
function polynomialScalarSign(x: unknown, comparison: '<' | '>'): -1 | 0 | 1 | null {
  let n: bigint | number;
  if (typeof x === 'bigint' || typeof x === 'number') n = x;
  else if (typeof x === 'boolean') n = Number(x);
  else if (x instanceof Integer) n = x.value;
  else if (x instanceof Rational) n = x.numerator;
  else if (x instanceof FiniteFieldElement) n = x.isZero() ? 0n : 1n;
  else if (
    x !== null &&
    typeof x === 'object' &&
    'toBigInt' in x &&
    typeof x.toBigInt === 'function'
  )
    n = x.toBigInt();
  else
    throw new TypeError(
      `'${comparison}' not supported between instances of '${polynomialScalarType(x)}' and 'int'`
    );
  return n < 0 ? -1 : n > 0 ? 1 : n === 0 || n === 0n ? 0 : null;
}

/** Preserve scalar-parent arithmetic in reverse's degree + 1 expression. */
function polynomialScalarAddOne(x: unknown): unknown {
  if (
    x !== null &&
    typeof x === 'object' &&
    !(x instanceof Rational) &&
    !(x instanceof Integer) &&
    'add' in x &&
    typeof x.add === 'function'
  )
    return x.add(1n);
  return null;
}

/** CoercionModel canonical parent selection, restricted to implemented coefficient rings. */
function polynomialCommonBase(
  a: CoefficientRing<RingElement>,
  b: CoefficientRing<RingElement>
): CoefficientRing<RingElement> | null {
  if (a === b) return a;
  if (a instanceof FractionField_generic || b instanceof FractionField_generic) {
    const frac = (a instanceof FractionField_generic ? a : b) as FractionField_generic<RingElement>;
    const other = a instanceof FractionField_generic ? b : a;
    if (other instanceof PolynomialRing && polynomialCanCoerce(other, frac)) return other;
    if (polynomialCanCoerce(frac.ring(), other)) return frac;
    if (
      other instanceof PolynomialRing &&
      other.base_ring instanceof PolynomialRing &&
      polynomialCanCoerce(frac.ring(), other.base_ring)
    )
      return new PolynomialRing(frac, other.variable_name);
    return null;
  }
  // Polynomial constructions are functorial in the coefficient base. Constants
  // embed into an existing polynomial base; incompatible variables do not.
  if (a instanceof PolynomialRing || b instanceof PolynomialRing) {
    if (a instanceof PolynomialRing && b instanceof PolynomialRing) {
      if (a.variable_name !== b.variable_name) {
        if (polynomialCanCoerce(a.base_ring, b)) return a;
        if (polynomialCanCoerce(b.base_ring, a)) return b;
        const variables = (ring: CoefficientRing<RingElement>): string[] => {
          const names: string[] = [];
          while (ring instanceof PolynomialRing) {
            names.push(ring.variable_name);
            ring = ring.base_ring;
          }
          return names;
        };
        const av = variables(a),
          bv = variables(b);
        if (
          a.base_ring instanceof PolynomialRing &&
          av.length > bv.length &&
          bv.every((v, i) => v === av[av.length - bv.length + i])
        ) {
          const coefficientParent = polynomialCommonBase(a.base_ring, b);
          if (coefficientParent instanceof PolynomialRing)
            return new PolynomialRing(coefficientParent, a.variable_name);
        }
        if (
          b.base_ring instanceof PolynomialRing &&
          bv.length > av.length &&
          av.every((v, i) => v === bv[bv.length - av.length + i])
        ) {
          const coefficientParent = polynomialCommonBase(a, b.base_ring);
          if (coefficientParent instanceof PolynomialRing)
            return new PolynomialRing(coefficientParent, b.variable_name);
        }
        return null;
      }
      const common = polynomialCommonBase(a.base_ring, b.base_ring);
      return common === null
        ? null
        : common === a.base_ring
          ? a
          : common === b.base_ring
            ? b
            : new PolynomialRing(common, a.variable_name);
    }
    const poly = (a instanceof PolynomialRing ? a : b) as PolynomialRing<RingElement>;
    const scalar = a instanceof PolynomialRing ? b : a;
    const common = polynomialCommonBase(poly.base_ring, scalar);
    return common === null
      ? null
      : common === poly.base_ring
        ? poly
        : new PolynomialRing(common, poly.variable_name);
  }
  if (a === RDF || b === RDF) {
    const other = a === RDF ? b : a;
    return other === RDF || ['Integer Ring', 'Rational Field'].includes(other.toString())
      ? RDF
      : null;
  }
  if (a.toString() === 'Integer Ring') return b.toString() === 'Integer Ring' ? a : b;
  if (b.toString() === 'Integer Ring') return a;
  if (a.toString() === 'Rational Field' || b.toString() === 'Rational Field')
    return a.toString() === 'Rational Field' && b.toString() === 'Rational Field' ? a : null;
  const a0 = a.zero(),
    b0 = b.zero();
  const left = a0 instanceof GF2Element ? new PrimeField(2n).zero() : a0;
  const right = b0 instanceof GF2Element ? new PrimeField(2n).zero() : b0;
  const finite = (
    x: RingElement
  ): x is IntegerMod | PrimeFieldElement | LegacyPrimeElement | FiniteFieldElement =>
    x instanceof IntegerMod ||
    x instanceof PrimeFieldElement ||
    x instanceof LegacyPrimeElement ||
    x instanceof FiniteFieldElement;
  if (finite(left) && finite(right)) {
    try {
      const [coerced] = canonicalFiniteOperands(left, right, '+');
      return coerced.parent === left.parent
        ? a
        : coerced.parent === right.parent
          ? b
          : coerced.parent;
    } catch (e) {
      if (e instanceof TypeError) return null;
      throw e;
    }
  }
  return null;
}
/** PolynomialRing._coerce_map_from_: same variable and a canonical coefficient map. */
function polynomialCommonOperands<C extends RingElement>(
  a: Polynomial<C>,
  b: Polynomial<C>,
  operation: string
): [Polynomial<C>, Polynomial<C>] {
  const common = polynomialCommonBase(a.parent, b.parent);
  if (
    common === null &&
    (operation === 'quo_rem' ||
      operation === 'gcd' ||
      operation === 'xgcd' ||
      operation === 'sylvester_matrix' ||
      operation === 'resultant' ||
      operation === 'multiplication_trunc' ||
      operation === 'pow')
  )
    throw new TypeError(
      `no common canonical parent for objects with parents: '${a.parent}' and '${b.parent}'`
    );
  if (common === null)
    throw new TypeError(
      `unsupported operand parent(s) for ${operation}: '${a.parent}' and '${b.parent}'`
    );
  const parent = common as unknown as PolynomialRingBase<C>;
  const convert = (f: Polynomial<C>): Polynomial<C> => {
    if (f.parent === parent) return f;
    const target = parent.base_ring;
    if (f.parent.variable_name !== parent.variable_name)
      return new Polynomial([target.__call__(f)], parent);
    const binary = target.zero() instanceof GF2Element;
    const coeffs = f.coeffs.map((c) => {
      const value =
        c instanceof GF2Element
          ? c.toBigInt()
          : binary &&
              (c instanceof Integer ||
                c instanceof IntegerMod ||
                c instanceof PrimeFieldElement ||
                c instanceof LegacyPrimeElement)
            ? c.value
            : c;
      return target.__call__(value) as C;
    });
    return new Polynomial(coeffs, parent as PolynomialRingBase<C>);
  };
  return [convert(a), convert(b)];
}

function polynomialCanCoerce(
  target: CoefficientRing<RingElement>,
  source: CoefficientRing<RingElement>
): boolean {
  if (target === source) return true;
  if (target instanceof FractionField_generic)
    return polynomialCanCoerce(
      target.ring(),
      source instanceof FractionField_generic ? source.ring() : source
    );
  if (target instanceof PolynomialRing) {
    if (source instanceof PolynomialRing && target.variable_name === source.variable_name)
      return polynomialCanCoerce(target.base_ring, source.base_ring);
    return polynomialCanCoerce(target.base_ring, source);
  }
  if (source instanceof PolynomialRing) return false;
  return polynomialCommonBase(target, source) === target;
}
function polynomialRealBase(base: CoefficientRing<RingElement>): boolean {
  if (base instanceof PolynomialRing) return polynomialRealBase(base.base_ring);
  return base === RDF || base.toString() === 'Integer Ring' || base.toString() === 'Rational Field';
}
function polynomialRealCoefficient(c: RingElement): number | null {
  if (c instanceof RealDoubleElement) return c.value;
  if (c instanceof Rational) return c.toNumber();
  if (c instanceof Integer) return Number(c.value);
  if (c instanceof Polynomial) {
    if (!polynomialRealBase(c.parent.base_ring)) return null;
    if (!c.coeffs.slice(1).every((x: RingElement) => polynomialRealCoefficient(x) === 0))
      return null;
    return c.coeffs.length ? polynomialRealCoefficient(c.coeffs[0]!) : 0;
  }
  return null;
}

/** Common numerator storage for the native QQ polynomial backend. */
function polynomialRationalData(f: { readonly coeffs: readonly unknown[] }): [bigint[], bigint] {
  const cs = f.coeffs as unknown as Rational[];
  let den = 1n;
  for (const c of cs) den = (den / gcdBigInt(den, c.denominator)) * c.denominator;
  return [cs.map((c) => c.numerator * (den / c.denominator)), den];
}

function polynomialElementParent(x: unknown): CoefficientRing<RingElement> {
  if (
    typeof x === 'bigint' ||
    x instanceof Integer ||
    (typeof x === 'number' && Number.isInteger(x))
  )
    return ZZ as unknown as CoefficientRing<RingElement>;
  if (x instanceof RealDoubleElement || isFractionElement(x)) return x.parent;
  if (x instanceof Rational) return QQ as unknown as CoefficientRing<RingElement>;
  if (
    x instanceof Polynomial ||
    x instanceof IntegerMod ||
    x instanceof PrimeFieldElement ||
    x instanceof LegacyPrimeElement ||
    x instanceof FiniteFieldElement ||
    x instanceof GF2Element
  )
    return x.parent as CoefficientRing<RingElement>;
  throw new AttributeError(`'${polynomialScalarType(x)}' object has no attribute 'parent'`);
}
function polynomialTruth(x: unknown): boolean {
  if (x === undefined || x === null) return false;
  if (typeof x === 'boolean') return x;
  if (typeof x === 'number' || typeof x === 'bigint') return x !== 0 && x !== 0n;
  if (typeof x === 'string' || Array.isArray(x)) return x.length > 0;
  if (typeof (x as RingElement).isZero === 'function') return !(x as RingElement).isZero();
  return true;
}
function polynomialSylvesterEntries<C extends RingElement>(
  a: Polynomial<C>,
  b: Polynomial<C>
): C[][] {
  const m = a.degree();
  const n = b.degree();
  const size = m + n;
  const matrix: C[][] = [];

  // Initialize matrix with zeros
  for (let i = 0; i < size; i++) {
    const row: C[] = [];
    for (let j = 0; j < size; j++) {
      row.push(a.parent.base_ring.zero() as C);
    }
    matrix.push(row);
  }

  // Fill in rows for f (n rows)
  // Row i contains coefficients of x^i * f
  for (let i = 0; i < n; i++) {
    for (let j = 0; j <= m; j++) {
      matrix[i]![i + (m - j)] = a.getCoeff(j);
    }
  }

  // Fill in rows for g (m rows)
  // Row n+i contains coefficients of x^i * g
  for (let i = 0; i < m; i++) {
    for (let j = 0; j <= n; j++) {
      matrix[n + i]![i + (n - j)] = b.getCoeff(j);
    }
  }

  return matrix;
}

/** Internal canonical-parent helpers shared with the polynomial constructor. */
export const _polynomial_coercion = { canCoerce: polynomialCanCoerce };

function polynomialRdfBase(base: CoefficientRing<RingElement>): boolean {
  return base === RDF || (base instanceof PolynomialRing && polynomialRdfBase(base.base_ring));
}

/** The original method checks scalar protocols before coefficient conversion. */
function polynomialPseudoType(x: unknown): string {
  if (
    x instanceof IntegerMod ||
    x instanceof PrimeFieldElement ||
    x instanceof LegacyPrimeElement ||
    x instanceof GF2Element
  ) {
    const p =
      x instanceof GF2Element
        ? 2n
        : x instanceof IntegerMod
          ? x.modulus
          : getRingCharacteristic(x.parent as CoefficientRing<RingElement>)!;
    return `sage.rings.finite_rings.integer_mod.IntegerMod_${p < 46341n ? 'int' : p < 2147483648n ? 'int64' : 'gmp'}`;
  }
  return polynomialScalarType(x);
}
/** Parent.__contains__: explicit conversion followed by equality, not a canonical map. */
function polynomialCoefficientContains(base: CoefficientRing<RingElement>, x: unknown): boolean {
  try {
    if (polynomialElementParent(x) === base) return true;
    if (
      x instanceof Polynomial &&
      x.degree() > 0 &&
      (x.parent.base_ring === base || base.zero() instanceof GF2Element)
    )
      return false;
    let value = x instanceof Polynomial && x.degree() <= 0 ? x.getCoeff(0) : x;
    if (!(base instanceof PolynomialRing))
      while (value instanceof Polynomial && value.degree() <= 0) value = value.getCoeff(0);
    const converted =
      base.zero() instanceof GF2Element
        ? new PolynomialRing(base).__call__(value).getCoeff(0)
        : base.__call__(value);
    if (x instanceof Polynomial) {
      if (
        base.zero() instanceof IntegerMod &&
        (base.zero() as IntegerMod).modulus === 1n &&
        x.parent.base_ring.zero() instanceof IntegerMod &&
        (x.parent.base_ring.zero() as IntegerMod).modulus !== 1n
      )
        return false;
      return x.eq(converted);
    }
    if (converted instanceof Polynomial) return converted.eq(x);
    if (converted instanceof Integer) {
      if (x instanceof Rational) return x.denominator === 1n && x.numerator === converted.value;
      if (x instanceof Integer) return x.value === converted.value;
      return (x as RingElement).eq(converted);
    }
    if (converted instanceof GF2Element) {
      const proxy = new PrimeField(2n).__call__(converted.toBigInt());
      return proxy.eq(x);
    }
    return converted.eq(x as RingElement);
  } catch (e) {
    if (
      e instanceof TypeError ||
      e instanceof ValueError ||
      e instanceof ArithmeticError ||
      e instanceof ZeroDivisionError
    )
      return false;
    throw e;
  }
}
function polynomialDomain(r: CoefficientRing<RingElement>): boolean {
  return r instanceof PolynomialRing
    ? polynomialDomain(r.base_ring)
    : r.is_field?.() === true || r.toString() === 'Integer Ring';
}
function polynomialCoefficientPower(
  x: RingElement,
  exponent: number,
  context?: CoefficientRing<RingElement>
): RingElement {
  if (x instanceof Polynomial && exponent < 0) {
    const backend = polynomialBackend(x.parent.base_ring);
    if ((backend === 'generic' || backend === 'large') && x.degree() <= 0) {
      const inverse = polynomialCoefficientPower(x.getCoeff(0), exponent);
      return x.parent.__call__(inverse);
    }
    if (!polynomialDomain(x.parent))
      throw new TypeError(`unsupported operand parent(s) for /: '${x.parent}' and '${x.parent}'`);
    return FractionField(x.parent as PolynomialRing<RingElement>).__call__(1n, x.pow(-exponent));
  }
  const pow = (x as unknown as { pow?: (e: bigint) => RingElement }).pow;
  if (typeof pow !== 'function') {
    if (!context)
      throw new NotImplementedError('SAGE_NOT_IMPLEMENTED: coefficient power in pseudo-division');
    if (exponent === 0) return context.one();
    let n = Math.abs(exponent),
      a = x;
    if (exponent < 0) {
      if (a.eq(context.one())) return a;
      const inv = (a as unknown as { inv?: () => RingElement }).inv;
      if (!inv) throw new NotImplementedError('SAGE_NOT_IMPLEMENTED: custom coefficient inverse');
      a = inv.call(a);
    }
    while (n % 2 === 0) {
      a = a.mul(a);
      n /= 2;
    }
    let result = a;
    n = Math.floor(n / 2);
    while (n) {
      a = a.mul(a);
      if (n % 2) result = a.mul(result);
      n = Math.floor(n / 2);
    }
    return result;
  }
  return pow.call(x, BigInt(exponent));
}
/** Scalar multiplication retains scalar parent and operand order in coercion errors. */
function polynomialScalarProduct(
  f: Polynomial<RingElement>,
  c: RingElement,
  left: boolean
): Polynomial<RingElement> {
  let cp: CoefficientRing<RingElement>;
  try {
    cp = polynomialElementParent(c);
  } catch (e) {
    if (!(e instanceof AttributeError) || c.constructor !== f.parent.base_ring.one().constructor)
      throw e;
    cp = f.parent.base_ring;
  }
  let common = c instanceof Integer ? f.parent : polynomialCommonBase(f.parent, cp);
  // ModuleAction forms the pushout of the actor and module base. A finite
  // field actor can base-extend QQ-polynomials by reducing their coefficients;
  // QQ itself cannot act in characteristic p (there is no connecting map).
  if (common === null && cp.is_field?.() && (getRingCharacteristic(cp) ?? 0n) > 0n) {
    const extend = (base: CoefficientRing<RingElement>): CoefficientRing<RingElement> | null => {
      if (base.toString() === 'Rational Field') return cp;
      if (base instanceof PolynomialRing) {
        const inner = extend(base.base_ring);
        return inner === null ? null : new PolynomialRing(inner, base.variable_name);
      }
      return null;
    };
    common = extend(f.parent);
  }
  if (!(common instanceof PolynomialRing))
    throw new TypeError(
      `unsupported operand parent(s) for *: '${left ? cp : f.parent}' and '${left ? f.parent : cp}'`
    );
  const ff = common === f.parent ? f : common.__call__(f);
  const cc = common.__call__(c).getCoeff(0);
  if (!ff.coeffs.length && polynomialBackend(common.base_ring) === 'generic') return ff;
  return new Polynomial(
    ff.coeffs.map((v) => (left ? cc.mul(v) : v.mul(cc))),
    common
  );
}

export function generic_power_trunc<C extends RingElement>(
  p: Polynomial<C>,
  n: bigint,
  prec: number
): Polynomial<C> {
  if (n < 0n) throw new ValueError('n must be a nonnegative integer');
  if (prec <= 0) return p.parent.zero();
  if (n === 0n) return p.parent.one();
  if (n === 1n) return p.truncate(prec);
  if (n === 2n) return p._mul_trunc_(p, prec);
  if (n === 3n) return p._mul_trunc_(p, prec)._mul_trunc_(p, prec);
  const a = p.truncate(prec),
    aa = a._mul_trunc_(a, prec);
  if (aa.eq(a)) return a;
  let apow = aa,
    i = 1n;
  while (!(n & (1n << i))) {
    apow = apow._mul_trunc_(apow, prec);
    i++;
  }
  let power = apow;
  i++;
  if (n & 1n) power = power._mul_trunc_(a, prec);
  while (1n << i <= n) {
    apow = apow._mul_trunc_(apow, prec);
    if (n & (1n << i)) power = power._mul_trunc_(apow, prec);
    i++;
  }
  return power;
}

function do_karatsuba_product<C extends RingElement>(
  a: readonly C[],
  b: readonly C[],
  threshold: number
): C[] {
  const n = a.length,
    m = b.length;
  if (!n || !m) return [];
  if (Math.min(n, m) <= Math.max(1, threshold)) return do_schoolbook_product(a, b, -1);
  if (n !== m) {
    // The longer input is split without reversing coefficient multiplication.
    const length = Math.min(n, m),
      out = do_karatsuba_product(a.slice(0, length), b.slice(0, length), threshold);
    for (let offset = length; offset < Math.max(n, m); offset += length) {
      const carry =
        n > m
          ? do_karatsuba_product(a.slice(offset, offset + length), b, threshold)
          : do_karatsuba_product(a, b.slice(offset, offset + length), threshold);
      for (let i = 0; i < carry.length; i++) {
        const k = offset + i;
        if (k < out.length) out[k] = out[k]!.add(carry[i]!);
        else out.push(carry[i]!);
      }
    }
    return out;
  }
  if (n === 2) {
    const bd = a[0]!.mul(b[0]!),
      ac = a[1]!.mul(b[1]!);
    return [bd, a[0]!.add(a[1]!).mul(b[0]!.add(b[1]!)).sub(ac).sub(bd), ac];
  }
  const e = Math.floor(n / 2),
    ac = do_karatsuba_product(a.slice(e), b.slice(e), threshold),
    bd = do_karatsuba_product(a.slice(0, e), b.slice(0, e), threshold);
  const A = a.slice(e),
    B = b.slice(e);
  for (let i = 0; i < e; i++) {
    A[i] = A[i]!.add(a[i]!);
    B[i] = B[i]!.add(b[i]!);
  }
  const middle = do_karatsuba_product(A, B, threshold);
  for (let i = 0; i < middle.length; i++)
    middle[i] = middle[i]!.sub(i < bd.length ? ac[i]!.add(bd[i]!) : ac[i]!);
  for (let i = 0; i < e - 1; i++) bd[e + i] = bd[e + i]!.add(middle[i]!);
  bd.push(middle[e - 1]!);
  for (let i = 0; i < ac.length - e; i++) ac[i] = ac[i]!.add(middle[e + i]!);
  return bd.concat(ac);
}

function polynomialModularPower<C extends RingElement>(
  f: Polynomial<C>,
  e: bigint,
  modulus: unknown,
  genericFallback = false
): Polynomial<C> | FractionElement<C> {
  const backend = polynomialBackend(f.parent.base_ring),
    base = f.parent.base_ring;
  const generic = genericFallback || backend === 'generic' || backend === 'large';
  if (backend === 'binary' && !generic && (e >= 1n << 63n || e < -(1n << 63n)))
    return polynomialModularPower(f, e, modulus, true);
  if (generic && (f.degree() <= 0 || e < 0n)) return f.pow(e);
  const zeroMod = !polynomialTruth(modulus);
  if (generic && zeroMod) return f.pow(e);
  if (
    backend === 'extension' &&
    !generic &&
    (modulus === 0n ||
      modulus === 0 ||
      (typeof modulus === 'object' &&
        modulus !== null &&
        'isZero' in modulus &&
        typeof modulus.isZero === 'function' &&
        modulus.isZero()))
  )
    throw new ZeroDivisionError('modulus must be nonzero');
  let m: Polynomial<C>;
  if (modulus instanceof Polynomial) {
    const [a, b] = polynomialCommonOperands(f, modulus, 'pow');
    if (a !== f)
      return polynomialModularPower(a, e, b, genericFallback) as Polynomial<C> | FractionElement<C>;
    m = b as Polynomial<C>;
  } else {
    if (
      backend === 'extension' &&
      !generic &&
      !(
        typeof modulus === 'bigint' ||
        modulus instanceof Integer ||
        modulus instanceof Rational ||
        modulus instanceof IntegerMod ||
        modulus instanceof PrimeFieldElement ||
        modulus instanceof LegacyPrimeElement ||
        modulus instanceof GF2Element ||
        modulus instanceof FiniteFieldElement ||
        modulus instanceof RealDoubleElement
      )
    )
      throw new AttributeError(
        `'${polynomialScalarType(modulus)}' object has no attribute 'is_zero'`
      );
    const primitiveFloat = typeof modulus === 'number' && !Number.isInteger(modulus),
      invalid = typeof modulus === 'string' || Array.isArray(modulus);
    const input =
      typeof modulus === 'boolean'
        ? BigInt(modulus)
        : primitiveFloat
          ? RDF.__call__(modulus as number)
          : modulus;
    const rightParent = invalid ? null : polynomialElementParent(input),
      zero = base.zero();
    const common =
      rightParent === null ||
      (modulus instanceof IntegerMod &&
        modulus.modulus === 1n &&
        zero instanceof IntegerMod &&
        zero.modulus !== 1n)
        ? null
        : polynomialCommonBase(f.parent, rightParent);
    if (common === null) {
      const label =
        invalid || primitiveFloat
          ? `<class '${polynomialScalarType(modulus)}'>`
          : String(rightParent);
      throw new TypeError(
        `no common canonical parent for objects with parents: '${f.parent}' and '${label}'`
      );
    }
    if (!(common instanceof PolynomialRing))
      throw new TypeError(
        `Cannot convert ${String(common)} to sage.rings.polynomial.polynomial_element.Polynomial`
      );
    const a = common.__call__(f),
      b = common.__call__(input);
    if (a !== f) return a.pow(e, b) as Polynomial<C> | FractionElement<C>;
    m = b as Polynomial<C>;
  }
  if (generic) {
    if (
      e > 0n &&
      m.parent === f.parent &&
      m.coeffs.filter((c) => !c.isZero()).length === 1 &&
      m.leading_coefficient().eq(base.one())
    )
      return f.power_trunc(e, m.degree());
    const a = f.mod(m);
    if (m.eq(1)) return a.parent.zero();
    if (e === 0n) return a.parent.one();
    let power = a,
      n = e;
    while (!(n & 1n)) {
      power = power.mul(power).mod(m);
      n >>= 1n;
    }
    let out = power;
    n >>= 1n;
    while (n) {
      power = power.mul(power).mod(m);
      if (n & 1n) out = out.mul(power).mod(m);
      n >>= 1n;
    }
    return out;
  }
  if (backend === 'word' || backend === 'extension') f = f.mod(m);
  const large = e > 0n && e.toString(2).length >= 32;
  if ((backend === 'binary' || !large) && f.isZero())
    return e === 0n ? f.parent.one() : f.parent.zero();
  if (m.isZero()) throw new ZeroDivisionError('modulus must be nonzero');
  if ((backend === 'binary' || !large) && m.eq(1)) return f.parent.zero();
  let out: Polynomial<C>;
  const exponent =
    e === -(1n << 63n) && (backend === 'binary' || backend === 'extension') ? e : e < 0n ? -e : e;
  try {
    if (backend === 'word') {
      const p = getRingCharacteristic(base)!,
        A = extractIntegerCoeffs(f),
        M = extractIntegerCoeffs(m);
      let result: bigint[];
      if (large) {
        const minv = _nmod_poly_inv_series_newton(M.slice().reverse(), M.length, p);
        result =
          f.degree() === 1 && f.getCoeff(0).isZero() && f.leading_coefficient().eq(base.one())
            ? _nmod_poly_powmod_x_fmpz_preinv(exponent, M, minv, p)
            : _nmod_poly_powmod_fmpz_binexp_preinv(A, exponent, M, minv, p);
      } else result = _nmod_poly_powmod_ui_binexp(A, exponent, M, p);
      out = new Polynomial(
        result.map((c) => base.__call__(c)),
        f.parent
      );
    } else if (backend === 'binary') {
      const pack = (v: Polynomial<C>) =>
        new GF2X(v.coeffs.map((c) => Number((c as unknown as { value: bigint | number }).value)));
      const M = pack(m),
        A = pack(f).rem(M),
        rep = GF2X.PowerMod(A, exponent, M).rep();
      out = new Polynomial(
        rep === 0n
          ? []
          : rep
              .toString(2)
              .split('')
              .reverse()
              .map((c) => base.__call__(BigInt(c))),
        f.parent
      );
    } else {
      const B = base as unknown as FiniteFieldElement['parent'],
        F = B.modulus.coeffs.map((c) => c.value);
      const pack = (v: Polynomial<C>) =>
        v.coeffs.map((c) => (c as unknown as FiniteFieldElement).lift.coeffs.map((v) => v.value));
      const A = pack(f),
        M = pack(m);
      let result: bigint[][];
      if (m.degree() === 1 && !large) result = ZZ_pEX_power(A, exponent, F, B.characteristic);
      else if (
        large &&
        f.degree() === 1 &&
        f.getCoeff(0).isZero() &&
        f.leading_coefficient().eq(base.one())
      )
        result = ZZ_pEX_PowerXMod(exponent, M, F, B.characteristic);
      else result = ZZ_pEX_PowerMod(A, exponent, M, F, B.characteristic);
      out = new Polynomial(
        result.map((c) => B.__call__(c) as unknown as C),
        f.parent
      );
    }
  } catch (error) {
    if (
      (backend === 'binary' || backend === 'extension') &&
      error instanceof Error &&
      error.name === 'Error'
    )
      throw new NTLError(error.message);
    throw error;
  }
  if (e >= 0n) return out;
  if (!polynomialDomain(out.parent))
    throw new TypeError(`unsupported operand parent(s) for /: '${out.parent}' and '${out.parent}'`);
  return FractionField(out.parent as PolynomialRing<C>).__call__(1n, out);
}

/** Sage caches the compiled instruction graph on the polynomial object. */
const polynomialGenerators = new WeakSet<Polynomial<RingElement>>();
const polynomialIrreducibility = new WeakMap<object, boolean>();
const polynomialEvaluationCache = new WeakMap<
  Polynomial<RingElement>,
  CompiledPolynomialFunction<RingElement>
>();
function polynomialEvaluate(
  f: Polynomial<RingElement>,
  input: unknown
): RingElement | number | EvaluationMatrix {
  let base = f.parent.base_ring;
  const backend = polynomialBackend(base);
  const inputParent = (x: unknown): CoefficientRing<RingElement> | null => {
    if (typeof x === 'number') return RDF;
    if (typeof x === 'boolean') return ZZ as unknown as CoefficientRing<RingElement>;
    try {
      return polynomialElementParent(x);
    } catch (e) {
      if (e instanceof AttributeError) {
        // Explicit coefficient adapters need not expose Sage's runtime parent.
        const zero = base.zero();
        if (
          x !== null &&
          typeof x === 'object' &&
          x.constructor !== Object &&
          x.constructor === zero.constructor &&
          'isZero' in x
        )
          return base;
        return null;
      }
      throw e;
    }
  };
  const coerce = (R: CoefficientRing<RingElement>, x: unknown): RingElement => {
    if (typeof x === 'boolean') x = BigInt(x);
    const value = R.__call__(x);
    return typeof value === 'bigint' ? (new Integer(value) as unknown as RingElement) : value;
  };
  const rawParent = inputParent(input);
  const rawInteger =
    typeof input === 'bigint' || typeof input === 'boolean' || input instanceof Integer;
  if (backend === 'integer' && rawInteger)
    return coerce(
      base,
      _fmpz_poly_evaluate_fmpz(
        extractIntegerCoeffs(f),
        input instanceof Integer ? input.value : BigInt(input as bigint | boolean)
      )
    );
  if (backend === 'rational') {
    const [a, den] = polynomialRationalData(f);
    if (input instanceof Polynomial && polynomialBackend(input.parent.base_ring) === 'rational') {
      const [b, db] = polynomialRationalData(input);
      const [out, d] = _fmpq_poly_compose(a, den, b, db);
      return new Polynomial(
        out.map((c) =>
          divideCoeffs(coerce(input.parent.base_ring, c), coerce(input.parent.base_ring, d))
        ),
        input.parent
      );
    }
    if (rawInteger || input instanceof Rational) {
      const [n, d] =
        input instanceof Rational
          ? _fmpq_poly_evaluate_fmpq(a, den, input.numerator, input.denominator)
          : _fmpq_poly_evaluate_fmpz(
              a,
              den,
              input instanceof Integer ? input.value : BigInt(input as bigint | boolean)
            );
      return divideCoeffs(coerce(base, n), coerce(base, d));
    }
  }
  if (
    (backend === 'word' ||
      backend === 'extension' ||
      (backend === 'large' && base.is_field?.() !== true)) &&
    rawParent &&
    polynomialCanCoerce(base, rawParent)
  ) {
    const point = coerce(base, input);
    if (backend === 'extension') {
      const B = base as unknown as FiniteFieldElement['parent'];
      const pack = (c: RingElement) => (c as FiniteFieldElement).lift.coeffs.map((v) => v.value);
      return B.__call__(
        ZZ_pEX_eval(
          f.coeffs.map(pack),
          pack(point),
          B.modulus.coeffs.map((v) => v.value),
          B.characteristic
        )
      );
    }
    const p = getRingCharacteristic(base)!;
    const a = f.coeffs.map((c) => (c as PrimeFieldElement).value);
    const x = (point as PrimeFieldElement).value;
    return coerce(
      base,
      backend === 'word' ? _nmod_poly_evaluate_nmod(a, x, p) : ZZ_pX_evaluate(a, x, p)
    );
  }
  const generator = (p: Polynomial<RingElement>) =>
    p.degree() === 1 && p.getCoeff(0).isZero() && p.getCoeff(1).eq(p.parent.base_ring.one());
  if (backend === 'word' && rawParent && polynomialCanCoerce(f.parent, rawParent)) {
    const g = f.parent.__call__(input);
    if (generator(g)) return f;
    return new Polynomial(
      _nmod_poly_compose(
        f.coeffs.map((c) => (c as PrimeFieldElement).value),
        g.coeffs.map((c) => (c as PrimeFieldElement).value),
        getRingCharacteristic(base)!
      ).map((c) => coerce(base, c)),
      f.parent
    );
  }
  let point = input;
  if (Array.isArray(point)) {
    const args = point;
    point = args[0];
    if (args.length > 1) {
      const top = f.degree() < 0 ? base.one() : f.leading_coefficient();
      if (!(top instanceof Polynomial)) throw new TypeError('Wrong number of arguments');
      let evaluated: RingElement | number | EvaluationMatrix;
      try {
        evaluated = polynomialEvaluate(top, args.slice(1));
      } catch (e) {
        if (e instanceof TypeError) throw new TypeError('Wrong number of arguments');
        throw e;
      }
      const newBase = inputParent(evaluated)!;
      f = new Polynomial(
        f.coeffs.map((c) =>
          coerce(newBase, polynomialEvaluate(c as Polynomial<RingElement>, args.slice(1)))
        ),
        new PolynomialRing(newBase, f.parent.variable_name)
      );
      base = newBase;
    }
  }
  if (point === null || point === undefined) point = f.parent.gen();
  const matrixResult = polynomialMatrixEvaluation(f, point);
  if (matrixResult !== undefined) return matrixResult;
  let R = inputParent(point);
  if (point instanceof Polynomial && point.parent.base_ring === base) {
    if (point.is_gen() || generator(point))
      return point.parent.variable_name === f.parent.variable_name
        ? f
        : new Polynomial([...f.coeffs], point.parent);
    if (point.degree() < 0) return point.parent.__call__(f.getCoeff(0));
    if (point.degree() === 0)
      return point.parent.__call__(polynomialEvaluate(f, point.getCoeff(0)));
    if (
      point.leading_coefficient().eq(base.one()) &&
      point.coeffs.slice(0, -1).every((c) => c.isZero())
    ) {
      const cs = Array.from({ length: Math.max(0, f.degree() * point.degree() + 1) }, () =>
        base.zero()
      );
      for (let i = 0; i < f.coeffs.length; i++) cs[i * point.degree()] = f.coeffs[i]!;
      return new Polynomial(cs, point.parent);
    }
  }
  const common = R && polynomialCommonBase(base, R);
  if (!common) {
    const label =
      typeof point === 'number'
        ? "<class 'float'>"
        : R
          ? String(R)
          : `<class '${polynomialScalarType(point)}'>`;
    throw new TypeError(
      `no common canonical parent for objects with parents: '${base}' and '${label}'`
    );
  }
  const primitiveFloat =
    typeof point === 'number' && ['Integer Ring', 'Rational Field'].includes(base.toString());
  const convert = (x: unknown): RingElement | number =>
    primitiveFloat ? (typeof x === 'number' ? x : RDF.__call__(x).value) : coerce(common, x);
  point = convert(point);
  const cst = convert(f.getCoeff(0));
  const exact = (B: CoefficientRing<RingElement>): boolean =>
    B instanceof PolynomialRing ? exact(B.base_ring) : B !== RDF;
  if (
    f.degree() <= 0 ||
    (typeof point !== 'number' && exact(common) && (point as RingElement).isZero())
  )
    return cst;
  const arithmetic = {
    multiply: (a: unknown, b: unknown): RingElement | number => {
      const A = convert(a),
        B = convert(b);
      return typeof A === 'number' ? A * (B as number) : A.mul(B as RingElement);
    },
    add: (a: unknown, b: unknown): RingElement | number => {
      const A = convert(a),
        B = convert(b);
      return typeof A === 'number' ? A + (B as number) : A.add(B as RingElement);
    },
  };
  if (f.degree() < 4 || f.degree() > 50000) {
    let result: unknown = f.leading_coefficient();
    for (let i = f.degree() - 1; i >= 0; i--)
      result = arithmetic.add(arithmetic.multiply(result, point), f.getCoeff(i));
    return result as RingElement | number;
  }
  let compiled = polynomialEvaluationCache.get(f);
  if (!compiled) {
    compiled = new CompiledPolynomialFunction(f.coeffs);
    polynomialEvaluationCache.set(f, compiled);
  }
  return compiled.eval(point, arithmetic) as RingElement | number;
}

function polynomialMatrixEvaluation(
  f: Polynomial<RingElement>,
  point: unknown
): EvaluationMatrix | undefined {
  const isMatrix = (a: unknown): a is EvaluationMatrix =>
    a instanceof EvaluationGenericMatrix ||
    a instanceof EvaluationIntegerMatrix ||
    a instanceof Matrix_mod2_dense ||
    a instanceof Matrix_modn_dense;
  if (!isMatrix(point)) return undefined;
  const integerBase = {
    zero: () => new Integer(0n),
    one: () => new Integer(1n),
    __call__: (x: unknown) =>
      x instanceof Integer ? x : new Integer(ZZ.__call__(x as IntegerInput)),
    is_field: () => false,
    toString: () => 'Integer Ring',
  } as unknown as CoefficientRing<RingElement>;
  const matrixBase = (a: EvaluationMatrix): CoefficientRing<RingElement> =>
    a instanceof EvaluationIntegerMatrix
      ? integerBase
      : a instanceof Matrix_mod2_dense
        ? (GF2 as CoefficientRing<RingElement>)
        : a instanceof Matrix_modn_dense
          ? (Zmod(a.modulus) as unknown as CoefficientRing<RingElement>)
          : a.base_ring;
  const base = f.parent.base_ring,
    pointBase = matrixBase(point),
    rows = point.nrows,
    cols = point.ncols;
  const nativeDefault = (B: CoefficientRing<RingElement>): boolean => {
    if (B === RDF || ['Integer Ring', 'Rational Field'].includes(String(B))) return true;
    if (B instanceof PolynomialRing && B.base_ring.is_field?.() === true) return true;
    const zero = B.zero(),
      p = getRingCharacteristic(B);
    if (zero instanceof FiniteFieldElement) {
      const q = zero.parent.cardinality();
      return (p === 2n && q <= 65536n) || q < 256n;
    }
    return p !== null && p > 0n && p < 94906266n;
  };
  const metadata = (a: EvaluationMatrix) =>
    polynomialMatrixParents.get(a) ?? {
      base: matrixBase(a),
      generic: a instanceof EvaluationGenericMatrix,
    };
  const describe = (a: EvaluationMatrix) => {
    const P = metadata(a);
    return `Full MatrixSpace of ${a.nrows} by ${a.ncols} dense matrices over ${P.base}${P.generic && nativeDefault(P.base) ? ' (using Matrix_generic_dense)' : ''}`;
  };
  const integerZero = String(base) === 'Integer Ring' && f.getCoeff(0).isZero();
  let common = polynomialCanCoerce(pointBase, base)
    ? pointBase
    : polynomialCommonBase(base, pointBase);
  // MatrixFunctor changes the base through pushout. QuotientFunctor.merge
  // rejects a trivial quotient intersection, even when a scalar coercion
  // into the zero ring exists (categories/pushout.py).
  if (
    !polynomialCanCoerce(pointBase, base) &&
    common &&
    getRingCharacteristic(common) === 1n &&
    (getRingCharacteristic(pointBase) ?? 0n) > 0n
  )
    common = null;
  // CoercionModel.canonical_coercion has a universal Integer(0) fallback;
  // rectangular MatrixSpaces have no scalar coercion map, but accept zero.
  if ((rows !== cols && !integerZero) || !common)
    throw new TypeError(
      `no common canonical parent for objects with parents: '${base}' and '${describe(point)}'`
    );
  type Parent = { base: CoefficientRing<RingElement>; generic: boolean };
  const pointParent = metadata(point);
  const target: Parent = polynomialCanCoerce(pointBase, base)
    ? pointParent
    : { base: common, generic: false };
  const scalar = (B: CoefficientRing<RingElement>, a: unknown): RingElement => {
    const c = B.__call__(a);
    return typeof c === 'bigint' ? (new Integer(c) as unknown as RingElement) : c;
  };
  const remember = (a: EvaluationMatrix, P: Parent): EvaluationMatrix => {
    polynomialMatrixParents.set(a, P);
    return a;
  };
  const construct = (P: Parent, entries: RingElement[][]): EvaluationMatrix => {
    const zero = P.base.zero();
    // The default MeatAxe converter requires Givaro even when its parent is PARI.
    if (
      !P.generic &&
      zero instanceof FiniteFieldElement &&
      zero.parent.characteristic !== 2n &&
      zero.parent.cardinality() < 256n
    )
      throw new AttributeError(
        "'FiniteField_pari_ffelt_with_category' object has no attribute '_cache'"
      );
    if (!P.generic && String(P.base) === 'Integer Ring')
      return remember(
        new EvaluationIntegerMatrix(
          rows,
          cols,
          entries.map((row) => row.map((c) => ZZ.__call__(c as unknown as IntegerInput)))
        ),
        P
      );
    if (!P.generic && zero instanceof IntegerMod)
      return remember(
        new Matrix_modn_dense(
          rows,
          cols,
          zero.modulus,
          entries.map((row) => row.map((c) => (c as IntegerMod).value))
        ),
        P
      );
    if (!P.generic && String(P.base) === 'Finite Field of size 2')
      return remember(
        new Matrix_mod2_dense(
          rows,
          cols,
          entries.map((row) => row.map((c) => Number(ZZ.__call__(c as unknown as IntegerInput))))
        ),
        P
      );
    return remember(new EvaluationGenericMatrix(P.base, rows, cols, entries), P);
  };
  const convert = (a: EvaluationMatrix, P: Parent): EvaluationMatrix => {
    const A = metadata(a);
    if (A.base === P.base && A.generic === P.generic) return a;
    if (a instanceof EvaluationIntegerMatrix && P.base === RDF && !P.generic)
      return remember(integer_to_real_double_dense(a), P);
    return construct(
      P,
      Array.from({ length: rows }, (_, i) =>
        Array.from({ length: cols }, (_, j) => scalar(P.base, a.get(i, j)))
      )
    );
  };
  const constant = (a: unknown, P: Parent): EvaluationMatrix => {
    const c = scalar(P.base, a),
      z = scalar(P.base, 0n);
    return construct(
      P,
      Array.from({ length: rows }, (_, i) =>
        Array.from({ length: cols }, (_, j) => (i === j ? c : z))
      )
    );
  };
  const commonParent = (a: EvaluationMatrix, b: EvaluationMatrix): Parent => {
    const A = metadata(a),
      B = metadata(b);
    if (polynomialCanCoerce(A.base, B.base) && polynomialCanCoerce(B.base, A.base))
      return { base: A.base, generic: A.generic && B.generic };
    return { base: polynomialCommonBase(A.base, B.base)!, generic: false };
  };
  const arithmetic = {
    multiply: (a: unknown, b: unknown): EvaluationMatrix => {
      if (isMatrix(a) && isMatrix(b)) {
        if (a.ncols !== b.nrows)
          throw new TypeError(
            `unsupported operand parent(s) for *: '${describe(a)}' and '${describe(b)}'`
          );
        const P = commonParent(a, b),
          A = convert(a, P),
          B = convert(b, P);
        const out =
          A instanceof EvaluationIntegerMatrix
            ? A.mul(B as EvaluationIntegerMatrix)
            : A instanceof Matrix_mod2_dense
              ? A.mul(B as Matrix_mod2_dense)
              : A instanceof Matrix_modn_dense
                ? A.mul(B as Matrix_modn_dense)
                : A.mul(B as EvaluationGenericMatrix<RingElement>);
        return remember(out, P);
      }
      const matrix = (isMatrix(a) ? a : b) as EvaluationMatrix,
        c = isMatrix(a) ? b : a;
      const M = metadata(matrix),
        P = polynomialCanCoerce(M.base, base)
          ? M
          : { base: polynomialCommonBase(base, M.base)!, generic: false };
      const A = convert(matrix, P),
        C = scalar(P.base, c);
      const out =
        A instanceof EvaluationIntegerMatrix
          ? A.scalar_mul(ZZ.__call__(C as unknown as IntegerInput))
          : A instanceof Matrix_mod2_dense
            ? construct(
                P,
                Array.from({ length: rows }, (_, i) =>
                  Array.from({ length: cols }, (_, j) => scalar(P.base, A.get(i, j)).mul(C))
                )
              )
            : A instanceof Matrix_modn_dense
              ? A.scalar_mul((C as IntegerMod).value)
              : A.scalar_mul(C);
      return remember(out, P);
    },
    add: (a: unknown, b: unknown): EvaluationMatrix => {
      const matrix = (isMatrix(a) ? a : b) as EvaluationMatrix;
      if (rows !== cols && !(isMatrix(a) && isMatrix(b))) {
        const c = isMatrix(a) ? b : a;
        if (!(c instanceof Integer && c.isZero()))
          throw new TypeError(
            `unsupported operand parent(s) for +: '${isMatrix(a) ? describe(a) : base}' and '${isMatrix(b) ? describe(b) : base}'`
          );
      }
      const M = metadata(matrix),
        P =
          isMatrix(a) && isMatrix(b)
            ? commonParent(a, b)
            : polynomialCanCoerce(M.base, base)
              ? M
              : { base: polynomialCommonBase(base, M.base)!, generic: false };
      const A = isMatrix(a) ? convert(a, P) : constant(a, P),
        B = isMatrix(b) ? convert(b, P) : constant(b, P);
      const out =
        A instanceof EvaluationIntegerMatrix
          ? A.add(B as EvaluationIntegerMatrix)
          : A instanceof Matrix_mod2_dense
            ? A.add(B as Matrix_mod2_dense)
            : A instanceof Matrix_modn_dense
              ? A.add(B as Matrix_modn_dense)
              : A.add(B as EvaluationGenericMatrix<RingElement>);
      return remember(out, P);
    },
  };
  const cst = constant(f.getCoeff(0), target);
  const exact = (B: CoefficientRing<RingElement>): boolean =>
    B instanceof PolynomialRing ? exact(B.base_ring) : B !== RDF;
  const isZero = Array.from({ length: rows }, (_, i) =>
    Array.from({ length: cols }, (_, j) => scalar(pointBase, point.get(i, j)).isZero())
  ).every((row) => row.every(Boolean));
  if (f.degree() <= 0 || (exact(pointBase) && isZero)) return cst;
  if (f.degree() < 4 || f.degree() > 50000) {
    let value: unknown = f.leading_coefficient();
    for (let i = f.degree() - 1; i >= 0; i--)
      value = arithmetic.add(arithmetic.multiply(value, point), f.getCoeff(i));
    return value as EvaluationMatrix;
  }
  let compiled = polynomialEvaluationCache.get(f);
  if (!compiled) {
    compiled = new CompiledPolynomialFunction(f.coeffs);
    polynomialEvaluationCache.set(f, compiled);
  }
  return compiled.eval(point, arithmetic) as EvaluationMatrix;
}

const polynomialMatrixParents = new WeakMap<
  EvaluationMatrix,
  { base: CoefficientRing<RingElement>; generic: boolean }
>();


/** Integer object adapter for native ZZ polynomial results. */
const polynomialIntegerCoefficientRing = {
  zero: () => new Integer(0n),
  one: () => new Integer(1n),
  __call__: (x: unknown) => x instanceof Integer ? x : new Integer(ZZ.__call__(x as IntegerInput)),
  is_field: () => false,
  toString: () => 'Integer Ring',
} as unknown as CoefficientRing<RingElement>;
const polynomialIntegerNumeratorParents = new Map<string, PolynomialRing<RingElement>>();
