import { idealhnf_principal, idealadd } from '@sagemath-ts/parigp-ts/src/base4.js';
import { type NfPrimeIdeal } from '@sagemath-ts/parigp-ts/src/base2.js';
import { generic_power_pos } from '../../arith/power.js';
import { type IntegerLike, type RationalLike, toBigInt } from '../../types/coercion.js';
import { Integer } from '../integer_ring.js';
import {
  idealHNF_inv,
  idealmul,
  mat_ideal_two_elt,
  idealdiv,
  idealintersect,
  idealval,
  idealfactor,
  idealismaximal,
} from '@sagemath-ts/parigp-ts/src/base4.js';
/**
 * @module sage/rings/number_field/number_field_ideal
 * @description Ideals of number fields
 *
 * Port of: sage/rings/number_field/number_field_ideal.py
 * Reference: reference/sage/src/sage/rings/number_field/number_field_ideal.py
 *
 * Ideals in number fields are represented using Hermite Normal Form (HNF)
 * with respect to an integral basis. Operations use PARI's ideal arithmetic.
 */

import { gcd as intGcd, lcm as intLcm, is_prime } from '../../arith/misc.js';
import {
  AttributeError,
  NotImplementedError,
  TypeError,
  ValueError,
  ZeroDivisionError,
} from '../../errors.js';
import { Rational } from '../rational.js';
import type { NumberField, NumberFieldInput } from './number_field.js';
import { NumberFieldElement } from './number_field_element.js';
import { ratInverse } from './pari_nf.js';

/** Scalar generators accepted by the number-field ideal constructors. */
export type NumberFieldIdealGenerator = RationalLike | number | NumberFieldElement;
/** A generator, generator list, or existing ideal. */
export type NumberFieldIdealInput =
  | NumberFieldIdealGenerator
  | NumberFieldIdealGenerator[]
  | NumberFieldIdeal;

/**
 * Hermite Normal Form representation of an ideal.
 * The HNF is an upper triangular matrix with respect to the power basis.
 */
export interface HNFMatrix {
  /**
   * Lower-triangular `n x n` integer matrix; row `i` gives the coordinates,
   * in the integral basis of the field, of the `i`-th element of a Z-basis of
   * `denominator * I`.  This is PARI's `idealhnf` shape, so `entries[0]` is a
   * multiple of `1` and `entries[0][0]/denominator` generates `I \cap Q`.
   */
  entries: bigint[][];
  /** The denominator (for fractional ideals) */
  denominator: bigint;
}

/**
 * An ideal of a number field (or its ring of integers).
 *
 * @see Reference: sage/rings/number_field/number_field_ideal.py:NumberFieldIdeal
 * @see Deviation: Number-field ideal intersection and construction adapters
 */
export class NumberFieldIdeal {
  protected readonly _number_field: NumberField;
  protected readonly _gens: NumberFieldElement[];
  protected _cachedHNF?: HNFMatrix;
  protected _cachedPrimeData?: NfPrimeIdeal;
  protected _cachedFactorization?: Array<[NumberFieldIdeal, bigint]>;
  protected _cachedNorm?: Rational;
  protected _cachedDenominator?: NumberFieldIdeal;
  protected _cachedNumerator?: NumberFieldIdeal;
  protected _cachedIsPrime?: boolean;
  protected _cachedTwoGenerators?: [NumberFieldElement, NumberFieldElement];
  protected _cachedFreeModule?: { basis: NumberFieldElement[]; rank: number };

  constructor(
    number_field: NumberField,
    gens: (NumberFieldIdealGenerator | NumberFieldIdealGenerator[])[]
  ) {
    this._number_field = number_field;
    const values =
      gens.length === 1 && Array.isArray(gens[0]) ? gens[0] : (gens as NumberFieldIdealGenerator[]);
    if (values.length === 0) {
      throw new ValueError(
        'gens must have length at least 1 (zero ideal is not a fractional ideal)'
      );
    }

    this._gens = values.map((g) => {
      if (g instanceof NumberFieldIdeal)
        throw new TypeError(`unable to convert ${g} to ${number_field}`);
      return number_field.__call__(g);
    });
  }

  /**
   * Return the number field this ideal belongs to.
   * @see Reference: sage/rings/number_field/number_field_ideal.py:number_field
   */
  number_field(): NumberField {
    return this._number_field;
  }

  /**
   * Return the generators of this ideal.
   * @see Reference: sage/rings/number_field/number_field_ideal.py:gens
   */
  gens(): NumberFieldElement[] {
    return [...this._gens];
  }

  /**
   * Return the number of generators.
   * @see Reference: sage/rings/number_field/number_field_ideal.py:ngens
   */
  ngens(): number {
    return this._gens.length;
  }

  /**
   * Return the i-th generator.
   */
  gen(i: number): NumberFieldElement {
    if (i < 0 || i >= this._gens.length) {
      throw new ValueError(`generator index ${i} out of range`);
    }
    return this._gens[i]!;
  }

  /**
   * Return the norm of this ideal as a rational number.
   *
   * The norm of an ideal I is [O_K : I], the index of I in the ring of integers.
   * For a principal ideal (a), N(I) = |N(a)|.
   *
   * @see Reference: sage/rings/number_field/number_field_ideal.py:norm
   */
  norm(): Rational {
    if (this._cachedNorm !== undefined) {
      return this._cachedNorm;
    }
    if (this.is_zero()) {
      this._cachedNorm = Rational.zero();
      return this._cachedNorm;
    }
    const n = BigInt(this._number_field.degree());
    const hnf = this._computeHNF();
    let det = 1n;
    for (let i = 0; i < hnf.entries.length; i++) {
      det *= hnf.entries[i]![i]!;
    }
    if (det < 0n) det = -det;
    this._cachedNorm = new Rational(det, hnf.denominator ** n);
    return this._cachedNorm;
  }

  /**
   * Return the absolute norm.
   * @see Reference: sage/rings/number_field/number_field_ideal.py:absolute_norm
   */
  absolute_norm(): Rational {
    return this.norm();
  }

  /**
   * Return the relative norm.
   * @see Reference: sage/rings/number_field/number_field_ideal.py:relative_norm
   */
  relative_norm(): Rational {
    // For absolute number fields, relative norm equals absolute norm
    return this.norm();
  }

  /**
   * Return two generators; the first generates the intersection with QQ.
   * Both generators are elements of the ambient number field, including when
   * the first is fractional. The result is cached as in Sage.
   *
   * @see Reference: sage/rings/number_field/number_field_ideal.py:gens_two
   * @see Deviation: Number-field ideal backend adapters
   */
  gens_two(): [NumberFieldElement, NumberFieldElement] {
    if (this._cachedTwoGenerators) return this._cachedTwoGenerators;
    const K = this._number_field;
    if (this.is_zero()) return (this._cachedTwoGenerators = [K.zero(), K.zero()]);
    const H = this._computeHNF(),
      value = H.entries[0]![0]!;
    const first = K.__call__(new Rational(value, H.denominator));
    // Sage returns zero as the second generator exactly for rational ideals.
    if (H.entries.every((column, j) => column.every((v, i) => v === (i === j ? value : 0n))))
      return (this._cachedTwoGenerators = [first, K.zero()]);
    const [, coordinates] = mat_ideal_two_elt(K._pari_ideal_data().multiplication, H.entries);
    const basis = K._pari_integral_basis();
    let second = K.zero();
    for (let i = 0; i < coordinates.length; i++)
      second = second.add(basis[i]!.mul(new Rational(coordinates[i]!, H.denominator)));
    return (this._cachedTwoGenerators = [first, second]);
  }

  /**
   * Return the smallest positive integer in this ideal.
   *
   * This is the smallest n > 0 such that n is in I.
   *
   * @see Reference: sage/rings/number_field/number_field_ideal.py:smallest_integer
   */
  smallest_integer(): bigint {
    if (this.is_zero()) {
      return 0n;
    }
    // Sage: ZZ(self.pari_hnf()[0,0].numerator())
    return this._intersectionWithQ().numerator;
  }

  /**
   * Check if this is a prime ideal.
   *
   * An ideal P is prime if for all a,b in O_K, ab in P implies a in P or b in P.
   *
   * @see Reference: sage/rings/number_field/number_field_ideal.py:is_prime
   */
  is_prime(): boolean {
    if (this._cachedIsPrime !== undefined) {
      return this._cachedIsPrime;
    }
    this._cachedIsPrime = this._computeIsPrime();
    return this._cachedIsPrime;
  }

  private _computeIsPrime(): boolean {
    if (this.is_zero()) return false;
    const H = this._computeHNF(),
      candidate = idealismaximal(this._number_field._pari_ideal_data(), H.entries, H.denominator);
    if (!candidate || !is_prime(candidate.p)) return false;
    this._cachedPrimeData = candidate;
    return true;
  }

  /**
   * Check if this is a maximal ideal.
   * @see Reference: sage/rings/number_field/number_field_ideal.py:is_maximal
   */
  is_maximal(): boolean {
    return this.is_prime() && !this.is_zero();
  }

  /**
   * Check if this is a principal ideal.
   *
   * An ideal is principal if it can be generated by a single element.
   * For class number 1 fields, all ideals are principal.
   *
   * @see Reference: sage/rings/number_field/number_field_ideal.py:is_principal
   */
  is_principal(): boolean {
    // Trivially principal
    if (this._gens.length === 1) {
      return true;
    }

    // For quadratic fields, check if class number is 1
    const K = this._number_field;
    if (K.degree() === 2) {
      try {
        const h = K.class_number();
        if (h === 1n) {
          return true;
        }
      } catch {
        // Can't compute class number, fall through
      }
    }

    // Check if all generators are scalar multiples of a single element
    if (this._gens.length === 2) {
      const g1 = this._gens[0]!;
      const g2 = this._gens[1]!;

      if (g1.is_zero()) return true;
      if (g2.is_zero()) return true;

      // Try to see if g2 = c * g1 for some c
      try {
        const ratio = g2.div(g1);
        if (ratio.is_integral()) {
          return true;
        }
        const ratioInv = g1.div(g2);
        if (ratioInv.is_integral()) {
          return true;
        }
      } catch {
        // Division failed
      }
    }

    // For general case, would need class group computation
    throw new NotImplementedError(
      'is_principal for general ideals requires class group computation'
    );
  }

  /**
   * Return a generator if this is principal.
   * @see Reference: sage/rings/number_field/number_field_ideal.py:gens_reduced
   */
  gens_reduced(): NumberFieldElement[] {
    if (this._gens.length === 1) {
      return [this._gens[0]!];
    }

    // Try to find a single generator
    if (this._gens.length === 2) {
      const g1 = this._gens[0]!;
      const g2 = this._gens[1]!;

      if (g1.is_zero()) return [g2];
      if (g2.is_zero()) return [g1];

      // Check if one divides the other
      try {
        const ratio = g2.div(g1);
        if (ratio.is_integral()) {
          // g2 = ratio * g1, so (g1, g2) = (g1)
          return [g1];
        }
      } catch {
        // Division failed
      }

      try {
        const ratioInv = g1.div(g2);
        if (ratioInv.is_integral()) {
          return [g2];
        }
      } catch {
        // Division failed
      }
    }

    throw new NotImplementedError('gens_reduced requires LLL reduction for general ideals');
  }

  /**
   * Factorize this ideal into prime ideals.
   *
   * Returns a list of pairs (P, e) where P is a prime ideal and e is its multiplicity.
   * Fractional ideals have negative exponents at denominator primes.
   * Factors follow Sage's lexicographic comparison of integral-basis HNF matrices.
   *
   * @see Reference: sage/rings/number_field/number_field_ideal.py:factor
   * @see Deviation: Number-field ideal class method adapters
   */
  factor(): Array<[NumberFieldIdeal, bigint]> {
    this._requireFractionalMethod('factor');
    if (this._cachedFactorization) return this._cachedFactorization;
    const K = this._number_field,
      H = this._computeHNF(),
      basis = K._pari_integral_basis(),
      factors = idealfactor(K._pari_ideal_data(), H.entries, H.denominator);
    const result = factors.map(([P, e]): [NumberFieldIdeal, bigint] => {
      const generator = P.generator.reduce(
          (a, c, i) => a.add(basis[i]!.scalarMul(new Rational(c))),
          K.zero()
        ),
        I = K.ideal([K.__call__(P.p), generator]);
      I._cachedPrimeData = P;
      I._cachedIsPrime = true;
      return [I, e];
    });
    // Factorization.sort falls back to NumberFieldIdeal._richcmp_: compare the
    // Sage matrices row by row, while our PARI HNFs are stored by columns.
    result.sort(([P], [Q]) => {
      const A = P._computeHNF(),
        B = Q._computeHNF();
      for (let i = 0; i < K.degree(); i++)
        for (let j = 0; j < K.degree(); j++) {
          const a = A.entries[j]![i]! * B.denominator;
          const b = B.entries[j]![i]! * A.denominator;
          if (a !== b) return a < b ? -1 : 1;
        }
      return 0;
    });
    return (this._cachedFactorization = result);
  }

  /**
   * Return the prime lying below (for prime ideals).
   *
   * If this is a prime ideal P lying above p, return p.
   *
   * Port alias for the rational prime in Sage pari_prime(), equal to P intersect ZZ.
   * @see Reference: sage/rings/number_field/number_field_ideal.py:smallest_integer
   * @see Deviation: Number-field ideal valuation adapters
   */
  prime_below(): bigint {
    if (!this.is_prime()) {
      throw new ValueError('ideal is not prime');
    }
    return this.smallest_integer();
  }

  /**
   * Return the ramification index of this prime ideal `P` over the rational
   * prime below it: the exponent of `P` in the factorisation of `p O_K`.
   *
   * @see Reference: sage/rings/number_field/number_field_ideal.py:ramification_index
   * @see Deviation: Number-field ideal class method adapters
   */
  ramification_index(): bigint {
    this._requireFractionalMethod('ramification_index');
    if (!this.is_prime()) {
      throw new ValueError(`${this} is not a prime ideal`);
    }
    return this._cachedPrimeData!.e;
  }

  /**
   * Return the residue class degree `f = [O_K/P : Z/pZ]`.
   *
   * @see Reference: sage/rings/number_field/number_field_ideal.py:residue_class_degree
   * @see Deviation: Number-field ideal class method adapters
   */
  residue_class_degree(): bigint {
    this._requireFractionalMethod('residue_class_degree');
    if (!this.is_prime()) {
      throw new ValueError(`${this} is not a prime ideal`);
    }
    return this._cachedPrimeData!.f;
  }

  /**
   * Return the residue field (for prime ideals).
   *
   * O_K/P is a finite field with p^f elements where f is the residue class degree.
   * Returns an object describing the finite field.
   *
   * @see Reference: sage/rings/number_field/number_field_ideal.py:residue_field
   * @see Deviation: Number-field ideal class method adapters
   */
  residue_field(): { characteristic: bigint; order: bigint; degree: bigint } {
    this._requireFractionalMethod('residue_field');
    if (!this.is_prime()) {
      throw new ValueError('The ideal must be prime');
    }

    const p = this.prime_below();
    const f = this.residue_class_degree();
    const order = p ** f;

    return {
      characteristic: p,
      order,
      degree: f,
    };
  }

  /**
   * Check if this ideal is integral.
   *
   * An ideal is integral if it is contained in the ring of integers.
   *
   * @see Reference: sage/rings/number_field/number_field_ideal.py:is_integral
   */
  is_integral(): boolean {
    if (this.is_zero()) return true;
    // I is integral iff its HNF with respect to the integral basis of O_K has
    // denominator 1.
    return this._computeHNF().denominator === 1n;
  }

  /**
   * Check if this is the zero ideal.
   */
  is_zero(): boolean {
    return this._gens.every((g) => g.is_zero());
  }

  /** Match methods defined only on Sage's NumberFieldFractionalIdeal class.
   * @see Deviation: Number-field ideal class method adapters
   */
  private _requireFractionalMethod(method: string): void {
    if (!(this instanceof NumberFieldFractionalIdeal))
      throw new AttributeError(`'NumberFieldIdeal' object has no attribute '${method}'`);
  }

  /**
   * Return the denominator ideal D in the coprime integral decomposition I = N/D.
   * @see Reference: sage/rings/number_field/number_field_ideal.py:denominator
   * @see Deviation: Number-field ideal class method adapters
   */
  denominator(): NumberFieldIdeal {
    this._requireFractionalMethod('denominator');
    if (!this._cachedDenominator)
      this._cachedDenominator = this.add(this._number_field.ideal(1n)).inverse();
    return this._cachedDenominator;
  }

  /**
   * Return the numerator ideal N in the coprime integral decomposition I = N/D.
   * @see Reference: sage/rings/number_field/number_field_ideal.py:numerator
   * @see Deviation: Number-field ideal class method adapters
   */
  numerator(): NumberFieldIdeal {
    this._requireFractionalMethod('numerator');
    if (!this._cachedNumerator) this._cachedNumerator = this.mul(this.denominator());
    return this._cachedNumerator;
  }

  /**
   * Check if an element is in this ideal.
   *
   * x is in I if x can be written as a linear combination of the generators.
   *
   * @see Reference: sage/rings/number_field/number_field_ideal.py:__contains__
   * @see Deviation: Number-field ideal coercion and centered integral bases
   */
  contains(value: NumberFieldInput | NumberFieldIdeal): boolean {
    let x: NumberFieldElement;
    try {
      x = this._number_field.__call__(value as NumberFieldInput);
    } catch (error) {
      // Ideal_generic.__contains__ catches TypeError only; malformed vector
      // lengths and nonfinite rational coefficients retain their ValueErrors.
      if (error instanceof TypeError) return false;
      throw error;
    }
    if (this.is_zero()) {
      return x.is_zero();
    }
    if (x.is_zero()) {
      return true;
    }
    const K = this._number_field;
    const n = K.degree();
    const hnf = this._computeHNF();
    // Coordinates of denominator * x in the integral basis.
    const basis = K._pari_integral_basis();
    const W: Rational[][] = basis.map((b) => b.list());
    const Winv = ratInverse(W);
    const xs = x.list();
    const v: Rational[] = [];
    for (let k = 0; k < n; k++) {
      let acc = Rational.zero();
      for (let l = 0; l < n; l++) acc = acc.add(xs[l]!.mul(Winv[l]![k]!));
      v.push(acc.mul(new Rational(hnf.denominator)));
    }
    // Solve v = t * H with H lower triangular: back-substitute from the last
    // coordinate downwards.
    const t: Rational[] = new Array(n).fill(Rational.zero());
    const rem = [...v];
    for (let i = n - 1; i >= 0; i--) {
      const d = hnf.entries[i]![i]!;
      const ti = rem[i]!.div(new Rational(d));
      if (ti.denominator !== 1n) return false;
      t[i] = ti;
      if (ti.isZero()) continue;
      for (let j = 0; j <= i; j++) {
        rem[j] = rem[j]!.sub(ti.mul(new Rational(hnf.entries[i]![j]!)));
      }
    }
    return rem.every((c) => c.isZero());
  }

  /**
   * Return the valuation of this ideal at the prime ideal p.
   * Scalar/list arguments are first converted through the field ideal factory.
   * The zero ideal returns 'Infinity'; fractional ideals can have negative values.
   * @see Reference: sage/rings/number_field/number_field_ideal.py:valuation
   * @see Deviation: Number-field ideal valuation adapters
   */
  valuation(value: NumberFieldIdealInput): bigint | 'Infinity' {
    const K = this._number_field,
      p = value instanceof NumberFieldIdeal ? value : K.ideal(value);
    if (p.is_zero()) throw new ValueError(`p (= ${p}) must be nonzero`);
    if (!p.is_prime()) throw new ValueError(`p (= ${p}) must be a prime`);
    if (p.number_field() !== K) {
      // The bundled Python format expression has two placeholders and one operand.
      throw new TypeError('not enough arguments for format string');
    }
    const prime = p._primeData(),
      H = this.is_zero() ? { entries: [], denominator: 1n } : this._computeHNF();
    return idealval(K._pari_ideal_data(), H.entries, H.denominator, prime);
  }

  /** Native prime finalization in the current maximal-order basis.
   * @see Deviation: Number-field ideal valuation adapters
   */
  private _primeData(): NfPrimeIdeal {
    if (!this.is_prime()) throw new ValueError(`${this} is not a prime ideal`);
    return this._cachedPrimeData!;
  }

  /**
   * Return an integral basis for this ideal.
   *
   * Returns a list of elements that form a Z-basis for this ideal.
   *
   * @see Reference: sage/rings/number_field/number_field_ideal.py:integral_basis
   * @see Deviation: Ideal basis and generator representations
   */
  integral_basis(): NumberFieldElement[] {
    return this.zk_basis();
  }

  /**
   * Return a free module representation.
   *
   * Returns an object describing the ideal as a Z-module.
   * The HNF basis spans the ideal in the maximal order, including fractional
   * and nonprincipal ideals. The zero module has empty basis and rank zero.
   * @see Deviation: Ideal basis and generator representations
   *
   * @see Reference: sage/rings/number_field/number_field_ideal.py:free_module
   */
  free_module(): { basis: NumberFieldElement[]; rank: number } {
    if (this._cachedFreeModule) return this._cachedFreeModule;
    const basis = this.integral_basis();
    return (this._cachedFreeModule = { basis, rank: basis.length });
  }

  /**
   * Return the inverse of this fractional ideal.
   *
   * I^{-1} = {x in K : x*I ⊆ O_K}
   *
   * Delegate HNF inversion to PARI's trace-dual ideal arithmetic.
   *
   * @see Reference: sage/rings/number_field/number_field_ideal.py:__invert__
   * @see Deviation: Number-field ideal backend adapters
   * @see Deviation: Number-field ideal class method adapters
   */
  inverse(): NumberFieldIdeal {
    if (!(this instanceof NumberFieldFractionalIdeal))
      throw new TypeError("bad operand type for unary ~: 'NumberFieldIdeal'");
    const H = this._computeHNF();
    const [entries, denominator] = idealHNF_inv(
      this._number_field._pari_ideal_data(),
      H.entries,
      H.denominator
    );
    return NumberFieldIdeal.fromHNF(this._number_field, { entries, denominator });
  }

  /**
   * Multiply two ideals.
   *
   * I * J is generated by all products ab where a in I, b in J.
   *
   * @see Reference: sage/rings/number_field/number_field_ideal.py:__mul__
   * @see Deviation: Number-field ideal coercion and centered integral bases
   */
  mul(other: NumberFieldIdealInput): NumberFieldIdeal {
    if (!(other instanceof NumberFieldIdeal)) other = this._number_field.ideal(other);

    if (this._gens.length === 1 && other._gens.length === 1) {
      if (this._number_field !== other._number_field)
        throw new TypeError(
          `unsupported operand parent(s) for *: '${this._number_field}' and '${other._number_field}'`
        );
      return this._number_field.ideal(this._gens[0]!.mul(other._gens[0]!));
    }
    // PARI returns an empty matrix for zero; Sage's ideal conversion rejects it.
    if (this.is_zero() || other.is_zero())
      throw new TypeError('[;] has unsupported PARI type t_MAT');
    const A = this._computeHNF(),
      B = other._computeHNF();
    const [entries, denominator] = idealmul(
      this._number_field._pari_ideal_data(),
      A.entries,
      A.denominator,
      B.entries,
      B.denominator
    );
    return NumberFieldIdeal.fromHNF(this._number_field, {
      entries,
      denominator,
    });
  }

  /** Build an ideal from a verified native HNF without expanding its generators. */
  private static fromHNF(K: NumberField, input: HNFMatrix): NumberFieldIdeal {
    let common = input.denominator;
    for (const column of input.entries) for (const x of column) common = intGcd(common, x);
    const entries = input.entries.map((column) => column.map((x) => x / common));
    const denominator = input.denominator / common,
      basis = K._pari_integral_basis();
    const generators = entries.map((column) => {
      let element = K.zero();
      for (let i = 0; i < column.length; i++)
        element = element.add(basis[i]!.mul(new Rational(column[i]!, denominator)));
      return element;
    });
    const result = K.ideal(...generators);
    result._cachedHNF = { entries, denominator };
    return result;
  }

  /**
   * Divide two ideals.
   *
   * I / J = I * J^{-1}
   *
   * @see Reference: sage/rings/number_field/number_field_ideal.py:_div_
   * @see Deviation: Number-field ideal coercion and centered integral bases
   * @see Deviation: Fractional ideal arithmetic and remaining intersection routing
   */
  div(other: NumberFieldIdealInput): NumberFieldIdeal {
    // Division uses the monoid's coercion model, unlike Ideal_generic.__mul__.
    const className =
      this instanceof NumberFieldFractionalIdeal
        ? 'NumberFieldFractionalIdeal'
        : 'NumberFieldIdeal';
    // canonical_coercion admits native numeric zero even without a real map.
    // Python's fallback masks a TypeError for a base zero ideal and native scalar.
    if (
      Array.isArray(other) ||
      (typeof other === 'number' && (other !== 0 || this.is_zero())) ||
      (this.is_zero() && typeof other === 'bigint')
    ) {
      const type = Array.isArray(other) ? 'list' : typeof other === 'number' ? 'float' : 'int';
      throw new TypeError(`unsupported operand type(s) for /: '${className}' and '${type}'`);
    }
    const foreign =
      other instanceof NumberFieldIdeal
        ? other.number_field()
        : other instanceof NumberFieldElement
          ? other.parent()
          : this._number_field;
    if (foreign !== this._number_field) {
      const right =
        other instanceof NumberFieldIdeal ? `Monoid of ideals of ${foreign}` : String(foreign);
      throw new TypeError(
        `unsupported operand parent(s) for /: 'Monoid of ideals of ${this._number_field}' and '${right}'`
      );
    }
    if (!(other instanceof NumberFieldIdeal)) other = this._number_field.ideal(other);
    // The zero ideal is a MonoidElement, not a MultiplicativeGroupElement.
    if (this.is_zero()) return this.mul(other.inverse());
    if (this._gens.length === 1 && other._gens.length === 1)
      return this._number_field.ideal(this._gens[0]!.div(other._gens[0]!));
    const A = this._computeHNF();
    const B = other.is_zero() ? { entries: [], denominator: 1n } : other._computeHNF();
    const [entries, denominator] = idealdiv(
      this._number_field._pari_ideal_data(),
      A.entries,
      A.denominator,
      B.entries,
      B.denominator
    );
    return NumberFieldIdeal.fromHNF(this._number_field, { entries, denominator });
  }

  /**
   * Return this ideal raised to a power.
   *
   * @see Reference: sage/rings/number_field/number_field_ideal.py:__pow__
   * @see Deviation: Number-field ideal class method adapters
   */
  pow(exponent: IntegerLike): NumberFieldIdeal {
    const n = toBigInt(exponent);
    if (!(this instanceof NumberFieldFractionalIdeal) && n < 0n) {
      // Element.__pow__ catches TypeError for Python-int exponents, while the
      // Sage Integer coercion path exposes the failed inverse directly.
      if (exponent instanceof Integer) return this.inverse();
      throw new TypeError(
        "unsupported operand type(s) for ** or pow(): 'NumberFieldIdeal' and 'int'"
      );
    }
    if (n === 0n) return this._number_field.ideal(1n);
    return n < 0n
      ? generic_power_pos(this.inverse(), -n)
      : generic_power_pos<NumberFieldIdeal>(this, n);
  }

  /**
   * Return the sum (GCD) of two ideals.
   *
   * I + J is the smallest ideal containing both I and J.
   *
   * @see Reference: sage/rings/number_field/number_field_ideal.py:__add__
   * @see Deviation: Number-field ideal coercion and centered integral bases
   */
  add(other: NumberFieldIdealInput): NumberFieldIdeal {
    if (!(other instanceof NumberFieldIdeal)) other = this._number_field.ideal(other);

    // Combine generators
    const newGens = [...this._gens, ...other._gens];
    return this._number_field.ideal(...newGens);
  }

  /**
   * Return the intersection (LCM) of two ideals.
   *
   * I ∩ J is the largest ideal contained in both I and J.
   * For principal ideals, (a) ∩ (b) is related to lcm(a, b) in the ring of integers.
   *
   * Delegate the rational HNF intersection to PARI's LLL-kernel route.
   *
   * @see Reference: sage/rings/number_field/number_field_ideal.py:intersection
   * @see Deviation: Number-field ideal intersection and construction adapters
   */
  intersection(other: NumberFieldIdealInput): NumberFieldIdeal {
    const K = this._number_field,
      J = K.ideal(other);
    const hnfOf = (I: NumberFieldIdeal): HNFMatrix =>
      I.is_zero() ? { entries: [], denominator: 1n } : I._computeHNF();
    const A = hnfOf(this),
      B = hnfOf(J);
    const [entries, denominator] = idealintersect(
      K._pari_ideal_data(),
      A.entries,
      A.denominator,
      B.entries,
      B.denominator
    );
    return entries.length ? NumberFieldIdeal.fromHNF(K, { entries, denominator }) : K.ideal(0n);
  }

  /**
   * Check equality of ideals.
   * @see Reference: sage/rings/number_field/number_field_ideal.py:__eq__
   */
  eq(other: NumberFieldIdeal): boolean {
    if (this._number_field !== other._number_field) {
      return false;
    }
    if (this.is_zero() || other.is_zero()) {
      return this.is_zero() && other.is_zero();
    }
    // The HNF with respect to the integral basis is a canonical form for the
    // underlying lattice, so it decides equality outright.
    return hnfKey(this._computeHNF()) === hnfKey(other._computeHNF());
  }

  /**
   * Check if this ideal divides another.
   *
   * I | J iff J ⊆ I
   *
   * @see Reference: sage/rings/number_field/number_field_ideal.py:divides
   * @see Deviation: Number-field ideal coercion and centered integral bases
   * @see Deviation: Number-field ideal class method adapters
   */
  divides(other: NumberFieldIdealInput): boolean {
    this._requireFractionalMethod('divides');
    if (!(other instanceof NumberFieldIdeal)) other = this._number_field.ideal(other);

    // I | J iff J/I is integral
    const quotient = other.div(this);
    return quotient.is_integral();
  }

  /**
   * Check if this ideal is coprime to another.
   *
   * I and J are coprime iff I + J = O_K.
   *
   * @see Reference: sage/rings/number_field/number_field_ideal.py:is_coprime
   * @see Deviation: Number-field ideal coercion and centered integral bases
   * @see Deviation: Number-field ideal class method adapters
   */
  is_coprime(other: NumberFieldIdealInput): boolean {
    this._requireFractionalMethod('is_coprime');
    other = this._number_field.ideal(other);
    const one = this._number_field.ideal(1n);
    if (this.is_integral() && other.is_integral()) {
      if (intGcd(this.norm().numerator, other.norm().numerator) === 1n) return true;
      return this.add(other).eq(one);
    }
    if (other.is_zero()) return this.eq(one);
    const D1 = this.denominator(),
      N1 = this.numerator();
    const D2 = other.denominator(),
      N2 = other.numerator();
    return N1.add(N2).eq(one) && N1.add(D2).eq(one) && D1.add(N2).eq(one) && D1.add(D2).eq(one);
  }

  /**
   * Return the ideal class of this ideal in the class group.
   *
   * Two ideals are in the same class if their quotient is a principal ideal.
   *
   * @see Reference: sage/rings/number_field/number_field_ideal.py:ideal_class
   */
  ideal_class(): unknown {
    // Get the class group
    const classGroup = this._number_field.class_group();

    // If class number is 1, all ideals are principal
    if (classGroup.order() === 1n) {
      return classGroup.identity();
    }

    // For principal ideals, return the identity
    if (this._gens.length === 1) {
      return classGroup.identity();
    }

    // General case requires computing the discrete log in the class group
    throw new NotImplementedError(
      'ideal_class requires PARI bnfisprincipal for non-principal ideals'
    );
  }

  /**
   * Return the ideal class in the narrow class group.
   *
   * The narrow class group uses totally positive generators for principal ideals.
   *
   * @see Reference: sage/rings/number_field/number_field_ideal.py:ideal_class_narrow
   */
  ideal_class_narrow(): unknown {
    // For imaginary quadratic fields, narrow = ordinary class group
    if (this._number_field.degree() === 2) {
      const disc = this._number_field.discriminant();
      if (disc < 0n) {
        return this.ideal_class();
      }
    }

    throw new NotImplementedError('ideal_class_narrow requires PARI for real fields');
  }

  /**
   * A Z-basis of this ideal as a lattice in `O_K`, read off the Hermite normal
   * form.  This is PARI's `idealhnf` shape (base4.c:1015), which SageMath
   * exposes as `I.pari_hnf()`.
   *
   * @see Reference: sage/rings/number_field/number_field_ideal.py:pari_hnf
   */
  zk_basis(): NumberFieldElement[] {
    if (this.is_zero()) return [];
    const K = this._number_field;
    const n = K.degree();
    const hnf = this._computeHNF();
    const zk = K._pari_integral_basis();
    const den = new Rational(hnf.denominator);
    const out: NumberFieldElement[] = [];
    for (let i = 0; i < n; i++) {
      let acc = K.zero();
      for (let j = 0; j < n; j++) {
        const c = hnf.entries[i]![j]!;
        if (c === 0n) continue;
        acc = acc.add(zk[j]!.scalarMul(new Rational(c).div(den)));
      }
      out.push(acc);
    }
    return out;
  }

  /**
   * Compute the Hermite normal form of this ideal with respect to the integral
   * basis of the field.
   *
   * SageMath obtains this from PARI (`nf.idealhnf`); here the Z-module
   * generated by `{g_i w_j}` is put in lower-triangular HNF, which is PARI's
   * shape and makes `entries[0][0]/denominator` the generator of `I \cap Q`.
   *
   * @see Reference: sage/rings/number_field/number_field_ideal.py:pari_hnf
   */
  protected _computeHNF(): HNFMatrix {
    if (this._cachedHNF) {
      return this._cachedHNF;
    }

    const K = this._number_field;
    const n = K.degree();
    const basis = K._pari_integral_basis();
    // W[i][j] = coefficient of alpha^j in w_i
    const W: Rational[][] = basis.map((b) => b.list());
    const Winv = ratInverse(W);

    // Sage builds each principal HNF before combining it with idealadd. The
    // modular bound from zkmultable_capZ prevents unbounded intermediate growth.
    const nf = K._pari_ideal_data();
    let entries: bigint[][] = [],
      denom = 1n;
    for (const generator of this._gens) {
      const power = generator.list();
      const coordinates = Array.from({ length: n }, (_, i) =>
        power.reduce((sum, c, j) => sum.add(c.mul(Winv[j]![i]!)), Rational.zero())
      );
      const d = coordinates.reduce((den, c) => intLcm(den, c.denominator), 1n);
      const [H, hDen] = idealhnf_principal(
        nf,
        coordinates.map((c) => c.numerator * (d / c.denominator)),
        d
      );
      [entries, denom] = idealadd(nf, entries, denom, H, hDen);
    }
    if (!entries.length) throw new ValueError('the zero ideal has no Hermite normal form');

    this._cachedHNF = { entries, denominator: denom };
    return this._cachedHNF;
  }

  /** The rational number `q` with `I \cap Q = q Z`. */
  private _intersectionWithQ(): Rational {
    const hnf = this._computeHNF();
    return new Rational(hnf.entries[0]![0]!, hnf.denominator);
  }

  protected _checkSameField(other: NumberFieldIdeal): void {
    if (this._number_field !== other._number_field) {
      throw new ValueError('ideals must be in the same number field');
    }
  }

  toString(): string {
    const generators = `(${this._gens.join(', ')})`;
    return this instanceof NumberFieldFractionalIdeal
      ? `Fractional ideal ${generators}`
      : `Ideal ${generators} of ${this._number_field}`;
  }
}

/**
 * A fractional ideal of a number field.
 * @see Reference: sage/rings/number_field/number_field_ideal.py:NumberFieldFractionalIdeal
 * @see Deviation: Number-field ideal intersection and construction adapters
 */
export class NumberFieldFractionalIdeal extends NumberFieldIdeal {
  constructor(
    number_field: NumberField,
    gens: (NumberFieldIdealGenerator | NumberFieldIdealGenerator[])[]
  ) {
    if (!gens.length)
      throw new ValueError(
        'gens must have length at least 1 (zero ideal is not a fractional ideal)'
      );
    const values =
      gens.length === 1 && Array.isArray(gens[0]) ? gens[0] : (gens as NumberFieldIdealGenerator[]);
    if (values.every(idealGeneratorIsZero))
      throw new ValueError(
        'gens must have a nonzero element (zero ideal is not a fractional ideal)'
      );
    super(number_field, values);
  }
}

/** Native truth test before field coercion in fractional-ideal construction. */
function idealGeneratorIsZero(value: NumberFieldIdealGenerator | NumberFieldIdeal): boolean {
  if (typeof value === 'bigint') return value === 0n;
  if (typeof value === 'number') return value === 0;
  if (value instanceof Integer) return value.value === 0n;
  if (value instanceof Rational) return value.numerator === 0n;
  return value.is_zero();
}

// Helper functions

/** Canonical string for an ideal HNF, used to decide equality of ideals. */
function hnfKey(h: HNFMatrix): string {
  return `${h.denominator}|${h.entries.map((r) => r.join(',')).join(';')}`;
}

/**
 * Find the smallest prime factor of n.
 */
