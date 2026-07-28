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

import { gcd as intGcd, lcm as intLcm, is_prime_power } from '../../arith/misc.js';
import { NotImplementedError, ValueError, ZeroDivisionError } from '../../errors.js';
import { Rational } from '../rational.js';
import type { NumberField } from './number_field.js';
import { NumberFieldElement } from './number_field_element.js';
import { hnfLower, ratInverse } from './pari_nf.js';

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
 */
export class NumberFieldIdeal {
  protected readonly _number_field: NumberField;
  protected readonly _gens: NumberFieldElement[];
  protected _cachedHNF?: HNFMatrix;
  protected _cachedNorm?: Rational;
  protected _cachedIsPrime?: boolean;

  constructor(number_field: NumberField, gens: NumberFieldElement[]) {
    this._number_field = number_field;

    if (gens.length === 0) {
      throw new ValueError(
        'gens must have length at least 1 (zero ideal is not a fractional ideal)'
      );
    }

    this._gens = gens;
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
   * Return a two-element representation (p, alpha).
   *
   * Every ideal in a number field can be written as (p, alpha) where
   * p is a rational integer and alpha is an algebraic integer.
   *
   * @see Reference: sage/rings/number_field/number_field_ideal.py:gens_two
   */
  gens_two(): [bigint, NumberFieldElement] {
    // For principal ideals
    if (this._gens.length === 1) {
      const gen = this._gens[0]!;
      return [this.smallest_integer(), gen];
    }

    const smallestInt = this.smallest_integer();
    if (smallestInt === 0n) {
      return [0n, this._number_field.zero()];
    }

    const isRational = (coeffs: Rational[]) => coeffs.slice(1).every((c) => c.isZero());
    for (const g of this._gens) {
      if (!isRational(g.list())) {
        return [smallestInt, g];
      }
    }
    return [smallestInt, this._number_field.__call__(smallestInt)];
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
    if (this.is_zero()) {
      // The zero ideal is prime in an integral domain, but Sage's is_prime
      // consults idealismaximal, which rejects it.
      return false;
    }
    if (!this.is_integral()) {
      return false;
    }
    const norm = this.norm();
    if (norm.denominator !== 1n) return false;
    const N = norm.numerator;
    if (N === 1n) return false;
    const data = is_prime_power(N, true);
    if (data[1] === 0n) return false; // norm is not a prime power
    const p = data[0];
    // Compare against the actual prime decomposition of p (Dedekind-Kummer).
    const decomposition = this._number_field.decomposition(p);
    const key = hnfKey(this._computeHNF());
    for (const [P] of decomposition) {
      if (hnfKey(P._computeHNF()) === key) return true;
    }
    return false;
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
   * For principal ideals generated by a prime power, returns the factorization.
   *
   * @see Reference: sage/rings/number_field/number_field_ideal.py:factor
   */
  factor(): Array<[NumberFieldIdeal, bigint]> {
    if (this.is_zero()) {
      return [];
    }

    // For principal ideals, factor the generator's norm
    if (this._gens.length === 1) {
      const gen = this._gens[0]!;
      const norm = gen.norm();
      const normInt = norm.numerator / norm.denominator;
      const normAbs = normInt < 0n ? -normInt : normInt;

      if (normAbs === 1n) {
        // Unit, so ideal is (1)
        return [];
      }

      // Factor the norm to find the primes involved
      // Then for each prime p, factor pO_K
      throw new NotImplementedError('factor requires prime decomposition via PARI');
    }

    throw new NotImplementedError('factor for non-principal ideals requires PARI idealfactor');
  }

  /**
   * Return the prime lying below (for prime ideals).
   *
   * If this is a prime ideal P lying above p, return p.
   *
   * @see Reference: sage/rings/number_field/number_field_ideal.py:prime_below
   */
  prime_below(): bigint {
    if (!this.is_prime()) {
      throw new ValueError('ideal is not prime');
    }
    return smallestPrimeFactor(this.norm().numerator);
  }

  /**
   * Return the ramification index of this prime ideal `P` over the rational
   * prime below it: the exponent of `P` in the factorisation of `p O_K`.
   *
   * @see Reference: sage/rings/number_field/number_field_ideal.py:ramification_index
   */
  ramification_index(): bigint {
    if (!this.is_prime()) {
      throw new ValueError('ramification index only defined for prime ideals');
    }
    const p = this.prime_below();
    const key = hnfKey(this._computeHNF());
    for (const [P, e] of this._number_field.decomposition(p)) {
      if (hnfKey(P._computeHNF()) === key) return e;
    }
    throw new ValueError('prime ideal not found in the decomposition of the prime below');
  }

  /**
   * Return the residue class degree `f = [O_K/P : Z/pZ]`.
   *
   * @see Reference: sage/rings/number_field/number_field_ideal.py:residue_class_degree
   */
  residue_class_degree(): bigint {
    if (!this.is_prime()) {
      throw new ValueError('residue class degree only defined for prime ideals');
    }
    const p = this.prime_below();
    let f = 0n;
    let temp = this.norm().numerator;
    while (temp % p === 0n) {
      temp /= p;
      f++;
    }
    if (temp !== 1n) {
      throw new ValueError('norm is not a power of the prime below');
    }
    return f;
  }

  /**
   * Return the residue field (for prime ideals).
   *
   * O_K/P is a finite field with p^f elements where f is the residue class degree.
   * Returns an object describing the finite field.
   *
   * @see Reference: sage/rings/number_field/number_field_ideal.py:residue_field
   */
  residue_field(): { characteristic: bigint; order: bigint; degree: bigint } {
    if (!this.is_prime()) {
      throw new ValueError('residue_field only defined for prime ideals');
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

  /**
   * Return the denominator.
   *
   * The denominator is the smallest positive integer d such that d*I is integral.
   *
   * @see Reference: sage/rings/number_field/number_field_ideal.py:denominator
   */
  denominator(): bigint {
    let denom = 1n;

    for (const g of this._gens) {
      denom = intLcm(denom, g.denominator());
    }

    return denom;
  }

  /**
   * Return the numerator.
   *
   * If I = J/d where J is integral and d is the denominator, return J.
   *
   * @see Reference: sage/rings/number_field/number_field_ideal.py:numerator
   */
  numerator(): NumberFieldIdeal {
    const d = this.denominator();
    if (d === 1n) {
      return this;
    }

    // Multiply all generators by d
    const newGens = this._gens.map((g) => g.scalarMul(new Rational(d)));
    return new NumberFieldIdeal(this._number_field, newGens);
  }

  /**
   * Check if an element is in this ideal.
   *
   * x is in I if x can be written as a linear combination of the generators.
   *
   * @see Reference: sage/rings/number_field/number_field_ideal.py:__contains__
   */
  contains(x: NumberFieldElement): boolean {
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
   * Return the valuation of x at this ideal.
   *
   * v_P(x) is the largest n such that x is in P^n.
   *
   * @see Reference: sage/rings/number_field/number_field_ideal.py:valuation
   */
  valuation(x: NumberFieldElement): bigint {
    if (!this.is_prime()) {
      throw new ValueError('valuation only defined at prime ideals');
    }

    if (x.is_zero()) {
      // Infinity for zero
      throw new ValueError('valuation of zero is infinity');
    }

    // For principal prime ideals generated by pi, v_P(x) is the largest k
    // such that pi^k | x
    if (this._gens.length === 1) {
      const pi = this._gens[0]!;
      let v = 0n;
      let current = x;

      // While pi divides current, increment v
      while (true) {
        try {
          const quotient = current.div(pi);
          if (!quotient.is_integral()) break;
          current = quotient;
          v++;

          // Safety limit
          if (v > 1000n) {
            throw new Error('valuation computation exceeded limit');
          }
        } catch {
          break;
        }
      }

      return v;
    }

    // For non-principal ideals, use the norm
    // v_P(x) can be computed from the factorization of the principal ideal (x)
    // N(P)^{v_P(x)} divides N(x)

    const norm = x.norm();
    const normAbs = norm.numerator < 0n ? -norm.numerator : norm.numerator;
    const p = this.prime_below();
    const f = this.residue_class_degree();

    // v_P(x) * f = v_p(N(x)) for degree 1 primes
    // More generally, this gives a lower bound

    let v = 0n;
    let remaining = normAbs / norm.denominator;

    while (remaining % p === 0n) {
      remaining /= p;
      v++;
    }

    // v is now v_p(N(x)), and v_P(x) >= v/f
    // For exact computation, need PARI
    if (f === 1n) {
      return v;
    }

    throw new NotImplementedError('valuation for non-degree-1 primes requires PARI');
  }

  /**
   * Return an integral basis for this ideal.
   *
   * Returns a list of elements that form a Z-basis for this ideal.
   *
   * @see Reference: sage/rings/number_field/number_field_ideal.py:integral_basis
   */
  integral_basis(): NumberFieldElement[] {
    // For a principal ideal (a), basis is {a * b_i} where {b_i} is basis of O_K
    if (this._gens.length === 1) {
      const a = this._gens[0]!;
      const fieldBasis = this._number_field.power_basis();
      return fieldBasis.map((b) => a.mul(b));
    }

    // For two-generator ideals (p, alpha) in a quadratic field,
    // we can compute the basis directly
    if (this._gens.length === 2 && this._number_field.degree() === 2) {
      const g1 = this._gens[0]!;
      const g2 = this._gens[1]!;

      // Check if g1 is a rational integer
      const g1Coeffs = g1.list();
      const isG1Rational = g1Coeffs.slice(1).every((c) => c.isZero());

      if (isG1Rational) {
        const p = g1Coeffs[0]!.numerator / g1Coeffs[0]!.denominator;
        // For a prime ideal (p, alpha), a basis is {p, alpha} or {alpha, p*omega}
        // where omega is the second basis element of O_K
        const fieldBasis = this._number_field.power_basis();
        return [g1.mul(fieldBasis[0]!), g2];
      }
    }

    throw new NotImplementedError(
      'integral_basis for general non-principal ideals requires HNF computation'
    );
  }

  /**
   * Return a free module representation.
   *
   * Returns an object describing the ideal as a Z-module.
   * For a principal ideal (a), returns the coordinates of the basis {a, a*alpha, ...}
   * in terms of the field's integral basis.
   *
   * @see Reference: sage/rings/number_field/number_field_ideal.py:free_module
   */
  free_module(): { basis: NumberFieldElement[]; rank: number } {
    const basis = this.integral_basis();
    return {
      basis,
      rank: basis.length,
    };
  }

  /**
   * Return the inverse of this fractional ideal.
   *
   * I^{-1} = {x in K : x*I ⊆ O_K}
   *
   * For a two-element ideal (a, b), we use the formula:
   * I^{-1} = (1/N(I)) * conjugate(I) for prime ideals in quadratic fields.
   *
   * @see Reference: sage/rings/number_field/number_field_ideal.py:inverse
   */
  inverse(): NumberFieldIdeal {
    if (this.is_zero()) {
      throw new ZeroDivisionError('cannot invert zero ideal');
    }

    // For a principal ideal (a), the inverse is (1/a)
    if (this._gens.length === 1) {
      const aInv = this._gens[0]!.inv();
      return new NumberFieldIdeal(this._number_field, [aInv]);
    }

    // For two-generator ideals in quadratic fields, we can compute the inverse
    // using the formula I^(-1) = conjugate(I) / N(I)
    if (this._gens.length === 2 && this._number_field.degree() === 2) {
      const normInv = this.norm().inv();

      // Conjugate each generator (for quadratic fields, conjugation negates the sqrt(d) part)
      const conjGens = this._gens.map((g) => {
        const coeffs = g.list();
        if (coeffs.length >= 2) {
          return new NumberFieldElement(this._number_field, [coeffs[0]!, coeffs[1]!.neg()]);
        }
        return g;
      });

      // Scale by 1/N(I)
      const invGens = conjGens.map((g) => g.scalarMul(normInv));
      return new NumberFieldIdeal(this._number_field, invGens);
    }

    throw new NotImplementedError(
      'inverse for general non-principal ideals requires PARI idealinv'
    );
  }

  /**
   * Multiply two ideals.
   *
   * I * J is generated by all products ab where a in I, b in J.
   *
   * @see Reference: sage/rings/number_field/number_field_ideal.py:__mul__
   */
  mul(other: NumberFieldIdeal): NumberFieldIdeal {
    this._checkSameField(other);

    // Generate products of generators
    const newGens: NumberFieldElement[] = [];

    for (const a of this._gens) {
      for (const b of other._gens) {
        newGens.push(a.mul(b));
      }
    }

    return new NumberFieldIdeal(this._number_field, newGens);
  }

  /**
   * Divide two ideals.
   *
   * I / J = I * J^{-1}
   *
   * @see Reference: sage/rings/number_field/number_field_ideal.py:__truediv__
   */
  div(other: NumberFieldIdeal): NumberFieldIdeal {
    return this.mul(other.inverse());
  }

  /**
   * Return this ideal raised to a power.
   *
   * @see Reference: sage/rings/number_field/number_field_ideal.py:__pow__
   */
  pow(n: bigint): NumberFieldIdeal {
    if (n === 0n) {
      // I^0 = O_K = (1)
      return new NumberFieldIdeal(this._number_field, [this._number_field.one()]);
    }

    if (n < 0n) {
      return this.inverse().pow(-n);
    }

    // Binary exponentiation
    let result = new NumberFieldIdeal(this._number_field, [this._number_field.one()]);
    let base: NumberFieldIdeal = this;

    while (n > 0n) {
      if (n % 2n === 1n) {
        result = result.mul(base);
      }
      base = base.mul(base);
      n = n / 2n;
    }

    return result;
  }

  /**
   * Return the sum (GCD) of two ideals.
   *
   * I + J is the smallest ideal containing both I and J.
   *
   * @see Reference: sage/rings/number_field/number_field_ideal.py:__add__
   */
  add(other: NumberFieldIdeal): NumberFieldIdeal {
    this._checkSameField(other);

    // Combine generators
    const newGens = [...this._gens, ...other._gens];
    return new NumberFieldIdeal(this._number_field, newGens);
  }

  /**
   * Return the intersection (LCM) of two ideals.
   *
   * I ∩ J is the largest ideal contained in both I and J.
   * For principal ideals, (a) ∩ (b) is related to lcm(a, b) in the ring of integers.
   *
   * We use the formula: I ∩ J = I * J / (I + J)
   * where I + J is the GCD (sum) of ideals.
   *
   * @see Reference: sage/rings/number_field/number_field_ideal.py:intersection
   */
  intersection(other: NumberFieldIdeal): NumberFieldIdeal {
    this._checkSameField(other);

    // For principal ideals in number fields where both generators are rational integers
    if (this._gens.length === 1 && other._gens.length === 1) {
      const a = this._gens[0]!;
      const b = other._gens[0]!;

      // Check if both are rational integers
      const aCoeffs = a.list();
      const bCoeffs = b.list();
      const isARational = aCoeffs.slice(1).every((c) => c.isZero());
      const isBRational = bCoeffs.slice(1).every((c) => c.isZero());

      if (isARational && isBRational) {
        const aInt = aCoeffs[0]!.numerator / aCoeffs[0]!.denominator;
        const bInt = bCoeffs[0]!.numerator / bCoeffs[0]!.denominator;
        const lcmVal = intLcm(aInt, bInt);
        return new NumberFieldIdeal(this._number_field, [this._number_field.__call__(lcmVal)]);
      }
    }

    // Use formula: I ∩ J = (I * J) / (I + J)
    // This works because I + J is the GCD
    const product = this.mul(other);
    const sum = this.add(other);

    try {
      return product.div(sum);
    } catch {
      throw new NotImplementedError('intersection requires HNF computation for this case');
    }
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
   */
  divides(other: NumberFieldIdeal): boolean {
    this._checkSameField(other);

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
   */
  is_coprime(other: NumberFieldIdeal): boolean {
    this._checkSameField(other);

    // Coprime iff I + J = O_K.
    const sum = this.add(other);
    return sum.norm().eq(Rational.one());
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

    // Coordinates in the integral basis of every product g_i * w_j.
    const coords: Rational[][] = [];
    for (const g of this._gens) {
      if (g.is_zero()) continue;
      for (const w of basis) {
        const prod = g.mul(w).list();
        const row: Rational[] = [];
        for (let k = 0; k < n; k++) {
          let acc = Rational.zero();
          for (let l = 0; l < n; l++) {
            acc = acc.add(prod[l]!.mul(Winv[l]![k]!));
          }
          row.push(acc);
        }
        coords.push(row);
      }
    }

    if (coords.length === 0) {
      throw new ValueError('the zero ideal has no Hermite normal form');
    }

    let denom = 1n;
    for (const row of coords) {
      for (const c of row) denom = intLcm(denom, c.denominator);
    }
    const rows = coords.map((row) => row.map((c) => c.numerator * (denom / c.denominator)));
    const entries = hnfLower(rows, n);
    // Reduce by the common content.
    let g = denom;
    for (const row of entries) {
      for (const x of row) g = intGcd(g, x);
    }
    if (g > 1n) {
      denom /= g;
      for (const row of entries) {
        for (let j = 0; j < n; j++) row[j] = row[j]! / g;
      }
    }

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
    if (this._gens.length === 1) {
      return `Fractional ideal (${this._gens[0]})`;
    }
    return `Fractional ideal (${this._gens.join(', ')})`;
  }
}

/**
 * A fractional ideal of a number field.
 * @see Reference: sage/rings/number_field/number_field_ideal.py:NumberFieldFractionalIdeal
 */
export class NumberFieldFractionalIdeal extends NumberFieldIdeal {
  // Fractional ideals have the same interface but allow denominators
}

// Helper functions

/** Canonical string for an ideal HNF, used to decide equality of ideals. */
function hnfKey(h: HNFMatrix): string {
  return `${h.denominator}|${h.entries.map((r) => r.join(',')).join(';')}`;
}

/**
 * Find the smallest prime factor of n.
 */
function smallestPrimeFactor(n: bigint): bigint {
  if (n <= 1n) return n;

  if (n % 2n === 0n) return 2n;

  let i = 3n;
  while (i * i <= n) {
    if (n % i === 0n) return i;
    i += 2n;
  }

  return n;
}
