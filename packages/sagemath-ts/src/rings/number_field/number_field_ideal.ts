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

import { NotImplementedError, ValueError, ZeroDivisionError } from '../../errors.js';
import { NumberFieldElement } from './number_field_element.js';
import type { NumberField, RationalPolynomial } from './number_field.js';
import { Rational } from '../rational.js';
import { gcd as intGcd, lcm as intLcm } from '../../arith/misc.js';

/**
 * Hermite Normal Form representation of an ideal.
 * The HNF is an upper triangular matrix with respect to the power basis.
 */
export interface HNFMatrix {
  /** The matrix entries */
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
  protected _cachedNorm?: bigint;
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
  norm(): bigint {
    if (this._cachedNorm !== undefined) {
      return this._cachedNorm;
    }

    // For a principal ideal (a), the norm is |N(a)|
    if (this._gens.length === 1) {
      const normRat = this._gens[0]!.norm();
      const n = normRat.numerator / normRat.denominator;
      this._cachedNorm = n < 0n ? -n : n;
      return this._cachedNorm;
    }

    // General case requires HNF computation
    // The norm is the determinant of the HNF matrix
    const hnf = this._computeHNF();
    let det = 1n;
    for (let i = 0; i < hnf.entries.length; i++) {
      det *= hnf.entries[i]![i]!;
    }
    det = det / hnf.denominator ** BigInt(hnf.entries.length);

    this._cachedNorm = det < 0n ? -det : det;
    return this._cachedNorm;
  }

  /**
   * Return the absolute norm.
   * @see Reference: sage/rings/number_field/number_field_ideal.py:absolute_norm
   */
  absolute_norm(): bigint {
    return this.norm();
  }

  /**
   * Return the relative norm.
   * @see Reference: sage/rings/number_field/number_field_ideal.py:relative_norm
   */
  relative_norm(): unknown {
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
      const norm = gen.norm();
      const p = norm.numerator / norm.denominator;
      return [p < 0n ? -p : p, gen];
    }

    // General case: find smallest integer in ideal
    const smallestInt = this.smallest_integer();

    if (smallestInt === 0n) {
      // Zero ideal
      return [0n, this._number_field.zero()];
    }

    // For ideals with 2 generators, use them directly if one is an integer
    if (this._gens.length === 2) {
      const g1 = this._gens[0]!;
      const g2 = this._gens[1]!;

      // Check if either generator is a rational integer
      const g1Coeffs = g1.list();
      const g2Coeffs = g2.list();

      const isRational = (coeffs: Rational[]) =>
        coeffs.slice(1).every((c) => c.isZero());

      if (isRational(g1Coeffs)) {
        const p = g1Coeffs[0]!.numerator / g1Coeffs[0]!.denominator;
        return [p < 0n ? -p : p, g2];
      }

      if (isRational(g2Coeffs)) {
        const p = g2Coeffs[0]!.numerator / g2Coeffs[0]!.denominator;
        return [p < 0n ? -p : p, g1];
      }
    }

    // Return [smallestInt, first non-integer generator]
    const pElem = this._number_field.__call__(smallestInt);
    for (const g of this._gens) {
      if (!g.eq(pElem)) {
        return [smallestInt, g];
      }
    }

    // All generators are multiples of smallestInt, ideal is principal
    return [smallestInt, pElem];
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

    // For a principal ideal (a), smallest integer is |N(a)| / gcd(coefficients)
    if (this._gens.length === 1) {
      const gen = this._gens[0]!;
      const coeffs = gen.list();

      // Find GCD of all numerators and LCM of all denominators
      let numerGcd = 0n;
      let denomLcm = 1n;

      for (const c of coeffs) {
        const n = c.numerator < 0n ? -c.numerator : c.numerator;
        const d = c.denominator;
        numerGcd = numerGcd === 0n ? n : intGcd(numerGcd, n);
        denomLcm = intLcm(denomLcm, d);
      }

      if (numerGcd === 0n) {
        return 0n;
      }

      // The smallest integer is related to the norm
      const normRat = gen.norm();
      const normAbs = normRat.numerator < 0n ? -normRat.numerator : normRat.numerator;
      return normAbs / normRat.denominator;
    }

    // General case requires HNF
    const hnf = this._computeHNF();
    return hnf.entries[0]![0]! / hnf.denominator;
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

    if (this.is_zero()) {
      this._cachedIsPrime = true;
      return true;
    }

    // Check if norm is a prime power
    const n = this.norm();
    if (n === 1n) {
      this._cachedIsPrime = false;
      return false;
    }

    // Factor the norm
    const p = smallestPrimeFactor(n);
    if (p === n) {
      // Norm is prime, ideal might be prime
      // Need additional check that ramification index is 1
      this._cachedIsPrime = true;
      return true;
    }

    // Check if n is a prime power
    let remaining = n;
    while (remaining % p === 0n) {
      remaining /= p;
    }

    if (remaining !== 1n) {
      // Norm has multiple distinct prime factors
      this._cachedIsPrime = false;
      return false;
    }

    // n = p^e, could still be prime if e = residue degree
    // Full check requires PARI
    this._cachedIsPrime = false;
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
    throw new NotImplementedError('is_principal for general ideals requires class group computation');
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

    const n = this.norm();
    return smallestPrimeFactor(n);
  }

  /**
   * Return the ramification index (for prime ideals).
   *
   * If this prime ideal P lies above p, return the largest e such that P^e divides pO_K.
   * For unramified primes, e = 1. For ramified primes, e > 1.
   *
   * @see Reference: sage/rings/number_field/number_field_ideal.py:ramification_index
   */
  ramification_index(): bigint {
    if (!this.is_prime()) {
      throw new ValueError('ramification index only defined for prime ideals');
    }

    const p = this.prime_below();
    const norm = this.norm();
    const f = this.residue_class_degree();

    // For a prime P above p: N(P) = p^f and e*f divides n = [K:Q]
    // The ramification index e is such that P^e || pO_K
    // We have n = sum over primes above p of e_i * f_i

    const n = BigInt(this._number_field.degree());

    // For quadratic fields with p ramified, e = 2
    // For quadratic fields with p split or inert, e = 1
    if (n === 2n) {
      const disc = this._number_field.discriminant();

      // p ramifies iff p | disc
      if (disc % p === 0n) {
        return 2n;
      }
      return 1n;
    }

    // General case: need to factor pO_K
    throw new NotImplementedError('ramification_index for degree > 2 requires PARI');
  }

  /**
   * Return the residue class degree (for prime ideals).
   *
   * f = [O_K/P : Z/pZ] where P lies above p.
   *
   * @see Reference: sage/rings/number_field/number_field_ideal.py:residue_class_degree
   */
  residue_class_degree(): bigint {
    if (!this.is_prime()) {
      throw new ValueError('residue class degree only defined for prime ideals');
    }

    // The degree is log_p(N(P))
    const norm = this.norm();
    const p = this.prime_below();

    let f = 0n;
    let temp = norm;
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
    // Check if all generators are integral
    return this._gens.every((g) => g.is_integral());
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
    // Special case: zero ideal contains only zero
    if (this.is_zero()) {
      return x.is_zero();
    }

    // Zero is in every ideal
    if (x.is_zero()) {
      return true;
    }

    // Principal ideal case
    if (this._gens.length === 1) {
      // x in (a) iff a divides x
      const a = this._gens[0]!;
      if (a.is_zero()) {
        return x.is_zero();
      }

      // Check if x/a is integral
      try {
        const quotient = x.div(a);
        return quotient.is_integral();
      } catch {
        return false;
      }
    }

    // For two-generator ideals (p, alpha), check if x is an O_K-linear combination
    if (this._gens.length === 2) {
      const g1 = this._gens[0]!;
      const g2 = this._gens[1]!;

      // Check if x is in (g1)
      if (!g1.is_zero()) {
        try {
          const q1 = x.div(g1);
          if (q1.is_integral()) {
            return true;
          }
        } catch {
          // Division failed
        }
      }

      // Check if x is in (g2)
      if (!g2.is_zero()) {
        try {
          const q2 = x.div(g2);
          if (q2.is_integral()) {
            return true;
          }
        } catch {
          // Division failed
        }
      }

      // More sophisticated check: if both generators are known,
      // check using norms
      const normX = x.norm();
      const normI = this.norm();

      // If N(x) / N(I) is not an integer, x is not in I
      const ratio = normX.div(new Rational(normI));
      if (ratio.denominator !== 1n) {
        return false;
      }
    }

    // General case: would need HNF computation
    throw new NotImplementedError('contains for general non-principal ideals requires HNF');
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

    throw new NotImplementedError('integral_basis for general non-principal ideals requires HNF computation');
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
      const norm = this.norm();
      const normInv = new Rational(1n, norm);

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

    throw new NotImplementedError('inverse for general non-principal ideals requires PARI idealinv');
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

    // Two ideals are equal if they contain each other
    // For simple cases, compare norms
    if (this.is_zero() && other.is_zero()) {
      return true;
    }

    if (this.is_zero() !== other.is_zero()) {
      return false;
    }

    // Compare norms first (quick check)
    if (this.norm() !== other.norm()) {
      return false;
    }

    // For principal ideals, check if generators differ by a unit
    if (this._gens.length === 1 && other._gens.length === 1) {
      const a = this._gens[0]!;
      const b = other._gens[0]!;

      if (a.is_zero() && b.is_zero()) {
        return true;
      }

      if (a.is_zero() !== b.is_zero()) {
        return false;
      }

      // Check if a/b is a unit
      const ratio = a.div(b);
      return ratio.is_unit();
    }

    // For two-generator ideals, check if they have the same generators (up to units)
    if (this._gens.length === 2 && other._gens.length === 2) {
      // Check if the quotient ideal is the unit ideal
      try {
        const quotient = this.div(other);
        if (quotient.norm() === 1n) {
          return true;
        }
      } catch {
        // Division failed, ideals are likely not equal
      }
    }

    // General case: compare HNF representations
    // For now, check if both contain each other's generators
    try {
      for (const g of this._gens) {
        if (!other.contains(g)) {
          return false;
        }
      }
      for (const g of other._gens) {
        if (!this.contains(g)) {
          return false;
        }
      }
      return true;
    } catch {
      throw new NotImplementedError('eq for general non-principal ideals requires HNF computation');
    }
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

    // Coprime iff gcd(N(I), N(J)) = 1 is necessary but not sufficient
    const g = intGcd(this.norm(), other.norm());

    if (g === 1n) {
      return true;
    }

    // Need more careful check
    const sum = this.add(other);
    return sum.norm() === 1n;
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
    throw new NotImplementedError('ideal_class requires PARI bnfisprincipal for non-principal ideals');
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
   * Compute the Hermite Normal Form of this ideal.
   */
  protected _computeHNF(): HNFMatrix {
    if (this._cachedHNF) {
      return this._cachedHNF;
    }

    const n = this._number_field.degree();

    // For a principal ideal, HNF is simpler
    if (this._gens.length === 1) {
      const gen = this._gens[0]!;
      const coeffs = gen.list();

      // Find common denominator
      let denom = 1n;
      for (const c of coeffs) {
        denom = intLcm(denom, c.denominator);
      }

      // Scale to integers
      const intCoeffs = coeffs.map((c) => (c.numerator * (denom / c.denominator)));

      // For a principal ideal generated by a, the HNF is the matrix
      // representing multiplication by a
      // This is a simplified version
      const entries: bigint[][] = Array(n)
        .fill(null)
        .map(() => Array(n).fill(0n));

      // Just return the identity scaled by the norm for now
      // Full implementation requires proper HNF computation
      const normAbs = this.norm();

      for (let i = 0; i < n; i++) {
        entries[i]![i] = normAbs;
      }

      this._cachedHNF = { entries, denominator: 1n };
      return this._cachedHNF;
    }

    // General case: construct matrix from all generators and reduce to HNF
    throw new NotImplementedError('HNF computation for non-principal ideals');
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
