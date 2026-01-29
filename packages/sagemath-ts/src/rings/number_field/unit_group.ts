/**
 * @module sage/rings/number_field/unit_group
 * @description Unit groups of number fields
 *
 * Port of: sage/rings/number_field/unit_group.py
 * Reference: reference/sage/src/sage/rings/number_field/unit_group.py
 *
 * NOTE: Full unit group computation requires PARI's bnfinit.
 * This implementation provides the structure and basic operations.
 * For full PARI compatibility, see DEVIATIONS.md.
 */

import { NotImplementedError, ValueError } from '../../errors.js';
import type { NumberField } from './number_field.js';
import { NumberFieldElement } from './number_field.js';
import { Rational } from '../rational.js';
import { is_root_of_unity } from './number_field_element.js';

/**
 * The unit group of a number field.
 *
 * The unit group O_K^* of the ring of integers O_K is isomorphic to:
 *   O_K^* ≅ μ_K × Z^r
 *
 * where:
 * - μ_K is the finite group of roots of unity in K
 * - r = r1 + r2 - 1 is the rank (by Dirichlet's unit theorem)
 * - r1 = number of real embeddings
 * - r2 = number of pairs of complex conjugate embeddings
 *
 * @see Reference: sage/rings/number_field/unit_group.py:UnitGroup
 */
export class UnitGroup {
  private readonly _number_field: NumberField;
  private readonly _torsion_order: bigint;
  private readonly _torsion_generator?: NumberFieldElement;
  private readonly _fundamental_units: readonly NumberFieldElement[];

  constructor(
    number_field: NumberField,
    torsion_order: bigint = 2n,
    torsion_generator?: NumberFieldElement,
    fundamental_units: NumberFieldElement[] = []
  ) {
    this._number_field = number_field;
    this._torsion_order = torsion_order;
    this._torsion_generator = torsion_generator;
    this._fundamental_units = Object.freeze([...fundamental_units]);
  }

  /**
   * Return the number field this unit group belongs to.
   */
  number_field(): NumberField {
    return this._number_field;
  }

  /**
   * Return the rank of the free part of the unit group.
   *
   * By Dirichlet's unit theorem: rank = r1 + r2 - 1
   * where r1, r2 are the signature of K.
   */
  rank(): number {
    const [r1, r2] = this._number_field.signature();
    return r1 + r2 - 1;
  }

  /**
   * Return the order of the torsion subgroup (roots of unity).
   */
  torsion_order(): bigint {
    return this._torsion_order;
  }

  /**
   * Alias for torsion_order.
   */
  zeta_order(): bigint {
    return this._torsion_order;
  }

  /**
   * Return a generator of the torsion subgroup.
   */
  torsion_generator(): NumberFieldElement {
    if (!this._torsion_generator) {
      throw new NotImplementedError('torsion generator computation requires PARI');
    }
    return this._torsion_generator;
  }

  /**
   * Return the roots of unity in this number field.
   */
  roots_of_unity(): NumberFieldElement[] {
    const gen = this.torsion_generator();
    const order = Number(this._torsion_order);
    const roots: NumberFieldElement[] = [];
    let current = gen.parent().one();

    for (let i = 0; i < order; i++) {
      roots.push(current);
      current = current.mul(gen);
    }

    return roots;
  }

  /**
   * Return a primitive n-th root of unity, if it exists.
   */
  zeta(n?: bigint): NumberFieldElement {
    if (n === undefined || n === 1n) {
      return this._number_field.one();
    }

    if (n === 2n) {
      // -1 is always a root of unity
      return this._number_field.__call__(-1n);
    }

    if (this._torsion_order % n !== 0n) {
      throw new ValueError(`no ${n}-th root of unity in this field`);
    }

    const gen = this.torsion_generator();
    const exp = this._torsion_order / n;
    return gen.pow(exp);
  }

  /**
   * Return the fundamental units.
   *
   * These are generators for the free part of the unit group.
   * There are r = r1 + r2 - 1 fundamental units.
   */
  fundamental_units(): NumberFieldElement[] {
    if (this._fundamental_units.length === 0 && this.rank() > 0) {
      throw new NotImplementedError('fundamental units computation requires PARI');
    }
    return [...this._fundamental_units];
  }

  /**
   * Return all generators (torsion + fundamental units).
   */
  gens(): NumberFieldElement[] {
    const result: NumberFieldElement[] = [];
    if (this._torsion_generator) {
      result.push(this._torsion_generator);
    }
    result.push(...this._fundamental_units);
    return result;
  }

  /**
   * Return the number of generators.
   */
  ngens(): number {
    return (this._torsion_generator ? 1 : 0) + this._fundamental_units.length;
  }

  /**
   * Return the i-th generator.
   */
  gen(i: number): NumberFieldElement {
    const gens = this.gens();
    if (i < 0 || i >= gens.length) {
      throw new ValueError(`generator index ${i} out of range`);
    }
    return gens[i]!;
  }

  /**
   * Return the regulator of this unit group.
   *
   * The regulator is the absolute value of the determinant of the
   * (r x r) matrix formed from the log embeddings of the fundamental units,
   * where r = r1 + r2 - 1 is the rank.
   *
   * For quadratic fields:
   * - Imaginary: regulator = 1 (no fundamental units)
   * - Real: regulator = log|epsilon| where epsilon is the fundamental unit
   *
   * Note: This requires numerical computation.
   */
  regulator(): number {
    const rank = this.rank();

    if (rank === 0) {
      // No fundamental units, regulator is 1 by convention
      return 1;
    }

    const units = this.fundamental_units();
    if (units.length === 0) {
      throw new NotImplementedError('regulator: fundamental units not computed');
    }

    const K = this._number_field;
    const [r1, r2] = K.signature();

    // Build the log embedding matrix
    // For each fundamental unit u_i, compute:
    // - log|sigma_j(u_i)| for j = 1, ..., r1 (real embeddings)
    // - 2 * log|sigma_j(u_i)| for j = 1, ..., r2 (complex embeddings)

    // For quadratic fields, the regulator is simply log|epsilon|
    if (K.degree() === 2 && r1 === 2) {
      // Real quadratic field
      const u = units[0]!;
      // Get coefficients [a, b] where u = a + b*sqrt(d)
      const coeffs = u.list();
      const a = Number(coeffs[0]!.numerator) / Number(coeffs[0]!.denominator);
      const b = Number(coeffs[1]!.numerator) / Number(coeffs[1]!.denominator);

      // Compute sqrt(d) from discriminant
      const disc = K.discriminant();
      const d = disc > 0n ? disc : -disc;
      const sqrtD = Math.sqrt(Number(d));

      // u = a + b*sqrt(d), |u| for real embedding
      const absU = Math.abs(a + b * sqrtD);
      return Math.log(absU);
    }

    // For higher degree fields, we would need to compute all embeddings
    // and form the full log embedding matrix
    throw new NotImplementedError('regulator: requires embedding computation for degree > 2');
  }

  /**
   * Return the logarithmic embedding matrix.
   *
   * This is the matrix with entries:
   * - L[i][j] = log|sigma_j(u_i)| for real embeddings (j <= r1)
   * - L[i][j] = 2*log|sigma_j(u_i)| for complex embeddings (j > r1)
   *
   * where u_i are the fundamental units.
   */
  log_embedding(): number[][] {
    const units = this.fundamental_units();
    const rank = this.rank();

    if (rank === 0) {
      return [];
    }

    if (units.length === 0) {
      throw new NotImplementedError('log_embedding: fundamental units not computed');
    }

    const K = this._number_field;
    const [r1, r2] = K.signature();

    // For quadratic fields
    if (K.degree() === 2) {
      if (r1 === 2) {
        // Real quadratic: 2 real embeddings, rank 1
        const u = units[0]!;
        const coeffs = u.list();
        const a = Number(coeffs[0]!.numerator) / Number(coeffs[0]!.denominator);
        const b = Number(coeffs[1]!.numerator) / Number(coeffs[1]!.denominator);

        const disc = K.discriminant();
        const d = disc > 0n ? disc : -disc;
        const sqrtD = Math.sqrt(Number(d));

        // Two embeddings: a + b*sqrt(d) and a - b*sqrt(d)
        const embed1 = a + b * sqrtD;
        const embed2 = a - b * sqrtD;

        // Log embedding matrix is [log|embed1|, log|embed2|]
        // But we only need r1 + r2 - 1 = 1 column for the regulator
        return [[Math.log(Math.abs(embed1))]];
      }
      // Imaginary quadratic: no fundamental units
      return [];
    }

    throw new NotImplementedError('log_embedding: requires embedding computation for degree > 2');
  }

  /**
   * Check if u is a unit.
   */
  contains(u: NumberFieldElement): boolean {
    return u.is_unit();
  }

  /**
   * Return the exponent vector of u with respect to the generators.
   *
   * This expresses u = ζ^e0 * u1^e1 * u2^e2 * ... * ur^er
   * where ζ is the torsion generator and u1,...,ur are fundamental units.
   *
   * Algorithm: First find the exponents for fundamental units using the
   * log embedding, then compute the torsion exponent.
   */
  log(u: NumberFieldElement): bigint[] {
    if (!this.contains(u)) {
      throw new ValueError('element is not a unit');
    }

    const ngens = this.ngens();
    const result: bigint[] = Array(ngens).fill(0n);

    if (ngens === 0) {
      return result;
    }

    const K = this._number_field;
    const rank = this.rank();

    // Handle the simple case of torsion-only (rank 0)
    if (rank === 0) {
      // Only torsion generator
      const torsion = this._torsion_order;
      const gen = this.torsion_generator();

      // Find e such that u = gen^e
      let current = K.one();
      for (let e = 0n; e < torsion; e++) {
        if (u.eq(current)) {
          result[0] = e;
          return result;
        }
        current = current.mul(gen);
      }
      throw new ValueError('failed to find torsion exponent');
    }

    // For rank > 0, we need to solve for fundamental unit exponents
    // This is done via the log embedding and linear algebra
    throw new NotImplementedError('log: requires numerical linear algebra for rank > 0');
  }

  /**
   * Return a unit from an exponent vector.
   */
  exp(v: bigint[]): NumberFieldElement {
    if (v.length !== this.ngens()) {
      throw new ValueError(`expected ${this.ngens()} exponents, got ${v.length}`);
    }

    const gens = this.gens();
    let result = this._number_field.one();

    for (let i = 0; i < v.length; i++) {
      result = result.mul(gens[i]!.pow(v[i]!));
    }

    return result;
  }

  /**
   * Check if two units are equivalent modulo torsion.
   *
   * Two units u and v are equivalent if u/v is a root of unity.
   */
  is_equivalent(u: NumberFieldElement, v: NumberFieldElement): boolean {
    if (!this.contains(u) || !this.contains(v)) {
      return false;
    }

    // Check if u/v is a root of unity
    // A unit w is a root of unity if and only if |N(w)| = 1 and all
    // embeddings have absolute value 1

    // Compute u/v
    const quotient = u.div(v);

    // Check if the quotient is a root of unity
    return is_root_of_unity(quotient);
  }

  /**
   * Return the order of this unit group (infinite if rank > 0).
   */
  order(): bigint | 'infinity' {
    if (this.rank() > 0) {
      return 'infinity';
    }
    return this._torsion_order;
  }

  /**
   * Check if this unit group is finite.
   */
  is_finite(): boolean {
    return this.rank() === 0;
  }

  toString(): string {
    const rank = this.rank();
    const torsion = this._torsion_order;

    let structure = `C${torsion}`;
    if (rank > 0) {
      structure += ` x Z`;
      if (rank > 1) {
        structure = `C${torsion} x ` + Array(rank).fill('Z').join(' x ');
      }
    }

    return `Unit group with structure ${structure} of ${this._number_field}`;
  }
}

/**
 * The S-unit group of a number field.
 *
 * For a set S of prime ideals, the S-units are elements u such that
 * v_P(u) = 0 for all primes P not in S.
 *
 * @see Reference: sage/rings/number_field/unit_group.py:S_units
 */
export class S_UnitGroup extends UnitGroup {
  private readonly _S: unknown[]; // Set of prime ideals

  constructor(
    number_field: NumberField,
    S: unknown[],
    torsion_order: bigint = 2n,
    torsion_generator?: NumberFieldElement,
    fundamental_units: NumberFieldElement[] = []
  ) {
    super(number_field, torsion_order, torsion_generator, fundamental_units);
    this._S = [...S];
  }

  /**
   * Return the set S of primes.
   */
  primes(): unknown[] {
    return [...this._S];
  }

  /**
   * Return the S-unit rank.
   *
   * rank(S-units) = rank(units) + |S|
   */
  override rank(): number {
    return super.rank() + this._S.length;
  }

  /**
   * Return the S-units generators.
   *
   * The generators consist of:
   * 1. The torsion generator (if any)
   * 2. The fundamental units (from the ordinary unit group)
   * 3. One generator for each prime in S (a uniformizer)
   *
   * This requires PARI's bnfunits for a complete implementation.
   */
  override gens(): NumberFieldElement[] {
    // Get the ordinary unit group generators
    const ordinaryGens = super.gens();

    // For each prime P in S, we need a generator for the P-part
    // These are elements with valuation 1 at P and 0 at other primes in S
    // Full computation requires PARI

    if (this._S.length > 0) {
      throw new NotImplementedError('S_UnitGroup.gens: full S-unit computation requires PARI');
    }

    return ordinaryGens;
  }

  override toString(): string {
    const rank = this.rank();
    const torsion = this.torsion_order();

    let structure = `C${torsion}`;
    if (rank > 0) {
      const zParts = Array(rank).fill('Z').join(' x ');
      structure = `C${torsion} x ${zParts}`;
    }

    return `S-unit group with structure ${structure} of ${this.number_field()} with S = (${this._S.length} primes)`;
  }
}

/**
 * Create a unit group for a quadratic field.
 *
 * For quadratic fields, we can compute units directly:
 * - Imaginary quadratic: only roots of unity, determined by discriminant
 * - Real quadratic: fundamental unit from continued fraction expansion
 */
export function quadraticUnitGroup(K: NumberField): UnitGroup {
  if (K.degree() !== 2) {
    throw new ValueError('not a quadratic field');
  }

  const [r1, r2] = K.signature();

  if (r2 === 1) {
    // Imaginary quadratic field
    // Roots of unity: -1 is always there
    // Additional roots of unity only for d = -1 (4th roots) and d = -3 (6th roots)
    const disc = K.discriminant();

    let torsionOrder: bigint;
    let torsionGen: NumberFieldElement;

    if (disc === -4n) {
      // Q(i): 4th roots of unity
      torsionOrder = 4n;
      torsionGen = K.gen(); // i
    } else if (disc === -3n) {
      // Q(zeta_3): 6th roots of unity
      torsionOrder = 6n;
      // (1 + sqrt(-3))/2 is a primitive 6th root
      torsionGen = K.__call__(new Rational(1n, 2n)).add(K.gen().scalarMul(new Rational(1n, 2n)));
    } else {
      // Only +/-1
      torsionOrder = 2n;
      torsionGen = K.__call__(-1n);
    }

    return new UnitGroup(K, torsionOrder, torsionGen, []);
  } else {
    // Real quadratic field
    // Has rank 1, so one fundamental unit
    // Finding the fundamental unit requires solving Pell's equation
    // or using continued fractions - complex to implement without PARI

    // For now, just set up the structure with torsion = 2 (just +/-1)
    const torsionOrder = 2n;
    const torsionGen = K.__call__(-1n);

    // The fundamental unit would be computed via:
    // 1. Continued fraction expansion of sqrt(d)
    // 2. Or solving x^2 - d*y^2 = +/-1

    // We return an incomplete unit group - fundamental unit not computed
    return new UnitGroup(K, torsionOrder, torsionGen, []);
  }
}
