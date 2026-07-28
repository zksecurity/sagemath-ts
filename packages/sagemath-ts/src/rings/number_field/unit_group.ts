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

import { isqrt } from '../../arith/misc.js';
import { NotImplementedError, ValueError } from '../../errors.js';
import { Rational } from '../rational.js';
import type { NumberField } from './number_field.js';
import type { NumberFieldElement } from './number_field.js';
import { is_root_of_unity } from './number_field_element.js';
import { quadunit } from './pari_nf.js';

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
    const [r1] = K.signature();

    // For a real quadratic field the regulator is log|epsilon| at either real
    // embedding.
    if (K.degree() === 2 && r1 === 2) {
      return Math.abs(quadraticLogAbs(units[0]!));
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
    const [r1] = K.signature();

    if (K.degree() === 2) {
      if (r1 === 2) {
        // Real quadratic: two real embeddings, rank 1; only r1 + r2 - 1 = 1
        // column is needed.
        return [[quadraticLogAbs(units[0]!)]];
      }
      // Imaginary quadratic: no fundamental units
      return [];
    }

    throw new NotImplementedError('log_embedding: requires embedding computation for degree > 2');
  }

  /**
   * Check if u is a unit of the ring of integers, i.e. an algebraic integer of
   * norm +/-1.  (`NumberFieldElement.is_unit()` answers for the *field*, where
   * every nonzero element is a unit.)
   */
  contains(u: NumberFieldElement): boolean {
    return u.is_integral_unit();
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

    // Rank 1 (real quadratic): u = +/- eps^k.  The exponent is found by a
    // doubling/binary search driven by the *exact* sign test of
    // `quadraticSign` (a floating point log of |u| is useless here: for
    // u = eps^-8 in Q(sqrt 2), 577 - 408*sqrt(2) cancels to nothing in a
    // double).  The result is then verified exactly.
    if (rank === 1 && K.degree() === 2) {
      const eps = this.fundamental_units()[0]!;
      // |u| as a field element, and the sign that was stripped
      const su = quadraticSign(u);
      const abs = su < 0 ? u.neg() : u;
      const e0 = su < 0 ? 1n : 0n;
      const k = quadraticUnitExponent(abs, eps);
      if (k !== null) {
        const cand = eps.pow(k);
        if (cand.eq(abs)) {
          result[0] = e0;
          result[1] = k;
          return result;
        }
      }
      throw new ValueError('failed to write the unit as a power of the fundamental unit');
    }

    // For higher rank we need the full log embedding matrix and lattice
    // reduction over the r real/complex places.
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
      structure += ' x Z';
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

  const [, r2] = K.signature();

  if (r2 === 1) {
    // Imaginary quadratic field: the unit group is exactly the roots of unity,
    // of order 6 for disc -3, 4 for disc -4 and 2 otherwise.
    const disc = K.discriminant();
    const [torsionOrder, torsionGen] = quadraticTorsion(K, disc);
    return new UnitGroup(K, torsionOrder, torsionGen, []);
  }

  // Real quadratic field: rank 1, torsion {+/-1}, and one fundamental unit,
  // which PARI's `quadunit` reads off the continued fraction expansion of
  // sqrt(disc) (equivalently, the fundamental solution of Pell's equation).
  const torsionGen = K.__call__(-1n);
  return new UnitGroup(K, 2n, torsionGen, [realQuadraticFundamentalUnit(K)]);
}

/**
 * The fundamental unit of the ring of integers of a *real* quadratic field,
 * as an element of `K` in the power basis of its own generator.
 *
 * SageMath obtains this from `K.pari_bnf().bnf_get_fu()`; for a quadratic
 * field PARI's `bnfinit` in turn reduces to `quadunit(disc(K))`, which is
 * exact and needs no class-group machinery.  `quadunit` returns `[u, v]` with
 * `epsilon = u + v * w_D`; here we translate `w_D` into the power basis of the
 * field's own generator `alpha`.
 *
 * If `x^2 + b x + c` is the (monic) defining polynomial then
 * `sqrtDelta := 2*alpha + b` squares to `Delta = b^2 - 4c`, and
 * `Delta = f^2 * D` with `D = disc(O_K)` and `f` a positive rational, so
 * `sqrt(D) = sqrtDelta / f`.
 *
 * The unit returned is `> 1` at the real embedding `alpha |-> (-b+sqrt(Delta))/2`
 * (the larger root), matching PARI's normalisation; e.g. Sage's
 * `NumberField(1/2*x^2 - 1/6).units()` is `(3*a + 2,)`, i.e. `2 + sqrt(3)`.
 *
 * @see Reference: reference/pari/src/basemath/quad.c:476 (quadunit)
 * @see Reference: sage/rings/number_field/number_field.py:7207 (units)
 */
export function realQuadraticFundamentalUnit(K: NumberField): NumberFieldElement {
  if (K.degree() !== 2) {
    throw new ValueError('not a quadratic field');
  }
  const [r1] = K.signature();
  if (r1 !== 2) {
    throw new ValueError('not a real quadratic field');
  }

  const D = K.discriminant();
  const [u, v] = quadunit(D);

  const poly = K.defining_polynomial();
  const b = poly.getCoeff(1);
  const c = poly.getCoeff(0);
  const delta = b.mul(b).sub(c.mul(new Rational(4n)));
  const f = rationalSqrt(delta.div(new Rational(D)));
  if (f === null) {
    throw new ValueError('disc(K) does not differ from disc(f) by a rational square');
  }

  // sqrtD = (2*alpha + b) / f
  const sqrtD = K.gen().scalarMul(new Rational(2n)).add(K.__call__(b)).scalarMul(f.inv());

  // w_D = sqrt(D)/2 (D even) or (1 + sqrt(D))/2 (D odd); epsilon = u + v*w_D.
  const half = new Rational(1n, 2n);
  const w = (D & 1n) === 1n ? K.one().add(sqrtD).scalarMul(half) : sqrtD.scalarMul(half);
  const eps = K.__call__(u).add(w.scalarMul(new Rational(v)));

  // Certify: a unit of O_K has norm +/-1 and is an algebraic integer.
  const nrm = eps.norm();
  if (nrm.denominator !== 1n || (nrm.numerator !== 1n && nrm.numerator !== -1n)) {
    throw new ValueError(`computed fundamental unit has norm ${nrm}, expected +/-1`);
  }
  if (!eps.is_integral()) {
    throw new ValueError('computed fundamental unit is not an algebraic integer');
  }
  return eps;
}

/**
 * The group of roots of unity of an imaginary quadratic field, as
 * `[order, generator]`.
 *
 * `sqrt(disc) = (2*alpha + b) / f` where `x^2 + b x + c` is the defining
 * polynomial, `Delta = b^2 - 4c` its discriminant and `f = sqrt(Delta/disc)`
 * the conductor of `Z[alpha]`.  A primitive 4th root of unity is
 * `sqrt(disc)/2` (disc = -4) and a primitive 6th root is `(1 + sqrt(disc))/2`
 * (disc = -3).
 */
function quadraticTorsion(K: NumberField, disc: bigint): [bigint, NumberFieldElement] {
  if (disc !== -3n && disc !== -4n) {
    return [2n, K.__call__(-1n)];
  }
  const poly = K.defining_polynomial();
  const b = poly.getCoeff(1);
  const c = poly.getCoeff(0);
  const delta = b.mul(b).sub(c.mul(new Rational(4n)));
  const fSquared = delta.div(new Rational(disc));
  const f = rationalSqrt(fSquared);
  if (f === null) {
    throw new ValueError('the conductor of the equation order is not a perfect square');
  }
  // sqrtDisc = (2 alpha + b) / f
  const sqrtDisc = K.gen().scalarMul(new Rational(2n)).add(K.__call__(b)).scalarMul(f.inv());

  let order: bigint;
  let zeta: NumberFieldElement;
  if (disc === -4n) {
    order = 4n;
    zeta = sqrtDisc.scalarMul(new Rational(1n, 2n));
  } else {
    order = 6n;
    zeta = K.one().add(sqrtDisc).scalarMul(new Rational(1n, 2n));
  }
  // Verify: zeta has exact multiplicative order `order`.
  if (!zeta.pow(order).is_one()) {
    throw new ValueError('computed torsion generator is not a root of unity');
  }
  for (const p of [2n, 3n]) {
    if (order % p === 0n && zeta.pow(order / p).is_one()) {
      throw new ValueError('computed torsion generator does not have the expected order');
    }
  }
  return [order, zeta];
}

/**
 * The exact sign of `x` at the real embedding `alpha |-> (-b + sqrt(Delta))/2`
 * of a real quadratic field with monic defining polynomial `x^2 + b x + c`
 * and `Delta = b^2 - 4c > 0`.
 *
 * Writing `x = p + q*alpha = A + B*sqrt(Delta)` with `A = p - q b/2`,
 * `B = q/2`, the sign is decided by comparing `A^2` with `B^2 Delta` when the
 * two terms disagree -- entirely in exact rational arithmetic, so there is no
 * cancellation.
 */
function quadraticSign(x: NumberFieldElement): number {
  const K = x.parent();
  const poly = K.defining_polynomial();
  const b = poly.getCoeff(1);
  const c = poly.getCoeff(0);
  const delta = b.mul(b).sub(c.mul(new Rational(4n)));
  if (delta.sign <= 0n) {
    throw new ValueError('quadraticSign requires a real quadratic field');
  }
  const coeffs = x.list();
  const p = coeffs[0] ?? Rational.zero();
  const q = coeffs[1] ?? Rational.zero();
  const A = p.sub(q.mul(b).mul(new Rational(1n, 2n)));
  const B = q.mul(new Rational(1n, 2n));
  const sA = Number(A.sign);
  const sB = Number(B.sign);
  if (sB === 0) return sA;
  if (sA === 0) return sB;
  if (sA === sB) return sA;
  // Opposite signs: |A| vs |B| sqrt(Delta), i.e. A^2 vs B^2 Delta.
  const lhs = A.mul(A);
  const rhs = B.mul(B).mul(delta);
  const d = lhs.sub(rhs);
  if (d.sign === 0n) return 0; // impossible: sqrt(Delta) is irrational
  return d.sign > 0n ? sA : sB;
}

/** Exact comparison `x <=> y` at the same real embedding. */
function quadraticCmp(x: NumberFieldElement, y: NumberFieldElement): number {
  return quadraticSign(x.sub(y));
}

/**
 * The integer `k` with `eps^k = a`, for `a` a totally positive-at-the-chosen-
 * embedding unit and `eps > 1` the fundamental unit.  Returns `null` when no
 * such `k` exists.  Uses O(log|k|) multiplications (doubling then a greedy
 * descent), never a linear scan.
 */
function quadraticUnitExponent(a: NumberFieldElement, eps: NumberFieldElement): bigint | null {
  const K = a.parent();
  const one = K.one();
  if (a.eq(one)) return 0n;
  if (quadraticSign(a) <= 0) return null;
  // Work with b >= 1 so that k >= 0, remembering whether to negate at the end.
  const inverted = quadraticCmp(a, one) < 0;
  let b = inverted ? a.inv() : a;
  // Doubling: powers[i] = eps^(2^i), stop once it exceeds b.
  const powers: NumberFieldElement[] = [eps];
  while (quadraticCmp(powers[powers.length - 1]!, b) < 0) {
    const next = powers[powers.length - 1]!;
    powers.push(next.mul(next));
    if (powers.length > 4096) return null;
  }
  let k = 0n;
  for (let i = powers.length - 1; i >= 0; i--) {
    const pw = powers[i]!;
    if (quadraticCmp(pw, b) <= 0) {
      b = b.div(pw);
      k += 1n << BigInt(i);
    }
  }
  if (!b.eq(one)) return null;
  return inverted ? -k : k;
}

/** Exact square root of a nonnegative rational, or null if it is not a square. */
function rationalSqrt(x: Rational): Rational | null {
  if (x.numerator < 0n) return null;
  const sn = isqrt(x.numerator);
  const sd = isqrt(x.denominator);
  if (sn * sn !== x.numerator || sd * sd !== x.denominator) return null;
  return new Rational(sn, sd);
}

/** `log|n|` for a nonzero bigint, without overflowing a double. */
function logAbsBigInt(n: bigint): number {
  const a = n < 0n ? -n : n;
  if (a === 0n) return Number.NEGATIVE_INFINITY;
  const bits = a.toString(2).length;
  if (bits <= 53) return Math.log(Number(a));
  const shift = BigInt(bits - 53);
  return Math.log(Number(a >> shift)) + Number(shift) * Math.LN2;
}

/** `log|r|` for a nonzero rational, without overflowing a double. */
function logAbsRational(r: Rational): number {
  return logAbsBigInt(r.numerator) - logAbsBigInt(r.denominator);
}

/**
 * `log|sigma(x)|` at the real embedding `alpha |-> (-b + sqrt(Delta))/2` of a
 * real quadratic field, computed without overflow and without catastrophic
 * cancellation.
 *
 * `x = A + B sqrt(Delta)`; when `A` and `B sqrt(Delta)` have the same sign the
 * two logs are combined with the usual log-sum-exp.  When they have opposite
 * signs the sum cancels (for `eps^-8` in `Q(sqrt 2)`, `577 - 408 sqrt 2`
 * loses every significant digit in a double), so we use
 * `|sigma(x)| = |N(x)| / |sigma'(x)|` instead: the conjugate `A - B sqrt(Delta)`
 * then has same-sign terms.
 */
function quadraticLogAbs(x: NumberFieldElement): number {
  const K = x.parent();
  const poly = K.defining_polynomial();
  const b = poly.getCoeff(1);
  const c = poly.getCoeff(0);
  const delta = b.mul(b).sub(c.mul(new Rational(4n)));
  if (delta.sign <= 0n) {
    throw new ValueError('quadraticLogAbs requires a real quadratic field');
  }
  const coeffs = x.list();
  const p = coeffs[0] ?? Rational.zero();
  const q = coeffs[1] ?? Rational.zero();
  const A = p.sub(q.mul(b).mul(new Rational(1n, 2n)));
  const B = q.mul(new Rational(1n, 2n));

  const combine = (opposite: boolean): number | null => {
    // log|A| and log|B sqrt(Delta)|
    if (B.sign === 0n) return A.sign === 0n ? Number.NEGATIVE_INFINITY : logAbsRational(A);
    if (A.sign === 0n) return logAbsRational(B) + 0.5 * logAbsRational(delta);
    const lA = logAbsRational(A);
    const lB = logAbsRational(B) + 0.5 * logAbsRational(delta);
    const m = Math.max(lA, lB);
    const s = Math.exp(lA - m) + (opposite ? -1 : 1) * Math.exp(lB - m);
    if (opposite && Math.abs(s) < 1e-8) return null; // cancellation
    return m + Math.log(Math.abs(s));
  };

  const sA = Number(A.sign);
  const sB = Number(B.sign);
  const opposite = sA !== 0 && sB !== 0 && sA !== sB;
  const direct = combine(opposite);
  if (direct !== null) return direct;

  // Cancellation: log|x| = log|N(x)| - log|conjugate(x)|.
  const nrm = x.norm();
  if (nrm.sign === 0n) return Number.NEGATIVE_INFINITY;
  const lA = logAbsRational(A);
  const lB = logAbsRational(B) + 0.5 * logAbsRational(delta);
  const m = Math.max(lA, lB);
  const conj = m + Math.log(Math.exp(lA - m) + Math.exp(lB - m));
  return logAbsRational(nrm) - conj;
}
