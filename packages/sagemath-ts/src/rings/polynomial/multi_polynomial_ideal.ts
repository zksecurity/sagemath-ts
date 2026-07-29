/**
 * @module sage/rings/polynomial/multi_polynomial_ideal
 * @description Ideals in multivariate polynomial rings and Gröbner basis computation
 *
 * Port of: sage/rings/polynomial/multi_polynomial_ideal.py
 *
 * Implements Buchberger's algorithm for computing Gröbner bases.
 */

import { ArithmeticError, NotImplementedError, TypeError, ValueError } from '../../errors.js';
import {
  type Exponent,
  type MPolynomialElement,
  type TermOrder,
  exponentToKey,
  keyToExponent,
  totalDegree,
} from './multi_polynomial_element.js';
import type { MPolynomialRing } from './multi_polynomial_ring.js';
import type { CoefficientRing, RingElement } from './polynomial_element.js';

/**
 * Compute the least common multiple of two exponent tuples.
 * This gives the exponent of the LCM of the leading monomials.
 */
function lcmExponent(a: Exponent, b: Exponent): number[] {
  const n = Math.max(a.length, b.length);
  const result: number[] = [];
  for (let i = 0; i < n; i++) {
    result.push(Math.max(a[i] ?? 0, b[i] ?? 0));
  }
  return result;
}

/**
 * Check if exponent a is divisible by exponent b.
 * Returns true if for all i, a[i] >= b[i].
 */
function exponentDivisible(a: Exponent, b: Exponent): boolean {
  const n = Math.max(a.length, b.length);
  for (let i = 0; i < n; i++) {
    if ((a[i] ?? 0) < (b[i] ?? 0)) {
      return false;
    }
  }
  return true;
}

/**
 * Subtract exponent b from exponent a.
 * Assumes a is divisible by b.
 */
function exponentSub(a: Exponent, b: Exponent): number[] {
  const n = Math.max(a.length, b.length);
  const result: number[] = [];
  for (let i = 0; i < n; i++) {
    result.push((a[i] ?? 0) - (b[i] ?? 0));
  }
  return result;
}

/**
 * Compute the S-polynomial of two polynomials f and g.
 *
 * The S-polynomial is defined as:
 *   S(f, g) = (lcm(LM(f), LM(g)) / LT(f)) * f - (lcm(LM(f), LM(g)) / LT(g)) * g
 *
 * where LM is leading monomial and LT is leading term.
 */
export function sPolynomial<R extends CoefficientRing, E extends RingElement>(
  f: MPolynomialElement<R, E>,
  g: MPolynomialElement<R, E>
): MPolynomialElement<R, E> {
  if (f.isZero() || g.isZero()) {
    return f.parent.zero();
  }

  const ltF = f.leadingTerm();
  const ltG = g.leadingTerm();

  if (!ltF || !ltG) {
    return f.parent.zero();
  }

  const [lmF, lcF] = ltF;
  const [lmG, lcG] = ltG;

  // Compute LCM of leading monomials
  const lcmExp = lcmExponent(lmF, lmG);

  // Compute the multipliers
  const multF = exponentSub(lcmExp, lmF);
  const multG = exponentSub(lcmExp, lmG);

  // Create monomials for multiplication
  const monomF = f.parent.monomial(multF);
  const monomG = g.parent.monomial(multG);

  // S(f, g) = (lcm/LT(f)) * f - (lcm/LT(g)) * g
  // We need to multiply by inverse of leading coefficients to make them equal
  // S(f, g) = monomF * f / lcF - monomG * g / lcG

  // For fields, we can do: S(f, g) = lcG * monomF * f - lcF * monomG * g
  const term1 = monomF.mul(f).scalarMul(lcG);
  const term2 = monomG.mul(g).scalarMul(lcF);

  return term1.sub(term2);
}

/**
 * Raise SageMath's error when the base ring of `f` is not a field.
 *
 * SageMath: `multi_polynomial_element.py:2488` —
 * `if not k.is_field(): raise TypeError("Can only reduce polynomials over fields.")`.
 * Without this guard, multivariate division over e.g. ZZ truncates the
 * quotient coefficient to zero, the subtrahend is zero, `p` never changes and
 * the division loop spins forever (audit H14).
 */
function requireFieldBaseRing<R extends CoefficientRing, E extends RingElement>(
  f: MPolynomialElement<R, E>
): void {
  const k = f.parent.base_ring as { is_field?: () => boolean };
  if (typeof k.is_field === 'function' && !k.is_field()) {
    throw new TypeError('Can only reduce polynomials over fields.');
  }
}

/**
 * Reduce polynomial f with respect to a set of polynomials G.
 *
 * Returns the remainder when f is divided by the polynomials in G,
 * using multivariate polynomial division.
 *
 * SageMath: `MPolynomial_polydict.reduce` (`multi_polynomial_element.py:2478`).
 *
 * @throws {TypeError} If the base ring is not a field
 * @see Deviation: Polynomials — Printing, Factor Shape, Term Orders and Base Rings
 */
export function reduce<R extends CoefficientRing, E extends RingElement>(
  f: MPolynomialElement<R, E>,
  G: MPolynomialElement<R, E>[]
): MPolynomialElement<R, E> {
  requireFieldBaseRing(f);

  if (G.length === 0 || f.isZero()) {
    return f;
  }

  let r = f.parent.zero();
  let p = f;

  while (!p.isZero()) {
    let divisionOccurred = false;

    for (const g of G) {
      if (g.isZero()) continue;

      const ltP = p.leadingTerm();
      const ltG = g.leadingTerm();

      if (!ltP || !ltG) continue;

      const [lmP, lcP] = ltP;
      const [lmG, lcG] = ltG;

      // Check if LM(g) divides LM(p)
      if (exponentDivisible(lmP, lmG)) {
        // Compute quotient monomial
        const quotExp = exponentSub(lmP, lmG);
        const quotMon = f.parent.monomial(quotExp);

        // Compute quotient coefficient: lcP / lcG
        // For fields, this is straightforward
        const quotCoeff = lcP.div(lcG) as E;

        // Base rings without a working is_field() would otherwise loop
        // forever here: if the division is not exact the leading term of p
        // is not cancelled and p never shrinks (audit H14).
        if (!(quotCoeff.mul(lcG) as E).eq(lcP)) {
          throw new TypeError('Can only reduce polynomials over fields.');
        }

        // p = p - quotCoeff * quotMon * g
        const subtrahend = quotMon.mul(g).scalarMul(quotCoeff);
        p = p.sub(subtrahend);
        divisionOccurred = true;
        break;
      }
    }

    if (!divisionOccurred) {
      // No divisor found, move leading term to remainder
      const lt = p.leadingTerm();
      if (lt) {
        const [lm, lc] = lt;
        const ltPoly = f.parent.monomial(lm as number[]).scalarMul(lc);
        r = r.add(ltPoly);
        p = p.sub(ltPoly);
      } else {
        break;
      }
    }
  }

  return r;
}

/**
 * Compute a Gröbner basis for an ideal using Buchberger's algorithm.
 *
 * @param generators - The generators of the ideal
 * @param options - Options for the computation
 * @returns A Gröbner basis for the ideal
 *
 * @example
 * ```typescript
 * const [R, [x, y]] = MPolynomialRingConstructor(QQ, ['x', 'y']);
 * const f = x.pow(2).sub(y);
 * const g = x.mul(y).sub(x);
 * const gb = groebner_basis([f, g]);
 * ```
 *
 * @throws {TypeError} If the base ring is not a field
 * @throws {ArithmeticError} If `maxIterations` S-pairs are exhausted before
 *   the basis closes up. Buchberger's algorithm always terminates, so this
 *   only signals that the budget was too small; returning the partial set
 *   would hand back something that is *not* a Gröbner basis (audit M19).
 * @see Deviation: Polynomials — Printing, Factor Shape, Term Orders and Base Rings
 */
export function groebner_basis<R extends CoefficientRing, E extends RingElement>(
  generators: MPolynomialElement<R, E>[],
  options?: {
    /** Whether to interreduce the result (default: true) */
    interreduce?: boolean;
    /** Maximum number of iterations (default: 10000) */
    maxIterations?: number;
  }
): MPolynomialElement<R, E>[] {
  if (generators.length === 0) {
    return [];
  }

  requireFieldBaseRing(generators[0]!);

  const interreduce = options?.interreduce ?? true;
  const maxIterations = options?.maxIterations ?? 10000;

  // Filter out zero polynomials and make monic
  let G = generators
    .filter((f) => !f.isZero())
    .map((f) => {
      const lt = f.leadingTerm();
      if (lt) {
        const [, lc] = lt;
        // Make monic by dividing by leading coefficient
        try {
          const lcInv = lc.inv() as E;
          return f.scalarMul(lcInv);
        } catch {
          return f;
        }
      }
      return f;
    });

  if (G.length === 0) {
    return [];
  }

  // Buchberger's algorithm
  // Keep track of pairs we've already processed
  const processedPairs = new Set<string>();

  const makePairKey = (i: number, j: number): string => {
    return i < j ? `${i},${j}` : `${j},${i}`;
  };

  let iterations = 0;
  let completed = false;

  while (iterations < maxIterations) {
    iterations++;

    // Find a pair (f, g) that hasn't been processed
    let newPairFound = false;
    let fi = -1;
    let fj = -1;

    for (let i = 0; i < G.length && !newPairFound; i++) {
      for (let j = i + 1; j < G.length && !newPairFound; j++) {
        const key = makePairKey(i, j);
        if (!processedPairs.has(key)) {
          processedPairs.add(key);
          fi = i;
          fj = j;
          newPairFound = true;
        }
      }
    }

    if (!newPairFound) {
      // All pairs processed, we're done
      completed = true;
      break;
    }

    const f = G[fi]!;
    const g = G[fj]!;

    // Compute S-polynomial
    const s = sPolynomial(f, g);

    // Reduce S-polynomial with respect to G
    const r = reduce(s, G);

    // If remainder is non-zero, add it to the basis
    if (!r.isZero()) {
      // Make monic
      let rMonic = r;
      const lt = r.leadingTerm();
      if (lt) {
        const [, lc] = lt;
        try {
          const lcInv = lc.inv() as E;
          rMonic = r.scalarMul(lcInv);
        } catch {
          // If we can't invert, just use r as-is
        }
      }

      G.push(rMonic);
    }
  }

  if (!completed) {
    // Never return a truncated set: it is not a Gröbner basis, and callers
    // (contains/reduce/dimension) would silently produce wrong answers.
    throw new ArithmeticError(
      `groebner_basis: exhausted maxIterations (${maxIterations}) with ` +
        `${G.length} basis elements and unprocessed S-pairs remaining; ` +
        'the result would not be a Gröbner basis'
    );
  }

  // Interreduce the basis to get a minimal/reduced Gröbner basis
  if (interreduce) {
    G = interreduceBasis(G);
  }

  return G;
}

/**
 * Interreduce a Gröbner basis to get a reduced Gröbner basis.
 *
 * A reduced Gröbner basis has the properties:
 * 1. Each polynomial is monic
 * 2. No leading monomial divides any other leading monomial
 * 3. No monomial in any polynomial is divisible by the leading monomial of another polynomial
 */
function interreduceBasis<R extends CoefficientRing, E extends RingElement>(
  G: MPolynomialElement<R, E>[]
): MPolynomialElement<R, E>[] {
  if (G.length === 0) return [];

  // First pass: remove polynomials whose leading monomial is divisible by another's
  const result: MPolynomialElement<R, E>[] = [];

  for (let i = 0; i < G.length; i++) {
    const f = G[i]!;
    const ltF = f.leadingTerm();
    if (!ltF) continue;

    const [lmF] = ltF;
    let isRedundant = false;

    for (let j = 0; j < G.length; j++) {
      if (i === j) continue;
      const g = G[j]!;
      const ltG = g.leadingTerm();
      if (!ltG) continue;

      const [lmG] = ltG;

      // Check if LM(g) divides LM(f) and g is "smaller" (to avoid removing both)
      if (j < i && exponentDivisible(lmF, lmG)) {
        isRedundant = true;
        break;
      }
    }

    if (!isRedundant) {
      result.push(f);
    }
  }

  // Second pass: reduce each polynomial by the others
  const reduced: MPolynomialElement<R, E>[] = [];

  for (let i = 0; i < result.length; i++) {
    const f = result[i]!;
    const others = [...result.slice(0, i), ...result.slice(i + 1), ...reduced];

    let r = reduce(f, others);

    // Make monic
    const lt = r.leadingTerm();
    if (lt) {
      const [, lc] = lt;
      try {
        const lcInv = lc.inv() as E;
        r = r.scalarMul(lcInv);
      } catch {
        // If we can't invert, just use r as-is
      }
    }

    if (!r.isZero()) {
      reduced.push(r);
    }
  }

  return reduced;
}

/**
 * Enumerate the subsets of {0, ..., n-1} of a given size, in lexicographic
 * order. Used by `MPolynomialIdeal.dimension` to reproduce the order in which
 * SageMath's `Set(...).subsets()` iterator yields subsets (by increasing size).
 */
function* subsetsOfSize(n: number, size: number): Generator<number[]> {
  const current: number[] = [];
  function* rec(start: number): Generator<number[]> {
    if (current.length === size) {
      yield [...current];
      return;
    }
    for (let i = start; i < n; i++) {
      current.push(i);
      yield* rec(i + 1);
      current.pop();
    }
  }
  yield* rec(0);
}

/**
 * An ideal in a multivariate polynomial ring.
 */
export class MPolynomialIdeal<R extends CoefficientRing, E extends RingElement> {
  readonly ring: MPolynomialRing<R, E>;
  readonly generators: MPolynomialElement<R, E>[];
  private _groebnerBasis: MPolynomialElement<R, E>[] | null = null;

  constructor(ring: MPolynomialRing<R, E>, generators: MPolynomialElement<R, E>[]) {
    this.ring = ring;
    this.generators = generators.filter((g) => !g.isZero());
  }

  /**
   * Compute and cache the Gröbner basis of this ideal.
   *
   * The zero ideal has Gröbner basis `[0]`, matching SageMath
   * (`multi_polynomial_ideal.py:4586`: `P.ideal([]).groebner_basis()` and
   * `P.ideal([0]).groebner_basis()` are both `[0]`).
   *
   * @see Deviation: Polynomials — Printing, Factor Shape, Term Orders and Base Rings
   */
  groebner_basis(options?: { interreduce?: boolean; maxIterations?: number }): MPolynomialElement<
    R,
    E
  >[] {
    if (this._groebnerBasis === null) {
      const gb = groebner_basis(this.generators, options);
      this._groebnerBasis = gb.length === 0 ? [this.ring.zero()] : gb;
    }
    return this._groebnerBasis;
  }

  /**
   * Check if a polynomial is in this ideal.
   */
  contains(f: MPolynomialElement<R, E>): boolean {
    const gb = this.groebner_basis();
    const r = reduce(f, gb);
    return r.isZero();
  }

  /**
   * Reduce a polynomial modulo this ideal.
   */
  reduce(f: MPolynomialElement<R, E>): MPolynomialElement<R, E> {
    const gb = this.groebner_basis();
    return reduce(f, gb);
  }

  /**
   * Check if this ideal is the zero ideal.
   */
  isZero(): boolean {
    return this.generators.length === 0;
  }

  /**
   * Check if this ideal is the whole ring (contains 1).
   */
  isOne(): boolean {
    const gb = this.groebner_basis();
    for (const g of gb) {
      if (g.isConstant() && !g.isZero()) {
        return true;
      }
    }
    return false;
  }

  /**
   * The dimension of the ring modulo this ideal.
   *
   * Port of SageMath's toy fallback in `MPolynomialIdeal.dimension`
   * (`multi_polynomial_ideal.py:1128-1192`), which follows Chapter 9,
   * Section 1 of Cox, Little and O'Shea's *Ideals, Varieties, and Algorithms*:
   *
   * - the base ring must be a field (`NotImplementedError` otherwise);
   * - a principal ideal is handled by Theorem 3.5.1 of [Ger2008]: `-1` for a
   *   unit, `n - 1` for a nonzero non-unit, `n` for zero;
   * - otherwise, with `M_i` the set of variables occurring in the i-th leading
   *   monomial, find the smallest subset `J` of the variables meeting every
   *   `M_i` and return `n - |J|`.
   *
   * If the ideal is the total ring, the dimension is `-1` by convention.
   *
   * @returns An integer (never `Infinity`)
   */
  dimension(): number {
    const k = this.ring.base_ring as { is_field?: () => boolean };
    if (typeof k.is_field === 'function' && !k.is_field()) {
      throw new NotImplementedError('implemented only over fields');
    }

    const n = this.ring.ngens();

    // multi_polynomial_ideal.py:1133-1141 -- principal ideals.
    // Note that the constructor drops zero generators, so an ideal built from
    // [0] (or []) reaches this branch with zero generators, i.e. g == 0.
    if (this.generators.length <= 1) {
      const g = this.generators[0];
      if (g === undefined || g.isZero()) {
        return n;
      }
      // Over a field the units are exactly the nonzero constants.
      if (g.isConstant()) {
        return -1;
      }
      return n - 1;
    }

    const gb = this.groebner_basis().filter((g) => !g.isZero());

    // "if self.ring().one() in gb: return -1"
    for (const g of gb) {
      if (g.isConstant()) {
        return -1;
      }
    }

    // "compute M_j, denoted by var_lms": the set of variable indices that
    // occur in each leading monomial.
    const varLms: Set<number>[] = [];
    for (const g of gb) {
      const lt = g.leadingTerm();
      if (!lt) continue;
      const [lm] = lt;
      const s = new Set<number>();
      for (let j = 0; j < n; j++) {
        if ((lm[j] ?? 0) > 0) {
          s.add(j);
        }
      }
      varLms.push(s);
    }

    // Enumerate the subsets J of {0, ..., n-1} by increasing size (this is the
    // order Sage's Set.subsets() iterator uses) and stop at the first J that
    // intersects every M_i.
    let minDimension = -1;
    for (let size = 0; size <= n && minDimension === -1; size++) {
      for (const J of subsetsOfSize(n, size)) {
        let intersectsAll = true;
        for (const M of varLms) {
          let meets = false;
          for (const j of J) {
            if (M.has(j)) {
              meets = true;
              break;
            }
          }
          if (!meets) {
            intersectsAll = false;
            break;
          }
        }
        if (intersectsAll) {
          minDimension = J.length;
          break;
        }
      }
    }

    if (minDimension === -1) {
      minDimension = n;
    }

    return n - minDimension;
  }

  toString(): string {
    if (this.generators.length === 0) {
      return 'Ideal (0)';
    }
    return `Ideal (${this.generators.map((g) => g.toString()).join(', ')})`;
  }
}

/**
 * Create an ideal from generators.
 *
 * SageMath spells this `R.ideal(gens)`, so the ring may always be given
 * explicitly; that form also accepts the empty generator list (the zero
 * ideal), which the ring cannot be inferred from.
 *
 * ```
 * sage: P = PolynomialRing(QQ, 't', 0)
 * sage: P.ideal([]).groebner_basis()
 * [0]
 * ```
 *
 * @param ringOrGenerators - The parent ring, or the list of generators
 * @param maybeGenerators - The generators, when the ring is given
 */
export function ideal<R extends CoefficientRing, E extends RingElement>(
  ringOrGenerators: MPolynomialRing<R, E> | MPolynomialElement<R, E>[],
  maybeGenerators?: MPolynomialElement<R, E>[]
): MPolynomialIdeal<R, E> {
  if (Array.isArray(ringOrGenerators)) {
    const generators = ringOrGenerators;
    if (generators.length === 0) {
      throw new ValueError(
        'cannot determine the parent ring of the zero ideal from an empty list of ' +
          'generators; pass the ring as the first argument'
      );
    }
    return new MPolynomialIdeal(
      generators[0]!.parent as unknown as MPolynomialRing<R, E>,
      generators
    );
  }
  return new MPolynomialIdeal(ringOrGenerators, maybeGenerators ?? []);
}
