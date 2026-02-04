/**
 * @module sage/rings/polynomial/multi_polynomial_ideal
 * @description Ideals in multivariate polynomial rings and Gröbner basis computation
 *
 * Port of: sage/rings/polynomial/multi_polynomial_ideal.py
 *
 * Implements Buchberger's algorithm for computing Gröbner bases.
 */

import { ValueError } from '../../errors.js';
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
 * Reduce polynomial f with respect to a set of polynomials G.
 *
 * Returns the remainder when f is divided by the polynomials in G,
 * using multivariate polynomial division.
 */
export function reduce<R extends CoefficientRing, E extends RingElement>(
  f: MPolynomialElement<R, E>,
  G: MPolynomialElement<R, E>[]
): MPolynomialElement<R, E> {
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
   */
  groebner_basis(options?: { interreduce?: boolean; maxIterations?: number }): MPolynomialElement<
    R,
    E
  >[] {
    if (this._groebnerBasis === null) {
      this._groebnerBasis = groebner_basis(this.generators, options);
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
   * Compute the dimension of the quotient ring R/I.
   * Returns -1 if the ideal is the whole ring, Infinity if infinite dimensional.
   *
   * Note: This is a simplified implementation that only handles some cases.
   *
   * @see Deviation: Multivariate Ideal Dimension Approximation
   */
  dimension(): number {
    if (this.isOne()) {
      return -1;
    }

    const gb = this.groebner_basis();
    if (gb.length === 0) {
      return Number.POSITIVE_INFINITY;
    }

    // For univariate case, dimension is 0 if there's a non-zero polynomial
    if (this.ring.ngens() === 1) {
      return gb.length > 0 ? 0 : Number.POSITIVE_INFINITY;
    }

    // General case: count variables that don't appear as pure powers in leading terms
    // This is a rough approximation
    const n = this.ring.ngens();
    const purePowers = new Set<number>();

    for (const g of gb) {
      const lt = g.leadingTerm();
      if (!lt) continue;
      const [lm] = lt;

      // Check if this is a pure power of a single variable
      let singleVar = -1;
      for (let i = 0; i < lm.length; i++) {
        if ((lm[i] ?? 0) > 0) {
          if (singleVar === -1) {
            singleVar = i;
          } else {
            singleVar = -2; // Multiple variables
            break;
          }
        }
      }

      if (singleVar >= 0) {
        purePowers.add(singleVar);
      }
    }

    // Dimension is roughly n - number of eliminated variables
    const dim = n - purePowers.size;
    return dim >= 0 ? dim : 0;
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
 */
export function ideal<R extends CoefficientRing, E extends RingElement>(
  generators: MPolynomialElement<R, E>[]
): MPolynomialIdeal<R, E> {
  if (generators.length === 0) {
    throw new ValueError('ideal requires at least one generator');
  }
  return new MPolynomialIdeal(generators[0]!.parent, generators);
}
