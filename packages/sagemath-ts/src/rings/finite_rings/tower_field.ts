/**
 * @module sage/rings/finite_rings/tower_field
 * @description Binary tower fields as used in Binius
 *
 * The tower is constructed as:
 * - T0 = GF(2)
 * - T_{i+1} = T_i[X]/(X^2 + α*X + 1) where α is the generator of T_i
 *
 * This gives fields of size 2^(2^i).
 */

import {
  type CoefficientRing,
  Polynomial,
  type RingElement,
} from '../polynomial/polynomial_element.js';
import { PolynomialRing } from '../polynomial/polynomial_ring.js';
import { QuotientRing, QuotientRingElement } from '../polynomial/quotient_ring.js';
import { GF2, GF2Element, type GF2Field } from './gf2.js';

/**
 * Type for tower field elements at any level.
 */
export type TowerFieldElement = GF2Element | QuotientRingElement<TowerFieldElement>;

/**
 * Type for tower fields at any level.
 */
export type TowerField = GF2Field | QuotientRing<TowerFieldElement>;

/**
 * Construct the i-th level of the binary tower field.
 *
 * - Ti(0) = GF(2)
 * - Ti(1) = GF(4)
 * - Ti(2) = GF(16)
 * - Ti(3) = GF(256)
 * - Ti(i) = GF(2^(2^i))
 *
 * @param i - Tower level (0-indexed)
 * @returns The tower field at level i
 *
 * @example
 * ```typescript
 * const T0 = Ti(0);  // GF(2), 2 elements
 * const T1 = Ti(1);  // GF(4), 4 elements
 * const T2 = Ti(2);  // GF(16), 16 elements
 * const T5 = Ti(5);  // GF(2^32), 4294967296 elements
 * ```
 */
export function Ti(i: number): TowerField {
  if (i < 0) {
    throw new Error('Tower level must be non-negative');
  }

  if (i === 0) {
    return GF2;
  }

  // Build iteratively
  let T: TowerField = GF2;

  for (let j = 0; j < i; j++) {
    T = buildNextLevel(T, j);
  }

  return T;
}

/**
 * Build T_{j+1} from T_j.
 *
 * T_{j+1} = T_j[X]/(X^2 + α*X + 1)
 * where α is the generator of T_j.
 */
function buildNextLevel(T: TowerField, level: number): QuotientRing<TowerFieldElement> {
  // Create polynomial ring T[X]
  const varName = `x${level}`;
  const polyRing = new PolynomialRing<TowerFieldElement>(
    T as CoefficientRing<TowerFieldElement>,
    varName
  );

  // Get the generator of T
  const alpha = T.gen();

  // Build modulus: X^2 + α*X + 1
  const x = polyRing.gen();
  const one = polyRing.one();

  // α * X
  const alphaX = new Polynomial<TowerFieldElement>(
    [T.zero() as TowerFieldElement, alpha as TowerFieldElement],
    polyRing
  );

  // X^2 + α*X + 1
  const modulus = x.mul(x).add(alphaX).add(one);

  return new QuotientRing(polyRing, modulus);
}

/**
 * Get the cardinality of Ti(i).
 *
 * |Ti(i)| = 2^(2^i)
 */
export function towerCardinality(i: number): bigint {
  return 2n ** (2n ** BigInt(i));
}

/**
 * List all elements of a tower field.
 * Only practical for small tower levels (i <= 4 or so).
 */
export function listTowerElements(T: TowerField): TowerFieldElement[] {
  const elements: TowerFieldElement[] = [];
  for (const elem of T as Iterable<TowerFieldElement>) {
    elements.push(elem);
  }
  return elements;
}

/**
 * Format a tower field element for display, matching SageMath output.
 */
export function formatTowerElement(elem: TowerFieldElement): string {
  if (elem instanceof GF2Element) {
    return elem.value.toString();
  }

  if (elem instanceof QuotientRingElement) {
    return elem.lift.toString();
  }

  return String(elem);
}

/**
 * Create a specific extension field GF(2^n) using the tower construction.
 *
 * Note: This only works for n that are powers of 2.
 */
export function GF2n(n: number): TowerField {
  // Check if n is a power of 2
  if (n <= 0 || (n & (n - 1)) !== 0) {
    throw new Error(`n must be a power of 2, got ${n}`);
  }

  // Find the tower level
  let level = 0;
  let size = 1;
  while (size < n) {
    level++;
    size *= 2;
  }

  return Ti(level);
}
