/**
 * @module sage/rings/number_field/order
 * @description Orders in number fields
 *
 * Port of: sage/rings/number_field/order.py
 * Reference: reference/sage/src/sage/rings/number_field/order.py
 *
 * An order in a number field is a subring that is a free Z-module of full rank.
 * The maximal order is also called the ring of integers.
 */

import { gcd as intGcd, isqrt } from '../../arith/misc.js';
import { NotImplementedError, ValueError } from '../../errors.js';
import { Rational } from '../rational.js';
import type { ClassGroup } from './class_group.js';
import { quadraticClassNumber } from './class_group.js';
import type { NumberField, RationalPolynomial } from './number_field.js';
import { NumberFieldElement as NumberFieldElementImpl } from './number_field_element.js';
import type { NumberFieldElement } from './number_field_element.js';
import type { NumberFieldIdeal } from './number_field_ideal.js';
import type { UnitGroup } from './unit_group.js';

/**
 * Representation of a Z-module as a matrix.
 * Each row represents a basis element expressed in terms of the power basis.
 */
export interface ZModule {
  /** Basis matrix: rows are the coordinates of basis elements in the power basis */
  basisMatrix: Rational[][];
  /** The ambient dimension (degree of the number field) */
  dimension: number;
}

/**
 * An order in a number field (a subring that is a free Z-module of full rank).
 *
 * @see Reference: sage/rings/number_field/order.py:Order
 */
export class Order {
  protected readonly _number_field: NumberField;
  protected readonly _basis: NumberFieldElement[];
  protected _cachedDiscriminant?: bigint;

  constructor(number_field: NumberField, basis?: NumberFieldElement[]) {
    this._number_field = number_field;
    if (basis) {
      this._basis = basis;
    } else {
      // Default to power basis [1, alpha, alpha^2, ..., alpha^{n-1}]
      this._basis = number_field.power_basis();
    }
  }

  /**
   * Return the number field this order belongs to.
   * @see Reference: sage/rings/number_field/order.py:number_field
   */
  number_field(): NumberField {
    return this._number_field;
  }

  /**
   * Return the rank of this order as a Z-module.
   * @see Reference: sage/rings/number_field/order.py:rank
   */
  rank(): number {
    return this._number_field.degree();
  }

  /**
   * Return a Z-basis for this order.
   * @see Reference: sage/rings/number_field/order.py:basis
   */
  basis(): NumberFieldElement[] {
    return [...this._basis];
  }

  /**
   * Return a power basis for this order.
   * For the equation order Z[alpha], this is [1, alpha, alpha^2, ...].
   * @see Reference: sage/rings/number_field/order.py:power_basis
   */
  power_basis(): NumberFieldElement[] {
    return this._number_field.power_basis();
  }

  /**
   * Return the discriminant of this order.
   *
   * The discriminant of an order O with Z-basis {b_1, ..., b_n} is
   * det(Tr(b_i * b_j)).
   *
   * @see Reference: sage/rings/number_field/order.py:discriminant
   */
  discriminant(): bigint {
    if (this._cachedDiscriminant !== undefined) {
      return this._cachedDiscriminant;
    }

    const n = this.rank();
    const basis = this._basis;

    // Compute the trace matrix: T[i][j] = Tr(b_i * b_j)
    const traceMatrix: Rational[][] = [];
    for (let i = 0; i < n; i++) {
      const row: Rational[] = [];
      for (let j = 0; j < n; j++) {
        const product = basis[i]!.mul(basis[j]!);
        row.push(product.trace());
      }
      traceMatrix.push(row);
    }

    // Compute determinant
    const det = matrixDeterminant(traceMatrix);

    // The discriminant should be an integer
    if (det.denominator !== 1n) {
      throw new ValueError('discriminant computation error: non-integral result');
    }

    this._cachedDiscriminant = det.numerator;
    return this._cachedDiscriminant;
  }

  /**
   * Return the index of this order in the maximal order.
   *
   * For the equation order Z[alpha], this is related to the discriminant by:
   * disc(Z[alpha]) = [O_K : Z[alpha]]^2 * disc(O_K)
   *
   * @see Reference: sage/rings/number_field/order.py:index_in_maximal_order
   */
  index_in_maximal_order(): bigint {
    // Compute using disc(O) = [O_K : O]^2 * disc(O_K)
    const orderDisc = this.discriminant();
    const fieldDisc = this._number_field.discriminant();

    if (orderDisc === fieldDisc) {
      return 1n;
    }

    // The index squared equals the ratio of discriminants
    const ratio = orderDisc / fieldDisc;
    if (ratio < 0n) {
      throw new ValueError('invalid discriminant ratio');
    }

    // Take square root
    const index = isqrt(ratio);
    if (index * index !== ratio) {
      throw new ValueError('discriminant ratio is not a perfect square');
    }

    return index;
  }

  /**
   * Check if this is the maximal order.
   * @see Reference: sage/rings/number_field/order.py:is_maximal
   */
  is_maximal(): boolean {
    // Check if discriminant equals field discriminant
    // For now, we can only do this for simple cases
    const orderDisc = this.discriminant();
    const fieldDisc = this._number_field.discriminant();

    return orderDisc === fieldDisc;
  }

  /**
   * Return the conductor of this order as an integer.
   *
   * The conductor is the unique positive integer f such that the discriminant
   * of this order is f^2 times the discriminant of the maximal order.
   * For quadratic orders only.
   *
   * @see Reference: sage/rings/number_field/order.py:conductor
   */
  conductor(): bigint {
    // Currently only implemented for quadratic fields
    if (this._number_field.degree() !== 2) {
      throw new NotImplementedError('conductor is only implemented for quadratic fields');
    }

    const orderDisc = this.discriminant();
    const fieldDisc = this._number_field.discriminant();

    // disc(O) = f^2 * disc(O_K)
    const ratio = orderDisc / fieldDisc;
    if (ratio < 0n) {
      throw new ValueError('invalid discriminant ratio');
    }

    const f = isqrt(ratio);
    if (f * f !== ratio) {
      throw new ValueError('discriminant ratio is not a perfect square');
    }

    return f;
  }

  /**
   * Return the ring of integers (maximal order) containing this order.
   * @see Reference: sage/rings/number_field/order.py:ring_of_integers
   */
  ring_of_integers(): Order {
    // Return the maximal order of the number field
    return this._number_field.ring_of_integers();
  }

  /**
   * Return an ideal of this order with given generators.
   * @see Reference: sage/rings/number_field/order.py:ideal
   */
  ideal(...gens: (bigint | NumberFieldElement)[]): NumberFieldIdeal {
    // Delegate to the number field's ideal method
    return this._number_field.ideal(...gens);
  }

  /**
   * Return the unit group of this order.
   * @see Reference: sage/rings/number_field/order.py:unit_group
   */
  unit_group(): UnitGroup {
    // Delegate to the number field's unit group
    return this._number_field.unit_group();
  }

  /**
   * Return the class group of this order.
   *
   * Currently only implemented for maximal orders.
   *
   * @see Reference: sage/rings/number_field/order.py:class_group
   */
  class_group(): ClassGroup {
    if (this.is_maximal()) {
      return this._number_field.class_group();
    }
    throw new NotImplementedError('class_group: non-maximal orders are not yet supported');
  }

  /**
   * Return the class number of this order.
   *
   * For non-maximal orders in quadratic fields, uses the discriminant
   * to compute the class number via PARI's qfbclassno algorithm.
   *
   * @see Reference: sage/rings/number_field/order.py:class_number
   */
  class_number(): bigint {
    if (this.is_maximal()) {
      return this._number_field.class_number();
    }

    // Non-maximal orders are only supported for quadratic fields
    if (this._number_field.degree() !== 2) {
      throw new NotImplementedError(
        'class_number: non-maximal orders in fields of degree > 2 are not supported'
      );
    }

    // For quadratic orders, use the order discriminant
    const disc = this.discriminant();
    return this._computeQuadraticClassNumber(disc);
  }

  /**
   * Compute the class number of a quadratic order from its discriminant.
   *
   * SageMath calls PARI's `qfbclassno`; this counts classes of binary
   * quadratic forms of the same (possibly non-fundamental) discriminant, which
   * is Gauss's theorem in the same form and is valid for orders as well.
   */
  protected _computeQuadraticClassNumber(disc: bigint): bigint {
    return quadraticClassNumber(disc);
  }

  /**
   * Return the Picard group of this order.
   *
   * The Picard group is the group of invertible fractional ideals
   * modulo principal ideals. For maximal orders, this equals the class group.
   *
   * @see Reference: sage/rings/number_field/order.py:picard_group
   */
  picard_group(): ClassGroup | { order: bigint; structure: bigint[] } {
    if (this.is_maximal()) {
      // For maximal orders, Picard group = class group
      return this.class_group();
    }

    // For non-maximal orders, the Picard group computation is more complex
    // It involves the class group and the conductor
    throw new NotImplementedError(
      'picard_group: non-maximal orders require more complex computation'
    );
  }

  /**
   * Check if an element is in this order.
   *
   * An element is in the order if its coordinates with respect to the
   * order's basis are all integers.
   *
   * @see Reference: sage/rings/number_field/order.py:__contains__
   */
  contains(x: NumberFieldElement): boolean {
    // Check if x can be expressed as a Z-linear combination of the basis
    const coords = this.coordinates(x);
    if (coords === null) return false;

    // All coordinates must be integers
    return coords.every((c) => c.denominator === 1n);
  }

  /**
   * Compute the coordinates of x with respect to this order's basis.
   * Returns null if x cannot be expressed in terms of the basis.
   */
  coordinates(x: NumberFieldElement): Rational[] | null {
    const n = this.rank();
    const xCoeffs = x.list();

    // Build the matrix of basis element coefficients
    const basisMatrix: Rational[][] = this._basis.map((b) => b.list());

    // Solve the system: sum_i c_i * basis[i] = x
    // This is equivalent to: basisMatrix^T * c = xCoeffs

    // Transpose the basis matrix
    const A: Rational[][] = [];
    for (let j = 0; j < n; j++) {
      const row: Rational[] = [];
      for (let i = 0; i < n; i++) {
        row.push(basisMatrix[i]![j]!);
      }
      A.push(row);
    }

    // Solve A * c = xCoeffs using Gaussian elimination
    return solveLinearSystem(A, xCoeffs);
  }

  /**
   * Coerce an element to this order.
   * @see Reference: sage/rings/number_field/order.py:__call__
   */
  __call__(x: unknown): NumberFieldElement {
    // Convert x to a number field element first
    const elem = this._number_field.__call__(x as bigint | Rational | number | NumberFieldElement);

    if (!this.contains(elem)) {
      throw new ValueError('element is not in this order');
    }

    return elem;
  }

  /**
   * Return a random element of this order.
   * @see Reference: sage/rings/number_field/order.py:random_element
   */
  random_element(bound: bigint = 10n): NumberFieldElement {
    // Generate random integer coefficients and form a linear combination of basis
    const basis = this._basis;
    let result = this._number_field.zero();

    for (const b of basis) {
      // Random integer in [-bound, bound]
      const coeff = randomBigInt(-bound, bound);
      result = result.add(b.scalarMul(new Rational(coeff)));
    }

    return result;
  }

  /**
   * Return the characteristic (which is 0).
   * @see Reference: sage/rings/number_field/order.py:characteristic
   */
  characteristic(): bigint {
    return 0n;
  }

  /**
   * Return the degree of this order over Z (same as rank).
   */
  degree(): number {
    return this.rank();
  }

  /**
   * Return the absolute degree (same as rank for absolute orders).
   */
  absolute_degree(): number {
    return this.rank();
  }

  toString(): string {
    if (this.is_maximal()) {
      return `Maximal Order in ${this._number_field}`;
    }
    return `Order in ${this._number_field}`;
  }
}

/**
 * The maximal order (ring of integers) of a number field.
 * @see Reference: sage/rings/number_field/order.py:AbsoluteOrder
 */
export class AbsoluteOrder extends Order {
  private _isMaximal: boolean = true;

  constructor(number_field: NumberField, basis?: NumberFieldElement[]) {
    super(number_field, basis);
  }

  override is_maximal(): boolean {
    return this._isMaximal;
  }

  /**
   * Return the different ideal.
   *
   * The different is the inverse of the codifferent. It measures ramification.
   * For maximal orders, this is the ideal generated by f'(alpha) where
   * f is the minimal polynomial and alpha is the generator.
   *
   * @see Reference: sage/rings/number_field/order.py:different
   */
  different(): NumberFieldIdeal {
    // The different ideal is the ideal generated by f'(alpha)
    // where f is the defining polynomial and alpha is the generator
    const K = this._number_field;
    const gen = K.gen();

    // Get the defining polynomial and its derivative
    const f = K.defining_polynomial();
    const fPrime = f.derivative();

    // Evaluate f'(alpha) in K by reading off its coefficients in the power basis
    const n = K.degree();
    const coeffs = Array.from({ length: n }, (_, i) => fPrime.getCoeff(i));
    const fPrimeAtAlpha = new NumberFieldElementImpl(K, coeffs);

    // The different is the ideal generated by f'(alpha) in the maximal order
    return K.ideal(fPrimeAtAlpha);
  }

  /**
   * Return the codifferent (dual or inverse different).
   *
   * The codifferent is {x in K : Tr(x * O) ⊆ Z}.
   * It is the inverse of the different ideal.
   *
   * @see Reference: sage/rings/number_field/order.py:codifferent
   */
  codifferent(): NumberFieldIdeal {
    // The codifferent is the inverse of the different
    const diff = this.different();
    return diff.inverse();
  }

  /**
   * Return the prime ideals lying above a rational prime p.
   *
   * Uses Dedekind's theorem for factorization in simple extensions:
   * Factor f(x) mod p into irreducibles g_i(x), then the primes above p
   * are (p, g_i(alpha)).
   *
   * @see Reference: sage/rings/number_field/order.py:prime_above
   */
  primes_above(p: bigint): NumberFieldIdeal[] {
    return this._number_field.primes_above(p);
  }

  /**
   * Return one prime ideal of this order lying above the rational prime `p`.
   *
   * @see Reference: sage/rings/number_field/order.py:prime_above
   */
  prime_above(p: bigint): NumberFieldIdeal {
    return this._number_field.prime_above(p);
  }

  /**
   * Decompose a prime p in this order.
   *
   * Returns the factorization of p*O into prime ideals as an array of
   * [prime ideal, exponent] pairs.
   *
   * @see Reference: sage/rings/number_field/order.py:decomposition
   */
  decomposition(p: bigint): Array<[NumberFieldIdeal, bigint]> {
    return this._number_field.decomposition(p);
  }

  override toString(): string {
    return `Maximal Order in ${this._number_field}`;
  }
}

/**
 * A relative order in a number field extension.
 * @see Reference: sage/rings/number_field/order.py:RelativeOrder
 */
export class RelativeOrder extends Order {
  private readonly _base_order: Order;

  constructor(number_field: NumberField, base_order: Order, basis?: NumberFieldElement[]) {
    super(number_field, basis);
    this._base_order = base_order;
  }

  /**
   * Return the base order.
   * @see Reference: sage/rings/number_field/order.py:base_order
   */
  base_order(): Order {
    return this._base_order;
  }

  /**
   * Return the absolute order.
   *
   * For a relative order in a tower K/L/Q, this returns the corresponding
   * order viewed as an order over Z (absolute extension).
   *
   * @see Reference: sage/rings/number_field/order.py:absolute_order
   */
  absolute_order(): AbsoluteOrder {
    // For now, we just return the maximal order of the absolute field
    // A full implementation would need to compute the absolute basis
    const K = this._number_field;
    return new AbsoluteOrder(K, K.power_basis());
  }

  override toString(): string {
    return `Relative Order in ${this._number_field} over ${this._base_order}`;
  }
}

/**
 * Create the equation order Z[alpha] for a number field.
 *
 * The equation order is the order generated by the field's generator.
 * It has basis [1, alpha, alpha^2, ..., alpha^{n-1}].
 */
export function EquationOrder(K: NumberField): Order {
  return new Order(K, K.power_basis());
}

// Helper functions

/**
 * Compute the determinant of a rational matrix using Gaussian elimination.
 */
function matrixDeterminant(m: Rational[][]): Rational {
  const n = m.length;
  if (n === 0) return Rational.one();

  // Make a copy
  const matrix = m.map((row) => [...row]);

  let det = Rational.one();
  let sign = 1;

  for (let col = 0; col < n; col++) {
    // Find pivot
    let pivotRow = -1;
    for (let row = col; row < n; row++) {
      if (!matrix[row]![col]!.isZero()) {
        pivotRow = row;
        break;
      }
    }

    if (pivotRow === -1) {
      return Rational.zero(); // Singular matrix
    }

    // Swap rows if needed
    if (pivotRow !== col) {
      [matrix[col], matrix[pivotRow]] = [matrix[pivotRow]!, matrix[col]!];
      sign = -sign;
    }

    const pivot = matrix[col]![col]!;
    det = det.mul(pivot);

    // Eliminate column below pivot
    for (let row = col + 1; row < n; row++) {
      const factor = matrix[row]![col]!.div(pivot);
      for (let j = col; j < n; j++) {
        matrix[row]![j] = matrix[row]![j]!.sub(factor.mul(matrix[col]![j]!));
      }
    }
  }

  return sign === 1 ? det : det.neg();
}

/**
 * Solve a linear system Ax = b using Gaussian elimination.
 * Returns null if the system has no solution.
 */
function solveLinearSystem(A: Rational[][], b: Rational[]): Rational[] | null {
  const n = A.length;
  if (n === 0) return [];
  if (b.length !== n) return null;

  // Augment matrix with b
  const matrix = A.map((row, i) => [...row, b[i]!]);

  // Gaussian elimination with partial pivoting
  for (let col = 0; col < n; col++) {
    // Find pivot
    let maxRow = col;
    let maxVal = matrix[col]![col]!.abs();

    for (let row = col + 1; row < n; row++) {
      const val = matrix[row]![col]!.abs();
      if (val.gt(maxVal)) {
        maxRow = row;
        maxVal = val;
      }
    }

    if (maxVal.isZero()) {
      return null; // Singular matrix
    }

    // Swap rows
    if (maxRow !== col) {
      [matrix[col], matrix[maxRow]] = [matrix[maxRow]!, matrix[col]!];
    }

    // Eliminate
    const pivot = matrix[col]![col]!;
    for (let row = col + 1; row < n; row++) {
      const factor = matrix[row]![col]!.div(pivot);
      for (let j = col; j <= n; j++) {
        matrix[row]![j] = matrix[row]![j]!.sub(factor.mul(matrix[col]![j]!));
      }
    }
  }

  // Back substitution
  const x: Rational[] = Array(n).fill(Rational.zero());

  for (let i = n - 1; i >= 0; i--) {
    let sum = matrix[i]![n]!;
    for (let j = i + 1; j < n; j++) {
      sum = sum.sub(matrix[i]![j]!.mul(x[j]!));
    }
    x[i] = sum.div(matrix[i]![i]!);
  }

  return x;
}

/**
 * Generate a random bigint in [min, max].
 */
function randomBigInt(min: bigint, max: bigint): bigint {
  const range = max - min + 1n;
  const bits = range.toString(2).length;
  const bytes = Math.ceil(bits / 8);

  // Use crypto.getRandomValues if available, otherwise use Math.random
  let result: bigint;
  if (typeof crypto !== 'undefined' && crypto.getRandomValues) {
    const array = new Uint8Array(bytes);
    crypto.getRandomValues(array);
    result = array.reduce((acc, byte) => (acc << 8n) | BigInt(byte), 0n);
  } else {
    result = BigInt(Math.floor(Math.random() * Number(range)));
  }

  return min + (result % range);
}
