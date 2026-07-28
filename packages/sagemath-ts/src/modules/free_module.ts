/**
 * @module sage/modules/free_module
 * @description Free modules over commutative rings
 * @see Reference: sage/modules/free_module.py
 *
 * Sage supports computation with free modules over an arbitrary commutative ring.
 * Nontrivial functionality is available over ZZ, fields, and some principal
 * ideal domains (e.g. QQ[x] and rings of integers of number fields).
 */

import { ArithmeticError, NotImplementedError, ValueError } from '../errors.js';
import {
  IntegerMatrixFromEntries,
  hermite_normal_form,
  saturation as matrixSaturation,
} from '../matrix/matrix_integer.js';
import { Rational } from '../rings/rational.js';
import {
  type FreeModuleElement,
  FreeModuleElementDense,
  FreeModuleElementSparse,
  type FreeModuleParent,
  type RingLike,
} from './free_module_element.js';

/**
 * Options for creating a free module.
 */
export interface FreeModuleOptions {
  sparse?: boolean;
  innerProductMatrix?: unknown;
  withBasis?: 'standard' | null;
}

/**
 * Create a free module over the given base ring.
 *
 * @param baseRing - A commutative ring
 * @param rank - The rank (dimension) of the module
 * @param options - Additional options
 * @returns A free module of the given rank over the base ring
 *
 * @example
 * ```typescript
 * // Create ZZ^3
 * const M = FreeModule(ZZ, 3);
 *
 * // Create a sparse module
 * const S = FreeModule(ZZ, 100, { sparse: true });
 *
 * // Create with inner product matrix
 * const Q = FreeModule(ZZ, 2, { innerProductMatrix: [[1, 0], [0, -1]] });
 * ```
 *
 * @see Reference: sage/modules/free_module.py:FreeModule
 */
export function FreeModule(
  baseRing: RingLike,
  rank: number | bigint,
  options?: FreeModuleOptions
): FreeModuleGeneric {
  const rankNum = typeof rank === 'bigint' ? Number(rank) : rank;
  if (rankNum < 0) {
    throw new ValueError(`rank (=${rankNum}) must be nonnegative`);
  }

  const sparse = options?.sparse ?? false;
  const innerProductMatrix = options?.innerProductMatrix;

  // Determine what type of module to create based on base ring
  // For now, we use the generic ambient module
  if (isField(baseRing)) {
    return new FreeModuleAmbientField(baseRing, rankNum, sparse, innerProductMatrix);
  }
  if (isPID(baseRing)) {
    return new FreeModuleAmbientPID(baseRing, rankNum, sparse, innerProductMatrix);
  }

  return new FreeModuleAmbient(baseRing, rankNum, sparse, innerProductMatrix);
}

/**
 * Create a vector space over a field.
 *
 * @param field - A field
 * @param dimension - The dimension of the vector space
 * @param options - Additional options
 * @returns A vector space of the given dimension
 *
 * @see Reference: sage/modules/free_module.py:VectorSpace
 */
export function VectorSpace(
  field: RingLike,
  dimension: number | bigint,
  options?: FreeModuleOptions
): FreeModuleField {
  if (!isField(field)) {
    throw new TypeError(`Argument K (= ${field}) must be a field.`);
  }

  const dimensionNum = typeof dimension === 'bigint' ? Number(dimension) : dimension;
  if (dimensionNum < 0) {
    throw new ValueError(`dimension (=${dimensionNum}) must be nonnegative`);
  }

  const sparse = options?.sparse ?? false;
  const innerProductMatrix = options?.innerProductMatrix;

  return new FreeModuleAmbientField(field, dimensionNum, sparse, innerProductMatrix);
}

/**
 * Return the span of the vectors in gens using scalars from baseRing.
 *
 * @param gens - List of vectors or lists of ring elements
 * @param baseRing - Optional base ring (default: inferred from gens)
 * @param options - Additional options
 * @returns The R-span of the generators
 *
 * @example
 * ```typescript
 * // Create span over QQ
 * const V = span([[1, 2, 5], [2, 2, 2]], QQ);
 *
 * // Create span over ZZ
 * const M = span([[1, 2, 3], [4, 5, 6]], ZZ);
 * ```
 *
 * @see Reference: sage/modules/free_module.py:span
 */
export function span(
  gens: unknown[][],
  baseRing?: RingLike,
  options?: { check?: boolean; alreadyEchelonized?: boolean }
): FreeModuleGeneric {
  if (gens.length === 0) {
    if (!baseRing) {
      throw new ValueError('base_ring must be specified for empty generators');
    }
    return FreeModule(baseRing, 0);
  }

  // Infer degree from first generator
  const degree = gens[0]!.length;

  // Infer base ring if not provided
  const ring = baseRing ?? inferRing(gens[0]!);

  // Create the ambient module
  const ambient = FreeModule(ring, degree);

  // Create the span
  const vectors = gens.map((g) => ambient.createElement(g));

  return ambient.span(vectors, ring, options);
}

// ============================================================================
// Helper functions
// ============================================================================

/**
 * Check if a ring is a field.
 */
function isField(ring: RingLike): boolean {
  if ('is_field' in ring && typeof ring.is_field === 'function') {
    return ring.is_field();
  }
  return false;
}

/**
 * Check if a ring is a PID.
 */
function isPID(ring: RingLike): boolean {
  // Check for common PIDs
  if ('is_principal_ideal_domain' in ring) {
    return (ring as { is_principal_ideal_domain: () => boolean }).is_principal_ideal_domain();
  }
  // ZZ is a PID
  if (ring.toString?.() === 'Integer Ring') {
    return true;
  }
  return false;
}

/**
 * Infer a ring from a list of elements.
 */
function inferRing(elements: unknown[]): RingLike {
  if (elements.length === 0) {
    // Default to integers
    return {
      zero: () => 0n,
      one: () => 1n,
      is_field: () => false,
    };
  }

  const first = elements[0];

  // Check for bigint
  if (typeof first === 'bigint') {
    return {
      zero: () => 0n,
      one: () => 1n,
      is_field: () => false,
      toString: () => 'Integer Ring',
    };
  }

  // Check for number (assume rationals)
  if (typeof first === 'number') {
    return {
      zero: () => 0,
      one: () => 1,
      is_field: () => true,
      toString: () => 'Rational Field',
    };
  }

  // Check for objects with parent method
  if (typeof first === 'object' && first !== null && 'parent' in first) {
    const parent = (first as { parent: () => RingLike }).parent();
    return parent;
  }

  // Default to integers
  return {
    zero: () => 0n,
    one: () => 1n,
    is_field: () => false,
  };
}
// ============================================================================
// Exact arithmetic helpers
//
// SageMath performs all of the linear algebra below over the base ring (or its
// fraction field) with exact arithmetic: Hermite normal form over ZZ, reduced
// row echelon form over a field.  These helpers provide the same exactness for
// the loosely typed `RingLike` rings used by this port.
// ============================================================================

/**
 * Exact arithmetic in the fraction field of a base ring.
 *
 * `lift` maps an entry of the module into the fraction field, `lower` maps a
 * fraction field element back to the representation used for entries.
 */
interface FractionFieldArithmetic {
  /** Whether the base ring is its own fraction field. */
  readonly isField: boolean;
  /** Whether the base ring is ZZ (entries are bigints, echelon form is HNF). */
  readonly isIntegral: boolean;
  /** Whether exact arithmetic in the fraction field is available at all. */
  readonly exact: boolean;
  zero(): unknown;
  one(): unknown;
  add(a: unknown, b: unknown): unknown;
  sub(a: unknown, b: unknown): unknown;
  mul(a: unknown, b: unknown): unknown;
  div(a: unknown, b: unknown): unknown;
  neg(a: unknown): unknown;
  isZero(a: unknown): boolean;
  eq(a: unknown, b: unknown): boolean;
  lift(x: unknown): unknown;
  lower(x: unknown): unknown;
  /** Whether a fraction field element belongs to the base ring. */
  inBaseRing(x: unknown): boolean;
  /** Denominator of a fraction field element (1 unless the ring is ZZ-like). */
  denominator(x: unknown): bigint;
}

/**
 * Convert an arbitrary entry to an exact rational.
 */
function toRational(x: unknown): Rational {
  if (x instanceof Rational) {
    return x;
  }
  if (typeof x === 'bigint') {
    return new Rational(x);
  }
  if (typeof x === 'number') {
    return Rational.from(x);
  }
  if (typeof x === 'object' && x !== null) {
    const value = (x as { value?: unknown }).value;
    if (typeof value === 'bigint') {
      return new Rational(value);
    }
    const num = (x as { numerator?: unknown }).numerator;
    const den = (x as { denominator?: unknown }).denominator;
    if (typeof num === 'bigint' && typeof den === 'bigint') {
      return new Rational(num, den);
    }
  }
  return Rational.from(String(x));
}

/**
 * Exact arithmetic over QQ, used for base rings whose elements are bigints
 * (ZZ), JavaScript numbers, or {@link Rational}s.
 */
class RationalArithmetic implements FractionFieldArithmetic {
  readonly isField: boolean;
  readonly isIntegral: boolean;
  readonly exact = true;
  private readonly mode: 'bigint' | 'number' | 'rational';

  constructor(mode: 'bigint' | 'number' | 'rational', isField: boolean) {
    this.mode = mode;
    this.isField = isField;
    this.isIntegral = mode === 'bigint';
  }

  zero(): unknown {
    return Rational.zero();
  }
  one(): unknown {
    return Rational.one();
  }
  add(a: unknown, b: unknown): unknown {
    return (a as Rational).add(b as Rational);
  }
  sub(a: unknown, b: unknown): unknown {
    return (a as Rational).sub(b as Rational);
  }
  mul(a: unknown, b: unknown): unknown {
    return (a as Rational).mul(b as Rational);
  }
  div(a: unknown, b: unknown): unknown {
    return (a as Rational).div(b as Rational);
  }
  neg(a: unknown): unknown {
    return (a as Rational).neg();
  }
  isZero(a: unknown): boolean {
    return (a as Rational).isZero();
  }
  eq(a: unknown, b: unknown): boolean {
    return (a as Rational).eq(b as Rational);
  }
  lift(x: unknown): unknown {
    return toRational(x);
  }
  lower(x: unknown): unknown {
    const r = x as Rational;
    if (this.mode === 'number') {
      return r.toNumber();
    }
    if (this.mode === 'rational') {
      return r;
    }
    return r.isInteger() ? r.numerator : r;
  }
  inBaseRing(x: unknown): boolean {
    if (this.mode === 'bigint') {
      return (x as Rational).isInteger();
    }
    return true;
  }
  denominator(x: unknown): bigint {
    return (x as Rational).denominator;
  }
}

/**
 * Exact arithmetic using the ring elements themselves; used when the base ring
 * is a field whose elements provide `div` (or `inv`).
 */
class RingElementArithmetic implements FractionFieldArithmetic {
  readonly isField = true;
  readonly isIntegral = false;
  readonly exact = true;
  private readonly ring: RingLike;

  constructor(ring: RingLike) {
    this.ring = ring;
  }

  zero(): unknown {
    return this.ring.zero();
  }
  one(): unknown {
    return this.ring.one();
  }
  add(a: unknown, b: unknown): unknown {
    return (a as { add: (x: unknown) => unknown }).add(b);
  }
  sub(a: unknown, b: unknown): unknown {
    return (a as { sub: (x: unknown) => unknown }).sub(b);
  }
  mul(a: unknown, b: unknown): unknown {
    return (a as { mul: (x: unknown) => unknown }).mul(b);
  }
  div(a: unknown, b: unknown): unknown {
    const x = a as { div?: (y: unknown) => unknown; mul: (y: unknown) => unknown };
    if (typeof x.div === 'function') {
      return x.div(b);
    }
    const y = b as { inv?: () => unknown; inverse?: () => unknown };
    if (typeof y.inv === 'function') {
      return x.mul(y.inv());
    }
    if (typeof y.inverse === 'function') {
      return x.mul(y.inverse());
    }
    throw new NotImplementedError('base ring elements do not support division');
  }
  neg(a: unknown): unknown {
    return (a as { neg: () => unknown }).neg();
  }
  isZero(a: unknown): boolean {
    return (a as { isZero: () => boolean }).isZero();
  }
  eq(a: unknown, b: unknown): boolean {
    return (a as { eq: (x: unknown) => boolean }).eq(b);
  }
  lift(x: unknown): unknown {
    if (typeof x === 'object' && x !== null && typeof (x as { add?: unknown }).add === 'function') {
      return x;
    }
    if (this.ring.__call__) {
      return this.ring.__call__(x);
    }
    return x;
  }
  lower(x: unknown): unknown {
    return x;
  }
  inBaseRing(_x: unknown): boolean {
    return true;
  }
  denominator(_x: unknown): bigint {
    return 1n;
  }
}

/**
 * Arithmetic for rings for which no fraction field is available.  Echelon
 * forms and linear solving are not implemented over such rings, exactly as in
 * SageMath, where the corresponding module class stores its generators
 * verbatim (`Submodule_free_ambient`).
 */
class InexactArithmetic implements FractionFieldArithmetic {
  readonly isField = false;
  readonly isIntegral = false;
  readonly exact = false;
  private readonly ring: RingLike;

  constructor(ring: RingLike) {
    this.ring = ring;
  }

  private fail(): never {
    throw new NotImplementedError('exact linear algebra is not implemented over this base ring');
  }

  zero(): unknown {
    return this.ring.zero();
  }
  one(): unknown {
    return this.ring.one();
  }
  add(_a: unknown, _b: unknown): unknown {
    this.fail();
  }
  sub(_a: unknown, _b: unknown): unknown {
    this.fail();
  }
  mul(_a: unknown, _b: unknown): unknown {
    this.fail();
  }
  div(_a: unknown, _b: unknown): unknown {
    this.fail();
  }
  neg(_a: unknown): unknown {
    this.fail();
  }
  isZero(_a: unknown): boolean {
    this.fail();
  }
  eq(_a: unknown, _b: unknown): boolean {
    this.fail();
  }
  lift(_x: unknown): unknown {
    this.fail();
  }
  lower(x: unknown): unknown {
    return x;
  }
  inBaseRing(_x: unknown): boolean {
    return true;
  }
  denominator(_x: unknown): bigint {
    return 1n;
  }
}

/**
 * Return exact fraction field arithmetic for the given base ring.
 */
function arithmeticFor(ring: RingLike): FractionFieldArithmetic {
  let zero: unknown;
  try {
    zero = ring.zero();
  } catch {
    return new InexactArithmetic(ring);
  }

  const isFieldRing = isField(ring);

  if (typeof zero === 'bigint') {
    return new RationalArithmetic('bigint', isFieldRing);
  }
  if (typeof zero === 'number') {
    return new RationalArithmetic('number', isFieldRing);
  }
  if (zero instanceof Rational) {
    return new RationalArithmetic('rational', true);
  }
  if (typeof zero === 'object' && zero !== null) {
    const z = zero as {
      add?: unknown;
      mul?: unknown;
      div?: unknown;
      inv?: unknown;
      value?: unknown;
    };
    if (typeof z.value === 'bigint' && typeof z.add !== 'function') {
      return new RationalArithmetic('bigint', isFieldRing);
    }
    // Only fields get element arithmetic: over a ring that is neither ZZ-like
    // nor a field there is no echelon form, exactly as in SageMath, where the
    // generators of a submodule are then stored verbatim.
    if (isFieldRing && typeof z.add === 'function' && typeof z.mul === 'function') {
      return new RingElementArithmetic(ring);
    }
  }

  return new InexactArithmetic(ring);
}

/**
 * Coerce the entries of a row into the given ring.
 *
 * Fractions are mapped to `num/den` in the target ring when that ring is a
 * field, which is how SageMath coerces a QQ-vector into GF(p) in
 * `change_ring`.  Over a non-field they are left alone: the coordinate ring of
 * the resulting module is then the fraction field of the base ring.
 */
function coerceRow(ring: RingLike, row: unknown[]): unknown[] {
  const ar = arithmeticFor(ring);
  return row.map((e) => {
    if (e instanceof Rational && ring.__call__) {
      if (e.isInteger()) {
        return ring.__call__(e.numerator);
      }
      if (isField(ring)) {
        return ar.div(ring.__call__(e.numerator), ring.__call__(e.denominator));
      }
    }
    return e;
  });
}

/**
 * Compute the GCD of two bigints.
 */
function bigintGcd(a: bigint, b: bigint): bigint {
  a = a < 0n ? -a : a;
  b = b < 0n ? -b : b;
  while (b !== 0n) {
    const t = b;
    b = a % b;
    a = t;
  }
  return a;
}

/**
 * Compute the LCM of two bigints.
 */
function bigintLcm(a: bigint, b: bigint): bigint {
  if (a === 0n || b === 0n) return 0n;
  const g = bigintGcd(a, b);
  const l = (a / g) * b;
  return l < 0n ? -l : l;
}

/**
 * Lift a matrix of entries into the fraction field.
 */
function liftRows(rows: unknown[][], ar: FractionFieldArithmetic): unknown[][] {
  return rows.map((row) => row.map((e) => ar.lift(e)));
}

/**
 * Lower a matrix of fraction field elements back to entries.
 */
function lowerRows(rows: unknown[][], ar: FractionFieldArithmetic): unknown[][] {
  return rows.map((row) => row.map((e) => ar.lower(e)));
}

/**
 * Reduced row echelon form of a lifted matrix over a field.
 *
 * @returns The RREF (without zero rows) and the list of pivot columns
 */
function rrefLifted(
  rows: unknown[][],
  ar: FractionFieldArithmetic
): { rows: unknown[][]; pivots: number[] } {
  const M = rows.map((row) => [...row]);
  const m = M.length;
  const n = m === 0 ? 0 : M[0]!.length;
  const pivots: number[] = [];
  let r = 0;

  for (let col = 0; col < n && r < m; col++) {
    // Find a pivot in this column
    let pivotRow = -1;
    for (let i = r; i < m; i++) {
      if (!ar.isZero(M[i]![col])) {
        pivotRow = i;
        break;
      }
    }
    if (pivotRow === -1) {
      continue;
    }
    if (pivotRow !== r) {
      [M[r], M[pivotRow]] = [M[pivotRow]!, M[r]!];
    }

    // Scale the pivot row so that the pivot is 1
    const pivot = M[r]![col];
    for (let j = col; j < n; j++) {
      M[r]![j] = ar.div(M[r]![j], pivot);
    }

    // Eliminate the column from every other row
    for (let i = 0; i < m; i++) {
      if (i === r) continue;
      const factor = M[i]![col];
      if (ar.isZero(factor)) continue;
      for (let j = col; j < n; j++) {
        M[i]![j] = ar.sub(M[i]![j], ar.mul(factor, M[r]![j]));
      }
    }

    pivots.push(col);
    r++;
  }

  return { rows: M.slice(0, r), pivots };
}

/**
 * Return the echelon form of the given rows over the base ring.
 *
 * Over ZZ this is the Hermite normal form (delegated to
 * `matrix_integer.hermite_normal_form`); over a field it is the reduced row
 * echelon form.  Zero rows are dropped, so the number of rows returned is the
 * rank of the input.
 *
 * @see Reference: sage/modules/free_module.py:FreeModule_submodule_with_basis_pid._echelonized_basis
 */
function echelonRows(rows: unknown[][], ar: FractionFieldArithmetic): unknown[][] {
  if (rows.length === 0 || rows[0]!.length === 0) {
    return rows.map((row) => [...row]);
  }
  if (!ar.exact) {
    return rows.map((row) => [...row]);
  }

  const lifted = liftRows(rows, ar);

  if (ar.isIntegral) {
    // Clear denominators, take the Hermite normal form, restore denominators.
    let d = 1n;
    for (const row of lifted) {
      for (const e of row) {
        d = bigintLcm(d, ar.denominator(e));
      }
    }
    const scaled: bigint[][] = lifted.map((row) =>
      row.map((e) => {
        const r = (e as Rational).mul(new Rational(d));
        if (!r.isInteger()) {
          throw new ArithmeticError('failed to clear denominators of the basis matrix');
        }
        return r.numerator;
      })
    );

    const H = hermite_normal_form(IntegerMatrixFromEntries(scaled), 'default', false, false) as {
      nrows: number;
      ncols: number;
      get: (i: number, j: number) => { value: bigint };
    };

    const out: unknown[][] = [];
    const dr = new Rational(d);
    for (let i = 0; i < H.nrows; i++) {
      const row: unknown[] = [];
      let nonzero = false;
      for (let j = 0; j < H.ncols; j++) {
        const v = H.get(i, j).value;
        if (v !== 0n) nonzero = true;
        row.push(ar.lower(new Rational(v).div(dr)));
      }
      if (nonzero) {
        out.push(row);
      }
    }
    return out;
  }

  const { rows: E } = rrefLifted(lifted, ar);
  return lowerRows(E, ar);
}

/**
 * Return the rank of the given rows over the base ring.
 */
function rankOfRows(rows: unknown[][], ar: FractionFieldArithmetic): number {
  if (rows.length === 0 || rows[0]!.length === 0) {
    return 0;
  }
  if (!ar.exact) {
    return rows.length;
  }
  return rrefLifted(liftRows(rows, ar), ar).rows.length;
}

/**
 * Return a basis of the right kernel `{x : A x = 0}` of the given rows,
 * in echelon form (SageMath's default `basis='echelon'`).
 */
function rightKernelRows(
  rows: unknown[][],
  ncols: number,
  ar: FractionFieldArithmetic
): unknown[][] {
  if (ncols === 0) {
    return [];
  }
  const lifted = rows.length === 0 ? [] : liftRows(rows, ar);
  const { rows: E, pivots } = rrefLifted(lifted, ar);

  const isPivot = new Array<boolean>(ncols).fill(false);
  for (const p of pivots) {
    isPivot[p] = true;
  }

  const kernel: unknown[][] = [];
  for (let free = 0; free < ncols; free++) {
    if (isPivot[free]) continue;
    const v: unknown[] = new Array(ncols).fill(ar.zero());
    v[free] = ar.one();
    for (let r = 0; r < pivots.length; r++) {
      v[pivots[r]!] = ar.neg(E[r]![free]);
    }
    kernel.push(v);
  }

  if (kernel.length === 0) {
    return [];
  }
  return lowerRows(rrefLifted(kernel, ar).rows, ar);
}

/**
 * Solve `x * B = v` exactly, where the rows of `B` are the basis vectors.
 *
 * @returns The (lifted) coefficient vector, or `null` if there is no solution
 */
function solveLeftLifted(
  B: unknown[][],
  v: unknown[],
  ar: FractionFieldArithmetic
): unknown[] | null {
  const n = B.length; // number of unknowns
  const m = v.length; // number of equations
  if (n === 0) {
    return v.every((e) => ar.isZero(e)) ? [] : null;
  }

  // Augmented matrix of the transposed system: row j is
  //   B[0][j] x_0 + ... + B[n-1][j] x_{n-1} = v[j]
  const A: unknown[][] = [];
  for (let j = 0; j < m; j++) {
    const row: unknown[] = [];
    for (let i = 0; i < n; i++) {
      row.push(B[i]![j]);
    }
    row.push(v[j]);
    A.push(row);
  }

  const { rows: E, pivots } = rrefLifted(A, ar);

  // Inconsistent if the augmented column is a pivot
  if (pivots.includes(n)) {
    return null;
  }

  const x: unknown[] = new Array(n).fill(ar.zero());
  for (let r = 0; r < pivots.length; r++) {
    x[pivots[r]!] = E[r]![n];
  }

  // Verify (free variables were set to zero, which is only valid if the
  // resulting vector really is a solution)
  for (let j = 0; j < m; j++) {
    let acc = ar.zero();
    for (let i = 0; i < n; i++) {
      acc = ar.add(acc, ar.mul(x[i], B[i]![j]));
    }
    if (!ar.eq(acc, v[j])) {
      return null;
    }
  }

  return x;
}

/**
 * Turn a basis of a QQ-kernel into a basis of the corresponding ZZ-kernel.
 *
 * The rows are scaled to be integral and then saturated, which is what
 * SageMath's `integer_kernel` computes.
 *
 * @see Reference: sage/matrix/matrix2.pyx:integer_kernel
 */
function integralKernelRows(rows: unknown[][], ar: FractionFieldArithmetic): unknown[][] {
  if (rows.length === 0) {
    return rows;
  }
  const cleared: bigint[][] = [];
  for (const row of rows) {
    const lifted = row.map((e) => ar.lift(e) as Rational);
    let d = 1n;
    for (const e of lifted) {
      d = bigintLcm(d, e.denominator);
    }
    cleared.push(lifted.map((e) => e.mul(new Rational(d)).numerator));
  }

  const S = matrixSaturation(IntegerMatrixFromEntries(cleared));
  const out: unknown[][] = [];
  for (let i = 0; i < S.nrows; i++) {
    const row: unknown[] = [];
    for (let j = 0; j < S.ncols; j++) {
      row.push(ar.lower(new Rational(S.get(i, j).value)));
    }
    out.push(row);
  }
  return out;
}

/**
 * Determinant of a square matrix of lifted entries, computed exactly by
 * Gaussian elimination over the fraction field.
 */
function determinantLifted(M: unknown[][], ar: FractionFieldArithmetic): unknown {
  const n = M.length;
  if (n === 0) {
    return ar.one();
  }
  const A = M.map((row) => [...row]);
  let det = ar.one();

  for (let col = 0; col < n; col++) {
    let pivotRow = -1;
    for (let i = col; i < n; i++) {
      if (!ar.isZero(A[i]![col])) {
        pivotRow = i;
        break;
      }
    }
    if (pivotRow === -1) {
      return ar.zero();
    }
    if (pivotRow !== col) {
      [A[col], A[pivotRow]] = [A[pivotRow]!, A[col]!];
      det = ar.neg(det);
    }
    const pivot = A[col]![col];
    det = ar.mul(det, pivot);
    for (let i = col + 1; i < n; i++) {
      const factor = ar.div(A[i]![col], pivot);
      if (ar.isZero(factor)) continue;
      for (let j = col; j < n; j++) {
        A[i]![j] = ar.sub(A[i]![j], ar.mul(factor, A[col]![j]));
      }
    }
  }

  return det;
}

// ============================================================================
// Base classes
// ============================================================================

/**
 * Base class for modules with elements represented by elements of a free module.
 * @see Reference: sage/modules/free_module.py:Module_free_ambient
 */
export abstract class ModuleFreeAmbient implements FreeModuleParent {
  protected _baseRing: RingLike;
  protected _degree: number;
  protected _sparse: boolean;
  protected _innerProductMatrix: unknown;

  constructor(
    baseRing: RingLike,
    degree: number,
    sparse: boolean = false,
    innerProductMatrix?: unknown
  ) {
    if (degree < 0) {
      throw new ValueError(`degree (=${degree}) must be nonnegative`);
    }
    this._baseRing = baseRing;
    this._degree = degree;
    this._sparse = sparse;
    this._innerProductMatrix = innerProductMatrix;
  }

  /**
   * Return the degree of this free module. This is the dimension of the
   * ambient vector space in which it is embedded.
   */
  degree(): number {
    return this._degree;
  }

  /**
   * Return the base ring of this module.
   */
  baseRing(): RingLike {
    return this._baseRing;
  }

  /**
   * Return whether this module uses sparse representation.
   */
  isSparse(): boolean {
    return this._sparse;
  }

  /**
   * Return whether elements are represented exactly.
   */
  isExact(): boolean {
    if ('is_exact' in this._baseRing && typeof this._baseRing.is_exact === 'function') {
      return this._baseRing.is_exact();
    }
    return true; // Assume exact by default
  }

  /**
   * Return the inner product matrix for this module.
   */
  innerProductMatrix(): unknown {
    return this._innerProductMatrix;
  }

  /**
   * Create an element of this module from entries.
   */
  createElement(entries: unknown[]): FreeModuleElement {
    // Coerce entries to the base ring if needed
    const coercedEntries = entries.map((e) => {
      if (this._baseRing.__call__) {
        return this._baseRing.__call__(e);
      }
      return e;
    });

    if (this._sparse) {
      return new FreeModuleElementSparse(this, coercedEntries);
    }
    return new FreeModuleElementDense(this, coercedEntries);
  }

  /**
   * Return the zero vector in this module.
   */
  zeroVector(): FreeModuleElement {
    const zero = this._baseRing.zero();
    const entries = new Array(this._degree).fill(zero);

    if (this._sparse) {
      return new FreeModuleElementSparse(this, new Map());
    }
    return new FreeModuleElementDense(this, entries);
  }

  /**
   * Return the zero submodule of this module.
   */
  zeroSubmodule(): FreeModuleGeneric {
    return this.span([], this._baseRing);
  }

  /**
   * Return the R-span of gens, where R is the base_ring.
   * @param gens - List of vectors
   * @param baseRing - Optional base ring
   * @param options - Additional options
   */
  span(
    gens: FreeModuleElement[],
    baseRing?: RingLike,
    options?: { check?: boolean; alreadyEchelonized?: boolean }
  ): FreeModuleGeneric {
    const ring = baseRing ?? this._baseRing;
    const self = this as unknown as FreeModuleGeneric;

    // The span lives in the ambient module, not in self
    // (free_module.py:1586: self._submodule_class(self.ambient_module(), ...)).
    const ambient = self.ambientModule ? self.ambientModule() : self;

    if (ring !== this._baseRing) {
      // The base ring changed: re-span in the ambient module over the new ring
      const M = ambient.changeRing(ring);
      return M.span(
        gens.map((g) => M.createElement(g.list())),
        ring,
        options
      );
    }

    const opts = {
      check: options?.check ?? true,
      alreadyEchelonized: options?.alreadyEchelonized ?? false,
    };

    if (isField(ring)) {
      return new FreeModuleSubspace(ambient as FreeModuleField, gens, opts);
    }
    if (arithmeticFor(ring).exact) {
      return new FreeModuleSubmodulePID(ambient, gens, opts);
    }
    // Over a general ring the generators are stored verbatim, exactly as in
    // SageMath's Submodule_free_ambient.
    return new FreeModuleSubmodule(ambient, gens, opts);
  }

  /**
   * Create the R-submodule of the ambient module with given generators.
   * @param gens - List of vectors or a free module
   * @param options - Additional options
   */
  submodule(
    gens: FreeModuleElement[] | FreeModuleGeneric,
    options?: { check?: boolean; alreadyEchelonized?: boolean }
  ): FreeModuleGeneric {
    const list = Array.isArray(gens) ? gens : gens.gens();
    const V = this.span(list, this._baseRing, options);

    if (options?.check ?? true) {
      if (!V.isSubmodule(this)) {
        throw new ArithmeticError(
          `argument gens (= ${list.map((g) => g.toString()).join(', ')}) does not generate a submodule of self`
        );
      }
    }

    return V;
  }

  /**
   * Return the quotient of self by the given submodule.
   * @param sub - A submodule of self
   * @param check - Whether to check that sub is a submodule
   *
   * @see Reference: sage/modules/free_module.py:Module_free_ambient.quotient_module
   */
  quotientModule(sub: FreeModuleGeneric, check: boolean = true): FreeModuleGeneric {
    // Check that sub is a valid submodule
    if (check) {
      if (!this.isSubmodule.call(sub, this as unknown as ModuleFreeAmbient)) {
        throw new ArithmeticError('sub must be a submodule of self');
      }
    }

    // For now, return the quotient structure
    // This creates a quotient module by the submodule
    return new FreeModuleQuotient(this as unknown as FreeModuleGeneric, sub);
  }

  /**
   * Return whether self is a submodule of other.
   * @param other - Another module
   *
   * @see Reference: sage/modules/free_module.py:Module_free_ambient.is_submodule
   */
  isSubmodule(other: ModuleFreeAmbient): boolean {
    if (this === other) {
      return true;
    }

    if (this._baseRing !== other._baseRing) {
      return false;
    }

    if (this._degree !== other._degree) {
      return false;
    }

    const self = this as unknown as FreeModuleGeneric;
    const target = other as unknown as FreeModuleGeneric;

    if (target.rank() < self.rank()) {
      return false;
    }

    // The zero module is always a submodule
    const selfBasis = self.basis();
    if (selfBasis.length === 0) {
      return true;
    }

    const otherBasis = target.basis();
    if (otherBasis.length === 0) {
      return false;
    }

    // Solve  self.basis_matrix() = M * other.basis_matrix()  and require every
    // entry of M to lie in the base ring (free_module.py:2287).
    const ar = arithmeticFor(this._baseRing);
    if (!ar.exact) {
      throw new NotImplementedError(
        'could not determine whether this is a submodule over this base ring'
      );
    }
    const B = liftRows(
      otherBasis.map((b) => b.list()),
      ar
    );

    for (const gen of selfBasis) {
      const v = gen.list().map((e) => ar.lift(e));
      const x = solveLeftLifted(B, v, ar);
      if (x === null) {
        return false;
      }
      for (const c of x) {
        if (!ar.inBaseRing(c)) {
          return false;
        }
      }
    }
    return true;
  }
}

/**
 * Base class for all free modules.
 * @see Reference: sage/modules/free_module.py:FreeModule_generic
 */
export class FreeModuleGeneric extends ModuleFreeAmbient {
  protected _rank: number;
  protected _coordinateRing: RingLike;
  protected _basis: FreeModuleElement[] | null = null;

  constructor(
    baseRing: RingLike,
    rank: number,
    degree: number,
    sparse: boolean = false,
    coordinateRing?: RingLike,
    innerProductMatrix?: unknown
  ) {
    super(baseRing, degree, sparse, innerProductMatrix);
    this._rank = rank;
    this._coordinateRing = coordinateRing ?? baseRing;
  }

  /**
   * Return the rank of this free module.
   */
  rank(): number {
    return this._rank;
  }

  /**
   * Return the dimension of this free module (same as rank).
   */
  dimension(): number {
    return this._rank;
  }

  /**
   * Return the codimension of this free module.
   */
  codimension(): number {
    return this._degree - this._rank;
  }

  /**
   * Return the basis of this module.
   */
  basis(): FreeModuleElement[] {
    if (this._basis !== null) {
      return this._basis;
    }

    // Default: standard basis
    const zero = this._baseRing.zero();
    const one = this._baseRing.one();
    const basisVectors: FreeModuleElement[] = [];

    for (let i = 0; i < this._rank; i++) {
      const entries: unknown[] = new Array(this._degree).fill(zero);
      entries[i] = one;
      const v = this.createElement(entries);
      v.setImmutable();
      basisVectors.push(v);
    }

    this._basis = basisVectors;
    return basisVectors;
  }

  /**
   * Return a tuple of basis elements.
   */
  gens(): FreeModuleElement[] {
    return this.basis();
  }

  /**
   * Return the i-th generator.
   * @param i - The index (default: 0)
   */
  gen(i: number = 0): FreeModuleElement {
    const b = this.basis();
    if (i < 0 || i >= b.length) {
      throw new ValueError(`generator index ${i} out of range [0, ${b.length})`);
    }
    return b[i]!;
  }

  /**
   * Return the number of basis elements.
   */
  ngens(): number {
    return this._rank;
  }

  /**
   * Return the matrix whose rows are the basis for this free module.
   * @param ring - Optional ring for the matrix
   */
  basisMatrix(ring?: RingLike): unknown {
    // For ambient modules, this is the identity matrix
    // For submodules, it's the matrix whose rows are the basis vectors
    const b = this.basis();
    const entries: unknown[][] = [];

    for (const v of b) {
      entries.push(v.list());
    }

    return entries;
  }

  /**
   * Return the echelonized basis matrix.
   *
   * The echelonized basis matrix is the row echelon form of the basis matrix.
   * For ambient modules, this is the identity matrix.
   *
   * @see Reference: sage/modules/free_module.py:FreeModule_generic.echelonized_basis_matrix
   */
  echelonizedBasisMatrix(): unknown[][] {
    if (this.isAmbient()) {
      // For ambient modules, the echelonized basis is the identity
      const n = this._rank;
      const result: unknown[][] = [];
      const zero = this._baseRing.zero();
      const one = this._baseRing.one();

      for (let i = 0; i < n; i++) {
        const row: unknown[] = [];
        for (let j = 0; j < n; j++) {
          row.push(i === j ? one : zero);
        }
        result.push(row);
      }
      return result;
    }

    // For submodules: the Hermite normal form over ZZ, the reduced row
    // echelon form over a field.
    const basisMat = this.basisMatrix() as unknown[][];
    return echelonRows(basisMat, arithmeticFor(this._baseRing));
  }

  /**
   * Return the echelonized basis of this module.
   *
   * @see Reference: sage/modules/free_module.py:FreeModule_submodule_with_basis_pid.echelonized_basis
   */
  echelonizedBasis(): FreeModuleElement[] {
    const ambient = this.ambientModule();
    return this.echelonizedBasisMatrix().map((row) => {
      const v = ambient.createElement(row);
      v.setImmutable();
      return v;
    });
  }

  /**
   * Return whether this module equals other, i.e. whether they have the same
   * ambient space and the same echelonized basis.
   *
   * @see Reference: sage/modules/free_module.py:FreeModule_generic._eq
   */
  equals(other: FreeModuleGeneric): boolean {
    if (this === other) {
      return true;
    }
    if (this._degree !== other.degree() || this._baseRing !== other.baseRing()) {
      return false;
    }
    if (this.rank() !== other.rank()) {
      return false;
    }
    const A = this.echelonizedBasisMatrix();
    const B = other.echelonizedBasisMatrix();
    if (A.length !== B.length) {
      return false;
    }
    const ar = arithmeticFor(this._baseRing);
    for (let i = 0; i < A.length; i++) {
      for (let j = 0; j < this._degree; j++) {
        if (!ar.eq(ar.lift(A[i]![j]), ar.lift(B[i]![j]))) {
          return false;
        }
      }
    }
    return true;
  }

  /**
   * Return the basis matrix of this module.
   */
  matrix(): unknown {
    return this.basisMatrix();
  }

  /**
   * Return the Gram matrix B*A*B^T where A is the inner product matrix
   * and B is the basis matrix. If A is the identity (standard inner product),
   * this is just B*B^T.
   * @returns A 2D array representing the Gram matrix
   */
  gramMatrix(): unknown[][] {
    const zero = this._baseRing.zero();
    const one = this._baseRing.one();

    if (this.isAmbient()) {
      // The Gram matrix of an ambient module is its inner product matrix
      // (the identity when there is none).
      const A = this.innerProductMatrix();
      if (Array.isArray(A)) {
        return (A as unknown[][]).map((row) => [...row]);
      }
      const G: unknown[][] = [];
      for (let i = 0; i < this._degree; i++) {
        const row: unknown[] = [];
        for (let j = 0; j < this._degree; j++) {
          row.push(i === j ? one : zero);
        }
        G.push(row);
      }
      return G;
    }

    // G = B*A*B^t, where A is the inner product matrix of the ambient module
    // and B the basis matrix; the inner product of the basis vectors already
    // applies A.
    const b = this.basis();
    const n = b.length;
    const G: unknown[][] = [];
    for (let i = 0; i < n; i++) {
      const row: unknown[] = [];
      for (let j = 0; j < n; j++) {
        row.push(b[i]!.innerProduct(b[j]!));
      }
      G.push(row);
    }
    return G;
  }

  /**
   * Return the discriminant of this free module.
   *
   * This is the determinant of the Gram matrix.  When the module carries an
   * inner product matrix it is a free quadratic module, whose discriminant
   * carries the extra sign `(-1)^(rank/2)`.
   *
   * @see Reference: sage/modules/free_module.py:FreeModule_generic.discriminant
   * @see Reference: sage/modules/free_quadratic_module.py:FreeQuadraticModule_generic.discriminant
   */
  discriminant(): unknown {
    const G = this.gramMatrix();
    const n = G.length;

    if (n === 0) {
      return this._baseRing.one();
    }

    const ar = arithmeticFor(this._baseRing);
    let det = determinantLifted(liftRows(G, ar), ar);

    if (this._innerProductMatrix !== null && this._innerProductMatrix !== undefined) {
      const r = Math.floor(this.rank() / 2);
      if (r % 2 === 1) {
        det = ar.neg(det);
      }
    }

    return ar.lower(det);
  }

  /**
   * Return the cardinality of this module.
   *
   * @returns The number of elements in this module, or Infinity if infinite.
   *
   * @see Reference: sage/modules/free_module.py:FreeModule_generic.cardinality
   * @see Deviation: finite cardinalities are bigints (they routinely exceed
   *   2^53); the infinite cardinality is the JavaScript Infinity.
   */
  cardinality(): bigint | number {
    if (this._rank === 0) {
      return 1n;
    }

    // Check if base ring is finite
    if ('cardinality' in this._baseRing && typeof this._baseRing.cardinality === 'function') {
      const baseCard = this._baseRing.cardinality();
      if (typeof baseCard === 'bigint') {
        return baseCard ** BigInt(this._rank);
      }
      if (typeof baseCard === 'number') {
        if (!Number.isFinite(baseCard)) {
          return Number.POSITIVE_INFINITY;
        }
        return BigInt(baseCard) ** BigInt(this._rank);
      }
      if (typeof baseCard === 'object' && baseCard !== null) {
        const value = (baseCard as { value?: unknown }).value;
        if (typeof value === 'bigint') {
          return value ** BigInt(this._rank);
        }
      }
    }

    // Check if base ring has is_finite method
    if ('is_finite' in this._baseRing && typeof this._baseRing.is_finite === 'function') {
      if (this._baseRing.is_finite()) {
        // Finite but cardinality not directly available
        throw new NotImplementedError('cardinality not computable for this finite ring');
      }
    }

    // Infinite ring means infinite module (unless rank is 0)
    return Number.POSITIVE_INFINITY;
  }

  /**
   * Return whether this is an ambient module.
   */
  isAmbient(): boolean {
    return false;
  }

  /**
   * Return whether the inner product on this module is the one induced by the
   * ambient inner product.
   *
   * @see Reference: sage/modules/free_module.py:FreeModule_generic.uses_ambient_inner_product
   */
  usesAmbientInnerProduct(): boolean {
    return true;
  }

  /**
   * Return whether this module is dense.
   */
  isDense(): boolean {
    return !this._sparse;
  }

  /**
   * Return whether the rank equals the degree.
   */
  isFull(): boolean {
    return this._rank === this._degree;
  }

  /**
   * Return whether this module is finite.
   */
  isFinite(): boolean {
    if (this._rank === 0) {
      return true;
    }
    // Check if base ring is finite
    if ('is_finite' in this._baseRing) {
      return (this._baseRing as { is_finite: () => boolean }).is_finite();
    }
    return false;
  }

  /**
   * Return whether the given basis has been specified by the user.
   */
  hasUserBasis(): boolean {
    return false;
  }

  /**
   * Return the coordinate ring of this module.
   */
  coordinateRing(): RingLike {
    return this._coordinateRing;
  }

  /**
   * Return the base field (fraction field of the base ring).
   *
   * @see Reference: sage/modules/free_module.py:FreeModule_generic.base_field
   */
  baseField(): RingLike {
    if (isField(this._baseRing)) {
      return this._baseRing;
    }

    // Try to get the fraction field
    if ('fraction_field' in this._baseRing && typeof this._baseRing.fraction_field === 'function') {
      return this._baseRing.fraction_field();
    }

    // For ZZ, the fraction field is QQ
    if (this._baseRing.toString?.() === 'Integer Ring') {
      return {
        zero: () => 0,
        one: () => 1,
        is_field: () => true,
        toString: () => 'Rational Field',
      };
    }

    throw new NotImplementedError('fraction_field not available for this base ring');
  }

  /**
   * Return the ambient module.
   */
  ambientModule(): FreeModuleGeneric {
    return this;
  }

  /**
   * Return the ambient vector space.
   *
   * This is the vector space over the fraction field that contains this module.
   *
   * @see Reference: sage/modules/free_module.py:FreeModule_generic.ambient_vector_space
   */
  ambientVectorSpace(): FreeModuleField {
    const field = this.baseField();
    return new FreeModuleAmbientField(field, this._degree, this._sparse, this._innerProductMatrix);
  }

  /**
   * Write v in terms of the basis for self.
   *
   * The coordinates are computed exactly in the fraction field of the base
   * ring, so they may be rational even when the module is defined over ZZ
   * (SageMath returns them in `FreeModule(R.fraction_field(), rank)`).
   *
   * @param v - A vector
   * @param check - Whether to verify v is in self
   * @returns The list of coefficients c such that v = sum(c[i] * basis[i])
   * @throws {ArithmeticError} If v is not in the span of the basis
   *
   * @see Reference: sage/modules/free_module.py:FreeModule_generic.coordinates
   * @see Deviation: an integral coordinate is returned as a bigint and a
   *   non-integral one as a Rational; `check: false` still raises when v is
   *   outside the span, where SageMath returns a meaningless vector.
   */
  coordinates(v: FreeModuleElement, check: boolean = true): unknown[] {
    const basis = this.basis();
    const n = basis.length;

    if (v.degree() !== this._degree) {
      throw new ArithmeticError('vector is not in free module');
    }

    if (n === 0) {
      if (v.isZero()) {
        return [];
      }
      throw new ArithmeticError('vector is not in free module');
    }

    // For ambient modules the coordinates are just the entries
    if (this.isAmbient() && this._rank === this._degree) {
      return v.list();
    }

    // Solve  x * B = v  exactly over the fraction field of the base ring.
    const ar = arithmeticFor(this._baseRing);
    if (!ar.exact) {
      throw new NotImplementedError('coordinates are not implemented over this base ring');
    }
    const B = liftRows(
      basis.map((b) => b.list()),
      ar
    );
    const target = v.list().map((e) => ar.lift(e));

    // The exact solve fails exactly when v is not in the span of the basis,
    // which is the condition SageMath's `check` verifies.  We therefore always
    // raise, even when check is false.
    const x = solveLeftLifted(B, target, ar);
    if (x === null) {
      throw new ArithmeticError('vector is not in free module');
    }

    return x.map((c) => ar.lower(c));
  }

  /**
   * Return the coordinate vector for v with respect to the basis.
   * @param v - A vector
   * @param check - Whether to verify v is in self
   */
  coordinateVector(v: FreeModuleElement, check?: boolean): FreeModuleElement {
    const coords = this.coordinates(v, check);

    // Create a new free module for the coordinate space
    const coordModule = FreeModule(this._baseRing, coords.length);
    return coordModule.createElement(coords);
  }

  /**
   * Return the direct sum of self and other.
   * @param other - Another free module
   * @returns A new module of rank self.rank() + other.rank()
   */
  directSum(other: FreeModuleGeneric): FreeModuleGeneric {
    // Direct sum has rank = rank(self) + rank(other)
    // and degree = degree(self) + degree(other) (for embedded direct sum)
    // The basis vectors are (v, 0) for v in basis(self) and (0, w) for w in basis(other)

    const newRank = this._rank + other.rank();
    const newDegree = this._degree + other.degree();

    const newModule = FreeModule(this._baseRing, newDegree, { sparse: this._sparse });

    // Build basis vectors
    const zero = this._baseRing.zero();
    const basisVectors: FreeModuleElement[] = [];

    // Basis vectors from self: (v, 0)
    for (const v of this.basis()) {
      const entries: unknown[] = [];
      for (let i = 0; i < this._degree; i++) {
        entries.push(v.getItem(i));
      }
      for (let i = 0; i < other.degree(); i++) {
        entries.push(zero);
      }
      basisVectors.push(newModule.createElement(entries));
    }

    // Basis vectors from other: (0, w)
    for (const w of other.basis()) {
      const entries: unknown[] = [];
      for (let i = 0; i < this._degree; i++) {
        entries.push(zero);
      }
      for (let i = 0; i < other.degree(); i++) {
        entries.push(w.getItem(i));
      }
      basisVectors.push(newModule.createElement(entries));
    }

    return newModule.span(basisVectors);
  }

  /**
   * Return the product of this module by a scalar.
   * If M is a module with basis b_1, ..., b_n, then scale(c) has basis c*b_1, ..., c*b_n.
   * @param scalar - A scalar from the base ring
   * @returns A new module with scaled basis
   */
  scale(scalar: unknown): FreeModuleGeneric {
    const b = this.basis();
    const scaledVectors: FreeModuleElement[] = [];

    for (const v of b) {
      scaledVectors.push(v.mul(scalar));
    }

    return this.span(scaledVectors);
  }

  /**
   * Return a random element of this module.
   * @param prob - Probability each coefficient is non-zero (default: 1.0)
   * @param min - Minimum value for random coefficients (default: -10)
   * @param max - Maximum value for random coefficients (default: 10)
   */
  randomElement(prob: number = 1.0, min: number = -10, max: number = 10): FreeModuleElement {
    const entries: unknown[] = [];
    const zero = this._baseRing.zero();

    for (let i = 0; i < this._degree; i++) {
      if (Math.random() < prob) {
        // Generate random value
        if (typeof zero === 'bigint') {
          const val = BigInt(Math.floor(Math.random() * (max - min + 1)) + min);
          entries.push(val);
        } else if (typeof zero === 'number') {
          entries.push(Math.floor(Math.random() * (max - min + 1)) + min);
        } else if (this._baseRing.__call__) {
          const val = Math.floor(Math.random() * (max - min + 1)) + min;
          entries.push(this._baseRing.__call__(val));
        } else {
          entries.push(zero);
        }
      } else {
        entries.push(zero);
      }
    }

    return this.createElement(entries);
  }

  /**
   * Return whether the given vectors are linearly dependent.
   * @param vecs - A list of vectors
   */
  areLinearlyDependent(vecs: FreeModuleElement[]): boolean {
    if (vecs.length === 0) {
      return false; // Empty set is linearly independent
    }

    if (vecs.length > this._degree) {
      return true; // More vectors than the degree means dependent
    }

    // A = matrix(vecs); A.echelonize(); any zero row means dependence.
    const ar = arithmeticFor(this._baseRing);
    if (!ar.exact) {
      throw new NotImplementedError('linear dependence is not implemented over this base ring');
    }
    return (
      rankOfRows(
        vecs.map((v) => v.list()),
        ar
      ) < vecs.length
    );
  }

  /**
   * Return a linear combination of the basis vectors.
   * @param coefficients - The coefficients
   */
  linearCombinationOfBasis(coefficients: unknown[]): FreeModuleElement {
    const b = this.basis();
    if (coefficients.length !== b.length) {
      throw new ValueError(`coefficients must have length ${b.length}, got ${coefficients.length}`);
    }

    let result = this.zeroVector();
    for (let i = 0; i < b.length; i++) {
      result = result.add(b[i]!.mul(coefficients[i]));
    }
    return result;
  }

  /**
   * Return the dense version of this module.
   *
   * @see Reference: sage/modules/free_module.py:FreeModule_generic.dense_module
   */
  denseModule(): FreeModuleGeneric {
    if (!this._sparse) {
      return this;
    }

    // Create a dense ambient module and span the same vectors
    const denseAmbient = FreeModule(this._baseRing, this._degree, { sparse: false });
    const basis = this.basis();
    const denseBasis: FreeModuleElement[] = [];

    for (const v of basis) {
      denseBasis.push(denseAmbient.createElement(v.list()));
    }

    return denseAmbient.span(denseBasis);
  }

  /**
   * Return the sparse version of this module.
   *
   * @see Reference: sage/modules/free_module.py:FreeModule_generic.sparse_module
   */
  sparseModule(): FreeModuleGeneric {
    if (this._sparse) {
      return this;
    }

    // Create a sparse ambient module and span the same vectors
    const sparseAmbient = FreeModule(this._baseRing, this._degree, { sparse: true });
    const basis = this.basis();
    const sparseBasis: FreeModuleElement[] = [];

    for (const v of basis) {
      sparseBasis.push(sparseAmbient.createElement(v.list()));
    }

    return sparseAmbient.span(sparseBasis);
  }

  /**
   * Return the free module over `ring` obtained by coercing each element of
   * the basis of self into a vector over the fraction field of `ring`, then
   * taking the resulting module.
   *
   * @param ring - A principal ideal domain
   * @throws {TypeError} If the new ring is not a principal ideal domain
   *
   * @see Reference: sage/modules/free_module.py:FreeModule_submodule_with_basis_pid.change_ring
   * @see Reference: sage/modules/free_module.py:FreeModule_ambient.change_ring
   */
  changeRing(ring: RingLike): FreeModuleGeneric {
    if (ring === this._baseRing) {
      return this;
    }
    if (!isField(ring) && !isPID(ring)) {
      throw new TypeError(
        `the new ring ${ring.toString?.() ?? ring} should be a principal ideal domain`
      );
    }

    if (this.isAmbient()) {
      return FreeModule(ring, this._rank, { sparse: this._sparse });
    }

    // Re-span the basis, in the ambient module of the same degree, over R.
    const M = this.ambientModule().changeRing(ring);
    const B = this.basis().map((b) => M.createElement(coerceRow(ring, b.list())));
    if (this.hasUserBasis() && M instanceof FreeModulePID) {
      return M.spanOfBasis(B, ring);
    }
    return M.span(B, ring);
  }

  /**
   * Return an ambient free module isomorphic to this one.
   */
  nonembeddedFreeModule(): FreeModuleGeneric {
    return FreeModule(this._baseRing, this._rank, { sparse: this._sparse });
  }

  /**
   * String representation.
   */
  toString(): string {
    const ringName = this._baseRing.toString?.() ?? 'Ring';
    if (this._degree === this._rank) {
      if (isField(this._baseRing)) {
        return `Vector space of dimension ${this._rank} over ${ringName}`;
      }
      return `Free module of rank ${this._rank} over ${ringName}`;
    }
    if (isField(this._baseRing)) {
      return `Vector space of degree ${this._degree} and dimension ${this._rank} over ${ringName}`;
    }
    return `Free module of degree ${this._degree} and rank ${this._rank} over ${ringName}`;
  }
}

/**
 * Base class for free modules over an integral domain.
 * @see Reference: sage/modules/free_module.py:FreeModule_generic_domain
 */
export class FreeModuleDomain extends FreeModuleGeneric {
  constructor(
    baseRing: RingLike,
    rank: number,
    degree: number,
    sparse: boolean = false,
    coordinateRing?: RingLike,
    innerProductMatrix?: unknown
  ) {
    super(baseRing, rank, degree, sparse, coordinateRing, innerProductMatrix);
  }

  /**
   * Return the sum of self and other.
   * @param other - Another submodule
   */
  add(other: FreeModuleGeneric): FreeModuleGeneric {
    // Combine generators from both modules
    const gens = [...this.gens(), ...other.gens()];
    return this.span(gens);
  }
}

/**
 * Base class for free modules over a PID.
 * @see Reference: sage/modules/free_module.py:FreeModule_generic_pid
 */
export class FreeModulePID extends FreeModuleDomain {
  constructor(
    baseRing: RingLike,
    rank: number,
    degree: number,
    sparse: boolean = false,
    coordinateRing?: RingLike,
    innerProductMatrix?: unknown
  ) {
    super(baseRing, rank, degree, sparse, coordinateRing, innerProductMatrix);
  }

  /**
   * Return the lattice index [other:self] of self in other.
   *
   * When self is contained in other, the lattice index is the usual index.
   * If the index is infinite, this returns Infinity.
   *
   * @param other - Another module
   *
   * @see Reference: sage/modules/free_module.py:FreeModule_generic_pid.index_in
   * @see Deviation: an integral index is returned as a bigint and a
   *   non-integral one as a Rational.
   */
  indexIn(other: FreeModuleGeneric): unknown {
    if (this._baseRing !== other.baseRing()) {
      throw new NotImplementedError(
        'lattice index only defined for modules over the same base ring.'
      );
    }
    if (this._degree !== other.degree()) {
      throw new ArithmeticError('self and other must be embedded in the same ambient space.');
    }

    const ar = arithmeticFor(this._baseRing);

    if (ar.isField) {
      if (this.equals(other)) {
        return 1n;
      }
      if (this.isSubmodule(other as unknown as ModuleFreeAmbient)) {
        return Number.POSITIVE_INFINITY;
      }
      throw new ArithmeticError('self must be contained in the vector space spanned by other.');
    }

    // C = [other.coordinates(b) for b in self.basis()]
    const C: unknown[][] = [];
    for (const b of this.basis()) {
      C.push(other.coordinates(b).map((c) => ar.lift(c)));
    }

    if (this.rank() < other.rank()) {
      return Number.POSITIVE_INFINITY;
    }

    const det = determinantLifted(C, ar);
    // For ZZ the index is the absolute value of the determinant
    const r = det as Rational;
    const abs = ar.isIntegral && r instanceof Rational ? r.abs() : det;
    return ar.lower(abs);
  }

  /**
   * Return the intersection of self and other.
   *
   * @param other - Another module
   *
   * @see Reference: sage/modules/free_module.py:FreeModule_generic_pid.intersection
   */
  intersection(other: FreeModuleGeneric): FreeModuleGeneric {
    if (this._degree !== other.degree()) {
      throw new ArithmeticError('self and other must be embedded in the same ambient space.');
    }

    // Dispense with the easy cases
    if (this.rank() === 0 || other.rank() === 0) {
      return this.zeroSubmodule();
    }
    if (this === other) {
      return this;
    }
    if (other.isSubmodule(this as unknown as ModuleFreeAmbient)) {
      return other;
    }
    if (this.isSubmodule(other as unknown as ModuleFreeAmbient)) {
      return this;
    }

    // Standard algorithm: let S be A1 stacked on A2; the vectors v in the
    // (left) kernel of S give the intersection as (v[:n]) * A1.
    const [V1, V2] =
      this.rank() <= other.rank()
        ? [this as unknown as FreeModuleGeneric, other]
        : [other, this as unknown as FreeModuleGeneric];

    const ar = arithmeticFor(this._baseRing);
    const A1 = V1.basis().map((v) => v.list());
    const A2 = V2.basis().map((v) => v.list());
    const S = [...A1, ...A2];
    const n = A1.length;

    // Left kernel of S = right kernel of S^t
    const St: unknown[][] = [];
    for (let j = 0; j < this._degree; j++) {
      St.push(S.map((row) => row[j]));
    }
    let K = rightKernelRows(St, S.length, ar);

    if (!ar.isField && ar.isIntegral && K.length > 0) {
      // integer_kernel: clear denominators and saturate, so that the kernel is
      // the full ZZ-module of integral relations.
      K = integralKernelRows(K, ar);
    }

    const gens: FreeModuleElement[] = [];
    const ambient = this.ambientModule();
    for (const v of K) {
      const coeffs = v.slice(0, n).map((e) => ar.lift(e));
      const entries: unknown[] = [];
      for (let j = 0; j < this._degree; j++) {
        let acc = ar.zero();
        for (let i = 0; i < n; i++) {
          acc = ar.add(acc, ar.mul(coeffs[i], ar.lift(A1[i]![j])));
        }
        entries.push(ar.lower(acc));
      }
      const w = ambient.createElement(entries);
      if (!w.isZero()) {
        gens.push(w);
      }
    }

    if (gens.length === 0) {
      return this.zeroSubmodule();
    }
    return this.span(gens);
  }

  /**
   * Return the index of this module in its saturation.
   *
   * The saturation of a submodule M of a free module F is the largest
   * submodule of F containing M with the same rank.
   *
   * @see Reference: sage/modules/free_module.py:FreeModule_generic_pid.index_in_saturation
   */
  indexInSaturation(): unknown {
    return this.indexIn(this.saturation());
  }

  /**
   * Return the saturated submodule of R^n that spans the same vector space.
   *
   * The saturation of a submodule M of a free module F is the largest
   * submodule S of F containing M with the same rank. Equivalently,
   * S is the intersection of the vector space span of M with F.
   *
   * For a lattice L in ZZ^n, the saturation is obtained by computing
   * the Hermite normal form of the basis matrix and dividing by the GCD
   * of the entries.
   *
   * @returns The saturation of this module
   *
   * @example
   * ```typescript
   * // Create a non-saturated lattice
   * const L = span([[9, 9, 6]], ZZ);
   * // L.saturation() returns span([[3, 3, 2]], ZZ)
   * ```
   *
   * @see Reference: sage/modules/free_module.py:FreeModule_generic_pid.saturation
   */
  saturation(): FreeModuleGeneric {
    // If the base ring is a field, the module is already saturated
    if (isField(this._baseRing)) {
      return this;
    }
    if (this._rank === 0) {
      return this;
    }

    const basisMat = this.basisMatrix() as unknown[][];
    if (basisMat.length === 0) {
      return this;
    }

    const ar = arithmeticFor(this._baseRing);
    if (!ar.isIntegral) {
      throw new NotImplementedError('saturation is only implemented over ZZ');
    }

    // A, _ = self.basis_matrix()._clear_denom(); S = self.span(A.saturation())
    const lifted = liftRows(basisMat, ar);
    let d = 1n;
    for (const row of lifted) {
      for (const e of row) {
        d = bigintLcm(d, ar.denominator(e));
      }
    }
    const cleared: bigint[][] = lifted.map((row) =>
      row.map((e) => (e as Rational).mul(new Rational(d)).numerator)
    );

    const S = matrixSaturation(IntegerMatrixFromEntries(cleared));

    const ambient = this.ambientModule();
    const saturatedVectors: FreeModuleElement[] = [];
    for (let i = 0; i < S.nrows; i++) {
      const row: unknown[] = [];
      for (let j = 0; j < S.ncols; j++) {
        row.push(ar.lower(new Rational(S.get(i, j).value)));
      }
      saturatedVectors.push(ambient.createElement(row));
    }

    const sat = ambient.span(saturatedVectors);

    // Return exactly self if it is already saturated
    return this.equals(sat) ? this : sat;
  }

  /**
   * Return the denominator of the basis matrix.
   *
   * This is the LCM of the denominators of all entries in the basis matrix
   * when expressed in the ambient space coordinates.
   *
   * @returns The denominator (LCM of all entry denominators)
   *
   * @see Reference: sage/modules/free_module.py:FreeModule_generic_pid.denominator
   */
  denominator(): unknown {
    const basisMat = this.basisMatrix() as unknown[][];

    if (basisMat.length === 0) {
      return this._baseRing.one();
    }

    const ar = arithmeticFor(this._baseRing);
    if (!ar.exact) {
      return this._baseRing.one();
    }

    let d = 1n;
    for (const row of basisMat) {
      for (const entry of row) {
        d = bigintLcm(d, ar.denominator(ar.lift(entry)));
      }
    }

    return ar.lower(new Rational(d));
  }

  /**
   * Return the free R-module with the given basis.
   * @param basis - A list of vectors
   * @param baseRing - Optional base ring
   * @param options - Additional options
   */
  spanOfBasis(
    basis: FreeModuleElement[],
    baseRing?: RingLike,
    options?: { check?: boolean; alreadyEchelonized?: boolean }
  ): FreeModuleGeneric {
    const ring = baseRing ?? this._baseRing;
    const ambient = this.ambientModule();

    if (ring !== this._baseRing) {
      const M = ambient.changeRing(ring);
      return (M as FreeModulePID).spanOfBasis(
        basis.map((b) => M.createElement(b.list())),
        ring,
        options
      );
    }

    const opts = {
      check: options?.check ?? true,
      echelonize: false,
      alreadyEchelonized: options?.alreadyEchelonized ?? false,
    };

    if (isField(ring)) {
      return new FreeModuleSubspaceWithBasis(ambient as FreeModuleField, basis, opts);
    }
    return new FreeModuleWithBasis(ambient, basis, opts);
  }

  /**
   * Create the R-submodule with given basis.
   * @param basis - A list of linearly independent vectors
   * @param options - Additional options
   */
  submoduleWithBasis(
    basis: FreeModuleElement[],
    options?: { check?: boolean; alreadyEchelonized?: boolean }
  ): FreeModuleGeneric {
    return this.spanOfBasis(basis, this._baseRing, options);
  }

  /**
   * Create a vector subspace of the ambient vector space.
   *
   * This creates a vector space over the fraction field of the base ring
   * that is spanned by the given generators.
   *
   * @param gens - A list of vectors
   * @param check - Whether to check vectors are in ambient space
   * @returns A vector space over the fraction field
   *
   * @see Reference: sage/modules/free_module.py:FreeModule_generic_pid.vector_space_span
   */
  vectorSpaceSpan(gens: FreeModuleElement[], check?: boolean): FreeModuleField {
    const field = this.baseField();
    const ambient = new FreeModuleAmbientField(
      field,
      this._degree,
      this._sparse,
      this._innerProductMatrix
    );

    if (gens.length === 0) {
      return ambient.subspace([]);
    }

    // Convert generators to vectors in the ambient field
    const fieldGens: FreeModuleElement[] = [];
    for (const gen of gens) {
      const entries = gen.list();
      fieldGens.push(ambient.createElement(entries));
    }

    return ambient.subspace(fieldGens, { check: check ?? true });
  }

  /**
   * Create a vector subspace with given basis.
   *
   * This creates a vector space over the fraction field of the base ring
   * with the given vectors as its basis.
   *
   * @param gens - A list of vectors (must be linearly independent)
   * @param check - Whether to check vectors are in ambient space
   * @returns A vector space with the given basis
   *
   * @see Reference: sage/modules/free_module.py:FreeModule_generic_pid.vector_space_span_of_basis
   */
  vectorSpaceSpanOfBasis(gens: FreeModuleElement[], check?: boolean): FreeModuleField {
    const field = this.baseField();
    const ambient = new FreeModuleAmbientField(
      field,
      this._degree,
      this._sparse,
      this._innerProductMatrix
    );

    if (gens.length === 0) {
      return ambient.subspaceWithBasis([]);
    }

    // Convert generators to vectors in the ambient field
    const fieldGens: FreeModuleElement[] = [];
    for (const gen of gens) {
      const entries = gen.list();
      fieldGens.push(ambient.createElement(entries));
    }

    return ambient.subspaceWithBasis(fieldGens, { check: check ?? true });
  }

  /**
   * Return the tensor product self tensor other over R.
   *
   * The tensor product of free modules of ranks m and n has rank m*n.
   *
   * @param other - Another free module
   * @returns The tensor product module
   *
   * @see Reference: sage/modules/free_module.py:FreeModule_generic_pid.tensor_product
   */
  tensorProduct(other: FreeModuleGeneric): FreeModuleGeneric {
    // Tensor product of R^m and R^n is R^(m*n)
    const newRank = this._rank * other.rank();
    const newModule = FreeModule(this._baseRing, newRank, { sparse: this._sparse });

    // For now, return the ambient module
    // A full implementation would track the basis elements as pairs (e_i, f_j)
    return newModule;
  }
}

/**
 * Base class for vector spaces (free modules over a field).
 * @see Reference: sage/modules/free_module.py:FreeModule_generic_field
 */
export class FreeModuleField extends FreeModulePID {
  constructor(
    baseRing: RingLike,
    rank: number,
    degree: number,
    sparse: boolean = false,
    innerProductMatrix?: unknown
  ) {
    super(baseRing, rank, degree, sparse, undefined, innerProductMatrix);
  }

  /**
   * Return the vector space of which this is a subspace.
   */
  vectorSpace(): FreeModuleField {
    return this;
  }

  /**
   * Return the subspace spanned by gens.
   * @param gens - A list of vectors
   * @param options - Additional options
   */
  subspace(
    gens: FreeModuleElement[],
    options?: { check?: boolean; alreadyEchelonized?: boolean }
  ): FreeModuleField {
    return new FreeModuleSubspace(this.ambientModule() as FreeModuleField, gens, {
      check: options?.check ?? true,
      alreadyEchelonized: options?.alreadyEchelonized ?? false,
    });
  }

  /**
   * Return the subspace with given basis.
   * @param gens - A list of linearly independent vectors
   * @param options - Additional options
   */
  subspaceWithBasis(
    gens: FreeModuleElement[],
    options?: { check?: boolean; alreadyEchelonized?: boolean }
  ): FreeModuleField {
    return new FreeModuleSubspaceWithBasis(this.ambientModule() as FreeModuleField, gens, {
      check: options?.check ?? true,
      echelonize: false,
      alreadyEchelonized: options?.alreadyEchelonized ?? false,
    });
  }

  /**
   * Return whether self is a subspace of other.
   *
   * @param other - Another vector space
   * @returns true if self is a subspace of other
   *
   * @see Reference: sage/modules/free_module.py:FreeModule_generic_field.is_subspace
   */
  isSubspace(other: FreeModuleField): boolean {
    // Use the underlying isSubmodule method
    return this.isSubmodule(other as unknown as ModuleFreeAmbient);
  }

  /**
   * Return the sum of self and other.
   * @param other - Another subspace
   */
  override add(other: FreeModuleField): FreeModuleField {
    const gens = [...this.gens(), ...other.gens()];
    return this.subspace(gens);
  }

  /**
   * Return the intersection of self and other.
   *
   * For vector spaces over a field, this uses the kernel-based algorithm.
   *
   * @param other - Another subspace
   * @returns The intersection of the two subspaces
   *
   * @see Reference: sage/modules/free_module.py:FreeModule_generic_field.intersection
   */
  override intersection(other: FreeModuleField): FreeModuleField {
    // Handle trivial cases
    if (this.dimension() === 0 || other.dimension() === 0) {
      return this.subspace([]);
    }

    if (this === other) {
      return this;
    }

    // If one is contained in the other, return the smaller one
    if (this.isSubspace(other)) {
      return this;
    }
    if (other.isSubspace(this)) {
      return other;
    }

    // Use the parent class intersection method (from FreeModulePID)
    const pidResult = super.intersection(other as unknown as FreeModuleGeneric);

    // Convert result to field subspace
    const gens = pidResult.gens();
    if (gens.length === 0) {
      return this.subspace([]);
    }

    // Create vectors in the ambient field space
    const ambient = this.ambientVectorSpace();
    const fieldGens: FreeModuleElement[] = [];

    for (const gen of gens) {
      fieldGens.push(ambient.createElement(gen.list()));
    }

    return ambient.subspace(fieldGens);
  }

  /**
   * Return the orthogonal complement of this subspace.
   *
   * The orthogonal complement is the set of all vectors v such that
   * <v, w> = 0 for all w in self.
   *
   * @returns The orthogonal complement of this subspace
   *
   * @see Reference: sage/modules/free_module.py:FreeModule_generic_field.orthogonal_complement
   */
  orthogonalComplement(): FreeModuleField {
    const basisMat = this.basisMatrix() as unknown[][];
    const n = this._degree;

    if (basisMat.length === 0) {
      // The complement of the zero space is the whole ambient space
      return this.ambientVectorSpace();
    }
    if (this.dimension() === n) {
      // The complement of the ambient space is zero
      return this.subspace([]);
    }

    // basis_matrix().right_kernel()
    const ar = arithmeticFor(this._baseRing);
    const kernel = rightKernelRows(basisMat, n, ar);

    const ambient = this.ambientVectorSpace();
    if (kernel.length === 0) {
      return ambient.subspace([]);
    }

    return ambient.subspace(kernel.map((row) => ambient.createElement(row)));
  }

  /**
   * Return the quotient self/other.
   *
   * Returns a vector space isomorphic to self/other.
   *
   * @param other - A subspace of self
   * @returns The quotient space
   *
   * @see Reference: sage/modules/free_module.py:FreeModule_generic_field.quotient
   */
  quotient(other: FreeModuleField): FreeModuleField {
    // The quotient V/W has dimension dim(V) - dim(W)
    if (!other.isSubspace(this)) {
      throw new ArithmeticError('other must be a subspace of self');
    }

    const quotientDim = this.dimension() - other.dimension();

    if (quotientDim === 0) {
      return this.subspace([]);
    }

    // Find a complement of other in self
    // This is a set of vectors that together with other's basis spans self
    const selfBasis = this.basis();
    const otherBasis = other.basis();

    // Find vectors in self's basis not in other's span
    const complementVectors: FreeModuleElement[] = [];

    for (const v of selfBasis) {
      // Check if v is in other's span
      try {
        other.coordinates(v, true);
        // v is in other, skip it
      } catch {
        // v is not in other, add to complement
        complementVectors.push(v);
        if (complementVectors.length >= quotientDim) {
          break;
        }
      }
    }

    // Return the span of complement vectors
    return this.subspace(complementVectors);
  }

  /**
   * Return an iterator over subspaces of given dimension.
   *
   * Note: This is only implemented for finite fields where we can
   * enumerate all subspaces.
   *
   * @param dim - The dimension
   * @yields Subspaces of the given dimension
   *
   * @see Reference: sage/modules/free_module.py:FreeModule_generic_field.subspaces
   */
  *subspaces(dim: number): IterableIterator<FreeModuleField> {
    if (dim < 0 || dim > this.dimension()) {
      return;
    }

    if (dim === 0) {
      yield this.subspace([]);
      return;
    }

    if (dim === this.dimension()) {
      yield this;
      return;
    }

    // For general fields, enumerating all subspaces is complex
    // This requires enumerating Grassmannian, which needs more infrastructure
    throw new NotImplementedError(
      'subspaces iteration is only implemented for finite fields with enumerable elements'
    );
  }

  /**
   * Return a complement of this subspace.
   *
   * The complement is a subspace W such that self + W = ambient and
   * self intersection W = {0}.
   *
   * @returns A complement subspace
   *
   * @see Reference: sage/modules/free_module.py:FreeModule_generic_field.complement
   */
  complement(): FreeModuleField {
    // Simple cases
    if (this.dimension() === 0) {
      return this.ambientVectorSpace();
    }

    if (this.dimension() === this.ambientVectorSpace().dimension()) {
      return this.subspace([]);
    }

    // The orthogonal complement with respect to the standard inner product
    // is a complement (though not unique)
    return this.orthogonalComplement();
  }
}

// ============================================================================
// Ambient modules
// ============================================================================

/**
 * Ambient free module over a ring.
 * @see Reference: sage/modules/free_module.py:FreeModule_ambient
 */
export class FreeModuleAmbient extends FreeModuleGeneric {
  constructor(
    baseRing: RingLike,
    rank: number,
    sparse: boolean = false,
    innerProductMatrix?: unknown
  ) {
    super(baseRing, rank, rank, sparse, undefined, innerProductMatrix);
  }

  /**
   * Return True since this is an ambient module.
   */
  override isAmbient(): boolean {
    return true;
  }
}

/**
 * Ambient free module over a PID.
 *
 * `FreeModule_ambient_pid` derives from `FreeModule_generic_pid` in SageMath,
 * so the full PID interface (span_of_basis, saturation, index_in, ...) is
 * available on ZZ^n.
 *
 * @see Reference: sage/modules/free_module.py:FreeModule_ambient_pid
 */
export class FreeModuleAmbientPID extends FreeModulePID {
  constructor(
    baseRing: RingLike,
    rank: number,
    sparse: boolean = false,
    innerProductMatrix?: unknown
  ) {
    super(baseRing, rank, rank, sparse, undefined, innerProductMatrix);
  }

  /**
   * Return True since this is an ambient module.
   */
  override isAmbient(): boolean {
    return true;
  }
}

/**
 * Ambient vector space over a field.
 * @see Reference: sage/modules/free_module.py:FreeModule_ambient_field
 */
export class FreeModuleAmbientField extends FreeModuleField {
  constructor(
    baseRing: RingLike,
    dimension: number,
    sparse: boolean = false,
    innerProductMatrix?: unknown
  ) {
    super(baseRing, dimension, dimension, sparse, innerProductMatrix);
  }

  /**
   * Return True since this is an ambient module.
   */
  override isAmbient(): boolean {
    return true;
  }
}

// ============================================================================
// Submodules
// ============================================================================

/**
 * Create a vector with the given entries, without coercing them.
 *
 * The entries produced by the exact linear algebra above already lie in the
 * coordinate ring of the module.
 */
function makeVector(parent: FreeModuleGeneric, entries: unknown[]): FreeModuleElement {
  const v = parent.isSparse()
    ? new FreeModuleElementSparse(parent, entries)
    : new FreeModuleElementDense(parent, entries);
  v.setImmutable();
  return v;
}

/**
 * Compute the user basis and the echelonized basis of a submodule.
 *
 * @see Reference: sage/modules/free_module.py:FreeModule_submodule_with_basis_pid.__init__
 */
function submoduleBases(
  ambient: FreeModuleGeneric,
  basis: FreeModuleElement[],
  options?: { check?: boolean; echelonize?: boolean; alreadyEchelonized?: boolean }
): { user: unknown[][]; echelonized: unknown[][] | null } {
  const ar = arithmeticFor(ambient.baseRing());
  const rows = basis.map((b) => b.list());

  if (options?.alreadyEchelonized) {
    return { user: rows, echelonized: rows };
  }
  if (options?.echelonize) {
    const E = echelonRows(rows, ar);
    return { user: E, echelonized: E };
  }
  if ((options?.check ?? true) && ar.exact && rows.length > 0) {
    if (rankOfRows(rows, ar) !== rows.length) {
      throw new ValueError('the given basis vectors must be linearly independent');
    }
  }
  return { user: rows, echelonized: null };
}

/**
 * A submodule of a free module over a general ring.
 *
 * Over a ring that is not a PID no echelon form is available, so the
 * generators are stored verbatim, exactly as in SageMath's
 * `Submodule_free_ambient`.
 *
 * @see Reference: sage/modules/free_module.py:Submodule_free_ambient
 */
export class FreeModuleSubmodule extends FreeModuleGeneric {
  protected _ambient: FreeModuleGeneric;
  protected _userBasis: FreeModuleElement[];

  constructor(
    ambient: FreeModuleGeneric,
    gens: FreeModuleElement[],
    _options?: { check?: boolean; echelonize?: boolean; alreadyEchelonized?: boolean }
  ) {
    super(
      ambient.baseRing(),
      gens.length,
      ambient.degree(),
      ambient.isSparse(),
      undefined,
      ambient.innerProductMatrix()
    );

    this._ambient = ambient;
    this._userBasis = gens.map((g) => makeVector(this, g.list()));
    this._basis = this._userBasis;
  }

  override ambientModule(): FreeModuleGeneric {
    return this._ambient;
  }

  override hasUserBasis(): boolean {
    return true;
  }
}

/**
 * Submodule of a free module over a PID with a user-specified basis.
 * @see Reference: sage/modules/free_module.py:FreeModule_submodule_with_basis_pid
 */
export class FreeModuleWithBasis extends FreeModulePID {
  protected _ambient: FreeModuleGeneric;
  protected _userBasis: FreeModuleElement[];
  protected _echelonizedBasisMatrix: unknown[][] | null = null;

  constructor(
    ambient: FreeModuleGeneric,
    basis: FreeModuleElement[],
    options?: { check?: boolean; echelonize?: boolean; alreadyEchelonized?: boolean }
  ) {
    const { user, echelonized } = submoduleBases(ambient, basis, options);
    super(
      ambient.baseRing(),
      user.length,
      ambient.degree(),
      ambient.isSparse(),
      undefined,
      ambient.innerProductMatrix()
    );

    this._ambient = ambient;
    this._userBasis = user.map((row) => makeVector(this, row));
    this._basis = this._userBasis;
    this._echelonizedBasisMatrix = echelonized;
  }

  /**
   * Return True since this module has a user-specified basis.
   */
  override hasUserBasis(): boolean {
    return true;
  }

  override ambientModule(): FreeModuleGeneric {
    return this._ambient;
  }

  /**
   * Return the basis matrix for self in row echelon form.
   *
   * @see Reference: sage/modules/free_module.py:FreeModule_submodule_with_basis_pid.echelonized_basis_matrix
   */
  override echelonizedBasisMatrix(): unknown[][] {
    if (this._echelonizedBasisMatrix === null) {
      this._echelonizedBasisMatrix = echelonRows(
        this.basisMatrix() as unknown[][],
        arithmeticFor(this._baseRing)
      );
    }
    return this._echelonizedBasisMatrix;
  }
}

/**
 * An R-submodule of K^n where K is the fraction field of the PID R, given by
 * generators.  Its basis is the echelon form of the generating matrix.
 *
 * @see Reference: sage/modules/free_module.py:FreeModule_submodule_pid
 */
export class FreeModuleSubmodulePID extends FreeModuleWithBasis {
  constructor(
    ambient: FreeModuleGeneric,
    gens: FreeModuleElement[],
    options?: { check?: boolean; alreadyEchelonized?: boolean }
  ) {
    super(ambient, gens, {
      check: options?.check ?? true,
      echelonize: !(options?.alreadyEchelonized ?? false),
      alreadyEchelonized: options?.alreadyEchelonized ?? false,
    });
  }

  /**
   * Return False: the basis is the echelon form, not a user basis.
   */
  override hasUserBasis(): boolean {
    return false;
  }
}

/**
 * Subspace of a vector space with a user-specified basis.
 * @see Reference: sage/modules/free_module.py:FreeModule_submodule_with_basis_field
 */
export class FreeModuleSubspaceWithBasis extends FreeModuleField {
  protected _ambient: FreeModuleField;
  protected _userBasis: FreeModuleElement[];
  protected _echelonizedBasisMatrix: unknown[][] | null = null;

  constructor(
    ambient: FreeModuleField,
    basis: FreeModuleElement[],
    options?: { check?: boolean; echelonize?: boolean; alreadyEchelonized?: boolean }
  ) {
    const { user, echelonized } = submoduleBases(ambient, basis, options);
    super(
      ambient.baseRing(),
      user.length,
      ambient.degree(),
      ambient.isSparse(),
      ambient.innerProductMatrix()
    );

    this._ambient = ambient;
    this._userBasis = user.map((row) => makeVector(this, row));
    this._basis = this._userBasis;
    this._echelonizedBasisMatrix = echelonized;
  }

  override hasUserBasis(): boolean {
    return true;
  }

  override ambientModule(): FreeModuleGeneric {
    return this._ambient;
  }

  override echelonizedBasisMatrix(): unknown[][] {
    if (this._echelonizedBasisMatrix === null) {
      this._echelonizedBasisMatrix = echelonRows(
        this.basisMatrix() as unknown[][],
        arithmeticFor(this._baseRing)
      );
    }
    return this._echelonizedBasisMatrix;
  }
}

/**
 * A subspace of a vector space, given by generators; its basis is the reduced
 * row echelon form of the generating matrix.
 *
 * @see Reference: sage/modules/free_module.py:FreeModule_submodule_field
 */
export class FreeModuleSubspace extends FreeModuleSubspaceWithBasis {
  constructor(
    ambient: FreeModuleField,
    gens: FreeModuleElement[],
    options?: { check?: boolean; alreadyEchelonized?: boolean }
  ) {
    super(ambient, gens, {
      check: options?.check ?? true,
      echelonize: !(options?.alreadyEchelonized ?? false),
      alreadyEchelonized: options?.alreadyEchelonized ?? false,
    });
  }

  /**
   * Return False: the basis is the echelon form, not a user basis.
   */
  override hasUserBasis(): boolean {
    return false;
  }
}

// ============================================================================
// Quotient Modules
// ============================================================================

/**
 * Quotient of a free module by a submodule.
 *
 * The quotient M/N is represented by cosets v + N for v in M.
 *
 * @see Reference: sage/modules/free_module.py:FreeModule_quotient
 */
export class FreeModuleQuotient extends FreeModuleGeneric {
  protected _cover: FreeModuleGeneric;
  protected _submodule: FreeModuleGeneric;

  constructor(cover: FreeModuleGeneric, submodule: FreeModuleGeneric) {
    // The rank of the quotient is rank(cover) - rank(submodule)
    const quotientRank = cover.rank() - submodule.rank();

    super(cover.baseRing(), quotientRank >= 0 ? quotientRank : 0, cover.degree(), cover.isSparse());

    this._cover = cover;
    this._submodule = submodule;
  }

  /**
   * Return the covering module M.
   */
  coveringModule(): FreeModuleGeneric {
    return this._cover;
  }

  /**
   * Return the submodule N that we are quotienting by.
   */
  relations(): FreeModuleGeneric {
    return this._submodule;
  }

  /**
   * Return the lift of an element from the quotient to the covering module.
   * @param v - An element of the quotient
   * @returns A representative in the covering module
   */
  lift(v: FreeModuleElement): FreeModuleElement {
    // For now, just return the element as-is since we're working with
    // representatives
    return v;
  }

  /**
   * Return the projection of an element to the quotient.
   * @param v - An element of the covering module
   * @returns The coset representative in the quotient
   */
  project(v: FreeModuleElement): FreeModuleElement {
    // Return the element as the coset representative
    // A full implementation would reduce v modulo the submodule
    return v;
  }

  override toString(): string {
    const ringName = this._baseRing.toString?.() ?? 'Ring';
    return `Quotient of ${this._cover.toString()} by ${this._submodule.toString()}`;
  }
}
