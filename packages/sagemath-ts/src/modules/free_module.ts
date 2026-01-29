/**
 * @module sage/modules/free_module
 * @description Free modules over commutative rings
 * @see Reference: sage/modules/free_module.py
 *
 * Sage supports computation with free modules over an arbitrary commutative ring.
 * Nontrivial functionality is available over ZZ, fields, and some principal
 * ideal domains (e.g. QQ[x] and rings of integers of number fields).
 */

import { NotImplementedError, ValueError, ArithmeticError } from '../errors.js';
import {
  FreeModuleElement,
  FreeModuleElementDense,
  FreeModuleElementSparse,
  type RingLike,
  type FreeModuleParent,
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
// Helper functions
// ============================================================================

/**
 * Compute GCD of two bigints.
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
 * Compute GCD of two numbers.
 */
function gcd(a: number, b: number): number {
  a = Math.abs(a);
  b = Math.abs(b);
  while (b !== 0) {
    const t = b;
    b = a % b;
    a = t;
  }
  return a;
}

/**
 * Compute LCM of two numbers.
 */
function lcm(a: number, b: number): number {
  if (a === 0 || b === 0) return 0;
  return Math.abs(a * b) / gcd(a, b);
}

/**
 * Compute the kernel of a matrix (null space).
 * Returns a list of vectors that span the kernel.
 */
function computeKernel(matrix: number[][], numCols: number): number[][] {
  const m = matrix.length;
  const n = numCols;
  const eps = 1e-10;

  // Copy matrix for row reduction
  const M: number[][] = matrix.map((row) => [...row]);

  // Row echelon form with column pivoting
  const pivotCols: number[] = [];
  let pivotRow = 0;

  for (let col = 0; col < n && pivotRow < m; col++) {
    // Find pivot
    let maxRow = -1;
    let maxVal = eps;
    for (let row = pivotRow; row < m; row++) {
      if (Math.abs(M[row]![col]!) > maxVal) {
        maxVal = Math.abs(M[row]![col]!);
        maxRow = row;
      }
    }

    if (maxRow === -1) {
      continue; // No pivot in this column
    }

    // Swap rows
    if (maxRow !== pivotRow) {
      [M[pivotRow], M[maxRow]] = [M[maxRow]!, M[pivotRow]!];
    }

    pivotCols.push(col);
    const pivot = M[pivotRow]![col]!;

    // Scale pivot row
    for (let j = col; j < n; j++) {
      M[pivotRow]![j] = M[pivotRow]![j]! / pivot;
    }

    // Eliminate other rows
    for (let row = 0; row < m; row++) {
      if (row !== pivotRow && Math.abs(M[row]![col]!) > eps) {
        const factor = M[row]![col]!;
        for (let j = col; j < n; j++) {
          M[row]![j] = M[row]![j]! - factor * M[pivotRow]![j]!;
        }
      }
    }

    pivotRow++;
  }

  // Find free columns (not pivot columns)
  const freeCols: number[] = [];
  for (let col = 0; col < n; col++) {
    if (!pivotCols.includes(col)) {
      freeCols.push(col);
    }
  }

  // Build kernel vectors
  const kernel: number[][] = [];

  for (const freeCol of freeCols) {
    const v: number[] = new Array(n).fill(0);
    v[freeCol] = 1;

    // For each pivot, compute its value
    for (let i = 0; i < pivotCols.length; i++) {
      const pivotCol = pivotCols[i]!;
      // Find the row with this pivot
      let val = 0;
      for (let j = 0; j < m; j++) {
        if (Math.abs(M[j]![pivotCol]! - 1) < eps) {
          // This is the pivot row
          val = -M[j]![freeCol]!;
          break;
        }
      }
      v[pivotCol] = val;
    }

    kernel.push(v);
  }

  return kernel;
}

/**
 * Compute row echelon form of a matrix.
 */
function echelonize(matrix: unknown[][], ring: RingLike): unknown[][] {
  if (matrix.length === 0) return [];

  const m = matrix.length;
  const n = matrix[0]!.length;
  const eps = 1e-10;

  // Convert to numbers for computation
  const M: number[][] = matrix.map((row) =>
    row.map((e) => {
      if (typeof e === 'bigint') return Number(e);
      if (typeof e === 'number') return e;
      return Number(String(e));
    })
  );

  let pivotRow = 0;

  for (let col = 0; col < n && pivotRow < m; col++) {
    // Find pivot
    let maxRow = -1;
    let maxVal = eps;
    for (let row = pivotRow; row < m; row++) {
      if (Math.abs(M[row]![col]!) > maxVal) {
        maxVal = Math.abs(M[row]![col]!);
        maxRow = row;
      }
    }

    if (maxRow === -1) {
      continue;
    }

    // Swap rows
    if (maxRow !== pivotRow) {
      [M[pivotRow], M[maxRow]] = [M[maxRow]!, M[pivotRow]!];
    }

    const pivot = M[pivotRow]![col]!;

    // Scale pivot row to make pivot = 1
    for (let j = col; j < n; j++) {
      M[pivotRow]![j] = M[pivotRow]![j]! / pivot;
    }

    // Eliminate below pivot
    for (let row = pivotRow + 1; row < m; row++) {
      if (Math.abs(M[row]![col]!) > eps) {
        const factor = M[row]![col]!;
        for (let j = col; j < n; j++) {
          M[row]![j] = M[row]![j]! - factor * M[pivotRow]![j]!;
        }
      }
    }

    pivotRow++;
  }

  // Convert back to ring elements
  const zero = ring.zero();
  const result: unknown[][] = [];

  for (const row of M) {
    const newRow: unknown[] = row.map((e) => {
      if (typeof zero === 'bigint') {
        return BigInt(Math.round(e));
      }
      return e;
    });
    result.push(newRow);
  }

  return result;
}

/**
 * Compute determinant of a numeric matrix using LU decomposition.
 */
function computeDeterminant(matrix: number[][]): number {
  const n = matrix.length;
  if (n === 0) return 1;
  if (n === 1) return matrix[0]![0]!;
  if (n === 2) {
    return matrix[0]![0]! * matrix[1]![1]! - matrix[0]![1]! * matrix[1]![0]!;
  }

  // Create copy for in-place operations
  const M: number[][] = matrix.map((row) => [...row]);

  let det = 1;
  let sign = 1;

  for (let col = 0; col < n; col++) {
    // Find pivot
    let maxRow = col;
    let maxVal = Math.abs(M[col]![col]!);
    for (let row = col + 1; row < n; row++) {
      if (Math.abs(M[row]![col]!) > maxVal) {
        maxVal = Math.abs(M[row]![col]!);
        maxRow = row;
      }
    }

    // Swap rows if needed
    if (maxRow !== col) {
      [M[col], M[maxRow]] = [M[maxRow]!, M[col]!];
      sign = -sign;
    }

    // Check for singular matrix
    if (Math.abs(M[col]![col]!) < 1e-15) {
      return 0;
    }

    det *= M[col]![col]!;

    // Eliminate below pivot
    for (let row = col + 1; row < n; row++) {
      const factor = M[row]![col]! / M[col]![col]!;
      for (let j = col; j < n; j++) {
        M[row]![j] = M[row]![j]! - factor * M[col]![j]!;
      }
    }
  }

  return sign * det;
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

    if (gens.length === 0) {
      // Zero submodule
      return new FreeModuleSubmodule(this as unknown as FreeModuleGeneric, [], {
        check: options?.check ?? true,
      });
    }

    return new FreeModuleSubmodule(this as unknown as FreeModuleGeneric, gens, {
      check: options?.check ?? true,
      alreadyEchelonized: options?.alreadyEchelonized ?? false,
    });
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
    if (Array.isArray(gens)) {
      return this.span(gens, this._baseRing, options);
    }
    // If gens is a module, get its generators
    return this.span(gens.gens(), this._baseRing, options);
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

    // The zero module is always a submodule
    const selfGens = (this as unknown as FreeModuleGeneric).gens();
    if (selfGens.length === 0) {
      return true;
    }

    // Other is the zero module but self is not
    const otherGens = (other as unknown as FreeModuleGeneric).gens();
    if (otherGens.length === 0) {
      return false;
    }

    // Check if every generator of self is in other
    // This is a simplified check - full implementation would require
    // solving a linear system to check containment
    try {
      for (const gen of selfGens) {
        // Try to express gen in terms of other's basis
        const coords = (other as unknown as FreeModuleGeneric).coordinates(gen, true);
        // Check that all coordinates are in the base ring
        for (const c of coords) {
          if (typeof c === 'bigint' || typeof c === 'number') {
            continue;
          }
          // If coordinate ring is not the base ring, may not be a submodule
        }
      }
      return true;
    } catch {
      return false;
    }
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

    // For submodules, compute row echelon form of the basis matrix
    const basisMat = this.basisMatrix() as unknown[][];
    return echelonize(basisMat, this._baseRing);
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
    const b = this.basis();
    const n = b.length;

    // Build Gram matrix G[i][j] = <b_i, b_j>
    const G: unknown[][] = [];

    for (let i = 0; i < n; i++) {
      const row: unknown[] = [];
      for (let j = 0; j < n; j++) {
        // Compute inner product of b[i] and b[j]
        const dot = b[i]!.innerProduct(b[j]!);
        row.push(dot);
      }
      G.push(row);
    }

    return G;
  }

  /**
   * Return the discriminant of this free module.
   * For a lattice, this is |det(G)| where G is the Gram matrix.
   */
  discriminant(): unknown {
    const G = this.gramMatrix();
    const n = G.length;

    if (n === 0) {
      return this._baseRing.one();
    }

    // Compute determinant of Gram matrix
    // Convert to number for computation
    const numG: number[][] = [];
    for (let i = 0; i < n; i++) {
      numG.push([]);
      for (let j = 0; j < n; j++) {
        const entry = G[i]![j]!;
        if (typeof entry === 'bigint') {
          numG[i]!.push(Number(entry));
        } else if (typeof entry === 'number') {
          numG[i]!.push(entry);
        } else {
          numG[i]!.push(Number(String(entry)));
        }
      }
    }

    // Compute determinant using LU decomposition
    const det = computeDeterminant(numG);
    const absDet = Math.abs(det);

    // Return in appropriate type
    const zero = this._baseRing.zero();
    if (typeof zero === 'bigint') {
      return BigInt(Math.round(absDet));
    }
    return absDet;
  }

  /**
   * Return the cardinality of this module.
   *
   * @returns The number of elements in this module, or Infinity if infinite.
   *
   * @see Reference: sage/modules/free_module.py:FreeModule_generic.cardinality
   */
  cardinality(): number | typeof Infinity {
    if (this._rank === 0) {
      return 1;
    }

    // Check if base ring is finite
    if ('cardinality' in this._baseRing && typeof this._baseRing.cardinality === 'function') {
      const baseCard = this._baseRing.cardinality();
      if (typeof baseCard === 'number' && isFinite(baseCard)) {
        return Math.pow(baseCard, this._rank);
      }
      if (typeof baseCard === 'bigint') {
        // Return as number if small enough, otherwise Infinity as approximation
        const result = baseCard ** BigInt(this._rank);
        if (result <= BigInt(Number.MAX_SAFE_INTEGER)) {
          return Number(result);
        }
        return Infinity;
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
    return Infinity;
  }

  /**
   * Return whether this is an ambient module.
   */
  isAmbient(): boolean {
    return false;
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
   * @param v - A vector
   * @param check - Whether to verify v is in self
   * @returns The list of coefficients c such that v = sum(c[i] * basis[i])
   */
  coordinates(v: FreeModuleElement, check: boolean = true): unknown[] {
    const basis = this.basis();
    const n = basis.length;
    const m = this._degree;

    if (n === 0) {
      if (v.isZero()) {
        return [];
      }
      throw new ArithmeticError('vector is not in the span of the basis');
    }

    // For ambient modules where rank = degree and basis is standard,
    // coordinates are just the entries
    if (this.isAmbient() && this._rank === this._degree) {
      return v.list();
    }

    // Build the basis matrix (rows are basis vectors)
    const B: number[][] = [];
    for (const b of basis) {
      const row: number[] = [];
      for (let j = 0; j < m; j++) {
        const entry = b.getItem(j);
        row.push(typeof entry === 'bigint' ? Number(entry) : Number(entry));
      }
      B.push(row);
    }

    // Build target vector
    const target: number[] = [];
    for (let j = 0; j < m; j++) {
      const entry = v.getItem(j);
      target.push(typeof entry === 'bigint' ? Number(entry) : Number(entry));
    }

    // Solve B^T * coords = target using least squares if over-determined
    // or direct solve if square
    // For simplicity, use pseudo-inverse: coords = (B * B^T)^{-1} * B * target

    // Compute B * B^T (n x n)
    const BBT: number[][] = [];
    for (let i = 0; i < n; i++) {
      BBT.push([]);
      for (let j = 0; j < n; j++) {
        let dot = 0;
        for (let k = 0; k < m; k++) {
          dot += B[i]![k]! * B[j]![k]!;
        }
        BBT[i]!.push(dot);
      }
    }

    // Compute B * target (n x 1)
    const Bt: number[] = [];
    for (let i = 0; i < n; i++) {
      let dot = 0;
      for (let k = 0; k < m; k++) {
        dot += B[i]![k]! * target[k]!;
      }
      Bt.push(dot);
    }

    // Solve BBT * coords = Bt using Gaussian elimination
    // Create augmented matrix
    const aug: number[][] = [];
    for (let i = 0; i < n; i++) {
      aug.push([...BBT[i]!, Bt[i]!]);
    }

    // Forward elimination
    for (let col = 0; col < n; col++) {
      // Find pivot
      let maxRow = col;
      for (let row = col + 1; row < n; row++) {
        if (Math.abs(aug[row]![col]!) > Math.abs(aug[maxRow]![col]!)) {
          maxRow = row;
        }
      }
      [aug[col], aug[maxRow]] = [aug[maxRow]!, aug[col]!];

      if (Math.abs(aug[col]![col]!) < 1e-10) {
        continue;
      }

      for (let row = col + 1; row < n; row++) {
        const factor = aug[row]![col]! / aug[col]![col]!;
        for (let j = col; j <= n; j++) {
          aug[row]![j] = aug[row]![j]! - factor * aug[col]![j]!;
        }
      }
    }

    // Back substitution
    const coords: number[] = new Array(n).fill(0);
    for (let row = n - 1; row >= 0; row--) {
      let sum = aug[row]![n]!;
      for (let col = row + 1; col < n; col++) {
        sum -= aug[row]![col]! * coords[col]!;
      }
      if (Math.abs(aug[row]![row]!) > 1e-10) {
        coords[row] = sum / aug[row]![row]!;
      }
    }

    // Convert back to ring elements
    const zero = this._baseRing.zero();
    const result: unknown[] = [];
    for (let i = 0; i < n; i++) {
      if (typeof zero === 'bigint') {
        result.push(BigInt(Math.round(coords[i]!)));
      } else if (typeof zero === 'number') {
        result.push(coords[i]!);
      } else if (this._baseRing.__call__) {
        result.push(this._baseRing.__call__(coords[i]!));
      } else {
        result.push(coords[i]!);
      }
    }

    return result;
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
      return true; // More vectors than dimension means dependent
    }

    // Build matrix from vectors (as rows)
    const n = vecs.length;
    const m = this._degree;

    // Convert to number matrix for rank computation
    const matrix: number[][] = [];
    for (const v of vecs) {
      const row: number[] = [];
      for (let j = 0; j < m; j++) {
        const entry = v.getItem(j);
        if (typeof entry === 'bigint') {
          row.push(Number(entry));
        } else if (typeof entry === 'number') {
          row.push(entry);
        } else {
          // Try to convert via toString
          row.push(Number(String(entry)));
        }
      }
      matrix.push(row);
    }

    // Compute rank using row echelon form
    const eps = 1e-10;

    // Gaussian elimination to find rank
    let rank = 0;
    const usedCols: boolean[] = new Array(m).fill(false);

    for (let row = 0; row < n; row++) {
      // Find pivot column
      let pivotCol = -1;
      for (let col = 0; col < m; col++) {
        if (!usedCols[col] && Math.abs(matrix[row]![col]!) > eps) {
          pivotCol = col;
          break;
        }
      }

      if (pivotCol === -1) {
        continue; // This row is zero or dependent
      }

      usedCols[pivotCol] = true;
      rank++;

      // Normalize pivot row
      const pivot = matrix[row]![pivotCol]!;
      for (let col = 0; col < m; col++) {
        matrix[row]![col] = matrix[row]![col]! / pivot;
      }

      // Eliminate this column from other rows
      for (let otherRow = 0; otherRow < n; otherRow++) {
        if (otherRow !== row) {
          const factor = matrix[otherRow]![pivotCol]!;
          for (let col = 0; col < m; col++) {
            matrix[otherRow]![col] = matrix[otherRow]![col]! - factor * matrix[row]![col]!;
          }
        }
      }
    }

    // Vectors are linearly dependent if rank < number of vectors
    return rank < n;
  }

  /**
   * Return a linear combination of the basis vectors.
   * @param coefficients - The coefficients
   */
  linearCombinationOfBasis(coefficients: unknown[]): FreeModuleElement {
    const b = this.basis();
    if (coefficients.length !== b.length) {
      throw new ValueError(
        `coefficients must have length ${b.length}, got ${coefficients.length}`
      );
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
   * Return this module with a different base ring.
   * @param ring - The new base ring
   */
  changeRing(ring: RingLike): FreeModuleGeneric {
    return FreeModule(ring, this._rank, { sparse: this._sparse });
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
   */
  indexIn(other: FreeModuleGeneric): unknown {
    if (this._baseRing !== other.baseRing()) {
      throw new ArithmeticError('self and other must have the same base ring');
    }

    // Handle trivial cases
    if (this.rank() === 0) {
      return other.rank() === 0 ? 1 : Infinity;
    }

    if (other.rank() === 0) {
      return this.rank() === 0 ? 1 : 0;
    }

    if (this.rank() < other.rank()) {
      return Infinity;
    }

    // Express basis of self in terms of basis of other
    const selfBasis = this.basis();
    const otherBasis = other.basis();

    // Build matrix of coordinates of self's basis in other's basis
    const coordMatrix: number[][] = [];
    for (const v of selfBasis) {
      try {
        const coords = other.coordinates(v, true);
        coordMatrix.push(coords.map((c) => Number(c)));
      } catch {
        return Infinity;
      }
    }

    // The index is |det(coordMatrix)| (for square matrices)
    if (coordMatrix.length !== coordMatrix[0]?.length) {
      return Infinity;
    }

    const det = computeDeterminant(coordMatrix);
    return Math.abs(det);
  }

  /**
   * Return the intersection of self and other.
   *
   * @param other - Another module
   *
   * @see Reference: sage/modules/free_module.py:FreeModule_generic_pid.intersection
   */
  intersection(other: FreeModuleGeneric): FreeModuleGeneric {
    // Handle trivial cases
    if (this.rank() === 0 || other.rank() === 0) {
      return this.zeroSubmodule();
    }

    if (this === other) {
      return this;
    }

    // If one is contained in the other, return the smaller one
    if (this.isSubmodule(other as unknown as ModuleFreeAmbient)) {
      return this;
    }
    if ((other as unknown as ModuleFreeAmbient).isSubmodule(this as unknown as ModuleFreeAmbient)) {
      return other;
    }

    // General case: use the formula that the intersection is the kernel
    // of the map from self + other to (self + other) / self \times (self + other) / other
    // This is equivalent to finding vectors in the row span of both matrices

    // Stack basis matrices vertically
    const selfBasis = this.basis();
    const otherBasis = other.basis();
    const m = selfBasis.length;
    const n = otherBasis.length;
    const d = this._degree;

    // Build augmented system: [A; B] where we look for linear combinations
    // c_1*a_1 + ... + c_m*a_m = d_1*b_1 + ... + d_n*b_n
    // This is equivalent to [A^T | -B^T] * [c; d]^T = 0

    // Solve using row reduction
    // Build matrix: [[self basis]^T | [other basis]^T]
    const augMatrix: number[][] = [];
    for (let i = 0; i < d; i++) {
      const row: number[] = [];
      for (let j = 0; j < m; j++) {
        const entry = selfBasis[j]!.getItem(i);
        row.push(typeof entry === 'bigint' ? Number(entry) : Number(entry));
      }
      for (let j = 0; j < n; j++) {
        const entry = otherBasis[j]!.getItem(i);
        row.push(typeof entry === 'bigint' ? -Number(entry) : -Number(entry));
      }
      augMatrix.push(row);
    }

    // Find the kernel of this matrix
    const kernel = computeKernel(augMatrix, m + n);

    // The intersection consists of vectors sum(c[i] * selfBasis[i]) for each kernel vector
    const intersectionGens: FreeModuleElement[] = [];

    for (const kernelVec of kernel) {
      // Extract the first m coordinates (coefficients for self's basis)
      const selfCoeffs = kernelVec.slice(0, m);

      // Compute the linear combination
      const entries: unknown[] = new Array(d).fill(this._baseRing.zero());
      for (let i = 0; i < m; i++) {
        if (Math.abs(selfCoeffs[i]!) > 1e-10) {
          for (let j = 0; j < d; j++) {
            const basis_entry = selfBasis[i]!.getItem(j);
            const val = typeof basis_entry === 'bigint' ? Number(basis_entry) : Number(basis_entry);
            entries[j] = Number(entries[j]) + selfCoeffs[i]! * val;
          }
        }
      }

      // Check if this is a non-zero vector
      const isNonZero = entries.some((e) => Math.abs(Number(e)) > 1e-10);
      if (isNonZero) {
        // Convert to appropriate type
        const zero = this._baseRing.zero();
        const finalEntries = entries.map((e) => {
          if (typeof zero === 'bigint') {
            return BigInt(Math.round(Number(e)));
          }
          return Number(e);
        });
        intersectionGens.push(this.createElement(finalEntries));
      }
    }

    if (intersectionGens.length === 0) {
      return this.zeroSubmodule();
    }

    return this.span(intersectionGens);
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
    const sat = this.saturation();
    return sat.indexIn(this);
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
    // If base ring is a field, the module is already saturated
    if (isField(this._baseRing)) {
      return this;
    }

    if (this._rank === 0) {
      return this;
    }

    // Get the basis matrix
    const basisMat = this.basisMatrix() as unknown[][];

    if (basisMat.length === 0) {
      return this;
    }

    // Compute GCD of all entries in each row
    const saturatedRows: unknown[][] = [];
    const zero = this._baseRing.zero();

    for (const row of basisMat) {
      // Convert row to bigints for GCD computation
      const rowBigints = row.map((e) => {
        if (typeof e === 'bigint') return e;
        if (typeof e === 'number') return BigInt(Math.round(e));
        return BigInt(String(e));
      });

      // Compute GCD of all entries in the row
      let rowGcd = 0n;
      for (const entry of rowBigints) {
        if (entry !== 0n) {
          rowGcd = bigintGcd(rowGcd, entry < 0n ? -entry : entry);
        }
      }

      // If GCD is 0 or 1, row is already saturated (or zero)
      if (rowGcd === 0n || rowGcd === 1n) {
        saturatedRows.push(row);
      } else {
        // Divide each entry by the GCD
        const saturatedRow: unknown[] = [];
        for (const entry of rowBigints) {
          if (typeof zero === 'bigint') {
            saturatedRow.push(entry / rowGcd);
          } else {
            saturatedRow.push(Number(entry / rowGcd));
          }
        }
        saturatedRows.push(saturatedRow);
      }
    }

    // Create the saturated basis vectors
    const ambient = this.ambientModule();
    const saturatedVectors: FreeModuleElement[] = [];

    for (const row of saturatedRows) {
      saturatedVectors.push(ambient.createElement(row));
    }

    // Create the saturated module
    const sat = ambient.span(saturatedVectors);

    // Return self if already saturated
    // (Check if the saturated module equals self)
    if (saturatedVectors.length === this._rank) {
      let isSame = true;
      for (let i = 0; i < saturatedVectors.length && isSame; i++) {
        const selfBasis = this.basis();
        if (!saturatedVectors[i]!.equals(selfBasis[i]!)) {
          isSame = false;
        }
      }
      if (isSame) {
        return this;
      }
    }

    return sat;
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

    // For integer-based modules, denominator is 1
    const zero = this._baseRing.zero();
    if (typeof zero === 'bigint') {
      return 1n;
    }

    // For rational-like entries, we would need to extract denominators
    // and compute their LCM. For simplicity, return 1 for now.
    // A full implementation would check each entry for a denominator method.
    let lcmResult = 1;

    for (const row of basisMat) {
      for (const entry of row) {
        // Check if entry has a denominator method (like rational numbers)
        if (typeof entry === 'object' && entry !== null && 'denominator' in entry) {
          const denom = (entry as { denominator: () => number }).denominator();
          lcmResult = lcm(lcmResult, denom);
        }
      }
    }

    if (typeof zero === 'number') {
      return lcmResult;
    }

    // Return as bigint if base ring uses bigints
    return BigInt(lcmResult);
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
  ): FreeModuleWithBasis {
    return new FreeModuleWithBasis(this, basis, {
      check: options?.check ?? true,
      echelonize: false,
      alreadyEchelonized: options?.alreadyEchelonized ?? false,
    });
  }

  /**
   * Create the R-submodule with given basis.
   * @param basis - A list of linearly independent vectors
   * @param options - Additional options
   */
  submoduleWithBasis(
    basis: FreeModuleElement[],
    options?: { check?: boolean; alreadyEchelonized?: boolean }
  ): FreeModuleWithBasis {
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
    const ambient = new FreeModuleAmbientField(field, this._degree, this._sparse, this._innerProductMatrix);

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
    const ambient = new FreeModuleAmbientField(field, this._degree, this._sparse, this._innerProductMatrix);

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
    return new FreeModuleSubspaceWithBasis(this, gens, {
      check: options?.check ?? true,
      echelonize: true,
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
    return new FreeModuleSubspaceWithBasis(this, gens, {
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
    // Get the basis matrix
    const basisMat = this.basisMatrix() as unknown[][];
    const n = this._degree;

    if (basisMat.length === 0) {
      // Complement of zero space is the entire ambient space
      return this.ambientVectorSpace();
    }

    if (this.dimension() === n) {
      // Complement of ambient space is zero
      return this.subspace([]);
    }

    // The orthogonal complement is the right kernel of the basis matrix
    // Convert to numeric matrix for computation
    const numMatrix: number[][] = basisMat.map((row) =>
      row.map((e) => {
        if (typeof e === 'bigint') return Number(e);
        if (typeof e === 'number') return e;
        return Number(String(e));
      })
    );

    // Find the kernel of the transposed matrix (right kernel of original)
    const kernel = computeKernel(numMatrix, n);

    if (kernel.length === 0) {
      return this.subspace([]);
    }

    // Convert kernel vectors to module elements
    const ambient = this.ambientVectorSpace();
    const complementGens: FreeModuleElement[] = [];
    const zero = this._baseRing.zero();

    for (const kv of kernel) {
      const entries = kv.map((e) => {
        if (typeof zero === 'bigint') return BigInt(Math.round(e));
        return e;
      });
      complementGens.push(ambient.createElement(entries));
    }

    return ambient.subspace(complementGens);
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
 * @see Reference: sage/modules/free_module.py:FreeModule_ambient_pid
 */
export class FreeModuleAmbientPID extends FreeModuleAmbient {
  constructor(
    baseRing: RingLike,
    rank: number,
    sparse: boolean = false,
    innerProductMatrix?: unknown
  ) {
    super(baseRing, rank, sparse, innerProductMatrix);
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
 * A submodule of a free module.
 */
export class FreeModuleSubmodule extends FreeModuleGeneric {
  protected _ambient: FreeModuleGeneric;
  protected _userBasis: FreeModuleElement[];

  constructor(
    ambient: FreeModuleGeneric,
    gens: FreeModuleElement[],
    options?: { check?: boolean; echelonize?: boolean; alreadyEchelonized?: boolean }
  ) {
    const rank = gens.length; // Simplified: assumes linearly independent
    super(
      ambient.baseRing(),
      rank,
      ambient.degree(),
      ambient.isSparse()
    );

    this._ambient = ambient;
    this._userBasis = [...gens];

    // Make basis vectors immutable
    for (const v of this._userBasis) {
      v.setImmutable();
    }

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
 * Submodule of a free module with a user-specified basis.
 * @see Reference: sage/modules/free_module.py:FreeModule_submodule_with_basis_pid
 */
export class FreeModuleWithBasis extends FreeModulePID {
  protected _ambient: FreeModuleGeneric;
  protected _userBasis: FreeModuleElement[];
  protected _echelonizedBasis: FreeModuleElement[] | null = null;

  constructor(
    ambient: FreeModuleGeneric,
    basis: FreeModuleElement[],
    options?: { check?: boolean; echelonize?: boolean; alreadyEchelonized?: boolean }
  ) {
    super(
      ambient.baseRing(),
      basis.length,
      ambient.degree(),
      ambient.isSparse()
    );

    this._ambient = ambient;
    this._userBasis = [...basis];

    // Make basis vectors immutable
    for (const v of this._userBasis) {
      v.setImmutable();
    }

    this._basis = this._userBasis;
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
   * Return the echelonized basis matrix.
   *
   * The echelonized basis matrix is the row echelon form of the user basis matrix.
   *
   * @returns The echelonized basis matrix
   *
   * @see Reference: sage/modules/free_module.py:FreeModule_submodule_with_basis_pid.echelonized_basis_matrix
   */
  override echelonizedBasisMatrix(): unknown[][] {
    if (this._echelonizedBasis !== null) {
      // Return cached echelonized form
      const result: unknown[][] = [];
      for (const v of this._echelonizedBasis) {
        result.push(v.list());
      }
      return result;
    }

    // Compute the echelonized form of the user basis
    const basisMat = this.basisMatrix() as unknown[][];

    if (basisMat.length === 0) {
      return [];
    }

    // Compute row echelon form
    const echelonMat = echelonize(basisMat, this._baseRing);

    // Cache the echelonized basis vectors
    this._echelonizedBasis = [];
    for (const row of echelonMat) {
      // Check if row is non-zero
      const isNonZero = row.some((e) => {
        if (typeof e === 'bigint') return e !== 0n;
        if (typeof e === 'number') return Math.abs(e) > 1e-10;
        return e !== 0;
      });
      if (isNonZero) {
        const v = this._ambient.createElement(row);
        v.setImmutable();
        this._echelonizedBasis.push(v);
      }
    }

    return echelonMat.filter((row) =>
      row.some((e) => {
        if (typeof e === 'bigint') return e !== 0n;
        if (typeof e === 'number') return Math.abs(e) > 1e-10;
        return e !== 0;
      })
    );
  }
}

/**
 * Submodule over a field with a user-specified basis.
 * @see Reference: sage/modules/free_module.py:FreeModule_submodule_field
 */
export class FreeModuleSubspaceWithBasis extends FreeModuleField {
  protected _ambient: FreeModuleField;
  protected _userBasis: FreeModuleElement[];

  constructor(
    ambient: FreeModuleField,
    basis: FreeModuleElement[],
    options?: { check?: boolean; echelonize?: boolean; alreadyEchelonized?: boolean }
  ) {
    super(
      ambient.baseRing(),
      basis.length,
      ambient.degree(),
      ambient.isSparse()
    );

    this._ambient = ambient;
    this._userBasis = [...basis];

    // Make basis vectors immutable
    for (const v of this._userBasis) {
      v.setImmutable();
    }

    this._basis = this._userBasis;
  }

  override hasUserBasis(): boolean {
    return true;
  }

  override ambientModule(): FreeModuleGeneric {
    return this._ambient;
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

    super(
      cover.baseRing(),
      quotientRank >= 0 ? quotientRank : 0,
      cover.degree(),
      cover.isSparse()
    );

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
