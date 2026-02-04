/**
 * @module sage/matrix/matrix_special
 * @description Special matrix constructors and algorithms
 *
 * Port of: sage/matrix/special.py, sage/matrix/berlekamp_massey.py
 */

import { NotImplementedError, ValueError } from '../errors.js';
import type { CoefficientRing, RingElement } from '../rings/polynomial/polynomial_element.js';
import { Polynomial } from '../rings/polynomial/polynomial_element.js';
import { PolynomialRing } from '../rings/polynomial/polynomial_ring.js';
import { Matrix } from './matrix_generic.js';

// ============================================================================
// Random Matrices
// ============================================================================

/**
 * Return a random matrix over the given ring.
 *
 * Creates a matrix with random entries from the given ring. For finite rings,
 * entries are uniformly distributed. For infinite rings, the distribution
 * depends on the ring's random element generation.
 *
 * @param ring - The base ring (must support random element generation)
 * @param nrows - Number of rows
 * @param ncols - Number of columns (default: nrows)
 * @param algorithm - Algorithm: 'randomize' (default), 'echelonizable', 'unimodular', 'diagonalizable'
 * @param implementation - Matrix implementation to use (not yet used)
 * @param options - Additional options (density for sparse matrices, etc.)
 * @returns A random matrix
 * @see Reference: sage/matrix/special.py:random_matrix
 */
export function random_matrix<R extends RingElement>(
  ring: CoefficientRing<R>,
  nrows: number,
  ncols?: number,
  algorithm: 'randomize' | 'echelonizable' | 'unimodular' | 'diagonalizable' = 'randomize',
  implementation?: string,
  options?: { density?: number; seed?: number }
): Matrix<R> {
  const actualNcols = ncols ?? nrows;

  if (nrows < 0 || actualNcols < 0) {
    throw new ValueError('matrix dimensions must be non-negative');
  }

  // Check if ring has a random element method
  const ringWithRandom = ring as unknown as { random_element?: () => R };
  if (typeof ringWithRandom.random_element !== 'function') {
    throw new NotImplementedError('random_matrix requires a ring with random_element() method');
  }

  if (algorithm === 'randomize') {
    // Simple random matrix: fill with random elements
    const entries: R[][] = [];
    const density = options?.density ?? 1.0;

    for (let i = 0; i < nrows; i++) {
      entries.push([]);
      for (let j = 0; j < actualNcols; j++) {
        if (Math.random() < density) {
          entries[i]!.push(ringWithRandom.random_element());
        } else {
          entries[i]!.push(ring.zero());
        }
      }
    }

    return new Matrix(ring, nrows, actualNcols, entries);
  }

  if (algorithm === 'unimodular') {
    return random_unimodular_matrix(ring, nrows);
  }

  if (algorithm === 'diagonalizable') {
    return random_diagonalizable_matrix(ring, nrows);
  }

  if (algorithm === 'echelonizable') {
    return random_echelonizable_matrix(ring, nrows, actualNcols);
  }

  throw new ValueError(`Unknown algorithm: ${algorithm}`);
}

/**
 * Return a random matrix in reduced row echelon form.
 *
 * Creates a random matrix that is already in RREF form with the specified
 * number of pivots (rank). The pivot positions are chosen randomly.
 *
 * @param ring - The base ring
 * @param nrows - Number of rows
 * @param ncols - Number of columns
 * @param num_pivots - Number of pivot columns (rank)
 * @returns A random RREF matrix
 * @see Reference: sage/matrix/special.py:random_rref_matrix
 */
export function random_rref_matrix<R extends RingElement>(
  ring: CoefficientRing<R>,
  nrows: number,
  ncols: number,
  num_pivots?: number
): Matrix<R> {
  if (nrows < 0 || ncols < 0) {
    throw new ValueError('matrix dimensions must be non-negative');
  }

  const maxPivots = Math.min(nrows, ncols);
  const actualNumPivots = num_pivots ?? maxPivots;

  if (actualNumPivots < 0 || actualNumPivots > maxPivots) {
    throw new ValueError(`num_pivots must be between 0 and min(nrows, ncols) = ${maxPivots}`);
  }

  // Check if ring has a random element method
  const ringWithRandom = ring as unknown as { random_element?: () => R };
  if (typeof ringWithRandom.random_element !== 'function') {
    throw new NotImplementedError(
      'random_rref_matrix requires a ring with random_element() method'
    );
  }

  const result = new Matrix<R>(ring, nrows, ncols);

  // Choose random pivot columns
  const allCols: number[] = [];
  for (let j = 0; j < ncols; j++) {
    allCols.push(j);
  }
  // Fisher-Yates shuffle
  for (let i = allCols.length - 1; i > 0; i--) {
    const j = Math.floor(Math.random() * (i + 1));
    [allCols[i], allCols[j]] = [allCols[j]!, allCols[i]!];
  }
  const pivotCols = allCols.slice(0, actualNumPivots).sort((a, b) => a - b);

  // Build RREF matrix
  for (let i = 0; i < actualNumPivots; i++) {
    const pivotCol = pivotCols[i]!;

    // Set pivot to 1
    result.set(i, pivotCol, ring.one());

    // Fill in random values in non-pivot columns to the right of this pivot
    for (let j = pivotCol + 1; j < ncols; j++) {
      if (!pivotCols.includes(j)) {
        result.set(i, j, ringWithRandom.random_element());
      }
    }
  }

  return result;
}

/**
 * Return a random echelonizable matrix of given rank.
 *
 * Creates a matrix that when reduced to echelon form will have the specified rank.
 * The matrix is constructed by creating a random matrix in reduced row echelon
 * form (RREF) and then applying random elementary row and column operations.
 *
 * @param ring - The base ring
 * @param nrows - Number of rows
 * @param ncols - Number of columns (default: nrows)
 * @param rank - The desired rank (default: min(nrows, ncols))
 * @param upper_bound - Upper bound on entry size (for integer rings, not used for finite fields)
 * @param max_tries - Maximum attempts (not currently used)
 * @returns A random echelonizable matrix
 * @see Reference: sage/matrix/special.py:random_echelonizable_matrix
 */
export function random_echelonizable_matrix<R extends RingElement>(
  ring: CoefficientRing<R>,
  nrows: number,
  ncols?: number,
  rank?: number,
  upper_bound?: number,
  max_tries?: number
): Matrix<R> {
  const actualNcols = ncols ?? nrows;
  const actualRank = rank ?? Math.min(nrows, actualNcols);

  if (nrows < 0 || actualNcols < 0) {
    throw new ValueError('matrix dimensions must be non-negative');
  }

  if (actualRank < 0 || actualRank > Math.min(nrows, actualNcols)) {
    throw new ValueError(
      `rank must be between 0 and min(nrows, ncols) = ${Math.min(nrows, actualNcols)}`
    );
  }

  // Check if ring has a random element method
  const ringWithRandom = ring as unknown as { random_element?: () => R };
  if (typeof ringWithRandom.random_element !== 'function') {
    throw new NotImplementedError(
      'random_echelonizable_matrix requires a ring with random_element() method'
    );
  }

  // Start with a random RREF matrix of the given rank
  // RREF has 1's at pivot positions, 0's below and above pivots, and can have
  // arbitrary values in non-pivot columns

  const result = new Matrix<R>(ring, nrows, actualNcols);

  // Choose random pivot columns
  const allCols: number[] = [];
  for (let j = 0; j < actualNcols; j++) {
    allCols.push(j);
  }
  // Shuffle and take first actualRank columns as pivots
  for (let i = allCols.length - 1; i > 0; i--) {
    const j = Math.floor(Math.random() * (i + 1));
    [allCols[i], allCols[j]] = [allCols[j]!, allCols[i]!];
  }
  const pivotCols = allCols.slice(0, actualRank).sort((a, b) => a - b);

  // Fill in RREF structure
  for (let i = 0; i < actualRank; i++) {
    const pivotCol = pivotCols[i]!;

    // Set pivot to 1
    result.set(i, pivotCol, ring.one());

    // Fill in random values in non-pivot columns to the right
    for (let j = pivotCol + 1; j < actualNcols; j++) {
      if (!pivotCols.includes(j)) {
        result.set(i, j, ringWithRandom.random_element());
      }
    }
  }

  // Apply random elementary row operations to scramble the matrix
  // while preserving its rank
  const numOperations = (nrows + actualNcols) * 2;

  for (let op = 0; op < numOperations; op++) {
    const opType = Math.floor(Math.random() * 2);
    const i = Math.floor(Math.random() * nrows);
    let j = Math.floor(Math.random() * nrows);
    while (j === i && nrows > 1) {
      j = Math.floor(Math.random() * nrows);
    }

    if (opType === 0 && nrows > 1) {
      // Row swap
      for (let k = 0; k < actualNcols; k++) {
        const tmp = result.get(i, k);
        result.set(i, k, result.get(j, k));
        result.set(j, k, tmp);
      }
    } else if (nrows > 1) {
      // Row addition
      const scale = ringWithRandom.random_element();
      for (let k = 0; k < actualNcols; k++) {
        result.set(i, k, result.get(i, k).add(scale.mul(result.get(j, k))) as R);
      }
    }
  }

  return result;
}

/**
 * Return a random matrix whose row and column spaces have given rank.
 *
 * Creates a matrix with predictable subspace structure useful for teaching.
 *
 * @param ring - The base ring
 * @param nrows - Number of rows
 * @param ncols - Number of columns
 * @param rank - The desired rank (default: min(nrows, ncols))
 * @returns A random subspaces matrix
 * @see Reference: sage/matrix/special.py:random_subspaces_matrix
 */
export function random_subspaces_matrix<R extends RingElement>(
  ring: CoefficientRing<R>,
  nrows: number,
  ncols: number,
  rank?: number
): Matrix<R> {
  // This is similar to random_echelonizable_matrix but designed
  // to have nice integer entries when computing subspaces
  return random_echelonizable_matrix(ring, nrows, ncols, rank);
}

/**
 * Return a random unimodular (det = +/- 1) matrix over a ring.
 *
 * A unimodular matrix is an invertible square matrix with determinant +/- 1
 * (or a unit in the base ring). For finite fields, any invertible matrix
 * can be considered "unimodular" since all non-zero elements are units.
 *
 * The algorithm creates a product of random elementary matrices, which
 * guarantees the result is unimodular.
 *
 * @param ring - The base ring
 * @param n - Size of the matrix
 * @param upper_bound - Upper bound on entry size (for integer rings, not used for finite fields)
 * @param max_tries - Maximum attempts (not currently used)
 * @returns A random unimodular matrix
 * @see Reference: sage/matrix/special.py:random_unimodular_matrix
 */
export function random_unimodular_matrix<R extends RingElement>(
  ring: CoefficientRing<R>,
  n: number,
  upper_bound?: number,
  max_tries?: number
): Matrix<R> {
  if (n < 0) {
    throw new ValueError('matrix size must be non-negative');
  }

  if (n === 0) {
    return new Matrix(ring, 0, 0);
  }

  // Check if ring has a random element method
  const ringWithRandom = ring as unknown as { random_element?: () => R };
  if (typeof ringWithRandom.random_element !== 'function') {
    throw new NotImplementedError(
      'random_unimodular_matrix requires a ring with random_element() method'
    );
  }

  // Start with identity matrix
  const result = new Matrix<R>(ring, n, n);
  for (let i = 0; i < n; i++) {
    result.set(i, i, ring.one());
  }

  // Apply random elementary operations
  // Each elementary operation has determinant 1 (row/col addition) or -1 (row/col swap)
  // or a unit (row/col scaling by unit)

  const numOperations = n * 2 + Math.floor(Math.random() * n);

  for (let op = 0; op < numOperations; op++) {
    const opType = Math.floor(Math.random() * 3);
    const i = Math.floor(Math.random() * n);
    let j = Math.floor(Math.random() * n);
    while (j === i && n > 1) {
      j = Math.floor(Math.random() * n);
    }

    if (opType === 0 && n > 1) {
      // Row swap: swap rows i and j (det *= -1)
      for (let k = 0; k < n; k++) {
        const tmp = result.get(i, k);
        result.set(i, k, result.get(j, k));
        result.set(j, k, tmp);
      }
    } else if (opType === 1 && n > 1) {
      // Row addition: add multiple of row j to row i (det unchanged)
      const scale = ringWithRandom.random_element();
      for (let k = 0; k < n; k++) {
        result.set(i, k, result.get(i, k).add(scale.mul(result.get(j, k))) as R);
      }
    } else {
      // Row scaling by a unit element
      // For finite fields, any non-zero element is a unit
      let unit = ringWithRandom.random_element();
      while (unit.isZero()) {
        unit = ringWithRandom.random_element();
      }
      for (let k = 0; k < n; k++) {
        result.set(i, k, result.get(i, k).mul(unit) as R);
      }
    }
  }

  return result;
}

/**
 * Return a random unitary matrix.
 *
 * A unitary matrix U satisfies U * U^* = I where U^* is the conjugate transpose.
 * For real matrices, this is the same as an orthogonal matrix (U * U^T = I).
 *
 * Note: For finite fields, the concept of unitary matrices doesn't directly apply.
 * This function requires rings with complex conjugation and square roots.
 *
 * @param ring - The base ring (must be RDF, CDF, or similar)
 * @param n - Size of the matrix
 * @returns A random unitary matrix
 * @see Reference: sage/matrix/special.py:random_unitary_matrix
 */
export function random_unitary_matrix<R extends RingElement>(
  ring: CoefficientRing<R>,
  n: number
): Matrix<R> {
  // For a proper implementation, we would use QR decomposition of a random matrix
  // or Haar measure sampling. This requires sqrt and is complex for general rings.
  throw new NotImplementedError(
    'random_unitary_matrix requires rings with sqrt (e.g., RDF, CDF). ' +
      'For finite fields, consider random_unimodular_matrix instead.'
  );
}

/**
 * Return a random diagonalizable matrix.
 *
 * Creates a diagonalizable matrix with specified or random eigenvalues.
 * The matrix is constructed as P * D * P^{-1} where D is diagonal
 * and P is a random invertible matrix.
 *
 * @param ring - The base ring
 * @param n - Size of the matrix
 * @param eigenvalues - Desired eigenvalues (if not provided, uses random elements)
 * @param dimensions - Dimensions of eigenspaces (each eigenvalue's multiplicity)
 * @returns A random diagonalizable matrix
 * @see Reference: sage/matrix/special.py:random_diagonalizable_matrix
 */
export function random_diagonalizable_matrix<R extends RingElement>(
  ring: CoefficientRing<R>,
  n: number,
  eigenvalues?: R[],
  dimensions?: number[]
): Matrix<R> {
  if (n < 0) {
    throw new ValueError('matrix size must be non-negative');
  }

  if (n === 0) {
    return new Matrix(ring, 0, 0);
  }

  // Check if ring has a random element method
  const ringWithRandom = ring as unknown as { random_element?: () => R };
  if (typeof ringWithRandom.random_element !== 'function') {
    throw new NotImplementedError(
      'random_diagonalizable_matrix requires a ring with random_element() method'
    );
  }

  // Build diagonal matrix D
  const D = new Matrix<R>(ring, n, n);

  if (eigenvalues !== undefined && dimensions !== undefined) {
    // Use provided eigenvalues with specified multiplicities
    if (eigenvalues.length !== dimensions.length) {
      throw new ValueError('eigenvalues and dimensions must have the same length');
    }
    const totalDim = dimensions.reduce((a, b) => a + b, 0);
    if (totalDim !== n) {
      throw new ValueError(`sum of dimensions (${totalDim}) must equal n (${n})`);
    }

    let idx = 0;
    for (let i = 0; i < eigenvalues.length; i++) {
      for (let j = 0; j < dimensions[i]!; j++) {
        D.set(idx, idx, eigenvalues[i]!);
        idx++;
      }
    }
  } else if (eigenvalues !== undefined) {
    // Use provided eigenvalues (one per diagonal entry)
    if (eigenvalues.length !== n) {
      throw new ValueError(`eigenvalues length (${eigenvalues.length}) must equal n (${n})`);
    }
    for (let i = 0; i < n; i++) {
      D.set(i, i, eigenvalues[i]!);
    }
  } else {
    // Generate random eigenvalues
    for (let i = 0; i < n; i++) {
      D.set(i, i, ringWithRandom.random_element());
    }
  }

  // Generate random invertible matrix P
  // For simplicity, we use a random unimodular matrix
  const P = random_unimodular_matrix(ring, n);

  // Compute P * D * P^{-1}
  // For fields, we can compute the inverse of P
  // For now, we use a different approach: construct using elementary operations

  // P * D
  const PD = P.mul(D);

  // For the inverse of P, we need field operations
  // We'll check if we can use the inverse function
  const fieldRing = ring as unknown as {
    inverse?: (x: R) => R;
    one: () => R;
    zero: () => R;
  };

  if (
    typeof fieldRing.inverse !== 'function' &&
    typeof (ring.one() as unknown as { inverse?: () => R }).inverse !== 'function'
  ) {
    // Cannot compute inverse, just return D (still diagonalizable)
    return D;
  }

  // Compute P^{-1} using Gaussian elimination
  const Aug = new Matrix<R>(ring, n, 2 * n);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      Aug.set(i, j, P.get(i, j));
      Aug.set(i, n + j, i === j ? ring.one() : ring.zero());
    }
  }

  // Forward elimination with partial pivoting
  for (let col = 0; col < n; col++) {
    // Find pivot
    let pivotRow = -1;
    for (let i = col; i < n; i++) {
      if (!Aug.get(i, col).isZero()) {
        pivotRow = i;
        break;
      }
    }

    if (pivotRow === -1) {
      // P is not invertible (shouldn't happen for unimodular matrix)
      return D;
    }

    // Swap rows
    if (pivotRow !== col) {
      for (let j = 0; j < 2 * n; j++) {
        const tmp = Aug.get(col, j);
        Aug.set(col, j, Aug.get(pivotRow, j));
        Aug.set(pivotRow, j, tmp);
      }
    }

    // Scale pivot row
    const pivot = Aug.get(col, col);
    const pivotInv = (pivot as unknown as { inverse?: () => R }).inverse
      ? (pivot as unknown as { inverse: () => R }).inverse()
      : (pivot as unknown as { inv?: () => R }).inv
        ? (pivot as unknown as { inv: () => R }).inv()
        : null;

    if (pivotInv === null) {
      return D;
    }

    for (let j = col; j < 2 * n; j++) {
      Aug.set(col, j, Aug.get(col, j).mul(pivotInv) as R);
    }

    // Eliminate
    for (let i = 0; i < n; i++) {
      if (i !== col && !Aug.get(i, col).isZero()) {
        const factor = Aug.get(i, col);
        for (let j = col; j < 2 * n; j++) {
          Aug.set(i, j, Aug.get(i, j).sub(factor.mul(Aug.get(col, j)) as R) as R);
        }
      }
    }
  }

  // Extract P^{-1}
  const PInv = new Matrix<R>(ring, n, n);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      PInv.set(i, j, Aug.get(i, n + j));
    }
  }

  // Return P * D * P^{-1}
  return PD.mul(PInv);
}

// ============================================================================
// Special Matrix Constructors
// ============================================================================

/**
 * Return a column matrix from a list of column vectors.
 *
 * Given a list of vectors (or lists), creates a matrix where each input
 * becomes a column.
 *
 * @param ring - The base ring
 * @param columns - Array of column vectors
 * @returns A matrix with the given columns
 * @see Reference: sage/matrix/special.py:column_matrix
 */
export function column_matrix<R extends RingElement>(
  ring: CoefficientRing<R>,
  columns: R[][]
): Matrix<R> {
  if (columns.length === 0) {
    return new Matrix(ring, 0, 0);
  }

  const ncols = columns.length;
  const nrows = columns[0]!.length;

  // Verify all columns have the same length
  for (let j = 1; j < ncols; j++) {
    if (columns[j]!.length !== nrows) {
      throw new ValueError(
        `all columns must have the same length, got ${nrows} and ${columns[j]!.length}`
      );
    }
  }

  const entries: R[][] = [];
  for (let i = 0; i < nrows; i++) {
    entries.push([]);
    for (let j = 0; j < ncols; j++) {
      entries[i]!.push(columns[j]![i]!);
    }
  }

  return new Matrix(ring, nrows, ncols, entries);
}

/**
 * Return a ones matrix.
 *
 * @param ring - The base ring
 * @param nrows - Number of rows
 * @param ncols - Number of columns (default: nrows)
 * @param sparse - Whether to use sparse representation (not used in dense impl)
 * @returns A matrix of all ones
 * @see Reference: sage/matrix/special.py:ones_matrix
 */
export function ones_matrix<R extends RingElement>(
  ring: CoefficientRing<R>,
  nrows: number = 1,
  ncols?: number,
  sparse?: boolean
): Matrix<R> {
  const actualNcols = ncols ?? nrows;
  const one = ring.one();

  const entries: R[][] = [];
  for (let i = 0; i < nrows; i++) {
    entries.push([]);
    for (let j = 0; j < actualNcols; j++) {
      entries[i]!.push(one);
    }
  }

  return new Matrix(ring, nrows, actualNcols, entries);
}

/**
 * Return a Lehmer matrix.
 *
 * The Lehmer matrix is an n x n matrix where element (i, j) is
 * min(i+1, j+1) / max(i+1, j+1), using 1-based indexing for the formula.
 *
 * @param ring - The base ring (must support division, e.g., QQ or a field)
 * @param n - Size of the matrix
 * @returns The Lehmer matrix
 * @see Reference: sage/matrix/special.py:lehmer
 */
export function lehmer<R extends RingElement>(ring: CoefficientRing<R>, n?: number): Matrix<R> {
  // Handle case where first argument is actually n (for convenience like SageMath)
  let size = n;
  const actualRing = ring;

  if (typeof (ring as unknown) === 'number' || typeof (ring as unknown) === 'bigint') {
    size = Number(ring as unknown);
    // We need a ring that supports fractions
    throw new NotImplementedError(
      'lehmer requires explicit ring argument when called with just dimension'
    );
  }

  if (size === undefined || size < 0) {
    throw new ValueError('n must be a non-negative integer');
  }

  if (size === 0) {
    return new Matrix(actualRing, 0, 0);
  }

  // Check if ring supports division
  const ringWithDiv = actualRing as unknown as {
    __call__?: (x: number | bigint | unknown) => R;
    div?: (a: R, b: R) => R;
  };

  if (typeof ringWithDiv.__call__ !== 'function') {
    throw new NotImplementedError('lehmer requires a ring with __call__ method');
  }

  const entries: R[][] = [];
  for (let i = 0; i < size; i++) {
    entries.push([]);
    for (let j = 0; j < size; j++) {
      // Lehmer(i,j) = min(i+1, j+1) / max(i+1, j+1)
      const minVal = Math.min(i + 1, j + 1);
      const maxVal = Math.max(i + 1, j + 1);

      // Try to create the fraction in the ring
      // This works for rings like QQ that can construct from ratios
      const numerator = ringWithDiv.__call__(minVal);
      const denominator = ringWithDiv.__call__(maxVal);

      // Compute numerator / denominator
      const denom = denominator as unknown as { inverse?: () => R; inv?: () => R };
      if (typeof denom.inverse === 'function') {
        entries[i]!.push(numerator.mul(denom.inverse()) as R);
      } else if (typeof denom.inv === 'function') {
        entries[i]!.push(numerator.mul(denom.inv()) as R);
      } else {
        throw new NotImplementedError('lehmer requires a ring with division');
      }
    }
  }

  return new Matrix(actualRing, size, size, entries);
}

/**
 * Return an elementary matrix.
 *
 * Elementary matrices correspond to elementary row operations:
 * - Row swap: swap rows row1 and row2
 * - Row scaling: multiply row1 by scale
 * - Row addition: add scale * row2 to row1
 *
 * @param ring - The base ring
 * @param n - Size of the matrix
 * @param options - Operation specification
 * @returns An elementary matrix
 * @see Reference: sage/matrix/special.py:elementary_matrix
 */
export function elementary_matrix<R extends RingElement>(
  ring: CoefficientRing<R>,
  n: number,
  options: { row1: number; row2?: number; scale?: R }
): Matrix<R> {
  if (n <= 0) {
    throw new ValueError('n must be positive');
  }

  const { row1, row2, scale } = options;

  if (row1 < 0 || row1 >= n) {
    throw new ValueError('row1 out of bounds');
  }

  // Start with identity matrix
  const entries: R[][] = [];
  for (let i = 0; i < n; i++) {
    entries.push([]);
    for (let j = 0; j < n; j++) {
      entries[i]!.push(i === j ? ring.one() : ring.zero());
    }
  }

  if (row2 !== undefined) {
    if (row2 < 0 || row2 >= n) {
      throw new ValueError('row2 out of bounds');
    }

    if (scale !== undefined) {
      // Row addition: add scale * row2 to row1
      // E[row1][row2] = scale (if row1 != row2)
      if (row1 !== row2) {
        entries[row1]![row2] = scale;
      }
    } else {
      // Row swap: swap rows row1 and row2
      // E[row1][row1] = 0, E[row1][row2] = 1
      // E[row2][row2] = 0, E[row2][row1] = 1
      entries[row1]![row1] = ring.zero();
      entries[row1]![row2] = ring.one();
      entries[row2]![row2] = ring.zero();
      entries[row2]![row1] = ring.one();
    }
  } else if (scale !== undefined) {
    // Row scaling: multiply row1 by scale
    entries[row1]![row1] = scale;
  } else {
    throw new ValueError('must specify either row2 or scale');
  }

  return new Matrix(ring, n, n, entries);
}

/**
 * Return a circulant matrix.
 *
 * A circulant matrix is defined by its first row, and each subsequent row
 * is a cyclic shift of the previous row.
 *
 * @param ring - The base ring
 * @param v - First row
 * @param sparse - Whether to use sparse representation (not used in dense impl)
 * @returns A circulant matrix
 * @see Reference: sage/matrix/special.py:circulant
 */
export function circulant<R extends RingElement>(
  ring: CoefficientRing<R>,
  v: R[],
  sparse?: boolean
): Matrix<R> {
  const n = v.length;
  if (n === 0) {
    return new Matrix(ring, 0, 0);
  }

  const entries: R[][] = [];
  for (let i = 0; i < n; i++) {
    entries.push([]);
    for (let j = 0; j < n; j++) {
      // Entry (i, j) is v[(j - i) mod n]
      const idx = (((j - i) % n) + n) % n;
      entries[i]!.push(v[idx]!);
    }
  }

  return new Matrix(ring, n, n, entries);
}

/**
 * Return a block matrix.
 *
 * Constructs a larger matrix by concatenating submatrices in a grid.
 * For example, block_matrix([[A, B], [C, D]]) creates:
 *   [ A  B ]
 *   [ C  D ]
 *
 * @param ring - The base ring
 * @param blocks - 2D array of blocks, where each block is a Matrix or a scalar
 * @param options - Optional configuration: { subdivide?: boolean }
 * @returns A block matrix
 * @see Reference: sage/matrix/special.py:block_matrix
 */
export function block_matrix<R extends RingElement>(
  ring: CoefficientRing<R>,
  blocks: Array<Array<Matrix<R> | R | number | 0>>,
  options?: { subdivide?: boolean }
): Matrix<R> {
  if (blocks.length === 0) {
    return new Matrix(ring, 0, 0);
  }

  const nBlockRows = blocks.length;
  const nBlockCols = blocks[0]!.length;

  // Validate that all rows have the same number of blocks
  for (let i = 1; i < nBlockRows; i++) {
    if (blocks[i]!.length !== nBlockCols) {
      throw new ValueError(
        `all block rows must have the same number of columns, got ${nBlockCols} and ${blocks[i]!.length}`
      );
    }
  }

  // Determine row heights and column widths from the actual matrices
  const rowHeights: number[] = new Array(nBlockRows).fill(0);
  const colWidths: number[] = new Array(nBlockCols).fill(0);

  // First pass: determine dimensions from non-scalar blocks
  for (let i = 0; i < nBlockRows; i++) {
    for (let j = 0; j < nBlockCols; j++) {
      const block = blocks[i]![j];
      if (block instanceof Matrix) {
        if (rowHeights[i] === 0) {
          rowHeights[i] = block.nrows;
        } else if (rowHeights[i] !== block.nrows) {
          throw new ValueError(
            `incompatible block heights in row ${i}: ${rowHeights[i]} vs ${block.nrows}`
          );
        }
        if (colWidths[j] === 0) {
          colWidths[j] = block.ncols;
        } else if (colWidths[j] !== block.ncols) {
          throw new ValueError(
            `incompatible block widths in column ${j}: ${colWidths[j]} vs ${block.ncols}`
          );
        }
      }
    }
  }

  // Second pass: for scalar blocks, they become square matrices matching the dimensions
  // If dimension is still 0, use 1x1
  for (let i = 0; i < nBlockRows; i++) {
    if (rowHeights[i] === 0) {
      // Find the column width for this block row
      for (let j = 0; j < nBlockCols; j++) {
        if (colWidths[j] !== 0) {
          rowHeights[i] = colWidths[j];
          break;
        }
      }
      if (rowHeights[i] === 0) {
        rowHeights[i] = 1;
      }
    }
  }
  for (let j = 0; j < nBlockCols; j++) {
    if (colWidths[j] === 0) {
      // Find the row height for this block column
      for (let i = 0; i < nBlockRows; i++) {
        if (rowHeights[i] !== 0) {
          colWidths[j] = rowHeights[i];
          break;
        }
      }
      if (colWidths[j] === 0) {
        colWidths[j] = 1;
      }
    }
  }

  // Calculate total dimensions
  const totalRows = rowHeights.reduce((a, b) => a + b, 0);
  const totalCols = colWidths.reduce((a, b) => a + b, 0);

  // Create the result matrix
  const result = new Matrix<R>(ring, totalRows, totalCols);

  // Fill in the blocks
  let rowOffset = 0;
  for (let i = 0; i < nBlockRows; i++) {
    let colOffset = 0;
    for (let j = 0; j < nBlockCols; j++) {
      const block = blocks[i]![j];
      const blockHeight = rowHeights[i]!;
      const blockWidth = colWidths[j]!;

      if (block instanceof Matrix) {
        // Copy matrix block
        for (let bi = 0; bi < blockHeight; bi++) {
          for (let bj = 0; bj < blockWidth; bj++) {
            result.set(rowOffset + bi, colOffset + bj, block.get(bi, bj));
          }
        }
      } else if (block === 0 || (block as R)?.isZero?.()) {
        // Zero block - already initialized to zeros
      } else {
        // Scalar block - create scalar matrix (scalar * I)
        const scalar =
          typeof block === 'number'
            ? (ring as unknown as { __call__: (x: number) => R }).__call__(block)
            : (block as R);
        const size = Math.min(blockHeight, blockWidth);
        for (let k = 0; k < size; k++) {
          result.set(rowOffset + k, colOffset + k, scalar);
        }
      }

      colOffset += blockWidth;
    }
    rowOffset += rowHeights[i]!;
  }

  return result;
}

/**
 * Return a block diagonal matrix.
 *
 * Constructs a matrix with the given matrices along the diagonal
 * and zeros elsewhere.
 *
 * @param ring - The base ring
 * @param sub_matrices - The diagonal blocks
 * @returns A block diagonal matrix
 * @see Reference: sage/matrix/special.py:block_diagonal_matrix
 */
export function block_diagonal_matrix<R extends RingElement>(
  ring: CoefficientRing<R>,
  ...sub_matrices: Matrix<R>[]
): Matrix<R> {
  if (sub_matrices.length === 0) {
    return new Matrix(ring, 0, 0);
  }

  // Calculate total dimensions
  let totalRows = 0;
  let totalCols = 0;
  for (const M of sub_matrices) {
    totalRows += M.nrows;
    totalCols += M.ncols;
  }

  const result = new Matrix<R>(ring, totalRows, totalCols);

  // Copy each block to the diagonal
  let rowOffset = 0;
  let colOffset = 0;
  for (const M of sub_matrices) {
    for (let i = 0; i < M.nrows; i++) {
      for (let j = 0; j < M.ncols; j++) {
        result.set(rowOffset + i, colOffset + j, M.get(i, j));
      }
    }
    rowOffset += M.nrows;
    colOffset += M.ncols;
  }

  return result;
}

/**
 * Return a Jordan block matrix.
 *
 * A Jordan block of size n with eigenvalue lambda is:
 *   [ lambda  1      0    ...  0   ]
 *   [ 0      lambda  1    ...  0   ]
 *   [ 0       0     lambda ... 0   ]
 *   [ ...    ...    ...   ... ...  ]
 *   [ 0       0      0    ... lambda]
 *
 * @param ring - The base ring
 * @param eigenvalue - The eigenvalue
 * @param size - Size of the block
 * @param sparse - Whether to use sparse representation (not used in dense impl)
 * @returns A Jordan block
 * @see Reference: sage/matrix/special.py:jordan_block
 */
export function jordan_block<R extends RingElement>(
  ring: CoefficientRing<R>,
  eigenvalue: R,
  size: number,
  sparse?: boolean
): Matrix<R> {
  if (size < 0) {
    throw new ValueError('size must be non-negative');
  }

  const one = ring.one();
  const result = new Matrix<R>(ring, size, size);

  // Set diagonal entries to eigenvalue
  for (let i = 0; i < size; i++) {
    result.set(i, i, eigenvalue);
  }

  // Set superdiagonal entries to 1
  for (let i = 0; i < size - 1; i++) {
    result.set(i, i + 1, one);
  }

  return result;
}

/**
 * Return a companion matrix for a monic polynomial.
 *
 * For a monic polynomial p(x) = x^n - c_{n-1}x^{n-1} - ... - c_1x - c_0,
 * the companion matrix has the property that its characteristic polynomial
 * is p(x).
 *
 * @param ring - The base ring
 * @param coeffs - Coefficients [c_0, c_1, ..., c_{n-1}] of the polynomial
 *                 (excluding the leading coefficient which must be 1)
 * @param format - 'right' (default) or 'left' companion form
 * @returns The companion matrix
 * @see Reference: sage/matrix/special.py:companion_matrix
 */
export function companion_matrix<R extends RingElement>(
  ring: CoefficientRing<R>,
  coeffs: R[],
  format: 'right' | 'left' = 'right'
): Matrix<R> {
  const n = coeffs.length;
  if (n === 0) {
    return new Matrix(ring, 0, 0);
  }

  const entries: R[][] = [];
  for (let i = 0; i < n; i++) {
    entries.push([]);
    for (let j = 0; j < n; j++) {
      entries[i]!.push(ring.zero());
    }
  }

  if (format === 'right') {
    // Right companion matrix:
    // [0 0 ... 0 c_0  ]
    // [1 0 ... 0 c_1  ]
    // [0 1 ... 0 c_2  ]
    // [...      ...   ]
    // [0 0 ... 1 c_n-1]
    for (let i = 1; i < n; i++) {
      entries[i]![i - 1] = ring.one();
    }
    for (let i = 0; i < n; i++) {
      entries[i]![n - 1] = coeffs[i]!;
    }
  } else {
    // Left companion matrix:
    // [c_{n-1} c_{n-2} ... c_1 c_0]
    // [1       0       ... 0   0  ]
    // [0       1       ... 0   0  ]
    // [...                        ]
    // [0       0       ... 1   0  ]
    for (let j = 0; j < n; j++) {
      entries[0]![j] = coeffs[n - 1 - j]!;
    }
    for (let i = 1; i < n; i++) {
      entries[i]![i - 1] = ring.one();
    }
  }

  return new Matrix(ring, n, n, entries);
}

/**
 * Return a Hilbert matrix.
 *
 * The Hilbert matrix has entry (i,j) = 1/(i+j+1), using 0-based indices.
 * This is equivalent to 1/(i+j-1) with 1-based indices as in the standard definition.
 *
 * The Hilbert matrix is a classic example of an ill-conditioned matrix.
 *
 * @param dim - Dimension of the matrix
 * @param ring - The base ring (must support fractions, e.g., QQ)
 * @returns The Hilbert matrix
 * @see Reference: sage/matrix/special.py:hilbert
 */
export function hilbert<R extends RingElement>(dim: number, ring: CoefficientRing<R>): Matrix<R> {
  if (dim < 0) {
    throw new ValueError('dimension must be non-negative');
  }

  if (dim === 0) {
    return new Matrix(ring, 0, 0);
  }

  // Check if ring supports division
  const ringWithDiv = ring as unknown as {
    __call__?: (x: number | bigint | unknown) => R;
  };

  if (typeof ringWithDiv.__call__ !== 'function') {
    throw new NotImplementedError('hilbert requires a ring with __call__ method');
  }

  const entries: R[][] = [];
  for (let i = 0; i < dim; i++) {
    entries.push([]);
    for (let j = 0; j < dim; j++) {
      // Hilbert(i,j) = 1 / (i + j + 1)
      // Using 0-based indexing
      const denomValue = i + j + 1;

      const one = ring.one();
      const denom = ringWithDiv.__call__(denomValue);

      // Compute 1 / denomValue
      const denomInv = denom as unknown as { inverse?: () => R; inv?: () => R };
      if (typeof denomInv.inverse === 'function') {
        entries[i]!.push(one.mul(denomInv.inverse()) as R);
      } else if (typeof denomInv.inv === 'function') {
        entries[i]!.push(one.mul(denomInv.inv()) as R);
      } else {
        throw new NotImplementedError('hilbert requires a ring with division');
      }
    }
  }

  return new Matrix(ring, dim, dim, entries);
}

/**
 * Return a Vandermonde matrix.
 *
 * The Vandermonde matrix for [a_0, a_1, ..., a_{n-1}] is:
 *   [ 1    a_0    a_0^2   ...  a_0^{n-1}  ]
 *   [ 1    a_1    a_1^2   ...  a_1^{n-1}  ]
 *   [ ...  ...    ...     ...  ...        ]
 *   [ 1    a_{n-1} a_{n-1}^2 ... a_{n-1}^{n-1} ]
 *
 * @param ring - The base ring
 * @param v - Vector of elements
 * @returns The Vandermonde matrix
 * @see Reference: sage/matrix/special.py:vandermonde
 */
export function vandermonde<R extends RingElement>(ring: CoefficientRing<R>, v: R[]): Matrix<R> {
  const n = v.length;
  if (n === 0) {
    return new Matrix(ring, 0, 0);
  }

  const one = ring.one();
  const entries: R[][] = [];

  for (let i = 0; i < n; i++) {
    entries.push([]);
    let power = one;
    for (let j = 0; j < n; j++) {
      entries[i]!.push(power);
      power = power.mul(v[i]!) as R;
    }
  }

  return new Matrix(ring, n, n, entries);
}

/**
 * Return a Toeplitz matrix.
 *
 * A Toeplitz matrix is constant along diagonals.
 * Entry (i, j) depends only on i - j.
 * c gives the first column, r gives the first row.
 * Note: c[0] should equal r[0].
 *
 * @param ring - The base ring
 * @param c - First column
 * @param r - First row
 * @returns The Toeplitz matrix
 * @see Reference: sage/matrix/special.py:toeplitz
 */
export function toeplitz<R extends RingElement>(
  ring: CoefficientRing<R>,
  c: R[],
  r: R[]
): Matrix<R> {
  const m = c.length;
  const n = r.length;

  if (m === 0 || n === 0) {
    return new Matrix(ring, m, n);
  }

  const entries: R[][] = [];
  for (let i = 0; i < m; i++) {
    entries.push([]);
    for (let j = 0; j < n; j++) {
      if (i >= j) {
        // Below or on diagonal: use column c
        entries[i]!.push(c[i - j]!);
      } else {
        // Above diagonal: use row r
        entries[i]!.push(r[j - i]!);
      }
    }
  }

  return new Matrix(ring, m, n, entries);
}

/**
 * Return a Hankel matrix.
 *
 * A Hankel matrix is constant along anti-diagonals.
 * Entry (i, j) depends only on i + j.
 * c gives the first column, r gives the last row.
 * If r is not provided, it's taken to be the same length as c with zeros.
 *
 * @param ring - The base ring
 * @param c - First column
 * @param r - Last row (optional)
 * @returns The Hankel matrix
 * @see Reference: sage/matrix/special.py:hankel
 */
export function hankel<R extends RingElement>(
  ring: CoefficientRing<R>,
  c: R[],
  r?: R[]
): Matrix<R> {
  const m = c.length;
  if (m === 0) {
    return new Matrix(ring, 0, 0);
  }

  // If r is not provided, use zeros for the last row except where it overlaps with c
  const actualR = r ?? [];
  const n = actualR.length > 0 ? actualR.length : m;

  const entries: R[][] = [];
  for (let i = 0; i < m; i++) {
    entries.push([]);
    for (let j = 0; j < n; j++) {
      const idx = i + j;
      if (idx < m) {
        entries[i]!.push(c[idx]!);
      } else if (actualR.length > 0 && idx - m + 1 < actualR.length) {
        entries[i]!.push(actualR[idx - m + 1]!);
      } else {
        entries[i]!.push(ring.zero());
      }
    }
  }

  return new Matrix(ring, m, n, entries);
}

/**
 * Return a rotation matrix that rotates a vector onto an axis.
 *
 * Returns an orthogonal matrix R such that R*v is a scalar multiple of e_i,
 * where e_i is the i-th standard basis vector.
 *
 * Note: This requires square roots and trigonometric functions, which are
 * not available for general rings. Works best with RDF or CDF.
 *
 * @param v - The vector
 * @param i - The axis index (0-based)
 * @param ring - The base ring
 * @returns The rotation matrix
 * @see Reference: sage/matrix/special.py:vector_on_axis_rotation_matrix
 */
export function vector_on_axis_rotation_matrix<R extends RingElement>(
  v: R[],
  i: number,
  ring?: CoefficientRing<R>
): Matrix<R> {
  // This requires sqrt and trigonometric functions
  throw new NotImplementedError(
    'vector_on_axis_rotation_matrix requires rings with sqrt (e.g., RDF, CDF)'
  );
}

/**
 * Return a rotation matrix that zeros the i-th coordinate.
 *
 * Returns an orthogonal matrix R such that (R*v)[i] = 0.
 * This is done using a Givens rotation.
 *
 * Note: This requires square roots and trigonometric functions, which are
 * not available for general rings. Works best with RDF or CDF.
 *
 * @param v - The vector
 * @param i - The coordinate to zero (0-based)
 * @param ring - The base ring
 * @returns The rotation matrix
 * @see Reference: sage/matrix/special.py:ith_to_zero_rotation_matrix
 */
export function ith_to_zero_rotation_matrix<R extends RingElement>(
  v: R[],
  i: number,
  ring?: CoefficientRing<R>
): Matrix<R> {
  // This requires sqrt and trigonometric functions
  throw new NotImplementedError(
    'ith_to_zero_rotation_matrix requires rings with sqrt (e.g., RDF, CDF)'
  );
}

// ============================================================================
// Berlekamp-Massey Algorithm
// ============================================================================

/**
 * Use the Berlekamp-Massey algorithm to find the minimal polynomial
 * of a linear recurrence sequence.
 *
 * The minimal polynomial of a linear recurrence {a_r} is the unique monic
 * polynomial g such that if {a_r} satisfies a linear recurrence
 * a_{j+k} + b_{j-1} * a_{j-1+k} + ... + b_0 * a_k = 0 (for all k >= 0),
 * then g divides x^j + sum_{i=0}^{j-1} b_i * x^i.
 *
 * @param a - List of even length of elements of a field
 * @returns The minimal polynomial of the sequence as a Polynomial
 * @see Reference: sage/matrix/berlekamp_massey.py:berlekamp_massey
 */
export function berlekamp_massey<R extends RingElement>(a: R[]): Polynomial<R> {
  if (a.length === 0) {
    throw new ValueError('sequence must be non-empty');
  }

  if (a.length % 2 !== 0) {
    throw new ValueError('argument must have an even number of terms');
  }

  const M = Math.floor(a.length / 2);

  if (M === 0) {
    throw new ValueError('sequence must have at least 2 elements');
  }

  // Get the field/ring from the first element
  const ring = a[0]!.parent() as CoefficientRing<R>;
  const polyRing = new PolynomialRing(ring, 'x');

  // Create polynomial from sequence: f0 = a[0] + a[1]*x + ... + a[2M-1]*x^{2M-1}
  // Note: We need to reverse since the polynomial class expects coeffs[i] = coeff of x^i
  const f0 = new Polynomial<R>(a, polyRing);

  // f1 = x^{2M}
  const f1Coeffs: R[] = [];
  for (let i = 0; i < 2 * M; i++) {
    f1Coeffs.push(ring.zero());
  }
  f1Coeffs.push(ring.one());
  let f1 = new Polynomial<R>(f1Coeffs, polyRing);
  let f0Current = f0;

  // s0 = 1, s1 = 0
  let s0 = new Polynomial<R>([ring.one()], polyRing);
  let s1 = new Polynomial<R>([ring.zero()], polyRing);

  // Extended Euclidean algorithm until degree < M
  while (f1.degree() >= M) {
    // Swap f0 and f1, compute quotient and remainder
    const temp = f1;
    const [q, r] = f0Current.quo_rem(f1);
    f0Current = temp;
    f1 = r;

    // Update s0, s1
    const newS1 = s0.sub(q.mul(s1)) as Polynomial<R>;
    s0 = s1;
    s1 = newS1;
  }

  // Return s1.reverse().monic()
  // reverse: if s1 = c_0 + c_1*x + ... + c_n*x^n, then reverse = c_n + c_{n-1}*x + ... + c_0*x^n
  const s1Coeffs = s1.coeffs.slice();
  while (s1Coeffs.length > 0 && s1Coeffs[s1Coeffs.length - 1]!.isZero()) {
    s1Coeffs.pop();
  }
  s1Coeffs.reverse();

  const reversed = new Polynomial<R>(s1Coeffs, polyRing);

  // Make monic
  return reversed.monic();
}

// ============================================================================
// Matrix Operations
// ============================================================================

/**
 * Return the tensor (Kronecker) product of two matrices.
 *
 * The tensor product A tensor B is the matrix with block structure where
 * the (i,j) block is A[i,j] * B.
 *
 * @param A - First matrix (m x n)
 * @param B - Second matrix (p x q)
 * @param subdivide - Whether to subdivide the result (not used in this impl)
 * @returns The tensor product A tensor B (mp x nq)
 * @see Reference: sage/matrix/matrix2.pyx:tensor_product
 */
export function tensor_product<R extends RingElement>(
  A: Matrix<R>,
  B: Matrix<R>,
  subdivide?: boolean
): Matrix<R> {
  const m = A.nrows;
  const n = A.ncols;
  const p = B.nrows;
  const q = B.ncols;
  const ring = A.base_ring;

  const result = new Matrix<R>(ring, m * p, n * q);

  for (let i = 0; i < m; i++) {
    for (let j = 0; j < n; j++) {
      const aij = A.get(i, j);
      for (let k = 0; k < p; k++) {
        for (let l = 0; l < q; l++) {
          result.set(i * p + k, j * q + l, aij.mul(B.get(k, l)) as R);
        }
      }
    }
  }

  return result;
}

/**
 * Return the elementwise (Hadamard) product of two matrices.
 *
 * The elementwise product is the matrix where (A * B)[i,j] = A[i,j] * B[i,j].
 *
 * @param A - First matrix
 * @param B - Second matrix (must have same dimensions as A)
 * @returns The elementwise product
 * @see Reference: sage/matrix/matrix2.pyx:elementwise_product
 */
export function elementwise_product<R extends RingElement>(A: Matrix<R>, B: Matrix<R>): Matrix<R> {
  if (A.nrows !== B.nrows || A.ncols !== B.ncols) {
    throw new ValueError(
      'matrices must have the same dimensions for elementwise product, ' +
        `got ${A.nrows}x${A.ncols} and ${B.nrows}x${B.ncols}`
    );
  }

  const ring = A.base_ring;
  const result = new Matrix<R>(ring, A.nrows, A.ncols);

  for (let i = 0; i < A.nrows; i++) {
    for (let j = 0; j < A.ncols; j++) {
      result.set(i, j, A.get(i, j).mul(B.get(i, j)) as R);
    }
  }

  return result;
}

/**
 * Return the rook vector of the matrix.
 *
 * The rook vector (r_0, r_1, ..., r_k) of a 0-1 matrix M is where r_i is
 * the number of ways to place i non-attacking rooks on the 1's of M.
 * More generally, r_i equals the i-th permanental minor of M.
 *
 * For a 0-1 matrix, the coefficient r_k counts the number of ways to place
 * k non-attacking rooks on positions where M has a 1.
 *
 * @param matrix - The matrix
 * @param algorithm - Algorithm to use: 'ButeraPernici' (default), 'Ryser', or 'naive'
 * @param complement - Whether to compute rook vector of complement matrix
 * @param use_complement - Force complement usage (auto-determined if not specified)
 * @returns The rook vector as array of ring elements
 * @see Reference: sage/matrix/matrix2.pyx:rook_vector
 * @see Reference: sage/matrix/matrix_misc.py:permanental_minor_polynomial
 */
export function rook_vector<R extends RingElement>(
  matrix: Matrix<R>,
  algorithm: 'ButeraPernici' | 'Ryser' | 'naive' = 'ButeraPernici',
  complement?: boolean,
  use_complement?: boolean
): R[] {
  const m = matrix.nrows;
  const n = matrix.ncols;
  const ring = matrix.base_ring;
  const mn = Math.min(m, n);
  const zero = ring.zero();
  const one = ring.one();

  // Check if matrix is 0-1 and count ones (needed for complement logic)
  let isZ2 = true;
  let numOnes = 0;
  outer: for (let i = 0; i < m; i++) {
    for (let j = 0; j < n; j++) {
      const x = matrix.get(i, j);
      if (!x.isZero()) {
        if (!x.eq(one)) {
          isZ2 = false;
          break outer;
        }
        numOnes++;
      }
    }
  }

  // Validate complement usage
  if (complement && !isZ2) {
    throw new ValueError('complement requires a 0-1 matrix');
  }

  // Auto-determine whether to use complement
  if (use_complement === undefined) {
    use_complement = isZ2 && numOnes > 0.55 * m * n;
  }

  let b: R[];

  if (use_complement) {
    // Compute rook vector of complement matrix
    const B = new Matrix<R>(ring, m, n);
    for (let i = 0; i < m; i++) {
      for (let j = 0; j < n; j++) {
        B.set(i, j, one.sub(matrix.get(i, j)) as R);
      }
    }
    b = rook_vector(B, algorithm, false, false);
    complement = !complement;
  } else if (algorithm === 'Ryser') {
    b = _rook_vector_ryser(matrix);
  } else if (algorithm === 'ButeraPernici') {
    b = _rook_vector_butera_pernici(matrix);
  } else if (algorithm === 'naive') {
    b = _rook_vector_naive(matrix);
  } else {
    throw new ValueError(`algorithm must be one of "ButeraPernici", "Ryser", or "naive", got "${algorithm}"`);
  }

  // Apply inclusion-exclusion if computing complement
  if (complement) {
    const a: R[] = [one];
    let c1 = 1;
    for (let k = 1; k <= mn; k++) {
      // c1 = C(m, k) * C(n, k) * k! / C(m, k-1) / C(n, k-1) / (k-1)!
      //    = (m-k+1) * (n-k+1) / k
      c1 = (c1 * (m - k + 1) * (n - k + 1)) / k;
      let c = c1;
      // s = c * b[0] + (-1)^k * b[k]
      let s = _scalarMul(ring, c, b[0]!);
      const sign_k = k % 2 === 0 ? 1 : -1;
      s = s.add(_scalarMul(ring, sign_k, b[k]!)) as R;

      for (let j = 1; j < k; j++) {
        // c = -c * (k-j+1) / ((m-j+1) * (n-j+1))
        c = (-c * (k - j + 1)) / ((m - j + 1) * (n - j + 1));
        s = s.add(_scalarMul(ring, c, b[j]!)) as R;
      }
      a.push(s);
    }
    return a;
  }

  return b;
}

/**
 * Compute rook vector using Ryser's algorithm for permanental minors.
 * Complexity: O(mn * 2^n) for each k, so O(mn^2 * 2^n) total.
 */
function _rook_vector_ryser<R extends RingElement>(matrix: Matrix<R>): R[] {
  const m = matrix.nrows;
  const n = matrix.ncols;
  const mn = Math.min(m, n);
  const ring = matrix.base_ring;

  const result: R[] = [ring.one()];
  for (let k = 1; k <= mn; k++) {
    result.push(_permanental_minor_ryser(matrix, k));
  }
  return result;
}

/**
 * Compute k-th permanental minor using Ryser's algorithm.
 */
function _permanental_minor_ryser<R extends RingElement>(matrix: Matrix<R>, k: number): R {
  const m = matrix.nrows;
  const n = matrix.ncols;
  const ring = matrix.base_ring;

  if (k === 0) return ring.one();
  if (k > m || k > n) return ring.zero();

  let pm = ring.zero();

  // Iterate over all k-subsets of rows and columns
  for (const rows of _choose_indices(m, k)) {
    for (const cols of _choose_indices(n, k)) {
      // Compute permanent of submatrix using Ryser's formula
      pm = pm.add(_permanent_submatrix_ryser(matrix, rows, cols)) as R;
    }
  }

  return pm;
}

/**
 * Compute permanent of a submatrix using Ryser's formula.
 * For a k x k matrix, complexity is O(k * 2^k).
 */
function _permanent_submatrix_ryser<R extends RingElement>(
  matrix: Matrix<R>,
  rows: number[],
  cols: number[]
): R {
  const k = rows.length;
  const ring = matrix.base_ring;

  if (k === 0) return ring.one();

  // Build the submatrix values for efficient access
  const submatrix: R[][] = [];
  for (let i = 0; i < k; i++) {
    submatrix.push([]);
    for (let j = 0; j < k; j++) {
      submatrix[i]!.push(matrix.get(rows[i]!, cols[j]!));
    }
  }

  // Ryser's formula: perm(A) = (-1)^n * sum_{S subset of {1..n}} (-1)^|S| * prod_i (sum_{j in S} a_ij)
  let perm = ring.zero();
  const numSubsets = 1 << k;

  for (let mask = 0; mask < numSubsets; mask++) {
    // Count bits in mask to get |S|
    let bitCount = 0;
    let temp = mask;
    while (temp) {
      bitCount += temp & 1;
      temp >>= 1;
    }

    // Compute product of row sums for selected columns
    let prod = ring.one();
    for (let i = 0; i < k; i++) {
      let rowSum = ring.zero();
      for (let j = 0; j < k; j++) {
        if (mask & (1 << j)) {
          rowSum = rowSum.add(submatrix[i]![j]!) as R;
        }
      }
      prod = prod.mul(rowSum) as R;
    }

    // Apply sign: (-1)^(k - |S|)
    const sign = (k - bitCount) % 2 === 0 ? 1 : -1;
    if (sign > 0) {
      perm = perm.add(prod) as R;
    } else {
      perm = perm.sub(prod) as R;
    }
  }

  return perm;
}

/**
 * Compute rook vector using the Butera-Pernici algorithm.
 *
 * This algorithm computes all permanental minors simultaneously by working
 * in a quotient polynomial ring where variables are nilpotent of order 2.
 *
 * Complexity: O(2^n * m * n) where m = nrows, n = ncols.
 *
 * @see [BP2015] P. Butera, M. Pernici, "Sums of permanental minors using Grassmann algebra"
 */
function _rook_vector_butera_pernici<R extends RingElement>(matrix: Matrix<R>): R[] {
  const m = matrix.nrows;
  const n = matrix.ncols;
  const mn = Math.min(m, n);
  const ring = matrix.base_ring;

  // Special case: empty matrix
  if (m === 0 || n === 0) {
    return [ring.one()];
  }

  // Get matrix rows for efficient access
  const rows: R[][] = [];
  for (let i = 0; i < m; i++) {
    const row: R[] = [];
    for (let j = 0; j < n; j++) {
      row.push(matrix.get(i, j));
    }
    rows.push(row);
  }

  // p is a dictionary mapping bitmasks to polynomial coefficients
  // The bitmask represents which eta_j variables are present
  // The value is an array of coefficients [c_0, c_1, ..., c_k] for polynomial c_0 + c_1*t + ... + c_k*t^k
  let p: Map<number, R[]> = new Map();
  p.set(0, [ring.one()]); // Start with 1

  // Track which columns still have nonzero entries in remaining rows
  const varsToDo = new Set<number>();
  for (let j = 0; j < n; j++) {
    varsToDo.add(j);
  }

  for (let i = 0; i < m; i++) {
    const a = rows[i]!;

    // Build p1 = 1 + t * sum_j A[i,j] * eta_j
    const p1: Map<number, R[]> = new Map();
    p1.set(0, [ring.one()]); // Constant term 1

    for (let j = 0; j < n; j++) {
      if (!a[j]!.isZero()) {
        // eta_j corresponds to bit j, coefficient is A[i,j] * t (degree 1)
        p1.set(1 << j, [ring.zero(), a[j]!]);
      }
    }

    // Determine which variables can be "integrated" (set to 1)
    // A variable eta_j can be integrated if it doesn't appear in any remaining row
    let maskFree = 0;
    const toRemove: number[] = [];
    for (const j of varsToDo) {
      let appearsLater = false;
      for (let k = i + 1; k < m; k++) {
        if (!rows[k]![j]!.isZero()) {
          appearsLater = true;
          break;
        }
      }
      if (!appearsLater) {
        maskFree |= 1 << j;
        toRemove.push(j);
      }
    }
    for (const j of toRemove) {
      varsToDo.delete(j);
    }

    // Multiply p and p1, applying the integration
    p = _prm_mul(ring, p, p1, maskFree, mn + 1);
  }

  // After processing all rows, p should have a single entry at key 0
  // containing the polynomial whose coefficients are the rook vector
  if (p.size === 0) {
    // All zero matrix
    const result = [ring.one()];
    for (let k = 1; k <= mn; k++) {
      result.push(ring.zero());
    }
    return result;
  }

  const coeffs = p.get(0);
  if (!coeffs || p.size !== 1) {
    // This shouldn't happen if the algorithm is correct
    throw new Error('Internal error in Butera-Pernici algorithm');
  }

  // Pad with zeros if needed
  const result: R[] = [];
  for (let k = 0; k <= mn; k++) {
    result.push(coeffs[k] ?? ring.zero());
  }

  return result;
}

/**
 * Multiply two polynomials in the nilpotent quotient ring R[eta_1,...,eta_n][t]
 * where eta_i^2 = 0, and integrate (set to 1) the variables in maskFree.
 *
 * @param ring - The coefficient ring
 * @param p1 - First polynomial (Map from bitmask to coefficient array)
 * @param p2 - Second polynomial
 * @param maskFree - Bitmask of variables to integrate
 * @param prec - Maximum degree to keep (exclusive)
 */
function _prm_mul<R extends RingElement>(
  ring: CoefficientRing<R>,
  p1: Map<number, R[]>,
  p2: Map<number, R[]>,
  maskFree: number,
  prec: number
): Map<number, R[]> {
  const result: Map<number, R[]> = new Map();

  for (const [exp1, v1] of p1) {
    for (const [exp2, v2] of p2) {
      // Skip if monomials share any variables (product is 0 due to nilpotency)
      if (exp1 & exp2) continue;

      // Multiply the polynomial coefficients
      const prod = _polyMul(ring, v1, v2, prec);
      if (_polyIsZero(prod)) continue;

      // Combine exponents and integrate free variables
      let exp = (exp1 | exp2) ^ ((exp1 | exp2) & maskFree);

      if (result.has(exp)) {
        result.set(exp, _polyAdd(ring, result.get(exp)!, prod));
      } else {
        result.set(exp, prod);
      }
    }
  }

  return result;
}

/**
 * Add two polynomials represented as coefficient arrays.
 */
function _polyAdd<R extends RingElement>(ring: CoefficientRing<R>, a: R[], b: R[]): R[] {
  const maxLen = Math.max(a.length, b.length);
  const result: R[] = [];
  for (let i = 0; i < maxLen; i++) {
    const ai = a[i] ?? ring.zero();
    const bi = b[i] ?? ring.zero();
    result.push(ai.add(bi) as R);
  }
  return result;
}

/**
 * Multiply two polynomials represented as coefficient arrays, truncating at prec.
 */
function _polyMul<R extends RingElement>(ring: CoefficientRing<R>, a: R[], b: R[], prec: number): R[] {
  const result: R[] = [];
  const maxDeg = Math.min(a.length + b.length - 1, prec);

  for (let k = 0; k < maxDeg; k++) {
    let coeff = ring.zero();
    for (let i = 0; i <= k && i < a.length; i++) {
      const j = k - i;
      if (j < b.length) {
        coeff = coeff.add(a[i]!.mul(b[j]!)) as R;
      }
    }
    result.push(coeff);
  }

  return result;
}

/**
 * Check if a polynomial is zero.
 */
function _polyIsZero<R extends RingElement>(p: R[]): boolean {
  return p.every((c) => c.isZero());
}

/**
 * Multiply a ring element by a scalar (integer).
 */
function _scalarMul<R extends RingElement>(ring: CoefficientRing<R>, scalar: number, elem: R): R {
  if (scalar === 0) return ring.zero();
  if (scalar === 1) return elem;
  if (scalar === -1) return elem.neg() as R;

  let result = ring.zero();
  const absScalar = Math.abs(scalar);
  for (let i = 0; i < absScalar; i++) {
    result = result.add(elem) as R;
  }
  return scalar > 0 ? (result as R) : (result.neg() as R);
}

/**
 * Compute rook vector using naive enumeration.
 * Only suitable for small matrices (up to ~12x12 for sparse matrices).
 * Complexity: O(C(p, k) * k) where p is number of nonzero positions.
 */
function _rook_vector_naive<R extends RingElement>(matrix: Matrix<R>): R[] {
  const m = matrix.nrows;
  const n = matrix.ncols;
  const ring = matrix.base_ring;
  const k = Math.min(m, n);

  // Find positions where the matrix has 1 (non-zero)
  const positions: Array<[number, number]> = [];
  for (let i = 0; i < m; i++) {
    for (let j = 0; j < n; j++) {
      if (!matrix.get(i, j).isZero()) {
        positions.push([i, j]);
      }
    }
  }

  // For very large position sets, this is infeasible
  // C(p, k) can be huge even for moderate p
  const maxPositions = 50; // Rough limit
  if (positions.length > maxPositions && k > 5) {
    throw new NotImplementedError(
      `rook_vector naive algorithm infeasible for ${positions.length} positions; use ButeraPernici or Ryser`
    );
  }

  // r_0 = 1 (empty placement)
  const result: R[] = [ring.one()];

  // Count placements for each number of rooks
  for (let numRooks = 1; numRooks <= k; numRooks++) {
    let count = 0;
    // Generate all combinations of numRooks positions
    for (const combo of _generate_combinations(positions, numRooks)) {
      // Check if this is a valid placement (no two rooks in same row/col)
      const rows = new Set<number>();
      const cols = new Set<number>();
      let valid = true;
      for (const [r, c] of combo) {
        if (rows.has(r) || cols.has(c)) {
          valid = false;
          break;
        }
        rows.add(r);
        cols.add(c);
      }
      if (valid) {
        count++;
      }
    }
    result.push(_scalarMul(ring, count, ring.one()));
  }

  return result;
}

/**
 * Generator for combinations to avoid storing all in memory.
 */
function* _generate_combinations<T>(arr: T[], k: number): Generator<T[]> {
  if (k === 0) {
    yield [];
    return;
  }
  if (arr.length < k) return;

  const indices = Array.from({ length: k }, (_, i) => i);

  while (true) {
    yield indices.map((i) => arr[i]!);

    // Find rightmost index that can be incremented
    let i = k - 1;
    while (i >= 0 && indices[i] === arr.length - k + i) {
      i--;
    }

    if (i < 0) break;

    indices[i]!++;
    for (let j = i + 1; j < k; j++) {
      indices[j] = indices[j - 1]! + 1;
    }
  }
}

/**
 * Generate all k-element subsets of {0, 1, ..., n-1}.
 */
function _choose_indices(n: number, k: number): number[][] {
  if (k === 0) return [[]];
  if (k > n) return [];

  const result: number[][] = [];

  function backtrack(start: number, current: number[]) {
    if (current.length === k) {
      result.push([...current]);
      return;
    }
    for (let i = start; i < n; i++) {
      current.push(i);
      backtrack(i + 1, current);
      current.pop();
    }
  }

  backtrack(0, []);
  return result;
}

/**
 * Return the product of row sums for selected columns.
 *
 * For each row, compute the sum of entries in the selected columns,
 * then return the product of all these row sums.
 *
 * @param matrix - The matrix
 * @param cols - Column indices
 * @returns The product of row sums
 * @see Reference: sage/matrix/matrix2.pyx:prod_of_row_sums
 */
export function prod_of_row_sums<R extends RingElement>(matrix: Matrix<R>, cols: number[]): R {
  const ring = matrix.base_ring;
  let prod = ring.one();

  for (let i = 0; i < matrix.nrows; i++) {
    let rowSum = ring.zero();
    for (const j of cols) {
      rowSum = rowSum.add(matrix.get(i, j)) as R;
    }
    prod = prod.mul(rowSum) as R;
  }

  return prod;
}

/**
 * Return the Hadamard bound on the determinant.
 *
 * The Hadamard bound is the product of the Euclidean norms of the rows
 * (or columns). For integer matrices, this gives an upper bound on |det(A)|.
 *
 * Note: This requires square root computation which is not available for
 * general rings. For exact computation over finite fields, this doesn't apply.
 *
 * @param matrix - The matrix
 * @returns The Hadamard bound (requires sqrt)
 * @see Reference: sage/matrix/matrix2.pyx:hadamard_bound
 */
export function hadamard_bound<R extends RingElement>(matrix: Matrix<R>): R {
  // Hadamard bound requires computing ||row_i|| = sqrt(sum of squares)
  // This requires sqrt which is not available for general rings
  throw new NotImplementedError('hadamard_bound requires rings with sqrt (e.g., RDF, RealField)');
}

/**
 * Subdivide the matrix.
 *
 * Divide the matrix into logical submatrices which can then be queried
 * and extracted. The subdivision is stored as metadata on the matrix.
 *
 * Note: This implementation stores subdivision data as a property on the matrix object.
 * The subdivision data consists of sorted lists of row and column positions.
 *
 * @param matrix - The matrix (mutated to add subdivision data)
 * @param row_lines - Row subdivision positions (single number or array)
 * @param col_lines - Column subdivision positions (single number or array)
 * @see Reference: sage/matrix/matrix2.pyx:subdivide
 */
export function subdivide<R extends RingElement>(
  matrix: Matrix<R>,
  row_lines?: number | number[] | [number[], number[]],
  col_lines?: number | number[]
): void {
  // Handle tuple argument: subdivide(matrix, (row_lines, col_lines))
  if (col_lines === undefined && row_lines !== undefined && Array.isArray(row_lines)) {
    if (row_lines.length === 2 && Array.isArray(row_lines[0]) && Array.isArray(row_lines[1])) {
      // It's a tuple (row_lines, col_lines)
      const tuple = row_lines as [number[], number[]];
      row_lines = tuple[0];
      col_lines = tuple[1];
    }
  }

  // Normalize row_lines
  let rowList: number[];
  if (row_lines === undefined || row_lines === null) {
    rowList = [];
  } else if (typeof row_lines === 'number') {
    rowList = [row_lines];
  } else {
    rowList = row_lines as number[];
  }

  // Normalize col_lines
  let colList: number[];
  if (col_lines === undefined || col_lines === null) {
    colList = [];
  } else if (typeof col_lines === 'number') {
    colList = [col_lines];
  } else {
    colList = col_lines as number[];
  }

  // If both empty, clear subdivisions
  const matrixWithSub = matrix as Matrix<R> & {
    _subdivisions?: [number[], number[]] | null;
  };

  if (rowList.length === 0 && colList.length === 0) {
    matrixWithSub._subdivisions = null;
    return;
  }

  // Sort and store with boundaries [0, ...positions..., nrows/ncols]
  const sortedRows = [0, ...rowList.map((x) => Math.floor(x)).sort((a, b) => a - b), matrix.nrows];
  const sortedCols = [0, ...colList.map((x) => Math.floor(x)).sort((a, b) => a - b), matrix.ncols];

  matrixWithSub._subdivisions = [sortedRows, sortedCols];
}

/**
 * Get a subdivision of the matrix.
 *
 * Return the (i,j)th submatrix of the matrix, according to a previously
 * set subdivision. Before a subdivision is set, the only valid arguments
 * are (0,0) which returns the entire matrix.
 *
 * @param matrix - The matrix
 * @param i - Row block index (0-based)
 * @param j - Column block index (0-based)
 * @returns The (i,j) subdivision block
 * @see Reference: sage/matrix/matrix2.pyx:subdivision
 */
export function subdivision<R extends RingElement>(
  matrix: Matrix<R>,
  i: number,
  j: number
): Matrix<R> {
  const matrixWithSub = matrix as Matrix<R> & {
    _subdivisions?: [number[], number[]] | null;
  };

  // Get or create default subdivisions
  let subs = matrixWithSub._subdivisions;
  if (subs === undefined || subs === null) {
    subs = [
      [0, matrix.nrows],
      [0, matrix.ncols],
    ];
  }

  const [rowSubs, colSubs] = subs;

  // Validate indices
  if (i < 0 || i >= rowSubs.length - 1) {
    throw new ValueError(`row subdivision index ${i} out of bounds`);
  }
  if (j < 0 || j >= colSubs.length - 1) {
    throw new ValueError(`column subdivision index ${j} out of bounds`);
  }

  const rowStart = rowSubs[i]!;
  const rowEnd = rowSubs[i + 1]!;
  const colStart = colSubs[j]!;
  const colEnd = colSubs[j + 1]!;

  const blockRows = rowEnd - rowStart;
  const blockCols = colEnd - colStart;

  const result = new Matrix<R>(matrix.base_ring, blockRows, blockCols);

  for (let bi = 0; bi < blockRows; bi++) {
    for (let bj = 0; bj < blockCols; bj++) {
      result.set(bi, bj, matrix.get(rowStart + bi, colStart + bj));
    }
  }

  return result;
}

/**
 * Return the subdivisions of the matrix.
 *
 * Returns the current subdivision of the matrix as a pair of lists:
 * (row_subdivisions, col_subdivisions). If no subdivision has been set,
 * returns ([], []).
 *
 * @param matrix - The matrix
 * @returns Pair (row_subdivisions, col_subdivisions) - the positions (not including 0 and nrows/ncols)
 * @see Reference: sage/matrix/matrix2.pyx:subdivisions
 */
export function subdivisions<R extends RingElement>(matrix: Matrix<R>): [number[], number[]] {
  const matrixWithSub = matrix as Matrix<R> & {
    _subdivisions?: [number[], number[]] | null;
  };

  if (matrixWithSub._subdivisions === undefined || matrixWithSub._subdivisions === null) {
    return [[], []];
  }

  const [rowSubs, colSubs] = matrixWithSub._subdivisions;

  // Return the internal positions (excluding 0 and nrows/ncols boundaries)
  return [rowSubs.slice(1, -1), colSubs.slice(1, -1)];
}

/**
 * Set a block of the matrix.
 *
 * Copies the entries of `block` into `matrix` starting at position (row, col).
 *
 * @param matrix - The matrix (mutated)
 * @param row - Starting row (0-based)
 * @param col - Starting column (0-based)
 * @param block - The block to insert
 * @see Reference: sage/matrix/matrix2.pyx:set_block
 */
export function set_block<R extends RingElement>(
  matrix: Matrix<R>,
  row: number,
  col: number,
  block: Matrix<R>
): void {
  if (row < 0 || col < 0) {
    throw new ValueError('row and col must be non-negative');
  }

  if (row + block.nrows > matrix.nrows || col + block.ncols > matrix.ncols) {
    throw new ValueError(
      `block of size ${block.nrows}x${block.ncols} starting at (${row}, ${col}) ` +
        `does not fit in matrix of size ${matrix.nrows}x${matrix.ncols}`
    );
  }

  for (let i = 0; i < block.nrows; i++) {
    for (let j = 0; j < block.ncols; j++) {
      matrix.set(row + i, col + j, block.get(i, j));
    }
  }
}

/**
 * Find entries satisfying a condition.
 *
 * @param matrix - The matrix
 * @param f - Predicate function
 * @param indices - Whether to return indices (default: false)
 * @returns Matching entries (if indices=false) or their (row, col) positions (if indices=true)
 * @see Reference: sage/matrix/matrix2.pyx:find
 */
export function find<R extends RingElement>(
  matrix: Matrix<R>,
  f: (x: R) => boolean,
  indices: boolean = false
): R[] | Array<[number, number]> {
  if (indices) {
    const result: Array<[number, number]> = [];
    for (let i = 0; i < matrix.nrows; i++) {
      for (let j = 0; j < matrix.ncols; j++) {
        if (f(matrix.get(i, j))) {
          result.push([i, j]);
        }
      }
    }
    return result;
  } else {
    const result: R[] = [];
    for (let i = 0; i < matrix.nrows; i++) {
      for (let j = 0; j < matrix.ncols; j++) {
        const entry = matrix.get(i, j);
        if (f(entry)) {
          result.push(entry);
        }
      }
    }
    return result;
  }
}

/**
 * Apply a function to each entry.
 *
 * Creates a new matrix where each entry is the result of applying phi
 * to the corresponding entry of the original matrix.
 *
 * @param matrix - The matrix
 * @param phi - The function to apply to each entry
 * @param ring - Ring for the result (if not provided, inferred from first output)
 * @param sparse - Whether to use sparse representation (not used in dense impl)
 * @returns A new matrix with phi applied elementwise
 * @see Reference: sage/matrix/matrix2.pyx:apply_map
 */
export function apply_map<R extends RingElement, S extends RingElement>(
  matrix: Matrix<R>,
  phi: (x: R) => S,
  ring?: CoefficientRing<S>,
  sparse?: boolean
): Matrix<S> {
  const m = matrix.nrows;
  const n = matrix.ncols;

  if (m === 0 || n === 0) {
    if (!ring) {
      throw new ValueError('must specify ring for empty matrix');
    }
    return new Matrix(ring, m, n);
  }

  // Apply phi to get the entries
  const entries: S[][] = [];
  for (let i = 0; i < m; i++) {
    entries.push([]);
    for (let j = 0; j < n; j++) {
      entries[i]!.push(phi(matrix.get(i, j)));
    }
  }

  // Determine ring from first entry if not provided
  // Use parent property (not method) or check for parent() method
  const firstResult = entries[0]![0]!;
  const firstResultWithParent = firstResult as unknown as {
    parent?: CoefficientRing<S> | (() => CoefficientRing<S>);
  };

  let actualRing: CoefficientRing<S>;
  if (ring) {
    actualRing = ring;
  } else if (typeof firstResultWithParent.parent === 'function') {
    actualRing = firstResultWithParent.parent();
  } else if (firstResultWithParent.parent) {
    actualRing = firstResultWithParent.parent as CoefficientRing<S>;
  } else {
    throw new ValueError('cannot determine result ring - please provide ring parameter');
  }

  return new Matrix(actualRing, m, n, entries);
}

/**
 * Apply a morphism to the entries.
 *
 * This is similar to apply_map but specifically for ring morphisms.
 *
 * @param matrix - The matrix
 * @param phi - The morphism (ring homomorphism)
 * @returns A new matrix with morphism applied
 * @see Reference: sage/matrix/matrix2.pyx:apply_morphism
 */
export function apply_morphism<R extends RingElement, S extends RingElement>(
  matrix: Matrix<R>,
  phi: (x: R) => S
): Matrix<S> {
  // apply_morphism is essentially the same as apply_map for our purposes
  return apply_map(matrix, phi);
}

/**
 * Substitute values into the matrix entries.
 *
 * All the arguments are transmitted unchanged to the method `subs` of
 * each coefficient. This function requires that the elements of the matrix
 * have a `subs` method (e.g., polynomial ring elements or symbolic expressions).
 *
 * @param matrix - The matrix
 * @param substitutions - A map of variable names to values, or positional arguments
 * @returns A new matrix with substitutions applied to each entry
 * @see Reference: sage/matrix/matrix2.pyx:subs
 */
export function subs<R extends RingElement>(
  matrix: Matrix<R>,
  substitutions: Record<string, unknown> | unknown
): Matrix<R> {
  const m = matrix.nrows;
  const n = matrix.ncols;
  const ring = matrix.base_ring;

  // Check if entries support subs
  if (m > 0 && n > 0) {
    const firstEntry = matrix.get(0, 0) as unknown as { subs?: (...args: unknown[]) => R };
    if (typeof firstEntry.subs !== 'function') {
      throw new TypeError('subs not defined for elements of the base ring');
    }
  }

  const entries: R[][] = [];
  for (let i = 0; i < m; i++) {
    entries.push([]);
    for (let j = 0; j < n; j++) {
      const entry = matrix.get(i, j) as unknown as { subs: (arg: unknown) => R };
      entries[i]!.push(entry.subs(substitutions));
    }
  }

  return new Matrix(ring, m, n, entries);
}

/**
 * Return the derivative of the matrix.
 *
 * Computes the derivative of each entry of the matrix with respect to
 * the given variable(s). This requires that the elements of the matrix
 * have a `derivative` method (e.g., polynomial ring elements).
 *
 * @param matrix - The matrix
 * @param variable - The variable to differentiate with respect to (or multiple variables)
 * @returns The derivative matrix
 * @see Reference: sage/matrix/matrix2.pyx:derivative
 */
export function derivative<R extends RingElement>(
  matrix: Matrix<R>,
  variable?: unknown
): Matrix<R> {
  const m = matrix.nrows;
  const n = matrix.ncols;
  const ring = matrix.base_ring;

  // Check if entries support derivative
  if (m > 0 && n > 0) {
    const firstEntry = matrix.get(0, 0) as unknown as {
      derivative?: (...args: unknown[]) => R;
      diff?: (...args: unknown[]) => R;
    };
    if (typeof firstEntry.derivative !== 'function' && typeof firstEntry.diff !== 'function') {
      throw new TypeError('derivative not defined for elements of the base ring');
    }
  }

  const entries: R[][] = [];
  for (let i = 0; i < m; i++) {
    entries.push([]);
    for (let j = 0; j < n; j++) {
      const entry = matrix.get(i, j) as unknown as {
        derivative?: (v?: unknown) => R;
        diff?: (v?: unknown) => R;
      };
      if (typeof entry.derivative === 'function') {
        entries[i]!.push(variable !== undefined ? entry.derivative(variable) : entry.derivative());
      } else if (typeof entry.diff === 'function') {
        entries[i]!.push(variable !== undefined ? entry.diff(variable) : entry.diff());
      } else {
        throw new TypeError('derivative not defined for elements of the base ring');
      }
    }
  }

  return new Matrix(ring, m, n, entries);
}

/**
 * Return a numerical approximation of the matrix.
 *
 * Returns a numerical approximation of the matrix with `prec` bits
 * (or decimal `digits`) of precision. Each entry is converted to a
 * numerical approximation.
 *
 * Note: This implementation requires that matrix entries have a
 * `numerical_approx` or `n` method. For matrices over finite fields
 * or integers without such methods, this will throw an error.
 *
 * @param matrix - The matrix
 * @param prec - Precision in bits
 * @param digits - Precision in digits (used if prec is not given)
 * @param algorithm - Algorithm to use (ignored for matrices)
 * @returns A numerical approximation
 * @see Reference: sage/matrix/matrix2.pyx:numerical_approx
 */
export function numerical_approx<R extends RingElement>(
  matrix: Matrix<R>,
  prec?: number,
  digits?: number,
  algorithm?: string
): Matrix<R> {
  const m = matrix.nrows;
  const n = matrix.ncols;
  const ring = matrix.base_ring;

  // Convert digits to prec if needed (roughly 3.32 bits per decimal digit)
  const actualPrec = prec ?? (digits !== undefined ? Math.ceil(digits * Math.log2(10)) : 53);

  // Check if entries support numerical_approx
  if (m > 0 && n > 0) {
    const firstEntry = matrix.get(0, 0) as unknown as {
      numerical_approx?: (prec?: number) => R;
      n?: (prec?: number) => R;
    };
    if (typeof firstEntry.numerical_approx !== 'function' && typeof firstEntry.n !== 'function') {
      throw new TypeError('numerical_approx not defined for elements of the base ring');
    }
  }

  const entries: R[][] = [];
  for (let i = 0; i < m; i++) {
    entries.push([]);
    for (let j = 0; j < n; j++) {
      const entry = matrix.get(i, j) as unknown as {
        numerical_approx?: (prec?: number) => R;
        n?: (prec?: number) => R;
      };
      if (typeof entry.numerical_approx === 'function') {
        entries[i]!.push(entry.numerical_approx(actualPrec));
      } else if (typeof entry.n === 'function') {
        entries[i]!.push(entry.n(actualPrec));
      } else {
        throw new TypeError('numerical_approx not defined for elements of the base ring');
      }
    }
  }

  return new Matrix(ring, m, n, entries);
}

/**
 * Return the denominator of the matrix.
 *
 * Returns the least common multiple of the denominators of all the
 * elements of the matrix. If there is no denominator function for
 * the base ring, or no LCM function for the denominators, raises a TypeError.
 *
 * For an empty matrix, returns 1.
 *
 * @param matrix - The matrix
 * @returns The LCM of denominators of all entries
 * @see Reference: sage/matrix/matrix2.pyx:denominator
 */
export function denominator<R extends RingElement>(matrix: Matrix<R>): R {
  const m = matrix.nrows;
  const n = matrix.ncols;
  const ring = matrix.base_ring;

  // For empty matrices, return 1
  if (m === 0 || n === 0) {
    return ring.one();
  }

  // Get the first element and check if it has a denominator method
  const firstEntry = matrix.get(0, 0) as unknown as {
    denominator?: () => R;
  };

  if (typeof firstEntry.denominator !== 'function') {
    throw new TypeError('denominator not defined for elements of the base ring');
  }

  // Get the denominator of the first element
  let result = firstEntry.denominator();

  // Compute LCM of all denominators
  for (let i = 0; i < m; i++) {
    for (let j = 0; j < n; j++) {
      if (i === 0 && j === 0) continue; // Skip first element

      const entry = matrix.get(i, j) as unknown as {
        denominator?: () => R;
      };

      if (typeof entry.denominator !== 'function') {
        throw new TypeError('denominator not defined for elements of the base ring');
      }

      const d = entry.denominator();

      // Compute LCM using result and d
      // LCM(a, b) = a * b / GCD(a, b)
      const resultWithLcm = result as unknown as {
        lcm?: (other: R) => R;
      };

      if (typeof resultWithLcm.lcm === 'function') {
        result = resultWithLcm.lcm(d);
      } else {
        throw new TypeError('lcm function not defined for denominators');
      }
    }
  }

  return result;
}

/**
 * Randomize the matrix entries.
 *
 * Mutates the matrix in place, filling entries with random elements
 * from the base ring.
 *
 * @param matrix - The matrix (mutated)
 * @param density - Density of non-zero entries (default: 1.0 = all entries)
 * @param nonzero - Whether all random entries should be non-zero (default: false)
 * @param args - Additional arguments for random element (not used)
 * @see Reference: sage/matrix/matrix2.pyx:randomize
 */
export function randomize<R extends RingElement>(
  matrix: Matrix<R>,
  density: number = 1.0,
  nonzero: boolean = false,
  ...args: unknown[]
): void {
  const ring = matrix.base_ring;

  // Check if ring has a random element method
  const ringWithRandom = ring as unknown as { random_element?: () => R };
  if (typeof ringWithRandom.random_element !== 'function') {
    throw new NotImplementedError('randomize requires a ring with random_element() method');
  }

  for (let i = 0; i < matrix.nrows; i++) {
    for (let j = 0; j < matrix.ncols; j++) {
      if (Math.random() < density) {
        let value = ringWithRandom.random_element();
        if (nonzero) {
          // Keep generating until we get a nonzero element
          let attempts = 0;
          while (value.isZero() && attempts < 100) {
            value = ringWithRandom.random_element();
            attempts++;
          }
        }
        matrix.set(i, j, value);
      }
    }
  }
}

/**
 * Return the bandwidth of the matrix.
 *
 * The lower bandwidth is the maximum row index minus column index for non-zero
 * entries below the diagonal. The upper bandwidth is the maximum column index
 * minus row index for non-zero entries above the diagonal.
 *
 * @param matrix - The matrix
 * @returns Pair (lower_bandwidth, upper_bandwidth)
 * @see Reference: sage/matrix/matrix2.pyx:get_bandwidth
 */
export function get_bandwidth<R extends RingElement>(matrix: Matrix<R>): [number, number] {
  let lower = 0;
  let upper = 0;

  for (let i = 0; i < matrix.nrows; i++) {
    for (let j = 0; j < matrix.ncols; j++) {
      if (!matrix.get(i, j).isZero()) {
        if (i > j) {
          // Below diagonal
          lower = Math.max(lower, i - j);
        } else if (j > i) {
          // Above diagonal
          upper = Math.max(upper, j - i);
        }
      }
    }
  }

  return [lower, upper];
}

/**
 * Return the matrix as a bipartite graph representation.
 *
 * Constructs a bipartite graph B representing the matrix uniquely.
 * Vertices are labeled 1 to nrows on the left (representing rows)
 * and nrows + 1 to nrows + ncols on the right (representing columns).
 * Each row is connected to each column with an edge weighted by
 * the value of the corresponding matrix entry.
 *
 * Note: This implementation returns a simplified representation as an object
 * with left vertices, right vertices, and edges. Full BipartiteGraph class
 * is not yet implemented.
 *
 * @param matrix - The matrix
 * @returns A bipartite graph representation with {left, right, edges}
 * @see Reference: sage/matrix/matrix2.pyx:as_bipartite_graph
 */
export function as_bipartite_graph<R extends RingElement>(
  matrix: Matrix<R>
): {
  left: number[];
  right: number[];
  edges: Array<[number, number, R]>;
} {
  const m = matrix.nrows;
  const n = matrix.ncols;

  // Left vertices are 1 to nrows (representing rows)
  const left: number[] = [];
  for (let i = 1; i <= m; i++) {
    left.push(i);
  }

  // Right vertices are nrows + 1 to nrows + ncols (representing columns)
  const right: number[] = [];
  for (let j = m + 1; j <= m + n; j++) {
    right.push(j);
  }

  // Edges connect row i to column j with weight matrix[i-1, j-m-1]
  const edges: Array<[number, number, R]> = [];
  for (let i = 0; i < m; i++) {
    for (let j = 0; j < n; j++) {
      edges.push([i + 1, m + j + 1, matrix.get(i, j)]);
    }
  }

  return { left, right, edges };
}

/**
 * Return the automorphisms of rows and columns.
 *
 * Returns the automorphisms of the matrix under permutations of rows and columns
 * as a list of pairs of permutations (row_perm, col_perm).
 *
 * An automorphism is a pair (sigma, tau) of permutations such that
 * applying sigma to the rows and tau to the columns leaves the matrix unchanged.
 *
 * Note: This is a simplified implementation that only finds the identity
 * automorphism and simple row/column swaps that preserve the matrix.
 * Full implementation requires graph isomorphism algorithms.
 *
 * @param matrix - The matrix
 * @returns Array of [row_permutation, col_permutation] pairs
 * @see Reference: sage/matrix/matrix2.pyx:automorphisms_of_rows_and_columns
 */
export function automorphisms_of_rows_and_columns<R extends RingElement>(
  matrix: Matrix<R>
): Array<[number[], number[]]> {
  const m = matrix.nrows;
  const n = matrix.ncols;

  // Result list of automorphisms
  const result: Array<[number[], number[]]> = [];

  // Identity permutation is always an automorphism
  const identityRow: number[] = [];
  const identityCol: number[] = [];
  for (let i = 0; i < m; i++) identityRow.push(i);
  for (let j = 0; j < n; j++) identityCol.push(j);
  result.push([identityRow, identityCol]);

  // For small matrices, we can check all permutations
  // For larger matrices, this is expensive
  if (m > 6 || n > 6) {
    // Only return identity for large matrices
    // Full implementation would use graph automorphism algorithms
    return result;
  }

  // Generate all permutations of rows
  const rowPerms = _generatePermutations(m);
  const colPerms = _generatePermutations(n);

  // Check each combination
  for (const rowPerm of rowPerms) {
    for (const colPerm of colPerms) {
      // Skip identity (already added)
      if (_isIdentityPerm(rowPerm) && _isIdentityPerm(colPerm)) {
        continue;
      }

      // Check if this permutation preserves the matrix
      if (_checkAutomorphism(matrix, rowPerm, colPerm)) {
        result.push([rowPerm, colPerm]);
      }
    }
  }

  return result;
}

/**
 * Generate all permutations of [0, 1, ..., n-1]
 */
function _generatePermutations(n: number): number[][] {
  if (n === 0) return [[]];
  if (n === 1) return [[0]];

  const result: number[][] = [];

  function permute(arr: number[], start: number) {
    if (start === arr.length - 1) {
      result.push([...arr]);
      return;
    }

    for (let i = start; i < arr.length; i++) {
      // Swap
      [arr[start], arr[i]] = [arr[i]!, arr[start]!];
      permute(arr, start + 1);
      // Swap back
      [arr[start], arr[i]] = [arr[i]!, arr[start]!];
    }
  }

  const initial: number[] = [];
  for (let i = 0; i < n; i++) initial.push(i);
  permute(initial, 0);

  return result;
}

/**
 * Check if a permutation is the identity
 */
function _isIdentityPerm(perm: number[]): boolean {
  for (let i = 0; i < perm.length; i++) {
    if (perm[i] !== i) return false;
  }
  return true;
}

/**
 * Check if applying rowPerm to rows and colPerm to columns preserves the matrix
 */
function _checkAutomorphism<R extends RingElement>(
  matrix: Matrix<R>,
  rowPerm: number[],
  colPerm: number[]
): boolean {
  const m = matrix.nrows;
  const n = matrix.ncols;

  for (let i = 0; i < m; i++) {
    for (let j = 0; j < n; j++) {
      // Original entry at (i, j)
      const original = matrix.get(i, j);
      // Entry at permuted position
      const permuted = matrix.get(rowPerm[i]!, colPerm[j]!);

      if (!original.eq(permuted)) {
        return false;
      }
    }
  }

  return true;
}

/**
 * Return the permutation normal form of the matrix.
 *
 * Takes the set of matrices that are the given matrix permuted by any row
 * and column permutation, and returns the maximal one of the set where
 * matrices are ordered lexicographically going along each row.
 *
 * @param matrix - The matrix
 * @param check - If true, return a tuple (maximal_matrix, (row_perm, col_perm))
 * @returns The permutation normal form (or tuple if check=true)
 * @see Reference: sage/matrix/matrix2.pyx:permutation_normal_form
 */
export function permutation_normal_form<R extends RingElement>(
  matrix: Matrix<R>,
  check?: boolean
): Matrix<R> | [Matrix<R>, [number[], number[]]] {
  const m = matrix.nrows;
  const n = matrix.ncols;
  const ring = matrix.base_ring;

  if (m === 0 || n === 0) {
    if (check) {
      return [matrix.copy(), [[], []]];
    }
    return matrix.copy();
  }

  // For small matrices, enumerate all permutations
  // For large matrices, use a heuristic approach
  if (m > 6 || n > 6) {
    // Heuristic: sort rows and columns by their "value"
    // This won't find the true maximum but gives a canonical form
    return _heuristicNormalForm(matrix, check);
  }

  // Generate all permutations
  const rowPerms = _generatePermutations(m);
  const colPerms = _generatePermutations(n);

  let maxMatrix: Matrix<R> | null = null;
  let maxRowPerm: number[] = [];
  let maxColPerm: number[] = [];

  for (const rowPerm of rowPerms) {
    for (const colPerm of colPerms) {
      const permuted = _applyPermutation(matrix, rowPerm, colPerm);

      if (maxMatrix === null || _compareMatricesLex(permuted, maxMatrix) > 0) {
        maxMatrix = permuted;
        maxRowPerm = rowPerm;
        maxColPerm = colPerm;
      }
    }
  }

  if (check) {
    return [maxMatrix!, [maxRowPerm, maxColPerm]];
  }

  return maxMatrix!;
}

/**
 * Apply row and column permutation to a matrix
 */
function _applyPermutation<R extends RingElement>(
  matrix: Matrix<R>,
  rowPerm: number[],
  colPerm: number[]
): Matrix<R> {
  const m = matrix.nrows;
  const n = matrix.ncols;
  const ring = matrix.base_ring;

  const entries: R[][] = [];
  for (let i = 0; i < m; i++) {
    entries.push([]);
    for (let j = 0; j < n; j++) {
      entries[i]!.push(matrix.get(rowPerm[i]!, colPerm[j]!));
    }
  }

  return new Matrix(ring, m, n, entries);
}

/**
 * Compare two matrices lexicographically
 * Returns positive if A > B, negative if A < B, 0 if equal
 */
function _compareMatricesLex<R extends RingElement>(A: Matrix<R>, B: Matrix<R>): number {
  const m = A.nrows;
  const n = A.ncols;

  for (let i = 0; i < m; i++) {
    for (let j = 0; j < n; j++) {
      const a = A.get(i, j);
      const b = B.get(i, j);

      if (a.eq(b)) continue;

      // Compare using toString as a fallback
      // For proper comparison, elements should have a compare method
      const aWithCmp = a as unknown as {
        cmp?: (other: R) => number;
        __lt__?: (other: R) => boolean;
      };
      const bAsR = b as R;

      if (typeof aWithCmp.cmp === 'function') {
        return aWithCmp.cmp(bAsR);
      }

      if (typeof aWithCmp.__lt__ === 'function') {
        return aWithCmp.__lt__(bAsR) ? -1 : 1;
      }

      // Fallback: string comparison
      const aStr = a.toString();
      const bStr = b.toString();
      if (aStr < bStr) return -1;
      if (aStr > bStr) return 1;
    }
  }

  return 0;
}

/**
 * Heuristic normal form for large matrices
 */
function _heuristicNormalForm<R extends RingElement>(
  matrix: Matrix<R>,
  check?: boolean
): Matrix<R> | [Matrix<R>, [number[], number[]]] {
  // Simple heuristic: sort rows by their string representation (descending)
  // and columns similarly
  const m = matrix.nrows;
  const n = matrix.ncols;

  // Get row indices sorted by row values (descending)
  const rowIndices: number[] = [];
  for (let i = 0; i < m; i++) rowIndices.push(i);

  rowIndices.sort((i1, i2) => {
    for (let j = 0; j < n; j++) {
      const a = matrix.get(i1, j);
      const b = matrix.get(i2, j);
      if (!a.eq(b)) {
        const aStr = a.toString();
        const bStr = b.toString();
        // Descending order (we want max)
        if (aStr > bStr) return -1;
        if (aStr < bStr) return 1;
      }
    }
    return 0;
  });

  // Apply row permutation first
  const rowPermuted = _applyPermutation(
    matrix,
    rowIndices,
    Array.from({ length: n }, (_, i) => i)
  );

  // Get column indices sorted by column values (descending)
  const colIndices: number[] = [];
  for (let j = 0; j < n; j++) colIndices.push(j);

  colIndices.sort((j1, j2) => {
    for (let i = 0; i < m; i++) {
      const a = rowPermuted.get(i, j1);
      const b = rowPermuted.get(i, j2);
      if (!a.eq(b)) {
        const aStr = a.toString();
        const bStr = b.toString();
        if (aStr > bStr) return -1;
        if (aStr < bStr) return 1;
      }
    }
    return 0;
  });

  const result = _applyPermutation(
    rowPermuted,
    Array.from({ length: m }, (_, i) => i),
    colIndices
  );

  if (check) {
    return [result, [rowIndices, colIndices]];
  }

  return result;
}

/**
 * Check if a matrix is a permutation of another.
 *
 * Returns True if there exists a permutation of rows and columns
 * sending A to B and False otherwise.
 *
 * @param A - First matrix
 * @param B - Second matrix
 * @param check - If true, return [boolean, permutation] where permutation is
 *                (row_perm, col_perm) or null if not a permutation
 * @returns True if A is a permutation of B (or tuple if check=true)
 * @see Reference: sage/matrix/matrix2.pyx:is_permutation_of
 */
export function is_permutation_of<R extends RingElement>(
  A: Matrix<R>,
  B: Matrix<R>,
  check?: boolean
): boolean | [boolean, [number[], number[]] | null] {
  // Different dimensions means definitely not a permutation
  if (A.nrows !== B.nrows || A.ncols !== B.ncols) {
    if (check) {
      return [false, null];
    }
    return false;
  }

  const m = A.nrows;
  const n = A.ncols;

  // Empty matrices are permutations of each other
  if (m === 0 || n === 0) {
    if (check) {
      return [true, [[], []]];
    }
    return true;
  }

  // Quick check: multiset of entries must be the same
  const entriesA: string[] = [];
  const entriesB: string[] = [];
  for (let i = 0; i < m; i++) {
    for (let j = 0; j < n; j++) {
      entriesA.push(A.get(i, j).toString());
      entriesB.push(B.get(i, j).toString());
    }
  }
  entriesA.sort();
  entriesB.sort();

  for (let i = 0; i < entriesA.length; i++) {
    if (entriesA[i] !== entriesB[i]) {
      if (check) {
        return [false, null];
      }
      return false;
    }
  }

  // For small matrices, try all permutations
  if (m <= 6 && n <= 6) {
    const rowPerms = _generatePermutations(m);
    const colPerms = _generatePermutations(n);

    for (const rowPerm of rowPerms) {
      for (const colPerm of colPerms) {
        const permuted = _applyPermutation(A, rowPerm, colPerm);
        if (permuted.eq(B)) {
          if (check) {
            return [true, [rowPerm, colPerm]];
          }
          return true;
        }
      }
    }

    if (check) {
      return [false, null];
    }
    return false;
  }

  // For larger matrices, compare permutation normal forms
  const normalA = permutation_normal_form(A, false) as Matrix<R>;
  const normalB = permutation_normal_form(B, false) as Matrix<R>;

  const isPermutation = normalA.eq(normalB);

  if (check) {
    if (isPermutation) {
      // We found they're permutations but don't have the exact permutation
      // This would require more work to find
      return [true, null];
    }
    return [false, null];
  }

  return isPermutation;
}
