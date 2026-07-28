/**
 * @module sage/matrix/matrix_special
 * @description Special matrix constructors and algorithms
 *
 * Port of: sage/matrix/special.py, sage/matrix/berlekamp_massey.py
 */

import { NotImplementedError, ValueError } from '../errors.js';
import { current_randstate } from '../misc/randstate.js';
import type { CoefficientRing, RingElement } from '../rings/polynomial/polynomial_element.js';
import { Polynomial } from '../rings/polynomial/polynomial_element.js';
import { PolynomialRing } from '../rings/polynomial/polynomial_ring.js';
import { QQ } from '../rings/rational_field.js';
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
        if (density >= 1.0 || _randomFraction() < density) {
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

  if (actualNumPivots === 0) {
    return result;
  }

  // Mirror SageMath: column 0 is always a pivot, the remaining pivot columns
  // are a random subset of the columns 1, ..., ncols-1.
  const subset: number[] = [];
  for (let j = 1; j < ncols; j++) {
    subset.push(j);
  }
  _shuffle(subset);
  const pivotCols = [0, ...subset.slice(0, actualNumPivots - 1).sort((a, b) => a - b)];

  // Build the RREF matrix: leading ones at the pivot positions ...
  for (let pivotRow = 0; pivotRow < actualNumPivots; pivotRow++) {
    result.set(pivotRow, pivotCols[pivotRow]!, ring.one());
  }

  // ... and random entries in the non-pivot columns, above the pivot rows
  for (let pivotIndex = 0; pivotIndex < actualNumPivots - 1; pivotIndex++) {
    for (let j = pivotCols[pivotIndex]! + 1; j < pivotCols[pivotIndex + 1]!; j++) {
      for (let i = 0; i <= pivotIndex; i++) {
        result.set(i, j, ringWithRandom.random_element());
      }
    }
  }
  for (let j = pivotCols[actualNumPivots - 1]! + 1; j < ncols; j++) {
    for (let i = 0; i < actualNumPivots; i++) {
      result.set(i, j, ringWithRandom.random_element());
    }
  }

  return result;
}

/**
 * A uniform fraction in [0, 1) drawn from SageMath's global random state.
 *
 * Used for the `density` parameter, so that seeding with `set_random_seed`
 * makes the random matrix constructors reproducible (JS `Math.random()` cannot
 * be seeded).
 */
function _randomFraction(): number {
  return Number(current_randstate().random_bits(53)) / 2 ** 53;
}

/**
 * Fisher-Yates shuffle driven by SageMath's global random state.
 */
function _shuffle<T>(arr: T[]): void {
  const rs = current_randstate();
  for (let i = arr.length - 1; i > 0; i--) {
    const j = Number(rs.randint(0n, BigInt(i)));
    [arr[i], arr[j]] = [arr[j]!, arr[i]!];
  }
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

  if (upper_bound !== undefined) {
    // SageMath only supports size control over ZZ and QQ, where it repeatedly
    // rejects row operations that push an entry past the bound.  This port
    // works over arbitrary rings with random_element() and has no notion of
    // absolute value, so the option is refused rather than silently ignored.
    throw new NotImplementedError(
      'size control (upper_bound) for random_echelonizable_matrix is only ' +
        'implemented over ZZ and QQ in SageMath; not supported by this port'
    );
  }

  // Start with a random RREF matrix of the given rank.
  const result = random_rref_matrix(ring, nrows, actualNcols, actualRank);

  // Scramble it using *only* transvections (adding a multiple of one row to
  // another).  These preserve the determinant, so for a full-rank square
  // matrix the result has determinant one, which is what makes
  // random_unimodular_matrix work.
  const addMultipleOfRow = (i: number, j: number, s: R): void => {
    for (let k = 0; k < actualNcols; k++) {
      result.set(i, k, result.get(i, k).add(s.mul(result.get(j, k)) as R) as R);
    }
  };

  for (let pivots = actualRank - 1; pivots >= 0; pivots--) {
    let row_index = 0;
    while (row_index < nrows) {
      if (pivots === row_index) {
        row_index += 1;
      }
      if (pivots !== row_index && row_index !== nrows) {
        // The pivot of index `pivots` lives in row `pivots` of the RREF
        addMultipleOfRow(row_index, pivots, ringWithRandom.random_element());
        row_index += 1;
      }
    }
  }
  if (nrows > 1) {
    const j = Number(current_randstate().randint(1n, BigInt(nrows - 1)));
    addMultipleOfRow(0, j, ringWithRandom.random_element());
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
 * Generate a random unimodular (determinant **one**) matrix of a desired size
 * over a desired ring.
 *
 * As in SageMath this simply delegates to
 * {@link random_echelonizable_matrix} with full rank: that routine only ever
 * adds multiples of one row to another, so the determinant stays equal to the
 * determinant of the full-rank RREF matrix it starts from, namely one.
 *
 * @param ring - The base ring
 * @param n - Size of the matrix
 * @param upper_bound - Upper bound on entry size (SageMath: only over ZZ or QQ)
 * @param max_tries - Number of tries used to generate each new random row
 * @returns A random matrix of determinant one
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

  // random_echelonizable_matrix() always returns a determinant one matrix if
  // given full rank -- it only ever applies transvections to the identity.
  return random_echelonizable_matrix(ring, n, n, n, upper_bound, max_tries);
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
 * Create a random matrix that diagonalizes nicely.
 *
 * The eigenvalues are elements of the base ring; when they are not supplied
 * they are drawn as `ring(randint(-10, 10))` and grouped into eigenspaces by
 * multiplicity, exactly as in SageMath.
 *
 * The matrix is `E * D * E^{-1}` where `D` is the diagonal matrix of
 * eigenvalues and `E` is built from the identity by column and row
 * transvections, hence unimodular.  `E^{-1}` is accumulated alongside `E` from
 * the inverse operations, so no division is ever performed and the routine
 * works over any (commutative) base ring.
 *
 * @param ring - The base ring
 * @param n - Size of the matrix
 * @param eigenvalues - The list of desired eigenvalues (requires `dimensions`)
 * @param dimensions - The list of dimensions of the corresponding eigenspaces
 * @returns A square, diagonalizable matrix
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

  if (eigenvalues !== undefined && dimensions === undefined) {
    throw new ValueError(
      'the list of eigenvalues must have a list of dimensions corresponding to each eigenvalue.'
    );
  }
  if (eigenvalues === undefined && dimensions !== undefined) {
    throw new ValueError('the list of dimensions must have a list of corresponding eigenvalues.');
  }

  if (n === 0) {
    return new Matrix(ring, 0, 0);
  }

  let values: R[];
  let dims: number[];

  if (eigenvalues === undefined) {
    // Create a list with `n` random eigenvalues in [-10, 10], then collapse it
    // to the distinct values together with their multiplicities.
    const rs = current_randstate();
    const drawn: bigint[] = [];
    for (let i = 0; i < n; i++) {
      drawn.push(rs.randint(-10n, 10n));
    }
    drawn.sort((a, b) => (a < b ? -1 : a > b ? 1 : 0));
    const distinct: bigint[] = [];
    for (const v of drawn) {
      if (!distinct.includes(v)) {
        distinct.push(v);
      }
    }
    values = distinct.map((v) => ring.__call__(v));
    dims = distinct.map((v) => drawn.filter((x) => x === v).length);
  } else {
    values = eigenvalues.slice();
    dims = dimensions!.slice();
  }

  const size_check = dims.reduce((a, b) => a + b, 0);
  if (n !== size_check) {
    throw new ValueError('the size of the matrix must equal the sum of the dimensions.');
  }
  if (Math.min(...dims) < 1) {
    throw new ValueError('eigenspaces must have a dimension of at least 1.');
  }
  if (values.length !== dims.length) {
    throw new ValueError(
      'each eigenvalue must have a corresponding dimension and each dimension a corresponding eigenvalue.'
    );
  }

  // Sort the dimensions in order of increasing size, keeping the eigenvalues
  // in step (a stable sort on the dimension only, as in SageMath's
  // sorted(zip(dimensions, eigenvalues)) up to ties in the dimension).
  const order = dims.map((_, i) => i);
  order.sort((a, b) => dims[a]! - dims[b]!);
  dims = order.map((i) => dims[i]!);
  values = order.map((i) => values[i]!);

  // Create the matrix of eigenvalues on the diagonal.
  const diagonal_matrix = new Matrix<R>(ring, n, n);
  let up_bound = 0;
  let low_bound = 0;
  for (let row_index = 0; row_index < dims.length; row_index++) {
    up_bound += dims[row_index]!;
    for (let entry = low_bound; entry < up_bound; entry++) {
      diagonal_matrix.set(entry, entry, values[row_index]!);
    }
    low_bound += dims[row_index]!;
  }

  // Create a matrix to hold each of the eigenvectors as its columns, beginning
  // with the identity matrix so that after the row and column operations the
  // resulting matrix is unimodular.  `inv` tracks its inverse.
  const eigenvector_matrix = new Matrix<R>(ring, n, n);
  const inv = new Matrix<R>(ring, n, n);
  for (let i = 0; i < n; i++) {
    eigenvector_matrix.set(i, i, ring.one());
    inv.set(i, i, ring.one());
  }

  // E <- E * T with T = I + s * e_{c2} e_{c1}^T; correspondingly
  // E^{-1} <- T^{-1} * E^{-1}, i.e. row c2 of the inverse loses s * row c1.
  const add_multiple_of_column = (c1: number, c2: number, s: bigint): void => {
    if (s === 0n) return;
    const scalar = ring.__call__(s);
    for (let i = 0; i < n; i++) {
      eigenvector_matrix.set(
        i,
        c1,
        eigenvector_matrix.get(i, c1).add(eigenvector_matrix.get(i, c2).mul(scalar) as R) as R
      );
    }
    for (let j = 0; j < n; j++) {
      inv.set(c2, j, inv.get(c2, j).sub(inv.get(c1, j).mul(scalar) as R) as R);
    }
  };

  // E <- T * E with T = I + s * e_{r1} e_{r2}^T; correspondingly
  // E^{-1} <- E^{-1} * T^{-1}, i.e. column r2 of the inverse loses s * column r1.
  const add_multiple_of_row = (r1: number, r2: number, s: bigint): void => {
    if (s === 0n) return;
    const scalar = ring.__call__(s);
    for (let j = 0; j < n; j++) {
      eigenvector_matrix.set(
        r1,
        j,
        eigenvector_matrix.get(r1, j).add(eigenvector_matrix.get(r2, j).mul(scalar) as R) as R
      );
    }
    for (let i = 0; i < n; i++) {
      inv.set(i, r2, inv.get(i, r2).sub(inv.get(i, r1).mul(scalar) as R) as R);
    }
  };

  const rs = current_randstate();
  const max_dim = Math.max(...dims);
  const min_dim = Math.min(...dims);

  // Assign the "protected" ones: a one is placed `dimensions[k]` rows up from
  // the bottom row and then diagonally down to the right.  Because the
  // dimensions are sorted increasingly, the target row is always strictly
  // below the target column and the source row is still a unit vector, so
  // setting the entry to one is exactly the transvection
  // "row r += 1 * row c"; performing it as such keeps `inv` in step.
  let upper_limit = 0;
  let lower_limit = 0;
  for (let dimension_index = 0; dimension_index < dims.length - 1; dimension_index++) {
    upper_limit += dims[dimension_index]!;
    let lowest_index_row_with_one = n - dims[dimension_index]!;
    for (let eigen_ones = lower_limit; eigen_ones < upper_limit; eigen_ones++) {
      if (lowest_index_row_with_one <= eigen_ones) {
        throw new NotImplementedError(
          'random_diagonalizable_matrix: unexpected eigenvector layout'
        );
      }
      add_multiple_of_row(lowest_index_row_with_one, eigen_ones, 1n);
      lowest_index_row_with_one += 1;
    }
    lower_limit += dims[dimension_index]!;
  }

  // A list giving the eigenvalue dimension corresponding to each column.
  const dimension_check: number[] = [];
  for (let i = 0; i < dims.length; i++) {
    for (let k = 0; k < dims[i]!; k++) {
      dimension_check.push(dims[i]!);
    }
  }

  // Fill the entries in the rows spanned by the protected ones using column
  // multiples, then the remaining rows using scalar row addition.
  for (
    let dimension_multiplicity = max_dim;
    dimension_multiplicity > min_dim;
    dimension_multiplicity--
  ) {
    const highest_one_row = n - dimension_multiplicity;
    // Find the column with the protected one in the lowest indexed row
    let highest_one_column = 0;
    while (eigenvector_matrix.get(highest_one_row, highest_one_column).isZero()) {
      highest_one_column += 1;
    }
    // dimension_check determines if a column has a low enough eigenvalue
    // dimension to take a column multiple
    for (
      let bottom_entry_filler = 0;
      bottom_entry_filler < dimension_check.length;
      bottom_entry_filler++
    ) {
      if (
        dimension_check[bottom_entry_filler]! < dimension_multiplicity &&
        eigenvector_matrix.get(highest_one_row, bottom_entry_filler).isZero()
      ) {
        // range of the multiplier determined experimentally by SageMath
        add_multiple_of_column(bottom_entry_filler, highest_one_column, rs.randint(-4n, 4n));
      }
    }
  }

  for (let row = n - max_dim; row < n; row++) {
    for (let upper_row = 0; upper_row < n - max_dim; upper_row++) {
      add_multiple_of_row(upper_row, row, rs.randint(-4n, 4n));
    }
  }

  return eigenvector_matrix.mul(diagonal_matrix).mul(inv);
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
 * Create a square matrix that corresponds to a row operation or a column
 * operation.
 *
 * The row-operation forms are (with `E` the returned matrix, so that `E*A`
 * performs the operation on `A`):
 *
 * - `{row1: i, row2: j}` — the matrix which swaps rows `i` and `j`
 * - `{row1: i, scale: s}` — the matrix which multiplies row `i` by `s`
 * - `{row1: i, row2: j, scale: s}` — the matrix which multiplies row `j` by
 *   `s` and adds it to row `i`
 *
 * Column operations are obtained in the analogous way by replacing `row1` by
 * `col1` and `row2` by `col2`; the resulting matrix is the transpose of the
 * corresponding row-operation matrix, and `A*E` performs the column operation.
 *
 * @param ring - The base ring
 * @param n - Size of the matrix (must be 1 or greater)
 * @param options - Operation specification (`row1`/`row2` or `col1`/`col2`, and `scale`)
 * @returns An elementary matrix
 * @see Reference: sage/matrix/special.py:elementary_matrix
 */
export function elementary_matrix<R extends RingElement>(
  ring: CoefficientRing<R>,
  n: number,
  options: { row1?: number; row2?: number; col1?: number; col2?: number; scale?: R }
): Matrix<R> {
  if (n <= 0) {
    throw new ValueError(`size of elementary matrix must be 1 or greater, not ${n}`);
  }

  const { scale } = options;

  // row operations or column operations?
  // a column operation matrix is the transpose of a row operation matrix
  if (options.row1 === undefined && options.col1 === undefined) {
    throw new ValueError('row1 or col1 must be specified');
  }
  if (options.row1 !== undefined && options.col1 !== undefined) {
    throw new ValueError('cannot specify both row1 and col1');
  }

  const rowop = options.row1 !== undefined;
  const opstring = rowop ? 'row' : 'column';
  const row1 = rowop ? options.row1! : options.col1!;
  const row2 = rowop ? options.row2 : options.col2;

  if (row1 < 0 || row1 >= n) {
    throw new ValueError(
      `${opstring} of elementary matrix must be positive and smaller than ${n}, not ${row1}`
    );
  }
  if (row2 !== undefined && (row2 < 0 || row2 >= n)) {
    throw new ValueError(
      `${opstring} of elementary matrix must be positive and smaller than ${n}, not ${row2}`
    );
  }

  // Start with the identity matrix
  const entries: R[][] = [];
  for (let i = 0; i < n; i++) {
    entries.push([]);
    for (let j = 0; j < n; j++) {
      entries[i]!.push(i === j ? ring.one() : ring.zero());
    }
  }

  if (row2 === undefined && scale === undefined) {
    throw new ValueError('insufficient parameters provided to construct elementary matrix');
  } else if (row2 !== undefined && scale !== undefined) {
    // Add a multiple of one row to another
    if (row1 === row2) {
      throw new ValueError(`cannot add a multiple of a ${opstring} to itself`);
    }
    entries[row1]![row2] = scale;
  } else if (row2 !== undefined) {
    // Swap two rows
    entries[row1]![row1] = ring.zero();
    entries[row2]![row2] = ring.zero();
    entries[row1]![row2] = ring.one();
    entries[row2]![row1] = ring.one();
  } else {
    // Scale a row
    if (scale!.isZero()) {
      throw new ValueError(`scale parameter of ${opstring} of elementary matrix must be nonzero`);
    }
    entries[row1]![row1] = scale!;
  }

  const elem = new Matrix(ring, n, n, entries);
  return rowop ? elem : elem.transpose();
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
 * Determine the dimensions of the rows and columns when assembling the
 * matrices in `sub_matrices` in a rectangular grid.
 *
 * Non-zero scalars are considered to be square matrices of any size, and
 * zeroes are considered to be zero matrices of any size.  A ValueError is
 * raised if there is insufficient or conflicting information.
 *
 * @returns A pair `[row_heights, col_widths]`
 * @see Reference: sage/matrix/special.py:_determine_block_matrix_grid
 */
function _determine_block_matrix_grid<R extends RingElement>(
  sub_matrices: Array<Array<Matrix<R> | R | number | 0>>
): [number[], number[]] {
  const nrows = sub_matrices.length;
  if (nrows === 0) {
    return [[], []];
  }
  const ncols = sub_matrices[0]!.length;
  if (ncols === 0) {
    return [new Array<number>(nrows).fill(0), []];
  }

  const row_heights: Array<number | null> = new Array<number | null>(nrows).fill(null);
  const col_widths: Array<number | null> = new Array<number | null>(ncols).fill(null);

  let changing = true;
  while (changing) {
    changing = false;
    for (let i = 0; i < nrows; i++) {
      for (let j = 0; j < ncols; j++) {
        const M = sub_matrices[i]![j];
        let sub_width: number | null = null;
        let sub_height: number | null = null;
        if (M instanceof Matrix) {
          sub_width = M.ncols;
          sub_height = M.nrows;
        } else if (!_isZeroScalar(M)) {
          // nonzero scalar is interpreted as a square matrix
          sub_width = row_heights[i] === null ? col_widths[j]! : row_heights[i]!;
          sub_height = sub_width;
        }
        if (sub_width !== null) {
          if (col_widths[j] === null) {
            changing = true;
            col_widths[j] = sub_width;
          } else if (col_widths[j] !== sub_width) {
            throw new ValueError('incompatible submatrix widths');
          }
        }
        if (sub_height !== null) {
          if (row_heights[i] === null) {
            changing = true;
            row_heights[i] = sub_height;
          } else if (row_heights[i] !== sub_height) {
            throw new ValueError('incompatible submatrix heights');
          }
        }
      }
    }
  }

  if (row_heights.includes(null) || col_widths.includes(null)) {
    throw new ValueError('insufficient information to determine dimensions.');
  }

  return [row_heights as number[], col_widths as number[]];
}

/**
 * Test whether the matrices in `sub_matrices` fit in a rectangular matrix when
 * assembled a row at a time.
 *
 * @returns `[row_heights, zero_widths, total_width]`
 * @see Reference: sage/matrix/special.py:_determine_block_matrix_rows
 */
function _determine_block_matrix_rows<R extends RingElement>(
  sub_matrices: Array<Array<Matrix<R> | R | number | 0>>
): [number[], number[], number] {
  let total_width: number | null = null;
  const row_heights: Array<number | null> = new Array<number | null>(sub_matrices.length).fill(
    null
  );
  const zero_widths: number[] = new Array<number>(sub_matrices.length).fill(0);

  // We first do a pass to see if we can determine the width
  let unknowns = false;
  for (let i = 0; i < sub_matrices.length; i++) {
    const R_ = sub_matrices[i]!;
    let height: number | null = null;
    let found_zeroes = false;
    for (const M of R_) {
      if (M instanceof Matrix) {
        if (height === null) {
          height = M.nrows;
        } else if (height !== M.nrows) {
          throw new ValueError('incompatible submatrix heights');
        }
      } else if (_isZeroScalar(M)) {
        found_zeroes = true;
      }
    }
    if (R_.length === 0) {
      height = 0;
    }

    if (height !== null && !found_zeroes) {
      let width = 0;
      for (const M of R_) {
        width += M instanceof Matrix ? M.ncols : height;
      }
      if (total_width === null) {
        total_width = width;
      } else if (total_width !== width) {
        throw new ValueError('incompatible submatrix widths');
      }
      row_heights[i] = height;
    } else {
      unknowns = true;
    }
  }

  if (total_width === null) {
    throw new ValueError('insufficient information to determine submatrix widths');
  }

  if (unknowns) {
    for (let i = 0; i < sub_matrices.length; i++) {
      if (row_heights[i] !== null) continue;
      const R_ = sub_matrices[i]!;
      // 0: no zeroes found, 1: consecutive zeroes found,
      // 2: consecutive zeroes followed by nonzero found, 3: non-consecutive zeroes
      let zero_state = 0;
      let scalars = 0;
      let width = 0;
      let height: number | null = null;
      for (const M of R_) {
        if (M instanceof Matrix) {
          height = M.nrows;
          width += M.ncols;
          if (zero_state === 1) zero_state = 2;
        } else if (_isZeroScalar(M)) {
          if (zero_state === 0) zero_state = 1;
          else if (zero_state === 2) zero_state = 3;
        } else {
          scalars += 1;
        }
      }

      let remaining_width = total_width - width;
      if (height !== null) {
        remaining_width -= scalars * height;
        if (remaining_width < 0) {
          throw new ValueError('incompatible submatrix widths');
        }
        if (remaining_width > 0 && zero_state === 3) {
          throw new ValueError('insufficient information to determine submatrix widths');
        }
        if (remaining_width > 0 && zero_state === 0) {
          throw new ValueError('incompatible submatrix widths');
        }
        row_heights[i] = height;
        zero_widths[i] = remaining_width;
      } else if (zero_state !== 0) {
        throw new ValueError('insufficient information to determine submatrix heights');
      } else if (total_width % R_.length !== 0) {
        throw new ValueError('incompatible submatrix widths');
      } else {
        row_heights[i] = total_width / R_.length;
      }
    }
  }

  return [row_heights as number[], zero_widths, total_width];
}

/**
 * Test whether a non-matrix block is the zero scalar (Python's `not M`).
 */
function _isZeroScalar<R extends RingElement>(M: Matrix<R> | R | number | 0 | undefined): boolean {
  if (M === undefined) return true;
  if (M instanceof Matrix) return false;
  if (typeof M === 'number') return M === 0;
  if (typeof M === 'bigint') return M === 0n;
  const elem = M as unknown as { isZero?: () => boolean };
  return typeof elem.isZero === 'function' ? elem.isZero() : false;
}

/**
 * Return a larger matrix made by concatenating submatrices (rows first, then
 * columns).
 *
 * For example, `block_matrix(R, [[A, B], [C, D]])` creates:
 *
 *     [ A  B ]
 *     [ C  D ]
 *
 * Non-zero scalars are interpreted as square scalar matrices of a size that is
 * deduced from the other blocks, and zeroes as zero matrices of any size.  If
 * the block dimensions cannot be deduced a ValueError is raised.
 *
 * As in SageMath, subdivisions along the block boundaries are set by default;
 * pass `{subdivide: false}` to suppress them.
 *
 * @param ring - The base ring
 * @param blocks - 2D array of blocks, where each block is a Matrix or a scalar
 * @param options - Optional configuration: { subdivide?: boolean } (default: true)
 * @returns A block matrix
 * @see Reference: sage/matrix/special.py:block_matrix
 */
export function block_matrix<R extends RingElement>(
  ring: CoefficientRing<R>,
  blocks: Array<Array<Matrix<R> | R | number | 0>>,
  options?: { subdivide?: boolean }
): Matrix<R> {
  const doSubdivide = options?.subdivide ?? true;

  if (blocks.length === 0) {
    return new Matrix(ring, 0, 0);
  }

  const nBlockRows = blocks.length;
  const nBlockCols = blocks[0]!.length;

  // Validate that all rows have the same number of blocks.  SageMath calls
  // this "list of rows is not valid" and only allows ragged input when
  // subdivide is False (in which case the matrices are fitted row by row);
  // this port always requires a rectangular list of lists.
  for (let i = 1; i < nBlockRows; i++) {
    if (blocks[i]!.length !== nBlockCols) {
      throw new ValueError('list of rows is not valid (rows are wrong types or lengths)');
    }
  }

  let rowHeights: number[] | null = null;
  let colWidths: number[] | null = null;
  let zeroWidths: number[] | null = null;
  let totalWidth: number | null = null;

  // We first try to place the matrices in a rectangular grid
  try {
    [rowHeights, colWidths] = _determine_block_matrix_grid(blocks);
  } catch (e) {
    if (doSubdivide) throw e;
    rowHeights = null;
    colWidths = null;
  }

  if (colWidths === null) {
    // Try placing the matrices in rows instead (only if subdivide is false)
    [rowHeights, zeroWidths, totalWidth] = _determine_block_matrix_rows(blocks);
  }

  const heights = rowHeights!;
  const totalRows = heights.reduce((a, b) => a + b, 0);
  const totalCols =
    colWidths !== null ? colWidths.reduce((a, b) => a + b, 0) : (totalWidth as number);

  const result = new Matrix<R>(ring, totalRows, totalCols);
  const ringCall = ring as unknown as { __call__: (x: unknown) => R };

  let rowOffset = 0;
  for (let i = 0; i < nBlockRows; i++) {
    let colOffset = 0;
    for (let j = 0; j < nBlockCols; j++) {
      const block = blocks[i]![j];
      const blockHeight = heights[i]!;

      if (block instanceof Matrix) {
        for (let bi = 0; bi < block.nrows; bi++) {
          for (let bj = 0; bj < block.ncols; bj++) {
            result.set(rowOffset + bi, colOffset + bj, block.get(bi, bj));
          }
        }
        colOffset += block.ncols;
        continue;
      }

      let blockWidth: number;
      if (_isZeroScalar(block) && zeroWidths !== null) {
        // A zero block soaks up whatever width is left over in this row
        blockWidth = zeroWidths[i]!;
        zeroWidths[i] = 0;
        colOffset += blockWidth;
        continue;
      } else if (zeroWidths !== null) {
        blockWidth = blockHeight;
      } else {
        blockWidth = colWidths![j]!;
      }

      if (!_isZeroScalar(block)) {
        // Non-zero scalar: a scalar matrix, which must be square
        if (blockHeight !== blockWidth) {
          throw new ValueError('nonzero scalar matrix must be square');
        }
        const scalar =
          typeof block === 'number' || typeof block === 'bigint'
            ? ringCall.__call__(block)
            : (block as R);
        for (let k = 0; k < blockHeight; k++) {
          result.set(rowOffset + k, colOffset + k, scalar);
        }
      }
      // zero blocks are already zero

      colOffset += blockWidth;
    }
    rowOffset += heights[i]!;
  }

  if (doSubdivide) {
    const rowLines: number[] = [];
    let acc = 0;
    for (let i = 0; i < heights.length - 1; i++) {
      acc += heights[i]!;
      rowLines.push(acc);
    }
    const colLines: number[] = [];
    acc = 0;
    const widths = colWidths!;
    for (let j = 0; j < widths.length - 1; j++) {
      acc += widths[j]!;
      colLines.push(acc);
    }
    subdivide(result, rowLines, colLines);
  }

  return result;
}

/**
 * Return a block diagonal matrix.
 *
 * Constructs a matrix with the given matrices along the diagonal
 * and zeros elsewhere.
 *
 * As in SageMath (which routes this through `block_matrix`), subdivisions
 * along the block boundaries are set.
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
  const n = sub_matrices.length;
  if (n === 0) {
    return new Matrix(ring, 0, 0);
  }

  // Mirror SageMath: build the n x n grid of blocks with zeros off the
  // diagonal and delegate to block_matrix (which sets the subdivisions).
  const grid: Array<Array<Matrix<R> | R | number | 0>> = [];
  for (let i = 0; i < n; i++) {
    const row: Array<Matrix<R> | R | number | 0> = new Array<Matrix<R> | R | number | 0>(n).fill(0);
    row[i] = sub_matrices[i]!;
    grid.push(row);
  }

  return block_matrix(ring, grid);
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
 * Create a companion matrix from a monic polynomial.
 *
 * The polynomial is given by the list of **all** of its coefficients, with
 * low-degree coefficients first, exactly as in SageMath.  The leading
 * coefficient must be given and must be one; the returned matrix has size
 * equal to the degree of the polynomial, has ones above or below the
 * diagonal, and the **negatives** of the coefficients along the indicated
 * border of the matrix.
 *
 * For `poly = [-2, -3, -4, -5, -6, 1]`, SageMath gives::
 *
 *     'right'          'left'           'bottom'         'top'
 *     [0 0 0 0 2]      [6 1 0 0 0]      [0 1 0 0 0]      [6 5 4 3 2]
 *     [1 0 0 0 3]      [5 0 1 0 0]      [0 0 1 0 0]      [1 0 0 0 0]
 *     [0 1 0 0 4]      [4 0 0 1 0]      [0 0 0 1 0]      [0 1 0 0 0]
 *     [0 0 1 0 5]      [3 0 0 0 1]      [0 0 0 0 1]      [0 0 1 0 0]
 *     [0 0 0 1 6]      [2 0 0 0 0]      [2 3 4 5 6]      [0 0 0 1 0]
 *
 * @param ring - The base ring
 * @param poly - All coefficients `[a_0, a_1, ..., a_n]` of the (monic) polynomial,
 *               low degree first; `a_n` must be one
 * @param format - One of `'right'` (default), `'left'`, `'top'` or `'bottom'`,
 *                 indicating which border holds the negated coefficients
 * @returns The companion matrix, of size `poly.length - 1`
 * @see Reference: sage/matrix/special.py:companion_matrix
 */
export function companion_matrix<R extends RingElement>(
  ring: CoefficientRing<R>,
  poly: R[],
  format: 'right' | 'left' | 'top' | 'bottom' = 'right'
): Matrix<R> {
  if (format !== 'right' && format !== 'left' && format !== 'top' && format !== 'bottom') {
    throw new ValueError(
      `format must be 'right', 'left', 'top' or 'bottom', not ${String(format)}`
    );
  }

  const n = poly.length - 1;
  if (n === -1) {
    throw new ValueError('polynomial cannot be specified by an empty list');
  }
  const leading = poly[n]!;
  if (!leading.eq(ring.one())) {
    throw new ValueError(
      `polynomial (or the polynomial implied by coefficients) must be monic, not a leading coefficient of ${leading.toString()}`
    );
  }

  const entries: R[][] = [];
  for (let i = 0; i < n; i++) {
    entries.push([]);
    for (let j = 0; j < n; j++) {
      entries[i]!.push(ring.zero());
    }
  }

  // 1s below the diagonal, or above the diagonal
  if (format === 'right' || format === 'top') {
    for (let i = 0; i < n - 1; i++) {
      entries[i + 1]![i] = ring.one();
    }
  } else {
    for (let i = 0; i < n - 1; i++) {
      entries[i]![i + 1] = ring.one();
    }
  }

  // right side, left side (reversed), bottom edge, top edge (reversed)
  if (format === 'right') {
    for (let i = 0; i < n; i++) {
      entries[i]![n - 1] = poly[i]!.neg() as R;
    }
  } else if (format === 'left') {
    for (let i = 0; i < n; i++) {
      entries[n - 1 - i]![0] = poly[i]!.neg() as R;
    }
  } else if (format === 'bottom') {
    for (let i = 0; i < n; i++) {
      entries[n - 1]![i] = poly[i]!.neg() as R;
    }
  } else {
    for (let i = 0; i < n; i++) {
      entries[0]![n - 1 - i] = poly[i]!.neg() as R;
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
 * Return a Toeplitz matrix of given first column and first row.
 *
 * In a Toeplitz matrix each descending diagonal from left to right is
 * constant, i.e. `T_{i,j} = T_{i+1,j+1}`.
 *
 * Following SageMath, `r` is the first row **counting from the second
 * column**, so the resulting matrix has `len(c)` rows and `len(r) + 1`
 * columns and the entry `(i, j)` is `c[i-j]` for `i >= j` and
 * `r[j-i-1]` otherwise.  In particular `c[0]` is the diagonal entry and
 * `r[0]` is the entry in position `(0, 1)`.
 *
 * @example
 * ```typescript
 * // matrix.toeplitz([1..4], [5..6]) in SageMath
 * // [1 5 6]
 * // [2 1 5]
 * // [3 2 1]
 * // [4 3 2]
 * ```
 *
 * @param ring - The base ring
 * @param c - First column
 * @param r - First row, starting at the second column
 * @returns The Toeplitz matrix (`c.length` by `r.length + 1`)
 * @see Reference: sage/matrix/special.py:toeplitz
 */
export function toeplitz<R extends RingElement>(
  ring: CoefficientRing<R>,
  c: R[],
  r: R[]
): Matrix<R> {
  const m = c.length;
  const n = r.length + 1;

  if (m === 0) {
    return new Matrix(ring, 0, n);
  }

  const entries: R[][] = [];
  for (let i = 0; i < m; i++) {
    entries.push([]);
    for (let j = 0; j < n; j++) {
      if (i >= j) {
        // Below or on the diagonal: use the column c
        entries[i]!.push(c[i - j]!);
      } else {
        // Above the diagonal: use the row r, which starts at column 1
        entries[i]!.push(r[j - i - 1]!);
      }
    }
  }

  return new Matrix(ring, m, n, entries);
}

/**
 * Return a Hankel matrix of given first column and last row.
 *
 * A Hankel matrix is constant along anti-diagonals: `H_{ij} = v_{i+j}` where
 * `v_i = c_i` for `i < len(c)` and `v_{len(c)+i} = r_i` otherwise.
 *
 * Following SageMath, `r` is the last row **from the second to the last
 * column**, so the resulting matrix has `len(c)` rows and `len(r) + 1`
 * columns.  If `r` is omitted it defaults to `len(c) - 1` zeros, which makes
 * the matrix square with zeros below the first anti-diagonal.
 *
 * @example
 * ```typescript
 * // matrix.hankel([1..3], [7..10]) in SageMath
 * // [ 1  2  3  7  8]
 * // [ 2  3  7  8  9]
 * // [ 3  7  8  9 10]
 * ```
 *
 * @param ring - The base ring
 * @param c - First column
 * @param r - Last row, from the second to the last column (optional)
 * @returns The Hankel matrix (`c.length` by `r.length + 1`)
 * @see Reference: sage/matrix/special.py:hankel
 */
export function hankel<R extends RingElement>(
  ring: CoefficientRing<R>,
  c: R[],
  r?: R[]
): Matrix<R> {
  const m = c.length;

  // Default last row: m - 1 zeros, which makes the matrix square
  const actualR: R[] = r ?? new Array<R>(Math.max(0, m - 1)).fill(ring.zero());
  const n = actualR.length;

  const entries: R[][] = [];
  for (let i = 0; i < m; i++) {
    entries.push([]);
    for (let j = 0; j <= n; j++) {
      const idx = i + j;
      entries[i]!.push(idx < m ? c[idx]! : actualR[idx - m]!);
    }
  }

  return new Matrix(ring, m, n + 1, entries);
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
 * Determine the field over which to run Berlekamp-Massey.
 *
 * Mirrors SageMath's `K = a[0].parent().fraction_field()`, with the
 * `AttributeError` fallback to the rational field: entries that carry a
 * `parent` (finite field elements, polynomials, ...) contribute their parent's
 * fraction field when it has one, and plain integers fall back to QQ.
 *
 * @see Reference: sage/matrix/berlekamp_massey.py:berlekamp_massey
 */
function _coefficient_field_of(
  a: Array<RingElement | number | bigint>
): CoefficientRing<RingElement> {
  const first = a[0];
  if (first === undefined || typeof first === 'number' || typeof first === 'bigint') {
    return QQ as unknown as CoefficientRing<RingElement>;
  }
  const parent = (first as unknown as { parent?: unknown }).parent;
  if (parent === undefined || parent === null) {
    return QQ as unknown as CoefficientRing<RingElement>;
  }
  const p = parent as {
    fraction_field?: () => CoefficientRing<RingElement>;
  };
  if (typeof p.fraction_field === 'function') {
    return p.fraction_field();
  }
  return parent as CoefficientRing<RingElement>;
}

/**
 * Use the Berlekamp-Massey algorithm to find the minimal polynomial
 * of a linear recurrence sequence.
 *
 * The minimal polynomial of a linear recurrence {a_r} is the unique monic
 * polynomial g such that if {a_r} satisfies a linear recurrence
 * a_{j+k} + b_{j-1} * a_{j-1+k} + ... + b_0 * a_k = 0 (for all k >= 0),
 * then g divides x^j + sum_{i=0}^{j-1} b_i * x^i.
 *
 * The sequence entries may be ring elements or plain integers; as in SageMath
 * the computation is carried out over the fraction field of their parent (QQ
 * for plain integers).
 *
 * @example
 * ```typescript
 * berlekamp_massey([1, 2, 1, 2, 1, 2]).toString();          // 'x^2 - 1'
 * berlekamp_massey([F7(1), F7(19), F7(1), F7(19)]);         // x^2 + 6 over GF(7)
 * ```
 *
 * @param a - List of even length of elements of a field (or of integers)
 * @returns The minimal polynomial of the sequence as a Polynomial
 * @see Reference: sage/matrix/berlekamp_massey.py:berlekamp_massey
 */
export function berlekamp_massey<R extends RingElement>(
  a: Array<R | number | bigint>
): Polynomial<R> {
  if (!Array.isArray(a)) {
    throw new TypeError('argument must be a list or tuple');
  }

  if (a.length % 2 !== 0) {
    throw new ValueError('argument must have an even number of terms');
  }

  const M = Math.floor(a.length / 2);

  // Determine the coefficient field, mirroring SageMath's
  // ``K = a[0].parent().fraction_field()`` with a fallback to QQ when the
  // entries are plain integers (which have no parent here).
  const ring = _coefficient_field_of(a) as CoefficientRing<R>;
  const polyRing = new PolynomialRing(ring, 'x');

  // Coerce the sequence into K, as SageMath's ``R(a)`` does
  const coeffs: R[] = a.map((x) =>
    typeof x === 'number' || typeof x === 'bigint' ? ring.__call__(x) : (x as R)
  );

  // Create polynomial from sequence: f0 = a[0] + a[1]*x + ... + a[2M-1]*x^{2M-1}
  const f0 = new Polynomial<R>(coeffs, polyRing);

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

  // Check if matrix is 0-1 and count ones (needed for complement logic).
  // SageMath starts the count at one; keep that so that the automatic
  // use_complement threshold matches exactly.
  let isZ2 = true;
  let numOnes = 1;
  let badEntry: { x: R; i: number; j: number } | null = null;
  outer: for (let i = 0; i < m; i++) {
    for (let j = 0; j < n; j++) {
      const x = matrix.get(i, j);
      if (!x.isZero()) {
        if (!x.eq(one)) {
          isZ2 = false;
          badEntry = { x, i, j };
          break outer;
        }
        numOnes++;
      }
    }
  }

  // Validate complement usage
  if (complement && !isZ2) {
    throw new ValueError(
      `coefficients must be zero or one, but we have '${badEntry!.x.toString()}' in position (${badEntry!.i},${badEntry!.j}).`
    );
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
    throw new ValueError(
      `algorithm must be one of "ButeraPernici", "Ryser", or "naive", got "${algorithm}"`
    );
  }

  // Apply inclusion-exclusion if computing complement.
  // The coefficients grow like C(m,k)*C(n,k)*k!, far past 2^53, so they are
  // computed in BigInt with Python's floor-division semantics.
  if (complement) {
    const a: R[] = [one];
    const mB = BigInt(m);
    const nB = BigInt(n);
    let c1 = 1n;
    for (let k = 1; k <= mn; k++) {
      const kB = BigInt(k);
      // c1 = C(m, k) * C(n, k) * k! / C(m, k-1) / C(n, k-1) / (k-1)!
      //    = (m-k+1) * (n-k+1) / k
      c1 = _floorDiv(c1 * (mB - kB + 1n) * (nB - kB + 1n), kB);
      let c = c1;
      // s = c * b[0] + (-1)^k * b[k]
      let s = _scalarMul(ring, c, b[0]!);
      const sign_k = k % 2 === 0 ? 1n : -1n;
      s = s.add(_scalarMul(ring, sign_k, b[k]!)) as R;

      for (let j = 1; j < k; j++) {
        // c = -c * (k-j+1) // ((m-j+1) * (n-j+1))
        const jB = BigInt(j);
        c = _floorDiv(-c * (kB - jB + 1n), (mB - jB + 1n) * (nB - jB + 1n));
        s = s.add(_scalarMul(ring, c, b[j]!)) as R;
      }
      a.push(s);
    }
    return a;
  }

  return b;
}

/**
 * Floor division on BigInt, matching Python's `//` operator (which rounds
 * towards negative infinity, unlike BigInt's truncating `/`).
 */
function _floorDiv(a: bigint, b: bigint): bigint {
  const q = a / b;
  return a % b !== 0n && a < 0n !== b < 0n ? q - 1n : q;
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

  // p is a dictionary mapping bitmasks to polynomial coefficients.
  // The bitmask represents which eta_j variables are present; SageMath uses
  // Python bignums for these masks, so they must be BigInt here -- a 32-bit
  // JS `1 << j` would alias column 32 onto column 0.
  // The value is an array of coefficients [c_0, c_1, ..., c_k] for polynomial c_0 + c_1*t + ... + c_k*t^k
  let p: Map<bigint, R[]> = new Map();
  p.set(0n, [ring.one()]); // Start with 1

  // Track which columns still have nonzero entries in remaining rows
  const varsToDo = new Set<number>();
  for (let j = 0; j < n; j++) {
    varsToDo.add(j);
  }

  for (let i = 0; i < m; i++) {
    const a = rows[i]!;

    // Build p1 = 1 + t * sum_j A[i,j] * eta_j
    const p1: Map<bigint, R[]> = new Map();
    p1.set(0n, [ring.one()]); // Constant term 1

    for (let j = 0; j < n; j++) {
      if (!a[j]!.isZero()) {
        // eta_j corresponds to bit j, coefficient is A[i,j] * t (degree 1)
        p1.set(1n << BigInt(j), [ring.zero(), a[j]!]);
      }
    }

    // Determine which variables can be "integrated" (set to 1)
    // A variable eta_j can be integrated if it doesn't appear in any remaining row
    let maskFree = 0n;
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
        maskFree |= 1n << BigInt(j);
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

  const coeffs = p.get(0n);
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
  p1: Map<bigint, R[]>,
  p2: Map<bigint, R[]>,
  maskFree: bigint,
  prec: number
): Map<bigint, R[]> {
  const result: Map<bigint, R[]> = new Map();

  for (const [exp1, v1] of p1) {
    if (_polyIsZero(v1)) continue;
    for (const [exp2, v2] of p2) {
      // Skip if monomials share any variables (product is 0 due to nilpotency)
      if ((exp1 & exp2) !== 0n) continue;

      // Multiply the polynomial coefficients
      const prod = _polyMul(ring, v1, v2, prec);
      if (_polyIsZero(prod)) continue;

      // Combine exponents and integrate free variables
      const exp = (exp1 | exp2) ^ ((exp1 | exp2) & maskFree);

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
function _polyMul<R extends RingElement>(
  ring: CoefficientRing<R>,
  a: R[],
  b: R[],
  prec: number
): R[] {
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
 * Multiply a ring element by an integer scalar using double-and-add.
 *
 * The inclusion-exclusion coefficients used by {@link rook_vector} are of the
 * order of `C(m,k) * C(n,k) * k!`, so repeated addition is not an option: this
 * runs in O(log |scalar|) ring additions.
 */
function _scalarMul<R extends RingElement>(ring: CoefficientRing<R>, scalar: bigint, elem: R): R {
  if (scalar === 0n) return ring.zero();
  if (scalar === 1n) return elem;
  if (scalar === -1n) return elem.neg() as R;

  let k = scalar < 0n ? -scalar : scalar;
  let acc: R | null = null;
  let addend = elem;
  while (k > 0n) {
    if (k & 1n) {
      acc = acc === null ? addend : (acc.add(addend) as R);
    }
    k >>= 1n;
    if (k > 0n) {
      addend = addend.add(addend) as R;
    }
  }
  const result = acc ?? ring.zero();
  return scalar > 0n ? result : (result.neg() as R);
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
    let count = 0n;
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
        count += 1n;
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

  // Check that the entries support substitution.  SageMath calls the entry's
  // `subs`; univariate polynomials in this port expose the same operation as
  // `evaluate`, so that is accepted too.
  type Substitutable = {
    subs?: (arg: unknown) => R;
    evaluate?: (arg: unknown) => R;
  };
  const applySubs = (entry: Substitutable): R => {
    if (typeof entry.subs === 'function') return entry.subs(substitutions);
    if (typeof entry.evaluate === 'function') return entry.evaluate(substitutions);
    throw new TypeError('subs not defined for elements of the base ring');
  };

  if (m > 0 && n > 0) {
    const firstEntry = matrix.get(0, 0) as unknown as Substitutable;
    if (typeof firstEntry.subs !== 'function' && typeof firstEntry.evaluate !== 'function') {
      throw new TypeError('subs not defined for elements of the base ring');
    }
  }

  const entries: R[][] = [];
  for (let i = 0; i < m; i++) {
    entries.push([]);
    for (let j = 0; j < n; j++) {
      entries[i]!.push(applySubs(matrix.get(i, j) as unknown as Substitutable));
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

  // The denominator of an entry: elements of QQ expose it as a bigint-valued
  // property, other rings may expose it as a method returning a ring element.
  const denomOf = (x: unknown): bigint => {
    const d = (x as { denominator?: unknown }).denominator;
    const raw = typeof d === 'function' ? (d as () => unknown).call(x) : d;
    if (typeof raw === 'bigint') return raw;
    if (typeof raw === 'number') return BigInt(raw);
    if (raw !== undefined && raw !== null) {
      const v = (raw as { value?: unknown }).value;
      if (typeof v === 'bigint') return v;
      const parsed = BigInt(String(raw));
      return parsed;
    }
    throw new TypeError('denominator not defined for elements of the base ring');
  };

  const gcdB = (a: bigint, b: bigint): bigint => {
    let x = a < 0n ? -a : a;
    let y = b < 0n ? -b : b;
    while (y !== 0n) {
      [x, y] = [y, x % y];
    }
    return x;
  };

  let result = 1n;
  for (let i = 0; i < m; i++) {
    for (let j = 0; j < n; j++) {
      const d = denomOf(matrix.get(i, j));
      if (d === 0n) continue;
      result = (result / gcdB(result, d)) * (d < 0n ? -d : d);
    }
  }

  return ring.__call__(result);
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
      if (density >= 1.0 || _randomFraction() < density) {
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

// ---------------------------------------------------------------------------
// Element comparison helpers
//
// SageMath compares matrix entries by *value*.  Ring elements in this port do
// not share a comparison interface, so these helpers recover the numeric value
// of an element (the lift, for residues and finite fields) and fall back to the
// element's own `cmp` where it has one.
// ---------------------------------------------------------------------------

/**
 * Return the integer value of a ring element, or null if it has none.
 *
 * Handles residues and prime finite field elements (which store a `value`
 * bigint) and finite field extension elements (whose `lift` is a polynomial
 * over the prime field, read as an integer in base p, matching SageMath's
 * `integer_representation`).
 */
function _elementValue<R extends RingElement>(x: R): bigint | null {
  const withValue = x as unknown as { value?: unknown };
  if (typeof withValue.value === 'bigint') {
    return withValue.value;
  }

  const withLift = x as unknown as {
    lift?: { coeffs?: ReadonlyArray<{ value?: unknown }> };
    parent?: { characteristic?: unknown };
  };
  const coeffs = withLift.lift?.coeffs;
  const p = withLift.parent?.characteristic;
  if (Array.isArray(coeffs) && typeof p === 'bigint') {
    let acc = 0n;
    for (let i = coeffs.length - 1; i >= 0; i--) {
      const c = coeffs[i]!.value;
      if (typeof c !== 'bigint') return null;
      acc = acc * p + c;
    }
    return acc;
  }

  return null;
}

/**
 * Compare two ring elements by value, returning a negative number, zero or a
 * positive number as `a` is smaller than, equal to, or greater than `b`.
 */
function _compareElements<R extends RingElement>(a: R, b: R): number {
  if (a.eq(b)) return 0;

  const withCmp = a as unknown as { cmp?: (other: R) => number };
  if (typeof withCmp.cmp === 'function') {
    return withCmp.cmp(b);
  }

  const va = _elementValue(a);
  const vb = _elementValue(b);
  if (va !== null && vb !== null) {
    return va < vb ? -1 : va > vb ? 1 : 0;
  }

  // Last resort: a stable, deterministic order on the printed representation
  const as = a.toString();
  const bs = b.toString();
  return as < bs ? -1 : as > bs ? 1 : 0;
}

/**
 * A key identifying an element up to equality, for use in multiset comparisons.
 */
function _elementKey<R extends RingElement>(x: R): string {
  const v = _elementValue(x);
  return v !== null ? v.toString() : x.toString();
}

/**
 * Compare two rows of ring elements lexicographically.
 */
function _compareRows<R extends RingElement>(a: R[], b: R[]): number {
  const len = Math.min(a.length, b.length);
  for (let i = 0; i < len; i++) {
    const c = _compareElements(a[i]!, b[i]!);
    if (c !== 0) return c;
  }
  return a.length - b.length;
}

// ---------------------------------------------------------------------------
// Exact search for a row/column permutation between two matrices
// ---------------------------------------------------------------------------

/**
 * Search for permutations `rowPerm`, `colPerm` with
 * `A[rowPerm[i]][colPerm[j]] === B[i][j]` for all `i, j`, where the entries are
 * given as comparison keys.
 *
 * This is the combinatorial core of {@link is_permutation_of}.  SageMath tests
 * the same property by deciding isomorphism of the edge-labelled bipartite
 * graphs of the two matrices; here it is decided by backtracking over the row
 * assignment with the column multisets of every prefix used to prune, which is
 * complete (it never reports a false negative).
 *
 * @returns The pair of permutations, or null if no such pair exists
 * @see Reference: sage/matrix/matrix2.pyx:is_permutation_of
 */
function _findPermutationOfKeys(A: string[][], B: string[][]): [number[], number[]] | null {
  const m = A.length;
  if (m !== B.length) return null;
  const n = m === 0 ? 0 : A[0]!.length;
  if (m === 0 || n === 0) {
    return [Array.from({ length: m }, (_, i) => i), Array.from({ length: n }, (_, j) => j)];
  }
  if (B[0]!.length !== n) return null;

  // Rows can only be matched to rows with the same multiset of entries, and
  // likewise for columns; both are necessary conditions and cheap to test.
  const sig = (v: string[]): string => [...v].sort().join(' ');
  const sigA = A.map(sig);
  const sigB = B.map(sig);
  if (sig(sigA) !== sig(sigB)) return null;

  const colOf = (M: string[][], j: number): string[] => M.map((row) => row[j]!);
  const colSigA: string[] = [];
  const colSigB: string[] = [];
  for (let j = 0; j < n; j++) {
    colSigA.push(sig(colOf(A, j)));
    colSigB.push(sig(colOf(B, j)));
  }
  if (sig(colSigA) !== sig(colSigB)) return null;

  // Assign the rows of B whose signature is rarest first: it keeps the search
  // tree narrow near the root.
  const countA = new Map<string, number>();
  for (const s of sigA) countA.set(s, (countA.get(s) ?? 0) + 1);
  const order = Array.from({ length: m }, (_, i) => i).sort(
    (i1, i2) => (countA.get(sigB[i1]!) ?? 0) - (countA.get(sigB[i2]!) ?? 0)
  );

  const rowPerm: number[] = new Array<number>(m).fill(-1);
  const used: boolean[] = new Array<boolean>(m).fill(false);

  // Column keys of B for the prefix of depth k, sorted; precomputed per depth
  const colKeysB: string[][] = [];
  {
    const running: string[] = new Array<string>(n).fill('');
    for (let k = 0; k < m; k++) {
      const bRow = B[order[k]!]!;
      for (let j = 0; j < n; j++) {
        running[j] = `${running[j]!} ${bRow[j]!}`;
      }
      colKeysB.push([...running].sort());
    }
  }

  const runningA: string[] = new Array<string>(n).fill('');

  const search = (depth: number): boolean => {
    if (depth === m) return true;
    const target = order[depth]!;
    const saved = [...runningA];
    for (let r = 0; r < m; r++) {
      if (used[r]) continue;
      if (sigA[r] !== sigB[target]) continue;
      for (let j = 0; j < n; j++) {
        runningA[j] = `${saved[j]!} ${A[r]![j]!}`;
      }
      // Prune: the column multisets of the two prefixes must agree
      const sortedA = [...runningA].sort();
      const targetB = colKeysB[depth]!;
      let ok = true;
      for (let j = 0; j < n; j++) {
        if (sortedA[j] !== targetB[j]) {
          ok = false;
          break;
        }
      }
      if (ok) {
        rowPerm[target] = r;
        used[r] = true;
        if (search(depth + 1)) return true;
        used[r] = false;
        rowPerm[target] = -1;
      }
      for (let j = 0; j < n; j++) {
        runningA[j] = saved[j]!;
      }
    }
    return false;
  };

  if (!search(0)) return null;

  // The full column keys agree as multisets, so read off the bijection
  const fullA: string[] = new Array<string>(n).fill('');
  const fullB: string[] = new Array<string>(n).fill('');
  for (let i = 0; i < m; i++) {
    for (let j = 0; j < n; j++) {
      fullA[j] = `${fullA[j]!} ${A[rowPerm[i]!]![j]!}`;
      fullB[j] = `${fullB[j]!} ${B[i]![j]!}`;
    }
  }
  const colPerm: number[] = new Array<number>(n).fill(-1);
  const usedCol: boolean[] = new Array<boolean>(n).fill(false);
  for (let j = 0; j < n; j++) {
    let found = -1;
    for (let c = 0; c < n; c++) {
      if (!usedCol[c] && fullA[c] === fullB[j]) {
        found = c;
        break;
      }
    }
    if (found === -1) return null;
    usedCol[found] = true;
    colPerm[j] = found;
  }

  return [rowPerm, colPerm];
}

/**
 * Convert a matrix into the array of comparison keys used by
 * {@link _findPermutationOfKeys}.
 */
function _matrixKeys<R extends RingElement>(M: Matrix<R>): string[][] {
  const keys: string[][] = [];
  for (let i = 0; i < M.nrows; i++) {
    const row: string[] = [];
    for (let j = 0; j < M.ncols; j++) {
      row.push(_elementKey(M.get(i, j)));
    }
    keys.push(row);
  }
  return keys;
}

// ---------------------------------------------------------------------------
// Permutation normal form
// ---------------------------------------------------------------------------

/**
 * A partial state of the search for the permutation normal form: the rows
 * placed so far, the rows still to place, the current column order, and the
 * sizes of the blocks of columns that are still interchangeable.
 */
interface _PNFState {
  placed: number[];
  remaining: number[];
  cols: number[];
  blockSizes: number[];
}

/**
 * Arrange row `r` of `matrix` inside the column blocks of `state`, sorting each
 * block descending, and return the resulting row values together with the
 * refined state.
 */
function _pnfExtend<R extends RingElement>(
  matrix: Matrix<R>,
  state: _PNFState,
  r: number
): { vals: R[]; next: _PNFState } {
  const cols: number[] = [];
  const blockSizes: number[] = [];
  const vals: R[] = [];

  let pos = 0;
  for (const size of state.blockSizes) {
    const block = state.cols.slice(pos, pos + size);
    pos += size;
    block.sort((c1, c2) => _compareElements(matrix.get(r, c2), matrix.get(r, c1)));
    // Split the block into runs of equal value; those columns stay interchangeable
    let runStart = 0;
    for (let k = 1; k <= block.length; k++) {
      if (k === block.length || !matrix.get(r, block[k]!).eq(matrix.get(r, block[runStart]!))) {
        blockSizes.push(k - runStart);
        runStart = k;
      }
    }
    for (const c of block) {
      cols.push(c);
      vals.push(matrix.get(r, c));
    }
  }

  return {
    vals,
    next: {
      placed: [...state.placed, r],
      remaining: state.remaining.filter((x) => x !== r),
      cols,
      blockSizes,
    },
  };
}

/**
 * Build the key matrix used to decide whether two search states are equivalent.
 *
 * The first row records the block a column belongs to (SageMath stores the same
 * information in the sentinel row `S`), which forces any permutation matching
 * two states to preserve the column blocks.
 */
function _pnfStateKeys<R extends RingElement>(matrix: Matrix<R>, state: _PNFState): string[][] {
  const blockRow: string[] = [];
  for (let b = 0; b < state.blockSizes.length; b++) {
    for (let k = 0; k < state.blockSizes[b]!; k++) {
      blockRow.push(`#${b}`);
    }
  }
  const keys: string[][] = [blockRow];
  for (const r of state.remaining) {
    keys.push(state.cols.map((c) => _elementKey(matrix.get(r, c))));
  }
  return keys;
}

/**
 * Take the set of matrices that are `matrix` permuted by any row and column
 * permutation, and return the maximal one of the set where matrices are ordered
 * lexicographically going along each row.
 *
 * The maximal matrix is found row by row: the first row must be the largest row
 * of the matrix sorted decreasingly, which pins the column order up to the
 * blocks of columns carrying equal entries; each further row is then chosen to
 * maximise its arrangement inside those blocks, refining them.  Whenever
 * several choices tie they are all kept, after discarding the ones that are
 * equivalent under a block-preserving permutation -- this is the structure of
 * SageMath's algorithm, whose sentinel row `S` encodes exactly those blocks.
 *
 * @param matrix - The matrix
 * @param check - If true, return a tuple (maximal_matrix, (row_perm, col_perm))
 *                with `maximal[i][j] === matrix[row_perm[i]][col_perm[j]]`
 * @returns The permutation normal form (or tuple if check=true)
 * @see Reference: sage/matrix/matrix2.pyx:permutation_normal_form
 */
export function permutation_normal_form<R extends RingElement>(
  matrix: Matrix<R>,
  check?: boolean
): Matrix<R> | [Matrix<R>, [number[], number[]]] {
  const m = matrix.nrows;
  const n = matrix.ncols;

  if (m === 0 || n === 0) {
    if (check) {
      return [
        matrix.copy(),
        [Array.from({ length: m }, (_, i) => i), Array.from({ length: n }, (_, j) => j)],
      ];
    }
    return matrix.copy();
  }

  let states: _PNFState[] = [
    {
      placed: [],
      remaining: Array.from({ length: m }, (_, i) => i),
      cols: Array.from({ length: n }, (_, j) => j),
      blockSizes: [n],
    },
  ];

  for (let l = 0; l < m; l++) {
    let best: R[] | null = null;
    let winners: _PNFState[] = [];

    for (const state of states) {
      for (const r of state.remaining) {
        const { vals, next } = _pnfExtend(matrix, state, r);
        const c = best === null ? 1 : _compareRows(vals, best);
        if (c > 0) {
          best = vals;
          winners = [next];
        } else if (c === 0) {
          winners.push(next);
        }
      }
    }

    // Discard states that are equivalent under a block-preserving permutation
    const kept: _PNFState[] = [];
    const keptKeys: string[][][] = [];
    for (const w of winners) {
      const keys = _pnfStateKeys(matrix, w);
      let duplicate = false;
      for (const other of keptKeys) {
        if (_findPermutationOfKeys(keys, other) !== null) {
          duplicate = true;
          break;
        }
      }
      if (!duplicate) {
        kept.push(w);
        keptKeys.push(keys);
      }
    }
    states = kept;
  }

  // Every surviving state realises the same (maximal) matrix
  const final = states[0]!;
  const result = _applyPermutation(matrix, final.placed, final.cols);

  if (check) {
    return [result, [final.placed, final.cols]];
  }
  return result;
}

/**
 * Apply row and column permutation to a matrix.
 *
 * The result satisfies `out[i][j] === matrix[rowPerm[i]][colPerm[j]]`.
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
 * Return true if there exists a permutation of rows and columns sending `A` to
 * `B`, and false otherwise.
 *
 * @param A - First matrix
 * @param B - Second matrix
 * @param check - If true, return `[boolean, permutation]` where permutation is
 *                `[row_perm, col_perm]` with
 *                `A[row_perm[i]][col_perm[j]] === B[i][j]`, or null
 * @returns True if A is a permutation of B (or a tuple if check=true)
 * @see Reference: sage/matrix/matrix2.pyx:is_permutation_of
 */
export function is_permutation_of<R extends RingElement>(
  A: Matrix<R>,
  B: Matrix<R>,
  check?: boolean
): boolean | [boolean, [number[], number[]] | null] {
  // Different dimensions means definitely not a permutation
  if (A.nrows !== B.nrows || A.ncols !== B.ncols) {
    return check ? [false, null] : false;
  }

  const perm = _findPermutationOfKeys(_matrixKeys(A), _matrixKeys(B));

  if (check) {
    return perm === null ? [false, null] : [true, perm];
  }
  return perm !== null;
}
