/**
 * @module sage/modules/bkz
 * @description Block Korkine-Zolotarev (BKZ) lattice reduction algorithm
 *
 * This module implements the BKZ lattice reduction algorithm, which produces
 * higher quality lattice basis reductions than LLL at the cost of increased
 * computation time. BKZ is parameterized by a block size beta, where larger
 * block sizes produce better reductions.
 *
 * Port of: sage/matrix/matrix_integer_dense.pyx (BKZ method)
 *
 * @see Reference: C.P. Schnorr and M. Euchner, "Lattice Basis Reduction:
 *                 Improved Practical Algorithms and Solving Subset Sum Problems"
 * @see Reference: sage/matrix/matrix_integer_dense.pyx:2779 (BKZ)
 */

import { ValueError } from '../errors.js';
import { IntegerMatrix, IntegerMatrixFromEntries } from '../matrix/matrix_integer.js';
import { type GramSchmidtResult, gramSchmidt, lllReduce } from './free_module_integer.js';

/**
 * Options for the BKZ algorithm.
 *
 * @see Reference: sage/matrix/matrix_integer_dense.pyx:2779
 */
export interface BKZOptions {
  /**
   * The delta parameter for LLL reduction within BKZ.
   * Must satisfy 0.25 < delta <= 1.0.
   * Default: 0.99
   */
  delta?: number;

  /**
   * The block size for BKZ reduction.
   * Must satisfy 2 <= blockSize <= rank.
   * Larger block sizes produce better reductions but take longer.
   * Default: 10
   */
  blockSize?: number;

  /**
   * Maximum number of BKZ tours (iterations over all blocks).
   * If not specified, continues until no improvement.
   * Default: 0 (unlimited)
   */
  maxLoops?: number;

  /**
   * Precision mode for floating-point arithmetic.
   * 'double' uses standard JavaScript numbers.
   * 'arbitrary' would use arbitrary precision (not yet implemented).
   * Default: 'double'
   */
  precision?: 'double' | 'arbitrary';

  /**
   * Pruning coefficients for enumeration (advanced).
   * If specified, should be an array of length blockSize with values in (0, 1].
   */
  pruning?: number[];

  /**
   * Whether to automatically abort if no progress is made.
   * Default: true
   */
  autoAbort?: boolean;
}

/**
 * Result of the BKZ algorithm.
 */
export interface BKZResult {
  /** The BKZ-reduced basis */
  basis: IntegerMatrix;
  /** Number of tours (complete passes) performed */
  tours: number;
  /** Whether the algorithm converged (no more improvements possible) */
  converged: boolean;
}

/**
 * Perform BKZ (Block Korkine-Zolotarev) lattice basis reduction.
 *
 * BKZ reduces a lattice basis to a form that is "more reduced" than LLL.
 * The algorithm works by:
 * 1. Starting with an LLL-reduced basis
 * 2. For each block of consecutive vectors of size blockSize:
 *    - Find the shortest vector in the projected sublattice (SVP)
 *    - If a shorter vector is found, insert it and re-LLL-reduce
 * 3. Repeat until no improvement is found
 *
 * The quality of BKZ reduction improves with larger block sizes:
 * - blockSize = 2: Equivalent to LLL
 * - blockSize = rank: Equivalent to HKZ (Hermite-Korkine-Zolotarev)
 *
 * @param basis - An IntegerMatrix whose rows form the lattice basis
 * @param options - BKZ parameters
 * @returns A BKZ-reduced IntegerMatrix
 *
 * @throws {ValueError} If delta is not in (0.25, 1.0]
 * @throws {ValueError} If blockSize is not in [2, rank]
 *
 * @example
 * ```typescript
 * const basis = IntegerMatrixFromEntries([
 *   [1n, 0n, 3n],
 *   [0n, 2n, 1n],
 *   [0n, 2n, 7n],
 * ]);
 *
 * // BKZ-20 reduction
 * const reduced = BKZ(basis, { blockSize: 20 });
 * ```
 *
 * @see Reference: sage/matrix/matrix_integer_dense.pyx:2779
 */
export function BKZ(basis: IntegerMatrix, options?: BKZOptions): IntegerMatrix {
  const delta = options?.delta ?? 0.99;
  const maxLoops = options?.maxLoops ?? 0;
  const autoAbort = options?.autoAbort ?? true;

  const n = basis.nrows;
  const m = basis.ncols;

  // Validate delta parameter first (before checking rank-dependent parameters)
  if (delta <= 0.25 || delta > 1.0) {
    throw new ValueError(`delta must be in (0.25, 1.0], got ${delta}`);
  }

  // Handle trivial cases before blockSize validation
  if (n === 0) {
    return new IntegerMatrix(0, m);
  }
  if (n === 1) {
    return basis.copy();
  }

  // Default blockSize is min(10, rank) to handle small matrices
  const blockSize = options?.blockSize ?? Math.min(10, n);

  // Validate blockSize parameter
  if (blockSize < 2) {
    throw new ValueError(`blockSize must be >= 2, got ${blockSize}`);
  }
  if (blockSize > n) {
    throw new ValueError(`blockSize must be <= rank (${n}), got ${blockSize}`);
  }

  // Step 1: Start with LLL-reduced basis
  let B = lllReduce(basis, { delta, eta: 0.501 });

  // If blockSize <= 2, LLL is sufficient
  if (blockSize <= 2) {
    return B;
  }

  // Main BKZ loop
  let tours = 0;
  let improved = true;
  const maxTours = maxLoops > 0 ? maxLoops : n * n * 10; // Safety limit

  while (improved && tours < maxTours) {
    improved = false;
    tours++;

    // Process each block starting position
    for (let kappa = 0; kappa <= n - 2; kappa++) {
      // Determine block size for this position
      const beta = Math.min(blockSize, n - kappa);

      if (beta < 2) {
        continue;
      }

      // Compute Gram-Schmidt for the current basis
      const gs = computeGramSchmidtFromMatrix(B);

      // Find shortest vector in the projected block [kappa, kappa + beta)
      const svpResult = enumerateSVP(B, gs, kappa, beta, options?.pruning);

      if (svpResult !== null) {
        const { coefficients, normSquared } = svpResult;

        // Check if we found a shorter vector than b*[kappa]
        const currentNormSq = gs.B[kappa]!;

        if (normSquared < currentNormSq * (1 - 1e-10)) {
          // Construct the new vector
          const newVector = constructVector(B, coefficients, kappa, beta);

          // Insert the new vector and re-reduce
          B = insertAndReduce(B, newVector, kappa, delta);

          improved = true;
        }
      }
    }

    // Auto-abort if no improvement in this tour
    if (autoAbort && !improved) {
      break;
    }
  }

  return B;
}

/**
 * Perform BKZ reduction and return detailed result.
 *
 * @param basis - The lattice basis matrix
 * @param options - BKZ parameters
 * @returns BKZResult with the reduced basis and metadata
 */
export function BKZWithInfo(basis: IntegerMatrix, options?: BKZOptions): BKZResult {
  const delta = options?.delta ?? 0.99;
  const maxLoops = options?.maxLoops ?? 0;
  const autoAbort = options?.autoAbort ?? true;

  const n = basis.nrows;
  const m = basis.ncols;

  // Validate delta first
  if (delta <= 0.25 || delta > 1.0) {
    throw new ValueError(`delta must be in (0.25, 1.0], got ${delta}`);
  }

  // Handle trivial cases before blockSize validation
  if (n === 0) {
    return { basis: new IntegerMatrix(0, m), tours: 0, converged: true };
  }
  if (n === 1) {
    return { basis: basis.copy(), tours: 0, converged: true };
  }

  // Default blockSize is min(10, rank) to handle small matrices
  const blockSize = options?.blockSize ?? Math.min(10, n);

  // Validate blockSize
  if (blockSize < 2) {
    throw new ValueError(`blockSize must be >= 2, got ${blockSize}`);
  }
  if (blockSize > n) {
    throw new ValueError(`blockSize must be <= rank (${n}), got ${blockSize}`);
  }

  // Start with LLL-reduced basis
  let B = lllReduce(basis, { delta, eta: 0.501 });

  if (blockSize <= 2) {
    return { basis: B, tours: 0, converged: true };
  }

  // Main BKZ loop
  let tours = 0;
  let improved = true;
  let converged = false;
  const maxTours = maxLoops > 0 ? maxLoops : n * n * 10;

  while (improved && tours < maxTours) {
    improved = false;
    tours++;

    for (let kappa = 0; kappa <= n - 2; kappa++) {
      const beta = Math.min(blockSize, n - kappa);

      if (beta < 2) {
        continue;
      }

      const gs = computeGramSchmidtFromMatrix(B);
      const svpResult = enumerateSVP(B, gs, kappa, beta, options?.pruning);

      if (svpResult !== null) {
        const { coefficients, normSquared } = svpResult;
        const currentNormSq = gs.B[kappa]!;

        if (normSquared < currentNormSq * (1 - 1e-10)) {
          const newVector = constructVector(B, coefficients, kappa, beta);
          B = insertAndReduce(B, newVector, kappa, delta);
          improved = true;
        }
      }
    }

    if (autoAbort && !improved) {
      converged = true;
      break;
    }
  }

  return { basis: B, tours, converged };
}

/**
 * Compute Gram-Schmidt orthogonalization from an IntegerMatrix.
 *
 * @param M - The matrix
 * @returns GramSchmidtResult
 */
function computeGramSchmidtFromMatrix(M: IntegerMatrix): GramSchmidtResult {
  const n = M.nrows;
  const m = M.ncols;
  const rows: number[][] = [];

  for (let i = 0; i < n; i++) {
    const row: number[] = [];
    for (let j = 0; j < m; j++) {
      row.push(Number(M.get(i, j).value));
    }
    rows.push(row);
  }

  return gramSchmidt(rows);
}

/**
 * Result of SVP enumeration.
 */
interface SVPResult {
  /** Coefficients of the shortest vector in terms of basis vectors */
  coefficients: number[];
  /** Squared norm of the shortest vector */
  normSquared: number;
}

/**
 * Schnorr-Euchner enumeration for finding the shortest vector in a projected block.
 *
 * This implements the enumeration algorithm from:
 * Schnorr and Euchner, "Lattice Basis Reduction: Improved Practical Algorithms"
 *
 * The algorithm enumerates lattice points in order of increasing length in the
 * projected sublattice, using the Gram-Schmidt orthogonalization to prune branches.
 *
 * @param B - The current basis matrix
 * @param gs - The Gram-Schmidt data for B
 * @param kappa - Starting index of the block
 * @param beta - Size of the block
 * @param pruning - Optional pruning coefficients
 * @returns The shortest vector found, or null if none found
 */
function enumerateSVP(
  B: IntegerMatrix,
  gs: GramSchmidtResult,
  kappa: number,
  beta: number,
  pruning?: number[]
): SVPResult | null {
  const n = B.nrows;
  const endIdx = Math.min(kappa + beta, n);
  const blockDim = endIdx - kappa;

  if (blockDim < 1) {
    return null;
  }

  // Extract the relevant portion of Gram-Schmidt data for this block
  const mu = gs.mu;
  const Bnorms = gs.B;

  // The enumeration bound is the norm of the first projected vector
  const enumBound = Bnorms[kappa]!;

  // Apply pruning if specified
  const pruningCoeffs = pruning ?? new Array(blockDim).fill(1.0);

  // Best solution found
  let bestNormSq = enumBound;
  let bestCoeffs: number[] | null = null;

  // Schnorr-Euchner enumeration
  // We use a recursive structure represented iteratively with a stack

  // State for each level: the coefficient value and partial sums
  const coeffs: number[] = new Array(blockDim).fill(0);
  const centers: number[] = new Array(blockDim).fill(0);
  const partialNorms: number[] = new Array(blockDim + 1).fill(0);
  const delta_coeffs: number[] = new Array(blockDim).fill(0); // For zig-zag pattern

  // Start enumeration from the last coordinate (highest index in block)
  let k = blockDim - 1;
  partialNorms[blockDim] = 0;

  // Compute center for the last coordinate
  centers[k] = 0;
  coeffs[k] = 0;
  delta_coeffs[k] = 1;

  // Main enumeration loop
  let iterations = 0;
  const maxIterations = 2 ** Math.min(blockDim, 20) * 1000; // Limit iterations

  while (k < blockDim && iterations < maxIterations) {
    iterations++;

    // Compute partial norm at this level
    const idx = kappa + k;
    const diff = coeffs[k]! - centers[k]!;
    const partialNorm = partialNorms[k + 1]! + diff * diff * Bnorms[idx]!;

    // Check pruning bound
    const prunedBound = bestNormSq * (pruningCoeffs[k] ?? 1.0);

    if (partialNorm >= prunedBound) {
      // Prune this branch and try next coefficient at this level
      if (delta_coeffs[k]! > 0) {
        coeffs[k] = coeffs[k]! + delta_coeffs[k]!;
        delta_coeffs[k] = -delta_coeffs[k]! - 1;
      } else {
        coeffs[k] = coeffs[k]! + delta_coeffs[k]!;
        delta_coeffs[k] = -delta_coeffs[k]! + 1;
      }

      // Check if we've exhausted this level
      // Use a simple bound based on the remaining norm budget
      const maxCoeff = Math.ceil(Math.sqrt(prunedBound / Bnorms[idx]!)) + 1;
      if (Math.abs(coeffs[k]!) > maxCoeff) {
        // Backtrack
        k++;
        if (k < blockDim) {
          // Move to next coefficient at parent level
          if (delta_coeffs[k]! > 0) {
            coeffs[k] = coeffs[k]! + delta_coeffs[k]!;
            delta_coeffs[k] = -delta_coeffs[k]! - 1;
          } else {
            coeffs[k] = coeffs[k]! + delta_coeffs[k]!;
            delta_coeffs[k] = -delta_coeffs[k]! + 1;
          }
        }
      }
    } else if (k === 0) {
      // Reached a leaf - we have a complete vector
      // Check if it's non-zero and shorter than current best
      let isZero = true;
      for (let i = 0; i < blockDim; i++) {
        if (coeffs[i] !== 0) {
          isZero = false;
          break;
        }
      }

      if (!isZero && partialNorm < bestNormSq) {
        bestNormSq = partialNorm;
        bestCoeffs = [...coeffs];
      }

      // Move to next coefficient at this level
      if (delta_coeffs[k]! > 0) {
        coeffs[k] = coeffs[k]! + delta_coeffs[k]!;
        delta_coeffs[k] = -delta_coeffs[k]! - 1;
      } else {
        coeffs[k] = coeffs[k]! + delta_coeffs[k]!;
        delta_coeffs[k] = -delta_coeffs[k]! + 1;
      }

      // Check bound
      const maxCoeff = Math.ceil(Math.sqrt(prunedBound / Bnorms[idx]!)) + 1;
      if (Math.abs(coeffs[k]!) > maxCoeff) {
        k++;
        if (k < blockDim) {
          if (delta_coeffs[k]! > 0) {
            coeffs[k] = coeffs[k]! + delta_coeffs[k]!;
            delta_coeffs[k] = -delta_coeffs[k]! - 1;
          } else {
            coeffs[k] = coeffs[k]! + delta_coeffs[k]!;
            delta_coeffs[k] = -delta_coeffs[k]! + 1;
          }
        }
      }
    } else {
      // Go deeper
      partialNorms[k] = partialNorm;
      k--;

      // Compute center for the next level
      let center = 0;
      for (let j = k + 1; j < blockDim; j++) {
        center -= coeffs[j]! * mu[kappa + j]![kappa + k]!;
      }
      centers[k] = center;
      coeffs[k] = Math.round(center);
      delta_coeffs[k] = coeffs[k]! >= center ? 1 : -1;
    }
  }

  if (bestCoeffs === null) {
    return null;
  }

  return {
    coefficients: bestCoeffs,
    normSquared: bestNormSq,
  };
}

/**
 * Construct a lattice vector from coefficients.
 *
 * @param B - The basis matrix
 * @param coeffs - The coefficients
 * @param kappa - Starting index
 * @param beta - Block size
 * @returns The resulting vector as an array of bigints
 */
function constructVector(
  B: IntegerMatrix,
  coeffs: number[],
  kappa: number,
  beta: number
): bigint[] {
  const m = B.ncols;
  const result: bigint[] = new Array(m).fill(0n);

  for (let i = 0; i < beta && kappa + i < B.nrows; i++) {
    const c = BigInt(Math.round(coeffs[i]!));
    for (let j = 0; j < m; j++) {
      result[j] = result[j]! + c * B.get(kappa + i, j).value;
    }
  }

  return result;
}

/**
 * Insert a new vector at position kappa and re-reduce the basis.
 *
 * The algorithm:
 * 1. Insert the new vector at position kappa
 * 2. Apply LLL reduction to restore the LLL property
 * 3. Remove any zero vectors that may result from linear dependencies
 *
 * @param B - The current basis
 * @param newVector - The new vector to insert
 * @param kappa - The position to insert at
 * @param delta - The LLL delta parameter
 * @returns The new reduced basis
 */
function insertAndReduce(
  B: IntegerMatrix,
  newVector: bigint[],
  kappa: number,
  delta: number
): IntegerMatrix {
  const n = B.nrows;
  const m = B.ncols;

  // Create new matrix with the vector inserted at position kappa
  // We swap the new vector with vector at kappa rather than inserting
  const newRows: bigint[][] = [];

  for (let i = 0; i < n; i++) {
    if (i === kappa) {
      // Insert the new vector at position kappa
      newRows.push([...newVector]);
    }
    const row: bigint[] = [];
    for (let j = 0; j < m; j++) {
      row.push(B.get(i, j).value);
    }
    newRows.push(row);
  }

  // Create matrix with n+1 rows
  let augmented = IntegerMatrixFromEntries(newRows);

  // Apply LLL to restore reduction
  augmented = lllReduce(augmented, { delta, eta: 0.501 });

  // Remove zero rows and extra rows to maintain original dimension
  const finalRows: bigint[][] = [];
  for (let i = 0; i < augmented.nrows && finalRows.length < n; i++) {
    let isZero = true;
    const row: bigint[] = [];
    for (let j = 0; j < m; j++) {
      const val = augmented.get(i, j).value;
      row.push(val);
      if (val !== 0n) {
        isZero = false;
      }
    }
    if (!isZero) {
      finalRows.push(row);
    }
  }

  // If we have fewer rows than expected, something went wrong
  // Just return the LLL-reduced version
  if (finalRows.length < n) {
    // Try to recover by using the first n non-zero rows
    const result: bigint[][] = [];
    for (let i = 0; i < augmented.nrows && result.length < n; i++) {
      let isZero = true;
      const row: bigint[] = [];
      for (let j = 0; j < m; j++) {
        const val = augmented.get(i, j).value;
        row.push(val);
        if (val !== 0n) {
          isZero = false;
        }
      }
      if (!isZero) {
        result.push(row);
      }
    }

    if (result.length === n) {
      return IntegerMatrixFromEntries(result);
    }

    // Fall back to returning original LLL-reduced matrix
    return lllReduce(B, { delta, eta: 0.501 });
  }

  return IntegerMatrixFromEntries(finalRows);
}

/**
 * Check if a matrix is BKZ-reduced with given parameters.
 *
 * A basis is BKZ-reduced with block size beta if:
 * 1. It is LLL-reduced
 * 2. For each block [i, i+beta), the first projected vector is the shortest
 *    vector in that projected sublattice
 *
 * Note: This is an approximate check that verifies LLL reduction and
 * basic BKZ properties. A full check would require SVP verification
 * for each block, which is computationally expensive.
 *
 * @param basis - The basis matrix to check
 * @param blockSize - The BKZ block size (default: 10)
 * @param delta - The LLL delta parameter (default: 0.99)
 * @returns True if the basis appears to be BKZ-reduced
 */
export function isBKZReduced(
  basis: IntegerMatrix,
  blockSize: number = 10,
  delta: number = 0.99
): boolean {
  const n = basis.nrows;

  if (n <= 1) {
    return true;
  }

  // First check LLL reduction
  const gs = computeGramSchmidtFromMatrix(basis);
  const { mu, B: Bnorms } = gs;
  const eta = 0.501;

  // Check size reduction
  for (let i = 1; i < n; i++) {
    for (let j = 0; j < i; j++) {
      if (Math.abs(mu[i]![j]!) > eta + 1e-10) {
        return false;
      }
    }
  }

  // Check Lovasz condition
  for (let k = 1; k < n; k++) {
    const muKK1 = mu[k]![k - 1]!;
    const lhs = delta * Bnorms[k - 1]!;
    const rhs = Bnorms[k]! + muKK1 * muKK1 * Bnorms[k - 1]!;
    if (lhs > rhs + 1e-10) {
      return false;
    }
  }

  // For a more complete check, we would verify SVP in each block
  // This is omitted for efficiency - the above checks are necessary but not sufficient

  return true;
}

/**
 * Perform HKZ (Hermite-Korkine-Zolotarev) reduction.
 *
 * HKZ reduction is equivalent to BKZ with block_size = rank.
 * An HKZ-reduced basis has the property that:
 * 1. The basis is size-reduced
 * 2. b_1 realizes the first minimum lambda_1(L)
 * 3. The projection of b_2, ..., b_n onto b_1^perp is HKZ reduced
 *
 * @param basis - The lattice basis matrix
 * @param options - BKZ options (blockSize is ignored and set to rank)
 * @returns An HKZ-reduced basis
 *
 * @example
 * ```typescript
 * const basis = IntegerMatrixFromEntries([
 *   [1n, 0n, 3n],
 *   [0n, 2n, 1n],
 *   [0n, 2n, 7n],
 * ]);
 *
 * const hkzReduced = HKZ(basis);
 * // hkzReduced[0] is a shortest vector in the lattice
 * ```
 *
 * @see Reference: sage/modules/free_module_integer.py:HKZ
 */
export function HKZ(basis: IntegerMatrix, options?: Omit<BKZOptions, 'blockSize'>): IntegerMatrix {
  return BKZ(basis, { ...options, blockSize: basis.nrows });
}

/**
 * Compute the Hermite factor of a reduced basis.
 *
 * The Hermite factor delta is defined as:
 * delta = (||b_1|| / det(L)^{1/n})^{1/n}
 *
 * This measures the quality of lattice reduction. For LLL-reduced bases,
 * delta is approximately 1.02 - 1.04. For BKZ-reduced bases with larger
 * block sizes, delta approaches 1.
 *
 * @param basis - A reduced basis matrix
 * @returns The root Hermite factor
 */
export function computeHermiteFactor(basis: IntegerMatrix): number {
  const n = basis.nrows;

  if (n === 0) {
    return 1;
  }

  // Compute ||b_1||
  let normSq = 0;
  for (let j = 0; j < basis.ncols; j++) {
    const v = Number(basis.get(0, j).value);
    normSq += v * v;
  }
  const norm = Math.sqrt(normSq);

  // Compute det(L)^{1/n}
  // For a basis matrix B, det(L) = sqrt(det(B * B^T))
  // We use the Gram matrix
  const gramDet = computeGramDeterminant(basis);
  const volume = Math.sqrt(Math.abs(gramDet));
  const volumeRoot = volume ** (1 / n);

  if (volumeRoot === 0) {
    return Number.POSITIVE_INFINITY;
  }

  const hermiteFactor = norm / volumeRoot;
  return hermiteFactor ** (1 / n);
}

/**
 * Compute the determinant of the Gram matrix B * B^T.
 */
function computeGramDeterminant(B: IntegerMatrix): number {
  const n = B.nrows;
  const m = B.ncols;

  // Compute Gram matrix entries
  const gram: number[][] = [];
  for (let i = 0; i < n; i++) {
    gram.push([]);
    for (let j = 0; j < n; j++) {
      let dot = 0;
      for (let k = 0; k < m; k++) {
        dot += Number(B.get(i, k).value) * Number(B.get(j, k).value);
      }
      gram[i]!.push(dot);
    }
  }

  // Compute determinant using LU decomposition
  // Simple implementation for small matrices
  return luDeterminant(gram);
}

/**
 * Compute determinant using LU decomposition.
 */
function luDeterminant(A: number[][]): number {
  const n = A.length;
  if (n === 0) return 1;
  if (n === 1) return A[0]![0]!;

  // Create a copy to avoid modifying the original
  const M = A.map((row) => [...row]);

  let det = 1;

  for (let i = 0; i < n; i++) {
    // Find pivot
    let maxVal = Math.abs(M[i]![i]!);
    let maxRow = i;
    for (let k = i + 1; k < n; k++) {
      if (Math.abs(M[k]![i]!) > maxVal) {
        maxVal = Math.abs(M[k]![i]!);
        maxRow = k;
      }
    }

    if (maxVal < 1e-15) {
      return 0; // Singular matrix
    }

    // Swap rows
    if (maxRow !== i) {
      [M[i], M[maxRow]] = [M[maxRow]!, M[i]!];
      det = -det;
    }

    det *= M[i]![i]!;

    // Eliminate below
    for (let k = i + 1; k < n; k++) {
      const factor = M[k]![i]! / M[i]![i]!;
      for (let j = i; j < n; j++) {
        M[k]![j] = M[k]![j]! - factor * M[i]![j]!;
      }
    }
  }

  return det;
}

/**
 * Estimate the BKZ block size needed to find vectors of a given length.
 *
 * Uses the Gaussian heuristic to estimate the required block size.
 *
 * @param dimension - The lattice dimension n
 * @param logVolume - log(volume) of the lattice
 * @param logTargetNorm - log(target vector norm)
 * @returns Estimated block size needed
 */
export function estimateBlockSize(
  dimension: number,
  logVolume: number,
  logTargetNorm: number
): number {
  // Using the BKZ root Hermite factor approximation:
  // delta_beta ≈ ((pi * beta)^{1/beta} * beta / (2 * pi * e))^{1/(2*(beta-1))}
  //
  // For BKZ-beta reduced basis: ||b_1|| ≈ delta_beta^n * vol(L)^{1/n}
  //
  // We want to find beta such that:
  // log(target) ≈ n * log(delta_beta) + logVolume / n

  const n = dimension;
  const targetLogDelta = (logTargetNorm - logVolume / n) / n;

  // Binary search for the required block size
  let lo = 2;
  let hi = n;

  while (lo < hi) {
    const mid = Math.floor((lo + hi) / 2);
    const logDelta = computeLogRootHermite(mid);

    if (logDelta > targetLogDelta) {
      lo = mid + 1;
    } else {
      hi = mid;
    }
  }

  return Math.min(lo, n);
}

/**
 * Compute log of the root Hermite factor for BKZ-beta.
 *
 * Uses the approximation: delta_beta ≈ (beta / (2 * pi * e))^{1/(2*beta)}
 */
function computeLogRootHermite(beta: number): number {
  if (beta <= 1) return Number.POSITIVE_INFINITY;

  // Hermite constant approximation for BKZ
  // gamma_beta ≈ beta / (2 * pi * e)
  const logGamma = Math.log(beta) - Math.log(2 * Math.PI * Math.E);
  return logGamma / (2 * beta);
}
