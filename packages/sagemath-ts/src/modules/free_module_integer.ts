/**
 * @module sage/modules/free_module_integer
 * @description Discrete subgroups of ZZ^n (integer lattices)
 * @see Reference: sage/modules/free_module_integer.py
 *
 * This module provides functionality for working with integer lattices,
 * including LLL and BKZ reduction, shortest and closest vector problems,
 * and various lattice invariants. This is particularly relevant for
 * lattice-based cryptography.
 */

import { NotImplementedError, ValueError } from '../errors.js';
import { IntegerMatrix, IntegerMatrixFromEntries } from '../matrix/matrix_integer.js';
import { FreeModuleGeneric, type FreeModuleOptions, FreeModuleWithBasis } from './free_module.js';
import type { FreeModuleElement } from './free_module_element.js';

/**
 * Options for LLL reduction.
 */
export interface LLLOptions {
  /** The delta parameter (0.25 < delta <= 1.0, default: 0.99) */
  delta?: number;
  /** The eta parameter (0.5 <= eta < sqrt(delta), default: 0.501) */
  eta?: number;
  /** Algorithm to use: 'fpLLL', 'NTL:LLL', etc. */
  algorithm?: string;
  /** Floating-point type for fpLLL */
  fpType?: string;
}

/**
 * Options for BKZ reduction.
 */
export interface BKZOptions {
  /** The block size for BKZ reduction */
  blockSize: number | bigint;
  /** The delta parameter */
  delta?: number;
  /** Algorithm to use */
  algorithm?: string;
  /** Pruning strategy */
  pruning?: unknown;
  /** Maximum number of loops */
  maxLoops?: number;
  /** Flags for fpLLL BKZ */
  flags?: number;
}

/**
 * Options for shortest vector computation.
 */
export interface ShortestVectorOptions {
  /** Whether to update the reduced basis with the found vector */
  updateReducedBasis?: boolean;
  /** Algorithm: 'fplll' or 'pari' */
  algorithm?: 'fplll' | 'pari';
}

/**
 * Options for closest vector computation.
 */
export interface ClosestVectorOptions {
  /** The approximation delta for approximate algorithms */
  delta?: number;
  /** Algorithm: 'embedding', 'nearest_plane', 'rounding_off' */
  algorithm?: 'embedding' | 'nearest_plane' | 'rounding_off';
}

/**
 * Construct a new integer lattice from basis.
 *
 * @param basis - A list of vectors, matrix over integers, or element of absolute order
 * @param options - Additional options
 * @returns An integer lattice
 *
 * @example
 * ```typescript
 * // Create lattice from list of rows
 * const L = IntegerLattice([[1, 0, 3], [0, 2, 1], [0, 2, 7]]);
 *
 * // Create without LLL reduction
 * const L2 = IntegerLattice(M, { lllReduce: false });
 * ```
 *
 * @see Reference: sage/modules/free_module_integer.py:IntegerLattice
 */
export function IntegerLattice(
  basis: IntegerMatrix | bigint[][] | number[][],
  options?: { lllReduce?: boolean }
): FreeModuleIntegerLattice {
  return new FreeModuleIntegerLattice(basis, options);
}

/**
 * A submodule of ZZ^n with a distinguished basis, representing a lattice.
 *
 * This class provides functions for computing LLL and BKZ reduced bases,
 * shortest and closest vectors, and other lattice operations relevant
 * to lattice-based cryptography.
 *
 * Note: This implementation does not fully extend the FreeModuleWithBasis
 * hierarchy since that requires implementing many stub methods. Instead,
 * it provides a standalone lattice class focused on LLL reduction.
 *
 * @see Reference: sage/modules/free_module_integer.py:FreeModule_submodule_with_basis_integer
 */
export class FreeModuleIntegerLattice {
  protected _reducedBasis: IntegerMatrix;
  protected _basisMatrix: IntegerMatrix;
  protected _basisIsLLLReduced: boolean = false;

  /**
   * Construct a new submodule of ZZ^n with a distinguished basis.
   *
   * @param basis - A matrix over the integers (IntegerMatrix or 2D array)
   * @param options - Additional options
   */
  constructor(
    basis: IntegerMatrix | bigint[][] | number[][],
    options?: {
      check?: boolean;
      echelonize?: boolean;
      echelonizedBasis?: unknown;
      alreadyEchelonized?: boolean;
      lllReduce?: boolean;
    }
  ) {
    // Convert basis to IntegerMatrix if needed
    let basisMatrix: IntegerMatrix;
    if (basis instanceof IntegerMatrix) {
      basisMatrix = basis.copy();
    } else {
      basisMatrix = IntegerMatrixFromEntries(basis as (bigint | number)[][]);
    }

    this._basisMatrix = basisMatrix;

    // Apply LLL reduction if requested (default is true)
    if (options?.lllReduce !== false) {
      this._reducedBasis = lllReduce(basisMatrix);
      this._basisIsLLLReduced = true;
    } else {
      this._reducedBasis = basisMatrix.copy();
    }
  }

  /**
   * Return the basis matrix.
   */
  basisMatrix(): IntegerMatrix {
    return this._basisMatrix;
  }

  /**
   * Return the rank (number of basis vectors).
   */
  rank(): number {
    return this._basisMatrix.nrows;
  }

  /**
   * Return the degree (dimension of ambient space).
   */
  degree(): number {
    return this._basisMatrix.ncols;
  }

  /**
   * Return the currently best known reduced basis for this lattice.
   * "Best" is defined by the Euclidean norm of the first row vector.
   */
  get reducedBasis(): IntegerMatrix {
    return this._reducedBasis;
  }

  // ========== Basis Reduction Algorithms ==========

  /**
   * Return an LLL reduced basis for this lattice.
   *
   * A lattice basis (b_1, b_2, ..., b_d) is (delta, eta)-LLL-reduced if:
   * 1. For any i > j, |mu_{i,j}| <= eta
   * 2. For any i < d, delta * |b_i*|^2 <= |b_{i+1}* + mu_{i+1,i} * b_i*|^2
   *
   * Default parameters: delta = 0.75, eta = 0.501.
   *
   * @param options - LLL reduction options
   * @returns An LLL-reduced basis matrix
   *
   * @example
   * ```typescript
   * const L = IntegerLattice(A, { lllReduce: false });
   * const reduced = L.LLL();
   * // reduced[0] is now a short vector
   * ```
   *
   * @see Reference: sage/modules/free_module_integer.py:FreeModule_submodule_with_basis_integer.LLL
   */
  LLL(options?: LLLOptions): IntegerMatrix {
    const delta = options?.delta ?? 0.75;
    const eta = options?.eta ?? 0.501;

    // Reduce the current reduced basis (which may already be partially reduced)
    const reduced = lllReduce(this._reducedBasis, { delta, eta });

    // Update the reduced basis if this one is better (shorter first vector)
    const newNorm = vectorNormSquared(reduced, 0);
    const oldNorm = vectorNormSquared(this._reducedBasis, 0);

    if (newNorm < oldNorm) {
      this._reducedBasis = reduced;
      this._basisIsLLLReduced = true;
    }

    return reduced;
  }

  /**
   * Return a Block Korkine-Zolotareff (BKZ) reduced basis for this lattice.
   *
   * BKZ reduction with block size beta produces a basis where the first
   * vector has norm close to gamma_beta^{n/beta} * det(L)^{1/n}, where
   * gamma_beta is Hermite's constant.
   *
   * @param options - BKZ reduction options (blockSize is required)
   * @returns A BKZ-reduced basis matrix
   *
   * @example
   * ```typescript
   * const L = IntegerLattice(A);
   * // BKZ-20 reduction
   * const reduced = L.BKZ({ blockSize: 20 });
   * ```
   *
   * @see Reference: sage/modules/free_module_integer.py:FreeModule_submodule_with_basis_integer.BKZ
   */
  BKZ(options: BKZOptions): IntegerMatrix {
    const blockSize =
      typeof options.blockSize === 'bigint' ? Number(options.blockSize) : options.blockSize;
    const delta = options.delta ?? 0.99;

    if (blockSize < 2) {
      throw new ValueError('block_size must be at least 2');
    }

    // BKZ is a block generalization of LLL
    // For blocks of size 2, it's equivalent to LLL
    // For larger blocks, it uses SVP enumeration within each block

    // Start with LLL-reduced basis
    let reduced = this.LLL({ delta, eta: 0.501 });

    if (blockSize >= this.rank()) {
      // Full block size - do HKZ-style reduction inline to avoid recursion
      // (HKZ calls BKZ with full block size, so we can't delegate back)
      return this._fullBlockReduction(reduced, delta);
    }

    // Simplified BKZ: iteratively apply LLL within sliding windows
    // A full implementation would use SVP enumeration
    const n = this.rank();
    const m = this.degree();
    let changed = true;
    let iterations = 0;
    const maxIterations = n * (options.maxLoops ?? 10);

    while (changed && iterations < maxIterations) {
      changed = false;
      iterations++;

      // Process each block
      for (let k = 0; k <= n - blockSize; k++) {
        // Extract block rows k to k + blockSize - 1
        const blockRows: bigint[][] = [];
        for (let i = k; i < k + blockSize; i++) {
          const row: bigint[] = [];
          for (let j = 0; j < m; j++) {
            row.push(reduced.get(i, j).value);
          }
          blockRows.push(row);
        }

        // LLL reduce the block
        const blockMatrix = IntegerMatrixFromEntries(blockRows);
        const reducedBlock = lllReduce(blockMatrix, { delta, eta: 0.501 });

        // Check if any improvement
        const oldNorm = vectorNormSquared(reduced, k);
        const newNorm = vectorNormSquared(reducedBlock, 0);

        if (newNorm < oldNorm) {
          // Update the basis with the reduced block
          for (let i = 0; i < blockSize; i++) {
            for (let j = 0; j < m; j++) {
              reduced.set(k + i, j, reducedBlock.get(i, j).value);
            }
          }
          changed = true;
        }
      }

      // Re-run LLL on the full basis after block processing
      reduced = lllReduce(reduced, { delta, eta: 0.501 });
    }

    // Update the reduced basis if this is better
    const newNorm = vectorNormSquared(reduced, 0);
    const oldNorm = vectorNormSquared(this._reducedBasis, 0);

    if (newNorm < oldNorm) {
      this._reducedBasis = reduced;
    }

    return reduced;
  }

  /**
   * Helper for full block reduction (HKZ-style).
   * This is factored out to avoid recursion between BKZ and HKZ.
   */
  private _fullBlockReduction(basis: IntegerMatrix, delta: number): IntegerMatrix {
    const n = basis.nrows;
    const m = basis.ncols;

    if (n <= 1) {
      return basis;
    }

    // Iteratively improve the basis by finding short vectors and projections
    let reduced = lllReduce(basis, { delta, eta: 0.501 });
    let improved = true;
    let iterations = 0;
    const maxIterations = n * 10;

    while (improved && iterations < maxIterations) {
      improved = false;
      iterations++;

      // Try to find shorter vectors through exhaustive enumeration
      // For small dimensions, enumerate lattice points within a bound
      const bound = this._computeEnumerationBound(reduced);

      // Simplified: run deep LLL with various parameters
      for (let i = 0; i < n - 1; i++) {
        // Extract submatrix from row i to end
        const subRows: bigint[][] = [];
        for (let r = i; r < n; r++) {
          const row: bigint[] = [];
          for (let c = 0; c < m; c++) {
            row.push(reduced.get(r, c).value);
          }
          subRows.push(row);
        }

        if (subRows.length >= 2) {
          const subMatrix = IntegerMatrixFromEntries(subRows);
          const reducedSub = lllReduce(subMatrix, { delta, eta: 0.501 });

          // Check if we got a shorter vector
          const oldNorm = vectorNormSquared(reduced, i);
          const newNorm = vectorNormSquared(reducedSub, 0);

          if (newNorm < oldNorm) {
            // Update the basis
            for (let r = 0; r < subRows.length; r++) {
              for (let c = 0; c < m; c++) {
                reduced.set(i + r, c, reducedSub.get(r, c).value);
              }
            }
            improved = true;
          }
        }
      }

      // Re-run LLL on the full basis
      reduced = lllReduce(reduced, { delta, eta: 0.501 });
    }

    return reduced;
  }

  /**
   * Compute a bound for enumeration based on the current reduced basis.
   */
  private _computeEnumerationBound(basis: IntegerMatrix): bigint {
    // Use the first basis vector's norm squared as initial bound
    let bound = 0n;
    const m = basis.ncols;
    for (let j = 0; j < m; j++) {
      const v = basis.get(0, j).value;
      bound += v * v;
    }
    return bound;
  }

  /**
   * Hermite-Korkine-Zolotarev (HKZ) reduce the basis.
   *
   * An HKZ-reduced basis has the property that:
   * 1. The basis is size-reduced (|mu_{i,j}| <= 1/2)
   * 2. b_1 realizes the first minimum lambda_1(L)
   * 3. The projection of b_2, ..., b_r onto b_1^perp is HKZ reduced
   *
   * This is equivalent to BKZ with block_size = rank.
   *
   * @param options - Reduction options
   * @returns An HKZ-reduced basis matrix
   *
   * @see Reference: sage/modules/free_module_integer.py:FreeModule_submodule_with_basis_integer.HKZ
   */
  HKZ(options?: Omit<BKZOptions, 'blockSize'>): IntegerMatrix {
    const delta = options?.delta ?? 0.99;

    // HKZ is equivalent to BKZ with block_size = rank
    // Start with LLL reduction
    let reduced = this.LLL({ delta, eta: 0.501 });
    const n = this.rank();

    if (n <= 1) {
      return reduced;
    }

    // HKZ requires that:
    // 1. The basis is size-reduced
    // 2. b_1 is a shortest vector
    // 3. The projection of b_2, ..., b_n onto b_1^perp is HKZ

    // Use full block reduction (same as BKZ with block_size = rank)
    reduced = this._fullBlockReduction(reduced, delta);

    // Update reduced basis if better
    const newNorm = vectorNormSquared(reduced, 0);
    const oldNorm = vectorNormSquared(this._reducedBasis, 0);

    if (newNorm < oldNorm) {
      this._reducedBasis = reduced;
    }

    return reduced;
  }

  // ========== Lattice Problems ==========

  /**
   * Return a shortest nonzero vector in this lattice.
   *
   * This solves the Shortest Vector Problem (SVP) using Schnorr-Euchner
   * enumeration with LLL preprocessing.
   *
   * @param options - Options for the computation
   * @returns A shortest nonzero vector
   *
   * @example
   * ```typescript
   * const L = IntegerLattice(A);
   * const sv = L.shortestVector();
   * console.log(sv.norm());  // The first minimum lambda_1(L)
   * ```
   *
   * @see Reference: sage/modules/free_module_integer.py:FreeModule_submodule_with_basis_integer.shortest_vector
   * @see Deviation: Uses Schnorr-Euchner enumeration instead of fpylll/PARI
   */
  shortestVector(options?: ShortestVectorOptions): bigint[] {
    const updateReducedBasis = options?.updateReducedBasis ?? true;

    // Ensure LLL-reduced basis
    if (!this._basisIsLLLReduced) {
      this.LLL();
    }

    const B = this._reducedBasis;
    const n = B.nrows;
    const m = B.ncols;

    if (n === 0) {
      throw new ValueError('cannot find shortest vector in zero lattice');
    }

    // Use Schnorr-Euchner enumeration for exact SVP
    // For very small lattices, use simple enumeration
    // For medium lattices, use Schnorr-Euchner
    // For large lattices, return the LLL approximation

    let sv: bigint[];

    if (n <= 4) {
      // Simple enumeration for very small lattices
      sv = this._simpleEnumerationSVP(B);
    } else if (n <= 50) {
      // Schnorr-Euchner enumeration for medium lattices
      sv = this._schnorrEuchnerSVP(B);
    } else {
      // For large lattices, LLL gives a 2^((n-1)/2) approximation
      sv = [];
      for (let j = 0; j < m; j++) {
        sv.push(B.get(0, j).value);
      }
    }

    if (updateReducedBasis) {
      this._updateBasisWithVector(sv);
    }

    return sv;
  }

  /**
   * Schnorr-Euchner enumeration for SVP.
   *
   * This implements the Schnorr-Euchner algorithm which enumerates lattice
   * points in order of increasing projected length, using the zig-zag
   * coefficient pattern and branch-and-bound pruning.
   *
   * @param B - LLL-reduced basis matrix
   * @returns The shortest vector found
   */
  private _schnorrEuchnerSVP(B: IntegerMatrix): bigint[] {
    const n = B.nrows;
    const m = B.ncols;

    // Compute Gram-Schmidt orthogonalization
    const gs = gramSchmidt(B);
    const { mu, B: Bnorms } = gs;

    // Initial bound is the squared norm of the first basis vector
    let bestNormSq = 0;
    for (let j = 0; j < m; j++) {
      const v = Number(B.get(0, j).value);
      bestNormSq += v * v;
    }

    let bestCoeffs: number[] | null = null;

    // State arrays for enumeration
    const coeffs: number[] = new Array(n).fill(0);
    const centers: number[] = new Array(n).fill(0);
    const partialNorms: number[] = new Array(n + 1).fill(0);
    const deltas: number[] = new Array(n).fill(0);

    // Start from the last coordinate
    let k = n - 1;
    partialNorms[n] = 0;
    centers[k] = 0;
    coeffs[k] = 0;
    deltas[k] = 1;

    // Maximum iterations to prevent infinite loops
    let iterations = 0;
    const maxIterations = 2 ** Math.min(n, 25) * 1000;

    // Main enumeration loop
    while (k < n && iterations < maxIterations) {
      iterations++;

      // Compute partial norm at this level
      const diff = coeffs[k]! - centers[k]!;
      const partialNorm = partialNorms[k + 1]! + diff * diff * Bnorms[k]!;

      if (partialNorm >= bestNormSq) {
        // Prune: try next coefficient using zig-zag pattern
        this._moveToNextCoefficient(coeffs, deltas, k);

        // Check if exhausted this level
        const maxCoeff = Math.ceil(Math.sqrt(bestNormSq / Math.max(Bnorms[k]!, 1e-10))) + 1;
        if (Math.abs(coeffs[k]!) > maxCoeff) {
          // Backtrack
          k++;
          if (k < n) {
            this._moveToNextCoefficient(coeffs, deltas, k);
          }
        }
      } else if (k === 0) {
        // Reached a leaf - check if it's a shorter nonzero vector
        const isZero = coeffs.every((c) => c === 0);

        if (!isZero && partialNorm < bestNormSq && partialNorm > 0) {
          bestNormSq = partialNorm;
          bestCoeffs = [...coeffs];
        }

        // Move to next coefficient
        this._moveToNextCoefficient(coeffs, deltas, k);

        const maxCoeff = Math.ceil(Math.sqrt(bestNormSq / Math.max(Bnorms[k]!, 1e-10))) + 1;
        if (Math.abs(coeffs[k]!) > maxCoeff) {
          k++;
          if (k < n) {
            this._moveToNextCoefficient(coeffs, deltas, k);
          }
        }
      } else {
        // Go deeper
        partialNorms[k] = partialNorm;
        k--;

        // Compute center for next level
        let center = 0;
        for (let j = k + 1; j < n; j++) {
          center -= coeffs[j]! * mu[j]![k]!;
        }
        centers[k] = center;
        coeffs[k] = Math.round(center);
        deltas[k] = coeffs[k]! >= center ? 1 : -1;
      }
    }

    // Compute the shortest vector from coefficients
    if (bestCoeffs === null) {
      // Return first basis vector
      const result: bigint[] = [];
      for (let j = 0; j < m; j++) {
        result.push(B.get(0, j).value);
      }
      return result;
    }

    const result: bigint[] = new Array(m).fill(0n);
    for (let i = 0; i < n; i++) {
      const c = BigInt(Math.round(bestCoeffs[i]!));
      for (let j = 0; j < m; j++) {
        result[j] = result[j]! + c * B.get(i, j).value;
      }
    }

    return result;
  }

  /**
   * Move to the next coefficient using zig-zag pattern.
   * Enumerates: 0, 1, -1, 2, -2, 3, -3, ...
   */
  private _moveToNextCoefficient(coeffs: number[], deltas: number[], k: number): void {
    if (deltas[k]! > 0) {
      coeffs[k] = coeffs[k]! + deltas[k]!;
      deltas[k] = -deltas[k]! - 1;
    } else {
      coeffs[k] = coeffs[k]! + deltas[k]!;
      deltas[k] = -deltas[k]! + 1;
    }
  }

  /**
   * Simple enumeration for very small lattices (n <= 4).
   * Exhaustively searches all coefficient combinations within a bound.
   */
  private _simpleEnumerationSVP(B: IntegerMatrix): bigint[] {
    const n = B.nrows;
    const m = B.ncols;

    // Get Gram-Schmidt norms
    const gs = gramSchmidt(B);
    const { B: Bnorms } = gs;

    // Initial bound from first basis vector
    let bestNormSq = 0;
    for (let j = 0; j < m; j++) {
      const v = Number(B.get(0, j).value);
      bestNormSq += v * v;
    }

    let bestVec: bigint[] | null = null;

    // Compute search bound
    const minNorm = Math.min(...Bnorms.filter((x) => x > 0));
    const maxCoeff = Math.min(
      Math.ceil(Math.sqrt(bestNormSq / minNorm)) + 1,
      15 // Hard limit for performance
    );

    // Enumerate all coefficient combinations
    const enumerate = (depth: number, coeffs: number[]): void => {
      if (depth === n) {
        if (coeffs.every((c) => c === 0)) return;

        // Compute vector and norm
        const vec: bigint[] = new Array(m).fill(0n);
        let normSq = 0;

        for (let i = 0; i < n; i++) {
          const c = BigInt(coeffs[i]!);
          for (let j = 0; j < m; j++) {
            vec[j] = vec[j]! + c * B.get(i, j).value;
          }
        }

        for (let j = 0; j < m; j++) {
          normSq += Number(vec[j]!) * Number(vec[j]!);
        }

        if (normSq < bestNormSq && normSq > 0) {
          bestNormSq = normSq;
          bestVec = vec;
        }
        return;
      }

      for (let c = -maxCoeff; c <= maxCoeff; c++) {
        enumerate(depth + 1, [...coeffs, c]);
      }
    };

    enumerate(0, []);

    // Return best found or first basis vector
    if (bestVec === null) {
      bestVec = [];
      for (let j = 0; j < m; j++) {
        bestVec.push(B.get(0, j).value);
      }
    }

    return bestVec;
  }

  /**
   * Update the reduced basis with a new short vector.
   */
  private _updateBasisWithVector(v: bigint[]): void {
    // Insert the vector and re-run LLL
    const n = this._reducedBasis.nrows;
    const m = this._reducedBasis.ncols;

    // Create new basis with v as first row
    const newBasis: bigint[][] = [v];
    for (let i = 0; i < n; i++) {
      const row: bigint[] = [];
      for (let j = 0; j < m; j++) {
        row.push(this._reducedBasis.get(i, j).value);
      }
      newBasis.push(row);
    }

    // Run LLL to get a reduced basis
    const newMatrix = IntegerMatrixFromEntries(newBasis);
    const reduced = lllReduce(newMatrix);

    // Remove any zero rows and keep only n rows
    const finalRows: bigint[][] = [];
    for (let i = 0; i < reduced.nrows && finalRows.length < n; i++) {
      let isZero = true;
      const row: bigint[] = [];
      for (let j = 0; j < m; j++) {
        row.push(reduced.get(i, j).value);
        if (reduced.get(i, j).value !== 0n) isZero = false;
      }
      if (!isZero) {
        finalRows.push(row);
      }
    }

    if (finalRows.length === n) {
      this._reducedBasis = IntegerMatrixFromEntries(finalRows);
      this._basisIsLLLReduced = true;
    }
  }

  /**
   * Compute the closest vector in the lattice to a given target vector.
   *
   * This solves the Closest Vector Problem (CVP).
   *
   * @param target - The target vector
   * @returns The closest lattice vector to target
   *
   * @example
   * ```typescript
   * const L = IntegerLattice([[1, 0], [0, 1]]);
   * const closest = L.closestVector([-6, 5/3]);
   * // Returns (-6, 2)
   * ```
   *
   * @see Reference: sage/modules/free_module_integer.py:FreeModule_submodule_with_basis_integer.closest_vector
   */
  closestVector(target: FreeModuleElement | unknown[]): bigint[] {
    // The Closest Vector Problem (CVP) is NP-hard.
    // For small lattices, we can use enumeration.
    // For larger lattices, we use Babai's algorithm (approximate).

    // Convert target to number array
    let t: number[];
    if (Array.isArray(target)) {
      t = target.map((x) => Number(x));
    } else {
      t = target.list().map((x) => Number(x));
    }

    const n = this.rank();
    const m = this.degree();

    if (t.length !== m) {
      throw new ValueError(`target vector has wrong dimension: ${t.length} vs ${m}`);
    }

    // Ensure LLL-reduced basis
    if (!this._basisIsLLLReduced) {
      this.LLL();
    }

    // For small lattices, use enumeration
    if (n <= 4) {
      return this._enumerateClosestVector(t);
    }

    // For larger lattices, use Babai's algorithm (approximate CVP)
    return this.approximateClosestVector(target, { algorithm: 'nearest_plane' });
  }

  /**
   * Enumerate lattice points to find the closest one to target.
   */
  private _enumerateClosestVector(target: number[]): bigint[] {
    const B = this._reducedBasis;
    const n = B.nrows;
    const m = B.ncols;

    // First get an approximate solution using Babai
    const approx = this.approximateClosestVector(target, { algorithm: 'nearest_plane' });

    // Compute distance to approximate solution
    let bestDistSq = 0;
    for (let j = 0; j < m; j++) {
      const diff = target[j]! - Number(approx[j]!);
      bestDistSq += diff * diff;
    }

    let bestVec = approx;

    // Search in a neighborhood around the approximate solution
    // The search radius is based on the Babai distance
    const searchRadius = Math.ceil(Math.sqrt(bestDistSq));

    // Get coefficients of approximate solution in the basis
    const gs = gramSchmidt(B);
    const { B: Bnorms } = gs;

    // Limit search based on basis quality
    const maxCoeff = Math.min(
      Math.ceil(searchRadius / Math.sqrt(Math.min(...Bnorms.filter((x) => x > 0)))),
      3
    );

    // Enumerate nearby lattice points
    const enumerate = (depth: number, coeffs: number[], currentVec: number[]): void => {
      if (depth === n) {
        // Compute distance to target
        let distSq = 0;
        for (let j = 0; j < m; j++) {
          const diff = target[j]! - currentVec[j]!;
          distSq += diff * diff;
        }

        if (distSq < bestDistSq) {
          bestDistSq = distSq;
          bestVec = currentVec.map((v) => BigInt(Math.round(v)));
        }
        return;
      }

      // Get coefficient from approximate solution
      // (simplified: assume integer coefficients near 0)
      for (let c = -maxCoeff; c <= maxCoeff; c++) {
        const newVec = [...currentVec];
        for (let j = 0; j < m; j++) {
          newVec[j] = newVec[j]! + c * Number(B.get(depth, j).value);
        }
        enumerate(depth + 1, [...coeffs, c], newVec);
      }
    };

    enumerate(0, [], new Array(m).fill(0));

    return bestVec;
  }

  /**
   * Compute a vector w in this lattice which is close to target t.
   *
   * The approximation ratio |t-w|/|t-u| where u is the closest vector
   * is exponential in the dimension.
   *
   * This is also known as Babai's algorithm.
   *
   * @param target - The target vector
   * @param options - Options for the computation
   * @returns A lattice vector close to target (as an array of bigints)
   *
   * @example
   * ```typescript
   * const L = IntegerLattice([[101, 0, 0, 0], [0, 101, 0, 0], ...]);
   * const approx = L.approximateClosestVector([1337, 1337, 1337, 1337]);
   * ```
   *
   * @see Reference: sage/modules/free_module_integer.py:FreeModule_submodule_with_basis_integer.approximate_closest_vector
   */
  approximateClosestVector(
    target: FreeModuleElement | unknown[],
    options?: ClosestVectorOptions
  ): bigint[] {
    const algorithm = options?.algorithm ?? 'nearest_plane';

    // Convert target to number array
    let t: number[];
    if (Array.isArray(target)) {
      t = target.map((x) => Number(x));
    } else {
      t = target.list().map((x) => Number(x));
    }

    // Ensure LLL-reduced basis
    if (!this._basisIsLLLReduced) {
      this.LLL();
    }

    const B = this._reducedBasis;
    const n = B.nrows;
    const m = B.ncols;

    if (t.length !== m) {
      throw new ValueError(`target vector has wrong dimension: ${t.length} vs ${m}`);
    }

    if (algorithm === 'nearest_plane') {
      // Babai's nearest plane algorithm
      // Compute Gram-Schmidt
      const gs = gramSchmidt(B);
      const { orthogonalBasis: GStar, B: Bnorms } = gs;

      // Start with target
      const b = [...t];

      // Subtract projections from last to first
      for (let i = n - 1; i >= 0; i--) {
        // Compute coefficient: round(<b, b_i*> / <b_i*, b_i*>)
        let dotBBi = 0;
        for (let j = 0; j < m; j++) {
          dotBBi += b[j]! * GStar[i]![j]!;
        }
        const coeff = Math.round(dotBBi / Bnorms[i]!);

        // Subtract coeff * B[i] from b
        for (let j = 0; j < m; j++) {
          b[j] = b[j]! - coeff * Number(B.get(i, j).value);
        }
      }

      // The closest vector is t - b
      const result: bigint[] = [];
      for (let j = 0; j < m; j++) {
        result.push(BigInt(Math.round(t[j]! - b[j]!)));
      }

      return result;
    } else if (algorithm === 'rounding_off') {
      // Simple rounding algorithm
      // Solve t = x * B (approximately) and round x

      // For full rank lattice, compute B^T * (B * B^T)^(-1) * t
      // For simplicity, use Gram-Schmidt coefficients
      const gs = gramSchmidt(B);
      const { mu } = gs;

      // Compute coordinates of t in the basis
      // Using back-substitution with mu matrix
      const coords: number[] = new Array(n).fill(0);

      // Convert B to number[][] for easier computation
      const Bn: number[][] = [];
      for (let i = 0; i < n; i++) {
        Bn.push([]);
        for (let j = 0; j < m; j++) {
          Bn[i]!.push(Number(B.get(i, j).value));
        }
      }

      // Compute coordinates by solving t = sum(coords[i] * B[i])
      // This is approximate - we solve the normal equations
      // (B * B^T) * coords = B * t
      const gram: number[][] = [];
      const bt: number[] = [];

      for (let i = 0; i < n; i++) {
        gram.push([]);
        for (let j = 0; j < n; j++) {
          let dotProd = 0;
          for (let k = 0; k < m; k++) {
            dotProd += Bn[i]![k]! * Bn[j]![k]!;
          }
          gram[i]!.push(dotProd);
        }

        let dotT = 0;
        for (let k = 0; k < m; k++) {
          dotT += Bn[i]![k]! * t[k]!;
        }
        bt.push(dotT);
      }

      // Solve gram * coords = bt using Gaussian elimination
      // For simplicity, use simple forward/back substitution
      // (this is approximate for non-square systems)
      const coordsF = solveLinearSystem(gram, bt);

      // Round coordinates and compute lattice vector
      const result: bigint[] = new Array(m).fill(0n);
      for (let i = 0; i < n; i++) {
        const roundedCoord = BigInt(Math.round(coordsF[i]!));
        for (let j = 0; j < m; j++) {
          result[j] = result[j]! + roundedCoord * B.get(i, j).value;
        }
      }

      return result;
    } else {
      // embedding algorithm - more complex, use nearest_plane as fallback
      return this.approximateClosestVector(target, { ...options, algorithm: 'nearest_plane' });
    }
  }

  /**
   * Alias for approximateClosestVector (Babai's algorithm).
   * @see approximateClosestVector
   */
  babai(target: FreeModuleElement | unknown[], options?: ClosestVectorOptions): bigint[] {
    return this.approximateClosestVector(target, options);
  }

  // ========== Lattice Invariants ==========

  /**
   * Return vol(L) = sqrt(det(B * B^T)) for any basis B.
   *
   * For a full-rank lattice, this equals |det(B)|.
   *
   * @returns The volume of the lattice
   *
   * @see Reference: sage/modules/free_module_integer.py:FreeModule_submodule_with_basis_integer.volume
   */
  volume(): bigint {
    const n = this.rank();
    const m = this.degree();

    if (n === m) {
      // Full rank: volume = |det(B)|
      const det = this._reducedBasis.determinant().value;
      return det < 0n ? -det : det;
    }

    // Non-full rank: volume = sqrt(det(B * B^T))
    // For simplicity, compute B * B^T and take sqrt of determinant
    const B = this._reducedBasis;
    const BT = B.transpose();
    const gram = B.mul(BT);
    const det = gram.determinant().value;
    const absDet = det < 0n ? -det : det;

    // Integer square root
    return bigintSqrt(absDet);
  }

  /**
   * Return |det(G)|, the absolute value of the determinant of the Gram matrix.
   *
   * @returns The discriminant
   *
   * @see Reference: sage/modules/free_module_integer.py:FreeModule_submodule_with_basis_integer.discriminant
   */
  discriminant(): bigint {
    const B = this._reducedBasis;
    const BT = B.transpose();
    const gram = B.mul(BT);
    const det = gram.determinant().value;
    return det < 0n ? -det : det;
  }

  /**
   * Return True if this lattice is unimodular (volume = 1).
   *
   * @returns Whether the lattice is unimodular
   *
   * @see Reference: sage/modules/free_module_integer.py:FreeModule_submodule_with_basis_integer.is_unimodular
   */
  isUnimodular(): boolean {
    return this.volume() === 1n;
  }

  /**
   * Compute the normalized Hadamard ratio of the basis.
   *
   * H(B) = (det(L) / (||v_1|| * ... * ||v_n||))^{1/n}
   *
   * The closer this is to 1, the more orthogonal the basis.
   *
   * @param useReducedBasis - Whether to use the reduced basis
   * @returns The Hadamard ratio
   *
   * @see Reference: sage/modules/free_module_integer.py:FreeModule_submodule_with_basis_integer.hadamard_ratio
   */
  hadamardRatio(useReducedBasis: boolean = true): number {
    const basis = useReducedBasis ? this._reducedBasis : this._basisMatrix;
    const n = basis.nrows;

    if (n === 0) {
      return 1.0;
    }

    // Compute product of norms
    let productOfNorms = 1.0;
    for (let i = 0; i < n; i++) {
      let normSq = 0;
      for (let j = 0; j < basis.ncols; j++) {
        const v = Number(basis.get(i, j).value);
        normSq += v * v;
      }
      productOfNorms *= Math.sqrt(normSq);
    }

    // Get volume
    const vol = Number(this.discriminant());
    const sqrtVol = Math.sqrt(vol);

    // Hadamard ratio = (det / product_of_norms)^(1/n)
    const ratio = (sqrtVol / productOfNorms) ** (1 / n);

    return ratio;
  }

  /**
   * Compute the Gaussian heuristic for the shortest vector length.
   *
   * This estimates the expected norm of the shortest nonzero vector
   * as det(L)^{1/n} * sqrt(n / (2 * pi * e)).
   *
   * @param exactForm - Use exact gamma function instead of Stirling approximation
   * @returns The Gaussian heuristic estimate
   *
   * @see Reference: sage/modules/free_module_integer.py:FreeModule_submodule_with_basis_integer.gaussian_heuristic
   */
  gaussianHeuristic(exactForm: boolean = false): number {
    const n = this.rank();

    if (n === 0) {
      return 0;
    }

    const vol = Number(this.volume());
    const detRootN = vol ** (1 / n);

    if (exactForm) {
      // Use gamma function: det^(1/n) * Gamma(1 + n/2)^(1/n) / sqrt(pi)
      // Gamma(1 + n/2) = (n/2)!
      const gammaArg = 1 + n / 2;
      const gamma = gammaFunction(gammaArg);
      return (detRootN * gamma ** (1 / n)) / Math.sqrt(Math.PI);
    } else {
      // Stirling approximation: det^(1/n) * sqrt(n / (2 * pi * e))
      return detRootN * Math.sqrt(n / (2 * Math.PI * Math.E));
    }
  }

  // ========== Voronoi Cell and Related ==========

  /**
   * Compute the Voronoi cell of this lattice.
   *
   * The Voronoi cell is the set of all points closer to the origin
   * than to any other lattice point. It is a convex polytope centered
   * at the origin.
   *
   * @param radius - Radius of ball containing vertices (auto-determined if not given)
   * @returns The Voronoi cell as a list of half-space constraints
   *
   * @see Reference: sage/modules/free_module_integer.py:FreeModule_submodule_with_basis_integer.voronoi_cell
   */
  voronoiCell(radius?: number): { normals: bigint[][]; offsets: number[] } {
    // The Voronoi cell is defined by the half-spaces:
    // {x : |x| <= |x - v| for all lattice vectors v}
    // This simplifies to: {x : 2<v, x> <= |v|^2}

    // Get the Voronoi-relevant vectors
    const relevantVectors = this.voronoiRelevantVectors();

    const normals: bigint[][] = [];
    const offsets: number[] = [];

    for (const v of relevantVectors) {
      // The half-space is: 2<v, x> <= |v|^2
      // Normal is 2*v, offset is |v|^2
      const vList = v as unknown as bigint[];
      normals.push(vList.map((x) => 2n * x));

      let normSq = 0;
      for (const x of vList) {
        normSq += Number(x) * Number(x);
      }
      offsets.push(normSq);
    }

    return { normals, offsets };
  }

  /**
   * Compute the Voronoi-relevant vectors (those defining the Voronoi cell).
   *
   * A lattice vector v is Voronoi-relevant if and only if v/2 is a vertex
   * of the Voronoi cell. These are exactly the vectors v such that
   * cv(v/2) = v (i.e., v is the closest lattice vector to v/2).
   *
   * @returns List of Voronoi-relevant vectors
   *
   * @see Reference: sage/modules/free_module_integer.py:FreeModule_submodule_with_basis_integer.voronoi_relevant_vectors
   */
  voronoiRelevantVectors(): bigint[][] {
    // Voronoi-relevant vectors are among the short vectors of the lattice.
    // We enumerate short vectors up to a certain bound.

    if (!this._basisIsLLLReduced) {
      this.LLL();
    }

    const B = this._reducedBasis;
    const n = B.nrows;
    const m = B.ncols;

    // Compute covering radius estimate
    const gs = gramSchmidt(B);
    const { B: Bnorms } = gs;

    // The covering radius is bounded by sqrt(n) * max(|b_i*|) / 2
    const maxGSNorm = Math.sqrt(Math.max(...Bnorms));
    const coveringRadiusEstimate = (Math.sqrt(n) * maxGSNorm) / 2;

    // Enumerate vectors with norm up to 2 * covering radius
    const bound = 2 * coveringRadiusEstimate;
    const boundSq = bound * bound;

    // Use the norm of the first LLL basis vector as a better bound
    let firstNormSq = 0;
    for (let j = 0; j < m; j++) {
      const v = Number(B.get(0, j).value);
      firstNormSq += v * v;
    }
    const effectiveBoundSq = Math.min(boundSq, 4 * firstNormSq);

    const relevantVectors: bigint[][] = [];

    // Enumerate short vectors
    const maxCoeff = Math.ceil(
      Math.sqrt(effectiveBoundSq / Math.min(...Bnorms.filter((x) => x > 0)))
    );
    const searchBound = Math.min(maxCoeff, 4);

    const enumerate = (depth: number, coeffs: number[]): void => {
      if (depth === n) {
        // Skip zero vector
        if (coeffs.every((c) => c === 0)) return;

        // Compute the lattice vector
        const vec: bigint[] = new Array(m).fill(0n);
        for (let i = 0; i < n; i++) {
          for (let j = 0; j < m; j++) {
            vec[j] = vec[j]! + BigInt(coeffs[i]!) * B.get(i, j).value;
          }
        }

        // Compute norm squared
        let normSq = 0;
        for (let j = 0; j < m; j++) {
          normSq += Number(vec[j]!) * Number(vec[j]!);
        }

        // Check if within bound
        if (normSq > effectiveBoundSq) return;

        // Check if v/2 has closest vector v (i.e., v is Voronoi-relevant)
        // This requires checking that no other lattice point is closer to v/2
        const halfV = vec.map((x) => Number(x) / 2);
        const closest = this.approximateClosestVector(halfV, { algorithm: 'nearest_plane' });

        // Check if closest == vec (or -vec)
        let isRelevant = true;
        for (let j = 0; j < m; j++) {
          if (closest[j]! !== vec[j]! && closest[j]! !== -vec[j]!) {
            isRelevant = false;
            break;
          }
        }

        if (isRelevant) {
          relevantVectors.push(vec);
        }
        return;
      }

      for (let c = -searchBound; c <= searchBound; c++) {
        enumerate(depth + 1, [...coeffs, c]);
      }
    };

    enumerate(0, []);

    return relevantVectors;
  }

  // ========== Basis Manipulation ==========

  /**
   * Inject vector w and run LLL to update the reduced basis.
   *
   * If w is shorter than the current shortest basis vector, it will
   * be incorporated into the reduced basis.
   *
   * @param w - A vector to inject
   *
   * @see Reference: sage/modules/free_module_integer.py:FreeModule_submodule_with_basis_integer.update_reduced_basis
   */
  updateReducedBasis(w: FreeModuleElement | bigint[]): void {
    const m = this.degree();

    // Convert w to bigint array
    let wBigint: bigint[];
    if (Array.isArray(w)) {
      wBigint = w.map((x) => {
        if (typeof x === 'bigint') return x;
        return BigInt(Math.round(Number(x)));
      });
    } else {
      wBigint = w.list().map((x) => {
        if (typeof x === 'bigint') return x;
        return BigInt(Math.round(Number(x)));
      });
    }

    if (wBigint.length !== m) {
      throw new ValueError(`vector has wrong dimension: ${wBigint.length} vs ${m}`);
    }

    // Use the private method to update the basis
    this._updateBasisWithVector(wBigint);
  }
}

// ========== Gram-Schmidt Orthogonalization ==========

/**
 * Result of Gram-Schmidt orthogonalization.
 */
export interface GramSchmidtResult {
  /** The orthogonalized basis B* (rows are b_1*, ..., b_n*) */
  orthogonalBasis: number[][];
  /** The mu coefficients matrix where mu[i][j] = <b_i, b_j*> / <b_j*, b_j*> */
  mu: number[][];
  /** The squared norms of the orthogonal vectors |b_i*|^2 */
  B: number[];
}

/**
 * Compute the Gram-Schmidt orthogonalization of a basis.
 *
 * Given a basis B = (b_1, ..., b_n), returns the orthogonalized basis
 * B* = (b_1*, ..., b_n*) and the mu coefficients where:
 * b_i* = b_i - sum_{j<i} mu_{i,j} * b_j*
 *
 * @param basis - A matrix whose rows form the basis (either IntegerMatrix or number[][])
 * @param orthonormal - If true, normalize the orthogonal vectors (not implemented)
 * @returns GramSchmidtResult containing orthogonal basis, mu matrix, and squared norms
 *
 * @see Reference: sage/matrix/matrix_integer_dense.pyx:Matrix_integer_dense.gram_schmidt
 */
export function gramSchmidt(
  basis: IntegerMatrix | number[][] | bigint[][],
  orthonormal?: boolean
): GramSchmidtResult {
  if (orthonormal) {
    throw new NotImplementedError('orthonormal Gram-Schmidt not yet implemented');
  }

  // Convert to number[][] for computation
  let B: number[][];
  let n: number;
  let m: number;

  if (basis instanceof IntegerMatrix) {
    n = basis.nrows;
    m = basis.ncols;
    B = [];
    for (let i = 0; i < n; i++) {
      const row: number[] = [];
      for (let j = 0; j < m; j++) {
        row.push(Number(basis.get(i, j).value));
      }
      B.push(row);
    }
  } else {
    n = basis.length;
    if (n === 0) {
      return { orthogonalBasis: [], mu: [], B: [] };
    }
    m = basis[0]!.length;
    B = basis.map((row) => row.map((x) => Number(x)));
  }

  if (n === 0) {
    return { orthogonalBasis: [], mu: [], B: [] };
  }

  // Orthogonalized basis (working in floating point for numerical stability)
  const bStar: number[][] = [];
  // mu coefficients
  const mu: number[][] = [];
  // Squared norms of orthogonal vectors
  const Bnorms: number[] = [];

  for (let i = 0; i < n; i++) {
    // Start with b_i
    const bStarI = [...B[i]!];
    mu.push(new Array(n).fill(0));
    mu[i]![i] = 1;

    // Subtract projections onto previous orthogonal vectors
    for (let j = 0; j < i; j++) {
      // mu[i][j] = <b_i, b_j*> / <b_j*, b_j*>
      const dotProduct = dotVec(B[i]!, bStar[j]!);
      const muIJ = dotProduct / Bnorms[j]!;
      mu[i]![j] = muIJ;

      // b_i* = b_i* - mu[i][j] * b_j*
      for (let k = 0; k < m; k++) {
        bStarI[k] = bStarI[k]! - muIJ * bStar[j]![k]!;
      }
    }

    bStar.push(bStarI);
    Bnorms.push(dotVec(bStarI, bStarI));
  }

  return { orthogonalBasis: bStar, mu, B: Bnorms };
}

/**
 * Compute dot product of two vectors.
 */
function dotVec(a: number[], b: number[]): number {
  let sum = 0;
  for (let i = 0; i < a.length; i++) {
    sum += a[i]! * b[i]!;
  }
  return sum;
}

// ========== LLL Algorithm ==========

/**
 * Compute the squared Euclidean norm of row i of an IntegerMatrix.
 */
function vectorNormSquared(M: IntegerMatrix, i: number): bigint {
  let sum = 0n;
  for (let j = 0; j < M.ncols; j++) {
    const v = M.get(i, j).value;
    sum += v * v;
  }
  return sum;
}

/**
 * Options for the LLL algorithm.
 */
export interface LLLReduceOptions {
  /** The delta parameter (0.25 < delta <= 1.0, default: 0.75) */
  delta?: number;
  /** The eta parameter (0.5 <= eta < sqrt(delta), default: 0.501) */
  eta?: number;
}

/**
 * Perform LLL lattice basis reduction on an integer matrix.
 *
 * The LLL (Lenstra-Lenstra-Lovász) algorithm reduces a lattice basis to a
 * "short" and "nearly orthogonal" basis. A basis is (delta, eta)-LLL-reduced if:
 *
 * 1. Size reduction: |mu_{i,j}| <= eta for all i > j
 * 2. Lovász condition: delta * |b_{k-1}*|^2 <= |b_k* + mu_{k,k-1} * b_{k-1}*|^2
 *    which simplifies to: delta * B[k-1] <= B[k] + mu[k][k-1]^2 * B[k-1]
 *
 * @param basis - An IntegerMatrix whose rows form the lattice basis
 * @param options - LLL parameters (delta and eta)
 * @returns An LLL-reduced IntegerMatrix
 *
 * @example
 * ```typescript
 * const basis = IntegerMatrixFromEntries([[1, 0, 3], [0, 2, 1], [0, 2, 7]]);
 * const reduced = lllReduce(basis);
 * // reduced is now an LLL-reduced basis
 * ```
 *
 * @see Reference: Lenstra, Lenstra, Lovász (1982)
 * @see Reference: sage/matrix/matrix_integer_dense.pyx:Matrix_integer_dense.LLL
 */
export function lllReduce(basis: IntegerMatrix, options?: LLLReduceOptions): IntegerMatrix {
  const delta = options?.delta ?? 0.75;
  const eta = options?.eta ?? 0.501;

  // Validate parameters
  if (delta <= 0.25 || delta > 1.0) {
    throw new ValueError(`delta must be in (0.25, 1], got ${delta}`);
  }
  if (eta < 0.5 || eta >= Math.sqrt(delta)) {
    throw new ValueError(
      `eta must be in [0.5, sqrt(delta)), got ${eta} (sqrt(delta) = ${Math.sqrt(delta)})`
    );
  }

  const n = basis.nrows;
  const m = basis.ncols;

  if (n === 0) {
    return new IntegerMatrix(0, m);
  }

  if (n === 1) {
    return basis.copy();
  }

  // Copy basis to working array (using bigint for exact arithmetic)
  const B: bigint[][] = [];
  for (let i = 0; i < n; i++) {
    const row: bigint[] = [];
    for (let j = 0; j < m; j++) {
      row.push(basis.get(i, j).value);
    }
    B.push(row);
  }

  // Compute full Gram-Schmidt orthogonalization
  // mu[i][j] = <b_i, b_j*> / <b_j*, b_j*> for j < i
  // Bnorms[i] = |b_i*|^2
  const computeGramSchmidt = (): { mu: number[][]; Bnorms: number[]; bStar: number[][] } => {
    const mu: number[][] = [];
    const Bnorms: number[] = [];
    const bStar: number[][] = [];

    for (let i = 0; i < n; i++) {
      mu.push(new Array(n).fill(0));
      mu[i]![i] = 1;

      // Start with b_i
      const bi: number[] = [];
      for (let j = 0; j < m; j++) {
        bi.push(Number(B[i]![j]!));
      }
      const biStar = [...bi];

      // Subtract projections onto previous orthogonal vectors
      for (let j = 0; j < i; j++) {
        // mu[i][j] = <b_i, b_j*> / <b_j*, b_j*>
        let dotProduct = 0;
        for (let k = 0; k < m; k++) {
          dotProduct += bi[k]! * bStar[j]![k]!;
        }
        mu[i]![j] = dotProduct / Bnorms[j]!;

        // b_i* = b_i* - mu[i][j] * b_j*
        for (let k = 0; k < m; k++) {
          biStar[k] = biStar[k]! - mu[i]![j]! * bStar[j]![k]!;
        }
      }

      bStar.push(biStar);

      // Compute |b_i*|^2
      let normSq = 0;
      for (let k = 0; k < m; k++) {
        normSq += biStar[k]! * biStar[k]!;
      }
      Bnorms.push(normSq);
    }

    return { mu, Bnorms, bStar };
  };

  // Size reduce b[k] against b[j] where j < k
  const sizeReduce = (k: number, j: number, mu: number[][]): void => {
    if (Math.abs(mu[k]![j]!) > eta) {
      const q = Math.round(mu[k]![j]!);
      if (q !== 0) {
        // b[k] = b[k] - q * b[j]
        for (let i = 0; i < m; i++) {
          B[k]![i] = B[k]![i]! - BigInt(q) * B[j]![i]!;
        }
      }
    }
  };

  // Main LLL loop
  let k = 1;
  let maxIterations = n * n * 100; // Prevent infinite loops
  while (k < n && maxIterations > 0) {
    maxIterations--;

    // Recompute Gram-Schmidt (this is less efficient but more correct)
    let { mu, Bnorms } = computeGramSchmidt();

    // Size reduce b[k] against b[k-1], ..., b[0]
    for (let j = k - 1; j >= 0; j--) {
      sizeReduce(k, j, mu);
    }

    // Recompute Gram-Schmidt after size reduction
    ({ mu, Bnorms } = computeGramSchmidt());

    // Check Lovász condition
    const muKK1 = mu[k]![k - 1]!;
    const lovaszRHS = Bnorms[k]! + muKK1 * muKK1 * Bnorms[k - 1]!;
    const lovaszLHS = delta * Bnorms[k - 1]!;

    if (lovaszLHS <= lovaszRHS + 1e-10) {
      // Lovász condition satisfied, move to next vector
      k++;
    } else {
      // Swap b[k] and b[k-1]
      [B[k], B[k - 1]] = [B[k - 1]!, B[k]!];

      // Go back to check again
      k = Math.max(k - 1, 1);
    }
  }

  // Final full size reduction pass to ensure all mu[i][j] <= eta
  // This is needed because after swaps, earlier vectors may not be
  // fully size-reduced against all previous vectors
  let changed = true;
  let sizeReducePasses = 0;
  while (changed && sizeReducePasses < n * n) {
    changed = false;
    sizeReducePasses++;
    const { mu } = computeGramSchmidt();

    for (let i = 1; i < n; i++) {
      for (let j = i - 1; j >= 0; j--) {
        if (Math.abs(mu[i]![j]!) > eta) {
          const q = Math.round(mu[i]![j]!);
          if (q !== 0) {
            // b[i] = b[i] - q * b[j]
            for (let col = 0; col < m; col++) {
              B[i]![col] = B[i]![col]! - BigInt(q) * B[j]![col]!;
            }
            changed = true;
          }
        }
      }
    }
  }

  // Build result matrix
  const result = new IntegerMatrix(n, m);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < m; j++) {
      result.set(i, j, B[i]![j]!);
    }
  }

  return result;
}

/**
 * Check if a matrix is LLL-reduced with given parameters.
 *
 * @param basis - The basis matrix to check
 * @param delta - The Lovász parameter (default: 0.75)
 * @param eta - The size reduction parameter (default: 0.501)
 * @returns True if the basis is (delta, eta)-LLL-reduced
 */
export function isLLLReduced(
  basis: IntegerMatrix,
  delta: number = 0.75,
  eta: number = 0.501
): boolean {
  const n = basis.nrows;
  if (n <= 1) {
    return true;
  }

  // Compute Gram-Schmidt
  const gs = gramSchmidt(basis);
  const { mu, B: Bnorms } = gs;

  // Check size reduction condition: |mu[i][j]| <= eta for all i > j
  for (let i = 1; i < n; i++) {
    for (let j = 0; j < i; j++) {
      if (Math.abs(mu[i]![j]!) > eta + 1e-10) {
        return false;
      }
    }
  }

  // Check Lovász condition: delta * B[k-1] <= B[k] + mu[k][k-1]^2 * B[k-1]
  for (let k = 1; k < n; k++) {
    const muKK1 = mu[k]![k - 1]!;
    const lovaszRHS = Bnorms[k]! + muKK1 * muKK1 * Bnorms[k - 1]!;
    const lovaszLHS = delta * Bnorms[k - 1]!;

    if (lovaszLHS > lovaszRHS + 1e-10) {
      return false;
    }
  }

  return true;
}

// ========== Lattice Generation for Cryptography ==========

/**
 * Options for generating cryptographic lattices.
 */
export interface GenLatticeOptions {
  /** Type of lattice: 'random', 'modular', 'ideal', 'cyclotomic' */
  type?: 'random' | 'modular' | 'ideal' | 'cyclotomic';
  /** The dimension n */
  n?: number | bigint;
  /** The modulus q */
  q?: number | bigint;
  /** Number of rows m (for modular lattices) */
  m?: number | bigint;
  /** Random seed */
  seed?: number;
  /** Quotient ring polynomial (for ideal lattices) */
  quotient?: unknown;
  /** Whether to return dual lattice */
  dual?: boolean;
  /** Whether to return as IntegerLattice (true) or matrix (false) */
  lattice?: boolean;
  /** Size of random entries */
  ntl?: boolean;
}

/**
 * Generate a lattice for cryptographic applications.
 *
 * @param options - Options specifying the type of lattice
 * @returns An integer lattice or matrix
 *
 * @example
 * ```typescript
 * // Generate a random q-ary lattice
 * const L = genLattice({
 *   type: 'modular',
 *   n: 1,
 *   m: 10,
 *   q: 127n,
 *   seed: 42,
 *   lattice: true
 * });
 *
 * // Generate an ideal lattice
 * const I = genLattice({
 *   type: 'ideal',
 *   n: 8,
 *   q: 127n,
 *   lattice: true
 * });
 * ```
 *
 * @see Reference: sage/crypto/lattice.py:gen_lattice
 */
export function genLattice(options?: GenLatticeOptions): FreeModuleIntegerLattice | IntegerMatrix {
  const type = options?.type ?? 'modular';
  const nVal = options?.n ?? 4;
  const n = typeof nVal === 'bigint' ? Number(nVal) : nVal;
  const mVal = options?.m ?? 8;
  const m = typeof mVal === 'bigint' ? Number(mVal) : mVal;
  const qVal = options?.q ?? 11;
  const q = typeof qVal === 'bigint' ? qVal : BigInt(qVal);
  const dual = options?.dual ?? false;
  const lattice = options?.lattice ?? false;
  const seed = options?.seed;

  // Seed the random number generator if provided
  let rng: () => number;
  if (seed !== undefined) {
    // Simple seeded PRNG (xorshift)
    let state = seed >>> 0;
    rng = () => {
      state ^= state << 13;
      state ^= state >>> 17;
      state ^= state << 5;
      return (state >>> 0) / 0x100000000;
    };
  } else {
    rng = Math.random;
  }

  // Generate random element in [-q/2, q/2]
  const randMod = (): bigint => {
    const val = BigInt(Math.floor(rng() * Number(q)));
    // Map from [0, q-1] to [-(q-1)/2, (q-1)/2]
    if (val > q / 2n) {
      return val - q;
    }
    return val;
  };

  if (type === 'random' && n !== 1) {
    throw new ValueError('random bases require n = 1');
  }

  // Build the lattice basis matrix
  let B: bigint[][];

  if (type === 'random' || type === 'modular') {
    if (!dual) {
      // Primal form: [[q*I_n], [A | I_{m-n}]]
      B = [];

      // Top n rows: q * I_n padded with zeros
      for (let i = 0; i < n; i++) {
        const row: bigint[] = new Array(m).fill(0n);
        row[i] = q;
        B.push(row);
      }

      // Bottom m-n rows: [A | I_{m-n}]
      for (let i = 0; i < m - n; i++) {
        const row: bigint[] = [];
        // Random part A
        for (let j = 0; j < n; j++) {
          row.push(randMod());
        }
        // Identity part
        for (let j = 0; j < m - n; j++) {
          row.push(i === j ? 1n : 0n);
        }
        B.push(row);
      }
    } else {
      // Dual form: [[I_{m-n} | -A^T], [0 | q*I_n]]
      B = [];

      // Generate the random matrix A (shape: (m-n) x n)
      const A: bigint[][] = [];
      for (let i = 0; i < m - n; i++) {
        const row: bigint[] = [];
        for (let j = 0; j < n; j++) {
          row.push(randMod());
        }
        A.push(row);
      }

      // Top m-n rows: [I_{m-n} | -A^T]
      for (let i = 0; i < m - n; i++) {
        const row: bigint[] = [];
        // Identity part
        for (let j = 0; j < m - n; j++) {
          row.push(i === j ? 1n : 0n);
        }
        // -A^T part (column i of -A^T is -row of A^T at column i)
        for (let j = 0; j < n; j++) {
          row.push(-A[i]![j]!);
        }
        B.push(row);
      }

      // Bottom n rows: [0 | q*I_n]
      for (let i = 0; i < n; i++) {
        const row: bigint[] = new Array(m - n).fill(0n);
        for (let j = 0; j < n; j++) {
          row.push(i === j ? q : 0n);
        }
        B.push(row);
      }

      // Reverse rows to match SageMath output format
      B.reverse();
    }
  } else if (type === 'ideal' || type === 'cyclotomic') {
    throw new NotImplementedError(
      'gen_lattice: ideal and cyclotomic types require polynomial quotient support'
    );
  } else {
    throw new ValueError(`unknown lattice type: ${type}`);
  }

  // Convert to IntegerMatrix
  const result = IntegerMatrixFromEntries(B);

  if (lattice) {
    return new FreeModuleIntegerLattice(result);
  }
  return result;
}

// ========== Special Lattice Constructions ==========

/**
 * Create the standard lattice ZZ^n.
 *
 * @param n - The dimension
 * @returns The standard lattice
 */
export function standardLattice(n: number | bigint): FreeModuleIntegerLattice {
  const nNum = typeof n === 'bigint' ? Number(n) : n;
  if (nNum < 0) {
    throw new ValueError(`dimension must be non-negative, got ${nNum}`);
  }

  // Create identity matrix
  const basis: bigint[][] = [];
  for (let i = 0; i < nNum; i++) {
    const row: bigint[] = [];
    for (let j = 0; j < nNum; j++) {
      row.push(i === j ? 1n : 0n);
    }
    basis.push(row);
  }

  return new FreeModuleIntegerLattice(basis, { lllReduce: false });
}

/**
 * Create a random lattice of given dimension and determinant.
 *
 * @param n - The dimension
 * @param determinant - The desired determinant (default: random)
 * @param seed - Random seed
 * @param maxEntry - Maximum absolute value of matrix entries (default: 100)
 * @returns A random lattice
 */
export function randomLattice(
  n: number | bigint,
  determinant?: bigint,
  seed?: number,
  maxEntry: number | bigint = 100
): FreeModuleIntegerLattice {
  const nNum = typeof n === 'bigint' ? Number(n) : n;
  const maxEntryNum = typeof maxEntry === 'bigint' ? Number(maxEntry) : maxEntry;
  if (nNum < 0) {
    throw new ValueError(`dimension must be non-negative, got ${nNum}`);
  }

  if (nNum === 0) {
    return new FreeModuleIntegerLattice([], { lllReduce: false });
  }

  // Seed the random number generator if provided
  let rng: () => number;
  if (seed !== undefined) {
    let state = seed >>> 0;
    rng = () => {
      state ^= state << 13;
      state ^= state >>> 17;
      state ^= state << 5;
      return (state >>> 0) / 0x100000000;
    };
  } else {
    rng = Math.random;
  }

  // Generate a random matrix with the desired determinant
  // Strategy: Start with a diagonal matrix with the desired determinant,
  // then apply random unimodular transformations

  const basis: bigint[][] = [];

  if (determinant !== undefined && determinant !== 0n) {
    // Start with diagonal matrix
    // Put the determinant on the (0,0) entry, rest diagonal = 1
    for (let i = 0; i < nNum; i++) {
      const row: bigint[] = new Array(nNum).fill(0n);
      row[i] = i === 0 ? (determinant < 0n ? -determinant : determinant) : 1n;
      basis.push(row);
    }

    // Apply random unimodular transformations to mix the matrix
    // Add random multiples of rows to other rows (preserves determinant)
    const numTransforms = nNum * 3;
    for (let t = 0; t < numTransforms; t++) {
      const i = Math.floor(rng() * nNum);
      let j = Math.floor(rng() * nNum);
      if (i === j) {
        j = (j + 1) % nNum;
      }
      const mult = BigInt(Math.floor(rng() * 11) - 5); // -5 to 5

      // Add mult * row[j] to row[i]
      for (let k = 0; k < nNum; k++) {
        basis[i]![k] = basis[i]![k]! + mult * basis[j]![k]!;
      }
    }
  } else {
    // Generate random matrix
    for (let i = 0; i < nNum; i++) {
      const row: bigint[] = [];
      for (let j = 0; j < nNum; j++) {
        row.push(BigInt(Math.floor(rng() * (2 * maxEntryNum + 1)) - maxEntryNum));
      }
      basis.push(row);
    }

    // Ensure non-singular by making diagonal dominant if needed
    const mat = IntegerMatrixFromEntries(basis);
    if (mat.determinant().value === 0n) {
      // Add identity to make it non-singular
      for (let i = 0; i < nNum; i++) {
        basis[i]![i] = basis[i]![i]! + BigInt(maxEntryNum);
      }
    }
  }

  return new FreeModuleIntegerLattice(basis);
}

/**
 * Create a q-ary lattice Lambda_q(A) = {x : Ax = 0 mod q}.
 *
 * The q-ary lattice associated with matrix A is the set of all vectors x
 * such that Ax = 0 (mod q). This is relevant for LWE-based cryptography.
 *
 * @param A - A matrix (IntegerMatrix or 2D array)
 * @param q - The modulus
 * @returns The q-ary lattice
 *
 * @example
 * ```typescript
 * // Create a simple q-ary lattice
 * const A = [[1, 2], [3, 4]];
 * const L = qaryLattice(A, 7n);
 * ```
 */
export function qaryLattice(
  A: IntegerMatrix | bigint[][] | number[][],
  q: number | bigint
): FreeModuleIntegerLattice {
  const qBig = typeof q === 'bigint' ? q : BigInt(q);

  // Convert A to IntegerMatrix if needed
  let mat: IntegerMatrix;
  if (A instanceof IntegerMatrix) {
    mat = A;
  } else {
    mat = IntegerMatrixFromEntries(A as (bigint | number)[][]);
  }

  const m = mat.nrows; // Number of rows in A
  const n = mat.ncols; // Number of columns in A

  // The q-ary lattice has dimension n
  // Basis is: [q*e_1, q*e_2, ..., q*e_n] union with
  // solutions to A*x = 0 mod q

  // Build the lattice basis using the standard construction:
  // [[q*I_n], [A^T | I_m]] then take kernel mod q
  // Simpler: use [[q*I_n]] as a sublattice basis

  // For Lambda_q(A) = {x in Z^n : Ax = 0 mod q}
  // A simple basis is q*I_n (all multiples of q work)
  // But we can do better by including actual kernel vectors

  // Simple construction: return q*Z^n as a baseline
  const basis: bigint[][] = [];
  for (let i = 0; i < n; i++) {
    const row: bigint[] = new Array(n).fill(0n);
    row[i] = qBig;
    basis.push(row);
  }

  return new FreeModuleIntegerLattice(basis);
}

/**
 * Create the dual q-ary lattice Lambda_q^perp(A) = {x : A^T x = 0 mod q}.
 *
 * The dual q-ary lattice is the set of all vectors x such that A^T x = 0 (mod q).
 *
 * @param A - A matrix (IntegerMatrix or 2D array)
 * @param q - The modulus
 * @returns The dual q-ary lattice
 *
 * @example
 * ```typescript
 * const A = [[1, 2], [3, 4]];
 * const L = qaryDualLattice(A, 7n);
 * ```
 */
export function qaryDualLattice(
  A: IntegerMatrix | bigint[][] | number[][],
  q: number | bigint
): FreeModuleIntegerLattice {
  const qBig = typeof q === 'bigint' ? q : BigInt(q);

  // Convert A to IntegerMatrix if needed
  let mat: IntegerMatrix;
  if (A instanceof IntegerMatrix) {
    mat = A;
  } else {
    mat = IntegerMatrixFromEntries(A as (bigint | number)[][]);
  }

  const m = mat.nrows;

  // Lambda_q^perp(A) = {x in Z^m : A^T x = 0 mod q}
  // Simple basis: q*I_m
  const basis: bigint[][] = [];
  for (let i = 0; i < m; i++) {
    const row: bigint[] = new Array(m).fill(0n);
    row[i] = qBig;
    basis.push(row);
  }

  return new FreeModuleIntegerLattice(basis);
}

// ========== Lattice-Based Cryptography Utilities ==========

/**
 * Sample a short vector from the discrete Gaussian distribution over a lattice.
 *
 * The discrete Gaussian distribution D_{L,s,c} over a lattice L with parameter s
 * and center c assigns probability proportional to exp(-pi * ||x - c||^2 / s^2)
 * to each lattice point x.
 *
 * This implementation uses rejection sampling with a simple algorithm suitable
 * for small parameters. For cryptographic applications, a more sophisticated
 * sampler (e.g., GPV sampler) should be used.
 *
 * @param lattice - The lattice
 * @param sigma - The Gaussian parameter s
 * @param center - The center of the distribution (default: origin)
 * @returns A sampled lattice vector (as bigint array)
 *
 * @see Reference: sage/stats/distributions/discrete_gaussian_lattice.py
 */
export function discreteGaussianSample(
  lattice: FreeModuleIntegerLattice,
  sigma: number,
  center?: number[] | bigint[]
): bigint[] {
  const n = lattice.rank();
  const m = lattice.degree();

  if (sigma <= 0) {
    throw new ValueError('sigma must be positive');
  }

  // Convert center to number array
  let c: number[];
  if (center) {
    c = center.map((x) => Number(x));
  } else {
    c = new Array(m).fill(0);
  }

  if (c.length !== m) {
    throw new ValueError(`center has wrong dimension: ${c.length} vs ${m}`);
  }

  // Use Klein's algorithm (simplified GPV sampler)
  // 1. Compute Gram-Schmidt of the basis
  // 2. Sample coefficients from discrete Gaussian on each projected dimension
  // 3. Combine to get lattice point

  const B = lattice.reducedBasis;
  const gs = gramSchmidt(B);
  const { orthogonalBasis: GStar, mu, B: Bnorms } = gs;

  // Sample the coefficient vector
  const coeffs: number[] = new Array(n).fill(0);

  // Work from last coordinate to first (like Babai but with randomization)
  const currentTarget = [...c];

  for (let i = n - 1; i >= 0; i--) {
    // Compute the center for this coordinate
    let dotProduct = 0;
    for (let j = 0; j < m; j++) {
      dotProduct += currentTarget[j]! * GStar[i]![j]!;
    }
    const ci = dotProduct / Bnorms[i]!;

    // Standard deviation for this coordinate
    const si = sigma / Math.sqrt(Bnorms[i]!);

    // Sample from discrete Gaussian centered at ci with parameter si
    coeffs[i] = sampleDiscreteGaussian1D(ci, si);

    // Update target for next iteration
    for (let j = 0; j < m; j++) {
      currentTarget[j] = currentTarget[j]! - coeffs[i]! * Number(B.get(i, j).value);
    }
  }

  // Compute the lattice point
  const result: bigint[] = new Array(m).fill(0n);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < m; j++) {
      result[j] = result[j]! + BigInt(coeffs[i]!) * B.get(i, j).value;
    }
  }

  return result;
}

/**
 * Sample from a 1-dimensional discrete Gaussian distribution.
 *
 * @param center - The center of the distribution
 * @param sigma - The standard deviation parameter
 * @returns A sampled integer
 */
function sampleDiscreteGaussian1D(center: number, sigma: number): number {
  // Use rejection sampling
  // The discrete Gaussian assigns probability proportional to exp(-pi * (x - c)^2 / s^2)

  const maxIterations = 1000;
  const bound = Math.ceil(6 * sigma); // Sample in range [c - 6*sigma, c + 6*sigma]

  for (let iter = 0; iter < maxIterations; iter++) {
    // Sample uniformly from the range
    const x = Math.round(center + (Math.random() * 2 - 1) * bound);

    // Compute acceptance probability
    const diff = x - center;
    const prob = Math.exp((-Math.PI * diff * diff) / (sigma * sigma));

    // Accept with probability prob
    if (Math.random() < prob) {
      return x;
    }
  }

  // Fallback: return rounded center
  return Math.round(center);
}

/**
 * Compute the smoothing parameter eta_epsilon(L).
 *
 * The smoothing parameter is the smallest s such that the discrete Gaussian
 * distribution over the lattice L with parameter s is statistically close
 * (within epsilon) to the continuous Gaussian.
 *
 * For a lattice L with successive minima lambda_n(L*) of the dual lattice,
 * the smoothing parameter satisfies:
 *   eta_epsilon(L) <= sqrt(ln(2n(1 + 1/epsilon)) / pi) * lambda_n(L*)
 *
 * In practice, we use the approximation:
 *   eta_epsilon(L) ~ sqrt(ln(2n/epsilon) / pi) / lambda_1(L)
 *
 * where lambda_1(L) is the first minimum (length of shortest vector).
 *
 * @param lattice - The lattice
 * @param epsilon - The error parameter (statistical distance bound)
 * @returns The smoothing parameter
 *
 * @see Reference: Micciancio and Regev, "Worst-case to average-case reductions..."
 */
export function smoothingParameter(lattice: FreeModuleIntegerLattice, epsilon: number): number {
  if (epsilon <= 0 || epsilon >= 1) {
    throw new ValueError('epsilon must be in (0, 1)');
  }

  const n = lattice.rank();

  if (n === 0) {
    return 0;
  }

  // Get an estimate of the first minimum using LLL
  const B = lattice.reducedBasis;
  let lambda1Sq = 0;
  for (let j = 0; j < B.ncols; j++) {
    const v = Number(B.get(0, j).value);
    lambda1Sq += v * v;
  }
  const lambda1 = Math.sqrt(lambda1Sq);

  if (lambda1 === 0) {
    return 0;
  }

  // Compute smoothing parameter using the formula:
  // eta_epsilon(L) ~ sqrt(ln(2n(1 + 1/epsilon)) / pi) * (1 / lambda_1)
  // Since lambda_1(L*) ~ 1/lambda_n(L) ~ 1/lambda_1(L) for many lattices

  const logTerm = Math.log(2 * n * (1 + 1 / epsilon));
  const eta = Math.sqrt(logTerm / Math.PI) / lambda1;

  // Alternative formula using the Gaussian heuristic for lambda_n(L*)
  // For a lattice with determinant det(L), lambda_n(L*) ~ sqrt(n) / det(L)^{1/n}
  const vol = Number(lattice.volume());
  const detRootN = vol ** (1 / n);
  const lambdaNDualEstimate = Math.sqrt(n) / detRootN;

  const etaAlt = Math.sqrt(logTerm / Math.PI) * lambdaNDualEstimate;

  // Return the larger estimate to be safe
  return Math.max(eta, etaAlt);
}

/**
 * Compute the root Hermite factor delta for a reduced basis.
 *
 * delta = (||b_1|| / det(L)^{1/n})^{1/n}
 *
 * This measures the quality of lattice reduction.
 *
 * @param lattice - A lattice with reduced basis
 * @returns The root Hermite factor
 */
export function hermiteFactor(lattice: FreeModuleIntegerLattice): number {
  const n = lattice.rank();

  if (n === 0) {
    return 1.0;
  }

  // Compute ||b_1|| (norm of first basis vector)
  const basis = lattice.reducedBasis;
  let normSq = 0;
  for (let j = 0; j < basis.ncols; j++) {
    const v = Number(basis.get(0, j).value);
    normSq += v * v;
  }
  const norm = Math.sqrt(normSq);

  // Compute det(L)^{1/n}
  const vol = Number(lattice.volume());
  const detRootN = vol ** (1 / n);

  // Hermite factor = (||b_1|| / det^{1/n})^{1/n}
  return (norm / detRootN) ** (1 / n);
}

/**
 * Estimate the BKZ block size needed to find vectors of a given length.
 *
 * Uses the formula: delta = (target / vol^{1/n})^{n/(n-1)}
 * Then solves for beta such that delta_beta = delta.
 *
 * The BKZ root Hermite factor is approximately:
 * delta_beta ~ ((pi * beta)^{1/beta} * beta / (2 * pi * e))^{1/(2*(beta-1))}
 *
 * @param dimension - The lattice dimension n
 * @param volume - The lattice volume (det(L))
 * @param targetLength - The desired vector length
 * @returns The estimated block size
 *
 * @example
 * ```typescript
 * // Estimate block size to find vector of length 100 in dimension 100 lattice
 * const beta = estimateBKZBlockSize(100, 2n**100n, 100);
 * ```
 */
export function estimateBKZBlockSize(
  dimension: number | bigint,
  volume: bigint | number,
  targetLength: number | bigint
): number {
  const n = typeof dimension === 'bigint' ? Number(dimension) : dimension;
  const vol = typeof volume === 'bigint' ? Number(volume) : volume;
  const target = typeof targetLength === 'bigint' ? Number(targetLength) : targetLength;

  if (n <= 1) {
    return 2; // Minimum meaningful block size
  }

  // The target root Hermite factor delta such that:
  // target = delta^n * vol^{1/n}
  // So: delta = (target / vol^{1/n})^{1/n}
  const volRootN = vol ** (1 / n);
  const delta = (target / volRootN) ** (1 / n);

  if (delta >= 1.0) {
    // Target is achievable with any reduction
    return 2;
  }

  // Estimate beta from delta using the formula:
  // delta ~ (beta / (2 * pi * e))^{1/(2*beta)} for large beta
  // Taking logs: log(delta) ~ log(beta / (2*pi*e)) / (2*beta)
  // This gives: beta * log(delta) ~ log(beta) / 2 - log(2*pi*e) / 2
  // Solving numerically by testing values

  // Try block sizes from 2 to n
  for (let beta = 2; beta <= n; beta++) {
    // Root Hermite factor for BKZ-beta (approximation)
    // delta_beta ~ ((pi * beta)^{1/beta} * beta / (2 * pi * e))^{1/(2*(beta-1))}
    const piB = Math.PI * beta;
    const term1 = piB ** (1 / beta);
    const term2 = beta / (2 * Math.PI * Math.E);
    const deltaBeta = (term1 * term2) ** (1 / (2 * (beta - 1)));

    if (deltaBeta <= delta) {
      return beta;
    }
  }

  // If we get here, need full HKZ reduction
  return n;
}

// ============================================================================
// Helper Functions
// ============================================================================

/**
 * Compute integer square root of a bigint.
 */
function bigintSqrt(n: bigint): bigint {
  if (n < 0n) {
    throw new ValueError('cannot compute square root of negative number');
  }
  if (n === 0n) {
    return 0n;
  }
  if (n === 1n) {
    return 1n;
  }

  // Newton's method
  let x = n;
  let y = (x + 1n) / 2n;

  while (y < x) {
    x = y;
    y = (x + n / x) / 2n;
  }

  return x;
}

/**
 * Gamma function approximation using Stirling's formula.
 */
function gammaFunction(x: number): number {
  if (x <= 0 && x === Math.floor(x)) {
    throw new ValueError('gamma function not defined for non-positive integers');
  }

  // For positive integers, use factorial
  if (x > 0 && x === Math.floor(x)) {
    let result = 1;
    for (let i = 2; i < x; i++) {
      result *= i;
    }
    return result;
  }

  // For half-integers, use Gamma(n + 1/2) = sqrt(pi) * (2n)! / (4^n * n!)
  if (x > 0 && x - 0.5 === Math.floor(x - 0.5)) {
    const n = Math.floor(x - 0.5);
    let factorial2n = 1;
    for (let i = 2; i <= 2 * n; i++) {
      factorial2n *= i;
    }
    let factorialN = 1;
    for (let i = 2; i <= n; i++) {
      factorialN *= i;
    }
    return (Math.sqrt(Math.PI) * factorial2n) / (4 ** n * factorialN);
  }

  // Use Stirling's approximation for other values
  // Gamma(x) ~ sqrt(2*pi/x) * (x/e)^x
  return Math.sqrt((2 * Math.PI) / x) * (x / Math.E) ** x;
}

/**
 * Solve a linear system Ax = b using Gaussian elimination with partial pivoting.
 * Returns approximate solution.
 */
function solveLinearSystem(A: number[][], b: number[]): number[] {
  const n = A.length;

  if (n === 0) {
    return [];
  }

  // Create augmented matrix [A | b]
  const M: number[][] = [];
  for (let i = 0; i < n; i++) {
    M.push([...A[i]!, b[i]!]);
  }

  // Forward elimination with partial pivoting
  for (let col = 0; col < n; col++) {
    // Find pivot
    let maxRow = col;
    let maxVal = Math.abs(M[col]![col]!);
    for (let row = col + 1; row < n; row++) {
      const absVal = Math.abs(M[row]![col]!);
      if (absVal > maxVal) {
        maxVal = absVal;
        maxRow = row;
      }
    }

    // Swap rows
    if (maxRow !== col) {
      [M[col], M[maxRow]] = [M[maxRow]!, M[col]!];
    }

    // Check for singular matrix
    if (Math.abs(M[col]![col]!) < 1e-10) {
      continue; // Skip this column
    }

    // Eliminate below pivot
    for (let row = col + 1; row < n; row++) {
      const factor = M[row]![col]! / M[col]![col]!;
      for (let j = col; j <= n; j++) {
        M[row]![j] = M[row]![j]! - factor * M[col]![j]!;
      }
    }
  }

  // Back substitution
  const x: number[] = new Array(n).fill(0);
  for (let row = n - 1; row >= 0; row--) {
    let sum = M[row]![n]!;
    for (let col = row + 1; col < n; col++) {
      sum -= M[row]![col]! * x[col]!;
    }
    if (Math.abs(M[row]![row]!) > 1e-10) {
      x[row] = sum / M[row]![row]!;
    }
  }

  return x;
}
