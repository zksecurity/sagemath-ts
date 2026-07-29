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
import {
  IntegerMatrix,
  IntegerMatrixFromEntries,
  hermite_normal_form,
} from '../matrix/matrix_integer.js';
import { FreeModuleGeneric, type FreeModuleOptions, FreeModuleWithBasis } from './free_module.js';
import type { FreeModuleElement } from './free_module_element.js';

/**
 * The exact square root of a positive integer that is **not** a perfect square.
 *
 * `FreeModuleIntegerLattice.volume()` on a lattice whose rank is smaller than
 * its degree returns `sqrt(det(B B^T))`, which is irrational in general.
 * SageMath returns a symbolic `sqrt(N)` there
 * (`free_module_integer.py:494-508`); we have no symbolic ring, so this tiny
 * value type carries the radicand exactly, prints exactly as SageMath does
 * (`sqrt(14)`) and coerces to a double through `valueOf` for the
 * floating-point consumers (`gaussianHeuristic`, `hermiteFactor`, ...).
 *
 * @see Deviation: Lattice Covolume Of Non-Full-Rank Lattices
 */
export class SqrtInteger {
  /** The radicand; never a perfect square. */
  readonly radicand: bigint;

  constructor(radicand: bigint) {
    if (radicand < 0n) {
      throw new ValueError('radicand must be nonnegative');
    }
    this.radicand = radicand;
  }

  /** `sqrt(N)`, matching SageMath's `str()` of the symbolic square root. */
  toString(): string {
    return `sqrt(${this.radicand})`;
  }

  /** The double-precision value, so `Number(v)` and `v ** x` behave. */
  valueOf(): number {
    return Math.sqrt(Number(this.radicand));
  }
}

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

    // Apply LLL reduction if requested (default is true).
    // SageMath: ``basis = matrix([v for v in basis.LLL() if v])`` --- the zero
    // rows LLL produces for a rank-deficient generating set are dropped, so
    // that the lattice always has a genuine basis
    // (free_module_integer.py:304).
    if (options?.lllReduce !== false) {
      const reduced = dropZeroRows(lllReduce(basisMatrix));
      this._basisMatrix = reduced;
      this._reducedBasis = reduced;
      this._basisIsLLLReduced = true;
    } else {
      // `FreeModule_submodule_with_basis_pid.__init__` is invoked with
      // `check=True` (`free_module_integer.py:305-311`) and rejects a
      // linearly dependent basis (`free_module.py:6737-6738`).  With
      // `lllReduce` on, LLL + `dropZeroRows` guarantees independence, so the
      // check can only ever fire on this branch.
      if ((options?.check ?? true) && basisMatrix.nrows > 0) {
        if (basisMatrix.rank() !== basisMatrix.nrows) {
          throw new ValueError('the given basis vectors must be linearly independent');
        }
      }
      this._basisMatrix = basisMatrix;
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
   *
   * Zero rows produced by LLL on a rank-deficient generating set are dropped
   * on construction, so this is the true rank of the lattice.
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
   * Default parameters: delta = 0.99, eta = 0.501 (as in SageMath).
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
    const delta = options?.delta ?? 0.99;
    const eta = options?.eta ?? 0.501;

    // Reduce the current reduced basis (which may already be partially reduced)
    // and drop the zero rows, exactly as SageMath does
    // (free_module_integer.py:395).
    const reduced = dropZeroRows(lllReduce(this._reducedBasis, { delta, eta }));

    // Update the reduced basis if this one is better (shorter first vector)
    if (reduced.nrows > 0) {
      const newNorm = vectorNormSquared(reduced, 0);
      const oldNorm = vectorNormSquared(this._reducedBasis, 0);

      if (newNorm < oldNorm) {
        this._reducedBasis = reduced;
        this._basisIsLLLReduced = true;
      }
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

    // fpLLL accepts block_size = 1 (LLL-only reduction) and rejects 0;
    // `matrix_integer_dense.pyx:2890,2980` passes the value straight through.
    if (blockSize <= 0) {
      throw new ValueError('block size must be > 0');
    }

    // BKZ is a block generalization of LLL
    // For blocks of size 2, it's equivalent to LLL
    // For larger blocks, it uses SVP enumeration within each block

    // Start with LLL-reduced basis (zero rows already dropped)
    let reduced = this.LLL({ delta, eta: 0.501 });

    if (blockSize >= reduced.nrows) {
      // Full block size - do HKZ-style reduction inline to avoid recursion
      // (HKZ calls BKZ with full block size, so we can't delegate back)
      return dropZeroRows(this._fullBlockReduction(reduced, delta));
    }

    // Schnorr-Euchner BKZ: for each kappa, replace `b_kappa` by a shortest
    // vector of the block `[kappa, kappa+beta)` projected away from
    // `b_0..b_{kappa-1}`, found by exact enumeration.  `blockSize === 1` is a
    // pure LLL reduction (fpLLL accepts it), which the LLL above already did.
    if (blockSize >= 2 && reduced.nrows > 1 && blockSize <= EXACT_SVP_MAX_RANK) {
      const rows = bkzTourExact(matrixRows(dropZeroRows(reduced)), blockSize, delta);
      reduced =
        rows.length === 0 ? new IntegerMatrix(0, this.degree()) : IntegerMatrixFromEntries(rows);
    }

    reduced = dropZeroRows(reduced);

    // Update the reduced basis if this is better
    if (reduced.nrows > 0) {
      const newNorm = vectorNormSquared(reduced, 0);
      const oldNorm = vectorNormSquared(this._reducedBasis, 0);

      if (newNorm < oldNorm) {
        this._reducedBasis = reduced;
      }
    }

    return reduced;
  }

  /**
   * Helper for full block reduction (HKZ-style).
   * This is factored out to avoid recursion between BKZ and HKZ.
   */
  private _fullBlockReduction(basis: IntegerMatrix, delta: number): IntegerMatrix {
    const n = basis.nrows;

    if (n <= 1) {
      return basis;
    }

    // A genuine BKZ tour with block size = rank, i.e. HKZ.  The block minima
    // are found by exact enumeration (`blockSVPExact`), so `b_1` really does
    // realise `lambda_1(L)` as `free_module_integer.py:451-478` promises; the
    // old body here only re-ran LLL on suffixes and therefore left strictly
    // shorter vectors sitting further down the basis.
    const reduced = lllReduce(basis, { delta, eta: 0.501 });
    if (n > EXACT_SVP_MAX_RANK) {
      // Above this rank a full enumeration is impractical; fall back to LLL.
      return reduced;
    }
    const rows = bkzTourExact(matrixRows(dropZeroRows(reduced)), n, delta);
    if (rows.length === 0) {
      return new IntegerMatrix(0, basis.ncols);
    }
    return IntegerMatrixFromEntries(rows);
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
    const n = reduced.nrows;

    if (n <= 1) {
      return reduced;
    }

    // HKZ requires that:
    // 1. The basis is size-reduced
    // 2. b_1 is a shortest vector
    // 3. The projection of b_2, ..., b_n onto b_1^perp is HKZ

    // Use full block reduction (same as BKZ with block_size = rank)
    reduced = dropZeroRows(this._fullBlockReduction(reduced, delta));

    // Update reduced basis if better
    if (reduced.nrows > 0) {
      const newNorm = vectorNormSquared(reduced, 0);
      const oldNorm = vectorNormSquared(this._reducedBasis, 0);

      if (newNorm < oldNorm) {
        this._reducedBasis = reduced;
      }
    }

    return reduced;
  }

  // ========== Lattice Problems ==========

  /**
   * Return a shortest nonzero vector in this lattice.
   *
   * This solves the Shortest Vector Problem (SVP) exactly by Fincke-Pohst
   * enumeration on the LLL-reduced basis, up to rank
   * `EXACT_SVP_MAX_RANK`; above that rank the first vector of the reduced
   * basis (an LLL approximation) is returned instead.
   *
   * @param options - Options for the computation
   * @returns A shortest nonzero vector
   *
   * @example
   * ```typescript
   * const L = IntegerLattice(A);
   * const sv = L.shortestVector();
   * ```
   *
   * @see Reference: sage/modules/free_module_integer.py:FreeModule_submodule_with_basis_integer.shortest_vector
   * @see Deviation: SageMath delegates to fpylll's `SVP.shortest_vector` or
   *   PARI's `qfminim`; we enumerate exactly in TypeScript, and fall back to
   *   the LLL approximation above rank EXACT_SVP_MAX_RANK.
   */
  shortestVector(options?: ShortestVectorOptions): bigint[] {
    const updateReducedBasis = options?.updateReducedBasis ?? true;
    const algorithm = options?.algorithm ?? 'fplll';
    if (algorithm !== 'fplll' && algorithm !== 'pari') {
      throw new ValueError(`algorithm '${algorithm}' unknown`);
    }

    // Ensure an LLL-reduced basis with linearly independent rows
    const B = this._basisIsLLLReduced ? this._reducedBasis : this.LLL();
    const rows = matrixRows(B);
    const n = rows.length;

    if (n === 0) {
      throw new ValueError('cannot find shortest vector in zero lattice');
    }

    let sv: bigint[];
    if (n <= EXACT_SVP_MAX_RANK) {
      sv = shortestVectorExact(rows);
    } else {
      // LLL gives a 2^((n-1)/2) approximation.
      sv = [...rows[0]!];
    }

    if (updateReducedBasis) {
      this.updateReducedBasis(sv);
    }

    return sv;
  }

  /**
   * Update the reduced basis with a new lattice vector.
   *
   * SageMath (`free_module_integer.py:609`):
   * ``L = w.stack(self.reduced_basis).LLL(); assert L[0] == 0;
   * self._reduced_basis = L.matrix_from_rows(range(1, L.nrows()))``.
   */
  private _updateBasisWithVector(v: bigint[]): void {
    const n = this._reducedBasis.nrows;
    const m = this._reducedBasis.ncols;

    // Stack w on top of the current basis and LLL-reduce; since w lies in the
    // lattice the result has exactly one zero row, which comes first.
    const stacked: bigint[][] = [[...v], ...matrixRows(this._reducedBasis)];
    const reduced = lllReduce(IntegerMatrixFromEntries(stacked));

    for (let j = 0; j < m; j++) {
      if (reduced.get(0, j).value !== 0n) {
        throw new ValueError('update_reduced_basis: the vector is not in the lattice');
      }
    }

    const rows = matrixRows(reduced).slice(1);
    this._reducedBasis = new IntegerMatrix(n, m, rows);
    this._basisIsLLLReduced = true;
  }

  /**
   * Compute the closest vector in the lattice to a given target vector.
   *
   * This solves the Closest Vector Problem (CVP) exactly.
   *
   * @param target - The target vector (integers, rationals or `"p/q"` strings)
   * @returns The closest lattice vector to target
   *
   * @example
   * ```typescript
   * const L = IntegerLattice([[1, 0], [0, 1]]);
   * const closest = L.closestVector([-6, '5/3']);
   * // Returns (-6, 2)
   * ```
   *
   * @see Reference: sage/modules/free_module_integer.py:FreeModule_submodule_with_basis_integer.closest_vector
   * @see Deviation: SageMath uses the Micciancio-Voulgaris algorithm on the
   *   diamond-cut Voronoi cell; we solve CVP exactly by Fincke-Pohst
   *   enumeration seeded with Babai's nearest plane.  Both return a closest
   *   vector; which one is returned may differ when several are equidistant.
   */
  closestVector(target: FreeModuleElement | unknown[]): bigint[] {
    const t = toRatVector(target);

    const m = this.degree();
    if (t.length !== m) {
      throw new ValueError(`target vector has wrong dimension: ${t.length} vs ${m}`);
    }

    // Ensure an LLL-reduced basis with linearly independent rows
    const B = this._basisIsLLLReduced ? this._reducedBasis : this.LLL();
    const rows = matrixRows(B);
    if (rows.length === 0) {
      return new Array(m).fill(0n);
    }

    // Enumeration minimises |x*B - t|^2 over the lattice; the component of t
    // orthogonal to the lattice span is a constant offset, so projecting t
    // first (as SageMath does) is unnecessary and would only lose exactness.
    const winners = closestVectorsExact(rows, t);
    return winners[0]!;
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
    const algorithm = options?.algorithm ?? 'embedding';
    const delta = options?.delta ?? 0.99;

    const t = toRatVector(target);

    const m = this.degree();
    if (t.length !== m) {
      throw new ValueError(`target vector has wrong dimension: ${t.length} vs ${m}`);
    }

    // Bound checks on delta are performed in isLLLReduced, as in SageMath.
    if (!isLLLReduced(this._reducedBasis, delta)) {
      this.LLL({ delta });
    }

    const B = this._reducedBasis;
    const rows = matrixRows(B);
    const n = rows.length;

    if (n === 0) {
      return new Array(m).fill(0n);
    }

    if (algorithm === 'embedding') {
      // SageMath builds
      //   L = [[B, 0], [t, weight]]  with weight = isqrt(<b_last, b_last>) + 1
      // over QQ, LLL-reduces it and reads off the row whose last entry is
      // +-weight.  A rational matrix is LLL-reduced by clearing denominators
      // (matrix_rational_dense.pyx:3061), and LLL is invariant under scaling
      // the whole matrix, so we scale by the common denominator of t.
      let den = 1n;
      for (const x of t) {
        den = (den / bigGcd(den, x.d)) * x.d;
      }

      const last = rows[n - 1]!;
      const weight = bigintSqrtFloor(dotBig(last, last)) + 1n;

      const entries: bigint[][] = [];
      for (let i = 0; i < n; i++) {
        const row: bigint[] = rows[i]!.map((v) => v * den);
        row.push(0n);
        entries.push(row);
      }
      const tRow: bigint[] = t.map((x) => (x.n * den) / x.d);
      tRow.push(weight * den);
      entries.push(tRow);

      const reduced = lllReduce(IntegerMatrixFromEntries(entries), { delta });

      for (let i = reduced.nrows - 1; i >= 0; i--) {
        const lastEntry = reduced.get(i, m).value;
        if (lastEntry === weight * den || lastEntry === -weight * den) {
          const sign = lastEntry > 0n ? 1n : -1n;
          const result: bigint[] = [];
          for (let j = 0; j < m; j++) {
            const v = ratSub(t[j]!, mkRat(sign * reduced.get(i, j).value, den));
            if (v.d !== 1n) {
              throw new ValueError('embedding CVP produced a non-integral vector');
            }
            result.push(v.n);
          }
          return result;
        }
      }
      throw new ValueError('No suitable vector found in basis.This is a bug, please report it.');
    } else if (algorithm === 'nearest_plane') {
      // Babai's nearest plane, exactly: b -= B[i] * ((b*G[i])/(G[i]*G[i])).round('even')
      return babaiNearestPlane(rows, t);
    } else if (algorithm === 'rounding_off') {
      // t = x*B may have no solution over QQ, so solve x*(B*B^T) = t*B^T
      // (SageMath: ``(B*B.T).solve_left(t*B.T)``) and round to even.
      const gram: Rat[][] = [];
      const rhs: Rat[] = [];
      for (let i = 0; i < n; i++) {
        const row: Rat[] = [];
        for (let j = 0; j < n; j++) {
          row.push(ratFromBigInt(dotBig(rows[i]!, rows[j]!)));
        }
        gram.push(row);
        rhs.push(dotRatBig(t, rows[i]!));
      }
      const sol = solveRationalSystem(gram, rhs);
      const result: bigint[] = new Array(m).fill(0n);
      for (let i = 0; i < n; i++) {
        const c = ratRoundEven(sol[i]!);
        if (c === 0n) continue;
        for (let j = 0; j < m; j++) {
          result[j] = result[j]! + c * rows[i]![j]!;
        }
      }
      return result;
    } else {
      throw new ValueError(
        "algorithm must be one of 'embedding', 'nearest_plane' or 'rounding_off'"
      );
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
  volume(): bigint | SqrtInteger {
    const n = this.rank();
    const m = this.degree();

    if (n === m) {
      // Full rank: volume = |det(B)|
      const det = this._reducedBasis.determinant().value;
      return det < 0n ? -det : det;
    }

    // Non-full rank: volume = sqrt(det(B * B^T)).  Upstream
    // (`free_module_integer.py:494-508`) returns
    // `gram_matrix().determinant().sqrt()`, which is an EXACT value: an
    // Integer when the Gram determinant is a perfect square and a symbolic
    // `sqrt(N)` otherwise.  Flooring the integer square root here (as this used
    // to) silently returned a wrong answer -- `<(1,1,0),(0,1,1)>` has covolume
    // sqrt(3), not 1.
    const B = this._reducedBasis;
    const BT = B.transpose();
    const gram = B.mul(BT);
    const det = gram.determinant().value;
    const absDet = det < 0n ? -det : det;

    const r = bigintSqrt(absDet);
    if (r * r === absDet) {
      return r;
    }
    return new SqrtInteger(absDet);
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
    // `volume() == 1` in Sage; a symbolic `sqrt(N)` with N not a perfect square
    // is never equal to 1, so a non-bigint volume is immediately False.
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
   * The Voronoi cell is the set of all points at least as close to the origin
   * as to any other lattice point.  It is the (bounded, for a full rank
   * lattice) polyhedron cut out by the half-spaces
   * `2<v, x> <= |v|^2` for the Voronoi-relevant vectors `v`; those
   * inequalities are exactly the facets of the cell.
   *
   * @param radius - unused; present for signature compatibility with SageMath
   * @returns The Voronoi cell as its (irredundant) H-representation, with each
   *   inequality `normals[i] . x <= offsets[i]` scaled down by the gcd of its
   *   coefficients
   *
   * @see Reference: sage/modules/free_module_integer.py:FreeModule_submodule_with_basis_integer.voronoi_cell
   * @see Deviation: SageMath returns a `Polyhedron` obtained by diamond
   *   cutting; we return the equivalent exact H-representation.
   */
  voronoiCell(radius?: number): { normals: bigint[][]; offsets: bigint[] } {
    void radius;

    // The Voronoi cell is {x : |x| <= |x - v| for all lattice vectors v},
    // i.e. {x : 2<v, x> <= |v|^2}.  Only the Voronoi-relevant vectors give
    // facets; every other lattice vector yields a redundant inequality.
    const relevantVectors = this.voronoiRelevantVectors();

    const normals: bigint[][] = [];
    const offsets: bigint[] = [];

    for (const v of relevantVectors) {
      const normal = v.map((x) => 2n * x);
      const offset = dotBig(v, v);

      // Normalize by the gcd, as a Polyhedron's H-representation is.
      let g = offset;
      for (const c of normal) {
        g = bigGcd(g, c);
      }
      if (g > 1n) {
        normals.push(normal.map((c) => c / g));
        offsets.push(offset / g);
      } else {
        normals.push(normal);
        offsets.push(offset);
      }
    }

    return { normals, offsets };
  }

  /**
   * Compute the Voronoi-relevant vectors (those defining the Voronoi cell).
   *
   * By Voronoi's theorem a nonzero lattice vector `v` is relevant if and only
   * if `+v` and `-v` are the *only* minimal length vectors of the coset
   * `v + 2L`.  There are at most `2(2^r - 1)` of them, one pair per nonzero
   * coset of `L / 2L`, and the set is closed under negation.
   *
   * @returns List of Voronoi-relevant vectors
   *
   * @example
   * ```typescript
   * IntegerLattice([[3, 0], [4, 0]]).voronoiRelevantVectors();
   * // [[-1, 0], [1, 0]]
   * ```
   *
   * @see Reference: sage/modules/free_module_integer.py:FreeModule_submodule_with_basis_integer.voronoi_relevant_vectors
   * @see Deviation: SageMath reads the relevant vectors off the diamond-cut
   *   Voronoi cell; we use Voronoi's characterization directly.  Both are
   *   exact and give the same set.
   */
  voronoiRelevantVectors(): bigint[][] {
    const B = this._basisIsLLLReduced ? this._reducedBasis : this.LLL();
    const rows = matrixRows(B);
    const r = rows.length;
    if (r === 0) {
      return [];
    }
    const m = B.ncols;

    // Basis of 2L
    const doubled = rows.map((row) => row.map((v) => 2n * v));

    const relevant: bigint[][] = [];

    // Walk over the 2^r - 1 nonzero cosets of L / 2L.
    const total = 1 << r;
    if (r > 24) {
      throw new NotImplementedError(
        `voronoi_relevant_vectors: rank ${r} is too large to enumerate L/2L`
      );
    }
    for (let mask = 1; mask < total; mask++) {
      const c: bigint[] = new Array(m).fill(0n);
      for (let i = 0; i < r; i++) {
        if ((mask >> i) & 1) {
          for (let j = 0; j < m; j++) {
            c[j] = c[j]! + rows[i]![j]!;
          }
        }
      }

      // Minimal length vectors of the coset c + 2L: c - w for the w in 2L
      // closest to c.
      const cRat = c.map((x) => ratFromBigInt(x));
      const closest = closestVectorsExact(doubled, cRat);
      if (closest.length === 2) {
        for (const w of closest) {
          relevant.push(c.map((x, j) => x - w[j]!));
        }
      }
    }

    // Deterministic (lexicographic) order, matching the order in which
    // SageMath's doctest reports them: [(-1, 0), (1, 0)].
    relevant.sort((a, b) => {
      for (let j = 0; j < a.length; j++) {
        if (a[j]! !== b[j]!) return a[j]! < b[j]! ? -1 : 1;
      }
      return 0;
    });

    return relevant;
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

// ========== Exact rational arithmetic ==========
//
// SageMath performs lattice reduction and Babai rounding over exact rings
// (ZZ / QQ, via fpLLL's exact wrapper and `Rational.round`).  Everything in
// this section is exact bigint arithmetic so that the algorithms below never
// go through IEEE doubles.

/** An exact rational number n/d with d > 0 and gcd(n, d) = 1. */
export interface Rat {
  n: bigint;
  d: bigint;
}

const RAT_ZERO: Rat = { n: 0n, d: 1n };
const RAT_ONE: Rat = { n: 1n, d: 1n };

function bigAbs(a: bigint): bigint {
  return a < 0n ? -a : a;
}

function bigGcd(a: bigint, b: bigint): bigint {
  let x = bigAbs(a);
  let y = bigAbs(b);
  while (y !== 0n) {
    const t = x % y;
    x = y;
    y = t;
  }
  return x;
}

/** Integer square root (floor) of a nonnegative bigint. */
function bigintSqrtFloor(n: bigint): bigint {
  if (n < 0n) {
    throw new ValueError('cannot compute square root of negative number');
  }
  if (n < 2n) {
    return n;
  }
  let x = n;
  let y = (x + 1n) / 2n;
  while (y < x) {
    x = y;
    y = (x + n / x) / 2n;
  }
  return x;
}

function mkRat(n: bigint, d: bigint): Rat {
  if (d === 0n) {
    throw new ValueError('rational with zero denominator');
  }
  let nn = n;
  let dd = d;
  if (dd < 0n) {
    nn = -nn;
    dd = -dd;
  }
  if (nn === 0n) {
    return RAT_ZERO;
  }
  const g = bigGcd(nn, dd);
  return { n: nn / g, d: dd / g };
}

function ratFromBigInt(n: bigint): Rat {
  return { n, d: 1n };
}

function ratAdd(a: Rat, b: Rat): Rat {
  return mkRat(a.n * b.d + b.n * a.d, a.d * b.d);
}

function ratSub(a: Rat, b: Rat): Rat {
  return mkRat(a.n * b.d - b.n * a.d, a.d * b.d);
}

function ratMul(a: Rat, b: Rat): Rat {
  return mkRat(a.n * b.n, a.d * b.d);
}

function ratDiv(a: Rat, b: Rat): Rat {
  if (b.n === 0n) {
    throw new ValueError('division by zero rational');
  }
  return mkRat(a.n * b.d, a.d * b.n);
}

function ratAbs(a: Rat): Rat {
  return a.n < 0n ? { n: -a.n, d: a.d } : a;
}

function ratIsZero(a: Rat): boolean {
  return a.n === 0n;
}

/** Return -1, 0 or 1 according to a < b, a == b, a > b. */
function ratCmp(a: Rat, b: Rat): number {
  const lhs = a.n * b.d;
  const rhs = b.n * a.d;
  return lhs < rhs ? -1 : lhs > rhs ? 1 : 0;
}

/** Floor division for bigints (Python semantics). */
function floorDivBig(a: bigint, b: bigint): bigint {
  let q = a / b;
  if (a % b !== 0n && a < 0n !== b < 0n) {
    q -= 1n;
  }
  return q;
}

/** floor(x) for an exact rational. */
function ratFloor(a: Rat): bigint {
  return floorDivBig(a.n, a.d);
}

/**
 * Round to nearest integer, ties away from the origin never occur: ties go to
 * the even integer.  This is SageMath's ``Rational.round('even')``.
 */
function ratRoundEven(a: Rat): bigint {
  const q = ratFloor(a);
  const r = a.n - q * a.d; // 0 <= r < d
  const twice = 2n * r;
  if (twice < a.d) {
    return q;
  }
  if (twice > a.d) {
    return q + 1n;
  }
  // Exact tie: pick the even one.
  return ((q % 2n) + 2n) % 2n === 0n ? q : q + 1n;
}

/** Round to nearest integer, ties up: floor(x + 1/2).  Used by LLL's RED step. */
function ratRoundHalfUp(a: Rat): bigint {
  return floorDivBig(2n * a.n + a.d, 2n * a.d);
}

/**
 * Convert a JavaScript number to the exact rational it represents.
 *
 * IEEE doubles are dyadic rationals, so this conversion is lossless: no
 * information is invented and none is dropped.
 */
function ratFromNumber(x: number): Rat {
  if (!Number.isFinite(x)) {
    throw new ValueError(`cannot convert ${x} to a rational number`);
  }
  if (Number.isInteger(x)) {
    return ratFromBigInt(BigInt(x));
  }
  let v = x;
  let d = 1n;
  while (!Number.isInteger(v)) {
    v *= 2;
    d *= 2n;
  }
  return mkRat(BigInt(v), d);
}

/**
 * Coerce an arbitrary scalar to an exact rational.
 *
 * Accepts bigint, number, a `"p/q"` string, `{n, d}` pairs and any object
 * exposing a bigint `value` field (e.g. `Integer`).
 */
export function toRat(x: unknown): Rat {
  if (typeof x === 'bigint') {
    return ratFromBigInt(x);
  }
  if (typeof x === 'number') {
    return ratFromNumber(x);
  }
  if (typeof x === 'string') {
    const parts = x.split('/');
    if (parts.length === 1) {
      return ratFromBigInt(BigInt(parts[0]!.trim()));
    }
    if (parts.length === 2) {
      return mkRat(BigInt(parts[0]!.trim()), BigInt(parts[1]!.trim()));
    }
    throw new ValueError(`cannot convert ${x} to a rational number`);
  }
  if (x !== null && typeof x === 'object') {
    const o = x as { n?: unknown; d?: unknown; value?: unknown };
    if (typeof o.n === 'bigint' && typeof o.d === 'bigint') {
      return mkRat(o.n, o.d);
    }
    if (typeof o.value === 'bigint') {
      return ratFromBigInt(o.value);
    }
    if (typeof o.value === 'number') {
      return ratFromNumber(o.value);
    }
  }
  throw new ValueError(`cannot convert ${String(x)} to a rational number`);
}

/** Exact dot product of two integer vectors. */
function dotBig(a: bigint[], b: bigint[]): bigint {
  let s = 0n;
  for (let i = 0; i < a.length; i++) {
    s += a[i]! * b[i]!;
  }
  return s;
}

/** Exact dot product of a rational vector with an integer vector. */
function dotRatBig(a: Rat[], b: bigint[]): Rat {
  let s = RAT_ZERO;
  for (let i = 0; i < a.length; i++) {
    s = ratAdd(s, ratMul(a[i]!, ratFromBigInt(b[i]!)));
  }
  return s;
}

/** Exact dot product of two rational vectors. */
function dotRat(a: Rat[], b: Rat[]): Rat {
  let s = RAT_ZERO;
  for (let i = 0; i < a.length; i++) {
    s = ratAdd(s, ratMul(a[i]!, b[i]!));
  }
  return s;
}

/**
 * Exact Gram-Schmidt data of an integer basis.
 *
 * `mu[i][j] = <b_i, b_j*> / <b_j*, b_j*>` for `j < i` and `B[i] = |b_i*|^2`,
 * computed from the integer Gram matrix so that no floating point is involved.
 * If `b_j*` vanishes (linearly dependent rows) then `<b_i, b_j*> = 0` and
 * `mu[i][j]` is defined to be 0.
 */
export interface ExactGSO {
  mu: Rat[][];
  B: Rat[];
}

export function exactGramSchmidt(rows: bigint[][]): ExactGSO {
  const n = rows.length;
  const mu: Rat[][] = [];
  const B: Rat[] = [];

  for (let i = 0; i < n; i++) {
    const muI: Rat[] = new Array(n).fill(RAT_ZERO);
    muI[i] = RAT_ONE;
    let Bi = ratFromBigInt(dotBig(rows[i]!, rows[i]!));

    for (let j = 0; j < i; j++) {
      // <b_i, b_j*> = <b_i, b_j> - sum_{k<j} mu[j][k] * <b_i, b_k*>
      //             = <b_i, b_j> - sum_{k<j} mu[j][k] * mu[i][k] * B[k]
      let s = ratFromBigInt(dotBig(rows[i]!, rows[j]!));
      for (let k = 0; k < j; k++) {
        s = ratSub(s, ratMul(mu[j]![k]!, ratMul(muI[k]!, B[k]!)));
      }
      muI[j] = ratIsZero(B[j]!) ? RAT_ZERO : ratDiv(s, B[j]!);
      Bi = ratSub(Bi, ratMul(ratMul(muI[j]!, muI[j]!), B[j]!));
    }

    mu.push(muI);
    B.push(Bi);
  }

  return { mu, B };
}

/**
 * Exact Fincke-Pohst enumeration.
 *
 * Enumerate every lattice vector `v = sum x_i b_i` of the lattice spanned by
 * `rows` with `|v - t|^2 <= bound`, calling `visit` for each one.  All
 * arithmetic is exact.
 *
 * @param rows - basis of the lattice (linearly independent rows)
 * @param t - target vector (exact rationals); use the zero vector for SVP
 * @param bound - squared radius (exact rational)
 * @param visit - callback receiving the vector and its squared distance to `t`
 */
function enumerateBall(
  rows: bigint[][],
  t: Rat[],
  bound: Rat,
  visit: (v: bigint[], distSq: Rat) => void
): void {
  const n = rows.length;
  if (n === 0) {
    return;
  }
  const m = rows[0]!.length;
  const { mu, B } = exactGramSchmidt(rows);
  for (let i = 0; i < n; i++) {
    if (ratIsZero(B[i]!)) {
      // Without this the enumeration at level i would never terminate.
      throw new ValueError('linearly dependent input for module version of Gram-Schmidt');
    }
  }

  // Coordinates of t with respect to the Gram-Schmidt basis:
  // tau_i = <t, b_i*> / B_i, with <t, b_i*> = <t, b_i> - sum_{k<i} mu[i][k] * <t, b_k*>
  const tau: Rat[] = [];
  const tDotStar: Rat[] = [];
  for (let i = 0; i < n; i++) {
    let s = dotRatBig(t, rows[i]!);
    for (let k = 0; k < i; k++) {
      s = ratSub(s, ratMul(mu[i]![k]!, tDotStar[k]!));
    }
    tDotStar.push(s);
    tau.push(ratIsZero(B[i]!) ? RAT_ZERO : ratDiv(s, B[i]!));
  }

  // Squared length of the component of t orthogonal to the lattice span.
  let ortho = dotRat(t, t);
  for (let i = 0; i < n; i++) {
    ortho = ratSub(ortho, ratMul(ratMul(tau[i]!, tau[i]!), B[i]!));
  }
  if (ratCmp(ortho, bound) > 0) {
    return; // nothing can be within the bound
  }

  const x: bigint[] = new Array(n).fill(0n);

  const recurse = (i: number, remaining: Rat): void => {
    // center_i = tau_i - sum_{j>i} mu[j][i] * x_j
    let center = tau[i]!;
    for (let j = i + 1; j < n; j++) {
      center = ratSub(center, ratMul(mu[j]![i]!, ratFromBigInt(x[j]!)));
    }

    const start = ratRoundEven(center);
    // Walk outwards from the centre in both directions while the projected
    // contribution B_i * (x_i - center)^2 still fits in the remaining budget.
    for (let dir = 0; dir < 2; dir++) {
      let step = dir === 0 ? 0n : -1n;
      for (;;) {
        const xi = start + step;
        const diff = ratSub(ratFromBigInt(xi), center);
        const contrib = ratMul(ratMul(diff, diff), B[i]!);
        if (ratCmp(contrib, remaining) > 0) {
          break;
        }
        x[i] = xi;
        if (i === 0) {
          const v: bigint[] = new Array(m).fill(0n);
          for (let r = 0; r < n; r++) {
            const c = x[r]!;
            if (c === 0n) continue;
            for (let col = 0; col < m; col++) {
              v[col] = v[col]! + c * rows[r]![col]!;
            }
          }
          // Exact squared distance |v - t|^2
          let distSq = RAT_ZERO;
          for (let col = 0; col < m; col++) {
            const d = ratSub(ratFromBigInt(v[col]!), t[col]!);
            distSq = ratAdd(distSq, ratMul(d, d));
          }
          visit(v, distSq);
        } else {
          recurse(i - 1, ratSub(remaining, contrib));
        }
        step = dir === 0 ? step + 1n : step - 1n;
      }
    }
    x[i] = 0n;
  };

  recurse(n - 1, ratSub(bound, ortho));
}

/**
 * Exact solution of the closest vector problem by enumeration.
 *
 * @param rows - basis rows (linearly independent)
 * @param t - target (exact rationals)
 * @returns every lattice vector at minimal distance from `t`
 */
function closestVectorsExact(rows: bigint[][], t: Rat[]): bigint[][] {
  const n = rows.length;
  if (n === 0) {
    return [];
  }

  // Babai nearest plane gives a starting radius.
  const start = babaiNearestPlane(rows, t);
  let best = RAT_ZERO;
  for (let j = 0; j < t.length; j++) {
    const d = ratSub(ratFromBigInt(start[j]!), t[j]!);
    best = ratAdd(best, ratMul(d, d));
  }

  let winners: bigint[][] = [start];

  enumerateBall(rows, t, best, (v, distSq) => {
    if (ratCmp(distSq, best) < 0) {
      best = distSq;
      winners = [v];
    } else if (ratCmp(distSq, best) === 0) {
      if (!winners.some((w) => w.every((c, idx) => c === v[idx]))) {
        winners.push(v);
      }
    }
  });

  return winners;
}

/**
 * Largest rank for which `shortestVector` runs a full exact enumeration.
 *
 * Above this rank the enumeration tree becomes impractical in TypeScript and
 * the LLL approximation is returned instead.
 */
const EXACT_SVP_MAX_RANK = 30;

/**
 * Exact shortest vector by Fincke-Pohst enumeration.
 *
 * The squared radius starts at the shortest row of the (LLL-reduced) basis,
 * so the enumeration is guaranteed to see a shortest vector.
 */
function shortestVectorExact(rows: bigint[][]): bigint[] {
  const n = rows.length;
  const m = rows[0]!.length;

  // Track the argmin row alongside the minimum, so that `best` and `bound`
  // describe the SAME vector.  Seeding `best` with row 0 while `bound` came
  // from a different row made `enumerateBall` (which only replaces `best` on a
  // strictly shorter hit) return row 0 whenever the minimum was already
  // optimal, i.e. return a vector that is not shortest.
  let bound = dotBig(rows[0]!, rows[0]!);
  let boundIdx = 0;
  for (let i = 1; i < n; i++) {
    const norm = dotBig(rows[i]!, rows[i]!);
    if (norm !== 0n && (bound === 0n || norm < bound)) {
      bound = norm;
      boundIdx = i;
    }
  }

  let best: bigint[] = [...rows[boundIdx]!];
  let bestNorm = ratFromBigInt(bound);
  const zero: Rat[] = new Array(m).fill(RAT_ZERO);

  enumerateBall(rows, zero, ratFromBigInt(bound), (v, distSq) => {
    if (ratIsZero(distSq)) {
      return; // the zero vector
    }
    if (ratCmp(distSq, bestNorm) < 0) {
      bestNorm = distSq;
      best = v;
    }
  });

  return best;
}

/**
 * Babai's nearest plane algorithm over exact rationals.
 *
 * Mirrors `free_module_integer.py:approximate_closest_vector` with
 * ``algorithm='nearest_plane'``: `b -= B[i] * ((b*G[i])/(G[i]*G[i])).round('even')`.
 */
function babaiNearestPlane(rows: bigint[][], t: Rat[]): bigint[] {
  const n = rows.length;
  if (n === 0) {
    return t.map(() => 0n);
  }
  const m = rows[0]!.length;
  const { mu, B } = exactGramSchmidt(rows);

  // G[i] = b_i* expressed exactly: b_i* = b_i - sum_{j<i} mu[i][j] b_j*
  const G: Rat[][] = [];
  for (let i = 0; i < n; i++) {
    const gi: Rat[] = rows[i]!.map((v) => ratFromBigInt(v));
    for (let j = 0; j < i; j++) {
      for (let k = 0; k < m; k++) {
        gi[k] = ratSub(gi[k]!, ratMul(mu[i]![j]!, G[j]![k]!));
      }
    }
    G.push(gi);
  }

  const b: Rat[] = [...t];
  for (let i = n - 1; i >= 0; i--) {
    if (ratIsZero(B[i]!)) {
      continue;
    }
    const c = ratRoundEven(ratDiv(dotRat(b, G[i]!), B[i]!));
    if (c !== 0n) {
      for (let k = 0; k < m; k++) {
        b[k] = ratSub(b[k]!, ratMul(ratFromBigInt(c), ratFromBigInt(rows[i]![k]!)));
      }
    }
  }

  // The lattice vector is t - b, which is integral.
  const result: bigint[] = [];
  for (let k = 0; k < b.length; k++) {
    const v = ratSub(t[k]!, b[k]!);
    if (v.d !== 1n) {
      throw new ValueError('Babai nearest plane produced a non-integral vector');
    }
    result.push(v.n);
  }
  return result;
}

/**
 * Solve the (square, symmetric) rational system `x * A = b` exactly by
 * Gaussian elimination with exact rational pivoting.
 *
 * Used for SageMath's ``rounding_off`` Babai variant, which solves
 * `x * (B B^T) = t B^T`.
 */
function solveRationalSystem(A: Rat[][], b: Rat[]): Rat[] {
  const n = A.length;
  if (n === 0) {
    return [];
  }
  const M: Rat[][] = A.map((row, i) => [...row, b[i]!]);

  let row = 0;
  const pivotOf: number[] = new Array(n).fill(-1);
  for (let col = 0; col < n && row < n; col++) {
    let pivot = -1;
    for (let r = row; r < n; r++) {
      if (!ratIsZero(M[r]![col]!)) {
        pivot = r;
        break;
      }
    }
    if (pivot === -1) {
      continue;
    }
    if (pivot !== row) {
      const t = M[row]!;
      M[row] = M[pivot]!;
      M[pivot] = t;
    }
    const p = M[row]![col]!;
    for (let j = col; j <= n; j++) {
      M[row]![j] = ratDiv(M[row]![j]!, p);
    }
    for (let r = 0; r < n; r++) {
      if (r === row) continue;
      const f = M[r]![col]!;
      if (ratIsZero(f)) continue;
      for (let j = col; j <= n; j++) {
        M[r]![j] = ratSub(M[r]![j]!, ratMul(f, M[row]![j]!));
      }
    }
    pivotOf[col] = row;
    row++;
  }

  const x: Rat[] = new Array(n).fill(RAT_ZERO);
  for (let col = 0; col < n; col++) {
    if (pivotOf[col]! >= 0) {
      x[col] = M[pivotOf[col]!]![n]!;
    }
  }
  return x;
}

/** Extract the rows of an IntegerMatrix as bigint arrays. */
function matrixRows(M: IntegerMatrix): bigint[][] {
  const rows: bigint[][] = [];
  for (let i = 0; i < M.nrows; i++) {
    const row: bigint[] = [];
    for (let j = 0; j < M.ncols; j++) {
      row.push(M.get(i, j).value);
    }
    rows.push(row);
  }
  return rows;
}

/**
 * Exact shortest vector of a **projected** block lattice.
 *
 * Let `pi_kappa` be the orthogonal projection away from `b_0, ..., b_{kappa-1}`.
 * This enumerates every nonzero integer combination
 * `x = (x_kappa, ..., x_{endIdx-1})` with
 *
 * ```
 * || sum_i x_i pi_kappa(b_i) ||^2  <  |b_kappa*|^2
 * ```
 *
 * and returns the minimising `x` (indexed from `kappa`), or `null` when no
 * combination beats `|b_kappa*|^2` — i.e. when the block is already
 * Korkine-Zolotarev reduced at `kappa`.
 *
 * All arithmetic is exact: the projected Gram-Schmidt data of the block is
 * literally `mu[i][j]` and `B[i]` for `kappa <= j <= i < endIdx`, which
 * {@link exactGramSchmidt} produces over the rationals.  The previous
 * double-precision enumeration in `bkz.ts` cut branches off using
 * `|x_k| > ceil(sqrt(bound / B[k])) + 1`, a test that ignores both the
 * enumeration centre and the norm already accumulated at deeper levels, and so
 * missed block minima outright.
 */
function blockSVPExact(
  rows: bigint[][],
  gso: ExactGSO,
  kappa: number,
  endIdx: number
): bigint[] | null {
  const d = endIdx - kappa;
  if (d < 2) {
    return null;
  }

  const { mu, B } = gso;
  const bound = B[kappa]!;
  if (ratCmp(bound, RAT_ZERO) <= 0) {
    return null;
  }

  const x: bigint[] = new Array(d).fill(0n);
  let best: bigint[] | null = null;
  let bestNorm = bound;

  // rho[i] = squared norm accumulated by levels i..d-1
  const rho: Rat[] = new Array(d + 1).fill(RAT_ZERO);

  // The enumeration radius is `|b_kappa*|`, which after LLL is tight, so the
  // tree is small in practice.  The budget only guards against a pathological
  // input; when it is hit the best vector found so far is still a genuine
  // improvement, it merely may not be the block minimum.
  let nodes = 0;
  const NODE_BUDGET = 5_000_000;

  /** Squared contribution of level i for the coefficient `k` about centre `c`. */
  const levelNorm = (i: number, c: Rat, k: bigint): Rat => {
    const diff = ratSub(ratFromBigInt(k), c);
    return ratMul(ratMul(diff, diff), B[kappa + i]!);
  };

  const descend = (i: number): void => {
    // Centre: c_i = -sum_{l>i} x_l * mu[kappa+l][kappa+i]
    let c: Rat = RAT_ZERO;
    for (let l = i + 1; l < d; l++) {
      if (x[l] !== 0n) {
        c = ratSub(c, ratMul(ratFromBigInt(x[l]!), mu[kappa + l]![kappa + i]!));
      }
    }
    const centre = ratRoundHalfUp(c);

    // `levelNorm` is a strictly convex quadratic in the integer coefficient and
    // `centre` is the nearest integer to `c`, so |k - c| grows monotonically as
    // we walk outward in either direction: the first coefficient that busts the
    // bound proves every further one in that direction busts it too.
    const visit = (k: bigint): boolean => {
      if (nodes++ > NODE_BUDGET) {
        return false;
      }
      const total = ratAdd(rho[i + 1]!, levelNorm(i, c, k));
      if (ratCmp(total, bestNorm) >= 0) {
        return false;
      }
      x[i] = k;
      if (i === 0) {
        if (x.some((v) => v !== 0n)) {
          bestNorm = total;
          best = [...x];
        }
      } else {
        rho[i] = total;
        descend(i - 1);
      }
      return true;
    };

    visit(centre);
    for (let k = centre + 1n; visit(k); k++) {
      // walk up
    }
    for (let k = centre - 1n; visit(k); k--) {
      // walk down
    }
    x[i] = 0n;
  };

  rho[d] = RAT_ZERO;
  descend(d - 1);

  return best;
}

/**
 * A `d x d` unimodular integer matrix whose first row is `x`.
 *
 * Such a matrix exists exactly when `gcd(x) == 1`, which holds for the
 * coefficient vector of any shortest vector of a block (otherwise `v / gcd(x)`
 * would be a strictly shorter lattice vector).
 *
 * Built by clearing `x` to `(1, 0, ..., 0)` with unimodular COLUMN operations
 * `x . T = e_1`, while accumulating `T^{-1}` through the mirrored ROW
 * operations; the answer is `T^{-1}`, whose first row is `e_1 . T^{-1} = x`.
 */
function unimodularWithFirstRow(x: bigint[]): bigint[][] {
  const d = x.length;
  const Uinv: bigint[][] = [];
  for (let i = 0; i < d; i++) {
    Uinv.push(Array.from({ length: d }, (_, j) => (i === j ? 1n : 0n)));
  }

  const y = [...x];
  for (let i = 1; i < d; i++) {
    if (y[i] === 0n) {
      continue;
    }
    const [g, s, t] = xgcdBig(y[0]!, y[i]!);
    const a = y[0]! / g;
    const b = y[i]! / g;
    // Column op: [col0 col_i] <- [col0 col_i] * [[s, -b], [t, a]] (det 1).
    // Its inverse acts on the rows of `Uinv` as [[a, b], [-t, s]].
    const row0 = Uinv[0]!;
    const rowI = Uinv[i]!;
    for (let j = 0; j < d; j++) {
      const p = row0[j]!;
      const q = rowI[j]!;
      row0[j] = a * p + b * q;
      rowI[j] = -t * p + s * q;
    }
    y[0] = g;
    y[i] = 0n;
  }

  if (y[0] !== 1n) {
    if (y[0] !== -1n) {
      throw new ValueError('block minimum has non-primitive coefficient vector');
    }
    // gcd came out as -1: negate the first row (and one other) to keep det ±1
    // while making the first row exactly `x`.
    for (let j = 0; j < d; j++) {
      Uinv[0]![j] = -Uinv[0]![j]!;
    }
  }

  return Uinv;
}

/** Extended gcd on bigints, returning `[g, s, t]` with `s*a + t*b == g >= 0`. */
function xgcdBig(a: bigint, b: bigint): [bigint, bigint, bigint] {
  let [oldR, r] = [a, b];
  let [oldS, s] = [1n, 0n];
  let [oldT, t] = [0n, 1n];
  while (r !== 0n) {
    const q = (oldR - (((oldR % r) + r) % r)) / r;
    [oldR, r] = [r, oldR - q * r];
    [oldS, s] = [s, oldS - q * s];
    [oldT, t] = [t, oldT - q * t];
  }
  if (oldR < 0n) {
    return [-oldR, -oldS, -oldT];
  }
  return [oldR, oldS, oldT];
}

/**
 * One BKZ tour with exact block SVP, followed by LLL.
 *
 * This is the Schnorr-Euchner BKZ loop: for every `kappa`, find a shortest
 * vector of the block `[kappa, kappa+beta)` projected away from
 * `b_0..b_{kappa-1}`; if it is strictly shorter than `b_kappa*`, insert the
 * corresponding lattice vector at position `kappa` and LLL the resulting
 * (redundant) generating set back down to a basis.
 *
 * @returns `true` when the basis changed.
 */
function bkzTourExact(rows: bigint[][], beta: number, delta: number): bigint[][] {
  let current = rows;
  const n = current.length;
  if (n < 2) {
    return current;
  }

  const maxTours = 4 * n + 8;
  for (let tour = 0; tour < maxTours; tour++) {
    let changed = false;
    for (let kappa = 0; kappa + 1 < current.length; kappa++) {
      const endIdx = Math.min(kappa + beta, current.length);
      if (endIdx - kappa < 2) {
        continue;
      }
      const gso = exactGramSchmidt(current);
      const x = blockSVPExact(current, gso, kappa, endIdx);
      if (x === null) {
        continue;
      }

      // Rebuild the block as `U * block`, where `U` is unimodular with first
      // row `x`; that puts the block minimum at position `kappa` while spanning
      // exactly the same lattice.  Handing a redundant `(d+1)`-row generating
      // set to LLL instead (as a first attempt did) lets LLL choose any basis
      // it likes and simply loses the vector we just found.
      const U = unimodularWithFirstRow(x);
      const block = current.slice(kappa, endIdx);
      const m = current[0]!.length;
      const newBlock: bigint[][] = U.map((urow) => {
        const out: bigint[] = new Array(m).fill(0n);
        for (let i = 0; i < urow.length; i++) {
          const ui = urow[i]!;
          if (ui === 0n) {
            continue;
          }
          const row = block[i]!;
          for (let j = 0; j < m; j++) {
            out[j] = out[j]! + ui * row[j]!;
          }
        }
        return out;
      });

      const next = [...current.slice(0, kappa), ...newBlock, ...current.slice(endIdx)];

      // LLL cannot swap row `kappa` back out: `b_kappa*` now realises the
      // minimum of the projected block, so
      // `B[kappa+1] >= (1 - mu^2) B[kappa] >= (delta - mu^2) B[kappa]`.
      // Rows before `kappa` are untouched, so their conditions still hold.
      const reduced = dropZeroRows(
        lllReduce(IntegerMatrixFromEntries(next), { delta, eta: 0.501 })
      );
      current = matrixRows(reduced);
      changed = true;
      if (current.length <= kappa) {
        break;
      }
    }
    if (!changed) {
      break;
    }
  }

  return current;
}

/**
 * Drop the zero rows of a matrix.
 *
 * SageMath does this with ``matrix([v for v in basis.LLL() if v])``
 * (free_module_integer.py:304).
 */
function dropZeroRows(M: IntegerMatrix): IntegerMatrix {
  const rows = matrixRows(M).filter((row) => row.some((v) => v !== 0n));
  if (rows.length === 0) {
    return new IntegerMatrix(0, M.ncols);
  }
  if (rows.length === M.nrows) {
    return M;
  }
  return new IntegerMatrix(rows.length, M.ncols, rows);
}

/**
 * Convert a target vector (array, FreeModuleElement, ...) to exact rationals.
 */
function toRatVector(target: FreeModuleElement | unknown[]): Rat[] {
  const raw: unknown[] = Array.isArray(target) ? target : (target as FreeModuleElement).list();
  return raw.map((x) => toRat(x));
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
  /** The delta parameter (0.25 < delta <= 1.0, default: 0.99) */
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
  const delta = options?.delta ?? 0.99;
  const eta = options?.eta ?? 0.501;

  // Validate parameters (matrix_integer_dense.pyx:3302 rejects delta <= 1/4,
  // delta > 1 and eta < 1/2).
  if (delta <= 0.25) {
    throw new ValueError(`delta must be > 0.25, got ${delta}`);
  }
  if (delta > 1.0) {
    throw new ValueError(`delta must be <= 1, got ${delta}`);
  }
  if (eta < 0.5) {
    throw new ValueError(`eta must be >= 0.5, got ${eta}`);
  }
  if (eta >= Math.sqrt(delta)) {
    throw new ValueError(
      `eta must be in [0.5, sqrt(delta)), got ${eta} (sqrt(delta) = ${Math.sqrt(delta)})`
    );
  }

  const n = basis.nrows;
  const m = basis.ncols;

  if (n === 0 || m === 0) {
    return new IntegerMatrix(n === 0 ? 0 : n, m);
  }

  // Copy basis to working array (bigint: exact arithmetic throughout)
  const rows: bigint[][] = [];
  for (let i = 0; i < n; i++) {
    const row: bigint[] = [];
    for (let j = 0; j < m; j++) {
      row.push(basis.get(i, j).value);
    }
    rows.push(row);
  }

  const deltaR = ratFromNumber(delta);
  const etaR = ratFromNumber(eta);

  let reduced = lllReduceIndependent(rows, deltaR, etaR);

  if (reduced === null) {
    // The rows are linearly dependent.  fpLLL/NTL return the zero rows first
    // followed by an LLL-reduced basis of the (lower rank) lattice
    // (matrix_integer_dense.pyx:3143 "Example with a nonzero kernel").  We get
    // an honest basis of the same lattice from the Hermite normal form, which
    // is exact, and reduce that.
    const hnf = hermite_normal_form(basis, 'default', false, false) as IntegerMatrix;
    const independent: bigint[][] = [];
    for (let i = 0; i < hnf.nrows; i++) {
      const row: bigint[] = [];
      let nonZero = false;
      for (let j = 0; j < m; j++) {
        const v = hnf.get(i, j).value;
        row.push(v);
        if (v !== 0n) nonZero = true;
      }
      if (nonZero) independent.push(row);
    }
    const sub = lllReduceIndependent(independent, deltaR, etaR);
    if (sub === null) {
      // Cannot happen: the nonzero rows of a Hermite normal form are
      // linearly independent.
      throw new ValueError('LLL: Hermite normal form returned dependent rows');
    }
    reduced = [];
    for (let i = 0; i < n - sub.length; i++) {
      reduced.push(new Array(m).fill(0n));
    }
    reduced.push(...sub);
  }

  // Build result matrix
  const result = new IntegerMatrix(n, m);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < m; j++) {
      result.set(i, j, reduced[i]![j]!);
    }
  }

  return result;
}

/**
 * LLL reduction of linearly independent rows, exactly, in place on a copy.
 *
 * This is Algorithm 2.6.3 of H. Cohen, *A Course in Computational Algebraic
 * Number Theory* (the same algorithm fpLLL implements, with the
 * Gram-Schmidt data kept as exact rationals rather than floating point
 * approximations, so that there is no precision to run out of).
 *
 * @returns the reduced rows, or `null` if the input rows are linearly dependent
 */
function lllReduceIndependent(rows: bigint[][], deltaR: Rat, etaR: Rat): bigint[][] | null {
  const n = rows.length;
  if (n === 0) {
    return [];
  }
  const m = rows[0]!.length;
  const b = rows.map((r) => [...r]);

  if (n === 1) {
    return dotBig(b[0]!, b[0]!) === 0n ? null : b;
  }

  const mu: Rat[][] = [];
  for (let i = 0; i < n; i++) {
    mu.push(new Array(n).fill(RAT_ZERO));
  }
  const B: Rat[] = new Array(n).fill(RAT_ZERO);

  B[0] = ratFromBigInt(dotBig(b[0]!, b[0]!));
  if (ratIsZero(B[0]!)) {
    return null;
  }

  // RED(k, l): size-reduce b_k against b_l
  const red = (k: number, l: number): void => {
    if (ratCmp(ratAbs(mu[k]![l]!), etaR) <= 0) {
      return;
    }
    const q = ratRoundHalfUp(mu[k]![l]!);
    if (q === 0n) {
      return;
    }
    for (let j = 0; j < m; j++) {
      b[k]![j] = b[k]![j]! - q * b[l]![j]!;
    }
    mu[k]![l] = ratSub(mu[k]![l]!, ratFromBigInt(q));
    for (let i = 0; i < l; i++) {
      mu[k]![i] = ratSub(mu[k]![i]!, ratMul(ratFromBigInt(q), mu[l]![i]!));
    }
  };

  // SWAP(k): exchange b_k and b_{k-1} and update the Gram-Schmidt data
  const swap = (k: number, kmax: number): void => {
    const tmp = b[k]!;
    b[k] = b[k - 1]!;
    b[k - 1] = tmp;
    for (let j = 0; j <= k - 2; j++) {
      const t = mu[k]![j]!;
      mu[k]![j] = mu[k - 1]![j]!;
      mu[k - 1]![j] = t;
    }
    const mu0 = mu[k]![k - 1]!;
    const BB = ratAdd(B[k]!, ratMul(ratMul(mu0, mu0), B[k - 1]!));
    mu[k]![k - 1] = ratDiv(ratMul(mu0, B[k - 1]!), BB);
    B[k] = ratDiv(ratMul(B[k - 1]!, B[k]!), BB);
    B[k - 1] = BB;
    for (let i = k + 1; i <= kmax; i++) {
      const t = mu[i]![k]!;
      mu[i]![k] = ratSub(mu[i]![k - 1]!, ratMul(mu0, t));
      mu[i]![k - 1] = ratAdd(t, ratMul(mu[k]![k - 1]!, mu[i]![k]!));
    }
  };

  let k = 1;
  let kmax = 0;
  // LLL provably terminates; the budget only guards against a pathological
  // non-terminating case for delta == 1, where the basis returned is still a
  // basis of the same lattice.
  let budget = 1000 * n * n + 100000;

  while (k < n && budget-- > 0) {
    if (k > kmax) {
      kmax = k;
      let Bk = ratFromBigInt(dotBig(b[k]!, b[k]!));
      for (let j = 0; j < k; j++) {
        let s = ratFromBigInt(dotBig(b[k]!, b[j]!));
        for (let i = 0; i < j; i++) {
          s = ratSub(s, ratMul(mu[j]![i]!, ratMul(mu[k]![i]!, B[i]!)));
        }
        if (ratIsZero(B[j]!)) {
          return null; // linearly dependent input
        }
        mu[k]![j] = ratDiv(s, B[j]!);
        Bk = ratSub(Bk, ratMul(ratMul(mu[k]![j]!, mu[k]![j]!), B[j]!));
      }
      B[k] = Bk;
      if (ratIsZero(Bk)) {
        return null; // linearly dependent input
      }
    }

    red(k, k - 1);

    const muKK1 = mu[k]![k - 1]!;
    const threshold = ratMul(ratSub(deltaR, ratMul(muKK1, muKK1)), B[k - 1]!);
    if (ratCmp(B[k]!, threshold) < 0) {
      swap(k, kmax);
      k = Math.max(1, k - 1);
    } else {
      for (let l = k - 2; l >= 0; l--) {
        red(k, l);
      }
      k++;
    }
  }

  return b;
}

/**
 * Check if a matrix is LLL-reduced with given parameters.
 *
 * Mirrors `Matrix_integer_dense.is_LLL_reduced` (matrix_integer_dense.pyx:3403)
 * with its ``'sage'`` algorithm, but exactly: `|mu_{i,j}| <= eta` for all
 * `i > j` and `|b_i*|^2 >= (delta - mu_{i,i-1}^2) |b_{i-1}*|^2`.
 *
 * Leading zero rows (as produced by LLL on a rank-deficient basis) are
 * ignored; genuinely dependent nonzero rows raise the same error SageMath's
 * ``'sage'`` algorithm raises.
 *
 * @param basis - The basis matrix to check
 * @param delta - The Lovász parameter (default: 0.99, as in SageMath)
 * @param eta - The size reduction parameter (default: 0.501)
 * @returns True if the basis is (delta, eta)-LLL-reduced
 */
export function isLLLReduced(
  basis: IntegerMatrix,
  delta: number = 0.99,
  eta: number = 0.501
): boolean {
  if (delta <= 0.25) {
    throw new ValueError(`delta must be > 0.25, got ${delta}`);
  }
  if (delta > 1.0) {
    throw new ValueError(`delta must be <= 1, got ${delta}`);
  }
  if (eta < 0.5) {
    throw new ValueError(`eta must be >= 0.5, got ${eta}`);
  }

  const rows: bigint[][] = [];
  for (let i = 0; i < basis.nrows; i++) {
    const row: bigint[] = [];
    let nonZero = false;
    for (let j = 0; j < basis.ncols; j++) {
      const v = basis.get(i, j).value;
      row.push(v);
      if (v !== 0n) nonZero = true;
    }
    if (nonZero) rows.push(row);
  }

  const n = rows.length;
  if (n <= 1) {
    return true;
  }

  const deltaR = ratFromNumber(delta);
  const etaR = ratFromNumber(eta);
  const { mu, B } = exactGramSchmidt(rows);

  for (let i = 0; i < n; i++) {
    if (ratIsZero(B[i]!)) {
      throw new ValueError('linearly dependent input for module version of Gram-Schmidt');
    }
  }

  // For any i > j, |mu_{i,j}| <= eta
  for (let i = 1; i < n; i++) {
    for (let j = 0; j < i; j++) {
      if (ratCmp(ratAbs(mu[i]![j]!), etaR) > 0) {
        return false;
      }
    }
  }

  // For any i < d, delta |b_i*|^2 <= |b_{i+1}* + mu_{i+1,i} b_i*|^2
  for (let i = 1; i < n; i++) {
    const muII1 = mu[i]![i - 1]!;
    const rhs = ratMul(ratSub(deltaR, ratMul(muII1, muII1)), B[i - 1]!);
    if (ratCmp(B[i]!, rhs) < 0) {
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
      // Dual form (crypto/lattice.py:294):
      //   B = [[I_n, -A'^T], [0, q*I_{m-n}]]
      // where A' is the random (m-n) x n block, so that |det(B)| = q^(m-n).
      B = [];

      // Generate the random matrix A' (shape: (m-n) x n)
      const A: bigint[][] = [];
      for (let i = 0; i < m - n; i++) {
        const row: bigint[] = [];
        for (let j = 0; j < n; j++) {
          row.push(randMod());
        }
        A.push(row);
      }

      // Top n rows: [I_n | -A'^T]
      for (let i = 0; i < n; i++) {
        const row: bigint[] = [];
        // Identity part (n columns)
        for (let j = 0; j < n; j++) {
          row.push(i === j ? 1n : 0n);
        }
        // -A'^T part (m-n columns): entry (i, j) of A'^T is A'[j][i]
        for (let j = 0; j < m - n; j++) {
          row.push(-A[j]![i]!);
        }
        B.push(row);
      }

      // Bottom m-n rows: [0 | q*I_{m-n}]
      for (let i = 0; i < m - n; i++) {
        const row: bigint[] = new Array(n).fill(0n);
        for (let j = 0; j < m - n; j++) {
          row.push(i === j ? q : 0n);
        }
        B.push(row);
      }

      // SageMath swaps rows i and m-i-1 for i in range(m//2), i.e. reverses
      // the row order.
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
 * Compute a basis of the kernel lattice `{x in Z^n : A x = 0 (mod q)}`.
 *
 * The lattice is obtained as the integer kernel of `[A^T | q I_m]`, projected
 * onto the first `n` coordinates: `(x, y)` with `x A^T + y (q I) = 0` is
 * exactly `A x = 0 (mod q)`.  The kernel is read off a Hermite normal form,
 * so the computation is exact.
 */
function kernelModQ(mat: IntegerMatrix, q: bigint): bigint[][] {
  const m = mat.nrows;
  const n = mat.ncols;

  if (q === 0n) {
    throw new ValueError('q must be nonzero');
  }
  const qAbs = q < 0n ? -q : q;

  // Rows of the auxiliary matrix, in Z^(m + n):
  //   [A^T_i | e_i]   for i < n      (A^T_i is column i of A)
  //   [q e_j | 0]     for j < m
  const aux: bigint[][] = [];
  for (let i = 0; i < n; i++) {
    const row: bigint[] = [];
    for (let j = 0; j < m; j++) {
      row.push(mat.get(j, i).value);
    }
    for (let j = 0; j < n; j++) {
      row.push(i === j ? 1n : 0n);
    }
    aux.push(row);
  }
  for (let j = 0; j < m; j++) {
    const row: bigint[] = new Array(m + n).fill(0n);
    row[j] = qAbs;
    aux.push(row);
  }

  const H = hermite_normal_form(
    IntegerMatrixFromEntries(aux),
    'default',
    false,
    true
  ) as IntegerMatrix;

  const basis: bigint[][] = [];
  for (let i = 0; i < H.nrows; i++) {
    let leadingZero = true;
    for (let j = 0; j < m; j++) {
      if (H.get(i, j).value !== 0n) {
        leadingZero = false;
        break;
      }
    }
    if (!leadingZero) continue;
    const row: bigint[] = [];
    let nonZero = false;
    for (let j = 0; j < n; j++) {
      const v = H.get(i, m + j).value;
      row.push(v);
      if (v !== 0n) nonZero = true;
    }
    if (nonZero) basis.push(row);
  }

  return basis;
}

/**
 * Create a q-ary lattice Lambda_q(A) = {x : Ax = 0 mod q}.
 *
 * The q-ary lattice associated with an m x n matrix A is the set of all
 * vectors x in Z^n such that Ax = 0 (mod q). This is relevant for LWE-based
 * cryptography.  It always contains q*Z^n, so it has rank n and volume
 * dividing q^n.
 *
 * @param A - A matrix (IntegerMatrix or 2D array)
 * @param q - The modulus
 * @returns The q-ary lattice
 *
 * @example
 * ```typescript
 * // {x in Z^2 : x_0 + 2 x_1 = 0 mod 7}, of volume 7
 * const L = qaryLattice([[1, 2], [2, 4]], 7n);
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

  return new FreeModuleIntegerLattice(kernelModQ(mat, qBig));
}

/**
 * Create the dual q-ary lattice Lambda_q^perp(A) = {x : A^T x = 0 mod q}.
 *
 * The dual q-ary lattice of an m x n matrix A is the set of all vectors
 * x in Z^m such that A^T x = 0 (mod q).
 *
 * @param A - A matrix (IntegerMatrix or 2D array)
 * @param q - The modulus
 * @returns The dual q-ary lattice
 *
 * @example
 * ```typescript
 * const L = qaryDualLattice([[1, 2], [3, 4]], 7n);
 * ```
 */
export function qaryDualLattice(
  A: IntegerMatrix | bigint[][] | number[][],
  q: number | bigint
): FreeModuleIntegerLattice {
  const qBig = typeof q === 'bigint' ? q : BigInt(q);

  let mat: IntegerMatrix;
  if (A instanceof IntegerMatrix) {
    mat = A;
  } else {
    mat = IntegerMatrixFromEntries(A as (bigint | number)[][]);
  }

  return new FreeModuleIntegerLattice(kernelModQ(mat.transpose(), qBig));
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
 * Root Hermite factor achieved by BKZ with block size beta.
 *
 * `delta(beta)` is *decreasing* in beta: the larger the block size, the
 * shorter the first basis vector.  For beta >= 40 the standard asymptotic
 * estimate `delta = ((pi beta)^(1/beta) beta / (2 pi e))^(1/(2(beta-1)))`
 * is used; below that it badly underestimates delta, so the experimental
 * values of Gama-Nguyen / Chen's BKZ simulator are interpolated instead.
 *
 * @param beta - block size (>= 2)
 * @returns the root Hermite factor delta_beta
 */
export function bkzRootHermiteFactor(beta: number): number {
  const small: [number, number][] = [
    [2, 1.0219],
    [5, 1.01862],
    [10, 1.01616],
    [15, 1.01485],
    [20, 1.0142],
    [25, 1.01342],
    [28, 1.01331],
    [40, 1.01295],
  ];

  if (beta <= 2) {
    return small[0]![1];
  }
  if (beta < 40) {
    for (let i = 1; i < small.length; i++) {
      const [b1, d1] = small[i]!;
      if (beta <= b1) {
        const [b0, d0] = small[i - 1]!;
        // linear interpolation between the two tabulated block sizes
        return d0 + ((beta - b0) * (d1 - d0)) / (b1 - b0);
      }
    }
  }
  if (beta === 40) {
    return small[small.length - 1]![1];
  }
  return (
    ((Math.PI * beta) ** (1 / beta) * (beta / (2 * Math.PI * Math.E))) ** (1 / (2 * (beta - 1)))
  );
}

/**
 * Estimate the BKZ block size needed to find vectors of a given length.
 *
 * A basis reduced with root Hermite factor delta has
 * `||b_1|| = delta^(n-1) * vol^(1/n)`, so the required delta is
 * `delta = (target / vol^(1/n))^(1/(n-1))`.  Since `delta_beta` decreases in
 * beta, the answer is the smallest beta with `delta_beta <= delta`.
 *
 * @param dimension - The lattice dimension n
 * @param volume - The lattice volume (det(L))
 * @param targetLength - The desired vector length
 * @returns The estimated block size (2 if LLL already suffices, n if even
 *   full HKZ reduction does not reach the target)
 *
 * @example
 * ```typescript
 * // A vector of length 100 in a dimension 100 lattice of volume 2^100 is
 * // far above the Gaussian heuristic, so LLL suffices.
 * estimateBKZBlockSize(100, 2n ** 100n, 100); // 2
 * ```
 */
export function estimateBKZBlockSize(
  dimension: number | bigint,
  volume: bigint | number,
  targetLength: number | bigint
): number {
  const n = typeof dimension === 'bigint' ? Number(dimension) : dimension;
  const target = typeof targetLength === 'bigint' ? Number(targetLength) : targetLength;

  if (n <= 1) {
    return 2; // Minimum meaningful block size
  }

  if (target <= 0) {
    throw new ValueError('targetLength must be positive');
  }

  // log(vol)/n computed without overflowing for cryptographic volumes.
  let logVol: number;
  if (typeof volume === 'bigint') {
    if (volume <= 0n) {
      throw new ValueError('volume must be positive');
    }
    logVol = bigLog(volume);
  } else {
    if (volume <= 0) {
      throw new ValueError('volume must be positive');
    }
    logVol = Math.log(volume);
  }

  // delta = (target / vol^(1/n))^(1/(n-1))
  const logDelta = (Math.log(target) - logVol / n) / (n - 1);
  const delta = Math.exp(logDelta);

  if (delta >= bkzRootHermiteFactor(2)) {
    // Target is achievable with plain LLL
    return 2;
  }

  for (let beta = 2; beta <= n; beta++) {
    if (bkzRootHermiteFactor(beta) <= delta) {
      return beta;
    }
  }

  // Even full HKZ reduction is not expected to reach the target.
  return n;
}

/**
 * Natural logarithm of a positive bigint, without overflowing to Infinity.
 */
function bigLog(n: bigint): number {
  if (n <= 0n) {
    throw new ValueError('logarithm of a nonpositive integer');
  }
  const bits = n.toString(2).length;
  if (bits <= 53) {
    return Math.log(Number(n));
  }
  const shift = BigInt(bits - 53);
  return Math.log(Number(n >> shift)) + Number(shift) * Math.LN2;
}

// ============================================================================
// Helper Functions
// ============================================================================

/**
 * Compute integer square root of a bigint.
 */
function bigintSqrt(n: bigint): bigint {
  return bigintSqrtFloor(n);
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
