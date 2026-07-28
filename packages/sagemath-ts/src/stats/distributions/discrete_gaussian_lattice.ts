/**
 * @module sage/stats/distributions/discrete_gaussian_lattice
 * @description Discrete Gaussian Distribution over Lattices
 *
 * This module implements the GPV (Gentry-Peikert-Vaikuntanathan) algorithm
 * for sampling from a discrete Gaussian distribution over a lattice, plus
 * Peikert's non-spherical sampler for a covariance matrix `Sigma`.
 *
 * The discrete Gaussian distribution D_{L,sigma,c} over a lattice L is defined by:
 *   P(x) = rho_{sigma,c}(x) / rho_{sigma,c}(L)
 *
 * where rho_{sigma,c}(x) = exp(-||x-c||^2 / (2*sigma^2)) is the Gaussian function.
 *
 * @see Reference: sage/stats/distributions/discrete_gaussian_lattice.py
 * @see [GPV2008] Gentry, Peikert, Vaikuntanathan, "Trapdoors for Hard Lattices...", 2008
 * @see [Pei2010] Peikert, "An Efficient and Parallel Gaussian Sampler for Lattices", 2010
 */

import {
  NotImplementedError,
  RuntimeError,
  TypeError as SageTypeError,
  ValueError,
} from '../../errors.js';
import { current_randstate } from '../../misc/randstate.js';
import { Rational } from '../../rings/rational.js';
import { type IntegerLike, toBigInt, toSafeNumber } from '../../types/coercion.js';
import {
  DiscreteGaussianDistributionIntegerSampler,
  type DiscreteGaussianOptions,
} from './discrete_gaussian_integer.js';

/**
 * Convert a basis/center entry to an exact rational.
 *
 * Sage's sampler takes a matrix over an exact ring (typically `ZZ` or `QQ`)
 * and does all lattice arithmetic exactly; only `sigma` lives in `RealField`.
 * We mirror that by keeping every basis and center entry as a {@link Rational}.
 */
function toRationalEntry(x: number | bigint | Rational): Rational {
  if (x instanceof Rational) {
    return x;
  }
  if (typeof x === 'bigint') {
    return new Rational(x, 1n);
  }
  if (typeof x !== 'number' || !Number.isFinite(x)) {
    throw new SageTypeError(`entries must be finite numbers, got ${x}`);
  }
  if (Number.isInteger(x)) {
    return new Rational(BigInt(x), 1n);
  }
  const str = x.toString();
  const eIndex = str.search(/[eE]/);
  if (eIndex < 0) {
    return Rational.from(str);
  }
  // Exponential notation, e.g. "1.5e-7": expand it exactly in base 10.
  const mantissa = Rational.from(str.slice(0, eIndex));
  const exponent = Number.parseInt(str.slice(eIndex + 1), 10);
  const power = new Rational(10n ** BigInt(Math.abs(exponent)), 1n);
  return exponent >= 0 ? mantissa.mul(power) : mantissa.div(power);
}

/**
 * Render a `RealField(53)` element the way Sage does: 15 significant digits.
 */
function rrRepr(x: number): string {
  if (x === 0) return '0.000000000000000';
  if (!Number.isFinite(x)) return String(x);
  return x.toPrecision(15);
}

/**
 * Render a matrix the way Sage does: one bracketed row per line, each column
 * padded to its widest entry and right-justified.
 */
function matrixRepr(rows: string[][]): string {
  if (rows.length === 0) return '[]';
  const ncols = rows[0]!.length;
  const widths: number[] = [];
  for (let j = 0; j < ncols; j++) {
    widths.push(Math.max(...rows.map((r) => r[j]!.length)));
  }
  return rows.map((r) => `[${r.map((e, j) => e.padStart(widths[j]!)).join(' ')}]`).join('\n');
}

/**
 * Result of Gram-Schmidt orthogonalization (exact).
 */
interface GramSchmidtResult {
  /** The orthogonalized basis B* (rows are b_1*, ..., b_n*) */
  bStar: Rational[][];
  /** The mu coefficients matrix where mu[i][j] = <b_i, b_j*> / <b_j*, b_j*> */
  mu: Rational[][];
  /** The squared norms of the orthogonal vectors |b_i*|^2 */
  bStarNormsSq: Rational[];
}

/**
 * Compute the exact Gram-Schmidt orthogonalization of a basis.
 *
 * This mirrors Sage's `B.gram_schmidt()` on an exact matrix
 * (`discrete_gaussian_lattice.py:587`), which returns rationals — not the
 * floating-point `gram_schmidt(..., orthonormal=True)` variant.
 *
 * @param basis - A matrix whose rows form the basis
 * @returns GramSchmidtResult containing orthogonal basis, mu matrix, and squared norms
 */
function gramSchmidt(basis: Rational[][]): GramSchmidtResult {
  const n = basis.length;
  if (n === 0) {
    return { bStar: [], mu: [], bStarNormsSq: [] };
  }
  const m = basis[0]!.length;

  const bStar: Rational[][] = [];
  const mu: Rational[][] = [];
  const bStarNormsSq: Rational[] = [];

  for (let i = 0; i < n; i++) {
    mu.push(new Array(n).fill(Rational.zero()));
    mu[i]![i] = Rational.one();

    // Start with b_i
    const bi = basis[i]!;
    const biStar = [...bi];

    // Subtract projections onto previous orthogonal vectors
    for (let j = 0; j < i; j++) {
      // mu[i][j] = <b_i, b_j*> / <b_j*, b_j*>
      let dotProduct = Rational.zero();
      for (let k = 0; k < m; k++) {
        dotProduct = dotProduct.add(bi[k]!.mul(bStar[j]![k]!));
      }
      if (bStarNormsSq[j]!.isZero()) {
        throw new ValueError('basis vectors must be linearly independent');
      }
      mu[i]![j] = dotProduct.div(bStarNormsSq[j]!);

      // b_i* = b_i* - mu[i][j] * b_j*
      for (let k = 0; k < m; k++) {
        biStar[k] = biStar[k]!.sub(mu[i]![j]!.mul(bStar[j]![k]!));
      }
    }

    bStar.push(biStar);

    // Compute |b_i*|^2
    let normSq = Rational.zero();
    for (let k = 0; k < m; k++) {
      normSq = normSq.add(biStar[k]!.mul(biStar[k]!));
    }
    bStarNormsSq.push(normSq);
  }

  return { bStar, mu, bStarNormsSq };
}

/**
 * Exact dot product of two rational vectors.
 */
function ratDot(a: Rational[], b: Rational[]): Rational {
  let sum = Rational.zero();
  for (let i = 0; i < a.length; i++) {
    sum = sum.add(a[i]!.mul(b[i]!));
  }
  return sum;
}

/** Exact matrix product. */
function ratMatMul(A: Rational[][], B: Rational[][]): Rational[][] {
  const n = A.length;
  const k = B.length;
  const m = B[0]!.length;
  const out: Rational[][] = [];
  for (let i = 0; i < n; i++) {
    const row: Rational[] = [];
    for (let j = 0; j < m; j++) {
      let s = Rational.zero();
      for (let t = 0; t < k; t++) {
        s = s.add(A[i]![t]!.mul(B[t]![j]!));
      }
      row.push(s);
    }
    out.push(row);
  }
  return out;
}

/** Exact transpose. */
function ratTranspose(A: Rational[][]): Rational[][] {
  const n = A.length;
  const m = A[0]!.length;
  const out: Rational[][] = [];
  for (let j = 0; j < m; j++) {
    const row: Rational[] = [];
    for (let i = 0; i < n; i++) row.push(A[i]![j]!);
    out.push(row);
  }
  return out;
}

/** Exact inverse by Gauss-Jordan elimination over QQ. */
function ratInverse(A: Rational[][]): Rational[][] {
  const n = A.length;
  if (n === 0 || A[0]!.length !== n) {
    throw new ValueError('matrix must be square to be inverted');
  }
  const M = A.map((r, i) => [
    ...r,
    ...Array.from({ length: n }, (_, j) => (i === j ? Rational.one() : Rational.zero())),
  ]);
  for (let col = 0; col < n; col++) {
    let pivot = -1;
    for (let row = col; row < n; row++) {
      if (!M[row]![col]!.isZero()) {
        pivot = row;
        break;
      }
    }
    if (pivot < 0) {
      throw new ValueError('matrix is not invertible');
    }
    if (pivot !== col) {
      const tmp = M[col]!;
      M[col] = M[pivot]!;
      M[pivot] = tmp;
    }
    const inv = M[col]![col]!.inv();
    M[col] = M[col]!.map((x) => x.mul(inv));
    for (let row = 0; row < n; row++) {
      if (row === col) continue;
      const factor = M[row]![col]!;
      if (factor.isZero()) continue;
      M[row] = M[row]!.map((x, j) => x.sub(factor.mul(M[col]![j]!)));
    }
  }
  return M.map((r) => r.slice(n));
}

/** Exact determinant by fraction-free-ish Gaussian elimination over QQ. */
function ratDet(A: Rational[][]): Rational {
  const n = A.length;
  if (n === 0) return Rational.one();
  const M = A.map((r) => [...r]);
  let det = Rational.one();
  for (let col = 0; col < n; col++) {
    let pivot = -1;
    for (let row = col; row < n; row++) {
      if (!M[row]![col]!.isZero()) {
        pivot = row;
        break;
      }
    }
    if (pivot < 0) return Rational.zero();
    if (pivot !== col) {
      const tmp = M[col]!;
      M[col] = M[pivot]!;
      M[pivot] = tmp;
      det = det.neg();
    }
    det = det.mul(M[col]![col]!);
    const inv = M[col]![col]!.inv();
    for (let row = col + 1; row < n; row++) {
      const factor = M[row]![col]!.mul(inv);
      if (factor.isZero()) continue;
      M[row] = M[row]!.map((x, j) => x.sub(factor.mul(M[col]![j]!)));
    }
  }
  return det;
}

/**
 * Solve `x * A = b` exactly (Sage's `A.solve_left(b)`), for a square
 * invertible `A`.
 */
function ratSolveLeft(A: Rational[][], b: Rational[]): Rational[] {
  const Ainv = ratInverse(A);
  return Ainv[0]!.map((_, j) =>
    b.reduce((acc, bi, i) => acc.add(bi.mul(Ainv[i]![j]!)), Rational.zero())
  );
}

/** Double-precision matrix product. */
function rrMatMul(A: number[][], B: number[][]): number[][] {
  const n = A.length;
  const k = B.length;
  const m = B[0]!.length;
  const out: number[][] = [];
  for (let i = 0; i < n; i++) {
    const row: number[] = [];
    for (let j = 0; j < m; j++) {
      let s = 0;
      for (let t = 0; t < k; t++) s += A[i]![t]! * B[t]![j]!;
      row.push(s);
    }
    out.push(row);
  }
  return out;
}

/** Double-precision matrix inverse via Gauss-Jordan with partial pivoting. */
function rrInverse(A: number[][]): number[][] {
  const n = A.length;
  const M = A.map((r, i) => [...r, ...Array.from({ length: n }, (_, j) => (i === j ? 1 : 0))]);
  for (let col = 0; col < n; col++) {
    let pivot = col;
    for (let row = col + 1; row < n; row++) {
      if (Math.abs(M[row]![col]!) > Math.abs(M[pivot]![col]!)) pivot = row;
    }
    if (M[pivot]![col] === 0) {
      throw new ValueError('matrix is not invertible');
    }
    if (pivot !== col) {
      const tmp = M[col]!;
      M[col] = M[pivot]!;
      M[pivot] = tmp;
    }
    const d = M[col]![col]!;
    M[col] = M[col]!.map((x) => x / d);
    for (let row = 0; row < n; row++) {
      if (row === col) continue;
      const f = M[row]![col]!;
      if (f === 0) continue;
      M[row] = M[row]!.map((x, j) => x - f * M[col]![j]!);
    }
  }
  return M.map((r) => r.slice(n));
}

/**
 * Cholesky decomposition `A = L L^T` with `L` lower triangular.
 *
 * Mirrors Sage's `Matrix.cholesky()` over `RealField`, which raises
 * `ValueError` when the matrix is not positive definite.
 */
function rrCholesky(A: number[][]): number[][] {
  const n = A.length;
  const L: number[][] = Array.from({ length: n }, () => new Array<number>(n).fill(0));
  for (let i = 0; i < n; i++) {
    for (let j = 0; j <= i; j++) {
      let s = A[i]![j]!;
      for (let k = 0; k < j; k++) s -= L[i]![k]! * L[j]![k]!;
      if (i === j) {
        if (s <= 0) {
          throw new ValueError('matrix is not positive definite');
        }
        L[i]![j] = Math.sqrt(s);
      } else {
        L[i]![j] = s / L[j]![j]!;
      }
    }
  }
  return L;
}

/**
 * Sage's `Matrix.is_positive_definite()` for a real symmetric matrix.
 */
function rrIsPositiveDefinite(A: number[][]): boolean {
  const n = A.length;
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      if (Math.abs(A[i]![j]! - A[j]![i]!) > 0) return false;
    }
  }
  try {
    rrCholesky(A);
    return true;
  } catch {
    return false;
  }
}

/**
 * LLL reduction of an integral row basis, mirroring `B.LLL()`.
 *
 * @see Deviation: Sage delegates to fpLLL with `delta = 0.99, eta = 0.501`;
 * this is the classical exact-rational LLL with the same `delta`.  The two
 * agree on the reduced *lattice* but may return different reduced bases for
 * inputs with several near-optimal reductions.  It is used only by
 * {@link DiscreteGaussianDistributionLatticeSampler._normalisation_factor_zz}'s
 * non-spherical branch, which needs it to pick an enumeration centre.
 */
function lllReduce(basis: Rational[][]): Rational[][] {
  const delta = new Rational(99n, 100n);
  const B = basis.map((r) => [...r]);
  const n = B.length;
  if (n <= 1) return B;
  let gs = gramSchmidt(B);
  let k = 1;
  let guard = 0;
  const maxIter = 1000 * n * n + 1000;
  while (k < n && guard++ < maxIter) {
    for (let j = k - 1; j >= 0; j--) {
      const mu = gs.mu[k]![j]!;
      const q = mu.round();
      if (q !== 0n) {
        const qr = new Rational(q, 1n);
        B[k] = B[k]!.map((x, t) => x.sub(qr.mul(B[j]![t]!)));
        gs = gramSchmidt(B);
      }
    }
    const mukk1 = gs.mu[k]![k - 1]!;
    const lhs = gs.bStarNormsSq[k]!;
    const rhs = delta.sub(mukk1.mul(mukk1)).mul(gs.bStarNormsSq[k - 1]!);
    if (lhs.cmp(rhs) >= 0) {
      k++;
    } else {
      const tmp = B[k]!;
      B[k] = B[k - 1]!;
      B[k - 1] = tmp;
      gs = gramSchmidt(B);
      k = Math.max(k - 1, 1);
    }
  }
  return B;
}

/**
 * Theta series of a positive definite integral quadratic form.
 *
 * `qfrep(Q, B)` returns, for `m = 1..floor(B)`, *half* the number of integer
 * vectors `x != 0` with `x Q x^T = m` — that is one representative of each
 * `{x, -x}` pair.  This is PARI's `qfrep(Q, B, 0)`, which Sage calls at
 * `discrete_gaussian_lattice.py:340` and `:351`.
 *
 * @see Deviation: Sage delegates to PARI (`Q.__pari__().qfrep(...)`).
 * `parigp-ts` does not expose `qfrep` yet, so this is a local Fincke-Pohst
 * enumeration over the Cholesky decomposition of `Q`.  The counts are exact
 * integers, so the result is identical to PARI's.
 */
export function qfrep(Q: number[][], bound: number): bigint[] {
  const b = Math.floor(bound);
  const counts = new Array<bigint>(Math.max(b, 0)).fill(0n);
  if (b <= 0) return counts;
  const n = Q.length;
  // Q = R^T R with R upper triangular: use the Cholesky factor of Q.
  const L = rrCholesky(Q);
  // Q(x) = sum_i ( sum_{j>=i} u[i][j] x_j )^2 with u = L^T scaled; use the
  // standard "quadratic completion" form q[i][j].
  // Write Q(x) = sum_i q_ii (x_i + sum_{j>i} q_ij x_j)^2 .
  const q: number[][] = Array.from({ length: n }, () => new Array<number>(n).fill(0));
  for (let i = 0; i < n; i++) {
    q[i]![i] = L[i]![i]! * L[i]![i]!;
    for (let j = 0; j < i; j++) {
      q[j]![i] = L[i]![j]! / L[j]![j]!;
    }
  }
  const x = new Array<number>(n).fill(0);
  // Enumerate from the last coordinate down.
  const rec = (i: number, remaining: number): void => {
    if (i < 0) {
      let s = 0;
      for (let a = 0; a < n; a++) {
        for (let c = 0; c < n; c++) s += x[a]! * Q[a]![c]! * x[c]!;
      }
      const m = Math.round(s);
      if (m >= 1 && m <= b) counts[m - 1] = counts[m - 1]! + 1n;
      return;
    }
    let centre = 0;
    for (let j = i + 1; j < n; j++) centre += q[i]![j]! * x[j]!;
    const width = Math.sqrt(Math.max(remaining, 0) / q[i]![i]!);
    const lo = Math.ceil(-centre - width - 1e-9);
    const hi = Math.floor(-centre + width + 1e-9);
    for (let v = lo; v <= hi; v++) {
      x[i] = v;
      const t = v + centre;
      rec(i - 1, remaining - q[i]![i]! * t * t);
    }
    x[i] = 0;
  };
  rec(n - 1, b);
  // The enumeration counted x and -x, and also the zero vector.
  return counts.map((c) => c / 2n);
}

/**
 * Options for constructing a discrete Gaussian lattice sampler.
 */
export interface DiscreteGaussianLatticeOptions {
  /**
   * Gaussian parameter, one of
   * - a real number `sigma > 0` (spherical),
   * - a positive definite matrix `Sigma` (non-spherical), or
   * - any matrix-like `S`, equivalent to `Sigma = S S^T`, when
   *   {@link sigma_basis} is set.
   */
  sigma: number | number[][] | Rational[][];

  /**
   * Center of the distribution (default: origin).
   * Should be a vector of the same dimension as the lattice.
   */
  c?: number[] | bigint[] | Rational[];

  /**
   * Tail cutoff parameter tau >= 1 (default: 6).
   *
   * @see Deviation: extra knob, forwarded to the per-coordinate integer
   * samplers.  Sage's lattice sampler always leaves `tau` at its default.
   */
  tau?: IntegerLike;

  /**
   * Rounding parameter `r` as defined in [Pei2010]_; ignored for a spherical
   * Gaussian parameter.  If not provided, set to the maximal possible value
   * such that `Sigma - r^2 B B^T` is positive definite.
   */
  r?: number;

  /** Bit precision `>= 53`. */
  precision?: number;

  /**
   * When set, `sigma` is treated as a (row) basis, i.e. the covariance matrix
   * is `Sigma = S S^T`.
   */
  sigma_basis?: boolean;
}

/**
 * Discrete Gaussian Distribution Sampler over a Lattice.
 *
 * Implements the GPV algorithm for a scalar `sigma` and Peikert's sampler for
 * a covariance matrix `Sigma`.
 *
 * @example
 * ```typescript
 * const D = new DiscreteGaussianDistributionLatticeSampler([[1, 0], [0, 1]], { sigma: 3 });
 * const v = D.sample(); // e.g., [2n, -1n]
 * ```
 *
 * @see Reference: sage/stats/distributions/discrete_gaussian_lattice.py
 */
export class DiscreteGaussianDistributionLatticeSampler {
  /**
   * The lattice basis (rows are basis vectors), as floating point.
   *
   * This is a lossy view kept for convenience; all sampling arithmetic uses
   * {@link basisExact}.
   */
  public readonly basis: number[][];

  /**
   * The lattice basis (rows are basis vectors), exactly.  Sage calls this `B`.
   */
  public readonly basisExact: Rational[][];

  /** Sage's `self.B`: the exact (row) basis. */
  public readonly B: Rational[][];

  /** Sage's `self.Q = B * B.T` (`discrete_gaussian_lattice.py:586`). */
  public readonly Q: Rational[][];

  /** Sage's `self.n = B.ncols()`. */
  public readonly n: number;

  /**
   * Whether every basis entry is an integer, i.e. the lattice sits inside Z^d.
   *
   * When false, {@link sample} cannot return integer coordinates and
   * {@link sampleExact} must be used instead.
   */
  public readonly isIntegral: boolean;

  /**
   * Whether the Gaussian parameter is a scalar (`true`) or a covariance
   * matrix (`false`).  Sage: `self.is_spherical`.
   */
  public readonly is_spherical: boolean;

  /**
   * Tail cutoff parameter.
   */
  public readonly tau: number;

  /**
   * Dimension of the lattice (number of basis vectors).
   */
  public readonly rank: number;

  /**
   * Dimension of the ambient space.
   */
  public readonly degree: number;

  /** Working precision, Sage's `compute_precision(precision, sigma)`. */
  public readonly precision: number;

  /**
   * Offline samples of `B^{-1} D_1`, used by the non-spherical sampler
   * (`discrete_gaussian_lattice.py:874-892`).
   */
  public offline_samples: number[][] = [];

  /** Peikert's rounding parameter (non-spherical only). */
  public r: number | null;

  /** Gaussian parameter: a scalar, or a covariance matrix in `RealField`. */
  private _sigma: number | number[][];

  /** Inverse of the covariance matrix (non-spherical only). */
  private sigma_inv: number[][] | null = null;

  /** Center. */
  private _c: Rational[];

  /** `c * B^{-1}` (non-spherical only), kept exact. */
  private _c_mul_B_inv: Rational[] | null = null;

  /** Exact inverse of the basis (non-spherical only). */
  private B_inv: Rational[][] | null = null;

  /** `Sigma_2 = Sigma - r^2 Q` Cholesky factor, transposed. */
  private B2: number[][] | null = null;

  /** `B2 * B^{-1}` (non-spherical only). */
  private B2_B_inv: number[][] | null = null;

  /** Whether `c` is in the lattice and `B^* == 1`. */
  private _c_in_lattice_and_lattice_trivial = false;

  /** Per-coordinate samplers used by `_call_simple`. */
  private D: DiscreteGaussianDistributionIntegerSampler[] | null = null;

  /**
   * Gram-Schmidt orthogonalization data (Sage's `self._G`, plus norms).
   */
  private readonly gs: GramSchmidtResult;

  /** Cached `_maximal_r()` (Sage caches it with `@cached_method`). */
  private _maximal_r_cache: number | null = null;

  /**
   * Construct a new discrete Gaussian lattice sampler.
   *
   * @param basis - The lattice basis (rows are basis vectors)
   * @param options - Configuration options
   */
  constructor(
    basis: number[][] | bigint[][] | Rational[][],
    options: DiscreteGaussianLatticeOptions
  ) {
    // Validate basis
    if (!Array.isArray(basis) || basis.length === 0) {
      throw new ValueError('basis must be a non-empty array of vectors');
    }

    // Keep the basis exactly; the number[][] view is only for reporting.
    this.basisExact = (basis as (number | bigint | Rational)[][]).map((row) =>
      row.map((x) => toRationalEntry(x))
    );
    this.B = this.basisExact;
    this.basis = this.basisExact.map((row) => row.map((x) => x.toNumber()));
    this.isIntegral = this.basisExact.every((row) => row.every((x) => x.isInteger()));
    this.rank = this.basisExact.length;
    this.degree = this.basisExact[0]!.length;
    this.n = this.degree;

    // Validate all rows have the same length
    for (const row of this.basisExact) {
      if (row.length !== this.degree) {
        throw new ValueError('all basis vectors must have the same dimension');
      }
    }

    // --- sigma: scalar or covariance matrix (discrete_gaussian_lattice.py:553-571)
    if (options.sigma === undefined || options.sigma === null) {
      throw new SageTypeError('sigma is required');
    }
    this.precision = DiscreteGaussianDistributionLatticeSampler.compute_precision(
      options.precision,
      options.sigma
    );
    if (typeof options.sigma === 'number') {
      if (!Number.isFinite(options.sigma)) {
        throw new SageTypeError(`sigma must be a finite number, got ${options.sigma}`);
      }
      if (options.sigma <= 0) {
        throw new ValueError(`sigma must be > 0, got ${options.sigma}`);
      }
      this._sigma = options.sigma;
      this.is_spherical = true;
    } else {
      let S = (options.sigma as (number | bigint | Rational)[][]).map((row) =>
        row.map((x) => toRationalEntry(x).toNumber())
      );
      // A scaled identity is handled as a scalar (py:564-565).
      const s00 = S[0]![0]!;
      const isScalarMatrix = S.every((row, i) => row.every((v, j) => v === (i === j ? s00 : 0)));
      if (isScalarMatrix) {
        this._sigma = s00;
        this.is_spherical = true;
      } else {
        if (options.sigma_basis) {
          S = rrMatMul(
            S,
            S[0]!.map((_, j) => S.map((r) => r[j]!))
          );
        }
        if (!rrIsPositiveDefinite(S)) {
          // Sage raises RuntimeError here (discrete_gaussian_lattice.py:570).
          throw new RuntimeError(
            `Sigma(=${matrixRepr(S.map((r) => r.map(rrRepr)))}) is not positive definite`
          );
        }
        this._sigma = S;
        this.is_spherical = false;
      }
    }

    // Validate and store tau (default: 6)
    const tauValue = options.tau !== undefined ? toSafeNumber(toBigInt(options.tau)) : 6;
    if (tauValue < 1) {
      throw new ValueError(`tau must be >= 1, got ${tauValue}`);
    }
    this.tau = tauValue;

    // Q = B * B^T and the exact Gram-Schmidt orthogonalization.
    this.Q = ratMatMul(this.basisExact, ratTranspose(this.basisExact));
    this.gs = gramSchmidt(this.basisExact);

    this.r = options.r ?? null;
    this._c = new Array(this.degree).fill(Rational.zero());
    this.set_c(options.c ?? null);
  }

  /**
   * Compute precision to use.
   *
   * Reference: `discrete_gaussian_lattice.py:155-188`.
   */
  static compute_precision(precision: number | undefined | null, _sigma: unknown): number {
    if (precision === undefined || precision === null) {
      // Our reals are doubles, so `sigma.precision()` is always 53.
      return 53;
    }
    return Math.max(53, precision);
  }

  /**
   * Gaussian parameter `sigma`.
   *
   * Reference: `discrete_gaussian_lattice.py:723-736`.
   */
  sigma(): number | number[][] {
    return this._sigma;
  }

  /**
   * Center `c`.
   *
   * Reference: `discrete_gaussian_lattice.py:738-756`.
   */
  c(): Rational[] {
    return this._c;
  }

  /**
   * Center `c` as a floating-point vector (convenience view).
   */
  cNumeric(): number[] {
    return this._c.map((x) => x.toNumber());
  }

  /**
   * Modify the center `c`, recomputing the cached data.
   *
   * Reference: `discrete_gaussian_lattice.py:758-792`.
   */
  set_c(c: number[] | bigint[] | Rational[] | null): void {
    if (c === null || c === undefined) {
      this._c = new Array(this.degree).fill(Rational.zero());
    } else {
      if (!Array.isArray(c) || c.length !== this.degree) {
        throw new ValueError(`c must be a vector of dimension ${this.degree}`);
      }
      this._c = (c as (number | bigint | Rational)[]).map((x) => toRationalEntry(x));
    }
    this._precompute_data();
  }

  /**
   * Precompute basis data.  Do not call directly.
   *
   * Reference: `discrete_gaussian_lattice.py:597-660`.
   */
  private _precompute_data(): void {
    this.D = null;
    this._c_in_lattice_and_lattice_trivial = false;
    if (this.is_spherical) {
      const sigma = this._sigma as number;
      const cIsZero = this._c.every((x) => x.isZero());
      const gIsIdentity = this._gramSchmidtIsIdentity();
      if (cIsZero && gIsIdentity) {
        this._c_in_lattice_and_lattice_trivial = true;
        const D = new DiscreteGaussianDistributionIntegerSampler({ sigma, tau: this.tau });
        this.D = new Array(this.rank).fill(D);
      } else if (gIsIdentity) {
        // w = B.solve_left(c); trivial iff w is integral.
        let w: Rational[] | null = null;
        try {
          w = ratSolveLeft(this.basisExact, this._c);
        } catch {
          w = null;
        }
        if (w?.every((x) => x.isInteger())) {
          this._c_in_lattice_and_lattice_trivial = true;
          const D: DiscreteGaussianDistributionIntegerSampler[] = [];
          for (let i = 0; i < this.rank; i++) {
            const sigma_ = sigma / Math.sqrt(this.gs.bStarNormsSq[i]!.toNumber());
            D.push(
              new DiscreteGaussianDistributionIntegerSampler({ sigma: sigma_, tau: this.tau })
            );
          }
          this.D = D;
        }
      }
      return;
    }

    // Non-spherical: [Pei2010] (discrete_gaussian_lattice.py:636-660).
    const Sigma = this._sigma as number[][];
    this.offline_samples = [];
    this.B_inv = ratInverse(this.basisExact);
    this.sigma_inv = rrInverse(Sigma);
    this._c_mul_B_inv = this._c.map((_, j) =>
      this._c.reduce((acc, ci, i) => acc.add(ci.mul(this.B_inv![i]![j]!)), Rational.zero())
    );

    if (this.r === null) {
      this.r = this._maximal_r() * 0.9999;
    }
    const r = this.r;
    const Qd = this.Q.map((row) => row.map((x) => x.toNumber()));
    const Sigma2 = Sigma.map((row, i) => row.map((v, j) => v - r * r * Qd[i]![j]!));
    let L: number[][];
    try {
      L = rrCholesky(Sigma2);
    } catch {
      throw new ValueError(
        `Σ₂ is not positive definite. Is your r(=${rrRepr(r)}) too large? It should be at most ${rrRepr(this._maximal_r())}`
      );
    }
    // Sage: self.B2 = Sigma2.cholesky().T
    this.B2 = L[0]!.map((_, j) => L.map((row) => row[j]!));
    this.B2_B_inv = rrMatMul(
      this.B2,
      this.B_inv.map((row) => row.map((x) => x.toNumber()))
    );
  }

  /** Whether the exact Gram-Schmidt matrix `_G` equals the identity. */
  private _gramSchmidtIsIdentity(): boolean {
    if (this.rank !== this.degree) return false;
    for (let i = 0; i < this.rank; i++) {
      for (let j = 0; j < this.degree; j++) {
        const want = i === j ? 1n : 0n;
        const v = this.gs.bStar[i]![j]!;
        if (!v.isInteger() || v.numerator !== want) return false;
      }
    }
    return true;
  }

  /**
   * The largest `r > 0` such that `Sigma - r^2 B B^T` is positive definite,
   * found by power iteration.
   *
   * Reference: `discrete_gaussian_lattice.py:357-389`.
   */
  _maximal_r(): number {
    if (!this.is_spherical) {
      if (this._maximal_r_cache !== null) return this._maximal_r_cache;
      const Sigma = this._sigma as number[][];
      const Qd = this.Q.map((row) => row.map((x) => x.toNumber()));
      // Q = self.Q / self._sigma  (i.e. Q * Sigma^{-1})
      const Qh = rrMatMul(Qd, rrInverse(Sigma));
      let v = [...Qh[0]!];
      const apply = (m: number[][], x: number[]): number[] =>
        m.map((row) => row.reduce((acc, e, j) => acc + e * x[j]!, 0));
      const norm = (x: number[]) => Math.sqrt(x.reduce((a, e) => a + e * e, 0));
      for (let cnt = 0; cnt < 10000; cnt++) {
        const w = apply(Qh, v);
        const nw = norm(w);
        const nv = w.map((e) => e / nw);
        const diff = Math.sqrt(nv.reduce((a, e, i) => a + (e - v[i]!) ** 2, 0));
        if (diff < 1e-12) {
          v = nv;
          break;
        }
        v = nv;
      }
      const res = Math.sqrt(v[0]! / apply(Qh, v)[0]!);
      this._maximal_r_cache = res;
      return res;
    }
    throw new ValueError('_maximal_r is only defined for non-spherical samplers');
  }

  /**
   * Compute the Gaussian `rho_{Lambda, c, Sigma}`.
   *
   * Reference: `discrete_gaussian_lattice.py:698-721`.
   */
  f(x: number[] | bigint[] | Rational[]): number {
    const xr = (x as (number | bigint | Rational)[]).map((v) => toRationalEntry(v));
    if (xr.length !== this.n) {
      throw new ValueError(`expected a vector of dimension ${this.n}`);
    }
    const d = xr.map((v, i) => v.sub(this._c[i]!));
    if (this.is_spherical) {
      const sigma = this._sigma as number;
      const normSq = ratDot(d, d).toNumber();
      return Math.exp(-normSq / (2 * sigma * sigma));
    }
    const dn = d.map((v) => v.toNumber());
    const S = this.sigma_inv!;
    let q = 0;
    for (let i = 0; i < dn.length; i++) {
      for (let j = 0; j < dn.length; j++) q += dn[i]! * S[i]![j]! * dn[j]!;
    }
    return Math.exp(-q / 2);
  }

  /**
   * Approximate `sum_{x in B} exp(-|x|^2/(2 sigma^2))`, the normalisation
   * factor, via Poisson summation.
   *
   * Reference: `discrete_gaussian_lattice.py:190-355`.
   *
   * @see Deviation: the sum is evaluated in double precision.  Sage evaluates
   * it in `RealField(prec)`, so a `prec > 53` request cannot be honoured
   * digit-for-digit — the value returned here agrees with Sage's to
   * approximately 15 significant digits.
   */
  _normalisation_factor_zz(tau?: number, prec?: number): number {
    void prec;
    if (!this.is_spherical) {
      // discrete_gaussian_lattice.py:295-313
      console.warn(
        'Note: `_normalisation_factor_zz` has not been properly implemented for non-spherical distributions.'
      );
      const basis = lllReduce(this.basisExact);
      const solved = ratSolveLeft(basis, this._c);
      const base = solved.map((v) => new Rational(v.round(), 1n));
      let BOUND = Math.max(1, Math.floor((Math.ceil(10 ** (4 / this.n)) - 1) / 2));
      BOUND = Math.min(BOUND, 10);
      let total = 0;
      const u = new Array<number>(this.n).fill(-BOUND);
      const rec = (i: number): void => {
        if (i === this.n) {
          const coords = u.map((v, j) => new Rational(BigInt(v), 1n).add(base[j]!));
          const pt = coords.map((_, j) =>
            coords.reduce((acc, ck, k) => acc.add(ck.mul(this.basisExact[k]![j]!)), Rational.zero())
          );
          total += this.f(pt);
          return;
        }
        for (let v = -BOUND; v <= BOUND; v++) {
          u[i] = v;
          rec(i + 1);
        }
      };
      rec(0);
      return total;
    }

    if (this.rank !== this.degree) {
      throw new NotImplementedError('basis must be a square matrix');
    }
    if (!this.isIntegral) {
      throw new NotImplementedError('lattice must be integral');
    }
    if (!this._c_in_lattice_and_lattice_trivial) {
      throw new NotImplementedError('center must be at zero and basis must be trivial');
    }

    const sigma = this._sigma as number;
    // f_or_hat (py:285-293)
    const f_or_hat = (x: number): number =>
      sigma > 1
        ? Math.exp(-Math.PI * Math.PI * (2 * sigma * sigma) * x)
        : Math.exp(-x / (2 * sigma * sigma));

    let det = 1;
    let norm_factor = 1;
    if (sigma > 1) {
      det = ratDet(this.basisExact).toNumber();
      norm_factor = (sigma * Math.sqrt(2 * Math.PI)) ** this.n / det;
    }

    const Qd = this.Q.map((row) => row.map((x) => x.toNumber()));
    if (tau !== undefined && tau !== null) {
      const freq = qfrep(Qd, tau * sigma);
      let res = 1;
      for (let x = 0; x < freq.length; x++) {
        res += 2 * Number(freq[x]!) * f_or_hat((x + 1) / det ** this.n);
      }
      return norm_factor * res;
    }

    let res = 1;
    let bound = 0;
    for (;;) {
      bound += 1;
      const cnt = Number(qfrep(Qd, bound)[bound - 1]!);
      const inc = 2 * cnt * f_or_hat(bound / det ** this.n);
      if (cnt > 0 && res === res + inc) {
        return norm_factor * res;
      }
      res += inc;
      if (bound > 100000) {
        throw new NotImplementedError(
          'SAGE_NOT_IMPLEMENTED: _normalisation_factor_zz did not converge; pass tau'
        );
      }
    }
  }

  /**
   * Randomly round to the lattice coset `Z + v` with Gaussian parameter `r`.
   *
   * Reference: `discrete_gaussian_lattice.py:391-407`.
   */
  _randomise(v: number[]): bigint[] {
    return v.map((vi) =>
      new DiscreteGaussianDistributionIntegerSampler({ sigma: this.r!, c: vi }).sample()
    );
  }

  /**
   * Precompute samples from `B^{-1} D_1` for {@link _call_non_spherical}.
   *
   * Reference: `discrete_gaussian_lattice.py:874-892`.
   */
  add_offline_samples(cnt = 1): void {
    if (this.is_spherical) {
      throw new ValueError('offline samples are only used by the non-spherical sampler');
    }
    const p = current_randstate().python_random();
    for (let i = 0; i < cnt; i++) {
      const coord = Array.from({ length: this.n }, () => p.normalvariate(0, 1));
      const M = this.B2_B_inv!;
      // vector(coord) * B2_B_inv
      this.offline_samples.push(
        M[0]!.map((_, j) => coord.reduce((acc, ci, i2) => acc + ci * M[i2]![j]!, 0))
      );
    }
  }

  /**
   * Sample from the discrete Gaussian distribution over the lattice, exactly.
   *
   * Dispatches exactly as Sage's `__call__`
   * (`discrete_gaussian_lattice.py:662-696`).
   */
  sampleExact(): Rational[] {
    if (!this.is_spherical) {
      return this._call_non_spherical();
    }
    if (this._c_in_lattice_and_lattice_trivial) {
      return this._call_simple();
    }
    return this._call();
  }

  /**
   * `_call_simple`: `c` is in the lattice and `B^* == 1`.
   *
   * Reference: `discrete_gaussian_lattice.py:822-840`.
   */
  _call_simple(): Rational[] {
    const w = this.D!.map((d) => new Rational(d.sample(), 1n));
    return this._c.map((ci, j) =>
      w.reduce((acc, wi, i) => acc.add(wi.mul(this.basisExact[i]![j]!)), ci)
    );
  }

  /**
   * `_call`: the general GPV loop.
   *
   * Reference: `discrete_gaussian_lattice.py:842-872`.
   */
  _call(): Rational[] {
    let c = [...this._c];
    let v: Rational[] = new Array(this.degree).fill(Rational.zero());
    const sigma = this._sigma as number;

    for (let i = this.rank - 1; i >= 0; i--) {
      const bStarI = this.gs.bStar[i]!;
      const bStarNormSq = this.gs.bStarNormsSq[i]!;
      const ci = ratDot(c, bStarI).div(bStarNormSq);
      const sigmaI = sigma / Math.sqrt(bStarNormSq.toNumber());
      const Di = new DiscreteGaussianDistributionIntegerSampler({
        sigma: sigmaI,
        c: ci.toNumber(),
        tau: this.tau,
        // Sage passes this explicitly (py:869).
        algorithm: 'uniform+online',
      });
      const z = Di.sample();

      const zRat = new Rational(z, 1n);
      const bi = this.basisExact[i]!;
      c = c.map((x, j) => x.sub(zRat.mul(bi[j]!)));
      v = v.map((x, j) => x.add(zRat.mul(bi[j]!)));
    }

    return v;
  }

  /**
   * `_call_non_spherical`: Peikert's sampler.
   *
   * Reference: `discrete_gaussian_lattice.py:894-915`.
   */
  _call_non_spherical(): Rational[] {
    if (this.offline_samples.length === 0) {
      this.add_offline_samples();
    }
    const off = this.offline_samples.pop()!;
    const vec = this._c_mul_B_inv!.map((x, i) => x.toNumber() - off[i]!);
    const z = this._randomise(vec);
    return this.basisExact[0]!.map((_, j) =>
      z.reduce(
        (acc, zi, i) => acc.add(new Rational(zi, 1n).mul(this.basisExact[i]![j]!)),
        Rational.zero()
      )
    );
  }

  /**
   * Sample from the discrete Gaussian distribution over the lattice.
   *
   * @returns A sample as a bigint vector
   * @throws ValueError if the sample is not integral; use {@link sampleExact}
   */
  sample(): bigint[] {
    if (!this.isIntegral) {
      throw new ValueError(
        'lattice basis is not integral; use sampleExact() for exact rational samples'
      );
    }
    const v = this.sampleExact();
    if (!v.every((x) => x.isInteger())) {
      // Possible with a non-integral centre on an integral lattice.
      throw new ValueError('sample is not integral; use sampleExact() for exact rational samples');
    }
    return v.map((x) => x.numerator);
  }

  /**
   * Alias for sample() - makes the sampler callable like a function.
   *
   * @returns A sample from the discrete Gaussian distribution
   */
  call(): bigint[] {
    return this.sample();
  }

  /**
   * Generate multiple samples.
   *
   * @param n - Number of samples to generate
   * @returns Array of samples
   */
  samples(n: number): bigint[][] {
    const result: bigint[][] = [];
    for (let i = 0; i < n; i++) {
      result.push(this.sample());
    }
    return result;
  }

  /**
   * Generate multiple exact samples.
   *
   * @param n - Number of samples to generate
   * @returns Array of exact samples
   */
  samplesExact(n: number): Rational[][] {
    const result: Rational[][] = [];
    for (let i = 0; i < n; i++) {
      result.push(this.sampleExact());
    }
    return result;
  }

  /**
   * Compute the smoothing parameter eta_epsilon(L) for this lattice.
   *
   * For practical purposes, we use the approximation:
   *   eta_epsilon(L) >= ||b_n*|| * sqrt(ln(2n(1 + 1/epsilon)) / pi)
   *
   * @see Deviation: not present in Sage; a convenience estimate.
   *
   * @param epsilon - The statistical parameter (default: 2^{-n})
   * @returns An estimate of the smoothing parameter
   */
  smoothingParameter(epsilon?: number): number {
    const n = this.rank;
    const eps = epsilon ?? 2 ** -n;
    const maxBStarNorm = Math.sqrt(Math.max(...this.gs.bStarNormsSq.map((x) => x.toNumber())));
    return maxBStarNorm * Math.sqrt(Math.log(2 * n * (1 + 1 / eps)) / Math.PI);
  }

  /**
   * Check if sigma is above the smoothing parameter.
   *
   * @param epsilon - The statistical parameter (default: 2^{-n})
   * @returns True if sigma >= eta_epsilon(L); always false for a covariance
   *   matrix, for which the scalar comparison is meaningless.
   */
  isAboveSmoothingParameter(epsilon?: number): boolean {
    if (!this.is_spherical) return false;
    return (this._sigma as number) >= this.smoothingParameter(epsilon);
  }

  /**
   * String representation.
   *
   * Reference: `discrete_gaussian_lattice.py:794-820`.
   */
  repr(): string {
    const sigma_str = this.is_spherical
      ? `σ = ${rrRepr(this._sigma as number)}`
      : `Σ =\n${matrixRepr((this._sigma as number[][]).map((r) => r.map(rrRepr)))}`;
    const cStr = `(${this._c.map((x) => x.toString()).join(', ')})`;
    const bStr = matrixRepr(this.basisExact.map((r) => r.map((x) => x.toString())));
    return `Discrete Gaussian sampler with Gaussian parameter ${sigma_str}, c=${cStr} over lattice with basis\n\n${bStr}`;
  }

  /**
   * String representation for console output.
   */
  toString(): string {
    return this.repr();
  }
}

/**
 * Factory function to create a discrete Gaussian lattice sampler.
 *
 * @param basis - The lattice basis
 * @param sigma - Gaussian parameter (scalar or covariance matrix)
 * @param c - Center (default: origin)
 * @param tau - Tail cutoff (default: 6)
 * @returns A discrete Gaussian lattice sampler
 */
export function DiscreteGaussianLattice(
  basis: number[][] | bigint[][] | Rational[][],
  sigma: number | number[][] | Rational[][],
  c?: number[] | bigint[] | Rational[],
  tau: IntegerLike = 6n
): DiscreteGaussianDistributionLatticeSampler {
  return new DiscreteGaussianDistributionLatticeSampler(basis, { sigma, c, tau });
}

/**
 * Sample a short vector from a lattice using discrete Gaussian sampling.
 *
 * This is a convenience function that creates a sampler and returns one sample.
 *
 * @param basis - The lattice basis
 * @param sigma - Standard deviation
 * @returns A short lattice vector
 */
export function sampleShortVector(
  basis: number[][] | bigint[][] | Rational[][],
  sigma: number
): bigint[] {
  const D = new DiscreteGaussianDistributionLatticeSampler(basis, { sigma });
  return D.sample();
}

/**
 * Preimage sampling: given a target t, sample a lattice vector close to t.
 *
 * This samples from D_{L, sigma, c} where c is chosen such that the sample
 * is statistically close to the nearest lattice vector to t.
 *
 * @param basis - The lattice basis
 * @param sigma - Standard deviation
 * @param target - The target vector
 * @returns A lattice vector close to target
 */
export function samplePreimage(
  basis: number[][] | bigint[][] | Rational[][],
  sigma: number,
  target: number[] | bigint[] | Rational[]
): bigint[] {
  const D = new DiscreteGaussianDistributionLatticeSampler(basis, { sigma, c: target });
  return D.sample();
}

/**
 * Discrete Gaussian sampler for polynomials (coefficient-wise).
 *
 * Samples a polynomial where each coefficient is drawn from a 1D discrete Gaussian.
 *
 * @see Deviation: SageMath's `DiscreteGaussianDistributionPolynomialSampler`
 * lives in `sage.crypto.lwe` (ported at `crypto/lwe.ts`), takes
 * `(P, n, sigma)` and returns an element of the polynomial ring `P`. This
 * class is an extra convenience wrapper living in the lattice module with
 * signature `(n, options)` that returns a coefficient array. It has no
 * SageMath counterpart; prefer `crypto/lwe.ts`'s class for parity.
 */
export class DiscreteGaussianDistributionPolynomialSampler {
  /**
   * The integer sampler used for each coefficient.
   */
  public readonly integerSampler: DiscreteGaussianDistributionIntegerSampler;

  /**
   * Number of coefficients.
   */
  public readonly n: number;

  /**
   * Construct a polynomial sampler.
   *
   * @param n - Number of coefficients
   * @param options - Options for the underlying integer sampler
   */
  constructor(n: number, options: DiscreteGaussianOptions) {
    if (n <= 0) {
      throw new ValueError(`n must be > 0, got ${n}`);
    }
    this.n = n;
    this.integerSampler = new DiscreteGaussianDistributionIntegerSampler(options);
  }

  /**
   * Sample a polynomial (as an array of coefficients).
   *
   * @returns Array of coefficients [a_0, a_1, ..., a_{n-1}]
   */
  sample(): bigint[] {
    return this.integerSampler.samples(this.n);
  }

  /**
   * Alias for sample().
   */
  call(): bigint[] {
    return this.sample();
  }

  /**
   * Generate multiple polynomial samples.
   *
   * @param count - Number of polynomials to sample
   * @returns Array of polynomial coefficient arrays
   */
  samples(count: number): bigint[][] {
    const result: bigint[][] = [];
    for (let i = 0; i < count; i++) {
      result.push(this.sample());
    }
    return result;
  }

  /**
   * Get the sigma parameter.
   */
  get sigma(): number {
    return this.integerSampler.sigma;
  }

  /**
   * Get the c parameter.
   */
  get c(): number {
    return this.integerSampler.c;
  }

  /**
   * String representation.
   */
  repr(): string {
    return `DiscreteGaussianDistributionPolynomialSampler(n=${this.n}, sigma=${this.sigma}, c=${this.c})`;
  }

  toString(): string {
    return this.repr();
  }
}

/**
 * Factory function to create a discrete Gaussian polynomial sampler.
 *
 * @param n - Number of coefficients
 * @param sigma - Standard deviation
 * @param c - Center (default: 0)
 * @param tau - Tail cutoff (default: 6)
 * @returns A discrete Gaussian polynomial sampler
 */
export function DiscreteGaussianPolynomial(
  n: number,
  sigma: number,
  c: IntegerLike = 0n,
  tau: IntegerLike = 6n
): DiscreteGaussianDistributionPolynomialSampler {
  return new DiscreteGaussianDistributionPolynomialSampler(n, { sigma, c, tau });
}
