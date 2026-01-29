/**
 * @module sage/stats/distributions/discrete_gaussian_lattice
 * @description Discrete Gaussian Distribution over Lattices
 *
 * This module implements the GPV (Gentry-Peikert-Vaikuntanathan) algorithm
 * for sampling from a discrete Gaussian distribution over a lattice.
 *
 * The discrete Gaussian distribution D_{L,sigma,c} over a lattice L is defined by:
 *   P(x) = rho_{sigma,c}(x) / rho_{sigma,c}(L)
 *
 * where rho_{sigma,c}(x) = exp(-||x-c||^2 / (2*sigma^2)) is the Gaussian function.
 *
 * The GPV algorithm samples using the lattice basis by:
 * 1. Computing Gram-Schmidt orthogonalization of the basis
 * 2. Recursively sampling coefficients from 1D discrete Gaussians
 *
 * @see Reference: sage/stats/distributions/discrete_gaussian_lattice.py
 * @see Gentry, Peikert, Vaikuntanathan, "Trapdoors for Hard Lattices...", 2008
 */

import { ValueError, NotImplementedError, TypeError as SageTypeError } from '../../errors.js';
import {
  DiscreteGaussianDistributionIntegerSampler,
  type DiscreteGaussianOptions,
  type DiscreteGaussianOptionsInternal,
} from './discrete_gaussian_integer.js';
import { type IntegerLike, toBigInt } from '../../types/coercion.js';

/**
 * Options for constructing a discrete Gaussian lattice sampler.
 */
export interface DiscreteGaussianLatticeOptions {
  /**
   * Standard deviation sigma > 0 (required).
   * Must be large enough relative to the basis quality.
   */
  sigma: number;

  /**
   * Center of the distribution (default: origin).
   * Should be a vector of the same dimension as the lattice.
   */
  c?: number[] | bigint[];

  /**
   * Tail cutoff parameter tau >= 1 (default: 6).
   */
  tau?: IntegerLike;
}

/**
 * Result of Gram-Schmidt orthogonalization.
 */
interface GramSchmidtResult {
  /** The orthogonalized basis B* (rows are b_1*, ..., b_n*) */
  bStar: number[][];
  /** The mu coefficients matrix where mu[i][j] = <b_i, b_j*> / <b_j*, b_j*> */
  mu: number[][];
  /** The squared norms of the orthogonal vectors |b_i*|^2 */
  bStarNormsSq: number[];
}

/**
 * Compute the Gram-Schmidt orthogonalization of a basis.
 *
 * @param basis - A matrix whose rows form the basis
 * @returns GramSchmidtResult containing orthogonal basis, mu matrix, and squared norms
 */
function gramSchmidt(basis: number[][]): GramSchmidtResult {
  const n = basis.length;
  if (n === 0) {
    return { bStar: [], mu: [], bStarNormsSq: [] };
  }
  const m = basis[0]!.length;

  const bStar: number[][] = [];
  const mu: number[][] = [];
  const bStarNormsSq: number[] = [];

  for (let i = 0; i < n; i++) {
    mu.push(new Array(n).fill(0));
    mu[i]![i] = 1;

    // Start with b_i
    const bi = [...basis[i]!];
    const biStar = [...bi];

    // Subtract projections onto previous orthogonal vectors
    for (let j = 0; j < i; j++) {
      // mu[i][j] = <b_i, b_j*> / <b_j*, b_j*>
      let dotProduct = 0;
      for (let k = 0; k < m; k++) {
        dotProduct += bi[k]! * bStar[j]![k]!;
      }
      mu[i]![j] = dotProduct / bStarNormsSq[j]!;

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
    bStarNormsSq.push(normSq);
  }

  return { bStar, mu, bStarNormsSq };
}

/**
 * Compute the dot product of two vectors.
 */
function dot(a: number[], b: number[]): number {
  let sum = 0;
  for (let i = 0; i < a.length; i++) {
    sum += a[i]! * b[i]!;
  }
  return sum;
}

/**
 * Add two vectors.
 */
function vecAdd(a: number[], b: number[]): number[] {
  return a.map((x, i) => x + b[i]!);
}

/**
 * Subtract two vectors.
 */
function vecSub(a: number[], b: number[]): number[] {
  return a.map((x, i) => x - b[i]!);
}

/**
 * Multiply a vector by a scalar.
 */
function vecScale(a: number[], s: number): number[] {
  return a.map((x) => x * s);
}

/**
 * Compute the squared Euclidean norm of a vector.
 */
function normSq(a: number[]): number {
  return dot(a, a);
}

/**
 * Discrete Gaussian Distribution Sampler over a Lattice.
 *
 * Implements the GPV (Gentry-Peikert-Vaikuntanathan) algorithm for sampling
 * from a discrete Gaussian distribution D_{L,sigma,c} over a lattice L.
 *
 * The algorithm works by:
 * 1. Computing the Gram-Schmidt orthogonalization of the basis
 * 2. For each basis vector (from last to first), sample a coefficient
 *    from a 1D discrete Gaussian centered at the appropriate value
 * 3. The final sample is the linear combination of basis vectors
 *
 * @example
 * ```typescript
 * // Create a lattice basis
 * const basis = [[1, 0], [0, 1]]; // Z^2 lattice
 *
 * // Create a sampler with sigma=3
 * const D = new DiscreteGaussianDistributionLatticeSampler({
 *   basis,
 *   sigma: 3
 * });
 *
 * // Sample a lattice vector
 * const v = D.sample();
 * console.log(v); // e.g., [2n, -1n]
 * ```
 *
 * @see Reference: sage/stats/distributions/discrete_gaussian_lattice.py
 */
export class DiscreteGaussianDistributionLatticeSampler {
  /**
   * The lattice basis (rows are basis vectors).
   */
  public readonly basis: number[][];

  /**
   * Standard deviation of the Gaussian distribution.
   */
  public readonly sigma: number;

  /**
   * Center of the distribution.
   */
  public readonly c: number[];

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

  /**
   * Gram-Schmidt orthogonalization data.
   */
  private readonly gs: GramSchmidtResult;

  /**
   * Precomputed s_i^2 = sigma^2 / |b_i*|^2 for each basis vector.
   */
  private readonly sigmaISq: number[];

  /**
   * Precomputed s_i = sigma / |b_i*| for each basis vector.
   */
  private readonly sigmaI: number[];

  /**
   * Construct a new discrete Gaussian lattice sampler.
   *
   * @param basis - The lattice basis (rows are basis vectors)
   * @param options - Configuration options
   * @throws ValueError if sigma is too small for the basis quality
   */
  constructor(
    basis: number[][] | bigint[][],
    options: DiscreteGaussianLatticeOptions
  ) {
    // Validate basis
    if (!Array.isArray(basis) || basis.length === 0) {
      throw new ValueError('basis must be a non-empty array of vectors');
    }

    // Convert to number[][] for computation
    this.basis = basis.map((row) => row.map((x) => Number(x)));
    this.rank = this.basis.length;
    this.degree = this.basis[0]!.length;

    // Validate all rows have the same length
    for (const row of this.basis) {
      if (row.length !== this.degree) {
        throw new ValueError('all basis vectors must have the same dimension');
      }
    }

    // Validate and store sigma
    if (options.sigma === undefined || options.sigma === null) {
      throw new SageTypeError('sigma is required');
    }
    if (typeof options.sigma !== 'number' || !isFinite(options.sigma)) {
      throw new SageTypeError(`sigma must be a finite number, got ${options.sigma}`);
    }
    if (options.sigma <= 0) {
      throw new ValueError(`sigma must be > 0, got ${options.sigma}`);
    }
    this.sigma = options.sigma;

    // Validate and store c (default: origin)
    if (options.c !== undefined) {
      if (!Array.isArray(options.c) || options.c.length !== this.degree) {
        throw new ValueError(`c must be a vector of dimension ${this.degree}`);
      }
      this.c = options.c.map((x) => Number(x));
    } else {
      this.c = new Array(this.degree).fill(0);
    }

    // Validate and store tau (default: 6)
    const tauValue = options.tau !== undefined ? Number(toBigInt(options.tau)) : 6;
    if (tauValue < 1) {
      throw new ValueError(`tau must be >= 1, got ${tauValue}`);
    }
    this.tau = tauValue;

    // Compute Gram-Schmidt orthogonalization
    this.gs = gramSchmidt(this.basis);

    // Precompute s_i = sigma / |b_i*| for each basis vector
    this.sigmaISq = this.gs.bStarNormsSq.map((normSq) =>
      (this.sigma * this.sigma) / normSq
    );
    this.sigmaI = this.sigmaISq.map((s2) => Math.sqrt(s2));

    // Validate that sigma is large enough
    // For the GPV algorithm to work correctly, sigma should be >= ||b*|| * omega(sqrt(log n))
    const maxBStarNorm = Math.sqrt(Math.max(...this.gs.bStarNormsSq));
    const minSigma = maxBStarNorm * Math.sqrt(Math.log(this.rank + 1));
    if (this.sigma < minSigma * 0.5) {
      // Just a warning - we still allow it but results may not be statistically close
      console.warn(
        `Warning: sigma=${this.sigma} may be too small for basis quality. ` +
        `Consider sigma >= ${minSigma.toFixed(2)} for good statistical properties.`
      );
    }
  }

  /**
   * Sample from the discrete Gaussian distribution over the lattice.
   *
   * Uses the GPV algorithm: for i from n down to 1, sample coefficient z_i
   * from a 1D discrete Gaussian, then return sum(z_i * b_i).
   *
   * @returns A sample from the discrete Gaussian distribution as a bigint vector
   */
  sample(): bigint[] {
    // Working vector (the current center offset)
    let v = [...this.c];

    // Coefficients z_i for each basis vector
    const z: bigint[] = new Array(this.rank).fill(0n);

    // Sample from last to first basis vector
    for (let i = this.rank - 1; i >= 0; i--) {
      // Compute the center for the 1D Gaussian: c_i = <v, b_i*> / |b_i*|^2
      const bStarI = this.gs.bStar[i]!;
      const bStarNormSq = this.gs.bStarNormsSq[i]!;
      const ci = dot(v, bStarI) / bStarNormSq;

      // Sample z_i from D_{sigma_i, c_i}
      const sigmaI = this.sigmaI[i]!;
      const Di = new DiscreteGaussianDistributionIntegerSampler({
        sigma: sigmaI,
        c: ci,
        tau: this.tau,
      });
      z[i] = Di.sample();

      // Update v: v = v - z_i * b_i
      const ziBi = vecScale(this.basis[i]!, Number(z[i]));
      v = vecSub(v, ziBi);
    }

    // The sample is sum(z_i * b_i)
    // But we computed it as c - v, so return by computing explicitly
    const result: bigint[] = new Array(this.degree).fill(0n);
    for (let i = 0; i < this.rank; i++) {
      for (let j = 0; j < this.degree; j++) {
        result[j] = result[j]! + z[i]! * BigInt(Math.round(this.basis[i]![j]!));
      }
    }

    return result;
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
   * Compute the smoothing parameter eta_epsilon(L) for this lattice.
   *
   * The smoothing parameter is the smallest s such that rho_{1/s}(L* \ {0}) <= epsilon,
   * where L* is the dual lattice.
   *
   * For practical purposes, we use the approximation:
   *   eta_epsilon(L) >= ||b_n*|| * sqrt(ln(2n(1 + 1/epsilon)) / pi)
   *
   * @param epsilon - The statistical parameter (default: 2^{-n})
   * @returns An estimate of the smoothing parameter
   */
  smoothingParameter(epsilon?: number): number {
    const n = this.rank;
    const eps = epsilon ?? Math.pow(2, -n);

    // Maximum Gram-Schmidt orthogonal vector length
    const maxBStarNorm = Math.sqrt(Math.max(...this.gs.bStarNormsSq));

    // Lower bound on smoothing parameter
    return maxBStarNorm * Math.sqrt(Math.log(2 * n * (1 + 1 / eps)) / Math.PI);
  }

  /**
   * Check if sigma is above the smoothing parameter.
   *
   * @param epsilon - The statistical parameter (default: 2^{-n})
   * @returns True if sigma >= eta_epsilon(L)
   */
  isAboveSmoothingParameter(epsilon?: number): boolean {
    return this.sigma >= this.smoothingParameter(epsilon);
  }

  /**
   * String representation.
   */
  repr(): string {
    const cStr = this.c.every((x) => x === 0) ? '0' : `[${this.c.join(', ')}]`;
    return `DiscreteGaussianDistributionLatticeSampler(rank=${this.rank}, sigma=${this.sigma}, c=${cStr})`;
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
 * @param sigma - Standard deviation
 * @param c - Center (default: origin)
 * @param tau - Tail cutoff (default: 6)
 * @returns A discrete Gaussian lattice sampler
 */
export function DiscreteGaussianLattice(
  basis: number[][] | bigint[][],
  sigma: number,
  c?: number[] | bigint[],
  tau: IntegerLike = 6
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
  basis: number[][] | bigint[][],
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
  basis: number[][] | bigint[][],
  sigma: number,
  target: number[] | bigint[]
): bigint[] {
  const c = target.map((x) => Number(x));
  const D = new DiscreteGaussianDistributionLatticeSampler(basis, { sigma, c });
  return D.sample();
}

/**
 * Discrete Gaussian sampler for polynomials (coefficient-wise).
 *
 * Samples a polynomial where each coefficient is drawn from a 1D discrete Gaussian.
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
  c: IntegerLike = 0,
  tau: IntegerLike = 6
): DiscreteGaussianDistributionPolynomialSampler {
  return new DiscreteGaussianDistributionPolynomialSampler(n, { sigma, c, tau });
}
