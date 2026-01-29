/**
 * @module sage/stats/distributions/discrete_gaussian_integer
 * @description Discrete Gaussian Distribution over the Integers
 *
 * This module implements a sampler for the discrete Gaussian distribution
 * over the integers. The discrete Gaussian distribution D_{sigma,c} is defined
 * by the probability density function:
 *
 *   rho_{sigma,c}(x) = exp(-(x-c)^2 / (2*sigma^2))
 *   P(x) = rho_{sigma,c}(x) / sum_{y in Z} rho_{sigma,c}(y)
 *
 * This is commonly used in lattice-based cryptography for error sampling.
 *
 * @see Reference: sage/stats/distributions/discrete_gaussian_integer.pyx
 */

import { ValueError, TypeError as SageTypeError } from '../../errors.js';
import { current_randstate } from '../../misc/randstate.js';
import { type IntegerLike, toBigInt } from '../../types/coercion.js';

/**
 * Algorithm choices for discrete Gaussian sampling.
 */
export type DiscreteGaussianAlgorithm = 'uniform+table' | 'uniform+online';

/**
 * Options for constructing a discrete Gaussian sampler.
 * Accepts IntegerLike for c and tau for user-facing API.
 */
export interface DiscreteGaussianOptions {
  /**
   * Standard deviation sigma > 0 (required).
   */
  sigma: number;

  /**
   * Center of the distribution (default: 0).
   */
  c?: IntegerLike;

  /**
   * Tail cutoff parameter tau >= 1 (default: 6).
   * Samples are drawn from [floor(c) - ceil(sigma*tau), floor(c) + ceil(sigma*tau)].
   */
  tau?: IntegerLike;

  /**
   * Algorithm to use:
   * - 'uniform+table': Precompute probability table (faster for repeated sampling)
   * - 'uniform+online': Compute probabilities on-the-fly (less memory)
   *
   * Default: 'uniform+table' when sigma*tau <= 10^6, otherwise 'uniform+online'
   */
  algorithm?: DiscreteGaussianAlgorithm;
}

/**
 * Internal options interface that allows number for c (needed by GPV algorithm).
 * @internal
 */
export interface DiscreteGaussianOptionsInternal {
  sigma: number;
  c?: number;
  tau?: number;
  algorithm?: DiscreteGaussianAlgorithm;
}

/**
 * Discrete Gaussian Distribution Sampler over the Integers.
 *
 * A sampler for the discrete Gaussian distribution D_{sigma,c} over the integers.
 * The distribution is characterized by:
 * - sigma: standard deviation (determines the spread)
 * - c: center of the distribution
 * - tau: tail cutoff (samples lie within tau standard deviations of center)
 *
 * @example
 * ```typescript
 * // Create a sampler with sigma=3, centered at 0
 * const D = new DiscreteGaussianDistributionIntegerSampler({ sigma: 3 });
 *
 * // Draw samples
 * const sample = D.sample();
 * console.log(sample); // e.g., 2n, -1n, 0n, etc.
 *
 * // Create with custom center and tail cutoff
 * const D2 = new DiscreteGaussianDistributionIntegerSampler({
 *   sigma: 5,
 *   c: 10,
 *   tau: 4
 * });
 * ```
 *
 * @see Reference: sage/stats/distributions/discrete_gaussian_integer.pyx:DiscreteGaussianDistributionIntegerSampler
 */
export class DiscreteGaussianDistributionIntegerSampler {
  /**
   * Standard deviation of the Gaussian distribution.
   */
  public readonly sigma: number;

  /**
   * Center of the distribution.
   */
  public readonly c: number;

  /**
   * Tail cutoff parameter.
   */
  public readonly tau: number;

  /**
   * The algorithm used for sampling.
   */
  public readonly algorithm: DiscreteGaussianAlgorithm;

  /**
   * Lower bound of the sampling range (inclusive).
   */
  public readonly lowerBound: bigint;

  /**
   * Upper bound of the sampling range (inclusive).
   */
  public readonly upperBound: bigint;

  /**
   * Precomputed probability table for 'uniform+table' algorithm.
   * Maps x -> rho(x) where rho(x) = exp(-(x-c)^2 / (2*sigma^2))
   */
  private readonly rhoTable: Map<bigint, number> | null;

  /**
   * Maximum value in the rho table (for rejection sampling).
   */
  private readonly rhoMax: number;

  /**
   * Precomputed value: -1 / (2 * sigma^2)
   */
  private readonly negHalfInvSigmaSq: number;

  /**
   * Construct a new discrete Gaussian sampler.
   *
   * @param options - Configuration options
   * @throws ValueError if sigma <= 0 or tau < 1
   *
   * @example
   * ```typescript
   * const sampler = new DiscreteGaussianDistributionIntegerSampler({
   *   sigma: 3.0,
   *   c: 0,
   *   tau: 6
   * });
   * ```
   */
  constructor(options: DiscreteGaussianOptions | DiscreteGaussianOptionsInternal) {
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

    // Validate and store c (default: 0)
    // Support both IntegerLike (public API) and number (internal/GPV algorithm)
    let cValue: number;
    if (options.c === undefined) {
      cValue = 0;
    } else if (typeof options.c === 'number') {
      if (!isFinite(options.c)) {
        throw new SageTypeError(`c must be a finite number, got ${options.c}`);
      }
      cValue = options.c;
    } else {
      // IntegerLike (bigint or Integer)
      cValue = Number(toBigInt(options.c));
    }
    this.c = cValue;

    // Validate and store tau (default: 6)
    // Support both IntegerLike (public API) and number (internal)
    let tauValue: number;
    if (options.tau === undefined) {
      tauValue = 6;
    } else if (typeof options.tau === 'number') {
      if (!isFinite(options.tau)) {
        throw new SageTypeError(`tau must be a finite number, got ${options.tau}`);
      }
      tauValue = options.tau;
    } else {
      // IntegerLike (bigint or Integer)
      tauValue = Number(toBigInt(options.tau));
    }
    if (tauValue < 1) {
      throw new ValueError(`tau must be >= 1, got ${tauValue}`);
    }
    this.tau = tauValue;

    // Precompute negHalfInvSigmaSq = -1 / (2 * sigma^2)
    this.negHalfInvSigmaSq = -1.0 / (2.0 * this.sigma * this.sigma);

    // Compute sampling range: [floor(c) - ceil(sigma*tau), floor(c) + ceil(sigma*tau)]
    const floorC = Math.floor(this.c);
    const halfWidth = Math.ceil(this.sigma * this.tau);
    this.lowerBound = BigInt(floorC - halfWidth);
    this.upperBound = BigInt(floorC + halfWidth);

    // Choose algorithm based on range size
    const rangeSize = Number(this.upperBound - this.lowerBound + 1n);
    if (options.algorithm) {
      this.algorithm = options.algorithm;
    } else {
      // Default: use table when range is reasonable (sigma*tau <= 10^6)
      this.algorithm = rangeSize <= 1_000_000 ? 'uniform+table' : 'uniform+online';
    }

    // The maximum value of rho(x) occurs at x closest to c
    // For computational stability, we use rho(x) without normalizing
    this.rhoMax = this._rho(this.c);

    // Precompute table if using 'uniform+table'
    if (this.algorithm === 'uniform+table') {
      this.rhoTable = new Map();
      for (let x = this.lowerBound; x <= this.upperBound; x++) {
        const rhoX = this._rho(Number(x));
        this.rhoTable.set(x, rhoX);
      }
    } else {
      this.rhoTable = null;
    }
  }

  /**
   * Compute the Gaussian weight rho_{sigma,c}(x) = exp(-(x-c)^2 / (2*sigma^2)).
   *
   * @param x - The point to evaluate
   * @returns The Gaussian weight at x
   */
  private _rho(x: number): number {
    const diff = x - this.c;
    return Math.exp(diff * diff * this.negHalfInvSigmaSq);
  }

  /**
   * Sample from the discrete Gaussian distribution.
   *
   * Uses rejection sampling: draw x uniformly from the range, accept with
   * probability rho(x)/rho_max.
   *
   * @returns A sample from the discrete Gaussian distribution
   *
   * @example
   * ```typescript
   * const D = new DiscreteGaussianDistributionIntegerSampler({ sigma: 3 });
   * const x = D.sample();
   * console.log(x); // e.g., -2n
   * ```
   */
  sample(): bigint {
    const range = this.upperBound - this.lowerBound + 1n;
    const rstate = current_randstate();

    // Rejection sampling loop
    while (true) {
      // Draw x uniformly from [lowerBound, upperBound]
      const x = this._uniformSample(range);

      // Get rho(x)
      let rhoX: number;
      if (this.algorithm === 'uniform+table' && this.rhoTable !== null) {
        rhoX = this.rhoTable.get(x)!;
      } else {
        rhoX = this._rho(Number(x));
      }

      // Accept with probability rho(x) / rhoMax
      const acceptProb = rhoX / this.rhoMax;
      if (rstate.random() < acceptProb) {
        return x;
      }
    }
  }

  /**
   * Sample a uniform random integer in [lowerBound, lowerBound + range - 1].
   */
  private _uniformSample(range: bigint): bigint {
    const rstate = current_randstate();
    const randomOffset = rstate.random_below(range);
    return this.lowerBound + randomOffset;
  }

  /**
   * Alias for sample() - makes the sampler callable like a function.
   * Matches SageMath's __call__ interface.
   *
   * @returns A sample from the discrete Gaussian distribution
   */
  call(): bigint {
    return this.sample();
  }

  /**
   * Generate multiple samples.
   *
   * @param n - Number of samples to generate
   * @returns Array of samples
   *
   * @example
   * ```typescript
   * const D = new DiscreteGaussianDistributionIntegerSampler({ sigma: 3 });
   * const samples = D.samples(100);
   * ```
   */
  samples(n: number): bigint[] {
    const result: bigint[] = [];
    for (let i = 0; i < n; i++) {
      result.push(this.sample());
    }
    return result;
  }

  /**
   * Compute the (unnormalized) probability mass at a given point.
   *
   * Returns rho(x) = exp(-(x-c)^2 / (2*sigma^2)).
   *
   * @param x - The point to evaluate
   * @returns The unnormalized probability mass at x
   *
   * @example
   * ```typescript
   * const D = new DiscreteGaussianDistributionIntegerSampler({ sigma: 3 });
   * const p0 = D.rho(0);  // Highest at center
   * const p5 = D.rho(5);  // Lower away from center
   * ```
   */
  rho(x: number | bigint): number {
    const xNum = typeof x === 'bigint' ? Number(x) : x;
    return this._rho(xNum);
  }

  /**
   * Compute the support (the set of integers with non-negligible probability).
   *
   * @returns Array of integers in the support
   */
  support(): bigint[] {
    const result: bigint[] = [];
    for (let x = this.lowerBound; x <= this.upperBound; x++) {
      result.push(x);
    }
    return result;
  }

  /**
   * Compute the normalization constant (sum of rho(x) over support).
   *
   * @returns The normalization constant
   */
  normalizationConstant(): number {
    let sum = 0;
    for (let x = this.lowerBound; x <= this.upperBound; x++) {
      sum += this._rho(Number(x));
    }
    return sum;
  }

  /**
   * Compute the normalized probability at a given point.
   *
   * @param x - The point to evaluate
   * @returns P(X = x) for the discrete Gaussian distribution
   */
  probability(x: number | bigint): number {
    const xBig = typeof x === 'bigint' ? x : BigInt(Math.round(x));
    if (xBig < this.lowerBound || xBig > this.upperBound) {
      return 0;
    }
    return this.rho(x) / this.normalizationConstant();
  }

  /**
   * Compute the expected value (mean) of the distribution.
   *
   * @returns The expected value (approximately c)
   */
  mean(): number {
    const Z = this.normalizationConstant();
    let sum = 0;
    for (let x = this.lowerBound; x <= this.upperBound; x++) {
      sum += Number(x) * this._rho(Number(x));
    }
    return sum / Z;
  }

  /**
   * Compute the variance of the distribution.
   *
   * @returns The variance (approximately sigma^2)
   */
  variance(): number {
    const mu = this.mean();
    const Z = this.normalizationConstant();
    let sum = 0;
    for (let x = this.lowerBound; x <= this.upperBound; x++) {
      const diff = Number(x) - mu;
      sum += diff * diff * this._rho(Number(x));
    }
    return sum / Z;
  }

  /**
   * Compute the standard deviation of the distribution.
   *
   * @returns The standard deviation (approximately sigma)
   */
  stddev(): number {
    return Math.sqrt(this.variance());
  }

  /**
   * String representation.
   */
  repr(): string {
    return `DiscreteGaussianDistributionIntegerSampler(sigma=${this.sigma}, c=${this.c}, tau=${this.tau})`;
  }

  /**
   * String representation for console output.
   */
  toString(): string {
    return this.repr();
  }

  /**
   * Create a copy of this sampler with modified parameters.
   *
   * @param options - Parameters to override
   * @returns A new sampler with the specified parameters
   */
  withOptions(options: Partial<DiscreteGaussianOptions>): DiscreteGaussianDistributionIntegerSampler {
    return new DiscreteGaussianDistributionIntegerSampler({
      sigma: options.sigma ?? this.sigma,
      c: options.c ?? this.c,
      tau: options.tau ?? this.tau,
      algorithm: options.algorithm ?? this.algorithm,
    });
  }
}

/**
 * Factory function to create a discrete Gaussian sampler.
 *
 * @param sigma - Standard deviation
 * @param c - Center (default: 0)
 * @param tau - Tail cutoff (default: 6)
 * @param algorithm - Sampling algorithm (default: auto-select)
 * @returns A discrete Gaussian sampler
 *
 * @example
 * ```typescript
 * const D = DiscreteGaussianInteger(3);
 * const sample = D.sample();
 * ```
 */
export function DiscreteGaussianInteger(
  sigma: number,
  c: IntegerLike = 0,
  tau: IntegerLike = 6,
  algorithm?: DiscreteGaussianAlgorithm
): DiscreteGaussianDistributionIntegerSampler {
  return new DiscreteGaussianDistributionIntegerSampler({ sigma, c, tau, algorithm });
}

/**
 * Compute the statistical distance between two discrete Gaussian distributions.
 *
 * The statistical distance (total variation distance) between distributions P and Q is:
 *   delta(P, Q) = (1/2) * sum_x |P(x) - Q(x)|
 *
 * @param D1 - First distribution
 * @param D2 - Second distribution
 * @returns The statistical distance in [0, 1]
 */
export function statisticalDistance(
  D1: DiscreteGaussianDistributionIntegerSampler,
  D2: DiscreteGaussianDistributionIntegerSampler
): number {
  // Compute the union of supports
  const minX = D1.lowerBound < D2.lowerBound ? D1.lowerBound : D2.lowerBound;
  const maxX = D1.upperBound > D2.upperBound ? D1.upperBound : D2.upperBound;

  let sum = 0;
  for (let x = minX; x <= maxX; x++) {
    const p1 = D1.probability(x);
    const p2 = D2.probability(x);
    sum += Math.abs(p1 - p2);
  }

  return sum / 2;
}

/**
 * Compute the Kullback-Leibler divergence from D1 to D2.
 *
 * KL(D1 || D2) = sum_x P1(x) * log(P1(x) / P2(x))
 *
 * @param D1 - First distribution
 * @param D2 - Second distribution
 * @returns The KL divergence (may be Infinity if supports don't overlap)
 */
export function klDivergence(
  D1: DiscreteGaussianDistributionIntegerSampler,
  D2: DiscreteGaussianDistributionIntegerSampler
): number {
  let sum = 0;

  for (let x = D1.lowerBound; x <= D1.upperBound; x++) {
    const p1 = D1.probability(x);
    if (p1 === 0) continue;

    const p2 = D2.probability(x);
    if (p2 === 0) return Infinity;

    sum += p1 * Math.log(p1 / p2);
  }

  return sum;
}

/**
 * Sample from the convolution of two discrete Gaussian distributions.
 *
 * If X ~ D_{sigma1, c1} and Y ~ D_{sigma2, c2} independently, then
 * X + Y ~ D_{sqrt(sigma1^2 + sigma2^2), c1 + c2} (approximately).
 *
 * @param D1 - First distribution
 * @param D2 - Second distribution
 * @returns A sample from the convolution
 */
export function sampleConvolution(
  D1: DiscreteGaussianDistributionIntegerSampler,
  D2: DiscreteGaussianDistributionIntegerSampler
): bigint {
  return D1.sample() + D2.sample();
}
