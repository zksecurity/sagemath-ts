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

import { TypeError as SageTypeError, ValueError } from '../../errors.js';
import { Rational } from '../../rings/rational.js';
import { type IntegerLike, toBigInt, toSafeNumber } from '../../types/coercion.js';
import {
  DiscreteGaussianDistributionIntegerSampler,
  type DiscreteGaussianOptions,
  type DiscreteGaussianOptionsInternal,
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
  c?: number[] | bigint[] | Rational[];

  /**
   * Tail cutoff parameter tau >= 1 (default: 6).
   */
  tau?: IntegerLike;
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
 * (`discrete_gaussian_lattice.py:585`), which returns rationals — not the
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
   * The lattice basis (rows are basis vectors), as floating point.
   *
   * This is a lossy view kept for convenience; all sampling arithmetic uses
   * {@link basisExact}.
   */
  public readonly basis: number[][];

  /**
   * The lattice basis (rows are basis vectors), exactly.
   */
  public readonly basisExact: Rational[][];

  /**
   * Whether every basis entry is an integer, i.e. the lattice sits inside Z^d.
   *
   * When false, {@link sample} cannot return integer coordinates and
   * {@link sampleExact} must be used instead.
   */
  public readonly isIntegral: boolean;

  /**
   * Standard deviation of the Gaussian distribution.
   */
  public readonly sigma: number;

  /**
   * Center of the distribution (floating-point view).
   */
  public readonly c: number[];

  /**
   * Center of the distribution, exactly.
   */
  public readonly cExact: Rational[];

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
   *
   * Sage computes this in `RealField(precision)`
   * (`discrete_gaussian_lattice.py:866`), so a double is faithful here.
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
    this.basis = this.basisExact.map((row) => row.map((x) => x.toNumber()));
    this.isIntegral = this.basisExact.every((row) => row.every((x) => x.isInteger()));
    this.rank = this.basisExact.length;
    this.degree = this.basisExact[0]!.length;

    // Validate all rows have the same length
    for (const row of this.basisExact) {
      if (row.length !== this.degree) {
        throw new ValueError('all basis vectors must have the same dimension');
      }
    }

    // Validate and store sigma
    if (options.sigma === undefined || options.sigma === null) {
      throw new SageTypeError('sigma is required');
    }
    if (typeof options.sigma !== 'number' || !Number.isFinite(options.sigma)) {
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
      this.cExact = (options.c as (number | bigint | Rational)[]).map((x) => toRationalEntry(x));
    } else {
      this.cExact = new Array(this.degree).fill(Rational.zero());
    }
    this.c = this.cExact.map((x) => x.toNumber());

    // Validate and store tau (default: 6)
    const tauValue = options.tau !== undefined ? toSafeNumber(toBigInt(options.tau)) : 6;
    if (tauValue < 1) {
      throw new ValueError(`tau must be >= 1, got ${tauValue}`);
    }
    this.tau = tauValue;

    // Compute the exact Gram-Schmidt orthogonalization
    this.gs = gramSchmidt(this.basisExact);

    // Precompute s_i = sigma / |b_i*| for each basis vector
    this.sigmaISq = this.gs.bStarNormsSq.map(
      (normSq) => (this.sigma * this.sigma) / normSq.toNumber()
    );
    this.sigmaI = this.sigmaISq.map((s2) => Math.sqrt(s2));

    // Validate that sigma is large enough
    // For the GPV algorithm to work correctly, sigma should be >= ||b*|| * omega(sqrt(log n))
    const maxBStarNorm = Math.sqrt(Math.max(...this.gs.bStarNormsSq.map((x) => x.toNumber())));
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
   * Sample from the discrete Gaussian distribution over the lattice, exactly.
   *
   * This is Sage's `_call` (`discrete_gaussian_lattice.py:842-872`), which
   * runs entirely over the exact base ring of the basis:
   *
   * ```
   * v = 0
   * for i in range(m - 1, -1, -1):
   *     b_ = self._G[i]
   *     c_ = c.dot_product(b_) / b_.dot_product(b_)
   *     sigma_ = sigma / b_.norm()
   *     z = DiscreteGaussianDistributionIntegerSampler(sigma=sigma_, c=c_,
   *                                                    algorithm='uniform+online')()
   *     c = c - z * B[i]
   *     v = v + z * B[i]
   * ```
   *
   * @returns A sample from the discrete Gaussian distribution as an exact
   *   rational vector
   */
  sampleExact(): Rational[] {
    // c is updated in place as in Sage's `_call`.
    let c = [...this.cExact];
    let v: Rational[] = new Array(this.degree).fill(Rational.zero());

    // Sample from last to first basis vector
    for (let i = this.rank - 1; i >= 0; i--) {
      // c_ = <c, b_i*> / <b_i*, b_i*>  (exact)
      const bStarI = this.gs.bStar[i]!;
      const bStarNormSq = this.gs.bStarNormsSq[i]!;
      const ci = ratDot(c, bStarI).div(bStarNormSq);

      // sigma_ = sigma / |b_i*|  (RealField in Sage)
      const sigmaI = this.sigmaI[i]!;
      const Di = new DiscreteGaussianDistributionIntegerSampler({
        sigma: sigmaI,
        c: ci.toNumber(),
        tau: this.tau,
        // Sage passes this explicitly, avoiding a sigma*tau-sized rho table
        // for every coordinate of every sample.
        algorithm: 'uniform+online',
      });
      const z = Di.sample();

      // c = c - z * B[i];  v = v + z * B[i]   (exact)
      const zRat = new Rational(z, 1n);
      const bi = this.basisExact[i]!;
      c = c.map((x, j) => x.sub(zRat.mul(bi[j]!)));
      v = v.map((x, j) => x.add(zRat.mul(bi[j]!)));
    }

    return v;
  }

  /**
   * Sample from the discrete Gaussian distribution over the lattice.
   *
   * Uses the GPV algorithm: for i from n down to 1, sample coefficient z_i
   * from a 1D discrete Gaussian, then return sum(z_i * b_i).
   *
   * @returns A sample from the discrete Gaussian distribution as a bigint vector
   * @throws ValueError if the lattice basis is not integral; use
   *   {@link sampleExact} in that case
   */
  sample(): bigint[] {
    if (!this.isIntegral) {
      throw new ValueError(
        'lattice basis is not integral; use sampleExact() for exact rational samples'
      );
    }
    return this.sampleExact().map((x) => x.numerator);
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
    const eps = epsilon ?? 2 ** -n;

    // Maximum Gram-Schmidt orthogonal vector length
    const maxBStarNorm = Math.sqrt(Math.max(...this.gs.bStarNormsSq.map((x) => x.toNumber())));

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
  basis: number[][] | bigint[][] | Rational[][],
  sigma: number,
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
