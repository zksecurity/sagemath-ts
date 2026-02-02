/**
 * @module sage/crypto/lwe
 * @description (Ring-)LWE oracle generators
 *
 * The Learning with Errors problem (LWE) is solving linear systems of equations
 * where the right hand side has been disturbed 'slightly' where 'slightly' is made
 * precise by a noise distribution - typically a discrete Gaussian distribution.
 *
 * The Ring Learning with Errors problem (Ring-LWE) is solving a set of univariate
 * polynomial equations - typically in a cyclotomic field - where the right hand
 * side was disturbed 'slightly'.
 *
 * This module implements generators of LWE samples where parameters are chosen
 * following proposals in the cryptographic literature.
 *
 * @see Reference: sage/crypto/lwe.py
 */

import { euler_phi, isqrt, next_prime } from '../arith/misc.js';
import { NotImplementedError, TypeError as SageTypeError, ValueError } from '../errors.js';
import { current_randstate } from '../misc/randstate.js';
import type { IntegerMod } from '../rings/finite_rings/integer_mod.js';
import { type IntegerModRing, Zmod } from '../rings/finite_rings/integer_mod_ring.js';
import type { Polynomial, RingElement } from '../rings/polynomial/polynomial_element.js';
import { PolynomialRing, PolynomialRingConstructor } from '../rings/polynomial/polynomial_ring.js';
import { QuotientRing, type QuotientRingElement } from '../rings/polynomial/quotient_ring.js';
import { DiscreteGaussianDistributionIntegerSampler } from '../stats/distributions/discrete_gaussian_integer.js';
import { type IntegerLike, toBigInt, toSafeNumber } from '../types/coercion.js';

/**
 * Distribution type for secret generation.
 */
export type SecretDistribution = 'uniform' | 'noise' | [number, number];

/**
 * Instance type for UniformNoiseLWE.
 */
export type UniformNoiseLWEInstance = 'key' | 'encrypt';

/**
 * Interface for any distribution sampler (can be used for errors or secrets).
 */
export interface DistributionSampler {
  /**
   * Sample from the distribution.
   */
  call(): bigint;
}

/**
 * Uniform sampling in a range of integers.
 *
 * @example
 * const sampler = new UniformSampler(-2, 2);
 * const sample = sampler.call(); // returns integer in range [-2, 2]
 *
 * @see Reference: sage/crypto/lwe.py:UniformSampler
 */
export class UniformSampler implements DistributionSampler {
  public readonly lower_bound: bigint;
  public readonly upper_bound: bigint;

  /**
   * Construct a uniform sampler with bounds `lower_bound` and
   * `upper_bound` (both endpoints inclusive).
   *
   * @param lower_bound - Lower bound (inclusive)
   * @param upper_bound - Upper bound (inclusive)
   */
  constructor(lower_bound: IntegerLike, upper_bound: IntegerLike) {
    const lb = toBigInt(lower_bound);
    const ub = toBigInt(upper_bound);

    if (lb > ub) {
      throw new SageTypeError('lower bound must be <= upper bound.');
    }

    this.lower_bound = lb;
    this.upper_bound = ub;
  }

  /**
   * Return a new sample.
   *
   * @returns A random integer in [lower_bound, upper_bound]
   */
  call(): bigint {
    const range = this.upper_bound - this.lower_bound + 1n;
    const rstate = current_randstate();

    const randomOffset = rstate.random_below(range);
    return this.lower_bound + randomOffset;
  }

  /**
   * String representation.
   */
  repr(): string {
    return `UniformSampler(${this.lower_bound}, ${this.upper_bound})`;
  }

  toString(): string {
    return this.repr();
  }
}

/**
 * Uniform sampler for polynomials.
 *
 * Samples univariate polynomials of degree n-1 where coefficients are
 * drawn uniformly at random between lower_bound and upper_bound.
 *
 * @example
 * const [R, x] = PolynomialRingConstructor(Zmod(101), 'x');
 * const sampler = new UniformPolynomialSampler(R, 8, -2, 2);
 * const poly = sampler.call();
 *
 * @see Reference: sage/crypto/lwe.py:UniformPolynomialSampler
 */
export class UniformPolynomialSampler<C extends RingElement = IntegerMod> {
  public readonly n: number;
  public readonly P: PolynomialRing<C>;
  public readonly lower_bound: bigint;
  public readonly upper_bound: bigint;
  public readonly D: UniformSampler;

  /**
   * Construct a sampler for univariate polynomials of degree n-1 where
   * coefficients are drawn uniformly at random between lower_bound and
   * upper_bound (both endpoints inclusive).
   *
   * @param P - A univariate polynomial ring
   * @param n - Number of coefficients to be sampled
   * @param lower_bound - Lower bound for coefficients (inclusive)
   * @param upper_bound - Upper bound for coefficients (inclusive)
   */
  constructor(
    P: PolynomialRing<C>,
    n: IntegerLike,
    lower_bound: IntegerLike,
    upper_bound: IntegerLike
  ) {
    this.n = toSafeNumber(toBigInt(n));
    this.P = P;

    const lb = toBigInt(lower_bound);
    const ub = toBigInt(upper_bound);

    if (lb > ub) {
      throw new SageTypeError('lower bound must be <= upper bound.');
    }

    this.lower_bound = lb;
    this.upper_bound = ub;
    this.D = new UniformSampler(lb, ub);
  }

  /**
   * Return a new sample polynomial.
   *
   * @returns A polynomial with random coefficients
   */
  call(): Polynomial<C> {
    const coeffs: C[] = [];
    for (let i = 0; i < this.n; i++) {
      const coeff = this.D.call();
      coeffs.push(this.P.base_ring.__call__(coeff) as C);
    }
    return this.P.__call__(coeffs);
  }

  /**
   * String representation.
   */
  repr(): string {
    return `UniformPolynomialSampler(${this.n}, ${this.lower_bound}, ${this.upper_bound})`;
  }

  toString(): string {
    return this.repr();
  }
}

/**
 * Discrete Gaussian sampler for polynomials.
 *
 * Samples univariate polynomials of degree n-1 where each coefficient
 * is sampled independently from a discrete Gaussian distribution.
 *
 * @example
 * const [R, x] = PolynomialRingConstructor(Zmod(101), 'x');
 * const sampler = new DiscreteGaussianDistributionPolynomialSampler(R, 8, 3.0);
 * const poly = sampler.call();
 *
 * @see Reference: sage/stats/distributions/discrete_gaussian_polynomial.py
 */
export class DiscreteGaussianDistributionPolynomialSampler<C extends RingElement = IntegerMod> {
  public readonly n: number;
  public readonly P: PolynomialRing<C>;
  public readonly D: DiscreteGaussianDistributionIntegerSampler;
  public readonly sigma: number;

  /**
   * Construct a sampler for univariate polynomials of degree n-1
   * where coefficients are drawn independently with standard deviation sigma.
   *
   * @param P - A univariate polynomial ring
   * @param n - Number of coefficients to be sampled
   * @param sigma - Standard deviation for each coefficient, or a DiscreteGaussianDistributionIntegerSampler
   */
  constructor(
    P: PolynomialRing<C>,
    n: number | bigint,
    sigma: number | DiscreteGaussianDistributionIntegerSampler
  ) {
    this.n = typeof n === 'bigint' ? Number(n) : n;
    this.P = P;

    if (typeof sigma === 'number') {
      this.sigma = sigma;
      this.D = new DiscreteGaussianDistributionIntegerSampler({ sigma });
    } else {
      this.D = sigma;
      this.sigma = sigma.sigma;
    }
  }

  /**
   * Return a new sample polynomial.
   *
   * @returns A polynomial with Gaussian-distributed coefficients
   */
  call(): Polynomial<C> {
    const coeffs: C[] = [];
    for (let i = 0; i < this.n; i++) {
      const coeff = this.D.call();
      coeffs.push(this.P.base_ring.__call__(coeff) as C);
    }
    return this.P.__call__(coeffs);
  }

  /**
   * String representation.
   */
  repr(): string {
    return `Discrete Gaussian sampler for polynomials of degree < ${this.n} with σ=${this.sigma.toFixed(6)} in each component`;
  }

  toString(): string {
    return this.repr();
  }
}

/**
 * A vector over Z/qZ for LWE samples.
 */
export class LWEVector {
  public readonly entries: IntegerMod[];
  public readonly ring: IntegerModRing;

  constructor(entries: IntegerMod[], ring: IntegerModRing) {
    this.entries = entries;
    this.ring = ring;
  }

  /**
   * Get the dimension of the vector.
   */
  get length(): number {
    return this.entries.length;
  }

  /**
   * Get element at index i.
   */
  get(i: number): IntegerMod {
    return this.entries[i]!;
  }

  /**
   * Compute the dot product with another vector.
   */
  dotProduct(other: LWEVector): IntegerMod {
    if (this.length !== other.length) {
      throw new ValueError('Vectors must have the same length for dot product');
    }

    let result = this.ring.zero();
    for (let i = 0; i < this.length; i++) {
      result = result.add(this.entries[i]!.mul(other.entries[i]!));
    }
    return result;
  }

  /**
   * Return entries as a tuple of bigints.
   */
  toTuple(): bigint[] {
    return this.entries.map((e) => e.value);
  }

  /**
   * String representation.
   */
  toString(): string {
    return `(${this.entries.map((e) => e.value).join(', ')})`;
  }
}

/**
 * Generate a random vector over Z/qZ.
 */
function randomVector(ring: IntegerModRing, n: number): LWEVector {
  const entries: IntegerMod[] = [];
  for (let i = 0; i < n; i++) {
    entries.push(ring.random_element());
  }
  return new LWEVector(entries, ring);
}

/**
 * Generate a vector from a distribution.
 */
function vectorFromDistribution(
  ring: IntegerModRing,
  n: number,
  D: DistributionSampler
): LWEVector {
  const entries: IntegerMod[] = [];
  for (let i = 0; i < n; i++) {
    entries.push(ring.__call__(D.call()));
  }
  return new LWEVector(entries, ring);
}

/**
 * Learning with Errors (LWE) oracle.
 *
 * The basic LWE oracle that generates samples of the form (a, <a,s> + e)
 * where s is a secret vector and e is drawn from a noise distribution.
 *
 * @example
 * // First, construct a noise distribution with uniform errors in [-2, 2]
 * const D = new UniformSampler(-2, 2);
 *
 * // Next, construct the oracle
 * const lwe = new LWE(20, 401n, D);
 *
 * // Sample from the oracle
 * const [a, c] = lwe.call();
 *
 * @see Reference: sage/crypto/lwe.py:LWE
 */
export class LWE {
  public readonly n: number;
  public readonly m: number | null;
  public readonly K: IntegerModRing;
  public readonly D: DistributionSampler;
  public readonly secret_dist: SecretDistribution;
  protected _secret: LWEVector;
  protected _sampleCount: number = 0;

  /**
   * Construct an LWE oracle in dimension n over a ring of order q
   * with noise distribution D.
   *
   * @param n - Dimension (integer > 0)
   * @param q - Modulus typically > n (integer > 0)
   * @param D - An error distribution such as an instance of
   *            DiscreteGaussianDistributionIntegerSampler or UniformSampler
   * @param secret_dist - Distribution of the secret:
   *   - 'uniform': secret follows the uniform distribution in Z/qZ
   *   - 'noise': secret follows the noise distribution
   *   - [lb, ub]: the secret is chosen uniformly from [lb,...,ub]
   * @param m - Number of allowed samples or null if no limit (default: null)
   */
  constructor(
    n: IntegerLike,
    q: IntegerLike,
    D: DistributionSampler,
    secret_dist: SecretDistribution = 'uniform',
    m: IntegerLike | null = null
  ) {
    this.n = toSafeNumber(toBigInt(n));
    this.m = m === null ? null : toSafeNumber(toBigInt(m));
    this.K = Zmod(toBigInt(q));
    this.D = D;
    this.secret_dist = secret_dist;

    // Generate the secret based on secret_dist
    if (secret_dist === 'uniform') {
      this._secret = randomVector(this.K, this.n);
    } else if (secret_dist === 'noise') {
      this._secret = vectorFromDistribution(this.K, this.n, D);
    } else if (Array.isArray(secret_dist)) {
      const [lb, ub] = secret_dist;
      const secretSampler = new UniformSampler(lb, ub);
      this._secret = vectorFromDistribution(this.K, this.n, secretSampler);
    } else {
      throw new SageTypeError(`Parameter secret_dist=${secret_dist} not understood.`);
    }
  }

  /**
   * Return the secret vector (for testing purposes).
   * In real usage this would be private.
   */
  get secret(): LWEVector {
    return this._secret;
  }

  /**
   * Return a new LWE sample (a, c) where c = <a,s> + e.
   *
   * @returns A tuple [a, c] where a is a vector and c is a scalar
   * @throws IndexError if number of available samples is exhausted
   */
  call(): [LWEVector, IntegerMod] {
    if (this.m !== null && this._sampleCount >= this.m) {
      throw new ValueError('Number of available samples exhausted.');
    }
    this._sampleCount++;

    // Generate random vector a
    const a = randomVector(this.K, this.n);

    // Compute c = <a, s> + e
    const as = a.dotProduct(this._secret);
    const e = this.K.__call__(this.D.call());
    const c = as.add(e);

    return [a, c];
  }

  /**
   * Generate m LWE samples.
   *
   * @param count - Number of samples to generate
   * @returns Array of [a, c] pairs
   */
  samples(count: IntegerLike): Array<[LWEVector, IntegerMod]> {
    const countNum = toSafeNumber(toBigInt(count));
    const result: Array<[LWEVector, IntegerMod]> = [];
    for (let i = 0; i < countNum; i++) {
      result.push(this.call());
    }
    return result;
  }

  /**
   * String representation.
   */
  repr(): string {
    const secretDistStr =
      typeof this.secret_dist === 'string'
        ? `'${this.secret_dist}'`
        : `(${this.secret_dist[0]}, ${this.secret_dist[1]})`;
    return `LWE(${this.n}, ${this.K.modulus}, ${this.D}, ${secretDistStr}, ${this.m})`;
  }

  toString(): string {
    return this.repr();
  }
}

/**
 * LWE oracle with parameters as in Regev's paper [Reg09].
 *
 * The modulus q and the standard deviation of the noise are chosen
 * as recommended in the original LWE paper.
 *
 * In [Reg09], q is the smallest prime >= n^2, and the standard deviation
 * is approximately q / (sqrt(n) * log(n)^2 * sqrt(2*pi)).
 *
 * Since we don't have DiscreteGaussian yet, we use a uniform distribution
 * with a range that approximates the same standard deviation.
 *
 * @example
 * const regev = new Regev(20);
 * const [a, c] = regev.call();
 *
 * @see Reference: sage/crypto/lwe.py:Regev
 */
export class Regev extends LWE {
  /**
   * Construct LWE instance parameterised by security parameter n where
   * the modulus q and the stddev of the noise are chosen as in [Reg09].
   *
   * @param n - Security parameter (integer > 0)
   * @param secret_dist - Distribution of the secret (default: 'uniform')
   * @param m - Number of allowed samples or null if no limit (default: null)
   */
  constructor(
    n: IntegerLike,
    secret_dist: SecretDistribution = 'uniform',
    m: IntegerLike | null = null
  ) {
    const nBig = toBigInt(n);
    const nNum = Number(nBig);

    // q is the smallest prime >= n^2
    const q = next_prime(nBig * nBig);

    // Standard deviation sigma = q / (sqrt(n) * log2(n)^2 * sqrt(2*pi))
    // For uniform distribution U(-bound, bound), stddev ~= bound / sqrt(3)
    // So bound ~= sigma * sqrt(3)
    const sqrtN = Math.sqrt(nNum);
    const log2N = Math.log2(nNum);
    const sigma = Number(q) / (sqrtN * log2N * log2N * Math.sqrt(2 * Math.PI));

    // For uniform distribution, use range that gives similar spread
    // We use a range of roughly [-3*sigma, 3*sigma] scaled down
    const bound = Math.max(1, Math.ceil(sigma));
    const D = new UniformSampler(-bound, bound);

    super(nNum, q, D, secret_dist, m === null ? null : toBigInt(m));
  }
}

/**
 * LWE oracle with parameters as in Lindner-Peikert [LP2011].
 *
 * The modulus q and the standard deviation of the noise are chosen
 * following the Lindner-Peikert recommendations.
 *
 * @example
 * const lp = new LindnerPeikert(20);
 * const [a, c] = lp.call();
 *
 * @see Reference: sage/crypto/lwe.py:LindnerPeikert
 */
export class LindnerPeikert extends LWE {
  /**
   * Construct LWE instance parameterised by security parameter n where
   * the modulus q and the stddev of the noise is chosen as in [LP2011].
   *
   * @param n - Security parameter (integer > 0)
   * @param delta - Error probability per symbol (default: 0.01)
   * @param m - Number of allowed samples or null, in which case m=2*n+128
   *            as in [LP2011] (default: null)
   */
  constructor(n: IntegerLike, delta: number = 0.01, m: IntegerLike | null = null) {
    const nNum = toSafeNumber(toBigInt(n));
    let mNum: number | null = m === null ? null : toSafeNumber(toBigInt(m));

    if (mNum === null) {
      mNum = 2 * nNum + 128;
    }

    // Simplified parameter selection based on [LP2011]
    // Find c >= 1 such that (c * exp((1-c^2)/2))^(2n) == 2^(-40)
    // We use numerical approximation
    let c = 1.0;
    for (let iter = 0; iter < 100; iter++) {
      const f = 2 * nNum * Math.log(c) + nNum * (1 - c * c) + 40 * Math.log(2);
      if (Math.abs(f) < 1e-10) break;
      const df = (2 * nNum) / c - 2 * nNum * c;
      c = c - f / df;
      if (c < 1) c = 1;
    }

    // Upper bound on s^2/t
    const s_t_bound = (Math.sqrt(2) * Math.PI) / c / Math.sqrt(2 * nNum * Math.log(2 / delta));

    // Choose q just large enough to allow for Gaussian parameter s >= 8
    const qApprox = (256 / s_t_bound) ** 1;
    const qBits = Math.ceil(Math.log2(qApprox));
    const q = next_prime(BigInt(Math.floor(2 ** qBits)));

    // Gaussian parameter and stddev
    const s = Math.sqrt(s_t_bound * Math.floor(Number(q) / 4));
    const stddev = s / Math.sqrt(2 * Math.PI);

    // Use uniform distribution approximating the Gaussian
    const bound = Math.max(1, Math.ceil(stddev * 2));
    const D = new UniformSampler(-bound, bound);

    // Call parent with noise secret distribution as in [LP2011]
    super(nNum, q, D, 'noise', mNum);
  }
}

/**
 * LWE oracle with uniform secret with parameters as in [CGW2013].
 *
 * @example
 * const lwe = new UniformNoiseLWE(89);
 * const [a, c] = lwe.call();
 *
 * @example
 * const lwe = new UniformNoiseLWE(89, 'encrypt');
 * const [a, c] = lwe.call();
 *
 * @see Reference: sage/crypto/lwe.py:UniformNoiseLWE
 */
export class UniformNoiseLWE extends LWE {
  /**
   * Construct LWE instance parameterised by security parameter n where
   * all other parameters are chosen as in [CGW2013].
   *
   * @param n - Security parameter (integer >= 89)
   * @param instance - One of:
   *   - 'key': the LWE-instance that hides the secret key is generated
   *   - 'encrypt': the LWE-instance that hides the message is generated
   * @param m - Number of allowed samples or null, in which case m is
   *            chosen as in [CGW2013] (default: null)
   */
  constructor(
    n: IntegerLike,
    instance: UniformNoiseLWEInstance = 'key',
    m: IntegerLike | null = null
  ) {
    const nNum = toSafeNumber(toBigInt(n));
    const mNum = m === null ? null : toSafeNumber(toBigInt(m));

    if (nNum < 89) {
      throw new SageTypeError('Parameter too small');
    }

    // Parameters from [CGW2013]
    const n2 = nNum;
    const C = 4 / Math.sqrt(2 * Math.PI);
    const kk = Math.floor((n2 - 2 * Math.log2(n2) ** 2) / 5);
    const n1 = Math.floor((3 * n2 - 5 * kk) / 2);
    const ke = Math.floor((n1 - 2 * Math.log2(n1) ** 2) / 5);
    const l = Math.floor((3 * n1 - 5 * ke) / 2) - n2;
    const sk = Math.ceil((C * (n1 + n2)) ** (3 / 2));
    const se = Math.ceil((C * (n1 + n2 + l)) ** (3 / 2));

    const q = next_prime(
      BigInt(
        Math.max(
          Math.ceil((4 * sk) ** ((n1 + n2) / n1)),
          Math.ceil((4 * se) ** ((n1 + n2 + l) / (n2 + l))),
          Math.ceil(4 * (n1 + n2) * se * sk + 4 * se + 1)
        )
      )
    );

    if (kk <= 0) {
      throw new SageTypeError('Parameter too small');
    }

    let D: UniformSampler;
    let effectiveN: number;
    let effectiveM: number | null;

    if (instance === 'key') {
      D = new UniformSampler(0, sk - 1);
      effectiveM = mNum === null ? n1 : mNum;
      effectiveN = n2;
    } else if (instance === 'encrypt') {
      D = new UniformSampler(0, se - 1);
      effectiveM = mNum === null ? n2 + l : mNum;
      effectiveN = n1;
    } else {
      throw new SageTypeError(`Parameter instance=${instance} not understood.`);
    }

    super(effectiveN, q, D, 'noise', effectiveM);
  }
}

/**
 * Interface for polynomial distribution samplers (used for Ring-LWE).
 */
export interface PolynomialDistributionSampler<C extends RingElement = IntegerMod> {
  /**
   * Number of coefficients.
   */
  n: number;

  /**
   * Sample from the distribution.
   */
  call(): Polynomial<C>;
}

/**
 * Ring Learning with Errors (Ring-LWE) oracle.
 *
 * Ring-LWE operates over polynomial rings, typically cyclotomic fields,
 * which allows for more compact representations and efficient operations.
 *
 * @example
 * const K = Zmod(257n);
 * const [R, x] = PolynomialRingConstructor(K, 'x');
 * const D = new DiscreteGaussianDistributionPolynomialSampler(R, 8, 3.0);
 * const ringlwe = new RingLWE(16, 257n, D, null, 'uniform');
 * const [a, c] = ringlwe.call();
 *
 * @see Reference: sage/crypto/lwe.py:RingLWE
 */
export class RingLWE {
  public readonly N: number;
  public readonly n: number;
  public readonly m: number | null;
  public readonly K: IntegerModRing;
  public readonly D: PolynomialDistributionSampler<IntegerMod>;
  public readonly q: bigint;
  public readonly poly: Polynomial<IntegerMod>;
  public readonly R_q: QuotientRing<IntegerMod>;
  public readonly secret_dist: SecretDistribution;
  protected _secret: QuotientRingElement<IntegerMod>;
  protected _sampleCount: number = 0;

  /**
   * Construct a Ring-LWE oracle in dimension n=phi(N) over a ring of order
   * q with noise distribution D.
   *
   * @param N - Index of cyclotomic polynomial (integer > 0)
   * @param q - Modulus typically > N (integer > 0)
   * @param D - An error distribution such as DiscreteGaussianDistributionPolynomialSampler
   * @param poly - A polynomial of degree phi(N). If null, the cyclotomic polynomial
   *               is used (default: null)
   * @param secret_dist - Distribution of the secret (default: 'uniform')
   * @param m - Number of allowed samples or null if no limit (default: null)
   */
  constructor(
    N: IntegerLike,
    q: IntegerLike,
    D: PolynomialDistributionSampler<IntegerMod>,
    poly: Polynomial<IntegerMod> | null = null,
    secret_dist: SecretDistribution = 'uniform',
    m: IntegerLike | null = null
  ) {
    this.N = toSafeNumber(toBigInt(N));
    this.n = Number(euler_phi(BigInt(this.N)));
    this.m = m === null ? null : toSafeNumber(toBigInt(m));
    const qBig = toBigInt(q);
    this.q = qBig;
    this.K = Zmod(qBig);

    // Check that noise distribution has correct dimension
    if (D.n !== this.n) {
      throw new ValueError(`Noise distribution has dimensions ${D.n} != ${this.n}`);
    }

    this.D = D;

    // Create polynomial ring over K
    const polyRing = new PolynomialRing(this.K, 'x');

    // Use cyclotomic polynomial if poly is not provided
    if (poly !== null) {
      this.poly = poly;
    } else {
      this.poly = polyRing.cyclotomic_polynomial(this.N);
    }

    // Create quotient ring R_q = K[x] / <poly>
    this.R_q = new QuotientRing(polyRing, this.poly);

    // Generate the secret
    this.secret_dist = secret_dist;
    if (secret_dist === 'uniform') {
      // Uniform sampling of secret from R_q
      this._secret = this._randomElement();
    } else if (secret_dist === 'noise') {
      // Secret from noise distribution
      const secretPoly = this.D.call();
      this._secret = this.R_q.__call__(secretPoly);
    } else {
      throw new SageTypeError(`Parameter secret_dist=${secret_dist} not understood.`);
    }
  }

  /**
   * Generate a random element of R_q.
   */
  protected _randomElement(): QuotientRingElement<IntegerMod> {
    const coeffs: IntegerMod[] = [];
    for (let i = 0; i < this.n; i++) {
      coeffs.push(this.K.random_element());
    }
    const poly = this.R_q.polynomial_ring.__call__(coeffs);
    return this.R_q.__call__(poly);
  }

  /**
   * Return the secret (for testing purposes).
   */
  get secret(): QuotientRingElement<IntegerMod> {
    return this._secret;
  }

  /**
   * Return a new Ring-LWE sample.
   *
   * Returns (a, a*s + e) where a is random, s is the secret, and e is noise.
   * Both a and a*s+e are returned as vectors of coefficients.
   *
   * @returns A tuple [a, c] where both are vectors
   * @throws IndexError if number of available samples is exhausted
   */
  call(): [bigint[], bigint[]] {
    if (this.m !== null && this._sampleCount >= this.m) {
      throw new ValueError('Number of available samples exhausted.');
    }
    this._sampleCount++;

    // Generate random a in R_q
    const a = this._randomElement();

    // Compute a * s
    const as = a.mul(this._secret);

    // Sample error from noise distribution
    const e = this.D.call();

    // Compute c = a*s + e
    const eInRq = this.R_q.__call__(e);
    const c = as.add(eInRq);

    // Return as vectors of coefficients
    return [this._toVector(a), this._toVector(c)];
  }

  /**
   * Convert a quotient ring element to a vector of bigints.
   */
  protected _toVector(elem: QuotientRingElement<IntegerMod>): bigint[] {
    const result: bigint[] = [];
    for (let i = 0; i < this.n; i++) {
      const coeff = elem.lift.getCoeff(i);
      result.push(coeff.value);
    }
    return result;
  }

  /**
   * String representation.
   */
  repr(): string {
    const secretDistStr =
      typeof this.secret_dist === 'string'
        ? `'${this.secret_dist}'`
        : `(${this.secret_dist[0]}, ${this.secret_dist[1]})`;
    return `RingLWE(${this.N}, ${this.q}, ${this.D}, ${this.poly}, ${secretDistStr}, ${this.m})`;
  }

  toString(): string {
    return this.repr();
  }
}

/**
 * Ring-LWE oracle with parameters as in Lindner-Peikert [LP2011].
 *
 * @example
 * const rlp = new RingLindnerPeikert(16);
 * const [a, c] = rlp.call();
 *
 * @see Reference: sage/crypto/lwe.py:RingLindnerPeikert
 */
export class RingLindnerPeikert extends RingLWE {
  /**
   * Construct a Ring-LWE oracle in dimension n=phi(N) where the modulus q
   * and the stddev of the noise is chosen as in [LP2011].
   *
   * @param N - Index of cyclotomic polynomial (integer > 0)
   * @param delta - Error probability per symbol (default: 0.01)
   * @param m - Number of allowed samples or null, in which case 3*n is used
   *            (default: null)
   */
  constructor(N: IntegerLike, delta: number = 0.01, m: IntegerLike | null = null) {
    const NNum = toSafeNumber(toBigInt(N));
    const n = Number(euler_phi(BigInt(NNum)));
    let mNum: number | null = m === null ? null : toSafeNumber(toBigInt(m));

    if (mNum === null) {
      mNum = 3 * n;
    }

    // Find c >= 1 such that (c * exp((1-c^2)/2))^(2n) == 2^(-40)
    // i.e. c >= 1 such that 2*n*log(c) + n*(1-c^2) + 40*log(2) == 0
    // Use Newton's method to find root
    let c = 1.0;
    for (let iter = 0; iter < 100; iter++) {
      const f = 2 * n * Math.log(c) + n * (1 - c * c) + 40 * Math.log(2);
      if (Math.abs(f) < 1e-10) break;
      const df = (2 * n) / c - 2 * n * c;
      c = c - f / df;
      if (c < 1) c = 1;
    }

    // Upper bound on s^2/t
    const s_t_bound = (Math.sqrt(2) * Math.PI) / c / Math.sqrt(2 * n * Math.log(2 / delta));

    // Choose q just large enough to allow for Gaussian parameter s >= 8
    const qApprox = 256 / s_t_bound;
    const qBits = Math.ceil(Math.log2(qApprox));
    const q = next_prime(BigInt(Math.floor(2 ** qBits)));

    // Gaussian parameter as defined in [LP2011]
    const s = Math.sqrt(s_t_bound * Math.floor(Number(q) / 4));

    // Transform s into stddev
    const stddev = s / Math.sqrt(2 * Math.PI);

    // Create the coefficient ring and polynomial ring
    const K = Zmod(q);
    const polyRing = new PolynomialRing(K, 'x');

    // Create the Discrete Gaussian polynomial sampler
    const D = new DiscreteGaussianDistributionPolynomialSampler(polyRing, n, stddev);

    // Call parent constructor
    super(NNum, q, D, null, 'noise', mNum);
  }
}

/**
 * Wrapper callable to convert Ring-LWE oracles into LWE oracles by
 * disregarding the additional structure.
 *
 * This allows using Ring-LWE instances in contexts that expect standard LWE.
 * Each Ring-LWE sample generates n LWE-style samples by rotating the polynomial.
 *
 * @example
 * const K = Zmod(257n);
 * const [R, x] = PolynomialRingConstructor(K, 'x');
 * const D = new DiscreteGaussianDistributionPolynomialSampler(R, 8, 5);
 * const rlwe = new RingLWE(16, 257n, D, null, 'uniform');
 * const lwe = new RingLWEConverter(rlwe);
 * const [a, c] = lwe.call();
 *
 * @see Reference: sage/crypto/lwe.py:RingLWEConverter
 */
export class RingLWEConverter {
  public readonly ringlwe: RingLWE;
  public readonly n: number;
  private _i: number = 0;
  private _ac: [bigint[], bigint[]] | null = null;

  /**
   * Construct a converter that wraps a Ring-LWE oracle.
   *
   * @param ringlwe - An instance of RingLWE
   */
  constructor(ringlwe: RingLWE) {
    this.ringlwe = ringlwe;
    this.n = ringlwe.n;
  }

  /**
   * Return a new LWE-style sample from the Ring-LWE oracle.
   *
   * Each call rotates through the coefficients of the Ring-LWE polynomial.
   * A new Ring-LWE sample is generated every n calls.
   *
   * @returns A tuple [a, c] where a is a tuple of bigints and c is a scalar bigint
   */
  call(): [bigint[], bigint] {
    // Get a new Ring-LWE sample every n calls
    if (this._i % this.n === 0) {
      this._ac = this.ringlwe.call();
    }

    const [a, c] = this._ac!;

    // Rotate the polynomial 'a' by x^i in R_q
    // This is multiplication by x^i in the quotient ring R_q = K[x]/<Phi_N(x)>
    // For cyclotomic polynomials of index N (power of 2), x^n = -1, so
    // x^i * (a_0 + a_1*x + ... + a_{n-1}*x^{n-1}) rotates coefficients with sign changes
    const rotatedA = this._rotatePolynomial(a, this._i % this.n);

    // Get the corresponding c coefficient
    const cScalar = c[this._i % this.n]!;

    this._i++;

    return [rotatedA, cScalar];
  }

  /**
   * Rotate a polynomial's coefficient vector by multiplying by x^i in R_q.
   * For Phi_N(x) with N a power of 2, x^n = -1.
   */
  private _rotatePolynomial(coeffs: bigint[], shift: number): bigint[] {
    const n = this.n;
    const q = this.ringlwe.q;
    const result: bigint[] = new Array(n);

    for (let j = 0; j < n; j++) {
      const newPos = (j + shift) % n;
      // When we wrap around (j + shift >= n), multiply by -1 due to x^n = -1
      const wraps = Math.floor((j + shift) / n);
      if (wraps % 2 === 0) {
        result[newPos] = coeffs[j]!;
      } else {
        // Negate in Z_q
        result[newPos] = coeffs[j] === 0n ? 0n : q - coeffs[j]!;
      }
    }

    return result;
  }

  /**
   * String representation.
   */
  repr(): string {
    return `RingLWEConverter(${this.ringlwe.repr()})`;
  }

  toString(): string {
    return this.repr();
  }
}

/**
 * Options for the samples function.
 */
export interface SamplesOptions {
  /**
   * Seed to be used for generation or null if no specific seed
   * shall be set.
   */
  seed?: number | null;

  /**
   * Use balance_sample to return balanced representations of
   * finite field elements (integers between -q//2 and q//2 instead
   * of 0 and q-1).
   * @default false
   */
  balanced?: boolean;
}

/**
 * Type for LWE oracle specification: either a class, instance, or name string.
 */
export type LWEOracle =
  | typeof LWE
  | typeof Regev
  | typeof LindnerPeikert
  | typeof UniformNoiseLWE
  | typeof RingLWE
  | typeof RingLindnerPeikert
  | LWE
  | RingLWE
  | 'Regev'
  | 'LindnerPeikert'
  | 'UniformNoiseLWE'
  | 'RingLindnerPeikert';

/**
 * Map of oracle names to classes.
 */
const ORACLE_CLASSES: Record<string, new (n: number, ...args: unknown[]) => LWE> = {
  Regev: Regev as unknown as new (n: number, ...args: unknown[]) => LWE,
  LindnerPeikert: LindnerPeikert as unknown as new (n: number, ...args: unknown[]) => LWE,
  UniformNoiseLWE: UniformNoiseLWE as unknown as new (n: number, ...args: unknown[]) => LWE,
};

/**
 * Return m LWE samples.
 *
 * @param m - The number of samples (integer > 0)
 * @param n - The security parameter (integer > 0)
 * @param lwe - Either:
 *   - A subclass of LWE such as Regev or LindnerPeikert
 *   - An instance of LWE or any subclass
 *   - The name of any such class (e.g., "Regev", "LindnerPeikert")
 * @param options - Additional options including seed and balanced flag
 * @returns Array of LWE samples
 *
 * @example
 * const S = samples(30, 20, 'Regev');
 *
 * @example
 * const S = samples(30, 20, LindnerPeikert);
 *
 * @example
 * const lwe = new LindnerPeikert(20);
 * const S = samples(30, 20, lwe);
 *
 * @example
 * // With balanced representation
 * const S = samples(30, 20, 'Regev', { balanced: true, seed: 1337 });
 *
 * @see Reference: sage/crypto/lwe.py:samples
 */
export function samples(
  m: IntegerLike,
  n: IntegerLike,
  lwe: LWEOracle,
  options: SamplesOptions = {}
): Array<[LWEVector | bigint[], IntegerMod | bigint]> {
  const { balanced = false } = options;
  const mNum = toSafeNumber(toBigInt(m));
  const nNum = toSafeNumber(toBigInt(n));

  let oracle: LWE;

  if (typeof lwe === 'string') {
    // Look up by name
    const OracleClass = ORACLE_CLASSES[lwe];
    if (!OracleClass) {
      throw new ValueError(`Unknown LWE oracle: ${lwe}`);
    }
    oracle = new OracleClass(nNum);
  } else if (lwe instanceof LWE) {
    // Already an instance
    if (lwe.n !== nNum) {
      throw new ValueError(
        `Passed LWE instance has n=${lwe.n}, but n=${nNum} was passed to this function.`
      );
    }
    oracle = lwe;
  } else if (typeof lwe === 'function') {
    // It's a class constructor
    oracle = new (lwe as new (n: number) => LWE)(nNum);
  } else {
    throw new ValueError('Invalid LWE oracle specification');
  }

  const result: Array<[LWEVector | bigint[], IntegerMod | bigint]> = [];

  for (let i = 0; i < mNum; i++) {
    const sample = oracle.call();
    if (balanced) {
      result.push(balance_sample(sample, oracle.K.modulus));
    } else {
      result.push(sample);
    }
  }

  return result;
}

/**
 * Given (a,c) = s return a tuple (a',c') where a' is an integer vector
 * with entries between -q//2 and q//2 and c is also within these bounds.
 *
 * This function is useful to convert between Sage's standard representation
 * of elements in Z/qZ as integers between 0 and q-1 and the usual
 * representation of such elements in lattice cryptography as integers
 * between -q//2 and q//2.
 *
 * @param s - Sample of the form [a, c] where a is a vector and c is a scalar
 * @param q - Modulus. If not provided, it is inferred from the sample's ring.
 * @returns Balanced representation of the sample
 *
 * @example
 * const balanced = balance_sample([a, c]);
 * // All elements in balanced[0] and balanced[1] are in range [-q//2, q//2]
 *
 * @see Reference: sage/crypto/lwe.py:balance_sample
 */
export function balance_sample(
  s: [LWEVector, IntegerMod],
  q?: IntegerLike | null
): [bigint[], bigint] {
  const [a, c] = s;

  // Get modulus
  const modulus = q != null ? toBigInt(q) : a.ring.modulus;
  const q2 = modulus / 2n;

  // Balance the vector
  const balancedA: bigint[] = a.entries.map((e) => {
    const v = e.value;
    return v <= q2 ? v : v - modulus;
  });

  // Balance the scalar
  const cVal = c.value;
  const balancedC = cVal <= q2 ? cVal : cVal - modulus;

  return [balancedA, balancedC];
}

export type { DiscreteGaussianDistributionIntegerSampler };

/**
 * Discrete Gaussian distribution sampler over polynomials.
 *
 * This is a placeholder type. The actual implementation should be in
 * a stats/distributions module.
 *
 * @see Reference: sage/stats/distributions/discrete_gaussian_polynomial.py
 */
export interface DiscreteGaussianDistributionPolynomialSampler {
  /**
   * Polynomial ring.
   */
  P: unknown;

  /**
   * Number of coefficients / degree bound.
   */
  n: number;

  /**
   * Standard deviation (sigma) of the Gaussian distribution for each coefficient.
   */
  sigma: number;

  /**
   * Sample from the distribution.
   */
  call(): unknown;
}
