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
import { current_randstate, set_random_seed } from '../misc/randstate.js';
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
 * Numerically find a root of `f` in the bracketing interval `[a, b]`.
 *
 * Mirrors `sage.numerical.optimize.find_root(f, a, b)`, which requires a sign
 * change on the interval and returns the root to full double precision.  We use
 * bisection, which is derivative-free and cannot fail on a valid bracket (a
 * Newton iteration started at an endpoint with `f'(x) = 0` diverges).
 *
 * @param f - Continuous function with `f(a)` and `f(b)` of opposite sign
 * @param a - Left endpoint of the bracket
 * @param b - Right endpoint of the bracket
 * @returns A point where `f` changes sign, to double precision
 */
function find_root(f: (x: number) => number, a: number, b: number): number {
  let lo = a;
  let hi = b;
  let flo = f(lo);
  const fhi = f(hi);

  if (flo === 0) return lo;
  if (fhi === 0) return hi;
  if (flo * fhi > 0) {
    throw new ValueError(`f appears to have no zero on the interval [${a}, ${b}]`);
  }

  // 200 halvings takes the bracket well below double precision.
  for (let i = 0; i < 200; i++) {
    const mid = (lo + hi) / 2;
    if (mid === lo || mid === hi) break;
    const fmid = f(mid);
    if (fmid === 0) return mid;
    if (flo * fmid < 0) {
      hi = mid;
    } else {
      lo = mid;
      flo = fmid;
    }
  }
  return (lo + hi) / 2;
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
    // Sage prints the Python `None` for an unbounded sample count.
    const mStr = this.m === null ? 'None' : String(this.m);
    return `LWE(${this.n}, ${this.K.modulus}, ${this.D}, ${secretDistStr}, ${mStr})`;
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
 * is q / (sqrt(n) * log(n, 2)^2 * sqrt(2*pi)).
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

    // Sage (lwe.py:397-400):
    //   q = ZZ(next_prime(n**2))
    //   s = RR(1/(RR(n).sqrt() * log(n, 2)**2) * q)
    //   D = DiscreteGaussianDistributionIntegerSampler(s/sqrt(2*pi.n()), q)
    const q = next_prime(nBig * nBig);

    const s = (1 / (Math.sqrt(nNum) * Math.log2(nNum) ** 2)) * Number(q);
    const D = new DiscreteGaussianDistributionIntegerSampler({
      sigma: s / Math.sqrt(2 * Math.PI),
      c: Number(q),
    });

    super(BigInt(nNum), q, D, secret_dist, m === null ? null : toBigInt(m));
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

    // Sage (lwe.py:438-448):
    //   c = find_root(2*n*log(c) + n*(1-c**2) + 40*log(2) == 0, 1, 10)
    // The function is positive at c = 1 (40*log(2)) and negative at c = 10, so
    // it is bracketed on [1, 10]; a Newton step from c = 1 is degenerate
    // (f'(1) = 0), hence the bracketed solver.
    const f = (x: number) => 2 * nNum * Math.log(x) + nNum * (1 - x * x) + 40 * Math.log(2);
    const c = find_root(f, 1, 10);

    // Upper bound on s**2/t
    const s_t_bound = (Math.SQRT2 * Math.PI) / c / Math.sqrt(2 * nNum * Math.log(2 / delta));

    // "choose q just large enough to allow for a Gaussian parameter s >= 8"
    //   q = next_prime(floor(2**round(log(256 / s_t_bound, 2))))
    const q = next_prime(BigInt(Math.floor(2 ** Math.round(Math.log2(256 / s_t_bound)))));

    // Gaussian parameter as defined in [LP2011], transformed into a stddev
    const s = Math.sqrt(s_t_bound * Math.floor(Number(q) / 4));
    const stddev = s / Math.sqrt(2 * Math.PI);
    const D = new DiscreteGaussianDistributionIntegerSampler({ sigma: stddev });

    // Call parent with noise secret distribution as in [LP2011]
    super(BigInt(nNum), q, D, 'noise', mNum === null ? null : BigInt(mNum));
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
      D = new UniformSampler(0n, BigInt(sk - 1));
      effectiveM = mNum === null ? n1 : mNum;
      effectiveN = n2;
    } else if (instance === 'encrypt') {
      D = new UniformSampler(0n, BigInt(se - 1));
      effectiveM = mNum === null ? n2 + l : mNum;
      effectiveN = n1;
    } else {
      throw new SageTypeError(`Parameter instance=${instance} not understood.`);
    }

    super(BigInt(effectiveN), q, D, 'noise', effectiveM === null ? null : BigInt(effectiveM));
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
    const mStr = this.m === null ? 'None' : String(this.m);
    return `RingLWE(${this.N}, ${this.q}, ${this.D}, ${this.poly}, ${secretDistStr}, ${mStr})`;
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

    // Sage (lwe.py:641-648): find c >= 1 such that
    // 2*n*log(c) + n*(1-c^2) + 40*log(2) == 0, via find_root on [1, 10].
    const f = (x: number) => 2 * n * Math.log(x) + n * (1 - x * x) + 40 * Math.log(2);
    const c = find_root(f, 1, 10);

    // Upper bound on s**2/t
    const s_t_bound = (Math.SQRT2 * Math.PI) / c / Math.sqrt(2 * n * Math.log(2 / delta));

    // "choose q just large enough to allow for a Gaussian parameter s >= 8"
    //   q = next_prime(floor(2**round(log(256 / s_t_bound, 2))))
    const q = next_prime(BigInt(Math.floor(2 ** Math.round(Math.log2(256 / s_t_bound)))));

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
    super(BigInt(NNum), q, D, null, 'noise', mNum === null ? null : BigInt(mNum));
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

    // Sage: r = vector((x**(self._i % self.n) * R_q(a.list())).list()), c[...]
    // i.e. the multiplication happens in R_q = K[x]/<Phi_N(x)>, which is *not*
    // a signed rotation unless Phi_N(x) = x^n + 1.
    const rotatedA = this._rotatePolynomial(a, this._i % this.n);

    // Get the corresponding c coefficient
    const cScalar = c[this._i % this.n]!;

    this._i++;

    return [rotatedA, cScalar];
  }

  /**
   * Multiply the polynomial with coefficient vector `coeffs` by x^shift in
   * `R_q = K[x]/<Phi_N(x)>` and return the resulting coefficient vector.
   */
  private _rotatePolynomial(coeffs: bigint[], shift: number): bigint[] {
    const R_q = this.ringlwe.R_q;
    const K = this.ringlwe.K;

    // R_q(a.list())
    let elem = R_q.__call__(coeffs.map((v) => K.__call__(v)));

    // x^shift * a, computed in R_q
    const x = R_q.gen();
    for (let i = 0; i < shift; i++) {
      elem = elem.mul(x);
    }

    const result: bigint[] = new Array(this.n);
    for (let j = 0; j < this.n; j++) {
      result[j] = elem.lift.getCoeff(j).value;
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
export interface OracleKeywords {
  /**
   * Secret distribution, passed through to `Regev`.
   */
  secret_dist?: SecretDistribution;

  /**
   * Error probability per symbol, passed through to `LindnerPeikert` and
   * `RingLindnerPeikert`.
   */
  delta?: number;

  /**
   * Instance kind, passed through to `UniformNoiseLWE`.
   */
  instance?: UniformNoiseLWEInstance;
}

export interface SamplesOptions extends OracleKeywords {
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
 * A sample as produced by an LWE, Ring-LWE or RingLWEConverter oracle.
 *
 * The right-hand side is a scalar for LWE oracles and a coefficient vector for
 * Ring-LWE oracles.
 */
export type LWESample = [LWEVector, IntegerMod] | [bigint[], bigint] | [bigint[], bigint[]];

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
  | RingLWEConverter
  | 'Regev'
  | 'LindnerPeikert'
  | 'UniformNoiseLWE'
  | 'RingLindnerPeikert';

/**
 * Map of oracle names to classes.
 */
const ORACLE_CLASSES: Record<string, LWEOracleClass> = {
  Regev,
  LindnerPeikert,
  UniformNoiseLWE,
  RingLindnerPeikert,
};

/**
 * Any oracle class that `samples` knows how to instantiate.
 */
type LWEOracleClass =
  | typeof Regev
  | typeof LindnerPeikert
  | typeof UniformNoiseLWE
  | typeof RingLindnerPeikert
  | typeof LWE
  | typeof RingLWE;

/**
 * Instantiate an oracle class the way Sage does: `lwe(n, m=m, **kwds)`.
 *
 * The `m` argument is the third positional parameter of every parameterised
 * oracle, so it has to be routed past the class-specific second parameter.
 */
function instantiateOracle(
  cls: LWEOracleClass,
  n: bigint,
  m: bigint,
  kwds: OracleKeywords
): LWE | RingLWE {
  if (cls === Regev) {
    return new Regev(n, kwds.secret_dist ?? 'uniform', m);
  }
  if (cls === LindnerPeikert) {
    return new LindnerPeikert(n, kwds.delta ?? 0.01, m);
  }
  if (cls === UniformNoiseLWE) {
    return new UniformNoiseLWE(n, kwds.instance ?? 'key', m);
  }
  if (cls === RingLindnerPeikert) {
    return new RingLindnerPeikert(n, kwds.delta ?? 0.01, m);
  }
  throw new ValueError(
    'Only parameterised LWE oracles can be constructed from a class; pass an instance instead.'
  );
}

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
): LWESample[] {
  const { balanced = false, seed = null, ...kwds } = options;
  const mNum = toSafeNumber(toBigInt(m));
  const nNum = toSafeNumber(toBigInt(n));

  // Sage: `if seed is not None: set_random_seed(seed)`
  if (seed !== null && seed !== undefined) {
    set_random_seed(seed);
  }

  let oracle: LWE | RingLWE | RingLWEConverter;

  if (typeof lwe === 'string') {
    // Sage: `lwe = eval(lwe)`, then instantiated below
    const OracleClass = ORACLE_CLASSES[lwe];
    if (!OracleClass) {
      throw new ValueError(`Unknown LWE oracle: ${lwe}`);
    }
    oracle = instantiateOracle(OracleClass, BigInt(nNum), BigInt(mNum), kwds);
  } else if (typeof lwe === 'function') {
    // Sage: `lwe = lwe(n, m=m, **kwds)`
    oracle = instantiateOracle(lwe as LWEOracleClass, BigInt(nNum), BigInt(mNum), kwds);
  } else if (lwe instanceof LWE || lwe instanceof RingLWE || lwe instanceof RingLWEConverter) {
    // Sage only checks `lwe.n != n` here, so Ring-LWE oracles are accepted too.
    if (lwe.n !== nNum) {
      throw new ValueError(
        `Passed LWE instance has n=${lwe.n}, but n=${nNum} was passed to this function.`
      );
    }
    oracle = lwe;
  } else {
    throw new ValueError('Invalid LWE oracle specification');
  }

  const modulus =
    oracle instanceof LWE
      ? oracle.K.modulus
      : oracle instanceof RingLWE
        ? oracle.q
        : oracle.ringlwe.q;

  const result: LWESample[] = [];

  for (let i = 0; i < mNum; i++) {
    const sample = oracle.call() as LWESample;
    result.push(balanced ? balance_sample(sample, modulus) : sample);
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
  s: LWESample,
  q?: IntegerLike | null
): [bigint[], bigint] | [bigint[], bigint[]] {
  const [a, c] = s;

  // Get modulus
  let modulus: bigint;
  if (q != null) {
    modulus = toBigInt(q);
  } else if (a instanceof LWEVector) {
    modulus = a.ring.modulus;
  } else {
    throw new ValueError('modulus q is required to balance an integer sample');
  }
  const q2 = modulus / 2n;

  const balance = (v: bigint): bigint => {
    const r = ((v % modulus) + modulus) % modulus;
    return r <= q2 ? r : r - modulus;
  };

  // Balance the vector
  const aValues = a instanceof LWEVector ? a.entries.map((e) => e.value) : a;
  const balancedA = aValues.map(balance);

  // Sage returns a vector when `c` is a vector (Ring-LWE) and a scalar
  // otherwise.
  if (Array.isArray(c)) {
    return [balancedA, c.map(balance)];
  }
  const cVal = typeof c === 'bigint' ? c : c.value;
  return [balancedA, balance(cVal)];
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
