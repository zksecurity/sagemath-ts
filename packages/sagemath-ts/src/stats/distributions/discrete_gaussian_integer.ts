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
 * All four of Sage's algorithms are implemented as ports of the vendored
 * `dgs` sources, which Sage compiles into the extension module:
 *
 * - `uniform+table`    -> `dgs_disc_gauss_mp_call_uniform_table` /
 *                         `..._uniform_table_offset` (`dgs_gauss_mp.c:306-331`)
 * - `uniform+online`   -> `dgs_disc_gauss_mp_call_uniform_online` (`dgs_gauss_mp.c:333-347`)
 * - `uniform+logtable` -> `dgs_disc_gauss_mp_call_uniform_logtable` (`dgs_gauss_mp.c:350-358`)
 * - `sigma2+logtable`  -> `dgs_disc_gauss_mp_call_sigma2_logtable` (`dgs_gauss_mp.c:360-382`)
 *
 * together with the Bernoulli machinery of `dgs_bern.c` that the two
 * `logtable` variants need.
 *
 * @see Reference: sage/stats/distributions/discrete_gaussian_integer.pyx
 * @see Reference: sage/stats/distributions/dgs_gauss_mp.c, dgs_bern.c
 */

import { TypeError as SageTypeError, ValueError } from '../../errors.js';
import { type RandState, current_randstate } from '../../misc/randstate.js';
import { type IntegerLike, toBigInt, toSafeNumber } from '../../types/coercion.js';

/**
 * Algorithm choices for discrete Gaussian sampling.
 *
 * These are the four names accepted by
 * `sage.stats.distributions.discrete_gaussian_integer` (`discrete_gaussian_integer.pyx:360-373`).
 */
export type DiscreteGaussianAlgorithm =
  | 'uniform+table'
  | 'uniform+online'
  | 'uniform+logtable'
  | 'sigma2+logtable';

/**
 * Working precision.
 *
 * Sage passes `sigma` through `RealField()` when it is not already a
 * `RealNumber` (`discrete_gaussian_integer.pyx:381-386`), and `RealField()`
 * defaults to 53 bits — the precision `dgs` then uses for every `mpfr_t`.
 */
const MPFR_PREC = 53;
const TWO_POW_PREC = 2 ** MPFR_PREC;

/**
 * Draw a uniform value in [0, 1), as MPFR's `mpfr_urandomb` does at
 * {@link MPFR_PREC} bits: it fills the significand with `prec` random bits and
 * normalises, which makes the result `k / 2^prec` for a uniform
 * `k in [0, 2^prec)`.
 */
function mpfrUrandomb(rstate: RandState): number {
  return Number(rstate.random_bits(MPFR_PREC)) / TWO_POW_PREC;
}

/**
 * Round to nearest integer, ties to even.
 *
 * This is MPFR's `MPFR_RNDN`, which `dgs_disc_gauss_mp_init` uses to split the
 * center into its integral part `c_z` and its fractional remainder `c_r`
 * (`dgs_gauss_mp.c:161-165`).
 */
function roundHalfToEven(x: number): number {
  const floor = Math.floor(x);
  const frac = x - floor;
  if (frac < 0.5) return floor;
  if (frac > 0.5) return floor + 1;
  // Exactly halfway: pick the even neighbour.
  return floor % 2 === 0 ? floor : floor + 1;
}

/**
 * Number of bits in the binary representation of a nonnegative bigint, with
 * `sizeinbase(0, 2) == 1` as GMP's `mpz_sizeinbase` has it.
 */
function sizeInBase2(x: bigint): number {
  return x === 0n ? 1 : x.toString(2).length;
}

/**
 * Balanced Bernoulli distribution with a 32-bit cache.
 *
 * Port of `dgs_bern_uniform_t` (`dgs_bern.h:82-144`, `dgs_bern.c:42-61`):
 * bits are drawn 32 at a time via `mpz_urandomb` and consumed from the low end
 * of the pool.
 */
class BernUniform {
  /** Number of bits sampled in each go (`DGS_BERN_UNIFORM_DEFAULT_LENGTH`). */
  private readonly length: number;
  /** Number of bits consumed so far. */
  private count: number;
  /** The pool of random bits. */
  private pool: bigint;

  constructor(length = 32) {
    this.length = length;
    this.count = length;
    this.pool = 0n;
  }

  /** `dgs_bern_uniform_call` (`dgs_bern.h:130-144`). */
  call(rstate: RandState): number {
    if (this.count === this.length) {
      this.pool = rstate.random_bits(this.length);
      this.count = 0;
    }
    const b = Number(this.pool & 1n);
    this.pool >>= 1n;
    this.count++;
    return b;
  }

  /** `dgs_bern_uniform_flush_cache` (`dgs_bern.h:173-175`). */
  flush_cache(): void {
    this.count = this.length;
  }
}

/**
 * Bernoulli distribution returning 1 with probability `p`.
 *
 * Port of `dgs_bern_mp_t` (`dgs_bern.c:67-93`): compare an `mpfr_urandomb`
 * draw against `p`.
 */
class BernMp {
  constructor(private readonly p: number) {}

  /** `dgs_bern_mp_call` (`dgs_bern.c:79-87`). */
  call(rstate: RandState): number {
    return mpfrUrandomb(rstate) < this.p ? 1 : 0;
  }
}

/**
 * Family of Bernoulli samplers with `p = exp(-x/f)` for nonnegative integers `x`.
 *
 * Port of `dgs_bern_exp_mp_t` (`dgs_bern.c:99-158`): precompute
 * `p[i] = exp(-2^i/f)` and, for a given `x`, return 1 only if every
 * sub-sampler indexed by a set bit of `x` returns 1.
 */
class BernExpMp {
  /** Supported inputs are `0 <= x < 2^l`. */
  readonly l: number;
  private readonly B: BernMp[];

  constructor(f: number, l: number) {
    const B: BernMp[] = [];
    let tmp = -1.0 / f;
    let len = l;
    for (let i = 0; i < l; i++) {
      const p = Math.exp(tmp);
      if (p === 0) {
        // dgs stops here and leaves `p[i]` uninitialised (`dgs_bern.c:121-124`);
        // storing the exact 0 makes the sub-sampler deterministically reject,
        // which is what an underflowed probability means.
        B.push(new BernMp(0));
        len = i + 1;
        break;
      }
      B.push(new BernMp(p));
      tmp = 2 * tmp;
    }
    this.B = B;
    this.l = len;
  }

  /** `dgs_bern_exp_mp_call` (`dgs_bern.c:146-158`). */
  call(x: bigint, rstate: RandState): number {
    const size = sizeInBase2(x);
    const start = size < this.l ? size : this.l;
    for (let i = start - 1; i >= 0; i--) {
      if ((x >> BigInt(i)) & 1n) {
        if (this.B[i]!.call(rstate) === 0) {
          return 0;
        }
      }
    }
    return 1;
  }
}

/**
 * Sampler for `D_{sigma_2}` with `sigma_2 = sqrt(1/(2 log 2))`, restricted to
 * nonnegative integers.
 *
 * Port of `dgs_disc_gauss_sigma2p_t` (`dgs_gauss_mp.c:40-69`): it needs only
 * fair coin flips, no exponentials.
 */
class DiscGaussSigma2p {
  readonly B = new BernUniform();

  /** `dgs_disc_gauss_sigma2p_mp_call` (`dgs_gauss_mp.c:47-69`). */
  call(rstate: RandState): bigint {
    for (;;) {
      if (!this.B.call(rstate)) {
        return 0n;
      }
      let dobreak = false;
      for (let i = 1; ; i++) {
        for (let j = 0; j < 2 * i - 2; j++) {
          if (this.B.call(rstate)) {
            dobreak = true;
            break;
          }
        }
        if (dobreak) {
          break;
        }
        if (!this.B.call(rstate)) {
          return BigInt(i);
        }
      }
    }
  }
}

/**
 * `sigma_2 = sqrt(1 / (2 * log 2))`.
 *
 * `dgs_gauss_mp.c:273-279` computes `mpfr_log(2)` — the *natural* logarithm,
 * despite the source comment writing it as `log_2 2`.
 */
const SIGMA2 = Math.sqrt(1 / (2 * Math.log(2)));

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
   * Samples are drawn from [round(c) - ceil(sigma*tau), round(c) + ceil(sigma*tau)],
   * where round() is round-half-to-even (MPFR's `MPFR_RNDN`).
   */
  tau?: IntegerLike;

  /**
   * Algorithm to use:
   * - 'uniform+table': Precompute probability table (faster for repeated sampling)
   * - 'uniform+online': Compute probabilities on-the-fly (less memory)
   * - 'uniform+logtable': Bernoulli rejection, no `exp` per sample; integer `c` only
   * - 'sigma2+logtable': sample `k*sigma_2` then reject; adjusts sigma; integer `c` only
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
   * We use tables for sigma*tau <= table_cutoff.
   *
   * Reference: `discrete_gaussian_integer.pyx:163`.
   */
  public static readonly table_cutoff = 10 ** 6;

  /**
   * Standard deviation of the Gaussian distribution.
   *
   * For `'sigma2+logtable'` this is the *adjusted* value `k*sigma_2`, exactly
   * as Sage reports it (`discrete_gaussian_integer.pyx:389-391`).
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
   *
   * For `'sigma2+logtable'` this describes `dgs`'s `upper_bound`, which that
   * algorithm uses only to size its Bernoulli table: its output is *not*
   * truncated to this window (`dgs_gauss_mp.c:360-382`).
   */
  public readonly lowerBound: bigint;

  /**
   * Upper bound of the sampling range (inclusive).
   *
   * @see {@link lowerBound} for the `'sigma2+logtable'` caveat.
   */
  public readonly upperBound: bigint;

  /** Integral part of the center, `mpfr_get_z(c, MPFR_RNDN)` (`dgs_gauss_mp.c:163`). */
  private readonly cZ: number;

  /** Fractional remainder of the center, `c - c_z` (`dgs_gauss_mp.c:165`). */
  private readonly cR: number;

  /** `upper_bound = ceil(sigma*tau + 1)` (`dgs_gauss_mp.c:120-122`). */
  private readonly upper_bound: bigint;

  /** `upper_bound - 1` (`dgs_gauss_mp.c:123`). */
  private readonly upper_bound_minus_one: bigint;

  /** `2*upper_bound - 1` (`dgs_gauss_mp.c:124-125`). */
  private readonly two_upper_bound_minus_one: bigint;

  /**
   * Precomputed table of `rho` values for `'uniform+table'`.
   *
   * For integer `c` this is `rho[x] = exp(-x^2/(2 sigma^2))` for
   * `0 <= x < upper_bound` with `rho[0]` halved (`dgs_gauss_mp.c:203-214`).
   * For non-integer `c` it is `rho[x + absmax] = exp(-(x - c_r)^2/(2 sigma^2))`
   * for `-absmax <= x <= absmax` (`dgs_gauss_mp.c:229-240`).
   */
  private readonly rhoTable: Float64Array | null;

  /** Bit cache used for the sign flip (`dgs_gauss_mp.c:188`, `:292`). */
  private readonly B: BernUniform | null;

  /** Bernoulli-exp family used by both `logtable` algorithms. */
  private readonly Bexp: BernExpMp | null;

  /** `D_{sigma_2}` sub-sampler used by `'sigma2+logtable'`. */
  private readonly D2: DiscGaussSigma2p | null;

  /** `k` with `sigma = k*sigma_2` (`dgs_gauss_mp.c:280-282`). */
  private readonly k: bigint;

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
    if (typeof options.sigma !== 'number' || !Number.isFinite(options.sigma)) {
      throw new SageTypeError(`sigma must be a finite number, got ${options.sigma}`);
    }
    if (options.sigma <= 0) {
      // Message from discrete_gaussian_integer.pyx:349.
      throw new ValueError(`sigma must be > 0.0 but got ${options.sigma.toFixed(6)}`);
    }
    let sigmaValue = options.sigma;

    // Validate and store c (default: 0)
    // Support both IntegerLike (public API) and number (internal/GPV algorithm)
    let cValue: number;
    if (options.c === undefined) {
      cValue = 0;
    } else if (typeof options.c === 'number') {
      if (!Number.isFinite(options.c)) {
        throw new SageTypeError(`c must be a finite number, got ${options.c}`);
      }
      cValue = options.c;
    } else {
      // IntegerLike (bigint or Integer) - validate safe range
      cValue = toSafeNumber(toBigInt(options.c));
    }
    this.c = cValue;

    // Validate and store tau (default: 6)
    // Support both IntegerLike (public API) and number (internal)
    let tauValue: number;
    if (options.tau === undefined) {
      tauValue = 6;
    } else if (typeof options.tau === 'number') {
      if (!Number.isFinite(options.tau)) {
        throw new SageTypeError(`tau must be a finite number, got ${options.tau}`);
      }
      tauValue = options.tau;
    } else {
      // IntegerLike (bigint or Integer) - validate safe range
      tauValue = toSafeNumber(toBigInt(options.tau));
    }
    if (tauValue < 1) {
      // Message from discrete_gaussian_integer.pyx:352.
      throw new ValueError(`tau must be >= 1 but got ${Math.trunc(tauValue)}`);
    }
    this.tau = tauValue;

    // Select the algorithm.  Sage uses a table for sigma*tau <= table_cutoff
    // (discrete_gaussian_integer.pyx:355-373).
    if (options.algorithm === undefined || options.algorithm === null) {
      this.algorithm =
        sigmaValue * this.tau <= DiscreteGaussianDistributionIntegerSampler.table_cutoff
          ? 'uniform+table'
          : 'uniform+online';
    } else if (
      options.algorithm === 'uniform+table' ||
      options.algorithm === 'uniform+online' ||
      options.algorithm === 'uniform+logtable' ||
      options.algorithm === 'sigma2+logtable'
    ) {
      this.algorithm = options.algorithm;
    } else {
      throw new ValueError(
        `Algorithm '${options.algorithm}' not supported by class 'DiscreteGaussianDistributionIntegerSampler'`
      );
    }

    // The two Bernoulli algorithms only support integer centers.  Sage raises
    // before even reaching dgs, and uses the 'uniform+logtable' wording for
    // both (discrete_gaussian_integer.pyx:365-372).
    if (this.algorithm === 'uniform+logtable' || this.algorithm === 'sigma2+logtable') {
      if (this.c % 1 !== 0) {
        throw new ValueError("algorithm 'uniform+logtable' requires c%1 == 0");
      }
    }

    // dgs splits the center into c_z = round_to_nearest_even(c) (MPFR_RNDN)
    // and the remainder c_r = c - c_z (dgs_gauss_mp.c:161-165).
    this.cZ = roundHalfToEven(this.c);
    this.cR = this.c - this.cZ;

    // 'sigma2+logtable' replaces sigma by k*sigma_2 *before* the bounds are
    // computed (dgs_gauss_mp.c:262-291).
    if (this.algorithm === 'sigma2+logtable') {
      const kNum = roundHalfToEven(sigmaValue / SIGMA2);
      if (kNum < 1) {
        // dgs would build a degenerate sampler (mpz_urandomm by 0).
        throw new ValueError(
          `sigma must be at least sigma_2 = ${SIGMA2.toFixed(6)} for algorithm 'sigma2+logtable' but got ${sigmaValue.toFixed(6)}`
        );
      }
      this.k = BigInt(kNum);
      sigmaValue = kNum * SIGMA2;
    } else {
      this.k = 0n;
    }
    this.sigma = sigmaValue;
    this.negHalfInvSigmaSq = -1.0 / (2.0 * this.sigma * this.sigma);

    // upper_bound = ceil(sigma*tau + 1) (dgs_gauss_mp.c:118-125); samples of
    // the uniform algorithms lie in [-(upper_bound-1), upper_bound-1] before
    // c_z is added.
    const halfWidth = Math.ceil(this.sigma * this.tau);
    this.upper_bound = BigInt(halfWidth) + 1n;
    this.upper_bound_minus_one = BigInt(halfWidth);
    this.two_upper_bound_minus_one = 2n * this.upper_bound - 1n;
    this.lowerBound = BigInt(this.cZ - halfWidth);
    this.upperBound = BigInt(this.cZ + halfWidth);

    // Per-algorithm precomputation, mirroring dgs_disc_gauss_mp_init.
    this.rhoTable = null;
    this.B = null;
    this.Bexp = null;
    this.D2 = null;

    if (this.algorithm === 'uniform+table') {
      this.B = new BernUniform();
      if (this.cR === 0) {
        // dgs_gauss_mp.c:191-215
        const size = toSafeNumber(this.upper_bound);
        const rho = new Float64Array(size);
        for (let x = 0; x < size; x++) {
          rho[x] = Math.exp(x * x * this.negHalfInvSigmaSq);
        }
        rho[0] = rho[0]! / 2;
        this.rhoTable = rho;
      } else {
        // dgs_gauss_mp.c:216-241
        const size = toSafeNumber(this.two_upper_bound_minus_one);
        const absmax = halfWidth;
        const rho = new Float64Array(size);
        for (let x = -absmax; x <= absmax; x++) {
          const d = x - this.cR;
          rho[x + absmax] = Math.exp(d * d * this.negHalfInvSigmaSq);
        }
        this.rhoTable = rho;
      }
    } else if (this.algorithm === 'uniform+logtable') {
      // _dgs_disc_gauss_mp_init_bexp (dgs_gauss_mp.c:129-136)
      this.Bexp = new BernExpMp(2 * this.sigma * this.sigma, 2 * sizeInBase2(this.upper_bound));
    } else if (this.algorithm === 'sigma2+logtable') {
      this.Bexp = new BernExpMp(2 * this.sigma * this.sigma, 2 * sizeInBase2(this.upper_bound));
      this.B = new BernUniform();
      this.D2 = new DiscGaussSigma2p();
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
   * Flush the internal cache of random bits.
   *
   * Port of `dgs_disc_gauss_mp_flush_cache` (`dgs_gauss.h:599-601`), exposed by
   * Sage as `_flush_cache()` (`discrete_gaussian_integer.pyx:395-441`).  Only
   * the algorithms that keep a `dgs_bern_uniform_t` have a cache; for the
   * others this is a no-op (dgs would dereference a NULL pointer).
   */
  _flush_cache(): void {
    if (this.B !== null) {
      this.B.flush_cache();
    }
  }

  /**
   * Sample from the discrete Gaussian distribution.
   *
   * Dispatches to the port of the corresponding `dgs_disc_gauss_mp_call_*`
   * function.
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
    const rstate = current_randstate();
    switch (this.algorithm) {
      case 'uniform+table':
        return this.cR === 0
          ? this._call_uniform_table(rstate)
          : this._call_uniform_table_offset(rstate);
      case 'uniform+online':
        return this._call_uniform_online(rstate);
      case 'uniform+logtable':
        return this._call_uniform_logtable(rstate);
      case 'sigma2+logtable':
        return this._call_sigma2_logtable(rstate);
    }
  }

  /**
   * `dgs_disc_gauss_mp_call_uniform_table` (`dgs_gauss_mp.c:306-318`).
   *
   * Integer center: draw `x` uniformly in `[0, upper_bound)`, accept with
   * probability `rho[x]` (with `rho[0]` halved so that `0` is not
   * double-counted by the sign flip), then flip the sign.
   */
  private _call_uniform_table(rstate: RandState): bigint {
    const rho = this.rhoTable!;
    let x: number;
    do {
      x = Number(rstate.random_below(this.upper_bound));
    } while (mpfrUrandomb(rstate) >= rho[x]!);
    let rop = BigInt(x);
    if (this.B!.call(rstate)) {
      rop = -rop;
    }
    return rop + BigInt(this.cZ);
  }

  /**
   * `dgs_disc_gauss_mp_call_uniform_table_offset` (`dgs_gauss_mp.c:320-331`).
   *
   * Non-integer center: the symmetry is lost, so the table covers the whole
   * window and there is no sign flip.
   */
  private _call_uniform_table_offset(rstate: RandState): bigint {
    const rho = this.rhoTable!;
    let x: number;
    do {
      x = Number(rstate.random_below(this.two_upper_bound_minus_one));
    } while (mpfrUrandomb(rstate) >= rho[x]!);
    return BigInt(x) - this.upper_bound_minus_one + BigInt(this.cZ);
  }

  /**
   * `dgs_disc_gauss_mp_call_uniform_online` (`dgs_gauss_mp.c:333-347`).
   */
  private _call_uniform_online(rstate: RandState): bigint {
    let x: bigint;
    for (;;) {
      x = rstate.random_below(this.two_upper_bound_minus_one) - this.upper_bound_minus_one;
      const d = Number(x) - this.cR;
      const z = Math.exp(d * d * this.negHalfInvSigmaSq);
      if (mpfrUrandomb(rstate) < z) {
        break;
      }
    }
    return x + BigInt(this.cZ);
  }

  /**
   * `dgs_disc_gauss_mp_call_uniform_logtable` (`dgs_gauss_mp.c:350-358`).
   *
   * Draw `x` uniformly from the window and accept with probability
   * `exp(-x^2/(2 sigma^2))`, evaluated with logarithmically many Bernoulli
   * trials instead of a call to `exp`.
   */
  private _call_uniform_logtable(rstate: RandState): bigint {
    let x: bigint;
    do {
      x = rstate.random_below(this.two_upper_bound_minus_one) - this.upper_bound_minus_one;
    } while (this.Bexp!.call(x * x, rstate) === 0);
    return x + BigInt(this.cZ);
  }

  /**
   * `dgs_disc_gauss_mp_call_sigma2_logtable` (`dgs_gauss_mp.c:360-382`).
   *
   * Sample `x` from the easily samplable `D_{sigma_2}`, `y` uniformly from
   * `[0, k)`, and accept `k*x + y` with probability
   * `exp(-((2kx + y)y)/(2 sigma^2))`, again through Bernoulli trials.
   */
  private _call_sigma2_logtable(rstate: RandState): bigint {
    let rop: bigint;
    for (;;) {
      let y = 0n;
      let x = 0n;
      for (;;) {
        x = this.D2!.call(rstate);
        y = rstate.random_below(this.k);
        // x2 = (2*k*x + y) * y
        const x2 = (2n * this.k * x + y) * y;
        if (this.Bexp!.call(x2, rstate) !== 0) {
          break;
        }
      }
      rop = this.k * x + y;
      if (rop === 0n) {
        if (this.B!.call(rstate)) {
          break;
        }
      } else {
        break;
      }
    }
    if (this.B!.call(rstate)) {
      rop = -rop;
    }
    return rop + BigInt(this.cZ);
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
      // Safe: range validated in constructor
      const xNum = toSafeNumber(x);
      sum += xNum * this._rho(xNum);
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
      // Safe: range validated in constructor
      const xNum = toSafeNumber(x);
      const diff = xNum - mu;
      sum += diff * diff * this._rho(xNum);
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
   *
   * Sage's `_repr_` (`discrete_gaussian_integer.pyx:487-497`):
   * `f"Discrete Gaussian sampler over the Integers with sigma = {self.sigma:.6f} and c = {self.c:.6f}"`.
   */
  repr(): string {
    return `Discrete Gaussian sampler over the Integers with sigma = ${this.sigma.toFixed(6)} and c = ${this.c.toFixed(6)}`;
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
  withOptions(
    options: Partial<DiscreteGaussianOptions>
  ): DiscreteGaussianDistributionIntegerSampler {
    const c = options.c === undefined ? this.c : toSafeNumber(toBigInt(options.c));
    const tau = options.tau === undefined ? this.tau : toSafeNumber(toBigInt(options.tau));
    return new DiscreteGaussianDistributionIntegerSampler({
      sigma: options.sigma ?? this.sigma,
      c,
      tau,
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
  c: IntegerLike = 0n,
  tau: IntegerLike = 6n,
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
    if (p2 === 0) return Number.POSITIVE_INFINITY;

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
