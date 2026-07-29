/**
 * sagemath-ts side of the `rand_stats` property-test area.
 *
 * Cases: tests/property/cases/rand_stats.cases.json
 * SageMath counterpart: tests/property/python/areas/rand_stats.py
 *
 * Scope
 * -----
 * - `sage.misc.randstate`        -> `src/misc/randstate.ts`
 * - `sage.stats.distributions.*` -> `src/stats/distributions/*`
 * - the `sage.rings.integer_ring` entry points that draw from randstate
 *   (`ZZ.random_element`), because those are what make a randstate defect
 *   observable in user code.
 *
 * Every function here seeds the generator explicitly and emits the whole draw
 * sequence: the port's RNG once had period 2, which a distribution-shape test
 * would happily pass.
 *
 * The formatting helpers at the top MUST stay byte-identical to the `_ints` /
 * `_vecs` / `_f` helpers in the Python module.
 */

import {
  current_randstate,
  initial_seed,
  random as sage_random,
  set_random_seed,
} from '../../../../packages/sagemath-ts/src/misc/randstate.js';
import { ZZ } from '../../../../packages/sagemath-ts/src/rings/integer_ring.js';
import { Rational } from '../../../../packages/sagemath-ts/src/rings/rational.js';
import {
  DiscreteGaussianDistributionIntegerSampler,
  type DiscreteGaussianOptionsInternal,
} from '../../../../packages/sagemath-ts/src/stats/distributions/discrete_gaussian_integer.js';
import {
  DiscreteGaussianDistributionLatticeSampler,
  DiscreteGaussianDistributionPolynomialSampler,
  RealNumberMP,
} from '../../../../packages/sagemath-ts/src/stats/distributions/discrete_gaussian_lattice.js';

const DGI = DiscreteGaussianDistributionIntegerSampler;
const DGL = DiscreteGaussianDistributionLatticeSampler;

// ---------------------------------------------------------------------------
// Formatting helpers — byte-identical to the Python side
// ---------------------------------------------------------------------------

/** `[1, 2, 3]` */
function ints(xs: readonly bigint[]): string {
  return `[${xs.map((x) => x.toString()).join(', ')}]`;
}

/** `[(1, 2), (3, 4)]` */
function vecs(vs: readonly (readonly bigint[])[]): string {
  return `[${vs.map((v) => `(${v.map((x) => x.toString()).join(', ')})`).join(', ')}]`;
}

/**
 * Canonical 12-significant-digit form of a double, e.g. `1.11439293741e-1`.
 *
 * 12 significant digits sits below the ~15.95 an IEEE-754 double carries, so
 * it absorbs last-ulp differences between MPFR's and libm's `exp`/`sqrt`
 * while still pinning the value hard.
 */
function f(x: number): string {
  if (Number.isNaN(x)) return 'nan';
  if (!Number.isFinite(x)) return x > 0 ? 'inf' : '-inf';
  const [mantissa, exponent] = x.toExponential(11).split('e') as [string, string];
  return `${mantissa}e${Number.parseInt(exponent, 10)}`;
}

function fs(xs: readonly number[]): string {
  return `[${xs.map(f).join(', ')}]`;
}

/** A double built from two integers, identically on both sides. */
function rat(num: bigint, den: bigint): number {
  return Number(num) / Number(den);
}

type Alg = 'uniform+table' | 'uniform+online' | 'uniform+logtable' | 'sigma2+logtable';

const ALGORITHMS: Record<string, Alg | undefined> = {
  '0': undefined,
  '1': 'uniform+table',
  '2': 'uniform+online',
  '3': 'uniform+logtable',
  '4': 'sigma2+logtable',
};

function alg(code: bigint): Alg | undefined {
  const key = code.toString();
  if (!(key in ALGORITHMS)) {
    throw new Error(`unknown algorithm code ${key}`);
  }
  return ALGORITHMS[key];
}

function square(flat: readonly bigint[], n: bigint): bigint[][] {
  const size = Number(n);
  const rows: bigint[][] = [];
  for (let i = 0; i < size; i++) {
    rows.push(flat.slice(i * size, (i + 1) * size) as bigint[]);
  }
  return rows;
}

function identity(n: bigint): number[][] {
  const size = Number(n);
  return Array.from({ length: size }, (_, i) =>
    Array.from({ length: size }, (_, j) => (i === j ? 1 : 0))
  );
}

function ratVec(flat: readonly bigint[], den: bigint): Rational[] {
  return flat.map((x) => new Rational(x, den));
}

function repeat<T>(count: bigint, fn: () => T): T[] {
  const out: T[] = [];
  for (let i = 0n; i < count; i++) out.push(fn());
  return out;
}

// ===========================================================================
// sage.misc.randstate
// ===========================================================================

function rs_random(seed: bigint, count: bigint): string {
  set_random_seed(seed);
  return ints(repeat(count, () => sage_random()));
}

/**
 * Seeded at `2^19937 - 20027 + offset`.
 *
 * `2^19937 - 20027` is the modulus GMP's `randseed_mt` reduces the seed by
 * (`rand/randmts.c:126`), so `offset == 0` drives `seed1` to 0 and exercises
 * the `seed1 += 2` fixup.
 */
function rs_random_modseed(offset: bigint, count: bigint): string {
  set_random_seed(2n ** 19937n - 20027n + offset);
  return ints(repeat(count, () => sage_random()));
}

function rs_initial_seed(seed: bigint): string {
  set_random_seed(seed);
  for (let i = 0; i < 3; i++) sage_random();
  return initial_seed().toString();
}

/**
 * `mpz_urandomm(state, n)` as reached through `ZZ.random_element(n)`.
 *
 * `ZZ._randomize_mpz` unconditionally burns one `c_random()` on its `den`
 * local before dispatching (`integer_ring.pyx:801`), so each draw is one
 * 31-bit discard followed by `mpz_urandomm`.  Composed here from the exported
 * randstate primitives so that this case stays a pure randstate test.
 */
function rs_urandomm(seed: bigint, n: bigint, count: bigint): string {
  set_random_seed(seed);
  const rstate = current_randstate();
  return ints(
    repeat(count, () => {
      rstate.c_random();
      return rstate.random_below(n);
    })
  );
}

/** `randstate.ZZ_seed()` == `ZZ.random_element(1 << 128)` (`randstate.pyx:629-642`). */
function rs_zz_seed(seed: bigint, count: bigint): string {
  set_random_seed(seed);
  const rstate = current_randstate();
  return ints(
    repeat(count, () => {
      rstate.c_random();
      return rstate.random_below(1n << 128n);
    })
  );
}

function rs_python_random(seed: bigint, count: bigint): string {
  set_random_seed(seed);
  const rnd = current_randstate().python_random();
  return fs(repeat(count, () => rnd.random()));
}

function rs_python_getrandbits(seed: bigint, k: bigint, count: bigint): string {
  set_random_seed(seed);
  const rnd = current_randstate().python_random();
  return ints(repeat(count, () => rnd.getrandbits(Number(k))));
}

function rs_python_randrange(seed: bigint, lo: bigint, hi: bigint, count: bigint): string {
  set_random_seed(seed);
  const rnd = current_randstate().python_random();
  return ints(repeat(count, () => rnd.randrange(lo, hi)));
}

function rs_python_randint(seed: bigint, a: bigint, b: bigint, count: bigint): string {
  set_random_seed(seed);
  const rnd = current_randstate().python_random();
  return ints(repeat(count, () => rnd.randint(a, b)));
}

function rs_python_normalvariate(
  seed: bigint,
  muN: bigint,
  muD: bigint,
  sigmaN: bigint,
  sigmaD: bigint,
  count: bigint
): string {
  set_random_seed(seed);
  const rnd = current_randstate().python_random();
  const mu = rat(muN, muD);
  const sigma = rat(sigmaN, sigmaD);
  return fs(repeat(count, () => rnd.normalvariate(mu, sigma)));
}

/**
 * `python_random()` first, then `python_random(seed=...)`.
 *
 * Upstream returns the *cached* object and ignores the new seed
 * (`randstate.pyx:617-620`: the `type(self._python_random) is cls` early
 * return happens before `seed` is looked at), so the second stream continues
 * the first.
 */
function rs_python_random_reseed(seed: bigint, explicitSeed: bigint, count: bigint): string {
  set_random_seed(seed);
  const rstate = current_randstate();
  const first = repeat(count, () => rstate.python_random().random());
  const second = repeat(count, () => rstate.python_random(explicitSeed).random());
  return `${fs(first)} ${fs(second)}`;
}

/** `python_random(seed=...)` as the *first* call on a fresh randstate. */
function rs_python_random_fresh_seed(seed: bigint, explicitSeed: bigint, count: bigint): string {
  set_random_seed(seed);
  const rnd = current_randstate().python_random(explicitSeed);
  return fs(repeat(count, () => rnd.random()));
}

// ===========================================================================
// sage.rings.integer_ring — the randstate consumers
// ===========================================================================

function zz_random_element_1n(seed: bigint, count: bigint): string {
  set_random_seed(seed);
  return ints(repeat(count, () => ZZ.random_element()));
}

function zz_random_element_range(seed: bigint, lo: bigint, hi: bigint, count: bigint): string {
  set_random_seed(seed);
  return ints(repeat(count, () => ZZ.random_element(lo, hi)));
}

function zz_random_element_rrandomb(seed: bigint, bits: bigint, count: bigint): string {
  set_random_seed(seed);
  return ints(repeat(count, () => ZZ.random_element(bits, undefined, 'mpz_rrandomb')));
}

function zz_random_element_gaussian(
  seed: bigint,
  sigmaN: bigint,
  sigmaD: bigint,
  count: bigint
): string {
  set_random_seed(seed);
  const sigma = rat(sigmaN, sigmaD);
  return ints(repeat(count, () => ZZ.random_element(sigma, undefined, 'gaussian')));
}

/**
 * Argument validation of `ZZ.random_element` (`integer_ring.pyx:768-781`).
 *
 * `nargs` is 0, 1 or 2 and says how many of `x`/`y` are actually passed.  Only
 * the *error surface* is compared here.
 */
function zz_error(nargs: bigint, x: bigint, y: bigint, dist: bigint): string {
  const distribution = (
    ['uniform', 'uniform', '1/n', 'mpz_rrandomb', 'gaussian', 'bogus'] as const
  )[Number(dist)];
  const passDistribution = dist === 0n ? undefined : distribution;
  const n = Number(nargs);
  set_random_seed(1n);
  try {
    if (n === 0) {
      ZZ.random_element(undefined, undefined, passDistribution as never);
    } else if (n === 1) {
      ZZ.random_element(x, undefined, passDistribution as never);
    } else {
      ZZ.random_element(x, y, passDistribution as never);
    }
  } catch (e) {
    const err = e as Error;
    return `${err.name}: ${err.message}`;
  }
  return 'ok';
}

/** `randrange` on an empty range. Both sides must raise `ValueError`. */
function rs_python_error(seed: bigint, lo: bigint, hi: bigint): string {
  set_random_seed(seed);
  const rnd = current_randstate().python_random();
  rnd.randrange(lo, hi);
  return 'ok';
}

// ===========================================================================
// sage.stats.distributions.discrete_gaussian_integer
// ===========================================================================

function makeDGI(options: DiscreteGaussianOptionsInternal) {
  return new DGI(options);
}

/**
 * Constructor validation of `DiscreteGaussianDistributionIntegerSampler`.
 *
 * Returns the exception's class and message as text so that a *message*
 * mismatch fails the case; `compare.ts` scores "both raised" as a pass
 * regardless of what they raised.
 */
function dgi_error(
  sigmaN: bigint,
  sigmaD: bigint,
  tau: bigint,
  cN: bigint,
  cD: bigint,
  algCode: bigint
): string {
  try {
    const D = makeDGI({
      sigma: rat(sigmaN, sigmaD),
      c: rat(cN, cD),
      tau: Number(tau),
      algorithm: alg(algCode),
    });
    return `ok: ${D.repr()}`;
  } catch (e) {
    const err = e as Error;
    return `${err.name}: ${err.message}`;
  }
}

function dgi_samples(
  seed: bigint,
  sigmaN: bigint,
  sigmaD: bigint,
  tau: bigint,
  algCode: bigint,
  count: bigint
): string {
  set_random_seed(seed);
  const D = makeDGI({
    sigma: rat(sigmaN, sigmaD),
    tau: Number(tau),
    algorithm: alg(algCode),
  });
  return ints(repeat(count, () => D.call()));
}

function dgi_samples_c(
  seed: bigint,
  sigmaN: bigint,
  sigmaD: bigint,
  cN: bigint,
  cD: bigint,
  algCode: bigint,
  count: bigint
): string {
  set_random_seed(seed);
  const D = makeDGI({
    sigma: rat(sigmaN, sigmaD),
    c: rat(cN, cD),
    algorithm: alg(algCode),
  });
  return ints(repeat(count, () => D.call()));
}

/**
 * Which algorithm the `algorithm=None` default picks:
 * `sigma*tau <= table_cutoff (10^6)` -> `uniform+table`
 * (`discrete_gaussian_integer.pyx:352-355`).
 */
function dgi_default_algorithm(sigmaN: bigint, sigmaD: bigint, tau: bigint): string {
  const D = makeDGI({ sigma: rat(sigmaN, sigmaD), tau: Number(tau) });
  return D.algorithm;
}

function dgi_repr(sigmaN: bigint, sigmaD: bigint, cN: bigint, cD: bigint): string {
  return makeDGI({ sigma: rat(sigmaN, sigmaD), c: rat(cN, cD) }).repr();
}

/**
 * The three streams of the `_flush_cache` doctest
 * (`discrete_gaussian_integer.pyx:404-434`).
 */
function dgi_stream_mode(
  seed: bigint,
  sigmaN: bigint,
  sigmaD: bigint,
  count: bigint,
  mode: bigint
): string {
  const m = Number(mode);
  set_random_seed(seed);
  const D = makeDGI({ sigma: rat(sigmaN, sigmaD) });
  const out: bigint[] = [];
  for (let i = 0n; i < count; i++) {
    if (m >= 1) set_random_seed(seed);
    if (m === 2) D._flush_cache();
    out.push(D.call());
  }
  return ints(out);
}

function dgi_histogram(
  seed: bigint,
  sigmaN: bigint,
  sigmaD: bigint,
  algCode: bigint,
  count: bigint,
  lo: bigint,
  hi: bigint
): string {
  set_random_seed(seed);
  const D = makeDGI({ sigma: rat(sigmaN, sigmaD), algorithm: alg(algCode) });
  const counts = new Map<bigint, bigint>();
  for (let v = lo; v <= hi; v++) counts.set(v, 0n);
  let outside = 0n;
  for (let i = 0n; i < count; i++) {
    const v = D.call();
    if (v >= lo && v <= hi) {
      counts.set(v, counts.get(v)! + 1n);
    } else {
      outside += 1n;
    }
  }
  const ordered: bigint[] = [];
  for (let v = lo; v <= hi; v++) ordered.push(counts.get(v)!);
  return `${ints(ordered)} outside=${outside}`;
}

/**
 * `(min, max, ceil(sigma*tau), every |x| within it)`.
 *
 * `dgs` bounds the support at `ceil(sigma*tau)` around `round(c)`
 * (`dgs_gauss_mp.c:161-200`).
 */
function dgi_support(
  seed: bigint,
  sigmaN: bigint,
  sigmaD: bigint,
  tau: bigint,
  algCode: bigint,
  count: bigint
): string {
  const sigma = rat(sigmaN, sigmaD);
  set_random_seed(seed);
  const D = makeDGI({ sigma, tau: Number(tau), algorithm: alg(algCode) });
  const xs = repeat(count, () => D.call());
  const bound = BigInt(Math.ceil(sigma * Number(tau)));
  let lo = xs[0]!;
  let hi = xs[0]!;
  let ok = true;
  for (const x of xs) {
    if (x < lo) lo = x;
    if (x > hi) hi = x;
    if ((x < 0n ? -x : x) > bound) ok = false;
  }
  return `(${lo}, ${hi}, ${bound}, ${ok ? 'True' : 'False'})`;
}

/** `rho_{sigma,c}(x) = exp(-(x-c)^2 / (2 sigma^2))`. */
function dgi_rho(sigmaN: bigint, sigmaD: bigint, cN: bigint, cD: bigint, x: bigint): string {
  const D = makeDGI({ sigma: rat(sigmaN, sigmaD), c: rat(cN, cD) });
  return f(D.rho(Number(x)));
}

// ===========================================================================
// sage.stats.distributions.discrete_gaussian_lattice
// ===========================================================================

/** Identity basis, centre 0 -> `_call_simple`. */
function dgl_samples(
  seed: bigint,
  n: bigint,
  sigmaN: bigint,
  sigmaD: bigint,
  count: bigint
): string {
  set_random_seed(seed);
  const D = new DGL(identity(n), { sigma: rat(sigmaN, sigmaD) });
  return vecs(repeat(count, () => D.call()));
}

/**
 * Identity basis with centre `c_flat / c_den`.
 *
 * `c_den == 1` keeps `c` in the lattice (`_call_simple`); anything else pushes
 * it off the lattice and into the GPV recursion `_call`.
 */
function dgl_samples_c(
  seed: bigint,
  n: bigint,
  sigmaN: bigint,
  sigmaD: bigint,
  cFlat: bigint[],
  cDen: bigint,
  count: bigint
): string {
  set_random_seed(seed);
  const D = new DGL(identity(n), {
    sigma: rat(sigmaN, sigmaD),
    c: ratVec(cFlat, cDen),
  });
  return vecs(repeat(count, () => D.call()));
}

/** Arbitrary (non-orthonormal) basis -> `_call`, the GPV recursion. */
function dgl_samples_basis(
  seed: bigint,
  basisFlat: bigint[],
  n: bigint,
  sigmaN: bigint,
  sigmaD: bigint,
  count: bigint
): string {
  set_random_seed(seed);
  const D = new DGL(square(basisFlat, n), { sigma: rat(sigmaN, sigmaD) });
  return vecs(repeat(count, () => D.call()));
}

function dgl_repr(n: bigint, sigmaN: bigint, sigmaD: bigint, cFlat: bigint[]): string {
  return new DGL(identity(n), {
    sigma: rat(sigmaN, sigmaD),
    c: ratVec(cFlat, 1n),
  }).repr();
}

/**
 * `__repr__` for a covariance matrix, and the `RuntimeError` upstream raises
 * when it is not positive definite (`discrete_gaussian_lattice.py:569-570`).
 */
function dgl_repr_sigma(sigmaFlat: bigint[], n: bigint): string {
  try {
    return new DGL(identity(n), {
      sigma: square(sigmaFlat, n) as unknown as number[][],
    }).repr();
  } catch (e) {
    const err = e as Error;
    return `${err.name}: ${err.message}`;
  }
}

/**
 * `compute_precision(precision, sigma)` (`discrete_gaussian_lattice.py:155`).
 *
 * `prec == -1` means `None`; `sigmaPrec == -1` means a plain integer sigma
 * (no `.precision()`, so the `AttributeError` branch).
 */
function dgl_compute_precision(prec: bigint, sigmaPrec: bigint): string {
  const p = Number(prec);
  const sp = Number(sigmaPrec);
  // `RealField(sp)(3)` == 3 with an `sp`-bit mantissa.
  const sigma: unknown = sp < 0 ? 3n : new RealNumberMP(1, 3n << BigInt(sp - 2), -(sp - 2), sp);
  return String(DGL.compute_precision(p < 0 ? undefined : p, sigma));
}

/** `_normalisation_factor_zz` on `ZZ^n` (Poisson summation via qfrep). */
function dgl_normalisation(
  n: bigint,
  sigmaN: bigint,
  sigmaD: bigint,
  tau: bigint,
  prec: bigint
): string {
  const D = new DGL(identity(n), { sigma: rat(sigmaN, sigmaD) });
  const value = D._normalisation_factor_zz(
    tau < 0n ? undefined : Number(tau),
    prec < 0n ? undefined : Number(prec)
  );
  return f(value.toNumber());
}

/** `round(_normalisation_factor_zz(prec=prec))` — exact integer readout. */
function dgl_normalisation_round(n: bigint, sigmaN: bigint, sigmaD: bigint, prec: bigint): string {
  const D = new DGL(identity(n), { sigma: rat(sigmaN, sigmaD) });
  return D._normalisation_factor_zz(undefined, Number(prec)).round().toString();
}

/**
 * `_normalisation_factor_zz` on a general basis / centre.
 *
 * Upstream raises `NotImplementedError` for a non-square basis, a non-integral
 * lattice, or a centre that is not the origin of a trivial lattice
 * (`discrete_gaussian_lattice.py:313-322`).
 */
function dgl_normalisation_basis(
  basisFlat: bigint[],
  rows: bigint,
  cols: bigint,
  basisDen: bigint,
  sigmaN: bigint,
  sigmaD: bigint,
  cFlat: bigint[],
  cDen: bigint
): string {
  const r = Number(rows);
  const c = Number(cols);
  const B: Rational[][] = [];
  for (let i = 0; i < r; i++) {
    B.push(basisFlat.slice(i * c, (i + 1) * c).map((x) => new Rational(x, basisDen)));
  }
  try {
    const D = new DGL(B, {
      sigma: rat(sigmaN, sigmaD),
      c: cFlat.length > 0 ? ratVec(cFlat, cDen) : undefined,
    });
    return f(D._normalisation_factor_zz().toNumber());
  } catch (e) {
    // Returned as text (not re-raised) so that a *message* mismatch fails the
    // case; compare.ts scores "both raised" as a pass regardless.
    const err = e as Error;
    return `${err.name}: ${err.message}`;
  }
}

/**
 * `_normalisation_factor_zz` for a covariance matrix (the approximate branch,
 * `discrete_gaussian_lattice.py:297-312`).
 */
function dgl_normalisation_nonspherical(sigmaFlat: bigint[], n: bigint, cFlat: bigint[]): string {
  const D = new DGL(identity(n), {
    sigma: square(sigmaFlat, n) as unknown as number[][],
    c: ratVec(cFlat, 1n),
  });
  return f(D._normalisation_factor_zz().toNumber());
}

/** `D.f(x)` — the (unnormalised) Gaussian weight. */
function dgl_f(
  n: bigint,
  sigmaN: bigint,
  sigmaD: bigint,
  cFlat: bigint[],
  xFlat: bigint[]
): string {
  const D = new DGL(identity(n), {
    sigma: rat(sigmaN, sigmaD),
    c: ratVec(cFlat, 1n),
  });
  return f(D.f(xFlat));
}

/** `D.f(x)` for a covariance matrix: `exp(-x Sigma^-1 x / 2)`. */
function dgl_f_sigma(sigmaFlat: bigint[], n: bigint, xFlat: bigint[]): string {
  const D = new DGL(identity(n), { sigma: square(sigmaFlat, n) as unknown as number[][] });
  return f(D.f(xFlat));
}

/** `_maximal_r()`: the largest `r` with `Sigma - r^2 Q` positive definite. */
function dgl_maximal_r(sigmaFlat: bigint[], n: bigint): string {
  const D = new DGL(identity(n), { sigma: square(sigmaFlat, n) as unknown as number[][] });
  return f(D._maximal_r());
}

/**
 * `_call_non_spherical` — Peikert's two-stage sampler.  Its offline stage
 * draws `n` `normalvariate(0, 1)` per sample from `randstate.python_random()`.
 */
function dgl_nonspherical(
  seed: bigint,
  sigmaFlat: bigint[],
  n: bigint,
  cFlat: bigint[],
  cDen: bigint,
  count: bigint
): string {
  set_random_seed(seed);
  const D = new DGL(identity(n), {
    sigma: square(sigmaFlat, n) as unknown as number[][],
    c: ratVec(cFlat, cDen),
  });
  return vecs(repeat(count, () => D.call()));
}

/** `DiscreteGaussianDistributionPolynomialSampler` coefficient streams. */
function dgl_poly_sampler(
  seed: bigint,
  n: bigint,
  sigmaN: bigint,
  sigmaD: bigint,
  count: bigint
): string {
  set_random_seed(seed);
  const sampler = new DiscreteGaussianDistributionPolynomialSampler(Number(n), {
    sigma: rat(sigmaN, sigmaD),
  });
  return vecs(repeat(count, () => sampler.sample()));
}

export const functions = {
  // randstate
  rs_random,
  rs_random_modseed,
  rs_initial_seed,
  rs_urandomm,
  rs_zz_seed,
  rs_python_random,
  rs_python_getrandbits,
  rs_python_randrange,
  rs_python_randint,
  rs_python_normalvariate,
  rs_python_random_reseed,
  rs_python_random_fresh_seed,
  // integer_ring consumers
  zz_random_element_1n,
  zz_random_element_range,
  zz_random_element_rrandomb,
  zz_random_element_gaussian,
  zz_error,
  rs_python_error,
  // discrete_gaussian_integer
  dgi_error,
  dgi_samples,
  dgi_samples_c,
  dgi_default_algorithm,
  dgi_repr,
  dgi_stream_mode,
  dgi_histogram,
  dgi_support,
  dgi_rho,
  // discrete_gaussian_lattice
  dgl_samples,
  dgl_samples_c,
  dgl_samples_basis,
  dgl_repr,
  dgl_repr_sigma,
  dgl_compute_precision,
  dgl_normalisation,
  dgl_normalisation_round,
  dgl_normalisation_basis,
  dgl_normalisation_nonspherical,
  dgl_f,
  dgl_f_sigma,
  dgl_maximal_r,
  dgl_nonspherical,
  dgl_poly_sampler,
};
