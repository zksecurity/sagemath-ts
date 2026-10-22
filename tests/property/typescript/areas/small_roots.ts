/**
 * sagemath-ts side of the `small_roots` property-test area.
 *
 * Covers Coppersmith's `small_roots` (`sage/rings/polynomial/polynomial_modn_dense_ntl.pyx`)
 * and the `Integer ** RealNumber` power that sets its default bound `X`.
 *
 * Cases: tests/property/cases/small_roots.cases.json
 * SageMath counterpart: tests/property/python/areas/small_roots.py
 *
 * `beta` and `epsilon` are indices into the tables below, which must match the
 * SageMath counterpart.
 */

import {
  mpfr_init2,
  mpfr_pow,
  mpfr_set_d,
  mpfr_set_z,
} from '../../../../packages/mpfr-ts/src/index.js';
import { GF, Zmod, next_prime } from '../../../../packages/sagemath-ts/src/index.js';
import {
  type CoefficientRing,
  type ModnCoefficient,
  type Polynomial,
  PolynomialRing,
  type SmallRootsOptions,
} from '../../../../packages/sagemath-ts/src/rings/polynomial/index.js';

/** `beta` values. The last two are rejected by `small_roots`. */
const BETAS = [1.0, 0.5, 0.25, 0.75, 0.4, 0.9, 0.0, 1.5];

/** `epsilon` values; `null` selects the default `beta / 8`. */
const EPSILONS: (number | null)[] = [null, 0.1, 0.05, -0.1, 0.0, 0.02];

type ModnPolynomial = Polynomial<ModnCoefficient>;

function fmt(roots: ModnCoefficient[]): string {
  return `[${roots.map(String).join(', ')}]`;
}

function options(
  xMode: bigint,
  x: bigint,
  betaIndex: bigint,
  epsilonIndex: bigint
): SmallRootsOptions {
  const result: SmallRootsOptions = { beta: BETAS[Number(betaIndex)]! };
  if (xMode) result.X = x;
  const epsilon = EPSILONS[Number(epsilonIndex)]!;
  if (epsilon !== null) result.epsilon = epsilon;
  return result;
}

function ring(base: unknown): PolynomialRing<ModnCoefficient> {
  return new PolynomialRing(base as CoefficientRing<ModnCoefficient>, 'x');
}

/** `prod(x - r) * (x^k + sum_{i<k} c^(i+1) x^i)` over `R`. */
function planted(
  R: PolynomialRing<ModnCoefficient>,
  roots: bigint[],
  c: bigint,
  k: bigint
): ModnPolynomial {
  const x = R.gen();
  let f: ModnPolynomial = R.one();
  for (const r of roots) f = f.mul(x.sub(R.__call__(r)));
  const cc = R.base_ring.__call__(c);
  let g: ModnPolynomial = x.pow(k) as ModnPolynomial;
  let power = R.base_ring.one();
  for (let i = 0n; i < k; i++) {
    power = power.mul(cc);
    g = g.add(R.__call__(power).mul(x.pow(i) as ModnPolynomial));
  }
  return f.mul(g);
}

export const functions: Record<string, (...args: never[]) => unknown> = {
  small_roots_monic: (
    N: bigint,
    lower: bigint[],
    xMode: bigint,
    x: bigint,
    betaIndex: bigint,
    epsilonIndex: bigint
  ) =>
    fmt(
      ring(Zmod(N))
        .__call__([...lower, 1n])
        .small_roots(options(xMode, x, betaIndex, epsilonIndex))
    ),
  small_roots_coeffs: (
    N: bigint,
    coeffs: bigint[],
    xMode: bigint,
    x: bigint,
    betaIndex: bigint,
    epsilonIndex: bigint
  ) =>
    fmt(
      ring(Zmod(N))
        .__call__(coeffs)
        .small_roots(options(xMode, x, betaIndex, epsilonIndex))
    ),
  small_roots_planted: (
    N: bigint,
    roots: bigint[],
    c: bigint,
    k: bigint,
    betaIndex: bigint,
    epsilonIndex: bigint
  ) =>
    fmt(planted(ring(Zmod(N)), roots, c, k).small_roots(options(0n, 0n, betaIndex, epsilonIndex))),
  small_roots_field: (a: bigint, roots: bigint[], c: bigint, k: bigint) =>
    fmt(planted(ring(GF(next_prime(a))), roots, c, k).small_roots()),
  small_roots_stereotyped: (a: bigint, b: bigint, kbits: bigint, K: bigint) => {
    const N = next_prime(a) * next_prime(b);
    const R = ring(Zmod(N));
    const padding = (1n << 512n) - (1n << kbits);
    const C = (padding + K) ** 3n % N;
    const f = (R.__call__(padding).add(R.gen()).pow(3n) as ModnPolynomial).sub(R.__call__(C));
    return fmt(f.small_roots());
  },
  small_roots_factor: (a: bigint, b: bigint, hidden: bigint, noise: bigint) => {
    const q = next_prime(b);
    const R = ring(Zmod(next_prime(a) * q));
    const f = R.gen().sub(R.__call__(q + noise));
    return fmt(f.small_roots({ X: (1n << hidden) - 1n, beta: 0.5 }));
  },
  small_roots_method_error: (kind: bigint) => {
    const base = [GF(2n), Zmod(2n), GF(4n, 'a'), GF(9n, 'a')][Number(kind)];
    return fmt(ring(base).gen().small_roots());
  },
  integer_pow_rr: (N: bigint, yNum: bigint, yDen: bigint) => {
    const x = mpfr_init2(53);
    const y = mpfr_init2(53);
    const z = mpfr_init2(53);
    mpfr_set_z(x, N);
    mpfr_set_d(y, Number(yNum) / Number(yDen));
    mpfr_pow(z, x, y);
    if (z.kind === 'zero') return '0';
    const shift = BigInt(z.exponent - z.precision);
    let num = BigInt(z.sign) * z.mantissa;
    let den = 1n;
    if (shift >= 0n) num <<= shift;
    else den <<= -shift;
    while ((num & 1n) === 0n && den > 1n) {
      num >>= 1n;
      den >>= 1n;
    }
    return den === 1n ? String(num) : `${num}/${den}`;
  },
};
