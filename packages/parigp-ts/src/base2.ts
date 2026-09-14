/** PARI base2.c:723–867, reduced resultants for integral-basis denominators. */
import { trimPolynomial } from './_polynomial_packing.js';
import { residue } from './_polynomial_division.js';
import { Flx_rem } from './Flx.js';
import { FpX_rem } from './FpX.js';
import { ZpM_echelon, zlm_echelon } from './hnf_snf.js';

/** Reduced resultant modulo p^m; the first polynomial has positive degree
 * and leading coefficient invertible modulo p (normally monic).
 */
export function ZpX_reduced_resultant(x: bigint[], y: bigint[], p: bigint, pm: bigint): bigint {
  const word = pm > 0n && pm < 1n << 64n;
  const f = trimPolynomial(x),
    n = f.length - 1;
  const rem = (a: bigint[]) =>
    word
      ? Flx_rem(
          a.map((c) => residue(c, pm)),
          f.map((c) => residue(c, pm)),
          pm
        )
      : FpX_rem(a, f, pm);
  let h = rem(y);
  // Native zero divisors raise in the remainder above. Positive-degree input
  // is required by the original Sylvester allocator.
  if (n <= 0)
    throw new RangeError('ZpX_reduced_resultant requires a polynomial of positive degree');
  const columns: bigint[][] = [];
  for (let j = 0; j < n; j++) {
    columns.push(Array.from({ length: n }, (_, i) => h[i] ?? 0n));
    if (j + 1 < n) h = rem([0n, ...h]);
  }
  const z = word ? zlm_echelon(columns, false, p, pm) : ZpM_echelon(columns, false, p, pm);
  return z?.[0]?.[0] ?? 0n;
}

/** Native doubling of p-adic precision, capped at M. */
export function ZpX_reduced_resultant_fast(f: bigint[], g: bigint[], p: bigint, M: number): bigint {
  let m =
    p === 2n
      ? 16
      : p === 3n
        ? 10
        : p === 5n
          ? 6
          : p === 7n
            ? 5
            : p === 11n || p === 13n
              ? 4
              : p < 41n
                ? 3
                : p < 257n
                  ? 2
                  : 1;
  let q: bigint | null = null;
  for (; M >= 2 * m; m *= 2) {
    q = q === null ? p ** BigInt(m) : q * q;
    const r = ZpX_reduced_resultant(f, g, p, q);
    if (r) return r;
  }
  q = p ** BigInt(M);
  return ZpX_reduced_resultant(f, g, p, q) || q;
}

import type { NfIdealData } from './base1.js';
import { zk_multable, pr_hnf } from './base4.js';
import { FpM_ker } from './alglin1.js';
import { FpM_mul } from './FpV.js';
import { inverseCoefficient } from './_polynomial_division.js';
import { ZX_resultant } from './galconj.js';
import { ZC_nfval } from './base3.js';
import { PariError } from './errors.js';

/** Native five-component prime data; tau=1 denotes an inert prime.
 * @see Deviation: Number-field ideal valuation adapters
 */
export interface NfPrimeIdeal {
  p: bigint;
  generator: bigint[];
  e: bigint;
  f: bigint;
  tau: bigint[][] | 1n;
}
const fpKernel = (A: bigint[][], p: bigint) =>
  FpM_ker([[], ...A.map((c) => [0n, ...c])], p)
    .slice(1)
    .map((c) => c.slice(1));
const fpProduct = (A: bigint[][], B: bigint[][], p: bigint) =>
  FpM_mul([[], ...A.map((c) => [0n, ...c])], [[], ...B.map((c) => [0n, ...c])], p)
    .slice(1)
    .map((c) => c.slice(1));
const columnProduct = (A: bigint[][], x: bigint[]) =>
  A[0]!.map((_, i) => x.reduce((s, c, j) => s + c * A[j]![i]!, 0n));

/** Native FpM_intersect_i, with null standing for the full ambient space. */
function mul_intersect(A: bigint[][] | null, B: bigint[][] | null, p: bigint): bigint[][] {
  if (A === null) return B!;
  if (B === null) return A;
  if (!A.length || !B.length) return [];
  return fpProduct(
    A,
    fpKernel([...A, ...B], p).map((c) => c.slice(0, A.length)),
    p
  );
}

/** Native get_LV prefix/suffix intersections of the other primes. */
function get_LV(L: bigint[][][], p: bigint, n: number): bigint[][][] {
  if (L.length === 1)
    return [
      Array.from({ length: n }, (_, j) => Array.from({ length: n }, (_, i) => (i === j ? 1n : 0n))),
    ];
  const A: (bigint[][] | null)[] = [null],
    B: (bigint[][] | null)[] = Array(L.length).fill(null);
  for (let i = 0; i < L.length - 1; i++) A[i + 1] = mul_intersect(A[i]!, L[i]!, p);
  for (let i = L.length - 1; i >= 1; i--) B[i - 1] = mul_intersect(L[i]!, B[i]!, p);
  return L.map((_, i) => mul_intersect(A[i]!, B[i]!, p));
}

/** Native uniformizer search, testing norms with the exact resultant fallback. */
function uniformizer(
  nf: NfIdealData,
  P: bigint[][],
  V: bigint[][],
  p: bigint,
  ramif: boolean
): bigint[] {
  const n = nf.multiplication.length,
    f = n - P.length,
    q = p ** BigInt(f + 1);
  const one = Array<bigint>(n).fill(0n);
  one[0] = 1n;
  // FpM_FpC_invimage: use the final dependency of [P | V | 1].
  const relations = fpKernel([...P, ...V, one], p),
    relation = relations.at(-1);
  if (!relation || !relation.at(-1)) throw new PariError('not a prime number in idealprimedec');
  const scale = inverseCoefficient(residue(-relation.at(-1)!, p), p, p < 1n << 64n);
  const center = (x: bigint) => {
    const r = x % p;
    return r > p / 2n ? r - p : r < -(p / 2n) ? r + p : r;
  };
  const u = columnProduct(
    P,
    relation.slice(0, P.length).map((x) => residue(x * scale, p))
  ).map(center);
  const isUniformizer = (a: bigint[]) => {
    const power = Array.from({ length: n }, (_, k) =>
      a.reduce((s, c, i) => s + c * nf.basis[i]![k]!, 0n)
    );
    const norm = ZX_resultant(nf.polynomial, power) / nf.basisDenominator ** BigInt(n);
    return norm % q !== 0n;
  };
  if (isUniformizer(u)) return u;
  u[0]! += u[0]! <= 0n ? p : -p;
  if (!ramif || isUniformizer(u)) return u;
  const Mv = zk_multable(
    nf.multiplication,
    u.map((x, i) => (i === 0 ? 1n : 0n) - x)
  );
  for (const c of P) {
    const x = columnProduct(Mv, c).map((v, i) => center(v + u[i]!));
    if (isUniformizer(x)) return x;
  }
  throw new PariError('not a prime number in idealprimedec');
}

/** Finalize certified prime HNFs using native primedec_end/get_pr (matrix route).
 * Input/output order is preserved. Inputs are integral prime HNFs above p.
 * @see Deviation: Number-field ideal valuation adapters
 */
export function primedec_end(nf: NfIdealData, ideals: bigint[][][], p: bigint): NfPrimeIdeal[] {
  return finalize_primes(
    nf,
    ideals.map((H) => H.filter((c, j) => c[j] === 1n).map((c) => c.slice())),
    p
  );
}
/** Native primedec_end, accepting the original mixed Kummer/finite-basis list. */
function finalize_primes(
  nf: NfIdealData,
  L: (bigint[][] | NfPrimeIdeal)[],
  p: bigint,
  limit = 0n
): NfPrimeIdeal[] {
  const n = nf.multiplication.length,
    spaces = L.map((P) => (Array.isArray(P) ? P : pr_hnf(nf, P).filter((c, j) => c[j] === 1n))),
    LV = get_LV(spaces, p, n),
    ramif = nf.discriminant % p === 0n;
  return L.flatMap((P, i): NfPrimeIdeal[] => {
    const f = Array.isArray(P) ? BigInt(n - P.length) : P.f;
    if (limit > 0n && f > limit) return [];
    if (!Array.isArray(P)) return [P];
    if (f === BigInt(n))
      return [
        {
          p,
          generator: Array.from({ length: n }, (_, j) => (j === 0 ? p : 0n)),
          e: 1n,
          f,
          tau: 1n,
        },
      ];
    const generator = uniformizer(nf, P, LV[i]!, p, ramif),
      t = fpKernel(zk_multable(nf.multiplication, generator), p)[0]!,
      pr: NfPrimeIdeal = { p, generator, e: 0n, f, tau: zk_multable(nf.multiplication, t) };
    pr.e = ramif ? 1n + ZC_nfval(t, pr) : 1n;
    return [pr];
  });
}

import { FpM_inv, FpM_image, FpM_suppl } from './alglin1.js';
import { FpX_factor, FpX_roots } from './FpX_factor.js';
import { FpX_divrem, FpX_mul, FpX_red } from './ffinit.js';
import { ZX_mul } from './ZX.js';
import { ZX_sub } from './galconj.js';
import { integerMatrixInverse } from './_matrix_inverse.js';
import { gen_pow_fold } from './bb_group.js';
const basisInverse = new WeakMap<NfIdealData, [bigint[][], bigint]>();
/** Native cached nf.invzk conversion, supplied the exact integral basis. */
function poltobasis(nf: NfIdealData, a: bigint[]): bigint[] {
  let inverse = basisInverse.get(nf);
  if (!inverse) {
    inverse = integerMatrixInverse(nf.basis)!;
    basisInverse.set(nf, inverse);
  }
  return nf.basis.map(
    (_, i) =>
      a.reduce((s, c, j) => s + c * inverse![0][j]![i]! * nf.basisDenominator, 0n) / inverse![1]
  );
}
function basis_norm(nf: NfIdealData, a: bigint[]): bigint {
  const n = nf.multiplication.length,
    power = Array.from({ length: n }, (_, j) =>
      a.reduce((s, c, i) => s + c * nf.basis[i]![j]!, 0n)
    );
  return ZX_resultant(nf.polynomial, power) / nf.basisDenominator ** BigInt(n);
}
const centerPrime = (a: bigint, p: bigint) => {
  const r = a % p;
  return r > p / 2n ? r - p : r < -(p / 2n) ? r + p : r;
};
const primeColumns = (A: bigint[][]) => [[], ...A.map((c) => [0n, ...c])];
const stripPrimeColumns = (A: bigint[][]) => A.slice(1).map((c) => c.slice(1));
const primeImage = (A: bigint[][], p: bigint) => stripPrimeColumns(FpM_image(primeColumns(A), p));

/** Native Dedekind-Kummer prime generator, e/f and anti-uniformizer construction.
 * @see Deviation: Number-field full prime decomposition adapters
 */
export function idealprimedec_kummer(
  nf: NfIdealData,
  u: bigint[],
  e: bigint,
  p: bigint
): NfPrimeIdeal {
  const n = nf.multiplication.length,
    f = BigInt(u.length - 1);
  if (f === BigInt(n))
    return { p, generator: Array.from({ length: n }, (_, i) => (i === 0 ? p : 0n)), e, f, tau: 1n };
  const t = poltobasis(nf, FpX_divrem(nf.polynomial, u, p)[0]).map((x) => centerPrime(x, p)),
    generator = poltobasis(nf, u).map((x) => centerPrime(x, p));
  if (e === 1n && basis_norm(nf, generator) % p ** (f + 1n) === 0n)
    generator[0]! += generator[0]! > 0n ? -p : p;
  return { p, generator, e, f, tau: zk_multable(nf.multiplication, t) };
}

/** Native pradical: Frobenius power kernel and x^p-x on the maximal order. */
function pradical(nf: NfIdealData, p: bigint): [bigint[][], bigint[][]] {
  const n = nf.multiplication.length,
    square = (a: bigint[]) =>
      columnProduct(zk_multable(nf.multiplication, a), a).map((c) => residue(c, p)),
    frob = Array.from({ length: n }, (_, i) => {
      const ei = Array.from({ length: n }, (_, j) => (i === j ? 1n : 0n));
      if (i === 0) return ei;
      const mul = zk_multable(nf.multiplication, ei);
      return gen_pow_fold(ei, p, square, (a) =>
        columnProduct(mul, square(a)).map((c) => residue(c, p))
      );
    });
  let m = frob,
    q = p;
  while (q < BigInt(n)) {
    q *= p;
    m = fpProduct(m, frob, p);
  }
  return [fpKernel(m, p), frob.map((c, j) => c.map((v, i) => v - (i === j ? 1n : 0n)))];
}
/** Native pol_min, taking the first dependency between successive powers. */
function pol_min(mul: bigint[][], p: bigint): bigint[] {
  const n = mul.length,
    powers: bigint[][] = [Array.from({ length: n }, (_, i) => (i === 0 ? 1n : 0n))];
  let z = mul[0]!;
  for (let i = 1; i <= n + 1; i++) {
    powers.push(z);
    if (i < n + 1) z = fpProduct(mul, [z], p)[0]!;
  }
  const polynomial = fpKernel(powers, p)[0]!;
  while (polynomial.at(-1) === 0n) polynomial.pop();
  return polynomial;
}

/** Native cmp_prime_over_p/ideal, using residue degree then generator coordinates. */
export function cmp_prime_ideal(P: NfPrimeIdeal, Q: NfPrimeIdeal): number {
  if (P.p !== Q.p) return P.p < Q.p ? -1 : 1;
  if (P.f !== Q.f) return P.f < Q.f ? -1 : 1;
  for (let i = 0; i < P.generator.length; i++)
    if (P.generator[i] !== Q.generator[i]) return P.generator[i]! < Q.generator[i]! ? -1 : 1;
  return 0;
}

/** Full native idealprimedec, including Kummer and Buchmann-Lenstra splitting.
 * @see Deviation: Number-field full prime decomposition adapters
 */
export function idealprimedec(nf: NfIdealData, p: bigint): NfPrimeIdeal[] {
  return idealprimedec_limit_f(nf, p, 0n);
}
/** Native residue-degree limit; zero requests all primes.
 * @see Deviation: Number-field full prime decomposition adapters
 */
export function idealprimedec_limit_f(nf: NfIdealData, p: bigint, limit: bigint): NfPrimeIdeal[] {
  if (limit < 0n) throw new PariError('domain error in idealprimedec: f < 0');
  const T = nf.polynomial,
    n = nf.multiplication.length,
    F = FpX_factor(T, p).map(([u, e]) => [u, BigInt(e)] as [bigint[], bigint]),
    index = nf.index ?? nf.basisDenominator ** BigInt(n) / absoluteBasisDeterminant(nf.basis);
  if (index % p !== 0n)
    return F.filter(([u]) => !limit || BigInt(u.length - 1) <= limit)
      .map(([u, e]) => idealprimedec_kummer(nf, u, e, p))
      .sort(cmp_prime_ideal);
  const g = F.reduce((a, [b]) => FpX_mul(a, b, p), [1n]),
    h = FpX_divrem(T, g, p)[0],
    f = FpX_red(
      ZX_sub(ZX_mul(g, h), T).map((c) => c / p),
      p
    ),
    L: (bigint[][] | NfPrimeIdeal)[] = [],
    bad: bigint[][] = [];
  for (const [u, e] of F) {
    if (e === 1n || FpX_divrem(f, u, p)[1].length) L.push(idealprimedec_kummer(nf, u, e, p));
    else bad.push(u);
  }
  const [Ip, phi] = pradical(nf, p),
    work: bigint[][][] = [];
  if (L.length) {
    const b = bad.reduce((a, b) => FpX_mul(a, b, p), [1n]),
      coordinates = poltobasis(nf, b).map((c) => residue(c, p));
    work.push(primeImage([...zk_multable(nf.multiplication, coordinates), ...Ip], p));
  } else work.push(Ip);
  const one = Array.from({ length: n }, (_, i) => (i === 0 ? 1n : 0n));
  while (work.length) {
    const H = work.pop()!,
      r = H.length,
      M = stripPrimeColumns(FpM_suppl(primeColumns([...H, one]), p)),
      Mi = stripPrimeColumns(FpM_inv(primeColumns(M), p)!),
      M2 = M.slice(r),
      Mi2 = Mi.map((c) => c.slice(r)),
      phi2 = fpProduct(Mi2, fpProduct(phi, M2, p), p),
      mat1 = fpKernel(phi2, p),
      dim = mat1.length;
    if (dim <= 1) {
      L.push(H);
      continue;
    }
    const a = fpProduct(M2, [mat1[1]!], p)[0]!,
      mula = zk_multable(nf.multiplication, a).map((c) => c.map((v) => residue(v, p))),
      mul2 = fpProduct(Mi2, fpProduct(mula, M2, p), p),
      roots = FpX_roots(pol_min(mul2, p), p);
    for (const root of roots) {
      const I = mula.map((c, j) => c.map((v, i) => v - (i === j ? root : 0n)));
      work.push(primeImage([...H, ...I], p));
    }
    if (roots.length === dim) for (let i = 0; i < roots.length; i++) L.push(work.pop()!);
  }
  return finalize_primes(nf, L, p, limit).sort(cmp_prime_ideal);
}
import { ZM_det as primeBasisDeterminant } from './qfrep.js';
function absoluteBasisDeterminant(A: bigint[][]): bigint {
  const d = primeBasisDeterminant(A);
  return d < 0n ? -d : d;
}
