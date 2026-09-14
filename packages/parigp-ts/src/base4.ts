/** PARI base4.c ideal HNF multiplication and trace-dual inversion.
 * Matrices use zero-indexed columns, matching the port's ideal lattice storage.
 * @see Deviation: Number-field ideal backend adapters
 */
import type { NfIdealData } from './base1.js';
import { ZM_hnfmodid, hnf_divscale } from './hnf_snf.js';
import { gcd } from './ff.js';
import { randomi } from './random.js';
import { residue } from './_polynomial_division.js';
import { PariError } from './errors.js';
import { ZM_lll, LLL_KER } from './lll.js';
import { RgM_mul } from './RgV.js';

type Matrix = bigint[][];
const identity = (n: number): Matrix =>
  Array.from({ length: n }, (_, j) => Array.from({ length: n }, (_, i) => (i === j ? 1n : 0n)));
const transpose = (A: Matrix): Matrix => (A.length ? A[0]!.map((_, i) => A.map((c) => c[i]!)) : []);
const hnf = (A: Matrix, mod: bigint, n: number): Matrix =>
  ZM_hnfmodid(A, Array<bigint>(n).fill(mod));
const equal = (A: Matrix, B: Matrix) =>
  A.length === B.length &&
  A.every((c, j) => c.length === B[j]!.length && c.every((v, i) => v === B[j]![i]));

/** zk_multable: columns are a*w_j in the supplied integral multiplication table. */
export function zk_multable(mul: bigint[][][], a: bigint[]): Matrix {
  const n = mul.length;
  return Array.from({ length: n }, (_, j) =>
    Array.from({ length: n }, (_, k) => {
      let z = 0n;
      for (let i = 0; i < n; i++) z += a[i]! * mul[i]![j]![k]!;
      return z;
    })
  );
}

/** Native get_random_a, with the original first-column scan and rejection loop. */
function get_random_a(mul: bigint[][][], I: Matrix, IZ: bigint): bigint[] {
  const n = mul.length,
    beta: bigint[][] = [],
    tables: Matrix[] = [];
  for (let j = 1; j < n; j++) {
    const T = zk_multable(mul, I[j]!).map((c) => c.map((v) => residue(v, IZ)));
    if (T.every((c) => c.every((v) => v === 0n))) continue;
    if (equal(I, hnf(T, IZ, n))) return I[j]!.slice();
    beta.push(I[j]!);
    tables.push(T);
  }
  for (;;) {
    const z = beta.map(() => randomi(IZ)),
      T = identity(n).map((c) => c.map(() => 0n));
    for (let k = 0; k < z.length; k++)
      for (let j = 0; j < n; j++)
        for (let i = 0; i < n; i++) T[j]![i] += z[k]! * tables[k]![j]![i]!;
    if (equal(I, hnf(T, IZ, n)))
      return Array.from({ length: n }, (_, i) => z.reduce((s, c, j) => s + c * beta[j]![i]!, 0n));
  }
}

/** Native two-element reduction for the matrix multiplication route.
 * The higher-degree smooth-norm approximation specialization remains pending.
 */
export function mat_ideal_two_elt(mul: bigint[][][], input: Matrix): [bigint, bigint[]] {
  const n = mul.length;
  if (n === 2) return [input[0]![0]!, input[1]!.slice()];
  let content = 0n;
  for (const c of input) for (const x of c) content = gcd(content, x);
  const IZ = input[0]![0]!;
  if (IZ === content) return [content, Array<bigint>(n).fill(0n)];
  const I = content === 1n ? input : input.map((c) => c.map((x) => x / content)),
    mod = IZ / content;
  // Reuse PARI's certified randomized route at higher degrees as well; the
  // smooth-norm ideal-approximation specialization is a separate open dependency.
  const a = get_random_a(mul, I, mod).map((x) => {
    const r = residue(x, mod);
    return (r > mod / 2n ? r - mod : r) * content;
  });
  return [IZ, a];
}

function idealHNF_mul_two(mul: bigint[][][], I: Matrix, [a, b]: [bigint, bigint[]]): Matrix {
  const n = mul.length,
    alpha = zk_multable(mul, b),
    columns: Matrix = [];
  for (const c of I)
    columns.push(
      Array.from({ length: n }, (_, i) => c.reduce((s, v, j) => s + alpha[j]![i]! * v, 0n))
    );
  for (const c of I) columns.push(c.map((v) => v * a));
  return hnf(columns, a * I[0]![0]!, n);
}
export function idealHNF_mul(nf: NfIdealData, I: Matrix, J: Matrix): Matrix {
  const i = I[0]![0]!,
    j = J[0]![0]!;
  if (i < j)
    return i === 1n
      ? J.map((c) => c.slice())
      : idealHNF_mul_two(nf.multiplication, J, mat_ideal_two_elt(nf.multiplication, I));
  return j === 1n
    ? I.map((c) => c.slice())
    : idealHNF_mul_two(nf.multiplication, I, mat_ideal_two_elt(nf.multiplication, J));
}

/** Native idealmul matrix route, including primitive content and dimension checks.
 * @see Deviation: Number-field ideal coercion and centered integral bases
 */
export function idealmul(
  nf: NfIdealData,
  I: Matrix,
  iDen: bigint,
  J: Matrix,
  jDen: bigint
): [Matrix, bigint] {
  if (iDen <= 0n || jDen <= 0n) throw new RangeError('idealmul requires positive denominators');
  if (!I.length || !J.length) return [[], 1n];
  const n = nf.multiplication.length;
  if (I.length !== n || J.length !== n || [...I, ...J].some((c) => c.length !== n))
    throw new PariError('inconsistent dimensions in idealmul');
  const primitive = (A: Matrix): [Matrix, bigint] => {
    let content = 0n;
    for (const c of A) for (const v of c) content = gcd(content, v);
    return [A.map((c) => c.map((v) => v / content)), content];
  };
  const [A, a] = primitive(I),
    [B, b] = primitive(J);
  const H = idealHNF_mul(nf, A, B).map((c) => c.map((v) => v * a * b));
  let den = iDen * jDen,
    content = den;
  for (const c of H) for (const v of c) content = gcd(content, v);
  den /= content;
  return [H.map((c) => c.map((v) => v / content)), den];
}

/** Return (I intersect Z)*I^-1 for an integral HNF ideal. */
export function idealHNF_inv_Z(nf: NfIdealData, I: Matrix): Matrix {
  const n = I.length,
    IZ = I[0]![0]!;
  if (IZ === 1n) return identity(n);
  const J = idealHNF_mul_two(nf.multiplication, I, nf.codifferentTwo);
  const dual = transpose(hnf_divscale(J, nf.traceInverse, IZ));
  return hnf(dual, IZ, n);
}
export function idealHNF_inv(nf: NfIdealData, I: Matrix, den: bigint): [Matrix, bigint] {
  const J = idealHNF_inv_Z(nf, I).map((c) => c.map((x) => x * den));
  let d = I[0]![0]!,
    g = d;
  for (const c of J) for (const x of c) g = gcd(g, x);
  d /= g;
  return [J.map((c) => c.map((x) => x / g)), d];
}

/** PARI idealdiv on rational HNF ideals, represented by columns and denominators. */
export function idealdiv(
  nf: NfIdealData,
  I: Matrix,
  iDen: bigint,
  J: Matrix,
  jDen: bigint
): [Matrix, bigint] {
  // idealdiv inverts the divisor before checking whether the dividend is zero.
  if (J.length === 0) throw new PariError('impossible inverse in ginv: 0');
  const [inverse, denominator] = idealHNF_inv(nf, J, jDen);
  if (I.length === 0) return [[], 1n];
  const product = idealHNF_mul(nf, I, inverse);
  let den = iDen * denominator,
    content = den;
  for (const column of product) for (const x of column) content = gcd(content, x);
  den /= content;
  return [product.map((column) => column.map((x) => x / content)), den];
}

/** Native idealintersect on rational HNF columns and positive denominators.
 * The caller supplies canonical ideal HNFs in the field's integral basis.
 * @see Deviation: Number-field ideal intersection and construction adapters
 */
export function idealintersect(
  nf: NfIdealData,
  I: Matrix,
  iDen: bigint,
  J: Matrix,
  jDen: bigint
): [Matrix, bigint] {
  const n = nf.multiplication.length;
  if (
    iDen <= 0n ||
    jDen <= 0n ||
    [I, J].some((A) => A.length && (A.length !== n || A.some((c) => c.length !== n)))
  )
    throw new RangeError('idealintersect requires square ideal HNFs and positive denominators');
  if (!I.length || !J.length) return [[], 1n];
  // Q_remove_denom returns the least denominator of each rational matrix.
  const remove = (A: Matrix, d: bigint): [Matrix, bigint] => {
    let g = d;
    for (const c of A) for (const v of c) g = gcd(g, v);
    return [A.map((c) => c.map((v) => v / g)), d / g];
  };
  [I, iDen] = remove(I, iDen);
  [J, jDen] = remove(J, jDen);
  const X = I.map((c) => c.map((v) => v * jDen)),
    Y = J.map((c) => c.map((v) => v * iDen));
  const K = (ZM_lll([...X, ...Y], 0.99, LLL_KER) as Matrix).map((c) => c.slice(0, n));
  const xZ = X[0]![0]!,
    yZ = Y[0]![0]!;
  const H = hnf(RgM_mul(X, K) as Matrix, (xZ / gcd(xZ, yZ)) * yZ, n);
  return remove(H, iDen * jDen);
}

import type { NfPrimeIdeal } from './base2.js';
import { Z_pvalrem } from './gen2.js';

/** Native idealHNF_val, with norm/intersection bounds and modular tau products. */
function idealHNF_val(A: Matrix, P: NfPrimeIdeal, Nval: bigint, Zval: bigint): bigint {
  if (Nval < P.f) return 0n;
  const e = P.e,
    p = P.p,
    vmax = Zval * e < Nval / P.f ? Zval * e : Nval / P.f,
    mul = P.tau as Matrix,
    n = mul.length,
    B: Matrix = [[]],
    vals: bigint[] = [0n];
  for (let j = 1; j < n; j++) {
    const y = Array.from({ length: n }, (_, i) =>
      A[j]!.slice(0, j + 1).reduce((s, c, k) => s + c * mul[k]![i]!, 0n)
    );
    if (y.some((c) => c % p !== 0n)) return 0n;
    B[j] = y.map((c) => c / p);
  }
  for (let j = 1; j < n; j++) {
    const content = B[j]!.reduce((g, c) => gcd(g, c), 0n);
    B[j] = B[j]!.map((c) => c / content);
    vals[j] = 1n + e * BigInt(Z_pvalrem(content, p)[0]);
  }
  let pk = p ** ((vmax + e - 1n) / e);
  for (let v = 1n; v < vmax; v++) {
    if (e === 1n || (vmax - v) % e === 0n) pk /= p;
    const cutoff = 1n << BigInt(64 * Math.ceil(pk.toString(2).length / 64));
    for (let j = 1; j < n; j++) {
      if (v < vals[j]!) continue;
      const y = Array.from({ length: n }, (_, i) =>
        B[j]!.reduce((s, c, k) => s + c * mul[k]![i]!, 0n)
      );
      if (y.some((c) => c % p !== 0n)) return v;
      B[j] = y.map((c) => {
        const a = c / p;
        return a >= cutoff || a <= -cutoff ? a % pk : a;
      });
    }
  }
  return vmax;
}

/** Native idealval/gpidealval matrix route, with exact fractional content.
 * Zero HNF returns the port's standard infinity sentinel.
 * @see Deviation: Number-field ideal valuation adapters
 */
export function idealval(
  nf: NfIdealData,
  I: Matrix,
  den: bigint,
  P: NfPrimeIdeal
): bigint | 'Infinity' {
  if (den <= 0n) throw new RangeError('idealval requires a positive denominator');
  if (!I.length) return 'Infinity';
  const n = nf.multiplication.length;
  if (I.length !== n || I.some((c) => c.length !== n))
    throw new RangeError('idealval requires a square ideal HNF');
  const content = I.reduce((g, c) => c.reduce((h, x) => gcd(h, x), g), 0n),
    A = I.map((c) => c.map((x) => x / content)),
    vc = BigInt(Z_pvalrem(content, P.p)[0] - Z_pvalrem(den, P.p)[0]);
  if (P.tau === 1n) return vc;
  const Zval = BigInt(Z_pvalrem(A[0]![0]!, P.p)[0]);
  if (!Zval) return vc * P.e;
  const Nval = A.reduce((v, c, j) => v + BigInt(Z_pvalrem(c[j]!, P.p)[0]), 0n);
  return idealHNF_val(A, P, Nval, Zval) + vc * P.e;
}

import { ZM_hnfmodprime } from './hnf_snf.js';
/** Native prime-ideal HNF, with the inert scalar shortcut.
 * @see Deviation: Number-field full prime decomposition adapters
 */
export function pr_hnf(nf: NfIdealData, P: NfPrimeIdeal): Matrix {
  if (P.tau === 1n) return identity(nf.multiplication.length).map((c) => c.map((x) => x * P.p));
  return ZM_hnfmodprime(zk_multable(nf.multiplication, P.generator), P.p);
}

import { idealprimedec, idealprimedec_limit_f, cmp_prime_ideal } from './base2.js';
import { Z_factor, isPrime } from './ifactor.js';
import { FpX_is_irred } from './ffinit.js';

/** Native idealnumden matrix route, returning coprime integral HNFs.
 * A null denominator denotes the unit ideal; an empty numerator denotes zero.
 * @see Deviation: Number-field full prime decomposition adapters
 */
export function idealnumden(nf: NfIdealData, I: Matrix, den: bigint): [Matrix, Matrix | null] {
  if (den <= 0n) throw new RangeError('idealnumden requires a positive denominator');
  if (!I.length) return [[], null];
  const n = nf.multiplication.length;
  if (I.length !== n || I.some((c) => c.length !== n))
    throw new PariError('inconsistent dimensions in idealnumden');
  const content = I.reduce((g, c) => c.reduce((h, x) => gcd(h, x), g), den),
    x = I.map((c) => c.map((v) => v / content)),
    d = den / content;
  if (d === 1n) return [x, null];
  const J = hnf(x, d, n),
    c = J[0]![0]!,
    B = idealHNF_inv_Z(nf, J).map((col) => col.map((v) => v * (d / c))),
    A = idealHNF_mul(nf, B, x).map((col) => col.map((v) => v / d));
  return [A, B];
}
/** Native integral HNF factorization: primitive part, content, norm and residue limits. */
function idealHNF_factor(nf: NfIdealData, I: Matrix): [NfPrimeIdeal, bigint][] {
  let content = I.reduce((g, c) => c.reduce((h, x) => gcd(h, x), g), 0n);
  const A = I.map((c) => c.map((v) => v / content)),
    out: [NfPrimeIdeal, bigint][] = [];
  for (const [p, z] of Z_factor(A[0]![0]!)) {
    let Nval = A.reduce((v, c, j) => v + BigInt(Z_pvalrem(c[j]!, p)[0]), 0n);
    const [vc, next] = Z_pvalrem(content, p);
    content = next;
    const primes = vc ? idealprimedec(nf, p) : idealprimedec_limit_f(nf, p, Nval);
    let j = 0;
    for (; Nval && j < primes.length; j++) {
      const P = primes[j]!,
        v = idealHNF_val(A, P, Nval, z);
      Nval -= v * P.f;
      const e = v + BigInt(vc) * P.e;
      if (e) out.push([P, e]);
    }
    if (vc) for (; j < primes.length; j++) out.push([primes[j]!, BigInt(vc) * primes[j]!.e]);
  }
  if (content !== 1n)
    for (const [p, e] of Z_factor(content))
      for (const P of idealprimedec(nf, p)) out.push([P, e * P.e]);
  return out;
}
/** Native idealfactor on rational HNF ideals; no norm-based unit shortcut.
 * @see Deviation: Number-field full prime decomposition adapters
 */
export function idealfactor(nf: NfIdealData, I: Matrix, den: bigint): [NfPrimeIdeal, bigint][] {
  if (!I.length) throw new PariError('domain error in idealfactor: ideal = 0');
  const [numerator, denominator] = idealnumden(nf, I, den),
    out = idealHNF_factor(nf, numerator);
  if (denominator)
    out.push(
      ...idealHNF_factor(nf, denominator).map(([P, e]) => [P, -e] as [NfPrimeIdeal, bigint])
    );
  return out.sort(([P], [Q]) => cmp_prime_ideal(P, Q));
}
/** Native idealismaximal_int, with the irreducibility shortcut when p avoids the index. */
function idealismaximal_int(nf: NfIdealData, p: bigint): NfPrimeIdeal | null {
  if (!isPrime(p)) return null;
  if (nf.index !== undefined && nf.index % p !== 0n && !FpX_is_irred(nf.polynomial, p)) return null;
  const L = idealprimedec(nf, p);
  return L.length === 1 && L[0]!.e === 1n ? L[0]! : null;
}
/** Native idealismaximal matrix route, avoiding factorization of the ideal norm.
 * @see Deviation: Number-field full prime decomposition adapters
 */
export function idealismaximal(nf: NfIdealData, I: Matrix, den: bigint): NfPrimeIdeal | null {
  if (!I.length) return null;
  const content = I.reduce((g, c) => c.reduce((h, x) => gcd(h, x), g), 0n),
    A = I.map((c) => c.map((v) => v / content)),
    p = A[0]![0]!;
  if (content !== den)
    return content % den === 0n && p === 1n ? idealismaximal_int(nf, content / den) : null;
  if (!isPrime(p)) return null;
  let f = 1n;
  for (let i = 1; i < A.length; i++) {
    const c = A[i]![i]!;
    if (c === p) f++;
    else if (c !== 1n) return null;
  }
  const L = idealprimedec_limit_f(nf, p, f);
  for (let i = L.length - 1; i >= 0; i--) {
    const P = L[i]!;
    if (P.f !== f) break;
    if (idealval(nf, A, 1n, P) === 1n) return P;
  }
  return null;
}

import { zkmultable_capZ } from './base3.js';
/** Native idealhnf_principal, using primitive content and a modular HNF bound.
 * Input is integral-basis coordinates over a positive common denominator.
 * @see Deviation: Number-field full prime decomposition adapters
 */
export function idealhnf_principal(
  nf: NfIdealData,
  coordinates: bigint[],
  den: bigint
): [Matrix, bigint] {
  if (den <= 0n) throw new RangeError('idealhnf_principal requires a positive denominator');
  const n = nf.multiplication.length;
  if (coordinates.length !== n)
    throw new RangeError('idealhnf_principal requires integral-basis coordinates');
  const content = coordinates.reduce((g, c) => gcd(g, c), 0n);
  if (!content) return [[], 1n];
  const divisor = gcd(content, den),
    scale = content / divisor,
    d = den / divisor;
  if (coordinates.slice(1).every((c) => c === 0n))
    return [identity(n).map((c) => c.map((x) => x * scale)), d];
  const M = zk_multable(
      nf.multiplication,
      coordinates.map((c) => c / content)
    ),
    H = hnf(M, zkmultable_capZ(M), n);
  return [H.map((c) => c.map((x) => x * scale)), d];
}
/** Native idealadd on rational HNFs, using gcd of their intersections with Z.
 * @see Deviation: Number-field full prime decomposition adapters
 */
export function idealadd(
  nf: NfIdealData,
  I: Matrix,
  iDen: bigint,
  J: Matrix,
  jDen: bigint
): [Matrix, bigint] {
  if (!I.length) return [J.map((c) => c.slice()), jDen];
  if (!J.length) return [I.map((c) => c.slice()), iDen];
  const d = (iDen / gcd(iDen, jDen)) * jDen,
    A = I.map((c) => c.map((x) => x * (d / iDen))),
    B = J.map((c) => c.map((x) => x * (d / jDen))),
    mod = gcd(A[0]![0]!, B[0]![0]!);
  const H =
      mod === 1n
        ? identity(nf.multiplication.length)
        : hnf([...A, ...B], mod, nf.multiplication.length),
    content = H.reduce((g, c) => c.reduce((h, x) => gcd(h, x), g), d);
  return [H.map((c) => c.map((x) => x / content)), d / content];
}
