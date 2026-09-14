/** NTL ZZXFactoring.cpp:101 SquareFreeDecomp (primitive input with positive leading coefficient).
 * @see Deviation: NTL integer polynomial squarefree decomposition
 */
import { GCD, _ZZX_kernels as k } from './ZZX1.js';
import { mul as modularProduct, type PolynomialProductState } from './ZZ_pX.js';
export function SquareFreeDecomp(ff: readonly bigint[], state?: PolynomialProductState): Array<[bigint[], number]> {
  const f = k.normalized(ff), out: Array<[bigint[], number]> = [];
  if (f.length <= 1) return out;
  const productState = state ?? {
    context: new FFTPrimeContext(), stream: new RandomStream(new Uint8Array(32)),
  };
  const diff = (a: bigint[]) => k.normalized(a.slice(1).map((x, i) => x * BigInt(i + 1)));
  const derivative = diff(f);
  let d = GCD(f, derivative, productState);
  if (d.length === 1) return [[f, 1]];
  let v = k.divide(f, d, productState)!, w = k.divide(derivative, d, productState)!, i = 0;
  for (;;) {
    i++;
    const t = diff(v), s = k.normalized(Array.from({length: Math.max(w.length, t.length)}, (_, j) => (w[j] ?? 0n) - (t[j] ?? 0n)));
    if (!s.length) { if (v.length > 1) out.push([v, i]); return out; }
    d = GCD(v, s, productState);
    v = k.divide(v, d, productState)!; w = k.divide(s, d, productState)!;
    if (d.length > 1) out.push([d, i]);
  }
}

const mod = (a: bigint, p: bigint) => ((a % p) + p) % p;
const bits = (a: bigint) => a === 0n ? 0 : (a < 0n ? -a : a).toString(2).length;
/** NTL ZZXFactoring.cpp ComputeTrace; returns a copied, updated trace vector.
 * @see Deviation: NTL integer factorization trace adapters
 */
export function ComputeTrace(Tr: readonly bigint[], f: readonly bigint[], d: number, P: bigint): bigint[] {
  let n = f.length - 1;
  while (n >= 0 && f[n] === 0n) n--;
  if (n <= 0 || f[n] !== 1n) throw new Error('ComputeTrace: internal error (1)');
  if (d <= 0) throw new Error('ComputeTrace: internal error (2)');
  if (Tr.length < d) throw new Error('ComputeTrace: internal error (3)');
  if (P <= 1n) throw new Error('ComputeTrace: internal error (4)');
  let t = 0n;
  if (d > n) {
    for (let i = 1; i <= n; i++) t += Tr[i + d - n - 2]! * f[i - 1]!;
  } else {
    t = f[n - d]! * BigInt(d);
    for (let i = 1; i < d; i++) t += Tr[i - 1]! * f[n - d + i]!;
  }
  const out = Tr.slice();
  out[d - 1] = mod(-t, P);
  return out;
}
/** NTL ZZXFactoring.cpp ChopTraces; output dimensions and unused suffix are retained.
 * @see Deviation: NTL integer factorization trace adapters
 */
export function ChopTraces(C: readonly bigint[], Tr: readonly bigint[], d: number, pb: readonly bigint[], pdelta: bigint, P: bigint, lc: bigint): bigint[] {
  if (d <= 0) throw new Error('ChopTraces: internal error (1)');
  if (C.length < d) throw new Error('ChopTraces: internal error (2)');
  if (Tr.length < d) throw new Error('ChopTraces: internal error (3)');
  if (pb.length < d) throw new Error('ChopTraces: internal error (4)');
  if (P <= 1n) throw new Error('ChopTraces: internal error (5)');
  const out = C.slice(), half = pdelta >> 1n, lcred = mod(lc, P);
  let lcpow = 1n;
  for (let i = 0; i < d; i++) {
    lcpow = lcpow * lcred % P;
    let t = mod((mod(lcpow * Tr[i]!, P) + (pb[i]! >> 1n)) / pb[i]!, pdelta);
    if (t > half) t -= pdelta;
    out[i] = t;
  }
  return out;
}
/** NTL dense trace combination; A uses rows with no dummy entry.
 * @see Deviation: NTL integer factorization trace adapters
 */
export function DenseChopTraces(C: readonly bigint[], Tr: readonly bigint[], d: number, d1: number, pb_eff: bigint, pdelta: bigint, P: bigint, lc: bigint, A: readonly (readonly bigint[])[]): bigint[] {
  const out = C.slice(), half = pdelta >> 1n, halfb = pb_eff >> 1n, lcred = mod(lc, P);
  for (let i = 0; i < d1; i++) {
    let lcpow = 1n, acc = 0n;
    for (let j = 0; j < d; j++) {
      lcpow = lcpow * lcred % P;
      const t = mod(lcpow * Tr[j]!, P) * mod(A[i]![j]!, P) % P;
      acc = (acc + t) % P;
    }
    let t = mod((acc + halfb) / pb_eff, pdelta);
    if (t > half) t -= pdelta;
    out[i] = t;
  }
  return out;
}
/** NTL trace precision update, returning [b, pb].
 * @see Deviation: NTL integer factorization trace adapters
 */
export function Compute_pb(b: readonly number[], pb: readonly bigint[], p: number | bigint, d: number, root_bound: bigint, n: number): [number[], bigint[]] {
  const bound = 2n * root_bound ** BigInt(d) * BigInt(n);
  let i = d === 1 ? 0 : b[d - 2]!, t = d === 1 ? 1n : pb[d - 2]!;
  while (t <= bound) { i++; t *= BigInt(p); }
  const bb = b.slice(0, d), pp = pb.slice(0, d);
  bb[d - 1] = i; pp[d - 1] = t;
  return [bb, pp];
}
/** NTL extra precision update, returning [delta, pdelta].
 * @see Deviation: NTL integer factorization trace adapters
 */
export function Compute_pdelta(delta: number, pdelta: bigint, p: number | bigint, bit_delta: number): [number, bigint] {
  while (bits(pdelta) <= bit_delta) { delta++; pdelta *= BigInt(p); }
  return [delta, pdelta];
}
/** NTL lattice construction, returning [M, C]; all matrices use rows.
 * @see Deviation: NTL integer factorization trace adapters
 */
export function BuildReductionMatrix(r: number, d: number, pdelta: bigint, chop_vec: readonly (readonly bigint[])[], B_L: readonly (readonly bigint[])[]): [bigint[][], number] {
  const s = B_L.length, C = Math.trunc(Math.sqrt(d * r) / 2) + 1;
  const M = Array.from({length:s+d}, () => Array<bigint>(r+d).fill(0n));
  for (let i = 0; i < s; i++) {
    for (let j = 0; j < r; j++) M[i]![j] = B_L[i]![j]! * BigInt(C);
    for (let j = 0; j < d; j++) {
      let t = 0n;
      for (let k = 0; k < r; k++) t += B_L[i]![k]! * chop_vec[k]![j]!;
      t = mod(t, pdelta);
      if (t > (pdelta >> 1n)) t -= pdelta;
      M[i]![j+r] = t;
    }
  }
  for (let i = 0; i < d; i++) M[i+s]![i+r] = pdelta;
  return [M, C];
}
/** NTL compressed-trace bound, returning [b_eff, pb_eff].
 * @see Deviation: NTL integer factorization trace adapters
 */
export function Compute_pb_eff(p: number | bigint, d: number, root_bound: bigint, n: number, ran_bits: number): [number, bigint] {
  const bound = root_bound === 1n ? (BigInt(d) * BigInt(n)) << BigInt(ran_bits + 1) : (root_bound ** BigInt(d) * BigInt(n)) << BigInt(ran_bits + 2);
  let i = 0, t = 1n;
  while (t <= bound) { i++; t *= BigInt(p); }
  return [i, t];
}
/** NTL number of compressed traces.
 * @see Deviation: NTL integer factorization trace adapters
 */
export function d1_val(bit_delta: number, r: number, s: number): number {
  return Math.trunc(0.30 * r * s / bit_delta) + 1;
}

import { _ZZ_pX_euclidean_kernels } from './ZZ_pX1.js';
import { XGCD as wordXGCD } from './lzz_pX1.js';
import { mul as wordMul, type zz_pXOptions } from './lzz_pX.js';
import { ZZ_pXModulus, ZZ_pXMultiplier, rem as quotientRem, MulMod as quotientMul } from './ZZ_pX.js';
import { mul as ZZX_mul } from './ZZX1.js';
const add = (a: bigint[], b: bigint[]) =>
  k.normalized(
    Array.from({ length: Math.max(a.length, b.length) }, (_, i) => (a[i] ?? 0n) + (b[i] ?? 0n))
  );
const sub = (a: bigint[], b: bigint[]) =>
  k.normalized(
    Array.from({ length: Math.max(a.length, b.length) }, (_, i) => (a[i] ?? 0n) - (b[i] ?? 0n))
  );
const exactScalarQuotient = (a: bigint[], p: bigint) => {
  if (a.some((c) => c % p !== 0n)) throw new Error('inexact division');
  return k.normalized(a.map((c) => c / p));
};
/** NTL MultiLift: balanced factor tree and staged quadratic Hensel lifting.
 * Coefficients are constant-first; p replaces the active native word-prime context.
 * Optional word context settings and shared product state survive all lifting stages.
 * @see Deviation: NTL multifactor Hensel lifting and integer products
 * @see Deviation: NTL stateful Hensel lifting
 */
export function MultiLift(
  a: readonly (readonly bigint[])[],
  f: readonly bigint[],
  e: number,
  p: bigint,
  options?: zz_pXOptions
): bigint[][] {
  const arithmetic = new zz_pXModulus(null, p, options).arithmetic,
    factors = a.map((x) => arithmetic.norm(x)),
    F = k.normalized(f),
    count = factors.length;
  if (count < 2 || e < 1 || e >= 2 ** 60 || F.at(-1) !== 1n || factors.some((x) => x.at(-1) !== 1n))
    throw new Error('MultiLift: bad args');
  if (e === 1) return factors;
  const E = [e];
  while (e > 1) {
    e = Math.floor((e + 1) / 2);
    E.push(e);
  }
  const v = Array<bigint[]>(2 * count - 2),
    w = Array<bigint[]>(2 * count - 2),
    link = Array<number>(2 * count - 2);
  for (let i = 0; i < count; i++) {
    v[i] = factors[i]!;
    link[i] = -(i + 1);
  }
  for (let j = 0, i = count; j < 2 * count - 4; j += 2) {
    let min = j;
    for (let s = j + 1; s < i; s++) if (v[s]!.length < v[min]!.length) min = s;
    [v[j], v[min]] = [v[min]!, v[j]!];
    [link[j], link[min]] = [link[min]!, link[j]!];
    min = j + 1;
    for (let s = j + 2; s < i; s++) if (v[s]!.length < v[min]!.length) min = s;
    [v[j + 1], v[min]] = [v[min]!, v[j + 1]!];
    [link[j + 1], link[min]] = [link[min]!, link[j + 1]!];
    v[i] = wordMul(v[j]!, v[j + 1]!, p, options);
    link[i] = j;
    i++;
  }
  for (let j = 0; j < 2 * count - 2; j += 2) {
    const [g, A, B] = wordXGCD(v[j]!, v[j + 1]!, p, options);
    if (g.length !== 1 || g[0] !== 1n) throw new Error('relatively prime polynomials expected');
    w[j] = A;
    w[j + 1] = B;
  }
  for (let i = E.length - 1; i > 0; i--) {
    const p0 = p ** BigInt(E[i]!),
      p1 = p ** BigInt(E[i - 1]! - E[i]!),
      ar = new ZZ_pXModulus(null, p1, options?.state).arithmetic;
    const lift = (target: bigint[], j: number): void => {
      if (j < 0) return;
      const g = v[j]!,
        h = v[j + 1]!,
        A = w[j]!,
        B = w[j + 1]!;
      const cc = ar.norm(exactScalarQuotient(sub(target, ZZX_mul(g, h, options?.state)), p0)),
        gg = ar.norm(g),
        hh = ar.norm(h),
        aa = ar.norm(A),
        bb = ar.norm(B);
      // HenselLift builds both quotient objects before either multiplier.
      // The final HenselLift1 stage uses raw products and skips inverse caches.
      const GG = new ZZ_pXModulus(gg, p1, options?.state);
      const HH = new ZZ_pXModulus(hh, p1, options?.state);
      const AA = i !== 1 ? new ZZ_pXMultiplier(aa, HH) : aa;
      const BB = i !== 1 ? new ZZ_pXMultiplier(bb, GG) : bb;
      const g1 = quotientMul(quotientRem(cc, GG), BB, GG),
        h1 = quotientMul(quotientRem(cc, HH), AA, HH);
      const G = add(
          g,
          g1.map((c) => c * p0)
        ),
        H = add(
          h,
          h1.map((c) => c * p0)
        );
      if (i !== 1) {
        const rr = ar.norm(
          exactScalarQuotient(
            sub([1n], add(ZZX_mul(A, G, options?.state), ZZX_mul(B, H, options?.state))),
            p0
          )
        );
        const a1 = quotientMul(quotientRem(rr, HH), AA, HH),
          b1 = quotientMul(quotientRem(rr, GG), BB, GG);
        w[j] = add(
          A,
          a1.map((c) => c * p0)
        );
        w[j + 1] = add(
          B,
          b1.map((c) => c * p0)
        );
      }
      v[j] = G;
      v[j + 1] = H;
      lift(G, link[j]!);
      lift(H, link[j + 1]!);
    };
    lift(F, v.length - 2);
  }
  const out = Array<bigint[]>(count);
  for (let i = 0; i < v.length; i++) if (link[i]! < 0) out[-link[i]! - 1] = v[i]!;
  return out;
}

import { image as ntl_image, LLL_plus as ntl_LLL_plus } from './LLL.js';
/** Native NTL ZZXFactoring.cpp PolyEval.
 * @see Deviation: NTL factor recovery bounds and lifting adapters
 */
export function PolyEval(f: readonly bigint[], a: bigint): bigint {
  const F = k.normalized(f);
  if (!F.length) throw new Error('PolyEval: internal error');
  let acc = F.at(-1)!;
  for (let i = F.length - 2; i >= 0; i--) acc = acc * a + F[i]!;
  return acc;
}
/** Native NTL ZZXFactoring.cpp RootBound.
 * @see Deviation: NTL factor recovery bounds and lifting adapters
 */
export function RootBound(f: readonly bigint[]): bigint {
  const F = k.normalized(f);
  if (!F.length || F[0] === 0n) throw new Error('RootBound: internal error');
  const n = F.length - 1,
    g = F.map((c, i) => (i === n ? (c < 0n ? -c : c) : c > 0n ? -c : c));
  let lb = 0n,
    ub = 1n;
  while (PolyEval(g, ub) < 0n) ub *= 2n;
  while (ub - lb > 1n) {
    const mb = (ub + lb) / 2n;
    if (PolyEval(g, mb) < 0n) lb = mb;
    else ub = mb;
  }
  return ub * g[n]!;
}
/** Native NTL ZZXFactoring.cpp CutAway.
 * @see Deviation: NTL factor recovery bounds and lifting adapters
 */
export function CutAway(
  D: readonly bigint[],
  M: readonly (readonly bigint[])[],
  C: number,
  r: number,
  d: number
): bigint[][] {
  let count = M.length;
  const c = BigInt(C),
    rr = BigInt(r),
    bound = 4n * c * c * rr + BigInt(d) * rr * rr;
  while (count >= 1 && 4n * D[count]! > bound * D[count - 1]!) count--;
  const B = Array.from({ length: count }, (_, i) =>
    Array.from({ length: r }, (_, j) => {
      if (c === 0n) throw new Error('division by zero in _ntl_gsdiv');
      const v = M[i]![j]!,
        q = v / c,
        rem = v % c;
      return rem !== 0n && rem < 0n !== c < 0n ? q - 1n : q;
    })
  );
  const [rank, , basis] = ntl_image(B);
  return basis.slice(count - rank);
}
/** Native NTL ZZXFactoring.cpp AdditionalLifting.
 * @see Deviation: NTL stateful Hensel lifting
 * @see Deviation: NTL factor recovery bounds and lifting adapters
 */
export function AdditionalLifting(
  P1: bigint,
  e1: number,
  w1: readonly (readonly bigint[])[],
  p: bigint,
  new_bound: number,
  f: readonly bigint[],
  doubling: boolean,
  state?: PolynomialProductState
): [bigint, number, bigint[][]] {
  const next = doubling ? Math.max(2 * e1, new_bound) : new_bound;
  if (next < 0) throw new Error('negative exponent in _ntl_zexps');
  if (
    next > 0 &&
    p !== 0n &&
    BigInt((p < 0n ? -p : p).toString(2).length) * BigInt(next) > (1n << 63n) - 64n
  )
    throw new Error('overflow in _ntl_gexps');
  const modulus = p ** BigInt(next),
    F = k.normalized(f),
    lc = F.at(-1) ?? 0n;
  let target: bigint[];
  if (lc === 1n) target = F;
  else if (lc === -1n) target = F.map((c) => -c);
  else {
    if (modulus === 0n) throw new Error('division by zero in _ntl_gdiv');
    if (modulus <= 1n) throw new Error('InvMod: second input <= 1');
    const ar = _ZZ_pX_euclidean_kernels(modulus);
    let inv: bigint;
    try {
      inv = ar.inverse(ar.mod(lc));
    } catch {
      throw new Error('InvMod: inverse undefined');
    }
    target = F.map((c) => ar.mod(c * inv));
  }
  const degree = F.length - 1;
  const maxroot = (degree <= 1 ? 0 : BigInt(degree - 1).toString(2).length) + 1;
  const factors = MultiLift(w1, target, next, p, { maxroot, state });
  return [modulus, next, factors];
}

/** NTL ZZXFactoring.cpp inplace_rev; returns a normalized independent polynomial.
 * @see Deviation: NTL factor selection adapters
 */
export function inplace_rev(f: readonly bigint[]): bigint[] {
  return k.normalized(k.normalized(f).reverse());
}
/** NTL distinct-degree pattern recording; repeated degree entries overwrite.
 * @see Deviation: NTL factor selection adapters
 */
export function RecordPattern(
  pat: readonly number[],
  fac: readonly (readonly [readonly bigint[], number])[],
  p: bigint
): number[] {
  const out = Array<number>(pat.length).fill(0),
    kp = _ZZ_pX_euclidean_kernels(p);
  for (const [f, d] of fac) out[d] = Math.trunc((kp.norm(f).length - 1) / d);
  return out;
}
/** NTL sum of the factor counts, including index zero.
 * @see Deviation: NTL factor selection adapters
 */
export function NumFactors(pat: readonly number[]): number {
  let n = 0;
  for (const x of pat) n += x;
  return n;
}
/** NTL subset-degree bitsets; the factor overload retains native suffix slots.
 * @see Deviation: NTL factor selection adapters
 */
export function CalcPossibleDegrees(pat: readonly number[]): bigint;
export function CalcPossibleDegrees(
  fac: readonly (readonly bigint[])[],
  count: number,
  p: bigint
): bigint[];
export function CalcPossibleDegrees(
  input: readonly number[] | readonly (readonly bigint[])[],
  count?: number,
  p?: bigint
): bigint | bigint[] {
  if (count === undefined) {
    let pd = 1n;
    const pat = input as readonly number[];
    for (let d = 1; d < pat.length; d++) for (let j = 0; j < pat[d]!; j++) pd |= pd << BigInt(d);
    return pd;
  }
  const fac = input as readonly (readonly bigint[])[],
    r = fac.length,
    S = Array<bigint>(r).fill(0n);
  if (!r) return S;
  if (count < 1 || count > r) throw new Error('CalcPossibleDegrees: bad args');
  const kp = _ZZ_pX_euclidean_kernels(p!),
    degree = fac.map((f) => BigInt(kp.norm(f).length - 1));
  S[r - 1] = 1n << degree[r - 1]!;
  for (let i = r - 2; i >= 0; i--) S[i] = (1n << degree[i]!) | S[i + 1]!;
  for (let l = 2; l <= count; l++) {
    let old = S[r - l]!;
    S[r - l] = S[r - l + 1]! << degree[r - l]!;
    for (let i = r - l - 1; i >= 0; i--) {
      const t = old << degree[i]!;
      old = S[i]!;
      S[i] = S[i + 1]! | t;
    }
  }
  return S;
}
/** NTL constant-term divisibility test with copied prefix-product cache.
 * @see Deviation: NTL factor selection adapters
 */
export function ConstTermTest(
  W: readonly (readonly bigint[])[],
  I: readonly number[],
  ct: bigint,
  lc: bigint,
  prod: readonly bigint[],
  ProdLen: number,
  p: bigint
): [number, bigint[], number] {
  const out = prod.map((x) => mod(x, p)),
    count = I.length;
  if (ProdLen === 0) {
    out[0] = mod(lc * (W[I[0]!]![0] ?? 0n), p);
    ProdLen++;
  }
  for (let i = ProdLen; i < count; i++) out[i] = mod(out[i - 1]! * (W[I[i]!]![0] ?? 0n), p);
  let t = out[count - 1]!;
  if (t > p >> 1n) t -= p;
  return [Number(t === 0n ? ct === 0n : ct % t === 0n), out, count - 1];
}
/** NTL balanced coefficient copy, with the positive midpoint retained.
 * @see Deviation: NTL factor selection adapters
 */
export function BalCopy(G: readonly bigint[], p: bigint): bigint[] {
  return _ZZ_pX_euclidean_kernels(p)
    .norm(G)
    .map((x) => (x > p >> 1n ? x - p : x));
}
/** NTL minimum-degree product tree, optionally selecting indices in caller order.
 * @see Deviation: NTL stateful modular polynomial products
 * @see Deviation: NTL factor selection adapters
 */
export function mul(
  W: readonly (readonly bigint[])[],
  p: bigint,
  I?: readonly number[],
  state?: PolynomialProductState
): bigint[] {
  const kp = _ZZ_pX_euclidean_kernels(p),
    a = (I === undefined ? W : I.map((i) => W[i]!)).map((f) => kp.norm(f));
  if (!a.length) return [1n];
  for (let i = 1; i < a.length; i++)
    for (let j = 0; j < a.length - i; j++)
      if (a[j]!.length < a[j + 1]!.length) [a[j], a[j + 1]] = [a[j + 1]!, a[j]!];
  while (a.length > 1) {
    const last = a.pop()!;
    a[a.length - 1] = state
      ? modularProduct(a[a.length - 1]!, last, p, state)
      : kp.mul(a[a.length - 1]!, last);
    let i = a.length - 1;
    while (i > 0 && a[i - 1]!.length < a[i]!.length) {
      [a[i - 1], a[i]] = [a[i]!, a[i - 1]!];
      i--;
    }
  }
  return a[0]!;
}

/** NTL product of the factors complementary to sorted, distinct indices I.
 * @see Deviation: NTL factor selection adapters
 * @see Deviation: NTL factor recombination cache state
 */
export function InvMul(
  W: readonly (readonly bigint[])[],
  I: readonly number[],
  p: bigint,
  state?: PolynomialProductState
): bigint[] {
  const remaining: bigint[][] = [];
  for (let j = 0, i = 0; j < W.length; j++) {
    if (i < I.length && j === I[i]) i++;
    else remaining.push(W[j]!.slice());
  }
  return mul(remaining, p, undefined, state);
}
/** NTL stable removal of sorted, distinct factor indices, returning copies.
 * @see Deviation: NTL factor selection adapters
 */
export function RemoveFactors(
  W: readonly (readonly bigint[])[],
  I: readonly number[],
  p: bigint
): bigint[][] {
  const kp = _ZZ_pX_euclidean_kernels(p),
    out: bigint[][] = [];
  for (let j = 0, i = 0; j < W.length; j++) {
    if (i < I.length && j === I[i]) i++;
    else out.push(kp.norm(W[j]!));
  }
  return out;
}
/** NTL bit unpacking uses the absolute value, including for negative integers.
 * @see Deviation: NTL factor selection adapters
 */
export function unpack(a: bigint, n: number): number[] {
  if (n < -1) throw new Error('negative length in vector::SetLength');
  const v = a < 0n ? -a : a;
  return Array.from({ length: n + 1 }, (_, i) => Number((v >> BigInt(i)) & 1n));
}
/** NTL subtraction of factor-count patterns, returning an independent vector.
 * @see Deviation: NTL factor selection adapters
 */
export function SubPattern(p1: readonly number[], p2: readonly number[]): number[] {
  if (p1.length !== p2.length) throw new Error('SubPattern: bad args');
  const out = p1.slice();
  for (let i = 0; i < out.length; i++) {
    out[i] = out[i]! - p2[i]!;
    if (out[i]! < 0) throw new Error('SubPattern: internal error');
  }
  return out;
}

import { PrimeSeq, RandomStream, RandomBits, RandomBnd } from './ZZ.js';
import { zz_pXModulus } from './lzz_pX.js';
import { GCD as wordGCD } from './lzz_pX1.js';
import { SFCanZass1, SFCanZass2 } from './lzz_pXFactoring.js';
interface LocalInfoStorage {
  primes: (bigint | undefined)[];
  primeLength: number;
  patterns: number[][];
  patternLength: number;
}
const localInfoStorage = new WeakMap<LocalInfoT, LocalInfoStorage>();
/** State retained by NTL's integer-polynomial small-prime selection.
 * Array accessors return snapshots; assignment copies the visible prefix while
 * retaining previously initialized storage beyond its new logical length.
 * @see Deviation: NTL stateful small-prime factor selection
 */
export class LocalInfoT {
  n = -1;
  NumPrimes = 0;
  NumFactors = 0;
  PossibleDegrees = 0n;
  readonly s = new PrimeSeq();
  context: Readonly<{ p: bigint; maxroot: number }> | null = null;
  private readonly storage: LocalInfoStorage;
  constructor() {
    this.storage = { primes: [], primeLength: 0, patterns: [], patternLength: 0 };
    localInfoStorage.set(this, this.storage);
  }
  get p(): (bigint | undefined)[] {
    return this.storage.primes.slice(0, this.storage.primeLength);
  }
  set p(values: readonly (bigint | undefined)[]) {
    for (let i = 0; i < values.length; i++) this.storage.primes[i] = values[i];
    this.storage.primeLength = values.length;
  }
  get pattern(): number[][] {
    return this.storage.patterns.slice(0, this.storage.patternLength).map((v) => v.slice());
  }
  set pattern(values: readonly (readonly number[])[]) {
    for (let i = 0; i < values.length; i++) this.storage.patterns[i] = values[i]!.slice();
    this.storage.patternLength = values.length;
  }
}
/** NTL tuning globals represented as options local to this operation. */
export interface SmallPrimeFactorizationOptions {
  InitNumPrimes?: number;
  MaxNumPrimes?: number;
  context?: FFTPrimeContext;
}
/** NTL ZZXFactoring.cpp small-prime selection and ordered factor recovery.
 * Input is a primitive squarefree integer polynomial. Returns null when the
 * modular patterns prove irreducibility, otherwise factors over info.context.p.
 * The prime sequence and admissibility data persist across calls.
 * @see Deviation: NTL stateful small-prime factor selection
 * @see Deviation: NTL cold word-context initialization
 */
export function SmallPrimeFactorization(
  info: LocalInfoT,
  input: readonly bigint[],
  stream: RandomStream,
  options?: SmallPrimeFactorizationOptions
): bigint[][] | null {
  const f = k.normalized(input),
    n = f.length - 1;
  info.n = n;
  info.NumPrimes = 0;
  info.NumFactors = 0;
  const initial = options?.InitNumPrimes ?? 7,
    maximum = options?.MaxNumPrimes ?? 50;
  if (initial < 1 || initial > 10000) throw new Error('bad ZZXFac_InitNumPrimes');
  if (maximum < initial || maximum > 10000) throw new Error('bad ZZXFac_MaxNumPrimes');
  const storage = localInfoStorage.get(info)!;
  storage.primeLength = initial;
  storage.primes.length = Math.max(storage.primes.length, initial);
  storage.patternLength = initial;
  while (storage.patterns.length < initial) storage.patterns.push([]);
  // NTL SetBit edits the magnitude and retains the integer's sign.
  const negativeDegrees = info.PossibleDegrees < 0n;
  const withBit =
    (negativeDegrees ? -info.PossibleDegrees : info.PossibleDegrees) | (1n << BigInt(n + 1));
  info.PossibleDegrees = (negativeDegrees ? -withBit : withBit) - 1n;
  let minr = n + 1,
    bestIndex = -1;
  let bestFactors: [bigint[], number][] = [],
    bestH: bigint[] = [];
  let bestContext: LocalInfoT['context'] = null;
  const maxroot = (n <= 1 ? 0 : BigInt(n - 1).toString(2).length) + 1;
  const wordOptions = {
    maxroot,
    state: options?.context ? { context: options.context, stream } : undefined,
  };
  while (info.NumPrimes < initial) {
    const p = info.s.next();
    if (p === 0n) throw new Error('out of small primes');
    if ((f[n] ?? 0n) % p === 0n) continue;
    const F = new zz_pXModulus(null, p, wordOptions),
      arithmetic = F.arithmetic;
    info.context = Object.freeze({ p, maxroot });
    let ff = arithmetic.norm(f);
    if (ff[ff.length - 1] !== 1n) {
      const scale = arithmetic.inverse(ff[ff.length - 1]!);
      ff = arithmetic.norm(ff.map((v) => v * scale));
    }
    const derivative = arithmetic.norm(ff.slice(1).map((v, i) => v * BigInt(i + 1)));
    const gcd = wordGCD(derivative, ff, p, wordOptions);
    if (gcd.length !== 1 || gcd[0] !== 1n) continue;
    const [factors, h] = SFCanZass1(ff, p, wordOptions);
    storage.primes[info.NumPrimes] = p;
    const pattern = RecordPattern(Array<number>(n + 1), factors, p);
    storage.patterns[info.NumPrimes] = pattern;
    const r = NumFactors(pattern);
    if (r === 1) return null;
    // NTL bit_and returns the intersection of magnitudes, without a sign.
    info.PossibleDegrees =
      (info.PossibleDegrees < 0n ? -info.PossibleDegrees : info.PossibleDegrees) &
      CalcPossibleDegrees(pattern);
    const withoutOne = info.PossibleDegrees & (info.PossibleDegrees - 1n);
    if (withoutOne !== 0n && (withoutOne & (withoutOne - 1n)) === 0n) return null;
    if (r < minr) {
      minr = r;
      bestFactors = factors;
      bestH = h;
      bestContext = info.context;
      bestIndex = info.NumPrimes;
    }
    info.NumPrimes++;
  }
  const last = info.NumPrimes - 1;
  [storage.patterns[bestIndex], storage.patterns[last]] = [
    storage.patterns[last]!,
    storage.patterns[bestIndex]!,
  ];
  storage.primes[bestIndex] = storage.primes[last];
  info.NumPrimes--;
  info.context = bestContext;
  return SFCanZass2(bestFactors, bestH, bestContext!.p, stream, wordOptions);
}

/** NTL tuning globals used while refining integer-factorization information. */
export interface UpdateLocalInfoOptions {
  van_Hoeij?: number;
  MaxNumPrimes?: number;
  state?: PolynomialProductState;
}
/** Refine retained prime patterns after integer factors have been recovered.
 * Return recomputed subset-degree suffixes, or null if that cache was unchanged.
 * Pattern updates preceding an exception remain visible; the saved word context
 * is preserved on both normal and exceptional exits.
 * @see Deviation: NTL local factor-information updates
 * @see Deviation: NTL factor recombination cache state
 */
export function UpdateLocalInfo(
  info: LocalInfoT,
  W: readonly (readonly bigint[])[],
  factors: readonly (readonly bigint[])[],
  input: readonly bigint[],
  count: number,
  modulus: bigint,
  options?: UpdateLocalInfoOptions
): bigint[] | null {
  const storage = localInfoStorage.get(info)!;
  let degrees: bigint[] | null = null;
  if (info.NumFactors < factors.length) {
    if (info.n + 1 < 0) throw new Error('negative length in vector::SetLength');
    const wordOptions = {
      maxroot: (info.n <= 1 ? 0 : BigInt(info.n - 1).toString(2).length) + 1,
      state: options?.state,
    };
    for (let i = 0; i < info.NumPrimes; i++) {
      const p = storage.primes[i]!;
      const F = new zz_pXModulus(null, p, wordOptions);
      for (let j = info.NumFactors; j < factors.length; j++) {
        let ff = F.arithmetic.norm(factors[j]!);
        if (ff.length && ff[ff.length - 1] !== 1n) {
          const inverse = F.arithmetic.inverse(ff[ff.length - 1]!);
          ff = F.arithmetic.norm(ff.map((x) => x * inverse));
        }
        const [fac] = SFCanZass1(ff, p, wordOptions);
        const pattern = RecordPattern(Array<number>(info.n + 1), fac, p);
        const retained = storage.patterns[i]!;
        if (retained.length !== pattern.length) throw new Error('SubPattern: bad args');
        // Native SubPattern retains each subtraction, including the failing slot.
        for (let d = 0; d < retained.length; d++) {
          retained[d] = retained[d]! - pattern[d]!;
          if (retained[d]! < 0) throw new Error('SubPattern: internal error');
        }
      }
      info.PossibleDegrees =
        (info.PossibleDegrees < 0n ? -info.PossibleDegrees : info.PossibleDegrees) &
        CalcPossibleDegrees(storage.patterns[i]!);
    }
    info.NumFactors = factors.length;
    degrees = CalcPossibleDegrees(W, count, modulus);
  }
  if (!(options?.van_Hoeij ?? 1) && info.NumPrimes + 1 < (options?.MaxNumPrimes ?? 50)) {
    const f = k.normalized(input),
      n = f.length - 1;
    const wordOptions = {
      maxroot: (n <= 1 ? 0 : BigInt(n - 1).toString(2).length) + 1,
      state: options?.state,
    };
    for (;;) {
      const p = info.s.next();
      if (!p) throw new Error('UpdateLocalInfo: out of primes');
      if ((f[n] ?? 0n) % p === 0n) continue;
      const F = new zz_pXModulus(null, p, wordOptions),
        arithmetic = F.arithmetic;
      let ff = arithmetic.norm(f);
      if (ff.length && ff[ff.length - 1] !== 1n) {
        const inverse = arithmetic.inverse(ff[ff.length - 1]!);
        ff = arithmetic.norm(ff.map((x) => x * inverse));
      }
      const derivative = arithmetic.norm(ff.slice(1).map((x, i) => x * BigInt(i + 1)));
      const gcd = wordGCD(derivative, ff, p, wordOptions);
      if (gcd.length !== 1 || gcd[0] !== 1n) continue;
      const [fac] = SFCanZass1(ff, p, wordOptions);
      const length = info.NumPrimes + 1;
      if (length < 0) throw new Error('negative length in vector::SetLength');
      storage.primeLength = length;
      storage.primes.length = Math.max(storage.primes.length, length);
      storage.patternLength = length;
      while (storage.patterns.length < length) storage.patterns.push([]);
      storage.primes[info.NumPrimes] = p;
      if (info.n + 1 < 0) throw new Error('negative length in vector::SetLength');
      const pattern = RecordPattern(Array<number>(info.n + 1), fac, p);
      storage.patterns[info.NumPrimes] = pattern;
      info.PossibleDegrees =
        (info.PossibleDegrees < 0n ? -info.PossibleDegrees : info.PossibleDegrees) &
        CalcPossibleDegrees(pattern);
      info.NumPrimes++;
      break;
    }
  }
  return degrees;
}

/** Native integer-factor recombination by subset cardinality.
 * The returned tuple contains recovered factors, the remaining integer
 * polynomial and the remaining modular factors; inputs are copied.
 * @see Deviation: NTL integer-factor cardinality search
 * @see Deviation: NTL factor recombination cache state
 */
export function CardinalitySearch(
  recovered: readonly (readonly bigint[])[],
  input: readonly bigint[],
  modular: readonly (readonly bigint[])[],
  info: LocalInfoT,
  count: number,
  bound: number,
  modulus: bigint,
  options?: UpdateLocalInfoOptions
): [bigint[][], bigint[], bigint[][]] {
  if (count < 0) throw new Error('negative length in vector::SetLength');
  if (count >= 2 ** 57) throw new Error('excessive length in vector::SetLength');
  const arithmetic = _ZZ_pX_euclidean_kernels(modulus);
  const factors = recovered.map(k.normalized);
  let f = k.normalized(input),
    W = modular.map(arithmetic.norm);
  const I = Array<number>(count).fill(0),
    D = Array<number>(count).fill(0);
  let prod = Array<bigint>(count).fill(0n),
    pdeg = CalcPossibleDegrees(W, count, modulus);
  if (count === 0 && W.length === 0)
    throw new Error(
      'CardinalitySearch: zero cardinality with no factors has no defined native result'
    );
  let r = W.length,
    counter = 0;
  let ct = (f[0] ?? 0n) * (f.at(-1) ?? 0n),
    lc = mod(f.at(-1) ?? 0n, modulus);
  I[0] = 0;
  done: while (I[0]! <= r - count) {
    let pd =
      pdeg[I[0]!]! & (info.PossibleDegrees < 0n ? -info.PossibleDegrees : info.PossibleDegrees);
    if (!pd) break;
    let allowed = unpack(pd, info.n);
    D[0] = W[I[0]!]!.length - 1;
    let i = 1,
      state = 0,
      ProdLen = 0;
    for (;;) {
      if (i < ProdLen) ProdLen = i;
      if (i === count) {
        if (counter > 2000000) {
          counter = 0;
          pdeg = UpdateLocalInfo(info, W, factors, f, count, modulus, options) ?? pdeg;
          pd =
            pdeg[I[0]!]! &
            (info.PossibleDegrees < 0n ? -info.PossibleDegrees : info.PossibleDegrees);
          if (!pd) break done;
          allowed = unpack(pd, info.n);
        }
        state = 1;
        if (!allowed[D[count - 1]!]) {
          i--;
          counter++;
          continue;
        }
        let passes: number;
        [passes, prod, ProdLen] = ConstTermTest(W, I, ct, lc, prod, ProdLen, modulus);
        if (!passes) {
          i--;
          counter += 100;
          continue;
        }
        counter += 1000;
        const selected = 2 * D[count - 1]! <= f.length - 1;
        const gg = selected
          ? mul(W, modulus, I, options?.state)
          : InvMul(W, I, modulus, options?.state);
        let g = BalCopy(arithmetic.norm(gg.map((x) => x * lc)), modulus);
        if (g.reduce((b, x) => Math.max(b, bits(x)), 0) > bound) {
          i--;
          continue;
        }
        g = k.primitive(g);
        const h = k.divide(f, g, options?.state);
        if (h === null) {
          i--;
          continue;
        }
        factors.push(selected ? g : h);
        f = selected ? h : g;
        ct = (f[0] ?? 0n) * (f.at(-1) ?? 0n);
        lc = mod(f.at(-1) ?? 0n, modulus);
        W = RemoveFactors(W, I, modulus);
        r = W.length;
        counter = 0;
        if (2 * count > r) break done;
        break;
      } else if (state === 0) {
        I[i] = I[i - 1]! + 1;
        D[i] = D[i - 1]! + W[I[i]!]!.length - 1;
        i++;
      } else {
        I[i] = I[i]! + 1;
        if (i === 0) break;
        if (I[i]! > r - count + i) i--;
        else {
          D[i] = D[i - 1]! + W[I[i]!]!.length - 1;
          i++;
          state = 0;
        }
      }
    }
  }
  return [factors, f, W];
}

// Original ZZXFactoring.cpp1251-1431,1620-1658, native unsigned 64-bit profile.
const word64 = (x: bigint) => BigInt.asUintN(64, x);
type LookupTable = bigint[][][];
type ShiftTable = number[][];
function pruning_bnd(r: number, k: number): number {
  let x = 0;
  for (let i = 0; i < k; i++) x += Math.log((r - i) / (k - i));
  return Math.trunc((x / Math.log(2)) * 0.75);
}
function shamt_tab_init(pos: number, card: number, pruning: number, thresholdBits: number): number {
  let x = 1;
  for (let i = 0; i < card; i++) x *= (pos - i) / (card - i);
  x *= pruning;
  if (pos <= 6) x *= 2;
  const t = Math.min(Math.max(Math.ceil(Math.log(x) / Math.log(2)), 6), 64 - thresholdBits);
  return 64 - t;
}
function DoInitTab(
  table: LookupTable,
  limit: number,
  ratio: bigint[],
  r: number,
  k: number,
  threshold: bigint,
  shifts: ShiftTable
): void {
  const sums = Array<bigint>(limit + 1).fill(0n),
    cards = Array<number>(limit + 1).fill(0),
    locations = Array<number>(limit + 1).fill(0);
  let j = 0,
    location = 0;
  while (j >= 0) {
    const sum = sums[j]!,
      card = cards[j]!;
    if (location === 0) {
      if (j >= limit || card >= k - 1) {
        if (card > 1) {
          const shift = BigInt(shifts[limit]![card]!);
          const first = word64(-sum) >> shift;
          const bucket = Number(first >> 6n);
          table[limit]![card]![bucket] = table[limit]![card]![bucket]! | (1n << (first & 63n));
          const second = word64(-sum + threshold) >> shift;
          if (first !== second) {
            const nextBucket = Number(second >> 6n);
            table[limit]![card]![nextBucket] =
              table[limit]![card]![nextBucket]! | (1n << (second & 63n));
          }
        }
        location = locations[j]!;
        j--;
        continue;
      }
      sums[j + 1] = sum;
      cards[j + 1] = card;
      locations[j + 1] = 1;
      j++;
      location = 0;
    } else if (location === 1) {
      sums[j + 1] = word64(sum + ratio[r - 1 - j]!);
      cards[j + 1] = card + 1;
      locations[j + 1] = 2;
      j++;
      location = 0;
    } else {
      location = locations[j]!;
      j--;
    }
  }
}
function InitTab(
  table: LookupTable,
  ratio: bigint[],
  r: number,
  k: number,
  threshold: bigint,
  shifts: ShiftTable,
  pruning: number
): void {
  if (!pruning) return;
  for (let i = 2; i <= pruning; i++) {
    for (let j = 2; j <= Math.min(k - 1, i); j++) {
      const length = Number(((1n << BigInt(64 - shifts[i]![j]!)) + 63n) >> 6n);
      for (let t = 0; t < length; t++) table[i]![j]![t] = 0n;
    }
    DoInitTab(table, i, ratio, r, k, threshold, shifts);
  }
}
function RatioInit1(
  W: bigint[][],
  lc: bigint,
  pruning: number,
  table: LookupTable,
  k: number,
  threshold: bigint,
  shifts: ShiftTable,
  modulus: bigint
): [bigint[], bigint[][]] {
  const scaled = (x: bigint) => (mod(x, modulus) << 64n) / modulus;
  const ratio = W.map((f) => scaled(f[f.length - 2]! * lc));
  InitTab(table, ratio, W.length, k, threshold, shifts, pruning);
  const pair = W.map(() => Array<bigint>(W.length).fill(0n));
  for (let i = 0; i < W.length; i++) {
    for (let j = 0; j < i; j++)
      pair[i]![j] = scaled(W[i]![W[i]!.length - 2]! * W[j]![W[j]!.length - 2]! * lc);
    pair[i]![i] = W[i]!.length >= 3 ? scaled(W[i]![W[i]!.length - 3]! * lc) : 0n;
  }
  return [ratio, pair];
}
function SecondOrderTest(
  I: number[],
  pair: bigint[][],
  stack: bigint[],
  length: number
): [boolean, number] {
  const k = I.length;
  let sum: bigint, threshold: bigint;
  if (length === 0) {
    const delta = (BigInt(k) * BigInt(k + 1)) >> 1n;
    sum = 1n + delta;
    threshold = 2n + delta;
    stack[k] = threshold;
  } else {
    sum = stack[length - 1]!;
    threshold = stack[k]!;
  }
  for (let i = length; i < k; i++) {
    const p = pair[I[i]!]!;
    for (let j = 0; j <= i; j++) sum = word64(sum + p[I[j]!]!);
    stack[i] = sum;
  }
  return [sum <= threshold, k - 1];
}

/** Operation-local NTL factor-recombination tuning values. */
export interface FactorRecombinationOptions extends UpdateLocalInfoOptions {
  MaxPrune?: number;
}
/** Native ratio-pruned integer-factor recombination for cardinalities above one.
 * Returns copied recovered factors, remaining integer polynomial and modular factors.
 * @see Deviation: NTL integer-factor cardinality search
 * @see Deviation: NTL factor recombination cache state
 */
export function CardinalitySearch1(
  recovered: readonly (readonly bigint[])[],
  input: readonly bigint[],
  modular: readonly (readonly bigint[])[],
  info: LocalInfoT,
  count: number,
  bound: number,
  modulus: bigint,
  options?: FactorRecombinationOptions
): [bigint[][], bigint[], bigint[][]] {
  if (count <= 1) throw new Error('internal error: call CardinalitySearch');
  if (bits(BigInt(count)) > 30) throw new Error('Cardinality Search: k too large...');
  const arithmetic = _ZZ_pX_euclidean_kernels(modulus);
  const factors = recovered.map(k.normalized);
  let f = k.normalized(input),
    W = modular.map(arithmetic.norm);
  let pdeg = CalcPossibleDegrees(W, count, modulus);
  if (!W.length)
    throw new Error('CardinalitySearch1: empty factor vector has no defined native result');
  const intersection = (x: bigint) =>
    x & (info.PossibleDegrees < 0n ? -info.PossibleDegrees : info.PossibleDegrees);
  if (!intersection(pdeg[0]!)) return [factors, f, W];
  let r = W.length;
  const I = Array<number>(count).fill(0),
    D = Array<number>(count).fill(0),
    ratioSum = Array<bigint>(count).fill(0n);
  const threshold = BigInt(count) + 1n,
    threshold1 = BigInt(count) + 2n;
  const thresholdBits = bits(threshold1);
  let pruning = Math.min(
    Math.floor(r / 2),
    options?.MaxPrune ?? 10,
    pruning_bnd(r, count),
    63 - thresholdBits
  );
  if (pruning <= 4) pruning = 0;
  const shifts: ShiftTable = Array.from({ length: pruning + 1 }, () => []);
  const table: LookupTable = Array.from({ length: pruning + 1 }, () => []);
  for (let i = 2; i <= pruning; i++) {
    shifts[i] = Array<number>(Math.min(count - 1, i) + 1).fill(0);
    table[i] = Array.from({ length: Math.min(count - 1, i) + 1 }, () => []);
    for (let j = 2; j <= Math.min(count - 1, i); j++) {
      shifts[i]![j] = shamt_tab_init(i, j, pruning, thresholdBits);
      table[i]![j] = Array<bigint>(Number(((1n << BigInt(64 - shifts[i]![j]!)) + 63n) >> 6n)).fill(
        0n
      );
    }
  }
  let prod = Array<bigint>(count).fill(0n),
    prod1 = Array<bigint>(count).fill(0n);
  const sumStack = Array<bigint>(count + 1).fill(0n);
  let counter = 0;
  let ct = (f[0] ?? 0n) * (f.at(-1) ?? 0n),
    lc = mod(f.at(-1) ?? 0n, modulus);
  let [ratio, pair] = RatioInit1(W, lc, pruning, table, count, threshold1, shifts, modulus);
  // These native filters are initialized once and retained after factor removal.
  const c1 = f.reduce((a, b) => a + b, 0n) * (f.at(-1) ?? 0n);
  let sums = W.map((a) => [
    mod(
      a.reduce((x, y) => x + y, 0n),
      modulus
    ),
  ]);
  let degree = W.map((a) => a.length - 1);
  I[0] = 0;
  done: while (I[0]! <= r - count) {
    let pd = intersection(pdeg[I[0]!]!);
    if (!pd) break;
    let allowed = unpack(pd, info.n);
    D[0] = degree[I[0]!]!;
    ratioSum[0] = word64(ratio[I[0]!]! + threshold);
    let i = 1,
      state = 0,
      ProdLen = 0,
      ProdLen1 = 0,
      SumLen = 0;
    let restart = false;
    for (;;) {
      counter++;
      if (counter > 2000000) {
        counter = 0;
        pdeg = UpdateLocalInfo(info, W, factors, f, count, modulus, options) ?? pdeg;
        pd = intersection(pdeg[I[0]!]!);
        if (!pd) break done;
        allowed = unpack(pd, info.n);
      }
      if (i === count - 1) {
        const previousRatio = ratioSum[count - 2]!,
          previousIndex = I[count - 2]!,
          previousDegree = D[count - 2]!;
        for (let next = previousIndex + 1; next < r; next++) {
          const rs = word64(previousRatio + ratio[next]!);
          if (rs > threshold1) {
            counter++;
            continue;
          }
          const candidateDegree = previousDegree + degree[next]!;
          if (!allowed[candidateDegree]) {
            counter++;
            continue;
          }
          I[count - 1] = next;
          let secondOrder: boolean;
          [secondOrder, SumLen] = SecondOrderTest(I, pair, sumStack, SumLen);
          if (!secondOrder) {
            counter += 2;
            continue;
          }
          let passes: number;
          [passes, prod1, ProdLen1] = ConstTermTest(sums, I, c1, lc, prod1, ProdLen1, modulus);
          if (!passes) {
            counter += 100;
            continue;
          }
          D[count - 1] = candidateDegree;
          [passes, prod, ProdLen] = ConstTermTest(W, I, ct, lc, prod, ProdLen, modulus);
          if (!passes) {
            counter += 100;
            continue;
          }
          counter += 1000;
          const selected = 2 * candidateDegree <= f.length - 1;
          const gg = selected
            ? mul(W, modulus, I, options?.state)
            : InvMul(W, I, modulus, options?.state);
          let g = BalCopy(arithmetic.norm(gg.map((x) => x * lc)), modulus);
          if (g.reduce((b, x) => Math.max(b, bits(x)), 0) > bound) continue;
          g = k.primitive(g);
          const h = k.divide(f, g, options?.state);
          if (h === null) continue;
          factors.push(selected ? g : h);
          f = selected ? h : g;
          ct = (f[0] ?? 0n) * (f.at(-1) ?? 0n);
          lc = mod(f.at(-1) ?? 0n, modulus);
          W = RemoveFactors(W, I, modulus);
          const removed = new Set(I);
          degree = degree.filter((_, j) => !removed.has(j));
          sums = sums.filter((_, j) => !removed.has(j));
          ratio = ratio.filter((_, j) => !removed.has(j));
          pair = pair
            .filter((_, j) => !removed.has(j))
            .map((row) => row.filter((_, j) => !removed.has(j)));
          r = W.length;
          counter = 0;
          pruning = Math.min(pruning, Math.floor(r / 2));
          if (pruning <= 4) pruning = 0;
          InitTab(table, ratio, r, count, threshold1, shifts, pruning);
          if (2 * count > r) break done;
          restart = true;
          break;
        }
        if (restart) break;
        i--;
        state = 1;
      } else {
        const pruned = (index: number) => {
          if (!pruning || r - index > pruning) return false;
          const pos = r - index,
            encoded = ratioSum[i - 1]! >> BigInt(shifts[pos]![count - i]!);
          return !(table[pos]![count - i]![Number(encoded >> 6n)]! & (1n << (encoded & 63n)));
        };
        if (state === 0) {
          const index = I[i - 1]! + 1;
          I[i] = index;
          if (pruned(index)) {
            i--;
            state = 1;
          } else {
            D[i] = D[i - 1]! + degree[index]!;
            ratioSum[i] = word64(ratioSum[i - 1]! + ratio[index]!);
            i++;
          }
        } else {
          if (i < ProdLen) ProdLen = i;
          if (i < ProdLen1) ProdLen1 = i;
          if (i < SumLen) SumLen = i;
          const index = I[i]! + 1;
          I[i] = index;
          if (i === 0) break;
          if (index > r - count + i) i--;
          else if (pruned(index)) i--;
          else {
            D[i] = D[i - 1]! + degree[index]!;
            ratioSum[i] = word64(ratioSum[i - 1]! + ratio[index]!);
            i++;
            state = 0;
          }
        }
      }
    }
  }
  return [factors, f, W];
}

/** Native cardinality-ordered integer-factor recovery from lifted factors.
 * @see Deviation: NTL integer-factor cardinality search
 * @see Deviation: NTL factor recombination cache state
 */
export function FindTrueFactors(
  input: readonly bigint[],
  lifted: readonly (readonly bigint[])[],
  modulus: bigint,
  info: LocalInfoT,
  bound: number,
  options?: FactorRecombinationOptions
): bigint[][] {
  const arithmetic = _ZZ_pX_euclidean_kernels(modulus);
  let factors: bigint[][] = [],
    f = k.normalized(input),
    W = lifted.map(arithmetic.norm);
  for (let count = 1; 2 * count <= W.length; count++) {
    [factors, f, W] =
      count <= 1
        ? CardinalitySearch(factors, f, W, info, count, bound, modulus, options)
        : CardinalitySearch1(factors, f, W, info, count, bound, modulus, options);
  }
  factors.push(f);
  return factors;
}

import { GenPrime_long } from './ZZ.js';
import { FFTPrimeContext, UseFFTPrime, GetFFTPrime } from './FFT.js';
import { gauss as wordGauss } from './mat_lzz_p.js';
import { inv as integerInverse } from './mat_ZZ.js';
/** Native certified row reduction for nonempty matrices with independent rows.
 * Returns [d, R] so R/d is the reduced row echelon form.
 * @see Deviation: NTL certified integer elimination
 */
export function gauss(
  M: readonly (readonly bigint[])[],
  context: FFTPrimeContext,
  stream: RandomStream,
  options: { columns?: number } = {}
): [bigint, bigint[][]] {
  const n = M.length,
    m = options.columns ?? M[0]?.length ?? 0;
  if (m < 0) throw new Error('SetDims: bad args');
  if (M.some((row) => row.length !== m)) throw new Error('nonrectangular matrix');
  if (!n || !m) throw new Error('gauss: internal error');
  const A = M.map((row) => [...row]);
  for (;;) {
    const p = GenPrime_long(60, stream);
    // Observable part of zz_p::init(p): initialize the shared CRT prime cache.
    // Exact BigInt arithmetic does not use native word preconditioner metadata.
    const bound = (p * p) << 29n;
    let product = 1n;
    for (let i = 0; product <= bound; i++) {
      UseFFTPrime(i, context, stream);
      product *= GetFFTPrime(i, context);
    }
    const [rank, MM] = wordGauss(A, p, { columns: m });
    if (rank < n) continue;
    const pos: number[] = [];
    for (let i = 0, j = 0; i < n; i++) {
      while (MM[i]![j] === 0n) j++;
      pos.push(j++);
    }
    const S = A.map((row) => pos.map((j) => row[j]!));
    const [d, S_inv] = integerInverse(S, context, stream, { status: true });
    if (!d) continue;
    const R = S_inv.map((row) =>
      Array.from({ length: m }, (_, j) => row.reduce((s, v, k) => s + v * A[k]![j]!, 0n))
    );
    let ok = true;
    for (let i = 0; i < n && ok; i++) {
      for (let j = 0; j < pos[i]! && ok; j++) if (R[i]![j] !== 0n) ok = false;
      if (R[i]![pos[i]!] !== d) ok = false;
      for (let j = 0; j < i && ok; j++) if (R[j]![pos[i]!] !== 0n) ok = false;
    }
    if (!ok) continue;
    return [d, R];
  }
}

/** Native van Hoeij recombination acceptance check with retained factor prefix.
 * @see Deviation: NTL van Hoeij factor acceptance
 */
export function GotThem(
  previous: readonly (readonly bigint[])[],
  B_L: readonly (readonly bigint[])[],
  lifted: readonly (readonly bigint[])[],
  input: readonly bigint[],
  bound: number,
  modulus: bigint,
  context: FFTPrimeContext,
  stream: RandomStream,
  options: { columns?: number } = {}
): [number, bigint[][]] {
  if (!Number.isInteger(bound) || bound < -(2 ** 63) || bound >= 2 ** 63)
    throw new RangeError('GotThem bound must be a signed native integer');
  if (modulus <= 1n) throw new Error('ZZ_pContext: p must be > 1');
  const arithmetic = _ZZ_pX_euclidean_kernels(modulus),
    f = k.normalized(input),
    W = lifted.map(arithmetic.norm),
    factors = previous.map(k.normalized);
  const [det, R] = gauss(B_L, context, stream, options),
    s = B_L.length,
    r = options.columns ?? B_L[0]?.length ?? 0;
  for (let j = 0; j < r; j++) {
    let count = 0;
    for (let i = 0; i < s; i++) {
      if (R[i]![j] === 0n) continue;
      if (R[i]![j] !== det) return [0, factors];
      count++;
    }
    if (count !== 1) return [0, factors];
  }
  const groups: number[][] = [];
  for (let i = 0; i < s; i++) {
    const I: number[] = [];
    for (let j = 0; j < r; j++) if (R[i]![j] !== 0n) I.push(j);
    groups.push(I);
  }
  for (const I of groups) if (I.length <= 3) return [0, factors];
  // The original adds every W[j] degree to every row's degree key, outside the
  // membership test. All keys are equal: its sort preserves certified row order.
  const ct = (f.at(-1) ?? 0n) * (f[0] ?? 0n),
    half = modulus >> 1n,
    lc = mod(f.at(-1) ?? 0n, modulus);
  for (const I of groups) {
    let prod = lc;
    for (const j of I) prod = (prod * (W[j]![0] ?? 0n)) % modulus;
    if (prod > half) prod -= modulus;
    if (prod === 0n ? ct !== 0n : ct % prod !== 0n) return [0, factors];
  }
  const fac: bigint[][] = [];
  for (let i = 0; i < s - 1; i++) {
    const product = mul(W, modulus, groups[i], { context, stream });
    let g = BalCopy(arithmetic.norm(product.map((x) => x * lc)), modulus);
    if (g.reduce((out, x) => Math.max(out, bits(x)), 0) > bound) return [0, factors];
    g = k.primitive(g);
    fac.push(g);
  }
  let remaining = f;
  for (const candidate of fac) {
    const quotient = k.divide(remaining, candidate, { context, stream });
    if (quotient === null) {
      console.error('X');
      return [0, factors];
    }
    remaining = quotient;
  }
  factors.push(...fac, remaining);
  return [1, factors];
}

/** Native tuning globals, scoped to one integer factorization operation.
 * @see Deviation: NTL complete integer factorization driver
 */
export interface IntegerFactorizationOptions extends FactorRecombinationOptions {
  InitNumPrimes?: number;
  PowerHack?: number;
  /** @internal Native early-abandon flag for intermediate deflation stages. */
  ok_to_abandon?: number;
}

/** Native adaptive van Hoeij recombination after bounded cardinality search.
 * @see Deviation: NTL complete integer factorization driver
 */
export function FindTrueFactors_vH(
  input: readonly bigint[],
  lifted: readonly (readonly bigint[])[],
  modulus: bigint,
  p: bigint,
  exponent: number,
  info: LocalInfoT,
  bound: number,
  options?: IntegerFactorizationOptions
): bigint[][] {
  const state = options?.state ?? {
    context: new FFTPrimeContext(),
    stream: new RandomStream(new Uint8Array(32)),
  };
  const settings = { ...options, state };
  const arithmetic = _ZZ_pX_euclidean_kernels(modulus);
  let factors: bigint[][] = [],
    f = k.normalized(input),
    W = lifted.map(arithmetic.norm);
  let count = 1;
  for (; 2 * count <= W.length && (count <= 3 || W.length <= 12); count++) {
    [factors, f, W] =
      count <= 1
        ? CardinalitySearch(factors, f, W, info, count, bound, modulus, settings)
        : CardinalitySearch1(factors, f, W, info, count, bound, modulus, settings);
  }
  if (2 * count > W.length) {
    factors.push(f);
    return factors;
  }
  let P1 = modulus,
    e1 = exponent,
    w1 = W.map((a) => a.slice());
  const r = W.length,
    n = f.length - 1,
    root_bound = RootBound(f);
  let B_L: bigint[][] = Array.from({ length: r }, (_, i) =>
    Array.from({ length: r }, (_, j) => (i === j ? 1n : 0n))
  );
  let d = 0,
    bit_delta = 0,
    b: number[] = [],
    pb: bigint[] = [],
    delta = 0,
    pdelta = 1n;
  let trace_vec: bigint[][] = Array.from({ length: r }, () => []);
  const chop_vec: bigint[][] = Array.from({ length: r }, () => []);
  let dense = false,
    s = r;
  const ran_bits = 32;
  for (;;) {
    if (
      settings.ok_to_abandon &&
      ((d >= 2 && s > 128) || (d >= 3 && s > 32) || (d >= 4 && s > 8) || d >= 5)
    ) {
      factors.push(f);
      break;
    }
    const d_last = d,
      d_inc = Math.min(1 + Math.trunc(d / (dense ? 4 : 8)), n - 1 - d);
    d += d_inc;
    if (bit_delta === 0) bit_delta = 2 * r;
    else
      bit_delta +=
        !dense || d_inc === 0
          ? 1 + Math.trunc(bit_delta / 8)
          : d1_val(bit_delta, r, s) > 1
            ? 1 + Math.trunc(bit_delta / 16)
            : 0;
    if (d > d1_val(bit_delta, r, s)) dense = true;
    [delta, pdelta] = Compute_pdelta(delta, pdelta, p, bit_delta);
    let d1: number, b_eff: number, pb_eff: bigint;
    if (!dense) {
      for (let index = d_last + 1; index <= d; index++)
        [b, pb] = Compute_pb(b, pb, p, index, root_bound, n);
      d1 = d;
      b_eff = b[d - 1]!;
      pb_eff = pb[d - 1]!;
    } else {
      d1 = d1_val(bit_delta, r, s);
      [b_eff, pb_eff] = Compute_pb_eff(p, d, root_bound, n, ran_bits);
    }
    if (b_eff + delta > e1) {
      [P1, e1, w1] = AdditionalLifting(P1, e1, w1, p, b_eff + delta, f, true, state);
      trace_vec = Array.from({ length: r }, (_, i) => {
        let trace = Array<bigint>(d_last).fill(0n);
        for (let index = 1; index <= d_last; index++)
          trace = ComputeTrace(trace, w1[i]!, index, P1);
        return trace;
      });
    }
    const A = dense
      ? Array.from({ length: d1 }, () =>
          Array.from({ length: d }, () => {
            const value = RandomBits(32, state.stream);
            return RandomBnd(2n, state.stream, { word: true }) ? -value : value;
          })
        )
      : [];
    for (let i = 0; i < r; i++) {
      let trace = Array.from({ length: d }, (_, j) => trace_vec[i]![j] ?? 0n);
      for (let index = d_last + 1; index <= d; index++)
        trace = ComputeTrace(trace, w1[i]!, index, P1);
      trace_vec[i] = trace;
      const chop = Array.from({ length: d1 }, (_, j) => chop_vec[i]![j] ?? 0n);
      chop_vec[i] = dense
        ? DenseChopTraces(chop, trace, d, d1, pb_eff, pdelta, P1, f.at(-1)!, A)
        : ChopTraces(chop, trace, d, pb, pdelta, P1, f.at(-1)!);
    }
    const [M, C] = BuildReductionMatrix(r, d1, pdelta, chop_vec, B_L);
    const [rank, D, reduced] = ntl_LLL_plus(M);
    if (rank !== s + d1) throw new Error('van Hoeij -- bad rank');
    const B1 = CutAway(D, reduced, C, r, d1);
    if (B1.length >= s) continue;
    B_L = B1;
    s = B_L.length;
    if (s === 0) throw new Error('oops! s == 0 should not happen!');
    if (s === 1) {
      factors.push(f);
      break;
    }
    if (s > r / 4) continue;
    const [done, next] = GotThem(factors, B_L, W, f, bound, modulus, state.context, state.stream);
    factors = next;
    if (done) break;
  }
  return factors;
}

/** Native primitive, squarefree factor driver before the deflation wrapper.
 * @see Deviation: NTL complete integer factorization driver
 */
export function ll_SFFactor(
  input: readonly bigint[],
  bound = 0,
  options?: IntegerFactorizationOptions
): bigint[][] {
  const ff = k.normalized(input);
  if (ff.length <= 2) return [ff];
  const state = options?.state ?? {
    context: new FFTPrimeContext(),
    stream: new RandomStream(new Uint8Array(32)),
  };
  const settings = { ...options, state };
  const xfac = ff[0] === 0n;
  let f = xfac ? ff.slice(1) : ff,
    sum = f.reduce((a, b) => a + b, 0n);
  const x1fac = sum === 0n;
  if (x1fac) {
    // ZZX1.cpp PlainPseudoDivRem with db=1, LC=1 and b[0]=-1.
    const remainder = f.slice(), quotient = Array<bigint>(f.length - 1);
    for (let i = quotient.length - 1; i >= 0; i--) {
      quotient[i] = remainder[i + 1]!;
      remainder[i] = remainder[i]! + quotient[i]!;
    }
    f = quotient;
  }
  sum = f.reduce((a, b) => a + b, 0n);
  const appendTrivial = (factors: bigint[][]): bigint[][] => {
    if (xfac) factors.push([0n, 1n]);
    if (x1fac) factors.push([-1n, 1n]);
    return factors;
  };
  if (f.length <= 2) return appendTrivial(f.length > 1 ? [f] : []);
  const abs = (x: bigint) => (x < 0n ? -x : x);
  const rev = abs(f.at(-1)!) > abs(f[0]!);
  if (rev) f = inplace_rev(f);
  const info = new LocalInfoT();
  const small = SmallPrimeFactorization(info, f, state.stream, {
    InitNumPrimes: settings.InitNumPrimes,
    MaxNumPrimes: settings.MaxNumPrimes,
    context: state.context,
  });
  if (small === null) return appendTrivial([rev ? inplace_rev(f) : f]);
  const n = f.length - 1;
  const bnd1 =
    f.reduce((best, c) => Math.max(best, bits(c)), 0) + Math.trunc((bits(BigInt(n + 1)) + 1) / 2);
  if (!bound || bnd1 < bound) bound = bnd1;
  let i = Math.trunc(n / 2);
  while (!((info.PossibleDegrees >> BigInt(i)) & 1n)) i--;
  const lc_bnd = bits(f.at(-1)!),
    coeff_bnd = bound + lc_bnd + i;
  const lift_bnd =
    Math.max(coeff_bnd + 15, bound + lc_bnd + 2 * bits(BigInt(n)) + 64, lc_bnd + bits(sum)) + 2;
  const p = info.context!.p;
  let exponent = Math.trunc(lift_bnd / (Math.log(Number(p)) / Math.log(2))),
    modulus = p ** BigInt(exponent);
  while (bits(modulus) <= lift_bnd) {
    modulus *= p;
    exponent++;
  }
  let target: bigint[];
  if (f.at(-1) === 1n) target = f;
  else if (f.at(-1) === -1n) target = f.map((x) => -x);
  else {
    const ar = _ZZ_pX_euclidean_kernels(modulus),
      inverse = ar.inverse(ar.mod(f.at(-1)!));
    target = f.map((x) => ar.mod(x * inverse));
  }
  const lifted = MultiLift(small, target, exponent, p, { maxroot: info.context!.maxroot, state });
  const factors =
    (settings.van_Hoeij ?? 1) && lifted.length > 12
      ? FindTrueFactors_vH(f, lifted, modulus, p, exponent, info, coeff_bnd, settings)
      : FindTrueFactors(f, lifted, modulus, info, coeff_bnd, settings);
  if (rev) {
    for (let j = 0; j < factors.length; j++) {
      factors[j] = inplace_rev(factors[j]!);
      if (factors[j]!.at(-1)! < 0n) factors[j] = factors[j]!.map((x) => -x);
    }
  }
  return appendTrivial(factors);
}

/** Native squarefree integer factorization, including recursive deflation.
 * @see Deviation: NTL complete integer factorization driver
 */
export function SFFactor(
  input: readonly bigint[],
  bound = 0,
  options?: IntegerFactorizationOptions
): bigint[][] {
  const f = k.normalized(input);
  if (!f.length) throw new Error('SFFactor: bad args');
  if (f.length <= 1) return [];
  const state = options?.state ?? {
    context: new FFTPrimeContext(),
    stream: new RandomStream(new Uint8Array(32)),
  };
  const settings = { ...options, state };
  if (!(settings.PowerHack ?? 1)) return ll_SFFactor(f, bound, { ...settings, ok_to_abandon: 0 });
  let m = 0;
  for (let i = 1; i < f.length && m !== 1; i++)
    if (f[i] !== 0n) {
      let a = m,
        b = i;
      while (b) {
        const t = a % b;
        a = b;
        b = t;
      }
      m = a;
    }
  if (m === 1) return ll_SFFactor(f, bound, { ...settings, ok_to_abandon: 0 });
  const v: number[] = [];
  let rest = m;
  for (let p = 2; rest > 1; p++)
    while (rest % p === 0) {
      v.push(p);
      rest /= p;
    }
  let result = [Array.from({ length: Math.trunc((f.length - 1) / m) + 1 }, (_, i) => f[i * m]!)];
  for (let j = v.length - 1; ; j--) {
    const next: bigint[][] = [];
    for (const part of result)
      next.push(
        ...ll_SFFactor(part, j < 0 ? bound : 0, { ...settings, ok_to_abandon: j < 0 ? 0 : 1 })
      );
    if (j < 0) return next;
    result = next.map((part) => {
      const out = Array<bigint>((part.length - 1) * v[j]! + 1).fill(0n);
      for (let i = 0; i < part.length; i++) out[i * v[j]!] = part[i]!;
      return out;
    });
  }
}

/** Native content and multiplicity-preserving integer factorization.
 * @see Deviation: NTL complete integer factorization driver
 */
export function factor(
  input: readonly bigint[],
  bound = 0,
  options?: IntegerFactorizationOptions
): [bigint, Array<[bigint[], number]>] {
  const f = k.normalized(input);
  if (f.length <= 1) return [f[0] ?? 0n, []];
  const state = options?.state ?? {
    context: new FFTPrimeContext(),
    stream: new RandomStream(new Uint8Array(32)),
  };
  const settings = { ...options, state },
    content = k.content(f),
    ff = f.map((x) => x / content);
  const bnd1 =
    ff.reduce((best, c) => Math.max(best, bits(c)), 0) +
    Math.trunc((bits(BigInt(ff.length)) + 1) / 2);
  if (!bound || bound > bnd1) bound = bnd1;
  const result: Array<[bigint[], number]> = [];
  for (const [part, multiplicity] of SquareFreeDecomp(ff, state))
    for (const factor of SFFactor(part, bound, settings)) result.push([factor, multiplicity]);
  return [content, result];
}
