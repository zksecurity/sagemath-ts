import type { zz_pXOptions } from './lzz_pX.js';
import { RandomBnd, type RandomStream } from './ZZ.js';
import {
  zz_pXModulus,
  zz_pXMultiplier,
  MulByXMod,
  mul,
  MulMod,
  rem,
  _zz_pX_transform_kernels,
} from './lzz_pX.js';
import { mul as matrix_mul } from './mat_lzz_p.js';
/** NTL block-composition argument, populated by build.
 * @see Deviation: NTL word polynomial composition adapters
 */
export class zz_pXNewArgument {
  mat: bigint[][] = [];
  poly: bigint[] = [];
}
const sqrtCount = (n: number): number => {
  const a = BigInt(n);
  if (a < 2n) return n;
  let x = 1n << BigInt(Math.ceil(a.toString(2).length / 2));
  for (;;) {
    const y = (x + a / x) >> 1n;
    if (y >= x) return Number(x);
    x = y;
  }
};
/** NTL matrix of powers and the next block power; commits cache state on success.
 * @see Deviation: NTL word polynomial composition adapters
 */
export function build(H: zz_pXNewArgument, h: readonly bigint[], F: zz_pXModulus, m: number): void {
  const k = F.arithmetic,
    base = k.norm(h),
    n = F.n,
    dh = base.length - 1;
  if (m <= 0 || dh >= n) throw new Error('build: bad args');
  if (m >= 2 ** 60) throw new Error('zz_pXNewArgument:build: m too big');
  let width: number;
  if (dh <= 0 || m === 1) width = 1;
  else if (dh <= Math.floor(n / (m - 1))) width = Math.min(n, dh * (m - 1) + 1);
  else width = n;
  const M = new zz_pXMultiplier(base, F);
  const mat = Array.from({ length: m }, () => Array<bigint>(width).fill(0n));
  let poly = [1n];
  for (let i = 0; i < m; i++) {
    for (let j = 0; j < width; j++) mat[i]![j] = poly[j] ?? 0n;
    poly = MulMod(poly, M, F);
  }
  H.mat = mat;
  H.poly = poly;
}
/** NTL modular composition, accepting a polynomial or a prebuilt block argument.
 * @see Deviation: NTL word polynomial composition adapters
 */
export function CompMod(
  g: readonly bigint[],
  h: readonly bigint[] | zz_pXNewArgument,
  F: zz_pXModulus
): bigint[] {
  const k = F.arithmetic,
    G = k.norm(g);
  if (!(h instanceof zz_pXNewArgument)) {
    const m = sqrtCount(G.length);
    if (m === 0) return [];
    const H = new zz_pXNewArgument();
    build(H, h, F, m);
    return CompMod(G, H, F);
  }
  if (G.length <= 1) return G;
  const m = h.mat.length;
  if (m === 0) throw new Error('CompMod: uninitialized argument');
  const l = Math.floor((G.length + m - 1) / m),
    width = h.mat[0]!.length;
  const gmat = Array.from({ length: l }, (_, i) =>
    Array.from({ length: m }, (_, j) => G[i * m + j] ?? 0n)
  );
  const xmat = matrix_mul(gmat, h.mat, F.p, width, m);
  let t = k.norm(xmat[l - 1]!);
  if (l > 1) {
    const M = new zz_pXMultiplier(h.poly, F);
    for (let i = l - 2; i >= 0; i--) t = k.add(MulMod(t, M, F), k.norm(xmat[i]!));
  }
  return t;
}
/** NTL rebuild of a block-composition argument after shrinking the modulus.
 * @see Deviation: NTL word polynomial composition adapters
 */
export function reduce(H: zz_pXNewArgument, F: zz_pXModulus): void {
  const m = H.mat.length;
  if (m === 0) throw new Error('reduce: uninitialized argument');
  const h = m > 1 ? H.mat[1]! : H.poly;
  build(H, rem(h, F), F, m);
}
/** NTL two compositions sharing one power matrix.
 * @see Deviation: NTL word polynomial composition adapters
 */
export function Comp2Mod(
  g1: readonly bigint[],
  g2: readonly bigint[],
  h: readonly bigint[],
  F: zz_pXModulus
): [bigint[], bigint[]] {
  const G1 = F.arithmetic.norm(g1),
    G2 = F.arithmetic.norm(g2),
    m = sqrtCount(G1.length + G2.length);
  if (m === 0) return [[], []];
  const H = new zz_pXNewArgument();
  build(H, h, F, m);
  return [CompMod(G1, H, F), CompMod(G2, H, F)];
}
/** NTL three compositions sharing one power matrix.
 * @see Deviation: NTL word polynomial composition adapters
 */
export function Comp3Mod(
  g1: readonly bigint[],
  g2: readonly bigint[],
  g3: readonly bigint[],
  h: readonly bigint[],
  F: zz_pXModulus
): [bigint[], bigint[], bigint[]] {
  const G1 = F.arithmetic.norm(g1),
    G2 = F.arithmetic.norm(g2),
    G3 = F.arithmetic.norm(g3),
    m = sqrtCount(G1.length + G2.length + G3.length);
  if (m === 0) return [[], [], []];
  const H = new zz_pXNewArgument();
  build(H, h, F, m);
  return [CompMod(G1, H, F), CompMod(G2, H, F), CompMod(G3, H, F)];
}
/** Shared native half-GCD transformation used by GCD and sequence reconstruction. */
function wordHalfGCDKernel(p: bigint, options?: zz_pXOptions) {
  const F = new zz_pXModulus(null, p, options),
    k = F.arithmetic;
  const halfCut = F.modCrossover === 45 ? 90 : F.modCrossover === 90 ? 180 : 350;
  const gcdCut = F.modCrossover === 45 ? 400 : F.modCrossover === 90 ? 800 : 1400;
  type M = [bigint[], bigint[], bigint[], bigint[]];
  const identity = (): M => [[1n], [], [], [1n]];
  const fftSize = (n: number): void => {
    const bits = n <= 1 ? 0 : BigInt(n - 1).toString(2).length;
    if (bits > F.MaxRoot) throw new Error('Polynomial too big for FFT');
  };
  const divrem = (a: bigint[], b: bigint[]): [bigint[], bigint[]] => {
    const n = b.length - 1,
      m = a.length - 1;
    const cutoff = F.PrimeCnt <= 1 ? 180 : F.PrimeCnt === 2 ? 350 : 750;
    if (n > cutoff && m - n > cutoff) {
      // FFTDivRem first builds a modulus for large quotients; otherwise it
      // inverts the reversed divisor before either FFT-size check.
      if (m >= 3 * n) new zz_pXModulus(b, p, options);
      else {
        k.inverse(b[n]!);
        fftSize(2 * (m - n) + 1);
        fftSize(n);
      }
    }
    return k.divrem(a, b);
  };
  const apply = (m: M, a: bigint[], b: bigint[]): [bigint[], bigint[]] => {
    // Native matrix application uses a cyclic FFT even for small entries.
    fftSize(a.length - m[3].length - 1);
    return [k.add(k.mul(m[0], a), k.mul(m[1], b)), k.add(k.mul(m[2], a), k.mul(m[3], b))];
  };
  const matrixProduct = (a: M, b: M): M => {
    fftSize(a[3].length + b[3].length - 1);
    return [
      k.add(k.mul(a[0], b[0]), k.mul(a[1], b[2])),
      k.add(k.mul(a[0], b[1]), k.mul(a[1], b[3])),
      k.add(k.mul(a[2], b[0]), k.mul(a[3], b[2])),
      k.add(k.mul(a[2], b[1]), k.mul(a[3], b[3])),
    ];
  };
  const multiply = (a: bigint[], b: bigint[]) => mul(a, b, p, options);
  const step = (m: M, q: bigint[]): M => [
    m[2],
    m[3],
    k.sub(m[0], multiply(q, m[2])),
    k.sub(m[1], multiply(q, m[3])),
  ];
  const half = (a: bigint[], b: bigint[], reduction: number): M => {
    if (!b.length || b.length <= a.length - reduction) return identity();
    const shift = Math.max(a.length - 1 - 2 * reduction + 2, 0);
    let A = a.slice(shift),
      B = b.slice(shift);
    if (reduction <= halfCut) {
      let m = identity();
      const goal = A.length - 1 - reduction;
      while (B.length - 1 > goal) {
        const [q, r] = k.divrem(A, B);
        [A, B] = [B, r];
        m = step(m, q);
      }
      return m;
    }
    const d1 = Math.max(1, Math.min(Math.floor((reduction + 1) / 2), reduction - 1));
    let m = half(A, B, d1);
    [A, B] = apply(m, A, B);
    const d2 = B.length - a.length + shift + reduction;
    if (!B.length || d2 <= 0) return m;
    const [q, r] = divrem(A, B);
    [A, B] = [B, r];
    return matrixProduct(half(A, B, d2), step(m, q));
  };
  const extendedHalf = (a: bigint[], b: bigint[], reduction: number): [M, bigint[], bigint[]] => {
    if (!b.length || b.length <= a.length - reduction) return [identity(), a, b];
    if (reduction <= halfCut) {
      let m = identity();
      const goal = a.length - 1 - reduction;
      while (b.length - 1 > goal) {
        const [q, r] = k.divrem(a, b);
        [a, b] = [b, r];
        m = step(m, q);
      }
      return [m, a, b];
    }
    const degree = a.length - 1;
    const d1 = Math.max(1, Math.min(Math.floor((reduction + 1) / 2), reduction - 1));
    const m1 = half(a, b, d1);
    [a, b] = apply(m1, a, b);
    const d2 = b.length - 1 - degree + reduction;
    if (!b.length || d2 <= 0) return [m1, a, b];
    const [q, r] = divrem(a, b);
    const [m2, u, v] = extendedHalf(b, r, d2);
    return [matrixProduct(m2, step(m1, q)), u, v];
  };
  return { k, half, apply, gcdCut, divrem, extendedHalf, multiply };
}
/** NTL word-polynomial GCD, retaining its half-GCD and Euclidean crossovers.
 * @see Deviation: NTL word polynomial composition adapters
 * @see Deviation: NTL degree-dependent word contexts
 */
export function GCD(
  u: readonly bigint[],
  v: readonly bigint[],
  p: bigint,
  options?: zz_pXOptions
): bigint[] {
  const { k, half, apply, gcdCut, divrem } = wordHalfGCDKernel(p, options);
  const halve = (a: bigint[], b: bigint[]): [bigint[], bigint[]] => {
    const reduction = Math.floor(a.length / 2);
    if (!b.length || b.length <= a.length - reduction) return [a, b];
    const degree = a.length - 1,
      d1 = Math.max(1, Math.min(Math.floor((reduction + 1) / 2), reduction - 1));
    [a, b] = apply(half(a, b, d1), a, b);
    const d2 = b.length - 1 - degree + reduction;
    if (!b.length || d2 <= 0) return [a, b];
    const r = divrem(a, b)[1];
    [a, b] = [b, r];
    return apply(half(a, b, d2), a, b);
  };
  let a = k.norm(u),
    b = k.norm(v);
  if (a.length === b.length) {
    if (!a.length) return [];
    b = k.divrem(b, a)[1];
  } else if (a.length < b.length) [a, b] = [b, a];
  while (a.length - 1 > gcdCut && b.length) {
    [a, b] = halve(a, b);
    if (b.length) {
      const r = divrem(a, b)[1];
      [a, b] = [b, r];
    }
  }
  while (b.length) {
    const r = k.divrem(a, b)[1];
    [a, b] = [b, r];
  }
  if (!a.length || a[a.length - 1] === 1n) return a;
  const scale = k.inverse(a[a.length - 1]!);
  return k.norm(a.map((c) => c * scale));
}

/** Native Berlekamp–Massey sequence reconstruction; uses the first 2*m entries. */
function BerlekampMassey(
  a: readonly bigint[],
  m: number,
  p: bigint,
  options?: zz_pXOptions
): bigint[] {
  const k = new zz_pXModulus(null, p, options).arithmetic;
  let Lambda = [1n],
    Sigma: bigint[] = [],
    L = 0,
    Delta = 1n,
    shift = 0;
  for (let r = 1; r <= 2 * m; r++) {
    let discrepancy = 0n;
    for (let i = 0; i < Lambda.length; i++)
      discrepancy = k.mod(discrepancy + Lambda[i]! * a[r - i - 1]!);
    if (discrepancy === 0n) shift++;
    else if (2 * L < r) {
      const scale = k.mod(discrepancy * k.inverse(Delta)),
        temp = k.norm(Sigma.map((x) => x * scale));
      Sigma = Lambda;
      Lambda = k.sub(Lambda, temp.length ? [...Array<bigint>(shift + 1).fill(0n), ...temp] : []);
      shift = 0;
      L = r - L;
      Delta = discrepancy;
    } else {
      shift++;
      const scale = k.mod(discrepancy * k.inverse(Delta)),
        temp = k.norm(Sigma.map((x) => x * scale));
      Lambda = k.sub(Lambda, temp.length ? [...Array<bigint>(shift).fill(0n), ...temp] : []);
    }
  }
  return Array.from({ length: L + 1 }, (_, i) => Lambda[L - i] ?? 0n);
}
/** Native half-GCD sequence reconstruction above the Berlekamp–Massey crossover. */
function GCDMinPolySeq(
  a: readonly bigint[],
  m: number,
  p: bigint,
  options?: zz_pXOptions
): bigint[] {
  const { k, half } = wordHalfGCDKernel(p, options);
  const reversed = k.norm(a.slice(0, 2 * m).reverse());
  const power = [...Array<bigint>(2 * m).fill(0n), 1n];
  const denominator = half(power, reversed, m + 1)[3];
  const scale = k.inverse(denominator[denominator.length - 1]!);
  return k.norm(denominator.map((x) => x * scale));
}
/** NTL minimum polynomial of a linearly generated sequence; m bounds its degree.
 * @see Deviation: NTL word polynomial sequence adapters
 * @see Deviation: NTL degree-dependent word contexts
 */
export function MinPolySeq(
  a: readonly bigint[],
  m: number,
  p: bigint,
  options?: zz_pXOptions
): bigint[] {
  const F = new zz_pXModulus(null, p, options);
  if (m < 0 || m >= 2 ** 60) throw new Error('MinPoly: bad args');
  if (a.length < 2 * m) throw new Error('MinPoly: sequence too short');
  const cutoff =
    F.PrimeCnt === 0 ? 400 : F.modCrossover === 45 ? 480 : F.modCrossover === 90 ? 900 : 1600;
  return m > cutoff ? GCDMinPolySeq(a, m, p, options) : BerlekampMassey(a, m, p, options);
}

function PlainUpdateMap(a: bigint[], b: bigint[], F: zz_pXModulus): bigint[] {
  if (!b.length) return [];
  const k = F.arithmetic,
    n = F.n,
    m = n - b.length;
  const dot = (c: bigint[], offset = 0): bigint => {
    let sum = 0n;
    for (let j = 0; j < c.length && j + offset < a.length; j++)
      sum = k.mod(sum + c[j]! * a[j + offset]!);
    return sum;
  };
  const x = Array<bigint>(n).fill(0n);
  for (let i = 0; i <= m; i++) x[i] = dot(b, i);
  if (b.length !== 1) {
    let c = m < 0 ? b.slice(-m) : [...Array<bigint>(m).fill(0n), ...b];
    for (let i = m + 1; i < n; i++) {
      c = MulByXMod(c, F.f, F.p);
      x[i] = dot(c);
    }
  }
  return x;
}
/** NTL transposed modular multiplication by a cached multiplier.
 * @see Deviation: NTL word projection rebuild state
 * @see Deviation: NTL word polynomial projection adapters
 */
export function UpdateMap(a: readonly bigint[], B: zz_pXMultiplier, F: zz_pXModulus): bigint[] {
  const k = F.arithmetic,
    A = k.norm(a),
    n = F.n;
  if (A.length > n) throw new Error('UpdateMap: bad args');
  if (!B.UseFFT) {
    if (B.b.length - 1 > n)
      throw new RangeError('UpdateMap: multiplier degree exceeds native buffer bounds');
    return k.norm(PlainUpdateMap(A, B.b, F));
  }
  const fft = _zz_pX_transform_kernels(F, B);
  let R1 = fft.RevTofftRep(A, fft.k, 0);
  let R2 = fft.mul(R1, 'FRep');
  const V1 = fft.RevFromfftRep(R2, 0, n - 2).map((x) => k.mod(-x));
  R2 = fft.RevTofftRep(V1, fft.l, n - 1);
  R2 = fft.mul(R2, 'B1');
  R1 = fft.mul(R1, 'B2');
  R2 = fft.AddExpand(R2, R1);
  return k.norm(fft.RevFromfftRep(R2, 0, n - 1));
}
/** NTL transposed block composition: inner products against successive powers.
 * @see Deviation: NTL word polynomial projection adapters
 */
export function ProjectPowers(
  a: readonly bigint[],
  count: number,
  h: readonly bigint[] | zz_pXNewArgument,
  F: zz_pXModulus
): bigint[] {
  if (a.length > F.n || count < 0) throw new Error('ProjectPowers: bad args');
  if (!(h instanceof zz_pXNewArgument)) {
    if (count === 0) return [];
    const H = new zz_pXNewArgument();
    build(H, h, F, sqrtCount(count));
    return ProjectPowers(a, count, H, F);
  }
  if (count >= 2 ** 60) throw new Error('ProjectPowers: excessive args');
  const m = h.mat.length;
  if (m === 0) throw new Error('CompMod: uninitialized argument');
  // Native prepared count=0 indexes an empty matrix and crashes; reject this undefined input.
  if (count === 0) throw new Error('ProjectPowers: prepared argument requires a positive count');
  const width = h.mat[0]!.length,
    l = Math.ceil(count / m),
    k = F.arithmetic;
  const transpose = Array.from({ length: width }, (_, i) => h.mat.map((row) => row[i]!));
  let s = k.norm(a);
  const amat: bigint[][] = [Array.from({ length: width }, (_, i) => s[i] ?? 0n)];
  if (l > 1) {
    const M = new zz_pXMultiplier(h.poly, F);
    for (let i = 1; i < l; i++) {
      s = UpdateMap(s, M, F);
      amat.push(Array.from({ length: width }, (_, j) => s[j] ?? 0n));
    }
  }
  const xmat = matrix_mul(amat, transpose, F.p, m, width);
  return Array.from({ length: count }, (_, i) => xmat[Math.floor(i / m)]![i % m]!);
}

/** NTL projected sequence reconstruction for a prescribed linear functional.
 * @see Deviation: NTL word polynomial element minimum polynomials
 */
export function DoMinPolyMod(
  g: readonly bigint[],
  F: zz_pXModulus,
  m: number,
  R: readonly bigint[]
): bigint[] {
  return MinPolySeq(ProjectPowers(R, 2 * m, g, F), m, F.p, F.options);
}
function elementMinPolyContext(
  F: zz_pXModulus,
  boundOrStream: number | RandomStream,
  context?: RandomStream
): [number, RandomStream] {
  return typeof boundOrStream === 'number' ? [boundOrStream, context!] : [F.n, boundOrStream];
}
/** NTL random-projection minimum-polynomial divisor, with an explicit stream.
 * @see Deviation: NTL word polynomial element minimum polynomials
 */
export function ProbMinPolyMod(
  g: readonly bigint[],
  F: zz_pXModulus,
  stream: RandomStream
): bigint[];
export function ProbMinPolyMod(
  g: readonly bigint[],
  F: zz_pXModulus,
  m: number,
  stream: RandomStream
): bigint[];
export function ProbMinPolyMod(
  g: readonly bigint[],
  F: zz_pXModulus,
  boundOrStream: number | RandomStream,
  context?: RandomStream
): bigint[] {
  const [m, stream] = elementMinPolyContext(F, boundOrStream, context);
  if (m < 1 || m > F.n) throw new Error('ProbMinPoly: bad args');
  const R = Array.from({ length: F.n }, () => RandomBnd(F.p, stream, { word: true }));
  return DoMinPolyMod(g, F, m, R);
}
/** NTL certified minimum polynomial, preserving retries and random draws.
 * @see Deviation: NTL word polynomial element minimum polynomials
 */
export function MinPolyMod(g: readonly bigint[], F: zz_pXModulus, stream: RandomStream): bigint[];
export function MinPolyMod(
  g: readonly bigint[],
  F: zz_pXModulus,
  m: number,
  stream: RandomStream
): bigint[];
export function MinPolyMod(
  g: readonly bigint[],
  F: zz_pXModulus,
  boundOrStream: number | RandomStream,
  context?: RandomStream
): bigint[] {
  const [m, stream] = elementMinPolyContext(F, boundOrStream, context),
    k = F.arithmetic;
  if (m < 1 || m > F.n) throw new Error('MinPoly: bad args');
  let h = ProbMinPolyMod(g, F, m, stream);
  if (h.length - 1 === m) return h;
  let h1 = CompMod(h, g, F);
  if (!h1.length) return h;
  for (;;) {
    let R = Array.from({ length: F.n }, () => RandomBnd(F.p, stream, { word: true }));
    const H1 = new zz_pXMultiplier(h1, F);
    R = UpdateMap(R, H1, F);
    const h2 = DoMinPolyMod(g, F, m - (h.length - 1), R);
    h = k.mul(h, h2);
    if (h.length - 1 === m) return h;
    const h3 = CompMod(h2, g, F);
    h1 = MulMod(h3, H1, F);
    if (!h1.length) return h;
  }
}
/** NTL deterministic minimum polynomial when the element's minimum polynomial is irreducible.
 * @see Deviation: NTL word polynomial element minimum polynomials
 */
export function IrredPolyMod(g: readonly bigint[], F: zz_pXModulus, m = F.n): bigint[] {
  if (m < 1 || m > F.n) throw new Error('IrredPoly: bad args');
  return DoMinPolyMod(g, F, m, [1n]);
}

/** Native small-input product; the caller supplies at least one root. */
function IterBuild(a: bigint[], p: bigint): void {
  const mod = (x: bigint) => ((x % p) + p) % p,
    n = a.length;
  a[0] = mod(-a[0]!);
  for (let k = 1; k < n; k++) {
    const b = mod(-a[k]!);
    a[k] = mod(b + a[k - 1]!);
    for (let i = k - 1; i >= 1; i--) a[i] = mod(a[i]! * b + a[i - 1]!);
    a[0] = mod(a[0]! * b);
  }
}
/** NTL monic polynomial with prescribed roots, preserving multiplicities and zero roots.
 * @see Deviation: NTL cached word samplers and root products
 * @see Deviation: NTL degree-dependent word contexts
 */
export function BuildFromRoots(a: readonly bigint[], p: bigint, options?: zz_pXOptions): bigint[] {
  const F = new zz_pXModulus(null, p, options),
    k = F.arithmetic,
    n = a.length;
  if (n === 0) return [1n];
  const cutoff = F.modCrossover === 45 ? 150 : F.modCrossover === 90 ? 300 : 500;
  if (n <= cutoff) {
    const x = a.map(k.mod);
    IterBuild(x, p);
    x.push(1n);
    return x;
  }
  if (BigInt(n - 1).toString(2).length > F.MaxRoot) throw new Error('Polynomial too big for FFT');
  let m = 1;
  while (m < n) m *= 2;
  // The native tree pads with zero roots to the next power of two. Its FFT
  // stages store the monic coefficient implicitly; full coefficient arrays
  // give the same balanced products at the portable multiplication boundary.
  let tree = Array.from({ length: m }, (_, i) => [k.mod(-(a[i] ?? 0n)), 1n]);
  while (tree.length > 1) {
    const next: bigint[][] = [];
    for (let i = 0; i < tree.length; i += 2) next.push(k.mul(tree[i]!, tree[i + 1]!));
    tree = next;
  }
  return tree[0]!.slice(m - n);
}

/** Native word polynomial extended GCD, returning [d,s,t] with d=a*s+b*t.
 * @see Deviation: NTL word polynomial composition adapters
 * @see Deviation: NTL degree-dependent word contexts
 */
export function XGCD(
  a: readonly bigint[],
  b: readonly bigint[],
  p: bigint,
  options?: zz_pXOptions
): [bigint[], bigint[], bigint[]] {
  const { k, divrem, extendedHalf, multiply } = wordHalfGCDKernel(p, options);
  let U = k.norm(a),
    V = k.norm(b),
    Q: bigint[] = [],
    flag = 0;
  if (!U.length && !V.length) return [[], [1n], []];
  if (U.length === V.length) {
    let r: bigint[];
    [Q, r] = divrem(U, V);
    [U, V] = [V, r];
    flag = 1;
  } else if (U.length < V.length) {
    [U, V] = [V, U];
    flag = 2;
  }
  const [M, d] = extendedHalf(U, V, U.length);
  const s = flag === 0 ? M[0] : M[1];
  const t = flag === 0 ? M[1] : flag === 1 ? k.sub(M[0], multiply(Q, M[1])) : M[0];
  const w = k.inverse(d[d.length - 1]!);
  return [d, s, t].map((v) => k.norm(v.map((c) => c * w))) as [bigint[], bigint[], bigint[]];
}
