import { FFTRoundUp } from './FFT_impl.js';
import { VectorRandomBnd, RandomStream } from './ZZ.js';
import {
  CalcMaxRoot,
  GetFFTPrime,
  UseFFTPrime,
  FFTPrimeContext,
  FFTFwd_trunc,
  FFTFwd_trans,
  FFTRev1_trans,
  type FFTPrimeInfo,
  FFTRev1_trunc,
} from './FFT.js';
import type { PolynomialProductState } from './ZZ_pX.js';
import { _ZZ_pX_euclidean_kernels } from './ZZ_pX1.js';
import { XGCD } from './lzz_pX1.js';
type K = ReturnType<typeof _ZZ_pX_euclidean_kernels>;
// lzz_pX.cpp ordinary multiplication and squaring crossovers by PrimeCnt.
const zz_pX_mul_crossover = [150, 150, 300, 500, 500] as const;
function wordCall<T>(run: () => T): T {
  try {
    return run();
  } catch (e) {
    if (e instanceof Error && e.message === "ZZ_pX InvMod: can't compute multiplicative inverse")
      throw new Error('InvMod: inverse undefined');
    if (e instanceof Error && e.message === 'ZZ_pX: division by zero')
      throw new Error('zz_pX: division by zero');
    throw e;
  }
}
function wordArithmetic(p: bigint): K {
  if (p <= 1n) throw new Error('zz_pContext: p must be > 1');
  if (p >= 1n << 60n) throw new Error('zz_pContext: modulus too big');
  const k = _ZZ_pX_euclidean_kernels(p);
  return {
    ...k,
    inverse: (a: bigint) => wordCall(() => k.inverse(a)),
    divrem: (a: bigint[], b: bigint[]) => wordCall(() => k.divrem(a, b)),
  };
}
/** Explicit native ordinary or FFT-prime coefficient context.
 * @see Deviation: NTL FFT-prime word coefficient contexts
 */
export interface zz_pXOptions {
  /** Native FFTInit context: p must equal this cached FFT prime. */
  fftPrime?: number;
  maxroot?: number;
  state?: PolynomialProductState;
}
// lzz_p.cpp chooses enough FFT primes that their product exceeds this bound.
// All four primes have CalcMaxRoot=25 on the bundled 60-bit word profile.
function modulusPrimeCount(p: bigint, maxroot: number, state?: PolynomialProductState): number {
  if (state) {
    const bound = (p * p) << BigInt(maxroot + 4);
    let product = 1n,
      count = 0;
    // Native initializes every required prime before reporting the four-prime limit.
    while (product <= bound) {
      UseFFTPrime(count, state.context, state.stream);
      product *= GetFFTPrime(count++, state.context);
    }
    if (count > 4) throw new Error('zz_pInit: too many primes');
    return count;
  }
  if (maxroot > 240) throw new Error('zz_pInit: too many primes');
  const bound = (p * p) << BigInt(maxroot + 4);
  let product = 1n;
  const primes = [
    882705526964617217n,
    891712726219358209n,
    855683929200394241n,
    1139410705724735489n,
  ];
  for (let i = 0; i < primes.length; i++) {
    product *= primes[i]!;
    if (product > bound) return i + 1;
  }
  throw new Error('zz_pInit: too many primes');
}
/** NTL word polynomial modulus and cached reciprocal; coefficients use bigint.
 * @see Deviation: NTL word modulus rebuild state
 * @see Deviation: NTL word polynomial quotient adapters
 * @see Deviation: NTL degree-dependent word contexts
 * @see Deviation: NTL FFT-prime word coefficient contexts
 * @see Deviation: NTL cold word-context initialization
 */
export class zz_pXModulus {
  f: bigint[] = [];
  n = -1;
  readonly p: bigint;
  /** @internal Requested context settings survive nested algorithm calls. */
  readonly options: Readonly<{
    maxroot: number;
    state?: PolynomialProductState;
    fftPrime?: number;
  }>;
  /** @internal Native context profile, distinct from requested maxroot. */
  readonly PrimeCnt: number;
  /** @internal Transform rows; unlike PrimeCnt, this is one for an FFT prime. */
  readonly NumPrimes: number;
  readonly MaxRoot: number;
  /** @internal Native construction threshold for the explicit word profile. */
  readonly modCrossover: number;
  /** @internal Portable arithmetic and reciprocal caches. */
  readonly arithmetic: K;
  /** @internal Reverse-series reciprocal for products of reduced operands. */
  reciprocal: bigint[] = [];
  constructor(f: readonly bigint[] | null, p: bigint, options?: zz_pXOptions) {
    // lzz_p.cpp INIT_FFT uses the selected prime directly and has no maxroot
    // argument. Capture the selected table for every nested polynomial call.
    if (options?.fftPrime !== undefined) {
      const index = options.fftPrime;
      if (!Number.isSafeInteger(index))
        throw new RangeError('FFT prime index must be a safe integer');
      if (index < 0) throw new Error('bad FFT prime index');
      const state = options.state ?? {
        context: new FFTPrimeContext(),
        stream: new RandomStream(new Uint8Array(32)),
      };
      UseFFTPrime(index, state.context, state.stream);
      const info = state.context.get(index);
      if (p !== info.q) throw new RangeError('coefficient modulus does not match FFT prime');
      this.p = p;
      this.arithmetic = wordArithmetic(p);
      this.PrimeCnt = 0;
      this.NumPrimes = 1;
      this.MaxRoot = CalcMaxRoot(info.q);
      this.modCrossover = 45;
      this.options = Object.freeze({ maxroot: this.MaxRoot, fftPrime: index, state });
      if (f !== null) build(this, f);
      return;
    }
    const maxroot = options?.maxroot ?? 25;
    if (!Number.isSafeInteger(maxroot))
      throw new RangeError('zz_pContext: maxroot must be a safe integer');
    if (maxroot < 0) throw new Error('zz_pContext: maxroot may not be negative');
    this.p = p;
    this.arithmetic = wordArithmetic(p);
    this.options = Object.freeze(options?.state ? { maxroot, state: options.state } : { maxroot });
    this.PrimeCnt = modulusPrimeCount(p, maxroot, options?.state);
    this.NumPrimes = this.PrimeCnt;
    this.MaxRoot = Math.min(
      maxroot,
      options?.state ? CalcMaxRoot(GetFFTPrime(this.PrimeCnt - 1, options.state.context)) : 25
    );
    this.modCrossover = this.PrimeCnt === 1 ? 45 : this.PrimeCnt === 2 ? 90 : 180;
    if (f !== null) build(this, f);
  }
}
// Native FFT representations survive failed builds independently of public f/n.
// Store their coefficient equivalents privately, including their transform sizes.
interface ModulusFFTCache {
  fallback?: PolynomialProductState;
  useFFT: boolean;
  k: number;
  l: number;
  fk: number;
  hk: number;
  f: bigint[];
  h: bigint[];
}
const modulusFFTCaches = new WeakMap<zz_pXModulus, ModulusFFTCache>();
function modulusFFTCache(F: zz_pXModulus): ModulusFFTCache {
  let cache = modulusFFTCaches.get(F);
  if (!cache) {
    cache = { useFFT: false, k: 0, l: 0, fk: -1, hk: -1, f: [], h: [] };
    modulusFFTCaches.set(F, cache);
  }
  return cache;
}
function fftSize(degree: number): number {
  return degree <= 1 ? 0 : BigInt(degree - 1).toString(2).length;
}
function checkFFTSize(size: number, F: zz_pXModulus): void {
  if (size > F.MaxRoot) throw new Error('Polynomial too big for FFT');
}
// Coefficient representation of a native cyclic transform, modulo X^(2^size)-1.
function foldFFT(a: readonly bigint[], size: number, arithmetic: K): bigint[] {
  const length = 2 ** size;
  if (a.length <= length) return arithmetic.norm(a);
  const result = Array<bigint>(length).fill(0n);
  for (let i = 0; i < a.length; i++) result[i % length] += a[i]!;
  return arithmetic.norm(result);
}
// In a failed rebuild an old reciprocal may exceed the declared product degree.
// Native inverse TFT then interpolates only its initialized frequency prefix.
function wordFFTPrimes(F: zz_pXModulus) {
  if (F.options.fftPrime !== undefined) return [F.options.state!.context.get(F.options.fftPrime)];
  const cache = modulusFFTCache(F);
  const state =
    F.options.state ??
    (cache.fallback ??= {
      context: new FFTPrimeContext(),
      stream: new RandomStream(new Uint8Array(32)),
    });
  return Array.from({ length: F.NumPrimes }, (_, i) => {
    UseFFTPrime(i, state.context, state.stream);
    return state.context.get(i);
  });
}
function truncatedProduct(
  a: bigint[],
  b: bigint[],
  size: number,
  requested: number,
  F: zz_pXModulus
): bigint[] {
  const arithmetic = F.arithmetic,
    len = FFTRoundUp(requested, size);
  // Use unreduced degrees: a leading coefficient may cancel modulo p while the
  // FFT-prime transform still sees a term outside the inverse TFT prefix.
  if (len === 2 ** size || a.length + b.length - 1 <= len)
    return foldFFT(arithmetic.mul(a, arithmetic.norm(b)), size, arithmetic);
  const primes = wordFFTPrimes(F);
  const residues = primes.map((info) => {
    const forward = (input: readonly bigint[], retained: boolean) => {
      // A retained HRep/B1 belongs to its original coefficient modulus. Native
      // multiplication reads its FFT-prime residues without re-normalizing it.
      const values = retained ? input.slice() : foldFFT(input, size, arithmetic),
        xn = FFTRoundUp(values.length, size);
      while (values.length < xn) values.push(0n);
      return FFTFwd_trunc(values, size, info, len, xn);
    };
    const A = forward(a, false),
      B = forward(b, true);
    return FFTRev1_trunc(
      A.map((x, i) => (x * B[i]!) % info.q),
      size,
      info,
      len
    );
  });
  return arithmetic.norm(FromModularRep(residues, F));
}
function fftRemainder(a: bigint[], F: zz_pXModulus): bigint[] {
  const cache = modulusFFTCache(F),
    arithmetic = F.arithmetic,
    n = F.n;
  checkFFTSize(cache.l, F);
  if (cache.l !== cache.hk) throw new Error('FFT rep mismatch');
  const highProduct = truncatedProduct(a.slice(n), cache.h, cache.l, 2 * n - 3, F);
  const q = highProduct.slice(n - 2, 2 * n - 3);
  checkFFTSize(cache.k, F);
  if (cache.k !== cache.fk) throw new Error('FFT rep mismatch');
  const product = foldFFT(arithmetic.mul(q, cache.f), cache.k, arithmetic);
  const input = foldFFT(a, cache.k, arithmetic);
  return arithmetic.norm(arithmetic.sub(input, product).slice(0, n));
}
/** NTL modulus or cached multiplier construction with native overwrite timing.
 * @see Deviation: NTL word modulus rebuild state
 * @see Deviation: NTL word polynomial quotient adapters
 * @see Deviation: NTL word polynomial projection adapters
 */
export function build(F: zz_pXModulus, f: readonly bigint[]): void;
export function build(B: zz_pXMultiplier, b: readonly bigint[], F: zz_pXModulus): void;
export function build(
  out: zz_pXModulus | zz_pXMultiplier,
  f: readonly bigint[],
  context?: zz_pXModulus
): void {
  if (out instanceof zz_pXMultiplier) return buildMultiplier(out, f, context!);
  const F = out;
  const k = F.arithmetic;
  F.f = k.norm(f);
  F.n = F.f.length - 1;
  const cache = modulusFFTCache(F);
  if (F.n <= 0) throw new Error('build: deg(f) must be at least 1');
  if (F.n <= F.modCrossover + 1) {
    cache.useFFT = false;
    return;
  }
  cache.useFFT = true;
  cache.k = fftSize(F.n);
  cache.l = fftSize(2 * F.n - 3);
  // FRep is overwritten before inversion; HRep remains intact if inversion fails.
  checkFFTSize(cache.k, F);
  cache.f = foldFFT(F.f, cache.k, k);
  cache.fk = cache.k;
  const reversed = F.f.slice().reverse();
  let a = [k.inverse(reversed[0]!)];
  while (a.length < F.n - 1) {
    const size = Math.min(F.n - 1, 2 * a.length);
    const correction = k.sub([2n], k.mul(reversed.slice(0, size), a).slice(0, size));
    a = k.mul(a, correction).slice(0, size);
    while (a.length < size) a.push(0n);
  }
  checkFFTSize(cache.l, F);
  cache.h = a.slice().reverse();
  cache.hk = cache.l;
  F.reciprocal = a;
}
/** NTL block reduction using the modulus' cached reversed reciprocal.
 * @see Deviation: NTL word modulus rebuild state
 * @see Deviation: NTL word polynomial quotient adapters
 */
export function rem(a: readonly bigint[], F: zz_pXModulus): bigint[] {
  if (F.n < 0) throw new Error('rem: uninitialized modulus');
  const k = F.arithmetic,
    A = k.norm(a),
    n = F.n,
    cache = modulusFFTCache(F);
  const reduce = (v: bigint[]): bigint[] => {
    if (v.length <= n) return v.slice();
    if (!cache.useFFT || v.length - 1 - n <= F.modCrossover) return k.divrem(v, F.f)[1];
    return fftRemainder(v, F);
  };
  if (A.length <= 2 * n - 1 || !cache.useFFT || A.length - 1 - n <= F.modCrossover)
    return reduce(A);
  // Native constructs a buffer of length 2*n-1 before entering this loop.
  if (n === 0) throw new Error('negative length in vector::SetLength');
  let buf: bigint[] = [],
    left = A.length;
  while (left > 0) {
    const amount = Math.min(2 * n - 1 - buf.length, left);
    buf = reduce(k.norm([...A.slice(left - amount, left), ...buf]));
    left -= amount;
  }
  return buf;
}
/** NTL reduced polynomial product, optionally using a cached multiplier.
 * @see Deviation: NTL word modulus rebuild state
 * @see Deviation: NTL word polynomial quotient adapters
 * @see Deviation: NTL word polynomial projection adapters
 */
export function MulMod(
  a: readonly bigint[],
  b: readonly bigint[] | zz_pXMultiplier,
  F: zz_pXModulus
): bigint[] {
  if (b instanceof zz_pXMultiplier) {
    const k = F.arithmetic,
      A = k.norm(a);
    if (A.length - 1 >= F.n)
      throw new Error(' bad args to MulMod(zz_pX,zz_pX,zz_pXMultiplier,zz_pXModulus)');
    if (!A.length) return [];
    const cache = modulusFFTCache(F);
    if (!b.UseFFT || !cache.useFFT || A.length - 1 <= F.modCrossover) {
      // A failed build can leave a plain multiplier flag and a new, oversized b.
      // Native ordinary mul still performs its FFT dispatch before reduction,
      // even when leading coefficients (or the entire product) cancel modulo p.
      const cutoff = zz_pX_mul_crossover[F.PrimeCnt]!;
      if (A.length - 1 > cutoff && b.b.length - 1 > cutoff)
        checkFFTSize(fftSize(A.length + b.b.length - 1), F);
      return rem(k.mul(A, b.b), F);
    }
    const multiplier = multiplierFFTCache(b);
    checkFFTSize(cache.l, F);
    if (cache.l !== multiplier.b1k) throw new Error('FFT rep mismatch');
    const high = truncatedProduct(
      A,
      multiplier.b1,
      cache.l,
      Math.max(2 ** cache.k, 2 * F.n - 2),
      F
    );
    const q = high.slice(F.n - 1, 2 * F.n - 2);
    if (cache.k !== multiplier.b2k) throw new Error('FFT rep mismatch');
    const product = foldFFT(k.mul(A, k.norm(multiplier.b2)), cache.k, k);
    checkFFTSize(cache.k, F);
    if (cache.k !== cache.fk) throw new Error('FFT rep mismatch');
    const correction = foldFFT(k.mul(q, cache.f), cache.k, k);
    return k.norm(k.sub(product, correction).slice(0, F.n));
  }
  if (F.n < 0) throw new Error('MulMod: uninitialized modulus');
  const k = F.arithmetic,
    A = k.norm(a),
    B = k.norm(b);
  if (A.length > F.n || B.length > F.n)
    throw new Error('bad args to MulMod(zz_pX,zz_pX,zz_pX,zz_pXModulus)');
  if (!A.length || !B.length) return [];
  const cache = modulusFFTCache(F);
  const cutoff = zz_pX_mul_crossover[F.PrimeCnt]!;
  if (cache.useFFT && A.length - 1 > cutoff && B.length - 1 > cutoff) {
    checkFFTSize(Math.max(fftSize(A.length + B.length - 1), cache.k), F);
    return fftRemainder(k.mul(A, B), F);
  }
  return rem(k.mul(A, B), F);
}
/** NTL reduced polynomial square with a prebuilt modulus.
 * @see Deviation: NTL word modulus rebuild state
 * @see Deviation: NTL word polynomial quotient adapters
 */
export function SqrMod(a: readonly bigint[], F: zz_pXModulus): bigint[] {
  if (F.n < 0) throw new Error('SqrMod: uninitialized modulus');
  const A = F.arithmetic.norm(a);
  if (A.length > F.n) throw new Error('bad args to SqrMod(zz_pX,zz_pX,zz_pXModulus)');
  const cache = modulusFFTCache(F);
  const cutoff = zz_pX_mul_crossover[F.PrimeCnt]!;
  if (cache.useFFT && A.length - 1 > cutoff) {
    checkFFTSize(Math.max(fftSize(2 * A.length - 1), cache.k), F);
    return fftRemainder(F.arithmetic.mul(A, A), F);
  }
  return rem(F.arithmetic.mul(A, A), F);
}
/** NTL multiply-by-X reduction using one scalar multiple of the modulus.
 * @see Deviation: NTL word polynomial quotient adapters
 */
export function MulByXMod(a: readonly bigint[], f: readonly bigint[], p: bigint): bigint[] {
  const k = wordArithmetic(p),
    A = k.norm(a),
    B = k.norm(f),
    n = B.length - 1,
    m = A.length - 1;
  if (m >= n || n === 0) throw new Error('MulByXMod: bad args');
  if (m < 0) return [];
  if (m < n - 1) return [0n, ...A];
  const z = k.mod(-A[n - 1]! * k.inverse(B[n]!));
  return k.norm(Array.from({ length: n }, (_, i) => (i ? A[i - 1]! : 0n) + z * B[i]!));
}
/** NTL inverse, or the monic gcd when no inverse exists.
 * @see Deviation: NTL word polynomial quotient adapters
 */
export function InvModStatus(
  a: readonly bigint[],
  f: readonly bigint[],
  p: bigint,
  options?: zz_pXOptions
): [number, bigint[]] {
  const k = new zz_pXModulus(null, p, options).arithmetic,
    A = k.norm(a),
    B = k.norm(f);
  if (A.length >= B.length || B.length === 1) throw new Error('InvModStatus: bad args');
  const [d, s] = wordCall(() => XGCD(A, B, p, options));
  return d.length === 1 && d[0] === 1n ? [0, s] : [1, d];
}
/** NTL polynomial inversion through the shared half-GCD implementation.
 * @see Deviation: NTL word polynomial quotient adapters
 */
export function InvMod(
  a: readonly bigint[],
  f: readonly bigint[],
  p: bigint,
  options?: zz_pXOptions
): bigint[] {
  const k = new zz_pXModulus(null, p, options).arithmetic,
    A = k.norm(a),
    B = k.norm(f);
  if (A.length >= B.length || B.length === 1) throw new Error('InvMod: bad args');
  const [d, s] = wordCall(() => XGCD(A, B, p, options));
  if (d.length !== 1 || d[0] !== 1n)
    throw new Error("zz_pX InvMod: can't compute multiplicative inverse");
  return s;
}
/** NTL left-to-right binary powering of X, inverting after negative powers.
 * @see Deviation: NTL word polynomial quotient adapters
 */
export function PowerXMod(e: bigint, F: zz_pXModulus): bigint[] {
  if (F.n < 0) throw new Error('PowerXMod: uninitialized modulus');
  if (e === 0n) return [1n];
  let h = [1n];
  for (const bit of (e < 0n ? -e : e).toString(2)) {
    h = SqrMod(h, F);
    if (bit === '1') h = MulByXMod(h, F.f, F.p);
  }
  return e < 0n ? InvMod(h, F.f, F.p, F.options) : h;
}
/** NTL binary powering of X+a, retaining the specialized linear multiply.
 * @see Deviation: NTL word polynomial quotient adapters
 */
export function PowerXPlusAMod(a: bigint, e: bigint, F: zz_pXModulus): bigint[] {
  if (F.n < 0) throw new Error('PowerXPlusAMod: uninitialized modulus');
  if (e === 0n) return [1n];
  let h = [1n];
  const k = F.arithmetic;
  for (const bit of (e < 0n ? -e : e).toString(2)) {
    h = SqrMod(h, F);
    if (bit === '1') h = k.add(MulByXMod(h, F.f, F.p), k.norm(h.map((x) => x * a)));
  }
  return e < 0n ? InvMod(h, F.f, F.p, F.options) : h;
}
/** NTL binary powering of a reduced polynomial, with native argument ordering.
 * @see Deviation: NTL word polynomial quotient adapters
 */
export function PowerMod(g: readonly bigint[], e: bigint, F: zz_pXModulus): bigint[] {
  const G = F.arithmetic.norm(g);
  if (G.length > F.n) throw new Error('PowerMod: bad args');
  if (e === 0n) return [1n];
  let h = [1n];
  for (const bit of (e < 0n ? -e : e).toString(2)) {
    h = SqrMod(h, F);
    if (bit === '1') h = MulMod(h, G, F);
  }
  return e < 0n ? InvMod(h, F.f, F.p, F.options) : h;
}

/** NTL cached polynomial multiplier; construct empty or from a reduced polynomial.
 * @see Deviation: NTL word polynomial projection adapters
 */
export class zz_pXMultiplier {
  b: bigint[] = [];
  UseFFT = 0;
  /** @internal Portable copy of the product represented by native FFT caches. */
  cached: bigint[] = [];
  /** @internal Reversed B2 polynomial and reversed B1 correction polynomial. */
  reverse: bigint[] = [];
  correction: bigint[] = [];
  constructor();
  constructor(b: readonly bigint[], F: zz_pXModulus);
  constructor(b?: readonly bigint[], F?: zz_pXModulus) {
    if (b !== undefined) build(this, b, F!);
  }
  /** Return an independent copy of the multiplier polynomial.
   * @see Deviation: NTL word polynomial projection adapters
   */
  val(): bigint[] {
    return this.b.slice();
  }
}
interface MultiplierFFTCache {
  primes: number;
  b1k: number;
  b2k: number;
  b1: bigint[];
  b2: bigint[];
}
const multiplierFFTCaches = new WeakMap<zz_pXMultiplier, MultiplierFFTCache>();
function multiplierFFTCache(B: zz_pXMultiplier): MultiplierFFTCache {
  let cache = multiplierFFTCaches.get(B);
  if (!cache) {
    cache = { primes: 0, b1k: -1, b2k: -1, b1: [], b2: [] };
    multiplierFFTCaches.set(B, cache);
  }
  return cache;
}
function buildMultiplier(B: zz_pXMultiplier, b: readonly bigint[], F: zz_pXModulus): void {
  if (F.n < 0) throw new Error('build zz_pXMultiplier: uninitialized modulus');
  B.b = F.arithmetic.norm(b);
  const degree = B.b.length - 1;
  if (degree >= F.n) throw new Error('build zz_pXMultiplier: deg(b) >= deg(f)');
  const cache = modulusFFTCache(F),
    multiplier = multiplierFFTCache(B),
    k = F.arithmetic;
  B.UseFFT = cache.useFFT && degree > F.modCrossover ? 1 : 0;
  if (B.UseFFT) {
    checkFFTSize(cache.l, F);
    // fftRep::DoSetSize retains its allocated prime count across plain builds.
    if (multiplier.primes !== 0 && multiplier.primes !== F.NumPrimes)
      throw new Error('fftRep: inconsistent use');
    multiplier.primes = F.NumPrimes;
    // Native writes B2 first; failure against HRep must retain the previous B1.
    multiplier.b2 = foldFFT(B.b, cache.k, k);
    multiplier.b2k = cache.k;
    B.cached = B.b.slice();
    B.reverse = Array.from({ length: F.n }, (_, i) => B.cached[F.n - 1 - i] ?? 0n);
    if (cache.l !== cache.hk) throw new Error('FFT rep mismatch');
    const product = truncatedProduct(B.b, cache.h, cache.l, 2 * F.n - 2, F);
    multiplier.b1 = product.slice(F.n - 1, 2 * F.n - 2);
    multiplier.b1k = cache.l;
    B.correction = Array.from({ length: F.n - 1 }, (_, i) => multiplier.b1[F.n - 2 - i] ?? 0n);
  }
}

/** NTL random word polynomial with degree below n and an explicit stream.
 * @see Deviation: NTL cached word samplers and root products
 */
export function random(n: number, p: bigint, stream: RandomStream): bigint[] {
  const F = new zz_pXModulus(null, p);
  if (n < 0) throw new Error('negative length in vector::SetLength');
  if (n >= 2 ** 57) throw new Error('excessive length in vector::SetLength');
  return F.arithmetic.norm(VectorRandomBnd(n, p, stream));
}

// The bundled AArch64 build contracts each CRT estimate update into FMADD.
// All operands here are nonnegative normals (or zero), with result <= 4.
// Align exact dyadic integers, round once with Number, then scale exactly.
const crtDoubleView = new DataView(new ArrayBuffer(8));
function fusedCRTEstimate(a: number, b: number, c: number): number {
  if (c === 0) return a * b;
  if (a === 0) return c;
  const parts = (x: number): [bigint, number] => {
    crtDoubleView.setFloat64(0, x, false);
    const bits = crtDoubleView.getBigUint64(0, false);
    return [(bits & ((1n << 52n) - 1n)) | (1n << 52n), Number((bits >> 52n) & 2047n) - 1075];
  };
  const [am, ae] = parts(a),
    [bm, be] = parts(b),
    [cm, ce] = parts(c);
  const exponent = Math.min(ae + be, ce);
  const exact = ((am * bm) << BigInt(ae + be - exponent)) + (cm << BigInt(ce - exponent));
  return Number(exact) * 2 ** exponent;
}
/** Native coefficient reconstruction from FFT-prime residue rows.
 * @see Deviation: NTL word CRT reconstruction
 */
export function FromModularRep(
  residues: readonly (readonly bigint[])[],
  F: zz_pXModulus
): bigint[] {
  if (residues.length !== F.NumPrimes)
    throw new RangeError('FromModularRep: incorrect number of prime rows');
  const len = residues[0]!.length;
  if (residues.some((row) => row.length !== len))
    throw new RangeError('FromModularRep: inconsistent coefficient counts');
  if (!len) return [];
  // FromfftRep bypasses CRT for p_info contexts and copies coefficient residues.
  if (F.PrimeCnt === 0) return residues[0]!.map((c) => F.arithmetic.mod(c));
  const arithmetic = F.arithmetic,
    primes = wordFFTPrimes(F);
  // The native vector routine uses exact centering for a single FFT prime.
  if (primes.length === 1) {
    const q = primes[0]!.q,
      half = q >> 1n;
    return residues[0]!.map((x) => {
      const value = ((x % q) + q) % q;
      return arithmetic.mod(value > half ? value - q : value);
    });
  }
  const modulus = primes.reduce((m, info) => m * info.q, 1n);
  const cofactors = primes.map((info) => modulus / info.q);
  const inverses = primes.map((info, i) =>
    _ZZ_pX_euclidean_kernels(info.q).inverse(cofactors[i]! % info.q)
  );
  // Match FromModularRep's nearest-integer CRT correction and binary64 sum.
  return Array.from({ length: len }, (_, j) => {
    let sum = 0n,
      approximate = 0;
    for (let i = 0; i < primes.length; i++) {
      const info = primes[i]!,
        value = (((residues[i]![j]! % info.q) + info.q) * inverses[i]!) % info.q;
      approximate = fusedCRTEstimate(Number(value), info.qrecip, approximate);
      sum += value * cofactors[i]!;
    }
    return arithmetic.mod(sum - BigInt(Math.floor(approximate + 0.5)) * modulus);
  });
}
interface WordFFTRep {
  k: number;
  rows: bigint[][];
}
// A build replaces the coefficient array only when the corresponding native
// representation is overwritten. Seeded prime identity retains its root table.
const forwardRepresentations = new WeakMap<
  bigint[],
  WeakMap<FFTPrimeInfo, Map<number, bigint[]>>
>();
/** @internal Native transposed FFT operations with retained quotient caches.
 * @see Deviation: NTL word projection rebuild state
 */
export function _zz_pX_transform_kernels(F: zz_pXModulus, B: zz_pXMultiplier) {
  const modulus = modulusFFTCache(F),
    multiplier = multiplierFFTCache(B);
  const primes = wordFFTPrimes(F),
    arithmetic = F.arithmetic;
  const RevTofftRep = (a: readonly bigint[], k: number, offset: number): WordFFTRep => {
    checkFFTSize(k, F);
    const n = 2 ** k,
      input = Array<bigint>(n).fill(0n);
    offset = ((offset % n) + n) % n;
    for (let i = 0; i < a.length; i++) {
      const j = (i + offset) % n;
      input[j] = arithmetic.mod(input[j]! + a[i]!);
    }
    return { k, rows: primes.map((info) => FFTRev1_trans(input, k, info)) };
  };
  const mul = (rep: WordFFTRep, name: 'FRep' | 'B1' | 'B2'): WordFFTRep => {
    const [k, coefficients] =
      name === 'FRep'
        ? ([modulus.fk, modulus.f] as const)
        : name === 'B1'
          ? ([multiplier.b1k, multiplier.b1] as const)
          : ([multiplier.b2k, multiplier.b2] as const);
    if (rep.k !== k) throw new Error('FFT rep mismatch');
    let cache = forwardRepresentations.get(coefficients);
    if (!cache) {
      cache = new WeakMap();
      forwardRepresentations.set(coefficients, cache);
    }
    return {
      k,
      rows: primes.map((info, i) => {
        let sizes = cache!.get(info);
        if (!sizes) {
          sizes = new Map();
          cache!.set(info, sizes);
        }
        let values = sizes.get(k);
        if (!values) {
          const input = foldFFT(coefficients, k, arithmetic),
            n = 2 ** k;
          while (input.length < n) input.push(0n);
          values = FFTFwd_trunc(input, k, info, n, n);
          sizes.set(k, values);
        }
        return rep.rows[i]!.map((x, j) => (x * values![j]!) % info.q);
      }),
    };
  };
  const RevFromfftRep = (rep: WordFFTRep, lo: number, hi: number): bigint[] => {
    const n = 2 ** rep.k;
    if (rep.rows.some((row) => row.length !== n)) throw new Error('RevFromfftRep: bad len');
    const count = Math.max(0, Math.min(hi, n - 1) - lo + 1);
    return FromModularRep(
      rep.rows.map((row, i) => FFTFwd_trans(row, rep.k, primes[i]!).slice(lo, lo + count)),
      F
    );
  };
  const AddExpand = (target: WordFFTRep, small: WordFFTRep): WordFFTRep => {
    if (target.k < small.k) throw new Error('AddExpand: bad args');
    const n = 2 ** small.k;
    if (target.rows.some((row) => row.length < n)) throw new Error('AddExpand: bad len');
    return {
      k: target.k,
      rows: target.rows.map((row, i) =>
        row.map((x, j) => (j < n ? (x + small.rows[i]![j]!) % primes[i]!.q : x))
      ),
    };
  };
  return { k: modulus.k, l: modulus.l, RevTofftRep, mul, RevFromfftRep, AddExpand };
}

/** Native ordinary word polynomial product.
 * @see Deviation: NTL word polynomial quotient adapters
 * @see Deviation: NTL degree-dependent word contexts
 */
export function mul(
  a: readonly bigint[],
  b: readonly bigint[],
  p: bigint,
  options?: zz_pXOptions
): bigint[] {
  if (a === b) return sqr(a, p, options);
  const F = new zz_pXModulus(null, p, options),
    k = F.arithmetic;
  const A = k.norm(a),
    B = k.norm(b),
    cutoff = zz_pX_mul_crossover[F.PrimeCnt]!;
  if (A.length - 1 > cutoff && B.length - 1 > cutoff)
    checkFFTSize(fftSize(A.length + B.length - 1), F);
  return k.mul(A, B);
}
/** Native ordinary word polynomial square.
 * @see Deviation: NTL word polynomial quotient adapters
 * @see Deviation: NTL degree-dependent word contexts
 */
export function sqr(a: readonly bigint[], p: bigint, options?: zz_pXOptions): bigint[] {
  const F = new zz_pXModulus(null, p, options),
    k = F.arithmetic,
    A = k.norm(a);
  if (A.length - 1 > zz_pX_mul_crossover[F.PrimeCnt]!) checkFFTSize(fftSize(2 * A.length - 1), F);
  return k.mul(A, A);
}
