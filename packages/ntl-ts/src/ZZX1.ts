/** Dense integer-polynomial kernels from NTL ZZX1.cpp.
 * @see Deviation: NTL stateful integer polynomial division
 * @see Deviation: NTL integer polynomial squarefree decomposition
 */
import { _ZZ_pX_euclidean_kernels } from './ZZ_pX1.js';
import type { PolynomialProductState } from './ZZ_pX.js';
import { UseFFTPrime, GetFFTPrime, FFTPrimeContext } from './FFT.js';
import { RandomStream } from './ZZ.js';
import { GCD as wordGCD } from './lzz_pX1.js';

const normalized = (a: readonly bigint[]): bigint[] => {
  const r = a.slice();
  while (r.length && r.at(-1) === 0n) r.pop();
  return r;
};
const gcd = (a: bigint, b: bigint): bigint => {
  while (b) [a, b] = [b, a % b];
  return a < 0n ? -a : a;
};
const content = (a: readonly bigint[]): bigint => {
  let d = 0n;
  for (const c of a) { d = gcd(d, c); if (d === 1n) break; }
  return (a.at(-1) ?? 0n) < 0n ? -d : d;
};
const primitive = (a: readonly bigint[]): bigint[] => {
  const f = normalized(a), c = content(f);
  return c ? f.map(x => x / c) : [];
};
const bits = (a: bigint): number => (a < 0n ? -a : a).toString(2).length;
const maxBits = (a: readonly bigint[]): number => a.reduce((m, x) => Math.max(m, x === 0n ? 0 : bits(x)), 0);
const dividesInteger = (a: bigint, b: bigint): boolean => b === 0n ? a === 0n : a % b === 0n;
const divideKernel = (
  a: readonly bigint[],
  b: readonly bigint[],
  state?: PolynomialProductState
): bigint[] | null => {
  const [status, q] = divide(a, b, { state });
  return status ? q : null;
};
const mod = (a: bigint, p: bigint): bigint => ((a % p) + p) % p;
const power = (a: bigint, e: bigint, p: bigint): bigint => {
  let r = 1n;
  while (e) { if (e & 1n) r = r * a % p; a = a * a % p; e >>= 1n; }
  return r;
};
// FFT.cpp candidate order with NTL_SP_NBITS=60. Only primality is needed here;
// the deterministic 64-bit witnesses replace randomized root-of-unity selection.
const primes: bigint[] = [];
let primeM = 59n, primeK = 0n;
function prime(index: number): bigint {
  while (primes.length <= index) {
    if (primeK === 0n) { primeM--; if (primeM < 5n) throw new RangeError('ran out of FFT primes'); primeK = 1n << (58n - primeM); }
    primeK--;
    const n = (1n << 59n) + (primeK << (primeM + 1n)) + (1n << primeM) + 1n;
    let d = n - 1n, s = 0;
    while (!(d & 1n)) { d >>= 1n; s++; }
    let probable = true;
    for (const a of [2n, 325n, 9375n, 28178n, 450775n, 9780504n, 1795265022n]) {
      let x = power(a % n, d, n);
      if (x === 1n || x === n - 1n) continue;
      let passed = false;
      for (let j = 1; j < s; j++) { x = x * x % n; if (x === n - 1n) { passed = true; break; } }
      if (!passed) { probable = false; break; }
    }
    if (probable) primes.push(n);
  }
  return primes[index]!;
}

/** Content removal, modular GCD, balanced CRT and exact-division certification.
 * @see Deviation: NTL integer polynomial squarefree decomposition
 */
export function GCD(a: readonly bigint[], b: readonly bigint[], state?: PolynomialProductState): bigint[] {
  const A = normalized(a), B = normalized(b);
  if (!A.length) return B.at(-1)! < 0n ? B.map(x => -x) : B;
  if (!B.length) return A.at(-1)! < 0n ? A.map(x => -x) : A;
  const ca = content(A), cb = content(B), c = gcd(ca, cb);
  const f1 = A.map(x => x / ca), f2 = B.map(x => x / cb), ld = gcd(f1.at(-1)!, f2.at(-1)!);
  const productState = state ?? {
    context: new FFTPrimeContext(), stream: new RandomStream(new Uint8Array(32)),
  };
  let product = 1n, g: bigint[] | null = null;
  for (let i = 0; ; i++) {
    // Native FFTInit precedes the leading-coefficient rejection and uses the
    // selected prime directly, without constructing an ordinary CRT context.
    UseFFTPrime(i, productState.context, productState.stream);
    const p = GetFFTPrime(i, productState.context);
    if (f1.at(-1)! % p === 0n || f2.at(-1)! % p === 0n) continue;
    const G = wordGCD(f1, f2, p, {fftPrime: i, state: productState}).map(x => x * (ld % p) % p);
    if (G.length === 1) return [c];
    if (g === null || G.length < g.length) {
      product = p; g = G.map(x => x > p / 2n ? x - p : x); continue;
    }
    if (G.length > g.length) continue;
    const modulus = product * p, inverse = power(product % p, p - 2n, p);
    let changed = false;
    g = g.map((x, j) => {
      let y = mod(x + product * mod((G[j]! - x) * inverse, p), modulus);
      if (y > modulus / 2n) y -= modulus;
      if (y !== x) changed = true;
      return y;
    });
    product = modulus;
    if (!changed) {
      const result = primitive(g);
      if (divideKernel(f1, result, productState) !== null && divideKernel(f2, result, productState) !== null) return result.map(x => x * c);
    }
  }
}
/** @internal Shared by the native squarefree decomposition port. */
export const _ZZX_kernels = { normalized, content, primitive, divide: divideKernel };

import { _ZZ_pX_power_kernels } from './ZZ_pX.js';
/** Native integer-polynomial transform selection from coefficient sizes and degrees.
 * @see Deviation: NTL stateful integer polynomial products
 */
export function ChooseSS(da: number, maxbitsa: number, db: number, maxbitsb: number): boolean {
  const ratio = SSRatio(da, maxbitsa, db, maxbitsb);
  const size = BigInt(maxbitsa) + BigInt(maxbitsb);
  // Native adds 64 before subtracting one: both intermediate sums must fit.
  if (size < -(1n << 63n) || size >= (1n << 63n) - 64n)
    throw new RangeError('ChooseSS coefficient sizes overflow');
  const k = (size + 63n) / 64n / 2n;
  return (
    (k >= 13n && ratio < 1.15) ||
    (k >= 26n && ratio < 1.3) ||
    (k >= 53n && ratio < 1.6) ||
    (k >= 106n && ratio < 1.8) ||
    (k >= 212n && ratio < 2.0)
  );
}
/** Exact dense integer product using the shared portable NTL product boundary.
 * @see Deviation: NTL multifactor Hensel lifting and integer products
 * @see Deviation: NTL stateful integer polynomial products
 */
export function mul(
  a: readonly bigint[],
  b: readonly bigint[],
  state?: PolynomialProductState
): bigint[] {
  const A = normalized(a),
    B = a === b ? A : normalized(b);
  if (!A.length || !B.length) return [];
  if (state) {
    const maxa = maxBits(A),
      maxb = maxBits(B),
      count = Math.min(A.length, B.length),
      minSize = Math.min(Math.ceil(maxa / 64), Math.ceil(maxb / 64));
    // Plain and Karatsuba branches never initialize native FFT-prime tables.
    if (count >= 80 && (minSize >= 30 || count >= 150)) {
      if (!ChooseSS(A.length - 1, maxa, B.length - 1, maxb)) {
        // NewFastCRTHelper adds two certification bits to HomMul/HomSqr's bound.
        const bound = BigInt(count).toString(2).length + maxa + maxb + 2;
        let product = 1n,
          index = 0;
        while (product.toString(2).length <= bound) {
          UseFFTPrime(index, state.context, state.stream);
          product *= GetFFTPrime(index++, state.context);
        }
      }
    }
  }
  const width =
    maxBits(A) + maxBits(B) + BigInt(Math.min(A.length, B.length)).toString(2).length + 1;
  const P = 1n << BigInt(width),
    half = P >> 1n,
    positive = (x: bigint[]) => x.map((c) => ((c % P) + P) % P);
  const X = positive(A),
    Y = A === B ? X : positive(B);
  return _ZZ_pX_power_kernels.multiply(X, Y, P).map((c) => (c > half ? c - P : c));
}
/** Exact dense integer square; preserves the shared alias-square path.
 * @see Deviation: NTL multifactor Hensel lifting and integer products
 * @see Deviation: NTL stateful integer polynomial products
 */
export function sqr(a: readonly bigint[], state?: PolynomialProductState): bigint[] {
  return mul(a, a, state);
}

/** NTL's integer-polynomial Schoenhage--Strassen padding ratio.
 * @see Deviation: NTL stateful modular polynomial products
 */
export function SSRatio(na: number, maxa: number, nb: number, maxb: number): number {
  const word = (x: number): bigint => {
    if (!Number.isInteger(x) || x < -(2 ** 63) || x >= 2 ** 63)
      throw new RangeError('SSRatio arguments must be signed native integers');
    return BigInt(x);
  };
  const a = word(na),
    ma = word(maxa),
    b = word(nb),
    mb = word(maxb);
  if (a <= 0n || b <= 0n) return 0;
  const n = a + b,
    bound = 2n + BigInt((a < b ? a : b).toString(2).length) + ma + mb;
  if (n >= (1n << 63n) - 1n || bound <= 0n || bound >= 1n << 63n)
    throw new RangeError('SSRatio dimensions or coefficient bound overflow');
  if (n >= 1n << 62n) throw new Error('NextPowerOfTwo: overflow');
  const l = BigInt(n.toString(2).length - 1);
  let mr = ((bound >> l) + 1n) << l;
  if (l >= 3n) {
    const alt = ((bound >> (l - 1n)) + 1n) << (l - 1n);
    if (alt < mr - mr / 8n) mr = alt;
  }
  if (mr >= (1n << 63n) - 1n) throw new RangeError('SSRatio padded coefficient bound overflow');
  return Number(mr + 1n) / Number(bound);
}

/** Native retained-output exact-division adapters.
 * @see Deviation: NTL stateful integer polynomial division
 */
export interface IntegerPolynomialDivisionOptions {
  previous?: readonly bigint[];
  state?: PolynomialProductState;
}
type DivisionResult = [0 | 1, bigint[]];
const rejectedDivision = (options: IntegerPolynomialDivisionOptions): DivisionResult => [
  0,
  normalized(options.previous ?? []),
];
const scalarDivision = (
  a: bigint[],
  b: bigint,
  options: IntegerPolynomialDivisionOptions
): DivisionResult => {
  if (!b) return a.length ? rejectedDivision(options) : [1, []];
  const q: bigint[] = [];
  for (const x of a) {
    if (x % b) return rejectedDivision(options);
    q.push(x / b);
  }
  return [1, q];
};
/** Native degree-first exact division with a coefficient bound.
 * @see Deviation: NTL stateful integer polynomial division
 */
export function PlainDivide(
  aa: readonly bigint[],
  bb: readonly bigint[],
  options: IntegerPolynomialDivisionOptions = {}
): DivisionResult {
  let a = normalized(aa),
    b = normalized(bb);
  if (!b.length) return a.length ? rejectedDivision(options) : [1, []];
  if (b.length === 1) return scalarDivision(a, b[0]!, options);
  // PlainDivide deliberately checks degree before its content reduction. In particular,
  // unlike HomDivide, a zero numerator and a nonconstant divisor report failure.
  if (a.length < b.length) return rejectedDivision(options);
  const ca = content(a),
    cb = content(b);
  if (ca % cb) return rejectedDivision(options);
  const cq = ca / cb;
  a = a.map((x) => x / ca);
  b = b.map((x) => x / cb);
  if (!dividesInteger(a.at(-1)!, b.at(-1)!) || !dividesInteger(a[0]!, b[0]!))
    return rejectedDivision(options);
  const bound = maxBits(a) + Math.floor((bits(BigInt(a.length)) + 1) / 2) + a.length - b.length;
  const q = Array<bigint>(a.length - b.length + 1).fill(0n),
    lead = b.at(-1)!;
  for (let i = q.length - 1; i >= 0; i--) {
    const x = a[i + b.length - 1]!;
    if (x % lead) return rejectedDivision(options);
    const t = x / lead;
    if (t !== 0n && bits(t) > bound) return rejectedDivision(options);
    q[i] = t;
    for (let j = b.length - 2; j >= 0; j--) a[i + j] = a[i + j]! - t * b[j]!;
  }
  if (a.slice(0, b.length - 1).some((x) => x !== 0n)) return rejectedDivision(options);
  return [1, normalized(q.map((x) => x * cq))];
}
/** Native modular quotient reconstruction and exact bit-bound certification.
 * @see Deviation: NTL stateful integer polynomial division
 */
export function HomDivide(
  aa: readonly bigint[],
  bb: readonly bigint[],
  options: IntegerPolynomialDivisionOptions = {}
): DivisionResult {
  let a = normalized(aa),
    b = normalized(bb);
  if (!b.length) return a.length ? rejectedDivision(options) : [1, []];
  if (!a.length) return [1, []];
  if (b.length === 1) return scalarDivision(a, b[0]!, options);
  if (a.length < b.length) return rejectedDivision(options);
  const ca = content(a),
    cb = content(b);
  if (ca % cb) return rejectedDivision(options);
  const cq = ca / cb;
  a = a.map((x) => x / ca);
  b = b.map((x) => x / cb);
  if (!dividesInteger(a.at(-1)!, b.at(-1)!) || !dividesInteger(a[0]!, b[0]!))
    return rejectedDivision(options);
  const aBound = maxBits(a),
    bBound = maxBits(b);
  let q: bigint[] = [],
    product = 1n,
    unstable = true;
  for (let i = 0; ; i++) {
    let p: bigint;
    if (options.state) {
      UseFFTPrime(i, options.state.context, options.state.stream);
      p = GetFFTPrime(i, options.state.context);
    } else p = prime(i);
    if (b.at(-1)! % p === 0n) continue;
    const k = _ZZ_pX_euclidean_kernels(p),
      A = k.norm(a),
      B = k.norm(b);
    if (!unstable) {
      const r = k.sub(A, k.mul(B, k.norm(q)));
      if (r.length >= B.length) unstable = true;
      else if (r.length) return rejectedDivision(options);
      else product *= p;
    }
    if (unstable) {
      const [Q, R] = k.divrem(A, B);
      if (R.length) return rejectedDivision(options);
      const modulus = product * p,
        inverse = k.inverse(product % p);
      unstable = false;
      q = Array.from({ length: Math.max(q.length, Q.length) }, (_, j) => {
        const x = q[j] ?? 0n;
        let y = mod(x + product * mod(((Q[j] ?? 0n) - x) * inverse, p), modulus);
        if (y > modulus / 2n) y -= modulus;
        if (y !== x) unstable = true;
        return y;
      });
      q = normalized(q);
      product = modulus;
    }
    const size = Math.min(b.length, q.length),
      sizeBits = size ? bits(BigInt(size)) : 0;
    if (!unstable && bits(product) > Math.max(aBound, bBound + maxBits(q) + sizeBits) + 3)
      return [1, q.map((x) => x * cq)];
  }
}
/** Native degree dispatch; scalar divisors preserve integer scalar division semantics.
 * @see Deviation: NTL stateful integer polynomial division
 */
export function divide(
  a: readonly bigint[],
  b: readonly bigint[] | bigint,
  options: IntegerPolynomialDivisionOptions = {}
): DivisionResult {
  const A = normalized(a);
  if (typeof b === 'bigint') return scalarDivision(A, b, options);
  const B = normalized(b);
  return B.length <= 9 || A.length - B.length <= 8
    ? PlainDivide(A, B, options)
    : HomDivide(A, B, options);
}
