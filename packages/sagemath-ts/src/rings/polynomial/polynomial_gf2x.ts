/**
 * @module sage/rings/polynomial/polynomial_gf2x
 * @description Univariate Polynomials over GF(2) with bit-packed storage
 *
 * Port of: sage/rings/polynomial/polynomial_gf2x.pyx
 * Reference: NTL's GF2X library
 */

import { ValueError, ZeroDivisionError } from '../../errors.js';
import { current_randstate } from '../../misc/randstate.js';

export class GF2X {
  readonly bits: bigint;

  constructor(bits: bigint) {
    this.bits = bits < 0n ? -bits : bits;
  }

  static fromBigInt(n: bigint): GF2X {
    return new GF2X(n < 0n ? -n : n);
  }

  static fromCoeffs(coeffs: (0 | 1 | number | boolean)[]): GF2X {
    let bits = 0n;
    for (let i = 0; i < coeffs.length; i++) {
      const c = coeffs[i];
      const bit = typeof c === 'boolean' ? (c ? 1n : 0n) : BigInt(c) & 1n;
      if (bit === 1n) bits |= 1n << BigInt(i);
    }
    return new GF2X(bits);
  }

  static fromHex(hex: string): GF2X {
    const normalized = hex.toLowerCase().replace(/^0x/, '');
    return new GF2X(BigInt('0x' + normalized));
  }

  static zero(): GF2X {
    return new GF2X(0n);
  }
  static one(): GF2X {
    return new GF2X(1n);
  }
  static x(): GF2X {
    return new GF2X(2n);
  }

  static monomial(n: number): GF2X {
    if (n < 0) throw new ValueError('monomial degree must be non-negative');
    return new GF2X(1n << BigInt(n));
  }

  static random(n: number): GF2X {
    if (n <= 0) return GF2X.zero();
    const randstate = current_randstate();
    let bits = 0n;
    for (let i = 0; i < n; i++) {
      if (randstate.random() < 0.5) bits |= 1n << BigInt(i);
    }
    return new GF2X(bits);
  }

  degree(): number {
    if (this.bits === 0n) return -1;
    let d = 0;
    let b = this.bits;
    while (b > 1n) {
      b >>= 1n;
      d++;
    }
    return d;
  }

  numBits(): number {
    return this.degree() + 1;
  }

  getCoeff(i: number): 0 | 1 {
    if (i < 0) return 0;
    return ((this.bits >> BigInt(i)) & 1n) === 1n ? 1 : 0;
  }

  setCoeff(i: number, c: 0 | 1): GF2X {
    if (i < 0) throw new ValueError('coefficient index must be non-negative');
    const mask = 1n << BigInt(i);
    return c === 1 ? new GF2X(this.bits | mask) : new GF2X(this.bits & ~mask);
  }

  leadingCoefficient(): 0 | 1 {
    return this.bits === 0n ? 0 : 1;
  }
  constantTerm(): 0 | 1 {
    return (this.bits & 1n) === 1n ? 1 : 0;
  }

  weight(): number {
    let count = 0;
    let b = this.bits;
    while (b > 0n) {
      if ((b & 1n) === 1n) count++;
      b >>= 1n;
    }
    return count;
  }

  isZero(): boolean {
    return this.bits === 0n;
  }
  isOne(): boolean {
    return this.bits === 1n;
  }
  isX(): boolean {
    return this.bits === 2n;
  }
  isConstant(): boolean {
    return this.bits <= 1n;
  }

  eq(other: GF2X | bigint): boolean {
    const otherBits = other instanceof GF2X ? other.bits : other;
    return this.bits === otherBits;
  }

  add(other: GF2X): GF2X {
    return new GF2X(this.bits ^ other.bits);
  }
  sub(other: GF2X): GF2X {
    return this.add(other);
  }
  neg(): GF2X {
    return new GF2X(this.bits);
  }

  mul(other: GF2X): GF2X {
    if (this.isZero() || other.isZero()) return GF2X.zero();
    let multiplicand: bigint;
    let multiplier: bigint;
    if (this.degree() < other.degree()) {
      multiplicand = other.bits;
      multiplier = this.bits;
    } else {
      multiplicand = this.bits;
      multiplier = other.bits;
    }
    let result = 0n;
    let shift = 0n;
    while (multiplier > 0n) {
      if ((multiplier & 1n) === 1n) result ^= multiplicand << shift;
      multiplier >>= 1n;
      shift += 1n;
    }
    return new GF2X(result);
  }

  sqr(): GF2X {
    if (this.isZero()) return GF2X.zero();
    let result = 0n;
    let bits = this.bits;
    let pos = 0n;
    while (bits > 0n) {
      if ((bits & 1n) === 1n) result |= 1n << (pos * 2n);
      bits >>= 1n;
      pos += 1n;
    }
    return new GF2X(result);
  }

  mulByX(): GF2X {
    return new GF2X(this.bits << 1n);
  }
  leftShift(n: number): GF2X {
    return n < 0 ? this.rightShift(-n) : new GF2X(this.bits << BigInt(n));
  }
  rightShift(n: number): GF2X {
    return n < 0 ? this.leftShift(-n) : new GF2X(this.bits >> BigInt(n));
  }

  divRem(other: GF2X): [GF2X, GF2X] {
    if (other.isZero()) throw new ZeroDivisionError('polynomial division by zero');
    const otherDeg = other.degree();
    if (this.degree() < otherDeg) return [GF2X.zero(), new GF2X(this.bits)];
    let remainder = this.bits;
    let quotient = 0n;
    const divisorBits = other.bits;
    let remainderDeg = this.degree();
    while (remainderDeg >= otherDeg) {
      const shift = remainderDeg - otherDeg;
      quotient |= 1n << BigInt(shift);
      remainder ^= divisorBits << BigInt(shift);
      if (remainder === 0n) break;
      remainderDeg = new GF2X(remainder).degree();
    }
    return [new GF2X(quotient), new GF2X(remainder)];
  }

  div(other: GF2X): GF2X {
    return this.divRem(other)[0];
  }
  mod(other: GF2X): GF2X {
    return this.divRem(other)[1];
  }
  isDivisibleBy(other: GF2X): boolean {
    return !other.isZero() && this.mod(other).isZero();
  }

  pow(n: number | bigint): GF2X {
    let exp = typeof n === 'bigint' ? n : BigInt(n);
    if (exp < 0n) throw new ValueError('negative exponent not supported');
    if (exp === 0n) return GF2X.one();
    if (this.isZero()) return GF2X.zero();
    let result = GF2X.one();
    let base = new GF2X(this.bits);
    while (exp > 0n) {
      if ((exp & 1n) === 1n) result = result.mul(base);
      base = base.sqr();
      exp >>= 1n;
    }
    return result;
  }

  powMod(n: number | bigint, modulus: GF2X): GF2X {
    if (modulus.isZero()) throw new ZeroDivisionError('modulus cannot be zero');
    let exp = typeof n === 'bigint' ? n : BigInt(n);
    if (exp < 0n) return this.invMod(modulus).powMod(-exp, modulus);
    if (exp === 0n) return GF2X.one();
    let result = GF2X.one();
    let base = this.mod(modulus);
    while (exp > 0n) {
      if ((exp & 1n) === 1n) result = result.mul(base).mod(modulus);
      base = base.sqr().mod(modulus);
      exp >>= 1n;
    }
    return result;
  }

  trunc(n: number): GF2X {
    if (n <= 0) return GF2X.zero();
    return new GF2X(this.bits & ((1n << BigInt(n)) - 1n));
  }

  gcd(other: GF2X): GF2X {
    let a = new GF2X(this.bits);
    let b = new GF2X(other.bits);
    while (!b.isZero()) {
      const r = a.mod(b);
      a = b;
      b = r;
    }
    return a;
  }

  xgcd(other: GF2X): [GF2X, GF2X, GF2X] {
    let a = new GF2X(this.bits);
    let b = new GF2X(other.bits);
    let s0 = GF2X.one();
    let s1 = GF2X.zero();
    let t0 = GF2X.zero();
    let t1 = GF2X.one();
    while (!b.isZero()) {
      const [q, r] = a.divRem(b);
      a = b;
      b = r;
      const newS = s0.sub(q.mul(s1));
      s0 = s1;
      s1 = newS;
      const newT = t0.sub(q.mul(t1));
      t0 = t1;
      t1 = newT;
    }
    return [a, s0, t0];
  }

  invMod(m: GF2X): GF2X {
    if (m.isZero()) throw new ZeroDivisionError('modulus cannot be zero');
    const [g, s] = this.xgcd(m);
    if (!g.isOne()) throw new ValueError('polynomial is not invertible');
    return s.mod(m);
  }

  is_irreducible(): boolean {
    const deg = this.degree();
    if (deg <= 0) return false;
    if (deg === 1) return true;
    if (this.constantTerm() === 0) return false;
    // NTL IterIrredTest: compute x^{2^i} mod f for i = 1 to deg-1
    // At each step where i divides deg, check gcd(h - x, f) = 1
    // After deg iterations, h should equal x
    const x = GF2X.x();
    let h = x;
    for (let i = 1; i <= deg; i++) {
      h = h.sqr().mod(this);
      // For i < deg: if i divides deg, check that gcd(h - x, f) = 1
      if (i < deg && deg % i === 0) {
        const g = h.sub(x).gcd(this);
        if (!g.isOne()) return false;
      }
    }
    // After deg iterations, x^{2^deg} should equal x (mod f)
    return h.eq(x);
  }

  derivative(): GF2X {
    if (this.isConstant()) return GF2X.zero();
    let result = 0n;
    const deg = this.degree();
    for (let i = 1; i <= deg; i += 2) {
      if (this.getCoeff(i) === 1) result |= 1n << BigInt(i - 1);
    }
    return new GF2X(result);
  }

  reverse(hi?: number): GF2X {
    const n = hi ?? this.degree();
    if (n < 0) return GF2X.zero();
    let result = 0n;
    for (let i = 0; i <= n; i++) {
      if (this.getCoeff(i) === 1) result |= 1n << BigInt(n - i);
    }
    return new GF2X(result);
  }

  toCoeffs(): (0 | 1)[] {
    if (this.isZero()) return [0];
    const deg = this.degree();
    const result: (0 | 1)[] = [];
    for (let i = 0; i <= deg; i++) result.push(this.getCoeff(i));
    return result;
  }

  toBigInt(): bigint {
    return this.bits;
  }
  toHex(): string {
    return '0x' + this.bits.toString(16);
  }

  toString(variable: string = 'x'): string {
    if (this.isZero()) return '0';
    const deg = this.degree();
    const terms: string[] = [];
    for (let i = deg; i >= 0; i--) {
      if (this.getCoeff(i) === 1) {
        if (i === 0) terms.push('1');
        else if (i === 1) terms.push(variable);
        else terms.push(variable + '^' + i);
      }
    }
    return terms.join(' + ') || '0';
  }

  repr(): string {
    return this.toString();
  }
}

export class GF2XRing {
  private static instance: GF2XRing;
  readonly variableName: string = 'x';
  private constructor() {}
  static getInstance(): GF2XRing {
    if (!GF2XRing.instance) GF2XRing.instance = new GF2XRing();
    return GF2XRing.instance;
  }
  __call__(x: bigint | number | (0 | 1)[] | GF2X): GF2X {
    if (x instanceof GF2X) return x;
    if (Array.isArray(x)) return GF2X.fromCoeffs(x);
    return GF2X.fromBigInt(BigInt(x));
  }
  zero(): GF2X {
    return GF2X.zero();
  }
  one(): GF2X {
    return GF2X.one();
  }
  gen(): GF2X {
    return GF2X.x();
  }
  /**
   * Return a random polynomial of the given degree (bounds).
   *
   * As in Sage, `degree` is either an exact degree or a `(min, max)` pair,
   * and the result has degree **exactly** `degree` when a single integer is
   * given (`R.random_element(6).degree() == 6`).  The zero polynomial has
   * degree -1, so it is only ever returned when the minimum degree is -1.
   *
   * @see Reference: sage/rings/polynomial/polynomial_ring.py:1344 (random_element)
   */
  random_element(degree: number | [number, number] = [-1, 2], monic: boolean = false): GF2X {
    let lo: number;
    let hi: number;
    if (Array.isArray(degree)) {
      if (degree.length !== 2) {
        throw new ValueError(
          'degree argument must be an integer or a tuple of 2 integers (min_degree, max_degree)'
        );
      }
      [lo, hi] = degree;
      if (lo > hi) {
        throw new ValueError('minimum degree must be less or equal than maximum degree');
      }
      if (hi < -1) {
        throw new ValueError(`maximum degree (=${hi}) must be at least -1`);
      }
    } else {
      if (degree < -1) {
        throw new ValueError(`degree (=${degree}) must be at least -1`);
      }
      lo = degree;
      hi = degree;
    }

    if (lo <= -2) {
      lo = -1;
    }

    if (lo === -1 && hi === -1) {
      return GF2X.zero();
    }

    if (lo === -1 && monic) {
      if (hi === 0) return GF2X.one();
      lo = 0;
    }

    const randstate = current_randstate();
    const randomCoeff = (): bigint => (randstate.random() < 0.5 ? 0n : 1n);

    if (lo === -1) {
      let bits = 0n;
      for (let i = 0; i <= hi; i++) {
        if (randomCoeff() === 1n) bits |= 1n << BigInt(i);
      }
      return new GF2X(bits);
    }

    const coefs = new Array<bigint>(hi + 1).fill(0n);
    let nonzero = false;
    while (!nonzero) {
      for (let i = 0; i <= hi - lo; i++) {
        const c = randomCoeff();
        if (monic) {
          coefs[hi - i] = c;
          if (!nonzero && c !== 0n) {
            coefs[hi - i] = 1n;
            nonzero = true;
          }
        } else {
          coefs[hi - i] = c;
          nonzero = nonzero || c !== 0n;
        }
      }
    }
    for (let i = hi - lo + 1; i <= hi; i++) {
      coefs[hi - i] = randomCoeff();
    }

    let bits = 0n;
    for (let i = 0; i <= hi; i++) {
      if (coefs[i] === 1n) bits |= 1n << BigInt(i);
    }
    return new GF2X(bits);
  }
  is_field(): boolean {
    return false;
  }
  is_finite(): boolean {
    return false;
  }
  toString(): string {
    return 'Univariate Polynomial Ring in x over Finite Field of size 2';
  }
}

export const GF2X_Ring = GF2XRing.getInstance();

export function buildIrred(n: number): GF2X {
  if (n < 1) throw new ValueError('degree must be at least 1');
  // NTL returns x (not x+1) for n = 1: `SetX(f)`
  // (`GF2XFactoring.cpp:481-484`).
  if (n === 1) return new GF2X(2n);
  const highBit = 1n << BigInt(n);
  for (let low = 1n; low < highBit; low += 2n) {
    const candidate = new GF2X(highBit | low);
    if (candidate.is_irreducible()) return candidate;
  }
  throw new Error('Failed to find irreducible polynomial of degree ' + n);
}

export function buildSparseIrred(n: number): GF2X {
  if (n < 1) throw new ValueError('degree must be at least 1');
  // NTL: `SetX(f)` for n = 1 (`GF2XFactoring.cpp:912-915`).
  if (n === 1) return new GF2X(2n);
  const highBit = 1n << BigInt(n);
  for (let k = 1; k < n; k++) {
    const candidate = new GF2X(highBit | (1n << BigInt(k)) | 1n);
    if (candidate.is_irreducible()) return candidate;
  }
  for (let k3 = 3; k3 < n; k3++) {
    for (let k2 = 2; k2 < k3; k2++) {
      for (let k1 = 1; k1 < k2; k1++) {
        const candidate = new GF2X(
          highBit | (1n << BigInt(k3)) | (1n << BigInt(k2)) | (1n << BigInt(k1)) | 1n
        );
        if (candidate.is_irreducible()) return candidate;
      }
    }
  }
  return buildIrred(n);
}

export function buildRandomIrred(n: number): GF2X {
  if (n < 1) throw new ValueError('degree must be at least 1');
  // The only monic irreducible of degree 1 that NTL's BuildIrred returns is x.
  if (n === 1) return new GF2X(2n);
  const highBit = 1n << BigInt(n);
  const randstate = current_randstate();
  while (true) {
    let bits = highBit | 1n;
    for (let i = 1; i < n; i++) {
      if (randstate.random() < 0.5) bits |= 1n << BigInt(i);
    }
    const candidate = new GF2X(bits);
    if (candidate.is_irreducible()) return candidate;
  }
}

export function GF2X_BuildIrred_list(n: number): (0 | 1)[] {
  return buildIrred(n).toCoeffs();
}
export function GF2X_BuildSparseIrred_list(n: number): (0 | 1)[] {
  return buildSparseIrred(n).toCoeffs();
}
export function GF2X_BuildRandomIrred_list(n: number): (0 | 1)[] {
  return buildRandomIrred(n).toCoeffs();
}

export function squareFreeDecomp(f: GF2X): Array<[GF2X, number]> {
  if (f.isZero()) throw new ValueError('square-free decomposition requires non-zero polynomial');
  if (f.isConstant()) return [];

  // Helper to compute square root (only valid when f' = 0 in characteristic 2)
  function sqrt(poly: GF2X): GF2X {
    const coeffs: (0 | 1)[] = [];
    for (let i = 0; i <= poly.degree(); i += 2) {
      coeffs.push(poly.getCoeff(i));
    }
    return GF2X.fromCoeffs(coeffs);
  }

  // Standard square-free factorization for characteristic 2
  // Algorithm from Knuth TAOCP or any algebraic number theory text
  const result: Array<[GF2X, number]> = [];

  function decompose(g: GF2X, baseMultiplicity: number): void {
    if (g.isConstant()) return;

    const gp = g.derivative();

    if (gp.isZero()) {
      // g is a perfect square (in characteristic 2, if g' = 0 then g = h^2 for some h)
      decompose(sqrt(g), baseMultiplicity * 2);
      return;
    }

    // Standard Yun's algorithm variant for characteristic p
    let c = g.gcd(gp);
    let w = g.div(c);
    let i = 1;

    while (!w.isOne()) {
      const y = w.gcd(c);
      const z = w.div(y);
      if (!z.isOne()) {
        result.push([z, i * baseMultiplicity]);
      }
      w = y;
      c = c.div(y);
      i++;
    }

    // If c is not 1, it's a perfect square; recurse
    if (!c.isOne()) {
      decompose(sqrt(c), baseMultiplicity * 2);
    }
  }

  decompose(f, 1);

  // Merge factors with same base (shouldn't happen with correct algorithm, but just in case)
  const merged = new Map<bigint, number>();
  for (const [poly, mult] of result) {
    const existing = merged.get(poly.bits) || 0;
    merged.set(poly.bits, existing + mult);
  }
  const finalResult: Array<[GF2X, number]> = [];
  for (const [bits, mult] of merged) {
    finalResult.push([new GF2X(bits), mult]);
  }

  finalResult.sort((a, b) => a[1] - b[1]);
  return finalResult;
}

export function distinctDegreeFactorization(f: GF2X): Array<[GF2X, number]> {
  if (f.isZero()) throw new ValueError('DDF requires non-zero polynomial');
  if (f.degree() <= 0) return [];
  const result: Array<[GF2X, number]> = [];
  const x = GF2X.x();
  let h = x.mod(f);
  let degree = 1;
  let remaining = new GF2X(f.bits);
  while (remaining.degree() >= 2 * degree) {
    h = h.sqr().mod(remaining);
    const g = h.sub(x).gcd(remaining);
    if (!g.isOne()) {
      result.push([g, degree]);
      const [q] = remaining.divRem(g);
      remaining = q;
      h = h.mod(remaining);
    }
    degree++;
  }
  if (!remaining.isOne()) result.push([remaining, remaining.degree()]);
  return result;
}

export function equalDegreeFactorization(f: GF2X, d: number): GF2X[] {
  if (f.isZero()) throw new ValueError('EDF requires non-zero polynomial');
  const deg = f.degree();
  if (deg === 0) return [];
  if (deg === d) return [f];
  if (deg % d !== 0)
    throw new ValueError('polynomial degree ' + deg + ' is not a multiple of factor degree ' + d);
  const numFactors = deg / d;
  if (numFactors === 1) return [f];
  let attempts = 0;
  const maxAttempts = 100 * deg;
  const factors: GF2X[] = [f];
  while (factors.length < numFactors && attempts < maxAttempts) {
    attempts++;
    const a = GF2X.random(deg);
    if (a.isZero()) continue;
    let t = a.mod(f);
    let ap = a;
    for (let i = 1; i < d; i++) {
      ap = ap.sqr().mod(f);
      t = t.add(ap);
    }
    for (let i = factors.length - 1; i >= 0; i--) {
      const fi = factors[i]!;
      if (fi.degree() === d) continue;
      const tmod = t.mod(fi);
      const g = tmod.gcd(fi);
      if (!g.isOne() && !g.eq(fi)) {
        factors.splice(i, 1);
        factors.push(g);
        const [q] = fi.divRem(g);
        factors.push(q);
      }
    }
  }
  return factors;
}

export function factor(f: GF2X): Array<[GF2X, number]> {
  if (f.isZero()) throw new ValueError('cannot factor zero polynomial');
  if (f.isConstant()) return [];
  const result: Array<[GF2X, number]> = [];
  const sqfree = squareFreeDecomp(f);
  for (const [sqf, mult] of sqfree) {
    const ddf = distinctDegreeFactorization(sqf);
    for (const [prod, d] of ddf) {
      const irreducibles = equalDegreeFactorization(prod, d);
      for (const irr of irreducibles) result.push([irr, mult]);
    }
  }
  result.sort((a, b) => {
    const degDiff = a[0].degree() - b[0].degree();
    if (degDiff !== 0) return degDiff;
    if (a[0].bits < b[0].bits) return -1;
    if (a[0].bits > b[0].bits) return 1;
    return 0;
  });
  return result;
}
