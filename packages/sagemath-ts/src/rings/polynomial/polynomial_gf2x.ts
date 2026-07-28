/**
 * @module sage/rings/polynomial/polynomial_gf2x
 * @description Univariate Polynomials over GF(2) with bit-packed storage
 *
 * Port of: sage/rings/polynomial/polynomial_gf2x.pyx
 *
 * SageMath's `Polynomial_GF2X` is a thin wrapper around **NTL's `GF2X`**: the
 * whole arithmetic layer comes from `sage/libs/ntl/ntl_GF2X_linkage.pxi`,
 * `is_irreducible` is `GF2X_IterIrredTest`, and the three
 * `GF2X_Build*Irred_list` helpers are `GF2X_BuildIrred` /
 * `GF2X_BuildSparseIrred` / `GF2X_BuildRandomIrred`
 * (`polynomial_gf2x.pyx:262-336`).
 *
 * Following CLAUDE.md ("when SageMath delegates to an external library, we
 * MUST also delegate to our port of that library"), every operation below
 * that `@sagemath-ts/ntl-ts` provides is delegated to it rather than
 * reimplemented here.  This class stays as the Sage-level element type: it
 * keeps the bit-packed `bits` field, Sage's error types, and the handful of
 * conveniences NTL does not expose.  See {@link toNTL}/{@link fromNTL}.
 *
 * The operations still implemented locally, with the reason:
 *
 * - `GF2X.random` / {@link buildRandomIrred}: `ntl-ts`'s `GF2X.random` throws
 *   `NTL_NOT_IMPLEMENTED` and it has no `BuildRandomIrred` (which needs NTL's
 *   `IrredPolyMod`/`GF2XModulus`), so randomness is drawn from Sage's own
 *   `current_randstate()` here.
 * - `reverse(hi)`: `ntl-ts`'s `reverse()` has no `hi` argument (NTL's
 *   `reverse(c, a, hi)`, `ntl/include/NTL/GF2X.h`, does).
 * - {@link squareFreeDecomp}, {@link distinctDegreeFactorization},
 *   {@link equalDegreeFactorization}, {@link factor}: `ntl-ts`'s
 *   `SquareFreeDecomp`/`DistinctDegFactor`/`EqualDegFactor`/`factor` all throw
 *   `NTL_NOT_IMPLEMENTED`.
 *
 * Reference: reference/ntl/src/GF2X.cpp, GF2X1.cpp, GF2XFactoring.cpp
 */

import {
  GF2X as NTL_GF2X,
  GF2X_BuildIrred,
  GF2X_BuildSparseIrred,
  GF2X_GCD,
  GF2X_IterIrredTest,
  GF2X_PowerMod,
  GF2X_XGCD,
} from '@sagemath-ts/ntl-ts';
import { ValueError, ZeroDivisionError } from '../../errors.js';
import { current_randstate } from '../../misc/randstate.js';

/** Wrap the bit-packed representation as an NTL `GF2X`. */
function toNTL(f: GF2X): NTL_GF2X {
  return new NTL_GF2X(f.bits);
}

/** Unwrap an NTL `GF2X` back into this module's element type. */
function fromNTL(f: NTL_GF2X): GF2X {
  return new GF2X(f.rep());
}

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

  /** NTL `add` (`GF2X.cpp`). */
  add(other: GF2X): GF2X {
    return fromNTL(toNTL(this).add(toNTL(other)));
  }
  /** NTL `sub` — identical to `add` in characteristic 2. */
  sub(other: GF2X): GF2X {
    return fromNTL(toNTL(this).sub(toNTL(other)));
  }
  /** NTL `negate` — the identity in characteristic 2. */
  neg(): GF2X {
    return fromNTL(toNTL(this).negate());
  }

  /** NTL `mul` (carry-less multiplication). */
  mul(other: GF2X): GF2X {
    return fromNTL(toNTL(this).mul(toNTL(other)));
  }

  /** NTL `sqr`. */
  sqr(): GF2X {
    return fromNTL(toNTL(this).sqr());
  }

  mulByX(): GF2X {
    return this.leftShift(1);
  }
  /** NTL `LeftShift` (negative shifts go right, as in NTL). */
  leftShift(n: number): GF2X {
    return fromNTL(toNTL(this).LeftShift(n));
  }
  /** NTL `RightShift`. */
  rightShift(n: number): GF2X {
    return fromNTL(toNTL(this).RightShift(n));
  }

  /**
   * NTL `DivRem`.
   *
   * NTL reports division by zero as a bare `ArithmeticError`; Sage raises
   * `ZeroDivisionError`, so the zero divisor is rejected here first.
   */
  divRem(other: GF2X): [GF2X, GF2X] {
    if (other.isZero()) throw new ZeroDivisionError('polynomial division by zero');
    const [q, r] = toNTL(this).DivRem(toNTL(other));
    return [fromNTL(q), fromNTL(r)];
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

  /**
   * NTL `PowerMod` (`GF2X1.cpp:1743`).
   *
   * NTL requires `deg(a) < deg(f)`, so the base is reduced first; NTL reports
   * a zero modulus as an `ArithmeticError` from its division, while Sage
   * raises `ZeroDivisionError`.
   */
  powMod(n: number | bigint, modulus: GF2X): GF2X {
    if (modulus.isZero()) throw new ZeroDivisionError('modulus cannot be zero');
    const exp = typeof n === 'bigint' ? n : BigInt(n);
    if (exp < 0n) return this.invMod(modulus).powMod(-exp, modulus);
    if (exp === 0n) return GF2X.one();
    return fromNTL(GF2X_PowerMod(toNTL(this.mod(modulus)), exp, toNTL(modulus)));
  }

  /** NTL `trunc`. */
  trunc(n: number): GF2X {
    if (n <= 0) return GF2X.zero();
    return fromNTL(toNTL(this).trunc(n));
  }

  /** NTL `GCD` (`GF2X1.cpp`). */
  gcd(other: GF2X): GF2X {
    return fromNTL(GF2X_GCD(toNTL(this), toNTL(other)));
  }

  /** NTL `XGCD` (`GF2X1.cpp:3625`): returns `[d, s, t]` with `d = s*a + t*b`. */
  xgcd(other: GF2X): [GF2X, GF2X, GF2X] {
    const [d, s, t] = GF2X_XGCD(toNTL(this), toNTL(other));
    return [fromNTL(d), fromNTL(s), fromNTL(t)];
  }

  /**
   * Inverse modulo `m`.
   *
   * NTL's `InvMod` insists on `deg(a) < deg(f)` and raises
   * `InvMod: inverse undefined`; Sage's `Polynomial.inverse_mod` reduces first
   * and raises on a non-unit, so the reduction and the error type are applied
   * here around NTL's `XGCD`.
   */
  invMod(m: GF2X): GF2X {
    if (m.isZero()) throw new ZeroDivisionError('modulus cannot be zero');
    const [g, s] = this.xgcd(m);
    if (!g.isOne()) throw new ValueError('polynomial is not invertible');
    return s.mod(m);
  }

  /**
   * NTL `IterIrredTest` (`GF2XFactoring.cpp:8`), which is exactly what Sage's
   * `Polynomial_GF2X.is_irreducible` calls (`polynomial_gf2x.pyx:281`).
   */
  is_irreducible(): boolean {
    return GF2X_IterIrredTest(toNTL(this));
  }

  /** NTL `diff` (the formal derivative). */
  derivative(): GF2X {
    return fromNTL(toNTL(this).diff());
  }

  /**
   * `x^hi * f(1/x)`, i.e. NTL's `reverse(c, a, hi)`
   * (`ntl/include/NTL/GF2X.h`); `hi` defaults to `deg(f)`.
   *
   * `ntl-ts` only ports the no-argument `reverse()`, so the `hi` form is done
   * here.
   */
  reverse(hi?: number): GF2X {
    if (hi === undefined) return fromNTL(toNTL(this).reverse());
    if (hi < 0) return GF2X.zero();
    let result = 0n;
    for (let i = 0; i <= hi; i++) {
      if (this.getCoeff(i) === 1) result |= 1n << BigInt(hi - i);
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

/**
 * The lexicographically smallest irreducible polynomial of degree `n`.
 *
 * Delegates to `ntl-ts`'s `GF2X_BuildIrred` (NTL `BuildIrred`,
 * `GF2XFactoring.cpp:472`), which returns `x` — not `x + 1` — for n = 1
 * (`GF2XFactoring.cpp:481-484`: `if (n == 1) { SetX(f); return; }`).
 */
export function buildIrred(n: number): GF2X {
  if (n < 1) throw new ValueError('degree must be at least 1');
  return fromNTL(GF2X_BuildIrred(n));
}

/**
 * An irreducible polynomial of degree `n` of minimal weight.
 *
 * Delegates to `ntl-ts`'s `GF2X_BuildSparseIrred` (NTL `BuildSparseIrred`,
 * `GF2XFactoring.cpp:900`), i.e. NTL's precomputed minimal-weight
 * trinomial/pentanomial table for n <= 2048 and its `FindTrinom`/`FindPent`
 * searches above that.  The table is *not* reproduced by "search trinomials in
 * increasing k, then pentanomials in increasing (k3,k2,k1)", which is what
 * this function used to do.
 */
export function buildSparseIrred(n: number): GF2X {
  if (n < 1) throw new ValueError('degree must be at least 1');
  return fromNTL(GF2X_BuildSparseIrred(n));
}

/**
 * A random irreducible polynomial of degree `n`.
 *
 * SageMath calls NTL's `BuildRandomIrred(f, BuildSparseIrred(n))`
 * (`polynomial_gf2x.pyx:325-336`), which takes the minimal polynomial of a
 * random element of GF(2)[x]/(g) via `IrredPolyMod`.  `ntl-ts` ports neither
 * `BuildRandomIrred` nor `IrredPolyMod`/`GF2XModulus`, so this draws monic
 * candidates with a nonzero constant term from Sage's `current_randstate()`
 * and keeps the first irreducible one (NTL's own `IterIrredTest` decides).
 * Both procedures return a uniformly-distributed-enough random irreducible;
 * the concrete polynomial for a given seed differs from Sage's.
 */
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

/**
 * Sage's `GF2X_BuildIrred_list` (`polynomial_gf2x.pyx:262`): the coefficient
 * list, constant term first, padded to `n + 1` entries.
 */
export function GF2X_BuildIrred_list(n: number): (0 | 1)[] {
  return padCoeffs(buildIrred(n), n);
}
/** Sage's `GF2X_BuildSparseIrred_list` (`polynomial_gf2x.pyx:285`). */
export function GF2X_BuildSparseIrred_list(n: number): (0 | 1)[] {
  return padCoeffs(buildSparseIrred(n), n);
}
/** Sage's `GF2X_BuildRandomIrred_list` (`polynomial_gf2x.pyx:311`). */
export function GF2X_BuildRandomIrred_list(n: number): (0 | 1)[] {
  return padCoeffs(buildRandomIrred(n), n);
}

/** `[GF2(f[i]) for i in range(n + 1)]` — Sage's list shape. */
function padCoeffs(f: GF2X, n: number): (0 | 1)[] {
  const out: (0 | 1)[] = [];
  for (let i = 0; i <= n; i++) out.push(f.getCoeff(i));
  return out;
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
