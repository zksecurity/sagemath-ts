/**
 * NTL polynomials over GF(2) (GF2X).
 * @see Reference: ntl/src/GF2X.cpp, ntl/src/GF2X1.cpp, ntl/src/GF2XFactoring.cpp
 *
 * GF2X represents polynomials with coefficients in GF(2).
 * These are used extensively in cryptography and coding theory.
 */

import { GF2 } from './GF2.js';
import { GF2X_irred_tab, GF2X_IRRED_TAB_MAX } from './GF2X_irred_tab.js';

// ============================================
// Internal bit-packed helpers
// ============================================

/**
 * Number of bits in a nonnegative bigint (0 for zero).
 */
function bitLength(x: bigint): number {
  if (x === 0n) return 0;
  return x.toString(2).length;
}

/**
 * Carry-less (GF(2)) multiplication of two bit-packed polynomials.
 */
function clmul(a: bigint, b: bigint): bigint {
  if (a === 0n || b === 0n) return 0n;
  // Iterate over the operand with fewer bits.
  let x = a;
  let y = b;
  if (bitLength(y) > bitLength(x)) {
    const t = x;
    x = y;
    y = t;
  }
  let r = 0n;
  while (y !== 0n) {
    if (y & 1n) r ^= x;
    x <<= 1n;
    y >>= 1n;
  }
  return r;
}

/**
 * Squaring of a bit-packed polynomial over GF(2): interleave the bits with 0.
 */
function clsqr(a: bigint): bigint {
  let r = 0n;
  let x = a;
  let i = 0n;
  while (x !== 0n) {
    if (x & 1n) r |= 1n << (2n * i);
    x >>= 1n;
    i++;
  }
  return r;
}

/**
 * Division with remainder of bit-packed polynomials.
 * @throws ArithmeticError "GF2X: division by zero" when b is zero
 */
function cldivrem(a: bigint, b: bigint): [bigint, bigint] {
  if (b === 0n) throw new Error('GF2X: division by zero');
  const db = bitLength(b) - 1;
  let r = a;
  let q = 0n;
  let dr = bitLength(r) - 1;
  while (dr >= db) {
    const shift = BigInt(dr - db);
    q ^= 1n << shift;
    r ^= b << shift;
    dr = bitLength(r) - 1;
  }
  return [q, r];
}

/**
 * Polynomials over GF(2).
 * TypeScript port of NTL's GF2X class.
 *
 * Coefficients are bits: bit `i` of the internal representation is the
 * coefficient of `x^i`.
 */
export class GF2X {
  private _xrep: bigint; // Bit-packed representation, bit i = coefficient of x^i

  /**
   * Whether to use hex output format.
   */
  static HexOutput: boolean = false;

  /**
   * Creates a new GF2X polynomial.
   * @param coeffs - Bit coefficients (constant term first), or a bit-packed
   *                 bigint whose bit `i` is the coefficient of `x^i`
   */
  constructor(coeffs?: number[] | bigint);
  /**
   * Creates the monomial `c * x^i` (NTL's `GF2X(i, c)` constructor).
   * @param i - Exponent
   * @param c - Coefficient (0 or 1)
   */
  constructor(i: number, c: number | GF2);
  constructor(a?: number[] | bigint | number, c?: number | GF2) {
    if (a === undefined) {
      this._xrep = 0n;
      return;
    }
    if (typeof a === 'number') {
      // NTL GF2X.h: GF2X(i, c) == c*X^i
      if (a < 0) throw new Error('SetCoeff: negative index');
      const bit = c === undefined ? 1 : c instanceof GF2 ? c.rep() : Number(c) & 1;
      this._xrep = bit ? 1n << BigInt(a) : 0n;
      return;
    }
    if (typeof a === 'bigint') {
      if (a < 0n) throw new Error('GF2X: negative bit representation');
      this._xrep = a;
      return;
    }
    let rep = 0n;
    for (let i = 0; i < a.length; i++) {
      if (Number(a[i]) & 1) rep |= 1n << BigInt(i);
    }
    this._xrep = rep;
  }

  /**
   * Wraps a bit-packed representation without copying (internal).
   */
  private static _fromRep(rep: bigint): GF2X {
    const f = new GF2X();
    f._xrep = rep;
    return f;
  }

  /**
   * The bit-packed representation: bit `i` is the coefficient of `x^i`.
   * @returns The bit-packed value
   */
  rep(): bigint {
    return this._xrep;
  }

  // ============================================
  // Basic Properties
  // ============================================

  /**
   * Returns the zero polynomial.
   */
  static zero(): GF2X {
    return new GF2X();
  }

  /**
   * Returns the one polynomial.
   */
  static one(): GF2X {
    return GF2X._fromRep(1n);
  }

  /**
   * Returns the polynomial X.
   */
  static X(): GF2X {
    return GF2X._fromRep(2n);
  }

  /**
   * Returns the degree (-1 for zero polynomial).
   * @returns The degree
   */
  deg(): number {
    return bitLength(this._xrep) - 1;
  }

  /**
   * Normalizes the polynomial (removes leading zeros).
   *
   * The bit-packed representation is always normalized, so this is a no-op;
   * it is kept for API compatibility with NTL.
   */
  normalize(): void {
    // no-op: bigint representation carries no leading zero words
  }

  /**
   * Sets maximum length (preallocates storage).
   *
   * Storage is managed by the bigint representation, so this only validates
   * its argument the way NTL does.
   * @param n - Maximum degree + 1
   * @throws LogicError "GF2X::SetMaxLength: negative length"
   */
  SetMaxLength(n: number): void {
    if (n < 0) throw new Error('GF2X::SetMaxLength: negative length');
  }

  /**
   * Sets length (number of coefficients); truncates when shortening.
   * @param n - Length
   * @throws LogicError "SetLength: negative index"
   */
  SetLength(n: number): void {
    if (n < 0) throw new Error('SetLength: negative index');
    this._xrep &= (1n << BigInt(n)) - 1n;
  }

  /**
   * Gets coefficient at index i.
   * @param i - Index
   * @returns The coefficient (0 or 1); 0 for i < 0 or i > deg
   */
  coeff(i: number): GF2 {
    if (i < 0) return GF2.zero();
    return new GF2(Number((this._xrep >> BigInt(i)) & 1n));
  }

  /**
   * Sets coefficient at index i to 1.
   * @param i - Index
   */
  SetCoeff(i: number): void;
  /**
   * Sets coefficient at index i.
   * @param i - Index
   * @param a - Value (0 or 1)
   */
  SetCoeff(i: number, a: number | GF2): void;
  SetCoeff(i: number, a?: number | GF2): void {
    if (i < 0) throw new Error('SetCoeff: negative index');
    const bit = a === undefined ? 1 : a instanceof GF2 ? a.rep() : Number(a) & 1;
    const mask = 1n << BigInt(i);
    if (bit) this._xrep |= mask;
    else this._xrep &= ~mask;
  }

  /**
   * Returns the leading coefficient.
   * @returns Leading coefficient (always 1 for nonzero, 0 for the zero polynomial)
   */
  LeadCoeff(): GF2 {
    return new GF2(this._xrep === 0n ? 0 : 1);
  }

  /**
   * Returns the constant term.
   * @returns Constant term
   */
  ConstTerm(): GF2 {
    return new GF2(Number(this._xrep & 1n));
  }

  // ============================================
  // Arithmetic Operations
  // ============================================

  /**
   * Adds two polynomials (XOR).
   * @param other - Polynomial to add
   * @returns The sum
   */
  add(other: GF2X): GF2X {
    return GF2X._fromRep(this._xrep ^ other._xrep);
  }

  /**
   * Subtracts a polynomial (same as add).
   * @param other - Polynomial to subtract
   * @returns The difference
   */
  sub(other: GF2X): GF2X {
    return GF2X._fromRep(this._xrep ^ other._xrep);
  }

  /**
   * Multiplies two polynomials.
   * @param other - Polynomial to multiply
   * @returns The product
   */
  mul(other: GF2X): GF2X {
    return GF2X._fromRep(clmul(this._xrep, other._xrep));
  }

  /**
   * Squares the polynomial.
   * @returns The square
   */
  sqr(): GF2X {
    return GF2X._fromRep(clsqr(this._xrep));
  }

  /**
   * Negates the polynomial (identity in characteristic 2).
   * @returns The negation
   */
  negate(): GF2X {
    return GF2X._fromRep(this._xrep);
  }

  // ============================================
  // Division Operations
  // ============================================

  /**
   * Polynomial division with remainder.
   * @param b - Divisor
   * @returns [quotient, remainder]
   * @throws ArithmeticError "GF2X: division by zero"
   */
  DivRem(b: GF2X): [GF2X, GF2X] {
    const [q, r] = cldivrem(this._xrep, b._xrep);
    return [GF2X._fromRep(q), GF2X._fromRep(r)];
  }

  /**
   * Polynomial division (quotient only).
   * @param b - Divisor
   * @returns The quotient
   * @throws ArithmeticError "GF2X: division by zero"
   */
  div(b: GF2X): GF2X {
    return GF2X._fromRep(cldivrem(this._xrep, b._xrep)[0]);
  }

  /**
   * Polynomial remainder.
   * @param b - Divisor
   * @returns The remainder
   * @throws ArithmeticError "GF2X: division by zero"
   */
  rem(b: GF2X): GF2X {
    return GF2X._fromRep(cldivrem(this._xrep, b._xrep)[1]);
  }

  // ============================================
  // GCD and Related
  // ============================================

  /**
   * Computes GCD of two polynomials.
   * @param a - First polynomial
   * @param b - Second polynomial
   * @returns The GCD (0 when both arguments are 0)
   */
  static GCD(a: GF2X, b: GF2X): GF2X {
    let u = a._xrep;
    let v = b._xrep;
    while (v !== 0n) {
      const r = cldivrem(u, v)[1];
      u = v;
      v = r;
    }
    return GF2X._fromRep(u);
  }

  /**
   * Extended GCD.
   *
   * Ports NTL's `XGCD` (ntl/src/GF2X1.cpp:3625).  NTL uses a half-GCD
   * recursion for large inputs; the (unique) `d`, `s`, `t` it returns are the
   * ones produced by the plain extended Euclidean algorithm used here, with
   * `deg(s) < deg(b) - deg(d)` and `deg(t) < deg(a) - deg(d)`.
   * @param a - First polynomial
   * @param b - Second polynomial
   * @returns [d, s, t] where d = s*a + t*b
   */
  static XGCD(a: GF2X, b: GF2X): [GF2X, GF2X, GF2X] {
    // NTL BaseXGCD: XGCD(a, 0) = (a, 1, 0)
    let r0 = a._xrep;
    let r1 = b._xrep;
    let s0 = 1n;
    let s1 = 0n;
    let t0 = 0n;
    let t1 = 1n;
    while (r1 !== 0n) {
      const [q, r] = cldivrem(r0, r1);
      const s = s0 ^ clmul(q, s1);
      const t = t0 ^ clmul(q, t1);
      r0 = r1;
      r1 = r;
      s0 = s1;
      s1 = s;
      t0 = t1;
      t1 = t;
    }
    return [GF2X._fromRep(r0), GF2X._fromRep(s0), GF2X._fromRep(t0)];
  }

  // ============================================
  // Modular Arithmetic
  // ============================================

  /**
   * Computes (a * b) mod f.
   * @param a - First polynomial
   * @param b - Second polynomial
   * @param f - Modulus polynomial
   * @returns (a * b) mod f
   */
  static MulMod(a: GF2X, b: GF2X, f: GF2X): GF2X {
    return GF2X._fromRep(cldivrem(clmul(a._xrep, b._xrep), f._xrep)[1]);
  }

  /**
   * Computes a^2 mod f.
   * @param a - Polynomial
   * @param f - Modulus polynomial
   * @returns a^2 mod f
   */
  static SqrMod(a: GF2X, f: GF2X): GF2X {
    return GF2X._fromRep(cldivrem(clsqr(a._xrep), f._xrep)[1]);
  }

  /**
   * Computes modular inverse.
   * @param a - Polynomial to invert
   * @param f - Modulus polynomial
   * @returns a^(-1) mod f
   * @throws LogicError "InvMod: bad args" if deg(a) >= deg(f) or deg(f) == 0,
   *         InvModError "InvMod: inverse undefined" if gcd(a, f) != 1
   */
  static InvMod(a: GF2X, f: GF2X): GF2X {
    // NTL GF2X1.cpp:2368 BaseInvMod / :2389 InvMod
    if (a.deg() >= f.deg() || f.deg() === 0) throw new Error('InvMod: bad args');
    const [d, s] = GF2X.XGCD(a, f);
    if (!d.IsOne()) throw new Error('InvMod: inverse undefined');
    return s;
  }

  /**
   * Computes a^e mod f.
   * @param a - Base polynomial
   * @param e - Exponent (may be negative, which inverts a)
   * @param f - Modulus polynomial
   * @returns a^e mod f
   * @throws LogicError "PowerMod: bad args" if deg(a) >= deg(f)
   */
  static PowerMod(a: GF2X, e: bigint, f: GF2X): GF2X {
    // NTL GF2X1.cpp:1743 -- sliding-window powering; identical result to
    // square-and-multiply.
    if (a.deg() >= f.deg()) throw new Error('PowerMod: bad args');
    if (e === 0n) return GF2X.one();
    let base = a;
    let n = e;
    if (n < 0n) {
      base = GF2X.InvMod(a, f);
      n = -n;
    }
    let res = GF2X.one().rem(f);
    while (n > 0n) {
      if (n & 1n) res = GF2X.MulMod(res, base, f);
      n >>= 1n;
      if (n > 0n) base = GF2X.SqrMod(base, f);
    }
    return res;
  }

  // ============================================
  // Factoring
  // ============================================

  /**
   * Factors a polynomial.
   * @returns Array of [factor, multiplicity] pairs
   */
  factor(): Array<[GF2X, number]> {
    throw new Error('NTL_NOT_IMPLEMENTED: GF2X.factor');
  }

  /**
   * Computes squarefree factorization.
   * @returns Array of [factor, multiplicity] pairs
   */
  SquareFreeDecomp(): Array<[GF2X, number]> {
    throw new Error('NTL_NOT_IMPLEMENTED: GF2X.SquareFreeDecomp');
  }

  /**
   * Distinct degree factorization.
   * @returns Array of [product of degree-d irreducibles, d] pairs
   */
  DistinctDegFactor(): Array<[GF2X, number]> {
    throw new Error('NTL_NOT_IMPLEMENTED: GF2X.DistinctDegFactor');
  }

  /**
   * Equal degree factorization.
   * @param d - Degree of factors
   * @returns Array of irreducible factors
   */
  EqualDegFactor(d: number): GF2X[] {
    throw new Error('NTL_NOT_IMPLEMENTED: GF2X.EqualDegFactor');
  }

  /**
   * Berlekamp factorization.
   * @returns Array of irreducible factors
   */
  BerlekampFactor(): GF2X[] {
    throw new Error('NTL_NOT_IMPLEMENTED: GF2X.BerlekampFactor');
  }

  // ============================================
  // Irreducibility
  // ============================================

  /**
   * Tests if polynomial is irreducible.
   *
   * Ports NTL's `IterIrredTest` (ntl/src/GF2XFactoring.cpp:8): it computes
   * `x^(2^d) - x mod f` iteratively for `d = 1 .. deg(f)/2`, batching the GCD
   * checks, and is deterministic.
   * @returns True if irreducible (false for constants and the zero polynomial)
   */
  isIrreducible(): boolean {
    const f = this;
    const df = f.deg();

    if (df <= 0) return false;
    if (df === 1) return true;

    // h = X^2 mod f
    let h = GF2X.SqrMod(GF2X.X(), f);

    const X = GF2X.X();
    let i = 0;
    let g = h;
    let d = 1;
    let limit = 2;
    let limit_sqr = limit * limit;
    let prod = GF2X.one();

    while (2 * d <= df) {
      const t = g.add(X);
      prod = GF2X.MulMod(prod, t, f);
      i++;
      if (i === limit_sqr) {
        if (!GF2X.GCD(f, prod).IsOne()) return false;
        prod = GF2X.one();
        limit++;
        limit_sqr = limit * limit;
        i = 0;
      }

      d = d + 1;
      if (2 * d <= df) {
        g = GF2X.SqrMod(g, f);
      }
    }

    if (i > 0) {
      if (!GF2X.GCD(f, prod).IsOne()) return false;
    }

    return true;
  }

  /**
   * Builds the lexicographically smallest irreducible polynomial of degree n.
   *
   * Deterministic: ports NTL's `BuildIrred` (ntl/src/GF2XFactoring.cpp:472),
   * which enumerates `x^n + ConvertBits(2i+1)` for i = 0, 1, 2, ... and
   * returns the first one passing `IterIrredTest`.
   * @param n - Degree (must be positive)
   * @returns The lexicographically smallest irreducible polynomial of degree n
   * @throws LogicError "BuildIrred: n must be positive"
   */
  static BuildIrred(n: number): GF2X {
    if (n <= 0) throw new Error('BuildIrred: n must be positive');

    if (n === 1) return GF2X.X();

    let i = 0n;
    for (;;) {
      if (i >> 63n) throw new Error('BuildIrred: limit exceeded');
      // ConvertBits(g, 2*i+1); SetCoeff(g, n)
      const g = GF2X._fromRep((2n * i + 1n) | (1n << BigInt(n)));
      i++;
      if (g.isIrreducible()) return g;
    }
  }

  /**
   * Builds an irreducible polynomial of degree n of minimal weight.
   *
   * Deterministic: ports NTL's `BuildSparseIrred`
   * (ntl/src/GF2XFactoring.cpp:900).  For n <= 2048 the answer comes from
   * NTL's precomputed table of minimal-weight trinomials/pentanomials; above
   * that bound NTL searches for an irreducible trinomial and then for an
   * irreducible pentanomial, falling back to `BuildIrred`.
   * @param n - Degree (must be positive)
   * @returns An irreducible polynomial of degree n of minimal weight
   * @throws LogicError "SparseIrred: n <= 0"
   */
  static BuildSparseIrred(n: number): GF2X {
    if (n <= 0) throw new Error('SparseIrred: n <= 0');

    if (n === 1) return GF2X.X();

    if (n <= GF2X_IRRED_TAB_MAX) {
      const row = GF2X_irred_tab[n]!;
      const f = new GF2X();
      f.SetCoeff(n);
      f.SetCoeff(row[0]);
      if (row[1] !== 0) {
        f.SetCoeff(row[1]);
        f.SetCoeff(row[2]);
      }
      f.SetCoeff(0);
      return f;
    }

    let k3 = FindTrinom(n);
    if (k3) {
      const f = new GF2X();
      f.SetCoeff(n);
      f.SetCoeff(k3);
      f.SetCoeff(0);
      return f;
    }

    const pent = FindPent(n);
    if (pent) {
      const [kk3, kk2, kk1] = pent;
      const f = new GF2X();
      f.SetCoeff(n);
      f.SetCoeff(kk3);
      f.SetCoeff(kk2);
      f.SetCoeff(kk1);
      f.SetCoeff(0);
      return f;
    }

    // The following is probably of only theoretical value: it is reasonable to
    // conjecture that for all n >= 2 there is either an irreducible trinomial
    // or pentanomial of degree n.
    return GF2X.BuildIrred(n);
  }

  // ============================================
  // Special Operations
  // ============================================

  /**
   * Computes derivative.
   * @returns The derivative
   */
  diff(): GF2X {
    // Over GF(2) only the odd-index coefficients survive, shifted down by one.
    let r = 0n;
    let x = this._xrep >> 1n;
    let i = 0n;
    while (x !== 0n) {
      // coefficient of x^(2i+1) contributes to x^(2i)
      if (x & 1n) r |= 1n << (2n * i);
      x >>= 2n;
      i++;
    }
    return GF2X._fromRep(r);
  }

  /**
   * Reverses coefficients (x^deg * f(1/x)).
   * @returns Reversed polynomial
   */
  reverse(): GF2X {
    const d = this.deg();
    if (d < 0) return new GF2X();
    let r = 0n;
    for (let i = 0; i <= d; i++) {
      if ((this._xrep >> BigInt(i)) & 1n) r |= 1n << BigInt(d - i);
    }
    return GF2X._fromRep(r);
  }

  /**
   * Left shift (multiply by X^n). Negative n shifts right.
   * @param n - Shift amount
   * @returns Shifted polynomial
   */
  LeftShift(n: number): GF2X {
    if (n < 0) return this.RightShift(-n);
    return GF2X._fromRep(this._xrep << BigInt(n));
  }

  /**
   * Right shift (divide by X^n). Negative n shifts left.
   * @param n - Shift amount
   * @returns Shifted polynomial
   */
  RightShift(n: number): GF2X {
    if (n < 0) return this.LeftShift(-n);
    return GF2X._fromRep(this._xrep >> BigInt(n));
  }

  /**
   * Truncates to degree < n.
   * @param n - Maximum degree + 1
   * @returns Truncated polynomial
   * @throws LogicError "trunc: bad args"
   */
  trunc(n: number): GF2X {
    if (n < 0) throw new Error('trunc: bad args');
    return GF2X._fromRep(this._xrep & ((1n << BigInt(n)) - 1n));
  }

  // ============================================
  // Comparison and Utility
  // ============================================

  /**
   * Checks if polynomial is zero.
   * @returns True if zero
   */
  IsZero(): boolean {
    return this._xrep === 0n;
  }

  /**
   * Checks if polynomial is one.
   * @returns True if one
   */
  IsOne(): boolean {
    return this._xrep === 1n;
  }

  /**
   * Checks if polynomial is X.
   * @returns True if X
   */
  IsX(): boolean {
    return this._xrep === 2n;
  }

  /**
   * Checks equality.
   * @param other - Polynomial to compare
   * @returns True if equal
   */
  equals(other: GF2X): boolean {
    return this._xrep === other._xrep;
  }

  /**
   * Converts to string, in NTL's format: the coefficient list from the
   * constant term up, e.g. `[1 1 1]` for x^2 + x + 1 (`[]` for zero).
   * Honours the static `HexOutput` flag.
   * @returns String representation
   */
  toString(): string {
    if (GF2X.HexOutput) return this.toHex();
    const d = this.deg();
    const bits: string[] = [];
    for (let i = 0; i <= d; i++) {
      bits.push(((this._xrep >> BigInt(i)) & 1n) === 1n ? '1' : '0');
    }
    return '[' + bits.join(' ') + ']';
  }

  /**
   * Converts to hex string in NTL's format (ntl/src/GF2X.cpp:378): coefficients
   * are packed four at a time into nibbles, least significant coefficient in
   * the least significant bit, and the nibbles are printed in increasing order.
   * @returns Hex representation, e.g. "0x31" for x^4 + x + 1
   */
  toHex(): string {
    const d = this.deg();
    if (d < 0) return '0x0';
    let out = '0x';
    let val = 0;
    let n = 0;
    for (let i = 0; i <= d; i++) {
      val = val | (Number((this._xrep >> BigInt(i)) & 1n) << n);
      n++;
      if (n === 4) {
        out += val.toString(16);
        val = 0;
        n = 0;
      }
    }
    if (val) out += val.toString(16);
    return out;
  }

  /**
   * Creates from hex string, inverting {@link GF2X.toHex} (ntl/src/GF2X.cpp:280).
   * @param hex - Hex string (e.g., "0x31")
   * @returns GF2X polynomial
   */
  static fromHex(hex: string): GF2X {
    let s = hex.trim();
    if (s.startsWith('0x') || s.startsWith('0X')) s = s.slice(2);
    let rep = 0n;
    let n = 0;
    for (const ch of s) {
      const val = parseInt(ch, 16);
      if (Number.isNaN(val)) throw new Error('bad GF2X input');
      for (let i = 0; i < 4; i++) {
        if (val & (1 << i)) rep |= 1n << BigInt(n + i);
      }
      n += 4;
    }
    return GF2X._fromRep(rep);
  }

  // ============================================
  // Random
  // ============================================

  /**
   * Generates a random polynomial of degree < n (n random coefficients),
   * using NTL's pseudo-random stream (ntl/src/GF2X.cpp:438).
   * @param n - Number of random coefficients; the result has degree < n
   * @returns Random polynomial of degree < n
   */
  static random(n: number): GF2X {
    throw new Error('NTL_NOT_IMPLEMENTED: GF2X.random');
  }
}

// ============================================
// Sparse irreducible searches (GF2XFactoring.cpp:868-898)
// ============================================

/**
 * Finds the smallest k with x^n + x^k + 1 irreducible, or 0 if none exists.
 * @see Reference: ntl/src/GF2XFactoring.cpp:868
 */
function FindTrinom(n: number): number {
  if (n < 2) throw new Error('tri--bad n');

  for (let k = 1; k <= n / 2; k++) {
    const f = new GF2X([]);
    f.SetCoeff(0);
    f.SetCoeff(k);
    f.SetCoeff(n);
    if (f.isIrreducible()) return k;
  }

  return 0;
}

/**
 * Finds [k3, k2, k1] with x^n + x^k3 + x^k2 + x^k1 + 1 irreducible,
 * or null if none exists.
 * @see Reference: ntl/src/GF2XFactoring.cpp:882
 */
function FindPent(n: number): [number, number, number] | null {
  if (n < 4) throw new Error('pent--bad n');

  for (let k3 = 3; k3 < n; k3++) {
    for (let k2 = 2; k2 < k3; k2++) {
      for (let k1 = 1; k1 < k2; k1++) {
        const f = new GF2X([]);
        f.SetCoeff(0);
        f.SetCoeff(k1);
        f.SetCoeff(k2);
        f.SetCoeff(k3);
        f.SetCoeff(n);
        if (f.isIrreducible()) return [k3, k2, k1];
      }
    }
  }

  return null;
}

// ============================================
// Standalone Functions (NTL-style API)
// ============================================

/**
 * Returns degree.
 */
export function deg(f: GF2X): number {
  return f.deg();
}

/**
 * Returns coefficient.
 */
export function coeff(f: GF2X, i: number): GF2 {
  return f.coeff(i);
}

/**
 * Returns leading coefficient.
 */
export function LeadCoeff(f: GF2X): GF2 {
  return f.LeadCoeff();
}

/**
 * Returns constant term.
 */
export function ConstTerm(f: GF2X): GF2 {
  return f.ConstTerm();
}

/**
 * Computes GCD.
 */
export function GCD(a: GF2X, b: GF2X): GF2X {
  return GF2X.GCD(a, b);
}

/**
 * Extended GCD.
 */
export function XGCD(a: GF2X, b: GF2X): [GF2X, GF2X, GF2X] {
  return GF2X.XGCD(a, b);
}

/**
 * Checks if zero.
 */
export function IsZero(f: GF2X): boolean {
  return f.IsZero();
}

/**
 * Checks if one.
 */
export function IsOne(f: GF2X): boolean {
  return f.IsOne();
}

/**
 * Checks if X.
 */
export function IsX(f: GF2X): boolean {
  return f.IsX();
}

/**
 * Computes modular product.
 */
export function MulMod(a: GF2X, b: GF2X, f: GF2X): GF2X {
  return GF2X.MulMod(a, b, f);
}

/**
 * Computes modular square.
 */
export function SqrMod(a: GF2X, f: GF2X): GF2X {
  return GF2X.SqrMod(a, f);
}

/**
 * Computes modular inverse.
 */
export function InvMod(a: GF2X, f: GF2X): GF2X {
  return GF2X.InvMod(a, f);
}

/**
 * Computes modular power.
 */
export function PowerMod(a: GF2X, e: bigint, f: GF2X): GF2X {
  return GF2X.PowerMod(a, e, f);
}

/**
 * Deterministic irreducibility test (NTL `IterIrredTest`).
 */
export function IterIrredTest(f: GF2X): boolean {
  return f.isIrreducible();
}

/**
 * Builds the lexicographically smallest irreducible polynomial of degree n.
 */
export function BuildIrred(n: number): GF2X {
  return GF2X.BuildIrred(n);
}

/**
 * Builds an irreducible polynomial of degree n of minimal weight.
 */
export function BuildSparseIrred(n: number): GF2X {
  return GF2X.BuildSparseIrred(n);
}
