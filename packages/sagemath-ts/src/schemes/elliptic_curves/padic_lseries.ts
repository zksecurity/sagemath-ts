/**
 * @module sage/schemes/elliptic_curves/padic_lseries
 * @description p-adic L-series of elliptic curves
 *
 * Port of: sage/schemes/elliptic_curves/padic_lseries.py
 *
 * This module provides p-adic L-series for elliptic curves over Q.
 * The p-adic L-function interpolates special values of the complex L-function.
 *
 * IMPORTANT: This is a partial implementation. The full SageMath implementation
 * requires:
 * - Modular symbols (sage.modular.modsym)
 * - p-adic rings and fields (sage.rings.padics)
 * - Power series rings (sage.rings.power_series_ring)
 * - Hyperelliptic curve computations (sage.schemes.hyperelliptic_curves.monsky_washnitzer)
 *
 * Many of these dependencies are not yet available in sagemath-ts.
 *
 * @see Deviation: Elliptic Curve p-adic L-series and Isogeny Class Partial Implementation
 */

import {
  binomial,
  gcd,
  is_prime,
  is_squarefree,
  kronecker_symbol,
  valuation,
} from '../../arith/misc.js';
import { NotImplementedError, ValueError } from '../../errors.js';
import { Qp, pAdicField } from '../../rings/padics/padic_generic.js';
import { pAdicGenericElement } from '../../rings/padics/padic_generic_element.js';
import {
  type CoefficientRing,
  LaurentSeriesRing,
  type PowerSeriesElement,
  PowerSeriesRing,
  type RingElement,
} from '../../rings/power_series_ring.js';
import type { EllipticCurveGeneric } from './ell_generic.js';
import type { EllipticCurveFormalGroup } from './formal_group.js';
import type { FieldElement } from './types.js';

/**
 * Rational number representation for modular symbols.
 */
export interface Rational {
  num: bigint;
  den: bigint;
}

/**
 * Create a rational number.
 */
export function rational(num: bigint, den: bigint = 1n): Rational {
  if (den === 0n) {
    throw new ValueError('denominator cannot be zero');
  }
  // Normalize sign
  if (den < 0n) {
    num = -num;
    den = -den;
  }
  // Reduce
  const g = gcd(num < 0n ? -num : num, den);
  return { num: num / g, den: den / g };
}

/**
 * A rational number element (for power series coefficients).
 */
export class RationalElement implements RingElement {
  readonly num: bigint;
  readonly den: bigint;

  constructor(num: bigint, den: bigint = 1n) {
    if (den === 0n) {
      throw new ValueError('denominator cannot be zero');
    }
    // Normalize
    if (den < 0n) {
      num = -num;
      den = -den;
    }
    const g = gcd(num < 0n ? -num : num, den);
    this.num = num / g;
    this.den = den / g;
  }

  add(other: RingElement): RationalElement {
    const o = other as RationalElement;
    return new RationalElement(this.num * o.den + o.num * this.den, this.den * o.den);
  }

  sub(other: RingElement): RationalElement {
    const o = other as RationalElement;
    return new RationalElement(this.num * o.den - o.num * this.den, this.den * o.den);
  }

  mul(other: RingElement): RationalElement {
    const o = other as RationalElement;
    return new RationalElement(this.num * o.num, this.den * o.den);
  }

  div(other: RingElement): RationalElement {
    const o = other as RationalElement;
    if (o.num === 0n) {
      throw new ValueError('division by zero');
    }
    return new RationalElement(this.num * o.den, this.den * o.num);
  }

  neg(): RationalElement {
    return new RationalElement(-this.num, this.den);
  }

  eq(other: RingElement | number | bigint): boolean {
    if (typeof other === 'number') {
      return this.num === BigInt(other) * this.den;
    }
    if (typeof other === 'bigint') {
      return this.num === other * this.den;
    }
    const o = other as RationalElement;
    return this.num * o.den === o.num * this.den;
  }

  isZero(): boolean {
    return this.num === 0n;
  }

  isOne(): boolean {
    return this.num === this.den;
  }

  isUnit(): boolean {
    return this.num !== 0n;
  }

  inv(): RationalElement {
    if (this.num === 0n) {
      throw new ValueError('cannot invert zero');
    }
    return new RationalElement(this.den, this.num);
  }

  toString(): string {
    if (this.den === 1n) {
      return this.num.toString();
    }
    return `${this.num}/${this.den}`;
  }
}

/**
 * The field of rational numbers Q.
 */
export class RationalRing implements CoefficientRing<RationalElement> {
  zero(): RationalElement {
    return new RationalElement(0n, 1n);
  }

  one(): RationalElement {
    return new RationalElement(1n, 1n);
  }

  __call__(x: unknown): RationalElement {
    if (x instanceof RationalElement) {
      return x;
    }
    if (typeof x === 'bigint') {
      return new RationalElement(x, 1n);
    }
    if (typeof x === 'number') {
      return new RationalElement(BigInt(x), 1n);
    }
    if (typeof x === 'object' && x !== null) {
      if ('num' in x && 'den' in x) {
        const r = x as { num: bigint; den: bigint };
        return new RationalElement(BigInt(r.num), BigInt(r.den));
      }
    }
    throw new ValueError(`cannot coerce ${typeof x} to rational`);
  }

  is_field(): boolean {
    return true;
  }

  characteristic(): bigint {
    return 0n;
  }
}

/**
 * Options for p-adic L-series construction.
 */
export interface pAdicLseriesOptions {
  /** Implementation to use: 'eclib' or 'sage' */
  implementation?: 'eclib' | 'sage' | 'num';
  /** Normalization: 'L_ratio', 'period', or 'none' */
  normalize?: 'L_ratio' | 'period' | 'none';
}

/**
 * A 2x2 matrix representing Frobenius on the Dieudonne module.
 */
export interface FrobeniusMatrix {
  /** The matrix entries [[a, b], [c, d]] */
  entries: [[pAdicGenericElement, pAdicGenericElement], [pAdicGenericElement, pAdicGenericElement]];
  /** The precision of the entries */
  precision: number;
  /** String representation */
  toString(): string;
}

/**
 * Check if D is a fundamental discriminant.
 *
 * A fundamental discriminant is either:
 * - D = 1 (mod 4) and squarefree
 * - D = 0 (mod 4), D/4 squarefree, and D/4 = 2 or 3 (mod 4)
 *
 * @param D - Integer to check
 * @returns true if D is a fundamental discriminant
 */
function isFundamentalDiscriminant(D: bigint): boolean {
  if (D === 1n) return true; // Trivial case

  const mod4 = ((D % 4n) + 4n) % 4n; // Proper mod for negative numbers

  if (mod4 === 1n) {
    // D = 1 (mod 4): D must be squarefree
    return is_squarefree(D);
  } else if (mod4 === 0n) {
    // D = 0 (mod 4): D/4 must be squarefree and D/4 = 2 or 3 (mod 4)
    const d = D / 4n;
    const dMod4 = ((d % 4n) + 4n) % 4n;
    return is_squarefree(d) && (dMod4 === 2n || dMod4 === 3n);
  }

  return false; // D = 2 or 3 (mod 4) is not fundamental
}

/** p-adic valuation of a nonzero bigint; `Infinity` for 0. */
function _vp(n: bigint, p: bigint): number {
  if (n === 0n) return Number.POSITIVE_INFINITY;
  let v = 0;
  let m = n < 0n ? -n : n;
  while (m % p === 0n) {
    m /= p;
    v++;
  }
  return v;
}

/** Modular inverse of a p-adic unit modulo p^k. */
function _unitInverse(u: bigint, p: bigint, k: number): bigint {
  const m = p ** BigInt(k);
  let [oldR, r] = [((u % m) + m) % m, m];
  let [oldS, s] = [1n, 0n];
  while (r !== 0n) {
    const q = oldR / r;
    [oldR, r] = [r, oldR - q * r];
    [oldS, s] = [s, oldS - q * s];
  }
  if (oldR !== 1n && oldR !== -1n) {
    throw new ValueError('not a p-adic unit');
  }
  const inv = oldR === 1n ? oldS : -oldS;
  return ((inv % m) + m) % m;
}

/**
 * The totally ramified quadratic extension `A = K[x]/(f)` of a p-adic field
 * `K = Qp(p, prec)` by an Eisenstein polynomial `f = x^2 - a_p*x + p`.
 *
 * This is the object SageMath builds with `K.extension(f, names='alpha')` in
 * `pAdicLseries.alpha` for a prime of supersingular reduction (there
 * `v_p(a_p) >= 1`, so `f` is Eisenstein and `A/K` is totally ramified of
 * degree 2 with uniformizer `alpha`).
 *
 * @see Reference: sage/schemes/elliptic_curves/padic_lseries.py:513-518
 * @see Reference: sage/rings/padics/padic_extension_leaves.py (EisensteinExtensionGeneric)
 * @see Deviation: the port's `sage/rings/padics` has a `pAdicExtension` shell with
 *   no element type, so the Eisenstein extension needed by `alpha()` is built here.
 */
export class pAdicEisensteinQuadraticExtension {
  private readonly _p: bigint;
  private readonly _ap: bigint;
  private readonly _prec: number;
  private readonly _name: string;

  /**
   * @param p - the residue characteristic
   * @param ap - the trace of Frobenius; must satisfy `p | ap` (Eisenstein)
   * @param prec - the precision cap of the base field `K = Qp(p, prec)`
   * @param name - the name of the generator (SageMath uses 'alpha')
   */
  constructor(p: bigint, ap: bigint, prec: number, name: string = 'alpha') {
    if (ap % p !== 0n) {
      throw new ValueError(`x^2 - ${ap}*x + ${p} is not Eisenstein at ${p}`);
    }
    this._p = p;
    this._ap = ap;
    this._prec = prec;
    this._name = name;
  }

  prime(): bigint {
    return this._p;
  }

  /** The `a_p` of the defining polynomial `x^2 - a_p*x + p`. */
  ap(): bigint {
    return this._ap;
  }

  variable_name(): string {
    return this._name;
  }

  /** The degree of `A` over `K`. */
  degree(): number {
    return 2;
  }

  /** The ramification index `e(A/K) = 2`. */
  e(): number {
    return 2;
  }

  /** The residue degree `f(A/K) = 1`. */
  f(): number {
    return 1;
  }

  /**
   * The precision cap, in powers of the uniformizer: `e * prec`.
   * (SageMath: an extension of ramification index `e` of a field of precision
   * cap `prec` has precision cap `e*prec`.)
   */
  precision_cap(): number {
    return 2 * this._prec;
  }

  /** The coefficients `[p, -a_p, 1]` of the defining polynomial. */
  defining_polynomial(): [bigint, bigint, bigint] {
    return [this._p, -this._ap, 1n];
  }

  /** The uniformizer `alpha`, i.e. the root of the defining polynomial. */
  gen(): pAdicEisensteinQuadraticElement {
    return new pAdicEisensteinQuadraticElement(this, 0n, 1n, 0, 1 + this.precision_cap());
  }

  zero(): pAdicEisensteinQuadraticElement {
    return new pAdicEisensteinQuadraticElement(this, 0n, 0n, 0, Number.POSITIVE_INFINITY);
  }

  one(): pAdicEisensteinQuadraticElement {
    return this.__call__(1n);
  }

  /**
   * Coerce a rational integer (or a `pAdicGenericElement` of `K`) into `A`.
   * An exact integer gets the full relative precision cap, as in SageMath.
   */
  __call__(x: bigint | number | pAdicGenericElement): pAdicEisensteinQuadraticElement {
    if (typeof x === 'number') return this.__call__(BigInt(x));
    if (typeof x === 'bigint') {
      if (x === 0n) return this.zero();
      const v = 2 * _vp(x, this._p);
      return new pAdicEisensteinQuadraticElement(this, x, 0n, 0, v + this.precision_cap());
    }
    // A capped-relative element of K: relative precision doubles in A.
    const lifted = x.lift();
    if (lifted === 0n) return this.zero();
    const v = 2 * _vp(lifted, this._p);
    const relprec = 2 * x.precision_relative();
    return new pAdicEisensteinQuadraticElement(
      this,
      lifted,
      0n,
      0,
      v + Math.min(relprec, this.precision_cap())
    );
  }

  toString(): string {
    return (
      `Eisenstein Extension in ${this._name} defined by x^2 - ${this._ap}*x + ${this._p} ` +
      `of ${this._p}-adic Field with capped relative precision ${this._prec}`
    );
  }
}

/**
 * An element `(c0 + c1*alpha) / p^den` of a
 * {@link pAdicEisensteinQuadraticExtension}, with an absolute precision
 * measured in powers of the uniformizer `alpha` (SageMath prints these as
 * `O(alpha^n)`).
 *
 * @see Deviation: see {@link pAdicEisensteinQuadraticExtension}.
 */
export class pAdicEisensteinQuadraticElement {
  private readonly _parent: pAdicEisensteinQuadraticExtension;
  private readonly _c0: bigint;
  private readonly _c1: bigint;
  private readonly _den: number;
  private readonly _absprec: number;

  constructor(
    parent: pAdicEisensteinQuadraticExtension,
    c0: bigint,
    c1: bigint,
    den: number,
    absprec: number
  ) {
    const p = parent.prime();
    // Normalize: cancel common powers of p between (c0, c1) and p^den.
    while (den > 0 && c0 % p === 0n && c1 % p === 0n) {
      c0 /= p;
      c1 /= p;
      den--;
    }
    this._parent = parent;
    this._c0 = c0;
    this._c1 = c1;
    this._den = den;
    this._absprec = absprec;
  }

  parent(): pAdicEisensteinQuadraticExtension {
    return this._parent;
  }

  prime(): bigint {
    return this._parent.prime();
  }

  /** The valuation in powers of `alpha` (`v(p) = 2`). */
  valuation(): number {
    if (this._c0 === 0n && this._c1 === 0n) return this._absprec;
    const p = this.prime();
    const v = Math.min(2 * _vp(this._c0, p), 2 * _vp(this._c1, p) + 1) - 2 * this._den;
    return Math.min(v, this._absprec);
  }

  /** The absolute precision, in powers of `alpha`. */
  precision_absolute(): number {
    return this._absprec;
  }

  /** The relative precision, in powers of `alpha`. */
  precision_relative(): number {
    if (this._absprec === Number.POSITIVE_INFINITY) return Number.POSITIVE_INFINITY;
    if (this._c0 === 0n && this._c1 === 0n) return 0;
    return this._absprec - this.valuation();
  }

  is_zero(): boolean {
    return this.valuation() >= this._absprec;
  }

  /** Cap the relative precision at the ring's precision cap. */
  private _capped(c0: bigint, c1: bigint, den: number, absprec: number) {
    const e = new pAdicEisensteinQuadraticElement(this._parent, c0, c1, den, absprec);
    if (absprec === Number.POSITIVE_INFINITY) return e;
    const cap = this._parent.precision_cap();
    const v = e.valuation();
    if (v !== Number.POSITIVE_INFINITY && absprec - v > cap) {
      return new pAdicEisensteinQuadraticElement(this._parent, c0, c1, den, v + cap);
    }
    return e;
  }

  add(other: pAdicEisensteinQuadraticElement): pAdicEisensteinQuadraticElement {
    const p = this.prime();
    const den = Math.max(this._den, other._den);
    const s = p ** BigInt(den - this._den);
    const t = p ** BigInt(den - other._den);
    return this._capped(
      this._c0 * s + other._c0 * t,
      this._c1 * s + other._c1 * t,
      den,
      Math.min(this._absprec, other._absprec)
    );
  }

  neg(): pAdicEisensteinQuadraticElement {
    return new pAdicEisensteinQuadraticElement(
      this._parent,
      -this._c0,
      -this._c1,
      this._den,
      this._absprec
    );
  }

  sub(other: pAdicEisensteinQuadraticElement): pAdicEisensteinQuadraticElement {
    return this.add(other.neg());
  }

  mul(other: pAdicEisensteinQuadraticElement): pAdicEisensteinQuadraticElement {
    const p = this.prime();
    const ap = this._parent.ap();
    // (a + b*alpha)(c + d*alpha) = (ac - p*bd) + (ad + bc + a_p*bd)*alpha,
    // using alpha^2 = a_p*alpha - p.
    const a = this._c0;
    const b = this._c1;
    const c = other._c0;
    const d = other._c1;
    const c0 = a * c - p * b * d;
    const c1 = a * d + b * c + ap * b * d;
    const va = this.valuation();
    const vb = other.valuation();
    const absprec = Math.min(va + other._absprec, vb + this._absprec);
    return this._capped(c0, c1, this._den + other._den, absprec);
  }

  /**
   * The norm `N_{A/K}(c0 + c1*alpha) = c0^2 + a_p*c0*c1 + p*c1^2` of the
   * numerator (an element of `K`).
   */
  private _numeratorNorm(): bigint {
    const p = this.prime();
    const ap = this._parent.ap();
    return this._c0 * this._c0 + ap * this._c0 * this._c1 + p * this._c1 * this._c1;
  }

  inv(): pAdicEisensteinQuadraticElement {
    if (this.is_zero()) {
      throw new ValueError('cannot invert zero');
    }
    const p = this.prime();
    const ap = this._parent.ap();
    // 1/x = p^den * conj / N, with conj = (c0 + a_p*c1) - c1*alpha.
    const N = this._numeratorNorm();
    const w = _vp(N, p);
    const U = N / p ** BigInt(w);
    const relprec = Math.min(this.precision_relative(), this._parent.precision_cap());
    // Enough p-adic digits to carry `relprec` powers of alpha (v(p) = 2).
    const k = Math.ceil(relprec / 2) + 2 + w;
    const Uinv = _unitInverse(U < 0n ? -U : U, p, k) * (U < 0n ? -1n : 1n);
    const mod = p ** BigInt(k);
    const c0 = (((this._c0 + ap * this._c1) * Uinv) % mod) as bigint;
    const c1 = ((-this._c1 * Uinv) % mod) as bigint;
    const v = -this.valuation();
    return this._capped(c0, c1, w - this._den, v + relprec);
  }

  div(other: pAdicEisensteinQuadraticElement): pAdicEisensteinQuadraticElement {
    return this.mul(other.inv());
  }

  pow(n: number): pAdicEisensteinQuadraticElement {
    if (n < 0) return this.pow(-n).inv();
    let result = this._parent.one();
    let base: pAdicEisensteinQuadraticElement = this;
    let e = n;
    while (e > 0) {
      if (e & 1) result = result.mul(base);
      base = base.mul(base);
      e >>= 1;
    }
    return result;
  }

  eq(other: pAdicEisensteinQuadraticElement): boolean {
    return this.sub(other).is_zero();
  }

  /**
   * The `alpha`-adic expansion: `[[exponent, digit], ...]` with digits in
   * `[0, p)`, from the valuation up to (but excluding) the absolute precision.
   */
  expansion(): Array<[number, bigint]> {
    const p = this.prime();
    const ap = this._parent.ap();
    if (this._absprec === Number.POSITIVE_INFINITY) {
      throw new ValueError('cannot expand an exact element to infinite precision');
    }
    // Clear the denominator: alpha^(2*den) / p^den = (s*alpha - 1)^den with
    // s = a_p/p (an integer because f is Eisenstein), so multiplying by that
    // integral element shifts every exponent by 2*den.
    let c0 = this._c0;
    let c1 = this._c1;
    const shift = 2 * this._den;
    if (this._den > 0) {
      const s = ap / p;
      let f0 = 1n;
      let f1 = 0n;
      for (let i = 0; i < this._den; i++) {
        // (f0 + f1*alpha) * (s*alpha - 1)
        const n0 = -f0 - p * f1 * s;
        const n1 = f0 * s - f1 + ap * f1 * s;
        f0 = n0;
        f1 = n1;
      }
      const n0 = c0 * f0 - p * c1 * f1;
      const n1 = c0 * f1 + c1 * f0 + ap * c1 * f1;
      c0 = n0;
      c1 = n1;
    }
    const out: Array<[number, bigint]> = [];
    for (let k = 0; k < this._absprec + shift; k++) {
      const d = ((c0 % p) + p) % p;
      if (d !== 0n) out.push([k - shift, d]);
      // (x - d)/alpha, using p/alpha = a_p - alpha
      const m = (c0 - d) / p;
      const n0 = m * ap + c1;
      const n1 = -m;
      c0 = n0;
      c1 = n1;
    }
    return out;
  }

  /** SageMath's `series` print mode for a ramified extension. */
  toString(): string {
    const name = this._parent.variable_name();
    if (this._absprec === Number.POSITIVE_INFINITY) {
      return this._c0 === 0n && this._c1 === 0n ? '0' : `${this._c0} + ${this._c1}*${name}`;
    }
    const parts: string[] = [];
    for (const [k, d] of this.expansion()) {
      const mon = k === 0 ? '' : k === 1 ? name : `${name}^${k}`;
      if (mon === '') parts.push(`${d}`);
      else parts.push(d === 1n ? mon : `${d}*${mon}`);
    }
    parts.push(`O(${name}^${this._absprec})`);
    return parts.join(' + ');
  }
}

/**
 * The p-adic L-series of an elliptic curve.
 *
 * This class computes the p-adic L-function attached to an elliptic curve E/Q.
 * The p-adic L-function is an element of the Iwasawa algebra that interpolates
 * the special values L(E, chi, 1) for Dirichlet characters chi of p-power conductor.
 *
 * For an elliptic curve E over Q and a prime p of good reduction, the p-adic
 * L-function L_p(E, s) is a p-adic analytic function that satisfies:
 *
 *   L_p(E, 1) = (1 - 1/alpha)^2 * L(E, 1) / Omega_E
 *
 * where alpha is the unit root of x^2 - a_p*x + p and Omega_E is the Neron period.
 *
 * @see Reference: sage/schemes/elliptic_curves/padic_lseries.py:pAdicLseries
 * @see Deviation: Elliptic Curve p-adic L-series and Isogeny Class Partial Implementation
 */
export class pAdicLseries<F extends FieldElement = FieldElement> {
  /** The elliptic curve */
  protected _E: EllipticCurveGeneric<F>;

  /** The prime p */
  protected _p: bigint;

  /** Implementation to use for modular symbols */
  protected _implementation: 'eclib' | 'sage' | 'num';

  /** Normalization method */
  protected _normalize: 'L_ratio' | 'period' | 'none';

  /** Cached alpha values by precision */
  protected _alpha: Map<number, pAdicGenericElement | pAdicEisensteinQuadraticElement> = new Map();

  /** Cached order of vanishing */
  protected __ord?: number;

  /** Cached series by (n, prec, D, eta) */
  protected __series: Map<string, PowerSeriesElement<pAdicGenericElement>> = new Map();

  /** Cached measure data */
  protected __measure_data: Map<
    string,
    [
      bigint,
      pAdicGenericElement,
      pAdicGenericElement,
      bigint,
      (r: Rational | number | bigint) => Rational,
    ]
  > = new Map();

  /**
   * Create the p-adic L-series for an elliptic curve.
   *
   * INPUT:
   * - E: an elliptic curve over Q
   * - p: a prime number of good or multiplicative reduction
   * - implementation: 'eclib', 'num', or 'sage' (default 'eclib')
   * - normalize: 'L_ratio' (default), 'period', or 'none'
   *
   * The prime p must be of semi-stable reduction, meaning the conductor
   * of E is not divisible by p^2.
   *
   * @example
   * ```typescript
   * const E = EllipticCurve([0, 0, 1, -1, 0]); // 11a1
   * const L = new pAdicLseries(E, 5n);
   * ```
   *
   * @see Reference: sage/schemes/elliptic_curves/padic_lseries.py:__init__
   */
  constructor(E: EllipticCurveGeneric<F>, p: bigint | number, options: pAdicLseriesOptions = {}) {
    const { implementation = 'eclib', normalize = 'L_ratio' } = options;

    this._E = E;
    this._p = BigInt(p);
    this._normalize = normalize;

    if (!['eclib', 'sage', 'num'].includes(implementation)) {
      throw new ValueError("Implementation should be one of 'eclib', 'num' or 'sage'");
    }
    this._implementation = implementation as 'eclib' | 'sage' | 'num';

    // Validate prime using the arith is_prime function
    if (!is_prime(this._p)) {
      throw new ValueError(`p (=${p}) must be a prime`);
    }

    // Check for semi-stable reduction: the conductor must not be divisible by p^2
    // (padic_lseries.py:182-183).  `_c_bound()` and the multiplicative-reduction
    // branch of `alpha()` both rely on this invariant.
    const N = this._conductor();
    if (N !== null && N % (this._p * this._p) === 0n) {
      throw new NotImplementedError(`p (=${p}) must be a prime of semi-stable reduction`);
    }
  }

  /**
   * Return the conductor of the curve, or null when the curve does not expose
   * a usable `conductor()` method.
   */
  private _conductor(): bigint | null {
    const E = this._E as unknown as { conductor?: () => bigint | number };
    if (typeof E.conductor !== 'function') {
      return null;
    }
    try {
      const N = E.conductor();
      return typeof N === 'bigint' ? N : BigInt(N);
    } catch {
      return null;
    }
  }

  /**
   * Return the elliptic curve to which this p-adic L-series is associated.
   *
   * @example
   * ```typescript
   * const L = E.padic_lseries(5n);
   * L.elliptic_curve() === E; // true
   * ```
   *
   * @see Reference: sage/schemes/elliptic_curves/padic_lseries.py:elliptic_curve
   */
  elliptic_curve(): EllipticCurveGeneric<F> {
    return this._E;
  }

  /**
   * Return the prime p.
   *
   * @example
   * ```typescript
   * const L = E.padic_lseries(5n);
   * L.prime(); // 5n
   * ```
   *
   * @see Reference: sage/schemes/elliptic_curves/padic_lseries.py:prime
   */
  prime(): bigint {
    return this._p;
  }

  /**
   * Compute the modular symbol for r with given sign.
   *
   * The modular symbol [r]^{sign} is defined as an integral of the
   * differential omega_E over a path from r to infinity.
   *
   * INPUT:
   * - r: a rational number (as Rational or { num, den })
   * - sign: +1 or -1 (default +1)
   * - quadratic_twist: a squarefree integer (default +1)
   *
   * OUTPUT: the modular symbol [r]^{sign}
   *
   * NOTE: This requires modular symbol computation which is not yet
   * implemented in sagemath-ts. The full implementation would require:
   * - Modular symbol space construction
   * - Boundary symbol computation
   * - Manin symbol representation
   *
   * @see Reference: sage/schemes/elliptic_curves/padic_lseries.py:modular_symbol
   */
  modular_symbol(
    r: Rational | number | bigint,
    sign: 1 | -1 = 1,
    quadratic_twist: bigint | number = 1
  ): Rational {
    // Reference: sage/schemes/elliptic_curves/padic_lseries.py:modular_symbol
    // This requires the modular symbol space computation
    // which depends on sage.modular.modsym
    //
    // The algorithm:
    // 1. If quadratic_twist == 1 and sign == +1, use self._modular_symbol(r)
    // 2. If sign == -1, need negative modular symbol space
    // 3. For D > 0, sum over kronecker(D, u) * m(r + u/D) for u = 1..D-1
    // 4. For D < 0, use negative symbol with sum over kronecker(D, u)

    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: pAdicLseries.modular_symbol - requires modular symbols (sage.modular.modsym)'
    );
  }

  /**
   * Compute the measure of a/p^n.
   *
   * The measure mu(a + p^n Z_p) is used in the definition of the p-adic L-function.
   * It is computed using modular symbols:
   *
   *   mu(a + p^n Z_p) = (1/alpha^n) * [a/p^n]^{sign} - (1/alpha^{n+1}) * [a/p^{n-1}]^{sign}
   *
   * INPUT:
   * - a: an integer
   * - n: a nonnegative integer
   * - prec: p-adic precision
   * - quadratic_twist: a squarefree integer (default +1)
   * - sign: +1 or -1 (default +1)
   *
   * OUTPUT: the measure mu(a + p^n Z_p) as a p-adic number
   *
   * @see Reference: sage/schemes/elliptic_curves/padic_lseries.py:measure
   */
  measure(
    a: bigint | number,
    n: number,
    prec: number,
    quadratic_twist: bigint | number = 1,
    sign: 1 | -1 = 1
  ): pAdicGenericElement {
    // Reference: sage/schemes/elliptic_curves/padic_lseries.py:measure
    // Algorithm:
    // 1. Validate sign is +/-1 and quadratic_twist constraints
    // 2. Get cached measure data or compute: p, alpha, z=1/alpha^n, w=p^(n-1), f=modular_symbol
    // 3. If quadratic_twist == 1:
    //    - If bad reduction at p: return z * f(a/(p*w))
    //    - Else: return z * (f(a/(p*w)) - f(a/w) / alpha)
    // 4. For quadratic_twist D != 1:
    //    - Compute chip = kronecker(D, p) (1 for supersingular)
    //    - Sum over u in range(1, |D|) with kronecker(D, u) weights

    const s = BigInt(sign);
    if (s !== 1n && s !== -1n) {
      throw new ValueError('Sign must be +- 1');
    }
    if (quadratic_twist !== 1 && s !== 1n) {
      throw new NotImplementedError('Quadratic twists not implemented for sign -1');
    }

    // This requires modular symbols and alpha computation
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: pAdicLseries.measure - requires modular symbols and p-adic arithmetic'
    );
  }

  /**
   * Compute the unit root of the characteristic polynomial of Frobenius.
   *
   * For a prime p of good ordinary reduction, the characteristic polynomial
   * of Frobenius is x^2 - a_p*x + p, which has two roots. The unit root alpha
   * is the one with p-adic valuation 0.
   *
   * For supersingular primes, alpha is a root in a quadratic extension of Q_p.
   *
   * INPUT:
   * - prec: p-adic precision (default 20)
   *
   * OUTPUT: the unit root alpha of x^2 - a_p*x + p
   *
   * ALGORITHM:
   * 1. For multiplicative reduction (p | conductor), alpha = a_p
   * 2. For good ordinary reduction, factor x^2 - a_p*x + p over Q_p
   *    and return the root with valuation < 1
   * 3. For supersingular reduction, return a root in a quadratic
   *    extension of Q_p
   *
   * @see Reference: sage/schemes/elliptic_curves/padic_lseries.py:alpha
   */
  alpha(prec: number = 20): pAdicGenericElement | pAdicEisensteinQuadraticElement {
    // Check cache
    const cached = this._alpha.get(prec);
    if (cached) {
      return cached;
    }

    const p = this._p;
    const K = Qp(p, prec);

    // Get a_p (trace of Frobenius)
    // This requires the ap() method on the elliptic curve
    let a_p: bigint;
    try {
      // Try to get a_p from the elliptic curve
      // The curve needs to have an ap() method
      if (
        'ap' in this._E &&
        typeof (this._E as unknown as { ap: (p: bigint) => bigint }).ap === 'function'
      ) {
        a_p = (this._E as unknown as { ap: (p: bigint) => bigint }).ap(p);
      } else {
        throw new NotImplementedError('Curve does not have ap() method');
      }
    } catch {
      throw new NotImplementedError(
        'SAGE_NOT_IMPLEMENTED: pAdicLseries.alpha - requires trace of Frobenius computation (ap)'
      );
    }

    // Multiplicative reduction (p divides the conductor): alpha = a_p
    // (padic_lseries.py:512-514).
    const N = this._conductor();
    if (N !== null && N % p === 0n) {
      const res = K.__call__(a_p);
      this._alpha.set(prec, res);
      return res;
    }

    // For good ordinary reduction: factor x^2 - a_p*x + p
    // The polynomial is x^2 - a_p*x + p
    // Roots are (a_p +/- sqrt(a_p^2 - 4p)) / 2

    // Check if ordinary: a_p not divisible by p
    const isOrdinary = a_p % p !== 0n;

    if (isOrdinary) {
      // Ordinary case: find the unit root using Hensel lifting
      // Start with a_p mod p (which is a unit)
      // The unit root satisfies alpha = a_p - p/alpha
      // Newton iteration: alpha' = alpha - (alpha^2 - a_p*alpha + p) / (2*alpha - a_p)

      // Initial approximation: solve x^2 - a_p*x + p = 0 mod p
      // This gives x = a_p mod p (since p = 0 mod p)
      let alpha_lift = ((a_p % p) + p) % p;
      if (alpha_lift === 0n) {
        alpha_lift = 1n; // Fallback
      }

      // Hensel lift to precision prec
      let modulus = p;
      for (let i = 1; i < prec; i++) {
        const newModulus = modulus * p;
        // f(x) = x^2 - a_p*x + p
        // f'(x) = 2x - a_p
        const fx = (alpha_lift * alpha_lift - a_p * alpha_lift + p) % newModulus;
        const fpx = (2n * alpha_lift - a_p) % newModulus;

        // Newton: x' = x - f(x)/f'(x)
        // Need to compute f(x)/f'(x) mod p^i
        const fpx_inv = this._modInverse(fpx, newModulus);
        if (fpx_inv === null) {
          throw new ValueError('Hensel lifting failed: derivative not invertible');
        }
        alpha_lift = (((alpha_lift - fx * fpx_inv) % newModulus) + newModulus) % newModulus;
        modulus = newModulus;
      }

      const result = K.__call__(alpha_lift);
      this._alpha.set(prec, result);
      return result;
    } else {
      // Supersingular case (padic_lseries.py:513-518):
      //     f = f.change_ring(K)
      //     A = K.extension(f, names='alpha')
      //     a = A.gen()
      // f = x^2 - a_p*x + p is Eisenstein here (p | a_p), so A/K is totally
      // ramified of degree 2 and alpha is a uniformizer.
      const A = new pAdicEisensteinQuadraticExtension(p, a_p, prec, 'alpha');
      const a = A.gen();
      this._alpha.set(prec, a);
      return a;
    }
  }

  /**
   * Compute modular inverse of a mod m.
   */
  private _modInverse(a: bigint, m: bigint): bigint | null {
    let [oldR, r] = [a % m, m];
    let [oldS, s] = [1n, 0n];

    while (r !== 0n) {
      const q = oldR / r;
      [oldR, r] = [r, oldR - q * r];
      [oldS, s] = [s, oldS - q * s];
    }

    if (oldR > 1n && oldR !== -1n) {
      return null; // No inverse
    }
    if (oldS < 0n) {
      oldS += m;
    }
    return oldS;
  }

  /**
   * Return the order of vanishing of the p-adic L-function at s=0.
   *
   * By a theorem of Kato, the order of vanishing is at least the rank of E(Q).
   * The p-adic BSD conjecture predicts equality (with a correction term for
   * split multiplicative reduction).
   *
   * OUTPUT: a non-negative integer
   *
   * NOTE: Currently only implemented for ordinary primes.
   *
   * ALGORITHM:
   * Compute successive approximations to the series and find the valuation.
   * By Kato's theorem, if the series vanishes to order v < r (where r is the rank),
   * we get a contradiction.
   *
   * @see Reference: sage/schemes/elliptic_curves/padic_lseries.py:order_of_vanishing
   */
  order_of_vanishing(): number {
    if (this.__ord !== undefined) {
      return this.__ord;
    }

    if (!this.is_ordinary()) {
      throw new NotImplementedError('order_of_vanishing only implemented for ordinary primes');
    }

    // This requires computing the series and finding its valuation
    // Also needs the rank of E(Q)
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: pAdicLseries.order_of_vanishing - requires series computation and rank'
    );
  }

  /**
   * Return Teichmuller lifts to the given precision.
   *
   * The Teichmuller character tau: (Z/pZ)* -> Z_p* is the unique multiplicative
   * section of the reduction map Z_p* -> (Z/pZ)*.
   *
   * INPUT:
   * - prec: positive integer
   *
   * OUTPUT: list of p-adic numbers; the cached Teichmuller lifts [0, tau(1), tau(2), ..., tau(p-1)]
   *
   * ALGORITHM:
   * For a in {1, ..., p-1}, the Teichmuller lift tau(a) is computed by iterating
   * the p-th power map: tau(a) = lim_{n->infinity} a^{p^n}
   * In practice, we compute a^{p^n} mod p^{prec} for large enough n.
   *
   * @see Reference: sage/schemes/elliptic_curves/padic_lseries.py:teichmuller
   */
  teichmuller(prec: number = 20): bigint[] {
    const p = this._p;
    const K = Qp(p, prec);
    const result: bigint[] = [0n];

    // Compute Teichmuller lifts for 1, 2, ..., p-1
    const teichSystem = K.teichmuller_system();
    for (const t of teichSystem) {
      // Get the residue at the given precision
      result.push(t.residue(prec));
    }

    return result;
  }

  /**
   * Compute the bounds on p-adic precision for series coefficients.
   *
   * This computes the valuations of the coefficients of omega_n = (1+T)^{p^n} - 1.
   *
   * INPUT:
   * - n: positive integer
   * - prec: number of terms
   *
   * OUTPUT: list of valuations [infinity, e_1, e_2, ..., e_{prec-1}]
   *
   * @see Reference: sage/schemes/elliptic_curves/padic_lseries.py:_e_bounds
   */
  protected _e_bounds(n: number, prec: number): number[] {
    // Reference: sage/schemes/elliptic_curves/padic_lseries.py:_e_bounds
    // Computes the valuations of binomial(p^n, j) for j = 0, 1, ..., prec-1
    // The sequence must be decreasing

    const pn = this._p ** BigInt(n);
    let enj = Number.POSITIVE_INFINITY;
    const res: number[] = [enj];

    for (let j = 1; j < prec; j++) {
      const bino = binomial(pn, BigInt(j));
      // binomial(pn, j) = 0 when j > pn, which has infinite valuation
      // For j <= pn, we compute the p-adic valuation
      let binoVal: number;
      if (bino === 0n) {
        binoVal = Number.POSITIVE_INFINITY;
      } else {
        binoVal = Number(valuation(bino, this._p));
      }
      enj = Math.min(binoVal, enj);
      res.push(enj);
    }

    return res;
  }

  /**
   * Compute the p-adic L-series as a power series.
   *
   * The p-adic L-series is computed as a power series in T, where T = gamma - 1
   * and gamma = 1 + p is a topological generator of 1 + pZ_p.
   *
   * INPUT:
   * - n: number of terms (default 2)
   * - quadratic_twist: a squarefree integer (default +1)
   * - prec: p-adic precision (default 5)
   * - eta: Teichmuller character twist (default 0)
   *
   * OUTPUT: a power series in T with coefficients in Q_p
   *
   * @see Reference: sage/schemes/elliptic_curves/padic_lseries.py:series
   */
  series(
    n: number = 2,
    quadratic_twist: bigint | number = 1,
    prec: number = 5,
    eta: number = 0
  ): PowerSeriesElement<pAdicGenericElement> {
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: pAdicLseries.series - requires p-adic fields and power series'
    );
  }

  /**
   * Return True if E has ordinary reduction at p.
   *
   * An elliptic curve has ordinary reduction at p if a_p is not divisible by p.
   * This is equivalent to the p-torsion E[p] having a non-trivial p-torsion point
   * over the algebraic closure of F_p.
   *
   * @see Reference: sage/schemes/elliptic_curves/padic_lseries.py:is_ordinary
   */
  is_ordinary(): boolean {
    // Check if E is ordinary at p by computing a_p mod p
    // This requires the ap() method on the elliptic curve
    try {
      if (
        'ap' in this._E &&
        typeof (this._E as unknown as { ap: (p: bigint) => bigint }).ap === 'function'
      ) {
        const a_p = (this._E as unknown as { ap: (p: bigint) => bigint }).ap(this._p);
        return a_p % this._p !== 0n;
      }
      // Also check is_ordinary method if available
      if (
        'is_ordinary' in this._E &&
        typeof (this._E as unknown as { is_ordinary: (p: bigint) => boolean }).is_ordinary ===
          'function'
      ) {
        return (this._E as unknown as { is_ordinary: (p: bigint) => boolean }).is_ordinary(this._p);
      }
    } catch {
      // Fall through to error
    }
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: pAdicLseries.is_ordinary - requires trace of Frobenius computation'
    );
  }

  /**
   * Return True if E has supersingular reduction at p.
   *
   * An elliptic curve has supersingular reduction at p if a_p is divisible by p.
   * This is the complement of ordinary reduction.
   *
   * @see Reference: sage/schemes/elliptic_curves/padic_lseries.py:is_supersingular
   */
  is_supersingular(): boolean {
    return !this.is_ordinary();
  }

  /**
   * Compute the Frobenius endomorphism on the formal group.
   *
   * Returns the matrix of Frobenius with respect to the basis {omega, eta}
   * of the Dieudonne module D_p(E) = H^1_dR(E/Q_p).
   *
   * INPUT:
   * - prec: precision (default 20)
   * - algorithm: 'mw' (Monsky-Washnitzer, default) or 'approx'
   *
   * OUTPUT: a 2x2 matrix over Q_p
   *
   * ALGORITHM:
   * The 'mw' algorithm uses Monsky-Washnitzer cohomology, which involves:
   * 1. Converting E to integral short Weierstrass form
   * 2. Computing the matrix of Frobenius on the MW cohomology
   * 3. Adjusting for the change of coordinates
   *
   * The 'approx' algorithm uses the Bernardi-Perrin-Riou approach,
   * which is slower but works for all primes.
   *
   * @see Reference: sage/schemes/elliptic_curves/padic_lseries.py:frobenius
   */
  frobenius(prec: number = 20, algorithm: 'mw' | 'approx' = 'mw'): FrobeniusMatrix {
    // Reference: sage/schemes/elliptic_curves/padic_lseries.py:frobenius
    // The Frobenius matrix phi satisfies phi^2 - (a_p/p)*phi + 1/p = 0
    //
    // For the 'mw' algorithm, we would use Monsky-Washnitzer cohomology
    // which is not yet implemented. For 'approx', we use Bernardi-Perrin-Riou.
    //
    // The matrix is 2x2 with entries in Q_p, representing the action of
    // Frobenius on the Dieudonne module D_p(E).

    const p = this._p;
    const E = this._E;

    if (algorithm !== 'mw' && algorithm !== 'approx') {
      throw new ValueError(`Unknown algorithm ${algorithm}.`);
    }

    // Get a_p
    let a_p: bigint;
    try {
      if ('ap' in E && typeof (E as unknown as { ap: (p: bigint) => bigint }).ap === 'function') {
        a_p = (E as unknown as { ap: (p: bigint) => bigint }).ap(p);
      } else {
        throw new Error('No ap method');
      }
    } catch {
      throw new NotImplementedError(
        'SAGE_NOT_IMPLEMENTED: frobenius - requires trace of Frobenius (ap) computation'
      );
    }

    // For the 'mw' algorithm, we need Monsky-Washnitzer cohomology
    // which is not available. We can only provide a partial result.
    if (algorithm === 'mw') {
      // The full MW algorithm requires:
      // 1. Convert E to integral short Weierstrass form
      // 2. Compute the Frobenius matrix using MW cohomology
      // 3. Change basis back to the original curve
      //
      // For now, throw an error indicating this is not implemented
      throw new NotImplementedError(
        'SAGE_NOT_IMPLEMENTED: frobenius with algorithm=mw - requires Monsky-Washnitzer cohomology'
      );
    }

    // For 'approx' algorithm (Bernardi-Perrin-Riou):
    // This uses the formal group to approximate Frobenius
    // Reference: sage/schemes/elliptic_curves/padic_lseries.py:__phi_bpr

    // The Frobenius satisfies phi^2 - (a_p/p)*phi + 1/p = 0
    // This means det(phi) = 1/p and tr(phi) = a_p/p
    //
    // For supersingular primes (a_p = 0), we have:
    // phi^2 + 1/p = 0, so phi = [[a, b], [c, d]] with a + d = 0 and ad - bc = 1/p
    //
    // For ordinary primes, the eigenvalues are alpha/p and beta/p
    // where alpha, beta are the roots of x^2 - a_p*x + p

    // The 'approx' algorithm is complex and slow
    // For now, provide a stub that shows the matrix structure
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: frobenius with algorithm=approx - requires formal group differential equation solver'
    );
  }

  /**
   * Compute the Bernardi sigma function.
   *
   * The Bernardi sigma function is a p-adic analogue of the Weierstrass sigma
   * function. It is used in the computation of p-adic heights and the p-adic
   * regulator.
   *
   * INPUT:
   * - prec: precision (default 20)
   *
   * OUTPUT: the sigma function as a power series in z = log(t), of absolute
   * precision ``prec + 5`` (exactly as in SageMath).
   *
   * ALGORITHM (verbatim from SageMath, padic_lseries.py:1626-1641):
   *
   *     Eh = E.formal()
   *     lo = Eh.log(prec + 5)
   *     F  = lo.reverse()(z)
   *     xofF = Eh.x(prec + 2)(F)
   *     g = (1/z^2 - xofF).power_series()
   *     h = g.integral().integral()
   *     sigma_of_z = z * h.exp()
   *
   * The only spelling difference is the Laurent bookkeeping: the port has no
   * Laurent-series arithmetic, so ``x(t) = t^-2 * u(t)`` is carried as the
   * power series ``u(t) = t^2 x(t)`` (whose coefficient list is exactly
   * ``Eh.x_list(prec+2)``) and the two Laurent divisions by ``z^2`` are done
   * by shifting.  Writing ``F = z*w`` with ``w(0) = 1``,
   *
   *     1/z^2 - x(F) = z^-2 * (1 - w^-2 * u(F))
   *
   * and the bracket has valuation >= 2, so the shift is exact.
   *
   * @see Reference: sage/schemes/elliptic_curves/padic_lseries.py:1613-1641
   */
  bernardi_sigma_function(prec: number = 20): PowerSeriesElement<RationalElement> {
    // Reference: sage/schemes/elliptic_curves/padic_lseries.py:1626
    //   E = self._E; Eh = E.formal()
    const E = this._E as unknown as {
      formal_group?: () => EllipticCurveFormalGroup;
      formal?: () => EllipticCurveFormalGroup;
    };
    const formal = E.formal_group ?? E.formal;
    if (typeof formal !== 'function') {
      throw new NotImplementedError(
        'SAGE_NOT_IMPLEMENTED: pAdicLseries.bernardi_sigma_function - curve does not have formal() method'
      );
    }
    const Eh = formal.call(E);

    if (prec < 1) {
      throw new ValueError('The precision must be positive.');
    }

    // Sage works in LaurentSeriesRing(QQ, 'z'); we do the same, over this
    // module's QQ.  Everything below is the image of Sage's computation.
    const QQr = new RationalRing();
    const S = new PowerSeriesRing<RationalElement>(QQr, 'z', prec + 5);

    // padic_lseries.py:1629  lo = Eh.log(prec + 5)
    const loCoeffs = Eh.log(prec + 5)
      .list()
      .map((c) => this._toRational(c));
    const lo = S.__call__(loCoeffs, prec + 5);

    // padic_lseries.py:1630/1634  F = lo.reverse()(z)
    const F = lo.reversion(prec + 5);

    // padic_lseries.py:1635  xofF = Eh.x(prec + 2)(F)
    // x(prec+2) has valuation -2 and absolute precision prec+2, so
    // u(t) = t^2*x(t) has absolute precision prec+4.
    const u = S.__call__(
      Eh.x_list(prec + 2).map((c) => this._toRational(c)),
      prec + 4
    );

    // 1/z^2 - x(F) = z^-2 * (1 - w^-2 * u(F))  with  w = F/z.
    const w = F._shiftRight(1);
    const bracket = S.one().sub(u.__call__(F).mul(w.pow(-2)));

    // The z^0 and z^1 coefficients of the bracket must vanish for the result
    // to be a power series (Sage relies on the same fact when it calls
    // ``.power_series()`` on g).
    for (let i = 0; i < 2; i++) {
      if (!bracket.__getitem__(i).isZero()) {
        throw new ValueError(
          `1/z^2 - x(F) is not a power series: z^${i - 2} coefficient is ${bracket.__getitem__(i)}`
        );
      }
    }

    // padic_lseries.py:1637-1639
    const g = bracket._shiftRight(2);
    const h = g.integral().integral();
    return S.gen().mul(h.exp());
  }

  /**
   * Coerce a coefficient of the formal group (an element of the curve's base
   * ring, which for `padic_lseries` is always QQ) into this module's
   * {@link RationalElement}.
   */
  private _toRational(c: unknown): RationalElement {
    if (c === undefined || c === null) return new RationalElement(0n, 1n);
    if (c instanceof RationalElement) return c;
    if (typeof c === 'bigint') return new RationalElement(c, 1n);
    if (typeof c === 'number') {
      if (!Number.isInteger(c)) {
        throw new ValueError(`cannot coerce ${c} to a rational number`);
      }
      return new RationalElement(BigInt(c), 1n);
    }
    if (typeof c === 'object') {
      const o = c as Record<string, unknown>;
      // sage/rings/rational.ts exposes numerator/denominator as getters;
      // other rings expose them as methods.
      for (const [n, d] of [
        ['numerator', 'denominator'],
        ['numer', 'denom'],
        ['num', 'den'],
      ] as const) {
        if (n in o && d in o) {
          const nv = typeof o[n] === 'function' ? (o[n] as () => unknown)() : o[n];
          const dv = typeof o[d] === 'function' ? (o[d] as () => unknown)() : o[d];
          if (
            (typeof nv === 'bigint' || typeof nv === 'number') &&
            (typeof dv === 'bigint' || typeof dv === 'number')
          ) {
            return new RationalElement(BigInt(nv), BigInt(dv));
          }
        }
      }
      // Last resort: the printed form "a" or "a/b".
      const s = String(c);
      const m = /^(-?\d+)(?:\/(\d+))?$/.exec(s);
      if (m) {
        return new RationalElement(BigInt(m[1]!), m[2] === undefined ? 1n : BigInt(m[2]));
      }
    }
    throw new ValueError(`cannot coerce ${c} to a rational number`);
  }

  /**
   * Helper to convert a power series element to have rational coefficients.
   */
  private _convertToRational(
    ps: PowerSeriesElement<RingElement>,
    targetRing: PowerSeriesRing<RationalElement>,
    prec: number
  ): PowerSeriesElement<RationalElement> {
    const QQ = targetRing.base_ring() as RationalRing;
    const coeffs: RationalElement[] = [];
    const psList = ps.list();
    for (let i = 0; i < Math.min(psList.length, prec); i++) {
      const c = psList[i];
      // Try to convert to rational
      if (c !== undefined) {
        const val = this._toBigInt(c);
        coeffs.push(QQ.__call__(val));
      } else {
        coeffs.push(QQ.zero());
      }
    }
    return targetRing.__call__(coeffs, prec);
  }

  /**
   * Helper to extract bigint from a ring element.
   */
  private _toBigInt(x: unknown): bigint {
    if (typeof x === 'bigint') return x;
    if (typeof x === 'number') return BigInt(x);
    if (x && typeof x === 'object') {
      if ('value' in x) return BigInt((x as { value: unknown }).value as bigint | number);
      if ('lift' in x && typeof (x as { lift: () => unknown }).lift === 'function') {
        return BigInt((x as { lift: () => unknown }).lift() as bigint | number);
      }
      if ('toString' in x) return BigInt((x as { toString: () => string }).toString());
    }
    return 0n;
  }

  /**
   * Get a cached series if available.
   * @see Reference: sage/schemes/elliptic_curves/padic_lseries.py:_get_series_from_cache
   */
  protected _get_series_from_cache(
    n: number,
    prec: number,
    D: bigint,
    eta: number
  ): PowerSeriesElement<pAdicGenericElement> | null {
    // Try exact match
    const key = `${n},${prec},${D},${eta}`;
    const cached = this.__series.get(key);
    if (cached) {
      return cached;
    }

    // Try to find a higher precision version
    for (const [cachedKey, cachedSeries] of this.__series) {
      const [_n, _prec, _D, _eta] = cachedKey
        .split(',')
        .map((x, i) => (i === 2 ? BigInt(x) : Number(x)));
      if (_n === n && _D === D && _eta === eta && Number(_prec) >= prec) {
        return cachedSeries.add_bigoh(prec);
      }
    }

    return null;
  }

  /**
   * Store a series in the cache.
   * @see Reference: sage/schemes/elliptic_curves/padic_lseries.py:_set_series_in_cache
   */
  protected _set_series_in_cache(
    n: number,
    prec: number,
    D: bigint,
    eta: number,
    f: PowerSeriesElement<pAdicGenericElement>
  ): void {
    const key = `${n},${prec},${D},${eta}`;
    this.__series.set(key, f);
  }

  /**
   * Compute the quotient of periods for a quadratic twist.
   *
   * For a fundamental discriminant D, computes the constant eta such that
   * sqrt(|D|) * Omega_{E_D}^+ = eta * Omega_E^{sign(D)}.
   *
   * According to [MTT1986] page 40, this is either 1 or 2 unless the
   * condition on the twist is not satisfied.
   *
   * NOTE: This computation requires period lattice computations which
   * are not yet fully implemented. For D = 1, returns 1. For other D,
   * throws NotImplementedError.
   *
   * @see Reference: sage/schemes/elliptic_curves/padic_lseries.py:_quotient_of_periods_to_twist
   */
  protected _quotient_of_periods_to_twist(D: bigint): bigint {
    // Reference: sage/schemes/elliptic_curves/padic_lseries.py:_quotient_of_periods_to_twist
    // This does not depend on p and could be moved elsewhere.
    //
    // Algorithm from SageMath:
    // 1. Compute the quadratic twist Et = E.quadratic_twist(D)
    // 2. For D > 1: qt = Et.period_lattice().basis()[0] / E.period_lattice().basis()[0]
    //    Then multiply by sqrt(D)
    // 3. For D < 0: qt = Et.period_lattice().basis()[1].imag() / E.period_lattice().basis()[0]
    //    Adjust for number of real components, multiply by sqrt(-D)
    // 4. Return QQ((8 * qt).round()) / 8

    if (D === 1n) {
      return 1n;
    }

    // The full computation requires:
    // - quadratic_twist() method on the curve
    // - period_lattice() computation
    // - real_components() method
    //
    // These are complex computations involving numerical integration
    // and lattice basis computation.

    // For many common cases, the result is 1 or 2
    // Check if the curve has the necessary methods
    const E = this._E;

    if ('quadratic_twist' in E && 'period_lattice' in E) {
      // Attempt the computation
      throw new NotImplementedError(
        'SAGE_NOT_IMPLEMENTED: _quotient_of_periods_to_twist - requires period lattice computation'
      );
    }

    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: _quotient_of_periods_to_twist - requires period lattice computation'
    );
  }

  /**
   * String representation.
   */
  toString(): string {
    let s = `${this._p}-adic L-series of ${this._E}`;
    if (this._normalize !== 'L_ratio') {
      s += ' (not normalized)';
    }
    return s;
  }
}

/**
 * p-adic L-series for ordinary primes.
 *
 * For ordinary primes, the p-adic L-function is an element of the Iwasawa algebra
 * Lambda(Z_p^*) with bounded coefficients.
 *
 * @see Reference: sage/schemes/elliptic_curves/padic_lseries.py:pAdicLseriesOrdinary
 * @see Deviation: Elliptic Curve p-adic L-series and Isogeny Class Partial Implementation
 */
export class pAdicLseriesOrdinary<F extends FieldElement = FieldElement> extends pAdicLseries<F> {
  /**
   * Create the p-adic L-series for an ordinary prime.
   *
   * @see Reference: sage/schemes/elliptic_curves/padic_lseries.py:pAdicLseriesOrdinary
   */
  constructor(E: EllipticCurveGeneric<F>, p: bigint | number, options: pAdicLseriesOptions = {}) {
    super(E, p, options);
  }

  /**
   * Return True - this class is for ordinary primes.
   */
  override is_ordinary(): boolean {
    return true;
  }

  /**
   * Return False - this class is for ordinary primes.
   */
  override is_supersingular(): boolean {
    return false;
  }

  /**
   * Compute the p-adic L-series as a power series.
   *
   * For ordinary primes, each coefficient is a p-adic number whose precision
   * is provably correct.
   *
   * The normalization is chosen such that:
   *   L_p(E, 1) = (1 - 1/alpha)^2 * L(E, 1) / Omega_E
   *
   * INPUT:
   * - n: a positive integer (number of p^n summands to use)
   * - quadratic_twist: a fundamental discriminant of a quadratic field (default +1)
   * - prec: maximal number of terms of the series (default 5)
   * - eta: power of Teichmuller character (default 0)
   *
   * ALGORITHM:
   * The series is computed as:
   *   L_p(E, T) = sum_{j=0}^{p^{n-1}-1} (sum_{a=1}^{p-1} tau(a)^eta * mu(tau(a) * gamma^j))
   *              * (1+T)^j
   *
   * where gamma = 1 + p, tau is the Teichmuller character, and mu is the measure.
   *
   * @see Reference: sage/schemes/elliptic_curves/padic_lseries.py:pAdicLseriesOrdinary.series
   */
  override series(
    n: number = 2,
    quadratic_twist: bigint | number = 1,
    prec: number = 5,
    eta: number = 0
  ): PowerSeriesElement<pAdicGenericElement> {
    // Input validation
    if (n < 1) {
      throw new ValueError(`n (=${n}) must be a positive integer`);
    }
    if (this._p === 2n && n === 1) {
      throw new ValueError(`n (=${n}) must be at least 2 if p is 2`);
    }
    if (prec < 1) {
      throw new ValueError(`Insufficient precision (${prec})`);
    }

    // eta only matters modulo p-1 (modulo 2 for p = 2): padic_lseries.py:868
    const modulus = this._p !== 2n ? Number(this._p - 1n) : 2;
    eta = ((eta % modulus) + modulus) % modulus;

    const D = BigInt(quadratic_twist);
    if (D !== 1n && eta !== 0) {
      throw new NotImplementedError(
        'quadratic twists only implemented for the 0th Teichmueller component'
      );
    }

    // Validate quadratic twist is a fundamental discriminant
    if (D !== 1n) {
      if (!isFundamentalDiscriminant(D)) {
        throw new ValueError(
          `quadratic_twist (=${D}) must be a fundamental discriminant of a quadratic field`
        );
      }
      if (gcd(D < 0n ? -D : D, this._p) !== 1n) {
        throw new ValueError(`quadratic twist (=${D}) must be coprime to p (=${this._p})`);
      }
    }

    // Full implementation requires modular symbols and p-adic arithmetic
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: pAdicLseriesOrdinary.series - requires modular symbols and p-adic arithmetic'
    );
  }

  /**
   * Compute bounds on coefficient precision.
   *
   * @see Reference: sage/schemes/elliptic_curves/padic_lseries.py:_prec_bounds
   */
  protected _prec_bounds(n: number, prec: number, sign: 1 | -1 = 1): number[] {
    // Reference: sage/schemes/elliptic_curves/padic_lseries.py:_prec_bounds
    // Uses _e_bounds and _c_bound

    let e: number[];
    if (this._p === 2n) {
      e = this._e_bounds(n - 2, prec);
    } else {
      e = this._e_bounds(n - 1, prec);
    }

    const c = this._c_bound(sign);
    return e.map((ej) => (ej === Number.POSITIVE_INFINITY ? Number.POSITIVE_INFINITY : ej - c));
  }

  /**
   * Compute upper bound on denominators of modular symbols.
   *
   * @see Reference: sage/schemes/elliptic_curves/padic_lseries.py:_c_bound
   */
  protected _c_bound(_sign: 1 | -1 = 1): number {
    // Sage first tests `E.galois_representation().is_irreducible(p)` (returning
    // 0 in that case) and otherwise bounds the p-adic valuation of the modular
    // symbol denominators using the modular symbol space or the X_0-optimal
    // curve of the isogeny class.  Neither `gal_reps.py` nor the modular symbol
    // machinery is ported, so there is no justified value to return here.
    //
    // Returning 0 unconditionally is *not* safe: c is subtracted from the
    // e-bounds in `_prec_bounds`, so a too-small c over-reports precision.
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: pAdicLseries._c_bound - requires E.galois_representation() and modular symbol denominators'
    );
  }

  /**
   * Alias for series().
   */
  power_series = this.series;
}

/**
 * p-adic L-series for supersingular primes.
 *
 * For supersingular primes, the p-adic L-function has coefficients in a quadratic
 * extension of Q_p, and the coefficients are unbounded.
 *
 * @see Reference: sage/schemes/elliptic_curves/padic_lseries.py:pAdicLseriesSupersingular
 * @see Deviation: Elliptic Curve p-adic L-series and Isogeny Class Partial Implementation
 */
export class pAdicLseriesSupersingular<
  F extends FieldElement = FieldElement,
> extends pAdicLseries<F> {
  /**
   * Create the p-adic L-series for a supersingular prime.
   *
   * @see Reference: sage/schemes/elliptic_curves/padic_lseries.py:pAdicLseriesSupersingular
   */
  constructor(E: EllipticCurveGeneric<F>, p: bigint | number, options: pAdicLseriesOptions = {}) {
    super(E, p, options);
  }

  /**
   * Return False - this class is for supersingular primes.
   */
  override is_ordinary(): boolean {
    return false;
  }

  /**
   * Return True - this class is for supersingular primes.
   */
  override is_supersingular(): boolean {
    return true;
  }

  /**
   * Compute the p-adic L-series as a power series.
   *
   * For supersingular primes, each coefficient is an element of a quadratic
   * extension of Q_p, and the coefficients are unbounded.
   *
   * INPUT:
   * - n: a positive integer (default 3)
   * - quadratic_twist: a fundamental discriminant (default +1)
   * - prec: maximal number of terms (default 5)
   * - eta: power of Teichmuller character (default 0)
   *
   * OUTPUT:
   * A power series with coefficients in a quadratic ramified extension of
   * the p-adic numbers generated by a root alpha of the characteristic
   * polynomial of Frobenius on T_pE.
   *
   * @see Reference: sage/schemes/elliptic_curves/padic_lseries.py:pAdicLseriesSupersingular.series
   */
  override series(
    n: number = 3,
    quadratic_twist: bigint | number = 1,
    prec: number = 5,
    eta: number = 0
  ): PowerSeriesElement<pAdicGenericElement> {
    // Input validation
    if (n < 1) {
      throw new ValueError(`n (=${n}) must be a positive integer`);
    }
    if (this._p === 2n && n === 1) {
      throw new ValueError(`n (=${n}) must be at least 2 when p=2`);
    }
    if (prec < 1) {
      throw new ValueError(`Insufficient precision (${prec})`);
    }

    const D = BigInt(quadratic_twist);
    if (D !== 1n && eta !== 0) {
      throw new NotImplementedError(
        'quadratic twists only implemented for the 0th Teichmueller component'
      );
    }

    // Validate quadratic twist is a fundamental discriminant
    if (D !== 1n) {
      if (!isFundamentalDiscriminant(D)) {
        throw new ValueError(
          `quadratic_twist (=${D}) must be a fundamental discriminant of a quadratic field`
        );
      }
    }

    // Full implementation requires modular symbols and p-adic arithmetic
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: pAdicLseriesSupersingular.series - requires modular symbols and p-adic arithmetic'
    );
  }

  /**
   * Compute bounds on coefficient precision (alpha-adic).
   *
   * @see Reference: sage/schemes/elliptic_curves/padic_lseries.py:pAdicLseriesSupersingular._prec_bounds
   */
  protected _prec_bounds(n: number, prec: number): number[] {
    // Reference: sage/schemes/elliptic_curves/padic_lseries.py
    // For supersingular, the bounds are alpha-adic (not p-adic)

    let e: number[];
    if (this._p === 2n) {
      e = this._e_bounds(n - 2, prec);
    } else {
      e = this._e_bounds(n - 1, prec);
    }

    const c0 = n + 2;
    const result: number[] = [Number.POSITIVE_INFINITY];
    for (let j = 1; j < e.length; j++) {
      result.push(2 * e[j]! - c0);
    }
    return result;
  }

  /**
   * Extract polynomial representation of a Qp[alpha] element.
   *
   * Given an element a in Qp[alpha], returns [v0, v1] such that a = v0 + v1*alpha.
   * Here alpha is a root of x^2 - a_p*x + p, where a_p = 0 for supersingular primes.
   *
   * INPUT:
   * - a: an element in the quadratic extension Q_p[alpha]
   *
   * OUTPUT: a tuple [v0, v1] of elements in Q_p such that a = v0 + v1*alpha
   *
   * ALGORITHM:
   * For elements of Eisenstein extensions (supersingular case), we use
   * the internal NTL representation to extract the polynomial coefficients.
   *
   * @see Reference: sage/schemes/elliptic_curves/padic_lseries.py:_poly
   */
  protected _poly(a: unknown): [pAdicGenericElement, pAdicGenericElement] {
    // Reference: sage/schemes/elliptic_curves/padic_lseries.py:_poly
    // The SageMath implementation uses:
    // v, k = a._ntl_rep_abs()
    // K = a.base_ring()
    // pi = K.uniformiser()
    // v0 = K(v[0]._sage_()) * pi**k
    // v1 = K(v[1]._sage_()) * pi**k
    // return [v0, v1]
    //
    // This requires the element to be from an Eisenstein extension
    // of Q_p, which we don't have full support for.

    // Check if a is already a base p-adic element (not in extension)
    if (a instanceof pAdicGenericElement) {
      // Element is in Q_p, so a = a + 0*alpha
      const K = Qp(this._p, 20);
      return [a, K.zero()];
    }

    // Check if a is zero
    if (
      a === 0 ||
      a === 0n ||
      (a && typeof a === 'object' && 'is_zero' in a && (a as { is_zero: () => boolean }).is_zero())
    ) {
      const K = Qp(this._p, 20);
      return [K.zero(), K.zero()];
    }

    // For proper Eisenstein extension elements, we need the _ntl_rep_abs method
    // This is specific to SageMath's p-adic extension implementation
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: _poly - requires p-adic Eisenstein extension elements with _ntl_rep_abs method'
    );
  }

  /**
   * Return a vector of two p-adic power series for the D_p valued L-series.
   *
   * The result v satisfies:
   *   (1 - phi)^{-2} * L_p(E, T) = v[0] * omega + v[1] * phi(omega)
   *
   * where omega is the invariant differential and phi is Frobenius.
   *
   * INPUT:
   * - n: positive integer for approximation level (default 3)
   * - quadratic_twist: fundamental discriminant (default +1)
   * - prec: number of terms (default 5)
   *
   * OUTPUT: a tuple (G, H) of two power series such that
   *   L_p(E) = G + H * alpha
   * where alpha is a root of the Frobenius characteristic polynomial.
   *
   * ALGORITHM:
   * 1. Compute the p-adic L-series via series()
   * 2. Split into components using the representation in Q_p[alpha]
   * 3. Apply the Frobenius transformation
   *
   * @see Reference: sage/schemes/elliptic_curves/padic_lseries.py:Dp_valued_series
   */
  Dp_valued_series(
    n: number = 3,
    quadratic_twist: bigint | number = 1,
    prec: number = 5
  ): [PowerSeriesElement<pAdicGenericElement>, PowerSeriesElement<pAdicGenericElement>] {
    // Reference: sage/schemes/elliptic_curves/padic_lseries.py:Dp_valued_series
    // Algorithm:
    // 1. lps = self.series(n, quadratic_twist=quadratic_twist, prec=prec)
    // 2. Split lps into G + H*alpha where G, H have coefficients in Q_p
    // 3. Compute phi = [[0, -1/p], [1, a_p/p]]
    // 4. lpv = [G + a_p*H, -p*H]
    // 5. eps = (1 - phi)^{-2}
    // 6. return lpv * eps.transpose()

    // This requires the series computation which depends on modular symbols
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: Dp_valued_series - requires series() which needs modular symbols'
    );
  }

  /**
   * Compute the D_p valued p-adic height.
   *
   * Returns the canonical p-adic height with values in the Dieudonne module D_p(E).
   * The height is defined as:
   *   h_eta * omega - h_omega * eta
   * where h_eta uses the Bernardi sigma function and h_omega = log_E^2.
   *
   * INPUT:
   * - prec: precision (default 20)
   *
   * OUTPUT: a function that takes a point P on E and returns a vector (v1, v2)
   * such that the D_p-valued height of P is v1*omega + v2*eta.
   *
   * ALGORITHM:
   * 1. Get the formal group of E
   * 2. Compute the formal log and Bernardi sigma function
   * 3. For a point P, compute z = log(t) where t is the local parameter
   * 4. h_omega = -z^2 / n^2 where n is a multiple making P reduce well
   * 5. h_eta = 2 * log(sigma(z)/e_Q) / n^2
   * 6. Return the vector [-h_eta, h_omega]
   *
   * @see Reference: sage/schemes/elliptic_curves/padic_lseries.py:Dp_valued_height
   */
  Dp_valued_height(prec: number = 20): (P: unknown) => [pAdicGenericElement, pAdicGenericElement] {
    // Reference: sage/schemes/elliptic_curves/padic_lseries.py:Dp_valued_height
    // This is a complex computation that requires:
    // 1. formal() method on the curve to get the formal group
    // 2. formal_log computation
    // 3. bernardi_sigma_function (which we've partially implemented)
    // 4. A way to compute n = multiple to make good reduction
    // 5. Point arithmetic on the curve

    const E = this._E;
    const p = this._p;

    // Check if the curve has the necessary infrastructure
    if (!('formal' in E)) {
      throw new NotImplementedError(
        'SAGE_NOT_IMPLEMENTED: Dp_valued_height - curve does not have formal() method'
      );
    }

    // The full implementation would return a function like:
    // (P) => {
    //   if (P.is_finite_order()) return [Qp.zero(), Qp.zero()];
    //   const Q = n * P;
    //   const t = -Q.x / Q.y;
    //   const z = elog(t);
    //   const h_omega = -z^2 / n^2;
    //   const sigma = this.bernardi_sigma_function(prec);
    //   const h_eta = 2 * log(sigma(z) / e_Q) / n^2;
    //   return [-h_eta, h_omega];
    // }

    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: Dp_valued_height - requires formal group and point arithmetic'
    );
  }

  /**
   * Compute the D_p valued p-adic regulator.
   *
   * Returns the canonical p-adic regulator with values in the Dieudonne module D_p(E)
   * as defined by Perrin-Riou.
   *
   * The result is written in the basis omega, phi(omega), and hence the
   * coordinates are independent of the chosen Weierstrass equation.
   *
   * INPUT:
   * - prec: precision (default 20)
   * - v1, v2: optional vectors for the computation (default: standard choice)
   *
   * OUTPUT: a vector (r1, r2) representing the regulator in the basis omega, phi(omega)
   *
   * ALGORITHM:
   * 1. Compute the D_p-valued height function h
   * 2. For each generator g_i of E(Q), compute h(g_i)
   * 3. Form the height pairing matrix
   * 4. Compute determinant and transform to the omega, phi(omega) basis
   *
   * NOTE: The definition here is corrected with respect to Perrin-Riou's article.
   * See [SW2013] for details.
   *
   * @see Reference: sage/schemes/elliptic_curves/padic_lseries.py:Dp_valued_regulator
   */
  Dp_valued_regulator(prec: number = 20): [pAdicGenericElement, pAdicGenericElement] {
    // Reference: sage/schemes/elliptic_curves/padic_lseries.py:Dp_valued_regulator
    // This requires:
    // 1. Dp_valued_height (which we've stubbed)
    // 2. The rank of E(Q) and its generators
    // 3. Frobenius matrix computation
    // 4. Matrix operations over Q_p

    // The algorithm from SageMath:
    // 1. Get the height function h = self.Dp_valued_height(prec)
    // 2. Define hv(vec, P) = -vec[0]*h(P)[1] + vec[1]*h(P)[0]
    // 3. Choose v1 = [0, 1] (eta) and v2 = [-1, 1] (eta - omega)
    // 4. For each pair of generators, compute the height pairing
    // 5. Compute regv(vec) = det of the height pairing matrix with vec
    // 6. reg_oe = (reg1 * v2 - reg2 * v1) / Dp_pairing(v2, v1)
    // 7. Transform using Frobenius to get coordinates in omega, phi(omega) basis

    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: Dp_valued_regulator - requires height computation, rank, and Frobenius'
    );
  }

  /**
   * Alias for series().
   */
  power_series = this.series;
}
