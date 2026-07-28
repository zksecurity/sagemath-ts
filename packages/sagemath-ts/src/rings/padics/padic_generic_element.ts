/**
 * @module sage/rings/padics/padic_generic_element
 * @description Elements of p-adic rings and fields
 *
 * Port of: sage/rings/padics/padic_generic_element.pyx
 * Reference: reference/sage/src/sage/rings/padics/padic_generic_element.pyx
 */

import { Z_factor } from '@sagemath-ts/parigp-ts';
import {
  NotImplementedError,
  PrecisionError,
  ValueError,
  ZeroDivisionError,
} from '../../errors.js';
import { Integer } from '../integer_ring.js';
import type { pAdicGeneric } from './padic_generic.js';

// Re-export PrecisionError for backwards compatibility
export { PrecisionError } from '../../errors.js';

/**
 * Sentinel used where SageMath returns ``+Infinity`` (valuations of exact zero,
 * infinite multiplicative/additive orders).
 *
 * SageMath's `sage.rings.infinity.PlusInfinity` has no bigint analogue, so we
 * use the IEEE `+Infinity` double.  Relational comparisons against `bigint`
 * behave as expected in JavaScript (`Infinity > 0n` is `true`).
 */
export type InfiniteOr<T> = T | number;

/**
 * Compute the p-adic valuation of n with respect to prime p.
 * Returns the largest k such that p^k divides n.
 */
function padic_valuation(n: bigint, p: bigint): bigint {
  if (n === 0n) {
    throw new ValueError('valuation of zero is infinity');
  }
  let v = 0n;
  let absN = n < 0n ? -n : n;
  while (absN % p === 0n) {
    absN = absN / p;
    v++;
  }
  return v;
}

/**
 * Compute p^n for bigint p and bigint n.
 */
function bigPow(p: bigint, n: bigint): bigint {
  if (n < 0n) {
    throw new ValueError('negative exponent');
  }
  let result = 1n;
  let base = p;
  let exp = n;
  while (exp > 0n) {
    if (exp % 2n === 1n) {
      result *= base;
    }
    base *= base;
    exp = exp / 2n;
  }
  return result;
}

/**
 * Compute the modular inverse of a mod m using extended Euclidean algorithm.
 */
function modInverse(a: bigint, m: bigint): bigint {
  let [oldR, r] = [a % m, m];
  let [oldS, s] = [1n, 0n];

  while (r !== 0n) {
    const q = oldR / r;
    [oldR, r] = [r, oldR - q * r];
    [oldS, s] = [s, oldS - q * s];
  }

  if (oldR > 1n) {
    throw new ValueError('modular inverse does not exist');
  }
  if (oldS < 0n) {
    oldS += m;
  }
  return oldS;
}

/**
 * Positive mod - always returns a value in [0, m)
 */
function posMod(a: bigint, m: bigint): bigint {
  const r = a % m;
  return r < 0n ? r + m : r;
}

/**
 * Modular exponentiation.
 */
function modPow(base: bigint, exp: bigint, mod: bigint): bigint {
  let result = 1n;
  let b = posMod(base, mod);
  let e = exp;
  while (e > 0n) {
    if (e % 2n === 1n) {
      result = (result * b) % mod;
    }
    e = e / 2n;
    b = (b * b) % mod;
  }
  return result;
}

/**
 * Return the multiplicative order of ``a`` in the residue field GF(p).
 *
 * This is the port of `sage.rings.finite_rings.integer_mod.IntegerMod.multiplicative_order`
 * specialised to a prime field: the group is cyclic of order `p - 1`, so the
 * order is obtained by removing prime factors of `p - 1` one at a time.
 * Factorisation is delegated to PARI, as SageMath does.
 *
 * @see Reference: sage/rings/finite_rings/integer_mod.pyx:multiplicative_order
 */
function multiplicativeOrderModP(a: bigint, p: bigint): bigint {
  const res = posMod(a, p);
  if (res === 0n) {
    throw new ValueError('multiplicative order of 0 is not defined');
  }
  let order = p - 1n;
  if (p === 2n) {
    return 1n;
  }
  for (const [q, _e] of Z_factor(p - 1n)) {
    if (q < 0n) continue;
    while (order % q === 0n && modPow(res, order / q, p) === 1n) {
      order = order / q;
    }
  }
  return order;
}

/**
 * An element of a p-adic ring or field.
 *
 * This implementation uses "capped-relative" precision, which stores:
 * - The valuation (ord_p of the element)
 * - The unit part (the element divided by p^valuation)
 * - The relative precision (number of significant p-adic digits)
 *
 * An element x is represented as: x = p^valuation * unit_part + O(p^(valuation + relprec))
 *
 * @see Reference: sage/rings/padics/padic_generic_element.pyx:pAdicGenericElement
 */
export class pAdicGenericElement {
  protected readonly _parent: pAdicGeneric;
  /** The valuation (ord_p) of this element */
  protected _valuation: bigint;
  /** The unit part (reduced mod p^precision) */
  protected _unit: bigint;
  /** The relative precision */
  protected _relprec: number;
  /** Whether this is an exact zero */
  protected _exactZero: boolean;

  constructor(parent: pAdicGeneric, value: bigint = 0n, absprec?: number, relprec?: number) {
    this._parent = parent;
    this._exactZero = false;
    this._relprec = relprec ?? parent.precision_cap();

    if (value === 0n) {
      // Zero element.  An inexact zero O(p^n) is stored as valuation n with
      // relative precision 0, so that precision_absolute() == n.
      this._valuation = BigInt(absprec ?? parent.precision_cap());
      this._unit = 0n;
      this._exactZero = absprec === undefined && relprec === undefined;
      if (!this._exactZero) {
        this._relprec = 0;
      }
    } else {
      const p = parent.prime();
      // Compute valuation
      this._valuation = padic_valuation(value, p);
      // Extract unit part
      const unit = value / bigPow(p, this._valuation);
      // Reduce unit mod p^relprec
      const prec = this._relprec;
      const modulus = bigPow(p, BigInt(prec));
      this._unit = posMod(unit, modulus);
    }
  }

  /**
   * Create a p-adic element from a valuation and unit part.
   */
  static fromValuationUnit(
    parent: pAdicGeneric,
    valuation: bigint,
    unit: bigint,
    relprec: number
  ): pAdicGenericElement {
    const elem = new pAdicGenericElement(parent, 0n);
    elem._valuation = valuation;
    elem._unit = unit;
    elem._relprec = relprec;
    elem._exactZero = false;
    return elem;
  }

  /**
   * Create an exact zero element.
   */
  static exactZero(parent: pAdicGeneric): pAdicGenericElement {
    const elem = new pAdicGenericElement(parent, 0n);
    elem._exactZero = true;
    return elem;
  }

  /**
   * Return the parent ring or field.
   * @see Reference: sage/rings/padics/padic_generic_element.pyx:parent
   */
  parent(): pAdicGeneric {
    return this._parent;
  }

  /**
   * Return the prime p.
   * @see Reference: sage/rings/padics/padic_generic_element.pyx:prime
   */
  prime(): bigint {
    return this._parent.prime();
  }

  /**
   * Return the precision of this element.
   * @see Reference: sage/rings/padics/padic_generic_element.pyx:precision_absolute
   */
  precision_absolute(): number {
    if (this._exactZero) {
      return Number.POSITIVE_INFINITY;
    }
    return Number(this._valuation) + this._relprec;
  }

  /**
   * Return the relative precision.
   * @see Reference: sage/rings/padics/padic_generic_element.pyx:precision_relative
   */
  precision_relative(): number {
    if (this._exactZero || this._unit === 0n) {
      return 0;
    }
    return this._relprec;
  }

  /**
   * Return the valuation of this element.
   * @see Reference: sage/rings/padics/padic_generic_element.pyx:valuation
   */
  valuation(): InfiniteOr<bigint> {
    if (this._exactZero) {
      // SageMath returns +Infinity for the valuation of an exact zero.
      return Number.POSITIVE_INFINITY;
    }
    if (this._unit === 0n) {
      // Inexact zero - valuation is the absolute precision
      return BigInt(this.precision_absolute());
    }
    return this._valuation;
  }

  /**
   * Alias for valuation.
   * @see Reference: sage/rings/padics/padic_generic_element.pyx:ordp
   */
  ordp(): InfiniteOr<bigint> {
    return this.valuation();
  }

  /**
   * Return the unit part of this element.
   * @see Reference: sage/rings/padics/padic_generic_element.pyx:unit_part
   */
  unit_part(): pAdicGenericElement {
    if (this.is_zero()) {
      throw new ValueError('unit part of zero is not defined');
    }
    return pAdicGenericElement.fromValuationUnit(this._parent, 0n, this._unit, this._relprec);
  }

  /**
   * Return the normalized valuation.
   * This is the valuation divided by the ramification index.
   * @see Reference: sage/rings/padics/padic_generic_element.pyx:normalized_valuation
   */
  normalized_valuation(): InfiniteOr<bigint> {
    const v = this.valuation();
    if (typeof v === 'number') {
      return v;
    }
    const e = this._parent.ramification_index();
    return v / e;
  }

  /**
   * Return the expansion as a list of coefficients.
   *
   * If this is a field element the expansion starts at `p^valuation` and has
   * `precision_relative()` entries; for a ring element it starts at `p^0` and
   * has `precision_absolute()` entries (so the first `valuation` entries are
   * zero).  This mirrors SageMath's `shift` bookkeeping.
   *
   * @see Reference: sage/rings/padics/padic_template_element.pxi:expansion
   */
  expansion(): bigint[] {
    const prec = this.precision_relative();
    if (prec === 0) {
      // Exact and inexact zeros have empty expansions in SageMath.
      return [];
    }
    const p = this.prime();
    const result: bigint[] = [];

    // For rings SageMath shifts the expansion by the valuation; for fields it
    // does not (the expansion is read relative to p^valuation).
    const shift = this._parent.is_field() ? 0 : Number(this._valuation);
    for (let i = 0; i < shift; i++) {
      result.push(0n);
    }

    // Extract digits from unit part
    let unit = this._unit;
    for (let i = 0; i < prec; i++) {
      result.push(unit % p);
      unit = unit / p;
    }
    return result;
  }

  /**
   * Return the list of coefficients (alias for expansion).
   * @see Reference: sage/rings/padics/padic_generic_element.pyx:list
   */
  list(): bigint[] {
    return this.expansion();
  }

  /**
   * Return the i-th coefficient.
   * @see Reference: sage/rings/padics/padic_generic_element.pyx:__getitem__
   */
  __getitem__(i: number): bigint {
    if (i < 0) {
      // For fields, negative indices are allowed
      if (!this._parent.is_field()) {
        throw new ValueError('negative index not allowed for ring elements');
      }
    }
    if (i >= this.precision_absolute()) {
      throw new PrecisionError('coefficient beyond precision');
    }

    const val = Number(this._valuation);
    if (i < val) {
      return 0n;
    }

    const p = this.prime();
    const shift = BigInt(i - val);
    const unit = this._unit / bigPow(p, shift);
    return unit % p;
  }

  /**
   * Return the sum of the `p^(i + l*k)` terms of the series expansion of this
   * element, for `i + l*k` between `i` and `j-1` inclusive.
   *
   * `i === null` starts at the valuation, `j === null` (or `+Infinity`) stops
   * at the absolute precision.
   *
   * @see Reference: sage/rings/padics/local_generic_element.pyx:slice
   */
  slice(
    i: number | null = null,
    j: number | null = null,
    k: number | null = 1
  ): pAdicGenericElement {
    if (k === null) {
      k = 1;
    }
    if (k <= 0) {
      throw new ValueError('slice step must be positive');
    }
    let iStart: number;
    if (i === null) {
      iStart = Number(this.valuation());
    } else {
      iStart = i;
    }
    let jStop: number;
    if (j === null || j === Number.POSITIVE_INFINITY) {
      jStop = this.precision_absolute();
      if (jStop === Number.POSITIVE_INFINITY) {
        return this._parent.zero();
      }
    } else {
      jStop = j;
    }

    let start = iStart;
    let stop = jStop;

    // For fields, expansion() contains only the coefficients starting from the
    // valuation, so shift the indices to make up for this.
    if (this._parent.is_field()) {
      const v = Number(this.valuation());
      start -= v;
      stop -= v;
    }

    // Make sure start and stop are nonnegative
    if (start < 0) {
      iStart += -start; // fixes the p-power of the first kept term
      start = 0;
    }
    stop = Math.max(stop, 0);

    const p = this.prime();
    const cap = this._parent.precision_cap();
    const digits = this.expansion();

    let ans = this._parent.zero();
    let ppow = BigInt(iStart);
    for (let idx = start; idx < stop && idx < digits.length; idx += k) {
      const d = digits[idx]!;
      if (d !== 0n) {
        ans = ans.add(pAdicGenericElement.fromValuationUnit(this._parent, ppow, d, cap));
      }
      ppow += BigInt(k);
    }

    // Fix the precision of the return value
    const selfPrec = this.precision_absolute();
    const ansPrec = ans.precision_absolute();
    if (jStop < ansPrec || selfPrec < ansPrec) {
      ans = ans.add_bigoh(Math.min(jStop, selfPrec));
    }

    return ans;
  }

  /**
   * Return the residue modulo p^n.
   * @see Reference: sage/rings/padics/padic_generic_element.pyx:residue
   */
  residue(n: number = 1): bigint {
    if (this.valuation() < 0n) {
      throw new ValueError('element must have nonnegative valuation in order to compute residue');
    }
    if (n > this.precision_absolute()) {
      throw new PrecisionError('not enough precision');
    }

    const p = this.prime();
    const modulus = bigPow(p, BigInt(n));
    const val = Number(this._valuation);

    if (val >= n) {
      return 0n;
    }

    return (this._unit * bigPow(p, BigInt(val))) % modulus;
  }

  /**
   * Lift to an integer.
   * @see Reference: sage/rings/padics/padic_generic_element.pyx:lift
   */
  lift(): bigint {
    if (this._exactZero) {
      return 0n;
    }
    const p = this.prime();
    return this._unit * bigPow(p, this._valuation);
  }

  /**
   * Lift to the integers mod p^n.
   * @see Reference: sage/rings/padics/padic_generic_element.pyx:lift_to_precision
   */
  lift_to_precision(n: number): pAdicGenericElement {
    if (this._exactZero) {
      return pAdicGenericElement.exactZero(this._parent);
    }

    const currentPrec = this.precision_absolute();
    if (n <= currentPrec) {
      return this;
    }

    // Create a new element with higher precision
    // The unit part stays the same, but relprec increases
    const newRelprec = n - Number(this._valuation);
    return pAdicGenericElement.fromValuationUnit(
      this._parent,
      this._valuation,
      this._unit,
      newRelprec
    );
  }

  /**
   * Check if this element is zero.
   * @see Reference: sage/rings/padics/padic_generic_element.pyx:is_zero
   */
  is_zero(): boolean {
    return this._exactZero || this._unit === 0n;
  }

  /**
   * Check if this element is one.
   * @see Reference: sage/rings/padics/padic_generic_element.pyx:is_one
   */
  is_one(): boolean {
    if (this._exactZero || this._unit === 0n) {
      return false;
    }
    return this._valuation === 0n && this._unit === 1n;
  }

  /**
   * Check if this element is a unit.
   * @see Reference: sage/rings/padics/padic_generic_element.pyx:is_unit
   */
  is_unit(): boolean {
    if (this.is_zero()) {
      return false;
    }
    return this._valuation === 0n;
  }

  /**
   * Check if this element is a p-adic unit (valuation is zero).
   */
  is_padic_unit(): boolean {
    return this.is_unit();
  }

  /**
   * Check if this element is integral.
   * @see Reference: sage/rings/padics/padic_generic_element.pyx:is_integral
   */
  is_integral(): boolean {
    if (this.is_zero()) {
      return true;
    }
    return this._valuation >= 0n;
  }

  /**
   * Check if this element is a square.
   * @see Reference: sage/rings/padics/padic_generic_element.pyx:is_square
   */
  is_square(): boolean {
    if (this._exactZero) {
      return true;
    }
    const parent = this._parent;
    const p = this.prime();

    if (p !== 2n) {
      if (this.is_zero()) {
        throw new PrecisionError(
          'not enough precision to be sure that this element has a square root'
        );
      }
      if (this._valuation % 2n !== 0n) {
        return false;
      }
      // Euler criterion in the residue field: a is a QR mod p iff a^((p-1)/2) = 1
      const residue = this.unit_part().residue(1);
      return modPow(residue, (p - 1n) / 2n, p) === 1n;
    }

    // p == 2: SageMath simply attempts the square root at precision
    // valuation + 2*e + 1.
    const e = parent.absolute_e();
    try {
      this.add_bigoh(Number(this._valuation) + 2 * e + 1)._nth_root(2n);
    } catch (err) {
      if (err instanceof PrecisionError) {
        throw new PrecisionError(
          'not enough precision to be sure that this element has a square root'
        );
      }
      if (err instanceof ValueError) {
        return false;
      }
      throw err;
    }
    return true;
  }

  /**
   * Return the square root of this p-adic number.
   *
   * @param options.extend - if true (default) and no root exists in the parent,
   *   raise NotImplementedError (SageMath would move to an extension); if false
   *   raise a ValueError.
   * @param options.all - if true, return the list of all square roots.
   *
   * @see Reference: sage/rings/padics/padic_generic_element.pyx:square_root
   */
  square_root(options?: { extend?: boolean; all?: false }): pAdicGenericElement;
  square_root(options: { extend?: boolean; all: true }): pAdicGenericElement[];
  square_root(options?: {
    extend?: boolean;
    all?: boolean;
  }): pAdicGenericElement | pAdicGenericElement[];
  square_root(options?: {
    extend?: boolean;
    all?: boolean;
  }): pAdicGenericElement | pAdicGenericElement[] {
    return this._square_root(options);
  }

  /**
   * Alias for {@link square_root}.
   * @see Reference: sage/rings/padics/local_generic_element.pyx:sqrt
   */
  sqrt(options?: { extend?: boolean; all?: false }): pAdicGenericElement;
  sqrt(options: { extend?: boolean; all: true }): pAdicGenericElement[];
  sqrt(options?: {
    extend?: boolean;
    all?: boolean;
  }): pAdicGenericElement | pAdicGenericElement[];
  sqrt(options?: {
    extend?: boolean;
    all?: boolean;
  }): pAdicGenericElement | pAdicGenericElement[] {
    return this._square_root(options);
  }

  /**
   * Return all square roots of this element (possibly the empty list).
   * @see Reference: sage/rings/padics/padic_generic_element.pyx:square_root
   */
  square_root_all(options?: { extend?: boolean }): pAdicGenericElement[] {
    const ans = this._square_root({ ...options, all: true });
    return Array.isArray(ans) ? ans : [ans];
  }

  private _square_root(options?: {
    extend?: boolean;
    all?: boolean;
  }): pAdicGenericElement | pAdicGenericElement[] {
    const extend = options?.extend ?? true;
    const all = options?.all ?? false;

    // We first check trivial cases and precision
    if (this._exactZero) {
      return this;
    }
    const parent = this._parent;
    const p = parent.prime();
    if (this.is_zero() || (p === 2n && this.precision_relative() < 1 + 2 * parent.absolute_e())) {
      throw new PrecisionError(
        'not enough precision to be sure that this element has a square root'
      );
    }

    let ans: pAdicGenericElement | null = null;
    try {
      ans = this._nth_root(2n);
    } catch (err) {
      if (!(err instanceof ValueError) || err instanceof PrecisionError) {
        throw err;
      }
    }

    if (ans !== null) {
      let ans2 = ans.neg();
      // SageMath chooses the root deterministically: the one whose expansion is
      // smaller at the first index where the two differ.
      const E1 = ans.expansion();
      const E2 = ans2.expansion();
      let i = ans.parent().is_field() ? 0 : Number(ans.valuation());
      while (i < E1.length && i < E2.length) {
        const d1 = E1[i]!;
        const d2 = E2[i]!;
        if (d1 > d2) {
          const tmp = ans;
          ans = ans2;
          ans2 = tmp;
          break;
        }
        if (d1 < d2) {
          break;
        }
        i += 1;
      }
      return all ? [ans, ans2] : ans;
    }

    if (extend) {
      throw new NotImplementedError('extending using the sqrt function not yet implemented');
    }
    if (all) {
      return [];
    }
    throw new ValueError('element is not a square');
  }

  private _modPow(base: bigint, exp: bigint, mod: bigint): bigint {
    return modPow(base, exp, mod);
  }

  /**
   * Return an n-th root of this element.
   *
   * @param options.all - if true, return all n-th roots instead of just one.
   *
   * @see Reference: sage/rings/padics/padic_generic_element.pyx:nth_root
   */
  nth_root(n: bigint, options?: { all?: false }): pAdicGenericElement;
  nth_root(n: bigint, options: { all: true }): pAdicGenericElement[];
  nth_root(n: bigint, options?: { all?: boolean }): pAdicGenericElement | pAdicGenericElement[];
  nth_root(n: bigint, options?: { all?: boolean }): pAdicGenericElement | pAdicGenericElement[] {
    return this._nth_root_dispatch(n, options);
  }

  /**
   * Return the list of all n-th roots of this element.
   * @see Reference: sage/rings/padics/padic_generic_element.pyx:nth_root
   */
  nth_root_all(n: bigint): pAdicGenericElement[] {
    const ans = this._nth_root_dispatch(n, { all: true });
    return Array.isArray(ans) ? ans : [ans];
  }

  private _nth_root_dispatch(
    n: bigint,
    options?: { all?: boolean }
  ): pAdicGenericElement | pAdicGenericElement[] {
    const root = this._nth_root(n);
    if (!options?.all) {
      return root;
    }
    return this._parent.roots_of_unity(n).map((zeta) => root.mul(zeta));
  }

  /**
   * Core n-th root computation (returns a single root).
   *
   * Follows SageMath's decomposition `n = p^v * m`: the `m`-th root (with
   * `gcd(m, p) = 1`) is obtained by a Newton iteration on the inverse root
   * seeded by the residue field's `m`-th root, and the `p`-th root is then
   * extracted `v` times.  Each `p`-th root extraction loses exactly one digit
   * of relative precision, as in SageMath.
   *
   * @see Reference: sage/rings/padics/padic_generic_element.pyx:nth_root
   */
  private _nth_root(n: bigint): pAdicGenericElement {
    if (n === 0n) {
      throw new ValueError('n must be a nonzero integer');
    }
    if (n === 1n) {
      return this;
    }
    if (n < 0n) {
      return this.inv()._nth_root(-n);
    }

    // We first check trivial cases
    if (this._exactZero) {
      return this;
    }
    if (this.is_zero()) {
      throw new PrecisionError('not enough precision to be sure that this element is a nth power');
    }

    const p = this.prime();

    // n = p^v * m with gcd(m, p) = 1
    let v = 0n;
    let m = n;
    while (m % p === 0n) {
      m = m / p;
      v += 1n;
    }

    // We check the valuation
    const val = this._valuation;
    if (val % n !== 0n) {
      throw new ValueError('this element is not a nth power');
    }

    let N = this.precision_relative();
    let unit = this._unit;

    // The m-th root of the unit part (no precision loss).
    unit = this._unit_mth_root(unit, m, N);

    // Then extract p-th roots v times (one digit lost each time).
    for (let i = 0n; i < v; i += 1n) {
      if (N < 2 || (p === 2n && N < 3)) {
        throw new PrecisionError(
          'not enough precision to be sure that this element is a nth power'
        );
      }
      unit = this._unit_pth_root(unit, N);
      N -= 1;
    }

    return pAdicGenericElement.fromValuationUnit(
      this._parent,
      val / n,
      posMod(unit, bigPow(p, BigInt(N))),
      N
    );
  }

  /**
   * Return the m-th root of the unit `u` modulo `p^N`, where `gcd(m, p) = 1`.
   *
   * This is SageMath's Newton iteration on the *inverse* root:
   * `root <- root + (1/m) * root * (1 - u * root^m)`, which converges to
   * `u^(-1/m)`; the answer is its inverse.
   */
  private _unit_mth_root(u: bigint, m: bigint, N: number): bigint {
    if (m === 1n) {
      return u;
    }
    const p = this.prime();

    // The residue field root: SageMath calls abar.nth_root(m).
    const abar = posMod(u, p);
    let xbar: bigint;
    try {
      xbar = new Integer(abar).nth_root_mod(m, p).value;
    } catch (_err) {
      throw new ValueError('this element is not a nth power');
    }
    if (modPow(xbar, m, p) !== abar) {
      throw new ValueError('this element is not a nth power');
    }

    let root = modInverse(xbar, p);
    let curprec = 1;
    while (curprec < N) {
      curprec = Math.min(2 * curprec, N);
      const mod = bigPow(p, BigInt(curprec));
      const invm = modInverse(posMod(m, mod), mod);
      const defect = posMod(1n - posMod(u, mod) * modPow(root, m, mod), mod);
      root = posMod(root + invm * root * defect, mod);
    }
    return modInverse(root, bigPow(p, BigInt(N)));
  }

  /**
   * Return the p-th root of the unit `u` known modulo `p^N`, as a unit known
   * modulo `p^(N-1)`.
   *
   * Writes `u = omega * w` with `omega` the Teichmuller representative (whose
   * p-th root is `omega` itself, since `omega^p = omega`) and `w` a 1-unit;
   * `w` is a p-th power exactly when `w = 1 mod p^(1 + e/(p-1))`, and the root
   * is then lifted one digit at a time.
   */
  private _unit_pth_root(u: bigint, N: number): bigint {
    const p = this.prime();
    const modN = bigPow(p, BigInt(N));

    // Teichmuller part
    let omega = 1n;
    if (p !== 2n) {
      omega = posMod(u, p);
      let mod = p;
      for (let i = 1; i < N; i++) {
        mod = mod * p;
        omega = modPow(omega, p, mod);
      }
    }
    const w = posMod(u * modInverse(omega, modN), modN);

    // p-th powers among the 1-units are exactly 1 + p^(1 + e/(p-1)) Z_p, i.e.
    // 1 mod p^2 for p odd and 1 mod 8 for p = 2.
    const kStart = p === 2n ? 2 : 1;
    if (posMod(w - 1n, bigPow(p, BigInt(kStart + 1))) !== 0n) {
      throw new ValueError('this element is not a nth power');
    }

    // Invariant: y^p = w mod p^(k+1); each step gains one digit.
    let y = 1n;
    for (let k = kStart; k <= N - 2; k++) {
      const mod = bigPow(p, BigInt(k + 2));
      const defect = posMod(w - modPow(y, p, mod), mod);
      const c = posMod(defect / bigPow(p, BigInt(k + 1)), p);
      // y^(p-1) = 1 mod p because y = 1 mod p, so no inversion is needed.
      y = posMod(y + c * bigPow(p, BigInt(k)), mod);
    }

    return posMod(omega * y, bigPow(p, BigInt(N - 1)));
  }

  /**
   * Check if this is an n-th power.
   * @see Reference: sage/rings/padics/padic_generic_element.pyx:is_nth_power
   */
  is_nth_power(n: bigint): boolean {
    if (n === 0n) {
      throw new ValueError('n must be nonzero');
    }
    if (n === 1n || n === -1n) {
      return true;
    }
    if (n === 2n) {
      return this.is_square();
    }
    if (this._exactZero) {
      return true;
    }
    // The residue-field criterion alone is not sufficient when p divides n, so
    // we run the actual root extraction.
    try {
      this._nth_root(n);
    } catch (err) {
      if (err instanceof PrecisionError) {
        throw err;
      }
      if (err instanceof ValueError) {
        return false;
      }
      throw err;
    }
    return true;
  }

  /**
   * Return the logarithm.
   * The logarithm is defined for 1-units using the power series:
   *   log(1-x) = -x - x^2/2 - x^3/3 - ...
   *
   * For general units u = a * v where v is a 1-unit and a is a Teichmuller
   * representative, we define log(u) = log(v).
   *
   * @see Reference: sage/rings/padics/padic_generic_element.pyx:log
   */
  log(options?: {
    p_branch?: pAdicGenericElement;
    pi_branch?: pAdicGenericElement;
    aprec?: number;
  }): pAdicGenericElement {
    if (this.is_zero()) {
      throw new ValueError('logarithm is not defined at zero');
    }

    const p = this.prime();
    const R = this._parent;

    // For non-units, require a branch specification
    if (!this.is_padic_unit()) {
      if (!options?.p_branch && !options?.pi_branch) {
        throw new ValueError('you must specify a branch of the logarithm for non-units');
      }

      // Handle non-units: x = p^v * u where u is a unit
      // log(x) = v * log(p) + log(u)
      // The branch specifies what log(p) should be
      const v = this._valuation;
      const u = this.unit_part();

      // Get log of unit part (recursive call, but u is now a unit)
      const logU = u.log({ aprec: options?.aprec });

      // Determine log(p) from the branch
      let logP: pAdicGenericElement;
      if (options?.p_branch) {
        logP = options.p_branch;
      } else if (options?.pi_branch) {
        // pi_branch is for uniformizer, which for Zp/Qp is just p
        logP = options.pi_branch;
      } else {
        throw new ValueError('branch not specified');
      }

      // log(x) = v * log(p) + log(u)
      return logP.scalar_mul(v).add(logU);
    }

    // Get the 1-unit to take log of
    let y = this.unit_part();
    let x = R.one().sub(y);

    // If x has valuation 0, we need to raise to (q-1)th power first
    // to get a 1-unit
    let denom = 1n;
    if (x.is_zero()) {
      // y = 1, so log(y) = 0
      return R.zero();
    }

    if (x.valuation() <= 0n) {
      const q = p; // For unramified base field, q = p
      const qm1 = q - 1n;
      y = y.pow(qm1);
      x = R.one().sub(y);
      denom = qm1;
      if (x.is_zero()) {
        return R.zero();
      }
    }

    const aprec = options?.aprec ?? this.precision_relative();

    // Compute log using the series log(1-x) = -x - x^2/2 - x^3/3 - ...
    // This converges when v(x) > 0
    // We need to be careful about division by n when n is divisible by p

    let result = x.neg();
    let power = x;

    for (let n = 2; n <= aprec * 2; n++) {
      power = power.mul(x);

      // Divide by n, handling powers of p
      const nBig = BigInt(n);
      let nUnit = nBig;
      let nPVal = 0n;
      while (nUnit % p === 0n) {
        nUnit = nUnit / p;
        nPVal++;
      }

      // term = power / n = power * (1/nUnit) * p^(-nPVal)
      if (power.is_zero()) break;

      const modulus = bigPow(p, BigInt(aprec));
      const nUnitInv = modInverse(nUnit, modulus);
      const term = pAdicGenericElement.fromValuationUnit(
        R,
        power._valuation - nPVal,
        (power._unit * nUnitInv) % modulus,
        power._relprec
      );

      if (Number(term.valuation()) >= aprec) {
        break;
      }

      result = result.sub(term);
    }

    // Divide by denom if we had to raise to power first
    if (denom !== 1n) {
      // denom = p - 1, which is coprime to p
      const modulus = bigPow(p, BigInt(aprec));
      const denomInv = modInverse(denom, modulus);
      result = pAdicGenericElement.fromValuationUnit(
        R,
        result._valuation,
        (result._unit * denomInv) % modulus,
        result._relprec
      );
    }

    return result.add_bigoh(aprec);
  }

  /**
   * Return the exponential.
   * The exponential converges when v(x) > e/(p-1) where e is the ramification index.
   *
   * @see Reference: sage/rings/padics/padic_generic_element.pyx:exp
   */
  exp(options?: { aprec?: number }): pAdicGenericElement {
    const p = this.prime();
    const e = Number(this._parent.ramification_index());

    // Check convergence: need v(x) > e/(p-1)
    // For base rings, e = 1, so need v(x) * (p-1) > 1
    if (Number(this.valuation()) * (Number(p) - 1) <= e) {
      throw new ValueError('Exponential does not converge for that input.');
    }

    // The optimal absolute precision on exp(self) is the absolute precision
    // on self, capped by the precision cap of the parent.
    const maxprec = Math.min(this.precision_absolute(), this._parent.precision_cap());
    const aprec = options?.aprec === undefined || options.aprec > maxprec ? maxprec : options.aprec;

    if (this._exactZero || this.is_zero()) {
      return this._parent.one().add_bigoh(aprec);
    }

    const R = this._parent;

    // Use the generic algorithm from SageMath
    // We compute sum_{n=0}^N x^n/n! using Horner's method
    const xVal = Number(this.valuation());

    // N is the number of terms needed
    const N = Math.floor(aprec / (xVal - e / (Number(p) - 1)));

    // We compute x^N + N*x^(N-1) + N*(N-1)*x^(N-2) + ... + N!
    // Then divide by N! at the end
    // This avoids computing factorials that are divisible by p

    // First compute N!
    let nFactorial = 1n;
    let nFactorialVal = 0n;
    for (let i = 1; i <= N; i++) {
      const iBig = BigInt(i);
      // Factor out powers of p from i
      let iUnit = iBig;
      let iPVal = 0n;
      while (iUnit % p === 0n) {
        iUnit = iUnit / p;
        iPVal++;
      }
      nFactorial *= iUnit;
      nFactorialVal += iPVal;
      // Reduce mod p^(some large precision)
      const largeMod = bigPow(p, BigInt(aprec * 2));
      nFactorial = nFactorial % largeMod;
    }

    // Now use Horner's method for the series
    let result = this._parent.one();
    for (let n = N; n >= 1; n--) {
      // result = result * x / n + 1
      result = result.mul(this);

      // Divide by n (handling powers of p)
      const nBig = BigInt(n);
      let nUnit = nBig;
      let nPVal = 0n;
      while (nUnit % p === 0n) {
        nUnit = nUnit / p;
        nPVal++;
      }
      // result = result * (1/nUnit) * p^(-nPVal)
      const modulus = bigPow(p, BigInt(aprec));
      const nUnitInv = modInverse(nUnit, modulus);
      result = pAdicGenericElement.fromValuationUnit(
        R,
        result._valuation - nPVal,
        (result._unit * nUnitInv) % modulus,
        result._relprec
      );

      // Add 1
      result = result.add(this._parent.one());
    }

    return result.add_bigoh(aprec);
  }

  /**
   * Return the Teichmuller lift.
   * The Teichmuller lift of x is the unique (p-1)th root of unity t
   * such that t = x mod p.
   *
   * @see Reference: sage/rings/padics/padic_generic_element.pyx:teichmuller
   */
  teichmuller(): pAdicGenericElement {
    if (this.valuation() > 0n) {
      return pAdicGenericElement.exactZero(this._parent);
    }

    const p = this.prime();
    const prec = this._parent.precision_cap();

    // Start with the residue mod p
    let t = this._unit % p;

    // Newton iteration: t' = t - (t^p - t) / (p*t^(p-1) - 1)
    // Simplifies to: t' = t * (1 - (t^p - t) * (inverse of (p*t^(p-1) - 1)))
    // For Teichmuller, we use: t' = t^p (this is the correct iteration)
    // Actually, we solve t^(p-1) = 1, so iterate: t' = t - (t^(p-1) - 1) * t / ((p-1)*t^(p-1))
    // But the standard formula is: t' = t^p (valid because we're looking for fixed point of Frobenius)

    // Use the formula: t_{n+1} = t_n - t_n * (t_n^{p-1} - 1) / (p-1)
    // For p-adic lift, we use: t_{n+1} = t_n^p mod p^{n+1}
    let modulus = p;
    for (let i = 1; i < prec; i++) {
      const newModulus = modulus * p;
      t = this._modPow(t, p, newModulus);
      modulus = newModulus;
    }

    return pAdicGenericElement.fromValuationUnit(this._parent, 0n, t, prec);
  }

  /**
   * Return the Artin-Hasse exponential.
   * AH(x) = exp(x + x^p/p + x^{p^2}/p^2 + ...)
   * @see Reference: sage/rings/padics/padic_generic_element.pyx:artin_hasse_exp
   */
  artin_hasse_exp(prec?: number): pAdicGenericElement {
    // Reference: sage/rings/padics/padic_generic_element.pyx:artin_hasse_exp
    // The Artin-Hasse exponential converges for valuation >= 1

    if (this.valuation() < 1n) {
      throw new ValueError('Artin-Hasse exponential does not converge on this input');
    }

    const targetPrec = prec ?? Math.min(this.precision_absolute(), this._parent.precision_cap());
    if (targetPrec <= 1) {
      // For precision 1, AH(x) ≡ 1 (mod p)
      return pAdicGenericElement.fromValuationUnit(this._parent, 0n, 1n, 1);
    }

    const p = this.prime();

    // Compute the argument: x + x^p/p + x^{p^2}/p^2 + ...
    // This sum converges since each term has increasing valuation
    let arg: pAdicGenericElement = this;
    let xPower: pAdicGenericElement = this;

    // Add terms x^{p^i}/p^i until they become negligible
    for (let i = 1; i < targetPrec; i++) {
      xPower = xPower.pow(p);
      const term = xPower.scalar_div(bigPow(p, BigInt(i)));
      if (term.valuation() >= BigInt(targetPrec)) {
        break;
      }
      arg = arg.add(term);
    }

    // Now compute exp(arg) using the standard exponential
    return arg.exp();
  }

  /**
   * Return the norm to the base field.
   * @see Reference: sage/rings/padics/padic_generic_element.pyx:norm
   */
  norm(): pAdicGenericElement {
    // For Qp/Zp, the norm is just the element itself
    return this;
  }

  /**
   * Return the trace to the base field.
   * @see Reference: sage/rings/padics/padic_generic_element.pyx:trace
   */
  trace(): pAdicGenericElement {
    // For Qp/Zp, the trace is just the element itself
    return this;
  }

  /**
   * Return the minimal polynomial.
   * For elements of Qp or Zp, this is x - self.
   * @see Reference: sage/rings/padics/padic_generic_element.pyx:minimal_polynomial
   */
  minimal_polynomial(name: string = 'x'): {
    coefficients: pAdicGenericElement[];
    variable: string;
  } {
    // Reference: sage/rings/padics/padic_generic_element.pyx:minimal_polynomial
    // For elements of the base ring Qp or Zp, the minimal polynomial is just x - a
    // where a is the element
    return {
      coefficients: [this.neg(), this._parent.one()],
      variable: name,
    };
  }

  /**
   * Return the characteristic polynomial.
   * For elements of Qp or Zp, this is the same as the minimal polynomial.
   * @see Reference: sage/rings/padics/padic_generic_element.pyx:charpoly
   */
  charpoly(name: string = 'x'): { coefficients: pAdicGenericElement[]; variable: string } {
    // Reference: sage/rings/padics/padic_generic_element.pyx
    // For the base ring, charpoly = minpoly = x - self
    return this.minimal_polynomial(name);
  }

  /**
   * Return the multiplicative order.
   * @see Reference: sage/rings/padics/padic_generic_element.pyx:multiplicative_order
   */
  multiplicative_order(prec?: number): InfiniteOr<bigint> {
    // eslint-disable-next-line @typescript-eslint/no-this-alias
    let self: pAdicGenericElement = this;
    if (prec !== undefined) {
      self = self.add_bigoh(prec);
    }
    if (self.is_zero() || self.valuation() !== 0n) {
      return Number.POSITIVE_INFINITY;
    }

    const parent = this._parent;
    const p = parent.prime();

    // Compute the multiplicative order outside p
    let order = multiplicativeOrderModP(self.residue(1), p);
    const one = parent.one();
    self = self.div(self.teichmuller());
    if (self.eq(one)) {
      return order;
    }

    // Compute the multiplicative order at p.  Roots of unity of p-power order
    // exist only when (p-1) divides the absolute ramification index e.
    const e = BigInt(parent.absolute_e());
    if (e % (p - 1n) !== 0n) {
      return Number.POSITIVE_INFINITY;
    }
    let n = 0;
    let ee = e;
    while (ee % p === 0n) {
      ee = ee / p;
      n += 1;
    }
    for (let i = 0; i <= n; i++) {
      order *= p;
      self = self.pow(p);
      if (self.eq(one)) {
        return order;
      }
    }
    return Number.POSITIVE_INFINITY;
  }

  /**
   * Return the additive order truncated at the given precision.
   * If prec is given, returns 1 if the element is zero mod p^prec, otherwise infinity.
   * @see Reference: sage/rings/padics/padic_generic_element.pyx:additive_order
   */
  additive_order(prec?: number): InfiniteOr<bigint> {
    // Reference: sage/rings/padics/padic_generic_element.pyx:additive_order
    //   if self.is_zero(prec): return Integer(1) else: return infinity
    if (prec !== undefined) {
      // Check if zero at this precision
      if (this.is_zero() || this.valuation() >= BigInt(prec)) {
        return 1n;
      }
    } else if (this.is_zero()) {
      // Both exact and inexact zeros are indistinguishable from zero.
      return 1n;
    }

    return Number.POSITIVE_INFINITY;
  }

  // Arithmetic operations

  /**
   * Add two elements.
   * @see Reference: sage/rings/padics/padic_generic_element.pyx:__add__
   */
  add(other: pAdicGenericElement): pAdicGenericElement {
    if (this._exactZero) return other;
    if (other._exactZero) return this;

    const p = this.prime();

    // Determine common precision
    const minAbsPrec = Math.min(this.precision_absolute(), other.precision_absolute());

    if (this.is_zero() && other.is_zero()) {
      // Both are (inexact) zeros
      const newVal = BigInt(minAbsPrec);
      return pAdicGenericElement.fromValuationUnit(this._parent, newVal, 0n, 0);
    }

    // Align both operands at the common (possibly negative) valuation and add
    // the unit parts once, as `cadd` does in CR_template.pxi.
    const minVal = this._valuation < other._valuation ? this._valuation : other._valuation;
    const thisNormalized = this._unit * bigPow(p, this._valuation - minVal);
    const otherNormalized = other._unit * bigPow(p, other._valuation - minVal);

    const sum = thisNormalized + otherNormalized;

    if (sum === 0n) {
      return pAdicGenericElement.fromValuationUnit(this._parent, BigInt(minAbsPrec), 0n, 0);
    }

    const newVal = minVal + padic_valuation(sum, p);
    const unit = sum / bigPow(p, newVal - minVal);

    const newRelPrec = minAbsPrec - Number(newVal);
    const modulus = bigPow(p, BigInt(newRelPrec > 0 ? newRelPrec : 0));
    const reducedUnit = posMod(unit, modulus);

    return pAdicGenericElement.fromValuationUnit(
      this._parent,
      newVal,
      reducedUnit,
      newRelPrec > 0 ? newRelPrec : 0
    );
  }

  /**
   * Subtract two elements.
   * @see Reference: sage/rings/padics/padic_generic_element.pyx:__sub__
   */
  sub(other: pAdicGenericElement): pAdicGenericElement {
    return this.add(other.neg());
  }

  /**
   * Multiply two elements.
   * @see Reference: sage/rings/padics/padic_generic_element.pyx:__mul__
   */
  mul(other: pAdicGenericElement): pAdicGenericElement {
    if (this._exactZero || other._exactZero) {
      return pAdicGenericElement.exactZero(this._parent);
    }

    if (this.is_zero() || other.is_zero()) {
      const newPrec = Math.min(this.precision_absolute(), other.precision_absolute());
      return pAdicGenericElement.fromValuationUnit(this._parent, BigInt(newPrec), 0n, 0);
    }

    const p = this.prime();
    const newVal = this._valuation + other._valuation;
    const newRelPrec = Math.min(this._relprec, other._relprec);
    const modulus = bigPow(p, BigInt(newRelPrec));

    const product = (this._unit * other._unit) % modulus;

    return pAdicGenericElement.fromValuationUnit(this._parent, newVal, product, newRelPrec);
  }

  /**
   * Divide two elements.
   * @see Reference: sage/rings/padics/padic_generic_element.pyx:__truediv__
   */
  div(other: pAdicGenericElement): pAdicGenericElement {
    if (other._exactZero) {
      throw new ZeroDivisionError('cannot divide by zero');
    }
    if (other.is_zero()) {
      throw new PrecisionError('cannot divide by something indistinguishable from zero');
    }

    if (this._exactZero) {
      return pAdicGenericElement.exactZero(this._parent);
    }

    if (this.is_zero()) {
      const newPrec = Math.min(
        this.precision_absolute() - Number(other.valuation()),
        this._parent.precision_cap()
      );
      return pAdicGenericElement.fromValuationUnit(this._parent, BigInt(newPrec), 0n, 0);
    }

    const p = this.prime();
    const newVal = this._valuation - other._valuation;
    const newRelPrec = Math.min(this._relprec, other._relprec);
    const modulus = bigPow(p, BigInt(newRelPrec));

    const otherInv = modInverse(other._unit, modulus);
    const quotient = (this._unit * otherInv) % modulus;

    return pAdicGenericElement.fromValuationUnit(this._parent, newVal, quotient, newRelPrec);
  }

  /**
   * Return the negation.
   * @see Reference: sage/rings/padics/padic_generic_element.pyx:__neg__
   */
  neg(): pAdicGenericElement {
    if (this._exactZero) {
      return pAdicGenericElement.exactZero(this._parent);
    }
    if (this.is_zero()) {
      return pAdicGenericElement.fromValuationUnit(
        this._parent,
        this._valuation,
        0n,
        this._relprec
      );
    }

    const p = this.prime();
    const modulus = bigPow(p, BigInt(this._relprec));
    const negUnit = posMod(-this._unit, modulus);

    return pAdicGenericElement.fromValuationUnit(
      this._parent,
      this._valuation,
      negUnit,
      this._relprec
    );
  }

  /**
   * Multiply by a scalar (bigint).
   * @see Reference: sage/rings/padics/padic_generic_element.pyx
   */
  scalar_mul(c: bigint): pAdicGenericElement {
    if (c === 0n) {
      return pAdicGenericElement.exactZero(this._parent);
    }
    if (c === 1n) {
      return this;
    }
    const other = this._parent.__call__(c);
    return this.mul(other);
  }

  /**
   * Divide by a scalar (bigint).
   * @see Reference: sage/rings/padics/padic_generic_element.pyx
   */
  scalar_div(c: bigint): pAdicGenericElement {
    if (c === 0n) {
      throw new ZeroDivisionError('cannot divide by zero');
    }
    if (c === 1n) {
      return this;
    }
    const other = this._parent.__call__(c);
    return this.div(other);
  }

  /**
   * Return the multiplicative inverse.
   * @see Reference: sage/rings/padics/padic_generic_element.pyx:__invert__
   */
  inv(): pAdicGenericElement {
    return this._parent.one().div(this);
  }

  /**
   * Return self raised to power n.
   * @see Reference: sage/rings/padics/padic_generic_element.pyx:__pow__
   */
  pow(n: bigint): pAdicGenericElement {
    if (n === 0n) {
      return this._parent.one();
    }
    if (n < 0n) {
      return this.inv().pow(-n);
    }

    if (this._exactZero) {
      return pAdicGenericElement.exactZero(this._parent);
    }

    if (this.is_zero()) {
      const newPrec = Number(this._valuation) * Number(n) + (this._relprec > 0 ? this._relprec : 0);
      return pAdicGenericElement.fromValuationUnit(this._parent, BigInt(newPrec), 0n, 0);
    }

    const p = this.prime();
    const newVal = this._valuation * n;
    const modulus = bigPow(p, BigInt(this._relprec));
    const powUnit = this._modPow(this._unit, n, modulus);

    return pAdicGenericElement.fromValuationUnit(this._parent, newVal, powUnit, this._relprec);
  }

  /**
   * Return the absolute value (p^-valuation).
   * @see Reference: sage/rings/padics/padic_generic_element.pyx:abs
   */
  abs(): number {
    if (this.is_zero()) {
      return 0;
    }
    const p = Number(this.prime());
    const v = Number(this.valuation());
    return p ** -v;
  }

  /**
   * Compare two elements.
   * @see Reference: sage/rings/padics/padic_generic_element.pyx:__eq__
   */
  eq(other: pAdicGenericElement): boolean {
    // Handle exact zeros
    if (this._exactZero && other._exactZero) {
      return true;
    }
    if (this._exactZero) {
      return other.is_zero();
    }
    if (other._exactZero) {
      return this.is_zero();
    }

    // Compare up to the minimum precision
    const minPrec = Math.min(this.precision_absolute(), other.precision_absolute());

    // Both are zero up to this precision
    if (Number(this.valuation()) >= minPrec && Number(other.valuation()) >= minPrec) {
      return true;
    }

    // Different valuations (up to precision)
    if (this._valuation !== other._valuation) {
      return false;
    }

    // Same valuation, compare units up to common relative precision
    const relPrec = minPrec - Number(this._valuation);
    if (relPrec <= 0) {
      return true;
    }

    const p = this.prime();
    const modulus = bigPow(p, BigInt(relPrec));
    return this._unit % modulus === other._unit % modulus;
  }

  /**
   * Add O(p^n) - reduce precision.
   */
  add_bigoh(n: number): pAdicGenericElement {
    if (this._exactZero) {
      return pAdicGenericElement.fromValuationUnit(this._parent, BigInt(n), 0n, 0);
    }

    const currentAbsPrec = this.precision_absolute();
    if (n >= currentAbsPrec) {
      return this;
    }

    const newRelPrec = n - Number(this._valuation);
    if (newRelPrec <= 0) {
      return pAdicGenericElement.fromValuationUnit(this._parent, BigInt(n), 0n, 0);
    }

    const p = this.prime();
    const modulus = bigPow(p, BigInt(newRelPrec));
    const reducedUnit = this._unit % modulus;

    return pAdicGenericElement.fromValuationUnit(
      this._parent,
      this._valuation,
      reducedUnit,
      newRelPrec
    );
  }

  toString(): string {
    if (this._exactZero) {
      return '0';
    }

    const p = this.prime();
    const prec = this.precision_absolute();
    // SageMath's printer omits the exponent when the absolute precision is 1.
    // Reference: sage/rings/padics/padic_printing.pyx:1057-1066
    const bigOh = prec === 1 ? `O(${p})` : `O(${p}^${prec})`;

    if (this.is_zero()) {
      return bigOh;
    }

    // Build series representation.  For fields the expansion starts at
    // p^valuation, which may be negative (e.g. Qp(5)(1/25) is 5^-2 + O(5^18)).
    const expansion = this.expansion();
    const startExp = this._parent.is_field() ? Number(this._valuation) : 0;
    const terms: string[] = [];

    for (let i = 0; i < expansion.length; i++) {
      const coef = expansion[i]!;
      if (coef === 0n) {
        continue;
      }
      const e = startExp + i;
      if (e === 0) {
        terms.push(coef.toString());
      } else if (e === 1) {
        terms.push(coef === 1n ? `${p}` : `${coef}*${p}`);
      } else {
        terms.push(coef === 1n ? `${p}^${e}` : `${coef}*${p}^${e}`);
      }
    }

    if (terms.length === 0) {
      return bigOh;
    }

    return `${terms.join(' + ')} + ${bigOh}`;
  }
}
