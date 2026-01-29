/**
 * @module sage/rings/padics/padic_generic_element
 * @description Elements of p-adic rings and fields
 *
 * Port of: sage/rings/padics/padic_generic_element.pyx
 * Reference: reference/sage/src/sage/rings/padics/padic_generic_element.pyx
 */

import {
  NotImplementedError,
  PrecisionError,
  ValueError,
  ZeroDivisionError,
} from '../../errors.js';
import type { pAdicGeneric } from './padic_generic.js';

// Re-export PrecisionError for backwards compatibility
export { PrecisionError } from '../../errors.js';

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

  constructor(
    parent: pAdicGeneric,
    value: bigint = 0n,
    absprec?: number,
    relprec?: number
  ) {
    this._parent = parent;
    this._exactZero = false;
    this._relprec = relprec ?? parent.precision_cap();

    if (value === 0n) {
      // Zero element
      this._valuation = BigInt(absprec ?? parent.precision_cap());
      this._unit = 0n;
      this._exactZero = absprec === undefined && relprec === undefined;
    } else {
      const p = parent.prime();
      // Compute valuation
      this._valuation = padic_valuation(value, p);
      // Extract unit part
      let unit = value / bigPow(p, this._valuation);
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
      return Infinity;
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
  valuation(): bigint {
    if (this._exactZero) {
      throw new ValueError('valuation of exact zero is infinity');
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
  ordp(): bigint {
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
    return pAdicGenericElement.fromValuationUnit(
      this._parent,
      0n,
      this._unit,
      this._relprec
    );
  }

  /**
   * Return the normalized valuation.
   * This is the valuation divided by the ramification index.
   * @see Reference: sage/rings/padics/padic_generic_element.pyx:normalized_valuation
   */
  normalized_valuation(): bigint {
    const e = this._parent.ramification_index();
    return this.valuation() / e;
  }

  /**
   * Return the expansion as a list of coefficients.
   * The i-th element is the coefficient of p^i in the expansion.
   * @see Reference: sage/rings/padics/padic_generic_element.pyx:expansion
   */
  expansion(): bigint[] {
    if (this._exactZero) {
      return [];
    }
    const p = this.prime();
    const result: bigint[] = [];
    const val = Number(this._valuation);
    const prec = this.precision_absolute();

    // Fill with zeros up to valuation
    for (let i = 0; i < val && i < prec; i++) {
      result.push(0n);
    }

    // Extract digits from unit part
    let unit = this._unit;
    for (let i = val; i < prec; i++) {
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
   * Return the slice of coefficients.
   * @see Reference: sage/rings/padics/padic_generic_element.pyx:slice
   */
  slice(start: number, end: number): pAdicGenericElement {
    const p = this.prime();
    let value = 0n;

    for (let i = start; i < end && i < this.precision_absolute(); i++) {
      const coef = i >= 0 ? this.__getitem__(i) : 0n;
      value += coef * bigPow(p, BigInt(i - start));
    }

    const result = new pAdicGenericElement(this._parent, value);
    return result;
  }

  /**
   * Return the residue modulo p^n.
   * @see Reference: sage/rings/padics/padic_generic_element.pyx:residue
   */
  residue(n: number = 1): bigint {
    if (this.valuation() < 0n) {
      throw new ValueError(
        'element must have nonnegative valuation in order to compute residue'
      );
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
    if (this.is_zero()) {
      return true;
    }

    // A p-adic number is a square iff:
    // 1. Its valuation is even
    // 2. Its unit part is a quadratic residue mod p
    if (this._valuation % 2n !== 0n) {
      return false;
    }

    const p = this.prime();
    const unitMod = this._unit % p;

    if (p === 2n) {
      // For p = 2, need to check mod 8
      const unitMod8 = this._unit % 8n;
      return unitMod8 === 1n;
    }

    // Euler criterion: a is QR mod p iff a^((p-1)/2) = 1 mod p
    const exp = (p - 1n) / 2n;
    let result = 1n;
    let base = unitMod;
    let e = exp;
    while (e > 0n) {
      if (e % 2n === 1n) {
        result = (result * base) % p;
      }
      base = (base * base) % p;
      e = e / 2n;
    }
    return result === 1n;
  }

  /**
   * Return the square root if it exists.
   * Uses Hensel's lemma for lifting.
   * @see Reference: sage/rings/padics/padic_generic_element.pyx:sqrt
   */
  sqrt(): pAdicGenericElement {
    if (this.is_zero()) {
      return pAdicGenericElement.exactZero(this._parent);
    }

    if (!this.is_square()) {
      throw new ValueError('element is not a square');
    }

    const p = this.prime();
    const val = this._valuation;
    const halfVal = val / 2n;

    // Start with square root of unit part mod p
    const unitMod = this._unit % p;
    let sqrtMod = this._tonelliShanks(unitMod, p);

    // Hensel lift to full precision
    const prec = this._relprec;
    let modulus = p;
    let sqrt = sqrtMod;

    for (let i = 1; i < prec; i++) {
      const newModulus = modulus * p;
      // Hensel lifting: sqrt' = sqrt - (sqrt^2 - unit) / (2*sqrt)
      const sqrtSquared = (sqrt * sqrt) % newModulus;
      const unitModNew = this._unit % newModulus;
      const diff = posMod(sqrtSquared - unitModNew, newModulus);
      const twoSqrtInv = modInverse((2n * sqrt) % newModulus, newModulus);
      sqrt = posMod(sqrt - diff * twoSqrtInv, newModulus);
      modulus = newModulus;
    }

    return pAdicGenericElement.fromValuationUnit(
      this._parent,
      halfVal,
      sqrt,
      prec
    );
  }

  /**
   * Tonelli-Shanks algorithm for computing square roots mod p.
   */
  private _tonelliShanks(n: bigint, p: bigint): bigint {
    if (n === 0n) return 0n;
    if (p === 2n) return n % 2n;

    // Find Q and S such that p - 1 = Q * 2^S
    let q = p - 1n;
    let s = 0n;
    while (q % 2n === 0n) {
      q = q / 2n;
      s++;
    }

    // Find a quadratic non-residue z
    let z = 2n;
    while (this._legendreSymbol(z, p) !== -1n) {
      z++;
    }

    let m = s;
    let c = this._modPow(z, q, p);
    let t = this._modPow(n, q, p);
    let r = this._modPow(n, (q + 1n) / 2n, p);

    while (true) {
      if (t === 0n) return 0n;
      if (t === 1n) return r;

      // Find the least i such that t^(2^i) = 1
      let i = 1n;
      let temp = (t * t) % p;
      while (temp !== 1n) {
        temp = (temp * temp) % p;
        i++;
      }

      // Update variables
      const exp = 1n << (m - i - 1n);
      const b = this._modPow(c, exp, p);
      m = i;
      c = (b * b) % p;
      t = (t * c) % p;
      r = (r * b) % p;
    }
  }

  private _legendreSymbol(a: bigint, p: bigint): bigint {
    const result = this._modPow(a, (p - 1n) / 2n, p);
    if (result === p - 1n) return -1n;
    return result;
  }

  private _modPow(base: bigint, exp: bigint, mod: bigint): bigint {
    let result = 1n;
    base = base % mod;
    while (exp > 0n) {
      if (exp % 2n === 1n) {
        result = (result * base) % mod;
      }
      exp = exp / 2n;
      base = (base * base) % mod;
    }
    return result;
  }

  /**
   * Return the n-th root if it exists.
   * @see Reference: sage/rings/padics/padic_generic_element.pyx:nth_root
   */
  nth_root(n: bigint): pAdicGenericElement {
    // Reference: sage/rings/padics/padic_generic_element.pyx:nth_root
    if (n === 0n) {
      throw new ValueError('n must be a nonzero integer');
    }
    if (n === 1n) {
      return this;
    }
    if (n < 0n) {
      return this.inv().nth_root(-n);
    }
    if (n === 2n) {
      return this.sqrt();
    }

    // Check if this is zero
    if (this.is_zero()) {
      return this;
    }

    const p = this.prime();
    const val = this.valuation();

    // Check if valuation is divisible by n
    if (val % n !== 0n) {
      throw new ValueError('this element is not a nth power');
    }

    // Get the unit part
    const unit = this.unit_part();
    const unitLift = unit.residue();

    // Check if the residue is an n-th power mod p
    // Need to check if u^((p-1)/gcd(n, p-1)) ≡ 1 (mod p)
    const g = this._gcd(n, p - 1n);
    const exp = (p - 1n) / g;
    const residuePow = this._modPow(unitLift, exp, p);
    if (residuePow !== 1n) {
      throw new ValueError('this element is not a nth power');
    }

    // Find n-th root of residue mod p using Tonelli-Shanks style algorithm
    // For simplicity, we use brute force for small p
    let rootMod = 1n;
    for (let x = 1n; x < p; x++) {
      if (this._modPow(x, n, p) === unitLift) {
        rootMod = x;
        break;
      }
    }

    // Hensel lift to get n-th root
    const prec = this.precision_relative();
    let root = rootMod;
    let modulus = p;

    // Newton iteration: x_{n+1} = x_n - (x_n^n - a) / (n * x_n^{n-1})
    // = x_n * (1 - (1 - a/x_n^n) / n)
    // = (x_n * (n-1) + a/x_n^{n-1}) / n
    for (let i = 1; i < prec; i++) {
      const newModulus = modulus * p;
      const xPowN = this._modPow(root, n, newModulus);
      const diff = posMod(unitLift * bigPow(p, BigInt(i)) - xPowN, newModulus);

      // We need to solve: x' = x + delta where (x+delta)^n ≡ a (mod p^{i+1})
      // Using x^n + n*x^{n-1}*delta ≡ a (mod p^{i+1})
      // delta ≡ (a - x^n) / (n * x^{n-1}) (mod p)
      const xPowNm1 = this._modPow(root, n - 1n, p);
      const nInv = modInverse(n % p, p);
      const delta = (diff / modulus) * nInv * modInverse(xPowNm1, p);
      root = posMod(root + delta * modulus, newModulus);
      modulus = newModulus;
    }

    // Construct the result with the correct valuation
    return pAdicGenericElement.fromValuationUnit(this._parent, val / n, root, prec);
  }

  private _gcd(a: bigint, b: bigint): bigint {
    a = a < 0n ? -a : a;
    b = b < 0n ? -b : b;
    while (b !== 0n) {
      const t = b;
      b = a % b;
      a = t;
    }
    return a;
  }

  /**
   * Check if this is an n-th power.
   * @see Reference: sage/rings/padics/padic_generic_element.pyx:is_nth_power
   */
  is_nth_power(n: bigint): boolean {
    // Reference: sage/rings/padics/padic_generic_element.pyx
    if (n === 0n) {
      throw new ValueError('n must be nonzero');
    }
    if (n === 1n || n === -1n) {
      return true;
    }
    if (n < 0n) {
      n = -n;
    }
    if (n === 2n) {
      return this.is_square();
    }

    if (this.is_zero()) {
      return true;
    }

    const p = this.prime();
    const val = this.valuation();

    // Valuation must be divisible by n
    if (val % n !== 0n) {
      return false;
    }

    // Unit part must be an n-th power
    const unit = this.unit_part();
    const unitLift = unit.residue();

    // Check if u^((p-1)/gcd(n, p-1)) ≡ 1 (mod p)
    const g = this._gcd(n, p - 1n);
    const exp = (p - 1n) / g;
    const residuePow = this._modPow(unitLift, exp, p);

    return residuePow === 1n;
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
        throw new ValueError(
          'you must specify a branch of the logarithm for non-units'
        );
      }

      // Handle non-units: x = p^v * u where u is a unit
      // log(x) = v * log(p) + log(u)
      // The branch specifies what log(p) should be
      const v = this.valuation();
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

    if (this._exactZero || this.is_zero()) {
      return this._parent.one();
    }

    const aprec =
      options?.aprec ??
      Math.min(this.precision_absolute(), this._parent.precision_cap());
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
    let arg = this;
    let xPower = this;

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
  minimal_polynomial(name: string = 'x'): { coefficients: pAdicGenericElement[]; variable: string } {
    // Reference: sage/rings/padics/padic_generic_element.pyx:minimal_polynomial
    // For elements of the base ring Qp or Zp, the minimal polynomial is just x - a
    // where a is the element
    return {
      coefficients: [this.neg(), this._parent.one()],
      variable: name
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
  multiplicative_order(): bigint {
    if (this.is_zero()) {
      throw new ValueError('multiplicative order of zero is not defined');
    }
    if (!this.is_unit()) {
      throw new ValueError(
        'multiplicative order is only defined for units in a ring'
      );
    }

    // For Qp, if it's not a root of unity, order is infinite
    // A unit is a root of unity iff it's a Teichmuller lift
    const p = this.prime();
    const order = p - 1n;

    // Check if this^(p-1) = 1
    const power = this.pow(order);
    if (power.is_one()) {
      // Find the actual order by testing divisors of p-1
      for (let d = 1n; d < order; d++) {
        if (order % d === 0n && this.pow(d).is_one()) {
          return d;
        }
      }
      return order;
    }

    throw new ValueError('element has infinite multiplicative order');
  }

  /**
   * Return the additive order truncated at the given precision.
   * If prec is given, returns 1 if the element is zero mod p^prec, otherwise infinity.
   * @see Reference: sage/rings/padics/padic_generic_element.pyx:additive_order
   */
  additive_order(prec?: number): bigint | 'Infinity' {
    // Reference: sage/rings/padics/padic_generic_element.pyx:additive_order
    // p-adic numbers have infinite additive order unless zero at the given precision

    if (prec !== undefined) {
      // Check if zero at this precision
      if (this.is_zero() || this.valuation() >= BigInt(prec)) {
        return 1n;
      }
    } else if (this._exactZero) {
      return 1n;
    }

    return 'Infinity';
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
    const minAbsPrec = Math.min(
      this.precision_absolute(),
      other.precision_absolute()
    );

    if (this.is_zero() && other.is_zero()) {
      // Both are (inexact) zeros
      const newVal = BigInt(minAbsPrec);
      return pAdicGenericElement.fromValuationUnit(this._parent, newVal, 0n, 0);
    }

    // Lift both to integers, add, then create new element
    const thisLift =
      this._unit * bigPow(p, this._valuation >= 0n ? this._valuation : 0n);
    const otherLift =
      other._unit * bigPow(p, other._valuation >= 0n ? other._valuation : 0n);

    // Handle negative valuations (field elements)
    let minVal = this._valuation < other._valuation ? this._valuation : other._valuation;
    const thisNormalized = thisLift * bigPow(p, this._valuation - minVal);
    const otherNormalized = otherLift * bigPow(p, other._valuation - minVal);

    const sum = thisNormalized + otherNormalized;

    if (sum === 0n) {
      return pAdicGenericElement.fromValuationUnit(
        this._parent,
        BigInt(minAbsPrec),
        0n,
        0
      );
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
      const newPrec = Math.min(
        this.precision_absolute(),
        other.precision_absolute()
      );
      return pAdicGenericElement.fromValuationUnit(
        this._parent,
        BigInt(newPrec),
        0n,
        0
      );
    }

    const p = this.prime();
    const newVal = this._valuation + other._valuation;
    const newRelPrec = Math.min(this._relprec, other._relprec);
    const modulus = bigPow(p, BigInt(newRelPrec));

    const product = (this._unit * other._unit) % modulus;

    return pAdicGenericElement.fromValuationUnit(
      this._parent,
      newVal,
      product,
      newRelPrec
    );
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
      throw new PrecisionError(
        'cannot divide by something indistinguishable from zero'
      );
    }

    if (this._exactZero) {
      return pAdicGenericElement.exactZero(this._parent);
    }

    if (this.is_zero()) {
      const newPrec = Math.min(
        this.precision_absolute() - Number(other.valuation()),
        this._parent.precision_cap()
      );
      return pAdicGenericElement.fromValuationUnit(
        this._parent,
        BigInt(newPrec),
        0n,
        0
      );
    }

    const p = this.prime();
    const newVal = this._valuation - other._valuation;
    const newRelPrec = Math.min(this._relprec, other._relprec);
    const modulus = bigPow(p, BigInt(newRelPrec));

    const otherInv = modInverse(other._unit, modulus);
    const quotient = (this._unit * otherInv) % modulus;

    return pAdicGenericElement.fromValuationUnit(
      this._parent,
      newVal,
      quotient,
      newRelPrec
    );
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
      const newPrec =
        Number(this._valuation) * Number(n) +
        (this._relprec > 0 ? this._relprec : 0);
      return pAdicGenericElement.fromValuationUnit(
        this._parent,
        BigInt(newPrec),
        0n,
        0
      );
    }

    const p = this.prime();
    const newVal = this._valuation * n;
    const modulus = bigPow(p, BigInt(this._relprec));
    const powUnit = this._modPow(this._unit, n, modulus);

    return pAdicGenericElement.fromValuationUnit(
      this._parent,
      newVal,
      powUnit,
      this._relprec
    );
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
    return Math.pow(p, -v);
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
    const minPrec = Math.min(
      this.precision_absolute(),
      other.precision_absolute()
    );

    // Both are zero up to this precision
    if (
      Number(this.valuation()) >= minPrec &&
      Number(other.valuation()) >= minPrec
    ) {
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
      return pAdicGenericElement.fromValuationUnit(
        this._parent,
        BigInt(n),
        0n,
        0
      );
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

    if (this.is_zero()) {
      return `O(${p}^${prec})`;
    }

    // Build series representation
    const expansion = this.expansion();
    const terms: string[] = [];

    for (let i = 0; i < expansion.length; i++) {
      const coef = expansion[i];
      if (coef !== 0n) {
        if (i === 0) {
          terms.push(coef.toString());
        } else if (i === 1) {
          if (coef === 1n) {
            terms.push(`${p}`);
          } else {
            terms.push(`${coef}*${p}`);
          }
        } else {
          if (coef === 1n) {
            terms.push(`${p}^${i}`);
          } else {
            terms.push(`${coef}*${p}^${i}`);
          }
        }
      }
    }

    if (terms.length === 0) {
      terms.push('0');
    }

    return terms.join(' + ') + ` + O(${p}^${prec})`;
  }
}
