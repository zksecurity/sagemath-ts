/**
 * @module sage/rings/padics/padic_generic
 * @description Generic p-adic rings and fields
 *
 * Port of: sage/rings/padics/padic_generic.py
 * Reference: reference/sage/src/sage/rings/padics/padic_generic.py
 */

import { Z_factor } from '@sagemath-ts/parigp-ts';
import { NotImplementedError, ValueError } from '../../errors.js';
import { pAdicGenericElement } from './padic_generic_element.js';

/**
 * Check if a number is prime (simple trial division for small numbers).
 */
function isPrime(n: bigint): boolean {
  if (n < 2n) return false;
  if (n === 2n) return true;
  if (n % 2n === 0n) return false;
  for (let i = 3n; i * i <= n; i += 2n) {
    if (n % i === 0n) return false;
  }
  return true;
}

function gcdBig(a: bigint, b: bigint): bigint {
  let x = a < 0n ? -a : a;
  let y = b < 0n ? -b : b;
  while (y !== 0n) {
    const t = y;
    y = x % y;
    x = t;
  }
  return x;
}

function powMod(base: bigint, exp: bigint, mod: bigint): bigint {
  let result = 1n;
  let b = ((base % mod) + mod) % mod;
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
 * Return the smallest primitive root modulo the prime p.
 * Factorisation of p-1 is delegated to PARI, as SageMath does.
 */
function primitiveRootModP(p: bigint): bigint {
  if (p === 2n) return 1n;
  const primes = Z_factor(p - 1n)
    .filter(([q]) => q > 0n)
    .map(([q]) => q);
  for (let g = 2n; g < p; g += 1n) {
    let ok = true;
    for (const q of primes) {
      if (powMod(g, (p - 1n) / q, p) === 1n) {
        ok = false;
        break;
      }
    }
    if (ok) return g;
  }
  throw new ValueError(`no primitive root modulo ${p}`);
}

/**
 * A generic p-adic ring or field.
 *
 * p-adic numbers are elements of the completion of the rationals
 * with respect to the p-adic valuation.
 *
 * @see Reference: sage/rings/padics/padic_generic.py:pAdicGeneric
 */
export class pAdicGeneric {
  protected readonly _p: bigint;
  protected readonly _prec: number;

  constructor(p: bigint, prec: number = 20) {
    if (!isPrime(p)) {
      throw new ValueError('p must be prime');
    }
    if (prec <= 0) {
      throw new ValueError('precision must be positive');
    }
    this._p = p;
    this._prec = prec;
  }

  /**
   * Return the prime p.
   * @see Reference: sage/rings/padics/padic_generic.py:prime
   */
  prime(): bigint {
    return this._p;
  }

  /**
   * Return the precision.
   * @see Reference: sage/rings/padics/padic_generic.py:precision_cap
   */
  precision_cap(): number {
    return this._prec;
  }

  /**
   * Return the ramification index.
   * @see Reference: sage/rings/padics/padic_generic.py:ramification_index
   */
  ramification_index(): bigint {
    return 1n; // For Qp and Zp
  }

  /**
   * Return the residue class degree.
   * @see Reference: sage/rings/padics/padic_generic.py:residue_class_degree
   */
  residue_class_degree(): bigint {
    return 1n; // For Qp and Zp
  }

  /**
   * Return the absolute degree.
   * @see Reference: sage/rings/padics/padic_generic.py:absolute_degree
   */
  absolute_degree(): bigint {
    return 1n; // For Qp and Zp
  }

  /**
   * Return the absolute e (ramification index times residue degree).
   * @see Reference: sage/rings/padics/padic_generic.py:absolute_e
   */
  absolute_e(): number {
    return Number(this.ramification_index());
  }

  /**
   * Return the absolute f (residue class degree).
   * @see Reference: sage/rings/padics/padic_generic.py:absolute_f
   */
  absolute_f(): number {
    return Number(this.residue_class_degree());
  }

  /**
   * Return the uniformizer (a generator of the maximal ideal).
   * For Zp and Qp, this is just p.
   * @see Reference: sage/rings/padics/padic_generic.py:uniformizer
   */
  uniformizer(): pAdicGenericElement {
    return new pAdicGenericElement(this, this._p);
  }

  /**
   * Return the uniformizer power series.
   * @see Reference: sage/rings/padics/padic_generic.py:uniformizer_pow
   */
  uniformizer_pow(n: bigint): pAdicGenericElement {
    if (n < 0n) {
      // For fields, we can have negative powers
      if (!this.is_field()) {
        // `padic_capped_relative_element.pyx:97`
        throw new ValueError('p divides the denominator');
      }
    }
    // Build p^n directly from its valuation and unit so that negative n works
    // (`this._p ** n` throws RangeError for negative exponents).
    return pAdicGenericElement.fromValuationUnit(this, n, 1n, this._prec);
  }

  /**
   * Return the n-th roots of unity contained in this ring/field.
   *
   * For `Q_p` (`p` odd) the roots of unity are exactly the Teichmuller
   * representatives, so there are `gcd(n, p-1)` of them; for `Q_2` they are
   * `{1, -1}`.
   *
   * @see Reference: sage/rings/padics/padic_generic.py:roots_of_unity
   */
  roots_of_unity(n?: bigint): pAdicGenericElement[] {
    const p = this._p;
    const groupOrder = p === 2n ? 2n : p - 1n;
    let d = n === undefined ? groupOrder : gcdBig(n, groupOrder);
    if (d <= 0n) {
      d = 1n;
    }

    if (p === 2n) {
      return d === 1n ? [this.one()] : [this.one(), this.one().neg()];
    }

    // The subgroup of order d of GF(p)^*, lifted by Teichmuller.
    const g = primitiveRootModP(p);
    const zeta = powMod(g, (p - 1n) / d, p);
    const result: pAdicGenericElement[] = [];
    let cur = 1n;
    for (let i = 0n; i < d; i += 1n) {
      result.push(this.teichmuller(cur));
      cur = (cur * zeta) % p;
    }
    return result;
  }

  /**
   * Return the residue field.
   * For Zp and Qp, this is Z/pZ (represented as bigints mod p).
   * @see Reference: sage/rings/padics/padic_generic.py:residue_field
   */
  residue_field(): { prime: bigint; size: bigint } {
    return { prime: this._p, size: this._p };
  }

  /**
   * Return the residue class field (alias for residue_field).
   * @see Reference: sage/rings/padics/padic_generic.py:residue_class_field
   */
  residue_class_field(): { prime: bigint; size: bigint } {
    return this.residue_field();
  }

  /**
   * Return the characteristic of the residue field.
   * @see Reference: sage/rings/padics/padic_generic.py:residue_characteristic
   */
  residue_characteristic(): bigint {
    return this._p;
  }

  /**
   * Check if this is a field.
   * Subclasses (pAdicRing, pAdicField) should override this.
   * @see Reference: sage/rings/padics/padic_generic.py:is_field
   */
  is_field(): boolean {
    // Default implementation - subclasses override
    // For the generic class, we assume it's a ring (not a field)
    return false;
  }

  /**
   * Return the characteristic.
   * @see Reference: sage/rings/padics/padic_generic.py:characteristic
   */
  characteristic(): bigint {
    return 0n;
  }

  /**
   * Return the fraction field.
   * Subclasses should override for proper behavior.
   * @see Reference: sage/rings/padics/padic_generic.py:fraction_field
   */
  fraction_field(): pAdicField {
    // Default implementation returns a pAdicField with same p and prec
    // Subclasses can override for more specific behavior
    return new pAdicField(this._p, this._prec);
  }

  /**
   * Return the integer ring.
   * Subclasses should override for proper behavior.
   * @see Reference: sage/rings/padics/padic_generic.py:integer_ring
   */
  integer_ring(): pAdicRing {
    // Default implementation returns a pAdicRing with same p and prec
    // Subclasses can override for more specific behavior
    return new pAdicRing(this._p, this._prec);
  }

  /**
   * Coerce an element to this ring.
   * @see Reference: sage/rings/padics/padic_generic.py:__call__
   */
  __call__(x: unknown, absprec?: number): pAdicGenericElement {
    if (x instanceof pAdicGenericElement) {
      if (x.parent() === this) {
        return x;
      }
      // Coerce from another p-adic ring
      return new pAdicGenericElement(this, x.lift(), absprec);
    }

    if (typeof x === 'bigint') {
      return new pAdicGenericElement(this, x, absprec);
    }

    if (typeof x === 'number') {
      if (!Number.isInteger(x)) {
        throw new ValueError('cannot coerce non-integer number to p-adic');
      }
      return new pAdicGenericElement(this, BigInt(x), absprec);
    }

    throw new ValueError(`cannot coerce ${typeof x} to p-adic`);
  }

  /**
   * Return zero.
   * @see Reference: sage/rings/padics/padic_generic.py:zero
   */
  zero(): pAdicGenericElement {
    return pAdicGenericElement.exactZero(this);
  }

  /**
   * Return one.
   * @see Reference: sage/rings/padics/padic_generic.py:one
   */
  one(): pAdicGenericElement {
    return new pAdicGenericElement(this, 1n);
  }

  /**
   * Return a random element.
   * @see Reference: sage/rings/padics/padic_generic.py:random_element
   */
  random_element(): pAdicGenericElement {
    // Generate a random integer mod p^prec
    const modulus = this._p ** BigInt(this._prec);
    // Simple random generation - not cryptographically secure
    const randomBits: bigint[] = [];
    for (let i = 0; i < this._prec; i++) {
      randomBits.push(BigInt(Math.floor(Math.random() * Number(this._p))));
    }
    let value = 0n;
    for (let i = 0; i < this._prec; i++) {
      value += randomBits[i] * this._p ** BigInt(i);
    }
    return new pAdicGenericElement(this, value);
  }

  /**
   * Return the Teichmuller lift of an element.
   * The Teichmuller lift of x is the unique (p-1)th root of unity t
   * such that t = x mod p.
   * @see Reference: sage/rings/padics/padic_generic.py:teichmuller
   */
  teichmuller(x: unknown): pAdicGenericElement {
    const elem = this.__call__(x);
    // Handle exact zero specially
    if (elem.is_zero()) {
      return this.zero();
    }
    if (elem.precision_absolute() <= 0) {
      throw new ValueError('not enough precision to determine Teichmuller representative');
    }
    if (elem.valuation() > 0n) {
      return this.zero();
    }
    return elem.teichmuller();
  }

  /**
   * Return all Teichmuller representatives.
   * These are the (p-1)th roots of unity plus zero.
   * @see Reference: sage/rings/padics/padic_generic.py:teichmuller_system
   */
  teichmuller_system(): pAdicGenericElement[] {
    const result: pAdicGenericElement[] = [];
    // Return Teichmuller lifts for 1, 2, ..., p-1
    for (let i = 1n; i < this._p; i++) {
      result.push(this.teichmuller(i));
    }
    return result;
  }

  /**
   * Return the extension by the polynomial f.
   * This creates a p-adic extension ring/field defined by the polynomial f.
   * @see Reference: sage/rings/padics/padic_generic.py:extension
   */
  extension(f: unknown, names: string): pAdicExtension {
    // Reference: sage/rings/padics/padic_generic.py:extension
    // Extract degree from polynomial f
    // For now, assume f is an array of coefficients or has a degree property
    let degree: number;
    if (Array.isArray(f)) {
      degree = f.length - 1;
    } else if (typeof f === 'object' && f !== null && 'degree' in f) {
      degree = (f as { degree: number }).degree;
    } else {
      throw new ValueError(
        'f must be a polynomial (array of coefficients) or have a degree property'
      );
    }
    return new pAdicExtension(this, degree);
  }

  /**
   * Return an unramified extension of degree n.
   * An unramified extension of degree n is defined by an irreducible
   * polynomial of degree n over Z/pZ lifted to Zp.
   * @see Reference: sage/rings/padics/padic_generic.py:unramified_extension
   */
  unramified_extension(n: number, names: string): pAdicExtension {
    // Reference: sage/rings/padics/padic_generic.py
    // For an unramified extension of degree n, we need an irreducible polynomial
    // of degree n over the residue field GF(p).
    // This is done using Conway polynomials or similar.
    return new pAdicExtension(this, n);
  }

  /**
   * Return a totally ramified extension of degree n.
   * A totally ramified extension of degree n is defined by an Eisenstein
   * polynomial of degree n.
   * @see Reference: sage/rings/padics/padic_generic.py:totally_ramified_extension
   */
  totally_ramified_extension(n: number, names: string): pAdicExtension {
    // Reference: sage/rings/padics/padic_generic.py
    // For a totally ramified extension, we use an Eisenstein polynomial
    // x^n - p (or a similar polynomial with p dividing the constant term)
    return new pAdicExtension(this, n);
  }

  toString(): string {
    return `p-adic Ring/Field with p=${this._p} and precision ${this._prec}`;
  }
}

/**
 * The ring of p-adic integers Zp.
 * @see Reference: sage/rings/padics/padic_base_leaves.py:pAdicRingCappedRelative
 */
export class pAdicRing extends pAdicGeneric {
  constructor(p: bigint, prec: number = 20) {
    super(p, prec);
  }

  /**
   * Check if this is a field (always false for Zp).
   * @see Reference: sage/rings/padics/padic_generic.py:is_field
   */
  override is_field(): boolean {
    return false;
  }

  /**
   * Return the fraction field Qp.
   * @see Reference: sage/rings/padics/padic_generic.py:fraction_field
   */
  override fraction_field(): pAdicField {
    return new pAdicField(this._p, this._prec);
  }

  /**
   * Return the integer ring (self).
   * @see Reference: sage/rings/padics/padic_generic.py:integer_ring
   */
  override integer_ring(): pAdicRing {
    return this;
  }

  toString(): string {
    return `${this._p}-adic Ring with capped relative precision ${this._prec}`;
  }
}

/**
 * The field of p-adic numbers Qp.
 * @see Reference: sage/rings/padics/padic_base_leaves.py:pAdicFieldCappedRelative
 */
export class pAdicField extends pAdicGeneric {
  constructor(p: bigint, prec: number = 20) {
    super(p, prec);
  }

  /**
   * Check if this is a field (always true for Qp).
   * @see Reference: sage/rings/padics/padic_generic.py:is_field
   */
  override is_field(): boolean {
    return true;
  }

  /**
   * Return the fraction field (self).
   * @see Reference: sage/rings/padics/padic_generic.py:fraction_field
   */
  override fraction_field(): pAdicField {
    return this;
  }

  /**
   * Return the integer ring Zp.
   * @see Reference: sage/rings/padics/padic_generic.py:integer_ring
   */
  override integer_ring(): pAdicRing {
    return new pAdicRing(this._p, this._prec);
  }

  toString(): string {
    return `${this._p}-adic Field with capped relative precision ${this._prec}`;
  }
}

/**
 * An extension of a p-adic ring or field.
 * @see Reference: sage/rings/padics/padic_extension_generic.py:pAdicExtensionGeneric
 */
export class pAdicExtension extends pAdicGeneric {
  private readonly _base: pAdicGeneric;
  private readonly _degree: number;

  constructor(base: pAdicGeneric, degree: number, prec?: number) {
    super(base.prime(), prec ?? base.precision_cap());
    this._base = base;
    this._degree = degree;
  }

  /**
   * Return the base ring.
   * @see Reference: sage/rings/padics/padic_extension_generic.py:base_ring
   */
  base_ring(): pAdicGeneric {
    return this._base;
  }

  /**
   * Return the degree of the extension.
   * @see Reference: sage/rings/padics/padic_extension_generic.py:degree
   */
  degree(): number {
    return this._degree;
  }

  /**
   * Return the absolute degree.
   * @see Reference: sage/rings/padics/padic_extension_generic.py:absolute_degree
   */
  override absolute_degree(): bigint {
    return BigInt(this._degree) * this._base.absolute_degree();
  }
}

/**
 * Create a p-adic ring Zp.
 * @see Reference: sage/rings/padics/factory.py:Zp
 */
export function Zp(p: bigint, prec: number = 20): pAdicRing {
  return new pAdicRing(p, prec);
}

/**
 * Create a p-adic field Qp.
 * @see Reference: sage/rings/padics/factory.py:Qp
 */
export function Qp(p: bigint, prec: number = 20): pAdicField {
  return new pAdicField(p, prec);
}
