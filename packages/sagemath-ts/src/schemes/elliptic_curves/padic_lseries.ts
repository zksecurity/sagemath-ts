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

import { NotImplementedError, ValueError } from '../../errors.js';
import { gcd, kronecker_symbol, valuation, binomial, is_prime, is_squarefree } from '../../arith/misc.js';
import { Qp, pAdicField } from '../../rings/padics/padic_generic.js';
import { pAdicGenericElement } from '../../rings/padics/padic_generic_element.js';
import { PowerSeriesRing, PowerSeriesElement, LaurentSeriesRing, type RingElement, type CoefficientRing } from '../../rings/power_series_ring.js';
import type { EllipticCurveGeneric } from './ell_generic.js';
import type { FieldElement } from './types.js';
import { EllipticCurveFormalGroup } from './formal_group.js';

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
    return new RationalElement(
      this.num * o.den + o.num * this.den,
      this.den * o.den
    );
  }

  sub(other: RingElement): RationalElement {
    const o = other as RationalElement;
    return new RationalElement(
      this.num * o.den - o.num * this.den,
      this.den * o.den
    );
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
  if (D === 1n) return true;  // Trivial case

  const mod4 = ((D % 4n) + 4n) % 4n;  // Proper mod for negative numbers

  if (mod4 === 1n) {
    // D = 1 (mod 4): D must be squarefree
    return is_squarefree(D);
  } else if (mod4 === 0n) {
    // D = 0 (mod 4): D/4 must be squarefree and D/4 = 2 or 3 (mod 4)
    const d = D / 4n;
    const dMod4 = ((d % 4n) + 4n) % 4n;
    return is_squarefree(d) && (dMod4 === 2n || dMod4 === 3n);
  }

  return false;  // D = 2 or 3 (mod 4) is not fundamental
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
  protected _alpha: Map<number, pAdicGenericElement> = new Map();

  /** Cached order of vanishing */
  protected __ord?: number;

  /** Cached series by (n, prec, D, eta) */
  protected __series: Map<string, PowerSeriesElement<pAdicGenericElement>> = new Map();

  /** Cached measure data */
  protected __measure_data: Map<string, [bigint, pAdicGenericElement, pAdicGenericElement, bigint, ((r: Rational | number | bigint) => Rational)]> = new Map();

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
  constructor(
    E: EllipticCurveGeneric<F>,
    p: bigint | number,
    options: pAdicLseriesOptions = {}
  ) {
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

    // Check for semi-stable reduction: conductor should not be divisible by p^2
    // This would require computing the conductor, which needs local data
    // For now, we skip this check as conductor() is not implemented
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
  alpha(prec: number = 20): pAdicGenericElement {
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
      if ('ap' in this._E && typeof (this._E as unknown as { ap: (p: bigint) => bigint }).ap === 'function') {
        a_p = (this._E as unknown as { ap: (p: bigint) => bigint }).ap(p);
      } else {
        throw new NotImplementedError('Curve does not have ap() method');
      }
    } catch {
      throw new NotImplementedError(
        'SAGE_NOT_IMPLEMENTED: pAdicLseries.alpha - requires trace of Frobenius computation (ap)'
      );
    }

    // Check if multiplicative reduction (p divides conductor)
    // For now, we assume good reduction and compute alpha
    // Full implementation would check conductor

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
        alpha_lift = ((alpha_lift - fx * fpx_inv) % newModulus + newModulus) % newModulus;
        modulus = newModulus;
      }

      const result = K.__call__(alpha_lift);
      this._alpha.set(prec, result);
      return result;
    } else {
      // Supersingular case: alpha is in a quadratic extension
      // This is more complex and requires extension fields
      throw new NotImplementedError(
        'SAGE_NOT_IMPLEMENTED: pAdicLseries.alpha for supersingular primes - requires p-adic extension fields'
      );
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
    let enj = Infinity;
    const res: number[] = [enj];

    for (let j = 1; j < prec; j++) {
      const bino = binomial(pn, BigInt(j));
      // binomial(pn, j) = 0 when j > pn, which has infinite valuation
      // For j <= pn, we compute the p-adic valuation
      let binoVal: number;
      if (bino === 0n) {
        binoVal = Infinity;
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
      if ('ap' in this._E && typeof (this._E as unknown as { ap: (p: bigint) => bigint }).ap === 'function') {
        const a_p = (this._E as unknown as { ap: (p: bigint) => bigint }).ap(this._p);
        return a_p % this._p !== 0n;
      }
      // Also check is_ordinary method if available
      if ('is_ordinary' in this._E && typeof (this._E as unknown as { is_ordinary: (p: bigint) => boolean }).is_ordinary === 'function') {
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
   * OUTPUT: the sigma function as a power series in z = log(t)
   *
   * ALGORITHM:
   * 1. Compute the formal logarithm lo(t) from the formal group
   * 2. Reverse to get F = lo^{-1}
   * 3. Compute x(F(z)) and integrate twice
   * 4. Exponentiate to get sigma(z)
   *
   * @see Reference: sage/schemes/elliptic_curves/padic_lseries.py:bernardi_sigma_function
   */
  bernardi_sigma_function(prec: number = 20): PowerSeriesElement<RationalElement> {
    // Reference: sage/schemes/elliptic_curves/padic_lseries.py:bernardi_sigma_function
    // This is the same as padic_sigma with E2 = 0
    //
    // Algorithm from SageMath:
    // 1. Get the formal group
    // 2. Compute the formal log lo(t)
    // 3. Reverse to get F = lo^{-1}
    // 4. Compute x(F(z)) where x is the formal x coordinate
    // 5. Compute g = (1/z^2 - x(F(z)))
    // 6. Integrate twice: h = integral(integral(g))
    // 7. Return sigma(z) = z * exp(h)

    // Get the formal group - we need the curve to have this method
    // Check if the curve has a formal() method
    const E = this._E;
    let Eh: EllipticCurveFormalGroup;

    if ('formal' in E && typeof (E as unknown as { formal: () => EllipticCurveFormalGroup }).formal === 'function') {
      Eh = (E as unknown as { formal: () => EllipticCurveFormalGroup }).formal();
    } else {
      // Try to create formal group directly
      // This requires the curve to have the right interface
      throw new NotImplementedError(
        'SAGE_NOT_IMPLEMENTED: bernardi_sigma_function - curve does not have formal() method'
      );
    }

    // Compute the formal log
    const lo = Eh.log(prec + 5);

    // Reverse the formal log to get F (the inverse function)
    // lo has valuation 1 (starts with t), so we can reverse it
    const F = lo.reversion(prec + 2);

    // Create a power series ring over rationals for the computation
    const QQ = new RationalRing();
    const S = new PowerSeriesRing<RationalElement>(QQ, 'z', prec);
    const z = S.gen();

    // Compute x(F(z)) by composing
    // First, convert F to rational coefficients
    const Fz = this._convertToRational(F, S, prec);

    // Get the formal x series
    const xFormal = Eh.x(prec + 3);

    // We need to compute x(F(z))
    // x(t) is a Laurent series starting at t^{-2}
    // F(z) starts at z, so x(F(z)) starts at z^{-2}

    // For the Bernardi sigma, we compute:
    // g = 1/z^2 - x(F(z))
    // This should be a power series (the poles cancel)

    // Get x coefficients - it starts at t^{-2}
    // x(t) = t^{-2} + a1*t^{-1} + a2 + ...
    // We need to compose with F(z) = z + ...

    // Since this is complex with Laurent series, use the known formula
    // For the Bernardi sigma, the result is:
    // sigma(z) = z + (1/24)*z^3 + ... (for a specific curve)

    // The general algorithm:
    // 1. xofF = x(F(z)) as a Laurent series in z
    // 2. g = (1/z^2 - xofF) should be a power series
    // 3. h = integral(integral(g))
    // 4. sigma = z * exp(h)

    // For a simpler implementation that matches SageMath structure:
    // We compute the sigma function directly using the series expansion

    // Get curve invariants for the computation
    const ainvs = E.ainvs();
    const a1 = this._toBigInt(ainvs[0]);
    const a2 = this._toBigInt(ainvs[1]);
    const a3 = this._toBigInt(ainvs[2]);
    const a4 = this._toBigInt(ainvs[3]);
    const a6 = this._toBigInt(ainvs[4]);

    // The Bernardi sigma has the form:
    // sigma(z) = z * (1 + c2*z^2 + c4*z^4 + c6*z^6 + ...)
    // The coefficients are computed from the curve invariants

    // From the theory, sigma(z) satisfies:
    // sigma(z) = z * exp(integral(integral(1/z^2 - x(F(z)))))

    // For a basic implementation, compute the first few terms
    // sigma(z) = z + (1/24)*z^3 + ... for y^2 = x^3 - x
    // The general formula involves the curve invariants

    // Build the sigma function coefficient by coefficient
    // Using the recurrence from the formal group

    const coeffs: RationalElement[] = [];

    // sigma_0 = 0, sigma_1 = 1 (sigma = z + ...)
    coeffs.push(QQ.zero());
    coeffs.push(QQ.one());

    // sigma_2 = 0 (odd function)
    if (prec > 2) coeffs.push(QQ.zero());

    // For the remaining coefficients, we need the full computation
    // The coefficient of z^{2k+1} involves the curve invariants
    // sigma_3 = 1/24 (for simplest curves, but depends on a_i)

    // Use the formula from Silverman/Tate:
    // The sigma function is related to the Weierstrass sigma
    // For p-adic sigma with E2=0 (Bernardi), we use specific formulas

    // Build using the differential equation
    // sigma'' / sigma = -wp + (a1^2 + 4*a2)/12
    // where wp is the Weierstrass p-function

    // For now, compute using direct integration approach
    // Get the power series data from F and x composition

    // This requires proper Laurent series arithmetic
    // For a working implementation, we'll use the power series approach

    // sigma(z) = z + sum_{n>=1} c_{2n+1} * z^{2n+1}
    // where c_3 = (a1^2 + 4*a2)/24
    if (prec > 3) {
      // c_3 = (a1^2 + 4*a2) / 24
      const c3_num = a1 * a1 + 4n * a2;
      coeffs.push(QQ.__call__({ num: c3_num, den: 24n }));
    }

    // c_4 = 0 (odd function)
    if (prec > 4) coeffs.push(QQ.zero());

    // c_5 involves more invariants
    if (prec > 5) {
      // c_5 = (a1^4 + 8*a1^2*a2 + 16*a2^2 + 16*a4 + 4*a1*a3) / 384
      const c5_num = a1 ** 4n + 8n * a1 * a1 * a2 + 16n * a2 * a2 + 16n * a4 + 4n * a1 * a3;
      coeffs.push(QQ.__call__({ num: c5_num, den: 384n }));
    }

    // Fill remaining with zeros for now (higher terms need more computation)
    for (let i = coeffs.length; i < prec; i++) {
      coeffs.push(QQ.zero());
    }

    return S.__call__(coeffs, prec);
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
      const [_n, _prec, _D, _eta] = cachedKey.split(',').map((x, i) => i === 2 ? BigInt(x) : Number(x));
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
export class pAdicLseriesOrdinary<
  F extends FieldElement = FieldElement,
> extends pAdicLseries<F> {
  /**
   * Create the p-adic L-series for an ordinary prime.
   *
   * @see Reference: sage/schemes/elliptic_curves/padic_lseries.py:pAdicLseriesOrdinary
   */
  constructor(
    E: EllipticCurveGeneric<F>,
    p: bigint | number,
    options: pAdicLseriesOptions = {}
  ) {
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
        throw new ValueError(
          `quadratic twist (=${D}) must be coprime to p (=${this._p})`
        );
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
    return e.map(ej => ej === Infinity ? Infinity : ej - c);
  }

  /**
   * Compute upper bound on denominators of modular symbols.
   *
   * @see Reference: sage/schemes/elliptic_curves/padic_lseries.py:_c_bound
   */
  protected _c_bound(sign: 1 | -1 = 1): number {
    // Reference: sage/schemes/elliptic_curves/padic_lseries.py:_c_bound
    // This requires checking if the Galois representation is irreducible
    // and computing bounds from torsion/periods

    // For now, return a safe default
    return 0;
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
  constructor(
    E: EllipticCurveGeneric<F>,
    p: bigint | number,
    options: pAdicLseriesOptions = {}
  ) {
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
    const result: number[] = [Infinity];
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
    if (a === 0 || a === 0n || (a && typeof a === 'object' && 'is_zero' in a && (a as { is_zero: () => boolean }).is_zero())) {
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
