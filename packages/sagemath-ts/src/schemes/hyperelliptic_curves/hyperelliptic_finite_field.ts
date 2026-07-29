/**
 * @module sage/schemes/hyperelliptic_curves/hyperelliptic_finite_field
 * @description Hyperelliptic curves over a finite field
 *
 * Port of: `sage/schemes/hyperelliptic_curves/hyperelliptic_finite_field.py`
 *
 * AUTHORS (upstream): David Kohel (2006), Robert Bradshaw (2007),
 * Alyson Deines / Marina Gresham / Gagan Sekhon (2010), Daniel Krenn (2011),
 * Jean-Pierre Flori / Jan Tuitman (2013), Kiran Kedlaya (2016),
 * Dean Bisogno (2017).
 *
 * ## What is and is not ported
 *
 * Upstream computes the characteristic polynomial of Frobenius with three
 * algorithms:
 *
 * - `'matrix'` — the `hypellfrob` C++ library (Kedlaya's algorithm).  Not ported.
 * - `'pari'` — PARI's `hyperellcharpoly`.  Not available in `parigp-ts`.
 * - `'cardinalities'` — count points over the first `g` extensions and apply
 *   Newton's identities.  **Ported in full**, and used as the default here.
 *
 * The Frobenius polynomial is uniquely determined by the curve, so the values
 * returned agree with SageMath; only the running time differs.
 */

import { binomial, isqrt } from '../../arith/misc.js';
import { NotImplementedError, ValueError } from '../../errors.js';
import { GFpn } from '../../rings/finite_rings/finite_field_extension.js';
import type { Polynomial, RingElement } from '../../rings/polynomial/polynomial_element.js';
import { PolynomialRing } from '../../rings/polynomial/polynomial_ring.js';
import {
  type HyperellipticBaseRing,
  absolute_trace_is_zero,
  cardinality_of,
  characteristic_of,
  degree_of,
  div_elements,
  element_to_bigint,
  field_embedding,
  is_square_of,
  iterate_field,
  pow_element,
  sqrt_all_of,
} from './field_ops.js';
import {
  HyperellipticCurve_generic,
  type HyperellipticPoint,
  sage_poly_repr,
} from './hyperelliptic_generic.js';

/**
 * A polynomial over `ZZ` represented by its coefficient list in ascending
 * degree order.  Used for the characteristic polynomial of Frobenius, which
 * Sage returns as an element of `ZZ['x']`.
 *
 * @see Deviation: the port has no `ZZ[x]` polynomial type, so integer
 * polynomials are returned as `bigint[]`; {@link zz_poly_repr} renders them
 * exactly the way Sage prints them.
 */
export type ZZPoly = bigint[];

/** Render an integer polynomial the way Sage prints elements of `ZZ['x']`. */
export function zz_poly_repr(coeffs: ZZPoly, name = 'x'): string {
  let d = coeffs.length - 1;
  while (d >= 0 && coeffs[d] === 0n) {
    d--;
  }
  const strs: (string | null)[] = [];
  for (let i = 0; i <= d; i++) {
    const c = coeffs[i]!;
    strs.push(c === 0n ? null : c.toString());
  }
  return sage_poly_repr(strs, name);
}

/** A rational function `num / den` over `ZZ`, as returned by `zeta_function`. */
export class ZZRationalFunction {
  constructor(
    readonly numerator: ZZPoly,
    readonly denominator: ZZPoly,
    readonly variable = 'x'
  ) {}

  toString(): string {
    return `(${zz_poly_repr(this.numerator, this.variable)})/(${zz_poly_repr(
      this.denominator,
      this.variable
    )})`;
  }
}

/**
 * Hyperelliptic curve over a finite field.
 */
export class HyperellipticCurve_finite_field<
  C extends RingElement,
> extends HyperellipticCurve_generic<C> {
  private __points: HyperellipticPoint<C>[] | null = null;
  private __frobenius_polynomial: ZZPoly | null = null;
  private __cardinality = new Map<number, bigint>();
  private __cartier: CartierData<C> | null = null;

  /* ------------------------------------------------------------------ */
  /* p-adic precision bounds                                            */
  /* ------------------------------------------------------------------ */

  /**
   * Bound on the number of `p`-adic digits needed to recover the Frobenius
   * polynomial from the characteristic polynomial of the Frobenius matrix.
   *
   * Port of `hyperelliptic_finite_field.py:73-125`.  Sage evaluates
   * `M = 2*binomial(2g, g)*sqrt(q)^g` in 53-bit floating point; this port
   * evaluates `M^2 = 4*binomial(2g, g)^2*q^g` exactly in `ZZ`, which agrees
   * with Sage whenever Sage's floating point is accurate.
   */
  _frobenius_coefficient_bound_charpoly(): number {
    const K = this.base_ring();
    const p = characteristic_of(K);
    const q = cardinality_of(K);
    if (q === null) {
      throw new ValueError('base ring must be finite');
    }
    const g = BigInt(this.genus());
    const c = binomial(2n * g, g);
    const M2 = 4n * c * c * q ** g;
    return bound_from_square(M2, p);
  }

  /**
   * Bound on the number of `p`-adic digits needed to recover `N_1, ..., N_n`
   * from the traces of the powers of the Frobenius matrix.
   *
   * Port of `hyperelliptic_finite_field.py:127-180` (`M = 4 g sqrt(q)^n`).
   */
  _frobenius_coefficient_bound_traces(n = 1): number {
    const K = this.base_ring();
    const p = characteristic_of(K);
    const q = cardinality_of(K);
    if (q === null) {
      throw new ValueError('base ring must be finite');
    }
    const g = BigInt(this.genus());
    const M2 = 16n * g * g * q ** BigInt(n);
    return bound_from_square(M2, p);
  }

  /* ------------------------------------------------------------------ */
  /* Frobenius                                                          */
  /* ------------------------------------------------------------------ */

  /** `hyperelliptic_finite_field.py:182-255` — requires the `hypellfrob` library. */
  frobenius_matrix_hypellfrob(_N?: number): never {
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: frobenius_matrix_hypellfrob (the hypellfrob library is not ported)'
    );
  }

  /** `hyperelliptic_finite_field.py:257-308` — requires the `hypellfrob` library. */
  frobenius_matrix(_options?: { N?: number; algorithm?: string }): never {
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: frobenius_matrix (the hypellfrob library is not ported)'
    );
  }

  /** `hyperelliptic_finite_field.py:369-431` — requires the `hypellfrob` library. */
  _frobenius_polynomial_matrix(_options?: { algorithm?: string }): never {
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: frobenius_polynomial(algorithm="matrix") (the hypellfrob library is not ported)'
    );
  }

  /** `hyperelliptic_finite_field.py:433-491` — requires PARI's `hyperellcharpoly`. */
  _frobenius_polynomial_pari(): never {
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: frobenius_polynomial(algorithm="pari") (parigp-ts does not provide hyperellcharpoly)'
    );
  }

  /**
   * Characteristic polynomial of Frobenius from the point counts over the
   * first `g` extensions of the base field.
   *
   * Port of `hyperelliptic_finite_field.py:310-367`.
   */
  _frobenius_polynomial_cardinalities(a?: bigint[]): ZZPoly {
    const g = this.genus();
    const q = cardinality_of(this.base_ring());
    if (q === null) {
      throw new ValueError('base ring must be finite');
    }

    const counts = a ?? this.count_points(g);

    // s_i = a_i - q^(i+1) - 1
    const s = counts.map((ai, i) => ai - q ** BigInt(i + 1) - 1n);

    const coeffs: bigint[] = [1n];
    for (let i = 1; i <= g; i++) {
      let c = 0n;
      for (let j = 0; j < i; j++) {
        c += s[i - 1 - j]! * coeffs[j]!;
      }
      if (c % BigInt(i) !== 0n) {
        throw new ValueError(
          `frobenius polynomial coefficient ${c}/${i} is not an integer; the point counts are inconsistent`
        );
      }
      coeffs.push(c / BigInt(i));
    }
    for (let i = 1; i <= g; i++) {
      coeffs.push(coeffs[g - i]! * q ** BigInt(i));
    }

    // ZZ['x'](coeffs).reverse()
    return coeffs.slice().reverse();
  }

  /**
   * Compute the characteristic polynomial of Frobenius, as an element of
   * `ZZ[x]` (returned as an ascending coefficient list).
   *
   * Port of `hyperelliptic_finite_field.py:506-630`.
   *
   * @see Deviation: with `hypellfrob` and PARI's `hyperellcharpoly` absent,
   * the default algorithm is always `'cardinalities'`.  The result is the same
   * polynomial; only the running time differs.
   */
  frobenius_polynomial(options?: { algorithm?: 'cardinalities' | 'matrix' | 'pari' }): ZZPoly {
    const algorithm = options?.algorithm;
    if (algorithm === 'matrix') {
      return this._frobenius_polynomial_matrix();
    }
    if (algorithm === 'pari') {
      return this._frobenius_polynomial_pari();
    }
    if (algorithm !== undefined && algorithm !== 'cardinalities') {
      throw new ValueError(`unknown algorithm ${algorithm}`);
    }
    if (algorithm === 'cardinalities') {
      return this._frobenius_polynomial_cardinalities();
    }
    if (this.__frobenius_polynomial === null) {
      this.__frobenius_polynomial = this._frobenius_polynomial_cardinalities();
    }
    return this.__frobenius_polynomial;
  }

  /**
   * `hyperelliptic_finite_field.py:1435-1471`
   *
   * `Z(x) = P^rev(x) / ((1 - x)(1 - q x))`.
   */
  zeta_function(): ZZRationalFunction {
    const q = cardinality_of(this.base_ring());
    if (q === null) {
      throw new ValueError('base ring must be finite');
    }
    const P = this.frobenius_polynomial();
    // P.reverse(): coefficient of x^k becomes P[deg - k]
    const num = P.slice().reverse();
    // (1-x)*(1-q*x) = 1 - (1+q) x + q x^2
    const den = [1n, -(1n + q), q];
    return new ZZRationalFunction(num, den);
  }

  /* ------------------------------------------------------------------ */
  /* Enumerating points                                                 */
  /* ------------------------------------------------------------------ */

  /**
   * List points by enumerating over `x` and solving the resulting quadratic
   * for `y`.
   *
   * Port of `hyperelliptic_finite_field.py:632-752`.
   */
  _points_fast_sqrt(): HyperellipticPoint<C>[] {
    const K = this.base_ring();
    const [f, h] = this.hyperelliptic_polynomials();
    const one = K.one() as C;
    const zero = K.zero() as C;
    const two = K.__call__(2n) as C;
    const char = characteristic_of(K);

    const points: HyperellipticPoint<C>[] = [];

    // start with the points at infinity
    if (this.defining_polynomial_eval(zero, one, zero).isZero()) {
      points.push(this.point([zero, one, zero], { check: true }));
    }
    const degP = this.defining_polynomial_degree();
    if (degP > 2) {
      // P(1, y, 0) = r*y + s
      const s = this.defining_polynomial_eval(one, zero, zero);
      const r = this.defining_polynomial_eval(one, one, zero).sub(s) as C;
      if (!r.isZero()) {
        points.push(this.point([one, div_elements(s.neg() as C, r), zero], { check: true }));
      }
    } else if (char === 2n) {
      for (const y of iterate_field(K)) {
        if (this.defining_polynomial_eval(one, y, zero).isZero()) {
          points.push(this.point([one, y, zero], { check: true }));
        }
      }
    } else {
      // P(1, y, 0) = y^2 + r*y + s
      const s = f.getCoeff(2).neg() as C;
      const r = h.getCoeff(1);
      const d = div_elements(r.mul(r) as C, two.mul(two) as C).sub(s) as C;
      const halfR = div_elements(r.neg() as C, two);
      if (d.isZero()) {
        points.push(this.point([one, halfR, zero], { check: true }));
      } else if (is_square_of(K, d)) {
        const sqrtd = sqrt_all_of(K, d)[0]!;
        points.push(this.point([one, halfR.add(sqrtd) as C, zero], { check: true }));
        points.push(this.point([one, halfR.sub(sqrtd) as C, zero], { check: true }));
      }
    }

    if (char === 2n) {
      if (h.isZero()) {
        for (const x of iterate_field(K)) {
          points.push(this.point([x, sqrt_all_of(K, f.evaluate(x))[0]!, one], { check: true }));
        }
      } else {
        // Artin-Schreier 2-roots: a_sqrts[x^2 + x] = x
        const roots = new Map<string, C>();
        for (const x of iterate_field(K)) {
          roots.set(key_of(x.mul(x).add(x) as C), x);
        }
        for (const x of iterate_field(K)) {
          const b = h.evaluate(x);
          const c = f.evaluate(x);
          if (!b.isZero()) {
            const r = roots.get(key_of(div_elements(c, b.mul(b) as C)));
            if (r !== undefined) {
              points.push(this.point([x, r.mul(b) as C, one], { check: true }));
              points.push(this.point([x, r.mul(b).add(b) as C, one], { check: true }));
            }
          } else {
            points.push(this.point([x, sqrt_all_of(K, c)[0]!, one], { check: true }));
          }
        }
      }
    } else if (h.isZero()) {
      for (const x of iterate_field(K)) {
        const y2 = f.evaluate(x);
        if (y2.isZero()) {
          points.push(this.point([x, y2, one], { check: true }));
        } else if (is_square_of(K, y2)) {
          const y = sqrt_all_of(K, y2)[0]!;
          points.push(this.point([x, y, one], { check: true }));
          points.push(this.point([x, y.neg() as C, one], { check: true }));
        }
      }
    } else {
      const b = h.neg().scalar_mul(inv_of(K, two));
      const D = b.mul(b).add(f);
      for (const x of iterate_field(K)) {
        const Dval = D.evaluate(x);
        if (Dval.isZero()) {
          points.push(this.point([x, b.evaluate(x), one], { check: true }));
        } else if (is_square_of(K, Dval)) {
          const sqrtD = sqrt_all_of(K, Dval)[0]!;
          const v = b.evaluate(x);
          points.push(this.point([x, v.add(sqrtD) as C, one], { check: true }));
          points.push(this.point([x, v.sub(sqrtD) as C, one], { check: true }));
        }
      }
    }
    return points;
  }

  /**
   * List points, caching all square roots ahead of time by squaring every
   * element of the field.
   *
   * Port of `hyperelliptic_finite_field.py:754-839`.  Sage builds
   * `square_roots[x*x] = x` while `x` runs over `K` in increasing order, so
   * the *last* (largest) root wins; this is reproduced exactly, because it
   * determines the order in which the two points over an `x`-coordinate are
   * listed.
   */
  _points_cache_sqrt(options?: { brute_force?: boolean }): HyperellipticPoint<C>[] {
    const brute_force = options?.brute_force ?? false;
    const K = this.base_ring();
    const char = characteristic_of(K);
    const [f, h] = this.hyperelliptic_polynomials();
    const one = K.one() as C;
    const zero = K.zero() as C;
    const two = K.__call__(2n) as C;

    let square_roots: Map<string, C> | null = null;
    if (char !== 2n) {
      square_roots = new Map<string, C>();
      for (const x of iterate_field(K)) {
        square_roots.set(key_of(x.mul(x) as C), x);
      }
    }

    const points: HyperellipticPoint<C>[] = [];

    if (this.defining_polynomial_eval(zero, one, zero).isZero()) {
      points.push(this.point([zero, one, zero], { check: true }));
    }
    const degP = this.defining_polynomial_degree();
    if (degP > 2) {
      const s = this.defining_polynomial_eval(one, zero, zero);
      const r = this.defining_polynomial_eval(one, one, zero).sub(s) as C;
      if (!r.isZero()) {
        points.push(this.point([one, div_elements(s.neg() as C, r), zero], { check: true }));
      }
    } else if (char === 2n) {
      for (const y of iterate_field(K)) {
        if (this.defining_polynomial_eval(one, y, zero).isZero()) {
          points.push(this.point([one, y, zero], { check: true }));
        }
      }
    } else {
      const s = f.getCoeff(2).neg() as C;
      const r = h.getCoeff(1);
      const d = div_elements(r.mul(r) as C, two.mul(two) as C).sub(s) as C;
      const sqrtd = square_roots!.get(key_of(d));
      const halfR = div_elements(r.neg() as C, two);
      if (d.isZero()) {
        points.push(this.point([one, halfR, zero], { check: true }));
      } else if (sqrtd !== undefined) {
        points.push(this.point([one, halfR.add(sqrtd) as C, zero], { check: true }));
        points.push(this.point([one, halfR.sub(sqrtd) as C, zero], { check: true }));
      }
    }

    if (char === 2n || brute_force) {
      for (const x of iterate_field(K)) {
        for (const y of iterate_field(K)) {
          if (this.defining_polynomial_eval(x, y, one).isZero()) {
            points.push(this.point([x, y, one], { check: true }));
          }
        }
      }
    } else if (h.isZero()) {
      for (const x of iterate_field(K)) {
        const y2 = f.evaluate(x);
        const y = square_roots!.get(key_of(y2));
        if (y2.isZero()) {
          points.push(this.point([x, y2, one], { check: true }));
        } else if (y !== undefined) {
          points.push(this.point([x, y, one], { check: true }));
          points.push(this.point([x, y.neg() as C, one], { check: true }));
        }
      }
    } else {
      const b = h.neg().scalar_mul(inv_of(K, two));
      const D = b.mul(b).add(f);
      for (const x of iterate_field(K)) {
        const Dval = D.evaluate(x);
        const sqrtD = square_roots!.get(key_of(Dval));
        if (Dval.isZero()) {
          points.push(this.point([x, b.evaluate(x), one], { check: true }));
        } else if (sqrtD !== undefined) {
          const v = b.evaluate(x);
          points.push(this.point([x, v.add(sqrtD) as C, one], { check: true }));
          points.push(this.point([x, v.sub(sqrtD) as C, one], { check: true }));
        }
      }
    }
    return points;
  }

  /**
   * All the points on the plane projective model of this curve.
   *
   * Port of `hyperelliptic_finite_field.py:841-909`.
   */
  points(): HyperellipticPoint<C>[] {
    if (this.__points !== null) {
      return this.__points;
    }
    const K = this.base_ring();
    if (degree_of(K) === 1) {
      this.__points = this._points_cache_sqrt();
    } else {
      const order = cardinality_of(K)!;
      // `zech_log_bound = 2**16` (finite_field_constructor.py:885)
      this.__points = order < 65536n ? this._points_fast_sqrt() : this._points_cache_sqrt();
    }
    return this.__points;
  }

  /* ------------------------------------------------------------------ */
  /* Counting points                                                    */
  /* ------------------------------------------------------------------ */

  /** `hyperelliptic_finite_field.py:911-958` — requires the `hypellfrob` library. */
  count_points_matrix_traces(_n = 1, _M?: unknown, _N?: number): never {
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: count_points_matrix_traces (the hypellfrob library is not ported)'
    );
  }

  /** `hyperelliptic_finite_field.py:1071-1156` — requires the `hypellfrob` library. */
  count_points_hypellfrob(_n = 1, _options?: { N?: number; algorithm?: string }): never {
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: count_points_hypellfrob (the hypellfrob library is not ported)'
    );
  }

  /** `hyperelliptic_finite_field.py:1355-1378` — requires the `hypellfrob` library. */
  cardinality_hypellfrob(_extension_degree = 1, _options?: { algorithm?: string }): never {
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: cardinality_hypellfrob (the hypellfrob library is not ported)'
    );
  }

  /**
   * Count the number of points over the first `n` extensions of the base field
   * from the Frobenius polynomial.
   *
   * Port of `hyperelliptic_finite_field.py:960-1011`.  Sage takes the
   * logarithm of the reciprocal polynomial as a power series over `QQ` and
   * reads off `N_i = q^i + 1 + i * flog[i]`.  Writing
   * `frev(t) = prod_j (1 - a_j t) = sum_k c_k t^k` we have
   * `i * flog[i] = -p_i` with `p_i = sum_j a_j^i`, and Newton's identity
   * `p_n = -n c_n - sum_{j<n} c_{n-j} p_j` computes the same values entirely
   * in `ZZ`.
   */
  count_points_frobenius_polynomial(n = 1, f?: ZZPoly): bigint[] {
    const P = f ?? this.frobenius_polynomial();
    const q = cardinality_of(this.base_ring());
    if (q === null) {
      throw new ValueError('base ring must be finite');
    }

    // frev: coefficient of t^k is P[deg - k]
    let deg = P.length - 1;
    while (deg >= 0 && P[deg] === 0n) {
      deg--;
    }
    const c: bigint[] = [];
    for (let k = 0; k <= deg; k++) {
      c.push(P[deg - k]!);
    }
    if (c[0] !== 1n) {
      throw new ValueError('the reciprocal Frobenius polynomial must have constant term 1');
    }

    const cAt = (k: number): bigint => (k < c.length ? c[k]! : 0n);

    const p: bigint[] = [0n]; // p[0] unused
    for (let k = 1; k <= n; k++) {
      let acc = -BigInt(k) * cAt(k);
      for (let j = 1; j < k; j++) {
        acc -= cAt(k - j) * p[j]!;
      }
      p.push(acc);
    }

    const out: bigint[] = [];
    for (let i = 1; i <= n; i++) {
      out.push(q ** BigInt(i) + 1n - p[i]!);
    }
    return out;
  }

  /**
   * Count the number of points over the first `n` extensions by exhaustive
   * search.
   *
   * Port of `hyperelliptic_finite_field.py:1013-1069`.
   */
  count_points_exhaustive(n = 1, options?: { naive?: boolean }): bigint[] {
    const naive = options?.naive ?? false;
    const g = this.genus();
    const a: bigint[] = [];
    for (let i = 1; i <= Math.min(n, g); i++) {
      a.push(this.cardinality_exhaustive(i));
    }
    if (n <= g) {
      return a;
    }
    if (naive) {
      for (let i = g + 1; i <= n; i++) {
        a.push(this.cardinality_exhaustive(i));
      }
      return a;
    }
    const f = this._frobenius_polynomial_cardinalities(a);
    return this.count_points_frobenius_polynomial(n, f);
  }

  /**
   * Count points over finite fields.
   *
   * Port of `hyperelliptic_finite_field.py:1158-1240`.
   *
   * @see Deviation: upstream prefers `hypellfrob` when the base field is prime
   * and large enough; that library is not ported, so exhaustive search (plus
   * Newton's identities beyond degree `g`) is always used.  The values agree.
   */
  count_points(n = 1): bigint[] {
    return this.count_points_exhaustive(n);
  }

  /**
   * Count points on a single extension of the base field by enumerating over
   * `x` and solving the resulting quadratic equation for `y`.
   *
   * Port of `hyperelliptic_finite_field.py:1242-1353`.
   */
  cardinality_exhaustive(extension_degree = 1): bigint {
    const K = this.base_ring();
    const g = this.genus();
    const n = extension_degree;
    const q = cardinality_of(K);
    if (q === null) {
      throw new ValueError('base ring must be finite');
    }

    if (g === 0) {
      // here is the projective line
      return q ** BigInt(n) + 1n;
    }

    const [f, h] = this.hyperelliptic_polynomials();
    let a = 0n;

    // Work in L = GF(q^n).
    let L: HyperellipticBaseRing<RingElement>;
    let fext: Polynomial<RingElement>;
    let hext: Polynomial<RingElement>;
    if (n === 1) {
      L = K as HyperellipticBaseRing<RingElement>;
      fext = f as unknown as Polynomial<RingElement>;
      hext = h as unknown as Polynomial<RingElement>;
    } else {
      const p = characteristic_of(K);
      const m = degree_of(K);
      // L = GF(K.cardinality()^n) = GF(p^(m*n))  (`hyperelliptic_finite_field.py:1310`)
      const ext = GFpn(p, m * n, undefined, 'z');
      L = ext as unknown as HyperellipticBaseRing<RingElement>;
      const P = new PolynomialRing<RingElement>(L, 't');
      const emb = field_embedding<C>(K, L);
      fext = P.__call__(coeff_list(f, K).map(emb));
      hext = P.__call__(coeff_list(h, K).map(emb));
    }

    // We solve equations of the form y^2 + r*y - s == 0.
    if (characteristic_of(K) === 2n) {
      // points at infinity: y^2 + h[g+1] y = f[2g+2]
      const r = h.getCoeff(g + 1);
      if (r.isZero()) {
        a += 1n;
      } else if (
        n % 2 === 0 ||
        absolute_trace_is_zero(K, div_elements(f.getCoeff(2 * g + 2), r.mul(r) as C))
      ) {
        a += 2n;
      }
      // affine points
      for (const x of iterate_field(L)) {
        const rr = hext.evaluate(x);
        if (rr.isZero()) {
          a += 1n;
        } else if (absolute_trace_is_zero(L, div_elements(fext.evaluate(x), rr.mul(rr)))) {
          a += 2n;
        }
      }
    } else {
      // points at infinity: d = h[g+1]^2 + 4 f[2g+2]
      const hg1 = h.getCoeff(g + 1);
      const d = hg1.mul(hg1).add((K.__call__(4n) as C).mul(f.getCoeff(2 * g + 2))) as C;
      if (d.isZero()) {
        a += 1n;
      } else if (n % 2 === 0 || is_square_of(K, d)) {
        a += 2n;
      }
      // affine points
      const four = L.__call__(4n);
      for (const x of iterate_field(L)) {
        const hx = hext.evaluate(x);
        const dv = hx.mul(hx).add(four.mul(fext.evaluate(x)));
        if (dv.isZero()) {
          a += 1n;
        } else if (is_square_of(L, dv)) {
          a += 2n;
        }
      }
    }

    return a;
  }

  /**
   * Count points on a single extension of the base field.
   *
   * Port of `hyperelliptic_finite_field.py:1380-1433`.
   *
   * @see Deviation: `hypellfrob` is not ported, so exhaustive search is always
   * used.  The values agree with SageMath.
   */
  cardinality(extension_degree = 1): bigint {
    const cached = this.__cardinality.get(extension_degree);
    if (cached !== undefined) {
      return cached;
    }
    const value = this.cardinality_exhaustive(extension_degree);
    this.__cardinality.set(extension_degree, value);
    return value;
  }

  /* ------------------------------------------------------------------ */
  /* Cartier / Hasse-Witt                                               */
  /* ------------------------------------------------------------------ */

  /**
   * Port of `hyperelliptic_finite_field.py:1475-1608` (`_Cartier_matrix_cached`).
   */
  _Cartier_matrix_cached(): CartierData<C> {
    if (this.__cartier !== null) {
      return this.__cartier;
    }
    const Fq = this.base_ring();
    const p = characteristic_of(Fq);

    if (p === 2n) {
      throw new ValueError('p must be odd');
    }

    const g = this.genus();
    const [f, h] = this.hyperelliptic_polynomials();
    if (!h.isZero()) {
      throw new ValueError('E must be of the form y^2 = f(x)');
    }
    const d = f.degree();
    if (d % 2 === 0) {
      throw new ValueError('In this implementation the degree of f must be odd');
    }
    const df = f.derivative();
    if (df.resultant(f).isZero()) {
      throw new ValueError('curve is not smooth');
    }

    // F = f^((p-1)/2)
    const F = f.pow((p - 1n) / 2n);
    const Coeff: C[] = coeff_list(F, Fq);
    const pn = Number(p);
    while (Coeff.length < pn * g) {
      Coeff.push(Fq.zero() as C);
    }

    const M: C[][] = [];
    for (let j = 1; j <= g; j++) {
      const H: C[] = [];
      for (let i = pn * j - 1; i > pn * j - g - 1; i--) {
        H.push(i >= 0 && i < Coeff.length ? Coeff[i]! : (Fq.zero() as C));
      }
      M.push(H);
    }
    this.__cartier = { M, Coeff, g, Fq, p };
    return this.__cartier;
  }

  /**
   * The Cartier matrix `M = (c_{p i - j})` where
   * `f(x)^{(p-1)/2} = sum c_i x^i`.
   *
   * Port of `hyperelliptic_finite_field.py:1611-1701`.
   */
  Cartier_matrix(): C[][] {
    return this._Cartier_matrix_cached().M;
  }

  /**
   * The Hasse-Witt matrix `N = M M^p ... M^{p^{g-1}}`.
   *
   * Port of `hyperelliptic_finite_field.py:1703-1868`.
   */
  Hasse_Witt(): C[][] {
    const { M, Coeff, g, Fq, p } = this._Cartier_matrix_cached();
    const pn = Number(p);

    const frob_mat = (k: number): C[][] => {
      const aExp = p ** BigInt(k);
      const CoeffPow = Coeff.map((c) => pow_element(c, aExp));
      const mat: C[][] = [];
      for (let i = 1; i <= g; i++) {
        const H: C[] = [];
        for (let j = pn * i - 1; j > pn * i - g - 1; j--) {
          H.push(j >= 0 && j < CoeffPow.length ? CoeffPow[j]! : (Fq.zero() as C));
        }
        mat.push(H);
      }
      return mat;
    };

    const Mall: C[][][] = [M];
    for (let k = 1; k < g; k++) {
      Mall.push(frob_mat(k));
    }
    Mall.reverse();

    let N = identity_matrix(Fq, g) as C[][];
    for (const l of Mall) {
      N = matrix_mul(Fq, N, l);
    }
    return N;
  }

  /** `hyperelliptic_finite_field.py:1870-1908` */
  a_number(): number {
    const { M, g, Fq } = this._Cartier_matrix_cached();
    return g - matrix_rank(Fq, M);
  }

  /** `hyperelliptic_finite_field.py:1910-1947` */
  p_rank(): number {
    const N = this.Hasse_Witt();
    return matrix_rank(this.base_ring(), N);
  }
}

/** Cached data returned by `_Cartier_matrix_cached`. */
export interface CartierData<C extends RingElement> {
  M: C[][];
  Coeff: C[];
  g: number;
  Fq: HyperellipticBaseRing<C>;
  p: bigint;
}

/* -------------------------------------------------------------------------- */
/* Helpers                                                                    */
/* -------------------------------------------------------------------------- */

/**
 * `B` such that `p^B >= M` where `M = sqrt(M2)`, matching Sage's
 * `B = ZZ(M.ceil()).exact_log(p); if p**B < M: B += 1`.
 */
function bound_from_square(M2: bigint, p: bigint): number {
  const r = isqrt(M2);
  const Mceil = r * r === M2 ? r : r + 1n;
  // exact_log: largest B with p^B <= Mceil
  let B = 0;
  let pw = 1n;
  while (pw * p <= Mceil) {
    pw *= p;
    B += 1;
  }
  // if p**B < M  <=>  p**(2B) < M2
  if (pw * pw < M2) {
    B += 1;
  }
  return B;
}

/** `f.list()`: the coefficients of `f` in ascending order (length `deg+1`). */
function coeff_list<C extends RingElement>(f: Polynomial<C>, K: HyperellipticBaseRing<C>): C[] {
  const out: C[] = [];
  for (let i = 0; i <= f.degree(); i++) {
    out.push(f.getCoeff(i));
  }
  if (out.length === 0) {
    out.push(K.zero() as C);
  }
  return out;
}

/** A hashable key for a base-ring element. */
function key_of<C extends RingElement>(a: C): string {
  const v = element_to_bigint(a);
  return v === null ? a.toString() : v.toString();
}

function inv_of<C extends RingElement>(K: HyperellipticBaseRing<C>, a: C): C {
  return div_elements(K.one() as C, a);
}

function identity_matrix<C extends RingElement>(K: HyperellipticBaseRing<C>, n: number): C[][] {
  const out: C[][] = [];
  for (let i = 0; i < n; i++) {
    const row: C[] = [];
    for (let j = 0; j < n; j++) {
      row.push((i === j ? K.one() : K.zero()) as C);
    }
    out.push(row);
  }
  return out;
}

function matrix_mul<C extends RingElement>(K: HyperellipticBaseRing<C>, A: C[][], B: C[][]): C[][] {
  const n = A.length;
  const m = B[0]?.length ?? 0;
  const k = B.length;
  const out: C[][] = [];
  for (let i = 0; i < n; i++) {
    const row: C[] = [];
    for (let j = 0; j < m; j++) {
      let s = K.zero() as C;
      for (let t = 0; t < k; t++) {
        s = s.add(A[i]![t]!.mul(B[t]![j]!) as C) as C;
      }
      row.push(s);
    }
    out.push(row);
  }
  return out;
}

/** Rank of a matrix over a field, by Gaussian elimination. */
function matrix_rank<C extends RingElement>(K: HyperellipticBaseRing<C>, A: C[][]): number {
  const m = A.map((row) => row.slice());
  const rows = m.length;
  const cols = rows === 0 ? 0 : m[0]!.length;
  let rank = 0;
  let row = 0;
  for (let col = 0; col < cols && row < rows; col++) {
    let pivot = -1;
    for (let i = row; i < rows; i++) {
      if (!m[i]![col]!.isZero()) {
        pivot = i;
        break;
      }
    }
    if (pivot === -1) {
      continue;
    }
    const tmp = m[row]!;
    m[row] = m[pivot]!;
    m[pivot] = tmp;
    const inv = inv_of(K, m[row]![col]!);
    for (let j = col; j < cols; j++) {
      m[row]![j] = m[row]![j]!.mul(inv) as C;
    }
    for (let i = 0; i < rows; i++) {
      if (i === row) {
        continue;
      }
      const factor = m[i]![col]!;
      if (factor.isZero()) {
        continue;
      }
      for (let j = col; j < cols; j++) {
        m[i]![j] = m[i]![j]!.sub(factor.mul(m[row]![j]!) as C) as C;
      }
    }
    row += 1;
    rank += 1;
  }
  return rank;
}
