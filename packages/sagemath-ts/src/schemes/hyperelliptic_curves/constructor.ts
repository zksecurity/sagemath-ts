/**
 * @module sage/schemes/hyperelliptic_curves/constructor
 * @description Hyperelliptic curve constructor
 *
 * Port of: `sage/schemes/hyperelliptic_curves/constructor.py`
 *
 * AUTHORS (upstream):
 *
 * - David Kohel (2006): initial version
 * - Anna Somoza (2019-04): dynamic class creation
 */

import { NotImplementedError, ValueError } from '../../errors.js';
import type { Polynomial, RingElement } from '../../rings/polynomial/polynomial_element.js';
import { cardinality_of, characteristic_of } from './field_ops.js';
import { HyperellipticCurve_finite_field } from './hyperelliptic_finite_field.js';
import {
  HyperellipticCurve_g2,
  HyperellipticCurve_g2_FiniteField,
  HyperellipticCurve_g2_RationalField,
} from './hyperelliptic_g2.js';
import {
  HyperellipticCurve_generic,
  _register_constructor_module,
  floorDiv,
  poly_repr,
} from './hyperelliptic_generic.js';
import { HyperellipticCurve_rational_field } from './hyperelliptic_rational_field.js';
// Loaded for its side effect: it registers itself so that `curve.jacobian()`
// works for any curve built through this constructor.
import './jacobian_generic.js';

export interface HyperellipticCurveOptions {
  /** Names of the coordinate functions; default `["x", "y"]`. */
  names?: string | [string, string];
  /**
   * Test that the input defines a hyperelliptic curve (default `true`).
   *
   * .. WARNING:: with `check_squarefree = false` the output curve is not to be
   * trusted; in particular `is_singular` always returns `false`.
   */
  check_squarefree?: boolean;
}

/**
 * Parse a defining equation `g(x, y) = y^2 + h(x) y - f(x)` given as a
 * multivariate polynomial.
 *
 * Port target: `constructor.py:32-91`.  Not ported: this port has no
 * `MPolynomial` type that the hyperelliptic module can consume.
 */
export function _parse_multivariate_defining_equation(_g: unknown): never {
  throw new NotImplementedError(
    'SAGE_NOT_IMPLEMENTED: _parse_multivariate_defining_equation ' +
      '(defining a hyperelliptic curve from a bivariate polynomial requires MPolynomial support)'
  );
}

/**
 * Return the hyperelliptic curve `y^2 + h y = f`, for univariate polynomials
 * `h` and `f`.  If `h` is not given, then it defaults to 0.
 *
 * Port of `constructor.py:94-368`.
 *
 * @example
 * ```typescript
 * const K = GF(37n);
 * const R = new PolynomialRing(K, 'x');
 * const x = R.gen();
 * const H = HyperellipticCurve(x.pow(5).add(x).add(R.__call__(2n)));
 * ```
 */
export function HyperellipticCurve<C extends RingElement>(
  f: Polynomial<C>,
  h?: Polynomial<C> | null,
  options?: HyperellipticCurveOptions
): HyperellipticCurve_generic<C> {
  const check_squarefree = options?.check_squarefree ?? true;
  const names = normalize_names(options?.names);

  const P = f.parent;
  const hh: Polynomial<C> = h ?? P.zero();
  const K = P.base_ring;

  // F is the discriminant; use this for the type check rather than f and h,
  // one of which might be constant (`constructor.py:269-286`).
  const four = P.__call__(K.__call__(4n) as C);
  const F = hh.mul(hh).add(four.mul(f));

  const df = f.degree();
  const dh_2 = 2 * hh.degree();
  const g = dh_2 < df ? floorDiv(df - 1, 2) : floorDiv(dh_2 - 1, 2);

  if (check_squarefree) {
    // Assuming we are working over a field, this checks that after resolving
    // the singularity at infinity, we get a smooth double cover of P^1.
    let should_be_coprime: [Polynomial<C>, Polynomial<C>];
    if (characteristic_of(K as never) === 2n) {
      // characteristic 2
      if (hh.isZero()) {
        throw new ValueError(
          `for characteristic 2, argument h = ${poly_repr(hh, names[0])} must be nonzero`
        );
      }
      const hg1 = hh.getCoeff(g + 1);
      if (
        hg1.isZero() &&
        f
          .getCoeff(2 * g + 1)
          .mul(f.getCoeff(2 * g + 1))
          .eq(
            f
              .getCoeff(2 * g + 2)
              .mul(hh.getCoeff(g))
              .mul(hh.getCoeff(g)) as C
          )
      ) {
        throw new ValueError('not a hyperelliptic curve: highly singular at infinity');
      }
      const hp = hh.derivative();
      const fp = f.derivative();
      should_be_coprime = [hh, f.mul(hp).mul(hp).add(fp.mul(fp))];
    } else {
      // characteristic not 2
      const dF = F.degree();
      if (dF !== 2 * g + 1 && dF !== 2 * g + 2) {
        throw new ValueError('not a hyperelliptic curve: highly singular at infinity');
      }
      should_be_coprime = [F, F.derivative()];
    }

    let smooth: boolean;
    try {
      smooth = should_be_coprime[0].gcd(should_be_coprime[1]).degree() === 0;
    } catch (_e) {
      try {
        smooth = !should_be_coprime[0].resultant(should_be_coprime[1]).isZero();
      } catch (_e2) {
        throw new NotImplementedError(
          `cannot determine whether polynomials [${poly_repr(
            should_be_coprime[0],
            names[0]
          )}, ${poly_repr(
            should_be_coprime[1],
            names[0]
          )}] have a common root, use check_squarefree=False to skip this check`
        );
      }
    }
    if (!smooth) {
      throw new ValueError('not a hyperelliptic curve: singularity in the provided affine patch');
    }
  }

  // Specialise to subclasses, mirroring `constructor.py:335-368`.  TypeScript
  // has no multiple inheritance, so the four combinations of
  // {generic, g2} x {generic, FiniteField, RationalField} are explicit classes.
  const isFinite = cardinality_of(K as never) !== null;
  const isRational = !isFinite && characteristic_of(K as never) === 0n && ring_is_QQ(K);

  if (g === 2) {
    if (isFinite) {
      return new HyperellipticCurve_g2_FiniteField<C>(f, hh, names, g);
    }
    if (isRational) {
      return new HyperellipticCurve_g2_RationalField<C>(f, hh, names, g);
    }
    return new HyperellipticCurve_g2<C>(f, hh, names, g);
  }
  if (isFinite) {
    return new HyperellipticCurve_finite_field<C>(f, hh, names, g);
  }
  if (isRational) {
    return new HyperellipticCurve_rational_field<C>(f, hh, names, g);
  }
  return new HyperellipticCurve_generic<C>(f, hh, names, g);
}

/** Sage's `normalize_names(2, names)`. */
function normalize_names(names?: string | [string, string]): [string, string] {
  if (names === undefined || names === null) {
    return ['x', 'y'];
  }
  if (typeof names === 'string') {
    const parts = names
      .split(',')
      .map((s) => s.trim())
      .filter((s) => s.length > 0);
    if (parts.length !== 2) {
      throw new ValueError(`variable names must be a list of 2 names, got ${names}`);
    }
    return [parts[0]!, parts[1]!];
  }
  return names;
}

/** Detect Sage's `RationalField` structurally. */
function ring_is_QQ(K: unknown): boolean {
  return String(K) === 'Rational Field';
}

_register_constructor_module({ HyperellipticCurve });
