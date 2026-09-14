/**
 * @module sage/schemes/elliptic_curves/weierstrass_morphism
 * @description Isomorphisms between Weierstrass models of elliptic curves
 *
 * Port of: sage/schemes/elliptic_curves/weierstrass_morphism.py
 *
 * Weierstrass isomorphisms are specified by lists of the form [u, r, s, t]
 * (with u != 0) which specifies a transformation (x, y) -> (x', y') where
 *
 *   (x, y) = (u^2 * x' + r, u^3 * y' + s * u^2 * x' + t)
 *
 * @see Reference: sage/schemes/elliptic_curves/weierstrass_morphism.py
 */

import { NotImplementedError, TypeError, ValueError } from '../../errors.js';
import { cmp_universal } from '@sagemath-ts/parigp-ts/src/gen2.js';
import { PariType } from '@sagemath-ts/parigp-ts/src/types.js';
import {
  FiniteFieldElement as ExtensionElement,
  FiniteFieldExtension,
  PrimeField,
} from '../../rings/finite_rings/finite_field_extension.js';
import { FiniteFieldPrime } from '../../rings/finite_rings/finite_field_prime.js';
import { GF2Element, GF2Field } from '../../rings/finite_rings/gf2.js';
import { QQ } from '../../rings/rational_field.js';
import { PolynomialRing } from '../../rings/polynomial/polynomial_ring.js';
import type { CoefficientRing, RingElement } from '../../rings/polynomial/polynomial_element.js';
import { EllipticCurve } from './constructor.js';
import { EllipticCurveGeneric } from './ell_generic.js';
import type { EllipticCurvePoint, FieldElement } from './ell_point.js';
import { _same_base_ring } from './types.js';

/** Scalar ordering used by the native morphism sorting key. */
function _compare_isomorphism_parameters(a: FieldElement, b: FieldElement): number {
  if (a.parent !== b.parent && a.parent && b.parent && _same_base_ring(a.parent, b.parent)) {
    b = a.parent.__call__(b);
  }
  const cmp = (a as unknown as { cmp?: (other: unknown) => number }).cmp;
  if (typeof cmp === 'function') return cmp.call(a, b);
  if (a instanceof ExtensionElement && b instanceof ExtensionElement) {
    return cmp_universal(
      {
        type: PariType.t_FFELT,
        p: a.parent.characteristic,
        degree: a.parent.degree,
        value: a.coefficients().map((c) => c.value),
        definingPoly: a.parent.modulus.coeffs.map((c) => c.value),
      },
      {
        type: PariType.t_FFELT,
        p: b.parent.characteristic,
        degree: b.parent.degree,
        value: b.coefficients().map((c) => c.value),
        definingPoly: b.parent.modulus.coeffs.map((c) => c.value),
      }
    );
  }
  const av = (a as unknown as { value?: unknown }).value;
  const bv = (b as unknown as { value?: unknown }).value;
  if (typeof av === 'bigint' && typeof bv === 'bigint') return av < bv ? -1 : av > bv ? 1 : 0;
  if (a instanceof GF2Element && b instanceof GF2Element) return a.value - b.value;
  if (a.eq(b)) return 0;
  throw new NotImplementedError(
    'SAGE_NOT_IMPLEMENTED: ordering isomorphism parameters over this field'
  );
}

/**
 * This class implements the basic arithmetic of isomorphisms between
 * Weierstrass models of elliptic curves.
 *
 * These are specified by lists of the form [u, r, s, t] (with u != 0)
 * which specifies a transformation (x, y) -> (x', y') where
 *
 *   (x, y) = (u^2 * x' + r, u^3 * y' + s * u^2 * x' + t)
 *
 * @see Reference: sage/schemes/elliptic_curves/weierstrass_morphism.py:baseWI
 */
export class baseWI<F extends FieldElement = FieldElement> {
  /** The u parameter */
  u: F;
  /** The r parameter */
  r: F;
  /** The s parameter */
  s: F;
  /** The t parameter */
  t: F;

  /**
   * Construct a baseWI with given parameters (defaults to identity).
   *
   * @param u - The u parameter (default: 1)
   * @param r - The r parameter (default: 0)
   * @param s - The s parameter (default: 0)
   * @param t - The t parameter (default: 0)
   * @see Reference: sage/schemes/elliptic_curves/weierstrass_morphism.py:baseWI.__init__
   */
  constructor(u: F, r: F, s: F, t: F) {
    if (u.isZero()) {
      throw new ValueError('u!=0 required for baseWI');
    }
    this.u = u;
    this.r = r;
    this.s = s;
    this.t = t;
  }

  /**
   * Return the parameters u, r, s, t as a tuple.
   *
   * @returns The parameters [u, r, s, t]
   * @see Reference: sage/schemes/elliptic_curves/weierstrass_morphism.py:baseWI.tuple
   */
  tuple(): [F, F, F, F] {
    return [this.u, this.r, this.s, this.t];
  }

  /**
   * Return the composition of this isomorphism and another.
   *
   * @param other - Another baseWI
   * @returns The composition self * other
   * @see Reference: sage/schemes/elliptic_curves/weierstrass_morphism.py:baseWI.__mul__
   */
  mul(other: baseWI<F>): baseWI<F> {
    const [u1, r1, s1, t1] = other.tuple();
    const [u2, r2, s2, t2] = this.tuple();

    // u' = u1 * u2
    const uNew = u1.mul(u2) as F;

    // r' = u1^2 * r2 + r1
    const rNew = u1.mul(u1).mul(r2).add(r1) as F;

    // s' = u1 * s2 + s1
    const sNew = u1.mul(s2).add(s1) as F;

    // t' = u1^3 * t2 + s1 * u1^2 * r2 + t1
    const u1Sq = u1.mul(u1) as F;
    const u1Cubed = u1Sq.mul(u1) as F;
    const tNew = u1Cubed.mul(t2).add(s1.mul(u1Sq).mul(r2)).add(t1) as F;

    return new baseWI<F>(uNew, rNew, sNew, tNew);
  }

  /**
   * Return the inverse of this isomorphism.
   *
   * @returns The inverse baseWI
   * @see Reference: sage/schemes/elliptic_curves/weierstrass_morphism.py:baseWI.__invert__
   */
  invert(): baseWI<F> {
    const [u, r, s, t] = this.tuple();
    const uInv = u.inv() as F;
    const uSq = u.mul(u) as F;
    const uCubed = uSq.mul(u) as F;

    // u' = 1/u
    const uNew = uInv;

    // r' = -r/u^2
    const rNew = r.neg().div(uSq) as F;

    // s' = -s/u
    const sNew = s.neg().div(u) as F;

    // t' = (r*s - t)/u^3
    const tNew = r.mul(s).sub(t).div(uCubed) as F;

    return new baseWI<F>(uNew, rNew, sNew, tNew);
  }

  /**
   * Return True if this is the identity isomorphism.
   *
   * @returns True if (u, r, s, t) == (1, 0, 0, 0)
   * @see Reference: sage/schemes/elliptic_curves/weierstrass_morphism.py:baseWI.is_identity
   */
  is_identity(): boolean {
    const K = this.u.parent;
    const one = K.one();
    return this.u.eq(one) && this.r.isZero() && this.s.isZero() && this.t.isZero();
  }

  /**
   * Base application of isomorphisms to curves and points.
   *
   * A baseWI w may be applied to a list [a1, a2, a3, a4, a6] representing
   * the a-invariants of an elliptic curve E, returning the a-invariants
   * of w(E); or to P=[x,y] or P=[x,y,z] representing a point in A^2 or P^2,
   * returning the transformed point.
   *
   * @param EorP - Either a-invariants or point coordinates
   * @returns The transformed curve or point
   * @see Reference: sage/schemes/elliptic_curves/weierstrass_morphism.py:baseWI.__call__
   */
  call(EorP: F[]): F[] {
    const [u, r, s, t] = this.tuple();

    if (EorP.length === 5) {
      // Transform a-invariants [a1, a2, a3, a4, a6]
      let [a1, a2, a3, a4, a6] = EorP;

      // a6' = a6 + r*(a4 + r*(a2 + r)) - t*(a3 + r*a1 + t)
      a6 = a6.add(r.mul(a4.add(r.mul(a2.add(r))))).sub(t.mul(a3.add(r.mul(a1)).add(t))) as F;

      // a4' = a4 - s*a3 + 2*r*a2 - (t + r*s)*a1 + 3*r^2 - 2*s*t
      // `Rational` deliberately has no `parent` property in this port.
      // Build the small constants from the nonzero isomorphism scale instead
      // of assuming every field element carries its parent.
      const one = u.div(u) as F;
      const two = one.add(one) as F;
      const three = two.add(one) as F;
      a4 = a4
        .sub(s.mul(a3))
        .add(two.mul(r).mul(a2))
        .sub(t.add(r.mul(s)).mul(a1))
        .add(three.mul(r).mul(r))
        .sub(two.mul(s).mul(t)) as F;

      // a3' = a3 + r*a1 + 2*t
      a3 = a3.add(r.mul(a1)).add(two.mul(t)) as F;

      // a2' = a2 - s*a1 + 3*r - s^2
      a2 = a2.sub(s.mul(a1)).add(three.mul(r)).sub(s.mul(s)) as F;

      // a1' = a1 + 2*s
      a1 = a1.add(two.mul(s)) as F;

      // Divide by appropriate powers of u
      const uSq = u.mul(u) as F;
      const uCubed = uSq.mul(u) as F;
      const u4 = uSq.mul(uSq) as F;
      const u6 = uCubed.mul(uCubed) as F;

      return [
        a1.div(u) as F,
        a2.div(uSq) as F,
        a3.div(uCubed) as F,
        a4.div(u4) as F,
        a6.div(u6) as F,
      ];
    }

    if (EorP.length === 2) {
      // Transform affine point (x, y)
      let [x, y] = EorP;

      // x' = (x - r) / u^2
      x = x.sub(r) as F;

      // y' = (y - s*x - t) / u^3  (note: x is already x - r here)
      // Actually: y' = (y - s*(x-r) - t) / u^3
      y = y.sub(s.mul(x)).sub(t) as F;

      const uSq = u.mul(u) as F;
      const uCubed = uSq.mul(u) as F;

      return [x.div(uSq) as F, y.div(uCubed) as F];
    }

    if (EorP.length === 3) {
      // Transform projective point (x, y, z)
      let [x, y, z] = EorP;

      // x' = (x - r*z) / u^2
      x = x.sub(r.mul(z)) as F;

      // y' = (y - s*x - t*z) / u^3  (note: x is already x - r*z)
      y = y.sub(s.mul(x)).sub(t.mul(z)) as F;

      const uSq = u.mul(u) as F;
      const uCubed = uSq.mul(u) as F;

      return [x.div(uSq) as F, y.div(uCubed) as F, z];
    }

    throw new ValueError('baseWI(a) only for a=(x,y), (x:y:z) or (a1,a2,a3,a4,a6)');
  }

  /**
   * Return the string representation of this isomorphism.
   *
   * @returns String representation
   * @see Reference: sage/schemes/elliptic_curves/weierstrass_morphism.py:baseWI.__repr__
   */
  toString(): string {
    return `(${this.u}, ${this.r}, ${this.s}, ${this.t})`;
  }
}

/**
 * Enumerate all isomorphisms between two elliptic curves, as a generator.
 *
 * @param E - First elliptic curve
 * @param F - Second elliptic curve
 * @returns Generator producing 4-tuples (u, r, s, t) representing an isomorphism
 * @see Deviation: Weierstrass Isomorphism Root Dispatch
 * @see Reference: sage/schemes/elliptic_curves/weierstrass_morphism.py:_isomorphisms
 */
export function* _isomorphisms<FE extends FieldElement>(
  E: EllipticCurveGeneric<FE>,
  F: EllipticCurveGeneric<FE>
): Generator<[FE, FE, FE, FE]> {
  if (!(E instanceof EllipticCurveGeneric) || !(F instanceof EllipticCurveGeneric)) {
    throw new ValueError('arguments are not elliptic curves');
  }
  // Several port classes represent Sage's same cached finite-field parent.
  // Put equivalent-parent coefficients in E's implementation before arithmetic.
  if (E.base_ring !== F.base_ring && _same_base_ring(E.base_ring, F.base_ring)) {
    F = new EllipticCurveGeneric(
      E.base_ring,
      F.a_invariants().map((v) => E.base_ring.__call__(v)) as [FE, FE, FE, FE, FE]
    );
  }
  // Quick check: j-invariants must match
  const jE = E.j_invariant();
  const jF = F.j_invariant();
  if (!jE.eq(jF)) {
    return;
  }

  const K = E.base_ring;
  const char = K.characteristic;
  const R = new PolynomialRing(K as unknown as CoefficientRing<FE & RingElement>, 'x');
  const zero = K.zero();
  const one = K.one();

  const [a1E, a2E, a3E, a4E, a6E] = E.a_invariants();
  const [a1F, a2F, a3F, a4F, a6F] = F.a_invariants();

  // Handle characteristic 2
  if (char === 2n) {
    const j = jE;
    if (j.isZero()) {
      // j = 0 in characteristic 2
      // Find u such that u^3 = a3E/a3F
      const uCubed = a3E.div(a3F) as FE;
      const uList = R.__call__([uCubed.neg(), zero, zero, one]).roots({ multiplicities: false });

      for (const u of uList) {
        const uSq = u.mul(u) as FE;
        const uCub = uSq.mul(u) as FE;

        // s satisfies: s^4 + a3E*s + (a2F^2 + a4F)*u^4 + a2E^2 + a4E = 0
        const sCoeff = uSq.mul(uSq) as FE; // u^4
        const constTerm = a2F.mul(a2F).add(a4F).mul(sCoeff).add(a2E.mul(a2E)).add(a4E) as FE;

        // Find roots of s^4 + a3E*s + constTerm = 0
        const sRoots = R.__call__([constTerm, a3E, zero, zero, one]).roots({
          multiplicities: false,
        });

        for (const s of sRoots) {
          // r = s^2 + a2E + a2F*u^2
          const r = s.mul(s).add(a2E).add(a2F.mul(uSq)) as FE;

          // Find t such that t^2 + a3E*t + r^3 + a2E*r^2 + a4E*r + a6E + a6F*u^6 = 0
          const tConstTerm = r
            .mul(r)
            .mul(r)
            .add(a2E.mul(r).mul(r))
            .add(a4E.mul(r))
            .add(a6E)
            .add(a6F.mul(uCub).mul(uCub)) as FE;
          const tRoots = R.__call__([tConstTerm, a3E, one]).roots({ multiplicities: false });

          for (const t of tRoots) {
            yield [u, r, s, t];
          }
        }
      }
    } else {
      // j != 0 in characteristic 2
      const u = a1E.div(a1F) as FE;
      const uSq = u.mul(u) as FE;
      const uCubed = uSq.mul(u) as FE;

      const r = a3E.add(a3F.mul(uCubed)).div(a1E) as FE;

      // Find s such that s^2 + a1E*s + r + a2E + a2F*u^2 = 0
      const sConstTerm = r.add(a2E).add(a2F.mul(uSq)) as FE;
      const sRoots = R.__call__([sConstTerm, a1E, one]).roots({ multiplicities: false });

      for (const s of sRoots) {
        const u4 = uSq.mul(uSq) as FE;
        const t = a4E
          .add(a4F.mul(u4))
          .add(s.mul(a3E))
          .add(r.mul(s).mul(a1E))
          .add(r.mul(r))
          .div(a1E) as FE;
        yield [u, r, s, t];
      }
    }
    return;
  }

  // Get b-invariants
  const [b2E, b4E, b6E] = E.b_invariants();
  const [b2F, b4F, b6F] = F.b_invariants();

  // Handle characteristic 3
  if (char === 3n) {
    const j = jE;
    if (j.isZero()) {
      // j = 0 in characteristic 3
      // Find u such that u^4 = b4E/b4F
      const u4 = b4E.div(b4F) as FE;
      const uList = R.__call__([u4.neg(), zero, zero, zero, one]).roots({ multiplicities: false });

      for (const u of uList) {
        const uSq = u.mul(u) as FE;
        const uCubed = uSq.mul(u) as FE;

        const s = a1E.sub(a1F.mul(u)) as FE;
        const t = a3E.sub(a3F.mul(uCubed)) as FE;

        // Find r such that r^3 - b4E*r + b6E - b6F*u^6 = 0
        const u6 = uCubed.mul(uCubed) as FE;
        const rConstTerm = b6E.sub(b6F.mul(u6)) as FE;
        const rRoots = R.__call__([rConstTerm, b4E.neg(), zero, one]).roots({
          multiplicities: false,
        });

        for (const r of rRoots) {
          yield [u, r, s, t.add(r.mul(a1E)) as FE];
        }
      }
    } else {
      // j != 0 in characteristic 3
      // Find u such that u^2 = b2E/b2F
      const uSqVal = b2E.div(b2F) as FE;
      const uList = R.__call__([uSqVal.neg(), zero, one]).roots({ multiplicities: false });

      for (const u of uList) {
        const uSq = u.mul(u) as FE;
        const uCubed = uSq.mul(u) as FE;
        const u4 = uSq.mul(uSq) as FE;

        const r = b4F.mul(u4).sub(b4E).div(b2E) as FE;
        const s = a1E.sub(a1F.mul(u)) as FE;
        const t = a3E.sub(a3F.mul(uCubed)).add(a1E.mul(r)) as FE;
        yield [u, r, s, t];
      }
    }
    return;
  }

  // Characteristic not 2 or 3
  const [c4E, c6E] = E.c_invariants();
  const [c4F, c6F] = F.c_invariants();

  const j = jE;
  let m: number;
  let um: FE;

  if (j.isZero()) {
    // j = 0: u^6 = c6E/c6F
    m = 6;
    um = c6E.div(c6F) as FE;
  } else if (j.eq(K.__call__(1728n))) {
    // j = 1728: u^4 = c4E/c4F
    m = 4;
    um = c4E.div(c4F) as FE;
  } else {
    // General case: u^2 = (c6E*c4F)/(c6F*c4E)
    m = 2;
    um = c6E.mul(c4F).div(c6F.mul(c4E)) as FE;
  }

  // Sage delegates every characteristic branch to distinct polynomial roots.
  const coefficients = Array.from({ length: m + 1 }, () => zero);
  coefficients[0] = um.neg();
  coefficients[m] = one;
  const uList = R.__call__(coefficients).roots({ multiplicities: false });

  const two = K.__call__(2n) as FE;
  const three = K.__call__(3n) as FE;

  for (const u of uList) {
    const uSq = u.mul(u) as FE;
    const uCubed = uSq.mul(u) as FE;

    // s = (a1F*u - a1E)/2
    const s = a1F.mul(u).sub(a1E).div(two) as FE;

    // r = (a2F*u^2 + a1E*s + s^2 - a2E)/3
    const r = a2F.mul(uSq).add(a1E.mul(s)).add(s.mul(s)).sub(a2E).div(three) as FE;

    // t = (a3F*u^3 - a1E*r - a3E)/2
    const t = a3F.mul(uCubed).sub(a1E.mul(r)).sub(a3E).div(two) as FE;

    yield [u, r, s, t];
  }
}

/**
 * Class representing a Weierstrass isomorphism between two elliptic curves.
 *
 * Given two Elliptic Curves E and F (represented by Weierstrass models as usual),
 * and a transformation urst from E to F, construct an isomorphism from E to F.
 * An exception is raised if urst(E) != F. At most one of E, F, urst can be None.
 *
 * Explicitly, the isomorphism defined by (u, r, s, t) maps a point (x, y) to
 *
 *   ((x - r) / u^2, (y - s(x-r) - t) / u^3)
 *
 * @see Reference: sage/schemes/elliptic_curves/weierstrass_morphism.py:WeierstrassIsomorphism
 */
export class WeierstrassIsomorphism<F extends FieldElement = FieldElement> extends baseWI<F> {
  private _domain: EllipticCurveGeneric<F>;
  private _codomain: EllipticCurveGeneric<F>;
  private _degree: bigint = 1n;

  /**
   * Construct a WeierstrassIsomorphism.
   *
   * @param E - An EllipticCurve, or null
   * @param urst - A 4-tuple [u, r, s, t], a baseWI object, or null
   * @param codomain - An EllipticCurve, or null
   * @see Reference: sage/schemes/elliptic_curves/weierstrass_morphism.py:WeierstrassIsomorphism.__init__
   */
  constructor(
    E?: EllipticCurveGeneric<F> | null,
    urst?: [F, F, F, F] | baseWI<F> | null,
    codomain?: EllipticCurveGeneric<F> | null
  ) {
    // Check that at most one argument is null/undefined
    const argCount = [E, urst, codomain].filter((x) => x !== null && x !== undefined).length;
    if (argCount < 2) {
      throw new ValueError('at most 1 argument can be None');
    }

    // Get the parameters
    let u: F;
    let r: F;
    let s: F;
    let t: F;

    if (urst instanceof baseWI) {
      [u, r, s, t] = urst.tuple();
    } else if (Array.isArray(urst)) {
      [u, r, s, t] = urst;
    } else if (E !== null && E !== undefined && codomain !== null && codomain !== undefined) {
      // Find isomorphism from E to codomain
      const isoGen = _isomorphisms(E, codomain);
      const first = isoGen.next();
      if (first.done) {
        throw new ValueError('elliptic curves not isomorphic');
      }
      [u, r, s, t] = first.value;
    } else {
      throw new ValueError('invalid parameters to WeierstrassIsomorphism constructor');
    }

    // Call baseWI constructor
    super(u, r, s, t);

    // Compute domain and codomain
    if (E !== null && E !== undefined) {
      this._domain = E;

      if (codomain === null || codomain === undefined) {
        // Compute codomain from urst(E)
        const newAinvs = this.call(E.a_invariants()) as [F, F, F, F, F];
        this._codomain = EllipticCurve(E.base_ring, newAinvs);
      } else {
        // Verify that urst(E) = codomain
        const newAinvs = this.call(E.a_invariants()) as [F, F, F, F, F];
        const computedCodomain = EllipticCurve(E.base_ring, newAinvs);
        if (!_curves_equal(computedCodomain, codomain)) {
          throw new ValueError(
            'second argument is not an isomorphism from first argument to third argument'
          );
        }
        this._codomain = codomain;
      }
    } else {
      // E is null, codomain is given
      // Compute domain by inverting urst
      this._codomain = codomain!;
      const invWI = new baseWI(u, r, s, t).invert();
      const domainAinvs = invWI.call(codomain!.a_invariants()) as [F, F, F, F, F];
      this._domain = EllipticCurve(codomain!.base_ring, domainAinvs);
    }
  }

  /**
   * Return the domain curve of this isomorphism.
   *
   * @returns The domain elliptic curve
   * @see Reference: sage/schemes/elliptic_curves/weierstrass_morphism.py:WeierstrassIsomorphism (domain is inherited from EllipticCurveHom)
   */
  domain(): EllipticCurveGeneric<F> {
    return this._domain;
  }

  /**
   * Return the codomain curve of this isomorphism.
   *
   * @returns The codomain elliptic curve
   * @see Reference: sage/schemes/elliptic_curves/weierstrass_morphism.py:WeierstrassIsomorphism (codomain is inherited from EllipticCurveHom)
   */
  codomain(): EllipticCurveGeneric<F> {
    return this._codomain;
  }

  /**
   * Return the degree of this isomorphism.
   *
   * For isomorphisms, the degree is always 1.
   *
   * @returns Always 1n
   * @see Reference: sage/schemes/elliptic_curves/weierstrass_morphism.py:WeierstrassIsomorphism (degree is always 1)
   */
  degree(): bigint {
    return 1n;
  }

  /**
   * Compare an isomorphism to another elliptic-curve morphism.
   *
   * @param left - Left morphism
   * @param right - Right morphism
   * @param op - Comparison operator
   * @returns Comparison result
   * @see Reference: sage/schemes/elliptic_curves/weierstrass_morphism.py:WeierstrassIsomorphism._comparison_impl
   * @see Deviation: Generic Curve Isomorphism Ordering
   */
  static _comparison_impl<F extends FieldElement>(
    left: WeierstrassIsomorphism<F>,
    right: WeierstrassIsomorphism<F>,
    op: string
  ): boolean;
  static _comparison_impl(left: unknown, right: unknown, op: string): boolean | null;
  static _comparison_impl<F extends FieldElement>(
    left: unknown,
    right: unknown,
    op: string
  ): boolean | null {
    if (!(left instanceof WeierstrassIsomorphism) || !(right instanceof WeierstrassIsomorphism)) {
      // The port uses null for Python's NotImplemented, as in _composition_impl.
      return null;
    }
    const operators: Record<string, string> = {
      eq: '==',
      ne: '!=',
      lt: '<',
      le: '<=',
      gt: '>',
      ge: '>=',
      '==': '==',
      '!=': '!=',
      '<': '<',
      '<=': '<=',
      '>': '>',
      '>=': '>=',
    };
    const symbol = Object.hasOwn(operators, op) ? operators[op] : undefined;
    if (!symbol) throw new ValueError(`unsupported comparison operator: ${op}`);
    for (const [a, b] of [
      [left._domain, right._domain],
      [left._codomain, right._codomain],
    ]) {
      if (!_curves_equal(a!, b!)) {
        if (symbol === '==') return false;
        if (symbol === '!=') return true;
        // Sage's WithEqualityById curve parents do not implement ordering.
        const className = (curve: EllipticCurveGeneric<F>) => {
          if (curve.base_ring === (QQ as unknown))
            return 'EllipticCurve_rational_field_with_category';
          const K = curve.base_ring;
          if (
            K instanceof FiniteFieldExtension ||
            K instanceof PrimeField ||
            K instanceof FiniteFieldPrime ||
            K instanceof GF2Field
          ) {
            return 'EllipticCurve_finite_field_with_category';
          }
          throw new NotImplementedError(
            'SAGE_NOT_IMPLEMENTED: ordering curves over this base ring'
          );
        };
        throw new TypeError(
          `'${symbol}' not supported between instances of '${className(a!)}' and '${className(b!)}'`
        );
      }
    }
    const v = left.tuple(),
      w = right.tuple();
    if (symbol === '==' || symbol === '!=') {
      const equal = v.every((x, i) => x.eq(left._domain.base_ring.__call__(w[i]!)));
      return symbol === '==' ? equal : !equal;
    }
    const compareTuple = (a: [F, F, F, F], b: [F, F, F, F]): number => {
      for (let i = 0; i < 4; i++) {
        const c = _compare_isomorphism_parameters(a[i]!, b[i]!);
        if (c) return c;
      }
      return 0;
    };
    const key = (iso: WeierstrassIsomorphism<F>) => {
      const v = iso.tuple(),
        w = iso.neg().tuple(),
        one = iso._domain.base_ring.one();
      const identity = (t: [F, F, F, F]) => t[0].eq(one) && t.slice(1).every((x) => x.isZero());
      return {
        i: identity(v) || identity(w) ? 0 : 1,
        minimum: compareTuple(v, w) <= 0 ? v : w,
        j: v[0].eq(one) ? 0 : w[0].eq(one) ? 1 : 2,
        v,
      };
    };
    const a = key(left),
      b = key(right);
    const c =
      a.i - b.i || compareTuple(a.minimum, b.minimum) || a.j - b.j || compareTuple(a.v, b.v);
    return symbol === '<' ? c < 0 : symbol === '<=' ? c <= 0 : symbol === '>' ? c > 0 : c >= 0;
  }

  /**
   * Less strict evaluation method for internal use.
   *
   * @param P - A sequence of 3 coordinates defining a point on self
   * @returns The result of evaluating self at the given point
   * @see Reference: sage/schemes/elliptic_curves/weierstrass_morphism.py:WeierstrassIsomorphism._eval
   */
  _eval(P: EllipticCurvePoint<F>): EllipticCurvePoint<F> {
    if (P.is_zero()) {
      return this._codomain.zero();
    }

    const coords = this.call([P.x(), P.y()]) as [F, F];
    return this._codomain.point(coords, false);
  }

  /**
   * Call function for WeierstrassIsomorphism class.
   *
   * @param P - A point on the domain curve
   * @returns The transformed point on the codomain curve
   * @see Reference: sage/schemes/elliptic_curves/weierstrass_morphism.py:WeierstrassIsomorphism._call_
   */
  _call_(P: EllipticCurvePoint<F>): EllipticCurvePoint<F> {
    if (P.is_zero()) {
      return this._codomain.zero();
    }

    const coords = this.call([P.x(), P.y()]) as [F, F];
    return this._codomain.point(coords, false);
  }

  /**
   * Return the inverse of this WeierstrassIsomorphism.
   *
   * @returns The inverse isomorphism
   * @see Reference: sage/schemes/elliptic_curves/weierstrass_morphism.py:WeierstrassIsomorphism.__invert__
   */
  override invert(): WeierstrassIsomorphism<F> {
    const invWI = super.invert();
    return new WeierstrassIsomorphism<F>(this._codomain, invWI.tuple(), this._domain);
  }

  /**
   * Return the composition of a WeierstrassIsomorphism with another elliptic-curve morphism.
   *
   * @param left - Left morphism
   * @param right - Right morphism
   * @returns The composition
   * @see Reference: sage/schemes/elliptic_curves/weierstrass_morphism.py:WeierstrassIsomorphism._composition_impl
   */
  static _composition_impl<F extends FieldElement>(
    left: WeierstrassIsomorphism<F>,
    right: WeierstrassIsomorphism<F>
  ): WeierstrassIsomorphism<F> | null {
    if (left instanceof WeierstrassIsomorphism && right instanceof WeierstrassIsomorphism) {
      if (!_curves_equal(left._domain, right._codomain)) {
        throw new ValueError('Domain of first argument must equal codomain of second');
      }
      const w = left.mul(right);
      return new WeierstrassIsomorphism<F>(right._domain, w.tuple(), left._codomain);
    }

    return null; // NotImplemented
  }

  /**
   * Return the pair of rational maps defining this isomorphism.
   *
   * @returns A pair [X_map, Y_map] of rational functions
   * @see Reference: sage/schemes/elliptic_curves/weierstrass_morphism.py:WeierstrassIsomorphism.rational_maps
   */
  rational_maps(): [string, string] {
    // For a WeierstrassIsomorphism with (u, r, s, t):
    // x' = (x - r) / u^2
    // y' = (y - s*(x - r) - t) / u^3
    const [u, r, s, t] = this.tuple();
    const uSq = u.mul(u);
    const uCubed = uSq.mul(u);

    // Return as strings for now (full rational function support would need polynomial rings)
    const xMap = r.isZero() ? `x/${uSq}` : `(x - ${r})/${uSq}`;

    const yMap = `(y - ${s}*(x - ${r}) - ${t})/${uCubed}`;

    return [xMap, yMap];
  }

  /**
   * Return the x-coordinate rational map of this isomorphism.
   *
   * @returns The x-coordinate rational function
   * @see Reference: sage/schemes/elliptic_curves/weierstrass_morphism.py:WeierstrassIsomorphism.x_rational_map
   */
  x_rational_map(): string {
    const [u, r] = this.tuple();
    const uSq = u.mul(u);
    return r.isZero() ? `x/${uSq}` : `(x - ${r})/${uSq}`;
  }

  /**
   * Return the kernel polynomial of this isomorphism.
   *
   * Isomorphisms have trivial kernel by definition, hence this method always returns 1.
   *
   * @returns The constant polynomial 1
   * @see Reference: sage/schemes/elliptic_curves/weierstrass_morphism.py:WeierstrassIsomorphism.kernel_polynomial
   */
  kernel_polynomial(): string {
    return '1';
  }

  /**
   * Return the dual isogeny of this isomorphism.
   *
   * For isomorphisms, the dual is just the inverse.
   *
   * @returns The inverse isomorphism
   * @see Reference: sage/schemes/elliptic_curves/weierstrass_morphism.py:WeierstrassIsomorphism.dual
   */
  dual(): WeierstrassIsomorphism<F> {
    return this.invert();
  }

  /**
   * Return the negative of this isomorphism, i.e., its composition with [-1].
   *
   * @returns The negated isomorphism
   * @see Reference: sage/schemes/elliptic_curves/weierstrass_morphism.py:WeierstrassIsomorphism.__neg__
   */
  neg(): WeierstrassIsomorphism<F> {
    const [a1, , a3] = this._domain.a_invariants();
    const K = this._domain.base_ring;
    const negOne = K.one().neg() as F;
    const zero = K.zero() as F;
    const negA1 = a1.neg() as F;
    const negA3 = a3.neg() as F;

    const negWI = new baseWI<F>(negOne, zero, negA1, negA3);
    const composed = this.mul(negWI);

    return new WeierstrassIsomorphism<F>(this._domain, composed.tuple(), this._codomain);
  }

  /**
   * Return the Weierstrass scaling factor associated to this isomorphism.
   *
   * The scaling factor is the constant u (in the base field) such that
   * phi^* omega_2 = u omega_1.
   *
   * @returns The scaling factor u
   * @see Reference: sage/schemes/elliptic_curves/weierstrass_morphism.py:WeierstrassIsomorphism.scaling_factor
   */
  scaling_factor(): F {
    return this.u;
  }

  /**
   * Return the inseparable degree of this Weierstrass isomorphism.
   *
   * For isomorphisms, this method always returns one.
   *
   * @returns Always 1n
   * @see Reference: sage/schemes/elliptic_curves/weierstrass_morphism.py:WeierstrassIsomorphism.inseparable_degree
   */
  inseparable_degree(): bigint {
    return 1n;
  }

  /**
   * Check if this Weierstrass isomorphism is the identity.
   *
   * @returns True if (u, r, s, t) == (1, 0, 0, 0)
   * @see Reference: sage/schemes/elliptic_curves/weierstrass_morphism.py:WeierstrassIsomorphism.is_identity
   */
  override is_identity(): boolean {
    return super.is_identity();
  }

  /**
   * Compute the order of this Weierstrass isomorphism if it is an automorphism.
   *
   * A ValueError is raised if the domain is not equal to the codomain.
   * A NotImplementedError is raised if the order is not 1, 2, 3, 4 or 6.
   *
   * @returns The order of the automorphism
   * @see Reference: sage/schemes/elliptic_curves/weierstrass_morphism.py:WeierstrassIsomorphism.order
   */
  order(): bigint {
    // Check if it is an actual endomorphism
    if (!_curves_equal(this._domain, this._codomain)) {
      throw new ValueError('the domain is different from the codomain');
    }

    if (this.is_identity()) {
      return 1n;
    }

    // Compute powers
    const ws2 = WeierstrassIsomorphism._composition_impl(this, this);
    if (ws2?.is_identity()) {
      return 2n;
    }

    const ws3 = ws2 ? WeierstrassIsomorphism._composition_impl(this, ws2) : null;
    if (ws3?.is_identity()) {
      return 3n;
    }

    const ws4 = ws2 ? WeierstrassIsomorphism._composition_impl(ws2, ws2) : null;
    if (ws4?.is_identity()) {
      return 4n;
    }

    const ws6 = ws2 && ws4 ? WeierstrassIsomorphism._composition_impl(ws2, ws4) : null;
    if (ws6?.is_identity()) {
      return 6n;
    }

    throw new NotImplementedError('the order of the endomorphism is not 1, 2, 3, 4 or 6');
  }

  /**
   * Check if this isomorphism equals another.
   *
   * @param other - Another isomorphism
   * @returns True if equal
   */
  eq(other: WeierstrassIsomorphism<F>): boolean {
    return WeierstrassIsomorphism._comparison_impl(this, other, 'eq');
  }

  /**
   * Return the string representation of this WeierstrassIsomorphism.
   *
   * @returns String representation
   * @see Reference: sage/schemes/elliptic_curves/weierstrass_morphism.py:WeierstrassIsomorphism.__repr__
   */
  override toString(): string {
    return `Elliptic-curve morphism:\n  From: ${this._domain}\n  To:   ${this._codomain}\n  Via:  (u,r,s,t) = ${super.toString()}`;
  }
}

/**
 * Helper: Check if two curves have the same a-invariants.
 */
function _curves_equal<F extends FieldElement>(
  E1: EllipticCurveGeneric<F>,
  E2: EllipticCurveGeneric<F>
): boolean {
  if (!_same_base_ring(E1.base_ring, E2.base_ring)) return false;
  const [a1, a2, a3, a4, a6] = E1.a_invariants();
  const [b1, b2, b3, b4, b6] = E2.a_invariants().map((v) => E1.base_ring.__call__(v));
  return a1.eq(b1) && a2.eq(b2) && a3.eq(b3) && a4.eq(b4) && a6.eq(b6);
}

/**
 * Given an elliptic curve E, return the identity morphism on E as a WeierstrassIsomorphism.
 *
 * @param E - An elliptic curve
 * @returns The identity morphism
 * @see Reference: sage/schemes/elliptic_curves/weierstrass_morphism.py:identity_morphism
 */
export function identity_morphism<F extends FieldElement>(
  E: EllipticCurveGeneric<F>
): WeierstrassIsomorphism<F> {
  const K = E.base_ring;
  const one = K.one() as F;
  const zero = K.zero() as F;
  return new WeierstrassIsomorphism<F>(E, [one, zero, zero, zero]);
}

/**
 * Given an elliptic curve E, return the negation endomorphism [-1] of E as a WeierstrassIsomorphism.
 *
 * @param E - An elliptic curve
 * @returns The negation morphism
 * @see Reference: sage/schemes/elliptic_curves/weierstrass_morphism.py:negation_morphism
 */
export function negation_morphism<F extends FieldElement>(
  E: EllipticCurveGeneric<F>
): WeierstrassIsomorphism<F> {
  const K = E.base_ring;
  const negOne = K.one().neg() as F;
  const zero = K.zero() as F;
  const negA1 = E.a1().neg() as F;
  const negA3 = E.a3().neg() as F;
  return new WeierstrassIsomorphism<F>(E, [negOne, zero, negA1, negA3]);
}
