import { EllipticCurveGeneric } from '../../../../packages/sagemath-ts/src/schemes/elliptic_curves/ell_generic.js';
import { EllipticCurve as OptimizedCurve } from '../../../../packages/sagemath-ts/src/schemes/elliptic_curves/ell_finite_field.js';
/**
 * sagemath-ts side of the `ec_advanced` property-test area.
 *
 * Differential oracle for the *advanced* elliptic-curve modules of the port:
 *
 * - `schemes/elliptic_curves/ell_curve_isogeny.ts`  (Velu, Kohel, duals, ...)
 * - `schemes/elliptic_curves/weierstrass_morphism.ts`
 * - `schemes/elliptic_curves/isogeny_class.ts`
 * - `schemes/elliptic_curves/formal_group.ts`
 * - `schemes/elliptic_curves/ell_torsion.ts`
 * - `schemes/elliptic_curves/cm.ts`
 *
 * Cases: tests/property/cases/ec_advanced.cases.json
 * SageMath counterpart: tests/property/python/areas/ec_advanced.py
 *
 * Conventions shared with the SageMath side
 * -----------------------------------------
 * - `p === 0n` selects `QQ` as the base field, otherwise `GF(p)`.
 * - Curves are given by a full list of a-invariants `[a1, a2, a3, a4, a6]` or a
 *   short list `[a4, a6]`.
 * - Every function returns an **already formatted string**; the generic
 *   `formatResult` in `runner.ts` is never relied upon.
 * - Every function is wrapped by `guard()`, which turns a thrown error into the
 *   string `ERR:<ErrorClass>:<message>`. This is deliberate: `compare.ts`
 *   scores "both sides raised" as a pass *without looking at the messages*, so
 *   letting the runner catch the error would hide a disagreement about *why* a
 *   call fails. Returning the message as an ordinary result makes it compared
 *   byte for byte like everything else.
 */

import { GF } from '../../../../packages/sagemath-ts/src/rings/finite_rings/finite_field_constructor.js';
import { QQ } from '../../../../packages/sagemath-ts/src/rings/rational_field.js';
import {
  cm_j_invariants,
  cm_orders,
  discriminants_with_bounded_class_number,
  hilbert_class_polynomial,
  is_cm_j_invariant,
  largest_disc_with_class_number,
  largest_fundamental_disc_with_class_number,
} from '../../../../packages/sagemath-ts/src/schemes/elliptic_curves/cm.js';
import { EllipticCurve } from '../../../../packages/sagemath-ts/src/schemes/elliptic_curves/constructor.js';
import {
  EllipticCurveIsogeny,
  compute_isogeny_bmss,
  compute_isogeny_kernel_polynomial,
  compute_vw_kohel_even_deg1,
  compute_vw_kohel_even_deg3,
  compute_vw_kohel_odd,
  fill_isogeny_matrix,
  isogeny_codomain_from_kernel,
  two_torsion_part,
  unfill_isogeny_matrix,
} from '../../../../packages/sagemath-ts/src/schemes/elliptic_curves/ell_curve_isogeny.js';
import { order_from_multiple } from '../../../../packages/sagemath-ts/src/schemes/elliptic_curves/ell_torsion.js';
import { EllipticCurveTorsionSubgroup } from '../../../../packages/sagemath-ts/src/schemes/elliptic_curves/ell_torsion.js';
import {
  Frobenius_filter,
  IsogenyClassRational,
} from '../../../../packages/sagemath-ts/src/schemes/elliptic_curves/isogeny_class.js';
import {
  WeierstrassIsomorphism,
  _isomorphisms,
  baseWI,
} from '../../../../packages/sagemath-ts/src/schemes/elliptic_curves/weierstrass_morphism.js';

/* eslint-disable @typescript-eslint/no-explicit-any */
// The elliptic-curve stack is generic over a `FieldElement` interface that
// neither `QQ`'s `Rational` nor the finite-field element type satisfies
// structurally (`characteristic` is a method on one and a property on the
// other), exactly as the existing formal-group unit tests note. The runtime
// paths are what this area tests, so the plumbing goes through `any`.
type Any = any;

// ---------------------------------------------------------------------------
// Formatting helpers (mirrored one-for-one in the SageMath area module)
// ---------------------------------------------------------------------------

/** Base field for the encoded characteristic `p` (0 meaning `QQ`). */
function field(p: bigint): Any {
  return p === 0n ? (QQ as Any) : (GF(p) as Any);
}

/** Elliptic curve over `field(p)` from a 2- or 5-element a-invariant list. */
function curve(p: bigint, ainvs: bigint[]): Any {
  const K = field(p);
  return EllipticCurve(K, ainvs.map((a) => K.__call__(a)) as Any) as Any;
}

/** `(a, b, c)` -- Python's own tuple rendering. */
function tup(seq: unknown[]): string {
  return `(${seq.map((x) => fmt(x)).join(', ')})`;
}

/** `[a, b, c]` -- Python's own list rendering. */
function lst(seq: unknown[]): string {
  return `[${seq.map((x) => fmt(x)).join(', ')}]`;
}

/** `str(x)` for the scalar shapes this area produces. */
function fmt(x: unknown): string {
  if (typeof x === 'boolean') return x ? 'True' : 'False';
  if (x === null || x === undefined) return 'None';
  return String(x);
}

function ainvsOf(E: Any): string {
  return tup(E.a_invariants());
}

/** Ascending coefficient list of a univariate polynomial given as an array. */
function polyList(coeffs: unknown[]): string {
  return lst(coeffs);
}

/**
 * All points of `E`, ordered by their string representation.
 *
 * Both runners enumerate the curve independently and in a different internal
 * order, so every table keyed by a point is sorted by `str(P)` -- a pure ASCII
 * comparison that Python and JavaScript agree on.
 */
function sortedPoints(E: Any): Any[] {
  const pts = E.torsion_points() as Any[];
  return pts
    .slice()
    .sort((a, b) => (a.toString() < b.toString() ? -1 : a.toString() > b.toString() ? 1 : 0));
}

/** Coefficients `f[lo] .. f[hi-1]` of a power/Laurent series. */
function series(f: Any, lo: number, hi: number): string {
  const out: unknown[] = [];
  for (let i = lo; i < hi; i++) {
    out.push(f.__getitem__(i).toString());
  }
  return lst(out);
}

/** Coefficients of a bivariate power series, by total degree then by `i`. */
function bivariate(F: Any, prec: number): string {
  const out: string[] = [];
  for (let n = 0; n < prec; n++) {
    for (let i = 0; i <= n; i++) {
      const j = n - i;
      out.push(`(${i},${j})=${F.coefficient(i, j).toString()}`);
    }
  }
  return `[${out.join(', ')}]`;
}

/** Wrap `fn` so thrown errors become comparable `ERR:<class>:<message>` results. */
function guard<A extends unknown[]>(fn: (...args: A) => string): (...args: A) => string {
  return (...args: A) => {
    try {
      return fn(...args);
    } catch (e) {
      const name = e instanceof Error ? e.constructor.name : 'Error';
      const message = e instanceof Error ? e.message : String(e);
      return `ERR:${name}:${message}`;
    }
  };
}

function isogeny(p: bigint, ainvs: bigint[], kx: bigint, ky: bigint): [Any, Any] {
  const E = curve(p, ainvs);
  const K = E.base_ring;
  const ker = E.point([K.__call__(kx), K.__call__(ky)]);
  return [E, new EllipticCurveIsogeny(E, ker)];
}

/** A kernel polynomial as the ascending `bigint[]` the port's Kohel entry points take. */
function kernelPoly(p: bigint, coeffs: bigint[]): bigint[] {
  const m = p === 0n ? null : p;
  return coeffs.map((c) => (m === null ? c : ((c % m) + m) % m));
}

// ===========================================================================
// ell_curve_isogeny.ts -- Velu / Kohel
// ===========================================================================

const raw: Record<string, (...args: Any[]) => string> = {
  /** a-invariants of the codomain of the isogeny with kernel <(kx, ky)>. */
  iso_codomain: (p: bigint, ainvs: bigint[], kx: bigint, ky: bigint) =>
    ainvsOf(isogeny(p, ainvs, kx, ky)[1].codomain()),

  iso_degree: (p: bigint, ainvs: bigint[], kx: bigint, ky: bigint) =>
    isogeny(p, ainvs, kx, ky)[1].degree().toString(),

  iso_kernel_poly: (p: bigint, ainvs: bigint[], kx: bigint, ky: bigint) =>
    polyList(isogeny(p, ainvs, kx, ky)[1].kernel_polynomial()),

  /**
   * `P |-> phi(P)` for **every** point of the domain, sorted by `str(P)`.
   *
   * This is the test that a wrong Velu y-coordinate cannot survive: an image
   * with the wrong sign is still a point of the right x-coordinate, so codomain
   * invariants and degrees stay correct while half of the table is wrong.
   */
  iso_image_table: (p: bigint, ainvs: bigint[], kx: bigint, ky: bigint) => {
    const [E, phi] = isogeny(p, ainvs, kx, ky);
    return sortedPoints(E)
      .map((P) => `${P}|->${phi.call(P)}`)
      .join(' ');
  },

  /** `True` iff every image point actually satisfies the codomain equation. */
  iso_images_on_codomain: (p: bigint, ainvs: bigint[], kx: bigint, ky: bigint) => {
    const [E, phi] = isogeny(p, ainvs, kx, ky);
    const E2 = phi.codomain();
    for (const P of sortedPoints(E)) {
      const Q = phi.call(P);
      if (Q.is_zero()) continue;
      if (!E2.is_on_curve(Q.x(), Q.y())) return `False at ${P}`;
    }
    return 'True';
  },

  /** `True` iff `phi(P + Q) == phi(P) + phi(Q)` for all P, Q on the domain. */
  iso_is_homomorphism: (p: bigint, ainvs: bigint[], kx: bigint, ky: bigint) => {
    const [E, phi] = isogeny(p, ainvs, kx, ky);
    const pts = sortedPoints(E);
    for (const P of pts) {
      for (const Q of pts) {
        if (!phi.call(P.add(Q)).eq(phi.call(P).add(phi.call(Q)))) {
          return `False at ${P} , ${Q}`;
        }
      }
    }
    return 'True';
  },

  /** Sorted list of the points that `phi` sends to infinity. */
  iso_kernel_is_kernel: (p: bigint, ainvs: bigint[], kx: bigint, ky: bigint) => {
    const [E, phi] = isogeny(p, ainvs, kx, ky);
    const ker = sortedPoints(E)
      .filter((P) => phi.call(P).is_zero())
      .map((P) => P.toString())
      .sort();
    return lst(ker);
  },

  iso_scaling_factor: (p: bigint, ainvs: bigint[], kx: bigint, ky: bigint) =>
    isogeny(p, ainvs, kx, ky)[1].scaling_factor().toString(),

  /** (separable, normalized, inseparable_degree, injective, surjective). */
  iso_flags: (p: bigint, ainvs: bigint[], kx: bigint, ky: bigint) => {
    const phi = isogeny(p, ainvs, kx, ky)[1];
    return tup([
      phi.is_separable(),
      phi.is_normalized(),
      phi.inseparable_degree(),
      phi.is_injective(),
      phi.is_surjective(),
    ]);
  },

  /** `str(phi)` -- pins the Weierstrass equation rendering of both curves. */
  iso_repr: (p: bigint, ainvs: bigint[], kx: bigint, ky: bigint) =>
    isogeny(p, ainvs, kx, ky)[1].toString(),

  /** The x-rational map evaluated at every element of the base field. */
  iso_x_rational_map_table: (p: bigint, ainvs: bigint[], kx: bigint, ky: bigint) => {
    const [E, phi] = isogeny(p, ainvs, kx, ky);
    const K = E.base_ring;
    const X = phi.x_rational_map();
    const out: string[] = [];
    for (let v = 0n; v < p; v++) {
      const a = K.__call__(v);
      let val: Any = null;
      try {
        val = X.evaluate(a);
      } catch {
        val = null;
      }
      out.push(`${a}->${val === null ? 'pole' : val}`);
    }
    return out.join(' ');
  },

  /** (degree, codomain a-invariants, kernel polynomial) of the dual isogeny. */
  iso_dual: (p: bigint, ainvs: bigint[], kx: bigint, ky: bigint) => {
    const phiHat = isogeny(p, ainvs, kx, ky)[1].dual();
    return tup([
      phiHat.degree(),
      tup(phiHat.codomain().a_invariants()),
      polyList(phiHat.kernel_polynomial()),
    ]);
  },

  /** `True` iff `phi_hat(phi(P)) == deg(phi) * P` for every P on the domain. */
  iso_dual_is_multiplication: (p: bigint, ainvs: bigint[], kx: bigint, ky: bigint) => {
    const [E, phi] = isogeny(p, ainvs, kx, ky);
    const phiHat = phi.dual();
    const d = phi.degree();
    for (const P of sortedPoints(E)) {
      if (!phiHat.call(phi.call(P)).eq(P.mul(d))) return `False at ${P}`;
    }
    return 'True';
  },

  /** (dual domain == phi codomain, dual codomain == phi domain). */
  iso_dual_domain_codomain: (p: bigint, ainvs: bigint[], kx: bigint, ky: bigint) => {
    const phi = isogeny(p, ainvs, kx, ky)[1];
    const phiHat = phi.dual();
    const same = (A: Any, B: Any) =>
      A.a_invariants().every((a: Any, i: number) => a.eq(B.a_invariants()[i]));
    return tup([same(phiHat.domain(), phi.codomain()), same(phiHat.codomain(), phi.domain())]);
  },

  /** Kohel: build the isogeny from a kernel polynomial, report its codomain. */
  iso_from_kernel_poly: (p: bigint, ainvs: bigint[], coeffs: bigint[]) => {
    const E = curve(p, ainvs);
    const phi = new EllipticCurveIsogeny(E, kernelPoly(p, coeffs) as Any);
    return tup([phi.degree(), tup(phi.codomain().a_invariants())]);
  },

  /** Kohel: the full image table of an isogeny built from a kernel polynomial. */
  iso_from_kernel_poly_images: (p: bigint, ainvs: bigint[], coeffs: bigint[]) => {
    const E = curve(p, ainvs);
    const phi = new EllipticCurveIsogeny(E, kernelPoly(p, coeffs) as Any);
    return sortedPoints(E)
      .map((P) => `${P}|->${phi.call(P)}`)
      .join(' ');
  },

  /** The standalone `isogeny_codomain_from_kernel` on a kernel polynomial. */
  iso_codomain_from_kernel: (p: bigint, ainvs: bigint[], coeffs: bigint[]) =>
    ainvsOf(isogeny_codomain_from_kernel(curve(p, ainvs), kernelPoly(p, coeffs) as Any)),

  iso_two_torsion_part: (p: bigint, ainvs: bigint[], coeffs: bigint[]) =>
    polyList(two_torsion_part(curve(p, ainvs), kernelPoly(p, coeffs))),

  /** `compute_isogeny_kernel_polynomial(E1, E2, ell)` -- domain+codomain only. */
  iso_kernel_poly_from_curves: (p: bigint, ainvs1: bigint[], ainvs2: bigint[], ell: bigint) =>
    polyList(compute_isogeny_kernel_polynomial(curve(p, ainvs1), curve(p, ainvs2), Number(ell))),

  iso_bmss: (p: bigint, ainvs1: bigint[], ainvs2: bigint[], ell: bigint) =>
    polyList(compute_isogeny_bmss(curve(p, ainvs1), curve(p, ainvs2), Number(ell))),

  /** Formal expansion of the isogeny as a power series in `t = -x/y`. */
  iso_formal: (p: bigint, ainvs: bigint[], kx: bigint, ky: bigint, prec: bigint) => {
    const coeffs = isogeny(p, ainvs, kx, ky)[1].formal(Number(prec)) as bigint[];
    const out: unknown[] = [];
    for (let i = 0; i < Number(prec); i++) out.push((coeffs[i] ?? 0n).toString());
    return lst(out);
  },

  /** Velu from an explicit *list* of kernel generators. */
  iso_kernel_list: (p: bigint, ainvs: bigint[], kxs: bigint[], kys: bigint[]) => {
    const E = curve(p, ainvs);
    const K = E.base_ring;
    const kernel = kxs.map((x, i) => E.point([K.__call__(x), K.__call__(kys[i] as bigint)]));
    const phi = new EllipticCurveIsogeny(E, kernel);
    return tup([phi.degree(), tup(phi.codomain().a_invariants())]);
  },

  /** Codomain a-invariants after requesting a particular Weierstrass model. */
  iso_model: (p: bigint, ainvs: bigint[], kx: bigint, ky: bigint, model: bigint) => {
    const models = [null, 'minimal', 'short_weierstrass', 'montgomery'] as const;
    const E = curve(p, ainvs);
    const K = E.base_ring;
    const ker = E.point([K.__call__(kx), K.__call__(ky)]);
    const phi = new EllipticCurveIsogeny(E, ker, null, null, models[Number(model)] as Any, true);
    return ainvsOf(phi.codomain());
  },

  /** `P |-> (-phi)(P)` for every point of the domain. */
  iso_neg_image_table: (p: bigint, ainvs: bigint[], kx: bigint, ky: bigint) => {
    const [E, phi] = isogeny(p, ainvs, kx, ky);
    const neg = phi.neg();
    return sortedPoints(E)
      .map((P) => `${P}|->${neg.call(P)}`)
      .join(' ');
  },

  /** `compute_vw_kohel_odd(b2, b4, b6, s1, s2, s3, n)` on raw field elements. */
  iso_vw_kohel_odd: (
    p: bigint,
    b2: bigint,
    b4: bigint,
    b6: bigint,
    s1: bigint,
    s2: bigint,
    s3: bigint,
    n: bigint
  ) => {
    const K = field(p);
    return tup(
      compute_vw_kohel_odd(
        K.__call__(b2),
        K.__call__(b4),
        K.__call__(b6),
        K.__call__(s1),
        K.__call__(s2),
        K.__call__(s3),
        Number(n)
      )
    );
  },

  iso_vw_kohel_even_deg1: (
    p: bigint,
    x0: bigint,
    y0: bigint,
    a1: bigint,
    a2: bigint,
    a4: bigint
  ) => {
    const K = field(p);
    return tup(
      compute_vw_kohel_even_deg1(
        K.__call__(x0),
        K.__call__(y0),
        K.__call__(a1),
        K.__call__(a2),
        K.__call__(a4)
      )
    );
  },

  iso_vw_kohel_even_deg3: (
    p: bigint,
    b2: bigint,
    b4: bigint,
    s1: bigint,
    s2: bigint,
    s3: bigint
  ) => {
    const K = field(p);
    return tup(
      compute_vw_kohel_even_deg3(
        K.__call__(b2),
        K.__call__(b4),
        K.__call__(s1),
        K.__call__(s2),
        K.__call__(s3)
      )
    );
  },

  /** `fill_isogeny_matrix` on an `n x n` matrix given row-major. */
  iso_fill_matrix: (n: bigint, flat: bigint[]) => {
    const size = Number(n);
    const M: bigint[][] = [];
    for (let i = 0; i < size; i++) {
      M.push(flat.slice(i * size, (i + 1) * size));
    }
    return lst(fill_isogeny_matrix(M).map((row) => tup(row)));
  },

  iso_unfill_matrix: (n: bigint, flat: bigint[]) => {
    const size = Number(n);
    const M: bigint[][] = [];
    for (let i = 0; i < size; i++) {
      M.push(flat.slice(i * size, (i + 1) * size));
    }
    return lst(unfill_isogeny_matrix(M).map((row) => tup(row)));
  },

  // =========================================================================
  // weierstrass_morphism.ts
  // =========================================================================

  /** Every `(u, r, s, t)` with `E1 --(u,r,s,t)--> E2`, in SageMath's own order. */
  wm_isomorphisms: (p: bigint, ainvs1: bigint[], ainvs2: bigint[]) =>
    lst([...(_isomorphisms(curve(p, ainvs1), curve(p, ainvs2)) as Any)].map((u: Any) => tup(u))),

  wm_automorphism_count: (p: bigint, ainvs: bigint[]) => {
    const E = curve(p, ainvs);
    return [...(_isomorphisms(E, E) as Any)].length.toString();
  },

  /** (codomain a-invariants, image of the point (px, py)). */
  wm_apply: (p: bigint, ainvs: bigint[], urst: bigint[], px: bigint, py: bigint) => {
    const E = curve(p, ainvs);
    const K = E.base_ring;
    const w = new WeierstrassIsomorphism(E, urst.map((c) => K.__call__(c)) as Any);
    const P = E.point([K.__call__(px), K.__call__(py)]);
    return tup([tup(w.codomain().a_invariants()), w._call_(P).toString()]);
  },

  /** `P |-> w(P)` for every point of the domain. */
  wm_apply_table: (p: bigint, ainvs: bigint[], urst: bigint[]) => {
    const E = curve(p, ainvs);
    const K = E.base_ring;
    const w = new WeierstrassIsomorphism(E, urst.map((c) => K.__call__(c)) as Any);
    return sortedPoints(E)
      .map((P) => `${P}|->${w._call_(P)}`)
      .join(' ');
  },

  /** `baseWI(urst1) * baseWI(urst2)` as a raw `(u, r, s, t)` tuple. */
  wm_compose: (p: bigint, urst1: bigint[], urst2: bigint[]) => {
    const K = field(p);
    const w1 = new baseWI(...(urst1.map((c) => K.__call__(c)) as [Any, Any, Any, Any]));
    const w2 = new baseWI(...(urst2.map((c) => K.__call__(c)) as [Any, Any, Any, Any]));
    return tup(w1.mul(w2).tuple());
  },

  wm_invert: (p: bigint, urst: bigint[]) => {
    const K = field(p);
    const w = new baseWI(...(urst.map((c) => K.__call__(c)) as [Any, Any, Any, Any]));
    return tup(w.invert().tuple());
  },

  wm_order: (p: bigint, ainvs: bigint[], urst: bigint[]) => {
    const E = curve(p, ainvs);
    const K = E.base_ring;
    return new WeierstrassIsomorphism(E, urst.map((c) => K.__call__(c)) as Any).order().toString();
  },

  wm_scaling_factor: (p: bigint, ainvs: bigint[], urst: bigint[]) => {
    const E = curve(p, ainvs);
    const K = E.base_ring;
    const w = new WeierstrassIsomorphism(E, urst.map((c) => K.__call__(c)) as Any);
    return tup([w.scaling_factor(), w.degree(), w.inseparable_degree()]);
  },

  wm_repr: (p: bigint, ainvs: bigint[], urst: bigint[]) => {
    const E = curve(p, ainvs);
    const K = E.base_ring;
    return new WeierstrassIsomorphism(E, urst.map((c) => K.__call__(c)) as Any).toString();
  },

  /** `(dual tuple, negation tuple)` of a Weierstrass isomorphism. */
  wm_dual_and_neg: (p: bigint, ainvs: bigint[], urst: bigint[]) => {
    const E = curve(p, ainvs);
    const K = E.base_ring;
    const w = new WeierstrassIsomorphism(E, urst.map((c) => K.__call__(c)) as Any);
    return tup([tup(w.dual().tuple()), tup(w.neg().tuple())]);
  },

  // =========================================================================
  // ell_torsion.ts
  // =========================================================================

  tors_invariants: (p: bigint, ainvs: bigint[]) =>
    tup(new EllipticCurveTorsionSubgroup(curve(p, ainvs)).invariants()),

  tors_order: (p: bigint, ainvs: bigint[]) =>
    new EllipticCurveTorsionSubgroup(curve(p, ainvs)).order().toString(),

  tors_point_order: (p: bigint, ainvs: bigint[], px: bigint, py: bigint) => {
    const E = curve(p, ainvs);
    const K = E.base_ring;
    return E.point([K.__call__(px), K.__call__(py)])
      .order()
      .toString();
  },

  tors_order_from_multiple: (p: bigint, ainvs: bigint[], px: bigint, py: bigint, m: bigint) => {
    const E = curve(p, ainvs);
    const K = E.base_ring;
    return order_from_multiple(E.point([K.__call__(px), K.__call__(py)]), m).toString();
  },

  /** `P |-> order(P)` for every point, sorted by `str(P)`. */
  tors_point_orders_table: (p: bigint, ainvs: bigint[]) =>
    sortedPoints(curve(p, ainvs))
      .map((P) => `${P}|->${P.order()}`)
      .join(' '),

  // =========================================================================
  // formal_group.ts
  // =========================================================================

  fg_w: (p: bigint, ainvs: bigint[], prec: bigint) =>
    series(curve(p, ainvs).formal_group().w(Number(prec)), 0, Number(prec)),

  fg_x: (p: bigint, ainvs: bigint[], prec: bigint) =>
    series(curve(p, ainvs).formal_group().x(Number(prec)), -2, Number(prec)),

  fg_y: (p: bigint, ainvs: bigint[], prec: bigint) =>
    series(curve(p, ainvs).formal_group().y(Number(prec)), -3, Number(prec)),

  fg_log: (p: bigint, ainvs: bigint[], prec: bigint) =>
    series(curve(p, ainvs).formal_group().log(Number(prec)), 0, Number(prec)),

  fg_inverse: (p: bigint, ainvs: bigint[], prec: bigint) =>
    series(curve(p, ainvs).formal_group().inverse(Number(prec)), 0, Number(prec)),

  fg_differential: (p: bigint, ainvs: bigint[], prec: bigint) =>
    series(curve(p, ainvs).formal_group().differential(Number(prec)), 0, Number(prec)),

  fg_sigma: (p: bigint, ainvs: bigint[], prec: bigint) =>
    series(curve(p, ainvs).formal_group().sigma(Number(prec)), 0, Number(prec)),

  fg_mult_by_n: (p: bigint, ainvs: bigint[], n: bigint, prec: bigint) =>
    series(curve(p, ainvs).formal_group().mult_by_n(n, Number(prec)), 0, Number(prec)),

  fg_group_law: (p: bigint, ainvs: bigint[], prec: bigint) =>
    bivariate(curve(p, ainvs).formal_group().group_law(Number(prec)), Number(prec)),

  // =========================================================================
  // cm.ts
  // =========================================================================

  cm_hilbert_class_polynomial: (D: bigint) =>
    polyList((hilbert_class_polynomial(D) as Any).coeffs.map((c: Any) => c.toString())),

  cm_orders_list: (h: bigint) => lst(cm_orders(h).map((o) => tup(o))),

  cm_is_cm_j_invariant: (j: bigint) => {
    const [flag, order] = is_cm_j_invariant(j);
    return tup([flag, order === null ? null : tup(order)]);
  },

  cm_largest_fundamental_disc: (h: bigint) => tup(largest_fundamental_disc_with_class_number(h)),

  cm_largest_disc: (h: bigint) => tup(largest_disc_with_class_number(h)),

  cm_discriminants_with_bounded_class_number: (hmax: bigint) => {
    const d = discriminants_with_bounded_class_number(hmax);
    const keys = [...d.keys()].sort((a, b) => (a < b ? -1 : a > b ? 1 : 0));
    return `{${keys.map((h) => `${h}:${lst((d.get(h) as Array<[bigint, bigint]>).map((o) => tup(o)))}`).join(', ')}}`;
  },

  cm_j_invariants_QQ: () => lst(cm_j_invariants(QQ as Any) as unknown[]),

  // =========================================================================
  // isogeny_class.ts
  // =========================================================================

  ic_frobenius_filter: (ainvs: bigint[], primes: bigint[]) =>
    lst(Frobenius_filter(curve(0n, ainvs), primes)),

  ic_rational_class: (ainvs: bigint[]) => {
    const C = new IsogenyClassRational(curve(0n, ainvs));
    const curves = C.curves.map((E) => tup(E.a_invariants()));
    const matrix = C.matrix().map((row) => tup(row));
    const primeMatrix = C.matrix(false).map((row) => tup(row));
    return tup([lst(curves), lst(matrix), lst(primeMatrix)]);
  },
};

export const functions: Record<string, (...args: Any[]) => string> = Object.fromEntries(
  Object.entries(raw).map(([name, fn]) => [name, guard(fn)])
);

import { GFpn } from '../../../../packages/sagemath-ts/src/rings/finite_rings/finite_field_extension.js';
functions.wm_polynomial_root_isomorphisms = (p: bigint, degree: bigint, modulus: bigint[], left: bigint[], right: bigint[], transform: bigint[]) => {
  try {
    const K: Any = degree > 1n ? GFpn(p, Number(degree), modulus.slice(0, -1).map(Number), 'a') : field(p);
    const decode = (v: bigint) => degree > 1n ? K.fromInteger(v) : K.__call__(v);
    const encode = (v: Any) => degree > 1n ? String(v.integer_representation()) : String(v);
    const E = EllipticCurve(K, left.map(decode) as Any);
    const F = EllipticCurve(K, transform.length ? new baseWI(...transform.map(decode) as [Any,Any,Any,Any]).call(E.a_invariants()) as Any : right.map(decode) as Any);
    const rows = [..._isomorphisms(E, F)];
    const value = rows.map(row => [row.map(encode), row.every(v => p === 0n || v.parent === K),
      new baseWI(...row).call(E.a_invariants()).every((v, i) => v.eq(F.a_invariants()[i]!))]);
    let first;
    try { first = {value: new WeierstrassIsomorphism(E, null, F).tuple().map(encode)}; }
    catch (e) { first = {error: (e as Error).name, message: (e as Error).message}; }
    return JSON.stringify({value: [value, first]});
  } catch (e) { return JSON.stringify({error: (e as Error).name, message: (e as Error).message}); }
};
functions.wm_isomorphism_argument_errors = (left: bigint, right: bigint) => {
  try {
    const E = curve(5n, [0n, 1n]), values = [E, null, 7, [], {}];
    return JSON.stringify({value: [..._isomorphisms(values[Number(left)], values[Number(right)])].length});
  } catch (e) { return JSON.stringify({error: (e as Error).name, message: (e as Error).message}); }
};
functions.ec_generic_isomorphism_order = (p: bigint, degree: bigint, modulus: bigint[], left: bigint[], right: bigint[], transform: bigint[]) => {
  try {
    const K: Any = degree > 1n ? GFpn(p, Number(degree), modulus.slice(0, -1).map(Number), 'a') : field(p);
    const decode = (v: bigint) => degree > 1n ? K.fromInteger(v) : K.__call__(v);
    const encode = (v: Any) => degree > 1n ? String(v.integer_representation()) : String(v);
    const E = EllipticCurve(K, left.map(decode) as Any);
    const F = EllipticCurve(K, transform.length ? new baseWI(...transform.map(decode) as [Any,Any,Any,Any]).call(E.a_invariants()) as Any : right.map(decode) as Any);
    let first;
    try {first = {value: E.isomorphism_to(F).map(encode)};}
    catch (e) {first = {error: (e as Error).name, message: (e as Error).message};}
    const ordered = E.isomorphisms(F).map(t => [t.map(encode),t.every(v => p === 0n || v.parent === K)]);
    const same = E.a_invariants().every((v,i) => v.eq(F.a_invariants()[i]!));
    const autos = same ? E.automorphisms().map(t => t.map(encode)) : null;
    return JSON.stringify({value: [first, ordered, autos, E.is_isomorphic(F)]});
  } catch (e) {return JSON.stringify({error: (e as Error).name, message: (e as Error).message});}
};
functions.wm_isomorphism_comparisons = (p: bigint, degree: bigint, modulus: bigint[], left: bigint[], transform: bigint[]) => {
  try {
    const K: Any = degree > 1n ? GFpn(p, Number(degree), modulus.slice(0, -1).map(Number), 'a') : field(p);
    const decode = (v: bigint) => degree > 1n ? K.fromInteger(v) : K.__call__(v);
    const E = EllipticCurve(K, left.map(decode) as Any);
    const F = transform.length ? EllipticCurve(K, new baseWI(...transform.map(decode) as [Any,Any,Any,Any]).call(E.a_invariants()) as Any) : E;
    const morphisms = [..._isomorphisms(E,F)].map(t => new WeierstrassIsomorphism(E,t,F));
    const compare = (a: Any,b: Any) => ['lt','le','eq','ne','gt','ge'].map(op => {
      try {return WeierstrassIsomorphism._comparison_impl(a,b,op);}
      catch(e) {return {error:(e as Error).name,message:(e as Error).message};}
    });
    const matrix=morphisms.map(a => morphisms.map(b => compare(a,b)));
    const identityE=new WeierstrassIsomorphism(E,[K.one(),K.zero(),K.zero(),K.zero()],E);
    const identityF=new WeierstrassIsomorphism(F,[K.one(),K.zero(),K.zero(),K.zero()],F);
    const domains=compare(identityE,identityF),codomains=morphisms.length?compare(identityE,morphisms[0]):null;
    const invalid=[compare(null,identityE),compare(identityE,null),compare(null,7)];
    return JSON.stringify({value:[matrix,domains,codomains,invalid]});
  } catch(e) {return JSON.stringify({error:(e as Error).name,message:(e as Error).message});}
};
functions.ec_is_isomorphic_arguments = (kind: bigint) => {
  try {
    const E=curve(5n,[0n,1n]),other=[null,7,[],{},E][Number(kind)];
    return JSON.stringify({value:E.is_isomorphic(other)});
  } catch(e) {return JSON.stringify({error:(e as Error).name,message:(e as Error).message});}
};
import { PrimeField } from '../../../../packages/sagemath-ts/src/rings/finite_rings/finite_field_extension.js';
import { GF2 } from '../../../../packages/sagemath-ts/src/rings/finite_rings/gf2.js';
functions.ec_isomorphism_parent_guards = (left: bigint,right: bigint) => {
  try {
    const field = (kind: number): Any => {
      if(kind===0)return QQ;
      if(kind>=16){const [p,g]:[bigint,bigint]=({16:[7n,2n],17:[7n,3n],18:[7n,2n],19:[7n,1n],20:[2n,0n],21:[2n,1n]} as Any)[kind];const R=new IsomorphismParentPolynomialRing(new PrimeField(p),'x');return new PrimeField(p,{modulus:R.__call__([-g,1n])});}
      if(kind===12)return new PrimeField(5n);
      if(kind===13)return GF2;
      const primes:Record<number,bigint>={1:2n,2:3n,3:5n,4:7n};
      if(primes[kind])return GF(primes[kind]);
      const [p,d,m,name]:Any=({5:[2n,2,[1,1],'a'],6:[2n,2,[1,1],'b'],7:[3n,2,[1,0],'a'],8:[3n,2,[2,1],'a'],9:[3n,2,[1,0],'b'],10:[5n,2,[2,0],'a'],11:[3n,2,[1,0],'a'],14:[2n,3,[1,1,0],'a'],15:[3n,3,[1,2,0],'a']} as Any)[kind];
      return GFpn(p,d,m,name);
    };
    const curve=(K:Any)=>EllipticCurve(K,K.characteristic===2n?[0n,0n,1n,0n,0n]:K.characteristic===3n?[0n,0n,0n,1n,0n]:[0n,0n,0n,0n,1n]);
    const E=curve(field(Number(left))),F=curve(field(Number(right)));
    const identity=(E:Any)=>new WeierstrassIsomorphism(E,[E.base_ring.one(),E.base_ring.zero(),E.base_ring.zero(),E.base_ring.zero()],E);
    const a=identity(E),b=identity(F);
    const comparison=['lt','le','eq','ne','gt','ge'].map(op=>{try{return WeierstrassIsomorphism._comparison_impl(a,b,op);}catch(e){return {error:(e as Error).name,message:(e as Error).message};}});
    let isomorphic;try{isomorphic={value:E.is_isomorphic(F)};}catch(e){isomorphic={error:(e as Error).name,message:(e as Error).message};}
    return JSON.stringify({value:[isomorphic,comparison]});
  }catch(e){return JSON.stringify({error:(e as Error).name,message:(e as Error).message});}
};
import { Polynomial as IsomorphismTracePolynomial } from '../../../../packages/sagemath-ts/src/rings/polynomial/polynomial_element.js';
functions.ec_isomorphism_root_trace = (p: bigint, degree: bigint, modulus: bigint[], coefficients: bigint[], operation: bigint) => {
  try {
    const K:Any=degree>1n?GFpn(p,Number(degree),modulus.slice(0,-1).map(Number),'a'):field(p);
    const decode=(v:bigint)=>degree>1n?K.fromInteger(v):K.__call__(v),encode=(v:Any)=>degree>1n?String(v.integer_representation()):String(v);
    const E=EllipticCurve(K,coefficients.map(decode) as Any),trace:Any[]=[];
    const prototype:Any=IsomorphismTracePolynomial.prototype,old=prototype.roots;
    try {
      prototype.roots=function(this:Any,options:Any){trace.push([this.coeffs.map(encode),options?.multiplicities??true]);return old.call(this,options);};
      const value=operation===0n?E.is_isomorphic(E):E.isomorphism_to(E).map(encode);
      return JSON.stringify({value:[value,trace]});
    }finally{prototype.roots=old;}
  }catch(e){return JSON.stringify({error:(e as Error).name,message:(e as Error).message});}
};

import { PolynomialRing as IsomorphismParentPolynomialRing } from '../../../../packages/sagemath-ts/src/rings/polynomial/polynomial_ring.js';

import { getrand as coordinateGetrand, setrand as coordinateSetrand } from '@sagemath-ts/parigp-ts';
functions.ec_coordinate_roots = (p: bigint, degree: bigint, modulus: bigint[], coefficients: Any[], coordinate: Any, operation: bigint, seed: bigint, optimized = false) => {
  const trace: Any[] = [];
  let result: Any;
  coordinateSetrand(seed);
  try {
    const K: Any = degree > 1n ? GFpn(p, Number(degree), modulus.slice(0,-1) as Any, 'a') : field(p);
    const decode = (v: Any) => degree > 1n ? K.fromInteger(BigInt(v)) : Array.isArray(v) ? K.__call__(v[0]).div(K.__call__(v[1])) : K.__call__(v);
    const encode = (v: Any) => degree > 1n ? String(v.integer_representation()) : String(v);
    const E: Any = optimized ? OptimizedCurve(K, [coefficients[3], coefficients[4]]) : EllipticCurve(K, coefficients.map(decode) as Any);
    const proto: Any = Object.getPrototypeOf(K.zero()), poly: Any = IsomorphismTracePolynomial.prototype;
    const square = proto.is_square, sqrt = proto.sqrt, roots = poly.roots;
    let depth = 0;
    try {
      proto.is_square = function(this: Any) {
        if (!depth) trace.push(['is_square', encode(this)]);
        depth++; try { return square.call(this); } finally { depth--; }
      };
      proto.sqrt = function(this: Any, options: Any) {
        if (!depth) trace.push(['sqrt', encode(this), options?.all ?? false]);
        depth++; try { return sqrt.call(this, options); } finally { depth--; }
      };
      poly.roots = function(this: Any, options: Any) {
        if (!depth) trace.push(['roots', this.coeffs.map(encode), options?.multiplicities ?? true]);
        depth++; try { return roots.call(this, options); } finally { depth--; }
      };
      let value: Any;
      if (operation === 0n) value = E.is_x_coord(decode(coordinate));
      else if (operation === 4n) value = String(E);
      else if (operation === 3n) value = (E.montgomery_model() as Any).a_invariants().map(encode);
      else {
        let pts: Any = operation === 1n ? E.lift_x(decode(coordinate), true) : [E.lift_x(decode(coordinate))];
        value = pts.map((P: Any) => [encode(optimized ? P.x : P.x()), encode(optimized ? P.y : P.y()), P.curve === E,
          (optimized ? [P.x, P.y, K.one()] : P.xyz()).every((v: Any) => p === 0n || v.parent === K)]);
      }
      result = {value};
    } finally { proto.is_square = square; proto.sqrt = sqrt; poly.roots = roots; }
  } catch(e) { result = {error: (e as Error).name, message: (e as Error).message}; }
  result.calls = trace;
  if (degree > 1n && p !== 2n && operation < 3n) result.state = String(coordinateGetrand());
  return JSON.stringify(result);
};

import { IntegerModRing as CoordinateModRing } from '../../../../packages/sagemath-ts/src/rings/finite_rings/integer_mod_ring.js';
import { Integer as CoordinateInteger } from '../../../../packages/sagemath-ts/src/rings/integer_ring.js';
const coordinateField = (kind: number): Any => {
  if (kind === 0) return QQ;
  if (kind === 12) return new PrimeField(3n);
  if (kind === 13) return GF2;
  if (kind === 14) return GF(257n);
  if (kind === 15) return GF(65537n);
  const primes: Record<number, bigint> = { 1: 2n, 2: 3n, 3: 5n, 4: 7n };
  if (primes[kind]) return GF(primes[kind]);
  const [p, d, T, name]: Any = (
    {
      5: [3n, 2, [1, 0], 'a'],
      6: [3n, 2, [1, 0], 'b'],
      7: [3n, 2, [2, 1], 'a'],
      8: [2n, 2, [1, 1], 'a'],
      9: [2n, 2, [1, 1], 'b'],
      10: [2n, 3, [1, 1, 0], 'a'],
      11: [3n, 3, [1, 2, 0], 'a'],
    } as Any
  )[kind];
  return GFpn(p, d, T, name);
};
function coordinateCoercionResult(
  target: bigint,
  source: bigint,
  numerator: bigint,
  denominator: bigint,
  operation: bigint,
  optimized = false
) {
  try {
    const K = coordinateField(Number(target));
    const E: Any = optimized ? OptimizedCurve(K, [1n, 0n]) : EllipticCurve(
      K,
      K.characteristic === 2n ? [1n, 0n, 0n, 0n, 1n] : [0n, 0n, 0n, 1n, 0n]
    );
    let x: Any;
    if (source === 14n) x = new CoordinateInteger(numerator);
    else if (source === 15n) x = numerator;
    else if (source === 16n) x = Number(numerator) / Number(denominator);
    else if (source === 17n) x = String(numerator);
    else if (source === 18n) x = null;
    else if (source === 19n) x = {};
    else if (source === 20n) x = Boolean(numerator);
    else if (source === 21n) x = NaN;
    else if (source === 22n) x = Infinity;
    else if (source === 23n) x = -Infinity;
    else if (source === 24n) x = 'not an integer';
    else if (source === 25n) x = '1/2';
    else if (source === 26n) x = undefined;
    else if (source === 27n) x = new CoordinateModRing(9n).__call__(numerator);
    else if (source === 28n) x = new CoordinateModRing(6n).__call__(numerator);
    else {
      const L = coordinateField(Number(source));
      x =
        L === QQ
          ? L.__call__(numerator).div(L.__call__(denominator))
          : L.fromInteger
            ? L.fromInteger(((numerator % L.order) + L.order) % L.order)
            : L.__call__(numerator);
    }
    let value: Any;
    if (operation === 4n) value = String(new CoordinateInteger(x));
    else if (operation === 3n) value = String(K.__call__(x));
    else if (operation === 0n) value = E.is_x_coord(x);
    else {
      const pts: Any =
        operation === 7n
          ? E.lift_x(x, true, true)
          : operation === 8n
            ? [E.lift_x(x, false, true)]
            : operation === 1n
              ? E.lift_x(x, true)
              : [E.lift_x(x)];
      value = pts.map((P: Any) => [
        String(typeof P.x === 'function' ? P.x() : P.x),
        String(typeof P.y === 'function' ? P.y() : P.y),
        String(P.curve.base_ring ?? P.curve.field),
        P.curve.a_invariants ? P.curve.a_invariants().map(String) : ['0', '0', '0', String(P.curve.a), String(P.curve.b)],
        P.curve === E,
      ]);
    }
    return JSON.stringify({ value });
  } catch (e) {
    return JSON.stringify({ error: (e as Error).name, message: (e as Error).message });
  }
};

functions.ec_coordinate_coercion = coordinateCoercionResult;
functions.ec_finite_coordinates = (target: bigint, source: bigint, numerator: bigint, denominator: bigint, operation: bigint) =>
  coordinateCoercionResult(target, source, numerator, denominator, operation, true);

functions.ec_curve_base_change = (
  source: bigint,
  target: bigint,
  numerator: bigint,
  denominator: bigint,
  operation: bigint
) => {
  const calls: string[] = [];
  let result: Any;
  try {
    const K = coordinateField(Number(source)),
      L = coordinateField(Number(target));
    const c =
      K === QQ
        ? K.__call__(numerator).div(K.__call__(denominator))
        : K.fromInteger
          ? K.fromInteger(((numerator % K.order) + K.order) % K.order)
          : K.__call__(numerator);
    const E = EllipticCurve(
      K,
      K.characteristic === 2n
        ? [K.one(), K.zero(), K.zero(), c, K.one()]
        : [K.zero(), K.zero(), K.zero(), c, K.one()]
    );
    const convert = L.__call__;
    if (operation >= 2n)
      L.__call__ = function (value: Any) {
        if (calls.length < 5)
          calls.push(
            value?.parent
              ? String(value.parent)
              : typeof value?.denominator === 'bigint'
                ? 'Rational Field'
                : 'Integer Ring'
          );
        return convert.call(this, value);
      };
    try {
      const changed = operation % 2n === 0n ? E.base_extend(L) : E.change_ring(L);
      result = {
        value: [String(changed.base_ring), changed.a_invariants().map(String), changed === E],
      };
    } finally {
      L.__call__ = convert;
    }
  } catch (e) {
    result = { error: (e as Error).name, message: (e as Error).message };
  }
  if (operation >= 2n) result.calls = calls;
  return JSON.stringify(result);
};

functions.ec_lift_extension = (
  kind: bigint,
  coefficients: bigint[],
  numerator: bigint,
  denominator: bigint,
  allPoints: bigint,
  extend: bigint
) => {
  try {
    const K = coordinateField(Number(kind));
    const decode = (n: bigint): Any =>
      K.fromInteger ? K.fromInteger(((n % K.order) + K.order) % K.order) : K.__call__(n);
    const E: Any = EllipticCurve(K, coefficients.map(decode));
    const x = K === QQ ? decode(numerator).div(K.__call__(denominator)) : decode(numerator);
    let points = E.lift_x(x, Boolean(allPoints), Boolean(extend));
    if (!allPoints) points = [points];
    const value = points.map((P: Any) => {
      const M = P.curve.base_ring;
      return [
        String(P.x()),
        String(P.y()),
        String(M),
        P.curve.a_invariants().map(String),
        P.curve === E,
        P.curve.is_on_curve(P.x(), P.y()),
        P.xyz().every(
          (c: Any) => (typeof c.parent === 'function' ? c.parent() : (c.parent ?? QQ)) === M
        ),
        P.mul(2n).xyz().map(String),
        P.neg().xyz().map(String),
      ];
    });
    return JSON.stringify({ value });
  } catch (e) {
    return JSON.stringify({ error: (e as Error).name, message: (e as Error).message });
  }
};

functions.ec_finite_coordinate_roots = (p: bigint, a: bigint, coordinate: bigint, operation: bigint, seed: bigint) =>
  functions.ec_coordinate_roots(p, 1n, [], [0n, 0n, 0n, a, 0n], coordinate, operation, seed, true);


import * as torsionGroups from '../../../../packages/sagemath-ts/src/groups/generic.js';
import * as torsionPoints from '../../../../packages/sagemath-ts/src/schemes/elliptic_curves/ell_point.js';
import { spyOn as torsionSpyOn } from 'bun:test';
functions.ec_primary_torsion = (
  p: bigint,
  T: bigint[],
  coefficients: bigint[],
  prime: bigint,
  bound: bigint
) => {
  const calls: string[][] = [];
  const originalDivision = torsionPoints.division_points,
    originalRelation = torsionGroups.linear_relation;
  const division = torsionSpyOn(torsionPoints, 'division_points').mockImplementation(((
    P: Any,
    n: Any,
    poly: Any
  ) => {
    calls.push(['division_points', String(P), String(n)]);
    return originalDivision(P, n, poly);
  }) as typeof originalDivision);
  const relation = torsionSpyOn(torsionGroups, 'linear_relation').mockImplementation(
    (P: Any, Q: Any, operation: Any) => {
      calls.push(['linear_relation', String(P), String(Q)]);
      return originalRelation(P, Q, operation);
    }
  );
  let result: Any;
  try {
    const K: Any = T.length ? GFpn(p, T.length - 1, T.slice(0, -1) as Any, 'a') : field(p);
    const decode = (n: bigint): Any => (T.length ? K.fromInteger(n) : K.__call__(n));
    const E = EllipticCurve(K, coefficients.map(decode) as Any);
    const basis = E._p_primary_torsion_basis(prime, bound === -99n ? undefined : bound);
    result = { value: basis.map(([P, k]) => [String(P), k]) };
  } catch (e) {
    result = { error: (e as Error).name, message: (e as Error).message };
  } finally {
    division.mockRestore();
    relation.mockRestore();
  }
  return JSON.stringify({ ...result, calls });
};

functions.ec_division_points = (
  p: bigint,
  T: bigint[],
  coefficients: bigint[],
  target: bigint[],
  m: bigint,
  polyOnly: bigint,
  knownOrder: bigint
) => {
  try {
    const K: Any = T.length ? GFpn(p, T.length - 1, T.slice(0, -1) as Any, 'a') : field(p);
    const decode = (n: bigint): Any => (T.length ? K.fromInteger(n) : K.__call__(n));
    const E: Any = EllipticCurve(K, coefficients.map(decode) as Any);
    const P: Any = target.length ? E.point(target.map(decode)) : E.zero();
    if (knownOrder) P._order = P.order();
    else delete P._order;
    coordinateSetrand(1n);
    const result: Any = torsionPoints.division_points(P, m, Boolean(polyOnly) as true);
    const value = Array.isArray(result)
      ? { points: result.map((Q) => [String(Q), Q._order === undefined ? null : String(Q._order)]) }
      : { polynomial: result.coeffs.map(String) };
    return JSON.stringify({
      value,
      target_order: P._order === undefined ? null : String(P._order),
    });
  } catch (e) {
    return JSON.stringify({ error: (e as Error).name, message: (e as Error).message });
  }
};


functions.ec_scalar_order = (modulus: bigint, m: bigint, known: bigint) => {
  const K: Any = new CoordinateModRing(modulus);
  const E: Any = EllipticCurve(K, [1n, 0n]);
  const P: Any = E.point([K.zero(), K.zero()]);
  if (known) P._order = 2n;
  const Q: Any = P.mul(m);
  return JSON.stringify([String(Q), Q._order === undefined ? null : String(Q._order)]);
};


functions.ec_finite_points = (p: bigint, a: bigint, b: bigint) => {
  const E: Any = OptimizedCurve(field(p), [a, b]);
  const calls: string[] = [];
  const original = E.abelian_group;
  if (original) E.abelian_group = () => { calls.push('abelian_group'); return original.call(E); };
  const points = E.points();
  const same = points === E.points();
  let mutation = null;
  try { points[0] = E.zero(); }
  catch (e) { mutation = [(e as Error).name, (e as Error).message]; }
  if (original) E.abelian_group = original;
  const group = original ? E.abelian_group() : null;
  return JSON.stringify({value: points.map(String), calls, same, mutation,
    group_cached: group !== null && group === E.abelian_group(),
    gens_updated: group !== null && E.gens().every((P:Any,i:number)=>P.eq(group.generators[i]))});
};


functions.ec_constructor_model = (p: bigint, T: bigint[], coefficients: bigint[], encoding: bigint, entry: bigint, x: bigint) => {
  try {
    const K: Any = T.length ? GFpn(p, T.length-1, T.slice(0,-1) as Any, 'a') : field(p);
    const encode=(n:bigint):Any => encoding===1n ? (T.length ? K.fromInteger(n) : K.__call__(n))
      : encoding===2n ? new CoordinateInteger(n) : encoding===3n ? QQ.__call__(n).div(QQ.__call__(2n))
      : encoding===4n ? String(n) : n;
    const cs:Any=coefficients.map(encode);
    const E:Any=entry===2n ? new EllipticCurveGeneric(K,cs) : entry===1n ? EllipticCurve(K,cs) : OptimizedCurve(K,cs);
    if(cs.length) cs[cs.length-1]=K.zero();
    const ainvs=E.a_invariants ? E.a_invariants() : [K.zero(),K.zero(),K.zero(),E.a,E.b];
    const value:Any[]=[String(E),ainvs.map(String),String(E.discriminant()),String(E.j_invariant()),
      ainvs.every((c:Any)=>(typeof c.parent==='function'?c.parent():(c.parent??QQ))===K)];
    if(entry!==2n) value.push(E.lift_x(K.__call__(x),true).map((P:Any)=>[String(P),String(P.mul(2n)),String(P.neg())]));
    return JSON.stringify({value});
  } catch(e) {return JSON.stringify({error:(e as Error).name,message:(e as Error).message});}
};
