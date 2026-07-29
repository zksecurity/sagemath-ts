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
import { Frobenius_filter } from '../../../../packages/sagemath-ts/src/schemes/elliptic_curves/isogeny_class.js';
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
};

export const functions: Record<string, (...args: Any[]) => string> = Object.fromEntries(
  Object.entries(raw).map(([name, fn]) => [name, guard(fn)])
);
