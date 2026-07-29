"""SageMath side of the ``ec_advanced`` property-test area.

Differential oracle for the *advanced* elliptic-curve modules of the port:

* ``schemes/elliptic_curves/ell_curve_isogeny.ts``  (Velu, Kohel, duals, ...)
* ``schemes/elliptic_curves/weierstrass_morphism.ts``
* ``schemes/elliptic_curves/isogeny_class.ts``
* ``schemes/elliptic_curves/formal_group.ts``
* ``schemes/elliptic_curves/ell_torsion.ts``
* ``schemes/elliptic_curves/cm.ts``

Cases: tests/property/cases/ec_advanced.cases.json
TypeScript counterpart: tests/property/typescript/areas/ec_advanced.ts

Conventions shared with the TypeScript side
-------------------------------------------
* ``p == 0`` selects ``QQ`` as the base field, otherwise ``GF(p)``.
* Curves are given by a full list of a-invariants ``[a1, a2, a3, a4, a6]`` or a
  short list ``[a4, a6]``.
* Every function returns an **already formatted string**; the generic
  ``format_result`` in ``runner.py`` is never relied upon.
* Every function is wrapped by :func:`_guard`, which turns an exception into the
  string ``ERR:<ExceptionClass>:<message>``.  This is deliberate: ``compare.ts``
  scores "both sides raised" as a pass *without looking at the messages*, so
  letting the runner catch the exception would hide a disagreement about *why*
  a call fails.  Returning the message as an ordinary result makes it compared
  byte for byte like everything else.
"""

from sage.all import *

from sage.schemes.elliptic_curves.ell_curve_isogeny import (
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
)
from sage.schemes.elliptic_curves.weierstrass_morphism import (
    WeierstrassIsomorphism,
    _isomorphisms,
    baseWI,
)
from sage.schemes.elliptic_curves.gal_reps_number_field import Frobenius_filter
from sage.schemes.elliptic_curves.cm import (
    cm_j_invariants,
    cm_orders,
    discriminants_with_bounded_class_number,
    hilbert_class_polynomial,
    is_cm_j_invariant,
    largest_disc_with_class_number,
    largest_fundamental_disc_with_class_number,
)
from sage.groups.generic import order_from_multiple


# ---------------------------------------------------------------------------
# Formatting helpers (mirrored one-for-one in the TypeScript area module)
# ---------------------------------------------------------------------------


def _field(p):
    """Base field for the encoded characteristic ``p`` (0 meaning ``QQ``)."""
    p = int(p)
    return QQ if p == 0 else GF(p)


def _curve(p, ainvs):
    """Elliptic curve over ``_field(p)`` from a 2- or 5-element a-invariant list."""
    K = _field(p)
    return EllipticCurve(K, [K(int(a)) for a in ainvs])


def _tup(seq):
    """``(a, b, c)`` -- Python's own tuple rendering."""
    return '(' + ', '.join(str(x) for x in seq) + ')'


def _lst(seq):
    """``[a, b, c]`` -- Python's own list rendering."""
    return '[' + ', '.join(str(x) for x in seq) + ']'


def _ainvs(E):
    return _tup(E.a_invariants())


def _poly_list(f):
    """Ascending coefficient list of a univariate polynomial."""
    return _lst(f.list())


def _sorted_points(E):
    """All points of ``E``, ordered by their string representation.

    Both runners enumerate the curve independently and in a different internal
    order, so every table keyed by a point is sorted by ``str(P)`` -- a pure
    ASCII comparison that Python and JavaScript agree on.
    """
    return sorted(E.points(), key=str)


def _series(f, lo, hi):
    """Coefficients ``f[lo] .. f[hi-1]`` of a power/Laurent series."""
    return _lst([f[i] for i in range(int(lo), int(hi))])


def _bivariate(F, prec):
    """Coefficients of a bivariate power series, by total degree then by ``i``."""
    coeffs = {}
    for mon, c in F.coefficients().items():
        e = mon.exponents()[0]
        coeffs[(int(e[0]), int(e[1]))] = str(c)
    out = []
    for n in range(int(prec)):
        for i in range(n + 1):
            j = n - i
            out.append('(%d,%d)=%s' % (i, j, coeffs.get((i, j), '0')))
    return '[' + ', '.join(out) + ']'


def _guard(fn):
    """Wrap ``fn`` so exceptions become comparable ``ERR:<class>:<message>`` results."""

    def wrapper(*args):
        try:
            return fn(*args)
        except Exception as e:  # noqa: BLE001 -- the message *is* the observation
            return 'ERR:%s:%s' % (type(e).__name__, e)

    return wrapper


def _isogeny(p, ainvs, kx, ky):
    E = _curve(p, ainvs)
    K = E.base_field()
    return E, EllipticCurveIsogeny(E, E(K(int(kx)), K(int(ky))))


def _kernel_polynomial_ring(E):
    return PolynomialRing(E.base_field(), 'x')


# ===========================================================================
# ell_curve_isogeny.ts -- Velu / Kohel
# ===========================================================================


def iso_codomain(p, ainvs, kx, ky):
    """a-invariants of the codomain of the isogeny with kernel <(kx, ky)>."""
    _, phi = _isogeny(p, ainvs, kx, ky)
    return _ainvs(phi.codomain())


def iso_degree(p, ainvs, kx, ky):
    _, phi = _isogeny(p, ainvs, kx, ky)
    return str(phi.degree())


def iso_kernel_poly(p, ainvs, kx, ky):
    _, phi = _isogeny(p, ainvs, kx, ky)
    return _poly_list(phi.kernel_polynomial())


def iso_image_table(p, ainvs, kx, ky):
    """``P |-> phi(P)`` for **every** point of the domain, sorted by ``str(P)``.

    This is the test that a wrong Velu y-coordinate cannot survive: an image
    with the wrong sign is still a point of the right x-coordinate, so codomain
    invariants and degrees stay correct while half of the table is wrong.
    """
    E, phi = _isogeny(p, ainvs, kx, ky)
    return ' '.join('%s|->%s' % (P, phi(P)) for P in _sorted_points(E))


def iso_images_on_codomain(p, ainvs, kx, ky):
    """``True`` iff every image point actually satisfies the codomain equation."""
    E, phi = _isogeny(p, ainvs, kx, ky)
    E2 = phi.codomain()
    for P in E.points():
        Q = phi(P)
        if Q.is_zero():
            continue
        if not E2.is_on_curve(Q[0], Q[1]):
            return 'False at %s' % P
    return 'True'


def iso_is_homomorphism(p, ainvs, kx, ky):
    """``True`` iff ``phi(P + Q) == phi(P) + phi(Q)`` for all P, Q on the domain."""
    E, phi = _isogeny(p, ainvs, kx, ky)
    pts = E.points()
    for P in pts:
        for Q in pts:
            if phi(P + Q) != phi(P) + phi(Q):
                return 'False at %s , %s' % (P, Q)
    return 'True'


def iso_kernel_is_kernel(p, ainvs, kx, ky):
    """Sorted list of the points that ``phi`` sends to infinity."""
    E, phi = _isogeny(p, ainvs, kx, ky)
    ker = sorted((str(P) for P in E.points() if phi(P).is_zero()))
    return _lst(ker)


def iso_scaling_factor(p, ainvs, kx, ky):
    _, phi = _isogeny(p, ainvs, kx, ky)
    return str(phi.scaling_factor())


def iso_flags(p, ainvs, kx, ky):
    """(separable, normalized, inseparable_degree, injective, surjective)."""
    _, phi = _isogeny(p, ainvs, kx, ky)
    return _tup(
        [
            phi.is_separable(),
            phi.is_normalized(),
            phi.inseparable_degree(),
            phi.is_injective(),
            phi.is_surjective(),
        ]
    )


def iso_repr(p, ainvs, kx, ky):
    """``str(phi)`` -- pins the Weierstrass equation rendering of both curves."""
    _, phi = _isogeny(p, ainvs, kx, ky)
    return str(phi)


def iso_x_rational_map_table(p, ainvs, kx, ky):
    """The x-rational map evaluated at every element of the base field."""
    E, phi = _isogeny(p, ainvs, kx, ky)
    K = E.base_field()
    X = phi.x_rational_map()
    out = []
    for v in range(int(K.characteristic())):
        a = K(v)
        try:
            out.append('%s->%s' % (a, X(a)))
        except (ZeroDivisionError, ValueError):
            out.append('%s->pole' % a)
    return ' '.join(out)


def iso_dual(p, ainvs, kx, ky):
    """(degree, codomain a-invariants, kernel polynomial) of the dual isogeny."""
    _, phi = _isogeny(p, ainvs, kx, ky)
    phi_hat = phi.dual()
    return _tup(
        [
            phi_hat.degree(),
            _tup(phi_hat.codomain().a_invariants()),
            _poly_list(phi_hat.kernel_polynomial()),
        ]
    )


def iso_dual_is_multiplication(p, ainvs, kx, ky):
    """``True`` iff ``phi_hat(phi(P)) == deg(phi) * P`` for every P on the domain."""
    E, phi = _isogeny(p, ainvs, kx, ky)
    phi_hat = phi.dual()
    d = phi.degree()
    for P in E.points():
        if phi_hat(phi(P)) != d * P:
            return 'False at %s' % P
    return 'True'


def iso_dual_domain_codomain(p, ainvs, kx, ky):
    """(dual domain == phi codomain, dual codomain == phi domain)."""
    _, phi = _isogeny(p, ainvs, kx, ky)
    phi_hat = phi.dual()
    return _tup(
        [
            phi_hat.domain().a_invariants() == phi.codomain().a_invariants(),
            phi_hat.codomain().a_invariants() == phi.domain().a_invariants(),
        ]
    )


def iso_from_kernel_poly(p, ainvs, coeffs):
    """Kohel: build the isogeny from a kernel polynomial, report its codomain."""
    E = _curve(p, ainvs)
    R = _kernel_polynomial_ring(E)
    f = R([E.base_field()(int(c)) for c in coeffs])
    phi = EllipticCurveIsogeny(E, f)
    return _tup([phi.degree(), _tup(phi.codomain().a_invariants())])


def iso_from_kernel_poly_images(p, ainvs, coeffs):
    """Kohel: the full image table of an isogeny built from a kernel polynomial."""
    E = _curve(p, ainvs)
    R = _kernel_polynomial_ring(E)
    f = R([E.base_field()(int(c)) for c in coeffs])
    phi = EllipticCurveIsogeny(E, f)
    return ' '.join('%s|->%s' % (P, phi(P)) for P in _sorted_points(E))


def iso_codomain_from_kernel(p, ainvs, coeffs):
    """The standalone ``isogeny_codomain_from_kernel`` on a kernel polynomial."""
    E = _curve(p, ainvs)
    R = _kernel_polynomial_ring(E)
    f = R([E.base_field()(int(c)) for c in coeffs])
    return _ainvs(isogeny_codomain_from_kernel(E, f))


def iso_two_torsion_part(p, ainvs, coeffs):
    E = _curve(p, ainvs)
    R = _kernel_polynomial_ring(E)
    f = R([E.base_field()(int(c)) for c in coeffs])
    return _poly_list(two_torsion_part(E, f))


def iso_kernel_poly_from_curves(p, ainvs1, ainvs2, ell):
    """``compute_isogeny_kernel_polynomial(E1, E2, ell)`` -- domain+codomain only."""
    E1 = _curve(p, ainvs1)
    E2 = _curve(p, ainvs2)
    return _poly_list(compute_isogeny_kernel_polynomial(E1, E2, int(ell)))


def iso_bmss(p, ainvs1, ainvs2, ell):
    E1 = _curve(p, ainvs1)
    E2 = _curve(p, ainvs2)
    return _poly_list(compute_isogeny_bmss(E1, E2, int(ell)))


def iso_formal(p, ainvs, kx, ky, prec):
    """Formal expansion of the isogeny as a power series in ``t = -x/y``."""
    _, phi = _isogeny(p, ainvs, kx, ky)
    return _series(phi.formal(int(prec)), 0, int(prec))


def iso_kernel_list(p, ainvs, kxs, kys):
    """Velu from an explicit *list* of kernel generators."""
    E = _curve(p, ainvs)
    K = E.base_field()
    kernel = [E(K(int(x)), K(int(y))) for x, y in zip(kxs, kys)]
    phi = EllipticCurveIsogeny(E, kernel)
    return _tup([phi.degree(), _tup(phi.codomain().a_invariants())])


_MODELS = {0: None, 1: 'minimal', 2: 'short_weierstrass', 3: 'montgomery'}


def iso_model(p, ainvs, kx, ky, model):
    """Codomain a-invariants after requesting a particular Weierstrass model."""
    E = _curve(p, ainvs)
    K = E.base_field()
    phi = EllipticCurveIsogeny(E, E(K(int(kx)), K(int(ky))), model=_MODELS[int(model)])
    return _ainvs(phi.codomain())


def iso_neg_image_table(p, ainvs, kx, ky):
    """``P |-> (-phi)(P)`` for every point of the domain."""
    E, phi = _isogeny(p, ainvs, kx, ky)
    neg = -phi
    return ' '.join('%s|->%s' % (P, neg(P)) for P in _sorted_points(E))


def iso_vw_kohel_odd(p, b2, b4, b6, s1, s2, s3, n):
    """``compute_vw_kohel_odd(b2, b4, b6, s1, s2, s3, n)`` on raw field elements."""
    K = _field(p)
    return _tup(
        compute_vw_kohel_odd(
            K(int(b2)), K(int(b4)), K(int(b6)), K(int(s1)), K(int(s2)), K(int(s3)), int(n)
        )
    )


def iso_vw_kohel_even_deg1(p, x0, y0, a1, a2, a4):
    K = _field(p)
    return _tup(
        compute_vw_kohel_even_deg1(K(int(x0)), K(int(y0)), K(int(a1)), K(int(a2)), K(int(a4)))
    )


def iso_vw_kohel_even_deg3(p, b2, b4, s1, s2, s3):
    K = _field(p)
    return _tup(
        compute_vw_kohel_even_deg3(K(int(b2)), K(int(b4)), K(int(s1)), K(int(s2)), K(int(s3)))
    )


def iso_fill_matrix(n, flat):
    """``fill_isogeny_matrix`` on an ``n x n`` matrix given row-major."""
    n = int(n)
    rows = [[Integer(int(flat[i * n + j])) for j in range(n)] for i in range(n)]
    M = Matrix(ZZ, rows)
    return str(fill_isogeny_matrix(M).rows())


def iso_unfill_matrix(n, flat):
    n = int(n)
    rows = [[Integer(int(flat[i * n + j])) for j in range(n)] for i in range(n)]
    M = Matrix(ZZ, rows)
    return str(unfill_isogeny_matrix(M).rows())


# ===========================================================================
# weierstrass_morphism.ts
# ===========================================================================


def wm_isomorphisms(p, ainvs1, ainvs2):
    """Every ``(u, r, s, t)`` with ``E1 --(u,r,s,t)--> E2``, in Sage's own order."""
    E1 = _curve(p, ainvs1)
    E2 = _curve(p, ainvs2)
    return _lst([_tup(urst) for urst in _isomorphisms(E1, E2)])


def wm_automorphism_count(p, ainvs):
    E = _curve(p, ainvs)
    return str(len(list(_isomorphisms(E, E))))


def wm_apply(p, ainvs, urst, px, py):
    """(codomain a-invariants, image of the point (px, py))."""
    E = _curve(p, ainvs)
    K = E.base_field()
    w = WeierstrassIsomorphism(E, tuple(K(int(c)) for c in urst))
    P = E(K(int(px)), K(int(py)))
    return _tup([_tup(w.codomain().a_invariants()), str(w(P))])


def wm_apply_table(p, ainvs, urst):
    """``P |-> w(P)`` for every point of the domain."""
    E = _curve(p, ainvs)
    K = E.base_field()
    w = WeierstrassIsomorphism(E, tuple(K(int(c)) for c in urst))
    return ' '.join('%s|->%s' % (P, w(P)) for P in _sorted_points(E))


def wm_compose(p, urst1, urst2):
    """``baseWI(urst1) * baseWI(urst2)`` as a raw ``(u, r, s, t)`` tuple."""
    K = _field(p)
    w1 = baseWI(*[K(int(c)) for c in urst1])
    w2 = baseWI(*[K(int(c)) for c in urst2])
    return _tup((w1 * w2).tuple())


def wm_invert(p, urst):
    K = _field(p)
    return _tup((~baseWI(*[K(int(c)) for c in urst])).tuple())


def wm_order(p, ainvs, urst):
    E = _curve(p, ainvs)
    K = E.base_field()
    return str(WeierstrassIsomorphism(E, tuple(K(int(c)) for c in urst)).order())


def wm_scaling_factor(p, ainvs, urst):
    E = _curve(p, ainvs)
    K = E.base_field()
    w = WeierstrassIsomorphism(E, tuple(K(int(c)) for c in urst))
    return _tup([w.scaling_factor(), w.degree(), w.inseparable_degree()])


def wm_repr(p, ainvs, urst):
    E = _curve(p, ainvs)
    K = E.base_field()
    return str(WeierstrassIsomorphism(E, tuple(K(int(c)) for c in urst)))


def wm_dual_and_neg(p, ainvs, urst):
    """``(dual tuple, negation tuple)`` of a Weierstrass isomorphism."""
    E = _curve(p, ainvs)
    K = E.base_field()
    w = WeierstrassIsomorphism(E, tuple(K(int(c)) for c in urst))
    return _tup([_tup(w.dual().tuple()), _tup((-w).tuple())])


# ===========================================================================
# ell_torsion.ts  (the port implements this for finite fields; SageMath's
# oracle for the group structure of E(F_q) is E.abelian_group())
# ===========================================================================


def tors_invariants(p, ainvs):
    E = _curve(p, ainvs)
    return _tup(E.abelian_group().invariants())


def tors_order(p, ainvs):
    E = _curve(p, ainvs)
    return str(E.cardinality())


def tors_point_order(p, ainvs, px, py):
    E = _curve(p, ainvs)
    K = E.base_field()
    return str(E(K(int(px)), K(int(py))).order())


def tors_order_from_multiple(p, ainvs, px, py, m):
    E = _curve(p, ainvs)
    K = E.base_field()
    return str(order_from_multiple(E(K(int(px)), K(int(py))), Integer(int(m))))


def tors_point_orders_table(p, ainvs):
    """``P |-> order(P)`` for every point, sorted by ``str(P)``."""
    E = _curve(p, ainvs)
    return ' '.join('%s|->%s' % (P, P.order()) for P in _sorted_points(E))


# ===========================================================================
# formal_group.ts
# ===========================================================================


def _formal(p, ainvs):
    return _curve(p, ainvs).formal_group()


def fg_w(p, ainvs, prec):
    return _series(_formal(p, ainvs).w(int(prec)), 0, int(prec))


def fg_x(p, ainvs, prec):
    return _series(_formal(p, ainvs).x(int(prec)), -2, int(prec))


def fg_y(p, ainvs, prec):
    return _series(_formal(p, ainvs).y(int(prec)), -3, int(prec))


def fg_log(p, ainvs, prec):
    return _series(_formal(p, ainvs).log(int(prec)), 0, int(prec))


def fg_inverse(p, ainvs, prec):
    return _series(_formal(p, ainvs).inverse(int(prec)), 0, int(prec))


def fg_differential(p, ainvs, prec):
    return _series(_formal(p, ainvs).differential(int(prec)), 0, int(prec))


def fg_sigma(p, ainvs, prec):
    return _series(_formal(p, ainvs).sigma(int(prec)), 0, int(prec))


def fg_mult_by_n(p, ainvs, n, prec):
    return _series(_formal(p, ainvs).mult_by_n(Integer(int(n)), int(prec)), 0, int(prec))


def fg_group_law(p, ainvs, prec):
    return _bivariate(_formal(p, ainvs).group_law(int(prec)), int(prec))


# ===========================================================================
# cm.ts
# ===========================================================================


def cm_hilbert_class_polynomial(D):
    return _poly_list(hilbert_class_polynomial(Integer(int(D))))


def cm_orders_list(h):
    return _lst([_tup(o) for o in cm_orders(Integer(int(h)))])


def cm_is_cm_j_invariant(j):
    return str(is_cm_j_invariant(QQ(int(j))))


def cm_largest_fundamental_disc(h):
    return _tup(largest_fundamental_disc_with_class_number(Integer(int(h))))


def cm_largest_disc(h):
    return _tup(largest_disc_with_class_number(Integer(int(h))))


def cm_discriminants_with_bounded_class_number(hmax):
    d = discriminants_with_bounded_class_number(Integer(int(hmax)))
    parts = []
    for h in sorted(d):
        parts.append('%s:%s' % (h, _lst([_tup(o) for o in d[h]])))
    return '{' + ', '.join(parts) + '}'


def cm_j_invariants_QQ():
    return _lst(cm_j_invariants(QQ))


# ===========================================================================
# isogeny_class.ts
# ===========================================================================


def ic_frobenius_filter(ainvs, primes):
    E = EllipticCurve(QQ, [Integer(int(a)) for a in ainvs])
    return _lst(Frobenius_filter(E, [Integer(int(l)) for l in primes]))


_RAW = {
    # --- ell_curve_isogeny -------------------------------------------------
    'iso_codomain': iso_codomain,
    'iso_degree': iso_degree,
    'iso_kernel_poly': iso_kernel_poly,
    'iso_image_table': iso_image_table,
    'iso_images_on_codomain': iso_images_on_codomain,
    'iso_is_homomorphism': iso_is_homomorphism,
    'iso_kernel_is_kernel': iso_kernel_is_kernel,
    'iso_scaling_factor': iso_scaling_factor,
    'iso_flags': iso_flags,
    'iso_repr': iso_repr,
    'iso_x_rational_map_table': iso_x_rational_map_table,
    'iso_dual': iso_dual,
    'iso_dual_is_multiplication': iso_dual_is_multiplication,
    'iso_dual_domain_codomain': iso_dual_domain_codomain,
    'iso_from_kernel_poly': iso_from_kernel_poly,
    'iso_from_kernel_poly_images': iso_from_kernel_poly_images,
    'iso_codomain_from_kernel': iso_codomain_from_kernel,
    'iso_two_torsion_part': iso_two_torsion_part,
    'iso_kernel_poly_from_curves': iso_kernel_poly_from_curves,
    'iso_bmss': iso_bmss,
    'iso_formal': iso_formal,
    'iso_kernel_list': iso_kernel_list,
    'iso_model': iso_model,
    'iso_neg_image_table': iso_neg_image_table,
    'iso_vw_kohel_odd': iso_vw_kohel_odd,
    'iso_vw_kohel_even_deg1': iso_vw_kohel_even_deg1,
    'iso_vw_kohel_even_deg3': iso_vw_kohel_even_deg3,
    'iso_fill_matrix': iso_fill_matrix,
    'iso_unfill_matrix': iso_unfill_matrix,
    # --- weierstrass_morphism ---------------------------------------------
    'wm_isomorphisms': wm_isomorphisms,
    'wm_automorphism_count': wm_automorphism_count,
    'wm_apply': wm_apply,
    'wm_apply_table': wm_apply_table,
    'wm_compose': wm_compose,
    'wm_invert': wm_invert,
    'wm_order': wm_order,
    'wm_scaling_factor': wm_scaling_factor,
    'wm_repr': wm_repr,
    'wm_dual_and_neg': wm_dual_and_neg,
    # --- ell_torsion -------------------------------------------------------
    'tors_invariants': tors_invariants,
    'tors_order': tors_order,
    'tors_point_order': tors_point_order,
    'tors_order_from_multiple': tors_order_from_multiple,
    'tors_point_orders_table': tors_point_orders_table,
    # --- formal_group ------------------------------------------------------
    'fg_w': fg_w,
    'fg_x': fg_x,
    'fg_y': fg_y,
    'fg_log': fg_log,
    'fg_inverse': fg_inverse,
    'fg_differential': fg_differential,
    'fg_sigma': fg_sigma,
    'fg_mult_by_n': fg_mult_by_n,
    'fg_group_law': fg_group_law,
    # --- cm ----------------------------------------------------------------
    'cm_hilbert_class_polynomial': cm_hilbert_class_polynomial,
    'cm_orders_list': cm_orders_list,
    'cm_is_cm_j_invariant': cm_is_cm_j_invariant,
    'cm_largest_fundamental_disc': cm_largest_fundamental_disc,
    'cm_largest_disc': cm_largest_disc,
    'cm_discriminants_with_bounded_class_number': cm_discriminants_with_bounded_class_number,
    'cm_j_invariants_QQ': cm_j_invariants_QQ,
    # --- isogeny_class -----------------------------------------------------
    'ic_frobenius_filter': ic_frobenius_filter,
}

FUNCTIONS = {name: _guard(fn) for name, fn in _RAW.items()}
