"""SageMath side of the ``matrix_extended`` property-test area.

Cases: tests/property/cases/matrix_extended.cases.json
TypeScript counterpart: tests/property/typescript/areas/matrix_extended.ts

This is the *oracle*: everything here is plain SageMath, so whatever it prints is
by definition the expected value.  The TypeScript module mirrors these functions
one-for-one (same names, same argument order, same output strings).

Areas covered: ``sage/matrix/matrix2.pyx`` (charpoly, minpoly, kernels, ...),
``sage/matrix/matrix_integer_dense.pyx`` (HNF, SNF, LLL, saturation, Frobenius),
``sage/matrix/matrix_modn_dense_template.pxi`` and
``sage/matrix/matrix_mod2_dense.pyx``.
"""

from sage.all import *


# ---------------------------------------------------------------------------
# Formatting helpers (mirrored exactly in the TypeScript module)
# ---------------------------------------------------------------------------

def _fmt_mat(M):
    """``[[1, 2], [3, 4]]``; a matrix with no rows renders as ``[]``."""
    rows = []
    for i in range(M.nrows()):
        rows.append('[' + ', '.join(str(M[i, j]) for j in range(M.ncols())) + ']')
    return '[' + ', '.join(rows) + ']'


def _fmt_list(xs):
    return '[' + ', '.join(str(x) for x in xs) + ']'


def _fmt_poly(p):
    """Coefficient list, constant term first."""
    return _fmt_list(p.list())


def _fmt_bool(b):
    return 'True' if b else 'False'


def _guard(f):
    """Run ``f``, returning ``ErrorName: message`` instead of raising.

    Only used where the port's exception text is known to be byte-identical, so
    that "both sides refuse this input, for the same reason" is a real assertion
    rather than the harness' vacuous both-errored pass.
    """
    try:
        return f()
    except Exception as e:
        return '%s: %s' % (type(e).__name__, e)


# ---------------------------------------------------------------------------
# Matrix builders
# ---------------------------------------------------------------------------

def _reshape(entries, nrows, ncols):
    entries = [Integer(e) for e in entries]
    if len(entries) != nrows * ncols:
        raise ValueError('expected %d entries, got %d' % (nrows * ncols, len(entries)))
    return [entries[i * ncols:(i + 1) * ncols] for i in range(nrows)]


def _gf(p, nrows, ncols, entries):
    return matrix(GF(Integer(p)), int(nrows), int(ncols), _reshape(entries, int(nrows), int(ncols)))


def _qq(nrows, ncols, entries):
    return matrix(QQ, int(nrows), int(ncols), _reshape(entries, int(nrows), int(ncols)))


def _zz(nrows, ncols, entries):
    return matrix(ZZ, int(nrows), int(ncols), _reshape(entries, int(nrows), int(ncols)))


def _modn(n, nrows, ncols, entries):
    return matrix(Zmod(Integer(n)), int(nrows), int(ncols),
                  _reshape(entries, int(nrows), int(ncols)))


def _mod2(nrows, ncols, entries):
    return matrix(GF(2), int(nrows), int(ncols), _reshape(entries, int(nrows), int(ncols)))


def _degenerate(rows, ncols):
    """Force a pivot-gapped, rank-deficient shape onto random entries.

    See the TypeScript counterpart for the rationale: uniformly random matrices
    over a small field essentially always have pivot columns ``0..r-1``, which is
    precisely the shape under which a transpose-echelonizing kernel routine
    accidentally agrees with SageMath.
    """
    out = [list(r) for r in rows]
    for row in out:
        if ncols > 0:
            row[0] = Integer(0)
        if ncols > 2:
            row[2] = Integer(3) * row[1]
    if out:
        out[-1] = [Integer(0)] * ncols
    return out


def _flatten(rows):
    return [x for row in rows for x in row]


# ---------------------------------------------------------------------------
# matrix2.pyx: charpoly / minpoly
# ---------------------------------------------------------------------------

def gf_charpoly(p, n, entries):
    return _fmt_poly(_gf(p, n, n, entries).charpoly())


def gf_charpoly_hessenberg(p, n, entries):
    """Oracle for the port's ``algorithm='hessenberg'`` characteristic polynomial.

    SageMath's specialised matrix classes do not expose a ``'hessenberg'``
    algorithm (``Matrix_modn_dense_template.charpoly`` accepts only ``'linbox'``
    and ``'generic'``, and ``Matrix_rational_dense.charpoly`` raises
    ``ValueError: no algorithm 'hessenberg'``).  The characteristic polynomial is
    unique, so the oracle for the port's *second* algorithm is simply SageMath's
    characteristic polynomial: ``_charpoly_hessenberg`` must reproduce it exactly.
    """
    return _fmt_poly(_gf(p, n, n, entries).charpoly())


def qq_charpoly(n, entries):
    return _fmt_poly(_qq(n, n, entries).charpoly())


def zz_charpoly(n, entries):
    """Oracle is SageMath's *integer* charpoly (FLINT), a different code path
    from ``matrix(QQ, ...).charpoly()``; the coefficient lists must agree."""
    return _fmt_poly(_zz(n, n, entries).charpoly())


def gf_minpoly(p, n, entries):
    return _fmt_poly(_gf(p, n, n, entries).minpoly())


def qq_minpoly(n, entries):
    return _fmt_poly(_qq(n, n, entries).minpoly())


# ---------------------------------------------------------------------------
# matrix2.pyx: kernels
# ---------------------------------------------------------------------------

def gf_rkm(p, r, c, entries):
    return _fmt_mat(_gf(p, r, c, entries).right_kernel_matrix())


def gf_rkm_pivot(p, r, c, entries):
    return _fmt_mat(_gf(p, r, c, entries).right_kernel_matrix(basis='pivot'))


def gf_lkm(p, r, c, entries):
    return _fmt_mat(_gf(p, r, c, entries).left_kernel_matrix())


def qq_rkm(r, c, entries):
    return _fmt_mat(_qq(r, c, entries).right_kernel_matrix(algorithm='generic'))


def qq_rkm_pivot(r, c, entries):
    return _fmt_mat(_qq(r, c, entries).right_kernel_matrix(algorithm='generic', basis='pivot'))


def gf_rkm_degenerate(p, r, c, entries):
    rows = _degenerate(_reshape(entries, int(r), int(c)), int(c))
    return _fmt_mat(_gf(p, r, c, _flatten(rows)).right_kernel_matrix())


def gf_rref_degenerate(p, r, c, entries):
    rows = _degenerate(_reshape(entries, int(r), int(c)), int(c))
    return _fmt_mat(_gf(p, r, c, _flatten(rows)).rref())


# ---------------------------------------------------------------------------
# matrix2.pyx: misc
# ---------------------------------------------------------------------------

def gf_rank(p, r, c, entries):
    return str(_gf(p, r, c, entries).rank())


def gf_right_nullity(p, r, c, entries):
    return str(_gf(p, r, c, entries).right_nullity())


def gf_det(p, n, entries):
    return str(_gf(p, n, n, entries).determinant())


def qq_det(n, entries):
    return str(_qq(n, n, entries).determinant())


def gf_inverse(p, n, entries):
    return _guard(lambda: _fmt_mat(_gf(p, n, n, entries).inverse()))


def qq_inverse(n, entries):
    return _guard(lambda: _fmt_mat(_qq(n, n, entries).inverse()))


def gf_adjugate(p, n, entries):
    return _fmt_mat(_gf(p, n, n, entries).adjugate())


def gf_permanent(p, r, c, entries):
    """Guarded: SageMath rejects ``m > n`` with
    ``ValueError: must have m <= n, but m (=..) and n (=..)`` (matrix2.pyx:1645).
    Without the guard both sides would merely raise and ``compare.ts`` would
    score the case as a vacuous pass without ever comparing the messages.
    """
    return _guard(lambda: str(_gf(p, r, c, entries).permanent()))


def gf_minors(p, r, c, k, entries):
    return _fmt_list(_gf(p, r, c, entries).minors(int(k)))


def gf_eigenvalues(p, n, entries):
    vals = sorted(Integer(v) for v in _gf(p, n, n, entries).eigenvalues(extend=False))
    return _fmt_list(vals)


def gf_is_diagonalizable(p, n, entries):
    return _fmt_bool(_gf(p, n, n, entries).is_diagonalizable())


def gf_solve_right(p, n, entries_a, entries_b):
    def go():
        A = _gf(p, n, n, entries_a)
        b = _gf(p, n, 1, entries_b)
        return _fmt_mat(A.solve_right(b))
    return _guard(go)


# ---------------------------------------------------------------------------
# matrix2.pyx / matrix_rational_dense: decompositions
# ---------------------------------------------------------------------------

def gf_rref(p, r, c, entries):
    return _fmt_mat(_gf(p, r, c, entries).rref())


def qq_rref(r, c, entries):
    return _fmt_mat(_qq(r, c, entries).rref())


def gf_pivots(p, r, c, entries):
    return _fmt_list(_gf(p, r, c, entries).pivots())


def gf_hessenberg(p, n, entries):
    return _fmt_mat(_gf(p, n, n, entries).hessenberg_form())


def qq_hessenberg(n, entries):
    return _fmt_mat(_qq(n, n, entries).hessenberg_form())


def gf_smith_form(p, r, c, entries):
    return _fmt_mat(_gf(p, r, c, entries).smith_form()[0])


def qq_smith_form(r, c, entries):
    return _fmt_mat(_qq(r, c, entries).smith_form()[0])


def gf_elementary_divisors(p, r, c, entries):
    return _fmt_list(_gf(p, r, c, entries).elementary_divisors())


def gf_hermite_form(p, r, c, entries):
    return _fmt_mat(_gf(p, r, c, entries).hermite_form())


def gf_jordan_form(p, n, entries):
    return _guard(lambda: _fmt_mat(_gf(p, n, n, entries).jordan_form(subdivide=False)))


def qq_jordan_form(n, entries):
    return _guard(lambda: _fmt_mat(_qq(n, n, entries).jordan_form(subdivide=False)))


# ---------------------------------------------------------------------------
# matrix_integer_dense.pyx
# ---------------------------------------------------------------------------

def zz_hnf(r, c, entries):
    return _fmt_mat(_zz(r, c, entries).hermite_form())


def zz_snf(r, c, entries):
    return _fmt_mat(_zz(r, c, entries).smith_form()[0])


def zz_elementary_divisors(r, c, entries):
    return _fmt_list(_zz(r, c, entries).elementary_divisors())


def zz_rank(r, c, entries):
    return str(_zz(r, c, entries).rank())


def zz_det(n, entries):
    return str(_zz(n, n, entries).determinant())


def zz_pivots(r, c, entries):
    return _fmt_list(_zz(r, c, entries).pivots())


def zz_kernel_echelon(r, c, entries):
    return _fmt_mat(_zz(r, c, entries).right_kernel_matrix(basis='echelon'))


def zz_left_kernel_echelon(r, c, entries):
    return _fmt_mat(_zz(r, c, entries).left_kernel_matrix(basis='echelon'))


def zz_lll(r, c, entries):
    return _fmt_mat(_zz(r, c, entries).LLL())


def zz_is_lll_reduced(r, c, entries):
    return _fmt_bool(_zz(r, c, entries).is_LLL_reduced())


def zz_p_minimal(n, p, s_max, entries):
    values = _zz(n, n, entries).p_minimal_polynomials(Integer(p), s_max=Integer(s_max))
    return '[' + ', '.join('[%s, %s]' % (s, _fmt_poly(f)) for s, f in values.items()) + ']'


def zz_null_ideal(n, b, entries):
    generators = _zz(n, n, entries).null_ideal(Integer(b)).gens()
    return '[' + ', '.join(_fmt_poly(f) for f in generators) + ']'


def zz_integer_valued_polynomials(n, entries):
    mu, generators = _zz(n, n, entries).integer_valued_polynomials_generators()
    return _fmt_poly(mu) + '|' + '[' + ', '.join(_fmt_poly(f) for f in generators) + ']'


def zz_height(r, c, entries):
    return str(_zz(r, c, entries).height())


def zz_gcd(r, c, entries):
    """``Matrix_integer_dense.gcd`` (matrix_integer_dense.pyx): gcd of all entries."""
    return str(_zz(r, c, entries).gcd())


def zz_is_primitive(r, c, entries):
    """``Matrix_integer_dense.is_primitive`` (matrix_integer_dense.pyx:1145).

    Note this is *Perron-Frobenius* primitivity -- all entries nonnegative and
    ``A^n`` entrywise positive for some ``n > 0`` -- not "the rows have gcd 1".
    """
    return _fmt_bool(_zz(r, c, entries).is_primitive())


def zz_saturation(r, c, entries):
    return _fmt_mat(_zz(r, c, entries).saturation().hermite_form())


def zz_index_in_saturation(r, c, entries):
    return str(_zz(r, c, entries).index_in_saturation())


def zz_frobenius(n, entries):
    return _fmt_mat(_zz(n, n, entries).frobenius_form())


# ---------------------------------------------------------------------------
# matrix_modn_dense_template.pxi
# ---------------------------------------------------------------------------

def modn_echelon(n, r, c, entries):
    return _guard(lambda: _fmt_mat(_modn(n, r, c, entries).echelon_form()))


def modn_pivots(n, r, c, entries):
    return _guard(lambda: _fmt_list(_modn(n, r, c, entries).pivots()))


def modn_rank(n, r, c, entries):
    return _guard(lambda: str(_modn(n, r, c, entries).rank()))


def modn_det(n, k, entries):
    return _guard(lambda: str(_modn(n, k, k, entries).determinant()))


def modn_charpoly(n, k, entries):
    return _guard(lambda: _fmt_list(_modn(n, k, k, entries).charpoly().list()))


def modn_minpoly(n, k, entries):
    return _guard(lambda: _fmt_list(_modn(n, k, k, entries).minpoly().list()))


def modn_inverse(n, k, entries):
    return _guard(lambda: _fmt_mat(_modn(n, k, k, entries).inverse()))


def modn_rkm(n, r, c, entries):
    return _guard(lambda: _fmt_mat(_modn(n, r, c, entries).right_kernel_matrix()))


def modn_rkm_computed(n, r, c, entries):
    return _guard(lambda: _fmt_mat(_modn(n, r, c, entries).right_kernel_matrix(basis='computed')))


# ---------------------------------------------------------------------------
# matrix_mod2_dense.pyx
# ---------------------------------------------------------------------------

def mod2_echelon(r, c, entries):
    return _fmt_mat(_mod2(r, c, entries).echelon_form())


def mod2_pivots(r, c, entries):
    return _fmt_list(_mod2(r, c, entries).pivots())


def mod2_rank(r, c, entries):
    return str(_mod2(r, c, entries).rank())


def mod2_det(n, entries):
    return str(_mod2(n, n, entries).determinant())


def mod2_inverse(n, entries):
    return _guard(lambda: _fmt_mat(_mod2(n, n, entries).inverse()))


def mod2_rkm(r, c, entries):
    return _fmt_mat(_mod2(r, c, entries).right_kernel_matrix())


def mod2_charpoly(n, entries):
    return _fmt_poly(_mod2(n, n, entries).charpoly())


FUNCTIONS = {
    # matrix_operations: charpoly
    'gf_charpoly': gf_charpoly,
    'gf_charpoly_hessenberg': gf_charpoly_hessenberg,
    'qq_charpoly': qq_charpoly,
    'zz_charpoly': zz_charpoly,
    # matrix_operations: minpoly
    'gf_minpoly': gf_minpoly,
    'qq_minpoly': qq_minpoly,
    # matrix_operations: kernels
    'gf_rkm': gf_rkm,
    'gf_rkm_pivot': gf_rkm_pivot,
    'gf_lkm': gf_lkm,
    'qq_rkm': qq_rkm,
    'qq_rkm_pivot': qq_rkm_pivot,
    'gf_rkm_degenerate': gf_rkm_degenerate,
    'gf_rref_degenerate': gf_rref_degenerate,
    # matrix_operations: misc
    'gf_rank': gf_rank,
    'gf_right_nullity': gf_right_nullity,
    'gf_det': gf_det,
    'qq_det': qq_det,
    'gf_inverse': gf_inverse,
    'qq_inverse': qq_inverse,
    'gf_adjugate': gf_adjugate,
    'gf_permanent': gf_permanent,
    'gf_minors': gf_minors,
    'gf_eigenvalues': gf_eigenvalues,
    'gf_is_diagonalizable': gf_is_diagonalizable,
    'gf_solve_right': gf_solve_right,
    # matrix_decompositions
    'gf_rref': gf_rref,
    'qq_rref': qq_rref,
    'gf_pivots': gf_pivots,
    'gf_hessenberg': gf_hessenberg,
    'qq_hessenberg': qq_hessenberg,
    'gf_smith_form': gf_smith_form,
    'qq_smith_form': qq_smith_form,
    'gf_elementary_divisors': gf_elementary_divisors,
    'gf_hermite_form': gf_hermite_form,
    'gf_jordan_form': gf_jordan_form,
    'qq_jordan_form': qq_jordan_form,
    # matrix_integer
    'zz_hnf': zz_hnf,
    'zz_snf': zz_snf,
    'zz_elementary_divisors': zz_elementary_divisors,
    'zz_rank': zz_rank,
    'zz_det': zz_det,
    'zz_pivots': zz_pivots,
    'zz_kernel_echelon': zz_kernel_echelon,
    'zz_left_kernel_echelon': zz_left_kernel_echelon,
    'zz_lll': zz_lll,
    'zz_is_lll_reduced': zz_is_lll_reduced,
    'zz_p_minimal': zz_p_minimal,
    'zz_null_ideal': zz_null_ideal,
    'zz_integer_valued_polynomials': zz_integer_valued_polynomials,
    'zz_height': zz_height,
    'zz_gcd': zz_gcd,
    'zz_is_primitive': zz_is_primitive,
    'zz_saturation': zz_saturation,
    'zz_index_in_saturation': zz_index_in_saturation,
    'zz_frobenius': zz_frobenius,
    # matrix_modn
    'modn_echelon': modn_echelon,
    'modn_pivots': modn_pivots,
    'modn_rank': modn_rank,
    'modn_det': modn_det,
    'modn_charpoly': modn_charpoly,
    'modn_minpoly': modn_minpoly,
    'modn_inverse': modn_inverse,
    'modn_rkm': modn_rkm,
    'modn_rkm_computed': modn_rkm_computed,
    # matrix_mod2
    'mod2_echelon': mod2_echelon,
    'mod2_pivots': mod2_pivots,
    'mod2_rank': mod2_rank,
    'mod2_det': mod2_det,
    'mod2_inverse': mod2_inverse,
    'mod2_rkm': mod2_rkm,
    'mod2_charpoly': mod2_charpoly,
}
