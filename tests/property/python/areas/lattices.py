"""SageMath side of the ``lattices`` property-test area.

Covers ``sage/modules/free_module.py``, ``free_module_element.pyx``,
``free_module_integer.py`` and the BKZ entry points that
``sage/modules/free_module_integer.py`` exposes.

Cases: tests/property/cases/lattices.cases.json
TypeScript counterpart: tests/property/typescript/areas/lattices.ts

Every function returns an **already formatted string** so that the two
languages land on byte-identical transcripts without touching the shared
runners.  The formatting helpers below are duplicated verbatim (modulo
syntax) in the TypeScript area module.

Matrices are passed in flattened, row-major, together with the column count,
because the shared argument generators can only produce integers and flat
lists of integers.
"""

from sage.all import (
    FreeModule,
    GF,
    Integer,
    QQ,
    VectorSpace,
    ZZ,
    matrix,
    vector,
)
from sage.modules.free_module_integer import IntegerLattice

# ---------------------------------------------------------------------------
# Formatting helpers (mirrored in the TypeScript area module)
# ---------------------------------------------------------------------------


def _fmt_bool(b):
    return 'True' if b else 'False'


def _fmt_list(xs):
    return '[' + ', '.join(str(x) for x in xs) + ']'


def _fmt_mat(rows):
    return '[' + ', '.join(_fmt_list(row) for row in rows) + ']'


def _mat_rows(M):
    return [[M[i, j] for j in range(M.ncols())] for i in range(M.nrows())]


def _fmt_matrix(M):
    return _fmt_mat(_mat_rows(M))


def _fmt_module(W):
    """rank + echelon basis matrix; the canonical fingerprint of a submodule."""
    return 'rank=%d basis=%s' % (W.rank(), _fmt_matrix(W.basis_matrix()))


# ---------------------------------------------------------------------------
# Input reshaping
# ---------------------------------------------------------------------------


def _reshape(flat, ncols):
    """Turn a flat row-major list into a list of rows of length ``ncols``."""
    ncols = int(ncols)
    flat = list(flat)
    if ncols <= 0:
        raise ValueError('ncols must be positive')
    if len(flat) % ncols != 0:
        raise ValueError('flat length %d is not a multiple of ncols %d' % (len(flat), ncols))
    return [flat[i * ncols:(i + 1) * ncols] for i in range(len(flat) // ncols)]


def _reshape_qq(flat, ncols):
    """``flat`` holds (numerator, denominator) pairs, row-major."""
    ncols = int(ncols)
    flat = list(flat)
    if len(flat) % 2 != 0:
        raise ValueError('rational payload must have even length')
    entries = [QQ(flat[2 * i]) / QQ(flat[2 * i + 1]) for i in range(len(flat) // 2)]
    return _reshape(entries, ncols)


# ---------------------------------------------------------------------------
# sage/modules/free_module.py -- spans, echelon bases, rank
# ---------------------------------------------------------------------------


def zz_span(flat, ncols):
    """Rank and Hermite echelon basis of the ZZ-span of the given rows."""
    V = FreeModule(ZZ, int(ncols))
    W = V.span(_reshape(flat, ncols))
    return _fmt_module(W)


def zz_span_qq_entries(flat, ncols):
    """ZZ-span of *non-integral* generators (rationals given as num/den pairs)."""
    V = FreeModule(ZZ, int(ncols))
    W = V.span(_reshape_qq(flat, ncols))
    return _fmt_module(W)


def zz_span_of_basis(flat, ncols):
    """User basis vs. echelonized basis of ``span_of_basis``."""
    V = FreeModule(ZZ, int(ncols))
    W = V.span_of_basis(_reshape(flat, ncols))
    return 'basis=%s echelon=%s' % (
        _fmt_matrix(W.basis_matrix()),
        _fmt_matrix(W.echelonized_basis_matrix()),
    )


def qq_span(flat, ncols):
    """Dimension and RREF basis of a QQ subspace."""
    V = VectorSpace(QQ, int(ncols))
    W = V.subspace(_reshape(flat, ncols))
    return 'dim=%d basis=%s' % (W.dimension(), _fmt_matrix(W.basis_matrix()))


def qq_span_qq_entries(flat, ncols):
    V = VectorSpace(QQ, int(ncols))
    W = V.subspace(_reshape_qq(flat, ncols))
    return 'dim=%d basis=%s' % (W.dimension(), _fmt_matrix(W.basis_matrix()))


def gf_span(flat, ncols, p):
    """Dimension and echelon basis of a GF(p) subspace."""
    V = VectorSpace(GF(Integer(p)), int(ncols))
    W = V.subspace(_reshape(flat, ncols))
    return 'dim=%d basis=%s' % (W.dimension(), _fmt_matrix(W.basis_matrix()))


def zz_intersection(a, b, ncols):
    V = FreeModule(ZZ, int(ncols))
    A = V.span(_reshape(a, ncols))
    B = V.span(_reshape(b, ncols))
    return '%s | %s' % (_fmt_module(A.intersection(B)), _fmt_module(B.intersection(A)))


def qq_intersection(a, b, ncols):
    V = VectorSpace(QQ, int(ncols))
    A = V.subspace(_reshape(a, ncols))
    B = V.subspace(_reshape(b, ncols))
    return '%s | %s' % (_fmt_module(A.intersection(B)), _fmt_module(B.intersection(A)))


def gf_intersection(a, b, ncols, p):
    V = VectorSpace(GF(Integer(p)), int(ncols))
    A = V.subspace(_reshape(a, ncols))
    B = V.subspace(_reshape(b, ncols))
    return '%s | %s' % (_fmt_module(A.intersection(B)), _fmt_module(B.intersection(A)))


def zz_sum(a, b, ncols):
    V = FreeModule(ZZ, int(ncols))
    A = V.span(_reshape(a, ncols))
    B = V.span(_reshape(b, ncols))
    return _fmt_module(A + B)


def zz_discriminant(flat, ncols):
    V = FreeModule(ZZ, int(ncols))
    return str(V.span(_reshape(flat, ncols)).discriminant())


def zz_saturation(flat, ncols):
    V = FreeModule(ZZ, int(ncols))
    return _fmt_module(V.span(_reshape(flat, ncols)).saturation())


def zz_index_in(a, b, ncols):
    """Index of span(a) in span(b); ``+Infinity`` when the index is infinite."""
    V = FreeModule(ZZ, int(ncols))
    A = V.span(_reshape(a, ncols))
    B = V.span(_reshape(b, ncols))
    return str(A.index_in(B))


def zz_index_in_ambient(flat, ncols):
    V = FreeModule(ZZ, int(ncols))
    return str(V.span(_reshape(flat, ncols)).index_in(V))


def zz_coordinates(flat, ncols, target):
    V = FreeModule(ZZ, int(ncols))
    W = V.span(_reshape(flat, ncols))
    return _fmt_list(W.coordinates(V(list(target))))


def zz_coordinates_user_basis(flat, ncols, target):
    V = FreeModule(ZZ, int(ncols))
    W = V.span_of_basis(_reshape(flat, ncols))
    return _fmt_list(W.coordinates(V(list(target))))


def zz_is_submodule(a, b, ncols):
    V = FreeModule(ZZ, int(ncols))
    A = V.span(_reshape(a, ncols))
    B = V.span(_reshape(b, ncols))
    return '%s %s' % (_fmt_bool(A.is_submodule(B)), _fmt_bool(B.is_submodule(A)))


def zz_module_eq(a, b, ncols):
    V = FreeModule(ZZ, int(ncols))
    A = V.span(_reshape(a, ncols))
    B = V.span(_reshape(b, ncols))
    return _fmt_bool(A == B)


def zz_quotient_invariants(flat, ncols):
    V = FreeModule(ZZ, int(ncols))
    W = V.span(_reshape(flat, ncols))
    return _fmt_list(list((V / W).invariants()))


def qq_complement(flat, ncols):
    V = VectorSpace(QQ, int(ncols))
    W = V.subspace(_reshape(flat, ncols))
    return _fmt_matrix(W.complement().basis_matrix())


def gf_complement(flat, ncols, p):
    V = VectorSpace(GF(Integer(p)), int(ncols))
    W = V.subspace(_reshape(flat, ncols))
    return _fmt_matrix(W.complement().basis_matrix())


def gf_cardinality(n, p):
    return str(VectorSpace(GF(Integer(p)), int(n)).cardinality())


# ---------------------------------------------------------------------------
# sage/modules/free_module_element.pyx
# ---------------------------------------------------------------------------


def vector_ops_zz(a, b):
    u = vector(ZZ, list(a))
    v = vector(ZZ, list(b))
    return 'add=%s sub=%s neg=%s dot=%s pairwise=%s hw=%d support=%s norm2=%s' % (
        _fmt_list(u + v),
        _fmt_list(u - v),
        _fmt_list(-u),
        str(u.dot_product(v)),
        _fmt_list(u.pairwise_product(v)),
        u.hamming_weight(),
        _fmt_list(u.support()),
        str(u.dot_product(u)),
    )


def vector_cross_zz(a, b):
    u = vector(ZZ, list(a))
    v = vector(ZZ, list(b))
    return _fmt_list(u.cross_product(v))


def vector_ops_gf(a, b, p):
    K = GF(Integer(p))
    u = vector(K, list(a))
    v = vector(K, list(b))
    return 'add=%s dot=%s pairwise=%s is_zero=%s hw=%d' % (
        _fmt_list(u + v),
        str(u.dot_product(v)),
        _fmt_list(u.pairwise_product(v)),
        _fmt_bool(u.is_zero()),
        u.hamming_weight(),
    )


def vector_scalar_zz(a, c):
    u = vector(ZZ, list(a))
    return _fmt_list(Integer(c) * u)


# ---------------------------------------------------------------------------
# Exact LLL-reducedness predicate
#
# Deliberately independent of the implementation under test: it re-derives the
# Gram-Schmidt data over QQ and checks Definition (delta, eta) exactly.  delta
# and eta are the SageMath defaults 0.99 / 0.501, as exact rationals.
# ---------------------------------------------------------------------------

_DELTA = QQ(99) / QQ(100)
_ETA = QQ(501) / QQ(1000)


def _is_lll_reduced_exact(rows):
    n = len(rows)
    if n <= 1:
        return True
    B = matrix(QQ, rows)
    G, mu = B.gram_schmidt()
    norms = [G[i].dot_product(G[i]) for i in range(n)]
    for i in range(n):
        for j in range(i):
            if abs(mu[i, j]) > _ETA:
                return False
    for i in range(1, n):
        lhs = _DELTA * norms[i - 1]
        rhs = norms[i] + mu[i, i - 1] ** 2 * norms[i - 1]
        if lhs > rhs:
            return False
    return True


def _hnf_rows(rows, ncols):
    M = matrix(ZZ, len(rows), int(ncols), [list(r) for r in rows]) if rows else matrix(ZZ, 0, int(ncols))
    return _mat_rows(M.hermite_form(include_zero_rows=False))


def _lattice_fingerprint(rows, ncols):
    """rank + HNF of the row span + exact LLL-reducedness + first squared norm."""
    n = len(rows)
    hnf = _hnf_rows(rows, ncols)
    first = sum(x * x for x in rows[0]) if n else 0
    return 'rows=%d hnf=%s reduced=%s first_norm2=%s' % (
        n,
        _fmt_mat(hnf),
        _fmt_bool(_is_lll_reduced_exact(rows)),
        str(first),
    )


# ---------------------------------------------------------------------------
# sage/modules/free_module_integer.py + BKZ
# ---------------------------------------------------------------------------


def lattice_rank_degree(flat, ncols):
    L = IntegerLattice(_reshape(flat, ncols))
    return 'rank=%d degree=%d' % (L.rank(), L.degree())


def lattice_basis(flat, ncols):
    """Basis matrix after the constructor's LLL reduction (exact comparison)."""
    L = IntegerLattice(_reshape(flat, ncols))
    return _fmt_matrix(L.basis_matrix())


def lattice_lll_fingerprint(flat, ncols):
    """Canonical fingerprint of the constructor's LLL-reduced basis.

    Unlike :func:`lattice_basis` this does not pin the reduced basis itself:
    an LLL-reduced basis is not unique (``v`` and ``-v``, or ``v`` and
    ``v - b_0`` when ``mu`` is exactly 1/2, are equally valid), so only the
    lattice it spans (its HNF), the reducedness predicate and the first
    squared norm are compared.
    """
    L = IntegerLattice(_reshape(flat, ncols))
    return _lattice_fingerprint(_mat_rows(L.basis_matrix()), ncols)


def lll_exact(flat, ncols):
    """LLL of a *non*-reduced lattice, exact matrix."""
    L = IntegerLattice(_reshape(flat, ncols), lll_reduce=False)
    return _fmt_matrix(L.LLL())


def lll_invariants(flat, ncols):
    L = IntegerLattice(_reshape(flat, ncols), lll_reduce=False)
    R = L.LLL()
    return _lattice_fingerprint(_mat_rows(R), ncols)


def bkz_exact(flat, ncols, block_size):
    L = IntegerLattice(_reshape(flat, ncols), lll_reduce=False)
    return _fmt_matrix(L.BKZ(block_size=int(block_size)))


def bkz_invariants(flat, ncols, block_size):
    L = IntegerLattice(_reshape(flat, ncols), lll_reduce=False)
    R = L.BKZ(block_size=int(block_size))
    return _lattice_fingerprint(_mat_rows(R), ncols)


def hkz_first_norm2(flat, ncols):
    """HKZ puts a shortest vector first, so ||b_1||^2 == lambda_1^2 is canonical."""
    L = IntegerLattice(_reshape(flat, ncols), lll_reduce=False)
    R = L.HKZ()
    row = [R[0, j] for j in range(R.ncols())]
    return str(sum(x * x for x in row))


def shortest_vector_norm2(flat, ncols):
    L = IntegerLattice(_reshape(flat, ncols))
    v = L.shortest_vector()
    return str(v.dot_product(v))


def closest_vector_dist2(flat, ncols, target):
    L = IntegerLattice(_reshape(flat, ncols))
    t = vector(ZZ, list(target))
    v = L.closest_vector(t)
    d = v - t
    return str(d.dot_product(d))


def lattice_volume(flat, ncols):
    L = IntegerLattice(_reshape(flat, ncols))
    return str(L.volume())


def lattice_discriminant(flat, ncols):
    L = IntegerLattice(_reshape(flat, ncols))
    return str(L.discriminant())


def lattice_is_unimodular(flat, ncols):
    L = IntegerLattice(_reshape(flat, ncols))
    return _fmt_bool(bool(L.is_unimodular()))


FUNCTIONS = {
    # free_module.py
    'zz_span': zz_span,
    'zz_span_qq_entries': zz_span_qq_entries,
    'zz_span_of_basis': zz_span_of_basis,
    'qq_span': qq_span,
    'qq_span_qq_entries': qq_span_qq_entries,
    'gf_span': gf_span,
    'zz_intersection': zz_intersection,
    'qq_intersection': qq_intersection,
    'gf_intersection': gf_intersection,
    'zz_sum': zz_sum,
    'zz_discriminant': zz_discriminant,
    'zz_saturation': zz_saturation,
    'zz_index_in': zz_index_in,
    'zz_index_in_ambient': zz_index_in_ambient,
    'zz_coordinates': zz_coordinates,
    'zz_coordinates_user_basis': zz_coordinates_user_basis,
    'zz_is_submodule': zz_is_submodule,
    'zz_module_eq': zz_module_eq,
    'zz_quotient_invariants': zz_quotient_invariants,
    'qq_complement': qq_complement,
    'gf_complement': gf_complement,
    'gf_cardinality': gf_cardinality,
    # free_module_element.pyx
    'vector_ops_zz': vector_ops_zz,
    'vector_cross_zz': vector_cross_zz,
    'vector_ops_gf': vector_ops_gf,
    'vector_scalar_zz': vector_scalar_zz,
    # free_module_integer.py / bkz
    'lattice_rank_degree': lattice_rank_degree,
    'lattice_basis': lattice_basis,
    'lattice_lll_fingerprint': lattice_lll_fingerprint,
    'lll_exact': lll_exact,
    'lll_invariants': lll_invariants,
    'bkz_exact': bkz_exact,
    'bkz_invariants': bkz_invariants,
    'hkz_first_norm2': hkz_first_norm2,
    'shortest_vector_norm2': shortest_vector_norm2,
    'closest_vector_dist2': closest_vector_dist2,
    'lattice_volume': lattice_volume,
    'lattice_discriminant': lattice_discriminant,
    'lattice_is_unimodular': lattice_is_unimodular,
}
