"""SageMath side of the ``coding_crypto`` property-test area.

Cases: tests/property/cases/coding_crypto.cases.json
TypeScript counterpart: tests/property/typescript/areas/coding_crypto.ts

Covers ``sage.coding.*`` (Reed-Muller, BCH, generalized Reed-Solomon, Goppa)
and ``sage.crypto.*`` (S-boxes, Boolean functions, LWE parameter derivation,
hard-lattice generation).

Conventions (the TypeScript module documents the same ones -- both must agree):

- every function returns an **already formatted string**, so the generic
  ``format_result`` never has to guess a shape;
- lists render as ``[a, b, c]``, matrices as ``[[a, b], [c, d]]``, booleans as
  ``True`` / ``False``;
- a call that is expected to raise on both sides is wrapped in ``_guard``,
  which renders any exception as the single token ``ERROR`` (exception
  *messages* are deliberately not compared);
- finite field elements are passed in and out as the integer whose base-p
  digits are the element's coefficients in the polynomial basis
  (``F.from_integer`` / ``x.to_integer()``).

Port defects this area currently pins (28 of 926 cases fail; every expected
value below was produced by running ``sage``, so the fix belongs in the port,
never in the case):

======================================  =====  =========================================
symptom                                 cases  root cause
======================================  =====  =========================================
``sbox_*`` on ``SBox([0,0,0,0])`` and   15     ``crypto/sbox.ts:228`` computes the output
``SBox([0, 2**49])``                           size as ``max(1, ceil(log2(max+1)))``;
                                               SageMath uses ``ZZ(max(S)).nbits()``
                                               (``sbox.pyx:208``), which is 0 for the
                                               constant S-box and exact above 2^53
``sbox_is_involution`` on non-          5      ``crypto/sbox.ts:383`` returns ``False``;
permutations                                   SageMath's ``is_involution`` is
                                               ``self == self.inverse()``
                                               (``sbox.pyx:1914``) and ``inverse``
                                               raises ``TypeError`` (``:1811``)
``rm_parameters_formula(2, 60)``,       2      ``coding/reed_muller_code.ts:232`` keeps
``(0, 63)``                                    the length in a JS ``number``, so it
                                               renders as 1152921504606847000
``goppa_generator_matrix``              4      ``coding/goppa_code.ts:451`` returns a
                                               non-echelon basis; SageMath returns
                                               ``from_parity_check_matrix(H)
                                               .generator_matrix()``
                                               (``goppa_code.py:434-437``)
``lattice_gen_random`` with q >= 2^31   2      ``crypto/lattice.ts:196`` always draws a
                                               single 31-bit ``c_random()``; SageMath
                                               only does that for the small-modulus
                                               dense matrix templates
======================================  =====  =========================================
"""

from sage.all import *
from sage.crypto.boolean_function import BooleanFunction
from sage.crypto.sbox import SBox, feistel_construction, misty_construction
from sage.crypto.lattice import gen_lattice
from sage.crypto.lwe import LindnerPeikert, Regev, RingLindnerPeikert


# ---------------------------------------------------------------------------
# formatting helpers (mirrored in the TypeScript area module)
# ---------------------------------------------------------------------------

def _lst(xs):
    return '[' + ', '.join(str(x) for x in xs) + ']'


def _mat(rows):
    return '[' + ', '.join(_lst(r) for r in rows) + ']'


def _bool(b):
    return 'True' if b else 'False'


def _guard(f):
    try:
        return f()
    except Exception:
        return 'ERROR'


def _rows(M):
    return [[str(M[i, j]) for j in range(M.ncols())] for i in range(M.nrows())]


def _ints(xs):
    return [int(x) for x in xs]


# ---------------------------------------------------------------------------
# sage.crypto.sbox
# ---------------------------------------------------------------------------

def _sbox(S, big_endian=1):
    return SBox(_ints(S), big_endian=bool(int(big_endian)))


def sbox_sizes(S):
    return _guard(lambda: _lst([_sbox(S).input_size(), _sbox(S).output_size()]))


def sbox_ddt(S):
    return _guard(lambda: _mat(_rows(_sbox(S).difference_distribution_table())))


def sbox_differential_uniformity(S):
    return _guard(lambda: str(_sbox(S).differential_uniformity()))


def sbox_is_apn(S):
    return _guard(lambda: _bool(_sbox(S).is_apn()))


_LAT_SCALES = ['bias', 'correlation', 'absolute_bias', 'fourier_coefficient']


def sbox_lat(S, scale):
    return _guard(lambda: _mat(_rows(
        _sbox(S).linear_approximation_table(scale=_LAT_SCALES[int(scale)]))))


def sbox_linearity(S):
    return _guard(lambda: str(_sbox(S).linearity()))


def sbox_max_linear_bias(S):
    return _guard(lambda: str(_sbox(S).maximal_linear_bias_absolute()))


def sbox_nonlinearity(S):
    return _guard(lambda: str(_sbox(S).nonlinearity()))


def sbox_max_degree(S):
    return _guard(lambda: str(_sbox(S).max_degree()))


def sbox_min_degree(S):
    return _guard(lambda: str(_sbox(S).min_degree()))


def sbox_is_permutation(S):
    return _guard(lambda: _bool(_sbox(S).is_permutation()))


def sbox_is_involution(S):
    return _guard(lambda: _bool(_sbox(S).is_involution()))


def sbox_is_balanced(S):
    return _guard(lambda: _bool(_sbox(S).is_balanced()))


def sbox_inverse(S):
    return _guard(lambda: _lst(list(_sbox(S).inverse())))


def sbox_derivative(S, u):
    return _guard(lambda: _lst(list(_sbox(S).derivative(int(u)))))


def sbox_component_function(S, b):
    return _guard(lambda: _lst(list(
        _sbox(S).component_function(int(b)).truth_table(format='int'))))


def sbox_fixed_points(S):
    return _guard(lambda: _lst(_sbox(S).fixed_points()))


def sbox_differential_branch_number(S):
    return _guard(lambda: str(_sbox(S).differential_branch_number()))


def sbox_linear_branch_number(S):
    return _guard(lambda: str(_sbox(S).linear_branch_number()))


def sbox_call_bits(S, x, big_endian):
    def go():
        s = _sbox(S, big_endian)
        return _lst(s(s.to_bits(int(x), s.input_size())))
    return _guard(go)


def sbox_to_bits(S, x, n, big_endian):
    return _guard(lambda: _lst(_sbox(S, big_endian).to_bits(int(x), int(n))))


def sbox_from_bits(S, bits, big_endian):
    return _guard(lambda: str(_sbox(S, big_endian).from_bits(_ints(bits))))


def _boxes(a, b, c):
    out = [_sbox(a), _sbox(b)]
    if len(c) > 0:
        out.append(_sbox(c))
    return out


def sbox_feistel(a, b, c):
    return _guard(lambda: _lst(list(feistel_construction(_boxes(a, b, c)))))


def sbox_misty(a, b, c):
    return _guard(lambda: _lst(list(misty_construction(_boxes(a, b, c)))))


def sbox_misty_stats(a, b, c):
    def go():
        S = misty_construction(_boxes(a, b, c))
        return _lst([S.differential_uniformity(), S.linearity()])
    return _guard(go)


# ---------------------------------------------------------------------------
# sage.crypto.boolean_function
# ---------------------------------------------------------------------------

def _bf(T):
    return BooleanFunction(_ints(T))


def bf_nvariables(T):
    return _guard(lambda: str(_bf(T).nvariables()))


def bf_truth_table_hex(T):
    return _guard(lambda: str(_bf(T).truth_table(format='hex')))


def bf_from_hex(value, width):
    def go():
        hexstr = format(int(value), '0%dx' % int(width))
        return _lst(list(BooleanFunction(hexstr).truth_table(format='int')))
    return _guard(go)


def bf_walsh(T):
    return _guard(lambda: _lst(_bf(T).walsh_hadamard_transform()))


def bf_absolute_walsh_spectrum(T):
    def go():
        spec = _bf(T).absolute_walsh_spectrum()
        return _lst(['(%s, %s)' % (k, spec[k]) for k in sorted(spec)])
    return _guard(go)


def bf_nonlinearity(T):
    return _guard(lambda: str(_bf(T).nonlinearity()))


def bf_is_bent(T):
    return _guard(lambda: _bool(_bf(T).is_bent()))


def bf_is_balanced(T):
    return _guard(lambda: _bool(_bf(T).is_balanced()))


def bf_correlation_immunity(T):
    return _guard(lambda: str(_bf(T).correlation_immunity()))


def bf_resiliency_order(T):
    return _guard(lambda: str(_bf(T).resiliency_order()))


def bf_autocorrelation(T):
    return _guard(lambda: _lst(_bf(T).autocorrelation()))


def bf_absolute_indicator(T):
    return _guard(lambda: str(_bf(T).absolute_indicator()))


def bf_sum_of_square_indicator(T):
    return _guard(lambda: str(_bf(T).sum_of_square_indicator()))


def bf_is_plateaued(T):
    return _guard(lambda: _bool(_bf(T).is_plateaued()))


def bf_algebraic_degree(T):
    return _guard(lambda: str(_bf(T).algebraic_degree()))


def bf_algebraic_immunity(T):
    return _guard(lambda: str(_bf(T).algebraic_immunity()))


def bf_anf_coefficients(T):
    """ANF as a coefficient vector indexed by monomial support (x_i <-> bit i)."""
    def go():
        f = _bf(T)
        n = f.nvariables()
        coeffs = [0] * (1 << n)
        p = f.algebraic_normal_form()
        for mono in p.monomials():
            idx = 0
            for v in mono.variables():
                idx += 1 << v.index()
            coeffs[idx] = 1
        return _lst(coeffs)
    return _guard(go)


def bf_derivative(T, u):
    return _guard(lambda: _lst(list(_bf(T).derivative(int(u)).truth_table(format='int'))))


def bf_is_linear_structure(T, u):
    return _guard(lambda: _bool(_bf(T).is_linear_structure(int(u))))


def bf_is_linear_structure_vec(T, u):
    return _guard(lambda: _bool(_bf(T).is_linear_structure(_ints(u))))


def bf_has_linear_structure(T):
    return _guard(lambda: _bool(_bf(T).has_linear_structure()))


def bf_is_symmetric(T):
    return _guard(lambda: _bool(_bf(T).is_symmetric()))


def bf_complement(T):
    return _guard(lambda: _lst(list((~_bf(T)).truth_table(format='int'))))


def bf_add(A, B):
    return _guard(lambda: _lst(list((_bf(A) + _bf(B)).truth_table(format='int'))))


def bf_mul(A, B):
    return _guard(lambda: _lst(list((_bf(A) * _bf(B)).truth_table(format='int'))))


def bf_concatenate(A, B):
    return _guard(lambda: _lst(list((_bf(A) | _bf(B)).truth_table(format='int'))))


def bf_zero_function(n):
    def go():
        f = BooleanFunction(int(n))
        return _lst([f.nvariables()] + list(f.truth_table(format='int')))
    return _guard(go)


# ---------------------------------------------------------------------------
# sage.coding.reed_muller_code
# ---------------------------------------------------------------------------

def _rm(r, m):
    return codes.ReedMullerCode(GF(2), int(r), int(m))


def rm_parameters(r, m):
    def go():
        C = _rm(r, m)
        return _lst([C.length(), C.dimension(), C.minimum_distance()])
    return _guard(go)


def rm_parameters_formula(r, m):
    """(length, dimension, minimum distance) for a Reed-Muller code too long to build.

    SageMath 10.3's ``AbstractLinearCode.__init__`` eagerly builds an element of
    the ambient `\\GF{2}^{2^m}`, so ``codes.ReedMullerCode(GF(2), 16, 32)``
    *segfaults* on this installation.  Upstream fixed exactly that (issue #33229,
    see ``reference/sage/src/sage/coding/reed_muller_code.py:406-414``, whose
    doctest gives ``C.dimension(), C.length() == (2448023843, 4294967296)``).

    So for these lengths the oracle is the arithmetic that upstream's own class
    performs -- ``2**num_of_var`` (``:448``), ``_binomial_sum(num_of_var, order)``
    (``:451``, reproduced verbatim from ``:46-66``) and ``2**(m - r)`` (``:492``)
    -- rather than the code object, and the port is still exercised through its
    ``ReedMullerCode`` class.  This is the M124 regression: the port used to
    compute the length with a 32-bit ``1 << m``.
    """
    def go():
        r_, m_ = int(r), int(m)
        if m_ < r_:
            raise ValueError('The order must be less than or equal to %s' % m_)
        s, nCi = Integer(1), Integer(1)
        for i in range(r_):
            nCi = ((m_ - i) * nCi) // (i + 1)
            s = nCi + s
        return _lst([Integer(2) ** m_, s, Integer(2) ** (m_ - r_)])
    return _guard(go)


def rm_generator_matrix(r, m):
    return _guard(lambda: _mat(_rows(_rm(r, m).generator_matrix())))


def rm_encode(r, m, msg):
    def go():
        C = _rm(r, m)
        return _lst(list(C.encode(vector(GF(2), _ints(msg)))))
    return _guard(go)


def _rm_received(C, msg, errs):
    v = list(C.encode(vector(GF(2), _ints(msg))))
    for e in _ints(errs):
        v[e] += 1
    return vector(GF(2), v)


def _rm_decode(C, v):
    """Nearest codeword to ``v`` inside the unique-decoding radius.

    ``decode_to_code`` would default to the ``Syndrome`` decoder, which
    tabulates all `2^{n-k}` coset leaders: 89 s per call for RM(2, 5) (measured;
    the table is rebuilt for every code object, so it does not amortise).  The
    ``InformationSet`` decoder is one of the three SageMath offers for this code
    (``C.decoders_available() == ['InformationSet', 'NearestNeighbor',
    'Syndrome']``) and returns the same word in milliseconds: inside the radius
    ``t = (d - 1) / 2`` the closest codeword is unique, so the choice of decoder
    cannot change the answer, only the running time.  Verified equal to the
    Syndrome decoder's output on RM(2, 4), RM(2, 5), RM(3, 5), RM(1, 4), RM(0, 3)
    and RM(1, 3) for every error weight 0 <= w <= t.
    """
    t = (C.minimum_distance() - 1) // 2
    return C.decode_to_code(v, decoder_name='InformationSet', number_errors=t)


def rm_decode_to_message(r, m, msg, errs):
    def go():
        C = _rm(r, m)
        # decode_to_message() forwards **kwargs to the *encoder* as well
        # (abstract_code.py:1001), so number_errors cannot be passed through it;
        # unencode() of the decoded codeword is what it does internally anyway.
        return _lst(list(C.unencode(_rm_decode(C, _rm_received(C, msg, errs)))))
    return _guard(go)


def rm_decode_to_code(r, m, msg, errs):
    def go():
        C = _rm(r, m)
        return _lst(list(_rm_decode(C, _rm_received(C, msg, errs))))
    return _guard(go)


# ---------------------------------------------------------------------------
# sage.coding.bch_code
# ---------------------------------------------------------------------------

def _bch(q, n, delta, b, l):
    return codes.BCHCode(GF(int(q), 'a'), int(n), int(delta),
                         offset=int(b), jump_size=int(l))


def bch_generator_polynomial(q, n, delta, b, l):
    return _guard(lambda: _lst([str(c) for c in _bch(q, n, delta, b, l)
                                .generator_polynomial().list()]))


def bch_defining_set(q, n, delta, b, l):
    return _guard(lambda: _lst(sorted(_bch(q, n, delta, b, l).defining_set())))


def bch_dimension(q, n, delta, b, l):
    return _guard(lambda: str(_bch(q, n, delta, b, l).dimension()))


def bch_decode(q, n, delta, b, l, msg, errs):
    def go():
        F = GF(int(q))
        C = _bch(q, n, delta, b, l)
        g = C.generator_polynomial()
        R = g.parent()
        cw = (R([F(int(c)) for c in msg]) * g).list()
        cw = cw + [F(0)] * (int(n) - len(cw))
        errs_i = _ints(errs)
        for i in range(0, len(errs_i), 2):
            cw[errs_i[i]] += F(errs_i[i + 1])
        return _lst([int(x) for x in C.decode_to_code(vector(F, cw))])
    return _guard(go)


# ---------------------------------------------------------------------------
# sage.coding.grs_code
# ---------------------------------------------------------------------------

def _grs(p, n, k):
    F = GF(int(p))
    pts = [F(i) for i in range(int(n))]
    return F, codes.GeneralizedReedSolomonCode(pts, int(k))


def grs_encode(p, n, k, msg):
    def go():
        F, C = _grs(p, n, k)
        return _lst([int(x) for x in C.encode(vector(F, _ints(msg)))])
    return _guard(go)


def grs_parity_column_multipliers(p, n, k):
    def go():
        _F, C = _grs(p, n, k)
        return _lst([int(x) for x in C.parity_check_matrix().row(0)])
    return _guard(go)


def grs_syndrome(p, n, k, word):
    def go():
        F, C = _grs(p, n, k)
        H = C.parity_check_matrix()
        return _lst([int(x) for x in H * vector(F, _ints(word))])
    return _guard(go)


def _grs_received(F, C, msg, errs):
    v = list(C.encode(vector(F, _ints(msg))))
    errs_i = _ints(errs)
    for i in range(0, len(errs_i), 2):
        v[errs_i[i]] += F(errs_i[i + 1])
    return vector(F, v)


def grs_decode_to_message(p, n, k, msg, errs):
    def go():
        F, C = _grs(p, n, k)
        return _lst([int(x) for x in C.decode_to_message(_grs_received(F, C, msg, errs))])
    return _guard(go)


def grs_decode_to_code(p, n, k, msg, errs):
    def go():
        F, C = _grs(p, n, k)
        return _lst([int(x) for x in C.decode_to_code(_grs_received(F, C, msg, errs))])
    return _guard(go)


# ---------------------------------------------------------------------------
# sage.coding.goppa_code
# ---------------------------------------------------------------------------

def _goppa(q, gcoeffs, support):
    F = GF(int(q), 'a')
    R = PolynomialRing(F, 'x')
    g = R([F.from_integer(int(c)) for c in gcoeffs])
    L = [F.from_integer(int(i)) for i in support]
    return F, codes.GoppaCode(g, L)


def goppa_parity_check_matrix(q, gcoeffs, support):
    return _guard(lambda: _mat(_rows(_goppa(q, gcoeffs, support)[1].parity_check_matrix())))


def goppa_dimension(q, gcoeffs, support):
    def go():
        _F, C = _goppa(q, gcoeffs, support)
        return _lst([C.length(), C.dimension()])
    return _guard(go)


def goppa_distance_bound(q, gcoeffs, support):
    return _guard(lambda: str(_goppa(q, gcoeffs, support)[1].distance_bound()))


def goppa_generator_matrix(q, gcoeffs, support):
    return _guard(lambda: _mat(_rows(_goppa(q, gcoeffs, support)[1].generator_matrix())))


def goppa_decode(q, gcoeffs, support, codeword, errs):
    def go():
        _F, C = _goppa(q, gcoeffs, support)
        word = _ints(codeword)
        for e in _ints(errs):
            word[e] = 1 - word[e]
        return _lst([int(x) for x in C.decode_to_code(vector(GF(2), word))])
    return _guard(go)


# ---------------------------------------------------------------------------
# sage.crypto.lwe
# ---------------------------------------------------------------------------

def _lwe_params(l, with_secret_dist=True):
    parts = [l.K.order(), '%.6f' % float(l.D.sigma),
             'None' if l.m is None else l.m]
    if with_secret_dist:
        parts.append(l.secret_dist)
    return _lst(parts)


def lwe_regev_params(n):
    return _guard(lambda: _lwe_params(Regev(int(n))))


def lwe_lindner_peikert_params(n):
    return _guard(lambda: _lwe_params(LindnerPeikert(int(n))))


def lwe_ring_lindner_peikert_params(n):
    def go():
        l = RingLindnerPeikert(int(n))
        return _lst([l.K.order(), '%.6f' % float(l.D.D.sigma),
                     'None' if l.m is None else l.m])
    return _guard(go)


# ---------------------------------------------------------------------------
# sage.crypto.lattice
# ---------------------------------------------------------------------------

def lattice_gen_modular(n, m, q, seed):
    return _guard(lambda: _mat(_rows(gen_lattice(
        type='modular', n=int(n), m=int(m), q=int(q), seed=int(seed)))))


def lattice_gen_modular_dual(n, m, q, seed):
    return _guard(lambda: _mat(_rows(gen_lattice(
        type='modular', n=int(n), m=int(m), q=int(q), seed=int(seed), dual=True))))


def lattice_gen_random(m, q, seed):
    return _guard(lambda: _mat(_rows(gen_lattice(
        type='random', n=1, m=int(m), q=int(q), seed=int(seed)))))


def _basis_invariants(B, q):
    qq = int(q)
    q_rows = 0
    for i in range(B.nrows()):
        if all(int(x) in (0, qq, -qq) for x in B.row(i)):
            q_rows += 1
    return _lst([B.nrows(), B.ncols(), abs(B.determinant()), q_rows])


def lattice_gen_invariants(n, m, q, seed, dual):
    return _guard(lambda: _basis_invariants(gen_lattice(
        type='modular', n=int(n), m=int(m), q=int(q), seed=int(seed),
        dual=bool(int(dual))), q))


def lattice_gen_cyclotomic_invariants(n, m, q, seed):
    return _guard(lambda: _basis_invariants(gen_lattice(
        type='cyclotomic', n=int(n), m=int(m), q=int(q), seed=int(seed)), q))


FUNCTIONS = {
    # sage.crypto.sbox
    'sbox_sizes': sbox_sizes,
    'sbox_ddt': sbox_ddt,
    'sbox_differential_uniformity': sbox_differential_uniformity,
    'sbox_is_apn': sbox_is_apn,
    'sbox_lat': sbox_lat,
    'sbox_linearity': sbox_linearity,
    'sbox_max_linear_bias': sbox_max_linear_bias,
    'sbox_nonlinearity': sbox_nonlinearity,
    'sbox_max_degree': sbox_max_degree,
    'sbox_min_degree': sbox_min_degree,
    'sbox_is_permutation': sbox_is_permutation,
    'sbox_is_involution': sbox_is_involution,
    'sbox_is_balanced': sbox_is_balanced,
    'sbox_inverse': sbox_inverse,
    'sbox_derivative': sbox_derivative,
    'sbox_component_function': sbox_component_function,
    'sbox_fixed_points': sbox_fixed_points,
    'sbox_differential_branch_number': sbox_differential_branch_number,
    'sbox_linear_branch_number': sbox_linear_branch_number,
    'sbox_call_bits': sbox_call_bits,
    'sbox_to_bits': sbox_to_bits,
    'sbox_from_bits': sbox_from_bits,
    'sbox_feistel': sbox_feistel,
    'sbox_misty': sbox_misty,
    'sbox_misty_stats': sbox_misty_stats,
    # sage.crypto.boolean_function
    'bf_nvariables': bf_nvariables,
    'bf_truth_table_hex': bf_truth_table_hex,
    'bf_from_hex': bf_from_hex,
    'bf_walsh': bf_walsh,
    'bf_absolute_walsh_spectrum': bf_absolute_walsh_spectrum,
    'bf_nonlinearity': bf_nonlinearity,
    'bf_is_bent': bf_is_bent,
    'bf_is_balanced': bf_is_balanced,
    'bf_correlation_immunity': bf_correlation_immunity,
    'bf_resiliency_order': bf_resiliency_order,
    'bf_autocorrelation': bf_autocorrelation,
    'bf_absolute_indicator': bf_absolute_indicator,
    'bf_sum_of_square_indicator': bf_sum_of_square_indicator,
    'bf_is_plateaued': bf_is_plateaued,
    'bf_algebraic_degree': bf_algebraic_degree,
    'bf_algebraic_immunity': bf_algebraic_immunity,
    'bf_anf_coefficients': bf_anf_coefficients,
    'bf_derivative': bf_derivative,
    'bf_is_linear_structure': bf_is_linear_structure,
    'bf_is_linear_structure_vec': bf_is_linear_structure_vec,
    'bf_has_linear_structure': bf_has_linear_structure,
    'bf_is_symmetric': bf_is_symmetric,
    'bf_complement': bf_complement,
    'bf_add': bf_add,
    'bf_mul': bf_mul,
    'bf_concatenate': bf_concatenate,
    'bf_zero_function': bf_zero_function,
    # sage.coding.reed_muller_code
    'rm_parameters': rm_parameters,
    'rm_parameters_formula': rm_parameters_formula,
    'rm_generator_matrix': rm_generator_matrix,
    'rm_encode': rm_encode,
    'rm_decode_to_message': rm_decode_to_message,
    'rm_decode_to_code': rm_decode_to_code,
    # sage.coding.bch_code
    'bch_generator_polynomial': bch_generator_polynomial,
    'bch_defining_set': bch_defining_set,
    'bch_dimension': bch_dimension,
    'bch_decode': bch_decode,
    # sage.coding.grs_code
    'grs_encode': grs_encode,
    'grs_parity_column_multipliers': grs_parity_column_multipliers,
    'grs_syndrome': grs_syndrome,
    'grs_decode_to_message': grs_decode_to_message,
    'grs_decode_to_code': grs_decode_to_code,
    # sage.coding.goppa_code
    'goppa_parity_check_matrix': goppa_parity_check_matrix,
    'goppa_dimension': goppa_dimension,
    'goppa_distance_bound': goppa_distance_bound,
    'goppa_generator_matrix': goppa_generator_matrix,
    'goppa_decode': goppa_decode,
    # sage.crypto.lwe
    'lwe_regev_params': lwe_regev_params,
    'lwe_lindner_peikert_params': lwe_lindner_peikert_params,
    'lwe_ring_lindner_peikert_params': lwe_ring_lindner_peikert_params,
    # sage.crypto.lattice
    'lattice_gen_modular': lattice_gen_modular,
    'lattice_gen_modular_dual': lattice_gen_modular_dual,
    'lattice_gen_random': lattice_gen_random,
    'lattice_gen_invariants': lattice_gen_invariants,
    'lattice_gen_cyclotomic_invariants': lattice_gen_cyclotomic_invariants,
}
