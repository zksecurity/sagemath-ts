"""SageMath side of the ``mpfr`` property-test area.

Covers ``sage.rings.real_mpfr`` (``RealField`` / ``RealNumber``) and
``sage.rings.complex_mpfr`` (``ComplexField`` / ``ComplexNumber``).

Cases: tests/property/cases/mpfr.cases.json
TypeScript counterpart: tests/property/typescript/areas/mpfr.ts

Design
------
The shared property-test runners can only generate *integer* arguments, so
every case here takes **indices** into the two literal tables below
(``REAL_VALUES`` / ``COMPLEX_VALUES``).  The tables must stay byte-identical
with the ones in the TypeScript counterpart -- both files list them in the
same order and both build their operands by *parsing the same decimal
string*, which is correctly rounded to the same 53-bit value on both sides.

Every function returns an **already formatted string** (see the README's
"Result Formatting Rules"), produced by one of:

``_fmt(x, EXACT)``
    17 significant decimal digits -- enough to round-trip an IEEE double
    exactly, so the comparison is bit-for-bit.  Used for operations MPFR
    computes with a correctly-rounded primitive that JavaScript also has
    (``+ - * /``, ``sqrt``, ``hypot``-free decompositions, ...).

``_fmt(x, DISPLAY)``
    15 significant decimal digits -- exactly what SageMath *prints* for a
    53-bit ``RealField`` element.  Used for the transcendental and special
    functions, where SageMath is correctly rounded (MPFR/PARI) but the port
    is limited to the host libm.  A disagreement at 15 digits is a real
    defect, not last-bit noise.

Exceptions are caught and rendered as ``ERROR: <message>``, so that the
comparison harness checks the *message* too (``compare.ts`` scores
"both sides raised" as a pass regardless of why, which would hide
divergent error behaviour).
"""

from decimal import ROUND_HALF_UP, Context, Decimal
from math import copysign, isinf, isnan

from sage.all import ComplexField, RealField
from sage.rings.complex_mpfr import ComplexNumber as _ComplexNumber

# --------------------------------------------------------------------------
# Operand tables.  KEEP IN SYNC (same order!) with typescript/areas/mpfr.ts.
# --------------------------------------------------------------------------

#: Real operands, given as decimal strings that both MPFR's ``mpfr_set_str``
#: and JavaScript's ``parseFloat`` round to the same 53-bit value.
REAL_VALUES = [
    '0',                        # 0
    '-0',                       # 1   negative zero
    '1',                        # 2
    '-1',                       # 3
    '0.5',                      # 4
    '-0.5',                     # 5
    '2',                        # 6
    '-2',                       # 7
    '3',                        # 8
    '0.1',                      # 9   not exactly representable
    '0.3333333333333333',       # 10
    '9007199254740992',         # 11  2^53
    '9007199254740993',         # 12  2^53 + 1  (rounds to 2^53)
    '9007199254740991',         # 13  2^53 - 1
    '4503599627370496.5',       # 14  2^52 + 1/2, an exact tie
    '0.49999999999999994',      # 15  the double just below 1/2
    '-0.49999999999999994',     # 16
    '1e-320',                   # 17  subnormal
    '5e-324',                   # 18  smallest subnormal
    '1.7976931348623157e308',   # 19  largest finite double
    '3.141592653589793',        # 20  pi
    '2.718281828459045',        # 21  e
    '1e-16',                    # 22
    '100',                      # 23
    '1e100',                    # 24
    '-1e100',                   # 25
    '0.0001',                   # 26
    '1e16',                     # 27
    '1e17',                     # 28
    '1.5',                      # 29
    '-1.5',                     # 30
    '2.5',                      # 31  tie for round()
    '-2.5',                     # 32  tie for round()
    '1e-8',                     # 33
    '1000000',                  # 34
    '6.02e23',                  # 35
    'NaN',                      # 36
    '+infinity',                # 37
    '-infinity',                # 38
    '-8',                       # 39
    '0.75',                     # 40
    '10',                       # 41
    '-0.75',                    # 42
    '1e-300',                   # 43
    '1e300',                    # 44
    '-3.5',                     # 45
    '20.333333333333332',       # 46  61/3 at 53 bits
    '0.9999999999999999',       # 47  the double just below 1
    '1.0000000000000002',       # 48  the double just above 1
    '-0.9999999999999999',      # 49
    '6',                        # 50
    '0.999',                    # 51
    '1e-5',                     # 52
    '-1e-5',                    # 53
    '12',                       # 54
]

#: Complex operands as (real, imaginary) decimal-string pairs.
COMPLEX_VALUES = [
    ('0', '0'),                 # 0
    ('1', '0'),                 # 1
    ('-1', '0'),                # 2   on the sqrt/log branch cut, from above
    ('-1', '-0'),               # 3   on the branch cut, from below (-0.0)
    ('0', '1'),                 # 4
    ('0', '-1'),                # 5
    ('1', '1'),                 # 6
    ('1', '-1'),                # 7
    ('-1', '1'),                # 8
    ('-1', '-1'),               # 9
    ('-4', '0'),                # 10
    ('-4', '-0'),               # 11
    ('2', '0'),                 # 12  arccos/arcsin branch cut, |z| > 1
    ('-2', '0'),                # 13
    ('0.5', '0'),               # 14
    ('-0.5', '0'),              # 15
    ('2', '1'),                 # 16
    ('-3', '0.5'),              # 17  sqrt "avoid_branch" path, im > 0
    ('-3', '-0.5'),             # 18  sqrt "avoid_branch" path, im < 0
    ('-3', '4'),                # 19  |im| > |re|, re < 0
    ('-3', '-4'),               # 20
    ('1e-100', '1e-100'),       # 21
    ('1', '1e-100'),            # 22
    ('-1', '1e-100'),           # 23  just above the branch cut
    ('-1', '-1e-100'),          # 24  just below the branch cut
    ('0', '2'),                 # 25
    ('0', '-2'),                # 26
    ('3', '-4'),                # 27
    ('1e300', '1e300'),         # 28  a^2 + b^2 overflows a double
    ('1e-200', '1e-200'),       # 29  a^2 + b^2 underflows a double
    ('0.5', '0.5'),             # 30
    ('2', '-0'),                # 31
    ('0', '0.5'),               # 32
    ('0', '-0.5'),              # 33
    ('1.5', '0'),               # 34
    ('-1.5', '0'),              # 35
    ('0.0001', '-0.0001'),      # 36
    ('27', '0'),                # 37
    ('-27', '0'),               # 38
    ('3.141592653589793', '0'), # 39
    ('0', '3.141592653589793'), # 40
    ('1e-8', '1'),              # 41
    ('-1e-8', '1'),             # 42
    ('100', '0.001'),           # 43
    ('-100', '0.001'),          # 44
    ('-100', '-0.001'),         # 45
    ('NaN', '1'),               # 46
    ('1', 'NaN'),               # 47
]

#: Rounding-mode names, indexed to match ``RoundingMode`` in the port.
RND_NAMES = ['RNDN', 'RNDZ', 'RNDD', 'RNDU', 'RNDA']

#: Significant digits for the two formatting regimes (see module docstring).
EXACT = 17
DISPLAY = 15

_R = RealField(53)
_C = ComplexField(53)


# --------------------------------------------------------------------------
# Operand construction and formatting
# --------------------------------------------------------------------------


def _rv(i):
    """The real operand with index ``i``."""
    return _R(REAL_VALUES[int(i)])


def _cv(i):
    """The complex operand with index ``i``."""
    re, im = COMPLEX_VALUES[int(i)]
    return _C(_R(re), _R(im))


#: Rounding used to cut a double down to ``digits`` significant decimals.
#: It must be *ties away from zero*, because that is what ECMA-262
#: ``Number.prototype.toExponential`` mandates ("if there are two such sets of
#: e and n, pick the e and n for which n * 10^(e-f) is larger").  Python's
#: ``'%.*e'`` would round ties to even instead, which really does differ:
#: ``cot(RR(pi))`` is exactly -8165619676597685, whose 15-digit rounding is a
#: tie, and the two rules disagree there.  Getting this wrong makes the two
#: runners report different strings for *identical* doubles.
_DEC = Context(prec=1, rounding=ROUND_HALF_UP)


def _fmt(x, digits):
    """Format a real value with ``digits`` significant decimal digits.

    Signed zeros, NaN and both infinities get their own literal spellings so
    that a sign flip on a branch cut cannot hide inside "0".
    """
    x = float(x)
    if isnan(x):
        return 'NaN'
    if isinf(x):
        return '+infinity' if x > 0 else '-infinity'
    if x == 0:
        return '-0' if copysign(1.0, x) < 0 else '0'
    _DEC.prec = digits
    # Decimal(float) is exact, so this rounds the true value of the double.
    rounded = _DEC.plus(Decimal(x))
    mant, _, exp = format(rounded, '.%de' % (digits - 1)).partition('e')
    return mant + 'e' + str(int(exp))


def _fmtc(z, digits):
    """Format a complex value as ``(re, im)``."""
    return '(%s, %s)' % (_fmt(z.real(), digits), _fmt(z.imag(), digits))


def _fmtv(v, digits):
    """Format a real *or* complex result.

    SageMath sometimes widens the codomain (``RR(-1).sqrt()`` is a
    ``ComplexNumber``); rendering the two shapes differently is what makes
    that visible instead of silently coercing.
    """
    if isinstance(v, _ComplexNumber):
        return _fmtc(v, digits)
    try:
        return _fmt(v, digits)
    except (TypeError, ValueError):
        return str(v)


def _fmtlist(zs, digits):
    return '[' + ', '.join(_fmtc(z, digits) for z in zs) + ']'


def _fmtrat(q):
    """Format a rational the way ``str()`` does: ``p/q``, or ``p`` if q == 1."""
    return str(q)


def _fmtint(n):
    return str(n)


def _fmtbool(b):
    return 'True' if b else 'False'


def _fmtorder(n):
    """Format a multiplicative/additive order (``1``, ``2``, ... or infinite)."""
    return str(n)


# --------------------------------------------------------------------------
# Real field / real number
# --------------------------------------------------------------------------


def real_field_str(rnd):
    return str(RealField(53, rnd=RND_NAMES[int(rnd)]))


def real_field_prec(prec):
    return str(RealField(int(prec)).precision())


def real_field_bad_prec(prec):
    return str(RealField(int(prec)).precision())


def real_field_rounding_mode(rnd):
    return str(RealField(53, rnd=RND_NAMES[int(rnd)]).rounding_mode())


def real_field_characteristic():
    return str(RealField(53).characteristic())


def real_rnd_div(rnd, i, j):
    """``a / b`` in a field with a directed rounding mode."""
    F = RealField(53, rnd=RND_NAMES[int(rnd)])
    return _fmt(F(REAL_VALUES[int(i)]) / F(REAL_VALUES[int(j)]), EXACT)


def real_rnd_sqrt(rnd, i):
    F = RealField(53, rnd=RND_NAMES[int(rnd)])
    return _fmtv(F(REAL_VALUES[int(i)]).sqrt(), EXACT)


def real_str(i):
    return str(_rv(i))


def real_add(i, j):
    return _fmt(_rv(i) + _rv(j), EXACT)


def real_sub(i, j):
    return _fmt(_rv(i) - _rv(j), EXACT)


def real_mul(i, j):
    return _fmt(_rv(i) * _rv(j), EXACT)


def real_div(i, j):
    return _fmt(_rv(i) / _rv(j), EXACT)


def real_neg(i):
    return _fmt(-_rv(i), EXACT)


def real_abs(i):
    return _fmt(abs(_rv(i)), EXACT)


def real_sign(i):
    return str(_rv(i).sign())


def real_pow(i, j):
    return _fmtv(_rv(i) ** _rv(j), EXACT)


def real_sqrt(i):
    return _fmtv(_rv(i).sqrt(), EXACT)


def real_cube_root(i):
    return _fmt(_rv(i).cube_root(), DISPLAY)


def real_nth_root(i, n):
    return _fmt(_rv(i).nth_root(int(n)), DISPLAY)


def real_floor(i):
    return str(_rv(i).floor())


def real_ceil(i):
    return str(_rv(i).ceil())


def real_round(i):
    return str(_rv(i).round())


def real_trunc(i):
    return str(_rv(i).trunc())


def real_frac(i):
    return _fmt(_rv(i).frac(), EXACT)


def real_exact_rational(i):
    return _fmtrat(_rv(i).exact_rational())


def real_simplest_rational(i):
    return _fmtrat(_rv(i).simplest_rational())


def real_nearby_rational_denom(i, d):
    return _fmtrat(_rv(i).nearby_rational(max_denominator=int(d)))


def real_nearby_rational_error(i, j):
    return _fmtrat(_rv(i).nearby_rational(max_error=_rv(j)))


def real_nearby_rational_both(i):
    return _fmtrat(_rv(i).nearby_rational(max_error=_R('0.1'), max_denominator=10))


def real_nearby_rational_none(i):
    return _fmtrat(_rv(i).nearby_rational())


def real_sign_mantissa_exponent(i):
    s, m, e = _rv(i).sign_mantissa_exponent()
    return '(%s, %s, %s)' % (s, m, e)


def real_fp_rank(i):
    return str(_rv(i).fp_rank())


def real_ulp(i):
    return _fmt(_rv(i).ulp(), EXACT)


def real_epsilon(i):
    return _fmt(_rv(i).epsilon(), EXACT)


def real_nextabove(i):
    return _fmt(_rv(i).nextabove(), EXACT)


def real_nextbelow(i):
    return _fmt(_rv(i).nextbelow(), EXACT)


def real_nexttoward(i, j):
    return _fmt(_rv(i).nexttoward(_rv(j)), EXACT)


def real_exp(i):
    return _fmt(_rv(i).exp(), DISPLAY)


def real_exp2(i):
    return _fmt(_rv(i).exp2(), DISPLAY)


def real_exp10(i):
    return _fmt(_rv(i).exp10(), DISPLAY)


def real_expm1(i):
    return _fmt(_rv(i).expm1(), DISPLAY)


def real_log(i):
    return _fmtv(_rv(i).log(), DISPLAY)


def real_log_base(i, b):
    return _fmtv(_rv(i).log(int(b)), DISPLAY)


def real_log2(i):
    return _fmtv(_rv(i).log2(), DISPLAY)


def real_log10(i):
    return _fmtv(_rv(i).log10(), DISPLAY)


def real_log1p(i):
    return _fmtv(_rv(i).log1p(), DISPLAY)


def real_sin(i):
    return _fmt(_rv(i).sin(), DISPLAY)


def real_cos(i):
    return _fmt(_rv(i).cos(), DISPLAY)


def real_tan(i):
    return _fmt(_rv(i).tan(), DISPLAY)


def real_arcsin(i):
    return _fmt(_rv(i).arcsin(), DISPLAY)


def real_arccos(i):
    return _fmt(_rv(i).arccos(), DISPLAY)


def real_arctan(i):
    return _fmt(_rv(i).arctan(), DISPLAY)


def real_sinh(i):
    return _fmt(_rv(i).sinh(), DISPLAY)


def real_cosh(i):
    return _fmt(_rv(i).cosh(), DISPLAY)


def real_tanh(i):
    return _fmt(_rv(i).tanh(), DISPLAY)


def real_arcsinh(i):
    return _fmt(_rv(i).arcsinh(), DISPLAY)


def real_arccosh(i):
    return _fmt(_rv(i).arccosh(), DISPLAY)


def real_arctanh(i):
    return _fmt(_rv(i).arctanh(), DISPLAY)


def real_cot(i):
    return _fmt(_rv(i).cot(), DISPLAY)


def real_sec(i):
    return _fmt(_rv(i).sec(), DISPLAY)


def real_csc(i):
    return _fmt(_rv(i).csc(), DISPLAY)


def real_coth(i):
    return _fmt(_rv(i).coth(), DISPLAY)


def real_sech(i):
    return _fmt(_rv(i).sech(), DISPLAY)


def real_csch(i):
    return _fmt(_rv(i).csch(), DISPLAY)


def real_gamma(i):
    return _fmt(_rv(i).gamma(), DISPLAY)


def real_log_gamma(i):
    return _fmtv(_rv(i).log_gamma(), DISPLAY)


def real_zeta(i):
    return _fmt(_rv(i).zeta(), DISPLAY)


def real_erf(i):
    return _fmt(_rv(i).erf(), DISPLAY)


def real_erfc(i):
    return _fmt(_rv(i).erfc(), DISPLAY)


def real_eint(i):
    return _fmt(_rv(i).eint(), DISPLAY)


def real_agm(i, j):
    return _fmt(_rv(i).agm(_rv(j)), DISPLAY)


def real_j0(i):
    return _fmt(_rv(i).j0(), DISPLAY)


def real_j1(i):
    return _fmt(_rv(i).j1(), DISPLAY)


def real_jn(i, n):
    return _fmt(_rv(i).jn(int(n)), DISPLAY)


def real_y0(i):
    return _fmt(_rv(i).y0(), DISPLAY)


def real_y1(i):
    return _fmt(_rv(i).y1(), DISPLAY)


def real_yn(i, n):
    return _fmt(_rv(i).yn(int(n)), DISPLAY)


def real_is_NaN(i):
    return _fmtbool(_rv(i).is_NaN())


def real_is_infinity(i):
    return _fmtbool(_rv(i).is_infinity())


def real_is_integer(i):
    return _fmtbool(_rv(i).is_integer())


def real_is_square(i):
    return _fmtbool(_rv(i).is_square())


def real_multiplicative_order(i):
    return _fmtorder(_rv(i).multiplicative_order())


# --------------------------------------------------------------------------
# Complex field / complex number
# --------------------------------------------------------------------------


def complex_field_str(prec):
    return str(ComplexField(int(prec)))


def complex_field_prec(prec):
    return str(ComplexField(int(prec)).prec())


def complex_field_is_exact():
    return _fmtbool(ComplexField(53).is_exact())


def complex_field_characteristic():
    return str(ComplexField(53).characteristic())


def complex_field_gen():
    return _fmtc(ComplexField(53).gen(), EXACT)


def complex_field_ngens():
    return str(ComplexField(53).ngens())


def complex_field_zeta(n):
    return _fmtc(ComplexField(53).zeta(int(n)), DISPLAY)


def complex_str(i):
    return str(_cv(i))


def complex_add(i, j):
    return _fmtc(_cv(i) + _cv(j), EXACT)


def complex_sub(i, j):
    return _fmtc(_cv(i) - _cv(j), EXACT)


def complex_mul(i, j):
    return _fmtc(_cv(i) * _cv(j), EXACT)


def complex_div(i, j):
    return _fmtc(_cv(i) / _cv(j), EXACT)


def complex_inv(i):
    return _fmtc(~_cv(i), EXACT)


def complex_neg(i):
    return _fmtc(-_cv(i), EXACT)


def complex_conjugate(i):
    return _fmtc(_cv(i).conjugate(), EXACT)


def complex_abs(i):
    return _fmt(abs(_cv(i)), EXACT)


def complex_norm(i):
    return _fmt(_cv(i).norm(), EXACT)


def complex_argument(i):
    return _fmt(_cv(i).argument(), DISPLAY)


def complex_sqrt(i):
    return _fmtc(_cv(i).sqrt(), EXACT)


def complex_sqrt_all(i):
    return _fmtlist(_cv(i).sqrt(all=True), EXACT)


def complex_nth_root(i, n):
    return _fmtc(_cv(i).nth_root(int(n)), DISPLAY)


def complex_nth_root_all(i, n):
    return _fmtlist(_cv(i).nth_root(int(n), all=True), DISPLAY)


def complex_exp(i):
    return _fmtc(_cv(i).exp(), DISPLAY)


def complex_log(i):
    return _fmtc(_cv(i).log(), DISPLAY)


def complex_log_base(i, b):
    return _fmtc(_cv(i).log(int(b)), DISPLAY)


def complex_cos(i):
    return _fmtc(_cv(i).cos(), DISPLAY)


def complex_sin(i):
    return _fmtc(_cv(i).sin(), DISPLAY)


def complex_tan(i):
    return _fmtc(_cv(i).tan(), DISPLAY)


def complex_cosh(i):
    return _fmtc(_cv(i).cosh(), DISPLAY)


def complex_sinh(i):
    return _fmtc(_cv(i).sinh(), DISPLAY)


def complex_tanh(i):
    return _fmtc(_cv(i).tanh(), DISPLAY)


def complex_arccos(i):
    return _fmtc(_cv(i).arccos(), DISPLAY)


def complex_arcsin(i):
    return _fmtc(_cv(i).arcsin(), DISPLAY)


def complex_arctan(i):
    return _fmtc(_cv(i).arctan(), DISPLAY)


def complex_arccosh(i):
    return _fmtc(_cv(i).arccosh(), DISPLAY)


def complex_arcsinh(i):
    return _fmtc(_cv(i).arcsinh(), DISPLAY)


def complex_arctanh(i):
    return _fmtc(_cv(i).arctanh(), DISPLAY)


def complex_cot(i):
    return _fmtc(_cv(i).cot(), DISPLAY)


def complex_sec(i):
    return _fmtc(_cv(i).sec(), DISPLAY)


def complex_csc(i):
    return _fmtc(_cv(i).csc(), DISPLAY)


def complex_coth(i):
    return _fmtc(_cv(i).coth(), DISPLAY)


def complex_sech(i):
    return _fmtc(_cv(i).sech(), DISPLAY)


def complex_csch(i):
    return _fmtc(_cv(i).csch(), DISPLAY)


def complex_gamma(i):
    return _fmtv(_cv(i).gamma(), DISPLAY)


def complex_zeta(i):
    return _fmtv(_cv(i).zeta(), DISPLAY)


def complex_dilog(i):
    return _fmtv(_cv(i).dilog(), DISPLAY)


def complex_eta(i):
    return _fmtv(_cv(i).eta(), DISPLAY)


def complex_agm(i, j):
    return _fmtv(_cv(i).agm(_cv(j)), DISPLAY)


def complex_is_real(i):
    return _fmtbool(_cv(i).is_real())


def complex_is_imaginary(i):
    return _fmtbool(_cv(i).is_imaginary())


def complex_is_integer(i):
    return _fmtbool(_cv(i).is_integer())


def complex_is_square(i):
    return _fmtbool(_cv(i).is_square())


def complex_is_NaN(i):
    return _fmtbool(_cv(i).is_NaN())


def complex_multiplicative_order(i):
    return _fmtorder(_cv(i).multiplicative_order())


def complex_additive_order(i):
    return _fmtorder(_cv(i).additive_order())


def complex_algdep(i, n):
    return '[' + ', '.join(str(c) for c in _cv(i).algdep(int(n)).list()) + ']'


# --------------------------------------------------------------------------
# Dispatch
# --------------------------------------------------------------------------

_RAW = {
    'real_field_str': real_field_str,
    'real_field_prec': real_field_prec,
    'real_field_bad_prec': real_field_bad_prec,
    'real_field_rounding_mode': real_field_rounding_mode,
    'real_field_characteristic': real_field_characteristic,
    'real_rnd_div': real_rnd_div,
    'real_rnd_sqrt': real_rnd_sqrt,
    'real_str': real_str,
    'real_add': real_add,
    'real_sub': real_sub,
    'real_mul': real_mul,
    'real_div': real_div,
    'real_neg': real_neg,
    'real_abs': real_abs,
    'real_sign': real_sign,
    'real_pow': real_pow,
    'real_sqrt': real_sqrt,
    'real_cube_root': real_cube_root,
    'real_nth_root': real_nth_root,
    'real_floor': real_floor,
    'real_ceil': real_ceil,
    'real_round': real_round,
    'real_trunc': real_trunc,
    'real_frac': real_frac,
    'real_exact_rational': real_exact_rational,
    'real_simplest_rational': real_simplest_rational,
    'real_nearby_rational_denom': real_nearby_rational_denom,
    'real_nearby_rational_error': real_nearby_rational_error,
    'real_nearby_rational_both': real_nearby_rational_both,
    'real_nearby_rational_none': real_nearby_rational_none,
    'real_sign_mantissa_exponent': real_sign_mantissa_exponent,
    'real_fp_rank': real_fp_rank,
    'real_ulp': real_ulp,
    'real_epsilon': real_epsilon,
    'real_nextabove': real_nextabove,
    'real_nextbelow': real_nextbelow,
    'real_nexttoward': real_nexttoward,
    'real_exp': real_exp,
    'real_exp2': real_exp2,
    'real_exp10': real_exp10,
    'real_expm1': real_expm1,
    'real_log': real_log,
    'real_log_base': real_log_base,
    'real_log2': real_log2,
    'real_log10': real_log10,
    'real_log1p': real_log1p,
    'real_sin': real_sin,
    'real_cos': real_cos,
    'real_tan': real_tan,
    'real_arcsin': real_arcsin,
    'real_arccos': real_arccos,
    'real_arctan': real_arctan,
    'real_sinh': real_sinh,
    'real_cosh': real_cosh,
    'real_tanh': real_tanh,
    'real_arcsinh': real_arcsinh,
    'real_arccosh': real_arccosh,
    'real_arctanh': real_arctanh,
    'real_cot': real_cot,
    'real_sec': real_sec,
    'real_csc': real_csc,
    'real_coth': real_coth,
    'real_sech': real_sech,
    'real_csch': real_csch,
    'real_gamma': real_gamma,
    'real_log_gamma': real_log_gamma,
    'real_zeta': real_zeta,
    'real_erf': real_erf,
    'real_erfc': real_erfc,
    'real_eint': real_eint,
    'real_agm': real_agm,
    'real_j0': real_j0,
    'real_j1': real_j1,
    'real_jn': real_jn,
    'real_y0': real_y0,
    'real_y1': real_y1,
    'real_yn': real_yn,
    'real_is_NaN': real_is_NaN,
    'real_is_infinity': real_is_infinity,
    'real_is_integer': real_is_integer,
    'real_is_square': real_is_square,
    'real_multiplicative_order': real_multiplicative_order,
    'complex_field_str': complex_field_str,
    'complex_field_prec': complex_field_prec,
    'complex_field_is_exact': complex_field_is_exact,
    'complex_field_characteristic': complex_field_characteristic,
    'complex_field_gen': complex_field_gen,
    'complex_field_ngens': complex_field_ngens,
    'complex_field_zeta': complex_field_zeta,
    'complex_str': complex_str,
    'complex_add': complex_add,
    'complex_sub': complex_sub,
    'complex_mul': complex_mul,
    'complex_div': complex_div,
    'complex_inv': complex_inv,
    'complex_neg': complex_neg,
    'complex_conjugate': complex_conjugate,
    'complex_abs': complex_abs,
    'complex_norm': complex_norm,
    'complex_argument': complex_argument,
    'complex_sqrt': complex_sqrt,
    'complex_sqrt_all': complex_sqrt_all,
    'complex_nth_root': complex_nth_root,
    'complex_nth_root_all': complex_nth_root_all,
    'complex_exp': complex_exp,
    'complex_log': complex_log,
    'complex_log_base': complex_log_base,
    'complex_cos': complex_cos,
    'complex_sin': complex_sin,
    'complex_tan': complex_tan,
    'complex_cosh': complex_cosh,
    'complex_sinh': complex_sinh,
    'complex_tanh': complex_tanh,
    'complex_arccos': complex_arccos,
    'complex_arcsin': complex_arcsin,
    'complex_arctan': complex_arctan,
    'complex_arccosh': complex_arccosh,
    'complex_arcsinh': complex_arcsinh,
    'complex_arctanh': complex_arctanh,
    'complex_cot': complex_cot,
    'complex_sec': complex_sec,
    'complex_csc': complex_csc,
    'complex_coth': complex_coth,
    'complex_sech': complex_sech,
    'complex_csch': complex_csch,
    'complex_gamma': complex_gamma,
    'complex_zeta': complex_zeta,
    'complex_dilog': complex_dilog,
    'complex_eta': complex_eta,
    'complex_agm': complex_agm,
    'complex_is_real': complex_is_real,
    'complex_is_imaginary': complex_is_imaginary,
    'complex_is_integer': complex_is_integer,
    'complex_is_square': complex_is_square,
    'complex_is_NaN': complex_is_NaN,
    'complex_multiplicative_order': complex_multiplicative_order,
    'complex_additive_order': complex_additive_order,
    'complex_algdep': complex_algdep,
}


def _guard(f):
    """Render an exception as ``ERROR: <message>`` instead of propagating it.

    ``compare.ts`` treats "both sides raised" as a pass without looking at
    the messages, so a raise has to become a *value* for the harness to check
    that SageMath and the port reject the same inputs for the same reason.
    """

    def wrapper(*args):
        try:
            return f(*args)
        except Exception as exc:  # noqa: BLE001 - deliberate: message is the result
            return 'ERROR: %s' % (exc,)

    wrapper.__name__ = getattr(f, '__name__', 'wrapper')
    return wrapper


FUNCTIONS = {name: _guard(f) for name, f in _RAW.items()}
