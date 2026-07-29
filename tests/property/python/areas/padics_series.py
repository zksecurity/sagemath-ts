"""SageMath side of the ``padics_series`` property-test area.

Covers ``sage.rings.padics`` (``Zp``/``Qp`` and their elements),
``sage.rings.power_series_ring`` and ``sage.rings.laurent_series_ring``.

Cases: tests/property/cases/padics_series.cases.json
sagemath-ts counterpart: tests/property/typescript/areas/padics_series.ts

Argument convention
-------------------
Every function takes one, two or three *lists of integers* (produced by
``fixedValue([...])`` / ``randomList(...)`` in the cases file), because those
are the only shapes the shared runners can generate identically on both sides.
The lists are decoded below; ``NONE`` (= -999999) is the sentinel for a Python
``None`` / TypeScript ``undefined`` argument.

Output convention
-----------------
Every function returns a plain ``str``.  Both sides build that string with the
same helpers so that any difference is a genuine behavioural difference:

* p-adic element  -> ``<repr> [v=<valuation>, rp=<rel prec>, ap=<abs prec>]``
* (Laurent) power series -> ``<repr> [v=..., prec=..., deg=...]``
* an exception    -> ``ERROR:<ExceptionClassName>: <message>``

The exception form is deliberate: ``compare.ts`` scores "both sides raised" as
a pass without looking at *why*, so raising is folded into the *result* string
instead.  That also pins SageMath's exact error messages, which CLAUDE.md
requires the port to preserve.  No case in this area is a vacuous pass: the
runner never reports an error for any of them, so every one of them compares a
value.

Two constraints shaped the case list; do not undo them without re-checking:

* ``randomBigint(min, max)`` is only trustworthy while ``max - min + 1`` fits in
  32 bits.  ``tests/property/typescript/mersenne-twister.ts:getrandbits`` packs
  the 32-bit words in the opposite order to CPython's ``_random_getrandbits``,
  so wider ranges make the two runners generate *different* arguments (which
  ``compare.ts`` then silently reports as "not found in TypeScript results").
  Large inputs -- values astride 2^53, 2^64, 2^80 -- are therefore pinned with
  ``fixedValue`` instead.
* The oracle is whatever ``sage`` is installed (10.3 here), which is older than
  the vendored ``reference/sage`` (10.9.beta4).  ``LaurentSeries.is_square``
  exists in the vendored source and in the port but not in 10.3, so it has no
  cases here; add them when the installed oracle catches up.
"""

from sage.all import (
    GF,
    QQ,
    Integer,
    LaurentSeriesRing,
    PowerSeriesRing,
    Qp,
    Zp,
)

# Sentinel meaning "argument omitted" (Python ``None``).
NONE = -999999


# ---------------------------------------------------------------------------
# formatting helpers
# ---------------------------------------------------------------------------


def _guard(fn):
    """Run ``fn`` and turn any exception into a comparable result string."""
    try:
        return fn()
    except Exception as e:  # noqa: BLE001 - the point is to compare the exception
        return 'ERROR:%s: %s' % (type(e).__name__, e)


def _opt(v):
    """Decode an optional integer argument."""
    v = int(v)
    return None if v == NONE else v


def _fmt_padic(x):
    return '%s [v=%s, rp=%s, ap=%s]' % (
        x,
        x.valuation(),
        x.precision_relative(),
        x.precision_absolute(),
    )


def _fmt_series(f):
    return '%s [v=%s, prec=%s, deg=%s]' % (f, f.valuation(), f.prec(), f.degree())


def _fmt_list(items):
    return '[' + ', '.join(str(x) for x in items) + ']'


# ---------------------------------------------------------------------------
# p-adic helpers
# ---------------------------------------------------------------------------


def _padic_ring(p, prec, is_field):
    p = Integer(p)
    prec = int(prec)
    return Qp(p, prec) if int(is_field) else Zp(p, prec)


def _padic_elt(spec):
    """Decode ``[p, prec, is_field, num, den]`` into ``(R, R(num)/R(den))``.

    Building the element as a quotient of two ring elements (rather than
    coercing a rational) keeps both sides on exactly the same code path.
    """
    p, prec, is_field, num, den = (Integer(x) for x in spec[:5])
    R = _padic_ring(p, prec, is_field)
    if den == 1:
        return R, R(num)
    return R, R(num) / R(den)


# ---------------------------------------------------------------------------
# p-adic rings
# ---------------------------------------------------------------------------


def padic_ring_info(spec):
    """``[p, prec, is_field]`` -> the ring's own invariants."""

    def go():
        p, prec, is_field = spec[:3]
        R = _padic_ring(p, prec, is_field)
        return ' | '.join(
            [
                str(R),
                'p=%s' % R.prime(),
                'cap=%s' % R.precision_cap(),
                'field=%s' % R.is_field(),
                'char=%s' % R.characteristic(),
                'res_char=%s' % R.residue_characteristic(),
                'e=%s' % R.absolute_e(),
                'f=%s' % R.absolute_f(),
                'unif=%s' % R.uniformizer(),
                'zero=%s' % _fmt_padic(R.zero()),
                'one=%s' % _fmt_padic(R.one()),
            ]
        )

    return _guard(go)


def padic_uniformizer_pow(spec):
    """``[p, prec, is_field, k]`` -> ``R.uniformizer_pow(k)``."""

    def go():
        R = _padic_ring(spec[0], spec[1], spec[2])
        return _fmt_padic(R.uniformizer_pow(Integer(spec[3])))

    return _guard(go)


def padic_teichmuller(spec):
    """``[p, prec, is_field, x]`` -> ``R.teichmuller(x)``."""

    def go():
        R = _padic_ring(spec[0], spec[1], spec[2])
        return _fmt_padic(R.teichmuller(Integer(spec[3])))

    return _guard(go)


def padic_teichmuller_system(spec):
    """``[p, prec, is_field]`` -> ``R.teichmuller_system()``."""

    def go():
        R = _padic_ring(spec[0], spec[1], spec[2])
        return _fmt_list(R.teichmuller_system())

    return _guard(go)


def padic_roots_of_unity(spec):
    """``[p, prec, is_field, n]`` -> ``R.roots_of_unity(n)`` (n = NONE: all)."""

    def go():
        R = _padic_ring(spec[0], spec[1], spec[2])
        n = _opt(spec[3])
        roots = R.roots_of_unity() if n is None else R.roots_of_unity(Integer(n))
        return _fmt_list(roots)

    return _guard(go)


# ---------------------------------------------------------------------------
# p-adic elements
# ---------------------------------------------------------------------------


def padic_repr(spec):
    """``[p, prec, is_field, num, den]`` -> the element itself."""
    return _guard(lambda: _fmt_padic(_padic_elt(spec)[1]))


def padic_parent(spec):
    """``[p, prec, is_field, num, den]`` -> the *parent* of the element.

    Division in ``Zp`` lands in ``Qp`` upstream, so this pins coercion.
    """
    return _guard(lambda: str(_padic_elt(spec)[1].parent()))


def padic_from_int(p, prec, value):
    """Three scalar arguments (used with ``randomBigint``)."""

    def go():
        R = Zp(Integer(p), int(prec))
        return _fmt_padic(R(Integer(value)))

    return _guard(go)


def padic_from_int_field(p, prec, value):
    def go():
        K = Qp(Integer(p), int(prec))
        return _fmt_padic(K(Integer(value)))

    return _guard(go)


def padic_expansion(spec):
    """``[p, prec, is_field, num, den]`` -> ``list(x.expansion())``."""

    def go():
        _R, x = _padic_elt(spec)
        return _fmt_list(list(x.expansion()))

    return _guard(go)


def padic_flags(spec):
    """``[p, prec, is_field, num, den]`` -> the boolean predicates."""

    def go():
        _R, x = _padic_elt(spec)
        return ' | '.join(
            [
                'zero=%s' % x.is_zero(),
                'one=%s' % x.is_one(),
                'unit=%s' % x.is_unit(),
                'integral=%s' % x.is_integral(),
                'val=%s' % x.valuation(),
                'lift=%s' % _guard(lambda: str(x.lift())),
            ]
        )

    return _guard(go)


def padic_unit_part(spec):
    return _guard(lambda: _fmt_padic(_padic_elt(spec)[1].unit_part()))


def padic_add(spec):
    """``[p, prec, is_field, n1, d1, n2, d2]``."""

    def go():
        _R, a = _padic_elt(spec)
        _R2, b = _padic_elt([spec[0], spec[1], spec[2], spec[5], spec[6]])
        return _fmt_padic(a + b)

    return _guard(go)


def padic_sub(spec):
    def go():
        _R, a = _padic_elt(spec)
        _R2, b = _padic_elt([spec[0], spec[1], spec[2], spec[5], spec[6]])
        return _fmt_padic(a - b)

    return _guard(go)


def padic_mul(spec):
    def go():
        _R, a = _padic_elt(spec)
        _R2, b = _padic_elt([spec[0], spec[1], spec[2], spec[5], spec[6]])
        return _fmt_padic(a * b)

    return _guard(go)


def padic_div(spec):
    def go():
        _R, a = _padic_elt(spec)
        _R2, b = _padic_elt([spec[0], spec[1], spec[2], spec[5], spec[6]])
        return _fmt_padic(a / b)

    return _guard(go)


def padic_pow(spec):
    """``[p, prec, is_field, num, den, e]``."""

    def go():
        _R, a = _padic_elt(spec)
        return _fmt_padic(a ** Integer(spec[5]))

    return _guard(go)


def padic_inv(spec):
    def go():
        _R, a = _padic_elt(spec)
        return _fmt_padic(a**-1)

    return _guard(go)


def padic_residue(spec):
    """``[p, prec, is_field, num, den, n]``."""

    def go():
        _R, a = _padic_elt(spec)
        return str(a.residue(int(spec[5])))

    return _guard(go)


def padic_getitem(spec):
    """``[p, prec, is_field, num, den, i]``."""

    def go():
        _R, a = _padic_elt(spec)
        return str(a[int(spec[5])])

    return _guard(go)


def padic_slice(spec):
    """``[p, prec, is_field, num, den, i, j, k]`` (NONE-able i, j, k)."""

    def go():
        _R, a = _padic_elt(spec)
        i, j, k = _opt(spec[5]), _opt(spec[6]), _opt(spec[7])
        return _fmt_padic(a.slice(i, j, k))

    return _guard(go)


def padic_add_bigoh(spec):
    """``[p, prec, is_field, num, den, n]``."""

    def go():
        _R, a = _padic_elt(spec)
        return _fmt_padic(a.add_bigoh(int(spec[5])))

    return _guard(go)


def padic_lift_to_precision(spec):
    """``[p, prec, is_field, num, den, cut, target]``."""

    def go():
        _R, a = _padic_elt(spec)
        return _fmt_padic(a.add_bigoh(int(spec[5])).lift_to_precision(int(spec[6])))

    return _guard(go)


def padic_inexact_zero(spec):
    """``[p, prec, is_field, k]`` -> everything about ``R(0, k)``."""

    def go():
        R = _padic_ring(spec[0], spec[1], spec[2])
        z = R(0, int(spec[3]))
        return ' | '.join(
            [
                _fmt_padic(z),
                'is_zero=%s' % z.is_zero(),
                'exp=%s' % _fmt_list(list(z.expansion())),
                'add_order=%s' % _guard(lambda: str(z.additive_order())),
                'plus_one=%s' % _fmt_padic(z + R(1)),
            ]
        )

    return _guard(go)


def padic_sqrt(spec):
    return _guard(lambda: _fmt_padic(_padic_elt(spec)[1].sqrt()))


def padic_sqrt_all(spec):
    return _guard(lambda: _fmt_list(_padic_elt(spec)[1].sqrt(all=True)))


def padic_is_square(spec):
    return _guard(lambda: str(_padic_elt(spec)[1].is_square()))


def padic_nth_root(spec):
    """``[p, prec, is_field, num, den, n]``."""

    def go():
        _R, a = _padic_elt(spec)
        return _fmt_padic(a.nth_root(Integer(spec[5])))

    return _guard(go)


def padic_nth_root_all(spec):
    def go():
        _R, a = _padic_elt(spec)
        return _fmt_list(a.nth_root(Integer(spec[5]), all=True))

    return _guard(go)


def padic_log(spec):
    return _guard(lambda: _fmt_padic(_padic_elt(spec)[1].log()))


def padic_exp(spec):
    return _guard(lambda: _fmt_padic(_padic_elt(spec)[1].exp()))


def padic_artin_hasse_exp(spec):
    return _guard(lambda: _fmt_padic(_padic_elt(spec)[1].artin_hasse_exp()))


def padic_multiplicative_order(spec):
    return _guard(lambda: str(_padic_elt(spec)[1].multiplicative_order()))


def padic_additive_order(spec):
    return _guard(lambda: str(_padic_elt(spec)[1].additive_order()))


# ---------------------------------------------------------------------------
# power series helpers
# ---------------------------------------------------------------------------


def _base_ring(code):
    """0 -> QQ, otherwise the prime field ``GF(code)``."""
    code = int(code)
    return QQ if code == 0 else GF(Integer(code))


def _ps_ring(params):
    """``params = [base, default_prec, ...]``."""
    return PowerSeriesRing(_base_ring(params[0]), 'x', default_prec=int(params[1]))


def _ps(params, coeffs, prec_index):
    """Build ``R(coeffs, prec)`` with ``prec = params[prec_index]``."""
    R = _ps_ring(params)
    prec = _opt(params[prec_index])
    lst = [Integer(c) for c in coeffs]
    if prec is None:
        return R, R(lst)
    return R, R(lst, prec)


# ---------------------------------------------------------------------------
# power series
# ---------------------------------------------------------------------------


def ps_repr(coeffs, params):
    """``params = [base, default_prec, prec]``."""
    return _guard(lambda: _fmt_series(_ps(params, coeffs, 2)[1]))


def ps_info(coeffs, params):
    def go():
        _R, f = _ps(params, coeffs, 2)
        return ' | '.join(
            [
                _fmt_series(f),
                'ap=%s' % f.precision_absolute(),
                'rp=%s' % f.precision_relative(),
                'zero=%s' % f.is_zero(),
                'one=%s' % f.is_one(),
                'unit=%s' % f.is_unit(),
                'monomial=%s' % f.is_monomial(),
                'list=%s' % _fmt_list(f.list()),
            ]
        )

    return _guard(go)


def ps_add(c1, c2, params):
    """``params = [base, default_prec, prec1, prec2]``."""

    def go():
        _R, f = _ps(params, c1, 2)
        _R2, g = _ps(params, c2, 3)
        return _fmt_series(f + g)

    return _guard(go)


def ps_sub(c1, c2, params):
    def go():
        _R, f = _ps(params, c1, 2)
        _R2, g = _ps(params, c2, 3)
        return _fmt_series(f - g)

    return _guard(go)


def ps_mul(c1, c2, params):
    def go():
        _R, f = _ps(params, c1, 2)
        _R2, g = _ps(params, c2, 3)
        return _fmt_series(f * g)

    return _guard(go)


def ps_div(c1, c2, params):
    def go():
        _R, f = _ps(params, c1, 2)
        _R2, g = _ps(params, c2, 3)
        return _fmt_series(f / g)

    return _guard(go)


def ps_compose(c1, c2, params):
    def go():
        _R, f = _ps(params, c1, 2)
        _R2, g = _ps(params, c2, 3)
        return _fmt_series(f(g))

    return _guard(go)


def ps_inv(coeffs, params):
    return _guard(lambda: _fmt_series(~_ps(params, coeffs, 2)[1]))


def ps_pow(coeffs, params):
    """``params = [base, default_prec, prec, e]``."""

    def go():
        _R, f = _ps(params, coeffs, 2)
        return _fmt_series(f ** Integer(params[3]))

    return _guard(go)


def ps_derivative(coeffs, params):
    return _guard(lambda: _fmt_series(_ps(params, coeffs, 2)[1].derivative()))


def ps_integral(coeffs, params):
    return _guard(lambda: _fmt_series(_ps(params, coeffs, 2)[1].integral()))


def ps_exp(coeffs, params):
    """``params = [base, default_prec, prec, target]``."""

    def go():
        _R, f = _ps(params, coeffs, 2)
        target = _opt(params[3])
        return _fmt_series(f.exp() if target is None else f.exp(target))

    return _guard(go)


def ps_log(coeffs, params):
    def go():
        _R, f = _ps(params, coeffs, 2)
        target = _opt(params[3])
        return _fmt_series(f.log() if target is None else f.log(target))

    return _guard(go)


def ps_sqrt(coeffs, params):
    def go():
        _R, f = _ps(params, coeffs, 2)
        target = _opt(params[3])
        return _fmt_series(f.sqrt() if target is None else f.sqrt(target))

    return _guard(go)


def ps_nth_root(coeffs, params):
    """``params = [base, default_prec, prec, n, target]``."""

    def go():
        _R, f = _ps(params, coeffs, 2)
        n = int(params[3])
        target = _opt(params[4])
        return _fmt_series(f.nth_root(n) if target is None else f.nth_root(n, target))

    return _guard(go)


def ps_pade(coeffs, params):
    """``params = [base, default_prec, prec, m, n]``."""

    def go():
        _R, f = _ps(params, coeffs, 2)
        return str(f.pade(int(params[3]), int(params[4])))

    return _guard(go)


def ps_reverse(coeffs, params):
    """``params = [base, default_prec, prec, precision]``."""

    def go():
        _R, f = _ps(params, coeffs, 2)
        precision = _opt(params[3])
        return _fmt_series(f.reverse() if precision is None else f.reverse(precision))

    return _guard(go)


def ps_V(coeffs, params):
    """``params = [base, default_prec, prec, n]``."""

    def go():
        _R, f = _ps(params, coeffs, 2)
        return _fmt_series(f.V(int(params[3])))

    return _guard(go)


def ps_shift(coeffs, params):
    def go():
        _R, f = _ps(params, coeffs, 2)
        return _fmt_series(f.shift(int(params[3])))

    return _guard(go)


def ps_truncate(coeffs, params):
    def go():
        _R, f = _ps(params, coeffs, 2)
        n = _opt(params[3])
        return str(f.truncate() if n is None else f.truncate(n))

    return _guard(go)


def ps_truncate_powerseries(coeffs, params):
    def go():
        _R, f = _ps(params, coeffs, 2)
        return _fmt_series(f.truncate_powerseries(int(params[3])))

    return _guard(go)


def ps_add_bigoh(coeffs, params):
    def go():
        _R, f = _ps(params, coeffs, 2)
        return _fmt_series(f.add_bigoh(int(params[3])))

    return _guard(go)


def ps_O(coeffs, params):
    def go():
        _R, f = _ps(params, coeffs, 2)
        return _fmt_series(f.O(int(params[3])))

    return _guard(go)


def ps_valuation_zero_part(coeffs, params):
    return _guard(lambda: _fmt_series(_ps(params, coeffs, 2)[1].valuation_zero_part()))


def ps_getitem(coeffs, params):
    def go():
        _R, f = _ps(params, coeffs, 2)
        return str(f[int(params[3])])

    return _guard(go)


def ps_is_square(coeffs, params):
    return _guard(lambda: str(_ps(params, coeffs, 2)[1].is_square()))


# ---------------------------------------------------------------------------
# Laurent series
# ---------------------------------------------------------------------------


def _ls_ring(params):
    return LaurentSeriesRing(_base_ring(params[0]), 'x', default_prec=int(params[1]))


def _ls(params, coeffs, shift_index, prec_index):
    """``L(R(coeffs), n)`` truncated at ``prec`` (both NONE-able)."""
    L = _ls_ring(params)
    R = L.power_series_ring()
    n = int(params[shift_index])
    prec = _opt(params[prec_index])
    f = R([Integer(c) for c in coeffs])
    g = L(f, n)
    return L, (g if prec is None else g.add_bigoh(prec))


def ls_repr(coeffs, params):
    """``params = [base, default_prec, n, prec]``."""
    return _guard(lambda: _fmt_series(_ls(params, coeffs, 2, 3)[1]))


def ls_info(coeffs, params):
    def go():
        _L, f = _ls(params, coeffs, 2, 3)
        return ' | '.join(
            [
                _fmt_series(f),
                'ap=%s' % f.precision_absolute(),
                'rp=%s' % f.precision_relative(),
                'zero=%s' % f.is_zero(),
                'unit=%s' % f.is_unit(),
                'monomial=%s' % f.is_monomial(),
                'exponents=%s' % _fmt_list(f.exponents()),
                'coefficients=%s' % _fmt_list(f.coefficients()),
                'residue=%s' % f.residue(),
                'vzp=%s' % f.valuation_zero_part(),
            ]
        )

    return _guard(go)


def ls_add(c1, c2, params):
    """``params = [base, default_prec, n1, prec1, n2, prec2]``."""

    def go():
        _L, f = _ls(params, c1, 2, 3)
        _L2, g = _ls(params, c2, 4, 5)
        return _fmt_series(f + g)

    return _guard(go)


def ls_sub(c1, c2, params):
    def go():
        _L, f = _ls(params, c1, 2, 3)
        _L2, g = _ls(params, c2, 4, 5)
        return _fmt_series(f - g)

    return _guard(go)


def ls_mul(c1, c2, params):
    def go():
        _L, f = _ls(params, c1, 2, 3)
        _L2, g = _ls(params, c2, 4, 5)
        return _fmt_series(f * g)

    return _guard(go)


def ls_div(c1, c2, params):
    def go():
        _L, f = _ls(params, c1, 2, 3)
        _L2, g = _ls(params, c2, 4, 5)
        return _fmt_series(f / g)

    return _guard(go)


def ls_inv(coeffs, params):
    return _guard(lambda: _fmt_series(_ls(params, coeffs, 2, 3)[1].inverse()))


def ls_pow(coeffs, params):
    """``params = [base, default_prec, n, prec, e]``."""

    def go():
        _L, f = _ls(params, coeffs, 2, 3)
        return _fmt_series(f ** Integer(params[4]))

    return _guard(go)


def ls_derivative(coeffs, params):
    return _guard(lambda: _fmt_series(_ls(params, coeffs, 2, 3)[1].derivative()))


def ls_integral(coeffs, params):
    return _guard(lambda: _fmt_series(_ls(params, coeffs, 2, 3)[1].integral()))


def ls_shift(coeffs, params):
    """``params = [base, default_prec, n, prec, k]``."""

    def go():
        _L, f = _ls(params, coeffs, 2, 3)
        return _fmt_series(f.shift(int(params[4])))

    return _guard(go)


def ls_truncate(coeffs, params):
    def go():
        _L, f = _ls(params, coeffs, 2, 3)
        return str(f.truncate(int(params[4])))

    return _guard(go)


def ls_truncate_laurentseries(coeffs, params):
    def go():
        _L, f = _ls(params, coeffs, 2, 3)
        return _fmt_series(f.truncate_laurentseries(int(params[4])))

    return _guard(go)


def ls_truncate_neg(coeffs, params):
    def go():
        _L, f = _ls(params, coeffs, 2, 3)
        return _fmt_series(f.truncate_neg(int(params[4])))

    return _guard(go)


def ls_add_bigoh(coeffs, params):
    def go():
        _L, f = _ls(params, coeffs, 2, 3)
        return _fmt_series(f.add_bigoh(int(params[4])))

    return _guard(go)


def ls_verschiebung(coeffs, params):
    def go():
        _L, f = _ls(params, coeffs, 2, 3)
        return _fmt_series(f.verschiebung(int(params[4])))

    return _guard(go)


def ls_reverse(coeffs, params):
    def go():
        _L, f = _ls(params, coeffs, 2, 3)
        precision = _opt(params[4])
        return _fmt_series(f.reverse() if precision is None else f.reverse(precision))

    return _guard(go)


def ls_nth_root(coeffs, params):
    """``params = [base, default_prec, n, prec, k, target]``."""

    def go():
        _L, f = _ls(params, coeffs, 2, 3)
        k = int(params[4])
        target = _opt(params[5])
        return _fmt_series(f.nth_root(k) if target is None else f.nth_root(k, target))

    return _guard(go)


def ls_is_square(coeffs, params):
    return _guard(lambda: str(_ls(params, coeffs, 2, 3)[1].is_square()))


def ls_power_series(coeffs, params):
    return _guard(lambda: _fmt_series(_ls(params, coeffs, 2, 3)[1].power_series()))


def ls_lift_to_precision(coeffs, params):
    def go():
        _L, f = _ls(params, coeffs, 2, 3)
        absprec = _opt(params[4])
        return _fmt_series(
            f.lift_to_precision() if absprec is None else f.lift_to_precision(absprec)
        )

    return _guard(go)


def ls_getitem(coeffs, params):
    def go():
        _L, f = _ls(params, coeffs, 2, 3)
        return str(f[int(params[4])])

    return _guard(go)


def ls_compose(c1, c2, params):
    def go():
        _L, f = _ls(params, c1, 2, 3)
        _L2, g = _ls(params, c2, 4, 5)
        return _fmt_series(f(g))

    return _guard(go)


def ls_ring_info(params):
    def go():
        L = _ls_ring(params)
        return ' | '.join(
            [
                str(L),
                'field=%s' % L.is_field(),
                'exact=%s' % L.is_exact(),
                'char=%s' % L.characteristic(),
                'default_prec=%s' % L.default_prec(),
                'gen=%s' % L.gen(),
                'unif=%s' % _guard(lambda: str(L.uniformizer())),
                'zero=%s' % _fmt_series(L.zero()),
                'one=%s' % _fmt_series(L.one()),
            ]
        )

    return _guard(go)


FUNCTIONS = {
    # p-adic rings
    'padic_ring_info': padic_ring_info,
    'padic_uniformizer_pow': padic_uniformizer_pow,
    'padic_teichmuller': padic_teichmuller,
    'padic_teichmuller_system': padic_teichmuller_system,
    'padic_roots_of_unity': padic_roots_of_unity,
    # p-adic elements
    'padic_repr': padic_repr,
    'padic_parent': padic_parent,
    'padic_from_int': padic_from_int,
    'padic_from_int_field': padic_from_int_field,
    'padic_expansion': padic_expansion,
    'padic_flags': padic_flags,
    'padic_unit_part': padic_unit_part,
    'padic_add': padic_add,
    'padic_sub': padic_sub,
    'padic_mul': padic_mul,
    'padic_div': padic_div,
    'padic_pow': padic_pow,
    'padic_inv': padic_inv,
    'padic_residue': padic_residue,
    'padic_getitem': padic_getitem,
    'padic_slice': padic_slice,
    'padic_add_bigoh': padic_add_bigoh,
    'padic_lift_to_precision': padic_lift_to_precision,
    'padic_inexact_zero': padic_inexact_zero,
    'padic_sqrt': padic_sqrt,
    'padic_sqrt_all': padic_sqrt_all,
    'padic_is_square': padic_is_square,
    'padic_nth_root': padic_nth_root,
    'padic_nth_root_all': padic_nth_root_all,
    'padic_log': padic_log,
    'padic_exp': padic_exp,
    'padic_artin_hasse_exp': padic_artin_hasse_exp,
    'padic_multiplicative_order': padic_multiplicative_order,
    'padic_additive_order': padic_additive_order,
    # power series
    'ps_repr': ps_repr,
    'ps_info': ps_info,
    'ps_add': ps_add,
    'ps_sub': ps_sub,
    'ps_mul': ps_mul,
    'ps_div': ps_div,
    'ps_compose': ps_compose,
    'ps_inv': ps_inv,
    'ps_pow': ps_pow,
    'ps_derivative': ps_derivative,
    'ps_integral': ps_integral,
    'ps_exp': ps_exp,
    'ps_log': ps_log,
    'ps_sqrt': ps_sqrt,
    'ps_nth_root': ps_nth_root,
    'ps_pade': ps_pade,
    'ps_reverse': ps_reverse,
    'ps_V': ps_V,
    'ps_shift': ps_shift,
    'ps_truncate': ps_truncate,
    'ps_truncate_powerseries': ps_truncate_powerseries,
    'ps_add_bigoh': ps_add_bigoh,
    'ps_O': ps_O,
    'ps_valuation_zero_part': ps_valuation_zero_part,
    'ps_getitem': ps_getitem,
    'ps_is_square': ps_is_square,
    # Laurent series
    'ls_ring_info': ls_ring_info,
    'ls_repr': ls_repr,
    'ls_info': ls_info,
    'ls_add': ls_add,
    'ls_sub': ls_sub,
    'ls_mul': ls_mul,
    'ls_div': ls_div,
    'ls_inv': ls_inv,
    'ls_pow': ls_pow,
    'ls_derivative': ls_derivative,
    'ls_integral': ls_integral,
    'ls_shift': ls_shift,
    'ls_truncate': ls_truncate,
    'ls_truncate_laurentseries': ls_truncate_laurentseries,
    'ls_truncate_neg': ls_truncate_neg,
    'ls_add_bigoh': ls_add_bigoh,
    'ls_verschiebung': ls_verschiebung,
    'ls_reverse': ls_reverse,
    'ls_nth_root': ls_nth_root,
    'ls_is_square': ls_is_square,
    'ls_power_series': ls_power_series,
    'ls_lift_to_precision': ls_lift_to_precision,
    'ls_getitem': ls_getitem,
    'ls_compose': ls_compose,
}
