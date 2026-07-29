"""SageMath side of the ``quadratic_forms`` property-test area.

Cases: tests/property/cases/quadratic_forms.cases.json
TypeScript counterpart: tests/property/typescript/areas/quadratic_forms.ts

Covers ``sage.quadratic_forms.binary_qf``: Gauss composition (``__mul__`` ->
PARI ``qfbcompraw``), reduction (``reduced_form`` -> PARI ``qfbred`` /
``qfbredsl2`` / Sage's own ``_reduce_indef``), the class-group enumerators
(``BinaryQF_reduced_representatives``), cycles, proper/improper equivalence and
``solve_integer`` (PARI ``qfbsolve`` / ``qfbcornacchia`` plus Sage's elementary
algorithm for square discriminants).

Argument convention
-------------------
Every function takes flat integer lists, because the property-test runners can
only generate ``bigint`` and ``bigint[]``.  A "form list" is a flat list of
triples ``[a1, b1, c1, a2, b2, c2, ...]``; a "pair list" is a flat list of
sextuples; a "solve list" is a flat list of quadruples ``[a, b, c, n]``.

Every function returns an already-formatted newline-separated string so that
the two runners are compared byte-for-byte without relying on the generic
``format_result``.

Oracle caveats (see the unit report)
------------------------------------
The installed oracle is SageMath 10.3 / PARI 2.15.4 while the port targets the
vendored ``reference/sage`` (10.9.beta4) / ``reference/pari`` (2.18.1).  Three
places where those two upstreams disagree with *each other* are deliberately
kept out of the case list, because a difference there says nothing about the
port:

1. ``is_reduced`` on a **singular** form (``D == 0``).  10.3 evaluates a
   floating-point criterion and returns ``False``; 10.9 raises
   ``ValueError('the quadratic form must be non-singular')``
   (``reference/sage/src/sage/quadratic_forms/binary_qf.py:1489``, Sage issue
   #37635).  ``D == 0`` is therefore exercised through every *other* method.
2. Plain ``reduced_form()`` of an **indefinite non-square** form.  PARI 2.16
   "changed qfbred to use standard normalization (same as qfbredsl2)"
   (``reference/pari/CHANGES-2.16:142``), so 2.15.4 and 2.18.1 return different
   representatives *of the same class*.  Those forms are compared through the
   class invariant ``sorted(reduced_form().cycle(proper=True))`` instead, and
   exactly through ``reduced_form(transformation=True)`` (``qfbredsl2``, which
   is unchanged).
3. Which representative ``qfbsolve`` returns.  2.15.4 and 2.18.1 both return a
   correct solution but not the same one (compare the doctest at
   ``reference/.../binary_qf.py:1650`` with 10.3's).  Those cases therefore
   compare the *complete* brute-forced solution set plus the validity and
   existence of PARI's answer, which is version independent.
"""

from sage.all import *


# --------------------------------------------------------------------------
# formatting helpers (mirrored verbatim in the TypeScript area module)
# --------------------------------------------------------------------------

def _F(f):
    """``(a,b,c)`` for a BinaryQF or a 3-tuple."""
    return '(%s,%s,%s)' % (f[0], f[1], f[2])


def _T(x):
    return 'True' if x else 'False'


def _try(fn):
    """Run ``fn`` returning its string, or ``!<message>`` if it raised.

    Error messages are compared byte-for-byte; ``compare.ts`` would otherwise
    score "both sides raised" as a pass no matter how different the reasons.
    """
    try:
        return fn()
    except Exception as e:
        return '!' + str(e)


def _chunks(flat, k):
    return [tuple(ZZ(x) for x in flat[i:i + k]) for i in range(0, len(flat), k)]


def _forms(flat):
    return [BinaryQF(a, b, c) for (a, b, c) in _chunks(flat, 3)]


def _lines(items):
    return '\n'.join(items)


def _is_square(n):
    return n >= 0 and ZZ(n).is_square()


def _pari_qfb_ok(f):
    """Whether PARI's ``Qfb`` constructor accepts this form.

    ``reference/pari/src/basemath/Qfb.c:174-176``: a negative definite form
    raises ``pari_err_IMPL`` and a square (hence also zero) discriminant raises
    ``pari_err_DOMAIN``.  Sage reaches that constructor through
    ``_pari_init_``/``__pari__`` on both operands of ``__mul__``.
    """
    D = f.discriminant()
    if D < 0:
        return f[0] > 0
    return not _is_square(D)


def _canonical(f):
    """A version-independent name for the *class* of ``f``.

    For an indefinite non-square form the reduced representative returned by
    PARI depends on the PARI version (CHANGES-2.16 item 45), but the proper
    cycle it lives in does not, so we sort the cycle.

    ``D == 0`` short-circuits to ``singular``: reduction of a singular form
    raises on both sides but with a message that changed between SageMath 10.3
    and the vendored 10.9 (see the module docstring, caveat 1).
    """
    D = f.discriminant()
    if D == 0:
        return 'singular'
    if D > 0 and not _is_square(D):
        return ','.join(sorted(_F(g) for g in f.reduced_form().cycle(proper=True)))
    return _F(f.reduced_form())


def _sort_forms(forms):
    return sorted(forms, key=lambda f: (f[0], f[1], f[2]))


# --------------------------------------------------------------------------
# predicates / accessors
# --------------------------------------------------------------------------

def qf_predicates(flat):
    """``discriminant``, ``content`` and every ``is_*`` predicate except
    ``is_reduced`` (which is version-dependent for ``D == 0``; see module doc).
    """
    out = []
    for f in _forms(flat):
        flags = [
            _T(f.is_primitive()), _T(f.is_zero()), _T(f.is_reducible()),
            _T(f.is_positive_definite()), _T(f.is_negative_definite()),
            _T(f.is_indefinite()), _T(f.is_singular()), _T(f.is_nonsingular()),
        ]
        out.append('%s D=%s content=%s %s'
                   % (_F(f), f.discriminant(), f.content(), ' '.join(flags)))
    return _lines(out)


def qf_is_reduced(flat):
    """``is_reduced`` on non-singular forms."""
    out = []
    for f in _forms(flat):
        out.append('%s %s' % (_F(f), _try(lambda f=f: _T(f.is_reduced()))))
    return _lines(out)


def qf_evaluate(flat):
    """``Q(x, y)`` for a flat list of ``[a, b, c, x, y]`` quintuples."""
    out = []
    for (a, b, c, x, y) in _chunks(flat, 5):
        f = BinaryQF(a, b, c)
        out.append('%s(%s,%s)=%s' % (_F(f), x, y, f(x, y)))
    return _lines(out)


def qf_matrix_action(flat):
    """Right and left action of a 2x2 matrix on a form.

    Flat list of ``[a, b, c, m00, m01, m10, m11]`` septuples.
    """
    out = []
    for (a, b, c, p, q, r, s) in _chunks(flat, 7):
        f = BinaryQF(a, b, c)
        M = Matrix(ZZ, [[p, q], [r, s]])
        out.append('%s [%s,%s;%s,%s] right=%s left=%s'
                   % (_F(f), p, q, r, s,
                      _try(lambda: _F(f.matrix_action_right(M))),
                      _try(lambda: _F(f.matrix_action_left(M)))))
    return _lines(out)


def qf_principal(Ds):
    """``BinaryQF.principal(D)`` including the ``D % 4 not in (0, 1)`` error."""
    out = []
    for D in Ds:
        out.append('%s %s' % (D, _try(lambda D=D: _F(BinaryQF.principal(D)))))
    return _lines(out)


# --------------------------------------------------------------------------
# reduction
# --------------------------------------------------------------------------

def qf_reduced_form(flat):
    """``reduced_form()`` under all three ``algorithm`` values.

    Restricted by the case list to definite and square-discriminant forms; the
    indefinite non-square case goes through :func:`qf_reduced_class` and
    :func:`qf_reduced_transformation` (see module doc, caveat 2).
    """
    out = []
    for f in _forms(flat):
        for alg in ('default', 'pari', 'sage'):
            out.append('%s %s -> %s'
                       % (_F(f), alg,
                          _try(lambda f=f, alg=alg: _F(f.reduced_form(algorithm=alg)))))
    return _lines(out)


def qf_reduced_transformation(flat):
    """``reduced_form(transformation=True)`` -- form, base change and a check
    that ``f.matrix_action_right(M)`` really is the reduced form."""
    out = []
    for f in _forms(flat):
        for alg in ('default', 'pari', 'sage'):
            def go(f=f, alg=alg):
                g, M = f.reduced_form(transformation=True, algorithm=alg)
                ok = _T(f.matrix_action_right(M) == g)
                det = M.determinant()
                return '%s [%s,%s;%s,%s] det=%s action=%s' % (
                    _F(g), M[0][0], M[0][1], M[1][0], M[1][1], det, ok)
            out.append('%s %s -> %s' % (_F(f), alg, _try(go)))
    return _lines(out)


def qf_reduced_transformation_default(flat):
    """``reduced_form(transformation=True)`` under the default algorithm only.

    Used for large indefinite discriminants: ``algorithm='sage'`` runs Sage's
    Python ``_reduce_indef`` loop, whose termination test is the 53-bit
    floating-point ``is_reduced`` in SageMath 10.3 (fixed in the vendored 10.9,
    Sage issue #37635), so it does not terminate there.  The default algorithm
    is PARI's ``qfbredsl2``, which is exact in both.
    """
    out = []
    for f in _forms(flat):
        def go(f=f):
            g, M = f.reduced_form(transformation=True)
            ok = _T(f.matrix_action_right(M) == g)
            return '%s [%s,%s;%s,%s] det=%s action=%s' % (
                _F(g), M[0][0], M[0][1], M[1][0], M[1][1], M.determinant(), ok)
        out.append('%s default -> %s' % (_F(f), _try(go)))
    return _lines(out)


def qf_reduced_class(flat):
    """Class invariant of ``reduced_form()``: the sorted proper cycle.

    Version-independent even for indefinite non-square discriminants.
    """
    out = []
    for f in _forms(flat):
        out.append('%s -> %s' % (_F(f), _try(lambda f=f: _canonical(f))))
    return _lines(out)


def qf_cycle(flat):
    """``cycle()`` and ``cycle(proper=True)``."""
    out = []
    for f in _forms(flat):
        for proper in (False, True):
            def go(f=f, proper=proper):
                return ','.join(_F(g) for g in f.cycle(proper=proper))
            out.append('%s proper=%s -> %s' % (_F(f), _T(proper), _try(go)))
    return _lines(out)


# --------------------------------------------------------------------------
# composition
# --------------------------------------------------------------------------

def qf_compose(flat):
    """Raw Gauss composition ``f * g`` (PARI ``qfbcompraw``), not reduced."""
    out = []
    for (a1, b1, c1, a2, b2, c2) in _chunks(flat, 6):
        f = BinaryQF(a1, b1, c1)
        g = BinaryQF(a2, b2, c2)
        out.append('%s * %s = %s' % (_F(f), _F(g), _try(lambda: _F(f * g))))
    return _lines(out)


def qf_compose_domain(flat):
    """Composition on inputs PARI's ``Qfb`` constructor rejects.

    Sage builds both operands with ``Qfb(a, b, c)`` (``_pari_init_``,
    ``binary_qf.py:158-178``) before calling ``qfbcompraw``, so a negative
    definite operand or a square/zero discriminant raises a ``PariError``
    rather than returning a form
    (``reference/pari/src/basemath/Qfb.c:174-176``).
    """
    return qf_compose(flat)


def qf_compose_table(D):
    """Full composition table of the reduced representatives of ``D``.

    ``primitive_only=False``, so imprimitive classes (which only exist for
    non-fundamental ``D``) are included -- those are the ``gcd(a, b) > 1``
    inputs on which the port's hand-written Gauss composition used to be wrong.
    Each entry lists the raw composite and the class it lands in.
    """
    D = ZZ(D)
    R = BinaryQF_reduced_representatives(D, primitive_only=False, proper=True)
    out = ['D=%s n=%s reps=%s' % (D, len(R), ','.join(_F(f) for f in R)),
           'contents=' + ','.join(str(f.content()) for f in R)]
    for f in R:
        for g in R:
            def go(f=f, g=g):
                h = f * g
                return '%s ~ %s' % (_F(h), _canonical(h))
            out.append('%s * %s = %s' % (_F(f), _F(g), _try(go)))
    return _lines(out)


def qf_compose_powers(flat):
    """Iterated composition ``f, f*f, f*f*f, ...`` reduced at each step.

    ``[a, b, c, k]`` quadruples: the order of the class in the form class group
    shows up as the first ``i`` with ``f^i`` principal.
    """
    out = []
    for (a, b, c, k) in _chunks(flat, 4):
        f = BinaryQF(a, b, c)

        def go(f=f, k=k):
            acc = f
            parts = []
            for i in range(1, int(k) + 1):
                if i > 1:
                    acc = acc * f
                parts.append('%s:%s' % (i, _canonical(acc)))
            return ' '.join(parts)
        out.append('%s -> %s' % (_F(f), _try(go)))
    return _lines(out)


# --------------------------------------------------------------------------
# class groups
# --------------------------------------------------------------------------

def qf_reduced_representatives(Ds):
    """``BinaryQF_reduced_representatives`` under all four flag combinations."""
    out = []
    for D in Ds:
        for primitive_only in (False, True):
            for proper in (True, False):
                def go(D=D, po=primitive_only, pr=proper):
                    R = BinaryQF_reduced_representatives(D, primitive_only=po,
                                                         proper=pr)
                    return '%s %s' % (len(R), ','.join(_F(f) for f in R))
                out.append('D=%s primitive_only=%s proper=%s -> %s'
                           % (D, _T(primitive_only), _T(proper), _try(go)))
    return _lines(out)


def qf_class_number(Ds):
    """``class_number(D)``.

    The port exports this as ``quadratic_forms/binary_qf.ts:class_number``;
    upstream has no such helper, so the oracle is its definition:
    ``len(BinaryQF_reduced_representatives(D, primitive_only=True,
    proper=True))``.
    """
    out = []
    for D in Ds:
        def go(D=D):
            return str(len(BinaryQF_reduced_representatives(
                D, primitive_only=True, proper=True)))
        out.append('%s %s' % (D, _try(go)))
    return _lines(out)


def qf_is_equivalent(flat):
    """``is_equivalent`` with ``proper=True`` and ``proper=False``."""
    out = []
    for (a1, b1, c1, a2, b2, c2) in _chunks(flat, 6):
        f = BinaryQF(a1, b1, c1)
        g = BinaryQF(a2, b2, c2)
        for proper in (True, False):
            out.append('%s ~ %s proper=%s -> %s'
                       % (_F(f), _F(g), _T(proper),
                          _try(lambda f=f, g=g, p=proper: _T(f.is_equivalent(g, proper=p)))))
    return _lines(out)


def qf_class_group_closure(D):
    """Check that reduction really is a class invariant for ``D``.

    Every reduced representative is hit by exactly one class, every rep is
    reduced, and composing with the principal form is the identity.
    """
    D = ZZ(D)
    R = BinaryQF_reduced_representatives(D, primitive_only=False, proper=True)
    P = BinaryQF.principal(D)
    out = ['D=%s principal=%s reps=%s' % (D, _F(P), len(R))]
    for f in R:
        def go(f=f, P=P):
            return 'reduced=%s canon=%s id=%s inv=%s' % (
                _T(f.is_reduced()), _canonical(f), _canonical(f * P),
                _canonical(f * BinaryQF(f[0], -f[1], f[2])))
        out.append('%s %s' % (_F(f), _try(go)))
    return _lines(out)


# --------------------------------------------------------------------------
# solve_integer
# --------------------------------------------------------------------------

def _brute_force_solutions(a, b, c, n):
    """All ``(x, y)`` with ``a x^2 + b x y + c y^2 == n``; ``a > 0``, ``D < 0``.

    ``4 a n = (2 a x + b y)^2 - D y^2`` bounds ``|y|``, so the search is finite
    and exhaustive.  Pure integer arithmetic, identical on both sides.
    """
    D = b * b - 4 * a * c
    sols = []
    if n < 0:
        return sols
    ymax = ZZ(4 * a * n // (-D)).isqrt()
    for y in range(-ymax, ymax + 1):
        disc = D * y * y + 4 * a * n
        if disc < 0:
            continue
        s = ZZ(disc).isqrt()
        if s * s != disc:
            continue
        roots = [s] if s == 0 else [s, -s]
        for r in roots:
            num = -b * y + r
            if num % (2 * a) == 0:
                sols.append((num // (2 * a), ZZ(y)))
    return sorted(set(sols))


def qf_solve_integer_definite(flat):
    """``solve_integer`` on positive definite forms, compared exhaustively.

    Which of several solutions PARI returns is version dependent, so we compare
    the complete brute-forced solution set (identical integer code on both
    sides) plus whether PARI's answer exists, is valid and is one of them.
    """
    out = []
    for (a, b, c, n) in _chunks(flat, 4):
        f = BinaryQF(a, b, c)

        def go(f=f, a=a, b=b, c=c, n=n):
            sols = _brute_force_solutions(a, b, c, n)
            xy = f.solve_integer(n)
            found = xy is not None
            valid = found and f(xy[0], xy[1]) == n
            member = found and (ZZ(xy[0]), ZZ(xy[1])) in sols
            agree = (found == (len(sols) > 0)) and (not found or member)
            first = '(%s,%s)' % sols[0] if sols else 'None'
            return 'nsols=%s min=%s found=%s valid=%s agree=%s' % (
                len(sols), first, _T(found), _T(valid), _T(agree))
        out.append('%s n=%s -> %s' % (_F(f), n, _try(go)))
    return _lines(out)


def qf_solve_integer_existence(flat):
    """``solve_integer`` on indefinite non-square forms: existence + validity.

    Indefinite forms have infinitely many representations, so only the
    version-independent facts are compared.
    """
    out = []
    for (a, b, c, n) in _chunks(flat, 4):
        f = BinaryQF(a, b, c)

        def go(f=f, n=n):
            xy = f.solve_integer(n)
            if xy is None:
                return 'found=False'
            return 'found=True valid=%s' % _T(f(xy[0], xy[1]) == n)
        out.append('%s n=%s -> %s' % (_F(f), n, _try(go)))
    return _lines(out)


def qf_solve_integer_square_disc(flat):
    """``solve_integer`` on square-discriminant forms.

    PARI has no ``t_QFB`` of square discriminant, so this is Sage's own
    elementary algorithm (``binary_qf.py:1751-1791``) -- fully deterministic,
    so the exact pair is compared.  Covers ``a = 0``, ``c = 0``, ``b = 0``,
    ``D = 0`` and the ``ZeroDivisionError`` on the zero form.
    """
    out = []
    for (a, b, c, n) in _chunks(flat, 4):
        f = BinaryQF(a, b, c)

        def go(f=f, n=n):
            xy = f.solve_integer(n)
            if xy is None:
                return 'None'
            return '(%s,%s) check=%s' % (xy[0], xy[1], _T(f(xy[0], xy[1]) == n))
        out.append('%s n=%s -> %s' % (_F(f), n, _try(go)))
    return _lines(out)


def qf_solve_integer_cornacchia(flat):
    """``solve_integer(n, algorithm='cornacchia')``.

    ``[c, n]`` pairs for the form ``x^2 + c y^2``, plus the ``a != 1`` /
    ``b != 0`` / ``c <= 0`` rejection.  ``n`` is a prime or four times a prime
    (Sage documents that this is *not* checked, and PARI's behaviour on
    composite ``n`` differs between versions, so composites are not pinned).
    """
    out = []
    for (a, b, c, n) in _chunks(flat, 4):
        f = BinaryQF(a, b, c)

        def go(f=f, n=n):
            xy = f.solve_integer(n, algorithm='cornacchia')
            if xy is None:
                return 'None'
            return '(%s,%s) check=%s' % (xy[0], xy[1], _T(f(xy[0], xy[1]) == n))
        out.append('%s n=%s -> %s' % (_F(f), n, _try(go)))
    return _lines(out)


def qf_solve_integer_domain(flat):
    """``solve_integer`` on inputs the port has historically mishandled.

    In particular ``n = 0`` with a non-square discriminant: PARI's ``ifactor``
    returns ``0^1`` for ``0`` (``reference/pari/src/basemath/ifactor1.c:4459``)
    and ``qfbsolve`` then answers "no solution", so Sage returns ``None``.
    """
    out = []
    for (a, b, c, n) in _chunks(flat, 4):
        f = BinaryQF(a, b, c)

        def go(f=f, n=n):
            xy = f.solve_integer(n)
            return 'None' if xy is None else '(%s,%s)' % (xy[0], xy[1])
        out.append('%s n=%s -> %s' % (_F(f), n, _try(go)))
    return _lines(out)


# --------------------------------------------------------------------------
# structured / random batteries
# --------------------------------------------------------------------------

def _forms_of_discriminant(D, bound):
    """Every ``(a, b, c)`` of discriminant ``D`` with ``|a| <= bound`` and
    ``|b| <= 2*bound`` (including ``a = 0`` when ``D`` is a square).

    Deliberately *structured*: it produces negative-``a``, ``a = 0``,
    ``|a| > |c|`` and imprimitive shapes that uniform random sampling of small
    triples almost never hits.
    """
    D = ZZ(D)
    out = []
    for a in range(-bound, bound + 1):
        for b in range(-2 * bound, 2 * bound + 1):
            num = b * b - D
            if a == 0:
                if num == 0:
                    for c in range(-bound, bound + 1):
                        out.append((ZZ(0), ZZ(b), ZZ(c)))
                continue
            den = 4 * a
            if num % den == 0:
                out.append((ZZ(a), ZZ(b), num // den))
    return out


def qf_sweep_discriminant(args):
    """Sweep every form of a given discriminant through the whole API.

    ``args = [D, bound]``.  This is the anti-"small dense random" case: it
    enumerates the degenerate shapes exhaustively instead of sampling them.
    """
    D, bound = ZZ(args[0]), int(args[1])
    out = []
    for (a, b, c) in _forms_of_discriminant(D, bound):
        f = BinaryQF(a, b, c)

        def trans(f=f):
            g, M = f.reduced_form(transformation=True)
            return '%s|%s,%s,%s,%s' % (_F(g), M[0][0], M[0][1], M[1][0], M[1][1])
        parts = [
            'content=%s' % f.content(),
            'reduced=%s' % _try(lambda f=f: _T(f.is_reduced())),
            'canon=%s' % _try(lambda f=f: _canonical(f)),
            'trans=%s' % _try(trans),
            'sq=%s' % (_try(lambda f=f: _F(f * f)) if _pari_qfb_ok(f) else 'domain'),
        ]
        out.append('%s %s' % (_F(f), ' '.join(parts)))
    return _lines(out)


def qf_random_battery(values):
    """Random forms built from the shared Mersenne-Twister stream.

    ``values`` is a flat list of random integers; consecutive triples become
    forms.  A deliberate mixture: raw random triples (arbitrary discriminant,
    frequently singular or square) plus their "same discriminant" partners.
    """
    out = []
    for (a, b, c) in _chunks(values, 3):
        f = BinaryQF(a, b, c)
        D = f.discriminant()
        ok = _pari_qfb_ok(f)
        parts = [
            'D=%s' % D,
            'content=%s' % f.content(),
            'primitive=%s' % _T(f.is_primitive()),
            'reducible=%s' % _T(f.is_reducible()),
            'posdef=%s' % _T(f.is_positive_definite()),
            'canon=%s' % _try(lambda f=f: _canonical(f)),
            'sq=%s' % (_try(lambda f=f: _F(f * f)) if ok else 'domain'),
            'sqcanon=%s' % (_try(lambda f=f: _canonical(f * f)) if ok else 'domain'),
            'equiv=%s' % (_try(lambda f=f: _T(f.is_equivalent(f)))
                          if D != 0 else 'singular'),
        ]
        out.append('%s %s' % (_F(f), ' '.join(parts)))
    return _lines(out)


def qf_random_class_group(values):
    """Random walks in the class group of a random negative discriminant.

    For each random value ``v`` the discriminant ``D = -(|v| + 3)`` is rounded
    down to ``0`` or ``1`` mod 4, its reduced representatives are enumerated,
    and a deterministic product of representatives is reduced.  Exercises
    composition on non-fundamental and imprimitive inputs without any
    randomness that the two runners could disagree on.
    """
    out = []
    for v in values:
        D = -(abs(ZZ(v)) + 3)
        while D % 4 not in (0, 1):
            D -= 1

        def go(D=D, v=v):
            R = BinaryQF_reduced_representatives(D, primitive_only=False,
                                                 proper=True)
            if not R:
                return 'empty'
            acc = R[0]
            names = []
            for i in range(len(R)):
                acc = acc * R[(int(abs(v)) + i) % len(R)]
                names.append(_canonical(acc))
            return 'h=%s %s' % (len(R), ' '.join(names))
        out.append('D=%s %s' % (D, _try(go)))
    return _lines(out)


FUNCTIONS = {
    'qf_predicates': qf_predicates,
    'qf_is_reduced': qf_is_reduced,
    'qf_evaluate': qf_evaluate,
    'qf_matrix_action': qf_matrix_action,
    'qf_principal': qf_principal,
    'qf_reduced_form': qf_reduced_form,
    'qf_reduced_transformation': qf_reduced_transformation,
    'qf_reduced_transformation_default': qf_reduced_transformation_default,
    'qf_reduced_class': qf_reduced_class,
    'qf_cycle': qf_cycle,
    'qf_compose': qf_compose,
    'qf_compose_domain': qf_compose_domain,
    'qf_compose_table': qf_compose_table,
    'qf_compose_powers': qf_compose_powers,
    'qf_reduced_representatives': qf_reduced_representatives,
    'qf_class_number': qf_class_number,
    'qf_is_equivalent': qf_is_equivalent,
    'qf_class_group_closure': qf_class_group_closure,
    'qf_solve_integer_definite': qf_solve_integer_definite,
    'qf_solve_integer_existence': qf_solve_integer_existence,
    'qf_solve_integer_square_disc': qf_solve_integer_square_disc,
    'qf_solve_integer_cornacchia': qf_solve_integer_cornacchia,
    'qf_solve_integer_domain': qf_solve_integer_domain,
    'qf_sweep_discriminant': qf_sweep_discriminant,
    'qf_random_battery': qf_random_battery,
    'qf_random_class_group': qf_random_class_group,
}
