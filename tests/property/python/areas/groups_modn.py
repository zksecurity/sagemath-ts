"""SageMath side of the ``groups_modn`` property-test area.

Differential oracle for the port's

* ``packages/sagemath-ts/src/groups/generic.ts``          (``sage.groups.generic``)
* ``packages/sagemath-ts/src/rings/finite_rings/integer_mod.ts``      (``sage.rings.finite_rings.integer_mod``)
* ``packages/sagemath-ts/src/rings/finite_rings/integer_mod_ring.ts`` (``sage.rings.finite_rings.integer_mod_ring``)
* ``packages/sagemath-ts/src/rings/finite_rings/roots_of_unity.ts``
* ``packages/sagemath-ts/src/rings/finite_rings/tower_field.ts``

Cases: tests/property/cases/groups_modn.cases.json
TypeScript counterpart: tests/property/typescript/areas/groups_modn.ts

Conventions
-----------
Every function here returns an **already formatted string**, so the runner's
generic ``format_result`` is a no-op and both sides are compared byte for byte.

Exceptions are *not* left to propagate: ``compare.ts`` scores "both sides
raised" as a pass no matter what was raised, which would make every error case
vacuous.  Instead every function funnels through :func:`_run`, which turns an
exception into the string ``"<ExceptionClass>: <message>"``.  A wrong exception
class or a drifted message is therefore a **failure**, not a silent pass.

Version note
------------
The port targets the vendored SageMath ``reference/sage`` (10.9.beta4) while the
installed oracle is SageMath 10.3.  Three behaviours changed in between; each is
transcribed from the vendored source at the call site and marked ``VENDORED``:

1. ``order_from_multiple``'s ``check`` failure: 10.3 uses a bare ``assert``,
   10.9 raises ``ValueError(f"The order of P(={P}) does not divide {M}")``
   (reference/sage/src/sage/groups/generic.py:1361).
2. ``order_from_bounds(P, None, ...)``: the auto-widening bounds loop only
   exists from 10.9 (reference/sage/src/sage/groups/generic.py:1471-1480).
   Its answer is by definition the exact order of ``P``, so the oracle used
   here is Sage's own ``multiplicative_order``.
3. ``IntegerMod.log``'s ``order``/``check`` parameters do not exist in 10.3
   (reference/sage/src/sage/rings/finite_rings/integer_mod.pyx:795-798).
"""

from sage.all import *
from sage.groups.generic import (
    bsgs,
    discrete_log,
    discrete_log_lambda,
    discrete_log_rho,
    has_order,
    multiple,
    multiples,
    order_from_multiple,
    order_from_bounds,
)

Z = Integer


# ---------------------------------------------------------------- formatting

def _run(f):
    """Call ``f`` and format its result, or the exception it raised."""
    try:
        return f()
    except Exception as e:  # noqa: BLE001 - the class name is part of the result
        return '%s: %s' % (type(e).__name__, e)


def _fmt_list(xs):
    return '[' + ', '.join(str(x) for x in xs) + ']'


def _fmt_tuple(xs):
    return '(' + ', '.join(str(x) for x in xs) + ')'


# ------------------------------------------------------- vendored transcripts

def _order_from_multiple(P, m, operation):
    """``order_from_multiple`` with the vendored 10.9 ``check`` behaviour.

    VENDORED: reference/sage/src/sage/groups/generic.py:1360-1362 ::

        if check and _multiple(P, M) != identity:
            raise ValueError(f"The order of P(={P}) does not divide {M}")

    Sage 10.3 spells the same test ``assert multiple(...) == identity``, which
    raises a message-less ``AssertionError``; everything after the check is
    identical in the two versions.
    """
    identity = P.parent().one() if operation == '*' else P.parent().zero()
    if multiple(P, m, operation=operation) != identity:
        raise ValueError('The order of P(=%s) does not divide %s' % (P, m))
    return order_from_multiple(P, m, operation=operation, check=False)


# ------------------------------------------------------------ groups.generic

def gg_bsgs_mul(p, a, b, lb, ub):
    """bsgs in (Z/pZ)*: find lb <= n <= ub with a^n = b."""
    return _run(lambda: bsgs(Mod(a, p), Mod(b, p), (Z(lb), Z(ub)), operation='*'))


def gg_bsgs_add(n, a, b, lb, ub):
    """bsgs in (Z/nZ, +): find lb <= k <= ub with k*a = b."""
    return _run(lambda: bsgs(Mod(a, n), Mod(b, n), (Z(lb), Z(ub)), operation='+'))


def gg_multiple_mul(p, a, k):
    return _run(lambda: multiple(Mod(a, p), Z(k), operation='*'))


def gg_multiple_add(n, a, k):
    return _run(lambda: multiple(Mod(a, n), Z(k), operation='+'))


def gg_dlog_mul(p, a, b, ord_):
    """discrete_log(a, b, ord) in (Z/pZ)*; ``ord`` may be a proper multiple."""
    return _run(lambda: discrete_log(Mod(a, p), Mod(b, p), Z(ord_), operation='*'))


def gg_dlog_mul_noord(p, a, b):
    return _run(lambda: discrete_log(Mod(a, p), Mod(b, p), operation='*'))


def gg_dlog_add(n, a, b, ord_):
    return _run(lambda: discrete_log(Mod(a, n), Mod(b, n), Z(ord_), operation='+'))


def gg_order_from_multiple_mul(n, a, m):
    return _run(lambda: _order_from_multiple(Mod(a, n), Z(m), '*'))


def gg_order_from_multiple_default(n, a, m):
    """No ``operation`` argument: exercises the default (``'+'`` in Sage)."""
    return _run(lambda: _order_from_multiple(Mod(a, n), Z(m), '+'))


def gg_order_from_bounds_mul(n, a, lb, ub):
    return _run(lambda: order_from_bounds(Mod(a, n), (Z(lb), Z(ub)), operation='*'))


def gg_order_from_bounds_d(n, a, lb, ub, d):
    return _run(lambda: order_from_bounds(Mod(a, n), (Z(lb), Z(ub)), Z(d), operation='*'))


def gg_order_from_bounds_add(n, a, lb, ub):
    return _run(lambda: order_from_bounds(Mod(a, n), (Z(lb), Z(ub)), operation='+'))


def gg_order_from_bounds_nobounds(n, a):
    """``order_from_bounds(P, None, operation='*')``.

    VENDORED: reference/sage/src/sage/groups/generic.py:1471-1480 widens
    ``(1, 256)`` by 16x until ``bsgs`` succeeds and then returns
    ``order_from_multiple``, i.e. the exact multiplicative order of ``P``.
    Sage 10.3 has no such branch, so the oracle is the order itself.
    """
    return _run(lambda: Mod(a, n).multiplicative_order())


def gg_multiple_of_order(n, a):
    """The port's ``multiple_of_order`` (no Sage counterpart) returns the exact
    order, being ``order_from_bounds`` with no bounds -- oracle as above."""
    return _run(lambda: Mod(a, n).multiplicative_order())


def gg_has_order_mul(n, a, m):
    return _run(lambda: has_order(Mod(a, n), Z(m), '*'))


def gg_has_order_default(n, a, m):
    """No ``operation`` argument: exercises the default (``'+'`` in Sage)."""
    return _run(lambda: has_order(Mod(a, n), Z(m)))


def gg_has_order_other(n, a, m):
    """``operation='other'``: Sage raises ``ValueError('unknown group operation')``."""
    return _run(lambda: has_order(Mod(a, n), Z(m), 'other'))


def gg_multiples_default(n, a, k):
    """No ``operation``/``indexed``: Sage yields ``0, a, 2a, ...`` (additive)."""
    return _run(lambda: _fmt_list(list(multiples(Mod(a, n), Z(k)))))


def gg_multiples_mul(p, a, k):
    return _run(lambda: _fmt_list(list(multiples(Mod(a, p), Z(k), operation='*'))))


def gg_multiples_indexed(n, a, k):
    return _run(lambda: _fmt_list([_fmt_tuple(t) for t in multiples(Mod(a, n), Z(k), indexed=True)]))


def gg_dlog_lambda(p, base, x, lb, ub):
    """Pollard lambda for the known logarithm ``x`` of ``base^x``."""
    return _run(lambda: discrete_log_lambda(
        Mod(base, p) ** Z(x), Mod(base, p), (Z(lb), Z(ub)), operation='*'))


def gg_dlog_rho(p, base, x, ord_):
    """Pollard rho; Sage requires ``ord`` to be prime."""
    return _run(lambda: discrete_log_rho(
        Mod(base, p) ** Z(x), Mod(base, p), Z(ord_), operation='*'))


# --------------------------------------------------------------- integer_mod

def im_log(n, a, b):
    return _run(lambda: Mod(a, n).log(Mod(b, n)))


def im_log_nobase(n, a):
    """``log()`` with no base: Sage uses ``parent.multiplicative_generator()``."""
    return _run(lambda: Mod(a, n).log())


def im_log_check(n, a, b, order):
    """``log(b, order, check=True)``.

    VENDORED: reference/sage/src/sage/rings/finite_rings/integer_mod.pyx:795-798 ::

        if check:
            if not has_order(b, order, '*'):
                raise ValueError('base does not have the provided order')

    Sage 10.3's ``log`` accepts only the base argument, so the guard is
    transcribed here and the rest delegated to Sage.
    """
    def go():
        x = Mod(a, n)
        base = Mod(b, n)
        if not x.is_unit():
            raise ValueError('logarithm of %s is not defined since it is not a unit modulo %s' % (x, n))
        if not base.is_unit():
            raise ValueError('logarithm with base %s is not defined since it is not a unit modulo %s' % (base, n))
        if not has_order(base, Z(order), '*'):
            raise ValueError('base does not have the provided order')
        return x.log(base)
    return _run(go)


def im_mult_order(n, a):
    return _run(lambda: Mod(a, n).multiplicative_order())


def im_inv(n, a):
    return _run(lambda: Mod(a, n) ** (-1))


def im_div(n, a, b):
    return _run(lambda: Mod(a, n) / Mod(b, n))


def im_pow(n, a, e):
    return _run(lambda: Mod(a, n) ** Z(e))


def im_add(n, a, b):
    return _run(lambda: Mod(a, n) + Mod(b, n))


def im_sub(n, a, b):
    return _run(lambda: Mod(a, n) - Mod(b, n))


def im_mul(n, a, b):
    return _run(lambda: Mod(a, n) * Mod(b, n))


def im_is_unit(n, a):
    return _run(lambda: Mod(a, n).is_unit())


def im_lift(n, a):
    return _run(lambda: Mod(a, n).lift())


# ---------------------------------------------------------- integer_mod_ring

def imr_unit_gens(n):
    return _run(lambda: _fmt_list(list(Integers(n).unit_gens())))


def imr_mult_gen(n):
    return _run(lambda: Integers(n).multiplicative_generator())


def imr_is_cyclic(n):
    return _run(lambda: Integers(n).multiplicative_group_is_cyclic())


def imr_is_field(n):
    return _run(lambda: Integers(n).is_field())


def imr_cardinality(n):
    return _run(lambda: Integers(n).cardinality())


def imr_units(n):
    """All units, ascending.  Sage: ``list_of_elements_of_multiplicative_group``."""
    return _run(lambda: _fmt_list(sorted(Z(u) for u in Integers(n).list_of_elements_of_multiplicative_group())))


def imr_list(n):
    return _run(lambda: _fmt_list(list(Integers(n))))


# ------------------------------------------------------------ roots_of_unity

def ru_has_primitive_root(p, n):
    """``n``-th roots of unity exist in GF(p) iff n | p-1 (n = 1 always)."""
    def go():
        if Z(n) <= 0:
            return False
        if Z(n) == 1:
            return True
        return (GF(p).order() - 1) % Z(n) == 0
    return _run(go)


def ru_primitive_nth_root(p, n):
    """``g^((p-1)/n)`` with ``g`` Sage's multiplicative generator of GF(p).

    This is the port's documented construction (roots_of_unity.ts:132), so the
    comparison also pins the choice of generator.
    """
    def go():
        F = GF(p)
        if Z(n) <= 0:
            raise ValueError('n must be positive')
        if Z(n) == 1:
            return F.one()
        if (F.order() - 1) % Z(n) != 0:
            raise ValueError(
                'No primitive %s-th root of unity exists in the field. '
                '%s does not divide %s = |F*|.' % (n, n, F.order() - 1))
        return F.multiplicative_generator() ** ((F.order() - 1) // Z(n))
    return _run(go)


def ru_roots_of_unity(p, n):
    """``[1, w, ..., w^(n-1)]`` for the primitive root of :func:`ru_primitive_nth_root`."""
    def go():
        F = GF(p)
        if Z(n) <= 0:
            raise ValueError('n must be positive')
        if Z(n) == 1:
            return _fmt_list([F.one()])
        if (F.order() - 1) % Z(n) != 0:
            raise ValueError(
                'No primitive %s-th root of unity exists in the field. '
                '%s does not divide %s = |F*|.' % (n, n, F.order() - 1))
        w = F.multiplicative_generator() ** ((F.order() - 1) // Z(n))
        out = []
        cur = F.one()
        for _ in range(Z(n)):
            out.append(cur)
            cur = cur * w
        return _fmt_list(out)
    return _run(go)


def ru_mult_order(p, a):
    def go():
        x = GF(p)(a)
        if x == 0:
            raise ArithmeticError('Multiplicative order of 0 not defined.')
        return x.multiplicative_order()
    return _run(go)


def ru_find_mult_gen(p):
    return _run(lambda: GF(p).multiplicative_generator())


def ru_element_of_order(p, n):
    def go():
        F = GF(p)
        if Z(n) <= 0:
            raise ValueError('order must be positive')
        if (F.order() - 1) % Z(n) != 0:
            raise ValueError('No element of order %s exists. %s does not divide %s.'
                             % (n, n, F.order() - 1))
        return F.multiplicative_generator() ** ((F.order() - 1) // Z(n))
    return _run(go)


def ru_cyclotomic(n):
    return _run(lambda: _fmt_list(cyclotomic_polynomial(Z(n)).list()))


def ru_max_fft_size(p):
    return _run(lambda: Z(2) ** ((GF(p).order() - 1).valuation(2)))


def ru_two_adicity(p):
    return _run(lambda: (GF(p).order() - 1).valuation(2))


def ru_valid_fft_sizes(p):
    return _run(lambda: _fmt_list(
        [Z(2) ** i for i in range((GF(p).order() - 1).valuation(2) + 1)]))


def ru_fft_domain(p, n):
    """The FFT domain ``[1, w, w^2, ...]`` of size ``n`` (a power of two)."""
    def go():
        F = GF(p)
        nn = Z(n)
        if nn <= 0:
            raise ValueError('domain size must be positive')
        if nn & (nn - 1) != 0:
            raise ValueError('domain size %s must be a power of 2' % nn)
        if (F.order() - 1) % nn != 0:
            raise ValueError('No %s-th roots of unity in the field. %s does not divide %s.'
                             % (nn, nn, F.order() - 1))
        w = F.multiplicative_generator() ** ((F.order() - 1) // nn)
        out = []
        cur = F.one()
        for _ in range(nn):
            out.append(cur)
            cur *= w
        return _fmt_list(out)
    return _run(go)


def ru_coset(p, n, offset):
    """``[c, c*w, c*w^2, ...]`` -- the FFT domain shifted by ``offset``."""
    def go():
        F = GF(p)
        nn = Z(n)
        w = F.multiplicative_generator() ** ((F.order() - 1) // nn)
        out = []
        cur = F(offset)
        for _ in range(nn):
            out.append(cur)
            cur *= w
        return _fmt_list(out)
    return _run(go)


def ru_coset_fold(p, n, offset):
    """One FRI fold: size n/2, generator w^2, offset c^2."""
    def go():
        F = GF(p)
        nn = Z(n)
        w = F.multiplicative_generator() ** ((F.order() - 1) // nn)
        w2 = w * w
        out = []
        cur = F(offset) ** 2
        for _ in range(nn // 2):
            out.append(cur)
            cur *= w2
        return _fmt_list(out)
    return _run(go)


# --------------------------------------------------------------- tower_field

def _tower(i):
    """The binary tower ``T_0 = GF(2)``, ``T_{j+1} = T_j[X]/(X^2 + a_j X + 1)``
    with ``a_j = T_j.gen()`` -- the construction of tower_field.ts:79-101."""
    T = GF(2)
    for j in range(int(i)):
        R = PolynomialRing(T, 'x%d' % j)
        x = R.gen()
        T = R.quotient(x ** 2 + T.gen() * x + 1)
    return T


def _tower_decode(T, i, v):
    """Integer -> element of ``T_i``, little-endian in the tower basis."""
    i = int(i)
    v = Z(v)
    if i == 0:
        return T(v % 2)
    half = Z(2) ** (Z(2) ** (i - 1))
    S = T.base_ring()
    lo = _tower_decode(S, i - 1, v % half)
    hi = _tower_decode(S, i - 1, v // half)
    return T(lo + hi * T.gen())


def _tower_encode(i, e):
    """Element of ``T_i`` -> integer (inverse of :func:`_tower_decode`)."""
    i = int(i)
    if i == 0:
        return Z(e)
    half = Z(2) ** (Z(2) ** (i - 1))
    c = e.lift().list()
    c = c + [e.parent().base_ring().zero()] * (2 - len(c))
    return _tower_encode(i - 1, c[0]) + half * _tower_encode(i - 1, c[1])


def tf_cardinality(i):
    return _run(lambda: _tower(i).cardinality())


def tf_add(i, a, b):
    def go():
        T = _tower(i)
        return _tower_encode(i, _tower_decode(T, i, a) + _tower_decode(T, i, b))
    return _run(go)


def tf_mul(i, a, b):
    def go():
        T = _tower(i)
        return _tower_encode(i, _tower_decode(T, i, a) * _tower_decode(T, i, b))
    return _run(go)


def tf_inv(i, a):
    def go():
        T = _tower(i)
        x = _tower_decode(T, i, a)
        if x == 0:
            raise ZeroDivisionError('division by zero')
        return _tower_encode(i, x ** (-1))
    return _run(go)


def tf_pow(i, a, e):
    def go():
        T = _tower(i)
        return _tower_encode(i, _tower_decode(T, i, a) ** Z(e))
    return _run(go)


def tf_gen(i):
    return _run(lambda: _tower_encode(i, _tower(i).gen()))


def tf_elements(i):
    """All elements of ``T_i`` in the ring's own iteration order, encoded."""
    return _run(lambda: _fmt_list([_tower_encode(i, e) for e in _tower(i)]))


def tf_mult_order(i, a):
    """Multiplicative order of the encoded element (brute force: |T_i| <= 256)."""
    def go():
        T = _tower(i)
        x = _tower_decode(T, i, a)
        if x == 0:
            raise ArithmeticError('Multiplicative order of 0 not defined.')
        k = 1
        y = x
        while y != T.one():
            y *= x
            k += 1
        return Z(k)
    return _run(go)


FUNCTIONS = {
    # groups/generic.ts
    'gg_bsgs_mul': gg_bsgs_mul,
    'gg_bsgs_add': gg_bsgs_add,
    'gg_multiple_mul': gg_multiple_mul,
    'gg_multiple_add': gg_multiple_add,
    'gg_dlog_mul': gg_dlog_mul,
    'gg_dlog_mul_noord': gg_dlog_mul_noord,
    'gg_dlog_add': gg_dlog_add,
    'gg_order_from_multiple_mul': gg_order_from_multiple_mul,
    'gg_order_from_multiple_default': gg_order_from_multiple_default,
    'gg_order_from_bounds_mul': gg_order_from_bounds_mul,
    'gg_order_from_bounds_d': gg_order_from_bounds_d,
    'gg_order_from_bounds_add': gg_order_from_bounds_add,
    'gg_order_from_bounds_nobounds': gg_order_from_bounds_nobounds,
    'gg_multiple_of_order': gg_multiple_of_order,
    'gg_has_order_mul': gg_has_order_mul,
    'gg_has_order_default': gg_has_order_default,
    'gg_has_order_other': gg_has_order_other,
    'gg_multiples_default': gg_multiples_default,
    'gg_multiples_mul': gg_multiples_mul,
    'gg_multiples_indexed': gg_multiples_indexed,
    'gg_dlog_lambda': gg_dlog_lambda,
    'gg_dlog_rho': gg_dlog_rho,
    # integer_mod.ts
    'im_log': im_log,
    'im_log_nobase': im_log_nobase,
    'im_log_check': im_log_check,
    'im_mult_order': im_mult_order,
    'im_inv': im_inv,
    'im_div': im_div,
    'im_pow': im_pow,
    'im_add': im_add,
    'im_sub': im_sub,
    'im_mul': im_mul,
    'im_is_unit': im_is_unit,
    'im_lift': im_lift,
    # integer_mod_ring.ts
    'imr_unit_gens': imr_unit_gens,
    'imr_mult_gen': imr_mult_gen,
    'imr_is_cyclic': imr_is_cyclic,
    'imr_is_field': imr_is_field,
    'imr_cardinality': imr_cardinality,
    'imr_units': imr_units,
    'imr_list': imr_list,
    # roots_of_unity.ts
    'ru_has_primitive_root': ru_has_primitive_root,
    'ru_primitive_nth_root': ru_primitive_nth_root,
    'ru_roots_of_unity': ru_roots_of_unity,
    'ru_mult_order': ru_mult_order,
    'ru_find_mult_gen': ru_find_mult_gen,
    'ru_element_of_order': ru_element_of_order,
    'ru_cyclotomic': ru_cyclotomic,
    'ru_max_fft_size': ru_max_fft_size,
    'ru_two_adicity': ru_two_adicity,
    'ru_valid_fft_sizes': ru_valid_fft_sizes,
    'ru_fft_domain': ru_fft_domain,
    'ru_coset': ru_coset,
    'ru_coset_fold': ru_coset_fold,
    # tower_field.ts
    'tf_cardinality': tf_cardinality,
    'tf_add': tf_add,
    'tf_mul': tf_mul,
    'tf_inv': tf_inv,
    'tf_pow': tf_pow,
    'tf_gen': tf_gen,
    'tf_elements': tf_elements,
    'tf_mult_order': tf_mult_order,
}


def gg_lambda_trace(seed, mode, p, b, lb, ub, x, hash_mode):
    """Live Sage Pollard lambda, with explicit signed hash and full walk trace.

    Modes 0/1/2 are addition/multiplication/custom addition. Width zero remains
    the separately documented upstream-crash deviation and is excluded here.
    """
    set_random_seed(Z(seed))
    trace = []
    base = Mod(b, p)
    operation = ['+', '*', 'other'][int(mode)]
    def hs(v):
        trace.append(int(v))
        z = int(v)**2 + 1
        return [z, -z, -1, -(1 << 100) + z, 0][int(hash_mode)]
    kw = dict(operation=operation, hash_function=hs)
    if operation == 'other':
        kw.update(identity=Mod(0, p), inverse=lambda v: -v, op=lambda u, v: u+v)
    try:
        target = base**Z(x) if operation == '*' else Z(x)*base
        result = str(discrete_log_lambda(target, base, (Z(lb), Z(ub)), **kw))
        error = None
    except Exception as e:
        result = None
        error = type(e).__name__ + ': ' + str(e)
    import json
    return json.dumps([result, error, trace], separators=(',', ':'))

FUNCTIONS['gg_lambda_trace'] = gg_lambda_trace


_BUNDLED_GROUPS_GENERIC = None


def gg_rho_state(seed, mode, p, b, order, target):
    """Run the bundled original with the port's documented DJB2 string hash.

    Loading the original Python module avoids the installed 10.3 BSGS identity
    error, which changed in the bundled version. No algorithm is transcribed.
    The next 64 CPython bits detect draw order, generator and shortcut drift.
    """
    native = _bundled_group_source()
    from sage.misc.prandom import getrandbits
    import json
    set_random_seed(Z(seed))
    base = Mod(b, p)
    operation = ['+', '*', 'other'][int(mode)]
    def hs(v):
        h = 5381
        for c in str(v):
            h = ((h * 33) ^ ord(c)) & 0xffffffff
        return h
    kw = dict(operation=operation, hash_function=hs)
    if operation == 'other':
        kw.update(identity=Mod(0, p), inverse=lambda v: -v, op=lambda u, v: u+v)
    try:
        result = str(native.discrete_log_rho(Mod(target, p), base, Z(order), **kw))
        error = None
    except Exception as e:
        result = None
        error = type(e).__name__ + ': ' + str(e)
    return json.dumps([result, error, str(getrandbits(64))], separators=(',', ':'))

FUNCTIONS['gg_rho_state'] = gg_rho_state


def _bundled_group_source():
    global _BUNDLED_GROUPS_GENERIC
    if _BUNDLED_GROUPS_GENERIC is None:
        import importlib.util
        from pathlib import Path
        source = Path(__file__).resolve().parents[4] / 'reference/sage/src/sage/groups/generic.py'
        spec = importlib.util.spec_from_file_location('_audit_bundled_groups_generic', source)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        _BUNDLED_GROUPS_GENERIC = module
    return _BUNDLED_GROUPS_GENERIC


def gg_group_parse(fn, mode, mask, variant):
    import json
    g = _bundled_group_source()
    set_random_seed(0)
    fn, mask, v = int(fn), int(mask), int(variant)
    operation = ['+', '*', 'other', 'custom'][int(mode)]
    a, b = Mod(2, 11), Mod(4, 11)
    kw = dict(operation=operation, identity=Mod(0, 11) if mask & 1 else None,
              inverse=(lambda x: -x) if mask & 2 else None,
              op=(lambda x, y: x+y) if mask & 4 else None)
    bounds = [(Z(0), Z(8)), (Z(-1), Z(8)), (Z(8), Z(0))][v]
    try:
        if fn == 0: r = g.multiple(a, [-1, 0, 3][v], **kw)
        elif fn == 1: r = g.bsgs(a, b, bounds, **kw)
        elif fn == 2: r = g.discrete_log(b, a, [10 if operation == '*' else 11, 0, -1][v], **kw)
        elif fn == 3: r = g.order_from_multiple(a, [110, 3, 0][v], **kw)
        elif fn == 4: r = g.order_from_bounds(a, [(1, 22), (-1, 22), (22, 1)][v], **kw)
        elif fn == 5: r = g.discrete_log_lambda(b, a, bounds, hash_function=lambda x: int(x)**2, **kw)
        elif fn == 6: r = g.discrete_log_rho(b, a, [11, 0, -1][v], **kw)
        else:
            _, identity, inverse, op = g._parse_group_def(a.parent(), **kw)
            r = json.dumps([int(identity), int(inverse(a)), int(op(a, b))], separators=(',', ':'))
        result, error = str(r), None
    except Exception as e:
        result, error = None, type(e).__name__ + ': ' + str(e)
    return json.dumps([result, error], separators=(',', ':'))


def gg_bounds_trace(value, lb, ub, d):
    import json
    trace = []
    def op(x, y):
        trace.append([int(x), int(y)])
        return x+y
    try:
        r = _bundled_group_source().order_from_bounds(Mod(value, 11), (Z(lb), Z(ub)), Z(d),
            operation='other', identity=Mod(0, 11), inverse=lambda x: -x, op=op)
        result, error = str(r), None
    except Exception as e:
        result, error = None, type(e).__name__ + ': ' + str(e)
    return json.dumps([result, error, trace], separators=(',', ':'))


def gg_multiple_trace(p, value, mode, n):
    import json
    trace = []
    def op(x, y):
        trace.append([int(x), int(y)])
        return x+y if mode == 0 else x*y if mode == 1 else min(x, y)
    identity = Mod(1 if mode == 1 else p-1 if mode == 2 else 0, p)
    inverse = (lambda x: ~x) if mode == 1 else (lambda x: -x)
    try:
        r = _bundled_group_source().multiple(Mod(value, p), Z(n), operation='other',
            identity=identity, inverse=inverse, op=op)
        result, error = str(r), None
    except Exception as e:
        result, error = None, type(e).__name__ + ': ' + str(e)
    return json.dumps([result, error, trace], separators=(',', ':'))

FUNCTIONS.update(gg_group_parse=gg_group_parse, gg_bounds_trace=gg_bounds_trace, gg_multiple_trace=gg_multiple_trace)


def gg_order_list_trace(modulus, value, multiple_value, factor_mode, packed_factors, plist_mode, plist):
    """Integer-valued factor/prime lists; compare all custom operation calls.

    Prime entries use Sage Integer, the target of the IntegerLike adapter.
    Zero-multiple cases whose original loop does not terminate are excluded.
    """
    import json
    trace = []
    def op(x, y):
        trace.append([int(x), int(y)])
        return x+y
    factors = None if factor_mode == 0 else [(Z(packed_factors[i]), Z(packed_factors[i+1])) for i in range(0, len(packed_factors), 2)]
    primes = None if plist_mode == 0 else [Z(p) for p in plist]
    try:
        r = _bundled_group_source().order_from_multiple(Mod(value, modulus), Z(multiple_value),
            factorization=factors, plist=primes, operation='other', identity=Mod(0, modulus),
            inverse=lambda x: -x, op=op)
        result, error = str(r), None
    except Exception as e:
        result, error = None, type(e).__name__ + ': ' + str(e)
    return json.dumps([result, error, trace], separators=(',', ':'))

FUNCTIONS['gg_order_list_trace'] = gg_order_list_trace


from group_iterator import iterator_state, iterator_standard, parent_field


def gg_iterator_state(*args): return iterator_state(_bundled_group_source(), *args)
def gg_iterator_standard(*args): return iterator_standard(_bundled_group_source(), *args)
def gg_parent_field(*args): return parent_field(_bundled_group_source(), *args)

FUNCTIONS.update(gg_iterator_state=gg_iterator_state, gg_iterator_standard=gg_iterator_standard, gg_parent_field=gg_parent_field)

from group_number_field import nf_scalar, nf_group

def gg_nf_group(*args): return nf_group(_bundled_group_source(), *args)
FUNCTIONS.update(gg_nf_scalar=nf_scalar, gg_nf_group=gg_nf_group)

from integer_rational import integer_group
def gg_integer(*args): return integer_group(_bundled_group_source(), *args)
FUNCTIONS["gg_integer"] = gg_integer

from real_double_group import rdf_multiple
def gg_rdf_multiple(*args):return rdf_multiple(_bundled_group_source(),*args)
FUNCTIONS["gg_rdf_multiple"] = gg_rdf_multiple

from group_schedule import group_schedule, rdf_schedule
def gg_schedule(*args):return group_schedule(_bundled_group_source(),*args)
def gg_rdf_schedule(*args):return rdf_schedule(_bundled_group_source(),*args)
FUNCTIONS.update(gg_schedule=gg_schedule,gg_rdf_schedule=gg_rdf_schedule)

from group_parent import group_parent
def gg_parent(*args):return group_parent(_bundled_group_source(),*args)
FUNCTIONS['gg_parent']=gg_parent

from group_pollard import pollard,pollard_matrix
def gg_pollard(*args):return pollard(_bundled_group_source(),*args)
def gg_pollard_matrix(*args):return pollard_matrix(_bundled_group_source(),*args)
FUNCTIONS.update(gg_pollard=gg_pollard,gg_pollard_matrix=gg_pollard_matrix)


def gg_linear_relation(modulus, a, b, mode, hints):
    import json
    P,Q=Mod(a,modulus),Mod(b,modulus)
    mode,hints=int(mode),int(hints)
    kw=dict(operation=['+','*','other'][mode])
    if mode==2: kw.update(identity=Mod(0,modulus),inverse=lambda x:-x,op=lambda x,y:x+y)
    if hints & 1: kw['ord_p']=P.multiplicative_order() if mode==1 else P.additive_order()
    if hints & 2: kw['ord_q']=Q.multiplicative_order() if mode==1 else Q.additive_order()
    try: result={'value':list(map(str,_bundled_group_source().linear_relation(P,Q,**kw)))}
    except Exception as e: result={'error':type(e).__name__,'message':str(e)}
    return json.dumps(result,separators=(',',':'))
FUNCTIONS['gg_linear_relation']=gg_linear_relation

FUNCTIONS['im_additive_order'] = lambda modulus,a: str(Mod(a,modulus).additive_order())
