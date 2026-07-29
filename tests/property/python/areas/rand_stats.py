"""SageMath side of the ``rand_stats`` property-test area.

Cases: tests/property/cases/rand_stats.cases.json
TypeScript counterpart: tests/property/typescript/areas/rand_stats.ts

Scope
-----
* ``sage.misc.randstate``            -> ``src/misc/randstate.ts``
* ``sage.stats.distributions.*``     -> ``src/stats/distributions/*``
* the ``sage.rings.integer_ring`` entry points that draw from ``randstate``
  (``ZZ.random_element``), because those are what make a randstate defect
  observable in user code.

Everything here is a *seeded stream* comparison: both sides call
``set_random_seed(s)`` and then emit the exact sequence of draws.  That is the
only way to catch the class of bug this area exists for -- the port's RNG once
had period 2, and a distribution-shape test would not have noticed.

Oracle-version note (IMPORTANT)
-------------------------------
``reference/sage`` is SageMath 10.9.beta4; the ``sage`` binary on this machine
is 10.3.  For ``sage.misc.randstate`` and
``sage.stats.distributions.discrete_gaussian_integer`` the two agree (verified
by running both), so those functions use the *installed* Sage.

``sage.stats.distributions.discrete_gaussian_lattice`` was rewritten between
10.3 and 10.9 (non-spherical support, Poisson-summation
``_normalisation_factor_zz``, ``c`` became a method, different ``__repr__``,
different ``_c_in_lattice`` predicate).  The port targets 10.9, so using the
10.3 module here would report version skew as port defects.  Instead we import
the **vendored 10.9 source file itself** and run it under the installed Sage --
it is pure Python and loads cleanly.  ``_load_vendored_dgl`` below asserts that
it reproduces the vendored docstring values, so the oracle is self-checking.
"""

import importlib.util
import os
import warnings

from sage.all import *
from sage.misc.randstate import current_randstate, initial_seed
from sage.misc.randstate import random as sage_random
from sage.stats.distributions.discrete_gaussian_integer import (
    DiscreteGaussianDistributionIntegerSampler as DGI,
)

# ---------------------------------------------------------------------------
# Vendored (10.9) discrete_gaussian_lattice
# ---------------------------------------------------------------------------

_REPO_ROOT = os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..')
)
_VENDORED_DGL_PATH = os.path.join(
    _REPO_ROOT, 'reference', 'sage', 'src', 'sage', 'stats', 'distributions',
    'discrete_gaussian_lattice.py',
)

_DGL_MODULE = None


def _load_vendored_dgl():
    """Import ``reference/sage/.../discrete_gaussian_lattice.py`` as a module.

    The module is cached.  On first load we check three values straight out of
    its own docstrings, so that a silently-wrong import can never masquerade as
    an oracle.
    """
    global _DGL_MODULE
    if _DGL_MODULE is not None:
        return _DGL_MODULE
    spec = importlib.util.spec_from_file_location(
        'sagemath_ts_vendored_dgl', _VENDORED_DGL_PATH
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    DGL = module.DiscreteGaussianDistributionLatticeSampler
    # discrete_gaussian_lattice.py:239-244 and :211
    check = DGL(ZZ**8, 0.5)
    assert str(check._normalisation_factor_zz(tau=3))[:6] == '3.1653', \
        'vendored discrete_gaussian_lattice.py did not reproduce its own doctest'
    assert str(check._normalisation_factor_zz())[:6] == '6.8249', \
        'vendored discrete_gaussian_lattice.py did not reproduce its own doctest'
    assert str(DGL(ZZ**3, 1.0)._normalisation_factor_zz())[:7] == '15.7496', \
        'vendored discrete_gaussian_lattice.py did not reproduce its own doctest'

    _DGL_MODULE = module
    return module


def _DGL():
    return _load_vendored_dgl().DiscreteGaussianDistributionLatticeSampler


# ---------------------------------------------------------------------------
# Formatting helpers (must be byte-identical to the TypeScript versions)
# ---------------------------------------------------------------------------

def _ints(xs):
    """``[1, 2, 3]``."""
    return '[' + ', '.join(str(Integer(x)) for x in xs) + ']'


def _vecs(vs):
    """``[(1, 2), (3, 4)]``."""
    return '[' + ', '.join(
        '(' + ', '.join(str(Integer(x)) for x in v) + ')' for v in vs
    ) + ']'


def _f(x):
    """Canonical 12-significant-digit form of a double, e.g. ``1.11439293741e-1``.

    12 significant digits is below the ~15.95 digits an IEEE-754 double
    carries, so it absorbs last-ulp differences between MPFR's and libm's
    ``exp``/``sqrt`` while still pinning the value hard.
    """
    x = float(x)
    if x != x:
        return 'nan'
    if x in (float('inf'), float('-inf')):
        return 'inf' if x > 0 else '-inf'
    mantissa, exponent = ('%.11e' % x).split('e')
    return '%se%d' % (mantissa, int(exponent))


def _fs(xs):
    return '[' + ', '.join(_f(x) for x in xs) + ']'


def _rat(num, den):
    """A double built from two integers, identically on both sides."""
    return float(Integer(num)) / float(Integer(den))


_ALGORITHMS = {
    0: None,
    1: 'uniform+table',
    2: 'uniform+online',
    3: 'uniform+logtable',
    4: 'sigma2+logtable',
}


def _alg(code):
    code = int(code)
    if code not in _ALGORITHMS:
        raise ValueError('unknown algorithm code %d' % code)
    return _ALGORITHMS[code]


def _square(flat, n):
    n = int(n)
    return [[Integer(flat[i * n + j]) for j in range(n)] for i in range(n)]


# ===========================================================================
# sage.misc.randstate
# ===========================================================================

def rs_random(seed, count):
    """``set_random_seed(seed)`` then ``count`` draws of the 31-bit generator.

    ``sage.misc.randstate.random`` is ``randstate.c_random`` == GMP
    ``gmp_urandomb_ui(state, 31)`` (``randstate.pyx:881``).
    """
    set_random_seed(Integer(seed))
    return _ints([sage_random() for _ in range(int(count))])


def rs_random_modseed(offset, count):
    """Same, but seeded at ``2^19937 - 20027 + offset``.

    ``2^19937 - 20027`` is the modulus GMP's ``randseed_mt`` reduces the seed
    by (``rand/randmts.c:126``), so ``offset == 0`` drives ``seed1`` to 0 and
    exercises the ``seed1 += 2`` fixup; ``offset == -1`` and ``offset == 1``
    straddle it.  Passing the seed as an offset keeps the transcript's ``args``
    field from carrying a 6000-digit literal.
    """
    seed = Integer(2) ** 19937 - 20027 + Integer(offset)
    set_random_seed(seed)
    return _ints([sage_random() for _ in range(int(count))])


def rs_initial_seed(seed):
    """``initial_seed()`` -- must not move as draws are taken."""
    set_random_seed(Integer(seed))
    for _ in range(3):
        sage_random()
    return str(initial_seed())


def rs_urandomm(seed, n, count):
    """``mpz_urandomm(state, n)`` as reached through ``ZZ.random_element(n)``.

    ``ZZ._randomize_mpz`` unconditionally burns one ``c_random()`` on its
    ``den`` local before dispatching (``integer_ring.pyx:801``), so each draw
    here is *one* 31-bit discard followed by ``mpz_urandomm``.
    """
    set_random_seed(Integer(seed))
    n = Integer(n)
    return _ints([ZZ.random_element(n) for _ in range(int(count))])


def rs_zz_seed(seed, count):
    """``randstate.ZZ_seed()`` == ``ZZ.random_element(1 << 128)``.

    Reference: ``randstate.pyx:629-642``.  This is the 128-bit draw Sage uses
    to seed every foreign generator, including ``python_random``.
    """
    set_random_seed(Integer(seed))
    rstate = current_randstate()
    return _ints([rstate.ZZ_seed() for _ in range(int(count))])


def rs_python_random(seed, count):
    """CPython ``random.Random.random()`` off ``randstate.python_random()``."""
    set_random_seed(Integer(seed))
    rnd = current_randstate().python_random()
    return _fs([rnd.random() for _ in range(int(count))])


def rs_python_getrandbits(seed, k, count):
    set_random_seed(Integer(seed))
    rnd = current_randstate().python_random()
    return _ints([rnd.getrandbits(int(k)) for _ in range(int(count))])


def rs_python_randrange(seed, lo, hi, count):
    set_random_seed(Integer(seed))
    rnd = current_randstate().python_random()
    return _ints([rnd.randrange(int(lo), int(hi)) for _ in range(int(count))])


def rs_python_randint(seed, a, b, count):
    set_random_seed(Integer(seed))
    rnd = current_randstate().python_random()
    return _ints([rnd.randint(int(a), int(b)) for _ in range(int(count))])


def rs_python_normalvariate(seed, mu_n, mu_d, sigma_n, sigma_d, count):
    set_random_seed(Integer(seed))
    rnd = current_randstate().python_random()
    mu = _rat(mu_n, mu_d)
    sigma = _rat(sigma_n, sigma_d)
    return _fs([rnd.normalvariate(mu, sigma) for _ in range(int(count))])


def rs_python_random_reseed(seed, explicit_seed, count):
    """``python_random()`` first, then ``python_random(seed=...)``.

    Upstream returns the *cached* object and ignores the new seed
    (``randstate.pyx:617-620``: the ``type(self._python_random) is cls`` early
    return happens before ``seed`` is looked at), so the second stream simply
    continues the first.
    """
    set_random_seed(Integer(seed))
    rstate = current_randstate()
    first = [rstate.python_random().random() for _ in range(int(count))]
    second = [
        rstate.python_random(seed=int(explicit_seed)).random()
        for _ in range(int(count))
    ]
    return _fs(first) + ' ' + _fs(second)


def rs_python_random_fresh_seed(seed, explicit_seed, count):
    """``python_random(seed=...)`` as the *first* call on a fresh randstate."""
    set_random_seed(Integer(seed))
    rnd = current_randstate().python_random(seed=int(explicit_seed))
    return _fs([rnd.random() for _ in range(int(count))])


# ===========================================================================
# sage.rings.integer_ring -- the randstate consumers
# ===========================================================================

def zz_random_element_1n(seed, count):
    """``ZZ.random_element()`` (the default ``1/n`` distribution)."""
    set_random_seed(Integer(seed))
    return _ints([ZZ.random_element() for _ in range(int(count))])


def zz_random_element_range(seed, lo, hi, count):
    """``ZZ.random_element(lo, hi)`` -- uniform on ``[lo, hi)``."""
    set_random_seed(Integer(seed))
    lo, hi = Integer(lo), Integer(hi)
    return _ints([ZZ.random_element(lo, hi) for _ in range(int(count))])


def zz_random_element_rrandomb(seed, bits, count):
    """``ZZ.random_element(bits, distribution='mpz_rrandomb')``.

    ``mpz_rrandomb`` is GMP's *runs* generator (long strings of equal bits),
    not ``mpz_urandomb``.
    """
    set_random_seed(Integer(seed))
    return _ints([
        ZZ.random_element(int(bits), distribution='mpz_rrandomb')
        for _ in range(int(count))
    ])


def zz_random_element_gaussian(seed, sigma_n, sigma_d, count):
    """``ZZ.random_element(sigma, distribution='gaussian')``.

    Routes into ``DiscreteGaussianDistributionIntegerSampler(algorithm=
    'uniform+logtable')`` and caches the sampler in a module global
    (``integer_ring.pyx:826-834``).
    """
    set_random_seed(Integer(seed))
    sigma = _rat(sigma_n, sigma_d)
    return _ints([
        ZZ.random_element(sigma, distribution='gaussian')
        for _ in range(int(count))
    ])


# ===========================================================================
# sage.stats.distributions.discrete_gaussian_integer
# ===========================================================================

def zz_error(nargs, x, y, dist):
    """Argument validation of ``ZZ.random_element`` (``integer_ring.pyx:768-781``).

    ``nargs`` is 0, 1 or 2 and says how many of ``x``/``y`` are actually
    passed.  Only the *error surface* is compared here -- the drawn value is
    covered by the ``zz_random_element_*`` functions.
    """
    distribution = {
        0: None,
        1: 'uniform',
        2: '1/n',
        3: 'mpz_rrandomb',
        4: 'gaussian',
        5: 'bogus',
    }[int(dist)]
    nargs = int(nargs)
    args = []
    if nargs >= 1:
        args.append(Integer(x))
    if nargs >= 2:
        args.append(Integer(y))
    set_random_seed(1)
    try:
        ZZ.random_element(*args, distribution=distribution)
    except Exception as e:
        return '%s: %s' % (type(e).__name__, e)
    return 'ok'


def rs_python_error(seed, lo, hi):
    """``randrange`` on an empty range.  Both sides must raise ``ValueError``."""
    set_random_seed(Integer(seed))
    rnd = current_randstate().python_random()
    rnd.randrange(int(lo), int(hi))
    return 'ok'


def dgi_samples(seed, sigma_n, sigma_d, tau, alg, count):
    set_random_seed(Integer(seed))
    D = DGI(sigma=_rat(sigma_n, sigma_d), tau=int(tau), algorithm=_alg(alg))
    return _ints([D() for _ in range(int(count))])


def dgi_samples_c(seed, sigma_n, sigma_d, c_n, c_d, alg, count):
    set_random_seed(Integer(seed))
    D = DGI(sigma=_rat(sigma_n, sigma_d), c=_rat(c_n, c_d), algorithm=_alg(alg))
    return _ints([D() for _ in range(int(count))])


def dgi_error(sigma_n, sigma_d, tau, c_n, c_d, alg):
    """Constructor validation of ``DiscreteGaussianDistributionIntegerSampler``.

    Returns the exception's class and message as text so that a *message*
    mismatch fails the case; ``compare.ts`` scores "both raised" as a pass
    regardless of what they raised.

    .. WARNING::

        Do not add ``sigma2+logtable`` cases with ``sigma < 0.43`` here.
        ``round(sigma/sigma_2)`` is then 0 and upstream's ``dgs`` calls
        ``mpz_urandomm`` with a zero modulus, which **segfaults the whole Sage
        runner** (verified on 10.3) and would wipe out every other case's
        transcript.  The port raises a ``ValueError`` there instead; that
        deviation is deliberately not exercised as a differential case.
    """
    try:
        D = DGI(
            sigma=_rat(sigma_n, sigma_d),
            c=_rat(c_n, c_d),
            tau=int(tau),
            algorithm=_alg(alg),
        )
    except Exception as e:
        return '%s: %s' % (type(e).__name__, e)
    return 'ok: ' + repr(D)


def dgi_default_algorithm(sigma_n, sigma_d, tau):
    """Which algorithm the ``algorithm=None`` default picks.

    ``sigma*tau <= table_cutoff (10^6)`` -> ``uniform+table``
    (``discrete_gaussian_integer.pyx:352-355``).
    """
    D = DGI(sigma=_rat(sigma_n, sigma_d), tau=int(tau))
    return D.algorithm


def dgi_repr(sigma_n, sigma_d, c_n, c_d):
    D = DGI(sigma=_rat(sigma_n, sigma_d), c=_rat(c_n, c_d))
    return repr(D)


def dgi_stream_mode(seed, sigma_n, sigma_d, count, mode):
    """The three streams of the ``_flush_cache`` doctest
    (``discrete_gaussian_integer.pyx:404-434``).

    mode 0 -- one seeding, ``count`` consecutive draws
    mode 1 -- reseed before every draw (shows the 32-bit sign-bit cache)
    mode 2 -- reseed *and* ``_flush_cache()`` before every draw
    """
    mode = int(mode)
    count = int(count)
    sigma = _rat(sigma_n, sigma_d)
    set_random_seed(Integer(seed))
    D = DGI(sigma=sigma)
    out = []
    for _ in range(count):
        if mode >= 1:
            set_random_seed(Integer(seed))
        if mode == 2:
            D._flush_cache()
        out.append(D())
    return _ints(out)


def dgi_histogram(seed, sigma_n, sigma_d, alg, count, lo, hi):
    """Counts of each value in ``[lo, hi]`` over ``count`` seeded draws.

    Exact (same seed both sides) *and* a distributional readout: the counts
    have to fall off like ``exp(-x^2/2 sigma^2)``.
    """
    set_random_seed(Integer(seed))
    D = DGI(sigma=_rat(sigma_n, sigma_d), algorithm=_alg(alg))
    lo, hi = int(lo), int(hi)
    counts = {v: 0 for v in range(lo, hi + 1)}
    outside = 0
    for _ in range(int(count)):
        v = int(D())
        if lo <= v <= hi:
            counts[v] += 1
        else:
            outside += 1
    return _ints([counts[v] for v in range(lo, hi + 1)]) + ' outside=%d' % outside


def dgi_support(seed, sigma_n, sigma_d, tau, alg, count):
    """``(min, max, |x| <= ceil(sigma*tau) for every draw)``.

    ``dgs`` bounds the support at ``ceil(sigma*tau)`` around ``round(c)``
    (``dgs_gauss_mp.c:161-200``); this pins that invariant on both sides.
    """
    import math
    sigma = _rat(sigma_n, sigma_d)
    set_random_seed(Integer(seed))
    D = DGI(sigma=sigma, tau=int(tau), algorithm=_alg(alg))
    xs = [Integer(D()) for _ in range(int(count))]
    bound = Integer(math.ceil(sigma * int(tau)))
    ok = all(abs(x) <= bound for x in xs)
    return '(%s, %s, %s, %s)' % (min(xs), max(xs), bound, 'True' if ok else 'False')


def dgi_rho(sigma_n, sigma_d, c_n, c_d, x):
    """``rho_{sigma,c}(x) = exp(-(x-c)^2 / (2 sigma^2))``."""
    sigma = _rat(sigma_n, sigma_d)
    c = _rat(c_n, c_d)
    x = float(Integer(x))
    return _f(RR(exp(RR(-((x - c) ** 2) / (2 * sigma * sigma)))))


# ===========================================================================
# sage.stats.distributions.discrete_gaussian_lattice
# ===========================================================================

def _identity(n):
    return identity_matrix(ZZ, int(n))


def dgl_samples(seed, n, sigma_n, sigma_d, count):
    """Identity basis, centre 0 -> ``_call_simple``."""
    DGL = _DGL()
    set_random_seed(Integer(seed))
    D = DGL(_identity(n), _rat(sigma_n, sigma_d))
    return _vecs([D() for _ in range(int(count))])


def dgl_samples_c(seed, n, sigma_n, sigma_d, c_flat, c_den, count):
    """Identity basis with centre ``c_flat / c_den``.

    ``c_den == 1`` keeps ``c`` in the lattice (``_call_simple``); anything else
    pushes it off the lattice and into the GPV recursion ``_call``.
    """
    DGL = _DGL()
    set_random_seed(Integer(seed))
    c_den = Integer(c_den)
    c = vector(QQ, [Integer(x) / c_den for x in c_flat])
    D = DGL(_identity(n), _rat(sigma_n, sigma_d), c=c)
    return _vecs([D() for _ in range(int(count))])


def dgl_samples_basis(seed, basis_flat, n, sigma_n, sigma_d, count):
    """Arbitrary (non-orthonormal) basis -> ``_call``, the GPV recursion."""
    DGL = _DGL()
    set_random_seed(Integer(seed))
    B = matrix(ZZ, _square(basis_flat, n))
    D = DGL(B, _rat(sigma_n, sigma_d))
    return _vecs([D() for _ in range(int(count))])


def dgl_repr(n, sigma_n, sigma_d, c_flat):
    DGL = _DGL()
    D = DGL(_identity(n), _rat(sigma_n, sigma_d), c=vector(ZZ, [Integer(x) for x in c_flat]))
    return repr(D)


def dgl_repr_sigma(sigma_flat, n):
    """``__repr__`` for a covariance matrix, and the ``RuntimeError`` upstream
    raises when it is not positive definite
    (``discrete_gaussian_lattice.py:569-570``)."""
    DGL = _DGL()
    try:
        return repr(DGL(_identity(n), matrix(ZZ, _square(sigma_flat, n))))
    except Exception as e:
        return '%s: %s' % (type(e).__name__, e)


def dgl_compute_precision(prec, sigma_prec):
    """``compute_precision(precision, sigma)`` (``discrete_gaussian_lattice.py:155``).

    ``prec == -1`` means ``None``; ``sigma_prec == -1`` means a plain
    ``Integer`` sigma (no ``.precision()``, so the ``AttributeError`` branch).
    """
    DGL = _DGL()
    prec = int(prec)
    sigma_prec = int(sigma_prec)
    sigma = Integer(3) if sigma_prec < 0 else RealField(sigma_prec)(3)
    return str(DGL.compute_precision(None if prec < 0 else prec, sigma))


def dgl_normalisation(n, sigma_n, sigma_d, tau, prec):
    """``_normalisation_factor_zz`` on ``ZZ^n`` (Poisson summation via qfrep)."""
    DGL = _DGL()
    D = DGL(_identity(n), _rat(sigma_n, sigma_d))
    tau = int(tau)
    prec = int(prec)
    value = D._normalisation_factor_zz(
        tau=None if tau < 0 else tau,
        prec=None if prec < 0 else prec,
    )
    return _f(value)


def dgl_normalisation_round(n, sigma_n, sigma_d, prec):
    """``round(_normalisation_factor_zz(prec=prec))`` -- exact integer readout."""
    DGL = _DGL()
    D = DGL(_identity(n), _rat(sigma_n, sigma_d))
    return str(Integer(round(D._normalisation_factor_zz(prec=int(prec)))))


def dgl_normalisation_basis(
    basis_flat, rows, cols, basis_den, sigma_n, sigma_d, c_flat, c_den
):
    """``_normalisation_factor_zz`` on a general basis / centre.

    Upstream raises ``NotImplementedError`` for a non-square basis, a
    non-integral lattice, or a centre that is not the origin of a trivial
    lattice (``discrete_gaussian_lattice.py:313-322``).
    """
    DGL = _DGL()
    rows, cols = int(rows), int(cols)
    basis_den = Integer(basis_den)
    entries = [
        [Integer(basis_flat[i * cols + j]) / basis_den for j in range(cols)]
        for i in range(rows)
    ]
    # Upstream's "lattice must be integral" test is ``B.base_ring() != ZZ``,
    # i.e. it is about the matrix's *parent*, not its entries: a QQ matrix of
    # integers still fails it.  The port has no parents and tests the entries,
    # so build over ZZ exactly when the entries are integral -- otherwise the
    # comparison would be about Sage's coercion model rather than about the
    # lattice sampler.
    B = matrix(ZZ if basis_den == 1 else QQ, entries)
    # An empty ``c_flat`` means "do not pass ``c`` at all", so that upstream's
    # default ``c=0`` is used.  Passing ``c=None`` explicitly is a different
    # thing: ``set_c(None)`` returns before ``_precompute_data()``
    # (``discrete_gaussian_lattice.py:778-780``) and leaves the sampler half
    # built, which the port's constructor has no way to express.
    kwargs = {}
    if len(c_flat) > 0:
        cd = Integer(c_den)
        kwargs['c'] = vector(QQ, [Integer(x) / cd for x in c_flat])
    try:
        D = DGL(B, _rat(sigma_n, sigma_d), **kwargs)
        return _f(D._normalisation_factor_zz())
    except Exception as e:
        # Returned as text (not re-raised) so that a *message* mismatch fails
        # the case; compare.ts scores "both raised" as a pass regardless.
        return '%s: %s' % (type(e).__name__, e)


def dgl_normalisation_nonspherical(sigma_flat, n, c_flat):
    """``_normalisation_factor_zz`` for a covariance matrix (the approximate
    branch, ``discrete_gaussian_lattice.py:297-312``)."""
    DGL = _DGL()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        D = DGL(
            _identity(n),
            matrix(ZZ, _square(sigma_flat, n)),
            c=vector(ZZ, [Integer(x) for x in c_flat]),
        )
        return _f(D._normalisation_factor_zz())


def dgl_f(n, sigma_n, sigma_d, c_flat, x_flat):
    """``D.f(x)`` -- the (unnormalised) Gaussian weight."""
    DGL = _DGL()
    D = DGL(
        _identity(n),
        _rat(sigma_n, sigma_d),
        c=vector(ZZ, [Integer(v) for v in c_flat]),
    )
    return _f(D.f([Integer(v) for v in x_flat]))


def dgl_f_sigma(sigma_flat, n, x_flat):
    """``D.f(x)`` for a covariance matrix: ``exp(-x Sigma^-1 x / 2)``."""
    DGL = _DGL()
    D = DGL(_identity(n), matrix(ZZ, _square(sigma_flat, n)))
    return _f(D.f([Integer(v) for v in x_flat]))


def dgl_maximal_r(sigma_flat, n):
    """``_maximal_r()``: the largest ``r`` with ``Sigma - r^2 Q`` positive definite."""
    DGL = _DGL()
    D = DGL(_identity(n), matrix(ZZ, _square(sigma_flat, n)))
    return _f(D._maximal_r())


def dgl_nonspherical(seed, sigma_flat, n, c_flat, c_den, count):
    """``_call_non_spherical`` -- Peikert's two-stage sampler.

    Its offline stage draws ``n`` ``normalvariate(0, 1)`` per sample from
    ``sage.misc.prandom``, i.e. from ``randstate.python_random()``.
    """
    DGL = _DGL()
    set_random_seed(Integer(seed))
    c_den = Integer(c_den)
    c = vector(QQ, [Integer(x) / c_den for x in c_flat])
    D = DGL(_identity(n), matrix(ZZ, _square(sigma_flat, n)), c=c)
    return _vecs([D() for _ in range(int(count))])


def dgl_poly_sampler(seed, n, sigma_n, sigma_d, count):
    """``DiscreteGaussianDistributionPolynomialSampler`` coefficient streams."""
    from sage.stats.distributions.discrete_gaussian_polynomial import (
        DiscreteGaussianDistributionPolynomialSampler as DGP,
    )
    set_random_seed(Integer(seed))
    n = int(n)
    P = PolynomialRing(ZZ, 'x')
    sampler = DGP(P, n, _rat(sigma_n, sigma_d))
    out = []
    for _ in range(int(count)):
        p = sampler()
        coeffs = p.list()
        coeffs = coeffs + [Integer(0)] * (n - len(coeffs))
        out.append(coeffs)
    return _vecs(out)


FUNCTIONS = {
    # randstate
    'rs_random': rs_random,
    'rs_random_modseed': rs_random_modseed,
    'rs_initial_seed': rs_initial_seed,
    'rs_urandomm': rs_urandomm,
    'rs_zz_seed': rs_zz_seed,
    'rs_python_random': rs_python_random,
    'rs_python_getrandbits': rs_python_getrandbits,
    'rs_python_randrange': rs_python_randrange,
    'rs_python_randint': rs_python_randint,
    'rs_python_normalvariate': rs_python_normalvariate,
    'rs_python_random_reseed': rs_python_random_reseed,
    'rs_python_random_fresh_seed': rs_python_random_fresh_seed,
    # integer_ring consumers
    'zz_random_element_1n': zz_random_element_1n,
    'zz_random_element_range': zz_random_element_range,
    'zz_random_element_rrandomb': zz_random_element_rrandomb,
    'zz_random_element_gaussian': zz_random_element_gaussian,
    'zz_error': zz_error,
    'rs_python_error': rs_python_error,
    # discrete_gaussian_integer
    'dgi_error': dgi_error,
    'dgi_samples': dgi_samples,
    'dgi_samples_c': dgi_samples_c,
    'dgi_default_algorithm': dgi_default_algorithm,
    'dgi_repr': dgi_repr,
    'dgi_stream_mode': dgi_stream_mode,
    'dgi_histogram': dgi_histogram,
    'dgi_support': dgi_support,
    'dgi_rho': dgi_rho,
    # discrete_gaussian_lattice
    'dgl_samples': dgl_samples,
    'dgl_samples_c': dgl_samples_c,
    'dgl_samples_basis': dgl_samples_basis,
    'dgl_repr': dgl_repr,
    'dgl_repr_sigma': dgl_repr_sigma,
    'dgl_compute_precision': dgl_compute_precision,
    'dgl_normalisation': dgl_normalisation,
    'dgl_normalisation_round': dgl_normalisation_round,
    'dgl_normalisation_basis': dgl_normalisation_basis,
    'dgl_normalisation_nonspherical': dgl_normalisation_nonspherical,
    'dgl_f': dgl_f,
    'dgl_f_sigma': dgl_f_sigma,
    'dgl_maximal_r': dgl_maximal_r,
    'dgl_nonspherical': dgl_nonspherical,
    'dgl_poly_sampler': dgl_poly_sampler,
}
