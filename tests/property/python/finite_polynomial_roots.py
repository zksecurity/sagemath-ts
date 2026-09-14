"""Ordered finite-field roots using the bundled Sage field method.

Extension coefficients are base-p integer representations with an explicit
modulus and the documented PARI element backend. No roots are re-sorted.
"""
from pathlib import Path
import textwrap
from sage.all import GF, PolynomialRing
from polynomial_factor_dispatch import _bundled_word_factor

_path = Path(__file__).resolve().parents[3] / 'reference/sage/src/sage/rings/finite_rings/finite_field_base.pyx'
_source = _path.read_text()
_start = _source.index('    def _roots_univariate_polynomial(')
_end = _source.index('\n\ndef unpickle_FiniteField_ext', _start)
exec(compile(textwrap.dedent(_source[_start:_end]).replace('pow(x, p, f)', '_bundled_power_mod(x, p, f)'), str(_path), 'exec'))


def _bundled_power_mod(x, exponent, modulus):
    # The bundled wrappers reject zero before native NTL/FLINT and accept large
    # Integer exponents. Installed 10.3 predates those changes. Execute the
    # original zero guard; use Sage's independent arithmetic power_mod for the
    # value of the modular power (no PARI randomness is consumed).
    from sage.arith.misc import power_mod
    K = modulus.base_ring()
    if not modulus:
        if K.degree() > 1 or K.characteristic() == 2:
            filename = 'polynomial_zz_pex.pyx' if K.degree() > 1 else 'polynomial_template.pxi'
            path = _path.parents[1] / 'polynomial' / filename
            line = next(line for line in path.read_text().splitlines()
                        if 'raise ZeroDivisionError("modulus must be nonzero")' in line)
            exec(compile(line.strip(), str(path), 'exec'))
        return x % modulus  # Bundled word/large-prime wrappers reduce first.
    return power_mod(x % modulus, exponent, modulus)


def finite_polynomial_roots(kind, p, degree, modulus, coefficients, factor_only=False, distinct=False):
    if kind == 3:
        K = GF(p**degree, 'a', modulus=PolynomialRing(GF(p), 't')(modulus), impl='pari_ffelt')
        def coefficient(c):
            digits = []
            for _ in range(degree):
                digits.append(c % p)
                c //= p
            return K(digits)
    else:
        K = GF(p)
        coefficient = K
    f = PolynomialRing(K, 'x')([coefficient(c) for c in coefficients])

    class FactorProxy:
        def __getattr__(self, name):
            return getattr(f, name)

        def factor(self):
            # Installed Sage 10.3 omits the bundled word-factor zero guard.
            return _bundled_word_factor(f) if degree == 1 and 2 < p < 2**63 else f.factor()

    if factor_only:
        F = FactorProxy().factor()
        pairs = [(g, e) for g, e in F]
        if F.unit() != 1:
            pairs.append((f.parent()(F.unit()), 1))
        pairs.sort(key=lambda t: (t[0].degree(), t[1], t[0]))
        return [[[str(c.integer_representation() if kind == 3 else int(c)) for c in g.list()], str(e)] for g,e in pairs]
    result = _roots_univariate_polynomial(K, f if distinct else FactorProxy(), None, not distinct)
    if distinct:
        return [[str(r.integer_representation() if kind == 3 else int(r)), r.parent() is K] for r in result]
    return [[str(r.integer_representation() if kind == 3 else int(r)), str(m), r.parent() is K]
            for r, m in result]


def finite_polynomial_roots_state(kind, p, coefficients):
    from pari_factor import native_pari_factor
    value = finite_polynomial_roots(kind, p, 1, [], coefficients)
    state = native_pari_factor(0, 1, p, coefficients).rsplit('|', 1)[1]
    return [value, state]


def finite_polynomial_distinct_roots_state(kind, p, coefficients):
    from pari_factor import native_pari_factor
    value = finite_polynomial_roots(kind, p, 1, [], coefficients, distinct=True)
    K=GF(p);R=PolynomialRing(K,'x');f=R(coefficients);x=R.gen()
    g=f.gcd(_bundled_power_mod(x,p,f)-x)
    state=native_pari_factor(0,1,p,[int(c) for c in g.list()]).rsplit('|',1)[1]
    return [value,state]
