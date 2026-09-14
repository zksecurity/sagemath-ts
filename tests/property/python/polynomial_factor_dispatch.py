"""Sage public factor shape and irreducibility comparisons.

Only the documented array/unit adaptation is normalized. Factor coefficients,
exponents and native order remain observable. For word finite fields execute
the bundled factor wrapper: installed Sage 10.3 lacks its zero guard.
"""
from pathlib import Path
import textwrap
from sage.all import ZZ, QQ, GF, PolynomialRing

_source_path = Path(__file__).resolve().parents[3] / 'reference/sage/src/sage/rings/polynomial/polynomial_zmod_flint.pyx'
_source = _source_path.read_text()
_start = _source.index('    def factor(self):')
_end = _source.index('    def monic(self):', _start)
factor_helper = lambda f: f.factor()
exec(compile(textwrap.dedent(_source[_start:_end]), str(_source_path), 'exec'))
_bundled_word_factor = factor


def polynomial_factor_fidelity(op, kind, coefficients, denominator):
    B = ZZ if kind == 0 else QQ if kind == 1 else GF(7) if kind == 2 else GF(2) if kind == 3 else GF(2**89-1)
    R = PolynomialRing(B, 'x')
    f = R([QQ(c)/denominator for c in coefficients])
    if op == 1:
        return f.is_irreducible()
    F = _bundled_word_factor(f) if kind == 2 else f.factor()
    pairs = [(R(g), e) for g, e in F]
    if F.unit() != 1:
        pairs.append((R(F.unit()), 1))
    # Inserting the separately held unit into the port's array representation
    # uses the same degree/exponent/polynomial order as Factorization.sort.
    pairs.sort(key=lambda t: (t[0].degree(), t[1], t[0]))
    return [[[str(c) for c in g.list()], str(e)] for g, e in pairs]


def polynomial_factor_state(op, kind, coefficients, denominator):
    """Native values plus bundled PARI effects, including cached false results."""
    import json
    from pari_qx_factor import pari_qx_factor
    R = PolynomialRing(ZZ if kind == 0 else QQ, 'x')
    f = R([QQ(c)/denominator for c in coefficients])
    if op == 0:
        value = polynomial_factor_fidelity(0, kind, coefficients, denominator)
        a = int(kind == 0 and 30 <= f.degree() <= 300)
        states, _ = json.loads(pari_qx_factor(19, 2, a, int(kind == 1), 1,
                                            denominator, coefficients, [], []))
        return json.dumps([value, states[0]], separators=(',', ':'))
    value = f.is_irreducible()
    assert f.is_irreducible() is value
    # QQ's specialized method factors the positive-leading primitive numerator
    # in ZZ[x], unlike QQ.factor(), which always delegates to PARI.
    g = f.numerator() if kind == 1 else f
    if kind == 1:
        content = abs(g.content())
        if g.leading_coefficient() < 0:
            content = -content
        g = g // content
    a = int(30 <= g.degree() <= 300)
    states, _ = json.loads(pari_qx_factor(19, 2, a, 0, 2, 1,
                                        [int(c) for c in g.list()], [], []))
    return json.dumps([[value, states[0]], [value, states[0]],
                       [R(f.list()).is_irreducible(), states[1]]], separators=(',', ':'))
