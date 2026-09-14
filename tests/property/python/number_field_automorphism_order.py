"""Sage automorphism ordering with the bundled integral-model implementation.

Installed Sage 10.3 uses polredbest here. Execute the bundled replacement from
number_field.py so rational/negative scaling tests use the source being ported.
The embeddings(self) body below is the shared source's exact conversion/order.
"""
import ast
import json
from pathlib import Path
from sage.all import QQ, ZZ, PolynomialRing, NumberField
from pari_qx_factor import pari_qx_factor
from sage.rings.number_field.number_field import put_natural_embedding_first

_source = Path(__file__).resolve().parents[3] / 'reference/sage/src/sage/rings/number_field/number_field.py'
_tree = ast.parse(_source.read_text())
_method = next(node for node in ast.walk(_tree) if isinstance(node, ast.FunctionDef) and node.name == '_pari_absolute_structure')
_method.decorator_list = []
exec(compile(ast.Module(body=[_method], type_ignores=[]), str(_source), 'exec'))


def nf_automorphism_order(coefficients, denominator):
    R = PolynomialRing(QQ, 'x')
    K = NumberField(R([QQ(c) / denominator for c in coefficients]), 'a')
    f, alpha, _ = _pari_absolute_structure(K)
    # The public bundled galoisconj kernel accepts the integral polynomial;
    # an nf object supplies a denominator bound but has the same sorted roots.
    raw, _ = json.loads(pari_qx_factor(15, 2, 0, 0, 0, 1,
                                      [int(f.polcoef(i)) for i in range(K.degree() + 1)], [], []))
    coefficient = lambda c: QQ(c[0]) / QQ(c[1]) if isinstance(c, list) else QQ(c)
    conj = [sum((f.variable()**i) * coefficient(c) for i, c in enumerate(h)) for h in raw]
    P = alpha.lift()
    # K's installed PARI conversion map still uses polredbest. Apply the
    # bundled beta map explicitly: its integral generator is scalar*K.gen().
    original = K.absolute_polynomial()
    scalar = ZZ((original * original.denominator()).leading_coefficient())
    values = []
    for g in conj:
        image = P(g.Mod(f)).lift()
        values.append(K(R([QQ(image.polcoef(i)) * scalar**i for i in range(K.degree())])))
    values = sorted(values)
    maps = [K.hom([e]) for e in values]
    put_natural_embedding_first(maps)
    return json.dumps([[str(c) for c in sigma(K.gen()).list()] for sigma in maps], separators=(',', ':'))
