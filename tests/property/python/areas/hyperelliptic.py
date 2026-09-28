"""Live-SageMath oracle for the ported hyperelliptic-curve surface."""

from sage.all import *


def _poly(K, coeffs):
    R = PolynomialRing(K, 'x')
    return R([K(c) for c in coeffs])


def _ints(xs):
    return ','.join(str(ZZ(x)) for x in xs)


def hyp_summary(p, f_coeffs, h_coeffs, extensions):
    K = GF(ZZ(p))
    f = _poly(K, f_coeffs)
    h = _poly(K, h_coeffs) if h_coeffs else 0
    H = HyperellipticCurve(f, h)
    counts = H.count_points(int(extensions))
    frob = H.frobenius_polynomial()
    return 'g=%s counts=%s frob=%s jac=%s' % (
        H.genus(),
        _ints(counts),
        _ints(frob.list()),
        frob(1),
    )


def hyp_cartier(p, f_coeffs):
    K = GF(ZZ(p))
    H = HyperellipticCurve(_poly(K, f_coeffs))
    cartier = ';'.join(_ints(row) for row in H.Cartier_matrix().rows())
    hasse_witt = ';'.join(_ints(row) for row in H.Hasse_Witt().rows())
    return 'C=%s HW=%s a=%s p=%s' % (
        cartier,
        hasse_witt,
        H.a_number(),
        H.p_rank(),
    )


FUNCTIONS = {
    'hyp_summary': hyp_summary,
    'hyp_cartier': hyp_cartier,
}


def hyp_base_cardinality(p,f,h,op):
    k=QQ if p==0 else GF(p)
    if op==0:return 'null'if k.cardinality()is infinity else str(k.cardinality())
    # Execute the bundled constructor to retain its updated validation messages.
    # The installed generic curve initializer matches the bundled homogenization.
    import ast,importlib
    from pathlib import Path
    global _audit_hyper_constructor
    if '_audit_hyper_constructor'not in globals():
        path=Path(__file__).resolve().parents[4]/'reference/sage/src/sage/schemes/hyperelliptic_curves/constructor.py'
        tree=ast.parse(path.read_text());ns=dict(vars(importlib.import_module('sage.schemes.hyperelliptic_curves.constructor')))
        tree.body=[n for n in tree.body if isinstance(n,(ast.FunctionDef,ast.Import,ast.ImportFrom))]
        exec(compile(tree,str(path),'exec'),ns)
        _audit_hyper_constructor=ns['HyperellipticCurve']
    R=PolynomialRing(k,'x');H=_audit_hyper_constructor(R(f),R(h))
    return str(H)+'|'+str(H.genus())+'|'+str(H.base_ring())
FUNCTIONS['hyp_base_cardinality']=hyp_base_cardinality

from hyperelliptic_root_callers import hyp_root_callers
FUNCTIONS['hyp_root_callers'] = hyp_root_callers
