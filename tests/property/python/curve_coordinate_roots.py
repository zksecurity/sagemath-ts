"""Execute bundled curve callers with native scalar and polynomial dependencies.

Only dependency bindings are adapted: AST calls to scalar predicates/roots
are traced, and extension roots use the bundled PARI version. No caller
arithmetic, branch, root ordering or error text is rewritten.
"""
import ast
import json
from pathlib import Path
from sage.all import GF, QQ, EllipticCurve, PolynomialRing
from sage.structure.coerce import py_scalar_to_element
from finite_polynomial_roots import _roots_univariate_polynomial
from pari_field_predicates import pari_field_predicates
from pari_ff_square_root import pari_ff_square_root

_path = Path(__file__).resolve().parents[3] / 'reference/sage/src/sage/schemes/elliptic_curves/ell_generic.py'
_tree = ast.parse(_path.read_text())
_methods = {n.name: n for c in _tree.body if isinstance(c, ast.ClassDef)
            for n in c.body if isinstance(n, ast.FunctionDef)}


class _ScalarBindings(ast.NodeTransformer):
    def visit_Call(self, node):
        self.generic_visit(node)
        if isinstance(node.func, ast.Attribute) and node.func.attr in ('is_square', 'sqrt'):
            node.args.insert(0, node.func.value)
            node.func = ast.Name(id='_scalar_' + node.func.attr, ctx=ast.Load())
        return node


def ec_coordinate_roots(p, degree, modulus, coefficients, coordinate, operation, seed):
    trace = []
    state = None
    try:
        if degree > 1:
            K = GF(p**degree, 'a', modulus=PolynomialRing(GF(p), 't')(modulus), impl='pari_ffelt')
            def decode(v):
                v = int(v); ds = []
                for _ in range(degree): ds.append(v % p); v //= p
                return K(ds)
            encode = lambda v: str(v.integer_representation())
            state = json.loads(pari_ff_square_root(p, modulus, [], seed))['state']
        else:
            K = GF(p) if p else QQ
            decode = lambda v: K(v[0])/K(v[1]) if isinstance(v,list) else K(v)
            encode = str
        E = EllipticCurve(K, list(map(decode, coefficients)))

        def predicate(x):
            trace.append(['is_square', encode(x)])
            if degree == 1: return x.is_square()
            return json.loads(pari_field_predicates(7, p, [int(c) for c in x.polynomial().list()], modulus, int(state)))['value']

        def sqrt(x, **options):
            nonlocal state
            trace.append(['sqrt', encode(x), bool(options.get('all', False))])
            if degree == 1: return x.sqrt(**options)
            raw = json.loads(pari_ff_square_root(p, modulus, [int(c) for c in x.polynomial().list()], int(state)))
            state = raw['state']
            if raw['value'] is None: return []
            r = K(list(map(int, raw['value'])))
            roots = [r] if not r or p == 2 else [r, -r]
            assert set(roots) == set(x.sqrt(all=True))
            return roots if options.get('all', False) else r

        class Poly:
            def __init__(self, f): self.f = f
            def roots(self, ring=None, multiplicities=True):
                trace.append(['roots', list(map(encode, self.f.list())), bool(multiplicities)])
                if K is QQ:
                    from polynomial_roots import roots
                    return roots(self.f, multiplicities=multiplicities)
                return _roots_univariate_polynomial(K, self.f, ring, multiplicities)

        namespace = dict(py_scalar_to_element=py_scalar_to_element,
                         PolynomialRing=lambda k, name: lambda cs: Poly(PolynomialRing(k, name)(cs)),
                         _scalar_is_square=predicate, _scalar_sqrt=sqrt)
        name = ['is_x_coord', 'lift_x', 'lift_x', 'montgomery_model', '_equation_string'][operation]
        import copy
        node = _ScalarBindings().visit(copy.deepcopy(_methods[name]))
        node.decorator_list = []
        exec(compile(ast.fix_missing_locations(ast.Module(body=[node], type_ignores=[])), str(_path), 'exec'), namespace)
        call = namespace[name]
        if operation == 0: value = bool(call(E, decode(coordinate)))
        elif operation in (1, 2):
            pts = call(E, decode(coordinate), all=operation == 1)
            if operation == 2: pts = [pts]
            value = [[encode(P[0]), encode(P[1]), P.curve() is E, all(c.parent() is K for c in P)] for P in pts]
        elif operation == 3:
            value = list(map(encode, call(E).ainvs()))
        else:
            value = "Elliptic Curve defined by " + call(E) + " over " + str(K)
        result = {'value': value}
    except Exception as error:
        result = {'error': type(error).__name__, 'message': str(error)}
    result['calls'] = trace
    # Polynomial dependencies have separate state oracles. This caller test
    # compares scalar state only where there are no polynomial root calls.
    if degree > 1 and p != 2 and operation < 3: result['state'] = state
    return json.dumps(result, separators=(',', ':'))
