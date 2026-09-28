"""Bundled hyperelliptic callers; preserve native scalar/polynomial root order."""
import ast
import copy
import importlib
import json
from pathlib import Path
from sage.all import QQ, GF, PolynomialRing, HyperellipticCurve
from finite_polynomial_roots import _roots_univariate_polynomial
from pari_ff_square_root import pari_ff_square_root
from curve_coordinate_coercion import _field

_root = Path(__file__).resolve().parents[3] / 'reference/sage/src/sage/schemes/hyperelliptic_curves'

def _method(filename, name):
    tree = ast.parse((_root / filename).read_text())
    return next(n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef) and n.name == name)

class _Bindings(ast.NodeTransformer):
    def visit_Call(self, node):
        self.generic_visit(node)
        if isinstance(node.func, ast.Attribute) and node.func.attr == 'roots':
            node.args.insert(0, node.func.value)
            node.func = ast.Name(id='_roots', ctx=ast.Load())
        return node

def hyp_root_callers(kind, coefficients, coordinate, operation, seed):
    trace = []
    try:
        K = _field(int(kind))
        def decode(n):
            if K is QQ or K.degree() == 1: return K(n)
            n = int(n) % int(K.order()); ds = []
            for _ in range(K.degree()): ds.append(n % K.characteristic()); n //= K.characteristic()
            return K(ds)
        encode = lambda v: str(v)
        R = PolynomialRing(K, 'x')
        f = R(list(map(decode, coefficients)))
        def roots(f, ring=None, multiplicities=True):
            trace.append(['roots', list(map(encode, f.list())), bool(multiplicities)])
            if K is QQ:
                from polynomial_roots import roots as rational_roots
                return rational_roots(f, multiplicities=multiplicities)
            return _roots_univariate_polynomial(K, f, None, multiplicities)
        if operation == 0:
            a = decode(coordinate)
            trace.append(['sqrt', encode(a), True, False])
            if K is not QQ and K.degree() > 1:
                modulus = list(map(int, K.modulus().list()))
                raw = json.loads(pari_ff_square_root(K.characteristic(), modulus, list(map(int, a.polynomial().list())), seed))
                if raw['value'] is None: ys = []
                else:
                    r = K(list(map(int, raw['value'])))
                    ys = [r] if not r or K.characteristic() == 2 else [r, -r]
                assert set(ys) == set(a.sqrt(all=True, extend=False))
            else: ys = a.sqrt(all=True, extend=False)
            result = {'value': list(map(encode, sorted(ys)))}
            if K is not QQ and K.degree() > 1: result['state'] = raw['state']
        elif operation == 3:
            result = {'value': list(map(encode, sorted([decode(coefficients[0]), decode(coordinate)])))}
        else:
            namespace = dict(_roots=roots, __name__='sage.schemes.hyperelliptic_curves.hyperelliptic_generic')
            name = 'lift_x' if operation == 4 else 'odd_degree_model' if operation == 1 else 'cantor_reduction'
            filename = 'hyperelliptic_generic.py' if operation in (1, 4) else 'jacobian_morphism.py'
            node = _Bindings().visit(copy.deepcopy(_method(filename, name)))
            node.decorator_list = []
            exec(compile(ast.fix_missing_locations(ast.Module(body=[node], type_ignores=[])), str(_root / filename), 'exec'), namespace)
            if operation in (1, 4):
                # Use the bundled constructor's validation and installed parent classes.
                ctorpath = _root / 'constructor.py'
                tree = ast.parse(ctorpath.read_text())
                tree.body = [n for n in tree.body if isinstance(n, (ast.FunctionDef, ast.Import, ast.ImportFrom))]
                ns = dict(vars(importlib.import_module('sage.schemes.hyperelliptic_curves.constructor')))
                exec(compile(tree, str(ctorpath), 'exec'), ns)
                H = ns['HyperellipticCurve'](f, R.one() if operation == 4 else R.zero())
                if operation == 4:
                    points = namespace[name](H, decode(coordinate), all=True)
                    result = {'value': [[list(map(encode, P)), P.codomain() is H] for P in points]}
                else:
                    changed = namespace[name](H)
                    result = {'value': [list(map(encode, changed.hyperelliptic_polynomials()[0].list())), changed is H]}
            else:
                # Construct a valid even-degree Mumford relation, k = a*q.
                x = R.gen(); t = decode(coordinate)
                a = x**3 + x + 1; b = R(t); h = R(1)
                r = decode(coefficients[0])
                q = r*r*x**3 + x + 1
                f = a*q + h*b + b*b
                aa, bb = namespace[name](a,b,f,h,2)
                result = {'value': [list(map(encode, aa.list())), list(map(encode, bb.list()))]}
    except Exception as error:
        result = {'error': type(error).__name__, 'message': str(error)}
    result['calls'] = trace
    return json.dumps(result, separators=(',', ':'))
