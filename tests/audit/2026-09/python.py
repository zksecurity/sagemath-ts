"""Live SageMath side of the September audit; see AUDIT-2026-09.md for source anchors."""
import json
from pathlib import Path
from sage.all import GF, ZZ, PolynomialRing, CRT, discrete_log


def run(operation, args):
    if operation == 'modulus':
        p = args[0]
        F = GF(p**2, 'a', modulus=[p - 1, 0, 1])
        u, v = F.gen() - 1, F.gen() + 1
        return [bool(F.is_field()), u == 0, v == 0, u*v == 0]
    if operation.startswith('field_'):
        a, b = args
        x, y = GF(a)(1), GF(b)(1)
        if operation == 'field_add': return str(x + y)
        if operation == 'field_div': return str(x / GF(b)(2))
        if operation == 'field_eq': return bool(x == y)
    if operation == 'frobenius':
        q, k = args
        return str(GF(q, 'a').gen().frobenius(k).integer_representation())
    if operation.startswith('polynomial_'):
        F = GF(args[0])
        x, y = PolynomialRing(F, 'x').gen(), PolynomialRing(F, 'y').gen()
        if operation == 'polynomial_eq': return bool(x == y)
        return bool(x - y == 0)
    if operation == 'squarefree':
        p, unit = args
        R = PolynomialRing(GF(p), 'x')
        f = unit*(R.gen() + 1)**2
        # Factorization.value() includes the unit: compare mathematical reconstruction,
        # independent of Sage's Factorization vs the port's array representation.
        return bool(f.squarefree_decomposition().value() == f)
    if operation == 'crt': return str(CRT(*map(ZZ, args)))
    if operation == 'additive_log':
        q, p = args
        x = GF(q, 'a').gen()
        return str(discrete_log(x + x, x, ZZ(p), operation='+'))
    if operation == 'integer_pow': return str(ZZ(args[0])**ZZ(args[1]))
    if operation == 'exact_log': return str(ZZ(args[0]).exact_log(ZZ(args[1])))
    raise ValueError('Unknown operation: ' + operation)


results = []
for case in json.loads(Path(__file__).with_name('cases.json').read_text()):
    try:
        results.append(dict(id=case['id'], value=run(case['operation'], case['args'])))
    except Exception as error:
        results.append(dict(id=case['id'], error=type(error).__name__, message=str(error)))
print(json.dumps(results))
