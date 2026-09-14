"""Live SageMath side of the September audit; see AUDIT-2026-09.md for source anchors."""
import json
from sage.all import GF, ZZ, PolynomialRing, CRT, discrete_log


def run(operation, args):
    if operation == 'modulus':
        p = args[0]
        F = GF(p**2, 'a', modulus=[p - 1, 0, 1])
        u, v = F.gen() - 1, F.gen() + 1
        return [bool(F.is_field()), u == 0, v == 0, u*v == 0]
    if operation.startswith('legacy_field_'): operation = operation.removeprefix('legacy_')
    if operation == 'extension_eq': return bool(GF(args[0], 'a').gen() == GF(args[1], 'b').gen())
    if operation == 'extension_scalar_eq': return bool(GF(args[0], 'a')(args[1]) == ZZ(args[1]))
    if operation == 'constant_polynomial_eq':
        a,b,c,d=args
        return bool(PolynomialRing(GF(a),'x')(c) == PolynomialRing(GF(b),'x' if d == 0 else 'y')(c))
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
    if operation in ('squarefree','factor_product','squarefree_inseparable'):
        p, unit = args
        R = PolynomialRing(GF(p), 'x')
        f = unit*(R.gen() + 1)**(p if operation == 'squarefree_inseparable' else 2)
        # Factorization.value() includes the unit: compare mathematical reconstruction,
        # independent of Sage's Factorization vs the port's array representation.
        return bool((f.factor() if operation == 'factor_product' else f.squarefree_decomposition()).value() == f)
    if operation == 'crt': return str(CRT(*map(ZZ, args)))
    if operation == 'additive_log':
        q, p = args
        x = GF(q, 'a').gen()
        return str(discrete_log(x + x, x, ZZ(p), operation='+'))
    if operation == 'invert_error': return str(~ZZ(args[0]))
    if operation in ('integer_pow','pow_error'): return str(ZZ(args[0])**ZZ(args[1]))
    if operation == 'exact_log': return str(ZZ(args[0]).exact_log(ZZ(args[1])))
    raise ValueError('Unknown operation: ' + operation)


# Encode exceptions as values, avoiding the comparison harness's any-two-errors rule.
def compare(operation, args):
    try: return json.dumps(dict(value=run(operation, args)), separators=(',', ':'))
    except Exception as error: return json.dumps(dict(error=type(error).__name__, **({'message': str(error)} if operation.endswith('_error') else {})), separators=(',', ':'))
FUNCTIONS = {op: (lambda *args, op=op: compare(op, args)) for op in ['invert_error','pow_error',
    'legacy_field_eq','legacy_field_add','extension_eq','extension_scalar_eq','constant_polynomial_eq','factor_product','squarefree_inseparable',
    'modulus', 'field_add', 'field_div', 'field_eq', 'frobenius', 'polynomial_eq',
    'polynomial_sub', 'squarefree', 'crt', 'additive_log', 'integer_pow', 'exact_log']}
