"""SageMath side of the ``finite_fields`` property-test area.

Cases: tests/property/cases/finite_fields.cases.json
"""

from sage.all import *

FUNCTIONS = {
    'ff_add': lambda p, a, b: GF(p)(a) + GF(p)(b),
    'ff_mul': lambda p, a, b: GF(p)(a) * GF(p)(b),
    'ff_inv': lambda p, a: GF(p)(a)**(-1),
    'ff_pow': lambda p, a, n: GF(p)(a)**n,
    'ff_pow_neg': lambda p, a, n: GF(p)(a)**n,
    'sqrt_mod': lambda a, p: (lambda F, x: (lambda r: p - r if r > (p-1)//2 else r)(Integer(x.sqrt())) if x.is_square() else None)(GF(p), GF(p)(a)),
    'sqrt_mod_p3mod4': lambda a, p: (lambda F, x: (lambda r: p - r if r > (p-1)//2 else r)(Integer(x.sqrt())) if x.is_square() else None)(GF(p), GF(p)(a)),
    'sqrt_mod_p5mod8': lambda a, p: (lambda F, x: (lambda r: p - r if r > (p-1)//2 else r)(Integer(x.sqrt())) if x.is_square() else None)(GF(p), GF(p)(a)),
    'sqrt_mod_general': lambda a, p: (lambda F, x: (lambda r: p - r if r > (p-1)//2 else r)(Integer(x.sqrt())) if x.is_square() else None)(GF(p), GF(p)(a)),
    'sqrt_mod_nonresidue': lambda a, p: (lambda F, x: (lambda r: p - r if r > (p-1)//2 else r)(Integer(x.sqrt())) if x.is_square() else None)(GF(p), GF(p)(a)),
    'primitive_root': primitive_root,
    'ff_multiplicative_generator': lambda p: Integer(GF(p).multiplicative_generator()),
    'discrete_log': lambda p, base, target: GF(p)(target).log(GF(p)(base)),
    'ff_ext_add': lambda p, n, a, b: (GF(p**n, 'a').fetch_int(Integer(a)) + GF(p**n, 'a').fetch_int(Integer(b))).integer_representation(),
    'ff_ext_mul': lambda p, n, a, b: (GF(p**n, 'a').fetch_int(Integer(a)) * GF(p**n, 'a').fetch_int(Integer(b))).integer_representation(),
    'ff_ext_inv': lambda p, n, a: (GF(p**n, 'a').fetch_int(Integer(a))**(-1)).integer_representation(),
    'ff_ext_pow': lambda p, n, a, e: (GF(p**n, 'a').fetch_int(Integer(a))**e).integer_representation(),
    'ff_ext_order': lambda p, n: p**n,
    'ff_ext_frobenius': lambda p, n, a: (GF(p**n, 'a').fetch_int(Integer(a))**p).integer_representation(),
    'ff_generator_order': lambda p, n: GF(p**n, 'a').multiplicative_generator().multiplicative_order(),
}
