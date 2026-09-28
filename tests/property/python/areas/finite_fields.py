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
    'ff_ext_order': lambda p, n: GF(p**n, 'a').cardinality(),
    'ff_ext_frobenius': lambda p, n, a: GF(p**n, 'a').fetch_int(Integer(a)).frobenius().integer_representation(),
    'ff_generator_order': lambda p, n: GF(p**n, 'a').multiplicative_generator().multiplicative_order(),
}

import json
import math

def comparison(run):
    def normalized(x):
        if isinstance(x, (list, tuple)): return [normalized(y) for y in x]
        if isinstance(x, bool) or x is None: return x
        return str(x)
    try: return json.dumps({'value': normalized(run())}, separators=(',', ':'), ensure_ascii=False)
    except Exception as e: return json.dumps({'error': type(e).__name__, 'message': str(e)}, separators=(',', ':'), ensure_ascii=False)

def extension(p,n):
    return GF(p**n, 'a', impl='pari_ffelt', modulus=GF(p**n,'a').modulus())

def from_integer(F,a):
    # Keep the historical fetch_int spelling available for older Sage runtimes.
    return F.from_integer(a) if hasattr(F, 'from_integer') else F.fetch_int(a)

def extension_unary(method,p,n,a):
    x=from_integer(extension(p,n),a)
    if method=='neg': return -x
    if method=='isZero': return x==0
    if method=='isOne': return x==1
    if method=='coefficients': return list(vector(x))
    if method in ('toString','repr'): return str(x)
    if method in ('minpoly','minimal_polynomial','minimalPolynomial'): return x.minpoly().list()
    return getattr(x,method)()

for method in ['_integer_','neg','isZero','isOne','trace','norm','coefficients','integer_representation','toString','repr','minpoly','minimal_polynomial','minimalPolynomial']:
    FUNCTIONS['ff_ext_'+method] = lambda p,n,a,m=method: comparison(lambda: extension_unary(m,p,n,a))
FUNCTIONS['ff_ext_fromInteger'] = lambda p,n,a: comparison(lambda: from_integer(extension(p,n),a).integer_representation())
FUNCTIONS['ff_ext_iteration'] = lambda p,n: comparison(lambda: [x.integer_representation() for x in extension(p,n)])
def extension_parent(p,n):
    F=extension(p,n)
    return [F.zero(),F.one(),F.gen(),F.cardinality(),F.is_field(),str(F),F.characteristic(),F.degree(),F.order(),
      list(F),F.multiplicative_generator(),F.primitive_element(),F.multiplicative_generator()]
FUNCTIONS['ff_ext_parent'] = lambda p,n: comparison(lambda: extension_parent(p,n))

def scalar_constructor(p,n,mode,a,b):
    F=GF(p) if n<=1 else extension(p,n)
    if mode==0: x=Integer(a)
    elif mode==1: x=QQ(a)/b
    elif mode==2: x=bool(a)
    elif mode==3: x=None
    elif mode==4: x=float(a)/float(b) if b else (float('nan') if not a else float('inf') if a>0 else -float('inf'))
    elif mode==5: x=str(a)
    elif mode==6: x=[QQ(a)/b,1]
    elif mode==7: x=[]
    # TypeScript has a single number type; integral numbers denote integers
    # throughout this port (DESIGN.md), including an integral division result.
    if mode==4 and math.isfinite(x) and x.is_integer(): x=Integer(x)
    return str(F(x))
FUNCTIONS['ff_scalar_constructor'] = lambda *args: comparison(lambda: scalar_constructor(*args))

def prime_unary(method,p,a):
    x=GF(p)(a)
    if method=='neg': return -x
    if method=='inv': return ~x
    if method=='isZero': return x==0
    if method=='isOne': return x==1
    if method=='isUnit': return x.is_unit()
    if method in ('toString','repr'): return str(x)
    if method=='toBigInt': return Integer(x)
    if method=='sqrt': return x.sqrt(extend=False)
    return getattr(x,method)()
for method in ['_integer_','neg','inv','isZero','isOne','is_square','sqrt','toString','repr','toBigInt','isUnit','lift','multiplicative_order']:
    FUNCTIONS['ff_prime_'+method] = lambda p,a,*rest,m=method: comparison(lambda: prime_unary(m,p,a))
import operator
for method in ['add','sub','mul','div','eq','pow']:
    FUNCTIONS['ff_prime_'+method] = lambda p,a,b,legacy,m=method: comparison(lambda: getattr(operator, 'truediv' if m=='div' else m)(GF(p)(a),Integer(b)))

def extension_conversion(p,n,a,mode):
    F=extension(p,n)
    if mode==0: x=from_integer(extension(p,n),a)
    elif mode==1: x=from_integer(GF(p**n,'b',impl='pari_ffelt',modulus=F.modulus()),a)
    elif mode==2: x=from_integer(GF(p**(n+1),'b',impl='pari_ffelt',modulus=GF(p**(n+1),'b').modulus()),a)
    elif mode==3: x=GF(3 if p==2 else 2)(a)
    elif mode in (4,6): x=PolynomialRing(GF(p),'z')([a,1])
    elif mode==5: x=PolynomialRing(GF(3 if p==2 else 2),'z')([a,1])
    if mode==6:
        from sage.rings.finite_rings.element_pari_ffelt import FiniteFieldElement_pari_ffelt
        return str(FiniteFieldElement_pari_ffelt(F,x))
    return str(F(x))
FUNCTIONS['ff_ext_conversion'] = lambda *args: comparison(lambda: extension_conversion(*args))

from sage.rings.finite_rings.finite_field_prime_modn import FiniteField_prime_modn
FUNCTIONS['ff_prime_constructor'] = lambda p,legacy: comparison(lambda: str(FiniteField_prime_modn(p)))
def prime_parent(p,legacy):
    F=GF(p)
    return [F.zero(),F.one(),F.gen(),F.cardinality(),F.is_field(),str(F),F.characteristic(),F.degree(),F.order(),
      list(F),list(F),F.multiplicative_generator(),list(F) if legacy else F.primitive_element()]
FUNCTIONS['ff_prime_parent'] = lambda p,legacy: comparison(lambda: prime_parent(p,legacy))
def direct_element(p,a,b,mode,legacy):
    from sage.rings.finite_rings.integer_mod import IntegerMod
    value=Integer(a) if mode==0 else QQ(a)/b if mode==1 else None if mode==2 else bool(a)
    return str(IntegerMod(GF(p),value))
FUNCTIONS['ff_direct_element'] = lambda *args: comparison(lambda: direct_element(*args))
def field_random(p,n,seed):
    F=GF(p) if n<=1 else extension(p,n)
    set_random_seed(seed)
    if n<=1: return [str(F.random_element()) for _ in range(20)]
    # Bundled Sage 10.9 finite_field_base.pyx:1058-1059 changed the default
    # from the vector randomizer used by installed 10.3 to one randrange draw.
    from sage.misc.prandom import randrange
    return [str(from_integer(F,randrange(F.order()))) for _ in range(20)]
FUNCTIONS['ff_random'] = lambda *args: comparison(lambda: field_random(*args))
FUNCTIONS['ff_pari_order'] = lambda p,a,multiple: pari(Mod(a,p)).znorder((p-1)*multiple).sage()
for method in ['add','sub','mul','div','eq']:
    FUNCTIONS['ff_scalar_'+method] = lambda p,n,a,b,d,mode,m=method: comparison(lambda: getattr(operator, 'truediv' if m=='div' else m)((GF(p) if n<=1 else extension(p,n))(a),Integer(b) if mode==0 else QQ(b)/d if mode==1 else bool(b)))

FUNCTIONS['ff_prime_p'] = lambda p,a: GF(p)(a).parent().characteristic()
# Port-only convenience: first nonsquare in Sage's ordered prime-field iteration.
FUNCTIONS['ff_prime_quadratic_non_residue'] = lambda p: Integer(next(x for x in GF(p) if not x.is_square()))
for method in ['sub','div','eq']:
    FUNCTIONS['ff_ext_'+method] = lambda p,n,a,b,m=method: comparison(lambda: getattr(operator,'truediv' if m=='div' else m)(from_integer(extension(p,n),a),from_integer(extension(p,n),b)))

factory_algorithms=['conway','first_lexicographic','minimal_weight','adleman-lenstra','primitive','unknown']
def field_factory(q,mode,codes,coeffs):
    name=''.join(chr(c) for c in codes)
    # The bundled newer constructor now validates names for prime fields too.
    # Execute that source branch with Sage's own name validator on 10.3.
    if mode!=7 and q>=2:
        # 10.3 perfect_power incorrectly misses small even squares (issue 40846).
        # The bundled integer.pyx fixes this; its PARI fallback gives the correct pair.
        exponent, prime = pari(q).ispower()
        from sage.structure.category_object import certify_names, normalize_names
        if exponent==1: certify_names([name])
        else: normalize_names(1,name)
        if not Integer(prime).is_prime():
            raise ValueError('the order of a finite field must be a prime power')
    if mode in (0,5,6): F=GF(q,name)
    elif mode in (1,2): F=GF(q,name,modulus=list(coeffs))
    elif mode==3: F=GF(q,name,modulus=factory_algorithms[coeffs[0]])
    elif mode==4: F=GF(q,name,modulus=PolynomialRing(GF(coeffs[0]),'z')(list(coeffs[1:])))
    elif mode==7: F=GF(q,name,check=False)
    return [str(F),str(F.gen()),str(F.gen()**2),F.cardinality(),[] if F.degree()==1 else F.modulus().list()]
FUNCTIONS['ff_factory'] = lambda *args: comparison(lambda: field_factory(*args))

FUNCTIONS['ff_conway'] = lambda p,n: comparison(lambda: [exists_conway_polynomial(p,n),conway_polynomial(p,n).list()[:-1]])
FUNCTIONS['ff_perfect_power_backend'] = lambda n: comparison(lambda: [Integer(pari(abs(n)).ispower()[1]),Integer(pari(abs(n)).ispower()[0])] if abs(n)>1 else [abs(n),1])

# Port convenience helpers expose the same prime-power decomposition as Sage's factory.
def field_order_helpers(q):
    exponent, prime = pari(q).ispower() if q>=2 else (1,q)
    valid = q>=2 and Integer(prime).is_prime()
    return [valid, [q,Integer(prime),Integer(exponent),exponent==1] if valid else None]
FUNCTIONS['ff_order_helpers'] = lambda q: comparison(lambda: field_order_helpers(q))

def field_sqrt_options(p,a,mode,legacy):
    x=GF(p)(a)
    options={} if mode==0 else {'extend':False} if mode==1 else {'extend':True} if mode==2 else {'all':True} if mode==3 else {'extend':False,'all':True} if mode==4 else {'extend':True,'all':True} if mode==5 else {'unknown':True}
    try:
        result=x.sqrt(**options)
    except NotImplementedError as error:
        # Bundled integer_mod.pyx adds an explanatory message to the same
        # all-roots extension branch that raised an empty error in Sage 10.3.
        if str(error)=='' and options.get('all') and options.get('extend',True) and not x.is_square():
            raise NotImplementedError('Finding all square roots in extensions is not implemented; try extend=False to find only roots in the base ring Zmod(n).')
        raise
    if isinstance(result,list): return list(map(str,result))
    return [str(result),str(result**2),str(result.parent())]
FUNCTIONS['ff_sqrt_options'] = lambda *args: comparison(lambda: field_sqrt_options(*args))

def field_string_expression(p,n,codes,polynomial):
    F=extension(p,n)
    if p>31:
        # Match the port's documented PARI fallback outside its Conway table.
        modulus=PolynomialRing(GF(p),'x').irreducible_element(n,algorithm='adleman-lenstra')
        F=GF(p**n,'a',impl='pari_ffelt',modulus=modulus)
    source=''.join(chr(c) for c in codes)
    if polynomial==2: return str(F.polynomial_ring())
    return list(map(str,F.polynomial_ring()(source).list())) if polynomial else str(F(source))
FUNCTIONS['ff_string_expression'] = lambda *args: comparison(lambda: field_string_expression(*args))

def prime_power_scalar(p,a,mode,b,d,legacy):
    exponent=Integer(b) if mode==0 else QQ(b)/d if mode==1 else bool(b) if mode==2 else None if mode==3 else str(b) if mode==4 else float(b)/float(d) if d else float('nan') if not b else float('inf') if b>0 else -float('inf')
    if mode==6: exponent=[]
    if mode==5 and math.isfinite(exponent) and exponent.is_integer(): exponent=Integer(exponent)
    return GF(p)(a)**exponent
FUNCTIONS['ff_prime_power_scalar']=lambda *args:comparison(lambda:prime_power_scalar(*args))

def mixed_binary(method,kind,p,n,a,other_kind,q,m,b):
    F=extension(p,n) if kind==2 else Zmod(p) if kind==3 else GF(p)
    x=from_integer(F,a) if kind==2 else F(a)
    if other_kind==0: y=Integer(b)
    elif other_kind==1: y=QQ(b)/m
    elif other_kind==2: y=bool(b)
    elif other_kind==3: y=None
    elif other_kind==4: y=str(b)
    elif other_kind==5:
        y=float(b)/float(m) if m else float('nan') if not b else float('inf')
        if math.isfinite(y) and y.is_integer(): y=Integer(y)
    elif other_kind==6: y=Zmod(q)(b)
    elif other_kind in (7,10): y=GF(q)(b)
    elif other_kind==8: y=from_integer(extension(q,m),b)
    elif other_kind==9: y=[] if b==0 else [Integer(b)]
    result=getattr(operator,'truediv' if method=='div' else method)(x,y)
    return result if isinstance(result,bool) else [str(result),str(parent(result))]
for method in ['add','sub','mul','div','eq']:
    FUNCTIONS['ff_mixed_'+method]=lambda *args,meth=method:comparison(lambda:mixed_binary(meth,*args))


def sequence_multiply(kind,p,a,mode):
    R=Zmod(p) if kind==3 else extension(p,2) if kind==2 else GF(p)
    sequence=['', 'ab', '𝄞é', [], [Integer(2),Integer(7)]][int(mode)]
    result=R(a)*sequence
    return [result,str(parent(result))]
FUNCTIONS['ff_sequence_multiply']=lambda *args:comparison(lambda:sequence_multiply(*args))


def parent_scalar(mode,a,b):
    if mode==0: return Integer(a)
    if mode==1: return QQ(a)/b
    if mode==2: return bool(a)
    if mode==3: return None
    if mode==4: return ''
    if mode==5: return str(a)
    if mode==6:
        value=float(a)/float(b) if b else float('nan') if not a else float('inf') if a>0 else -float('inf')
        return Integer(value) if value.is_integer() else value
    if mode==7: return [] if not a else [Integer(a)]
    if mode==8: return Zmod(b)(a)
    if mode==9: return GF(b)(a)
    if mode==10: return extension(b,2).from_integer(a)

def generator_index(kind,p,mode,a,b):
    R=QQ if kind==4 else Zmod(p) if kind==3 else extension(p,2) if kind==2 else GF(p)
    return R.gen(parent_scalar(mode,a,b))

def parent_constructor(kind,mode,a,b,check):
    from sage.rings.finite_rings.finite_field_prime_modn import FiniteField_prime_modn
    R=FiniteField_prime_modn(parent_scalar(mode,a,b),check=bool(check) if kind else True)
    return [str(R),R.characteristic(),R.order(),R.zero(),R.one(),R.is_field()]

FUNCTIONS['ff_generator_index']=lambda *args:comparison(lambda:generator_index(*args))
FUNCTIONS['ff_parent_constructor']=lambda *args:comparison(lambda:parent_constructor(*args))


def integer_conversion(kind,p,n,a):
    x=from_integer(extension(p,n),a) if kind==2 else GF(p)(a)
    return [ZZ(x),Integer(x)]
FUNCTIONS['ff_integer_conversion']=lambda *args:comparison(lambda:integer_conversion(*args))

# Invoke the exported PARI C kernels directly; quotient-polynomial power is not
# restricted to irreducible moduli. Preserve zero/constant polynomial GEN types.
def pari_quotient_kernel(method,p,coeffs,modulus,exponent):
    from pari_prime_quotient_cache import native_pari_prime_quotient_cache
    result=native_pari_prime_quotient_cache(1 if method=='inv' else 0,int(p),int(exponent),list(map(int,modulus)),list(map(int,coeffs)))
    return json.loads(result)

def quotient_comparison(method,*args):
    return comparison(lambda:pari_quotient_kernel(method,*args))
FUNCTIONS['ff_pari_quotient_pow']=lambda *args:quotient_comparison('pow',*args)
FUNCTIONS['ff_pari_quotient_powBig']=lambda *args:quotient_comparison('pow',*args)
FUNCTIONS['ff_pari_quotient_inv']=lambda *args:quotient_comparison('inv',*args)

def extension_power_scalar(p,n,a,mode,b,d):
    x=from_integer(extension(p,n),a)
    return x**parent_scalar(mode,b,d)
FUNCTIONS['ff_extension_power_scalar']=lambda *args:comparison(lambda:extension_power_scalar(*args))

def quotient_inverse_lift(p,degree,coeffs,modulus):
    from pari_prime_quotient_cache import native_pari_prime_quotient_cache
    return json.loads(native_pari_prime_quotient_cache(3,int(p),int(degree),list(map(int,modulus)),list(map(int,coeffs))))
FUNCTIONS['ff_pari_inverse_lift']=lambda *args:comparison(lambda:quotient_inverse_lift(*args))

from pari_ff_square_root import pari_ff_square_root, ff_extension_sqrt
FUNCTIONS['pari_ff_square_root'] = pari_ff_square_root
FUNCTIONS['ff_extension_sqrt'] = ff_extension_sqrt

from pari_field_predicates import pari_field_predicates, ff_extension_is_square
FUNCTIONS['pari_field_predicates'] = pari_field_predicates
FUNCTIONS['ff_extension_is_square'] = ff_extension_is_square
from pari_field_predicates import pari_field_record_scalar
FUNCTIONS['pari_field_record_scalar'] = pari_field_record_scalar

from pari_field_predicates import pari_field_predicates as pari_field_trace
FUNCTIONS['pari_field_trace'] = pari_field_trace
from pari_field_predicates import ff_extension_trace
FUNCTIONS['ff_extension_trace'] = ff_extension_trace

def pari_field_trace_scalar(p, a, seed):
    return pari_field_trace(8, p, [a], [0, 1], seed)
FUNCTIONS['pari_field_trace_scalar'] = pari_field_trace_scalar


def pari_field_charpoly(op, p, a, T, seed):
    return pari_field_trace(op, p, a, T, seed)

def pari_bivariate_resultant(op, p, T, coefficients, width, seed):
    Q = [coefficients[i:i+int(width)] for i in range(0,len(coefficients),int(width))]
    return pari_field_trace(op, p, T, Q, seed)

FUNCTIONS.update(pari_field_charpoly=pari_field_charpoly, pari_bivariate_resultant=pari_bivariate_resultant)

from pari_field_predicates import ff_extension_norm
FUNCTIONS['ff_extension_norm'] = ff_extension_norm

FUNCTIONS['pari_field_charpoly_scalar'] = lambda p,a,seed: pari_field_trace(14,p,[a],[0,1],seed)
