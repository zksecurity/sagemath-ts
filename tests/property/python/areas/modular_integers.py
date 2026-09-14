"""Direct residue-ring comparisons, including result parents and exact errors."""
from sage.all import *
from sage.rings.finite_rings.integer_mod import IntegerMod
from sage.rings.finite_rings.integer_mod_ring import IntegerModRing_generic
import json, operator

def comparison(run):
    def normalize(x):
        if isinstance(x,(list,tuple)): return [normalize(y) for y in x]
        if isinstance(x,bool) or x is None: return x
        return str(x)
    try: return json.dumps({'value':normalize(run())},separators=(',',':'),ensure_ascii=False)
    except Exception as e: return json.dumps({'error':type(e).__name__,'message':str(e)},separators=(',',':'),ensure_ascii=False)

def scalar(mode,a,b):
    if mode==0: return Integer(a)
    if mode==1: return QQ(a)/b
    if mode==2: return bool(a)
    if mode==3: return None
    if mode==4: return str(a)
    if mode==5: return float(a)/float(b) if b else float('nan') if not a else float('inf') if a>0 else -float('inf')
    if mode==6: return Zmod(b)(a)
    if mode==7: return GF(b)(a)
    if mode==8: return []

def binary(method,n,a,mode,b,d):
    x=Zmod(n)(a); y=scalar(mode,b,d)
    if mode==5 and y.is_integer(): y=Integer(y) # Integral JS numbers denote integers.
    z=getattr(operator,'truediv' if method=='div' else method)(x,y)
    return z if isinstance(z,bool) else [z,str(parent(z))]

def constructor(n,mode,a,b,direct):
    x=scalar(mode,a,b)
    if mode==5 and isinstance(x,float) and x.is_integer(): x=Integer(x)
    return IntegerMod(Zmod(n),x) if direct else Zmod(n)(x)

def unary(method,n,a):
    x=Zmod(n)(a)
    if method=='neg': return -x
    if method=='inv': return ~x
    if method=='isZero': return x==0
    if method=='isOne': return x==1
    if method=='isUnit': return x.is_unit()
    if method=='modulus': return x.modulus()
    if method=='toBigInt': return Integer(x)
    if method=='repr': return 'Mod(%s, %s)'%(x,n) # Documented debugging convenience.
    if method=='toString': return str(x)
    return getattr(x,method)()

def ring(n):
    R=Zmod(n)
    return [str(R),R.characteristic(),R.order(),R.zero(),R.one(),R.gen(),R.is_field(),R.cardinality(),list(R),list(R),R.list_of_elements_of_multiplicative_group(),R.multiplicative_group_is_cyclic(),R.unit_gens()]

def random(n,seed,bound,mode):
    set_random_seed(seed); R=Zmod(n)
    return [R.random_element() if mode==0 else R.random_element(bound) for _ in range(20)]

FUNCTIONS={
 'mi_constructor':lambda *args:comparison(lambda:constructor(*args)),
 'mi_ring':lambda n:comparison(lambda:ring(n)),
 'mi_random':lambda *args:comparison(lambda:random(*args)),
 'mi_ring_direct':lambda n:comparison(lambda:str(IntegerModRing_generic(n))),
}
for method in ['add','sub','mul','div','eq','pow']:
    FUNCTIONS['mi_'+method]=lambda *args,m=method:comparison(lambda:binary(m,*args))
for method in ['neg','inv','isZero','isOne','isUnit','modulus','lift','toBigInt','toString','repr','multiplicative_order']:
    FUNCTIONS['mi_'+method]=lambda *args,m=method:comparison(lambda:unary(m,*args))

FUNCTIONS['mi_pari_order']=lambda a,n,order:comparison(lambda:pari(Mod(a,n)).znorder(None if order==0 else order).sage())

FUNCTIONS['mi_string']=lambda n,codes:comparison(lambda:Zmod(n)(''.join(chr(c) for c in codes)))


def factory(mode,kind,a,b):
    order=scalar(kind,a,b)
    factory=[Zmod,Integers,IntegerModRing][mode] if mode<2 else Zmod
    R=factory() if kind==9 else factory(order)
    return [str(R),R is ZZ,R.characteristic(),R(12)]

def mod_factory(n,kind,a,b,parent_order,with_parent):
    from sage.rings.finite_rings.integer_mod import Mod
    value=scalar(kind,a,b)
    # Integral JavaScript numbers denote integers throughout this port.
    if isinstance(value,float) and value.is_integer(): value=Integer(value)
    result=Mod(value,n,Zmod(parent_order)) if with_parent else Mod(value,n)
    if n==0: return [result,result is value]
    return [result,str(result.parent())]
FUNCTIONS['mi_factory']=lambda *args:comparison(lambda:factory(*args))
FUNCTIONS['mi_mod_factory']=lambda *args:comparison(lambda:mod_factory(*args))

FUNCTIONS['mi_factory_identity']=lambda n:comparison(lambda:[Zmod(n) is Zmod(n),Zmod(n) is Zmod(-n),Zmod(n) is Zmod(QQ(n)),Zmod(n) is Integers(n),Mod(1,n).parent() is Zmod(n)])


def mod_prime_parent(n,p,kind,a,b,legacy):
    result=Mod(scalar(kind,a,b),n,GF(p))
    if n==0: return result
    return [str(result),str(result.parent()),result.is_square()]
FUNCTIONS['mi_mod_prime_parent']=lambda *args:comparison(lambda:mod_prime_parent(*args))

from modular_polynomial_roots import mi_polynomial_roots
FUNCTIONS['mi_polynomial_roots']=lambda *args:comparison(lambda:mi_polynomial_roots(*args))
from modular_polynomial_roots import mi_field,mi_factored_order,mi_residue_root_lift
for _name,_function in [('mi_field',mi_field),('mi_factored_order',mi_factored_order),('mi_residue_root_lift',mi_residue_root_lift)]:
    FUNCTIONS[_name]=lambda *args,f=_function:comparison(lambda:f(*args))

FUNCTIONS['mi_modular_roots_hook']=lambda *args:comparison(lambda:mi_polynomial_roots(*args))
