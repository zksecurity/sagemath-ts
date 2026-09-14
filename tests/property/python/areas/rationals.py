"""Sage rational methods and explicit operator/convenience adapters."""
import json
import operator
import struct
from sage.all import QQ, ZZ, matrix

UNARY = ['abs','neg','inv','toString','floor','ceil','round','trunc','isZero','isOne',
'isInteger','isPositive','isNegative','asIntegerRatio','height','continued_fraction_list',
'continued_fraction','is_square','period','real','imag','conjugate','norm','relative_norm',
'absolute_norm','trace','support','is_unit','is_integral','additive_order','multiplicative_order',
'simplify','list','signAsInt']
BINARY = ['add','sub','mul','div','eq','lt','le','gt','ge','cmp','rational_gcd','rational_lcm','content']
INTEGER = ['pow','valuation','ord','is_nth_power','nth_root','val_unit','padic_valuation',
'denominator_valuation','numerator_valuation','str']
SET = ['prime_to_S_part','is_S_unit','is_S_integral']

ALIASES = {k:getattr(operator,v) for k,v in dict(add='add',sub='sub',mul='mul',div='truediv',
neg='neg',pow='pow',eq='eq',lt='lt',le='le',gt='gt',ge='ge').items()}
ALIASES.update(dict(inv=lambda x:~x, cmp=lambda x,y:ZZ((x>y)-(x<y)), toString=str,
    isZero=lambda x:x==0,isOne=lambda x:x==1,isInteger=lambda x:x.denominator()==1,
    isPositive=lambda x:x>0,isNegative=lambda x:x<0,
    asIntegerRatio=lambda x:[x.numerator(),x.denominator()],
    rational_gcd=lambda x,y:x.gcd(y),rational_lcm=lambda x,y:x.lcm(y),
    # Vendored rational.pyx defaults to even; installed 10.3 defaults to away.
    round=lambda x:x.round('even'), continued_fraction=lambda x:list(x.continued_fraction()),
    simplify=lambda x:x,signAsInt=lambda x:x.sign(),padic_valuation=lambda x,p:x.valuation(p),
    denominator_valuation=lambda x,p:x.denominator().valuation(p),
    numerator_valuation=lambda x,p:x.numerator().valuation(p)))

def normalize(x):
    if isinstance(x,bool):return x
    if isinstance(x,(list,tuple)):return [normalize(v) for v in x]
    return str(x).replace('+Infinity','Infinity')

def run(method,n,d,args):
    x=QQ((n,d))
    if method in ALIASES:return ALIASES[method](x,*args)
    return getattr(x,method)(*args)

def compare(method,n,d,args):
    try:return json.dumps(dict(value=normalize(run(method,n,d,args))),separators=(',',':'))
    except Exception as e:return json.dumps(dict(error=type(e).__name__),separators=(',',':'))

FUNCTIONS={f'rational_{m}':lambda n,d,m=m:compare(m,n,d,[]) for m in UNARY}
FUNCTIONS.update({f'rational_{m}':lambda n,d,a,b,m=m:compare(m,n,d,[QQ((a,b))]) for m in BINARY})
FUNCTIONS.update({f'rational_{m}':lambda n,d,a,m=m:compare(m,n,d,[a]) for m in INTEGER+SET})

# Binary encodings compare exact IEEE-754 output, including signed underflow.
def numeric(method,n,d):
    try:
        x=QQ((n,d))
        if method == 'matrix_norm': value=struct.pack('>d',float(matrix(QQ,[[x]]).norm(2))).hex()
        elif method in ('toNumber','n','numerical_approx'):
            value=struct.pack('>d',float(x if method=='toNumber' else x.n())).hex()
        else:value=str(ZZ((x.global_height()*10**8).round()))
        return json.dumps(dict(value=value),separators=(',',':'))
    except Exception as e:return json.dumps(dict(error=type(e).__name__),separators=(',',':'))
for method in ['toNumber','n','numerical_approx','global_height','matrix_norm']:
    FUNCTIONS['rational_'+method]=lambda n,d,method=method:numeric(method,n,d)
FUNCTIONS['rational_sqrt']=lambda n,d: compare('sqrt',n,d,[]) if QQ((n,d)).is_square() else compare('sqrt_no_extend',n,d,[])
ALIASES['sqrt_no_extend']=lambda x:x.sqrt(extend=False)
FUNCTIONS['rational_gamma']=lambda n,d:compare('gamma',n,d,[])

# Remaining conveniences and factory spellings, mapped to Sage coercion.
def extra(method,n,d,arg=0):
    try:
        if method == 'fromString': value=QQ(''.join(chr(int(c)) for c in n))
        elif method in ('from','qq_from'):
            number=float(n)/float(d) if d else float('inf' if n>0 else '-inf' if n<0 else 'nan')
            value=QQ(number)
        elif method in ('constructor','fromTuple'): value=QQ((n,d))
        elif method == 'zero': value=QQ.zero()
        elif method == 'one': value=QQ.one()
        else:
            x=QQ((n,d))
            if method in ('minpoly','charpoly'): value=list(getattr(x,method)())
            elif method == 'factorial': value=ZZ(x).factorial()
            elif method == 'roundToRational':
                scale=ZZ(10)**arg
                value=(x*scale).round('even')/scale
            elif method == 'ndigits': value=len(abs(x.numerator()).str(arg))+len(x.denominator().str(arg))
            elif method == 'nbits': value=x.numerator().nbits()+x.denominator().nbits()
            elif method == 'local_height': value=ZZ((x.local_height(arg)*10**8).round())
            elif method == 'round_mode': value=x.round(['even','away','toward','up','down','odd'][int(arg)])
            elif method == 'continued_fraction_hj': value=x.continued_fraction_list('hj')
        return json.dumps(dict(value=normalize(value)),separators=(',',':'))
    except Exception as e:return json.dumps(dict(error=type(e).__name__),separators=(',',':'))
for method in ['from','qq_from','fromTuple','zero','one','constructor','minpoly','charpoly','factorial','roundToRational','ndigits','nbits','local_height','round_mode','continued_fraction_hj']:
    FUNCTIONS['rational_'+method]=lambda n,d,*args,method=method:extra(method,n,d,*args)
FUNCTIONS['rational_fromString']=lambda codes:extra('fromString',codes,1)

for method in ['numerator','numer','denominator','denom','sign']:
    FUNCTIONS['rational_'+method]=lambda n,d,method=method:compare(method,n,d,[])
FUNCTIONS['rational_from_rational']=lambda n,d:extra('constructor',n,d)
FUNCTIONS['rational_from_bigint']=lambda n,d:extra('constructor',n,1)
for method in ['qq_fromString','from_string']:
    FUNCTIONS['rational_'+method]=FUNCTIONS['rational_fromString']

def valuation_error(n,d,p,method):
    names=['valuation','ord','val_unit','numerator_valuation','denominator_valuation','padic_valuation']
    try:return json.dumps(dict(value=normalize(run(names[int(method)],n,d,[ZZ(p)]))),separators=(',',':'))
    except Exception as e:return json.dumps(dict(error=type(e).__name__,message=str(e)),separators=(',',':'))
FUNCTIONS['rational_valuation_error']=valuation_error

import itertools
from sage.rings.rational_field import RationalField
QQ_ALIASES={'getInstance':lambda:QQ,'toString':lambda:str(QQ),'some_elements':lambda:list(QQ.some_elements())}
def qq_result(run):
    try:return json.dumps(dict(value=normalize(run())),separators=(',',':'))
    except Exception as e:return json.dumps(dict(error=type(e).__name__),separators=(',',':'))
QQ_METHODS=['getInstance', 'zero', 'one', 'is_field', 'characteristic', 'is_ring', 'is_integral_domain', 'is_prime_field', 'is_absolute', 'is_finite', 'degree', 'absolute_degree', 'ngens', 'gens', 'an_element', 'some_elements', 'discriminant', 'absolute_discriminant', 'relative_discriminant', 'class_number', 'signature', 'order', 'maximal_order', 'ring_of_integers', 'number_field', 'power_basis', 'toString', '_latex_']

for method in QQ_METHODS:
    FUNCTIONS['qq_'+method]=lambda method=method:qq_result(QQ_ALIASES.get(method,lambda:getattr(QQ,method)()))
for method in ['gen','zeta','primes_of_bounded_norm_iter','range_by_height']:
    def qq_bound(a,*args,method=method):
        def run():
            x=getattr(QQ,method)(a,*args)
            return list(x) if method in ['primes_of_bounded_norm_iter','range_by_height'] else x
        return qq_result(run)
    FUNCTIONS['qq_'+method]=qq_bound
FUNCTIONS['qq_iterator']=lambda n:qq_result(lambda:list(itertools.islice(QQ,int(n))))
for method in ['selmer_generators','selmer_group_iterator']:
    FUNCTIONS['qq_'+method]=lambda S,m,method=method:qq_result(lambda:QQ.selmer_generators(S,m,orders=True) if method=='selmer_generators' else list(QQ.selmer_group_iterator(S,m)))
FUNCTIONS['qq_quadratic_defect']=lambda n,d,p:qq_result(lambda:QQ.quadratic_defect(QQ((n,d)),p))
def ieee_divide(n,d):
    return float(n)/float(d) if d else float('inf' if n>0 else '-inf' if n<0 else 'nan')
FUNCTIONS['qq_pair']=lambda n,d,a,b,mode:qq_result(lambda:QQ((QQ((n,d)),QQ((a,b)))) if mode==0 else QQ((ieee_divide(n,d),ieee_divide(a,b))))
FUNCTIONS['qq_prime_fraction_bound']=lambda n,d:qq_result(lambda:list(QQ.primes_of_bounded_norm_iter(QQ((n,d)))))
FUNCTIONS['qq___call__']=lambda n:qq_result(lambda:QQ(n))

def rational_root_constructed(base,k,offset):
    x=QQ((ZZ(base)**k+offset,ZZ(3)**k))
    exact=x.is_nth_power(k)
    return json.dumps(dict(value=[bool(exact),str(x.nth_root(k)) if exact else 'none']),separators=(',',':'))
FUNCTIONS['rational_root_constructed']=rational_root_constructed

def rational_factory_scalar(mode):
    cases=[lambda:QQ(),lambda:QQ(ZZ(7)),lambda:QQ(True),lambda:QQ(None),lambda:QQ(),lambda:QQ(None),lambda:QQ(True),lambda:QQ((None,1)),lambda:QQ((1,None)),lambda:QQ((True,2)),lambda:QQ((False,2)),lambda:QQ((1,False)),lambda:QQ([None]),lambda:QQ(False),lambda:QQ(ZZ(-3))]
    return qq_result(cases[int(mode)])
FUNCTIONS['rational_factory_scalar']=rational_factory_scalar

FUNCTIONS['qq_tuple_string']=lambda codes,den:qq_result(lambda:QQ((''.join(chr(int(c)) for c in codes),den)))


def rational_eq_float(n, d, bits):
    import struct
    value = struct.unpack('>d', struct.pack('>Q', int(bits)))[0]
    return bool(QQ(n)/d == value)
FUNCTIONS['rational_eq_float'] = rational_eq_float
