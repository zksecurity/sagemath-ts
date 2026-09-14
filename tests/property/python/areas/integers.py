"""Live integer method oracle, with no substitutions for wrong port results."""
import json
from sage.all import Integer, power_mod, Mod, RR, frac, pari, ZZ, QQ, euler_phi, moebius, number_of_divisors, sigma, bell_number, catalan_number, fibonacci, lucas_number2, number_of_partitions, primitive_root, nth_prime, prime_pi, carmichael_lambda, legendre_symbol
import operator

unary = ['abs','sign','is_prime','is_unit','isqrt','is_square','nbits',
  'prime_divisors','divisors','squarefree_part','next_prime','next_prime_power','is_prime_power',
  'is_perfect_power','is_irreducible','is_pseudoprime','is_squarefree','is_discriminant',
  'is_fundamental_discriminant','factorial','radical','popcount','multiplicative_order',
  'sqrtrem','bits','bit_length','numerator','denominator','floor','ceil','round','trunc',
  '__invert__','previous_prime', 'class_number', 'factor', 'euler_phi', 'moebius', 'number_of_divisors', 'hamming_weight', 'square', 'cube', 'content', 'primitive_part', 'continued_fraction', 'frac', 'is_power_of_two', 'is_even', 'is_odd', 'bell_number', 'catalan_number', 'fibonacci', 'lucas_number', 'number_of_partitions', 'primitive_root', 'nth_prime', 'prime_pi', 'carmichael_lambda', 'neg', 'isZero', 'valueOf', 'toString']
binary = ['gcd','lcm','xgcd','mod','quo_rem','ndigits','valuation','nth_root','exact_log',
  'prime_to_m_part','jacobi','kronecker','binomial','digits','__mod__','__floordiv__',
  'divides','inverse_mod','ord','trial_division', 'sigma', 'legendre_symbol', 'is_coprime', 'add', 'sub', 'mul', 'div', 'pow', 'eq', 'lt', 'le', 'gt', 'ge']

def normalize(value):
    if isinstance(value, bool): return value
    if isinstance(value, (list,tuple)): return [normalize(x) for x in value]
    # The port represents positive infinity by the string 'Infinity'.
    return str(value).replace('+Infinity', 'Infinity')

# Port-only conveniences map to Sage's free functions/operators rather than
# pretending Sage Integer has methods with these names (see DEVIATIONS.md).
aliases = {m: globals()[m] for m in ['euler_phi','moebius','number_of_divisors','sigma',
    'bell_number','catalan_number','fibonacci','number_of_partitions','primitive_root',
    'nth_prime','prime_pi','carmichael_lambda','legendre_symbol']}
aliases.update({m: getattr(operator,op) for m,op in dict(add='add',sub='sub',mul='mul',
    div='floordiv',pow='pow',eq='eq',lt='lt',le='le',gt='gt',ge='ge',neg='neg').items()})
# Installed primecountpy crashes even on small inputs on this host. PARI's
# independent exact primepi supplies the same counting function for these cases.
aliases['prime_pi'] = lambda n: pari(n).primepi().sage()
aliases.update(dict(hamming_weight=lambda n:n.popcount(), square=lambda n:n**2,
    cube=lambda n:n**3, content=abs, primitive_part=lambda n:n.sign(),
    continued_fraction=lambda n:list(QQ(n).continued_fraction()), frac=frac,
    is_coprime=lambda n,a:n.gcd(a)==1, is_power_of_two=lambda n:n>0 and n.prime_to_m_part(2)==1,
    is_even=lambda n:n%2==0, is_odd=lambda n:n%2==1, lucas_number=lambda n:lucas_number2(n,1,-1),
    isZero=lambda n:n==0, valueOf=lambda n:n, toString=lambda n:str(n)))

def evaluate(method,n,args):
    if method=='factor':
        f=n.factor()
        return ([[-1,1]] if f.unit()==-1 else []) + list(f)
    if method in aliases: return aliases[method](n,*args)
    return getattr(n,method)(*args)

def compare(method,n,args):
    try:
        # VENDORED 10.9.beta4: categories/unique_factorization_domains.py:283-284
        # returns self for zero; Sage 10.3 predates that zero guard.
        value = ZZ(0) if method == 'radical' and n == 0 else evaluate(method,ZZ(n),args)
        return json.dumps(dict(value=normalize(value)),separators=(',',':'))
    except Exception as error: return json.dumps(dict(error=type(error).__name__),separators=(',',':'))
FUNCTIONS = {f'integer_{m}': (lambda n, m=m: compare(m,n,[])) for m in unary}
FUNCTIONS.update({f'integer_{m}': (lambda n,a,m=m: compare(m,n,[a])) for m in binary})


# The following methods are port conveniences. Root certificates compare with
# Sage's root existence and verify the returned root; the chosen nth root is
# intentionally not canonical (DEVIATIONS.md: Port-only conveniences).
def special(method, n, *args):
    try:
        if method == 'digit_sum': value = sum(abs(d) for d in ZZ(n).digits(args[0]))
        elif method == 'sqrt_mod':
            roots = Mod(n,args[0]).sqrt(all=True)
            value = min(ZZ(r) for r in roots) if roots else 'None'
        elif method == 'nth_root_mod':
            value = bool(Mod(n,args[1]).nth_root(args[0],all=True))
        elif method == 'is_primitive_root':
            x = Mod(n,args[0])
            value = bool(x.is_unit() and x.multiplicative_order() == euler_phi(args[0]))
        elif method == 'powermod':
            # integer.pyx:3656 calls GMP mpz_powm without checking inversion;
            # a negative exponent on a non-unit aborts installed Sage. The
            # port deliberately uses Sage's safe arith.power_mod error path.
            if args[0] < 0 and ZZ(n).gcd(args[1]) != 1:
                value = power_mod(n,*args)
            else: value = ZZ(n).powermod(*args)
        elif method == 'log': value = '-Infinity' if n == 0 else ZZ(RR(n).log().floor())
        elif method == 'real_log': value = ZZ((RR(n).log()*10**8).round())
        elif method == 'global_height': value = ZZ((ZZ(n).global_height()*10**8).round())
        elif method == 'is_quadratic_residue': value = bool(Mod(n,args[0]).is_square())
        elif method == 'bit': value = ZZ(n).test_bit(args[0])
        elif method == 'core':
            # Port-only t-th-power-free part; a Sage factorization certificate.
            value = ZZ(n).sign()
            if n != 0:
                for p,e in abs(ZZ(n)).factor(): value *= p**(e % args[0])
        elif method == 'is_strong_pseudoprime':
            from gmpy2 import is_strong_prp
            value = bool(is_strong_prp(int(n),int(args[0])))
        return json.dumps(dict(value=normalize(value)),separators=(',',':'))
    except Exception as error: return json.dumps(dict(error=type(error).__name__),separators=(',',':'))
for name in ['digit_sum','sqrt_mod','nth_root_mod','is_primitive_root','powermod','log',
             'real_log','global_height','is_quadratic_residue','bit','core','is_strong_pseudoprime']:
    FUNCTIONS['integer_'+name] = lambda n,*args,name=name:special(name,n,*args)

FUNCTIONS['integer_root_constructed']=lambda base,k,offset,sign:compare('nth_root',ZZ(sign)*(ZZ(base)**k+offset),[k,True])
FUNCTIONS['integer_root_truncated']=lambda n,k:compare('nth_root',n,[k,True])

def integer_conversion(value,base,wrapper):
    try:
        result=(Integer if wrapper else ZZ)(value,base=base)
        return json.dumps(dict(value=str(result)),separators=(',',':'))
    except Exception as error:return json.dumps(dict(error=type(error).__name__,message=str(error)),separators=(',',':'))
FUNCTIONS['zz_string']=lambda codes,base,wrapper:integer_conversion(''.join(chr(int(c)) for c in codes),base,wrapper)
FUNCTIONS['zz_number']=lambda n,d,wrapper:integer_conversion(float(n)/float(d) if d else float('inf' if n>0 else '-inf' if n<0 else 'nan'),0,wrapper)
FUNCTIONS['zz_rational']=lambda n,d,wrapper:integer_conversion(QQ((n,d)),0,wrapper)
FUNCTIONS['zz_wrapped_integer']=lambda n,wrapper:integer_conversion(ZZ(n),0,wrapper)
FUNCTIONS['zz_boolean']=lambda n,wrapper:integer_conversion(bool(n),0,wrapper)
FUNCTIONS['zz_empty']=lambda mode,wrapper:integer_conversion(None,0,wrapper)
FUNCTIONS['zz_digits']=lambda digits,base,wrapper:integer_conversion(list(digits),base,wrapper)
for method in ['getInstance','zero','one','characteristic','is_field','is_ring','is_integral_domain','toString']:
    FUNCTIONS['zz_'+method]=lambda method=method:json.dumps(dict(value=normalize(ZZ if method=='getInstance' else str(ZZ) if method=='toString' else getattr(ZZ,method)())),separators=(',',':'))

def zz_rational_no_base(n,d,wrapper):
    try:return json.dumps(dict(value=str((Integer if wrapper else ZZ)(QQ((n,d))))),separators=(',',':'))
    except Exception as error:return json.dumps(dict(error=type(error).__name__,message=str(error)),separators=(',',':'))
FUNCTIONS['zz_rational_no_base']=zz_rational_no_base

FUNCTIONS['zz_large_string']=lambda n,base,wrapper:integer_conversion(ZZ(n).str(base),base,wrapper)
def zz_modular(n,m,wrapper):
    try:return json.dumps(dict(value=str((Integer if wrapper else ZZ)(Mod(n,m)))),separators=(',',':'))
    except Exception as error:return json.dumps(dict(error=type(error).__name__,message=str(error)),separators=(',',':'))
FUNCTIONS['zz_modular']=zz_modular
def zz_hook(n,wrapper):
    class Hook:
        def _integer_(self,parent):
            assert parent is ZZ
            return ZZ(n)
    return json.dumps(dict(value=str((Integer if wrapper else ZZ)(Hook()))),separators=(',',':'))
FUNCTIONS['zz_hook']=zz_hook

from integer_rational import integer_rational
FUNCTIONS["integer_rational"] = integer_rational
