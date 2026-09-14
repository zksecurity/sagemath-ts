"""SageMath side of the ``arith_special`` property-test area.

Cases: tests/property/cases/arith_special.cases.json
"""

from sage.all import *
from sage.arith.misc import dedekind_psi, subfactorial, primes, differences, primes_first_n, dedekind_sum, integer_floor, integer_ceil, integer_trunc, algebraic_dependency


def arith_binomial(n, k):
    """Compute binomial coefficient C(n,k)."""
    return binomial(n, k)


def arith_fibonacci(n):
    """Compute the n-th Fibonacci number."""
    return fibonacci(n)


def arith_lucas_number(n):
    """Compute the n-th Lucas number using lucas_number2."""
    return lucas_number2(n, 1, -1)


def arith_factorial(n):
    """Compute n!."""
    return factorial(n)


def arith_bernoulli_numerator(n):
    """Return the numerator of the n-th Bernoulli number."""
    return bernoulli(n).numerator()


def arith_bernoulli_denominator(n):
    """Return the denominator of the n-th Bernoulli number."""
    return bernoulli(n).denominator()


def arith_multinomial(*args):
    """Compute multinomial coefficient."""
    return multinomial(*args)


def arith_primitive_root(n):
    """Find a primitive root modulo n."""
    return primitive_root(n)


def arith_nth_prime(n):
    """Return the n-th prime number."""
    # Bundled 10.9 renamed "non-positive" to "nonpositive"; installed 10.3 differs only in text.
    if n <= 0:
        raise ValueError("nth prime meaningless for nonpositive n (=%s)" % n)
    return nth_prime(n)


def arith_subfactorial(n):
    """Compute !n (number of derangements)."""
    return subfactorial(n)


def arith_carmichael_lambda(n):
    """Compute Carmichael's lambda function."""
    return carmichael_lambda(n)


def arith_dedekind_psi(n):
    """Compute Dedekind's psi function."""
    # Use the imported dedekind_psi from sage.arith.misc
    return int(dedekind_psi(n))


def arith_primes_iterator(mode, start, stop, count, proof, wrapped):
    p = None if proof == 0 else proof == 1
    iterator = primes() if mode == 0 else primes(start) if mode == 1 else primes(start, infinity if mode == 3 else stop, p)
    values, done = [], False
    for _ in range(count):
        try:
            values.append(next(iterator))
        except StopIteration:
            done = True
            break
    return ','.join(map(str, values)) + '|' + str(int(done))


def arith_differences(length, degree, offset, order, wrapped):
    return differences([(Integer(i) + offset)**degree for i in range(length)], order)


def arith_dedekind_sum(h, k, mode, wrapped):
    # Bundled 10.9 calls cypari2 directly; installed 10.3 uses the GP subprocess
    # interface, whose transport-generated TypeError is not the original backend error.
    if mode == 3:
        from sage.libs.pari import pari as native_pari
        value = QQ(native_pari.sumdedekind(h, k))
    elif mode == 0:
        value = dedekind_sum(h, k)
    else:
        value = dedekind_sum(h, k, algorithm=['default', 'default', 'flint', 'pari', 'invalid', None][mode])
    return [value.numerator(), value.denominator()]


ROUNDING_TEXTS = ['', 'foo', '-2.5', '+2.5', '1_000.5', ' nan ', '-Infinity', '0e1000', '1e1000', '-1e-1000', '１２.５', '١.٥']
class RoundingMethodError:
    # Protocol adapter for a port Integer with overridden rounding methods.
    # A Sage Integer subclass constructor does not retain these Python overrides.
    def __init__(self, value, error):
        self.value, self.error = value, error
    def floor(self): raise self.error('rounding sentinel')
    def ceil(self): raise self.error('rounding sentinel')
    def __float__(self): return float(self.value)
    def __ge__(self, other): return self.value >= other

def arith_integer_rounding(kind, a, b, precision, operation):
    import struct
    if kind == 0 or kind == 1: value = ZZ(a)
    elif kind == 2: value = QQ(a)/b
    elif kind == 3: value = float(a)/float(b)
    elif kind == 4: value = RealField(int(precision))(str(a)+'e'+str(b))
    elif kind == 5: value = RealField(int(precision))(['NaN','inf','-inf','-0'][a])
    elif kind == 6: value = struct.unpack('>d',struct.pack('>Q',int(a)))[0]
    elif kind == 7: value = ROUNDING_TEXTS[a]
    elif kind == 8: value = [True,False,None,b'2.5',b'x',b'nan'][a]
    elif kind == 12: value = RoundingMethodError(a, AttributeError)
    elif kind == 13: value = RoundingMethodError(a, ValueError)
    elif kind == 9: value = Mod(a,7)
    else: value = GF(7 if kind == 10 else 2)(a)
    return [integer_floor,integer_ceil,integer_trunc][operation](value)


def arith_algdep_dispatch(bits, degree, mode, bound, first, second):
    import struct
    z = struct.unpack('>d', struct.pack('>Q', int(bits)))[0]
    options = {} if mode == 0 else {'height_bound': bound} if mode == 1 else \
        {'height_bound': bound, 'proof': True} if mode == 2 else {'proof': True} if mode == 3 else \
        {'known_bits': first} if mode == 4 else {'use_bits': first} if mode == 5 else \
        {'known_digits': first} if mode == 6 else {'use_digits': first} if mode == 7 else \
        {'known_bits': first, 'use_bits': second, 'known_digits': first, 'use_digits': second}
    result = algebraic_dependency(z, degree, **options)
    return None if result is None else list(result)


FUNCTIONS = {
    'algdep_dispatch': arith_algdep_dispatch,
    'integer_rounding': arith_integer_rounding,
    'dedekind_native': lambda h, k, backend: arith_dedekind_sum(h, k, 2 if backend == 0 else 3, 0),
    'dedekind_sum': arith_dedekind_sum,
    'prime_native': lambda n: pari.prime(n).sage(),
    'first_primes_nonfinite': lambda mode: primes_first_n(float('nan') if mode == 0 else float('inf') if mode == 1 else -float('inf')),
    'primes_iterator': arith_primes_iterator,
    'differences': arith_differences,
    'binomial': arith_binomial,
    'fibonacci': arith_fibonacci,
    'lucas_number': arith_lucas_number,
    'factorial': arith_factorial,
    'bernoulli_numerator': arith_bernoulli_numerator,
    'bernoulli_denominator': arith_bernoulli_denominator,
    'multinomial': arith_multinomial,
    'primitive_root': arith_primitive_root,
    'nth_prime': arith_nth_prime,
    'nth_prime_wrapped': lambda n, wrapped: arith_nth_prime(n),
    'primes_first_n': lambda n, mode: primes_first_n(float(n)/1000 if mode == 3 else float(n) if mode == 0 else n),
    'subfactorial': arith_subfactorial,
    'carmichael_lambda': arith_carmichael_lambda,
    'dedekind_psi': arith_dedekind_psi,
}

ALGDEP_TEXTS = ['0','-0','1','-1','1.5','-1.5','1.4142135623730950488016887242096980785696718753769','3.1415926535897932384626433832795028841971693993751','1e-30','1e30','NaN','inf','-inf','1.6180339887498948482045868343656381177203091798058']
def algdep_real(kind, value, denominator, precision, degree, mode, bound, first, second):
    z = ZZ(value) if kind < 2 else QQ(value)/denominator if kind == 2 else RealField(int(precision))(ALGDEP_TEXTS[value])
    options = {} if mode == 0 else {'height_bound': bound} if mode == 1 else {'height_bound': bound, 'proof': True} if mode == 2 else {'proof': True} if mode == 3 else {'known_bits': first} if mode == 4 else {'use_bits': first} if mode == 5 else {'known_digits': first} if mode == 6 else {'use_digits': first} if mode == 7 else {'known_bits': first, 'use_bits': second, 'known_digits': first, 'use_digits': second}
    result = z.algebraic_dependency(degree) if mode == 9 else z.algdep(degree) if mode == 10 else algebraic_dependency(z,degree,**options)
    return None if result is None else list(result)
FUNCTIONS['algdep_real'] = algdep_real

def native_algdep(bits,degree):
    # Select the bundled PARI version, whose LLL strategy differs from 10.3's.
    import json
    from pari_lll_dependents import pari_lll_dependents
    return [ZZ(c) for c in json.loads(pari_lll_dependents(5,degree,0,[bits]))]
FUNCTIONS['native_algdep'] = native_algdep

def prime_traversal(value,operation,wrapped,data):
    functions=[next_prime_power,previous_prime_power,next_probable_prime,is_power_of_two,is_pseudoprime_power,is_prime_power,next_prime,previous_prime]
    if operation in (4,5):
        result=functions[operation](value,get_data=bool(data))
        return list(result) if data else result
    return functions[operation](value)
def prime_power_range(start,stop,mode,wrapped):
    from sage.arith.misc import prime_powers,eratosthenes
    return prime_powers(start) if mode==0 else prime_powers(start,stop) if mode==1 else eratosthenes(start)
def native_prime_traversal(value,operation):
    from sage.libs.pari import pari as native_pari
    if operation==0:return ZZ(native_pari(value).nextprime())
    if operation==1:return ZZ(native_pari(value).precprime())
    answer=ZZ(value).is_prime_power(proof=False,get_data=True)
    return list(answer) if answer[1] else None
FUNCTIONS['prime_traversal']=prime_traversal
FUNCTIONS['prime_power_range']=prime_power_range
FUNCTIONS['native_prime_traversal']=native_prime_traversal

def hilbert_dispatch(an,ad,bn,bd,p,mode,wrapped):
    from sage.arith.misc import hilbert_symbol
    a=QQ(an)/ad;b=QQ(bn)/bd
    if mode==0:return hilbert_symbol(a,b,p)
    return hilbert_symbol(a,b,p,algorithm=['pari','direct','all','invalid',None,''][mode-1])
def hilbert_conductor_dispatch(a,b,operation,wrapped):
    from sage.arith.misc import hilbert_conductor,hilbert_conductor_inverse
    return list(hilbert_conductor_inverse(a)) if operation else hilbert_conductor(a,b)
def native_hilbert(a,b,p):
    from sage.libs.pari import pari as native_pari
    return ZZ(native_pari(a).hilbert(b,p))
FUNCTIONS['hilbert_dispatch']=hilbert_dispatch
FUNCTIONS['hilbert_conductor_dispatch']=hilbert_conductor_dispatch
FUNCTIONS['native_hilbert']=native_hilbert

def native_valuation_unit(n,p):
    from sage.libs.pari import pari as native_pari
    v=ZZ(native_pari(n).valuation(p))
    return [v,n//p**v]
FUNCTIONS['native_valuation_unit']=native_valuation_unit

# Free integer-valued APIs map Sage int/Integer to JavaScript bigint.
# The frame checks that documented return representation as well as the value.
def arith_scalar_edge(n,m,operation,wrapped):
    from sage.arith.misc import quadratic_residues,fundamental_discriminant,odd_part,prime_to_m_part
    import json
    if operation==0:return json.dumps([str(x) for x in quadratic_residues(n)],separators=(',',':'))
    result=fundamental_discriminant(n) if operation==1 else odd_part(n) if operation==2 else prime_to_m_part(n,m)
    return json.dumps([str(result),'integer'],separators=(',',':'))
def arith_continuant(length,power,offset,order,mode,wrapped):
    from sage.arith.misc import continuant
    import json
    values=[ZZ(i)**power+offset for i in range(length)]
    result=continuant(values) if mode==0 else continuant(values,order)
    return json.dumps([str(result),'integer'],separators=(',',':'))
def arith_squarefree_prefix(n,take,wrapped):
    from sage.arith.misc import squarefree_divisors
    import json
    iterator=squarefree_divisors(n);values=[];done=False
    for _ in range(take):
        try:values.append(str(next(iterator)))
        except StopIteration:done=True;break
    return json.dumps([values,done],separators=(',',':'))
FUNCTIONS['arith_scalar_edge']=arith_scalar_edge
FUNCTIONS['arith_continuant']=arith_continuant
FUNCTIONS['arith_squarefree_prefix']=arith_squarefree_prefix

def arithmetic_factory(order,a,b,operation,wrapped):
    from sage.arith.misc import get_gcd,get_inverse_mod
    try:
        return [get_gcd,get_inverse_mod][operation](order)(a,b)
    except SystemError:
        # Cython's -1 sentinel: CPython may include an unstable object address,
        # or emit the shorter diagnostic after specializing the call site.
        raise SystemError('error return without exception set') from None
def arithmetic_factory_selection(order,operation,wrapped):
    from sage.arith.misc import get_gcd,get_inverse_mod
    maker=[get_gcd,get_inverse_mod][operation]
    f=maker(order);g=maker(order)
    return f.__name__+'|'+str(int(f is g))
FUNCTIONS['arithmetic_factory']=arithmetic_factory
FUNCTIONS['arithmetic_factory_selection']=arithmetic_factory_selection

# Sage 10.3 predates CRT_basis's non-coprime option and CRT_list's singleton
# change. Execute the bundled 10.9 bodies unchanged with their original globals.
def _bundled_crt():
    import ast
    from pathlib import Path
    import sage.arith.misc as original
    if not hasattr(_bundled_crt,'functions'):
        source=Path(__file__).resolve().parents[4]/'reference/sage/src/sage/arith/misc.py'
        names={'crt','CRT_list','CRT_basis','CRT_vectors','smooth_part','coprime_part','two_squares','three_squares','four_squares','sum_of_k_squares'}
        nodes=[n for n in ast.parse(source.read_text()).body if isinstance(n,ast.FunctionDef) and n.name in names]
        namespace=dict(vars(original))
        exec(compile(ast.Module(body=nodes,type_ignores=[]),str(source),'exec'),namespace)
        namespace['CRT']=namespace['crt']
        _bundled_crt.functions=namespace
    return _bundled_crt.functions

_crt_bases=[[],[3],[3,5],[60,90,150],[7,6,10],[2,4,8,16],[-3,5],[-3,-5],[0],[0,3],[1,-1,3],[2,3,5,7,11],[4,6,9,10]]
def crt_scalar_dispatch(a,b,m,n,wrapped):
    return _bundled_crt()['crt'](a,b,m,n)
def crt_list_dispatch(code,offset,delta,mode,wrapped):
    funcs=_bundled_crt();moduli=list(map(ZZ,_crt_bases[code]))
    values=[offset+i*(i+1) for i in range(max(0,len(moduli)+delta))]
    if mode in (3,7,8):
        values=[Mod(v,m) if mode==3 else GF(m)(v) for v,m in zip(values,moduli)]
        result=funcs['CRT_list'](values)
        return [result.lift(),result.modulus(),len(values)==1 and result is values[0]]
    if mode==4:return funcs['CRT_list']('bad',moduli)
    if mode==5:return funcs['CRT_list'](values,'bad')
    if mode==6:return funcs['CRT_list'](values)
    return funcs[['CRT_list','crt','CRT'][mode]](values,moduli)
def crt_basis_dispatch(length,encoding,required,wrapped):
    moduli=[]
    for _ in range(length):
        moduli.append(ZZ(encoding%7)-3);encoding//=7
    if length<0:moduli=list(map(ZZ,_crt_bases[encoding]))
    result=_bundled_crt()['CRT_basis'](moduli,require_coprime_moduli=bool(required))
    # Avoid the generic transcript formatter's two-entry factorization heuristic.
    import json
    def frame(x):return [frame(v) for v in x] if isinstance(x,(list,tuple)) else x if isinstance(x,bool) else str(x)
    return json.dumps(frame(result),separators=(',',':'))

def crt_vector_dispatch(code,offset,columns,shape,wrapped):
    moduli=list(map(ZZ,_crt_bases[code]))
    values=[[offset+i*j+j for j in range(columns)] for i in range(len(moduli))]
    if shape==1:values.append([offset]*columns)
    if shape==2 and values:values[-1]=values[-1][:-1]
    if shape==3 and values:values[-1].append(offset)
    return _bundled_crt()['CRT_vectors'](values,moduli)
FUNCTIONS['crt_scalar_dispatch']=crt_scalar_dispatch
FUNCTIONS['crt_list_dispatch']=crt_list_dispatch
FUNCTIONS['crt_basis_dispatch']=crt_basis_dispatch
FUNCTIONS['crt_vector_dispatch']=crt_vector_dispatch

# ProductTree/prod_with_derivative's installed bodies equal the bundled source
# after stripping docstrings (verified structurally with ast).
_product_bases=[[], [2, 3], [3, 2], [6], [6, 14], [-6, 14], [6, -14], [2, 2], [0], [1], [-1], [-2, 3], [5, 7, 11, 13], [4, 9, 25], [2, 3, 5, 7, 11], [2, 3, 5, 7, 11, 13, 17], [2, 3, 5, 7, 11, 13, 17, 19], [2, 3, 5, 7, 11, 13, 17, 19, 23], [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47], [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53], [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59], [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127], [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131], [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137], [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307], [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311], [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313], [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709], [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719], [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727]]
def _product_frame(result):
    import json
    def frame(x):return [frame(v) for v in x] if isinstance(x,(list,tuple)) else x if isinstance(x,(bool,str)) or x is None else str(x)
    return json.dumps(frame(result),separators=(',',':'))
def product_tree_dispatch(code,x,operation,delta,wrapped):
    from sage.rings.generic import ProductTree
    leaves=list(map(ZZ,_product_bases[code]));tree=ProductTree(iter(leaves) if wrapped&2 else leaves)
    if operation==0:return _product_frame([len(tree),list(tree),tree.leaves(),tree.layers,tree.leaves() is tree.layers[0]])
    if operation==1:return tree.root()
    if operation==2:return tree.remainders(x)
    xs=[x+i*i for i in range(max(0,len(leaves)+delta))]
    if operation==3:
        first=tree.interpolation(xs);cache=tree._crt_bases
        return _product_frame([first,tree.interpolation(xs),cache is tree._crt_bases])
    try:tree.interpolation(xs)
    except Exception as e:error=[type(e).__name__,str(e)]
    else:error=None
    try:second=tree.interpolation([x]*len(leaves));answer=[str(second)]
    except Exception as e:answer=[type(e).__name__,str(e)]
    return _product_frame([error,tree._crt_bases,answer])
def smooth_part_dispatch(x,code,operation,mode):
    from sage.rings.generic import ProductTree
    leaves=list(map(ZZ,_product_bases[code]));base=ProductTree(leaves) if mode&1 else iter(leaves) if mode&4 else leaves
    result=_bundled_crt()[['smooth_part','coprime_part'][operation]](x,base)
    return _product_frame(list(result)) if operation==0 else result
def smooth_nontermination(x,code,operation):
    import signal
    def timed_out(signum,frame):raise TimeoutError()
    previous=signal.signal(signal.SIGALRM,timed_out)
    try:
        signal.setitimer(signal.ITIMER_REAL,0.1)
        return smooth_part_dispatch(x,code,operation,0)
    except TimeoutError:
        # Only fixed source-proven nonterminating inputs use this adapter.
        raise NotImplementedError('SAGE_NOT_IMPLEMENTED: smooth_part: original does not terminate for this factor base') from None
    finally:
        signal.setitimer(signal.ITIMER_REAL,0)
        signal.signal(signal.SIGALRM,previous)
def product_derivative_dispatch(length,x,wrapped):
    from sage.rings.generic import prod_with_derivative
    pairs=[(x+i,ZZ(-1)**i*i) for i in range(length)]
    return list(prod_with_derivative(iter(pairs) if wrapped&2 else pairs))
FUNCTIONS['product_tree_dispatch']=product_tree_dispatch
FUNCTIONS['smooth_part_dispatch']=smooth_part_dispatch
FUNCTIONS['smooth_nontermination']=smooth_nontermination
FUNCTIONS['product_derivative_dispatch']=product_derivative_dispatch

def xlcm_dispatch(m,n,wrapped):
    from sage.arith.misc import xlcm
    return list(xlcm(m,n))
def factor_scalar_dispatch(n,operation,wrapped):
    from sage.arith.misc import dedekind_psi,carmichael_lambda,euler_phi,radical,number_of_divisors,prime_factors
    # Same bundled zero guard as the integer-method oracle; installed 10.3
    # predates categories/unique_factorization_domains.py's radical(0) behavior.
    result=ZZ(0) if operation==3 and n==0 else [dedekind_psi,carmichael_lambda,euler_phi,radical,number_of_divisors,prime_factors][operation](n)
    return _product_frame(result)
FUNCTIONS['xlcm_dispatch']=xlcm_dispatch
FUNCTIONS['factor_scalar_dispatch']=factor_scalar_dispatch

def native_factor_arithmetic(n,operation):
    from sage.libs.pari import pari as native_pari
    value=native_pari(n)
    if operation==0:return ZZ(value.eulerphi())
    if operation==1:return ZZ(value.numdiv())
    factors=value.factor()
    return _product_frame([[ZZ(factors[i,0]),ZZ(factors[i,1])] for i in range(factors.nrows())])
def native_factor_power(p,exponent,unit,operation):
    return native_factor_arithmetic(unit*p**exponent,operation)
def native_quadratic_roots(p,exponent,valuation,shift,sign):
    from sage.libs.pari import pari as native_pari
    if not hasattr(native_quadratic_roots,'function'):
        # Bind the actual installed PARI symbol. Separate calls prevent GP from
        # declaring the new function name as a polynomial variable before install.
        native_pari('install("Zn_quad_roots","GGG","audit_native_quadroots")')
        native_quadratic_roots.function=native_pari('audit_native_quadroots')
    # Construct solvable equations (x+shift)^2=p^valuation, valuation even.
    # The raw C entry point returns NULL for unsolvable equations, which is not
    # a GP value; those domains are covered through the existing qfbsolve oracle.
    result=native_quadratic_roots.function(sign*p**exponent,2*shift,shift**2-p**valuation)
    return _product_frame([ZZ(result[0]),[ZZ(x) for x in result[1]]])
FUNCTIONS['native_factor_arithmetic']=native_factor_arithmetic
FUNCTIONS['native_factor_power']=native_factor_power
FUNCTIONS['native_quadratic_roots']=native_quadratic_roots

def square_decomposition(op,n,k,mode):
    from sage.rings import sum_of_squares as native
    functions=_bundled_crt()
    if op==3:
        if mode//2==3:k=float(k)+(0.5 if k>=0 else -0.5)
        elif mode//2==4:k=[float('nan'),float('inf'),float('-inf')][k]
        elif mode//2==0:k=int(k)
        return list(functions['sum_of_k_squares'](k,n))
    if op<3:return list(functions[['two_squares','three_squares','four_squares'][op]](n))
    result=getattr(native,['two_squares_pyx','three_squares_pyx','four_squares_pyx','is_sum_of_two_squares_pyx'][op-4])(n)
    return result if isinstance(result,bool) else list(result)
FUNCTIONS['square_decomposition']=square_decomposition

FUNCTIONS['square_decomposition_power']=lambda op,n,e,mode: square_decomposition(op,n*4**e,0,mode)

def mqrr_dispatch(u,m,t,wrapped):
    result=mqrr_rational_reconstruction(u,m,t)
    return None if result is None else list(result)
FUNCTIONS['mqrr_dispatch']=mqrr_dispatch

import json,struct
from sage.arith.misc import sort_complex_numbers_for_display

def complex_display_order(base,step,imag,mode,order,size):
 base,step,imag,mode,order,size=map(int,[base,step,imag,mode,order,size])
 ims=[[0.,-0.,1.,-1.,1e-12,-1e-12,2.],[1.,-1.,2.,-2.,3.,-3.,4.],[1e-12,-1e-12,0.,1.,-1.,2.,-2.],[1e-10,-1e-10,1e-11,-1e-11,1.,0.,-0.],[1e-300,-1e-300,0.,1e300,-1e300,1.,-1.]][imag]
 vals=[]
 for i in range(size):
  r=struct.unpack('>d',struct.pack('>Q',(base+((i%7)-3)*step)%(2**64)))[0]
  z=CDF(r,ims[i%7]);vals.append((z,i) if mode else z)
 if order==1:vals.reverse()
 elif order==2:vals=vals[::2]+vals[1::2]
 before=vals.copy();out=sort_complex_numbers_for_display(vals)
 return json.dumps([[str(next(i for i,x in enumerate(before) if x is z)) for z in out],out is vals,all(x is y for x,y in zip(vals,before))],separators=(",",":"))

def python_sort_dispatch(seed,n,pattern):
 MASK=2**64-1
 seed,n,pattern=map(int,[seed,n,pattern]);values=[];x=seed
 for i in range(n):
  x=(1664525*x+1013904223)%(2**32)
  if pattern==0:v=x%101-50
  elif pattern==1:v=i+x%3
  elif pattern==2:v=n-i-x%3
  elif pattern==3:v=float('nan') if x%11==0 else x%101-50
  elif pattern==4:v=float('nan') if i%4==0 else x%101-50
  elif pattern==5:v=(i//64)*128+63-i%64
  else:v=(i%4)*n+i//4
  values.append(v)
 trace=14695981039346656037
 class Entry:
  def __init__(self,index):self.index=index
  def __lt__(self,other):
   nonlocal trace
   trace=(((trace^self.index)*1099511628211)^other.index)&MASK
   return values[self.index]<values[other.index]
 arr=[Entry(i) for i in range(n)];out=sorted(arr)
 return json.dumps([[str(x.index) for x in out],str(trace)],separators=(",",":"))

FUNCTIONS['complex_display_order']=complex_display_order
FUNCTIONS['python_sort_dispatch']=python_sort_dispatch

def valuation_dispatch(op,n,d,p,mode):
 if op<3:
  z=ZZ(n) if mode%2 else int(n)
  if op==0:result=valuation(z,p)
  else:result=getattr(ZZ(n),['valuation','ord'][op-1])(p)
 else:
  z=QQ(n)/QQ(d)
  if op==3:result=valuation(z,p)
  elif op in [4,5]:result=getattr(z,['valuation','ord'][op-4])(p)
  elif op==6:result=z.val_unit(p)
  elif op==7:result=z.numerator().valuation(p)
  elif op==8:result=z.denominator().valuation(p)
  else:result=z.valuation(p)
 def scalar(x):return 'Infinity' if x==Infinity else str(x)
 return json.dumps(dict(value=[scalar(x) for x in result] if isinstance(result,tuple) else scalar(result)),separators=(",",":"))

def valuation_protocol(n,p,mode):
    calls=[];saved=[n,p]
    class V:
        def _integer_(self,parent):calls.append('integer');return n
        def valuation(self,p):
            calls.append('valuation')
            if mode==1:raise AttributeError('inside')
            if mode==2:raise ValueError('valuation failed')
            return saved if mode==4 else n+p
    class Getter(V):
        @property
        def valuation(self):calls.append('get');raise AttributeError('missing')
    class Coercible:
        def _integer_(self,parent):
            calls.append('integer')
            if mode==6:raise AttributeError('integer failed')
            if mode==7:raise ValueError('integer failed')
            return n
    obj=Coercible() if mode>=5 else Getter() if mode==3 else V()
    try:
        value=valuation(obj,p)
        result=dict(value=[str(x) for x in value] if isinstance(value,list) else ('Infinity' if value==Infinity else str(value)),same=value is saved,calls=calls)
    except Exception as e:result=dict(error=type(e).__name__,message=str(e),calls=calls)
    return json.dumps(result,separators=(',',':'))
FUNCTIONS['valuation_protocol']=valuation_protocol

def gmp_remove(n,p):
    import gmpy2
    unit,exponent=gmpy2.remove(n,p)
    return json.dumps([str(exponent),str(unit)],separators=(',',':'))
FUNCTIONS['gmp_remove']=gmp_remove
FUNCTIONS['valuation_dispatch']=valuation_dispatch
FUNCTIONS['valuation_power']=lambda op,p,e,unit,d: gmp_remove(unit*p**e,p) if op==10 else valuation_dispatch(op,unit*p**e,d,p,3)
FUNCTIONS['integer_bit_size']=lambda bits,delta,sign,op: getattr(ZZ(sign*((ZZ(1)<<bits)+delta)),['nbits','bit_length'][op])()

def gmp_factorial(op,n):
    """Native 64-bit GMP kernels; compare full values by length and SHA-256 bytes."""
    import ctypes, hashlib, gmpy2
    n = int(n)
    if op == 3:
        # Internal symbol is exported by the GMP 6.3.0 already loaded into Sage.
        lib = ctypes.CDLL(None)
        native = lib.__gmp_primesieve
        native.argtypes = [ctypes.POINTER(ctypes.c_ulong),ctypes.c_ulong]
        native.restype = ctypes.c_ulong
        assert ctypes.sizeof(ctypes.c_ulong) == 8
        words = ((n-5)|1)//3//64+1
        buf = (ctypes.c_ulong*words)()
        count = native(buf,n)
        payload = b''.join(int(word).to_bytes(8,'little') for word in buf)
        return json.dumps([str(count),hashlib.sha256(payload).hexdigest()],separators=(',',':'))
    if op >= 5:
        from sage.arith.misc import factorial as original_factorial
        selector = op-5 if op >= 10 else op
        arg = ZZ(n) if op >= 10 else n
        value = int(original_factorial(arg) if selector == 5 else original_factorial(arg, {6:'gmp',7:'unknown',8:'pari',9:None}[selector]))
    elif op == 0:
        try:
            value = int(ZZ(n).factorial())
        except ValueError as error:
            # Bundled integer.pyx only changes the installed 10.3 error spelling.
            if str(error) == 'factorial only defined for non-negative integers':
                raise ValueError('factorial only defined for nonnegative integers') from None
            raise
    else:
        value = int(gmpy2.double_fac(n if n%2 else n-1)) if op == 4 else int(gmpy2.fac(n))
        if op == 2:
            value >>= (value & -value).bit_length()-1
    bits = value.bit_length()
    payload = value.to_bytes((bits+7)//8,'big')
    return json.dumps([str(bits),hashlib.sha256(payload).hexdigest()],separators=(',',':'))
FUNCTIONS['gmp_factorial'] = gmp_factorial

def pari_real_add(op,px,py,ex,ey,mx,my,sx,sy):
    """Native PARI add.c, with controlled zero padding after real operands.

    addrr_sign reads x[lx] past the logical mantissa at whole-word gaps when
    extending. Zero padding makes that upstream undefined read reproducible;
    see DEVIATIONS.md: PARI real-addition padding. No result is normalized.
    """
    import ctypes
    op,px,py,ex,ey,mx,my,sx,sy=map(int,[op,px,py,ex,ey,mx,my,sx,sy])
    if not hasattr(pari_real_add,'native'):
        pari(1)
        assert ctypes.sizeof(ctypes.c_ulong)==8
        lib=ctypes.CDLL(None);pointer=ctypes.POINTER(ctypes.c_ulong)
        rr=lib.addrr_sign;rr.argtypes=[pointer,ctypes.c_long,pointer,ctypes.c_long];rr.restype=pointer
        ir=lib.addir_sign;ir.argtypes=rr.argtypes;ir.restype=pointer
        sr=lib.addsr;sr.argtypes=[ctypes.c_long,pointer];sr.restype=pointer
        get=lib.get_avma;get.argtypes=[];get.restype=ctypes.c_ulong
        restore=lib.set_avma;restore.argtypes=[ctypes.c_ulong];restore.restype=None
        pari_real_add.native=(rr,ir,sr,get,restore)
    rr,ir,sr,get,restore=pari_real_add.native
    mask64=(1<<64)-1
    def mantissa(p,mode):
        low=1<<(p-1);mask=low-1
        return low+[0,mask,1,mask-1,mask//3,mask//3*2,
                    mask-((1<<min(64,p-1))-1),mask//7,low//2-1][mode]
    def real(p,e,mode,sign):
        words=p//64+2 if sign else 2
        # Keep one physically allocated zero guard after the logical GEN.
        a=(ctypes.c_ulong*(words+1))();a[0]=(2<<57)|words
        a[1]=((sign&3)<<62)+(1<<61)+e;m=mantissa(p,mode)
        for i in range(words-2):a[i+2]=(m>>(p-64*(i+1)))&mask64
        return a
    y=real(py,ey,my,sy);saved=get()
    try:
        if op<2:
            x=real(px,ex,mx,sx);z=rr(x,sx,y,sy if op==0 else -sy)
        else:
            m=mantissa(px,mx);shift=ex-px+1;n=sx*(m<<shift if shift>=0 else m>>-shift)
            if op<4:
                mag=abs(n);words=(mag.bit_length()+63)//64+2
                x=(ctypes.c_ulong*words)();x[0]=(1<<57)|words
                sign=1 if n>0 else -1 if n<0 else 0;x[1]=((sign&3)<<62)|words
                for i in range(words-2):x[i+2]=(mag>>(i*64))&mask64
                z=ir(x,sign,y,sy if op==2 else -sy)
            else:z=sr(n if op==4 else -n,y)
        words=z[0]&((1<<56)-1);sign=z[1]>>62;sign=-1 if sign==3 else sign
        exponent=(z[1]&((1<<62)-1))-(1<<61);mag=0
        for i in range(2,words):mag=(mag<<64)+z[i]
        return json.dumps([sign,str(exponent),str(mag),(words-2)*64],separators=(',',':'))
    finally:restore(saved)
FUNCTIONS['pari_real_add']=pari_real_add

def pari_real_trans(op,p,n,shift,k):
    """Native PARI kernels; tune the installed word-precision ABI to bundled limits."""
    import ctypes,hashlib
    op,p,n,shift,k=map(int,[op,p,n,shift,k])
    pari(1);lib=ctypes.CDLL(None);P=ctypes.POINTER(ctypes.c_ulong)
    exp_limit=ctypes.c_long.in_dll(lib,'EXPNEWTON_LIMIT')
    log_limit=ctypes.c_long.in_dll(lib,'LOGAGM_LIMIT')
    inv_limit=ctypes.c_long.in_dll(lib,'INVNEWTON_LIMIT')
    limits=(exp_limit.value,log_limit.value,inv_limit.value)
    # Installed PARI uses word lengths including two header words; bundled PARI
    # uses bits. Match kernel/gmp/tune.h's 4224-bit and 384-bit thresholds.
    exp_limit.value=4224//64+2;log_limit.value=384//64+2;inv_limit.value=4800//64+2
    get=lib.get_avma;get.argtypes=[];get.restype=ctypes.c_ulong
    restore=lib.set_avma;restore.argtypes=[ctypes.c_ulong];restore.restype=None
    def native(name,args,restype=P):
        f=getattr(lib,name);f.argtypes=args;f.restype=restype;return f
    def intgen(value):
        mag=abs(value);words=(mag.bit_length()+63)//64+2
        a=(ctypes.c_ulong*words)();a[0]=(1<<57)|words
        sign=1 if value>0 else -1 if value<0 else 0;a[1]=((sign&3)<<62)|words
        for i in range(words-2):a[i+2]=(mag>>(64*i))&((1<<64)-1)
        return a
    def raw_real(x):
        prec=int(x.bitprecision()) if x else 0;exponent=int(x.exponent())
        words=prec//64+2;a=(ctypes.c_ulong*(words+1))();a[0]=(2<<57)|words
        sign=int(x.sign());a[1]=((sign&3)<<62)+(1<<61)+exponent
        mag=int(abs((x<<int(prec-1-exponent)).sage())) if sign else 0
        for i in range(words-2):a[i+2]=(mag>>(prec-64*(i+1)))&((1<<64)-1)
        return a
    def raw_frame(z):
        words=z[0]&((1<<56)-1);sign=z[1]>>62;sign=-1 if sign==3 else sign
        exponent=(z[1]&((1<<62)-1))-(1<<61);mag=0
        if not sign:return [0,str(exponent),'0',(words-2)*64]
        for i in range(2,words):mag=(mag<<64)+z[i]
        return [sign,str(exponent),str(mag),(words-2)*64]
    def gen_frame(z):
        if not z:return [0,str(z.exponent()),'0',0]
        prec=int(z.bitprecision());exponent=int(z.exponent())
        return [int(z.sign()),str(exponent),str(int(abs((z<<int(prec-1-exponent)).sage()))),prec]
    def integer_frame(z):
        words=z[1]&((1<<56)-1);sign=z[1]>>62;mag=0
        for i in range(words-1,1,-1):mag=(mag<<64)+z[i]
        return -mag if sign==3 else mag
    def hash_integer(value):
        mag=abs(value);bits=mag.bit_length()
        return [str(-1 if value<0 else 1 if value>0 else 0),str(bits),hashlib.sha256(mag.to_bytes((bits+7)//8,'big')).hexdigest()]
    try:
        if op in [3,7,8]:
            if op==3:z=pari.factorial(n,precision=p)
            else:
                from sage.arith.misc import factorial as original_factorial
                z=original_factorial(ZZ(n) if op==8 else n,'pari')
            result=gen_frame(z)
        elif op==10:
            z=pari.bernfrac(n);result=[str(z.numerator()),str(z.denominator())]
        elif op==12:
            result=str(native('quadratic_prec_mask',[ctypes.c_long],ctypes.c_ulong)(n))
        elif op in [9,11,13]:
            saved=get()
            try:
                if op==9:z=native('mpfact',[ctypes.c_long])(n)
                elif op==13:z=native('mulu_interval_step',[ctypes.c_ulong]*3)(n,shift,k)
                else:
                    numerator=intgen(n<<shift if shift>=0 else n>>-shift);denominator=intgen(k)
                    z=native('rdivii',[P,P,ctypes.c_long])(numerator,denominator,p//64+2)
                result=raw_frame(z) if op==11 else hash_integer(integer_frame(z))
            finally:restore(saved)
        elif op==0:result=gen_frame(pari.Pi(precision=p))
        else:
            unit=pari.factorial(0,precision=p)
            x=((pari(n)*unit) if n else unit-unit)<<int(shift-(not n))
            if op==1:result=gen_frame(x.exp())
            elif op==4:result=gen_frame(abs(x).log())
            elif op==5:result=gen_frame(x.agm(1))
            else:
                raw=raw_real(x);saved=get()
                try:
                    if op==2:z=native('exp1r_abs',[P])(raw)
                    elif op==5:z=native('agm1r_abs',[P])(raw)
                    elif op==6:z=native('powru',[P,ctypes.c_ulong])(raw,k)
                    elif op==14:z=native('invr',[P])(raw)
                    else:raise ValueError('invalid real kernel operation')
                    result=raw_frame(z)
                finally:restore(saved)
        return json.dumps(result,separators=(',',':'))
    finally:exp_limit.value,log_limit.value,inv_limit.value=limits
FUNCTIONS['pari_real_trans']=pari_real_trans

def pari_product_trace(op,size,seed):
    """Execute original product/powering schedules with instrumented callbacks."""
    import ctypes,hashlib
    op,size,seed=map(int,[op,size,seed]);pari(1);lib=ctypes.CDLL(None)
    P=ctypes.POINTER(ctypes.c_ulong);owners=[];trace=[];mask=(1<<64)-1
    def integer(value):
        words=(value.bit_length()+63)//64+2;a=(ctypes.c_ulong*words)()
        a[0]=(1<<57)|words;a[1]=((1 if value else 0)<<62)|words
        for i in range(words-2):a[2+i]=(value>>(64*i))&mask
        owners.append(a);return a
    def value(a):
        words=a[1]&((1<<56)-1);out=0
        for i in range(words-1,1,-1):out=(out<<64)+a[i]
        return out
    def product(a,b):
        trace.append(f'M:{a},{b}')
        return a*b%((1<<61)-1) if op==2 else ((a*65599)^b)&mask
    MUL=ctypes.CFUNCTYPE(ctypes.c_void_p,ctypes.c_void_p,P,P)
    SQR=ctypes.CFUNCTYPE(ctypes.c_void_p,ctypes.c_void_p,P)
    @MUL
    def mul(ctx,a,b):return ctypes.addressof(integer(product(value(a),value(b))))
    @SQR
    def sqr(ctx,a):
        x=value(a);trace.append(f'S:{x}')
        return ctypes.addressof(integer(x*x%((1<<61)-1)))
    get=lib.get_avma;get.argtypes=[];get.restype=ctypes.c_ulong
    restore=lib.set_avma;restore.argtypes=[ctypes.c_ulong];restore.restype=None
    saved=get()
    try:
        if op==2:
            base=integer(seed);f=lib.gen_powu_i;f.argtypes=[P,ctypes.c_ulong,ctypes.c_void_p,SQR,MUL];f.restype=P
            z=f(base,size,None,sqr,mul)
        else:
            vec=(ctypes.c_ulong*(size+1))();vec[0]=(17<<57)|(size+1);state=seed
            for i in range(size):
                state=(1664525*state+1013904223)&((1<<32)-1)
                vec[i+1]=ctypes.addressof(integer(state))
            if op==0:
                f=lib.ZV_prod;f.argtypes=[P];f.restype=P;z=f(vec)
            else:
                f=lib.gen_product;f.argtypes=[P,ctypes.c_void_p,MUL];f.restype=P;z=f(vec,None,mul)
        out=value(z);bits=out.bit_length()
        return json.dumps([str(bits),hashlib.sha256(out.to_bytes((bits+7)//8,'big')).hexdigest(),hashlib.sha256('|'.join(trace).encode()).hexdigest()],separators=(',',':'))
    finally:restore(saved)
FUNCTIONS['pari_product_trace']=pari_product_trace

def pari_real_div(op,px,py,ex,ey,mx,my,sx,sy):
    """Native C quotients with deterministic operand padding and bundled tuning."""
    import ctypes
    op,px,py,ex,ey,mx,my,sx,sy=map(int,[op,px,py,ex,ey,mx,my,sx,sy])
    pari(1);lib=ctypes.CDLL(None);P=ctypes.POINTER(ctypes.c_ulong)
    div_limit=ctypes.c_long.in_dll(lib,'DIVRR_GMP_LIMIT')
    inv_limit=ctypes.c_long.in_dll(lib,'INVNEWTON_LIMIT')
    limits=(div_limit.value,inv_limit.value)
    div_limit.value=256//64+2;inv_limit.value=4800//64+2
    get=lib.get_avma;get.argtypes=[];get.restype=ctypes.c_ulong
    restore=lib.set_avma;restore.argtypes=[ctypes.c_ulong];restore.restype=None
    mask=(1<<64)-1
    def mantissa(p,mode):
        if mode<0:return -mode
        low=1<<(p-1);tail=low-1
        return low+[0,tail,1,tail-1,tail//3,tail//3*2,tail-((1<<min(64,p-1))-1),tail//7][mode]
    def integer_value(p,e,mode,sign):
        m=mantissa(p,mode);shift=e-p+1
        return sign*(m<<shift if shift>=0 else m>>-shift)
    def integer(p,e,mode,sign):
        n=integer_value(p,e,mode,sign);mag=abs(n);words=(mag.bit_length()+63)//64+2
        a=(ctypes.c_ulong*words)();a[0]=(1<<57)|words;a[1]=((sign&3)<<62)|words if n else words
        for i in range(words-2):a[i+2]=(mag>>(64*i))&mask
        return a
    def real(p,e,mode,sign):
        words=p//64+2;a=(ctypes.c_ulong*(words+1))();a[0]=(2<<57)|words
        a[1]=((sign&3)<<62)+(1<<61)+e;m=mantissa(p,mode)
        for i in range(words-2):a[i+2]=(m>>(p-64*(i+1)))&mask
        return a
    saved=get()
    try:
        x=(integer if op==1 else real)(px,ex,mx,sx)
        y=integer_value(py,ey,my,sy) if op==3 else (integer if op==2 else real)(py,ey,my,sy)
        f=getattr(lib,['divrr','divir','divri','divru'][op]);f.argtypes=[P,ctypes.c_ulong if op==3 else P];f.restype=P
        z=f(x,y);words=z[0]&((1<<56)-1);sign=z[1]>>62;sign=-1 if sign==3 else sign
        exponent=(z[1]&((1<<62)-1))-(1<<61);m=0
        if not sign:return json.dumps([0,str(exponent),'0',0],separators=(',',':'))
        for i in range(2,words):m=(m<<64)+z[i]
        return json.dumps([sign,str(exponent),str(m),(words-2)*64],separators=(',',':'))
    finally:
        restore(saved);div_limit.value,inv_limit.value=limits
FUNCTIONS['pari_real_div']=pari_real_div


_pari_zero_native=None
def pari_real_zero(op,p,q,e,n):
    """Preserve zero allocation; compare inverse origin under documented error mapping."""
    global _pari_zero_native
    if _pari_zero_native is None:
        import os,sys
        from pathlib import Path
        from sage.misc.cython import cython_import
        pari(1);sys.stdout.flush();saved=os.dup(1)
        try:
            os.dup2(2,1)
            _pari_zero_native=cython_import(str(Path(__file__).parent.parent/'pari_real_zero_native.pyx'))
        finally:sys.stdout.flush();os.dup2(saved,1);os.close(saved)
    try:result=_pari_zero_native.native_zero(*map(int,[op,p,q,e,n]))
    except Exception as error:
        if type(error).__name__!='PariError' or not str(error).startswith('impossible inverse in '):raise
        # The dependency already documents separate error subclasses and no GEN rendering.
        result=['inverse',str(error).split(':',1)[0]]
    return json.dumps(result,separators=(',',':'))
FUNCTIONS['pari_real_zero']=pari_real_zero


def pari_real_mul(op,px,py,ex,ey,mx,my,sx,sy):
    """Native C products with deterministic operand padding and bundled tuning."""
    import ctypes
    op,px,py,ex,ey,mx,my,sx,sy=map(int,[op,px,py,ex,ey,mx,my,sx,sy])
    pari(1);lib=ctypes.CDLL(None);P=ctypes.POINTER(ctypes.c_ulong)
    mul_limit=ctypes.c_long.in_dll(lib,'MULRR_MULII_LIMIT')
    sqr_limit=ctypes.c_long.in_dll(lib,'SQRR_SQRI_LIMIT')
    limits=(mul_limit.value,sqr_limit.value)
    mul_limit.value=3520//64+2;sqr_limit.value=768//64+2
    get=lib.get_avma;get.argtypes=[];get.restype=ctypes.c_ulong
    restore=lib.set_avma;restore.argtypes=[ctypes.c_ulong];restore.restype=None
    mask=(1<<64)-1
    def mantissa(p,mode):
        if mode<0:return -mode
        low=1<<(p-1);tail=low-1
        return low+[0,tail,1,tail-1,tail//3,tail//3*2,tail-((1<<min(64,p-1))-1),tail//7][mode]
    def integer_value(p,e,mode,sign):
        m=mantissa(p,mode);shift=e-p+1
        return sign*(m<<shift if shift>=0 else m>>-shift)
    def integer(p,e,mode,sign):
        n=integer_value(p,e,mode,sign);mag=abs(n);words=(mag.bit_length()+63)//64+2
        a=(ctypes.c_ulong*words)();a[0]=(1<<57)|words;a[1]=((sign&3)<<62)|words if n else words
        for i in range(words-2):a[i+2]=(mag>>(64*i))&mask
        return a
    def real(p,e,mode,sign):
        words=p//64+2;a=(ctypes.c_ulong*(words+1))();a[0]=(2<<57)|words
        a[1]=((sign&3)<<62)+(1<<61)+e;m=mantissa(p,mode)
        for i in range(words-2):a[i+2]=(m>>(p-64*(i+1)))&mask
        return a
    saved=get()
    try:
        x=integer_value(px,ex,mx,sx) if op==5 else (integer if op==1 else real)(px,ex,mx,sx)
        y=integer_value(py,ey,my,sy) if op==6 else (integer if op==2 else real)(py,ey,my,sy)
        f=getattr(lib,['mulrr','mulir','mulri','sqrr','mulrr','mulsr','mulrs'][op]);f.restype=P
        f.argtypes=[P] if op==3 else [ctypes.c_long if op==5 else P,ctypes.c_long if op==6 else P]
        z=f(x) if op==3 else f(x,x) if op==4 else f(x,y)
        words=z[0]&((1<<56)-1);sign=z[1]>>62;sign=-1 if sign==3 else sign
        exponent=(z[1]&((1<<62)-1))-(1<<61);m=0
        if not sign:return json.dumps([0,str(exponent),'0',(words-2)*64],separators=(',',':'))
        for i in range(2,words):m=(m<<64)+z[i]
        return json.dumps([sign,str(exponent),str(m),(words-2)*64],separators=(',',':'))
    finally:
        restore(saved);mul_limit.value,sqr_limit.value=limits
FUNCTIONS['pari_real_mul']=pari_real_mul


def pari_real_sqrt(op,p,e,m,sign):
    """Native GMP-backed real roots, complex signs, exact integer root/remainder."""
    import ctypes,hashlib
    op,p,e,m,sign=map(int,[op,p,e,m,sign])
    pari(1);lib=ctypes.CDLL(None);P=ctypes.POINTER(ctypes.c_ulong)
    get=lib.get_avma;get.argtypes=[];get.restype=ctypes.c_ulong
    restore=lib.set_avma;restore.argtypes=[ctypes.c_ulong];restore.restype=None
    def frame(z):
        typ=z[0]>>57;words=z[0]&((1<<56)-1)
        if typ==6:return ['complex',frame(ctypes.cast(z[1],P)),frame(ctypes.cast(z[2],P))]
        if typ==1:
            words=z[1]&((1<<56)-1)
            n=0
            for i in range(words-1,1,-1):n=(n<<64)+z[i]
            return ['integer',str(-n if z[1]>>62==3 else n)]
        s=z[1]>>62;s=-1 if s==3 else s;ex=(z[1]&((1<<62)-1))-(1<<61);mantissa=0
        if s:
            for i in range(2,words):mantissa=(mantissa<<64)+z[i]
        return ['real',s,str(ex),str(mantissa),(words-2)*64]
    def integer_root(n,root_only=False):
        mag=abs(n);words=(mag.bit_length()+63)//64+2;x=(ctypes.c_ulong*words)()
        x[0]=(1<<57)|words;x[1]=((3 if n<0 else 1 if n else 0)<<62)|words
        for i in range(words-2):x[i+2]=(mag>>(64*i))&((1<<64)-1)
        f=lib.sqrti if root_only else lib.sqrtremi
        f.argtypes=[P] if root_only else [P,ctypes.POINTER(P)];f.restype=P;remainder=P();saved=get()
        try:
            root=f(x) if root_only else f(x,ctypes.byref(remainder))
            return [frame(root)[1]] if root_only else [frame(root)[1],frame(remainder)[1]]
        finally:restore(saved)
    if op in [4,5]:return json.dumps(['integer-rootonly']+integer_root((m<<e)+sign,True),separators=(',',':'))
    if op==2:return json.dumps(['integer-root']+integer_root((m<<e)+sign),separators=(',',':'))
    if op==3:
        digest=hashlib.sha256()
        for i in range(m):digest.update((','.join(integer_root(e+i*sign))+';').encode())
        return digest.hexdigest()
    words=p//64+2;x=(ctypes.c_ulong*(words+1))();x[0]=(2<<57)|words
    x[1]=((sign&3)<<62)+(1<<61)+e
    for i in range(words-2):x[i+2]=(m>>(p-64*(i+1)))&((1<<64)-1)
    f=getattr(lib,'sqrtr' if op else 'sqrtr_abs');f.argtypes=[P];f.restype=P;saved=get()
    try:return json.dumps(frame(f(x)),separators=(',',':'))
    finally:restore(saved)
FUNCTIONS['pari_real_sqrt']=pari_real_sqrt


def pari_real_integer_error(op,p,e,mode,s):
    """gcvtoi's error exponent is the native subtraction exponent, including zero."""
    import ctypes
    op,p,e,mode,s=map(int,[op,p,e,mode,s])
    pari(1);lib=ctypes.CDLL(None);P=ctypes.POINTER(ctypes.c_ulong)
    get=lib.get_avma;get.argtypes=[];get.restype=ctypes.c_ulong
    restore=lib.set_avma;restore.argtypes=[ctypes.c_ulong];restore.restype=None
    words=p//64+2;x=(ctypes.c_ulong*(words+1))();x[0]=(2<<57)|words
    x[1]=((s&3)<<62)+(1<<61)+e;low=1<<max(0,p-1);tail=low-1
    m=low+[0,tail,1,max(0,tail-1),tail//3,tail//3*2,tail//7,tail//5][mode] if s else 0
    for i in range(words-2):x[i+2]=(m>>(p-64*(i+1)))&((1<<64)-1)
    f=lib.gcvtoi;f.argtypes=[P,ctypes.POINTER(ctypes.c_long)];f.restype=P
    error=ctypes.c_long();saved=get()
    try:
        z=f(x,ctypes.byref(error));length=z[1]&((1<<56)-1);n=0
        for i in range(length-1,1,-1):n=(n<<64)+z[i]
        if z[1]>>62==3:n=-n
        return json.dumps([str(n),str(error.value)],separators=(',',':'))
    finally:restore(saved)
FUNCTIONS['pari_real_integer_error']=pari_real_integer_error


_pari_double_native=None
def pari_real_double(op,p,e,m,sign):
    global _pari_double_native
    if _pari_double_native is None:
        import os,sys
        from pathlib import Path
        from sage.misc.cython import cython_import
        pari(1);sys.stdout.flush();saved=os.dup(1)
        try:
            os.dup2(2,1)
            _pari_double_native=cython_import(str(Path(__file__).parent.parent/'pari_real_double_native.pyx'))
        finally:sys.stdout.flush();os.dup2(saved,1);os.close(saved)
    return json.dumps(_pari_double_native.native_double(*map(int,[op,p,e,m,sign])),separators=(',',':'))
FUNCTIONS['pari_real_double']=pari_real_double


_pari_compare_native=None
def pari_real_compare(p,e,m,s,q,f,n,t):
    global _pari_compare_native
    if _pari_compare_native is None:
        import os,sys
        from pathlib import Path
        from sage.misc.cython import cython_import
        pari(1);sys.stdout.flush();saved=os.dup(1)
        try:
            os.dup2(2,1)
            _pari_compare_native=cython_import(str(Path(__file__).parent.parent/'pari_real_compare_native.pyx'))
        finally:sys.stdout.flush();os.dup2(saved,1);os.close(saved)
    return _pari_compare_native.native_compare(*map(int,[p,e,m,s,q,f,n,t]))
FUNCTIONS['pari_real_compare']=pari_real_compare

FUNCTIONS['pari_buch_real_add']=pari_real_add
FUNCTIONS['pari_buch_real_div']=pari_real_div
FUNCTIONS['pari_buch_real_mul']=pari_real_mul

FUNCTIONS['pari_buch_real_zero']=pari_real_zero

_pari_buch_convert_native=None
def pari_buch_convert(op,p,e,m,s,q,n):
    global _pari_buch_convert_native
    if _pari_buch_convert_native is None:
        import os,sys
        from pathlib import Path
        from sage.misc.cython import cython_import
        pari(1);sys.stdout.flush();saved=os.dup(1)
        try:
            os.dup2(2,1)
            _pari_buch_convert_native=cython_import(str(Path(__file__).parent.parent/'pari_buch_convert_native.pyx'))
        finally:sys.stdout.flush();os.dup2(saved,1);os.close(saved)
    return json.dumps(_pari_buch_convert_native.native_convert(*map(int,[op,p,e,m,s,q,n])),separators=(',',':'))
FUNCTIONS['pari_buch_convert']=pari_buch_convert

FUNCTIONS['pari_real_truncate']=pari_buch_convert

FUNCTIONS['pari_buch_real_sqrt']=pari_real_sqrt

_pari_logexp_native=None
def pari_real_logexp(op,p,e,m,s,q):
    global _pari_logexp_native
    if _pari_logexp_native is None:
        import os,sys
        from pathlib import Path
        from sage.misc.cython import cython_import
        pari(1);sys.stdout.flush();saved=os.dup(1)
        try:
            os.dup2(2,1)
            _pari_logexp_native=cython_import(str(Path(__file__).parent.parent/'pari_real_logexp_native.pyx'))
        finally:sys.stdout.flush();os.dup2(saved,1);os.close(saved)
    return json.dumps(_pari_logexp_native.native_logexp(*map(int,[op,p,e,m,s,q])),separators=(',',':'))
FUNCTIONS['pari_real_logexp']=pari_real_logexp

_pari_split_native=None
def pari_binary_split(op,n,a,b):
    if op:
        import ctypes
        pari(1);lib=ctypes.CDLL(None);P=ctypes.POINTER(ctypes.c_long)
        class AStruct(ctypes.Structure):_fields_=[(name,ctypes.POINTER(P)) for name in ['a','b','p','q']]
        class RStruct(ctypes.Structure):_fields_=[(name,P) for name in ['P','Q','B','T']]
        get=lib.get_avma;get.argtypes=[];get.restype=ctypes.c_ulong
        restore=lib.set_avma;restore.argtypes=[ctypes.c_ulong];restore.restype=None
        init=lib.abpq_init;init.argtypes=[ctypes.POINTER(AStruct),ctypes.c_long];init.restype=None
        run=lib.abpq_sum;run.argtypes=[ctypes.POINTER(RStruct),ctypes.c_long,ctypes.c_long,ctypes.POINTER(AStruct)];run.restype=None
        make=lib.stoi;make.argtypes=[ctypes.c_long];make.restype=P
        saved=get()
        try:
            A=AStruct();R=RStruct();size=int(n);offset=int(b);state=int(a)
            init(ctypes.byref(A),size+offset)
            for i in range(size+offset+1):
                for field,mod,shift in [('a',17,8),('b',19,9),('p',23,11),('q',29,14)]:
                    state=(1664525*state+1013904223)&0xffffffff
                    getattr(A,field)[i]=make(state%mod-shift)
            run(ctypes.byref(R),offset,offset+size,ctypes.byref(A))
            def integer(x):
                words=x[1]&((1<<56)-1);value=0
                for i in range(words-1,1,-1):value=(value<<64)+(x[i]&((1<<64)-1))
                return str(-value if x[1]<0 else value)
            return json.dumps([integer(getattr(R,k)) for k in ['P','Q','B','T']],separators=(',',':'))
        finally:restore(saved)
    global _pari_split_native
    if _pari_split_native is None:
        import os,sys
        from pathlib import Path
        from sage.misc.cython import cython_import
        pari(1);sys.stdout.flush();saved=os.dup(1)
        try:
            os.dup2(2,1)
            _pari_split_native=cython_import(str(Path(__file__).parent.parent/'pari_split_native.pyx'))
        finally:sys.stdout.flush();os.dup2(saved,1);os.close(saved)
    return json.dumps(_pari_split_native.native_split(*map(int,[op,n,a,b])),separators=(',',':'))
FUNCTIONS['pari_binary_split']=pari_binary_split


_pari_exp1_native=None
def pari_real_exp1(p,e,m,s):
    global _pari_exp1_native
    if _pari_exp1_native is None:
        import os,sys
        from pathlib import Path
        from sage.misc.cython import cython_import
        pari(1);sys.stdout.flush();saved=os.dup(1)
        try:
            os.dup2(2,1)
            _pari_exp1_native=cython_import(str(Path(__file__).parent.parent/'pari_exp1_native.pyx'))
        finally:sys.stdout.flush();os.dup2(saved,1);os.close(saved)
    return json.dumps(_pari_exp1_native.native_exp1(*map(int,[p,e,m,s])),separators=(',',':'))
FUNCTIONS['pari_real_exp1']=pari_real_exp1
