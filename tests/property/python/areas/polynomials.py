"""SageMath side of the ``polynomials`` property-test area (GF(p)[x]).

Cases: tests/property/cases/polynomials.cases.json
"""

from sage.all import *

from ._helpers import (
    poly_add,
    poly_derivative,
    poly_eval,
    poly_factor,
    poly_gcd,
    poly_is_irreducible,
    poly_mod,
    poly_mul,
    poly_pow,
    poly_quo_rem,
    poly_roots,
)

FUNCTIONS = {
    'poly_add': poly_add,
    'poly_mul': poly_mul,
    'poly_quo_rem': poly_quo_rem,
    'poly_mod': poly_mod,
    'poly_gcd': poly_gcd,
    'poly_eval': poly_eval,
    'poly_factor': poly_factor,
    'poly_derivative': poly_derivative,
    'poly_is_irreducible': poly_is_irreducible,
    'poly_roots': poly_roots,
    'poly_pow': poly_pow,
}


import json

def comparison(run):
    def normalized(x):
        if isinstance(x, (list, tuple)): return [normalized(a) for a in x]
        if isinstance(x, bool) or x is None: return x
        return str(x)
    try: return json.dumps({'value': normalized(run())}, separators=(',', ':'))
    except Exception as e: return json.dumps({'error': type(e).__name__, 'message': str(e)}, separators=(',', ':'))

def polynomial_base(kind):
    return QQ if kind == 0 else GF(7) if kind in (1, 2) else Zmod(14) if kind == 3 else GF(49, 'a', impl='pari_ffelt') if kind == 4 else GF(2)

def polynomial_input(kind, n):
    if kind == 26: return Zmod(1)(n)
    if kind == 24: return float(1)/3
    if kind == 25: return float.fromhex("0x0.0000000000001p-1022")
    if kind == 21: return GF(49, 'a', impl='pari_ffelt')(n)
    if kind == 22: return GF(2)(n)
    if kind == 23: return GF(7)(n)
    return Integer(n) if kind in (0, 1) else QQ(n)/2 if kind == 2 else bool(n) if kind == 3 else str(n) if kind == 4 else None if kind == 5 else float(n)+0.5 if kind == 6 else GF(7)(n) if kind == 7 else Zmod(14)(n) if kind == 8 else GF(5)(n) if kind == 9 else [] if kind == 10 else [n] if kind == 11 else GF(49, 'a', impl='pari_ffelt').gen()+n if kind == 12 else [n, 1] if kind == 13 else [[n, 1]] if kind == 14 else [[n]] if kind == 15 else [None] if kind == 16 else [[]] if kind == 17 else float('nan')

def polynomial_constructor_input(base, kind, n, shape):
    R = PolynomialRing(polynomial_base(base), 'x')
    c = polynomial_input(kind, n)
    x = [c] if shape == 0 else [c, 0, 1, 0] if shape == 1 else [] if shape == 2 else None
    f = R() if shape == 4 else R(x)
    return [list(map(str, f.list())), f.degree(), str(f.parent()), all(a == R.base_ring()(a) for a in f.list())]
FUNCTIONS['polynomial_constructor_input'] = lambda *args: comparison(lambda: polynomial_constructor_input(*args))
FUNCTIONS['polynomial_gen_input'] = lambda base, kind, n: comparison(lambda: str(PolynomialRing(polynomial_base(base), 'x').gen(polynomial_input(kind,n))))

FUNCTIONS['ff_prime__rational_'] = lambda p, n, legacy: comparison(lambda: GF(p)(n)._rational_())
FUNCTIONS['modular__rational_'] = lambda p, n: comparison(lambda: Zmod(p)(n)._rational_())
def polynomial_coefficient_conversion(base, source, n):
    R = polynomial_base(base)
    F = GF(49, 'a', impl='pari_ffelt')
    x = Zmod(14)(n) if source == 0 else GF(7)(n) if source == 1 else F(list(Integer(n).digits(7))) if source == 2 else Zmod(5)(n)
    return [str(R(x)), str(R)]
FUNCTIONS['polynomial_coefficient_conversion'] = lambda *args: comparison(lambda: polynomial_coefficient_conversion(*args))

def polynomial_parent_constructor(base, nested):
    R = PolynomialRing(polynomial_base(base), 'y')
    f = R([2, 1])
    if not nested: return R(f) is f
    S = PolynomialRing(R, 'x')
    g = S(f)
    return [g.degree(), list(map(str, g.list())), str(g.parent())]
FUNCTIONS['polynomial_parent_constructor'] = lambda *args: comparison(lambda: polynomial_parent_constructor(*args))

def unary_polynomial_base(kind):
    if kind in (13,14):return Zmod(2**127-(2 if kind==13 else 1))
    return GF(2) if kind == 10 else ZZ if kind == 6 else Zmod(1) if kind == 7 else GF(2**127-1) if kind == 8 else GF(8, 'a', impl='pari_ffelt') if kind == 9 else polynomial_base(kind)
def polynomial_monic(base, coeffs, denominator):
    R = PolynomialRing(unary_polynomial_base(base), 'x')
    f = R([c if denominator == 1 else QQ(c)/denominator for c in coeffs])
    g = f.monic()
    return [list(map(str,g.list())), str(g.parent()), g is f]
FUNCTIONS['polynomial_monic'] = lambda *args: comparison(lambda: polynomial_monic(*args))
FUNCTIONS['polynomial_is_monic'] = lambda base,coeffs: comparison(lambda: PolynomialRing(unary_polynomial_base(base), 'x')(coeffs).is_monic())
def polynomial_monic_extension(p, degree, coeffs):
    F = GF(p**degree, 'a', impl='pari_ffelt')
    R = PolynomialRing(F, 'x')
    f = R([F(list(Integer(c).digits(p))) for c in coeffs])
    g = f.monic()
    return [list(map(str,g.list())), str(g.parent()), g is f]
FUNCTIONS['polynomial_monic_extension'] = lambda *args: comparison(lambda: polynomial_monic_extension(*args))

FUNCTIONS['gf2_zero_protocol'] = lambda operation,a,b: comparison(lambda: str(GF(2)(a)/GF(2)(b) if operation == 0 else ~GF(2)(a) if operation == 1 else GF(2)(a)**(float(b)+0.5 if operation == 3 else Integer(b))))

def polynomial_index_operation(base, operation, coeffs, kind, n):
    B = PolynomialRing(QQ, 'y') if base == 11 else unary_polynomial_base(base)
    R = PolynomialRing(B, 'x')
    f = R(coeffs)
    index = float('inf') if kind == 19 else float('-inf') if kind == 20 else polynomial_input(kind,n)
    if operation == 3: return str(f[index])
    if operation == 0 and base in (0,6) and f and kind in (2,6,7,8,9,12,19,21) and index > 0:
        # Sage 10.3's failed [QQ.zero()] * non-index scalar decrefs the cached
        # zero and eventually segfaults. Empty-list multiplication executes the
        # same rejected sequence dispatch without corrupting a coefficient.
        # Bundled Polynomial.shift (polynomial_element.pyx:10244) uses this
        # sequence-repetition protocol, which rejects even integral Rationals.
        [] * index
    if operation == 2 and base in (5,7,8,10,11) and kind == 0 and 2**63 <= index < 2**64:
        [] * (index + 1 - len(f.list()))
    if operation == 2 and base in (5,7,8,10,11) and kind == 21 and len(f.list()) < index + 1:
        [] * (index + 1 - len(f.list()))
    if operation == 2 and base in (5,7,8,10,11) and kind == 2 and index.denominator() == 1 and 0 <= index < 2**64 and index + 1 > len(f.list()):
        # Generic reverse has the same compiled error-cleanup bug when padding
        # with a Rational count; the unsigned-long check has already succeeded.
        return [] * (index + 1 - len(f.list()))
    g = f.shift(index) if operation == 0 else f.truncate(index) if operation == 1 else f.reverse(index)
    return None if g is None else [list(map(str,g.list())), str(g.parent()), g is f]
def polynomial_index_comparison(*args):
    result = comparison(lambda: polynomial_index_operation(*args))
    # Bundled polynomial_element.pyx:8143 and specialized reverse methods use
    # "nonnegative"; installed Sage 10.3 uses the older "non-negative" spelling.
    if args[1] == 2:
        payload = json.loads(result)
        if payload.get('error') == 'ValueError' and payload.get('message','').startswith('degree argument must be a non-negative integer, got '):
            payload['message'] = payload['message'].replace('non-negative', 'nonnegative', 1)
            result = json.dumps(payload, separators=(',', ':'))
    return result
FUNCTIONS['polynomial_index_operation'] = polynomial_index_comparison

def polynomial_mixed_arithmetic(left, right, operation, a, b, different_variable, shape=0, variable_kind=-1, variable_n=0):
    def base(k): return RDF if k == 17 else Zmod([6,10,15][k-12]) if 12 <= k <= 14 else GF(5) if k == 15 else GF(11) if k == 16 else PolynomialRing(QQ, 't') if k == 11 else unary_polynomial_base(k)
    R = PolynomialRing(base(left), 't' if different_variable == 3 else 'x')
    S = PolynomialRing(base(right), 't' if different_variable == 2 else 'y' if different_variable == 1 else 'x')
    def coeffs(B, k, cs): return [B.gen()+c if shape == 3 and k in (4,9) else B.gen()+c if shape == 1 and k == 11 else QQ(c)/2 if shape == 2 and k in (0,11) else c for c in cs]
    f = R(coeffs(R.base_ring(), left, a))
    h = f if different_variable == 4 else S(coeffs(S.base_ring(), right, b))
    if operation == 11:
        g=f.multiplication_trunc(h,variable_n)
        return [list(map(str,g.list())),str(g.parent()),g is f,g is h]
    if operation == 3: return f == h
    if operation == 10:
        result=f.pseudo_quo_rem(h)
        return [[str(g),str(g.parent()),g is f,g is h,[g is k for k in result]] for g in result]
    if operation == 9:
        g = f.resultant(h) if variable_kind == -1 else f.resultant(h,proof=polynomial_input(variable_kind,variable_n))
        return [str(g),str(g.parent())]
    if operation == 8:
        variable = R.gen() if variable_kind == 26 else S.gen() if variable_kind == 27 else PolynomialRing(QQ, 'z').gen() if variable_kind == 28 else polynomial_input(variable_kind,variable_n)
        matrix = f.sylvester_matrix(h) if variable_kind == -1 else f.sylvester_matrix(h,variable)
        return [list(map(str,row)) for row in matrix.rows()]
    if operation == 7:
        result = polynomial_reference_gcd(f,h,'xgcd')
        from sage.rings.polynomial.polynomial_element import Polynomial
        return [[list(map(str,g.list())) if isinstance(g,Polynomial) else str(g), str(g.parent()), g is f, g is h, [g is x for x in result]] for g in result]
    if operation == 4:
        q, r = f.quo_rem(h)
        return [[list(map(str, g.list())), str(g.parent()), g is f, g is h] for g in (q,r)]
    g = polynomial_reference_gcd(f, h) if operation == 6 else f % h if operation == 5 else f + h if operation == 0 else f - h if operation == 1 else f * h
    return [list(map(str, g.list())), str(g.parent()), g is f, g is h]
FUNCTIONS['polynomial_mixed_arithmetic'] = lambda *args: comparison(lambda: polynomial_mixed_arithmetic(*args))

def polynomial_repr(base, coeffs, shape):
    B = PolynomialRing(QQ, 't') if base == 11 else unary_polynomial_base(base)
    R = PolynomialRing(B, 'x')
    values = [QQ(c)/2 if shape == 1 else (B.gen()**(i+1)+1)*c if shape == 2 else c for i,c in enumerate(coeffs)]
    return str(R(values))
FUNCTIONS['polynomial_repr'] = lambda *args: comparison(lambda: polynomial_repr(*args))

def polynomial_scalar_equal(base, coeffs, kind, n, denominator=1):
    B = PolynomialRing(QQ, 't') if base == 11 else unary_polynomial_base(base)
    f = PolynomialRing(B, 'x')([c if denominator == 1 else QQ(c)/denominator for c in coeffs])
    scalar = float('inf') if kind == 19 else float('-inf') if kind == 20 else polynomial_input(kind,n)
    return f == scalar
FUNCTIONS['polynomial_scalar_equal'] = lambda *args: comparison(lambda: polynomial_scalar_equal(*args))

_polynomial_gcd_guard = None

def polynomial_reference_gcd(f, h, operation='gcd'):
    # Bundled polynomial_template.pxi now short-circuits equal nonzero inputs
    # through monic(). Installed Sage 10.3 predates that fix and its composite
    # native GCD can abort outside sig_on(). Execute the actual native operation
    # inside a signal guard, matching nmod_poly_linkage.pxi's current error guard.
    from sage.structure.element import canonical_coercion
    a, b = canonical_coercion(f, h)
    template = type(a).__name__ in ('Polynomial_zmod_flint', 'Polynomial_GF2X', 'Polynomial_ZZ_pEX')
    if operation == 'gcd' and template and a and b and a == b:
        return a.monic()
    if type(a).__name__ != 'Polynomial_zmod_flint' or a.base_ring().is_field():
        return getattr(a, operation)(b)
    global _polynomial_gcd_guard
    import os, sys, ctypes
    from pathlib import Path
    sys.stdout.flush()
    libc = ctypes.CDLL(None)
    libc.fflush(None)
    saved = os.dup(1)
    try:
        # Keep compiler and FLINT native diagnostics off the JSON transcript.
        os.dup2(2, 1)
        if _polynomial_gcd_guard is None:
            from sage.misc.cython import cython_import
            _polynomial_gcd_guard = cython_import(str(Path(__file__).parent.parent / 'polynomial_gcd_guard.pyx'))
        try:
            return getattr(_polynomial_gcd_guard, 'guarded_' + operation)(a, b)
        except RuntimeError:
            if operation == 'xgcd':
                raise ValueError('non-invertible elements encountered during XGCD') from None
            raise RuntimeError('FLINT gcd calculation failed') from None
    finally:
        sys.stdout.flush()
        libc.fflush(None)
        os.dup2(saved, 1)
        os.close(saved)

def polynomial_sylvester_operand(base, kind, n, coeffs):
    R = PolynomialRing(PolynomialRing(QQ, 't') if base == 11 else unary_polynomial_base(base), 'x')
    return [list(map(str,row)) for row in R(coeffs).sylvester_matrix(polynomial_input(kind,n)).rows()]
FUNCTIONS['polynomial_sylvester_operand'] = lambda *args: comparison(lambda: polynomial_sylvester_operand(*args))

def polynomial_constructor_scalar(base, kind, n):
    R = PolynomialRing(PolynomialRing(QQ, 't') if base == 11 else unary_polynomial_base(base), 'x')
    f = R(polynomial_input(kind,n))
    return [f.degree(),list(map(str,f.list())),str(f.parent())]
FUNCTIONS['polynomial_constructor_scalar'] = lambda *args: comparison(lambda: polynomial_constructor_scalar(*args))

def polynomial_sylvester_degree_one(p, a, b):
    R = PolynomialRing(GF(p), 'x')
    return [list(map(str,row)) for row in R(a).sylvester_matrix(R(b)).rows()]
FUNCTIONS['polynomial_sylvester_degree_one'] = lambda *args: comparison(lambda: polynomial_sylvester_degree_one(*args))

def polynomial_resultant_operand(base, kind, n, coeffs):
    R = PolynomialRing(PolynomialRing(QQ, 't') if base == 11 else unary_polynomial_base(base), 'x')
    g = R(coeffs).resultant(polynomial_input(kind,n))
    return [str(g),str(g.parent())]
FUNCTIONS['polynomial_resultant_operand'] = lambda *args: comparison(lambda: polynomial_resultant_operand(*args))

def rdf_coefficient(operation, a, b):
    import struct
    from sage.all import RDF
    x = RDF(struct.unpack('>d',int(a%(2**64)).to_bytes(8,'big'))[0])
    y = RDF(struct.unpack('>d',int(b%(2**64)).to_bytes(8,'big'))[0])
    if operation == 6: return x == y
    if operation == 9: return [str(RDF),RDF.characteristic(),RDF.is_field(),RDF.zero() is RDF.zero(),RDF.one() is RDF.one(),RDF(x) is x]
    r = x if operation == 0 else x+y if operation == 1 else x-y if operation == 2 else x*y if operation == 3 else x/y if operation == 4 else -x if operation == 5 else RDF.zero() if operation == 7 else RDF.one() if operation == 8 else RDF(QQ(a)/b) if operation == 10 else RDF(a) if operation == 11 else RDF(bool(a))
    return ['NaN' if r.is_NaN() else str(int.from_bytes(struct.pack('>d',float(r)),'big')),str(r),r.is_zero(),str(r.parent())]
FUNCTIONS['rdf_coefficient'] = lambda *args: comparison(lambda: rdf_coefficient(*args))

def rdf_polynomial_resultant(a,b):
    from sage.all import RDF
    import struct
    cv=lambda raw: RDF(struct.unpack('>d',int(raw%(2**64)).to_bytes(8,'big'))[0])
    R=PolynomialRing(RDF,'x')
    r=R([cv(x) for x in a]).resultant(R([cv(x) for x in b]))
    return ['NaN' if r.is_NaN() else str(int.from_bytes(struct.pack('>d',float(r)),'big')),str(r),str(r.parent())]
FUNCTIONS['rdf_polynomial_resultant'] = lambda *args: comparison(lambda: rdf_polynomial_resultant(*args))

def derivative_arguments(R, form, n):
    x=R.gen()
    B=R.base_ring()
    z=B.gen() if hasattr(B,'gen') else 0
    return [] if form==0 else [n] if form==1 else [Integer(n)] if form==2 else [QQ(n)/2] if form==3 else [bool(n)] if form==4 else [str(n)] if form==5 else [None] if form==6 else [x] if form==7 else [n*x] if form==8 else [PolynomialRing(B,'y').gen()] if form==9 else [PolynomialRing(QQ,'x').gen()] if form==10 else [z] if form==11 else [[]] if form==12 else [[x]] if form==13 else [[x,n]] if form==14 else [x,n] if form==15 else [x,n,n] if form==16 else [n,x,n] if form==17 else [n,None,n] if form==18 else [[None,x]] if form==19 else [[None,[]]] if form==20 else [float('nan')] if form==21 else [float(n)+0.5]

def polynomial_derivative_protocol(base, method, form, n, coeffs, shape):
    B=RDF if base==17 else PolynomialRing(ZZ,'t') if base==18 else PolynomialRing(RDF,'t') if base==19 else PolynomialRing(QQ,'t') if base==11 else unary_polynomial_base(base)
    R=PolynomialRing(B,'x')
    if shape==3:
        import struct
        values=[RDF(struct.unpack('>d',int(c%(2**64)).to_bytes(8,'big'))[0]) for c in coeffs]
    else: values=[B.gen()+c if shape==2 else QQ(c)/2 if shape==1 else c for c in coeffs]
    f=R(values)
    if method==5:return [f.diff==f.derivative,f.differentiate==f.derivative]
    args=derivative_arguments(R,form,n)
    if method==6:
        from sage.misc.derivative import multi_derivative
        result=derivative_reference_call(lambda: multi_derivative(f,args))
    else: result=derivative_reference_call(lambda: getattr(f,['derivative','_derivative','diff','differentiate','gradient'][method])(*args))
    frame=lambda g:[list(map(str,g.list())),str(g.parent()),g is f]
    return [frame(g) for g in result] if method==4 else frame(result)
FUNCTIONS['polynomial_derivative_protocol'] = lambda *args: comparison(lambda: polynomial_derivative_protocol(*args))


def derivative_reference_call(run):
    # Installed Sage 10.3 says "non-negative"; bundled misc/derivative.pyx:163
    # says "nonnegative". Normalize only this historical wording difference.
    try: return run()
    except ValueError as e:
        if str(e)=='derivative counts must be non-negative':
            raise ValueError('derivative counts must be nonnegative') from None
        raise

def derivative_parse_protocol(form,n):
    from sage.misc.derivative import derivative_parse
    R=PolynomialRing(QQ,'x')
    args=derivative_arguments(R,form,n)
    out=derivative_reference_call(lambda: derivative_parse(args))
    def frame(x):return [frame(y) for y in x] if isinstance(x,list) else None if x is None else str(x)
    return [frame(out),len(args)==1 and out is args[0]]
FUNCTIONS['derivative_parse_protocol'] = lambda *args: comparison(lambda: derivative_parse_protocol(*args))


def polynomial_pseudo_operand(base,kind,n,coeffs):
    B=PolynomialRing(QQ,'t') if base==11 else unary_polynomial_base(base)
    R=PolynomialRing(B,'x');f=R(coeffs);h=polynomial_input(kind,n)
    result=f.pseudo_quo_rem(h)
    return [[str(g),str(g.parent()),g is f,g is h,[g is k for k in result]] for g in result]
FUNCTIONS['polynomial_pseudo_operand']=lambda *args:comparison(lambda:polynomial_pseudo_operand(*args))

def polynomial_constant_fraction(base,a,b,d,op):
    R=PolynomialRing(unary_polynomial_base(base),'x');K=R.fraction_field()
    f=K(R(a))/K(d);g=K(R(b))/K(d)
    if op==6:return f==g
    h=f+g if op==1 else f-g if op==2 else f*g if op==3 else f/g if op==4 else -f if op==5 else f
    return [str(h),str(h.numerator()),str(h.denominator()),str(h.parent()),h.is_zero(),K.is_field(),K.ring() is R,str(K.zero()),str(K.one()),K(f) is f]
FUNCTIONS['polynomial_constant_fraction']=lambda *args:comparison(lambda:polynomial_constant_fraction(*args))

def polynomial_fraction_protocol(base,op,a,ad,b,bd,flavor,mode,shape=0):
    from sage.rings.fraction_field import FractionField_generic
    B=unary_polynomial_base(base);R=PolynomialRing(B,'x')
    K=FractionField_generic(R) if flavor else R.fraction_field()
    def coefficients(values):return [B.gen()+c if shape and base in (4,9) else QQ(c)/2 if shape and base==0 else c for c in values]
    A,D,E,F=map(lambda v:R(coefficients(v)),[a,ad,b,bd])
    f=K(A,D) if mode==0 else K._element_class(K,A,D,reduce=(mode==1))
    g=K(E,F) if mode==0 else K._element_class(K,E,F,reduce=(mode==1))
    if op==6:return f==g
    # FpTElement.__pow__ reads the leading coefficient of an empty denominator
    # for zero**negative and segfaults (reproduced on the original). Compare
    # this singular input through its checked, original inverse instead.
    if op==10 and type(f).__name__=='FpTElement' and f.numerator().is_zero(): return (~f)**2
    if op==11:f.reduce();h=f
    else:h=f+g if op==1 else f-g if op==2 else f*g if op==3 else f/g if op==4 else -f if op==5 else ~f if op==7 else f**2 if op==8 else f**0 if op==9 else f**(-2) if op==10 else f
    n=h.numerator();d=h.denominator()
    return [str(h),list(map(str,n.list())),list(map(str,d.list())),str(h.parent()),type(h).__name__,h.parent() is K,h is f,h is g,n is h.numerator(),d is h.denominator(),K.zero() is K.zero(),K.one() is K.one(),h.is_zero(),str(K.base_ring()),str(K.characteristic()),K.is_exact(),str(K.gen()),K.ngens(),R.fraction_field() is R.fraction_field(),FractionField(R) is R.fraction_field(),type(h).__name__!='FpTElement' or (h.numer()==n and h.denom()==d)]
FUNCTIONS['polynomial_fraction_protocol']=lambda *args:comparison(lambda:polynomial_fraction_protocol(*args))


def polynomial_fraction_parent(base,flavor,index):
    from sage.rings.fraction_field import FractionField_generic
    from sage.rings.fraction_field_FpT import FpT
    R=PolynomialRing(unary_polynomial_base(base),'x')
    # Bundled FpT.__init__ rejects non-Polynomial_zmod_flint rings after its
    # characteristic range check; installed Sage 10.3 lacks this validation.
    if flavor==2 and 2<R.characteristic()<FpT.INTEGER_LIMIT:
        from sage.rings.polynomial.polynomial_zmod_flint import Polynomial_zmod_flint
        if not issubclass(R.element_class,Polynomial_zmod_flint):raise TypeError('unsupported polynomial ring')
    try: K=FpT(R) if flavor==2 else FractionField_generic(R) if flavor==1 else FractionField(R)
    except TypeError as e:
        # Sage 10.3's factory has a trailing period, removed in bundled source.
        if str(e)=='R must be an integral domain.':raise TypeError('R must be an integral domain') from None
        raise
    return [str(K),type(K).__name__.removesuffix('_with_category'),str(K.gen(index)),K.ngens(),str(K.base_ring()),str(K.characteristic()),K.is_exact()]
FUNCTIONS['polynomial_fraction_parent']=lambda *args:comparison(lambda:polynomial_fraction_parent(*args))


def polynomial_fraction_span(base,n,c):
    R=PolynomialRing(unary_polynomial_base(base),'x');x=R.gen();A=FreeModule(R,3)
    L=A.span([A([x,x**n+c,0]),A([0,0,x])]);V=L.vector_space_span(L.basis())
    return [str(A.base_field()),int(V.dimension()),[[str(e) for e in row] for row in V.basis_matrix().rows()]]
FUNCTIONS['polynomial_fraction_span']=lambda *args:comparison(lambda:polynomial_fraction_span(*args))


def polynomial_fraction_conversion(source,target,mode,operation,a,d,variable=0):
    from sage.rings.fraction_field import FractionField_generic
    B=unary_polynomial_base(source);R=PolynomialRing(B,'x')
    T=R if source==target and variable==0 else PolynomialRing(unary_polynomial_base(target),'y' if variable else 'x')
    K=FractionField_generic(R) if mode in (1,2) else R.fraction_field()
    f=K._element_class(K,R(a),R(d),reduce=(mode in (0,1)))
    n=f.numerator();den=f.denominator()
    try:
        h=T(f) if operation==0 else T.base_ring()(f) if operation==1 else f._conversion(T.base_ring()) if operation==3 else T.fraction_field()(f)
        out=['value',str(h),str(h.parent()),h is n]
    except Exception as e:
        message=str(e)
        if operation in (1,3) and isinstance(e,TypeError) and message=='not a constant polynomial':
            offending=f.numerator() if f.numerator().degree()>0 else f.denominator()
            message=f'{offending} is not a constant polynomial'
        out=['error',type(e).__name__,message]
    return [out,str(f),list(map(str,f.numerator().list())),list(map(str,f.denominator().list())),n is f.numerator(),den is f.denominator()]
FUNCTIONS['polynomial_fraction_conversion']=lambda *args:comparison(lambda:polynomial_fraction_conversion(*args))


def polynomial_fraction_import(mode):
    from importlib import import_module
    file,names=[('fraction_field',['FractionField','FractionField_generic','FractionField_1poly_field']),('fraction_field_element',['FractionFieldElement','FractionFieldElement_1poly_field']),('fraction_field_FpT',['FpT','FpTElement'])][mode]
    module=import_module('sage.rings.'+file)
    return [callable(getattr(module,name)) for name in names]
FUNCTIONS['polynomial_fraction_import']=lambda *args:comparison(lambda:polynomial_fraction_import(*args))

FUNCTIONS['polynomial_ring_coercion']=lambda source,target,variable:comparison(lambda:PolynomialRing(unary_polynomial_base(target),'y' if variable else 'x').has_coerce_map_from(PolynomialRing(unary_polynomial_base(source),'x')))


def polynomial_fraction_tower(base,mode,operation,a,d):
    from sage.rings.fraction_field import FractionField_generic
    U=PolynomialRing(unary_polynomial_base(base),'t')
    F=FractionField_generic(U) if mode in (1,2) else U.fraction_field()
    f=F._element_class(F,U(a),U(d),reduce=(mode in (0,1)))
    n=f.numerator();den=f.denominator()
    # The TypeScript `new` constructor retains its supplied base parent. Use
    # Sage's corresponding non-interned class here: the global factory can
    # otherwise substitute an equal fraction-parent subclass from an earlier case.
    from sage.rings.polynomial.polynomial_ring import PolynomialRing_field, PolynomialRing_integral_domain
    R=PolynomialRing_integral_domain(U,'x') if operation in (1,3,5,7) else PolynomialRing_field(F,'x')
    try:
        if operation<2:h=R(f);result=[str(h),str(h.parent()),list(map(str,h.list()))]
        elif operation==10:
            result=[R(f)==f,R.gen()==f,U.gen()==f,R(f)==R(f),f==U.gen()]
        elif operation>=6:
            z=R.gen() if operation in (6,7) else R(f)
            h=F(z) if operation in (6,8) else U(z)
            result=[str(h),str(h.parent()),h is f]
        else:
            g=R([U.gen()]) if operation<4 else R([U.one(),U.zero(),U.gen()])
            dividend=R.zero() if operation<4 else R.one()
            values=dividend.pseudo_quo_rem(g)
            result=[[str(v),str(v.parent()),[v is w for w in values]] for v in values]
        out=['value',result]
    except Exception as e:
        message=str(e)
        if isinstance(e,TypeError) and message=='not a constant polynomial' and operation in (6,7):
            message=f'{R.gen()} is not a constant polynomial'
        out=['error',type(e).__name__,message]
    return [out,str(f),list(map(str,f.numerator().list())),list(map(str,f.denominator().list())),n is f.numerator(),den is f.denominator()]
FUNCTIONS['polynomial_fraction_tower']=lambda *args:comparison(lambda:polynomial_fraction_tower(*args))


def polynomial_scalar_conversion(source,target,operation,coeffs,denominator):
    R=PolynomialRing(unary_polynomial_base(source),'x');f=R([QQ(c)/denominator for c in coeffs]);T=unary_polynomial_base(target)
    try:
        h=T(f) if operation==0 else f._scalar_conversion(T) if operation==1 else f._integer_(ZZ) if operation==2 else f._rational_()
    except TypeError as e:
        if operation==0 and str(e)=='not a constant polynomial':raise TypeError(f'{f} is not a constant polynomial') from None
        raise
    return [str(h),str(ZZ if operation==2 else QQ if operation==3 else T)]
FUNCTIONS['polynomial_scalar_conversion']=lambda *args:comparison(lambda:polynomial_scalar_conversion(*args))
FUNCTIONS['gf2_scalar_hook']=lambda value,operation:comparison(lambda:str(GF(2)(value)._integer_() if operation==0 else GF(2)(value)._rational_()))


def bundled_plain_power(f,e):
    h=f**e
    # Sage 10.3 inherits Polynomial.__pow__; bundled Polynomial_dense_mod_p
    # adds an NTL wrapper that reconstructs a new polynomial even for exponent 1.
    # Execute that native path, rather than replacing an expected identity bit.
    if type(f).__name__=='Polynomial_dense_mod_p' and f.degree()>0 and e==1:
        return f.parent()(f.ntl_ZZ_pX()**ZZ(e),construct=True)
    return h

def polynomial_power_integer(base,kind,exponent,coeffs):
    from sage.rings.polynomial.polynomial_element import Polynomial
    R=PolynomialRing(unary_polynomial_base(base),'x');f=R(coeffs)
    e=int(exponent) if kind==1 else ZZ(exponent)
    if exponent == -(2**63) and base in (5,8,10):
        # Signed C-long negation overflows. The original short-polynomial NTL
        # wrapper omits sig_on and aborts; add the same signal guard used by
        # its long-polynomial branch so the original NTL error is catchable.
        from pathlib import Path
        from sage.misc.cython import cython_import
        global _safe_power_module
        if '_safe_power_module' not in globals():
            _safe_power_module=cython_import(str(Path(__file__).parent.parent/'polynomial_power_native.pyx'))
        h=_safe_power_module.safe_power(f,e)
    else:h=bundled_plain_power(f,e)
    return [str(h),str(h.parent()),h is f,['polynomial',list(map(str,h.list()))] if isinstance(h,Polynomial) else ['fraction',list(map(str,h.numerator().list())),list(map(str,h.denominator().list()))]]
FUNCTIONS['polynomial_power_integer']=lambda *args:comparison(lambda:polynomial_power_integer(*args))

def polynomial_power_strict_rational(coeffs,denominator,exponent):
    R=PolynomialRing(QQ,'x');h=R([QQ(c)/denominator for c in coeffs])**exponent
    return [str(h),list(map(str,h.list())),all(c.parent() is QQ for c in h.list())]
FUNCTIONS['polynomial_power_strict_rational']=lambda *args:comparison(lambda:polynomial_power_strict_rational(*args))

def polynomial_power_exponent(base,kind,numerator,denominator,coeffs):
    from sage.rings.polynomial.polynomial_element import Polynomial
    R=PolynomialRing(unary_polynomial_base(base),'x');f=R(coeffs)
    exponent=numerator if kind==0 else (int(numerator//denominator) if numerator%denominator==0 else float(numerator)/float(denominator)) if kind==1 else QQ(numerator)/denominator if kind==2 else bool(numerator) if kind==3 else str(numerator) if kind==4 else None if kind==5 else [numerator] if kind==6 else R(numerator) if kind==7 else GF(7)(numerator) if kind==8 else GF(2)(numerator) if kind==9 else RDF(numerator)/denominator if kind==10 else float('nan') if kind==11 else float('inf') if kind==12 else float('-inf') if kind==13 else R([numerator,denominator])
    h=bundled_plain_power(f,exponent)
    return [str(h),str(h.parent()),h is f,['polynomial',list(map(str,h.list()))] if isinstance(h,Polynomial) else ['fraction',list(map(str,h.numerator().list())),list(map(str,h.denominator().list()))]]
FUNCTIONS['polynomial_power_exponent']=lambda *args:comparison(lambda:polynomial_power_exponent(*args))

def polynomial_root_series(base,operation,coeffs,denominator,n,precision,other,start_kind):
    from sage.rings.polynomial.polynomial_element import generic_power_trunc
    R=PolynomialRing(unary_polynomial_base(base),'x');f=R([QQ(c)/denominator for c in coeffs]);g=R(other)
    h=generic_power_trunc(f,ZZ(n),precision) if operation==7 else f.nth_root(n) if operation==0 else f._nth_root_series(n,precision,None if start_kind==0 else R.one() if start_kind==1 else g) if operation==1 else f.inverse_series_trunc(precision) if operation==2 else f.power_trunc(n,precision) if operation==3 else f._power_trunc(n,precision) if operation==4 else f._mul_trunc_(g,precision) if operation==5 else f.multiplication_trunc(g,precision)
    return [str(h),str(h.parent()),h is f,h is g,list(map(str,h.list()))]
def polynomial_root_series_oracle(*args):
    try:return polynomial_root_series(*args)
    except ValueError as e:
        # Sage 10.3 spells this with a hyphen; bundled generic_power_trunc does not.
        if str(e)=='n must be a non-negative integer':raise ValueError('n must be a nonnegative integer') from None
        raise
FUNCTIONS['polynomial_root_series']=lambda *args:comparison(lambda:polynomial_root_series_oracle(*args))

def finite_coefficient_nth_root(kind,p,degree,n,value):
    B=Zmod(p) if kind==0 else GF(p) if kind in (1,2) else GF(2) if kind==3 else GF(p**degree,'a')
    c=B(list(ZZ(value).digits(p))) if kind==4 else B(value)
    return [str(c.nth_root(n)),str(B)]
FUNCTIONS['finite_coefficient_nth_root']=lambda *args:comparison(lambda:finite_coefficient_nth_root(*args))

def polynomial_factor_order(base, factors, multiplicities):
    R=PolynomialRing(unary_polynomial_base(base),'x');f=R.one()
    for j,m in enumerate(multiplicities):
        f*=R(factors[3*j:3*j+3])**m
    F=f.factor()
    answer=[[str(g),str(e)] for g,e in F]
    if F.unit()!=1: answer.insert(0,[str(F.unit()),'1'])
    return answer
FUNCTIONS['polynomial_factor_order']=lambda *args:comparison(lambda:polynomial_factor_order(*args))

def polynomial_series_input(base,operation,slot,kind,value,coeffs):
    R=PolynomialRing(unary_polynomial_base(base),'x');f=R(coeffs)
    v=polynomial_input(kind,value);n=v if slot==0 else 2;precision=v if slot==1 else 3
    h=f.nth_root(n) if operation==0 else f._nth_root_series(n,precision) if operation==1 else f.inverse_series_trunc(precision) if operation==2 else f.power_trunc(n,precision) if operation==3 else f._power_trunc(n,precision) if operation==4 else f._mul_trunc_(f,precision) if operation==5 else f.multiplication_trunc(f,precision)
    return [str(h),str(h.parent()),h is f]
FUNCTIONS['polynomial_series_input']=lambda *args:comparison(lambda:polynomial_series_input(*args))


def polynomial_truncated_operand(base,kind,value,coeffs,precision,internal,precision_kind=-1):
    R=PolynomialRing(RDF if base==11 else unary_polynomial_base(base),'x');f=R(coeffs)
    other=RDF(value)/2 if kind==27 else R(value) if kind==28 else R.fraction_field()(1)/(R.gen()+1) if kind==29 else R([value,1]) if kind==30 else polynomial_input(kind,value)
    precision=precision if precision_kind==-1 else polynomial_input(precision_kind,precision)
    g=f._mul_trunc_(other,precision) if internal else f.multiplication_trunc(other,precision)
    return [str(g),str(g.parent()),g is f,g is other]
FUNCTIONS['polynomial_truncated_operand']=lambda *args:comparison(lambda:polynomial_truncated_operand(*args))

def polynomial_real_truncated_product(a,b,n,public):
    import struct,math
    def unpack(x):return RDF(struct.unpack('>d',int(x).to_bytes(8,'big'))[0])
    def pack(x):return 'nan' if math.isnan(float(x)) else struct.pack('>d',float(x)).hex()
    R=PolynomialRing(RDF,'x');f=R(list(map(unpack,a)));g=f if public in (3,4,5) else R(list(map(unpack,b)))
    h=f*g if public in (2,3) else f.multiplication_trunc(g,n) if public in (1,5) else f._mul_trunc_(g,n)
    return [list(map(pack,h.list())),str(h)]
FUNCTIONS['polynomial_real_truncated_product']=lambda *args:comparison(lambda:polynomial_real_truncated_product(*args))

def polynomial_exact_product(base,a,b,same,shape,denominator):
    import hashlib
    B=PolynomialRing(QQ,'t') if base==11 else PolynomialRing(QQ,'t').fraction_field() if base==12 else unary_polynomial_base(base)
    R=PolynomialRing(B,'x')
    def coefficient(c):
        if base in (4,9) and shape:return B(list(ZZ(c).digits(B.characteristic())))
        if base in (11,12):return (B.gen()+QQ(c)/denominator) if shape else B(QQ(c)/denominator)
        return QQ(c)/denominator if base==0 else c
    f=R(list(map(coefficient,a)));g=f if same else R(list(map(coefficient,b)));h=f*g
    data=','.join(map(str,h.list())).encode()
    return [len(h.list()),hashlib.sha256(data).hexdigest(),str(h.parent()),h is f,h is g]
FUNCTIONS['polynomial_exact_product']=lambda *args:comparison(lambda:polynomial_exact_product(*args))

def polynomial_matrix_product(a,b,same):
    M=MatrixSpace(ZZ,2);R=PolynomialRing(M,'x')
    def poly(v):return R([M(v[i:i+4]) for i in range(0,len(v),4)])
    f=poly(a);g=f if same else poly(b);h=f*g
    return [[str(c) for c in m.list()] for m in h.list()]
FUNCTIONS['polynomial_matrix_product']=lambda *args:comparison(lambda:polynomial_matrix_product(*args))

def polynomial_modular_power(base,a,kind,value,m,exponent_kind,exponent):
    B=PolynomialRing(QQ,'t') if base==11 else RDF if base==12 else unary_polynomial_base(base)
    R=PolynomialRing(B,'x');f=R(a)
    modulus=PolynomialRing(QQ,'x')(m) if kind==29 else PolynomialRing(GF(5),'x')(m) if kind==30 else PolynomialRing(B,'y')(m) if kind==31 else PolynomialRing(PolynomialRing(QQ,'t'),'x')(m) if kind==32 else R(m) if kind==0 else f if kind==1 else polynomial_input(kind-2,value) if kind>=2 else None
    e=polynomial_input(exponent_kind,exponent)
    h=bundled_modular_power(f,e,modulus)
    return [str(h),str(h.parent()),h is f,h is modulus]


def bundled_modular_power(f,e,modulus):
    # Only wrapper-level differences from installed Sage 10.3 are transcribed;
    # arithmetic continues through the original native implementation.
    if modulus is None:return bundled_plain_power(f,e)
    name=type(f).__name__
    if name in ('Polynomial_integer_dense_flint','Polynomial_rational_flint'):
        raise NotImplementedError('pow() with a modulus is not implemented for this ring')
    extension=name=='Polynomial_ZZ_pEX'
    word=name=='Polynomial_zmod_flint'
    if extension or word:
        e=ZZ(e)
        if extension and modulus.is_zero():raise ZeroDivisionError('modulus must be nonzero')
        if f.parent() is not parent(modulus):
            from sage.structure.element import get_coercion_model
            a,m=get_coercion_model().canonical_coercion(f,modulus)
            if a is not f:return bundled_modular_power(a,e,m)
            modulus=m
        f=f % modulus
    if extension and e>0 and e.nbits()>=32 and modulus.degree()<=0:
        # New wrapper enters the native big-exponent path before the template's
        # zero shortcut. Guard only the NTL build, not PARI-backed arithmetic.
        from .polynomial_ops import native_modular_module
        B=f.base_ring();p=B.characteristic()
        context=ntl.ZZ_pEContext(ntl.ZZ_pX(list(B.modulus()),p))
        def poly(v):return ntl.ZZ_pEX([ntl.ZZ_pE(list(c.polynomial()),context) for c in v.list()],context)
        return native_modular_module().native_extension_modular_power(poly(f),e,poly(modulus),0)
    if name in ('Polynomial_ZZ_pEX','Polynomial_GF2X') and f and -2**63<=e<2**63 and not (extension and e>0 and e.nbits()>=32):
        if modulus.is_zero():raise ZeroDivisionError('modulus must be nonzero')
        if modulus.is_one():return f.parent().zero()
    if name=='Polynomial_GF2X' and e>=2**63 and f.is_gen() and not modulus:
        # Sage 10.3's overflowing list repetition damages cached GF(2) values
        # after repeated calls. Execute the original once in an isolated process.
        import subprocess,sys,json,builtins
        script="""from sage.all import *
import json
R=PolynomialRing(GF(2),'x')
try: pow(R.gen(),ZZ(%r),R.zero())
except Exception as err: print(json.dumps([type(err).__name__,str(err)]))
""" % str(e)
        result=subprocess.run([sys.executable,'-c',script],capture_output=True,text=True,check=True)
        kind,message=json.loads(result.stdout)
        raise getattr(builtins,kind)(message)
    if name=='Polynomial_GF2X' and e == -2**63 and f:
        # The template's C-long negation wraps at LONG_MIN. Run just the
        # original native kernel under a signal guard, then construct/invert
        # outside it to avoid Sage 10.3's error-path coefficient-cache damage.
        from .polynomial_ops import native_modular_module
        h=native_modular_module().native_binary_modular_power(ntl.GF2X(list(f)),e,ntl.GF2X(list(modulus)))
        return ~f.parent()([ZZ(str(c)) for c in h.list()])
    return pow(f,e,modulus)
FUNCTIONS['polynomial_modular_power']=lambda *args:comparison(lambda:polynomial_modular_power(*args))

def polynomial_real_integer_power(bits,n,method):
    import struct,math
    x=RDF(struct.unpack('>d',int(bits).to_bytes(8,'big'))[0])
    if method>=3:
        from .polynomial_ops import native_power_module
        h=native_power_module().native_gsl_power_value(x,int(n),int(method-3))
    elif method==0:h=x**n
    elif method==1:h=~x
    else:h=PolynomialRing(RDF,'x')([x])**n;h=h[0]
    value=float(h)
    return 'nan' if math.isnan(value) else struct.pack('>d',value).hex()
FUNCTIONS['polynomial_real_integer_power']=lambda *args:comparison(lambda:polynomial_real_integer_power(*args))

def polynomial_real_quotient(a,b,same):
    import struct,math
    def unpack(x):return RDF(struct.unpack('>d',int(x).to_bytes(8,'big'))[0])
    def pack(x):return 'nan' if math.isnan(float(x)) else struct.pack('>d',float(x)).hex()
    R=PolynomialRing(RDF,'x');f=R(list(map(unpack,a)));g=f if same else R(list(map(unpack,b)))
    return [[list(map(pack,h.list())),str(h)] for h in f.quo_rem(g)]
FUNCTIONS['polynomial_real_quotient']=lambda *args:comparison(lambda:polynomial_real_quotient(*args))

def polynomial_real_evaluate(a,point):
    import struct,math
    def unpack(v):return RDF(struct.unpack('>d',int(v).to_bytes(8,'big'))[0])
    value=float(PolynomialRing(RDF,'x')(list(map(unpack,a)))(unpack(point)))
    return 'nan' if math.isnan(value) else struct.pack('>d',value).hex()
FUNCTIONS['polynomial_real_evaluate']=lambda *args:comparison(lambda:polynomial_real_evaluate(*args))

def polynomial_compiled_evaluation(a,point,algorithm,method):
    import struct,math
    from sage.rings.polynomial.polynomial_compiled import CompiledPolynomialFunction
    def unpack(v):return RDF(struct.unpack('>d',int(v).to_bytes(8,'big'))[0])
    def pack(v):return 'nan' if math.isnan(float(v)) else struct.pack('>d',float(v)).hex()
    f=CompiledPolynomialFunction(list(map(unpack,a)),['binary','pippenger','other',''][int(algorithm)])
    x=unpack(point)
    return str(f) if method==1 else [pack(f(x)),pack(f(-x)),pack(f(x))] if method==2 else pack(f(x))
FUNCTIONS['polynomial_compiled_evaluation']=lambda *args:comparison(lambda:polynomial_compiled_evaluation(*args))

def bundled_polynomial_evaluation(f,x):
    h=f(x)
    # Installed Sage 10.3's large-composite NTL scalar wrapper constructs its
    # result in the polynomial parent. Bundled __call__ constructs in _base.
    if type(f).__name__=='Polynomial_dense_modn_ntl_ZZ' and f.base_ring().has_coerce_map_from(parent(x)):
        return f.base_ring()(h)
    return h

def polynomial_evaluation_input(base,a,kind,n):
    from sage.rings.polynomial.polynomial_element import Polynomial
    B=PolynomialRing(QQ,'t') if base==11 else RDF if base==12 else unary_polynomial_base(base)
    R=PolynomialRing(B,'x');f=R(a)
    x=B(n) if kind==-1 else RDF(n) if kind==27 else RDF('NaN') if kind==28 else RDF('-infinity' if n<0 else '+infinity') if kind==29 else f if kind==30 else R.gen() if kind==31 else R.zero() if kind==32 else PolynomialRing(QQ,'y')([n,1]) if kind==33 else PolynomialRing(B,'y')([n,1]) if kind==34 else R([n,1]) if kind==35 else R([n,0,1]) if kind==36 else polynomial_input(kind,n)
    h=bundled_polynomial_evaluation(f,x)
    import struct,math
    value=('float:' + ('nan' if math.isnan(h) else struct.pack('>d',h).hex())) if isinstance(h,float) else str(h)
    return [value,str(parent(h)),h is f,h is x if isinstance(h,Polynomial) else False]
FUNCTIONS['polynomial_evaluation_input']=lambda *args:comparison(lambda:polynomial_evaluation_input(*args))

def polynomial_generator_identity(base,kind):
    B=PolynomialRing(QQ,'t') if base==11 else RDF if base==12 else unary_polynomial_base(base)
    R=PolynomialRing(B,'x');x=R.gen()
    values=[x,R([0,1]),R(x),x+R.zero(),x*R.one(),-(-x),bundled_plain_power(x,1),x**0,R.zero(),R.one(),x(x),R([0,1])(x),x(R([0,1]))]
    f=values[int(kind)]
    return [f.is_gen(),f is x,x is R.gen(),str(f)]
FUNCTIONS['polynomial_generator_identity']=lambda *args:comparison(lambda:polynomial_generator_identity(*args))

def polynomial_composition_input(base,a,b,target,method):
    import struct,math
    B=PolynomialRing(QQ,'t') if base==11 else RDF if base==12 else unary_polynomial_base(base)
    R=PolynomialRing(B,'x');S=R if target==0 else PolynomialRing(B if target==1 else QQ if target==2 else RDF,'y' if target==1 else 'x')
    f=R(a);g=S(b);h=f(g)
    def pack(c):return ('nan' if math.isnan(c) else struct.pack('>d',float(c)).hex()) if parent(c) is RDF else str(c)
    return [str(h),str(parent(h)),h is f,h is g,list(map(pack,h.list()))]
FUNCTIONS['polynomial_composition_input']=lambda *args:comparison(lambda:polynomial_composition_input(*args))

def polynomial_evaluation_stress(base,length,pattern,den,num,xden,kind):
    import struct,math
    B=RDF if base==12 else unary_polynomial_base(base);N=int(length)
    coeffs=[]
    for i in range(N):
        n=0 if (pattern==1 and i!=0 and i!=N-1) or (pattern==2 and i%13!=0) else ((i%7)-3)*(2**300 if pattern==3 else 1)
        coeffs.append(B(n if base==6 else QQ(n)/den))
    f=PolynomialRing(B,'x')(coeffs);q=QQ(num)/xden
    x=q if kind==0 else RDF(q) if kind==1 else RDF('-infinity' if num<0 else '+infinity') if kind==2 else float(q) if kind==3 else float('nan')
    def pack(h):
        if isinstance(h,float) or parent(h) is RDF:return ['nan' if math.isnan(h) else struct.pack('>d',float(h)).hex(),str(parent(h))]
        return [str(h),str(parent(h))]
    return [pack(f(x)),pack(f(x))]
FUNCTIONS['polynomial_evaluation_stress']=lambda *args:comparison(lambda:polynomial_evaluation_stress(*args))

def polynomial_generator_constructor(base,a,flag):
    B=PolynomialRing(QQ,'t') if base==11 else RDF if base==12 else unary_polynomial_base(base)
    R=PolynomialRing(B,'x');f=R.element_class(R,a,is_gen=bool(flag))
    return [str(f),f.degree(),f.is_gen(),list(map(str,f.list())),str(bundled_polynomial_evaluation(f,ZZ(2)))]
FUNCTIONS['polynomial_generator_constructor']=lambda *args:comparison(lambda:polynomial_generator_constructor(*args))

def polynomial_matrix_evaluation(base,target,rows,cols,pattern,a):
    B=PolynomialRing(QQ,'t') if base==11 else RDF if base==12 else unary_polynomial_base(base)
    T=ZZ if target==0 else QQ if target==1 else RDF if target==2 else Zmod(7) if target==3 else Zmod(14) if target==4 else GF(7) if target==5 else unary_polynomial_base(4) if target==6 else GF(2) if target in (7,9) else Zmod(1)
    m=int(rows);n=int(cols)
    entries=[0 if pattern==0 else (1 if i==j else 0) if pattern==1 else (i*n+j)%7-2 for i in range(m) for j in range(n)]
    # GenericMatrix mirrors matrix_generic_dense.pyx, so request that actual
    # original backend; IntegerMatrix and Matrix_modn_dense use native defaults.
    point=MatrixSpace(T,m,n,implementation=None if target in (0,3,4,9) else 'generic')(entries);f=PolynomialRing(B,'x')(a);h=bundled_matrix_evaluation(f,point)
    if type(h).__name__=='Matrix_gf2e_dense' and not hasattr(h.base_ring(),'_cache'):
        # Serialize the original native value without its unrelated Givaro-only
        # getter. Evaluation itself already returned this matrix successfully.
        values=[[h.base_ring()(ZZ(c).digits(2)) for c in row] for row in native_matrix_evaluation_module().native_binary_matrix_rows(h)]
    else:values=h.rows()
    return [[list(map(str,row)) for row in values],str(h.base_ring()),h.nrows(),h.ncols(),h is point]
FUNCTIONS['polynomial_matrix_evaluation']=lambda *args:comparison(lambda:polynomial_matrix_evaluation(*args))

_bundled_matrix_is_exact=None
def bundled_matrix_evaluation(f,point):
    # Sage 10.3 inherits Parent.is_exact() == True. Execute the bundled
    # MatrixSpace method verbatim; it delegates exactness to the base ring.
    from sage.matrix.matrix_space import MatrixSpace as MatrixSpaceClass
    if 'is_exact' in MatrixSpaceClass.__dict__:return f(point)
    global _bundled_matrix_is_exact
    if _bundled_matrix_is_exact is None:
        import ast
        from pathlib import Path
        source=Path(__file__).parents[4]/'reference/sage/src/sage/matrix/matrix_space.py'
        cls=next(n for n in ast.parse(source.read_text()).body if isinstance(n,ast.ClassDef) and n.name=='MatrixSpace')
        method=next(n for n in cls.body if isinstance(n,ast.FunctionDef) and n.name=='is_exact')
        namespace={}
        exec(compile(ast.Module(body=[method],type_ignores=[]),str(source),'exec'),namespace)
        _bundled_matrix_is_exact=namespace['is_exact']
    MatrixSpaceClass.is_exact=_bundled_matrix_is_exact
    try:return f(point)
    finally:del MatrixSpaceClass.is_exact

_native_matrix_evaluation_module=None
_native_matrix_evaluation_error=None
def native_matrix_evaluation_module():
    global _native_matrix_evaluation_module,_native_matrix_evaluation_error
    if _native_matrix_evaluation_error is not None:raise _native_matrix_evaluation_error
    if _native_matrix_evaluation_module is None:
        import os,sys
        from pathlib import Path
        from sage.misc.cython import cython_import
        sys.stdout.flush();saved=os.dup(1)
        try:
            os.dup2(2,1)
            _native_matrix_evaluation_module=cython_import(str(Path(__file__).parent.parent/'polynomial_matrix_native.pyx'))
        except Exception as error:
            _native_matrix_evaluation_error=error
            raise
        finally:sys.stdout.flush();os.dup2(saved,1);os.close(saved)
    return _native_matrix_evaluation_module

def polynomial_matrix_actions(base,target,n,wrap,degree,den,a,b):
    import struct,math
    def real(v):return RDF(struct.unpack('>d',int(v).to_bytes(8,'big'))[0])
    B=PolynomialRing(QQ,'t') if base==11 else RDF if base==12 else unary_polynomial_base(base)
    T=ZZ if target==0 else QQ if target==1 else RDF
    entries=[real(b[i%len(b)]) if target==2 else T(b[i%len(b)])/den if target==1 else T(b[i%len(b)]) for i in range(n*n)]
    point=MatrixSpace(T,n,implementation=None if target==0 else 'generic')(entries)
    cs={int(i*degree//max(1,len(a)-1)):real(v) if base==12 else B([QQ(v)/den,QQ(i+1)/den]) if base==11 else B(QQ(v)/den) for i,v in enumerate(a)}
    f=PolynomialRing(B,'x')(cs)
    def frame(h):
        return [[[('NaN' if math.isnan(float(c)) else str(int.from_bytes(struct.pack('>d',float(c)),'big'))) if h.base_ring() is RDF else str(c) for c in row] for row in h.rows()],str(h.base_ring()),h.nrows(),h.ncols(),h is point]
    first=frame(bundled_matrix_evaluation(f,[point] if wrap else point))
    # Reuse the same compiled plan across scalar and matrix calls.
    f(B(2))
    second=frame(bundled_matrix_evaluation(f,[point] if wrap else point))
    return [first,second]
FUNCTIONS['polynomial_matrix_actions']=lambda *args:comparison(lambda:polynomial_matrix_actions(*args))

def polynomial_matrix_double_conversion(rows,cols,a):
    import struct
    from sage.matrix.change_ring import integer_to_real_double_dense
    h=integer_to_real_double_dense(matrix(ZZ,rows,cols,[a[i%len(a)] for i in range(rows*cols)]))
    return [[[str(int.from_bytes(struct.pack('>d',float(c)),'big')) for c in row] for row in h.rows()],str(h.base_ring()),h.nrows(),h.ncols()]
FUNCTIONS['polynomial_matrix_double_conversion']=lambda *args:comparison(lambda:polynomial_matrix_double_conversion(*args))

from real_double_group import rdf_scalar
FUNCTIONS["rdf_scalar"] = rdf_scalar

from real_double_group import rdf_alias
FUNCTIONS["rdf_alias"] = rdf_alias

from flint_remainder import native_rem
FUNCTIONS.update(flint_remainder=native_rem)

# Preserve bundled root order, multiplicities, and ring-specific zero errors.
from polynomial_roots import polynomial_roots_fidelity
FUNCTIONS['polynomial_roots_fidelity'] = lambda *args: comparison(lambda: polynomial_roots_fidelity(*args))

from polynomial_factor_dispatch import polynomial_factor_fidelity
FUNCTIONS['polynomial_factor_fidelity'] = lambda *args: comparison(lambda: polynomial_factor_fidelity(*args))
from polynomial_factor_dispatch import polynomial_factor_state
FUNCTIONS['polynomial_factor_state'] = polynomial_factor_state

from finite_polynomial_roots import finite_polynomial_roots
FUNCTIONS['finite_polynomial_roots'] = lambda *args: comparison(lambda: finite_polynomial_roots(*args))

from pari_ffelt_compare import pari_ffelt_compare
FUNCTIONS['pari_ffelt_compare'] = pari_ffelt_compare

FUNCTIONS['finite_polynomial_factors'] = lambda *args: comparison(lambda: finite_polynomial_roots(*args, factor_only=True))

FUNCTIONS['pari_ffelt_compare_scalar'] = pari_ffelt_compare

from finite_polynomial_roots import finite_polynomial_roots_state
FUNCTIONS['finite_polynomial_roots_state'] = lambda *args: comparison(lambda: finite_polynomial_roots_state(*args))


FUNCTIONS['finite_polynomial_distinct_roots'] = lambda *args: comparison(lambda: finite_polynomial_roots(*args, distinct=True))
from polynomial_roots import polynomial_distinct_roots_fidelity
FUNCTIONS['polynomial_distinct_roots_fidelity'] = lambda *args: comparison(lambda: polynomial_distinct_roots_fidelity(*args))

from finite_polynomial_roots import finite_polynomial_distinct_roots_state
FUNCTIONS['finite_polynomial_distinct_roots_state'] = lambda *args: comparison(lambda: finite_polynomial_distinct_roots_state(*args))
