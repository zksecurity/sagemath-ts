"""Live-SageMath oracle for rational function fields over prime fields."""

from sage.all import *


def _fraction_from_strings(F,n,d=None):
    # Fraction-field fallback calls sage_eval. The bundled preparser removed the
    # old backslash operator, so execute that exact source instead of Sage 10.3's.
    import importlib,importlib.util
    from pathlib import Path
    global _audit_preparser
    if '_audit_preparser' not in globals():
        path=Path(__file__).resolve().parents[4]/'reference/sage/src/sage/repl/preparse.py'
        spec=importlib.util.spec_from_file_location('_audit_preparser',path)
        _audit_preparser=importlib.util.module_from_spec(spec);spec.loader.exec_module(_audit_preparser)
    module=importlib.import_module('sage.misc.sage_eval');old=module.preparser
    try:
        module.preparser=_audit_preparser
        return F(n)if d is None else F(n,d)
    finally:module.preparser=old



def _ints(xs):
    return ','.join(str(x) for x in xs)


def _element(K, coeffs):
    return K(K._ring(coeffs))


def _parts(f):
    return '%s/%s' % (_ints(f.numerator().list()), _ints(f.denominator().list()))


def ff_arithmetic(p, numerator, denominator, exponent):
    K = FunctionField(QQ if p==0 else GF(ZZ(p)), 'x')
    f = _element(K, numerator) / _element(K, denominator)
    g = f ** ZZ(exponent)
    return 'f=%s pow=%s degree=%s vx=%s square=%s' % (
        _parts(f),
        _parts(g),
        f.degree(),
        f.valuation(K.gen()),
        f.is_square(),
    )


def ff_factor(p, numerator, denominator):
    K = FunctionField(QQ if p==0 else GF(ZZ(p)), 'x')
    f = _element(K, numerator) / _element(K, denominator)
    # Bundled Sage (polynomial_zmod_flint.pyx:820 and
    # polynomial_element.pyx:5045) rejects zero before factoring. Installed
    # Sage 10.3 omits the guard for its word-modular backend.
    if f.is_zero():
        raise ArithmeticError('factorization of 0 is not defined')
    F = f.factor()
    factors = ';'.join('%s:%s' % (_parts(q), e) for q, e in F)
    return 'unit=%s factors=%s' % (_parts(K(F.unit())), factors)


def ff_places(p, degree):
    K = FunctionField(QQ if p==0 else GF(ZZ(p)), 'x')
    places = K.places(int(degree))
    return ';'.join(str(P) for P in places)


FUNCTIONS = {
    'ff_arithmetic': ff_arithmetic,
    'ff_factor': ff_factor,
    'ff_places': ff_places,
}


def ff_parent(p,op,n):
    K=FunctionField(QQ if p==0 else GF(p),'x');k=K.constant_base_field()
    if op==0:return '|'.join(map(str,[K,K.variable_name(),','.join(K.variable_names()),K.ngens(),K.degree(),K.genus(),K.characteristic(),K.is_finite(),K.is_global(),K.is_perfect(),K.is_field(),K.base_field() is K,K.rational_function_field() is K,K.constant_field() is k,K.constant_base_field() is k]))
    if op==1:return _parts(K.gen(n))
    if op==2:return str(K.degree(None if n==0 else K if n==1 else FunctionField(k,'x' if n==2 else 'y')))
    if op==3:return str(K is FunctionField(k,'x' if n==0 else ('x',)))
    if op==4:
        L,from_L,to_L=K.change_variable_name('x' if n==0 else 'y')
        return '|'.join(map(str,[L,L.change_variable_name('x')[0] is K,_parts(from_L(L.gen())),_parts(to_L(K.gen())),from_L(to_L(K.gen()))==K.gen()]))
    if op==5:return str(K.different())
    if op==6:return _parts(K.maximal_order().gen(n))
    raise ValueError('unknown ff_parent operation')


def ff_element(p,op,num,den,num2,den2,n):
    K=FunctionField(QQ if p==0 else GF(p),'x')
    f=_element(K,num)/_element(K,den);g=_element(K,num2)/_element(K,den2)
    if op==0:return '|'.join(map(str,[_parts(f),f.is_zero(),f.is_one(),f.degree(),_parts(f.trace()),_parts(f.norm()),_parts(f.matrix()[0,0]),_parts(f.list()[0]),_parts(f.element()),f.parent() is K]))
    if op==1:return _parts(-f)
    if op==2:return _parts(~f)
    if op==3:return _parts(f**n)
    if op==4:return _parts(f.sqrt())
    if op==5:return ';'.join(_parts(x)for x in f.sqrt(all=True))
    if op==6:return str(f.is_square())
    if op==7:return str((f>g)-(f<g))
    if op==8:return _parts(f+g)
    if op==9:return _parts(f-g)
    if op==10:return _parts(f*g)
    if op==11:return _parts(f/g)
    if op==12:return _parts(f.nth_root(n))
    if op==13:return str(f.is_nth_power(n))
    if op==14:return str(f)
    if op==15:return _parts(f*K.constant_base_field()(n))
    raise ValueError('unknown ff_element operation')
FUNCTIONS.update(ff_parent=ff_parent,ff_element=ff_element)


def ff_valuation(p,num,den,vnum,vden,mode):
    K=FunctionField(QQ if p==0 else GF(p),'x')
    f=_element(K,num)/_element(K,den);v=_element(K,vnum)/_element(K,vden)
    return str(f.valuation(v if mode==0 else v.numerator()))
FUNCTIONS['ff_valuation']=ff_valuation


def ff_zero_ideal(p,infinite,op):
    K=FunctionField(QQ if p==0 else GF(p),'x')
    O=K.maximal_order_infinite() if infinite else K.maximal_order();I=O.ideal(K.zero())
    if op==0:return str(I.is_zero())
    # The finite ideal's _factor calls its generator's factor(), reaching the
    # bundled polynomial zero guard; the infinite ideal uses degree arithmetic.
    if not infinite:raise ArithmeticError('factorization of 0 is not defined')
    if op==1:return ';'.join('%s:%s'%(_parts(P.gens()[0]),e)for P,e in I.factor())
    return str(I.divisor() if op==2 else I.divisor_of_zeros() if op==3 else I.divisor_of_poles())
FUNCTIONS['ff_zero_ideal']=ff_zero_ideal


def ff_bridge(p,q,op,mode,num,den,codes):
    K=FunctionField(QQ if p==0 else GF(p),'x')
    if op==0:return str(K.field())
    if op in (5,6):
        f=_fraction_from_strings(K,''.join(chr(int(c))for c in codes))
        if op==6:return '|'.join(map(str,[f.element().parent() is K.field(),f.numerator().parent() is K._ring,f.denominator().parent() is K._ring]))
        return _parts(f)+'|'+str(f)
    f=_element(K,num)/_element(K,den)
    if op==1:
        value=f.element();return '|'.join(map(str,[value.parent(),value.parent() is K.field(),value is f.element(),type(value).__name__]))
    if op==2:
        value=K(f.element());return _parts(value)+'|'+str(value)
    if op==3:
        R=K._ring;value=K(R.fraction_field()(R(num),R(den)));return _parts(value)+'|'+str(value)
    L=FunctionField(K.constant_base_field()if q==p else QQ if q==0 else GF(q),'x'if mode==0 else'y')
    value=K(_element(L,num)/_element(L,den));return _parts(value)+'|'+str(value)
FUNCTIONS['ff_bridge']=ff_bridge


def ff_fraction_strings(p,mode,ncodes,dcodes):
    K=FunctionField(QQ if p==0 else GF(p),'x');F=K.field()
    n=''.join(chr(int(c))for c in ncodes);d=''.join(chr(int(c))for c in dcodes)
    f=_fraction_from_strings(F,n,None if mode==0 else d)
    return _parts(f)+'|'+str(f)
FUNCTIONS['ff_fraction_strings']=ff_fraction_strings


def ff_wrapped(p,op,num,den):
    K=FunctionField(QQ if p==0 else GF(p),'x');F=K.field();R=K._ring
    raw=F._element_class(F,R(num),R(den),reduce=False);f=K.element_class(K,raw)
    before=_parts(f)
    if op==0:return '|'.join(map(str,[before,f.is_zero(),f.is_one(),f==K.zero(),f==K.one(),str(f),f.element() is raw]))
    if op==1:
        try:R(raw);err='ok'
        except Exception as e:err=type(e).__name__+':'+str(e)
        return before+'|'+err+'|'+_parts(f)+'|'+str(f)
    if op==2:
        if hasattr(raw,'reduce'):raw.reduce()
        return before+'|'+_parts(f)+'|'+str(f)
    g=-f if op==3 else f+K.one()if op==4 else f-K.one()if op==5 else f*K.gen()if op==6 else f/K.gen()if op==7 else f**(op-8)
    return _parts(g)+'|'+str(g)
FUNCTIONS['ff_wrapped']=ff_wrapped


def ff_inverse(p,mode,num,den,inum,iden):
    K=FunctionField(QQ if p==0 else GF(p),'x');f=_element(K,num)/_element(K,den)
    g=_element(K,inum)/_element(K,iden)
    if mode in (4,5):
        F=K.field();R=K._ring
        raw=F._element_class(F,R(num if mode==4 else inum),R(den if mode==4 else iden),reduce=False)
        if mode==4:f=K.element_class(K,raw)
        else:g=K.element_class(K,raw)
    if mode<2:I=(K.maximal_order()if mode==0 else K.maximal_order_infinite()).ideal(g)
    else:
        from types import SimpleNamespace
        I=SimpleNamespace(gens=lambda:(g,g)if mode==2 else (g,))
    v=f.inverse_mod(I);return _parts(v)+'|'+str(v)
FUNCTIONS['ff_inverse']=ff_inverse


def ff_order(p,infinite,op,mode,num,den,codes):
    K=FunctionField(QQ if p==0 else GF(p),'x');O=K.maximal_order_infinite()if infinite else K.maximal_order()
    if op==0:
        M=O.ideal_monoid()
        from sage.rings.function_field.ideal import IdealMonoid
        return '|'.join(map(str,[O,O.function_field() is K,O.fraction_field() is K,O.is_field(),O.is_noetherian(),','.join(map(str,(1,)if infinite else O.basis())),M,M.ring() is O,M is O.ideal_monoid(),M is IdealMonoid(O)]))
    if op==7:
        # Bundled order_rational.py:516 returns the literal tuple (1,).
        # Sage 10.3 incorrectly returns the function-field element 1/x.
        basis=(1,)if infinite else O.basis()
        return ';'.join(type(x).__name__+':'+str(x)for x in basis)
    if op==1:
        f=_fraction_from_strings(O,''.join(chr(int(c))for c in codes));return _parts(f)+'|'+str(f)
    if op==5:return str(O.is_subring(K if mode==0 else O if mode==1 else K.maximal_order()if infinite else K.maximal_order_infinite()))
    f=_element(K,num)/_element(K,den)
    if op==2:return _parts(O(f))+'|'+str(O(f))
    if op==6:return str(f in O)
    P=O if mode==0 else K.maximal_order()if infinite else K.maximal_order_infinite();I=P.ideal(f)
    if op==3:return str(O.ideal(I))
    M=O.ideal_monoid()
    if op==4:
        if mode==2:
            from types import SimpleNamespace
            I=SimpleNamespace(gens=lambda:(f,))
        v=M(I);return str(v)+'|'+str(v is I)
    raise ValueError('unknown ff_order operation')
FUNCTIONS['ff_order']=ff_order


def ff_flint_word(op,a,p,e,sign):
    from pathlib import Path
    from sage.misc.cython import cython_import
    global _flint_word_native
    if '_flint_word_native' not in globals():
        _flint_word_native=cython_import(str(Path(__file__).parent.parent/'flint_word_sqrt.pyx'))
    return _flint_word_native.native_word(int(op),a,p,e,sign)
FUNCTIONS['ff_flint_word']=ff_flint_word


def ff_flint_poly(p,op,n,cutoff,a,b):
    from pathlib import Path
    from sage.misc.cython import cython_import
    global _flint_poly_sqrt_native
    if '_flint_poly_sqrt_native' not in globals():
        _flint_poly_sqrt_native=cython_import(str(Path(__file__).parent.parent/'flint_poly_sqrt.pyx'))
    operation=0 if op==0 else 1 if op in (1,2,3)else 2 if op==4 else 3 if op==5 else 4 if op==6 else 5
    # gr_poly/inv_series_basecase.c returns success immediately for len == 0;
    # the older public nmod wrapper instead truncates its input and aborts.
    if op==5 and n==0:return ''
    import os,ctypes
    libc=ctypes.CDLL(None);libc.fflush(None);saved=os.dup(1)
    try:
        os.dup2(2,1)
        value=_flint_poly_sqrt_native.native_poly(operation,a,b,p,n)
    finally:
        libc.fflush(None);os.dup2(saved,1);os.close(saved)
    return 'None'if value is None else ','.join(map(str,value))
FUNCTIONS['ff_flint_poly']=ff_flint_poly


def ff_square(p,op,raw,all_roots,extend,num,den):
    K=FunctionField(QQ if p==0 else GF(p),'x');F=K.field();R=K._ring
    x=F._element_class(F,R(num),R(den),reduce=False)if raw else F(R(num),R(den))
    f=K.element_class(K,x)
    if op==0:return str(f.is_square())
    if op==2:return str(x.is_square())
    if op==4:
        v=x._sqrt_or_None()
        return 'None'if v is None else _parts(v)+'|'+str(v)+'|'+str(v is x)
    if op==5:
        ok,v=x.is_square(root=True)
        return str(ok)+'|'+('None'if v is None else _parts(v))
    v=f.sqrt(all=bool(all_roots))if op==1 else x.sqrt(extend=bool(extend),all=bool(all_roots))
    if isinstance(v,list):return ';'.join(_parts(z)+'|'+str(z)for z in v)
    return _parts(v)+'|'+str(v)+'|'+str(v is (f if op==1 else x))
FUNCTIONS['ff_square']=ff_square


def ff_polynomial_square(p,root,a):
    R=PolynomialRing(ZZ if p==-1 else QQ if p==0 else GF(p),'x');f=R(a)
    if not root:return str(f.is_square())
    ok,g=f.is_square(root=True)
    return str(ok)+'|'+('None'if g is None else _ints(g.list())+'|'+str(g)+'|'+str(g is f)+'|'+str(g.parent() is R))
FUNCTIONS['ff_polynomial_square']=ff_polynomial_square


def ff_ntl_poly(op,a,b):
    from sage.libs.ntl.ntl_ZZX import ntl_ZZX
    A=ntl_ZZX(a)
    if op==0:return _ints(A.gcd(ntl_ZZX(b)).list())
    return ';'.join(_ints(f.list())+':'+str(e)for f,e in A.squarefree_decomposition())
FUNCTIONS['ff_ntl_poly']=ff_ntl_poly


def ff_squarefree(p,a):
    R=PolynomialRing(ZZ if p==-1 else QQ if p==0 else Zmod(p),'x');f=R(a)
    # Bundled polynomial_zmod_flint.pyx checks zero before the field guard.
    # Sage 10.3 instead returns a zero factorization (or the composite error).
    if 2 < p < 2**63 and not f:
        raise ArithmeticError('square-free decomposition of 0 is not defined')
    fac=f.squarefree_decomposition()
    return str(fac.unit())+'|'+ ';'.join(_ints(g.list())+':'+str(e)+':'+str(g.parent() is R)for g,e in fac)
FUNCTIONS['ff_squarefree']=ff_squarefree


def ff_flint_squarefree(p,a):
    global _native_flint_squarefree
    if '_native_flint_squarefree'not in globals():
        from sage.misc.cython import cython_import
        from pathlib import Path
        _native_flint_squarefree=cython_import(str(Path(__file__).resolve().parents[1]/'flint_poly_sqrt.pyx')).native_squarefree
    return ';'.join(_ints(f)+':'+str(e)for f,e in _native_flint_squarefree(a,p))
FUNCTIONS['ff_flint_squarefree']=ff_flint_squarefree

def ff_squarefree_field(p,n,op,e,a):
    if n==1:k=GF(p)
    else:
        mod={2:[1,1,0,1],3:[1,0,1],5:[2,0,1],7:[1,0,1]}[int(p)]
        k=GF(p**n,'a',modulus=PolynomialRing(GF(p),'z')(mod))
    R=PolynomialRing(k,'x')
    def coeff(v):return sum(k((v//p**i)%p)*k.gen()**i for i in range(n))
    f=R([coeff(v)for v in a])
    if op:f=f**e
    if op==2:f*= (R.gen()+k.gen())**(e+1)
    if n==1 and 2<p<2**63 and not f:raise ArithmeticError('square-free decomposition of 0 is not defined')
    fac=f.squarefree_decomposition()
    enc=lambda c:_ints(c.polynomial().list())
    return enc(fac.unit())+'|'+ ';'.join('/'.join(enc(c)for c in g.list())+':'+str(m)+':'+str(g.parent() is R)for g,m in fac)
FUNCTIONS['ff_squarefree_field']=ff_squarefree_field


def ff_gf2_roots(op,a,b,raw,all_roots,extend):
    k=GF(2);x=type(k(ZZ(0)))(k,a)if raw else k(a)
    if op==13:return str(x.is_one())+'|'+repr(x)
    if op==0:return str(x.is_square())
    if op==2:
        R=PolynomialRing(k,'x');f=R([a,b,a]);ok,r=f.is_square(root=True)
        return str(ok)+'|'+('None'if r is None else _ints(r.list())+'|'+str(r.parent()is R))
    if op==3:return str(k.zero()is k(ZZ(0)))+'|'+str(k.one()is k(ZZ(1)))+'|'+str(k.gen()is k(ZZ(1)))
    if op==12:return x.sqrt(oops=True)
    if op==1:v=x.sqrt(all=bool(all_roots),extend=bool(extend))
    elif op==4:v=x+k(b)
    elif op==5:v=x-k(b)
    elif op==6:v=x*k(b)
    elif op==7:v=x/k(b)
    elif op==8:v=-x
    elif op==9:v=~x
    elif op==10:v=x**b
    elif op==11:v=k(x)
    def enc(v):return str(v)+'|'+str(v is x)+'|'+str(v is k(ZZ(v)))
    return ';'.join(enc(z)for z in v)if isinstance(v,list)else enc(v)
FUNCTIONS['ff_gf2_roots']=ff_gf2_roots


def ff_pari_squarefree(p,a):
    global _pari_squarefree_installed
    if '_pari_squarefree_installed'not in globals():
        pari('install("FpX_factor_squarefree","GG","audit_squarefree")')
        _pari_squarefree_installed=True
    try:
        v=pari('audit_squarefree')(pari(a).Polrev(),pari(p))
        return 'value|'+ ';'.join(','.join(str(g.polcoef(i))for i in range(int(g.poldegree())+1))for g in v)
    except Exception as e:
        return 'error|'+type(e).__name__+'|'+str(e).removeprefix('audit_squarefree: ')
FUNCTIONS['ff_pari_squarefree']=ff_pari_squarefree


def ff_place_maps(p,infinite,a,op,mode,num,den):
    from sage.rings.function_field.place import PlaceSet
    from sage.rings.function_field.valuation_ring import FunctionFieldValuationRing
    K=FunctionField(QQ if p==0 else GF(p),'x');x=K.gen();O=K.maximal_order_infinite()if infinite else K.maximal_order()
    I=O.ideal(1/x if infinite else x-a);P=I.place();S=K.place_set();V=P.valuation_ring()
    if op==0:return '|'.join(map(str,[S,S.function_field()is K,S is PlaceSet(K),P,P.degree(),P.is_infinite_place(),_parts(P.local_uniformizer()),P.parent()is S,P.function_field()is K,P.prime_ideal()is I,V,V.place()==P]))
    if op==1:
        Q=S(I);W=FunctionFieldValuationRing(K,Q)
        return '|'.join(map(str,[P is Q,P==Q,V is P.valuation_ring(),V is Q.valuation_ring(),V is W,W.place()is V.place()]))
    if op==2:
        A=P.residue_field();B=V.residue_field();C=P.residue_field(name='a');D=V.residue_field(name='a')
        return '|'.join(map(str,[A[0],A is B,A[1]is B[1],A[2]is B[2],C is D,A is C,A[0]is K.constant_base_field()]))
    if op==6:return str(S(P)is P)
    if op==7:
        Q=K.maximal_order().ideal(x-(a+1)).place();return str((P>Q)-(P<Q))
    f=_element(K,num)/_element(K,den)
    if op==5:return _parts(V(f))+'|'+str(V(f)is f)
    A=P._residue_field()if mode==0 else P.residue_field()if mode==1 else V.residue_field(name='b')
    k,fr,to=A
    if op==3:
        c=K.constant_base_field()(a);g=fr(c);return _parts(g)+'|'+str(g.parent()is K)
    if op==4:return str(to(f))
    if op==8:return str(to(num[0]if num else ZZ(0)))
    raise ValueError('unknown ff_place_maps operation')
FUNCTIONS['ff_place_maps']=ff_place_maps


def ff_poly_numden(p,a,b):
    k=ZZ if p==-1 else QQ if p==0 else GF(p);R=PolynomialRing(k,'x')
    f=R([k(n)/k(d)if p!=-1 else k(n)for n,d in zip(a,b)])
    n=f.numerator();d=f.denominator()
    return '|'.join(map(str,[','.join(map(str,n.list())),d,n.parent().base_ring(),n.parent()is f.parent(),n.parent()is f.numerator().parent()]))
FUNCTIONS['ff_poly_numden']=ff_poly_numden


def ff_poly_nested_numden(mode,a,b,e):
    T=PolynomialRing(QQ if mode!=2 else GF(5),'t');t=T.gen();k=T if mode==0 else T.fraction_field();R=PolynomialRing(k,'x')
    cs=[k(T([n/ZZ(d),1]))if mode==0 else k(n)/(k(t+d)**j)for n,d,j in zip(a,b,e)]
    f=R(cs);n=f.numerator();d=f.denominator()
    return '|'.join(map(str,[n,d,n.parent()is f.parent()]))
FUNCTIONS['ff_poly_nested_numden']=ff_poly_nested_numden


def ff_poly_lcm(p,a,b):
    if p in (-3,-2):
        from pathlib import Path
        from sage.misc.cython import cython_import
        global _flint_poly_sqrt_native
        if '_flint_poly_sqrt_native'not in globals():_flint_poly_sqrt_native=cython_import(str(Path(__file__).parent.parent/'flint_poly_sqrt.pyx'))
        r=_flint_poly_sqrt_native.native_lcm(0 if p==-3 else 1,a,b)
        return ','.join(map(str,r))if p==-3 else ','.join(map(str,r[0]))+'/'+str(r[1])
    k=ZZ if p==-1 else QQ if p==0 else GF(p);R=PolynomialRing(k,'x');f=R(a);g=R(b)
    return ','.join(map(str,f.lcm(g).list()))
FUNCTIONS['ff_poly_lcm']=ff_poly_lcm


def ff_divisor(p,op,keys,ms,ns,scalar,num,den,vec):
    from sage.rings.function_field.divisor import divisor,prime_divisor,DivisorGroup
    K=FunctionField(QQ if p==0 else GF(p),'x');x=K.gen();G=K.divisor_group()
    def place(a):return K.maximal_order_infinite().ideal(1/x).place()if a==-99 else K.maximal_order().ideal(x-a).place()
    ps=[place(a)for a in keys];D=divisor(K,dict(zip(ps,ms)));E=divisor(K,dict(zip(ps,ns)));P=ps[0]if ps else place(0)
    def ds(d):return str(d)+'|'+ ';'.join(str(q)+':'+str(m)for q,m in d.list())
    if op==0:return '|'.join(map(str,[ds(D),D.degree(),D.is_effective(),','.join(map(str,D.support())),D.parent()is G,D.parent()._field is K,D.dict()is D.dict(),G,G is DivisorGroup(K),G.function_field()is K]))
    if op==1:return ds(-D)
    if op==2:return ds(D+E)
    if op==3:return ds(D-E)
    if op==4:return ds(scalar*D)
    if op==5:return str((D>E)-(D<E))
    if op==6:return ds(D.numerator())+'|'+ds(D.denominator())
    if op==7:return str(D.multiplicity(place(scalar)))+'|'+str(D.valuation(place(scalar)))
    if op==8:return ds(G(D))+'|'+str(G(D)is D)
    if op==9:return ds(G(P))
    if op==10:
        v=ZZ(scalar)if not vec else K(scalar);return ds(G(v))
    if op==11:return str(D.dimension())+'|'+ ';'.join(_parts(f)for f in D._basis())+'|'+ ';'.join(_parts(f)for f in D.basis_function_space())
    if op==12:
        f=_element(K,num)/_element(K,den);return ','.join(map(str,D.function_space()[2](f)))
    if op==13:
        v=D.function_space()[1]([K.constant_base_field()(n)for n in vec]);return ('ff|'+_parts(v))if v in K and hasattr(v,'parent')and v.parent()is K else 'scalar|'+str(v)
    if op==14:return D._format(lambda x:'<'+str(x)+'>',' @ ',' / ')
    if op==15:
        A=D._function_space();B=D._function_space();C=D.function_space();E=D.function_space();return '|'.join(map(str,[A is B,A[0]is D.basis_function_space(),A[1]is B[1],C is E,C[1]is E[1],C[2]is E[2]]))
    if op==16:return ';'.join(_parts(f)for f in D._differential_space()[0])
    if op==17:return ds(prime_divisor(K,P,scalar))
    if op==18:
        D.dict()[P]=scalar;return ds(D)
    if op==19:return '|'.join(map(str,[G.zero()is G.zero(),G(ZZ(0))is G(ZZ(0)),G(G.zero())is G.zero()]))
    raise ValueError('unknown ff_divisor operation')
FUNCTIONS['ff_divisor']=ff_divisor


def ff_divisor_echelon(p,modes,num,den):
    K=FunctionField(QQ if p==0 else GF(p),'x');x=K.gen();D=K.divisor_group().zero()
    choices=[K.zero(),K.one(),x,x+1,2*x,1/x,1/(x+1),x/(x+1)]
    bs,co=D._echelon_basis([choices[int(m)]for m in modes]);f=_element(K,num)/_element(K,den)
    return ';'.join(_parts(b)for b in bs)+'|'+','.join(map(str,co(f)))
FUNCTIONS['ff_divisor_echelon']=ff_divisor_echelon


def ff_ideal(p,infinite,raw,op,num,den,num2,den2,e):
    from sage.rings.function_field.ideal_rational import FunctionFieldIdeal_rational,FunctionFieldIdealInfinite_rational
    K=FunctionField(QQ if p==0 else GF(p),'x');O=K.maximal_order_infinite()if infinite else K.maximal_order()
    f=_element(K,num)/_element(K,den);g=_element(K,num2)/_element(K,den2)
    cls=FunctionFieldIdealInfinite_rational if infinite else FunctionFieldIdeal_rational
    I=cls(O,f)if raw else O.ideal(f);J=cls(O,g)if raw else O.ideal(g)
    def fmt(I):return str(I)+'|'+_parts(I.gen())
    if op==0:return '|'.join(map(str,[fmt(I),I._repr_short(),I.is_zero(),I.is_prime(),I.ring()is O,I.base_ring()is O,I.parent()is O.ideal_monoid(),I.gen()is I.gens()[0],I.gen()is I.gens_over_base()[0],I.gen()is f]))
    if op==1:return str((I>J)-(I<J))
    if op==2:return fmt(I+J)
    if op==3:return fmt(I*J)
    if op==4:return fmt(I/J)
    if op==5:return fmt(~I)
    if op==6:
        v=I**e;return fmt(v)+'|'+str(v is I)
    if op==7:return str(g in I)
    if op==8:return str(I.valuation(J))
    if op==9:return str(I.place())
    if op==10:return fmt(g*I)
    if op==11:return str(ZZ(e)in I)
    if op==12:return ','.join(map(str,I.denominator().list()))
    if op==13:return fmt(O.ideal([f,g]))
    if op in (14,15,16,17):
        if not infinite and I.gen().is_zero():raise ArithmeticError('factorization of 0 is not defined')
        if op==14:return ';'.join(_parts(P.gen())+':'+str(m)for P,m in I.factor())
        return str(I.divisor()if op==15 else I.divisor_of_zeros()if op==16 else I.divisor_of_poles())
    if op==18:return ';'.join(_parts(g)for g in I.gens_reduced())
    if op==19:return '|'.join(map(str,[I**0 is I**0,I**1 is I]))
    if op==20:return '|'.join(map(str,[I.parent().one(),I.parent().one()is I.parent().one(),I**0 is I.parent().one()]))
    raise ValueError('unknown ff_ideal operation')
FUNCTIONS['ff_ideal']=ff_ideal


def ff_enumeration(p,kind,op,mode,d,limit):
    from itertools import islice
    k=QQ if p==0 else Zmod(p)if kind==1 else GF(p**2,'a',modulus=([1,1,1]if p==2 else [1,0,1]if p==3 else [2,0,1]),impl='pari_ffelt')if kind==2 else GF(p)
    if op==10:return str(k.cardinality())
    if op==11:
        # The port's finite-cardinality adapter explicitly rejects infinite fields.
        if not k.is_finite():raise TypeError(f'{k} is not a finite field')
        return str(k.cardinality())
    R=PolynomialRing(k,'x');d=int(d);limit=int(limit)
    opts={'of_degree':d}if mode==0 else {'max_degree':d}if mode==1 else {}if mode==2 else {'of_degree':d,'max_degree':d}
    if op==6:
        # Negative QQ/word FLINT dictionary exponents can segfault or request
        # exabytes in the original unchecked C kernels. Compare the documented
        # safe guard instead of invoking undefined native memory behavior.
        if d<0 and (p==0 or kind!=2 and 2<p<2**64):
            class RangeError(Exception):pass
            raise RangeError('negative monomial degree is outside the native FLINT contract')
        return ','.join(map(str,R.monomial(d).list()))
    if op>=7:
        K=FunctionField(k,'x')
        if op==8:return str(K.get_place(d))
        return ';'.join(map(str,islice(K._places_finite(d),limit)))if op==7 else ';'.join(map(str,K.places(d)))
    g=R.polynomials(**opts)if op==0 else R.monics(**opts)if op==1 else R._polys_degree(d)if op==2 else R._polys_max(d)if op==3 else R._monics_degree(d)if op==4 else R._monics_max(d)
    return ';'.join(','.join(map(str,f.list()))for f in islice(g,limit))
FUNCTIONS['ff_enumeration']=ff_enumeration


def ff_mrange(sizes,kind,mode,limit):
    from sage.misc.mrange import _xmrange_iter
    from itertools import islice
    sizes=list(map(int,sizes));counts=[0]*len(sizes)
    class Input:
        def __init__(self,n,i):self.n=n;self.i=i
        def is_finite(self):
            if kind in (4,5,7,10):raise AttributeError()
            if kind==3:raise ValueError('finiteness unknown')
            if kind==6:raise RuntimeError('finiteness unknown')
            return False if kind in (8,9) else self.n>=0
        def cardinality(self):
            if kind==4:return infinity if self.n<0 else self.n
            if kind==5:raise NotImplementedError()
            if kind==8:raise TypeError('cardinality unknown')
            if kind in (9,10):raise ValueError('cardinality unknown')
            raise AttributeError()
        def __len__(self):
            if self.n<0:raise TypeError()
            return self.n
        def __iter__(self):
            v=0
            while self.n<0 or v<self.n:
                counts[self.i]+=1;yield v;v+=1
    ins=[list(range(n))if kind==0 else Input(n,i)for i,n in enumerate(sizes)]
    if kind==2:ins=list(map(iter,ins))
    def convert(*args):
        if not args:return 777 if mode==1 else [999]
        return sum(args[0])if mode==1 else args[0]
    rows=list(islice(_xmrange_iter(ins,list if mode==0 else convert),int(limit)))
    same=len(rows)>1 and rows[0]is rows[-1]if mode!=1 else False
    return ';'.join(','.join(map(str,r))if mode!=1 else str(r)for r in rows)+'|'+str(same)+'|'+','.join(map(str,counts))
FUNCTIONS['ff_mrange']=ff_mrange


def ff_mpoly_repr(p,order,cs):
    R=PolynomialRing(QQ if p==0 else GF(p),['x0','x1','x2'],order=['lex','deglex','degrevlex'][int(order)])
    x,y,z=R.gens();ms=[R.one(),x**2,y*z,x*y,z**2]
    return str(sum((R(c)*m for c,m in zip(cs,ms)),R.zero()))
FUNCTIONS['ff_mpoly_repr']=ff_mpoly_repr

FUNCTIONS['qq_cardinality']=lambda:str(QQ.cardinality())


def ff_element_observe(p,raw,op,num,den,a):
    K=FunctionField(QQ if p==0 else GF(p),'x');F=K.field();R=K._ring
    f=K.element_class(K,F._element_class(F,R(num),R(den),reduce=False))if raw else _element(K,num)/_element(K,den)
    if op==0:return '|'.join(map(str,[_parts(f),_parts(f.matrix()[0,0]),_parts(f.trace()),_parts(f.norm()),f.matrix()[0,0]is f,f.trace()is f,f.norm()is f,f.matrix()is f.matrix()]))
    if op==1:
        D=f.divisor();return str(D)+'|'+str(D is f.divisor())+'|'+str(D.degree())
    if op==2:return str(f.divisor_of_zeros())
    if op==3:return str(f.divisor_of_poles())
    if op==4:return ';'.join(map(str,f.zeros()))
    if op==5:return ';'.join(map(str,f.poles()))
    if op in (6,7):
        P=K.maximal_order().ideal(K.gen()-a).place()if op==6 else K.maximal_order_infinite().prime_ideal().place()
        r=f.evaluate(P);return str(r)+'|'+str(r.parent()is K.constant_base_field())
    if op==8:
        M=f.matrix();return str(M.is_immutable())+'|'+str(M.nrows())+'|'+str(M.ncols())
    if op==9:
        M=f.matrix();M[0,0]=K(a);return _parts(M[0,0])
    if op==10:
        from copy import copy
        M=f.matrix();C=copy(M);C[0,0]=K(a)
        return '|'.join(map(str,[_parts(M[0,0]),_parts(C[0,0]),_parts(f),C.is_immutable()]))
    if op==11:
        before=_parts(f);M=f.matrix();t=f.trace();n=f.norm()
        return '|'.join([before,_parts(f),_parts(M[0,0]),_parts(t),_parts(n)])
    raise ValueError('unknown ff_element_observe operation')
FUNCTIONS['ff_element_observe']=ff_element_observe


def ff_predicate(p,op,kind):
    from sage.rings.function_field.element import FunctionFieldElement
    from sage.rings.function_field.function_field import FunctionField as FunctionFieldBase
    from sage.categories.function_fields import FunctionFields
    from sage.categories.fields import Fields
    K=FunctionField(QQ if p==0 else GF(p),'x');counts=[0,0]
    class CategoryProbe:pass
    class Probe:
        def parent(self):
            counts[0]+=1
            if kind==21:raise ValueError('parent unavailable')
            return K if kind==17 or kind==19 and counts[0]>1 else QQ if kind in (18,19) else self
        def category(self):
            counts[1]+=1
            if kind==25:raise ValueError('category unavailable')
            if kind==26:raise AttributeError('category unavailable')
            if kind in (27,40):return CategoryProbe()
            if kind==41:return 0
            if kind==28:return None
            return FunctionFields()if kind==23 else Fields()
    vals=[K.gen(),ZZ(1),QQ(1)/2,K._ring.gen(),K.field().gen(),K,K.maximal_order(),K.maximal_order().ideal(K.gen()),K.maximal_order().ideal(K.gen()).place(),QQ,GF(3),[ZZ(1)],{},None,True,float(1.5),'x']
    if kind<len(vals):x=vals[int(kind)]
    else:
        x=Probe()
        if kind==20:del Probe.parent
        if kind==29:x.category=1
        if kind==30:x.parent=1
        if kind==31:x=lambda:ZZ(1)
        if kind==32:x.category=ZZ(1)
        if kind==33:x.parent=ZZ(1)
        if kind==34:x.parent={}
        if kind==35:x.parent=[]
        if kind==36:x.parent=None
        if kind==37:x.parent=False
        if kind==38:x.parent='x'
        if kind==39:x.category=None
        if kind==40:CategoryProbe.is_subcategory=ZZ(1)
        if kind==42:x=ZZ(1)
        if kind==43:x=GF(3 if p==0 else p)(1)
        if kind==44:x=GF(2)(1)
        if kind==45:x=matrix(QQ,[[1]])
    # Execute the bundled predicate body. Sage 10.3 predates its category fallback;
    # omit only the warning hook, which is not part of the result transcript.
    if op==0:
        if isinstance(x,FunctionFieldElement):r=True
        elif isinstance(x.parent(),FunctionFieldBase):r=True
        else:r=x.parent()in FunctionFields()
    else:r=isinstance(x,FunctionFieldBase)or x in FunctionFields()
    return str(r)+'|'+','.join(map(str,counts))
FUNCTIONS['ff_predicate']=ff_predicate


def ff_predicate_access(p,op,kind):
    from sage.rings.function_field.element import FunctionFieldElement
    from sage.rings.function_field.function_field import FunctionField as FunctionFieldBase
    from sage.categories.function_fields import FunctionFields
    K=FunctionField(QQ if p==0 else GF(p),'x');counts=[0]*6
    def parent_get():
        counts[0]+=1;n=counts[0]
        if kind==2:raise AttributeError('parent unavailable')
        def call():
            counts[1]+=1
            return (QQ if n==1 else K)if kind==7 else x
        return call
    def category_get():
        counts[2]+=1
        if kind==3:raise AttributeError('category unavailable')
        if kind==4:raise ValueError('category unavailable')
        if kind==6:return None
        def call():counts[3]+=1;return c
        return call
    def sub_get():
        counts[4]+=1
        if kind==5:raise AttributeError('subcategory unavailable')
        def call(C):counts[5]+=1;return C is FunctionFields()
        return call
    class Probe:
        parent=property(lambda self:parent_get())
        category=property(lambda self:category_get())
    class CategoryProbe:
        is_subcategory=property(lambda self:sub_get())
    class DynamicProbe:
        def __getattr__(self,name):
            if name=='parent':return parent_get()
            if name=='category':return category_get()
            raise AttributeError(name)
    class DynamicCategory:
        def __getattr__(self,name):
            if name=='is_subcategory':return sub_get()
            raise AttributeError(name)
    x=DynamicProbe()if kind==1 else Probe();c=DynamicCategory()if kind==1 else CategoryProbe()
    if op==0:
        if isinstance(x,FunctionFieldElement):r=True
        elif isinstance(x.parent(),FunctionFieldBase):r=True
        else:r=x.parent()in FunctionFields()
    else:r=isinstance(x,FunctionFieldBase)or x in FunctionFields()
    return str(r)+'|'+','.join(map(str,counts))
FUNCTIONS['ff_predicate_access']=ff_predicate_access


def ff_matrix_base(p,raw,warm,kind,num,den):
    K=FunctionField(QQ if p==0 else GF(p),'x');F=K.field();R=K._ring
    # Factory parents survive between cases. Start each scenario with a cold
    # native free-module cache; subsequent calls deliberately share it.
    K.free_module.clear_cache()
    f=K.element_class(K,F._element_class(F,R(num),R(den),reduce=False))if raw else _element(K,num)/_element(K,den)
    bases=[None,K,QQ,FunctionField(K.constant_base_field(),'y'),F,R,False,ZZ(0),[],{},'x',K.gen(),K.maximal_order(),ZZ,GF(3)]
    class MatrixBaseProbe:pass
    bases += [set(),{},MatrixBaseProbe(),{},lambda:None]
    base=bases[int(kind)]
    def observe(g,b):
        try:return _parts(g.matrix(b)[0,0])
        except Exception as e:return type(e).__name__+':'+str(e)
    if warm==1:f.matrix()
    if warm==2:f.trace()
    if warm==3:f.norm()
    if warm==4:K.one().matrix()
    if warm==5:observe(f,QQ)
    if warm==6:FunctionField(K.constant_base_field(),'z').gen().matrix()
    first=observe(f,base);second=observe(f,base);success=observe(f,K);last=observe(K.one(),base)
    return '|'.join([first,second,success,last,_parts(f)])
FUNCTIONS['ff_matrix_base']=ff_matrix_base


def ff_variable_map(p,same,op,raw,kind,num,den):
    K=FunctionField(QQ if p==0 else GF(p),'x');L,from_L,to_L=K.change_variable_name('x'if same else'y')
    source=L if op==0 else K;target=K if op==0 else L;map=from_L if op==0 else to_L
    R=source._ring;F=source.field()
    f=source.element_class(source,F._element_class(F,R(num),R(den),reduce=False))if raw else _element(source,num)/_element(source,den)
    vals=[f,ZZ(1),QQ(1)/2,source.constant_base_field()(1),R(num),target._ring(num),f.element(),num,None,'1','badname',{},True,float(1.5),target.gen()]
    vals += [{2:ZZ(3),0:ZZ(1)},{'bad':ZZ(1)},[],[[ZZ(1),ZZ(2)]],ZZ(1)]
    vals += [{key:ZZ(1)}for key in ["a'b",'a\"b',"a'b\"c",'a\nb','a\\b']]
    x=vals[int(kind)];result=map(x)
    return '|'.join([_parts(result),str(result.parent()is target),str(result is x),_parts(f)])
FUNCTIONS['ff_variable_map']=ff_variable_map


def ff_poly_dict(p,keykind,valuekind,degrees,values):
    k=ZZ if p==-1 else QQ if p==0 else GF(p);R=PolynomialRing(k,'x')
    def value(v):
        if valuekind==1:return [v,ZZ(1)]
        if valuekind==2:return R([v,1])
        if valuekind==3:return QQ(v)/2
        if valuekind==4:return str(v)
        if valuekind==5:return 'x'
        if valuekind==6:return '1/2'
        if valuekind==7:return None
        if valuekind==8:return True
        if valuekind==9:return float(v)+0.5
        if valuekind==10:return []
        if valuekind==11:return [[]]
        if valuekind==12:return R(v)
        if valuekind==13:return R.zero()
        return v
    d={(tuple([e])if keykind==2 else str(e)if keykind==3 else float(e)+0.5 if keykind==4 else ()if keykind==5 else (e,7)if keykind==6 else int(e)if keykind==7 else e):value(v)for e,v in zip(degrees,values)}
    # QQ/word-FLINT negative exponents invoke undefined native kernels. Apply
    # the existing documented safety guard instead of running those kernels.
    if keykind in (0,1) and any(e<0 for e in degrees) and (p==0 or 2<p<2**64):
        class RangeError(Exception):pass
        raise RangeError('negative monomial degree is outside the native FLINT contract')
    # JS Records enumerate nonnegative integer keys in ascending order.
    if keykind==1:d=dict(sorted(d.items(),key=lambda kv:(0,int(kv[0]))if 0<=kv[0]<2**32-1 else(1,list(d).index(kv[0]))))
    f=R(d);return ','.join(map(str,f.list()))
FUNCTIONS['ff_poly_dict']=ff_poly_dict


def ff_poly_dict_order(p,kind):
    k=ZZ if p==-1 else QQ if p==0 else GF(p);R=PolynomialRing(k,'x')
    items=[(ZZ(0),[ZZ(1),ZZ(2)]),('bad',ZZ(1))]
    if kind==1:items.reverse()
    if kind==2:items=[(ZZ(0),'bad'),(ZZ(2),ZZ(1))]
    if kind==3:items=[(ZZ(0),[ZZ(1),ZZ(2)]),((ZZ(2),),ZZ(1))]
    return ','.join(map(str,R(dict(items)).list()))
FUNCTIONS['ff_poly_dict_order']=ff_poly_dict_order


def ff_flint_rng(op,seed1,seed2,limit,count):
    from flint_random import native_rng
    return native_rng(op,seed1,seed2,limit,count)
FUNCTIONS['ff_flint_rng']=ff_flint_rng


def ff_flint_factor(op,p,d,a,b,c,seed1,seed2):
    from flint_factor import native_factor
    return native_factor(op,p,d,a,b,c,seed1,seed2)
FUNCTIONS['ff_flint_factor']=ff_flint_factor


def ff_word_factor(p,a):
    R=PolynomialRing(GF(p),'x');f=R(a)
    # Bundled Sage added the zero check to polynomial_zmod_flint.factor;
    # installed Sage 10.3 predates it. Other backends execute their own checks.
    if not f and 2<p<2**63: raise ArithmeticError('factorization of 0 is not defined')
    factors=f.factor()
    out=list(factors)
    if f.degree()==0 or factors.unit()!=1:out.insert(0,(R(factors.unit()),1))
    return ';'.join(','.join(map(str,g.list()))+':'+str(e)for g,e in out)
FUNCTIONS['ff_word_factor']=ff_word_factor

def ff_flint_matrix(p,rows,inner,cols):
    a=[(i*17+(i//inner)*13+3)%p for i in range(int(rows*inner))]
    b=[(i*23+(i//cols)*7+5)%p for i in range(int(inner*cols))]
    return ff_flint_factor(18,p,inner,a,b,[],rows,cols)
FUNCTIONS['ff_flint_matrix']=ff_flint_matrix


def ff_pari_rng(op,seed,limit,count,length):
    from pari_random import native_pari_rng
    return native_pari_rng(op,seed,limit,count,length)
FUNCTIONS['ff_pari_rng']=ff_pari_rng


def ff_pari_binary(op,a,b,c,d,e):
    from pari_binary import native_pari_binary
    return native_pari_binary(op,a,b,c,d,e)
FUNCTIONS['ff_pari_binary']=ff_pari_binary


def ff_binary_factor(bits,kind):
    k=Zmod(2)if kind==2 else GF(2);R=PolynomialRing(k,'x')
    f=R([(bits>>i)&1 for i in range(int(bits).bit_length())]);factors=f.factor();out=list(factors)
    if f.degree()==0 or factors.unit()!=1:out.insert(0,(R(factors.unit()),1))
    return ';'.join(','.join(map(str,g.list()))+':'+str(e)for g,e in out)
FUNCTIONS['ff_binary_factor']=ff_binary_factor


def ff_pari_quotient(op,p,n,x,T):
    from pari_quotient import native_pari_quotient
    return native_pari_quotient(op,p,n,x,T)
FUNCTIONS['ff_pari_quotient']=ff_pari_quotient


def ff_pari_matrix(op,p,m,n,k,ab,bb,seed):
    from pari_matrix import native_pari_matrix
    return native_pari_matrix(op,p,m,n,k,ab,bb,seed)
FUNCTIONS['ff_pari_matrix']=ff_pari_matrix


def ff_pari_matrix_guard(op,p,rows,x,y):
    # Preserve unsafe word/allocation guards; generic empty matrices return before p.
    class RangeError(Exception):pass
    if op==0 and not 1<=p<2**64:raise RangeError('modulus must be a positive word integer')
    if op==2:
        if not 0<=rows<=2**53-1:raise RangeError('row count must be nonnegative')
        if any(c<0 or c>>int(rows)for c in x)or any(c<0 or c>>len(x)for c in y):raise RangeError('column exceeds matrix row count')
    else:
        from pari_matrix import native_pari_matrix
        native_pari_matrix(1 if op==0 else 2,p,0,0,0,0,0,0)
    return True
FUNCTIONS['ff_pari_matrix_guard']=ff_pari_matrix_guard


def ff_pari_composition(op,p,n,Q,x,T):
    from pari_composition import native_pari_composition
    return native_pari_composition(op,p,n,Q,x,T)
FUNCTIONS['ff_pari_composition']=ff_pari_composition


def ff_pari_polymul(op,p,a,b):
    from pari_polymul import native_pari_polymul
    return native_pari_polymul(op,p,a,b)
FUNCTIONS['ff_pari_polymul']=ff_pari_polymul


def ff_pari_polymul_pattern(op,p,n,m,ea,eb,kind):
    def poly(length,e,salt):
        out=[]
        for i in range(int(length)):
            c=Integer(0)if kind==0 else (Integer(1)<<int(e))+(Integer(i)*17+salt*31)%(Integer(1)<<int(e))
            if kind==2 and i%2:c=-c
            if kind==3 and i%5==0:c=Integer(0)
            if op in (4,5):c=p-1 if kind==4 else c%p
            out.append(c)
        return out
    return ff_pari_polymul(op,p,poly(n,ea,1),poly(m,eb,2))
FUNCTIONS['ff_pari_polymul_pattern']=ff_pari_polymul_pattern


def ff_pari_polydiv(op,p,a,b):
    from pari_polydiv import native_pari_polydiv
    return native_pari_polydiv(op,p,a,b)
FUNCTIONS['ff_pari_polydiv']=ff_pari_polydiv


def ff_pari_polydiv_pattern(op,p,n,m,kind):
    def poly(length,divisor):
        out=[]
        for i in range(int(length)):
            if divisor and i==length-1:c=p-1 if kind==1 and p>2 else 1
            else:
                c=(i*i*17+i*(13 if divisor else 31)+7)%p
                if kind==0:c=0
                if kind==2 and op<3:c=c-p if i%2 else c+p
                if kind==3 and i%5!=0:c=0
            out.append(c)
        return out
    return ff_pari_polydiv(op,p,poly(n,op in (2,5)),poly(m,True))
FUNCTIONS['ff_pari_polydiv_pattern']=ff_pari_polydiv_pattern


def ff_pari_polygcd(op,p,a,b):
    from pari_polygcd import native_pari_polygcd
    return native_pari_polygcd(op,p,a,b)
FUNCTIONS['ff_pari_polygcd']=ff_pari_polygcd


def ff_pari_polygcd_pattern(op,p,n,m,kind):
    def poly(length,salt):
        state=salt;mask=(1<<128)-1;out=[]
        for i in range(int(length)):
            state=(state*6364136223846793005+1442695040888963407)&mask
            c=state%p
            if kind==0 or(kind==3 and i%5!=0):c=0
            out.append(c)
        if kind==2 and out:
            v=[0]*(len(out)+2)
            for i,c in enumerate(out):v[i]=(v[i]+c)%p;v[i+2]=(v[i+2]+c)%p
            return v
        return out
    return ff_pari_polygcd(op,p,poly(n,17),poly(m,43))
FUNCTIONS['ff_pari_polygcd_pattern']=ff_pari_polygcd_pattern

def ff_pari_minpoly(op,seed,p,a,b):
    from pari_minpoly import native_pari_minpoly
    return native_pari_minpoly(op,seed,p,a,b)
FUNCTIONS['ff_pari_minpoly']=ff_pari_minpoly

def ff_pari_ddf(op,p,a):
    from pari_ddf import native_pari_ddf
    return native_pari_ddf(op,p,a)
FUNCTIONS['ff_pari_ddf']=ff_pari_ddf

def ff_pari_sqrt(op,p,a,y):
    from pari_sqrt import native_pari_sqrt
    return native_pari_sqrt(op,p,a,y)
FUNCTIONS['ff_pari_sqrt']=ff_pari_sqrt

def ff_pari_normalize(op,p,a):
    from pari_normalize import native_pari_normalize
    return native_pari_normalize(op,p,a)
FUNCTIONS['ff_pari_normalize']=ff_pari_normalize

def ff_pari_factor(op,seed,p,a):
    from pari_factor import native_pari_factor
    return native_pari_factor(op,seed,p,a)
FUNCTIONS['ff_pari_factor']=ff_pari_factor

def ff_pari_roots(op,p,a):
    from pari_roots import native_pari_roots
    return native_pari_roots(int(op),int(p),list(map(int,a)))
FUNCTIONS['ff_pari_roots']=ff_pari_roots

def ff_pari_extension(mode,op,p,T,alen,aflat,blen,bflat):
    from pari_extension import native_pari_extension
    def decode(lengths,flat):
        values=[];offset=0
        for length in lengths:
            n=int(length)
            if n<0:values.append(int(flat[offset]));offset+=1
            else:values.append(list(map(int,flat[offset:offset+n])));offset+=n
        return values
    mode=int(mode)
    return native_pari_extension(mode,int(op),int(p),int(T[0]) if mode==2 else list(map(int,T)),decode(alen,aflat),decode(blen,bflat))
FUNCTIONS['ff_pari_extension']=ff_pari_extension

def ff_pari_extension_large(degree,zero,host):
    import hashlib
    from pari_extension import native_pari_extension
    d=int(degree)
    result=native_pari_extension(1,0,7,[1]+[0]*(d-1)+[1],[[1]*d],[] if zero else [[1]])
    return hashlib.sha256(result.encode()).hexdigest()
FUNCTIONS['ff_pari_extension_large']=ff_pari_extension_large

def ff_pari_extension_division(mode,op,p,T,alen,aflat,blen,bflat):
    from pari_extension_division import native_pari_extension_division
    def decode(lengths,flat):
        values=[];offset=0
        for length in lengths:
            n=int(length)
            if n<0:values.append(int(flat[offset]));offset+=1
            else:values.append(list(map(int,flat[offset:offset+n])));offset+=n
        return values
    mode=int(mode)
    return native_pari_extension_division(mode,int(op),int(p),int(T[0]) if mode==2 else list(map(int,T)),decode(alen,aflat),decode(blen,bflat))
FUNCTIONS['ff_pari_extension_division']=ff_pari_extension_division

def ff_pari_extension_division_pattern(mode,op,degree):
    import hashlib
    from pari_extension_division import native_pari_extension_division
    m=int(mode);n=int(degree);p=2 if m==2 else 17 if m==1 else 2**64+13
    T=7 if m==2 else [3,0,1]
    a=[i%4 for i in range(n+1)] if m==2 else [[i%17,(i//17)%17] for i in range(n+1)]
    a[-1]=1 if m==2 else [1]
    b=[3,1,1] if m==2 else [[1,1],[1],[1]]
    result=native_pari_extension_division(m,int(op),p,T,a,b)
    return hashlib.sha256(result.encode()).hexdigest()
FUNCTIONS['ff_pari_extension_division_pattern']=ff_pari_extension_division_pattern

def ff_pari_extension_gcd(mode,op,p,T,alen,aflat,blen,bflat):
    from pari_extension_gcd import native_pari_extension_gcd
    def decode(lengths,flat):
        values=[];offset=0
        for length in lengths:
            n=int(length)
            if n<0:values.append(int(flat[offset]));offset+=1
            else:values.append(list(map(int,flat[offset:offset+n])));offset+=n
        return values
    mode=int(mode)
    return native_pari_extension_gcd(mode,int(op),int(p),int(T[0]) if mode==2 else list(map(int,T)),decode(alen,aflat),decode(blen,bflat))
FUNCTIONS['ff_pari_extension_gcd']=ff_pari_extension_gcd

def ff_pari_extension_gcd_termination(degree,op):
    from pari_extension_gcd import native_binary_extension_gcd_termination
    return native_binary_extension_gcd_termination(int(degree),int(op))
FUNCTIONS['ff_pari_extension_gcd_termination']=ff_pari_extension_gcd_termination

def ff_pari_extension_quotient(mode,op,p,n,T,alen,aflat,blen,bflat,slen,sflat):
    from pari_extension_quotient import native_pari_extension_quotient
    def decode(lengths,flat):
        values=[];offset=0
        for length in lengths:
            size=int(length)
            if size<0:values.append(int(flat[offset]));offset+=1
            else:values.append(list(map(int,flat[offset:offset+size])));offset+=size
        return values
    mode=int(mode)
    return native_pari_extension_quotient(mode,int(op),int(p),int(n),int(T[0]) if mode==2 else list(map(int,T)),
                                         decode(alen,aflat),decode(blen,bflat),decode(slen,sflat))
FUNCTIONS['ff_pari_extension_quotient']=ff_pari_extension_quotient

def _extension_decode(lengths,flat):
    out=[];offset=0
    for length in lengths:
        n=int(length)
        if n<0:out.append(int(flat[offset]));offset+=1
        else:out.append(list(map(int,flat[offset:offset+n])));offset+=n
    return out

def _extension_columns(shape,lengths,flat):
    coefficients=_extension_decode(lengths,flat);out=[];offset=0
    for length in shape:
        n=int(length);out.append(coefficients[offset:offset+n]);offset+=n
    return out

def ff_pari_extension_matrix(cache,p,T,ashape,alen,aflat,bshape,blen,bflat):
    from pari_extension_matrix import native_pari_extension_matrix
    return native_pari_extension_matrix(int(cache),int(p),list(map(int,T)),_extension_columns(ashape,alen,aflat),_extension_columns(bshape,blen,bflat))
FUNCTIONS['ff_pari_extension_matrix']=ff_pari_extension_matrix

def ff_pari_extension_composition(mode,code,p,T,qlen,qflat,xlen,xflat,vshape,vlen,vflat,slen,sflat):
    from pari_extension_composition import native_pari_extension_composition
    mode=int(mode)
    return native_pari_extension_composition(mode,int(code),int(p),int(T[0])if mode==2 else list(map(int,T)),
      _extension_decode(qlen,qflat),_extension_decode(xlen,xflat),_extension_columns(vshape,vlen,vflat),_extension_decode(slen,sflat))
FUNCTIONS['ff_pari_extension_composition']=ff_pari_extension_composition

def ff_pari_inner_cache(op,p,n,Q,x,vlen,vflat,T):
    from pari_inner_cache import native_pari_inner_cache
    return native_pari_inner_cache(int(op),int(p),int(n),list(map(int,Q)),list(map(int,x)),_extension_decode(vlen,vflat),list(map(int,T)))
FUNCTIONS['ff_pari_inner_cache']=ff_pari_inner_cache

def ff_pari_aut_cache(op,p,n,x,a,T):
    from pari_aut_cache import native_pari_aut_cache
    return native_pari_aut_cache(int(op),int(p),int(n),list(map(int,x)),list(map(int,a)),list(map(int,T)))
FUNCTIONS['ff_pari_aut_cache']=ff_pari_aut_cache

def ff_pari_coefficient_composition(mode,op,p,T,plen,pflat,x,vlen,vflat):
    from pari_coefficient_composition import native_pari_coefficient_composition
    mode=int(mode)
    return native_pari_coefficient_composition(mode,int(op),int(p),int(T[0])if mode>=2 else list(map(int,T)),
      _extension_decode(plen,pflat),int(x[0])if mode>=2 else list(map(int,x)),_extension_decode(vlen,vflat))
FUNCTIONS['ff_pari_coefficient_composition']=ff_pari_coefficient_composition

def ff_pari_extension_aut(mode,code,p,n,T,phi,slen,sflat,blen,bflat,alen,aflat):
    from pari_extension_aut import native_pari_extension_aut
    mode=int(mode)
    return native_pari_extension_aut(mode,int(code),int(p),int(n),int(T[0])if mode==2 else list(map(int,T)),
      _extension_decode(slen,sflat),int(phi[0])if mode==2 else list(map(int,phi)),_extension_decode(blen,bflat),_extension_decode(alen,aflat))
FUNCTIONS['ff_pari_extension_aut']=ff_pari_extension_aut

def ff_pari_extension_projection(mode,code,p,seed,n,T,alen,aflat,blen,bflat):
    from pari_extension_projection import native_pari_extension_projection
    return native_pari_extension_projection(int(mode),int(code),int(p),int(seed),int(n),list(map(int,T)),_extension_decode(alen,aflat),_extension_decode(blen,bflat))
FUNCTIONS['ff_pari_extension_projection']=ff_pari_extension_projection

def ff_pari_extension_minpoly(mode,code,p,seed,T,slen,sflat,xlen,xflat):
    from pari_extension_minpoly import native_pari_extension_minpoly
    return native_pari_extension_minpoly(int(mode),int(code),int(p),int(seed),list(map(int,T)),_extension_decode(slen,sflat),_extension_decode(xlen,xflat))
FUNCTIONS['ff_pari_extension_minpoly']=ff_pari_extension_minpoly

def ff_pari_prime_quotient_cache(op,p,n,T,a):
    from pari_prime_quotient_cache import native_pari_prime_quotient_cache
    return native_pari_prime_quotient_cache(int(op),int(p),int(n),list(map(int,T)),list(map(int,a)))
FUNCTIONS['ff_pari_prime_quotient_cache']=ff_pari_prime_quotient_cache

def ff_pari_extension_frobenius(mode,code,p,T,slen,sflat,alen,aflat):
    from pari_extension_frobenius import native_pari_extension_frobenius
    t=int(T[0])if mode==2 else list(map(int,T))
    return native_pari_extension_frobenius(int(mode),int(code),int(p),t,_extension_decode(slen,sflat),_extension_decode(alen,aflat))
FUNCTIONS['ff_pari_extension_frobenius']=ff_pari_extension_frobenius

def ff_pari_prime_linear(op,p,c,a,b):
    from pari_prime_linear import native_pari_prime_linear
    return native_pari_prime_linear(int(op),int(p),int(c),list(map(int,a)),list(map(int,b)))
FUNCTIONS['ff_pari_prime_linear']=ff_pari_prime_linear

def ff_pari_extension_root_count(mode,op,p,T,flen,fflat):
    from pari_extension_root_count import native_pari_extension_root_count
    t=int(T[0])if mode==2 else list(map(int,T))
    return native_pari_extension_root_count(int(mode),int(op),int(p),t,_extension_decode(flen,fflat))
FUNCTIONS['ff_pari_extension_root_count']=ff_pari_extension_root_count

def ff_pari_scalar_inverse(op,p,a):
    from pari_scalar_inverse import native_pari_scalar_inverse
    return native_pari_scalar_inverse(int(op),int(p),int(a))
FUNCTIONS['ff_pari_scalar_inverse']=ff_pari_scalar_inverse

def ff_pari_prime_observation(op,p,x,f):
    from pari_prime_observations import native_pari_prime_observation
    return native_pari_prime_observation(int(op),int(p),int(x),list(map(int,f)))
FUNCTIONS['ff_pari_prime_observation']=ff_pari_prime_observation

def ff_pari_fp_scalar(op,p,a,b,c):
    from pari_fp_scalars import native_pari_fp_scalar
    return native_pari_fp_scalar(int(op),int(p),int(a),int(b),int(c))
FUNCTIONS['ff_pari_fp_scalar']=ff_pari_fp_scalar

def ff_pari_fp_power(op,p,a,e):
    from pari_fp_power import native_pari_fp_power
    return native_pari_fp_power(int(op),int(p),int(a),int(e))
FUNCTIONS['ff_pari_fp_power']=ff_pari_fp_power

def ff_pari_fp_predicate(op,p,a,e):
    from pari_fp_predicates import native_pari_fp_predicate
    return native_pari_fp_predicate(int(op),int(p),int(a),int(e))
FUNCTIONS['ff_pari_fp_predicate']=ff_pari_fp_predicate

def ff_pari_qpoly(op,p,half,den,perm,L,M):
    from pari_qpoly import native_pari_qpoly
    return native_pari_qpoly(int(op),int(p),int(half),int(den),[int(x)for x in perm],[int(x)for x in L],[int(x)for x in M])
FUNCTIONS['ff_pari_qpoly']=ff_pari_qpoly

def ff_pari_galois_integer(op,a,b,k):
    from pari_galois_integer import native_pari_galois_integer
    return native_pari_galois_integer(int(op),int(a),int(b),int(k))
FUNCTIONS['ff_pari_galois_integer']=ff_pari_galois_integer

def ff_pari_zx_index(op,a,b,c,p):
    from pari_zx_index import native_pari_zx_index
    return native_pari_zx_index(int(op),[int(v)for v in a], [int(v)for v in b]if isinstance(b,list)else int(b),int(c),int(p))
FUNCTIONS['ff_pari_zx_index']=ff_pari_zx_index

def ff_pari_zp_lift(op,f,a,t,Q,p,e):
    from pari_zp_lift import native_pari_zp_lift
    return native_pari_zp_lift(int(op),*[list(map(int,v))for v in [f,a,t,Q]],int(p),int(e))
FUNCTIONS['ff_pari_zp_lift']=ff_pari_zp_lift

def ff_pari_vandermonde(op,L,d,p):
    from pari_vandermonde import native_pari_vandermonde
    return native_pari_vandermonde(int(op),list(map(int,L)),int(d),int(p))
FUNCTIONS['ff_pari_vandermonde']=ff_pari_vandermonde

def ff_pari_perm_vector(op,a,b,e):
    from pari_perm_vector import native_pari_perm_vector
    return native_pari_perm_vector(int(op),list(map(int,a)),list(map(int,b)),int(e))
FUNCTIONS['ff_pari_perm_vector']=ff_pari_perm_vector

def ff_pari_group_structure(op,n,g,go,h,ho,p):
    from pari_group_structure import native_pari_group_structure
    return native_pari_group_structure(int(op),int(n),*[list(map(int,v)) for v in (g,go,h,ho,p)])
FUNCTIONS['ff_pari_group_structure']=ff_pari_group_structure

def ff_pari_sympol(op,n,a,v,w,p):
    from pari_sympol import native_pari_sympol
    return native_pari_sympol(int(op),int(n),*[list(map(int,x)) for x in (a,v,w)],int(p))
FUNCTIONS['ff_pari_sympol']=ff_pari_sympol

def ff_pari_subgroup(op,n,a,b,p):
    from pari_subgroup import native_pari_subgroup
    return native_pari_subgroup(int(op),int(n),list(map(int,a)),list(map(int,b)),int(p))
FUNCTIONS['ff_pari_subgroup']=ff_pari_subgroup

def ff_pari_subgroup_bounds(op,n,a,b,p):
    from pari_subgroup_bounds import native_pari_subgroup_bounds
    return native_pari_subgroup_bounds(int(op),int(n),list(map(int,a)),list(map(int,b)),int(p))
FUNCTIONS['ff_pari_subgroup_bounds']=ff_pari_subgroup_bounds

def ff_pari_galois_validation(op,a,hasden,den):
    from pari_galois_validation import native_pari_galois_validation
    return native_pari_galois_validation(int(op),list(map(int,a)),int(hasden),int(den))
FUNCTIONS['ff_pari_galois_validation']=ff_pari_galois_validation

def ff_pari_galois_actions(*args):
    from pari_galois_actions import native_pari_galois_actions
    return native_pari_galois_actions(*args)
FUNCTIONS['ff_pari_galois_actions']=ff_pari_galois_actions

def ff_pari_symmetric_search(width,a,w,p):
    from pari_symmetric_search import native_pari_symmetric_search
    return native_pari_symmetric_search(int(width),list(map(int,a)),list(map(int,w)),int(p))
FUNCTIONS['ff_pari_symmetric_search']=ff_pari_symmetric_search

def ff_pari_findpsi(*args):
    from pari_findpsi import native_pari_findpsi
    return native_pari_findpsi(*args)
FUNCTIONS['ff_pari_findpsi']=ff_pari_findpsi

def ff_pari_galois_filter(*args):
    from pari_galois_filter import native_pari_galois_filter
    return native_pari_galois_filter(*args)
FUNCTIONS['ff_pari_galois_filter']=ff_pari_galois_filter

def ff_pari_galois_helpers(*args):
    from pari_galois_helpers import native_pari_galois_helpers
    return native_pari_galois_helpers(*args)
FUNCTIONS['ff_pari_galois_helpers']=ff_pari_galois_helpers

def ff_pari_galois_kernel(op,nr,nc,m,a,b):
    from pari_galois_kernel import native_pari_galois_kernel
    return native_pari_galois_kernel(int(op),int(nr),int(nc),list(map(int,m)),int(a),int(b))
FUNCTIONS['ff_pari_galois_kernel']=ff_pari_galois_kernel

def ff_pari_fixedfield_reprime(*args):
    from pari_fixedfield_reprime import native_pari_fixedfield_reprime
    return native_pari_fixedfield_reprime(*args)
FUNCTIONS['ff_pari_fixedfield_reprime']=ff_pari_fixedfield_reprime

def ff_pari_factor_stages(op,n,rounds,seed):
    from pari_factor_stages import native_pari_factor_stages
    return native_pari_factor_stages(int(op),int(n),int(rounds),int(seed))
FUNCTIONS['ff_pari_factor_stages']=ff_pari_factor_stages

def ff_pari_power_helpers(op,n,a,b,c):
    from pari_power_helpers import native_pari_power_helpers
    return native_pari_power_helpers(int(op),int(n),int(a),int(b),int(c))
FUNCTIONS['ff_pari_power_helpers']=ff_pari_power_helpers

def ff_pari_ecm(op,n,nbc,seed,b,rounds):
    from pari_ecm import native_pari_ecm
    return native_pari_ecm(int(op),int(n),int(nbc),int(seed),int(b),int(rounds))
FUNCTIONS['ff_pari_ecm']=ff_pari_ecm

def ff_pari_mpqs_relations(n,q,y1,r1,y2,r2,mode):
    from pari_mpqs_relations import native_pari_mpqs_relations
    return native_pari_mpqs_relations(int(n),int(q),int(y1),list(map(int,r1)),int(y2),list(map(int,r2)),int(mode))
FUNCTIONS['ff_pari_mpqs_relations']=ff_pari_mpqs_relations

def ff_pari_mpqs_sqrt(a,p):
    from pari_mpqs_sqrt import native_pari_mpqs_sqrt
    return native_pari_mpqs_sqrt(int(a),int(p))
FUNCTIONS['ff_pari_mpqs_sqrt']=ff_pari_mpqs_sqrt

def ff_pari_sparse_kernel(rows,seed,packed):
    M=[];i=0
    while i<len(packed):
        n=int(packed[i]);i+=1;M.append(packed[i:i+n]);i+=n
    from pari_sparse_kernel import native_pari_sparse_kernel
    return native_pari_sparse_kernel(int(rows),int(seed),[[int(v)for v in c]for c in M])
FUNCTIONS['ff_pari_sparse_kernel']=ff_pari_sparse_kernel

def ff_pari_mpqs_symbols(op,x,y):
    from pari_mpqs_symbols import native_pari_mpqs_symbols
    return native_pari_mpqs_symbols(int(op),int(x),int(y))
FUNCTIONS['ff_pari_mpqs_symbols']=ff_pari_mpqs_symbols

def ff_pari_mpqs_debug(op,n,y,relp,q,mode):
    from pari_mpqs_debug import native_pari_mpqs_debug
    return native_pari_mpqs_debug(int(op),int(n),int(y),[int(v)for v in relp],int(q),int(mode))
FUNCTIONS['ff_pari_mpqs_debug']=ff_pari_mpqs_debug

def ff_pari_mpqs_inverse(a,p):
    from pari_mpqs_inverse import native_pari_mpqs_inverse
    return native_pari_mpqs_inverse(int(a),int(p))
FUNCTIONS['ff_pari_mpqs_inverse']=ff_pari_mpqs_inverse

def ff_pari_mpqs_init(D,L,rounds,missing):
    from pari_mpqs_init import native_pari_mpqs_init
    return native_pari_mpqs_init(int(D),int(L),int(rounds),[int(v)for v in missing])
FUNCTIONS['ff_pari_mpqs_init']=ff_pari_mpqs_init

def ff_pari_mpqs_fb(n,size,want):
    from pari_mpqs_fb import native_pari_mpqs_fb
    return native_pari_mpqs_fb(int(n),int(size),int(want))
FUNCTIONS['ff_pari_mpqs_fb']=ff_pari_mpqs_fb

def ff_pari_mpqs_candidates(M,t,p):
    from pari_mpqs_candidates import native_pari_mpqs_candidates
    return native_pari_mpqs_candidates(int(M),int(t),int(p))
FUNCTIONS['ff_pari_mpqs_candidates']=ff_pari_mpqs_candidates

def ff_pari_mpqs_driver(n,seed):
    from pari_mpqs_driver import native_pari_mpqs_driver
    return native_pari_mpqs_driver(int(n),int(seed))
FUNCTIONS['ff_pari_mpqs_driver']=ff_pari_mpqs_driver

def ff_pari_mpqs_hash(size,packed):
    from pari_mpqs_hash import native_pari_mpqs_hash
    return native_pari_mpqs_hash(int(size),[int(v)for v in packed])
FUNCTIONS['ff_pari_mpqs_hash']=ff_pari_mpqs_hash

def ff_pari_mpqs_warning(n,packed):
    from pari_mpqs_warning import native_pari_mpqs_warning
    return native_pari_mpqs_warning(int(n),[int(v)for v in packed])
FUNCTIONS['ff_pari_mpqs_warning']=ff_pari_mpqs_warning

def ff_pari_mpqs_class_candidates(D,L,rounds,missing):
    from pari_mpqs_class_candidates import native_pari_mpqs_class_candidates
    return native_pari_mpqs_class_candidates(int(D),int(L),int(rounds),[int(v)for v in missing])
FUNCTIONS['ff_pari_mpqs_class_candidates']=ff_pari_mpqs_class_candidates


def ff_pari_zp_precision(op,f,a,t,p,e):
    from pari_zp_lift import native_pari_zp_lift
    return native_pari_zp_lift(int(op),list(map(int,f)),list(map(int,a)),list(map(int,t)),[],int(p),int(e))
FUNCTIONS['ff_pari_zp_precision']=ff_pari_zp_precision
