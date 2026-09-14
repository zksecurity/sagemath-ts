"""Public curve ordering through bundled generator and rich-comparison bodies."""
from pathlib import Path
import json,textwrap
from sage.all import GF,QQ,PolynomialRing,EllipticCurve
import sage.schemes.elliptic_curves.weierstrass_morphism as wm
from weierstrass_roots import _isomorphisms

_path=Path(__file__).resolve().parents[3]/'reference/sage/src/sage/schemes/elliptic_curves/weierstrass_morphism.py'
_source=_path.read_text();_start=_source.index('    def _comparison_impl(');_end=_source.index('    def _eval(',_start)
_namespace=dict(vars(wm));exec(compile(textwrap.dedent(_source[_start:_end]),str(_path),'exec'),_namespace)
_comparison_impl=_namespace['_comparison_impl']

def ec_generic_isomorphism_order(p, degree, modulus, left, right, transform):
    try:
        if degree>1:
            K=GF(p**degree,'a',modulus=PolynomialRing(GF(p),'t')(modulus),impl='pari_ffelt')
            def decode(v):
                ds=[]
                for _ in range(degree):ds.append(v%p);v//=p
                return K(ds)
            encode=lambda v:str(v.integer_representation())
        else:
            K=GF(p) if p else QQ;decode=K;encode=lambda v:str(v)
        E=EllipticCurve(K,list(map(decode,left)))
        F=EllipticCurve(K,wm.baseWI(*map(decode,transform))(E.ainvs())) if transform else EllipticCurve(K,list(map(decode,right)))
        rows=list(_isomorphisms(E,F))
        old_gen=wm._isomorphisms;old_cmp=wm.WeierstrassIsomorphism._comparison_impl
        try:
            # Reuse only this call's freshly computed native roots. The public
            # methods still execute their own construction/sorting/error paths.
            wm._isomorphisms=lambda e,f:iter(rows) if e is E and f is F else _isomorphisms(e,f)
            wm.WeierstrassIsomorphism._comparison_impl=staticmethod(_comparison_impl)
            try:first={'value':list(map(encode,E.isomorphism_to(F).tuple()))}
            except Exception as e:first={'error':type(e).__name__,'message':str(e)}
            ordered=[[list(map(encode,w.tuple())),all(v.parent() is K for v in w.tuple())] for w in E.isomorphisms(F)]
            autos=[list(map(encode,w.tuple())) for w in E.automorphisms()] if E is F else None
            isomorphic=bool(E.is_isomorphic(F))
            return json.dumps({'value':[first,ordered,autos,isomorphic]},separators=(',',':'))
        finally:
            wm._isomorphisms=old_gen;wm.WeierstrassIsomorphism._comparison_impl=staticmethod(old_cmp)
    except Exception as e:return json.dumps({'error':type(e).__name__,'message':str(e)},separators=(',',':'))

def wm_isomorphism_comparisons(p,degree,modulus,left,transform):
    try:
        if degree>1:
            K=GF(p**degree,'a',modulus=PolynomialRing(GF(p),'t')(modulus),impl='pari_ffelt')
            def decode(v):
                ds=[]
                for _ in range(degree):ds.append(v%p);v//=p
                return K(ds)
        else:K=GF(p) if p else QQ;decode=K
        E=EllipticCurve(K,list(map(decode,left)))
        F=EllipticCurve(K,wm.baseWI(*map(decode,transform))(E.ainvs())) if transform else E
        morphisms=[wm.WeierstrassIsomorphism(E,t,F) for t in _isomorphisms(E,F)]
        def compare(a,b):
            results=[]
            for op in range(6):
                try:
                    v=_comparison_impl(a,b,op)
                    results.append(None if v is NotImplemented else bool(v))
                except Exception as e:results.append({'error':type(e).__name__,'message':str(e)})
            return results
        matrix=[[compare(a,b) for b in morphisms] for a in morphisms]
        identityE=wm.WeierstrassIsomorphism(E,(K.one(),K.zero(),K.zero(),K.zero()),E)
        identityF=wm.WeierstrassIsomorphism(F,(K.one(),K.zero(),K.zero(),K.zero()),F)
        domains=compare(identityE,identityF)
        codomains=compare(identityE,morphisms[0]) if morphisms else None
        invalid=[compare(None,identityE),compare(identityE,None),compare(None,7)]
        return json.dumps({'value':[matrix,domains,codomains,invalid]},separators=(',',':'))
    except Exception as e:return json.dumps({'error':type(e).__name__,'message':str(e)},separators=(',',':'))

def ec_is_isomorphic_arguments(kind):
    try:
        E=EllipticCurve(GF(5),[0,1])
        other=[None,7,[],object(),E][kind]
        return json.dumps({'value':bool(E.is_isomorphic(other))},separators=(',',':'))
    except Exception as e:return json.dumps({'error':type(e).__name__,'message':str(e)},separators=(',',':'))

def ec_isomorphism_parent_guards(left,right):
    try:
        def field(kind):
            if kind==0:return QQ
            if kind>=16:
                p,g={16:(7,2),17:(7,3),18:(7,2),19:(7,1),20:(2,0),21:(2,1)}[kind]
                return GF(p,modulus=PolynomialRing(GF(p),'x')([-g,1]))
            primes={1:2,2:3,3:5,4:7,12:5,13:2}
            if kind in primes:return GF(primes[kind])
            p,d,m,name={5:(2,2,[1,1,1],'a'),6:(2,2,[1,1,1],'b'),
                7:(3,2,[1,0,1],'a'),8:(3,2,[2,1,1],'a'),9:(3,2,[1,0,1],'b'),
                10:(5,2,[2,0,1],'a'),11:(3,2,[1,0,1],'a'),14:(2,3,[1,1,0,1],'a'),15:(3,3,[1,2,0,1],'a')}[kind]
            return GF(p**d,name,modulus=PolynomialRing(GF(p),'t')(m),impl='pari_ffelt')
        def curve(K):return EllipticCurve(K,[0,0,1,0,0] if K.characteristic()==2 else [0,0,0,1,0] if K.characteristic()==3 else [0,0,0,0,1])
        E,F=curve(field(left)),curve(field(right))
        a=wm.WeierstrassIsomorphism(E,(1,0,0,0),E);b=wm.WeierstrassIsomorphism(F,(1,0,0,0),F)
        comparison=[]
        for op in range(6):
            try:comparison.append(bool(_comparison_impl(a,b,op)))
            except Exception as e:comparison.append({'error':type(e).__name__,'message':str(e)})
        try:isomorphic={'value':bool(E.is_isomorphic(F))}
        except Exception as e:isomorphic={'error':type(e).__name__,'message':str(e)}
        return json.dumps({'value':[isomorphic,comparison]},separators=(',',':'))
    except Exception as e:return json.dumps({'error':type(e).__name__,'message':str(e)},separators=(',',':'))

def ec_isomorphism_root_trace(p,degree,modulus,coefficients,operation):
    from weierstrass_roots import _Polynomial
    try:
        if degree>1:
            K=GF(p**degree,'a',modulus=PolynomialRing(GF(p),'t')(modulus),impl='pari_ffelt')
            def decode(v):
                ds=[]
                for _ in range(degree):ds.append(v%p);v//=p
                return K(ds)
            encode=lambda v:str(v.integer_representation())
        else:K=GF(p);decode=K;encode=lambda v:str(v)
        E=EllipticCurve(K,list(map(decode,coefficients)));trace=[]
        old_gen=wm._isomorphisms;old_roots=_Polynomial.roots
        def roots(poly,multiplicities=False):
            trace.append([list(map(encode,poly.f.list())),bool(multiplicities)])
            return old_roots(poly,multiplicities=multiplicities)
        try:
            wm._isomorphisms=_isomorphisms;_Polynomial.roots=roots
            value=bool(E.is_isomorphic(E)) if operation==0 else list(map(encode,E.isomorphism_to(E).tuple()))
            return json.dumps({'value':[value,trace]},separators=(',',':'))
        finally:wm._isomorphisms=old_gen;_Polynomial.roots=old_roots
    except Exception as e:return json.dumps({'error':type(e).__name__,'message':str(e)},separators=(',',':'))
