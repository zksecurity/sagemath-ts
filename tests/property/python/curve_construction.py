"""Live constructor comparisons, preserving supplied coordinate models."""
import json
from sage.all import EllipticCurve, Integer, QQ

def ec_constructor_model(p,T,coefficients,encoding,entry,x):
    try:
        # Reuse the explicit field adapter without constructing a test curve.
        from sage.all import GF, PolynomialRing
        K=GF(p**(len(T)-1),'a',modulus=PolynomialRing(GF(p),'t')(T),impl='pari_ffelt') if T else (GF(p) if p else QQ)
        def encode(n):
            if encoding==1:
                if T:
                    ds=[]
                    for _ in range(len(T)-1): ds.append(n%p); n//=p
                    return K(ds)
                return K(n)
            if encoding==2: return Integer(n)
            if encoding==3: return QQ(n)/2
            if encoding==4: return str(n)
            return n
        cs=list(map(encode,coefficients))
        if entry==2:
            from sage.schemes.elliptic_curves.ell_generic import EllipticCurve_generic
            E=EllipticCurve_generic(K,cs)
        else: E=EllipticCurve(K,cs)
        # Caller-owned input is not part of Sage's immutable defining tuple.
        if cs: cs[-1]=K(0)
        value=[str(E),list(map(str,E.a_invariants())),str(E.discriminant()),str(E.j_invariant()),
               all(c.parent() is K for c in E.a_invariants())]
        if entry!=2:
            points=E.lift_x(K(x),all=True)
            value.append([[str(P),str(2*P),str(-P)] for P in points])
        return json.dumps({'value':value},separators=(',',':'))
    except Exception as e: return json.dumps({'error':type(e).__name__,'message':str(e)},separators=(',',':'))
