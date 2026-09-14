"""Native Sage ideal lattices and cached generator/module results."""
import json
from sage.all import NumberField, PolynomialRing, QQ, ZZ, lcm, matrix

def nf_ideal_basis(op,cs,fn,fd,flat,gd,count,_exponent):
    R=PolynomialRing(QQ,'x'); K=NumberField(R([QQ(c)*fn/fd for c in cs]),'a'); n=K.degree()
    I=K.ideal([K(R([QQ(c)/gd for c in flat[j*n:(j+1)*n]])) for j in range(count)])
    def lattice(basis):
        rows=[z.list() for z in basis if z]
        if not rows:return []
        d=lcm([c.denominator() for row in rows for c in row])
        H=matrix(ZZ,[[ZZ(c*d) for c in row] for row in rows]).hermite_form()
        return [[str(QQ(c)/d) for c in row] for row in H.rows() if any(row)]
    try:
        if op==12:result=I.free_module() is I.free_module()
        elif op in [9,10]:
            basis=I.integral_basis() if op==9 else [K(v) for v in I.free_module().basis()]
            result=[len(basis),lattice(basis)]
        else:
            pair=I.gens_two();a,b=pair
            result=[str(a),lattice(K.ideal(a,b).basis())==lattice(I.basis()),b.is_zero(),a.parent() is K,b.parent() is K,I.gens_two() is pair]
        error=None
    except Exception as e:result=None;error=type(e).__name__+': '+str(e)
    return json.dumps([result,error],separators=(',',':'))
