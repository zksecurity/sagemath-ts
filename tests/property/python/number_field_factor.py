"""Sage ideal factor values and ordering; order uses a common nf_NOLLL basis."""
import json
from sage.all import *
from number_field_ideal_inputs import _field,_lattice
from pari_full_ideals import pari_full_ideal

def nf_ideal_factor(op,cs,flat,d,count):
    K=_field(cs);n=K.degree()
    gens=[K([QQ(c)/d for c in flat[i*n:(i+1)*n]])for i in range(int(count))];I=K.ideal(gens)
    if op==3:
        # Retain Sage Factorization.sort and NumberFieldIdeal._richcmp_. Only
        # replace cached HNFs to use the bundled PARI basis on both sides.
        rows=json.loads(pari_full_ideal(1,cs,flat,d,count));pairs=[]
        for prime,H,e in rows:
            P=K.ideal(1)
            P._NumberFieldIdeal__pari_hnf=pari(matrix(ZZ,H).transpose())
            pairs.append((P,ZZ(e)))
        F=Factorization(pairs)
        value=[[[[str(c)for c in col]for col in P.pari_hnf()],str(e)]for P,e in F]
    else:
        F=I.factor()if op==0 else K.factor(I)if op==1 else K.factor(gens)
        value=[[_lattice(P),str(e)]for P,e in F]
        value.sort(key=lambda row:json.dumps(row[0],separators=(',',':')))
    return json.dumps(value,separators=(',',':'))

def nf_factor_properties(cs,flat,d,count):
    from sage.rings.number_field.number_field_ideal import NumberFieldFractionalIdeal
    K=_field(cs);n=K.degree();I=K.ideal([K([QQ(c)/d for c in flat[i*n:(i+1)*n]])for i in range(int(count))])
    F=I.factor();product=K.ideal(1);rows=[]
    for P,e in F:
        product*=P**e
        rows.append([_lattice(P),str(e),P.is_prime(),isinstance(P,NumberFieldFractionalIdeal),str(P.ramification_index()),str(P.residue_class_degree()),str(I.valuation(P))])
    rows.sort(key=lambda row:json.dumps(row[0],separators=(',',':')))
    return json.dumps([F is I.factor(),F is K.factor(I),_lattice(product),rows],separators=(',',':'))
