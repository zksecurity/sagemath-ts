import json
from sage.all import *

def nf_ideal_audit(op,cs,fn,fd,flat,gd,count,exponent):
 R=PolynomialRing(QQ,'x');K=NumberField(R([QQ(c)*fn/fd for c in cs]),'a');n=K.degree();gens=[K(R([QQ(c)/gd for c in flat[j*n:(j+1)*n]]))for j in range(count)];I=K.ideal(gens)
 def ideal(J):
  if J.is_zero():return ['0',[]]
  v=[z.list()for z in J.basis()];d=lcm([c.denominator()for r in v for c in r]);H=matrix(ZZ,[[ZZ(c*d)for c in row]for row in v]).hermite_form();return [str(J.norm()),[[str(QQ(c)/d)for c in row]for row in H.rows()]]
 try:
  if op==0:r=ideal(I)
  elif op==1:r=ideal(~I)
  elif op==2:r=ideal(I*(~I))
  elif op==3:r=ideal(I**int(exponent))
  elif op==4:r=ideal(K.different()) # port order convenience maps to ambient maximal-order different
  elif op==5:r=ideal(~K.different())
  elif op==6:r=[ideal(K.different()),K.different() is K.different()]
  else:
   result=I**(int(exponent) if op==7 else ZZ(exponent));r=[ideal(result),result is I]
  error=None
 except Exception as e:r=None;error=type(e).__name__+': '+str(e)
 return json.dumps([r,error],separators=(',',':'))
