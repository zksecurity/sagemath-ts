import json
from sage.all import *
def nf_ideal_construction(cs,fn,fd,flat,den,count):
 R=PolynomialRing(QQ,'x');K=NumberField(R([QQ(c)*fn/fd for c in cs]),'a');n=K.degree();gens=[K(R([QQ(c)/den for c in flat[j*n:(j+1)*n]]))for j in range(count)]
 try:
  I=K.ideal(*gens);r=[[[str(c)for c in a.list()]for a in I.gens()],I.ngens(),I is K.ideal(*gens),I is K.ideal(0)];error=None
 except Exception as e:r=None;error=type(e).__name__+': '+str(e)
 return json.dumps([r,error],separators=(',',':'))
