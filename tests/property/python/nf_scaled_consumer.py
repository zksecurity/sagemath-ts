import json
from sage.all import *
from sage.schemes.elliptic_curves.gal_reps_number_field import Frobenius_filter

def nf_scaled_consumer(op,cs,num,den):
 R=PolynomialRing(QQ,'x');K=NumberField(R([QQ(c)*num/den for c in cs]),'a');a=K.gen()
 try:
  if op==0:r=[str(K.different().norm()),str((~K.different()).norm())]
  else:r=list(map(str,Frobenius_filter(EllipticCurve(K,[0,0,0,a+1,a]),[2,3,5,7,11,13])))
  error=None
 except Exception as e:r=None;error=type(e).__name__+': '+str(e)
 return json.dumps([r,error],separators=(',',':'))
