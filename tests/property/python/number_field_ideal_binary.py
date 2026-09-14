import json
from sage.all import *
def nf_ideal_binary(op,cs,fn,fd,aa,ad,ac,bb,bd,bc):
 R=PolynomialRing(QQ,'x');K=NumberField(R([QQ(c)*fn/fd for c in cs]),'a');n=K.degree()
 def ideal(flat,d,count):return K.ideal([K(R([QQ(c)/d for c in flat[j*n:(j+1)*n]]))for j in range(count)])
 I=ideal(aa,ad,ac);J=ideal(bb,bd,bc)
 def fmt(L):
  if L.is_zero():return ['0',[]]
  rows=[a.list()for a in L.basis()];d=lcm([c.denominator()for row in rows for c in row]);H=matrix(ZZ,[[ZZ(c*d)for c in row]for row in rows]).hermite_form();return [str(L.norm()),[[str(QQ(c)/d)for c in row]for row in H.rows()]]
 try:
  if op==8:r=fmt(I*J)
  elif op==9:
   nf=K.pari_nf();H=nf.idealdiv(I.pari_hnf(),J.pari_hnf());r=fmt(K.ideal(0)) if H.ncols()==0 else fmt(K.ideal([K(nf.nfbasistoalg(H[j])) for j in range(H.ncols())]))
  elif op==0:r=fmt(I+J)
  elif op==1:r=fmt(I.intersection(J))
  elif op==2:r=I.is_coprime(J)
  elif op==3:r=I.divides(J)
  elif op==4:r=fmt(I/J)
  elif op==5:
   D=I.denominator();r=[fmt(D),D is I.denominator(),D.number_field() is K]
  elif op==6:
   N=I.numerator();r=[fmt(N),N is I.numerator(),N is I]
  else:r=[fmt(I.numerator()/I.denominator()),I.numerator().is_coprime(I.denominator())]
  error=None
 except Exception as e:r=None;error=type(e).__name__+': '+str(e)
 return json.dumps([r,error],separators=(',',':'))
