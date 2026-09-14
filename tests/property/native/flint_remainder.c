#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"
#include <stdio.h>
#include <string.h>
char * audit_remainder(int op,const char *as,const char *bs,const char *ads,const char *bds){
 fmpz_poly_t A,B,Q,R;fmpz_poly_init(A);fmpz_poly_init(B);fmpz_poly_init(Q);fmpz_poly_init(R);
 fmpz_poly_set_str(A,as);fmpz_poly_set_str(B,bs);ulong d=0;
 fmpz_t den,da,db;fmpz_init(den);fmpz_init(da);fmpz_init(db);fmpz_set_str(da,ads,10);fmpz_set_str(db,bds,10);
 if(op==0)fmpz_poly_pseudo_divrem_divconquer(Q,R,&d,A,B);
 else if(op==1)fmpz_poly_pseudo_divrem_basecase(Q,R,&d,A,B);
 else if(op==2)fmpz_poly_pseudo_rem(R,&d,A,B);
 else{
  fmpq_poly_t X,Y,Z;fmpq_poly_init(X);fmpq_poly_init(Y);fmpq_poly_init(Z);
  fmpq_poly_set_fmpz_poly(X,A);fmpq_poly_scalar_div_fmpz(X,X,da);
  fmpq_poly_set_fmpz_poly(Y,B);fmpq_poly_scalar_div_fmpz(Y,Y,db);
  fmpq_poly_rem(Z,X,Y);fmpq_poly_get_numerator(R,Z);fmpq_poly_get_denominator(den,Z);
  fmpq_poly_clear(X);fmpq_poly_clear(Y);fmpq_poly_clear(Z);
 }
 if(op!=3)fmpz_set_ui(den,d);
 char *q=fmpz_poly_get_str(Q),*r=fmpz_poly_get_str(R),*ds=fmpz_get_str(NULL,10,den);
 char *out=flint_malloc(strlen(q)+strlen(r)+strlen(ds)+3);sprintf(out,"%s\n%s\n%s",q,r,ds);
 flint_free(q);flint_free(r);flint_free(ds);fmpz_clear(den);fmpz_clear(da);fmpz_clear(db);
 fmpz_poly_clear(A);fmpz_poly_clear(B);fmpz_poly_clear(Q);fmpz_poly_clear(R);return out;
}
void audit_remainder_free(char *p){flint_free(p);}
