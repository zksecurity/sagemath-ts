/* Bundled PARI matrix products and integer rescaling. Precision uses the bit ABI. */
#include "pari.h"
#include "paripriv.h"
#include <stdio.h>
#ifndef LONG_IS_64BIT
#error "Real matrix oracle requires 64-bit PARI"
#endif
static void value(GEN x) {
  if (!x) { printf("null"); return; }
  if (typ(x)==t_INT) { pari_printf("\"%Ps\"",x); return; }
  if (typ(x)==t_REAL) {
    GEN m=gen_0;
    if(signe(x))for(long i=2;i<lg(x);i++)m=addii(shifti(m,64),utoi((ulong)x[i]));
    printf("[%ld,\"%ld\",",signe(x),expo(x));pari_printf("\"%Ps\",%ld]",m,(lg(x)-2)*64);return;
  }
  printf("[");for(long i=1;i<lg(x);i++){if(i>1)printf(",");value(gel(x,i));}printf("]");
}
static void error(const char *s){
 printf("ERROR ");putchar(34);
 for(;*s;s++){unsigned char c=*s;if(c==34||c==92){putchar(92);putchar(c);}else if(c<32)printf("\\u%04x",c);else putchar(c);}
 putchar(34);putchar(10);
}


static GEN matrix(GEN v,long m,long n,long mode,long p,long shift){
 GEN A=cgetg(n+1,t_MAT);
 for(long j=1;j<=n;j++){
  gel(A,j)=cgetg(m+1,t_COL);
  for(long i=1;i<=m;i++){
   GEN z=gel(v,(i-1)*n+j);
   if(mode==1||(mode==2&&(i+j)%2))z=shiftr(itor(z,p),shift);
   if(mode==3 && signe(z))z=shiftr(itor(z,p),shift);
   gcoeff(A,i,j)=z;
  }
 }
 return A;
}
int main(void){
 long op,mode,mode2,m,k,n,p,q,shift;static char data[4000000],other[4000000];
 pari_init(256000000,500000);
 while(scanf("%ld %ld %ld %ld %ld %ld %ld %ld %ld %3999999s %3999999s",&op,&mode,&mode2,&m,&k,&n,&p,&q,&shift,data,other)==11){
  pari_sp av=avma;
  pari_CATCH(CATCH_ALL){char *s=pari_err2str(pari_err_last());error(s);pari_free(s);}
  pari_TRY{
   GEN A=matrix(gp_read_str(data),m,k,mode,p,shift),B=matrix(gp_read_str(other),k,n,mode2,q,-shift);
   GEN z;
   if(op==0)z=RgM_rescale_to_int(A);
   else if(op==1)z=gram_matrix(A);
   else if(op==2)z=RgM_mul(A,B);
   else if(op==3)z=RgV_dotsquare(gel(A,1));
   else if(op==4)z=RgV_dotproduct(gel(A,1),gel(A,2));
   else z=RgV_dotproduct(gel(A,1),gel(A,1));
   value(z);printf("\n");
  }pari_ENDCATCH;set_avma(av);fflush(stdout);
 }
 pari_close();return 0;
}
