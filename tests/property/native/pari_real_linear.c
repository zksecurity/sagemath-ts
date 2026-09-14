/* Bundled PARI real matrix and numeric branches. Precision uses the bit ABI. */
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

int main(void) {
 long op,mode,n,p,q,shift;static char data[4000000];pari_init(256000000,500000);
 while(scanf("%ld %ld %ld %ld %ld %ld %3999999s",&op,&mode,&n,&p,&q,&shift,data)==7){
  pari_sp av=avma;
  pari_CATCH(CATCH_ALL){char *s=pari_err2str(pari_err_last());error(s);pari_free(s);}
  pari_TRY{
   GEN v=gp_read_str(data),A=cgetg(n+1,t_MAT);
   for(long j=1;j<=n;j++){
    gel(A,j)=cgetg(n+1,t_COL);
    for(long i=1;i<=n;i++){
     GEN z=gel(v,(i-1)*n+j);
     if(!(mode&1) || signe(z)) {
      long bits=(mode&2)&&j==1?p:q;
      z=(!signe(z)&&(mode&4))?real_0_bit(shift-bits):shiftr(itor(z,bits),shift);
     }
     gcoeff(A,i,j)=z;
    }
   }
   if(op<3)value(op==0?qfgaussred_positive(A):op==1?RgM_Cholesky(A,p):RgM_inv_upper(A));
   else {
    GEN x=gcoeff(A,1,1),y=gcoeff(A,1,2);
    value(op==3?gmul(x,y):op==4?gsub(x,y):op==5?gneg(x):op==6?gdiv(x,y):ginv(x));
   }
   printf("\n");
  }pari_ENDCATCH;set_avma(av);fflush(stdout);
 }pari_close();return 0;
}
