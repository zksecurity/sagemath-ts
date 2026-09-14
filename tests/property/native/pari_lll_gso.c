/* Bundled PARI adaptive QR/Cholesky precision. Internal functions from lll.c. */
#include "lll.c"
#include <stdio.h>
#ifndef LONG_IS_64BIT
#error "Adaptive Gram-Schmidt oracle requires 64-bit PARI"
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

int main(void){
 long op,mode,m,n,p,shift;static char data[4000000];pari_init(256000000,500000);
 while(scanf("%ld %ld %ld %ld %ld %ld %3999999s",&op,&mode,&m,&n,&p,&shift,data)==7){
  pari_sp av=avma;
  pari_CATCH(CATCH_ALL){char *s=pari_err2str(pari_err_last());error(s);pari_free(s);}
  pari_TRY{
   GEN v=gp_read_str(data),A=cgetg(n+1,t_MAT);
   for(long j=1;j<=n;j++){
    gel(A,j)=cgetg(m+1,t_COL);
    for(long i=1;i<=m;i++)gcoeff(A,i,j)=gel(v,(i-1)*n+j);
   }
   if(op==0)value(gramschmidt_upper(A));
   else if(op==1)value(gramschmidt_dynprec(A));
   else if(op==2)value(RgM_Cholesky_dynprec(gram_matrix(A)));
   else{
    if(mode)A=RgM_gtofp(A,p);
    if(mode&&shift)A=gmul2n(A,shift);
    value(mkvecn(7,stoi(drop(A)),stoi(potential(A)),stoi(spread(A)),stoi(condition_bound(A,0)),stoi(condition_bound(A,1)),stoi(GS_extraprec(A,0)),stoi(GS_extraprec(A,1))));
   }
   printf("\n");
  }pari_ENDCATCH;set_avma(av);fflush(stdout);
 }
 pari_close();return 0;
}
