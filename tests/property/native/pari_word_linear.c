/* Bundled PARI word matrix pivots and solves. */
#include "pari.h"
#include "paripriv.h"
#include <stdio.h>
#ifndef LONG_IS_64BIT
#error "Modular matrix oracle requires 64-bit PARI"
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

static GEN matrix(GEN v,long m,long n){
 GEN A=cgetg(n+1,t_MAT);
 for(long j=1;j<=n;j++){
  gel(A,j)=cgetg(m+1,t_COL);
  for(long i=1;i<=m;i++)gcoeff(A,i,j)=gel(v,(i-1)*n+j);
 }
 return A;
}
int main(void){
 long op,m,n,k;ulong p;static char data[4000000],other[4000000];pari_init(256000000,500000);
 while(scanf("%ld %ld %ld %ld %lu %3999999s %3999999s",&op,&m,&n,&k,&p,data,other)==7){
  pari_sp av=avma;
  pari_CATCH(CATCH_ALL){char *s=pari_err2str(pari_err_last());error(s);pari_free(s);}
  pari_TRY{
   GEN A=ZM_to_Flm(matrix(gp_read_str(data),m,n),p),B=ZM_to_Flm(matrix(gp_read_str(other),m,k),p);
   if(op==0){long r;GEN d=Flm_pivots(A,p,&r,0);printf("[");if(!d)printf("null");else {printf("[");for(long i=1;i<lg(d);i++){if(i>1)printf(",");printf("%ld",d[i]);}printf("]");}printf(",%ld]",r);}
   else{GEN X=Flm_gauss(A,B,p);value(X?Flm_to_ZM(X):NULL);}
   printf("\n");
  }pari_ENDCATCH;set_avma(av);fflush(stdout);
 }
 pari_close();return 0;
}
