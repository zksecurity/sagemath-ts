/* Bundled PARI integer matrix pivots and solves. */
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
 long op,m,n,k;static char data[4000000],other[4000000];pari_init(256000000,500000);
 while(scanf("%ld %ld %ld %ld %3999999s %3999999s",&op,&m,&n,&k,data,other)==6){
  pari_sp av=avma;
  pari_CATCH(CATCH_ALL){char *s=pari_err2str(pari_err_last());error(s);pari_free(s);}
  pari_TRY{
   GEN A=matrix(gp_read_str(data),m,n),B=matrix(gp_read_str(other),m,k);
   if(op==0){long r;GEN d=ZM_pivots(A,&r);printf("[");if(!d)printf("null");else {printf("[");for(long i=1;i<lg(d);i++){if(i>1)printf(",");printf("%ld",d[i]);}printf("]");}printf(",%ld]",r);}
   else if(op==1)printf("%ld",ZM_rank(A));
   else{GEN X=ZM_gauss(A,B),den=NULL;if(!X)printf("null");else{X=Q_remove_denom(X,&den);printf("[");value(X);printf(",");value(den?den:gen_1);printf("]");}}
   printf("\n");
  }pari_ENDCATCH;set_avma(av);fflush(stdout);
 }
 pari_close();return 0;
}
