/* Unmodified bundled PARI QR routines. Real precision uses the bit ABI. */
#include "pari.h"
#include "paripriv.h"
#include <stdio.h>
#ifndef LONG_IS_64BIT
#error "QR oracle requires 64-bit PARI"
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
 long op,mode,m,n,p,q,shift;static char data[4000000];pari_init(256000000,500000);
 while(scanf("%ld %ld %ld %ld %ld %ld %ld %3999999s",&op,&mode,&m,&n,&p,&q,&shift,data)==8){
  pari_sp av=avma;
  pari_CATCH(CATCH_ALL){char *s=pari_err2str(pari_err_last());error(s);pari_free(s);}
  pari_TRY{
   GEN v=gp_read_str(data),A=cgetg(n+1,t_MAT);
   for(long j=1;j<=n;j++){
    gel(A,j)=cgetg(m+1,t_COL);
    for(long i=1;i<=m;i++){
     GEN z=gel(v,(i-1)*n+j);
     if(mode==1||(mode==2&&(i+j)%2))z=shiftr(itor(z,q),shift);
     gcoeff(A,i,j)=z;
    }
   }
   if(op==0){GEN B=NULL,Q=NULL,L=NULL;long s=QR_init(A,&B,&Q,&L,p);
    printf("[%ld,",s);
    /* QR_init allocates Q with n-1 public entries. For tall matrices its
       final scratch write occupies B's header: print B by the known count. */
    if(B){printf("[");for(long i=1;i<=n;i++){if(i>1)printf(",");value(gel(B,i));}printf("]");}else printf("null");
    printf(",");value(Q);printf(",");value(L);printf("]");
   }else value(op==1?R_from_QR(A,p):gaussred_from_QR(A,p));
   printf("\n");
  }pari_ENDCATCH;set_avma(av);fflush(stdout);
 }pari_close();return 0;
}
