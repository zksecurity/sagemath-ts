/* Bundled PARI permuted HNF and knapsack lattice recognition. */
#include "pari.h"
#include "paripriv.h"
#include <stdio.h>
static GEN matrix_from(GEN v,long m,long n){GEN A=cgetg(n+1,t_MAT);for(long j=1;j<=n;j++){gel(A,j)=cgetg(m+1,t_COL);for(long i=1;i<=m;i++)gcoeff(A,i,j)=gel(v,(i-1)*n+j);}return A;}
static void value(GEN x){if(!x){printf("null");return;}if(typ(x)==t_INT){pari_printf("\"%Ps\"",x);return;}printf("[");for(long i=1;i<lg(x);i++){if(i>1)printf(",");if(typ(x)==t_VECSMALL)printf("\"%ld\"",x[i]);else value(gel(x,i));}printf("]");}
static void error(const char*s){printf("ERROR ");putchar(34);for(;*s;s++){unsigned char c=*s;if(c==34||c==92){putchar(92);putchar(c);}else if(c<32)printf("\\u%04x",c);else putchar(c);}putchar(34);putchar(10);}
int main(void){static char data[4000000];long op,m,n;pari_init(256000000,500000);
while(scanf("%ld %ld %ld %3999999s",&op,&m,&n,data)==4){pari_sp av=avma;
pari_CATCH(CATCH_ALL){char*s=pari_err2str(pari_err_last());error(s);pari_free(s);}
pari_TRY{GEN v=gp_read_str(data),A=matrix_from(v,m,n);
 if(op==4)value(ZM_hnf_knapsack(A));else if(op==5)value(hnfperm(A));
 else {GEN U=NULL,P=NULL,H=ZM_hnfperm(A,op&1?&U:NULL,op&2?&P:NULL);printf("[");value(H);printf(",");value(U);printf(",");value(P);printf("]");}
 printf("\n");}pari_ENDCATCH;set_avma(av);fflush(stdout);}pari_close();return 0;}
