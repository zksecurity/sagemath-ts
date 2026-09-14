#include "pari.h"
#include "paripriv.h"
#include <stdio.h>
static GEN read_matrix(GEN v,long m,long n,long start){GEN A=cgetg(n+1,t_MAT);for(long j=1;j<=n;j++){gel(A,j)=cgetg(m+1,t_COL);for(long i=1;i<=m;i++)gcoeff(A,i,j)=gel(v,start+(i-1)*n+j);}return A;}
static void print_matrix(GEN A){printf("[");if(lg(A)>1)for(long i=1;i<lgcols(A);i++){if(i>1)printf(",");printf("[");for(long j=1;j<lg(A);j++){if(j>1)printf(",");pari_printf("\"%Ps\"",gcoeff(A,i,j));}printf("]");}printf("]");}
int main(void){long op,m,n;unsigned long p;static char data[2000000];pari_init(256000000,500000);
while(scanf("%ld %lu %ld %ld %1999999s",&op,&p,&m,&n,data)==5){pari_sp av=avma;
pari_CATCH(CATCH_ALL){printf("ERROR\n");}pari_TRY{
 GEN v=gp_read_str(data),r,d=gen_1;
 if(op==2){r=halfgcdii(gel(v,1),gel(v,2));printf("[");print_matrix(gel(r,1));pari_printf(",[\"%Ps\",\"%Ps\"]]",gmael(r,2,1),gmael(r,2,2));}
 else if(op==3){GEN a,b;if(Fp_ratlift(gel(v,1),gel(v,2),gel(v,3),gel(v,3),&a,&b))pari_printf("[\"%Ps\",\"%Ps\"]",a,b);else printf("null");}
 else{GEN A=read_matrix(v,m,n,0);
  if(op==0)r=Flm_to_ZM(Flm_adjoint(ZM_to_Flm(A,p),p));
  else if(op==1)r=ZM_inv(A,&d);
  else r=hnf_divscale(A,read_matrix(v,m,n,m*n),gel(v,2*m*n+1));
  if(!r)printf("null");else{printf("[");print_matrix(r);pari_printf(",\"%Ps\"]",d);}
 }
 putchar(10);
}pari_ENDCATCH;set_avma(av);fflush(stdout);}pari_close();return 0;}
