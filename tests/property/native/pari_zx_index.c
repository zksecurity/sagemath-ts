#include "pari.h"
#include "paripriv.h"
#include <stdio.h>
static void print_error(const char *s){
 printf("ERROR ");putchar(34);
 for(;*s;s++){unsigned char c=*s;if(c==34||c==92){putchar(92);putchar(c);}else if(c<32)printf("\\u%04x",c);else putchar(c);}
 putchar(34);putchar(10);
}


static void poly(GEN f){printf("[");for(long i=2;i<lg(f);i++){if(i>2)printf(",");pari_printf("%Ps",gel(f,i));}printf("]");}
int main(void){long op;static char sa[2000000],sb[2000000],sc[30000],sp[30000];pari_init(256000000,500000);
while(scanf("%ld %1999999s %1999999s %29999s %29999s",&op,sa,sb,sc,sp)==5){pari_sp av=avma;
pari_CATCH(CATCH_ALL){char*e=pari_err2str(pari_err_last());print_error(e);pari_free(e);}
pari_TRY{GEN a=gtopolyrev(gp_read_str(sa),0),b=gtopolyrev(gp_read_str(sb),0),c=gp_read_str(sc),p=gp_read_str(sp),z;
switch(op){case 0:z=ZX_neg(a);break;case 1:z=ZX_Z_mul(a,c);break;case 2:z=ZX_deriv(a);break;case 3:z=ZX_resultant(a,b);break;case 4:z=ZX_disc(a);break;case 5:z=stoi(ZX_is_squarefree(a));break;case 6:z=indexpartial(a,NULL);break;case 7:z=indexpartial(a,c);break;case 8:z=ZpX_reduced_resultant(a,b,p,powiu(p,itou(c)));break;case 9:z=ZpX_reduced_resultant_fast(a,b,p,itos(c));break;case 10:{GEN U;GEN F=absZ_factor_limit_strict(c,0,&U);printf("OK [");for(long i=1;i<lg(gel(F,1));i++){if(i>1)printf(",");pari_printf("[%Ps,%Ps]",gcoeff(F,i,1),gcoeff(F,i,2));}printf("],");if(U)pari_printf("%Ps",U);else printf("0");printf("\n");goto done;}
case 11:z=ZX_gcd(a,b);break;
case 12:case 13:case 14:case 15:{GEN v=gp_read_str(sa);long nr=itos(c),nc=lg(v)==1?0:(lg(v)-1)/nr;GEN mat=cgetg(nc+1,t_MAT);for(long j=1;j<=nc;j++){GEN col=cgetg(nr+1,(op%2==0)?t_COL:t_VECSMALL);for(long i=1;i<=nr;i++){GEN e=gel(v,(j-1)*nr+i);if(op%2==0)gel(col,i)=e;else uel(col,i)=itou(e);}gel(mat,j)=col;}GEN pm=gp_read_str(sb);z=(op%2==0)?ZpM_echelon(mat,op>=14,p,pm):zlm_echelon(mat,op>=14,itou(p),itou(pm));if(!z){printf("OK null\n");goto done;}printf("OK [");for(long j=1;j<lg(z);j++){if(j>1)printf(",");printf("[");for(long i=1;i<lg(gel(z,j));i++){if(i>1)printf(",");if(op%2==0)pari_printf("%Ps",gcoeff(z,i,j));else printf("%lu",uel(gel(z,j),i));}printf("]");}printf("]\n");goto done;}
}
printf("OK ");if(op<=2||op==11)poly(z);else pari_printf("%Ps",z);printf("\n");done:;}
pari_ENDCATCH;set_avma(av);fflush(stdout);}pari_close();return 0;}
