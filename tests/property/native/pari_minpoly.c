/* Bundled PARI minimal-polynomial result and consumed random state. */
#include "pari.h"
#include "paripriv.h"
#include <stdio.h>
static void print_poly(GEN f){printf("[");for(long i=2;i<lg(f);i++){if(i>2)printf(",");pari_printf("%Ps",gel(f,i));}printf("]");}
int main(void){long op;static char ss[16000],sp[16000],sx[2000000],st[2000000];pari_init(256000000,500000);
 while(scanf("%ld %15999s %15999s %1999999s %1999999s",&op,ss,sp,sx,st)==5){pari_sp av=avma;
  pari_CATCH(CATCH_ALL){char *error=pari_err2str(pari_err_last());printf("ERROR %s\n",error);pari_free(error);}
  pari_TRY{GEN p=gp_read_str(sp),x=gtopolyrev(gp_read_str(sx),0),T=gtopolyrev(gp_read_str(st),0),r;pari_init_rand();setrand(gp_read_str(ss));
   if(op==0)r=FpXQ_minpoly(x,T,p);else{ulong q=itou(p);r=Flx_to_ZX(Flxq_minpoly(ZX_to_Flx(x,q),ZX_to_Flx(T,q),q));}
   printf("OK ");print_poly(r);pari_printf("|%Ps\n",getrand());
  }pari_ENDCATCH;set_avma(av);fflush(stdout);
 }pari_close();return 0;}
