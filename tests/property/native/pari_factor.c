#include "pari.h"
#include "paripriv.h"
#include <stdio.h>
static void poly(GEN f){printf("[");for(long i=2;i<lg(f);i++){if(i>2)printf(",");pari_printf("%Ps",gel(f,i));}printf("]");}
int main(void){long op;static char ss[16000],sp[16000],sf[2000000];pari_init(256000000,500000);
 while(scanf("%ld %15999s %15999s %1999999s",&op,ss,sp,sf)==4){pari_sp av=avma;
  pari_CATCH(CATCH_ALL){char *e=pari_err2str(pari_err_last());if(op==4)pari_printf("OK ERROR %s|%Ps\n",e,getrand());else printf("ERROR %s\n",e);pari_free(e);}
  pari_TRY{GEN p=gp_read_str(sp),f=gtopolyrev(gp_read_str(sf),0),F;pari_init_rand();setrand(gp_read_str(ss));
   if(op<=2||op==4){F=op!=1?FpX_factor(f,p):Flx_factor(ZX_to_Flx(f,itou(p)),itou(p));printf("OK [");
    for(long i=1;i<lg(gel(F,1));i++){if(i>1)printf(",");printf("[");poly(op!=1?gmael(F,1,i):Flx_to_ZX(gmael(F,1,i)));printf(",%ld]",uel(gel(F,2),i));}printf("]");}
   else{GEN r=FpX_factor_squarefree(f,p);printf("OK [");for(long i=1;i<lg(r);i++){if(i>1)printf(",");poly(gel(r,i));}printf("]");}
   pari_printf("|%Ps\n",getrand());
  }pari_ENDCATCH;set_avma(av);fflush(stdout);
 }pari_close();return 0;}
