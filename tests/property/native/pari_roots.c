#include "pari.h"
#include "paripriv.h"
#include <stdio.h>
int main(void){long op;static char sp[16000],sf[2000000];pari_init(128000000,500000);
while(scanf("%ld %15999s %1999999s",&op,sp,sf)==3){pari_sp av=avma;
pari_CATCH(CATCH_ALL){char *e=pari_err2str(pari_err_last());printf("ERROR %s\n",e);pari_free(e);}
pari_TRY{GEN p=gp_read_str(sp),f=gtopolyrev(gp_read_str(sf),0),r;
if(op<2){r=op?Flx_roots(ZX_to_Flx(f,itou(p)),itou(p)):FpX_roots(f,p);printf("OK [");for(long i=1;i<lg(r);i++){if(i>1)printf(",");if(op)printf("%lu",uel(r,i));else pari_printf("%Ps",gel(r,i));}printf("]\n");}
else if(op==2)printf("OK %d\n",FpX_is_totally_split(f,p));
else if(op==3)printf("OK %ld\n",FpX_nbroots(f,p));
else if(op==4)printf("OK %d\n",Flx_is_totally_split(ZX_to_Flx(f,itou(p)),itou(p)));
else printf("OK %ld\n",Flx_nbroots(ZX_to_Flx(f,itou(p)),itou(p)));
}pari_ENDCATCH;set_avma(av);fflush(stdout);}pari_close();return 0;}
