#include "pari.h"
#include "paripriv.h"
#include <stdio.h>
static void poly(GEN f){printf("[");for(long i=2;i<lg(f);i++){if(i>2)printf(",");pari_printf("%Ps",gel(f,i));}printf("]");}
int main(void){long op;static char sp[16000],sf[2000000];pari_init(256000000,500000);
 while(scanf("%ld %15999s %1999999s",&op,sp,sf)==3){pari_sp av=avma;
  pari_CATCH(CATCH_ALL){char *e=pari_err2str(pari_err_last());printf("ERROR %s\n",e);pari_free(e);}
  pari_TRY{GEN p=gp_read_str(sp),f=gtopolyrev(gp_read_str(sf),0);
   if(op<2){GEN F=op==0?FpX_ddf(f,p):Flx_ddf(ZX_to_Flx(f,itou(p)),itou(p));printf("OK [");
    for(long i=1;i<lg(gel(F,1));i++){if(i>1)printf(",");printf("[%ld,",uel(gel(F,2),i));poly(op==0?gmael(F,1,i):Flx_to_ZX(gmael(F,1,i)));printf("]");}printf("]\n");}
   else if(op==5){GEN F=FpX_ddf(f,p),D=const_vecsmall(degpol(f),0);long nb=0;
    for(long i=1;i<lg(gel(F,1));i++){long d=uel(gel(F,2),i),count=degpol(gmael(F,1,i))/d;D[d]=count;nb+=count;}
    printf("OK [[0");for(long i=1;i<lg(D);i++)printf(",%ld",uel(D,i));printf("],%ld]\n",nb);}
   else if(op==2){long n=FpX_nbfact(f,p);printf("OK %ld\n",n);}
   else{long nb;GEN D=Flx_nbfact_by_degree(ZX_to_Flx(f,itou(p)),&nb,itou(p));printf("OK [[0");for(long i=1;i<lg(D);i++)printf(",%ld",uel(D,i));printf("],%ld]\n",nb);}
  }pari_ENDCATCH;set_avma(av);fflush(stdout);
 }pari_close();return 0;}
