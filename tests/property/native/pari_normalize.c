#include "pari.h"
#include <stdio.h>
static void poly(GEN f){printf("[");for(long i=2;i<lg(f);i++){if(i>2)printf(",");pari_printf("%Ps",gel(f,i));}printf("]");}
int main(void){long op;static char sp[16000],sf[200000],sx[16000];pari_init(64000000,500000);
 while(scanf("%ld %15999s %199999s %15999s",&op,sp,sf,sx)==4){pari_sp av=avma;
 pari_CATCH(CATCH_ALL){char*e=pari_err2str(pari_err_last());printf("ERROR %s\n",e);pari_free(e);}
 pari_TRY{GEN p=gp_read_str(sp),f=gtopolyrev(gp_read_str(sf),0),x=gp_read_str(sx),r;
 if(op==0)r=FpX_normalize(f,p);else if(op==2)r=FpX_Fp_mul_to_monic(f,x,p);else{ulong q=itou(p);f=ZX_to_Flx(f,q);r=Flx_to_ZX(op==1?Flx_normalize(f,q):Flx_Fl_mul_to_monic(f,itou(x),q));}
 printf("OK ");poly(r);printf("\n");}pari_ENDCATCH;set_avma(av);fflush(stdout);
 }pari_close();return 0;}
