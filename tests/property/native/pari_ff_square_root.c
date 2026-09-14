/* All root selection and random-state transitions execute bundled PARI. */
#include "pari.h"
#include "paripriv.h"
#include <stdio.h>
static void print_value(GEN r) {
  GEN v=gel(r,2);
  if (r[1]==t_FF_F2xq) v=F2x_to_ZX(v);
  else if (r[1]==t_FF_Flxq) v=Flx_to_ZX(v);
  printf("[");for(long i=2;i<lg(v);i++){if(i>2)printf(",");pari_printf("\"%Ps\"",gel(v,i));}printf("]");
}
int main(void){static char sp[16000],st[200000],sx[200000],ss[16000];
 pari_init(256000000,500000);
 while(scanf("%15999s %199999s %199999s %15999s",sp,st,sx,ss)==4){pari_sp av=avma;
  pari_CATCH(CATCH_ALL){char *e=pari_err2str(pari_err_last());printf("ERROR %s\n",e);pari_free(e);}
  pari_TRY{
   GEN p=gp_read_str(sp),T=gtopolyrev(gp_read_str(st),0),a=gtopolyrev(gp_read_str(sx),0);
   GEN g=ffgen(FpX_to_mod(T,p),0),x=Fq_to_FF(a,g),r;
   setrand(gp_read_str(ss)); long ok=FF_issquareall(x,&r);
   printf("{\"value\":");if(ok)print_value(r);else printf("null");
   pari_printf(",\"state\":\"%Ps\"}\n",getrand());
  }pari_ENDCATCH;set_avma(av);fflush(stdout);
 }pari_close();return 0;
}
