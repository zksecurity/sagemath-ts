#include "pari.h"
#include "paripriv.h"
#include <stdio.h>
int main(void) {
  static char sp[16000], st[200000], sx[200000], sq[16000], su[200000], sy[200000];
  pari_init(128000000, 500000);
  while (scanf("%15999s %199999s %199999s %15999s %199999s %199999s", sp,st,sx,sq,su,sy)==6) {
    pari_sp av=avma;
    pari_CATCH(CATCH_ALL) {
      char *e=pari_err2str(pari_err_last()); printf("ERROR %s\n",e); pari_free(e);
    } pari_TRY {
      GEN p=gp_read_str(sp),q=gp_read_str(sq);
      GEN T=gtopolyrev(gp_read_str(st),0),U=gtopolyrev(gp_read_str(su),0);
      GEN ff=ffgen(FpX_to_mod(T,p),0),gg=ffgen(FpX_to_mod(U,q),0);
      GEN x=Fq_to_FF(gtopolyrev(gp_read_str(sx),0),ff);
      GEN y=Fq_to_FF(gtopolyrev(gp_read_str(sy),0),gg);
      printf("%d\n",cmp_universal(x,y));
    } pari_ENDCATCH;
    set_avma(av); fflush(stdout);
  }
  pari_close(); return 0;
}
