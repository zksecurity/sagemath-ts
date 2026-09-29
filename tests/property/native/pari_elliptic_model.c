/* Execute the bundled PARI curve coordinate kernels and scalar caller. */
#include "pari.h"
#include "paripriv.h"
#include <stdio.h>
int main(void) {
  char *line=NULL; size_t capacity=0;
  pari_init(64000000,500000);
  while (getline(&line,&capacity,stdin)>0) {
    pari_sp av=avma;
    pari_CATCH(CATCH_ALL) {
      char *error=pari_err2str(pari_err_last());
      printf("ERROR %s\n",error); pari_free(error);
    } pari_TRY {
      GEN v=gp_read_str(line), p=gel(v,1), cs=gel(v,2), P=gel(v,3), ch=gel(v,4), n=gel(v,5), result;
      long op=itos(gel(v,6));
      if (lg(P)==1) P=ellinf();
      if (op==0) result=FpE_changepoint(P,ch,p);
      else if(op==1) result=FpE_changepointinv(P,ch,p);
      else {
        GEN E=ellinit(cs,p,DEFAULTPREC);
        if(op==2) result=liftall(ellmul(E,P,n));
        else result=ellff_get_a4a6(E);
      }
      pari_printf("%Ps\n",result);
    } pari_ENDCATCH;
    set_avma(av); fflush(stdout);
  }
  free(line); pari_close(); return 0;
}
