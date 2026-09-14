/* Bundled PARI modular square roots and fused powering schedules. */
#include "pari.h"
#include "paripriv.h"
#include <stdio.h>
struct audit_fold {GEN p; ulong squares, fused;};
static GEN fold_square(void *data,GEN x){struct audit_fold *s=data;s->squares++;return modii(addis(mulis(x,3),1),s->p);}
static GEN fold_fused(void *data,GEN x){struct audit_fold *s=data;s->fused++;return modii(addis(mulis(x,7),2),s->p);}
int main(void){long op;static char sp[16000],sa[16000],sy[16000];pari_init(64000000,500000);
 while(scanf("%ld %15999s %15999s %15999s",&op,sp,sa,sy)==4){pari_sp av=avma;
  pari_CATCH(CATCH_ALL){char *error=pari_err2str(pari_err_last());printf("ERROR %s\n",error);pari_free(error);}
  pari_TRY{GEN p=gp_read_str(sp),a=gp_read_str(sa),y=gp_read_str(sy),r;
   if(op==3){struct audit_fold state={p,0,0};r=gen_pow_fold(a,y,&state,fold_square,fold_fused);pari_printf("OK [%Ps,%lu,%lu]\n",r,state.squares,state.fused);}
   else{if(op==0)r=Fp_sqrt(a,p);else if(op==1)r=Fp_sqrt_i(a,signe(y)?y:NULL,p);else{ulong z=Fl_sqrt(itou(a),itou(p));r=z==~0UL?NULL:utoi(z);}printf("OK ");if(r)pari_printf("%Ps",r);else printf("null");printf("\n");}
  }pari_ENDCATCH;set_avma(av);fflush(stdout);
 }pari_close();return 0;}
