/* Bundled PARI resultant, norm and square predicate; retained RNG is observable. */
#include "pari.h"
#include "paripriv.h"
#include <stdio.h>
static void print_error(const char *s){printf("ERROR \"");for(;*s;s++){unsigned char c=*s;if(c==34||c==92){putchar(92);putchar(c);}else if(c<32)printf("\\u%04x",c);else putchar(c);}printf("\"\n");}
int main(void){long op;static char sp[16000],sa[2000000],sb[2000000],ss[16000];pari_init(512000000,500000);
 while(scanf("%ld %15999s %1999999s %1999999s %15999s",&op,sp,sa,sb,ss)==5){pari_sp av=avma;
  pari_CATCH(CATCH_ALL){char *e=pari_err2str(pari_err_last());print_error(e);pari_free(e);}
  pari_TRY{GEN p=gp_read_str(sp),a=gtopolyrev(gp_read_str(sa),0),b=gtopolyrev(gp_read_str(sb),0),r=NULL,x=NULL;long yes=0;
   if(op==6||op==7||op==8)x=Fq_to_FF(a,ffgen(FpX_to_mod(b,p),0));
   setrand(gp_read_str(ss));
   switch(op){
    case 0:r=FpX_resultant(a,b,p);break;
    case 1:r=utoi(Flx_resultant(ZX_to_Flx(a,itou(p)),ZX_to_Flx(b,itou(p)),itou(p)));break;
    case 2:r=FpXQ_norm(a,b,p);break;
    case 3:r=utoi(Flxq_norm(ZX_to_Flx(a,itou(p)),ZX_to_Flx(b,itou(p)),itou(p)));break;
    case 4:yes=FpXQ_issquare(a,b,p);break;
    case 5:yes=Flxq_issquare(ZX_to_Flx(a,itou(p)),ZX_to_Flx(b,itou(p)),itou(p));break;
    case 6:r=FF_norm(x);break;
    case 7:yes=FF_issquare(x);break;
    case 8:r=FF_trace(x);break;
    case 9:r=FpXQ_trace(a,b,p);break;
    case 10:r=utoi(Flxq_trace(ZX_to_Flx(a,itou(p)),ZX_to_Flx(b,itou(p)),itou(p)));break;
    case 11:r=utoi(F2xq_trace(ZX_to_F2x(a),ZX_to_F2x(b)));break;
    default:pari_err_BUG("unknown field operation");
   }
   printf("{\"value\":");if(r)pari_printf("\"%Ps\"",r);else printf(yes?"true":"false");
   pari_printf(",\"state\":\"%Ps\"}\n",getrand());
  }pari_ENDCATCH;set_avma(av);fflush(stdout);
 }pari_close();return 0;
}
