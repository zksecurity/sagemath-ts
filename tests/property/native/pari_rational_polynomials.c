/* Original QQ polynomial reduction, product, norm and trace. */
#include "pari.h"
#include "paripriv.h"
#include <stdio.h>
static void value(GEN x){if(typ(x)==t_INT||typ(x)==t_FRAC){pari_printf("\"%Ps\"",x);return;}printf("[");for(long i=1;i<lg(x);i++){if(i>1)printf(",");value(gel(x,i));}printf("]");}
static void error(const char*s){printf("ERROR ");putchar(34);for(;*s;s++){unsigned char c=*s;if(c==34||c==92){putchar(92);putchar(c);}else if(c<32)printf("\\u%04x",c);else putchar(c);}putchar(34);putchar(10);}
int main(void){static char data[4000000];pari_init(256000000,500000);
while(scanf("%3999999s",data)==1){pari_sp av=avma;
pari_CATCH(CATCH_ALL){char*s=pari_err2str(pari_err_last());error(s);pari_free(s);}
pari_TRY{GEN v=gp_read_str(data),a=gdiv(RgV_to_RgX(gel(v,2),0),gel(v,3)),b=gdiv(RgV_to_RgX(gel(v,4),0),gel(v,5)),T=gdiv(RgV_to_RgX(gel(v,6),0),gel(v,7)),z;long op=itos(gel(v,1));
 switch(op){case 0:z=RgX_rem(a,b);break;case 1:z=QX_mul(a,b);break;case 2:z=QX_ZX_rem(a,b);break;case 3:z=QXQ_mul(a,b,T);break;case 4:z=RgXQ_mul(a,b,T);break;case 5:z=QXQ_norm(a,T);break;case 6:z=RgXQ_norm(a,T);break;case 7:z=RgXQ_trace(a,T);break;case 8:z=gnorm(gmodulo(a,T));break;case 9:z=gtrace(gmodulo(a,T));break;default:z=ZX_rem(a,b);}
 if(op<=4||op==10){GEN d,w=Q_remove_denom(z,&d);value(mkvec2(RgX_to_RgV(w,lgpol(w)),d?d:gen_1));}else pari_printf("\"%Ps\"",z);printf("\n");
}pari_ENDCATCH;set_avma(av);fflush(stdout);}pari_close();return 0;}
