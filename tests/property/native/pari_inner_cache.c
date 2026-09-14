#include "pari.h"
#include "paripriv.h"
#include <stdio.h>
static GEN dense(GEN v,long var){GEN x=cgetg(lg(v)+1,t_POL);x[1]=evalsigne(1)|evalvarn(var);for(long i=1;i<lg(v);i++)gel(x,i+1)=typ(gel(v,i))==t_VEC?gtopolyrev(gel(v,i),1):gel(v,i);return FpXX_renormalize(x,lg(x));}
static GEN packed(GEN n){if(!signe(n))return pol0_F2x(evalvarn(1));long d=expi(n);GEN v=cgetg(d+2,t_VEC);for(long i=0;i<=d;i++)gel(v,i+1)=bittest(n,i)?gen_1:gen_0;return ZX_to_F2x(gtopolyrev(v,1));}
static GEN binary(GEN v){GEN x=cgetg(lg(v)+1,t_POL);x[1]=evalsigne(1)|evalvarn(0);for(long i=1;i<lg(v);i++)gel(x,i+1)=packed(gel(v,i));return F2xX_renormalize(x,lg(x));}
static void print_coefficient(GEN c){if(typ(c)==t_INT){pari_printf("%Ps",c);return;}printf("[");for(long j=2;j<lg(c);j++){if(j>2)printf(",");pari_printf("%Ps",gel(c,j));}printf("]");}
static void print_outer(GEN R,long mode){printf("[");for(long i=2;i<lg(R);i++){if(i>2)printf(",");if(mode==2)pari_printf("%Ps",ZX_Z_eval(F2x_to_ZX(gel(R,i)),gen_2));else print_coefficient(mode==1?Flx_to_ZX(gel(R,i)):gel(R,i));}printf("]");}


static void print_error(const char *s){
 printf("ERROR ");putchar(34);
 for(;*s;s++){unsigned char c=*s;if(c==34||c==92){putchar(92);putchar(c);}else if(c<32)printf("\\u%04x",c);else putchar(c);}
 putchar(34);putchar(10);
}

static GEN powers(GEN a){GEN V=cgetg(lg(a),t_VEC);for(long i=1;i<lg(a);i++)gel(V,i)=gtopolyrev(gel(a,i),0);return V;}
int main(void){long op,n;static char sp[16000],st[2000000],sq[2000000],sx[2000000],sv[2000000];pari_init(256000000,500000);
while(scanf("%ld %15999s %ld %1999999s %1999999s %1999999s %1999999s",&op,sp,&n,sq,sx,sv,st)==7){pari_sp av=avma;
pari_CATCH(CATCH_ALL){char *e=pari_err2str(pari_err_last());print_error(e);pari_free(e);}
pari_TRY{GEN p=gp_read_str(sp),T=gtopolyrev(gp_read_str(st),0),Q=gtopolyrev(gp_read_str(sq),0),x=gtopolyrev(gp_read_str(sx),0),V=powers(gp_read_str(sv)),R;
if(op==0)R=FpXQ_powers(x,n,T,p);
else if(op==1)R=FlxV_to_ZXV(Flxq_powers(ZX_to_Flx(x,itou(p)),n,ZX_to_Flx(T,itou(p)),itou(p)));
else if(op==2)R=FpX_FpXQV_eval(Q,V,T,p);
else if(op==3)R=Flx_to_ZX(Flx_FlxqV_eval(ZX_to_Flx(Q,itou(p)),ZXV_to_FlxV(V,itou(p)),ZX_to_Flx(T,itou(p)),itou(p)));
else if(op==4)R=FpX_FpXQ_eval(Q,x,T,p);
else R=Flx_to_ZX(Flx_Flxq_eval(ZX_to_Flx(Q,itou(p)),ZX_to_Flx(x,itou(p)),ZX_to_Flx(T,itou(p)),itou(p)));
printf("OK ");if(op<2){putchar('[');for(long i=1;i<lg(R);i++){if(i>1)putchar(',');print_coefficient(gel(R,i));}putchar(']');}else print_coefficient(R);putchar(10);}
pari_ENDCATCH;set_avma(av);fflush(stdout);}pari_close();return 0;}
