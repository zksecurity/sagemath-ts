#include <math.h>
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

static GEN table(GEN v,long mode,ulong p){GEN V=cgetg(lg(v),t_VEC);for(long i=1;i<lg(v);i++){GEN a=mode==2?binary(gel(v,i)):dense(gel(v,i),0);gel(V,i)=mode==1?ZXX_to_FlxX(a,p,1):a;}return V;}
int main(void){long mode,op,cache,inner;static char sp[16000],st[2000000],sq[2000000],sx[2000000],sv[2000000],ss[2000000];pari_init(256000000,500000);
while(scanf("%ld %ld %15999s %1999999s %1999999s %1999999s %1999999s %1999999s",&mode,&op,sp,st,sq,sx,sv,ss)==8){cache=(op%100)/10;inner=op/100;op%=10;pari_sp av=avma;
pari_CATCH(CATCH_ALL){char *e=pari_err2str(pari_err_last());print_error(e);pari_free(e);}
pari_TRY{GEN p=gp_read_str(sp),T,Q,X,V,S,R;
if(mode==2){T=packed(gp_read_str(st));Q=binary(gp_read_str(sq));X=binary(gp_read_str(sx));V=table(gp_read_str(sv),mode,2);S=binary(gp_read_str(ss));if(inner)T=F2x_get_red(T);if(cache)S=cache==1?mkvec2(F2xqX_invBarrett(S,T),S):F2xqX_get_red(S,T);
R=op?F2xqX_F2xqXQV_eval(Q,V,S,T):F2xqX_F2xqXQ_eval(Q,X,S,T);}
else{T=gtopolyrev(gp_read_str(st),1);Q=dense(gp_read_str(sq),0);X=dense(gp_read_str(sx),0);S=dense(gp_read_str(ss),0);V=table(gp_read_str(sv),mode,mode==1?itou(p):0);
if(mode==1){ulong q=itou(p);T=ZX_to_Flx(T,q);Q=ZXX_to_FlxX(Q,q,1);X=ZXX_to_FlxX(X,q,1);S=ZXX_to_FlxX(S,q,1);if(inner)T=Flx_get_red(T,q);if(cache)S=cache==1?mkvec2(FlxqX_invBarrett(S,T,q),S):FlxqX_get_red(S,T,q);
if(op==8)R=pol_0(0);
else if(op==9){T=Flx_get_red(T,q);S=FlxqX_get_red(S,T,q);(void)FlxqXQ_powers(X,(long)sqrt((double)degpol(Q)),S,T,q);R=pol_0(0);}
else R=op?FlxqX_FlxqXQV_eval(Q,V,S,T,q):FlxqX_FlxqXQ_eval(Q,X,S,T,q);}
else{if(inner)T=FpX_get_red(T,p);if(cache)S=cache==1?mkvec2(FpXQX_invBarrett(S,T,p),S):FpXQX_get_red(S,T,p);
R=op?FpXQX_FpXQXQV_eval(Q,V,S,T,p):FpXQX_FpXQXQ_eval(Q,X,S,T,p);}}
printf("OK ");print_outer(R,mode);putchar(10);}
pari_ENDCATCH;set_avma(av);fflush(stdout);}pari_close();return 0;}
