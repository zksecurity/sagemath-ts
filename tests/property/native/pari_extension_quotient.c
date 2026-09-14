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
int main(void){long mode,op,cache,inner;static char sp[16000],sn[16000],st[2000000],sa[2000000],sb[2000000],ss[2000000];pari_init(256000000,500000);
while(scanf("%ld %ld %15999s %15999s %1999999s %1999999s %1999999s %1999999s",&mode,&op,sp,sn,st,sa,sb,ss)==8){cache=(op%100)/10;inner=op/100;op%=10;pari_sp av=avma;
pari_CATCH(CATCH_ALL){char *e=pari_err2str(pari_err_last());print_error(e);pari_free(e);}
pari_TRY{GEN p=gp_read_str(sp),n=gp_read_str(sn),T,A,B,S,R;
if(mode==2){T=packed(gp_read_str(st));A=binary(gp_read_str(sa));B=binary(gp_read_str(sb));S=binary(gp_read_str(ss));if(inner)T=F2x_get_red(T);if(cache)S=cache==1?mkvec2(F2xqX_invBarrett(S,T),S):F2xqX_get_red(S,T);
R=op==8?pol_0(0):op==0?F2xqXQ_mul(A,B,S,T):op==1?F2xqXQ_sqr(A,S,T):op==2?F2xqXQ_invsafe(A,S,T):op==3?F2xqXQ_inv(A,S,T):op==4?F2xqXQ_mul(A,F2xqXQ_inv(B,S,T),S,T):op==5?F2xqXQ_pow(A,n,S,T):F2xqXQ_powers(A,itos(n),S,T);}
else{T=gtopolyrev(gp_read_str(st),1);A=dense(gp_read_str(sa),0);B=dense(gp_read_str(sb),0);S=dense(gp_read_str(ss),0);
if(mode==1){ulong q=itou(p);T=ZX_to_Flx(T,q);A=ZXX_to_FlxX(A,q,1);B=ZXX_to_FlxX(B,q,1);S=ZXX_to_FlxX(S,q,1);if(inner)T=Flx_get_red(T,q);if(cache)S=cache==1?mkvec2(FlxqX_invBarrett(S,T,q),S):FlxqX_get_red(S,T,q);
R=op==8?pol_0(0):op==0?FlxqXQ_mul(A,B,S,T,q):op==1?FlxqXQ_sqr(A,S,T,q):op==2?FlxqXQ_invsafe(A,S,T,q):op==3?FlxqXQ_inv(A,S,T,q):op==4?FlxqXQ_div(A,B,S,T,q):op==5?FlxqXQ_pow(A,n,S,T,q):op==7?FlxqXQ_powu(A,itou(n),S,T,q):FlxqXQ_powers(A,itos(n),S,T,q);}
else {if(inner)T=FpX_get_red(T,p);if(cache)S=cache==1?mkvec2(FpXQX_invBarrett(S,T,p),S):FpXQX_get_red(S,T,p);
R=op==8?pol_0(0):op==0?FpXQXQ_mul(A,B,S,T,p):op==1?FpXQXQ_sqr(A,S,T,p):op==2?FpXQXQ_invsafe(A,S,T,p):op==3?FpXQXQ_inv(A,S,T,p):op==4?FpXQXQ_div(A,B,S,T,p):op==5?FpXQXQ_pow(A,n,S,T,p):FpXQXQ_powers(A,itos(n),S,T,p);}}
printf("OK ");if(!R)printf("null");else if(op==6){printf("[");for(long i=1;i<lg(R);i++){if(i>1)printf(",");print_outer(gel(R,i),mode);}printf("]");}else print_outer(R,mode);printf("\n");
}pari_ENDCATCH;set_avma(av);fflush(stdout);}pari_close();return 0;}
