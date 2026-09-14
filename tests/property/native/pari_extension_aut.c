#include <math.h>
#include <errno.h>
#include <stdlib.h>
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

int main(void){
 long mode,op;static char sp[16000],sn[16000],st[2000000],ss[2000000],sphi[2000000],sb[2000000],sa[2000000];pari_init(256000000,500000);
 while(scanf("%ld %ld %15999s %15999s %1999999s %1999999s %1999999s %1999999s %1999999s",&mode,&op,sp,sn,st,ss,sphi,sb,sa)==9){
  pari_sp av=avma;
  pari_CATCH(CATCH_ALL){char *e=pari_err2str(pari_err_last());print_error(e);pari_free(e);}
  pari_TRY{
   long cache=(op%100)/10,inner=op/100;op%=10;
   GEN p=gp_read_str(sp),N=gp_read_str(sn),T,S,phi,B,a,aut,R;
   long n;
   if(mode==1&&op==1)n=(long)itou(N);
   else{
    /* itos excludes LONG_MIN, although the C entry accepts that signed word. */
    char *end;errno=0;n=strtol(sn,&end,10);
    if(errno==ERANGE||*end)pari_err_OVERFLOW("t_INT-->long assignment");
   }
   if(mode==2){T=packed(gp_read_str(st));S=binary(gp_read_str(ss));phi=packed(gp_read_str(sphi));B=binary(gp_read_str(sb));a=binary(gp_read_str(sa));
    if(cache)S=cache==1?mkvec2(F2xqX_invBarrett(S,T),S):F2xqX_get_red(S,T);
    aut=op?mkvec3(phi,B,a):mkvec2(phi,B);
    R=op?F2xqXQ_auttrace(aut,n,S,T):F2xqXQ_autpow(aut,n,S,T);
   }else{
    T=gtopolyrev(gp_read_str(st),1);S=dense(gp_read_str(ss),0);phi=gtopolyrev(gp_read_str(sphi),1);B=dense(gp_read_str(sb),0);a=dense(gp_read_str(sa),0);
    if(mode==1){ulong q=itou(p);T=ZX_to_Flx(T,q);S=ZXX_to_FlxX(S,q,1);phi=ZX_to_Flx(phi,q);B=ZXX_to_FlxX(B,q,1);a=ZXX_to_FlxX(a,q,1);
     if(inner)T=Flx_get_red(T,q);if(cache)S=cache==1?mkvec2(FlxqX_invBarrett(S,T,q),S):FlxqX_get_red(S,T,q);
     aut=op==0?mkvec2(phi,B):op==1?mkvec2(B,a):mkvec3(phi,B,a);
     R=op==0?FlxqXQ_autpow(aut,n,S,T,q):op==1?FlxqXQ_auttrace(aut,(ulong)n,S,T,q):FlxqXQ_autsum(aut,n,S,T,q);
    }else{
     if(inner)T=FpX_get_red(T,p);if(cache)S=cache==1?mkvec2(FpXQX_invBarrett(S,T,p),S):FpXQX_get_red(S,T,p);
     aut=op==0?mkvec2(phi,B):op==1?mkvec2(B,a):mkvec3(phi,B,a);
     R=op==0?FpXQXQ_autpow(aut,n,S,T,p):op==1?FpXQXQ_auttrace(aut,n,S,T,p):FpXQXQ_autsum(aut,n,S,T,p);
    }
   }
   printf("OK [");
   if(op==1&&mode!=2)print_outer(gel(R,1),mode);
   else if(mode==2)pari_printf("%Ps",ZX_Z_eval(F2x_to_ZX(gel(R,1)),gen_2));
   else print_coefficient(mode==1?Flx_to_ZX(gel(R,1)):gel(R,1));
   for(long i=2;i<lg(R);i++){putchar(',');print_outer(gel(R,i),mode);}puts("]");
  }pari_ENDCATCH;
  set_avma(av);fflush(stdout);
 }
 pari_close();return 0;
}
