/* Original nf_NOLLL prime decomposition, HNF and factorization. */
#include "pari.h"
#include "paripriv.h"
#include <stdio.h>
static void value(GEN x){if(!x){printf("null");return;}if(typ(x)==t_INT||typ(x)==t_FRAC){pari_printf("\"%Ps\"",x);return;}printf("[");for(long i=1;i<lg(x);i++){if(i>1)printf(",");value(gel(x,i));}printf("]");}
static void error(const char*s){printf("ERROR ");putchar(34);for(;*s;s++){unsigned char c=*s;if(c==34||c==92){putchar(92);putchar(c);}else if(c<32)printf("\\u%04x",c);else putchar(c);}putchar(34);putchar(10);}
int main(void){static char data[4000000];pari_init(256000000,500000);
while(scanf("%3999999s",data)==1){pari_sp av=avma;
 pari_CATCH(CATCH_ALL){char*s=pari_err2str(pari_err_last());error(s);pari_free(s);}
 pari_TRY{GEN v=gp_read_str(data),scale,T=ZX_to_monic(RgV_to_RgX(gel(v,2),0),&scale),nf=nfinit0(T,nf_NOLLL,DEFAULTPREC),out;long op=itos(gel(v,1));
 if(op==0){GEN L=idealprimedec_limit_f(nf,gel(v,3),itos(gel(v,4)));out=cgetg(lg(L),t_VEC);for(long i=1;i<lg(L);i++)gel(out,i)=mkvec2(gel(L,i),pr_hnf(nf,gel(L,i)));}
 else{GEN generators=gel(v,3),I=zeromat(0,0);for(long i=1;i<lg(generators);i++)I=idealadd(nf,I,gdiv(gsubst(RgV_to_RgX(gel(generators,i),0),0,gdiv(pol_x(0),scale)),gel(v,4)));
 if(op==1){GEN F=idealfactor(nf,I);out=cgetg(lg(gel(F,1)),t_VEC);for(long i=1;i<lg(out);i++){GEN P=gmael(F,1,i);gel(out,i)=mkvec3(P,pr_hnf(nf,P),gmael(F,2,i));}}
 else if(op==2)out=idealnumden(nf,I);
 else if(op==3)out=idealismaximal(nf,I);
 else if(op==4){GEN d,H=Q_remove_denom(I,&d);out=mkvec2(H,d?d:gen_1);}
 else{GEN d,a=algtobasis(nf,gdiv(gsubst(RgV_to_RgX(gel(generators,1),0),0,gdiv(pol_x(0),scale)),gel(v,4))),M=zk_multable(nf,Q_remove_denom(a,&d));out=mkvec2(zkmultable_inv(M),zkmultable_capZ(M));}
 }
 value(out);printf("\n");}pari_ENDCATCH;set_avma(av);fflush(stdout);}pari_close();return 0;}
