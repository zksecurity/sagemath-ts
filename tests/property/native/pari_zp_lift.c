#include "pari.h"
#include "paripriv.h"
#include <stdio.h>
static void record(GEN x);
struct iso_trace {void *E;GEN (*eval)(void*,GEN,GEN);GEN (*invd)(void*,GEN,GEN,GEN,long);};
static GEN iso_eval(void *E,GEN x,GEN q){
 struct iso_trace *d=E;GEN v=d->eval(d->E,x,q);record(mkvec4(gen_0,q,x,v));return v;
}
static GEN iso_invd(void *E,GEN V,GEN v,GEN q,long M){
 struct iso_trace *d=E;GEN r=d->invd(d->E,V,v,q,M);record(mkvecn(6,gen_1,q,stoi(M),V,v,r));return r;
}
static GEN audit_iso_newton(GEN x,GEN p,long n,void *E,GEN (*eval)(void*,GEN,GEN),GEN (*invd)(void*,GEN,GEN,GEN,long)){
 struct iso_trace d={E,eval,invd};return gen_ZpX_Newton(x,p,n,&d,iso_eval,iso_invd);
}
#include "pari_padic_frobenius.h"
/* Keep the upstream zero-exponent regression deterministic: stop immediately
 * before the original generic power's nonzero-exponent precondition is broken. */
static GEN audit_nonzero_power(GEN x,GEN n,void *E,GEN (*sqr)(void*,GEN),GEN (*mul)(void*,GEN,GEN)) {
 if(!signe(n))pari_err_BUG("zero exponent in Fp_pow2n");
 return gen_pow(x,n,E,sqr,mul);
}
#define gen_pow audit_nonzero_power
#include "pari_guarded_root.h"
#undef gen_pow
static void print_error(const char *s){
 printf("ERROR ");putchar(34);
 for(;*s;s++){unsigned char c=*s;if(c==34||c==92){putchar(92);putchar(c);}else if(c<32)printf("\\u%04x",c);else putchar(c);}
 putchar(34);putchar(10);
}



static void poly(GEN f){printf("[");for(long i=2;i<lg(f);i++){if(i>2)printf(",");pari_printf("%Ps",gel(f,i));}printf("]");}
/* Observe the public native lifting drivers without replacing their bodies. */
static GEN trace[32768];static long trace_count;
static GEN lift_a,lift_T,lift_ai;
static void record(GEN x){if(trace_count==32768)pari_err_BUG("precision trace capacity");trace[trace_count++]=gclone(x);}
static void tree(GEN x){
 if(typ(x)==t_POL){poly(x);return;}
 if(typ(x)==t_INT){pari_printf("%Ps",x);return;}
 printf("[");for(long i=1;i<lg(x);i++){if(i>1)printf(",");tree(gel(x,i));}printf("]");
}
static GEN newton_eval(void *E,GEN x,GEN q){
 (void)E;GEN f=FpX_Fp_sub(FpXQ_mul(x,FpX_red(lift_a,q),FpX_red(lift_T,q),q),gen_1,q);
 record(mkvec4(gen_0,q,x,f));return mkvec2(f,x);
}
static GEN newton_invd(void *E,GEN V,GEN v,GEN q,long M){
 (void)E;GEN r=FpXQ_mul(V,gel(v,2),FpX_red(lift_T,q),q);
 record(mkvecn(6,gen_1,q,stoi(M),V,gel(v,2),r));return r;
}
static GEN dixon_lin(void *E,GEN F,GEN d,GEN q){
 (void)E;GEN r=FpXQ_mul(gmael(F,1,1),d,gel(F,2),q);
 record(mkvec5(gen_0,q,F,d,r));return r;
}
static GEN dixon_invl(void *E,GEN d){
 GEN p=(GEN)E;GEN r=FpXQ_mul(lift_ai,d,FpX_red(lift_T,p),p);
 record(mkvec3(gen_1,d,r));return r;
}
static void audit_precision(long op,GEN f,GEN a,GEN T,GEN p,long e){
 GEN r,x;lift_a=f;lift_T=T;
 if((op>=13 && op<=15)||op==21){
   GEN aa=constant_coeff(a),exponent=(op==13||op==21)?constant_coeff(T):gen_2;
   GEN b=op==15?constant_coeff(f):addii(Fp_pow(aa,exponent,powiu(p,e)),mulii(p,constant_coeff(f)));
   r=op==21?audit_Zp_sqrtnlift(b,exponent,aa,p,e):op==13?Zp_sqrtnlift(b,exponent,aa,p,e):op==14?Zp_sqrtlift(b,aa,p,e):Zp_sqrt(b,p,e);
 }
 else if(op==22)r=Zp_exp(mulii(constant_coeff(f),equaliu(p,2)?stoi(4):p),p,e);
 else if(op==23)r=ZpXQ_log(ZX_Z_add(ZX_Z_mul(f,p),gen_1),T,p,e);
 else if(op==24)r=Zp_inv(constant_coeff(f),p,e);
 else if(op==25)r=Zp_div(constant_coeff(a),constant_coeff(f),p,e);
 else if(op==26)r=Zp_invlift(constant_coeff(f),e==1?constant_coeff(a):Zp_inv(constant_coeff(f),p,1),p,e);
 else if(op==27||op==28){
   x=cvtop(constant_coeff(f),p,e);
   if(op==28){setvalp(x,valp(x)+itos(remii(constant_coeff(a),stoi(5))));x=Qp_exp(x);}
   r=mkvec4(padic_u(x),stoi(valp(x)),stoi(precp(x)),padic_pd(x));
 }
 else if(op==29)r=ZpXQ_sqrtnorm(ZX_Z_add(ZX_Z_mul(f,p),gen_1),T,powiu(p,e),p,e);
 else if(op==32){
   GEN q=powiu(p,e),u=constant_coeff(f),v=polcoef_i(f,1,0),w=constant_coeff(a);
   GEN P=mkmat3(mkcol3(u,v,w),mkcol3(addis(w,1),subii(u,w),v),mkcol3(v,subis(w,1),u));
   GEN xp=FpXQ_powers(FpX_rem(FpX_red(f,q),T,q),2,T,q),yp=FpXQ_powers(FpX_rem(FpX_red(a,q),T,q),2,T,q);
   r=mkvec4(FpM_FpXV_bilinear(P,xp,yp,q),FpM_FpXQV_bilinear(P,xp,yp,T,q),FpXC_powderiv(xp,q),FpXV_FpC_mul(xp,gel(P,1),q));
 }
 else if(op==33)r=getc2(mkmat2(RgX_to_RgV(f,lgpol(f)),RgX_to_RgV(a,lgpol(a))),FpX_rem(pol_x(0),T,powiu(p,e)),T,powiu(p,e),itou(p),e);
 else if(op==35){GEN act;long dj;GEN phi=get_Kohel_polynomials(itou(p),&act,&dj);r=mkvec3(phi?phi:gen_0,act?act:gen_0,stoi(dj));}
 else if(op==37){GEN act;long dj;GEN phi=get_Kohel_polynomials(itou(p),&act,&dj),q=powiu(p,e);r=phi?getc2(act,FpX_rem(FpX_red(f,q),T,q),T,q,itou(p),e):gen_0;}
 else if(op==34||op==36){
   ulong pp=itou(p),pi=SMALL_ULONG(pp)?0:get_Fl_red(pp);GEN phi,act;long dj;
   if(op==36)phi=get_Kohel_polynomials(pp,&act,&dj);
   else {
     phi=zeromatcopy(pp+1,pp+1);
     for(long c=1;c<=pp+1;c++)for(long i=1;i<=pp+1;i++)gmael(phi,c,i)=mulii(p,addis(constant_coeff(a),c*i));
     gmael(phi,1,pp+1)=subis(gmael(phi,1,pp+1),1);gmael(phi,2,1)=addis(gmael(phi,2,1),1);
   }
   if(!phi)r=gen_0;
   else {
     GEN Tp=ZX_to_Flx(T,pp),q=powiu(p,e),lr=Flxq_lroot_pre(polx_Flx(0),Tp,pp,pi);
     GEN sqx=Flxq_powers_pre(lr,pp-1,Tp,pp,pi);
     x=e==1?f:Flx_to_ZX(Flx_rem(ZX_to_Flx(f,pp),Tp,pp));T=Flx_Teichmuller(Tp,pp,e);
     GEN Xm=FpXQ_powers(pol_xn(degpol(T),0),pp-1,T,q);
     r=lift_isogeny(phi,x,e,Xm,T,sqx,Tp,pp,pi);
   }
 }
 else if(op==30||op==31){
   ulong pp=itou(p),pi=SMALL_ULONG(pp)?0:get_Fl_red(pp);
   GEN Tp=ZX_to_Flx(T,pp),q=powiu(p,e);
   GEN lr=Flxq_lroot_pre(polx_Flx(0),Tp,pp,pi),sqx=Flxq_powers_pre(lr,pp-1,Tp,pp,pi);
   x=Flx_to_ZX(Flx_rem(ZX_to_Flx(f,pp),Tp,pp));
   if(op==31)r=Flx_to_ZX(Flxq_lroot_fast_pre(ZX_to_Flx(x,pp),sqx,Tp,pp,pi));
   else {
     T=Flx_Teichmuller(Tp,pp,e);
     GEN Xm=FpXQ_powers(pol_xn(degpol(T),0),pp-1,T,q);
     r=Teichmuller_lift(x,Xm,T,sqx,Tp,pp,pi,e);
   }
 }
 else if(op>=16 && op<=20){
   GEN q=powiu(p,e);
   if(op==17)T=Flx_Teichmuller(ZX_to_Flx(T,itou(p)),itou(p),e);
   if(!(op==19&&degpol(T)==1))f=FpX_rem(FpX_red(f,q),T,q);
   if(op==16)r=ZpXQ_frob_cyc(f,T,q,itou(p));
   else if(op==17)r=ZpXQ_frob(f,FpXQ_powers(pol_xn(degpol(T),0),itos(p)-1,T,q),T,q,itou(p));
   else if(op==18)r=ZpXQ_frob(f,cgetg(1,t_VEC),T,q,itou(p));
   else {
     if(lgpol(FpX_red(f,p))==0)f=pol_1(0);
     if(op==19)r=ZpXQ_norm_pcyc(f,T,q,p);
     else r=ZpXQ_sqrtnorm_pcyc(FpXQ_sqr(f,T,q),T,q,p,e);
   }
 }
 else if(op==11)r=Flx_Teichmuller(ZX_to_Flx(f,itou(p)),itou(p),e);
 else if(op==6)r=ZpXQ_inv(f,T,p,e);
 else if(op==7){x=e==1?a:ZpXQ_inv(f,T,p,1);r=ZpXQ_invlift(f,x,T,p,e);}
 else if(op==8)r=ZpXQ_div(a,f,T,powiu(p,e),p,e);
 else {
  lift_ai=ZpXQ_inv(f,T,p,1);
  if(op==9||op==12)r=gen_ZpX_Newton(lift_ai,p,e,NULL,newton_eval,newton_invd);
  else r=gen_ZpX_Dixon(mkvec2(mkvec(f),T),a,powiu(p,e),p,e,(void*)p,dixon_lin,dixon_invl);
 }
 printf("OK [");if(r)tree(r);else printf("null");printf(",[");
 for(long i=0;i<trace_count;i++){if(i)printf(",");tree(trace[i]);}printf("]]\n");
}
int main(void){long op,e;static char sf[2000000],sa[2000000],st[2000000],sq[2000000],sp[30000];pari_init(256000000,500000);
while(scanf("%ld %1999999s %1999999s %1999999s %1999999s %29999s %ld",&op,sf,sa,st,sq,sp,&e)==7){pari_sp av=avma;
pari_CATCH(CATCH_ALL){char*err=pari_err2str(pari_err_last());print_error(err);pari_free(err);}
pari_TRY{GEN f=gtopolyrev(gp_read_str(sf),0),a=gp_read_str(sa),t=gtopolyrev(gp_read_str(st),0),flat=gp_read_str(sq),p=gp_read_str(sp),z,Q;long count=0;
for(long i=1;i<lg(flat);){long n=itos(gel(flat,i));i+=n+1;count++;}Q=cgetg(count+1,t_VEC);
for(long i=1,j=1;i<lg(flat);j++){long n=itos(gel(flat,i++));GEN v=cgetg(n+1,t_VEC);for(long k=1;k<=n;k++)gel(v,k)=gel(flat,i++);gel(Q,j)=gtopolyrev(v,0);}
/* Match the initial native evaluation/inverse phases, but stop before
 * unchecked exact division when the root precondition is violated. */
if(op==5 && e>1){GEN S=gtopolyrev(a,0),TT=FpX_get_red(t,powiu(p,e)),q=sqri(p),Tq=FpXT_red(TT,q),Tq2=FpXT_red(Tq,p);
(void)FpXQ_inv(FpX_FpXQ_eval(FpX_deriv(f,p),S,Tq2,p),Tq2,p);
GEN value=FpX_FpXQ_eval(FpX_red(f,q),S,Tq,q);
for(long i=2;i<lg(value);i++)if(signe(modii(gel(value,i),p))){printf("GUARD Hensel lifting requires exact polynomial division\n");goto done;}}
if(op>=6){audit_precision(op,f,gtopolyrev(a,0),t,p,e);goto done;}
switch(op){case 0:z=ZpX_liftroot(f,gel(a,1),p,e);break;case 1:z=ZpX_liftroots(f,a,p,e);break;case 2:z=ZpX_roots(f,p,e);break;case 3:z=ZpX_liftfact(f,Q,powiu(p,e),p,e);break;case 4:z=bezout_lift_fact(f,Q,p,e);break;default:z=ZpX_ZpXQ_liftroot(f,gtopolyrev(a,0),t,p,e);}
printf("OK ");if(op==0)pari_printf("%Ps",z);else if(op==5)poly(z);else{printf("[");for(long i=1;i<lg(z);i++){if(i>1)printf(",");if(op==3||op==4)poly(gel(z,i));else pari_printf("%Ps",gel(z,i));}printf("]");}printf("\n");done:;}
pari_ENDCATCH;for(long i=0;i<trace_count;i++)gunclone(trace[i]);trace_count=0;set_avma(av);fflush(stdout);}pari_close();return 0;}
