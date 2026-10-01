/* Direct bundled PARI FlxqE/FpXQE calls; base-p indices are transport only. */
#include "pari.h"
#include "paripriv.h"
#include <stdio.h>
/* Compile the actual bundled Fp_ellcard body with observed backend calls. */
static GEN dispatch_calls;
static long dispatch_count,dispatch_cm;
static void dispatch_log(long tag,GEN a4,GEN a6,GEN p,long smallfact) {
  gel(dispatch_calls,++dispatch_count)=mkvec5(stoi(tag),a4,a6,p,stoi(smallfact));
}
static GEN dispatch_CM(GEN a4,GEN a6,GEN p) {
  dispatch_log(0,a4,a6,p,0);return dispatch_cm?stoi(777):NULL;
}
static GEN dispatch_SEA(GEN a4,GEN a6,GEN p,long smallfact) {
  dispatch_log(1,a4,a6,p,smallfact);return stoi(123);
}
static long dispatch_word(ulong a4,ulong a6,ulong p) {
  dispatch_log(2,utoi(a4),utoi(a6),utoi(p),0);return 456;
}
static GEN dispatch_big(GEN a4,GEN a6,GEN p) {
  dispatch_log(3,a4,a6,p,0);return stoi(789);
}
#define Fl_elltrace_naive Fl_elltrace
#define Fp_ellcard_CM dispatch_CM
#define Fp_ellcard_SEA dispatch_SEA
#define Fl_ellcard_Shanks dispatch_word
#define Fp_ellcard_Shanks dispatch_big
#include "pari_cardinality_dispatch.h"
#undef Fl_elltrace_naive
#undef Fp_ellcard_CM
#undef Fp_ellcard_SEA
#undef Fl_ellcard_Shanks
#undef Fp_ellcard_Shanks
#include "pari_constant_j.h"
static GEN polynomial(GEN n,GEN p,int word) {
  GEN v=cgetg(2+expi(addiu(n,1)),t_VEC);
  long k=1;
  while(signe(n)){GEN r; n=dvmdii(n,p,&r);gel(v,k++)=r;}
  setlg(v,k);GEN x=gtopolyrev(v,1);
  return word?ZX_to_Flx(x,itou(p)):x;
}
static GEN point(GEN v,GEN p,int word) {
  return lg(v)==1?ellinf():mkvec2(polynomial(gel(v,1),p,word),polynomial(gel(v,2),p,word));
}
static GEN multiply(GEN P,GEN n,GEN a,GEN T,GEN p,int word) {
  return word?FlxqE_mul(P,n,a,T,itou(p)):FpXQE_mul(P,n,a,T,p);
}
static void printpoint(GEN P,GEN p,int word) {
  if(ell_is_inf(P)){printf("[]");return;}
  GEN x=gel(P,1),y=gel(P,2);
  if(word){x=Flx_to_ZX(x);y=Flx_to_ZX(y);}
  pari_printf("[%Ps,%Ps]",ZX_Z_eval(x,p),ZX_Z_eval(y,p));
}
static GEN ffvalue(GEN x,GEN p,GEN fg,int integer) {
  return integer?x:Fq_to_FF(polynomial(x,p,0),fg);
}
static void ffprint(GEN x,GEN p) {pari_printf("%Ps",ZX_Z_eval(FF_to_FpXQ(x),p));}
static void modeltest(GEN v) {
  GEN p=gel(v,1),T=gtopolyrev(gel(v,2),1),fg=Tp_to_FF(T,p),cs=gel(v,3),rawP=gel(v,4);
  long mode=itos(gel(v,8)),encoding=itos(gel(gel(v,6),1));
  GEN a=cgetg(lg(cs),t_VEC),P=lg(rawP)==1?ellinf():mkvec2(ffvalue(gel(rawP,1),p,fg,encoding&2),ffvalue(gel(rawP,2),p,fg,encoding&2));
  for(long i=1;i<lg(cs);i++)gel(a,i)=ffvalue(gel(cs,i),p,fg,encoding&1);
  if(mode==2 || mode==3 || (mode>=6 && mode<=8 && lg(rawP)>1)) {
    if(mode==3)gel(a,2)=gneg(gdiv(gsqr(gel(a,1)),stoi(4)));
    GEN x=gel(P,1),y=gel(P,2),x2=gsqr(x);
    gel(a,5)=gsub(gadd(gsqr(y),gadd(gmul(gmul(gel(a,1),x),y),gmul(gel(a,3),y))),
      gadd(gmul(x2,x),gadd(gmul(gel(a,2),x2),gmul(gel(a,4),x))));
  }
  GEN E;
  if(mode==1) {
    E=cgetg(17,t_VEC);for(long i=1;i<=12;i++)gel(E,i)=gel(a,i);
    for(long i=13;i<=16;i++)gel(E,i)=gen_0;
    E=FF_ellinit(E,fg);
  }else E=ellinit(a,mode==4?p:fg,DEFAULTPREC);
  if(lg(E)==1){printf("OK null\n");return;}
  GEN R=ellinf(),R2=ellinf(),order=NULL;
  if(mode>=6 && mode<=8) {
    GEN o=gel(v,7);if(mode==7)o=Z_factor(o);else if(mode==8)o=mkvec2(o,Z_factor(o));
    order=FF_ellorder(E,P,o);
  }else if(mode!=1){R=ellmul(E,P,gel(v,7));R2=FF_ellmul(E,P,gel(v,7));}
  printf("OK [[");
  for(long i=1;i<=13;i++){if(i>1)printf(",");ffprint(gel(E,i),p);}
  GEN model=ellff_get_a4a6(E),ma=gel(model,1),ch=gel(model,3);
  int word=lgefint(p)==3;
  printf("],[");
  if(typ(ma)==t_VEC){printf("[");pari_printf("%Ps",ZX_Z_eval(Flx_to_ZX(gel(ma,1)),p));printf("]");}
  else pari_printf("%Ps",ZX_Z_eval(word?Flx_to_ZX(ma):ma,p));
  pari_printf(",%Ps,[",ZX_Z_eval(word?Flx_to_ZX(gel(model,2)):gel(model,2),p));
  for(long i=1;i<=4;i++){if(i>1)printf(",");GEN c=gel(ch,i);pari_printf("%Ps",ZX_Z_eval(word?Flx_to_ZX(c):c,p));}
  printf("]]");
  if(order) pari_printf(",%Ps",order);
  else if(mode!=1) {
    GEN results[2]={R,R2};
    for(int j=0;j<2;j++) {
      printf(",");R=results[j];
      if(ell_is_inf(R))printf("[]");else{printf("[");ffprint(gel(R,1),p);printf(",");ffprint(gel(R,2),p);printf("]");}
    }
  }
  pari_printf(",[%ld,%Ps,%d]",ell_get_type(E),p,!FF_equal0(ell_get_disc(E)));
  printf("]\n");
}
/* Generic order is the shared dependency of all three extension kernels. */
static GEN order_clones[32768], order_mod;
static long order_count;
static void order_log(long code,GEN a,GEN n) {
  if(order_count>=32768) pari_err_BUG("order oracle trace capacity");
  order_clones[order_count++]=gclone(mkvec3(stoi(code),a,n));
}
static GEN order_pow(void *unused,GEN a,GEN n) {
  (void)unused;order_log(1,a,n);return modii(mulii(a,n),order_mod);
}
static int order_identity(GEN a) {order_log(0,a,gen_0);return !signe(a);}
static void ordertest(GEN v) {
  GEN rows=gel(v,4),n=gel(v,7),o=n;
  long encoding=itos(gel(gel(v,6),1));
  if(encoding==1 || encoding==2) {
    GEN ps=cgetg(lg(rows),t_COL),es=cgetg(lg(rows),t_COL);
    for(long i=1;i<lg(rows);i++){gel(ps,i)=gel(gel(rows,i),1);gel(es,i)=gel(gel(rows,i),2);}
    o=mkmat2(ps,es);if(encoding==2)o=mkvec2(n,o);
  }else if(encoding==3)o=NULL;
  order_count=0;order_mod=gel(v,1);
  const struct bb_group group={NULL,order_pow,NULL,NULL,NULL,order_identity,NULL};
  GEN r=gen_order(modii(gel(v,3),order_mod),o,NULL,&group);
  GEN trace=cgetg(order_count+1,t_VEC);
  for(long i=0;i<order_count;i++)gel(trace,i+1)=order_clones[i];
  pari_printf("OK [%Ps,%Ps]\n",r,trace);
}
static void error(const char *s) {
  printf("ERROR \"");
  for(;*s;s++){unsigned char c=*s;if(c==34||c==92){putchar(92);putchar(c);}else if(c<32)printf("\\u%04x",c);else putchar(c);}
  printf("\"\n");
}
int main(void) {
  char *line=NULL;size_t capacity=0;pari_init(64000000,500000);
  while(getline(&line,&capacity,stdin)>0) {
    pari_sp av=avma;
    pari_CATCH(CATCH_ALL){char *s=pari_err2str(pari_err_last());error(s);pari_free(s);}
    pari_TRY {
      GEN v=gp_read_str(line),p=gel(v,1),T=gtopolyrev(gel(v,2),1),a,P,Q,ch,n=gel(v,7),R;
      int word=itos(gel(v,9));long op=itos(gel(v,10));int paired=op>=10;
      if(op==30)ordertest(v);
      else if(op==31){
        GEN E=ellinit(gel(v,3),p,DEFAULTPREC),PP=gel(v,4);
        GEN r=ellorder(E,lg(PP)==1?ellinf():PP,n);pari_printf("OK %Ps\n",r);
      }
      else if(op==32){GEN r=elltrace_extension(gel(v,3),itos(n),p);pari_printf("OK %Ps\n",r);}
      else if(op==33){GEN cs=gel(v,3),q=powiu(p,itos(n));GEN r=Fp_ffellcard(gel(cs,1),gel(cs,2),q,itos(n),p);pari_printf("OK %Ps\n",r);}
      else if(op==34){
        GEN cs=gel(v,3);dispatch_cm=itos(n);dispatch_calls=cgetg(8,t_VEC);dispatch_count=0;
        GEN r=audit_Fp_ellcard(modii(gel(cs,1),p),modii(gel(cs,2),p),p);
        setlg(dispatch_calls,dispatch_count+1);pari_printf("OK [%Ps,%Ps]\n",r,dispatch_calls);
      }
      else if(op==35){
        ulong pp=itou(p),aa=itou(gel(v,3)),pi=get_Fl_red(pp);long mode=itos(gel(v,8));
        GEN rawP=gel(v,4),rawQ=gel(v,5);
        GEN P=lg(rawP)==1?ellinf():ZV_to_Flv(rawP,pp),Q=lg(rawQ)==1?ellinf():ZV_to_Flv(rawQ,pp),r;
        if(mode==0)r=Flj_dbl_pre(P,aa,pp,pi);
        else if(mode==1)r=Flj_add_pre(P,Q,aa,pp,pi);
        else if(mode==2)r=Flj_neg(P,pp);
        else if(mode==3)r=Flj_mulu_pre(P,itou(n),aa,pp,pi);
        else if(mode==4)r=Fle_to_Flj(P);
        else if(mode==5)r=Flj_to_Fle_pre(P,pp,pi);
        else if(mode==6)r=Fle_dbl(P,aa,pp);
        else if(mode==7)r=Fle_add(P,Q,aa,pp);
        else if(mode==8)r=Fle_mulu(P,itou(n),aa,pp);
        else r=Fle_order(P,n,aa,pp);
        if(mode==9)pari_printf("OK %Ps\n",r);
        else if(mode>=5 && ell_is_inf(r))printf("OK []\n");
        else pari_printf("OK %Ps\n",Flv_to_ZV(r));
      }
      else if(op==36){
        GEN cs=gel(v,3);setrand(n);
        GEN a=modii(gel(cs,1),p),b=modii(gel(cs,2),p);
        /* Generated inputs satisfy the counter's nonsingularity precondition. */
        while(!signe(modii(addii(mului(4,powiu(a,3)),mului(27,sqri(b))),p))) b=modii(addiu(b,1),p);
        GEN r=Fp_ellcard(a,b,p);
        pari_printf("OK [%Ps,%Ps]\n",r,modii(getrand(),subiu(shifti(gen_1,127),1)));
      }
      else if(op==37){
        ulong pp=itou(p);long mode=itos(gel(v,8));setrand(n);
        GEN TT=ZX_to_Flx(T,pp),z=ZX_to_Flx(gtopolyrev(gel(v,3),1),pp),r;
        if(mode==0)r=Flxq_sqrt(z,TT,pp);
        else if(mode==1)r=Flxq_sqrt_pre(z,TT,pp,get_Fl_red(pp));
        else if(mode==2){GEN zz=gel(v,3);r=Fl2_sqrt_pre(mkvecsmall2(itou(gel(zz,1)),itou(gel(zz,2))),itou(gel(gel(v,2),1)),pp,get_Fl_red(pp));}
        else {
          GEN b=Flx_to_F2x(z),t=Flx_to_F2x(TT);
          if(mode==3)r=F2xq_sqrt(b,t);
          else if(mode==4)r=F2xq_autpow(b,itos(gel(gel(v,6),1)),t);
          else r=F2xq_sqrt_fast(b,F2xq_sqrt(polx_F2x(1),t),t);
          r=F2x_to_Flx(r);
        }
        if(r) {
          GEN out=mode==2?Flv_to_ZV(r):Flx_to_ZX(r);
          if(mode!=2)out=RgX_to_RgV(out,lg(out)-2);
          pari_printf("OK [%Ps,%Ps]\n",out,modii(getrand(),subiu(shifti(gen_1,127),1)));
        }else pari_printf("OK [null,%Ps]\n",modii(getrand(),subiu(shifti(gen_1,127),1)));
      }
      else if(op==38){
        GEN cs=gel(v,3);ulong pp=itou(p);setrand(n);
        GEN aa=ZX_to_Flx(gtopolyrev(gel(cs,1),1),pp),bb=ZX_to_Flx(gtopolyrev(gel(cs,2),1),pp);
        if(itos(gel(v,8)))aa=mkvec(aa);
        GEN r=Flxq_ellcard(aa,bb,ZX_to_Flx(T,pp),pp);
        pari_printf("OK [%Ps,%Ps]\n",r,modii(getrand(),subiu(shifti(gen_1,127),1)));
      }
      else if(op==39){
        long kind=itos(gel(v,8)),d=degpol(T);GEN j,aa,bb,z=FpX_red(gtopolyrev(gel(v,3),1),p);
        setrand(n);if(degpol(z)<0)z=pol_1(1);
        if(kind==0){j=gen_0;aa=pol_0(1);bb=z;}
        else if(kind==1){j=modsi(1728,p);aa=z;bb=pol_0(1);}
        else {
          j=modsi(kind,p);GEN j1728=modsi(1728,p);
          while(!signe(j)||equalii(j,j1728))j=modii(addiu(j,1),p);
          GEN g=Fp_div(j,Fp_sub(utoi(1728),j,p),p);
          aa=FpX_Fp_mul(FpXQ_sqr(z,T,p),Fp_mulu(g,3,p),p);
          bb=FpX_Fp_mul(FpXQ_powu(z,3,T,p),Fp_mulu(g,2,p),p);
        }
        GEN q=powiu(p,d),r;
        if(word)r=Flxq_ellcardj(ZX_to_Flx(aa,itou(p)),ZX_to_Flx(bb,itou(p)),itou(j),ZX_to_Flx(T,itou(p)),q,itou(p),d);
        else r=FpXQ_ellcardj(aa,bb,j,T,q,p,d);
        pari_printf("OK [%Ps,%Ps]\n",r,modii(getrand(),subiu(shifti(gen_1,127),1)));
      }
      else if(op>=20)modeltest(v);
      else {
      if(word)T=ZX_to_Flx(T,itou(p));
      a=polynomial(gel(v,3),p,word);
      if(itos(gel(v,8)))a=mkvec(a);
      P=point(gel(v,4),p,word);
      Q=paired?multiply(P,gel(v,5),a,T,p,word):point(gel(v,5),p,word);
      if(paired)op-=10;
      ch=cgetg(5,t_VEC);for(long i=1;i<=4;i++)gel(ch,i)=polynomial(gel(gel(v,6),i),p,word);
      if(op==0)R=word?FlxqE_add(P,Q,a,T,itou(p)):FpXQE_add(P,Q,a,T,p);
      else if(op==1)R=word?FlxqE_dbl(P,a,T,itou(p)):FpXQE_dbl(P,a,T,p);
      else if(op==2)R=word?FlxqE_neg(P,T,itou(p)):FpXQE_neg(P,T,p);
      else if(op==3)R=word?FlxqE_sub(P,Q,a,T,itou(p)):FpXQE_sub(P,Q,a,T,p);
      else if(op==4)R=multiply(P,n,a,T,p,word);
      else if(op==5)R=word?FlxqE_changepoint(P,ch,T,itou(p)):FpXQE_changepoint(P,ch,T,p);
      else R=word?FlxqE_changepointinv(P,ch,T,itou(p)):FpXQE_changepointinv(P,ch,T,p);
      printf("OK ");
      if(paired){printf("[");printpoint(Q,p,word);printf(",");printpoint(R,p,word);printf("]");}
      else printpoint(R,p,word);
      printf("\n");
      }
    } pari_ENDCATCH;
    for(long i=0;i<order_count;i++)gunclone(order_clones[i]);
    order_count=0;set_avma(av);fflush(stdout);
  }
  free(line);pari_close();return 0;
}
