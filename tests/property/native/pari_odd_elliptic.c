/* Direct bundled PARI FlxqE/FpXQE calls; base-p indices are transport only. */
#include "pari.h"
#include "paripriv.h"
#include <stdio.h>
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
  if(mode==2 || mode==3) {
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
  GEN R=mode==1?ellinf():ellmul(E,P,gel(v,7)),R2=mode==1?ellinf():FF_ellmul(E,P,gel(v,7));
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
  if(mode!=1) {
    GEN results[2]={R,R2};
    for(int j=0;j<2;j++) {
      printf(",");R=results[j];
      if(ell_is_inf(R))printf("[]");else{printf("[");ffprint(gel(R,1),p);printf(",");ffprint(gel(R,2),p);printf("]");}
    }
  }
  pari_printf(",[%ld,%Ps,%d]",ell_get_type(E),p,!FF_equal0(ell_get_disc(E)));
  printf("]\n");
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
      if(op>=20)modeltest(v);
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
    set_avma(av);fflush(stdout);
  }
  free(line);pari_close();return 0;
}
