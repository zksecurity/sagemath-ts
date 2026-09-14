/* Exercise original static prime finalization and public idealval, with native nf data. */
#include "base2.c"
#include <stdio.h>
static void value(GEN x){
 if(typ(x)==t_INT||typ(x)==t_FRAC||typ(x)==t_INFINITY){pari_printf("\"%Ps\"",x);return;}
 printf("[");for(long i=1;i<lg(x);i++){if(i>1)printf(",");value(gel(x,i));}printf("]");
}
static void error(const char*s){printf("ERROR ");putchar(34);for(;*s;s++){unsigned char c=*s;if(c==34||c==92){putchar(92);putchar(c);}else if(c<32)printf("\\u%04x",c);else putchar(c);}putchar(34);putchar(10);}
static int compare_hnf(void *E,GEN a,GEN b){(void)E;for(long j=1;j<lg(a);j++)for(long i=1;i<lg(gel(a,j));i++){int c=cmpii(gcoeff(a,i,j),gcoeff(b,i,j));if(c)return c;}return 0;}
int main(void){static char data[4000000];pari_init(256000000,500000);
 while(scanf("%3999999s",data)==1){pari_sp av=avma;
  pari_CATCH(CATCH_ALL){char*s=pari_err2str(pari_err_last());error(s);pari_free(s);}
  pari_TRY{
   GEN v=gp_read_str(data),nf=nfinit0(RgV_to_RgX(gel(v,2),0),nf_NOLLL,DEFAULTPREC),p=gel(v,3);
   GEN primes=idealprimedec(nf,p),L=cgetg(lg(primes),t_VEC),out=cgetg(lg(primes),t_VEC);
   for(long i=1;i<lg(primes);i++)gel(L,i)=idealhnf(nf,gel(primes,i));
   L=gen_sort(L,NULL,compare_hnf);
   if(itos(gel(v,1))==0){
    GEN bases=cgetg(lg(L),t_VEC);for(long i=1;i<lg(L);i++){
     GEN H=gel(L,i),B=cgetg(lg(H),t_MAT);long k=1;
     for(long j=1;j<lg(H);j++)if(equali1(gcoeff(H,j,j)))gel(B,k++)=gel(H,j);
     setlg(B,k);gel(bases,i)=B;
    }
    GEN P=primedec_end(nf,bases,p,0);
    for(long i=1;i<lg(P);i++)gel(out,i)=mkvec2(gel(L,i),gel(P,i));
   }else{
    GEN generators=gel(v,4),I=zeromat(0,0);
    for(long i=1;i<lg(generators);i++)I=idealadd(nf,I,gdiv(RgV_to_RgX(gel(generators,i),0),gel(v,5)));
    for(long i=1;i<lg(L);i++){
     GEN P=NULL;for(long j=1;j<lg(primes);j++)if(gequal(idealhnf(nf,gel(primes,j)),gel(L,i))){P=gel(primes,j);break;}
     gel(out,i)=mkvec2(gel(L,i),itos(gel(v,1))==2?stoi(ZC_nfval(shallowtrans(gel(generators,1)),P)):gpidealval(nf,I,P));
    }
   }
   value(out);printf("\n");
  }pari_ENDCATCH;set_avma(av);fflush(stdout);
 }pari_close();return 0;
}
