/* Bundled PARI ideal intersection and independent ideal-HNF fixture construction. */
#include "pari.h"
#include "paripriv.h"
#include <stdio.h>
#ifndef LONG_IS_64BIT
#error "Ideal intersection oracle requires 64-bit PARI"
#endif
static void value(GEN x) {
  if (!x) { printf("null"); return; }
  if (typ(x)==t_INT) { pari_printf("\"%Ps\"",x); return; }
  if (typ(x)==t_REAL) {
    GEN m=gen_0;
    if(signe(x))for(long i=2;i<lg(x);i++)m=addii(shifti(m,64),utoi((ulong)x[i]));
    printf("[%ld,\"%ld\",",signe(x),expo(x));pari_printf("\"%Ps\",%ld]",m,(lg(x)-2)*64);return;
  }
  printf("[");for(long i=1;i<lg(x);i++){if(i>1)printf(",");value(gel(x,i));}printf("]");
}
static void error(const char *s){
 printf("ERROR ");putchar(34);
 for(;*s;s++){unsigned char c=*s;if(c==34||c==92){putchar(92);putchar(c);}else if(c<32)printf("\\u%04x",c);else putchar(c);}
 putchar(34);putchar(10);
}

int main(void){
 long op,n;static char data[4000000];GEN fields[65]={0};pari_init(256000000,500000);
 while(scanf("%ld %ld %3999999s",&op,&n,data)==3){
  pari_sp av=avma;
  pari_CATCH(CATCH_ALL){char *s=pari_err2str(pari_err_last());error(s);pari_free(s);}
  pari_TRY{
   if(!fields[n])fields[n]=gclone(nfinit0(gsubgs(pol_xn(n,0),2),0,DEFAULTPREC));
   GEN nf=fields[n],v=gp_read_str(data),I,J,z,den;
   if(op){
    I=idealhnf0(nf,powiu(gel(v,1),itos(gel(v,3))),gpowgs(gsub(pol_x(0),gel(v,2)),itos(gel(v,3))));
    J=idealhnf0(nf,powiu(gel(v,4),itos(gel(v,6))),gpowgs(gsub(pol_x(0),gel(v,5)),itos(gel(v,6))));
    value(mkvec2(I,J));
   }else{
    long ni=itos(gel(v,4)),nj=itos(gel(v,5));I=cgetg(ni+1,t_MAT);J=cgetg(nj+1,t_MAT);
    for(long j=1;j<=ni;j++){gel(I,j)=cgetg(n+1,t_COL);for(long i=1;i<=n;i++)gcoeff(I,i,j)=gdiv(gel(v,5+(i-1)*ni+j),gel(v,2));}
    for(long j=1;j<=nj;j++){gel(J,j)=cgetg(n+1,t_COL);for(long i=1;i<=n;i++)gcoeff(J,i,j)=gdiv(gel(v,5+n*ni+(i-1)*nj+j),gel(v,3));}
    setrand(gel(v,1));z=Q_remove_denom(idealintersect(nf,I,J),&den);
    value(mkvec3(z,den?den:gen_1,randomi(int2n(128))));
   }
   printf("\n");
  }pari_ENDCATCH;set_avma(av);fflush(stdout);
 }
 for(long i=0;i<65;i++)if(fields[i])gunclone(fields[i]);pari_close();return 0;
}
