/* Original PARI integer rounding, HNF centering, nfbasis and rational idealmul. */
#include "pari.h"
#include "paripriv.h"
#include <stdio.h>
static void value(GEN x){
 if(typ(x)==t_INT||typ(x)==t_FRAC){pari_printf("\"%Ps\"",x);return;}
 printf("[");for(long i=1;i<lg(x);i++){if(i>1)printf(",");value(gel(x,i));}printf("]");
}
static void error(const char*s){printf("ERROR ");putchar(34);for(;*s;s++){unsigned char c=*s;if(c==34||c==92){putchar(92);putchar(c);}else if(c<32)printf("\\u%04x",c);else putchar(c);}putchar(34);putchar(10);}
static GEN matrix_from(GEN v,long start,long n,long m,GEN den){
 GEN M=cgetg(m+1,t_MAT);for(long j=1;j<=m;j++){gel(M,j)=cgetg(n+1,t_COL);for(long i=1;i<=n;i++)gcoeff(M,i,j)=gdiv(gel(v,start+(i-1)*m+j),den);}return M;
}
int main(void){
 long op,n;static char data[4000000];GEN fields[65]={0};pari_init(256000000,500000);
 while(scanf("%ld %ld %3999999s",&op,&n,data)==3){pari_sp av=avma;
  pari_CATCH(CATCH_ALL){char*s=pari_err2str(pari_err_last());error(s);pari_free(s);}
  pari_TRY{
   GEN v=gp_read_str(data),result;
   if(op==0)result=ZM_hnfcenter(matrix_from(v,0,n,n,gen_1));
   else if(op==1)result=diviiround(gel(v,1),gel(v,2));
   else if(op==3){
    GEN f=RgV_to_RgX(v,0),basis=nfbasis(f,NULL);long degree=degpol(f);result=cgetg(degree+1,t_VEC);
    for(long j=1;j<=degree;j++){GEN c=cgetg(degree+1,t_VEC);gel(result,j)=c;for(long i=0;i<degree;i++)gel(c,i+1)=polcoef(gel(basis,j),i,-1);}
   }else{
    if(!fields[n])fields[n]=gclone(nfinit0(gsubgs(pol_xn(n,0),2),0,DEFAULTPREC));
    long ni=itos(gel(v,3)),nj=itos(gel(v,4));GEN I=matrix_from(v,4,ni,ni,gel(v,1)),J=matrix_from(v,4+ni*ni,nj,nj,gel(v,2)),den;
    setrand(utoi(42));GEN H=Q_remove_denom(idealmul(fields[n],I,J),&den);result=mkvec3(H,den?den:gen_1,randomi(int2n(128)));
   }
   value(result);printf("\n");
  }pari_ENDCATCH;set_avma(av);fflush(stdout);
 }
 for(long i=0;i<65;i++)if(fields[i])gunclone(fields[i]);pari_close();return 0;
}
