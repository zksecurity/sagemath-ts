/* Original finite-field image, basis completion, inverse and prime HNF routines. */
#include "pari.h"
#include "paripriv.h"
#include <stdio.h>
static GEN matrix_from(GEN v,long m,long n){GEN A=cgetg(n+1,t_MAT);for(long j=1;j<=n;j++){gel(A,j)=cgetg(m+1,t_COL);for(long i=1;i<=m;i++)gcoeff(A,i,j)=gel(v,(i-1)*n+j);}return A;}
static void value(GEN x){if(!x){printf("null");return;}if(typ(x)==t_INT){pari_printf("\"%Ps\"",x);return;}printf("[");for(long i=1;i<lg(x);i++){if(i>1)printf(",");if(typ(x)==t_VECSMALL)printf("\"%ld\"",x[i]);else value(gel(x,i));}printf("]");}
static void error(const char*s){printf("ERROR ");putchar(34);for(;*s;s++){unsigned char c=*s;if(c==34||c==92){putchar(92);putchar(c);}else if(c<32)printf("\\u%04x",c);else putchar(c);}putchar(34);putchar(10);}
int main(void){static char data[4000000],prime[20000];long op,m,n;pari_init(256000000,500000);
while(scanf("%ld %ld %ld %19999s %3999999s",&op,&m,&n,prime,data)==5){pari_sp av=avma;
pari_CATCH(CATCH_ALL){char*s=pari_err2str(pari_err_last());error(s);pari_free(s);}
pari_TRY{GEN p=gp_read_str(prime),v=gp_read_str(data),A=matrix_from(v,m,n),out;
 if(op==0)out=FpM_inv(A,p);else if(op==1)out=FpM_image(A,p);else if(op==2)out=FpM_suppl(A,p);else if(op==3)out=ZM_hnfmodprime(A,p);
 else if(op==4){long r;GEN d=F2m_gauss_pivot(ZM_to_F2m(A),&r);out=mkvec2(d?zv_to_ZV(d):gen_0,stoi(r));}
 else {GEN B=matid(m);out=F2m_gauss(ZM_to_F2m(A),ZM_to_F2m(B));if(out)out=F2m_to_ZM(out);}
 value(out);printf("\n");}pari_ENDCATCH;set_avma(av);fflush(stdout);}pari_close();return 0;}
