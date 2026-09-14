#include "pari.h"
#include "paripriv.h"
#include <stdio.h>
static void print_error(const char *s){
 printf("ERROR ");putchar(34);
 for(;*s;s++){unsigned char c=*s;if(c==34||c==92){putchar(92);putchar(c);}else if(c<32)printf("\\u%04x",c);else putchar(c);}
 putchar(34);putchar(10);
}



int main(void){long op;static char sl[2000000],sd[30000],sp[30000];pari_init(256000000,500000);
while(scanf("%ld %1999999s %29999s %29999s",&op,sl,sd,sp)==4){pari_sp av=avma;
pari_CATCH(CATCH_ALL){char*err=pari_err2str(pari_err_last());print_error(err);pari_free(err);}
pari_TRY{GEN L=gp_read_str(sl),den=gp_read_str(sd),p=gp_read_str(sp),z=op==0?FpV_invVandermonde(L,den,p):op==1?FpV_inv(L,p):producttree_scheme(itos(den));printf("OK [");if(op==0){for(long i=1;i<lg(L);i++){if(i>1)printf(",");printf("[");for(long j=1;j<lg(z);j++){if(j>1)printf(",");pari_printf("%Ps",gcoeff(z,i,j));}printf("]");}}else for(long i=1;i<lg(z);i++){if(i>1)printf(",");if(op==2)printf("%ld",z[i]);else pari_printf("%Ps",gel(z,i));}printf("]\n");}
pari_ENDCATCH;set_avma(av);fflush(stdout);}pari_close();return 0;}
