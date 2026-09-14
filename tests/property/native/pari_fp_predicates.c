#include "pari.h"
#include "paripriv.h"
#include <stdio.h>
static void print_error(const char *s){
 printf("ERROR ");putchar(34);
 for(;*s;s++){unsigned char c=*s;if(c==34||c==92){putchar(92);putchar(c);}else if(c<32)printf("\\u%04x",c);else putchar(c);}
 putchar(34);putchar(10);
}

int main(void){long op;static char sp[30000],sa[30000],se[30000];pari_init(256000000,500000);
while(scanf("%ld %29999s %29999s %29999s",&op,sp,sa,se)==4){pari_sp av=avma;
pari_CATCH(CATCH_ALL){char*e=pari_err2str(pari_err_last());print_error(e);pari_free(e);}
pari_TRY{GEN p=gp_read_str(sp),a=gp_read_str(sa),e=gp_read_str(se),z; if(op==0)z=stoi(Fp_issquare(a,p));else if(op==1)z=stoi(kronecker(a,p));else if(op==2)z=gcdii(a,p);else if(op==3){GEN u,v,d=bezout(a,p,&u,&v);z=mkvec3(d,u,v);}else if(op==4)z=Fp_order(a,e,p);else z=znorder(gmodulo(a,p),op==5?NULL:e);pari_printf("OK %Ps\n",z);}
pari_ENDCATCH;set_avma(av);fflush(stdout);}pari_close();return 0;}
