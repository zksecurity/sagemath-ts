#include "pari.h"
#include "paripriv.h"
#include <stdio.h>
static void print_error(const char *s){
 printf("ERROR ");putchar(34);
 for(;*s;s++){unsigned char c=*s;if(c==34||c==92){putchar(92);putchar(c);}else if(c<32)printf("\\u%04x",c);else putchar(c);}
 putchar(34);putchar(10);
}

int main(void){long op;static char sp[30000],sa[30000],sb[30000],sc[30000];pari_init(256000000,500000);
while(scanf("%ld %29999s %29999s %29999s %29999s",&op,sp,sa,sb,sc)==5){pari_sp av=avma;
pari_CATCH(CATCH_ALL){char*e=pari_err2str(pari_err_last());print_error(e);pari_free(e);}
pari_TRY{GEN p=gp_read_str(sp),a=gp_read_str(sa),b=gp_read_str(sb),c=gp_read_str(sc),z;
switch(op){case 0:z=Fp_red(a,p);break;case 1:z=Fp_add(a,b,p);break;case 2:z=Fp_sub(a,b,p);break;
case 3:z=Fp_neg(a,p);break;case 4:z=Fp_mul(a,b,p);break;case 5:z=Fp_sqr(a,p);break;
case 6:z=Fp_inv(a,p);break;case 7:z=Fp_div(a,b,p);break;case 8:z=Fp_addmul(a,b,c,p);break;
case 9:z=Fp_double(a,p);break;case 10:z=Fp_halve(a,p);break;case 11:z=Fp_center(a,p,shifti(p,-1));break;
case 12:z=Fp_mulu(a,itou(b),p);break;default:z=stoi(equalii(a,b));}
pari_printf("OK %Ps\n",z);}
pari_ENDCATCH;set_avma(av);fflush(stdout);}pari_close();return 0;}
