#include "pari.h"
#include "paripriv.h"
#include <stdio.h>
static void poly(GEN f){printf("[");for(long i=2;i<lg(f);i++){if(i>2)printf(",");pari_printf("%Ps",gel(f,i));}printf("]");}
static void print_error(const char *s){
 printf("ERROR ");putchar(34);
 for(;*s;s++){unsigned char c=*s;if(c==34||c==92){putchar(92);putchar(c);}else if(c<32)printf("\\u%04x",c);else putchar(c);}
 putchar(34);putchar(10);
}

int main(void){long op;static char sp[30000],sx[30000],sf[2000000];pari_init(256000000,500000);
while(scanf("%ld %29999s %29999s %1999999s",&op,sp,sx,sf)==4){pari_sp av=avma;
pari_CATCH(CATCH_ALL){char*e=pari_err2str(pari_err_last());print_error(e);pari_free(e);}
pari_TRY{GEN p=gp_read_str(sp),x=gp_read_str(sx),f=gtopolyrev(gp_read_str(sf),0),z;
if(op==0)z=FpX_deriv(f,p);else if(op==1)z=FpX_eval(f,x,p);else if(op==2)z=FpX_center(f,p,x);else if(op==3)z=FpX_div_by_X_x(f,x,p,NULL);else z=stoi(FpX_is_squarefree(f,p));
printf("OK ");if(op==1||op==4)pari_printf("%Ps",z);else poly(z);printf("\n");}
pari_ENDCATCH;set_avma(av);fflush(stdout);}pari_close();return 0;}
