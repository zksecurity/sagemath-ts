#include "../../../reference/pari/src/basemath/galconj.c"
#include <stdio.h>
static void print_error(const char *s){
 printf("ERROR ");putchar(34);
 for(;*s;s++){unsigned char c=*s;if(c==34||c==92){putchar(92);putchar(c);}else if(c<32)printf("\\u%04x",c);else putchar(c);}
 putchar(34);putchar(10);
}


static void print_poly(GEN f){printf("[");for(long i=2;i<lg(f);i++){if(i>2)printf(",");pari_printf("%Ps",gel(f,i));}printf("]");}
int main(void){long op;static char sp[30000],sh[30000],sd[30000],ss[30000],sl[2000000],sm[2000000];pari_init(256000000,500000);
while(scanf("%ld %29999s %29999s %29999s %29999s %1999999s %1999999s",&op,sp,sh,sd,ss,sl,sm)==7){pari_sp av=avma;
pari_CATCH(CATCH_ALL){char*e=pari_err2str(pari_err_last());print_error(e);pari_free(e);}
pari_TRY{GEN p=gp_read_str(sp),half=gp_read_str(sh),den=gp_read_str(sd),perm=gtovecsmall(gp_read_str(ss)),L=gp_read_str(sl),v=gp_read_str(sm),z;
if(op==3){long n=lg(L)-1;GEN M=cgetg(n+1,t_MAT);for(long j=1;j<=n;j++){gel(M,j)=cgetg(n+1,t_COL);for(long i=1;i<=n;i++)gcoeff(M,i,j)=gel(v,(i-1)*n+j);}z=permtopol(perm,L,M,den,p,half,0);}
else {z=gdiv(gtopolyrev(L,0),den);if(op==2)z=RgX_to_FpX(z,p);}
printf("OK ");if(op==0){GEN d,num=Q_remove_denom(z,&d);print_poly(num);pari_printf("/%Ps",d?d:gen_1);}else print_poly(z);printf("\n");}
pari_ENDCATCH;set_avma(av);fflush(stdout);}pari_close();return 0;}
