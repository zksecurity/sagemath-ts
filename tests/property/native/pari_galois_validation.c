#include "../../../reference/pari/src/basemath/galconj.c"
#include <stdio.h>
static void print_error(const char *s){
 printf("ERROR ");putchar(34);
 for(;*s;s++){unsigned char c=*s;if(c==34||c==92){putchar(92);putchar(c);}else if(c<32)printf("\\u%04x",c);else putchar(c);}
 putchar(34);putchar(10);
}




static void print_gen(GEN v){long t=typ(v);if(t==t_INT){pari_printf("%Ps",v);return;}printf("[");for(long i=1;i<lg(v);i++){if(i>1)printf(",");if(t==t_VECSMALL)printf("%ld",v[i]);else print_gen(gel(v,i));}printf("]");}
static GEN matrix(GEN flat,long n){long k=(lg(flat)-1)/n;GEN M=cgetg(k+1,t_MAT);for(long i=1;i<=k;i++){GEN c=cgetg(n+1,t_COL);for(long j=1;j<=n;j++)gel(c,j)=gel(flat,(i-1)*n+j);gel(M,i)=c;}return M;}
static void poly(GEN z){GEN d,num=Q_remove_denom(z,&d);if(typ(num)!=t_POL)num=scalarpol(num,0);printf("[");for(long j=2;j<lg(num);j++){if(j>2)printf(",");pari_printf("%Ps",gel(num,j));}pari_printf("]/%Ps",d?d:gen_1);}
int main(void){long op,hasden;static char sa[2000000],sd[30000];pari_init(256000000,500000);
while(scanf("%ld %1999999s %ld %29999s",&op,sa,&hasden,sd)==4){pari_sp av=avma;
pari_CATCH(CATCH_ALL){char*err=pari_err2str(pari_err_last());print_error(err);pari_free(err);}
pari_TRY{GEN T=gtopolyrev(gp_read_str(sa),0),d=hasden?gp_read_str(sd):NULL,z;
if(op==0){z=galoisinit(T,d);printf("OK %d\n",typ(z)==t_VEC?1:0);}
else {z=galoisconj4(T,d);printf("OK [");for(long i=1;i<lg(z);i++){if(i>1)printf(",");poly(gel(z,i));}printf("]\n");}}
pari_ENDCATCH;set_avma(av);fflush(stdout);}pari_close();return 0;}
