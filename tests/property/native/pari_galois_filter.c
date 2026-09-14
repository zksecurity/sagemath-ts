#include "../../../reference/pari/src/basemath/galconj.c"
#include <stdio.h>
static void print_error(const char *s){
 printf("ERROR ");putchar(34);
 for(;*s;s++){unsigned char c=*s;if(c==34||c==92){putchar(92);putchar(c);}else if(c<32)printf("\\u%04x",c);else putchar(c);}
 putchar(34);putchar(10);
}




static void print_gen(GEN v){long t=typ(v);if(t==t_POL){printf("[");for(long i=2;i<lg(v);i++){if(i>2)printf(",");print_gen(gel(v,i));}printf("]");return;}if(t==t_INT){pari_printf("%Ps",v);return;}printf("[");for(long i=1;i<lg(v);i++){if(i>1)printf(",");if(t==t_VECSMALL)printf("%ld",v[i]);else print_gen(gel(v,i));}printf("]");}
static GEN matrix(GEN flat,long n){long k=(lg(flat)-1)/n;GEN M=cgetg(k+1,t_MAT);for(long i=1;i<=k;i++){GEN c=cgetg(n+1,t_COL);for(long j=1;j<=n;j++)gel(c,j)=gel(flat,(i-1)*n+j);gel(M,i)=c;}return M;}

int main(void){long n;static char sl[200000],sm[2000000],sb[30000],sq[30000],sp[200000];pari_init(256000000,500000);
while(scanf("%ld %199999s %1999999s %29999s %29999s %199999s",&n,sl,sm,sb,sq,sp)==6){pari_sp av=avma;
pari_CATCH(CATCH_ALL){char*err=pari_err2str(pari_err_last());print_error(err);pari_free(err);}
pari_TRY{struct galois_test td;GEN L=gp_read_str(sl),M=matrix(gp_read_str(sm),n),A=gp_read_str(sp);inittest(L,M,gp_read_str(sb),gp_read_str(sq),&td);printf("OK [");for(long i=1;i<lg(A);i+=n){GEN p=cgetg(n+1,t_VECSMALL);for(long j=1;j<=n;j++)p[j]=itos(gel(A,i+j-1));long ok=galois_test_perm(&td,p);if(i>1)printf(",");printf("[%ld,",ok);print_gen(td.order);printf(",[");for(long k=1;k<=n;k++){if(k>1)printf(",");if(td.PV[k])print_gen(gel(td.PV,k));else printf("0");}printf("]]");}printf("]\n");freetest(&td);}
pari_ENDCATCH;set_avma(av);fflush(stdout);}pari_close();return 0;}
