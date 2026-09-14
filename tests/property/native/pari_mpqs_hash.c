#include "../../../reference/pari/src/basemath/mpqs.c"
#include <stdio.h>
static void print_error(const char *s){
 printf("ERROR ");putchar(34);
 for(;*s;s++){unsigned char c=*s;if(c==34||c==92){putchar(92);putchar(c);}else if(c<32)printf("\\u%04x",c);else putchar(c);}
 putchar(34);putchar(10);
}




static void print_gen(GEN v){long t=typ(v);if(t==t_POL){printf("[");for(long i=2;i<lg(v);i++){if(i>2)printf(",");print_gen(gel(v,i));}printf("]");return;}if(t==t_INT){pari_printf("%Ps",v);return;}printf("[");for(long i=1;i<lg(v);i++){if(i>1)printf(",");if(t==t_VECSMALL)printf("%ld",v[i]);else print_gen(gel(v,i));}printf("]");}
static GEN matrix(GEN flat,long n){long k=(lg(flat)-1)/n;GEN M=cgetg(k+1,t_MAT);for(long i=1;i<=k;i++){GEN c=cgetg(n+1,t_COL);for(long j=1;j<=n;j++)gel(c,j)=gel(flat,(i-1)*n+j);gel(M,i)=c;}return M;}




int main(){char*line=NULL;size_t cap=0;pari_init(128000000,500000);while(getline(&line,&cap,stdin)>0){pari_sp av=avma;pari_CATCH(CATCH_ALL){char*e=pari_err2str(pari_err_last());print_error(e);pari_free(e);}pari_TRY{GEN a=gp_read_str(line),packed=gel(a,2);hashtable h;hash_init_GEN(&h,itou(gel(a,1)),gequal,1);printf("OK [[");long first=1;for(long i=1;i<lg(packed);){GEN Y=gel(packed,i++);long n=itos(gel(packed,i++));GEN r=cgetg(n+1,t_VECSMALL);for(long j=1;j<=n;j++)r[j]=itos(gel(packed,i++));GEN R=mkvec2(Y,r);if(!first)printf(",");first=0;printf("%lu",hash_GEN(R));frel_add(&h,R);}printf("],");print_gen(hash_keys_GEN(&h));printf(",%lu]\n",h.nb);}pari_ENDCATCH;set_avma(av);fflush(stdout);}pari_close();}
