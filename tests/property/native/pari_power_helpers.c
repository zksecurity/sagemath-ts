#include "../../../reference/pari/src/basemath/ifactor1.c"
#include <stdio.h>
static void print_error(const char *s){
 printf("ERROR ");putchar(34);
 for(;*s;s++){unsigned char c=*s;if(c==34||c==92){putchar(92);putchar(c);}else if(c<32)printf("\\u%04x",c);else putchar(c);}
 putchar(34);putchar(10);
}




static void print_gen(GEN v){long t=typ(v);if(t==t_POL){printf("[");for(long i=2;i<lg(v);i++){if(i>2)printf(",");print_gen(gel(v,i));}printf("]");return;}if(t==t_INT){pari_printf("%Ps",v);return;}printf("[");for(long i=1;i<lg(v);i++){if(i>1)printf(",");if(t==t_VECSMALL)printf("%ld",v[i]);else print_gen(gel(v,i));}printf("]");}
static GEN matrix(GEN flat,long n){long k=(lg(flat)-1)/n;GEN M=cgetg(k+1,t_MAT);for(long i=1;i<=k;i++){GEN c=cgetg(n+1,t_COL);for(long j=1;j<=n;j++)gel(c,j)=gel(flat,(i-1)*n+j);gel(M,i)=c;}return M;}


int main(void){long op,a,b,c;static char sn[30000];pari_init(256000000,500000);
while(scanf("%ld %29999s %ld %ld %ld",&op,sn,&a,&b,&c)==5){pari_sp av=avma;
pari_CATCH(CATCH_ALL){char*err=pari_err2str(pari_err_last());print_error(err);pari_free(err);}
pari_TRY{GEN n=gp_read_str(sn),root=gen_0,z;long e;
switch(op){
case 0:{ulong mask=a;e=is_357_power(n,&root,&mask);z=mkvec3(stoi(e),e?root:gen_0,utoi(mask));break;}
case 1:e=is_kth_power(n,a,&root);z=e?root:gen_0;break;
case 2:e=Z_isanypower(n,&root);z=mkvec2(stoi(e),e?root:n);break;
case 3:{forprime_t T;u_forprime_init(&T,a,b);e=is_pth_power(n,&root,&T,c);z=mkvec3(stoi(e),e?root:n,utoi(T.p>(ulong)b?0:T.p));break;}
case 4:z=utoi(tridiv_bound(n));break;
case 5:e=Z_issquareall(n,&root);z=e?mkvec(root):cgetg(1,t_VEC);break;
case 6:e=isprimepower(n,&root);z=e?mkvec2(root,stoi(e)):cgetg(1,t_VEC);break;
default:z=gen_0;break;}
printf("OK ");print_gen(z);printf("\n");}
pari_ENDCATCH;set_avma(av);fflush(stdout);}pari_close();return 0;}
