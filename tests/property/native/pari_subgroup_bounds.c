#include "../../../reference/pari/src/basemath/galconj.c"
#include <stdio.h>
static void print_error(const char *s){
 printf("ERROR ");putchar(34);
 for(;*s;s++){unsigned char c=*s;if(c==34||c==92){putchar(92);putchar(c);}else if(c<32)printf("\\u%04x",c);else putchar(c);}
 putchar(34);putchar(10);
}




static void print_gen(GEN v){long t=typ(v);if(t==t_INT){pari_printf("%Ps",v);return;}printf("[");for(long i=1;i<lg(v);i++){if(i>1)printf(",");if(t==t_VECSMALL)printf("%ld",v[i]);else print_gen(gel(v,i));}printf("]");}
static GEN matrix(GEN flat,long n){long k=(lg(flat)-1)/n;GEN M=cgetg(k+1,t_MAT);for(long i=1;i<=k;i++){GEN c=cgetg(n+1,t_COL);for(long j=1;j<=n;j++)gel(c,j)=gel(flat,(i-1)*n+j);gel(M,i)=c;}return M;}
int main(void){long op,n;static char sa[2000000],sb[30000],sp[30000];pari_init(256000000,500000);
while(scanf("%ld %ld %1999999s %29999s %29999s",&op,&n,sa,sb,sp)==5){pari_sp av=avma;
pari_CATCH(CATCH_ALL){char*err=pari_err2str(pari_err_last());print_error(err);pari_free(err);}
pari_TRY{GEN a=gp_read_str(sa),b=gp_read_str(sb),p=gp_read_str(sp),z;
if(op==0){z=znstar(p);GEN v=gel(z,3);for(long i=1;i<lg(v);i++)gel(v,i)=lift(gel(v,i));}
else if(op==1)z=subgrouplist(a,mkvec(p));
else if(op==2)z=ZM_hnfmodid(matrix(a,n),b);
else if(op==3)z=listznstarelts(n,itos(p));
else if(op==4)z=pgener_Zp(p);
else if(op==5)z=subgrouplist(a,NULL);
else if(op==6)z=subgrouplist(a,p);
else if(op==8)z=subgrouplist(a,b);
else {GEN U,F=absZ_factor_limit_strict(p,itou(gel(a,1)),&U);z=mkvec2(F,U?U:gen_0);}
printf("OK ");print_gen(z);printf("\n");}
pari_ENDCATCH;set_avma(av);fflush(stdout);}pari_close();return 0;}
