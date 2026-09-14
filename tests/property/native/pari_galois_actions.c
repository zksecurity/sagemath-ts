#include "../../../reference/pari/src/basemath/galconj.c"
#include <stdio.h>
static void print_error(const char *s){
 printf("ERROR ");putchar(34);
 for(;*s;s++){unsigned char c=*s;if(c==34||c==92){putchar(92);putchar(c);}else if(c<32)printf("\\u%04x",c);else putchar(c);}
 putchar(34);putchar(10);
}




static void print_gen(GEN v){long t=typ(v);if(t==t_POL){printf("[");for(long i=2;i<lg(v);i++){if(i>2)printf(",");print_gen(gel(v,i));}printf("]");return;}if(t==t_INT){pari_printf("%Ps",v);return;}printf("[");for(long i=1;i<lg(v);i++){if(i>1)printf(",");if(t==t_VECSMALL)printf("%ld",v[i]);else print_gen(gel(v,i));}printf("]");}
static GEN matrix(GEN flat,long n){long k=(lg(flat)-1)/n;GEN M=cgetg(k+1,t_MAT);for(long i=1;i<=k;i++){GEN c=cgetg(n+1,t_COL);for(long j=1;j<=n;j++)gel(c,j)=gel(flat,(i-1)*n+j);gel(M,i)=c;}return M;}

static GEN perms(GEN flat,long n){GEN v=cgetg((lg(flat)-1)/n+1,t_VEC);for(long j=1;j<lg(v);j++){GEN p=cgetg(n+1,t_VECSMALL);for(long i=1;i<=n;i++)p[i]=itos(gel(flat,(j-1)*n+i));gel(v,j)=p;}return v;}
static void qpoly(GEN z){GEN d,num=Q_remove_denom(z,&d);if(typ(num)!=t_POL)num=scalarpol(num,0);print_gen(num);pari_printf("/%Ps",d?d:gen_1);}
int main(void){long op,e,kind,flag;static char st[30000],sp[30000],sl[200000],sm[2000000],sd[30000],sg[200000],ss[200000],so[30000],sa[200000],sb[30000];pari_init(256000000,500000);
while(scanf("%ld %29999s %29999s %ld %199999s %1999999s %29999s %199999s %199999s %29999s %ld %199999s %29999s %ld",&op,st,sp,&e,sl,sm,sd,sg,ss,so,&kind,sa,sb,&flag)==14){pari_sp av=avma;
pari_CATCH(CATCH_ALL){char*err=pari_err2str(pari_err_last());print_error(err);pari_free(err);}
pari_TRY{GEN T=gtopolyrev(gp_read_str(st),0),p=gp_read_str(sp),gal=cgetg(9,t_VEC),a=gp_read_str(sa),b=gp_read_str(sb),P,z;long n=degpol(T);
gel(gal,1)=T;gel(gal,2)=mkvec3(p,stoi(e),powiu(p,e));gel(gal,3)=gp_read_str(sl);gel(gal,4)=matrix(gp_read_str(sm),n);gel(gal,5)=gp_read_str(sd);gel(gal,6)=perms(gp_read_str(sg),n);gel(gal,7)=perms(gp_read_str(ss),n);gel(gal,8)=ZV_to_zv(gp_read_str(so));
P=kind==0?ZV_to_zv(a):kind==1?perms(a,n):mkvec2(perms(a,n),ZV_to_zv(b));
if(op<=1)z=galoispermtopol(gal,P);else if(op==2)z=galoisfixedfield(gal,P,flag,1);else if(op==3)z=galois_group(gal);else z=galoissubgroups(op==4?gal:P);
printf("OK ");if(op==0)qpoly(z);else if(op==1){printf("[");for(long i=1;i<lg(z);i++){if(i>1)printf(",");qpoly(gel(z,i));}printf("]");}
else if(op==2){if(flag==1)print_gen(z);else{printf("[");print_gen(gel(z,1));printf(",");qpoly(lift(gel(z,2)));if(flag==2){GEN F=gel(z,3);printf(",[");for(long i=1;i<lg(F);i++){GEN f=gel(F,i);if(i>1)printf(",");printf("[");for(long j=2;j<lg(f);j++){if(j>2)printf(",");qpoly(gel(f,j));}printf("]");}printf("]");}printf("]");}}
else print_gen(z);printf("\n");}
pari_ENDCATCH;set_avma(av);fflush(stdout);}pari_close();return 0;}
