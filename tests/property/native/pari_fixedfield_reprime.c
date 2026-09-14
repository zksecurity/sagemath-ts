#include "../../../reference/pari/src/basemath/galconj.c"
#include <stdio.h>
static void print_error(const char *s){
 printf("ERROR ");putchar(34);
 for(;*s;s++){unsigned char c=*s;if(c==34||c==92){putchar(92);putchar(c);}else if(c<32)printf("\\u%04x",c);else putchar(c);}
 putchar(34);putchar(10);
}




static void print_gen(GEN v){long t=typ(v);if(t==t_POL){printf("[");for(long i=2;i<lg(v);i++){if(i>2)printf(",");print_gen(gel(v,i));}printf("]");return;}if(t==t_INT){pari_printf("%Ps",v);return;}printf("[");for(long i=1;i<lg(v);i++){if(i>1)printf(",");if(t==t_VECSMALL)printf("%ld",v[i]);else print_gen(gel(v,i));}printf("]");}
static GEN matrix(GEN flat,long n){long k=(lg(flat)-1)/n;GEN M=cgetg(k+1,t_MAT);for(long i=1;i<=k;i++){GEN c=cgetg(n+1,t_COL);for(long j=1;j<=n;j++)gel(c,j)=gel(flat,(i-1)*n+j);gel(M,i)=c;}return M;}

int main(void){long e,o,start;static char st[30000],sl[30000],sL[200000],ss[30000],sd[30000],sp[200000],sr[30000];pari_init(256000000,500000);
while(scanf("%29999s %29999s %ld %199999s %29999s %29999s %199999s %ld %ld %29999s",st,sl,&e,sL,ss,sd,sp,&o,&start,sr)==10){pari_sp av=avma;
pari_CATCH(CATCH_ALL){char*err=pari_err2str(pari_err_last());print_error(err);pari_free(err);}
pari_TRY{GEN T=gtopolyrev(gp_read_str(st),0),l=gp_read_str(sl),L=gp_read_str(sL),S=gdiv(gtopolyrev(gp_read_str(ss),0),gp_read_str(sd)),O=perm_cycles(ZV_to_zv(gp_read_str(sp))),V=NULL;struct galois_borne gb={0};struct galois_frobenius gf={0};gb.l=l;gb.valabs=gb.valsol=e;gb.ladicabs=gb.ladicsol=powiu(l,e);gb.dis=ZX_disc(T);gb.bornesol=gen_0;gf.deg=o;gf.p=start;gf.Tmod=cgetg(1,t_VEC);gf.psi=cgetg(1,t_VECSMALL);setrand(gp_read_str(sr));GEN PG=galoisgenfixedfield0(O,L,S,T,NULL,&V,&gf,&gb);printf("OK [");if(PG)print_gen(mkvec2(gg_get_std(gel(PG,1)),gel(PG,2)));else printf("0");printf(",");if(V)print_gen(V);else printf("0");printf(",%ld,",gf.p);print_gen(gf.Tmod);printf(",");print_gen(gf.psi);printf("]|");print_gen(getrand());printf("\n");}
pari_ENDCATCH;set_avma(av);fflush(stdout);}pari_close();return 0;}
