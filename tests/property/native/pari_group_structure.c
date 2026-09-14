#include "pari.h"
#include "paripriv.h"
#include <stdio.h>
static void print_error(const char *s){
 printf("ERROR ");putchar(34);
 for(;*s;s++){unsigned char c=*s;if(c==34||c==92){putchar(92);putchar(c);}else if(c<32)printf("\\u%04x",c);else putchar(c);}
 putchar(34);putchar(10);
}



static GEN group(GEN flat,GEN ord,long n){GEN gen=cgetg(lg(ord),t_VEC);for(long j=1;j<lg(ord);j++){GEN v=cgetg(n+1,t_VECSMALL);for(long i=1;i<=n;i++)v[i]=itos(gel(flat,(j-1)*n+i));gel(gen,j)=v;}return mkvec2(gen,ZV_to_zv(ord));}
static void vec(GEN v){printf("[");for(long i=1;i<lg(v);i++){if(i>1)printf(",");if(typ(v)==t_VEC)vec(gel(v,i));else printf("%ld",v[i]);}printf("]");}
int main(void){long op,n;static char sa[2000000],sb[30000],sc[2000000],sd[30000],sp[2000000];pari_init(256000000,500000);
while(scanf("%ld %ld %1999999s %29999s %1999999s %29999s %1999999s",&op,&n,sa,sb,sc,sd,sp)==7){pari_sp av=avma;
pari_CATCH(CATCH_ALL){char*err=pari_err2str(pari_err_last());print_error(err);pari_free(err);}
pari_TRY{GEN G=group(gp_read_str(sa),gp_read_str(sb),n),H=group(gp_read_str(sc),gp_read_str(sd),n),p=ZV_to_zv(gp_read_str(sp)),z=NULL,C;long scalar=-1;
switch(op){case 0:z=group_elts(G,n);break;case 1:z=group_leftcoset(G,p);break;case 2:z=group_rightcoset(G,p);break;case 3:z=group_set(G,n);break;case 4:z=groupelts_set(group_elts(G,n),n);break;case 5:scalar=group_perm_normalize(G,p);break;case 6:z=group_subgroups(G);break;case 7:scalar=group_isA4S4(G);break;case 8:z=group_quotient(G,H);break;case 9:z=groupelts_quotient(group_elts(G,n),H);break;case 10:C=groupelts_quotient(group_elts(G,n),H);z=quotient_perm(C,p);break;case 11:C=groupelts_quotient(group_elts(G,n),H);z=quotient_group(C,G);break;case 12:C=groupelts_quotient(group_elts(G,n),H);z=quotient_subgroup_lift(C,H,trivialgroup());break;case 13:C=groupelts_quotient(group_elts(G,n),H);z=quotient_subgroup_lift(C,H,quotient_group(C,G));break;case 14:scalar=perm_relorder(p,group_set(H,n));break;case 15:{GEN el=group_elts(H,n);z=perm_generate(p,el,perm_relorder(p,groupelts_set(el,n)));break;}case 16:z=cyclicgroup(p,perm_orderu(p));break;case 17:z=dicyclicgroup(gmael(G,1,1),gmael(G,1,2),mael(G,2,1),mael(G,2,2));break;case 18:z=trivialgroup();break;}
printf("OK ");if(scalar>=0)printf("%ld",scalar);else if(op==3||op==4){printf("[");for(long i=1;i<=n;i++){if(i>1)printf(",");printf("%ld",F2v_coeff(z,i));}printf("]");}else vec(z);printf("\n");}
pari_ENDCATCH;set_avma(av);fflush(stdout);}pari_close();return 0;}
