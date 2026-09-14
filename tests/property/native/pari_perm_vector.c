#include "pari.h"
#include "paripriv.h"
#include <stdio.h>
static void print_error(const char *s){
 printf("ERROR ");putchar(34);
 for(;*s;s++){unsigned char c=*s;if(c==34||c==92){putchar(92);putchar(c);}else if(c<32)printf("\\u%04x",c);else putchar(c);}
 putchar(34);putchar(10);
}



static void vec(GEN v){printf("[");for(long i=1;i<lg(v);i++){if(i>1)printf(",");if(typ(v)==t_VEC)vec(gel(v,i));else printf("%.0f",(double)v[i]);}printf("]");}
int main(void){long op;static char sa[2000000],sb[2000000],se[30000];pari_init(256000000,500000);
while(scanf("%ld %1999999s %1999999s %29999s",&op,sa,sb,se)==4){pari_sp av=avma;
pari_CATCH(CATCH_ALL){char*err=pari_err2str(pari_err_last());print_error(err);pari_free(err);}
pari_TRY{GEN a=ZV_to_zv(gp_read_str(sa)),b=ZV_to_zv(gp_read_str(sb)),E=gp_read_str(se),z=NULL;long c=0;ulong u=0;int scalar=0;
/* zv_prod explicitly assumes no signed overflow: certify each prefix
 * with native arbitrary-precision arithmetic before calling it. */
if(op==12){GEN product=gen_1,max=shifti(gen_1,63),min=negi(max);for(long i=1;i<lg(a);i++){product=mulis(product,a[i]);if(cmpii(product,min)<0||cmpii(product,max)>=0){printf("GUARD zv_prod requires a product within the signed-word range\n");goto done;}}}
switch(op){case 0:z=perm_mul(a,b);break;case 1:z=perm_sqr(a);break;case 2:z=perm_inv(a);break;case 3:z=perm_conj(a,b);break;case 4:c=perm_commute(a,b);scalar=1;break;case 5:z=perm_powu(a,itou(E));break;case 6:z=perm_cycles(a);break;case 7:u=perm_orderu(a);scalar=2;break;case 8:z=cyc_pow(perm_cycles(a),itos(E));break;case 9:c=vecsmall_lexcmp(a,b);scalar=1;break;case 10:c=zv_equal(a,b);scalar=1;break;case 11:z=vecsmall_uniq(a);break;case 12:c=group_order(mkvec2(cgetg(1,t_VEC),a));scalar=1;break;case 13:c=group_domain(mkvec2(signe(E)?mkvec(a):cgetg(1,t_VEC),cgetg(1,t_VECSMALL)));scalar=1;break;case 14:z=vecperm_orbits(itos(E)==0?cgetg(1,t_VEC):itos(E)==1?mkvec(a):mkvec2(a,b),lg(a)-1);break;}
printf("OK ");if(scalar==1)printf("%.0f",(double)c);else if(scalar==2)printf("%.0f",(double)u);else vec(z);printf("\n");done:;}
pari_ENDCATCH;set_avma(av);fflush(stdout);}pari_close();return 0;}
