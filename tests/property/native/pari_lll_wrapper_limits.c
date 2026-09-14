#include "lll.c"
#include <stdio.h>
int main(int argc,char **argv){pari_init(64000000,500000);GEN x=gp_read_str(argv[2]),r,N;int op=atoi(argv[1]);pari_CATCH(CATCH_ALL){pari_printf("ERROR %Ps\n",pari_err_last());}pari_TRY{if(op==0)r=ZM_lll_norms(x,.99,LLL_IM,&N);else if(op==1)r=ZM_flatter_rank(x,1,LLL_INPLACE);else if(op==2)r=ZM_flattergram_rank(x,1,0);else r=ZM_flattergram(x,0);pari_printf("%Ps\n",r);}pari_ENDCATCH;return 0;}
