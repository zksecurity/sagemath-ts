/* Direct ABI adapter for the bundled FLINT. Algorithm implementations are linked unchanged. */
#include "flint.h"
#include "ulong_extras.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "nmod_poly_factor.h"
#include "nmod_mat.h"

static int emit_word(ulong *out, slong cap, slong *len, ulong value) {
    if (*len >= cap) return 0;
    out[(*len)++] = value; return 1;
}
static int emit_vec(ulong *out, slong cap, slong *len, const ulong *a, slong n, ulong exponent) {
    if (!emit_word(out,cap,len,exponent) || !emit_word(out,cap,len,n)) return 0;
    for(slong i=0;i<n;i++) if(!emit_word(out,cap,len,a[i]))return 0;
    return 1;
}

slong audit_flint_factor(int op, ulong p, slong d,
 const ulong *a, slong an, const ulong *b, slong bn, const ulong *c, slong cn,
 int seeded, ulong seed1, ulong seed2, ulong *out, slong cap) {
    nmod_poly_t A,B,C,R,V;
    nmod_poly_factor_t factors;
    flint_rand_t state;
    slong len=0, chunks=1;
    ulong scalar=0;
    nmod_poly_init(A,p); nmod_poly_init(B,p); nmod_poly_init(C,p); nmod_poly_init(R,p);nmod_poly_init(V,p);
    nmod_poly_factor_init(factors);flint_rand_init(state);
    if(seeded)flint_rand_set_seed(state,seed1,seed2);
    for(slong i=0;i<an;i++)nmod_poly_set_coeff_ui(A,i,a[i]);
    for(slong i=0;i<bn;i++)nmod_poly_set_coeff_ui(B,i,b[i]);
    for(slong i=0;i<cn;i++)nmod_poly_set_coeff_ui(C,i,c[i]);
    slong *degrees=flint_malloc(sizeof(slong)*(an+2));
    ulong *vector=NULL;
    nmod_poly_struct *polys=NULL,*results=NULL;
    nmod_mat_t matrix;
    int matrix_initialized=0;
    switch(op) {
      case 0:scalar=nmod_poly_deflation(A);chunks=0;break;
      case 1:nmod_poly_deflate(R,A,d);break;
      case 2:nmod_poly_inflate(R,A,d);break;
      case 3:scalar=nmod_poly_remove(A,B);nmod_poly_set(R,A);break;
      case 4:case 5:
        vector=flint_malloc(sizeof(ulong)*(d+1));
        if(op==4)_nmod_vec_rand(vector,state,d,A->mod);else _nmod_vec_randtest(vector,state,d,A->mod);
        break;
      case 6:nmod_poly_randtest(R,state,d);break;
      case 7:scalar=nmod_poly_factor_equal_deg_prob(R,state,A,d);if(!scalar)chunks=0;break;
      case 8:nmod_poly_factor_equal_deg(factors,A,d);chunks=factors->num;break;
      case 9:nmod_poly_factor_distinct_deg(factors,A,&degrees);chunks=factors->num;break;
      case 10:nmod_poly_factor_cantor_zassenhaus(factors,A);chunks=factors->num;break;
      case 11:nmod_poly_factor_kaltofen_shoup(factors,A);chunks=factors->num;break;
      case 12:scalar=nmod_poly_factor(factors,A);chunks=factors->num;break;
      case 13:scalar=nmod_poly_factor_with_cantor_zassenhaus(factors,A);chunks=factors->num;break;
      case 14:scalar=nmod_poly_factor_with_kaltofen_shoup(factors,A);chunks=factors->num;break;
      case 15:
        nmod_poly_reverse(V,B,B->length);nmod_poly_inv_series(V,V,B->length);
        nmod_mat_init(matrix,n_sqrt(B->length-1)+1,B->length-1,p);matrix_initialized=1;
        nmod_poly_precompute_matrix(matrix,A,B,V);chunks=matrix->r;break;
      case 16:
        nmod_poly_reverse(V,C,C->length);nmod_poly_inv_series(V,V,C->length);
        nmod_mat_init(matrix,n_sqrt(C->length-1)+1,C->length-1,p);matrix_initialized=1;
        nmod_poly_precompute_matrix(matrix,B,C,V);
        nmod_poly_compose_mod_brent_kung_precomp_preinv(R,A,matrix,C,V);break;
      case 17:
        nmod_poly_reverse(V,C,C->length);nmod_poly_inv_series(V,V,C->length);
        polys=flint_malloc(sizeof(nmod_poly_struct)*d);results=flint_malloc(sizeof(nmod_poly_struct)*d);
        for(slong i=0;i<d;i++) {
          nmod_poly_init(polys+i,p);nmod_poly_init(results+i,p);
          for(slong j=0;j<C->length-1;j++)nmod_poly_set_coeff_ui(polys+i,j,a[i*(C->length-1)+j]);
        }
        nmod_poly_compose_mod_brent_kung_vec_preinv(results,polys,d,d,B,C,V);chunks=d;break;
      case 18: {
        nmod_mat_t left,right;
        nmod_mat_init(left,seed1,d,p);nmod_mat_init(right,d,seed2,p);
        nmod_mat_init(matrix,seed1,seed2,p);matrix_initialized=1;
        for(slong i=0;i<(slong)seed1;i++)for(slong j=0;j<d;j++)nmod_mat_entry(left,i,j)=a[i*d+j];
        for(slong i=0;i<d;i++)for(slong j=0;j<(slong)seed2;j++)nmod_mat_entry(right,i,j)=b[i*seed2+j];
        nmod_mat_mul(matrix,left,right);chunks=seed1;
        nmod_mat_clear(left);nmod_mat_clear(right);break;
      }
      case 19: {
        slong n=C->length-1;
        nmod_poly_reverse(V,C,C->length);nmod_poly_inv_series(V,V,C->length);
        nmod_mat_init(matrix,d,n,p);matrix_initialized=1;
        for(slong i=0;i<d;i++)for(slong j=0;j<n;j++)nmod_mat_entry(matrix,i,j)=a[i*n+j];
        nmod_poly_fit_length(B,n);for(slong i=B->length;i<n;i++)B->coeffs[i]=0;
        nmod_poly_fit_length(R,n);
        _nmod_poly_mod_matrix_rows_evaluate(R->coeffs,matrix,B->coeffs,n,C->coeffs,C->length,V->coeffs,V->length,A->mod);
        R->length=n;_nmod_poly_normalise(R);break;
      }
      case 20: {
        nmod_mat_t before;
        nmod_poly_reverse(V,B,B->length);nmod_poly_inv_series(V,V,B->length);
        nmod_mat_init(before,n_sqrt(B->length-1)+1,B->length-1,p);
        nmod_poly_precompute_matrix(before,A,B,V);
        _nmod_poly_reduce_matrix_mod_poly(matrix,before,C);matrix_initialized=1;chunks=matrix->r;
        nmod_mat_clear(before);break;
      }
    }
    if (!emit_word(out,cap,&len,scalar) || !emit_word(out,cap,&len,state->__randval)
        || !emit_word(out,cap,&len,state->__randval2) || !emit_word(out,cap,&len,chunks)) {len=-1;goto cleanup;}
    if(op==4 || op==5) {if(!emit_vec(out,cap,&len,vector,d,1))len=-1;}
    else if(op>=8 && op<=14) {
      for(slong i=0;i<chunks;i++)if(!emit_vec(out,cap,&len,factors->p[i].coeffs,factors->p[i].length,op==9?degrees[i]:factors->exp[i])){len=-1;break;}
    } else if(op==15 || op==18 || op==20) {
      for(slong i=0;i<chunks;i++)if(!emit_vec(out,cap,&len,nmod_mat_entry_ptr(matrix,i,0),matrix->c,1)){len=-1;break;}
    } else if(op==17) {
      for(slong i=0;i<chunks;i++)if(!emit_vec(out,cap,&len,results[i].coeffs,results[i].length,1)){len=-1;break;}
    } else if(chunks) {if(!emit_vec(out,cap,&len,R->coeffs,R->length,1))len=-1;}
cleanup:
    if(matrix_initialized)nmod_mat_clear(matrix);
    if(polys)for(slong i=0;i<d;i++){nmod_poly_clear(polys+i);nmod_poly_clear(results+i);}
    flint_free(polys);flint_free(results);flint_free(vector);flint_free(degrees);
    nmod_poly_clear(A);nmod_poly_clear(B);nmod_poly_clear(C);nmod_poly_clear(R);nmod_poly_clear(V);
    nmod_poly_factor_clear(factors);flint_rand_clear(state);return len;
}
