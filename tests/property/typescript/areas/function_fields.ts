import {mpqsInternals} from '../../../../packages/parigp-ts/src/mpqs.js';
import {ECM} from '../../../../packages/parigp-ts/src/_ecm.js';
import * as pariPower from '../../../../packages/parigp-ts/src/ifactor.js';
import {squfof,pollardbrent,Z_pollardbrent} from '../../../../packages/parigp-ts/src/ifactor.js';
import { eulerphiu as nativeWordTotient } from '../../../../packages/parigp-ts/src/arith2.js';
import { absZ_factor_limit_strict as pariStrictFactor } from '../../../../packages/parigp-ts/src/ifactor1.js';
import { znstar as pariZnstar } from '../../../../packages/parigp-ts/src/char.js';
import { subgrouplist as pariSubgroups } from '../../../../packages/parigp-ts/src/subgroup.js';
import { ZM_hnfmodid as pariHnfModId } from '../../../../packages/parigp-ts/src/hnf_snf.js';
import { pgener_Zp as pariPrimeGenerator } from '../../../../packages/parigp-ts/src/arith1.js';
import { Fp_powu as pariUnsignedPower } from '../../../../packages/parigp-ts/src/arith1.js';
import * as pariQPoly from '../../../../packages/parigp-ts/src/galconj.js';
import * as pariFpScalars from '../../../../packages/parigp-ts/src/ff.js';
import * as pariObservations from '../../../../packages/parigp-ts/src/galconj.js';
import { residue as pariScalarResidue, inverseCoefficient as pariScalarInverse } from '../../../../packages/parigp-ts/src/_polynomial_division.js';
import * as pariPrimeLinear from '../../../../packages/parigp-ts/src/ffinit.js';
import * as pariExtensionFactor from '../../../../packages/parigp-ts/src/FpXQX_factor.js';
import { Flx_Frobenius as pariWordFrobenius } from '../../../../packages/parigp-ts/src/Flx.js';
import { FpXQ_powBig as pariPrimePowerAlias, FpX_Frobenius as pariPrimeFrobenius } from '../../../../packages/parigp-ts/src/galconj.js';
import { extensionMinimalPolynomial } from '../../../../packages/parigp-ts/src/_extension_minpoly.js';
import { extensionProjection } from '../../../../packages/parigp-ts/src/_extension_projection.js';
import { extensionTruncatedPolynomial } from '../../../../packages/parigp-ts/src/_extension_polynomial.js';
import { FlxqX_dotproduct as extensionWordDot } from '../../../../packages/parigp-ts/src/Flx.js';
import { extensionAutomorphism } from '../../../../packages/parigp-ts/src/_extension_automorphism.js';
import { coefficientSubstitution } from '../../../../packages/parigp-ts/src/_coefficient_composition.js';
import { gen_pow_i, gen_pow_fold } from '../../../../packages/parigp-ts/src/bb_group.js';
import { FpXQ_pow as pari_gcd_power } from '../../../../packages/parigp-ts/src/ffinit.js';
import { FpX_split_part, FpX_nbroots } from '../../../../packages/parigp-ts/src/galconj.js';
import { FpX_halfgcd_all, FpX_halfgcd } from '../../../../packages/parigp-ts/src/FpX.js';
import { Flx_gcd, Flx_extgcd, Flx_halfgcd_all, Flx_halfgcd } from '../../../../packages/parigp-ts/src/Flx.js';
import { FpXQ_inv as pari_gcd_inverse, ZpX_liftfact, bezout_lift_fact } from '../../../../packages/parigp-ts/src/galconj.js';
import { FpX_gcd as pari_FpX_gcd } from '../../../../packages/parigp-ts/src/ffinit.js';
import { FpX_extgcd as pari_FpX_extgcd } from '../../../../packages/parigp-ts/src/galconj.js';
import { FpX_invBarrett } from '../../../../packages/parigp-ts/src/FpX.js';
import { Flx_divrem, Flx_rem, Flx_invBarrett } from '../../../../packages/parigp-ts/src/Flx.js';
import { FpX_divrem as pari_FpX_divrem, FpX_rem as pari_FpX_rem } from '../../../../packages/parigp-ts/src/ffinit.js';
import { ZX_sqr } from '../../../../packages/parigp-ts/src/ZX.js';
import { FpX_sqr } from '../../../../packages/parigp-ts/src/FpX.js';
import { Flx_mul, Flx_sqr } from '../../../../packages/parigp-ts/src/Flx.js';
import { ZX_mul as pari_ZX_mul } from '../../../../packages/parigp-ts/src/galconj.js';
import { FpX_FpXQV_eval as pari_FpXQV_eval, FpXQ_auttrace } from '../../../../packages/parigp-ts/src/FpX.js';
import { Flxq_powers, Flx_Flxq_eval, Flx_FlxqV_eval } from '../../../../packages/parigp-ts/src/Flx.js';
import { brent_kung_optpow } from '../../../../packages/parigp-ts/src/RgX.js';
import { FpX_FpXQ_eval as pari_FpX_eval } from '../../../../packages/parigp-ts/src/galconj.js';
import { FpX_mul as pari_FpX_mul, FpX_Fp_mul as pari_FpX_scalar } from '../../../../packages/parigp-ts/src/ffinit.js';
import { Flm_mul, FpM_mul } from '../../../../packages/parigp-ts/src/FpV.js';
import { ZM_mul as pari_ZM_mul } from '../../../../packages/parigp-ts/src/buch.js';
import { FpXQ_powers, FpXQ_autpow, FpXQ_autpowers } from '../../../../packages/parigp-ts/src/galconj.js';
import { F2x_factor, F2x_factor_squarefree, F2x_ddf } from '../../../../packages/parigp-ts/src/FpX_factor.js';
import { F2x_degree, F2x_add, F2x_mul, F2x_sqr, F2x_sqrt, F2x_rem, F2x_divrem, F2x_gcd, F2x_deriv, F2x_valrem, F2xq_mul, F2xq_sqr, F2xq_powers, F2x_Frobenius, F2x_matFrobenius } from '../../../../packages/parigp-ts/src/F2x.js';
import { F2m_mul, F2m_ker_sp, F2m_ker } from '../../../../packages/parigp-ts/src/F2v.js';
import { pari_init_rand, pari_rand, setrand, getrand, random_bits, random_Fl, randomi, random_F2x, random_zv } from '../../../../packages/parigp-ts/src/random.js';
import { random_Flx } from '../../../../packages/parigp-ts/src/Flx.js';
import { random_FpX } from '../../../../packages/parigp-ts/src/FpX.js';
import { nmod_mat_mul } from '../../../../packages/flint-ts/src/nmod_mat/mul.js';
import { _nmod_poly_mod_matrix_rows_evaluate } from '../../../../packages/flint-ts/src/nmod_poly/mod_matrix_rows_evaluate.js';
import { nmod_poly_factor_distinct_deg } from '../../../../packages/flint-ts/src/nmod_poly_factor/factor_distinct_deg.js';
import { nmod_poly_factor_kaltofen_shoup } from '../../../../packages/flint-ts/src/nmod_poly_factor/factor_kaltofen_shoup.js';
import { nmod_poly_factor, nmod_poly_factor_with_cantor_zassenhaus, nmod_poly_factor_with_kaltofen_shoup } from '../../../../packages/flint-ts/src/nmod_poly_factor/factor.js';
import { nmod_poly_factor_equal_deg_prob } from '../../../../packages/flint-ts/src/nmod_poly_factor/factor_equal_deg_prob.js';
import { nmod_poly_factor_equal_deg } from '../../../../packages/flint-ts/src/nmod_poly_factor/factor_equal_deg.js';
import { nmod_poly_factor_cantor_zassenhaus } from '../../../../packages/flint-ts/src/nmod_poly_factor/factor_cantor_zassenhaus.js';
import { _nmod_poly_inv_series_newton } from '../../../../packages/flint-ts/src/nmod_poly/inv_series_newton.js';
import { nmod_poly_precompute_matrix, nmod_poly_compose_mod_brent_kung_precomp_preinv, _nmod_poly_reduce_matrix_mod_poly } from '../../../../packages/flint-ts/src/nmod_poly/compose_mod_brent_kung_precomp_preinv.js';
import { nmod_poly_compose_mod_brent_kung_vec_preinv } from '../../../../packages/flint-ts/src/nmod_poly/compose_mod_brent_kung_vec_preinv.js';
import { nmod_poly_deflation } from '../../../../packages/flint-ts/src/nmod_poly/deflation.js';
import { nmod_poly_deflate } from '../../../../packages/flint-ts/src/nmod_poly/deflate.js';
import { nmod_poly_inflate } from '../../../../packages/flint-ts/src/nmod_poly/inflate.js';
import { nmod_poly_remove } from '../../../../packages/flint-ts/src/nmod_poly/remove.js';
import { nmod_poly_randtest } from '../../../../packages/flint-ts/src/nmod_poly/randtest.js';
import { _nmod_vec_rand, _nmod_vec_randtest } from '../../../../packages/flint-ts/src/nmod_vec/rand.js';
import { flint_rand_init, flint_rand_set_seed, flint_rand_get_seed, flint_rand_clear, n_randlimb, n_randint, n_urandint, n_randbits, n_randtest_bits, n_randtest, n_randtest_not_zero } from '../../../../packages/flint-ts/src/index.js';
import { RationalFunctionField } from '../../../../packages/sagemath-ts/src/rings/function_field/function_field_rational.js';
import { Matrix } from '../../../../packages/sagemath-ts/src/matrix/matrix_generic.js';
import { is_FunctionField, is_FunctionFieldElement } from '../../../../packages/sagemath-ts/src/rings/function_field/index.js';
import { Rational } from '../../../../packages/sagemath-ts/src/rings/rational.js';
import { MPolynomialRing } from '../../../../packages/sagemath-ts/src/rings/polynomial/multi_polynomial_ring.js';
import { constant_field_cardinality } from '../../../../packages/sagemath-ts/src/rings/function_field/constant_field.js';
import { _xmrange_iter } from '../../../../packages/sagemath-ts/src/misc/mrange.js';
import { FunctionFieldIdeal_rational, FunctionFieldIdealInfinite_rational } from '../../../../packages/sagemath-ts/src/rings/function_field/index.js';
import { divisor, prime_divisor, DivisorGroup } from '../../../../packages/sagemath-ts/src/rings/function_field/index.js';
import { _fmpz_poly_lcm, _fmpq_poly_lcm } from '../../../../packages/flint-ts/src/index.js';
import { PlaceSet, FunctionFieldValuationRing } from '../../../../packages/sagemath-ts/src/rings/function_field/index.js';
import { FpX_factor_squarefree } from '../../../../packages/parigp-ts/src/index.js';
import { GF2, GF2Element } from '../../../../packages/sagemath-ts/src/rings/finite_rings/gf2.js';
import { nmod_poly_factor_squarefree } from '../../../../packages/flint-ts/src/index.js';
import { GFpn } from '../../../../packages/sagemath-ts/src/rings/finite_rings/finite_field_extension.js';
import { Zmod } from '../../../../packages/sagemath-ts/src/rings/finite_rings/integer_mod_ring.js';
import { ZZX_GCD, ZZX_SquareFreeDecomp } from '../../../../packages/ntl-ts/src/index.js';
import { Integer, ZZ } from '../../../../packages/sagemath-ts/src/rings/integer_ring.js';
import { PolynomialRing } from '../../../../packages/sagemath-ts/src/rings/polynomial/polynomial_ring.js';
import { _nmod_poly_sqrt } from '../../../../packages/flint-ts/src/nmod_poly/sqrt.js';
import { _nmod_poly_sqrt_series } from '../../../../packages/flint-ts/src/nmod_poly/sqrt_series.js';
import { _gr_poly_sqrt_series_basecase } from '../../../../packages/flint-ts/src/gr_poly/sqrt_series_basecase.js';
import { _gr_poly_sqrt_series_newton } from '../../../../packages/flint-ts/src/gr_poly/sqrt_series_newton.js';
import { _gr_poly_rsqrt_series_basecase } from '../../../../packages/flint-ts/src/gr_poly/rsqrt_series_basecase.js';
import { _gr_poly_inv_series_basecase } from '../../../../packages/flint-ts/src/gr_poly/inv_series_basecase.js';
import { _nmod_poly_mulhigh } from '../../../../packages/flint-ts/src/nmod_poly/mulhigh.js';
import { _nmod_poly_mulhigh_classical } from '../../../../packages/flint-ts/src/nmod_poly/mulhigh_classical.js';
import { _n_jacobi_unsigned, n_jacobi_unsigned, n_jacobi, n_is_square, n_preinvert_limb, n_powmod2_ui_preinv, n_sqrtmod } from '../../../../packages/flint-ts/src/index.js';
/**
 * sagemath-ts side of the live rational-function-field differential area.
 */

import { QQ } from '../../../../packages/sagemath-ts/src/rings/rational_field.js';
import { GF } from '../../../../packages/sagemath-ts/src/rings/finite_rings/index.js';
import type {
  ConstantField,
  ConstantFieldElement,
} from '../../../../packages/sagemath-ts/src/rings/function_field/constant_field.js';
import { FunctionFieldElement_rational } from '../../../../packages/sagemath-ts/src/rings/function_field/element_rational.js';
import type { RationalFunctionField_global } from '../../../../packages/sagemath-ts/src/rings/function_field/function_field_rational.js';
import { IdealMonoid } from '../../../../packages/sagemath-ts/src/rings/function_field/ideal.js';
import { FunctionField } from '../../../../packages/sagemath-ts/src/rings/function_field/index.js';

type CE = ConstantFieldElement;
type FF = FunctionFieldElement_rational<CE>;
const bool = (value: boolean): string => (value ? 'True' : 'False');

function field(p: bigint): RationalFunctionField_global<CE> {
  return FunctionField(
    (p === 0n ? QQ : GF(p)) as unknown as ConstantField<CE>,
    'x'
  ) as RationalFunctionField_global<CE>;
}

function element(K: RationalFunctionField_global<CE>, coeffs: bigint[]): FF {
  const k = K.constant_base_field();
  return K.__call__(coeffs.map((c) => k.__call__(c))) as FF;
}

const ints = (xs: Array<{ value?: bigint; toString(): string }>): string =>
  xs.map((x) => (typeof x.value === 'bigint' ? x.value : x.toString())).join(',');

function parts(f: {numerator(): FF['_num']; denominator(): FF['_den']}): string {
  return `${ints(f.numerator().coeffs)}/${ints(f.denominator().coeffs)}`;
}

function ff_arithmetic(
  p: bigint,
  numerator: bigint[],
  denominator: bigint[],
  exponent: bigint
): string {
  const K = field(p);
  const f = element(K, numerator).div(element(K, denominator));
  const g = f.pow(exponent);
  return `f=${parts(f)} pow=${parts(g)} degree=${f.degree()} vx=${f.valuation(
    K.gen()
  )} square=${bool(f.is_square())}`;
}

function ff_factor(p: bigint, numerator: bigint[], denominator: bigint[]): string {
  const K = field(p);
  const f = element(K, numerator).div(element(K, denominator));
  const F = f.factor();
  const factors = F.factors.map(([q, e]) => `${parts(q)}:${e}`).join(';');
  return `unit=${parts(F.unit)} factors=${factors}`;
}

function ff_places(p: bigint, degree: bigint): string {
  return field(p).places(Number(degree)).map(String).join(';');
}

export const functions = {
  ff_arithmetic,
  ff_factor,
  ff_places,
  ff_parent,
  ff_element,
  ff_valuation,
  ff_zero_ideal,
  ff_bridge,
  ff_square,
  ff_pari_squarefree,
  ff_place_maps,
  ff_poly_numden,
  ff_poly_nested_numden,
  ff_poly_lcm,
  ff_divisor,
  ff_divisor_echelon,
  ff_ideal,
  ff_enumeration,
  ff_mrange,
  qq_cardinality: () => QQ.cardinality() === 'Infinity' ? '+Infinity' : String(QQ.cardinality()),
  ff_mpoly_repr,
  ff_element_observe,
  ff_predicate,
  ff_predicate_access, ff_matrix_base, ff_variable_map, ff_poly_dict, ff_poly_dict_order, ff_flint_rng, ff_flint_factor, ff_word_factor, ff_flint_matrix, ff_pari_rng, ff_pari_binary, ff_binary_factor, ff_pari_quotient, ff_pari_matrix, ff_pari_matrix_guard, ff_pari_composition, ff_pari_polymul, ff_pari_polymul_pattern, ff_pari_polydiv, ff_pari_polydiv_pattern, ff_pari_polygcd, ff_pari_polygcd_pattern, ff_pari_minpoly, ff_pari_ddf, ff_pari_sqrt, ff_pari_normalize, ff_pari_factor, ff_pari_roots, ff_pari_extension, ff_pari_extension_large, ff_pari_extension_division, ff_pari_extension_division_pattern, ff_pari_extension_gcd, ff_pari_extension_gcd_termination, ff_pari_extension_quotient, ff_pari_extension_matrix, ff_pari_extension_composition, ff_pari_inner_cache, ff_pari_aut_cache, ff_pari_coefficient_composition, ff_pari_extension_aut, ff_pari_extension_projection, ff_pari_extension_minpoly, ff_pari_prime_quotient_cache, ff_pari_extension_frobenius, ff_pari_prime_linear, ff_pari_extension_root_count, ff_pari_scalar_inverse, ff_pari_prime_observation, ff_pari_fp_scalar, ff_pari_fp_power, ff_pari_fp_predicate, ff_pari_qpoly, ff_pari_galois_integer, ff_pari_zx_index, ff_pari_zp_lift, ff_pari_vandermonde, ff_pari_perm_vector, ff_pari_group_structure, ff_pari_sympol, ff_pari_subgroup, ff_pari_subgroup_bounds, ff_pari_galois_validation, ff_pari_galois_actions, ff_pari_symmetric_search, ff_pari_findpsi, ff_pari_galois_filter, ff_pari_galois_helpers, ff_pari_galois_kernel, ff_pari_fixedfield_reprime, ff_pari_factor_stages, ff_pari_power_helpers, ff_pari_ecm, ff_pari_mpqs_relations, ff_pari_mpqs_sqrt, ff_pari_sparse_kernel, ff_pari_mpqs_symbols, ff_pari_mpqs_debug, ff_pari_mpqs_inverse, ff_pari_mpqs_init, ff_pari_mpqs_fb, ff_pari_mpqs_candidates, ff_pari_mpqs_driver, ff_pari_mpqs_hash, ff_pari_mpqs_warning, ff_pari_mpqs_class_candidates,
  ff_gf2_roots,
  ff_flint_squarefree,
  ff_squarefree_field,
  ff_squarefree,
  ff_ntl_poly,
  ff_polynomial_square,
  ff_flint_poly,
  ff_flint_word,
  ff_order,
  ff_inverse,
  ff_wrapped,
  ff_fraction_strings,
};


function ff_parent(p: bigint, op: bigint, n: bigint): string {
  const K=field(p), k=K.constant_base_field();
  if(op===0n)return [K.toString(),K.variable_name(),K.variable_names().join(','),K.ngens(),K.degree(),K.genus(),K.characteristic(),bool(K.is_finite()),bool(K.is_global()),bool(K.is_perfect()),bool(K.is_field()),bool(K.base_field()===K),bool(K.rational_function_field()===K),bool(K.constant_field()===k),bool(K.constant_base_field()===k)].join('|');
  if(op===1n)return parts(K.gen(Number(n)) as FF);
  if(op===2n)return String(K.degree(n===0n?undefined:n===1n?K:FunctionField(k,n===2n?'x':'y')));
  if(op===3n)return bool(K===FunctionField(k,n===0n?'x':['x']));
  if(op===4n){const [L,from,to]=K.change_variable_name(n===0n?'x':'y');return [String(L),bool(L.change_variable_name('x')[0]===K),parts(from(L.gen())),parts(to(K.gen())),bool(from(to(K.gen())).eq(K.gen()))].join('|');}
  if(op===5n)return String(K.different());
  if(op===6n)return parts(K.maximal_order().gen(Number(n)));
  throw new Error('unknown ff_parent operation');
}

function ff_element(p: bigint,op: bigint,num: bigint[],den: bigint[],num2: bigint[],den2: bigint[],n: bigint): string {
  const K=field(p), f=element(K,num).div(element(K,den)), g=element(K,num2).div(element(K,den2));
  if(op===0n)return [parts(f),bool(f.is_zero()),bool(f.is_one()),f.degree(),parts(f.trace()),parts(f.norm()),parts(f.matrix()[0]![0]!),parts(f.list()[0]!),parts(f.element()),bool(f.parent===K)].join('|');
  if(op===1n)return parts(f.neg());
  if(op===2n)return parts(f.inv());
  if(op===3n)return parts(f.pow(n));
  if(op===4n)return parts(f.sqrt() as FF);
  if(op===5n)return (f.sqrt(true) as FF[]).map(parts).join(';');
  if(op===6n)return bool(f.is_square());
  if(op===7n)return String(f.cmp(g));
  if(op===8n)return parts(f.add(g));
  if(op===9n)return parts(f.sub(g));
  if(op===10n)return parts(f.mul(g));
  if(op===11n)return parts(f.div(g));
  if(op===12n)return parts(f.nth_root(n));
  if(op===13n)return bool(f.is_nth_power(n));
  if(op===14n)return f.toString();
  if(op===15n)return parts(f.scalar_mul(K.constant_base_field().__call__(n)));
  throw new Error('unknown ff_element operation');
}


function ff_valuation(p: bigint,num: bigint[],den: bigint[],vnum: bigint[],vden: bigint[],mode: bigint): string {
  const K=field(p), f=element(K,num).div(element(K,den));
  const v=element(K,vnum).div(element(K,vden));
  const result=f.valuation(mode===0n?v:v.numerator());
  return result===Number.POSITIVE_INFINITY?'+Infinity':String(result);
}


function ff_zero_ideal(p: bigint,infinite: bigint,op: bigint): string {
  const K=field(p), O=infinite?K.maximal_order_infinite():K.maximal_order(), I=O.ideal(K.zero());
  if(op===0n)return bool(I.is_zero());
  if(op===1n)return I.factor().map(([P,e])=>`${parts(P.gens()[0] as FF)}:${e}`).join(';');
  return String(op===2n?I.divisor():op===3n?I.divisor_of_zeros():I.divisor_of_poles());
}


function ff_bridge(p: bigint,q: bigint,op: bigint,mode: bigint,num: bigint[],den: bigint[],codes: bigint[]): string {
  const K=field(p);
  if(op===0n)return String(K.field());
  if(op===5n||op===6n){const f=K.__call__(String.fromCodePoint(...codes.map(Number)));if(op===6n)return [bool(f.element().parent===K.field()),bool(f.numerator().parent===K._ring),bool(f.denominator().parent===K._ring)].join('|');return parts(f)+'|'+String(f);}
  const f=element(K,num).div(element(K,den));
  if(op===1n){const value=f.element();return [String(value.parent),bool(value.parent===K.field()),bool(value===f.element()),value.constructor.name].join('|');}
  if(op===2n){const value=K.__call__(f.element());return parts(value)+'|'+String(value);}
  if(op===3n){const R=K._ring, k=K.constant_base_field();const value=K.__call__(R.fraction_field().__call__(R.__call__(num.map(x=>k.__call__(x))),R.__call__(den.map(x=>k.__call__(x)))));return parts(value)+'|'+String(value);}
  const L=FunctionField(q===p?K.constant_base_field():(q===0n?QQ:GF(q)) as unknown as ConstantField<CE>,mode===0n?'x':'y');
  const k=L.constant_base_field(), value=K.__call__(L.__call__(num.map(x=>k.__call__(x))).div(L.__call__(den.map(x=>k.__call__(x)))));
  return parts(value)+'|'+String(value);
}

function ff_fraction_strings(p: bigint,mode: bigint,ncodes: bigint[],dcodes: bigint[]): string {
  const F=field(p).field(),n=String.fromCodePoint(...ncodes.map(Number)),d=String.fromCodePoint(...dcodes.map(Number));
  const f=mode===0n?F.__call__(n):F.__call__(n,d);return parts(f)+'|'+String(f);
}

function ff_wrapped(p: bigint,op: bigint,num: bigint[],den: bigint[]): string {
  const K=field(p),F=K.field(),R=K._ring,k=K.constant_base_field();
  const raw=new F._element_class(F,R.__call__(num.map(x=>k.__call__(x))),R.__call__(den.map(x=>k.__call__(x))),{reduce:false});
  const f=new FunctionFieldElement_rational(K,raw),before=parts(f);
  if(op===0n)return [before,bool(f.is_zero()),bool(f.is_one()),bool(f.eq(K.zero())),bool(f.eq(K.one())),String(f),bool(f.element()===raw)].join('|');
  if(op===1n){let err='ok';try{R.__call__(raw);}catch(e){err=(e as Error).name+':'+(e as Error).message;}return before+'|'+err+'|'+parts(f)+'|'+String(f);}
  if(op===2n){if('reduce'in raw)(raw as typeof raw & {reduce():void}).reduce();return before+'|'+parts(f)+'|'+String(f);}
  const g=op===3n?f.neg():op===4n?f.add(K.one()):op===5n?f.sub(K.one()):op===6n?f.mul(K.gen()):op===7n?f.div(K.gen()):f.pow(op-8n);
  return parts(g)+'|'+String(g);
}

function ff_inverse(p: bigint,mode: bigint,num: bigint[],den: bigint[],inum: bigint[],iden: bigint[]): string {
  const K=field(p);let f=element(K,num).div(element(K,den)),g=element(K,inum).div(element(K,iden));
  if(mode===4n||mode===5n){const F=K.field(),R=K._ring,k=K.constant_base_field(),raw=new F._element_class(F,R.__call__((mode===4n?num:inum).map(x=>k.__call__(x))),R.__call__((mode===4n?den:iden).map(x=>k.__call__(x))),{reduce:false});if(mode===4n)f=new FunctionFieldElement_rational(K,raw);else g=new FunctionFieldElement_rational(K,raw);}
  const I=mode<2n?(mode===0n?K.maximal_order():K.maximal_order_infinite()).ideal(g):{gens:()=>mode===2n?[g,g]:[g]};
  const v=f.inverse_mod(I as Parameters<FF['inverse_mod']>[0]);return parts(v)+'|'+String(v);
}

function ff_order(p: bigint,infinite: bigint,op: bigint,mode: bigint,num: bigint[],den: bigint[],codes: bigint[]): string {
  const K=field(p),O=infinite?K.maximal_order_infinite():K.maximal_order();
  if(op===0n){const M=O.ideal_monoid();return [String(O),bool(O.function_field()===K),bool(O.fraction_field()===K),bool(O.is_field()),bool(O.is_noetherian()),O.basis().map(String).join(','),String(M),bool(M.ring()===O),bool(M===O.ideal_monoid()),bool(M===new IdealMonoid(O))].join('|');}
  if(op===7n)return O.basis().map(x=>(typeof x==='bigint'?'int':x.constructor.name)+':'+String(x)).join(';');
  if(op===1n){const f=O.__call__(String.fromCodePoint(...codes.map(Number)));return parts(f)+'|'+String(f);}
  if(op===5n)return bool(O.is_subring(mode===0n?K:mode===1n?O:infinite?K.maximal_order():K.maximal_order_infinite()));
  const f=element(K,num).div(element(K,den));
  if(op===2n)return parts(O.__call__(f))+'|'+String(O.__call__(f));
  if(op===6n)return bool(O.contains(f));
  const P=mode===0n?O:infinite?K.maximal_order():K.maximal_order_infinite(),I=P.ideal(f);
  if(op===3n)return String(O.ideal(I));
  if(op===4n){const input=mode===2n?{gens:()=>[f]}:I,v=O.ideal_monoid().__call__(input);return String(v)+'|'+bool(v===input);}
  throw new Error('unknown ff_order operation');
}

function ff_flint_word(op: bigint,a: bigint,p: bigint,e: bigint,sign: bigint): string {
  if(op===0n)return String(n_sqrtmod(a,p));
  if(op===1n)return String(n_jacobi_unsigned(a,p));
  if(op===2n)return String(_n_jacobi_unsigned(a,p,Number(sign)));
  if(op===3n)return bool(n_is_square(a));
  if(op===4n)return String(n_powmod2_ui_preinv(a,e,p,n_preinvert_limb(p)));
  if(op===5n)return String(n_preinvert_limb(p));
  if(op===6n)return String(n_jacobi(a,p));
  throw new Error('unknown native word operation');
}

function ff_flint_poly(p: bigint,op: bigint,n: bigint,cutoff: bigint,a: bigint[],b: bigint[]): string {
  const value=op===0n?_nmod_poly_sqrt(a,p):op===1n?_nmod_poly_sqrt_series(a,Number(n),p):op===2n?_gr_poly_sqrt_series_basecase(a,Number(n),p):op===3n?_gr_poly_sqrt_series_newton(a,Number(n),Number(cutoff),p):op===4n?_gr_poly_rsqrt_series_basecase(a,Number(n),p):op===5n?_gr_poly_inv_series_basecase(a,Number(n),p):op===6n?_nmod_poly_mulhigh(a,b,Number(n),p):_nmod_poly_mulhigh_classical(a,b,Number(n),p);
  return value===null?'None':value.join(',');
}

function ff_square(p: bigint,op: bigint,raw: bigint,allRoots: bigint,extend: bigint,num: bigint[],den: bigint[]): string {
  const K=field(p),F=K.field(),R=K._ring,k=K.constant_base_field();
  const a=R.__call__(num.map(x=>k.__call__(x))),b=R.__call__(den.map(x=>k.__call__(x)));
  const x=(raw?new F._element_class(F,a,b,{reduce:false}):F.__call__(a,b)) as ReturnType<typeof F.__call__> & {is_square(root?:boolean):boolean|[boolean,ReturnType<typeof F.__call__>|null];sqrt(extend:boolean,all:boolean):ReturnType<typeof F.__call__>|ReturnType<typeof F.__call__>[];_sqrt_or_None():ReturnType<typeof F.__call__>|null};
  const f=new FunctionFieldElement_rational(K,x);
  if(op===0n)return bool(f.is_square());
  if(op===2n)return bool(x.is_square() as boolean);
  if(op===4n){const v=x._sqrt_or_None();return v===null?'None':parts(v)+'|'+String(v)+'|'+bool(v===x);}
  if(op===5n){const [ok,v]=x.is_square(true) as [boolean,ReturnType<typeof F.__call__>|null];return bool(ok)+'|'+(v===null?'None':parts(v));}
  const v=op===1n?f.sqrt(Boolean(allRoots)):x.sqrt(Boolean(extend),Boolean(allRoots));
  if(Array.isArray(v))return v.map(z=>parts(z)+'|'+String(z)).join(';');
  return parts(v)+'|'+String(v)+'|'+bool(v===(op===1n?f:x));
}

function ff_polynomial_square(p: bigint,root: bigint,a: bigint[]): string {
  const R=new PolynomialRing<any>(p===-1n?{zero:()=>new Integer(0n),one:()=>new Integer(1n),__call__:(x:unknown)=>new Integer(ZZ.__call__(x as bigint)),is_field:()=>false,characteristic:()=>0n,toString:()=> 'Integer Ring'}:p===0n?QQ:GF(p),'x'),f=R.__call__(a);
  if(!root)return bool(f.is_square());
  const [ok,g]=f.is_square(true);
  return bool(ok)+'|'+(g===null?'None':ints(g.coeffs)+'|'+String(g)+'|'+bool(g===f)+'|'+bool(g.parent===R));
}

function ff_ntl_poly(op: bigint,a: bigint[],b: bigint[]): string {
  return op===0n?ZZX_GCD(a,b).join(','):ZZX_SquareFreeDecomp(a).map(([f,e])=>f.join(',')+':'+e).join(';');
}

function ff_squarefree(p: bigint,a: bigint[]): string {
  const base=p===-1n?{zero:()=>new Integer(0n),one:()=>new Integer(1n),__call__:(x:unknown)=>new Integer(ZZ.__call__(x as bigint)),is_field:()=>false,characteristic:()=>0n,toString:()=> 'Integer Ring'}:p===0n?QQ:Zmod(p);
  const R=new PolynomialRing<any>(base,'x'),f=R.__call__(a),fac=f.squarefree_decomposition();
  let unit=R.base_ring.one();
  for(const [g,e]of fac)if(g.degree()<=0)unit=unit.mul(g.getCoeff(0).pow(BigInt(e)));
  return String(unit)+'|'+fac.filter(([g])=>g.degree()>0).map(([g,e])=>ints(g.coeffs)+':'+e+':'+bool(g.parent===R)).join(';');
}

function ff_flint_squarefree(p: bigint,a: bigint[]): string {
  return nmod_poly_factor_squarefree(a,p).map(([f,e])=>f.join(',')+':'+e).join(';');
}
function ff_squarefree_field(p: bigint,n: bigint,op: bigint,e: bigint,a: bigint[]): string {
  const mod=Number(p)===2?[1n,1n,0n]:Number(p)===3?[1n,0n]:Number(p)===5?[2n,0n]:[1n,0n];
  const k:any=n===1n?GF(p):GFpn(p,Number(n),mod,'a'),R=new PolynomialRing<any>(k,'x');
  let f=R.__call__(a.map(v=>n===1n?k.__call__(v):k.fromInteger(v%(p**n))));
  if(op)f=f.pow(Number(e)) as typeof f;
  if(op===2n)f=f.mul(R.gen().add(R.__call__(k.gen())).pow(Number(e+1n)) as typeof f);
  const fac=f.squarefree_decomposition();let unit=k.one();
  for(const [g,m]of fac)if(g.degree()<=0)unit=unit.mul(g.getCoeff(0).pow(BigInt(m)));
  const enc=(c:any)=>{const xs:bigint[]=typeof c.value==='bigint'?[c.value]:c.coefficients().map((z:any)=>z.value);while(xs.at(-1)===0n)xs.pop();return xs.join(',');};
  return enc(unit)+'|'+fac.filter(([g])=>g.degree()>0).map(([g,m])=>g.coeffs.map(enc).join('/')+':'+m+':'+bool(g.parent===R)).join(';');
}

function ff_gf2_roots(op: bigint,a: bigint,b: bigint,raw: bigint,allRoots: bigint,extend: bigint): string {
  const x:any=raw?new GF2Element(a,GF2):GF2.__call__(a);
  if(op===13n)return bool(x.isOne())+'|'+x.repr();
  if(op===0n)return bool(x.is_square());
  if(op===2n){const R=new PolynomialRing<any>(GF2,'x'),f=R.__call__([a,b,a]),[ok,r]=f.is_square(true);return bool(ok)+'|'+(r===null?'None':ints(r.coeffs)+'|'+bool(r.parent===R));}
  if(op===3n)return bool(GF2.zero()===GF2.__call__(0n))+'|'+bool(GF2.one()===GF2.__call__(1n))+'|'+bool(GF2.gen()===GF2.__call__(1n));
  if(op===12n)return x.sqrt({oops:true});
  let v:any;
  if(op===1n)v=x.sqrt({all:Boolean(allRoots),extend:Boolean(extend)});
  else if(op===4n)v=x.add(GF2.__call__(b));else if(op===5n)v=x.sub(GF2.__call__(b));else if(op===6n)v=x.mul(GF2.__call__(b));else if(op===7n)v=x.div(GF2.__call__(b));else if(op===8n)v=x.neg();else if(op===9n)v=x.inv();else if(op===10n)v=x.pow(b);else v=GF2.__call__(x);
  const enc=(v:GF2Element)=>String(v)+'|'+bool(v===x)+'|'+bool(v===GF2.__call__(BigInt(v.value)));
  return Array.isArray(v)?v.map(enc).join(';'):enc(v);
}

function ff_pari_squarefree(p: bigint,a: bigint[]): string {
  try{return 'value|'+FpX_factor_squarefree(a,p).map(f=>f.join(',')).join(';');}
  catch(e){return 'error|'+(e as Error).name+'|'+(e as Error).message;}
}

function ff_place_maps(p:bigint,infinite:bigint,a:bigint,op:bigint,mode:bigint,num:bigint[],den:bigint[]):string{
  const K=field(p),x=K.gen(),O=infinite?K.maximal_order_infinite():K.maximal_order();
  const I=O.ideal(infinite?x.inv():x.sub(K.__call__(a))),P=I.place(),S=K.place_set(),V=P.valuation_ring();
  const join=(xs:unknown[])=>xs.map(x=>typeof x==='boolean'?bool(x):String(x)).join('|');
  if(op===0n)return join([S,S.function_field()===K,S===new PlaceSet(K),P,P.degree(),P.is_infinite_place(),parts(P.local_uniformizer()as FF),P.parent()===S,P.function_field()===K,P.prime_ideal()===I,V,V.place().eq(P)]);
  if(op===1n){const Q=S.__call__(I),W=new FunctionFieldValuationRing(K,Q);return join([P===Q,P.eq(Q),V===P.valuation_ring(),V===Q.valuation_ring(),V===W,W.place()===V.place()]);}
  if(op===2n){const A=P.residue_field(),B=V.residue_field(),C=P.residue_field('a'),D=V.residue_field('a');return join([A[0],A===B,A[1]===B[1],A[2]===B[2],C===D,A===C,A[0]===K.constant_base_field()]);}
  if(op===6n)return bool(S.__call__(P)===P);
  if(op===7n){const Q=K.maximal_order().ideal(x.sub(K.__call__(a+1n))).place();return String(P.cmp(Q));}
  const f=element(K,num).div(element(K,den));
  if(op===5n){const g=V.__call__(f);return parts(g as FF)+'|'+bool(g===f);}
  const [k,fr,to]=mode===0n?P._residue_field():mode===1n?P.residue_field():V.residue_field('b');
  if(op===3n){const g=fr(K.constant_base_field().__call__(a));return parts(g as FF)+'|'+bool(g.parent===K);}
  if(op===4n)return String(to(f));
  if(op===8n)return String(to((num[0]??0n)as unknown as FF));
  throw new Error('unknown ff_place_maps operation');
}

function ff_poly_numden(p:bigint,a:bigint[],b:bigint[]):string{
  const k=(p===-1n?{zero:()=>new Integer(0n),one:()=>new Integer(1n),__call__:(x:unknown)=>new Integer(ZZ.__call__(x as bigint)),is_field:()=>false,characteristic:()=>0n,toString:()=> 'Integer Ring'}:p===0n?QQ:GF(p))as unknown as ConstantField<CE>,R=new PolynomialRing(k,'x');
  const f=R.__call__(a.map((n,i)=>p===-1n?k.__call__(n):(k.__call__(n)as CE&{div(x:CE):CE}).div(k.__call__(b[i]!))));
  const n=f.numerator(),d=f.denominator();
  return [n.coeffs.map(String).join(','),String(d),String(n.parent.base_ring),bool(n.parent===f.parent),bool(n.parent===f.numerator().parent)].join('|');
}

function ff_poly_nested_numden(mode:bigint,a:bigint[],b:bigint[],e:bigint[]):string{
  const T=new PolynomialRing<any>(mode!==2n?QQ:GF(5n),'t'),t=T.gen(),k:any=mode===0n?T:T.fraction_field(),R=new PolynomialRing<any>(k,'x');
  const cs=a.map((n,i)=>mode===0n?k.__call__(T.__call__([QQ.__call__(n).div(QQ.__call__(b[i]!)),QQ.one()])):k.__call__(n).div(k.__call__(t.add(T.__call__(b[i]!))).pow(e[i]!)));
  const f=R.__call__(cs),n=f.numerator(),d=f.denominator();return [String(n),String(d),bool(n.parent===f.parent)].join('|');
}

function ff_poly_lcm(p:bigint,a:bigint[],b:bigint[]):string{
  const trim=(input:bigint[])=>{const f=input.slice();while(f.length&&f[f.length-1]===0n)f.pop();return f;};
  if(p===-3n)return _fmpz_poly_lcm(trim(a),trim(b)).join(',');
  if(p===-2n){const[n,d]=_fmpq_poly_lcm(trim(a),trim(b));return n.join(',')+'/'+d;}
  const base=p===-1n?{zero:()=>new Integer(0n),one:()=>new Integer(1n),__call__:(x:unknown)=>new Integer(ZZ.__call__(x as bigint)),is_field:()=>false,characteristic:()=>0n,toString:()=> 'Integer Ring'}:p===0n?QQ:GF(p);
  const R=new PolynomialRing<any>(base,'x');return R.__call__(a).lcm(R.__call__(b)).coeffs.map(String).join(',');
}

function ff_divisor(p:bigint,op:bigint,keys:bigint[],ms:bigint[],ns:bigint[],scalar:bigint,num:bigint[],den:bigint[],vec:bigint[]):string{
  const K=field(p),x=K.gen(),G=K.divisor_group();
  const place=(a:bigint)=>a===-99n?K.maximal_order_infinite().ideal(x.inv()).place():K.maximal_order().ideal(x.sub(K.__call__(a))).place();
  const ps=keys.map(place),D=divisor(K,ps.map((q,i)=>[q,ms[i]!])),E=divisor(K,ps.slice(0,ns.length).map((q,i)=>[q,ns[i]!])),P=ps[0]??place(0n);
  const ds=(d:ReturnType<typeof G.zero>)=>String(d)+'|'+d.list().map(([q,m])=>q+':'+m).join(';');
  const join=(xs:unknown[])=>xs.map(x=>typeof x==='boolean'?bool(x):String(x)).join('|');
  if(op===0n)return join([ds(D),D.degree(),D.is_effective(),D.support().map(String).join(','),D.parent()===G,D.function_field()===K,D.dict()===D.dict(),G,G===new DivisorGroup(K),G.function_field()===K]);
  if(op===1n)return ds(D.neg());if(op===2n)return ds(D.add(E));if(op===3n)return ds(D.sub(E));if(op===4n)return ds(D.scalar_mul(scalar));
  if(op===5n)return String(D.cmp(E));if(op===6n)return ds(D.numerator())+'|'+ds(D.denominator());if(op===7n)return D.multiplicity(place(scalar))+'|'+D.valuation(place(scalar));
  if(op===8n)return ds(G.__call__(D))+'|'+bool(G.__call__(D)===D);if(op===9n)return ds(G.__call__(P));if(op===10n)return ds(G.__call__(vec.length?K.__call__(scalar):scalar));
  if(op===11n)return D.dimension()+'|'+D._basis().map(parts).join(';')+'|'+D.basis_function_space().map(parts).join(';');
  if(op===12n)return D.function_space()[2](element(K,num).div(element(K,den))).map(String).join(',');
  if(op===13n){const v=D.function_space()[1](vec.map(n=>K.constant_base_field().__call__(n)))as unknown;return v instanceof FunctionFieldElement_rational?'ff|'+parts(v):'scalar|'+v;}
  if(op===14n)return (D._format as unknown as (formatter:(x:unknown)=>string,mul:string,cr:string)=>string)(x=>'<'+x+'>',' @ ',' / ');
  if(op===15n){const A=D._function_space(),B=D._function_space(),C=D.function_space(),E=D.function_space();return join([A===B,A[0]===D.basis_function_space(),A[1]===B[1],C===E,C[1]===E[1],C[2]===E[2]]);}
  if(op===16n)return D._differential_space()[0].map(parts).join(';');if(op===17n)return ds(prime_divisor(K,P,scalar));if(op===18n){D.dict().set(P._key(),[P,scalar]);return ds(D);}
  if(op===19n)return join([G.zero()===G.zero(),G.__call__(0n)===G.__call__(0n),G.__call__(G.zero())===G.zero()]);
  throw new Error('unknown ff_divisor operation');
}

function ff_divisor_echelon(p:bigint,modes:bigint[],num:bigint[],den:bigint[]):string{
  const K=field(p),x=K.gen(),D=K.divisor_group().zero(),choices=[K.zero(),K.one(),x,x.add(K.one()),x.mul(K.__call__(2n)),x.inv(),x.add(K.one()).inv(),x.div(x.add(K.one()))];
  const[bs,co]=D._echelon_basis(modes.map(m=>choices[Number(m)]!));return bs.map(parts).join(';')+'|'+co(element(K,num).div(element(K,den))).map(String).join(',');
}

function ff_ideal(p:bigint,infinite:bigint,raw:bigint,op:bigint,num:bigint[],den:bigint[],num2:bigint[],den2:bigint[],e:bigint):string{
  const K=field(p),O=infinite?K.maximal_order_infinite():K.maximal_order(),f=element(K,num).div(element(K,den)),g=element(K,num2).div(element(K,den2));
  const cls=infinite?FunctionFieldIdealInfinite_rational:FunctionFieldIdeal_rational;
  const I:any=raw?new cls(O,f):O.ideal(f),J:any=raw?new cls(O,g):O.ideal(g);
  const fmt=(I:any)=>String(I)+'|'+parts(I.gen()),join=(xs:unknown[])=>xs.map(x=>typeof x==='boolean'?bool(x):String(x)).join('|');
  if(op===0n)return join([fmt(I),I._repr_short(),I.is_zero(),I.is_prime(),I.ring()===O,I.base_ring()===O,I.parent()===O.ideal_monoid(),I.gen()===I.gens()[0],I.gen()===I.gens_over_base()[0],I.gen()===f]);
  if(op===1n)return String(I.cmp(J));if(op===2n)return fmt(I.add(J));if(op===3n)return fmt(I.mul(J));if(op===4n)return fmt(I.div(J));if(op===5n)return fmt(I.inv());
  if(op===6n){const v=I.pow(e);return fmt(v)+'|'+bool(v===I);}if(op===7n)return bool(I.contains(g));if(op===8n){const v=I.valuation(J);return v===Infinity?'+Infinity':String(v);}if(op===9n)return String(I.place());if(op===10n)return fmt(I.acted_upon(g));if(op===11n)return bool(I.contains(e));if(op===12n)return I.denominator().coeffs.map(String).join(',');if(op===13n)return fmt(O.ideal([f,g]));
  if(op===14n)return I.factor().map(([P,m]:any)=>parts(P.gen())+':'+m).join(';');if(op===15n)return String(I.divisor());if(op===16n)return String(I.divisor_of_zeros());if(op===17n)return String(I.divisor_of_poles());if(op===18n)return I.gens_reduced().map(parts).join(';');if(op===19n)return join([I.pow(0n)===I.pow(0n),I.pow(1n)===I]);
  if(op===20n)return join([I.parent().one(),I.parent().one()===I.parent().one(),I.pow(0n)===I.parent().one()]);
  throw new Error('unknown ff_ideal operation');
}

function ff_enumeration(p:bigint,kind:bigint,op:bigint,mode:bigint,d:bigint,limit:bigint):string{
 const k:any=p===0n?QQ:kind===1n?Zmod(p):kind===2n?GFpn(p,2,p===2n?[1n,1n]:p===3n?[1n,0n]:[2n,0n],'a'):GF(p),R:any=new PolynomialRing(k,'x');
 if(op===11n)return String(constant_field_cardinality(k));
 if(op===10n){const n=k.cardinality();return n==='Infinity'?'+Infinity':String(n);}
 const opts=mode===0n?{of_degree:Number(d)}:mode===1n?{max_degree:Number(d)}:mode===2n?{}:{of_degree:Number(d),max_degree:Number(d)};
 if(op===6n)return R.monomial(Number(d)).coeffs.map(String).join(',');
 const take=(g:Iterator<any>)=>{const out:any[]=[];for(let i=0n;i<limit;i++){const r=g.next();if(r.done)break;out.push(r.value);}return out;};
 if(op>=7n){const K:any=FunctionField(k,'x');if(op===8n)return String(K.get_place(Number(d)));return (op===7n?take(K._places_finite(Number(d))):K.places(Number(d))).map(String).join(';');}
 const g=op===0n?R.polynomials(opts):op===1n?R.monics(opts):op===2n?R._polys_degree(Number(d)):op===3n?R._polys_max(Number(d)):op===4n?R._monics_degree(Number(d)):R._monics_max(Number(d));
 return take(g).map(f=>f.coeffs.map(String).join(',')).join(';');
}

function ff_mrange(sizes:bigint[],kind:bigint,mode:bigint,limit:bigint):string{
 const counts=sizes.map(()=>0);
 const fail=(name:string,message=''):never=>{const e=new Error(message);e.name=name;throw e;};
 let ins:any[]=sizes.map((n,i)=>kind===0n?Array.from({length:Number(n)},(_,v)=>v):{is_finite:()=>{if([4n,5n,7n,10n].includes(kind))return fail('AttributeError');if(kind===3n)return fail('ValueError','finiteness unknown');if(kind===6n)return fail('RuntimeError','finiteness unknown');return [8n,9n].includes(kind)?false:n>=0n;},cardinality:()=>{if(kind===4n)return n<0n?Infinity:n;if(kind===5n)return fail('NotImplementedError');if(kind===8n)return fail('TypeError','cardinality unknown');if([9n,10n].includes(kind))return fail('ValueError','cardinality unknown');return fail('AttributeError');},get length(){return n<0n?undefined:Number(n);},*[Symbol.iterator](){for(let v=0;n<0n||BigInt(v)<n;v++){counts[i]++;yield v;}}});
 if(kind===2n)ins=ins.map(x=>x[Symbol.iterator]());
 const g=_xmrange_iter<number,any>(ins,mode===0n?(r=[])=>r.slice():mode===1n?r=>r===undefined?777:r.reduce((a,b)=>a+b,0):r=>r===undefined?[999]:r),rows:any[]=[];
 for(let i=0n;i<limit;i++){const r=g.next();if(r.done)break;rows.push(r.value);}
 return rows.map(r=>mode===1n?String(r):r.join(',')).join(';')+'|'+bool(mode!==1n&&rows.length>1&&rows[0]===rows.at(-1))+'|'+counts.join(',');
}

function ff_mpoly_repr(p:bigint,order:bigint,cs:bigint[]):string{
 const k:any=p===0n?QQ:GF(p),R=new MPolynomialRing<any>(k,['x0','x1','x2'],(['lex','deglex','degrevlex']as const)[Number(order)]),es=[[0,0,0],[2,0,0],[0,1,1],[1,1,0],[0,0,2]];
 return String(R.fromTerms(cs.map((c,i)=>[k.__call__(c),es[i]!] as [any,number[]])));
}

function ff_element_observe(p:bigint,raw:bigint,op:bigint,num:bigint[],den:bigint[],a:bigint):string{
 const K=field(p),F=K.field(),R=K._ring,k=K.constant_base_field();
 const f=raw?new FunctionFieldElement_rational(K,new F._element_class(F,R.__call__(num.map(x=>k.__call__(x))),R.__call__(den.map(x=>k.__call__(x))),{reduce:false})):element(K,num).div(element(K,den));
 const join=(xs:unknown[])=>xs.map(x=>typeof x==='boolean'?bool(x):String(x)).join('|');
 if(op===0n)return join([parts(f),parts(f.matrix()[0]![0]!),parts(f.trace()),parts(f.norm()),f.matrix()[0]![0]===f,f.trace()===f,f.norm()===f,f.matrix()===f.matrix()]);
 if(op===1n){const D=f.divisor();return String(D)+'|'+bool(D===f.divisor())+'|'+D.degree();}
 if(op===2n)return String(f.divisor_of_zeros());if(op===3n)return String(f.divisor_of_poles());if(op===4n)return f.zeros().map(String).join(';');if(op===5n)return f.poles().map(String).join(';');
 if(op===6n||op===7n){const P=op===6n?K.maximal_order().ideal(K.gen().sub(K.__call__(a))).place():K.maximal_order_infinite().prime_ideal().place(),r=f.evaluate(P)as any;return String(r)+'|'+bool((r instanceof Rational?QQ:typeof r.parent==='function'?r.parent():r.parent)===k);}
 if(op===8n){const M=f.matrix();return join([Object.isFrozen(M)&&Object.isFrozen(M[0]),M.length,M[0]!.length]);}
 if(op===9n){const M=f.matrix();M[0]![0]=K.__call__(a);return parts(M[0]![0]!);}
 if(op===10n){const M=f.matrix(),C=M.map(row=>row.slice());C[0]![0]=K.__call__(a);return join([parts(M[0]![0]!),parts(C[0]![0]!),parts(f),Object.isFrozen(C)]);}
 if(op===11n){const before=parts(f),M=f.matrix(),t=f.trace(),n=f.norm();return [before,parts(f),parts(M[0]![0]!),parts(t),parts(n)].join('|');}
 throw new Error('unknown ff_element_observe operation');
}

function ff_predicate(p:bigint,op:bigint,kind:bigint):string{
 const K=field(p),counts=[0,0],fail=(name:string,message:string):never=>{const e=new Error(message);e.name=name;throw e;};
 class CategoryProbe{}
 class Probe{
  parent(){counts[0]++;if(kind===21n)return fail('ValueError','parent unavailable');return kind===17n||kind===19n&&counts[0]!>1?K:kind===18n||kind===19n?QQ:this;}
  category():any{counts[1]++;if(kind===25n)return fail('ValueError','category unavailable');if(kind===26n)return fail('AttributeError','category unavailable');if(kind===27n||kind===40n)return new CategoryProbe();if(kind===41n)return 0;if(kind===28n)return null;return {is_subcategory:(c:unknown)=>kind===23n&&String(c)==='Category of function fields'};}
 }
 const vals:any[]=[K.gen(),1n,QQ.__call__(1n).div(QQ.__call__(2n)),K._ring.gen(),K.field().gen(),K,K.maximal_order(),K.maximal_order().ideal(K.gen()),K.maximal_order().ideal(K.gen()).place(),QQ,GF(3n),[1n],{},null,true,1.5,'x'];
 let x:any;if(kind<BigInt(vals.length))x=vals[Number(kind)];else{x=new Probe();if(kind===20n)delete (Probe.prototype as any).parent;if(kind===29n)x.category=1;if(kind===30n)x.parent=1;if(kind===31n)x=()=>1n;if(kind===32n)x.category=1n;if(kind===33n)x.parent=1n;if(kind===34n)x.parent={};if(kind===35n)x.parent=[];if(kind===36n)x.parent=null;if(kind===37n)x.parent=false;if(kind===38n)x.parent='x';if(kind===39n)x.category=null;if(kind===40n)(CategoryProbe.prototype as any).is_subcategory=1n;if(kind===42n)x=new Integer(1n);if(kind===43n)x=GF(p===0n?3n:p).__call__(1n);if(kind===44n)x=GF2.one();if(kind===45n)x=new Matrix(QQ,1,1,[[QQ.one()]]);}
 return bool(op===0n?is_FunctionFieldElement(x):is_FunctionField(x))+'|'+counts.join(',');
}

function ff_predicate_access(p:bigint,op:bigint,kind:bigint):string{
 const K=field(p),counts=[0,0,0,0,0,0],fail=(name:string,message:string):never=>{const e=new Error(message);e.name=name;throw e;};let x:any,c:any;
 const parentGet=()=>{const n=++counts[0]!;if(kind===2n)return fail('AttributeError','parent unavailable');return ()=>{counts[1]++;return kind===7n?(n===1?QQ:K):x;};};
 const categoryGet=()=>{counts[2]++;if(kind===3n)return fail('AttributeError','category unavailable');if(kind===4n)return fail('ValueError','category unavailable');if(kind===6n)return null;return ()=>{counts[3]++;return c;};};
 const subGet=()=>{counts[4]++;if(kind===5n)return fail('AttributeError','subcategory unavailable');return (C:unknown)=>{counts[5]++;return String(C)==='Category of function fields';};};
 class Probe{get parent(){return parentGet();}get category(){return categoryGet();}}
 class CategoryProbe{get is_subcategory(){return subGet();}}
 class DynamicProbe{}class DynamicCategory{}
 x=kind===1n?new Proxy(new DynamicProbe(),{get:(target,key,receiver)=>key==='parent'?parentGet():key==='category'?categoryGet():Reflect.get(target,key,receiver)}):new Probe();
 c=kind===1n?new Proxy(new DynamicCategory(),{get:(target,key,receiver)=>key==='is_subcategory'?subGet():Reflect.get(target,key,receiver)}):new CategoryProbe();
 return bool(op===0n?is_FunctionFieldElement(x):is_FunctionField(x))+'|'+counts.join(',');
}

function ff_matrix_base(p:bigint,raw:bigint,warm:bigint,kind:bigint,num:bigint[],den:bigint[]):string{
 // Use a fresh parent so each scenario starts cold, like clear_cache() in Sage.
 const k:any=p===0n?QQ:GF(p),K=new RationalFunctionField(k,'x'),F=K.field(),R=K._ring;
 const f=raw?new FunctionFieldElement_rational(K,new F._element_class(F,R.__call__(num.map(x=>k.__call__(x))),R.__call__(den.map(x=>k.__call__(x))),{reduce:false})):element(K,num).div(element(K,den));
 class MatrixBaseProbe{}
 const bases:any[]=[null,K,QQ,FunctionField(k,'y'),F,R,false,0n,[],{},'x',K.gen(),K.maximal_order(),ZZ,GF(3n),new Set(),new Map(),new MatrixBaseProbe(),Object.create(null),()=>null],base=bases[Number(kind)];
 const observe=(g:any,b:any)=>{try{return parts(g.matrix(b)[0][0]);}catch(e:any){return e.name+':'+e.message;}};
 if(warm===1n)f.matrix();if(warm===2n)f.trace();if(warm===3n)f.norm();if(warm===4n)K.one().matrix();if(warm===5n)observe(f,QQ);if(warm===6n)FunctionField(k,'z').gen().matrix();
 return [observe(f,base),observe(f,base),observe(f,K),observe(K.one(),base),parts(f)].join('|');
}

function ff_variable_map(p:bigint,same:bigint,op:bigint,raw:bigint,kind:bigint,num:bigint[],den:bigint[]):string{
 const K=field(p),[L,fromL,toL]=K.change_variable_name(same?'x':'y'),source=op===0n?L:K,target=op===0n?K:L,map:any=op===0n?fromL:toL,R=source._ring,F=source.field(),k=source.constant_base_field();
 const f=raw?new FunctionFieldElement_rational(source,new F._element_class(F,R.__call__(num.map(x=>k.__call__(x))),R.__call__(den.map(x=>k.__call__(x))),{reduce:false})):element(source as any,num).div(element(source as any,den));
 const vals:any[]=[f,1n,QQ.__call__(1n).div(QQ.__call__(2n)),k.__call__(1n),R.__call__(num.map(x=>k.__call__(x))),target._ring.__call__(num.map(x=>k.__call__(x))),f.element(),num,null,'1','badname',{},true,1.5,target.gen(),{2:3n,0:1n},new Map([['bad',1n]]),[],[[1n,2n]],new Integer(1n),...["a'b",'a\"b',"a'b\"c",'a\nb','a\\b'].map(key=>new Map([[key,1n]]))],x=vals[Number(kind)],result=map(x);
 return [parts(result),bool(result.parent===target),bool(result===x),parts(f)].join('|');
}

function ff_poly_dict(p:bigint,keykind:bigint,valuekind:bigint,degrees:bigint[],values:bigint[]):string{
 const k:any=p===-1n?{zero:()=>new Integer(0n),one:()=>new Integer(1n),__call__:(x:unknown)=>new Integer(ZZ.__call__(x as bigint)),is_field:()=>false,characteristic:()=>0n,toString:()=> 'Integer Ring'}:p===0n?QQ:GF(p),R=new PolynomialRing<any>(k,'x');
 const value=(v:bigint):any=>valuekind===1n?[v,1n]:valuekind===2n?R.__call__([v,1n]):valuekind===3n?QQ.__call__(v).div(QQ.__call__(2n)):valuekind===4n?String(v):valuekind===5n?'x':valuekind===6n?'1/2':valuekind===7n?null:valuekind===8n?true:valuekind===9n?Number(v)+0.5:valuekind===10n?[]:valuekind===11n?[[]]:valuekind===12n?R.__call__(v):valuekind===13n?R.zero():v;
 const entries=degrees.map((e,i)=>[keykind===2n?[e]:keykind===3n?String(e):keykind===4n?Number(e)+0.5:keykind===5n?[]:keykind===6n?[e,7n]:keykind===7n?Number(e):e,value(values[i]!) ]as [any,any]);
 const d=keykind===1n?Object.fromEntries(entries):new Map(entries);return R.__call__(d).coeffs.map(String).join(',');
}

function ff_poly_dict_order(p:bigint,kind:bigint):string{
 const k:any=p===-1n?{zero:()=>new Integer(0n),one:()=>new Integer(1n),__call__:(x:unknown)=>new Integer(ZZ.__call__(x as bigint)),is_field:()=>false,characteristic:()=>0n,toString:()=> 'Integer Ring'}:p===0n?QQ:GF(p),R=new PolynomialRing<any>(k,'x');
 let items:any[]=[[0n,[1n,2n]],['bad',1n]];if(kind===1n)items.reverse();if(kind===2n)items=[[0n,'bad'],[2n,1n]];if(kind===3n)items=[[0n,[1n,2n]],[[2n],1n]];
 return R.__call__(new Map(items)).coeffs.map(String).join(',');
}

function ff_flint_rng(op:bigint,seed1:bigint,seed2:bigint,limit:bigint,count:bigint):string{
 const state=flint_rand_init();if(seed1>=0n)flint_rand_set_seed(state,seed1,seed2);const out:bigint[]=[];
 for(let i=0n;i<count;i++)out.push(op===0n?n_randlimb(state):op===1n?n_randint(state,limit):op===2n?n_urandint(state,limit):op===3n?n_randbits(state,Number(limit)):op===4n?n_randtest_bits(state,Number(limit)):op===5n?n_randtest(state):n_randtest_not_zero(state));
 flint_rand_clear(state);return [...out,...flint_rand_get_seed(state)].join(',');
}

function ff_flint_factor(op:bigint,p:bigint,d:bigint,a:bigint[],b:bigint[],c:bigint[],seed1:bigint,seed2:bigint):string{
 const state=flint_rand_init();if(seed1>=0n)flint_rand_set_seed(state,seed1,seed2);
 let scalar=0n, chunks:Array<[bigint[],number]>=[];
 if(op===0n)scalar=BigInt(nmod_poly_deflation(a,p));
 else if(op===1n)chunks=[[nmod_poly_deflate(a,Number(d),p),1]];
 else if(op===2n)chunks=[[nmod_poly_inflate(a,Number(d),p),1]];
 else if(op===3n){const [v,e]=nmod_poly_remove(a,b,p);scalar=BigInt(e);chunks=[[v,1]];}
 else if(op===4n)chunks=[[_nmod_vec_rand(state,Number(d),p),1]];
 else if(op===5n)chunks=[[_nmod_vec_randtest(state,Number(d),p),1]];
 else if(op===6n)chunks=[[nmod_poly_randtest(state,Number(d),p),1]];
 else if(op===7n){const f=nmod_poly_factor_equal_deg_prob(state,a,Number(d),p);if(f){scalar=1n;chunks=[[f,1]];}}
 else if(op===8n)chunks=nmod_poly_factor_equal_deg(a,Number(d),p);
 else if(op===9n)chunks=nmod_poly_factor_distinct_deg(a,p);
 else if(op===11n)chunks=nmod_poly_factor_kaltofen_shoup(a,p);
 else if(op===12n)[scalar,chunks]=nmod_poly_factor(a,p);
 else if(op===13n)[scalar,chunks]=nmod_poly_factor_with_cantor_zassenhaus(a,p);
 else if(op===14n)[scalar,chunks]=nmod_poly_factor_with_kaltofen_shoup(a,p);
 else if(op===10n)chunks=nmod_poly_factor_cantor_zassenhaus(a,p);
 else if(op>=15n&&op<=17n){
  const modulus=op===15n?b:c,inv=_nmod_poly_inv_series_newton(modulus.slice().reverse(),modulus.length,p);
  if(op===15n)chunks=nmod_poly_precompute_matrix(a,b,inv,p).map(v=>[v,1]);
  else if(op===16n)chunks=[[nmod_poly_compose_mod_brent_kung_precomp_preinv(a,nmod_poly_precompute_matrix(b,c,inv,p),c,inv,p),1]];
  else {const polys=Array.from({length:Number(d)},(_,i)=>a.slice(i*(c.length-1),(i+1)*(c.length-1)));chunks=nmod_poly_compose_mod_brent_kung_vec_preinv(polys,Number(d),b,c,inv,p).map(v=>[v,1]);}
 }
 else if(op===18n){const left=Array.from({length:Number(seed1)},(_,i)=>a.slice(i*Number(d),(i+1)*Number(d))),right=Array.from({length:Number(d)},(_,i)=>b.slice(i*Number(seed2),(i+1)*Number(seed2)));chunks=nmod_mat_mul(left,right,p,Number(seed2)).map(v=>[v,1]);}
 else if(op===19n){const rows=Array.from({length:Number(d)},(_,i)=>a.slice(i*(c.length-1),(i+1)*(c.length-1))),inv=_nmod_poly_inv_series_newton(c.slice().reverse(),c.length,p);chunks=[[_nmod_poly_mod_matrix_rows_evaluate(rows,b,c,inv,p),1]];}
 else if(op===20n){const inv=_nmod_poly_inv_series_newton(b.slice().reverse(),b.length,p);chunks=_nmod_poly_reduce_matrix_mod_poly(nmod_poly_precompute_matrix(a,b,inv,p),c,p).map(v=>[v,1]);}
 else throw new Error('unknown ff_flint_factor operation');
 return [scalar,...flint_rand_get_seed(state),chunks.length,...chunks.flatMap(([f,e])=>[e,f.length,...f])].join(',');
}

function ff_word_factor(p:bigint,a:bigint[]):string{
 const R=new PolynomialRing(GF(p),'x');return R.__call__(a).factor().map(([f,e])=>f.coeffs.map(String).join(',')+':'+e).join(';');
}
function ff_flint_matrix(p:bigint,rows:bigint,inner:bigint,cols:bigint):string{
 const a=Array.from({length:Number(rows*inner)},(_,i)=>(BigInt(i)*17n+BigInt(Math.floor(i/Number(inner)))*13n+3n)%p);
 const b=Array.from({length:Number(inner*cols)},(_,i)=>(BigInt(i)*23n+BigInt(Math.floor(i/Number(cols)))*7n+5n)%p);
 return ff_flint_factor(18n,p,inner,a,b,[],rows,cols);
}

function ff_pari_rng(op:bigint,seed:bigint,limit:bigint,count:bigint,length:bigint):string{
 pari_init_rand();if(seed!==-1n)setrand(seed);if(op===9n)pari_init_rand();
 const out:string[]=[];
 for(let i=0n;i<count;i++){
  if(op===0n||op===9n)out.push(String(pari_rand()));
  else if(op===1n)out.push(String(random_bits(Number(limit))));
  else if(op===2n)out.push(String(random_Fl(limit)));
  else if(op===3n)out.push(String(randomi(limit)));
  else if(op===4n)out.push(String(random_F2x(Number(limit))));
  else if(op===5n)out.push(random_zv(Number(limit)).join(','));
  else if(op===6n)out.push(random_Flx(Number(length),limit).join(','));
  else if(op===7n)out.push(random_FpX(Number(length),limit).join(','));
  else if(op===8n){const old=getrand(),a=pari_rand();setrand(old);out.push(a+','+pari_rand());}
 }
 return out.join(';')+'|'+getrand();
}

function ff_pari_binary(op:bigint,a:bigint,b:bigint,c:bigint,d:bigint,e:bigint):string{
 if(op===0n)return String(F2x_mul(a,b));if(op===1n)return String(F2x_sqr(a));if(op===2n)return String(F2x_sqrt(a));
 if(op===3n)return String(F2x_gcd(a,b));if(op===4n)return F2x_divrem(a,b).join(',');if(op===5n)return String(F2x_rem(a,b));
 if(op===6n)return String(F2x_deriv(a));if(op===7n)return String(F2x_Frobenius(a));if(op===8n)return F2x_matFrobenius(a).join(',');
 if(op===9n){
  if(d<0n||d>BigInt(Number.MAX_SAFE_INTEGER))return String(F2m_ker_sp([],Number(d),Number(e)));
  const mask=(1n<<d)-1n,columns=Array.from({length:Number(b)},(_,i)=>i===Number(b)-1?a>>(BigInt(i)*d):(a>>(BigInt(i)*d))&mask);
  const result=e===2n?F2m_ker(columns,Number(d)):F2m_ker_sp(columns,Number(d),Number(e));
  return (result===null?'null':Array.isArray(result)?'m:'+result.join(','):'v:'+result)+'|'+columns.join(',')+'|'+Number(result===columns);
 }
 if(op===10n)return F2xq_powers(a,Number(d),b).join(',');if(op===11n)return String(F2x_degree(a));if(op===12n)return F2x_valrem(a).join(',');
 if(op===13n)return String(F2xq_mul(a,b,c));if(op===14n)return String(F2xq_sqr(a,b));if(op===15n)return String(F2x_add(a,b));
 if(op>=16n&&op<=18n){setrand(b);const result=op===16n?F2x_factor(a).flat():op===17n?F2x_factor_squarefree(a):F2x_ddf(a).flat();return result.join(',')+'|'+getrand();}
 throw new Error('unknown binary kernel operation');
}

function ff_binary_factor(bits:bigint,kind:bigint):string{
 const k:any=kind===0n?GF(2n):kind===1n?GF2:Zmod(2n),R=new PolynomialRing<any>(k,'x');
 const coeffs=Array.from({length:bits?bits.toString(2).length:0},(_,i)=>k.__call__((bits>>BigInt(i))&1n));
 return R.__call__(coeffs).factor().map(([f,e])=>f.coeffs.map(String).join(',')+':'+e).join(';');
}

function ff_pari_quotient(op:bigint,p:bigint,n:bigint,x:bigint[],T:bigint[]):string{
 const r=op===0n?FpXQ_powers(x,Number(n),T,p):op===1n?FpXQ_autpow(x,Number(n),T,p):FpXQ_autpowers(x,Number(n),T,p).slice(1);
 const poly=(f:bigint[])=>'['+f.join(',')+']';return op===1n?poly(r as bigint[]):'['+(r as bigint[][]).map(poly).join(',')+']';
}

function ff_pari_matrix(op:bigint,p:bigint,m:bigint,n:bigint,k:bigint,ab:bigint,bb:bigint,seed:bigint):string{
 const entry=(i:bigint,j:bigint,bits:bigint,s:bigint):bigint=>{
  if(!bits||!s||(s%4n===3n&&(i+3n*j)%4n!==0n))return 0n;
  const v=((1n<<(bits-1n))+(i*17n+j*23n+s*31n)%257n)%(1n<<bits);return s%2n===0n&&(i+j)%2n!==0n?-v:v;
 };
 const matrix=(rows:bigint,cols:bigint,bits:bigint,s:bigint):bigint[][]=>[[],...Array.from({length:Number(cols)},(_,j)=>[0n,...Array.from({length:Number(rows)},(_,i)=>entry(BigInt(i+1),BigInt(j+1),bits,s))])];
 const a=matrix(m,n,ab,seed),b=matrix(n,k,bb,seed?seed+5n:0n);
 const red=(v:bigint[][])=>[[],...v.slice(1).map(c=>[0n,...c.slice(1).map(z=>((z%p)+p)%p)])];
 let r:bigint[][];
 if(op===0n)r=pari_ZM_mul(a,b);
 else if(op===1n)r=Flm_mul(red(a),red(b),p);
 else if(op===2n)r=FpM_mul(a,b,p);
 else {
  const pack=(v:bigint[][])=>v.slice(1).map(c=>c.slice(1).reduce((z,v,i)=>z|((v&1n)<<BigInt(i)),0n));
  const rows=n?m:0n;r=[[],...F2m_mul(pack(a),pack(b),Number(m)).map(v=>[0n,...Array.from({length:Number(rows)},(_,i)=>(v>>BigInt(i))&1n)])];
 }
 return '['+r.slice(1).map(c=>'['+c.slice(1).join(',')+']').join(',')+']';
}

function ff_pari_matrix_guard(op:bigint,p:bigint,rows:bigint,x:bigint[],y:bigint[]):boolean{
 if(op===0n)Flm_mul([[]],[[]],p);else if(op===1n)FpM_mul([[]],[[]],p);else F2m_mul(x,y,Number(rows));return true;
}

function ff_pari_composition(op:bigint,p:bigint,n:bigint,Q:bigint[],x:bigint[],T:bigint[]):string{
 if(op===8n||op===9n){const r=op===8n?pari_FpXQV_eval(Q,[],T,p):Flx_FlxqV_eval(Q,[],T,p);return '['+r.join(',')+']';}
 if(op===4n)return String(brent_kung_optpow(Number(p),Number(n),Number(Q[0]??0n)));
 if(op===5n)return '['+FpXQ_auttrace([Q,x],n,T,p).map(f=>'['+f.join(',')+']').join(',')+']';
 const r=op===0n?pari_FpX_eval(Q,x,T,p):op===1n?pari_FpXQV_eval(Q,FpXQ_powers(x,Number(n),T,p),T,p):op===2n?Flx_Flxq_eval(Q,x,T,p):op===3n?Flx_FlxqV_eval(Q,Flxq_powers(x,Number(n),T,p),T,p):op===6n?pari_FpX_mul(Q,x,p):pari_FpX_scalar(Q,x[0]??0n,p);
 return '['+r.join(',')+']';
}

function ff_pari_polymul(op:bigint,p:bigint,a:bigint[],b:bigint[]):string{
 const r=op===0n?pari_FpX_mul(a,b,p):op===1n?FpX_sqr(a,p):op===2n?pari_ZX_mul(a,b):op===3n?ZX_sqr(a):op===4n?Flx_mul(a,b,p):op===5n?Flx_sqr(a,p):op===6n?pari_ZX_mul(a,a):pari_FpX_mul(a,a,p);
 return '['+r.join(',')+']';
}

function ff_pari_polymul_pattern(op:bigint,p:bigint,n:bigint,m:bigint,ea:bigint,eb:bigint,kind:bigint):string{
 const poly=(length:bigint,e:bigint,salt:bigint):bigint[]=>Array.from({length:Number(length)},(_,i)=>{
  let c=kind===0n?0n:(1n<<e)+(BigInt(i)*17n+salt*31n)%(1n<<e);
  if(kind===2n&&i%2)c=-c;if(kind===3n&&i%5===0)c=0n;
  if(op===4n||op===5n)c=kind===4n?p-1n:((c%p)+p)%p;return c;
 });
 return ff_pari_polymul(op,p,poly(n,ea,1n),poly(m,eb,2n));
}

function ff_pari_polydiv(op:bigint,p:bigint,a:bigint[],b:bigint[]):string{
 const r=op===0n?pari_FpX_divrem(a,b,p):op===1n?pari_FpX_rem(a,b,p):op===2n?FpX_invBarrett(a,p):op===3n?Flx_divrem(a,b,p):op===4n?Flx_rem(a,b,p):Flx_invBarrett(a,p);
 return JSON.stringify(r,(_,v)=>typeof v==='bigint'?v.toString():v).replace(/"/g,'');
}

function ff_pari_polydiv_pattern(op:bigint,p:bigint,n:bigint,m:bigint,kind:bigint):string{
 const poly=(length:bigint,divisor:boolean):bigint[]=>Array.from({length:Number(length)},(_,i)=>{
  if(divisor&&i===Number(length)-1)return kind===1n&&p>2n?p-1n:1n;
  let c=(BigInt(i)*BigInt(i)*17n+BigInt(i)*(divisor?13n:31n)+7n)%p;
  if(kind===0n)c=0n;
  if(kind===2n&&op<3n)c=i%2?c-p:c+p;
  if(kind===3n&&i%5!==0)c=0n;
  return c;
 });
 return ff_pari_polydiv(op,p,poly(n,op===2n||op===5n),poly(m,true));
}

function ff_pari_polygcd(op:bigint,p:bigint,a:bigint[],b:bigint[]):string{
 if(op===15n){let squares=0,multiplies=0;const value=gen_pow_i(a[0]??0n,b[0]??0n,x=>{squares++;return (3n*x+1n)%p;},(x,y)=>{multiplies++;return (5n*x+7n*y+1n)%p;});return '['+[value,squares,multiplies].join(',')+']';}
 const r=op===14n?pari_gcd_power(a,b[0]??0n,b.slice(1),p):op===11n?FpX_split_part(a,p):op===12n?FpX_nbroots(a,p):op===0n?pari_FpX_gcd(a,b,p):op===1n?pari_FpX_extgcd(a,b,p):op===2n?FpX_halfgcd_all(a,b,p):op===3n?Flx_gcd(a,b,p):op===4n?Flx_extgcd(a,b,p):op===5n?Flx_halfgcd_all(a,b,p):op===9n?FpX_halfgcd(a,b,p):op===10n?Flx_halfgcd(a,b,p):op===6n?pari_gcd_inverse(a,b,p):op===7n?ZpX_liftfact(a,[[],b,pari_FpX_divrem(a,b,p)[0]],p,3).slice(1):bezout_lift_fact(a,[[],b,pari_FpX_divrem(a,b,p)[0]],p,3).slice(1);
 return JSON.stringify(r,(_,v)=>typeof v==='bigint'?v.toString():v).replace(/"/g,'');
}

function ff_pari_polygcd_pattern(op:bigint,p:bigint,n:bigint,m:bigint,kind:bigint):string{
 const poly=(length:bigint,salt:bigint):bigint[]=>{
  let state=salt;const mask=(1n<<128n)-1n,out:bigint[]=[];
  for(let i=0;i<Number(length);i++){
   state=(state*6364136223846793005n+1442695040888963407n)&mask;
   let c=state%p;if(kind===0n||kind===3n&&i%5!==0)c=0n;out.push(c);
  }
  if(kind===2n&&out.length){const v=new Array<bigint>(out.length+2).fill(0n);for(let i=0;i<out.length;i++){v[i]=(v[i]!+out[i]!)%p;v[i+2]=(v[i+2]!+out[i]!)%p;}return v;}
  return out;
 };
 return ff_pari_polygcd(op,p,poly(n,17n),poly(m,43n));
}

import { FpXQ_minpoly as pari_minpoly } from '../../../../packages/parigp-ts/src/galconj.js';
function ff_pari_minpoly(op:bigint,seed:bigint,p:bigint,a:bigint[],b:bigint[]):string{
 const saved=getrand();try{setrand(seed);const r=op===1n?pari_word_minpoly(a,b,p):pari_minpoly(a,b,p);return '['+r.join(',')+']|'+getrand();}finally{setrand(saved);}
}

import { Flxq_minpoly as pari_word_minpoly } from '../../../../packages/parigp-ts/src/Flx.js';

import { FpX_ddf as pari_ddf, FpX_nbfact_by_degree as pari_ddf_count } from '../../../../packages/parigp-ts/src/galconj.js';
function ff_pari_ddf(op:bigint,p:bigint,a:bigint[]):string{
 const counts=()=>{const {D,nb}=pari_ddf_count(a,p);return [D,nb];};
 const wordCounts=()=>{const {D,nb}=pari_word_ddf_count(a,p);return [D,nb];};
 const r=op===4n?wordCounts():op===3n||op===5n?counts():op===2n?pari_nbfact(a,p):op===1n?pari_word_ddf(a,p).map(([f,d])=>[d,f]):[...pari_ddf(a,p)];
 return JSON.stringify(r,(_,v)=>typeof v==='bigint'?v.toString():v).replaceAll('"','');
}

import { Flx_ddf as pari_word_ddf, FpX_nbfact as pari_nbfact } from '../../../../packages/parigp-ts/src/FpX_factor.js';

import { Flx_nbfact_by_degree as pari_word_ddf_count } from '../../../../packages/parigp-ts/src/FpX_factor.js';

import { Fp_sqrt as native_pari_sqrt, Fp_sqrt_i, Fl_sqrt } from '../../../../packages/parigp-ts/src/ff.js';
function ff_pari_sqrt(op:bigint,p:bigint,a:bigint,y:bigint):string{
 if(op===3n){let squares=0,fused=0;const mod=(x:bigint)=>((x%p)+p)%p;const r=gen_pow_fold(a,y,x=>{squares++;return mod(3n*x+1n);},x=>{fused++;return mod(7n*x+2n);});return `[${r},${squares},${fused}]`;}
 const r=op===1n?Fp_sqrt_i(a,y===0n?null:y,p):op===2n?Fl_sqrt(a,p):native_pari_sqrt(a,p);return r===null?'null':String(r);
}

import { Flx_normalize as pariNormalizeWord } from '../../../../packages/parigp-ts/src/Flx.js';
import { FpX_normalize as pariNormalizePolynomial } from '../../../../packages/parigp-ts/src/ffinit.js';
function ff_pari_normalize(op:bigint,p:bigint,a:bigint[]):string{
 return '['+(op===1n?pariNormalizeWord(a,p):pariNormalizePolynomial(a,p)).join(',')+']';
}

import { _galconj_factor_squarefree_irreducibles as pariFactorSquarefree } from '../../../../packages/parigp-ts/src/galconj.js';
function ff_pari_factor(op:bigint,seed:bigint,p:bigint,a:bigint[]):string{
 const saved=getrand();
 try{
  setrand(seed);
  try{
   const r=op===3n?FpX_factor_squarefree(a,p):op===2n?pariFactorSquarefree(a,p).map(f=>[f,1]):op===1n?pariWordFactor(a,p):pariFactor(a,p);
   return JSON.stringify(r,(_,v)=>typeof v==='bigint'?v.toString():v).replaceAll('"','')+'|'+getrand();
  }catch(error){
   if(op===4n)return 'ERROR '+(error as Error).message+'|'+getrand();
   throw error;
  }
 }finally{setrand(saved);}
}

import { FpX_factor as pariFactor, Flx_factor as pariWordFactor } from '../../../../packages/parigp-ts/src/FpX_factor.js';


import { FpX_roots as pariRoots, FpX_is_totally_split as pariTotallySplit, FpX_nbroots as pariRootCount } from '../../../../packages/parigp-ts/src/galconj.js';
function ff_pari_roots(op:bigint,p:bigint,a:bigint[]):string{
 const r=op===0n?pariRoots(a,p):op===1n?pariWordRoots(a,p):op===2n?Number(pariTotallySplit(a,p)):op===4n?Number(pariWordTotallySplit(a,p)):op===5n?pariWordRootCount(a,p):pariRootCount(a,p);
 return JSON.stringify(r,(_,v)=>typeof v==='bigint'?v.toString():v).replaceAll('"','');
}

import {Flx_roots as pariWordRoots, Flx_nbroots as pariWordRootCount, Flx_is_totally_split as pariWordTotallySplit} from '../../../../packages/parigp-ts/src/FpX_factor.js';


import {FpXQX_mul,FpXQX_sqr,FpXQX_red} from '../../../../packages/parigp-ts/src/FpXX.js';
import {FpXQX_normalize} from '../../../../packages/parigp-ts/src/polarit3.js';
import {FlxqX_mul,FlxqX_sqr,FlxqX_red,FlxqX_normalize} from '../../../../packages/parigp-ts/src/FlxX.js';
import {F2xqX_mul,F2xqX_sqr,F2xqX_red,F2xqX_normalize} from '../../../../packages/parigp-ts/src/F2x.js';
function ff_pari_extension(mode:bigint,op:bigint,p:bigint,T:bigint[],alen:bigint[],aflat:bigint[],blen:bigint[],bflat:bigint[]):string{
 const decode=(lengths:bigint[],flat:bigint[]):(bigint|bigint[])[]=>{let offset=0;return lengths.map(l=>{const n=Number(l);if(n<0)return flat[offset++]!;const c=flat.slice(offset,offset+n);offset+=n;return c;});};
 const a=decode(alen,aflat),b=decode(blen,bflat);let r;
 if(mode===0n)r=op===0n?FpXQX_mul(a,b,T,p):op===1n?FpXQX_sqr(a,T,p):op===2n?FpXQX_red(a,T,p):FpXQX_normalize(a,T,p);
 else if(mode===1n){const x=a as bigint[][],y=b as bigint[][];r=op===0n?FlxqX_mul(x,y,T,p):op===1n?FlxqX_sqr(x,T,p):op===2n?FlxqX_red(x,T,p):FlxqX_normalize(x,T,p);}
 else{const x=a as bigint[],y=b as bigint[];r=op===0n?F2xqX_mul(x,y,T[0]!):op===1n?F2xqX_sqr(x,T[0]!):op===2n?F2xqX_red(x,T[0]!):F2xqX_normalize(x,T[0]!);}
 return JSON.stringify(r,(_,v)=>typeof v==='bigint'?v.toString():v).replaceAll('"','');
}

import {createHash as extensionHash} from 'node:crypto';
import {nodeWordExtensionProduct} from '../pari_extension_node.js';
function ff_pari_extension_large(degree:bigint,zero:bigint,host:bigint):string{
 const d=Number(degree);
 if(host)return nodeWordExtensionProduct(d,Boolean(zero));
 const result=FlxqX_mul([new Array<bigint>(d).fill(1n)],zero?[]:[[1n]],[1n,...new Array<bigint>(d-1).fill(0n),1n],7n);
 const text=JSON.stringify(result,(_,v)=>typeof v==='bigint'?v.toString():v).replaceAll('"','');
 return extensionHash('sha256').update(text).digest('hex');
}

import * as pariExtensionGeneric from '../../../../packages/parigp-ts/src/FpXX.js';
import * as pariExtensionWord from '../../../../packages/parigp-ts/src/FlxX.js';
import * as pariExtensionBinary from '../../../../packages/parigp-ts/src/F2x.js';
import {extensionDivision as pariExtensionDivision, type ExtensionModulus} from '../../../../packages/parigp-ts/src/_extension_division.js';
import type {ExtensionPolynomial} from '../../../../packages/parigp-ts/src/_extension_polynomial.js';
function ff_pari_extension_division(mode:bigint,operation:bigint,p:bigint,T:bigint[],alen:bigint[],aflat:bigint[],blen:bigint[],bflat:bigint[]):string{
 const decode=(lengths:bigint[],flat:bigint[]):ExtensionPolynomial=>{let offset=0;return lengths.map(l=>{const n=Number(l);if(n<0)return flat[offset++]!;const c=flat.slice(offset,offset+n);offset+=n;return c;});};
 const a=decode(alen,aflat),b=decode(blen,bflat),op=Number(operation%10n),cache=Number(operation/10n);
 let r;
 if(mode===0n){
  const api=pariExtensionGeneric;
  const S:ExtensionModulus=op===2||!cache?b:cache===1?{polynomial:b,inverse:api.FpXQX_invBarrett(b,T,p)}:api.FpXQX_get_red(b,T,p);
  r=op===0?api.FpXQX_divrem(a,S,T,p):op===1?api.FpXQX_rem(a,S,T,p):op===2?api.FpXQX_invBarrett(a,T,p):op===3?pariExtensionDivision(0,3,p,T,a,S):api.FpXQX_div(a,S,T,p);
 }else if(mode===1n){
  const api=pariExtensionWord,x=a as bigint[][],y=b as bigint[][];
  const S:ExtensionModulus<bigint[][]>=op===2||!cache?y:cache===1?{polynomial:y,inverse:api.FlxqX_invBarrett(y,T,p)}:api.FlxqX_get_red(y,T,p);
  r=op===0?api.FlxqX_divrem(x,S,T,p):op===1?api.FlxqX_rem(x,S,T,p):op===2?api.FlxqX_invBarrett(x,T,p):op===3?pariExtensionDivision(1,3,p,T,x,S):api.FlxqX_div(x,S,T,p);
 }else{
  const api=pariExtensionBinary,x=a as bigint[],y=b as bigint[],t=T[0]!;
  const S:ExtensionModulus<bigint[]>=op===2||!cache?y:cache===1?{polynomial:y,inverse:api.F2xqX_invBarrett(y,t)}:api.F2xqX_get_red(y,t);
  r=op===0?api.F2xqX_divrem(x,S,t):op===1?api.F2xqX_rem(x,S,t):op===2?api.F2xqX_invBarrett(x,t):op===3?pariExtensionDivision(2,3,2n,t,x,S):api.F2xqX_div(x,S,t);
 }
 return JSON.stringify(r,(_,v)=>typeof v==='bigint'?v.toString():v).replaceAll('"','');
}

function ff_pari_extension_division_pattern(mode:bigint,op:bigint,degree:bigint):string{
 const m=Number(mode) as 0|1|2,n=Number(degree),p=m===2?2n:m===1?17n:(1n<<64n)+13n,T=m===2?7n:[3n,0n,1n];
 const a:ExtensionPolynomial=m===2?Array.from({length:n+1},(_,i)=>BigInt(i%4)):Array.from({length:n+1},(_,i)=>[BigInt(i%17),BigInt(Math.floor(i/17)%17)]);
 a[n]=m===2?1n:[1n];const b=m===2?[3n,1n,1n]:[[1n,1n],[1n],[1n]];
 const result=pariExtensionDivision(m,Number(op) as 0|1,p,T,a,b);
 const text=JSON.stringify(result,(_,v)=>typeof v==='bigint'?v.toString():v).replaceAll('"','');
 return extensionHash('sha256').update(text).digest('hex');
}

import { extensionGcd } from '../../../../packages/parigp-ts/src/_extension_gcd.js';
function ff_pari_extension_gcd(mode:bigint,op:bigint,p:bigint,T:bigint[],alen:bigint[],aflat:bigint[],blen:bigint[],bflat:bigint[]):string{
 const decode=(lengths:bigint[],flat:bigint[])=>{
  let offset=0;return lengths.map(length=>{
   const n=Number(length);if(n<0)return flat[offset++]!;
   const c=flat.slice(offset,offset+n);offset+=n;return c;
  });
 };
 const m=Number(mode)as 0|1|2,o=Number(op)as 0|1|2|3,a=decode(alen,aflat),b=decode(blen,bflat);
 const compact=(v:any):string=>Array.isArray(v)?'['+v.map(compact).join(',')+']':String(v);
 const result=o===3?extensionGcd(m,o,p,m===2?T[0]!:T,a,b)
  :m===0?[pariExtensionGeneric.FpXQX_gcd,pariExtensionGeneric.FpXQX_extgcd,pariExtensionGeneric.FpXQX_halfgcd][o]!(a,b,T,p)
  :m===1?[pariExtensionWord.FlxqX_gcd,pariExtensionWord.FlxqX_extgcd,pariExtensionWord.FlxqX_halfgcd][o]!(a as bigint[][],b as bigint[][],T,p)
  :[pariExtensionBinary.F2xqX_gcd,pariExtensionBinary.F2xqX_extgcd,pariExtensionBinary.F2xqX_halfgcd][o]!(a as bigint[],b as bigint[],T[0]!);
 return compact(result);
}

import {binaryExtensionGcdTermination} from '../pari_extension_gcd_termination.js';
function ff_pari_extension_gcd_termination(degree:bigint,op:bigint):string{
 return binaryExtensionGcdTermination(Number(degree),Number(op));
}

import {extensionQuotient} from '../../../../packages/parigp-ts/src/_extension_quotient.js';
import {polynomialQuotient as extensionInnerQuotient} from '../../../../packages/parigp-ts/src/_polynomial_quotient.js';
import {extensionGetRed} from '../../../../packages/parigp-ts/src/_extension_division.js';
import {validateExtensionInputs as extensionInputGuard} from '../../../../packages/parigp-ts/src/_extension_polynomial.js';
function ff_pari_extension_quotient(mode:bigint,operation:bigint,p:bigint,n:bigint,T:bigint[],alen:bigint[],aflat:bigint[],blen:bigint[],bflat:bigint[],slen:bigint[],sflat:bigint[]):string{
 const decode=(lengths:bigint[],flat:bigint[]):ExtensionPolynomial=>{
  let offset=0;return lengths.map(length=>{
   const size=Number(length);if(size<0)return flat[offset++]!;
   const c=flat.slice(offset,offset+size);offset+=size;return c;
  });
 };
 const m=Number(mode)as 0|1|2,op=Number(operation%10n)as 0|1|2|3|4|5|6|7,cache=Number(operation%100n/10n);
 const a=decode(alen,aflat),b=decode(blen,bflat),raw=decode(slen,sflat);
 let t=m===2?T[0]!:T;
 extensionInputGuard(m,p,t,a,b);extensionInputGuard(m,p,t,raw);
 if(op===6&&(n<0n||n>=0xffffffffn))throw new RangeError('power count must be a nonnegative array length');
 if(op===7&&(n<0n||n>=1n<<64n))throw new RangeError('exponent must be an unsigned word integer');
 if(Array.isArray(t)){t=t.slice();while(t.length&&t.at(-1)===0n)t.pop();}
 const inner=operation>=100n&&m!==2?extensionInnerQuotient(t as bigint[],p,m===1).inverse:undefined;
 const S:ExtensionModulus=cache===1?{polynomial:raw,inverse:pariExtensionDivision(m,2,p,t,raw,[],inner)as ExtensionPolynomial}
  :cache===2?extensionGetRed(m,raw,t,p,inner):raw;
 let result;
 if(inner!==undefined || m===2&&op===4)result=extensionQuotient(m,op,p,n,t,a,b,S,inner);
 else{
  const api:any=m===0?pariExtensionGeneric:m===1?pariExtensionWord:pariExtensionBinary;
  const prefix=m===0?'FpXQXQ_':m===1?'FlxqXQ_':'F2xqXQ_';
  const tail=m===2?[S,t]:[S,t,p];
  const name=['mul','sqr','invsafe','inv','div','pow','powers','powu'][op]!;
  result=op===0||op===4?api[prefix+name](a,b,...tail)
   :op===5||op===7?api[prefix+name](a,n,...tail)
   :op===6?api[prefix+name](a,Number(n),...tail):api[prefix+name](a,...tail);
 }
 const compact=(v:any):string=>Array.isArray(v)?'['+v.map(compact).join(',')+']':String(v);
 return compact(result);
}

import {FlxqM_mul} from '../../../../packages/parigp-ts/src/alglin1.js';
import {wordExtensionMatrix} from '../../../../packages/parigp-ts/src/_extension_matrix.js';
import {extensionComposition} from '../../../../packages/parigp-ts/src/_extension_composition.js';
function extensionDecode(lengths:bigint[],flat:bigint[]):ExtensionPolynomial{
 let offset=0;return lengths.map(length=>{
  const n=Number(length);if(n<0)return flat[offset++]!;
  const v=flat.slice(offset,offset+n);offset+=n;return v;
 });
}
function extensionColumns(shape:bigint[],lengths:bigint[],flat:bigint[]):ExtensionPolynomial[]{
 const coefficients=extensionDecode(lengths,flat);let offset=0;
 return shape.map(length=>{const n=Number(length),v=coefficients.slice(offset,offset+n);offset+=n;return v;});
}
const extensionCompact=(v:any):string=>Array.isArray(v)?'['+v.map(extensionCompact).join(',')+']':String(v);
function ff_pari_extension_matrix(cache:bigint,p:bigint,T:bigint[],ashape:bigint[],alen:bigint[],aflat:bigint[],bshape:bigint[],blen:bigint[],bflat:bigint[]):string{
 const A=extensionColumns(ashape,alen,aflat)as bigint[][][],B=extensionColumns(bshape,blen,bflat)as bigint[][][];
 if(!cache){
  const sentinel=(M:bigint[][][])=>[[],...M.map(c=>[[],...c])];
  return extensionCompact(FlxqM_mul(sentinel(A),sentinel(B),T,p).slice(1).map(c=>c.slice(1)));
 }
 extensionInputGuard(1,p,T,[]);
 for(const c of [...A,...B])extensionInputGuard(1,p,T,c);
 if(A.some(c=>c.length!==A[0]!.length)||B.some(c=>c.length!==A.length))throw new RangeError('extension matrices must be rectangular and dimension-compatible');
 T=T.slice();while(T.length&&T.at(-1)===0n)T.pop();
 const inverse=extensionInnerQuotient(T,p,true).inverse;
 return extensionCompact(wordExtensionMatrix(A,B,T,p,inverse));
}
function ff_pari_extension_composition(mode:bigint,code:bigint,p:bigint,T:bigint[],qlen:bigint[],qflat:bigint[],xlen:bigint[],xflat:bigint[],vshape:bigint[],vlen:bigint[],vflat:bigint[],slen:bigint[],sflat:bigint[]):string{
 const m=Number(mode)as 0|1|2,op=Number(code%10n)as 0|1,cache=Number(code%100n/10n);
 const Q=extensionDecode(qlen,qflat),x=extensionDecode(xlen,xflat),V=extensionColumns(vshape,vlen,vflat),raw=extensionDecode(slen,sflat);
 let t=m===2?T[0]!:T;
 extensionInputGuard(m,p,t,Q,x);extensionInputGuard(m,p,t,raw);
 for(const v of V)extensionInputGuard(m,p,t,v);
 if(Array.isArray(t)){t=t.slice();while(t.length&&t.at(-1)===0n)t.pop();}
 const inverse=code>=100n&&m!==2?extensionInnerQuotient(t as bigint[],p,m===1).inverse:undefined;
 const S:ExtensionModulus=cache===1?{polynomial:raw,inverse:pariExtensionDivision(m,2,p,t,raw,[],inverse)as ExtensionPolynomial}
  :cache===2?extensionGetRed(m,raw,t,p,inverse):raw;
 if(inverse!==undefined)return extensionCompact(extensionComposition(m,op,p,t,Q,x,V,S,inverse));
 const api:any=m===0?pariExtensionGeneric:m===1?pariExtensionWord:pariExtensionBinary;
 const name=(m===0?'FpXQX_FpXQXQ':m===1?'FlxqX_FlxqXQ':'F2xqX_F2xqXQ')+(op?'V_eval':'_eval');
 return extensionCompact(m===2?api[name](Q,op?V:x,S,t):api[name](Q,op?V:x,S,t,p));
}

function ff_pari_inner_cache(op:bigint,p:bigint,n:bigint,Q:bigint[],x:bigint[],vlen:bigint[],vflat:bigint[],T:bigint[]):string{
 const V=extensionDecode(vlen,vflat)as bigint[][];
 return extensionCompact(op===0n?FpXQ_powers(x,Number(n),T,p):op===1n?Flxq_powers(x,Number(n),T,p)
  :op===2n?pari_FpXQV_eval(Q,V,T,p):op===3n?Flx_FlxqV_eval(Q,V,T,p)
  :op===4n?pari_FpX_eval(Q,x,T,p):Flx_Flxq_eval(Q,x,T,p));
}

function ff_pari_aut_cache(op:bigint,p:bigint,n:bigint,x:bigint[],a:bigint[],T:bigint[]):string{
 return extensionCompact(op===0n?FpXQ_autpow(x,Number(n),T,p):op===1n?FpXQ_autpowers(x,Number(n),T,p).slice(1):FpXQ_auttrace([x,a],n,T,p));
}

function ff_pari_coefficient_composition(mode:bigint,code:bigint,p:bigint,T:bigint[],plen:bigint[],pflat:bigint[],x:bigint[],vlen:bigint[],vflat:bigint[]):string{
 const m=Number(mode),op=Number(code%10n),P=extensionDecode(plen,pflat),V=extensionDecode(vlen,vflat);
 if(m===3)return extensionCompact(op?pariExtensionBinary.F2x_F2xqV_eval(P[0]as bigint,V as bigint[],T[0]!):pariExtensionBinary.F2x_F2xq_eval(P[0]as bigint,x[0]!,T[0]!));
 let t:any=m===2?T[0]!:T.slice();if(Array.isArray(t)){while(t.length&&t.at(-1)===0n)t.pop();}
 const inverse=code>=10n&&m!==2?extensionInnerQuotient(t,p,m===1).inverse:undefined;
 if(inverse!==undefined)return extensionCompact(coefficientSubstitution(m as 0|1|2,op as 0|1,p,t,P,m===2?x[0]!:x,V as any,inverse));
 const api:any=m===0?pariExtensionGeneric:m===1?pariExtensionWord:pariExtensionBinary;
 const name=(m===0?'FpXY_FpXQ':m===1?'FlxY_Flxq':'F2xY_F2xq')+(op?'V_evalx':'_evalx');
 return extensionCompact(m===2?api[name](P,op?V:x[0]!,t):api[name](P,op?V:x,t,p));
}

function ff_pari_extension_aut(mode:bigint,code:bigint,p:bigint,n:bigint,T:bigint[],phi:bigint[],slen:bigint[],sflat:bigint[],blen:bigint[],bflat:bigint[],alen:bigint[],aflat:bigint[]):string{
 const m=Number(mode),op=Number(code%10n),cache=Number(code%100n/10n),raw=extensionDecode(slen,sflat),B=extensionDecode(blen,bflat),a=extensionDecode(alen,aflat);
 let t:any=m===2?T[0]!:T.slice();if(Array.isArray(t)){while(t.length&&t.at(-1)===0n)t.pop();}
 const inner=code>=100n&&m!==2?extensionInnerQuotient(t,p,m===1).inverse:undefined;
 const S:ExtensionModulus=cache===1?{polynomial:raw,inverse:pariExtensionDivision(m as 0|1|2,2,p,t,raw,[],inner)as ExtensionPolynomial}:cache===2?extensionGetRed(m as 0|1|2,raw,t,p,inner):raw;
 const aut:any=op===0?[m===2?phi[0]!:phi,B]:op===1&&m!==2?[B,a]:[m===2?phi[0]!:phi,B,a];
 if(inner!==undefined)return extensionCompact(extensionAutomorphism(m as 0|1|2,op as 0|1|2,aut,n,S,t,p,inner));
 const api:any=m===0?pariExtensionGeneric:m===1?pariExtensionWord:pariExtensionBinary;
 const name=(m===0?'FpXQXQ':m===1?'FlxqXQ':'F2xqXQ')+['_autpow','_auttrace','_autsum'][op];
 return extensionCompact(m===2?api[name](aut,n,S,t):api[name](aut,n,S,t,p));
}

function ff_pari_extension_projection(mode:bigint,code:bigint,p:bigint,seed:bigint,n:bigint,T:bigint[],alen:bigint[],aflat:bigint[],blen:bigint[],bflat:bigint[]):string {
 const m=Number(mode)as 0|1,op=Number(code%10n),a=extensionDecode(alen,aflat),b=extensionDecode(blen,bflat),saved=getrand();
 try {
  pari_init_rand();setrand(seed);
  const t=T.slice();while(t.length&&t.at(-1)===0n)t.pop();
  const inner=code>=10n?extensionInnerQuotient(t,p,m===1).inverse:undefined;
  let result:any;
  if(inner!==undefined)result=op<2?extensionProjection(m,op as 0|1,Number(n),t,p,a,b,inner):extensionTruncatedPolynomial(m,op===3,Number(n),t,p,a,b,inner);
  else if(op===0)result=m===0?pariExtensionGeneric.random_FpXQX(Number(n),t,p):pariExtensionWord.random_FlxqX(Number(n),t,p);
  else if(op===1)result=m===0?pariExtensionGeneric.FpXQX_dotproduct(a,b,t,p):extensionWordDot(a as bigint[][],b as bigint[][],t,p);
  else {const api:any=m===0?pariExtensionGeneric:pariExtensionWord,name=(m===0?'FpXQXn':'FlxqXn')+(op===2?'_mul':'_sqr');result=op===2?api[name](a,b,Number(n),t,p):api[name](a,Number(n),t,p);}
  return extensionCompact(result)+'|'+getrand();
 }finally{setrand(saved);}
}

function ff_pari_extension_minpoly(mode:bigint,code:bigint,p:bigint,seed:bigint,T:bigint[],slen:bigint[],sflat:bigint[],xlen:bigint[],xflat:bigint[]):string {
 const m=Number(mode)as 0|1,raw=extensionDecode(slen,sflat),x=extensionDecode(xlen,xflat),saved=getrand();
 try {
  pari_init_rand();setrand(seed);
  const t=T.slice();while(t.length&&t.at(-1)===0n)t.pop();
  const inner=code>=100n?extensionInnerQuotient(t,p,m===1).inverse:undefined;
  const S:ExtensionModulus=code%100n===10n?{polynomial:raw,inverse:pariExtensionDivision(m,2,p,t,raw,[],inner)as ExtensionPolynomial}:code%100n===20n?extensionGetRed(m,raw,t,p,inner):raw;
  const result=inner!==undefined?extensionMinimalPolynomial(m,x,S,t,p,inner):m===0?pariExtensionGeneric.FpXQXQ_minpoly(x,S,t,p):pariExtensionWord.FlxqXQ_minpoly(x as bigint[][],S as ExtensionModulus<bigint[][]>,t,p);
  return extensionCompact(result)+'|'+getrand();
 }finally{setrand(saved);}
}

function ff_pari_prime_quotient_cache(op:bigint,p:bigint,n:bigint,T:bigint[],a:bigint[]):string {
 const r=op===0n?pari_gcd_power(a,n,T,p):op===1n?pari_gcd_inverse(a,T,p):op===2n?pariPrimeFrobenius(T,p):op===3n?pari_gcd_inverse(a,T,p**n,p):op===4n?pariPrimePowerAlias(a,n,T,p):pari_gcd_inverse(a,T,p,undefined);
 return '['+r.join(',')+']';
}

function ff_pari_extension_frobenius(mode:bigint,code:bigint,p:bigint,T:bigint[],slen:bigint[],sflat:bigint[],alen:bigint[],aflat:bigint[]):string {
 const m=Number(mode)as 0|1|2,op=Number(code%10n),raw=extensionDecode(slen,sflat),a=extensionDecode(alen,aflat),t=m===2?T[0]!:T;
 if(op===2)return extensionCompact(pariWordFrobenius(T,p));
 // Reject invalid converted coefficient contexts before the fixture builds a cache.
 if(m===0&&p>=2n&&p<1n<<64n){let n=T.length;while(n&&T[n-1]!%p===0n)n--;if(n<2)throw new RangeError('extension modulus must have positive degree');}
 const S:ExtensionModulus=code/10n===1n?{polynomial:raw,inverse:pariExtensionDivision(m,2,p,t,raw,[])as ExtensionPolynomial}:code/10n===2n?extensionGetRed(m,raw,t,p):raw;
 const api:any=pariExtensionFactor,name=(m===0?'FpXQX':m===1?'FlxqX':'F2xqX')+(op===0?'_Frobenius':'Q_halfFrobenius');
 return extensionCompact(m===2?api[name](S,t):op===0?api[name](S,t,p):api[name](a,S,t,p));
}

function ff_pari_prime_linear(op:bigint,p:bigint,c:bigint,a:bigint[],b:bigint[]):string {
 const v=op===0n?pariPrimeLinear.FpX_red(a,p):op===1n?pariPrimeLinear.FpX_add(a,b,p):op===2n?pariPrimeLinear.FpX_sub(a,b,p):op===3n?pariPrimeLinear.FpX_neg(a,p):pariPrimeLinear.FpX_Fp_mul(a,c,p);
 return '['+v.join(',')+']';
}

function ff_pari_extension_root_count(mode:bigint,op:bigint,p:bigint,T:bigint[],flen:bigint[],fflat:bigint[]):string {
 const f=extensionDecode(flen,fflat),R=pariExtensionFactor;
 const v=mode===3n?R.FqX_nbroots(f,null,p):mode===4n?R.FqX_nbroots(f,T,p):mode===2n?R.F2xqX_nbroots(f as bigint[],T[0]!):op===0n?R.FpXQX_split_part(f,T,p):op===1n?mode===0n?R.FpXQX_nbroots(f,T,p):R.FlxqX_nbroots(f as bigint[][],T,p):op===2n?mode===0n?R.FpXQX_is_squarefree(f,T,p):R.FlxqX_is_squarefree(f as bigint[][],T,p):mode===0n?pariExtensionGeneric.FpXX_deriv(f,p):pariExtensionWord.FlxX_deriv(f as bigint[][],p);
 return typeof v==='boolean'?v?'1':'0':extensionCompact(v);
}

function ff_pari_scalar_inverse(op:bigint,p:bigint,a:bigint):string {
 return String(op===1n?pariScalarResidue(a,p):pariScalarInverse(a,p,op===2n));
}

function ff_pari_prime_observation(op:bigint,p:bigint,x:bigint,f:bigint[]):string {
 const v=op===0n?pariObservations.FpX_deriv(f,p):op===1n?pariObservations.FpX_eval(f,x,p):op===2n?pariObservations.FpX_center(f,p,x):op===3n?pariObservations.FpX_div_by_X_x(f,x,p):pariObservations.FpX_is_squarefree(f,p);
 return typeof v==='boolean'?v?'1':'0':Array.isArray(v)?'['+v.join(',')+']':String(v);
}

function ff_pari_fp_scalar(op:bigint,p:bigint,a:bigint,b:bigint,c:bigint):string {
 const P=pariFpScalars;
 const v=op===0n?P.Fp_red(a,p):op===1n?P.Fp_add(a,b,p):op===2n?P.Fp_sub(a,b,p):op===3n?P.Fp_neg(a,p):op===4n?P.Fp_mul(a,b,p):op===5n?P.Fp_sqr(a,p):op===6n?P.Fp_inv(a,p):op===7n?P.Fp_div(a,b,p):op===8n?P.Fp_addmul(a,b,c,p):op===9n?P.Fp_double(a,p):op===10n?P.Fp_halve(a,p):op===11n?P.Fp_center(a,p):op===12n?P.Fp_mulu(a,Number(b),p):P.Fp_eq(a,b);
 return typeof v==='boolean'?v?'1':'0':String(v);
}

function ff_pari_fp_power(op:bigint,p:bigint,a:bigint,e:bigint):string {
 return String(pariFpScalars.Fp_pow(a,e,p));
}

function ff_pari_fp_predicate(op:bigint,p:bigint,a:bigint,e:bigint):string {
 const P=pariFpScalars;
 const v=op===0n?P.Fp_issquare(a,p):op===1n?P.kronecker(a,p):op===2n?P.gcd(a,p):op===3n?P.xgcd(a,p):op===4n?P.Fp_order(a,e,p):P.znorder(a,p,op===5n?undefined:e);
 return typeof v==='boolean'?v?'1':'0':Array.isArray(v)?'['+v.join(', ')+']':String(v);
}

function ff_pari_qpoly(op:bigint,p:bigint,half:bigint,den:bigint,perm:bigint[],L:bigint[],flat:bigint[]):string {
 const P=pariQPoly;const q={num:L,den};
 if(op===0n){const z=P.QPoly_normalize(q);return '['+z.num.join(',')+']/'+z.den;}
 if(op===2n)return '['+P.QPoly_to_FpX(q,p).join(',')+']';
 let z=q;if(op===3n){const n=L.length,M=[Array<bigint>(n+1).fill(0n)];for(let i=0;i<n;i++)M.push([0n,...flat.slice(i*n,(i+1)*n)]);z=P.permtopol([0,...perm.map(Number)],[0n,...L],M,den,p,half);}
 return '['+P.QPoly_to_fractions(z).map(([a,b])=>b===1n?String(a):a+'/'+b).join(',')+']';
}

function ff_pari_galois_integer(op:bigint,a:bigint,b:bigint,k:bigint):string {
 const P=pariQPoly,word=(x:number)=>String(BigInt(x));
 if(op===0n)return word(P.ugcd(Number(a),Number(b)));
 if(op===1n)return word(P.ulcm(Number(a),Number(b)));
 if(op===2n)return '['+P.factoru_small(Number(a)).map(p=>'['+p.map(word).join(',')+']').join(',')+']';
 if(op===3n){const p=new P.Forprime(Number(a));return '['+Array.from({length:Number(k)},()=>word(p.next())).join(',')+']';}
 return String(P.logint(a,b));
}

import {
  ZpX_reduced_resultant,
  ZpX_reduced_resultant_fast,
} from '../../../../packages/parigp-ts/src/base2.js';
import { ZpM_echelon, zlm_echelon } from '../../../../packages/parigp-ts/src/hnf_snf.js';
import { ZX_gcd } from '../../../../packages/parigp-ts/src/QX_factor.js';
import { absZ_factor_limit_strict_default } from '../../../../packages/parigp-ts/src/ifactor.js';
function ff_pari_zx_index(
  op: bigint,
  a: bigint[],
  b: bigint[] | bigint,
  c: bigint,
  p: bigint
): string {
  const P = pariQPoly,
    poly = (v: bigint[]) => '[' + v.join(',') + ']',
    B = b as bigint[];
  switch (Number(op)) {
    case 0:
      return poly(P.ZX_neg(a));
    case 1:
      return poly(P.ZX_Z_mul(a, c));
    case 2:
      return poly(P.ZX_deriv(a));
    case 3:
      return String(P.ZX_resultant(a, B));
    case 4:
      return String(P.ZX_disc(a));
    case 5:
      return P.ZX_is_squarefree(a) ? '1' : '0';
    case 6:
      return String(P.indexpartial(a));
    case 7:
      return String(P.indexpartial(a, c));
    case 8:
      return String(ZpX_reduced_resultant(a, B, p, p ** c));
    case 9:
      return String(ZpX_reduced_resultant_fast(a, B, p, Number(c)));
    case 10: {
      const [f, u] = absZ_factor_limit_strict_default(c);
      return '[' + f.map(poly).join(',') + '],' + (u ? '[' + u.join(', ') + ']' : '0');
    }
    case 11:
      return poly(ZX_gcd(a, B));
    default: {
      const nr = Number(c),
        mat = nr
          ? Array.from({ length: a.length / nr }, (_, j) => a.slice(j * nr, (j + 1) * nr))
          : [];
      const z = (op % 2n === 0n ? ZpM_echelon : zlm_echelon)(mat, op >= 14n, p, b as bigint);
      return z === null ? 'null' : '[' + z.map(poly).join(',') + ']';
    }
  }
}

function ff_pari_zp_lift(
  op: bigint,
  f: bigint[],
  a: bigint[],
  t: bigint[],
  flat: bigint[],
  p: bigint,
  e: bigint
): string {
  const P = pariQPoly,
    Q: bigint[][] = [[]],
    fmt = (v: bigint | bigint[] | bigint[][]): string =>
      Array.isArray(v) ? '[' + v.map((x) => fmt(x)).join(',') + ']' : String(v);
  for (let i = 0; i < flat.length; ) {
    const n = Number(flat[i++]);
    Q.push(flat.slice(i, i + n));
    i += n;
  }
  switch (Number(op)) {
    case 0:
      return fmt(P.ZpX_liftroot(f, a[0]!, p, Number(e)));
    case 1:
      return fmt(P.ZpX_liftroots(f, [0n, ...a], p, Number(e)).slice(1));
    case 2:
      return fmt(P.ZpX_roots(f, p, Number(e)).slice(1));
    case 3:
      return fmt(P.ZpX_liftfact(f, Q, p, Number(e)).slice(1));
    case 4:
      return fmt(P.bezout_lift_fact(f, Q, p, Number(e)).slice(1));
    default:
      return fmt(P.ZpX_ZpXQ_liftroot(f, a, t, p, Number(e)));
  }
}

import { FpV_inv as pariBatchInverse } from '../../../../packages/parigp-ts/src/FpX.js';
import { producttree_scheme as pariProductScheme } from '../../../../packages/parigp-ts/src/bb_group.js';
function ff_pari_vandermonde(op: bigint, L: bigint[], d: bigint, p: bigint): string {
  if (op === 1n) return '[' + pariBatchInverse(L, p).join(',') + ']';
  if (op === 2n) return '[' + pariProductScheme(Number(d)).join(',') + ']';
  const M = pariQPoly
    .FpV_invVandermonde([0n, ...L], d, p)
    .slice(1)
    .map((row) => row.slice(1));
  return '[' + M.map((row) => '[' + row.join(',') + ']').join(',') + ']';
}

function ff_pari_perm_vector(op: bigint, A: bigint[], B: bigint[], exponent: bigint): string {
  const P = pariQPoly,
    a = [0, ...A.map(Number)],
    b = [0, ...B.map(Number)],
    e = Number(exponent);
  const fmt = (v: number | number[] | number[][]): string =>
    Array.isArray(v) ? '[' + v.map((x) => fmt(x)).join(',') + ']' : String(BigInt(v));
  switch (Number(op)) {
    case 0:
      return fmt(P.perm_mul(a, b).slice(1));
    case 1:
      return fmt(P.perm_sqr(a).slice(1));
    case 2:
      return fmt(P.perm_inv(a).slice(1));
    case 3:
      return fmt(P.perm_conj(a, b).slice(1));
    case 4:
      return P.perm_commute(a, b) ? '1' : '0';
    case 5:
      return fmt(P.perm_powu(a, e).slice(1));
    case 6:
      return fmt(
        P.perm_cycles(a)
          .slice(1)
          .map((v) => v.slice(1))
      );
    case 7:
      return fmt(P.perm_orderu(a));
    case 8:
      return fmt(
        P.cyc_pow(P.perm_cycles(a), e)
          .slice(1)
          .map((v) => v.slice(1))
      );
    case 9:
      return fmt(P.vecsmall_lexcmp(a, b));
    case 10:
      return P.zv_equal(a, b) ? '1' : '0';
    case 11:
      return fmt(P.vecsmall_uniq(a).slice(1));
    case 12:
      return fmt(P.group_order({ gen: [[]], ord: a }));
    case 13:
      return fmt(P.group_domain({ gen: e ? [[], a] : [[]], ord: [0] }));
    default:
      return fmt(
        P.vecperm_orbits(e === 0 ? [] : e === 1 ? [a] : [a, b], a.length - 1)
          .slice(1)
          .map((v) => v.slice(1))
      );
  }
}
function ff_pari_group_structure(
  opArg: bigint,
  nArg: bigint,
  g: bigint[],
  go: bigint[],
  h: bigint[],
  ho: bigint[],
  perm: bigint[]
): string {
  const P = pariQPoly,
    op = Number(opArg),
    n = Number(nArg);
  const group = (flat: bigint[], ord: bigint[]): pariQPoly.Group => ({
    gen: [[], ...ord.map((_, i) => [0, ...flat.slice(i * n, (i + 1) * n).map(Number)])],
    ord: [0, ...ord.map(Number)],
  });
  const G = group(g, go),
    H = group(h, ho),
    p = [0, ...perm.map(Number)];
  const fmt = (v: number | number[] | number[][] | number[][][]): string =>
    Array.isArray(v) ? '[' + v.map((x) => fmt(x)).join(',') + ']' : String(v);
  const perms = (v: pariQPoly.Perm[]) => v.slice(1).map((p) => p.slice(1));
  const Gfmt = (g: pariQPoly.Group) => [perms(g.gen), g.ord.slice(1)];
  const Cfmt = (c: pariQPoly.Quotient) => [perms(c.gen), c.coset.slice(1)];
  let v: any;
  const C = () => P.groupelts_quotient(P.group_elts(G, n), H);
  switch (op) {
    case 0:
      v = perms(P.group_elts(G, n));
      break;
    case 1:
      v = perms(P.group_leftcoset(G, p));
      break;
    case 2:
      v = perms(P.group_rightcoset(G, p));
      break;
    case 3:
      v = P.group_set(G, n)
        .slice(1)
        .map((v) => (v ? 1 : 0));
      break;
    case 4:
      v = P.groupelts_set(P.group_elts(G, n), n)
        .slice(1)
        .map((v) => (v ? 1 : 0));
      break;
    case 5:
      v = P.group_perm_normalize(G, p) ? 1 : 0;
      break;
    case 6:
      v = P.group_subgroups(G).map(Gfmt);
      break;
    case 7:
      v = P.group_isA4S4(G);
      break;
    case 8:
      v = Cfmt(P.group_quotient(G, H));
      break;
    case 9:
      v = Cfmt(C());
      break;
    case 10:
      v = P.quotient_perm(C(), p).slice(1);
      break;
    case 11:
      v = Gfmt(P.quotient_group(C(), G));
      break;
    case 12:
      v = Gfmt(P.quotient_subgroup_lift(C(), H, P.trivialgroup()));
      break;
    case 13: {
      const c = C();
      v = Gfmt(P.quotient_subgroup_lift(c, H, P.quotient_group(c, G)));
      break;
    }
    case 14:
      v = P.perm_relorder(p, P.group_set(H, n));
      break;
    case 15:
      v = perms(P.perm_generate(p, P.group_elts(H, n), P.perm_relorder(p, P.group_set(H, n))));
      break;
    case 16:
      v = Gfmt(P.cyclicgroup(p, P.perm_orderu(p)));
      break;
    case 17:
      v = Gfmt(P.dicyclicgroup(G.gen[1]!, G.gen[2]!, G.ord[1]!, G.ord[2]!));
      break;
    case 18:
      v = Gfmt(P.trivialgroup());
      break;
  }
  return fmt(v);
}

function ff_pari_sympol(
  op: bigint,
  width: bigint,
  a: bigint[],
  v: bigint[],
  w: bigint[],
  p: bigint
): string {
  const P = pariQPoly,
    n = Number(width);
  const fmt = (x: any): string => (Array.isArray(x) ? '[' + x.map(fmt).join(',') + ']' : String(x));
  if (op === 4n) return String(pariUnsignedPower(a[0]!, w[0]!, p));
  const M: bigint[][] = [
    [],
    ...Array.from({ length: a.length / n }, (_, i) => [0n, ...a.slice(i * n, (i + 1) * n)]),
  ];
  if (op === 1n) {
    const z = P.sympol_eval({ v: [0, ...v.map(Number)], w: [0, ...w.map(Number)] }, M, p);
    return fmt(Array.isArray(z) ? z.slice(1) : z);
  }
  if (op === 2n) {
    const z = P.fixedfieldsympol(M, p);
    return fmt([z.v.slice(1), z.w.slice(1)]);
  }
  return fmt(
    P.fixedfieldorbits(
      M.map((x) => x.map(Number)),
      [0n, ...v]
    )
      .slice(1)
      .map((x) => x.slice(1))
  );
}

function ff_pari_subgroup(op: bigint, width: bigint, a: bigint[], b: bigint[], p: bigint): string {
  const fmt = (x: any): string => (Array.isArray(x) ? '[' + x.map(fmt).join(',') + ']' : String(x));
  if (op === 0n) return fmt(pariZnstar(p));
  if (op === 1n) return fmt(pariSubgroups(a, [p]));
  if (op === 2n) {
    const n = Number(width);
    return fmt(
      pariHnfModId(
        Array.from({ length: a.length / n }, (_, i) => a.slice(i * n, (i + 1) * n)),
        b
      )
    );
  }
  if (op === 3n)
    return fmt(pariQPoly.listznstarelts(Number(width), Number(p)).map((v) => v.slice(1)));
  return String(pariPrimeGenerator(p));
}


function ff_pari_subgroup_bounds(
  op: bigint,
  width: bigint,
  a: bigint[],
  b: bigint[],
  p: bigint
): string {
  const fmt = (x: any): string => (Array.isArray(x) ? '[' + x.map(fmt).join(',') + ']' : String(x));
  if (op === 7n) {
    const [f, u] = pariStrictFactor(p, a[0]!);
    return fmt([[f.map((x) => x[0]), f.map((x) => x[1])], u ?? 0n]);
  }
  return fmt(pariSubgroups(a, op === 5n ? undefined : op === 8n ? b : p));
}

function ff_pari_galois_validation(op: bigint, a: bigint[], hasden: bigint, den: bigint): string {
  const d = hasden ? den : null;
  if (op === 0n) return pariQPoly.galoisinit(a, d) ? '1' : '0';
  return (
    '[' +
    pariQPoly
      .galoisconj4(a, d)
      .map((p) => '[' + p.num.join(',') + ']/' + p.den)
      .join(',') +
    ']'
  );
}

function ff_pari_galois_actions(
  op: bigint,
  pol: bigint[],
  p: bigint,
  e: bigint,
  L: bigint[],
  m: bigint[],
  den: bigint,
  G: bigint[],
  gen: bigint[],
  ord: bigint[],
  kind: bigint,
  a: bigint[],
  b: bigint[],
  flag: bigint
): string {
  const n = pol.length - 1;
  const perms = (v: bigint[]): pariQPoly.Perm[] => [
    [],
    ...Array.from({ length: v.length / n }, (_, i) => [
      0,
      ...v.slice(i * n, (i + 1) * n).map(Number),
    ]),
  ];
  const gal: pariQPoly.GaloisInit = {
    pol,
    p,
    e: Number(e),
    mod: p ** e,
    roots: [0n, ...L],
    invvdm: [
      [],
      ...Array.from({ length: n }, (_, r) => [
        0n,
        ...Array.from({ length: n }, (_, c) => m[c * n + r]!),
      ]),
    ],
    den,
    group: perms(G),
    gen: perms(gen),
    orders: [0, ...ord.map(Number)],
  };
  const H: pariQPoly.Group = { gen: perms(a), ord: [0, ...b.map(Number)] },
    perm = kind === 0n ? [0, ...a.map(Number)] : kind === 1n ? perms(a) : H;
  const q = (x: pariQPoly.QPoly) => '[' + x.num.join(',') + ']/' + x.den;
  const fmt = (x: unknown): string =>
    Array.isArray(x) ? '[' + x.map(fmt).join(',') + ']' : String(x);
  const gf = (g: pariQPoly.Group) => [g.gen.slice(1).map((p) => p.slice(1)), g.ord.slice(1)];
  if (op === 0n) return q(pariQPoly.galoispermtopol(gal, perm as pariQPoly.Perm));
  if (op === 1n)
    return (
      '[' +
      pariQPoly
        .galoisvecpermtopol(gal, (perm as pariQPoly.Perm[]).slice(1))
        .map(q)
        .join(',') +
      ']'
    );
  if (op === 2n) {
    const z = pariQPoly.galoisfixedfield(gal, perm, Number(flag) as 0);
    return flag === 1n
      ? fmt(z.P)
      : '[' +
          fmt(z.P) +
          ',' +
          q(z.S!) +
          (flag === 2n
            ? ',[' + z.factors!.map((f) => '[' + f.map(q).join(',') + ']').join(',') + ']'
            : '') +
          ']';
  }
  if (op === 3n) return fmt(gf(pariQPoly.galois_group(gal)));
  return fmt(pariQPoly.galoissubgroups(op === 4n ? gal : H).map(gf));
}

function ff_pari_symmetric_search(width: bigint, a: bigint[], w: bigint[], p: bigint): string {
  const n = Number(width),
    NS = [
      [],
      ...Array.from({ length: a.length / n }, (_, i) => [0n, ...a.slice(i * n, (i + 1) * n)]),
    ];
  const s = pariQPoly.fixedfieldsurmer(p, NS, [0, ...w.map(Number)]);
  return s ? JSON.stringify([s.v.slice(1), s.w.slice(1)]) : '0';
}

function ff_pari_findpsi(
  T: bigint[],
  num: bigint[],
  den: bigint,
  o: bigint,
  start: bigint,
  bad: bigint,
  seed: bigint
): string {
  const saved = getrand();
  try {
    setrand(seed);
    const out = { Tmod: [] as pariQPoly.ZX[], psi: [] as number[], p: 0 };
    const p = pariQPoly.findpsi(
      pariQPoly.ZX_disc(T) * den * bad,
      Number(start),
      T,
      { num, den },
      Number(o),
      out
    );
    return (
      '[' +
      p +
      ',[' +
      out.Tmod.slice(1)
        .map((t) => '[' + t.join(',') + ']')
        .join(',') +
      '],[' +
      out.psi.slice(1).join(',') +
      ']]|' +
      getrand()
    );
  } finally {
    setrand(saved);
  }
}
function ff_pari_galois_filter(
  sz: bigint,
  l: bigint[],
  m: bigint[],
  b: bigint,
  q: bigint,
  a: bigint[]
): string {
  const n = Number(sz),
    td = pariQPoly.inittest(
      [0n, ...l],
      [
        [],
        ...Array.from({ length: n }, (_, r) => [
          0n,
          ...Array.from({ length: n }, (_, c) => m[c * n + r]!),
        ]),
      ],
      b,
      q
    ),
    seq = [];
  for (let i = 0; i < a.length; i += n) {
    const ok = pariQPoly.galois_test_perm(td, [0, ...a.slice(i, i + n).map(Number)]);
    seq.push([
      Number(ok),
      td.order.slice(1),
      td.PV.slice(1).map((M) =>
        M ? M.slice(1).map((c) => c.slice(1).map((x) => BigInt.asIntN(64, x))) : 0
      ),
    ]);
  }
  const fmt = (x: unknown): string =>
    Array.isArray(x) ? '[' + x.map(fmt).join(',') + ']' : String(x);
  return fmt(seq);
}
function ff_pari_galois_helpers(op: bigint, sz: bigint, a: bigint[], b: bigint, f: bigint): string {
  const n = Number(sz);
  let z: unknown;
  if (op === 0n) z = pariQPoly.ZX_sub(a, b);
  else if (op === 1n) z = Number(pariQPoly.ZX_equal(a, b));
  else if (op === 2n) z = pariQPoly.ZX_xn(n);
  else if (op === 3n) {
    const lo = Array.from({ length: a.length / n }, (_, i) => [
      0,
      ...a.slice(i * n, (i + 1) * n).map(Number),
    ]);
    z = pariQPoly.galoisfindgroups(lo, [0, ...b.map(Number)], Number(f)).map((x) => x.slice(1));
  } else if (op === 6n) z = nativeWordTotient(sz);
  else z = BigInt(op === 4n ? pariQPoly.radicalu(n) : pariQPoly.eulerphiu(n));
  const fmt = (x: unknown): string =>
    Array.isArray(x) ? '[' + x.map(fmt).join(',') + ']' : String(x);
  return fmt(z);
}

function ff_pari_galois_kernel(
  op: bigint,
  rows: bigint,
  cols: bigint,
  a: bigint[],
  p: bigint,
  b: bigint
): string {
  if (op !== 0n) return String(pariQPoly.lcmBig(p, b));
  const nr = Number(rows),
    nc = Number(cols),
    M = Array.from({ length: nr }, (_, i) => a.slice(i * nc, (i + 1) * nc));
  return (
    '[' +
    pariQPoly
      .FpM_ker(M, nr, nc, p)
      .map((v) => '[' + v.join(',') + ']')
      .join(',') +
    ']'
  );
}

function ff_pari_fixedfield_reprime(
  T: bigint[],
  prime: bigint,
  e: bigint,
  roots: bigint[],
  num: bigint[],
  den: bigint,
  perm: bigint[],
  o: bigint,
  start: bigint,
  seed: bigint
): string {
  const saved = getrand();
  try {
    const val = Number(e),
      mod = prime ** e,
      gb = {
        l: prime,
        valabs: val,
        valsol: val,
        ladicabs: mod,
        ladicsol: mod,
        dis: pariQPoly.ZX_disc(T),
        bornesol: 0n,
      },
      gf = {
        p: Number(start),
        fp: 0,
        deg: Number(o),
        Tmod: [] as pariQPoly.ZX[],
        psi: [] as number[],
      };
    setrand(seed);
    const z = pariQPoly.galoisgenfixedfield0(
      pariQPoly.perm_cycles([0, ...perm.map(Number)]),
      [0n, ...roots],
      { num, den },
      T,
      null,
      gf,
      gb
    );
    const fmt = (x: unknown): string =>
      Array.isArray(x) ? '[' + x.map(fmt).join(',') + ']' : String(x);
    return (
      fmt([
        z ? [[z.PG.gen.slice(1).map((p) => p.slice(1)), z.PG.orders.slice(1)], z.Pg.slice(1)] : 0,
        z ? [[z.V.sym.v.slice(1), z.V.sym.w.slice(1)], z.V.PL.slice(1), z.V.P] : 0,
        gf.p,
        gf.Tmod.slice(1),
        gf.psi.slice(1),
      ]) +
      '|' +
      getrand()
    );
  } finally {
    setrand(saved);
  }
}

function ff_pari_factor_stages(op: bigint, n: bigint, rounds: bigint, seed: bigint): string {
  const value =
    op === 0n
      ? squfof(n)
      : op === 1n
        ? pollardbrent(n)
        : Z_pollardbrent(n, Number(rounds), Number(seed));
  return value ? '[' + value.join(',') + ']' : '0';
}

function ff_pari_power_helpers(op: bigint, n: bigint, a: bigint, b: bigint, c: bigint): string {
  let z: unknown;
  switch (Number(op)) {
    case 0:
      z = pariPower.is_357_power(n, Number(a));
      break;
    case 1:
      z = pariPower.is_kth_power(n, Number(a)) ?? 0n;
      break;
    case 2:
      z = pariPower.Z_isanypower(n);
      break;
    case 3: {
      let last = 0;
      const it = pariPower.forprime(Number(a), Number(b));
      const t = {
        next: () => {
          const v = it.next();
          last = v.done ? 0 : v.value;
          return v.done ? null : v.value;
        },
      };
      const [e, root] = pariPower.is_pth_power(n, t, Number(c));
      z = [e, root, last];
      break;
    }
    case 4:
      z = pariPower.tridiv_bound(n);
      break;
    case 5: {
      const root = pariPower.Z_issquareall(n);
      z = root === null ? [] : [root];
      break;
    }
    case 6:
      z = pariPower.isprimepower(n) ?? [];
      break;
    default:
      throw new Error('unknown power helper');
  }
  const fmt = (x: unknown): string =>
    Array.isArray(x) ? '[' + x.map(fmt).join(',') + ']' : String(x);
  return fmt(z);
}

function ff_pari_ecm(
  op: bigint,
  n: bigint,
  nbc: bigint,
  seed: bigint,
  b: bigint,
  rounds: bigint
): string {
  let factor: bigint | null = null;
  if (op === 2n) factor = pariPower.ellfacteur(n, false);
  else if (op === 3n) factor = pariPower.ellfacteur(n, true, Number(rounds));
  else {
    const state = new ECM(n, Number(nbc), Number(seed));
    for (let i = 0n; i < rounds; i++) {
      factor = state.round(Number(b));
      if (factor) break;
    }
  }
  return String(factor ?? 0n);
}

function ff_pari_mpqs_relations(
  n: bigint,
  q: bigint,
  y1: bigint,
  r1: bigint[],
  y2: bigint,
  r2: bigint[],
  mode: bigint
): string {
  const h = mpqsInternals.newHandle();
  h.N = n;
  h.size_of_FB = 6;
  mpqsInternals.mpqs_FB_ctor(h).p.set([0, 0, 2, 3, 5, 7, 11, 13]);
  const z = mpqsInternals.combine_large_primes(
    h,
    Number(q),
    { Y: y1, relp: r1.map(Number) },
    { Y: y2, relp: r2.map(Number) },
    Number(mode)
  );
  return z === null
    ? '0'
    : typeof z === 'bigint'
      ? String(z)
      : '[' + z.Y + ',[' + z.relp.join(',') + ']]';
}

function ff_pari_mpqs_sqrt(a: bigint, p: bigint): string {
  return String(mpqsInternals.Fl_sqrt(Number(a), Number(p)));
}

import { F2Ms_ker, F2Ms_colelim } from '../../../../packages/parigp-ts/src/F2v.js';
function ff_pari_sparse_kernel(rows: bigint, seed: bigint, packed: bigint[]): string {
  const M: bigint[][] = [];
  for (let i = 0; i < packed.length;) { const n = Number(packed[i++]!); M.push(packed.slice(i, i+n)); i += n; }
  const saved = getrand();
  try {
    setrand(seed);
    const matrix = M.map((c) => c.map(Number)),
      n = Number(rows);
    const p = F2Ms_colelim(matrix, n),
      K = F2Ms_ker(matrix, n);
    // Exercise the MPQS word adapter with the same initial random state.
    const state = getrand();
    setrand(seed);
    const adapted = mpqsInternals.F2Ms_ker(matrix, n).map((v) => {
      let bits = 0n;
      for (let i = v.length - 1; i >= 0; i--) bits = (bits << 32n) | BigInt(v[i]!);
      return bits >> 1n;
    });
    if (adapted.join(',') !== K.join(',') || getrand() !== state)
      throw new Error('MPQS sparse-kernel adapter mismatch');
    return '[' + '[' + p.join(',') + '],[' + K.join(',') + '],' + state + ']';
  } finally {
    setrand(saved);
  }
}

function ff_pari_mpqs_symbols(op: bigint, x: bigint, y: bigint): string {
  return String(
    op === 1n ? mpqsInternals.kroiu(x, Number(y)) : mpqsInternals.krouu(Number(x), Number(y))
  );
}
function ff_pari_mpqs_debug(
  op: bigint,
  n: bigint,
  y: bigint,
  relp: bigint[],
  q: bigint,
  mode: bigint
): string {
  const h = mpqsInternals.newHandle();
  h.N = n;
  h.size_of_FB = 6;
  mpqsInternals.mpqs_FB_ctor(h).p.set([0, 0, 2, 3, 5, 7, 11, 13]);
  const r = relp.map(Number);
  if (op === 0n) return String(mpqsInternals.mpqs_factorback(h, r));
  mpqsInternals.mpqs_check_rel(h, { Y: y, relp: r }, Number(q), Number(mode));
  return '1';
}

function ff_pari_mpqs_inverse(a: bigint, p: bigint): string {
  return String(mpqsInternals.Fl_inv(Number(a), Number(p)));
}
function ff_pari_mpqs_fb(n: bigint, size: bigint, want: bigint): string {
  const S = mpqsInternals,
    h = S.newHandle();
  h.N = h.kN = n;
  h.size_of_FB = Number(size);
  h._k = S.cand_multipliers[0]!;
  h.index0_FB = 3;
  h.pmin_index1 = 3;
  const f = S.mpqs_create_FB(h, !!want);
  return JSON.stringify(
    f
      ? [f]
      : [
          0,
          h.index1_FB,
          h.largest_FB_p,
          Array.from({ length: h.size_of_FB - 1 }, (_, i) => [
            h.FB.p[i + 3],
            h.FB.sqrt_kN[i + 3],
            h.FB.flags[i + 3],
          ]),
        ]
  );
}
function ff_pari_mpqs_candidates(M: bigint, t: bigint, p: bigint): string {
  const S = mpqsInternals,
    h = S.newHandle();
  h.M = Number(M);
  h.sieve_threshold = Number(t);
  S.mpqs_sieve_array_ctor(h);
  for (let i = 0; i < 2 * h.M; i++)
    h.sieve_array[i] =
      p === 0n
        ? 0
        : p === 1n
          ? 255
          : p === 2n
            ? i % 2
              ? 0
              : 255
            : p === 3n
              ? i % 256
              : p === 4n
                ? (i * 73 + 19) % 256
                : i % 8 === 7
                  ? 127
                  : 255;
  const n = S.mpqs_eval_sieve(h);
  return JSON.stringify([n, [...h.candidates.subarray(0, n + 1)]]);
}
function ff_pari_mpqs_init(D: bigint, L: bigint, R: bigint, missing: bigint[]): string {
  const S = mpqsInternals,
    h = S.newHandle();
  h.N = h.kN = D;
  const d = Number(D & 15n);
  h.two_is_norm = Number(![0, 4, 5, 13].includes(d));
  h.two_is_bad = Number([0, 4].includes(d));
  h._k = S.cand_multipliers[0]!;
  S.mpqs_set_parameters(h);
  h.size_of_FB = Number(L) + Number(!h.two_is_norm);
  S.mpqs_create_FB(h, false);
  S.mpqs_sieve_array_ctor(h);
  S.mpqs_poly_ctor(h);
  h.lp_bound = h.largest_FB_p;
  h.dkN = Number(-D);
  S.mpqs_set_sieve_threshold(h);
  const ok = S.mpqs_locate_A_range(h),
    rounds: unknown[] = [];
  const out = [
    Number(ok),
    h.two_is_norm,
    h.two_is_bad,
    h.size_of_FB,
    h.index1_FB,
    h.largest_FB_p,
    h.sieve_threshold,
    Array.from({ length: h.size_of_FB - 1 }, (_, i) => {
      const j = i + 3;
      return [h.FB.p[j], h.FB.sqrt_kN[j], h.FB.flags[j]];
    }),
    rounds,
  ];
  if (ok) {
    h.index_j = 0xffffffff;
    for (let i = 0; i < Number(R); i++) {
      if (!S.mpqs_self_init(h, missing.length ? missing.map(Number) : null)) {
        rounds.push(0);
        break;
      }
      S.mpqs_sieve(h);
      const snapshot = [
        h.A,
        h.B,
        h.bin_index,
        h.index_i,
        h.index_j,
        h.index2_FB,
        Array.from({ length: h.omega_A }, (_, j) => [h.per_A_i[j], h.per_A_H[j]]),
        Array.from({ length: h.size_of_FB - 1 }, (_, i) => {
          const j = i + 3;
          return [h.FB.p[j], h.FB.start1[j], h.FB.start2[j], h.FB.logval[j], h.FB.flags[j]];
        }),
        Buffer.from(h.sieve_array.subarray(0, 2 * h.M)).toString('hex'),
        [...h.candidates.subarray(0, S.mpqs_eval_sieve(h))],
      ];
      rounds.push(snapshot);
    }
  }
  const fmt = (v: unknown): string =>
    Array.isArray(v)
      ? '[' + v.map(fmt).join(',') + ']'
      : typeof v === 'bigint'
        ? String(v)
        : JSON.stringify(v);
  return fmt(out);
}

import { mpqs as pariMpqs } from '../../../../packages/parigp-ts/src/mpqs.js';
import { RelationTable, relationHash } from '../../../../packages/parigp-ts/src/_mpqs_hash.js';
function ff_pari_mpqs_driver(n: bigint, seed: bigint): string {
  const saved = getrand();
  try {
    setrand(seed);
    const r = pariMpqs(n);
    return (
      '[' +
      (r === null ? 'null' : '[' + r.map((v) => '[' + v.join(',') + ']').join(',') + ']') +
      ',' +
      getrand() +
      ']'
    );
  } finally {
    setrand(saved);
  }
}
function ff_pari_mpqs_hash(size: bigint, packed: bigint[]): string {
  const table = new RelationTable(Number(size)),
    hashes: bigint[] = [];
  for (let i = 0; i < packed.length; ) {
    const Y = packed[i++]!,
      n = Number(packed[i++]!),
      relp = packed.slice(i, i + n).map(Number);
    i += n;
    const row = { Y, relp };
    hashes.push(relationHash(row));
    table.add(row);
  }
  return (
    '[[' +
    hashes.join(',') +
    '],[' +
    [...table.values()].map((r) => '[' + r.Y + ',[' + r.relp.join(',') + ']]').join(',') +
    '],' +
    table.size +
    ']'
  );
}

function ff_pari_mpqs_warning(n: bigint, packed: bigint[]): string {
  const P = mpqsInternals,
    h = P.newHandle();
  h.N = n;
  h.debug = true;
  h.size_of_FB = 6;
  P.mpqs_FB_ctor(h).p.set([0, 0, 2, 3, 5, 7, 11, 13]);
  const table = new RelationTable(20);
  for (let i = 0; i < packed.length; ) {
    const Y = packed[i++]!,
      len = Number(packed[i++]!),
      relp = packed.slice(i, i + len).map(Number);
    i += len;
    table.add({ Y, relp });
  }
  let warnings = 0;
  const saved = console.warn;
  console.warn = (message, ...rest) => {
    if (message !== 'MPQS: wrong relation found after Gauss' || rest.length)
      throw new Error('unexpected MPQS warning');
    warnings++;
  };
  try {
    const r = P.mpqs_solve_linear_system(h, table);
    return (
      '[' +
      (r === null ? 'null' : '[' + r.map((v) => '[' + v.join(',') + ']').join(',') + ']') +
      ',' +
      warnings +
      ']'
    );
  } finally {
    console.warn = saved;
  }
}
function ff_pari_mpqs_class_candidates(D: bigint, L: bigint, R: bigint, missing: bigint[]): string {
  const S = mpqsInternals,
    h = S.newHandle();
  h.N = h.kN = D;
  h.debug = true;
  const d = Number(D & 15n);
  h.two_is_norm = Number(![0, 4, 5, 13].includes(d));
  h.two_is_bad = Number([0, 4].includes(d));
  h._k = S.cand_multipliers[0]!;
  S.mpqs_set_parameters(h);
  h.size_of_FB = Number(L) + Number(!h.two_is_norm);
  S.mpqs_create_FB(h, false);
  S.mpqs_sieve_array_ctor(h);
  S.mpqs_poly_ctor(h);
  h.lp_bound = h.largest_FB_p;
  h.dkN = Number(-D);
  S.mpqs_set_sieve_threshold(h);
  const ok = S.mpqs_locate_A_range(h),
    table = new RelationTable(h.target_rels),
    lprel = new Map<number, { Y: bigint; relp: number[] }>(),
    rounds: unknown[] = [];
  if (ok) {
    h.index_j = 0xffffffff;
    for (let i = 0; i < Number(R); i++) {
      if (!S.mpqs_self_init(h, missing.length ? missing.map(Number) : null)) {
        rounds.push([0]);
        break;
      }
      S.mpqs_sieve(h);
      const tc = S.mpqs_eval_sieve(h),
        factor = S.mpqs_eval_cand(h, tc, table, lprel, 1);
      rounds.push([1, tc, table.size, factor ?? 0]);
    }
  }
  const fmt = (v: unknown): string =>
    Array.isArray(v) ? '[' + v.map(fmt).join(',') + ']' : String(v);
  return fmt([Number(ok), rounds, [...table.values()].map((r) => [r.Y, r.relp])]);
}


import * as precisionLift from '../../../../packages/parigp-ts/src/Zp.js';
import {
  FpX_red as precisionRed,
  FpX_sub as precisionSub,
  FpXQ_mul as precisionMul,
} from '../../../../packages/parigp-ts/src/ffinit.js';
functions.ff_pari_zp_precision = (
  op: bigint,
  f: bigint[],
  a: bigint[],
  T: bigint[],
  p: bigint,
  e: bigint
): string => {
  const trace: unknown[] = [],
    n = Number(e),
    q = op === 27n || op === 28n ? 1n : p ** e;
  let r: bigint[] | bigint | null;
  while (a.length && a.at(-1) === 0n) a = a.slice(0, -1);
  if ((op >= 13n && op <= 15n) || op === 21n) {
    const aa = a[0] ?? 0n,
      exponent = op === 13n || op === 21n ? T[0]! : 2n;
    const b = op === 15n ? (f[0] ?? 0n) : pariFpScalars.Fp_pow(aa, exponent, q) + p * (f[0] ?? 0n);
    r =
      (op === 13n || op === 21n)
        ? precisionLift.Zp_sqrtnlift(b, exponent, aa, p, n)
        : op === 14n
          ? precisionLift.Zp_sqrtlift(b, aa, p, n)
          : precisionRoots.Zp_sqrt(b, p, n);
  } else if (op === 22n) r = precisionLift.Zp_exp((f[0] ?? 0n) * (p === 2n ? 4n : p), p, n);
  else if (op === 23n) {
    const x = f.map(c => c * p);
    x[0] = (x[0] ?? 0n) + 1n;
    r = precisionLift.ZpXQ_log(x, T, p, n);
  } else if (op === 24n) r = precisionLift.Zp_inv(f[0] ?? 0n, p, n);
  else if (op === 25n) r = precisionLift.Zp_div(a[0] ?? 0n, f[0] ?? 0n, p, n);
  else if (op === 26n) r = precisionLift.Zp_invlift(f[0] ?? 0n,
    n === 1 ? (a[0] ?? 0n) : precisionLift.Zp_inv(f[0] ?? 0n, p, 1), p, n);
  else if (op === 27n || op === 28n) {
    let x = precisionCvtop(f[0] ?? 0n, p, n);
    if (op === 28n) x = precisionExp({ ...x, valuation: x.valuation + Number((a[0] ?? 0n) % 5n) });
    r = [x.unit, BigInt(x.valuation), BigInt(x.precision), x.modulus];
  } else if (op === 29n) {
    const x = f.map(c => c * p);
    x[0] = (x[0] ?? 0n) + 1n;
    r = precisionFrob.ZpXQ_sqrtnorm(x, T, q, p, n);
  } else if (op === 30n || op === 31n) {
    const Tp = precisionRed(T, p);
    const root = precisionPow([0n,1n], p ** BigInt(Tp.length - 2), Tp, p);
    const sqx = precisionWordPowers(root, Number(p) - 1, Tp, p);
    const x = precisionRem(precisionRed(f, p), Tp, p);
    if (op === 31n) r = precisionLroot(x, sqx, Tp, p, 0n);
    else {
      T = precisionLift.Flx_Teichmuller(Tp, p, n);
      const Xm = precisionPowers([...Array<bigint>(T.length - 1).fill(0n), 1n], Number(p) - 1, T, q);
      r = precisionFrob.Teichmuller_lift(x, Xm, T, sqx, Tp, p, 0n, n);
    }
  }
  else if (op >= 16n && op <= 20n) {
    if (op === 17n) T = precisionLift.Flx_Teichmuller(precisionRed(T, p), p, n);
    let x = op === 19n && T.length === 2 ? f : precisionRem(precisionRed(f, q), T, q);
    if (op === 16n) r = precisionFrob.ZpXQ_frob_cyc(x, T, q, p);
    else if (op === 17n)
      r = precisionFrob.ZpXQ_frob(
        x,
        precisionPowers([...Array<bigint>(T.length - 1).fill(0n), 1n], Number(p) - 1, T, q),
        T,
        q,
        p
      );
    else if (op === 18n) r = precisionFrob.ZpXQ_frob(x, [], T, q, p);
    else {
      if (!precisionRed(x, p).length) x = [1n];
      r =
        op === 19n
          ? precisionFrob.ZpXQ_norm_pcyc(x, T, q, p)
          : precisionFrob.ZpXQ_sqrtnorm_pcyc(precisionMul(x, x, T, q), T, q, p, n);
    }
  } else if (op === 12n) {
    const driver = precisionLift.gen_ZpX_Newton;
    const spy = precisionSpyOn(precisionLift, 'gen_ZpX_Newton').mockImplementation(
      (x, p, n, evaluate, invd) =>
        driver(
          x,
          p,
          n,
          (x, q) => {
            const v = evaluate(x, q);
            trace.push([0n, q, x, v[0]]);
            return v;
          },
          (V, v, q, M) => {
            const r = invd(V, v, q, M);
            trace.push([1n, q, BigInt(M), V, v[1], r]);
            return r;
          }
        )
    );
    try {
      r = pariQPoly.FpXQ_inv(f, T, q, p);
    } finally {
      spy.mockRestore();
    }
  } else if (op === 11n) r = precisionLift.Flx_Teichmuller(precisionRed(f, p), p, n);
  else if (op === 6n) r = precisionLift.ZpXQ_inv(f, T, p, n);
  else if (op === 7n)
    r = precisionLift.ZpXQ_invlift(f, n === 1 ? a : precisionLift.ZpXQ_inv(f, T, p, 1), T, p, n);
  else if (op === 8n) r = precisionLift.ZpXQ_div(a, f, T, q, p, n);
  else {
    const ai = precisionLift.ZpXQ_inv(f, T, p, 1);
    if (op === 9n)
      r = precisionLift.gen_ZpX_Newton<[bigint[], bigint[]]>(
        ai,
        p,
        n,
        (x, q) => {
          const f1 = precisionSub(
            precisionMul(x, precisionRed(f, q), precisionRed(T, q), q),
            [1n],
            q
          );
          trace.push([0n, q, x, f1]);
          return [f1, x];
        },
        (V, v, q, M) => {
          const r = precisionMul(V, v[1], precisionRed(T, q), q);
          trace.push([1n, q, BigInt(M), V, v[1], r]);
          return r;
        }
      );
    else
      r = precisionLift.gen_ZpX_Dixon<[bigint[][], bigint[]]>(
        [[f], T],
        a,
        q,
        p,
        n,
        (F, d, q) => {
          const r = precisionMul(F[0][0]!, d, F[1], q);
          trace.push([0n, q, F, d, r]);
          return r;
        },
        (d) => {
          const r = precisionMul(ai, d, precisionRed(T, p), p);
          trace.push([1n, d, r]);
          return r;
        }
      );
  }
  return JSON.stringify([r, trace], (_, v) => (typeof v === 'bigint' ? String(v) : v)).replace(
    /"(-?\d+)"/g,
    '$1'
  );
};

import { spyOn as precisionSpyOn } from 'bun:test';

import * as precisionRoots from '../../../../packages/parigp-ts/src/qfb.js';
import { cvtop as precisionCvtop } from '../../../../packages/parigp-ts/src/gen2.js';
import { Qp_exp as precisionExp } from '../../../../packages/parigp-ts/src/trans1.js';
import { FpXQ_pow as precisionPow } from '../../../../packages/parigp-ts/src/ffinit.js';
import { Flxq_lroot_fast_pre as precisionLroot, Flxq_powers as precisionWordPowers } from '../../../../packages/parigp-ts/src/Flx.js';
import * as precisionFrob from '../../../../packages/parigp-ts/src/FlxqE.js';
import {
  FpX_rem as precisionRem,
  FpXQ_powers as precisionPowers,
} from '../../../../packages/parigp-ts/src/FpX.js';

functions.ff_pari_zp_binary_linear = (b:bigint,e:bigint):string =>
  String(precisionLift.Zp_sqrtnlift(2n*b+1n,1n,1n,2n,Number(e)));
