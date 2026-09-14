/** PARI basemath/FpXX.c: polynomials over a prime-field extension.
 * @see Deviation: PARI extension-polynomial coefficient adapters
 */
import { extensionPolynomial, type ExtensionPolynomial } from './_extension_polynomial.js';
/** Native FpXQX_mul; ascending outer and inner coefficients.
 * @see Deviation: PARI extension-polynomial coefficient adapters
 */
export function FpXQX_mul(
  x: ExtensionPolynomial,
  y: ExtensionPolynomial,
  T: bigint[],
  p: bigint
): ExtensionPolynomial {
  return extensionPolynomial(0, 0, p, T, x, y) as ExtensionPolynomial;
}
/** Native FpXQX_sqr; ascending outer and inner coefficients.
 * @see Deviation: PARI extension-polynomial coefficient adapters
 */
export function FpXQX_sqr(x: ExtensionPolynomial, T: bigint[], p: bigint): ExtensionPolynomial {
  return extensionPolynomial(0, 1, p, T, x) as ExtensionPolynomial;
}
/** Native FpXQX_red; ascending outer and inner coefficients.
 * @see Deviation: PARI extension-polynomial coefficient adapters
 */
export function FpXQX_red(x: ExtensionPolynomial, T: bigint[], p: bigint): ExtensionPolynomial {
  return extensionPolynomial(0, 2, p, T, x) as ExtensionPolynomial;
}

import {extensionDivision,extensionGetRed,type ExtensionModulus} from './_extension_division.js';
/** Native FpXQX_divrem, with PARI's coefficient and dispatch rules.
 * @see Deviation: PARI extension-polynomial division adapters
 */
export function FpXQX_divrem(x:ExtensionPolynomial, S:ExtensionModulus<ExtensionPolynomial>, T:bigint[], p:bigint):[ExtensionPolynomial,ExtensionPolynomial]{
 return extensionDivision(0,0,p,T,x,S) as [ExtensionPolynomial,ExtensionPolynomial];
}
/** Native FpXQX_rem, with PARI's coefficient and dispatch rules.
 * @see Deviation: PARI extension-polynomial division adapters
 */
export function FpXQX_rem(x:ExtensionPolynomial, S:ExtensionModulus<ExtensionPolynomial>, T:bigint[], p:bigint):ExtensionPolynomial{
 return extensionDivision(0,1,p,T,x,S) as ExtensionPolynomial;
}
/** Native FpXQX_div, with PARI's coefficient and dispatch rules.
 * @see Deviation: PARI extension-polynomial division adapters
 */
export function FpXQX_div(x:ExtensionPolynomial, S:ExtensionModulus<ExtensionPolynomial>, T:bigint[], p:bigint):ExtensionPolynomial{
 return extensionDivision(0,4,p,T,x,S) as ExtensionPolynomial;
}
/** Native FpXQX_invBarrett, with PARI's coefficient and dispatch rules.
 * @see Deviation: PARI extension-polynomial division adapters
 */
export function FpXQX_invBarrett(x:ExtensionPolynomial, T:bigint[], p:bigint):ExtensionPolynomial{
 return extensionDivision(0,2,p,T,x) as ExtensionPolynomial;
}
/** Native cached reciprocal selection; existing reduction objects are retained.
 * @see Deviation: PARI extension-polynomial division adapters
 */
export function FpXQX_get_red(S:ExtensionModulus<ExtensionPolynomial>,T:bigint[],p:bigint):ExtensionModulus<ExtensionPolynomial>{
 return extensionGetRed(0,S,T,p);
}

import { extensionGcd, type ExtensionMatrix } from './_extension_gcd.js';
/** Native unscaled greatest common divisor.
 * @see Deviation: PARI extension-polynomial GCD adapters
 */
export function FpXQX_gcd(x:ExtensionPolynomial, y:ExtensionPolynomial, T:bigint[], p:bigint):ExtensionPolynomial{
  return extensionGcd(0,0,p,T,x,y) as ExtensionPolynomial;
}
/** Native [gcd, U, V] with U*x + V*y = gcd.
 * @see Deviation: PARI extension-polynomial GCD adapters
 */
export function FpXQX_extgcd(x:ExtensionPolynomial, y:ExtensionPolynomial, T:bigint[], p:bigint):[ExtensionPolynomial,ExtensionPolynomial,ExtensionPolynomial]{
  return extensionGcd(0,1,p,T,x,y) as [ExtensionPolynomial,ExtensionPolynomial,ExtensionPolynomial];
}
/** Native half-GCD transformation matrix in row-major order.
 * @see Deviation: PARI extension-polynomial GCD adapters
 */
export function FpXQX_halfgcd(x:ExtensionPolynomial, y:ExtensionPolynomial, T:bigint[], p:bigint):ExtensionMatrix<ExtensionPolynomial>{
  return extensionGcd(0,2,p,T,x,y) as ExtensionMatrix<ExtensionPolynomial>;
}

import {extensionQuotient} from './_extension_quotient.js';
/** Product in the extension quotient.
 * @see Deviation: PARI extension quotient adapters
 */
export function FpXQXQ_mul(x:ExtensionPolynomial, y:ExtensionPolynomial, S:ExtensionModulus<ExtensionPolynomial>, T:bigint[], p:bigint):ExtensionPolynomial{
  return extensionQuotient(0,0,p,0n,T,x,y,S)as ExtensionPolynomial;
}
/** Square in the extension quotient.
 * @see Deviation: PARI extension quotient adapters
 */
export function FpXQXQ_sqr(x:ExtensionPolynomial, S:ExtensionModulus<ExtensionPolynomial>, T:bigint[], p:bigint):ExtensionPolynomial{
  return extensionQuotient(0,1,p,0n,T,x,[],S)as ExtensionPolynomial;
}
/** Inverse or null when no inverse exists.
 * @see Deviation: PARI extension quotient adapters
 */
export function FpXQXQ_invsafe(x:ExtensionPolynomial, S:ExtensionModulus<ExtensionPolynomial>, T:bigint[], p:bigint):ExtensionPolynomial|null{
  return extensionQuotient(0,2,p,0n,T,x,[],S)as ExtensionPolynomial|null;
}
/** Native quotient inverse, with exact GEN error display.
 * @see Deviation: PARI extension quotient adapters
 */
export function FpXQXQ_inv(x:ExtensionPolynomial, S:ExtensionModulus<ExtensionPolynomial>, T:bigint[], p:bigint):ExtensionPolynomial{
  return extensionQuotient(0,3,p,0n,T,x,[],S)as ExtensionPolynomial;
}
/** Divide in the extension quotient.
 * @see Deviation: PARI extension quotient adapters
 */
export function FpXQXQ_div(x:ExtensionPolynomial, y:ExtensionPolynomial, S:ExtensionModulus<ExtensionPolynomial>, T:bigint[], p:bigint):ExtensionPolynomial{
  return extensionQuotient(0,4,p,0n,T,x,y,S)as ExtensionPolynomial;
}
/** Native signed-integer powering schedule.
 * @see Deviation: PARI extension quotient adapters
 */
export function FpXQXQ_pow(x:ExtensionPolynomial, n:bigint, S:ExtensionModulus<ExtensionPolynomial>, T:bigint[], p:bigint):ExtensionPolynomial{
  return extensionQuotient(0,5,p,n,T,x,[],S)as ExtensionPolynomial;
}
/** Powers from zero through l, retaining the raw first power.
 * @see Deviation: PARI extension quotient adapters
 */
export function FpXQXQ_powers(x:ExtensionPolynomial, l:number, S:ExtensionModulus<ExtensionPolynomial>, T:bigint[], p:bigint):ExtensionPolynomial[]{
  if(!Number.isSafeInteger(l)||l<0||l>=0xffffffff)throw new RangeError('power count must be a nonnegative array length');
  return extensionQuotient(0,6,p,BigInt(l),T,x,[],S)as ExtensionPolynomial[];
}

import {extensionComposition} from './_extension_composition.js';
/** Native extension-polynomial composition.
 * @see Deviation: PARI extension composition adapters
 */
export function FpXQX_FpXQXQ_eval(Q:ExtensionPolynomial, x:ExtensionPolynomial, S:ExtensionModulus<ExtensionPolynomial>, T:bigint[], p:bigint):ExtensionPolynomial{
 return extensionComposition(0,0,p,T,Q,x,[],S)as ExtensionPolynomial;
}
/** Native extension-polynomial composition from a supplied power table.
 * @see Deviation: PARI extension composition adapters
 */
export function FpXQX_FpXQXQV_eval(Q:ExtensionPolynomial, V:ExtensionPolynomial[], S:ExtensionModulus<ExtensionPolynomial>, T:bigint[], p:bigint):ExtensionPolynomial{
 return extensionComposition(0,1,p,T,Q,[],V,S)as ExtensionPolynomial;
}

import { coefficientSubstitution } from './_coefficient_composition.js';

/** Native substitution in polynomial coefficients.
 * @see Deviation: PARI coefficient-substitution adapters
 */
export function FpXY_FpXQ_evalx(
  P: ExtensionPolynomial,
  x: bigint[],
  T: bigint[],
  p: bigint
): ExtensionPolynomial {
  return coefficientSubstitution(0, 0, p, T, P, x, []);
}

/** Native substitution in polynomial coefficients from a supplied power table.
 * @see Deviation: PARI coefficient-substitution adapters
 */
export function FpXY_FpXQV_evalx(
  P: ExtensionPolynomial,
  V: bigint[][],
  T: bigint[],
  p: bigint
): ExtensionPolynomial {
  return coefficientSubstitution(0, 1, p, T, P, [], V);
}

import { extensionAutomorphism } from './_extension_automorphism.js';

/** Native extension automorphism power.
 * @see Deviation: PARI extension automorphism adapters
 */
export function FpXQXQ_autpow(
  aut: [bigint[], ExtensionPolynomial],
  n: bigint,
  S: ExtensionModulus<ExtensionPolynomial>,
  T: bigint[],
  p: bigint
): [bigint[], ExtensionPolynomial] {
  return extensionAutomorphism(0, 0, aut, n, S, T, p) as [bigint[], ExtensionPolynomial];
}

/** Native additive trace.
 * @see Deviation: PARI extension automorphism adapters
 */
export function FpXQXQ_auttrace(
  aut: [ExtensionPolynomial, ExtensionPolynomial],
  n: bigint,
  S: ExtensionModulus<ExtensionPolynomial>,
  T: bigint[],
  p: bigint
): [ExtensionPolynomial, ExtensionPolynomial] {
  return extensionAutomorphism(0, 1, aut, n, S, T, p) as [ExtensionPolynomial, ExtensionPolynomial];
}

/** Native multiplicative aggregate.
 * @see Deviation: PARI extension automorphism adapters
 */
export function FpXQXQ_autsum(
  aut: [bigint[], ExtensionPolynomial, ExtensionPolynomial],
  n: bigint,
  S: ExtensionModulus<ExtensionPolynomial>,
  T: bigint[],
  p: bigint
): [bigint[], ExtensionPolynomial, ExtensionPolynomial] {
  return extensionAutomorphism(0, 2, aut, n, S, T, p) as [
    bigint[],
    ExtensionPolynomial,
    ExtensionPolynomial,
  ];
}


import { extensionProjection } from './_extension_projection.js';
import { extensionTruncatedPolynomial } from './_extension_polynomial.js';
/** Native random_FpXQX; ascending coefficient arrays, without variable metadata.
 * @see Deviation: PARI extension projection adapters
 */
export function random_FpXQX(length: number, T: bigint[], p: bigint): ExtensionPolynomial {
  return extensionProjection(0, 0, length, T, p) as ExtensionPolynomial;
}
/** Native FpXQX_dotproduct; ascending coefficient arrays, without variable metadata.
 * @see Deviation: PARI extension projection adapters
 */
export function FpXQX_dotproduct(
  x: ExtensionPolynomial,
  y: ExtensionPolynomial,
  T: bigint[],
  p: bigint
): bigint | bigint[] {
  return extensionProjection(0, 1, 0, T, p, x, y) as bigint | bigint[];
}
/** Native FpXQXn_mul; ascending coefficient arrays, without variable metadata.
 * @see Deviation: PARI extension projection adapters
 */
export function FpXQXn_mul(
  x: ExtensionPolynomial,
  y: ExtensionPolynomial,
  n: number,
  T: bigint[],
  p: bigint
): ExtensionPolynomial {
  return extensionTruncatedPolynomial(0, false, n, T, p, x, y) as ExtensionPolynomial;
}
/** Native FpXQXn_sqr; ascending coefficient arrays, without variable metadata.
 * @see Deviation: PARI extension projection adapters
 */
export function FpXQXn_sqr(
  x: ExtensionPolynomial,
  n: number,
  T: bigint[],
  p: bigint
): ExtensionPolynomial {
  return extensionTruncatedPolynomial(0, true, n, T, p, x) as ExtensionPolynomial;
}


import { extensionMinimalPolynomial } from './_extension_minpoly.js';
/** Native FpXQXQ_minpoly: Shoup projection over the extension coefficient field.
 * @see Deviation: PARI extension minimal-polynomial adapters
 */
export function FpXQXQ_minpoly(
  x: ExtensionPolynomial,
  S: ExtensionModulus<ExtensionPolynomial>,
  T: bigint[],
  p: bigint
): ExtensionPolynomial {
  return extensionMinimalPolynomial(0, x, S, T, p) as ExtensionPolynomial;
}

import { extensionDerivative } from './_extension_derivative.js';
/** Native outer derivative, retaining coefficient tags.
 * @see Deviation: PARI extension root-count and derivative adapters
 */
export function FpXX_deriv(x: ExtensionPolynomial, p: bigint): ExtensionPolynomial {
  return extensionDerivative(false, x, p);
}
