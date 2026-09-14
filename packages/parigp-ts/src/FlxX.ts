/** PARI basemath/FlxX.c: polynomials over a prime-field extension.
 * @see Deviation: PARI extension-polynomial coefficient adapters
 */
import { extensionPolynomial } from './_extension_polynomial.js';
/** Native FlxqX_mul; ascending outer and inner coefficients.
 * @see Deviation: PARI extension-polynomial coefficient adapters
 */
export function FlxqX_mul(x: bigint[][], y: bigint[][], T: bigint[], p: bigint): bigint[][] {
  return extensionPolynomial(1, 0, p, T, x, y) as bigint[][];
}
/** Native FlxqX_sqr; ascending outer and inner coefficients.
 * @see Deviation: PARI extension-polynomial coefficient adapters
 */
export function FlxqX_sqr(x: bigint[][], T: bigint[], p: bigint): bigint[][] {
  return extensionPolynomial(1, 1, p, T, x) as bigint[][];
}
/** Native FlxqX_red; ascending outer and inner coefficients.
 * @see Deviation: PARI extension-polynomial coefficient adapters
 */
export function FlxqX_red(x: bigint[][], T: bigint[], p: bigint): bigint[][] {
  return extensionPolynomial(1, 2, p, T, x) as bigint[][];
}
/** Native word normalization inverts even a leading constant one.
 * @see Deviation: PARI extension-polynomial coefficient adapters
 */
export function FlxqX_normalize(x: bigint[][], T: bigint[], p: bigint): bigint[][] {
  return extensionPolynomial(1, 3, p, T, x) as bigint[][];
}

import {extensionDivision,extensionGetRed,type ExtensionModulus} from './_extension_division.js';
/** Native FlxqX_divrem, with PARI's coefficient and dispatch rules.
 * @see Deviation: PARI extension-polynomial division adapters
 */
export function FlxqX_divrem(x:bigint[][], S:ExtensionModulus<bigint[][]>, T:bigint[], p:bigint):[bigint[][],bigint[][]]{
 return extensionDivision(1,0,p,T,x,S) as [bigint[][],bigint[][]];
}
/** Native FlxqX_rem, with PARI's coefficient and dispatch rules.
 * @see Deviation: PARI extension-polynomial division adapters
 */
export function FlxqX_rem(x:bigint[][], S:ExtensionModulus<bigint[][]>, T:bigint[], p:bigint):bigint[][]{
 return extensionDivision(1,1,p,T,x,S) as bigint[][];
}
/** Native FlxqX_div, with PARI's coefficient and dispatch rules.
 * @see Deviation: PARI extension-polynomial division adapters
 */
export function FlxqX_div(x:bigint[][], S:ExtensionModulus<bigint[][]>, T:bigint[], p:bigint):bigint[][]{
 return extensionDivision(1,4,p,T,x,S) as bigint[][];
}
/** Native FlxqX_invBarrett, with PARI's coefficient and dispatch rules.
 * @see Deviation: PARI extension-polynomial division adapters
 */
export function FlxqX_invBarrett(x:bigint[][], T:bigint[], p:bigint):bigint[][]{
 return extensionDivision(1,2,p,T,x) as bigint[][];
}
/** Native cached reciprocal selection; existing reduction objects are retained.
 * @see Deviation: PARI extension-polynomial division adapters
 */
export function FlxqX_get_red(S:ExtensionModulus<bigint[][]>,T:bigint[],p:bigint):ExtensionModulus<bigint[][]>{
 return extensionGetRed(1,S,T,p);
}

import { extensionGcd, type ExtensionMatrix } from './_extension_gcd.js';
/** Native unscaled greatest common divisor.
 * @see Deviation: PARI extension-polynomial GCD adapters
 */
export function FlxqX_gcd(x:bigint[][], y:bigint[][], T:bigint[], p:bigint):bigint[][]{
  return extensionGcd(1,0,p,T,x,y) as bigint[][];
}
/** Native [gcd, U, V] with U*x + V*y = gcd.
 * @see Deviation: PARI extension-polynomial GCD adapters
 */
export function FlxqX_extgcd(x:bigint[][], y:bigint[][], T:bigint[], p:bigint):[bigint[][],bigint[][],bigint[][]]{
  return extensionGcd(1,1,p,T,x,y) as [bigint[][],bigint[][],bigint[][]];
}
/** Native half-GCD transformation matrix in row-major order.
 * @see Deviation: PARI extension-polynomial GCD adapters
 */
export function FlxqX_halfgcd(x:bigint[][], y:bigint[][], T:bigint[], p:bigint):ExtensionMatrix<bigint[][]>{
  return extensionGcd(1,2,p,T,x,y) as ExtensionMatrix<bigint[][]>;
}

import {extensionQuotient} from './_extension_quotient.js';
/** Product in the extension quotient.
 * @see Deviation: PARI extension quotient adapters
 */
export function FlxqXQ_mul(x:bigint[][], y:bigint[][], S:ExtensionModulus<bigint[][]>, T:bigint[], p:bigint):bigint[][]{
  return extensionQuotient(1,0,p,0n,T,x,y,S)as bigint[][];
}
/** Square in the extension quotient.
 * @see Deviation: PARI extension quotient adapters
 */
export function FlxqXQ_sqr(x:bigint[][], S:ExtensionModulus<bigint[][]>, T:bigint[], p:bigint):bigint[][]{
  return extensionQuotient(1,1,p,0n,T,x,[],S)as bigint[][];
}
/** Inverse or null when no inverse exists.
 * @see Deviation: PARI extension quotient adapters
 */
export function FlxqXQ_invsafe(x:bigint[][], S:ExtensionModulus<bigint[][]>, T:bigint[], p:bigint):bigint[][]|null{
  return extensionQuotient(1,2,p,0n,T,x,[],S)as bigint[][]|null;
}
/** Native quotient inverse, with exact GEN error display.
 * @see Deviation: PARI extension quotient adapters
 */
export function FlxqXQ_inv(x:bigint[][], S:ExtensionModulus<bigint[][]>, T:bigint[], p:bigint):bigint[][]{
  return extensionQuotient(1,3,p,0n,T,x,[],S)as bigint[][];
}
/** Divide in the extension quotient.
 * @see Deviation: PARI extension quotient adapters
 */
export function FlxqXQ_div(x:bigint[][], y:bigint[][], S:ExtensionModulus<bigint[][]>, T:bigint[], p:bigint):bigint[][]{
  return extensionQuotient(1,4,p,0n,T,x,y,S)as bigint[][];
}
/** Native signed-integer powering schedule.
 * @see Deviation: PARI extension quotient adapters
 */
export function FlxqXQ_pow(x:bigint[][], n:bigint, S:ExtensionModulus<bigint[][]>, T:bigint[], p:bigint):bigint[][]{
  return extensionQuotient(1,5,p,n,T,x,[],S)as bigint[][];
}
/** Powers from zero through l, retaining the raw first power.
 * @see Deviation: PARI extension quotient adapters
 */
export function FlxqXQ_powers(x:bigint[][], l:number, S:ExtensionModulus<bigint[][]>, T:bigint[], p:bigint):bigint[][][]{
  if(!Number.isSafeInteger(l)||l<0||l>=0xffffffff)throw new RangeError('power count must be a nonnegative array length');
  return extensionQuotient(1,6,p,BigInt(l),T,x,[],S)as bigint[][][];
}
/** Native unsigned-word powering shortcuts.
 * @see Deviation: PARI extension quotient adapters
 */
export function FlxqXQ_powu(x:bigint[][], n:bigint, S:ExtensionModulus<bigint[][]>, T:bigint[], p:bigint):bigint[][]{
  return extensionQuotient(1,7,p,n,T,x,[],S)as bigint[][];
}

import {extensionComposition} from './_extension_composition.js';
/** Native extension-polynomial composition.
 * @see Deviation: PARI extension composition adapters
 */
export function FlxqX_FlxqXQ_eval(Q:bigint[][], x:bigint[][], S:ExtensionModulus<bigint[][]>, T:bigint[], p:bigint):bigint[][]{
 return extensionComposition(1,0,p,T,Q,x,[],S)as bigint[][];
}
/** Native extension-polynomial composition from a supplied power table.
 * @see Deviation: PARI extension composition adapters
 */
export function FlxqX_FlxqXQV_eval(Q:bigint[][], V:bigint[][][], S:ExtensionModulus<bigint[][]>, T:bigint[], p:bigint):bigint[][]{
 return extensionComposition(1,1,p,T,Q,[],V,S)as bigint[][];
}

import { coefficientSubstitution } from './_coefficient_composition.js';

/** Native substitution in polynomial coefficients.
 * @see Deviation: PARI coefficient-substitution adapters
 */
export function FlxY_Flxq_evalx(P: bigint[][], x: bigint[], T: bigint[], p: bigint): bigint[][] {
  return coefficientSubstitution(1, 0, p, T, P, x, []) as bigint[][];
}

/** Native substitution in polynomial coefficients from a supplied power table.
 * @see Deviation: PARI coefficient-substitution adapters
 */
export function FlxY_FlxqV_evalx(P: bigint[][], V: bigint[][], T: bigint[], p: bigint): bigint[][] {
  return coefficientSubstitution(1, 1, p, T, P, [], V) as bigint[][];
}

import { extensionAutomorphism } from './_extension_automorphism.js';

/** Native extension automorphism power.
 * @see Deviation: PARI extension automorphism adapters
 */
export function FlxqXQ_autpow(
  aut: [bigint[], bigint[][]],
  n: bigint,
  S: ExtensionModulus<bigint[][]>,
  T: bigint[],
  p: bigint
): [bigint[], bigint[][]] {
  return extensionAutomorphism(1, 0, aut, n, S, T, p) as [bigint[], bigint[][]];
}

/** Native additive trace.
 * @see Deviation: PARI extension automorphism adapters
 */
export function FlxqXQ_auttrace(
  aut: [bigint[][], bigint[][]],
  n: bigint,
  S: ExtensionModulus<bigint[][]>,
  T: bigint[],
  p: bigint
): [bigint[][], bigint[][]] {
  return extensionAutomorphism(1, 1, aut, n, S, T, p) as [bigint[][], bigint[][]];
}

/** Native multiplicative aggregate.
 * @see Deviation: PARI extension automorphism adapters
 */
export function FlxqXQ_autsum(
  aut: [bigint[], bigint[][], bigint[][]],
  n: bigint,
  S: ExtensionModulus<bigint[][]>,
  T: bigint[],
  p: bigint
): [bigint[], bigint[][], bigint[][]] {
  return extensionAutomorphism(1, 2, aut, n, S, T, p) as [bigint[], bigint[][], bigint[][]];
}


import { extensionProjection } from './_extension_projection.js';
import { extensionTruncatedPolynomial } from './_extension_polynomial.js';
/** Native random_FlxqX; ascending coefficient arrays, without variable metadata.
 * @see Deviation: PARI extension projection adapters
 */
export function random_FlxqX(length: number, T: bigint[], p: bigint): bigint[][] {
  return extensionProjection(1, 0, length, T, p) as bigint[][];
}
/** Native FlxqXn_mul; ascending coefficient arrays, without variable metadata.
 * @see Deviation: PARI extension projection adapters
 */
export function FlxqXn_mul(
  x: bigint[][],
  y: bigint[][],
  n: number,
  T: bigint[],
  p: bigint
): bigint[][] {
  return extensionTruncatedPolynomial(1, false, n, T, p, x, y) as bigint[][];
}
/** Native FlxqXn_sqr; ascending coefficient arrays, without variable metadata.
 * @see Deviation: PARI extension projection adapters
 */
export function FlxqXn_sqr(x: bigint[][], n: number, T: bigint[], p: bigint): bigint[][] {
  return extensionTruncatedPolynomial(1, true, n, T, p, x) as bigint[][];
}


import { extensionMinimalPolynomial } from './_extension_minpoly.js';
/** Native FlxqXQ_minpoly: Shoup projection over the extension coefficient field.
 * @see Deviation: PARI extension minimal-polynomial adapters
 */
export function FlxqXQ_minpoly(
  x: bigint[][],
  S: ExtensionModulus<bigint[][]>,
  T: bigint[],
  p: bigint
): bigint[][] {
  return extensionMinimalPolynomial(1, x, S, T, p) as bigint[][];
}

import { extensionDerivative } from './_extension_derivative.js';
/** Native outer derivative, retaining coefficient tags.
 * @see Deviation: PARI extension root-count and derivative adapters
 */
export function FlxX_deriv(x: bigint[][], p: bigint): bigint[][] {
  return extensionDerivative(true, x, p) as bigint[][];
}
