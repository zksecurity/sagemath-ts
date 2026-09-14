/** PARI basemath/F2x.c, with packed bigint coefficient bits.
 * @see Deviation: PARI packed binary-polynomial kernels
 */
import { PariError } from './errors.js';
const WORD_MASK = (1n << 64n) - 1n;
function requirePolynomial(x: bigint): void {
  if (x < 0n) throw new RangeError('polynomial bits must be nonnegative');
}
export function F2x_degree(x: bigint): number {
  requirePolynomial(x);
  return x === 0n ? -1 : x.toString(2).length - 1;
}
export function F2x_add(x: bigint, y: bigint): bigint {
  requirePolynomial(x);
  requirePolynomial(y);
  return x ^ y;
}
/** Native 64-bit GMP cutoffs: 11 words for Karatsuba, 41 for integer packing. */
export function F2x_mul(x: bigint, y: bigint): bigint {
  requirePolynomial(x);
  requirePolynomial(y);
  if (!x || !y) return 0n;
  let shift = 0n;
  while (!(x & WORD_MASK)) {
    x >>= 64n;
    shift += 64n;
  }
  while (!(y & WORD_MASK)) {
    y >>= 64n;
    shift += 64n;
  }
  let nx = Math.ceil((F2x_degree(x) + 1) / 64),
    ny = Math.ceil((F2x_degree(y) + 1) / 64);
  if (nx < ny) {
    [x, y] = [y, x];
    [nx, ny] = [ny, nx];
  }
  if (nx < 11) {
    let out = 0n;
    while (y) {
      if (y & 1n) out ^= x;
      y >>= 1n;
      x <<= 1n;
    }
    return out << shift;
  }
  if (ny >= 41) {
    const dx = F2x_degree(x),
      dy = F2x_degree(y),
      bs = (1 + Math.max(dx, dy)).toString(2).length;
    const zeros = '0'.repeat(bs - 1);
    const encode = (f: bigint) => BigInt('0b' + f.toString(2).split('').join(zeros));
    const product = (encode(x) * encode(y)).toString(2);
    const bits = Array.from(
      { length: dx + dy + 1 },
      (_, i) => product[product.length - 1 - (dx + dy - i) * bs] ?? '0'
    ).join('');
    return BigInt('0b' + bits) << shift;
  }
  const width = BigInt(Math.ceil(nx / 2) * 64),
    mask = (1n << width) - 1n;
  const lo = x & mask,
    hi = x >> width;
  if (ny <= Math.ceil(nx / 2)) return (F2x_mul(lo, y) ^ (F2x_mul(hi, y) << width)) << shift;
  const ylo = y & mask,
    yhi = y >> width,
    a = F2x_mul(lo, ylo),
    b = F2x_mul(hi, yhi);
  const middle = F2x_mul(lo ^ hi, ylo ^ yhi) ^ a ^ b;
  return (a ^ (middle << width) ^ (b << (2n * width))) << shift;
}
export function F2x_sqr(x: bigint): bigint {
  requirePolynomial(x);
  const table = [
    '00',
    '01',
    '04',
    '05',
    '10',
    '11',
    '14',
    '15',
    '40',
    '41',
    '44',
    '45',
    '50',
    '51',
    '54',
    '55',
  ];
  return BigInt('0x' + Array.from(x.toString(16), (c) => table[Number.parseInt(c, 16)]).join(''));
}
/** Native sqrt assumes a square; nonsquares would index outside its lookup table. */
export function F2x_sqrt(x: bigint): bigint {
  requirePolynomial(x);
  if (F2x_deriv(x) !== 0n) throw new RangeError('polynomial must be a square');
  let hex = x.toString(16);
  if (hex.length % 2) hex = '0' + hex;
  const out: string[] = [];
  for (let i = 0; i < hex.length; i += 2) {
    const b = Number.parseInt(hex.slice(i, i + 2), 16);
    out.push(((b & 1) | ((b >> 1) & 2) | ((b >> 2) & 4) | ((b >> 3) & 8)).toString(16));
  }
  return BigInt('0x' + out.join(''));
}
export function F2x_rem(x: bigint, y: bigint): bigint {
  requirePolynomial(x);
  const dy = F2x_degree(y);
  if (dy < 0) throw new RangeError('polynomial divisor must be nonzero');
  if (dy === 0) return 0n;
  for (let dx = F2x_degree(x); dx >= dy; dx = F2x_degree(x)) x ^= y << BigInt(dx - dy);
  return x;
}
export function F2x_divrem(x: bigint, y: bigint): [bigint, bigint] {
  requirePolynomial(x);
  const dy = F2x_degree(y);
  if (dy < 0) throw new PariError('impossible inverse in F2x_divrem: Vecsmall([0]).');
  if (dy === 0) return [x, 0n];
  let quotient = 0n;
  for (let dx = F2x_degree(x); dx >= dy; dx = F2x_degree(x)) {
    const shift = BigInt(dx - dy);
    quotient |= 1n << shift;
    x ^= y << shift;
  }
  return [quotient, x];
}
export function F2x_gcd(a: bigint, b: bigint): bigint {
  requirePolynomial(a);
  requirePolynomial(b);
  if (F2x_degree(b) > F2x_degree(a)) [a, b] = [b, a];
  while (b) [a, b] = [b, F2x_rem(a, b)];
  return a;
}
export function F2x_deriv(x: bigint): bigint {
  requirePolynomial(x);
  return BigInt(
    '0x' +
      Array.from(x.toString(16), (c) => ((Number.parseInt(c, 16) >> 1) & 5).toString(16)).join('')
  );
}
/** Returns [valuation, quotient]; zero has native LONG_MAX valuation. */
export function F2x_valrem(x: bigint): [bigint, bigint] {
  requirePolynomial(x);
  if (!x) return [(1n << 63n) - 1n, 0n];
  const bits = x.toString(2),
    v = bits.length - 1 - bits.lastIndexOf('1');
  return [BigInt(v), x >> BigInt(v)];
}
export function F2xq_mul(x: bigint, y: bigint, f: bigint): bigint {
  return F2x_rem(F2x_mul(x, y), f);
}
export function F2xq_sqr(x: bigint, f: bigint): bigint {
  return F2x_rem(F2x_sqr(x), f);
}
/** Includes the zeroth power and the original first power, as gen_powers does. */
export function F2xq_powers(x: bigint, l: number, f: bigint): bigint[] {
  requirePolynomial(x);
  requirePolynomial(f);
  if (!Number.isSafeInteger(l) || l < 0) throw new RangeError('power count must be nonnegative');
  const out = [1n];
  if (!l) return out;
  out.push(x);
  if (l === 1) return out;
  out.push(F2xq_sqr(x, f));
  const useSquare = 2 * F2x_degree(x) >= F2x_degree(f);
  for (let i = 3; i <= l; i++)
    out.push(useSquare && i % 2 === 0 ? F2xq_sqr(out[i / 2]!, f) : F2xq_mul(out[i - 1]!, x, f));
  return out;
}
export function F2x_Frobenius(f: bigint): bigint {
  return F2xq_sqr(2n, f);
}
/** Packed columns of the native matrix of Frobenius. */
export function F2x_matFrobenius(f: bigint): bigint[] {
  const degree = F2x_degree(f);
  if (degree <= 0) throw new RangeError('modulus must have positive degree');
  return F2xq_powers(F2x_Frobenius(f), degree - 1, f);
}

import {extensionPolynomial} from './_extension_polynomial.js';
/** Native F2xqX_mul; each outer coefficient is a packed bit polynomial.
 * @see Deviation: PARI extension-polynomial coefficient adapters
 */
export function F2xqX_mul(x:bigint[], y:bigint[], T:bigint):bigint[]{
 return extensionPolynomial(2,0,2n,T,x,y) as bigint[];
}
/** Native F2xqX_sqr; each outer coefficient is a packed bit polynomial.
 * @see Deviation: PARI extension-polynomial coefficient adapters
 */
export function F2xqX_sqr(x:bigint[], T:bigint):bigint[]{
 return extensionPolynomial(2,1,2n,T,x) as bigint[];
}
/** Native F2xqX_red; each outer coefficient is a packed bit polynomial.
 * @see Deviation: PARI extension-polynomial coefficient adapters
 */
export function F2xqX_red(x:bigint[], T:bigint):bigint[]{
 return extensionPolynomial(2,2,2n,T,x) as bigint[];
}
/** Native F2xqX_normalize; each outer coefficient is a packed bit polynomial.
 * @see Deviation: PARI extension-polynomial coefficient adapters
 */
export function F2xqX_normalize(x:bigint[], T:bigint):bigint[]{
 return extensionPolynomial(2,3,2n,T,x) as bigint[];
}

import {extensionDivision,extensionGetRed,type ExtensionModulus} from './_extension_division.js';
/** Native F2xqX_divrem, with PARI's coefficient and dispatch rules.
 * @see Deviation: PARI extension-polynomial division adapters
 */
export function F2xqX_divrem(x:bigint[], S:ExtensionModulus<bigint[]>, T:bigint):[bigint[],bigint[]]{
 return extensionDivision(2,0,2n,T,x,S) as [bigint[],bigint[]];
}
/** Native F2xqX_rem, with PARI's coefficient and dispatch rules.
 * @see Deviation: PARI extension-polynomial division adapters
 */
export function F2xqX_rem(x:bigint[], S:ExtensionModulus<bigint[]>, T:bigint):bigint[]{
 return extensionDivision(2,1,2n,T,x,S) as bigint[];
}
/** Native F2xqX_div, with PARI's coefficient and dispatch rules.
 * @see Deviation: PARI extension-polynomial division adapters
 */
export function F2xqX_div(x:bigint[], S:ExtensionModulus<bigint[]>, T:bigint):bigint[]{
 return extensionDivision(2,4,2n,T,x,S) as bigint[];
}
/** Native F2xqX_invBarrett, with PARI's coefficient and dispatch rules.
 * @see Deviation: PARI extension-polynomial division adapters
 */
export function F2xqX_invBarrett(x:bigint[], T:bigint):bigint[]{
 return extensionDivision(2,2,2n,T,x) as bigint[];
}
/** Native cached reciprocal selection; existing reduction objects are retained.
 * @see Deviation: PARI extension-polynomial division adapters
 */
export function F2xqX_get_red(S:ExtensionModulus<bigint[]>,T:bigint):ExtensionModulus<bigint[]>{
 return extensionGetRed(2,S,T,2n);
}

import { extensionGcd, type ExtensionMatrix } from './_extension_gcd.js';
/** Native unscaled greatest common divisor.
 * @see Deviation: PARI extension-polynomial GCD adapters
 */
export function F2xqX_gcd(x:bigint[], y:bigint[], T:bigint):bigint[]{
  return extensionGcd(2,0,2n,T,x,y) as bigint[];
}
/** Native [gcd, U, V] with U*x + V*y = gcd.
 * @see Deviation: PARI extension-polynomial GCD adapters
 */
export function F2xqX_extgcd(x:bigint[], y:bigint[], T:bigint):[bigint[],bigint[],bigint[]]{
  return extensionGcd(2,1,2n,T,x,y) as [bigint[],bigint[],bigint[]];
}
/** Native half-GCD transformation matrix in row-major order.
 * @see Deviation: PARI extension-polynomial GCD adapters
 */
export function F2xqX_halfgcd(x:bigint[], y:bigint[], T:bigint):ExtensionMatrix<bigint[]>{
  return extensionGcd(2,2,2n,T,x,y) as ExtensionMatrix<bigint[]>;
}

import {extensionQuotient} from './_extension_quotient.js';
/** Product in the extension quotient.
 * @see Deviation: PARI extension quotient adapters
 */
export function F2xqXQ_mul(x:bigint[], y:bigint[], S:ExtensionModulus<bigint[]>, T:bigint):bigint[]{
  return extensionQuotient(2,0,2n,0n,T,x,y,S)as bigint[];
}
/** Square in the extension quotient.
 * @see Deviation: PARI extension quotient adapters
 */
export function F2xqXQ_sqr(x:bigint[], S:ExtensionModulus<bigint[]>, T:bigint):bigint[]{
  return extensionQuotient(2,1,2n,0n,T,x,[],S)as bigint[];
}
/** Inverse or null when no inverse exists.
 * @see Deviation: PARI extension quotient adapters
 */
export function F2xqXQ_invsafe(x:bigint[], S:ExtensionModulus<bigint[]>, T:bigint):bigint[]|null{
  return extensionQuotient(2,2,2n,0n,T,x,[],S)as bigint[]|null;
}
/** Native quotient inverse, with exact GEN error display.
 * @see Deviation: PARI extension quotient adapters
 */
export function F2xqXQ_inv(x:bigint[], S:ExtensionModulus<bigint[]>, T:bigint):bigint[]{
  return extensionQuotient(2,3,2n,0n,T,x,[],S)as bigint[];
}
/** Native signed-integer powering schedule.
 * @see Deviation: PARI extension quotient adapters
 */
export function F2xqXQ_pow(x:bigint[], n:bigint, S:ExtensionModulus<bigint[]>, T:bigint):bigint[]{
  return extensionQuotient(2,5,2n,n,T,x,[],S)as bigint[];
}
/** Powers from zero through l, retaining the raw first power.
 * @see Deviation: PARI extension quotient adapters
 */
export function F2xqXQ_powers(x:bigint[], l:number, S:ExtensionModulus<bigint[]>, T:bigint):bigint[][]{
  if(!Number.isSafeInteger(l)||l<0||l>=0xffffffff)throw new RangeError('power count must be a nonnegative array length');
  return extensionQuotient(2,6,2n,BigInt(l),T,x,[],S)as bigint[][];
}

import {extensionComposition} from './_extension_composition.js';
/** Native extension-polynomial composition.
 * @see Deviation: PARI extension composition adapters
 */
export function F2xqX_F2xqXQ_eval(Q:bigint[], x:bigint[], S:ExtensionModulus<bigint[]>, T:bigint):bigint[]{
 return extensionComposition(2,0,2n,T,Q,x,[],S)as bigint[];
}
/** Native extension-polynomial composition from a supplied power table.
 * @see Deviation: PARI extension composition adapters
 */
export function F2xqX_F2xqXQV_eval(Q:bigint[], V:bigint[][], S:ExtensionModulus<bigint[]>, T:bigint):bigint[]{
 return extensionComposition(2,1,2n,T,Q,[],V,S)as bigint[];
}

import { coefficientSubstitution } from './_coefficient_composition.js';
import { binaryComposition } from './_binary_composition.js';

/** Native substitution in polynomial coefficients.
 * @see Deviation: PARI coefficient-substitution adapters
 */
export function F2xY_F2xq_evalx(P: bigint[], x: bigint, T: bigint): bigint[] {
  return coefficientSubstitution(2, 0, 2n, T, P, x, []) as bigint[];
}

/** Native substitution in polynomial coefficients from a supplied power table.
 * @see Deviation: PARI coefficient-substitution adapters
 */
export function F2xY_F2xqV_evalx(P: bigint[], V: bigint[], T: bigint): bigint[] {
  return coefficientSubstitution(2, 1, 2n, T, P, 0n, V) as bigint[];
}

/** Native binary scalar composition.
 * @see Deviation: PARI packed binary-polynomial kernels
 */
export function F2x_F2xq_eval(Q: bigint, x: bigint, T: bigint): bigint {
  return binaryComposition(Q, x, [], T, true);
}

/** Native binary scalar composition from a supplied power table.
 * @see Deviation: PARI packed binary-polynomial kernels
 */
export function F2x_F2xqV_eval(Q: bigint, V: bigint[], T: bigint): bigint {
  return binaryComposition(Q, 0n, V, T, false);
}

import { extensionAutomorphism } from './_extension_automorphism.js';

/** Native extension automorphism power.
 * @see Deviation: PARI extension automorphism adapters
 */
export function F2xqXQ_autpow(
  aut: [bigint, bigint[]],
  n: bigint,
  S: ExtensionModulus<bigint[]>,
  T: bigint
): [bigint, bigint[]] {
  return extensionAutomorphism(2, 0, aut, n, S, T, 2n) as [bigint, bigint[]];
}

/** Native additive trace.
 * @see Deviation: PARI extension automorphism adapters
 */
export function F2xqXQ_auttrace(
  aut: [bigint, bigint[], bigint[]],
  n: bigint,
  S: ExtensionModulus<bigint[]>,
  T: bigint
): [bigint, bigint[], bigint[]] {
  return extensionAutomorphism(2, 1, aut, n, S, T, 2n) as [bigint, bigint[], bigint[]];
}
