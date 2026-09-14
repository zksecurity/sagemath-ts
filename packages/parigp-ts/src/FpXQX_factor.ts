/** PARI basemath/FpXQX_factor.c: finite extension factorization dependencies. */
import { extensionFrobenius } from './_extension_frobenius.js';
import { type ExtensionPolynomial as P } from './_extension_polynomial.js';
import { type ExtensionModulus as M } from './_extension_division.js';
/** Native FpXQX_Frobenius.
 * @see Deviation: PARI extension Frobenius adapters
 */
export function FpXQX_Frobenius(S: M<P>, T: bigint[], p: bigint): P {
  return extensionFrobenius(0, false, [], S, T, p);
}
/** Native FlxqX_Frobenius.
 * @see Deviation: PARI extension Frobenius adapters
 */
export function FlxqX_Frobenius(S: M<bigint[][]>, T: bigint[], p: bigint): bigint[][] {
  return extensionFrobenius(1, false, [], S, T, p) as bigint[][];
}
/** Native F2xqX_Frobenius.
 * @see Deviation: PARI extension Frobenius adapters
 */
export function F2xqX_Frobenius(S: M<bigint[]>, T: bigint): bigint[] {
  return extensionFrobenius(2, false, [], S, T, 2n) as bigint[];
}
/** Native FpXQXQ_halfFrobenius.
 * @see Deviation: PARI extension Frobenius adapters
 */
export function FpXQXQ_halfFrobenius(a: P, S: M<P>, T: bigint[], p: bigint): P {
  return extensionFrobenius(0, true, a, S, T, p);
}
/** Native FlxqXQ_halfFrobenius.
 * @see Deviation: PARI extension Frobenius adapters
 */
export function FlxqXQ_halfFrobenius(
  a: bigint[][],
  S: M<bigint[][]>,
  T: bigint[],
  p: bigint
): bigint[][] {
  return extensionFrobenius(1, true, a, S, T, p) as bigint[][];
}

import { extensionSplitPart } from './_extension_roots.js';
import { extensionGcd } from './_extension_gcd.js';
import { FpXX_deriv } from './FpXX.js';
import { FlxX_deriv } from './FlxX.js';
import { FpX_nbroots } from './galconj.js';
/** Native FpXQX_split_part.
 * @see Deviation: PARI extension root-count and derivative adapters
 */
export function FpXQX_split_part(f: P, T: bigint[], p: bigint): P {
  return extensionSplitPart(0, f, T, p);
}
/** Native FpXQX_nbroots.
 * @see Deviation: PARI extension root-count and derivative adapters
 */
export function FpXQX_nbroots(f: P, T: bigint[], p: bigint): number {
  return extensionSplitPart(0, f, T, p).length - 1;
}
/** Native FlxqX_nbroots.
 * @see Deviation: PARI extension root-count and derivative adapters
 */
export function FlxqX_nbroots(f: bigint[][], T: bigint[], p: bigint): number {
  return extensionSplitPart(1, f, T, p).length - 1;
}
/** Native F2xqX_nbroots.
 * @see Deviation: PARI extension root-count and derivative adapters
 */
export function F2xqX_nbroots(f: bigint[], T: bigint): number {
  return extensionSplitPart(2, f, T, 2n).length - 1;
}
/** Native FqX_nbroots.
 * @see Deviation: PARI extension root-count and derivative adapters
 */
export function FqX_nbroots(f: P, T: bigint[] | null, p: bigint): number {
  if (T !== null) return FpXQX_nbroots(f, T, p);
  if (f.some((c) => typeof c !== 'bigint'))
    throw new TypeError('prime-field polynomial coefficients must be integers');
  return FpX_nbroots(f as bigint[], p);
}
/** Native FpXQX_is_squarefree.
 * @see Deviation: PARI extension root-count and derivative adapters
 */
export function FpXQX_is_squarefree(f: P, T: bigint[], p: bigint): boolean {
  const derivative = FpXX_deriv(f, p);
  return (extensionGcd(0, 0, p, T, f, derivative) as P).length === 1;
}
/** Native FlxqX_is_squarefree.
 * @see Deviation: PARI extension root-count and derivative adapters
 */
export function FlxqX_is_squarefree(f: bigint[][], T: bigint[], p: bigint): boolean {
  const derivative = FlxX_deriv(f, p);
  return (extensionGcd(1, 0, p, T, f, derivative) as P).length === 1;
}
