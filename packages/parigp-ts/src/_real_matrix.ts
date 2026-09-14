/** Numeric branches of PARI gen1.c used by inexact matrix dependencies.
 * @see Deviation: PARI real Cholesky and triangular inverse adapters
 */
import { type MpReal, mulrr, subrr, negr, divrr } from './qfb.js';
import { invr } from './kernel/none/mp_indep.js';
import { PariError } from './errors.js';
/** Exact gen_0 remains distinct from a finite-accuracy real zero. */
export type MatrixReal = MpReal | 0n;
export const matrixRealMul = (x: MatrixReal, y: MatrixReal): MatrixReal =>
  x === 0n || y === 0n ? 0n : mulrr(x, y);
export const matrixRealNeg = (x: MatrixReal): MatrixReal => (x === 0n ? 0n : negr(x));
export const matrixRealSub = (x: MatrixReal, y: MatrixReal): MatrixReal =>
  y === 0n ? x : x === 0n ? negr(y) : subrr(x, y);
function zeroError(name: string, y: MatrixReal): never {
  const zero =
    y === 0n ? '0' : '0.E' + (Math.trunc(y.e * 0.30102999566398119521 - (y.e < 0 ? 1 : 0)) + 1);
  throw new PariError(`impossible inverse in ${name}: ${zero}`);
}
export function matrixRealInv(x: MatrixReal): MpReal {
  if (x === 0n) zeroError('ginv', x);
  if (!x.s) {
    if (!x.p) zeroError('invr', x);
    if (x.p > 4800) return invr(x);
    zeroError('divrr', x);
  }
  return invr(x);
}
export function matrixRealDiv(x: MatrixReal, y: MatrixReal): MatrixReal {
  if (y === 0n || !y.s)
    zeroError(x === 0n && y === 0n ? 'dvmdii' : x === 0n || y === 0n ? 'gdiv' : 'divrr', y);
  return x === 0n ? 0n : divrr(x, y);
}
