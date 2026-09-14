/** PARI Flv.c word-matrix kernel. Both matrix axes have unused slot zero. */
import { fieldKernel } from './_matrix_kernel.js';
import { fromColumns } from './_matrix_mul.js';
import { residue } from './_polynomial_division.js';

export function Flm_ker(x: bigint[][], p: bigint): bigint[][] {
  const a = fromColumns(x).map((row) => row.map((v) => residue(v, p)));
  return [
    [],
    ...fieldKernel(a, a.length, Math.max(0, x.length - 1), p, true).map((v) => [0n, ...v]),
  ];
}

import { wordMatrixPivots, wordMatrixSolve } from './_matrix_inverse.js';
import { toColumns } from './_matrix_mul.js';
/** Native word-matrix pivot rows/nullity; matrices and nonnull pivot vectors have unused slot zero.
 * @see Deviation: PARI word matrix pivot and solve adapters
 */
export function Flm_pivots(x: bigint[][], p: bigint): [number[] | null, number] {
  const a = wordColumns(x, p),
    n = Math.max(0, x.length - 1);
  const [d, r] = wordMatrixPivots(a, a.length, n, p);
  return [d === null ? null : [0, ...d], r];
}
/** Native word-matrix solve; both axes have unused slot zero, null on rank failure.
 * @see Deviation: PARI word matrix pivot and solve adapters
 */
export function Flm_gauss(a: bigint[][], b: bigint[][], p: bigint): bigint[][] | null {
  const A = wordColumns(a, p),
    n = Math.max(0, a.length - 1),
    k = Math.max(0, b.length - 1);
  const B = k ? wordColumns(b, p) : A.map(() => []);
  const u = wordMatrixSolve(A, B, A.length, n, k, p);
  return u === null ? null : !n ? [[]] : toColumns(u, k);
}

function wordColumns(x: bigint[][], p: bigint): bigint[][] {
  if (p <= 1n || p >= 1n << 64n)
    throw new RangeError('word matrix operations require a prime in the unsigned 64-bit range');
  if (x.length > 1 && (!x[1]!.length || x.slice(1).some((c) => c.length !== x[1]!.length)))
    throw new RangeError('word matrices require rectangular columns with unused slot zero');
  return fromColumns(x);
}
