import { PariError } from './errors.js';
import { type MpReal, mulrr, negr, real_1 } from './qfb.js';
import { addrr } from './kernel/none/add.js';
import { divrr } from './kernel/gmp/mp.js';

type RealCell = MpReal | null;
// null is PARI's exact gen_0, distinct from a t_REAL zero with finite accuracy.
const mul = (a: RealCell, b: RealCell): RealCell => (a === null || b === null ? null : mulrr(a, b));
const sub = (a: RealCell, b: RealCell): RealCell =>
  b === null ? a : a === null ? negr(b) : addrr(a, negr(b));
const exponent = (a: RealCell): number => (a === null ? -Infinity : a.e);

/**
 * PARI det/RgM_det2/det_simple_gauss for single-word real matrices, in row order.
 * Null cells model exact integer zeros used to pad the Sylvester matrix.
 * @see Deviation: Native Single-Word Real Elimination
 */
export function det(a: readonly (readonly (MpReal | null)[])[]): MpReal | null {
  const n = a.length;
  if (a.some((row) => row.length !== n)) throw new RangeError('inconsistent dimensions in det');
  if (!n) return real_1(64);
  if (n === 1) return a[0]![0]!;
  if (n === 2) return sub(mul(a[0]![0]!, a[1]![1]!), mul(a[0]![1]!, a[1]![0]!));
  const work = a.map((row) => [...row]);
  let product: MpReal | undefined,
    sign = 1;
  for (let i = 0; i < n - 1; i++) {
    // gauss_get_pivot_max uses exponent alone, preserving the first tie.
    let k = -1,
      e = -Infinity;
    for (let j = i; j < n; j++)
      if (exponent(work[j]![i]!) > e) {
        e = exponent(work[j]![i]!);
        k = j;
      }
    if (k < 0) return work[i]![i]!;
    const pivot = work[k]![i]!;
    const reference = a[k]![i]!;
    const originalExponent =
      reference === null ? Math.max(...a.map((row) => exponent(row[i]!))) : exponent(reference);
    // cx_approx0 compares the evolving pivot to the original input column.
    if (pivot === null || !pivot.s || originalExponent - pivot.e > pivot.p) return work[i]![i]!;
    if (k !== i) {
      for (let j = i; j < n; j++) [work[i]![j], work[k]![j]] = [work[k]![j]!, work[i]![j]!];
      sign = -sign;
    }
    const p = work[i]![i]! as MpReal;
    product = product === undefined ? p : mulrr(product, p);
    // Native elimination updates columns, leaving the eliminated cells alone.
    for (let col = i + 1; col < n; col++) {
      const entry = work[i]![col]!;
      if (entry === null || !entry.s) continue;
      const ratio = divrr(entry, p);
      for (let row = i + 1; row < n; row++)
        work[row]![col] = sub(work[row]![col]!, mul(ratio, work[row]![i]!));
    }
  }
  return mul(sign < 0 ? negr(product!) : product!, work[n - 1]![n - 1]!);
}

import { wordExtensionMatrix } from './_extension_matrix.js';
import { validateExtensionInputs } from './_extension_polynomial.js';
/** Word extension matrix product; columns and rows have unused slot zero.
 * @see Deviation: PARI extension composition adapters
 */
export function FlxqM_mul(A: bigint[][][], B: bigint[][][], T: bigint[], p: bigint): bigint[][][] {
  const a = A.slice(1).map((c) => c.slice(1)),
    b = B.slice(1).map((c) => c.slice(1));
  validateExtensionInputs(1, p, T, []);
  for (const c of [...a, ...b]) validateExtensionInputs(1, p, T, c);
  if (a.some((c) => c.length !== a[0]!.length) || b.some((c) => c.length !== a.length))
    throw new RangeError('extension matrices must be rectangular and dimension-compatible');
  return [[], ...wordExtensionMatrix(a, b, T, p).map((c) => [[], ...c])];
}

import { fieldKernel } from './_matrix_kernel.js';
import { fromColumns } from './_matrix_mul.js';
import { residue } from './_polynomial_division.js';
import { F2m_ker } from './F2v.js';
import { F3m_ker } from './F3v.js';
import { Flm_ker } from './Flv.js';

/** Native FpM_ker dispatch; columns and rows have unused slot zero. */
export function FpM_ker(x: bigint[][], p: bigint): bigint[][] {
  if (x.length <= 1) return [[]];
  const rows = x[1]!.length - 1,
    n = x.length - 1;
  if (p === 2n || p === 3n) {
    const bits = p === 2n ? 1n : 2n;
    const packed = x
      .slice(1)
      .map((c) =>
        c.slice(1).reduce((v, a, i) => v | (residue(a, p) << BigInt(Number(bits) * i)), 0n)
      );
    const kernel = p === 2n ? F2m_ker(packed, rows) : F3m_ker(packed, rows);
    return [
      [],
      ...kernel.map((v) => [
        0n,
        ...Array.from(
          { length: n },
          (_, i) => (v >> BigInt(Number(bits) * i)) & ((1n << bits) - 1n)
        ),
      ]),
    ];
  }
  if (p < 1n << 64n) return Flm_ker(x, p);
  return [[], ...fieldKernel(fromColumns(x), rows, n, p, false).map((v) => [0n, ...v])];
}

import {
  type MatrixReal,
  matrixRealInv,
  matrixRealMul,
  matrixRealNeg,
  matrixRealSub,
  matrixRealDiv,
} from './_real_matrix.js';
/** Native back substitution for an invertible upper-triangular real matrix.
 * Columns contain MpReal or exact integer zero cells.
 * @see Deviation: PARI real Cholesky and triangular inverse adapters
 */
export function RgM_inv_upper(A: MatrixReal[][]): MatrixReal[][] {
  const n = A.length;
  if (A.some((c) => c.length !== n)) throw new RangeError('RgM_inv_upper requires square columns');
  return A.map((_, index) => {
    const u = Array<MatrixReal>(n).fill(0n);
    u[index] = matrixRealInv(A[index]![index]!);
    for (let i = index - 1; i >= 0; i--) {
      let m = matrixRealNeg(matrixRealMul(A[i + 1]![i]!, u[i + 1]!));
      for (let j = i + 2; j < n; j++) m = matrixRealSub(m, matrixRealMul(A[j]![i]!, u[j]!));
      u[i] = matrixRealDiv(m, A[i]![i]!);
    }
    return u;
  });
}

import { integerMatrixPivots, integerMatrixSolve } from './_matrix_inverse.js';
import { toColumns as integerColumnsFromRows } from './_matrix_mul.js';
/** Native integer pivot profile/nullity; matrices and nonnull pivot vectors have unused slot zero.
 * @see Deviation: PARI modular integer rank and solve adapters
 */
export function ZM_pivots(x: bigint[][]): [number[] | null, number] {
  const a = checkedIntegerColumns(x),
    n = Math.max(0, x.length - 1);
  const [d, r] = integerMatrixPivots(a, a.length, n);
  return [d === null ? null : [0, ...d], r];
}
/** Native modular integer rank, including exact certification and native errors.
 * @see Deviation: PARI modular integer rank and solve adapters
 */
export function ZM_rank(x: bigint[][]): number {
  return Math.max(0, x.length - 1) - ZM_pivots(x)[1];
}
/** Native modular solve; return integer numerator columns/common denominator, or null on rank failure.
 * Both matrix axes have unused slot zero; requires rows >= columns.
 * @see Deviation: PARI modular integer rank and solve adapters
 */
export function ZM_gauss(a: bigint[][], b: bigint[][]): [bigint[][], bigint] | null {
  const A = checkedIntegerColumns(a),
    n = Math.max(0, a.length - 1),
    k = Math.max(0, b.length - 1);
  const B = k ? checkedIntegerColumns(b) : A.map(() => []);
  const result = integerMatrixSolve(A, B, A.length, n, k);
  return result === null ? null : [!n ? [[]] : integerColumnsFromRows(result[0], k), result[1]];
}
function checkedIntegerColumns(x: bigint[][]): bigint[][] {
  if (x.length > 1 && (!x[1]!.length || x.slice(1).some((c) => c.length !== x[1]!.length)))
    throw new RangeError('integer matrices require rectangular columns with unused slot zero');
  return fromColumns(x);
}

import { fieldMatrixPivots, fieldMatrixSolve } from './_matrix_kernel.js';
import { wordMatrixPivots, wordMatrixSolve } from './_matrix_inverse.js';
import { F2m_gauss_pivot, F2m_gauss } from './F2v.js';

function FpM_pivot_data(A: bigint[][], p: bigint): [number[] | null, number] {
  const n = A.length - 1,
    m = n > 0 ? A[1]!.length - 1 : 0;
  if (p === 2n)
    return F2m_gauss_pivot(
      A.slice(1).map((c) => c.slice(1).reduce((v, x, i) => v | ((x & 1n) << BigInt(i)), 0n)),
      m
    );
  return (p < 1n << 64n ? wordMatrixPivots : fieldMatrixPivots)(fromColumns(A), m, n, p);
}
/** Native FpM_image, retaining original columns at the native pivot indices.
 * @see Deviation: PARI prime-decomposition matrix adapters
 */
export function FpM_image(A: bigint[][], p: bigint): bigint[][] {
  const [d] = FpM_pivot_data(A, p);
  return [
    [],
    ...A.slice(1)
      .filter((_, i) => d?.[i])
      .map((c) => c.slice()),
  ];
}
/** Native FpM_suppl, retaining independent columns then appending missing unit columns.
 * @see Deviation: PARI prime-decomposition matrix adapters
 */
export function FpM_suppl(A: bigint[][], p: bigint): bigint[][] {
  if (A.length <= 1) throw new PariError('sorry, suppl [empty matrix] is not yet implemented');
  const [d] = FpM_pivot_data(A, p),
    m = A[1]!.length - 1,
    used = new Set(d!.filter(Boolean)),
    out = [
      [],
      ...A.slice(1)
        .filter((_, i) => d?.[i])
        .map((c) => c.slice()),
    ];
  for (let i = 1; i <= m; i++)
    if (!used.has(i))
      out.push([0n, ...Array.from({ length: m }, (_, j) => (j === i - 1 ? 1n : 0n))]);
  return out;
}
/** Native FpM_inv dispatches binary, word-CUP and generic-CUP solvers.
 * @see Deviation: PARI prime-decomposition matrix adapters
 */
export function FpM_inv(A: bigint[][], p: bigint): bigint[][] | null {
  if (A.length <= 1) return [[]];
  const m = A[1]!.length - 1,
    n = A.length - 1;
  if (p === 2n) {
    const x = A.slice(1).map((c) =>
        c.slice(1).reduce((v, x, i) => v | ((x & 1n) << BigInt(i)), 0n)
      ),
      result = F2m_gauss(
        x,
        Array.from({ length: m }, (_, i) => 1n << BigInt(i)),
        m
      );
    return result
      ? [
          [],
          ...result.map((c) => [0n, ...Array.from({ length: n }, (_, i) => (c >> BigInt(i)) & 1n)]),
        ]
      : null;
  }
  const I = Array.from({ length: m }, (_, i) =>
      Array.from({ length: m }, (_, j) => (i === j ? 1n : 0n))
    ),
    result = (p < 1n << 64n ? wordMatrixSolve : fieldMatrixSolve)(fromColumns(A), I, m, n, m, p);
  return result ? integerColumnsFromRows(result, m) : null;
}
