/**
 * PARI algebraic dependencies for native double inputs.
 * Reference: reference/pari/src/basemath/bibli1.c:algdep0/lindep2/lindepfull_bit.
 * @see Deviation: PARI algebraic dependencies
 */
import { dbltor, trunc2nr } from './kernel/none/mp_indep.js';
import { mulrr } from './qfb.js';
import { ZM_lll, LLL_INPLACE } from './lll.js';

/** Return ascending coefficients of PARI's algdep polynomial for a double input. */
export function algdep(z: number, degree: bigint): bigint[] {
  // cypari2 converts the value before the C-long degree argument.
  const x = dbltor(z);
  if (degree < -(1n << 63n) || degree >= 1n << 63n) {
    const error = new Error('Python int too large to convert to C long');
    error.name = 'OverflowError';
    throw error;
  }
  if (!x.s) return [0n, 1n];
  if (degree < 0n) {
    const error = new Error('domain error in algdep: degree < 0');
    error.name = 'PariError';
    throw error;
  }
  if (degree === 0n) return [1n];
  // dbltor creates one 64-bit real word. lindep2 uses floor(0.8 * 64) bits.
  const bit = 51;
  const scaled = [1n << BigInt(bit)];
  let power = x;
  for (let i = 1n; i <= degree; i++) {
    scaled.push(trunc2nr(power, bit));
    // This is PARI real multiplication, with a 64-bit mantissa at every step.
    if (i < degree) power = mulrr(power, x);
  }
  // lindepfull_bit reduces [e_i, scaled_i] directly. Gram-mode reduction can
  // choose another valid relation because it follows a different LLL path.
  const basis = scaled.map((b, j) => [...scaled.map((_, i) => (i === j ? 1n : 0n)), b]);
  const reduced = ZM_lll(basis, 0.99, LLL_INPLACE) as bigint[][];
  let coefficients = reduced[0]!.slice(0, scaled.length);
  while (coefficients.length > 1 && coefficients[coefficients.length - 1] === 0n)
    coefficients.pop();
  if (coefficients[coefficients.length - 1]! < 0n) coefficients = coefficients.map((c) => -c);
  return coefficients;
}

import {
  type MpReal,
  addir,
  addrr,
  subir,
  subrr,
  mulir,
  sqrr,
  sqrtr_abs,
  itor,
  rtor,
  negr,
} from './qfb.js';
import { invr } from './kernel/none/mp_indep.js';

export type QrScalar = bigint | MpReal;
export type Householder = [MpReal, QrScalar[]];
const qrSign = (x: QrScalar): number =>
  typeof x === 'bigint' ? (x < 0n ? -1 : x > 0n ? 1 : 0) : x.s;
const qrAdd = (x: QrScalar, y: QrScalar): QrScalar =>
  typeof x === 'bigint'
    ? typeof y === 'bigint'
      ? x + y
      : addir(x, y)
    : typeof y === 'bigint'
      ? addir(y, x)
      : addrr(x, y);
const qrSub = (x: QrScalar, y: QrScalar): QrScalar =>
  typeof x === 'bigint'
    ? typeof y === 'bigint'
      ? x - y
      : subir(x, y)
    : typeof y === 'bigint'
      ? addir(-y, x)
      : subrr(x, y);
const qrMul = (x: QrScalar, y: QrScalar): QrScalar =>
  typeof x === 'bigint'
    ? typeof y === 'bigint'
      ? x * y
      : mulir(x, y)
    : typeof y === 'bigint'
      ? mulir(y, x)
      : mulrr(x, y);
const qrSqr = (x: QrScalar): QrScalar => (typeof x === 'bigint' ? x * x : sqrr(x));

function ApplyQ([beta, v]: Householder, r: QrScalar[]): void {
  const offset = r.length - v.length;
  let s = qrMul(v[0]!, r[offset]!);
  for (let i = 1; i < v.length; i++) s = qrAdd(s, qrMul(v[i]!, r[offset + i]!));
  s = qrMul(beta, s);
  for (let i = 0; i < v.length; i++)
    if (qrSign(v[i]!)) r[offset + i] = qrSub(r[offset + i]!, qrMul(s, v[i]!));
}

/** Native Householder QR; zero-indexed columns, integer/real cells, bit precision.
 * Returns [success, squared norms, reflectors, lower triangular L].
 * @see Deviation: PARI Householder QR adapters
 */
export function QR_init(
  x: QrScalar[][],
  precision = 64
): [number, QrScalar[] | null, Householder[] | null, QrScalar[][] | null] {
  const n = x.length,
    rows = x[0]?.length ?? 0;
  if (!n || rows < n || x.some((c) => c.length !== rows))
    throw new RangeError('QR_init requires nonempty rectangular columns with rows >= columns');
  if (!Number.isSafeInteger(precision) || precision <= 0 || precision % 64)
    throw new RangeError('QR_init precision must be a positive multiple of 64 bits');
  const B: QrScalar[] = [],
    Q: Householder[] = [];
  const L: QrScalar[][] = Array.from({ length: n }, () => Array<QrScalar>(n).fill(0n));
  for (let k = 0; k < n; k++) {
    const r = [...x[k]!];
    for (let j = 0; j < k; j++) ApplyQ(Q[j]!, r);
    const x1 = r[k]!;
    let x2 = qrSqr(x1);
    if (k < rows - 1) {
      const v = r.slice(k);
      for (let i = k + 1; i < rows; i++) x2 = qrAdd(x2, qrSqr(r[i]!));
      if (!qrSign(x2)) return [0, null, null, null];
      let Nx = sqrtr_abs(typeof x2 === 'bigint' ? itor(x2, precision) : x2);
      if (qrSign(x1) < 0) Nx = negr(Nx);
      v[0] = qrAdd(x1, Nx);
      const beta = !qrSign(x1)
        ? typeof x2 === 'bigint'
          ? itor(x2, precision)
          : rtor(x2, precision)
        : (qrAdd(x2, qrMul(Nx, x1)) as MpReal);
      Q[k] = [invr(beta as MpReal), v];
      L[k]![k] = negr(Nx);
    } else L[k]![k] = x1;
    B[k] = x2;
    for (let i = 0; i < k; i++) L[i]![k] = r[i]!;
    if (typeof x2 !== 'bigint' && x2.p <= 64 && x2.e >= 32) return [0, null, null, null];
  }
  // Native Q has n-1 public entries even for a tall matrix. The additional
  // last reflector is internal scratch and is not needed by later columns.
  return [1, B, Q.slice(0, n - 1), L];
}

/** Native upper R from Householder QR, or null on precision failure.
 * @see Deviation: PARI Householder QR adapters
 */
export function R_from_QR(x: QrScalar[][], precision = 64): QrScalar[][] | null {
  const [, , , L] = QR_init(x, precision);
  return L === null ? null : L.map((_, j) => L.map((c) => c[j]!));
}

/** Native Gaussian reduction of x^T*x through QR; x must be square and full rank.
 * @see Deviation: PARI Householder QR adapters
 */
export function gaussred_from_QR(x: QrScalar[][], precision = 64): QrScalar[][] | null {
  if (x.some((c) => c.length !== x.length))
    throw new RangeError('gaussred_from_QR requires square columns');
  const [, B, , L] = QR_init(x, precision);
  if (L === null) return null;
  for (let j = 0; j < L.length - 1; j++) {
    const m = L[j]!,
      invNx = invr(m[j] as MpReal);
    m[j] = B![j]!;
    for (let i = j + 1; i < L.length; i++) m[i] = qrMul(invNx, m[i]!);
  }
  L[L.length - 1]![L.length - 1] = B![L.length - 1]!;
  return L.map((_, j) => L.map((c) => c[j]!));
}
