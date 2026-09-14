import { CRTInRange } from './ZZ.js';
export interface IntegerMatrixOptions {
  columns?: number;
}
export interface IntegerMatrixCRTOptions extends IntegerMatrixOptions {
  residueColumns?: number;
}
function integerMatrix(A: readonly (readonly bigint[])[], columns: number): bigint[][] {
  if (columns < 0) throw new Error('SetDims: bad args');
  if (A.some((row) => row.length !== columns)) throw new Error('nonrectangular matrix');
  return A.map((row) => [...row]);
}
function integerSqrt(n: bigint): bigint {
  if (n < 2n) return n;
  let x = 1n << BigInt(Math.ceil(n.toString(2).length / 2));
  for (;;) {
    const y = (x + n / x) >> 1n;
    if (y >= x) return x;
    x = y;
  }
}
/** Native Hadamard bound used by modular determinant and inverse reconstruction.
 * @see Deviation: NTL integer CRT and modular determinant adapters
 */
export function DetBound(
  A: readonly (readonly bigint[])[],
  options: IntegerMatrixOptions = {}
): number {
  const M = integerMatrix(A, options.columns ?? A[0]?.length ?? 0);
  let res = 1n;
  for (const row of M) {
    let t = row.reduce((s, x) => s + x * x, 0n);
    if (t > 1n) t = integerSqrt(t) + 1n;
    res *= t;
  }
  return res ? res.toString(2).length : 0;
}
/** Native matrix CRT with one shared inverse; [modified, matrix, new modulus].
 * @see Deviation: NTL integer CRT and modular determinant adapters
 */
export function CRT(
  gg: readonly (readonly bigint[])[],
  a: bigint,
  G: readonly (readonly bigint[])[],
  p: bigint,
  options: IntegerMatrixCRTOptions = {}
): [number, bigint[][], bigint] {
  if (p <= 1n) throw new Error('zz_pContext: p must be > 1');
  if (p >= 1n << 60n) throw new Error('zz_pContext: modulus too big');
  const columns = options.columns ?? gg[0]?.length ?? 0;
  const residueColumns = options.residueColumns ?? G[0]?.length ?? 0;
  const M = integerMatrix(gg, columns),
    R = integerMatrix(G, residueColumns);
  if (M.length !== R.length || columns !== residueColumns)
    throw new Error('CRT: dimension mismatch');
  if (a <= 0n) throw new RangeError('CRT: require a > 0');
  const norm = (x: bigint, m: bigint) => ((x % m) + m) % m;
  let u = norm(a, p),
    v = p,
    s = 1n,
    t = 0n;
  while (v) {
    const q = u / v;
    [u, v] = [v, u - q * v];
    [s, t] = [t, s - q * t];
  }
  if (u !== 1n) throw new Error('InvMod: inverse undefined');
  const aInv = norm(s, p),
    p1 = p / 2n,
    a1 = a / 2n;
  let modified = 0;
  for (let i = 0; i < M.length; i++)
    for (let j = 0; j < columns; j++) {
      let g = M[i]![j]!;
      if (!CRTInRange(g, a)) {
        modified = 1;
        g = norm(g, a);
        if (g > a1) g -= a;
      }
      let h = norm((norm(R[i]![j]!, p) - norm(g, p)) * aInv, p);
      if (h > p1) h -= p;
      if (h) {
        modified = 1;
        if (!(p & 1n) && g > 0n && h === p1) g -= a * h;
        else g += a * h;
      }
      M[i]![j] = g;
    }
  return [modified, M, a * p];
}

import { CRT as scalarCRT, GenPrime, type RandomStream } from './ZZ.js';
import { FFTPrimeContext, UseFFTPrime, GetFFTPrime } from './FFT.js';
import { determinant as wordDeterminant, relaxed_inv as wordInverse } from './mat_lzz_p.js';
import { determinant as modularDeterminant } from './mat_ZZ_p.js';
export interface IntegerMatrixDeterminantOptions extends IntegerMatrixOptions {
  deterministic?: boolean;
}
export interface IntegerMatrixInverseOptions extends IntegerMatrixOptions {
  status?: false;
}
export interface IntegerMatrixStatusInverseOptions extends IntegerMatrixDeterminantOptions {
  status: true;
  previous?: readonly (readonly bigint[])[];
}
function reconstructionBits(n: bigint): number {
  return n === 0n ? 0 : (n < 0n ? -n : n).toString(2).length;
}
/** Native multimodular determinant, including probabilistic stabilization checks.
 * @see Deviation: NTL integer matrix reconstruction contexts
 */
export function determinant(
  A: readonly (readonly bigint[])[],
  context: FFTPrimeContext,
  stream: RandomStream,
  options: IntegerMatrixDeterminantOptions = {}
): bigint {
  const columns = options.columns ?? A[0]?.length ?? 0;
  const M = integerMatrix(A, columns),
    n = M.length;
  if (columns !== n) throw new Error('determinant: nonsquare matrix');
  if (!n) return 1n;
  const bound = 2 + DetBound(M);
  let res = 0n,
    prod = 1n,
    instable = 1,
    gpCount = 0;
  for (let i = 0; ; i++) {
    if (reconstructionBits(prod) > bound) break;
    if (
      !options.deterministic &&
      !instable &&
      bound > 1000 &&
      reconstructionBits(prod) < 0.25 * bound
    ) {
      const plen = 90 + reconstructionBits(BigInt(Math.max(bound, reconstructionBits(res))));
      const P = GenPrime(plen, stream, { err: 90 + 2 * reconstructionBits(BigInt(gpCount++)) });
      const dd = modularDeterminant(M, P);
      let changed: number;
      [changed, res, prod] = scalarCRT(res, prod, dd, P);
      if (changed) instable = 1;
      else break;
    }
    UseFFTPrime(i, context, stream);
    const p = GetFFTPrime(i, context),
      dd = wordDeterminant(M, p);
    [instable, res, prod] = scalarCRT(res, prod, dd, p, { word: true });
  }
  return res;
}
/** Native integral inverse, or status/adjugate output when status is true.
 * @see Deviation: NTL integer matrix reconstruction contexts
 */
export function inv(
  A: readonly (readonly bigint[])[],
  context: FFTPrimeContext,
  stream: RandomStream,
  options: IntegerMatrixStatusInverseOptions
): [bigint, bigint[][]];
export function inv(
  A: readonly (readonly bigint[])[],
  context: FFTPrimeContext,
  stream: RandomStream,
  options?: IntegerMatrixInverseOptions
): bigint[][];
export function inv(
  A: readonly (readonly bigint[])[],
  context: FFTPrimeContext,
  stream: RandomStream,
  options: IntegerMatrixInverseOptions | IntegerMatrixStatusInverseOptions = {}
): bigint[][] | [bigint, bigint[][]] {
  const columns = options.columns ?? A[0]?.length ?? 0;
  const M = integerMatrix(A, columns),
    n = M.length;
  const previous = options.status ? (options.previous ?? []) : [];
  const retained = integerMatrix(previous, previous[0]?.length ?? 0);
  if (columns !== n) throw new Error('solve: nonsquare matrix');
  if (!n) return options.status ? [1n, []] : [];
  let x = Array.from({ length: n }, () => Array<bigint>(n).fill(0n));
  let d = 0n,
    d1 = 0n,
    dProd = 1n,
    xProd = 1n,
    dInstable = 1,
    xInstable = 1,
    gpCount = 0,
    check = false;
  const bound = 2 + DetBound(M),
    deterministic = options.status && options.deterministic;
  for (let i = 0; ; i++) {
    if ((check || d === 0n) && !dInstable) {
      if (reconstructionBits(dProd) > bound) break;
      if (!deterministic && bound > 1000 && reconstructionBits(dProd) < 0.25 * bound) {
        const plen = 90 + reconstructionBits(BigInt(Math.max(bound, reconstructionBits(d))));
        const P = GenPrime(plen, stream, { err: 90 + 2 * reconstructionBits(BigInt(gpCount++)) });
        const dd = modularDeterminant(M, P);
        let changed: number;
        [changed, d, dProd] = scalarCRT(d, dProd, dd, P);
        if (changed) dInstable = 1;
        else break;
      }
    }
    UseFFTPrime(i, context, stream);
    const p = GetFFTPrime(i, context);
    if (!check) {
      const [dd, xx] = wordInverse(M, p, { relax: false });
      [dInstable, d, dProd] = scalarCRT(d, dProd, dd, p, { word: true });
      if (dd) {
        const adjugate = xx.map((row) => row.map((v) => (v * dd) % p));
        [xInstable, x, xProd] = CRT(x, xProd, adjugate, p);
      } else xInstable = 1;
      if (!dInstable && !xInstable) {
        // Original exact matrix-product certificate, after both reconstructions stabilize.
        const y = x.map((row) =>
          Array.from({ length: n }, (_, j) => row.reduce((s, v, k) => s + v * M[k]![j]!, 0n))
        );
        if (y.every((row, r) => row.every((v, c) => v === (r === c ? d : 0n)))) {
          d1 = d;
          check = true;
        }
      }
    } else {
      const dd = wordDeterminant(M, p);
      [dInstable, d, dProd] = scalarCRT(d, dProd, dd, p, { word: true });
    }
  }
  if (check && d1 !== d) {
    x = x.map((row) =>
      row.map((v) => {
        const numerator = v * d;
        if (numerator % d1 !== 0n) throw new Error('inexact division');
        return numerator / d1;
      })
    );
  }
  if (options.status) return [d, check ? x : retained];
  if (d === -1n) return x.map((row) => row.map((v) => -v));
  if (d !== 1n) throw new Error('inv: non-invertible matrix');
  return x;
}
