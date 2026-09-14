/*
 * The Strassen schedule below follows NTL mat_lzz_p.cpp, derived from FLINT.
 * Copyright (C) 2008, Martin Albrecht
 * Copyright (C) 2008, 2009 William Hart.
 * Copyright (C) 2010, Fredrik Johansson
 * This file is part of FLINT.
 * FLINT is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (LGPL) as published
 * by the Free Software Foundation; either version 2.1 of the License, or
 * (at your option) any later version. See reference/ntl/src/mat_lzz_p.cpp.
 */
/** NTL matrix product with native Strassen–Winograd scheduling.
 * @see Deviation: NTL word matrix multiplication adapter
 */
export function mul(
  A: readonly (readonly bigint[])[],
  B: readonly (readonly bigint[])[],
  p: bigint,
  columns = B[0]?.length ?? 0,
  inner = A[0]?.length ?? B.length
): bigint[][] {
  if (p <= 1n) throw new Error('zz_pContext: p must be > 1');
  if (p >= 1n << 60n) throw new Error('zz_pContext: modulus too big');
  if (columns < 0 || inner < 0) throw new Error('SetDims: bad args');
  if (
    inner !== B.length ||
    A.some((row) => row.length !== inner) ||
    B.some((row) => row.length !== columns)
  )
    throw new Error('matrix mul: dimension mismatch');
  const norm = (x: bigint) => ((x % p) + p) % p;
  const aa = A.map((row) => row.map(norm)),
    bb = B.map((row) => row.map(norm));
  const zeros = (r: number, c: number) =>
    Array.from({ length: r }, () => Array<bigint>(c).fill(0n));
  const block = (M: bigint[][], r: number, c: number, nr: number, nc: number) =>
    M.slice(r, r + nr).map((row) => row.slice(c, c + nc));
  const combine = (X: bigint[][], Y: bigint[][], sign: bigint) =>
    X.map((row, i) =>
      row.map((x, j) => {
        const v = x + sign * Y[i]![j]!;
        return v < 0n ? v + p : v >= p ? v - p : v;
      })
    );
  const product = (X: bigint[][], Y: bigint[][], a: number, b: number, c: number): bigint[][] => {
    const C = zeros(a, c);
    if (!a || !b || !c) return C;
    if (a <= 448 || b <= 448 || c <= 448) {
      // Encode each row/column as a nonnegative polynomial. The middle product
      // coefficient is the exact dot product; the slot bound prevents carries.
      const width = BigInt(((p - 1n) * (p - 1n) * BigInt(b)).toString(2).length),
        mask = (1n << width) - 1n,
        shift = BigInt(b - 1) * width;
      const rows = X.map((row) => {
        let v = 0n;
        for (let k = b - 1; k >= 0; k--) v = (v << width) | row[k]!;
        return v;
      });
      const cols = Array.from({ length: c }, (_, j) => {
        let v = 0n;
        for (let k = 0; k < b; k++) v = (v << width) | Y[k]![j]!;
        return v;
      });
      for (let i = 0; i < a; i++)
        for (let j = 0; j < c; j++) C[i]![j] = (((rows[i]! * cols[j]!) >> shift) & mask) % p;
      return C;
    }
    const ar = Math.floor(a / 2),
      br = Math.floor(b / 2),
      cc = Math.floor(c / 2);
    const A11 = block(X, 0, 0, ar, br),
      A12 = block(X, 0, br, ar, br),
      A21 = block(X, ar, 0, ar, br),
      A22 = block(X, ar, br, ar, br);
    const B11 = block(Y, 0, 0, br, cc),
      B12 = block(Y, 0, cc, br, cc),
      B21 = block(Y, br, 0, br, cc),
      B22 = block(Y, br, cc, br, cc);
    let X1 = combine(A11, A21, -1n),
      X2 = combine(B22, B12, -1n);
    let C21 = product(X1, X2, ar, br, cc);
    X1 = combine(A21, A22, 1n);
    X2 = combine(B12, B11, -1n);
    let C22 = product(X1, X2, ar, br, cc);
    X1 = combine(X1, A11, -1n);
    X2 = combine(B22, X2, -1n);
    let C12 = product(X1, X2, ar, br, cc);
    X1 = combine(A12, X1, -1n);
    let C11 = product(X1, B22, ar, br, cc);
    const first = product(A11, B11, ar, br, cc);
    C12 = combine(first, C12, 1n);
    C21 = combine(C12, C21, 1n);
    C12 = combine(C12, C22, 1n);
    C22 = combine(C21, C22, 1n);
    C12 = combine(C12, C11, 1n);
    X2 = combine(X2, B21, -1n);
    C11 = product(A22, X2, ar, br, cc);
    C21 = combine(C21, C11, -1n);
    C11 = product(A12, B21, ar, br, cc);
    C11 = combine(first, C11, 1n);
    for (let i = 0; i < ar; i++) {
      for (let j = 0; j < cc; j++) {
        C[i]![j] = C11[i]![j]!;
        C[i]![j + cc] = C12[i]![j]!;
        C[i + ar]![j] = C21[i]![j]!;
        C[i + ar]![j + cc] = C22[i]![j]!;
      }
    }
    if (c > 2 * cc) {
      const tail = product(X, block(Y, 0, 2 * cc, b, c - 2 * cc), a, b, c - 2 * cc);
      for (let i = 0; i < a; i++) for (let j = 2 * cc; j < c; j++) C[i]![j] = tail[i]![j - 2 * cc]!;
    }
    if (a > 2 * ar) {
      const tail = product(block(X, 2 * ar, 0, a - 2 * ar, b), Y, a - 2 * ar, b, c);
      for (let i = 2 * ar; i < a; i++) C[i] = tail[i - 2 * ar]!;
    }
    if (b > 2 * br) {
      const tail = product(
        block(X, 0, 2 * br, 2 * ar, b - 2 * br),
        block(Y, 2 * br, 0, b - 2 * br, 2 * cc),
        2 * ar,
        b - 2 * br,
        2 * cc
      );
      for (let i = 0; i < 2 * ar; i++)
        for (let j = 0; j < 2 * cc; j++) C[i]![j] = (C[i]![j]! + tail[i]![j]!) % p;
    }
    return C;
  };
  return product(aa, bb, A.length, inner, columns);
}

// Source: mat_lzz_p.cpp basic_inv, basic_tri, elim_basic and their public dispatchers.
// The exact scalar Gaussian schedule has the same cubic arithmetic complexity as
// native packed/blocked kernels; modular reductions replace bounded word accumulators.
export interface WordMatrixOptions {
  columns?: number;
}
export interface WordMatrixEliminationOptions extends WordMatrixOptions {
  w?: number;
}
export interface WordMatrixInverseOptions extends WordMatrixOptions {
  relax?: boolean;
  previous?: readonly (readonly bigint[])[];
}
export interface WordMatrixSolveOptions extends WordMatrixOptions {
  relax?: boolean;
  left?: boolean;
  previous?: readonly bigint[];
}
function linearContext(p: bigint): void {
  if (p <= 1n) throw new Error('zz_pContext: p must be > 1');
  if (p >= 1n << 60n) throw new Error('zz_pContext: modulus too big');
}
const linearNorm = (x: bigint, p: bigint): bigint => ((x % p) + p) % p;
function linearMatrix(A: readonly (readonly bigint[])[], p: bigint, columns: number): bigint[][] {
  if (columns < 0) throw new Error('SetDims: bad args');
  if (A.some((row) => row.length !== columns)) throw new Error('nonrectangular matrix');
  return A.map((row) => row.map((x) => linearNorm(x, p)));
}
function linearInverse(a: bigint, p: bigint, relax: boolean): bigint | null {
  let u = a,
    v = p,
    s = 1n,
    t = 0n;
  while (v) {
    const q = u / v;
    [u, v] = [v, u - q * v];
    [s, t] = [t, s - q * t];
  }
  if (u !== 1n) {
    if (relax) return null;
    throw new Error('InvMod: inverse undefined');
  }
  return linearNorm(s, p);
}
function linearEliminate(
  A: readonly (readonly bigint[])[],
  p: bigint,
  columns: number,
  w: number,
  full: boolean,
  wantImage: boolean,
  wantKernel: boolean
): [number, bigint[][], bigint[][]] {
  linearContext(p);
  const M = linearMatrix(A, p, columns),
    n = A.length;
  if (w < 0 || w > columns) throw new Error('elim: bad args');
  if (!n) return [0, [], []];
  const P = Array.from({ length: n }, (_, i) => i),
    pcol: number[] = [];
  let rank = 0;
  for (let k = 0; k < w; k++) {
    let pos = rank;
    while (pos < n && M[pos]![k] === 0n) pos++;
    if (pos === n) continue;
    const inverse = linearInverse(M[pos]![k]!, p, false)!;
    if (pos !== rank) {
      [M[pos], M[rank]] = [M[rank]!, M[pos]!];
      P[rank] = pos;
    }
    for (let i = rank + 1; i < n; i++) {
      const t = linearNorm(-M[i]![k]! * inverse, p);
      M[i]![k] = t;
      if (!t) continue;
      for (let j = k + 1; j < columns; j++) M[i]![j] = (M[i]![j]! + t * M[rank]![j]!) % p;
    }
    pcol.push(k);
    rank++;
  }
  const im = wantImage
    ? M.slice(0, full ? n : rank).map((row, i) =>
        row.map((x, j) => (j < (i < rank ? pcol[i]! : w) ? 0n : x))
      )
    : [];
  const ker: bigint[][] = [];
  if (wantKernel)
    for (let i = 0; i < n - rank; i++) {
      const row = Array<bigint>(n).fill(0n);
      row[rank + i] = 1n;
      for (let k = rank - 1; k >= 0; k--) {
        let acc = M[rank + i]![pcol[k]!]!;
        for (let j = k + 1; j < rank; j++) acc = (acc + row[j]! * M[j]![pcol[k]!]!) % p;
        row[k] = acc;
      }
      for (let k = n - 1; k >= 0; k--) {
        const pos = P[k]!;
        if (pos !== k) [row[pos], row[k]] = [row[k]!, row[pos]!];
      }
      ker.push(row);
    }
  return [rank, im, ker];
}
/** Native Gaussian elimination with copied full matrix output.
 * @see Deviation: NTL word matrix linear algebra
 */
export function gauss(
  A: readonly (readonly bigint[])[],
  p: bigint,
  options?: WordMatrixEliminationOptions
): [number, bigint[][]] {
  const columns = options?.columns ?? A[0]?.length ?? 0;
  const [rank, M] = linearEliminate(A, p, columns, options?.w ?? columns, true, true, false);
  return [rank, M];
}
/** Native image basis in pivot-row order, without pivot normalization.
 * @see Deviation: NTL word matrix linear algebra
 */
export function image(
  A: readonly (readonly bigint[])[],
  p: bigint,
  options?: WordMatrixOptions
): bigint[][] {
  const columns = options?.columns ?? A[0]?.length ?? 0;
  return linearEliminate(A, p, columns, columns, false, true, false)[1];
}
/** Native left-kernel basis, with the original reverse row-swap order.
 * @see Deviation: NTL word matrix linear algebra
 */
export function kernel(
  A: readonly (readonly bigint[])[],
  p: bigint,
  options?: WordMatrixOptions
): bigint[][] {
  const columns = options?.columns ?? A[0]?.length ?? 0;
  return linearEliminate(A, p, columns, columns, false, false, true)[2];
}
/** Native determinant/inverse status; a singular result retains copied previous output.
 * @see Deviation: NTL word matrix linear algebra
 */
export function relaxed_inv(
  A: readonly (readonly bigint[])[],
  p: bigint,
  options?: WordMatrixInverseOptions
): [bigint, bigint[][]] {
  linearContext(p);
  const columns = options?.columns ?? A[0]?.length ?? 0,
    M = linearMatrix(A, p, columns),
    n = A.length;
  if (columns !== n) throw new Error('inv: nonsquare matrix');
  const P = Array.from({ length: n }, (_, i) => i);
  let det = 1n;
  for (let k = 0; k < n; k++) {
    let pos = -1,
      inverse = 0n;
    for (let i = k; i < n; i++)
      if (M[i]![k]) {
        const value = linearInverse(M[i]![k]!, p, options?.relax ?? true);
        if (value !== null) {
          pos = i;
          inverse = value;
          break;
        }
      }
    if (pos === -1)
      return [0n, (options?.previous ?? []).map((row) => row.map((x) => linearNorm(x, p)))];
    if (pos !== k) {
      [M[pos], M[k]] = [M[k]!, M[pos]!];
      det = linearNorm(-det, p);
      P[k] = pos;
    }
    det = (det * M[k]![k]!) % p;
    for (let j = 0; j < n; j++) M[k]![j] = (M[k]![j]! * inverse) % p;
    M[k]![k] = inverse;
    for (let i = 0; i < n; i++) {
      if (i === k) continue;
      const t = linearNorm(-M[i]![k]!, p);
      M[i]![k] = 0n;
      if (!t) continue;
      for (let j = 0; j < n; j++) M[i]![j] = (M[i]![j]! + t * M[k]![j]!) % p;
    }
  }
  for (const row of M)
    for (let k = n - 1; k >= 0; k--) {
      const pos = P[k]!;
      if (pos !== k) [row[pos], row[k]] = [row[k]!, row[pos]!];
    }
  return [det, M];
}
/** Native strict inverse, raising on nonunits and singular matrices.
 * @see Deviation: NTL word matrix linear algebra
 */
export function inv(
  A: readonly (readonly bigint[])[],
  p: bigint,
  options?: WordMatrixOptions
): bigint[][] {
  const [d, X] = relaxed_inv(A, p, { ...options, relax: false });
  if (!d) throw new Error('inv: non-invertible matrix');
  return X;
}
function linearTri(
  A: readonly (readonly bigint[])[],
  p: bigint,
  b: readonly bigint[] | null,
  options?: WordMatrixSolveOptions
): [bigint, bigint[]] {
  linearContext(p);
  const columns = options?.columns ?? A[0]?.length ?? 0;
  let M = linearMatrix(A, p, columns);
  const n = A.length;
  // The public tri dispatcher uses "inv", before reaching basic_tri's own checks.
  if (columns !== n) throw new Error('inv: nonsquare matrix');
  if (b && b.length !== n) throw new Error('tri: dimension mismatch');
  if (options?.left) M = Array.from({ length: n }, (_, i) => M.map((row) => row[i]!));
  const bv = b?.map((x) => linearNorm(x, p)) ?? [];
  let det = 1n;
  for (let k = 0; k < n; k++) {
    let pos = -1,
      inverse = 0n;
    for (let i = k; i < n; i++)
      if (M[i]![k]) {
        const value = linearInverse(M[i]![k]!, p, options?.relax ?? true);
        if (value !== null) {
          pos = i;
          inverse = value;
          break;
        }
      }
    if (pos === -1) return [0n, (options?.previous ?? []).map((x) => linearNorm(x, p))];
    if (pos !== k) {
      [M[pos], M[k]] = [M[k]!, M[pos]!];
      det = linearNorm(-det, p);
      if (b) [bv[pos], bv[k]] = [bv[k]!, bv[pos]!];
    }
    det = (det * M[k]![k]!) % p;
    for (let j = k + 1; j < n; j++) M[k]![j] = (M[k]![j]! * inverse) % p;
    if (b) bv[k] = (bv[k]! * inverse) % p;
    for (let i = k + 1; i < n; i++) {
      const t = linearNorm(-M[i]![k]!, p);
      if (!t) continue;
      for (let j = k + 1; j < n; j++) M[i]![j] = (M[i]![j]! + t * M[k]![j]!) % p;
      if (b) bv[i] = (bv[i]! + t * bv[k]!) % p;
    }
  }
  const X = Array<bigint>(b ? n : 0).fill(0n);
  if (b)
    for (let i = n - 1; i >= 0; i--) {
      let sum = 0n;
      for (let j = i + 1; j < n; j++) sum = (sum + X[j]! * M[i]![j]!) % p;
      X[i] = linearNorm(bv[i]! - sum, p);
    }
  return [det, X];
}
/** Native determinant with optional nonunit-pivot skipping.
 * @see Deviation: NTL word matrix linear algebra
 */
export function relaxed_determinant(
  A: readonly (readonly bigint[])[],
  p: bigint,
  options?: WordMatrixOptions & { relax?: boolean }
): bigint {
  return linearTri(A, p, null, options)[0];
}
/** Native strict determinant.
 * @see Deviation: NTL word matrix linear algebra
 */
export function determinant(
  A: readonly (readonly bigint[])[],
  p: bigint,
  options?: WordMatrixOptions
): bigint {
  return relaxed_determinant(A, p, { ...options, relax: false });
}
/** Native solve status; left selects x*A=b, otherwise A*x=b.
 * @see Deviation: NTL word matrix linear algebra
 */
export function relaxed_solve(
  A: readonly (readonly bigint[])[],
  b: readonly bigint[],
  p: bigint,
  options?: WordMatrixSolveOptions
): [bigint, bigint[]] {
  return linearTri(A, p, b, options);
}
/** Native strict solve status, preserving previous output on a singular matrix.
 * @see Deviation: NTL word matrix linear algebra
 */
export function solve(
  A: readonly (readonly bigint[])[],
  b: readonly bigint[],
  p: bigint,
  options?: Omit<WordMatrixSolveOptions, 'relax'>
): [bigint, bigint[]] {
  return relaxed_solve(A, b, p, { ...options, relax: false });
}
