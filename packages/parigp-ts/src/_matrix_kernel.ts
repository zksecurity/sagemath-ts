/** PARI alglin1.c/Flv.c kernel elimination, in zero-indexed row order.
 * Generic fields keep the native unreduced additions and negations; word fields
 * reduce them. The recursive echelon and triangular solves share that schedule.
 */
import { inverseCoefficient, residue } from './_polynomial_division.js';
import { FpM_mul } from './FpV.js';
import { fromColumns, toColumns, type Rows } from './_matrix_mul.js';

type Field = { p: bigint; word: boolean; cutoff: number };
const zero = (m: number, n: number): Rows =>
  Array.from({ length: m }, () => Array<bigint>(n).fill(0n));
const red = (x: bigint, f: Field) => residue(x, f.p);
const add = (x: bigint, y: bigint, f: Field) => (f.word ? red(x + y, f) : x + y);
const neg = (x: bigint, f: Field) => (f.word ? red(-x, f) : -x);
const sub = (a: Rows, b: Rows, f: Field): Rows =>
  a.map((r, i) => r.map((x, j) => add(x, neg(b[i]![j]!, f), f)));
const slice = (a: Rows, r: number, m: number, c: number, n: number): Rows =>
  a.slice(r, r + m).map((row) => row.slice(c, c + n));
const product = (a: Rows, b: Rows, m: number, n: number, k: number, f: Field): Rows => {
  if (!m || !n || !k) return zero(m, k);
  return fromColumns(FpM_mul(toColumns(a, n), toColumns(b, k), f.p));
};
const complement = (r: number[], n: number): number[] => {
  const pivots = new Set(r);
  return Array.from({ length: n }, (_, i) => i).filter((i) => !pivots.has(i));
};

type Echelon = { R: number[]; C: Rows };
/** Native gen_CUP_basecase / Flm_CUP_basecase, requesting only R and C. */
function gen_CUP_basecase(input: Rows, m: number, n: number, f: Field): Echelon {
  const a = input.map((r) => r.slice()),
    R: number[] = [];
  let pr = -1;
  for (let j = 0; j < n; j++) {
    let pc = -1;
    for (pr++; pr < m; pr++) {
      for (let k = j; k < n; k++) {
        a[pr]![k] = red(a[pr]![k]!, f);
        if (pc < 0 && a[pr]![k] !== 0n) pc = k;
      }
      if (pc >= 0) break;
    }
    if (pc < 0) break;
    R.push(pr);
    if (pc !== j) for (const row of a) [row[j], row[pc]] = [row[pc]!, row[j]!];
    const u = inverseCoefficient(a[pr]![j]!, f.p, f.word);
    for (let i = pr + 1; i < m; i++) {
      const v = red(a[i]![j]! * u, f);
      a[i]![j] = v;
      for (let k = j + 1; k < n; k++) a[i]![k] = add(a[i]![k]!, red(a[pr]![k]! * neg(v, f), f), f);
    }
  }
  return { R, C: a.map((r) => r.slice(0, R.length)) };
}

/** Solve L*X=A with implicit unit diagonal, using native block recursion. */
function gen_rsolve_lower_unit(L: Rows, A: Rows, m: number, k: number, f: Field): Rows {
  if (m === 0) return [];
  if (m === 1) return [A[0]!.slice()];
  if (m === 2)
    return [A[0]!.slice(), A[1]!.map((x, j) => add(x, neg(red(A[0]![j]! * L[1]![0]!, f), f), f))];
  const m1 = Math.ceil(m / 2),
    m2 = m - m1;
  const X1 = gen_rsolve_lower_unit(slice(L, 0, m1, 0, m1), slice(A, 0, m1, 0, k), m1, k, f);
  const A2 = sub(slice(A, m1, m2, 0, k), product(slice(L, m1, m2, 0, m1), X1, m2, m1, k, f), f);
  const X2 = gen_rsolve_lower_unit(slice(L, m1, m2, m1, m2), A2, m2, k, f);
  return [...X1, ...X2];
}

/** Native gen_lsolve_lower_unit solves X*L=A, despite its source comment. */
function gen_lsolve_lower_unit(L: Rows, A: Rows, m: number, k: number, f: Field): Rows {
  if (m <= 1) return A.map((r) => r.slice());
  if (m === 2) return A.map((r) => [add(r[0]!, neg(red(r[1]! * L[1]![0]!, f), f), f), r[1]!]);
  const m1 = Math.ceil(m / 2),
    m2 = m - m1;
  const X2 = gen_lsolve_lower_unit(slice(L, m1, m2, m1, m2), slice(A, 0, k, m1, m2), m2, k, f);
  const A1 = sub(slice(A, 0, k, 0, m1), product(X2, slice(L, m1, m2, 0, m1), k, m2, m1, f), f);
  const X1 = gen_lsolve_lower_unit(slice(L, 0, m1, 0, m1), A1, m1, k, f);
  return X1.map((row, i) => row.concat(X2[i]!));
}

/** Column echelon form: same splits, row profile and merge as gen_echelon. */
function gen_echelon(A: Rows, m: number, n: number, f: Field): Echelon {
  if (m < f.cutoff || n < f.cutoff) return gen_CUP_basecase(A, m, n, f);
  const n1 = Math.ceil(n / 2),
    n2 = n - n1;
  const first = gen_echelon(slice(A, 0, m, 0, n1), m, n1, f),
    r1 = first.R.length;
  const A2 = slice(A, 0, m, n1, n2);
  if (!r1) return gen_echelon(A2, m, n2, f);
  if (r1 === m) return first;
  const Rc = complement(first.R, m);
  const C11 = first.R.map((i) => first.C[i]!),
    C21 = Rc.map((i) => first.C[i]!);
  const M12 = gen_rsolve_lower_unit(
    C11,
    first.R.map((i) => A2[i]!),
    r1,
    n2,
    f
  );
  const B2 = sub(
    Rc.map((i) => A2[i]!),
    product(C21, M12, m - r1, r1, n2, f),
    f
  );
  const second = gen_echelon(B2, m - r1, n2, f),
    r2 = second.R.length;
  if (!r2) return first;
  const R2 = second.R.map((i) => Rc[i]!),
    C2 = zero(m, r2);
  for (let i = 0; i < Rc.length; i++) C2[Rc[i]!] = second.C[i]!;
  const R: number[] = [],
    C = zero(m, 0);
  let j1 = 0,
    j2 = 0;
  while (j1 < r1 || j2 < r2) {
    const useFirst = j2 === r2 || (j1 < r1 && first.R[j1]! < R2[j2]!);
    const index = useFirst ? j1++ : j2++,
      source = useFirst ? first.C : C2;
    R.push((useFirst ? first.R : R2)[index]!);
    for (let i = 0; i < m; i++) C[i]!.push(source[i]![index]!);
  }
  return { R, C };
}

/** Native small-matrix Gaussian kernel, preserving the first available pivot. */
function gen_ker(A: Rows, m: number, n: number, f: Field): Rows {
  const x = A.map((r) => r.slice()),
    used = Array<boolean>(m).fill(false),
    d = Array<number>(n).fill(-1);
  for (let k = 0; k < n; k++) {
    let j = 0;
    for (; j < m; j++)
      if (!used[j]) {
        x[j]![k] = red(x[j]![k]!, f);
        if (x[j]![k] !== 0n) break;
      }
    if (j === m) continue;
    const pivot = neg(inverseCoefficient(x[j]![k]!, f.p, f.word), f);
    used[j] = true;
    d[k] = j;
    x[j]![k] = -1n;
    for (let i = k + 1; i < n; i++) x[j]![i] = red(pivot * x[j]![i]!, f);
    for (let t = 0; t < m; t++) {
      if (t === j) continue;
      const c = red(x[t]![k]!, f);
      if (!c) continue;
      x[t]![k] = 0n;
      for (let i = k + 1; i < n; i++) x[t]![i] = red(x[t]![i]! + c * x[j]![i]!, f);
    }
  }
  const out: Rows = [];
  for (let k = 0; k < n; k++)
    if (d[k] === -1) {
      const v = Array<bigint>(n).fill(0n);
      v[k] = 1n;
      for (let i = 0; i < k; i++) if (d[i] !== -1) v[i] = red(x[d[i]!]![k]!, f);
      out.push(v);
    }
  return out;
}

/** Kernel columns; native word and generic echelon thresholds are 8 and 5. */
export function fieldKernel(A: Rows, m: number, n: number, p: bigint, word: boolean): Rows {
  const f = { p, word, cutoff: word ? 8 : 5 };
  if (m < f.cutoff || n < f.cutoff) return gen_ker(A, m, n, f);
  const transpose = Array.from({ length: n }, (_, i) =>
    Array.from({ length: m }, (_, j) => A[j]![i]!)
  );
  const { R, C } = gen_echelon(transpose, n, m, f),
    r = R.length,
    Rc = complement(R, n);
  const S = gen_lsolve_lower_unit(
    R.map((i) => C[i]!),
    Rc.map((i) => C[i]!),
    r,
    n - r,
    f
  );
  return Rc.map((free, i) => {
    const v = Array<bigint>(n).fill(0n);
    v[free] = 1n;
    for (let j = 0; j < r; j++) v[R[j]!] = neg(S[i]![j]!, f);
    return v;
  });
}

/** Native upper-triangular solve, with the generic field's reduction schedule. */
function gen_upper(U: Rows, B: Rows, f: Field, left = false): Rows {
  const n = U.length,
    k = left ? B.length : (B[0]?.length ?? 0);
  if (!n) return B.map((r) => r.slice());
  if (n === 1) {
    const u = inverseCoefficient(U[0]![0]!, f.p, f.word);
    return B.map((r) => r.map((x) => red(x * u, f)));
  }
  if (n === 2) {
    const a = U[0]![0]!,
      b = U[0]![1]!,
      d = U[1]![1]!,
      di = inverseCoefficient(red(a * d, f), f.p, f.word),
      ai = red(d * di, f),
      bi = red(a * di, f);
    if (left)
      return B.map((r) => {
        const x = red(r[0]! * ai, f);
        return [x, red(add(r[1]!, neg(red(x * b, f), f), f) * bi, f)];
      });
    const x2 = B[1]!.map((x) => red(x * bi, f)),
      x1 = B[0]!.map((x, j) => red(add(x, neg(red(b * x2[j]!, f), f), f) * ai, f));
    return [x1, x2];
  }
  const n1 = Math.ceil(n / 2),
    n2 = n - n1,
    U11 = slice(U, 0, n1, 0, n1),
    U12 = slice(U, 0, n1, n1, n2),
    U22 = slice(U, n1, n2, n1, n2);
  if (left) {
    const X1 = gen_upper(U11, slice(B, 0, k, 0, n1), f, true),
      B2 = sub(slice(B, 0, k, n1, n2), product(X1, U12, k, n1, n2, f), f),
      X2 = gen_upper(U22, B2, f, true);
    return X1.map((r, i) => r.concat(X2[i]!));
  }
  const X2 = gen_upper(U22, slice(B, n1, n2, 0, k), f),
    B1 = sub(slice(B, 0, n1, 0, k), product(U12, X2, n1, n2, k, f), f);
  return [...gen_upper(U11, B1, f), ...X2];
}
type CUP = { R: number[]; C: Rows; U: Rows; P: number[] };
/** Native gen_CUP, preserving row profiles, column permutations and five-row cutoff. */
function gen_CUP(A: Rows, m: number, n: number, f: Field): CUP {
  if (m < f.cutoff || n < f.cutoff) {
    const a = A.map((r) => r.slice()),
      R: number[] = [],
      P = Array.from({ length: n }, (_, i) => i);
    let pr = -1;
    for (let j = 0; j < n; j++) {
      let pc = -1;
      for (pr++; pr < m; pr++) {
        for (let k = j; k < n; k++) {
          a[pr]![k] = red(a[pr]![k]!, f);
          if (pc < 0 && a[pr]![k] !== 0n) pc = k;
        }
        if (pc >= 0) break;
      }
      if (pc < 0) break;
      R.push(pr);
      if (pc !== j) {
        for (const r of a) [r[j], r[pc]] = [r[pc]!, r[j]!];
        [P[j], P[pc]] = [P[pc]!, P[j]!];
      }
      const u = inverseCoefficient(a[pr]![j]!, f.p, f.word);
      for (let i = pr + 1; i < m; i++) {
        const v = red(a[i]![j]! * u, f);
        a[i]![j] = v;
        for (let k = j + 1; k < n; k++)
          a[i]![k] = add(a[i]![k]!, red(a[pr]![k]! * neg(v, f), f), f);
      }
    }
    return { R, C: a.map((r) => r.slice(0, R.length)), U: R.map((i) => a[i]!.slice()), P };
  }
  const m1 = Math.ceil(Math.min(m, n) / 2),
    m2 = m - m1,
    first = gen_CUP(A.slice(0, m1), m1, n, f),
    r1 = first.R.length,
    A2 = A.slice(m1);
  if (!r1) {
    const second = gen_CUP(A2, m2, n, f);
    return {
      ...second,
      R: second.R.map((i) => i + m1),
      C: [...zero(m1, second.R.length), ...second.C],
    };
  }
  const U11 = slice(first.U, 0, r1, 0, r1),
    U12 = slice(first.U, 0, r1, r1, n - r1),
    T21 = A2.map((r) => first.P.slice(0, r1).map((j) => r[j]!)),
    T22 = A2.map((r) => first.P.slice(r1).map((j) => r[j]!)),
    C21 = gen_upper(U11, T21, f, true),
    B2 = sub(T22, product(C21, U12, m2, r1, n - r1, f), f),
    second = gen_CUP(B2, m2, n - r1, f),
    r2 = second.R.length;
  return {
    R: first.R.concat(second.R.map((i) => i + m1)),
    C: [
      ...first.C.map((r) => r.concat(Array<bigint>(r2).fill(0n))),
      ...C21.map((r, i) => r.concat(second.C[i]!)),
    ],
    U: [
      ...U11.map((r, i) => r.concat(second.P.map((j) => U12[i]![j]!))),
      ...second.U.map((r) => Array<bigint>(r1).fill(0n).concat(r)),
    ],
    P: first.P.slice(0, r1).concat(second.P.map((j) => first.P[r1 + j]!)),
  };
}

/** Native generic Fp matrix pivot rows (one-based, zero for dependent columns).
 * @see Deviation: PARI prime-decomposition matrix adapters
 */
export function fieldMatrixPivots(
  A: Rows,
  m: number,
  n: number,
  p: bigint
): [number[] | null, number] {
  if (!n) return [null, 0];
  const f: Field = { p, word: false, cutoff: 5 },
    d = Array<number>(n).fill(0);
  if (m >= 5 && n >= 5) {
    const { R, P } = gen_CUP(A, m, n, f);
    R.forEach((i, j) => {
      d[P[j]!] = i + 1;
    });
    return [d, n - R.length];
  }
  const a = A.map((r) => r.slice()),
    used = Array<boolean>(m).fill(false);
  let r = 0;
  for (let k = 0; k < n; k++) {
    let j = 0;
    for (; j < m; j++)
      if (!used[j]) {
        a[j]![k] = red(a[j]![k]!, f);
        if (a[j]![k] !== 0n) break;
      }
    if (j === m) {
      r++;
      continue;
    }
    const piv = -inverseCoefficient(a[j]![k]!, p, false);
    used[j] = true;
    d[k] = j + 1;
    for (let i = k + 1; i < n; i++) a[j]![i] = red(piv * a[j]![i]!, f);
    for (let t = 0; t < m; t++)
      if (!used[t]) {
        const c = red(a[t]![k]!, f);
        if (!c) continue;
        a[t]![k] = 0n;
        for (let i = k + 1; i < n; i++) a[t]![i] = red(a[t]![i]! + c * a[j]![i]!, f);
      }
    for (let i = k; i < n; i++) a[j]![i] = 0n;
  }
  return [d, r];
}

/** Native generic gen_gauss, including the five-column CUP solve.
 * @see Deviation: PARI prime-decomposition matrix adapters
 */
export function fieldMatrixSolve(
  A: Rows,
  B: Rows,
  m: number,
  n: number,
  k: number,
  p: bigint
): Rows | null {
  if (!n) return [];
  const f: Field = { p, word: false, cutoff: 5 };
  if (n >= 5) {
    if (m < n) return null;
    const { R, C, U, P } = gen_CUP(A, m, n, f);
    if (R.length < n) return null;
    const Y = gen_rsolve_lower_unit(
        R.map((i) => C[i]!),
        R.map((i) => B[i]!),
        n,
        k,
        f
      ),
      X = gen_upper(U, Y, f),
      out: Rows = [];
    P.forEach((j, i) => {
      out[j] = X[i]!;
    });
    return out;
  }
  const a = A.map((r) => r.slice()),
    b = B.map((r) => r.slice());
  for (let i = 0; i < n; i++) {
    let j = i;
    for (; j < m; j++) {
      const v = red(a[j]![i]!, f);
      if (v) {
        a[j]![i] = inverseCoefficient(v, p, false);
        break;
      }
      a[j]![i] = 0n;
    }
    if (j === m) return null;
    if (j !== i) {
      for (let col = i; col < n; col++) [a[i]![col], a[j]![col]] = [a[j]![col]!, a[i]![col]!];
      [b[i], b[j]] = [b[j]!, b[i]!];
    }
    if (i === n - 1) break;
    for (let row = i + 1; row < m; row++) {
      let v = red(a[row]![i]!, f);
      a[row]![i] = 0n;
      if (!v) continue;
      v = red(-v * a[i]![i]!, f);
      for (let col = i + 1; col < n; col++) {
        a[i]![col] = red(a[i]![col]!, f);
        a[row]![col]! += v * a[i]![col]!;
      }
      for (let col = 0; col < k; col++) {
        b[i]![col] = red(b[i]![col]!, f);
        b[row]![col]! += v * b[i]![col]!;
      }
    }
  }
  const out = zero(n, k);
  for (let col = 0; col < k; col++)
    for (let i = n - 1; i >= 0; i--) {
      let v = b[i]![col]!;
      for (let j = i + 1; j < n; j++) v -= a[i]![j]! * out[j]![col]!;
      out[i]![col] = red(red(v, f) * a[i]![i]!, f);
    }
  return out;
}
