/** Native PARI extension half-GCD recursion, matrix schedules and unscaled Bézout outputs.
 * @see Deviation: PARI extension-polynomial GCD adapters
 */
import { type ExtensionPolynomial as P, validateExtensionInputs } from './_extension_polynomial.js';
import { extensionField, trimExtension as trim } from './_extension_field.js';
import { extensionDivision } from './_extension_division.js';
import { FpX_red } from './ffinit.js';
export type ExtensionMatrix<A extends P = P> = [[A, A], [A, A]];
type PolynomialMatrix = ExtensionMatrix;
type HalfGcdResult = [PolynomialMatrix, P, P];
interface GcdArithmetic {
  add(a: P, b: P): P;
  sub(a: P, b: P): P;
  mul(a: P, b: P): P;
  divrem(a: P, b: P): [P, P];
  rem(a: P, b: P): P;
  halfLimit: number;
  gcdLimit: number;
  extendedLimit: number;
  z: bigint | bigint[];
  u: bigint | bigint[];
  div(a: P, b: P): P;
}
const identity = (c: GcdArithmetic): PolynomialMatrix => [
  [[c.u], []],
  [[], [c.u]],
];
const swapMatrix = (c: GcdArithmetic): PolynomialMatrix => [
  [[], [c.u]],
  [[c.u], []],
];
const low = (a: P, n: number): P => trim(a.slice(0, n));
const shift = (a: P, n: number, c: GcdArithmetic): P =>
  !a.length ? [] : n < 0 ? a.slice(-n) : new Array(n).fill(c.z).concat(a);
function apply(M: PolynomialMatrix, a: P, b: P, c: GcdArithmetic): [P, P] {
  return [c.add(c.mul(M[0][0], a), c.mul(M[0][1], b)), c.add(c.mul(M[1][0], a), c.mul(M[1][1], b))];
}
/** Native Strassen seven-product 2x2 polynomial matrix multiplication. */
function multiplyMatrix(
  A: PolynomialMatrix,
  B: PolynomialMatrix,
  c: GcdArithmetic
): PolynomialMatrix {
  const [[a, b], [d, e]] = A,
    [[f, g], [h, i]] = B;
  const m1 = c.mul(c.add(a, e), c.add(f, i)),
    m2 = c.mul(c.add(d, e), f),
    m3 = c.mul(a, c.sub(g, i)),
    m4 = c.mul(e, c.sub(h, f)),
    m5 = c.mul(c.add(a, b), i),
    m6 = c.mul(c.sub(d, a), c.add(f, g)),
    m7 = c.mul(c.sub(b, e), c.add(h, i));
  return [
    [c.add(c.add(m1, m4), c.sub(m7, m5)), c.add(m3, m5)],
    [c.add(m2, m4), c.add(c.sub(m1, m2), c.add(m3, m6))],
  ];
}
function quotientMatrix(q: P, M: PolynomialMatrix, c: GcdArithmetic): PolynomialMatrix {
  return [
    [M[1][0], M[1][1]],
    [c.sub(M[0][0], c.mul(M[1][0], q)), c.sub(M[0][1], c.mul(M[1][1], q))],
  ];
}
function polynomialHalfGcdBasecase(a: P, b: P, c: GcdArithmetic): HalfGcdResult {
  const n = Math.floor(a.length / 2);
  let M = identity(c);
  while (b.length > n) {
    const [q, r] = c.divrem(a, b);
    a = b;
    b = r;
    M = quotientMatrix(q, M, c);
  }
  return [M, a, b];
}
function halfInternal(x: P, y: P, c: GcdArithmetic): HalfGcdResult {
  if (x.length < c.halfLimit) return polynomialHalfGcdBasecase(x, y, c);
  const n = Math.floor(x.length / 2);
  if (y.length <= n) return [identity(c), x.slice(), y.slice()];
  const [R, a, b] = halfInternal(shift(x, -n, c), shift(y, -n, c), c);
  const [v1, v2] = apply(R, low(x, n), low(y, n), c);
  const x1 = c.add(shift(a, n, c), v1),
    y1 = c.add(shift(b, n, c), v2);
  if (y1.length <= n) return [R, x1, y1];
  const k = 2 * n - (y1.length - 1),
    [q, r] = c.divrem(x1, y1);
  const [S, a2, b2] = halfInternal(shift(y1, -k, c), shift(r, -k, c), c);
  const T = multiplyMatrix(S, quotientMatrix(q, R, c), c);
  const [w1, w2] = apply(S, low(y1, k), low(r, k), c);
  return [T, c.add(shift(a2, k, c), w1), c.add(shift(b2, k, c), w2)];
}
/** Native initial swap/remainder conventions; half-GCD does not reduce inputs first. */
function polynomialHalfGcd(x: P, y: P, c: GcdArithmetic): HalfGcdResult {
  if (!x.length) return [swapMatrix(c), y.slice(), x.slice()];
  if (y.length < x.length) return halfInternal(x, y, c);
  const [q, r] = c.divrem(y, x),
    [R, a, b] = halfInternal(x, r, c);
  R[0][0] = c.sub(R[0][0], c.mul(q, R[0][1]));
  R[1][0] = c.sub(R[1][0], c.mul(q, R[1][1]));
  return [R, a, b];
}
function polynomialGcd(x: P, y: P, c: GcdArithmetic): P {
  if (!x.length) return y.slice();
  while (y.length >= c.gcdLimit) {
    if (y.length <= Math.floor(x.length / 2)) {
      const r = c.rem(x, y);
      x = y;
      y = r;
    }
    [, x, y] = polynomialHalfGcd(x, y, c);
  }
  while (y.length) {
    const r = c.rem(x, y);
    x = y;
    y = r;
  }
  return x.slice();
}
function extendedBasecase(a: P, b: P, c: GcdArithmetic, needU: boolean): [P, P, P] {
  const A = a,
    B = b;
  let v: P = [],
    v1: P = [c.u];
  while (b.length) {
    const [q, r] = c.divrem(a, b);
    [v, v1] = [v1, c.sub(v, c.mul(q, v1))];
    a = b;
    b = r;
  }
  const u = needU ? c.div(c.sub(a, c.mul(B, v)), A) : [];
  return [a, u, v];
}
/** Collect native step matrices, then replay them backwards onto final cofactors. */
function polynomialExtendedGcd(x: P, y: P, c: GcdArithmetic, needU = true): [P, P, P] {
  if (y.length < c.extendedLimit) return extendedBasecase(x, y, c, needU);
  const matrices: PolynomialMatrix[] = [];
  while (y.length >= c.extendedLimit) {
    if (y.length <= Math.floor(x.length / 2)) {
      const [q, r] = c.divrem(x, y);
      x = y;
      y = r;
      matrices.push([
        [[], [c.u]],
        [[c.u], c.sub([], q)],
      ]);
    } else {
      const [M, a, b] = polynomialHalfGcd(x, y, c);
      matrices.push(M);
      x = a;
      y = b;
    }
  }
  let [d, u, v] = extendedBasecase(x, y, c, true);
  for (let i = matrices.length - 1; i >= 0; i--) {
    const R = matrices[i]!;
    const uu = i === 0 && !needU ? [] : c.add(c.mul(u, R[0][0]), c.mul(v, R[1][0]));
    const vv = c.add(c.mul(u, R[0][1]), c.mul(v, R[1][1]));
    u = uu;
    v = vv;
  }
  return [d, u, v];
}

/** Internal modes: GCD, both Bézout outputs, half-GCD matrix, V-only Bézout. */
export function extensionGcd(
  mode: 0 | 1 | 2,
  op: 0 | 1 | 2 | 3,
  p: bigint,
  T: bigint | bigint[],
  a: P,
  b: P,
  innerInverse?: bigint[]
): P | [P, P, P] | PolynomialMatrix | [P, P] {
  validateExtensionInputs(mode, p, T, a, b);
  a = trim(a);
  b = trim(b);
  if (mode === 0 && p < 1n << 64n) {
    const convert = (v: P) => trim(v.map((c) => FpX_red(typeof c === 'bigint' ? [c] : c, p)));
    const collapse = (v: P) =>
      v.map((c) => ((c as bigint[]).length < 2 ? ((c as bigint[])[0] ?? 0n) : c));
    const r = extensionGcd(1, op, p, FpX_red(T as bigint[], p), convert(a), convert(b), innerInverse === undefined ? undefined : FpX_red(innerInverse,p));
    return op === 0
      ? collapse(r as P)
      : op === 2
        ? ((r as PolynomialMatrix).map((row) => row.map(collapse)) as PolynomialMatrix)
        : ((r as P[]).map(collapse) as [P, P, P] | [P, P]);
  }
  const f = extensionField(mode, T, p, innerInverse);
  const add = (a: P, b: P) =>
    trim(
      Array.from({ length: Math.max(a.length, b.length) }, (_, i) =>
        i < a.length ? (i < b.length ? f.add(a[i]!, b[i]!) : a[i]!) : b[i]!
      )
    );
  const c: GcdArithmetic = {
    z: f.z,
    u: f.u,
    add,
    sub: f.psub,
    mul: f.pmul,
    divrem: (a, b) => extensionDivision(mode, 0, p, T, a, b, innerInverse) as [P, P],
    rem: (a, b) => extensionDivision(mode, 1, p, T, a, b, innerInverse) as P,
    div: (a, b) => extensionDivision(mode, 4, p, T, a, b, innerInverse) as P,
    halfLimit: mode === 0 ? 29 : mode === 1 ? 30 : 89,
    gcdLimit: mode === 0 ? 101 : mode === 1 ? 116 : 310,
    extendedLimit: mode === 1 ? 21 : 10,
  };
  if (op === 2) return polynomialHalfGcd(a, b, c)[0];
  a = f.pred(a);
  b = f.pred(b);
  if (op === 0) return polynomialGcd(a, b, c);
  const r = polynomialExtendedGcd(a, b, c, op === 1);
  return op === 3 ? [r[0], r[2]] : r;
}
