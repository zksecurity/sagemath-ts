/** Native PARI half-GCD recursion and polynomial matrix schedules. */
type P = bigint[];
export type PolynomialMatrix = [[P, P], [P, P]];
export type HalfGcdResult = [PolynomialMatrix, P, P];
export interface GcdArithmetic {
  add(a: P, b: P): P;
  sub(a: P, b: P): P;
  mul(a: P, b: P): P;
  divrem(a: P, b: P): [P, P];
  rem(a: P, b: P): P;
  halfLimit: number;
  gcdLimit: number;
  extendedLimit: number;
  word: boolean;
}
import { Fp_powu } from './arith1.js';
import { Fp_mul, Fp_neg } from './ff.js';
interface ResultantState {
  value: bigint;
  leading: bigint;
  degree0: number;
  degree1: number;
  offset: number;
  p: bigint;
}
/** FpX_halfres_update / Flx_halfres_update_pre, including deferred degree drops. */
function updateResultant(da: number, db: number, dr: number, state: ResultantState): void {
  if (dr >= 0) {
    if (state.leading !== 1n) {
      state.leading = Fp_powu(state.leading, BigInt(da - dr), state.p);
      state.value = Fp_mul(state.value, state.leading, state.p);
    }
    if ((da + state.offset) % 2 !== 0 && (db + state.offset) % 2 !== 0)
      state.value = Fp_neg(state.value, state.p);
  } else if (db === 0) {
    if (state.leading !== 1n) {
      state.leading = Fp_powu(state.leading, BigInt(da), state.p);
      state.value = Fp_mul(state.value, state.leading, state.p);
    }
  } else state.value = 0n;
}
const identity = (): PolynomialMatrix => [
  [[1n], []],
  [[], [1n]],
];
const swapMatrix = (): PolynomialMatrix => [
  [[], [1n]],
  [[1n], []],
];
const low = (a: P, n: number): P => {
  const b = a.slice(0, n);
  while (b.length && b[b.length - 1] === 0n) b.pop();
  return b;
};
const shift = (a: P, n: number): P =>
  !a.length ? [] : n < 0 ? a.slice(-n) : new Array<bigint>(n).fill(0n).concat(a);
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
export function polynomialHalfGcdBasecase(a: P, b: P, c: GcdArithmetic): HalfGcdResult {
  return halfBasecase(a, b, c);
}
function halfBasecase(a: P, b: P, c: GcdArithmetic, state?: ResultantState): HalfGcdResult {
  const n = Math.floor(a.length / 2);
  let M = identity();
  while (b.length > n) {
    const [q, r] = c.divrem(a, b);
    if (state) {
      const da = a.length - 1,
        db = b.length - 1,
        dr = r.length - 1;
      state.leading = b[db]!;
      if (dr >= n) updateResultant(da, db, dr, state);
      else {
        state.degree0 = da;
        state.degree1 = db;
      }
    }
    a = b;
    b = r;
    M = quotientMatrix(q, M, c);
  }
  return [M, a, b];
}
function halfInternal(x: P, y: P, c: GcdArithmetic, state?: ResultantState): HalfGcdResult {
  if (x.length < c.halfLimit)
    return state ? halfBasecase(x, y, c, state) : polynomialHalfGcdBasecase(x, y, c);
  const n = Math.floor(x.length / 2);
  if (y.length <= n) return [identity(), x.slice(), y.slice()];
  if (state) {
    state.leading = y[y.length - 1]!;
    state.degree0 -= n;
    state.degree1 -= n;
    state.offset += n;
  }
  const [R, a, b] = halfInternal(shift(x, -n), shift(y, -n), c, state);
  if (state) {
    state.offset -= n;
    state.degree0 += n;
    state.degree1 += n;
  }
  const [v1, v2] = apply(R, low(x, n), low(y, n), c);
  const x1 = c.add(shift(a, n), v1),
    y1 = c.add(shift(b, n), v2);
  if (y1.length <= n) return [R, x1, y1];
  const k = 2 * n - (y1.length - 1),
    [q, r] = c.divrem(x1, y1);
  if (state) {
    const dx1 = x1.length - 1,
      dy1 = y1.length - 1,
      dr = r.length - 1;
    if (dy1 < y.length - 1) updateResultant(state.degree0, state.degree1, dy1, state);
    state.leading = y1[dy1]!;
    state.degree0 = dx1;
    state.degree1 = dy1;
    if (dr >= n) {
      updateResultant(dx1, dy1, dr, state);
      state.degree0 = dy1;
      state.degree1 = dr;
    }
    state.degree0 -= k;
    state.degree1 -= k;
    state.offset += k;
  }
  const [S, a2, b2] = halfInternal(shift(y1, -k), shift(r, -k), c, state);
  if (state) {
    state.degree0 += k;
    state.degree1 += k;
    state.offset -= k;
  }
  const T = multiplyMatrix(S, quotientMatrix(q, R, c), c);
  const [w1, w2] = apply(S, low(y1, k), low(r, k), c);
  return [T, c.add(shift(a2, k), w1), c.add(shift(b2, k), w2)];
}
/** FpX/Flx_halfgcd_all: native initial swap/remainder conventions. */
export function polynomialHalfGcd(x: P, y: P, c: GcdArithmetic): HalfGcdResult {
  if (!x.length) return [swapMatrix(), y.slice(), x.slice()];
  if (y.length < x.length) return halfInternal(x, y, c);
  const [q, r] = c.divrem(y, x),
    [R, a, b] = halfInternal(x, r, c);
  R[0][0] = c.sub(R[0][0], c.mul(q, R[0][1]));
  R[1][0] = c.sub(R[1][0], c.mul(q, R[1][1]));
  return [R, a, b];
}
export function polynomialGcd(x: P, y: P, c: GcdArithmetic): P {
  if (!x.length) return y.slice();
  while (y.length >= c.gcdLimit) {
    if (y.length <= Math.floor(x.length / 2)) {
      const r = c.rem(x, y);
      x = y;
      y = r;
    }
    [, x, y] = polynomialHalfGcd(x, y, c);
  }
  if (c.word && y.length > x.length) [x, y] = [y, x];
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
  let u: P = [1n],
    u1: P = [],
    v: P = [],
    v1: P = [1n];
  if (!b.length) return [a.slice(), u, v];
  while (b.length) {
    const [q, r] = c.divrem(a, b);
    a = b;
    b = r;
    if (c.word && needU) {
      [u, u1] = [u1, u];
      u1 = c.sub(u1, c.mul(u, q));
    }
    [v, v1] = [v1, v];
    if (!c.word && !b.length) break;
    v1 = c.sub(v1, c.mul(v, q));
  }
  if (!c.word && needU) u = c.divrem(c.sub(a, c.mul(B, v)), A)[0];
  return [a, u, v];
}
/** Collect native step matrices, then replay them backwards onto final cofactors. */
export function polynomialExtendedGcd(x: P, y: P, c: GcdArithmetic, needU = true): [P, P, P] {
  if (y.length < c.extendedLimit) return extendedBasecase(x, y, c, needU);
  const matrices: PolynomialMatrix[] = [];
  while (y.length >= c.extendedLimit) {
    if (y.length <= Math.floor(x.length / 2)) {
      const [q, r] = c.divrem(x, y);
      x = y;
      y = r;
      matrices.push([
        [[], [1n]],
        [[1n], c.sub([], q)],
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

/** Native half-resultant accumulation, sharing the half-GCD schedule. */
export function polynomialResultant(x: P, y: P, c: GcdArithmetic, p: bigint): bigint {
  if (!x.length || !y.length) return 0n;
  let result = 1n;
  if (x.length < y.length) {
    if (x.length % 2 === 0 && y.length % 2 === 0) result = Fp_neg(result, p);
    [x, y] = [y, x];
  }
  while (y.length >= c.gcdLimit) {
    if (y.length <= Math.floor(x.length / 2)) {
      const r = c.rem(x, y),
        dx = x.length - 1,
        dy = y.length - 1,
        dr = r.length - 1;
      const leading = y[dy]!;
      if (leading !== 1n) result = Fp_mul(result, Fp_powu(leading, BigInt(dx - dr), p), p);
      if (dx % 2 && dy % 2) result = Fp_neg(result, p);
      x = y;
      y = r;
    }
    const state: ResultantState = {
      value: result,
      leading: y[y.length - 1] ?? 0n,
      degree0: x.length - 1,
      degree1: y.length - 1,
      offset: 0,
      p,
    };
    const dy = y.length - 1;
    [, x, y] = halfInternal(x, y, c, state);
    if (y.length - 1 < dy) updateResultant(state.degree0, state.degree1, y.length - 1, state);
    result = state.value;
  }
  if (!x.length || !y.length) return 0n;
  // Native Euclidean basecase: the accumulated leading powers and signs are
  // multiplied by the final constant to the degree of the preceding divisor.
  while (y.length > 1) {
    const dx = x.length - 1,
      dy = y.length - 1,
      leading = y[dy]!,
      r = c.rem(x, y);
    if (!r.length) return 0n;
    if (dx % 2 && dy % 2) result = Fp_neg(result, p);
    if (leading !== 1n)
      result = Fp_mul(result, Fp_powu(leading, BigInt(dx - (r.length - 1)), p), p);
    x = y;
    y = r;
  }
  return Fp_mul(result, Fp_powu(y[0]!, BigInt(x.length - 1), p), p);
}
