/** @see Deviation: Polynomial Roots and Truncated Series */
/** Dense NTL ZZ_pX1.cpp XGCD/half-GCD for extension coefficient inversion. */
import { _ZZ_pX_power_kernels } from './ZZ_pX.js';
type Matrix = [bigint[], bigint[], bigint[], bigint[]];
/** @internal Shared modular polynomial arithmetic used by ZZX reconstruction. */
export function _ZZ_pX_euclidean_kernels(p: bigint) {
  if (p <= 1n) throw new RangeError('modulus must be greater than one');
  const mod = (c: bigint) => ((c % p) + p) % p;
  const norm = (a: readonly bigint[]) => {
    const r = a.map(mod);
    while (r.length && !r[r.length - 1]) r.pop();
    return r;
  };
  const add = (a: bigint[], b: bigint[]) =>
    norm(
      Array.from({ length: Math.max(a.length, b.length) }, (_, i) => (a[i] ?? 0n) + (b[i] ?? 0n))
    );
  const sub = (a: bigint[], b: bigint[]) =>
    norm(
      Array.from({ length: Math.max(a.length, b.length) }, (_, i) => (a[i] ?? 0n) - (b[i] ?? 0n))
    );
  const mul = (a: bigint[], b: bigint[]) => _ZZ_pX_power_kernels.multiply(a, b, p);
  const inverse = (a: bigint) => {
    let r = p,
      s = mod(a),
      u = 0n,
      v = 1n;
    while (s) {
      const q = r / s;
      [r, s] = [s, r - q * s];
      [u, v] = [v, u - q * v];
    }
    if (r !== 1n) throw new Error("ZZ_pX InvMod: can't compute multiplicative inverse");
    return mod(u);
  };
  const divrem = (a: bigint[], b: bigint[]): [bigint[], bigint[]] => {
    if (!b.length) throw new Error('ZZ_pX: division by zero');
    if (a.length < b.length) return [[], a.slice()];
    const length = a.length - b.length + 1;
    if (Math.min(length, b.length) < 32) {
      const r = a.slice(),
        q = Array<bigint>(length).fill(0n),
        inv = inverse(b[b.length - 1]!);
      for (let k = length - 1; k >= 0; k--) {
        const c = mod(r[k + b.length - 1]! * inv);
        q[k] = c;
        for (let j = 0; j < b.length; j++) r[k + j] = mod(r[k + j]! - c * b[j]!);
      }
      return [norm(q), norm(r)];
    }
    // NTL NewtonDivRem: reversed quotient via a truncated reciprocal.
    const reverse = b.slice().reverse();
    let g = [inverse(reverse[0]!)],
      size = 1;
    while (size < length) {
      size = Math.min(2 * size, length);
      g = sub(add(g, g), mul(mul(reverse.slice(0, size), g).slice(0, size), g).slice(0, size));
    }
    const q = mul(a.slice().reverse().slice(0, length), g).slice(0, length);
    while (q.length < length) q.push(0n);
    q.reverse();
    return [norm(q), sub(a, mul(q, b))];
  };
  return { mod, norm, add, sub, mul, inverse, divrem };
}
export function XGCD(
  a: readonly bigint[],
  b: readonly bigint[],
  p: bigint
): [bigint[], bigint[], bigint[]] {
  const { mod, norm, add, sub, mul, inverse, divrem } = _ZZ_pX_euclidean_kernels(p);
  const identity = (): Matrix => [[1n], [], [], [1n]];
  const apply = (M: Matrix, u: bigint[], v: bigint[]): [bigint[], bigint[]] => [
    add(mul(M[0], u), mul(M[1], v)),
    add(mul(M[2], u), mul(M[3], v)),
  ];
  const matrixMul = (A: Matrix, B: Matrix): Matrix => [
    add(mul(A[0], B[0]), mul(A[1], B[2])),
    add(mul(A[0], B[1]), mul(A[1], B[3])),
    add(mul(A[2], B[0]), mul(A[3], B[2])),
    add(mul(A[2], B[1]), mul(A[3], B[3])),
  ];
  const step = (M: Matrix, q: bigint[]): Matrix => [
    M[2],
    M[3],
    sub(M[0], mul(q, M[2])),
    sub(M[1], mul(q, M[3])),
  ];
  const iterative = (u: bigint[], v: bigint[], reduction: number): [Matrix, bigint[], bigint[]] => {
    let M = identity();
    const goal = u.length - 1 - reduction;
    while (v.length - 1 > goal) {
      const [q, r] = divrem(u, v);
      [u, v] = [v, r];
      M = step(M, q);
    }
    return [M, u, v];
  };
  const half = (u: bigint[], v: bigint[], reduction: number): Matrix => {
    if (!v.length || v.length <= u.length - reduction) return identity();
    const shift = Math.max(0, u.length - 1 - 2 * reduction + 2);
    let U = u.slice(shift),
      V = v.slice(shift);
    if (reduction <= 25) return iterative(U, V, reduction)[0];
    const d1 = Math.max(1, Math.min(Math.ceil(reduction / 2), reduction - 1));
    let M = half(U, V, d1);
    [U, V] = apply(M, U, V);
    const d2 = V.length - u.length + shift + reduction;
    if (!V.length || d2 <= 0) return M;
    const [q, r] = divrem(U, V);
    [U, V] = [V, r];
    const N = half(U, V, d2);
    return matrixMul(N, step(M, q));
  };
  const xhalf = (u: bigint[], v: bigint[], reduction: number): [Matrix, bigint[], bigint[]] => {
    if (!v.length || v.length <= u.length - reduction) return [identity(), u, v];
    if (reduction <= 25) return iterative(u, v, reduction);
    const degree = u.length - 1,
      d1 = Math.max(1, Math.min(Math.ceil(reduction / 2), reduction - 1));
    let M = half(u, v, d1);
    [u, v] = apply(M, u, v);
    const d2 = v.length - 1 - degree + reduction;
    if (!v.length || d2 <= 0) return [M, u, v];
    const [q, r] = divrem(u, v);
    [u, v] = [v, r];
    const [N, U, V] = xhalf(u, v, d2);
    return [matrixMul(N, step(M, q)), U, V];
  };
  let U = norm(a),
    V = norm(b),
    Q: bigint[] = [],
    flag = 0;
  if (!U.length && !V.length) return [[], [1n], []];
  if (U.length === V.length) {
    [Q, U] = divrem(U, V);
    [U, V] = [V, U];
    flag = 1;
  } else if (U.length < V.length) {
    [U, V] = [V, U];
    flag = 2;
  }
  const [M, D] = xhalf(U, V, U.length);
  let S = flag === 0 ? M[0] : M[1],
    T = flag === 0 ? M[1] : flag === 1 ? sub(M[0], mul(Q, M[1])) : M[0];
  const w = inverse(D[D.length - 1]!);
  return [D.map((c) => mod(c * w)), S.map((c) => mod(c * w)), T.map((c) => mod(c * w))];
}
