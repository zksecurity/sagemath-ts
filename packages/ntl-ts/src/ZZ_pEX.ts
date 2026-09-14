/** @see Deviation: Polynomial Integer Powers and Portable Native Products */
/** Dense-array NTL extension-polynomial power; coefficients are reduced modulo a monic f. */
import { _ZZ_pX_power_kernels } from './ZZ_pX.js';
import { XGCD as coefficientXGCD } from './ZZ_pX1.js';
export function power(
  a: readonly (readonly bigint[])[],
  e: bigint,
  f: readonly bigint[],
  p: bigint
): bigint[][] {
  if (p <= 1n || f.length < 2 || f[f.length - 1] !== 1n)
    throw new RangeError('extension modulus must be monic with positive degree');
  if (e < 0n) throw new Error('power: negative exponent');
  if (e >= 1n << 63n) throw new RangeError('exponent must fit a signed word');
  const { coefficient, normalize, multiply } = extensionKernels(f, p);
  const A = normalize(a.map(coefficient));
  if (e === 0n) return [[1n]];
  if (!A.length || (A.length === 1 && A[0]!.length === 1 && A[0]![0] === 1n)) return A;
  if (BigInt(A.length - 1) > ((1n << 63n) - 2n) / e) throw new Error('overflow in power');
  let result: bigint[][] = [[1n]];
  for (const bit of e.toString(2)) {
    result = multiply(result, result);
    if (bit === '1') result = multiply(result, A);
  }
  return result;
}

/** @see Deviation: Polynomial Roots and Truncated Series */
export function InvTrunc(
  a: readonly (readonly bigint[])[],
  e: number,
  f: readonly bigint[],
  p: bigint
): bigint[][] {
  if (!Number.isSafeInteger(e) || e < 0) throw new Error('InvTrunc: bad args');
  if (e === 0) return [];
  if (p <= 1n || f.length < 2 || f[f.length - 1] !== 1n)
    throw new RangeError('extension modulus must be monic with positive degree');
  const { coefficient, normalize, multiply } = extensionKernels(f, p),
    A = normalize(a.map(coefficient));
  const [d, x] = coefficientXGCD(A[0] ?? [], f, p);
  if (d.length !== 1 || d[0] !== 1n)
    throw new Error("ZZ_pX InvMod: can't compute multiplicative inverse");
  let g = [x];
  const sizes = [e];
  while (sizes[sizes.length - 1]! > 1) sizes.push(Math.ceil(sizes[sizes.length - 1]! / 2));
  // ZZ_pEX.cpp NewtonInv: compute only the new high half at each lift.
  for (let i = sizes.length - 1; i > 0; i--) {
    const k = sizes[i]!,
      l = sizes[i - 1]! - k;
    const g1 = multiply(A.slice(0, k + l), g).slice(k, k + l);
    const g2 = multiply(g1, g).slice(0, l);
    while (g.length < k) g.push([]);
    for (let j = 0; j < g2.length; j++) g[k + j] = coefficient(g2[j]!.map((c) => -c));
    g = normalize(g);
  }
  return g;
}

function extensionKernels(f: readonly bigint[], p: bigint) {
  const canonical = (c: bigint) => ((c % p) + p) % p;
  const coefficient = (v: readonly bigint[]): bigint[] => {
    const r = v.map(canonical);
    while (r.length && r[r.length - 1] === 0n) r.pop();
    while (r.length >= f.length) {
      const c = r[r.length - 1]!,
        shift = r.length - f.length;
      for (let i = 0; i < f.length - 1; i++) r[shift + i] = canonical(r[shift + i]! - c * f[i]!);
      r.pop();
      while (r.length && r[r.length - 1] === 0n) r.pop();
    }
    return r;
  };
  const normalize = (v: bigint[][]): bigint[][] => {
    while (v.length && !v[v.length - 1]!.length) v.pop();
    return v;
  };
  const multiply = (x: bigint[][], y: bigint[][]): bigint[][] => {
    if (!x.length || !y.length) return [];
    if (x.length === 1)
      return normalize(y.map((c) => coefficient(_ZZ_pX_power_kernels.multiply(c, x[0]!, p))));
    if (y.length === 1)
      return normalize(x.map((c) => coefficient(_ZZ_pX_power_kernels.multiply(c, y[0]!, p))));
    // ZZ_pEX.cpp mul/sqr: embed at stride 2*extension_degree-1.
    const stride = 2 * (f.length - 1) - 1;
    const pack = (v: bigint[][]): bigint[] => {
      const out = Array<bigint>(v.length * stride).fill(0n);
      for (let i = 0; i < v.length; i++)
        for (let j = 0; j < v[i]!.length; j++) out[i * stride + j] = v[i]![j]!;
      while (out.length && out[out.length - 1] === 0n) out.pop();
      return out;
    };
    const X = pack(x),
      Y = x === y ? X : pack(y),
      product = _ZZ_pX_power_kernels.multiply(X, Y, p),
      out: bigint[][] = [];
    for (let i = 0; i < product.length; i += stride)
      out.push(coefficient(product.slice(i, i + stride)));
    return normalize(out);
  };
  return { coefficient, normalize, multiply };
}

export function mul(
  a: readonly (readonly bigint[])[],
  b: readonly (readonly bigint[])[],
  f: readonly bigint[],
  p: bigint
): bigint[][] {
  if (p <= 1n || f.length < 2 || f[f.length - 1] !== 1n)
    throw new RangeError('extension modulus must be monic with positive degree');
  const { coefficient, normalize, multiply } = extensionKernels(f, p);
  const A = normalize(a.map(coefficient)),
    B = a === b ? A : normalize(b.map(coefficient));
  return multiply(A, B);
}

/** NTL ZZ_pEX.cpp: modular products and reciprocal-preconditioned remainders.
 * @see Deviation: Polynomial Modular Powers
 */
function extensionPolynomialKernels(f: readonly bigint[], p: bigint) {
  if (p <= 1n || f.length < 2 || f[f.length - 1] !== 1n)
    throw new RangeError('extension modulus must be monic with positive degree');
  const k = extensionKernels(f, p),
    norm = (a: readonly (readonly bigint[])[]) => k.normalize(a.map(k.coefficient));
  const cmul = (a: readonly bigint[], b: readonly bigint[]) =>
    k.coefficient(_ZZ_pX_power_kernels.multiply(a, b, p));
  const csub = (a: readonly bigint[], b: readonly bigint[]) =>
    k.coefficient(
      Array.from({ length: Math.max(a.length, b.length) }, (_, i) => (a[i] ?? 0n) - (b[i] ?? 0n))
    );
  const add = (a: bigint[][], b: bigint[][]) =>
    norm(
      Array.from({ length: Math.max(a.length, b.length) }, (_, i) =>
        Array.from(
          { length: Math.max(a[i]?.length ?? 0, b[i]?.length ?? 0) },
          (_, j) => (a[i]?.[j] ?? 0n) + (b[i]?.[j] ?? 0n)
        )
      )
    );
  const sub = (a: bigint[][], b: bigint[][]) =>
    norm(
      Array.from({ length: Math.max(a.length, b.length) }, (_, i) => csub(a[i] ?? [], b[i] ?? []))
    );
  const inverse = (a: bigint[]) => {
    const [d, x] = coefficientXGCD(a, f, p);
    if (d.length !== 1 || d[0] !== 1n)
      throw new Error("ZZ_pX InvMod: can't compute multiplicative inverse");
    return x;
  };
  const mul = k.multiply;
  const reciprocalRemainder = (
    a: bigint[][],
    b: bigint[][],
    g: bigint[][]
  ): [bigint[][], bigint[][]] => {
    if (a.length < b.length) return [[], a.slice()];
    const n = a.length - b.length + 1,
      q = mul(a.slice().reverse().slice(0, n), g).slice(0, n);
    while (q.length < n) q.push([]);
    q.reverse();
    return [norm(q), sub(a, mul(q, b))];
  };
  const divrem = (a: bigint[][], b: bigint[][]): [bigint[][], bigint[][]] => {
    if (!b.length) throw new Error('ZZ_pEX: division by zero');
    if (a.length < b.length) return [[], a.slice()];
    const n = a.length - b.length + 1;
    if (Math.min(n, b.length) < 16) {
      const r = a.slice(),
        q: bigint[][] = Array.from({ length: n }, () => []),
        inv = inverse(b[b.length - 1]!);
      for (let i = n - 1; i >= 0; i--) {
        q[i] = cmul(r[i + b.length - 1]!, inv);
        for (let j = 0; j < b.length; j++) r[i + j] = csub(r[i + j]!, cmul(q[i]!, b[j]!));
      }
      return [norm(q), norm(r)];
    }
    return reciprocalRemainder(a, b, InvTrunc(b.slice().reverse(), n, f, p));
  };
  const reducer = (b: bigint[][]) => {
    if (b.length <= 1) throw new Error('build(ZZ_pEXModulus,ZZ_pEX): deg(f) <= 0');
    // NTL ModCross is 8; above it, build a reciprocal once for all products.
    const g = b.length - 1 < 8 ? null : InvTrunc(b.slice().reverse(), b.length - 2, f, p);
    return (a: bigint[][]) => {
      if (a.length < b.length) return a.slice();
      if (g && a.length <= 2 * b.length - 3) return reciprocalRemainder(a, b, g)[1];
      return divrem(a, b)[1];
    };
  };
  return { norm, add, sub, mul, inverse, cmul, divrem, reducer };
}

/** NTL sliding-window modular power, including inversion for negative exponents.
 * @see Deviation: Polynomial Modular Powers
 */
export function PowerMod(
  a: readonly (readonly bigint[])[],
  e: bigint,
  m: readonly (readonly bigint[])[],
  f: readonly bigint[],
  p: bigint
): bigint[][] {
  const K = extensionPolynomialKernels(f, p),
    A = K.norm(a),
    M = K.norm(m),
    reduce = K.reducer(M);
  if (A.length >= M.length) throw new Error('PowerMod: bad args');
  const inverse = (a: bigint[][]) => {
    const [d, x] = XGCD(a, M, f, p);
    if (d.length !== 1 || d[0]!.length !== 1 || d[0]![0] !== 1n)
      throw new Error("ZZ_pEX InvMod: can't compute multiplicative inverse");
    return x;
  };
  if (e === 0n) return [[1n]];
  if (e === 1n) return A;
  if (e === -1n) return inverse(A);
  const multiply = (a: bigint[][], b: bigint[][]) => reduce(K.mul(a, b));
  if (e === 2n) return multiply(A, A);
  if (e === -2n) return inverse(multiply(A, A));
  const bits = (e < 0n ? -e : e).toString(2);
  let out: bigint[][] = [[1n]];
  if (bits.length < 16) {
    for (const bit of bits) {
      out = multiply(out, out);
      if (bit === '1') out = multiply(out, A);
    }
  } else {
    let k = 1;
    while (k < 3 && BigInt(bits.length) > BigInt((k + 1) * (k + 2)) * (1n << BigInt(k - 1))) k++;
    const v = [A];
    if (k > 1) {
      const square = multiply(A, A);
      for (let i = 1; i < 1 << (k - 1); i++) v.push(multiply(v[i - 1]!, square));
    }
    let value = 0;
    for (let i = 0; i < bits.length; i++) {
      value = (value << 1) + (bits[i] === '1' ? 1 : 0);
      if (value === 0) out = multiply(out, out);
      else if (value >= 1 << (k - 1) || i === bits.length - 1) {
        let zeros = 0;
        while (!(value & 1)) {
          value >>= 1;
          zeros++;
        }
        for (let n = value; n > 0; n >>= 1) out = multiply(out, out);
        out = multiply(out, v[value >> 1]!);
        while (zeros-- > 0) out = multiply(out, out);
        value = 0;
      }
    }
  }
  return e < 0n ? inverse(out) : out;
}
/** NTL PowerXMod uses modular coefficient shifts instead of multiplying by x.
 * @see Deviation: Polynomial Modular Powers
 */
export function PowerXMod(
  e: bigint,
  m: readonly (readonly bigint[])[],
  f: readonly bigint[],
  p: bigint
): bigint[][] {
  const K = extensionPolynomialKernels(f, p),
    M = K.norm(m),
    reduce = K.reducer(M);
  let out: bigint[][] = [[1n]];
  for (const bit of (e < 0n ? -e : e).toString(2)) {
    out = reduce(K.mul(out, out));
    if (bit === '1') out = reduce([[], ...out]);
  }
  if (e < 0n) {
    const [d, x] = XGCD(out, M, f, p);
    if (d.length !== 1 || d[0]!.length !== 1 || d[0]![0] !== 1n)
      throw new Error("ZZ_pEX InvMod: can't compute multiplicative inverse");
    return x;
  }
  return out;
}

/** NTL ZZ_pEX.cpp XGCD/XHalfGCD, with ordered inverse-matrix recursion.
 * @see Deviation: Polynomial Modular Powers
 */
export function XGCD(
  a: readonly (readonly bigint[])[],
  b: readonly (readonly bigint[])[],
  f: readonly bigint[],
  p: bigint
): [bigint[][], bigint[][], bigint[][]] {
  const { norm, add, sub, mul, inverse, cmul, divrem } = extensionPolynomialKernels(f, p);
  type Matrix = [bigint[][], bigint[][], bigint[][], bigint[][]];
  const identity = (): Matrix => [[[1n]], [], [], [[1n]]];
  const apply = (M: Matrix, u: bigint[][], v: bigint[][]): [bigint[][], bigint[][]] => [
    add(mul(M[0], u), mul(M[1], v)),
    add(mul(M[2], u), mul(M[3], v)),
  ];
  const matrixMul = (A: Matrix, B: Matrix): Matrix => [
    add(mul(A[0], B[0]), mul(A[1], B[2])),
    add(mul(A[0], B[1]), mul(A[1], B[3])),
    add(mul(A[2], B[0]), mul(A[3], B[2])),
    add(mul(A[2], B[1]), mul(A[3], B[3])),
  ];
  const step = (M: Matrix, q: bigint[][]): Matrix => [
    M[2],
    M[3],
    sub(M[0], mul(q, M[2])),
    sub(M[1], mul(q, M[3])),
  ];
  const iterative = (
    u: bigint[][],
    v: bigint[][],
    reduction: number
  ): [Matrix, bigint[][], bigint[][]] => {
    let M = identity();
    const goal = u.length - 1 - reduction;
    while (v.length - 1 > goal) {
      const [q, r] = divrem(u, v);
      [u, v] = [v, r];
      M = step(M, q);
    }
    return [M, u, v];
  };
  const half = (u: bigint[][], v: bigint[][], reduction: number): Matrix => {
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
  const xhalf = (
    u: bigint[][],
    v: bigint[][],
    reduction: number
  ): [Matrix, bigint[][], bigint[][]] => {
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
    Q: bigint[][] = [],
    flag = 0;
  if (!U.length && !V.length) return [[], [[1n]], []];
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
  return [D.map((c) => cmul(c, w)), S.map((c) => cmul(c, w)), T.map((c) => cmul(c, w))];
}

/** NTL ZZ_pEX.cpp eval: Horner over the coefficient extension field.
 * @see Deviation: Polynomial Evaluation and Composition
 */
function evaluate(
  a: readonly (readonly bigint[])[],
  x: readonly bigint[],
  f: readonly bigint[],
  p: bigint
): bigint[] {
  if (p <= 1n || f.length < 2 || f[f.length - 1] !== 1n)
    throw new RangeError('extension modulus must be monic with positive degree');
  const k = extensionKernels(f, p),
    point = k.coefficient(x),
    A = k.normalize(a.map(k.coefficient));
  let value: bigint[] = [];
  for (let i = A.length - 1; i >= 0; i--) {
    const product = _ZZ_pX_power_kernels.multiply(value, point, p),
      c = A[i]!;
    value = k.coefficient(
      Array.from(
        { length: Math.max(product.length, c.length) },
        (_, j) => (product[j] ?? 0n) + (c[j] ?? 0n)
      )
    );
  }
  return value;
}
export { evaluate as eval };
