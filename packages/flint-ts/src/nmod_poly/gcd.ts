import { _nmod_poly_mul } from './mul.js';
/** Dense-array nmod GCD backend; source: nmod_poly/gcd.c and gr_poly/hgcd.c.
 * Returns the final (not necessarily monic) remainder. Composite moduli are accepted
 * until a required inverse fails; a nonzero constant terminates without inversion.
 * @see Deviation: Dense FLINT polynomial GCD kernels
 */
export function _nmod_poly_gcd(a: readonly bigint[], b: readonly bigint[], p: bigint): bigint[] {
  if (p < 2n) throw new RangeError('modulus must be at least 2');
  let A = normalized(a, p),
    B = normalized(b, p);
  if (A.length < B.length) [A, B] = [B, A];
  if (!B.length) return A;
  const index = Math.min(p.toString(2).length - 1, 63);
  if (B.length < GCD_CUTOFF[index]!) return euclidean(A, B, p);
  let R = divrem(A, B, p)[1];
  if (!R.length) return B;
  if (R.length === 1) return R;
  [, A, B] = halfGcd(B, R, p, INNER_CUTOFF[index]!);
  while (B.length) {
    R = divrem(A, B, p)[1];
    if (!R.length) return B;
    if (R.length === 1) return R;
    if (B.length < OUTER_CUTOFF[index]!) return euclidean(B, R, p);
    [, A, B] = halfGcd(B, R, p, INNER_CUTOFF[index]!);
  }
  return A;
}
// The 64-bit tuning tables from nmod_poly/gcd.c.
const GCD_CUTOFF = [
  470, 657, 626, 689, 689, 723, 796, 876, 835, 876, 919, 964, 964, 1012, 1170, 1228, 919, 1115,
  1115, 1115, 1491, 1491, 1643, 1725, 1811, 1725, 1811, 1725, 1901, 1901, 1901, 1170, 1115, 1170,
  1228, 1170, 1228, 1170, 1228, 1289, 1289, 1228, 1289, 1170, 1353, 1643, 1725, 1811, 1725, 1643,
  1901, 1901, 1901, 1811, 1901, 1901, 1901, 1996, 1901, 1901, 1996, 2199, 1289, 1725,
];
const INNER_CUTOFF = [
  68, 88, 60, 76, 64, 64, 76, 68, 72, 68, 80, 88, 96, 108, 68, 140, 116, 108, 96, 92, 104, 132, 140,
  136, 88, 104, 120, 112, 136, 200, 188, 176, 172, 112, 144, 156, 132, 172, 156, 148, 180, 200, 188,
  164, 160, 168, 168, 168, 168, 244, 180, 180, 200, 204, 192, 228, 212, 208, 168, 196, 188, 216,
  208, 76,
];
const OUTER_CUTOFF = [
  475, 432, 341, 412, 697, 697, 633, 633, 664, 603, 731, 931, 767, 845, 731, 731, 697, 697, 845,
  887, 697, 1129, 931, 1025, 1185, 1025, 1185, 1244, 1371, 1585, 1585, 1510, 375, 498, 432, 375,
  453, 412, 325, 412, 498, 475, 498, 393, 412, 310, 393, 432, 393, 575, 633, 633, 664, 548, 475,
  845, 845, 548, 767, 522, 603, 575, 731, 1510,
];
function normalized(a: readonly bigint[], p: bigint): bigint[] {
  const r = a.map((c) => mod(c, p));
  while (r.length && r[r.length - 1] === 0n) r.pop();
  return r;
}
function mod(a: bigint, p: bigint): bigint {
  const r = a % p;
  return r < 0n ? r + p : r;
}
function inverse(a: bigint, p: bigint): bigint {
  let r = p,
    s = mod(a, p),
    x = 0n,
    y = 1n;
  while (s) {
    const q = r / s;
    [r, s] = [s, r - q * s];
    [x, y] = [y, x - q * y];
  }
  if (r !== 1n) throw new RangeError('coefficient is not invertible');
  return mod(x, p);
}
function add(a: readonly bigint[], b: readonly bigint[], p: bigint): bigint[] {
  return normalized(
    Array.from({ length: Math.max(a.length, b.length) }, (_, i) => (a[i] ?? 0n) + (b[i] ?? 0n)),
    p
  );
}
function sub(a: readonly bigint[], b: readonly bigint[], p: bigint): bigint[] {
  return normalized(
    Array.from({ length: Math.max(a.length, b.length) }, (_, i) => (a[i] ?? 0n) - (b[i] ?? 0n)),
    p
  );
}
/** Share native modular products, including bounded packing for large operands.
 * @see Deviation: Polynomial Modular Powers
 */
function mul(a: readonly bigint[], b: readonly bigint[], p: bigint): bigint[] {
  return _nmod_poly_mul(a, b, p);
}
/** nmod_poly/divrem.c: small classical division, otherwise reversed series inversion. */
function divrem(a: readonly bigint[], b: readonly bigint[], p: bigint): [bigint[], bigint[]] {
  if (!b.length) throw new RangeError('polynomial division by zero');
  if (a.length < b.length) return [[], a.slice()];
  const n = a.length - b.length + 1,
    pb = p.toString(2).length;
  if (
    a.length <= 20 ||
    b.length <= 8 ||
    n <= 7 ||
    (pb <= 61 && a.length <= 40) ||
    (pb <= 29 && a.length <= 70)
  ) {
    const r = a.slice(),
      q = Array<bigint>(n).fill(0n),
      inv = inverse(b[b.length - 1]!, p);
    for (let i = n - 1; i >= 0; i--) {
      q[i] = mod(r[i + b.length - 1]! * inv, p);
      for (let j = 0; j < b.length; j++) r[i + j] = mod(r[i + j]! - q[i]! * b[j]!, p);
    }
    return [normalized(q, p), normalized(r, p)];
  }
  const br = b.slice().reverse();
  let g = [inverse(br[0]!, p)];
  for (let k = 1; k < n; ) {
    const size = Math.min(2 * k, n);
    const correction = sub([2n], mul(br.slice(0, size), g, p).slice(0, size), p);
    g = mul(g, correction, p).slice(0, size);
    k = size;
  }
  const qr = mul(a.slice().reverse().slice(0, n), g, p).slice(0, n);
  while (qr.length < n) qr.push(0n);
  const q = normalized(qr.reverse(), p);
  return [q, sub(a, mul(b, q, p), p)];
}
type Matrix = [bigint[], bigint[], bigint[], bigint[]];
function matrixMul(a: Matrix, b: Matrix, p: bigint): Matrix {
  return [
    add(mul(a[0], b[0], p), mul(a[1], b[2], p), p),
    add(mul(a[0], b[1], p), mul(a[1], b[3], p), p),
    add(mul(a[2], b[0], p), mul(a[3], b[2], p), p),
    add(mul(a[2], b[1], p), mul(a[3], b[3], p), p),
  ];
}
function apply(m: Matrix, a: bigint[], b: bigint[], p: bigint): [bigint[], bigint[]] {
  return [add(mul(m[0], a, p), mul(m[1], b, p), p), add(mul(m[2], a, p), mul(m[3], b, p), p)];
}
/** gr_poly/hgcd.c, with the matrix mapping inputs to remainders (inverse of C's matrix). */
function halfGcd(
  a: bigint[],
  b: bigint[],
  p: bigint,
  cutoff: number,
  res?: ResultantState
): [Matrix, bigint[], bigint[]] {
  const m = Math.floor(a.length / 2),
    identity: Matrix = [[1n], [], [], [1n]];
  if (b.length < m + 1) return [identity, a, b];
  if (a.length < cutoff) {
    let M = identity;
    while (b.length >= m + 1) {
      if (res) res.lc = b[b.length - 1]!;
      const [q, r] = divrem(a, b, p);
      if (res) {
        if (r.length >= m + 1) resultantStep(res, a.length, b.length, r.length, p);
        else {
          res.len0 = a.length;
          res.len1 = b.length;
        }
      }
      M = matrixMul([[], [1n], [1n], q.map((c) => mod(-c, p))], M, p);
      [a, b] = [b, r];
    }
    return [M, a, b];
  }
  if (res) {
    res.lc = b[b.length - 1]!;
    res.len0 -= m;
    res.len1 -= m;
    res.off += m;
  }
  const [N1] = halfGcd(a.slice(m), b.slice(m), p, cutoff, res);
  if (res) {
    res.len0 += m;
    res.len1 += m;
    res.off -= m;
  }
  const [a2, b2] = apply(N1, a, b, p);
  if (b2.length < m + 1) return [N1, a2, b2];
  const k = 2 * m - b2.length + 1;
  if (res) {
    if (b2.length < b.length) resultantStep(res, res.len0, res.len1, b2.length, p);
    res.lc = b2[b2.length - 1]!;
    res.len0 = a2.length;
    res.len1 = b2.length;
  }
  const [q, d] = divrem(a2, b2, p);
  if (res) {
    if (d.length >= m + 1) {
      resultantStep(res, a2.length, b2.length, d.length, p);
      res.len0 = b2.length;
      res.len1 = d.length;
    }
    res.len0 -= k;
    res.len1 -= k;
    res.off += k;
  }
  const [N2] = halfGcd(b2.slice(k), d.slice(k), p, cutoff, res);
  if (res) {
    res.len0 += k;
    res.len1 += k;
    res.off -= k;
  }
  const [A, B] = apply(N2, b2, d, p);
  return [matrixMul(N2, matrixMul([[], [1n], [1n], q.map((c) => mod(-c, p))], N1, p), p), A, B];
}
function euclidean(a: bigint[], b: bigint[], p: bigint): bigint[] {
  // FLINT stops at a nonzero constant without attempting to invert it.
  // This matters over composite moduli, where a constant may be a nonunit.
  while (b.length > 1) [a, b] = [b, divrem(a, b, p)[1]];
  return b.length ? b : a;
}

/** @internal Shared classical/Newton division kernel; public validation is in divrem.ts. */
export { divrem as _nmod_poly_divrem_kernel };

/** @internal Shared arithmetic and inverse-matrix HGCD representation for XGCD. */
export const _nmod_poly_xgcd_kernels = {
  normalized,
  inverse,
  add,
  sub,
  mul,
  divrem,
  halfGcd,
  apply,
  innerCutoff: (n: bigint) => INNER_CUTOFF[Math.min(n.toString(2).length - 1, 63)]!,
};

/** gr_poly/hgcd.c resultant accumulator, with deferred degree-drop factors. */
type ResultantState = { value: bigint; lc: bigint; len0: number; len1: number; off: number };
function power(a: bigint, e: number, n: bigint): bigint {
  let r = 1n,
    exponent = BigInt(e);
  for (; exponent; exponent >>= 1n, a = (a * a) % n) if (exponent & 1n) r = (r * a) % n;
  return r;
}
function resultantStep(
  res: ResultantState,
  len0: number,
  len1: number,
  lenR: number,
  n: bigint
): void {
  if (lenR) {
    res.lc = power(res.lc, len0 - lenR, n);
    res.value = (res.value * res.lc) % n;
    if ((len0 + res.off) % 2 === 0 && (len1 + res.off) % 2 === 0) res.value = mod(-res.value, n);
  } else if (len1 === 1) {
    res.lc = power(res.lc, len0 - 1, n);
    res.value = (res.value * res.lc) % n;
  } else res.value = 0n;
}
function resultantHalfGcd(
  a: bigint[],
  b: bigint[],
  n: bigint,
  value: bigint
): [bigint, bigint[], bigint[]] {
  const res: ResultantState = {
    value,
    lc: b[b.length - 1]!,
    len0: a.length,
    len1: b.length,
    off: 0,
  };
  const [, A, B] = halfGcd(a, b, n, INNER_CUTOFF[Math.min(n.toString(2).length - 1, 63)]!, res);
  if (B.length < b.length) resultantStep(res, res.len0, res.len1, B.length, n);
  return [res.value, A, B];
}
/** @internal Resultant accumulation shares the same recursive HGCD matrix operations. */
export const _nmod_poly_resultant_kernels = {
  normalized,
  divrem,
  mod,
  power,
  resultantStep,
  resultantHalfGcd,
  gcdCutoff: (n: bigint) => GCD_CUTOFF[Math.min(n.toString(2).length - 1, 63)]!,
  outerCutoff: (n: bigint) => OUTER_CUTOFF[Math.min(n.toString(2).length - 1, 63)]!,
};
