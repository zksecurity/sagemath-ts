import { _fmpz_poly_resultant_kernels } from './gcd.js';
/** @see Deviation: Polynomial Integer Powers and Portable Native Products */
/** Dense array port of fmpz_poly/divrem.c and its basecase/divide-and-conquer kernels. */
export function _fmpz_poly_divrem(
  a: readonly bigint[],
  b: readonly bigint[],
  exact?: false
): [bigint[], bigint[]];
export function _fmpz_poly_divrem(
  a: readonly bigint[],
  b: readonly bigint[],
  exact: true
): [bigint[], bigint[]] | null;
export function _fmpz_poly_divrem(
  a: readonly bigint[],
  b: readonly bigint[],
  exact: boolean
): [bigint[], bigint[]] | null;
export function _fmpz_poly_divrem(
  a: readonly bigint[],
  b: readonly bigint[],
  exact: boolean = false
): [bigint[], bigint[]] | null {
  const A = normalized(a),
    B = normalized(b);
  if (!B.length) throw new RangeError('division by zero polynomial');
  return divide(A, B, exact);
}
function normalized(a: readonly bigint[]): bigint[] {
  const r = a.slice();
  while (r.length && r[r.length - 1] === 0n) r.pop();
  return r;
}
function abs(a: bigint): bigint {
  return a < 0n ? -a : a;
}
function floorDivide(a: bigint, b: bigint): bigint {
  const q = a / b,
    r = a % b;
  return r !== 0n && r < 0n !== b < 0n ? q - 1n : q;
}
// The same FLINT multiplication boundary is shared with GCD and power kernels.
// Bounded packed products avoid engine temporary limits during divide-and-conquer.
function multiply(a: readonly bigint[], b: readonly bigint[]): bigint[] {
  return _fmpz_poly_resultant_kernels.multiply(a, b);
}
function subtract(a: readonly bigint[], b: readonly bigint[]): bigint[] {
  return normalized(
    Array.from({ length: Math.max(a.length, b.length) }, (_, i) => (a[i] ?? 0n) - (b[i] ?? 0n))
  );
}
function shifted(a: readonly bigint[], n: number): bigint[] {
  return a.length ? Array<bigint>(n).fill(0n).concat(a) : [];
}
/** The 2n-1 by n recursion in divrem_divconquer_recursive.c, cutoff 16. */
function balancedDivision(
  a: readonly bigint[],
  b: readonly bigint[],
  exact: boolean
): [bigint[], bigint[]] | null {
  const n = b.length;
  if (n <= 16) {
    const r = a.slice(),
      q = Array<bigint>(n).fill(0n),
      lead = b[n - 1]!;
    for (let i = n - 1; i >= 0; i--) {
      const c = r[i + n - 1] ?? 0n;
      if (exact && c % lead !== 0n) return null;
      // FLINT skips coefficients smaller in magnitude than the divisor lead,
      // then uses floor division (not truncation) for all other coefficients.
      q[i] = abs(c) < abs(lead) ? 0n : floorDivide(c, lead);
      for (let j = 0; j < n; j++) r[i + j] = (r[i + j] ?? 0n) - q[i]! * b[j]!;
    }
    return [normalized(q), multiply(q, b)];
  }
  const n2 = Math.floor(n / 2),
    n1 = n - n2;
  const first = balancedDivision(a.slice(2 * n2), b.slice(n2), exact);
  if (!first) return null;
  const q1 = first[0],
    dq1 = multiply(b, q1);
  const t = subtract(a.slice(n2), dq1);
  const p2 = Array<bigint>(n2 - 1)
    .fill(0n)
    .concat(Array.from({ length: n2 }, (_, i) => t[n1 - 1 + i] ?? 0n));
  const second = balancedDivision(p2, b.slice(n1), exact);
  if (!second) return null;
  const q2 = second[0];
  const q = shifted(q1, n2);
  for (let i = 0; i < q2.length; i++) q[i] = (q[i] ?? 0n) + q2[i]!;
  const product = shifted(dq1, n2),
    dq2 = multiply(b, q2);
  for (let i = 0; i < dq2.length; i++) product[i] = (product[i] ?? 0n) + dq2[i]!;
  return [normalized(q), normalized(product)];
}
/** fmpz_poly_divrem.c and divrem_divconquer.c, with the native exactness flag. */
function divide(
  a: readonly bigint[],
  b: readonly bigint[],
  exact: boolean
): [bigint[], bigint[]] | null {
  if (a.length < b.length) return [[], a.slice()];
  if (b.length < 6) {
    const r = a.slice(),
      q = Array<bigint>(a.length - b.length + 1).fill(0n),
      lead = b[b.length - 1]!;
    for (let i = q.length - 1; i >= 0; i--) {
      const c = r[i + b.length - 1]!;
      if (exact && c % lead !== 0n) return null;
      // FLINT skips coefficients smaller in magnitude than the divisor lead,
      // then uses floor division (not truncation) for all other coefficients.
      q[i] = abs(c) < abs(lead) ? 0n : floorDivide(c, lead);
      for (let j = 0; j < b.length; j++) r[i + j] = r[i + j]! - q[i]! * b[j]!;
    }
    return [normalized(q), normalized(r)];
  }
  if (a.length <= 2 * b.length - 1) {
    const n = a.length - b.length + 1,
      offset = b.length - n;
    const result = balancedDivision(a.slice(offset), b.slice(offset), exact);
    if (!result) return null;
    const product = shifted(result[1], offset),
      low = multiply(b.slice(0, offset), result[0]);
    for (let i = 0; i < low.length; i++) product[i] = (product[i] ?? 0n) + low[i]!;
    return [result[0], subtract(a, product)];
  }
  let r = a.slice(),
    length = a.length;
  const q = Array<bigint>(a.length - b.length + 1).fill(0n),
    block = 2 * b.length - 1;
  while (length >= block) {
    const shift = length - block;
    const result = balancedDivision(r.slice(shift, length), b, exact);
    if (!result) return null;
    for (let i = 0; i < result[0].length; i++) q[shift + i] = result[0][i]!;
    for (let i = 0; i < result[1].length; i++) r[shift + i] = (r[shift + i] ?? 0n) - result[1][i]!;
    length -= b.length;
  }
  if (length >= b.length) {
    const result = divide(r.slice(0, length), b, exact);
    if (!result) return null;
    for (let i = 0; i < result[0].length; i++) q[i] = result[0][i]!;
    for (let i = 0; i < length; i++) r[i] = result[1][i] ?? 0n;
  }
  return [normalized(q), normalized(r)];
}
