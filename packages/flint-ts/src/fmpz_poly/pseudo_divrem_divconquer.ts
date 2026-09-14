/** Native pseudo division with FLINT's degree splits and cutoff 16.
 * @see Deviation: Polynomial Integer Powers and Portable Native Products
 */
import { _fmpz_poly_pseudo_divrem_basecase } from './pseudo_divrem_basecase.js';
import { _fmpz_poly_mul } from './mul.js';
import { _fmpz_poly_resultant_kernels as k } from './gcd.js';

type Division = [bigint[], bigint[], bigint];

/** Public array adapter: L^d A = Q B + R, degree(R) < degree(B). */
export function _fmpz_poly_pseudo_divrem_divconquer(
  a: readonly bigint[],
  b: readonly bigint[]
): Division {
  const A = k.normalized(a),
    B = k.normalized(b);
  if (!B.length) throw new RangeError('division by zero polynomial');
  if (A.length < B.length) return [[], A, 0n];
  const [q, r, d] = divide(A, B);
  return [k.normalized(q), k.normalized(r), d];
}

/** Fixed-length buffers preserve the original recursion even on degree drops. */
function divide(A: readonly bigint[], B: readonly bigint[]): Division {
  const lenA = A.length,
    lenB = B.length;
  if (lenB <= 16 || (lenA > 2 * lenB - 1 && lenA < 128)) {
    return _fmpz_poly_pseudo_divrem_basecase(A, B);
  }
  const n2 = Math.floor(lenB / 2),
    n1 = lenB - n2,
    lead = B[lenB - 1]!;
  if (lenA <= lenB + n2 - 1) {
    const p1 = A.slice(n1);
    p1.fill(0n, 0, n2 - 1);
    const [q, , d] = divide(p1, B.slice(n1));
    // The native buffer moves implement this same low-part subtraction.
    const qb = _fmpz_poly_mul(q, B),
      scale = lead ** d;
    const r = Array.from({ length: lenB - 1 }, (_, i) => A[i]! * scale - (qb[i] ?? 0n));
    return [q, r, d];
  }
  if (lenA > 2 * lenB - 1) {
    const shift = lenA - 2 * lenB + 1,
      p1 = A.slice(shift);
    p1.fill(0n, 0, lenB - 1);
    const [q1, r1, s1] = divide(p1, B);
    const scale1 = lead ** s1,
      t = A.slice(0, lenA - lenB).map((c) => c * scale1);
    for (let i = 0; i < r1.length; i++) t[shift + i] = t[shift + i]! + r1[i]!;
    const [q2, r, s2] = divide(t, B);
    const q = Array<bigint>(lenA - lenB + 1).fill(0n),
      scale2 = lead ** s2;
    for (let i = 0; i < q2.length; i++) q[i] = q2[i]!;
    for (let i = 0; i < q1.length; i++) q[shift + i] = q[shift + i]! + q1[i]! * scale2;
    return [q, r, s1 + s2];
  }
  const p1 = A.slice(2 * n2);
  p1.fill(0n, 0, n1 - 1);
  const [q1, r1, s1] = divide(p1, B.slice(n2));
  const d2q1 = _fmpz_poly_mul(B.slice(0, n2), q1);
  const scale1 = lead ** s1,
    t = A.slice(0, lenB + n2 - 1).map((c) => c * scale1);
  for (let i = 0; i < r1.length; i++) t[2 * n2 + i] = t[2 * n2 + i]! + r1[i]!;
  for (let i = 0; i < d2q1.length; i++) t[n2 + i] = t[n2 + i]! - d2q1[i]!;
  const [q2, r, s2] = divide(t, B);
  const q = Array<bigint>(lenA - lenB + 1).fill(0n),
    scale2 = lead ** s2;
  for (let i = 0; i < q2.length; i++) q[i] = q2[i]!;
  for (let i = 0; i < q1.length; i++) q[n2 + i] = q[n2 + i]! + q1[i]! * scale2;
  return [q, r, s1 + s2];
}
