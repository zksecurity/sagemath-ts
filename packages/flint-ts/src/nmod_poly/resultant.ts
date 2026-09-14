import { _nmod_poly_resultant_kernels as k } from './gcd.js';

/** FLINT modular resultant, including small-formula and recursive half-GCD dispatch.
 * Same-array identity preserves the native pointer-alias shortcut, including constants.
 * @see Deviation: FLINT Resultant Kernels
 */
export function _nmod_poly_resultant(
  a: readonly bigint[],
  b: readonly bigint[],
  n: bigint
): bigint {
  if (n < 2n) throw new RangeError('modulus must be at least 2');
  let A = k.normalized(a, n),
    B = k.normalized(b, n);
  if (!A.length || !B.length || a === b) return 0n;
  const swapped = A.length < B.length;
  if (swapped) [A, B] = [B, A];
  let r: bigint;
  if (B.length === 1) r = k.power(B[0]!, A.length - 1, n);
  else if (B.length === 2) {
    // gr_poly/resultant_small.c evaluates the homogeneous form without inverses.
    const powers = [1n];
    for (let i = 1; i < A.length; i++) powers.push((powers[i - 1]! * B[0]!) % n);
    let factor = 1n;
    r = 0n;
    for (let i = A.length - 1; i >= 0; i--) {
      r = (r + A[i]! * powers[i]! * factor) % n;
      factor = k.mod(-factor * B[1]!, n);
    }
  } else if (A.length === 3 && B.length === 3) {
    const t0 = A[0]! * A[2]!,
      t1 = B[0]! * B[2]!,
      t2 = A[0]! * B[2]!,
      t3 = A[2]! * B[0]!,
      t7 = A[1]! * B[1]!;
    r = k.mod(t0 * (B[1]! ** 2n - 2n * t1) + t1 * A[1]! ** 2n + t2 * (t2 - t7) + t3 * (t3 - t7), n);
  } else if (B.length < k.gcdCutoff(n)) r = euclidean(A, B, n);
  else {
    let G = A,
      J = B;
    r = 1n;
    for (;;) {
      const R = k.divrem(G, J, n)[1];
      const state = { value: r, lc: J[J.length - 1]!, len0: G.length, len1: J.length, off: 0 };
      k.resultantStep(state, G.length, J.length, R.length, n);
      r = state.value;
      if (!R.length) break;
      // The initial reduction is always followed by HGCD in the native algorithm.
      if (G !== A && J.length < k.outerCutoff(n)) {
        r = (r * euclidean(J, R, n)) % n;
        break;
      }
      [r, G, J] = k.resultantHalfGcd(J, R, n, r);
      if (!J.length) break;
    }
  }
  return swapped && A.length % 2 === 0 && B.length % 2 === 0 ? k.mod(-r, n) : r;
}
function euclidean(A: bigint[], B: bigint[], n: bigint): bigint {
  if (B.length === 1) return k.power(B[0]!, A.length - 1, n);
  const state = { value: 1n, lc: 1n, len0: A.length, len1: B.length, off: 0 };
  while (B.length) {
    const R = k.divrem(A, B, n)[1];
    state.lc = B[B.length - 1]!;
    k.resultantStep(state, A.length, B.length, R.length, n);
    [A, B] = [B, R];
  }
  return state.value;
}
