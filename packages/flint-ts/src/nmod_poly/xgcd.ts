import { _nmod_poly_xgcd_kernels as k } from './gcd.js';

/** Dense modular extended GCD, with FLINT's Euclidean/half-GCD dispatch.
 * Returns monic [g,s,t] satisfying g = s*a + t*b. Inputs are not mutated.
 * @see Deviation: Polynomial Extended GCD Finite Backends
 * @see Reference: nmod_poly/xgcd.c and gr_poly/xgcd_hgcd.c
 */
export function _nmod_poly_xgcd(
  a: readonly bigint[],
  b: readonly bigint[],
  n: bigint
): [bigint[], bigint[], bigint[]] {
  if (n < 2n) throw new RangeError('modulus must be at least 2');
  let A = k.normalized(a, n),
    B = k.normalized(b, n);
  const swapped = A.length < B.length;
  if (swapped) [A, B] = [B, A];
  const finish = (g: bigint[], s: bigint[], t: bigint[]): [bigint[], bigint[], bigint[]] => {
    if (g.length) {
      const inv = k.inverse(g[g.length - 1]!, n);
      const scale = (v: bigint[]) => v.map((c) => (c * inv) % n);
      [g, s, t] = [scale(g), scale(s), scale(t)];
    }
    return swapped ? [g, t, s] : [g, s, t];
  };
  if (!A.length) return [[], [], []];
  if (!B.length) return finish(A, [1n], []);
  if (B.length === 1) return finish(B, [], [1n]);
  const cutoff = CUTOFF[Math.min(n.toString(2).length - 1, 63)]!;
  if (B.length < cutoff) return finish(...euclidean(A, B, n));
  const r0 = k.divrem(A, B, n)[1];
  if (!r0.length) return finish(B, [], [1n]);
  let [M, h, j] = k.halfGcd(B, r0, n, k.innerCutoff(n));
  // Our HGCD matrix maps (B,r0) to (h,j); these are their A-coefficients.
  let S = M[1],
    T = M[3],
    G = h;
  while (j.length) {
    const [q, r] = k.divrem(h, j, n);
    const v = k.mul(q, T, n);
    [S, T] = [T, S];
    if (!r.length) {
      G = j;
      break;
    }
    T = k.sub(T, v, n);
    if (j.length < cutoff) {
      const [g, u0, u1] = euclidean(j, r, n);
      G = g;
      S = k.add(k.mul(S, u0, n), k.mul(T, u1, n), n);
      break;
    }
    [M, h, j] = k.halfGcd(j, r, n, k.innerCutoff(n));
    [S, T] = k.apply(M, S, T, n);
    G = h;
  }
  T = k.divrem(k.sub(G, k.mul(S, A, n), n), B, n)[0];
  return finish(G, S, T);
}

const CUTOFF = [
  159, 89, 75, 85, 72, 78, 81, 89, 85, 89, 93, 97, 97, 93, 101, 97, 97, 106, 106, 121, 111, 101,
  101, 101, 106, 106, 106, 101, 106, 182, 191, 191, 200, 231, 200, 191, 220, 210, 210, 231, 231,
  220, 210, 220, 242, 210, 231, 220, 254, 254, 242, 242, 242, 254, 306, 306, 266, 306, 292, 306,
  321, 306, 388, 166,
];
function euclidean(A: bigint[], B: bigint[], n: bigint): [bigint[], bigint[], bigint[]] {
  if (B.length === 1) return [B, [], [1n]];
  const r = k.divrem(A, B, n)[1];
  if (!r.length) return [B, [], [1n]];
  let D = B,
    U: bigint[] = [],
    V1 = [1n],
    V3 = r;
  do {
    const [q, r] = k.divrem(D, V3, n);
    [U, V1] = [V1, k.sub(U, k.mul(V1, q, n), n)];
    [D, V3] = [V3, r];
  } while (V3.length);
  const T = k.divrem(k.sub(D, k.mul(A, U, n), n), B, n)[0];
  return [D, U, T];
}
