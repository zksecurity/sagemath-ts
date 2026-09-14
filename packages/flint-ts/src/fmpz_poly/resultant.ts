import { _fmpz_poly_resultant_kernels as k } from './gcd.js';
import { _nmod_poly_resultant } from '../nmod_poly/resultant.js';

/** FLINT integer resultant: signed content, subresultant PRS and balanced modular CRT.
 * @see Deviation: FLINT Resultant Kernels
 */
export function _fmpz_poly_resultant(a: readonly bigint[], b: readonly bigint[]): bigint {
  let A = k.normalized(a),
    B = k.normalized(b);
  if (!A.length || !B.length) return 0n;
  const swapped = A.length < B.length;
  if (swapped) [A, B] = [B, A];
  const r =
    B.length > 144 || B.length ** 3 * (k.maxBits(A) + k.maxBits(B)) > 6000000
      ? modular(A, B)
      : euclidean(A, B);
  return swapped && A.length > 1 && A.length % 2 === 0 && B.length % 2 === 0 ? -r : r;
}
function euclidean(A: bigint[], B: bigint[]): bigint {
  if (B.length === 1) return B[0]! ** BigInt(A.length - 1);
  const ca = k.content(A),
    cb = k.content(B),
    content = ca ** BigInt(B.length - 1) * cb ** BigInt(A.length - 1);
  A = A.map((c) => c / ca);
  B = B.map((c) => c / cb);
  let g = 1n,
    h = 1n,
    sign = 1n;
  do {
    const d = A.length - B.length;
    if (A.length % 2 === 0 && B.length % 2 === 0) sign = -sign;
    const R = k.pseudoRemainder(A, B);
    if (!R.length) return 0n;
    [A, B] = [B, R];
    const hd = h ** BigInt(d);
    B = B.map((c) => c / (g * hd));
    h = (h * A[A.length - 1]! ** BigInt(d)) / hd;
    g = A[A.length - 1]!;
  } while (B.length > 1);
  h = (h * B[0]! ** BigInt(A.length - 1)) / h ** BigInt(A.length - 1);
  return sign * content * h;
}
function modular(A: bigint[], B: bigint[]): bigint {
  if (B.length === 1) return B[0]! ** BigInt(A.length - 1);
  const ca = k.content(A),
    cb = k.content(B);
  A = A.map((c) => c / ca);
  B = B.map((c) => c / cb);
  const lead = A[A.length - 1]! * B[B.length - 1]!;
  const normA = A.reduce((s, c) => s + c * c, 0n),
    normB = B.reduce((s, c) => s + c * c, 0n);
  const norm = normA ** BigInt(B.length - 1) * normB ** BigInt(A.length - 1);
  // Exact bit length of floor(sqrt(norm)) + 1, without materializing its square root.
  const half = Math.ceil(k.bits(norm) / 2),
    edge = (1n << BigInt(half)) - 1n;
  const bound = half + (norm >= edge * edge ? 1 : 0) + 2;
  const residues: [bigint, bigint][] = [];
  for (const p of k.modularPrimes()) {
    if (lead % p === 0n) continue;
    residues.push([
      _nmod_poly_resultant(
        A.map((c) => k.mod(c, p)),
        B.map((c) => k.mod(c, p)),
        p
      ),
      p,
    ]);
    if (residues.length * 63 >= bound) break;
  }
  // Balanced product-tree reconstruction follows fmpz_multi_CRT's asymptotic shape.
  let nodes = residues;
  while (nodes.length > 1) {
    const next: [bigint, bigint][] = [];
    for (let i = 0; i < nodes.length; i += 2) {
      if (i + 1 === nodes.length) {
        next.push(nodes[i]!);
        continue;
      }
      const [a, m] = nodes[i]!,
        [b, n] = nodes[i + 1]!;
      next.push([a + m * k.mod((b - a) * k.inverse(m % n, n), n), m * n]);
    }
    nodes = next;
  }
  let [r, modulus] = nodes[0]!;
  if (r > modulus / 2n) r -= modulus;
  return r * ca ** BigInt(B.length - 1) * cb ** BigInt(A.length - 1);
}
