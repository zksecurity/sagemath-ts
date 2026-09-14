import { _ntl_gexteucl as xgcd } from './lip.js';
const exact = (a: bigint, b: bigint) => {
  if (a % b) throw new Error('ExactDiv: nonzero remainder');
  return a / b;
};
const dot = (a: bigint[], b: bigint[]) => a.reduce((s, x, i) => s + x * b[i]!, 0n);
const bal = (a: bigint, d: bigint) => {
  let q = a / d,
    r = a % d;
  if (r < 0n) {
    q--;
    r += d;
  }
  return r * 2n > d || (r * 2n === d && q < 0n) ? q + 1n : q;
};
const transform = (
  A: bigint[],
  B: bigint[],
  x: bigint,
  y: bigint,
  u: bigint,
  v: bigint
): [bigint[], bigint[]] => [A.map((a, i) => x * a + y * B[i]!), A.map((a, i) => u * a + v * B[i]!)];
function reduceKernel(
  input: readonly (readonly bigint[])[],
  a: bigint,
  b: bigint,
  wantU: boolean,
  imageOnly: boolean
): [number, bigint[], bigint[][], bigint[][] | null] {
  const m = input.length,
    B = [[], ...input.map((r) => r.slice())],
    U: bigint[][] | null = wantU
      ? [
          [],
          ...Array.from({ length: m }, (_, i) =>
            Array.from({ length: m }, (_, j) => (i === j ? 1n : 0n))
          ),
        ]
      : null;
  const P = Array<number>(m + 1).fill(0),
    D = Array<bigint>(m + 1).fill(0n),
    lam = Array.from({ length: m + 1 }, () => Array<bigint>(m + 1).fill(0n));
  D[0] = 1n;
  let s = 0,
    k = 1,
    maxK = 0,
    force = true;
  const gs = (k: number) => {
    for (let j = 1; j < k; j++) {
      const pos = P[j]!;
      if (!pos) continue;
      let u = dot(B[k]!, B[j]!);
      for (let i = 1; i < pos; i++) u = (D[i]! * u - lam[k]![i]! * lam[j]![i]!) / D[i - 1]!;
      lam[k]![pos] = u;
    }
    let u = dot(B[k]!, B[k]!);
    for (let i = 1; i <= s; i++) u = (D[i]! * u - lam[k]![i]! * lam[k]![i]!) / D[i - 1]!;
    if (u === 0n) P[k] = 0;
    else {
      P[k] = ++s;
      D[s] = u;
    }
  };
  const reduce = (k: number, l: number) => {
    const p = P[l]!;
    if (!p) return;
    const twice = 2n * lam[k]![p]!;
    if ((twice < 0n ? -twice : twice) <= D[p]!) return;
    const r = bal(lam[k]![p]!, D[p]!);
    B[k] = B[k]!.map((v, i) => v - r * B[l]![i]!);
    if (U) U[k] = U[k]!.map((v, i) => v - r * U[l]![i]!);
    for (let j = 1; j < l; j++) if (P[j]) lam[k]![P[j]!] = lam[k]![P[j]!]! - r * lam[l]![P[j]!]!;
    lam[k]![p] = lam[k]![p]! - r * D[p]!;
  };
  const swap = (k: number): boolean => {
    const pk = P[k]!,
      prev = P[k - 1]!;
    if (pk) {
      [B[k - 1], B[k]] = [B[k]!, B[k - 1]!];
      if (U) [U[k - 1], U[k]] = [U[k]!, U[k - 1]!];
      for (let j = 1; j <= k - 2; j++)
        if (P[j]) [lam[k - 1]![P[j]!], lam[k]![P[j]!]] = [lam[k]![P[j]!]!, lam[k - 1]![P[j]!]!];
      const l = lam[k]![pk - 1]!;
      for (let i = k + 1; i <= maxK; i++) {
        const c1 = lam[i]![pk - 1]!,
          c2 = lam[i]![pk]!;
        lam[i]![pk - 1] = exact(l * c1 + D[pk - 2]! * c2, D[pk - 1]!);
        lam[i]![pk] = exact(D[pk]! * c1 - l * c2, D[pk - 1]!);
      }
      D[pk - 1] = exact(D[pk - 2]! * D[pk]! + l * l, D[pk - 1]!);
      return false;
    }
    if (lam[k]![prev] !== 0n) {
      const [x, y, e] = xgcd(lam[k]![prev]!, D[prev]!);
      const t1 = exact(lam[k]![prev]!, e),
        t3 = exact(D[prev]!, e),
        t2 = -t3;
      [B[k - 1], B[k]] = transform(B[k - 1]!, B[k]!, t1, t2, y, x);
      if (U) [U[k - 1], U[k]] = transform(U[k - 1]!, U[k]!, t1, t2, y, x);
      for (let j = 1; j <= k - 2; j++)
        if (P[j]) {
          const p = P[j]!,
            c1 = lam[k - 1]![p]!,
            c2 = lam[k]![p]!;
          lam[k - 1]![p] = t1 * c1 + t2 * c2;
          lam[k]![p] = y * c1 + x * c2;
        }
      const denom = t2 * t2;
      D[prev] = exact(D[prev]!, denom);
      for (let i = k + 1; i <= maxK; i++)
        if (P[i]) {
          const p = P[i]!;
          D[p] = exact(D[p]!, denom);
          for (let j = i + 1; j <= maxK; j++) lam[j]![p] = exact(lam[j]![p]!, denom);
        }
      for (let i = k + 1; i <= maxK; i++) lam[i]![prev] = exact(lam[i]![prev]!, t3);
      [P[k - 1], P[k]] = [P[k]!, P[k - 1]!];
      return true;
    }
    [B[k - 1], B[k]] = [B[k]!, B[k - 1]!];
    if (U) [U[k - 1], U[k]] = [U[k]!, U[k - 1]!];
    for (let j = 1; j <= k - 2; j++)
      if (P[j]) [lam[k - 1]![P[j]!], lam[k]![P[j]!]] = [lam[k]![P[j]!]!, lam[k - 1]![P[j]!]!];
    [P[k - 1], P[k]] = [P[k]!, P[k - 1]!];
    return false;
  };
  while (k <= m) {
    if (k > maxK) {
      gs(k);
      maxK = k;
    }
    if (k === 1) {
      force = true;
      k++;
      continue;
    }
    if (force) for (let j = k - 1; j >= 1; j--) reduce(k, j);
    const p = P[k]!,
      mustSwap =
        P[k - 1] !== 0 &&
        (p === 0 ||
          (!imageOnly &&
            a * D[p - 1]! * D[p - 1]! >
              b * (D[p]! * D[p - 2]! + lam[k]![p - 1]! * lam[k]![p - 1]!)));
    if (mustSwap) {
      force = swap(k);
      k--;
    } else {
      force = true;
      k++;
    }
  }
  return [s, D.slice(0, s + 1), B.slice(1), U?.slice(1) ?? null];
}
/** NTL LLL.cpp LLL_plus; row matrices and copied output-reference results.
 * @see Deviation: NTL exact integer lattice adapters
 */
export function LLL_plus(
  B: readonly (readonly bigint[])[],
  a = 3,
  b = 4,
  transform = false
): [number, bigint[], bigint[][], bigint[][] | null] {
  if (a <= 0 || b <= 0 || a > b || Math.trunc(b / 4) >= a) throw new Error('LLL_plus: bad args');
  return reduceKernel(B, BigInt(a), BigInt(b), transform, false);
}
/** NTL LLL.cpp LLL; row matrices and copied output-reference results.
 * @see Deviation: NTL exact integer lattice adapters
 */
export function LLL(
  B: readonly (readonly bigint[])[],
  a = 3,
  b = 4,
  transform = false
): [number, bigint, bigint[][], bigint[][] | null] {
  if (a <= 0 || b <= 0 || a > b || Math.trunc(b / 4) >= a) throw new Error('LLL: bad args');
  const [r, D, M, U] = reduceKernel(B, BigInt(a), BigInt(b), transform, false);
  return [r, D[r]!, M, U];
}
/** NTL LLL.cpp image; row matrices and copied output-reference results.
 * @see Deviation: NTL exact integer lattice adapters
 */
export function image(
  B: readonly (readonly bigint[])[],
  transform = false
): [number, bigint, bigint[][], bigint[][] | null] {
  const [r, D, M, U] = reduceKernel(B, 3n, 4n, transform, true);
  return [r, D[r]!, M, U];
}
/** NTL LLL.cpp LatticeSolve; row matrices and copied output-reference results.
 * @see Deviation: NTL exact integer lattice adapters
 */
export function LatticeSolve(
  A: readonly (readonly bigint[])[],
  y: readonly bigint[],
  reduce = 0,
  initial: readonly bigint[] = [],
  columns = A[0]?.length ?? 0
): [number, bigint[]] {
  const n = A.length;
  if (y.length !== columns) throw new Error('LatticeSolve: dimension mismatch');
  if (reduce < 0 || reduce > 2) throw new Error('LatticeSolve: bad reduce parameter');
  if (y.every((c) => c === 0n)) return [1, Array<bigint>(n).fill(0n)];
  const [rank, , A1, U1] = image(A, true),
    ker = n - rank;
  const [newRank, , , U2] = image([...A1.slice(ker), y.slice()], true);
  if (newRank !== rank || (U2![0]![rank] !== 1n && U2![0]![rank] !== -1n))
    return [0, initial.slice()];
  let x1 = U2![0]!.slice(0, rank);
  if (U2![0]![rank] === 1n) x1 = x1.map((c) => -c);
  let x2 = Array<bigint>(n).fill(0n);
  for (let i = 0; i < rank; i++) x2 = x2.map((c, j) => c + U1![ker + i]![j]! * x1[i]!);
  if (reduce === 0) return [1, x2];
  let K = U1!.slice(0, ker);
  if (reduce === 2) K = LLL(K)[2];
  const result = image([...K, x2])[2];
  return [1, result[ker]!];
}
