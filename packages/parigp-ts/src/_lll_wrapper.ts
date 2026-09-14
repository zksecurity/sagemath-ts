/** PARI lll.c adaptive integer LLL and iterated compression.
 * @see Deviation: PARI adaptive LLL and FLATTER adapters
 */
import { fplll_fast, fplll_heuristic, fplll_dpe, fplll } from './lll.js';
import { gramschmidt_dynprec, RgM_Cholesky_dynprec, drop, potential, spread } from './_lll_gso.js';
import { type MpReal, itor, mkqfb, redimagsl2 } from './qfb.js';
import { gram_matrix, RgM_mul, type RgScalar } from './RgV.js';
import { RgM_rescale_to_int } from './polarit2.js';
import { RgM_inv_upper } from './alglin1.js';
import { RgM_Cholesky } from './alglin2.js';
import { QR_init } from './bibli1.js';
import { integerMatrixInverse, integerMatrixPivots } from './_matrix_inverse.js';
import { roundr_safe } from './gen3.js';
import { type MatrixReal } from './_real_matrix.js';
import { abscmprr } from './kernel/none/cmp.js';
import { PariError } from './errors.js';

type Matrix = bigint[][];
type RealMatrix = RgScalar[][];
export type LllResult = Matrix | [Matrix, Matrix];
export const LLL_KER = 1,
  LLL_IM = 2,
  LLL_ALL = 4,
  LLL_GRAM = 0x100,
  LLL_KEEP_FIRST = 0x200,
  LLL_INPLACE = 0x400,
  LLL_COMPATIBLE = 0x800,
  LLL_UPPER = 0x1000,
  LLL_NOCERTIFY = 0x2000,
  LLL_NOFLATTER = 0x4000;
const identity = (n: number): Matrix =>
  Array.from({ length: n }, (_, j) => Array.from({ length: n }, (_, i) => (i === j ? 1n : 0n)));
const copy = <T>(x: T[][]): T[][] => x.map((c) => c.slice());
const transpose = <T>(x: T[][]): T[][] => x[0]?.map((_, i) => x.map((c) => c[i]!)) ?? [];
const mul = (x: Matrix, y: Matrix): Matrix => RgM_mul(x, y) as Matrix;
const real = (x: Matrix, p: number): MpReal[][] => x.map((c) => c.map((v) => itor(v, p)));
const upper = (x: Matrix): boolean =>
  x.length === (x[0]?.length ?? 0) && x.every((c, j) => c.every((v, i) => i <= j || v === 0n));
const lower = (x: Matrix): boolean => upper(transpose(x));
const knapsack = (x: Matrix): boolean =>
  x.length === (x[0]?.length ?? 0) &&
  x.every((c, j) => c.every((v, i) => i === 0 || i === j || v === 0n));
const exp = (x: bigint): bigint =>
  x === 0n ? -(1n << 61n) : BigInt((x < 0n ? -x : x).toString(2).length - 1);
const normExp = (x: Matrix): bigint =>
  exp(x.reduce((s, c) => s + c.reduce((t, v) => t + v * v, 0n), 0n));
const qfApply = (x: Matrix, t: Matrix): Matrix => mul(transpose(t), mul(x, t));
function augment(x: Matrix, i: number): Matrix {
  const I = identity(x.length);
  return x.map((c, j) => [...c.map((v) => v << BigInt(i)), ...I[j]!]);
}
const shiftGram = (x: Matrix, i: number): Matrix =>
  x.map((c, j) => c.map((v, k) => (v << BigInt(i)) + (j === k ? 1n : 0n)));
function trivial<T extends RgScalar>(
  x: T[][],
  flag: number
): T[][] | Matrix | [Matrix, T[][] | Matrix] {
  const n = x.length;
  if (!n) return flag & LLL_ALL ? [[], []] : [];
  if (x[0]!.every((v) => (typeof v === 'bigint' ? v === 0n : !v.s))) {
    if (flag & LLL_KER) return identity(1);
    if (flag & (LLL_IM | LLL_INPLACE)) return [];
    return [identity(1), []];
  }
  if (flag & LLL_INPLACE) return copy(x);
  if (flag & LLL_KER) return [];
  if (flag & LLL_IM) return identity(1);
  return [[], flag & LLL_GRAM ? copy(x) : identity(1)];
}
function finish(h: Matrix, k: number, flag: number): LllResult {
  if (!(flag & (LLL_IM | LLL_KER | LLL_ALL | LLL_INPLACE))) return h;
  if (flag & (LLL_IM | LLL_INPLACE)) return h.slice(k);
  if (flag & LLL_KER) return h.slice(0, k);
  return [h.slice(0, k), h.slice(k)];
}
function sizered(T1: Matrix, T3: Matrix, R1: MatrixReal[][], R2: MatrixReal[][]): Matrix {
  const inverse = integerMatrixInverse(transpose(T1))!;
  const M = RgM_mul(transpose(inverse[0]), RgM_mul(RgM_mul(RgM_inv_upper(R1), R2), T3));
  return mul(
    T1.map((c) => c.map((v) => -v)),
    M.map((c) => c.map((v) => (typeof v === 'bigint' ? v : roundr_safe(v))))
  );
}
/** Native one-step FLATTER, including statistics before the final overlap.
 * @see Deviation: PARI adaptive LLL and FLATTER adapters
 */
export function flat(M: Matrix, flag = 0): [Matrix, Matrix | null, bigint, bigint] {
  if (M.length < 2) throw new RangeError('FLATTER requires at least two columns');
  const k = M.length,
    n = k >> 1,
    n2 = k - n,
    m = n >> 1;
  let R = gramschmidt_dynprec(M) as MatrixReal[][];
  const slice = (x: MatrixReal[][], row: number, rows: number, col: number, cols: number) =>
    x.slice(col, col + cols).map((c) => c.slice(row, row + rows));
  const subflag = LLL_IM | LLL_UPPER | LLL_NOCERTIFY;
  const R1 = slice(R, 0, n, 0, n),
    R2 = slice(R, 0, n, n, n2),
    R3 = slice(R, n, n2, n, n2);
  const T1 = lllfp(R1, 0.99, subflag | (flag & LLL_KEEP_FIRST)) as Matrix;
  let T3 = lllfp(R3, 0.99, subflag) as Matrix;
  const T2 = sizered(T1, T3, R1, R2);
  const T: Matrix = [
    ...T1.map((c) => [...c, ...Array<bigint>(n2).fill(0n)]),
    ...T3.map((c, j) => [...T2[j]!, ...c]),
  ];
  M = mul(M, T);
  R = gramschmidt_dynprec(M) as MatrixReal[][];
  T3 = lllfp(slice(R, m, n2, m, n2), 0.99, subflag) as Matrix;
  const S = identity(k);
  for (let j = 0; j < n2; j++) for (let i = 0; i < n2; i++) S[m + j]![m + i] = T3[j]![i]!;
  return [mul(M, S), flag & LLL_INPLACE ? null : mul(T, S), drop(R), potential(R)];
}
/** Native iterated FLATTER; a rejected final step is deliberately discarded.
 * @see Deviation: PARI adaptive LLL and FLATTER adapters
 */
export function ZM_flatter(M: Matrix, flag = 0): Matrix | null {
  let s = -1n,
    pot = (1n << 63n) - 1n,
    T: Matrix | null = null;
  for (let i = 1; ; i++) {
    const [M2, U, t, pot2] = flat(M, flag);
    if (t === 0n || (s >= 0n && ((s === t && pot >= pot2) || (s < t && i > 20)))) break;
    s = t;
    pot = pot2;
    M = M2;
    if (!(flag & LLL_INPLACE)) T = T ? mul(T, U!) : U;
  }
  return flag & LLL_INPLACE ? M : T;
}
/** Native rank-deficient augmentation and transformation growth certificate.
 * @see Deviation: PARI adaptive LLL and FLATTER adapters
 */
export function ZM_flatter_rank(M: Matrix, rank: number, flag = 0): Matrix | null {
  const n = M.length;
  if (rank === n) return ZM_flatter(M, flag);
  let T: Matrix | null = null,
    sm = (1n << 63n) - 1n;
  for (let i = 1; ; i++) {
    const S = ZM_flatter(augment(M, i), flag);
    if (!S) break;
    const s = normExp(S);
    if (s >= sm) break;
    sm = s;
    // ZM_mul reads only its left operand's number of columns from each RHS
    // column. The native rank-deficient INPLACE path passes augmented columns.
    const truncated = S.map((c) => c.slice(0, n));
    T = T ? mul(T, truncated) : S;
    M = mul(M, truncated);
  }
  return flag & LLL_INPLACE ? M : (T ?? identity(n));
}
function flattergramStep(M: Matrix, flag: number): Matrix {
  return lllfp(
    RgM_Cholesky_dynprec(M),
    0.99,
    LLL_IM | LLL_UPPER | LLL_NOCERTIFY | (flag & LLL_KEEP_FIRST)
  ) as Matrix;
}
/** Native Gram compression on a positive-definite form of dimension >= 2.
 * @see Deviation: PARI adaptive LLL and FLATTER adapters
 */
export function ZM_flattergram(M: Matrix, flag = 0): Matrix | null {
  if (M.length < 2) throw new RangeError('Gram FLATTER requires at least two columns');
  let T: Matrix | null = null,
    s = -1n;
  for (let i = 1; ; i++) {
    const S = flattergramStep(M, flag),
      t = normExp(S);
    if (t === 0n || (s !== 0n && (s === t || (s < t && i > 20)))) break;
    T = T ? mul(T, S) : S;
    M = qfApply(M, S);
    s = t;
  }
  return T;
}
function flatterGramRank(M: Matrix, rank: number, flag: number): Matrix | null {
  if (rank === M.length) return ZM_flattergram(M, flag);
  let T: Matrix | null = null;
  for (let i = 1; ; i++) {
    // The bundled loop has no successful exit for a nonsingular transform in
    // dimension > 1. Bound this resource failure without inventing a result;
    // native arithmetic exceptions before the limit remain observable.
    if (i > 64)
      throw new RangeError('rank-deficient Gram FLATTER exceeded 64 augmentation attempts');
    const S = ZM_flattergram(shiftGram(M, i), flag);
    if (!S) break;
    T = T ? mul(T, S) : S;
    M = qfApply(M, S);
  }
  return !T || T.every((c, j) => c.every((v, i) => v === (i === j ? 1n : 0n))) ? null : T;
}
function initialGso(M: Matrix, rank: number, gram: boolean): RealMatrix | null {
  const n = M.length,
    p = Math.ceil((3 * n + 30) / 64) * 64;
  if (gram) return RgM_Cholesky(real(rank < n ? shiftGram(M, 1) : M, p), p);
  const [ok, , , L] = QR_init(real(rank < n ? augment(M, 1) : M, p), p);
  return ok && L && L.every((c, j) => (typeof c[j] === 'bigint' ? c[j] !== 0n : c[j]!.s !== 0))
    ? L
    : null;
}
function supExponent(R: RealMatrix): bigint {
  let m: RgScalar = R[0]![0]!;
  for (const c of R)
    for (const v of c) {
      if (typeof m === 'bigint' && typeof v === 'bigint') {
        if ((v < 0n ? -v : v) > (m < 0n ? -m : m)) m = v;
      } else if (typeof m !== 'bigint' && typeof v !== 'bigint' && abscmprr(v, m) > 0) m = v;
    }
  return typeof m === 'bigint' ? exp(m) : BigInt(m.e);
}
// Literal tuning tables from the bundled lll.c; changing them changes U.
const thre = [
  31783, 34393, 20894, 22525, 13533, 1928, 672, 671, 422, 506, 315, 313, 222, 205, 167, 154, 139,
  138, 110, 120, 98, 94, 81, 75, 74, 64, 74, 74, 79, 96, 112, 111, 105, 104, 96, 86, 84, 78, 75, 70,
  66, 62, 62, 57, 56, 47, 45, 52, 50, 44, 48, 42, 36, 35, 35, 34, 40, 33, 34, 32, 36, 31, 38, 38,
  40, 38, 38, 37, 35, 31, 34, 36, 34, 32, 34, 32, 28, 27, 25, 31, 25, 27, 28, 26, 25, 21, 21, 25,
  25, 22, 21, 24, 24, 22, 21, 23, 22, 22, 22, 22, 21, 24, 21, 22, 19, 20, 19, 20, 19, 19, 19, 18,
  19, 18, 18, 20, 19, 20, 18, 19, 18, 21, 18, 20, 18, 18,
];
const thsn = [
  23280, 30486, 50077, 44136, 78724, 15690, 1801, 1611, 981, 1359, 978, 1042, 815, 866, 788, 775,
  726, 712, 626, 613, 548, 564, 474, 481, 504, 447, 453, 508, 705, 794, 1008, 946, 767, 898, 886,
  763, 842, 757, 725, 774, 639, 655, 705, 627, 635, 704, 511, 613, 583, 595, 568, 640, 541, 640,
  567, 540, 577, 584, 546, 509, 526, 572, 637, 746, 772, 743, 743, 742, 800, 708, 832, 768, 707,
  692, 692, 768, 696, 635, 709, 694, 768, 719, 655, 569, 590, 644, 685, 623, 627, 720, 633, 636,
  602, 635, 575, 631, 642, 647, 632, 656, 573, 511, 688, 640, 528, 616, 511, 559, 601, 620, 635,
  688, 608, 768, 658, 582, 644, 704, 555, 673, 600, 601, 641, 661, 601, 670,
];
/** Adaptive native integer LLL. Matrices are zero-based columns; input is copied.
 * @see Deviation: PARI adaptive LLL and FLATTER adapters
 */
function lllWithNorms(
  input: Matrix,
  delta: number,
  flag: number,
  wantNorms: boolean
): [LllResult, MpReal[] | null] {
  const n = input.length,
    m = input[0]?.length ?? 0,
    keep = !!(flag & LLL_KEEP_FIRST);
  if (input.some((c) => c.length !== m))
    throw new RangeError('ZM_lll requires rectangular columns');
  if (n <= 1) return [trivial(input, flag) as LllResult, null];
  if (!m) {
    if (flag & LLL_KER) return [identity(n), null];
    if (flag & (LLL_INPLACE | LLL_IM)) return [[], null];
    return [[identity(n), []], null];
  }
  if (flag & LLL_GRAM && m !== n) throw new RangeError('ZM_lll requires a square Gram matrix');
  if (n === 2 && m === 2 && flag & LLL_IM && !keep) {
    const G = flag & LLL_GRAM ? input : (gram_matrix(input) as Matrix);
    const a = G[0]![0]!,
      b = G[1]![0]! * 2n,
      c = G[1]![1]!,
      D = b * b - 4n * a * c;
    if (D < 0n) {
      // Bundled lll.c passes a t_QFB to RgM_gram_schmidt in this path.
      // Preserve the caught native error without causing an actual memory fault.
      if (wantNorms) throw new PariError('bug in PARI/GP (Segmentation Fault), please report');
      const U = transpose(redimagsl2(mkqfb(a, b, c, D)).U);
      return [flag & LLL_INPLACE ? mul(input, U) : U, null];
    }
  }
  let B: Matrix | null = flag & LLL_GRAM ? null : copy(input),
    G: Matrix | null = flag & LLL_GRAM ? copy(input) : null;
  let U: Matrix | null = flag & LLL_GRAM || !(flag & LLL_INPLACE) ? identity(n) : null;
  const isUpper = B ? !!(flag & LLL_UPPER) || upper(B) : false,
    isLower = B && !isUpper && !keep ? lower(B) : false;
  let L = isLower
    ? copy(B!)
        .reverse()
        .map((c) => c.reverse())
    : null;
  // The pivot helper returns nullity, while the native wrapper needs rank.
  const actualRank = flag & LLL_NOFLATTER ? 0 : n - integerMatrixPivots(transpose(input), m, n)[1];
  let useflatter = false;
  if (n > 2 && !(flag & LLL_NOFLATTER)) {
    const R = B
      ? isUpper
        ? B
        : isLower
          ? L
          : initialGso(B, actualRank, false)
      : initialGso(G!, actualRank, true);
    if (R) {
      const spr = spread(R);
      let sz = supExponent(R),
        thr: number;
      if ((isUpper && knapsack(B!)) || (isLower && knapsack(L!)))
        thr = thsn[Math.min(n - 3, thsn.length - 1)]!;
      else {
        thr = thre[Math.min(n - 3, thre.length - 1)]!;
        if (n >= 10) sz = spr;
      }
      useflatter = sz >= BigInt(thr);
    } else useflatter = true;
  }
  if (useflatter) {
    if (B) {
      let basis = isLower ? L! : B;
      const T = ZM_flatter_rank(basis, actualRank, flag | (isUpper || isLower ? LLL_UPPER : 0));
      if (T) {
        if (U) {
          U = mul(U, T);
          basis = mul(basis, T);
        } else basis = T;
      }
      if (isLower) {
        B = basis.map((c) => c.slice().reverse());
        if (U) U = U.map((c) => c.slice().reverse());
      } else B = basis;
    } else {
      for (let i = 0; i < n; i++) for (let j = 0; j < i; j++) G![j]![i] = G![i]![j]!;
      const T = flatterGramRank(G!, actualRank, flag);
      if (T) {
        if (U) U = mul(U, T);
        G = qfApply(G!, T);
      }
    }
  }
  let zeros = -1;
  let norms: MpReal[] | null = null;
  const applyU = (T: Matrix | null) => {
    if (U && T) U = mul(U, T);
  };
  if (B) {
    let stage = fplll_fast(B, delta, 0.51, keep, !!U);
    zeros = stage[0];
    B = stage[1];
    applyU(stage[2]);
    for (let p = 64, t = 0; zeros < 0 && t < (n < 100 ? 1 : 2); p += 64, t++) {
      stage = fplll_heuristic(B, delta, 0.51, keep, !!U, p, p);
      zeros = stage[0];
      B = stage[1];
      applyU(stage[2]);
    }
  }
  if (zeros < 0 || !(flag & LLL_NOCERTIFY)) {
    let stage = fplll_dpe(B, G, delta, 0.51, keep, !!U, wantNorms);
    zeros = stage[0];
    G = stage[1];
    B = stage[2];
    applyU(stage[3]);
    norms = stage[4];
    if (zeros < 0)
      for (let p = 64; ; p += 64) {
        stage = fplll(B, G, delta, 0.51, keep, !!U, wantNorms, p);
        zeros = stage[0];
        G = stage[1];
        B = stage[2];
        applyU(stage[3]);
        norms = stage[4];
        if (zeros >= 0) break;
      }
  }
  return [finish(U ?? B!, zeros, flag), norms];
}
/** Native ZM_lll, without requesting Gram--Schmidt norms. */
export function ZM_lll(input: Matrix, delta = 0.99, flag = LLL_IM): LllResult {
  return lllWithNorms(input, delta, flag, false)[0];
}
/** Native ZM_lll_norms. Null norms mean the selected native path leaves pN unset.
 * @see Deviation: PARI adaptive LLL norm output
 */
export function ZM_lll_norms(
  input: Matrix,
  delta = 0.99,
  flag = LLL_IM
): [LllResult, MpReal[] | null] {
  return lllWithNorms(input, delta, flag, true);
}
/** Integer/real lllfp adapter; full inexact Gram input must use real coefficients.
 * @see Deviation: PARI adaptive LLL and FLATTER adapters
 */
export function lllfp(
  input: RealMatrix,
  delta = 0.99,
  flag = LLL_IM
): LllResult | RealMatrix | [Matrix, RealMatrix] | null {
  if (input.length <= 1) return trivial(input, flag);
  let x = input;
  if (flag & LLL_GRAM) {
    if (x.some((c) => c.length !== x.length))
      throw new PariError('inconsistent dimensions in qflllgram');
    const reals = x.flat().filter((v): v is MpReal => typeof v !== 'bigint');
    if (reals.length) {
      if (x.some((c) => c.some((v) => typeof v === 'bigint' && v !== 0n)))
        throw new RangeError('inexact Gram LLL requires real coefficients or exact zero');
      const R = RgM_Cholesky(x as MatrixReal[][], Math.min(...reals.map((v) => v.p)));
      if (!R) return null;
      x = R;
      flag &= ~LLL_GRAM;
    }
  }
  return ZM_lll(RgM_rescale_to_int(x), delta, flag);
}
