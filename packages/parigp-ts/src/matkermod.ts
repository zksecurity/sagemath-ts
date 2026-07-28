/**
 * @module parigp-ts/matkermod
 * @description Linear algebra over Z/dZ for **arbitrary** (possibly composite) d,
 * via black-box Hermite rings and the Howell normal form.
 *
 * Direct port of `reference/pari/src/basemath/bb_hnf.c` (contributed to PARI by
 * Aurel Page, 2017), specialised to the `Fp_hermite` black box (the ring Z/dZ):
 *
 * - `bb_hnf.c:66-171`   `_Fp_add/_Fp_neg/_Fp_mul/_Fp_rann/_Fp_lquo/_Fp_unit/_Fp_extgcd`
 *                       plus the helpers `Z_split` and `Z_stab`.
 * - `bb_hnf.c:492-772`  `gen_howell_i` (triangularisation + reduced Howell form)
 * - `bb_hnf.c:826-885`  `gen_kernel_howell`, `gen_kernel_from_howell`, `gen_kernel`
 * - `bb_hnf.c:806-824`  `gen_matimage`
 * - `bb_hnf.c:1036-1050` `matkermod`
 * - `bb_hnf.c:981-990`  `matimagemod`
 * - `bb_hnf.c:1015-1033` `matdetmod`
 * - `bb_hnf.c:1053-1071` `matinvmod`
 *
 * SageMath delegates `Matrix_integer_dense.right_kernel_matrix()` over Z/nZ with
 * composite n to PARI's `matkermod` (`reference/sage/src/sage/matrix/matrix2.pyx`,
 * `right_kernel_matrix`), which is why this lives in `parigp-ts`.
 *
 * MATRIX LAYOUT (important): matrices use PARI's own column-major layout, i.e.
 * `A[j][i]` is the entry in row `i`, column `j`, with **both indices 0-based** at
 * the public API boundary. Use {@link zm_from_rows} / {@link zm_to_rows} to convert
 * from/to the more usual row-major `bigint[][]`.
 * The kernel/image bases are returned as the **columns** of the result, exactly as
 * in PARI.
 */

/**
 * Matrix over Z in PARI's column-major layout: `A[j][i]` = entry at row `i`,
 * column `j` (0-based). `A.length` is the number of columns.
 */
export type ZM = bigint[][];

/** Convert a row-major matrix (`rows[i][j]` = row i, column j) to a {@link ZM}. */
export function zm_from_rows(rows: readonly (readonly bigint[])[]): ZM {
  const m = rows.length;
  const n = m === 0 ? 0 : rows[0]!.length;
  for (const r of rows) {
    if (r.length !== n) throw new Error('zm_from_rows: ragged matrix');
  }
  const out: ZM = [];
  for (let j = 0; j < n; j++) {
    const col: bigint[] = new Array(m);
    for (let i = 0; i < m; i++) col[i] = rows[i]![j]!;
    out.push(col);
  }
  return out;
}

/** Convert a {@link ZM} back to row-major form. */
export function zm_to_rows(A: ZM): bigint[][] {
  const n = A.length;
  const m = n === 0 ? 0 : A[0]!.length;
  const out: bigint[][] = [];
  for (let i = 0; i < m; i++) {
    const row: bigint[] = new Array(n);
    for (let j = 0; j < n; j++) row[j] = A[j]![i]!;
    out.push(row);
  }
  return out;
}

/* ------------------------------------------------------------------ */
/* Errors (mirroring PARI's pari_err_*)                                */
/* ------------------------------------------------------------------ */

/*
 * NOTE: these classes are the shared PARI error kinds used by this module and
 * by `ffinit.ts` / `qfb.ts`. They live here only because `parigp-ts` has no
 * shared errors module yet; they should move to one (e.g. `src/errors.ts`) the
 * day such a file exists, with these names re-exported for compatibility.
 */

/** PARI `pari_err_TYPE` */
export class PariTypeError extends Error {
  constructor(fun: string, what: string) {
    super(`incorrect type in ${fun} (${what})`);
    this.name = 'PariTypeError';
  }
}
/** PARI `pari_err_DOMAIN` */
export class PariDomainError extends Error {
  constructor(fun: string, v: string, op: string, lim: string) {
    super(`domain error in ${fun}: ${v} ${op} ${lim}`);
    this.name = 'PariDomainError';
  }
}
/** PARI `pari_err_DIM` */
export class PariDimError extends Error {
  constructor(fun: string) {
    super(`inconsistent dimensions in ${fun}`);
    this.name = 'PariDimError';
  }
}
/** PARI `pari_err_INV` */
export class PariInvError extends Error {
  constructor(fun: string) {
    super(`impossible inverse in ${fun}`);
    this.name = 'PariInvError';
  }
}
/** PARI `pari_err_PRIME` */
export class PariPrimeError extends Error {
  constructor(fun: string, p: bigint) {
    super(`not a prime number in ${fun}: ${p}`);
    this.name = 'PariPrimeError';
  }
}
/** PARI `pari_err_SQRTN` */
export class PariSqrtnError extends Error {
  constructor(fun: string, what: string) {
    super(`not an n-th power residue in ${fun}: ${what}`);
    this.name = 'PariSqrtnError';
  }
}
/** PARI `pari_err_FLAG` */
export class PariFlagError extends Error {
  constructor(fun: string) {
    super(`incorrect flag in ${fun}`);
    this.name = 'PariFlagError';
  }
}

/* ------------------------------------------------------------------ */
/* Internal 1-indexed representation (mirrors PARI's GEN indexing)     */
/* ------------------------------------------------------------------ */

/** Column with a dummy slot at index 0; entries live at 1..m. */
type Col = bigint[];
/** Matrix with a dummy slot at index 0; columns live at 1..n. */
type Mat = Col[];

const DUMMY = 0n;

function toInternal(A: ZM): { H: Mat; m: number; n: number } {
  const n = A.length;
  const m = n === 0 ? 0 : A[0]!.length;
  const H: Mat = [[]];
  for (let j = 0; j < n; j++) {
    const src = A[j]!;
    if (src.length !== m) throw new PariDimError('matkermod');
    const col: Col = new Array(m + 1);
    col[0] = DUMMY;
    for (let i = 1; i <= m; i++) col[i] = src[i - 1]!;
    H.push(col);
  }
  return { H, m, n };
}

function toExternal(H: Mat, m: number): ZM {
  const out: ZM = [];
  for (let j = 1; j < H.length; j++) {
    const col = H[j]!;
    const dst: bigint[] = new Array(m);
    for (let i = 1; i <= m; i++) dst[i - 1] = col[i]!;
    out.push(dst);
  }
  return out;
}

/* ------------------------------------------------------------------ */
/* Integer helpers                                                     */
/* ------------------------------------------------------------------ */

function iabs(x: bigint): bigint {
  return x < 0n ? -x : x;
}

/** gcd >= 0 */
function gcdii(a: bigint, b: bigint): bigint {
  a = iabs(a);
  b = iabs(b);
  while (b) {
    const t = a % b;
    a = b;
    b = t;
  }
  return a;
}

/** PARI `bezout`: returns `[d,u,v]` with `d = gcd(x,y) >= 0` and `u*x + v*y = d`. */
function bezout(x: bigint, y: bigint): [bigint, bigint, bigint] {
  let old_r = x,
    r = y;
  let old_s = 1n,
    s = 0n;
  let old_t = 0n,
    t = 1n;
  while (r !== 0n) {
    const q = fdiv(old_r, r); // any consistent quotient works
    [old_r, r] = [r, old_r - q * r];
    [old_s, s] = [s, old_s - q * s];
    [old_t, t] = [t, old_t - q * t];
  }
  if (old_r < 0n) {
    old_r = -old_r;
    old_s = -old_s;
    old_t = -old_t;
  }
  return [old_r, old_s, old_t];
}

/** floor division */
function fdiv(a: bigint, b: bigint): bigint {
  let q = a / b;
  if (a % b !== 0n && (a < 0n) !== (b < 0n)) q -= 1n;
  return q;
}

/**
 * PARI `truedvmdii(x,y,&r)`: `x = q*y + r` with `0 <= r < |y|`.
 * `bb_hnf.c:81` `_Fp_lquo`.
 */
function truedvmdii(x: bigint, y: bigint): [bigint, bigint] {
  if (y === 0n) throw new Error('truedvmdii: division by zero');
  let q = x / y;
  let r = x - q * y;
  if (r < 0n) {
    if (y > 0n) {
      q -= 1n;
      r += y;
    } else {
      q += 1n;
      r -= y;
    }
  }
  return [q, r];
}

function Fp_red_(x: bigint, N: bigint): bigint {
  const r = x % N;
  return r < 0n ? r + N : r;
}

/** number of bits of |x| (PARI `expi(x)+1` for x != 0) */
function bitlen(x: bigint): number {
  x = iabs(x);
  return x === 0n ? 0 : x.toString(2).length;
}

/** PARI `Z_split` (`bb_hnf.c:85-96`): D = M*N with p|M => p∤a, p|N => p|a; return M. */
function Z_split(D: bigint, a: bigint): bigint {
  const e = bitlen(D) - 1; // expi(D)
  const n = e < 2 ? 1 : bitlen(BigInt(e)) - 1 + 1; // expu(e)+1
  for (let i = 1; i <= n; i++) a = Fp_red_(a * a, D);
  const N = gcdii(a, D);
  return D / N;
}

/** PARI `Z_stab` (`bb_hnf.c:99-108`): c s.t. gcd(a+c*b,N) = gcd(a,b,N), no factoring. */
function Z_stab(a: bigint, b: bigint, N: bigint): bigint {
  let g = gcdii(a, b);
  g = gcdii(g, N);
  const N2 = N / g;
  const a2 = a / g;
  return Z_split(N2, a2);
}

/* ------------------------------------------------------------------ */
/* The Fp black-box Hermite ring (`bb_hnf.c:66-171`)                   */
/* ------------------------------------------------------------------ */

class FpHermite {
  constructor(readonly N: bigint) {}

  add(x: bigint, y: bigint): bigint {
    return x + y;
  }
  neg(x: bigint): bigint {
    return -x;
  }
  mul(x: bigint, y: bigint): bigint {
    return x * y;
  }
  red(x: bigint): bigint {
    return Fp_red_(x, this.N);
  }
  equal0(x: bigint): boolean {
    return x === 0n;
  }
  equal1(x: bigint): boolean {
    return x === 1n;
  }
  s(x: number): bigint {
    if (x === 0) return 0n;
    if (x === 1) return 1n;
    return Fp_red_(BigInt(x), this.N);
  }
  /** `_Fp_rann`: b with b*R = ann(a) */
  rann(x: bigint): bigint {
    if (x === 0n) return 1n;
    const d = gcdii(x, this.N);
    return Fp_red_(this.N / d, this.N);
  }
  /** `_Fp_lquo`: q with r = x - y*q the canonical representative of x mod yR */
  lquo(x: bigint, y: bigint): [bigint, bigint] {
    return truedvmdii(x, y);
  }
  /** `_Fp_unit`: returns `[g,u]` with u a unit and x*u = g the canonical
   *  generator of xR, or `null` if x is already canonical (or zero). */
  unit(x: bigint): [bigint, bigint] | null {
    const N = this.N;
    if (x === 0n) return null;
    const [g, s0] = bezout(x, N);
    let s = s0;
    if (g === 1n || gcdii(s, N) === 1n) return [g, s];
    const N2 = N / g;
    for (let i = 0; i < 5; i++) {
      s = s + N2;
      if (gcdii(s, N) === 1n) return [g, s];
    }
    let d = Z_stab(s, N2, N);
    d = d * N2;
    const v = Fp_red_(s + d, N);
    if (v === 1n) return null;
    return [g, v];
  }
  /** `_Fp_extgcd`: `[d, U]` with U in GL_2 and `[x;y]*U = [0;d]`.
   *  U is returned row-major as `[[U11,U12],[U21,U22]]`.
   *  `smallop` is true when U needs no reduction afterwards. */
  extgcd(x: bigint, y: bigint): { d: bigint; U: bigint[][]; smallop: boolean } {
    if (y === 1n) {
      return {
        d: y,
        U: [
          [1n, 0n],
          [Fp_red_(-x, this.N), 1n],
        ],
        smallop: true,
      };
    }
    const [d, u, v] = bezout(x, y);
    if (d === 0n) {
      return {
        d,
        U: [
          [1n, 0n],
          [0n, 1n],
        ],
        smallop: false,
      };
    }
    return {
      d,
      U: [
        [y / d, u],
        [-(x / d), v],
      ],
      smallop: false,
    };
  }
}

/* ------------------------------------------------------------------ */
/* Column primitives (`bb_hnf.c:173-329`)                              */
/* ------------------------------------------------------------------ */

function gen_redcol(C: Col, lim: number, R: FpHermite): void {
  for (let i = 1; i <= lim; i++) if (!R.equal0(C[i]!)) C[i] = R.red(C[i]!);
}

/** `gen_rightmulcol`: C*a, or null if a == 0. */
function gen_rightmulcol(
  C: Col,
  a: bigint,
  lim: number,
  fillzeros: boolean,
  R: FpHermite
): Col | null {
  if (R.equal1(a)) return C;
  if (R.equal0(a)) return null;
  // NB: PARI leaves the slots past `lim` uninitialised unless `fillzeros`;
  // they are never read there, so zero-filling is behaviourally identical
  // and avoids `undefined` slots in JS.
  const Ca: Col = new Array(C.length).fill(0n);
  for (let i = 1; i <= lim; i++) Ca[i] = R.mul(C[i]!, a);
  void fillzeros;
  return Ca;
}

/** C1 <- C1 + C2 (C2 assumed 0 after lim) */
function gen_addcol(C1: Col, C2: Col, lim: number, R: FpHermite): void {
  for (let i = 1; i <= lim; i++) C1[i] = R.add(C1[i]!, C2[i]!);
}

/** H[,i] <- H[,i] + C*a */
function gen_addrightmul(
  H: Mat,
  C: Col,
  a: bigint,
  i: number,
  lim: number,
  R: FpHermite
): void {
  if (R.equal0(a)) return;
  const Ca = gen_rightmulcol(C, a, lim, false, R);
  if (Ca) gen_addcol(H[i]!, Ca, lim, R);
}

function gen_zerocol(n: number): Col {
  const C: Col = new Array(n + 1).fill(0n);
  return C;
}

function gen_zeromat(m: number, n: number): Mat {
  const M: Mat = [[]];
  for (let i = 1; i <= n; i++) M.push(gen_zerocol(m));
  return M;
}

function gen_colei(n: number, i: number): Col {
  const C = gen_zerocol(n);
  C[i] = 1n;
  return C;
}

function gen_matid_hermite(n: number): Mat {
  const M: Mat = [[]];
  for (let i = 1; i <= n; i++) M.push(gen_colei(n, i));
  return M;
}

function matDims(M: Mat): [number, number] {
  const n = M.length - 1;
  const m = n === 0 ? 0 : M[1]!.length - 1;
  return [m, n];
}

function gen_matmul_hermite(A: Mat, B: Mat, R: FpHermite): Mat {
  const [a, c] = matDims(A);
  const [c2, b] = matDims(B);
  if (c !== c2) throw new PariDimError('gen_matmul_hermite');
  const M: Mat = [[]];
  for (let j = 1; j <= b; j++) {
    const col: Col = new Array(a + 1);
    col[0] = DUMMY;
    for (let i = 1; i <= a; i++) {
      let sum = 0n;
      for (let k = 1; k <= c; k++) sum = R.add(sum, R.mul(A[k]![i]!, B[j]![k]!));
      col[i] = sum;
    }
    gen_redcol(col, a, R);
    M.push(col);
  }
  return M;
}

/** C <- A*u1 + B*u2 (all zero after lim) */
function gen_rightlincomb(
  A: Col,
  B: Col,
  u1: bigint,
  u2: bigint,
  lim: number,
  R: FpHermite
): Col {
  const Au1 = gen_rightmulcol(A, u1, lim, true, R);
  const Bu2 = gen_rightmulcol(B, u2, lim, true, R);
  if (!Au1 && !Bu2) return gen_zerocol(A.length - 1);
  if (!Au1) return Bu2!;
  if (!Bu2) return Au1;
  const out = Au1 === A ? A.slice() : Au1;
  gen_addcol(out, Bu2, lim, R);
  return out;
}

/** (H[,i] | H[,j]) <- (H[,i] | H[,j]) * U */
function gen_elem(
  H: Mat,
  U: bigint[][],
  i: number,
  j: number,
  lim: number,
  R: FpHermite
): void {
  const Hi = H[i]!.slice();
  const Hj = H[j]!.slice();
  H[i] = gen_rightlincomb(Hi, Hj, U[0]![0]!, U[1]![0]!, lim, R);
  H[j] = gen_rightlincomb(Hi, Hj, U[0]![1]!, U[1]![1]!, lim, R);
}

function gen_is_zerocol(C: Col, lim: number, R: FpHermite): boolean {
  for (let i = 1; i <= lim; i++) if (!R.equal0(C[i]!)) return false;
  return true;
}

/* ------------------------------------------------------------------ */
/* Elementary operations (`bb_hnf.c:331-460`)                          */
/* ------------------------------------------------------------------ */

type Op =
  | { k: 'perm'; perm: number[] } // 1-indexed permutation, vecpermute semantics
  | { k: 'swap'; i: number; j: number }
  | { k: 'mul'; i: number; u: bigint }
  | { k: 'transv'; i: number; j: number; a: bigint } // Ci <- Ci + Cj*a
  | { k: 'U'; i: number; j: number; U: bigint[][] };

function mkoptransv(i: number, j: number, a: bigint, R: FpHermite): Op | null {
  a = R.red(a);
  if (R.equal0(a)) return null;
  return { k: 'transv', i, j, a };
}

function mkopU(i: number, j: number, U: bigint[][], R: FpHermite): Op | null {
  if (R.equal1(U[0]![0]!) && R.equal0(U[0]![1]!) && R.equal1(U[1]![1]!))
    return mkoptransv(i, j, U[1]![0]!, R);
  return { k: 'U', i, j, U };
}

function mkopmul(i: number, u: bigint, R: FpHermite): Op | null {
  if (R.equal1(u)) return null;
  return { k: 'mul', i, u };
}

/** apply op to M by right multiplication */
function gen_rightapply(M: Mat, op: Op, R: FpHermite): void {
  const [m] = matDims(M);
  switch (op.k) {
    case 'perm': {
      const M2: Mat = [[]];
      for (let i = 1; i < op.perm.length; i++) M2.push(M[op.perm[i]!]!);
      for (let i = 1; i < M.length; i++) M[i] = M2[i]!;
      return;
    }
    case 'swap': {
      const t = M[op.i]!;
      M[op.i] = M[op.j]!;
      M[op.j] = t;
      return;
    }
    case 'mul': {
      const c = gen_rightmulcol(M[op.i]!, op.u, m, false, R);
      if (c) M[op.i] = c;
      gen_redcol(M[op.i]!, m, R);
      return;
    }
    case 'transv': {
      gen_addrightmul(M, M[op.j]!, op.a, op.i, m, R);
      gen_redcol(M[op.i]!, m, R);
      return;
    }
    case 'U': {
      gen_elem(M, op.U, op.i, op.j, m, R);
      gen_redcol(M[op.i]!, m, R);
      gen_redcol(M[op.j]!, m, R);
      return;
    }
  }
}

/** apply op to the column C by left multiplication */
function gen_leftapply(C: Col, op: Op, R: FpHermite): void {
  switch (op.k) {
    case 'perm': {
      // C2 = vecpermute(C, perm_inv(op)); C <- C2
      const inv: number[] = new Array(op.perm.length);
      for (let i = 1; i < op.perm.length; i++) inv[op.perm[i]!] = i;
      const C2: Col = new Array(C.length);
      C2[0] = DUMMY;
      for (let i = 1; i < C.length; i++) C2[i] = C[inv[i]!]!;
      for (let i = 1; i < C.length; i++) C[i] = C2[i]!;
      return;
    }
    case 'swap': {
      const t = C[op.i]!;
      C[op.i] = C[op.j]!;
      C[op.j] = t;
      return;
    }
    case 'mul': {
      C[op.i] = R.red(R.mul(op.u, C[op.i]!));
      return;
    }
    case 'transv': {
      // Ci <- Ci + Cj*a  ==>  transpose acts as Cj += a*Ci
      if (R.equal0(C[op.i]!)) return;
      C[op.j] = R.add(C[op.j]!, R.mul(op.a, C[op.i]!));
      return;
    }
    case 'U': {
      const X = op.U;
      const ci = C[op.i]!,
        cj = C[op.j]!;
      C[op.i] = R.red(R.add(R.mul(X[0]![0]!, ci), R.mul(X[0]![1]!, cj)));
      C[op.j] = R.red(R.add(R.mul(X[1]![0]!, ci), R.mul(X[1]![1]!, cj)));
      return;
    }
  }
}

function permSign(perm: number[]): number {
  const n = perm.length - 1;
  const seen = new Array(n + 1).fill(false);
  let sign = 1;
  for (let i = 1; i <= n; i++) {
    if (seen[i]) continue;
    let len = 0;
    let j = i;
    while (!seen[j]) {
      seen[j] = true;
      j = perm[j]!;
      len++;
    }
    if (len % 2 === 0) sign = -sign;
  }
  return sign;
}

/** `gen_detops` (`bb_hnf.c:462-509`) */
function gen_detops(ops: Op[], R: FpHermite): bigint {
  let d = 1n;
  for (const op of ops) {
    switch (op.k) {
      case 'perm':
        if (permSign(op.perm) < 0) d = R.neg(d);
        break;
      case 'swap':
        d = R.neg(d);
        break;
      case 'mul':
        d = R.red(R.mul(d, op.u));
        break;
      case 'U': {
        const A = op.U[0]![0]!,
          B = op.U[0]![1]!,
          C = op.U[1]![0]!,
          D = op.U[1]![1]!;
        d = R.red(R.mul(d, R.add(R.mul(A, D), R.neg(R.mul(B, C)))));
        break;
      }
      case 'transv':
        break;
    }
  }
  return d;
}

function gen_is_inv(x: bigint, R: FpHermite): boolean {
  const u = R.unit(x);
  if (!u) return R.equal1(x);
  return R.equal1(u[0]);
}

function gen_last_inv_diago(A: Mat, R: FpHermite): number {
  const [m, n] = matDims(A);
  for (let i = 1, j = n - m + 1; i <= m; i++, j++)
    if (!gen_is_inv(A[j]![i]!, R)) return i - 1;
  return m;
}

/* ------------------------------------------------------------------ */
/* Howell form (`bb_hnf.c:511-772` gen_howell_i)                        */
/* ------------------------------------------------------------------ */

interface HowellResult {
  H: Mat | null;
  ops?: Op[];
  /** set when early_abort triggered: the offending non-invertible pivot */
  abortPivot?: bigint;
}

/**
 * `gen_howell_i` (`bb_hnf.c:513-772`).
 * @param removeZerocols 0 none, 1 until square, 2 all
 */
function gen_howell_i(
  A: Mat,
  removeZerocols: 0 | 1 | 2,
  permuteZerocols: boolean,
  earlyAbort: boolean,
  onlyTriangular: boolean,
  wantOps: boolean,
  R: FpHermite
): HowellResult {
  let [m, n] = matDims(A);
  const ops: Op[] = [];
  let piv = 0n;
  let lastinv = 0;
  const zero = 0n;
  const one = 1n;

  if (earlyAbort && n < m) return { H: null, abortPivot: zero };

  let H: Mat;
  let extra: number;
  if (n < m + 1) {
    extra = m + 1 - n;
    H = gen_zeromat(m, extra);
    for (let j = 1; j <= n; j++) H.push(A[j]!.slice());
  } else {
    extra = 0;
    H = [[]];
    for (let j = 1; j <= n; j++) H.push(A[j]!.slice());
  }
  [m, n] = matDims(H);
  const s = n - m; // shift

  /* put in triangular form */
  for (let i = m, si = s + m; i > 0 && si > extra; i--, si--) {
    H[si]![i] = R.red(H[si]![i]!);
    for (let j = extra + 1; j < si; j++) {
      H[j]![i] = R.red(H[j]![i]!);
      if (R.equal0(H[j]![i]!)) continue;
      const eg = R.extgcd(H[j]![i]!, H[si]![i]!);
      let d = eg.d;
      const U = eg.U;
      if (n > 10) {
        const u = R.unit(d);
        if (u) {
          U[0]![1] = R.mul(U[0]![1]!, u[1]);
          U[1]![1] = R.mul(U[1]![1]!, u[1]);
          d = u[0];
        }
      }
      gen_elem(H, U, j, si, i - 1, R);
      if (wantOps) {
        const op = mkopU(j, si, U, R);
        if (op) ops.push(op);
      }
      H[j]![i] = zero;
      H[si]![i] = d;
      if (!eg.smallop) {
        gen_redcol(H[si]!, i - 1, R);
        gen_redcol(H[j]!, i - 1, R);
      }
    }

    if (earlyAbort) {
      let d = H[si]![i]!;
      const u = R.unit(d);
      if (u) d = u[0];
      if (!R.equal1(d)) return { H: null, abortPivot: d };
    }
  }

  if (!wantOps) lastinv = gen_last_inv_diago(H, R);

  /* put in reduced Howell form */
  if (!onlyTriangular) {
    for (let i = m, si = s + m; i > 0; i--, si--) {
      /* normalize diagonal coefficient */
      if (i <= lastinv) {
        H[si]![i] = one;
      } else {
        const u = R.unit(H[si]![i]!);
        if (u) {
          const c = gen_rightmulcol(H[si]!, u[1], i - 1, true, R);
          if (c) H[si] = c;
          H[si]![i] = u[0];
          gen_redcol(H[si]!, i - 1, R);
          if (wantOps) {
            const op = mkopmul(si, u[1], R);
            if (op) ops.push(op);
          }
        } else {
          H[si]![i] = R.red(H[si]![i]!);
        }
      }
      piv = H[si]![i]!;

      /* reduce above diagonal */
      if (!R.equal0(piv)) {
        const C = H[si]!;
        for (let j = si + 1; j <= n; j++) {
          if (i <= lastinv) {
            H[j]![i] = zero;
          } else {
            H[j]![i] = R.red(H[j]![i]!);
            let q: bigint, r: bigint;
            if (R.equal1(piv)) {
              q = H[j]![i]!;
              r = zero;
            } else {
              [q, r] = R.lquo(H[j]![i]!, piv);
            }
            q = R.neg(q);
            gen_addrightmul(H, C, q, j, i - 1, R);
            if (wantOps) {
              const op = mkoptransv(j, si, q, R);
              if (op) ops.push(op);
            }
            H[j]![i] = r;
          }
        }
      }

      /* ensure Howell property */
      if (i > 1) {
        const a = R.rann(piv);
        if (!R.equal0(a)) {
          const c = gen_rightmulcol(H[si]!, a, i - 1, true, R);
          H[1] = c === null ? gen_zerocol(m) : c === H[si] ? c.slice() : c;
          if (wantOps) {
            const op = mkoptransv(1, si, a, R);
            if (op) ops.push(op);
          }
          for (let i2 = i - 1, si2 = s + i2; i2 > 0; i2--, si2--) {
            H[1]![i2] = R.red(H[1]![i2]!);
            if (R.equal0(H[1]![i2]!)) continue;
            H[si2]![i2] = R.red(H[si2]![i2]!);
            if (R.equal0(H[si2]![i2]!)) {
              const t = H[1]!;
              H[1] = H[si2]!;
              H[si2] = t;
              if (wantOps) ops.push({ k: 'swap', i: 1, j: si2 });
              continue;
            }
            const eg = R.extgcd(H[1]![i2]!, H[si2]![i2]!);
            const d = eg.d;
            const U = eg.U;
            gen_elem(H, U, 1, si2, i2 - 1, R);
            if (wantOps) {
              const op = mkopU(1, si2, U, R);
              if (op) ops.push(op);
            }
            H[1]![i2] = zero;
            H[si2]![i2] = d;
            if (!eg.smallop) {
              gen_redcol(H[si2]!, i2, R);
              gen_redcol(H[1]!, i2 - 1, R);
            }
          }
        }
      }
    }
  }

  for (let j = 1; j <= n; j++) {
    const lim = Math.max(0, m - n + j);
    gen_redcol(H[j]!, lim, R);
    for (let i = lim + 1; i <= m; i++) H[j]![i] = zero;
  }

  /* put zero columns first */
  const iszero: boolean[] = new Array(n + 1).fill(false);
  let nbz = 0;
  for (let i = 1; i <= n; i++) {
    iszero[i] = gen_is_zerocol(H[i]!, Math.max(0, m - n + i), R);
    if (iszero[i]) nbz++;
  }

  let j = 1;
  let perm: number[];
  if (permuteZerocols) {
    perm = new Array(n + 1).fill(0);
    for (let i = 1; i <= n; i++)
      if (iszero[i]) {
        perm[j] = i;
        j++;
      }
  } else {
    perm = new Array(n - nbz + 1).fill(0);
  }
  for (let i = 1; i <= n; i++)
    if (!iszero[i]) {
      perm[j] = i;
      j++;
    }

  if (permuteZerocols || removeZerocols === 2) {
    const H2: Mat = [[]];
    for (let i = 1; i < perm.length; i++) H2.push(H[perm[i]!]!);
    H = H2;
  }
  if (permuteZerocols && removeZerocols === 2) H = vecslice(H, nbz + 1, n);
  if (removeZerocols === 1) H = vecslice(H, s + 1, n);
  if (permuteZerocols && wantOps) ops.push({ k: 'perm', perm });

  return wantOps ? { H, ops } : { H };
}

function vecslice(M: Mat, a: number, b: number): Mat {
  const out: Mat = [[]];
  for (let i = a; i <= b; i++) out.push(M[i]!);
  return out;
}

function matslice(M: Mat, r1: number, r2: number, c1: number, c2: number): Mat {
  const out: Mat = [[]];
  for (let j = c1; j <= c2; j++) {
    const col: Col = new Array(r2 - r1 + 2);
    col[0] = DUMMY;
    for (let i = r1; i <= r2; i++) col[i - r1 + 1] = M[j]![i]!;
    out.push(col);
  }
  return out;
}

/* ------------------------------------------------------------------ */
/* Kernel (`bb_hnf.c:826-885`)                                          */
/* ------------------------------------------------------------------ */

/** `gen_kernel_howell` (`bb_hnf.c:827-848`); H in true Howell form, no zero cols. */
function gen_kernel_howell(H: Mat, R: FpHermite): Mat {
  const [m, n] = matDims(H);
  const K = gen_zeromat(n, n);
  let i = m;
  for (let j = n; j > 0; j--) {
    while (R.equal0(H[j]![i]!)) i--;
    const piv = H[j]![i]!;
    if (R.equal0(piv)) continue;
    K[j]![j] = R.rann(piv);
    if (j < n) {
      const FK = gen_matmul_hermite(
        matslice(H, i, i, j + 1, n),
        matslice(K, j + 1, n, j + 1, n),
        R
      );
      for (let j2 = j + 1; j2 <= n; j2++)
        K[j2]![j] = R.neg(R.lquo(FK[j2 - j]![1]!, piv)[0]);
    }
  }
  return K;
}

/** `gen_kernel_from_howell` (`bb_hnf.c:850-877`) */
function gen_kernel_from_howell(H: Mat, ops: Op[], n: number, R: FpHermite): Mat {
  const [m, r] = matDims(H);
  if (!r) return gen_matid_hermite(n);
  const n2 = Math.max(n, m + 1);
  const extra = n2 - n;
  const nbz = n2 - r;
  const KH = gen_kernel_howell(H, R);
  const K: Mat = [[]];
  for (let i = 1; i <= nbz; i++) K.push(gen_colei(nbz + r, i));
  for (let i = 1; i <= r; i++) {
    const col: Col = new Array(nbz + r + 1);
    col[0] = DUMMY;
    for (let t = 1; t <= nbz; t++) col[t] = 0n;
    for (let t = 1; t <= r; t++) col[nbz + t] = KH[i]![t]!;
    K.push(col);
  }
  for (let i = 1; i < K.length; i++) {
    for (let o = ops.length - 1; o >= 0; o--) gen_leftapply(K[i]!, ops[o]!, R);
    gen_redcol(K[i]!, nbz + r, R);
  }
  /* rowpermute(K, cyclic_perm(n2, extra)) */
  const cperm: number[] = new Array(n2 + 1).fill(0);
  {
    let i = 1;
    for (; i <= n2 - extra; i++) cperm[i] = i + extra;
    for (; i <= n2; i++) cperm[i] = i - n2 + extra;
  }
  const K2: Mat = [[]];
  for (let j = 1; j < K.length; j++) {
    const col: Col = new Array(n2 + 1);
    col[0] = DUMMY;
    for (let i = 1; i <= n2; i++) col[i] = K[j]![cperm[i]!]!;
    K2.push(col);
  }
  const KK = gen_howell_i(K2, 2, false, false, false, false, R).H!;
  {
    let i = n2;
    for (let j = KK.length - 1; j > 0; j--) {
      while (R.equal0(KK[j]![i]!)) i--;
      if (i <= n) return matslice(KK, 1, n, 1, j);
    }
  }
  return [[]];
}

/** `gen_kernel` (`bb_hnf.c:879-889`) */
function gen_kernel(A: Mat, wantIm: boolean, R: FpHermite): { K: Mat; im?: Mat } {
  const n = A.length - 1;
  const { H, ops } = gen_howell_i(A, 2, true, false, false, true, R);
  const K = gen_kernel_from_howell(H!, ops!, n, R);
  return wantIm ? { K, im: H! } : { K };
}

/** `gen_matimage` (`bb_hnf.c:806-824`) */
function gen_matimage(A: Mat, wantU: boolean, R: FpHermite): { H: Mat; U?: Mat } {
  if (wantU) {
    const [m, n] = matDims(A);
    const { H, ops } = gen_howell_i(A, 2, true, false, false, true, R);
    const r = H!.length - 1;
    let U: Mat = gen_zeromat(n, Math.max(0, m - n + 1));
    const id = gen_matid_hermite(n);
    for (let i = 1; i <= n; i++) U.push(id[i]!);
    const n2 = U.length - 1;
    for (const op of ops!) gen_rightapply(U, op, R);
    if (r < n2) U = vecslice(U, n2 - r + 1, n2);
    return { H: H!, U };
  }
  return { H: gen_howell_i(A, 2, false, false, false, false, R).H! };
}

/** `gen_inv` (`bb_hnf.c:892-917`): left inverse of the transpose input. */
function gen_inv(A: Mat, R: FpHermite): Mat {
  const [m, n] = matDims(A);
  const res = gen_howell_i(A, 0, false, true, false, true, R);
  if (!res.H) throw new PariInvError('gen_inv');
  const H = res.H;
  const ops = res.ops!;
  const n2 = H.length - 1;
  let U = gen_zeromat(n2, m);
  for (let j = 1; j <= m; j++) U[j]![j + n2 - m] = 1n;
  for (let j = 1; j <= m; j++) {
    for (let o = ops.length - 1; o >= 0; o--) gen_leftapply(U[j]!, ops[o]!, R);
    gen_redcol(U[j]!, n2, R);
  }
  if (n2 > n) {
    const U2: Mat = [[]];
    for (let j = 1; j <= m; j++) {
      const col: Col = new Array(n + 1);
      col[0] = DUMMY;
      for (let i = 1; i <= n; i++) col[i] = U[j]![n2 - n + i]!;
      U2.push(col);
    }
    U = U2;
  }
  return U;
}

function transposeMat(A: Mat): Mat {
  const [m, n] = matDims(A);
  const out: Mat = [[]];
  for (let i = 1; i <= m; i++) {
    const col: Col = new Array(n + 1);
    col[0] = DUMMY;
    for (let j = 1; j <= n; j++) col[j] = A[j]![i]!;
    out.push(col);
  }
  return out;
}

/* ------------------------------------------------------------------ */
/* Public API                                                          */
/* ------------------------------------------------------------------ */

function checkd(fun: string, d: bigint): void {
  if (typeof d !== 'bigint') throw new PariTypeError(fun, 'not a t_INT');
  if (d <= 0n) throw new PariDomainError(fun, 'd', '<=', '0');
}

/**
 * PARI `matkermod(A, d, &im)` (`bb_hnf.c:1036-1050`).
 *
 * Returns a Howell basis of the kernel of `A` modulo `d` (an arbitrary positive
 * integer, **not necessarily prime**), as the columns of the returned matrix.
 *
 * @param A matrix in column-major {@link ZM} layout
 * @param d modulus, `d > 0`
 * @param wantIm if true, also return a basis of the image of `A`
 */
export function matkermod(
  A: ZM,
  d: bigint,
  wantIm = false
): { ker: ZM; im?: ZM } {
  checkd('matkermod', d);
  if (d === 1n) return wantIm ? { ker: [], im: [] } : { ker: [] };
  const R = new FpHermite(d);
  let { H, m, n } = toInternal(A);
  if (!wantIm && m > 2 * n) {
    // A = shallowtrans(matimagemod(shallowtrans(A), d, NULL))
    const T = transposeMat(H);
    const im = gen_matimage(T, false, R).H;
    H = transposeMat(im);
    [m, n] = matDims(H);
  }
  const res = gen_kernel(H, wantIm, R);
  const ker = toExternal(res.K, n);
  if (wantIm) return { ker, im: toExternal(res.im!, m) };
  return { ker };
}

/**
 * Convenience wrapper for callers using the usual row-major layout (e.g. Sage's
 * `right_kernel_matrix`): takes `A` as an array of **rows** and returns the
 * kernel basis as an array of **rows** (each returned row `v` satisfies
 * `A*v == 0 mod d`).
 */
export function matkermod_basis(
  A: readonly (readonly bigint[])[],
  d: bigint
): bigint[][] {
  const { ker } = matkermod(zm_from_rows(A), d);
  // ker is column-major: each ker[j] is a basis vector already.
  return ker.map((c) => c.slice());
}

/**
 * PARI `matimagemod(A, d, &U)` (`bb_hnf.c:981-990`).
 * Howell basis of the image of `A` mod `d`, as columns. If `wantU`, also returns
 * `U` with `A*U = H` (mod d).
 */
export function matimagemod(A: ZM, d: bigint, wantU = false): { im: ZM; U?: ZM } {
  checkd('matimagemod', d);
  if (d === 1n) return wantU ? { im: [], U: [] } : { im: [] };
  const R = new FpHermite(d);
  const { H, m } = toInternal(A);
  const res = gen_matimage(H, wantU, R);
  const im = toExternal(res.H, m);
  if (wantU) return { im, U: toExternal(res.U!, A.length) };
  return { im };
}

/**
 * PARI `matdetmod(A, d)` (`bb_hnf.c:1015-1033`).
 * Determinant of the square matrix `A` modulo `d`.
 */
export function matdetmod(A: ZM, d: bigint): bigint {
  checkd('matdetmod', d);
  const n = A.length;
  if (!n) return d === 1n ? 0n : 1n;
  const { H: A1, m } = toInternal(A);
  if (n !== m) throw new PariDimError('matdetmod');
  if (d === 1n) return 0n;
  const R = new FpHermite(d);
  const { H, ops } = gen_howell_i(A1, 1, false, false, true, true, R);
  let D = gen_detops(ops!, R);
  D = Fp_inv_(D, d);
  for (let i = 1; i <= n; i++) D = Fp_red_(D * H![i]![i]!, d);
  return D;
}

function Fp_inv_(a: bigint, p: bigint): bigint {
  const [g, u] = bezout(Fp_red_(a, p), p);
  if (g !== 1n) throw new PariInvError(`Fp_inv (Mod(${Fp_red_(a, p)}, ${p}))`);
  return Fp_red_(u, p);
}

/**
 * PARI `matinvmod(A, d)` (`bb_hnf.c:1053-1071`): a left inverse of `A` mod `d`.
 */
export function matinvmod(A: ZM, d: bigint): ZM {
  checkd('matinvmod', d);
  const { H, m, n } = toInternal(A);
  if (d === 1n) {
    if (m < n) throw new PariInvError('matinvmod');
    return zm_from_rows(
      Array.from({ length: n }, () => Array.from({ length: m }, () => 0n))
    );
  }
  const R = new FpHermite(d);
  const U = gen_inv(transposeMat(H), R);
  const Ut = transposeMat(U);
  return toExternal(Ut, matDims(Ut)[0]);
}
