/**
 * @module parigp-ts/buch
 * @description Class group and unit group of a quadratic field
 * (McCurley / Buchmann index calculus), a direct port of
 * `reference/pari/src/basemath/buch1.c` together with the pieces of
 * `reference/pari/src/basemath/buch2.c`, `hnf_snf.c`, `Qfb.c` and `alglin1.c`
 * that it calls.
 *
 * Ported functions, with their upstream location:
 *
 * - `buch1.c:1035-1241`  `Buchquad_i`, `buch1.c:1242-1248` `Buchquad`
 * - `buch1.c:1260-1295`  `quadclassunit0`, `quadclassno`, `quadclassnos`
 * - `buch1.c:189-197`    `bnf_increase_LIMC`
 * - `buch1.c:199-262`    `isless_iu`, `Z_isquasismooth_prod`, `factorquad`
 * - `buch1.c:264-302`    `largeprime` (large prime relations), `clearhash`
 * - `buch1.c:318-338`    `cache_prime_quad`
 * - `buch1.c:340-363`    `compute_invresquad`
 * - `buch1.c:365-390`    `is_bad`, `nthidealquad`
 * - `buch1.c:392-427`    `quadGRHchk`
 * - `buch1.c:430-505`    `FBquad`, `subFBquad`
 * - `buch1.c:507-536`    `powsubFBquad`
 * - `buch1.c:538-592`    `sub_fact`, `add_fact`
 * - `buch1.c:594-645`    `get_clgp`, `trivial_relations`
 * - `buch1.c:656-775`    `rel_to_col`, `imag_relations`, `imag_be_honest`
 * - `buch1.c:777-970`    `dist`, `real_relations`, `real_be_honest`
 * - `buch1.c:972-1032`   `crabs`, `gcdreal`, `get_R`, `quad_be_honest`
 * - `buch1.c:58-173`     `qfr3_canon`, `qfr5_canon`, `QFR3_comp`, `QFR5_comp`,
 *                        `qfr5_rho_pow`, `qfr5_pf`, `qfr3_pf`, `init_form`,
 *                        `random_form`
 * - `buch2.c:354-386`    `init_GRHcheck`, `GRHok`
 * - `Qfb.c:23-31`        `check_quaddisc`
 * - `Qfb.c:409-560`      `fix_expo`, `qfr5_dist`, `rho_get_BC`, `qfr3_rho`,
 *                        `qfr5_rho`, `qfr_to_qfr5`, `ab_isreduced`, `qfr5_red`,
 *                        `qfr3_red`, `qfr_data_init`, `qfr_1_fill`, `qfr5_1`
 * - `Qfb.c:1013-1071`    `qfb_sqr`, `qfb_comp` (on 3-component containers)
 * - `Qfb.c:1478-1505`    `qfr5_compraw`, `qfr5_comp`, `qfr3_compraw`, `qfr3_comp`
 * - `Qfb.c:1543-1580`    `qfr3_powraw`, `qfr3_pow`
 * - `Qfb.c:1661-1694`    `primeform_u`
 * - `hnf_snf.c:70-201`   `count`, `count2`, `hnffinal`
 * - `hnf_snf.c:222-239`  `ZM_rowrankprofile`
 * - `hnf_snf.c:254-554`  `hnfspec_i`
 * - `hnf_snf.c:610-658`  `hnfadd_i`
 * - `hnf_snf.c:671-746`  `ZC_elem`
 * - `hnf_snf.c:1546-1810` `Minus`, `findi`, `findi_normalize`, `reduce2`,
 *                        `hnfswap`, `reverse_rows`, `must_swap`, `ZM_hnflll`
 * - `hnf_snf.c:2330-2346` `bezout_step`
 * - `hnf_snf.c:2355-2380` `ZM_snf_no_divide`, `ZM_redpart`
 * - `hnf_snf.c:2386-2570` `ZM_snfall_i` (square HNF input, the only case
 *                        reachable from `get_clgp`)
 * - `hnf_snf.c:2820-2900` `ZM_snfclean`, `snf_group`, `ZM_snf_group`
 * - `alglin1.c:3695-3760` `ZM_pivots` (see the note on that function)
 * - `ZV.c:1526-1536`     `ZM_det_triangular`
 *
 * ## Conventions
 *
 * Matrices and vectors are **1-based** here: `V[0]` is a dummy slot, `V[i]` for
 * `i = 1..n` are the entries and `V.length = n+1` plays the role of PARI's
 * `lg()`.  Matrices are lists of columns, `M[j][i]` = entry `(i,j)`, exactly as
 * in PARI.  This keeps the transcription of `hnfspec`, whose index arithmetic
 * is intricate, literal.
 *
 * ## Real arithmetic
 *
 * `Buchquad` needs PARI's `t_REAL` (Shanks distances of indefinite forms, the
 * regulator, the residue of `zeta_K`).  `parigp-ts` has no multiprecision float
 * layer, so a minimal one is provided here ({@link Real}); it follows PARI's
 * representation (sign, normalized mantissa, binary exponent, precision in
 * bits) and the semantics of the `mp.c` primitives that `Buchquad` uses.  It is
 * a helper, not a transcription of `mp.c`: the elementary operations round to
 * nearest instead of reproducing PARI's exact rounding, so the last bits of a
 * regulator may differ from PARI's.  Accuracy is asserted by the algorithm
 * itself, exactly as upstream: `get_R` (`buch1.c:996`) accepts a tentative
 * regulator only when `h * R * invhr` lies in `(0.8, 1.3)`.
 */

import { Fp_sqrt, kronecker } from './ff.js';
import { NotImplementedError, Z_factor, isPrime } from './ifactor.js';
import { PariDomainError, PariSqrtnError, PariTypeError, matinvmod } from './matkermod.js';
import { type Qfb, mkqfb, primeform, qfbcomp, qfbpow, qfbred } from './qfb.js';

export { PariDomainError, PariSqrtnError, PariTypeError };
export { NotImplementedError };

/* ================================================================== */
/* Small integer helpers (PARI kernel)                                 */
/* ================================================================== */

const iabs = (x: bigint): bigint => (x < 0n ? -x : x);
const signe = (x: bigint): number => (x > 0n ? 1 : x < 0n ? -1 : 0);

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

/** extended gcd: returns `[d, u, v]` with `u*x + v*y = d >= 0`. PARI `bezout` */
function bezout(x: bigint, y: bigint): [bigint, bigint, bigint] {
  let a = x;
  let b = y;
  let u0 = 1n;
  let v0 = 0n;
  let u1 = 0n;
  let v1 = 1n;
  while (b !== 0n) {
    const q = a / b;
    const r = a - q * b;
    a = b;
    b = r;
    const nu = u0 - q * u1;
    const nv = v0 - q * v1;
    u0 = u1;
    v0 = v1;
    u1 = nu;
    v1 = nv;
  }
  if (a < 0n) {
    a = -a;
    u0 = -u0;
    v0 = -v0;
  }
  return [a, u0, v0];
}

/** PARI `truedvmdii`: Euclidean division with `0 <= r < |y|` */
function truedvmdii(x: bigint, y: bigint): [bigint, bigint] {
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

/** PARI `truedivii` */
const truedivii = (x: bigint, y: bigint): bigint => truedvmdii(x, y)[0];

/** PARI `diviiround`: round to nearest, ties away from zero */
function diviiround(x: bigint, y: bigint): bigint {
  const s = (x < 0n ? -1 : 1) * (y < 0n ? -1 : 1);
  const a = iabs(x);
  const b = iabs(y);
  const q = (2n * a + b) / (2n * b);
  return s < 0 ? -q : q;
}

function shifti(x: bigint, n: number): bigint {
  return n >= 0 ? x << BigInt(n) : x >> BigInt(-n);
}

/** number of bits of |x| minus 1 (PARI `expi`); `expi(0)` is undefined */
function expi(x: bigint): number {
  const a = iabs(x);
  if (a === 0n) return -1;
  return a.toString(2).length - 1;
}

/* PARI stores the sign apart, so `mod2`/`mod4`/... are residues of |x| */
const modk = (x: bigint, k: bigint): number => Number(iabs(x) % k);
const mod2 = (x: bigint): number => modk(x, 2n);
const mod4 = (x: bigint): number => modk(x, 4n);
const mod8 = (x: bigint): number => modk(x, 8n);
const mod16 = (x: bigint): number => modk(x, 16n);
const umodiu = (x: bigint, p: bigint): bigint => ((x % p) + p) % p;

/** integer square root */
export function sqrti(n: bigint): bigint {
  if (n < 0n) throw new PariDomainError('sqrti', 'n', '<', '0');
  if (n < 2n) return n;
  let x = 1n << BigInt((expi(n) >> 1) + 1);
  for (;;) {
    const y = (x + n / x) >> 1n;
    if (y >= x) break;
    x = y;
  }
  return x;
}

function Z_issquare(n: bigint): boolean {
  if (n < 0n) return false;
  const r = sqrti(n);
  return r * r === n;
}

/* ================================================================== */
/* 1-based vector / matrix helpers                                     */
/* ================================================================== */

/** 1-based column of integers (`V[0]` unused). */
export type ZC = bigint[];
/** 1-based matrix, list of columns: `M[j][i]` is the entry `(i,j)`. */
export type ZMat = ZC[];
/** 1-based vector of machine integers (PARI `t_VECSMALL`). */
export type zv = number[];

const lg = (v: readonly unknown[]): number => v.length;

function cgetg_int(n: number): ZC {
  const v: ZC = new Array<bigint>(n + 1).fill(0n);
  return v;
}
function cgetg_zv(n: number): zv {
  const v: zv = new Array<number>(n + 1).fill(0);
  return v;
}
function zero_zv(n: number): zv {
  return cgetg_zv(n);
}
function zerocol(n: number): ZC {
  return cgetg_int(n);
}
function identity_perm(n: number): zv {
  const v = cgetg_zv(n);
  for (let i = 1; i <= n; i++) v[i] = i;
  return v;
}
function matid(n: number): ZMat {
  const M: ZMat = new Array<ZC>(n + 1);
  M[0] = [];
  for (let j = 1; j <= n; j++) {
    const c = zerocol(n);
    c[j] = 1n;
    M[j] = c;
  }
  return M;
}
function zeromat(m: number, n: number): ZMat {
  const M: ZMat = new Array<ZC>(n + 1);
  M[0] = [];
  for (let j = 1; j <= n; j++) M[j] = zerocol(m);
  return M;
}
function nbrows(M: ZMat): number {
  return lg(M) === 1 ? 0 : lg(M[1]!) - 1;
}
function shallowcopyMat(M: ZMat): ZMat {
  const N: ZMat = new Array<ZC>(lg(M));
  N[0] = [];
  for (let j = 1; j < lg(M); j++) N[j] = M[j]!.slice();
  return N;
}
/** PARI `vecslice(v, a, b)` on a 1-based vector */
function vecslice<T>(v: T[], a: number, b: number): T[] {
  const out = new Array<T>(b - a + 2);
  out[0] = v[0]!;
  for (let i = a; i <= b; i++) out[i - a + 1] = v[i]!;
  return out;
}
function shallowconcat<T>(x: T[], y: T[]): T[] {
  const lx = lg(x);
  const ly = lg(y);
  const out = new Array<T>(lx + ly - 1);
  out[0] = x[0]!;
  for (let i = 1; i < lx; i++) out[i] = x[i]!;
  for (let i = 1; i < ly; i++) out[lx + i - 1] = y[i]!;
  return out;
}
function rowslice(M: ZMat, a: number, b: number): ZMat {
  const out: ZMat = new Array<ZC>(lg(M));
  out[0] = [];
  for (let j = 1; j < lg(M); j++) out[j] = vecslice(M[j]!, a, b);
  return out;
}
function rowpermute(M: ZMat, p: zv): ZMat {
  const out: ZMat = new Array<ZC>(lg(M));
  out[0] = [];
  for (let j = 1; j < lg(M); j++) {
    const c = M[j]!;
    const d = cgetg_int(lg(p) - 1);
    for (let i = 1; i < lg(p); i++) d[i] = c[p[i]!]!;
    out[j] = d;
  }
  return out;
}
function vecsmallpermute(v: zv, p: zv): zv {
  const out = cgetg_zv(lg(p) - 1);
  for (let i = 1; i < lg(p); i++) out[i] = v[p[i]!]!;
  return out;
}
function vconcat(A: ZMat, B: ZMat): ZMat {
  const out: ZMat = new Array<ZC>(lg(A));
  out[0] = [];
  for (let j = 1; j < lg(A); j++) {
    const a = A[j]!;
    const b = B[j]!;
    const c = cgetg_int(lg(a) - 1 + lg(b) - 1);
    for (let i = 1; i < lg(a); i++) c[i] = a[i]!;
    for (let i = 1; i < lg(b); i++) c[lg(a) - 1 + i] = b[i]!;
    out[j] = c;
  }
  return out;
}
/** PARI `ZM_mul` */
export function ZM_mul(A: ZMat, B: ZMat): ZMat {
  const n = lg(B) - 1;
  const m = lg(A) === 1 ? 0 : lg(A[1]!) - 1;
  const k = lg(A) - 1;
  const out: ZMat = new Array<ZC>(n + 1);
  out[0] = [];
  for (let j = 1; j <= n; j++) {
    const c = zerocol(m);
    const bj = B[j]!;
    for (let t = 1; t <= k; t++) {
      const v = bj[t]!;
      if (!v) continue;
      const at = A[t]!;
      for (let i = 1; i <= m; i++) c[i] += v * at[i]!;
    }
    out[j] = c;
  }
  return out;
}
function ZM_sub(A: ZMat, B: ZMat): ZMat {
  const out: ZMat = new Array<ZC>(lg(A));
  out[0] = [];
  for (let j = 1; j < lg(A); j++) {
    const a = A[j]!;
    const b = B[j]!;
    const c = cgetg_int(lg(a) - 1);
    for (let i = 1; i < lg(a); i++) c[i] = a[i]! - b[i]!;
    out[j] = c;
  }
  return out;
}
/** `A * B` where `B` is a small-integer matrix (PARI `ZM_zm_mul`) */
function ZM_zm_mul(A: ZMat, B: zv[]): ZMat {
  const n = lg(B) - 1;
  const m = lg(A) === 1 ? 0 : lg(A[1]!) - 1;
  const k = lg(A) - 1;
  const out: ZMat = new Array<ZC>(n + 1);
  out[0] = [];
  for (let j = 1; j <= n; j++) {
    const c = zerocol(m);
    const bj = B[j]!;
    for (let t = 1; t <= k; t++) {
      const v = bj[t]!;
      if (!v) continue;
      const bv = BigInt(v);
      const at = A[t]!;
      for (let i = 1; i <= m; i++) c[i] += bv * at[i]!;
    }
    out[j] = c;
  }
  return out;
}
function zm_to_ZM(M: zv[]): ZMat {
  const out: ZMat = new Array<ZC>(lg(M));
  out[0] = [];
  for (let j = 1; j < lg(M); j++) {
    const c = M[j]!;
    const d = cgetg_int(lg(c) - 1);
    for (let i = 1; i < lg(c); i++) d[i] = BigInt(c[i]!);
    out[j] = d;
  }
  return out;
}
/** rows `a..b` of `M` permuted by `perm` (PARI `rowslicepermute`) */
function rowslicepermute(M: zv[], perm: zv, a: number, b: number): zv[] {
  const out: zv[] = new Array<zv>(lg(M));
  out[0] = [];
  for (let j = 1; j < lg(M); j++) {
    const c = M[j]!;
    const d = cgetg_zv(b - a + 1);
    for (let i = a; i <= b; i++) d[i - a + 1] = c[perm[i]!]!;
    out[j] = d;
  }
  return out;
}
function ZV_togglesign(v: ZC): void {
  for (let i = 1; i < lg(v); i++) v[i] = -v[i]!;
}
/** PARI `ZC_lincomb1_inplace(X, Y, v)`: `X <- X + v Y` */
function ZC_lincomb1_inplace(X: ZC, Y: ZC, v: bigint): void {
  if (!v) return;
  for (let i = 1; i < lg(X); i++) X[i] = X[i]! + v * Y[i]!;
}
function ZV_equal0(v: ZC): boolean {
  for (let i = 1; i < lg(v); i++) if (v[i]) return false;
  return true;
}

/** PARI `ZM_det_triangular` (`ZV.c:1526`) */
export function ZM_det_triangular(M: ZMat): bigint {
  const l = lg(M);
  if (l < 3) return l < 2 ? 1n : M[1]![1]!;
  let s = M[1]![1]!;
  for (let i = 2; i < l; i++) s *= M[i]![i]!;
  return s;
}

/* ================================================================== */
/* Real arithmetic (helper: minimal PARI t_REAL)                       */
/* ================================================================== */

/**
 * A PARI `t_REAL`: value `s * m * 2^(e - p + 1)` with `2^(p-1) <= m < 2^p`
 * when `s != 0`.  `e` is PARI's `expo()`, `p` its bit precision.
 */
export interface Real {
  /** sign: -1, 0 or 1 */
  s: number;
  /** mantissa, `>= 0`, normalized to exactly `p` bits when `s != 0` */
  m: bigint;
  /** binary exponent: `|x| in [2^e, 2^(e+1))` */
  e: number;
  /** precision in bits */
  p: number;
}

/** PARI `DEFAULTPREC` (64-bit words): 64 bits of mantissa. */
export const DEFAULTPREC = 64;

export function real_0_bit(e: number, p = DEFAULTPREC): Real {
  return { s: 0, m: 0n, e, p };
}
export function real_0(p = DEFAULTPREC): Real {
  return real_0_bit(-p, p);
}
/** normalize a value `s * m * 2^k` to precision `p` bits (round to nearest) */
function mkreal(s: number, m: bigint, k: number, p: number): Real {
  if (s === 0 || m === 0n) return real_0(p);
  const bits = m.toString(2).length;
  if (bits > p) {
    const drop = bits - p;
    m = (m + (1n << BigInt(drop - 1))) >> BigInt(drop);
    k += drop;
    if (m.toString(2).length > p) {
      /* rounding overflowed to 2^p */
      m >>= 1n;
      k += 1;
    }
  } else if (bits < p) {
    m <<= BigInt(p - bits);
    k -= p - bits;
  }
  return { s: s > 0 ? 1 : -1, m, e: k + p - 1, p };
}
/** exact value of `x` as a scaled integer: `x = s * m * 2^shift(x)` */
const rshift = (x: Real): number => x.e - x.p + 1;

export function itor(x: bigint, p = DEFAULTPREC): Real {
  if (x === 0n) return real_0(p);
  return mkreal(x < 0n ? -1 : 1, iabs(x), 0, p);
}
export function real_1(p = DEFAULTPREC): Real {
  return { s: 1, m: 1n << BigInt(p - 1), e: 0, p };
}
export function real_neg(x: Real): Real {
  return { s: -x.s, m: x.m, e: x.e, p: x.p };
}
export function real_abs(x: Real): Real {
  return { s: x.s === 0 ? 0 : 1, m: x.m, e: x.e, p: x.p };
}
export function real_sign(x: Real): number {
  return x.s;
}
export function real_expo(x: Real): number {
  return x.e;
}
/** PARI `shiftr(x, n)`: exact multiplication by `2^n` */
export function shiftr(x: Real, n: number): Real {
  if (x.s === 0) return { s: 0, m: 0n, e: x.e + n, p: x.p };
  return { s: x.s, m: x.m, e: x.e + n, p: x.p };
}
export function setprec(x: Real, p: number): Real {
  if (x.p === p) return x;
  if (x.s === 0) return real_0_bit(x.e, p);
  return mkreal(x.s, x.m, rshift(x), p);
}

export function addrr(x: Real, y: Real): Real {
  const p = Math.max(x.p, y.p);
  if (x.s === 0) return setprec(y, p);
  if (y.s === 0) return setprec(x, p);
  const sx = rshift(x);
  const sy = rshift(y);
  const k = Math.min(sx, sy);
  const mx = (x.s < 0 ? -x.m : x.m) << BigInt(sx - k);
  const my = (y.s < 0 ? -y.m : y.m) << BigInt(sy - k);
  const s = mx + my;
  if (s === 0n) return real_0_bit(Math.min(x.e, y.e) - p, p);
  return mkreal(s < 0n ? -1 : 1, iabs(s), k, p);
}
export function subrr(x: Real, y: Real): Real {
  return addrr(x, real_neg(y));
}
export function mulrr(x: Real, y: Real): Real {
  const p = Math.max(x.p, y.p);
  if (x.s === 0 || y.s === 0) return real_0_bit(x.e + y.e, p);
  return mkreal(x.s * y.s, x.m * y.m, rshift(x) + rshift(y), p);
}
export function sqrr(x: Real): Real {
  return mulrr(x, x);
}
export function divrr(x: Real, y: Real): Real {
  if (y.s === 0) throw new PariDomainError('divrr', 'y', '=', '0');
  const p = Math.max(x.p, y.p);
  if (x.s === 0) return real_0_bit(x.e - y.e, p);
  /* compute floor(x.m * 2^(p+2) / y.m) then normalize */
  const shift = p + 2;
  const q = (x.m << BigInt(shift)) / y.m;
  return mkreal(x.s * y.s, q, rshift(x) - rshift(y) - shift, p);
}
export function mulir(x: bigint, y: Real): Real {
  if (x === 0n || y.s === 0) return real_0_bit(y.e, y.p);
  return mkreal((x < 0n ? -1 : 1) * y.s, iabs(x) * y.m, rshift(y), y.p);
}
export function mulur(x: number, y: Real): Real {
  return mulir(BigInt(x), y);
}
/** PARI `divri(x, y)`: real divided by integer */
export function divri(x: Real, y: bigint): Real {
  if (y === 0n) throw new PariDomainError('divri', 'y', '=', '0');
  if (x.s === 0) return real_0_bit(x.e, x.p);
  const p = x.p;
  const shift = p + 2 + expi(iabs(y)) + 1;
  const q = (x.m << BigInt(shift)) / iabs(y);
  return mkreal(x.s * (y < 0n ? -1 : 1), q, rshift(x) - shift, p);
}
export function divru(x: Real, y: number): Real {
  return divri(x, BigInt(y));
}
/** PARI `divir(x, y)`: integer divided by real */
export function divir(x: bigint, y: Real): Real {
  return divrr(itor(x, y.p), y);
}
export function cmprr(x: Real, y: Real): number {
  const d = subrr(x, y);
  return d.s;
}
/** PARI `truncr`: truncate towards 0 */
export function truncr(x: Real): bigint {
  if (x.s === 0) return 0n;
  const k = rshift(x);
  const v = k >= 0 ? x.m << BigInt(k) : x.m >> BigInt(-k);
  return x.s < 0 ? -v : v;
}
/**
 * PARI `gcvtoi(x, &e)` (`gen3.c:2668-2683`): truncate to an integer; `e` is the
 * number of error bits on the integral part (`e > 0` means the result is
 * meaningless).
 */
export function gcvtoi(x: Real): { z: bigint; e: number } {
  if (x.s === 0) return { z: 0n, e: x.e };
  if (x.e < 0) return { z: 0n, e: x.e };
  const e1 = x.e - x.p + 1;
  const z = truncr(x);
  if (e1 > 0) return { z, e: e1 };
  /* e = expo(x - y): exponent of the discarded fractional part */
  const frac = subrr(x, itor(z, x.p));
  return { z, e: frac.s === 0 ? -(1 << 30) : frac.e };
}
export function rtodbl(x: Real): number {
  if (x.s === 0) return 0;
  /* take 53 significant bits */
  const p = x.p;
  let m = x.m;
  let k = rshift(x);
  if (p > 60) {
    const drop = p - 60;
    m >>= BigInt(drop);
    k += drop;
  }
  return x.s * Number(m) * 2 ** k;
}
export function dbltor(d: number, p = DEFAULTPREC): Real {
  if (d === 0) return real_0(p);
  if (!Number.isFinite(d)) throw new PariDomainError('dbltor', 'd', '=', 'oo');
  const s = d < 0 ? -1 : 1;
  let a = Math.abs(d);
  let k = 0;
  while (a < 1) {
    a *= 2;
    k--;
  }
  while (a >= 2) {
    a /= 2;
    k++;
  }
  /* a in [1,2): 53 bits */
  const m = BigInt(Math.round(a * 2 ** 52));
  return mkreal(s, m, k - 52, p);
}

/** PARI `sqrtr` for `x > 0` */
export function sqrtr(x: Real): Real {
  if (x.s < 0) throw new PariDomainError('sqrtr', 'x', '<', '0');
  if (x.s === 0) return real_0_bit(x.e >> 1, x.p);
  const p = x.p;
  /* x = m * 2^k; want sqrt = sqrt(m * 2^(k+2t)) * 2^(-t) with enough bits */
  let m = x.m;
  let k = rshift(x);
  const extra = p + 4;
  m <<= BigInt(2 * extra);
  k -= 2 * extra;
  if (k % 2 !== 0) {
    m <<= 1n;
    k -= 1;
  }
  const r = sqrti(m);
  return mkreal(1, r, k / 2, p);
}

let LOG2_CACHE: Real | null = null;
/** PARI `mplog2(prec)` */
export function mplog2(p = DEFAULTPREC): Real {
  if (LOG2_CACHE && LOG2_CACHE.p >= p) return setprec(LOG2_CACHE, p);
  /* log 2 = 2 * atanh(1/3) */
  const w = p + 32;
  const one = real_1(w);
  const third = divru(one, 3);
  LOG2_CACHE = shiftr(atanh_small(third, w), 1);
  return setprec(LOG2_CACHE, p);
}

/** `atanh(t)` for `|t| <= 1/3` by its Taylor series, at `p` bits */
function atanh_small(t: Real, p: number): Real {
  const t2 = mulrr(setprec(t, p), setprec(t, p));
  let term = setprec(t, p);
  let sum = term;
  for (let n = 3; ; n += 2) {
    term = mulrr(term, t2);
    if (term.s === 0 || term.e < sum.e - p - 4) break;
    sum = addrr(sum, divru(term, n));
  }
  return sum;
}

/**
 * PARI `logr_abs(x)`: natural logarithm of `|x|`, `x != 0`.
 *
 * Argument reduction `x = m * 2^k` with `m in [1,2)`, then `sqrt` until
 * `m` is close to 1 and the `atanh` series converges fast.
 */
export function logr_abs(x: Real): Real {
  if (x.s === 0) throw new PariDomainError('logr_abs', 'x', '=', '0');
  const p = x.p;
  const w = p + 32;
  const k = x.e;
  let m = setprec({ s: 1, m: x.m, e: 0, p: x.p }, w); /* |x| / 2^k in [1,2) */
  let nsq = 0;
  /* sqrt until m - 1 <= 1/4, i.e. expo(m-1) <= -2 */
  for (;;) {
    const d = subrr(m, real_1(w));
    if (d.s === 0 || d.e <= -3) break;
    m = sqrtr(m);
    nsq++;
    if (nsq > 4 * w) break;
  }
  /* log m = 2 atanh((m-1)/(m+1)) */
  const num = subrr(m, real_1(w));
  const den = addrr(m, real_1(w));
  const t = num.s === 0 ? real_0(w) : divrr(num, den);
  let lm = shiftr(atanh_small(t, w), 1);
  lm = shiftr(lm, nsq); /* undo the nsq square roots */
  const res = k === 0 ? lm : addrr(lm, mulir(BigInt(k), mplog2(w)));
  return setprec(res, p);
}

/** `exp(x)` (only used by the tests as an inverse oracle for `logr_abs`) */
export function expr(x: Real): Real {
  const p = x.p;
  const w = p + 32;
  if (x.s === 0) return real_1(p);
  /* x = k log 2 + r, |r| <= log2/2 */
  const l2 = mplog2(w);
  const kk = gcvtoi(divrr(setprec(x, w), l2)).z;
  let r = subrr(setprec(x, w), mulir(kk, l2));
  /* halve r until small */
  let nh = 0;
  while (r.s !== 0 && r.e > -8) {
    r = shiftr(r, -1);
    nh++;
  }
  let term = real_1(w);
  let sum = real_1(w);
  for (let n = 1; ; n++) {
    term = divru(mulrr(term, r), n);
    if (term.s === 0 || term.e < sum.e - w - 4) break;
    sum = addrr(sum, term);
  }
  for (let i = 0; i < nh; i++) sum = mulrr(sum, sum);
  return setprec(shiftr(sum, Number(kk)), p);
}

/**
 * A PARI `t_COMPLEX` with `t_REAL` real part and `t_INT` imaginary part: the
 * archimedean component attached to a relation in `buch1.c` (a Shanks distance
 * plus a bit recording the sign of the corresponding quadratic number).
 */
export interface CReal {
  re: Real;
  im: bigint;
}
function creal(re: Real, im: bigint): CReal {
  return { re, im };
}
function cadd(x: CReal, y: CReal): CReal {
  return creal(addrr(x.re, y.re), x.im + y.im);
}
function csub(x: CReal, y: CReal): CReal {
  return creal(subrr(x.re, y.re), x.im - y.im);
}
function cmulint(x: CReal, n: bigint): CReal {
  return creal(mulir(n, x.re), n * x.im);
}
function cneg(x: CReal): CReal {
  return creal(real_neg(x.re), -x.im);
}

/* ================================================================== */
/* ZM_pivots / row rank profile (alglin1.c:3695)                       */
/* ================================================================== */

/**
 * PARI `ZM_pivots(M0, &rr)` (`alglin1.c:3695-3760`).
 *
 * `d[k]` is the row index of the pivot used in column `k`, or `0` when column
 * `k` is a linear combination of the previous ones; `rr` is the number of such
 * columns, i.e. `dim ker M0`.
 *
 * Upstream computes this modulo random word-size primes and then *certifies*
 * the answer by an exact linear-algebra check, so its output is the true (and
 * canonical: pivot = first unused row) rank profile.  We compute the same
 * object directly by one-step fraction-free (Bareiss) elimination, which is
 * exact by construction; the returned `d` and `rr` are identical.
 */
export function ZM_pivots(M0: ZMat): { d: zv | null; rr: number } {
  const n = lg(M0) - 1;
  if (n === 0) return { d: null, rr: 0 };
  const m = nbrows(M0);
  /* working copy, row-major, 1-based */
  const x: bigint[][] = new Array(m + 1);
  for (let i = 1; i <= m; i++) {
    const row = new Array<bigint>(n + 1).fill(0n);
    for (let j = 1; j <= n; j++) row[j] = M0[j]![i]!;
    x[i] = row;
  }
  const c = zero_zv(m); /* c[i] = column whose pivot sits on row i */
  const d = cgetg_zv(n);
  let r = 0;
  let prev = 1n;
  for (let k = 1; k <= n; k++) {
    let j = 0;
    for (let i = 1; i <= m; i++)
      if (!c[i] && x[i]![k]) {
        j = i;
        break;
      }
    if (!j) {
      r++;
      d[k] = 0;
      continue;
    }
    c[j] = k;
    d[k] = j;
    const piv = x[j]![k]!;
    for (let t = 1; t <= m; t++) {
      if (t === j || c[t]) continue;
      const xt = x[t]!;
      const v = xt[k]!;
      if (!v) {
        /* still need the Bareiss scaling to keep the invariant */
        for (let i = k + 1; i <= n; i++) {
          const num = piv * xt[i]!;
          xt[i] = num / prev;
        }
        xt[k] = 0n;
        continue;
      }
      const xj = x[j]!;
      for (let i = k + 1; i <= n; i++) {
        const num = piv * xt[i]! - v * xj[i]!;
        if (num % prev !== 0n) throw new Error('ZM_pivots: inexact Bareiss division');
        xt[i] = num / prev;
      }
      xt[k] = 0n;
    }
    prev = piv;
  }
  return { d, rr: r };
}

/**
 * PARI `ZM_rowrankprofile(x, &nlze)` (`hnf_snf.c:222-239`): permutation giving
 * `imagecompl(x') | image(x')` where `x'` is the transpose of `x`.
 */
function ZM_rowrankprofile(x: ZMat): { perm: zv; nlze: number } {
  const l = nbrows(x) + 1; /* lg of the transpose */
  /* transpose */
  const t: ZMat = new Array<ZC>(l);
  t[0] = [];
  for (let i = 1; i < l; i++) {
    const col = cgetg_int(lg(x) - 1);
    for (let j = 1; j < lg(x); j++) col[j] = x[j]![i]!;
    t[i] = col;
  }
  const { d, rr } = ZM_pivots(t);
  if (!d) return { perm: identity_perm(l - 1), nlze: rr };
  const y = cgetg_zv(l - 1);
  let j = 1;
  let k = rr + 1;
  for (let i = 1; i < l; i++) {
    if (d[i]) y[k++] = i;
    else y[j++] = i;
  }
  return { perm: y, nlze: rr };
}

/* ================================================================== */
/* HNFLLL (Havas, Majewski, Mathews) -- hnf_snf.c:1546-1810           */
/* ================================================================== */

/** `hnf_snf.c:1547-1554` */
function Minus(j: number, lambda: ZMat): void {
  const n = lg(lambda);
  for (let k = 1; k < j; k++) lambda[j]![k] = -lambda[j]![k]!;
  for (let k = j + 1; k < n; k++) lambda[k]![j] = -lambda[k]![j]!;
}

/** index of first nonzero entry (`hnf_snf.c:1557-1564`) */
function findi(M: ZC): number {
  const n = lg(M);
  for (let i = 1; i < n; i++) if (M[i]) return i;
  return 0;
}

/** `hnf_snf.c:1566-1576` */
function findi_normalize(Aj: ZC, B: ZMat | null, j: number, lambda: ZMat): number {
  const r = findi(Aj);
  if (r && Aj[r]! < 0n) {
    ZV_togglesign(Aj);
    if (B) ZV_togglesign(B[j]!);
    Minus(j, lambda);
  }
  return r;
}

/** `hnf_snf.c:1624-1670` */
function reduce2(
  A: ZMat,
  B: ZMat | null,
  k: number,
  j: number,
  lambda: ZMat,
  D: bigint[]
): { row0: number; row1: number } {
  let q: bigint;
  const row0 = findi_normalize(A[j]!, B, j, lambda);
  const row1 = findi_normalize(A[k]!, B, k, lambda);
  if (row0) {
    q = truedivii(A[k]![row0]!, A[j]![row0]!);
  } else {
    const cmp = iabs(shifti(lambda[k]![j]!, 1)) - iabs(D[j]!);
    if (cmp > 0n) q = diviiround(lambda[k]![j]!, D[j]!);
    else return { row0, row1 };
  }
  if (q) {
    const Lk = lambda[k]!;
    const Lj = lambda[j]!;
    q = -q;
    if (row0) ZC_lincomb1_inplace(A[k]!, A[j]!, q);
    if (B) ZC_lincomb1_inplace(B[k]!, B[j]!, q);
    Lk[j] = Lk[j]! + q * D[j]!;
    for (let i = 1; i < j; i++) if (Lj[i]) Lk[i] = Lk[i]! + q * Lj[i]!;
  }
  return { row0, row1 };
}

/** `hnf_snf.c:1698-1725` */
function hnfswap(A: ZMat, B: ZMat | null, k: number, lambda: ZMat, D: bigint[]): void {
  const n = lg(A);
  const Lk = lambda[k]!;
  const tmp = A[k]!;
  A[k] = A[k - 1]!;
  A[k - 1] = tmp;
  if (B) {
    const t = B[k]!;
    B[k] = B[k - 1]!;
    B[k - 1] = t;
  }
  for (let j = k - 2; j; j--) {
    const t = lambda[k - 1]![j]!;
    lambda[k - 1]![j] = Lk[j]!;
    Lk[j] = t;
  }
  for (let i = k + 1; i < n; i++) {
    const Li = lambda[i]!;
    if (Li[k - 1] === 0n && Li[k] === 0n) continue;
    const t = Li[k - 1]! * D[k]! - Li[k]! * Lk[k - 1]!;
    const u = Li[k]! * D[k - 2]! + Li[k - 1]! * Lk[k - 1]!;
    Li[k - 1] = u / D[k - 1]!;
    Li[k] = t / D[k - 1]!;
  }
  D[k - 1] = (D[k - 2]! * D[k]! + Lk[k - 1]! * Lk[k - 1]!) / D[k - 1]!;
}

/** `hnf_snf.c:1729-1744` */
function reverse_rows(A: ZMat): ZMat {
  const n = lg(A);
  if (n === 1) return A;
  const h = lg(A[1]!);
  for (let j = 1; j < n; j++) {
    const col = A[j]!;
    for (let i = (h - 1) >> 1; i; i--) {
      const t = col[i]!;
      col[i] = col[h - i]!;
      col[h - i] = t;
    }
  }
  return A;
}

/** `hnf_snf.c:1747-1753` */
function must_swap(k: number, lambda: ZMat, D: bigint[]): boolean {
  const z = D[k - 2]! * D[k]! + lambda[k]![k - 1]! * lambda[k]![k - 1]!;
  return z < D[k - 1]! * D[k - 1]!;
}

/**
 * PARI `ZM_hnflll(A, &B, remove)` (`hnf_snf.c:1755-1810`).
 * Returns the HNF of `A` (column HNF) and, if `wantB`, the transformation `B`
 * with `A_orig * B = A_hnf`.
 */
export function ZM_hnflll(A0: ZMat, wantB: boolean, remove: boolean): { H: ZMat; B: ZMat | null } {
  const n = lg(A0);
  let A = reverse_rows(shallowcopyMat(A0));
  let B: ZMat | null = wantB ? matid(n - 1) : null;
  const D: bigint[] = new Array<bigint>(n + 1).fill(1n); /* D[0..n-1] */
  const lambda = zeromat(n - 1, n - 1);
  let k = 2;
  while (k < n) {
    const { row0, row1 } = reduce2(A, B, k, k - 1, lambda, D);
    let do_swap: boolean;
    if (row0) do_swap = !row1 || row0 <= row1;
    else if (row1) do_swap = false;
    else do_swap = must_swap(k, lambda, D);
    if (do_swap) {
      hnfswap(A, B, k, lambda, D);
      if (k > 2) k--;
    } else {
      for (let i = k - 2; i; i--) reduce2(A, B, k, i, lambda, D);
      k++;
    }
  }
  if (n === 2) findi_normalize(A[1]!, B, 1, lambda);
  A = reverse_rows(A);
  if (remove) {
    let i = 1;
    for (; i < n; i++) if (!ZV_equal0(A[i]!)) break;
    /* remove_0cols(i-1, &A, &B, remove) */
    const t = i - 1;
    if (t) {
      A = vecslice(A, t + 1, n - 1);
      if (B) B = vecslice(B, t + 1, n - 1);
    }
  }
  return { H: A, B };
}

/* ================================================================== */
/* Special HNF for relation matrices -- hnf_snf.c:70-658               */
/* ================================================================== */

/** `hnf_snf.c:70-85`: number of nonzero entries in a row, `-1` if some |entry| > 1 */
function count(mat: zv[], row: number, len: number): { n: number; first: number } {
  let n = 0;
  let first = 0;
  for (let j = 1; j <= len; j++) {
    const p = mat[j]![row]!;
    if (p) {
      if (Math.abs(p) !== 1) return { n: -1, first };
      n++;
      first = j;
    }
  }
  return { n, first };
}

/** `hnf_snf.c:87-94`: last column with a `+/-1` on that row */
function count2(mat: zv[], row: number, len: number): number {
  for (let j = len; j; j--) if (Math.abs(mat[j]![row]!) === 1) return j;
  return 0;
}

/** the archimedean component vector attached to the relations */
type CVec = CReal[];

function RgV_ZM_mul(C: CVec, T: ZMat): CVec {
  const n = lg(T) - 1;
  const out: CVec = new Array<CReal>(n + 1);
  out[0] = C[0]!;
  for (let j = 1; j <= n; j++) {
    const t = T[j]!;
    let s: CReal | null = null;
    for (let i = 1; i < lg(t); i++) {
      const v = t[i]!;
      if (!v) continue;
      const z = cmulint(C[i]!, v);
      s = s ? cadd(s, z) : z;
    }
    out[j] = s ?? creal(real_0(C[1]!.re.p), 0n);
  }
  return out;
}
function RgV_zm_mul(C: CVec, T: zv[]): CVec {
  const n = lg(T) - 1;
  const out: CVec = new Array<CReal>(n + 1);
  out[0] = C[0]!;
  for (let j = 1; j <= n; j++) {
    const t = T[j]!;
    let s: CReal | null = null;
    for (let i = 1; i < lg(t); i++) {
      const v = t[i]!;
      if (!v) continue;
      const z = cmulint(C[i]!, BigInt(v));
      s = s ? cadd(s, z) : z;
    }
    out[j] = s ?? creal(real_0(C[1]!.re.p), 0n);
  }
  return out;
}
function CVec_sub(A: CVec, B: CVec): CVec {
  const out: CVec = new Array<CReal>(lg(A));
  out[0] = A[0]!;
  for (let i = 1; i < lg(A); i++) out[i] = csub(A[i]!, B[i]!);
  return out;
}

interface HnfState {
  dep: ZMat;
  B: ZMat;
  C: CVec;
}

/**
 * PARI `hnffinal(matgen, perm, &dep, &B, &C)` (`hnf_snf.c:96-201`).
 * `perm` is modified in place.
 */
function hnffinal(matgen: ZMat, perm: zv, st: HnfState): ZMat {
  const B = st.B;
  const C = st.C;
  let dep = st.dep;
  const co = lg(C);
  let col = lg(matgen) - 1;
  if (col === 0) return matgen;
  let lnz = nbrows(matgen);
  const nlze = nbrows(dep);
  let lig = nlze + lnz;

  const hl = ZM_hnflll(matgen, true, false);
  let H = hl.H;
  const U = hl.B!;
  /* H += lg(H)-1 - lnz : keep the last lnz columns */
  H = vecslice(H, lg(H) - lnz, lg(H) - 1);
  const zc = col - lnz;
  if (nlze) {
    dep = ZM_mul(dep, U);
    dep = vecslice(dep, zc + 1, lg(dep) - 1);
  }
  const diagH1: boolean[] = new Array<boolean>(lnz + 1).fill(false);

  const Cnew: CVec = new Array<CReal>(co);
  Cnew[0] = C[0]!;
  {
    const Ccut = vecslice(C, 1, col);
    const p1 = RgV_ZM_mul(Ccut, U);
    for (let j = 1; j <= col; j++) Cnew[j] = p1[j]!;
    for (let j = col + 1; j < co; j++) Cnew[j] = C[j]!;
  }

  /* Clean up B using new H */
  let s = 0;
  for (let i = lnz; i; i--) {
    const Di = dep[i];
    const Hi = H[i]!;
    let h: bigint | null = Hi[i]!;
    diagH1[i] = h === 1n || h === -1n;
    if (diagH1[i]) {
      h = null;
      s++;
    }
    for (let j = col + 1; j < co; j++) {
      const z = B[j - col]!;
      let p1 = z[i + nlze]!;
      if (h) p1 = truedivii(p1, h);
      if (!p1) continue;
      let k = 1;
      for (; k <= nlze; k++) z[k] = z[k]! - p1 * Di![k]!;
      for (; k <= lig; k++) z[k] = z[k]! - p1 * Hi[k - nlze]!;
      Cnew[j] = csub(Cnew[j]!, cmulint(Cnew[i + zc]!, p1));
    }
  }
  /* push the 1 rows down */
  {
    const p1 = cgetg_zv(lnz);
    let i1 = 0;
    let j1 = lnz - s;
    for (let i = 1; i <= lnz; i++) {
      if (diagH1[i]) p1[++j1] = perm[nlze + i]!;
      else perm[nlze + ++i1] = perm[nlze + i]!;
    }
    for (let i = i1 + 1; i <= lnz; i++) perm[nlze + i] = p1[i]!;
  }

  lig -= s;
  col -= s;
  const lnz0 = lnz;
  lnz -= s;
  const Hnew: ZMat = new Array<ZC>(lnz + 1);
  Hnew[0] = [];
  const depnew: ZMat = new Array<ZC>(lnz + 1);
  depnew[0] = [];
  const Bnew: ZMat = new Array<ZC>(co - col);
  Bnew[0] = [];
  const Cout: CVec = Cnew.slice();
  {
    let i1 = 0;
    let j1 = 0;
    for (let j = 1; j <= lnz0; j++) {
      const z = H[j]!;
      let p1: ZC;
      let off = 0;
      if (diagH1[j]) {
        i1++;
        Cout[i1 + col] = Cnew[j + zc]!;
        p1 = cgetg_int(lig);
        Bnew[i1] = p1;
        for (let i = 1; i <= nlze; i++) p1[i] = dep[j]![i]!;
        off = nlze;
      } else {
        j1++;
        Cout[j1 + zc] = Cnew[j + zc]!;
        p1 = cgetg_int(lnz);
        Hnew[j1] = p1;
        depnew[j1] = dep[j]!;
      }
      let k = 1;
      for (let i = 1; k <= lnz; i++) if (!diagH1[i]) p1[off + k++] = z[i]!;
    }
    for (let j = s + 1; j < co - col; j++) {
      const z = B[j - s]!;
      const p1 = cgetg_int(lig);
      Bnew[j] = p1;
      for (let i = 1; i <= nlze; i++) p1[i] = z[i]!;
      let k = 1;
      for (let i = 1; k <= lnz; i++) if (!diagH1[i]) p1[nlze + k++] = z[nlze + i]!;
    }
  }
  st.dep = depnew;
  st.C = Cout;
  st.B = Bnew;
  return Hnew;
}

/**
 * PARI `hnfspec_i(mat0, perm, &dep, &B, &C, k0)` (`hnf_snf.c:254-554`).
 *
 * HNF-reduce a relation matrix by column operations and row permutations.
 * `mat0` is a `t_VECSMALL` matrix, `perm` a permutation of its rows (modified
 * in place), `k0` the number of dense top rows.
 */
export function hnfspec_i(mat0: zv[], perm: zv, st: HnfState, k0: number): ZMat {
  const li = lg(perm);
  const CO = lg(mat0);
  let C = st.C;
  let co = CO;
  if (co > 300 && co > 1.5 * li) {
    /* treat the rest at the end; unreachable from Buchquad, whose relation
     * matrix has co ~ li */
    co = Math.floor(1.2 * li);
    C = vecslice(C, 1, co - 1);
  }
  const matt: ZMat = new Array<ZC>(co);
  matt[0] = [];
  const mat: zv[] = new Array<zv>(co);
  mat[0] = [];
  for (let j = 1; j < co; j++) {
    const matj = mat0[j]!.slice();
    const p1 = cgetg_int(k0);
    matt[j] = p1;
    mat[j] = matj;
    for (let i = 1; i <= k0; i++) p1[i] = BigInt(matj[perm[i]!]!);
  }

  let i = li - 1;
  let lig = li - 1;
  let col = co - 1;
  let lk0 = k0;
  const T: ZMat | null = k0 || lg(C) > 1 ? matid(col) : null;
  /* Look for lines with a single nonzero entry, equal to 1 in absolute value */
  while (i > lk0 && col) {
    const { n, first } = count(mat, perm[i]!, col);
    if (n === 0) {
      lk0++;
      const t = perm[i]!;
      perm[i] = perm[lk0]!;
      perm[lk0] = t;
      i = lig;
      continue;
    }
    if (n === 1) {
      const t = perm[i]!;
      perm[i] = perm[lig]!;
      perm[lig] = t;
      if (T) {
        const u = T[first]!;
        T[first] = T[col]!;
        T[col] = u;
      }
      const u = mat[first]!;
      mat[first] = mat[col]!;
      mat[col] = u;
      const p = mat[col]!;
      if (p[perm[lig]!]! < 0) {
        for (let h = lk0 + 1; h < lig; h++) p[perm[h]!] = -p[perm[h]!]!;
        if (T) {
          const p1 = T[col]!;
          for (let h = 1; ; h++)
            if (p1[h]) {
              p1[h] = -p1[h]!;
              break;
            }
        }
      }
      lig--;
      col--;
      i = lig;
      continue;
    }
    i--;
  }

  /* Get rid of all lines containing only 0 and +/- 1 */
  let s = 0;
  while (lig > lk0 && col && s < 0x40000000) {
    let ii = lig;
    let nn = 0;
    for (; ii > lk0; ii--) {
      const r = count(mat, perm[ii]!, col);
      if (r.n > 0) {
        nn = r.first;
        break;
      }
    }
    if (ii === lk0) break;
    {
      const t = perm[ii]!;
      perm[ii] = perm[lig]!;
      perm[lig] = t;
    }
    {
      const u = mat[nn]!;
      mat[nn] = mat[col]!;
      mat[col] = u;
    }
    const p = mat[col]!;
    if (T) {
      const u = T[nn]!;
      T[nn] = T[col]!;
      T[col] = u;
    }
    if (p[perm[lig]!]! < 0) {
      for (let h = lk0 + 1; h <= lig; h++) p[perm[h]!] = -p[perm[h]!]!;
      if (T) ZV_togglesign(T[col]!);
    }
    for (let j = 1; j < col; j++) {
      const matj = mat[j]!;
      const t = matj[perm[lig]!]!;
      if (!t) continue;
      if (t === 1) {
        for (let h = lk0 + 1; h <= lig; h++) {
          matj[perm[h]!] = matj[perm[h]!]! - p[perm[h]!]!;
          const z = Math.abs(matj[perm[h]!]!);
          if (z > s) s = z;
        }
      } else {
        for (let h = lk0 + 1; h <= lig; h++) {
          matj[perm[h]!] = matj[perm[h]!]! + p[perm[h]!]!;
          const z = Math.abs(matj[perm[h]!]!);
          if (z > s) s = z;
        }
      }
      if (T) ZC_lincomb1_inplace(T[j]!, T[col]!, BigInt(-t));
    }
    lig--;
    col--;
  }

  /* As above with lines containing a +/- 1 (no other assumption) */
  const vmax = cgetg_zv(co - 1);
  for (let j = 1; j <= col; j++) {
    const matj = mat[j]!;
    let sm = 0;
    for (let h = lk0 + 1; h <= lig; h++) {
      const z = Math.abs(matj[h]!);
      if (z > sm) sm = z;
    }
    vmax[j] = sm;
  }
  END2: while (lig > lk0 && col) {
    let ii = lig;
    let nn = 0;
    for (; ii > lk0; ii--) {
      nn = count2(mat, perm[ii]!, col);
      if (nn) break;
    }
    if (ii === lk0) break;
    {
      const t = vmax[nn]!;
      vmax[nn] = vmax[col]!;
      vmax[col] = t;
    }
    {
      const t = perm[ii]!;
      perm[ii] = perm[lig]!;
      perm[lig] = t;
    }
    {
      const u = mat[nn]!;
      mat[nn] = mat[col]!;
      mat[col] = u;
    }
    const p = mat[col]!;
    if (T) {
      const u = T[nn]!;
      T[nn] = T[col]!;
      T[col] = u;
    }
    if (p[perm[lig]!]! < 0) {
      for (let h = lk0 + 1; h <= lig; h++) p[perm[h]!] = -p[perm[h]!]!;
      if (T) ZV_togglesign(T[col]!);
    }
    for (let j = 1; j < col; j++) {
      const matj = mat[j]!;
      const t = matj[perm[lig]!]!;
      if (!t) continue;
      /* single precision guard, as upstream (HIGHBIT for 32-bit longs is
       * replaced by 2^52, the exact range of a JS number) */
      if (vmax[col] && Math.abs(t) >= (4503599627370496 - vmax[j]!) / vmax[col]!) break END2;
      let sm = 0;
      for (let h = lk0 + 1; h <= lig; h++) {
        matj[perm[h]!] = matj[perm[h]!]! - t * p[perm[h]!]!;
        const z = Math.abs(matj[perm[h]!]!);
        if (z > sm) sm = z;
      }
      vmax[j] = sm;
      if (T) ZC_lincomb1_inplace(T[j]!, T[col]!, BigInt(-t));
    }
    lig--;
    col--;
  }

  /* clean up mat: remove everything to the right of the 1s on diagonal */
  const matb: ZMat = new Array<ZC>(co);
  matb[0] = [];
  for (let j = 1; j < co; j++) {
    const matj = mat[j]!;
    const p1 = cgetg_int(li - 1 - k0); /* indices 1..li-1-k0 <-> rows k0+1..li-1 */
    matb[j] = p1;
    for (let h = k0 + 1; h < li; h++) p1[h - k0] = BigInt(matj[perm[h]!]!);
  }
  for (let ii = li - 2; ii > lig; ii--) {
    const i0 = ii - k0;
    const k = ii + co - li;
    const Bk = matb[k]!;
    for (let j = k + 1; j < co; j++) {
      const Bj = matb[j]!;
      const v = Bj[i0]!;
      if (!v) continue;
      Bj[i0] = 0n;
      for (let h = 1; h < i0; h++) Bj[h] = Bj[h]! - v * Bk[h]!;
      if (T) ZC_lincomb1_inplace(T[j]!, T[k]!, -v);
    }
  }

  const nlze = lk0 - k0;
  const lnz = lig - nlze + 1;
  const mattT = T ? ZM_mul(matt, T) : matt;
  const extramat: ZMat = new Array<ZC>(col + 1);
  extramat[0] = [];
  for (let j = 1; j <= col; j++) {
    const z = mattT[j]!;
    const t = matb[j]!;
    const p2 = cgetg_int(lnz - 1);
    extramat[j] = p2;
    let ii = 1;
    for (; ii <= k0; ii++) p2[ii] = z[ii]!;
    for (; ii < lnz; ii++) p2[ii] = t[ii + nlze - k0]!;
  }
  let permpro: zv;
  let nr: number;
  if (!col) {
    permpro = identity_perm(lnz - 1);
    nr = lnz - 1;
  } else {
    const rp = ZM_rowrankprofile(extramat);
    permpro = rp.perm;
    nr = rp.nlze;
  }
  if (nlze) {
    const p1 = cgetg_zv(lk0);
    for (let ii = 1; ii <= nlze; ii++) p1[ii] = perm[ii + k0]!;
    for (let ii = nlze + 1; ii <= lk0; ii++) p1[ii] = perm[ii - nlze]!;
    for (let ii = 1; ii <= lk0; ii++) perm[ii] = p1[ii]!;
  }
  {
    const p1 = cgetg_zv(lnz - 1);
    for (let ii = 1; ii < lnz; ii++) p1[ii] = perm[nlze + permpro[ii]!]!;
    for (let ii = 1; ii < lnz; ii++) perm[nlze + ii] = p1[ii]!;
  }

  const matbnew: ZMat = new Array<ZC>(col + 1);
  matbnew[0] = [];
  const dep: ZMat = new Array<ZC>(col + 1);
  dep[0] = [];
  for (let j = 1; j <= col; j++) {
    const z = extramat[j]!;
    const p1 = cgetg_int(nlze + nr);
    dep[j] = p1;
    const p2 = cgetg_int(lnz - 1 - nr);
    matbnew[j] = p2;
    let ii = 1;
    for (; ii <= nlze; ii++) p1[ii] = 0n;
    for (ii = 1; ii <= nr; ii++) p1[nlze + ii] = z[permpro[ii]!]!;
    for (; ii < lnz; ii++) p2[ii - nr] = z[permpro[ii]!]!;
  }

  const Bmat: ZMat = new Array<ZC>(co - col);
  Bmat[0] = [];
  for (let j = col + 1; j < co; j++) {
    const y = mattT[j]!;
    const z = matb[j]!;
    const p1 = cgetg_int(lig);
    Bmat[j - col] = p1;
    for (let ii = 1; ii <= nlze; ii++) p1[ii] = z[ii]!;
    for (let k = 1; k < lnz; k++) {
      const ii = permpro[k]!;
      /* upstream: z += nlze-k0 before this loop */
      p1[nlze + k] = ii <= k0 ? y[ii]! : z[ii + nlze - k0]!;
    }
  }
  if (T) C = RgV_ZM_mul(C, T);
  st.dep = dep;
  st.B = Bmat;
  st.C = C;
  let H = hnffinal(matbnew, perm, st);
  if (CO > co) {
    /* treat the rest, N columns at a time (unreachable from Buchquad) */
    const N = 300;
    const L = CO - co;
    let l = Math.min(L, N);
    let off = co - 1;
    for (let a = l; ; ) {
      const MAT: zv[] = new Array<zv>(l + 1);
      MAT[0] = [];
      const emb: CVec = new Array<CReal>(l + 1);
      emb[0] = st.C[0]!;
      for (let j = 1; j <= l; j++) {
        MAT[j] = mat0[off + j]!;
        emb[j] = st.C[off + j]!;
      }
      H = hnfadd_i(H, perm, st, MAT, emb);
      if (a === L) break;
      off += l;
      a += l;
      if (a > L) {
        l = L - (a - l);
        a = L;
      }
    }
  }
  return H;
}

/**
 * PARI `hnfadd_i(H, perm, &dep, &B, &C, extramat, extraC)`
 * (`hnf_snf.c:610-658`): add new relations to a matrix treated by `hnfspec`.
 */
export function hnfadd_i(H: ZMat, perm: zv, st: HnfState, extramat: zv[], extraC: CVec): ZMat {
  if (lg(extramat) === 1) return H;
  const B = st.B;
  const C = st.C;
  const dep = st.dep;
  const co = lg(C) - 1;
  const lH = lg(H) - 1;
  const lB = lg(B) - 1;
  const li = lg(perm) - 1;
  const lig = li - lB;
  const col = co - lB;

  let extratop = zm_to_ZM(rowslicepermute(extramat, perm, 1, lig));
  let extraCC = extraC;
  if (li !== lig) {
    const A = vecslice(C, col + 1, co);
    const c = rowslicepermute(extramat, perm, lig + 1, li);
    extraCC = CVec_sub(extraC, RgV_zm_mul(A, c));
    extratop = ZM_sub(extratop, ZM_zm_mul(B, c));
  }
  const extra = shallowconcat(extratop, vconcat(dep, H));
  let Cnew = shallowconcat(extraCC, vecslice(C, col - lH + 1, co));
  const rp = ZM_rowrankprofile(extra);
  const permpro0 = rp.perm;
  const nlze = rp.nlze;
  const extra2 = rowpermute(extra, permpro0);
  st.B = rowpermute(B, permpro0);
  const permpro = vecsmallpermute(perm, permpro0);
  for (let i = 1; i <= lig; i++) perm[i] = permpro[i]!;

  st.dep = rowslice(extra2, 1, nlze);
  const matb = rowslice(extra2, nlze + 1, lig);
  st.C = Cnew;
  const Hnew = hnffinal(matb, perm, st);
  Cnew = st.C;
  st.C = shallowconcat(vecslice(C, 1, col - lH), Cnew);
  return Hnew;
}

/* ================================================================== */
/* Smith normal form -- hnf_snf.c:2330-2900                            */
/* ================================================================== */

/** PARI `ZC_elem` (`hnf_snf.c:671-713`) */
function ZC_elem(aj: bigint, ak: bigint, A: ZMat, U: ZMat | null, j: number, k: number): void {
  if (!ak) {
    const t = A[j]!;
    A[j] = A[k]!;
    A[k] = t;
    if (U) {
      const u = U[j]!;
      U[j] = U[k]!;
      U[k] = u;
    }
    return;
  }
  const [d, u, v] = bezout(aj, ak);
  if (!u) {
    /* ak | aj */
    const p1 = -(aj / ak);
    ZC_lincomb1_inplace(A[j]!, A[k]!, p1);
    if (U) ZC_lincomb1_inplace(U[j]!, U[k]!, p1);
    return;
  }
  if (!v) {
    /* aj | ak */
    const p1 = -(ak / aj);
    ZC_lincomb1_inplace(A[k]!, A[j]!, p1);
    const t = A[j]!;
    A[j] = A[k]!;
    A[k] = t;
    if (U) {
      ZC_lincomb1_inplace(U[k]!, U[j]!, p1);
      const uu = U[j]!;
      U[j] = U[k]!;
      U[k] = uu;
    }
    return;
  }
  let aj2 = aj;
  let ak2 = ak;
  if (d !== 1n && d !== -1n) {
    aj2 = aj / d;
    ak2 = ak / d;
  }
  aj2 = -aj2;
  const lincomb = (x: bigint, y: bigint, X: ZC, Y: ZC): ZC => {
    const out = cgetg_int(lg(X) - 1);
    for (let i = 1; i < lg(X); i++) out[i] = x * X[i]! + y * Y[i]!;
    return out;
  };
  const p1 = A[k]!;
  A[k] = lincomb(u, v, A[j]!, p1);
  A[j] = lincomb(aj2, ak2, p1, A[j]!);
  if (U) {
    const q1 = U[k]!;
    U[k] = lincomb(u, v, U[j]!, q1);
    U[j] = lincomb(aj2, ak2, q1, U[j]!);
  }
}

/** PARI `bezout_step` (`hnf_snf.c:2329-2346`) */
function bezout_step(
  a: bigint,
  b: bigint
): { d: bigint; a: bigint; b: bigint; u: bigint; v: bigint } {
  if (iabs(a) === iabs(b)) {
    const sa = signe(a);
    const sb = signe(b);
    if (sb === sa) {
      if (sa > 0) return { d: a, a: 1n, b: 1n, u: 1n, v: 0n };
      return { d: iabs(a), a: 1n, b: 1n, u: -1n, v: 0n };
    }
    if (sa > 0) return { d: a, a: 1n, b: -1n, u: 1n, v: 0n };
    return { d: b, a: -1n, b: 1n, u: -1n, v: 0n };
  }
  const [d, u, v] = bezout(a, b);
  return { d, a: a / d, b: b / d, u, v };
}

/** PARI `ZM_snf_no_divide` (`hnf_snf.c:2356-2368`) */
function ZM_snf_no_divide(x: ZMat, i: number): number {
  const b = x[i]![i]!;
  if (b === 1n || b === -1n) return 0;
  for (let k = 1; k < i; k++) for (let j = 1; j < i; j++) if (x[j]![k]! % b !== 0n) return k;
  return 0;
}

/** PARI `ZM_redpart` (`hnf_snf.c:2371-2381`) */
function ZM_redpart(x: ZMat, p: bigint, I: number): void {
  for (let j = 1; j <= I; j++) {
    const col = x[j]!;
    for (let i = 1; i <= I; i++) {
      const c = col[i]!;
      if (iabs(c) > p) col[i] = c % p;
    }
  }
}

/** PARI `update` (`hnf_snf.c:775-792`): `(c1, c2) *= [u,-b; v,a]` */
function updateCols(
  u: bigint,
  v: bigint,
  a: bigint,
  b: bigint,
  U: ZMat,
  i: number,
  j: number
): void {
  const c1 = U[i]!;
  const c2 = U[j]!;
  const n = lg(c1);
  const p1 = cgetg_int(n - 1);
  const p2 = cgetg_int(n - 1);
  for (let t = 1; t < n; t++) {
    p1[t] = u * c1[t]! + v * c2[t]!;
    p2[t] = a * c2[t]! - b * c1[t]!;
  }
  U[i] = p1;
  U[j] = p2;
}

function ZM_ishnf(x: ZMat): boolean {
  const n = lg(x) - 1;
  if (n === 0) return true;
  const m = nbrows(x);
  if (m !== n) return false;
  for (let i = 1; i <= n; i++) {
    if (x[i]![i]! <= 0n) return false;
    for (let j = 1; j < i; j++) if (x[j]![i]! !== 0n) return false;
    for (let j = i + 1; j <= n; j++) {
      const c = x[j]![i]!;
      if (c < 0n || c >= x[i]![i]!) return false;
    }
  }
  return true;
}

function transposeMat(M: ZMat): ZMat {
  const n = lg(M) - 1;
  const m = nbrows(M);
  const out: ZMat = new Array<ZC>(m + 1);
  out[0] = [];
  for (let i = 1; i <= m; i++) {
    const c = cgetg_int(n);
    for (let j = 1; j <= n; j++) c[j] = M[j]![i]!;
    out[i] = c;
  }
  return out;
}

/**
 * PARI `ZM_snfall_i(x, NULL, &V, 1|2)` (`hnf_snf.c:2386-2566`), restricted to
 * the case reached from `get_clgp`: `x` square, in HNF, of nonzero determinant.
 * Returns the diagonal of the SNF and `V` with `U x V = D`.
 */
function ZM_snfall_hnf(x0: ZMat): { D: ZC; V: ZMat } {
  const n = lg(x0) - 1;
  if (!n) return { D: cgetg_int(0), V: [[]] };
  const m0 = nbrows(x0);
  if (!ZM_ishnf(x0))
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: ZM_snfall for a matrix that is not a square HNF ' +
        '(needs ZM_hnfperm/ZM_hnfmod, hnf_snf.c:2400-2430)'
    );
  let mdet = ZM_det_triangular(x0);
  const V0 = matid(n);
  const A = x0;
  const x = shallowcopyMat(x0);
  let U: ZMat | null = matid(n);
  for (let i = n; i > 1; i--) {
    for (;;) {
      let c = 0;
      for (let j = i - 1; j >= 1; j--) {
        const b = x[j]![i]!;
        if (!b) continue;
        const a = x[i]![i]!;
        ZC_elem(b, a, x, null, j, i);
      }
      for (let j = i - 1; j >= 1; j--) {
        let b = x[i]![j]!;
        if (!b) continue;
        let a = x[i]![i]!;
        const st = bezout_step(a, b);
        const d = st.d;
        a = st.a;
        b = st.b;
        const u = st.u;
        const v = st.v;
        for (let k = 1; k < i; k++) {
          const t = u * x[k]![i]! + v * x[k]![j]!;
          x[k]![j] = a * x[k]![j]! - b * x[k]![i]!;
          x[k]![i] = t;
        }
        x[i]![j] = 0n;
        x[i]![i] = d;
        if (U) updateCols(u, v, a, b, U, i, j);
        c = 1;
      }
      if (!c) {
        const k = ZM_snf_no_divide(x, i);
        if (!k) break;
        for (let j = 1; j <= i; j++) x[j]![i] = x[j]![i]! + x[j]![k]!;
        if (U) {
          const ci = U[i]!;
          const ck = U[k]!;
          for (let t = 1; t < lg(ci); t++) ci[t] = ci[t]! + ck[t]!;
        }
      }
      ZM_redpart(x, mdet, i);
      if (U) ZM_redpart(U, mdet, n);
    }
  }
  for (let k = n; k; k--) {
    const d = gcdii(x[k]![k]!, mdet);
    x[k]![k] = d;
    if (d !== 1n) mdet = mdet / d;
  }
  U = transposeMat(U!);
  /* U A V = D => D^(-1) U A = V^(-1) */
  let W = ZM_mul(U, A);
  for (let i = 1; i <= n; i++) {
    const c = x[i]![i]!;
    if (c === 1n || c === -1n) break;
    for (let j = 1; j < lg(W); j++) {
      const col = W[j]!;
      if (col[i]! % c !== 0n) throw new Error('ZM_snfall: inexact row division');
      col[i] = col[i]! / c;
    }
  }
  const d1 = x[1]![1]!;
  const Wred: ZMat = new Array<ZC>(lg(W));
  Wred[0] = [];
  for (let j = 1; j < lg(W); j++) {
    const col = W[j]!;
    const c = cgetg_int(lg(col) - 1);
    for (let i = 1; i < lg(col); i++) c[i] = umodiu(col[i]!, d1);
    Wred[j] = c;
  }
  W = fromExternalZM(matinvmod(toExternalZM(Wred), d1));
  const V = ZM_mul(V0, W);
  /* return_vec: diagonal */
  const dg = cgetg_int(n);
  for (let i = 1; i <= n; i++) dg[i] = x[i]![i]!;
  let Dvec = dg;
  if (m0 > n) {
    const z = cgetg_int(m0 - n);
    Dvec = shallowconcat(z, dg);
  }
  return { D: Dvec, V };
}

/* conversion to/from the 0-based column-major layout of matkermod.ts */
function toExternalZM(M: ZMat): bigint[][] {
  const out: bigint[][] = [];
  for (let j = 1; j < lg(M); j++) out.push(M[j]!.slice(1));
  return out;
}
function fromExternalZM(M: bigint[][]): ZMat {
  const out: ZMat = new Array<ZC>(M.length + 1);
  out[0] = [];
  for (let j = 0; j < M.length; j++) out[j + 1] = [0n, ...M[j]!];
  return out;
}

/**
 * PARI `ZM_snf_group(H, NULL, &newUi)` (`hnf_snf.c:2893-2899` and
 * `snf_group`, `hnf_snf.c:2857-2889`).
 *
 * `H` is a relation matrix among the rows of generators `g`, in HNF.  Returns
 * the diagonal of the SNF with the `1`s removed and `Ui` such that the new
 * generators are `G = g Ui`.
 */
export function ZM_snf_group(H: ZMat): { D: ZC; Ui: ZMat } {
  const r = ZM_snfall_hnf(H);
  let D = r.D;
  let V = r.V;
  /* ZM_snfclean(D, NULL, V) */
  let c = 1;
  const l = lg(D);
  for (; c < l; c++) {
    const t = D[c]!;
    if (t === 1n || t === -1n) break;
  }
  D = vecslice(D, 1, c - 1);
  V = vecslice(V, 1, c - 1);
  /* snf_group: Ui = H * (V mod D), then divide column i by D[i], reduce mod H */
  const ln = lg(D);
  for (let i = 1; i < ln; i++) {
    const col = V[i]!;
    for (let k = 1; k < lg(col); k++) col[k] = umodiu(col[k]!, D[i]!);
  }
  const Ui = ZM_mul(H, V);
  for (let i = 1; i < ln; i++) {
    const col = Ui[i]!;
    for (let k = 1; k < lg(col); k++) {
      if (col[k]! % D[i]! !== 0n) throw new Error('snf_group: inexact division');
      col[k] = col[k]! / D[i]!;
    }
  }
  /* ZM_hnfrem(Ui, H): reduce columns of Ui modulo the HNF lattice H */
  const n = lg(H) - 1;
  for (let j = 1; j < lg(Ui); j++) {
    const col = Ui[j]!;
    for (let i = n; i >= 1; i--) {
      const h = H[i]![i]!;
      const q = diviiround(col[i]!, h);
      if (!q) continue;
      for (let k = 1; k < lg(col); k++) col[k] = col[k]! - q * H[i]![k]!;
    }
  }
  return { D, Ui };
}

/* ================================================================== */
/* Binary quadratic forms: qfr3 / qfr5 containers (Qfb.c:396-620)      */
/* ================================================================== */

/** PARI `struct qfr_data` (`paripriv.h`), built by `qfr_data_init` */
export interface QfrData {
  D: bigint;
  sqrtD: Real;
  isqrtD: bigint;
}

/** a `qfr3` container: `[a,b,c]` */
export type Qfr3 = [bigint, bigint, bigint];
/** a `qfr5` container: `[a,b,c,e,d]`, `e` a binary exponent, `d` a `t_REAL` */
export interface Qfr5 {
  a: bigint;
  b: bigint;
  c: bigint;
  e: bigint;
  d: Real;
}

/** PARI `qfr_data_init` (`Qfb.c:553-559`) */
export function qfr_data_init(D: bigint, prec: number): QfrData {
  const sqrtD = sqrtr(itor(D, prec));
  return { D, sqrtD, isqrtD: truncr(sqrtD) };
}

const EMAX = 22;
/** PARI `fix_expo` (`Qfb.c:411-418`) */
function fix_expo(x: Qfr5): void {
  if (real_expo(x.d) >= 1 << EMAX) {
    x.e = x.e + 1n;
    x.d = shiftr(x.d, -(1 << EMAX));
  }
}

/** PARI `qfr5_dist(e, d, prec)` (`Qfb.c:423-431`): `(1/2) log(|d| 2^(e 2^EMAX))` */
export function qfr5_dist(e: bigint, d: Real, prec: number): Real {
  let t = logr_abs(setprec(d, prec));
  if (e) {
    const u = shiftr(mulir(e, mplog2(prec)), EMAX);
    t = addrr(t, u);
  }
  return shiftr(t, -1);
}

/** PARI `rho_get_BC` (`Qfb.c:433-441`) */
function rho_get_BC(a: bigint, b: bigint, c: bigint, S: QfrData): [bigint, bigint] {
  const t = iabs(S.isqrtD) >= iabs(c) ? S.isqrtD : iabs(c);
  const [q, u] = truedvmdii(t + b, c << 1n);
  const B = t - u;
  const C = a - q * (b - q * c);
  return [B, C];
}

/** PARI `qfr3_rho` (`Qfb.c:444-450`) */
export function qfr3_rho(x: Qfr3, S: QfrData): Qfr3 {
  const [B, C] = rho_get_BC(x[0], x[1], x[2], S);
  return [x[2], B, C];
}

/** PARI `qfr5_rho` (`Qfb.c:453-471`) */
export function qfr5_rho(x: Qfr5, S: QfrData): Qfr5 {
  const sb = signe(x.b);
  const [B, C] = rho_get_BC(x.a, x.b, x.c, S);
  const y: Qfr5 = { a: x.c, b: B, c: C, e: x.e, d: x.d };
  if (sb) {
    let t: Real;
    const num = x.b * x.b - S.D;
    if (sb < 0) {
      const den = subrr(itor(x.b, S.sqrtD.p), S.sqrtD);
      t = divir(num, mulrr(den, den));
    } else {
      const s = addrr(itor(x.b, S.sqrtD.p), S.sqrtD);
      t = divri(mulrr(s, s), num);
    }
    y.d = mulrr(t, y.d);
    fix_expo(y);
  } else {
    y.d = real_neg(y.d);
  }
  return y;
}

/** PARI `qfr_to_qfr5` (`Qfb.c:473-475`) */
function qfr_to_qfr5(x: Qfr3, prec: number): Qfr5 {
  return { a: x[0], b: x[1], c: x[2], e: 0n, d: real_1(prec) };
}

/** PARI `ab_isreduced` (`Qfb.c:505-513`) */
function ab_isreduced(a: bigint, b: bigint, isqrtD: bigint): boolean {
  if (signe(b) <= 0 || iabs(b) > iabs(isqrtD)) return false;
  const t = isqrtD - iabs(2n * a);
  return signe(t) < 0 ? iabs(b) >= iabs(t) : iabs(b) > iabs(t);
}

/** PARI `qfr5_red` (`Qfb.c:518-531`) */
export function qfr5_red(x: Qfr5, S: QfrData): Qfr5 {
  while (!ab_isreduced(x.a, x.b, S.isqrtD)) x = qfr5_rho(x, S);
  return x;
}

/** PARI `qfr3_red` (`Qfb.c:534-550`) */
export function qfr3_red(x: Qfr3, S: QfrData): Qfr3 {
  let [a, b, c] = x;
  while (!ab_isreduced(a, b, S.isqrtD)) {
    const [B, C] = rho_get_BC(a, b, c, S);
    a = c;
    b = B;
    c = C;
  }
  return [a, b, c];
}

/** PARI `qfb_sqr` (`Qfb.c:1013-1036`), on a 3-component container */
function qfb_sqr3(x: Qfr3): Qfr3 {
  const [d1, x2] = bezout(x[1], x[0]);
  let c = x[2];
  let m = c * x2;
  let v1: bigint;
  let v2: bigint;
  if (d1 === 1n) {
    v1 = v2 = x[0];
  } else {
    v1 = x[0] / d1;
    v2 = v1 * gcdii(d1, c);
    c = c * d1;
  }
  m = -m;
  const r = ((m % v2) + v2) % v2;
  const p1 = r * v1;
  const c3 = c + r * (x[1] + p1);
  return [v1 * v2, x[1] + (p1 << 1n), c3 / v2];
}

/** PARI `qfb_comp` (`Qfb.c:1038-1071`), on 3-component containers */
function qfb_comp3(x: Qfr3, y: Qfr3): Qfr3 {
  const n = (y[1] - x[1]) >> 1n;
  let v1 = x[0];
  let v2 = y[0];
  let c = y[2];
  const [d, y1] = bezout(v2, v1);
  let m: bigint;
  if (d === 1n) {
    m = y1 * n;
  } else {
    const s = y[1] - n;
    const [d1, x2, y2] = bezout(s, d);
    if (d1 !== 1n) {
      v1 = v1 / d1;
      v2 = v2 / d1;
      v1 = v1 * gcdii(c, gcdii(x[2], gcdii(d1, n)));
      c = c * d1;
    }
    m = y1 * y2 * n + y[2] * x2;
  }
  m = -m;
  const r = ((m % v1) + v1) % v1;
  const p1 = r * v2;
  const c3 = c + r * (y[1] + p1);
  return [v1 * v2, y[1] + (p1 << 1n), c3 / v1];
}

/** PARI `qfr5_compraw` (`Qfb.c:1477-1493`) */
function qfr5_compraw(x: Qfr5, y: Qfr5): Qfr5 {
  const same = x === y;
  const z3 = same ? qfb_sqr3([x.a, x.b, x.c]) : qfb_comp3([x.a, x.b, x.c], [y.a, y.b, y.c]);
  const z: Qfr5 = {
    a: z3[0],
    b: z3[1],
    c: z3[2],
    e: same ? x.e << 1n : x.e + y.e,
    d: same ? sqrr(x.d) : mulrr(x.d, y.d),
  };
  fix_expo(z);
  return z;
}
/** PARI `qfr5_comp` (`Qfb.c:1494-1495`) */
export function qfr5_comp(x: Qfr5, y: Qfr5, S: QfrData): Qfr5 {
  return qfr5_red(qfr5_compraw(x, y), S);
}
/** PARI `qfr3_compraw` (`Qfb.c:1497-1502`) */
function qfr3_compraw(x: Qfr3, y: Qfr3): Qfr3 {
  return x === y ? qfb_sqr3(x) : qfb_comp3(x, y);
}
/** PARI `qfr3_comp` (`Qfb.c:1503-1505`) */
export function qfr3_comp(x: Qfr3, y: Qfr3, S: QfrData): Qfr3 {
  return qfr3_red(qfr3_compraw(x, y), S);
}
/** PARI `qfr3_powraw` (`Qfb.c:1542-1553`), `m > 0` */
function qfr3_powraw(x: Qfr3, m: bigint): Qfr3 {
  let y: Qfr3 | null = null;
  let b = x;
  for (; m; m >>= 1n) {
    if (m & 1n) y = y ? qfr3_compraw(y, b) : b;
    if (m === 1n) break;
    b = qfb_sqr3(b);
  }
  return y!;
}
/** PARI `qfr3_pow` (`Qfb.c:1555-1580`) */
export function qfr3_pow(x: Qfr3, n: bigint, S: QfrData): Qfr3 {
  const s = signe(n);
  if (!s) return qfr3_1(S);
  let base = x;
  if (s < 0) base = [x[0], -x[1], x[2]];
  let m = s < 0 ? -n : n;
  let y: Qfr3 | null = null;
  for (; m; m >>= 1n) {
    if (m & 1n) y = y ? qfr3_comp(y, base, S) : base;
    if (m === 1n) break;
    base = qfr3_comp(base, base, S);
  }
  return y!;
}

/** PARI `qfr_1_fill` (`Qfb.c:1244-1253`) */
function qfr_1_fill(S: QfrData): Qfr3 {
  let y2 = S.isqrtD;
  if (mod2(S.D) !== mod2(y2)) y2 = y2 - 1n;
  return [1n, y2, (y2 * y2 - S.D) >> 2n];
}
/** PARI `qfr3_1` (`Qfb.c:1262-1267`) */
function qfr3_1(S: QfrData): Qfr3 {
  return qfr_1_fill(S);
}
/** PARI `qfr5_1` (`Qfb.c:1254-1261`) */
function qfr5_1(S: QfrData, prec: number): Qfr5 {
  const y = qfr_1_fill(S);
  return { a: y[0], b: y[1], c: y[2], e: 0n, d: real_1(prec) };
}

/**
 * PARI `primeform_u(x, p)` (`Qfb.c:1661-1694`): the prime form of
 * discriminant `x` above `p`, as a 3-component container.
 */
export function primeform_u3(x: bigint, p: bigint): Qfr3 {
  let s = mod8(x);
  if (signe(x) < 0 && s) s = 8 - s;
  if (s & 2) throw new PariDomainError('primeform', 'disc % 4', '>', '1');
  let b: bigint;
  let c: bigint;
  if (p === 2n) {
    switch (s) {
      case 0:
        b = 0n;
        break;
      case 1:
        b = 1n;
        break;
      case 4:
        b = 2n;
        break;
      default:
        throw new PariSqrtnError('primeform', `Mod(${x}, ${p})`);
    }
    c = (BigInt(s) - x) >> 3n;
  } else {
    const r = Fp_sqrt(umodiu(x, p), p);
    if (r === null) throw new PariSqrtnError('primeform', `Mod(${x}, ${p})`);
    b = r;
    if ((mod2(b) ^ (s & 1)) & 1) b = p - b;
    c = ((b * b - x) >> 2n) / p;
  }
  return [p, b, c];
}

/* ================================================================== */
/* GRH check (buch2.c:354-386, buch1.c:318-427)                        */
/* ================================================================== */

interface GRHprime {
  p: number;
  logp: number;
  /** kronecker(D, p) */
  dec: number;
}
interface GRHcheck {
  cN: number;
  cD: number;
  primes: GRHprime[] /* 0-based, in increasing order */;
  limp: number;
}

/** PARI `init_GRHcheck` (`buch2.c:353-367`) */
function init_GRHcheck(N: number, R1: number, LOGD: number): GRHcheck {
  const c1 = (Math.PI * Math.PI) / 2;
  const c2 = 3.663862376709;
  const c3 = 3.801387092431; /* Euler + log(8*Pi) */
  return {
    cN: R1 * c2 + N * c1,
    cD: LOGD - N * c3 - (R1 * Math.PI) / 2,
    primes: [],
    limp: 0,
  };
}

/** PARI `GRHok` (`buch2.c:381-386`) */
function GRHok(S: GRHcheck, L: number, SA: number, SB: number): boolean {
  return S.cD + (S.cN + 2 * SB) / L - 2 * SA < -1e-8;
}

/** simple prime iterator (PARI `u_forprime_next`) */
function nextprime_u(n: number): number {
  let p = n < 2 ? 2 : n + 1;
  for (;;) {
    if (isPrime(BigInt(p))) return p;
    p++;
  }
}

/** PARI `kroiu(D, p)`: Kronecker symbol `(D/p)` for `p` prime */
function kroiu(D: bigint, p: number): number {
  return kronecker(D, BigInt(p));
}

/** PARI `cache_prime_quad` (`buch1.c:318-338`) */
function cache_prime_quad(S: GRHcheck, LIM: number, D: bigint): void {
  if (S.limp >= LIM) return;
  let p = S.primes.length ? S.primes[S.primes.length - 1]!.p : 1;
  for (;;) {
    p = nextprime_u(p);
    S.primes.push({ p, logp: Math.log(p), dec: kroiu(D, p) });
    /* store up to nextprime(LIM) included */
    if (p >= LIM) {
      S.limp = p;
      break;
    }
  }
}

/** PARI `compute_invresquad` (`buch1.c:340-363`) */
function compute_invresquad(S: GRHcheck, LIMC: number): Real {
  let invres = real_1(DEFAULTPREC);
  const limp = Math.log(LIMC) / 2;
  for (const pr of S.primes) {
    const s = pr.dec;
    if (s) {
      const p = pr.p;
      if (s > 0 || pr.logp <= limp) invres = mulur(p - s, divru(invres, p));
      else if (s < 0) invres = mulur(p, divru(invres, p - 1));
    }
  }
  return invres;
}

/** PARI `is_bad` (`buch1.c:365-377`): `p | conductor` of the order of disc `D`? */
function is_bad(D: bigint, p: number): boolean {
  if (p === 2) {
    let r = mod16(D) >> 1;
    if (r && signe(D) < 0) r = 8 - r;
    return r < 4;
  }
  const P = BigInt(p);
  return D % (P * P) === 0n;
}

/** PARI `nthidealquad` (`buch1.c:379-390`) */
function nthidealquad(D: bigint, n: number): number {
  let p = 1;
  for (;;) {
    p = nextprime_u(p);
    if (n <= 0) return p;
    if (!is_bad(D, p) && kroiu(D, p) >= 0) n--;
  }
}

/** PARI `quadGRHchk` (`buch1.c:392-427`) */
function quadGRHchk(D: bigint, S: GRHcheck, LIMC: number): boolean {
  const logC = Math.log(LIMC);
  let SA = 0;
  let SB = 0;
  cache_prime_quad(S, LIMC, D);
  for (let i = 0; ; i++) {
    const pr = S.primes[i]!;
    const p = pr.p;
    if (p > LIMC) break;
    let logNP: number;
    let q: number;
    if (pr.dec < 0) {
      logNP = 2 * pr.logp;
      q = 1 / p;
    } else {
      logNP = pr.logp;
      q = 1 / Math.sqrt(p);
    }
    let A = logNP * q;
    let B = logNP * A;
    const M = Math.floor(logC / logNP);
    if (M > 1) {
      const inv1_q = 1 / (1 - q);
      A *= (1 - q ** M) * inv1_q;
      B *= (1 - q ** M * (M + 1 - M * q)) * inv1_q * inv1_q;
    }
    if (pr.dec > 0) {
      SA += 2 * A;
      SB += 2 * B;
    } else {
      SA += A;
      SB += B;
    }
    if (p === LIMC) break;
  }
  return GRHok(S, logC, SA, SB);
}

/**
 * PARI `bnf_increase_LIMC` (`buch1.c:189-197`): suggest a larger factor base
 * bound after a failure.
 */
export function bnf_increase_LIMC(LIMC: number, D: number): number {
  if (LIMC <= D / 13.333) LIMC *= 2;
  else LIMC += Math.max(1, Math.floor(D / 20));
  return Math.floor(LIMC);
}

/* ================================================================== */
/* Buchquad (buch1.c)                                                  */
/* ================================================================== */

/**
 * Deterministic pseudo-random source replacing PARI's `pari_rand`
 * (`random_bits`, `buch1.c:165`).  The algorithm certifies its own output
 * (`get_R`, `quad_be_honest`, the HNF rank test), so the stream only affects
 * how long the relation search takes, never the result.
 */
let RAND_STATE = 0x2545f491;
export function setBuchRandomSeed(seed: number): void {
  RAND_STATE = seed >>> 0 || 1;
}
function random_bits(k: number): number {
  /* xorshift32 */
  let x = RAND_STATE;
  x ^= x << 13;
  x >>>= 0;
  x ^= x >>> 17;
  x ^= x << 5;
  x >>>= 0;
  RAND_STATE = x;
  return x >>> (32 - k);
}

const RANDOM_BITS = 4;
const CBUCH = (1 << RANDOM_BITS) - 1;

/** entry of the large prime hash table (`buch1.c:265-288`) */
interface LPEntry {
  nrho: number;
  np: number;
  q: number;
  ex: number[];
}

interface BuchQuad {
  limhash: number;
  KC: number;
  KC2: number;
  PRECREG: number;
  primfact: number[];
  exprimfact: number[];
  hashtab: Map<number, LPEntry>;
  FB: zv;
  numFB: number[];
  prodFB: bigint;
  powsubFB: (Qfb[] | Qfr5[])[];
  vperm: zv;
  subFB: zv;
  badprim: bigint | null;
  q: QfrData;
}

/** PARI `Z_isquasismooth_prod` (`buch1.c:206-216`) */
function Z_isquasismooth_prod(N: bigint, P: bigint): bigint {
  P = gcdii(P, N);
  while (P !== 1n) {
    N = N / P;
    P = gcdii(N, P);
  }
  return N;
}

/**
 * PARI `factorquad` (`buch1.c:218-262`): try to factor the first coefficient of
 * the form `f` over the factor base.  Returns 0 on failure, 1 when fully
 * factored, and the (possibly composite) unfactored part otherwise.
 */
function factorquad(B: BuchQuad, a: bigint, nFB: number, limp: number): number {
  const P = B.primfact;
  const E = B.exprimfact;
  let lo = 0;
  let x = iabs(a);
  if (B.badprim && gcdii(x, B.badprim) !== 1n) return 0;
  const F = Z_isquasismooth_prod(x, B.prodFB);
  if (F > BigInt(B.limhash)) return 0;
  for (let i = 1; ; i++) {
    const p = BigInt(B.FB[i]!);
    let q = x / p;
    let r = x - q * p;
    if (!r) {
      let k = 0;
      do {
        k++;
        x = q;
        q = x / p;
        r = x - q * p;
      } while (!r);
      lo++;
      P[lo] = B.FB[i]!;
      E[lo] = k;
    }
    if (q <= p) break;
    if (i === nFB) return 0;
  }
  if (x === 1n) {
    P[0] = lo;
    return 1;
  }
  if (x > BigInt(B.limhash)) return 0;
  const X = Number(x);
  if (X !== 1 && X <= limp) {
    lo++;
    P[lo] = X;
    E[lo] = 1;
    P[0] = lo;
    return 1;
  }
  P[0] = lo;
  return X;
}

/** PARI `largeprime` (`buch1.c:264-288`) */
function largeprime(B: BuchQuad, q: number, ex: zv, np: number, nrho: number): LPEntry | null {
  const l = lg(B.subFB);
  const pt = B.hashtab.get(q);
  if (!pt) {
    const e: number[] = new Array<number>(l).fill(0);
    for (let i = 1; i < l; i++) e[i] = ex[i]!;
    B.hashtab.set(q, { nrho, np, q, ex: e });
    return null;
  }
  for (let i = 1; i < l; i++) if (pt.ex[i] !== ex[i]) return pt;
  return pt.np === np ? null : pt;
}

/** PARI `FBquad` (`buch1.c:430-474`) */
function FBquad(B: BuchQuad, C2: number, C1: number, S: GRHcheck): void {
  const D = B.q.D;
  cache_prime_quad(S, C2, D);
  const numFB: number[] = new Array<number>(C2 + 1).fill(0);
  const FB: zv = [0];
  let i = 0;
  B.KC = 0;
  let badprim = 1n;
  for (let idx = 0; ; idx++) {
    const pr = S.primes[idx]!;
    const p = pr.p;
    if (!B.KC && p > C1) B.KC = i;
    if (p > C2) break;
    /* upstream switch: case -1 inert (skip); case 0 ramified, skipped when
     * `is_bad`, otherwise falls through to the split case */
    if (pr.dec !== -1) {
      if (pr.dec === 0 && is_bad(D, p)) badprim = badprim * BigInt(p);
      else {
        i++;
        numFB[p] = i;
        FB[i] = p;
      }
    }
    if (p === C2) {
      if (!B.KC) B.KC = i;
      break;
    }
  }
  B.KC2 = i;
  B.FB = FB;
  B.numFB = numFB;
  B.badprim = badprim === 1n ? null : badprim;
  let prod = 1n;
  for (let j = 1; j <= B.KC2; j++) prod *= BigInt(FB[j]!);
  B.prodFB = prod;
}

/** PARI `subFBquad` (`buch1.c:476-505`) */
function subFBquad(B: BuchQuad, D: bigint, PROD: number, minSFB: number): zv {
  const lv = B.KC + 1;
  let lgsub = 1;
  let ino = 1;
  let prod = 1;
  const vperm = cgetg_zv(lv - 1);
  const no = cgetg_zv(lv - 1);
  let j = 1;
  for (; j < lv; j++) {
    const p = B.FB[j]!;
    if (umodiu(D, BigInt(p)) === 0n) no[ino++] = j;
    else {
      vperm[lgsub++] = j;
      prod *= p;
      if (lgsub > minSFB && prod > PROD) break;
    }
  }
  let i = lgsub;
  for (j = 1; j < ino; i++, j++) vperm[i] = no[j]!;
  for (; i < lv; i++) vperm[i] = i;
  B.vperm = vperm;
  return vecslice(vperm, 1, lgsub - 1);
}

/** PARI `powsubFBquad` (`buch1.c:507-536`) */
function powsubFBquad(B: BuchQuad, n: number): (Qfb[] | Qfr5[])[] {
  const l = lg(B.subFB);
  const x: (Qfb[] | Qfr5[])[] = new Array(l);
  const D = B.q.D;
  if (B.PRECREG) {
    for (let i = 1; i < l; i++) {
      const F = qfr5_pf(B.q, B.FB[B.subFB[i]!]!, B.PRECREG);
      const y: Qfr5[] = new Array(n + 1);
      y[1] = F;
      for (let j = 2; j <= n; j++) y[j] = QFR5_comp(y[j - 1]!, F, B.q);
      x[i] = y;
    }
  } else {
    for (let i = 1; i < l; i++) {
      const F = qfi_pf(D, B.FB[B.subFB[i]!]!);
      const y: Qfb[] = new Array(n + 1);
      y[1] = F;
      for (let j = 2; j <= n; j++) y[j] = qfbcomp(y[j - 1]!, F);
      x[i] = y;
    }
  }
  return x;
}

/* --- the qfr3/qfr5 wrappers of buch1.c:58-173 --- */

/** PARI `qfr3_canon` (`buch1.c:58-68`) */
function qfr3_canon(x: Qfr3, S: QfrData): Qfr3 {
  const [a, , c] = x;
  if (signe(a) < 0) {
    if (iabs(a) === iabs(c)) return qfr3_rho(x, S);
    return [-a, x[1], -c];
  }
  return x;
}
/** PARI `qfr5_canon` (`buch1.c:80-90`) */
function qfr5_canon(x: Qfr5, S: QfrData): Qfr5 {
  if (signe(x.a) < 0) {
    if (iabs(x.a) === iabs(x.c)) return qfr5_rho(x, S);
    return { a: -x.a, b: x.b, c: -x.c, e: x.e, d: x.d };
  }
  return x;
}
/** PARI `QFR5_comp` (`buch1.c:91-93`) */
function QFR5_comp(x: Qfr5, y: Qfr5, S: QfrData): Qfr5 {
  return qfr5_canon(qfr5_comp(x, y, S), S);
}
/** PARI `QFR3_comp` (`buch1.c:94-96`) */
function QFR3_comp(x: Qfr3, y: Qfr3, S: QfrData): Qfr3 {
  return qfr3_canon(qfr3_comp(x, y, S), S);
}
/** PARI `qfr5_rho_pow` (`buch1.c:99-114`) */
function qfr5_rho_pow(x: Qfr5, n: number, S: QfrData): Qfr5 {
  for (let i = 1; i <= n; i++) x = qfr5_rho(x, S);
  return x;
}
/** PARI `qfr5_pf` (`buch1.c:116-121`) */
function qfr5_pf(S: QfrData, p: number, prec: number): Qfr5 {
  const y = primeform_u3(S.D, BigInt(p));
  return qfr5_canon(qfr5_red(qfr_to_qfr5(y, prec), S), S);
}
/** PARI `qfr3_pf` (`buch1.c:123-128`) */
function qfr3_pf(S: QfrData, p: number): Qfr3 {
  const y = primeform_u3(S.D, BigInt(p));
  return qfr3_canon(qfr3_red(y, S), S);
}
/** PARI `qfi_pf` = `primeform_u` (`buch1.c:130`) */
function qfi_pf(D: bigint, p: number): Qfb {
  const y = primeform_u3(D, BigInt(p));
  return mkqfb(y[0], y[1], y[2], D);
}

/** PARI `init_form` (`buch1.c:133-146`) for qfr5 */
function qfr5_factorback(B: BuchQuad, ex: readonly number[]): Qfr5 | null {
  const l = lg(B.powsubFB);
  let F: Qfr5 | null = null;
  for (let i = 1; i < l; i++)
    if (ex[i]) {
      const t = (B.powsubFB[i] as Qfr5[])[ex[i]!]!;
      F = F ? QFR5_comp(F, t, B.q) : t;
    }
  return F;
}
/** PARI `init_form` for qfi (`buch1.c:150-154`) */
function qfi_factorback(B: BuchQuad, ex: readonly number[]): Qfb | null {
  const l = lg(B.powsubFB);
  let F: Qfb | null = null;
  for (let i = 1; i < l; i++)
    if (ex[i]) {
      const t = (B.powsubFB[i] as Qfb[])[ex[i]!]!;
      F = F ? qfbcomp(F, t) : t;
    }
  return F;
}
/** PARI `random_form` (`buch1.c:156-173`) */
function qfr3_random(B: BuchQuad, ex: zv): Qfr3 {
  const l = lg(ex);
  for (;;) {
    for (let i = 1; i < l; i++) ex[i] = random_bits(RANDOM_BITS);
    const F = qfr5_factorback_to3(B, ex);
    if (F) return F;
  }
}
function qfr5_factorback_to3(B: BuchQuad, ex: zv): Qfr3 | null {
  /* init_form with QFR3_comp: the qfr3 analogue of qfr5_factorback */
  const l = lg(B.powsubFB);
  let F: Qfr3 | null = null;
  for (let i = 1; i < l; i++)
    if (ex[i]) {
      const t5 = (B.powsubFB[i] as Qfr5[])[ex[i]!]!;
      const t: Qfr3 = [t5.a, t5.b, t5.c];
      F = F ? QFR3_comp(F, t, B.q) : t;
    }
  return F;
}
function qfi_random(B: BuchQuad, ex: zv): Qfb {
  const l = lg(ex);
  for (;;) {
    for (let i = 1; i < l; i++) ex[i] = random_bits(RANDOM_BITS);
    const F = qfi_factorback(B, ex);
    if (F) return F;
  }
}

/** PARI `sub_fact` (`buch1.c:538-550`) */
function sub_fact(B: BuchQuad, col: zv, b: bigint): void {
  for (let i = 1; i <= B.primfact[0]!; i++) {
    const p = B.primfact[i]!;
    const k = B.numFB[p]!;
    let e = B.exprimfact[i]!;
    if (umodiu(b, BigInt(p << 1)) > BigInt(p)) e = -e;
    col[k] = col[k]! - e;
  }
}
/** PARI `add_fact` (`buch1.c:580-592`) */
function add_fact(B: BuchQuad, col: zv, b: bigint): void {
  for (let i = 1; i <= B.primfact[0]!; i++) {
    const p = B.primfact[i]!;
    const k = B.numFB[p]!;
    let e = B.exprimfact[i]!;
    if (umodiu(b, BigInt(p << 1)) > BigInt(p)) e = -e;
    col[k] = col[k]! + e;
  }
}

/** PARI `get_clgp` (`buch1.c:594-629`) */
function get_clgp(B: BuchQuad, W: ZMat): { cyc: ZC; gen: Qfb[] } {
  const { D, Ui } = ZM_snf_group(W);
  const l = lg(W);
  const c = lg(D);
  const res: Qfb[] = new Array<Qfb>(c);
  const init: Qfb[] = new Array<Qfb>(l);
  for (let i = 1; i < l; i++) init[i] = qfi_pf(B.q.D, B.FB[B.vperm[i]!]!);
  for (let j = 1; j < c; j++) {
    let g: Qfb | null = null;
    if (signe(B.q.D) > 0) {
      let g3: Qfr3 | null = null;
      for (let i = 1; i < l; i++) {
        const u = Ui[j]![i]!;
        if (!u) continue;
        const f = init[i]!;
        const t = qfr3_pow([f.a, f.b, f.c], u, B.q);
        g3 = g3 ? qfr3_comp(g3, t, B.q) : t;
      }
      const r = qfr3_canon_safe(qfr3_red(g3!, B.q), B.q);
      g = mkqfb(r[0], r[1], r[2], B.q.D);
    } else {
      for (let i = 1; i < l; i++) {
        const u = Ui[j]![i]!;
        if (!u) continue;
        const t = qfbpow(init[i]!, u);
        g = g ? qfbcomp(g, t) : t;
      }
    }
    res[j] = g!;
  }
  return { cyc: D, gen: res };
}
/** PARI `qfr3_canon_safe` (`buch1.c:69-79`) */
function qfr3_canon_safe(x: Qfr3, S: QfrData): Qfr3 {
  const [a, b, c] = x;
  if (signe(a) < 0) {
    if (iabs(a) === iabs(c)) return qfr3_rho(x, S);
    return [-a, b, -c];
  }
  return x;
}

/** PARI `trivial_relations` (`buch1.c:631-645`) */
function trivial_relations(B: BuchQuad, mat: zv[], C: CVec, prec: number): number {
  let j = 0;
  const D = B.q.D;
  for (let i = 1; i <= B.KC; i++) {
    if (umodiu(D, BigInt(B.FB[i]!)) !== 0n) continue;
    const col = zero_zv(B.KC);
    col[i] = 2;
    j++;
    mat[j] = col;
    C[j] = creal(real_0(prec || DEFAULTPREC), 0n);
  }
  return j;
}

/** PARI `imag_relations` (`buch1.c:668-736`) */
function imag_relations(
  B: BuchQuad,
  need: number,
  pc: { v: number },
  LIMC: number,
  mat: zv[],
  off: number
): void {
  const lgsub = lg(B.subFB);
  let current = pc.v;
  let nbtest = 0;
  let s = 0;
  const ex = cgetg_zv(lgsub - 1);
  if (!current) current = 1;
  for (;;) {
    if (s >= need) break;
    let form = qfi_random(B, ex);
    form = qfbcomp(form, qfi_pf(B.q.D, B.FB[current]!));
    nbtest++;
    const fpc = factorquad(B, form.a, B.KC, LIMC);
    if (!fpc) {
      if ((nbtest & 0xff) === 0 && ++current > B.KC) current = 1;
      continue;
    }
    let col: zv;
    if (fpc > 1) {
      const fpd = largeprime(B, fpc, ex, current, 0);
      if (!fpd) continue;
      let form2 = qfi_factorback(B, fpd.ex)!;
      form2 = qfbcomp(form2, qfi_pf(B.q.D, B.FB[fpd.np]!));
      const p = BigInt(fpc << 1);
      const b1 = umodiu(form2.b, p);
      const b2 = umodiu(form.b, p);
      if (b1 !== b2 && b1 + b2 !== p) continue;
      col = mat[off + ++s]!;
      add_fact(B, col, form.b);
      factorquad(B, form2.a, B.KC, LIMC);
      if (b1 === b2) {
        for (let i = 1; i < lgsub; i++) col[B.subFB[i]!] = col[B.subFB[i]!]! + fpd.ex[i]! - ex[i]!;
        sub_fact(B, col, form2.b);
        col[fpd.np] = col[fpd.np]! + 1;
      } else {
        for (let i = 1; i < lgsub; i++) col[B.subFB[i]!] = col[B.subFB[i]!]! - fpd.ex[i]! - ex[i]!;
        add_fact(B, col, form2.b);
        col[fpd.np] = col[fpd.np]! - 1;
      }
    } else {
      col = mat[off + ++s]!;
      for (let i = 1; i < lgsub; i++) col[B.subFB[i]!] = -ex[i]!;
      add_fact(B, col, form.b);
    }
    col[current] = col[current]! - 1;
    if (++current > B.KC) current = 1;
  }
  pc.v = current;
}

/** PARI `imag_be_honest` (`buch1.c:757-775`) */
function imag_be_honest(B: BuchQuad): boolean {
  let s = B.KC;
  let nbtest = 0;
  const ex = cgetg_zv(lg(B.subFB) - 1);
  while (s < B.KC2) {
    const p = B.FB[s + 1]!;
    const F = qfbcomp(qfi_pf(B.q.D, p), qfi_random(B, ex));
    const fpc = factorquad(B, F.a, s, p - 1);
    if (fpc === 1) {
      nbtest = 0;
      s++;
    } else if (++nbtest > 40) return false;
  }
  return true;
}

/** PARI `dist` (`buch1.c:777-778`) */
function dist(e: bigint, d: Real, prec: number): { t: Real; d: Real } {
  return { t: qfr5_dist(e, d, prec), d };
}

/** PARI `real_relations` (`buch1.c:789-945`) */
function real_relations(
  B: BuchQuad,
  need: number,
  pc: { v: number },
  lim: number,
  LIMC: number,
  mat: zv[],
  C: CVec,
  off: number
): void {
  const lgsub = lg(B.subFB);
  const prec = B.PRECREG;
  let current = pc.v;
  let s = 0;
  let first = current === 0;
  const ex = cgetg_zv(lgsub - 1);
  if (!current) current = 1;
  if (lim > need) lim = need;

  OUTER: for (;;) {
    if (s >= need) break;
    if (first && s >= lim) first = false;
    let form: Qfr3 = qfr3_random(B, ex);
    if (!first) form = QFR3_comp(form, qfr3_pf(B.q, B.FB[current]!), B.q);
    const form0 = form;
    let form1: Qfr5 | null = null;
    let endcycle = false;
    let rhoacc = 0;
    let rho = -1;

    /* CYCLE */
    for (;;) {
      if (endcycle || rho > 5000) {
        if (++current > B.KC) current = 1;
        continue OUTER;
      }
      if (rho < 0) rho = 0;
      else {
        form = qfr3_rho(form, B.q);
        rho++;
        rhoacc++;
        if (first) {
          endcycle = iabs(form[0]) === iabs(form0[0]) && form[1] === form0[1];
        } else {
          if (iabs(form[0]) === iabs(form[2])) {
            if (iabs(form[0]) === iabs(form0[0]) && form[1] === form0[1]) continue;
            form = qfr3_rho(form, B.q);
            rho++;
            rhoacc++;
          } else {
            form = [iabs(form[0]), form[1], -iabs(form[2])];
          }
          if (form[0] === form0[0] && form[1] === form0[1]) continue;
        }
      }
      const fpc = factorquad(B, form[0], B.KC, LIMC);
      if (!fpc) continue;
      let col: zv;
      let d: { t: Real; d: Real };
      if (fpc > 1) {
        const fpd = largeprime(B, fpc, ex, first ? 0 : current, rhoacc);
        if (!fpd) continue;
        if (!form1) {
          form1 = qfr5_factorback(B, ex)!;
          if (!first) form1 = QFR5_comp(form1, qfr5_pf(B.q, B.FB[current]!, prec), B.q);
        }
        form1 = qfr5_rho_pow(form1, rho, B.q);
        rho = 0;

        let form2 = qfr5_factorback(B, fpd.ex)!;
        if (fpd.np) form2 = QFR5_comp(form2, qfr5_pf(B.q, B.FB[fpd.np]!, prec), B.q);
        form2 = qfr5_rho_pow(form2, fpd.nrho, B.q);
        if (iabs(form2.a) !== iabs(form2.c))
          form2 = { a: iabs(form2.a), b: form2.b, c: -iabs(form2.c), e: form2.e, d: form2.d };
        const p = BigInt(fpc << 1);
        const b1 = umodiu(form2.b, p);
        const b2 = umodiu(form1.b, p);
        if (b1 !== b2 && b1 + b2 !== p) continue;

        col = mat[off + ++s]!;
        add_fact(B, col, form1.b);
        factorquad(B, form2.a, B.KC, LIMC);
        if (b1 === b2) {
          for (let i = 1; i < lgsub; i++)
            col[B.subFB[i]!] = col[B.subFB[i]!]! + fpd.ex[i]! - ex[i]!;
          sub_fact(B, col, form2.b);
          if (fpd.np) col[fpd.np] = col[fpd.np]! + 1;
          d = dist(form1.e - form2.e, divrr(form1.d, form2.d), prec);
        } else {
          for (let i = 1; i < lgsub; i++)
            col[B.subFB[i]!] = col[B.subFB[i]!]! - fpd.ex[i]! - ex[i]!;
          add_fact(B, col, form2.b);
          if (fpd.np) col[fpd.np] = col[fpd.np]! - 1;
          d = dist(form1.e + form2.e, mulrr(form1.d, form2.d), prec);
        }
      } else {
        if (!form1) {
          form1 = qfr5_factorback(B, ex)!;
          if (!first) form1 = QFR5_comp(form1, qfr5_pf(B.q, B.FB[current]!, prec), B.q);
        }
        form1 = qfr5_rho_pow(form1, rho, B.q);
        rho = 0;
        col = mat[off + ++s]!;
        for (let i = 1; i < lgsub; i++) col[B.subFB[i]!] = -ex[i]!;
        add_fact(B, col, form1.b);
        d = dist(form1.e, form1.d, prec);
      }
      C[off + s] = creal(d.t, real_sign(d.d) < 0 ? 1n : 0n);
      if (first) {
        if (s >= need) break OUTER;
        if (s >= lim) continue OUTER;
        continue;
      }
      col[current] = col[current]! - 1;
      if (++current > B.KC) current = 1;
      continue OUTER;
    }
  }
  pc.v = current;
}

/** PARI `real_be_honest` (`buch1.c:947-970`) */
function real_be_honest(B: BuchQuad): boolean {
  let s = B.KC;
  let nbtest = 0;
  const ex = cgetg_zv(lg(B.subFB) - 1);
  while (s < B.KC2) {
    const p = B.FB[s + 1]!;
    let F = QFR3_comp(qfr3_random(B, ex), qfr3_pf(B.q, p), B.q);
    const F0 = F;
    for (;;) {
      const fpc = factorquad(B, F[0], s, p - 1);
      if (fpc === 1) {
        nbtest = 0;
        s++;
        break;
      }
      if (++nbtest > 40) return false;
      F = qfr3_canon(qfr3_rho(F, B.q), B.q);
      if (F[0] === F0[0] && F[1] === F0[1]) break;
    }
  }
  return true;
}

/** PARI `crabs` (`buch1.c:972-976`) */
function crabs(a: CReal): CReal {
  return real_sign(a.re) < 0 ? cneg(a) : a;
}

/** PARI `gcdreal` (`buch1.c:978-994`) */
function gcdreal(a: CReal, b: CReal): CReal | null {
  if (!real_sign(a.re)) return crabs(b);
  if (!real_sign(b.re)) return crabs(a);
  if (real_expo(a.re) < -5) return crabs(b);
  if (real_expo(b.re) < -5) return crabs(a);
  a = crabs(a);
  b = crabs(b);
  while (real_expo(b.re) >= -5 && real_sign(b.re)) {
    const { z: q, e } = gcvtoi(divrr(a.re, b.re));
    if (e > 0) return null;
    const r = csub(a, cmulint(b, q));
    a = b;
    b = r;
  }
  return crabs(a);
}

const fupb_NONE = 0;
const fupb_RELAT = 1;
const fupb_PRECI = 2;

/** PARI `get_R` (`buch1.c:996-1020`) */
function get_R(B: BuchQuad, C: CVec, sreg: number, z: Real): { code: number; R: CReal | null } {
  let R: CReal = creal(real_1(B.PRECREG || DEFAULTPREC), 0n);
  if (B.PRECREG) {
    R = crabs(C[1]!);
    for (let i = 2; i <= sreg; i++) {
      const t = gcdreal(C[i]!, R);
      if (!t) return { code: fupb_PRECI, R: null };
      R = t;
    }
    if (real_sign(R.re) === 0 || real_expo(R.re) <= -3) return { code: fupb_RELAT, R: null };
  }
  const c = rtodbl(mulrr(z, R.re));
  if (c < 0.8 || c > 1.3) return { code: fupb_RELAT, R: null };
  return { code: fupb_NONE, R };
}

/** PARI `quad_be_honest` (`buch1.c:1022-1032`) */
function quad_be_honest(B: BuchQuad): boolean {
  if (B.KC2 <= B.KC) return true;
  return B.PRECREG ? real_be_honest(B) : imag_be_honest(B);
}

/** PARI `check_quaddisc` (`Qfb.c:22-31`) */
function check_quaddisc(x: bigint, f: string): number {
  const s = signe(x);
  if (Z_issquare(x)) throw new PariDomainError(f, 'issquare(disc)', '=', '1');
  let r = mod4(x);
  if (s < 0 && r) r = 4 - r;
  if (r > 1) throw new PariDomainError(f, 'disc % 4', '>', '1');
  return s;
}

/** the result of {@link Buchquad} / {@link quadclassunit0} */
export interface QuadClassUnit {
  /** class number */
  no: bigint;
  /** cyclic structure of the class group (elementary divisors, decreasing) */
  cyc: bigint[];
  /** generators of the cyclic factors, as binary quadratic forms */
  gen: Qfb[];
  /** regulator (1 for imaginary fields) */
  reg: Real;
  /**
   * `-1` if the fundamental unit has norm `-1` (real fields only), `1`
   * otherwise; absent for imaginary fields, exactly as PARI returns a
   * 4-component vector there.
   */
  sign?: bigint;
}

/**
 * PARI `Buchquad(D, c, c2, prec)` (`buch1.c:1035-1248`): class group, class
 * number and regulator of the quadratic order of discriminant `D`.
 */
export function Buchquad(D: bigint, cbach = 0, cbach2 = 0, prec = DEFAULTPREC): QuadClassUnit {
  const MAXRELSUP = 20;
  const SFB_MAX = 3;
  const RELSUP = 5;
  const s = check_quaddisc(D, 'Buchquad');
  const BQ: BuchQuad = {
    limhash: 0,
    KC: 0,
    KC2: 0,
    PRECREG: 0,
    primfact: new Array<number>(100).fill(0),
    exprimfact: new Array<number>(100).fill(0),
    hashtab: new Map<number, LPEntry>(),
    FB: [0],
    numFB: [],
    prodFB: 1n,
    powsubFB: [],
    vperm: [0],
    subFB: [0],
    badprim: null,
    q: { D, sqrtD: real_0(), isqrtD: 0n },
  };
  if (s < 0) {
    if (iabs(D) <= 4n) return { no: 1n, cyc: [], gen: [], reg: real_1(DEFAULTPREC) };
    BQ.PRECREG = 0;
  } else {
    /* PARI: maxss(prec+EXTRAPREC64, nbits2prec(2*expi(D)+128)) */
    BQ.PRECREG = Math.max(prec + 64, 2 * expi(D) + 128);
  }
  const drc = Number(iabs(D));
  const LOGD = Math.log(drc);
  const LOGD2 = LOGD * LOGD;
  const sdrc = Math.sqrt(drc);
  let lim = sdrc;
  if (!BQ.PRECREG) lim /= Math.sqrt(3);
  let cp = Math.floor(Math.exp(Math.sqrt((LOGD * Math.log(LOGD)) / 8)));
  if (cp < 20) cp = 20;
  if (cbach > 6) {
    if (cbach2 < cbach) cbach2 = cbach;
    cbach = 6;
  }
  if (cbach < 0) throw new PariDomainError('Buchquad', 'Bach constant', '<', '0');

  const minSFB = expi(D) > 15 ? 3 : 2;
  const GRHcheckS = init_GRHcheck(2, BQ.PRECREG ? 2 : 0, LOGD);
  let LIMC0 = Math.max(Math.floor(cbach2 * LOGD2), 1);
  let low = LIMC0;
  let high = LIMC0;
  const LIMCMAX = Math.floor(4 * LOGD2);
  cache_prime_quad(GRHcheckS, expi(D) < 16 ? 97 : 1223, D);
  while (!quadGRHchk(D, GRHcheckS, high)) {
    low = high;
    high *= 2;
  }
  while (high - low > 1) {
    const test = Math.floor((low + high) / 2);
    if (quadGRHchk(D, GRHcheckS, test)) high = test;
    else low = test;
  }
  let LIMC2: number;
  if (high === LIMC0 + 1 && quadGRHchk(D, GRHcheckS, LIMC0)) LIMC2 = LIMC0;
  else LIMC2 = high;
  if (LIMC2 > LIMCMAX) LIMC2 = LIMCMAX;
  LIMC0 = Math.floor(cbach * LOGD2);
  let LIMC = cbach ? LIMC0 : LIMC2;
  LIMC = Math.max(LIMC, nthidealquad(D, 2));

  let FIRST = 1;
  let W: ZMat | null = null;
  let h = 0n;
  let R: CReal | null = null;
  let cyc: ZC = [0n];
  let gen: Qfb[] = [];
  let nsubFB = 0;
  const st: HnfState = { dep: [[]], B: [[]], C: [] as CVec };

  START: for (;;) {
    do {
      if (!FIRST) LIMC = bnf_increase_LIMC(LIMC, LIMCMAX);
      FIRST = 0;
      BQ.hashtab.clear();
      if (LIMC < cp) LIMC = cp;
      if (LIMC2 < LIMC) LIMC2 = LIMC;
      if (BQ.PRECREG) BQ.q = qfr_data_init(D, BQ.PRECREG);
      FBquad(BQ, LIMC2, LIMC, GRHcheckS);
      BQ.subFB = subFBquad(BQ, D, lim + 0.5, minSFB);
      nsubFB = lg(BQ.subFB) - 1;
    } while (nsubFB < (expi(D) > 15 ? 3 : 2));
    /* invhr = 2^r1 (2pi)^r2 / (sqrt(D) w) ~ L(chi,1) / hR */
    const invhr = mulrr(
      dbltor((BQ.PRECREG ? 2 : Math.PI) / sdrc),
      compute_invresquad(GRHcheckS, LIMC)
    );
    BQ.powsubFB = powsubFBquad(BQ, CBUCH + 1);
    /* PARI: (LIMC & HIGHMASK)? (HIGHBIT>>1): LIMC*LIMC (buch1.c:1145) */
    BQ.limhash = LIMC < 3037000499 ? LIMC * LIMC : 2 ** 62;

    let need = BQ.KC + RELSUP - 2;
    const pc = { v: 0 };
    W = null;
    let sfb_trials = 0;
    let nreldep = 0;
    let nrelsup = 0;
    const sReal = nsubFB + RELSUP;

    do {
      if ((nreldep & 3) === 1 || (nrelsup & 7) === 1) {
        BQ.subFB = vecslice(BQ.vperm, 1, nsubFB);
        BQ.powsubFB = powsubFBquad(BQ, CBUCH + 1);
        BQ.hashtab.clear();
      }
      need += 2;
      const mat: zv[] = new Array<zv>(need + 1);
      mat[0] = [];
      const extraC: CVec = new Array<CReal>(need + 1);
      extraC[0] = creal(real_0(), 0n);
      let triv = 0;
      if (!W) {
        st.C = extraC;
        triv = trivial_relations(BQ, mat, extraC, BQ.PRECREG);
      }
      for (let i = triv + 1; i <= need; i++) {
        mat[i] = zero_zv(BQ.KC);
        extraC[i] = creal(real_0(BQ.PRECREG || DEFAULTPREC), 0n);
      }
      if (BQ.PRECREG) real_relations(BQ, need - triv, pc, sReal, LIMC, mat, extraC, triv);
      else imag_relations(BQ, need - triv, pc, LIMC, mat, triv);

      if (!W) W = hnfspec_i(mat, BQ.vperm, st, nsubFB);
      else W = hnfadd_i(W, BQ.vperm, st, mat, extraC);
      need = BQ.KC - (lg(W) - 1) - (lg(st.B) - 1);
      if (need) {
        if (++nreldep > 15 && cbach < 1) continue START;
        continue;
      }
      h = ZM_det_triangular(W);
      const sreg = lg(st.C) - 1 - (lg(st.B) - 1) - (lg(W) - 1);
      const gr = get_R(BQ, st.C, sreg, mulir(h, invhr));
      if (gr.code === fupb_PRECI) {
        BQ.PRECREG *= 2;
        FIRST = 1;
        continue START;
      }
      if (gr.code === fupb_RELAT) {
        if (++nrelsup > MAXRELSUP) {
          if (++sfb_trials > SFB_MAX && cbach <= 1) continue START;
          if (nsubFB < Math.min(10, BQ.KC)) nsubFB++;
        }
        need = Math.min(BQ.KC, nrelsup);
      } else R = gr.R;
    } while (need);

    if (!quad_be_honest(BQ)) continue;
    BQ.hashtab.clear();
    const cg = get_clgp(BQ, W!);
    cyc = cg.cyc;
    gen = cg.gen;
    break;
  }
  const cycArr: bigint[] = [];
  for (let i = 1; i < lg(cyc); i++) cycArr.push(cyc[i]!);
  const genArr: Qfb[] = [];
  for (let i = 1; i < lg(cyc); i++) genArr.push(gen[i]!);
  if (BQ.PRECREG)
    return {
      no: h,
      cyc: cycArr,
      gen: genArr,
      reg: R!.re,
      sign: mod2(R!.im) ? -1n : 1n,
    };
  return { no: h, cyc: cycArr, gen: genArr, reg: real_1(DEFAULTPREC) };
}

/**
 * PARI `quadclassunit0(x, flag, data, prec)` (`buch1.c:1260-1281`).
 * `data = [c1, c2]` are the Bach constants.
 */
export function quadclassunit0(
  x: bigint,
  flag = 0,
  data: number[] | null = null,
  prec = DEFAULTPREC
): QuadClassUnit {
  let c1 = 0;
  let c2 = 0;
  if (data) {
    if (data.length > 6) throw new PariDomainError('quadclassunit', 'tech vector', '>', '6');
    if (data.length > 0) c1 = data[0]!;
    if (data.length > 1) c2 = data[1]!;
  }
  if (flag) throw new NotImplementedError('SAGE_NOT_IMPLEMENTED: narrow class group');
  return Buchquad(x, c1, c2, prec);
}

/** PARI `quadclassno(D)` (`buch1.c:1282-1288`) */
export function quadclassno(D: bigint): bigint {
  return Buchquad(D, 0, 0, DEFAULTPREC).no;
}

/* ================================================================== */
/* bnfinit (buch2.c)                                                   */
/* ================================================================== */

/**
 * The result of {@link bnfinit}: the class group and the unit group of a
 * number field, in the shape of the components of PARI's `bnf` that
 * `bnfinit` computes (`buch2.c:4374-4420`).
 */
export interface Bnf {
  /** degree of the defining polynomial */
  degree: number;
  /** field discriminant `disc(K)` */
  disc: bigint;
  /** signature */
  r1: number;
  r2: number;
  /** class group: order, elementary divisors, generators (as forms, degree 2) */
  clgp: { no: bigint; cyc: bigint[]; gen: Qfb[] };
  /** regulator */
  reg: Real;
  /** number of roots of unity in `K` */
  tu: bigint;
  /** `-1` if the fundamental unit has norm `-1` (real quadratic only) */
  unitNorm?: bigint;
}

/** discriminant of the monic polynomial `x^2 + T[1] x + T[0]` */
function quadraticDisc(T: readonly bigint[]): bigint {
  return T[1]! * T[1]! - 4n * T[0]!;
}

/**
 * Fundamental discriminant `D0` with `disc = f^2 D0`: the discriminant of the
 * maximal order of `Q(sqrt(disc))`.
 */
function fundamentalDisc(disc: bigint, factor: (n: bigint) => Array<[bigint, bigint]>): bigint {
  let D = disc;
  for (;;) {
    let changed = false;
    for (const [p, e] of factor(iabs(D))) {
      if (e < 2n) continue;
      const p2 = p * p;
      const cand = D / p2;
      const r = mod4(cand);
      const s = signe(cand) < 0 && r ? 4 - r : r;
      if (s <= 1) {
        D = cand;
        changed = true;
        break;
      }
    }
    if (!changed) return D;
  }
}

/**
 * PARI `bnfinit` = `Buchall(P, flag, prec)` (`buch2.c:3700-3701`, driver
 * `Buchall_param`, `buch2.c:3946-4440`): class group, class number, regulator
 * and unit group of the number field `Q[x]/(T)`.
 *
 * `T` is a **monic** polynomial over `Z`, given by its coefficients in
 * increasing degree order (`[c0, c1, ..., 1]`).
 *
 * Implemented cases:
 * - `deg T <= 1`: `Buchall_deg1` (`buch2.c:3900-3930`), the trivial `bnf` of
 *   `Q`;
 * - `deg T == 2`: the whole computation is that of the quadratic order of the
 *   field discriminant, i.e. {@link Buchquad} (`buch1.c`), which is a complete
 *   and exact index calculus for that case.  PARI runs its general algorithm
 *   here as well and gets the same class group and the same regulator (with
 *   `Nrelid = 0`, cf. `buch2.c:3972`); only the auxiliary data attached to the
 *   `bnf` (units as algebraic numbers, `Vbase`, `C`, `W`, `B`, ...) differ.
 *
 * Degree `>= 3` throws: see the message for the exact list of missing upstream
 * routines.
 */
export function bnfinit(T: readonly bigint[], prec = DEFAULTPREC): Bnf {
  const n = T.length - 1;
  if (n >= 1 && T[n] !== 1n) throw new PariTypeError('bnfinit', 'nonmonic polynomial');
  if (n <= 1) {
    /* PARI Buchall_deg1 (buch2.c:3900): Q has trivial class group, R = 1 */
    return {
      degree: 1,
      disc: 1n,
      r1: 1,
      r2: 0,
      clgp: { no: 1n, cyc: [], gen: [] },
      reg: real_1(DEFAULTPREC),
      tu: 2n,
    };
  }
  if (n === 2) {
    /* eslint-disable-next-line @typescript-eslint/no-var-requires */
    const disc = quadraticDisc(T);
    if (disc === 0n || Z_issquare(disc))
      throw new PariDomainError('bnfinit', 'issquare(disc)', '=', '1');
    const D = fundamentalDisc(disc, (m) => Z_factor(m).map(([p, e]) => [p, BigInt(e)]));
    const r = Buchquad(D, 0, 0, prec);
    const tu = D === -3n ? 6n : D === -4n ? 4n : 2n;
    return {
      degree: 2,
      disc: D,
      r1: D > 0n ? 2 : 0,
      r2: D > 0n ? 0 : 1,
      clgp: { no: r.no, cyc: r.cyc, gen: r.gen },
      reg: r.reg,
      tu,
      ...(D > 0n ? { unitNorm: r.sign } : {}),
    };
  }
  throw new NotImplementedError(
    'SAGE_NOT_IMPLEMENTED: bnfinit of a field of degree > 2. PARI Buchall_param ' +
      '(buch2.c:3946) needs the number field layer, none of which is ported to ' +
      'parigp-ts yet: nfinit_basic/nfinit_complete (base1.c:2104,2143) and ' +
      'nfmaxord/round 4 (base2.c:462), the archimedean embeddings and the T2 form ' +
      '(nfmaxord_to_nf, base1.c:1829, via QX_complex_roots), idealprimedec ' +
      '(base2.c:2386) for FBgen (buch2.c:650) and nthideal (buch2.c:600), ideal ' +
      'arithmetic in HNF (idealhnf/idealmul/idealred, base4.c) for factorgen ' +
      '(buch2.c:852) and rnd_rel (buch2.c:2860), the Fincke-Pohst pseudo-minima of ' +
      'small_norm (buch2.c:2540) under the LLL-reduced T2 form (lll.c), ' +
      'cleanarch/getfu/makeunits (buch2.c:898,1126,1238) for the unit group and ' +
      'nfrootsof1 (base3.c) for the torsion. The relation-matrix half of the ' +
      'algorithm IS ported here (hnfspec_i/hnfadd_i/ZM_snf_group, i.e. ' +
      'class_group_gen, buch2.c:3730) and is exercised by Buchquad.'
  );
}
