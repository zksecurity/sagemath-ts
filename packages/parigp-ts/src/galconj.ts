/**
 * @module parigp-ts/galconj
 * @description Port of PARI/GP's Galois group machinery,
 * `reference/pari/src/basemath/galconj.c`.
 *
 * This provides `galoisinit` (the Galois group of a Galois number field as a
 * permutation group, via Allombert's "subgroup by subgroup" p-adic algorithm),
 * together with `galoispermtopol`, `galoisvecpermtopol`, `galoisfixedfield`
 * and `galoissubgroups`.
 *
 * Everything is exact: `bigint` integer and rational arithmetic only.  Upstream
 * uses floating point in exactly one place (`galoisborne`, to size the p-adic
 * accuracy); we replace that by an exact Hadamard/Cramer bound, see
 * {@link galoisborne}.
 *
 * Structure and function names mirror upstream; each function carries its
 * `galconj.c:<line>` (or `perm.c`, `Zp.c`, `FpX.c`, ...) citation.
 *
 * Supporting routines that upstream takes from other files (`perm.c` group and
 * permutation handling, `Zp.c` Hensel lifting, `FpX.c` polynomial arithmetic
 * over `F_p`) are ported here as well, in the sections below, because
 * `parigp-ts` has no `perm`/`Zp` module yet.
 *
 * NOT ported (each throws a `NotImplementedError` naming itself when it would
 * be reached):
 *  - `s4galoisgen` (galconj.c:1519) and `f36galoisgen` (galconj.c:1698), the
 *    dedicated searches for `S4` (degree 24) and `3x3:4` (degree 36).  They
 *    need `FpX_ffisom`/`FpXQ_ffisom_inv`/`FpXV_chinese`, which parigp-ts does
 *    not have.  `a4galoisgen` (degree 12) *is* ported.
 *  - `galoisgenlift_nilp` and the `pc_*` polycyclic layer (galconj.c:2389-2744),
 *    the "central" shortcut.  It is unreachable below degree 105, because
 *    `galoisanalysis` sets `ga_easy` for `n <= 104` (galconj.c:1104, 2827).
 *  - the `poliscyclo` shortcut of `galoisconj4_main` (galconj.c:2998), i.e.
 *    `galoiscyclo`: cyclotomic fields simply go through the generic algorithm,
 *    which returns the same group (with a different `l` and root labelling).
 *  - `galoisinitfromaut` (galconj.c:2926) and `galoissplittinginit`
 *    (galconj.c:2976); the latter needs `nfsplitting0` (base1.c:1413), which
 *    lives outside galconj.c.
 */

import {
  type FpX,
  FpX_add,
  FpX_degree,
  FpX_divrem,
  FpX_Fp_mul,
  FpX_gcd,
  FpX_mul,
  FpX_normalize,
  FpX_red,
  FpX_rem,
  FpX_renormalize,
  FpX_sub,
  pol_xn,
} from './ffinit.js';
import { isqrt, isPrime, NotImplementedError, Z_factor } from './ifactor.js';
import { PariBugError, PariImplError } from './polmodular.js';
import {
  PariDomainError,
  PariFlagError,
  PariInvError,
  PariTypeError,
} from './matkermod.js';

export {
  NotImplementedError,
  PariBugError,
  PariDomainError,
  PariFlagError,
  PariImplError,
  PariInvError,
  PariTypeError,
};

/** PARI `pari_err_IRREDPOL` */
export class PariIrredpolError extends Error {
  constructor(fun: string) {
    super(`not an irreducible polynomial in ${fun}`);
    this.name = 'PariIrredpolError';
  }
}

/* ================================================================== */
/*  Small integer helpers                                             */
/* ================================================================== */

/** `BITS_IN_LONG`; the `intheadlong` machinery of galconj.c is 64-bit. */
const BITS_IN_LONG = 64;

function bmod(a: bigint, m: bigint): bigint {
  const r = a % m;
  return r < 0n ? r + m : r;
}

function babs(a: bigint): bigint {
  return a < 0n ? -a : a;
}

function bmax(a: bigint, b: bigint): bigint {
  return a > b ? a : b;
}

/** gcd of two nonnegative machine integers (PARI `ugcd`) */
export function ugcd(a: number, b: number): number {
  a = Math.abs(a);
  b = Math.abs(b);
  while (b) {
    const t = a % b;
    a = b;
    b = t;
  }
  return a;
}

/** lcm of two nonnegative machine integers (PARI `ulcm`/`clcm`) */
export function ulcm(a: number, b: number): number {
  if (!a || !b) return 0;
  return (a / ugcd(a, b)) * b;
}

/** PARI `expu`: floor(log2(n)) for n >= 1 */
function expu(n: number): number {
  return 31 - Math.clz32(n);
}

/** PARI `usqrt`: floor(sqrt(n)) */
function usqrt(n: number): number {
  let r = Math.floor(Math.sqrt(n));
  while (r * r > n) r--;
  while ((r + 1) * (r + 1) <= n) r++;
  return r;
}

/**
 * modular inverse of `a` mod `m`, `gcd(a,m)=1`.
 * NOT exported: `ff.ts` already exports a function of this name.
 */
function Fp_inv(a: bigint, m: bigint): bigint {
  let [r0, r1] = [bmod(a, m), m];
  let [s0, s1] = [1n, 0n];
  while (r1 !== 0n) {
    const q = r0 / r1;
    [r0, r1] = [r1, r0 - q * r1];
    [s0, s1] = [s1, s0 - q * s1];
  }
  if (r0 !== 1n) throw new PariInvError('Fp_inv');
  return bmod(s0, m);
}

/** modular inverse of a small integer mod a small modulus (PARI `Fl_inv`) */
function Fl_inv(a: number, m: number): number {
  return Number(Fp_inv(BigInt(a), BigInt(m)));
}

function Fl_mul(a: number, b: number, m: number): number {
  return Number((BigInt(a) * BigInt(b)) % BigInt(m));
}

/** PARI `Fl_powu` */
function Fl_powu(a: number, e: number, m: number): number {
  let r = 1;
  let b = a % m;
  while (e > 0) {
    if (e & 1) r = Fl_mul(r, b, m);
    b = Fl_mul(b, b, m);
    e >>= 1;
  }
  return r;
}

/** PARI `logint(x, l)`: floor(log_l x) for x >= 1 */
export function logint(x: bigint, l: bigint): number {
  if (x < 1n) throw new PariDomainError('logint', 'x', '<', '1');
  let n = 0;
  let q = l;
  while (q <= x) {
    q *= l;
    n++;
  }
  return n;
}

/** factorisation of a small positive integer, as [prime, exponent] pairs */
export function factoru_small(n: number): Array<[number, number]> {
  const res: Array<[number, number]> = [];
  let m = n;
  for (let p = 2; p * p <= m; p++) {
    if (m % p) continue;
    let e = 0;
    while (m % p === 0) {
      m /= p;
      e++;
    }
    res.push([p, e]);
  }
  if (m > 1) res.push([m, 1]);
  return res;
}

/**
 * PARI `factoru_pow`: returns `[Fp, Fe, Fpe]`, 1-indexed vectors of the
 * primes, exponents and prime powers dividing `n`.
 */
function factoru_pow(n: number): { Fp: number[]; Fe: number[]; Fpe: number[] } {
  const f = factoru_small(n);
  const Fp = [0];
  const Fe = [0];
  const Fpe = [0];
  for (const [p, e] of f) {
    Fp.push(p);
    Fe.push(e);
    Fpe.push(Math.pow(p, e));
  }
  return { Fp, Fe, Fpe };
}

/** PARI `radicalu`: product of the primes dividing n */
function radicalu(n: number): number {
  let r = 1;
  for (const [p] of factoru_small(n)) r *= p;
  return r;
}

/** PARI `uisprimepower` */
function uisprimepower(n: number): boolean {
  if (n < 2) return false;
  return factoru_small(n).length === 1;
}

/** Euler phi of a small integer */
function eulerphiu(n: number): number {
  let r = n;
  for (const [p] of factoru_small(n)) r = (r / p) * (p - 1);
  return r;
}

/** iterator over primes >= start (PARI `u_forprime_init`/`u_forprime_next`) */
export class Forprime {
  private cur: number;
  constructor(start: number) {
    this.cur = Math.max(2, start) - 1;
  }
  next(): number {
    let n = this.cur + 1;
    if (n <= 2) {
      this.cur = 2;
      return 2;
    }
    if (n % 2 === 0) n++;
    for (;;) {
      if (isPrime(BigInt(n))) {
        this.cur = n;
        return n;
      }
      n += 2;
    }
  }
}

/* ================================================================== */
/*  ZX: dense integer polynomials, `f[i]` = coefficient of x^i         */
/* ================================================================== */

/** A polynomial with integer coefficients, `f[i]` = coefficient of `x^i`. */
export type ZX = bigint[];

export function ZX_renormalize(f: ZX): ZX {
  let d = f.length;
  while (d > 0 && f[d - 1] === 0n) d--;
  return d === f.length ? f.slice() : f.slice(0, d);
}

export function ZX_degree(f: ZX): number {
  return ZX_renormalize(f).length - 1;
}

export function ZX_add(a: ZX, b: ZX): ZX {
  const n = Math.max(a.length, b.length);
  const r: ZX = new Array(n).fill(0n);
  for (let i = 0; i < n; i++) r[i] = (a[i] ?? 0n) + (b[i] ?? 0n);
  return ZX_renormalize(r);
}

export function ZX_sub(a: ZX, b: ZX): ZX {
  const n = Math.max(a.length, b.length);
  const r: ZX = new Array(n).fill(0n);
  for (let i = 0; i < n; i++) r[i] = (a[i] ?? 0n) - (b[i] ?? 0n);
  return ZX_renormalize(r);
}

export function ZX_neg(a: ZX): ZX {
  return a.map((c) => -c);
}

export function ZX_mul(a: ZX, b: ZX): ZX {
  if (a.length === 0 || b.length === 0) return [];
  const r: ZX = new Array(a.length + b.length - 1).fill(0n);
  for (let i = 0; i < a.length; i++) {
    const ai = a[i]!;
    if (ai === 0n) continue;
    for (let j = 0; j < b.length; j++) r[i + j] += ai * b[j]!;
  }
  return ZX_renormalize(r);
}

export function ZX_Z_mul(a: ZX, c: bigint): ZX {
  if (c === 0n) return [];
  return a.map((x) => x * c);
}

export function ZX_deriv(f: ZX): ZX {
  const r: ZX = [];
  for (let i = 1; i < f.length; i++) r.push(f[i]! * BigInt(i));
  return ZX_renormalize(r);
}

export function ZX_equal(a: ZX, b: ZX): boolean {
  const x = ZX_renormalize(a);
  const y = ZX_renormalize(b);
  if (x.length !== y.length) return false;
  for (let i = 0; i < x.length; i++) if (x[i] !== y[i]) return false;
  return true;
}

export function ZX_is_monic(f: ZX): boolean {
  const g = ZX_renormalize(f);
  return g.length > 0 && g[g.length - 1] === 1n;
}

/** `x^n` as a ZX */
export function ZX_xn(n: number): ZX {
  const r: ZX = new Array(n + 1).fill(0n);
  r[n] = 1n;
  return r;
}

/** Euclidean division by a *monic* divisor; exact over Z. */
export function ZX_divrem_monic(a: ZX, b: ZX): [ZX, ZX] {
  const B = ZX_renormalize(b);
  if (B.length === 0 || B[B.length - 1] !== 1n)
    throw new PariTypeError('ZX_divrem_monic', 'divisor is not monic');
  const A = ZX_renormalize(a);
  if (A.length < B.length) return [[], A];
  const r = A.slice();
  const q: ZX = new Array(A.length - B.length + 1).fill(0n);
  for (let i = A.length - B.length; i >= 0; i--) {
    const c = r[i + B.length - 1]!;
    q[i] = c;
    if (c === 0n) continue;
    for (let j = 0; j < B.length; j++) r[i + j] -= c * B[j]!;
  }
  return [ZX_renormalize(q), ZX_renormalize(r.slice(0, B.length - 1))];
}

/** sup norm of the coefficients */
function ZX_supnorm(f: ZX): bigint {
  let m = 0n;
  for (const c of f) {
    const a = babs(c);
    if (a > m) m = a;
  }
  return m;
}

/**
 * Resultant of two ZX, computed modulo several primes and reconstructed by
 * CRT.  The Goldstein--Graham/Hadamard bound `|res| <= |A|_2^deg B * |B|_2^deg A`
 * tells us when to stop.  (PARI: `ZX_resultant`, polarit3.c:2333.)
 */
export function ZX_resultant(a: ZX, b: ZX): bigint {
  const A = ZX_renormalize(a);
  const B = ZX_renormalize(b);
  if (A.length === 0 || B.length === 0) return 0n;
  const dA = A.length - 1;
  const dB = B.length - 1;
  if (dA === 0 && dB === 0) return 1n;
  // Hadamard-type bound on |res|
  let n2A = 0n;
  for (const c of A) n2A += c * c;
  let n2B = 0n;
  for (const c of B) n2B += c * c;
  const normA = isqrt(n2A) + 1n;
  const normB = isqrt(n2B) + 1n;
  const bound = 2n * normA ** BigInt(dB) * normB ** BigInt(dA) + 1n;
  let mod = 1n;
  let res = 0n;
  const it = new Forprime(1 << 20);
  while (mod <= 2n * bound) {
    const p = BigInt(it.next());
    if (bmod(A[dA]!, p) === 0n || bmod(B[dB]!, p) === 0n) continue;
    const rp = FpX_resultant(FpX_red(A, p), FpX_red(B, p), p);
    // CRT
    if (mod === 1n) {
      res = rp;
      mod = p;
    } else {
      const t = bmod((rp - res) * Fp_inv(mod, p), p);
      res += t * mod;
      mod *= p;
    }
    res = bmod(res, mod);
  }
  // centre
  return res * 2n > mod ? res - mod : res;
}

/** resultant over F_p by the Euclidean algorithm */
function FpX_resultant(a: FpX, b: FpX, p: bigint): bigint {
  let A = FpX_red(a, p);
  let B = FpX_red(b, p);
  let res = 1n;
  for (;;) {
    const dA = FpX_degree(A);
    const dB = FpX_degree(B);
    if (dB < 0) return 0n;
    if (dB === 0) return (res * B[0]! ** BigInt(dA)) % p;
    const R = FpX_rem(A, B, p);
    const dR = FpX_degree(R);
    if (dR < 0) return 0n;
    // res(A,B) = (-1)^{dA dB} lc(B)^{dA-dR} res(B,R)
    let s = (B[dB]! ** BigInt(dA - dR)) % p;
    if (dA % 2 === 1 && dB % 2 === 1) s = bmod(-s, p);
    res = (res * s) % p;
    A = B;
    B = R;
  }
}

/** `disc(T)` for a monic ZX (PARI `ZX_disc`) */
export function ZX_disc(T: ZX): bigint {
  const t = ZX_renormalize(T);
  const n = t.length - 1;
  if (n <= 0) return 0n;
  const r = ZX_resultant(t, ZX_deriv(t));
  const s = ((n * (n - 1)) / 2) % 2 === 0 ? 1n : -1n;
  const lc = t[n]!;
  return (s * r) / lc;
}

export function ZX_is_squarefree(T: ZX): boolean {
  return ZX_disc(T) !== 0n;
}

/**
 * PARI `indexpartial` (base2.c:1895): a multiple of the denominator of an
 * algebraic integer of `Q[X]/(T)` written in the power basis.
 *
 * Upstream refines each prime power `p^(e/2) || DT` by
 * `ZpX_reduced_resultant_fast`; we keep `p^(e/2)` itself, which is still a
 * multiple of the denominator (only the p-adic accuracy suffers).
 */
export function indexpartial(T: ZX, DT?: bigint): bigint {
  const D = DT === undefined ? ZX_disc(T) : DT;
  if (D === 0n) throw new PariDomainError('indexpartial', 'disc', '=', '0');
  let res = 1n;
  for (const [p, e] of Z_factor(babs(D))) {
    const e2 = e >> 1n;
    res *= e2 >= 2n ? p ** e2 : p;
  }
  return res;
}

/* ================================================================== */
/*  FpX additions (FpX.c)                                             */
/* ================================================================== */

export function FpX_deriv(f: FpX, p: bigint): FpX {
  const r: FpX = [];
  for (let i = 1; i < f.length; i++) r.push((f[i]! * BigInt(i)) % p);
  return FpX_red(r, p);
}

export function FpX_eval(f: FpX, x: bigint, p: bigint): bigint {
  let r = 0n;
  const xx = bmod(x, p);
  for (let i = f.length - 1; i >= 0; i--) r = (r * xx + f[i]!) % p;
  return bmod(r, p);
}

/**
 * PARI `Fp_center`: representative in (-p/2, p/2].
 * NOT exported: `ff.ts` already exports a 2-argument function of this name.
 */
function Fp_center(a: bigint, p: bigint, p2: bigint): bigint {
  const r = bmod(a, p);
  return r > p2 ? r - p : r;
}

/** PARI `FpX_center_i` */
export function FpX_center(f: FpX, p: bigint, p2: bigint): ZX {
  return ZX_renormalize(f.map((c) => Fp_center(c, p, p2)));
}

/** PARI `FpX_div_by_X_x`: quotient of T by (x - a) */
export function FpX_div_by_X_x(T: FpX, a: bigint, p: bigint): FpX {
  const n = T.length - 1;
  if (n < 0) return [];
  const q: FpX = new Array(n).fill(0n);
  let c = T[n]!;
  for (let i = n - 1; i >= 0; i--) {
    q[i] = c;
    c = bmod(T[i]! + c * a, p);
  }
  return FpX_renormalize(q);
}

/** extended gcd over F_p: returns [d, u, v] with u*a + v*b = d */
export function FpX_extgcd(a: FpX, b: FpX, p: bigint): [FpX, FpX, FpX] {
  let r0 = FpX_red(a, p);
  let r1 = FpX_red(b, p);
  let s0: FpX = [1n];
  let s1: FpX = [];
  let t0: FpX = [];
  let t1: FpX = [1n];
  while (r1.length !== 0) {
    const [q, r] = FpX_divrem(r0, r1, p);
    [r0, r1] = [r1, r];
    [s0, s1] = [s1, FpX_sub(s0, FpX_mul(q, s1, p), p)];
    [t0, t1] = [t1, FpX_sub(t0, FpX_mul(q, t1, p), p)];
  }
  if (r0.length === 0) return [r0, s0, t0];
  const inv = Fp_inv(r0[r0.length - 1]!, p);
  return [FpX_Fp_mul(r0, inv, p), FpX_Fp_mul(s0, inv, p), FpX_Fp_mul(t0, inv, p)];
}

/** PARI `FpXQ_powers`: [1, x, x^2, ..., x^n] mod (T,p) */
export function FpXQ_powers(x: FpX, n: number, T: FpX, p: bigint): FpX[] {
  const V: FpX[] = [FpX_rem([1n], T, p)];
  for (let i = 1; i <= n; i++) V.push(FpX_rem(FpX_mul(V[i - 1]!, x, p), T, p));
  return V;
}

/** PARI `FpXQ_pow` (exported from ffinit under the same semantics) */
export function FpXQ_powBig(a: FpX, e: bigint, T: FpX, p: bigint): FpX {
  let r: FpX = FpX_rem([1n], T, p);
  let b = FpX_rem(a, T, p);
  let k = e;
  while (k > 0n) {
    if (k & 1n) r = FpX_rem(FpX_mul(r, b, p), T, p);
    b = FpX_rem(FpX_mul(b, b, p), T, p);
    k >>= 1n;
  }
  return r;
}

/** PARI `FpX_FpXQ_eval`: P(S) mod (T,p) */
export function FpX_FpXQ_eval(P: FpX, S: FpX, T: FpX, p: bigint): FpX {
  let r: FpX = [];
  for (let i = P.length - 1; i >= 0; i--) {
    r = FpX_rem(FpX_mul(r, S, p), T, p);
    if (P[i] !== 0n) r = FpX_add(r, [P[i]!], p);
  }
  return r;
}

/** PARI `FpX_Frobenius` (FpX.c:2213): `x^p mod (T,p)` */
export function FpX_Frobenius(T: FpX, p: bigint): FpX {
  return FpXQ_powBig(pol_xn(1), p, T, p);
}

/** PARI `FpXQ_autpow` (FpX.c:2354): `n`-fold composition of the automorphism `x` */
export function FpXQ_autpow(x: FpX, n: number, T: FpX, p: bigint): FpX {
  if (n === 0) return FpX_rem(pol_xn(1), T, p);
  let r = FpX_rem(x, T, p);
  for (let i = 1; i < n; i++) r = FpX_FpXQ_eval(r, FpX_rem(x, T, p), T, p);
  return r;
}

/**
 * PARI `FpXQ_autpowers` (FpX.c:2322): `V[i]` = `aut^(i-1)` (composition),
 * for `i = 1 .. f+1`.  1-indexed (`V[0]` is unused).
 */
export function FpXQ_autpowers(aut: FpX, f: number, T: FpX, p: bigint): FpX[] {
  const V: FpX[] = [[]];
  V.push(FpX_rem(pol_xn(1), T, p));
  if (f === 0) return V;
  V.push(FpX_rem(aut, T, p));
  for (let i = 3; i <= f + 1; i++) V.push(FpX_FpXQ_eval(V[i - 1]!, V[2]!, T, p));
  return V;
}

/** PARI `FpV_roots_to_pol` (FpX.c:1404): prod_i (x - V[i]); V is 1-indexed */
export function FpV_roots_to_pol(V: bigint[], p: bigint): FpX {
  let r: FpX = [1n];
  for (let i = 1; i < V.length; i++) r = FpX_mul(r, [bmod(-V[i]!, p), 1n], p);
  return r;
}

/** kernel of a matrix over F_p (rows x cols, 0-indexed); returns a basis */
function FpM_ker(M: bigint[][], nrows: number, ncols: number, p: bigint): bigint[][] {
  const A = M.map((r) => r.slice());
  const pivotOf: number[] = new Array(ncols).fill(-1);
  let row = 0;
  for (let col = 0; col < ncols && row < nrows; col++) {
    let piv = -1;
    for (let i = row; i < nrows; i++)
      if (bmod(A[i]![col]!, p) !== 0n) {
        piv = i;
        break;
      }
    if (piv < 0) continue;
    [A[row], A[piv]] = [A[piv]!, A[row]!];
    const inv = Fp_inv(A[row]![col]!, p);
    for (let j = 0; j < ncols; j++) A[row]![j] = (A[row]![j]! * inv) % p;
    for (let i = 0; i < nrows; i++) {
      if (i === row) continue;
      const c = bmod(A[i]![col]!, p);
      if (c === 0n) continue;
      for (let j = 0; j < ncols; j++) A[i]![j] = bmod(A[i]![j]! - c * A[row]![j]!, p);
    }
    pivotOf[col] = row;
    row++;
  }
  const basis: bigint[][] = [];
  for (let col = 0; col < ncols; col++) {
    if (pivotOf[col] !== -1) continue;
    const v: bigint[] = new Array(ncols).fill(0n);
    v[col] = 1n;
    for (let c2 = 0; c2 < ncols; c2++) {
      const r2 = pivotOf[c2]!;
      if (r2 === -1) continue;
      v[c2] = bmod(-A[r2]![col]!, p);
    }
    basis.push(v);
  }
  return basis;
}

/**
 * PARI `FpXQ_minpoly` (FpX.c:3095): minimal polynomial of `x` in F_p[X]/(T).
 * We use the straightforward linear algebra version.
 */
export function FpXQ_minpoly(x: FpX, T: FpX, p: bigint): FpX {
  const n = FpX_degree(T);
  const pows: FpX[] = [];
  let cur: FpX = FpX_rem([1n], T, p);
  for (let i = 0; i <= n; i++) {
    pows.push(cur);
    cur = FpX_rem(FpX_mul(cur, x, p), T, p);
  }
  // find the least d with 1, x, ..., x^d linearly dependent
  for (let d = 1; d <= n; d++) {
    // matrix rows = coefficient index, columns = power index 0..d
    const rows: bigint[][] = [];
    for (let i = 0; i < n; i++) {
      const r: bigint[] = [];
      for (let j = 0; j <= d; j++) r.push(pows[j]![i] ?? 0n);
      rows.push(r);
    }
    const K = FpM_ker(rows, n, d + 1, p);
    if (K.length === 0) continue;
    const v = K[0]!;
    return FpX_normalize(FpX_red(v, p), p);
  }
  throw new PariBugError('FpXQ_minpoly');
}

/** PARI `FpX_is_squarefree` */
export function FpX_is_squarefree(f: FpX, p: bigint): boolean {
  const g = FpX_red(f, p);
  if (FpX_degree(g) <= 0) return FpX_degree(g) === 0;
  const d = FpX_deriv(g, p);
  if (d.length === 0) return false;
  return FpX_degree(FpX_gcd(g, d, p)) === 0;
}

/** PARI `FpX_split_part`: gcd(x^p - x, f), the product of the linear factors */
export function FpX_split_part(f: FpX, p: bigint): FpX {
  const T = FpX_normalize(FpX_red(f, p), p);
  if (FpX_degree(T) <= 0) return T;
  const xp = FpXQ_powBig(pol_xn(1), p, T, p);
  return FpX_gcd(FpX_sub(xp, pol_xn(1), p), T, p);
}

/** PARI `Flx_nbroots`: number of roots of f in F_p */
export function FpX_nbroots(f: FpX, p: bigint): number {
  return FpX_degree(FpX_split_part(f, p));
}

/** PARI `Flx_is_totally_split` */
export function FpX_is_totally_split(f: FpX, p: bigint): boolean {
  return FpX_nbroots(f, p) === FpX_degree(FpX_red(f, p));
}

/** deterministic-ish pseudo random generator, so that factoring is repeatable */
class Rand {
  private s: number;
  constructor(seed: number) {
    this.s = seed >>> 0 || 1;
  }
  next(): number {
    // xorshift32
    let x = this.s;
    x ^= x << 13;
    x >>>= 0;
    x ^= x >>> 17;
    x ^= x << 5;
    x >>>= 0;
    this.s = x;
    return x;
  }
  bigint(p: bigint): bigint {
    let r = 0n;
    let bound = 1n;
    while (bound < p) {
      r = r * 4294967296n + BigInt(this.next());
      bound *= 4294967296n;
    }
    return r % p;
  }
}

/**
 * Distinct-degree factorisation of a monic squarefree `f` over F_p: returns
 * `parts[d]` = product of the irreducible factors of degree `d`.
 */
export function FpX_ddf(f: FpX, p: bigint): Map<number, FpX> {
  let T = FpX_normalize(FpX_red(f, p), p);
  const out = new Map<number, FpX>();
  let xq = pol_xn(1);
  let d = 0;
  while (FpX_degree(T) > 0) {
    d++;
    if (2 * d > FpX_degree(T)) {
      out.set(FpX_degree(T), T);
      break;
    }
    xq = FpXQ_powBig(xq, p, T, p);
    const g = FpX_gcd(FpX_sub(xq, pol_xn(1), p), T, p);
    if (FpX_degree(g) > 0) {
      out.set(d, g);
      T = FpX_divrem(T, g, p)[0];
      xq = FpX_rem(xq, T, p);
    }
  }
  return out;
}

/** Cantor--Zassenhaus equal degree splitting of `f` (all factors of degree d) */
function FpX_edf(f: FpX, d: number, p: bigint, rnd: Rand): FpX[] {
  const n = FpX_degree(f);
  if (n === d) return [f];
  const q = (p ** BigInt(d) - 1n) / 2n;
  for (;;) {
    const a: FpX = [];
    for (let i = 0; i < n; i++) a.push(rnd.bigint(p));
    const A = FpX_red(a, p);
    if (FpX_degree(A) <= 0) continue;
    let g = FpX_gcd(A, f, p);
    if (FpX_degree(g) === 0) {
      if (p === 2n) {
        // trace map
        let t: FpX = A;
        let s: FpX = A;
        for (let i = 1; i < d; i++) {
          s = FpX_rem(FpX_mul(s, s, p), f, p);
          t = FpX_add(t, s, p);
        }
        g = FpX_gcd(t, f, p);
      } else {
        g = FpX_gcd(FpX_sub(FpXQ_powBig(A, q, f, p), [1n], p), f, p);
      }
    }
    const dg = FpX_degree(g);
    if (dg <= 0 || dg === n) continue;
    const h = FpX_divrem(f, g, p)[0];
    return [...FpX_edf(FpX_normalize(g, p), d, p, rnd), ...FpX_edf(FpX_normalize(h, p), d, p, rnd)];
  }
}

/**
 * PARI `Flx_factor` restricted to *squarefree* input: the irreducible factors
 * of `f` over F_p, sorted (by degree, then lexicographically on coefficients)
 * so that the result is deterministic.
 */
export function FpX_factor_squarefree(f: FpX, p: bigint): FpX[] {
  const T = FpX_normalize(FpX_red(f, p), p);
  if (FpX_degree(T) <= 0) return [];
  if (!FpX_is_squarefree(T, p))
    throw new PariDomainError('FpX_factor_squarefree', 'issquarefree(f)', '=', '0');
  const rnd = new Rand(1);
  const parts = FpX_ddf(T, p);
  const res: FpX[] = [];
  for (const [d, g] of parts) res.push(...FpX_edf(g, d, p, rnd));
  res.sort(cmp_FpX);
  return res;
}

/** PARI `cmp_Flx`: compare by degree then by coefficients */
export function cmp_FpX(a: FpX, b: FpX): number {
  if (a.length !== b.length) return a.length - b.length;
  for (let i = a.length - 1; i >= 0; i--) {
    if (a[i]! !== b[i]!) return a[i]! < b[i]! ? -1 : 1;
  }
  return 0;
}

/** PARI `FpX_roots`: the roots of `f` in F_p, sorted increasingly */
export function FpX_roots(f: FpX, p: bigint): bigint[] {
  const g = FpX_split_part(f, p);
  if (FpX_degree(g) <= 0) return [];
  const rnd = new Rand(1);
  const fac = FpX_edf(FpX_normalize(g, p), 1, p, rnd);
  const r = fac.map((h) => bmod(-h[0]!, p));
  r.sort((x, y) => (x < y ? -1 : x > y ? 1 : 0));
  return r;
}

/**
 * PARI `Flx_nbfact_by_degree` (Flx.c): `D[i]` = number of irreducible factors
 * of degree `i`; also returns the total number of factors.
 */
export function FpX_nbfact_by_degree(f: FpX, p: bigint): { D: number[]; nb: number } {
  const n = FpX_degree(f);
  const D: number[] = new Array(n + 1).fill(0);
  let nb = 0;
  for (const [d, g] of FpX_ddf(f, p)) {
    const k = FpX_degree(g) / d;
    D[d] = k;
    nb += k;
  }
  return { D, nb };
}

/* ================================================================== */
/*  Zp: Hensel lifting (Zp.c)                                         */
/* ================================================================== */

/** PARI `ZpX_liftroot` (Zp.c:820): lift a simple root `a` of `f` mod p to p^e */
export function ZpX_liftroot(f: ZX, a: bigint, p: bigint, e: number): bigint {
  let q = p;
  let x = bmod(a, p);
  const df = ZX_deriv(f);
  while (q < p ** BigInt(e)) {
    const q2 = q * q > p ** BigInt(e) ? p ** BigInt(e) : q * q;
    const fx = bmod(FpX_eval(FpX_red(f, q2), x, q2), q2);
    const d = FpX_eval(FpX_red(df, q2), x, q2);
    x = bmod(x - fx * Fp_inv(d, q2), q2);
    q = q2;
  }
  return bmod(x, p ** BigInt(e));
}

/**
 * PARI `ZpX_roots` (Zp.c:755) for a `T` that is totally split mod `p`
 * (the only case galconj.c needs: `p` is chosen totally split).
 * Returns the roots as a 1-indexed vector of integers mod `p^e`.
 */
export function ZpX_roots(F: ZX, p: bigint, e: number): bigint[] {
  const r = FpX_roots(FpX_red(F, p), p);
  const out: bigint[] = [0n];
  for (const a of r) out.push(ZpX_liftroot(F, a, p, e));
  return out;
}

/** PARI `ZpX_liftroots` (Zp.c:833) */
export function ZpX_liftroots(f: ZX, S: bigint[], p: bigint, e: number): bigint[] {
  const out: bigint[] = [0n];
  for (let i = 1; i < S.length; i++) out.push(ZpX_liftroot(f, S[i]!, p, e));
  return out;
}

/**
 * Two-factor Hensel lift: given monic `A0`, `B0` over F_p with
 * `T = A0*B0 mod p` and `gcd(A0,B0) = 1`, return `[A,B]` monic mod `p^e` with
 * `T = A*B mod p^e`.  (PARI `ZpX_liftfact`/`MultiLift`, Zp.c:640.)
 */
function ZpX_lift2(T: ZX, A0: FpX, B0: FpX, p: bigint, e: number): [FpX, FpX] {
  const [g, s, t] = FpX_extgcd(A0, B0, p);
  if (FpX_degree(g) !== 0) throw new PariBugError('ZpX_lift2: factors are not coprime');
  let A: FpX = A0.slice();
  let B: FpX = B0.slice();
  let q = p;
  for (let k = 1; k < e; k++) {
    const qq = q * p;
    const Tq = FpX_red(T, qq);
    const prod = FpX_mul(A, B, qq);
    const E = FpX_sub(Tq, prod, qq);
    // E is divisible by q
    const Eq = FpX_red(
      E.map((c) => c / q),
      p
    );
    // a*B + b*A = Eq (mod p), deg a < deg A
    const a = FpX_rem(FpX_mul(t, Eq, p), A, p);
    const rest = FpX_sub(Eq, FpX_mul(a, B, p), p);
    const [b, rem] = FpX_divrem(rest, A, p);
    if (rem.length !== 0) throw new PariBugError('ZpX_lift2: inexact division');
    A = FpX_add(A, FpX_Fp_mul(a, q, qq), qq);
    B = FpX_add(B, FpX_Fp_mul(b, q, qq), qq);
    q = qq;
  }
  return [A, B];
}

/**
 * PARI `ZpX_liftfact` (Zp.c:640): lift the factorisation `Q` of `pol` mod `p`
 * to `p^e`.  `Q` is a 1-indexed vector of monic pairwise-coprime factors.
 */
export function ZpX_liftfact(pol: ZX, Q: FpX[], p: bigint, e: number): FpX[] {
  const pe = p ** BigInt(e);
  const g = Q.length - 1;
  if (g <= 1) return [[], FpX_red(pol, pe)];
  let rest: ZX = FpX_red(pol, pe);
  const out: FpX[] = [[]];
  for (let i = 1; i < g; i++) {
    // split rest = Q[i] * (rest / Q[i]) mod p
    let B0: FpX = [1n];
    for (let j = i + 1; j <= g; j++) B0 = FpX_mul(B0, Q[j]!, p);
    const [A, B] = ZpX_lift2(rest, Q[i]!, B0, p, e);
    out.push(A);
    rest = B;
  }
  out.push(rest);
  return out;
}

/**
 * PARI `bezout_lift_fact` (Zp.c:689): the Bezout coefficients of the lifted
 * factorisation, `U[i] = 1 mod Qlift[i]`, `0 mod Qlift[j]` for `j != i`.
 * 1-indexed.
 */
export function bezout_lift_fact(pol: ZX, Q: FpX[], p: bigint, e: number): FpX[] {
  const k = Q.length - 1;
  const pe = p ** BigInt(e);
  if (k === 1) return [[], [1n]];
  const F = ZpX_liftfact(pol, Q, p, e);
  const Tq = FpX_red(pol, pe);
  const out: FpX[] = [[]];
  for (let i = 1; i <= k; i++) {
    const [G, rem] = FpX_divrem(Tq, F[i]!, pe);
    if (rem.length !== 0) throw new PariBugError('bezout_lift_fact: inexact division');
    // H = G^-1 mod (F[i], p^e), by Hensel from mod p
    const Fi = F[i]!;
    const Fip = FpX_red(Fi, p);
    const [g0, u0] = FpX_extgcd(FpX_red(G, p), Fip, p);
    if (FpX_degree(g0) !== 0) throw new PariBugError('bezout_lift_fact: not coprime');
    let H = u0;
    let q = p;
    while (q < pe) {
      const qq = q * q > pe ? pe : q * q;
      const Gq = FpX_red(G, qq);
      const Fq = FpX_red(Fi, qq);
      // H <- H*(2 - G*H)
      const t = FpX_sub([2n], FpX_rem(FpX_mul(Gq, H, qq), Fq, qq), qq);
      H = FpX_rem(FpX_mul(H, t, qq), Fq, qq);
      q = qq;
    }
    out.push(FpX_rem(FpX_mul(G, H, pe), Tq, pe));
  }
  return out;
}

/**
 * PARI `ZpX_ZpXQ_liftroot` (Zp.c:1333): Newton lift of a root `S` of `P` in
 * `Z_p[x]/(T)` from precision 1 to `n`.  (The `early` exit of
 * `ZpX_ZpXQ_liftroot_ea` is a pure speed-up: the value it returns is the same
 * lift, reduced mod `p^n`.)
 */
export function ZpX_ZpXQ_liftroot(P: ZX, S: FpX, T: ZX, p: bigint, n: number): FpX {
  const pn = p ** BigInt(n);
  if (n === 1) return FpX_red(S, p);
  const dP = ZX_deriv(P);
  let q = p;
  let x = FpX_red(S, p);
  while (q < pn) {
    const qq = q * q > pn ? pn : q * q;
    const Tq = FpX_red(T, qq);
    const Pq = FpX_red(P, qq);
    const dq = FpX_red(dP, qq);
    const fx = FpX_FpXQ_eval(Pq, x, Tq, qq);
    const dx = FpX_FpXQ_eval(dq, x, Tq, qq);
    const inv = FpXQ_inv(dx, Tq, qq, p);
    x = FpX_sub(x, FpX_rem(FpX_mul(fx, inv, qq), Tq, qq), qq);
    q = qq;
  }
  return FpX_red(x, pn);
}

/** inverse in `(Z/p^e)[x]/(T)`, by Hensel from the inverse mod p */
export function FpXQ_inv(a: FpX, T: FpX, q: bigint, p: bigint): FpX {
  const Tp = FpX_red(T, p);
  const [g, u] = FpX_extgcd(FpX_red(a, p), Tp, p);
  if (FpX_degree(g) !== 0) throw new PariInvError('FpXQ_inv');
  let H = FpX_rem(u, Tp, p);
  let m = p;
  while (m < q) {
    const mm = m * m > q ? q : m * m;
    const Tm = FpX_red(T, mm);
    const am = FpX_red(a, mm);
    const t = FpX_sub([2n], FpX_rem(FpX_mul(am, H, mm), Tm, mm), mm);
    H = FpX_rem(FpX_mul(H, t, mm), Tm, mm);
    m = mm;
  }
  return FpX_red(H, q);
}

/**
 * PARI `FpV_invVandermonde` (FpX.c:1865): `den` times the inverse of the
 * Vandermonde matrix of `L`.  `L` is 1-indexed; the result is the matrix
 * `M[k][i]`, `k, i = 1..n`, with `M[k][i]` = coefficient of `x^(k-1)` in
 * `den * T(x)/((x - L[i]) T'(L[i]))`.
 */
export function FpV_invVandermonde(L: bigint[], den: bigint, p: bigint): bigint[][] {
  const n = L.length - 1;
  const T = FpV_roots_to_pol(L, p);
  const dT = FpX_deriv(T, p);
  const M: bigint[][] = [];
  for (let k = 0; k <= n; k++) M.push(new Array(n + 1).fill(0n));
  for (let i = 1; i <= n; i++) {
    const R = (bmod(den, p) * Fp_inv(FpX_eval(dT, L[i]!, p), p)) % p;
    const P = FpX_Fp_mul(FpX_div_by_X_x(T, L[i]!, p), R, p);
    for (let k = 1; k <= n; k++) M[k]![i] = P[k - 1] ?? 0n;
  }
  return M;
}

/* ================================================================== */
/*  Permutations (perm.c)                                             */
/*                                                                    */
/*  As in PARI, all vectors here are 1-indexed: entry 0 is a dummy.    */
/* ================================================================== */

/** A permutation of `{1..n}`, 1-indexed (`p[0]` is a dummy 0). */
export type Perm = number[];

/** PARI `identity_perm` (pariinl.h:1137) */
export function identity_perm(n: number): Perm {
  const v: Perm = new Array(n + 1).fill(0);
  for (let i = 1; i <= n; i++) v[i] = i;
  return v;
}

/** PARI `perm_mul` (pariinl.h:1151): `(s*x)[i] = s[x[i]]` */
export function perm_mul(s: Perm, x: Perm): Perm {
  const v: Perm = new Array(x.length).fill(0);
  for (let i = 1; i < x.length; i++) v[i] = s[x[i]!]!;
  return v;
}

/** PARI `perm_sqr` */
export function perm_sqr(x: Perm): Perm {
  return perm_mul(x, x);
}

/** PARI `perm_inv` */
export function perm_inv(x: Perm): Perm {
  const v: Perm = new Array(x.length).fill(0);
  for (let i = 1; i < x.length; i++) v[x[i]!] = i;
  return v;
}

/** PARI `perm_conj`: `s*t*s^-1` */
export function perm_conj(s: Perm, t: Perm): Perm {
  const v: Perm = new Array(s.length).fill(0);
  for (let i = 1; i < s.length; i++) v[s[i]!] = s[t[i]!]!;
  return v;
}

/** PARI `perm_commute` (perm.c:766) */
export function perm_commute(s: Perm, t: Perm): boolean {
  for (let i = 1; i < t.length; i++) if (t[s[i]!] !== s[t[i]!]) return false;
  return true;
}

/** PARI `perm_powu` (perm.c:693) */
export function perm_powu(perm: Perm, exp: number): Perm {
  const r = perm.length - 1;
  const p: Perm = new Array(r + 1).fill(0);
  const v: Perm = new Array(r + 1).fill(0);
  for (let i = 1; i <= r; i++) {
    if (p[i]) continue;
    v[1] = i;
    let n = 1;
    for (let k = perm[i]!; k !== i; k = perm[k]!, n++) v[n + 1] = k;
    const e = exp % n;
    let l = e;
    for (let k = 1; k <= n; k++) {
      p[v[k]!] = v[l + 1]!;
      if (++l === n) l = 0;
    }
  }
  return p;
}

/** PARI `vecperm_orbits_i` (perm.c:417): orbits of `<v>` acting on `{1..n}` */
export function vecperm_orbits(v: Perm[], n: number): Perm[] {
  const cycle: Perm[] = [[]];
  const bit: number[] = new Array(n + 2).fill(0);
  let mj = 1;
  for (let k = 1; k <= n; ) {
    let m = 1;
    const cy: number[] = [0];
    for (; bit[mj]; mj++);
    k++;
    cy[m++] = mj;
    bit[mj++] = 1;
    for (;;) {
      const mold = m;
      for (let o = 0; o < v.length; o++) {
        const vo = v[o]!;
        for (let pp = 1; pp < m; pp++) {
          const j = vo[cy[pp]!]!;
          if (!bit[j]) cy[m++] = j;
          bit[j] = 1;
        }
      }
      if (m === mold) break;
      k += m - mold;
    }
    cycle.push(cy);
  }
  return cycle;
}

/** PARI `perm_cycles` (perm.c:478) */
export function perm_cycles(v: Perm): Perm[] {
  return vecperm_orbits([v], v.length - 1);
}

/** PARI `perm_orderu` (perm.c:493) */
export function perm_orderu(v: Perm): number {
  const c = perm_cycles(v);
  let d = 1;
  for (let i = 1; i < c.length; i++) d = ulcm(d, c[i]!.length - 1);
  return d;
}

/** PARI `cyc_pow` (perm.c:616) */
export function cyc_pow(cyc: Perm[], exp: number): Perm[] {
  const c: Perm[] = [[]];
  for (let j = 1; j < cyc.length; j++) {
    const v = cyc[j]!;
    const n = v.length - 1;
    let e = exp % n;
    if (e < 0) e += n;
    const g = ugcd(n, e);
    const m = n / g;
    for (let i = 0; i < g; i++) {
      const p: number[] = new Array(m + 1).fill(0);
      let l = i;
      for (let k = 1; k <= m; k++) {
        p[k] = v[l + 1]!;
        l += e;
        if (l >= n) l -= n;
      }
      c.push(p);
    }
  }
  return c;
}

/** PARI `vecpermute`: `A[x[i]]`, both 1-indexed */
export function vecpermute<T>(A: T[], x: number[]): T[] {
  const r: T[] = [A[0]!];
  for (let i = 1; i < x.length; i++) r.push(A[x[i]!]!);
  return r;
}

/** lexicographic comparison of 1-indexed integer vectors (PARI `vecsmall_lexcmp`) */
export function vecsmall_lexcmp(a: number[], b: number[]): number {
  const l = Math.min(a.length, b.length);
  for (let i = 1; i < l; i++) if (a[i] !== b[i]) return a[i]! < b[i]! ? -1 : 1;
  return a.length - b.length;
}

export function zv_equal(a: number[], b: number[]): boolean {
  return vecsmall_lexcmp(a, b) === 0;
}

/** PARI `vecsmall_uniq`: sorted, duplicates removed (1-indexed) */
export function vecsmall_uniq(v: number[]): number[] {
  const s = v.slice(1).sort((a, b) => a - b);
  const out: number[] = [0];
  for (const x of s) if (out.length === 1 || out[out.length - 1] !== x) out.push(x);
  return out;
}

/* ================================================================== */
/*  Groups (perm.c)                                                   */
/* ================================================================== */

/**
 * A group as PARI stores it: `[gen, ord]`, a polycyclic presentation.
 * `gen` and `ord` are 1-indexed.
 */
export interface Group {
  gen: Perm[];
  ord: number[];
}

/** PARI `group_order` (perm.c:818) */
export function group_order(G: Group): number {
  let r = 1;
  for (let i = 1; i < G.ord.length; i++) r *= G.ord[i]!;
  return r;
}

/** PARI `group_domain` (perm.c:824) */
export function group_domain(G: Group): number {
  if (G.gen.length < 2) throw new PariDomainError('group_domain', '#G', '=', '1');
  return G.gen[1]!.length - 1;
}

/** PARI `trivialgroup` (perm.c:945) */
export function trivialgroup(): Group {
  return { gen: [[]], ord: [0] };
}

/** PARI `cyclicgroup` (perm.c:948) */
export function cyclicgroup(g: Perm, s: number): Group {
  return { gen: [[], g.slice()], ord: [0, s] };
}

/** PARI `dicyclicgroup` (perm.c:953) */
export function dicyclicgroup(g1: Perm, g2: Perm, s1: number, s2: number): Group {
  return { gen: [[], g1.slice(), g2.slice()], ord: [0, s1, s2] };
}

/** PARI `group_elts` (perm.c:866): all elements, from the presentation */
export function group_elts(G: Group, n: number): Perm[] {
  const gen = G.gen;
  const ord = G.ord;
  const res: Perm[] = [[]];
  res.push(identity_perm(n));
  let k = 1;
  for (let i = 1; i < gen.length; i++) {
    const c = k * (ord[i]! - 1);
    res.push(gen[i]!.slice());
    k++;
    for (let j = 2; j <= c; j++) {
      res.push(perm_mul(res[j]!, gen[i]!));
      k++;
    }
  }
  return res;
}

/** PARI `group_leftcoset` (perm.c:836): `gG` */
export function group_leftcoset(G: Group, g: Perm): Perm[] {
  const gen = G.gen;
  const ord = G.ord;
  const res: Perm[] = [[], g.slice()];
  let k = 1;
  for (let i = 1; i < gen.length; i++) {
    const c = k * (ord[i]! - 1);
    for (let j = 1; j <= c; j++) {
      res.push(perm_mul(res[j]!, gen[i]!));
      k++;
    }
  }
  return res;
}

/** PARI `group_rightcoset` (perm.c:852): `Gg` */
export function group_rightcoset(G: Group, g: Perm): Perm[] {
  const gen = G.gen;
  const ord = G.ord;
  const res: Perm[] = [[], g.slice()];
  let k = 1;
  for (let i = 1; i < gen.length; i++) {
    const c = k * (ord[i]! - 1);
    for (let j = 1; j <= c; j++) {
      res.push(perm_mul(gen[i]!, res[j]!));
      k++;
    }
  }
  return res;
}

/** PARI `perm_generate` (perm.c:808) */
export function perm_generate(S: Perm, H: Perm[], o: number): Perm[] {
  const n = H.length - 1;
  const L: Perm[] = [[]];
  for (let i = 1; i <= n; i++) L.push(H[i]!.slice());
  for (let i = n + 1; i <= n * o; i++) L.push(perm_mul(L[i - n]!, S));
  return L;
}

/** PARI `groupelts_set` (perm.c:897): bitmap of `{g[1] : g in elts}` */
export function groupelts_set(elts: Perm[], n: number): boolean[] {
  const res: boolean[] = new Array(n + 1).fill(false);
  for (let i = 1; i < elts.length; i++) res[elts[i]![1]!] = true;
  return res;
}

/** PARI `group_set` (perm.c:906) */
export function group_set(G: Group, n: number): boolean[] {
  return groupelts_set(group_elts(G, n), n);
}

/** PARI `perm_relorder` (perm.c:798): order of `p` modulo the group `set` */
export function perm_relorder(p: Perm, set: boolean[]): number {
  let n = 1;
  let q = p[1]!;
  while (!set[q]) {
    q = p[q]!;
    n++;
  }
  return n;
}

/** PARI `group_perm_normalize` (perm.c:1081) */
export function group_perm_normalize(N: Group, g: Perm): boolean {
  const a = group_leftcoset(N, g).slice(1).sort(vecsmall_lexcmp);
  const b = group_rightcoset(N, g).slice(1).sort(vecsmall_lexcmp);
  if (a.length !== b.length) return false;
  for (let i = 0; i < a.length; i++) if (!zv_equal(a[i]!, b[i]!)) return false;
  return true;
}

/** A quotient map, PARI `[gen, coset]` (see perm.c:1010) */
export interface Quotient {
  gen: Perm[];
  coset: number[];
}

/** PARI `groupelts_quotient` (perm.c:977) */
export function groupelts_quotient(elt: Perm[], H: Group): Quotient {
  const n = elt[1]!.length - 1;
  const o = group_order(H);
  const le = elt.length - 1;
  const used: boolean[] = new Array(le + 2).fill(false);
  const l = le / o;
  const p2: Perm[] = [[]];
  const p3: number[] = new Array(n + 1).fill(0);
  const el: number[] = new Array(n + 1).fill(0);
  let a = 1;
  for (let i = 1; i <= le; i++) el[elt[i]![1]!] = i;
  for (let i = 1; i <= l; i++) {
    while (used[a]) a++;
    const V = group_leftcoset(H, elt[a]!);
    p2.push(V[1]!);
    for (let j = 1; j < V.length; j++) {
      const b = el[V[j]![1]!]!;
      if (b === 0) throw new PariImplError('group_quotient for a non-WSS group');
      used[b] = true;
    }
    for (let j = 1; j <= o; j++) p3[V[j]![1]!] = i;
  }
  return { gen: p2, coset: p3 };
}

/** PARI `group_quotient` (perm.c:1005) */
export function group_quotient(G: Group, H: Group): Quotient {
  return groupelts_quotient(group_elts(G, group_domain(G)), H);
}

/** PARI `quotient_perm` (perm.c:1001) */
export function quotient_perm(C: Quotient, p: Perm): Perm {
  const gen = C.gen;
  const coset = C.coset;
  const p3: Perm = new Array(gen.length).fill(0);
  for (let j = 1; j < gen.length; j++) {
    p3[j] = coset[p[gen[j]![1]!]!]!;
    if (p3[j] === 0) throw new PariImplError('quotient_perm for a non-WSS group');
  }
  return p3;
}

/** PARI `quotient_subgroup_lift` (perm.c:1021) */
export function quotient_subgroup_lift(C: Quotient, H: Group, S: Group): Group {
  const genH = H.gen;
  const genS = S.gen;
  const genC = C.gen;
  const L: Perm[] = [[]];
  for (let j = 1; j < genH.length; j++) L.push(genH[j]!);
  for (let j = 1; j < genS.length; j++) L.push(genC[genS[j]![1]!]!);
  const ord: number[] = [0];
  for (let j = 1; j < H.ord.length; j++) ord.push(H.ord[j]!);
  for (let j = 1; j < S.ord.length; j++) ord.push(S.ord[j]!);
  return { gen: L, ord };
}

/** PARI `quotient_group` (perm.c:1039) */
export function quotient_group(C: Quotient, G: Group): Group {
  const Cgen = C.gen;
  const Ggen = G.gen;
  const n = Cgen.length - 1;
  const Qord: number[] = [0];
  const Qgen: Perm[] = [[]];
  let Qelt: Perm[] = [[], identity_perm(n)];
  let Qset = groupelts_set(Qelt, n);
  for (let i = 1; i < Ggen.length; i++) {
    const g = quotient_perm(C, Ggen[i]!);
    const o = perm_relorder(g, Qset);
    if (o !== 1) {
      Qgen.push(g);
      Qord.push(o);
      Qelt = perm_generate(g, Qelt, o);
      Qset = groupelts_set(Qelt, n);
    }
  }
  return { gen: Qgen, ord: Qord };
}

/** PARI `trivialsubgroups` (perm.c:793) */
function trivialsubgroups(): Group[] {
  return [trivialgroup()];
}

/** PARI `liftlistsubgroups` (perm.c:1091) */
function liftlistsubgroups(L: Group[], C: Perm[], r: number): Group[] {
  const c = C.length - 1;
  if (L.length === 0) return [];
  const n = C[1]!.length - 1;
  const R: Group[] = [];
  for (const S of L) {
    const Selt = group_set(S, n);
    for (let j = 1; j <= c; j++) {
      const p = C[j]!;
      if (perm_relorder(p, Selt) === r && group_perm_normalize(S, p)) {
        const gen = S.gen.slice();
        gen.push(p);
        const ord = S.ord.slice();
        ord.push(r);
        R.push({ gen, ord });
      }
    }
  }
  return R;
}

/** PARI `liftsubgroup` (perm.c:1120) */
function liftsubgroup(C: Quotient, H: Group, S: Group): Group[] {
  let V: Group[] = trivialsubgroups();
  const Sgen = S.gen;
  const Sord = S.ord;
  const Cgen = C.gen;
  for (let i = 1; i < Sgen.length; i++) {
    const W = group_leftcoset(H, Cgen[Sgen[i]![1]!]!);
    V = liftlistsubgroups(V, W, Sord[i]!);
  }
  return V;
}

/** PARI `group_isA4S4` (perm.c:1137): 1:A4, 2:S4, 3:F36, 0: other */
export function group_isA4S4(G: Group): number {
  const elt = G.gen;
  const ord = G.ord;
  const n = ord.length;
  if (n !== 4 && n !== 5) return 0;
  if (n === 4 && ord[1] === 3 && ord[2] === 3 && ord[3] === 4) {
    const p = elt[1]!;
    const q = elt[2]!;
    const r = elt[3]!;
    for (let i = 1; i <= 36; i++) if (p[r[i]!] !== r[q[i]!]) return 0;
    return 3;
  }
  if (ord[1] !== 2 || ord[2] !== 2 || ord[3] !== 3) return 0;
  if (perm_commute(elt[1]!, elt[3]!)) return 0;
  if (n === 4) return 1;
  if (ord[4] !== 2) return 0;
  if (perm_commute(elt[3]!, elt[4]!)) return 0;
  return 2;
}

/** PARI `group_subgroups` (perm.c:1162): all subgroups of `G` */
export function group_subgroups(G: Group): Group[] {
  const gen = G.gen;
  const ord = G.ord;
  const n = gen.length;
  let H: Group;
  let sg3: Group[] | null;
  if (n === 1) return trivialsubgroups();
  const t = group_isA4S4(G);
  if (t === 3) {
    const HH: Group = {
      gen: [[], gen[1]!, gen[2]!, perm_sqr(gen[3]!)],
      ord: [0, 3, 3, 2],
    };
    const S = group_subgroups(HH);
    const V: Group[] = [];
    V.push(cyclicgroup(gen[3]!, 4));
    for (let i = 2; i < 10; i++)
      V.push(cyclicgroup(perm_mul(V[i - 2]!.gen[1]!, gen[i % 3 === 1 ? 2 : 1]!), 4));
    V.push(G);
    return [...S, ...V];
  } else if (t) {
    const s = gen[1]!;
    const tt = gen[2]!;
    const st = perm_mul(s, tt);
    H = dicyclicgroup(s, tt, 2, 2);
    sg3 = [cyclicgroup(s, 2), cyclicgroup(tt, 2), cyclicgroup(st, 2)];
    if (n === 5) {
      let u = gen[3]!;
      let v = gen[4]!;
      let w: Perm;
      let u2: Perm;
      if (zv_equal(perm_conj(u, s), tt)) u2 = perm_sqr(u);
      else {
        u2 = u;
        u = perm_sqr(u);
      }
      if (perm_orderu(v) === 2) {
        if (!perm_commute(s, v)) {
          v = perm_conj(u, v);
          if (!perm_commute(s, v)) v = perm_conj(u, v);
        }
        w = perm_mul(v, tt);
      } else {
        w = v;
        if (!zv_equal(perm_sqr(w), s)) {
          w = perm_conj(u, w);
          if (!zv_equal(perm_sqr(w), s)) w = perm_conj(u, w);
        }
        v = perm_mul(w, tt);
      }
      sg3.push(dicyclicgroup(s, v, 2, 2));
      sg3.push(dicyclicgroup(tt, perm_conj(u, v), 2, 2));
      sg3.push(dicyclicgroup(st, perm_conj(u2, v), 2, 2));
      sg3.push(dicyclicgroup(s, w, 2, 2));
      sg3.push(dicyclicgroup(tt, perm_conj(u, w), 2, 2));
      sg3.push(dicyclicgroup(st, perm_conj(u2, w), 2, 2));
    }
  } else {
    const osig = factoru_small(ord[1]!)[0]![0];
    const sig = perm_powu(gen[1]!, ord[1]! / osig);
    H = cyclicgroup(sig, osig);
    sg3 = null;
  }
  const C = group_quotient(G, H);
  const Q = quotient_group(C, G);
  const M = group_subgroups(Q);
  const sg1: Group[] = M.map((m) => quotient_subgroup_lift(C, H, m));
  const sg2: Group[] = [];
  for (const m of M) sg2.push(...liftsubgroup(C, H, m));
  let p1 = [...sg1, ...sg2];
  if (sg3) {
    p1 = [...p1, ...sg3];
    if (n === 5) {
      /* ensure that the D4 subgroups of S4 are in supersolvable format */
      for (let j = 2; j <= 4; j++) {
        const c = p1[j]!.gen;
        if (!perm_commute(c[1]!, c[3]!)) {
          if (perm_commute(c[2]!, c[3]!)) {
            const tmp = c[1]!;
            c[1] = c[2]!;
            c[2] = tmp;
          } else c[1] = perm_mul(c[2]!, c[1]!);
        }
      }
    }
  }
  return p1;
}

/* ================================================================== */
/*  galconj.c                                                         */
/* ================================================================== */

/** a polynomial with rational coefficients: `num[i]/den` is the coeff of x^i */
export interface QPoly {
  num: ZX;
  den: bigint;
}

/** reduce a QPoly to lowest terms (global denominator) */
export function QPoly_normalize(q: QPoly): QPoly {
  let g = q.den;
  for (const c of q.num) {
    let a = babs(c);
    while (a) {
      const t = g % a;
      g = a;
      a = t;
    }
    if (g === 1n) break;
  }
  if (g <= 1n) return { num: ZX_renormalize(q.num), den: q.den };
  return { num: ZX_renormalize(q.num.map((c) => c / g)), den: q.den / g };
}

/** the coefficients of a QPoly as reduced fractions `[num, den]` */
export function QPoly_to_fractions(q: QPoly): Array<[bigint, bigint]> {
  return q.num.map((c) => {
    let a = babs(c);
    let b = q.den;
    while (a) {
      const t = b % a;
      b = a;
      a = t;
    }
    return [c / b, q.den / b] as [bigint, bigint];
  });
}

/** PARI `RgX_to_FpX` for a QPoly */
export function QPoly_to_FpX(q: QPoly, p: bigint): FpX {
  const inv = Fp_inv(q.den, p);
  return FpX_red(
    q.num.map((c) => bmod(c * inv, p)),
    p
  );
}

export interface GaloisBorne {
  l: bigint;
  valsol: number;
  valabs: number;
  bornesol: bigint;
  ladicsol: bigint;
  ladicabs: bigint;
  dis: bigint;
}

interface GaloisLift {
  T: ZX;
  den: bigint;
  p: bigint;
  L: bigint[];
  Lden: bigint[];
  e: number;
  Q: bigint;
  TQ: FpX;
  gb: GaloisBorne;
}

interface GaloisTestlift {
  n: number;
  f: number;
  g: number;
  bezoutcoeff: FpX[];
  pauto: FpX[];
  C: Array<Array<FpX | null>>;
  Cd: bigint[][];
}

interface GaloisTest {
  order: number[];
  borne: bigint;
  lborne: bigint;
  ladic: bigint;
  PV: Array<bigint[][] | null>;
  M: bigint[][];
  L: bigint[];
}

/* result of the study of Frobenius degrees (galconj.c:117) */
const ga_all_normal = 1;
const ga_ext_2 = 2;
const ga_non_wss = 4;
const ga_all_nilpotent = 8;
const ga_easy = 16;

interface GaloisAnalysis {
  p: number;
  deg: number;
  ord: number;
  l: number;
  p4: number;
  group: number;
}

interface GaloisFrobenius {
  p: number;
  fp: number;
  deg: number;
  Tmod: ZX[];
  psi: number[];
}

/** a group in the "std" form PARI's galoisgen returns: generators + relative orders */
export interface GaloisGens {
  gen: Perm[];
  orders: number[];
}

/* ------------------------------------------------------------------ */
/* Coefficient bounds (galconj.c:200-282)                             */
/* ------------------------------------------------------------------ */

/**
 * Exact replacement for PARI's `initgaloisborne` (galconj.c:200), which
 * computes complex roots and `||den * V^-1||_oo` in floating point.
 *
 * We bound the same quantity exactly: if `c` is the coefficient vector of an
 * element of `Q[x]/(T)` whose conjugates are all bounded by `W`, then Cramer's
 * rule gives `c_k = det(V with column k replaced by w)/det(V)` with
 * `|det V| = sqrt|disc T|`, and Hadamard's inequality bounds the numerator by
 * `n^(n/2) B^(n(n-1)/2) W` where `B` bounds the roots (Cauchy's bound).  Hence
 *
 *     |den * c_k| <= borne * W,  borne = den n^(n/2) B^(n(n-1)/2) / sqrt|disc|.
 *
 * Returns `{den, borne, borneroots}` with `borneroots = B`.
 */
export function initgaloisborne(
  T: ZX,
  dn: bigint | null
): { den: bigint; borne: bigint; borneroots: bigint; D: bigint } {
  const n = ZX_degree(T);
  const D = ZX_disc(T);
  if (D === 0n) throw new PariDomainError('galoisinit', 'issquarefree(pol)', '=', '0');
  const den = dn === null ? indexpartial(T, D) : dn;
  /* Cauchy's bound on the roots of a monic T */
  const B = 1n + ZX_supnorm(T.slice(0, n));
  const nn = BigInt(n) ** BigInt(n);
  const num = den * (isqrt(nn) + 1n) * B ** BigInt((n * (n - 1)) / 2);
  const den2 = isqrt(babs(D));
  const borne = (num + den2 - 1n) / den2;
  return { den, borne, borneroots: B, D };
}

/** PARI `galoisborne` (galconj.c:244) */
export function galoisborne(T: ZX, dn: bigint | null, gb: GaloisBorne, d: number): bigint {
  const step = 3;
  const n = ZX_degree(T);
  const { den, borne, borneroots, D } = initgaloisborne(T, dn);
  gb.dis = D;
  const dnn = dn === null ? den : dn;
  const bornetrace =
    BigInt(Math.floor((2 * step * n) / d)) * borneroots ** BigInt(Math.min(n, step));
  const br = borne * borneroots;
  const borneabs = bmax(borne * bornetrace, bornetrace ** BigInt(d));
  /* We use d-1 tests, so we must overlift to 2^BITS_IN_LONG */
  gb.valsol = logint(br << BigInt(2 + BITS_IN_LONG), gb.l) + 1;
  gb.valabs = logint(borneabs << 2n, gb.l) + 1;
  gb.valabs = Math.max(gb.valsol, gb.valabs);
  gb.bornesol = br << 1n;
  gb.ladicsol = gb.l ** BigInt(gb.valsol);
  gb.ladicabs = gb.l ** BigInt(gb.valabs);
  return dnn;
}

/** PARI `makeLden` (galconj.c:284) */
function makeLden(L: bigint[], den: bigint, gb: GaloisBorne): bigint[] {
  const r: bigint[] = [0n];
  for (let i = 1; i < L.length; i++) r.push(bmod(L[i]! * den, gb.ladicsol));
  return r;
}

/** PARI `initlift` (galconj.c:289) */
function initlift(
  T: ZX,
  den: bigint,
  p: number,
  L: bigint[],
  Lden: bigint[],
  gb: GaloisBorne
): GaloisLift {
  const P = BigInt(p);
  let e = logint(gb.bornesol << BigInt(2 + BITS_IN_LONG), P) + 1;
  if (e < 2) e = 2;
  const Q = P ** BigInt(e);
  return {
    T,
    den: den === 1n ? 1n : den,
    p: P,
    L,
    Lden,
    e,
    Q,
    TQ: FpX_red(T, Q),
    gb,
  };
}

/**
 * PARI `poltopermtest` (galconj.c:311): check that `f` is (with high
 * probability) an automorphism and compute the permutation it induces.
 */
function poltopermtest(f: ZX, gl: GaloisLift, pf: Perm): boolean {
  const B = gl.gb.bornesol;
  for (let i = 0; i < f.length; i++) if (babs(f[i]!) > B) return false;
  const ll = gl.L.length;
  const fp: boolean[] = new Array(ll).fill(true);
  for (let i = 1; i < ll; i++) {
    const fx = FpX_eval(FpX_red(f, gl.gb.ladicsol), gl.L[i]!, gl.gb.ladicsol);
    let j = 1;
    for (; j < ll; j++)
      if (fp[j] && fx === gl.Lden[j]!) {
        pf[i] = j;
        fp[j] = false;
        break;
      }
    if (j === ll) return false;
  }
  return true;
}

/** PARI `galoisfrobeniustest` (galconj.c:337) */
function galoisfrobeniustest(aut: FpX, gl: GaloisLift, frob: Perm): boolean {
  let tlift = aut;
  if (gl.den !== 1n) tlift = FpX_Fp_mul(tlift, gl.den, gl.Q);
  const t = FpX_center(FpX_red(tlift, gl.Q), gl.Q, gl.Q >> 1n);
  return poltopermtest(t, gl, frob);
}

/** PARI `automorphismlift` (galconj.c:387) */
function automorphismlift(S: FpX, gl: GaloisLift): FpX {
  return ZpX_ZpXQ_liftroot(gl.T, S, gl.T, gl.p, gl.e);
}

/** PARI `galoisdolift` (galconj.c:393) */
function galoisdolift(gl: GaloisLift): FpX {
  const Tp = FpX_red(gl.T, gl.p);
  const S = FpX_Frobenius(Tp, gl.p);
  return automorphismlift(S, gl);
}

/** PARI `galoisdoliftn` (galconj.c:402) */
function galoisdoliftn(gl: GaloisLift, e: number): FpX {
  const Tp = FpX_red(gl.T, gl.p);
  const S = FpXQ_autpow(FpX_Frobenius(Tp, gl.p), e, Tp, gl.p);
  return automorphismlift(S, gl);
}

/** PARI `findpsi` (galconj.c:411) */
function findpsi(
  D: bigint,
  pstart: number,
  P: ZX,
  S: QPoly,
  o: number,
  out: { Tmod: ZX[]; psi: number[]; p: number }
): number {
  const n = ZX_degree(P);
  const g = n / o;
  const iter = new Forprime(pstart);
  for (;;) {
    const p = iter.next();
    const P_ = BigInt(p);
    if (bmod(D, P_) === 0n) continue;
    const F = FpX_factor_squarefree(FpX_red(P, P_), P_);
    if (F.length !== g) continue;
    const Fv: FpX[] = [[], ...F];
    const psi: number[] = new Array(g + 1).fill(0);
    const Sp = QPoly_to_FpX(S, P_);
    let gp = 0;
    let j = 1;
    for (; j <= g; j++) {
      const Fj = Fv[j]!;
      const Sj = FpX_rem(Sp, Fj, P_);
      const A = FpXQ_autpowers(FpX_Frobenius(Fj, P_), o, Fj, P_);
      let i = 1;
      for (; i <= o; i++)
        if (cmp_FpX(Sj, A[i + 1]!) === 0) {
          psi[j] = i;
          break;
        }
      if (i > o) break;
      if (gp === 0 && i === 1) gp = j;
    }
    if (gp && j > g) {
      if (gp !== 1) {
        const t = Fv[1]!;
        Fv[1] = Fv[gp]!;
        Fv[gp] = t;
        const u = psi[1]!;
        psi[1] = psi[gp]!;
        psi[gp] = u;
      }
      const inv = Fl_inv(psi[g]!, o);
      const Tpsi: number[] = [0];
      for (let i = 1; i <= g; i++) Tpsi.push(Fl_mul(psi[i]!, inv, o));
      out.Tmod = Fv;
      out.psi = Tpsi;
      out.p = p;
      return p;
    }
  }
}

/** PARI `inittestlift` (galconj.c:457) */
function inittestlift(plift: FpX, Tmod: ZX[], gl: GaloisLift): GaloisTestlift {
  const n = gl.L.length - 1;
  const g = Tmod.length - 1;
  const f = n / g;
  const bezoutcoeff = bezout_lift_fact(
    gl.T,
    Tmod.map((t, i) => (i === 0 ? [] : FpX_red(t, gl.p))),
    gl.p,
    gl.e
  );
  const pauto = FpXQ_autpowers(plift, f - 1, gl.TQ, gl.Q);
  return { n, g, f, bezoutcoeff, pauto, C: [], Cd: [] };
}

/* the `intheadlong` technique, galconj.c:471-511 */

function intheadlong(x: bigint, mod: bigint): bigint {
  return BigInt.asUintN(BITS_IN_LONG, (bmod(x, mod) << BigInt(BITS_IN_LONG)) / mod);
}

function polheadlong(P: FpX, n: number, mod: bigint): bigint {
  return P.length > n ? intheadlong(P[n]!, mod) : 0n;
}

function headlongisint(Z: bigint, n: number): boolean {
  return BigInt.asUintN(BITS_IN_LONG, -Z) <= BigInt(n);
}

function uadd(a: bigint, b: bigint): bigint {
  return BigInt.asUintN(BITS_IN_LONG, a + b);
}

/** PARI `frobeniusliftall` (galconj.c:513) */
function frobeniusliftall(
  sg: number[],
  el: number,
  gl: GaloisLift,
  gt: GaloisTestlift,
  frob: Perm
): number[] | null {
  const c = sg.length - 1;
  const n = gl.L.length - 1;
  const m = gt.g;
  const d = m / c;
  const ord = gt.f;
  let c_idx = gt.g - 1;
  const pf: number[] = new Array(m).fill(0);
  /* number of tests: m! / (c * (d!)^c) */
  let NN = 1n;
  for (let i = 2; i <= m; i++) NN *= BigInt(i);
  let df = 1n;
  for (let i = 2; i <= d; i++) df *= BigInt(i);
  NN /= BigInt(c) * df ** BigInt(c);
  if (NN > 10n ** 15n)
    throw new NotImplementedError(
      `SAGE_NOT_IMPLEMENTED: galconj frobeniusliftall would need ${NN} tests ` +
        '(PARI warns and returns a partial result; we refuse to guess)'
    );
  const C = gt.C;
  const Cd = gt.Cd;
  let v = FpX_rem(FpX_mul(gt.pauto[1 + (el % ord)]!, gt.bezoutcoeff[m]!, gl.Q), gl.TQ, gl.Q);
  if (gl.den !== 1n) v = FpX_Fp_mul(v, gl.den, gl.Q);
  const SG: number[] = [0];
  for (let i = 1; i <= c; i++) SG.push(((el * sg[i]!) % ord) + 1);
  const cache: bigint[] = new Array(m + 1).fill(0n);
  cache[m] = polheadlong(v, 1, gl.Q);
  const headcache = polheadlong(v, 2, gl.Q);
  for (let i = 1; i < m; i++) pf[i] = 1 + Math.floor(i / d);
  for (let iter = 0n; ; iter++) {
    for (let j = c_idx; j > 0; j--) {
      const h = SG[pf[j]!]!;
      if (!C[h]![j]) {
        let r = FpX_rem(FpX_mul(gt.pauto[h]!, gt.bezoutcoeff[j]!, gl.Q), gl.TQ, gl.Q);
        if (gl.den !== 1n) r = FpX_Fp_mul(r, gl.den, gl.Q);
        C[h]![j] = r;
        Cd[h]![j] = polheadlong(r, 1, gl.Q);
      }
      cache[j] = uadd(cache[j + 1]!, Cd[h]![j]!);
    }
    if (headlongisint(cache[1]!, n)) {
      let head = headcache;
      for (let j = 1; j < m; j++) head = uadd(head, polheadlong(C[SG[pf[j]!]!]![j]!, 2, gl.Q));
      if (headlongisint(head, n)) {
        let u: ZX = v;
        for (let j = 1; j < m; j++) u = ZX_add(u, C[SG[pf[j]!]!]![j]!);
        const uc = FpX_center(FpX_red(u, gl.Q), gl.Q, gl.Q >> 1n);
        if (poltopermtest(uc, gl, frob)) {
          const psi: number[] = new Array(m).fill(0);
          for (let i = 1; i < m; i++) psi[i] = pf[i]!;
          return psi;
        }
      }
    }
    if (iter >= NN - 1n) break;
    /* next multiset permutation (galconj.c:600-605) */
    let j = 2;
    for (; j < m && pf[j - 1]! >= pf[j]!; j++);
    for (let k = 1; k < j - k && pf[k] !== pf[j - k]; k++) {
      const t = pf[k]!;
      pf[k] = pf[j - k]!;
      pf[j - k] = t;
    }
    let k = j - 1;
    for (; pf[k]! >= pf[j]!; k--);
    const t = pf[j]!;
    pf[j] = pf[k]!;
    pf[k] = t;
    c_idx = j;
  }
  return null;
}

/** PARI `Vmatrix` (galconj.c:612) */
function Vmatrix(i: number, td: GaloisTest): bigint[][] {
  const n = td.L.length - 1;
  const W: bigint[][] = [[]];
  for (let j = 1; j <= n; j++) {
    const col: bigint[] = [0n];
    for (let k = 1; k <= n; k++) col.push(intheadlong(bmod(td.L[k]! * td.M[i]![j]!, td.ladic), td.ladic));
    W.push(col);
  }
  return W;
}

/** PARI `inittest` (galconj.c:621) */
function inittest(L: bigint[], M: bigint[][], borne: bigint, ladic: bigint): GaloisTest {
  const n = L.length - 1;
  const p: number[] = new Array(n + 1).fill(0);
  for (let i = 1; i <= n - 2; i++) p[i] = i + 2;
  p[n - 1] = 1;
  p[n] = 2;
  const td: GaloisTest = {
    order: p,
    borne,
    lborne: ladic - borne,
    ladic,
    L,
    M,
    PV: new Array(n + 1).fill(null),
  };
  td.PV[2] = Vmatrix(2, td);
  return td;
}

/** PARI `padicisint` (galconj.c:651) */
function padicisint(P: bigint, td: GaloisTest): boolean {
  const U = bmod(P, td.ladic);
  return U <= td.borne || U >= td.lborne;
}

/** PARI `galois_test_perm` (galconj.c:662) */
function galois_test_perm(td: GaloisTest, pf: Perm): boolean {
  const n = td.L.length - 1;
  let i = 1;
  for (; i < n; i++) {
    const ord = td.order[i]!;
    const PW = td.PV[ord];
    if (PW) {
      let head = PW[1]![pf[1]!]!;
      for (let j = 2; j <= n; j++) head = uadd(head, PW[j]![pf[j]!]!);
      if (!headlongisint(head, n)) break;
    } else {
      let V = 0n;
      for (let j = 1; j <= n; j++) V += td.M[ord]![j]! * td.L[pf[j]!]!;
      if (!padicisint(V, td)) {
        td.PV[ord] = Vmatrix(ord, td);
        break;
      }
    }
  }
  if (i === n) return true;
  if (i > 1) {
    const z = td.order[i]!;
    for (let j = i; j > 1; j--) td.order[j] = td.order[j - 1]!;
    td.order[1] = z;
  }
  return false;
}

/** PARI `testpermutation` (galconj.c:713) */
function testpermutation(
  F: Perm[],
  B: Perm[],
  x: number[],
  s: number,
  e: number,
  cut: number,
  td: GaloisTest
): Perm | null {
  const a = F.length - 1;
  const b = F[1]!.length - 1;
  const c = B.length - 1;
  const d = B[1]!.length - 1;
  const n = a * b;
  s = (b + s) % b;
  const pf: Perm = new Array(n + 1).fill(0);
  const ar: bigint[] = new Array(a + 2).fill(0n);
  const G: Perm[] = new Array(a + 1).fill(null).map(() => []);
  const W = td.PV[td.order[n]!]!;
  {
    let i = 1;
    let j = 1;
    for (let cx = 1; cx <= a; cx++, i++) {
      G[cx] = F[B[j]![i]!]!;
      if (i === d) {
        i = 0;
        j++;
      }
    }
  }
  let NN = BigInt(b) ** BigInt(c * (d - Math.floor(d / e)));
  NN /= BigInt(cut);
  if (NN > 10n ** 14n)
    throw new NotImplementedError(
      `SAGE_NOT_IMPLEMENTED: galconj testpermutation would need ${NN} tests ` +
        '(PARI warns and returns a partial result; we refuse to guess)'
    );
  let start = 0;
  for (let l1 = 0n; l1 < NN; l1++) {
    let i: number;
    if (start) {
      let j = e;
      i = 1;
      while (i < a) {
        if (++x[i]! !== b) break;
        x[i++] = 0;
        if (i === j) {
          i++;
          j += e;
        }
      }
    } else {
      start = 1;
      i = a - 1;
    }
    /* intheadlong test */
    let p5 = ((i + 1) % d) - 1;
    for (let p1 = i + 1; p1 >= 1; p1--, p5--) {
      let V = 0n;
      let p6: number;
      if (p5 === -1) {
        p5 = d - 1;
        p6 = p1 + 1 - d;
      } else p6 = p1 + 1;
      const G1 = G[p1]!;
      const G6 = G[p6]!;
      const p4 = p5 ? x[p1 - 1]! : 0;
      let p3 = 1 + x[p1]!;
      for (let p2 = 1 + p4; p2 <= b; p2++) {
        V = uadd(V, W[G6[p3]!]![G1[p2]!]!);
        p3 += s;
        if (p3 > b) p3 -= b;
      }
      p3 = 1 + x[p1]! - s;
      if (p3 <= 0) p3 += b;
      for (let p2 = p4; p2 >= 1; p2--) {
        V = uadd(V, W[G6[p3]!]![G1[p2]!]!);
        p3 -= s;
        if (p3 <= 0) p3 += b;
      }
      ar[p1] = uadd(ar[p1 + 1]!, V);
    }
    if (!headlongisint(ar[1]!, n)) continue;
    /* full computation */
    for (let p1 = 1, p5b = d; p1 <= a; p1++, p5b++) {
      let p4: number;
      if (p5b === d) {
        p5b = 0;
        p4 = 0;
      } else p4 = x[p1 - 1]!;
      const p6 = p5b === d - 1 ? p1 + 1 - d : p1 + 1;
      let p3 = 1 + x[p1]!;
      for (let p2 = 1 + p4; p2 <= b; p2++) {
        pf[G[p1]![p2]!] = G[p6]![p3]!;
        p3 += s;
        if (p3 > b) p3 -= b;
      }
      p3 = 1 + x[p1]! - s;
      if (p3 <= 0) p3 += b;
      for (let p2 = p4; p2 >= 1; p2--) {
        pf[G[p1]![p2]!] = G[p6]![p3]!;
        p3 -= s;
        if (p3 <= 0) p3 += b;
      }
    }
    if (galois_test_perm(td, pf)) return pf;
  }
  return null;
}

/**
 * PARI `listznstarelts` (galconj.c:830): the subgroups of `(Z/mZ)^*` whose
 * order divides `o`, as the (sorted) lists of their elements, by increasing
 * order.
 */
export function listznstarelts(m: number, o: number): number[][] {
  if (m === 2) return [[0, 1]];
  const units: number[] = [];
  for (let i = 1; i < m; i++) if (ugcd(i, m) === 1) units.push(i);
  const phi = units.length;
  o = ugcd(o, phi);
  /* enumerate all subgroups by closure of subsets of elements, keeping those
   * whose order divides o */
  const seen = new Map<string, number[]>();
  const closure = (gens: number[]): number[] => {
    const S = new Set<number>([1]);
    let added = true;
    while (added) {
      added = false;
      for (const g of gens)
        for (const s of Array.from(S)) {
          const t = (g * s) % m;
          if (!S.has(t)) {
            S.add(t);
            added = true;
          }
        }
    }
    return Array.from(S).sort((a, b) => a - b);
  };
  const queue: number[][] = [[]];
  seen.set('', [1]);
  const groups: number[][] = [[1]];
  while (queue.length) {
    const g = queue.pop()!;
    const cur = closure(g);
    for (const u of units) {
      if (cur.includes(u)) continue;
      const ng = [...g, u];
      const cl = closure(ng);
      if (cl.length > o) continue;
      const key = cl.join(',');
      if (seen.has(key)) continue;
      seen.set(key, cl);
      groups.push(cl);
      queue.push(ng);
    }
  }
  const res = groups.filter((g) => o % g.length === 0);
  res.sort((a, b) => a.length - b.length || vecsmall_lexcmp([0, ...a], [0, ...b]));
  return res.map((g) => [0, ...g]);
}

/* ---- symmetric polynomials (galconj.c:849-1035) ---- */

/** a sympol `sum_i v[i] * s_{w[i]}`; both 1-indexed (galconj.c:849) */
export interface SymPol {
  v: number[];
  w: number[];
}

/** PARI `Flm_newtonsum` (galconj.c:855) */
function Flm_newtonsum(M: bigint[][], e: number, p: bigint): bigint[] {
  const NS: bigint[] = [0n];
  for (let i = 1; i < M.length; i++) {
    let s = 0n;
    const Mi = M[i]!;
    for (let j = 1; j < Mi.length; j++) s = (s + bmod(Mi[j]!, p) ** BigInt(e)) % p;
    NS.push(s);
  }
  return NS;
}

/** PARI `Flv_sympol_eval` (galconj.c:871) */
function Flv_sympol_eval(v: number[], NS: bigint[][], p: bigint): bigint[] {
  const l = v.length;
  const n = NS[1]!.length;
  const S: bigint[] = new Array(n).fill(0n);
  for (let i = 1; i < l; i++) {
    if (!v[i]) continue;
    for (let j = 1; j < n; j++) S[j] = (S[j]! + BigInt(v[i]!) * NS[i]![j]!) % p;
  }
  return S;
}

/** PARI `sympol_eval_newtonsum` (galconj.c:882) */
function sympol_eval_newtonsum(e: number, O: bigint[][], mod: bigint): bigint[] {
  const PL: bigint[] = [0n];
  for (let i = 1; i < O.length; i++) {
    let s = 0n;
    for (let j = 1; j < O[i]!.length; j++) {
      let t = 1n;
      const b = bmod(O[i]![j]!, mod);
      for (let k = 0; k < e; k++) t = (t * b) % mod;
      s = (s + t) % mod;
    }
    PL.push(s);
  }
  return PL;
}

/** PARI `sympol_eval` (galconj.c:897) */
export function sympol_eval(sym: SymPol, O: bigint[][], mod: bigint): bigint[] {
  const n = O.length;
  const S: bigint[] = new Array(n).fill(0n);
  for (let i = 1; i < sym.v.length; i++) {
    if (!sym.v[i]) continue;
    const N = sympol_eval_newtonsum(sym.w[i]!, O, mod);
    for (let j = 1; j < n; j++) S[j] = bmod(S[j]! + BigInt(sym.v[i]!) * N[j]!, mod);
  }
  return S;
}

/** PARI `sympol_aut_evalmod` (galconj.c:912) */
function sympol_aut_evalmod(sym: SymPol, g: number, sigma: QPoly, Tp: FpX, p: bigint): FpX {
  const sig = QPoly_to_FpX(sigma, p);
  let f: FpX = pol_xn(1);
  let s: FpX = [];
  for (let i = 1; i <= g; i++) {
    if (i > 1) f = FpX_FpXQ_eval(f, sig, Tp, p);
    for (let j = 1; j < sym.v.length; j++) {
      if (!sym.v[j]) continue;
      s = FpX_add(s, FpX_Fp_mul(FpXQ_powBig(f, BigInt(sym.w[j]!), Tp, p), BigInt(sym.v[j]!), p), p);
    }
  }
  return s;
}

/** PARI `fixedfieldfactmod` (galconj.c:934) */
function fixedfieldfactmod(Sp: FpX, p: bigint, Tmod: ZX[]): FpX[] {
  const F: FpX[] = [[]];
  for (let i = 1; i < Tmod.length; i++) {
    const Ti = FpX_red(Tmod[i]!, p);
    F.push(FpXQ_minpoly(FpX_rem(Sp, Ti, p), Ti, p));
  }
  return F;
}

function vecsmall_is1to1(L: bigint[]): boolean {
  const s = new Set<bigint>();
  for (let i = 1; i < L.length; i++) {
    if (s.has(L[i]!)) return false;
    s.add(L[i]!);
  }
  return true;
}

/** PARI `fixedfieldsurmer` (galconj.c:947) */
function fixedfieldsurmer(l: bigint, NS: bigint[][], W: number[]): SymPol | null {
  const step = 3;
  const n = W.length - 1;
  const m = 1 << ((n - 1) << 1);
  const sym: number[] = new Array(n + 1).fill(0);
  for (let j = 1; j < n; j++) sym[j] = step;
  sym[n] = 0;
  for (let i = 0; i < m; i++) {
    let j = 1;
    for (; sym[j] === step; j++) sym[j] = 0;
    sym[j]!++;
    const L = Flv_sympol_eval(sym, NS, l);
    if (!vecsmall_is1to1(L)) continue;
    return { v: sym.slice(), w: W.slice() };
  }
  return null;
}

/** PARI `sympol_is1to1_lg` (galconj.c:971) */
function sympol_is1to1_lg(NS: bigint[][], n: number): boolean {
  const l = NS[1]!.length;
  for (let i = 1; i < l; i++)
    for (let j = i + 1; j < l; j++) {
      let k = 1;
      for (; k < n; k++) if (NS[k]![j]! !== NS[k]![i]!) break;
      if (k >= n) return false;
    }
  return true;
}

/** PARI `fixedfieldsympol` (galconj.c:989) */
export function fixedfieldsympol(O: bigint[][], l: bigint): SymPol {
  const n = (BITS_IN_LONG >> 1) - 1;
  const NS: bigint[][] = [[]];
  const W: number[] = new Array(n + 1).fill(0);
  let sym: SymPol | null = null;
  let e = 1;
  const Ol: bigint[][] = O.map((c, i) => (i === 0 ? [] : c.map((x, j) => (j === 0 ? 0n : bmod(x, l)))));
  for (let i = 1; !sym && i <= n; i++) {
    let L = Flm_newtonsum(Ol, e++, l);
    if (Ol.length > 2) {
      const isconst = (v: bigint[]) => {
        for (let k = 2; k < v.length; k++) if (v[k] !== v[1]) return false;
        return true;
      };
      while (isconst(L)) L = Flm_newtonsum(Ol, e++, l);
    }
    W[i] = e - 1;
    NS[i] = L;
    if (sympol_is1to1_lg(NS, i + 1)) sym = fixedfieldsurmer(l, NS, W.slice(0, i + 1));
  }
  if (!sym) throw new PariBugError('fixedfieldsympol [p too small]');
  return sym;
}

/** PARI `fixedfieldorbits` (galconj.c:1015) */
export function fixedfieldorbits(O: Perm[], L: bigint[]): bigint[][] {
  const S: bigint[][] = [[]];
  for (let i = 1; i < O.length; i++) S.push(vecpermute(L, O[i]!));
  return S;
}

/** PARI `fixedfieldinclusion` (galconj.c:1024) */
function fixedfieldinclusion(O: Perm[], PL: bigint[]): bigint[] {
  const f = O.length - 1;
  const g = O[1]!.length - 1;
  const S: bigint[] = new Array(f * g + 1).fill(0n);
  for (let i = 1; i <= f; i++) {
    const Oi = O[i]!;
    for (let j = 1; j <= g; j++) S[Oi[j]!] = PL[i]!;
  }
  return S;
}

/** PARI `vectopol` (galconj.c:1038) */
function vectopol(v: bigint[], M: bigint[][], den: bigint, mod: bigint, mod2: bigint): QPoly {
  const l = v.length;
  const num: ZX = [];
  for (let k = 1; k < l; k++) {
    let s = 0n;
    for (let i = 1; i < l; i++) s += M[k]![i]! * v[i]!;
    num.push(Fp_center(s, mod, mod2));
  }
  return QPoly_normalize({ num: ZX_renormalize(num), den });
}

/** PARI `permtopol` (galconj.c:1050) */
export function permtopol(
  p: Perm,
  L: bigint[],
  M: bigint[][],
  den: bigint,
  mod: bigint,
  mod2: bigint
): QPoly {
  if (p.length !== L.length) throw new PariTypeError('permtopol [permutation]', 'p');
  return vectopol(vecpermute(L, p), M, den, mod, mod2);
}

/* ---- galoisanalysis (galconj.c:1075-1235) ---- */

function notgalois(p: number, ga: GaloisAnalysis): void {
  ga.p = p;
  ga.deg = 0;
}

/** PARI `init_group` (galconj.c:1084) */
function init_group(n: number, np: number, Fp: number[], Fe: number[]): { group: number; order: number } {
  const prim_nonwss_orders = [48, 56, 60, 72, 75, 80, 196, 200, 216];
  let phi_order = 1;
  let order = 1;
  let group = 0;
  for (const o of prim_nonwss_orders)
    if (n % o === 0) {
      group |= ga_non_wss;
      break;
    }
  if (np === 2 && Fp[2] === 3 && Fe[2] === 1 && Fe[1]! > 2) group |= ga_ext_2;
  for (let i = np; i > 0; i--) {
    const p = Fp[i]!;
    if (phi_order % p === 0) {
      group |= ga_all_normal;
      break;
    }
    order *= p;
    phi_order *= p - 1;
    if (Fe[i]! > 1) break;
  }
  if (uisprimepower(n) || n === 135) group |= ga_all_nilpotent;
  if (n <= 104) group |= ga_easy;
  return { group, order };
}

/** PARI `improves` (galconj.c:1109) */
function improves(
  a: number,
  b: number,
  plift: number,
  p: number,
  n: number,
  karma: { v: number }
): boolean {
  if (!plift || a > b) {
    karma.v = ugcd(p - 1, n);
    return true;
  }
  if (a === b) {
    const k = ugcd(p - 1, n);
    if (k > karma.v) {
      karma.v = k;
      return true;
    }
  }
  return false;
}

/** PARI `galoisanalysis` (galconj.c:1121); returns false if not Galois or not WSS */
export function galoisanalysis(
  T: ZX,
  ga: GaloisAnalysis,
  calcul_l: boolean,
  bad: bigint | null
): boolean {
  const n = ZX_degree(T);
  const O: number[] = new Array(n + 1).fill(0);
  const { Fp, Fe, Fpe } = factoru_pow(n);
  const np = Fp.length - 1;
  let { group, order } = init_group(n, np, Fp, Fe);
  let deg = Fp[np]!;
  let plift = 0;
  let nbtest = 0;
  const karma = { v: 0 };
  const nbmax = 8 + (n >> 1);
  const S = new Forprime(n * Math.max(expu(n) - 3, 2));
  for (;;) {
    const cont =
      !plift ||
      (nbtest < nbmax && (nbtest <= 8 || order < n >> 1)) ||
      ((n === 24 || n === 36) && O[6] === 0 && O[4] === 0) ||
      (group & ga_non_wss && order === Fp[np]!);
    if (!cont) break;
    if (group & ga_non_wss && nbtest >= 3 * nbmax) break;
    nbtest++;
    const p = S.next();
    const P = BigInt(p);
    if (bad !== null && bmod(bad, P) === 0n) continue;
    const Tp = FpX_red(T, P);
    if (!FpX_is_squarefree(Tp, P)) {
      if (!--nbtest) nbtest = 1;
      continue;
    }
    const { D, nb: d } = FpX_nbfact_by_degree(Tp, P);
    const o = n / d;
    if (!Number.isInteger(o) || D[o] !== d) {
      notgalois(p, ga);
      return false;
    }
    if (!O[o]) O[o] = p;
    let norm_o = 1;
    ga_end: {
      if (o % deg) break ga_end;
      if (group & ga_all_normal && o < order) break ga_end;
      if (o * Fp[1]! >= n) norm_o = o;
      else {
        for (let i = np; i > 0; i--) {
          if (o % Fpe[i]!) break;
          norm_o *= Fpe[i]!;
        }
      }
      if (norm_o !== 1) {
        if (!(group & ga_all_normal) || o > order) karma.v = ugcd(p - 1, n);
        else if (!improves(norm_o, deg, plift, p, n, karma)) break ga_end;
        deg = norm_o;
        group |= ga_all_normal;
      } else if (group & ga_all_normal) break ga_end;
      else if (!improves(o, order, plift, p, n, karma)) break ga_end;
      order = o;
      plift = p;
    }
  }
  ga.p = plift;
  if (!plift || (group & ga_non_wss && order === Fp[np]!)) return false;
  const linf = 2 * n * usqrt(n);
  if (calcul_l && O[1]! <= linf) {
    const S2 = new Forprime(linf + 1);
    for (;;) {
      const p = S2.next();
      const P = BigInt(p);
      const Tp = FpX_red(T, P);
      const nb = FpX_nbroots(Tp, P);
      if (nb === n) {
        O[1] = p;
        break;
      }
      if (nb && FpX_is_squarefree(Tp, P)) {
        notgalois(p, ga);
        return false;
      }
    }
  }
  ga.group = group;
  ga.deg = deg;
  ga.ord = order;
  ga.l = O[1]!;
  ga.p4 = n >= 4 ? O[4]! : 0;
  return true;
}

/* ---- the lift of the Frobenius (galconj.c:1878-2110) ---- */

/** PARI `galoisfindgroups` (galconj.c:1879) */
function galoisfindgroups(lo: number[][], sg: number[], f: number): number[][] {
  const V: number[][] = [];
  for (const loi of lo) {
    const W: number[] = [0];
    for (let k = 1; k < loi.length; k++) W.push(loi[k]! % f);
    if (zv_equal(vecsmall_uniq(W), sg)) V.push(loi);
  }
  return V;
}

/** PARI `galoismakepsi` (galconj.c:1897) */
function galoismakepsi(g: number, sg: number[], pf: number[]): number[] {
  const psi: number[] = new Array(g + 1).fill(0);
  for (let i = 1; i < g; i++) psi[i] = sg[pf[i]!]!;
  psi[g] = sg[1]!;
  return psi;
}

/** PARI `galoisfrobeniuslift_nilp` (galconj.c:1906) */
function galoisfrobeniuslift_nilp(
  T: ZX,
  den: bigint,
  L: bigint[],
  Lden: bigint[],
  gf: GaloisFrobenius,
  gb: GaloisBorne
): Perm | null {
  let deg = 1;
  const g = gf.Tmod.length - 1;
  const res: Perm = new Array(L.length).fill(0);
  gf.psi = new Array(g + 1).fill(1);
  gf.psi[0] = 0;
  const gl = initlift(T, den, gf.p, L, Lden, gb);
  let aut = galoisdolift(gl);
  if (galoisfrobeniustest(aut, gl, res)) {
    gf.deg = gf.fp;
    return res;
  }
  const F = factoru_small(gf.fp);
  const frob: Perm = new Array(L.length).fill(0);
  for (let k = F.length - 1; k >= 0; k--) {
    const [Fp_, Fe_] = F[k]!;
    let fres: Perm | null = null;
    let el = gf.fp;
    let dg = 1;
    let dgf = 1;
    for (let e = 1; e <= Fe_; e++) {
      dg *= Fp_;
      el /= Fp_;
      if (el === 1) break;
      aut = galoisdoliftn(gl, el);
      if (!galoisfrobeniustest(aut, gl, frob)) break;
      dgf = dg;
      fres = frob.slice();
    }
    if (dgf === 1) continue;
    const pr = deg * dgf;
    if (deg === 1) {
      for (let i = 1; i < res.length; i++) res[i] = fres![i]!;
    } else {
      const cp = perm_mul(res, fres!);
      for (let i = 1; i < res.length; i++) res[i] = cp[i]!;
    }
    deg = pr;
  }
  if (deg === 1) return null;
  gf.deg = deg;
  return res;
}

/** PARI `galoisfrobeniuslift` (galconj.c:1966) */
function galoisfrobeniuslift(
  T: ZX,
  den: bigint,
  L: bigint[],
  Lden: bigint[],
  gf: GaloisFrobenius,
  gb: GaloisBorne
): Perm | null {
  const n = L.length - 1;
  let deg = 1;
  const g = gf.Tmod.length - 1;
  const res: Perm = new Array(L.length).fill(0);
  gf.psi = new Array(g + 1).fill(1);
  gf.psi[0] = 0;
  const gl = initlift(T, den, gf.p, L, Lden, gb);
  const aut = galoisdolift(gl);
  if (galoisfrobeniustest(aut, gl, res)) {
    gf.deg = gf.fp;
    return res;
  }
  const gt = inittestlift(aut, gf.Tmod, gl);
  gt.C = new Array(gf.fp + 1).fill(null).map(() => new Array(gt.g + 1).fill(null));
  gt.Cd = new Array(gf.fp + 1).fill(null).map(() => new Array(gt.g + 1).fill(0n));

  const F = factoru_small(gf.fp);
  const frob: Perm = new Array(L.length).fill(0);
  for (let k = F.length - 1; k >= 0; k--) {
    const [Fp_, Fe_] = F[k]!;
    let psi: number[] | null = null;
    let fres: Perm | null = null;
    let sg: number[] = [0, 1];
    let el = gf.fp;
    let dg = 1;
    let dgf = 1;
    for (let e = 1; e <= Fe_; e++) {
      dg *= Fp_;
      el /= Fp_;
      if (galoisfrobeniustest(gt.pauto[el + 1]!, gl, frob)) {
        psi = new Array(g + 1).fill(1);
        psi[0] = 0;
        dgf = dg;
        fres = frob.slice();
        continue;
      }
      let lo = listznstarelts(dg, n / gf.fp);
      if (e !== 1) lo = galoisfindgroups(lo, sg, dgf);
      let l = 0;
      for (; l < lo.length; l++) {
        if (lo[l]!.length <= 2) continue;
        const pf = frobeniusliftall(lo[l]!, el, gl, gt, frob);
        if (pf) {
          sg = lo[l]!.slice();
          psi = galoismakepsi(g, sg, pf);
          dgf = dg;
          fres = frob.slice();
          break;
        }
      }
      if (l === lo.length) break;
    }
    if (dgf === 1) continue;
    const pr = deg * dgf;
    if (deg === 1) {
      for (let i = 1; i < res.length; i++) res[i] = fres![i]!;
      for (let i = 1; i < psi!.length; i++) gf.psi[i] = psi![i]!;
    } else {
      const cp = perm_mul(res, fres!);
      for (let i = 1; i < res.length; i++) res[i] = cp[i]!;
      for (let i = 1; i < psi!.length; i++) gf.psi[i] = (dgf * gf.psi[i]! + deg * psi![i]!) % pr;
    }
    deg = pr;
  }
  if (deg === 1) return null;
  /* Normalize result so that psi[g] = 1 */
  const im = Fl_inv(gf.psi[g]!, deg);
  const cp = perm_powu(res, im);
  for (let i = 1; i < res.length; i++) res[i] = cp[i]!;
  for (let i = 1; i < gf.psi.length; i++) gf.psi[i] = Fl_mul(im, gf.psi[i]!, deg);
  gf.deg = deg;
  return res;
}


/**
 * PARI `a4galoisgen` (galconj.c:1237): the dedicated generator search for
 * `A4` (degree 12).  Verbatim transcription, including the two hard-coded
 * enumerations of 10395 and 60 pairings and the final 2x4 search for the
 * order-3 generator.
 */
function a4galoisgen(td: GaloisTest): GaloisGens | null {
  const n = 12;
  const pft: Perm = new Array(n + 1).fill(0);
  const pfu: Perm = new Array(n + 1).fill(0);
  const pfv: Perm = new Array(n + 1).fill(0);
  const ar: bigint[] = new Array(5).fill(0n);
  const mt = td.PV[td.order[n]!]!;
  /* t and u are 0-indexed (upstream shifts the t_VECSMALL by one) */
  const t: number[] = new Array(n);
  for (let k = 0; k < n; k++) t[k] = k + 1;
  const u: number[] = new Array(n).fill(0);
  /* MT(i,j) = mt(i,j) + mt(j,i), symmetric; 1-indexed */
  const MT: bigint[][] = new Array(n + 1).fill(null).map(() => new Array(n + 1).fill(0n));
  for (let j = 1; j <= n; j++)
    for (let i = 1; i < j; i++) {
      const v = uadd(mt[j]![i]!, mt[i]![j]!);
      MT[i]![j] = v;
      MT[j]![i] = v;
    }
  /* n = 2k = 12; N = (2k)!/(k! 2^k) = 10395 */
  let N = 10395;
  ar[4] = MT[11]![12]!;
  ar[3] = uadd(ar[4]!, MT[9]![10]!);
  ar[2] = uadd(ar[3]!, MT[7]![8]!);
  ar[1] = uadd(ar[2]!, MT[5]![6]!);
  let i = 0;
  for (i = 0; i < N; i++) {
    if (i) {
      let a: number;
      let x = i;
      let y = 1;
      do {
        y += 2;
        a = x % y;
        x = Math.floor(x / y);
      } while (!a);
      switch (y) {
        case 3:
          [t[2], t[2 - a]] = [t[2 - a]!, t[2]!];
          break;
        case 5: {
          const z = t[0]!;
          t[0] = t[2]!;
          t[2] = t[1]!;
          t[1] = z;
          [t[4], t[4 - a]] = [t[4 - a]!, t[4]!];
          ar[1] = uadd(ar[2]!, MT[t[4]!]![t[5]!]!);
          break;
        }
        case 7: {
          const z = t[0]!;
          t[0] = t[4]!;
          t[4] = t[3]!;
          t[3] = t[1]!;
          t[1] = t[2]!;
          t[2] = z;
          [t[6], t[6 - a]] = [t[6 - a]!, t[6]!];
          ar[2] = uadd(ar[3]!, MT[t[6]!]![t[7]!]!);
          ar[1] = uadd(ar[2]!, MT[t[4]!]![t[5]!]!);
          break;
        }
        case 9: {
          const z = t[0]!;
          t[0] = t[6]!;
          t[6] = t[5]!;
          t[5] = t[3]!;
          t[3] = z;
          [t[1], t[4]] = [t[4]!, t[1]!];
          [t[8], t[8 - a]] = [t[8 - a]!, t[8]!];
          ar[3] = uadd(ar[4]!, MT[t[8]!]![t[9]!]!);
          ar[2] = uadd(ar[3]!, MT[t[6]!]![t[7]!]!);
          ar[1] = uadd(ar[2]!, MT[t[4]!]![t[5]!]!);
          break;
        }
        case 11: {
          const z = t[0]!;
          t[0] = t[8]!;
          t[8] = t[7]!;
          t[7] = t[5]!;
          t[5] = t[1]!;
          t[1] = t[6]!;
          t[6] = t[3]!;
          t[3] = t[2]!;
          t[2] = t[4]!;
          t[4] = z;
          [t[10], t[10 - a]] = [t[10 - a]!, t[10]!];
          ar[4] = MT[t[10]!]![t[11]!]!;
          ar[3] = uadd(ar[4]!, MT[t[8]!]![t[9]!]!);
          ar[2] = uadd(ar[3]!, MT[t[6]!]![t[7]!]!);
          ar[1] = uadd(ar[2]!, MT[t[4]!]![t[5]!]!);
          break;
        }
      }
    }
    const g = uadd(uadd(ar[1]!, MT[t[0]!]![t[1]!]!), MT[t[2]!]![t[3]!]!);
    if (headlongisint(g, n)) {
      for (let k = 0; k < n; k += 2) {
        pft[t[k]!] = t[k + 1]!;
        pft[t[k + 1]!] = t[k]!;
      }
      if (galois_test_perm(td, pft)) break;
    }
  }
  if (i === N) return null;
  /* N = (k)!/(k/2)! / 2 = 60 */
  N = 60;
  for (let k = 0; k < n; k += 4) {
    u[k + 3] = t[k + 3]!;
    u[k + 2] = t[k + 1]!;
    u[k + 1] = t[k + 2]!;
    u[k] = t[k]!;
  }
  for (i = 0; i < N; i++) {
    let g = 0n;
    if (i) {
      let a: number;
      let x = i;
      let y = -2;
      do {
        y += 4;
        a = x % y;
        x = Math.floor(x / y);
      } while (!a);
      [u[0], u[2]] = [u[2]!, u[0]!];
      switch (y) {
        case 2:
          break;
        case 6:
          [u[4], u[6]] = [u[6]!, u[4]!];
          if (!(a & 1)) {
            a = 4 - (a >> 1);
            [u[6], u[a]] = [u[a]!, u[6]!];
            [u[4], u[a - 2]] = [u[a - 2]!, u[4]!];
          }
          break;
        case 10: {
          const z = u[6]!;
          u[6] = u[3]!;
          u[3] = u[2]!;
          u[2] = u[4]!;
          u[4] = u[1]!;
          u[1] = u[0]!;
          u[0] = z;
          if (a >= 3) a += 2;
          a = 8 - a;
          [u[10], u[a]] = [u[a]!, u[10]!];
          [u[8], u[a - 2]] = [u[a - 2]!, u[8]!];
          break;
        }
      }
    }
    for (let k = 0; k < n; k += 2) g = uadd(g, MT[u[k]!]![u[k + 1]!]!);
    if (headlongisint(g, n)) {
      for (let k = 0; k < n; k += 2) {
        pfu[u[k]!] = u[k + 1]!;
        pfu[u[k + 1]!] = u[k]!;
      }
      if (galois_test_perm(td, pfu)) break;
    }
  }
  if (i === N) return null;
  const orb = [pft, pfu];
  const O = vecperm_orbits(orb, 12);
  const O1 = O[1]!;
  const O2 = O[2]!;
  const O3 = O[3]!;
  for (let j = 0; j < 2; j++) {
    pfv[O1[1]!] = O2[1]!;
    pfv[O1[2]!] = O2[3 + j]!;
    pfv[O1[3]!] = O2[4 - (j << 1)]!;
    pfv[O1[4]!] = O2[2 + j]!;
    for (let i2 = 0; i2 < 4; i2++) {
      let g = 0n;
      switch (i2) {
        case 0:
          break;
        case 1:
          [O3[1], O3[2]] = [O3[2]!, O3[1]!];
          [O3[3], O3[4]] = [O3[4]!, O3[3]!];
          break;
        case 2:
          [O3[1], O3[4]] = [O3[4]!, O3[1]!];
          [O3[2], O3[3]] = [O3[3]!, O3[2]!];
          break;
        case 3:
          [O3[1], O3[2]] = [O3[2]!, O3[1]!];
          [O3[3], O3[4]] = [O3[4]!, O3[3]!];
          break;
      }
      pfv[O2[1]!] = O3[1]!;
      pfv[O2[3 + j]!] = O3[4 - j]!;
      pfv[O2[4 - (j << 1)]!] = O3[2 + (j << 1)]!;
      pfv[O2[2 + j]!] = O3[3 - j]!;
      pfv[O3[1]!] = O1[1]!;
      pfv[O3[4 - j]!] = O1[2]!;
      pfv[O3[2 + (j << 1)]!] = O1[3]!;
      pfv[O3[3 - j]!] = O1[4]!;
      for (let k = 1; k <= n; k++) g = uadd(g, mt[k]![pfv[k]!]!);
      if (headlongisint(g, n) && galois_test_perm(td, pfv))
        return { gen: [[], pft, pfu, pfv], orders: [0, 2, 2, 3] };
    }
  }
  return null;
}

/** PARI `galoisfindfrobenius` (galconj.c:2057); returns null if not Galois */
function galoisfindfrobenius(
  T: ZX,
  L: bigint[],
  den: bigint,
  bad: bigint | null,
  gf: GaloisFrobenius,
  gb: GaloisBorne,
  ga: GaloisAnalysis
): Perm | null {
  let Try = 0;
  const n = ZX_degree(T);
  /* Upstream can loop over primes for ever here when the group is one of the
   * three shapes it handles by a dedicated routine (a4/s4/f36galoisgen), since
   * it only ever reaches this point after that routine has failed.  We port
   * a4galoisgen but not s4/f36galoisgen, so we extend upstream's own
   * `Try > 3n/2` guard (galconj.c:2101) to those shapes. */
  const isSpecial = n === 12 && ga.ord === 3 && !ga.p4;
  const deg = (gf.deg = ga.deg);
  let gmask = ga.group & ga_ext_2 ? 3 : 1;
  const Lden = makeLden(L, den, gb);
  const is_nilpotent = ga.group & ga_all_nilpotent;
  const S = new Forprime(ga.p);
  for (;;) {
    gf.p = S.next();
    const P = BigInt(gf.p);
    const Tp = FpX_red(T, P);
    if (!FpX_is_squarefree(Tp, P)) continue;
    if (bad !== null && bmod(bad, P) === 0n) continue;
    const Ti = FpX_factor_squarefree(Tp, P);
    const nb = Ti.length;
    const d = FpX_degree(Ti[0]!);
    if (nb > 1 && FpX_degree(Ti[nb - 1]!) !== d) return null;
    if (((gmask & 1) === 0 || d % deg) && ((gmask & 2) === 0 || d % 2 === 1)) continue;
    gf.fp = d;
    gf.Tmod = [[], ...Ti];
    const frob = is_nilpotent
      ? galoisfrobeniuslift_nilp(T, den, L, Lden, gf, gb)
      : galoisfrobeniuslift(T, den, L, Lden, gf, gb);
    if (frob) return frob;
    if (is_nilpotent) continue;
    if (ga.group & ga_all_normal && d % deg === 0) gmask &= ~1;
    if (!gmask) return null;
    if ((ga.group & ga_non_wss || isSpecial) && ++Try > (3 * n) >> 1) return null;
  }
}

/** PARI `get_image` (galconj.c:2114) */
function get_image(tau: FpX, P: FpX, Pmod: FpX[], p: bigint): number {
  const gp = Pmod.length - 1;
  let t = FpX_FpXQ_eval(Pmod[gp]!, tau, P, p);
  t = FpX_normalize(FpX_gcd(P, t, p), p);
  for (let g = 1; g <= gp; g++) if (cmp_FpX(t, Pmod[g]!) === 0) return g;
  return 0;
}

/** PARI `galoisgenfixedfield` (galconj.c:2136) */
function galoisgenfixedfield(
  Pmod: FpX[],
  PL: bigint[],
  P: ZX,
  ip: bigint,
  bad: bigint | null,
  gb: GaloisBorne
): { PG: GaloisGens; Pg: number[] } | null {
  const Pp = FpX_red(P, ip);
  if (ZX_degree(P) === 2 && bad === null) {
    const PG: GaloisGens = { gen: [[], [0, 2, 1]], orders: [0, 2] };
    /* tau = -x - P[1] */
    const tau: ZX = ZX_renormalize([-(P[1] ?? 0n), -1n]);
    const g = get_image(FpX_red(tau, ip), Pp, Pmod, ip);
    if (!g) return null;
    return { PG, Pg: [0, g] };
  }
  const Pga: GaloisAnalysis = { p: 0, deg: 0, ord: 0, l: 0, p4: 0, group: 0 };
  if (!galoisanalysis(P, Pga, false, null)) return null;
  if (bad !== null) Pga.group &= ~ga_easy;
  const Pgb: GaloisBorne = {
    l: gb.l,
    valsol: 0,
    valabs: 0,
    bornesol: 0n,
    ladicsol: 0n,
    ladicabs: 0n,
    dis: 0n,
  };
  const Pden = galoisborne(P, null, Pgb, ZX_degree(P));
  if (Pgb.valabs > gb.valabs) PL = ZpX_liftroots(P, PL, gb.l, Pgb.valabs);
  else if (Pgb.valabs < gb.valabs)
    PL = PL.map((x, i) => (i === 0 ? 0n : bmod(x, Pgb.ladicabs)));
  const PM = FpV_invVandermonde(PL, Pden, Pgb.ladicabs);
  const PG = galoisgen(
    P,
    PL,
    PM,
    Pden,
    bad !== null ? lcmBig(Pgb.dis, bad) : null,
    Pgb,
    Pga
  );
  if (!PG) return null;
  const lP = PG.gen.length;
  const mod = Pgb.ladicabs;
  const mod2 = mod >> 1n;
  const Pg: number[] = [0];
  for (let j = 1; j < lP; j++) {
    const tau = permtopol(PG.gen[j]!, PL, PM, Pden, mod, mod2);
    const g = get_image(QPoly_to_FpX(tau, ip), Pp, Pmod, ip);
    if (!g) return null;
    Pg.push(g);
  }
  return { PG, Pg };
}

function lcmBig(a: bigint, b: bigint): bigint {
  let x = babs(a);
  let y = babs(b);
  if (x === 0n || y === 0n) return 0n;
  const g = (function gcd(u: bigint, v: bigint): bigint {
    while (v) {
      const t = u % v;
      u = v;
      v = t;
    }
    return u;
  })(x, y);
  return (x / g) * y;
}

/** PARI `galoisgenfixedfield0` (galconj.c:2195) */
function galoisgenfixedfield0(
  O: Perm[],
  L: bigint[],
  sigma: QPoly,
  T: ZX,
  bad: bigint | null,
  gf: GaloisFrobenius,
  gb: GaloisBorne
): { PG: GaloisGens; Pg: number[]; V: { sym: SymPol; PL: bigint[]; P: ZX } } | null {
  const mod = gb.ladicabs;
  const mod2 = mod >> 1n;
  const OL = fixedfieldorbits(O, L);
  const sym = fixedfieldsympol(OL, gb.l);
  const PL = sympol_eval(sym, OL, mod);
  const P = FpX_center(FpV_roots_to_pol(PL, mod), mod, mod2);
  if (!FpX_is_squarefree(FpX_red(P, BigInt(gf.p)), BigInt(gf.p))) {
    const badp = lcmBig(bad !== null ? bad : gb.dis, ZX_disc(P));
    const out = { Tmod: gf.Tmod, psi: gf.psi, p: gf.p };
    findpsi(badp, gf.p, T, sigma, gf.deg, out);
    gf.p = out.p;
    gf.Tmod = out.Tmod;
    gf.psi = out.psi;
  }
  const p = BigInt(gf.p);
  const Tp = FpX_red(T, p);
  const Sp = sympol_aut_evalmod(sym, gf.deg, sigma, Tp, p);
  const Pmod = fixedfieldfactmod(Sp, p, gf.Tmod);
  const PG = galoisgenfixedfield(Pmod, PL, P, p, bad, gb);
  if (PG === null) return null;
  return { PG: PG.PG, Pg: PG.Pg, V: { sym, PL, P } };
}

/** PARI `stpow` (galconj.c:2227) */
function stpow(s: number, e: number, m: number): number {
  let n = 1;
  for (let i = 1; i < e; i++) n = (1 + n * s) % m;
  return n;
}

/** PARI `wpow` (galconj.c:2235) */
function wpow(s: number, m: number, e: number, n: number): number[] {
  const w: number[] = new Array(n + 1).fill(0);
  let si = s;
  w[1] = 1;
  for (let i = 2; i <= n; i++) w[i] = w[i - 1]! * e;
  for (let i = n; i >= 1; i--) {
    si = Fl_powu(si, e, m);
    w[i] = Fl_mul(s - 1 < 0 ? s - 1 + m : s - 1, stpow(si, w[i]!, m), m);
  }
  return w;
}

/** PARI `galoisgenliftauto` (galconj.c:2251) */
function galoisgenliftauto(O: Perm[], gj: Perm, s: number, n: number, td: GaloisTest): Perm | null {
  const deg = O[1]!.length - 1;
  const X: number[] = new Array(O.length).fill(0);
  const oX: number[] = new Array(O.length).fill(0);
  const B = perm_cycles(gj);
  const oj = B[1]!.length - 1;
  const F = factoru_small(oj);
  let pf = identity_perm(n);
  for (let k = F.length - 1; k >= 0; k--) {
    const p = F[k]![0];
    const e = F[k]![1];
    const op = oj / Math.pow(p, e);
    let dg = 1;
    let el = oj;
    let osel = 1;
    let a = 0;
    let pf1: Perm | null = null;
    const Be: Perm[][] = new Array(e + 1).fill(null).map(() => []);
    Be[e] = cyc_pow(B, op);
    for (let i = e - 1; i >= 1; i--) Be[i] = cyc_pow(Be[i + 1]!, p);
    const w = wpow(Fl_powu(s, op, deg), deg, p, e);
    const wg: number[] = new Array(e + 2).fill(0);
    wg[e + 1] = deg;
    for (let i = e; i >= 1; i--) wg[i] = ugcd(wg[i + 1]!, w[i]!);
    for (let i = 1; i < O.length; i++) oX[i] = 0;
    for (let f = 1; f <= e; f++) {
      const Bel = Be[f]!;
      dg *= p;
      el /= p;
      const sel = Fl_powu(s, el, deg);
      const sr = ugcd(stpow(sel, p, deg), deg);
      pf1 = null;
      let t = 0;
      for (; t < sr; t++)
        if ((a + t * w[f]!) % wg[f + 1]! === 0) {
          for (let i = 1; i < X.length; i++) X[i] = 0;
          for (let i = 0; i < X.length - 1; i += dg)
            for (let j = 1, kk = p, st = t; kk <= dg; j++, kk += p) {
              X[kk + i] = (oX[j + i]! + st) % deg;
              st = (t + st * osel) % deg;
            }
          pf1 = testpermutation(O, Bel, X, sel, p, sr, td);
          if (pf1) break;
        }
      if (!pf1) return null;
      for (let i = 1; i < O.length; i++) oX[i] = X[i]!;
      osel = sel;
      a = (a + t * w[f]!) % deg;
    }
    pf = perm_mul(pf, perm_powu(pf1!, el));
  }
  return pf;
}

/** PARI `galoisgenlift` (galconj.c:2746) */
function galoisgenlift(
  PG: GaloisGens,
  Pg: number[],
  O: Perm[],
  L: bigint[],
  M: bigint[][],
  frob: Perm,
  gb: GaloisBorne,
  gf: GaloisFrobenius
): GaloisGens | null {
  const lP = PG.gen.length;
  const n = L.length - 1;
  const td = inittest(L, M, gb.bornesol, gb.ladicsol);
  const gen: Perm[] = [[], frob.slice()];
  const orders: number[] = [0, gf.deg, ...PG.orders.slice(1)];
  for (let j = 1; j < lP; j++) {
    const pf = galoisgenliftauto(O, PG.gen[j]!, gf.psi[Pg[j]!]!, n, td);
    if (!pf) return null;
    gen.push(pf);
  }
  return { gen, orders };
}

/** PARI `psi_order` (galconj.c:2770) */
function psi_order(psi: number[], d: number): number {
  let s = 1;
  for (let i = 1; i < psi.length; i++) s = ulcm(s, d / ugcd(psi[i]! - 1 < 0 ? psi[i]! - 1 + d : psi[i]! - 1, d));
  return s;
}

/** PARI `galoisgen` (galconj.c:2780) */
export function galoisgen(
  T: ZX,
  L: bigint[],
  M: bigint[][],
  den: bigint,
  bad: bigint | null,
  gb: GaloisBorne,
  ga: GaloisAnalysis
): GaloisGens | null {
  const n = ZX_degree(T);
  if (!ga.deg) return null;
  const special = (name: string) =>
    new NotImplementedError(
      `SAGE_NOT_IMPLEMENTED: galconj ${name} (the dedicated generator search ` +
        'upstream runs for this group is not ported; the generic weakly-super-solvable ' +
        'search does not apply to it)'
    );
  if (n === 12 && ga.ord === 3 && !ga.p4) {
    /* A4 is very probable: test it first (galconj.c:2794) */
    const td = inittest(L, M, gb.bornesol, gb.ladicsol);
    const PG = a4galoisgen(td);
    if (PG) return PG;
  }
  if (n === 24 && ga.ord === 3 && ga.p4) throw special('s4galoisgen');
  if (n === 36 && ga.ord === 3 && ga.p4) throw special('f36galoisgen');
  const gf: GaloisFrobenius = { p: 0, fp: 0, deg: 0, Tmod: [], psi: [] };
  const frob = galoisfindfrobenius(T, L, den, bad, gf, gb, ga);
  if (!frob) {
    if (n === 12 && ga.ord === 3 && !ga.p4) throw special('a4galoisgen');
    return null;
  }
  const po = psi_order(gf.psi, gf.deg);
  if (!(ga.group & ga_easy) && po < gf.deg && (gf.deg / radicalu(gf.deg)) % po === 0)
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: galconj galoisgenlift_nilp (the central/nilpotent ' +
        'shortcut, only reachable for degree > 104)'
    );
  const sigma = permtopol(frob, L, M, den, gb.ladicabs, gb.ladicabs >> 1n);
  if (gf.deg === n) return { gen: [[], frob], orders: [0, n] }; /* cyclic */
  const O = perm_cycles(frob);
  const PG = galoisgenfixedfield0(O, L, sigma, T, null, gf, gb);
  if (PG === null) return null;
  return galoisgenlift(PG.PG, PG.Pg, O, L, M, frob, gb, gf);
}

/* ================================================================== */
/*  Public entry points (galconj.c:2976-3470)                         */
/* ================================================================== */

/** The structure PARI's `galoisinit` returns (galconj.c:3044, an 8-component vector) */
export interface GaloisInit {
  /** gal[1]: the defining polynomial */
  pol: ZX;
  /** gal[2][1]: the prime `l` (totally split in the field) */
  p: bigint;
  /** gal[2][2]: the l-adic accuracy */
  e: number;
  /** gal[2][3]: `l^e` */
  mod: bigint;
  /** gal[3]: the l-adic roots, 1-indexed */
  roots: bigint[];
  /** gal[4]: `den` times the inverse Vandermonde matrix, `M[k][i]` */
  invvdm: bigint[][];
  /** gal[5]: a denominator */
  den: bigint;
  /** gal[6]: all the elements of the group, as permutations of the roots */
  group: Perm[];
  /** gal[7]: the generators */
  gen: Perm[];
  /** gal[8]: their relative orders */
  orders: number[];
}

function checkZXmonic(T: ZX, fun: string): void {
  const n = ZX_degree(T);
  if (n <= 0) throw new PariIrredpolError(fun);
  if (!ZX_is_monic(T)) throw new PariImplError(`${fun}(nonmonic)`);
  if (!ZX_is_squarefree(T))
    throw new PariDomainError(fun, 'issquarefree(pol)', '=', '0');
}

/**
 * PARI `galoisconj4_main` (galconj.c:2987) with `flag = 1`, i.e.
 * `galoisinit(T, den)` (galconj.c:3176).
 *
 * Returns `null` when the field is not Galois, or when its Galois group is not
 * weakly super solvable (PARI returns `gen_0` in both cases).
 */
export function galoisinit(T: ZX, den: bigint | null = null): GaloisInit | null {
  const pol = ZX_renormalize(T);
  const n = ZX_degree(pol);
  checkZXmonic(pol, 'galoisinit');
  const ga: GaloisAnalysis = { p: 0, deg: 0, ord: 0, l: 0, p4: 0, group: 0 };
  if (n === 1) {
    ga.l = 3;
    ga.deg = 1;
    den = 1n;
  } else if (!galoisanalysis(pol, ga, true, null)) return null;
  if (den !== null) den = babs(den);
  const gb: GaloisBorne = {
    l: BigInt(ga.l),
    valsol: 0,
    valabs: 0,
    bornesol: 0n,
    ladicsol: 0n,
    ladicabs: 0n,
    dis: 0n,
  };
  const d = galoisborne(pol, den, gb, n);
  const L = ZpX_roots(pol, gb.l, gb.valabs);
  if (L.length - 1 !== n) throw new PariBugError('galoisinit: l is not totally split');
  const M = FpV_invVandermonde(L, d, gb.ladicabs);
  let G: GaloisGens | null;
  if (n === 1) G = { gen: [[]], orders: [0] };
  else G = galoisgen(pol, L, M, d, null, gb, ga);
  if (!G) return null;
  const grp: GaloisInit = {
    pol,
    p: gb.l,
    e: gb.valabs,
    mod: gb.ladicabs,
    roots: L,
    invvdm: M,
    den: d,
    group: group_elts({ gen: G.gen, ord: G.orders }, n),
    gen: G.gen,
    orders: G.orders,
  };
  certify(grp);
  return grp;
}

/**
 * Sanity check that has no counterpart upstream: verify that the permutations
 * we return really do form a group of order `deg T` acting simply transitively
 * on the roots.  This guarantees we never return a plausible wrong answer
 * (upstream trusts its own p-adic bounds).
 */
function certify(grp: GaloisInit): void {
  const n = ZX_degree(grp.pol);
  const elts = grp.group;
  if (elts.length - 1 !== n) throw new PariBugError('galoisinit: wrong group order');
  const seen = new Set<string>();
  for (let i = 1; i < elts.length; i++) {
    const s = elts[i]!;
    if (s.length !== n + 1) throw new PariBugError('galoisinit: bad permutation');
    const key = s.slice(1).join(',');
    if (seen.has(key)) throw new PariBugError('galoisinit: repeated group element');
    seen.add(key);
    /* each element must permute the l-adic roots */
    const sig = permtopol(s, grp.roots, grp.invvdm, grp.den, grp.mod, grp.mod >> 1n);
    for (let i2 = 1; i2 <= n; i2++) {
      const v = bmod(
        FpX_eval(FpX_red(sig.num, grp.mod), grp.roots[i2]!, grp.mod) * Fp_inv(sig.den, grp.mod),
        grp.mod
      );
      if (v !== bmod(grp.roots[s[i2]!]!, grp.mod))
        throw new PariBugError('galoisinit: permutation is not an automorphism');
    }
  }
  /* closure */
  for (let i = 1; i < elts.length; i++)
    for (let j = 1; j < elts.length; j++)
      if (!seen.has(perm_mul(elts[i]!, elts[j]!).slice(1).join(',')))
        throw new PariBugError('galoisinit: the group is not closed');
}

/** PARI `galoisvecpermtopol` (galconj.c:1057) */
export function galoisvecpermtopol(gal: GaloisInit, vec: Perm[]): QPoly[] {
  const out: QPoly[] = [];
  for (const p of vec) {
    if (p.length !== gal.roots.length) throw new PariTypeError('galoispermtopol', 'perm');
    out.push(permtopol(p, gal.roots, gal.invvdm, gal.den, gal.mod, gal.mod >> 1n));
  }
  return out;
}

/** PARI `galoispermtopol` (galconj.c:3202) */
export function galoispermtopol(gal: GaloisInit, perm: Perm): QPoly {
  return permtopol(perm, gal.roots, gal.invvdm, gal.den, gal.mod, gal.mod >> 1n);
}

/** PARI `galoiscosets` (galconj.c:3213) */
function galoiscosets(O: Perm[], perm: Perm[]): number[] {
  const l = O.length;
  const C: number[] = new Array(l).fill(0);
  const o = O[1]!;
  const f = o.length;
  const u = o[1]!;
  const RC: boolean[] = new Array(perm[1]!.length).fill(false);
  for (let i = 1, j = 1; j < l; i++) {
    const p = perm[i]!;
    if (RC[p[u]!]) continue;
    for (let k = 1; k < f; k++) RC[p[o[k]!]!] = true;
    C[j++] = i;
  }
  return C;
}

/** PARI `fixedfieldfactor` (galconj.c:3230) */
function fixedfieldfactor(
  L: bigint[],
  O: Perm[],
  perm: Perm[],
  M: bigint[][],
  den: bigint,
  mod: bigint,
  mod2: bigint
): QPoly[][] {
  const l = O.length;
  const lo = O[1]!.length;
  const cosets = galoiscosets(O, perm);
  const res: QPoly[][] = [];
  for (let i = 1; i < l; i++) {
    const F: QPoly[] = new Array(lo + 1).fill(null);
    F[lo] = { num: [1n], den: 1n };
    const Lp = vecpermute(L, perm[cosets[i]!]!);
    const G: bigint[][] = [[]];
    for (let k = 1; k < l; k++) {
      const roots = vecpermute(Lp, O[k]!);
      G.push(FpV_roots_to_pol(roots, mod));
    }
    for (let j = 1; j < lo; j++) {
      const V: bigint[] = [0n];
      /* gmael(G,k,j+1) is the coefficient of x^(j-1) of the t_POL G[k] */
      for (let k = 1; k < l; k++) V.push(G[k]![j - 1] ?? 0n);
      F[j] = vectopol(V, M, den, mod, mod2);
    }
    res.push(F.slice(1, lo + 1));
  }
  return res;
}

/** the result of {@link galoisfixedfield} */
export interface FixedField {
  /** the polynomial defining the fixed field */
  P: ZX;
  /** a root of `P` expressed in `Q[x]/(T)` (`flag != 1`) */
  S?: QPoly;
  /** the factorisation of `T` over the fixed field (`flag = 2`) */
  factors?: QPoly[][];
}

/** PARI `galoisfixedfield` (galconj.c:3276) */
export function galoisfixedfield(
  gal: GaloisInit,
  perm: Perm | Perm[] | Group,
  flag: 0 | 1 | 2 = 0
): FixedField {
  if (flag < 0 || flag > 2) throw new PariFlagError('galoisfixedfield');
  const T = gal.pol;
  let L = gal.roots;
  const n = L.length - 1;
  let mod = gal.mod;
  let O: Perm[];
  const chk = (p: Perm) => {
    if (p.length !== n + 1) throw new PariTypeError('galoisfixedfield', 'perm');
  };
  if (Array.isArray(perm) && perm.length > 0 && Array.isArray(perm[0])) {
    const v = perm as Perm[];
    const gens = v[0]!.length === 0 ? v.slice(1) : v;
    for (const p of gens) chk(p);
    O = vecperm_orbits(gens, n);
  } else if (!Array.isArray(perm)) {
    const G = perm as Group;
    const gens = G.gen.slice(1);
    for (const p of gens) chk(p);
    O = vecperm_orbits(gens, n);
  } else {
    const p = perm as Perm;
    chk(p);
    O = perm_cycles(p);
  }
  let mod2 = mod >> 1n;
  const OL = fixedfieldorbits(O, L);
  const sym = fixedfieldsympol(OL, gal.p);
  let PL = sympol_eval(sym, OL, mod);
  const P = FpX_center(FpV_roots_to_pol(PL, mod), mod, mod2);
  if (flag === 1) return { P };
  const Sv = fixedfieldinclusion(O, PL);
  const S = vectopol(Sv, gal.invvdm, gal.den, mod, mod2);
  if (flag === 0) return { P, S };
  const Pgb: GaloisBorne = {
    l: gal.p,
    valsol: 0,
    valabs: 0,
    bornesol: 0n,
    ladicsol: 0n,
    ladicabs: 0n,
    dis: 0n,
  };
  const val = gal.e;
  const Pden = galoisborne(P, null, Pgb, ZX_degree(T) / ZX_degree(P));
  if (Pgb.valabs > val) {
    PL = ZpX_liftroots(P, PL, Pgb.l, Pgb.valabs);
    L = ZpX_liftroots(T, L, Pgb.l, Pgb.valabs);
    mod = Pgb.ladicabs;
    mod2 = mod >> 1n;
  }
  const PM = FpV_invVandermonde(PL, Pden, mod);
  const factors = fixedfieldfactor(L, O, gal.group, PM, Pden, mod, mod2);
  return { P, S, factors };
}

/** PARI `galois_group` (galconj.c:3337) */
export function galois_group(gal: GaloisInit): Group {
  return { gen: gal.gen, ord: gal.orders };
}

/** PARI `galoissubgroups` (galconj.c:3450) */
export function galoissubgroups(gal: GaloisInit | Group): Group[] {
  const G: Group = 'pol' in gal ? galois_group(gal) : gal;
  return group_subgroups(G);
}

/**
 * PARI `galoisconj4` (galconj.c:3104) / `nfgaloisconj(T, 4)`: the automorphisms
 * of `Q[x]/(T)` as polynomials, sorted.
 */
export function galoisconj4(T: ZX, den: bigint | null = null): QPoly[] | null {
  const gal = galoisinit(T, den);
  if (!gal) return null;
  const aut = galoisvecpermtopol(gal, gal.group.slice(1));
  aut.sort((a, b) => {
    const A = QPoly_to_fractions(a);
    const B = QPoly_to_fractions(b);
    if (A.length !== B.length) return A.length - B.length;
    for (let i = A.length - 1; i >= 0; i--) {
      const x = A[i]![0] * B[i]![1];
      const y = B[i]![0] * A[i]![1];
      if (x !== y) return x < y ? -1 : 1;
    }
    return 0;
  });
  return aut;
}
