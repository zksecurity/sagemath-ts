/**
 * @module parigp-ts/elliptic/ellsea
 * @description The Schoof-Elkies-Atkin point counting algorithm.
 *
 * Port of `reference/pari/src/basemath/ellsea.c` (2113 lines), restricted to
 * prime fields (`T = NULL` everywhere in upstream, i.e. `Fq_ellcard_SEA` with
 * `q = p`).  The structure, function names and control flow follow upstream
 * statement by statement; every function carries the upstream file:line it was
 * transcribed from.
 *
 * The modular equations come from our port of `polmodular.c`: upstream's
 * `get_modular_eqn` (ellsea.c:108-123) reads them from the optional `seadata`
 * package and, when that is absent, falls back to
 * `polmodular_ZXX(ell, ell==3? 0: 5, vx, vy)` -- exactly what we do (`seadata`
 * is not vendored, so the `'A'` and `'C'` branches of `find_isogenous` are
 * unreachable; see `find_isogenous`).
 *
 * Upstream reference points:
 * - `Fq_ellcard_SEA`          ellsea.c:1978-2109 (main loop, CRT of the trace)
 * - `find_trace`              ellsea.c:1505-1573 (per-prime dispatch)
 * - `find_trace_Elkies_power` ellsea.c:1396-1451 (Elkies, trace mod ell^k)
 * - `find_trace_Atkin`        ellsea.c:1455-1485 (Atkin, trace set)
 * - `match_and_sort`          ellsea.c:1815-1958 (BSGS over the Atkin data)
 * - `find_isogenous_from_J`   ellsea.c:1164-1222 (isogenous curve)
 * - `find_kernel`             ellsea.c:845-861   (kernel polynomial)
 * - `find_eigen_value_power`  ellsea.c:808-829   (eigenvalue of Frobenius)
 * - division polynomials      ellsea.c:225-520
 * - supersingularity test     FpE.c:698-805 (`Fp_elljissupersingular`)
 *
 * Not ported, because unreachable without `seadata`: the `'A'` (Atkin) and
 * `'C'` (canonical) modular-equation types of `find_isogenous`
 * (`find_isogenous_from_Atkin`, ellsea.c:900-960, and
 * `find_isogenous_from_canonical`, ellsea.c:964-1051).  `find_isogenous`
 * throws a `NotImplementedError` naming them.  Extension fields
 * (`T != NULL`) are likewise out of scope: the exported entry point is
 * `Fp_ellcard_SEA`.
 *
 * Timings (Apple silicon, bun 1.3.10), single process, cold modular-equation
 * cache: 65 bits 0.03 s, 100 bits 0.1 s, 128 bits 29 s, 159 bits 46 s,
 * 191 bits 30 s, 255 bits (Curve25519) 232 s, 256 bits (NIST P-256) 262 s.
 * Once the modular equations are cached a further 256-bit curve takes 13 s;
 * computing `polmodular_ZXX(ell, INV_G2)` is 80-90% of the cold time (which is
 * exactly what PARI's `seadata` package exists to avoid).
 */

import { Fp_ellcard_CM, Fp_ellj_get_CM } from './group.js';
import {
  INV_G2,
  INV_J,
  PariBugError,
  polclass0,
  polmodular_ZXX,
  polmodular_db_init,
} from '../polmodular.js';
import { PariPrimeError } from '../matkermod.js';
import { NotImplementedError } from '../ifactor.js';
import { Qfb, qfbsolve, Zp_sqrt } from '../qfb.js';
import { Fp_sqrt, kronecker } from '../ff.js';

/* ================================================================== */
/* Debugging (upstream's DEBUGLEVEL_ellsea traces, ellsea.c:37)        */
/* ================================================================== */

let DEBUGLEVEL = 0;
/** Mirror of PARI's `default(debug, n)` for the `ellsea` traces. */
export function setSeaDebugLevel(n: number): void {
  DEBUGLEVEL = n;
}
function err_printf(s: string): void {
  if (DEBUGLEVEL) process.stderr.write(s);
}

/* ================================================================== */
/* Z/N arithmetic ("Fq" with T = NULL; N is p or p^e)                 */
/* ================================================================== */

function mod(a: bigint, N: bigint): bigint {
  const r = a % N;
  return r < 0n ? r + N : r;
}

/** PARI `Fp_inv`; raises like `pari_err_INV` when `a` is not invertible. */
function Fp_inv(a: bigint, N: bigint): bigint {
  let [old_r, r] = [mod(a, N), N];
  let [old_s, s] = [1n, 0n];
  while (r !== 0n) {
    const q = old_r / r;
    [old_r, r] = [r, old_r - q * r];
    [old_s, s] = [s, old_s - q * s];
  }
  if (old_r !== 1n) throw new Error(`impossible inverse in Fp_inv: Mod(${mod(a, N)}, ${N})`);
  return mod(old_s, N);
}

const Fq_add = (a: bigint, b: bigint, N: bigint) => mod(a + b, N);
const Fq_sub = (a: bigint, b: bigint, N: bigint) => mod(a - b, N);
const Fq_neg = (a: bigint, N: bigint) => mod(-a, N);
const Fq_mul = (a: bigint, b: bigint, N: bigint) => mod(a * b, N);
const Fq_sqr = (a: bigint, N: bigint) => mod(a * a, N);
const Fq_mulu = (a: bigint, b: number | bigint, N: bigint) => mod(a * BigInt(b), N);
const Fq_div = (a: bigint, b: bigint, N: bigint) => mod(a * Fp_inv(b, N), N);
const Fq_halve = (a: bigint, N: bigint) => mod(a * Fp_inv(2n, N), N);

function Fq_powu(a: bigint, n: number, N: bigint): bigint {
  let r = 1n;
  let b = mod(a, N);
  let e = n;
  while (e > 0) {
    if (e & 1) r = mod(r * b, N);
    e >>= 1;
    if (e > 0) b = mod(b * b, N);
  }
  return r;
}

function Fp_pow(a: bigint, n: bigint, N: bigint): bigint {
  let r = 1n;
  let b = mod(a, N);
  let e = n;
  while (e > 0n) {
    if (e & 1n) r = mod(r * b, N);
    e >>= 1n;
    if (e > 0n) b = mod(b * b, N);
  }
  return r;
}

/** PARI `Zq_inv` (ellsea.c:559-563), `T = NULL`. */
function Zq_inv(b: bigint, pp: bigint, e: number): bigint {
  return e === 1 ? Fp_inv(b, pp) : Fp_inv(b, pp ** BigInt(e));
}
/** PARI `Zq_div` (ellsea.c:565-570). */
function Zq_div(a: bigint, b: bigint, N: bigint, pp: bigint, e: number): bigint {
  if (e === 1) return Fq_div(a, b, N);
  return Fq_mul(a, Zq_inv(b, pp, e), N);
}
/** PARI `Zq_sqrt` (ellsea.c:572-577). */
function Zq_sqrt(b: bigint, pp: bigint, e: number): bigint {
  if (e === 1) {
    const r = Fp_sqrt(mod(b, pp), pp);
    if (r === null) throw new Error('not a square in Zq_sqrt');
    return r;
  }
  const r = Zp_sqrt(mod(b, pp ** BigInt(e)), pp, e);
  if (r === null) throw new Error('not a square in Zq_sqrt');
  return r;
}
/** PARI `Zq_pval` (ellsea.c:583-585), `T = NULL`. */
function Zq_pval(a: bigint, p: bigint): number {
  if (a === 0n) return 0;
  let v = 0;
  let x = a;
  while (x % p === 0n) {
    x /= p;
    v++;
  }
  return v;
}
/** PARI `u_pvalrem`. */
function u_pvalrem(b: bigint, p: bigint): { v: number; n: bigint } {
  let v = 0;
  let n = b;
  while (n % p === 0n) {
    n /= p;
    v++;
  }
  return { v, n };
}
/**
 * PARI `Zq_divu_safe` (ellsea.c:587-601): `a / b` in `Z/q`, `q = p^e`,
 * returning `null` when the `p`-adic valuation of `b` exceeds that of `a`.
 */
function Zq_divu_safe(a: bigint, b: number, N: bigint, pp: bigint, e: number): bigint | null {
  if (e === 1) return Fq_div(a, BigInt(b), N);
  const { v, n } = u_pvalrem(BigInt(b), pp);
  let A = a;
  if (v > 0) {
    if (A === 0n) return 0n;
    const w = Zq_pval(A, pp);
    if (v > w) return null;
    A = A / pp ** BigInt(v);
  }
  return Fq_mul(A, Fp_inv(n, N), N);
}

/** PARI `ulogint(x, p)`: largest `k` with `p^k <= x`. */
function ulogint(x: number, p: bigint): number {
  let k = 0;
  let q = 1n;
  const X = BigInt(x);
  for (;;) {
    q *= p;
    if (q > X) return k;
    k++;
  }
}

/** PARI `quadratic_prec_mask` (Zp.c:37-47). */
function quadratic_prec_mask(n: number): bigint {
  let a = n;
  let mask = 0n;
  for (let i = 1; ; i++, mask <<= 1n) {
    mask |= BigInt(a & 1);
    a = (a + 1) >> 1;
    if (a === 1) return mask | (1n << BigInt(i));
  }
}

/* ================================================================== */
/* FpX: univariate polynomials over Z/N                               */
/* little-endian coefficient arrays, no trailing zero, zero poly = [] */
/* ================================================================== */

type FpX = bigint[];

function FpX_renormalize(a: FpX): FpX {
  let i = a.length;
  while (i > 0 && a[i - 1] === 0n) i--;
  a.length = i;
  return a;
}
function FpX_red(a: readonly bigint[], N: bigint): FpX {
  return FpX_renormalize(a.map((c) => mod(c, N)));
}
function FpX_deg(a: FpX): number {
  return a.length - 1;
}
function pol_x(): FpX {
  return [0n, 1n];
}
function FpX_add(a: FpX, b: FpX, N: bigint): FpX {
  const n = Math.max(a.length, b.length);
  const r = new Array<bigint>(n);
  for (let i = 0; i < n; i++) r[i] = mod((a[i] ?? 0n) + (b[i] ?? 0n), N);
  return FpX_renormalize(r);
}
function FpX_sub(a: FpX, b: FpX, N: bigint): FpX {
  const n = Math.max(a.length, b.length);
  const r = new Array<bigint>(n);
  for (let i = 0; i < n; i++) r[i] = mod((a[i] ?? 0n) - (b[i] ?? 0n), N);
  return FpX_renormalize(r);
}
function FpX_neg(a: FpX, N: bigint): FpX {
  return a.map((c) => (c === 0n ? 0n : N - c));
}
function FpX_Fp_mul(a: FpX, c: bigint, N: bigint): FpX {
  const cc = mod(c, N);
  if (cc === 0n) return [];
  return FpX_renormalize(a.map((x) => mod(x * cc, N)));
}
function FpX_Fp_add(a: FpX, c: bigint, N: bigint): FpX {
  const r = a.slice();
  if (r.length === 0) return FpX_renormalize([mod(c, N)]);
  r[0] = mod(r[0]! + c, N);
  return FpX_renormalize(r);
}
function FpX_mulu(a: FpX, c: number, N: bigint): FpX {
  return FpX_Fp_mul(a, BigInt(c), N);
}
function FpX_mul(a: FpX, b: FpX, N: bigint): FpX {
  if (a.length === 0 || b.length === 0) return [];
  const r = new Array<bigint>(a.length + b.length - 1).fill(0n);
  for (let i = 0; i < a.length; i++) {
    const ai = a[i]!;
    if (ai === 0n) continue;
    for (let j = 0; j < b.length; j++) {
      const bj = b[j]!;
      if (bj === 0n) continue;
      r[i + j] = (r[i + j]! + ai * bj) % N;
    }
  }
  return FpX_renormalize(r);
}
function FpX_sqr(a: FpX, N: bigint): FpX {
  return FpX_mul(a, a, N);
}
function FpX_halve(a: FpX, N: bigint): FpX {
  return FpX_Fp_mul(a, Fp_inv(2n, N), N);
}
/** Euclidean division; the leading coefficient of `b` must be invertible. */
function FpX_divrem(a: FpX, b: FpX, N: bigint): [FpX, FpX] {
  if (b.length === 0) throw new Error('impossible inverse in FpX_divrem: 0');
  const db = b.length - 1;
  if (a.length - 1 < db) return [[], a.slice()];
  const inv = Fp_inv(b[db]!, N);
  const r = a.slice();
  const q = new Array<bigint>(a.length - db).fill(0n);
  for (let i = r.length - 1; i >= db; i--) {
    const c = r[i]!;
    if (c === 0n) continue;
    const qc = mod(c * inv, N);
    q[i - db] = qc;
    const off = i - db;
    for (let j = 0; j <= db; j++) r[off + j] = mod(r[off + j]! - qc * b[j]!, N);
  }
  return [FpX_renormalize(q), FpX_renormalize(r)];
}
function FpX_rem(a: FpX, b: FpX, N: bigint): FpX {
  return FpX_divrem(a, b, N)[1]!;
}
function FpX_div(a: FpX, b: FpX, N: bigint): FpX {
  return FpX_divrem(a, b, N)[0]!;
}
/** PARI `FpX_normalize`: make monic. */
function FpX_normalize(a: FpX, N: bigint): FpX {
  if (a.length === 0) return a;
  const lc = a[a.length - 1]!;
  if (lc === 1n) return a.slice();
  return FpX_Fp_mul(a, Fp_inv(lc, N), N);
}
function FpX_gcd(a: FpX, b: FpX, N: bigint): FpX {
  let x = a.slice();
  let y = b.slice();
  while (y.length !== 0) {
    const r = FpX_rem(x, y, N);
    x = y;
    y = r;
  }
  return x;
}
function FpX_deriv(a: FpX, N: bigint): FpX {
  if (a.length <= 1) return [];
  const r = new Array<bigint>(a.length - 1);
  for (let i = 1; i < a.length; i++) r[i - 1] = mod(a[i]! * BigInt(i), N);
  return FpX_renormalize(r);
}
function FpX_eval(a: FpX, x: bigint, N: bigint): bigint {
  let s = 0n;
  for (let i = a.length - 1; i >= 0; i--) s = mod(s * x + a[i]!, N);
  return s;
}
/** PARI `RgX_shift_shallow`: multiply by `x^n` (drop low terms if `n < 0`). */
function FpX_shift(a: FpX, n: number): FpX {
  if (n === 0) return a.slice();
  if (n > 0) return a.length === 0 ? [] : new Array<bigint>(n).fill(0n).concat(a);
  const m = -n;
  return m >= a.length ? [] : a.slice(m);
}
/** PARI `RgXn_red_shallow`: truncate mod `x^n`. */
function FpXn_red(a: FpX, n: number): FpX {
  if (n <= 0) return [];
  return FpX_renormalize(a.slice(0, n));
}
function FpXn_mul(a: FpX, b: FpX, n: number, N: bigint): FpX {
  return FpXn_red(FpX_mul(a, b, N), n);
}
/** PARI `FpX_mulhigh_i` (FpX.c). */
function FpX_mulhigh_i(f: FpX, g: FpX, n: number, N: bigint): FpX {
  return FpX_shift(FpX_mul(f, g, N), -n);
}
/** PARI `RgX_blocks(f, n, 2)`. */
function FpX_blocks2(f: FpX, n: number): [FpX, FpX] {
  return [FpXn_red(f, n), FpX_shift(f, -n)];
}
/** PARI `FpXn_mulhigh` (FpX.c). */
function FpXn_mulhigh(f: FpX, g: FpX, n2: number, n: number, N: bigint): FpX {
  const [fl, fh] = FpX_blocks2(f, n2);
  return FpX_add(FpX_mulhigh_i(fl, g, n2, N), FpXn_mul(fh, g, n - n2, N), N);
}
/**
 * PARI `FpX_integXn` (FpX.c): `sum x_i/(n+i) X^i`, using the gcd trick so that
 * a coefficient divisible by `p` may still be divided by a multiple of `p`.
 */
function FpX_integXn(x: FpX, n: number, N: bigint): FpX {
  if (x.length === 0) return [];
  const y = new Array<bigint>(x.length);
  for (let i = 0; i < x.length; i++) {
    const xi = x[i]!;
    if (xi === 0n) {
      y[i] = 0n;
      continue;
    }
    const j = BigInt(n + i + 1); /* PARI: j = n + i - 1 with its 1-based index */
    let d = j;
    let a = xi;
    while (d !== 0n) {
      const t = a % d;
      a = d;
      d = t;
    }
    const g = a; /* gcd(xi, j) */
    y[i] = Fq_div(xi / g, j / g, N);
  }
  return FpX_renormalize(y);
}
/**
 * Inverse of the power series `f` modulo `x^e`.  PARI computes this with a
 * Newton iteration (`FpXn_div`, FpX.c); the object is the unique inverse power
 * series, so we use the (quadratically slower but identical) direct
 * recurrence -- the degrees involved here are `O(ell)`.
 */
function FpXn_inv(f: FpX, e: number, N: bigint): FpX {
  if (f.length === 0) throw new Error('impossible inverse in FpXn_inv: 0');
  const i0 = Fp_inv(f[0]!, N);
  const r = new Array<bigint>(e).fill(0n);
  r[0] = i0;
  for (let k = 1; k < e; k++) {
    let s = 0n;
    for (let i = 1; i <= k; i++) {
      const fi = f[i] ?? 0n;
      if (fi !== 0n) s += fi * r[k - i]!;
    }
    r[k] = mod(-s % N * i0, N);
  }
  return FpX_renormalize(r);
}
/** PARI `FpXn_expint` (FpX.c): `exp(int(h))` mod `x^e`. */
function FpXn_expint(h: FpX, e: number, N: bigint): FpX {
  let n = 1;
  let f: FpX = [1n];
  let g: FpX = [1n];
  let mask = quadratic_prec_mask(e);
  while (mask > 1n) {
    const n2 = n;
    n <<= 1;
    if (mask & 1n) n--;
    mask >>= 1n;
    let u = FpXn_mul(g, FpX_mulhigh_i(f, FpXn_red(h, n2 - 1), n2 - 1, N), n - n2, N);
    u = FpX_add(u, FpX_shift(FpXn_red(h, n - 1), 1 - n2), N);
    const w = FpXn_mul(f, FpX_integXn(u, n2 - 1, N), n - n2, N);
    f = FpX_add(f, FpX_shift(w, n2), N);
    if (mask <= 1n) break;
    u = FpXn_mul(g, FpXn_mulhigh(f, g, n2, n, N), n - n2, N);
    g = FpX_sub(g, FpX_shift(u, n2), N);
  }
  return f;
}
/** PARI `RgX_recip` (RgX.c): `x^deg(P) P(1/x)`. */
function FpX_recip(x: FpX): FpX {
  return FpX_renormalize(x.slice().reverse());
}
/** PARI `deg1pol_shallow(a, b)` = `a*x + b`. */
function deg1pol(a: bigint, b: bigint, N: bigint): FpX {
  return FpX_renormalize([mod(b, N), mod(a, N)]);
}
/** PARI `FpX_div_by_X_x` (FpX.c). */
function FpX_div_by_X_x(a: FpX, x: bigint, N: bigint): { q: FpX; r: bigint } {
  const l = a.length;
  if (l <= 1) return { q: [], r: l === 0 ? 0n : a[0]! };
  const z = new Array<bigint>(l - 1);
  z[l - 2] = a[l - 1]!;
  for (let i = l - 3; i >= 0; i--) z[i] = mod(a[i + 1]! + x * z[i + 1]!, N);
  return { q: FpX_renormalize(z), r: mod(a[0]! + x * z[0]!, N) };
}

/* --- modular arithmetic in Fp[x]/(T) ------------------------------ */

function FpXQ_mul(a: FpX, b: FpX, T: FpX, N: bigint): FpX {
  return FpX_rem(FpX_mul(a, b, N), T, N);
}
function FpXQ_sqr(a: FpX, T: FpX, N: bigint): FpX {
  return FpX_rem(FpX_sqr(a, N), T, N);
}
function FpXQ_pow(a: FpX, n: bigint, T: FpX, N: bigint): FpX {
  let r: FpX = FpX_rem([1n], T, N);
  let b = FpX_rem(a, T, N);
  let e = n;
  while (e > 0n) {
    if (e & 1n) r = FpXQ_mul(r, b, T, N);
    e >>= 1n;
    if (e > 0n) b = FpXQ_sqr(b, T, N);
  }
  return r;
}
/** Modular composition `P(A) mod T` (Brent-Kung, as PARI's `FpX_FpXQV_eval`). */
function FpXQ_compose(P: FpX, A: FpX, T: FpX, N: bigint): FpX {
  const d = P.length;
  if (d === 0) return [];
  const k = Math.max(1, Math.ceil(Math.sqrt(d)));
  /* baby steps A^0 .. A^k */
  const pow: FpX[] = [[1n]];
  for (let i = 1; i <= k; i++) pow.push(FpXQ_mul(pow[i - 1]!, A, T, N));
  const nb = Math.ceil(d / k);
  let res: FpX = [];
  for (let b = nb - 1; b >= 0; b--) {
    let block: FpX = [];
    for (let j = 0; j < k; j++) {
      const c = P[b * k + j];
      if (c === undefined || c === 0n) continue;
      block = FpX_add(block, FpX_Fp_mul(pow[j]!, c, N), N);
    }
    res = FpX_add(FpXQ_mul(res, pow[k]!, T, N), block, N);
  }
  return res;
}
/** PARI `FpX_Frobenius`: `x^p mod T`. */
function FpX_Frobenius(T: FpX, p: bigint): FpX {
  return FpXQ_pow(pol_x(), p, T, p);
}

/* --- root finding (FpX_factor.c) ---------------------------------- */

let rngState = 0x9e3779b97f4a7c15n;
function nextRand(): bigint {
  rngState ^= (rngState << 13n) & 0xffffffffffffffffn;
  rngState ^= rngState >> 7n;
  rngState ^= (rngState << 17n) & 0xffffffffffffffffn;
  return rngState;
}
function randomFp(p: bigint): bigint {
  let r = 0n;
  for (let i = 0; i < 5; i++) r = (r << 64n) | nextRand();
  return r % p;
}

/** Split a polynomial that is a product of distinct linear factors. */
function FpX_split_linear(g: FpX, p: bigint, out: bigint[]): void {
  const d = FpX_deg(g);
  if (d <= 0) return;
  if (d === 1) {
    out.push(mod(-g[0]! * Fp_inv(g[1]!, p), p));
    return;
  }
  for (;;) {
    const a = randomFp(p);
    const t = FpXQ_pow(FpX_renormalize([a, 1n]), (p - 1n) / 2n, g, p);
    const h = FpX_gcd(FpX_sub(t, [1n], p), g, p);
    const dh = FpX_deg(h);
    if (dh <= 0 || dh === d) continue;
    FpX_split_linear(FpX_normalize(h, p), p, out);
    FpX_split_linear(FpX_normalize(FpX_div(g, h, p), p), p, out);
    return;
  }
}
/** PARI `FpX_roots`: the roots in `Fp` of `f`, without multiplicity. */
function FpX_roots(f: FpX, p: bigint): bigint[] {
  const T = FpX_normalize(FpX_red(f, p), p);
  if (FpX_deg(T) <= 0) return [];
  const xp = FpX_Frobenius(T, p);
  const g = FpX_gcd(FpX_sub(xp, pol_x(), p), T, p);
  if (FpX_deg(g) <= 0) return [];
  const out: bigint[] = [];
  FpX_split_linear(FpX_normalize(g, p), p, out);
  out.sort((a, b) => (a < b ? -1 : a > b ? 1 : 0));
  return out;
}
/** PARI `FpX_oneroot`: one root in `Fp`, or `null`. */
function FpX_oneroot(f: FpX, p: bigint): bigint | null {
  const r = FpX_roots(f, p);
  return r.length ? r[0]! : null;
}
/** PARI `FpX_nbroots`. */
function FpX_nbroots(f: FpX, p: bigint): number {
  const T = FpX_normalize(FpX_red(f, p), p);
  if (FpX_deg(T) <= 0) return 0;
  const xp = FpX_Frobenius(T, p);
  const g = FpX_gcd(FpX_sub(xp, pol_x(), p), T, p);
  return Math.max(FpX_deg(g), 0);
}
/** PARI `FpX_is_squarefree`. */
function FpX_is_squarefree(f: FpX, p: bigint): boolean {
  const T = FpX_red(f, p);
  if (T.length <= 1) return true;
  const d = FpX_deriv(T, p);
  if (d.length === 0) return false;
  return FpX_deg(FpX_gcd(T, d, p)) === 0;
}
/**
 * PARI `FpX_ddf_degree` (FpX_factor.c): the smallest `r >= 1` with
 * `x^(p^r) = x mod T`, i.e. the common degree of the irreducible factors of a
 * squarefree `T` all of whose factors have the same degree.  Transcribed with
 * upstream's baby-step/giant-step structure.
 */
function FpX_ddf_degree(T: FpX, XP: FpX, p: bigint): number {
  const n = FpX_deg(T);
  const X = pol_x();
  if (FpX_deg(XP) === 1 && XP[0] === 0n && XP[1] === 1n) return 1;
  const B = n >> 1;
  const l = Math.max(1, Math.floor(Math.sqrt(B)));
  const m = Math.floor((B + l - 1) / l);
  const h = new Map<string, number>();
  const key = (P: FpX) => P.join(',');
  h.set(key(X), 0);
  h.set(key(XP), 1);
  let b = XP;
  for (let i = 3; i <= l + 1; i++) {
    b = FpXQ_compose(b, XP, T, p);
    if (b.length === 2 && b[0] === 0n && b[1] === 1n) return i - 1;
    h.set(key(b), i - 1);
  }
  let g = b;
  for (let i = 2; i <= m + 1; i++) {
    g = FpXQ_compose(g, b, T, p);
    const j = h.get(key(g));
    if (j !== undefined) return l * i - j;
  }
  return n;
}

/* ================================================================== */
/* FpXY: bivariate polynomials over Z/N                               */
/* `Q[a]` is the coefficient (a polynomial in the inner variable vy)   */
/* of `vx^a`; so `Q[a][b]` is the coefficient of `vx^a vy^b`. This is  */
/* PARI's t_POL in vx with FpX (in vy) coefficients.                  */
/* ================================================================== */

type FpXY = FpX[];

function FpXY_renormalize(Q: FpXY): FpXY {
  let i = Q.length;
  while (i > 0 && Q[i - 1]!.length === 0) i--;
  Q.length = i;
  return Q;
}
/** PARI `FpXY_evaly(Q, y, p, vx)` (FpXX.c:1472): substitute `y` for `vx`. */
function FpXY_evaly(Q: FpXY, y: bigint, N: bigint): FpX {
  if (Q.length === 0) return [];
  let z = Q[Q.length - 1]!.slice();
  for (let i = Q.length - 2; i >= 0; i--) z = FpX_add(Q[i]!, FpX_Fp_mul(z, y, N), N);
  return z;
}
/** Derivative with respect to the main (outer) variable `vx`. */
function FpXY_deriv_x(Q: FpXY, N: bigint): FpXY {
  if (Q.length <= 1) return [];
  const r: FpXY = [];
  for (let i = 1; i < Q.length; i++) r.push(FpX_mulu(Q[i]!, i, N));
  return FpXY_renormalize(r);
}
/** PARI `RgX_splitting(P, 3)` on the outer variable: `H_k[a] = P[3a+k]`. */
function FpXY_splitting3(P: FpXY): [FpXY, FpXY, FpXY] {
  const H: [FpXY, FpXY, FpXY] = [[], [], []];
  for (let i = 0; i < P.length; i++) H[(i % 3) as 0 | 1 | 2].push(P[i]!);
  return [FpXY_renormalize(H[0]), FpXY_renormalize(H[1]), FpXY_renormalize(H[2])];
}
/**
 * PARI `RgXY_deflatex(H, n, d)` (ellsea.c:1070-1082): shift each coefficient
 * (a polynomial in the inner variable) by `d` and deflate it by `n`.
 */
function FpXY_deflatex(H: FpXY, n: number, d: number): FpXY {
  const R: FpXY = [];
  for (const Hi of H) {
    const s = FpX_shift(Hi, d);
    const q: FpX = [];
    for (let k = 0; k * n < s.length; k++) q.push(s[k * n]!);
    R.push(FpX_renormalize(q));
  }
  return FpXY_renormalize(R);
}

/* ================================================================== */
/* Modular equations (ellsea.c:60-124)                                */
/* ================================================================== */

/**
 * PARI's `struct meqn` (ellsea.c:84-89).  Without the `seadata` package only
 * `type = 'J'` occurs (ellsea.c:118-122).
 */
interface Meqn {
  type: 'J' | 'A' | 'C';
  eq: FpXY;
  evalR: FpX | null;
  evaldR: FpX | null;
  evalddR: FpX | null;
}

/**
 * Cache of the integral modular equations, mirroring upstream's `modular_eqn`
 * global (ellsea.c:38): `Phi_ell` only depends on `ell`.
 */
const modular_eqn_cache = new Map<number, FpXY>();

/**
 * PARI `get_modular_eqn` (ellsea.c:107-123).  `seadata` is not vendored, so we
 * always take upstream's fallback branch:
 * `type = 'J'` and `eq = polmodular_ZXX(ell, ell==3? 0: 5)`.
 */
function get_modular_eqn(ell: number): Meqn {
  let eq = modular_eqn_cache.get(ell);
  if (!eq) {
    const M = ell === 3 ? polmodular_ZXX(3, INV_J) : polmodular_ZXX(ell, INV_G2);
    eq = FpXY_renormalize(M.map((col) => FpX_renormalize(col.slice())));
    modular_eqn_cache.set(ell, eq);
  }
  return { type: 'J', eq, evalR: null, evaldR: null, evalddR: null };
}

/**
 * PARI `Fq_polmodular_eval` (ellsea.c:1084-1151): from `Phi^g2_N(X, Y)` (or
 * `Phi_3(X, Y)` for `N = 3`) build `[R, dR, ddR]` where `R(Y) = Phi_N(Y, j)`
 * is the classical modular polynomial evaluated at the curve's `j`, and `dR`,
 * `ddR` are its first two derivatives with respect to `j`.
 */
function Fq_polmodular_eval(meqn: FpXY, j: bigint, N: number, p: bigint): [FpX, FpX, FpX] {
  const t0 = N % 3 === 1 ? 2 : 0;
  const t2 = N % 3 === 1 ? 0 : 2;
  if (N === 3) {
    const P = meqn.map((c) => FpX_red(c, p));
    const dP = FpXY_deriv_x(P, p);
    const ddP = FpXY_deriv_x(dP, p);
    return [FpXY_evaly(P, j, p), FpXY_evaly(dP, j, p), FpXY_evaly(ddP, j, p)];
  }
  const P5 = meqn.map((c) => FpX_red(c, p));
  const H = FpXY_splitting3(P5);
  const H0 = FpXY_deflatex(H[0], 3, -t0);
  const H1 = FpXY_deflatex(H[1], 3, -1);
  const H2 = FpXY_deflatex(H[2], 3, -t2);
  const h0 = FpXY_evaly(H0, j, p);
  const h1 = FpXY_evaly(H1, j, p);
  const h2 = FpXY_evaly(H2, j, p);
  const dH0 = FpXY_deriv_x(H0, p);
  const dH1 = FpXY_deriv_x(H1, p);
  const dH2 = FpXY_deriv_x(H2, p);
  const ddH0 = FpXY_deriv_x(dH0, p);
  const ddH1 = FpXY_deriv_x(dH1, p);
  const ddH2 = FpXY_deriv_x(dH2, p);
  const d0 = FpXY_evaly(dH0, j, p);
  const d1 = FpXY_evaly(dH1, j, p);
  const d2 = FpXY_evaly(dH2, j, p);
  const dd0 = FpXY_evaly(ddH0, j, p);
  const dd1 = FpXY_evaly(ddH1, j, p);
  const dd2 = FpXY_evaly(ddH2, j, p);

  const h02 = FpX_sqr(h0, p);
  const h12 = FpX_sqr(h1, p);
  const h22 = FpX_sqr(h2, p);
  const h03 = FpX_mul(h0, h02, p);
  const h13 = FpX_mul(h1, h12, p);
  const h23 = FpX_mul(h2, h22, p);
  const h012 = FpX_mul(FpX_mul(h0, h1, p), h2, p);
  const dh03 = FpX_mul(FpX_mulu(d0, 3, p), h02, p);
  const dh13 = FpX_mul(FpX_mulu(d1, 3, p), h12, p);
  const dh23 = FpX_mul(FpX_mulu(d2, 3, p), h22, p);
  const dh012 = FpX_add(
    FpX_add(FpX_mul(FpX_mul(d0, h1, p), h2, p), FpX_mul(FpX_mul(h0, d1, p), h2, p), p),
    FpX_mul(FpX_mul(h0, h1, p), d2, p),
    p
  );
  const R1 = FpX_sub(h13, FpX_mulu(h012, 3, p), p);
  const j2 = Fq_sqr(j, p);
  const R = FpX_add(
    FpX_add(
      FpX_Fp_mul(FpX_shift(h23, t2), j2, p),
      FpX_Fp_mul(FpX_shift(R1, 1), j, p),
      p
    ),
    FpX_shift(h03, t0),
    p
  );
  const dR1 = FpX_sub(dh13, FpX_mulu(dh012, 3, p), p);
  const dR = FpX_add(
    FpX_add(
      FpX_shift(
        FpX_add(FpX_Fp_mul(dh23, j2, p), FpX_Fp_mul(h23, Fq_mulu(j, 2, p), p), p),
        t2
      ),
      FpX_shift(FpX_add(FpX_Fp_mul(dR1, j, p), R1, p), 1),
      p
    ),
    FpX_shift(dh03, t0),
    p
  );
  const ddh03 = FpX_mulu(
    FpX_add(FpX_mul(dd0, h02, p), FpX_mul(FpX_mulu(FpX_sqr(d0, p), 2, p), h0, p), p),
    3,
    p
  );
  const ddh13 = FpX_mulu(
    FpX_add(FpX_mul(dd1, h12, p), FpX_mul(FpX_mulu(FpX_sqr(d1, p), 2, p), h1, p), p),
    3,
    p
  );
  const ddh23 = FpX_mulu(
    FpX_add(FpX_mul(dd2, h22, p), FpX_mul(FpX_mulu(FpX_sqr(d2, p), 2, p), h2, p), p),
    3,
    p
  );
  const ddh012 = FpX_add(
    FpX_add(
      FpX_add(FpX_mul(FpX_mul(dd0, h1, p), h2, p), FpX_mul(FpX_mul(h0, dd1, p), h2, p), p),
      FpX_mul(FpX_mul(h0, h1, p), dd2, p),
      p
    ),
    FpX_mulu(
      FpX_add(
        FpX_add(FpX_mul(FpX_mul(d0, d1, p), h2, p), FpX_mul(FpX_mul(d0, h1, p), d2, p), p),
        FpX_mul(FpX_mul(h0, d1, p), d2, p),
        p
      ),
      2,
      p
    ),
    p
  );
  const ddR1 = FpX_sub(ddh13, FpX_mulu(ddh012, 3, p), p);
  const ddR2 = FpX_add(
    FpX_add(FpX_Fp_mul(ddh23, j2, p), FpX_Fp_mul(dh23, Fq_mulu(j, 4, p), p), p),
    FpX_mulu(h23, 2, p),
    p
  );
  const ddR = FpX_add(
    FpX_add(
      FpX_shift(ddR2, t2),
      FpX_shift(FpX_add(FpX_mulu(dR1, 2, p), FpX_Fp_mul(ddR1, j, p), p), 1),
      p
    ),
    FpX_shift(ddh03, t0),
    p
  );
  return [R, dR, ddR];
}

/** PARI `meqn_j` (ellsea.c:1153-1163). */
function meqn_j(M: Meqn, j: bigint, ell: number, N: bigint): FpX {
  const [R, dR, ddR] = Fq_polmodular_eval(M.eq, j, ell, N);
  M.evalR = R;
  M.evaldR = dR;
  M.evalddR = ddR;
  return R;
}

/* ================================================================== */
/* n-division polynomials (ellsea.c:225-520)                          */
/* ================================================================== */

/** PARI's `struct bb_algebra`, specialised to `Fp[x]` and `Fp[x]/(h)`. */
interface BBAlg {
  mul(a: FpX, b: FpX): FpX;
  sqr(a: FpX): FpX;
  add(a: FpX, b: FpX): FpX;
  sub(a: FpX, b: FpX): FpX;
  red(a: FpX): FpX;
  one(): FpX;
  zero(): FpX;
}
function get_FpX_algebra(N: bigint): BBAlg {
  return {
    mul: (a, b) => FpX_mul(a, b, N),
    sqr: (a) => FpX_sqr(a, N),
    add: (a, b) => FpX_add(a, b, N),
    sub: (a, b) => FpX_sub(a, b, N),
    red: (a) => FpX_red(a, N),
    one: () => [1n],
    zero: () => [],
  };
}
function get_FpXQ_algebra(h: FpX, N: bigint): BBAlg {
  return {
    mul: (a, b) => FpXQ_mul(a, b, h, N),
    sqr: (a) => FpXQ_sqr(a, h, N),
    add: (a, b) => FpX_add(a, b, N),
    sub: (a, b) => FpX_sub(a, b, N),
    red: (a) => FpX_rem(a, h, N),
    one: () => FpX_rem([1n], h, N),
    zero: () => [],
  };
}

/** PARI's `struct divpolmod_red` (ellsea.c:400-405). */
interface DivPolMod {
  ff: BBAlg;
  /** `t[0][n] = f_n`, `t[1][n] = f_n^2`, `t[2][n] = f_n f_{n-2}` */
  t: [Array<FpX | null>, Array<FpX | null>, Array<FpX | null>];
  r2: FpX;
}

/** PARI `divpol` (ellsea.c:269-306). */
function divpol(d: DivPolMod, n: number): FpX {
  const { t, ff, r2 } = d;
  const m = Math.floor(n / 2);
  if (n === 0) return ff.zero();
  const c = t[0][n];
  if (c) return c;
  let f: FpX;
  switch (n) {
    case 1:
    case 2:
      f = ff.one();
      break;
    default:
      if (n % 2 === 1) {
        if (m % 2 === 1)
          f = ff.sub(
            ff.mul(divpol_ff(d, m + 2), divpol_f2(d, m)),
            ff.mul(r2, ff.mul(divpol_ff(d, m + 1), divpol_f2(d, m + 1)))
          );
        else
          f = ff.sub(
            ff.mul(r2, ff.mul(divpol_ff(d, m + 2), divpol_f2(d, m))),
            ff.mul(divpol_ff(d, m + 1), divpol_f2(d, m + 1))
          );
      } else
        f = ff.sub(
          ff.mul(divpol_ff(d, m + 2), divpol_f2(d, m - 1)),
          ff.mul(divpol_ff(d, m), divpol_f2(d, m + 1))
        );
  }
  f = ff.red(f);
  t[0][n] = f;
  return f;
}
/** PARI `divpol_f2` (ellsea.c:248-256): `f_n^2`. */
function divpol_f2(d: DivPolMod, n: number): FpX {
  if (n === 0) return d.ff.zero();
  if (n <= 2) return d.ff.one();
  const c = d.t[1][n];
  if (c) return c;
  const r = d.ff.sqr(divpol(d, n));
  d.t[1][n] = r;
  return r;
}
/** PARI `divpol_ff` (ellsea.c:258-267): `f_n f_{n-2}`. */
function divpol_ff(d: DivPolMod, n: number): FpX {
  if (n <= 2) return d.ff.zero();
  const c = d.t[2][n];
  if (c) return c;
  if (n <= 4) return divpol(d, n);
  const r = d.ff.mul(divpol(d, n), divpol(d, n - 2));
  d.t[2][n] = r;
  return r;
}

/** PARI `Fq_elldivpol34` (ellsea.c:352-379). */
function Fq_elldivpol34(n: number, a4: bigint, a6: bigint, h: FpX | null, N: bigint): FpX {
  let res: FpX;
  switch (n) {
    case 3:
      res = FpX_red([Fq_neg(Fq_sqr(a4, N), N), Fq_mulu(a6, 12, N), Fq_mulu(a4, 6, N), 0n, 3n], N);
      break;
    case 4: {
      const a42 = Fq_sqr(a4, N);
      res = FpX_red(
        [
          Fq_sub(Fq_mul(Fq_sqr(a6, N), mod(-8n, N), N), Fq_mul(a4, a42, N), N),
          Fq_mul(Fq_mul(a4, a6, N), mod(-4n, N), N),
          Fq_mul(a42, mod(-5n, N), N),
          Fq_mulu(a6, 20, N),
          Fq_mulu(a4, 5, N),
          0n,
          1n,
        ],
        N
      );
      res = FpX_mulu(res, 2, N);
      break;
    }
    default:
      throw new PariBugError('Fq_elldivpol34');
  }
  if (h) res = FpX_rem(res, h, N);
  return res;
}
/** PARI `rhs` (ellsea.c:383-387): `x^3 + a4 x + a6`. */
function rhs(a4: bigint, a6: bigint, N: bigint): FpX {
  return FpX_red([a6, a4, 0n, 1n], N);
}
/** PARI `Fq_elldivpolmod_init` (ellsea.c:418-433) + `divpolmod_init`. */
function Fq_elldivpolmod_init(
  a4: bigint,
  a6: bigint,
  n: number,
  h: FpX | null,
  N: bigint
): DivPolMod {
  const D3 = n >= 0 ? Fq_elldivpol34(3, a4, a6, h, N) : null;
  const D4 = n >= 1 ? Fq_elldivpol34(4, a4, a6, h, N) : null;
  let RHS = rhs(a4, a6, N);
  if (h) RHS = FpX_rem(RHS, h, N);
  RHS = FpX_mulu(RHS, 4, N);
  const ff = h ? get_FpXQ_algebra(h, N) : get_FpX_algebra(N);
  const k = n + 2;
  const t: [Array<FpX | null>, Array<FpX | null>, Array<FpX | null>] = [
    new Array<FpX | null>(k + 1).fill(null),
    new Array<FpX | null>(k + 1).fill(null),
    new Array<FpX | null>(k + 1).fill(null),
  ];
  if (k >= 3 && D3) t[0][3] = D3;
  if (k >= 4 && D4) t[0][4] = D4;
  return { ff, t, r2: ff.sqr(RHS) };
}
/**
 * PARI `Fq_elldivpolmod` (ellsea.c:466-482): the `n`-division polynomial
 * modulo `h`.
 */
export function Fq_elldivpolmod(
  a4: bigint,
  a6: bigint,
  n: number,
  h: FpX | null,
  N: bigint
): FpX {
  const d = Fq_elldivpolmod_init(a4, a6, n, h, N);
  return divpol(d, n);
}
/** PARI `Fq_ellyn` (ellsea.c:492-512). */
function Fq_ellyn(d: DivPolMod, k: number): [FpX, FpX] {
  const ff = d.ff;
  if (k === 1) return [ff.one(), ff.one()];
  const pn2 = divpol(d, k - 2);
  const pp2 = divpol(d, k + 2);
  const pn12 = divpol_f2(d, k - 1);
  const pp12 = divpol_f2(d, k + 1);
  const on = ff.red(ff.sub(ff.mul(pp2, pn12), ff.mul(pn2, pp12)));
  const f = divpol(d, k);
  const f2 = divpol_f2(d, k);
  let f3 = ff.mul(f, f2);
  if (k % 2 === 0) f3 = ff.mul(f3, d.r2);
  return [on, f3];
}
/** PARI `Fq_elldivpol2` (ellsea.c:518-519). */
function Fq_elldivpol2(a4: bigint, a6: bigint, N: bigint): FpX {
  return FpX_red([Fq_mulu(a6, 4, N), Fq_mulu(a4, 4, N), 0n, 4n], N);
}
/** PARI `Fq_elldivpol2d` (ellsea.c:522-523). */
function Fq_elldivpol2d(a4: bigint, N: bigint): FpX {
  return FpX_red([Fq_mulu(a4, 2, N), 0n, 6n], N);
}
/** PARI `FqX_numer_isog_abscissa` (ellsea.c:526-542). */
function FqX_numer_isog_abscissa(h: FpX, a4: bigint, a6: bigint, N: bigint): FpX {
  const m = FpX_deg(h);
  const mp1 = h[m - 1] ?? 0n; /* PARI: gel(h, m+1), coefficient of x^(m-1) */
  const dh = FpX_deriv(h, N);
  const ddh = FpX_deriv(dh, N);
  const t = Fq_elldivpol2(a4, a6, N);
  const u = Fq_elldivpol2d(a4, N);
  const t1 = FpX_sub(FpX_sqr(dh, N), FpX_mul(ddh, h, N), N);
  const t2 = FpX_mul(u, FpX_mul(h, dh, N), N);
  const t3 = FpX_mul(FpX_sqr(h, N), deg1pol(BigInt(2 * m), Fq_mulu(mp1, 2, N), N), N);
  const f0 = FpX_add(FpX_sub(FpX_mul(t, t1, N), t2, N), t3, N);
  const t4 = FpX_mul(pol_x(), FpX_sqr(h, N), N);
  return FpX_add(t4, f0, N);
}

/* ================================================================== */
/* Kernel polynomial (ellsea.c:604-720, 845-861)                      */
/* ================================================================== */

/** PARI `FqX_mulhigh` (ellsea.c:611-616). */
function FqX_mulhigh(f: FpX, g: FpX, n2: number, n: number, N: bigint): FpX {
  const [fl, fh] = FpX_blocks2(f, n2);
  return FpX_add(FpX_mulhigh_i(fl, g, n2, N), FpXn_mul(fh, g, n - n2, N), N);
}
/** PARI `FqX_invlift1` (ellsea.c:618-623). */
function FqX_invlift1(Q: FpX, P: FpX, t1: number, t2: number, N: bigint): FpX {
  const H = FpXn_mul(FqX_mulhigh(Q, P, t1, t2, N), Q, t2 - t1, N);
  return FpX_sub(Q, FpX_shift(H, t1), N);
}
/** PARI `FqX_invsqrtlift1` (ellsea.c:625-631). */
function FqX_invsqrtlift1(Q: FpX, P: FpX, t1: number, t2: number, N: bigint): FpX {
  const D = FqX_mulhigh(P, FpX_sqr(Q, N), t1, t2, N);
  const H = FpXn_mul(Q, FpX_halve(D, N), t2 - t1, N);
  return FpX_sub(Q, FpX_shift(H, t1), N);
}
/** PARI `ZqX_integ2Xn` (ellsea.c:634-648). */
function ZqX_integ2Xn(P: FpX, n: number, N: bigint, pp: bigint, e: number): FpX | null {
  const d = FpX_deg(P);
  if (d === -1) return [];
  const Q = new Array<bigint>(d + 1);
  for (let k = 0; k <= d; k++) {
    const q = Zq_divu_safe(P[k]!, 2 * (k + n) + 1, N, pp, e);
    if (q === null) return null;
    Q[k] = q;
  }
  return FpX_renormalize(Q);
}
/**
 * PARI `Zq_Weierstrass` (ellsea.c:651-687): the solution of
 * `G*(S'^2) = (S/x)*(H o S)` mod `x^m`, `G = 1 + a4 x^2 + a6 x^3`,
 * `H = 1 + b4 x^2 + b6 x^3`.
 */
function Zq_Weierstrass(
  a4: bigint,
  a6: bigint,
  b4: bigint,
  b6: bigint,
  m: number,
  N: bigint,
  pp: bigint,
  n: number
): FpX | null {
  let mask = quadratic_prec_mask(m);
  let iGdS2: FpX = [1n];
  const G = FpX_red([1n, 0n, a4, a6], N);
  let GdS2 = G;
  let S: FpX = pol_x();
  let sG: FpX = [1n];
  let isG: FpX = sG;
  let dS: FpX = sG;
  let NN = 1;
  for (; mask > 1n; ) {
    const N2 = NN;
    NN <<= 1;
    if (mask & 1n) NN--;
    mask >>= 1n;
    const d = NN - N2;
    const S2 = FpX_sqr(S, N);
    let HS = FpX_Fp_add(FpX_Fp_mul(S, b6, N), b4, N);
    HS = FpX_Fp_add(FpXn_mul(S2, HS, NN, N), 1n, N);
    HS = FpXn_mul(HS, FpX_shift(S, -1), NN, N);
    sG = FpXn_mul(G, isG, N2, N);
    const dK = FpXn_mul(
      FpX_shift(FpX_sub(HS, GdS2, N), -N2),
      FpXn_mul(iGdS2, isG, d, N),
      d,
      N
    );
    const K = ZqX_integ2Xn(dK, N2, N, pp, n);
    if (!K) return null;
    const E = FpXn_mul(FpXn_mul(K, sG, d, N), dS, d, N);
    S = FpX_add(S, FpX_shift(E, N2 + 1), N);
    if (mask <= 1n) break;
    isG = FqX_invsqrtlift1(isG, G, N2, NN, N);
    dS = FpX_deriv(S, N);
    GdS2 = FpX_mul(G, FpX_sqr(dS, N), N);
    iGdS2 = FqX_invlift1(iGdS2, GdS2, N2, NN, N);
  }
  return S;
}
/** PARI `ZqXn_WNewton` (ellsea.c:689-715). */
function ZqXn_WNewton(
  S: FpX,
  l: number,
  a4: bigint,
  a6: bigint,
  pp1: bigint,
  N: bigint,
  pp: bigint,
  e: number
): FpX | null {
  const d = FpX_deg(S);
  const Ge = new Array<bigint>(d + 1).fill(0n);
  Ge[0] = mod(pp1, N);
  const S_ = (i: number) => S[i] ?? 0n; /* PARI gel(S, 2+i) */
  if (d >= 2) {
    const g = Zq_divu_safe(Fq_sub(S_(2), Fq_mulu(a4, l - 1, N), N), 6, N, pp, e);
    if (g === null) return null;
    Ge[1] = g;
  }
  if (d >= 3) {
    const g = Zq_divu_safe(
      Fq_sub(Fq_sub(S_(3), Fq_mul(a4, Fq_mulu(pp1, 6, N), N), N), Fq_mulu(a6, (l - 1) * 2, N), N),
      10,
      N,
      pp,
      e
    );
    if (g === null) return null;
    Ge[2] = g;
  }
  /* PARI: gel(Ge,k+1) = (gel(S,k+2) - a4 (4k-6) gel(Ge,k-1)
   *                                 - a6 (4k-8) gel(Ge,k-2)) / (4k-2)
   * with gel(P,i) the coefficient of degree i-2, i.e. in 0-based terms
   * Ge[k-1] = (S[k] - a4 (4k-6) Ge[k-3] - a6 (4k-8) Ge[k-4]) / (4k-2). */
  for (let k = 4; k <= d; k++) {
    const g = Zq_divu_safe(
      Fq_sub(
        Fq_sub(S_(k), Fq_mul(a4, Fq_mulu(Ge[k - 3]!, 4 * k - 6, N), N), N),
        Fq_mul(a6, Fq_mulu(Ge[k - 4]!, 4 * k - 8, N), N),
        N
      ),
      4 * k - 2,
      N,
      pp,
      e
    );
    if (g === null) return null;
    Ge[k - 1] = g;
  }
  return FpX_renormalize(Ge);
}
/** PARI `find_kernel` (ellsea.c:845-861). */
function find_kernel(
  a4: bigint,
  a6: bigint,
  l: number,
  b4: bigint,
  b6: bigint,
  pp1: bigint,
  N: bigint,
  pp: bigint,
  e: number
): FpX | null {
  const d = ((l + 1) >> 1) + 1;
  if (l === 3) return deg1pol(1n, Fq_neg(pp1, N), N);
  let S = Zq_Weierstrass(a4, a6, b4, b6, d + 1, N, pp, e);
  if (S === null) return null;
  S = FpX_shift(S, -1);
  const Sd = FpXn_inv(S, d, N);
  let Ge = ZqXn_WNewton(Sd, l, a4, a6, pp1, N, pp, e);
  if (!Ge) return null;
  Ge = FpX_neg(Ge, N);
  Ge = FpXn_expint(Ge, d, N);
  Ge = FpX_recip(FpX_red(Ge, pp));
  if (FpX_deg(Ge) === (l - 1) >> 1) return Ge;
  err_printf(`[find_kernel l=${l}: deg Ge = ${FpX_deg(Ge)} != ${(l - 1) >> 1}]\n`);
  return null;
}

/* ================================================================== */
/* Isogenous curve (ellsea.c:718-1060, 1164-1249)                     */
/* ================================================================== */

/** PARI `Fq_ellj` (ellsea.c:721-729). */
function Fq_ellj(a4: bigint, a6: bigint, N: bigint): bigint {
  const a43 = Fq_mulu(Fq_powu(a4, 3, N), 4, N);
  return Fq_div(Fq_mulu(a43, 1728, N), Fq_add(a43, Fq_mulu(Fq_sqr(a6, N), 27, N), N), N);
}
/** PARI `Zq_ellj` (ellsea.c:731-739). */
function Zq_ellj(a4: bigint, a6: bigint, N: bigint, pp: bigint, e: number): bigint {
  const a43 = Fq_mulu(Fq_powu(a4, 3, N), 4, N);
  return Zq_div(Fq_mulu(a43, 1728, N), Fq_add(a43, Fq_mulu(Fq_sqr(a6, N), 27, N), N), N, pp, e);
}
/** PARI `corr` (ellsea.c:1053-1060). */
function corr(c4: bigint, c6: bigint, N: bigint, pp: bigint, e: number): bigint {
  const c46 = Zq_div(Fq_sqr(c4, N), c6, N, pp, e);
  const c64 = Zq_div(c6, c4, N, pp, e);
  const a = Fq_div(2n, 3n, N);
  return Fq_add(Fq_halve(c46, N), Fq_mul(a, c64, N), N);
}
/** PARI `a4a6t_from_J` (ellsea.c:876-883). */
function a4a6t_from_J(
  l: number,
  C4t: bigint,
  C6t: bigint,
  N: bigint
): { a4t: bigint; a6t: bigint } {
  const l2 = mod(BigInt(l) * BigInt(l), N);
  const l4 = Fq_sqr(l2, N);
  const l6 = Fq_mul(l4, l2, N);
  const v = Fp_inv(mod(-864n, N), N);
  const u = Fq_mulu(v, 18, N);
  return { a4t: Fq_mul(C4t, Fq_mul(u, l4, N), N), a6t: Fq_mul(C6t, Fq_mul(v, l6, N), N) };
}
/** PARI `ZpX_liftroot` (Zp.c:809-829), `T = NULL`. */
function ZpX_liftroot(f: FpX, a0: bigint, p: bigint, e: number): bigint {
  let q = p;
  let a = mod(a0, q);
  if (e === 1) return a;
  let mask = quadratic_prec_mask(e);
  let fr = FpX_red(f, q);
  let W = Fp_inv(FpX_eval(FpX_deriv(fr, q), a, q), q);
  for (;;) {
    q = q * q;
    if (mask & 1n) q = q / p;
    mask >>= 1n;
    fr = FpX_red(f, q);
    a = Fq_sub(a, Fq_mul(W, FpX_eval(fr, a, q), q), q);
    if (mask === 1n) return a;
    W = Fq_sub(mod(W * 2n, q), Fq_mul(Fq_sqr(W, q), FpX_eval(FpX_deriv(fr, q), a, q), q), q);
  }
}

/**
 * PARI `find_isogenous_from_J` (ellsea.c:1164-1222): the curve `ell`-isogenous
 * to `E` and its kernel polynomial, from the classical modular equation.
 */
function find_isogenous_from_J(
  a4: bigint,
  a6: bigint,
  ell: number,
  MEQN: Meqn,
  g: bigint,
  pp: bigint,
  e: number
): { a4t: bigint; a6t: bigint; h: FpX } | null {
  const N = e === 1 ? pp : pp ** BigInt(e);
  const R = MEQN.evalR!;
  const dR = MEQN.evaldR!;
  const ddR = MEQN.evalddR!;
  if (mod(g, N) === 0n || Fq_sub(g, 1728n, N) === 0n) return null;
  const C4 = Fq_mul(a4, mod(-48n, N), N);
  const C6 = Fq_mul(a6, mod(-864n, N), N);
  if (C4 === 0n || C6 === 0n) return null;
  const j = Zq_ellj(a4, a6, N, pp, e);
  const jp = Fq_mul(j, Zq_div(C6, C4, N, pp, e), N);
  const co = corr(C4, C6, N, pp, e);
  const Py = FpX_deriv(R, N);
  const Pxy = FpX_deriv(dR, N);
  const Pyy = FpX_deriv(Py, N);
  const Pxj = FpX_eval(dR, g, N);
  if (Pxj === 0n) return null;
  const Pyj = FpX_eval(Py, g, N);
  const Pxxj = FpX_eval(ddR, g, N);
  const Pxyj = FpX_eval(Pxy, g, N);
  const Pyyj = FpX_eval(Pyy, g, N);
  const jtp = Fq_div(Fq_mul(jp, Zq_div(Pxj, Pyj, N, pp, e), N), mod(-BigInt(ell), N), N);
  const jtp2 = Fq_sqr(jtp, N);
  const jtp3 = Fq_mul(jtp, jtp2, N);
  const den = Fq_mul(Fq_sqr(g, N), Fq_sub(g, 1728n, N), N);
  const D = Zq_inv(den, pp, e);
  const C4t = Fq_mul(jtp2, Fq_mul(g, D, N), N);
  const C6t = Fq_mul(jtp3, D, N);
  const s0 = Fq_mul(Fq_sqr(jp, N), Pxxj, N);
  const s1 = Fq_mul(Fq_mulu(Fq_mul(jp, jtp, N), 2 * ell, N), Pxyj, N);
  const s2 = Fq_mul(Fq_mulu(jtp2, ell * ell, N), Pyyj, N);
  const s3 = Zq_div(Fq_add(s0, Fq_add(s1, s2, N), N), Fq_mul(jp, Pxj, N), N, pp, e);
  const cot = corr(C4t, C6t, N, pp, e);
  const c0 = Fq_sub(co, Fq_mulu(cot, ell, N), N);
  const p_1 = Fq_div(Fq_mulu(Fq_add(s3, c0, N), ell, N), mod(-4n, N), N);
  const { a4t, a6t } = a4a6t_from_J(ell, C4t, C6t, N);
  const h = find_kernel(a4, a6, ell, a4t, a6t, p_1, N, pp, e);
  if (!h) return null;
  return { a4t, a6t, h };
}

/**
 * PARI `find_isogenous` (ellsea.c:1224-1249).
 *
 * The `'C'` (canonical) and `'A'` (Atkin) modular-equation types only occur
 * when the `seadata` package is installed (ellsea.c:107-123); it is not
 * vendored under `reference/`, `get_modular_eqn` therefore always produces
 * type `'J'` and those two branches are unreachable.
 */
function find_isogenous(
  a4: bigint,
  a6: bigint,
  ell: number,
  MEQN: Meqn,
  g0: bigint,
  p: bigint
): { a4t: bigint; a6t: bigint; h: FpX } | null {
  const pp = p < 1n << 64n ? p : 0n;
  const e =
    pp !== 0n ? ulogint(((ell + 1) >> 1) + 1, pp) + ulogint(2 * ell + 4, pp) + 1 : 1;
  if (a4 === 0n || a6 === 0n) return null;
  let g = g0;
  if (e > 1) {
    const pe = p ** BigInt(e);
    const meqnj = meqn_j(MEQN, Zq_ellj(a4, a6, pe, p, e), ell, pe);
    g = ZpX_liftroot(meqnj, g, p, e);
  }
  if (MEQN.type !== 'J')
    throw new NotImplementedError(
      `SAGE_NOT_IMPLEMENTED: find_isogenous for modular equation type ` +
        `'${MEQN.type}' (find_isogenous_from_${MEQN.type === 'A' ? 'Atkin' : 'canonical'}, ` +
        `ellsea.c:${MEQN.type === 'A' ? '900-960' : '964-1051'}); those types only ` +
        `occur with the PARI seadata package, which is not vendored`
    );
  return find_isogenous_from_J(a4, a6, ell, MEQN, g, p, e);
}

/** PARI `FqX_homogenous_eval` (ellsea.c:1251-1261). */
function FqX_homogenous_eval(P: FpX, A: FpX, B: FpX, N: bigint): FpX {
  const d = FpX_deg(P);
  if (d < 0) return [];
  let s: FpX = FpX_red([P[d]!], N);
  let Bn: FpX = [1n];
  for (let i = d - 1; i >= 0; i--) {
    Bn = FpX_mul(Bn, B, N);
    s = FpX_add(FpX_mul(s, A, N), FpX_Fp_mul(Bn, P[i]!, N), N);
  }
  return s;
}
/** PARI `FqX_homogenous_div` (ellsea.c:1263-1275). */
function FqX_homogenous_div(
  P: FpX,
  Q: FpX,
  A: FpX,
  B: FpX,
  N: bigint
): { num: FpX; den: FpX } {
  const d = FpX_deg(Q) - FpX_deg(P);
  let num = FqX_homogenous_eval(P, A, B, N);
  let den = FqX_homogenous_eval(Q, A, B, N);
  if (d > 0) num = FpX_mul(num, FpX_powu(B, d, N), N);
  else if (d < 0) den = FpX_mul(den, FpX_powu(B, -d, N), N);
  return { num, den };
}
function FpX_powu(a: FpX, n: number, N: bigint): FpX {
  let r: FpX = [1n];
  let b = a;
  let e = n;
  while (e > 0) {
    if (e & 1) r = FpX_mul(r, b, N);
    e >>= 1;
    if (e > 0) b = FpX_sqr(b, N);
  }
  return r;
}

/** PARI `find_kernel_power` (ellsea.c:1277-1310). */
function find_kernel_power(
  Eba4: bigint,
  Eba6: bigint,
  Eca4: bigint,
  Eca6: bigint,
  ell: number,
  MEQN: Meqn,
  kpoly: FpX,
  Ib: { num: FpX; den: FpX },
  p: bigint
): {
  a4t: bigint;
  a6t: bigint;
  kpoly_new: FpX;
  gtmp: FpX;
  Ic: { num: FpX; den: FpX };
} | null {
  const num_iso = FqX_numer_isog_abscissa(kpoly, Eba4, Eba6, p);
  const mpoly = meqn_j(MEQN, Fq_ellj(Eca4, Eca6, p), ell, p);
  const mroots = FpX_roots(mpoly, p);
  const kpoly2 = FpX_sqr(kpoly, p);
  for (let i = 0; i < mroots.length; i++) {
    const tmp = find_isogenous(Eca4, Eca6, ell, MEQN, mroots[i]!, p);
    if (!tmp) return null;
    const { a4t, a6t, h: gtmp } = tmp;
    /* check that the kernel kpoly is the good one */
    const h = FqX_homogenous_eval(gtmp, num_iso, kpoly2, p);
    if (Fq_elldivpolmod(Eba4, Eba6, ell, h, p).length !== 0) {
      const Ic = FqX_homogenous_div(num_iso, kpoly2, Ib.num, Ib.den, p);
      const kpoly_new = FqX_homogenous_eval(gtmp, Ic.num, Ic.den, p);
      return { a4t, a6t, kpoly_new, gtmp, Ic };
    }
  }
  return null;
}

/* ================================================================== */
/* Eigenvalue of the Frobenius (ellsea.c:741-835)                     */
/* ================================================================== */

/** PARI `Fq_find_eigen_Frobenius` (ellsea.c:750-757). */
function Fq_find_eigen_Frobenius(a4: bigint, a6: bigint, h: FpX, p: bigint): FpX {
  const RHS = FpX_rem(rhs(a4, a6, p), h, p);
  return FpXQ_pow(RHS, p >> 1n, h, p);
}
/** PARI `find_eigen_value_oneroot` (ellsea.c:762-778). */
function find_eigen_value_oneroot(
  a4: bigint,
  a6: bigint,
  ell: number,
  tr: number[],
  h: FpX,
  p: bigint
): number {
  const Gy = Fq_find_eigen_Frobenius(a4, a6, h, p);
  const L = BigInt(ell);
  let t = Number(mod(BigInt(tr[0]!) * Fp_inv(2n, L), L));
  if (t < ell >> 1) t = ell - t;
  const d = Fq_elldivpolmod_init(a4, a6, t, h, p);
  const f = Fq_ellyn(d, t);
  const Dy = FpXQ_mul(Gy, f[1], h, p);
  if (!FpX_equal(f[0], Dy)) t = ell - t;
  return t;
}
function FpX_equal(a: FpX, b: FpX): boolean {
  if (a.length !== b.length) return false;
  for (let i = 0; i < a.length; i++) if (a[i] !== b[i]) return false;
  return true;
}
/** PARI `Fq_find_eigen_value_power` (ellsea.c:808-829). */
function find_eigen_value_power(
  a4: bigint,
  a6: bigint,
  ell: number,
  k: number,
  lambda: number,
  h: FpX,
  p: bigint
): number {
  const ellk1 = Math.pow(ell, k - 1);
  const ellk = ell * ellk1;
  const Gy = Fq_find_eigen_Frobenius(a4, a6, h, p);
  const d = Fq_elldivpolmod_init(a4, a6, ellk, h, p);
  let t = lambda;
  for (; t < ellk; t += ellk1) {
    const f = Fq_ellyn(d, t);
    const Dr = FpXQ_mul(Gy, f[1], h, p);
    if (FpX_equal(f[0], Dr)) break;
    if (FpX_equal(f[0], FpX_neg(Dr, p))) {
      t = ellk - t;
      break;
    }
  }
  return t;
}

/* ================================================================== */
/* CM curves (ellsea.c:159-224)                                       */
/* ================================================================== */

/** PARI `factoru`, small `n`. */
function factoru(n: number): Array<[number, number]> {
  const f: Array<[number, number]> = [];
  let m = n;
  for (let d = 2; d * d <= m; d++) {
    if (m % d) continue;
    let e = 0;
    while (m % d === 0) {
      m /= d;
      e++;
    }
    f.push([d, e]);
  }
  if (m > 1) f.push([m, 1]);
  return f;
}
/** PARI `divisorsu_fact`. */
function divisorsu(n: number): number[] {
  const d: number[] = [1];
  for (const [q, e] of factoru(n)) {
    const cur = d.slice();
    let qk = 1;
    for (let i = 1; i <= e; i++) {
      qk *= q;
      for (const x of cur) d.push(x * qk);
    }
  }
  return d.sort((a, b) => a - b);
}
/**
 * PARI `coredisc2u_fact` for a positive `n`: writes `-n = d f^2` with `d` a
 * fundamental discriminant, and returns `|d|` together with the factorisation
 * of `f`.  (ellsea.c uses `coredisc2u_fact(F, -1, &P, &E)`.)
 */
function coredisc2u_neg(n: number): { d: number; f: number[] } {
  /* -n = D * f^2 with D fundamental (D < 0) */
  const fa = factoru(n);
  let core = 1;
  let f = 1;
  for (const [q, e] of fa) {
    if (e & 1) core *= q;
    f *= Math.pow(q, e >> 1);
  }
  /* -core is squarefree; D = -core if -core = 1 mod 4, else -4core */
  let D = -core;
  if (((D % 4) + 4) % 4 !== 1) {
    D = -4 * core;
    /* f must halve */
    if (f % 2 === 0) f = f / 2;
    else {
      /* -n = D f^2 forces f even; recompute */
      D = -core;
    }
  }
  return { d: Math.abs(D), f: divisorsu(f) };
}
/**
 * PARI `list_singular_discs` (ellsea.c:172-193): the discriminants `D` with
 * `polclass(D) | poldisc(polmodular(l))`, as a bit set.
 */
function list_singular_discs(l: number): Set<number> {
  const _4l2 = 4 * l * l;
  const V = new Set<number>();
  V.add(4); /* v = 0 */
  V.add(3); /* v = l */
  for (let v = 1; v < 2 * l; v++) {
    if (v === l) continue;
    const { d, f } = coredisc2u_neg(_4l2 - v * v);
    for (const c of f) V.add(d * c * c);
  }
  return V;
}
/** PARI `find_CM` (ellsea.c:195-212). */
function find_CM(l: number, j: bigint, p: bigint): number {
  const inv = INV_J;
  const v = list_singular_discs(l);
  const n = 4 * l * l;
  const db = polmodular_db_init(inv);
  for (let i = 1; i < n; i++) {
    if (!v.has(i)) continue;
    const C = polclass0(-i, inv, db);
    const F = FpX_eval(FpX_red(C, p), j, p);
    if (F === 0n) return -i;
  }
  return 0;
}
/** PARI `Fq_ellcard_CM` (ellsea.c:220-232), `T = NULL` so `d = 1`, `q = p`. */
function Fq_ellcard_CM(disc: number, a4: bigint, a6: bigint, p: bigint): bigint {
  const q1 = p + 1n;
  /* upstream passes the factorisation of 4 p^d = 2^2 p directly */
  const Q = qfbsolve(Qfb(1n, 0n, BigInt(-disc)), 4n * p, 3, [
    [2n, 2n],
    [p, 1n],
  ]) as Array<[bigint, bigint]>;
  if (Q.length === 0) return q1;
  const S = Q.map(([x, y]) => q1 + (y > 0n ? x : -x));
  return gen_select_order(S, a4, a6, p);
}

/* ================================================================== */
/* Elliptic curve group operations over Fp (Fle.c), used by the        */
/* match-and-sort stage; kept local so that `ellsea` does not create a */
/* module cycle with `group.ts`.                                       */
/* ================================================================== */

type Fle = { x: bigint; y: bigint } | null;

function Fle_dbl(P: Fle, a4: bigint, p: bigint): Fle {
  if (!P || P.y === 0n) return null;
  const lam = Fq_div(Fq_add(Fq_mulu(Fq_sqr(P.x, p), 3, p), a4, p), Fq_mulu(P.y, 2, p), p);
  const x = Fq_sub(Fq_sqr(lam, p), Fq_mulu(P.x, 2, p), p);
  return { x, y: Fq_sub(Fq_mul(lam, Fq_sub(P.x, x, p), p), P.y, p) };
}
function Fle_add(P: Fle, Q: Fle, a4: bigint, p: bigint): Fle {
  if (!P) return Q;
  if (!Q) return P;
  if (P.x === Q.x) {
    if (P.y === Q.y) return Fle_dbl(P, a4, p);
    return null;
  }
  const lam = Fq_div(Fq_sub(Q.y, P.y, p), Fq_sub(Q.x, P.x, p), p);
  const x = Fq_sub(Fq_sub(Fq_sqr(lam, p), P.x, p), Q.x, p);
  return { x, y: Fq_sub(Fq_mul(lam, Fq_sub(P.x, x, p), p), P.y, p) };
}
function Fle_neg(P: Fle, p: bigint): Fle {
  return P ? { x: P.x, y: Fq_neg(P.y, p) } : null;
}
function Fle_mul(P: Fle, n: bigint, a4: bigint, p: bigint): Fle {
  if (n === 0n || !P) return null;
  let m = n;
  let Q: Fle = P;
  if (m < 0n) {
    m = -m;
    Q = Fle_neg(P, p);
  }
  let R: Fle = null;
  let base = Q;
  while (m > 0n) {
    if (m & 1n) R = Fle_add(R, base, a4, p);
    m >>= 1n;
    if (m > 0n) base = Fle_dbl(base, a4, p);
  }
  return R;
}
function random_Fle(a4: bigint, a6: bigint, p: bigint): Fle {
  for (;;) {
    const x = randomFp(p);
    const rhsv = mod(x * (x * x + a4) + a6, p);
    if (rhsv === 0n && mod(3n * x * x + a4, p) === 0n) continue;
    if (kronecker(rhsv, p) < 0) continue;
    const y = Fp_sqrt(rhsv, p);
    if (y === null) throw new PariPrimeError('random_Fle', p);
    return { x, y };
  }
}
/** PARI `gen_select_order` (bb_group.c:762-794) for `FpE`. */
function gen_select_order(o: bigint[], a4: bigint, a6: bigint, p: bigint): bigint {
  const lo = o.length;
  if (lo === 0) throw new PariBugError('gen_select_order (no candidate)');
  if (lo === 1) return o[0]!;
  const so = o.map((_, i) => i).sort((i, j) => (o[i]! < o[j]! ? -1 : o[i]! > o[j]! ? 1 : 0));
  const vo = new Array<boolean>(lo).fill(false);
  let nbo = lo;
  let lastgood = o[so[lo - 1]!]!;
  for (;;) {
    let lasto = 0n;
    const P = random_Fle(a4, a6, p);
    let t: Fle = null;
    for (let i = 0; i < lo; i++) {
      const newo = o[so[i]!]!;
      if (vo[i]) continue;
      t = Fle_add(t, Fle_mul(P, newo - lasto, a4, p), a4, p);
      lasto = newo;
      if (t !== null) {
        if (--nbo === 1) return lastgood;
        vo[i] = true;
      } else lastgood = lasto;
    }
  }
}

/* ================================================================== */
/* Supersingularity (FpE.c:700-805)                                   */
/* ================================================================== */

/** F_{p^2} = F_p[w]/(w^2 - g); elements are `[a, b]` for `a + b w`. */
type Fp2 = [bigint, bigint];

function Fp2_mul(x: Fp2, y: Fp2, g: bigint, p: bigint): Fp2 {
  const a = mod(x[0] * y[0] + ((x[1] * y[1]) % p) * g, p);
  const b = mod(x[0] * y[1] + x[1] * y[0], p);
  return [a, b];
}
function Fp2_sqr(x: Fp2, g: bigint, p: bigint): Fp2 {
  return Fp2_mul(x, x, g, p);
}
function Fp2_sub(x: Fp2, y: Fp2, p: bigint): Fp2 {
  return [mod(x[0] - y[0], p), mod(x[1] - y[1], p)];
}
function Fp2_mulu(x: Fp2, n: bigint, p: bigint): Fp2 {
  return [mod(x[0] * n, p), mod(x[1] * n, p)];
}
function Fp2_halve(x: Fp2, p: bigint): Fp2 {
  const h = Fp_inv(2n, p);
  return [mod(x[0] * h, p), mod(x[1] * h, p)];
}
function Fp2_pow(x: Fp2, n: bigint, g: bigint, p: bigint): Fp2 {
  let r: Fp2 = [1n, 0n];
  let b = x;
  let e = n;
  while (e > 0n) {
    if (e & 1n) r = Fp2_mul(r, b, g, p);
    e >>= 1n;
    if (e > 0n) b = Fp2_sqr(b, g, p);
  }
  return r;
}
function Fp2_equal(x: Fp2, y: Fp2): boolean {
  return x[0] === y[0] && x[1] === y[1];
}
/** Tonelli-Shanks in F_{p^2} (PARI's `FpXQ_sqrt` -> `gen_Shanks_sqrtn`). */
function Fp2_sqrt(a: Fp2, g: bigint, p: bigint): Fp2 | null {
  if (a[0] === 0n && a[1] === 0n) return [0n, 0n];
  const q = p * p;
  let Q = q - 1n;
  let e = 0;
  while ((Q & 1n) === 0n) {
    Q >>= 1n;
    e++;
  }
  if (!Fp2_equal(Fp2_pow(a, (q - 1n) / 2n, g, p), [1n, 0n])) return null;
  /* find a non-residue */
  let z: Fp2 = [0n, 0n];
  for (let c = 1n; ; c++) {
    const cand: Fp2 = [c, 1n];
    if (Fp2_equal(Fp2_pow(cand, (q - 1n) / 2n, g, p), [mod(-1n, p), 0n])) {
      z = Fp2_pow(cand, Q, g, p);
      break;
    }
  }
  let M = e;
  let c = z;
  let t = Fp2_pow(a, Q, g, p);
  let R = Fp2_pow(a, (Q + 1n) / 2n, g, p);
  while (!Fp2_equal(t, [1n, 0n])) {
    let i = 0;
    let t2 = t;
    while (!Fp2_equal(t2, [1n, 0n])) {
      t2 = Fp2_sqr(t2, g, p);
      i++;
      if (i === M) return null;
    }
    let b = c;
    for (let k = 0; k < M - i - 1; k++) b = Fp2_sqr(b, g, p);
    M = i;
    c = Fp2_sqr(b, g, p);
    t = Fp2_mul(t, c, g, p);
    R = Fp2_mul(R, b, g, p);
  }
  return R;
}
/** PARI `Fp_2gener` (arith1.c:1067): a generator of the 2-Sylow of `Fp^*`. */
function Fp_2gener(p: bigint): bigint {
  let q = p - 1n;
  let e = 0;
  while ((q & 1n) === 0n) {
    q >>= 1n;
    e++;
  }
  if (e === 0) throw new PariBugError('Fp_2gener');
  for (let x = 2n; ; x++) {
    if (kronecker(x, p) < 0) return Fp_pow(x, q, p);
  }
}
/** PARI `FpX_quad_root` (FpX.c). */
function FpX_quad_root(x: FpX, p: bigint): bigint | null {
  const b = x[1] ?? 0n;
  const c = x[0] ?? 0n;
  const D = mod(b * b - 4n * c, p);
  if (kronecker(D, p) === -1) return null;
  const s = Fp_sqrt(D, p);
  if (s === null) return null;
  return Fq_halve(Fq_sub(s, b, p), p);
}
/** PARI `FqX_quad_root` over F_{p^2}. */
function Fp2X_quad_root(x: Fp2[], g: bigint, p: bigint): Fp2 | null {
  const b = x[1] ?? [0n, 0n];
  const c = x[0] ?? [0n, 0n];
  const D = Fp2_sub(Fp2_sqr(b, g, p), Fp2_mulu(c, 4n, p), p);
  const s = Fp2_sqrt(D, g, p);
  if (!s) return null;
  return Fp2_halve(Fp2_sub(s, b, p), p);
}
/** Evaluate the (symmetric) `Phi_2` at `x = j`, over F_p or F_{p^2}. */
function Phi2_evalx_Fp(Phi2: FpXY, j: bigint, p: bigint): FpX {
  return FpXY_evaly(Phi2, j, p);
}
function Phi2_evalx_Fp2(Phi2: FpXY, j: Fp2, g: bigint, p: bigint): Fp2[] {
  /* Horner over the outer variable, coefficients are FpX in the inner one */
  const out: Fp2[] = [];
  const n = Math.max(...Phi2.map((c) => c.length), 0);
  for (let b = 0; b < n; b++) {
    let s: Fp2 = [0n, 0n];
    for (let a = Phi2.length - 1; a >= 0; a--) {
      s = Fp2_mul(s, j, g, p);
      const c = Phi2[a]![b] ?? 0n;
      s = [mod(s[0] + c, p), s[1]];
    }
    out.push(s);
  }
  while (out.length && out[out.length - 1]![0] === 0n && out[out.length - 1]![1] === 0n)
    out.pop();
  return out;
}
function Fp2X_div_by_X_x(a: Fp2[], x: Fp2, g: bigint, p: bigint): Fp2[] {
  const l = a.length;
  if (l <= 1) return [];
  const z = new Array<Fp2>(l - 1);
  z[l - 2] = a[l - 1]!;
  for (let i = l - 3; i >= 0; i--) {
    const t = Fp2_mul(x, z[i + 1]!, g, p);
    z[i] = [mod(a[i + 1]![0] + t[0], p), mod(a[i + 1]![1] + t[1], p)];
  }
  return z;
}
/** PARI `Fq_path_extends_to_floor` (FpE.c:698-721), over F_{p^2}. */
function Fq_path_extends_to_floor(
  j_prev0: Fp2,
  j0: Fp2,
  g: bigint,
  p: bigint,
  Phi2: FpXY,
  max_len: number
): boolean {
  let j_prev = j_prev0;
  let j = j0;
  for (let d = 1; d <= max_len; d++) {
    const Phi2_j = Fp2X_div_by_X_x(Phi2_evalx_Fp2(Phi2, j, g, p), j_prev, g, p);
    const j_next = Fp2X_quad_root(Phi2_j, g, p);
    if (!j_next) return true;
    j_prev = j;
    j = j_next;
  }
  return false;
}
/** PARI `Fp_path_extends_to_floor` (FpE.c:723-750). */
function Fp_path_extends_to_floor(
  j_prev0: bigint[],
  j0: bigint[],
  p: bigint,
  Phi2: FpXY,
  max_len: number
): { found: boolean; j: bigint; j_prev: bigint } {
  const j_prev = j_prev0.slice();
  const j = j0.slice();
  const l = j.length;
  for (let d = 1; d <= max_len; d++) {
    for (let i = 0; i < l; i++) {
      const Phi2_j = FpX_div_by_X_x(Phi2_evalx_Fp(Phi2, j[i]!, p), j_prev[i]!, p).q;
      const j_next = FpX_quad_root(Phi2_j, p);
      if (j_next === null) return { found: true, j: j[i]!, j_prev: j_prev[i]! };
      j_prev[i] = j[i]!;
      j[i] = j_next;
    }
  }
  return { found: false, j: 0n, j_prev: 0n };
}
/** PARI `Fp_jissupersingular` (FpE.c:753-775). */
function Fp_jissupersingular(j: bigint, p: bigint): boolean {
  const max_path_len = p.toString(2).length; /* expi(p)+1 */
  const Phi2raw = polmodular_ZXX(2, INV_J);
  const Phi2: FpXY = FpXY_renormalize(Phi2raw.map((col) => FpX_red(col, p)));
  const Phi2_j = Phi2_evalx_Fp(Phi2, j, p);
  const roots = FpX_roots(Phi2_j, p);
  const nbroots = roots.length;
  if (nbroots === 0) return false;
  /* upstream: S = deg2pol(1, 0, -Fp_2gener(p)), i.e. F_p[w]/(w^2 - g) */
  const g = Fp_2gener(p);
  let jj: Fp2;
  let jprev: Fp2;
  if (nbroots === 1 && FpX_is_squarefree(Phi2_j, p)) {
    jprev = [j, 0n];
    const q = FpX_div_by_X_x(Phi2_j, roots[0]!, p).q;
    const r = Fp2X_quad_root(
      q.map((c) => [c, 0n] as Fp2),
      g,
      p
    );
    if (!r) return false;
    jj = r;
  } else {
    const r = Fp_path_extends_to_floor(
      new Array<bigint>(nbroots).fill(j),
      roots,
      p,
      Phi2,
      max_path_len
    );
    if (!r.found) return true;
    jj = [r.j, 0n];
    jprev = [r.j_prev, 0n];
  }
  return !Fq_path_extends_to_floor(jprev, jj, g, p, Phi2, max_path_len);
}
/** PARI `Fp_elljissupersingular` (FpE.c:794-805). */
export function Fp_elljissupersingular(j: bigint, p: bigint): boolean {
  if (p <= 5n) return mod(j, p) === 0n;
  const CM = Fp_ellj_get_CM(mod(j, p), 1n, p);
  if (CM < 0) return kronecker(BigInt(CM), p) < 0;
  return Fp_jissupersingular(mod(j, p), p);
}

/* ================================================================== */
/* TRACE (ellsea.c:1312-1517)                                         */
/* ================================================================== */

/** PARI `enum mod_type` (ellsea.c:1315). */
enum ModType {
  MTcm = 0,
  MTpathological = 1,
  MTAtkin = 2,
  MTElkies = 3,
  MTone_root = 4,
  MTroots = 5,
}

/** PARI `Fp_study_eqn` (ellsea.c:1328-1337). */
function Fp_study_eqn(mpoly: FpX, p: bigint): { g: bigint | null; dG: number; r: number } {
  const T = FpX_normalize(mpoly, p);
  const XP = FpX_Frobenius(T, p);
  const G = FpX_gcd(FpX_sub(XP, pol_x(), p), T, p);
  const dG = FpX_deg(G);
  if (dG <= 0) return { g: null, dG: 0, r: FpX_ddf_degree(T, XP, p) };
  return { g: FpX_oneroot(G, p), dG, r: 0 };
}

/** PARI `study_modular_eqn` (ellsea.c:1364-1392). */
function study_modular_eqn(
  ell: number,
  mpoly: FpX,
  p: bigint
): { g: bigint | null; mt: ModType; r: number } {
  if (!FpX_is_squarefree(mpoly, p)) return { g: null, mt: ModType.MTcm, r: 0 };
  const { g, dG, r } = Fp_study_eqn(mpoly, p);
  let mt: ModType;
  switch (dG) {
    case 0:
      mt = ModType.MTAtkin;
      break;
    case 1:
      mt = ModType.MTone_root;
      break;
    case 2:
      mt = ModType.MTElkies;
      break;
    default:
      mt = dG === ell + 1 ? ModType.MTroots : ModType.MTpathological;
  }
  return { g, mt, r };
}

/** PARI `Fl_div`, `Fl_add`, ... on small moduli, as plain numbers. */
function Fl_inv(a: number, m: number): number {
  let [old_r, r] = [((a % m) + m) % m, m];
  let [old_s, s] = [1, 0];
  while (r !== 0) {
    const q = Math.floor(old_r / r);
    [old_r, r] = [r, old_r - q * r];
    [old_s, s] = [s, old_s - q * s];
  }
  if (old_r !== 1) throw new Error(`impossible inverse modulo ${m}: ${a}`);
  return ((old_s % m) + m) % m;
}
function Fl_mul(a: number, b: number, m: number): number {
  return Number((BigInt(a) * BigInt(b)) % BigInt(m));
}
function Fl_div(a: number, b: number, m: number): number {
  return Fl_mul(a, Fl_inv(b, m), m);
}
function Fl_add(a: number, b: number, m: number): number {
  return (a + b) % m;
}
function Fl_sub(a: number, b: number, m: number): number {
  return ((a - b) % m + m) % m;
}
function Fl_sqrt(a: number, m: number): number {
  const r = Fp_sqrt(BigInt(a), BigInt(m));
  if (r === null) throw new Error(`Fl_sqrt: ${a} is not a square mod ${m}`);
  return Number(r);
}

/**
 * PARI `find_trace_Elkies_power` (ellsea.c:1396-1451): the trace modulo
 * `ell^k` when `ell` is an Elkies prime.
 */
function find_trace_Elkies_power(
  a4: bigint,
  a6: bigint,
  ell: number,
  ptk: { k: number },
  MEQN: Meqn,
  g: bigint,
  tr: number[] | null,
  q: bigint,
  p: bigint,
  smallfact: number
): number[] | null {
  let k = ptk.k;
  let ellk = Math.pow(ell, k);
  let pellk = Number(mod(q, BigInt(ellk)));
  let Eba4 = a4;
  let Eba6 = a6;
  const tmp = find_isogenous(a4, a6, ell, MEQN, g, p);
  if (!tmp) return null;
  let Eca4 = tmp.a4t;
  let Eca6 = tmp.a6t;
  let kpoly = tmp.h;
  let Ib: { num: FpX; den: FpX } = { num: pol_x(), den: [1n] };
  let lambda = tr
    ? find_eigen_value_oneroot(a4, a6, ell, tr, kpoly, p)
    : find_eigen_value_power(a4, a6, ell, 1, 1, kpoly, p);
  if (smallfact && smallfact % ell !== 0) {
    const pell = pellk % ell;
    const ap = Fl_add(lambda, Fl_div(pell, lambda, ell), ell);
    if (Fl_sub(pell, ap, ell) === ell - 1) return [ap];
    if (smallfact < 0 && Fl_add(pell, ap, ell) === ell - 1) return [ap];
  }
  for (let cnt = 2; cnt <= k; cnt++) {
    const t2 = find_kernel_power(Eba4, Eba6, Eca4, Eca6, ell, MEQN, kpoly, Ib, p);
    if (!t2) {
      k = cnt - 1;
      break;
    }
    lambda = find_eigen_value_power(a4, a6, ell, cnt, lambda, t2.kpoly_new, p);
    Eba4 = Eca4;
    Eba6 = Eca6;
    Eca4 = t2.a4t;
    Eca6 = t2.a6t;
    kpoly = t2.gtmp;
    Ib = t2.Ic;
  }
  ellk = Math.pow(ell, k);
  pellk = Number(mod(q, BigInt(ellk)));
  ptk.k = k;
  return [Fl_add(lambda, Fl_div(pellk, lambda, ellk), ellk)];
}

/**
 * PARI `find_trace_Atkin` (ellsea.c:1455-1485): the possible traces when `ell`
 * is an Atkin prime and the modular equation has splitting degree `r`.
 */
function find_trace_Atkin(ell: number, r: number, q: bigint): number[] {
  const val_pos: number[] = [];
  const pell = Number(mod(q, BigInt(ell)));
  const invp = Fl_inv(pell, ell);
  const P = factoru(r).map(([q0]) => q0);
  /* arithmetic in F_ell[x]/(x^2 - teta x + pell) */
  const mulmod = (a: [number, number], b: [number, number], teta: number): [number, number] => {
    /* (a0 + a1 x)(b0 + b1 x) with x^2 = teta x - pell */
    const c0 = Fl_mul(a[0], b[0], ell);
    const c1 = Fl_add(Fl_mul(a[0], b[1], ell), Fl_mul(a[1], b[0], ell), ell);
    const c2 = Fl_mul(a[1], b[1], ell);
    return [
      Fl_sub(c0, Fl_mul(c2, pell, ell), ell),
      Fl_add(c1, Fl_mul(c2, teta, ell), ell),
    ];
  };
  const powmod = (a: [number, number], n: number, teta: number): [number, number] => {
    let res: [number, number] = [1, 0];
    let b = a;
    let e = n;
    while (e > 0) {
      if (e & 1) res = mulmod(res, b, teta);
      e >>= 1;
      if (e > 0) b = mulmod(b, b, teta);
    }
    return res;
  };
  const is1 = (a: [number, number]) => a[0] === 1 && a[1] === 0;
  if (r === 2 && kronecker(BigInt(ell - pell), BigInt(ell)) < 0) val_pos.push(0);
  for (let teta = 1; teta < ell; teta++) {
    const disc = Fl_sub(Fl_mul(teta, teta, ell), Fl_mul(4, pell, ell), ell);
    if (kronecker(BigInt(disc), BigInt(ell)) >= 0) continue;
    const U: [number, number] = [ell - 1, Fl_mul(invp, teta, ell)];
    const a = powmod(U, Math.floor(r / P[0]!), teta);
    if (!is1(a) && is1(powmod(a, P[0]!, teta))) {
      let i = 1;
      for (; i < P.length; i++) if (is1(powmod(U, Math.floor(r / P[i]!), teta))) break;
      if (i === P.length) val_pos.push(teta);
    }
  }
  return val_pos;
}

/** PARI `find_trace_one_root` (ellsea.c:1488-1493). */
function find_trace_one_root(ell: number, q: bigint): number[] {
  const s = Fl_sqrt(Number(mod(q, BigInt(ell))), ell);
  const a = Fl_add(s, s, ell);
  return [a, ell - a];
}
/** PARI `find_trace_lp1_roots` (ellsea.c:1495-1502). */
function find_trace_lp1_roots(ell: number, q: bigint): number[] {
  const ell2 = ell * ell;
  const pell = Number(mod(q, BigInt(ell2)));
  const a = Fl_sqrt(pell % ell, ell);
  const pa = Fl_add(Fl_div(pell, a, ell2), a, ell2);
  return [pa, ell2 - pa];
}

/**
 * PARI `find_trace` (ellsea.c:1505-1573): the trace modulo `ell^kt`, as a list
 * of possibilities.  Returns `{card}` in the CM case (upstream returns a
 * t_INT), `null` when the prime must be skipped.
 */
function find_trace(
  a4: bigint,
  a6: bigint,
  j: bigint,
  ell: number,
  q: bigint,
  p: bigint,
  ptkt: { kt: number },
  smallfact: number
): { tr: number[] } | { card: bigint } | null {
  let kt = Math.max(
    Math.floor(Math.log(expi(q) * Math.LN2) / Math.log(ell)),
    1
  );
  err_printf(`SEA: Prime ${ell} `);
  const t_eqn = Date.now();
  const MEQN = get_modular_eqn(ell);
  const t_meqn = Date.now();
  const meqnj = meqn_j(MEQN, j, ell, p);
  const { g, mt, r } = study_modular_eqn(ell, meqnj, p);
  err_printf(
    `[eqn ${t_meqn - t_eqn}ms study ${Date.now() - t_meqn}ms] ` +
      ['CM', 'Pathological', 'Atkin', 'Elkies', 'One root', 'l+1 roots'][mt] +
      '\t'
  );
  let tr: number[] | null;
  switch (mt) {
    case ModType.MTone_root: {
      const tr2 = find_trace_one_root(ell, q);
      const ptk = { k: kt };
      tr = find_trace_Elkies_power(a4, a6, ell, ptk, MEQN, g!, tr2, q, p, smallfact);
      kt = ptk.k;
      if (!tr) {
        tr = tr2;
        kt = 1;
      }
      break;
    }
    case ModType.MTElkies: {
      /* Contrary to MTone_root, may look mod higher powers of ell */
      if (p <= BigInt(2 * ell + 3)) kt = 1;
      const ptk = { k: kt };
      tr = find_trace_Elkies_power(a4, a6, ell, ptk, MEQN, g!, null, q, p, smallfact);
      kt = ptk.k;
      if (!tr) return null;
      break;
    }
    case ModType.MTroots:
      tr = find_trace_lp1_roots(ell, q);
      kt = 2;
      break;
    case ModType.MTAtkin:
      tr = find_trace_Atkin(ell, r, q);
      if (tr.length === 0) throw new PariPrimeError('ellap', p);
      kt = 1;
      break;
    case ModType.MTcm: {
      const D = find_CM(ell, j, p);
      return { card: Fq_ellcard_CM(D, a4, a6, p) };
    }
    default: /* MTpathological */
      return null;
  }
  ptkt.kt = kt;
  err_printf(`${tr.length} trace(s) [${Date.now() - t_eqn}ms]\n`);
  return { tr };
}

/* ================================================================== */
/* MATCH AND SORT (ellsea.c:1576-1958)                                */
/* ================================================================== */

/** `floor(log2(n))`, PARI `expi`. */
function expi(n: bigint): number {
  return n.toString(2).length - 1;
}
/** PARI `sqrti`: `floor(sqrt(n))`. */
function sqrti(n: bigint): bigint {
  if (n < 0n) throw new Error('sqrti: negative argument');
  if (n < 2n) return n;
  let x = n;
  let y = (x + 1n) / 2n;
  while (y < x) {
    x = y;
    y = (x + n / x) / 2n;
  }
  return x;
}
/** Floor division for bigints (PARI `truedivii`). */
function truedivii(a: bigint, b: bigint): bigint {
  const q = a / b;
  return a - q * b !== 0n && a < 0n !== b < 0n ? q - 1n : q;
}

/** One entry of `compile_atkin`: the modulus and the possible traces. */
interface AtkinEntry {
  mod: bigint;
  traces: number[];
}

/** PARI `separation` (ellsea.c:1576-1604). */
function separation(cnt: number[]): number {
  const k = cnt.length;
  const l = (1 << k) - 1;
  let P = 1n;
  for (let j = 0; j < k; j++) P *= BigInt(cnt[j]!);
  const P3 = 3n * P;
  let best_i = 0;
  let best_r = P3;
  for (let i = 1; i < l; i++) {
    let p_b = 1n;
    for (let j = 0; j < k; j++) if (i & (1 << j)) p_b *= BigInt(cnt[j]!);
    const r = (p_b * p_b) * 4n - P3;
    if (r === 0n) return i;
    const ar = r < 0n ? -r : r;
    const abr = best_r < 0n ? -best_r : best_r;
    if (ar < abr) {
      best_i = i;
      best_r = r;
    }
  }
  return best_i;
}

/** PARI `multiple_crt` (ellsea.c:1610-1633). */
function multiple_crt(x: bigint[], y: number[], q: bigint, P: bigint): { x: bigint[]; P: bigint } {
  /* bezout(P, q, &u, &v): u P + v q = 1 */
  const { u, v } = bezout(P, q);
  const a1 = P * u;
  const a2 = q * v;
  const PQ = P * q;
  const out = new Array<bigint>(x.length * y.length);
  let k = 0;
  for (let i = 0; i < x.length; i++) {
    const a2x = mod(a2 * x[i]!, PQ);
    for (let j = 0; j < y.length; j++) out[k++] = mod(a1 * BigInt(y[j]!) + a2x, PQ);
  }
  return { x: out, P: PQ };
}
function bezout(a: bigint, b: bigint): { u: bigint; v: bigint } {
  let [old_r, r] = [a, b];
  let [old_s, s] = [1n, 0n];
  let [old_t, t] = [0n, 1n];
  while (r !== 0n) {
    const q = old_r / r;
    [old_r, r] = [r, old_r - q * r];
    [old_s, s] = [s, old_s - q * s];
    [old_t, t] = [t, old_t - q * t];
  }
  if (old_r < 0n) return { u: -old_s, v: -old_t };
  return { u: old_s, v: old_t };
}

/** PARI `possible_traces` (ellsea.c:1640-1669). */
function possible_traces(
  compile: AtkinEntry[],
  mask: number
): { V: bigint[]; P: bigint } {
  const C: AtkinEntry[] = [];
  for (let i = 0; i < compile.length; i++) if (mask & (1 << i)) C.push(compile[i]!);
  let V = C[0]!.traces.map((t) => BigInt(t));
  let Pfinal = C[0]!.mod;
  for (let i = 1; i < C.length; i++) {
    const r = multiple_crt(V, C[i]!.traces, C[i]!.mod, Pfinal);
    V = r.x;
    Pfinal = r.P;
  }
  return { V, P: Pfinal };
}

/** PARI `cost` (ellsea.c:1671-1681). */
function cost(mask: number, cost_vec: number[]): bigint {
  let c = 1n;
  for (let i = 0; i < cost_vec.length; i++) if (mask & (1 << i)) c *= BigInt(cost_vec[i]!);
  return c;
}
/** PARI `value` (ellsea.c:1683-1693). */
function value(mask: number, atkin: AtkinEntry[], k: number): bigint {
  let c = 1n;
  for (let i = 0; i < k; i++) if (mask & (1 << i)) c *= atkin[i]!.mod;
  return c;
}
/** PARI `get_lgatkin` (ellsea.c:1706-1713). */
function get_lgatkin(compile_atkin: AtkinEntry[], k: number): number[] {
  const v = new Array<number>(k);
  for (let j = 0; j < k; j++) v[j] = compile_atkin[j]!.traces.length;
  return v;
}
/** PARI `set_cost` (ellsea.c:1695-1704). */
function set_cost(B: number[], b: number, cost_vec: number[], pi: { i: number }): void {
  const costb = cost(b, cost_vec);
  let i = pi.i;
  while (i >= 1 && costb < cost(B[i]!, cost_vec)) --i;
  B[++i] = b;
  pi.i = i;
}
/** PARI `champion` (ellsea.c:1715-1758). */
function champion(
  atkin: AtkinEntry[],
  k: number,
  bound_champ: bigint
): { mask: number; cost: bigint } | null {
  const two_k = 1 << k;
  const cost_vec = get_lgatkin(atkin, k);
  if (k === 1) return { mask: 1, cost: BigInt(cost_vec[0]!) };
  const B = new Array<number>(two_k + 1).fill(0);
  const Bp = new Array<number>(two_k + 1).fill(0);
  Bp[2] = 1;
  let n = 2;
  for (let j = 2; j <= k; j++) {
    const pi = { i: 1 };
    let i1 = 2;
    let i2 = 1;
    for (; i1 <= n; ) {
      let b: number;
      const b1 = Bp[i1]!;
      const b2 = Bp[i2]! | (1 << (j - 1));
      if (value(b1, atkin, k) < value(b2, atkin, k)) {
        b = b1;
        i1++;
      } else {
        b = b2;
        i2++;
      }
      set_cost(B, b, cost_vec, pi);
    }
    for (; i2 <= n; i2++) {
      const b = Bp[i2]! | (1 << (j - 1));
      set_cost(B, b, cost_vec, pi);
    }
    n = pi.i;
    for (let i = 1; i <= n; i++) Bp[i] = B[i]!;
  }
  let res: { mask: number; cost: bigint } | null = null;
  for (let i = 1; i <= two_k; i++) {
    if (!B[i]) continue;
    const b = cost(B[i]!, cost_vec);
    const v = value(B[i]!, atkin, k);
    if (v <= bound_champ) continue;
    if (res && b >= res.cost) continue;
    res = { mask: B[i]!, cost: b };
  }
  return res;
}

/** PARI `compute_diff` (ellsea.c:1760-1767). */
function compute_diff(v: bigint[]): bigint[] {
  const s = new Set<string>();
  const out: bigint[] = [];
  for (let i = 0; i + 1 < v.length; i++) {
    const d = v[i + 1]! - v[i]!;
    const key = d.toString();
    if (!s.has(key)) {
      s.add(key);
      out.push(d);
    }
  }
  out.sort((a, b) => (a < b ? -1 : a > b ? 1 : 0));
  return out;
}
/** PARI `cmp_atkin` (ellsea.c:1769-1778). */
function cmp_atkin(a: AtkinEntry | null, b: AtkinEntry | null): number {
  const ta = a === null ? 1 : 0;
  const tb = b === null ? 1 : 0;
  if (ta || tb) return ta - tb;
  const c = a!.traces.length - b!.traces.length;
  if (c) return c;
  return b!.mod < a!.mod ? -1 : b!.mod > a!.mod ? 1 : 0;
}
/** PARI `gen_search` (bb_group.c). */
function gen_search(T: Array<AtkinEntry | null>, x: AtkinEntry): number {
  let u = T.length;
  if (!u) return -1;
  let l = 1;
  let i = 0;
  let s = 0;
  do {
    i = (l + u) >> 1;
    s = cmp_atkin(x, T[i - 1]!);
    if (!s) return i;
    if (s < 0) u = i - 1;
    else l = i + 1;
  } while (u >= l);
  return -(s < 0 ? i : i + 1);
}
/** PARI `add_atkin` (ellsea.c:1780-1790). */
function add_atkin(atkin: Array<AtkinEntry | null>, trace: AtkinEntry, nb: { n: number }): void {
  const l = atkin.length;
  let k = gen_search(atkin, trace);
  if (k > 0 || (k = -k) > l) return;
  for (let i = l; i > k; i--) atkin[i - 1] = atkin[i - 2]!;
  if (atkin[l - 1] === null) nb.n++;
  atkin[k - 1] = trace;
}

/** A hash of the x-coordinate, playing the role of PARI's `grp->hash`. */
function xhash(P: Fle): number {
  const x = P === null ? 0n : P.x;
  let h = 0;
  let v = x;
  while (v > 0n) {
    h = (h * 31 + Number(v & 0xffffn)) | 0;
    v >>= 16n;
  }
  return h;
}
function xcoord(P: Fle): bigint {
  return P === null ? 0n : P.x;
}

/** PARI `BSGS_pre` (ellsea.c:1792-1809). */
function BSGS_pre(
  V: bigint[],
  P: Fle,
  a4: bigint,
  p: bigint
): { diff: bigint[]; pre: Fle[]; index: Map<string, number> } {
  const diff = compute_diff(V);
  const index = new Map<string, number>();
  if (diff.length === 0) return { diff, pre: [], index };
  const pre = new Array<Fle>(diff.length);
  pre[0] = Fle_mul(P, diff[0]!, a4, p);
  for (let i = 1; i < diff.length; i++) {
    const d = diff[i]! - diff[i - 1]!;
    pre[i] = Fle_add(pre[i - 1]!, Fle_mul(P, d, a4, p), a4, p);
  }
  for (let i = 0; i < diff.length; i++) index.set(diff[i]!.toString(), i);
  return { diff, pre, index };
}

/**
 * PARI `match_and_sort` (ellsea.c:1815-1958): combine the Atkin information
 * with the trace `u` known modulo `Mu` by a baby-step/giant-step search.
 */
function match_and_sort(
  compile_atkin: AtkinEntry[],
  Mu: bigint,
  u: bigint,
  q: bigint,
  a4: bigint,
  a6: bigint
): bigint {
  const pp1 = q + 1n;
  const k = compile_atkin.length;
  const bound = sqrti(q << 2n);
  const p = q; /* prime field */
  if (k === 1) {
    /* only one Atkin prime, check the cardinality with random points */
    const r = compile_atkin[0]!;
    const C = Mu * r.mod;
    const Cs2 = C >> 1n;
    const card: bigint[] = [];
    for (const t2 of r.traces) {
      let t = crt2(u, Mu, BigInt(t2), r.mod);
      if (t > Cs2) t -= C;
      const at = t < 0n ? -t : t;
      if (at <= bound) card.push(pp1 - t);
    }
    return gen_select_order(card, a4, a6, p);
  }
  const best_i = separation(get_lgatkin(compile_atkin, k));
  const babyR = possible_traces(compile_atkin, best_i);
  const giantR = possible_traces(compile_atkin, (1 << k) - 1 - best_i);
  const baby = babyR.V;
  const giant = giantR.V;
  const Mb = babyR.P;
  const Mg = giantR.P;
  let den = Fp_inv(mod(Mu * Mb, Mg), Mg);
  for (let i = 0; i < giant.length; i++) giant[i] = mod(giant[i]! * den, Mg);
  giant.sort((x, y) => (x < y ? -1 : x > y ? 1 : 0));
  const Sg = mod(-u * den, Mg);
  den = Fp_inv(mod(Mu * Mg, Mb), Mb);
  let dec_inf = (Mb * (Mg + 2n * Sg)) / (2n * Mg);
  dec_inf = -dec_inf; /* ceil(-(Mb/2 + Sg Mb/Mg)) */
  const div = truedivii(dec_inf, Mb) * Mb;
  for (let i = 0; i < baby.length; i++) {
    let b = mod((baby[i]! - u) * den, Mb) + div;
    if (b < dec_inf) b += Mb;
    baby[i] = b;
  }
  baby.sort((x, y) => (x < y ? -1 : x > y ? 1 : 0));
  const SgMb = Sg * Mb;
  const lbaby = baby.length;
  const lgiant = giant.length;

  for (;;) {
    /* MATCH_RESTART */
    const card: bigint[] = [];
    const lcard = 100;
    let restart = false;
    const P = random_Fle(a4, a6, p);
    let point = Fle_mul(P, Mu, a4, p);
    const Pb = Fle_mul(point, Mg, a4, p);
    const Pg = Fle_mul(point, Mb, a4, p);
    /* Precomputation for babies */
    let pre = BSGS_pre(baby, Pb, a4, p);
    const table = new Array<number>(lbaby);
    point = Fle_mul(P, pp1 - u - Mu * (SgMb + Mg * baby[0]!), a4, p);
    table[0] = xhash(point);
    for (let i = 1; i < lbaby; i++) {
      const d = baby[i]! - baby[i - 1]!;
      const idx = pre.index.get(d.toString())!;
      point = Fle_add(point, Fle_neg(pre.pre[idx]!, p), a4, p);
      table[i] = xhash(point);
    }
    /* Precomputations for giants */
    pre = BSGS_pre(giant, Pg, a4, p);
    /* Look for a collision among the x-coordinates */
    const table_ind = table.map((_, i) => i).sort((i, j) => table[i]! - table[j]!);
    const tables = table_ind.map((i) => table[i]!);

    point = Fle_mul(Pg, giant[0]!, a4, p);
    for (let i = 0; ; i++) {
      const h = xhash(point);
      let s = zv_search(tables, h);
      if (s) {
        while (s > 1 && tables[s - 2] === h) s--;
        for (; s <= lbaby && tables[s - 1] === h; s++) {
          const B = baby[table_ind[s - 1]!]!;
          const G = giant[i]!;
          const GMb = G * Mb;
          const BMg = B * Mg;
          const Be = pp1 - u - Mu * (SgMb + BMg);
          const Bp = Fle_mul(P, Be, a4, p);
          if (xcoord(Bp) === xcoord(point)) {
            const card1 = Be - Mu * GMb;
            const card2 = card1 + 2n * Mu * GMb;
            const d1 = pp1 - card1;
            if ((d1 < 0n ? -d1 : d1) <= bound) card.push(card1);
            if (card.length >= lcard) {
              restart = true;
              break;
            }
            const d2 = pp1 - card2;
            if ((d2 < 0n ? -d2 : d2) <= bound) card.push(card2);
            if (card.length >= lcard) {
              restart = true;
              break;
            }
          }
        }
        if (restart) break;
      }
      if (i === lgiant - 1) break;
      const d = giant[i + 1]! - giant[i]!;
      const idx = pre.index.get(d.toString())!;
      point = Fle_add(point, pre.pre[idx]!, a4, p);
    }
    if (restart) continue;
    if (card.length === 0) throw new PariBugError('match_and_sort');
    return gen_select_order(card, a4, a6, p);
  }
}
/** PARI `zv_search`: 1-based index of `x` in the sorted array, or 0. */
function zv_search(T: number[], x: number): number {
  let l = 1;
  let u = T.length;
  while (u >= l) {
    const i = (l + u) >> 1;
    if (x < T[i - 1]!) u = i - 1;
    else if (x > T[i - 1]!) l = i + 1;
    else return i;
  }
  return 0;
}
/** CRT of `a mod A` and `b mod B` with `A`, `B` coprime, result in `[0, AB)`. */
function crt2(a: bigint, A: bigint, b: bigint, B: bigint): bigint {
  const { u, v } = bezout(A, B);
  const AB = A * B;
  return mod(a * (B * v) + b * (A * u), AB);
}

/** PARI `get_bound_bsgs` (ellsea.c:1960-1971).  Upstream is inexact here. */
function get_bound_bsgs(lp: number): number {
  let B: number;
  if (lp <= 160) B = Math.pow(1.048, lp) / 9;
  else if (lp <= 192) B = Math.pow(1.052, lp) / 16.65;
  else B = Math.pow(1.035, Math.min(lp, 307)) * 1.35;
  return B * 1000000;
}

/* ================================================================== */
/* The main entry point (ellsea.c:1973-2113)                          */
/* ================================================================== */

/** PARI `Z_incremental_CRT` (polarit3.c:1193) + `Fl_chinese_coprime`. */
function Z_incremental_CRT(
  H: bigint,
  Hp: number,
  q: bigint,
  p: number
): { H: bigint; q: bigint } {
  const P = BigInt(p);
  const qp = q * P;
  const qinv = Fp_inv(mod(q, P), P);
  const amod = mod(H, P);
  if (BigInt(Hp) === amod) return { H, q: qp };
  const d = mod((BigInt(Hp) - amod) * qinv, P);
  let ax: bigint;
  if (d >= 1n + (P >> 1n)) ax = H - (P - d) * q;
  else {
    ax = H + d * q;
    if (ax > qp >> 1n) ax -= qp;
  }
  return { H: ax, q: qp };
}

/** Odd primes, in increasing order (PARI `u_forprime_init(&TT, 3, ...)`). */
function* oddPrimes(): Generator<number> {
  const primes: number[] = [];
  for (let n = 3; ; n += 2) {
    let isP = true;
    for (const q of primes) {
      if (q * q > n) break;
      if (n % q === 0) {
        isP = false;
        break;
      }
    }
    if (isP) {
      primes.push(n);
      yield n;
    }
  }
}

/**
 * PARI `Fq_ellcard_SEA` (ellsea.c:1978-2109) with `T = NULL`, i.e.
 * `Fp_ellcard_SEA` (ellsea.c:2111-2113).
 *
 * @param a4 - curve coefficient, `0 <= a4 < p`
 * @param a6 - curve coefficient, `0 <= a6 < p`
 * @param p  - a prime `> 3`
 * @param smallfact - as upstream: stop and return `0` as soon as a prime
 *   factor of `#E(Fp)` not dividing `smallfact` is found (negative: also test
 *   the quadratic twist)
 */
export function Fp_ellcard_SEA(a4: bigint, a6: bigint, p: bigint, smallfact = 0): bigint {
  const MAX_ATKIN = 21;
  const q = p;
  const j = Fq_ellj(a4, a6, p);
  if (j === 0n || Fq_sub(j, 1728n, p) === 0n) {
    const c = Fp_ellcard_CM(a4, a6, p);
    if (c === null) throw new PariBugError('Fp_ellcard_SEA (j = 0 or 1728)');
    return c;
  }
  if (Fp_elljissupersingular(j, p)) return p + 1n;

  /* First compute the trace modulo 2 */
  let TR: bigint;
  let TR_mod: bigint;
  switch (FpX_nbroots(rhs(a4, a6, p), p)) {
    case 3: {
      /* bonus time: 4 | #E(Fq) = q+1 - t */
      let i = Number(mod(q, 4n)) + 1;
      if (i > 2) i -= 4;
      TR_mod = 4n;
      TR = BigInt(i);
      break;
    }
    case 1:
      TR_mod = 2n;
      TR = 0n;
      break;
    default:
      TR_mod = 2n;
      TR = 1n;
      break;
  }
  if (smallfact % 2 !== 0 && mod(TR, 2n) === 0n) return 0n;

  const bound = sqrti(q << 4n);
  const bound_bsgs = get_bound_bsgs(expi(q));
  const compile_atkin: Array<AtkinEntry | null> = new Array(MAX_ATKIN).fill(null);
  const nb = { n: 0 };
  let prod_atkin = 1n;
  let max_traces = 0n;
  let bound_gr = 1;
  const growth_factor = 1.26;

  for (const ell of oddPrimes()) {
    if (BigInt(ell) === p) continue;
    const ptkt = { kt: 1 };
    const res = find_trace(a4, a6, j, ell, q, p, ptkt, smallfact);
    if (!res) continue;
    if ('card' in res) return res.card;
    const trace_mod = res.tr;
    const kt = ptkt.kt;
    const nbtrace = trace_mod.length;
    const ellkt = Math.pow(ell, kt);
    if (nbtrace === 1) {
      const t_mod_ellkt = trace_mod[0]!;
      if (smallfact && smallfact % ell !== 0) {
        /* does ell divide q + 1 - t ? */
        const q_mod_ell_plus_one = Number(mod(q, BigInt(ell))) + 1;
        const card_mod_ell = ((q_mod_ell_plus_one - t_mod_ellkt) % ell + ell) % ell;
        let tcard_mod_ell = 1;
        if (card_mod_ell && smallfact < 0)
          tcard_mod_ell = ((q_mod_ell_plus_one + t_mod_ellkt) % ell + ell) % ell;
        if (!card_mod_ell || !tcard_mod_ell) return 0n;
      }
      const upd = Z_incremental_CRT(TR, t_mod_ellkt % ellkt, TR_mod, ellkt);
      TR = upd.H;
      TR_mod = upd.q;
    } else {
      add_atkin(compile_atkin, { mod: BigInt(ellkt), traces: trace_mod }, nb);
      prod_atkin = value(-1, compile_atkin.slice(0, nb.n) as AtkinEntry[], nb.n);
    }
    err_printf(`  missing ${expi(bound) - expi(TR_mod * prod_atkin)} bits\n`);
    if (TR_mod * prod_atkin > bound) {
      if (!nb.n) return q + 1n - TR;
      const bound_tr = bound_bsgs * bound_gr;
      bound_gr *= growth_factor;
      if (max_traces !== 0n) {
        max_traces = (max_traces * BigInt(nbtrace)) / BigInt(ellkt);
      }
      if (Number(max_traces) < bound_tr) {
        const bound_atkin = truedivii(bound, TR_mod);
        const champ = champion(compile_atkin.slice(0, nb.n) as AtkinEntry[], nb.n, bound_atkin);
        if (champ) {
          max_traces = champ.cost;
          if (Number(max_traces) < bound_tr) {
            const cat: AtkinEntry[] = [];
            for (let i = 0; i < nb.n; i++)
              if (champ.mask & (1 << i)) cat.push(compile_atkin[i]!);
            err_printf(`Match and sort for ${max_traces} possibilities.\n`);
            return match_and_sort(cat, TR_mod, TR, q, a4, a6);
          }
        }
      }
    }
  }
  /* not reached */
  throw new PariBugError('Fp_ellcard_SEA');
}

/** Test-only surface. */
export const _internal = {
  FpX_red,
  FpX_mul,
  FpX_rem,
  FpX_gcd,
  FpX_roots,
  FpX_ddf_degree,
  FpX_Frobenius,
  FpXn_inv,
  FpXn_expint,
  FpX_recip,
  Fq_polmodular_eval,
  get_modular_eqn,
  meqn_j,
  find_isogenous,
  find_kernel,
  Fq_elldivpolmod,
  find_trace_Atkin,
  find_trace_one_root,
  find_trace_lp1_roots,
  study_modular_eqn,
  separation,
  champion,
  possible_traces,
  Z_incremental_CRT,
  list_singular_discs,
  Fp_jissupersingular,
  Fp2_sqrt,
  sqrti,
  expi,
  ModType,
};

