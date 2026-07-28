/**
 * @module parigp-ts/qfb
 * @description Binary quadratic forms: composition, reduction and representation.
 *
 * Direct port of `reference/pari/src/basemath/Qfb.c`:
 *
 * - `Qfb.c:44-61`      `quadpoly_bc` (used by the principal form)
 * - `Qfb.c:148-176`    `Qfb0` (constructor + checks)
 * - `Qfb.c:187-271`    `dvmdii_round`, `REDB`, `REDBU`
 * - `Qfb.c:273-311`    `qfi_redsl2_basecase`
 * - `Qfb.c:376-395`    `qfi_rho`
 * - `Qfb.c:434-450`    `rho_get_BC`, `qfr3_rho`
 * - `Qfb.c:504-533`    `ab_isreduced`, `qfr3_red`
 * - `Qfb.c:594-654`    `qfr_red_basecase_i`, `qfr_rhosl2`, `qfr_redsl2_basecase`
 * - `Qfb.c:656-825`    Schoenhage fast reduction `pqfbred_1/is_minimal/pqfbred_rec`
 * - `Qfb.c:826-1006`   `qfr_redsl2`, `qfi_redsl2`, `qfbredsl2`, `qfr_red`, `qfi_red`, `qfbred0`
 * - `Qfb.c:1013-1071`  `qfb_sqr`, `qfb_comp` (Gauss/Dirichlet composition)
 * - `Qfb.c:1074-1180`  `qfb_comp_gen`, `qficomp0`, `qfrcomp0`, `qfbcomp`, `qfbcompraw`
 * - `Qfb.c:1296-1325`  `qfipowraw`, `qfipow` (and the qfr analogues at 1608-1632)
 * - `Qfb.c:1661-1749`  `primeform_u`, `primeform`
 * - `Qfb.c:1751-1783`  `normforms`
 * - `Qfb.c:1785-1895`  `qfisolvep_all`, `qfrsolve_normform`
 * - `Qfb.c:1916-1995`  `qfbsolve_primitive`, `qfbsolve_all`, `qfbsolve`
 * - `Qfb.c:1997-2017`  `cornacchia`, `Qfb.c:2019-2085` `cornacchia2`
 * - `quad.c:1149-1260` `Zn_quad_roots`
 *
 * Everything here is exact integer arithmetic. PARI additionally carries a
 * floating-point "Shanks distance" for *indefinite* forms, represented as a
 * `t_VEC [qfb, t_REAL]` (`qfr5_*`, `qfbred`/`qfbpow`/`qfbcomp` on such a pair).
 * That family is NOT ported: the functions here accept a plain {@link Qfb} only,
 * which is exactly PARI's `t_QFB` branch, and every PARI code path we reproduce
 * for it (`qfr3_*`) is integer-only. If you need Shanks distances, they are
 * still missing.
 */

import { Fp_sqrt, kronecker } from './ff.js';
import { isPrime, NotImplementedError, Z_factor, type Factorization } from './ifactor.js';
import {
  PariDomainError,
  PariFlagError,
  PariInvError,
  PariSqrtnError,
  PariTypeError,
} from './matkermod.js';

/* ------------------------------------------------------------------ */
/* Errors                                                              */
/* ------------------------------------------------------------------ */

/*
 * PARI error kinds are shared across the parigp-ts modules; they are currently
 * defined in `matkermod.ts` (see the note there) and re-exported here.
 */
export {
  PariDomainError,
  PariFlagError,
  PariSqrtnError,
  PariTypeError,
} from './matkermod.js';
export { NotImplementedError } from './ifactor.js';

/* ------------------------------------------------------------------ */
/* Integer helpers mirroring the PARI kernel                           */
/* ------------------------------------------------------------------ */

const iabs = (x: bigint): bigint => (x < 0n ? -x : x);
const isign = (x: bigint): number => (x > 0n ? 1 : x < 0n ? -1 : 0);

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

/** `[d,u,v]` with `d = gcd(x,y) >= 0` and `u*x+v*y = d` (PARI `bezout`) */
function bezout(x: bigint, y: bigint): [bigint, bigint, bigint] {
  let [r0, r1] = [x, y];
  let [s0, s1] = [1n, 0n];
  let [t0, t1] = [0n, 1n];
  while (r1 !== 0n) {
    const q = r0 / r1; // truncated; any consistent choice works
    [r0, r1] = [r1, r0 - q * r1];
    [s0, s1] = [s1, s0 - q * s1];
    [t0, t1] = [t1, t0 - q * t1];
  }
  if (r0 < 0n) return [-r0, -s0, -t0];
  return [r0, s0, t0];
}

/** PARI `dvmdii`: truncated division, remainder has the sign of `x`. */
function dvmdii(x: bigint, y: bigint): [bigint, bigint] {
  const q = x / y;
  return [q, x - q * y];
}

/** PARI `truedvmdii`: `x = q*y + r`, `0 <= r < |y|`. */
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

const truedivii = (x: bigint, y: bigint): bigint => truedvmdii(x, y)[0];

/** PARI `shifti(x,n)`: sign-magnitude shift (truncation toward zero for n<0). */
function shifti(x: bigint, n: number): bigint {
  if (n >= 0) return x << BigInt(n);
  const s = x < 0n;
  const v = (s ? -x : x) >> BigInt(-n);
  return s ? -v : v;
}

/** PARI `remi2n(x,n)`: `x` mod `2^n`, keeping the sign of `x`. */
function remi2n(x: bigint, n: number): bigint {
  const s = x < 0n;
  const v = (s ? -x : x) & ((1n << BigInt(n)) - 1n);
  return s ? -v : v;
}

/** PARI `expi(x)` = floor(log2 |x|); we use -1 for 0 (see the notes in Qfb.c usage). */
function expi(x: bigint): number {
  const v = iabs(x);
  return v === 0n ? -1 : v.toString(2).length - 1;
}

/*
 * PARI's `mod2`/`mod4`/`mod8`/`mod16` operate on the *magnitude* of a t_INT
 * (PARI stores integers in sign + magnitude form), e.g. `mod8(-23) == 7`, not 1.
 * Several algorithms here (primeform, Z2_sqrt) depend on that convention.
 */
const modk = (x: bigint, k: bigint): number => Number(iabs(x) % k);
const mod2 = (x: bigint): number => modk(x, 2n);
const mod4 = (x: bigint): number => modk(x, 4n);
const mod8 = (x: bigint): number => modk(x, 8n);
const mod16 = (x: bigint): number => modk(x, 16n);
/** PARI's capital `Mod4`/`Mod8`: the true residue, handling the sign. */
const Mod4 = (x: bigint): number => Number(((x % 4n) + 4n) % 4n);
const Mod8 = (x: bigint): number => Number(((x % 8n) + 8n) % 8n);

/** integer square root (floor) of `n >= 0`; PARI `sqrti` */
export function sqrti(n: bigint): bigint {
  if (n < 0n) throw new PariDomainError('sqrti', 'argument', '<', '0');
  if (n < 2n) return n;
  let x = 1n << BigInt(((n.toString(2).length + 1) >> 1) + 1);
  for (;;) {
    const y = (x + n / x) >> 1n;
    if (y >= x) break;
    x = y;
  }
  return x;
}

/** PARI `Z_issquareall`: `[true, sqrt]` if `n` is a perfect square */
export function Z_issquareall(n: bigint): [boolean, bigint] {
  if (n < 0n) return [false, 0n];
  const r = sqrti(n);
  return r * r === n ? [true, r] : [false, 0n];
}

/** PARI `modii(x,y)`: the representative of `x` mod `y` in `[0, |y|)`. */
function modii(x: bigint, y: bigint): bigint {
  const m = iabs(y);
  const r = x % m;
  return r < 0n ? r + m : r;
}

function Fp_red(a: bigint, p: bigint): bigint {
  const r = a % p;
  return r < 0n ? r + p : r;
}

function Fp_inv(a: bigint, p: bigint): bigint {
  const [g, u] = bezout(Fp_red(a, p), p);
  if (g !== 1n) throw new PariInvError(`Fp_inv (mod ${p})`);
  return Fp_red(u, p);
}

function Fp_pow(a: bigint, e: bigint, p: bigint): bigint {
  let r = 1n;
  let b = Fp_red(a, p);
  let k = e;
  while (k > 0n) {
    if (k & 1n) r = (r * b) % p;
    b = (b * b) % p;
    k >>= 1n;
  }
  return r;
}

/* ------------------------------------------------------------------ */
/* Type                                                                */
/* ------------------------------------------------------------------ */

/** A binary quadratic form `a x^2 + b x y + c y^2`, with cached discriminant. */
export interface Qfb {
  readonly a: bigint;
  readonly b: bigint;
  readonly c: bigint;
  /** `b^2 - 4ac` */
  readonly D: bigint;
}

export function mkqfb(a: bigint, b: bigint, c: bigint, D: bigint): Qfb {
  return { a, b, c, D };
}

export function qfb_disc3(a: bigint, b: bigint, c: bigint): bigint {
  return b * b - 4n * a * c;
}

export function qfb_disc(q: Qfb): bigint {
  return q.D;
}

/** definite (D < 0) forms are PARI's "qfi", indefinite ones its "qfr" */
export function qfb_is_qfi(q: Qfb): boolean {
  return q.D < 0n;
}

/** PARI `Qfb0` (`Qfb.c:148-176`): build a form with PARI's validity checks. */
export function Qfb(a: bigint, b: bigint, c: bigint): Qfb {
  const D = qfb_disc3(a, b, c);
  if (D < 0n) {
    if (a < 0n) throw new NotImplementedError('negative definite t_QFB');
  } else if (Z_issquareall(D)[0]) {
    throw new PariDomainError('Qfb', 'issquare(disc)', '=', '1');
  }
  return mkqfb(a, b, c, D);
}

export function qfb_equal(x: Qfb, y: Qfb): boolean {
  return x.a === y.a && x.b === y.b && x.c === y.c;
}

/** PARI `qfbinv` (`Qfb.c:145-147`) */
export function qfbinv(x: Qfb): Qfb {
  return mkqfb(x.a, -x.b, x.c, x.D);
}

/** PARI `qfi_1_by_disc` / `qfr_1_by_disc` (`Qfb.c:1213-1234`, `1268-1284`) */
export function qfb_1(x: Qfb): Qfb {
  const D = x.D;
  if (D < 0n) {
    /* quadpoly_bc(D, D odd): b = -1, c = (1-D)/4; qfi_1_by_disc flips b to 1 */
    const b = mod2(D) !== 0 ? 1n : 0n;
    const c = mod2(D) !== 0 ? (1n - D) / 4n : -(D / 4n);
    return mkqfb(1n, b, c, D);
  }
  /* qfr_1_by_disc: s = floor(sqrt(D)), decremented so that s = D mod 2 */
  let s = sqrti(D);
  if (mod2(s) !== mod2(D)) s = s - 1n;
  return mkqfb(1n, s, (s * s - D) / 4n, D);
}

/* ------------------------------------------------------------------ */
/* Composition (Qfb.c:1013-1071)                                       */
/* ------------------------------------------------------------------ */

/** PARI `qfb_sqr` (`Qfb.c:1013-1036`) */
function qfb_sqr(x: Qfb): Qfb {
  const [d1, x2] = bezout(x.b, x.a); /* usually d1 = 1 */
  let c = x.c;
  let m = c * x2;
  let v1: bigint, v2: bigint;
  if (d1 === 1n) {
    v1 = v2 = x.a;
  } else {
    v1 = x.a / d1;
    v2 = v1 * gcdii(d1, c);
    c = c * d1;
  }
  m = -m;
  const r = modii(m, v2);
  const p1 = r * v1;
  const c3 = c + r * (x.b + p1);
  return mkqfb(v1 * v2, x.b + 2n * p1, c3 / v2, x.D);
}

/** PARI `qfb_comp` (`Qfb.c:1038-1071`) */
function qfb_comp(x: Qfb, y: Qfb): Qfb {
  /* PARI checks pointer identity only (`Qfb.c:1043`); structurally equal but
   * distinct forms go through the general composition. */
  if (x === y) return qfb_sqr(x);
  const n = shifti(y.b - x.b, -1);
  let v1 = x.a;
  let v2 = y.a;
  let c = y.c;
  const [d, y1] = bezout(v2, v1);
  let m: bigint;
  if (d === 1n) {
    m = y1 * n;
  } else {
    const s = y.b - n;
    const [d1, x2, y2] = bezout(s, d); /* x2*s + y2*d = d1 */
    if (d1 !== 1n) {
      v1 = v1 / d1;
      v2 = v2 / d1;
      v1 = v1 * gcdii(c, gcdii(x.c, gcdii(d1, n)));
      c = c * d1;
    }
    m = y1 * y2 * n + y.c * x2;
  }
  m = -m;
  const r = modii(m, v1);
  const p1 = r * v2;
  const c3 = c + r * (y.b + p1);
  return mkqfb(v1 * v2, y.b + 2n * p1, c3 / v1, x.D);
}

function content3(a: bigint, b: bigint, c: bigint): bigint {
  return gcdii(a, gcdii(b, c));
}

/** extended gcd of three integers: `[g, u1, u2, u3]` with `sum ui*vi = g` */
function ZV3_extgcd(A: bigint, B: bigint, C: bigint): [bigint, bigint, bigint, bigint] {
  const [g1, u1, u2] = bezout(A, B);
  const [g, w1, w2] = bezout(g1, C);
  return [g, w1 * u1, w1 * u2, w2];
}

/** PARI `qfb_comp_gen` (`Qfb.c:1074-1116`): composition of forms of different discriminants */
function qfb_comp_gen(x: Qfb, y: Qfb): Qfb | null {
  let d1 = x.D;
  let d2 = y.D;
  let a1 = x.a,
    b1 = x.b,
    c1 = x.c;
  let a2 = y.a,
    b2 = y.b,
    c2 = y.c;
  let cx = content3(x.a, x.b, x.c);
  const cy = content3(y.a, y.b, y.c);
  if (cx !== 1n && cx !== -1n) {
    a1 /= cx;
    b1 /= cx;
    c1 /= cx;
    d1 /= cx * cx;
  }
  if (cy !== 1n && cy !== -1n) {
    a2 /= cy;
    b2 /= cy;
    c2 /= cy;
    d2 /= cy * cy;
  }
  let D = gcdii(d1, d2);
  if (d1 < 0n) D = -D;
  const [ok1, n1] = Z_issquareall(d1 / D);
  const [ok2, n2] = Z_issquareall(d2 / D);
  if (!ok1 || !ok2) return null;
  const A0 = a1 * n2;
  const B0 = a2 * n1;
  const C0 = shifti(b1 * n2 + b2 * n1, -1);
  const [m, U1, U2, U3] = ZV3_extgcd(A0, B0, C0);
  let A = ((a1 * b2) / m) * U1;
  let B = ((a2 * b1) / m) * U2;
  let C = b1 * b2 + D * n1 * n2;
  C = (shifti(C, -1) / m) * U3;
  B = A + B + C;
  const m2 = m * m;
  A = (a1 * a2) / m2;
  C = shifti(B * B - D, -2) / A;
  cx = cx * cy;
  if (cx !== 1n && cx !== -1n) {
    A = A * cx;
    B = B * cx;
    C = C * cx;
    D = D * cx * cx;
  }
  return mkqfb(A, B, C, D);
}

/* ------------------------------------------------------------------ */
/* Reduction: definite forms (Qfb.c:187-395, 947-989)                  */
/* ------------------------------------------------------------------ */

/** PARI `dvmdii_round` (`Qfb.c:188-198`): `b = q*2a + r`, `-a < r <= a`, `a > 0`. */
function dvmdii_round(b: bigint, a: bigint): [bigint, bigint] {
  const a2 = shifti(a, 1);
  let [q, r] = dvmdii(b, a2);
  if (b >= 0n) {
    if (iabs(r) > a) {
      q = q + 1n;
      r = r - a2;
    }
  } else {
    if (iabs(r) >= a) {
      q = q - 1n;
      r = r + a2;
    }
  }
  return [q, r];
}

/** PARI `REDB` (`Qfb.c:220-227`): reduce `b` mod `2a`, updating `b` and `c`. */
function REDB(a: bigint, b: bigint, c: bigint): [bigint, bigint] {
  const [q, r] = dvmdii_round(b, a);
  if (q === 0n) return [b, c];
  return [r, c - q * shifti(b + r, -1)];
}

/** PARI `REDBU` (`Qfb.c:263-271`) */
function REDBU(
  a: bigint,
  b: bigint,
  c: bigint,
  u1: bigint,
  u2: bigint
): [bigint, bigint, bigint] {
  const [q, r] = dvmdii_round(b, a);
  return [r, c - q * shifti(b + r, -1), u2 - q * u1];
}

/** PARI `qfi_red_basecase_av` (`Qfb.c:947-974`) */
function qfi_red_basecase(q: Qfb): Qfb {
  let a = q.a,
    b = q.b,
    c = q.c;
  let cmp = iabs(a) - iabs(b);
  if (cmp < 0n) [b, c] = REDB(a, b, c);
  else if (cmp === 0n && b < 0n) b = -b;
  for (;;) {
    cmp = iabs(a) - iabs(c);
    if (cmp <= 0n) break;
    [a, c] = [c, a];
    b = -b;
    [b, c] = REDB(a, b, c);
  }
  if (cmp === 0n && b < 0n) b = -b;
  return mkqfb(a, b, c, q.D);
}

/** PARI `qfi_rho` (`Qfb.c:376-395`) */
export function qfi_rho(x: Qfb): Qfb {
  let a = x.a,
    b = x.b,
    c = x.c;
  const fl = iabs(a) - iabs(c);
  if (fl <= 0n) {
    const fg = iabs(a) - iabs(b);
    if (fg >= 0n) {
      let bb = b;
      if ((fl === 0n || fg === 0n) && bb < 0n) bb = -bb;
      return mkqfb(a, bb, c, x.D);
    }
  }
  [a, c] = [c, a];
  b = -b;
  [b, c] = REDB(a, b, c);
  return mkqfb(a, b, c, x.D);
}

/** PARI `qfi_redsl2_basecase` (`Qfb.c:273-311`) */
function qfi_redsl2_basecase(q: Qfb): { Q: Qfb; U: bigint[][] } {
  let a = q.a,
    b = q.b,
    c = q.c;
  let u1 = 1n,
    u2 = 0n;
  let cmp = iabs(a) - iabs(b);
  if (cmp < 0n) {
    [b, c, u2] = REDBU(a, b, c, u1, u2);
  } else if (cmp === 0n && b < 0n) {
    b = -b;
    u2 = 1n;
  }
  for (;;) {
    cmp = iabs(a) - iabs(c);
    if (cmp <= 0n) break;
    [a, c] = [c, a];
    b = -b;
    const z = u1;
    u1 = u2;
    u2 = -z;
    [b, c, u2] = REDBU(a, b, c, u1, u2);
  }
  if (cmp === 0n && b < 0n) {
    b = -b;
    const z = u1;
    u1 = u2;
    u2 = -z;
  }
  let z = shifti(b - q.b, -1);
  let v1 = z * u1 - a * u2;
  v1 = v1 / q.c;
  z = z - b;
  let v2 = z * u2 + c * u1;
  v2 = v2 / q.c;
  return {
    Q: mkqfb(a, b, c, q.D),
    U: [
      [u1, u2],
      [v1, v2],
    ],
  };
}

/* ------------------------------------------------------------------ */
/* Reduction: indefinite forms (Qfb.c:434-533, 613-654)                */
/* ------------------------------------------------------------------ */

/** PARI `rho_get_BC` (`Qfb.c:434-442`) */
function rho_get_BC(
  a: bigint,
  b: bigint,
  c: bigint,
  isqrtD: bigint
): [bigint, bigint] {
  const t = iabs(isqrtD) >= iabs(c) ? isqrtD : iabs(c);
  const [q, u] = truedvmdii(t + b, shifti(c, 1));
  return [t - u, a - q * (b - q * c)];
}

/** PARI `ab_isreduced` (`Qfb.c:507-516`) */
function ab_isreduced(a: bigint, b: bigint, isqrtD: bigint): boolean {
  if (b <= 0n || iabs(b) > iabs(isqrtD)) return false;
  const t = isqrtD - iabs(shifti(a, 1));
  return t < 0n ? iabs(b) >= iabs(t) : iabs(b) > iabs(t);
}

/** PARI `qfr3_red` (`Qfb.c:518-533`) */
function qfr3_red(a: bigint, b: bigint, c: bigint, isqrtD: bigint): [bigint, bigint, bigint] {
  while (!ab_isreduced(a, b, isqrtD)) {
    const [B, C] = rho_get_BC(a, b, c, isqrtD);
    a = c;
    b = B;
    c = C;
  }
  return [a, b, c];
}

/** PARI `qfr3_rho` (`Qfb.c:444-450`) */
function qfr3_rho(a: bigint, b: bigint, c: bigint, isqrtD: bigint): [bigint, bigint, bigint] {
  const [B, C] = rho_get_BC(a, b, c, isqrtD);
  return [c, B, C];
}

/** PARI `qfr_rhosl2_i` (`Qfb.c:613-626`) */
function qfr_rhosl2_i(
  a: bigint,
  b: bigint,
  c: bigint,
  u1: bigint,
  u2: bigint,
  v1: bigint,
  v2: bigint,
  rd: bigint
): [bigint, bigint, bigint, bigint, bigint, bigint, bigint] {
  const C = iabs(c);
  const t = b + (rd > C ? rd : C);
  let [q, r] = truedvmdii(t, shifti(C, 1));
  if (c < 0n) q = -q;
  const na = c;
  const nb = t - (r + b);
  const nc = a - q * (b - q * c);
  const r1 = u1;
  const nu1 = v1;
  const nv1 = q * v1 - r1;
  const r2 = u2;
  const nu2 = v2;
  const nv2 = q * v2 - r2;
  return [na, nb, nc, nu1, nu2, nv1, nv2];
}

/** PARI `qfr_redsl2_basecase` (`Qfb.c:640-655`) */
function qfr_redsl2_basecase(V: Qfb, rd: bigint): { Q: Qfb; U: bigint[][] } {
  let u1 = 1n,
    u2 = 0n,
    v1 = 0n,
    v2 = 1n;
  let a = V.a,
    b = V.b,
    c = V.c;
  const d = V.D;
  while (!ab_isreduced(a, b, rd)) {
    [a, b, c, u1, u2, v1, v2] = qfr_rhosl2_i(a, b, c, u1, u2, v1, v2, rd);
  }
  return {
    Q: mkqfb(a, b, c, d),
    U: [
      [u1, v1],
      [u2, v2],
    ],
  };
}

/** PARI `qfr_rhosl2` (`Qfb.c:628-638`) */
function qfr_rhosl2(A: { Q: Qfb; U: bigint[][] }, rd: bigint): { Q: Qfb; U: bigint[][] } {
  const V = A.Q;
  const M = A.U;
  const [a, b, c, u1, u2, v1, v2] = qfr_rhosl2_i(
    V.a,
    V.b,
    V.c,
    M[0]![0]!,
    M[1]![0]!,
    M[0]![1]!,
    M[1]![1]!,
    rd
  );
  return {
    Q: mkqfb(a, b, c, V.D),
    U: [
      [u1, v1],
      [u2, v2],
    ],
  };
}

/* ------------------------------------------------------------------ */
/* Schoenhage fast reduction (Qfb.c:656-825)                           */
/* ------------------------------------------------------------------ */

/** PARI `lti2n` (`Qfb.c:668-669`) */
function lti2n(a: bigint, m: number): boolean {
  return a < 0n || expi(a) < m;
}

/** PARI `pqfbred_1` (`Qfb.c:671-704`) */
function pqfbred_1(Q: Qfb, m: number, U: bigint[][]): Qfb {
  let a = Q.a,
    b = Q.b,
    c = Q.c;
  const d = Q.D;
  if (iabs(a) < iabs(c)) {
    const r2 = shifti(a, m + 2) + d;
    const e2 = expi(r2);
    const r = 1n << BigInt(r2 < 0n || e2 < 2 * m + 2 ? m + 1 : e2 >> 1);
    const t = truedivii(b - r, shifti(a, 1));
    if (t === 0n) throw new Error('pqfbred_1: bug (t = 0)');
    const at = a * t;
    c = c - b * t + at * t;
    b = b - shifti(at, 1);
    U[0]![1] = U[0]![1]! - U[0]![0]! * t;
    U[1]![1] = U[1]![1]! - U[1]![0]! * t;
  } else {
    const r2 = shifti(c, m + 2) + d;
    const e2 = expi(r2);
    const r = 1n << BigInt(r2 < 0n || e2 < 2 * m + 2 ? m + 1 : e2 >> 1);
    const t = truedivii(b - r, shifti(c, 1));
    if (t === 0n) throw new Error('pqfbred_1: bug (t = 0)');
    const ct = c * t;
    a = a - b * t + ct * t;
    b = b - shifti(ct, 1);
    U[0]![0] = U[0]![0]! - U[0]![1]! * t;
    U[1]![0] = U[1]![0]! - U[1]![1]! * t;
  }
  return mkqfb(a, b, c, d);
}

/** PARI `is_minimal` (`Qfb.c:706-713`) */
function is_minimal(Q: Qfb, m: number): boolean {
  return (
    lti2n(Q.a - Q.b + Q.c, m) ||
    (lti2n(Q.b - shifti(Q.a, 1), m + 1) && lti2n(Q.b - shifti(Q.c, 1), m + 1))
  );
}

/** PARI `pqfbred_iter_1` (`Qfb.c:715-729`) */
function pqfbred_iter_1(Q: Qfb, m: number, U: bigint[][]): Qfb {
  while (!is_minimal(Q, m)) Q = pqfbred_1(Q, m, U);
  return Q;
}

function matid2(): bigint[][] {
  return [
    [1n, 0n],
    [0n, 1n],
  ];
}

function ZM2_mul(A: bigint[][], B: bigint[][]): bigint[][] {
  return [
    [
      A[0]![0]! * B[0]![0]! + A[0]![1]! * B[1]![0]!,
      A[0]![0]! * B[0]![1]! + A[0]![1]! * B[1]![1]!,
    ],
    [
      A[1]![0]! * B[0]![0]! + A[1]![1]! * B[1]![0]!,
      A[1]![0]! * B[0]![1]! + A[1]![1]! * B[1]![1]!,
    ],
  ];
}

const qfb_maxexpi = (Q: Qfb): number => 1 + Math.max(expi(Q.a), expi(Q.b), expi(Q.c));

function qfb_minexpi(Q: Qfb): number {
  const m = Math.min(expi(Q.a), expi(Q.b), expi(Q.c));
  return m < 0 ? 0 : m;
}

/** PARI `qfb3_SL2_apply` (`Qfb.c:760-776`) */
function qfb3_SL2_apply(a: bigint, b: bigint, c: bigint, M: bigint[][]): [bigint, bigint, bigint] {
  const x = M[0]![0]!,
    y = M[1]![0]!,
    z = M[0]![1]!,
    t = M[1]![1]!;
  const by = b * y,
    bt = b * t,
    bz = b * z;
  const a2 = shifti(a, 1),
    c2 = shifti(c, 1);
  const A1 = x * (a * x + by);
  const A2 = c * y * y;
  const B1 = x * (a2 * z + bt);
  const B2 = y * (c2 * t + bz);
  const C1 = z * (a * z + bt);
  const C2 = c * t * t;
  return [A1 + A2, B1 + B2, C1 + C2];
}

/** PARI `pqfbred_rec` (`Qfb.c:778-823`) */
function pqfbred_rec(Q: Qfb, m: number): { Q: Qfb; U: bigint[][] } {
  const d = Q.D;
  const n = qfb_maxexpi(Q) - m;
  let goingToR8 = false;
  if (n < 170) {
    const U = matid2();
    const R = pqfbred_iter_1(Q, m, U);
    return { Q: R, U };
  }
  let U: bigint[][];
  let QR: Qfb;
  if (qfb_minexpi(Q) <= m + 2) {
    U = matid2();
    QR = Q;
  } else {
    let p: number;
    let mm: number;
    let Q0: [bigint, bigint, bigint] | null = null;
    let Q1: Qfb;
    if (m <= n) {
      mm = m;
      p = 0;
      Q1 = Q;
    } else {
      mm = n;
      p = m + 1 - n;
      Q0 = [remi2n(Q.a, p), remi2n(Q.b, p), remi2n(Q.c, p)];
      const a1 = shifti(Q.a, -p),
        b1 = shifti(Q.b, -p),
        c1 = shifti(Q.c, -p);
      Q1 = mkqfb(a1, b1, c1, qfb_disc3(a1, b1, c1));
    }
    const h = mm + (n >> 1);
    if (qfb_minexpi(Q1) <= h) {
      U = matid2();
      QR = Q1;
    } else {
      const r = pqfbred_rec(Q1, h);
      QR = r.Q;
      U = r.U;
    }
    while (qfb_maxexpi(QR) > h) {
      if (is_minimal(QR, mm)) {
        goingToR8 = true;
        break;
      }
      QR = pqfbred_1(QR, mm, U);
    }
    if (!goingToR8) {
      const r = pqfbred_rec(QR, mm);
      QR = r.Q;
      U = ZM2_mul(U, r.U);
    }
    if (p > 0 && Q0) {
      const [x0, y0, z0] = qfb3_SL2_apply(Q0[0], Q0[1], Q0[2], U);
      QR = mkqfb(shifti(QR.a, p) + x0, shifti(QR.b, p) + y0, shifti(QR.c, p) + z0, d);
    }
  }
  QR = pqfbred_iter_1(QR, m, U);
  return { Q: QR, U };
}

/**
 * PARI `qfi_red_fast` (`Qfb.c:746-752`): should we use Schoenhage's
 * asymptotically fast reduction?  Exposed so the tests can exercise both paths.
 */
export function qfi_red_fast(Q: Qfb, limit = 9000): boolean {
  return 2 * qfb_maxexpi(Q) - expi(Q.D) > limit;
}

/* ------------------------------------------------------------------ */
/* Top-level reduction                                                 */
/* ------------------------------------------------------------------ */

/** conjugate: `(a,b,c) -> (a,-b,c)`, the inverse class */
function qfb_conj(q: Qfb): Qfb {
  return mkqfb(q.a, -q.b, q.c, q.D);
}

/**
 * PARI `qfi_red_av` (`Qfb.c:975-987`).
 *
 * DEVIATION (sign of b on the Schoenhage path). `pqfbred_rec` requires
 * non-negative coefficients, so PARI flips the sign of `b` before calling it —
 * but it never flips the result back, so above the `qfi_red_fast` threshold
 * upstream returns the reduced form of the *conjugate* (inverse) class whenever
 * `b < 0`. We conjugate the result back so that the fast path returns exactly
 * what the base case returns; the `*_schoenhage` tests check this on every
 * sample form with the threshold forced to -infinity.
 */
function qfi_red(Q: Qfb, limit = 9000): Qfb {
  if (qfi_red_fast(Q, limit)) {
    const flip = Q.b < 0n;
    let R = flip ? qfb_conj(Q) : Q;
    R = pqfbred_rec(R, 0).Q;
    if (flip) R = qfb_conj(R);
    return qfi_red_basecase(R);
  }
  return qfi_red_basecase(Q);
}

/** PARI `qfr_red_i` (`Qfb.c:914-940`) with `flag` restricted to t_QFB inputs */
function qfr_red_i(Q: Qfb, step: boolean, isqrtD: bigint | null, limit = 9000): Qfb {
  const rd = isqrtD ?? sqrti(Q.D);
  if (!qfi_red_fast(Q, limit)) {
    if (step) {
      const [a, b, c] = qfr3_rho(Q.a, Q.b, Q.c, rd);
      return mkqfb(a, b, c, Q.D);
    }
    const [a, b, c] = qfr3_red(Q.a, Q.b, Q.c, rd);
    return mkqfb(a, b, c, Q.D);
  }
  let a = Q.a,
    b = Q.b,
    c = Q.c;
  const d = Q.D;
  const sa = isign(a);
  if (sa < 0) {
    a = -a;
    b = -b;
    c = -c;
  }
  if (c < 0n) {
    const t = truedivii(rd - b, shifti(a, 1)) + 1n;
    const at = a * t;
    c = c - b * t + at * t;
    b = b - shifti(at, 1);
  }
  const flip = b < 0n;
  let Qr = pqfbred_rec(mkqfb(a, iabs(b), c, d), 0).Q;
  if (flip) Qr = qfb_conj(Qr); /* see the deviation note on qfi_red */
  if (sa < 0) Qr = mkqfb(-Qr.a, -Qr.b, -Qr.c, Qr.D);
  if (step) {
    const [a2, b2, c2] = qfr3_rho(Qr.a, Qr.b, Qr.c, rd);
    return mkqfb(a2, b2, c2, d);
  }
  const [a2, b2, c2] = qfr3_red(Qr.a, Qr.b, Qr.c, rd);
  return mkqfb(a2, b2, c2, d);
}

/**
 * PARI `qfbred0(x, flag, isqrtD, sqrtD)` (`Qfb.c:991-1002`) restricted to
 * `t_QFB` inputs (no Shanks distance).
 *
 * @param flag bit 0 (`qf_STEP`): perform a single reduction step (`rho`);
 *             bit 1 (`qf_NOD`) is implied for `t_QFB` inputs and ignored.
 */
export function qfbred(q: Qfb, flag = 0, isqrtD: bigint | null = null): Qfb {
  if (flag < 0 || flag > 3) throw new PariFlagError('qfbred');
  const step = (flag & 1) !== 0;
  if (qfb_is_qfi(q)) return step ? qfi_rho(q) : qfi_red(q);
  return qfr_red_i(q, step, isqrtD);
}

/** internal: reduce with a lowered Schoenhage threshold (test hook) */
export function qfbred_withLimit(q: Qfb, limit: number): Qfb {
  if (qfb_is_qfi(q)) return qfi_red(q, limit);
  return qfr_red_i(q, false, null, limit);
}

/**
 * PARI `qfi_redsl2` (`Qfb.c:857-880`).
 *
 * DEVIATION, same as in {@link qfi_red}: upstream negates `b` before the
 * Schoenhage reduction and then only negates the second *row* of `U`, which
 * yields `det U = -1` and the reduced form of the conjugate class. We instead
 * conjugate the intermediate form back and use `D*U*D` (`D = diag(1,-1)`),
 * which keeps `det U = 1` and `Q o U == result`, matching the base case.
 */
function qfi_redsl2(Q: Qfb, limit = 9000): { Q: Qfb; U: bigint[][] } {
  if (!qfi_red_fast(Q, limit)) return qfi_redsl2_basecase(Q);
  const flip = Q.b < 0n;
  const r = pqfbred_rec(flip ? qfb_conj(Q) : Q, 0);
  let RQ = r.Q;
  let RU = r.U;
  if (flip) {
    RQ = qfb_conj(RQ);
    RU = [
      [RU[0]![0]!, -RU[0]![1]!],
      [-RU[1]![0]!, RU[1]![1]!],
    ];
  }
  const w = qfi_redsl2_basecase(RQ);
  return { Q: w.Q, U: ZM2_mul(RU, w.U) };
}

/**
 * PARI `qfr_redsl2` (`Qfb.c:825-855`).
 *
 * DEVIATIONS on the Schoenhage path, same spirit as {@link qfi_redsl2}:
 *  - the `b < 0` conjugation is undone on the form as well as on `U`;
 *  - upstream also negates all three coefficients when `a < 0` *without*
 *    touching `U`, which breaks the `Q o U == result` contract; we simply fall
 *    back to the (always correct) base case in that situation.
 */
function qfr_redsl2(Q: Qfb, isqrtD: bigint, limit = 9000): { Q: Qfb; U: bigint[][] } {
  if (!qfi_red_fast(Q, limit) || Q.a < 0n) return qfr_redsl2_basecase(Q, isqrtD);
  let a = Q.a,
    b = Q.b,
    c = Q.c;
  const d = Q.D;
  const sa = isign(a);
  let t: bigint | null = null;
  if (sa < 0) {
    a = -a;
    b = -b;
    c = -c;
  }
  if (c < 0n) {
    t = truedivii(isqrtD - b, shifti(a, 1)) + 1n;
    const at = a * t;
    c = c - b * t + at * t;
    b = b - shifti(at, 1);
  }
  const sb = isign(b);
  const r = pqfbred_rec(mkqfb(a, sb < 0 ? -b : b, c, d), 0);
  let Qr = r.Q;
  let U = r.U;
  if (sa < 0) Qr = mkqfb(-Qr.a, -Qr.b, -Qr.c, Qr.D);
  if (sb < 0) {
    Qr = qfb_conj(Qr);
    U = [
      [U[0]![0]!, -U[0]![1]!],
      [-U[1]![0]!, U[1]![1]!],
    ];
  }
  if (t !== null) U = [[U[0]![0]! - U[1]![0]! * t, U[0]![1]! - U[1]![1]! * t], U[1]!];
  const W = qfr_redsl2_basecase(Qr, isqrtD);
  return { Q: W.Q, U: ZM2_mul(U, W.U) };
}

/**
 * PARI `qfbredsl2(q, isD)` (`Qfb.c:889-902`): return the reduced form together
 * with the base change `U` in SL_2(Z). `U` is returned row-major:
 * `U[0] = [u1, v1]`, `U[1] = [u2, v2]`, and `q o U == reduced form`.
 */
export function qfbredsl2(q: Qfb, isD: bigint | null = null): { Q: Qfb; U: bigint[][] } {
  if (qfb_is_qfi(q)) {
    if (isD !== null) throw new PariTypeError('qfbredsl2', 'isD given for a definite form');
    return qfi_redsl2(q);
  }
  return qfr_redsl2(q, isD ?? sqrti(q.D));
}

/** internal test hook: force the Schoenhage path by lowering the threshold */
export function qfbredsl2_withLimit(q: Qfb, limit: number): { Q: Qfb; U: bigint[][] } {
  if (qfb_is_qfi(q)) return qfi_redsl2(q, limit);
  return qfr_redsl2(q, sqrti(q.D), limit);
}

/** apply `M` in GL_2(Z) to the form `q` (`Qfb.c:760-776`) */
export function qfb_apply(q: Qfb, M: bigint[][]): Qfb {
  const [a, b, c] = qfb3_SL2_apply(q.a, q.b, q.c, M);
  return mkqfb(a, b, c, qfb_disc3(a, b, c));
}

/* ------------------------------------------------------------------ */
/* Public composition / powering                                       */
/* ------------------------------------------------------------------ */

/** PARI `qfbcompraw` (`Qfb.c:1165-1181`) */
export function qfbcompraw(x: Qfb, y: Qfb): Qfb {
  if (x.D !== y.D) {
    const z = qfb_comp_gen(x, y);
    if (!z) throw new PariDomainError('qfbcompraw', 'discriminants', '!=', 'equal');
    return z;
  }
  return qfb_comp(x, y);
}

/** PARI `qfbcomp` (`Qfb.c:1145-1161`) */
export function qfbcomp(x: Qfb, y: Qfb): Qfb {
  if (x.D !== y.D) {
    const z = qfb_comp_gen(x, y);
    if (!z) throw new PariDomainError('qfbcomp', 'discriminants', '!=', 'equal');
    return qfbred(z);
  }
  return qfbred(qfb_comp(x, y));
}

/** PARI `qfbsqr` (`Qfb.c:1208-1213`) */
export function qfbsqr(x: Qfb): Qfb {
  return qfbred(qfb_sqr(x));
}

/** PARI `qfbsqr` with `raw = 1` (`qfisqr0`/`qfrsqr0`, `Qfb.c:1183-1204`) */
export function qfbsqrraw(x: Qfb): Qfb {
  return qfb_sqr(x);
}

/**
 * PARI `qfbpowraw` (`Qfb.c:1634-1639`): `x^n` without reduction.
 *
 * NB the two branches use *different* addition chains in PARI, and raw
 * composition is not canonical, so the chain matters: definite forms go through
 * `gen_powu` (left-to-right binary, `bb_group.c:120-153`) while indefinite forms
 * go through `qfr3_powraw` (right-to-left, `Qfb.c:1508-1519`). We mirror both.
 */
export function qfbpowraw(x: Qfb, n: bigint): Qfb {
  if (n === 0n) return qfb_1(x);
  if (n === 1n) return x;
  if (n === -1n) return qfbinv(x);
  const base = n < 0n ? qfbinv(x) : x;
  const e = n < 0n ? -n : n;
  if (qfb_is_qfi(base)) return leftrightPow(base, e, qfb_sqr, qfb_comp);
  return rightleftPow(base, e, qfb_sqr, qfb_comp);
}

/** left-to-right binary powering, PARI `leftright_binary_powu` (`bb_group.c:120-144`) */
function leftrightPow(
  x: Qfb,
  n: bigint,
  sqr: (a: Qfb) => Qfb,
  mul: (a: Qfb, b: Qfb) => Qfb
): Qfb {
  if (n === 1n) return x;
  const bits = n.toString(2);
  let y = x;
  for (let i = 1; i < bits.length; i++) {
    y = sqr(y);
    if (bits[i] === '1') y = mul(y, x);
  }
  return y;
}

/** right-to-left binary powering, PARI `qfr3_pow`/`qfr3_powraw` (`Qfb.c:1508-1573`) */
function rightleftPow(
  x: Qfb,
  n: bigint,
  sqr: (a: Qfb) => Qfb,
  mul: (a: Qfb, b: Qfb) => Qfb
): Qfb {
  let y: Qfb | null = null;
  let b = x;
  let e = n;
  for (;;) {
    if (e & 1n) y = y === null ? b : mul(y, b);
    e >>= 1n;
    if (e === 0n) break;
    b = sqr(b);
  }
  return y!;
}

/**
 * PARI `qfbpow` (`Qfb.c:1643-1648`): `x^n` with reduction after every step.
 *
 * Definite forms: `qfipow` (`Qfb.c:1311-1321`) reduces the base first and uses
 * `gen_pow`; the reduced form of a definite class is unique so the addition
 * chain is irrelevant. Indefinite forms: `qfrpow` (`Qfb.c:1608-1632`) does
 * **not** reduce the base and uses the right-to-left chain of `qfr3_pow`; the
 * reduced representative of an indefinite class is not unique, so we reproduce
 * that chain exactly.
 */
export function qfbpow(x: Qfb, n: bigint): Qfb {
  if (n === 0n) return qfb_1(x);
  const base0 = n < 0n ? qfbinv(x) : x;
  const e = n < 0n ? -n : n;
  if (qfb_is_qfi(base0)) {
    const base = qfbred(base0);
    if (e === 1n) return base;
    return leftrightPow(
      base,
      e,
      (a) => qfbred(qfb_sqr(a)),
      (a, b) => qfbred(qfb_comp(a, b))
    );
  }
  if (e === 1n) return qfbred(base0);
  return rightleftPow(
    base0,
    e,
    (a) => qfbred(qfb_sqr(a)),
    (a, b) => qfbred(qfb_comp(a, b))
  );
}

/* ------------------------------------------------------------------ */
/* Prime forms (Qfb.c:1661-1749)                                       */
/* ------------------------------------------------------------------ */

/** PARI `primeform(x, p)` (`Qfb.c:1696-1749`) for `p > 0` */
export function primeform(x: bigint, p: bigint): Qfb {
  const f = 'primeform';
  if (x === 0n) throw new PariDomainError(f, 'D', '=', '0');
  if (p === 0n) throw new PariDomainError(f, 'p', '=', '0');
  if (p < 0n) throw new NotImplementedError('primeform with p < 0 (negative definite t_QFB)');
  if (p === 1n) return qfb_1(mkqfb(1n, 0n, 0n, x));
  let s = mod8(x);
  if (x < 0n && s) s = 8 - s;
  if (s & 2) throw new PariDomainError(f, 'disc % 4', '>', '1');
  if (p === 2n) {
    let b: bigint;
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
        throw new PariSqrtnError(f, `Mod(${x}, 2)`);
    }
    const c = shifti(BigInt(s) - x, -3);
    return mkqfb(2n, b, c, x);
  }
  const sq = Fp_sqrt(Fp_red(x, p), p);
  if (sq === null) throw new PariSqrtnError(f, `Mod(${x}, ${p})`);
  let b = sq;
  const s1 = s & 1;
  if ((b === 0n && s1) || mod2(b) !== s1) b = p - b;
  const c = shifti(b * b - x, -2) / p;
  return mkqfb(p, b, c, x);
}

/* ------------------------------------------------------------------ */
/* Cornacchia (Qfb.c:1997-2085)                                        */
/* ------------------------------------------------------------------ */

/**
 * The Euclidean remainder sequence of `(A,B)`, largest first. PARI reaches the
 * relevant remainders through `halfgcdii` (`kernel/none/halfgcd.c:385-397`),
 * which returns the first pair `(a,b)` of consecutive remainders with
 * `b^2 < max(|A|,|B|) <= a^2`; we simply walk the sequence, which produces the
 * same remainders.
 */
function remainderSequence(A: bigint, B: bigint): bigint[] {
  const out: bigint[] = [A, B];
  let a = A,
    b = B;
  while (b !== 0n) {
    const r = a % b;
    a = b;
    b = r;
    if (b !== 0n) out.push(b);
  }
  return out;
}

/**
 * PARI `cornacchia(d, p, &x, &y)` (`Qfb.c:1999-2017`): solve `x^2 + d y^2 = p`.
 * Assumes `d > 0` and `p` prime. Returns `null` when there is no solution.
 */
export function cornacchia(d: bigint, p: bigint): [bigint, bigint] | null {
  const bb = p - d;
  if (bb < 0n) return null;
  if (bb === 0n) return [0n, 1n];
  let b = Fp_sqrt(Fp_red(bb, p), p); /* sqrt(-d) mod p */
  if (b === null) return null;
  /* first remainder with square < p */
  const seq = remainderSequence(p, b);
  let r: bigint | null = null;
  for (const v of seq)
    if (v * v < p) {
      r = v;
      break;
    }
  if (r === null) return null;
  b = r;
  const num = p - b * b;
  if (num % d !== 0n) return null;
  const [ok, c] = Z_issquareall(num / d);
  if (!ok) return null;
  return [b, c];
}

/**
 * PARI `cornacchia2(d, p, &x, &y)` (`Qfb.c:2060-2081`): solve `x^2 + d y^2 = 4p`.
 * Assumes `d > 0` congruent to 0 or 3 mod 4 and `p` prime.
 */
export function cornacchia2(d: bigint, p: bigint): [bigint, bigint] | null {
  const p4 = shifti(p, 2);
  if (iabs(p4) < d) return null;
  if (iabs(p) === 2n) {
    if (d === 4n) return [2n, 1n];
    if (d === 7n) return [1n, 1n];
    return null;
  }
  let b = Fp_sqrt(Fp_red(-d, p), p);
  if (b === null) return null;
  if (b === 0n) {
    /* d = p, 2p, 3p, 4p */
    if (iabs(d) === iabs(p4)) return [0n, 1n];
    if (iabs(d) === iabs(p)) return [0n, 2n];
    return null;
  }
  if (mod2(b) !== mod2(d)) b = p - b;
  /* largest remainder of the (2p, b) sequence whose square is <= 4p */
  const seq = remainderSequence(shifti(p, 1), b);
  let r: bigint | null = null;
  for (const v of seq)
    if (v * v <= p4) {
      r = v;
      break;
    }
  if (r === null) return null;
  const num = p4 - r * r;
  if (num % d !== 0n) return null;
  const [ok, c] = Z_issquareall(num / d);
  if (!ok) return null;
  return [r, c];
}

/** PARI `qfbcornacchia(d, p)` (`Qfb.c:2087-2099`) */
export function qfbcornacchia(d: bigint, p: bigint): [bigint, bigint] | null {
  if (d <= 0n) throw new PariTypeError('qfbcornacchia', 'd <= 0');
  if (p < 2n) throw new PariTypeError('qfbcornacchia', 'p < 2');
  return mod4(p) ? cornacchia(d, p) : cornacchia2(d, shifti(p, -2));
}

/* ------------------------------------------------------------------ */
/* Zn_quad_roots (quad.c:1149-1260)                                    */
/* ------------------------------------------------------------------ */

/** PARI `Z2_sqrt(x, e)` (`trans1.c:1410-1445`): a square root of the unit `x` mod `2^e`. */
export function Z2_sqrt(x: bigint, e: number): bigint | null {
  const r = x >= 0n ? mod16(x) : 16 - mod16(x);
  switch (e) {
    case 1:
      return 1n;
    case 2:
      return (r & 3) === 1 ? 1n : null;
    case 3:
      return (r & 7) === 1 ? 1n : null;
    case 4:
      if (r === 1) return 1n;
      return r === 9 ? 3n : null;
    default:
      if ((r & 7) !== 1) return null;
  }
  let z = r === 1 ? 1n : 3n;
  let ez = 3;
  for (;;) {
    ez = (ez << 1) - 1;
    if (ez > e) ez = e;
    const m = 1n << BigInt(ez);
    z = z + remi2n(x * Fp_inv(z, m), ez);
    z = shifti(z, -1);
    if (e === ez) return z;
    if (ez < e) ez--;
  }
}

/** PARI `Zp_sqrt(x, p, e)` (`Zp.c:203-213`): Hensel lift of `Fp_sqrt`. */
export function Zp_sqrt(x: bigint, p: bigint, e: number): bigint | null {
  if (p === 2n) return Z2_sqrt(x, e);
  let z = Fp_sqrt(Fp_red(x, p), p);
  if (z === null) return null;
  let k = 1;
  let q = p;
  while (k < e) {
    const k2 = Math.min(2 * k, e);
    const q2 = p ** BigInt(k2);
    /* Newton: z <- z - (z^2 - x)/(2z) mod p^k2 */
    const inv = Fp_inv(2n * z, q2);
    z = Fp_red(z - (z * z - x) * inv, q2);
    k = k2;
    q = q2;
  }
  void q;
  return z;
}

/** CRT: solve `x = a[i] mod m[i]` for pairwise coprime moduli. */
function crt(a: bigint[], m: bigint[]): [bigint, bigint] {
  let x = 0n;
  let M = 1n;
  for (let i = 0; i < a.length; i++) {
    const mi = m[i]!;
    const [g, u] = bezout(M, mi);
    if (g !== 1n) throw new Error('crt: moduli not coprime');
    const t = Fp_red((a[i]! - x) * u, mi);
    x = x + M * t;
    M = M * mi;
    x = Fp_red(x, M);
  }
  return [x, M];
}

function Z_pvalrem(x: bigint, p: bigint): [number, bigint] {
  let v = 0;
  while (x % p === 0n) {
    x /= p;
    v++;
  }
  return [v, x];
}

/**
 * PARI `Zn_quad_roots(N, B, C)` (`quad.c:1149-1260`): return `[N', v]` where `v`
 * lists all `x` mod `N'` with `x^2 + B x + C == 0` mod `N`.
 * `N` is given by its factorisation.
 */
export function Zn_quad_roots(
  fa0: Factorization,
  B: bigint,
  C: bigint
): { Np: bigint; roots: bigint[] } | null {
  let N = 1n;
  for (const [p, e] of fa0) N *= p ** e;
  N = iabs(N);
  /* PARI: `fa = clean_Z_factor(fa)` (`quad.c:1158`, `arith2.c:273`) drops the
   * leading `-1` that `Z_factor` records for a negative argument.  Without this
   * the unit survives into the `4N` factorisation below, where `Z_pvalrem(D, -1)`
   * never terminates and `Np` would pick up a spurious sign. */
  const fa: Factorization = fa0.length > 0 && fa0[0]![0] === -1n ? fa0.slice(1) : fa0;
  const N4 = shifti(N, 2);
  const D = Fp_red(B * B - shifti(C, 2), N4);
  if (D === 0n) {
    /* (x + B/2)^2 = 0 (mod N) */
    let Np = 1n;
    for (const [p, e] of fa) Np *= p ** ((e + 1n) >> 1n);
    const Bh = shifti(B, -1);
    return { Np, roots: [Fp_red(-Bh, Np)] };
  }
  /* factorisation of 4N */
  const fa4: Array<[bigint, bigint]> = [];
  let has2 = false;
  for (const [p, e] of fa) {
    if (p === 2n) {
      fa4.push([2n, e + 2n]);
      has2 = true;
    } else fa4.push([p, e]);
  }
  if (!has2) fa4.push([2n, 2n]);
  fa4.sort((u, v) => (u[0] < v[0] ? -1 : u[0] > v[0] ? 1 : 0));

  const F: bigint[] = [];
  const mF: bigint[] = [];
  const Qm: bigint[] = [];
  let Q0 = 1n;
  let F0: bigint | null = 0n;
  let F0set = true;
  for (const [p, es] of fa4) {
    const s = Number(es);
    const [t, D0] = Z_pvalrem(D, p);
    const dd = s - t;
    if (dd <= 0) {
      Q0 = Q0 * p ** BigInt((s + 1) >> 1);
      continue;
    }
    if (t % 2 === 1) return null;
    const t2 = t >> 1;
    let q: bigint, f: bigint;
    if (p !== 2n) {
      if (kronecker(D0, p) === -1) return null;
      q = p ** BigInt(s - t2);
      const r = Zp_sqrt(D0, p, dd);
      if (r === null) return null;
      f = t2 ? p ** BigInt(t2) * r : r;
    } else {
      if (dd <= 3) {
        if (dd === 3 && Mod8(D0) !== 1) return null;
        if (dd === 2 && Mod4(D0) !== 1) return null;
        Q0 = 1n << BigInt(1 + t2);
        F0 = null;
        F0set = false;
        continue;
      }
      if (Mod8(D0) !== 1) return null;
      q = 1n << BigInt(dd - 1 + t2);
      const r = Z2_sqrt(D0, dd);
      if (r === null) return null;
      f = t2 ? shifti(r, t2) : r;
    }
    Qm.push(q);
    F.push(Fp_red(f, q));
    mF.push(Fp_red(-f, q));
  }

  const A: bigint[] = F.slice();
  const mods = Qm.slice();
  if (Q0 !== 1n) {
    if (!F0set || F0 === null) F0 = shifti(Q0, -1);
    A.push(F0);
    mF.push(F0);
    mods.push(Q0);
  }
  const nfree = F.length;
  const ct = 1 << nfree;
  let Np = 1n;
  for (const q of mods) Np *= q;
  let Bm: bigint | null = Fp_red(B, Np);
  if (Bm === 0n) Bm = null;
  const Nphalf = shifti(Np, -1);
  const roots: bigint[] = [];
  for (let j = 1; j <= ct; j++) {
    let m = j - 1;
    const AA = A.slice();
    for (let i = 0; i < nfree; i++) {
      AA[i] = m & 1 ? mF[i]! : F[i]!;
      m >>= 1;
    }
    let [u] = crt(AA, mods);
    if (Bm !== null) u = u - Bm;
    roots.push(Fp_red(shifti(u, -1), Nphalf));
  }
  return { Np: Nphalf, roots };
}

/* ------------------------------------------------------------------ */
/* normforms + qfbsolve (Qfb.c:1751-1995)                              */
/* ------------------------------------------------------------------ */

function factorback(fa: Factorization): bigint {
  let n = 1n;
  for (const [p, e] of fa) n *= p ** e;
  return n;
}

/**
 * PARI `normforms(D, fa)` (`Qfb.c:1751-1783`). `arg` mirrors PARI's polymorphic
 * second argument: the value `a` is read off without factoring when `n` is a
 * bare integer, so `n = 0` short-circuits before any factorisation (as in PARI).
 */
function normforms(D: bigint, arg: SolveArg): Qfb[] | null {
  const a = arg.n !== null ? arg.n : factorback(arg.fa!);
  const sa = isign(a);
  if (sa === 0 || (D < 0n && sa < 0)) return null;
  const fa = solveArgFa(arg);
  const D_odd = mod2(D) !== 0;
  const V = D_odd
    ? Zn_quad_roots(fa, 1n, shifti(1n - D, -2))
    : Zn_quad_roots(fa, 0n, -shifti(D, -2));
  if (!V) return null;
  const N = V.Np;
  const Bs = V.roots;
  const N2 = shifti(N, 1);
  /* PARI (`Qfb.c:1766`): `aN = itou(diviiexact(a, N))`, annotated there as
   * "|a|/N" - `itou` returns the unsigned mantissa, i.e. the absolute value.
   * `aN` is only ever a loop count and an exact divisor of `b + N`, so the sign
   * must be dropped here; the signed quotient makes the `++j === aN` loop below
   * run forever whenever `a < 0` (an indefinite form with `n < 0`). */
  const aN = iabs(a) / N;
  const L: Qfb[] = [];
  for (const B of Bs) {
    let b = shifti(B, 1);
    if (D_odd) b = b + 1n;
    let c = shifti(b * b - D, -2) / a;
    for (let j = 0; ; b = b + N2) {
      L.push(mkqfb(a, b, c, D));
      if (++j === Number(aN)) break;
      let Cc = b + N;
      if (aN > 1n) Cc = Cc / aN;
      c = sa > 0 ? c + Cc : c - Cc;
    }
  }
  return L;
}

/** PARI `SL2_div_mul_e1(N, M)` (`Qfb.c:1786-1794`): `(N*M^-1)[,1]` */
function SL2_div_mul_e1(N: bigint[][], M: bigint[][]): [bigint, bigint] {
  const b = M[1]![0]!,
    d = M[1]![1]!;
  const A = N[0]![0]! * d,
    B = N[0]![1]! * b;
  const C = N[1]![0]! * d,
    Dd = N[1]![1]! * b;
  return [A - B, C - Dd];
}

/** PARI `qfisolve_normform` (`Qfb.c:1795-1801`) */
function qfisolve_normform(
  Qr: { Q: Qfb; U: bigint[][] },
  P: Qfb
): [bigint, bigint] | null {
  const { Q: b, U: M } = qfi_redsl2_basecase(P);
  if (!qfb_equal(Qr.Q, b)) return null;
  return SL2_div_mul_e1(Qr.U, M);
}

/** PARI `GL2_qfb_equal` (`Qfb.c:1804-1810`) */
function GL2_qfb_equal(a: Qfb, b: Qfb): boolean {
  return a.a === b.a && iabs(a.b) === iabs(b.b) && a.c === b.c;
}

/** PARI `allsols` (`Qfb.c:1813-1829`) */
function allsols(
  Q: Qfb,
  s: number,
  u: bigint,
  v: bigint
): [bigint, bigint] | Array<[bigint, bigint]> {
  if (v < 0n) {
    u = -u;
    v = -v;
  }
  const w: [bigint, bigint] = [u, v];
  if (s < 0) return w;
  if (s === 0) return [w];
  const b = Q.b;
  if (b !== 0n) {
    const [t, r] = dvmdii(b * v, Q.a);
    if (r !== 0n) return [w];
    u = u + t;
  }
  return [w, [-u, v]];
}

/** PARI `qfisolvep_all` (`Qfb.c:1830-1876`) */
function qfisolvep_all(
  Q: Qfb,
  p: bigint,
  all: boolean
): [bigint, bigint] | Array<[bigint, bigint]> | null {
  const D = Q.D;
  let s = kronecker(D, p);
  if (s < 0) return null;
  if (!all) s = -1;
  if (Q.b === 0n) {
    if (Q.a === 1n) {
      const r = cornacchia(Q.c, p);
      if (!r) return null;
      return allsols(Q, s, r[0], r[1]);
    }
    if (Q.c === 1n) {
      const r = cornacchia(Q.a, p);
      if (!r) return null;
      return allsols(Q, s, r[1], r[0]);
    }
  }
  const RU = qfi_redsl2_basecase(Q);
  const R = RU.Q;
  const U = RU.U;
  if (R.a === 1n) {
    let x: [bigint, bigint];
    if (R.b === 0n) {
      const r = cornacchia(R.c, p);
      if (!r) return null;
      x = [r[0], r[1]];
    } else {
      const r = cornacchia2(-D, p);
      if (!r) return null;
      let t = r[0] - r[1];
      if (mod2(t) !== 0) return null;
      x = [shifti(t, -1), r[1]];
    }
    /* x <- U * x (column), then transposed to a row vector */
    const x0 = U[0]![0]! * x[0] + U[0]![1]! * x[1];
    const x1 = U[1]![0]! * x[0] + U[1]![1]! * x[1];
    return allsols(Q, s, x0, x1);
  }
  const qV = qfi_redsl2_basecase(primeform(D, p));
  const q = qV.Q;
  let V = qV.U;
  if (!GL2_qfb_equal(R, q)) return null;
  if (isign(R.b) !== isign(q.b)) V = [V[0]!, [-V[1]![0]!, V[1]![1]!]];
  const x = SL2_div_mul_e1(U, V);
  return allsols(Q, s, x[0], x[1]);
}

/** PARI `qfrsolve_normform` (`Qfb.c:1878-1895`) */
function qfrsolve_normform(
  N: { Q: Qfb; U: bigint[][] },
  Ps: Qfb,
  rd: bigint
): [bigint, bigint] | null {
  let M = N;
  const P = qfr_redsl2_basecase(Ps, rd);
  let Qc = P;
  for (;;) {
    if (qfb_equal(M.Q, P.Q)) return SL2_div_mul_e1(M.U, P.U);
    if (qfb_equal(N.Q, Qc.Q)) return SL2_div_mul_e1(N.U, Qc.U);
    M = qfr_rhosl2(M, rd);
    if (qfb_equal(M.Q, N.Q)) return null;
    Qc = qfr_rhosl2(Qc, rd);
    if (qfb_equal(P.Q, Qc.Q)) return null;
  }
}

/**
 * PARI `known_prime` (`Qfb.c:1897-1906`). When `n` is given as a plain integer
 * PARI only runs a BPSW test (`check_arith_all` returns NULL, so no factoring);
 * when it is given as a factorisation, the factorisation must be a single prime
 * to the first power.
 */
function known_prime(arg: SolveArg): bigint | null {
  if (!arg.fa) return arg.n !== null && isPrime(arg.n) ? arg.n : null;
  if (arg.fa.length !== 1) return null;
  const [p, e] = arg.fa[0]!;
  return e === 1n && p > 1n ? p : null;
}

/**
 * `n` as PARI's `check_arith_all` sees it: either a bare integer (`fa === null`,
 * no factorisation known yet) or an explicit factorisation.
 */
interface SolveArg {
  n: bigint | null;
  fa: Factorization | null;
}

/** the factorisation of `arg`, computing it if necessary */
function solveArgFa(arg: SolveArg): Factorization {
  if (arg.fa) return arg.fa;
  const fa = Z_factor(arg.n!);
  arg.fa = fa;
  return fa;
}

type SolveState = { Qr: { Q: Qfb; U: bigint[][] } | null };

/** PARI `qfbsolve_primitive_i` (`Qfb.c:1908-1930`) */
function qfbsolve_primitive_i(
  Q: Qfb,
  rd: bigint | null,
  st: SolveState,
  arg: SolveArg,
  all: boolean
): [bigint, bigint] | Array<[bigint, bigint]> | null {
  if (rd === null) {
    const p = known_prime(arg);
    if (p !== null) return qfisolvep_all(Q, p, all);
  }
  const F = normforms(Q.D, arg);
  if (!F) return null;
  if (!st.Qr) st.Qr = qfbredsl2(Q, rd);
  const W: Array<[bigint, bigint]> = [];
  for (const f of F) {
    const x = rd !== null ? qfrsolve_normform(st.Qr, f, rd) : qfisolve_normform(st.Qr, f);
    if (x) {
      if (!all) return x;
      W.push(x);
    }
  }
  if (W.length === 0) return null;
  return lexsort(W);
}

function lexcmp(a: readonly bigint[], b: readonly bigint[]): number {
  for (let i = 0; i < Math.min(a.length, b.length); i++) {
    if (a[i]! < b[i]!) return -1;
    if (a[i]! > b[i]!) return 1;
  }
  return a.length - b.length;
}

function lexsort<T extends bigint[]>(v: T[]): T[] {
  return v.slice().sort(lexcmp);
}

/** PARI `qfb_initrd` (`Qfb.c:1932-1933`) */
function qfb_initrd(Q: Qfb): bigint | null {
  return Q.D > 0n ? sqrti(Q.D) : null;
}

/** divisors of `prod p^floor(e/2)`, with their factorisations (PARI `divisors_factored`) */
function divisorsFactored(fa: Factorization): Array<{ d: bigint; fa: Factorization }> {
  let out: Array<{ d: bigint; fa: Factorization }> = [{ d: 1n, fa: [] }];
  for (const [p, e] of fa) {
    const half = e >> 1n;
    const next: Array<{ d: bigint; fa: Factorization }> = [];
    for (const cur of out)
      for (let k = 0n; k <= half; k++)
        next.push({
          d: cur.d * p ** k,
          fa: k === 0n ? cur.fa : [...cur.fa, [p, k] as [bigint, bigint]],
        });
    out = next;
  }
  out.sort((u, v) => (u.d < v.d ? -1 : u.d > v.d ? 1 : 0));
  return out;
}

function famatDivSqr(fa: Factorization, g: Factorization): Factorization {
  const m = new Map<bigint, bigint>();
  for (const [p, e] of fa) m.set(p, (m.get(p) ?? 0n) + e);
  for (const [p, e] of g) m.set(p, (m.get(p) ?? 0n) - 2n * e);
  const out: Factorization = [];
  for (const [p, e] of m) if (e !== 0n) out.push([p, e]);
  out.sort((u, v) => (u[0] < v[0] ? -1 : u[0] > v[0] ? 1 : 0));
  return out;
}

/**
 * PARI `qfbsolve(Q, n, flag)` (`Qfb.c:1987-1995`): solve `Q(x,y) = n`.
 *
 * @param flag bit 0: return *all* solutions instead of a single one;
 *             bit 1: also look for imprimitive solutions.
 * @param fa   optional known factorisation of `n` (PARI accepts `[n, factor(n)]`);
 *             pass it to avoid an expensive factorisation, exactly as in PARI.
 * @returns `[x,y]` (or `[]` if there is none) when bit 0 is clear, otherwise a
 *          lexicographically sorted list of solutions.
 */
export function qfbsolve(
  Q: Qfb,
  n: bigint,
  flag = 0,
  fa: Factorization | null = null
): Array<bigint> | Array<[bigint, bigint]> {
  if (flag < 0 || flag > 3) throw new PariFlagError('qfbsolve');
  const all = (flag & 1) !== 0;
  if (flag & 2) return qfbsolve_all(Q, n, all, fa);
  return qfbsolve_primitive(Q, n, all, fa);
}

function qfbsolve_primitive(
  Q: Qfb,
  n: bigint,
  all: boolean,
  fa: Factorization | null
): Array<bigint> | Array<[bigint, bigint]> {
  const rdQ = qfb_initrd(Q);
  const st: SolveState = { Qr: null };
  const x = qfbsolve_primitive_i(Q, rdQ, st, { n, fa }, all);
  if (!x) return [];
  return x as Array<bigint> | Array<[bigint, bigint]>;
}

function qfbsolve_all(
  Q: Qfb,
  n: bigint,
  all: boolean,
  fa0: Factorization | null
): Array<bigint> | Array<[bigint, bigint]> {
  const fa = fa0 ?? Z_factor(n);
  const rdQ = qfb_initrd(Q);
  const st: SolveState = { Qr: null };
  const D = divisorsFactored(fa);
  const W: Array<[bigint, bigint]> = [];
  for (let i = 0; i < D.length; i++) {
    const d = D[i]!;
    const FA = d.d === 1n ? fa : famatDivSqr(fa, d.fa);
    const w = qfbsolve_primitive_i(Q, rdQ, st, { n: factorback(FA), fa: FA }, all);
    if (w) {
      if (!all) {
        const s = w as [bigint, bigint];
        return d.d === 1n ? s : [s[0] * d.d, s[1] * d.d];
      }
      for (const s of w as Array<[bigint, bigint]>)
        W.push(d.d === 1n ? s : [s[0] * d.d, s[1] * d.d]);
    }
  }
  if (W.length === 0) return [];
  return lexsort(W);
}
