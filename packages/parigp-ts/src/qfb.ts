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
 * PARI additionally carries Shanks' logarithmic distance for *indefinite*
 * forms, as a `t_VEC [t_QFB, t_REAL]` (`check_qfbext`, `Qfb.c:111-123`). That
 * family is ported too:
 *
 * - `Qfb.c:396-430`    `fix_expo`, `qfr5_dist` (the `(e, d)` distance coding)
 * - `Qfb.c:434-472`    `rho_get_BC`, `qfr3_rho`, `qfr5_rho`
 * - `Qfb.c:474-503`    `qfr_to_qfr5`, `qfr5_to_qfr` (`qfr3_to_qfr` is the
 *                      trailing `mkqfb` of `qfr5_to_qfr` here)
 * - `Qfb.c:505-550`    `ab_isreduced`, `qfr5_red`, `qfr3_red`
 * - `Qfb.c:552-591`    `qfr_data_init`, `qfr5_init`, `qfr3_init`
 * - `Qfb.c:594-611`    `qfr_red_basecase_i`
 * - `Qfb.c:1124-1137`  `qfrcomp0`, `Qfb.c:1194-1204` `qfrsqr0`
 * - `Qfb.c:1252-1277`  `qfr_1_fill`, `qfr5_1`, `qfr3_1`
 * - `Qfb.c:1470-1568`  `qfr5_compraw/comp/powraw/pow`, `qfr3_compraw/comp/powraw/pow`
 * - `Qfb.c:1570-1628`  `qfrinvraw`, `qfrpowraw`, `qfrpow`
 *
 * so {@link qfbred}, {@link qfbcomp}, {@link qfbcompraw}, {@link qfbsqr},
 * {@link qfbsqrraw}, {@link qfbpow} and {@link qfbpowraw} all accept a
 * {@link QfbExt} (`{ Q, d }`) wherever PARI accepts its `t_VEC`, and return one.
 * The `t_REAL` itself ({@link MpReal}) is a port of PARI's mp kernel; see the
 * long note on it below, in particular its rounding deviation.
 *
 * Everything else here is exact integer arithmetic.
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

/* ================================================================== */
/* t_REAL: PARI's multiprecision binary floating point                 */
/* ------------------------------------------------------------------ */
/*
 * The Shanks-distance family below (`qfr5_*`) carries the distance as a
 * PARI `t_REAL`, so we need PARI's float type.  A `t_REAL` is a
 * sign/exponent/mantissa triple: `signe(x)`, `expo(x)` and a mantissa of
 * `realprec(x)` bits whose top bit is set, so that
 *
 *     x = signe(x) * m * 2^(expo(x) + 1 - realprec(x)),   2^e <= |x| < 2^(e+1)
 *
 * and a "real zero" (`real_0_bit(e)`, `lg == 2`) which carries only the
 * exponent `e` = "0 known to absolute accuracy 2^e".  From PARI 2.16 on a
 * `prec` is measured in **bits** (`pariinl.h:1626-1640`: `prec2nbits(x) = x`,
 * `nbits2prec(x) = ceil(x/64)*64`), which is what we use here.
 *
 * Sources, ported function by function:
 *   `kernel/none/level1.h:435-450`   `real_0_bit`, `real_0`, `real_1`
 *   `kernel/none/mp_indep.c:139-600` `mul0r`, `div0r`, `mulsr`, `mulur`,
 *                                    `mulrr`, `sqrr`, `mulir`, `divir`,
 *                                    `divru`, `mulrrz_end` (the rounding rule)
 *   `kernel/none/mp.c:635-765`       `divrr`, `divri`
 *   `kernel/none/add.c:110-330`      `addir_sign`, `addsr_sign`, `addrr_sign`
 *   `kernel/none/mp.c:2063-2160`     `sqrtr_abs`
 *   `basemath/trans1.c:2841-2985`    `log2_split`, `constlog2`, `mplog2`,
 *                                    `logr_aux`, `logr_abs`
 *   `basemath/gen3.c:2669-2695`      `gcvtoi`
 *   `headers/parigen.h:122-140`      `lg`, `signe`, `expo`, `realprec`
 *
 * DEVIATION (rounding).  PARI's mp kernel is *nearly* correctly rounded: each
 * primitive keeps one guard word and rounds up when its top bit is set
 * (`mulrrz_end`, `mp_indep.c:216-222`), which can differ from the correctly
 * rounded result in the last bit.  We round the exact result to nearest with
 * ties away from zero -- the same rule, applied to the exact value instead of
 * to a one-word approximation of it.  Results therefore agree with PARI to
 * within a few units in the last place, not bit for bit.  PARI's own printed
 * output is already off by one unit in the 38th digit on the `test/in/qfb`
 * distance (see `qfb.test.ts`), so bit-compatibility is not even well defined.
 * The *output precision* of every operation is PARI's (see each function).
 */

/** PARI `t_REAL`. `s` is `signe`, `e` is `expo`, `p` is `realprec` in bits. */
export interface MpReal {
  /** `signe(x)` */
  readonly s: -1 | 0 | 1;
  /** `expo(x)`: `2^e <= |x| < 2^(e+1)` when `s != 0`; the accuracy when `s = 0` */
  readonly e: number;
  /** mantissa: exactly `p` bits with the top bit set (`0n` when `s = 0`) */
  readonly m: bigint;
  /** `realprec(x)`, in bits (`0` for a real zero, whose `lg` is 2) */
  readonly p: number;
}

/** bit length of `x > 0` (PARI `expi(x) + 1`) */
function bitlen(x: bigint): number {
  if (x <= 0n) return 0;
  let n = 0;
  let v = x;
  while (v >= 1n << 1024n) {
    v >>= 1024n;
    n += 1024;
  }
  while (v >= 1n << 32n) {
    v >>= 32n;
    n += 32;
  }
  let w = Number(v);
  while (w >= 1) {
    w = Math.floor(w / 2);
    n++;
  }
  return n;
}

/** PARI `nbits2prec` (`pariinl.h:1626`), in bits: round up to a whole word */
export function nbits2prec(n: number): number {
  return Math.ceil(n / 64) * 64;
}

/** PARI `realprec(x)` (`parigen.h:134`) */
export function realprec(x: MpReal): number {
  return x.p;
}
/** PARI `expo(x)` (`parigen.h:139`) */
export function expo(x: MpReal): number {
  return x.e;
}
/** PARI `signe(x)` for a t_REAL */
export function realsigne(x: MpReal): number {
  return x.s;
}

/** PARI `LOWDEFAULTPREC` (`parigen.h:47`) */
const LOWDEFAULTPREC = 64;

/**
 * PARI `precision(x)` for a t_REAL (`gen3.c:140-142`, `precREAL`/`prec0`):
 * the working precision `x` was computed at.  Unlike `realprec` this is
 * nonzero for a real zero, whose exponent records its absolute accuracy.
 */
export function precision(x: MpReal): number {
  if (x.s !== 0) return x.p;
  return x.e < 0 ? nbits2prec(-x.e) : LOWDEFAULTPREC;
}

/** PARI `real_0_bit(e)` (`level1.h:435`) */
export function real_0_bit(e: number): MpReal {
  return { s: 0, e, m: 0n, p: 0 };
}
/** PARI `real_0(prec)` (`level1.h:437`) */
export function real_0(prec: number): MpReal {
  return real_0_bit(-prec);
}
/** PARI `real_1(prec)` (`level1.h:441`) */
export function real_1(prec: number): MpReal {
  return { s: 1, e: 0, m: 1n << BigInt(prec - 1), p: prec };
}

/**
 * Build the t_REAL closest to `s * mag * 2^k` at `p` bits (`mag > 0`).
 * PARI rounds up when the first discarded bit is set (`mp_indep.c:216-222`).
 */
function realmk(s: -1 | 1, mag: bigint, k: number, p: number): MpReal {
  const n = bitlen(mag);
  let e = k + n - 1;
  let m: bigint;
  const sh = n - p;
  if (sh > 0) {
    m = mag >> BigInt(sh);
    if ((mag >> BigInt(sh - 1)) & 1n) {
      m += 1n;
      if (m >> BigInt(p)) {
        m >>= 1n;
        e += 1;
      }
    }
  } else {
    m = mag << BigInt(-sh);
  }
  return { s, e, m, p };
}

/** exact value of `x` as `num * 2^k` */
function realExact(x: MpReal): [bigint, number] {
  if (x.s === 0) return [0n, 0];
  return [x.s < 0 ? -x.m : x.m, x.e + 1 - x.p];
}

/**
 * absolute ulp exponent: the weight of the last mantissa bit.  For a real
 * zero PARI's `addrr_sign` (`add.c:200-215`) uses `expo + 1`, which is the
 * same quantity.
 */
function realulp(x: MpReal): number {
  return x.s === 0 ? x.e + 1 : x.e + 1 - x.p;
}

/** round the exact value `v * 2^k` at the ulp `u` (PARI's addition rule) */
function roundAtUlp(v: bigint, k: number, u: number): MpReal {
  if (v === 0n) return real_0_bit(u);
  const s: -1 | 1 = v < 0n ? -1 : 1;
  const mag = v < 0n ? -v : v;
  const e = k + bitlen(mag) - 1;
  const p = e + 1 - u;
  if (p <= 0) return real_0_bit(u);
  return realmk(s, mag, k, p);
}

/** PARI `negr` */
export function negr(x: MpReal): MpReal {
  return x.s === 0 ? x : { s: -x.s as -1 | 1, e: x.e, m: x.m, p: x.p };
}
/** PARI `absr` / `mpabs` */
export function absr(x: MpReal): MpReal {
  return x.s < 0 ? { s: 1, e: x.e, m: x.m, p: x.p } : x;
}
/** PARI `shiftr(x,n)` = `x * 2^n` (exact) */
export function shiftr(x: MpReal, n: number): MpReal {
  return { s: x.s, e: x.e + n, m: x.m, p: x.p };
}
/** PARI `setexpo(x,e)` (returns a copy: our t_REALs are immutable) */
export function setexpo(x: MpReal, e: number): MpReal {
  return { s: x.s, e, m: x.m, p: x.p };
}
/** PARI `gequal1` for a t_REAL: exactly 1 */
export function gequal1(x: MpReal): boolean {
  return x.s === 1 && x.e === 0 && x.m === 1n << BigInt(x.p - 1);
}

/** PARI `itor(x, prec)` (`mp.c`), rounding an integer to `prec` bits */
export function itor(x: bigint, prec: number): MpReal {
  if (x === 0n) return real_0(prec);
  return realmk(x < 0n ? -1 : 1, x < 0n ? -x : x, 0, prec);
}
/** PARI `rtor(x, prec)`: change the precision of a t_REAL */
export function rtor(x: MpReal, prec: number): MpReal {
  if (x.s === 0) return x;
  if (prec === x.p) return x;
  return realmk(x.s, x.m, x.e + 1 - x.p, prec);
}
/** PARI `truncr(x)`: the integer part, towards zero */
export function truncr(x: MpReal): bigint {
  if (x.s === 0 || x.e < 0) return 0n;
  const k = x.e + 1 - x.p;
  const v = k >= 0 ? x.m << BigInt(k) : x.m >> BigInt(-k);
  return x.s < 0 ? -v : v;
}
/**
 * PARI `gcvtoi(x, &e)` for a t_REAL (`gen3.c:2669-2681`): the integer part,
 * together with the exponent of the discarded part (or of the ulp when the
 * t_REAL does not even determine the integer part).
 */
export function gcvtoi(x: MpReal): [bigint, number] {
  if (x.s === 0) return [0n, x.e];
  if (x.e < 0) return [0n, x.e];
  let e1 = x.e - x.p + 1;
  const y = truncr(x);
  if (e1 <= 0) {
    /* e1 = expo(subri(x,y)): exponent of the fractional part */
    const [v, k] = realExact(x);
    const frac = v - (y << BigInt(-k)); /* k <= 0 here */
    e1 = frac === 0n ? -(1 << 30) : k + bitlen(frac < 0n ? -frac : frac) - 1;
  }
  return [y, e1];
}

/** PARI `mul0r` (`mp_indep.c:157-163`) */
function mul0r(x: MpReal): MpReal {
  const l = x.p;
  const e = l > 0 ? x.e - l : x.e < 0 ? 2 * x.e : 0;
  return real_0_bit(e);
}

/** PARI `addrr` (`add.c:181-330`): round the exact sum at `max(ulp x, ulp y)` */
export function addrr(x: MpReal, y: MpReal): MpReal {
  if (x.s === 0 && y.s === 0) return real_0_bit(Math.max(x.e, y.e));
  const [vx, kx] = realExact(x);
  const [vy, ky] = realExact(y);
  const k = Math.min(kx, ky);
  const v = (vx << BigInt(kx - k)) + (vy << BigInt(ky - k));
  return roundAtUlp(v, k, Math.max(realulp(x), realulp(y)));
}
/** PARI `subrr` */
export function subrr(x: MpReal, y: MpReal): MpReal {
  return addrr(x, negr(y));
}
/** PARI `addir` (`add.c:118-143`): exact integer + t_REAL */
export function addir(x: bigint, y: MpReal): MpReal {
  if (x === 0n) return y;
  const [vy, ky] = realExact(y);
  const k = Math.min(0, ky);
  const v = (x << BigInt(-k)) + (vy << BigInt(ky - k));
  return roundAtUlp(v, k, realulp(y));
}
/** PARI `subir(x,y) = x - y` */
export function subir(x: bigint, y: MpReal): MpReal {
  return addir(x, negr(y));
}
/** PARI `addsr` / `subrs` (`add.c:146-179`) */
export function addrs(x: MpReal, n: number): MpReal {
  return addir(BigInt(n), x);
}
export function subrs(x: MpReal, n: number): MpReal {
  return addir(BigInt(-n), x);
}

/** PARI `mulrr` (`mp_indep.c:391-405`): result precision `min(px, py)` */
export function mulrr(x: MpReal, y: MpReal): MpReal {
  if (x === y) return sqrr(x);
  if (x.s === 0 || y.s === 0) return real_0_bit(x.e + y.e);
  const p = Math.min(x.p, y.p);
  const s: -1 | 1 = x.s === y.s ? 1 : -1;
  return realmk(s, x.m * y.m, x.e + 1 - x.p + (y.e + 1 - y.p), p);
}
/** PARI `sqrr` (`mp_indep.c:409-418`): result precision `px` */
export function sqrr(x: MpReal): MpReal {
  if (x.s === 0) return real_0_bit(2 * x.e);
  return realmk(1, x.m * x.m, 2 * (x.e + 1 - x.p), x.p);
}
/** PARI `mulir` (`mp_indep.c:421-450`): result precision `py` */
export function mulir(x: bigint, y: MpReal): MpReal {
  if (x === 0n) return mul0r(y);
  if (y.s === 0) return real_0_bit(expi(x) + y.e);
  const s: -1 | 1 = (x < 0n ? -1 : 1) === y.s ? 1 : -1;
  return realmk(s, (x < 0n ? -x : x) * y.m, y.e + 1 - y.p, y.p);
}
/** PARI `mulri` (`mp_indep.c`): result precision `px` */
export function mulri(x: MpReal, y: bigint): MpReal {
  if (y === 0n) return mul0r(x);
  if (x.s === 0) return real_0_bit(expi(y) + x.e);
  const s: -1 | 1 = (y < 0n ? -1 : 1) === x.s ? 1 : -1;
  return realmk(s, (y < 0n ? -y : y) * x.m, x.e + 1 - x.p, x.p);
}
/** PARI `mulsr` (`mp_indep.c:171-188`) */
export function mulsr(n: number, y: MpReal): MpReal {
  return mulir(BigInt(n), y);
}
/** PARI `mulrs` */
export function mulrs(x: MpReal, n: number): MpReal {
  return mulri(x, BigInt(n));
}

/**
 * Divide the exact magnitudes `(na * 2^ka) / (nb * 2^kb)` and round to `p`
 * bits.  A sticky bit makes the rounding exactly the one `realmk` would do on
 * the exact quotient.
 */
function divmag(s: -1 | 1, na: bigint, ka: number, nb: bigint, kb: number, p: number): MpReal {
  /* enough shift for the quotient to carry at least p + 64 bits */
  const g = p + 64 + Math.max(0, bitlen(nb) - bitlen(na));
  const num = na << BigInt(g);
  let q = num / nb;
  if (num % nb !== 0n) q |= 1n;
  return realmk(s, q, ka - kb - g, p);
}

/** PARI `divrr` (`mp.c:635-745`): result precision `min(px, py)` */
export function divrr(x: MpReal, y: MpReal): MpReal {
  if (y.s === 0) throw new PariInvError('divrr');
  if (x.s === 0) return real_0_bit(x.e - y.e);
  const p = Math.min(x.p, y.p);
  const s: -1 | 1 = x.s === y.s ? 1 : -1;
  return divmag(s, x.m, x.e + 1 - x.p, y.m, y.e + 1 - y.p, p);
}
/** PARI `divir` (`mp_indep.c:555-570`): result precision `py` */
export function divir(x: bigint, y: MpReal): MpReal {
  if (y.s === 0) throw new PariInvError('divir');
  if (x === 0n) return real_0_bit(-y.p - y.e);
  const s: -1 | 1 = (x < 0n ? -1 : 1) === y.s ? 1 : -1;
  return divmag(s, x < 0n ? -x : x, 0, y.m, y.e + 1 - y.p, y.p);
}
/** PARI `divri` (`mp.c:748-765`): result precision `px` */
export function divri(x: MpReal, y: bigint): MpReal {
  if (y === 0n) throw new PariInvError('divri');
  if (x.s === 0) return real_0_bit(x.e - expi(y));
  const s: -1 | 1 = (y < 0n ? -1 : 1) === x.s ? 1 : -1;
  return divmag(s, x.m, x.e + 1 - x.p, y < 0n ? -y : y, 0, x.p);
}
/** PARI `divru` (`mp_indep.c:688`) */
export function divru(x: MpReal, n: number): MpReal {
  return divri(x, BigInt(n));
}

/**
 * PARI `sqrtr_abs` (`mp.c:2063-2160`).  PARI uses a Newton iteration on the
 * mantissa; we take the exact integer square root of the scaled mantissa and
 * round (with a sticky bit), which is the correctly rounded result.
 */
export function sqrtr_abs(x: MpReal): MpReal {
  if (x.s === 0) return real_0_bit(x.e >> 1);
  const p = x.p;
  let k = x.e + 1 - p;
  let mag = x.m;
  /* want (mag << g) to have >= 2p+4 bits and (k-g) even */
  let g = 2 * p + 4;
  if ((k - g) % 2 !== 0) g += 1;
  mag <<= BigInt(g);
  k -= g;
  let r = sqrti(mag);
  if (r * r !== mag) r |= 1n;
  return realmk(1, r, k / 2, p);
}
/** PARI `sqrtr` (`trans1.c`): errors on a negative argument */
export function sqrtr(x: MpReal): MpReal {
  if (x.s < 0) throw new PariDomainError('sqrtr', 'argument', '<', '0');
  return sqrtr_abs(x);
}

/**
 * `atanh(1/q) * 2^N`, rounded down, by the defining series
 * `sum_{k>=0} q^-(2k+1)/(2k+1)`.  PARI evaluates the same series by binary
 * splitting (`atanhuu`, `trans1.c`); only the speed differs.
 */
function atanhuu_scaled(q: bigint, N: number): bigint {
  const one = 1n << BigInt(N);
  const q2 = q * q;
  let S = 0n;
  let qp = q;
  let k = 1n;
  while (qp <= one) {
    S += one / (k * qp);
    qp *= q2;
    k += 2n;
  }
  return S;
}

let log2Cache: MpReal | null = null;
/**
 * PARI `mplog2` / `constlog2` / `log2_split` (`trans1.c:2841-2868`):
 * `log 2 = 18 atanh(1/26) - 2 atanh(1/4801) + 8 atanh(1/8749)`.
 */
export function mplog2(prec: number): MpReal {
  if (log2Cache === null || log2Cache.p < prec) {
    const N = nbits2prec(prec + 128) + 64;
    const L =
      18n * atanhuu_scaled(26n, N) - 2n * atanhuu_scaled(4801n, N) + 8n * atanhuu_scaled(8749n, N);
    log2Cache = realmk(1, L, -N, nbits2prec(prec + 128));
  }
  return rtor(log2Cache, prec);
}

/** approximate `log2 |x|` as a double (PARI `dbllog2r`) */
function dbllog2r(x: MpReal): number {
  if (x.s === 0) return -1e30;
  const sh = x.p - 53;
  const top = sh > 0 ? Number(x.m >> BigInt(sh)) : Number(x.m) * 2 ** -sh;
  return x.e + Math.log2(top / 2 ** 52);
}

/**
 * PARI `logr_aux` (`trans1.c:2892-2925`): `log(x)/2` where
 * `y = (x-1)/(x+1)` is close to 0, via `y * (1 + y^2/3 + y^4/5 + ...)`.
 * PARI raises the working precision as the loop advances; we run the whole
 * Horner recurrence at the input precision, which is at least as accurate.
 */
function logr_aux(y: MpReal): MpReal {
  const L = y.p;
  const d = -2 * dbllog2r(y);
  let k = Math.floor(2 * (L / d));
  k |= 1;
  if (k >= 3) {
    const y2 = sqrr(y);
    let S = divru(real_1(L), k);
    let T = S;
    for (k -= 2; ; k -= 2) {
      T = mulrr(S, y2);
      if (k === 1) break;
      S = addrr(divru(real_1(L), k), T);
    }
    return mulrr(y, addrs(T, 1));
  }
  return y;
}

/**
 * PARI `logr_abs(X)` (`trans1.c:2926-2985`): `log |X|`.
 *
 * PARI's `logagmr_abs` branch (taken when `realprec(X) > LOGAGM_LIMIT`, a
 * tuning constant far above the precisions this module works at) is not
 * ported: the series path below computes the same value.
 */
export function logr_abs(X: MpReal): MpReal {
  if (X.s === 0) throw new PariDomainError('logr_abs', 'argument', '=', '0');
  const p = X.p;
  let EX = X.e;
  /* choose the smaller of x-1 and 1-x/2 (`trans1.c:2937-2951`) */
  const u = X.m >> BigInt(p - 64);
  let D: bigint;
  if (u > 12297829382473034410n) {
    /* (~0UL/3)*2: x > 4/3, use 1 - x/2 */
    EX++;
    D = (1n << BigInt(p)) - 1n - X.m;
  } else {
    D = X.m - (1n << BigInt(p - 1));
  }
  if (D === 0n) return EX ? mulsr(EX, mplog2(p)) : real_0(p);
  const a = p - bitlen(D); /* ~ -log2 |1-x| */
  let L = p + 64; /* EXTRAPRECWORD */
  const b = L - 64 * Math.floor(a / 64);
  const dd = -a / 2;
  let m = Math.floor(dd + Math.sqrt(dd * dd + b / 6));
  if (m > b - a) m = b - a;
  if (m < 0.2 * a) m = 0;
  else L += nbits2prec(m);
  let x = shiftr(rtor(absr(X), L), -EX); /* 2/3 < x < 4/3 */
  for (let i = 1; i <= m; i++) x = sqrtr_abs(x);
  let y = divrr(subrs(x, 1), addrs(x, 1));
  y = logr_aux(y);
  y = shiftr(y, m + 1);
  if (EX) y = addrr(y, mulsr(EX, mplog2(p + 64)));
  const outp = EX ? p : Math.max(64, p - 64 * Math.floor(a / 64));
  return rtor(y, outp);
}

/**
 * The exact value of `x` as a fraction `[num, den]`, `den > 0` a power of two.
 * Not an upstream function: an exactness hook for the tests, which compare our
 * t_REALs with decimal literals by exact rational arithmetic.
 */
export function mpreal_to_frac(x: MpReal): [bigint, bigint] {
  const [v, k] = realExact(x);
  return k >= 0 ? [v << BigInt(k), 1n] : [v, 1n << BigInt(-k)];
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

/**
 * PARI's "extended t_QFB": the `t_VEC [t_QFB, t_REAL]` of an *indefinite*
 * form together with Shanks' logarithmic distance (`check_qfbext`,
 * `Qfb.c:111-123`). Every public entry point below accepts it wherever PARI
 * does, and returns one when given one.
 */
export interface QfbExt {
  readonly Q: Qfb;
  /** Shanks' distance, a t_REAL */
  readonly d: MpReal;
}

/** a `t_QFB` or an extended `t_QFB` */
export type QfbLike = Qfb | QfbExt;

/** is this PARI's `t_VEC [t_QFB, t_REAL]` rather than a bare `t_QFB`? */
export function is_qfbext(x: QfbLike): x is QfbExt {
  return (x as QfbExt).Q !== undefined;
}

/** PARI `check_qfbext(fun, x)` (`Qfb.c:111-123`): returns the underlying form */
function check_qfbext(fun: string, x: QfbLike): Qfb {
  if (!is_qfbext(x)) return x;
  if (x.Q.D < 0n) throw new PariTypeError(fun, 'definite form with a distance');
  return x.Q;
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
function REDBU(a: bigint, b: bigint, c: bigint, u1: bigint, u2: bigint): [bigint, bigint, bigint] {
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
function rho_get_BC(a: bigint, b: bigint, c: bigint, isqrtD: bigint): [bigint, bigint] {
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
    [A[0]![0]! * B[0]![0]! + A[0]![1]! * B[1]![0]!, A[0]![0]! * B[0]![1]! + A[0]![1]! * B[1]![1]!],
    [A[1]![0]! * B[0]![0]! + A[1]![1]! * B[1]![0]!, A[1]![0]! * B[0]![1]! + A[1]![1]! * B[1]![1]!],
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

/* ------------------------------------------------------------------ */
/* qfr3 / qfr5: Shanks' distance (Qfb.c:396-560, 1470-1640)            */
/* ------------------------------------------------------------------ */

/*
 * PARI's own comment (`Qfb.c:396-406`): "qfr3 / qfr5 routines take a container
 * of t_INTs as argument, at least 3 (resp. 5) components. A qfr3 [a,b,c]
 * contains the form coeffs, in a qfr5 [a,b,c,e,d] the t_INT e is a binary
 * exponent, d a t_REAL, coding the distance in multiplicative form: the true
 * distance is obtained from qfr5_dist."
 */

/** PARI's `qfr3` container `[a,b,c]` */
export interface Qfr3 {
  readonly a: bigint;
  readonly b: bigint;
  readonly c: bigint;
}
/** PARI's `qfr5` container `[a,b,c, e, d]`: `e` a binary exponent, `d` a t_REAL */
export interface Qfr5 extends Qfr3 {
  readonly e: bigint;
  readonly d: MpReal;
}
/** PARI `struct qfr_data` (`paripriv.h`): the discriminant and its square roots */
export interface QfrData {
  D: bigint;
  sqrtD: MpReal | null;
  isqrtD: bigint | null;
}

/** PARI `EMAX` (`Qfb.c:409`) */
const EMAX = 22;

/** PARI `fix_expo` (`Qfb.c:410-417`) */
function fix_expo(x: Qfr5): Qfr5 {
  if (expo(x.d) >= 1 << EMAX) {
    return { a: x.a, b: x.b, c: x.c, e: x.e + 1n, d: shiftr(x.d, -(1 << EMAX)) };
  }
  return x;
}

/**
 * PARI `qfr5_dist(e, d, prec)` (`Qfb.c:419-430`):
 * `(1/2) log(|d| * 2^{e * 2^EMAX})`, the true Shanks distance coded by the
 * `(e, d)` pair of a qfr5.
 */
export function qfr5_dist(e: bigint, d: MpReal, prec: number): MpReal {
  let t = logr_abs(d);
  if (e !== 0n) {
    let u = mulir(e, mplog2(prec));
    u = shiftr(u, EMAX);
    t = addrr(t, u);
  }
  return shiftr(t, -1);
}

/** PARI `qfr5_rho` (`Qfb.c:452-472`): `rho`, updating the distance */
function qfr5_rho(x: Qfr5, S: QfrData): Qfr5 {
  const { a, b, c } = x;
  const sb = isign(b);
  const [B, C] = rho_get_BC(a, b, c, S.isqrtD!);
  let y: Qfr5 = { a: c, b: B, c: C, e: x.e, d: x.d };
  if (sb) {
    const t0 = b * b - S.D;
    /* t = (b + sqrt(D)) / (b - sqrt(D)), evaluated stably */
    const t = sb < 0 ? divir(t0, sqrr(subir(b, S.sqrtD!))) : divri(sqrr(addir(b, S.sqrtD!)), t0);
    y = { a: y.a, b: y.b, c: y.c, e: y.e, d: mulrr(t, y.d) };
    y = fix_expo(y);
  } else {
    y = { a: y.a, b: y.b, c: y.c, e: y.e, d: negr(y.d) };
  }
  return y;
}

/** PARI `qfr5_red` (`Qfb.c:517-531`) */
function qfr5_red(x: Qfr5, S: QfrData): Qfr5 {
  while (!ab_isreduced(x.a, x.b, S.isqrtD!)) x = qfr5_rho(x, S);
  return x;
}

/** PARI `qfr3_red` (`Qfb.c:533-550`), qfr3 container flavour */
function qfr3_red3(x: Qfr3, S: QfrData): Qfr3 {
  const [a, b, c] = qfr3_red(x.a, x.b, x.c, S.isqrtD!);
  return { a, b, c };
}

/** PARI `qfr3_rho` (`Qfb.c:444-450`), qfr3 container flavour */
function qfr3_rho3(x: Qfr3, S: QfrData): Qfr3 {
  const [a, b, c] = qfr3_rho(x.a, x.b, x.c, S.isqrtD!);
  return { a, b, c };
}

/** PARI `qfr_to_qfr5` (`Qfb.c:474-477`) */
function qfr_to_qfr5(x: Qfb, prec: number): Qfr5 {
  return { a: x.a, b: x.b, c: x.c, e: 0n, d: real_1(prec) };
}

/**
 * PARI `qfr5_to_qfr(x, D, d0)` (`Qfb.c:479-503`): fold the multiplicative
 * distance `(e, d)` accumulated by the qfr5 routines into the initial
 * (logarithmic) distance `d0`, and rebuild a t_QFB.
 *
 * DEVIATION. Upstream reads `mplog2(lg(d0))` (`Qfb.c:495`), i.e. it passes the
 * *word length* of `d0` where a bit precision is expected -- a call site that
 * was missed when PARI 2.16 changed `prec` from words to bits, and which would
 * add `n * log 2` at 4 bits of accuracy. We pass `precision(d0)` (`gen3.c:142`). The branch is
 * reachable only when `fix_expo` has fired, i.e. `|d| >= 2^(2^22)`.
 */
function qfr5_to_qfr(x: Qfr3 | Qfr5, D: bigint, d0: MpReal | null): QfbLike {
  let dd = d0;
  if (dd !== null) {
    const X = x as Qfr5;
    let n = X.e;
    let d = absr(X.d);
    if (n !== 0n) {
      n = (n << BigInt(EMAX)) + BigInt(expo(d));
      d = setexpo(d, 0);
      d = logr_abs(d);
      if (n !== 0n) d = addrr(d, mulir(n, mplog2(precision(dd))));
      d = shiftr(d, -1);
      dd = addrr(dd, d);
    } else if (!gequal1(d)) {
      /* avoid loss of precision */
      d = logr_abs(d);
      d = shiftr(d, -1);
      dd = addrr(dd, d);
    }
  }
  const q = mkqfb(x.a, x.b, x.c, D);
  return dd !== null ? { Q: q, d: dd } : q;
}

/** PARI `qfr_data_init(D, prec, S)` (`Qfb.c:552-558`) */
export function qfr_data_init(D: bigint, prec: number): QfrData {
  const sqrtD = sqrtr(itor(D, prec));
  return { D, sqrtD, isqrtD: truncr(sqrtD) };
}

/** PARI `qfr5_init` (`Qfb.c:560-583`) */
function qfr5_init(x: Qfb, d: MpReal, S: QfrData): Qfr5 {
  let prec = realprec(d);
  let l = -expo(d);
  if (l < 64) l = 64;
  prec = Math.max(prec, nbits2prec(l));
  S.D = x.D;
  const y = qfr_to_qfr5(x, prec);
  if (S.sqrtD === null) S.sqrtD = sqrtr(itor(S.D, prec));
  if (S.isqrtD === null) {
    const [n, e] = gcvtoi(S.sqrtD);
    S.isqrtD = e > -2 ? sqrti(S.D) : n;
  }
  return y;
}

/** PARI `qfr3_init` (`Qfb.c:585-591`) */
function qfr3_init(x: Qfb, S: QfrData): Qfr3 {
  S.D = x.D;
  if (S.isqrtD === null) S.isqrtD = sqrti(S.D);
  return { a: x.a, b: x.b, c: x.c };
}

/**
 * PARI `qfr_1_fill` (`Qfb.c:1252-1261`), the principal form of `S.D`.
 *
 * DEVIATION. Upstream writes `y2 = subiu(y,1)` where `y` is the *container
 * being filled*, not the integer `y2` (`Qfb.c:1257`) -- a typo that would
 * produce garbage. We use `y2 - 1`, which is what `qfr_1_by_disc`
 * (`Qfb.c:1215-1232`) computes for the same discriminant. Unreachable from
 * the public entry points, which handle `n = 0` before calling `qfr*_pow`.
 */
function qfr_1_fill(S: QfrData): Qfr3 {
  let y2 = S.isqrtD!;
  if (mod2(S.D) !== mod2(y2)) y2 = y2 - 1n;
  return { a: 1n, b: y2, c: shifti(y2 * y2 - S.D, -2) };
}
/** PARI `qfr5_1` (`Qfb.c:1263-1270`) */
function qfr5_1(S: QfrData, prec: number): Qfr5 {
  const y = qfr_1_fill(S);
  return { a: y.a, b: y.b, c: y.c, e: 0n, d: real_1(prec) };
}
/** PARI `qfr3_1` (`Qfb.c:1272-1277`) */
function qfr3_1(S: QfrData): Qfr3 {
  return qfr_1_fill(S);
}

/**
 * `qfb_comp` on qfr containers. PARI dispatches to `qfb_sqr` on *pointer*
 * identity (`Qfb.c:1043`), which is what the `x === y` test reproduces: it is
 * how `qfr5_powraw`/`qfr5_pow` square (`x = qfr5_compraw(x,x)`).
 */
function qfr_compraw3(x: Qfr3, y: Qfr3, D: bigint): Qfr3 {
  const X = mkqfb(x.a, x.b, x.c, D);
  if (x === y) return qfb_sqr(X);
  return qfb_comp(X, mkqfb(y.a, y.b, y.c, D));
}

/** PARI `qfr5_compraw` (`Qfb.c:1470-1486`) */
function qfr5_compraw(x: Qfr5, y: Qfr5, D: bigint): Qfr5 {
  const z = qfr_compraw3(x, y, D);
  const e = x === y ? shifti(x.e, 1) : x.e + y.e;
  const d = x === y ? sqrr(x.d) : mulrr(x.d, y.d);
  return fix_expo({ a: z.a, b: z.b, c: z.c, e, d });
}
/** PARI `qfr5_comp` (`Qfb.c:1487-1489`) */
function qfr5_comp(x: Qfr5, y: Qfr5, S: QfrData): Qfr5 {
  return qfr5_red(qfr5_compraw(x, y, S.D), S);
}
/** PARI `qfr3_compraw` (`Qfb.c:1490-1496`) */
function qfr3_compraw(x: Qfr3, y: Qfr3, D: bigint): Qfr3 {
  return qfr_compraw3(x, y, D);
}
/** PARI `qfr3_comp` (`Qfb.c:1497-1499`) */
function qfr3_comp(x: Qfr3, y: Qfr3, S: QfrData): Qfr3 {
  return qfr3_red3(qfr3_compraw(x, y, S.D), S);
}

/** PARI `qfr5_powraw` (`Qfb.c:1501-1513`), `m > 0` */
function qfr5_powraw(x: Qfr5, m: bigint, D: bigint): Qfr5 {
  let y: Qfr5 | null = null;
  for (; m; m >>= 1n) {
    if (m & 1n) y = y ? qfr5_compraw(y, x, D) : x;
    if (m === 1n) break;
    x = qfr5_compraw(x, x, D);
  }
  return y!;
}

/**
 * PARI `qfr5_pow` (`Qfb.c:1515-1534`), assuming `n > 0` (see the note on
 * {@link qfrpow}).
 *
 * DEVIATION. Upstream loops over the *machine words* of `n`, least significant
 * first, with `if (m == 1 && i == 2) break;`, so a word with leading zero bits
 * ends its inner loop early and the squarings it owed are never done; and it
 * reads the word into a *signed* `long m`, so `m >>= 1` on a word whose top
 * bit is set is an arithmetic shift that never reaches 0. Both were confirmed
 * on the live PARI 2.15.4:
 *
 *     f = Qfb(-2020,879,1142);
 *     qfbpow(f, 2^64+1) == qfbpow(f, 3)   \\ both Qfb(418, 3025, -508)
 *     qfbpow(f, 2^63)                     \\ runs until the stack overflows
 *
 * and the 2.18.1 source we ported is unchanged. We use the plain right-to-left
 * binary chain over the whole exponent -- which is exactly what the word loop
 * computes for every single-word `n`, and what the `t_QFB` powering in this
 * file has always done. `qfb.test.ts` pins the consequence: our `f^(10^20)`
 * carries a distance congruent to its form's cycle distance mod the regulator,
 * PARI 2.15.4's does not (residual 1057.8, with R = 2641.55).
 */
function qfr5_pow(x: Qfr5, n: bigint, S: QfrData): Qfr5 {
  if (n === 0n) return qfr5_1(S, precision(x.d));
  let y: Qfr5 | null = null;
  let m = n < 0n ? -n : n;
  for (;;) {
    if (m & 1n) y = y === null ? x : qfr5_comp(y, x, S);
    m >>= 1n;
    if (m === 0n) break;
    x = qfr5_comp(x, x, S);
  }
  return y!;
}

/** PARI `qfr3_powraw` (`Qfb.c:1536-1547`), `m > 0` */
function qfr3_powraw(x: Qfr3, m: bigint, D: bigint): Qfr3 {
  let y: Qfr3 | null = null;
  for (; m; m >>= 1n) {
    if (m & 1n) y = y ? qfr3_compraw(y, x, D) : x;
    if (m === 1n) break;
    x = qfr3_compraw(x, x, D);
  }
  return y!;
}

/** PARI `qfr3_pow` (`Qfb.c:1549-1568`); same deviation as {@link qfr5_pow} */
function qfr3_pow(x: Qfr3, n: bigint, S: QfrData): Qfr3 {
  if (n === 0n) return qfr3_1(S);
  let y: Qfr3 | null = null;
  let m = n < 0n ? -n : n;
  for (;;) {
    if (m & 1n) y = y === null ? x : qfr3_comp(y, x, S);
    m >>= 1n;
    if (m === 0n) break;
    x = qfr3_comp(x, x, S);
  }
  return y!;
}

/** PARI `qfrinvraw` (`Qfb.c:1570-1575`) */
function qfrinvraw(x: QfbLike): QfbLike {
  if (is_qfbext(x)) return { Q: qfbinv(x.Q), d: negr(x.d) };
  return qfbinv(x);
}

/**
 * PARI `qfrpowraw(x, n)` (`Qfb.c:1577-1602`) for an extended t_QFB.
 *
 * DEVIATION. Upstream negates `n` in place before computing the final distance
 * (`if (n < 0) { x = qfb_inv(x); n = -n; } ... qfr5_to_qfr(x, S.D, mulrs(d0,n))`,
 * `Qfb.c:1594-1599`), so `qfbpowraw([f,d], -k)` gets the distance `+k*d`
 * instead of `-k*d` (verified on PARI 2.15.4: `qfbpowraw(f,-3)` and
 * `qfbpowraw(f,3)` both report `12.760258257204765447179324751139353799`).
 * We use the signed exponent, so that the distance of `x^-k` is `-k` times the
 * distance of `x`, as it must be. The *form* is upstream's.
 */
function qfrpowraw_ext(X: QfbExt, n: bigint): QfbLike {
  const d0 = X.d;
  let x = X.Q;
  if (n === 0n) return { Q: qfb_1(x), d: real_0(precision(d0)) };
  const n0 = n;
  if (n < 0n) {
    x = qfbinv(x);
    n = -n;
  }
  const S: QfrData = { D: 0n, sqrtD: null, isqrtD: null };
  let y = qfr5_init(x, d0, S);
  if (n !== 1n) y = qfr5_powraw(y, n, S.D);
  return qfr5_to_qfr(y, S.D, mulri(d0, n0));
}

/**
 * PARI `qfrpow(x, n)` (`Qfb.c:1604-1628`) for an extended t_QFB.
 *
 * DEVIATION. Upstream inverts the form for `n < 0` and then hands the *signed*
 * `n` to `qfr5_pow`, which inverts a second time (`Qfb.c:1521`); PARI therefore
 * returns `x^|n|` for `n <= -2` (verified on PARI 2.15.4: `qfbpow(f,-6)` and
 * `qfbpow(f,6)` return the same form, with opposite distances -- and the same
 * happens on the plain `t_QFB` path through `qfr3_pow`). We invert once and
 * raise to `|n|`, so `x^-n` really is the inverse of `x^n`.
 */
function qfrpow_ext(X: QfbExt, n: bigint): QfbLike {
  const d0 = X.d;
  let x = X.Q;
  const s = isign(n);
  if (s === 0) return { Q: qfb_1(x), d: real_0(precision(d0)) };
  if (s < 0) x = qfbinv(x);
  const S: QfrData = { D: 0n, sqrtD: null, isqrtD: null };
  let y = qfr5_init(x, d0, S);
  y = n === 1n || n === -1n ? qfr5_red(y, S) : qfr5_pow(y, n < 0n ? -n : n, S);
  return qfr5_to_qfr(y, S.D, mulri(d0, n));
}

/**
 * PARI `qfr_red_basecase_i(x, flag, isqrtD, sqrtD)` (`Qfb.c:594-611`), the
 * only reduction path PARI uses for an extended t_QFB (`qfr_red_i`,
 * `Qfb.c:911-913`, sends every `t_VEC` here).
 */
function qfr_red_basecase_i(
  x: QfbLike,
  flag: number,
  isqrtD: bigint | null,
  sqrtD: MpReal | null
): QfbLike {
  let d: MpReal | null = null;
  let q: Qfb;
  if (is_qfbext(x)) {
    d = x.d;
    q = x.Q;
  } else {
    q = x;
    flag |= qf_NOD;
  }
  const S: QfrData = { D: 0n, sqrtD, isqrtD };
  let y: Qfr3 | Qfr5;
  if (flag & qf_NOD) {
    y = qfr3_init(q, S);
    y = flag & qf_STEP ? qfr3_rho3(y, S) : qfr3_red3(y, S);
  } else {
    const y5 = qfr5_init(q, d!, S);
    y = flag & qf_STEP ? qfr5_rho(y5, S) : qfr5_red(y5, S);
  }
  return qfr5_to_qfr(y, q.D, d);
}

const qf_NOD = 2;
const qf_STEP = 1;

/**
 * PARI `qfbred0(x, flag, isqrtD, sqrtD)` (`Qfb.c:991-1002`).
 *
 * @param flag bit 0 (`qf_STEP`): perform a single reduction step (`rho`);
 *             bit 1 (`qf_NOD`): ignore the distance. `qf_NOD` is forced on for
 *             a bare `t_QFB` and forced off for an extended one, exactly as in
 *             `Qfb.c:997-998`.
 * @param sqrtD optional `sqrt(D)` as a t_REAL, only used for an extended form.
 */
export function qfbred(q: Qfb, flag?: number, isqrtD?: bigint | null, sqrtD?: MpReal | null): Qfb;
export function qfbred(
  q: QfbLike,
  flag?: number,
  isqrtD?: bigint | null,
  sqrtD?: MpReal | null
): QfbLike;
export function qfbred(
  q: QfbLike,
  flag = 0,
  isqrtD: bigint | null = null,
  sqrtD: MpReal | null = null
): QfbLike {
  if (flag < 0 || flag > 3) throw new PariFlagError('qfbred');
  const qq = check_qfbext('qfbred', q);
  if (qfb_is_qfi(qq)) {
    const step = (flag & qf_STEP) !== 0;
    return step ? qfi_rho(qq) : qfi_red(qq);
  }
  if (!is_qfbext(q)) return qfr_red_i(q, (flag & qf_STEP) !== 0, isqrtD);
  return qfr_red_basecase_i(q, flag & ~qf_NOD, isqrtD, sqrtD);
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

/**
 * PARI `qfrcomp0(x, y, raw)` (`Qfb.c:1124-1137`): composition of two forms of
 * the same (positive) discriminant, adding their distances.
 */
function qfrcomp0(x: QfbLike, y: QfbLike, raw: boolean): QfbLike {
  const dx = is_qfbext(x) ? x.d : null;
  const dy = is_qfbext(y) ? y.d : null;
  const X = is_qfbext(x) ? x.Q : x;
  const Y = is_qfbext(y) ? y.Q : y;
  const z = qfb_comp(X, Y);
  let w: QfbLike = z;
  if (dx) w = { Q: z, d: dy ? addrr(dx, dy) : dx };
  else if (dy) w = { Q: z, d: dy };
  if (raw) return w;
  return qfbred(w);
}

/** PARI `qfrsqr0(x, raw)` (`Qfb.c:1194-1204`) */
function qfrsqr0(x: QfbLike, raw: boolean): QfbLike {
  const dx = is_qfbext(x) ? x.d : null;
  const X = is_qfbext(x) ? x.Q : x;
  const z = qfb_sqr(X);
  const w: QfbLike = dx ? { Q: z, d: shiftr(dx, 1) } : z;
  if (raw) return w;
  return qfbred(w);
}

/** PARI `qfbcompraw` (`Qfb.c:1165-1181`) */
export function qfbcompraw(x: Qfb, y: Qfb): Qfb;
export function qfbcompraw(x: QfbLike, y: QfbLike): QfbLike;
export function qfbcompraw(x: QfbLike, y: QfbLike): QfbLike {
  const qx = check_qfbext('qfbcompraw', x);
  const qy = check_qfbext('qfbcompraw', y);
  if (qx.D !== qy.D) {
    const z = qfb_comp_gen(qx, qy);
    if (is_qfbext(x) || is_qfbext(y))
      throw new NotImplementedError("Shanks's distance in general composition");
    if (!z) throw new PariDomainError('qfbcompraw', 'discriminants', '!=', 'equal');
    return z;
  }
  if (qfb_is_qfi(qx)) return qfb_comp(qx, qy);
  return qfrcomp0(x, y, true);
}

/** PARI `qfbcomp` (`Qfb.c:1145-1161`) */
export function qfbcomp(x: Qfb, y: Qfb): Qfb;
export function qfbcomp(x: QfbLike, y: QfbLike): QfbLike;
export function qfbcomp(x: QfbLike, y: QfbLike): QfbLike {
  const qx = check_qfbext('qfbcomp', x);
  const qy = check_qfbext('qfbcomp', y);
  if (qx.D !== qy.D) {
    const z = qfb_comp_gen(qx, qy);
    if (is_qfbext(x) || is_qfbext(y))
      throw new NotImplementedError("Shanks's distance in general composition");
    if (!z) throw new PariDomainError('qfbcomp', 'discriminants', '!=', 'equal');
    return qfbred(z);
  }
  if (qfb_is_qfi(qx)) return qfbred(qfb_comp(qx, qy));
  return qfrcomp0(x, y, false);
}

/** PARI `qfbsqr` (`Qfb.c:1208-1213`) */
export function qfbsqr(x: Qfb): Qfb;
export function qfbsqr(x: QfbLike): QfbLike;
export function qfbsqr(x: QfbLike): QfbLike {
  const qx = check_qfbext('qfbsqr', x);
  if (qfb_is_qfi(qx)) return qfbred(qfb_sqr(qx));
  return qfrsqr0(x, false);
}

/** PARI `qfbsqr` with `raw = 1` (`qfisqr0`/`qfrsqr0`, `Qfb.c:1183-1204`) */
export function qfbsqrraw(x: Qfb): Qfb;
export function qfbsqrraw(x: QfbLike): QfbLike;
export function qfbsqrraw(x: QfbLike): QfbLike {
  const qx = check_qfbext('qfbsqr', x);
  if (qfb_is_qfi(qx)) return qfb_sqr(qx);
  return qfrsqr0(x, true);
}

/**
 * PARI `qfbpowraw` (`Qfb.c:1634-1639`): `x^n` without reduction.
 *
 * NB the two branches use *different* addition chains in PARI, and raw
 * composition is not canonical, so the chain matters: definite forms go through
 * `gen_powu` (left-to-right binary, `bb_group.c:120-153`) while indefinite forms
 * go through `qfr3_powraw` (right-to-left, `Qfb.c:1508-1519`). We mirror both.
 */
export function qfbpowraw(x: Qfb, n: bigint): Qfb;
export function qfbpowraw(x: QfbLike, n: bigint): QfbLike;
export function qfbpowraw(x: QfbLike, n: bigint): QfbLike {
  const q = check_qfbext('qfbpowraw', x);
  if (is_qfbext(x)) {
    if (n === 1n) return x;
    if (n === -1n) return qfrinvraw(x);
    return qfrpowraw_ext(x, n);
  }
  if (n === 0n) return qfb_1(q);
  if (n === 1n) return q;
  if (n === -1n) return qfbinv(q);
  const base = n < 0n ? qfbinv(q) : q;
  const e = n < 0n ? -n : n;
  if (qfb_is_qfi(base)) return leftrightPow(base, e, qfb_sqr, qfb_comp);
  const y = qfr3_powraw({ a: base.a, b: base.b, c: base.c }, e, base.D);
  return mkqfb(y.a, y.b, y.c, base.D);
}

/** left-to-right binary powering, PARI `leftright_binary_powu` (`bb_group.c:120-144`) */
function leftrightPow(x: Qfb, n: bigint, sqr: (a: Qfb) => Qfb, mul: (a: Qfb, b: Qfb) => Qfb): Qfb {
  if (n === 1n) return x;
  const bits = n.toString(2);
  let y = x;
  for (let i = 1; i < bits.length; i++) {
    y = sqr(y);
    if (bits[i] === '1') y = mul(y, x);
  }
  return y;
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
export function qfbpow(x: Qfb, n: bigint): Qfb;
export function qfbpow(x: QfbLike, n: bigint): QfbLike;
export function qfbpow(x: QfbLike, n: bigint): QfbLike {
  const q = check_qfbext('qfbpow', x);
  if (is_qfbext(x)) return qfrpow_ext(x, n);
  if (n === 0n) return qfb_1(q);
  const base0 = n < 0n ? qfbinv(q) : q;
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
  const S: QfrData = { D: base0.D, sqrtD: null, isqrtD: null };
  const x3 = qfr3_init(base0, S);
  const y = e === 1n ? qfr3_red3(x3, S) : qfr3_pow(x3, e, S);
  return mkqfb(y.a, y.b, y.c, base0.D);
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
function qfisolve_normform(Qr: { Q: Qfb; U: bigint[][] }, P: Qfb): [bigint, bigint] | null {
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
