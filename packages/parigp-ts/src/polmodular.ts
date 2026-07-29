/**
 * @module parigp-ts/polmodular
 * @description Modular polynomials `Phi_L(X, Y)` — a port of PARI's
 * `reference/pari/src/basemath/polmodular.c` (Sutherland's CM/isogeny-volcano
 * algorithm, contributed to PARI in 2014).
 *
 * Why this lives in `parigp-ts`: PARI's SEA implementation reads the modular
 * equation of level `ell` from the optional `seadata` package, but when that
 * file is absent `get_modular_eqn` (`reference/pari/src/basemath/ellsea.c:112-124`)
 * falls back to
 *
 * ```c
 *   M->type = 'J';
 *   M->eq = polmodular_ZXX(ell, ell==3? 0: 5, vx, vy);
 * ```
 *
 * i.e. PARI computes the modular equation itself. So a faithful port of SEA
 * needs a faithful port of `polmodular`.
 *
 * Upstream files ported here:
 * - `reference/pari/src/basemath/polmodular.c` (all of the code below unless
 *   stated otherwise)
 * - `reference/pari/src/basemath/polclass.c:2111-2147` (`check_modinv`)
 * - `reference/pari/src/basemath/mftrace.c:2489-2508` (`hclassno6_count`)
 *
 * MATRIX LAYOUT: as everywhere else in this package (see `matkermod.ts`),
 * a `ZM` is PARI column-major: `M[j][i]` is the entry in row `i`, column `j`,
 * both 0-based. `polmodular_ZM(L, inv)[j][i]` is therefore the coefficient of
 * `X^i Y^j` in `Phi_L`, which for the symmetric `Phi_L` is the same either way.
 */

import { Fp_inv, Fp_mul, Fp_pow, Fp_red, Fp_sqrt, kronecker } from './ff.js';
import {
  FpXQ_pow,
  FpX_degree,
  FpX_divrem,
  FpX_gcd,
  FpX_mul,
  FpX_normalize,
  FpX_renormalize,
  FpX_sub,
} from './ffinit.js';
import { NotImplementedError, Z_factor, factoru, isPrime, isqrt } from './ifactor.js';
import { PariDomainError, type ZM } from './matkermod.js';
import {
  type Qfb,
  mkqfb,
  primeform,
  qfb_1,
  qfb_equal,
  qfbcomp,
  qfbinv,
  qfbpow,
  qfbred,
  qfbsqr,
} from './qfb.js';

/* ================================================================== */
/* Errors                                                             */
/* ================================================================== */

/** PARI `pari_err_BUG` (`reference/pari/src/language/err.c`). */
export class PariBugError extends Error {
  constructor(what: string) {
    super(`Bug in ${what}, please report`);
    this.name = 'PariBugError';
  }
}

/** PARI `pari_err_IMPL`. */
export class PariImplError extends Error {
  constructor(what: string) {
    super(`sorry, ${what} is not yet implemented`);
    this.name = 'PariImplError';
  }
}

/** PARI `pari_err_PRIORITY`. */
export class PariPriorityError extends Error {
  constructor(fun: string, x: string, op: string, v: number) {
    super(`incorrect priority in ${fun}: variable ${x} ${op} ${v}`);
    this.name = 'PariPriorityError';
  }
}

/** PARI `pari_err(e_ARCH, ...)`. */
export class PariArchError extends Error {
  constructor(what: string) {
    super(`sorry, ${what} not available on this system`);
    this.name = 'PariArchError';
  }
}

/* ================================================================== */
/* Class invariants (`paripriv.h:208-238`)                            */
/* ================================================================== */

export const INV_J = 0;
export const INV_F = 1;
export const INV_F2 = 2;
export const INV_F3 = 3;
export const INV_F4 = 4;
export const INV_G2 = 5;
export const INV_W2W3 = 6;
export const INV_F8 = 8;
export const INV_W3W3 = 9;
export const INV_W2W5 = 10;
export const INV_W2W7 = 14;
export const INV_W3W5 = 15;
export const INV_W3W7 = 21;
export const INV_W2W3E2 = 23;
export const INV_W2W5E2 = 24;
export const INV_W2W13 = 26;
export const INV_W2W7E2 = 27;
export const INV_W3W3E2 = 28;
export const INV_W5W7 = 35;
export const INV_W3W13 = 39;
export const INV_ATKIN3 = 103;
export const INV_ATKIN5 = 105;
export const INV_ATKIN7 = 107;
export const INV_ATKIN11 = 111;
export const INV_ATKIN13 = 113;
export const INV_ATKIN17 = 117;
export const INV_ATKIN19 = 119;
export const INV_ATKIN23 = 123;
export const INV_ATKIN29 = 129;
export const INV_ATKIN31 = 131;
export const INV_LAST = 131;

/** PARI `check_modinv` (`polclass.c:2111-2147`). */
export function check_modinv(inv: number): void {
  switch (inv) {
    case INV_J:
    case INV_F:
    case INV_F2:
    case INV_F3:
    case INV_F4:
    case INV_G2:
    case INV_W2W3:
    case INV_F8:
    case INV_W3W3:
    case INV_W2W5:
    case INV_W2W7:
    case INV_W3W5:
    case INV_W3W7:
    case INV_W2W3E2:
    case INV_W2W5E2:
    case INV_W2W13:
    case INV_W2W7E2:
    case INV_W3W3E2:
    case INV_W5W7:
    case INV_W3W13:
    case INV_ATKIN3:
    case INV_ATKIN5:
    case INV_ATKIN7:
    case INV_ATKIN11:
    case INV_ATKIN13:
    case INV_ATKIN17:
    case INV_ATKIN19:
    case INV_ATKIN23:
    case INV_ATKIN29:
    case INV_ATKIN31:
      return;
    default:
      throw new PariDomainError('polmodular', 'inv', 'invalid invariant', String(inv));
  }
}

/**
 * PARI `modinv_level` (`polmodular.c:28-63`): the square-free part of the
 * level of the modular function `inv`.
 */
export function modinv_level(inv: number): number {
  switch (inv) {
    case INV_J:
      return 1;
    case INV_G2:
    case INV_W3W3E2:
      return 3;
    case INV_F:
    case INV_F2:
    case INV_F4:
    case INV_F8:
      return 6;
    case INV_F3:
      return 2;
    case INV_W3W3:
      return 6;
    case INV_W2W7E2:
    case INV_W2W7:
      return 14;
    case INV_W3W5:
      return 15;
    case INV_W2W3E2:
    case INV_W2W3:
      return 6;
    case INV_W2W5E2:
    case INV_W2W5:
      return 30;
    case INV_W2W13:
      return 26;
    case INV_W3W7:
      return 42;
    case INV_W5W7:
      return 35;
    case INV_W3W13:
      return 39;
    case INV_ATKIN3:
    case INV_ATKIN5:
    case INV_ATKIN7:
    case INV_ATKIN11:
    case INV_ATKIN13:
    case INV_ATKIN17:
    case INV_ATKIN19:
    case INV_ATKIN23:
    case INV_ATKIN29:
    case INV_ATKIN31:
      return inv - 100;
  }
  throw new PariBugError('modinv_level');
}

/**
 * PARI `modinv_degree` (`polmodular.c:69-98`): where applicable returns
 * `N = p1*p2` (possibly `p2 = 1`) such that two j's related to the same f are
 * `N`-isogenous, and 0 otherwise. `out` receives `[p1, p2]`.
 */
export function modinv_degree(out: { p1: number; p2: number }, inv: number): number {
  switch (inv) {
    case INV_W3W5:
      out.p1 = 3;
      out.p2 = 5;
      return 15;
    case INV_W2W3E2:
    case INV_W2W3:
      out.p1 = 2;
      out.p2 = 3;
      return 6;
    case INV_W2W5E2:
    case INV_W2W5:
      out.p1 = 2;
      out.p2 = 5;
      return 10;
    case INV_W2W7E2:
    case INV_W2W7:
      out.p1 = 2;
      out.p2 = 7;
      return 14;
    case INV_W2W13:
      out.p1 = 2;
      out.p2 = 13;
      return 26;
    case INV_W3W7:
      out.p1 = 3;
      out.p2 = 7;
      return 21;
    case INV_W3W3E2:
    case INV_W3W3:
      out.p1 = 3;
      out.p2 = 3;
      return 9;
    case INV_W5W7:
      out.p1 = 5;
      out.p2 = 7;
      return 35;
    case INV_W3W13:
      out.p1 = 3;
      out.p2 = 13;
      return 39;
    case INV_ATKIN3:
    case INV_ATKIN5:
    case INV_ATKIN7:
    case INV_ATKIN11:
    case INV_ATKIN13:
    case INV_ATKIN17:
    case INV_ATKIN19:
    case INV_ATKIN23:
    case INV_ATKIN29:
    case INV_ATKIN31:
      out.p1 = inv - 100;
      out.p2 = 1;
      return inv - 100;
  }
  out.p1 = 1;
  out.p2 = 1;
  return 0;
}

/** PARI `modinv_odd_conductor` (`polmodular.c:103-112`). */
export function modinv_odd_conductor(inv: number): boolean {
  switch (inv) {
    case INV_F:
    case INV_W3W3:
    case INV_W3W7:
      return true;
  }
  return false;
}

/** PARI `modinv_height_factor` (`polmodular.c:114-149`). */
export function modinv_height_factor(inv: number): number {
  switch (inv) {
    case INV_J:
      return 1;
    case INV_G2:
      return 3;
    case INV_F:
      return 72;
    case INV_F2:
      return 36;
    case INV_F3:
      return 24;
    case INV_F4:
      return 18;
    case INV_F8:
      return 9;
    case INV_W2W3:
      return 72;
    case INV_W3W3:
      return 36;
    case INV_W2W5:
      return 54;
    case INV_W2W7:
      return 48;
    case INV_W3W5:
      return 36;
    case INV_W2W13:
      return 42;
    case INV_W3W7:
      return 32;
    case INV_W2W3E2:
      return 36;
    case INV_W2W5E2:
      return 27;
    case INV_W2W7E2:
      return 24;
    case INV_W3W3E2:
      return 18;
    case INV_W5W7:
      return 24;
    case INV_W3W13:
      return 28;
    case INV_ATKIN3:
    case INV_ATKIN5:
    case INV_ATKIN7:
    case INV_ATKIN11:
    case INV_ATKIN13:
    case INV_ATKIN17:
    case INV_ATKIN19:
    case INV_ATKIN23:
    case INV_ATKIN29:
    case INV_ATKIN31:
      return (inv - 99) / 2;
    default:
      throw new PariBugError('modinv_height_factor');
  }
}

/** PARI `modinv_sparse_factor` (`polmodular.c:185-221`). */
export function modinv_sparse_factor(inv: number): number {
  switch (inv) {
    case INV_G2:
    case INV_F8:
    case INV_W3W5:
    case INV_W2W5E2:
    case INV_W3W3E2:
      return 3;
    case INV_F:
      return 24;
    case INV_F2:
    case INV_W2W3:
      return 12;
    case INV_F3:
      return 8;
    case INV_F4:
    case INV_W2W3E2:
    case INV_W2W5:
    case INV_W3W3:
      return 6;
    case INV_W2W7:
      return 4;
    case INV_W2W7E2:
    case INV_W2W13:
    case INV_W3W7:
      return 2;
  }
  return 1;
}

/* `polmodular.c:223-226` */
const IQ_FILTER_1MOD3 = 1;
const IQ_FILTER_1MOD4 = 4;

/** PARI `modinv_pfilter` (`polmodular.c:223-248`). */
export function modinv_pfilter(inv: number): number {
  switch (inv) {
    case INV_G2:
    case INV_W3W3:
    case INV_W3W3E2:
    case INV_W3W5:
    case INV_W2W5:
    case INV_W2W3E2:
    case INV_W2W5E2:
    case INV_W3W13:
      return IQ_FILTER_1MOD3; /* ensure unique cube roots */
    case INV_W2W7:
    case INV_F3:
      return IQ_FILTER_1MOD4; /* ensure at most two 4th/8th roots */
    case INV_F:
    case INV_F2:
    case INV_F4:
    case INV_F8:
    case INV_W2W3:
      return IQ_FILTER_1MOD3 | IQ_FILTER_1MOD4;
  }
  return 0;
}

/** PARI `modinv_good_prime` (`polmodular.c:249-272`). */
export function modinv_good_prime(inv: number, p: number): boolean {
  switch (inv) {
    case INV_G2:
    case INV_W2W3E2:
    case INV_W3W3:
    case INV_W3W3E2:
    case INV_W3W5:
    case INV_W2W5E2:
    case INV_W2W5:
      return p % 3 === 2;
    case INV_W2W7:
    case INV_F3:
      return (p & 3) !== 1;
    case INV_F2:
    case INV_F4:
    case INV_F8:
    case INV_F:
    case INV_W2W3:
      return p % 3 === 2 && (p & 3) !== 1;
  }
  return true;
}

/** PARI `prime_to_conductor` (`polmodular.c:274-283`). */
function prime_to_conductor(D: number, p: number): boolean {
  if (p > 2) return D % (p * p) !== 0;
  const b = D & 0xf;
  return b !== 0 && b !== 4; /* 2 | cond(D) <=> D = 0,4 mod 16 */
}

/** PARI `red_primeform` (`polmodular.c:284-294`); `null` if `p | cond(D)`. */
function red_primeform(D: number, p: number): Qfb | null {
  if (!prime_to_conductor(D, p)) return null;
  return qfbred(primeform(BigInt(D), BigInt(p)));
}

/**
 * PARI `qfb_nform` (`polmodular.c:295-313`): product of prime forms over the
 * primes appearing in the factorisation of `n` (with multiplicity).
 */
export function qfb_nform(D: number, n: number): Qfb | null {
  const fa = factoru(BigInt(n));
  let N: Qfb | null = null;
  for (let i = 0; i < fa.length; i++) {
    const Q = red_primeform(D, Number(fa[i]![0]));
    if (!Q) return null;
    const e = Number(fa[i]![1]);
    let j: number;
    if (i === 0) {
      N = Q;
      j = 1;
    } else {
      j = 0;
    }
    for (; j < e; ++j) N = qfbcomp(Q, N!);
  }
  return N;
}

/** PARI `qfb_is_two_torsion` (`polmodular.c:315-321`). */
function qfb_is_two_torsion(x: Qfb): boolean {
  return x.a === 1n || x.b === 0n || x.a === x.b || x.a === x.c;
}

/** PARI `qfb_distinct_prods` (`polmodular.c:323-338`). */
function qfb_distinct_prods(D: number, p1: number, p2: number): boolean {
  let P1 = red_primeform(D, p1);
  if (!P1) return false;
  P1 = qfbsqr(P1);
  let P2 = red_primeform(D, p2);
  if (!P2) return false;
  P2 = qfbsqr(P2);
  return !(P1.a === P2.a && (P1.b < 0n ? -P1.b : P1.b) === (P2.b < 0n ? -P2.b : P2.b));
}

/** PARI `modinv_double_eta_good_disc` (`polmodular.c:340-380`). */
function modinv_double_eta_good_disc(D: number, inv: number): boolean {
  const pp = { p1: 1, p2: 1 };
  const N = modinv_degree(pp, inv);
  if (!N) return false;
  const i1 = kross(D, pp.p1);
  if (i1 < 0) return false;
  /* Exclude ramified case for w_{p,p} */
  if (pp.p1 === pp.p2 && i1 === 0) return false;
  const i2 = kross(D, pp.p2);
  if (i2 < 0) return false;
  let P = red_primeform(D, pp.p1);
  if (!P || P.a === 1n || (i1 !== 0 && qfb_is_two_torsion(P))) return false;
  if (pp.p1 === pp.p2) return !qfb_is_two_torsion(qfbsqr(P));

  P = red_primeform(D, pp.p2);
  if (!P || P.a === 1n || (i2 !== 0 && qfb_is_two_torsion(P))) return false;

  if (i1 > 0 && i2 > 0 && !qfb_distinct_prods(D, pp.p1, pp.p2)) return false;
  if (i1 === 0 && i2 === 0) {
    const Q = qfb_nform(D, N);
    if (Q && Q.a === 1n) return false;
  }
  return true;
}

/** PARI `modinv_ramified` (`polmodular.c:388-396`). */
export function modinv_ramified(D: number, inv: number, out: { N: number }): boolean {
  const pp = { p1: 1, p2: 1 };
  out.N = modinv_degree(pp, inv);
  if (out.N <= 1) return false;
  return D % pp.p1 === 0 && D % pp.p2 === 0;
}

/** PARI `modinv_good_atkin` (`polmodular.c:398-412`). */
function modinv_good_atkin(L: number, D: number): boolean {
  const L2 = L * L;
  if (kross(D, L) < 0 || -D % L2 === 0) return false;
  if (-D > 4 * L2) return true;
  let q = red_primeform(D, L);
  if (!q || q.a === 1n) return false;
  if (D % L === 0) return true;
  q = qfbsqr(q);
  return q.a !== 1n;
}

/** PARI `modinv_good_disc` (`polmodular.c:413-469`). */
export function modinv_good_disc(inv: number, D: number): boolean {
  switch (inv) {
    case INV_J:
      return true;
    case INV_G2:
      return D % 3 !== 0;
    case INV_F3:
      return (-D & 7) === 7;
    case INV_F:
    case INV_F2:
    case INV_F4:
    case INV_F8:
      return (-D & 7) === 7 && D % 3 !== 0;
    case INV_W3W5:
      return D % 3 !== 0 && modinv_double_eta_good_disc(D, inv);
    case INV_W3W3E2:
      return D % 3 !== 0 && modinv_double_eta_good_disc(D, inv);
    case INV_W3W3:
      return (D & 1) !== 0 && D % 3 !== 0 && modinv_double_eta_good_disc(D, inv);
    case INV_W2W3E2:
      return D % 3 !== 0 && modinv_double_eta_good_disc(D, inv);
    case INV_W2W3:
      return (-D & 7) === 7 && D % 3 !== 0 && modinv_double_eta_good_disc(D, inv);
    case INV_W2W5:
      return -D % 80 !== 20 && D % 3 !== 0 && modinv_double_eta_good_disc(D, inv);
    case INV_W2W5E2:
      return D % 3 !== 0 && modinv_double_eta_good_disc(D, inv);
    case INV_W2W7E2:
      return -D % 112 !== 84 && modinv_double_eta_good_disc(D, inv);
    case INV_W2W7:
      return (-D & 7) === 7 && modinv_double_eta_good_disc(D, inv);
    case INV_W2W13:
      return -D % 208 !== 52 && modinv_double_eta_good_disc(D, inv);
    case INV_W3W7:
      return (D & 1) !== 0 && -D % 21 !== 0 && modinv_double_eta_good_disc(D, inv);
    case INV_W5W7:
      return modinv_double_eta_good_disc(D, inv);
    case INV_W3W13:
      return (D & 1) !== 0 && D % 3 !== 0 && modinv_double_eta_good_disc(D, inv);
    case INV_ATKIN3:
    case INV_ATKIN5:
    case INV_ATKIN7:
    case INV_ATKIN11:
    case INV_ATKIN13:
    case INV_ATKIN17:
    case INV_ATKIN19:
    case INV_ATKIN23:
    case INV_ATKIN29:
    case INV_ATKIN31:
      return modinv_good_atkin(inv - 100, D);
  }
  throw new PariBugError('modinv_good_disc');
}

/** PARI `disc_best_modinv` (`polmodular.c:152-184`). */
export function disc_best_modinv(D: number): number {
  const order = [
    INV_F,
    INV_W2W3,
    INV_W2W5,
    INV_W2W7,
    INV_W2W13,
    INV_W3W3,
    INV_W2W3E2,
    INV_W3W5,
    INV_W3W7,
    INV_W3W13,
    INV_W2W5E2,
    INV_F3,
    INV_W2W7E2,
    INV_W5W7,
    INV_W3W3E2,
    INV_ATKIN31,
    INV_ATKIN29,
    INV_ATKIN23,
    INV_ATKIN19,
    INV_ATKIN17,
    INV_ATKIN13,
    INV_ATKIN11,
    INV_ATKIN7,
    INV_ATKIN5,
    INV_G2,
    INV_ATKIN3,
  ];
  for (const inv of order) if (modinv_good_disc(inv, D)) return inv;
  return INV_J;
}

/** PARI `modinv_is_Weber` (`polmodular.c:470-475`). */
export function modinv_is_Weber(inv: number): boolean {
  return inv === INV_F || inv === INV_F2 || inv === INV_F3 || inv === INV_F4 || inv === INV_F8;
}

/** PARI `modinv_is_double_eta` (`polmodular.c:477-503`). */
export function modinv_is_double_eta(inv: number): boolean {
  switch (inv) {
    case INV_W2W3:
    case INV_W2W3E2:
    case INV_W2W5:
    case INV_W2W5E2:
    case INV_W2W7:
    case INV_W2W7E2:
    case INV_W2W13:
    case INV_W3W3:
    case INV_W3W3E2:
    case INV_W3W5:
    case INV_W3W7:
    case INV_W5W7:
    case INV_W3W13:
    case INV_ATKIN3:
    case INV_ATKIN5:
    case INV_ATKIN7:
    case INV_ATKIN11:
    case INV_ATKIN13:
    case INV_ATKIN17:
    case INV_ATKIN19:
    case INV_ATKIN23:
    case INV_ATKIN29:
    case INV_ATKIN31:
      return true;
  }
  return false;
}

/** PARI `modinv_parent` (`polmodular.c:2296-2310`). */
function modinv_parent(inv: number): number {
  switch (inv) {
    case INV_F2:
    case INV_F4:
    case INV_F8:
      return INV_F;
    case INV_W2W3E2:
      return INV_W2W3;
    case INV_W2W5E2:
      return INV_W2W5;
    case INV_W2W7E2:
      return INV_W2W7;
    case INV_W3W3E2:
      return INV_W3W3;
    default:
      throw new PariBugError('modinv_parent');
  }
}

/** PARI `modinv_parent_power` (`polmodular.c:2312-2326`). */
function modinv_parent_power(inv: number): number {
  switch (inv) {
    case INV_F4:
      return 4;
    case INV_F8:
      return 8;
    case INV_F2:
    case INV_W2W3E2:
    case INV_W2W5E2:
    case INV_W2W7E2:
    case INV_W3W3E2:
      return 2;
    default:
      throw new PariBugError('modinv_parent_power');
  }
}

/** PARI `modinv_max_internal_level` (`polmodular.c:1824-1860`). */
export function modinv_max_internal_level(inv: number): number {
  switch (inv) {
    case INV_J:
      return 5;
    case INV_G2:
      return 2;
    case INV_F:
    case INV_F2:
    case INV_F4:
    case INV_F8:
      return 5;
    case INV_W2W5:
    case INV_W2W5E2:
      return 7;
    case INV_W2W3:
    case INV_W2W3E2:
    case INV_W3W3:
    case INV_W3W7:
      return 5;
    case INV_W3W3E2:
      return 2;
    case INV_F3:
    case INV_W2W7:
    case INV_W2W7E2:
    case INV_W2W13:
      return 3;
    case INV_W3W5:
    case INV_W5W7:
    case INV_W3W13:
    case INV_ATKIN3:
    case INV_ATKIN5:
    case INV_ATKIN7:
    case INV_ATKIN11:
    case INV_ATKIN13:
    case INV_ATKIN17:
    case INV_ATKIN19:
    case INV_ATKIN23:
    case INV_ATKIN29:
    case INV_ATKIN31:
      return 2;
  }
  throw new PariBugError('modinv_max_internal_level');
}

/* ================================================================== */
/* Small helpers                                                      */
/* ================================================================== */

/** PARI `kross(long, long)` — Kronecker symbol on machine integers. */
function kross(x: number, y: number): number {
  return kronecker(BigInt(x), BigInt(y));
}

/** PARI `uu32toi(a, b)` = `a * 2^32 + b` (`paripriv`/`gen`). */
function uu32toi(a: number, b: number): bigint {
  return (BigInt(a) << 32n) + BigInt(b);
}

/** PARI `mkintn(n, w_{n-1}, ..., w_0)`: base-2^32 limbs, most significant first. */
function mkintn(...words: number[]): bigint {
  let r = 0n;
  for (const w of words) r = (r << 32n) + BigInt(w >>> 0);
  return r;
}

/** `ceildivuu(a, b)` */
function ceildivuu(a: number, b: number): number {
  return Math.floor((a + b - 1) / b);
}

/** PARI `u_lvalrem(x, L, &n)`: returns v_L(x), sets n = x / L^v. */
function u_lvalrem(x: bigint, L: bigint): { v: number; n: bigint } {
  let v = 0;
  let n = x;
  while (n % L === 0n) {
    n /= L;
    v++;
  }
  return { v, n };
}

/** PARI `z_lval(x, L)` */
function z_lval(x: number, L: number): number {
  let v = 0;
  let n = Math.abs(x);
  if (n === 0) return 0;
  while (n % L === 0) {
    n /= L;
    v++;
  }
  return v;
}

function ugcd(a: number, b: number): number {
  a = Math.abs(a);
  b = Math.abs(b);
  while (b) {
    const t = a % b;
    a = b;
    b = t;
  }
  return a;
}

/** PARI `pari_PRIMES[i]`, 1-indexed: 2, 3, 5, ... */
const PARI_PRIMES = [
  0, 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
  101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193,
  197, 199, 211, 223, 227, 229, 233, 239, 241, 251,
];

/** PARI `unextprime(n)`: smallest prime >= n (n >= 1). */
function unextprime(n: number): number {
  if (n <= 2) return 2;
  let m = n | 1;
  if (m < n) m += 2;
  for (;;) {
    if (isPrime(BigInt(m))) return m;
    m += 2;
  }
}

/* ================================================================== */
/* Modular polynomials of small level (`polmodular.c:2036-2455`)      */
/* ================================================================== */

/*  Phi2, the modular polynomial of level 2 (`polmodular.c:2039-2079`):
 *
 *  X^3 + X^2 * (-Y^2 + 1488*Y - 162000)
 *      + X * (1488*Y^2 + 40773375*Y + 8748000000)
 *      + Y^3 - 162000*Y^2 + 8748000000*Y - 157464000000000  */
function phi2_ZV(): bigint[] {
  return [-uu32toi(36662, 1908994048), uu32toi(2, 158065408), 40773375n, -162000n, 1488n, -1n];
}

/* `polmodular.c:2081-2113` */
function phi3_ZV(): bigint[] {
  return [
    0n,
    uu32toi(100, 2503270400) << 32n,
    -uu32toi(179476562, 2147483648),
    uu32toi(105468, 3221225472),
    uu32toi(2072, 1050738688),
    2587918086n,
    36864000n,
    -1069956n,
    2232n,
    -1n,
  ];
}

/* `polmodular.c:2115-2146` */
function phi5_ZV(): bigint[] {
  return [
    mkintn(0x18c2cc9c, 0x484382b2, 0xdc000000, 0x0, 0x0),
    mkintn(0x2638f, 0x2ff02690, 0x68026000, 0x0, 0x0),
    -mkintn(0x308, 0xac9d9a4, 0xe0fdab12, 0xc0000000, 0x0),
    mkintn(0x13, 0xaae09f9d, 0x1b5ef872, 0x30000000, 0x0),
    mkintn(0x1b802fa9, 0x77ba0653, 0xd2f78000, 0x0),
    mkintn(0xfbfd, 0x278e4756, 0xdf08a7c4, 0x40000000),
    mkintn(0x35f922, 0x62ccea6f, 0x153d0000, 0x0),
    -mkintn(0x97d, 0x29203faf, 0xc3036909, 0x80000000),
    mkintn(0x56e9e892, 0xd7781867, 0xf2ea0000),
    -mkintn(0x5d6d, 0xe0a58f4e, 0x9ee68c14),
    mkintn(0x1100d, 0x85cea769, 0x40000000),
    mkintn(0x1b38, 0x43cf461f, 0x3a900000),
    mkintn(0x14, 0xc45a616e, 0x4801680f),
    uu32toi(0x17f4350, 0x493ca3e0),
    uu32toi(0x183, 0xe54ce1f8),
    uu32toi(0x1c9, 0x18860000),
    -uu32toi(0x39, 0x6f7a2206),
    2028551200n,
    -4550940n,
    3720n,
    -1n,
  ];
}

/** A sparse vector of length `n` with the given `[index (1-based), value]` pairs. */
function sparseZV(n: number, entries: Array<[number, bigint]>): bigint[] {
  const v = new Array<bigint>(n).fill(0n);
  for (const [i, x] of entries) v[i - 1] = x;
  return v;
}

/* `polmodular.c:2148-2295` */
const phi5_f_ZV = (): bigint[] =>
  sparseZV(21, [
    [3, 4n],
    [21, -1n],
  ]);
const phi3_f3_ZV = (): bigint[] =>
  sparseZV(10, [
    [3, 8n],
    [10, -1n],
  ]);
const phi2_g2_ZV = (): bigint[] => [-54000n, 0n, 495n, 0n, 0n, -1n];
const phi5_w2w3_ZV = (): bigint[] =>
  sparseZV(21, [
    [3, -1n],
    [10, 5n],
    [21, -1n],
  ]);
const phi7_w2w5_ZV = (): bigint[] =>
  sparseZV(36, [
    [3, -1n],
    [15, 56n],
    [19, 42n],
    [24, 21n],
    [30, 7n],
    [36, -1n],
  ]);
const phi5_w3w3_ZV = (): bigint[] =>
  sparseZV(21, [
    [3, 9n],
    [6, -15n],
    [15, 5n],
    [21, -1n],
  ]);
const phi3_w2w7_ZV = (): bigint[] =>
  sparseZV(10, [
    [3, -1n],
    [6, 3n],
    [10, -1n],
  ]);
const phi2_w3w5_ZV = (): bigint[] => [0n, 0n, 1n, 0n, 0n, -1n];
const phi5_w3w7_ZV = (): bigint[] =>
  sparseZV(21, [
    [3, -1n],
    [6, 10n],
    [8, 5n],
    [10, 35n],
    [13, 20n],
    [15, 10n],
    [17, 5n],
    [19, 5n],
    [21, -1n],
  ]);
const phi3_w2w13_ZV = (): bigint[] =>
  sparseZV(10, [
    [3, -1n],
    [6, 3n],
    [8, 3n],
    [10, -1n],
  ]);
const phi2_w3w3e2_ZV = (): bigint[] => [0n, 0n, 3n, 0n, 0n, -1n];
const phi2_w5w7_ZV = (): bigint[] => [0n, 0n, 1n, 0n, 2n, -1n];
const phi2_w3w13_ZV = (): bigint[] => [0n, 0n, -1n, 0n, 2n, -1n];
const phi2_atkin3_ZV = (): bigint[] => [28166076n, 741474n, 17343n, 1566n, 0n, -1n];
const phi2_atkin5_ZV = (): bigint[] => [323456n, 24244n, 1519n, 268n, 0n, -1n];
const phi2_atkin7_ZV = (): bigint[] => [27100n, 3810n, 407n, 102n, 0n, -1n];
const phi2_atkin11_ZV = (): bigint[] => [1600n, 470n, 91n, 34n, 0n, -1n];
const phi2_atkin13_ZV = (): bigint[] => [656n, 240n, 55n, 24n, 0n, -1n];
const phi2_atkin17_ZV = (): bigint[] => [156n, 86n, 27n, 14n, 0n, -1n];
const phi2_atkin19_ZV = (): bigint[] => [100n, 60n, 19n, 12n, 0n, -1n];
const phi2_atkin23_ZV = (): bigint[] => [2n, 6n, 9n, 4n, 2n, -1n];
const phi2_atkin29_ZV = (): bigint[] => [0n, 0n, 3n, 2n, 2n, -1n];
const phi2_atkin31_ZV = (): bigint[] => [-2n, 0n, 1n, 2n, 2n, -1n];

/**
 * PARI `internal_db` (`polmodular.c:2397-2441`): the coefficient vector of the
 * modular polynomial of smallest level for `inv`, or `null` when the answer has
 * to be obtained through `polmodular0_powerup_ZM`.
 *
 * NOTE: upstream's `case INV_J:` falls through into `case INV_F:` when
 * `L` is not 2, 3 or 5. That path is unreachable because `polmodular_small_ZM`
 * is only ever called with `L <= modinv_max_internal_level(inv) = 5`; we mirror
 * the fall-through anyway so that behaviour is bit-identical.
 */
function internal_db(L: number, inv: number): bigint[] | null {
  if (inv === INV_J) {
    switch (L) {
      case 2:
        return phi2_ZV();
      case 3:
        return phi3_ZV();
      case 5:
        return phi5_ZV();
      default:
        return phi5_f_ZV(); /* upstream fall-through, unreachable */
    }
  }
  switch (inv) {
    case INV_F:
      return phi5_f_ZV();
    case INV_F2:
      return null;
    case INV_F3:
      return phi3_f3_ZV();
    case INV_F4:
      return null;
    case INV_G2:
      return phi2_g2_ZV();
    case INV_W2W3:
      return phi5_w2w3_ZV();
    case INV_F8:
      return null;
    case INV_W3W3:
      return phi5_w3w3_ZV();
    case INV_W2W5:
      return phi7_w2w5_ZV();
    case INV_W2W7:
      return phi3_w2w7_ZV();
    case INV_W3W5:
      return phi2_w3w5_ZV();
    case INV_W3W7:
      return phi5_w3w7_ZV();
    case INV_W2W3E2:
      return null;
    case INV_W2W5E2:
      return null;
    case INV_W2W13:
      return phi3_w2w13_ZV();
    case INV_W2W7E2:
      return null;
    case INV_W3W3E2:
      return phi2_w3w3e2_ZV();
    case INV_W5W7:
      return phi2_w5w7_ZV();
    case INV_W3W13:
      return phi2_w3w13_ZV();
    case INV_ATKIN3:
      return phi2_atkin3_ZV();
    case INV_ATKIN5:
      return phi2_atkin5_ZV();
    case INV_ATKIN7:
      return phi2_atkin7_ZV();
    case INV_ATKIN11:
      return phi2_atkin11_ZV();
    case INV_ATKIN13:
      return phi2_atkin13_ZV();
    case INV_ATKIN17:
      return phi2_atkin17_ZV();
    case INV_ATKIN19:
      return phi2_atkin19_ZV();
    case INV_ATKIN23:
      return phi2_atkin23_ZV();
    case INV_ATKIN29:
      return phi2_atkin29_ZV();
    case INV_ATKIN31:
      return phi2_atkin31_ZV();
  }
  throw new PariBugError('internal_db');
}

/** Allocate a zero `n x m` matrix in PARI column-major layout. */
function zero_ZM(rows: number, cols: number): ZM {
  const M: ZM = new Array(cols);
  for (let j = 0; j < cols; j++) M[j] = new Array<bigint>(rows).fill(0n);
  return M;
}

/**
 * PARI `sympol_to_ZM` (`polmodular.c:1810-1821`): desymmetrise the
 * `(L+1)(L+2)/2` coefficients of a symmetric `Phi_L` into an
 * `(L+2) x (L+2)` matrix (the coefficient of `X^{L+1}`/`Y^{L+1}`, always 1,
 * is not part of the input).
 */
export function sympol_to_ZM(phi: readonly bigint[], L: number): ZM {
  const res = zero_ZM(L + 2, L + 2);
  let c = 0;
  for (let i = 1; i <= L + 1; ++i)
    for (let j = 1; j <= i; ++j, ++c) {
      const v = phi[c]!;
      res[j - 1]![i - 1] = v; /* gcoeff(res, i, j) */
      res[i - 1]![j - 1] = v; /* gcoeff(res, j, i) */
    }
  res[0]![L + 1] = 1n;
  res[L + 1]![0] = 1n;
  return res;
}

/* ================================================================== */
/* The polmodular database (`polmodular.c:941-1029`)                  */
/* ================================================================== */

/**
 * PARI's polmodular database: a pair of tables indexed by level `L`, the first
 * for `INV_J`, the second for the requested invariant (`gen_0` when the
 * requested invariant *is* `INV_J`).
 */
export interface PolmodularDB {
  j: Array<ZM | null>;
  inv: Array<ZM | null> | null;
}

/** PARI `polmodular_db_init` (`polmodular.c:945-953`). */
export function polmodular_db_init(inv: number): PolmodularDB {
  const LEN = 32;
  return {
    j: new Array<ZM | null>(LEN + 1).fill(null),
    inv: inv === INV_J ? null : new Array<ZM | null>(LEN + 1).fill(null),
  };
}

/** PARI `polmodular_db_for_inv` (`polmodular.c:987-988`). */
export function polmodular_db_for_inv(db: PolmodularDB, inv: number): Array<ZM | null> {
  return inv === INV_J ? db.j : db.inv!;
}

/** PARI `polmodular_db_add_level` (`polmodular.c:955-979`). */
export function polmodular_db_add_level(db: PolmodularDB, L: number, inv: number): void {
  const key = inv === INV_J ? 'j' : 'inv';
  let tab = key === 'j' ? db.j : db.inv!;
  const max_L = tab.length - 1;
  if (L > max_L) {
    const newlen = 2 * L;
    const newtab = new Array<ZM | null>(newlen + 1).fill(null);
    for (let i = 1; i <= max_L; ++i) newtab[i] = tab[i]!;
    if (key === 'j') db.j = newtab;
    else db.inv = newtab;
    tab = newtab;
  }
  if (tab[L] == null) {
    /* may set tab[L] recursively */
    const x = polmodular0_ZM(L, inv, null, null, false, db);
    (key === 'j' ? db.j : db.inv!)[L] = x;
  }
}

/** PARI `polmodular_db_add_levels` (`polmodular.c:981-985`). */
export function polmodular_db_add_levels(
  db: PolmodularDB,
  levels: readonly number[],
  inv: number
): void {
  for (const L of levels) polmodular_db_add_level(db, L, inv);
}

/** PARI `ZM_to_Flm`. */
function ZM_to_Flm(M: ZM, p: bigint): bigint[][] {
  return M.map((col) => col.map((x) => Fp_red(x, p)));
}

/** PARI `polmodular_db_getp` (`polmodular.c:991-998`): `Phi_L mod p`. */
export function polmodular_db_getp(db: Array<ZM | null>, L: number, p: bigint): bigint[][] {
  const f = db[L];
  if (f == null) throw new PariBugError('polmodular_db_getp');
  return ZM_to_Flm(f, p);
}

/* ================================================================== */
/* Evaluation of small-level modular polynomials                      */
/* (`polmodular.c:1035-1170`)                                         */
/* ================================================================== */

/* `Flm` here is the same column-major convention: `M[j][i]` = row i, col j. */
type Flm = bigint[][];

function mcoeff(M: Flm, i: number, j: number): bigint {
  return M[j - 1]![i - 1]!;
}

/**
 * PARI `Fl_addmul2/3/4/5` (`polmodular.c:874-940`): `sum_i x_i y_{n-1-i} mod p`
 * (note the reversal of the second argument list).
 */
function Fl_addmulrev(xs: readonly bigint[], ys: readonly bigint[], p: bigint): bigint {
  const n = xs.length;
  let s = 0n;
  for (let i = 0; i < n; i++) s += Fp_mul(xs[i]!, ys[n - 1 - i]!, p);
  return s % p;
}

/** PARI `Flm_Fl_phi2_evalx` (`polmodular.c:1039-1064`), as a coefficient list. */
function Flm_Fl_phi2_evalx(phi2: Flm, j: bigint, p: bigint): bigint[] {
  const j2 = Fp_mul(j, j, p);
  const J = [j, j2];
  const res = new Array<bigint>(4);

  let t1 = (j + mcoeff(phi2, 3, 1)) % p;
  t1 = Fl_addmulrev(J, [t1, mcoeff(phi2, 2, 1)], p);
  res[0] = (t1 + mcoeff(phi2, 1, 1)) % p;

  t1 = Fl_addmulrev(J, [mcoeff(phi2, 3, 2), mcoeff(phi2, 2, 2)], p);
  res[1] = (t1 + mcoeff(phi2, 2, 1)) % p;

  t1 = Fp_mul(j, mcoeff(phi2, 3, 2), p);
  t1 = (t1 + mcoeff(phi2, 3, 1)) % p;
  res[2] = (t1 - j2 + p) % p;

  res[3] = 1n;
  return res;
}

/** PARI `Flm_Fl_phi3_evalx` (`polmodular.c:1066-1097`). */
function Flm_Fl_phi3_evalx(phi3: Flm, j: bigint, p: bigint): bigint[] {
  const j2 = Fp_mul(j, j, p);
  const j3 = Fp_mul(j, j2, p);
  const J = [j, j2, j3];
  const res = new Array<bigint>(5);

  let t1 = (j + mcoeff(phi3, 4, 1)) % p;
  t1 = Fl_addmulrev(J, [t1, mcoeff(phi3, 3, 1), mcoeff(phi3, 2, 1)], p);
  res[0] = (t1 + mcoeff(phi3, 1, 1)) % p;

  t1 = Fl_addmulrev(J, [mcoeff(phi3, 4, 2), mcoeff(phi3, 3, 2), mcoeff(phi3, 2, 2)], p);
  res[1] = (t1 + mcoeff(phi3, 2, 1)) % p;

  t1 = Fl_addmulrev(J, [mcoeff(phi3, 4, 3), mcoeff(phi3, 3, 3), mcoeff(phi3, 3, 2)], p);
  res[2] = (t1 + mcoeff(phi3, 3, 1)) % p;

  t1 = Fl_addmulrev([j, j2], [mcoeff(phi3, 4, 3), mcoeff(phi3, 4, 2)], p);
  t1 = (t1 + mcoeff(phi3, 4, 1)) % p;
  res[3] = (t1 - j3 + p) % p;

  res[4] = 1n;
  return res;
}

/** PARI `Flm_Fl_phi5_evalx` (`polmodular.c:1099-1154`). */
function Flm_Fl_phi5_evalx(phi5: Flm, j: bigint, p: bigint): bigint[] {
  const j2 = Fp_mul(j, j, p);
  const j3 = Fp_mul(j, j2, p);
  const j4 = Fp_mul(j2, j2, p);
  const j5 = Fp_mul(j, j4, p);
  const J = [j, j2, j3, j4, j5];
  const res = new Array<bigint>(7);

  let t1 = (j + mcoeff(phi5, 6, 1)) % p;
  t1 = Fl_addmulrev(
    J,
    [t1, mcoeff(phi5, 5, 1), mcoeff(phi5, 4, 1), mcoeff(phi5, 3, 1), mcoeff(phi5, 2, 1)],
    p
  );
  res[0] = (t1 + mcoeff(phi5, 1, 1)) % p;

  t1 = Fl_addmulrev(
    J,
    [
      mcoeff(phi5, 6, 2),
      mcoeff(phi5, 5, 2),
      mcoeff(phi5, 4, 2),
      mcoeff(phi5, 3, 2),
      mcoeff(phi5, 2, 2),
    ],
    p
  );
  res[1] = (t1 + mcoeff(phi5, 2, 1)) % p;

  t1 = Fl_addmulrev(
    J,
    [
      mcoeff(phi5, 6, 3),
      mcoeff(phi5, 5, 3),
      mcoeff(phi5, 4, 3),
      mcoeff(phi5, 3, 3),
      mcoeff(phi5, 3, 2),
    ],
    p
  );
  res[2] = (t1 + mcoeff(phi5, 3, 1)) % p;

  t1 = Fl_addmulrev(
    J,
    [
      mcoeff(phi5, 6, 4),
      mcoeff(phi5, 5, 4),
      mcoeff(phi5, 4, 4),
      mcoeff(phi5, 4, 3),
      mcoeff(phi5, 4, 2),
    ],
    p
  );
  res[3] = (t1 + mcoeff(phi5, 4, 1)) % p;

  t1 = Fl_addmulrev(
    J,
    [
      mcoeff(phi5, 6, 5),
      mcoeff(phi5, 5, 5),
      mcoeff(phi5, 5, 4),
      mcoeff(phi5, 5, 3),
      mcoeff(phi5, 5, 2),
    ],
    p
  );
  res[4] = (t1 + mcoeff(phi5, 5, 1)) % p;

  t1 = Fl_addmulrev(
    [j, j2, j3, j4],
    [mcoeff(phi5, 6, 5), mcoeff(phi5, 6, 4), mcoeff(phi5, 6, 3), mcoeff(phi5, 6, 2)],
    p
  );
  t1 = (t1 + mcoeff(phi5, 6, 1)) % p;
  res[5] = (t1 - j5 + p) % p;

  res[6] = 1n;
  return res;
}

/**
 * PARI `Flm_Fl_polmodular_evalx` (`polmodular.c:1155-1169`): `Phi_L(X, j) mod p`
 * as a coefficient list in `X` (index `i` = coefficient of `X^i`).
 */
export function Flm_Fl_polmodular_evalx(phi: Flm, L: number, j: bigint, p: bigint): bigint[] {
  switch (L) {
    case 2:
      return Flm_Fl_phi2_evalx(phi, j, p);
    case 3:
      return Flm_Fl_phi3_evalx(phi, j, p);
    case 5:
      return Flm_Fl_phi5_evalx(phi, j, p);
    default: {
      const jp = Fl_powers(j, L + 1, p);
      return Flm_Flc_mul(phi, jp, p);
    }
  }
}

/** PARI `Fl_powers_pre(x, n, p, pi)`: `[1, x, ..., x^n]`. */
function Fl_powers(x: bigint, n: number, p: bigint): bigint[] {
  const v = new Array<bigint>(n + 1);
  v[0] = 1n % p;
  for (let i = 1; i <= n; i++) v[i] = Fp_mul(v[i - 1]!, x, p);
  return v;
}

/** PARI `Flm_Flc_mul`: `M * c` (M column-major). */
function Flm_Flc_mul(M: Flm, c: readonly bigint[], p: bigint): bigint[] {
  const rows = M[0]!.length;
  const out = new Array<bigint>(rows).fill(0n);
  for (let jj = 0; jj < M.length && jj < c.length; jj++) {
    const cj = c[jj]!;
    if (cj === 0n) continue;
    const col = M[jj]!;
    for (let i = 0; i < rows; i++) out[i] = (out[i]! + Fp_mul(col[i]!, cj, p)) % p;
  }
  return out;
}

/* ================================================================== */
/* Top-level entry points (`polmodular.c:1864-2034`)                  */
/* ================================================================== */

/**
 * PARI `polmodular_small_ZM` (`polmodular.c:2442-2448`). Only valid when
 * `L <= modinv_max_internal_level(inv)`.
 */
function polmodular_small_ZM(L: number, inv: number, db: PolmodularDB): ZM {
  const f = internal_db(L, inv);
  if (!f) return polmodular0_powerup_ZM(L, inv, db);
  return sympol_to_ZM(f, L);
}

/**
 * PARI `polmodular0_ZM` (`polmodular.c:1865-1935`).
 *
 * When `J`/`Q` are given the result is the vector of `Phi_L(X, J) mod Q`
 * (and its first two derivatives if `compute_derivs`), packed as the columns
 * of the returned matrix.
 */
export function polmodular0_ZM(
  L: number,
  inv: number,
  J: bigint | null,
  Q: bigint | null,
  compute_derivs: boolean,
  db: PolmodularDB
): ZM {
  const lvl = modinv_level(inv);
  if (ugcd(L, lvl) !== 1)
    throw new PariDomainError('polmodular0_ZM', 'invariant', 'incompatible with', String(L));

  if (L <= modinv_max_internal_level(inv)) return polmodular_small_ZM(L, inv, db);

  return polmodular_CM_ZM(L, inv, J, Q, compute_derivs, db);
}

/** PARI `polmodular_ZM` (`polmodular.c:1942-1957`). */
export function polmodular_ZM(L: number, inv: number = INV_J): ZM {
  if (L < 2) throw new PariDomainError('polmodular_ZM', 'L', '<', '2');
  /* TODO (upstream too): handle nonprime L. */
  if (!isPrime(BigInt(L))) throw new PariImplError('composite level');
  const db = polmodular_db_init(inv);
  return polmodular0_ZM(L, inv, null, null, false, db);
}

/**
 * A bivariate polynomial over Z: `f[j][i]` is the coefficient of `X^i Y^j`
 * (this is exactly the `ZM` returned by {@link polmodular_ZM}).
 */
export type ZXX = ZM;

/**
 * PARI `polmodular_ZXX` (`polmodular.c:1959-1970`).
 *
 * PARI's `vx`/`vy` are variable *numbers*; the only thing they influence is the
 * priority check and which variable ends up outermost. We keep the check (so
 * the error behaviour matches) and always return the matrix with `X` indexing
 * rows and `Y` indexing columns, as `RgM_to_RgXX` does.
 */
export function polmodular_ZXX(L: number, inv: number = INV_J, vx = 0, vy = 1): ZXX {
  if (vx < 0) vx = 0;
  if (vy < 0) vy = 1;
  if (vx >= vy) throw new PariPriorityError('polmodular_ZXX', `x${vx}`, '<=', vy);
  return polmodular_ZM(L, inv);
}

/** PARI `FpV_deriv` (`polmodular.c:1972-1980`). */
function FpV_deriv(v: readonly bigint[], deg: number, P: bigint): bigint[] {
  const ln = v.length;
  const dv = new Array<bigint>(ln).fill(0n);
  let d = deg;
  for (let i = ln - 1; i > 0; i--, d--) dv[i] = Fp_mul(v[i - 1]!, BigInt(d) % P, P);
  dv[0] = 0n;
  return dv;
}

/**
 * PARI `Fp_polmodular_evalx` (`polmodular.c:1981-2010`): `Phi_L(X, J) mod P`,
 * as a coefficient list in `X`. With `compute_derivs`, returns
 * `[Phi, dPhi/dY, d^2Phi/dY^2]` evaluated at `Y = J`.
 */
export function Fp_polmodular_evalx(
  L: number,
  inv: number,
  J: bigint,
  P: bigint,
  compute_derivs = false
): bigint[] | [bigint[], bigint[], bigint[]] {
  if (L <= modinv_max_internal_level(inv)) {
    const phi = polmodular_ZM(L, inv).map((col) => col.map((x) => Fp_red(x, P)));
    let j_powers = Fp_powers(J, L + 1, P);
    const modpol = Flm_Flc_mul(phi, j_powers, P);
    if (compute_derivs) {
      j_powers = FpV_deriv(j_powers, L + 1, P);
      const d1 = Flm_Flc_mul(phi, j_powers, P);
      j_powers = FpV_deriv(j_powers, L + 1, P);
      const d2 = Flm_Flc_mul(phi, j_powers, P);
      return [modpol, d1, d2];
    }
    return modpol;
  }

  const db = polmodular_db_init(inv);
  const phi = polmodular0_ZM(L, inv, J, P, compute_derivs, db);
  return compute_derivs ? [phi[0]!, phi[1]!, phi[2]!] : phi[0]!;
}

/** PARI `Fp_powers(x, n, p)`. */
function Fp_powers(x: bigint, n: number, p: bigint): bigint[] {
  return Fl_powers(Fp_red(x, p), n, p);
}

/**
 * PARI `polmodular` (`polmodular.c:2011-2066`), restricted to the cases our
 * type system can express: `x = null` returns `Phi_L(X, Y)` over Z, otherwise
 * `x = {j, p}` is a `t_INTMOD` and we return `Phi_L(X, j) mod p`.
 */
export function polmodular(
  L: number,
  inv: number = INV_J,
  x: { j: bigint; p: bigint } | null = null,
  v = 1,
  compute_derivs = false
): ZXX | bigint[] | [bigint[], bigint[], bigint[]] {
  check_modinv(inv);
  if (!x) {
    if (compute_derivs) throw new Error('incorrect flag in polmodular');
    return polmodular_ZXX(L, inv, 0, v);
  }
  if (v < 0) v = 1;
  return Fp_polmodular_evalx(L, inv, x.j, x.p, compute_derivs);
}

/* ================================================================== */
/* The CM algorithm (`polmodular.c:1865-1935`, the `L > max` branch)   */
/* ================================================================== */

/* ------------------------------------------------------------------ */
/* Flx: polynomials over F_p (`Flx.c`, `FpX_factor.c`)                 */
/* `Flx` is `bigint[]` with `f[i]` the coefficient of `X^i`, trailing  */
/* zeros stripped; the zero polynomial is `[]`.                        */
/* ------------------------------------------------------------------ */

type Flx = bigint[];

function Flx_deg(f: Flx): number {
  return f.length - 1;
}

/** PARI `Flx_eval_pre`. */
function Flx_eval(f: Flx, x: bigint, p: bigint): bigint {
  let r = 0n;
  for (let i = f.length - 1; i >= 0; i--) r = (Fp_mul(r, x, p) + f[i]!) % p;
  return r;
}

/** PARI `Flx_deg1_root` (`Flx.c`): root of `a X + b`. */
function Flx_deg1_root(f: Flx, p: bigint): bigint {
  return Fp_mul(p - (f[0]! % p), Fp_inv(f[1]!, p), p) % p;
}

/** PARI `Flx_div_by_X_x(f, a, p, &rem)`: synthetic division by `X - a`. */
function Flx_div_by_X_x(f: Flx, a: bigint, p: bigint): { q: Flx; rem: bigint } {
  const n = f.length - 1;
  if (n < 0) return { q: [], rem: 0n };
  const q = new Array<bigint>(Math.max(n, 0));
  let carry = f[n]!;
  for (let i = n - 1; i >= 0; i--) {
    q[i] = carry;
    carry = (f[i]! + Fp_mul(carry, a, p)) % p;
  }
  return { q: FpX_renormalize(q), rem: carry };
}

/** PARI `Flx_remove_root` (`volcano.c:235-242`). */
function Flx_remove_root(f: Flx, a: bigint, p: bigint): Flx {
  const { q, rem } = Flx_div_by_X_x(f, a, p);
  if (rem !== 0n) throw new PariBugError('Flx_remove_root');
  return q;
}

/**
 * Distinct roots of `f` in `F_p` (PARI `Flx_roots_pre`, `FpX_factor.c`).
 * Cantor-Zassenhaus equal-degree splitting after `gcd(f, X^p - X)`.
 */
function Flx_roots(f: Flx, p: bigint): bigint[] {
  const g0 = FpX_normalize(FpX_renormalize(f.map((c) => Fp_red(c, p))), p);
  if (FpX_degree(g0) < 1) return [];
  /* g = gcd(f, X^p - X): the product of the distinct linear factors */
  const xp = FpXQ_pow([0n, 1n], p, g0, p);
  const g = FpX_normalize(FpX_gcd(g0, FpX_sub(xp, [0n, 1n], p), p), p);
  const out: bigint[] = [];
  Flx_split_linear(g, p, out);
  out.sort((a, b) => (a < b ? -1 : a > b ? 1 : 0));
  return out;
}

/** Equal-degree splitting of a squarefree product of linear factors. */
function Flx_split_linear(g: Flx, p: bigint, out: bigint[]): void {
  const d = FpX_degree(g);
  if (d < 1) return;
  if (d === 1) {
    out.push(Flx_deg1_root(g, p));
    return;
  }
  const e = (p - 1n) / 2n;
  for (;;) {
    const a = randomFl(p);
    let h = FpXQ_pow([a, 1n], e, g, p);
    h = FpX_sub(h, [1n], p);
    const c = FpX_normalize(FpX_gcd(g, h, p), p);
    const dc = FpX_degree(c);
    if (dc > 0 && dc < d) {
      Flx_split_linear(c, p, out);
      Flx_split_linear(FpX_normalize(FpX_divrem(g, c, p)[0], p), p, out);
      return;
    }
  }
}

/** PARI `Flx_oneroot_pre`: one root of `f`, or `null` if there is none. */
function Flx_oneroot(f: Flx, p: bigint): bigint | null {
  const r = Flx_roots(f, p);
  return r.length ? r[0]! : null;
}

/** PARI `Flx_nbroots`: number of distinct roots of `f` in `F_p`. */
function Flx_nbroots(f: Flx, p: bigint): number {
  return Flx_roots(f, p).length;
}

/** PARI `Flv_roots_to_pol(rts, p, 0)`: `prod_i (X - rts[i])`. */
function Flv_roots_to_pol(rts: readonly bigint[], p: bigint): Flx {
  let f: Flx = [1n % p];
  for (const r of rts) f = FpX_mul(f, [(p - r) % p, 1n], p);
  return f;
}

/** PARI `Flv_Flm_polint(x, M, p, 0)`: Lagrange interpolation, column by column. */
function Flv_Flm_polint(x: readonly bigint[], M: Flm, p: bigint): Flx[] {
  const n = x.length;
  /* Z(X) = prod (X - x_i) */
  const Z = Flv_roots_to_pol(x, p);
  /* basis_i = Z / (X - x_i) / Z'(x_i) */
  const basis: Flx[] = [];
  for (let i = 0; i < n; i++) {
    const { q } = Flx_div_by_X_x(Z, x[i]!, p);
    let d = 1n;
    for (let k = 0; k < n; k++) if (k !== i) d = Fp_mul(d, (x[i]! - x[k]! + p) % p, p);
    const di = Fp_inv(d, p);
    basis.push(q.map((c) => Fp_mul(c, di, p)));
  }
  const out: Flx[] = [];
  for (const col of M) {
    const acc = new Array<bigint>(n).fill(0n);
    for (let i = 0; i < n; i++) {
      const y = col[i]! % p;
      if (y === 0n) continue;
      const b = basis[i]!;
      for (let k = 0; k < b.length; k++) acc[k] = (acc[k]! + Fp_mul(b[k]!, y, p)) % p;
    }
    out.push(FpX_renormalize(acc));
  }
  return out;
}

/* ------------------------------------------------------------------ */
/* Deterministic PRNG                                                  */
/* ------------------------------------------------------------------ */

/*
 * DEVIATION: PARI uses its own global `pari_rand()`. Since every use of
 * randomness here is Las Vegas (the answer is independent of the random
 * choices), we use a deterministic xorshift so that runs are reproducible.
 */
let rngState = 0x2545f491n;
function nextRand(): bigint {
  rngState ^= (rngState << 13n) & 0xffffffffffffffffn;
  rngState ^= rngState >> 7n;
  rngState ^= (rngState << 17n) & 0xffffffffffffffffn;
  return rngState;
}
/** PARI `random_Fl(p)`. */
function randomFl(p: bigint): bigint {
  let r = 0n;
  for (let i = 0; i < 4; i++) r = (r << 64n) | nextRand();
  return r % p;
}

/* ------------------------------------------------------------------ */
/* Elliptic curves over F_p in short Weierstrass form (`Fle.c`)        */
/* ------------------------------------------------------------------ */

/** An affine point, or `null` for the point at infinity. */
type Fle = { x: bigint; y: bigint } | null;

function Fle_dbl(P: Fle, a4: bigint, p: bigint): Fle {
  if (!P || P.y === 0n) return null;
  const s = Fp_mul((3n * Fp_mul(P.x, P.x, p) + a4) % p, Fp_inv((2n * P.y) % p, p), p);
  const x = (Fp_mul(s, s, p) - 2n * P.x + 2n * p) % p;
  const y = (Fp_mul(s, (P.x - x + p) % p, p) - P.y + p) % p;
  return { x, y };
}

function Fle_add(P: Fle, Q: Fle, a4: bigint, p: bigint): Fle {
  if (!P) return Q;
  if (!Q) return P;
  if (P.x === Q.x) {
    if (P.y === Q.y) return Fle_dbl(P, a4, p);
    return null;
  }
  const s = Fp_mul((P.y - Q.y + p) % p, Fp_inv((P.x - Q.x + p) % p, p), p);
  const x = (Fp_mul(s, s, p) - P.x - Q.x + 2n * p) % p;
  const y = (Fp_mul(s, (P.x - x + p) % p, p) - P.y + p) % p;
  return { x, y };
}

function Fle_mulu(P: Fle, n: bigint, a4: bigint, p: bigint): Fle {
  if (n === 0n) return null;
  let R: Fle = null;
  let Q = P;
  let m = n;
  while (m > 0n) {
    if (m & 1n) R = Fle_add(R, Q, a4, p);
    Q = Fle_dbl(Q, a4, p);
    m >>= 1n;
  }
  return R;
}

/** PARI `random_Fle_pre(a4, a6, p, pi)`. */
function random_Fle(a4: bigint, a6: bigint, p: bigint): Fle {
  for (;;) {
    const x = randomFl(p);
    const rhs = (Fp_mul(x, (Fp_mul(x, x, p) + a4) % p, p) + a6) % p;
    const y = Fp_sqrt(rhs, p);
    if (y !== null) return { x, y };
  }
}

/** PARI `Fl_elltrace(a4, a6, p)` (`Fle.c`), naive Legendre-symbol count. */
function Fl_elltrace(a4: bigint, a6: bigint, p: bigint): bigint {
  let s = 0n;
  for (let x = 0n; x < p; x++) {
    const rhs = (Fp_mul(x, (Fp_mul(x, x, p) + a4) % p, p) + a6) % p;
    s += BigInt(kronecker(rhs, p));
  }
  return -s;
}

/** PARI `Fl_ellj_pre` (`Fle.c:504-516`). */
function Fl_ellj(a4: bigint, a6: bigint, p: bigint): bigint {
  const a43 = (4n * Fp_mul(a4, Fp_mul(a4, a4, p), p)) % p;
  const a62 = Fp_mul(Fp_mul(a6, a6, p), 27n % p, p);
  const z1 = Fp_mul(a43, 1728n % p, p);
  const z2 = (a43 + a62) % p;
  return Fp_mul(z1, Fp_inv(z2, p), p);
}

/** PARI `Fl_ellj_to_a4a6` (`Fle.c:518-532`). */
function Fl_ellj_to_a4a6(j: bigint, p: bigint): { a4: bigint; a6: bigint } {
  const zagier = 1728n % p;
  if (j === 0n) return { a4: 0n, a6: 1n % p };
  if (j === zagier) return { a4: 1n % p, a6: 0n };
  const k = (zagier - j + p) % p;
  const kj = Fp_mul(k, j, p);
  const k2j = Fp_mul(kj, k, p);
  return { a4: (3n * kj) % p, a6: (2n * k2j) % p };
}

/** PARI `Fl_elltwist_disc` (`Fle.c:558-564`). */
function Fl_elltwist_disc(
  a4: bigint,
  a6: bigint,
  D: bigint,
  p: bigint
): { a4: bigint; a6: bigint } {
  const D2 = Fp_mul(D, D, p);
  return { a4: Fp_mul(a4, D2, p), a6: Fp_mul(a6, Fp_mul(D, D2, p), p) };
}

/**
 * PARI `Flj_order_ufact` (`Fle.c:222-241`): order of `P` on `y^2 = x^3+a4x+a6`
 * given that it divides `n`, whose factorisation is `fa`. Returns 0 if `n P`
 * is not the point at infinity.
 */
function Fle_order_ufact(
  P: Fle,
  n: bigint,
  fa: Array<[bigint, bigint]>,
  a4: bigint,
  p: bigint
): bigint {
  let res = 1n;
  for (const [t, e] of fa) {
    let b = P;
    const te = t ** e;
    if (fa.length !== 1) b = Fle_mulu(b, n / te, a4, p);
    let j = 0n;
    for (; j < e && b !== null; j++) b = Fle_mulu(b, t, a4, p);
    if (b !== null) return 0n;
    res *= t ** j;
  }
  return res;
}

/* ------------------------------------------------------------------ */
/* Volcano traversal (`reference/pari/src/basemath/volcano.c`)         */
/* ------------------------------------------------------------------ */

/** `volcano.c:34-36` */
function is_j_exceptional(j: bigint, p: bigint): boolean {
  return j === 0n || j === 1728n % p;
}

/** `volcano.c:38-44` */
function node_degree(phi: Flm, L: number, j: bigint, p: bigint): number {
  return Flx_nbroots(Flm_Fl_polmodular_evalx(phi, L, j, p), p);
}

/** `volcano.c:52-62`: `Phi_L(X, path[d]) / (X - path[d-1])`. */
function nhbr_polynomial(path: bigint[], d: number, phi: Flm, p: bigint, L: number): Flx {
  const modpol = Flm_Fl_polmodular_evalx(phi, L, path[d]!, p);
  const { q, rem } = Flx_div_by_X_x(modpol, path[d - 1]!, p);
  if (rem !== 0n) throw new PariBugError('nhbr_polynomial: invalid preceding j');
  return q;
}

/** `volcano.c:70-84` */
function extend_path(path: bigint[], phi: Flm, p: bigint, L: number, max_len: number): number {
  let d = 1;
  for (; d < max_len; d++) {
    const nhbr = Flx_oneroot(nhbr_polynomial(path, d, phi, p, L), p);
    if (nhbr === null) break; /* no root: we are on the floor */
    path[d + 1] = nhbr;
  }
  return d;
}

/** PARI `ascend_volcano` (`volcano.c:86-137`), Sutherland 2009 Algorithm Ascend. */
function ascend_volcano(
  phi: Flm,
  j: bigint,
  p: bigint,
  level: number,
  L: number,
  depth: number,
  steps: number
): bigint {
  const path: bigint[] = new Array(depth + 2).fill(0n);
  let max_len = depth - level;
  let first_iter = true;
  if (steps <= 0 || max_len < 0) throw new PariBugError('ascend_volcano: bad params');
  while (steps--) {
    const nhbr_pol = first_iter
      ? Flm_Fl_polmodular_evalx(phi, L, j, p)
      : nhbr_polynomial(path, 1, phi, p, L);
    const nhbrs = Flx_roots(nhbr_pol, p);
    const nhbrs_len = nhbrs.length;
    path[0] = j;
    first_iter = false;

    j = nhbrs[nhbrs_len - 1]!;
    for (let i = 0; i < nhbrs_len - 1; i++) {
      const next_j = nhbrs[i]!;
      if (is_j_exceptional(next_j, p)) {
        if (steps) throw new PariBugError('ascend_volcano: Got to the top with more steps to go!');
        j = next_j;
        break;
      }
      path[1] = next_j;
      const len = extend_path(path, phi, p, L, max_len);
      const last_j = path[len]!;
      if (len === max_len && (is_j_exceptional(last_j, p) || node_degree(phi, L, last_j, p) > 1)) {
        j = next_j;
        break;
      }
    }
    path[1] = j;
    max_len++;
  }
  return j;
}

/** `volcano.c:139-153` */
function random_distinct_neighbours_of(
  phi: Flm,
  j: bigint,
  p: bigint,
  L: number,
  must_have_two: boolean
): [bigint, bigint | null] {
  let modpol = Flm_Fl_polmodular_evalx(phi, L, j, p);
  const n1 = Flx_oneroot(modpol, p);
  if (n1 === null) throw new PariBugError('random_distinct_neighbours_of [no neighbour]');
  modpol = Flx_div_by_X_x(modpol, n1, p).q;
  const n2 = Flx_oneroot(modpol, p);
  if (must_have_two && n2 === null)
    throw new PariBugError('random_distinct_neighbours_of [single neighbour]');
  return [n1, n2];
}

/** PARI `descend_volcano` (`volcano.c:155-204`), Sutherland 2009 Algorithm Descend. */
function descend_volcano(
  phi: Flm,
  j: bigint,
  p: bigint,
  level: number,
  L: number,
  depth: number,
  steps: number
): bigint {
  if (steps <= 0 || level + steps > depth) throw new PariBugError('descend_volcano');
  const max_len = depth - level;
  const path: bigint[] = new Array(max_len + 2).fill(0n);
  path[0] = j;
  if (!level) {
    const nhbrs = Flx_roots(Flm_Fl_polmodular_evalx(phi, L, j, p), p);
    let i: number;
    for (i = 0; i < 3; i++) {
      path[1] = nhbrs[i]!;
      const len = extend_path(path, phi, p, L, max_len);
      if (len < max_len || node_degree(phi, L, path[len]!, p) === 1) break;
    }
    if (i >= 3) throw new PariBugError('descend_volcano [2]');
  } else {
    const [nhbr1, nhbr2] = random_distinct_neighbours_of(phi, j, p, L, true);
    path[1] = nhbr1;
    const len = extend_path(path, phi, p, L, max_len);
    if (
      len === max_len &&
      (is_j_exceptional(path[len]!, p) || node_degree(phi, L, path[len]!, p) !== 1)
    ) {
      path[1] = nhbr2!;
      extend_path(path, phi, p, L, steps);
    }
  }
  return path[steps]!;
}

/** PARI `j_level_in_volcano` (`volcano.c:206-233`). */
function j_level_in_volcano(phi: Flm, j: bigint, p: bigint, L: number, depth: number): number {
  if (depth === 0 || is_j_exceptional(j, p)) return 0;
  const path1: bigint[] = new Array(depth + 2).fill(0n);
  const path2: bigint[] = new Array(depth + 2).fill(0n);
  path1[0] = path2[0] = j;
  const [n1, n2] = random_distinct_neighbours_of(phi, j, p, L, false);
  path1[1] = n1;
  if (n2 === null) return depth;
  path2[1] = n2;
  const path1_len = extend_path(path1, phi, p, L, depth);
  const path2_len = extend_path(path2, phi, p, L, path1_len);
  return depth - path2_len;
}

/** `volcano.c:244-251` */
function get_nbrs(phi: Flm, L: number, J: bigint, xJ: bigint | null, p: bigint): bigint[] {
  let f = Flm_Fl_polmodular_evalx(phi, L, J, p);
  if (xJ !== null) f = Flx_remove_root(f, xJ, p);
  return Flx_roots(f, p);
}

/**
 * PARI `surface_path` (`volcano.c:262-367`): a path of length `n` along the
 * surface of an `L`-volcano of height `h` starting at the surface node `J`.
 * `W`/`wo` is the output array and its base offset.
 */
function surface_path(
  W: bigint[],
  wo: number,
  n: number,
  phi: Flm,
  L: number,
  h: number,
  J: bigint,
  nJ: bigint | null,
  p: bigint
): number {
  const W0 = J;
  W[wo] = J;
  if (n === 1) return 1;

  const T: bigint[][] = new Array(h + 1);
  let v = get_nbrs(phi, L, J, nJ, p);
  if (nJ !== null) v = [nJ, ...v];
  T[0] = v;
  let k = v.length;

  switch (k) {
    case 0:
      throw new PariBugError('surface_path');
    case 1:
      if (h) throw new PariBugError('surface_path');
      W[wo + 1] = v[0]!;
      if (W[wo + 1] === W[wo]) return 1;
      return 2;
    case 2: {
      if (L === 2) {
        const u = get_nbrs(phi, L, v[0]!, J, p);
        const nn = u.length - (u.includes(J) ? 1 : 0);
        W[wo + 1] = nn === 1 ? v[0]! : v[1]!;
        return 2;
      }
      if (h) throw new PariBugError('surface_path');
      W[wo + 1] = v[0]!;
      for (let w = 2; w < n; w++) {
        v = get_nbrs(phi, L, W[wo + w - 1]!, W[wo + w - 2]!, p);
        if (v.length !== 1) throw new PariBugError('surface_path');
        W[wo + w] = v[0]!;
        if (W[wo + w] === W0) return w;
      }
      return n;
    }
  }
  if (!h) throw new PariBugError('surface_path');

  let w = 1;
  for (let x = 0; ; x++) {
    W[wo + w] = T[(w - 1) % h]![x]!;
    if (W[wo + w] === W0) return w;
    if (x === k - 1 && w === n - 1) return n;

    for (let j = w; ; ) {
      const vv = get_nbrs(phi, L, W[wo + j]!, W[wo + j - 1]!, p);
      const m = vv.length;
      if (!m) break; /* hit the floor */
      if (m !== L) throw new PariBugError('surface_path');
      T[j % h] = vv;
      W[wo + ++j] = vv[0]!;
      if (j === w + h) {
        ++w;
        if (W[wo + w] === W0) return w;
        x = 0;
        k = L;
      }
      if (w === n) return w;
    }
  }
}

/** PARI `next_surface_nbr` (`volcano.c:369-405`). */
function next_surface_nbr(
  phi: Flm,
  L: number,
  h: number,
  J: bigint,
  pJ: bigint | null,
  p: bigint
): bigint | null {
  const S = get_nbrs(phi, L, J, pJ, p);
  const k = S.length;
  if (!k) return null;
  if (k === 1 || (pJ === null && k === 2)) return S[0]!;
  if (!h) throw new PariBugError('next_surface_nbr');

  const P: bigint[] = new Array(h + 2).fill(0n);
  P[0] = J;
  let i: number;
  let j = 0;
  for (i = 0; i < k; i++) {
    P[1] = S[i]!;
    for (j = 1; j <= h; j++) {
      const T = get_nbrs(phi, L, P[j]!, P[j - 1]!, p);
      if (T.length === 0) break;
      P[j + 1] = T[0]!;
    }
    if (j < h) throw new PariBugError('next_surface_nbr');
    if (j > h) break;
  }
  if (i === k) throw new PariBugError('next_surf_nbr');
  return S[i]!;
}

/** PARI `common_nbr` (`volcano.c:407-428`). */
function common_nbr(
  J1: bigint,
  Phi1: Flm,
  L1: number,
  J2: bigint,
  Phi2: Flm,
  L2: number,
  p: bigint
): bigint[] {
  const g = Flm_Fl_polmodular_evalx(Phi1, L1, J1, p);
  const f = Flm_Fl_polmodular_evalx(Phi2, L2, J2, p);
  const d = FpX_normalize(FpX_gcd(f, g, p), p);
  if (FpX_degree(d) === 1) return [Flx_deg1_root(d, p)];
  if (FpX_degree(d) !== 2) throw new PariBugError('common_neighbour');
  const r = Flx_roots(d, p);
  if (!r.length) throw new PariBugError('common_neighbour');
  /* PARI returns `rlen`, the number of DISTINCT roots of the degree-2 gcd (so
   * 1 for a double root, 2 otherwise), and sets nbr[0] = r[1], nbr[1] =
   * r[rlen].  Callers branch on that count being 2, so the returned array must
   * have length `rlen`: returning [r0, r0] for a double root made every caller
   * take the ambiguous two-candidate branch.  In `surface_parallel_path` with
   * n == 2 that branch is an unconditional failure, which made `polclass0`
   * reject every j-invariant it drew (e.g. D = -288). */
  return r.length === 1 ? [r[0]!] : [r[0]!, r[1]!];
}

/** PARI `common_nbr_pred_poly` (`volcano.c:430-442`). */
function common_nbr_pred_poly(
  J1: bigint,
  Phi1: Flm,
  L1: number,
  J2: bigint,
  Phi2: Flm,
  L2: number,
  J0: bigint,
  p: bigint
): Flx {
  let g = Flm_Fl_polmodular_evalx(Phi1, L1, J1, p);
  g = Flx_remove_root(g, J0, p);
  const f = Flm_Fl_polmodular_evalx(Phi2, L2, J2, p);
  return FpX_normalize(FpX_gcd(f, g, p), p);
}

/** PARI `common_nbr_pred` (`volcano.c:444-456`). */
function common_nbr_pred(
  J1: bigint,
  Phi1: Flm,
  L1: number,
  J2: bigint,
  Phi2: Flm,
  L2: number,
  J0: bigint,
  p: bigint
): bigint | null {
  const d = common_nbr_pred_poly(J1, Phi1, L1, J2, Phi2, L2, J0, p);
  return FpX_degree(d) === 1 ? Flx_deg1_root(d, p) : null;
}

/** PARI `common_nbr_verify` (`volcano.c:458-470`). */
function common_nbr_verify(
  J1: bigint,
  Phi1: Flm,
  L1: number,
  J2: bigint,
  Phi2: Flm,
  L2: number,
  J0: bigint,
  p: bigint
): bigint | null {
  const d = common_nbr_pred_poly(J1, Phi1, L1, J2, Phi2, L2, J0, p);
  if (FpX_degree(d) <= 0) return null;
  if (FpX_degree(d) > 1) throw new PariBugError('common_neighbour_verify');
  return Flx_deg1_root(d, p);
}

/** PARI `Flm_Fl_polmodular_evalxy` (`volcano.c:472-478`). */
function Flm_Fl_polmodular_evalxy(Phi: Flm, L: number, x: bigint, y: bigint, p: bigint): bigint {
  return Flx_eval(Flm_Fl_polmodular_evalx(Phi, L, x, p), y, p);
}

/** PARI `common_nbr_corner` (`volcano.c:480-501`). */
function common_nbr_corner(
  J1: bigint,
  Phi1: Flm,
  L1: number,
  h1: number,
  J2: bigint,
  Phi2: Flm,
  L2: number,
  J0: bigint,
  p: bigint
): bigint | null {
  const nbrs = common_nbr(J1, Phi1, L1, J2, Phi2, L2, p);
  if (nbrs.length === 2) {
    const nJ2 = next_surface_nbr(Phi1, L1, h1, J2, J0, p);
    if (nJ2 === null) return null;
    let nJ1 = next_surface_nbr(Phi1, L1, h1, nbrs[0]!, J1, p);
    if (nJ1 === null) return null;
    if (Flm_Fl_polmodular_evalxy(Phi2, L2, nJ1, nJ2, p) !== 0n) {
      nbrs[0] = nbrs[1]!;
    } else {
      nJ1 = next_surface_nbr(Phi1, L1, h1, nbrs[1]!, J1, p);
      if (nJ1 === null) return null;
      if (Flm_Fl_polmodular_evalxy(Phi2, L2, nJ1, nJ2, p) === 0n) return null;
    }
  }
  return nbrs[0]!;
}

/** PARI `surface_gcd_cycle` (`volcano.c:503-559`). */
function surface_gcd_cycle(
  W: bigint[],
  wo: number,
  V: bigint[],
  vo: number,
  n: number,
  Phi1: Flm,
  L1: number,
  Phi2: Flm,
  L2: number,
  e: number,
  p: bigint,
  same: boolean
): number {
  let i1 = 0;
  let j2 = 0;
  let i2 = e - 1;
  let j1 = e - 1;
  if (!same) {
    i1 = j1 + 1;
    i2 = n - 1;
  }
  do {
    let f = Flm_Fl_polmodular_evalx(Phi2, L2, V[vo + i1]!, p);
    let g = Flm_Fl_polmodular_evalx(Phi1, L1, W[wo + j1]!, p);
    g = Flx_remove_root(g, W[wo + j1 - 1]!, p);
    const h1 = FpX_normalize(FpX_gcd(f, g, p), p);
    if (FpX_degree(h1) !== 1) break;
    const h11 = h1[1]!;
    const h10 = h1[0]!;

    f = Flm_Fl_polmodular_evalx(Phi2, L2, V[vo + i2]!, p);
    g = Flm_Fl_polmodular_evalx(Phi1, L1, W[wo + j2]!, p);
    let k = j2 + 1;
    if (k === n) k = 0;
    g = Flx_remove_root(g, W[wo + k]!, p);
    const h2 = FpX_normalize(FpX_gcd(f, g, p), p);
    if (FpX_degree(h2) !== 1) break;
    const h21 = h2[1]!;
    const h20 = h2[0]!;

    i1++;
    i2--;
    if (i2 < 0) i2 = n - 1;
    j1++;
    j2--;
    if (j2 < 0) j2 = n - 1;

    const t0 = Fp_mul(h11, h21, p);
    const t1 = Fp_inv(t0, p);
    W[wo + j1] = (p - Fp_mul(Fp_mul(t1, h21, p), h10, p)) % p;
    W[wo + j2] = (p - Fp_mul(Fp_mul(t1, h11, p), h20, p)) % p;
  } while (j2 > j1 + 1);
  return n - j2 + j1 + 1;
}

/** PARI `surface_gcd_path` (`volcano.c:561-586`). */
function surface_gcd_path(
  W: bigint[],
  wo: number,
  V: bigint[],
  vo: number,
  n: number,
  Phi1: Flm,
  L1: number,
  Phi2: Flm,
  L2: number,
  e: number,
  p: bigint,
  same: boolean
): number {
  let i = 0;
  let j = e;
  if (!same) i = j;
  while (j < n) {
    const f = Flm_Fl_polmodular_evalx(Phi2, L2, V[vo + i]!, p);
    let g = Flm_Fl_polmodular_evalx(Phi1, L1, W[wo + j - 1]!, p);
    g = Flx_remove_root(g, W[wo + j - 2]!, p);
    const d = FpX_normalize(FpX_gcd(f, g, p), p);
    if (FpX_degree(d) !== 1) break;
    W[wo + j] = Flx_deg1_root(d, p);
    i++;
    j++;
  }
  return j;
}

/** PARI `surface_parallel_path` (`volcano.c:588-613`). */
function surface_parallel_path(
  W: bigint[],
  wo: number,
  V: bigint[],
  vo: number,
  n: number,
  Phi1: Flm,
  L1: number,
  Phi2: Flm,
  L2: number,
  p: bigint,
  cycle: boolean
): number {
  const nbrs = common_nbr(W[wo]!, Phi1, L1, V[vo + 1]!, Phi2, L2, p);
  if (nbrs.length === 2) {
    if (n <= 2) return 1;
    if (common_nbr_verify(nbrs[0]!, Phi1, L1, V[vo + 2]!, Phi2, L2, W[wo]!, p) === null)
      nbrs[0] = nbrs[1]!;
    else if (common_nbr_verify(nbrs[1]!, Phi1, L1, V[vo + 2]!, Phi2, L2, W[wo]!, p) !== null)
      return 1;
  }
  W[wo + 1] = nbrs[0]!;
  if (n <= 2) return n;
  return cycle
    ? surface_gcd_cycle(W, wo, V, vo, n, Phi1, L1, Phi2, L2, 2, p, false)
    : surface_gcd_path(W, wo, V, vo, n, Phi1, L1, Phi2, L2, 2, p, false);
}

/* ------------------------------------------------------------------ */
/* Norm equations and polycyclic presentations (`polclass.c`)          */
/* ------------------------------------------------------------------ */

/** PARI's `norm_eqn_t` (`paripriv.h`), set by `norm_eqn_set` (`polclass.c:1644`). */
interface NormEqn {
  D: number;
  u: number;
  t: number;
  v: number;
  faw: Array<[bigint, bigint]>;
  p: bigint;
  T: bigint /* a quadratic non-residue mod p */;
}

function norm_eqn_set(
  D: number,
  t: number,
  u: number,
  v: number,
  faw: Array<[bigint, bigint]>,
  p: bigint
): NormEqn {
  let T: bigint;
  do T = randomFl(p);
  while (T === 0n || kronecker(T, p) !== -1);
  return { D, u, t, v, faw, p, T };
}

/** PARI's polycyclic presentation of `cl(D)` (`polclass.c:1166-1282`). */
interface Pcp {
  L: number[] /* generator norms */;
  n: number[] /* relative orders */;
  o: number[] /* absolute orders */;
  m: number[] /* radices */;
  r: number[] /* power relations, `r[ri(i) + j]` */;
  L0: number;
  k: number;
  enum_cnt: number;
  h: number;
  inv: number;
  D: number;
  D0: number;
  u: number;
  fau: Array<[bigint, bigint]>;
}

const evec_ri = (i: number): number => (i * (i - 1)) >> 1;

/** `polclass.c:805-822` */
function evec_reduce(e: number[], n: number[], r: number[], k: number): void {
  if (!k) return;
  for (let i = k - 1; i > 0; i--) {
    if (e[i]! >= n[i]!) {
      const q = Math.floor(e[i]! / n[i]!);
      const ri = evec_ri(i);
      for (let j = 0; j < i; j++) e[j]! += q * r[ri + j]!;
      e[i]! -= q * n[i]!;
    }
  }
  e[0] = ((e[0]! % n[0]!) + n[0]!) % n[0]!;
}

/** `polclass.c:824-833` */
function evec_compose(
  e3: number[],
  e1: number[],
  e2: number[],
  n: number[],
  r: number[],
  k: number
): void {
  for (let i = 0; i < k; i++) e3[i] = e1[i]! + e2[i]!;
  evec_reduce(e3, n, r, k);
}

/** `polclass.c:860-884` */
function evec_inverse(e2: number[], e1: number[], n: number[], r: number[], k: number): void {
  const e3 = e1.slice(0, k);
  const e4 = new Array<number>(k).fill(0);
  for (let i = k - 1; i >= 0; i--)
    if (e3[i]) {
      e4[i]! += n[i]! - e3[i]!;
      evec_reduce(e4, n, r, k);
      e3[i] = n[i]!;
      evec_reduce(e3, n, r, k);
    }
  for (let i = 0; i < k; i++) e2[i] = e4[i]!;
}

/** `polclass.c:900-917` */
function evec_order(e: number[], n: number[], r: number[], k: number): number {
  const f = e.slice(0, k);
  let o = 1;
  for (let i = k - 1; i >= 0; i--)
    if (f[i]) {
      const m = n[i]! / ugcd(f[i]!, n[i]!);
      for (let j = 0; j < k; j++) f[j]! *= m;
      evec_reduce(f, n, r, k);
      o *= m;
    }
  return o;
}

/** `polclass.c:919-934` */
function evec_orders(o: number[], n: number[], r: number[], k: number): void {
  const e = new Array<number>(k).fill(0);
  for (let i = 0; i < k; i++) {
    e[i] = 1;
    if (i) e[i - 1] = 0;
    o[i] = evec_order(e, n, r, k);
  }
}

/** `polclass.c:956-962` */
function evec_n_to_m(m: number[], n: number[], k: number): void {
  m[0] = n[0]!;
  for (let i = 1; i < k; ++i) m[i] = m[i - 1]! * n[i]!;
}

/** `polclass.c:964-972` */
function logfac(n: number): number {
  /* upstream literal `0.57236494292470008707171367567653` (`polclass.c:969`),
   * which is log(pi)/2 rounded to the nearest double */
  const HALFLOGPI = 0.5723649429247001;
  return n * Math.log(n) - n + Math.log(n * (1.0 + 4.0 * n * (1.0 + 2.0 * n))) / 6.0 + HALFLOGPI;
}

/** PARI `upper_bound_on_classpoly_coeffs` (`polclass.c:974-994`), Sutherland Lemma 8. */
function upper_bound_on_classpoly_coeffs(D: number, h: number, qfinorms: number[]): number {
  const C = 2114.567;
  const t = Math.PI * Math.sqrt(-D);
  let B = 0.0;
  let maxak = 0;
  let lnMh = 0;
  for (let k = 1; k <= h; ++k) {
    const ak = qfinorms[k - 1]!;
    const tk = t / ak;
    const lnMk = tk + Math.log(1.0 + C * Math.exp(-tk));
    B += lnMk;
    if (ak > maxak) {
      maxak = ak;
      lnMh = lnMk;
    }
  }
  const m = Math.floor((h + 1) / (Math.exp(lnMh) + 1.0));
  const logbinom = m > 0 && m < h ? logfac(h) - logfac(m) - logfac(h - m) : 0;
  return (B + logbinom - m * lnMh) / Math.LN2 + 2.0;
}

/** PARI `classgp_pcp_check_generators` (`polclass.c:1126-1163`). */
function classgp_pcp_check_generators(n: number[], r: number[], k: number, L0: number): number {
  const s = L0 ? 1 : 0;
  const e1 = new Array<number>(k).fill(0);
  for (let i = s + 1; i < k; i++) {
    if (n[i] !== 2) continue;
    const ri = evec_ri(i);
    let j: number;
    for (j = s; j < i; j++) if (r[ri + j]) break;
    if (j === i) continue;
    for (let i0 = s; i0 < i; i0++) {
      if (4 % n[i0]!) continue;
      e1.fill(0);
      e1[i0] = 4;
      evec_reduce(e1, n, r, k);
      for (j = s; j < i; j++) if (e1[j]) break;
      if (j < i) continue;
      e1.fill(0);
      e1[i0] = 2;
      evec_reduce(e1, n, r, k);
      for (j = s; j < i; j++) if (e1[j] !== r[ri + j]) break;
      if (j === i) return i;
      evec_inverse(e1, e1, n, r, k);
      for (j = s; j < i; j++) if (e1[j] !== r[ri + j]) break;
      if (j === i) return i;
    }
  }
  return -1;
}

/** PARI `next_generator` (`polclass.c:777-795`). */
function next_generator(
  D: number,
  u: number,
  filter: number,
  P: { p: number }
): { gen: Qfb; red: Qfb } {
  let p = P.p;
  for (;;) {
    p = unextprime(p + 1);
    if (kross(D, p) !== -1 && u % p !== 0 && filter % p !== 0) {
      const gen = primeform(BigInt(D), BigInt(p));
      const red = qfbred(gen);
      if (red.a !== 1n) {
        P.p = p;
        return { gen, red };
      }
    }
  }
}

const qfbKey = (q: Qfb): string => `${q.a},${q.b}`;

/** PARI `classgp_make_pcp` (`polclass.c:1165-1282`), Sutherland 2009 Algorithm 2.2. */
function classgp_make_pcp(
  h: number,
  D: number,
  D0: number,
  u: number,
  fau: Array<[bigint, bigint]>,
  inv: number,
  Lfilter: number
): { G: Pcp; height: number } {
  const MAX_GENS = 16;
  const lvl = modinv_level(inv);
  const L = new Array<number>(MAX_GENS).fill(0);
  const m = new Array<number>(MAX_GENS).fill(0);
  const n = new Array<number>(MAX_GENS).fill(0);
  const o = new Array<number>(MAX_GENS).fill(0);
  const r = new Array<number>((MAX_GENS * (MAX_GENS - 1)) / 2).fill(0);

  const G: Pcp = {
    L,
    n,
    o,
    m,
    r,
    L0: 0,
    k: 0,
    enum_cnt: 0,
    h,
    inv,
    D,
    D0,
    u,
    fau,
  };
  if (h === 1) {
    return { G, height: upper_bound_on_classpoly_coeffs(D, h, [1]) };
  }
  /* L0 is only nonzero for ramified double-eta invariants, unsupported here */
  const L0 = 0;
  const enum_cnt = h;
  let GLfilter = (Lfilter * lvl) / ugcd(Lfilter, lvl);
  const DD = BigInt(D);
  let k = 0;
  let T: Qfb[] = [];
  for (;;) {
    k = 0;
    const tbl = new Map<string, number>();
    const ident = primeform(DD, 1n);
    tbl.set(qfbKey(ident), 0);
    T = [ident];
    let nelts = 1;
    const curr = { p: 1 };

    while (nelts < h) {
      if (k === MAX_GENS) throw new PariImplError('classgp_pcp');
      const { gen: gamma_i, red } = next_generator(D, u, GLfilter, curr);
      let beta = red;
      let ri = 1;
      let e = tbl.get(qfbKey(beta));
      let N = T.length;
      while (e === undefined) {
        for (let j = 0; j < N; ++j) {
          const t = qfbcomp(beta, T[j]!);
          tbl.set(qfbKey(t), T.length);
          T.push(t);
        }
        beta = qfbcomp(beta, gamma_i);
        ++ri;
        e = tbl.get(qfbKey(beta));
      }
      if (ri > 1) {
        L[k] = curr.p;
        n[k] = ri;
        nelts *= ri;
        let NN = 1;
        const si = e;
        for (let j = 0; j < k; ++j) {
          r[evec_ri(k) + j] = Math.floor(si / NN) % n[j]!;
          NN *= n[j]!;
        }
        ++k;
      }
      N = T.length;
    }

    const i = classgp_pcp_check_generators(n, r, k, L0);
    if (i < 0) {
      G.L0 = L0;
      G.k = k;
      G.enum_cnt = enum_cnt;
      evec_orders(o, n, r, k);
      evec_n_to_m(m, n, k);
      break;
    }
    GLfilter *= L[i]!;
  }
  const v = new Array<number>(h);
  v[0] = 1;
  for (let i = 1; i < h; ++i) v[i] = Number(T[i]!.a);
  const height = upper_bound_on_classpoly_coeffs(D, enum_cnt, v);

  const L1 = L[k - 1]!;
  const L2 = k > 1 ? L[k - 2]! : 1;
  if (2 * (1 + Math.log2(L1) + Math.log2(L2)) >= 64) throw new PariImplError('classgp_pcp');
  return { G, height };
}

/** PARI `enum_roots` (`volcano.c:615-696`). */
function enum_roots(
  J0: bigint,
  ne: NormEqn,
  fdb: Array<ZM | null>,
  G: Pcp,
  vshape: Array<[bigint, bigint]> | null
): bigint[] | null {
  const MAX_HEIGHT = 64;
  const s = G.L0 ? 1 : 0;
  const n = G.n.slice(s);
  const L = G.L.slice(s);
  const o = G.o.slice(s);
  const k = G.k - s;
  const N = G.enum_cnt;
  const p = ne.p;
  if (!k) return [J0];

  const roots: bigint[] = new Array(N + MAX_HEIGHT).fill(0n);
  const shape = vshape ?? factoru(BigInt(ne.v));

  const Phi: Flm[] = [];
  const e = new Array<number>(k).fill(0);
  const off = new Array<number>(k).fill(0);
  const poff = new Array<number>(k).fill(0);
  const h = new Array<number>(k).fill(0);
  for (let i = 0; i < k; ++i) {
    h[i] = 0;
    for (const [q, ee] of shape)
      if (Number(q) === L[i]) {
        h[i] = Number(ee);
        break;
      }
    Phi.push(polmodular_db_getp(fdb, L[i]!, p));
  }

  let t = surface_path(roots, 0, n[0]!, Phi[0]!, L[0]!, h[0]!, J0, null, p);
  if (t < n[0]!) return null;
  if (k === 1) return roots.slice(0, t);

  const M = new Array<number>(k).fill(0);
  M[0] = 1;
  for (let i = 1; i < k; ++i) M[i] = M[i - 1]! * n[i - 1]!;
  let i = 1;
  while (i < k) {
    let j: number;
    for (j = i + 1; j < k && !e[j]; ++j);
    let root: bigint | null;
    if (j < k) {
      if (e[i]) {
        root = common_nbr_pred(
          roots[off[i]!]!,
          Phi[i]!,
          L[i]!,
          roots[t - M[j]!]!,
          Phi[j]!,
          L[j]!,
          roots[poff[i]!]!,
          p
        );
      } else {
        root = common_nbr_corner(
          roots[off[i]!]!,
          Phi[i]!,
          L[i]!,
          h[i]!,
          roots[t - M[j]!]!,
          Phi[j]!,
          L[j]!,
          roots[poff[j]!]!,
          p
        );
      }
    } else {
      root = next_surface_nbr(
        Phi[i]!,
        L[i]!,
        h[i]!,
        roots[off[i]!]!,
        e[i] ? roots[poff[i]!]! : null,
        p
      );
    }
    if (root === null) break;
    roots[t] = root;
    if (roots[t] === roots[0]) break;

    poff[i] = off[i]!;
    off[i] = t;
    e[i]!++;
    for (j = i - 1; j > 0; --j) {
      e[j] = 0;
      off[j] = off[j + 1]!;
    }

    const t0 = surface_parallel_path(
      roots,
      t,
      roots,
      poff[i]!,
      n[0]!,
      Phi[0]!,
      L[0]!,
      Phi[i]!,
      L[i]!,
      p,
      n[0] === o[0]
    );
    if (t0 < n[0]!) break;
    t += n[0]!;
    for (i = 1; i < k && e[i] === n[i]! - 1; i++);
  }
  if (t !== N) return null;
  return roots.slice(0, t);
}

/* ------------------------------------------------------------------ */
/* Class numbers (`mftrace.c`, `quad.c`, `arith2.c`)                   */
/* ------------------------------------------------------------------ */

/** PARI `hclassno6_count` (`mftrace.c:2487-2508`); `D < 0` fundamental, `-D < 500000`. */
function hclassno6_count(D: number): number {
  const d = -D;
  let h = 0;
  let f = false;
  let b = d & 1;
  let b2 = (1 + d) >> 2;
  let a: number;
  if (!b) {
    for (a = 1; a * a < b2; a++) if (b2 % a === 0) h++;
    f = a * a === b2;
    b = 2;
    b2 = (4 + d) >> 2;
  }
  while (b2 * 3 < d) {
    if (b2 % b === 0) h++;
    for (a = b + 1; a * a < b2; a++) if (b2 % a === 0) h += 2;
    if (a * a === b2) h++;
    b += 2;
    b2 = (b * b + d) >> 2;
  }
  if (b2 * 3 === d) return 6 * h + 2;
  if (f) return 6 * h + 3;
  return 6 * h;
}

/** PARI `quadclassnos(D)` for `D < 0` fundamental (via `hclassno6_count`). */
function quadclassnos(D: number): number {
  if (D === -3 || D === -4) return 1;
  if (-D >= 500000)
    throw new NotImplementedError(
      `SAGE_NOT_IMPLEMENTED: quadclassnos(${D}); |D| >= 500000 needs Buchquad`
    );
  return hclassno6_count(D) / 6;
}

/** PARI `coredisc2u_fact(fa, -1, &P, &E)` (`arith2.c:693-717`). */
function coredisc2u_fact(fa: Array<[bigint, bigint]>): {
  D: number;
  P: number[];
  E: number[];
} {
  let D = 1;
  const P: number[] = [];
  const E: number[] = [];
  for (const [pp, ee] of fa) {
    const p = Number(pp);
    let e = Number(ee);
    if (e & 1) D *= p;
    e >>= 1;
    if (e) {
      P.push(p);
      E.push(e);
    }
  }
  if ((D & 3) !== 3) {
    D *= 4;
    if (!--E[0]!) {
      P.shift();
      E.shift();
    }
  }
  return { D, P, E };
}

/** PARI `uquadclassnoF_fact` (`quad.c:686-717`) with `s = -1`. */
function uquadclassnoF_fact(d: number, P: number[], E: number[]): number {
  let H = 1;
  for (let i = 0; i < P.length; i++) {
    const p = P[i]!;
    const e = E[i]!;
    let Dm = p === 2 ? d & 7 : d % p;
    Dm = -Dm;
    const a = kross(Dm, p);
    if (!a) H *= p ** e;
    else {
      H *= p - a;
      if (e >= 2) H *= p ** (e - 1);
    }
  }
  if (P.length === 0) return H;
  if (d === 4) H >>= 1;
  else if (d === 3) H /= 3;
  return H;
}

/** PARI `uhclassnoF_fact` (`quad.c:1066-1079`). */
function uhclassnoF_fact(faF: Array<[number, number]>, D: number): number {
  let t = 1;
  for (const [p, e] of faF) {
    const s = kross(D, p);
    if (e === 1) {
      t *= 1 + p - s;
      continue;
    }
    if (s === 1) {
      t *= p ** e;
      continue;
    }
    /* usumpow(p, e) = 1 + p + ... + p^(e-1) */
    let q = 0;
    for (let i = 0; i < e; i++) q += p ** i;
    t *= 1 + q * (p - s);
  }
  return t;
}

/** PARI `quadnegclassnou` (`polclass.c:1959-1965`). */
function quadnegclassnou(D: number): {
  h: number;
  D0: number;
  P: number[];
  E: number[];
} {
  const d = -D;
  const { D: d0, P, E } = coredisc2u_fact(factoru(BigInt(d)));
  const h = uquadclassnoF_fact(d0, P, E) * quadclassnos(-d0);
  return { h, D0: -d0, P, E };
}

/** PARI `hclassno_wrapper` (`polclass.c:1418-1431`). */
function hclassno_wrapper(h: number, D: number, Faf: Array<[number, number]>): number {
  if (Faf.length === 0) {
    switch (D) {
      case -3:
        return 1 / 3;
      case -4:
        return 1 / 2;
      default:
        return h;
    }
  }
  return h * uhclassnoF_fact(Faf, D);
}

/* ------------------------------------------------------------------ */
/* Prime selection for class polynomials (`polclass.c:1284-1583`)      */
/* ------------------------------------------------------------------ */

const SMALL_PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31];

/*
 * PARI stores `SMOOTH_INTS` and `HURWITZ_RATIO` as literal tables of 1200
 * entries (`polclass.c:1287-1416`); the source documents the GP code that
 * generates them, which is what we run here (identical values, no risk of a
 * transcription slip).
 */
const V_MAX_TABLE = 1200;
const SMOOTH_INTS: number[] = (() => {
  const t = new Array<number>(V_MAX_TABLE + 1).fill(-1);
  t[0] = -1;
  for (let v = 1; v <= V_MAX_TABLE; v++) {
    let m = v;
    let mask = 0;
    for (let i = 0; i < SMALL_PRIMES.length; i++) {
      const p = SMALL_PRIMES[i]!;
      if (m % p === 0) {
        mask |= 1 << i;
        while (m % p === 0) m /= p;
      }
    }
    t[v] = m === 1 ? mask : -1;
  }
  return t;
})();

/* ceil(128 * prod_{p | v} (p+1)/(p-1)), or 0 if v is not 31-smooth */
const HURWITZ_RATIO: number[] = (() => {
  const t = new Array<number>(V_MAX_TABLE + 1).fill(0);
  for (let v = 1; v <= V_MAX_TABLE; v++) {
    let m = v;
    let num = 1;
    let den = 1;
    let ok = true;
    for (let d = 2; d <= m; d++) {
      if (m % d !== 0) continue;
      if (d > 31) {
        ok = false;
        break;
      }
      num *= d + 1;
      den *= d - 1;
      while (m % d === 0) m /= d;
    }
    t[v] = ok ? Math.ceil((num / den) * 128 - 1e-9) : 0;
  }
  return t;
})();

/** PARI `factor_uv` (`polclass.c:1433-1451`): `factor(u*v)`. */
function factor_uv(
  fau: Array<[number, number]>,
  v: number,
  vfactors: number
): Array<[number, number]> {
  if (!vfactors) return fau;
  const map = new Map<number, number>(fau);
  let vv = v;
  for (let i = 0; vfactors; i++, vfactors >>= 1)
    if (vfactors & 1) {
      const p = SMALL_PRIMES[i]!;
      let e = 0;
      while (vv % p === 0) {
        vv /= p;
        e++;
      }
      map.set(p, (map.get(p) ?? 0) + e);
      if (vv === 1) break;
    }
  return [...map.entries()].sort((a, b) => a[0] - b[0]);
}

interface PoolEntry {
  p: bigint;
  t: number;
  v: number;
  rho_inv: number;
  vfactors: number;
  faw: Array<[number, number]>;
}

/** PARI `select_classpoly_prime_pool` (`polclass.c:1453-1531`), Sutherland Alg. 2.1. */
function select_classpoly_prime_pool(min_bits: number, delta: number, G: Pcp): PoolEntry[] {
  const V_MAX = 1200;
  const d = -G.D;
  const inv = G.inv;
  const fau: Array<[number, number]> = G.fau.map(([a, b]) => [Number(a), Number(b)]);
  const hurwitz = hclassno_wrapper(G.h, G.D0, fau);
  const res: PoolEntry[] = [];
  const t_min = new Array<number>(V_MAX).fill(2);
  let bits = 0.0;
  const BIL = 62; /* PARI's BITS_IN_LONG - 2 ceiling for products; see below */
  const t_size_lim = 2.0 * Math.sqrt(2 ** BIL - d / 4);

  for (let z = d / (2.0 * hurwitz); ; z *= delta + 1.0) {
    const v_bound_aux = 4.0 * z * hurwitz;
    for (let v = 1; v < V_MAX; v++) {
      const vfactors = SMOOTH_INTS[v]!;
      if (vfactors < 0) continue;
      const hurwitz_ratio_bound = HURWITZ_RATIO[v]! / 128.0;
      const vd = v * d;
      if (vd >= v_bound_aux * hurwitz_ratio_bound) break;
      const v2d = v * vd;
      const faw = factor_uv(fau, v, vfactors);
      const H = hclassno_wrapper(G.h, G.D0, faw);
      const max_p = z * v * hurwitz * hurwitz_ratio_bound;
      const t_max = 2.0 * Math.sqrt(Math.min(2 ** BIL - v2d / 4, max_p));
      let t = t_min[v]!;
      if ((t & 1) !== (v2d & 1)) t++;
      /* 4p = t^2 + v^2 d; p can exceed 2^53, so keep it in BigInt */
      let pp = (BigInt(t) * BigInt(t) + BigInt(v2d)) >> 2n;
      for (; t <= t_max; pp += BigInt(t + 1), t += 2) {
        if (pp > BigInt(Number.MAX_SAFE_INTEGER))
          throw new NotImplementedError(
            'SAGE_NOT_IMPLEMENTED: select_classpoly_prime_pool needs primes > 2^53'
          );
        const pn = Number(pp);
        if (modinv_good_prime(inv, pn) && isPrime(pp)) {
          res.push({ p: pp, t, v, rho_inv: Math.floor(pn / H), vfactors, faw });
          bits += Math.log2(pn);
        }
      }
      t_min[v] = t;
      if (bits > min_bits) return res;
    }
    if (t_min[1]! >= t_size_lim) throw new PariArchError(`class polynomial of discriminant ${G.D}`);
  }
}

/** PARI `height_margin` (`polclass.c:1537-1546`). */
function height_margin(inv: number): number {
  if (inv === INV_F) return 64;
  if (inv === INV_G2) return 5;
  if (inv !== INV_J) return 256;
  return 0;
}

/** PARI `select_classpoly_primes` (`polclass.c:1548-1583`). */
function select_classpoly_primes(
  delta: number,
  G: Pcp,
  height: number
): { primes: PoolEntry[]; vfactors: number; biggest_v: number } {
  const k = 2;
  const s = modinv_height_factor(G.inv);
  const b = height / s + height_margin(G.inv);
  const min_prime_bits = k * b;
  const pool = select_classpoly_prime_pool(min_prime_bits, delta, G);
  pool.sort((x, y) => x.rho_inv - y.rho_inv);
  let prime_bits = 0.0;
  let biggest_v = 0;
  let vfactors = 0;
  let i = 0;
  for (; i < pool.length; i++) {
    const q = pool[i]!;
    vfactors |= q.vfactors;
    prime_bits += Math.log2(Number(q.p));
    if (q.v > biggest_v) biggest_v = q.v;
    if (prime_bits > b) break;
  }
  return { primes: pool.slice(0, Math.min(i + 1, pool.length)), vfactors, biggest_v };
}

/* ------------------------------------------------------------------ */
/* Finding a j-invariant with a prescribed trace (`polclass.c:55-773`) */
/* ------------------------------------------------------------------ */

/** PARI `hasse_bounds` (`polclass.c:58-60`). */
function hasse_bounds(p: bigint): { low: bigint; high: bigint } {
  const u = isqrt(4n * p);
  return { low: p + 1n - u, high: p + 1n + u };
}

function famatsmall_divexact(
  a: Array<[bigint, bigint]>,
  b: Array<[bigint, bigint]>
): Array<[bigint, bigint]> {
  const map = new Map<bigint, bigint>(a);
  for (const [q, e] of b) {
    const cur = map.get(q);
    if (cur === undefined) continue;
    if (cur === e) map.delete(q);
    else map.set(q, cur - e);
  }
  return [...map.entries()].sort((x, y) => (x[0] < y[0] ? -1 : 1));
}

/** PARI `test_curve_order` (`polclass.c:82-131`), Sutherland 2009 TestCurveOrder. */
function test_curve_order(
  ne: NormEqn,
  a4in: bigint,
  a6in: bigint,
  N0in: bigint,
  N1in: bigint,
  n0in: Array<[bigint, bigint]>,
  n1in: Array<[bigint, bigint]>,
  hasse_low: bigint,
  hasse_high: bigint
): boolean {
  const p = ne.p;
  let a4 = a4in;
  let a6 = a6in;
  let N0 = N0in;
  let N1 = N1in;
  let n0 = n0in;
  let n1 = n1in;
  let m0 = 1n;
  let m1 = 1n;
  let swapped = false;
  let first = true;
  let a4t = 0n;
  let a6t = 0n;
  if (p <= 11n) {
    const card = p + 1n - Fl_elltrace(a4, a6, p);
    return card === N0 || card === N1;
  }
  for (;;) {
    let n_s: bigint;
    let Q = random_Fle(a4, a6, p);
    if (m0 === 1n) n_s = Fle_order_ufact(Q, N0, n0, a4, p);
    else if (N0 % m0 !== 0n) n_s = 0n;
    else {
      const fa0 = famatsmall_divexact(n0, factoru(m0));
      Q = Fle_mulu(Q, m0, a4, p);
      n_s = Fle_order_ufact(Q, N0 / m0, fa0, a4, p);
    }
    if (n_s === 0n) {
      if (swapped || N1 % m0 !== 0n || N0 % m1 !== 0n) return false;
      [n0, n1] = [n1, n0];
      [N0, N1] = [N1, N0];
      swapped = true;
      continue;
    }
    m0 *= n_s;
    const a1 = (2n * p + 2n) % m1;
    let x = ((hasse_low + m0 - 1n) / m0) * m0;
    for (; x <= hasse_high; x += m0) if (x % m1 === a1 && x !== N0 && x !== N1) break;
    if (x > hasse_high) return true;
    if (first) {
      ({ a4: a4t, a6: a6t } = Fl_elltwist_disc(a4, a6, ne.T, p));
      first = false;
    }
    [a4, a4t] = [a4t, a4];
    [a6, a6t] = [a6t, a6];
    [m0, m1] = [m1, m0];
  }
}

/**
 * PARI `find_j_inv_with_given_trace` (`polclass.c:710-773`), Sutherland 2009
 * Algorithm 1.1.
 *
 * DEVIATION: upstream picks a torsion constraint `m > 1` from Sutherland's
 * `torcosts.h` tables and draws curves from the corresponding `X_1(m)`
 * parameterisation (`reference/pari/src/basemath/crvwtors.c`). We use the
 * `m = 1`, `twist = 3` entry of those tables, i.e. plain uniformly random
 * curves and the two-sided test `(p+1)P == tP`. This is one of the choices the
 * upstream tables can return; it costs a constant factor in speed and cannot
 * change the result, which is a j-invariant of trace +/- t.
 */
function find_j_inv_with_given_trace(ne: NormEqn): bigint {
  const p = ne.p;
  const t = BigInt(ne.t);
  if (p === 2n || p === 3n) {
    if (t === 0n) throw new PariBugError('find_j_inv_with_given_trace');
    return t;
  }
  const p1 = p + 1n;
  const N0 = p1 - t;
  const N1 = p1 + t;
  const n0 = factoru(N0);
  const n1 = factoru(N1);
  const { low, high } = hasse_bounds(p);
  for (;;) {
    const a4 = randomFl(p);
    const a6 = randomFl(p);
    if (a4 === 0n || a6 === 0n) continue;
    if ((4n * Fp_mul(a4, Fp_mul(a4, a4, p), p) + 27n * Fp_mul(a6, a6, p)) % p === 0n) continue;
    const P = random_Fle(a4, a6, p);
    const Pp1 = Fle_mulu(P, p1, a4, p);
    const Pt = Fle_mulu(P, t, a4, p);
    const eq = (Pp1 === null && Pt === null) || (Pp1 !== null && Pt !== null && Pp1.x === Pt.x);
    if (!eq) continue;
    if (test_curve_order(ne, a4, a6, N0, N1, n0, n1, low, high)) return Fl_ellj(a4, a6, p);
  }
}

/** PARI `oneroot_of_classpoly` (`polclass.c:1585-1633`), Sutherland Alg. 1.2. */
function polclass_oneroot_of_classpoly(
  j0: bigint,
  ne: NormEqn,
  jdb: Array<ZM | null>
): { j: bigint; cert: boolean; ok: boolean } {
  const p = ne.p;
  let j = j0;
  if (j === 0n || j === 1728n % p) throw new PariBugError('oneroot_of_classpoly');
  let cert = true;
  const factors = ne.faw;
  if (factors.length === 0) return { j, cert, ok: true };
  const L_bound = Math.max(Math.log(-ne.D), ne.v);
  for (const [Lb, depthb] of factors) {
    const L = Number(Lb);
    const depth = Number(depthb);
    if (L > L_bound) {
      cert = false;
      break;
    }
    const phi = polmodular_db_getp(jdb, L, p);
    const jlvl = j_level_in_volcano(phi, j, p, L, depth);
    const lvl_diff = z_lval(ne.u, L) - jlvl;
    if (lvl_diff < 0) j = ascend_volcano(phi, j, p, jlvl, L, depth, -lvl_diff);
    else if (lvl_diff > 0) j = descend_volcano(phi, j, p, jlvl, L, depth, lvl_diff);
  }
  return { j, cert, ok: j !== 0n && j !== 1728n % p };
}

/** PARI `modfn_root` (`polmodular.c:800-824`), for the invariants we support. */
function modfn_root(j: bigint, ne: NormEqn, inv: number): bigint {
  switch (inv) {
    case INV_J:
      return j;
    case INV_G2:
      /* Fl_sqrtl(j, 3, p): p = 2 mod 3 here (modinv_good_prime), so the cube
       * root is unique and equal to j^((2p-1)/3). */
      return Fp_pow(j, (2n * ne.p - 1n) / 3n, ne.p);
  }
  throw new NotImplementedError(
    `SAGE_NOT_IMPLEMENTED: modfn_root for invariant ${inv} (Weber/double-eta` +
      ' class invariants are not ported)'
  );
}

/** PARI `modfn_preimage` (`polmodular.c:850-864`), for the invariants we support. */
function modfn_preimage(x: bigint, p: bigint, inv: number): bigint {
  switch (inv) {
    case INV_J:
      return x;
    case INV_G2:
      return Fp_pow(x, 3n, p);
  }
  throw new NotImplementedError(`SAGE_NOT_IMPLEMENTED: modfn_preimage for invariant ${inv}`);
}

/**
 * PARI `find_jinv` (`polclass.c:1716-1746`).
 *
 * Returns the root together with `oneroot_of_classpoly`'s endomorphism-ring
 * certificate (`*cert` upstream): `cert` is false exactly when the walk up the
 * volcano was cut short because a prime dividing `w` exceeded `L_bound`, so
 * that `j` is only *probably* on the right level.
 */
function find_jinv(
  ne: NormEqn,
  inv: number,
  jdb: Array<ZM | null>
): { j: bigint; cert: boolean } {
  for (;;) {
    const j_t = find_j_inv_with_given_trace(ne);
    if (j_t === 0n)
      throw new PariBugError("polclass0: Couldn't find j-invariant with given trace.");
    const { j, cert, ok } = polclass_oneroot_of_classpoly(j_t, ne, jdb);
    if (ok) return { j: modfn_root(j, ne, inv), cert };
  }
}

/**
 * PARI `vecsmall_isin_skip(v, x, k)` (`volcano.c:1636-1643`): index of the
 * first `i >= k` with `v[i] == x`, or 0.  Zero-based here: `k` is the first
 * index examined and the result is `-1` when `x` does not occur.
 */
function vecsmall_isin_skip(v: readonly bigint[], x: bigint, k: number): number {
  for (let i = k; i < v.length; ++i) if (v[i] === x) return i;
  return -1;
}

/** PARI `polclass_roots_modp` (`polclass.c:1748-1777`). */
function polclass_roots_modp(ne: NormEqn, G: Pcp, db: PolmodularDB): bigint[] {
  const inv = G.inv;
  const jdb = polmodular_db_for_inv(db, INV_J);
  const fdb = polmodular_db_for_inv(db, inv);
  const vshape = factoru(BigInt(ne.v));
  for (;;) {
    const { j, cert } = find_jinv(ne, inv, jdb);
    const res = enum_roots(j, ne, fdb, G, vshape);
    /* `enum_roots` only fails when j has the wrong endomorphism ring, which
     * `oneroot_of_classpoly` rules out when it returns a certificate.  PARI
     * raises e_BUG here rather than looping forever. */
    if (!res && cert) throw new PariBugError('polclass_roots_modp');
    /* Without a certificate, a repeat of roots[0] later in the cycle means j
     * was on the wrong level after all: draw again. */
    if (res && !cert && vecsmall_isin_skip(res, res[0]!, 1) >= 0) continue;
    if (res) return res;
  }
}

/* ------------------------------------------------------------------ */
/* Chinese remaindering (`ncV_chinese_center`, `nmV_chinese_center`)   */
/* ------------------------------------------------------------------ */

/**
 * PARI `ncV_chinese_center(V, plist, NULL)`: CRT the vectors `V[i] mod
 * plist[i]` into a single integer vector with centred residues.
 */
function ncV_chinese_center(V: bigint[][], plist: bigint[]): bigint[] {
  const len = V[0]!.length;
  let P = 1n;
  for (const p of plist) P *= p;
  const half = P >> 1n;
  const coef: bigint[] = plist.map((p) => {
    const Q = P / p;
    return Fp_mul(Q, Fp_inv(Q % p, p), P);
  });
  const out = new Array<bigint>(len).fill(0n);
  for (let k = 0; k < len; k++) {
    let x = 0n;
    for (let i = 0; i < V.length; i++) x = (x + Fp_mul(coef[i]!, V[i]![k]!, P)) % P;
    out[k] = x > half ? x - P : x;
  }
  return out;
}

/** As above for a list of matrices (PARI `nmV_chinese_center`). */
function nmV_chinese_center(V: Flm[], plist: bigint[]): ZM {
  const ncols = V[0]!.length;
  const out: ZM = [];
  for (let c = 0; c < ncols; c++)
    out.push(
      ncV_chinese_center(
        V.map((m) => m[c]!),
        plist
      )
    );
  return out;
}

/* ------------------------------------------------------------------ */
/* Hilbert class polynomials (`polclass.c:1980-2108`)                  */
/* ------------------------------------------------------------------ */

/** PARI `polclass_small_disc` (`polclass.c:1944-1957`). */
function polclass_small_disc(D: number, inv: number): Flx | null {
  if (D === -3) return [0n, 1n];
  if (D === -4) {
    switch (inv) {
      case INV_J:
        return [-1728n, 1n];
      case INV_G2:
        return [-12n, 1n];
      default:
        throw new PariBugError('polclass_small_disc');
    }
  }
  return null;
}

/**
 * PARI `polclass0` (`polclass.c:1980-2108`): the class polynomial
 * `H_D(X)` for the class invariant `inv`, as a coefficient list.
 */
export function polclass0(D: number, inv: number, db: PolmodularDB): bigint[] {
  const small = polclass_small_disc(D, inv);
  if (small) return small;
  if (modinv_is_double_eta(inv) || modinv_is_Weber(inv))
    throw new NotImplementedError(
      `SAGE_NOT_IMPLEMENTED: polclass0 for invariant ${inv}; only INV_J and INV_G2 ` +
        'are ported (Weber/double-eta need double_eta_raw + adjust_signs/orientation)'
    );

  const delta = 0.5;
  const filter = 1;
  const { h, D0, P: Pu, E: Eu } = quadnegclassnou(D);
  const u = D === D0 ? 1 : Math.round(Math.sqrt(D / D0));
  const fau: Array<[bigint, bigint]> = Pu.map((p, i) => [BigInt(p), BigInt(Eu[i]!)]);

  const { G, height } = classgp_make_pcp(h, D, D0, u, fau, inv, filter);
  const { primes, vfactors, biggest_v } = select_classpoly_primes(delta, G, height);

  /* Prepopulate db with all the modular polynomials we might need */
  if (u > 1) {
    const maxL = Math.max(Math.log(-D), biggest_v);
    for (const L of Pu) {
      if (L > maxL) break;
      polmodular_db_add_level(db, L, INV_J);
    }
  }
  let vf = vfactors;
  for (let i = 0; vf; ++i) {
    if (vf & 1) polmodular_db_add_level(db, SMALL_PRIMES[i]!, INV_J);
    vf >>= 1;
  }
  const s = G.L0 ? 1 : 0;
  polmodular_db_add_levels(db, G.L.slice(s, G.k), inv);

  const H: bigint[][] = [];
  const plist: bigint[] = [];
  for (const q of primes) {
    const faw: Array<[bigint, bigint]> = q.faw.map(([a, b]) => [BigInt(a), BigInt(b)]);
    const ne = norm_eqn_set(D, q.t, u, q.v, faw, q.p);
    const roots = polclass_roots_modp(ne, G, db);
    H.push(Flv_roots_to_pol(roots, q.p));
    plist.push(q.p);
  }
  /* pad all the reduced polynomials to the same length (they are monic of the
   * same degree, but `FpX` strips trailing zeros) */
  let deg = 0;
  for (const v of H) deg = Math.max(deg, v.length);
  for (const v of H) while (v.length < deg) v.push(0n);
  return ncV_chinese_center(H, plist);
}

/* ------------------------------------------------------------------ */
/* Modular polynomials by CM (`polmodular.c:1035-1935, 3663-4595`)     */
/* ------------------------------------------------------------------ */

/** PARI `modpoly_height_bound` (`polmodular.c:3675-3710`). */
function modpoly_height_bound(L: number, inv: number): number {
  let nbits = 6.0 * L * Math.log2(L) + (16 / Math.LN2) * L + 8.0 * Math.sqrt(L) * Math.log2(L);
  const nbits2a = 6.0 * L * Math.log2(L) + (17 / Math.LN2) * L;
  if (nbits2a < nbits) nbits = nbits2a;
  const hf = modinv_height_factor(inv);
  if (hf > 1) {
    const nbits2 = nbits - 4.01 * L - 3.0;
    nbits = nbits2 / hf + 4.01 * L + 3.0;
  }
  if (inv === INV_F) {
    let c: number;
    if (L < 30) c = 45;
    else if (L < 100) c = 36;
    else if (L < 300) c = 32;
    else if (L < 600) c = 26;
    else if (L < 1200) c = 24;
    else if (L < 2400) c = 22;
    else c = 20;
    nbits = (6.0 * L * Math.log2(L) + c * L) / hf;
  }
  return nbits;
}

const SMOOTH_PRIMES = 31; /* (BITS_IN_LONG >> 1) - 1 for 64-bit longs */
const MAX_L1 = 255;
const MODPOLY_MAX_DCNT = 64;
const MODPOLY_USE_L1 = 1;
const MODPOLY_IGNORE_SPARSE_FACTOR = 8;

interface DiscInfo {
  inv: number;
  L: number;
  D0: number;
  D1: number;
  L0: number;
  L1: number;
  n1: number;
  n2: number;
  dl1: number;
  dl2_0: number;
  dl2_1: number;
  nprimes: number;
  cost: number;
  bits: number;
  primes: bigint[];
  traces: number[];
}

/** Is the reduced form `Q` the principal (identity) class? */
function qfb_is_1(Q: Qfb): boolean {
  return Q.a === 1n;
}

/**
 * PARI `rec_order` (`bb_group.c:668-697`). Note that, exactly as upstream, the
 * per-prime loop is bounded by the exponent: if `a^o != 1` the result is
 * garbage rather than an infinite loop.
 */
function rec_order(a: Qfb, o: bigint, m: Array<[bigint, bigint]>, x: number, y: number): bigint {
  if (qfb_is_1(a)) return 1n;
  if (x === y) {
    const [p, e] = m[x]!;
    let b = a;
    for (let i = 0n; i < e; i++) {
      if (qfb_is_1(b)) return p ** i;
      b = qfbpow(b, p);
    }
    return p ** e;
  }
  const z = (x + y) >> 1;
  let cof = 1n;
  for (let i = x; i <= z; i++) cof *= m[i]![0] ** m[i]![1];
  let b = qfbpow(a, cof);
  const o1 = rec_order(b, o / cof, m, z + 1, y);
  b = qfbpow(a, o1);
  const o2 = rec_order(b, o / o1, m, x, z);
  return o1 * o2;
}

/**
 * PARI `qfi_order(Q, o)` = `gen_order(Q, o, NULL, &qfi_group)`
 * (`quad.c:582-584`, `bb_group.c:699-713`).
 */
function qfi_order(Q: Qfb, N: bigint): bigint {
  if (N === 1n) return 1n;
  const m = Z_factor(N);
  return rec_order(Q, N, m, 0, m.length - 1);
}

/** PARI `qfi_Shanks(R, Q, n)`: BSGS discrete log of `R` in `<Q>`, `|Q| | n`. */
function qfi_Shanks(R: Qfb, Q: Qfb, n: bigint): bigint | null {
  if (n <= 0n) return null;
  let m = isqrt(n);
  if (m * m < n) m += 1n;
  const tbl = new Map<string, bigint>();
  let cur = qfbred(qfb_1(Q));
  for (let i = 0n; i < m; i++) {
    const key = qfbKey(cur);
    if (!tbl.has(key)) tbl.set(key, i);
    cur = qfbcomp(cur, Q);
  }
  const Qm = qfbinv(qfbpow(Q, m)); /* Q^-m */
  let gamma = qfbred(R);
  for (let i = 0n; i <= m; i++) {
    const j = tbl.get(qfbKey(gamma));
    if (j !== undefined) {
      const x = (i * m + j) % n;
      if (qfb_equal(qfbred(qfbpow(Q, x)), qfbred(R))) return x;
    }
    gamma = qfbcomp(gamma, Qm);
  }
  return null;
}

/** PARI `qform_primeform2` (`polmodular.c:3752-3781`). */
function qform_primeform2(p: number, D: number): Qfb | null {
  const a = BigInt(p) * BigInt(p);
  const Dp2 = a * BigInt(D);
  const M = BigInt(p - 1);
  for (let k = D & 1; k <= p; k += 2) {
    const c = (k * k - D) / 4;
    if (c % p === 0) continue;
    const q = mkqfb(a, BigInt(k * p), BigInt(c), Dp2);
    const Q = qfbred(q);
    const ord = qfi_order(Q, M);
    if (ord === M) {
      if (qfbpow(Q, M).a === 1n) return q;
      break;
    }
  }
  return null;
}

/** PARI `primeform_discrete_log` (`polmodular.c:3783-3794`). */
function primeform_discrete_log(L0: number, L: number, n: number, D: number): number {
  const DD = BigInt(D);
  const Q = primeform(DD, BigInt(L0));
  const R = primeform(DD, BigInt(L));
  const X = qfi_Shanks(R, Q, BigInt(n));
  return X === null ? -1 : Number(X);
}

/** PARI `select_L0` (`polmodular.c:3796-3843`). */
function select_L0(L: number, inv: number, initial_L0: number): number {
  const modinv_N = modinv_level(inv);
  if (modinv_N % L === 0) throw new PariBugError('select_L0');

  if (
    inv === INV_F ||
    inv === INV_F2 ||
    inv === INV_F4 ||
    inv === INV_F8 ||
    inv === INV_W2W3 ||
    inv === INV_W2W3E2 ||
    inv === INV_W3W3
  ) {
    if (L === 19) return 13;
    if (L === 29) return 7;
  }
  if (inv === INV_W2W5 && L === 19) return 13;
  if (inv === INV_W2W5E2 && (L === 7 || L === 19)) return 13;
  if ((inv === INV_W2W7 || inv === INV_W2W7E2) && L === 11) return 13;
  if (inv === INV_W3W5) {
    if (L === 7) return 13;
    if (L === 17) return 7;
  }
  if (inv === INV_W3W7) {
    if (L === 29 || L === 101) return 11;
    if (L === 11 || L === 19) return 13;
  }
  let L0 = unextprime(initial_L0 + 1);
  while (L0 === L || modinv_N % L0 === 0) L0 = unextprime(L0 + 1);
  return L0;
}

/** PARI `primeform_exp_order` (`polmodular.c:3845-3852`). */
function primeform_exp_order(L: number, n: number, D: number, ord: number): number {
  const Q = qfbpow(primeform(BigInt(D), BigInt(L)), BigInt(n));
  return Number(qfi_order(Q, BigInt(ord)));
}

/** PARI `orientation_ambiguity` (`polmodular.c:3854-3901`). */
function orientation_ambiguity(
  D1: number,
  L0: number,
  modinv_p1: number,
  modinv_p2: number
): boolean {
  let ambiguity = false;
  let Q1 = red_primeform(D1, modinv_p1);
  if (!Q1) return true;
  let Q2: Qfb | null = null;
  if (modinv_p2 > 1) {
    if (modinv_p1 === modinv_p2) Q1 = qfbsqr(Q1);
    else {
      let P2 = red_primeform(D1, modinv_p2);
      if (!P2) return true;
      const Q = qfbsqr(P2);
      const R = qfbsqr(Q1);
      if (Q.a === R.a && (Q.b < 0n ? -Q.b : Q.b) === (R.b < 0n ? -R.b : R.b)) {
        ambiguity = true;
      } else {
        Q2 = qfbcomp(Q1, P2);
        P2 = qfbinv(P2);
        Q1 = qfbcomp(Q1, P2);
      }
    }
  }
  if (!ambiguity) {
    const P0 = red_primeform(D1, L0);
    if (!P0) return true;
    const P = qfbsqr(P0);
    if (P.a === Q1.a || (modinv_p2 > 1 && modinv_p1 !== modinv_p2 && Q2 !== null && P.a === Q2.a))
      ambiguity = true;
  }
  return ambiguity;
}

/** PARI `check_generators` (`polmodular.c:3903-3939`). */
function check_generators(
  out: { n1: number; m: number },
  D: number,
  h: number,
  n: number,
  subgrp_sz: number,
  L0: number,
  L1: number
): boolean {
  const m = primeform_exp_order(L0, n, D, h);
  out.m = m;
  const n1 = n * m;
  if (!n1) throw new PariBugError('check_generators');
  out.n1 = n1;
  if (n1 < subgrp_sz / 2 || (!L1 && n1 < subgrp_sz)) return false;
  if (n1 < subgrp_sz && !(n1 & 1)) {
    const D1 = BigInt(D);
    const Q = qfbpow(primeform(D1, BigInt(L0)), BigInt(n1 / 2));
    if (qfb_equal(Q, qfbred(primeform(D1, BigInt(L1))))) return false;
  }
  return true;
}

/** PARI `modpoly_pickD_primes` (`polmodular.c:3941-4062`). */
function modpoly_pickD_primes(
  collect: boolean,
  max: number,
  xprimes: bigint[],
  minbits: number,
  Dinfo: DiscInfo
): { n: number; totbits: number; primes: bigint[]; traces: number[] } {
  const D = Dinfo.D1;
  const absD = -D;
  const L0 = Dinfo.L0;
  const L1 = Dinfo.L1;
  const L = Dinfo.L;
  const inv = Dinfo.inv;
  const FF_BITS = 62; /* BITS_IN_LONG - 2 */
  const primes: bigint[] = [];
  const traces: number[] = [];

  const pfilter = modinv_pfilter(inv);
  if (pfilter & IQ_FILTER_1MOD3 && D % 3 === 0) return { n: 0, totbits: 0, primes, traces };
  if (pfilter & IQ_FILTER_1MOD4 && (D & 0xf) === 0) return { n: 0, totbits: 0, primes, traces };

  let one_prime = false;
  let totbits = 0;
  if (max <= 1) {
    const q = (pfilter & IQ_FILTER_1MOD3 ? 2 : 1) * (pfilter & IQ_FILTER_1MOD4 ? 2 : 1);
    one_prime =
      2 ** Math.floor((FF_BITS + 1) / 2) * (Math.log2(L * L * -D) - 1) >
      q * L * minbits * FF_BITS * Math.LN2;
    if (one_prime) totbits = minbits + 1; /* lie */
  }

  let m = 0;
  let n = 0;
  let bits = 0.0;
  outer: for (let v = 1; v < 100 && bits < minbits; v++) {
    if (ugcd(absD, v) !== 1) continue;
    if (v > 2 && modinv_is_double_eta(inv) && ugcd(modinv_level(inv), v) !== 1) continue;
    if (v & 1 && (D & 7) === 1) continue;
    if (L0 === 2 && !(v & 3)) continue;
    if (pfilter & IQ_FILTER_1MOD4 && (v * v * D) % 16 === 0) continue;
    if (pfilter & IQ_FILTER_1MOD3 && v % 3 === 0) continue;
    if (L0 !== 2 && v % L0 === 0) continue;
    if (L1 && v % L1 === 0) continue;
    if ((v * v * absD) / 4 > 2 ** FF_BITS / (L * L)) break;
    let a1_start: number;
    let a1_delta: number;
    if (v & 1 && D & 1) {
      a1_start = 1;
      a1_delta = 2;
    } else {
      a1_start = (v * v * D) & 7 ? 2 : 0;
      a1_delta = 4;
    }
    for (let a1 = a1_start; bits < minbits; a1 += a1_delta) {
      const A1 = BigInt(a1);
      const a2 = (A1 * A1 + BigInt(v * v * absD)) >> 2n;
      if (a2 % BigInt(L) === 0n) continue;
      const t = a1 * L + 2;
      const p = a2 * BigInt(L * L) + BigInt(t - 1);
      if (
        (p & 1n) === 0n ||
        BigInt(t) * BigInt(t) + BigInt(v * v * L * L) * BigInt(absD) !== 4n * p
      )
        throw new PariBugError('modpoly_pickD_primes');
      if (p > 2n ** BigInt(FF_BITS)) break;
      if (p > BigInt(Number.MAX_SAFE_INTEGER))
        throw new NotImplementedError(
          'SAGE_NOT_IMPLEMENTED: modpoly_pickD_primes needs primes > 2^53'
        );
      const pnum = Number(p);
      if (xprimes.length) {
        while (m < xprimes.length && xprimes[m]! < p) m++;
        if (m < xprimes.length && p === xprimes[m]!) continue;
      }
      if (!modinv_good_prime(inv, pnum) || !isPrime(p)) continue;
      if (collect) {
        if (n >= max) break outer;
        primes.push(p);
        traces.push(t);
      }
      n++;
      bits += Math.log2(pnum);
      if (one_prime) break outer;
    }
  }
  if (!n) return { n: 0, totbits: 0, primes, traces };
  if (!totbits) totbits = Math.floor(bits);
  return { n, totbits, primes, traces };
}

const MAX_VOLCANO_FLOOR_SIZE = 100000000;

/** PARI `calc_primes_for_discriminants` (`polmodular.c:4064-4108`). */
function calc_primes_for_discriminants(
  Ds: DiscInfo[],
  DcntIn: number,
  L: number,
  minbits: number
): number {
  let Dcnt = DcntIn;
  let D1 = Ds[0]!.D1;
  for (let i = 1; i < Dcnt; i++) if (Ds[i]!.D1 > D1) D1 = Ds[i]!.D1;
  const n = Math.ceil(minbits / (Math.log2(L * L * -D1) - 2)) + 1;
  let primes: bigint[] = [];
  let totbits = 0;
  for (let i = 0; i < Dcnt && totbits < minbits; i++) {
    const r = modpoly_pickD_primes(true, n, primes, minbits - totbits, Ds[i]!);
    Ds[i]!.bits = r.totbits;
    Ds[i]!.nprimes = r.n;
    Ds[i]!.primes = r.primes;
    Ds[i]!.traces = r.traces;
    totbits += Ds[i]!.bits;
    if (totbits >= minbits || i === Dcnt - 1) {
      Dcnt = i + 1;
      break;
    }
    primes = [...primes, ...r.primes].sort((a, b) => (a < b ? -1 : a > b ? 1 : 0));
  }
  if (totbits < minbits) return 0;
  return Dcnt;
}

interface DEntry {
  /** bit 0: `ord(L0) < h(D)`; bits `2j`, `2j+1`: `1 + (D | prime(j))`.
   * `bigint` because upstream packs 31 two-bit fields into a 64-bit long. */
  m: bigint;
  D: number;
  h: number;
}

/** PARI `modpoly_pickD` (`polmodular.c:4110-4386`). */
function modpoly_pickD(
  Ds: DiscInfo[],
  L: number,
  inv: number,
  L0: number,
  max_L1: number,
  minbits: number,
  flags: number,
  tab: DEntry[]
): number {
  const pp = { p1: 1, p2: 1 };
  const modinv_deg = modinv_degree(pp, inv);
  const modinv_p1 = pp.p1;
  const modinv_p2 = pp.p2;
  const pfilter = modinv_pfilter(inv);
  const modinv_N = modinv_level(inv);
  const L_bits = Math.log2(L);
  if (!(L & 1)) throw new PariBugError('modpoly_pickD');

  const d =
    flags & MODPOLY_IGNORE_SPARSE_FACTOR ? L + 2 : ceildivuu(L + 1, modinv_sparse_factor(inv)) + 1;
  const use_L1 = kross(L0, L) > 0 || !!(flags & MODPOLY_USE_L1);

  let Dcnt = 0;
  let best_cost = 0;
  let totbits = 0;

  for (let D0_i = 0; D0_i < tab.length; D0_i++) {
    const D0_entry = tab[D0_i]!;
    const D0 = D0_entry.D;
    if (!modinv_good_disc(inv, D0)) continue;
    if (
      kross(D0, L) < 1 ||
      (modinv_p1 > 1 && kross(D0, modinv_p1) < 1) ||
      (modinv_p2 > 1 && kross(D0, modinv_p2) < 1)
    )
      continue;
    const deg = D0_entry.h;
    const h0 = D0_entry.m & 2n ? 2 * deg : deg;
    const n0 = h0 / (Number(D0_entry.m & 1n) + 1);

    let L1 = 0;
    let i: number;
    for (i = 1; i <= SMOOTH_PRIMES; i++) {
      const p = PARI_PRIMES[i]!;
      if (p <= L0) continue;
      if (Number((D0_entry.m >> BigInt(2 * i)) & 3n) === 1) {
        if (p <= max_L1 && modinv_N % p !== 0 && kross(p, L) < 0) {
          L1 = p;
          break;
        }
      }
    }
    if (i > SMOOTH_PRIMES && (n0 < h0 || use_L1)) continue;

    if (totbits > minbits && best_cost && h0 * (L - 1) > 3 * best_cost) break;

    const D0_bits = Math.log2(-D0);
    if (D0_bits + 2 * L_bits > 63) continue;

    const m0 = primeform_exp_order(L0, n0, L * L * D0, n0 * (L - 1));
    if (m0 < (L - 1) / 2) continue;
    const H_cost = 2 * deg * deg;

    const twofactor = (D0 & 7) === 5 ? (D0_entry.m & 0xcn ? 1 : 3) : 0;

    for (i = 0; i <= SMOOTH_PRIMES; i++) {
      let h1: number;
      let D1: number;
      let p: number;
      let q: number;
      let j: number;
      let p_bits: number;
      if (i) {
        if (modinv_odd_conductor(inv) && i === 1) continue;
        p = PARI_PRIMES[i]!;
        if (p > max_L1) break;
        if (p === L0 || p === L1 || p === L || p === modinv_p1 || p === modinv_p2) continue;
        p_bits = Math.log2(p);
        h1 = h0 * (p - Number((D0_entry.m >> BigInt(2 * i)) & 0x3n) + 1);
        for (j = 1, q = p; h1 < d; j++, q *= p, h1 *= p);
        D1 = q * q * D0;
        if (pfilter & IQ_FILTER_1MOD4 && (D1 & 0xf) === 0) continue;
      } else {
        h1 = h0;
        D1 = D0;
        p = q = j = 1;
        p_bits = 0;
      }
      if (twofactor && q & 1) {
        if (modinv_odd_conductor(inv)) continue;
        D1 *= 4;
        h1 *= twofactor;
      }
      if (totbits > minbits && best_cost && h1 * (L - 1) > 2.2 * best_cost) continue;
      if (D0_bits + 2 * j * p_bits + 2 * L_bits + (twofactor && q & 1 ? 2.0 : 0.0) > 63) continue;

      const g1 = { n1: 0, m: 0 };
      if (!check_generators(g1, D1, h1, n0, d, L0, L1)) continue;
      const n1 = g1.n1;

      let dl1: number;
      if (n1 >= h1) dl1 = -1;
      else {
        dl1 = primeform_discrete_log(L0, L, n1, D1);
        if (dl1 < 0) continue;
      }
      if (modinv_deg && orientation_ambiguity(D1, L0, modinv_p1, modinv_p2)) continue;

      const D2 = L * L * D1;
      const h2 = h1 * (L - 1);
      const g2 = { n1: 0, m: 0 };
      if (!check_generators(g2, D2, h2, n1, d * (L - 1), L0, L1)) continue;
      const n2 = g2.n1;
      const m = g2.m;

      if (m < (L - 1) / 2 || (!L1 && m < L - 1)) continue;
      let dl20 = n1;
      let dl21 = 0;
      if (m < L - 1) {
        const Q1a = qform_primeform2(L, D1);
        if (!Q1a) throw new PariBugError('modpoly_pickD');
        let Q2 = primeform(BigInt(D2), BigInt(L1));
        Q2 = qfbcomp(Q1a, Q2);
        let Q1 = primeform(BigInt(D2), BigInt(L0));
        const kk = (n2 & 1 ? 2 * n2 : n2) / (L - 1);
        Q1 = qfbpow(Q1, BigInt(kk));
        const X = qfi_Shanks(Q2, Q1, BigInt(L - 1));
        if (X === null) continue;
        dl20 = Number(X) * kk;
        dl21 = 1;
      }
      if (!(m < L - 1 || n2 < d * (L - 1)) && n1 >= d && !use_L1) L1 = 0;
      if (!L1 && use_L1) continue;
      if (L1 && !dl21) continue;

      const enum_cost = n2 * (5 * L0 * L0 + 0.25 * L1 * L1);
      const cost = enum_cost + H_cost;
      if (best_cost && cost > 2.2 * best_cost) break;
      if (best_cost && cost >= 0.99 * best_cost) continue;

      const Dinfo: DiscInfo = {
        inv,
        L,
        D0,
        D1,
        L0,
        L1,
        n1,
        n2,
        dl1,
        dl2_0: dl20,
        dl2_1: dl21,
        nprimes: 0,
        cost,
        bits: 0,
        primes: [],
        traces: [],
      };
      const r = modpoly_pickD_primes(false, 0, [], minbits, Dinfo);
      if (!r.n) continue;
      Dinfo.bits = r.totbits;

      let jj: number;
      for (jj = 0; jj < Dcnt; jj++) if (Dinfo.cost < Ds[jj]!.cost) break;
      if (n2 > MAX_VOLCANO_FLOOR_SIZE && n2 * (L1 ? 2 : 1) > 1.2 * (d * (L - 1))) continue;
      if (jj === Dcnt && Dcnt === MODPOLY_MAX_DCNT) continue;
      totbits += Dinfo.bits;
      if (Dcnt === MODPOLY_MAX_DCNT) totbits -= Ds[Dcnt - 1]!.bits;
      if (Dcnt < MODPOLY_MAX_DCNT) Dcnt++;
      let kk: number;
      for (kk = Dcnt - 1; kk > jj; kk--) Ds[kk] = Ds[kk - 1]!;
      Ds[kk] = Dinfo;
      best_cost = totbits > minbits ? Ds[Dcnt - 1]!.cost : 0;
      if (!i) break;
    }
  }
  if (!Dcnt) return 0;

  Dcnt = calc_primes_for_discriminants(Ds, Dcnt, L, minbits);

  for (let i = 0; i < Dcnt; i++)
    if (Ds[i]!.dl1 < 0) {
      Ds[i]!.dl1 = primeform_discrete_log(L0, L, Ds[i]!.n1, Ds[i]!.D1);
      if (Ds[i]!.dl1 < 0) throw new PariBugError('modpoly_pickD');
    }
  return Dcnt;
}

/**
 * PARI `scanD0` (`polmodular.c:4405-4487`): fundamental discriminants
 * `minD <= |D| <= maxD` whose class group is cyclic and generated by a form of
 * norm `L0`, with class number at most `maxh`.
 */
function scanD0(state: { minD: number }, maxD: number, maxh: number, L0: number): DEntry[] {
  const tab: DEntry[] = [];
  const minD = state.minD;
  for (let dd = minD; dd <= maxD; dd++) {
    /* squarefree, odd, d = 3 mod 4 */
    if ((dd & 3) !== 3) continue;
    if (!(dd & 1)) continue;
    const fa = factoru(BigInt(dd));
    let squarefree = true;
    for (const [, e] of fa) if (e > 1n) squarefree = false;
    if (!squarefree) continue;
    if (fa.length > 2) continue; /* restrict to possibly cyclic class groups */
    const D = -dd;
    if (kross(D, L0) < 1) continue;

    let L1 = fa.length > 1 && Number(fa[0]![0]) <= MAX_L1 ? Number(fa[0]![0]) : 0;
    const h = hclassno6_count(D) / 6;
    if (h > 2 * maxh || (!L1 && h > maxh)) continue;

    const DD = BigInt(D);
    const f = primeform(DD, BigInt(L0));
    const n = Number(qfi_order(qfbred(f), BigInt(h)));
    if (n < h / 2 || (!L1 && n < h)) continue;

    const k = fa.length;
    for (let j = 1; ; j++) {
      if (n === h || (L1 && qfi_Shanks(primeform(DD, BigInt(L1)), f, BigInt(n)) === null)) break;
      if (!L1) break;
      L1 = j < k && k > 1 && Number(fa[j - 1]![0]) <= MAX_L1 ? Number(fa[j - 1]![0]) : 0;
    }
    let m = n < h ? 1n : 0n;
    for (let j = 1; j <= SMOOTH_PRIMES; j++) {
      const x = BigInt(1 + kross(D, PARI_PRIMES[j]!));
      m |= x << BigInt(2 * j);
    }
    tab.push({ D, h, m });
  }
  /* `_qsort_cmp` (`polmodular.c:4388-4403`): by class number, then by |D| */
  tab.sort((x, y) => {
    const u = x.h * ((x.m & 2n ? 1 : 0) + 1);
    const v = y.h * ((y.m & 2n ? 1 : 0) + 1);
    if (u !== v) return u - v;
    return y.D - x.D;
  });
  state.minD = maxD + 3 - (maxD & 3);
  return tab;
}

/** PARI `discriminant_with_classno_at_least` (`polmodular.c:4489-4595`). */
function discriminant_with_classno_at_least(
  bestD: DiscInfo[],
  L: number,
  inv: number,
  Q: bigint | null,
  ignore_sparse: boolean
): number {
  const SMALL_L_BOUND = 101;
  let max_max_D = 160000 * (inv ? 2 : 1);
  const s = modinv_sparse_factor(inv);
  let maxD = 10000;
  const maxh = L / s < SMALL_L_BOUND ? 10 * L : 1.5 * L;
  const flags = ignore_sparse ? MODPOLY_IGNORE_SPARSE_FACTOR : 0;
  const L0 = select_L0(L, inv, 0);
  const max_L1 = Math.floor(L / 2) + 2;
  let minbits = modpoly_height_bound(L, inv);
  if (Q) minbits += Q.toString(2).length;
  const state = { minD: 7 };
  let best_cost = -1.0;
  let best_cnt = 0;

  while (!best_cnt) {
    while (maxD <= max_max_D) {
      const tab = scanD0(state, maxD, maxh, L0);
      const Ds: DiscInfo[] = [];
      const Dcnt = modpoly_pickD(Ds, L, inv, L0, max_L1, minbits, flags, tab);
      let cost = 0.0;
      if (Dcnt) {
        let n1 = 0;
        for (let i = 0; i < Dcnt; i++) {
          n1 = Math.max(n1, Ds[i]!.n1);
          cost += Ds[i]!.cost;
        }
        const eps = (n1 * s - L) / L;
        if (best_cost < 0.0 || cost < best_cost) {
          for (let i = 0; i < Dcnt; i++) bestD[i] = Ds[i]!;
          best_cost = cost;
          best_cnt = Dcnt;
          if (L / s <= SMALL_L_BOUND || eps < 0.05) break;
        }
      } else if (Math.log2(maxD) > 64 - 2 * (Math.log2(L) + 2)) {
        throw new PariArchError(`modular polynomial of level ${L} and invariant ${inv}`);
      }
      maxD *= 2;
      state.minD += 4;
    }
    max_max_D *= 2;
  }
  return best_cnt;
}

/* --- the per-prime computation (`polmodular.c:1364-1760`) --- */

/** PARI `get_Lsqr_cycle` (`polmodular.c:1366-1381`). */
function get_Lsqr_cycle(dinfo: DiscInfo): number[] {
  const n1 = dinfo.n1;
  const L = dinfo.L;
  const cyc = new Array<number>(L).fill(0);
  cyc[1] = 0;
  let i: number;
  for (i = 2; i <= L / 2; ++i) cyc[i] = cyc[i - 1]! + n1;
  if (!dinfo.L1) {
    for (; i < L; ++i) cyc[i] = cyc[i - 1]! + n1;
  } else {
    cyc[L - 1] = 2 * dinfo.n2 - n1 / 2;
    for (i = L - 2; i > L / 2; --i) cyc[i] = cyc[i + 1]! - n1;
  }
  return cyc;
}

/** PARI `update_Lsqr_cycle` (`polmodular.c:1383-1392`). */
function update_Lsqr_cycle(cyc: number[], dinfo: DiscInfo): void {
  const L = dinfo.L;
  for (let i = 1; i < L; ++i) cyc[i]!++;
  if (dinfo.L1 && cyc[L - 1] === 2 * dinfo.n2) {
    const n1 = dinfo.n1;
    for (let i = Math.floor(L / 2) + 1; i < L; ++i) cyc[i]! -= n1;
  }
}

/** PARI `oneroot_of_classpoly` (`polmodular.c:1394-1418`). */
function modpoly_oneroot_of_classpoly(
  hilb: bigint[],
  factu: Array<[bigint, bigint]>,
  ne: NormEqn,
  jdb: Array<ZM | null>
): bigint {
  const p = ne.p;
  const hilbp = hilb.map((c) => Fp_red(c, p));
  let j0 = Flx_oneroot(FpX_renormalize(hilbp), p);
  if (j0 === null)
    throw new PariBugError("oneroot_of_classpoly: Didn't find a root of the class polynomial");
  for (const [Lb, valb] of factu) {
    const L = Number(Lb);
    const phi = polmodular_db_getp(jdb, L, p);
    const val = Number(valb) + z_lval(ne.v, L);
    j0 = descend_volcano(phi, j0, p, 0, L, val, val);
  }
  return j0;
}

/** PARI `make_pcp_surface` (`polmodular.c:1420-1427`). */
function make_pcp_surface(dinfo: DiscInfo): Pcp {
  return {
    L: [dinfo.L0],
    n: [dinfo.n1],
    o: [dinfo.n1],
    m: [dinfo.n1],
    r: [],
    L0: 0,
    k: 1,
    enum_cnt: dinfo.n1,
    h: dinfo.n1,
    inv: dinfo.inv,
    D: 0,
    D0: 0,
    u: 1,
    fau: [],
  };
}

/** PARI `make_pcp_floor` (`polmodular.c:1429-1446`). */
function make_pcp_floor(dinfo: DiscInfo): Pcp {
  const k = dinfo.L1 ? 2 : 1;
  const L = k === 1 ? [dinfo.L0] : [dinfo.L0, dinfo.L1];
  const n = k === 1 ? [dinfo.n2] : [dinfo.n2, 2];
  const o = k === 1 ? [dinfo.n2] : [dinfo.n2, 2];
  const m: number[] = [];
  evec_n_to_m(m, n, k);
  return {
    L,
    n,
    o,
    m,
    r: new Array<number>((k * (k - 1)) / 2).fill(0),
    L0: 0,
    k,
    enum_cnt: dinfo.n2 * k,
    h: dinfo.n2 * k,
    inv: dinfo.inv,
    D: 0,
    D0: 0,
    u: 1,
    fau: [],
  };
}

/** PARI `carray_reverse_inplace` (`polmodular.c:1471-1478`). */
function carray_reverse_inplace(arr: bigint[], from: number, n: number): void {
  const lim = n >> 1;
  for (let i = 0; i < lim; i++) {
    const t = arr[from + i]!;
    arr[from + i] = arr[from + n - 1 - i]!;
    arr[from + n - 1 - i] = t;
  }
}

/** PARI `append_neighbours` (`polmodular.c:1480-1487`). */
function append_neighbours(
  rts: bigint[],
  surface_js: bigint[],
  njs: number,
  L: number,
  m: number,
  i: number
): void {
  const r_idx = (i - 1 + m) % njs;
  const l_idx = (((i - 1 - m) % njs) + njs) % njs;
  rts[L - 1] = surface_js[l_idx]!;
  rts[L] = surface_js[r_idx]!;
}

/** PARI `roots_to_coeffs` (`polmodular.c:1488-1505`). */
function roots_to_coeffs(rts: bigint[][], p: bigint, L: number): Flm {
  const M: Flm = [];
  for (let k = 0; k < L + 2; k++) M.push(new Array<bigint>(rts.length).fill(0n));
  for (let i = 0; i < rts.length; ++i) {
    const modpol = Flv_roots_to_pol(rts[i]!, p);
    for (let k = 0; k < L + 2; ++k) M[k]![i] = modpol[k] ?? 0n;
  }
  return M;
}

/** PARI `root_matrix` (`polmodular.c:1514-1572`), for non-double-eta invariants. */
function root_matrix(
  L: number,
  dinfo: DiscInfo,
  njinvs: number,
  surface_js: bigint[],
  floor_js: bigint[],
  n: bigint,
  card: bigint,
  val: number,
  ne: NormEqn
): bigint[][] {
  const m = dinfo.dl1;
  const njs = surface_js.length;
  const inv = dinfo.inv;
  const p = ne.p;
  const rt_mat: bigint[][] = [];
  for (let i = 0; i < njinvs; i++) rt_mat.push(new Array<bigint>(L + 1).fill(0n));

  const cyc = get_Lsqr_cycle(dinfo);
  let rts = rt_mat[0]!;
  for (let i = 1; i < L; i++) rts[i - 1] = floor_js[cyc[i]!]!;
  append_neighbours(rts, surface_js, njs, L, m, 1);

  update_Lsqr_cycle(cyc, dinfo);
  rts = rt_mat[1]!;
  for (let i = 1; i < L; i++) rts[i - 1] = floor_js[cyc[i]!]!;

  if (modinv_is_double_eta(inv))
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: root_matrix for double-eta class invariants'
    );
  const j1pr = modfn_preimage(rts[0]!, p, inv);
  const j1 = compute_L_isogenous_curve(L, n, ne, j1pr, card, val, false);
  const rev = j1 !== modfn_preimage(surface_js[1]!, p, inv);
  if (rev) carray_reverse_inplace(surface_js, 1, njs - 1);
  append_neighbours(rts, surface_js, njs, L, m, 2);

  for (let i = 3; i <= njinvs; ++i) {
    update_Lsqr_cycle(cyc, dinfo);
    rts = rt_mat[i - 1]!;
    for (let k = 1; k < L; k++) rts[k - 1] = floor_js[cyc[k]!]!;
    append_neighbours(rts, surface_js, njs, L, m, i);
  }
  return rt_mat;
}

/** PARI `interpolate_coeffs` (`polmodular.c:1574-1586`). */
function interpolate_coeffs(phi_modp: Flm, p: bigint, j_invs: bigint[], coeff_mat: Flm): void {
  const pols = Flv_Flm_polint(j_invs, coeff_mat, p);
  for (let i = 0; i < pols.length; ++i) {
    const pol = pols[i]!;
    for (let k = 0; k < pol.length; ++k) phi_modp[i]![k] = pol[k]!;
  }
}

/** PARI `Flv_lastnonzero` (`polmodular.c:1588-1595`), 1-based. */
function Flv_lastnonzero(v: bigint[]): number {
  for (let i = v.length; i > 0; --i) if (v[i - 1]) return i;
  return 0;
}

/** PARI `inflate_polys` (`polmodular.c:1597-1620`). */
function inflate_polys(phi: Flm, L: number, s: number): void {
  const deg = L + 1;
  const maxr = phi[0]!.length;
  for (let k = 0; k <= deg; ) {
    const c = (((L * (1 - k) + 1) % s) + s) % s;
    ++k;
    for (let i = Flv_lastnonzero(phi[k - 1]!); i > 0; --i) {
      const r = c + (i - 1) * s + 1;
      if (r > maxr) {
        phi[k - 1]![i - 1] = 0n;
        continue;
      }
      if (r !== i) {
        phi[k - 1]![r - 1] = phi[k - 1]![i - 1]!;
        phi[k - 1]![i - 1] = 0n;
      }
    }
  }
}

/** PARI `normalise_coeffs` (`polmodular.c:1629-1662`). */
function normalise_coeffs(coeffs: Flm, js: bigint[], L: number, s: number, p: bigint): void {
  if (s <= 1) throw new PariBugError('normalise_coeffs');
  const pows: bigint[][] = [];
  pows.push(new Array<bigint>(js.length).fill(1n % p));
  const inv_js = js.map((x) => Fp_inv(x, p));
  pows.push(inv_js);
  for (let k = 3; k <= s; ++k) pows.push(inv_js.map((x) => Fp_pow(x, BigInt(k - 1), p)));
  for (let k = 0; k < coeffs.length; ++k) {
    const c = (((L * (1 - k) + 1) % s) + s) % s;
    const col = coeffs[k]!;
    const C = pows[c]!;
    for (let i = 0; i < col.length; ++i) col[i] = Fp_mul(col[i]!, C[i]!, p);
  }
}

/* --- Velu and L-isogenous curves (`polmodular.c:1192-1362`) --- */

/** PARI `Fle_quotient_from_kernel_generator` (`polmodular.c:1192-1224`). */
function Fle_quotient_from_kernel_generator(
  a4: bigint,
  a6: bigint,
  pt: { x: bigint; y: bigint },
  p: bigint
): { a4: bigint; a6: bigint } {
  let t = 0n;
  let w = 0n;
  let Q: Fle = { x: pt.x, y: pt.y };
  let xQ: bigint;
  do {
    xQ = Q!.x;
    const yQ = Q!.y;
    const tQ = (6n * Fp_mul(xQ, xQ, p) + 2n * a4) % p;
    const uQ = (4n * Fp_mul(yQ, yQ, p) + Fp_mul(tQ, xQ, p)) % p;
    t = (t + tQ) % p;
    w = (w + uQ) % p;
    Q = Fle_add(pt, Q, a4, p);
  } while (Q !== null && Q.x !== xQ);
  return { a4: (a4 - 5n * t + 5n * p) % p, a6: (a6 - 7n * w + 7n * p) % p };
}

/** PARI `find_L_tors_point` (`polmodular.c:1231-1254`). */
function find_L_tors_point(
  a4: bigint,
  a6: bigint,
  p: bigint,
  n: bigint,
  L: bigint,
  val: number
): { P: Fle; i: number } {
  let P: Fle = null;
  for (;;) {
    const Q = random_Fle(a4, a6, p);
    P = Fle_mulu(Q, n, a4, p);
    if (P !== null) break;
  }
  let i = 0;
  for (; i < val; ++i) {
    const Q = Fle_mulu(P, L, a4, p);
    if (Q === null) break;
    P = Q;
  }
  return { P, i };
}

/** PARI `select_curve_with_L_tors_point` (`polmodular.c:1256-1287`). */
function select_curve_with_L_tors_point(
  L: number,
  j: bigint,
  n: bigint,
  card: bigint,
  val: number,
  ne: NormEqn
): { a4: bigint; a6: bigint; P: Fle } {
  const p = ne.p;
  if (card % BigInt(L) !== 0n)
    throw new PariBugError('select_curve_with_L_tors_point: Cardinality not divisible by L');
  let { a4: A4, a6: A6 } = Fl_ellj_to_a4a6(j, p);
  let { a4: A4t, a6: A6t } = Fl_elltwist_disc(A4, A6, ne.T, p);
  for (;;) {
    const { P, i } = find_L_tors_point(A4, A6, p, n, BigInt(L), val);
    if (i < val) return { a4: A4, a6: A6, P };
    [A4, A4t] = [A4t, A4];
    [A6, A6t] = [A6t, A6];
  }
}

/** PARI `verify_L_sylow_is_cyclic` (`polmodular.c:1289-1310`). */
function verify_L_sylow_is_cyclic(e: bigint, a4: bigint, a6: bigint, p: bigint): boolean {
  const N_RETRIES = 3;
  for (let i = 0; i < N_RETRIES; ++i) {
    let P = random_Fle(a4, a6, p);
    P = Fle_mulu(P, e, a4, p);
    if (P !== null) return true;
  }
  return false;
}

/** PARI `find_noniso_L_isogenous_curve` (`polmodular.c:1312-1345`). */
function find_noniso_L_isogenous_curve(
  L: number,
  n: bigint,
  ne: NormEqn,
  e: bigint,
  val: number,
  a4: bigint,
  a6: bigint,
  init_pt: Fle,
  verify: boolean
): bigint {
  const p = ne.p;
  let pt = init_pt;
  for (;;) {
    const img = Fle_quotient_from_kernel_generator(a4, a6, pt!, p);
    if (!verify || verify_L_sylow_is_cyclic(e, img.a4, img.a6, p))
      return Fl_ellj(img.a4, img.a6, p);
    pt = find_L_tors_point(a4, a6, p, n, BigInt(L), val).P;
  }
}

/** PARI `compute_L_isogenous_curve` (`polmodular.c:1347-1364`). */
function compute_L_isogenous_curve(
  L: number,
  n: bigint,
  ne: NormEqn,
  j: bigint,
  card: bigint,
  val: number,
  verify: boolean
): bigint {
  if (ne.p < 5n || j === 0n || j === 1728n % ne.p)
    throw new PariBugError('compute_L_isogenous_curve');
  const { a4, a6, P } = select_curve_with_L_tors_point(L, j, n, card, val, ne);
  const e = card / BigInt(L);
  if (e * BigInt(L) !== card) throw new PariBugError('compute_L_isogenous_curve');
  return find_noniso_L_isogenous_curve(L, n, ne, e, val, a4, a6, P, verify);
}

/** PARI `polmodular_split_p_Flm` (`polmodular.c:1694-1760`), Sutherland Alg. 2.1. */
function polmodular_split_p_Flm(
  L: number,
  hilb: bigint[],
  factu: Array<[bigint, bigint]>,
  ne: NormEqn,
  db: PolmodularDB,
  G_surface: Pcp,
  G_floor: Pcp,
  dinfo: DiscInfo
): Flm {
  const p = ne.p;
  const inv = dinfo.inv;
  const s = modinv_sparse_factor(inv);
  const nj_selected = Math.ceil((L + 1) / s) + 1;
  const jdb = polmodular_db_for_inv(db, INV_J);
  const fdb = polmodular_db_for_inv(db, inv);

  const card = p + 1n - BigInt(ne.t);
  const { v: val, n } = u_lvalrem(card, BigInt(L));

  const j0 = modpoly_oneroot_of_classpoly(hilb, factu, ne, jdb);
  const j0pr = compute_L_isogenous_curve(L, n, ne, j0, card, val, true);
  if (modinv_is_double_eta(inv))
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: polmodular_split_p_Flm for double-eta class invariants'
    );
  const j0_rt = modfn_root(j0, ne, inv);
  const j0pr_rt = modfn_root(j0pr, ne, inv);

  const surface_js = enum_roots(j0_rt, ne, fdb, G_surface, null);
  if (!surface_js) throw new PariBugError('polmodular_split_p_Flm [surface]');
  /* L^2 D is the discriminant for the order R = Z + L OO */
  const eqn: NormEqn = {
    ...ne,
    D: L * L * ne.D,
    u: L * ne.u,
    v: L * ne.u * ne.v,
  };
  const floor_js = enum_roots(j0pr_rt, eqn, fdb, G_floor, null);
  if (!floor_js) throw new PariBugError('polmodular_split_p_Flm [floor]');

  const rts = root_matrix(L, dinfo, nj_selected, surface_js, floor_js, n, card, val, ne);

  let switched_signs = false;
  for (;;) {
    const coeffs = roots_to_coeffs(rts, p, L);
    const surf = surface_js.slice(0, nj_selected);
    if (s > 1) {
      normalise_coeffs(coeffs, surf, L, s, p);
      for (let i = 0; i < surf.length; i++) surf[i] = Fp_pow(surf[i]!, BigInt(s), p);
    }
    const phi_modp: Flm = [];
    for (let i = 0; i < L + 2; i++) phi_modp.push(new Array<bigint>(L + 2).fill(0n));
    interpolate_coeffs(phi_modp, p, surf, coeffs);
    if (s > 1) inflate_polys(phi_modp, L, s);

    if (phi_modp[L]![L] === p - 1n) return phi_modp;
    if (switched_signs) throw new PariBugError('polmodular_split_p_Flm');
    for (let i = 0; i < rts.length; ++i) {
      surface_js[i] = (p - surface_js[i]!) % p;
      rts[i]![L - 1] = (p - rts[i]![L - 1]!) % p;
      rts[i]![L] = (p - rts[i]![L]!) % p;
    }
    switched_signs = true;
  }
}

/** PARI `Flv_deriv_pre_inplace` (`polmodular.c:1762-1768`). */
function Flv_deriv_inplace(v: bigint[], deg: number, p: bigint): void {
  let d = BigInt(deg) % p;
  for (let i = v.length - 1; i > 0; --i, --d) v[i] = Fp_mul(v[i - 1]!, (d + p) % p, p);
  v[0] = 0n;
}

/** PARI `eval_modpoly_modp` (`polmodular.c:1770-1787`). */
function eval_modpoly_modp(Tp: Flm, j_powers: bigint[], p: bigint, compute_derivs: boolean): Flm {
  const L = j_powers.length - 2;
  const jp = j_powers.map((x) => Fp_red(x, p));
  const out: Flm = [Flm_Flc_mul(Tp, jp, p)];
  if (compute_derivs) {
    Flv_deriv_inplace(jp, L + 1, p);
    out.push(Flm_Flc_mul(Tp, jp, p));
    Flv_deriv_inplace(jp, L + 1, p);
    out.push(Flm_Flc_mul(Tp, jp, p));
  }
  return out;
}

/** PARI `polmodular0_powerup_ZM` (`polmodular.c:2327-2395`). */
function polmodular0_powerup_ZM(L: number, inv: number, db: PolmodularDB): ZM {
  const parent = modinv_parent(inv);
  const e = BigInt(modinv_parent_power(inv));
  const Ds: DiscInfo[] = [];
  const nDs = discriminant_with_classno_at_least(Ds, L, inv, null, true);
  if (nDs !== 1) throw new PariBugError('polmodular0_powerup_ZM');
  const D = Ds[0]!.D1;
  const nprimes = Ds[0]!.nprimes + 1;
  const mp = polmodular0_ZM(L, parent, null, null, false, db);
  const H = polclass0(D, parent, db);

  const N = L + 2;
  if (H.length - 1 < N) throw new PariBugError('polmodular0_powerup_ZM');

  const mats: Flm[] = [];
  const plist: bigint[] = [];
  for (let sidx = 1; sidx < nprimes; ++sidx) {
    const p = Ds[0]!.primes[sidx - 1]!;
    const phi_modp: Flm = [];
    for (let i = 0; i < L + 2; i++) phi_modp.push(new Array<bigint>(N).fill(0n));
    const Hp = FpX_renormalize(H.map((c) => Fp_red(c, p)));
    const Hrts = Flx_roots(Hp, p);
    if (Hrts.length < N) throw new PariBugError('polmodular0_powerup_ZM');
    const js = new Array<bigint>(N);
    for (let i = 0; i < N; ++i) js[i] = Fp_pow(Hrts[i]!, e, p);

    const Phip = ZM_to_Flm(mp, p);
    const coeff_mat: Flm = [];
    for (let i = 0; i < L + 2; i++) coeff_mat.push(new Array<bigint>(N).fill(0n));
    for (let i = 0; i < N; ++i) {
      let phi_at_ji = Flm_Fl_polmodular_evalx(Phip, L, Hrts[i]!, p);
      const mprts = Flx_roots(phi_at_ji, p);
      if (mprts.length !== L + 1) throw new PariBugError('polmodular0_powerup_ZM');
      const powered = mprts.map((x) => Fp_pow(x, e, p));
      phi_at_ji = Flv_roots_to_pol(powered, p);
      for (let k = 0; k < L + 2; ++k) coeff_mat[k]![i] = phi_at_ji[k] ?? 0n;
    }
    interpolate_coeffs(phi_modp, p, js, coeff_mat);
    mats.push(phi_modp);
    plist.push(p);
  }
  return nmV_chinese_center(mats, plist);
}

/** The CM branch of PARI `polmodular0_ZM` (`polmodular.c:1878-1934`). */
function polmodular_CM_ZM(
  L: number,
  inv: number,
  J: bigint | null,
  Q: bigint | null,
  compute_derivs: boolean,
  db: PolmodularDB
): ZM {
  if (inv !== INV_J && inv !== INV_G2)
    throw new NotImplementedError(
      `SAGE_NOT_IMPLEMENTED: polmodular for invariant ${inv} at level ${L}; only ` +
        'INV_J and INV_G2 are ported (Weber/double-eta class invariants need ' +
        'double_eta_raw / modinv_f_from_j from polmodular.c:500-870)'
    );
  const Ds: DiscInfo[] = [];
  const Dcnt = discriminant_with_classno_at_least(Ds, L, inv, Q, false);
  const mats: Flm[] = [];
  const plist: bigint[] = [];
  let j_powers: bigint[] | null = null;
  if (J !== null && Q !== null) j_powers = Fp_powers(J, L + 1, Q);

  for (let d = 0; d < Dcnt; d++) {
    const dinfo = Ds[d]!;
    const D = dinfo.D1;
    const DK = dinfo.D0;
    const cond = Math.round(Math.sqrt(D / DK));
    const factu = factoru(BigInt(cond));

    polmodular_db_add_level(db, dinfo.L0, inv);
    if (dinfo.L1) polmodular_db_add_level(db, dinfo.L1, inv);
    const hilb = polclass0(DK, INV_J, db);
    if (cond > 1)
      polmodular_db_add_levels(
        db,
        factu.map(([q]) => Number(q)),
        INV_J
      );
    const G_surface = make_pcp_surface(dinfo);
    const G_floor = make_pcp_floor(dinfo);

    for (let i = 0; i < dinfo.nprimes; i++) {
      const p = dinfo.primes[i]!;
      const t = dinfo.traces[i]!;
      const q = (4n * p - BigInt(t) * BigInt(t)) / BigInt(-D);
      const vL = isqrt(q);
      if (vL * vL !== q) throw new PariBugError('polmodular_worker');
      const ne = norm_eqn_set(D, t, cond, Number(vL), [], p);
      let Tp = polmodular_split_p_Flm(L, hilb, factu, ne, db, G_surface, G_floor, dinfo);
      if (j_powers) Tp = eval_modpoly_modp(Tp, j_powers, p, compute_derivs);
      mats.push(Tp);
      plist.push(p);
    }
  }
  let modpoly = nmV_chinese_center(mats, plist);
  if (Q !== null) modpoly = modpoly.map((col) => col.map((x) => Fp_red(x, Q)));
  return modpoly;
}

/* ------------------------------------------------------------------ */
/* Internals exposed for testing only                                  */
/* ------------------------------------------------------------------ */

/** @internal Not part of the public API; used by `polmodular.test.ts`. */
export const _internal = {
  Flx_roots,
  Flx_oneroot,
  Flx_nbroots,
  Flv_roots_to_pol,
  Flv_Flm_polint,
  Flx_div_by_X_x,
  Fle_add,
  Fle_mulu,
  random_Fle,
  Fl_ellj,
  Fl_ellj_to_a4a6,
  Fl_elltwist_disc,
  hclassno6_count,
  quadclassnos,
  quadnegclassnou,
  qfi_order,
  qfi_Shanks,
  classgp_make_pcp,
  select_classpoly_primes,
  select_L0,
  scanD0,
  discriminant_with_classno_at_least,
  modpoly_height_bound,
  ncV_chinese_center,
  descend_volcano,
  ascend_volcano,
  j_level_in_volcano,
  enum_roots,
  norm_eqn_set,
  find_j_inv_with_given_trace,
  SMOOTH_INTS,
  HURWITZ_RATIO,
};
