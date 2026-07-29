/**
 * Tests for the port of PARI's Galois machinery (`galconj.c`).
 *
 * Oracles, in decreasing order of authority:
 *
 *  1. **PARI itself**, run through the SageMath installed on this machine:
 *     `galoisinit(T)` for its relative orders, `nfgaloisconj(T,4)` for the
 *     automorphism polynomials (printed exactly as `gp` prints them),
 *     `#galoissubgroups(G)` and the degrees of `galoissubfields(G,1)`.
 *     Every string in {@link ORACLE} was produced by that PARI and is copied
 *     verbatim.
 *  2. **Exact algebraic identities**, independent of anything p-adic:
 *     - every permutation returned is a genuine automorphism, i.e.
 *       `T(sigma(x)) = 0 in Q[x]/(T)` computed in exact integer arithmetic;
 *     - the permutations are pairwise distinct and closed under composition;
 *     - `galoisfixedfield`: `deg P = [G:H]` (the Galois correspondence),
 *       `P(S) = 0 in Q[x]/(T)`, and (flag 2) the returned factorisation
 *       multiplies back to `T` over `Q[y]/(P)`.
 *  3. **Brute force** for the low-level layers (roots mod p, factorisation
 *     mod p, Hensel lifts, inverse Vandermonde, subgroup lattices).
 */

import { describe, expect, test } from 'bun:test';
import { FpX_add, FpX_degree, FpX_mul, FpX_red, FpX_rem, FpX_sub } from './ffinit.js';
import {
  type GaloisInit,
  type Group,
  type Perm,
  type QPoly,
  type ZX,
  FpV_invVandermonde,
  FpX_factor_squarefree,
  FpX_is_squarefree,
  FpX_nbfact_by_degree,
  FpX_roots,
  FpXQ_minpoly,
  NotImplementedError,
  PariDomainError,
  PariImplError,
  QPoly_to_fractions,
  ZX_add,
  ZX_degree,
  ZX_disc,
  ZX_divrem_monic,
  ZX_mul,
  ZX_renormalize,
  ZX_resultant,
  ZpX_liftfact,
  ZpX_roots,
  bezout_lift_fact,
  galoisconj4,
  galoisfixedfield,
  galoisinit,
  galoispermtopol,
  galoissubgroups,
  galoisvecpermtopol,
  group_elts,
  group_order,
  identity_perm,
  indexpartial,
  listznstarelts,
  perm_cycles,
  perm_mul,
  perm_powu,
  vecperm_orbits,
} from './galconj.js';

/* ------------------------------------------------------------------ */
/* helpers                                                             */
/* ------------------------------------------------------------------ */

/** print a QPoly the way `gp` prints a t_POL, so we can compare with PARI */
function gpstr(q: QPoly): string {
  const fr = QPoly_to_fractions(q);
  const parts: string[] = [];
  for (let i = fr.length - 1; i >= 0; i--) {
    const [n, d] = fr[i]!;
    if (n === 0n) continue;
    const neg = n < 0n;
    const an = neg ? -n : n;
    const coef = d === 1n ? (an === 1n && i > 0 ? '' : an.toString()) : `${an}/${d}`;
    const mon = i === 0 ? '' : i === 1 ? 'x' : `x^${i}`;
    const term = coef && mon ? `${coef}*${mon}` : coef + mon;
    parts.push((parts.length === 0 ? (neg ? '-' : '') : neg ? ' - ' : ' + ') + term);
  }
  return parts.length ? parts.join('') : '0';
}

/** `A(x)^k` reduced mod the monic `T` */
function powmodT(A: ZX, k: number, T: ZX): ZX {
  let r: ZX = [1n];
  for (let i = 0; i < k; i++) r = ZX_divrem_monic(ZX_mul(r, A), T)[1];
  return r;
}

/**
 * Exact test that `F(sigma) = 0` in `Q[x]/(T)`, where `sigma = A/d`:
 * `sum_k F_k A^k d^(m-k) = 0 mod T` in Z[x], `m = deg F`.
 */
function vanishesAt(F: ZX, sigma: QPoly, T: ZX): boolean {
  const m = ZX_degree(F);
  const A = sigma.num;
  const d = sigma.den;
  let acc: ZX = [];
  for (let k = 0; k <= m; k++) {
    if (F[k] === 0n || F[k] === undefined) continue;
    const t = ZX_mul(powmodT(A, k, T), [F[k]! * d ** BigInt(m - k)]);
    acc = ZX_add(acc, t);
  }
  return ZX_divrem_monic(acc, T)[1].length === 0;
}

/* rational-coefficient polynomials over Q[y]/(P), used to check the
 * factorisation returned by galoisfixedfield(..., 2) */
type QModP = { num: ZX; den: bigint };
function qmulP(a: QModP, b: QModP, P: ZX): QModP {
  return { num: ZX_divrem_monic(ZX_mul(a.num, b.num), P)[1], den: a.den * b.den };
}
function qaddP(a: QModP, b: QModP): QModP {
  return { num: ZX_add(ZX_mul(a.num, [b.den]), ZX_mul(b.num, [a.den])), den: a.den * b.den };
}
function qisZero(a: QModP): boolean {
  return ZX_renormalize(a.num).length === 0;
}
/** product of polynomials in x whose coefficients live in Q[y]/(P) */
function prodFactors(factors: QModP[][], P: ZX): QModP[] {
  let acc: QModP[] = [{ num: [1n], den: 1n }];
  for (const f of factors) {
    const out: QModP[] = new Array(acc.length + f.length - 1)
      .fill(null)
      .map(() => ({ num: [] as ZX, den: 1n }));
    for (let i = 0; i < acc.length; i++)
      for (let j = 0; j < f.length; j++)
        out[i + j] = qaddP(out[i + j]!, qmulP(acc[i]!, f[j]!, P));
    acc = out;
  }
  return acc;
}

function isPerm(p: Perm, n: number): boolean {
  if (p.length !== n + 1) return false;
  const seen = new Set<number>();
  for (let i = 1; i <= n; i++) {
    if (p[i]! < 1 || p[i]! > n || seen.has(p[i]!)) return false;
    seen.add(p[i]!);
  }
  return true;
}

function key(p: Perm): string {
  return p.slice(1).join(',');
}

/* ------------------------------------------------------------------ */
/* The PARI oracle                                                     */
/* ------------------------------------------------------------------ */

interface Oracle {
  /** T, little-endian */
  pol: bigint[];
  /** PARI `galoisinit(T)[8]`, the relative orders of the generators */
  ord: number[];
  /** PARI `#galoissubgroups(galoisinit(T))` */
  nsub: number;
  /** degrees of PARI `galoissubfields(galoisinit(T),1)`, sorted */
  degs: number[];
  /** PARI `nfgaloisconj(T,4)`, printed by gp */
  conj: string[];
}

const O = (pol: bigint[], ord: number[], nsub: number, degs: number[], conj: string[]): Oracle => ({
  pol,
  ord,
  nsub,
  degs,
  conj,
});

const ORACLE: Record<string, Oracle> = {
  /* Q(sqrt 2), C2 */
  'x^2-2': O([-2n, 0n, 1n], [2], 2, [1, 2], ['-x', 'x']),
  /* the cyclic cubic field of conductor 9, C3 */
  'x^3-3*x+1': O([1n, -3n, 0n, 1n], [3], 2, [1, 3], ['x', '-x^2 - x + 2', 'x^2 - 2']),
  /* Q(zeta_5), C4 */
  'polcyclo(5)': O([1n, 1n, 1n, 1n, 1n], [4], 3, [1, 2, 4], [
    'x',
    'x^2',
    '-x^3 - x^2 - x - 1',
    'x^3',
  ]),
  /* Q(zeta_7), C6 */
  'polcyclo(7)': O([1n, 1n, 1n, 1n, 1n, 1n, 1n], [6], 4, [1, 2, 3, 6], [
    'x',
    'x^2',
    'x^3',
    'x^4',
    '-x^5 - x^4 - x^3 - x^2 - x - 1',
    'x^5',
  ]),
  /* Q(zeta_8), C2 x C2 */
  'polcyclo(8)': O([1n, 0n, 0n, 0n, 1n], [2, 2], 5, [1, 2, 2, 2, 4], ['-x', 'x', '-x^3', 'x^3']),
  /* Q(zeta_9), C6 */
  'polcyclo(9)': O([1n, 0n, 0n, 1n, 0n, 0n, 1n], [6], 4, [1, 2, 3, 6], [
    'x',
    'x^2',
    '-x^4 - x',
    'x^4',
    '-x^5 - x^2',
    'x^5',
  ]),
  /* Q(zeta_12), C2 x C2 */
  'polcyclo(12)': O([1n, 0n, -1n, 0n, 1n], [2, 2], 5, [1, 2, 2, 2, 4], [
    '-x',
    'x',
    '-x^3 + x',
    'x^3 - x',
  ]),
  /* Q(zeta_15), C4 x C2 */
  'polcyclo(15)': O(
    [1n, -1n, 0n, 1n, -1n, 1n, 0n, -1n, 1n],
    [4, 2],
    8,
    [1, 2, 2, 2, 4, 4, 4, 8],
    [
      'x',
      'x^2',
      'x^4',
      '-x^6 - x',
      '-x^7 + x^5 - x^4 - x + 1',
      '-x^7 + x^6 - x^4 + x^3 - x^2 + 1',
      'x^7 - x^5 + x^4 - x^3 + x - 1',
      'x^7',
    ]
  ),
  /* Q(zeta_16), C4 x C2 */
  'polcyclo(16)': O(
    [1n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 1n],
    [4, 2],
    8,
    [1, 2, 2, 2, 4, 4, 4, 8],
    ['-x', 'x', '-x^3', 'x^3', '-x^5', 'x^5', '-x^7', 'x^7']
  ),
  /* the splitting field of x^3-2: S3 (PARI's own galoisinit example) */
  'x^6+108': O(
    [108n, 0n, 0n, 0n, 0n, 0n, 1n],
    [3, 2],
    6,
    [1, 2, 3, 3, 3, 6],
    [
      '-x',
      'x',
      '-1/12*x^4 - 1/2*x',
      '-1/12*x^4 + 1/2*x',
      '1/12*x^4 - 1/2*x',
      '1/12*x^4 + 1/2*x',
    ]
  ),
  /* the splitting field of x^4-2: D4 */
  'x^8+28*x^4+2500': O(
    [2500n, 0n, 0n, 0n, 28n, 0n, 0n, 0n, 1n],
    [2, 2, 2],
    10,
    [1, 2, 2, 2, 4, 4, 4, 4, 4, 8],
    [
      '-x',
      'x',
      '-1/48*x^5 - 7/24*x',
      '-1/60*x^5 + 11/30*x',
      '-1/80*x^5 - 39/40*x',
      '1/80*x^5 + 39/40*x',
      '1/60*x^5 - 11/30*x',
      '1/48*x^5 + 7/24*x',
    ]
  ),
  /* the quaternion field of discriminant 1296: Q8 */
  'x^8-12*x^6+36*x^4-36*x^2+9': O(
    [9n, 0n, -36n, 0n, 36n, 0n, -12n, 0n, 1n],
    [2, 2, 2],
    6,
    [1, 2, 2, 2, 4, 8],
    [
      '-x',
      'x',
      '-1/3*x^5 + 3*x^3 - 3*x',
      '1/3*x^5 - 3*x^3 + 3*x',
      '-2/3*x^7 + 22/3*x^5 - 17*x^3 + 10*x',
      '-1/3*x^7 + 10/3*x^5 - 5*x^3 - x',
      '1/3*x^7 - 10/3*x^5 + 5*x^3 + x',
      '2/3*x^7 - 22/3*x^5 + 17*x^3 - 10*x',
    ]
  ),
  /* the splitting field of x^5-5x+12: D5 */
  'x^10+10*x^8+125*x^6+500*x^4+2500*x^2+4000': O(
    [4000n, 0n, 2500n, 0n, 500n, 0n, 125n, 0n, 10n, 0n, 1n],
    [5, 2],
    8,
    [1, 2, 5, 5, 5, 5, 5, 10],
    [
      '-x',
      'x',
      '-1/960*x^9 - 1/600*x^8 - 1/120*x^7 - 1/40*x^6 - 101/960*x^5 - 1/4*x^4 - 29/96*x^3 - 5/4*x^2 - 23/12*x - 10/3',
      '-1/960*x^9 - 1/800*x^8 - 1/120*x^7 - 1/120*x^6 - 101/960*x^5 - 11/96*x^4 - 29/96*x^3 + 5/48*x^2 - 17/12*x - 5/6',
      '-1/960*x^9 + 1/800*x^8 - 1/120*x^7 + 1/120*x^6 - 101/960*x^5 + 11/96*x^4 - 29/96*x^3 - 5/48*x^2 - 17/12*x + 5/6',
      '-1/960*x^9 + 1/600*x^8 - 1/120*x^7 + 1/40*x^6 - 101/960*x^5 + 1/4*x^4 - 29/96*x^3 + 5/4*x^2 - 23/12*x + 10/3',
      '1/960*x^9 - 1/600*x^8 + 1/120*x^7 - 1/40*x^6 + 101/960*x^5 - 1/4*x^4 + 29/96*x^3 - 5/4*x^2 + 23/12*x - 10/3',
      '1/960*x^9 - 1/800*x^8 + 1/120*x^7 - 1/120*x^6 + 101/960*x^5 - 11/96*x^4 + 29/96*x^3 + 5/48*x^2 + 17/12*x - 5/6',
      '1/960*x^9 + 1/800*x^8 + 1/120*x^7 + 1/120*x^6 + 101/960*x^5 + 11/96*x^4 + 29/96*x^3 - 5/48*x^2 + 17/12*x + 5/6',
      '1/960*x^9 + 1/600*x^8 + 1/120*x^7 + 1/40*x^6 + 101/960*x^5 + 1/4*x^4 + 29/96*x^3 + 5/4*x^2 + 23/12*x + 10/3',
    ]
  ),
};


/**
 * PARI `galoissubfields(galoisinit(T), 1)` -- every subfield of every field
 * above, as gp prints them, sorted.  Copied verbatim from that PARI.
 */
const SUBFIELDS: Record<string, string[]> = {
  'x^2-2': ['x', 'x^2 - 2'],
  'polcyclo(5)': ['x + 1', 'x^2 + x - 1', 'x^4 + x^3 + x^2 + x + 1'],
  'polcyclo(7)': ['x + 1', 'x^2 + x + 2', 'x^3 + x^2 - 2*x - 1', 'x^6 + x^5 + x^4 + x^3 + x^2 + x + 1'],
  'polcyclo(8)': ['x', 'x^2 + 2', 'x^2 + 4', 'x^2 - 2', 'x^4 + 1'],
  'polcyclo(9)': ['x', 'x^2 + 3*x + 9', 'x^3 - 3*x + 1', 'x^6 + x^3 + 1'],
  'polcyclo(12)': ['x', 'x^2 + 1', 'x^2 - 2*x + 4', 'x^2 - 3', 'x^4 - x^2 + 1'],
  'polcyclo(15)': ['x - 1', 'x^2 - x + 1', 'x^2 - x + 4', 'x^2 - x - 1', 'x^4 - x^3 + 2*x^2 + x + 1', 'x^4 - x^3 + x^2 - x + 1', 'x^4 - x^3 - 4*x^2 + 4*x + 1', 'x^8 - x^7 + x^5 - x^4 + x^3 - x + 1'],
  'polcyclo(16)': ['x', 'x^2 + 16', 'x^2 + 8', 'x^2 - 8', 'x^4 + 16', 'x^4 + 4*x^2 + 2', 'x^4 - 4*x^2 + 2', 'x^8 + 1'],
  'x^6+108': ['x', 'x^2 + 972', 'x^3 + 54', 'x^3 + 864', 'x^3 - 54', 'x^6 + 108'],
  'x^8+28*x^4+2500': ['x', 'x^2 + 112*x + 40000', 'x^2 + 512', 'x^2 - 288', 'x^4 + 112*x^2 + 40000', 'x^4 + 648', 'x^4 + 8', 'x^4 - 32', 'x^4 - 512', 'x^8 + 28*x^4 + 2500'],
  'x^8-12*x^6+36*x^4-36*x^2+9': ['x', 'x^2 - 24*x + 120', 'x^2 - 24*x + 72', 'x^2 - 24*x + 96', 'x^4 - 24*x^3 + 144*x^2 - 288*x + 144', 'x^8 - 12*x^6 + 36*x^4 - 36*x^2 + 9'],
  'x^10+10*x^8+125*x^6+500*x^4+2500*x^2+4000': ['x', 'x^10 + 10*x^8 + 125*x^6 + 500*x^4 + 2500*x^2 + 4000', 'x^2 + 100000', 'x^5 + 100*x^2 + 275*x + 200', 'x^5 + 20*x^4 + 500*x^3 + 4000*x^2 + 40000*x + 128000', 'x^5 + 25*x^3 + 50*x^2 + 150*x + 900', 'x^5 + 25*x^3 - 50*x^2 + 150*x - 900', 'x^5 - 100*x^2 + 275*x - 200'],
};

/** print a ZX the way gp prints a t_POL */
function gpZX(p: ZX): string {
  const parts: string[] = [];
  for (let i = p.length - 1; i >= 0; i--) {
    const c = p[i]!;
    if (c === 0n) continue;
    const neg = c < 0n;
    const a = neg ? -c : c;
    const coef = a === 1n && i > 0 ? '' : a.toString();
    const mon = i === 0 ? '' : i === 1 ? 'x' : `x^${i}`;
    const term = coef && mon ? `${coef}*${mon}` : coef + mon;
    parts.push((parts.length === 0 ? (neg ? '-' : '') : neg ? ' - ' : ' + ') + term);
  }
  return parts.length ? parts.join('') : '0';
}

/* the A4 case is kept separate: its automorphisms are long */
const A4_POL: bigint[] = [
  331776n, 0n, 165888n, 0n, -16128n, 0n, 1664n, 0n, 96n, 0n, 0n, 0n, 1n,
];
const A4_CONJ = [
  '-x',
  'x',
  '-275/19650816*x^11 - 77/3275136*x^9 - 2947/1637568*x^7 - 15389/614088*x^5 + 22885/102348*x^3 - 6764/2843*x',
  '-275/39301632*x^11 - 245/26201088*x^10 - 77/6550272*x^9 + 125/1091712*x^8 - 2947/3275136*x^7 - 385/545856*x^6 - 15389/1228176*x^5 + 821/409392*x^4 + 22885/204696*x^3 + 3035/8529*x^2 - 9607/5686*x - 6122/2843',
  '-275/39301632*x^11 - 245/26201088*x^10 - 77/6550272*x^9 + 125/1091712*x^8 - 2947/3275136*x^7 - 385/545856*x^6 - 15389/1228176*x^5 + 821/409392*x^4 + 22885/204696*x^3 + 3035/8529*x^2 - 3921/5686*x - 6122/2843',
  '-275/39301632*x^11 + 245/26201088*x^10 - 77/6550272*x^9 - 125/1091712*x^8 - 2947/3275136*x^7 + 385/545856*x^6 - 15389/1228176*x^5 - 821/409392*x^4 + 22885/204696*x^3 - 3035/8529*x^2 - 9607/5686*x + 6122/2843',
  '-275/39301632*x^11 + 245/26201088*x^10 - 77/6550272*x^9 - 125/1091712*x^8 - 2947/3275136*x^7 + 385/545856*x^6 - 15389/1228176*x^5 - 821/409392*x^4 + 22885/204696*x^3 - 3035/8529*x^2 - 3921/5686*x + 6122/2843',
  '275/39301632*x^11 - 245/26201088*x^10 + 77/6550272*x^9 + 125/1091712*x^8 + 2947/3275136*x^7 - 385/545856*x^6 + 15389/1228176*x^5 + 821/409392*x^4 - 22885/204696*x^3 + 3035/8529*x^2 + 3921/5686*x - 6122/2843',
  '275/39301632*x^11 - 245/26201088*x^10 + 77/6550272*x^9 + 125/1091712*x^8 + 2947/3275136*x^7 - 385/545856*x^6 + 15389/1228176*x^5 + 821/409392*x^4 - 22885/204696*x^3 + 3035/8529*x^2 + 9607/5686*x - 6122/2843',
  '275/39301632*x^11 + 245/26201088*x^10 + 77/6550272*x^9 - 125/1091712*x^8 + 2947/3275136*x^7 + 385/545856*x^6 + 15389/1228176*x^5 - 821/409392*x^4 - 22885/204696*x^3 - 3035/8529*x^2 + 3921/5686*x + 6122/2843',
  '275/39301632*x^11 + 245/26201088*x^10 + 77/6550272*x^9 - 125/1091712*x^8 + 2947/3275136*x^7 + 385/545856*x^6 + 15389/1228176*x^5 - 821/409392*x^4 - 22885/204696*x^3 - 3035/8529*x^2 + 9607/5686*x + 6122/2843',
  '275/19650816*x^11 + 77/3275136*x^9 + 2947/1637568*x^7 + 15389/614088*x^5 - 22885/102348*x^3 + 6764/2843*x',
];

/* ================================================================== */
describe('ZX arithmetic', () => {
  test('disc(x^n + a) = (-1)^(n(n-1)/2) n^n a^(n-1)', () => {
    for (const [n, a] of [
      [2, 1n],
      [3, -2n],
      [4, -2n],
      [5, 3n],
      [6, 108n],
      [7, 5n],
    ] as Array<[number, bigint]>) {
      const T: ZX = new Array(n + 1).fill(0n);
      T[0] = a;
      T[n] = 1n;
      const s = ((n * (n - 1)) / 2) % 2 === 0 ? 1n : -1n;
      expect(ZX_disc(T)).toBe(s * BigInt(n) ** BigInt(n) * a ** BigInt(n - 1));
    }
  });

  test('disc of some familiar polynomials', () => {
    expect(ZX_disc([1n, -3n, 0n, 1n])).toBe(81n); // the cyclic cubic
    expect(ZX_disc([1n, 1n, 1n, 1n, 1n])).toBe(125n); // polcyclo(5)
    expect(ZX_disc([-1n, 0n, 1n])).toBe(4n);
    expect(ZX_disc([0n, -1n, 0n, 1n])).toBe(4n); // x^3-x, squarefree
    expect(ZX_disc([0n, 0n, 1n])).toBe(0n); // x^2, not squarefree
  });

  test('resultant matches the product of B over the roots of A', () => {
    // res(x^2-2, x^2-3) = prod (a_i^2-3) = (2-3)(2-3) = 1
    expect(ZX_resultant([-2n, 0n, 1n], [-3n, 0n, 1n])).toBe(1n);
    // res(x-1, f) = f(1)
    expect(ZX_resultant([-1n, 1n], [5n, 7n, 1n])).toBe(13n);
    // res(A,B) = (-1)^(deg A deg B) res(B,A)
    const A: ZX = [3n, 1n, 4n, 1n];
    const B: ZX = [5n, 9n, 2n];
    expect(ZX_resultant(A, B)).toBe(ZX_resultant(B, A));
  });

  test('indexpartial divides a power of the discriminant', () => {
    expect(indexpartial([-5n, 0n, 1n])).toBe(10n); // disc = 20 = 2^2*5
    expect(indexpartial([-2n, 0n, 1n])).toBe(2n); // disc = 8
    const d = indexpartial([1n, -3n, 0n, 1n]);
    expect(d).toBe(9n); // PARI's galoisinit(x^3-3*x+1)[5] = 9
  });
});

describe('FpX layer', () => {
  const primes = [2n, 3n, 5n, 7n, 11n, 13n, 31n, 101n];

  test('factorisation mod p multiplies back and has the right degrees', () => {
    const T: ZX = [108n, 0n, 0n, 0n, 0n, 0n, 1n];
    for (const p of primes) {
      const Tp = FpX_red(T, p);
      if (FpX_degree(Tp) !== 6 || !FpX_is_squarefree(Tp, p)) continue;
      const F = FpX_factor_squarefree(Tp, p);
      let prod: bigint[] = [1n];
      for (const f of F) prod = FpX_mul(prod, f, p);
      expect(FpX_sub(prod, Tp, p).length).toBe(0);
      const { D, nb } = FpX_nbfact_by_degree(Tp, p);
      expect(nb).toBe(F.length);
      const count = new Map<number, number>();
      for (const f of F) count.set(FpX_degree(f), (count.get(FpX_degree(f)) ?? 0) + 1);
      for (const [d, c] of count) expect(D[d]).toBe(c);
    }
  });

  test('roots mod p agree with exhaustive search', () => {
    const pols: ZX[] = [
      [108n, 0n, 0n, 0n, 0n, 0n, 1n],
      [1n, 0n, 0n, 0n, 1n],
      [1n, -3n, 0n, 1n],
      [-2n, 0n, 1n],
    ];
    for (const T of pols)
      for (const p of primes) {
        const Tp = FpX_red(T, p);
        if (FpX_degree(Tp) < 1) continue;
        const brute: bigint[] = [];
        for (let i = 0n; i < p; i++) {
          let v = 0n;
          for (let k = Tp.length - 1; k >= 0; k--) v = (v * i + Tp[k]!) % p;
          if (v % p === 0n) brute.push(i);
        }
        // FpX_roots only returns the *distinct* roots
        expect(FpX_roots(Tp, p)).toEqual(brute);
      }
  });

  test('FpXQ_minpoly', () => {
    const p = 7n;
    const T = FpX_red([1n, 1n, 0n, 1n], p); // x^3+x+1, irreducible mod 7
    expect(FpXQ_minpoly([0n, 1n], T, p)).toEqual(T);
    expect(FpXQ_minpoly([3n], T, p)).toEqual([4n, 1n]); // x - 3
    // x^2 generates the same field but has its own minimal polynomial;
    // check the defining property rather than a hard-coded value
    const m = FpXQ_minpoly([0n, 0n, 1n], T, p);
    expect(FpX_degree(m)).toBe(3);
    let v: bigint[] = [];
    for (let i = m.length - 1; i >= 0; i--) {
      v = FpX_rem(FpX_mul(v, [0n, 0n, 1n], p), T, p);
      if (m[i] !== 0n) v = FpX_add(v, [m[i]!], p);
    }
    expect(v).toEqual([]);
    expect(m).toEqual([6n, 1n, 2n, 1n]); // y^3 + 2y^2 + y + 6
  });

  test('FpV_invVandermonde inverts the Vandermonde matrix', () => {
    const p = 10007n;
    const L = [0n, 3n, 11n, 42n, 100n, 7n];
    const den = 6n;
    const M = FpV_invVandermonde(L, den, p);
    for (let k = 1; k <= 5; k++)
      for (let j = 1; j <= 5; j++) {
        let s = 0n;
        for (let i = 1; i <= 5; i++) {
          let pw = 1n;
          for (let t = 1; t < j; t++) pw = (pw * L[i]!) % p;
          s = (s + M[k]![i]! * pw) % p;
        }
        expect(s % p).toBe(k === j ? den % p : 0n);
      }
  });
});

describe('Hensel lifting (Zp.c)', () => {
  const T: ZX = [108n, 0n, 0n, 0n, 0n, 0n, 1n];

  test('ZpX_roots lifts every root exactly', () => {
    const l = 31n; // totally split in Q[x]/(x^6+108) (PARI picks it too)
    const e = 8;
    const mod = l ** BigInt(e);
    const L = ZpX_roots(T, l, e);
    expect(L.length).toBe(7);
    expect(new Set(L.slice(1).map(String)).size).toBe(6);
    for (let i = 1; i <= 6; i++) {
      let v = 0n;
      for (let k = T.length - 1; k >= 0; k--) v = (v * L[i]! + T[k]!) % mod;
      expect(((v % mod) + mod) % mod).toBe(0n);
    }
  });

  test('ZpX_liftfact and bezout_lift_fact', () => {
    for (const p of [5n, 7n, 11n, 13n]) {
      if (!FpX_is_squarefree(FpX_red(T, p), p)) continue;
      const F = FpX_factor_squarefree(FpX_red(T, p), p);
      if (F.length < 2) continue;
      const Q = [[], ...F];
      const e = 4;
      const pe = p ** BigInt(e);
      const FL = ZpX_liftfact(T, Q, p, e);
      let prod: bigint[] = [1n];
      for (let i = 1; i < FL.length; i++) prod = FpX_mul(prod, FL[i]!, pe);
      expect(FpX_sub(prod, FpX_red(T, pe), pe).length).toBe(0);
      // each lifted factor reduces to the original mod p
      for (let i = 1; i < FL.length; i++)
        expect(FpX_sub(FpX_red(FL[i]!, p), F[i - 1]!, p).length).toBe(0);
      const U = bezout_lift_fact(T, Q, p, e);
      for (let i = 1; i < U.length; i++)
        for (let j = 1; j < FL.length; j++)
          expect(FpX_rem(U[i]!, FL[j]!, pe)).toEqual(i === j ? [1n] : []);
    }
  });
});

describe('permutations and groups (perm.c)', () => {
  test('perm_mul is composition: (s*t)[i] = s[t[i]]', () => {
    const s: Perm = [0, 2, 3, 1];
    const t: Perm = [0, 2, 1, 3];
    expect(perm_mul(s, t)).toEqual([0, 3, 2, 1]);
    expect(perm_mul(s, identity_perm(3))).toEqual(s);
    expect(perm_powu(s, 3)).toEqual(identity_perm(3));
  });

  test('perm_cycles and vecperm_orbits', () => {
    const s: Perm = [0, 2, 1, 4, 5, 3];
    expect(perm_cycles(s)).toEqual([[], [0, 1, 2], [0, 3, 4, 5]]);
    expect(vecperm_orbits([s], 5)).toEqual([[], [0, 1, 2], [0, 3, 4, 5]]);
    // the group generated by (12) and (345) has the same orbits
    expect(vecperm_orbits([[0, 2, 1, 3, 4, 5], [0, 1, 2, 4, 5, 3]], 5)).toEqual([
      [],
      [0, 1, 2],
      [0, 3, 4, 5],
    ]);
  });

  test('listznstarelts lists the subgroups of (Z/mZ)^*', () => {
    // (Z/8Z)^* = {1,3,5,7} = C2 x C2: subgroups of order dividing 4
    const l8 = listznstarelts(8, 4).map((g) => g.slice(1));
    expect(l8).toEqual([[1], [1, 3], [1, 5], [1, 7], [1, 3, 5, 7]]);
    // order dividing 2 only
    expect(listznstarelts(8, 2).map((g) => g.slice(1))).toEqual([[1], [1, 3], [1, 5], [1, 7]]);
    // (Z/5Z)^* is cyclic of order 4
    expect(listznstarelts(5, 4).map((g) => g.slice(1))).toEqual([[1], [1, 4], [1, 2, 3, 4]]);
    expect(listznstarelts(2, 1)).toEqual([[0, 1]]);
  });
});

/* ================================================================== */
describe('galoisinit', () => {
  for (const [name, o] of Object.entries(ORACLE)) {
    test(`${name}: agrees with PARI`, () => {
      const T = o.pol;
      const n = ZX_degree(T);
      const gal = galoisinit(T);
      expect(gal).not.toBeNull();
      const G = gal!;
      /* --- PARI's own values --- */
      expect(G.orders.slice(1)).toEqual(o.ord);
      expect(G.group.length - 1).toBe(n);
      expect(o.ord.reduce((a, b) => a * b, 1)).toBe(n);
      expect(galoisconj4(T)!.map(gpstr)).toEqual(o.conj);

      /* --- exact, non-p-adic verification --- */
      const elts = G.group.slice(1);
      for (const s of elts) expect(isPerm(s, n)).toBe(true);
      expect(new Set(elts.map(key)).size).toBe(n); // pairwise distinct
      const set = new Set(elts.map(key));
      for (const a of elts) for (const b of elts) expect(set.has(key(perm_mul(a, b)))).toBe(true);
      // every element really is an automorphism of Q[x]/(T)
      for (const s of elts) expect(vanishesAt(T, galoispermtopol(G, s), T)).toBe(true);
      // galoispermtopol and galoisvecpermtopol agree
      expect(galoisvecpermtopol(G, elts).map(gpstr)).toEqual(elts.map((s) => gpstr(galoispermtopol(G, s))));
      // the permutation really is the action on the l-adic roots
      for (const s of elts) {
        const sig = galoispermtopol(G, s);
        for (let i = 1; i <= n; i++) {
          let v = 0n;
          for (let k = sig.num.length - 1; k >= 0; k--) v = (v * G.roots[i]! + sig.num[k]!) % G.mod;
          // v / sig.den == roots[s[i]]
          expect(((v - sig.den * G.roots[s[i]!]!) % G.mod + G.mod) % G.mod).toBe(0n);
        }
      }
    });
  }

  test('x^12+96*x^8+...: A4 (a4galoisgen)', () => {
    const gal = galoisinit(A4_POL);
    expect(gal).not.toBeNull();
    expect(gal!.orders.slice(1)).toEqual([2, 2, 3]);
    expect(gal!.group.length - 1).toBe(12);
    expect(galoisconj4(A4_POL)!.map(gpstr)).toEqual(A4_CONJ);
    for (const s of gal!.group.slice(1))
      expect(vanishesAt(A4_POL, galoispermtopol(gal!, s), A4_POL)).toBe(true);
    expect(galoissubgroups(gal!).length).toBe(10);
    // PARI galoissubfields(galoisinit(T),1), verbatim
    expect(
      galoissubgroups(gal!)
        .map((H) => gpZX(galoisfixedfield(gal!, H, 1).P))
        .sort()
    ).toEqual(['x', 'x^12 + 96*x^8 + 1664*x^6 - 16128*x^4 + 165888*x^2 + 331776', 'x^3 - 768*x + 4096', 'x^4 + 384*x^2 + 4096*x + 36864', 'x^4 + 384*x^2 + 4096*x + 36864', 'x^4 + 512*x + 3072', 'x^4 - 512*x + 3072', 'x^6 + 384*x^4 + 13312*x^3 - 258048*x^2 + 5308416*x + 21233664', 'x^6 - 768*x^2 - 4096', 'x^6 - 768*x^2 - 4096']);
  });

  test('degree 20: the splitting field of x^5-2 (F20)', () => {
    const T: ZX = [50000n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 0n, 2500n];
    const P: ZX = new Array(21).fill(0n);
    P[0] = 50000n;
    P[10] = 2500n;
    P[20] = 1n;
    const gal = galoisinit(P);
    expect(gal).not.toBeNull();
    expect(gal!.orders.slice(1)).toEqual([5, 4]); // PARI: Vecsmall([5,4])
    expect(gal!.group.length - 1).toBe(20);
    expect(galoissubgroups(gal!).length).toBe(14); // PARI: #galoissubgroups = 14
    expect(galoisconj4(P)!.map(gpstr)).toEqual([
      '-x',
      'x',
      '-1/1100*x^11 - 18/11*x',
      '1/1100*x^11 + 18/11*x',
      '-3/22000*x^16 - 19/55*x^6 - 1/2*x',
      '-3/22000*x^16 - 19/55*x^6 + 1/2*x',
      '-1/11000*x^16 - 1/2200*x^11 - 47/220*x^6 - 29/22*x',
      '-1/11000*x^16 - 1/2200*x^11 - 47/220*x^6 - 7/22*x',
      '-1/11000*x^16 + 1/2200*x^11 - 47/220*x^6 + 7/22*x',
      '-1/11000*x^16 + 1/2200*x^11 - 47/220*x^6 + 29/22*x',
      '-1/22000*x^16 - 1/2200*x^11 - 29/220*x^6 - 9/11*x',
      '-1/22000*x^16 + 1/2200*x^11 - 29/220*x^6 + 9/11*x',
      '1/22000*x^16 - 1/2200*x^11 + 29/220*x^6 - 9/11*x',
      '1/22000*x^16 + 1/2200*x^11 + 29/220*x^6 + 9/11*x',
      '1/11000*x^16 - 1/2200*x^11 + 47/220*x^6 - 29/22*x',
      '1/11000*x^16 - 1/2200*x^11 + 47/220*x^6 - 7/22*x',
      '1/11000*x^16 + 1/2200*x^11 + 47/220*x^6 + 7/22*x',
      '1/11000*x^16 + 1/2200*x^11 + 47/220*x^6 + 29/22*x',
      '3/22000*x^16 + 19/55*x^6 - 1/2*x',
      '3/22000*x^16 + 19/55*x^6 + 1/2*x',
    ]);
    for (const s of gal!.gen.slice(1)) expect(vanishesAt(P, galoispermtopol(gal!, s), P)).toBe(true);
    void T;
  });

  test('an explicit denominator (PARI\'s second argument) is accepted', () => {
    const T = ORACLE['x^6+108']!.pol;
    // PARI: galoisinit(x^6+108)[5] = 648
    const gal = galoisinit(T, 648n);
    expect(gal).not.toBeNull();
    expect(gal!.den).toBe(648n);
    expect(gal!.orders.slice(1)).toEqual([3, 2]);
    for (const s of gal!.group.slice(1))
      expect(vanishesAt(T, galoispermtopol(gal!, s), T)).toBe(true);
  });

  test('the trivial extension', () => {
    const gal = galoisinit([0n, 1n]);
    expect(gal).not.toBeNull();
    expect(gal!.group).toEqual([[], [0, 1]]);
    expect(gal!.orders.slice(1)).toEqual([]);
  });

  test('non-Galois and reducible input: PARI returns 0, we return null', () => {
    // x^4-5x^2+6 = (x^2-2)(x^2-3): PARI galoisinit -> 0
    expect(galoisinit([6n, 0n, -5n, 0n, 1n])).toBeNull();
    // x^3-2 is irreducible but not Galois
    expect(galoisinit([-2n, 0n, 0n, 1n])).toBeNull();
    // x^4-2
    expect(galoisinit([-2n, 0n, 0n, 0n, 1n])).toBeNull();
  });

  test('argument checks mirror PARI', () => {
    expect(() => galoisinit([1n, 0n, 0n, 2n])).toThrow(PariImplError); // nonmonic
    expect(() => galoisinit([0n, 0n, 1n])).toThrow(PariDomainError); // x^2 not squarefree
  });

  test('S4 (degree 24) reports exactly what is missing', () => {
    const S4: ZX = [
      1215289321n, 0n, 2074608720n, 0n, 756753632n, 0n, -61067260n, 0n, 48052400n, 0n, -462400n,
      0n, 1241022n, 0n, -23120n, 0n, 7520n, 0n, 340n, 0n, 80n, 0n, 0n, 0n, 1n,
    ];
    expect(() => galoisinit(S4)).toThrow(NotImplementedError);
    expect(() => galoisinit(S4)).toThrow(/s4galoisgen/);
  });
});

/* ================================================================== */
describe('galoissubgroups', () => {
  for (const [name, o] of Object.entries(ORACLE)) {
    test(`${name}: ${o.nsub} subgroups, all genuine`, () => {
      const gal = galoisinit(o.pol)!;
      const n = ZX_degree(o.pol);
      const subs = galoissubgroups(gal);
      expect(subs.length).toBe(o.nsub);
      const seen = new Set<string>();
      for (const H of subs) {
        const elts = group_elts(H, n).slice(1);
        expect(elts.length).toBe(group_order(H));
        expect(n % group_order(H)).toBe(0); // Lagrange
        const set = new Set(elts.map(key));
        expect(set.size).toBe(group_order(H)); // distinct
        for (const a of elts)
          for (const b of elts) expect(set.has(key(perm_mul(a, b)))).toBe(true); // closed
        const sig = Array.from(set).sort().join('|');
        expect(seen.has(sig)).toBe(false); // no repeats
        seen.add(sig);
      }
      // the trivial group and the whole group are among them
      expect(subs.some((H) => group_order(H) === 1)).toBe(true);
      expect(subs.some((H) => group_order(H) === n)).toBe(true);
    });
  }
});

/* ================================================================== */
describe('galoisfixedfield', () => {
  for (const [name, o] of Object.entries(ORACLE)) {
    test(`${name}: the Galois correspondence`, () => {
      const T = o.pol;
      const n = ZX_degree(T);
      const gal = galoisinit(T)!;
      const subs = galoissubgroups(gal);
      const degs: number[] = [];
      for (const H of subs) {
        const { P, S } = galoisfixedfield(gal, H, 0);
        const d = ZX_degree(P);
        // Galois correspondence: [fixed field : Q] = [G:H]
        expect(d).toBe(n / group_order(H));
        degs.push(d);
        expect(ZX_disc(P)).not.toBe(0n); // separable
        // S is a root of P inside Q[x]/(T)
        expect(vanishesAt(P, S!, T)).toBe(true);
        // flag 1 returns the same polynomial
        expect(galoisfixedfield(gal, H, 1).P).toEqual(P);
      }
      expect(degs.sort((a, b) => a - b)).toEqual(o.degs);
    });
  }

  test('the subfield polynomials are PARI\'s galoissubfields(G,1), verbatim', () => {
    for (const [name, expected] of Object.entries(SUBFIELDS)) {
      const gal = galoisinit(ORACLE[name]!.pol)!;
      const mine = galoissubgroups(gal)
        .map((H) => gpZX(galoisfixedfield(gal, H, 1).P))
        .sort();
      expect(mine).toEqual(expected.slice().sort());
    }
  });

  test('a single permutation, a list of permutations and a Group all work', () => {
    const T = ORACLE['x^6+108']!.pol;
    const gal = galoisinit(T)!;
    const s = gal.gen[1]!; // of order 3
    const byPerm = galoisfixedfield(gal, s, 1).P;
    const byList = galoisfixedfield(gal, [s], 1).P;
    const byGroup = galoisfixedfield(gal, { gen: [[], s], ord: [0, 3] }, 1).P;
    expect(byList).toEqual(byPerm);
    expect(byGroup).toEqual(byPerm);
    expect(ZX_degree(byPerm)).toBe(2); // index 2
  });

  test('flag 2: the factorisation of T over the fixed field multiplies back to T', () => {
    for (const name of ['x^2-2', 'x^6+108', 'polcyclo(8)']) {
      const T = ORACLE[name]!.pol;
      const n = ZX_degree(T);
      const gal = galoisinit(T)!;
      for (const H of galoissubgroups(gal)) {
        const r = galoisfixedfield(gal, H, 2);
        const P = r.P;
        const m = ZX_degree(P);
        expect(r.factors!.length).toBe(m);
        for (const f of r.factors!) expect(f.length).toBe(n / m + 1);
        const prod = prodFactors(
          r.factors!.map((f) => f.map((c) => ({ num: c.num, den: c.den }))),
          P
        );
        expect(prod.length).toBe(n + 1);
        for (let i = 0; i <= n; i++) {
          // prod[i] must equal the constant T[i]
          const diff = qaddP(prod[i]!, { num: [-(T[i] ?? 0n)], den: 1n });
          expect(qisZero(diff)).toBe(true);
        }
      }
    }
  });
});
