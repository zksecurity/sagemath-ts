/**
 * @module parigp-ts/ffinit
 * @description Construction of irreducible polynomials over F_p (PARI's `ffinit`),
 * following the Adleman-Lenstra method.
 *
 * Direct port of `reference/pari/src/basemath/polarit3.c`:
 *
 * - `polarit3.c:3352-3358`  `ffinit_rand`
 * - `polarit3.c:3360-3388`  `ffinit_Artin_Schreier_2`
 * - `polarit3.c:3390-3410`  `ffinit_Artin_Schreier`
 * - `polarit3.c:3412-3424`  `flinit_check` / `flinit`
 * - `polarit3.c:3455-3499`  `fpinit_check` / `fpinit`  (Adleman-Lenstra variant)
 * - `polarit3.c:3426-3437`  `ffinit_fact_Flx` / `polarit3.c:3484-3496` `ffinit_fact`
 * - `polarit3.c:3439-3460`  `init_Flxq_i` / `polarit3.c:3498-3517` `init_Fq_i`, `ffinit`
 * - `polarit3.c:3524-3538`  `ffnbirred`
 * - `polarit3.c:1983-2005`  `FpX_composedsum`, `polarit3.c:2029-2038` `FpXV_composedsum`
 * - `subcyclo.c:424-447`    `polsubcyclo_cyclic` / `subcyclo.c:843-874` `polsubcyclo_g`
 *
 * Sage uses `pari(p).ffinit(n)` as the default modulus for `GF(p^n)` whenever no
 * Conway polynomial is available
 * (`reference/sage/src/sage/rings/polynomial/polynomial_ring.py`, `irreducible_element`,
 * algorithm `"first_lexicographic"`/default -> `pari(...).ffinit(...)`).
 *
 * DEVIATIONS from PARI's implementation (identical results, different route):
 *
 * 1. `FpX_composedsum` / the resultant used by `ffinit_Artin_Schreier`. PARI
 *    computes the composed sum through Newton sums and a truncated power-series
 *    product (`Flx_composedsum`, `Flx.c:4310-4340`), falling back to a p-adic lift
 *    (`ZpX_invLaplace_init`) when `p <= deg P * deg Q`, because the Laplace
 *    transform divides by factorials. We instead evaluate the *defining*
 *    resultant `Res_y(P(y), Q(x-y))` exactly, as the determinant of the
 *    multiplication-by-`Q(x-y)` operator on `F_p[x][y]/(P(y))` (Bareiss
 *    fraction-free elimination over the integral domain `F_p[x]`). Both compute
 *    `prod_{i,j} (x - a_i - b_j)`, so the output polynomial is identical; only
 *    the asymptotic complexity differs (the degrees involved in `ffinit` are
 *    bounded by `n`, so this is irrelevant here).
 *
 * 2. `polsubcyclo(n, l, 0)` for prime `n`. PARI computes the Gaussian periods
 *    numerically first to get a size bound and then p-adically
 *    (`subcyclo.c:polsubcyclo_g`). We compute the very same periods exactly in
 *    `Z[x]/(x^n - 1)` (or in `F_p[x]/(x^n - 1)` when only the reduction mod p is
 *    needed). Same polynomial, no floating point anywhere.
 *
 * 3. `FpX_is_irred`. PARI routes through its factorisation machinery
 *    (`FpX_factor.c:2384`); we use Rabin's irreducibility test. Same boolean.
 */

import { isPrime, factoru } from './ifactor.js';
import { PariDomainError, PariPrimeError } from './matkermod.js';

/* ------------------------------------------------------------------ */
/* Errors                                                              */
/* ------------------------------------------------------------------ */

/*
 * PARI error kinds are shared across the parigp-ts modules; they are currently
 * defined in `matkermod.ts` (see the note there) and re-exported here so that
 * `ffinit` callers can `import { PariDomainError } from './ffinit.js'`.
 */
export { PariDomainError, PariPrimeError } from './matkermod.js';

/* ------------------------------------------------------------------ */
/* FpX: dense polynomials over F_p                                     */
/* ------------------------------------------------------------------ */

/**
 * A polynomial over F_p: little-endian dense coefficient array, `f[i]` is the
 * coefficient of `x^i`, every entry in `[0, p)`, no trailing zeros.
 * The zero polynomial is `[]`.
 */
export type FpX = bigint[];

function mod(a: bigint, p: bigint): bigint {
  const r = a % p;
  return r < 0n ? r + p : r;
}

/** strip trailing zeros */
export function FpX_renormalize(f: FpX): FpX {
  let d = f.length;
  while (d > 0 && f[d - 1] === 0n) d--;
  return d === f.length ? f : f.slice(0, d);
}

/** reduce every coefficient mod p and renormalize */
export function FpX_red(f: readonly bigint[], p: bigint): FpX {
  return FpX_renormalize(f.map((c) => mod(c, p)));
}

export function FpX_degree(f: FpX): number {
  return f.length - 1;
}

export function FpX_add(a: FpX, b: FpX, p: bigint): FpX {
  const n = Math.max(a.length, b.length);
  const out: FpX = new Array(n);
  for (let i = 0; i < n; i++) out[i] = mod((a[i] ?? 0n) + (b[i] ?? 0n), p);
  return FpX_renormalize(out);
}

export function FpX_sub(a: FpX, b: FpX, p: bigint): FpX {
  const n = Math.max(a.length, b.length);
  const out: FpX = new Array(n);
  for (let i = 0; i < n; i++) out[i] = mod((a[i] ?? 0n) - (b[i] ?? 0n), p);
  return FpX_renormalize(out);
}

export function FpX_neg(a: FpX, p: bigint): FpX {
  return FpX_renormalize(a.map((c) => mod(-c, p)));
}

export function FpX_mul(a: FpX, b: FpX, p: bigint): FpX {
  if (a.length === 0 || b.length === 0) return [];
  const out: FpX = new Array(a.length + b.length - 1).fill(0n);
  for (let i = 0; i < a.length; i++) {
    const ai = a[i]!;
    if (ai === 0n) continue;
    for (let j = 0; j < b.length; j++) out[i + j] = (out[i + j]! + ai * b[j]!) % p;
  }
  return FpX_renormalize(out);
}

/** multiply by a scalar */
export function FpX_Fp_mul(a: FpX, c: bigint, p: bigint): FpX {
  const cc = mod(c, p);
  if (cc === 0n) return [];
  return FpX_renormalize(a.map((x) => (x * cc) % p));
}

function Fp_inv(a: bigint, p: bigint): bigint {
  let [old_r, r] = [mod(a, p), p];
  let [old_s, s] = [1n, 0n];
  while (r !== 0n) {
    const q = old_r / r;
    [old_r, r] = [r, old_r - q * r];
    [old_s, s] = [s, old_s - q * s];
  }
  if (old_r !== 1n) throw new Error(`impossible inverse modulo ${p}: ${a}`);
  return mod(old_s, p);
}

/** Euclidean division: returns `[q, r]` with `a = q*b + r`, `deg r < deg b`. */
export function FpX_divrem(a: FpX, b: FpX, p: bigint): [FpX, FpX] {
  if (b.length === 0) throw new Error('FpX_divrem: division by zero');
  if (a.length < b.length) return [[], a.slice()];
  const inv = Fp_inv(b[b.length - 1]!, p);
  const r = a.slice();
  const q: FpX = new Array(a.length - b.length + 1).fill(0n);
  for (let i = a.length - b.length; i >= 0; i--) {
    const c = (r[i + b.length - 1]! * inv) % p;
    q[i] = c;
    if (c === 0n) continue;
    for (let j = 0; j < b.length; j++)
      r[i + j] = mod(r[i + j]! - c * b[j]!, p);
  }
  return [FpX_renormalize(q), FpX_renormalize(r.slice(0, b.length - 1))];
}

export function FpX_rem(a: FpX, b: FpX, p: bigint): FpX {
  return FpX_divrem(a, b, p)[1];
}

/** make monic */
export function FpX_normalize(f: FpX, p: bigint): FpX {
  if (f.length === 0) return f;
  const lc = f[f.length - 1]!;
  if (lc === 1n) return f;
  return FpX_Fp_mul(f, Fp_inv(lc, p), p);
}

export function FpX_gcd(a: FpX, b: FpX, p: bigint): FpX {
  let x = FpX_renormalize(a.slice());
  let y = FpX_renormalize(b.slice());
  while (y.length !== 0) {
    const r = FpX_rem(x, y, p);
    x = y;
    y = r;
  }
  return x.length === 0 ? x : FpX_normalize(x, p);
}

/** `x^n` as an FpX */
export function pol_xn(n: number): FpX {
  const out: FpX = new Array(n + 1).fill(0n);
  out[n] = 1n;
  return out;
}

export function FpXQ_mul(a: FpX, b: FpX, T: FpX, p: bigint): FpX {
  return FpX_rem(FpX_mul(a, b, p), T, p);
}

/** `x^e mod (T, p)` for `e >= 0` */
export function FpXQ_pow(a: FpX, e: bigint, T: FpX, p: bigint): FpX {
  if (e < 0n) throw new Error('FpXQ_pow: negative exponent');
  let result: FpX = FpX_rem([1n], T, p);
  let base = FpX_rem(a, T, p);
  let k = e;
  while (k > 0n) {
    if (k & 1n) result = FpXQ_mul(result, base, T, p);
    base = FpXQ_mul(base, base, T, p);
    k >>= 1n;
  }
  return result;
}

/**
 * Rabin's irreducibility test for a non-constant `f` over F_p.
 * @see the module-level deviation note (PARI uses `FpX_factor_i(f,p,2)`).
 */
export function FpX_is_irred(f: FpX, p: bigint): boolean {
  const g = FpX_normalize(FpX_red(f, p), p);
  const d = FpX_degree(g);
  if (d <= 0) return false;
  if (d === 1) return true;
  // x^(p^d) == x (mod g)
  let xq = pol_xn(1);
  for (let i = 0; i < d; i++) xq = FpXQ_pow(xq, p, g, p);
  if (FpX_sub(xq, pol_xn(1), p).length !== 0) return false;
  // gcd(x^(p^(d/r)) - x, g) == 1 for every prime r | d
  const primes = factoru(BigInt(d)).map(([q]) => Number(q));
  for (const r of primes) {
    let y = pol_xn(1);
    for (let i = 0; i < d / r; i++) y = FpXQ_pow(y, p, g, p);
    const h = FpX_gcd(FpX_sub(y, pol_xn(1), p), g, p);
    if (FpX_degree(h) !== 0) return false;
  }
  return true;
}

/**
 * PARI `ffinit_rand` (`polarit3.c:3352-3358`): a random monic irreducible
 * polynomial of degree `n` over F_p. Not used by {@link ffinit} itself (PARI
 * doesn't either); provided for completeness.
 */
export function ffinit_rand(p: bigint, n: number, rnd: () => number = Math.random): FpX {
  for (;;) {
    const pol = pol_xn(n);
    for (let i = 0; i < n; i++) pol[i] = BigInt(Math.floor(rnd() * Number(p))) % p;
    if (FpX_is_irred(pol, p)) return FpX_renormalize(pol);
  }
}

/* ------------------------------------------------------------------ */
/* Bivariate resultant Res_y(T(y), Q(x,y))                             */
/* ------------------------------------------------------------------ */

/**
 * A polynomial in `F_p[x][y]`: `Q[k]` is the coefficient of `y^k`, itself an
 * {@link FpX} in `x`.
 */
export type FpXY = FpX[];

function FpXY_renormalize(Q: FpXY): FpXY {
  let d = Q.length;
  while (d > 0 && Q[d - 1]!.length === 0) d--;
  return Q.slice(0, d);
}

/** exact division in F_p[x]; throws if the division is not exact */
function FpX_divexact(a: FpX, b: FpX, p: bigint): FpX {
  const [q, r] = FpX_divrem(a, b, p);
  if (r.length !== 0) throw new Error('FpX_divexact: inexact division');
  return q;
}

/**
 * `Res_y(T(y), Q(x,y))` for `T` **monic** in `F_p[y]`, returned as an
 * {@link FpX} in `x` of degree `deg(T) * deg_x(Q)`.
 *
 * This is exactly PARI's `Flx_FlxY_resultant(T, Q, p)`
 * (`polarit3.c:1916-1930`): `T` is the polynomial in the eliminated variable and
 * `Q` is bivariate; the result lives in the remaining variable.
 *
 * Since `T` is monic, `Res_y(T,Q) = prod_{T(a)=0} Q(x,a) = det(mult by Q on
 * F_p[x][y]/(T(y)))`. We build that matrix and take its determinant with
 * fraction-free (Bareiss) elimination over the integral domain `F_p[x]`.
 */
export function FpX_FpXY_resultant(T: FpX, Q: FpXY, p: bigint): FpX {
  const d = FpX_degree(T);
  if (d < 0) throw new Error('FpX_FpXY_resultant: T = 0');
  if (T[d] !== 1n) throw new Error('FpX_FpXY_resultant: T must be monic');
  if (d === 0) return [1n];
  const G = FpXY_renormalize(Q.map((c) => FpX_renormalize(c.slice())));
  if (G.length === 0) return [];

  // reduce a bivariate poly modulo T(y) (T monic, constant in x)
  const redmodT = (A: FpXY): FpXY => {
    const R = A.map((c) => c.slice());
    for (let k = R.length - 1; k >= d; k--) {
      const c = R[k]!;
      if (c.length === 0) continue;
      R[k] = [];
      for (let i = 0; i < d; i++) {
        if (T[i] === 0n) continue;
        R[k - d + i] = FpX_sub(R[k - d + i]!, FpX_Fp_mul(c, T[i]!, p), p);
      }
    }
    R.length = Math.min(R.length, d);
    while (R.length < d) R.push([]);
    return R;
  };

  // columns of the multiplication-by-G matrix: y^j * G mod T
  const cols: FpXY[] = [];
  let cur = redmodT(G);
  cols.push(cur);
  for (let j = 1; j < d; j++) {
    const shifted: FpXY = [[], ...cur.map((c) => c.slice())];
    cur = redmodT(shifted);
    cols.push(cur);
  }

  // M[i][j] = coefficient of y^i in column j
  const M: FpX[][] = [];
  for (let i = 0; i < d; i++) {
    const row: FpX[] = [];
    for (let j = 0; j < d; j++) row.push(cols[j]![i] ?? []);
    M.push(row);
  }

  /* Bareiss fraction-free determinant over F_p[x] */
  let sign = 1n;
  let prev: FpX = [1n];
  for (let k = 0; k < d - 1; k++) {
    if (M[k]![k]!.length === 0) {
      let r = -1;
      for (let i = k + 1; i < d; i++)
        if (M[i]![k]!.length !== 0) {
          r = i;
          break;
        }
      if (r < 0) return [];
      const t = M[k]!;
      M[k] = M[r]!;
      M[r] = t;
      sign = mod(-sign, p);
    }
    for (let i = k + 1; i < d; i++)
      for (let j = k + 1; j < d; j++) {
        const num = FpX_sub(
          FpX_mul(M[i]![j]!, M[k]![k]!, p),
          FpX_mul(M[i]![k]!, M[k]![j]!, p),
          p
        );
        M[i]![j] = FpX_divexact(num, prev, p);
      }
    prev = M[k]![k]!;
  }
  return FpX_Fp_mul(M[d - 1]![d - 1]!, sign, p);
}

/**
 * PARI `FpX_composedsum(P, Q, p)` (`polarit3.c:1983-2005`): the monic (up to
 * leading coefficients) polynomial whose roots are all sums `a_i + b_j` of a
 * root of `P` and a root of `Q`; equal to `Res_y(P(y), Q(x-y))`.
 *
 * @see module deviation note 1 for the algorithmic difference.
 */
export function FpX_composedsum(P: FpX, Q: FpX, p: bigint): FpX {
  const dP = FpX_degree(P);
  const dQ = FpX_degree(Q);
  if (dP < 0 || dQ < 0) return [];
  if (dP === 0 || dQ === 0) return [1n];
  const lead = (Fp_pow(P[dP]!, BigInt(dQ), p) * Fp_pow(Q[dQ]!, BigInt(dP), p)) % p;
  // G(x,y) = Q(x-y) expanded as a polynomial in y with FpX coefficients
  const G: FpXY = new Array(dQ + 1).fill(null).map(() => [] as FpX);
  // (x-y)^k = sum_j C(k,j) x^(k-j) (-y)^j
  const binom: bigint[][] = [];
  for (let k = 0; k <= dQ; k++) {
    binom.push(new Array(k + 1).fill(0n));
    binom[k]![0] = 1n;
    for (let j = 1; j <= k; j++)
      binom[k]![j] = mod((binom[k - 1]![j - 1] ?? 0n) + (binom[k - 1]![j] ?? 0n), p);
  }
  for (let k = 0; k <= dQ; k++) {
    const qk = Q[k] ?? 0n;
    if (qk === 0n) continue;
    for (let j = 0; j <= k; j++) {
      const sgn = j % 2 === 0 ? 1n : p - 1n;
      const c = (((qk * binom[k]![j]!) % p) * sgn) % p;
      if (c === 0n) continue;
      const term = pol_xn(k - j).map((v) => (v * c) % p);
      G[j] = FpX_add(G[j]!, FpX_renormalize(term), p);
    }
  }
  const Pm = FpX_normalize(P, p);
  const R = FpX_FpXY_resultant(Pm, G, p);
  return FpX_Fp_mul(R, lead, p);
}

function Fp_pow(a: bigint, e: bigint, p: bigint): bigint {
  let r = 1n;
  let b = mod(a, p);
  let k = e;
  while (k > 0n) {
    if (k & 1n) r = (r * b) % p;
    b = (b * b) % p;
    k >>= 1n;
  }
  return r;
}

/**
 * PARI `FpXV_composedsum(V, p)` (`polarit3.c:2029-2038`): iterated composed sum
 * over a vector of polynomials, computed with a product tree
 * (`bb_group.c:355-377` `gen_product`).
 */
export function FpXV_composedsum(V: readonly FpX[], p: bigint): FpX {
  if (V.length === 0) return [1n];
  let cur = V.map((f) => f.slice());
  while (cur.length > 1) {
    const next: FpX[] = [];
    for (let i = 0; i + 1 < cur.length; i += 2)
      next.push(FpX_composedsum(cur[i]!, cur[i + 1]!, p));
    if (cur.length % 2 === 1) next.push(cur[cur.length - 1]!);
    cur = next;
  }
  return cur[0]!;
}

/* ------------------------------------------------------------------ */
/* polsubcyclo for prime conductor                                     */
/* ------------------------------------------------------------------ */

/** multiplicative order of `a` mod `n`, knowing it divides `o` (PARI `Fl_order`) */
function Fl_order(a: bigint, o: bigint, n: bigint): bigint {
  let ord = o;
  for (const [r] of factoru(o)) {
    while (ord % r === 0n && Fp_pow(a, ord / r, n) === 1n) ord /= r;
  }
  return ord;
}

/** smallest primitive root mod the prime `n` */
function primitiveRoot(n: bigint): bigint {
  if (n === 2n) return 1n;
  const fac = factoru(n - 1n).map(([q]) => q);
  for (let g = 2n; g < n; g++) {
    let ok = true;
    for (const q of fac)
      if (Fp_pow(g, (n - 1n) / q, n) === 1n) {
        ok = false;
        break;
      }
    if (ok) return g;
  }
  throw new Error(`no primitive root mod ${n}`);
}

/**
 * `polsubcyclo(n, l, 0)` for a **prime** `n` with `l | n-1`: the degree-`l`
 * subfield of `Q(zeta_n)`, i.e. the minimal polynomial of the Gaussian periods.
 *
 * Port of `subcyclo.c:843-874` `polsubcyclo_g` restricted to prime conductor
 * (the only case `fpinit` needs, since `fpinit_check` requires `n` prime), with
 * `polsubcyclo_cyclic` (`subcyclo.c:424-447`) providing the periods
 *   eta_i = sum_{k=0}^{m-1} zeta_n^(g^i * g^(l*k)),  m = (n-1)/l,
 * and the answer being `prod_{i=0}^{l-1} (X - eta_i)`.
 *
 * The result does not depend on the choice of primitive root `g`: the periods
 * are indexed by the `l` cosets of the unique index-`l` subgroup of (Z/nZ)^*.
 *
 * @param modulus if given, the computation is carried out in `F_p[x]/(x^n-1)`
 *                and the result is returned reduced mod `p` (this is PARI's
 *                `fpinit`, which only ever needs `FpX_red(polsubcyclo(n,l,0),p)`).
 * @returns the coefficients, little-endian; without `modulus` these are the true
 *          (possibly negative) integers, matching PARI's `FpX_center`ed output.
 */
export function polsubcyclo_prime(n: number, l: number, modulus?: bigint): bigint[] {
  const N = BigInt(n);
  if (!isPrime(N)) throw new PariDomainError('polsubcyclo_prime', 'n', 'is not', 'prime');
  if ((n - 1) % l !== 0) throw new PariDomainError('polsubcyclo_prime', 'l', 'does not divide', 'n-1');
  const m = (n - 1) / l;
  const g = primitiveRoot(N);
  const gl = Fp_pow(g, BigInt(l), N);

  // element of Z[x]/(x^n - 1) (or F_p[x]/(x^n-1)); index = exponent of zeta
  const red = modulus ? (v: bigint): bigint => mod(v, modulus) : (v: bigint): bigint => v;
  const mulRing = (a: bigint[], b: bigint[]): bigint[] => {
    const out = new Array<bigint>(n).fill(0n);
    for (let i = 0; i < n; i++) {
      const ai = a[i]!;
      if (ai === 0n) continue;
      for (let j = 0; j < n; j++) {
        const bj = b[j]!;
        if (bj === 0n) continue;
        const k = i + j >= n ? i + j - n : i + j;
        out[k] = red(out[k]! + ai * bj);
      }
    }
    return out;
  };
  const subRing = (a: bigint[], b: bigint[]): bigint[] =>
    a.map((x, i) => red(x - b[i]!));

  // the l Gaussian periods
  const etas: bigint[][] = [];
  let base = 1n; // g^i
  for (let i = 0; i < l; i++) {
    const e = new Array<bigint>(n).fill(0n);
    let ex = base;
    for (let k = 0; k < m; k++) {
      const idx = Number(ex);
      e[idx] = red(e[idx]! + 1n);
      ex = (ex * gl) % N;
    }
    etas.push(e);
    base = (base * g) % N;
  }

  // prod (X - eta_i) in (Z[x]/(x^n-1))[X]; coefficients end up rational integers
  let poly: bigint[][] = [new Array<bigint>(n).fill(0n)];
  poly[0]![0] = 1n; // constant polynomial 1
  for (const e of etas) {
    const next: bigint[][] = new Array(poly.length + 1);
    for (let i = 0; i <= poly.length; i++) next[i] = new Array<bigint>(n).fill(0n);
    for (let i = 0; i < poly.length; i++) {
      // X * coeff
      next[i + 1] = poly[i]!.slice();
      // - eta * coeff
      next[i] = subRing(next[i]!, mulRing(poly[i]!, e));
    }
    poly = next;
  }

  // Each coefficient is a rational integer: reduce modulo Phi_n = 1+x+...+x^(n-1)
  // (i.e. kill the x^(n-1) coefficient) and read off the constant term.
  const out: bigint[] = [];
  for (const c of poly) {
    const t = c[n - 1]!;
    for (let i = 0; i < n; i++) c[i] = red(c[i]! - t);
    for (let i = 1; i < n - 1; i++)
      if (c[i] !== 0n)
        throw new Error('polsubcyclo_prime: period polynomial is not rational');
    out.push(c[0]!);
  }
  return out;
}

/* ------------------------------------------------------------------ */
/* Adleman-Lenstra: fpinit                                             */
/* ------------------------------------------------------------------ */

/**
 * PARI `fpinit_check` (`polarit3.c:3455-3462`) / `flinit_check`
 * (`polarit3.c:3412-3419`): is `polsubcyclo(n,l,0)` irreducible mod `p`?
 */
export function fpinit_check(p: bigint, n: number, l: number): boolean {
  const N = BigInt(n);
  if (!isPrime(N)) return false;
  const q = mod(p, N);
  if (q === 0n) return false;
  const ord = Fl_order(q, N - 1n, N);
  return gcdInt((N - 1n) / ord, BigInt(l)) === 1n;
}

function gcdInt(a: bigint, b: bigint): bigint {
  a = a < 0n ? -a : a;
  b = b < 0n ? -b : b;
  while (b) {
    const t = a % b;
    a = b;
    b = t;
  }
  return a;
}

/**
 * PARI `fpinit` (`polarit3.c:3464-3477`): an irreducible polynomial of degree
 * `l` over F_p, obtained as a subfield of a cyclotomic field. This is the
 * variant of Adleman-Lenstra, "Finding irreducible polynomials over finite
 * fields", ACM 1986 (5) 350--355, described in PARI's comment.
 *
 * PRECONDITION (PARI's comment: "assume k*p does not divide l", k = 2 if
 * p % 4 == 1 else 4): in particular `p` must not divide `l`. PARI's search loop
 * would spin forever otherwise -- e.g. `p = 2, l = 8`: any prime `n = 1 mod 8`
 * makes 2 a quadratic residue, so `(n-1)/ord_n(2)` is always even and
 * `fpinit_check` never succeeds. `ffinit` only ever calls this with `l` a power
 * of a prime different from `p` (the `p`-part goes through
 * {@link ffinit_Artin_Schreier}). We raise instead of looping.
 */
export function fpinit(p: bigint, l: number): FpX {
  if (BigInt(l) % p === 0n)
    throw new PariDomainError('fpinit', 'p', 'divides', `l = ${l}`);
  let n = 1 + l;
  while (!fpinit_check(p, n, l)) n += l;
  return FpX_red(polsubcyclo_prime(n, l, p), p);
}

/* ------------------------------------------------------------------ */
/* Artin-Schreier towers                                               */
/* ------------------------------------------------------------------ */

/**
 * PARI `ffinit_Artin_Schreier_2` (`polarit3.c:3360-3388`): an extension of
 * degree `2^l` of F_2, `l > 0`.
 */
export function ffinit_Artin_Schreier_2(l: number): FpX {
  if (l <= 0) throw new PariDomainError('ffinit_Artin_Schreier_2', 'l', '<=', '0');
  if (l === 1) return [1n, 1n, 1n]; /* x^2 + x + 1 */
  /* Q = x^2 + x + y(y^2+y), as a polynomial in y with F_2[x] coefficients */
  const Q: FpXY = [[0n, 1n, 1n], [], [1n], [1n]];
  /* T = x^4 + x + 1, minimal polynomial of a root of Q */
  let T: FpX = [1n, 1n, 0n, 0n, 1n];
  for (let i = 2; i < l; i++) T = FpX_FpXY_resultant(T, Q, 2n);
  return T;
}

/**
 * PARI `ffinit_Artin_Schreier` (`polarit3.c:3390-3410`): an extension of degree
 * `p^l` of F_p, `l > 0`.
 */
export function ffinit_Artin_Schreier(p: bigint, l: number): FpX {
  if (p === 2n) return ffinit_Artin_Schreier_2(l);
  if (l <= 0) throw new PariDomainError('ffinit_Artin_Schreier', 'l', '<=', '0');
  const pn = Number(p);
  /* T = x^p - x - 1 */
  let T: FpX = FpX_sub(pol_xn(pn), [1n, 1n], p);
  if (l === 1) return T;
  /* Q = x^p - x - (y^(2p-1) - y^p), as a polynomial in y over F_p[x] */
  const Q: FpXY = new Array(2 * pn).fill(null).map(() => [] as FpX);
  Q[0] = FpX_sub(pol_xn(pn), pol_xn(1), p);
  Q[pn] = [1n];
  Q[2 * pn - 1] = [mod(-1n, p)];
  for (let i = 2; i <= l; i++) T = FpX_FpXY_resultant(T, Q, p);
  return T;
}

/* ------------------------------------------------------------------ */
/* ffinit                                                              */
/* ------------------------------------------------------------------ */

/** PARI `factoru_pow` (`arith2.c:81`): `[primes, exponents, prime powers]` */
function factoru_pow(n: number): Array<[bigint, number, number]> {
  return factoru(BigInt(n)).map(([q, e]) => [q, Number(e), Number(q ** e)]);
}

/**
 * PARI `ffinit_fact` (`polarit3.c:3484-3496`) / `ffinit_fact_Flx`
 * (`polarit3.c:3426-3437`).
 */
function ffinit_fact(p: bigint, n: number): FpX {
  const F = factoru_pow(n);
  const P: FpX[] = F.map(([q, e, m]) =>
    p === q ? ffinit_Artin_Schreier(q, e) : fpinit(p, m)
  );
  return FpXV_composedsum(P, p);
}

/**
 * PARI `init_Fq` (`polarit3.c:3498-3519`) — the un-`Mod`ded version of
 * {@link ffinit}: a monic irreducible polynomial of degree `n` over F_p.
 *
 * PARI's `ffinit(p,n,v)` is `FpX_to_mod(init_Fq_i(p,n,v), p)`, i.e. the same
 * polynomial with `Mod(.,p)` coefficients; we return the plain coefficient
 * vector, so {@link ffinit} and {@link init_Fq} agree.
 */
export function init_Fq(p: bigint, n: number): FpX {
  if (n <= 0) throw new PariDomainError('ffinit', 'degree', '<=', '0');
  if (p < 2n) throw new PariPrimeError('ffinit', p);
  if (n === 1) return pol_xn(1);
  if (p !== 2n && p % 2n === 0n) throw new PariPrimeError('ffinit', p);
  /* if 1 + x + ... + x^n is irreducible mod p, use it (polcyclo(n+1)) */
  if (fpinit_check(p, n + 1, n)) return new Array(n + 1).fill(1n) as FpX;
  return ffinit_fact(p, n);
}

/**
 * PARI `ffinit(p, n)` (`polarit3.c:3520-3525`): a monic irreducible polynomial
 * of degree `n` over F_p, returned as an {@link FpX} (little-endian
 * coefficients in `[0,p)`). PARI wraps the same polynomial in `Mod(.,p)`.
 *
 * This is the default modulus Sage uses for `GF(p^n)` when no Conway polynomial
 * is available.
 */
export function ffinit(p: bigint, n: number): FpX {
  if (!isPrime(p)) throw new PariPrimeError('ffinit', p);
  return init_Fq(p, n);
}

/**
 * PARI `ffnbirred(p, n)` (`polarit3.c:3527-3540`): the number of monic
 * irreducible polynomials of degree `n` over F_p.
 */
export function ffnbirred(p: bigint, n: number): bigint {
  if (n <= 0) throw new PariDomainError('ffnbirred', 'degree', '<=', '0');
  // s = sum_{d | n} mu(d) p^(n/d), divided by n
  let s = 0n;
  const primes = factoru(BigInt(n)).map(([q]) => Number(q));
  const sub = (mask: number): [number, number] => {
    let d = 1;
    let mu = 1;
    for (let i = 0; i < primes.length; i++)
      if (mask & (1 << i)) {
        d *= primes[i]!;
        mu = -mu;
      }
    return [d, mu];
  };
  for (let mask = 0; mask < 1 << primes.length; mask++) {
    const [d, mu] = sub(mask);
    const pd = p ** BigInt(n / d);
    s += mu > 0 ? pd : -pd;
  }
  return s / BigInt(n);
}
