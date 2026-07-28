/**
 * @module sage/rings/number_field/pari_nf
 * @description The number-field kernel routines that SageMath delegates to PARI/GP.
 *
 * SageMath's `sage.rings.number_field` reaches PARI through cypari2 for
 *
 * - `nfbasis` / `nfdisc`  (`reference/pari/src/basemath/base2.c:maxord`),
 * - `idealprimedec`       (`reference/pari/src/basemath/base2.c:idealprimedec`),
 * - `nfgaloisconj`        (`reference/pari/src/basemath/galconj.c`),
 * - `polisirreducible`    (`reference/pari/src/basemath/polarit2.c`).
 *
 * `parigp-ts` currently has no `nf` module, so this file ports the underlying
 * algorithms directly:
 *
 * - `nfbasis`/`nfdisc`: the Pohst--Zassenhaus *Round 2* p-maximal order
 *   algorithm (Cohen, *A Course in Computational Algebraic Number Theory*,
 *   Algorithm 6.1.8), which is exactly what PARI's `maxord` runs when the
 *   Round-4 shortcut does not apply.
 * - `idealprimedec`: the Dedekind--Kummer theorem applied to `f mod p`
 *   (valid whenever `p` does not divide `[O_K : Z[theta]]`, which is the
 *   `p_2` branch of PARI's `idealprimedec`).
 * - `nfgaloisconj`: p-adic reconstruction of the conjugates of `theta`
 *   (PARI's "Allombert" method, with rational reconstruction in place of LLL).
 * - polynomial factorisation over `F_p` (Cantor--Zassenhaus) and over `Z`
 *   (Zassenhaus with a single big prime, so no Hensel lift is required).
 *
 * Everything here is exact integer/rational arithmetic.
 *
 * @see Deviation: Number Field Kernel Ported Locally Instead of parigp-ts
 */

import { factor as intFactor, is_prime, isqrt, xgcd } from '../../arith/misc.js';
import { NotImplementedError, ValueError } from '../../errors.js';
import { IntegerMatrix, LLL } from '../../matrix/matrix_integer.js';
import { Rational } from '../rational.js';

/** A dense polynomial with integer coefficients, `coeffs[i]` is the coefficient of `x^i`. */
export type ZPoly = bigint[];

// ---------------------------------------------------------------------------
// Small integer helpers
// ---------------------------------------------------------------------------

function babs(a: bigint): bigint {
  return a < 0n ? -a : a;
}

function bgcd(a: bigint, b: bigint): bigint {
  a = babs(a);
  b = babs(b);
  while (b) {
    [a, b] = [b, a % b];
  }
  return a;
}

/** Euclidean remainder in `[0, m)`. */
function mmod(a: bigint, m: bigint): bigint {
  const r = a % m;
  return r < 0n ? r + m : r;
}

/** Balanced representative of `a mod m`, in `(-m/2, m/2]`. */
function balanced(a: bigint, m: bigint): bigint {
  const r = mmod(a, m);
  return r * 2n > m ? r - m : r;
}

/** Modular inverse of `a` modulo `m` (requires `gcd(a, m) = 1`). */
export function invMod(a: bigint, m: bigint): bigint {
  const [g, u] = xgcd(mmod(a, m), m);
  if (g !== 1n) {
    throw new ValueError(`inverse of ${a} modulo ${m} does not exist`);
  }
  return mmod(u, m);
}

function powMod(base: bigint, exp: bigint, m: bigint): bigint {
  let r = 1n;
  let b = mmod(base, m);
  let e = exp;
  while (e > 0n) {
    if (e & 1n) r = (r * b) % m;
    b = (b * b) % m;
    e >>= 1n;
  }
  return r;
}

// ---------------------------------------------------------------------------
// Integer polynomial arithmetic
// ---------------------------------------------------------------------------

export function zpNorm(f: ZPoly): ZPoly {
  let l = f.length;
  while (l > 0 && f[l - 1] === 0n) l--;
  return f.slice(0, l);
}

export function zpDeg(f: ZPoly): number {
  return zpNorm(f).length - 1;
}

export function zpAdd(a: ZPoly, b: ZPoly): ZPoly {
  const n = Math.max(a.length, b.length);
  const r: ZPoly = new Array(n).fill(0n);
  for (let i = 0; i < n; i++) r[i] = (a[i] ?? 0n) + (b[i] ?? 0n);
  return zpNorm(r);
}

export function zpSub(a: ZPoly, b: ZPoly): ZPoly {
  const n = Math.max(a.length, b.length);
  const r: ZPoly = new Array(n).fill(0n);
  for (let i = 0; i < n; i++) r[i] = (a[i] ?? 0n) - (b[i] ?? 0n);
  return zpNorm(r);
}

export function zpMul(a: ZPoly, b: ZPoly): ZPoly {
  const A = zpNorm(a);
  const B = zpNorm(b);
  if (A.length === 0 || B.length === 0) return [];
  const r: ZPoly = new Array(A.length + B.length - 1).fill(0n);
  for (let i = 0; i < A.length; i++) {
    if (A[i] === 0n) continue;
    for (let j = 0; j < B.length; j++) {
      r[i + j] = r[i + j]! + A[i]! * B[j]!;
    }
  }
  return zpNorm(r);
}

export function zpScale(a: ZPoly, c: bigint): ZPoly {
  return zpNorm(a.map((x) => x * c));
}

/** Remainder of `a` on division by the *monic* integer polynomial `g`. */
export function zpRemMonic(a: ZPoly, g: ZPoly): ZPoly {
  const G = zpNorm(g);
  const dg = G.length - 1;
  if (dg < 0) throw new ValueError('division by the zero polynomial');
  if (G[dg] !== 1n) throw new ValueError('zpRemMonic requires a monic modulus');
  const r = [...zpNorm(a)];
  for (let i = r.length - 1; i >= dg; i--) {
    const c = r[i]!;
    if (c === 0n) continue;
    r[i] = 0n;
    for (let j = 0; j < dg; j++) {
      r[i - dg + j] = (r[i - dg + j] ?? 0n) - c * G[j]!;
    }
  }
  return zpNorm(r);
}

/** Product of `a` and `b` reduced modulo the monic integer polynomial `g`. */
export function zpMulMod(a: ZPoly, b: ZPoly, g: ZPoly): ZPoly {
  return zpRemMonic(zpMul(a, b), g);
}

export function zpDerivative(a: ZPoly): ZPoly {
  const r: ZPoly = [];
  for (let i = 1; i < a.length; i++) r.push(a[i]! * BigInt(i));
  return zpNorm(r);
}

/** Exact division `a / b` of integer polynomials, or `null` when it is not exact. */
export function zpExactDiv(a: ZPoly, b: ZPoly): ZPoly | null {
  const A = [...zpNorm(a)];
  const B = zpNorm(b);
  if (B.length === 0) return null;
  const db = B.length - 1;
  const lb = B[db]!;
  const q: ZPoly = new Array(Math.max(0, A.length - db)).fill(0n);
  for (let i = A.length - 1; i >= db; i--) {
    const c = A[i]!;
    if (c === 0n) continue;
    if (c % lb !== 0n) return null;
    const t = c / lb;
    q[i - db] = t;
    for (let j = 0; j <= db; j++) {
      A[i - db + j] = A[i - db + j]! - t * B[j]!;
    }
  }
  if (zpNorm(A).length !== 0) return null;
  return zpNorm(q);
}

/** Content (gcd of the coefficients), with the sign of the leading coefficient. */
export function zpContent(a: ZPoly): bigint {
  const A = zpNorm(a);
  if (A.length === 0) return 0n;
  let g = 0n;
  for (const c of A) g = bgcd(g, c);
  return A[A.length - 1]! < 0n ? -g : g;
}

export function zpPrimitive(a: ZPoly): ZPoly {
  const c = zpContent(a);
  if (c === 0n || c === 1n) return zpNorm(a);
  return zpNorm(a).map((x) => x / c);
}

/**
 * Discriminant of the integer polynomial `f`, computed as
 * `(-1)^(n(n-1)/2) res(f, f') / lc(f)` with a fraction-free (Euclidean)
 * resultant.
 */
export function zpDiscriminant(f: ZPoly): bigint {
  const F = zpNorm(f);
  const n = F.length - 1;
  if (n < 1) throw new ValueError('discriminant undefined for constant polynomials');
  if (n === 1) return 1n;
  const res = zpResultant(F, zpDerivative(F));
  const sign = ((n * (n - 1)) / 2) % 2 === 0 ? 1n : -1n;
  const lc = F[n]!;
  const d = (sign * res) / lc;
  if (sign * res !== d * lc) {
    throw new ValueError('discriminant is not integral');
  }
  return d;
}

/**
 * Resultant of two integer polynomials, by the Euclidean recursion
 * `res(A, B) = (-1)^(deg A * deg B) lc(B)^(deg A - deg R) res(B, R)` with
 * `R = A mod B`.  Exact rational arithmetic; the degrees involved here are
 * tiny (a number field's defining polynomial).
 */
export function zpResultant(a: ZPoly, b: ZPoly): bigint {
  const trim = (v: Rational[]): Rational[] => {
    let l = v.length;
    while (l > 0 && v[l - 1]!.isZero()) l--;
    return v.slice(0, l);
  };
  let A = trim(zpNorm(a).map((c) => new Rational(c)));
  let B = trim(zpNorm(b).map((c) => new Rational(c)));
  if (A.length === 0 || B.length === 0) return 0n;
  let res = Rational.one();
  for (;;) {
    const dA = A.length - 1;
    const dB = B.length - 1;
    if (dB === 0) {
      // res(A, b0) = b0^deg(A)
      for (let i = 0; i < dA; i++) res = res.mul(B[0]!);
      break;
    }
    // R = A mod B
    const R = [...A];
    const lb = B[dB]!;
    for (let i = dA; i >= dB; i--) {
      const c = R[i]!;
      if (c.isZero()) continue;
      const t = c.div(lb);
      for (let j = 0; j <= dB; j++) {
        R[i - dB + j] = R[i - dB + j]!.sub(t.mul(B[j]!));
      }
    }
    const Rt = trim(R);
    if (Rt.length === 0) return 0n;
    const dR = Rt.length - 1;
    if ((dA * dB) % 2 === 1) res = res.neg();
    for (let i = 0; i < dA - dR; i++) res = res.mul(lb);
    A = B;
    B = Rt;
  }
  if (res.denominator !== 1n) {
    throw new ValueError('resultant is not integral');
  }
  return res.numerator;
}

// ---------------------------------------------------------------------------
// Polynomial arithmetic over F_p
// ---------------------------------------------------------------------------

export function fpNorm(f: ZPoly, p: bigint): ZPoly {
  const r = f.map((c) => mmod(c, p));
  let l = r.length;
  while (l > 0 && r[l - 1] === 0n) l--;
  return r.slice(0, l);
}

function fpAdd(a: ZPoly, b: ZPoly, p: bigint): ZPoly {
  return fpNorm(zpAdd(a, b), p);
}

function fpSub(a: ZPoly, b: ZPoly, p: bigint): ZPoly {
  return fpNorm(zpSub(a, b), p);
}

function fpMul(a: ZPoly, b: ZPoly, p: bigint): ZPoly {
  const A = fpNorm(a, p);
  const B = fpNorm(b, p);
  if (A.length === 0 || B.length === 0) return [];
  const r: ZPoly = new Array(A.length + B.length - 1).fill(0n);
  for (let i = 0; i < A.length; i++) {
    if (A[i] === 0n) continue;
    for (let j = 0; j < B.length; j++) {
      r[i + j] = (r[i + j]! + A[i]! * B[j]!) % p;
    }
  }
  return fpNorm(r, p);
}

function fpDivmod(a: ZPoly, b: ZPoly, p: bigint): [ZPoly, ZPoly] {
  const A = fpNorm(a, p);
  const B = fpNorm(b, p);
  if (B.length === 0) throw new ValueError('division by the zero polynomial');
  const db = B.length - 1;
  const inv = invMod(B[db]!, p);
  const R = [...A];
  const q: ZPoly = new Array(Math.max(0, A.length - db)).fill(0n);
  for (let i = R.length - 1; i >= db; i--) {
    const c = R[i]!;
    if (c === 0n) continue;
    const t = (c * inv) % p;
    q[i - db] = t;
    for (let j = 0; j <= db; j++) {
      R[i - db + j] = mmod(R[i - db + j]! - t * B[j]!, p);
    }
  }
  return [fpNorm(q, p), fpNorm(R, p)];
}

function fpRem(a: ZPoly, b: ZPoly, p: bigint): ZPoly {
  return fpDivmod(a, b, p)[1];
}

function fpMonic(a: ZPoly, p: bigint): ZPoly {
  const A = fpNorm(a, p);
  if (A.length === 0) return A;
  const inv = invMod(A[A.length - 1]!, p);
  return fpNorm(
    A.map((c) => c * inv),
    p
  );
}

export function fpGcd(a: ZPoly, b: ZPoly, p: bigint): ZPoly {
  let A = fpNorm(a, p);
  let B = fpNorm(b, p);
  while (B.length > 0) {
    const R = fpRem(A, B, p);
    A = B;
    B = R;
  }
  if (A.length === 0) return A;
  return fpMonic(A, p);
}

/** `base^exp mod (modulus, p)`. */
function fpPowMod(base: ZPoly, exp: bigint, modulus: ZPoly, p: bigint): ZPoly {
  let result: ZPoly = [1n];
  let b = fpRem(base, modulus, p);
  let e = exp;
  while (e > 0n) {
    if (e & 1n) result = fpRem(fpMul(result, b, p), modulus, p);
    b = fpRem(fpMul(b, b, p), modulus, p);
    e >>= 1n;
  }
  return fpNorm(result, p);
}

/** `f(x) = g(x^p)`; return `g` (the p-th root in `F_p[x]`). */
function fpPthRoot(f: ZPoly, p: bigint): ZPoly {
  const F = fpNorm(f, p);
  const out: ZPoly = [];
  const step = Number(p);
  for (let i = 0; i < F.length; i += step) {
    out.push(F[i]!);
  }
  for (let i = 0; i < F.length; i++) {
    if (i % step !== 0 && F[i] !== 0n) {
      throw new ValueError('polynomial is not a p-th power');
    }
  }
  return fpNorm(out, p);
}

/**
 * Squarefree decomposition of a monic polynomial over `F_p`
 * (Cohen, Algorithm 3.4.2).
 */
export function fpSquarefree(f: ZPoly, p: bigint): Array<[ZPoly, number]> {
  const out: Array<[ZPoly, number]> = [];
  let T: ZPoly = fpMonic(f, p);
  let e = 1;
  while (T.length - 1 > 0) {
    const Tp = fpNorm(zpDerivative(T), p);
    if (Tp.length === 0) {
      T = fpPthRoot(T, p);
      e *= Number(p);
      continue;
    }
    let T0 = fpGcd(T, Tp, p);
    let V = fpDivmod(T, T0, p)[0];
    let k = 0;
    while (V.length - 1 > 0) {
      k++;
      if (BigInt(k) % p === 0n) {
        T0 = fpDivmod(T0, V, p)[0];
        k++;
      }
      const W = fpGcd(T0, V, p);
      const A = fpDivmod(V, W, p)[0];
      if (A.length - 1 > 0) out.push([fpMonic(A, p), e * k]);
      V = W;
      T0 = fpDivmod(T0, W, p)[0];
    }
    T = T0;
  }
  return out;
}

/** Distinct-degree factorisation of a squarefree monic `u` over `F_p`. */
function fpDistinctDegree(u: ZPoly, p: bigint): Array<[ZPoly, number]> {
  const out: Array<[ZPoly, number]> = [];
  let v = fpMonic(u, p);
  let w: ZPoly = [0n, 1n]; // x
  let d = 0;
  while (v.length - 1 >= 2 * (d + 1)) {
    d++;
    w = fpPowMod(w, p, v, p);
    const g = fpGcd(fpSub(w, [0n, 1n], p), v, p);
    if (g.length - 1 > 0) {
      out.push([g, d]);
      v = fpDivmod(v, g, p)[0];
      w = fpRem(w, v, p);
    }
  }
  if (v.length - 1 > 0) out.push([v, v.length - 1]);
  return out;
}

/**
 * Equal-degree splitting (Cantor--Zassenhaus): split `u`, a product of
 * distinct monic irreducible factors all of degree `d`, over `F_p`.
 *
 * Candidate elements are enumerated deterministically so results are
 * reproducible.
 */
function fpEqualDegree(u: ZPoly, d: number, p: bigint): ZPoly[] {
  const n = u.length - 1;
  if (n === d) return [fpMonic(u, p)];
  const factors: ZPoly[] = [];
  const stack: ZPoly[] = [fpMonic(u, p)];
  // Deterministic supply of candidate polynomials of degree < n.
  const candidates: ZPoly[] = [];
  for (let deg = 1; deg < Math.max(2, n); deg++) {
    for (let c = 0n; c < (p < 40n ? p : 40n); c++) {
      const poly: ZPoly = new Array(deg + 1).fill(0n);
      poly[deg] = 1n;
      poly[0] = c;
      candidates.push(poly);
      if (deg >= 2) {
        const poly2: ZPoly = new Array(deg + 1).fill(0n);
        poly2[deg] = 1n;
        poly2[1] = 1n;
        poly2[0] = c;
        candidates.push(poly2);
      }
    }
  }
  while (stack.length > 0) {
    const v = stack.pop()!;
    if (v.length - 1 === d) {
      factors.push(v);
      continue;
    }
    let split = false;
    for (const a of candidates) {
      const ar = fpRem(a, v, p);
      if (ar.length === 0) continue;
      let g = fpGcd(ar, v, p);
      if (g.length - 1 > 0 && g.length - 1 < v.length - 1) {
        stack.push(g, fpDivmod(v, g, p)[0]);
        split = true;
        break;
      }
      let b: ZPoly;
      if (p === 2n) {
        // Trace map: a + a^2 + a^4 + ... + a^(2^(d-1))
        b = ar;
        let t = ar;
        for (let i = 1; i < d; i++) {
          t = fpRem(fpMul(t, t, p), v, p);
          b = fpAdd(b, t, p);
        }
      } else {
        const pd = p ** BigInt(d);
        b = fpSub(fpPowMod(ar, (pd - 1n) / 2n, v, p), [1n], p);
      }
      g = fpGcd(b, v, p);
      if (g.length - 1 > 0 && g.length - 1 < v.length - 1) {
        stack.push(g, fpDivmod(v, g, p)[0]);
        split = true;
        break;
      }
    }
    if (!split) {
      throw new ValueError('equal-degree factorisation failed to split');
    }
  }
  return factors;
}

/**
 * Factor a polynomial over `F_p` into monic irreducibles with multiplicities.
 * The leading coefficient is dropped (the result is the factorisation of the
 * monic normalisation).
 */
export function fpFactor(f: ZPoly, p: bigint): Array<[ZPoly, number]> {
  const F = fpMonic(f, p);
  if (F.length - 1 <= 0) return [];
  const out: Array<[ZPoly, number]> = [];
  for (const [sqf, mult] of fpSquarefree(F, p)) {
    for (const [part, d] of fpDistinctDegree(sqf, p)) {
      for (const irr of fpEqualDegree(part, d, p)) {
        out.push([irr, mult]);
      }
    }
  }
  out.sort((a, b) => {
    if (a[0].length !== b[0].length) return a[0].length - b[0].length;
    for (let i = a[0].length - 1; i >= 0; i--) {
      if (a[0][i] !== b[0][i]) return a[0][i]! < b[0][i]! ? -1 : 1;
    }
    return 0;
  });
  return out;
}

/** Roots of `f` in `F_p`, each listed once. */
export function fpRoots(f: ZPoly, p: bigint): bigint[] {
  const F = fpNorm(f, p);
  if (F.length === 0) throw new ValueError('the zero polynomial has every root');
  const roots: bigint[] = [];
  // gcd(f, x^p - x) isolates the linear factors.
  const xp = fpPowMod([0n, 1n], p, F, p);
  const g = fpGcd(fpSub(xp, [0n, 1n], p), F, p);
  if (g.length - 1 <= 0) return roots;
  for (const [irr] of fpFactor(g, p)) {
    if (irr.length - 1 === 1) {
      roots.push(mmod(-irr[0]!, p));
    }
  }
  roots.sort((a, b) => (a < b ? -1 : a > b ? 1 : 0));
  return roots;
}

// ---------------------------------------------------------------------------
// Factorisation over Z (Zassenhaus with a single big prime)
// ---------------------------------------------------------------------------

/**
 * Mignotte's bound `2^n * ||f||_2` on the absolute value of any coefficient of
 * a factor of `f` (rounded up to an integer).
 */
function mignotteBound(f: ZPoly): bigint {
  const F = zpNorm(f);
  const n = F.length - 1;
  let sum = 0n;
  for (const c of F) sum += c * c;
  const norm = isqrt(sum) + 1n;
  return (1n << BigInt(n)) * norm;
}

/**
 * Factor a squarefree primitive integer polynomial into irreducible primitive
 * factors (Zassenhaus).  A single prime `p` larger than `2 * B * |lc|` is used,
 * which makes the Hensel lift unnecessary.
 */
export function zpFactorSquarefree(f: ZPoly): ZPoly[] {
  const F = zpPrimitive(f);
  const n = F.length - 1;
  if (n <= 0) return [];
  if (n === 1) return [F];
  const lc = F[n]!;
  const need = 2n * mignotteBound(F) * babs(lc) + 1n;
  if (need > 1n << 200n) {
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: integer polynomial factorisation for this size requires a Hensel lift'
    );
  }
  let p = need;
  if (p < 3n) p = 3n;
  for (;;) {
    p += 1n;
    if (!is_prime(p)) continue;
    if (lc % p === 0n) continue;
    const fp = fpNorm(F, p);
    const g = fpGcd(fp, fpNorm(zpDerivative(F), p), p);
    if (g.length - 1 !== 0) continue;
    break;
  }
  const modFactors = fpFactor(F, p).map(([h]) => h);
  const r = modFactors.length;
  if (r === 1) return [F];

  const out: ZPoly[] = [];
  let remaining = F;
  const used = new Set<number>();
  let size = 1;
  while (2 * size <= r - used.size) {
    const remLc = remaining[remaining.length - 1]!;
    const indices = [...Array(r).keys()].filter((i) => !used.has(i));
    let progressed = false;
    for (const combo of combinations(indices, size)) {
      let prod: ZPoly = [mmod(remLc, p)];
      for (const i of combo) prod = fpMul(prod, modFactors[i]!, p);
      const cand = zpPrimitive(prod.map((c) => balanced(c, p)));
      if (cand.length - 1 <= 0) continue;
      const q = zpExactDiv(remaining, cand);
      if (q !== null) {
        out.push(cand);
        remaining = zpPrimitive(q);
        for (const i of combo) used.add(i);
        progressed = true;
        break;
      }
    }
    if (!progressed) size++;
  }
  if (remaining.length - 1 > 0) out.push(remaining);
  return out;
}

function combinations(items: number[], k: number): number[][] {
  const out: number[][] = [];
  const cur: number[] = [];
  const rec = (start: number) => {
    if (cur.length === k) {
      out.push([...cur]);
      return;
    }
    for (let i = start; i < items.length; i++) {
      cur.push(items[i]!);
      rec(i + 1);
      cur.pop();
    }
  };
  rec(0);
  return out;
}

/**
 * Test whether the integer polynomial `f` (degree >= 1) is irreducible over `Q`.
 *
 * Mirrors PARI's `polisirreducible`: a polynomial with a repeated factor or
 * more than one Zassenhaus factor is reducible.  A non-primitive constant
 * factor does not matter over `Q`.
 */
export function zpIsIrreducibleOverQ(f: ZPoly): boolean {
  const F = zpPrimitive(f);
  const n = F.length - 1;
  if (n < 1) return false;
  if (n === 1) return true;
  const g = zpGcdPoly(F, zpDerivative(F));
  if (g.length - 1 > 0) return false; // not squarefree
  return zpFactorSquarefree(F).length === 1;
}

/** GCD of two integer polynomials, as a primitive integer polynomial. */
export function zpGcdPoly(a: ZPoly, b: ZPoly): ZPoly {
  let A = zpPrimitive(a);
  let B = zpPrimitive(b);
  while (B.length > 0) {
    // pseudo-remainder
    const db = B.length - 1;
    const lb = B[db]!;
    let R = [...A];
    let guard = 0;
    while (R.length - 1 >= db && zpNorm(R).length > 0) {
      const dr = R.length - 1;
      const c = R[dr]!;
      R = R.map((x) => x * lb);
      for (let j = 0; j <= db; j++) {
        R[dr - db + j] = R[dr - db + j]! - c * B[j]!;
      }
      R = zpNorm(R);
      if (guard++ > 10000) throw new ValueError('polynomial gcd did not terminate');
    }
    A = B;
    B = R.length === 0 ? [] : zpPrimitive(R);
  }
  return A.length === 0 ? [] : zpPrimitive(A);
}

// ---------------------------------------------------------------------------
// Integer matrices: HNF, kernels, inverses
// ---------------------------------------------------------------------------

/**
 * Row Hermite normal form of an integer matrix with `n` columns, assuming the
 * rows span a rank-`n` lattice.  Row `i` of the result has its pivot in column
 * `i` (upper triangular), positive diagonal, and entries above each pivot
 * reduced modulo it.
 */
export function hnf(rows: bigint[][], n: number): bigint[][] {
  const m = rows.map((r) => {
    const c = new Array<bigint>(n).fill(0n);
    for (let i = 0; i < n; i++) c[i] = r[i] ?? 0n;
    return c;
  });
  let r = 0;
  for (let c = 0; c < n; c++) {
    let pivot = -1;
    for (let i = r; i < m.length; i++) {
      if (m[i]![c] !== 0n) {
        pivot = i;
        break;
      }
    }
    if (pivot === -1) continue;
    [m[r], m[pivot]] = [m[pivot]!, m[r]!];
    for (let i = r + 1; i < m.length; i++) {
      if (m[i]![c] === 0n) continue;
      const [g, u, v] = xgcd(m[r]![c]!, m[i]![c]!);
      const a = m[r]![c]! / g;
      const b = m[i]![c]! / g;
      const rowR = m[r]!;
      const rowI = m[i]!;
      const newR = new Array<bigint>(n);
      const newI = new Array<bigint>(n);
      for (let j = 0; j < n; j++) {
        newR[j] = u * rowR[j]! + v * rowI[j]!;
        newI[j] = -b * rowR[j]! + a * rowI[j]!;
      }
      m[r] = newR;
      m[i] = newI;
    }
    if (m[r]![c]! < 0n) {
      m[r] = m[r]!.map((x) => -x);
    }
    const d = m[r]![c]!;
    for (let i = 0; i < r; i++) {
      let q = m[i]![c]! / d;
      if (m[i]![c]! - q * d < 0n) q -= 1n;
      if (q !== 0n) {
        for (let j = 0; j < n; j++) m[i]![j] = m[i]![j]! - q * m[r]![j]!;
      }
    }
    r++;
    if (r === n) break;
  }
  if (r < n) {
    throw new ValueError('HNF input does not have full rank');
  }
  return m.slice(0, n);
}

/**
 * "Lower" Hermite normal form: row `i` has its last nonzero entry in column
 * `i`.  This is the shape PARI uses for `nfbasis` and for ideal HNFs, so that
 * the first basis vector is a rational integer multiple of `1`.
 */
export function hnfLower(rows: bigint[][], n: number): bigint[][] {
  const rev = (v: bigint[]): bigint[] => {
    const out = new Array<bigint>(n);
    for (let i = 0; i < n; i++) out[i] = v[n - 1 - i] ?? 0n;
    return out;
  };
  const H = hnf(rows.map(rev), n);
  // Row i of H has pivot at reversed column i, i.e. original column n-1-i.
  const out = H.map(rev);
  out.reverse();
  return out;
}

/** Determinant of a square integer matrix (via the rational Gauss algorithm). */
export function intDet(m: bigint[][]): bigint {
  const n = m.length;
  const a = m.map((r) => r.map((x) => new Rational(x)));
  let det = Rational.one();
  let sign = 1n;
  for (let c = 0; c < n; c++) {
    let pivot = -1;
    for (let i = c; i < n; i++) {
      if (!a[i]![c]!.isZero()) {
        pivot = i;
        break;
      }
    }
    if (pivot === -1) return 0n;
    if (pivot !== c) {
      [a[c], a[pivot]] = [a[pivot]!, a[c]!];
      sign = -sign;
    }
    const pv = a[c]![c]!;
    det = det.mul(pv);
    for (let i = c + 1; i < n; i++) {
      const f = a[i]![c]!.div(pv);
      if (f.isZero()) continue;
      for (let j = c; j < n; j++) {
        a[i]![j] = a[i]![j]!.sub(f.mul(a[c]![j]!));
      }
    }
  }
  const out = det.mul(new Rational(sign));
  if (out.denominator !== 1n) throw new ValueError('determinant is not integral');
  return out.numerator;
}

/** Inverse of a square rational matrix. */
export function ratInverse(m: Rational[][]): Rational[][] {
  const n = m.length;
  const a = m.map((r, i) => [
    ...r,
    ...Array.from({ length: n }, (_, j) => (i === j ? Rational.one() : Rational.zero())),
  ]);
  for (let c = 0; c < n; c++) {
    let pivot = -1;
    for (let i = c; i < n; i++) {
      if (!a[i]![c]!.isZero()) {
        pivot = i;
        break;
      }
    }
    if (pivot === -1) throw new ValueError('matrix is singular');
    [a[c], a[pivot]] = [a[pivot]!, a[c]!];
    const pv = a[c]![c]!.inv();
    for (let j = 0; j < 2 * n; j++) a[c]![j] = a[c]![j]!.mul(pv);
    for (let i = 0; i < n; i++) {
      if (i === c) continue;
      const f = a[i]![c]!;
      if (f.isZero()) continue;
      for (let j = 0; j < 2 * n; j++) {
        a[i]![j] = a[i]![j]!.sub(f.mul(a[c]![j]!));
      }
    }
  }
  return a.map((r) => r.slice(n));
}

/** Basis of `{x in F_p^n : A x = 0}` for an `m x n` matrix `A`. */
export function nullspaceModP(A: bigint[][], p: bigint): bigint[][] {
  const m = A.length;
  const n = m === 0 ? 0 : A[0]!.length;
  const a = A.map((r) => r.map((x) => mmod(x, p)));
  const pivotCol: number[] = [];
  let row = 0;
  for (let c = 0; c < n && row < m; c++) {
    let pivot = -1;
    for (let i = row; i < m; i++) {
      if (a[i]![c] !== 0n) {
        pivot = i;
        break;
      }
    }
    if (pivot === -1) continue;
    [a[row], a[pivot]] = [a[pivot]!, a[row]!];
    const inv = invMod(a[row]![c]!, p);
    for (let j = 0; j < n; j++) a[row]![j] = (a[row]![j]! * inv) % p;
    for (let i = 0; i < m; i++) {
      if (i === row) continue;
      const f = a[i]![c]!;
      if (f === 0n) continue;
      for (let j = 0; j < n; j++) {
        a[i]![j] = mmod(a[i]![j]! - f * a[row]![j]!, p);
      }
    }
    pivotCol.push(c);
    row++;
  }
  const isPivot = new Array<boolean>(n).fill(false);
  for (const c of pivotCol) isPivot[c] = true;
  const basis: bigint[][] = [];
  for (let free = 0; free < n; free++) {
    if (isPivot[free]) continue;
    const v = new Array<bigint>(n).fill(0n);
    v[free] = 1n;
    for (let i = 0; i < pivotCol.length; i++) {
      v[pivotCol[i]!] = mmod(-a[i]![free]!, p);
    }
    basis.push(v);
  }
  return basis;
}

/** Basis of `{c in F_p^m : c A = 0}` for an `m x n` matrix `A`. */
function leftNullspaceModP(A: bigint[][], p: bigint): bigint[][] {
  const m = A.length;
  const n = m === 0 ? 0 : A[0]!.length;
  const T: bigint[][] = Array.from({ length: n }, (_, j) =>
    Array.from({ length: m }, (_, i) => A[i]![j]!)
  );
  return nullspaceModP(T, p);
}

// ---------------------------------------------------------------------------
// nfbasis / nfdisc (Round 2)
// ---------------------------------------------------------------------------

/** The result of `nfbasis`. */
export interface NfBasisResult {
  /** Monic integral polynomial defining the field, in the variable `theta`. */
  g: ZPoly;
  /**
   * `n x n` integer matrix; the i-th element of the integral basis is
   * `(1/den) * sum_j basis[i][j] * theta^j`.  Lower triangular (PARI's shape),
   * so `basis[0]` is a multiple of `1`.
   */
  basis: bigint[][];
  /** Common denominator of the integral basis. */
  den: bigint;
  /** `[O_K : Z[theta]]`. */
  index: bigint;
  /** `disc(O_K)`. */
  disc: bigint;
}

/**
 * The multiplication table of an order, together with the data needed to move
 * between the order basis and the power basis of `theta`.
 */
class OrderTable {
  readonly n: number;
  readonly g: ZPoly;
  readonly M: bigint[][];
  readonly den: bigint;
  private readonly Minv: Rational[][];
  private readonly sc: bigint[][][];

  constructor(g: ZPoly, M: bigint[][], den: bigint) {
    this.n = zpDeg(g);
    this.g = g;
    this.M = M;
    this.den = den;
    this.Minv = ratInverse(M.map((r) => r.map((x) => new Rational(x))));
    const n = this.n;
    this.sc = [];
    for (let i = 0; i < n; i++) {
      const rowI: bigint[][] = [];
      for (let j = 0; j < n; j++) {
        // (b_i b_j) * den^2 in the power basis
        const prod = zpMulMod(M[i]!, M[j]!, g);
        // coordinates in the order basis: (prod / den) * Minv
        const w: Rational[] = [];
        for (let k = 0; k < n; k++) {
          w.push(new Rational(prod[k] ?? 0n, den));
        }
        const t: bigint[] = [];
        for (let k = 0; k < n; k++) {
          let acc = Rational.zero();
          for (let l = 0; l < n; l++) {
            acc = acc.add(w[l]!.mul(this.Minv[l]![k]!));
          }
          if (acc.denominator !== 1n) {
            throw new ValueError('order basis is not closed under multiplication');
          }
          t.push(acc.numerator);
        }
        rowI.push(t);
      }
      this.sc.push(rowI);
    }
  }

  /** Structure constants: coordinates of `b_i * b_j` in the order basis. */
  structure(i: number, j: number): bigint[] {
    return this.sc[i]![j]!;
  }

  /** Product of two coordinate vectors modulo `p`. */
  mulModP(u: bigint[], v: bigint[], p: bigint): bigint[] {
    const n = this.n;
    const out = new Array<bigint>(n).fill(0n);
    for (let i = 0; i < n; i++) {
      if (u[i] === 0n) continue;
      for (let j = 0; j < n; j++) {
        if (v[j] === 0n) continue;
        const c = (u[i]! * v[j]!) % p;
        const s = this.sc[i]![j]!;
        for (let k = 0; k < n; k++) {
          out[k] = (out[k]! + c * s[k]!) % p;
        }
      }
    }
    return out.map((x) => mmod(x, p));
  }

  /**
   * Coordinates in the order basis of a rational vector given in the power
   * basis of `theta`.
   */
  fromPowerBasis(v: Rational[]): Rational[] {
    const n = this.n;
    // v = t * (M / den)  =>  t = (v * den) * M^{-1}
    const scaled = v.map((c) => c.mul(new Rational(this.den)));
    const out: Rational[] = [];
    for (let k = 0; k < n; k++) {
      let acc = Rational.zero();
      for (let l = 0; l < n; l++) {
        acc = acc.add(scaled[l]!.mul(this.Minv[l]![k]!));
      }
      out.push(acc);
    }
    return out;
  }
}

/**
 * The `p`-radical of the order `O`: `{x in O : x^(p^k) in pO}` for
 * `p^k >= n`.  Returned as the `n x n` matrix of a Z-basis in terms of the
 * order basis.
 */
function pRadical(T: OrderTable, p: bigint): bigint[][] {
  const n = T.n;
  let k = 0;
  let q = 1n;
  while (q < BigInt(n)) {
    q *= p;
    k++;
  }
  if (k === 0) k = 1;
  // Matrix of the Frobenius x -> x^p in the order basis: row i = coords(b_i^p).
  const one = oneCoords(T);
  const R: bigint[][] = [];
  for (let i = 0; i < n; i++) {
    const e = new Array<bigint>(n).fill(0n);
    e[i] = 1n;
    // b_i^p by square and multiply
    let result = one;
    let base = e;
    let exp = p;
    while (exp > 0n) {
      if (exp & 1n) result = T.mulModP(result, base, p);
      base = T.mulModP(base, base, p);
      exp >>= 1n;
    }
    R.push(result);
  }
  // Compose k times: the map is c |-> c R, so we need R^k.
  let Rk = R;
  for (let i = 1; i < k; i++) {
    Rk = matMulModP(Rk, R, p);
  }
  const kernel = leftNullspaceModP(Rk, p);
  const rows: bigint[][] = [];
  for (let i = 0; i < n; i++) {
    const e = new Array<bigint>(n).fill(0n);
    e[i] = p;
    rows.push(e);
  }
  for (const v of kernel) rows.push(v);
  return hnf(rows, n);
}

function matMulModP(A: bigint[][], B: bigint[][], p: bigint): bigint[][] {
  const m = A.length;
  const n = B[0]!.length;
  const k = B.length;
  const out: bigint[][] = [];
  for (let i = 0; i < m; i++) {
    const row = new Array<bigint>(n).fill(0n);
    for (let l = 0; l < k; l++) {
      const a = A[i]![l]!;
      if (a === 0n) continue;
      for (let j = 0; j < n; j++) {
        row[j] = (row[j]! + a * B[l]![j]!) % p;
      }
    }
    out.push(row.map((x) => mmod(x, p)));
  }
  return out;
}

/** Coordinates of `1` in the order basis. */
function oneCoords(T: OrderTable): bigint[] {
  const n = T.n;
  // 1 has power-basis coordinates (1, 0, ..., 0); in the order basis this is
  // (den, 0, ..., 0) * M^{-1}.
  const v: Rational[] = new Array(n).fill(Rational.zero());
  v[0] = new Rational(T.den);
  const inv = ratInverse(T.M.map((r) => r.map((x) => new Rational(x))));
  const out: bigint[] = [];
  for (let k = 0; k < n; k++) {
    let acc = Rational.zero();
    for (let l = 0; l < n; l++) acc = acc.add(v[l]!.mul(inv[l]![k]!));
    if (acc.denominator !== 1n) throw new ValueError('1 is not in the order');
    out.push(acc.numerator);
  }
  return out;
}

/**
 * One Round-2 step: given an order `O` (as an `OrderTable`) and its `p`-radical
 * `I`, return the ring of multipliers `{x in K : x I subset I}` as a new
 * `(matrix, denominator)` pair in the power basis, or `null` when `O` is
 * already `p`-maximal.
 */
function multiplierRing(
  T: OrderTable,
  A: bigint[][],
  p: bigint
): { M: bigint[][]; den: bigint } | null {
  const n = T.n;
  const Ainv = ratInverse(A.map((r) => r.map((x) => new Rational(x))));
  // Row i of L collects, for c = e_i, the vectors (y u_j) expressed in the
  // I-basis and reduced mod p.
  const L: bigint[][] = [];
  for (let i = 0; i < n; i++) {
    const row: bigint[] = [];
    for (let j = 0; j < n; j++) {
      // w = sum_k A[j][k] * (b_i b_k)
      const w = new Array<bigint>(n).fill(0n);
      for (let k = 0; k < n; k++) {
        const a = A[j]![k]!;
        if (a === 0n) continue;
        const s = T.structure(i, k);
        for (let l = 0; l < n; l++) w[l] = w[l]! + a * s[l]!;
      }
      // t = w * Ainv, integral because I is an ideal of O
      for (let l = 0; l < n; l++) {
        let acc = Rational.zero();
        for (let m = 0; m < n; m++) {
          acc = acc.add(new Rational(w[m]!).mul(Ainv[m]![l]!));
        }
        if (acc.denominator !== 1n) {
          throw new ValueError('the p-radical is not an ideal of the order');
        }
        row.push(mmod(acc.numerator, p));
      }
    }
    L.push(row);
  }
  const V = leftNullspaceModP(L, p);
  if (V.length === 0) return null;
  // New module: O + (1/p) * span(V), in the power basis with denominator den*p.
  const rows: bigint[][] = [];
  for (let i = 0; i < n; i++) {
    rows.push(T.M[i]!.map((x) => x * p));
  }
  let added = false;
  for (const v of V) {
    const w = new Array<bigint>(n).fill(0n);
    for (let i = 0; i < n; i++) {
      if (v[i] === 0n) continue;
      for (let j = 0; j < n; j++) w[j] = w[j]! + v[i]! * T.M[i]![j]!;
    }
    rows.push(w);
    added = true;
  }
  if (!added) return null;
  const H = hnfLower(rows, n);
  let den = T.den * p;
  // Reduce by the common content.
  let g = den;
  for (const r of H) for (const x of r) g = bgcd(g, x);
  if (g > 1n) {
    den /= g;
    for (const r of H) for (let j = 0; j < n; j++) r[j] = r[j]! / g;
  }
  return { M: H, den };
}

/**
 * Compute a `p`-maximal order containing `Z[theta]`.
 *
 * Port of Cohen, Algorithm 6.1.8 ("Round 2"), which is the algorithm behind
 * PARI's `maxord`.
 */
function pMaximalOrder(
  g: ZPoly,
  start: { M: bigint[][]; den: bigint },
  p: bigint
): { M: bigint[][]; den: bigint } {
  let cur = start;
  for (let guard = 0; guard < 200; guard++) {
    const T = new OrderTable(g, cur.M, cur.den);
    const radical = pRadical(T, p);
    const next = multiplierRing(T, radical, p);
    if (next === null) return cur;
    cur = next;
  }
  throw new ValueError('Round 2 did not terminate');
}

/**
 * `nfbasis(g)` for a monic integral polynomial `g`: an integral basis of the
 * maximal order of `Q[x]/(g)`, its index in `Z[theta]` and the field
 * discriminant.
 *
 * @see Reference: reference/pari/src/basemath/base2.c:nfbasis
 */
export function nfbasis(g: ZPoly): NfBasisResult {
  const G = zpNorm(g);
  const n = G.length - 1;
  if (n < 1) throw new ValueError('nfbasis requires a polynomial of degree at least 1');
  if (G[n] !== 1n) throw new ValueError('nfbasis requires a monic polynomial');
  const polyDisc = zpDiscriminant(G);
  if (polyDisc === 0n) throw new ValueError('polynomial is not squarefree');

  let cur: { M: bigint[][]; den: bigint } = {
    M: Array.from({ length: n }, (_, i) =>
      Array.from({ length: n }, (_, j): bigint => (i === j ? 1n : 0n))
    ),
    den: 1n,
  };

  if (n > 1) {
    for (const [p, e] of intFactor(babs(polyDisc))) {
      if (e < 2n) continue;
      cur = pMaximalOrder(G, cur, p);
    }
  }

  const det = intDet(cur.M);
  const index = cur.den ** BigInt(n) / babs(det);
  if (index * babs(det) !== cur.den ** BigInt(n)) {
    throw new ValueError('inconsistent index computation');
  }
  const disc = polyDisc / (index * index);
  if (disc * index * index !== polyDisc) {
    throw new ValueError('inconsistent discriminant computation');
  }
  return { g: G, basis: cur.M, den: cur.den, index, disc };
}

/**
 * `nfdisc(g)`: the discriminant of the maximal order of `Q[x]/(g)`.
 *
 * @see Reference: reference/pari/src/basemath/base2.c:nfdisc
 */
export function nfdisc(g: ZPoly): bigint {
  return nfbasis(g).disc;
}

// ---------------------------------------------------------------------------
// nfgaloisconj: conjugates of theta lying in Q[x]/(g)
// ---------------------------------------------------------------------------

/**
 * Rational reconstruction: given `a mod m`, find `x/y = a (mod m)` with
 * `|x| <= N` and `0 < y <= D`, or `null` if no such pair exists.
 */
export function rationalReconstruct(a: bigint, m: bigint, N: bigint, D: bigint): Rational | null {
  let r0 = m;
  let r1 = mmod(a, m);
  let s0 = 0n;
  let s1 = 1n;
  while (r1 > N) {
    if (r1 === 0n) break;
    const q = r0 / r1;
    [r0, r1] = [r1, r0 - q * r1];
    [s0, s1] = [s1, s0 - q * s1];
  }
  if (r1 > N) return null;
  const y = babs(s1);
  if (y === 0n || y > D) return null;
  if (bgcd(r1, y) !== 1n) return null;
  const x = s1 < 0n ? -r1 : r1;
  return new Rational(x, y);
}

/**
 * `numberofconjugates(T)`: a rigorous *upper bound* on `#Aut(K)` for
 * `K = Q[x]/(T)`, obtained from the factorisation shapes of `T` modulo several
 * unramified primes.
 *
 * If `m = #Aut(K)` and `p` is unramified, `Aut(K)` acts on the primes above
 * `p`; the stabiliser of `P` is the (cyclic) decomposition group, of order
 * dividing `f_P`, so the primes of residue degree `d` fall into orbits of size
 * `m / |D_P|` with `|D_P| | d`.  Summing `d` over them shows `m | L[d] * d`,
 * where `L[d]` is the number of degree-`d` factors of `T mod p`.  `m` also
 * divides `n`.  The gcd of all these is therefore a multiple of `m`.
 *
 * @see Reference: reference/pari/src/basemath/galconj.c:3113 (numberofconjugates)
 */
export function numberofconjugates(T: ZPoly, pinit = 2n): bigint {
  const n = zpDeg(T);
  if (n === 1) return 1n;
  const nbmax = n < 10 ? 20 : 2 * n + 1;
  const disc = zpDiscriminant(T);
  let nbtest = 0;
  let c = BigInt(n);
  for (let p = pinit; ; p++) {
    if (!is_prime(p)) continue;
    if (disc % p === 0n) continue; // unramified / squarefree mod p
    nbtest++;
    const L = new Array<number>(n + 1).fill(0);
    let nb = 0;
    for (const [gi, e] of fpFactor(T, p)) {
      L[zpDeg(gi)] = (L[zpDeg(gi)] ?? 0) + e;
      nb += e;
    }
    if (L[Math.floor(n / nb)] === nb) {
      // all factors have the same degree: no information, probably Galois
      if (c === BigInt(n) && nbtest > 10) break;
    } else {
      c = bgcd(c, BigInt(L[1]!));
      for (let i = 2; i <= n; i++) {
        if (L[i]) c = bgcd(c, BigInt(L[i]! * i));
      }
      if (c === 1n) break;
    }
    if (nbtest === nbmax) break;
  }
  return c;
}

/**
 * A rigorous bound on `|D * c_i|`, where `beta = sum c_i theta^i` is any root
 * of `g` lying in `K = Q[x]/(g)` and `D = [O_K : Z[theta]]`.
 *
 * `beta` is a root of the monic `g`, hence an algebraic integer whose every
 * archimedean conjugate is at most the Cauchy root bound `B = 1 + max|g_i|`.
 * Writing `c = V^{-1} w` with `V` the Vandermonde matrix of the roots and
 * `w` the vector of conjugates of `beta`, `|det V| = sqrt(|disc g|)` and
 * Hadamard bounds every cofactor by `(n-1)^((n-1)/2) B^((n-1)(n-2))`, so
 *
 *   `|c_i| <= n * ceil((n-1)^((n-1)/2)) * B^((n-1)(n-2)+1) / sqrt(|disc g|)`.
 *
 * `D * beta` lies in `Z[theta]` because `D * O_K subseteq Z[theta]`, so the
 * quantity returned bounds the *integer* coordinates the lattice search below
 * has to find.
 */
function conjugateCoeffBound(g: ZPoly, D: bigint): bigint {
  const n = zpDeg(g);
  let maxc = 0n;
  for (let i = 0; i < n; i++) {
    const c = babs(g[i] ?? 0n);
    if (c > maxc) maxc = c;
  }
  const B = maxc + 1n;
  const hadamard = isqrt(BigInt(n - 1) ** BigInt(n - 1)) + 1n;
  const num = D * BigInt(n) * hadamard * B ** BigInt((n - 1) * (n - 2) + 1);
  const den = isqrt(babs(zpDiscriminant(g))); // <= |det V|
  return num / den + 1n;
}

/**
 * All leading principal minors `d_1, ..., d_k` of a symmetric positive
 * definite integer matrix, by fraction-free (Bareiss) elimination.  Every
 * division is exact.
 */
function leadingPrincipalMinors(gram: bigint[][]): bigint[] {
  const k = gram.length;
  const a = gram.map((r) => [...r]);
  const d: bigint[] = [];
  let prev = 1n;
  for (let i = 0; i < k; i++) {
    d.push(a[i]![i]!);
    if (i === k - 1) break;
    for (let j = i + 1; j < k; j++) {
      for (let l = i + 1; l < k; l++) {
        a[j]![l] = (a[j]![l]! * a[i]![i]! - a[j]![i]! * a[i]![l]!) / prev;
      }
    }
    prev = a[i]![i]!;
  }
  return d;
}

/**
 * `nfgaloisconj(g)`: the roots of `g` that lie in `K = Q[x]/(g)`, returned as
 * coefficient vectors in the power basis of `theta`.  The identity `theta` is
 * always the first entry.  There is no degree restriction.
 *
 * Method (PARI's, `galconj.c`): pick a prime `p` unramified in `K` for which
 * `g` splits completely, and lift its roots `r_1, ..., r_n` to `Z/p^k`.  An
 * automorphism `sigma` is determined by `sigma(theta) = beta` with
 * `beta(r_1) = r_j` for a single `j`, because the `n` embeddings
 * `K -> Q_p, theta |-> r_i` are permuted simply transitively by `Aut(K)`.  The
 * coefficients of `D * beta` (`D = [O_K : Z[theta]]`) are therefore the short
 * vector of the lattice
 *
 *   `{ (a, e) in Z^(n+1) : sum a_i r_1^i = e * D * r_j  (mod p^k) }`
 *
 * with `e = 1`, and LLL finds it.  This replaces the earlier enumeration of
 * all `n!` permutations of the `p`-adic roots, which forced a degree-8 cap.
 *
 * Two independent guarantees make the answer *proved*, not heuristic:
 *
 * - no false positives: every candidate is checked exactly with `g(beta) = 0`
 *   in `Q[x]/(g)`;
 * - no false negatives: for every `j` that produced nothing, the reduced basis
 *   certifies that the lattice has *no* vector at all of norm `<= A*sqrt(n+1)`
 *   (via `min_i ||b*_i|| <= lambda_1`), so no such `beta` exists.  When the
 *   certificate is inconclusive the `p`-adic precision is squared and the
 *   search repeated.  `numberofconjugates` provides a cheap early exit.
 *
 * @see Reference: reference/pari/src/basemath/galconj.c:2988 (galoisconj4_main)
 * @see Reference: reference/pari/src/basemath/galconj.c:3113 (numberofconjugates)
 */
export function nfgaloisconj(g: ZPoly): Rational[][] {
  const G = zpNorm(g);
  const n = G.length - 1;
  if (n === 1) return [[new Rational(-(G[0] ?? 0n))]];
  const identity: Rational[] = new Array(n).fill(Rational.zero());
  identity[1] = Rational.one();
  if (n === 2) {
    // The conjugate of theta is -b - theta for g = x^2 + b x + c.
    const other: Rational[] = new Array(n).fill(Rational.zero());
    other[0] = new Rational(-(G[1] ?? 0n));
    other[1] = new Rational(-1n);
    return [identity, other];
  }

  const disc = zpDiscriminant(G);
  if (disc === 0n) throw new ValueError('nfgaloisconj requires a squarefree polynomial');
  const D = nfbasis(G).index;
  const upperBound = numberofconjugates(G);

  // Find a small prime for which g splits completely into distinct roots.
  let p = 0n;
  let roots0: bigint[] = [];
  for (let cand = 2n; cand < 100000n; cand++) {
    if (!is_prime(cand)) continue;
    if (disc % cand === 0n) continue;
    const r = fpRoots(G, cand);
    if (r.length === n) {
      p = cand;
      roots0 = r;
      break;
    }
  }
  if (p === 0n) {
    // By Chebotarev a totally split prime has positive density, so this only
    // happens if the first 100000 integers were unlucky.
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: nfgaloisconj could not find a totally split prime'
    );
  }

  const A = conjugateCoeffBound(G, D);
  const V2 = A * A * BigInt(n + 1); // squared norm bound of the sought vector
  const gp = zpDerivative(G);

  const liftRoots = (pk: bigint): bigint[] =>
    roots0.map((r0) => {
      let r = r0;
      let mod = p;
      while (mod < pk) {
        mod = mod * mod;
        if (mod > pk) mod = pk;
        const gv = zpEvalMod(G, r, mod);
        const gpv = zpEvalMod(gp, r, mod);
        r = mmod(r - gv * invMod(gpv, mod), mod);
      }
      return mmod(r, pk);
    });

  let target = A * A + 1n;
  for (let attempt = 0; attempt < MAX_GALOISCONJ_ATTEMPTS; attempt++) {
    let pk = p;
    while (pk <= target) pk *= p;
    const roots = liftRoots(pk);
    const r1 = roots[0]!;
    const pows: bigint[] = [1n];
    for (let i = 1; i < n; i++) pows.push(mmod(pows[i - 1]! * r1, pk));

    const found: Rational[][] = [];
    const seen = new Set<string>();
    const certified = new Array<boolean>(n).fill(false);

    for (let j = 0; j < n; j++) {
      // Lattice basis: the solutions of  sum a_i r_1^i = e * D * r_j (mod p^k),
      // with the extra coordinate e scaled by A so that the wanted vector
      // (a, A) is short.
      const rows: bigint[][] = [];
      const row0 = new Array<bigint>(n + 1).fill(0n);
      row0[0] = pk;
      rows.push(row0);
      for (let i = 1; i < n; i++) {
        const r = new Array<bigint>(n + 1).fill(0n);
        r[0] = -pows[i]!;
        r[i] = 1n;
        rows.push(r);
      }
      const lastRow = new Array<bigint>(n + 1).fill(0n);
      lastRow[0] = mmod(D * roots[j]!, pk);
      lastRow[n] = A;
      rows.push(lastRow);

      const reduced = LLL(new IntegerMatrix(n + 1, n + 1, rows)) as IntegerMatrix;
      const R: bigint[][] = [];
      for (let i = 0; i <= n; i++) {
        const row: bigint[] = [];
        for (let k = 0; k <= n; k++) row.push(reduced.get(i, k).value);
        R.push(row);
      }

      for (const row of R) {
        const e = row[n]!;
        if (e === 0n || e % A !== 0n) continue;
        const m = e / A;
        const a: bigint[] = [];
        let ok = true;
        for (let k = 0; k < n; k++) {
          const v = row[k]!;
          if (v % m !== 0n) {
            ok = false;
            break;
          }
          a.push(v / m);
        }
        if (!ok) continue;
        const coeffs = a.map((x) => new Rational(x, D));
        if (!isRootOfPoly(coeffs, G)) continue;
        const key = coeffs.map((r) => `${r.numerator}/${r.denominator}`).join(',');
        if (seen.has(key)) continue;
        seen.add(key);
        found.push(coeffs);
      }
      if (BigInt(found.length) >= upperBound) break;

      // Rigorous non-existence certificate for this j: every nonzero lattice
      // vector has norm >= min_i ||b*_i||, and ||b*_i||^2 = d_i / d_{i-1}.
      const gram: bigint[][] = [];
      for (let i = 0; i <= n; i++) {
        const grow: bigint[] = [];
        for (let k = 0; k <= n; k++) {
          let s = 0n;
          for (let l = 0; l <= n; l++) s += R[i]![l]! * R[k]![l]!;
          grow.push(s);
        }
        gram.push(grow);
      }
      let minOK = true;
      let prev = 1n;
      for (const di of leadingPrincipalMinors(gram)) {
        if (di <= V2 * prev) {
          minOK = false;
          break;
        }
        prev = di;
      }
      certified[j] = minOK;
    }

    let complete = BigInt(found.length) >= upperBound;
    if (!complete) {
      // Which r_j are realised by a conjugate we did find?
      const realised = new Array<boolean>(n).fill(false);
      for (const c of found) {
        let val = 0n;
        let pw = 1n;
        for (let i = 0; i < n; i++) {
          const ci = mmod(c[i]!.numerator * invMod(mmod(c[i]!.denominator, pk), pk), pk);
          val = mmod(val + ci * pw, pk);
          pw = mmod(pw * r1, pk);
        }
        const j = roots.indexOf(val);
        if (j >= 0) realised[j] = true;
      }
      complete = realised.every((r, j) => r || certified[j]!);
    }

    if (complete) {
      // Put the identity first, as PARI does.
      const idKey = identity.map((r) => `${r.numerator}/${r.denominator}`).join(',');
      found.sort((x, y) => {
        const kx = x.map((r) => `${r.numerator}/${r.denominator}`).join(',');
        const ky = y.map((r) => `${r.numerator}/${r.denominator}`).join(',');
        if (kx === idKey) return -1;
        if (ky === idKey) return 1;
        return kx < ky ? -1 : kx > ky ? 1 : 0;
      });
      return found;
    }
    target = target * target;
  }

  throw new NotImplementedError(
    `SAGE_NOT_IMPLEMENTED: nfgaloisconj could not certify the conjugates of a degree ${n} ` +
      `polynomial after ${MAX_GALOISCONJ_ATTEMPTS} precision doublings`
  );
}

/** Number of times `nfgaloisconj` squares the p-adic precision before giving up. */
const MAX_GALOISCONJ_ATTEMPTS = 12;

function zpEvalMod(f: ZPoly, x: bigint, m: bigint): bigint {
  let acc = 0n;
  for (let i = f.length - 1; i >= 0; i--) {
    acc = mmod(acc * x + f[i]!, m);
  }
  return acc;
}

/** Check that `sum c_i theta^i` is a root of `g` in `Q[x]/(g)`. */
function isRootOfPoly(c: Rational[], g: ZPoly): boolean {
  const n = zpDeg(g);
  // Work with a common denominator so that everything stays integral.
  let den = 1n;
  for (const r of c) den = (den / bgcd(den, r.denominator)) * r.denominator;
  const num = c.map((r) => r.numerator * (den / r.denominator));
  // Evaluate g at beta = num/den, i.e. compute den^n * g(num/den) mod g.
  let acc: ZPoly = [];
  let power: ZPoly = [1n];
  let denPow = 1n;
  const terms: Array<{ coeff: bigint; poly: ZPoly; denPow: bigint }> = [];
  for (let i = 0; i <= n; i++) {
    terms.push({ coeff: g[i] ?? 0n, poly: power, denPow });
    power = zpMulMod(power, num, g);
    denPow *= den;
  }
  const total = den ** BigInt(n);
  for (const t of terms) {
    if (t.coeff === 0n) continue;
    const scale = total / t.denPow;
    acc = zpAdd(acc, zpScale(t.poly, t.coeff * scale));
  }
  return zpNorm(acc).length === 0;
}

// ---------------------------------------------------------------------------
// idealprimedec: the Buchmann--Lenstra ("round 4") prime decomposition
// ---------------------------------------------------------------------------

/**
 * The multiplication table of `O_K` in a fixed `Z`-basis `w_1, ..., w_n`:
 * `mul[i][j]` are the coordinates of `w_i * w_j` in that basis.  `w_1` must
 * be `1`.
 */
export type MulTable = bigint[][][];

/** One prime of `O_K` above `p`, as produced by `primedec`. */
export interface PrimeDecEntry {
  /**
   * An `F_p`-basis (rows, coordinates in the order basis) of `P / p O_K`
   * inside `O_K / p O_K`, so that `P = p O_K + sum Z * gens[i]`.
   */
  gens: bigint[][];
  /** Ramification index `e = v_P(p)`. */
  e: bigint;
  /** Residue class degree `f = [O_K/P : F_p]`. */
  f: bigint;
}

/** Product of two elements of `O_K/p`, given in coordinates. */
function algMul(a: bigint[], b: bigint[], mul: MulTable, p: bigint): bigint[] {
  const n = a.length;
  const out = new Array<bigint>(n).fill(0n);
  for (let i = 0; i < n; i++) {
    if (a[i] === 0n) continue;
    for (let j = 0; j < n; j++) {
      if (b[j] === 0n) continue;
      const c = mmod(a[i]! * b[j]!, p);
      const row = mul[i]![j]!;
      for (let k = 0; k < n; k++) out[k] = mmod(out[k]! + c * row[k]!, p);
    }
  }
  return out;
}

/** `a^e` in `O_K/p`. */
function algPow(a: bigint[], e: bigint, mul: MulTable, p: bigint): bigint[] {
  const n = a.length;
  let result = new Array<bigint>(n).fill(0n);
  result[0] = 1n % p; // w_1 = 1
  let base = a;
  let k = e;
  while (k > 0n) {
    if (k & 1n) result = algMul(result, base, mul, p);
    base = algMul(base, base, mul, p);
    k >>= 1n;
  }
  return result;
}

/**
 * The matrix of multiplication by `a` acting on row vectors: row `i` is
 * `w_i * a`, so `v |-> v * M` is `x |-> x * a`.
 */
function algMultable(a: bigint[], mul: MulTable, p: bigint): bigint[][] {
  const n = a.length;
  const M: bigint[][] = [];
  for (let i = 0; i < n; i++) {
    const ei = new Array<bigint>(n).fill(0n);
    ei[i] = 1n;
    M.push(algMul(ei, a, mul, p));
  }
  return M;
}

/** Reduced row echelon basis of the row space of `rows` over `F_p`. */
function fpRowSpace(rows: bigint[][], n: number, p: bigint): bigint[][] {
  const a = rows.map((r) => {
    const c = new Array<bigint>(n).fill(0n);
    for (let i = 0; i < n; i++) c[i] = mmod(r[i] ?? 0n, p);
    return c;
  });
  let row = 0;
  for (let c = 0; c < n && row < a.length; c++) {
    let pivot = -1;
    for (let i = row; i < a.length; i++) {
      if (a[i]![c] !== 0n) {
        pivot = i;
        break;
      }
    }
    if (pivot === -1) continue;
    [a[row], a[pivot]] = [a[pivot]!, a[row]!];
    const inv = invMod(a[row]![c]!, p);
    for (let j = 0; j < n; j++) a[row]![j] = mmod(a[row]![j]! * inv, p);
    for (let i = 0; i < a.length; i++) {
      if (i === row) continue;
      const f = a[i]![c]!;
      if (f === 0n) continue;
      for (let j = 0; j < n; j++) a[i]![j] = mmod(a[i]![j]! - f * a[row]![j]!, p);
    }
    row++;
  }
  return a.slice(0, row);
}

/** Transpose. */
function fpTranspose(M: bigint[][]): bigint[][] {
  if (M.length === 0) return [];
  const rows = M.length;
  const cols = M[0]!.length;
  const T: bigint[][] = [];
  for (let j = 0; j < cols; j++) {
    const r: bigint[] = [];
    for (let i = 0; i < rows; i++) r.push(M[i]![j]!);
    T.push(r);
  }
  return T;
}

/** `{ v : v * M = 0 }` over `F_p`, as a list of row vectors. */
function fpLeftKernel(M: bigint[][], p: bigint): bigint[][] {
  if (M.length === 0) return [];
  return nullspaceModP(fpTranspose(M), p);
}

/** Matrix product over `F_p`. */
function fpMatMul(A: bigint[][], B: bigint[][], p: bigint): bigint[][] {
  const n = A.length;
  if (n === 0) return [];
  const m = B[0]!.length;
  const k = B.length;
  const C: bigint[][] = [];
  for (let i = 0; i < n; i++) {
    const row = new Array<bigint>(m).fill(0n);
    for (let t = 0; t < k; t++) {
      const a = A[i]![t]!;
      if (a === 0n) continue;
      for (let j = 0; j < m; j++) row[j] = mmod(row[j]! + a * B[t]![j]!, p);
    }
    C.push(row);
  }
  return C;
}

/** Inverse of a square matrix over `F_p`. */
function fpInverse(M: bigint[][], p: bigint): bigint[][] {
  const n = M.length;
  const a = M.map((r, i) => {
    const row = r.map((x) => mmod(x, p));
    const ext = new Array<bigint>(n).fill(0n);
    ext[i] = 1n;
    return row.concat(ext);
  });
  for (let c = 0; c < n; c++) {
    let pivot = -1;
    for (let i = c; i < n; i++) {
      if (a[i]![c] !== 0n) {
        pivot = i;
        break;
      }
    }
    if (pivot === -1) throw new ValueError('matrix is singular mod p');
    [a[c], a[pivot]] = [a[pivot]!, a[c]!];
    const inv = invMod(a[c]![c]!, p);
    for (let j = 0; j < 2 * n; j++) a[c]![j] = mmod(a[c]![j]! * inv, p);
    for (let i = 0; i < n; i++) {
      if (i === c) continue;
      const f = a[i]![c]!;
      if (f === 0n) continue;
      for (let j = 0; j < 2 * n; j++) a[i]![j] = mmod(a[i]![j]! - f * a[c]![j]!, p);
    }
  }
  return a.map((r) => r.slice(n));
}

/**
 * The minimal polynomial over `F_p` of the element whose multiplication matrix
 * (acting on row vectors) is `M` in an algebra whose first basis vector is 1.
 *
 * @see Reference: reference/pari/src/basemath/base2.c:2185 (pol_min)
 */
function algMinPoly(M: bigint[][], p: bigint): ZPoly {
  const d = M.length;
  const powers: bigint[][] = [];
  let v = new Array<bigint>(d).fill(0n);
  v[0] = 1n % p; // the identity
  for (let i = 0; i <= d; i++) {
    powers.push(v);
    v = fpMatMul([v], M, p)[0]!;
  }
  // First linear dependency among 1, a, a^2, ...
  for (let k = 1; k <= d; k++) {
    const ker = fpLeftKernel(powers.slice(0, k + 1), p);
    if (ker.length > 0) {
      const c = ker[0]!;
      // normalise: leading coefficient 1
      let lead = k;
      while (lead >= 0 && c[lead] === 0n) lead--;
      const inv = invMod(c[lead]!, p);
      return zpNorm(c.slice(0, lead + 1).map((x) => mmod(x * inv, p)));
    }
  }
  throw new ValueError('no minimal polynomial found');
}

/**
 * `primedec(mul, p)`: the primes of `O_K` above `p`, by the Buchmann--Lenstra
 * "round 4" algorithm, which does *not* need `Z[theta]` to be `p`-maximal and
 * therefore works at inessential discriminant divisors (where no generator
 * `gamma` with `p` prime to `[O_K : Z[gamma]]` exists, so the Dedekind--Kummer
 * theorem cannot be applied at all -- e.g. `p = 2` in
 * `Q[x]/(x^3 - x^2 - 2x - 8)`).
 *
 * The algorithm:
 *
 * 1. `I_p`, the `p`-radical of `O_K`, is the kernel of the `k`-th power of the
 *    `F_p`-linear Frobenius `x |-> x^p` on `O_K/p`, for `p^k >= n`.
 * 2. `A = (O_K/p) / I_p` is an etale (separable commutative) `F_p`-algebra,
 *    hence a product of finite fields; its maximal ideals are the `P/pO_K`.
 * 3. `A` is split by picking `a` in `ker(x |-> x^p - x)` that is not a scalar:
 *    its minimal polynomial divides `x^p - x`, so it splits into distinct
 *    linear factors, and the images of `a - lambda_i` cut `A` into pieces.
 *    Repeat until every piece is a field.
 * 4. `f = n - dim(P/pO_K)`; `e = v_P(p)` is obtained from the largest `k` with
 *    `p O_K` contained in `P^k`, using exact lattice arithmetic.
 *
 * @param mul - multiplication table of `O_K` in a basis whose first element is 1
 * @param p - a rational prime
 * @see Reference: reference/pari/src/basemath/base2.c:2248 (primedec_aux)
 * @see Reference: reference/pari/src/basemath/base2.c:2150 (pradical)
 */
export function primedec(mul: MulTable, p: bigint): PrimeDecEntry[] {
  const n = mul.length;
  if (n === 0) throw new ValueError('primedec: empty multiplication table');

  // 1. p-radical: kernel of Frobenius^k with p^k >= n.
  // Frobenius is F_p-linear because (sum v_i w_i)^p = sum v_i w_i^p mod p.
  const frob: bigint[][] = [];
  for (let i = 0; i < n; i++) {
    const ei = new Array<bigint>(n).fill(0n);
    ei[i] = 1n;
    frob.push(algPow(ei, p, mul, p));
  }
  let m = frob;
  let q = p;
  while (q < BigInt(n)) {
    q *= p;
    m = fpMatMul(m, frob, p);
  }
  const Ip = fpLeftKernel(m, p);
  // phi = x -> x^p - x
  const phi = frob.map((row, i) => row.map((x, j) => mmod(x - (i === j ? 1n : 0n), p)));

  // 2/3. Split the etale algebra (O_K/p)/I_p.
  const maximal: bigint[][][] = [];
  const worklist: bigint[][][] = [Ip];
  let guard = 0;
  while (worklist.length > 0) {
    if (++guard > 4 * n * n + 16) {
      throw new ValueError('primedec: splitting did not terminate');
    }
    const H = worklist.pop()!;
    const r = H.length;
    if (r === n) continue; // the unit ideal
    // Complete H to a basis of F_p^n whose (r+1)-st vector is 1 = w_1.
    const Mrows: bigint[][] = H.map((row) => [...row]);
    const cand: bigint[][] = [];
    const e1 = new Array<bigint>(n).fill(0n);
    e1[0] = 1n;
    cand.push(e1);
    for (let i = 0; i < n; i++) {
      const ei = new Array<bigint>(n).fill(0n);
      ei[i] = 1n;
      cand.push(ei);
    }
    for (const c of cand) {
      if (Mrows.length === n) break;
      if (fpRowSpace([...Mrows, c], n, p).length > Mrows.length) Mrows.push(c);
    }
    if (Mrows.length !== n) throw new ValueError('primedec: could not complete a basis');
    const Minv = fpInverse(Mrows, p);
    const M2 = Mrows.slice(r); // basis of the chosen complement of H
    const project = (X: bigint[][]): bigint[][] =>
      fpMatMul(fpMatMul(M2, X, p), Minv, p).map((row) => row.slice(r));

    const phi2 = project(phi);
    const kernel = fpLeftKernel(phi2, p);
    const dim = kernel.length; // A2 is a product of `dim` fields

    if (dim <= 1) {
      // A2 is a field: H is a maximal ideal of O_K/p.
      maximal.push(H);
      continue;
    }

    // Split A2 with an element of ker(x -> x^p - x) that is not a scalar.
    let didSplit = false;
    for (const u of kernel) {
      const a = fpMatMul([u], M2, p)[0]!;
      const mula = algMultable(a, mul, p);
      const mu = algMinPoly(project(mula), p);
      // a^p = a in A2, so mu divides x^p - x: it splits into distinct linear
      // factors and has deg(mu) roots.  deg(mu) = 1 means a is a scalar.
      if (zpDeg(mu) <= 1) continue;
      const roots = fpRoots(mu, p);
      if (roots.length !== zpDeg(mu)) {
        throw new ValueError('primedec: the minimal polynomial does not split');
      }
      const pieces = roots.map((lambda) =>
        fpRowSpace(
          [...H, ...mula.map((row, i) => row.map((x, j) => mmod(x - (i === j ? lambda : 0n), p)))],
          n,
          p
        )
      );
      // n roots == dim components means each piece is already maximal.
      for (const piece of pieces) (roots.length === dim ? maximal : worklist).push(piece);
      didSplit = true;
      break;
    }
    if (!didSplit) {
      throw new ValueError('primedec: no splitting element found in an etale algebra');
    }
  }

  // 4. Ramification indices, by exact lattice arithmetic in O_K.
  const out: PrimeDecEntry[] = [];
  for (const H of maximal) {
    const f = BigInt(n - H.length);
    const lattice = primeLattice(H, p, n);
    let e = 1n;
    let power = lattice;
    for (;;) {
      const next = latticeMul(power, lattice, mul);
      if (!latticeContainsP(next, p, n)) break;
      power = next;
      e += 1n;
      if (e > BigInt(n)) throw new ValueError('primedec: ramification index out of range');
    }
    out.push({ gens: H, e, f });
  }
  return out;
}

/** The `Z`-lattice `p O_K + span(H)` as an `n x n` HNF matrix. */
function primeLattice(H: bigint[][], p: bigint, n: number): bigint[][] {
  const rows: bigint[][] = [];
  for (const h of H) rows.push([...h]);
  for (let i = 0; i < n; i++) {
    const r = new Array<bigint>(n).fill(0n);
    r[i] = p;
    rows.push(r);
  }
  return hnf(rows, n);
}

/** Product of two full-rank `O_K`-lattices given by `n x n` bases. */
function latticeMul(A: bigint[][], B: bigint[][], mul: MulTable): bigint[][] {
  const n = A.length;
  const rows: bigint[][] = [];
  for (const a of A) {
    for (const b of B) {
      const out = new Array<bigint>(n).fill(0n);
      for (let i = 0; i < n; i++) {
        if (a[i] === 0n) continue;
        for (let j = 0; j < n; j++) {
          if (b[j] === 0n) continue;
          const c = a[i]! * b[j]!;
          const row = mul[i]![j]!;
          for (let k = 0; k < n; k++) out[k] = out[k]! + c * row[k]!;
        }
      }
      rows.push(out);
    }
  }
  return hnf(rows, n);
}

/** Is `p O_K` contained in the lattice `L` (given as an upper-triangular HNF)? */
function latticeContainsP(L: bigint[][], p: bigint, n: number): boolean {
  for (let i = 0; i < n; i++) {
    const v = new Array<bigint>(n).fill(0n);
    v[i] = p;
    if (!latticeContains(L, v, n)) return false;
  }
  return true;
}

/** Membership test in an upper-triangular full-rank HNF lattice. */
function latticeContains(L: bigint[][], v: bigint[], n: number): boolean {
  const x = [...v];
  for (let i = 0; i < n; i++) {
    // find the pivot column of row i
    let c = i;
    while (c < n && L[i]![c] === 0n) c++;
    if (c === n) return false;
    if (x[c]! % L[i]![c]! !== 0n) return false;
    const q = x[c]! / L[i]![c]!;
    if (q !== 0n) for (let j = 0; j < n; j++) x[j] = x[j]! - q * L[i]![j]!;
  }
  return x.every((y) => y === 0n);
}

// ---------------------------------------------------------------------------
// quadunit / quadunitnorm: fundamental unit of a real quadratic order
// ---------------------------------------------------------------------------

/**
 * `quadunit(D)`: the fundamental unit of the quadratic order of discriminant
 * `D > 0`, returned as the pair `[u, v]` with `epsilon = u + v * w_D`, where
 * `w_D` is PARI's `quadgen(D)`, i.e. the root of `quadpoly(D)`:
 *
 * - `D = 0 (mod 4)`: `quadpoly = x^2 - D/4`,     `w_D = sqrt(D)/2`;
 * - `D = 1 (mod 4)`: `quadpoly = x^2 - x - (D-1)/4`, `w_D = (1 + sqrt(D))/2`.
 *
 * This is a transcription of PARI's `quadunit_uv_basecase`: the continued
 * fraction expansion of `(P + sqrt(D))/Q` is run until the period closes (the
 * even-period test `p1 == p` and the odd-period test `q == q1`), while the
 * convergent numerators/denominators `u_i`, `v_i` are accumulated; the
 * fundamental solution of Pell's equation is then read off from the last two
 * convergents.  The returned unit is the smallest one `> 1` and has norm `+-1`.
 *
 * PARI switches to `quadunit_uv` (quad.c:429), which multiplies the same
 * elementary matrices with a product tree, once `D >= 2000000`; that is a
 * pure speed optimisation and returns the identical `[u, v]`.
 *
 * @param D - a positive discriminant, `D = 0, 1 (mod 4)`, not a perfect square
 * @see Reference: reference/pari/src/basemath/quad.c:281 (quadunit_uv_basecase)
 * @see Reference: reference/pari/src/basemath/quad.c:476 (quadunit)
 */
export function quadunit(D: bigint): [bigint, bigint] {
  if (D <= 0n) {
    throw new ValueError(`quadunit: disc <= 0: ${D}`);
  }
  if (mmod(D, 4n) > 1n) {
    throw new ValueError(`quadunit: disc % 4 > 1: ${D}`);
  }
  const d = isqrt(D); // floor(sqrt(D))
  if (d * d === D) {
    throw new ValueError(`quadunit: issquare(disc) = 1: ${D}`);
  }
  const rem = D - d * d; // sqrtremi(D, &a)
  const m = (D & 1n) === 1n; // mpodd(D)

  let p = d;
  let q1 = rem >> 1n;
  let q = 2n;
  if (((d & 1n) === 1n) !== m) {
    p = d - 1n;
    q1 = q1 + d; // q1 = (D - p^2)/2
  }
  let u1 = 2n;
  let u2 = p;
  let v1 = 0n;
  let v2 = 1n;
  let first = true;
  let a: bigint;
  let b: bigint;
  let c: bigint;
  for (;;) {
    const t = q;
    if (first) {
      first = false;
      q = q1;
    } else {
      const A = (p + d) / q;
      const r = (p + d) % q;
      const p1 = p;
      p = d - r;
      if (p1 === p) {
        // even period
        a = u2 * u2;
        b = v2 * v2;
        c = (u2 + v2) * (u2 + v2);
        break;
      }
      const nu = u1 + A * u2;
      u1 = u2;
      u2 = nu;
      const nv = v1 + A * v2;
      v1 = v2;
      v2 = nv;
      q = q1 - A * (p - p1);
    }
    q1 = t;
    if (q === t) {
      // odd period
      a = u1 * u2;
      b = v1 * v2;
      c = (u1 + v1) * (u2 + v2);
      break;
    }
  }
  let u = (a + D * b) / q;
  const v = (c - (a + b)) / q;
  if (m) u = u - v;
  return [u >> 1n, v];
}

/**
 * `quadunitnorm(D)`: the norm (`+1` or `-1`) of the fundamental unit of the
 * quadratic order of discriminant `D > 0`.
 *
 * Computed from `quadunit`'s own output rather than PARI's `quadunit_q`
 * shortcut; the two agree because `N(u + v w_D) = u^2 + u v - v^2 (D-1)/4`
 * (`D` odd) resp. `u^2 - v^2 D/4` (`D` even).
 *
 * @see Reference: reference/pari/src/basemath/quad.c:520 (quadunitnorm)
 */
export function quadunitnorm(D: bigint): bigint {
  const [u, v] = quadunit(D);
  return quadNormUV(D, u, v);
}

/** `N(u + v w_D)` where `w_D = quadgen(D)`. */
export function quadNormUV(D: bigint, u: bigint, v: bigint): bigint {
  if ((D & 1n) === 1n) {
    // w^2 = w + (D-1)/4, so N(u + v w) = u^2 + u v - v^2 (D-1)/4
    return u * u + u * v - v * v * ((D - 1n) / 4n);
  }
  // w^2 = D/4, so N(u + v w) = u^2 - v^2 D/4
  return u * u - v * v * (D / 4n);
}

// ---------------------------------------------------------------------------
// Non-monic / non-integral defining polynomials
// ---------------------------------------------------------------------------

/**
 * Given a monic polynomial `f` over `Q` of degree `n`, return the smallest
 * positive integer `d` such that `g(y) = d^n f(y/d)` has integer coefficients,
 * together with `g`.  Then `theta = d * alpha` is an algebraic integer
 * generating the same field.
 *
 * This mirrors what Sage does with `K.pari_polynomial()`.
 */
export function integralDefiningPolynomial(coeffs: Rational[]): { g: ZPoly; scale: bigint } {
  const n = coeffs.length - 1;
  if (n < 1) throw new ValueError('defining polynomial must have degree at least 1');
  if (!coeffs[n]!.eq(Rational.one())) {
    throw new ValueError('integralDefiningPolynomial requires a monic polynomial');
  }
  // d must satisfy: d^(n-i) * a_i is an integer for every i.
  let d = 1n;
  for (let i = 0; i < n; i++) {
    const den = coeffs[i]!.denominator;
    if (den === 1n) continue;
    // Grow d until d^(n-i) is divisible by den.
    for (const [p, e] of intFactor(den)) {
      const k = BigInt(n - i);
      // need v_p(d) * k >= e
      let needed = e / k;
      if (needed * k < e) needed += 1n;
      let have = 0n;
      let t = d;
      while (t % p === 0n) {
        t /= p;
        have += 1n;
      }
      if (have < needed) {
        d *= p ** (needed - have);
      }
    }
  }
  const g: ZPoly = [];
  for (let i = 0; i <= n; i++) {
    const v = coeffs[i]!.mul(new Rational(d ** BigInt(n - i)));
    if (v.denominator !== 1n) {
      throw new ValueError('failed to make the defining polynomial integral');
    }
    g.push(v.numerator);
  }
  return { g: zpNorm(g), scale: d };
}
