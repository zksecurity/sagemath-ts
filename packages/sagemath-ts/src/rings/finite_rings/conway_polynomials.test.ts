/**
 * Tests for sage/rings/finite_rings/conway_polynomials
 *
 * The whole point of the table is that every entry is monic, **irreducible**
 * and **primitive** (SageMath: `sage/rings/finite_rings/conway_polynomials.py`,
 * whose database is Frank Luebeck's `CPimport.txt`).  Nothing outside this file
 * checks that, so the table is verified here from first principles with
 * self-contained arithmetic:
 *
 *  - irreducibility via Rabin's test (`x^(p^n) == x mod f` and
 *    `gcd(x^(p^(n/q)) - x, f) == 1` for every prime `q | n`), the test
 *    NTL/FLINT use (`reference/flint/src/nmod_poly/is_irreducible_rabin.c`);
 *  - primitivity by checking that `x^((p^n-1)/r) != 1 mod f` for every prime
 *    `r | p^n - 1`.
 *
 * The factorisations needed for the primitivity test are computed with a
 * local Pollard-rho, so the test does not depend on `arith.factor` (which
 * currently gives up on some 19-digit composites such as 2^62 - 1).
 */

import { describe, test, expect } from 'bun:test';
import {
  CONWAY_POLYNOMIALS,
  conway_polynomial,
  has_conway_polynomial,
  available_degrees,
  available_characteristics,
} from './conway_polynomials.js';

// ---------------------------------------------------------------------------
// Integer helpers (local, so the table is verified independently of arith/)
// ---------------------------------------------------------------------------

function powmod(base: bigint, exp: bigint, m: bigint): bigint {
  let result = 1n;
  let b = base % m;
  let e = exp;
  while (e > 0n) {
    if (e & 1n) result = (result * b) % m;
    b = (b * b) % m;
    e >>= 1n;
  }
  return result;
}

/** Deterministic Miller-Rabin for n < 3.3 * 10^24. */
function isPrime(n: bigint): boolean {
  if (n < 2n) return false;
  for (const p of [2n, 3n, 5n, 7n, 11n, 13n, 17n, 19n, 23n, 29n, 31n, 37n]) {
    if (n === p) return true;
    if (n % p === 0n) return false;
  }
  let d = n - 1n;
  let s = 0n;
  while (d % 2n === 0n) {
    d /= 2n;
    s += 1n;
  }
  for (const a of [2n, 3n, 5n, 7n, 11n, 13n, 17n, 19n, 23n, 29n, 31n, 37n]) {
    let x = powmod(a, d, n);
    if (x === 1n || x === n - 1n) continue;
    let composite = true;
    for (let i = 1n; i < s; i++) {
      x = (x * x) % n;
      if (x === n - 1n) {
        composite = false;
        break;
      }
    }
    if (composite) return false;
  }
  return true;
}

function gcdBig(a: bigint, b: bigint): bigint {
  let x = a < 0n ? -a : a;
  let y = b < 0n ? -b : b;
  while (y) {
    [x, y] = [y, x % y];
  }
  return x;
}

/** Brent's variant of Pollard's rho. */
function pollardRho(n: bigint): bigint {
  if (n % 2n === 0n) return 2n;
  for (let c = 1n; ; c++) {
    let x = 2n;
    let y = 2n;
    let d = 1n;
    while (d === 1n) {
      x = (x * x + c) % n;
      y = (y * y + c) % n;
      y = (y * y + c) % n;
      d = gcdBig(x > y ? x - y : y - x, n);
    }
    if (d !== n) return d;
  }
}

/** Distinct prime divisors of n (n > 0). */
function primeDivisors(n: bigint): bigint[] {
  const out = new Set<bigint>();
  const stack: bigint[] = [n];
  while (stack.length > 0) {
    const m = stack.pop()!;
    if (m === 1n) continue;
    if (isPrime(m)) {
      out.add(m);
      continue;
    }
    // strip small factors first, rho does not split prime powers well
    let rest = m;
    let split = false;
    for (let p = 2n; p < 100000n; p += p === 2n ? 1n : 2n) {
      if (p * p > rest) break;
      if (rest % p === 0n) {
        out.add(p);
        while (rest % p === 0n) rest /= p;
        split = true;
      }
    }
    if (rest === 1n) continue;
    if (isPrime(rest)) {
      out.add(rest);
      continue;
    }
    if (split) {
      stack.push(rest);
      continue;
    }
    const d = pollardRho(rest);
    stack.push(d, rest / d);
  }
  return [...out].sort((a, b) => (a < b ? -1 : a > b ? 1 : 0));
}

// ---------------------------------------------------------------------------
// Polynomial arithmetic in GF(p)[x] / (f), f monic of degree n
// ---------------------------------------------------------------------------

type Poly = number[]; // increasing degree order, length <= n, no trailing zeros

function trim(a: Poly): Poly {
  let i = a.length;
  while (i > 0 && a[i - 1] === 0) i--;
  return a.slice(0, i);
}

/** Full coefficient list of the Conway polynomial: [c_0, ..., c_{n-1}, 1]. */
function monicCoeffs(p: number, n: number): Poly {
  return [...conway_polynomial(p, n), 1];
}

function polyRem(a: Poly, f: Poly, p: number): Poly {
  const r = [...a];
  const df = f.length - 1;
  for (let i = r.length - 1; i >= df; i--) {
    const c = r[i]!;
    if (c === 0) continue;
    for (let j = 0; j <= df; j++) {
      r[i - df + j] = (((r[i - df + j]! - c * f[j]!) % p) + p) % p;
    }
  }
  return trim(r);
}

function polyMulMod(a: Poly, b: Poly, f: Poly, p: number): Poly {
  if (a.length === 0 || b.length === 0) return [];
  const prod = new Array<number>(a.length + b.length - 1).fill(0);
  for (let i = 0; i < a.length; i++) {
    const ai = a[i]!;
    if (ai === 0) continue;
    for (let j = 0; j < b.length; j++) {
      prod[i + j] = (prod[i + j]! + ai * b[j]!) % p;
    }
  }
  return polyRem(prod, f, p);
}

function polyPowMod(a: Poly, e: bigint, f: Poly, p: number): Poly {
  let result: Poly = [1];
  let base = polyRem([...a], f, p);
  let exp = e;
  while (exp > 0n) {
    if (exp & 1n) result = polyMulMod(result, base, f, p);
    base = polyMulMod(base, base, f, p);
    exp >>= 1n;
  }
  return result;
}

function polySub(a: Poly, b: Poly, p: number): Poly {
  const len = Math.max(a.length, b.length);
  const out = new Array<number>(len).fill(0);
  for (let i = 0; i < len; i++) {
    out[i] = ((((a[i] ?? 0) - (b[i] ?? 0)) % p) + p) % p;
  }
  return trim(out);
}

function polyGcd(a: Poly, b: Poly, p: number): Poly {
  let x = trim([...a]);
  let y = trim([...b]);
  while (y.length > 0) {
    // make y monic so polyRem applies
    const inv = modInverse(y[y.length - 1]!, p);
    const ym = y.map((c) => (c * inv) % p);
    const r = polyRem(x, ym, p);
    x = ym;
    y = r;
  }
  return x;
}

function modInverse(a: number, p: number): number {
  let r = ((a % p) + p) % p;
  let result = 1;
  // p is prime and small: a^(p-2)
  let e = p - 2;
  while (e > 0) {
    if (e & 1) result = (result * r) % p;
    r = (r * r) % p;
    e >>= 1;
  }
  return result;
}

/** Rabin's irreducibility test for monic f of degree n over GF(p). */
function isIrreducible(f: Poly, p: number): boolean {
  const n = f.length - 1;
  if (n <= 0) return false;
  const P = BigInt(p);
  const x: Poly = [0, 1];
  // x^(p^n) == x
  if (polyPowMod(x, P ** BigInt(n), f, p).join(',') !== polyRem(x, f, p).join(',')) {
    return false;
  }
  for (const q of primeDivisors(BigInt(n))) {
    const m = BigInt(n) / q;
    const g = polyGcd(polySub(polyPowMod(x, P ** m, f, p), x, p), f, p);
    if (g.length !== 1) return false;
  }
  return true;
}

/** Does the class of x generate GF(p^n)^*? */
function isPrimitive(f: Poly, p: number): boolean {
  const n = f.length - 1;
  const order = BigInt(p) ** BigInt(n) - 1n;
  const x: Poly = [0, 1];
  if (polyPowMod(x, order, f, p).join(',') !== '1') return false;
  for (const r of primeDivisors(order)) {
    if (polyPowMod(x, order / r, f, p).join(',') === '1') return false;
  }
  return true;
}

// ---------------------------------------------------------------------------

function allEntries(): { p: number; n: number }[] {
  const out: { p: number; n: number }[] = [];
  for (const p of available_characteristics()) {
    for (const n of available_degrees(p)) out.push({ p, n });
  }
  return out;
}

describe('Conway polynomial database', () => {
  test('local helpers agree with known factorisations', () => {
    expect(primeDivisors(2n ** 62n - 1n)).toEqual([3n, 715827883n, 2147483647n]);
    expect(primeDivisors(2n ** 64n - 1n)).toEqual([3n, 5n, 17n, 257n, 641n, 65537n, 6700417n]);
    expect(isPrime(2n ** 61n - 1n)).toBe(true);
  });

  test('every entry has exactly n coefficients, all reduced mod p', () => {
    for (const { p, n } of allEntries()) {
      const c = conway_polynomial(p, n);
      expect(c.length).toBe(n);
      for (const ci of c) {
        expect(Number.isInteger(ci)).toBe(true);
        expect(ci).toBeGreaterThanOrEqual(0);
        expect(ci).toBeLessThan(p);
      }
      // A Conway polynomial always has a nonzero constant term (it is
      // irreducible of degree >= 2, so x does not divide it).
      expect(c[0]).not.toBe(0);
    }
  });

  test('every entry is irreducible over GF(p)', () => {
    for (const { p, n } of allEntries()) {
      const f = monicCoeffs(p, n);
      expect({ p, n, irreducible: isIrreducible(f, p) }).toEqual({ p, n, irreducible: true });
    }
  });

  test('every entry is primitive (x generates GF(p^n)^*)', () => {
    for (const { p, n } of allEntries()) {
      const f = monicCoeffs(p, n);
      expect({ p, n, primitive: isPrimitive(f, p) }).toEqual({ p, n, primitive: true });
    }
  });

  test('Conway compatibility: C_{p,m}(x^((p^n-1)/(p^m-1))) == 0 mod C_{p,n}', () => {
    // This is the defining property of the Conway family (Luebeck): for every
    // m | n the norm-compatible embedding GF(p^m) -> GF(p^n) sends the root of
    // C_{p,m} to x^((p^n-1)/(p^m-1)).  A merely irreducible+primitive
    // substitute fails it, so this pins the actual Conway values.
    let checked = 0;
    for (const p of available_characteristics()) {
      const degrees = new Set(available_degrees(p));
      for (const n of degrees) {
        const f = monicCoeffs(p, n);
        for (let m = 1; m < n; m++) {
          if (n % m !== 0 || !degrees.has(m)) continue;
          const e = (BigInt(p) ** BigInt(n) - 1n) / (BigInt(p) ** BigInt(m) - 1n);
          const y = polyPowMod([0, 1], e, f, p); // image of the GF(p^m) generator
          // evaluate C_{p,m} at y by Horner
          const cm = monicCoeffs(p, m);
          let acc: Poly = [];
          for (let i = cm.length - 1; i >= 0; i--) {
            acc = polyMulMod(acc, y, f, p);
            if (cm[i] !== 0) acc = polySub(acc, [p - cm[i]!], p);
          }
          expect({ p, n, m, value: acc.join(',') }).toEqual({ p, n, m, value: '' });
          checked++;
        }
      }
    }
    expect(checked).toBeGreaterThan(100);
  });

  test('constant term is (-1)^n times a primitive root of GF(p)', () => {
    // Luebeck's normalisation: f(0) = (-1)^n * g with g the smallest primitive
    // root mod p.  This is the standard sanity check on a Conway table.
    const smallestPrimitiveRoot = (p: number): number => {
      const qs = primeDivisors(BigInt(p - 1));
      for (let g = 2n; ; g++) {
        if (qs.every((q) => powmod(g, BigInt(p - 1) / q, BigInt(p)) !== 1n)) return Number(g);
      }
    };
    for (const p of available_characteristics()) {
      if (p === 2) continue; // GF(2): g = 1, sign irrelevant
      const g = smallestPrimitiveRoot(p);
      for (const n of available_degrees(p)) {
        const expected = n % 2 === 0 ? g : ((-g % p) + p) % p;
        expect({ p, n, c0: conway_polynomial(p, n)[0] }).toEqual({ p, n, c0: expected });
      }
    }
  });

  test('values reported by SageMath conway_polynomial() are reproduced', () => {
    // sage: conway_polynomial(p, n).list()  (leading 1 dropped)
    const pinned: [number, number, number[]][] = [
      [2, 2, [1, 1]],
      [2, 3, [1, 1, 0]],
      [2, 8, [1, 0, 1, 1, 1, 0, 0, 0]],
      // C2/M21 regressions: the port previously stored reducible or
      // merely-irreducible substitutes for these.
      [2, 14, [1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0]],
      [3, 8, [2, 2, 2, 0, 1, 2, 0, 0]],
      [5, 5, [3, 4, 0, 0, 0]],
      [5, 6, [2, 0, 1, 4, 1, 0]],
      [7, 5, [4, 1, 0, 0, 0]],
      [13, 4, [2, 12, 3, 0]],
      [17, 2, [3, 16]],
      [17, 3, [14, 1, 0]],
      [23, 3, [18, 2, 0]],
      [29, 2, [2, 24]],
      [31, 2, [3, 29]],
      [31, 3, [28, 1, 0]],
    ];
    for (const [p, n, coeffs] of pinned) {
      expect(conway_polynomial(p, n)).toEqual(coeffs);
    }
  });

  test('the fabricated GF(2^128) pentanomial is gone', () => {
    // x^128 + x^7 + x^2 + x + 1 was stored as "the Conway polynomial" for
    // GF(2^128); no Conway polynomial of that degree is known/tabulated, so
    // the entry must be absent rather than wrong.
    expect(has_conway_polynomial(2, 128)).toBe(false);
    expect(() => conway_polynomial(2, 128)).toThrow('No Conway polynomial in database for GF(2^128)');
  });

  test('lookup API', () => {
    expect(has_conway_polynomial(2, 2)).toBe(true);
    expect(has_conway_polynomial(2, 64)).toBe(true);
    expect(has_conway_polynomial(37, 2)).toBe(false);
    expect(() => conway_polynomial(37, 2)).toThrow(
      'No Conway polynomials in database for characteristic 37'
    );
    expect(available_characteristics()).toEqual([2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]);
    expect(available_degrees(2)[0]).toBe(2);
    expect(available_degrees(37)).toEqual([]);
    // CONWAY_POLYNOMIALS is the zero-padded view of the compact storage
    expect(CONWAY_POLYNOMIALS[2]![5]).toEqual([1, 0, 1, 0, 0]);
  });
});
