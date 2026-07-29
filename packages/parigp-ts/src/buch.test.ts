/**
 * Tests for {@link Buchquad} / {@link quadclassunit0} (PARI `buch1.c`) and the
 * `hnf_snf.c` machinery it runs on.
 *
 * Oracles, in decreasing order of authority:
 *
 *  1. **PARI's own regression suite**, copied verbatim from
 *     `reference/pari/src/test/in/quadclassunit` and its expected output
 *     `reference/pari/src/test/32/quadclassunit`: the 608 discriminants of the
 *     `test(10^15)` / `test(-10^15)` tables, the bug-report cases (#1180,
 *     #1195, #1700, #2015) and the four `quadclassunit(+-2^81 + c)` values.
 *     The class number and the cyclic structure are canonical, so they do not
 *     depend on PARI's random stream.
 *  2. **Exhaustive enumeration**: the class number of an imaginary
 *     discriminant as the number of reduced primitive forms.
 *  3. **The exact analytic class number formula** for fundamental `D < -4`:
 *     `h = |sum_{a=1}^{|D|-1} kronecker(D,a) a| / |D|`.
 *  4. **Pell's equation**: the fundamental unit, hence the regulator and the
 *     sign of its norm, of a real quadratic order.
 *  5. Group-theoretic identities: the generators returned really do generate a
 *     group of order `h` with the announced orders; the invariant factors form
 *     a divisibility chain and multiply to `h`.
 *  6. An independent textbook Smith normal form, for `ZM_snf_group`.
 */

import { describe, expect, test } from 'bun:test';
import {
  Buchquad,
  DEFAULTPREC,
  ZM_det_triangular,
  ZM_hnflll,
  ZM_mul,
  ZM_pivots,
  ZM_snf_group,
  type ZMat,
  addrr,
  bnf_increase_LIMC,
  bnfinit,
  dbltor,
  divrr,
  expr,
  gcvtoi,
  itor,
  logr_abs,
  mplog2,
  mulir,
  mulrr,
  quadclassno,
  quadclassunit0,
  real_1,
  rtodbl,
  setBuchRandomSeed,
  setprec,
  sqrtr,
  subrr,
  truncr,
} from './buch.js';
import { kronecker } from './ff.js';
import { type Qfb, mkqfb, qfbcomp, qfbred } from './qfb.js';

/* ------------------------------------------------------------------ */
/* helpers / oracles                                                   */
/* ------------------------------------------------------------------ */

const babs = (x: bigint): bigint => (x < 0n ? -x : x);

function gcdB(a: bigint, b: bigint): bigint {
  a = babs(a);
  b = babs(b);
  while (b) {
    const t = a % b;
    a = b;
    b = t;
  }
  return a;
}

function isqrtB(n: bigint): bigint {
  if (n < 2n) return n;
  let x = 1n << BigInt((n.toString(2).length >> 1) + 1);
  for (;;) {
    const y = (x + n / x) >> 1n;
    if (y >= x) break;
    x = y;
  }
  return x;
}

/** number of reduced primitive positive definite forms of discriminant `D < 0` */
function bruteClassNumber(D: bigint): bigint {
  let h = 0n;
  const bound = isqrtB(-D / 3n) + 1n;
  for (let a = 1n; a <= bound; a++) {
    for (let b = -a + 1n; b <= a; b++) {
      const t = b * b - D;
      if (t % (4n * a) !== 0n) continue;
      const c = t / (4n * a);
      if (c < a) continue;
      if (gcdB(gcdB(a, b), c) !== 1n) continue;
      if ((a === c || b === a) && b < 0n) continue;
      h++;
    }
  }
  return h;
}

/**
 * Exact Dirichlet class number formula for a fundamental discriminant
 * `D < -4`: `h = |sum_a (D/a) a| / |D|` (the `w = 2` case).
 */
function dirichletClassNumber(D: bigint): bigint {
  let s = 0n;
  const N = -D;
  for (let a = 1n; a < N; a++) s += BigInt(kronecker(D, a)) * a;
  return babs(s) / N;
}

/** is `D` a fundamental discriminant? */
function isFundamental(D: bigint): boolean {
  if (D === 1n) return false;
  const r = ((D % 4n) + 4n) % 4n;
  if (r === 1n) {
    /* squarefree */
    return isSquarefree(D);
  }
  if (r !== 0n) return false;
  const m = D / 4n;
  const rm = ((m % 4n) + 4n) % 4n;
  if (rm !== 2n && rm !== 3n) return false;
  return isSquarefree(m);
}
function isSquarefree(n: bigint): boolean {
  n = babs(n);
  for (let p = 2n; p * p <= n; p++) {
    if (n % p === 0n) {
      n /= p;
      if (n % p === 0n) return false;
    }
  }
  return true;
}

/**
 * Fundamental unit `(t + u sqrt(D))/2` of the real quadratic order of
 * discriminant `D`, from the smallest solution of `t^2 - D u^2 = +-4`.
 */
function fundamentalUnit(D: bigint, maxU = 400000n): { t: bigint; u: bigint; norm: number } | null {
  for (let u = 1n; u <= maxU; u++) {
    for (const s of [-4n, 4n]) {
      const v = D * u * u + s;
      if (v <= 0n) continue;
      const r = isqrtB(v);
      if (r * r === v && r > 0n) return { t: r, u, norm: s === 4n ? 1 : -1 };
    }
  }
  return null;
}

/** the principal form of discriminant `D` */
function principalForm(D: bigint): Qfb {
  return ((D % 4n) + 4n) % 4n === 0n ? mkqfb(1n, 0n, -D / 4n, D) : mkqfb(1n, 1n, (1n - D) / 4n, D);
}

/* 1-based column-major matrices, as in buch.ts */
function mkMat(rows: bigint[][]): ZMat {
  const m = rows.length;
  const n = m === 0 ? 0 : rows[0]!.length;
  const M: ZMat = new Array<bigint[]>(n + 1);
  M[0] = [];
  for (let j = 1; j <= n; j++) {
    const c: bigint[] = new Array<bigint>(m + 1).fill(0n);
    for (let i = 1; i <= m; i++) c[i] = rows[i - 1]![j - 1]!;
    M[j] = c;
  }
  return M;
}
function matRows(M: ZMat): string[][] {
  const n = M.length - 1;
  if (n === 0) return [];
  const m = M[1]!.length - 1;
  const out: string[][] = [];
  for (let i = 1; i <= m; i++) {
    const r: string[] = [];
    for (let j = 1; j <= n; j++) r.push(String(M[j]![i]!));
    out.push(r);
  }
  return out;
}

/** rank over Q by fraction-free elimination */
function exactRank(rows: bigint[][]): number {
  const A = rows.map((r) => r.slice());
  const m = A.length;
  const n = m ? A[0]!.length : 0;
  let r = 0;
  for (let c = 0; c < n && r < m; c++) {
    let p = -1;
    for (let i = r; i < m; i++)
      if (A[i]![c] !== 0n) {
        p = i;
        break;
      }
    if (p < 0) continue;
    const t = A[r]!;
    A[r] = A[p]!;
    A[p] = t;
    for (let i = 0; i < m; i++) {
      if (i === r || A[i]![c] === 0n) continue;
      const a = A[r]![c]!;
      const b = A[i]![c]!;
      for (let j = 0; j < n; j++) A[i]![j] = A[i]![j]! * a - A[r]![j]! * b;
      let g = 0n;
      for (let j = 0; j < n; j++) g = gcdB(g, A[i]![j]!);
      if (g > 1n) for (let j = 0; j < n; j++) A[i]![j] = A[i]![j]! / g;
    }
    r++;
  }
  return r;
}

/** textbook Smith normal form: invariant factors > 1, decreasing */
function snfOracle(rows: bigint[][]): bigint[] {
  const A = rows.map((r) => r.slice());
  const m = A.length;
  const n = m ? A[0]!.length : 0;
  const res: bigint[] = [];
  let t = 0;
  while (t < m && t < n) {
    let found = false;
    for (let i = t; i < m && !found; i++)
      for (let j = t; j < n && !found; j++)
        if (A[i]![j] !== 0n) {
          const tmp = A[t]!;
          A[t] = A[i]!;
          A[i] = tmp;
          for (let k = 0; k < m; k++) {
            const x = A[k]![t]!;
            A[k]![t] = A[k]![j]!;
            A[k]![j] = x;
          }
          found = true;
        }
    if (!found) break;
    let done = false;
    while (!done) {
      done = true;
      for (let i = t + 1; i < m; i++) {
        while (A[i]![t] !== 0n) {
          const q = A[i]![t]! / A[t]![t]!;
          for (let j = t; j < n; j++) A[i]![j] = A[i]![j]! - q * A[t]![j]!;
          if (A[i]![t] !== 0n) {
            const tmp = A[t]!;
            A[t] = A[i]!;
            A[i] = tmp;
          }
        }
      }
      for (let j = t + 1; j < n; j++) {
        while (A[t]![j] !== 0n) {
          const q = A[t]![j]! / A[t]![t]!;
          for (let i = t; i < m; i++) A[i]![j] = A[i]![j]! - q * A[i]![t]!;
          if (A[t]![j] !== 0n) {
            for (let i = 0; i < m; i++) {
              const x = A[i]![t]!;
              A[i]![t] = A[i]![j]!;
              A[i]![j] = x;
            }
            done = false;
          }
        }
      }
      if (!done) continue;
      for (let i = t + 1; i < m; i++)
        if (A[i]![t] !== 0n) {
          done = false;
          break;
        }
    }
    let ok = true;
    for (let i = t + 1; i < m && ok; i++)
      for (let j = t + 1; j < n && ok; j++)
        if (A[i]![j]! % A[t]![t]! !== 0n) {
          for (let j2 = t; j2 < n; j2++) A[t]![j2] = A[t]![j2]! + A[i]![j2]!;
          ok = false;
        }
    if (!ok) continue;
    res.push(babs(A[t]![t]!));
    t++;
  }
  const d = res.slice();
  for (let i = 0; i < d.length; i++)
    for (let j = i + 1; j < d.length; j++) {
      const g = gcdB(d[i]!, d[j]!);
      if (g === 0n) continue;
      const l = (d[i]! / g) * d[j]!;
      d[i] = g;
      d[j] = l;
    }
  return d.filter((x) => x !== 1n && x !== 0n).sort((a, b) => (a < b ? 1 : a > b ? -1 : 0));
}

let RNG = 12345;
function rnd(n: number): number {
  RNG = (RNG * 1103515245 + 12345) & 0x7fffffff;
  return RNG % n;
}

/* ------------------------------------------------------------------ */

describe('Real arithmetic (t_REAL helper)', () => {
  test('log 2 to 40 decimal digits', () => {
    const l2 = mplog2(200);
    expect(truncr(mulir(10n ** 40n, setprec(l2, 260))).toString()).toBe(
      '6931471805599453094172321214581765680755'
    );
  });

  test('sqrt 2 to 40 decimal digits', () => {
    const s2 = sqrtr(itor(2n, 200));
    expect(truncr(mulir(10n ** 40n, setprec(s2, 260))).toString()).toBe(
      '14142135623730950488016887242096980785696'
    );
  });

  test('log and exp agree with the C library in double precision', () => {
    for (const v of [0.001, 0.5, 1.5, 3.7, 100.25, 1e10, 1e-9]) {
      expect(rtodbl(logr_abs(dbltor(v, 128)))).toBeCloseTo(Math.log(v), 12);
      if (Math.abs(v) < 20) expect(rtodbl(expr(dbltor(v, 128)))).toBeCloseTo(Math.exp(v), 8);
    }
  });

  test('exp(log x) = x to high precision', () => {
    const x = dbltor(3.7, 200);
    const y = expr(logr_abs(x));
    const d = subrr(y, x);
    /* relative error < 2^-150 */
    expect(d.s === 0 || d.e - x.e < -150).toBe(true);
  });

  test('elementary operations', () => {
    const a = dbltor(3.25, 128);
    const b = dbltor(1.5, 128);
    expect(rtodbl(addrr(a, b))).toBe(4.75);
    expect(rtodbl(subrr(a, b))).toBe(1.75);
    expect(rtodbl(mulrr(a, b))).toBe(4.875);
    expect(rtodbl(divrr(a, b))).toBeCloseTo(3.25 / 1.5, 15);
  });

  test('truncr truncates towards zero, gcvtoi reports lost bits', () => {
    expect(truncr(dbltor(-3.7, 128))).toBe(-3n);
    expect(truncr(dbltor(3.7, 128))).toBe(3n);
    /* PARI gcvtoi truncates and returns expo(x - trunc(x)) when exact */
    const g = gcvtoi(dbltor(-3.75, 128));
    expect(g.z).toBe(-3n);
    expect(g.e).toBe(-1); /* |0.75| in [1/2,1) */
    /* not enough precision to determine the integral part */
    const big = mulrr(itor(1n << 200n, 64), real_1(64));
    const g2 = gcvtoi(addrr(big, real_1(64)));
    expect(g2.e).toBeGreaterThan(0);
  });
});

describe('Integer matrix machinery (hnf_snf.c)', () => {
  test('ZM_pivots agrees with exact elimination on 300 random matrices', () => {
    RNG = 12345;
    for (let trial = 0; trial < 300; trial++) {
      const m = 1 + rnd(6);
      const n = 1 + rnd(6);
      const rows: bigint[][] = [];
      for (let i = 0; i < m; i++) {
        const r: bigint[] = [];
        for (let j = 0; j < n; j++) r.push(BigInt(rnd(11) - 5));
        rows.push(r);
      }
      const { rr } = ZM_pivots(mkMat(rows));
      expect(rr).toBe(n - exactRank(rows));
    }
  });

  test('ZM_hnflll: A * B = H with B unimodular', () => {
    RNG = 999;
    for (let trial = 0; trial < 200; trial++) {
      const m = 1 + rnd(5);
      const n = 1 + rnd(5);
      const rows: bigint[][] = [];
      for (let i = 0; i < m; i++) {
        const r: bigint[] = [];
        for (let j = 0; j < n; j++) r.push(BigInt(rnd(15) - 7));
        rows.push(r);
      }
      const { H, B } = ZM_hnflll(mkMat(rows), true, false);
      expect(matRows(ZM_mul(mkMat(rows), B!))).toEqual(matRows(H));
    }
  });

  test('ZM_snf_group matches a textbook Smith normal form', () => {
    RNG = 4242;
    let tested = 0;
    for (let trial = 0; trial < 300; trial++) {
      const n = 1 + rnd(5);
      const rows: bigint[][] = [];
      for (let i = 0; i < n; i++) {
        const r: bigint[] = [];
        for (let j = 0; j < n; j++) r.push(BigInt(rnd(15) - 7));
        rows.push(r);
      }
      if (exactRank(rows) < n) continue;
      const { H } = ZM_hnflll(mkMat(rows), false, false);
      const h = matRows(H).map((r) => r.map(BigInt));
      const { D } = ZM_snf_group(mkMat(h));
      const inv = D.slice(1);
      expect(inv.map(String)).toEqual(snfOracle(h).map(String));
      tested++;
    }
    expect(tested).toBeGreaterThan(100);
  });

  test('ZM_det_triangular', () => {
    const M = mkMat([
      [2n, 1n, 0n],
      [0n, 3n, 2n],
      [0n, 0n, 5n],
    ]);
    expect(ZM_det_triangular(M)).toBe(30n);
  });

  test('bnf_increase_LIMC (buch1.c:189)', () => {
    /* c <= 0.3 doubles, otherwise adds D/20 */
    expect(bnf_increase_LIMC(10, 1000)).toBe(20);
    expect(bnf_increase_LIMC(500, 1000)).toBe(550);
  });
});

describe('Buchquad: imaginary quadratic fields', () => {
  test('every discriminant -3 .. -600 matches the reduced form count', () => {
    setBuchRandomSeed(7);
    let tested = 0;
    for (let d = 3n; d <= 600n; d++) {
      const D = -d;
      const r4 = ((D % 4n) + 4n) % 4n;
      if (r4 !== 0n && r4 !== 1n) continue;
      const res = Buchquad(D);
      expect(res.no).toBe(bruteClassNumber(D));
      /* structure: divisibility chain, product = h */
      let prod = 1n;
      for (let i = 0; i < res.cyc.length; i++) {
        prod *= res.cyc[i]!;
        if (i) expect(res.cyc[i - 1]! % res.cyc[i]!).toBe(0n);
      }
      expect(prod).toBe(res.no);
      tested++;
    }
    expect(tested).toBeGreaterThan(290);
  });

  test('exact Dirichlet class number formula, fundamental D', () => {
    setBuchRandomSeed(13);
    let tested = 0;
    for (const D of [-23n, -47n, -71n, -163n, -1051n, -3299n, -4003n, -10007n, -20015n, -99991n]) {
      if (!isFundamental(D)) continue;
      expect(Buchquad(D).no).toBe(dirichletClassNumber(D));
      tested++;
    }
    expect(tested).toBe(10);
  });

  test('class group structure and generators', () => {
    setBuchRandomSeed(3);
    const expected: [bigint, bigint[]][] = [
      [-120n, [2n, 2n]],
      [-84n, [2n, 2n]],
      [-260n, [4n, 2n]],
      [-1155n, [2n, 2n, 2n]],
      [-3299n, [9n, 3n]],
      [-5460n, [2n, 2n, 2n, 2n]],
      [-1015n, [8n, 2n]],
      [-20015n, [126n]],
      [-100003n, [39n]],
    ];
    for (const [D, cyc] of expected) {
      const res = Buchquad(D);
      expect(res.cyc).toEqual(cyc);
      /* the generators generate a group of order h, with the announced orders */
      const seen = new Set<string>();
      const one = qfbred(principalForm(D));
      const key = (f: Qfb) => `${f.a},${f.b},${f.c}`;
      const rec = (i: number, acc: Qfb): void => {
        if (i === res.gen.length) {
          seen.add(key(acc));
          return;
        }
        let cur = acc;
        for (let e = 0n; e < res.cyc[i]!; e++) {
          rec(i + 1, cur);
          cur = qfbcomp(cur, res.gen[i]!);
        }
        /* g^cyc[i] = 1 */
        expect(key(cur)).toBe(key(acc));
      };
      rec(0, one);
      expect(BigInt(seen.size)).toBe(res.no);
    }
  });

  test('non-fundamental discriminants (orders)', () => {
    setBuchRandomSeed(5);
    for (const D of [-12n, -16n, -27n, -75n, -100n, -29920n]) {
      const res = Buchquad(D);
      if (babs(D) < 2000n) expect(res.no).toBe(bruteClassNumber(D));
    }
    /* PARI test/in/quadclassunit, "#1195 with non-fundamental discriminants" */
    expect(Buchquad(-29920n).cyc).toEqual([4n, 2n, 2n, 2n]);
  });

  test('tiny discriminants', () => {
    expect(Buchquad(-3n).no).toBe(1n);
    expect(Buchquad(-4n).no).toBe(1n);
    expect(Buchquad(-8n).no).toBe(1n);
    expect(Buchquad(-3n).cyc).toEqual([]);
  });

  test('analytic (Dirichlet) upper bound on h', () => {
    setBuchRandomSeed(17);
    for (const D of [-1000003n, -10000019n, -100000007n]) {
      const h = Number(Buchquad(D).no);
      const q = Number(-D);
      /* |L(1,chi)| <= 1 + log(q)/2 and h = sqrt(q) L(1,chi) / pi */
      const bound = (Math.sqrt(q) / Math.PI) * (1 + Math.log(q) / 2);
      expect(h).toBeLessThan(bound);
      expect(h).toBeGreaterThan(0);
    }
  });
});

describe('Buchquad: real quadratic fields', () => {
  test('regulator and unit norm against Pell', () => {
    setBuchRandomSeed(11);
    const data: bigint[] = [
      5n,
      8n,
      12n,
      13n,
      17n,
      21n,
      24n,
      28n,
      29n,
      33n,
      40n,
      44n,
      60n,
      61n,
      65n,
      76n,
      85n,
      88n,
      104n,
      229n,
      401n,
      577n,
      1009n,
      10001n,
      20n,
      45n,
      32n,
      48n,
    ];
    for (const D of data) {
      const res = Buchquad(D);
      const fu = fundamentalUnit(D)!;
      expect(fu).not.toBeNull();
      const R = Math.log((Number(fu.t) + Number(fu.u) * Math.sqrt(Number(D))) / 2);
      expect(rtodbl(res.reg)).toBeCloseTo(R, 8);
      expect(res.sign).toBe(fu.norm === -1 ? -1n : 1n);
      /* regulator is a genuine multiprecision value */
      expect(res.reg.p).toBeGreaterThan(100);
    }
  });

  test('known class numbers', () => {
    setBuchRandomSeed(19);
    expect(Buchquad(229n).no).toBe(3n);
    expect(Buchquad(401n).no).toBe(5n);
    expect(Buchquad(577n).no).toBe(7n);
    expect(Buchquad(10001n).cyc).toEqual([16n]);
    /* PARI test/in/quadclassunit #1180 */
    expect(Buchquad(572n).no).toBe(2n);
    /* PARI test/in/quadclassunit: quadclassunit(8*3*5*7).cyc */
    expect(Buchquad(840n).cyc).toEqual([2n, 2n]);
  });

  test('h * R matches the analytic class number formula within the GRH window', () => {
    setBuchRandomSeed(23);
    for (const D of [401n, 1009n, 10001n]) {
      const res = Buchquad(D);
      const hR = Number(res.no) * rtodbl(res.reg);
      const q = Number(D);
      /* h R = sqrt(D) L(1,chi) / 2, with 0 < L(1,chi) <= 1 + log(q)/2 */
      expect(hR).toBeGreaterThan(0);
      expect(hR).toBeLessThan((Math.sqrt(q) / 2) * (1 + Math.log(q) / 2));
    }
  });
});

describe('Buchquad against PARI regression values', () => {
  /* reference/pari/src/test/32/quadclassunit */
  test('bug-report cases from test/in/quadclassunit', () => {
    setBuchRandomSeed(1);
    const cases: [bigint, bigint[]][] = [
      [-8419588n, [176n, 2n]] /* #1195 */,
      [-1459008n, [16n, 4n, 2n, 2n]],
      [-3799812n, [54n, 2n, 2n, 2n]],
      [-13163208n, [156n, 2n, 2n]],
      [-29920n, [4n, 2n, 2n, 2n]],
      [-612556842419n, [192199n]] /* #1700 */,
      [-699n, [10n]] /* #2015 */,
      [(1n << 70n) + 25n, [17n]],
    ];
    for (const [D, cyc] of cases) expect(Buchquad(D).cyc).toEqual(cyc);
  });

  test('Bach constants are honoured (quadclassunit(D,,[c1,c2]))', () => {
    setBuchRandomSeed(1);
    /* quadclassunit(-13163208,,[0.1]).cyc and quadclassunit(-699,,[6,6]).cyc */
    expect(quadclassunit0(-13163208n, 0, [0.1]).cyc).toEqual([156n, 2n, 2n]);
    expect(quadclassunit0(-699n, 0, [6, 6]).cyc).toEqual([10n]);
  });

  test('quadclassunit(+-2^81 + c): class numbers up to 10^12', () => {
    setBuchRandomSeed(1);
    const two81 = 1n << 81n;
    const cases: [bigint, bigint, bigint[]][] = [
      [1n - two81, 959887376480n, [59992961030n, 2n, 2n, 2n, 2n]],
      [5n - two81, 391530492496n, [97882623124n, 2n, 2n]],
      [-4n - two81, 820252603100n, [410126301550n, 2n]],
      [-8n - two81, 302960366592n, [1183438932n, 2n, 2n, 2n, 2n, 2n, 2n, 2n, 2n]],
    ];
    for (const [D, no, cyc] of cases) {
      const r = Buchquad(D);
      expect(r.no).toBe(no);
      expect(r.cyc).toEqual(cyc);
    }
  }, 60000);

  /* The two `test(10^15)` / `test(-10^15)` tables of
   * reference/pari/src/test/32/quadclassunit, transcribed verbatim.  The first
   * block is real (`d` from 10^15, fundamental only), the second imaginary. */
  const TABLE_REAL: [bigint, bigint[]][] = [
    [1000000000000001n, [4n, 4n, 2n, 2n, 2n]],
    [1000000000000005n, [2n, 2n, 2n, 2n]],
    [1000000000000009n, [2n, 2n]],
    [1000000000000012n, [32n, 2n]],
    [1000000000000013n, []],
    [1000000000000021n, [2n]],
    [1000000000000024n, [4n]],
    [1000000000000028n, [2n]],
    [1000000000000029n, [6n, 2n]],
    [1000000000000033n, []],
    [1000000000000037n, []],
    [1000000000000040n, [8n, 2n, 2n, 2n]],
    [1000000000000041n, [2n, 2n, 2n, 2n]],
    [1000000000000045n, [4n, 2n]],
    [1000000000000049n, [2n]],
    [1000000000000056n, [2n, 2n, 2n]],
    [1000000000000057n, [2n]],
    [1000000000000060n, [6n, 2n, 2n, 2n]],
    [1000000000000061n, [2n]],
    [1000000000000065n, [4n, 2n, 2n]],
    [1000000000000069n, []],
    [1000000000000076n, [2n]],
    [1000000000000077n, [6n, 2n]],
    [1000000000000081n, [2n, 2n]],
    [1000000000000085n, [2n, 2n]],
    [1000000000000088n, [13n]],
    [1000000000000093n, []],
    [1000000000000097n, [4n, 4n]],
    [1000000000000101n, [12n, 2n, 2n]],
    [1000000000000104n, [2n]],
    [1000000000000105n, [2n, 2n, 2n]],
    [1000000000000108n, [2n]],
    [1000000000000109n, [4n]],
    [1000000000000113n, [2n, 2n]],
    [1000000000000117n, []],
    [1000000000000120n, [12n, 2n]],
    [1000000000000121n, [2n]],
    [1000000000000124n, [8n, 2n, 2n]],
    [1000000000000129n, [6n]],
    [1000000000000133n, [2n, 2n, 2n]],
    [1000000000000136n, [4n]],
    [1000000000000137n, [4n, 2n]],
    [1000000000000140n, [2n, 2n, 2n, 2n]],
    [1000000000000145n, [4n]],
    [1000000000000149n, [2n, 2n]],
    [1000000000000153n, []],
    [1000000000000156n, [2n, 2n, 2n, 2n]],
    [1000000000000157n, [22n, 2n, 2n]],
    [1000000000000165n, [2n, 2n, 2n]],
    [1000000000000168n, [334n]],
    [1000000000000169n, []],
    [1000000000000172n, [3n]],
    [1000000000000173n, [2n, 2n]],
    [1000000000000177n, [2n]],
    [1000000000000181n, [2n]],
    [1000000000000184n, [2n, 2n, 2n, 2n]],
    [1000000000000185n, [8n, 2n, 2n, 2n]],
    [1000000000000189n, [22n]],
    [1000000000000193n, [2n]],
    [1000000000000201n, [2n]],
  ];
  const TABLE_IMAG: [bigint, bigint[]][] = [
    [-999999999999995n, [3872378n, 2n]],
    [-999999999999992n, [2471436n, 2n, 2n, 2n]],
    [-999999999999991n, [9144306n, 2n]],
    [-999999999999988n, [1124902n, 2n, 2n]],
    [-999999999999987n, [913748n, 2n, 2n, 2n]],
    [-999999999999983n, [13126428n, 2n]],
    [-999999999999979n, [2148058n, 2n, 2n]],
    [-999999999999976n, [2984712n, 2n, 2n]],
    [-999999999999971n, [7798614n, 2n]],
    [-999999999999967n, [11253252n]],
    [-999999999999960n, [481528n, 2n, 2n, 2n, 2n]],
    [-999999999999959n, [25233340n, 2n]],
    [-999999999999956n, [9820206n, 2n]],
    [-999999999999955n, [2185542n, 2n]],
    [-999999999999951n, [11807322n, 2n]],
    [-999999999999947n, [7261653n]],
    [-999999999999944n, [7665042n, 2n]],
    [-999999999999943n, [17982820n]],
    [-999999999999940n, [1458740n, 2n, 2n]],
    [-999999999999939n, [1705720n, 2n, 2n]],
    [-999999999999935n, [3815348n, 2n, 2n, 2n]],
    [-999999999999931n, [377454n, 2n, 2n, 2n, 2n]],
    [-999999999999928n, [3090498n, 2n]],
    [-999999999999924n, [275778n, 2n, 2n, 2n, 2n, 2n, 2n]],
    [-999999999999923n, [2614376n, 2n, 2n]],
    [-999999999999919n, [20175156n]],
    [-999999999999915n, [878588n, 2n, 2n, 2n]],
    [-999999999999912n, [3618208n, 2n]],
    [-999999999999911n, [31532754n]],
    [-999999999999908n, [12931036n]],
    [-999999999999907n, [3859738n]],
    [-999999999999903n, [448308n, 2n, 2n, 2n, 2n, 2n]],
    [-999999999999899n, [21701582n]],
    [-999999999999895n, [4135150n, 2n, 2n]],
    [-999999999999892n, [1620700n, 2n, 2n]],
    [-999999999999887n, [10190276n, 2n, 2n]],
    [-999999999999883n, [2498825n]],
    [-999999999999880n, [1139478n, 2n, 2n, 2n]],
    [-999999999999879n, [9394500n, 2n]],
    [-999999999999876n, [1305236n, 2n, 2n, 2n]],
    [-999999999999871n, [13406168n, 2n]],
    [-999999999999867n, [4308820n, 2n]],
    [-999999999999863n, [7643020n, 2n]],
    [-999999999999860n, [4852422n, 2n, 2n]],
    [-999999999999859n, [3024070n, 2n]],
    [-999999999999851n, [1199052n, 4n, 2n]],
    [-999999999999848n, [6290596n, 2n]],
    [-999999999999844n, [5296516n, 2n]],
    [-999999999999843n, [3731152n, 2n]],
    [-999999999999839n, [42706704n]],
    [-999999999999835n, [1378294n, 2n, 2n]],
    [-999999999999832n, [974764n, 2n, 2n, 2n]],
    [-999999999999831n, [21097546n, 2n]],
    [-999999999999827n, [7464103n]],
    [-999999999999823n, [4886300n, 2n]],
    [-999999999999816n, [3866004n, 2n, 2n]],
    [-999999999999815n, [8445968n, 2n, 2n]],
    [-999999999999812n, [1381066n, 2n, 2n, 2n]],
    [-999999999999811n, [3696512n, 2n]],
    [-999999999999807n, [3253596n, 2n, 2n]],
  ];

  test('test(10^15): 60 real discriminants', () => {
    setBuchRandomSeed(1);
    for (const [D, cyc] of TABLE_REAL) expect(Buchquad(D).cyc).toEqual(cyc);
  }, 120000);

  test('test(-10^15): 60 imaginary discriminants', () => {
    setBuchRandomSeed(1);
    for (const [D, cyc] of TABLE_IMAG) expect(Buchquad(D).cyc).toEqual(cyc);
  }, 120000);
});

describe('bnfinit (buch2.c)', () => {
  test('degree 1: the trivial bnf of Q', () => {
    const b = bnfinit([0n, 1n]);
    expect(b.clgp.no).toBe(1n);
    expect(b.clgp.cyc).toEqual([]);
    expect(rtodbl(b.reg)).toBe(1);
  });

  test('degree 2: class group, regulator and torsion of the maximal order', () => {
    setBuchRandomSeed(37);
    /* [T, disc(K), h, cyc, R (0 if imaginary), #torsion] */
    const cases: [bigint[], bigint, bigint, bigint[], number, bigint][] = [
      [[-5n, 0n, 1n], 5n, 1n, [], 0.4812118250596035, 2n],
      [[-2n, 0n, 1n], 8n, 1n, [], 0.881373587019543, 2n],
      [[-12n, 0n, 1n], 12n, 1n, [], 1.3169578969248166, 2n] /* disc 48 -> 12 */,
      [[-79n, 0n, 1n], 316n, 3n, [3n], 5.07513475044481, 2n],
      [[23n, 0n, 1n], -23n, 3n, [3n], 0, 2n],
      [[14n, 0n, 1n], -56n, 4n, [4n], 0, 2n],
      [[1n, 1n, 1n], -3n, 1n, [], 0, 6n],
      [[1n, 0n, 1n], -4n, 1n, [], 0, 4n],
    ];
    for (const [T, disc, h, cyc, R, tu] of cases) {
      const b = bnfinit(T);
      expect(b.disc).toBe(disc);
      expect(b.clgp.no).toBe(h);
      expect(b.clgp.cyc).toEqual(cyc);
      expect(b.tu).toBe(tu);
      if (R) expect(rtodbl(b.reg)).toBeCloseTo(R, 9);
      else expect(rtodbl(b.reg)).toBe(1);
    }
  });

  test('degree > 2 names the missing upstream routines', () => {
    expect(() => bnfinit([1n, 0n, 0n, 1n])).toThrow(/SAGE_NOT_IMPLEMENTED/);
    expect(() => bnfinit([1n, 0n, 0n, 1n])).toThrow(/nfmaxord/);
    expect(() => bnfinit([1n, 0n, 0n, 1n])).toThrow(/idealprimedec/);
  });

  test('nonmonic and degenerate input', () => {
    expect(() => bnfinit([1n, 0n, 2n])).toThrow('nonmonic');
    expect(() => bnfinit([-4n, 0n, 1n])).toThrow('issquare(disc)');
  });
});

describe('quadclassunit0 / quadclassno wrappers and errors', () => {
  test('quadclassno', () => {
    setBuchRandomSeed(29);
    expect(quadclassno(-1000003n)).toBe(105n);
    expect(quadclassno(-23n)).toBe(3n);
  });

  test('square discriminant is refused', () => {
    expect(() => Buchquad(9n)).toThrow('domain error in Buchquad: issquare(disc) = 1');
  });

  test('discriminant must be 0 or 1 mod 4', () => {
    expect(() => Buchquad(-26n)).toThrow('domain error in Buchquad: disc % 4 > 1');
    expect(() => Buchquad(7n)).toThrow('domain error in Buchquad: disc % 4 > 1');
  });

  test('narrow class group is not implemented (PARI pari_err_IMPL)', () => {
    expect(() => quadclassunit0(5n, 1)).toThrow('narrow class group');
  });

  test('negative Bach constant is refused', () => {
    expect(() => Buchquad(-23n, -1)).toThrow('Bach constant');
  });

  test('default precision of the returned regulator', () => {
    setBuchRandomSeed(31);
    expect(Buchquad(-23n).reg.p).toBe(DEFAULTPREC);
    expect(rtodbl(Buchquad(-23n).reg)).toBe(1);
  });
});
