/**
 * Tests for the SEA port (`ellsea.ts`).
 *
 * Oracles used here:
 * - exhaustive point counting over `F_p` for small `p` (an independent,
 *   obviously-correct implementation);
 * - `ellcard` (Shanks/Mestre, `group.ts`) and `Fp_ellcard_Schoof`
 *   (`advanced.ts`), which share no code with `ellsea.ts`;
 * - `Fp_polmodular_evalx` (`polmodular.ts`), an independent route to
 *   `Phi_L(X, j) mod p` and its derivatives;
 * - the published cardinalities of NIST P-256 and Curve25519;
 * - a *proof* of the answer for large `p`: exhibit a point `P` with `N P = O`
 *   and `ord(P) > 4 sqrt(p)`, so that `N` is the only multiple of `ord(P)` in
 *   the Hasse interval (`assertCardProven` below).
 */

import { describe, expect, test } from 'bun:test';
import {
  Fp_ellcard_SEA,
  Fp_elljissupersingular,
  Fq_elldivpolmod,
  _internal,
} from './ellsea.js';
import { ellcard_sea, Fp_ellcard_Schoof } from './advanced.js';
import {
  FpE_mul,
  FpE_random,
  ell_is_inf,
  ellcard,
  ellinit_Fp,
  ellorder,
} from './group.js';
import { Fp_polmodular_evalx, INV_J } from '../polmodular.js';
import { kronecker } from '../ff.js';

/* ------------------------------------------------------------------ */
/* helpers                                                            */
/* ------------------------------------------------------------------ */

function mod(a: bigint, p: bigint): bigint {
  const r = a % p;
  return r < 0n ? r + p : r;
}
function Fp_inv(a: bigint, p: bigint): bigint {
  let [r0, r1] = [mod(a, p), p];
  let [s0, s1] = [1n, 0n];
  while (r1 !== 0n) {
    const q = r0 / r1;
    [r0, r1] = [r1, r0 - q * r1];
    [s0, s1] = [s1, s0 - q * s1];
  }
  return mod(s0, p);
}
/** #E(F_p) by exhaustive enumeration of the x-coordinates. */
function bruteCard(a4: bigint, a6: bigint, p: bigint): bigint {
  let n = 1n;
  for (let x = 0n; x < p; x++) {
    const r = mod(((x * x) % p) * x + a4 * x + a6, p);
    if (r === 0n) n += 1n;
    else n += kronecker(r, p) > 0 ? 2n : 0n;
  }
  return n;
}
function ellj(a4: bigint, a6: bigint, p: bigint): bigint {
  const a43 = mod(((4n * a4 * a4) % p) * a4, p);
  return (((a43 * 1728n) % p) * Fp_inv(mod(a43 + 27n * a6 * a6, p), p)) % p;
}
function isqrt(n: bigint): bigint {
  if (n < 2n) return n;
  let x = n;
  let y = (x + 1n) / 2n;
  while (y < x) {
    x = y;
    y = (x + n / x) / 2n;
  }
  return x;
}
/**
 * Prove that `N = #E(F_p)`: find a point of order `> 4 sqrt(p)` killed by `N`;
 * then `N` is the only multiple of that order inside the Hasse interval.
 */
function assertCardProven(a4: bigint, a6: bigint, p: bigint, N: bigint): void {
  const E = ellinit_Fp(a4, a6, p);
  const s = isqrt(p);
  expect((p + 1n - N) * (p + 1n - N) <= 4n * p).toBe(true); /* Hasse */
  for (let k = 0; k < 30; k++) {
    const P = FpE_random(E);
    expect(ell_is_inf(FpE_mul(P, N, a4, p))).toBe(true);
    const o = ellorder(E, P, N);
    if (o <= 4n * s + 4n) continue;
    const lo = p + 1n - 2n * s - 1n;
    const hi = p + 1n + 2n * s + 1n;
    let cnt = 0n;
    for (let m = (lo / o) * o; m <= hi; m += o) if (m >= lo) cnt++;
    if (cnt === 1n) return; /* proven */
  }
  throw new Error('assertCardProven: no point of large enough order found');
}


/* ------------------------------------------------------------------ */

describe('FpX helpers (ellsea.c support routines)', () => {
  const p = 1000003n;

  test('FpXn_inv is the truncated inverse power series', () => {
    for (let k = 1; k <= 20; k++) {
      const f = [1n + BigInt(k), 3n, 5n, 7n, 11n, 13n].map((c) => mod(c * BigInt(k), p));
      f[0] = mod(f[0]! + 1n, p);
      const n = 6;
      const g = _internal.FpXn_inv(f, n, p);
      const prod = _internal.FpX_mul(f, g, p).slice(0, n);
      while (prod.length && prod[prod.length - 1] === 0n) prod.pop();
      expect(prod).toEqual([1n]);
    }
  });

  test('FpXn_expint(h, e) solves f’ = h f, f(0) = 1', () => {
    for (let k = 1; k <= 10; k++) {
      const h = [BigInt(k), 2n, 3n, 4n, 5n].map((c) => mod(c * 7n + BigInt(k), p));
      const e = 6;
      const f = _internal.FpXn_expint(h, e, p);
      /* f' - h f = 0 mod x^(e-1) */
      const fp: bigint[] = [];
      for (let i = 1; i < f.length; i++) fp.push(mod(f[i]! * BigInt(i), p));
      const hf = _internal.FpX_mul(h, f, p);
      const diff: bigint[] = [];
      for (let i = 0; i < e - 1; i++) diff.push(mod((fp[i] ?? 0n) - (hf[i] ?? 0n), p));
      while (diff.length && diff[diff.length - 1] === 0n) diff.pop();
      expect(diff).toEqual([]);
      expect(f[0]).toBe(1n);
    }
  });

  test('FpX_roots finds exactly the roots', () => {
    for (const q of [1009n, 1000003n]) {
      for (let k = 1; k <= 8; k++) {
        const f = [BigInt(k), 5n, 0n, 3n, 1n].map((c) => mod(c, q));
        const r = _internal.FpX_roots(f, q);
        const brute: bigint[] = [];
        for (let x = 0n; x < (q < 5000n ? q : 5000n); x++) {
          let v = 0n;
          for (let i = f.length - 1; i >= 0; i--) v = mod(v * x + f[i]!, q);
          if (v === 0n) brute.push(x);
        }
        if (q < 5000n) expect(r).toEqual(brute);
        else for (const x of brute) expect(r).toContain(x);
      }
    }
  });

  test('FpX_ddf_degree returns the order of the Frobenius', () => {
    const q = 1009n;
    /* products of irreducible polynomials of equal degree */
    for (const f of [
      [1n, 1n, 1n], // x^2+x+1: irreducible mod 1009? checked below by brute force
      [2n, 0n, 1n],
      [1n, 0n, 0n, 1n],
    ]) {
      const T = _internal.FpX_red(f, q);
      const XP = _internal.FpX_Frobenius(T, q);
      const r = _internal.FpX_ddf_degree(T, XP, q);
      /* brute force: smallest r with x^(p^r) = x mod T */
      let brute = 0;
      let cur = XP;
      for (let i = 1; i <= 6; i++) {
        if (cur.length === 2 && cur[0] === 0n && cur[1] === 1n) {
          brute = i;
          break;
        }
        /* cur <- cur^p */
        let acc: bigint[] = [];
        let pw: bigint[] = [1n];
        for (let d = 0; d < cur.length; d++) {
          acc = _internal.FpX_red(
            addPoly(acc, mulScalar(pw, cur[d]!, q), q),
            q
          );
          pw = _internal.FpX_rem(_internal.FpX_mul(pw, XP, q), T, q);
        }
        cur = acc;
      }
      expect(r).toBe(brute);
    }
    function addPoly(a: bigint[], b: bigint[], m: bigint): bigint[] {
      const n = Math.max(a.length, b.length);
      const out: bigint[] = [];
      for (let i = 0; i < n; i++) out.push(mod((a[i] ?? 0n) + (b[i] ?? 0n), m));
      return out;
    }
    function mulScalar(a: bigint[], c: bigint, m: bigint): bigint[] {
      return a.map((x) => mod(x * c, m));
    }
  });

  test('Fp2_sqrt returns a genuine square root in F_{p^2}', () => {
    const q = 1009n; /* g = Fp_2gener(1009) */
    const g = 11n; /* any non-residue works for the algebra; check first */
    expect(kronecker(g, q)).toBe(-1);
    for (let a = 1n; a <= 30n; a++)
      for (let b = 0n; b <= 3n; b++) {
        const x: [bigint, bigint] = [a, b];
        const sq: [bigint, bigint] = [
          mod(x[0] * x[0] + x[1] * x[1] * g, q),
          mod(2n * x[0] * x[1], q),
        ];
        const r = _internal.Fp2_sqrt(sq, g, q);
        expect(r).not.toBeNull();
        const back: [bigint, bigint] = [
          mod(r![0] * r![0] + r![1] * r![1] * g, q),
          mod(2n * r![0] * r![1], q),
        ];
        expect(back).toEqual(sq);
      }
  }, 60000);

  test('Z_incremental_CRT agrees with a direct CRT and centres the result', () => {
    let H = 1n;
    let Q = 4n;
    const mods = [27, 25, 49];
    const vals = [5, 13, 30];
    for (let i = 0; i < mods.length; i++) {
      const r = _internal.Z_incremental_CRT(H, vals[i]!, Q, mods[i]!);
      H = r.H;
      Q = r.q;
      expect(mod(H, BigInt(mods[i]!))).toBe(BigInt(vals[i]!));
      expect(H > -Q / 2n - 1n && H <= Q / 2n).toBe(true);
    }
    expect(mod(H, 4n)).toBe(1n);
    expect(Q).toBe(4n * 27n * 25n * 49n);
  });

  test('possible_traces enumerates every CRT combination', () => {
    const compile = [
      { mod: 5n, traces: [1, 2] },
      { mod: 7n, traces: [3, 4, 5] },
      { mod: 11n, traces: [0] },
    ];
    const { V, P } = _internal.possible_traces(compile, 0b111);
    expect(P).toBe(385n);
    expect(V.length).toBe(6);
    const seen = new Set(V.map(String));
    expect(seen.size).toBe(6);
    for (const v of V) {
      expect([1n, 2n]).toContain(mod(v, 5n));
      expect([3n, 4n, 5n]).toContain(mod(v, 7n));
      expect(mod(v, 11n)).toBe(0n);
    }
  });
});

describe('modular equations', () => {
  test('Fq_polmodular_eval reproduces Phi_L(X, j) and its j-derivatives', () => {
    const p = 1000003n;
    for (const ell of [3, 5, 7, 11, 13, 17, 19, 23]) {
      for (const j of [12345n, 777n, 999999n]) {
        const M = _internal.get_modular_eqn(ell);
        const [R, dR, ddR] = _internal.Fq_polmodular_eval(M.eq, j, ell, p);
        const [A, dA, ddA] = Fp_polmodular_evalx(ell, INV_J, j, p, true) as [
          bigint[],
          bigint[],
          bigint[],
        ];
        const norm = (v: bigint[]) => {
          const w = v.slice();
          while (w.length && w[w.length - 1] === 0n) w.pop();
          return w;
        };
        expect(norm(R)).toEqual(norm(A));
        expect(norm(dR)).toEqual(norm(dA));
        expect(norm(ddR)).toEqual(norm(ddA));
        expect(R.length - 1).toBe(ell + 1); /* monic of degree ell+1 */
        expect(R[ell + 1]).toBe(1n);
      }
    }
  }, 300000);

  test('the kernel polynomial divides the ell-division polynomial', () => {
    /* upstream's own check, ellsea.c:1300 */
    const cases: Array<[bigint, bigint, bigint]> = [
      [340282366920938463463374607431768211297n, 7932n, 104734n],
      [1000000007n, 12345n, 6789n],
      [1000003n, 3n, 5n],
    ];
    let elkies = 0;
    for (const [p, a4, a6] of cases) {
      const j = ellj(a4, a6, p);
      for (const ell of [3, 5, 7, 11, 13, 17, 19, 23]) {
        const M = _internal.get_modular_eqn(ell);
        const mj = _internal.meqn_j(M, j, ell, p);
        const st = _internal.study_modular_eqn(ell, mj, p);
        if (st.mt !== _internal.ModType.MTElkies && st.mt !== _internal.ModType.MTone_root)
          continue;
        const iso = _internal.find_isogenous(a4, a6, ell, M, st.g!, p);
        expect(iso).not.toBeNull();
        expect(iso!.h.length - 1).toBe((ell - 1) >> 1);
        expect(Fq_elldivpolmod(a4, a6, ell, iso!.h, p)).toEqual([]);
        elkies++;
      }
    }
    expect(elkies).toBeGreaterThan(8);
  }, 300000);

  test('find_trace_Atkin always contains the true trace', () => {
    const p = 1000003n;
    let atkin = 0;
    for (let a4 = 1n; a4 <= 6n; a4++)
      for (let a6 = 1n; a6 <= 6n; a6++) {
        if (mod(4n * a4 * a4 * a4 + 27n * a6 * a6, p) === 0n) continue;
        const t = p + 1n - ellcard({ a4, a6, p });
        const j = ellj(a4, a6, p);
        if (j === 0n || j === mod(1728n, p)) continue;
        for (const ell of [3, 5, 7, 11, 13]) {
          const M = _internal.get_modular_eqn(ell);
          const mj = _internal.meqn_j(M, j, ell, p);
          const st = _internal.study_modular_eqn(ell, mj, p);
          if (st.mt !== _internal.ModType.MTAtkin) continue;
          const tr = _internal.find_trace_Atkin(ell, st.r, p);
          expect(tr.length).toBeGreaterThan(0);
          expect(tr).toContain(Number(mod(t, BigInt(ell))));
          atkin++;
        }
      }
    expect(atkin).toBeGreaterThan(10);
  }, 300000);

  test('the l+1 roots branch (MTroots) is exercised and correct', () => {
    /* Phi_ell(X, j) splits completely: the trace is known modulo ell^2 */
    const p = 1000003n;
    let found = 0;
    for (const [a4, a6, ell] of [
      [12n, 8n, 3],
      [13n, 3n, 3],
      [16n, 2n, 3],
      [12n, 9n, 7],
      [18n, 3n, 7],
    ] as Array<[bigint, bigint, number]>) {
      const j = ellj(a4, a6, p);
      const M = _internal.get_modular_eqn(ell);
      const mj = _internal.meqn_j(M, j, ell, p);
      const st = _internal.study_modular_eqn(ell, mj, p);
      expect(st.mt).toBe(_internal.ModType.MTroots);
      const t = p + 1n - ellcard({ a4, a6, p });
      const tr = _internal.find_trace_lp1_roots(ell, p);
      expect(tr).toContain(Number(mod(t, BigInt(ell * ell))));
      expect(Fp_ellcard_SEA(a4, a6, p)).toBe(ellcard({ a4, a6, p }));
      found++;
    }
    expect(found).toBe(5);
  }, 300000);

  test('find_trace_one_root / find_trace_lp1_roots satisfy their defining congruence', () => {
    const q = 1000003n;
    for (const ell of [5, 7, 11, 13, 17]) {
      const L = BigInt(ell);
      if (kronecker(mod(q, L), L) < 0) continue;
      const tr = _internal.find_trace_one_root(ell, q);
      for (const t of tr) expect(mod(BigInt(t) * BigInt(t) - 4n * q, L)).toBe(0n);
      const tr2 = _internal.find_trace_lp1_roots(ell, q);
      for (const t of tr2) expect(mod(BigInt(t) * BigInt(t) - 4n * q, L * L)).toBe(0n);
    }
  }, 60000);
});

describe('supersingularity (FpE.c:753-805)', () => {
  test('Fp_elljissupersingular matches exhaustive point counting', () => {
    for (const p of [1009n, 2003n]) {
      const ssTrue: bigint[] = [];
      const ssPred: bigint[] = [];
      for (let j = 0n; j < p; j++) {
        if (j === 0n || j === mod(1728n, p)) continue;
        const k = mod(1728n - j, p);
        const a4 = (((3n * j) % p) * k) % p;
        const a6 = ((((2n * j) % p) * k) % p) * k % p;
        if (mod(4n * a4 * a4 * a4 + 27n * a6 * a6, p) === 0n) continue;
        if (bruteCard(a4, a6, p) === p + 1n) ssTrue.push(j);
        if (Fp_elljissupersingular(j, p)) ssPred.push(j);
      }
      expect(ssPred).toEqual(ssTrue);
      expect(ssTrue.length).toBeGreaterThan(0);
    }
  }, 120000);

  test('SEA returns p+1 on supersingular curves', () => {
    const p = 1009n;
    for (const j of [149n, 155n, 157n, 529n]) {
      const k = mod(1728n - j, p);
      const a4 = (((3n * j) % p) * k) % p;
      const a6 = ((((2n * j) % p) * k) % p) * k % p;
      expect(Fp_ellcard_SEA(a4, a6, p)).toBe(p + 1n);
    }
  }, 120000);
});

describe('Fp_ellcard_SEA against exhaustive point counting', () => {
  for (const p of [1009n, 2003n, 4001n]) {
    test(`p = ${p}`, () => {
      let n = 0;
      for (let a4 = 0n; a4 <= 11n; a4++)
        for (let a6 = 0n; a6 <= 11n; a6++) {
          if (mod(4n * a4 * a4 * a4 + 27n * a6 * a6, p) === 0n) continue;
          expect(Fp_ellcard_SEA(a4, a6, p)).toBe(bruteCard(a4, a6, p));
          n++;
        }
      expect(n).toBeGreaterThan(100);
    }, 600000);
  }
});

describe('Fp_ellcard_SEA against Shanks/Mestre (group.ts) and Schoof (advanced.ts)', () => {
  let st = 12345n;
  function rnd(m: bigint): bigint {
    st = (st * 6364136223846793005n + 1442695040888963407n) & ((1n << 64n) - 1n);
    return st % m;
  }
  for (const [p, count] of [
    [1000003n, 250],
    [1000000007n, 250],
    [4294967311n, 150],
    [281474976710677n, 40],
  ] as Array<[bigint, number]>) {
    test(`${count} random curves over F_${p}`, () => {
      let n = 0;
      for (let i = 0; i < count; i++) {
        const a4 = rnd(p);
        const a6 = rnd(p);
        if (mod(4n * a4 * a4 * a4 + 27n * a6 * a6, p) === 0n) continue;
        const N = Fp_ellcard_SEA(a4, a6, p);
        expect(N).toBe(ellcard({ a4, a6, p }));
        const t = p + 1n - N;
        expect(t * t <= 4n * p).toBe(true); /* Hasse */
        n++;
      }
      expect(n).toBeGreaterThan(count - 5);
    }, 600000);
  }

  test('agrees with Fp_ellcard_Schoof (independent implementation)', () => {
    for (const [p, a4, a6] of [
      [4294967311n, 1234567n, 7654321n],
      [281474976710677n, 3n, 5n],
    ] as Array<[bigint, bigint, bigint]>) {
      expect(Fp_ellcard_SEA(a4, a6, p)).toBe(Fp_ellcard_Schoof(a4, a6, p));
    }
  }, 600000);

  test('ellcard_sea (advanced.ts) delegates to SEA', () => {
    const E = { a4: 239810037n, a6: 543121245n, p: 1000000007n };
    expect(ellcard_sea(E)).toBe(1000047980n);
    expect(ellcard_sea(E)).toBe(ellcard({ ...E }));
  }, 120000);
});

describe('special branches', () => {
  test('j = 0 and j = 1728 are handled by the CM formulas', () => {
    for (const p of [1000003n, 1000000007n]) {
      expect(Fp_ellcard_SEA(0n, 5n, p)).toBe(bruteCardFast(0n, 5n, p));
      expect(Fp_ellcard_SEA(7n, 0n, p)).toBe(bruteCardFast(7n, 0n, p));
    }
    function bruteCardFast(a4: bigint, a6: bigint, p: bigint): bigint {
      return ellcard({ a4, a6, p });
    }
  }, 120000);

  test('the CM branch (non-squarefree modular equation) gives the right answer', () => {
    /* every curve over F_1009 whose modular equation is not squarefree for
     * some small ell (including a4=1,a6=10, whose CM branch reaches the
     * discriminant -288) */
    const p = 1009n;
    let cm = 0;
    for (let a4 = 1n; a4 <= 12n; a4++)
      for (let a6 = 1n; a6 <= 12n; a6++) {
        if (mod(4n * a4 * a4 * a4 + 27n * a6 * a6, p) === 0n) continue;
        const j = ellj(a4, a6, p);
        if (j === 0n || j === mod(1728n, p)) continue;
        let isCM = false;
        for (const ell of [3, 5, 7, 11, 13, 17, 19]) {
          const M = _internal.get_modular_eqn(ell);
          const mj = _internal.meqn_j(M, j, ell, p);
          if (_internal.study_modular_eqn(ell, mj, p).mt === _internal.ModType.MTcm)
            isCM = true;
        }
        if (!isCM) continue;
        expect(Fp_ellcard_SEA(a4, a6, p)).toBe(bruteCard(a4, a6, p));
        cm++;
      }
    expect(cm).toBeGreaterThan(15);
  }, 300000);

  test('smallfact aborts when a small prime divides the order', () => {
    const p = 1000000007n;
    for (let a4 = 1n; a4 <= 8n; a4++) {
      const a6 = 3n * a4 + 1n;
      const N = Fp_ellcard_SEA(a4, a6, p);
      /* smallfact = 1: abort as soon as any prime factor is found */
      const r = Fp_ellcard_SEA(a4, a6, p, 1);
      if (r === 0n) {
        /* upstream returns 0 only when it has *proved* a small factor */
        expect(N % 2n === 0n || smallestPrimeFactorBelow(N, 200) !== null).toBe(true);
      } else {
        expect(r).toBe(N);
      }
    }
    function smallestPrimeFactorBelow(n: bigint, B: number): bigint | null {
      for (let d = 2n; d <= BigInt(B); d++) if (n % d === 0n) return d;
      return null;
    }
  }, 300000);
});

describe("PARI's own ellsea regression vectors", () => {
  /* reference/pari/src/test/in/ellsea, expected values from
   * reference/pari/src/test/32/ellsea (lines 2-12 for `ellap`, 23-38 for
   * `ellsea(E, smallfact)`).  Copied verbatim. */
  test('ellap on the curves of the `do(i, v)` table', () => {
    const V: Array<[bigint, bigint, bigint, bigint]> = [
      [
        202600005603433095160409308644759862837n,
        25496852782325453225973142890909600552n,
        129550610797481291887769966647995045232n,
        -18627161351017007203n,
      ],
      [
        173327739907566197112155895875385467119n,
        52716988591102938437323369716512206005n,
        43087597392844950895070462564402654315n,
        18827282990304904850n,
      ],
      [
        523583591335747530615071369664554118036421993253n,
        25679429559575246581833628827363203226862930934n,
        78220741356817481602535950765825830003112603824n,
        -311256626765211726406998n,
      ],
      [
        1319450668936329467137913739322239157303860926441n,
        807652438980115949692649657326438677571309575087n,
        1017125626316888896817395440127041355136940446205n,
        -1156815323986765479761266n,
      ],
      [
        439581010348913995032270658729785287035964480323270935583n,
        323922016281172901245123590380881598241426088528431020005n,
        54496426275749371996644207660602248980615186517525561222n,
        8021839135157401454666601928n,
      ],
      [
        2979720374579183569554262247145622188470249961843364603751n,
        428654869348535206084607029945439317783967748844874233571n,
        1182279475380088064870625220629639405548336474256523329003n,
        69384671472347162238655401774n,
      ],
      [
        6243380271698146227966925307851825694742847655729810693741n,
        4068721281680536125235363885580194460678653324971583338307n,
        2519148351962491328666249705249360758373031631978108875818n,
        -28652256072001057705168347198n,
      ],
      [
        1606938044258990275541962092341162602522202993782792835304761n,
        1n,
        252199199707645577897249048746397012330572101453777389069968n,
        1271547588042840381566950172346n,
      ],
      [
        1267650600228229401496703205953n,
        1n,
        417990942431022911086532367249n,
        1854715558584444n,
      ],
      [590295810358705651741n, 1n, 3n, 20420247695n],
      [18446744073709551629n, 1n, 42n, -4742075250n],
    ];
    /* all 11 entries of upstream's table */
    expect(V.length).toBe(11);
    for (const [p, a4, a6, ap] of V) expect(p + 1n - Fp_ellcard_SEA(a4, a6, p)).toBe(ap);
  }, 2400000);

  test('ellsea(E, smallfact)', () => {
    const S: Array<[bigint, bigint, bigint, number, bigint]> = [
      [1048845330395786101209709n, 1n, 56n, 1, 0n],
      [1048845330395786101209709n, 1n, 56n, 2, 1048845330393999479019082n],
      [523n, 1n, 519n, 0, 486n],
      [
        1267650600228229401496703205653n,
        1122988618244467583984567614936n,
        429172847969450664478514342664n,
        -1,
        1267650600228229268303269105757n,
      ],
      [
        1267650600228229401496703205653n,
        1213812743793711191989251498394n,
        677975617584150034841507871840n,
        2,
        1267650600228228204181264322228n,
      ],
      [
        1267650600228229401496703205653n,
        1213812743793711191989251498394n,
        677975617584150034841507871840n,
        -2,
        0n,
      ],
      [1048845330395786101209839n, 1n, 15n, 1, 0n],
      [1048845330395786101209839n, 1n, 15n, -1, 0n],
      [1048845330395786101209839n, 1n, -15n, 1, 1048845330396436420140511n],
      [1048845330395786101209839n, 1n, -15n, -1, 0n],
      [1048845330395786101209839n, 1n, 15n, 7, 1048845330395135782279169n],
      [1048845330395786101209839n, 1n, 15n, -7, 1048845330395135782279169n],
      [1048845330395786101209839n, 1n, -15n, 7, 1048845330396436420140511n],
      [1048845330395786101209839n, 1n, -15n, -7, 1048845330396436420140511n],
    ];
    for (const [p, a4, a60, sf, want] of S)
      expect(Fp_ellcard_SEA(a4, mod(a60, p), p, sf)).toBe(want);
  }, 600000);
});

describe('cryptographic sizes', () => {
  test('128-bit prime, verified by the point-order proof', () => {
    const p = 340282366920938463463374607431768211297n;
    const a4 = 7932n;
    const a6 = 104734n;
    const N = Fp_ellcard_SEA(a4, a6, p);
    expect(N).toBe(340282366920938463430152381949486749310n);
    assertCardProven(a4, a6, p, N);
  }, 900000);

  test('NIST P-256 has its published cardinality', () => {
    const p = 0xffffffff00000001000000000000000000000000ffffffffffffffffffffffffn;
    const a4 = p - 3n;
    const a6 = 0x5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604bn;
    const expected =
      0xffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551n;
    expect(Fp_ellcard_SEA(a4, a6, p)).toBe(expected);
  }, 1800000);

  test('Curve25519 has cardinality 8 * l', () => {
    const p = (1n << 255n) - 19n;
    const A = 486662n;
    const a4 = (mod(3n - A * A, p) * Fp_inv(3n, p)) % p;
    const a6 = (mod(2n * A * A * A - 9n * A, p) * Fp_inv(27n, p)) % p;
    const expected =
      8n * ((1n << 252n) + 27742317777372353535851937790883648493n);
    expect(Fp_ellcard_SEA(a4, a6, p)).toBe(expected);
  }, 1800000);
});
