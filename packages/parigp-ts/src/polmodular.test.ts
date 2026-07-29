/**
 * @module polmodular.test
 * @description Tests for modular polynomials `Phi_L(X, Y)` (PARI `polmodular.c`).
 *
 * ORACLES (all of them independent of this port):
 *
 * 1. **PARI's own regression suite.** `reference/pari/src/test/in/polmodular`
 *    hashes `polmodular(L, inv)` with a small DJB-style hash and compares
 *    against a table of golden values. We reimplement that hash verbatim
 *    (`hash_Z`/`hash_ZX`/`hash_ZXX` below, from the GP source of that file) and
 *    check the same golden values: the nine `modpoly_hashes` for `inv = 0` and
 *    `L = 2, 3, 5, 7, 11, 13, 17, 19, 23`, the `modfn_in` entries for
 *    `inv = INV_G2`, and the four `check_eval_modpoly` cases.
 *
 * 2. **Published closed forms.** `Phi_2` and `Phi_3` are checked coefficient by
 *    coefficient against the classical expressions (also reproduced in the
 *    comments of `polmodular.c:2039-2113`).
 *
 * 3. **Algebraic identities.** `Phi_L` is symmetric and has degree `L+1` in
 *    each variable; `Phi_L(j, j')` vanishes exactly on `L`-isogenous pairs, so
 *    `Phi_L(X, j)` splits into `L+1` roots over the CM primes used internally.
 *
 * 4. **PARI `polclass`** for the Hilbert class polynomials computed on the way.
 *
 * The golden values were produced with the PARI shipped inside the local
 * SageMath (2.15.4); the hashes in item 1 come from the vendored 2.18.1 test
 * file. They agree, as they must: `polmodular` is deterministic.
 */

import { describe, expect, it } from 'bun:test';
import { NotImplementedError } from './ifactor.js';
import { PariDomainError, type ZM } from './matkermod.js';
import {
  Flm_Fl_polmodular_evalx,
  Fp_polmodular_evalx,
  _internal as I,
  INV_F,
  INV_G2,
  INV_J,
  INV_W2W3,
  PariBugError,
  PariImplError,
  PariPriorityError,
  check_modinv,
  disc_best_modinv,
  modinv_good_disc,
  modinv_height_factor,
  modinv_level,
  modinv_max_internal_level,
  modinv_sparse_factor,
  polclass0,
  polmodular,
  polmodular_ZM,
  polmodular_ZXX,
  polmodular_db_add_level,
  polmodular_db_getp,
  polmodular_db_init,
  sympol_to_ZM,
} from './polmodular.js';

/* ================================================================== */
/* PARI's own hash from reference/pari/src/test/in/polmodular          */
/* ================================================================== */

const MASK64 = (1n << 64n) - 1n;
const INIT_H = 5381n;

/** GP `glue(h, a) = bitand(((h << 5) + h) + a, 2^64 - 1)` */
const glue = (h: bigint, a: bigint): bigint => ((h << 5n) + h + a) & MASK64;

/** GP `hash_Z(n)`: DJB2 over the base-256 limbs of `|n|`, then `sign(n) % 256`. */
function hash_Z(n: bigint): bigint {
  const sgn = n === 0n ? 0n : n > 0n ? 1n : 255n; /* GP's % is Euclidean */
  let m = n < 0n ? -n : n;
  let h = INIT_H;
  while (m !== 0n) {
    h = glue(h, m & 255n);
    m >>= 8n;
  }
  return glue(h, sgn);
}

/** GP `hash_ZX(pol)`; `v` is `Vec(pol)`, i.e. highest degree first. */
function hash_ZX(v: readonly bigint[]): bigint {
  let h = INIT_H;
  for (const c of v) h = glue(h, hash_Z(c));
  return h;
}

/** `Vec(f)` for a coefficient list `f` indexed by degree: strip and reverse. */
function toVec(f: readonly bigint[]): bigint[] {
  let d = -1;
  for (let i = f.length - 1; i >= 0; i--)
    if (f[i] !== 0n) {
      d = i;
      break;
    }
  if (d < 0) return [0n];
  const out: bigint[] = [];
  for (let i = d; i >= 0; i--) out.push(f[i]!);
  return out;
}

/** GP `hash_ZXX(pol)` for our column-major `Phi[j][i] = [X^i Y^j]`. */
function hash_ZXX(M: ZM, n: number): bigint {
  let h = INIT_H;
  for (let i = n - 1; i >= 0; i--) {
    const col: bigint[] = [];
    for (let j = 0; j < n; j++) col.push(M[j]![i]!);
    h = glue(h, hash_ZX(toVec(col)));
  }
  return h;
}

/* ================================================================== */

describe('class invariant helpers (polmodular.c:28-503)', () => {
  it('modinv_level matches upstream', () => {
    expect(modinv_level(INV_J)).toBe(1);
    expect(modinv_level(INV_G2)).toBe(3);
    expect(modinv_level(INV_F)).toBe(6);
    expect(modinv_level(INV_W2W3)).toBe(6);
    expect(modinv_level(131)).toBe(31); /* INV_ATKIN31 */
  });

  it('modinv_sparse_factor / height_factor / max_internal_level', () => {
    expect(modinv_sparse_factor(INV_J)).toBe(1);
    expect(modinv_sparse_factor(INV_G2)).toBe(3);
    expect(modinv_sparse_factor(INV_F)).toBe(24);
    expect(modinv_height_factor(INV_J)).toBe(1);
    expect(modinv_height_factor(INV_G2)).toBe(3);
    expect(modinv_height_factor(131)).toBe(16);
    expect(modinv_max_internal_level(INV_J)).toBe(5);
    expect(modinv_max_internal_level(INV_G2)).toBe(2);
  });

  it('modinv_good_disc / disc_best_modinv', () => {
    /* INV_J accepts everything; INV_G2 needs 3 not dividing D */
    expect(modinv_good_disc(INV_J, -71)).toBe(true);
    expect(modinv_good_disc(INV_G2, -71)).toBe(true);
    expect(modinv_good_disc(INV_G2, -3)).toBe(false);
    /* INV_F needs -D = 7 mod 8 and 3 not dividing D */
    expect(modinv_good_disc(INV_F, -71)).toBe(true);
    expect(modinv_good_disc(INV_F, -8)).toBe(false);
    expect(disc_best_modinv(-71)).toBe(INV_F);
  });

  it('check_modinv rejects invalid invariants (polclass.c:2111-2147)', () => {
    expect(() => check_modinv(7)).toThrow(PariDomainError);
    expect(() => check_modinv(INV_J)).not.toThrow();
    expect(() => check_modinv(131)).not.toThrow();
  });
});

describe('SMOOTH_INTS / HURWITZ_RATIO tables (polclass.c:1287-1416)', () => {
  /*
   * We regenerate the two 1200-entry tables with the GP code documented in
   * `polclass.c`; these are the first rows of the literal tables, copied
   * verbatim from the vendored source, as a transcription check.
   */
  it('SMOOTH_INTS agrees with the upstream literal', () => {
    const head = [
      0, 1, 2, 1, 4, 3, 8, 1, 2, 5, 16, 3, 32, 9, 6, 1, 64, 3, 128, 5, 10, 17, 256, 3, 4, 33, 2, 9,
      512, 7, 1024, 1, 18, 65, 12, 3, -1, 129, 34, 5, -1, 11, -1, 17, 6, 257, -1, 3, 8, 5, 66, 33,
      -1, 3, 20, 9, 130, 513, -1, 7,
    ];
    for (let v = 1; v <= head.length; v++) expect(I.SMOOTH_INTS[v]).toBe(head[v - 1]!);
    expect(I.SMOOTH_INTS[1200]).toBe(7); /* last entry of the upstream table */
  });

  it('HURWITZ_RATIO agrees with the upstream literal', () => {
    const head = [
      128, 384, 256, 384, 192, 768, 171, 384, 256, 576, 154, 768, 150, 512, 384, 384, 144, 768, 143,
      576, 342, 461, 140, 768, 192, 448, 256, 512, 138, 1152, 137, 384, 308, 432, 256, 768, 0, 427,
      299, 576, 0, 1024, 0, 461, 384, 419, 0, 768, 171, 576, 288, 448, 0, 768, 231, 512, 285, 412,
      0, 1152,
    ];
    for (let v = 1; v <= head.length; v++) expect(I.HURWITZ_RATIO[v]).toBe(head[v - 1]!);
    expect(I.HURWITZ_RATIO[1200]).toBe(1152);
  });
});

describe('Flx helpers used by the volcano walk', () => {
  const p = 1009n;

  it('Flv_roots_to_pol and Flx_roots are inverse', () => {
    const rts = [1n, 3n, 5n, 800n];
    const f = I.Flv_roots_to_pol(rts, p);
    expect(f[f.length - 1]).toBe(1n);
    expect(I.Flx_roots(f, p)).toEqual(rts);
    expect(I.Flx_nbroots(f, p)).toBe(4);
    /* X^2 + 1 is irreducible mod 1009? 1009 = 1 mod 4, so it has roots */
    expect(I.Flx_roots([1n, 0n, 1n], p).length).toBe(2);
    /* a genuinely rootless quadratic: X^2 - t for a non-residue t */
    let t = 2n;
    while (I.Flx_roots([(p - t) % p, 0n, 1n], p).length) t++;
    expect(I.Flx_oneroot([(p - t) % p, 0n, 1n], p)).toBeNull();
  });

  it('Flv_Flm_polint interpolates', () => {
    /* y = x^2 through (1,1), (2,4), (3,9) */
    const pols = I.Flv_Flm_polint([1n, 2n, 3n], [[1n, 4n, 9n]], p);
    expect(pols[0]).toEqual([0n, 0n, 1n]);
  });

  it('Flx_div_by_X_x performs synthetic division', () => {
    const f = I.Flv_roots_to_pol([2n, 7n], p);
    const { q, rem } = I.Flx_div_by_X_x(f, 2n, p);
    expect(rem).toBe(0n);
    expect(q).toEqual([(p - 7n) % p, 1n]);
  });
});

describe('class numbers (mftrace.c:2487-2508, quad.c, arith2.c)', () => {
  it('quadclassnos reproduces h(D) for fundamental D', () => {
    /* h(-3)=h(-4)=h(-7)=h(-8)=h(-11)=1, h(-15)=2, h(-23)=3, h(-47)=5,
     * h(-71)=7, h(-163)=1 */
    const table: Array<[number, number]> = [
      [-3, 1],
      [-4, 1],
      [-7, 1],
      [-8, 1],
      [-11, 1],
      [-15, 2],
      [-23, 3],
      [-47, 5],
      [-71, 7],
      [-163, 1],
      [-10007, 77],
    ];
    for (const [D, h] of table) expect(I.quadclassnos(D)).toBe(h);
  });

  it('quadnegclassnou splits D = u^2 D0 and returns h(D)', () => {
    /* -108 = 6^2 * (-3): h(-108) = 3 */
    const r = I.quadnegclassnou(-108);
    expect(r.D0).toBe(-3);
    expect(r.h).toBe(3);
    /* -104 is fundamental */
    expect(I.quadnegclassnou(-104).D0).toBe(-104);
    expect(I.quadnegclassnou(-104).h).toBe(6);
  });
});

describe('polclass0: Hilbert class polynomials (polclass.c:1980-2108)', () => {
  /* Golden values: PARI `Vecrev(polclass(D))`. */
  const cases: Array<[number, string]> = [
    [-7, '3375,1'],
    [-15, '-121287375,191025,1'],
    [-23, '12771880859375,-5151296875,3491750,1'],
    [
      -47,
      '16042929600623870849609375,-14982472850828613281250,5115161850595703125,' +
        '-9987963828125,2257834125,1',
    ],
    [
      -71,
      '737707086760731113357714241006081263,-425319473946139603274605151187659,' +
        '5138800366453976780323726329446,-823534263439730779968091389,' +
        '98394038810047812049302,-3091990138604570,313645809715,1',
    ],
    [
      -104,
      '65437179730333545242323676123103232,-25735039642229334200564710375424,' +
        '1378339984770204584193868955648,31013571054009020830449664,' +
        '739545196164376195072,-82028232174464,1',
    ],
    [
      -299,
      '-18273883965326272223717626628647422907813731016193733558272,' +
        '45797528808215150136248975363201860724351225694802411520,' +
        '-19207839443594488822936988943836177115227877227364352,' +
        '6417141278133218665289808655954275181523718111232,' +
        '-186547260770756829961971675685151791296544768,' +
        '2094055410006322146651491130721133658112,' +
        '-28635280874816126174326167699456,391086320728105978429440,1',
    ],
    [
      -1051,
      '8743043565409016736756337526364572620943665551429688944132020955512832,' +
        '609478936308852703138855519171409018549012233989963733448483930112,' +
        '84072685629795186930092590274005854385288302697479331815358464,' +
        '74010210743744060468842095314484372058081554805555200,' +
        '170566836806239391545422096319774885221433344,1',
    ],
  ];

  for (const [D, want] of cases) {
    it(`H_{${D}}(X) matches PARI polclass(${D})`, () => {
      const db = polmodular_db_init(INV_J);
      const H = polclass0(D, INV_J, db);
      expect(H.map(String).join(',')).toBe(want);
    }, 120000);
  }

  /*
   * Non-fundamental discriminants (u > 1): these exercise `enum_roots`'s
   * multi-generator walk, which the fundamental cases above never reach with a
   * cyclic factor of order 2.  Regression for `common_nbr`, which used to
   * report a double root of the degree-2 gcd as two distinct candidates and so
   * made `surface_parallel_path` fail for every j-invariant when n[0] == 2
   * (D = -288 never terminated).  Golden values: PARI `Vecrev(polclass(D))`.
   */
  const nonFundamental: Array<[number, string]> = [
    [-32, '12167000000,-52250000,1'],
    [-75, '5209253090426880,654403829760,1'],
    [-108, '-1879994705688000000000,224179462188000000,-151013228706000,1'],
    [
      -128,
      '-345363656226658026765625000000,-55499520947716391500000000,' +
        '-395258439243352250000,-2729960418308000,1',
    ],
    [-147, '11356800389480448000000,34848505552896000,1'],
    [
      -256,
      '-1064410681181869521037208505239142408,26925623396663008311375890966784,' +
        '-1826592673506207200904172752,-6761166974781862161312,1',
    ],
    [
      -288,
      '40994594700208456153393000000000000,136478143044657426076564000000000,' +
        '-87330008255955399131086000000,-142637765058468510772000,1',
    ],
    [
      -432,
      '42889619864187195342544128412237640625000000000000,' +
        '3869372376492639837782614434923625000000000000,' +
        '34904627315764077727184412247908187500000000,' +
        '1007059405271040783775694468925000000000,' +
        '280179539493990596285512318134750000,-22804995243537595825782822000,1',
    ],
  ];

  for (const [D, want] of nonFundamental) {
    it(`H_{${D}}(X) matches PARI polclass(${D}) (non-fundamental)`, () => {
      const db = polmodular_db_init(INV_J);
      const H = polclass0(D, INV_J, db);
      expect(H.map(String).join(',')).toBe(want);
    }, 120000);
  }

  it('small discriminants (polclass.c:1944-1957)', () => {
    const db = polmodular_db_init(INV_J);
    expect(polclass0(-3, INV_J, db)).toEqual([0n, 1n]);
    expect(polclass0(-4, INV_J, db)).toEqual([-1728n, 1n]);
    expect(polclass0(-4, INV_G2, db)).toEqual([-12n, 1n]);
  });
});

describe('Phi_2 and Phi_3 against their published closed forms', () => {
  /*
   * Phi_2(X, Y) = X^3 + Y^3 - X^2 Y^2
   *             + 1488 (X^2 Y + X Y^2) - 162000 (X^2 + Y^2)
   *             + 40773375 X Y + 8748000000 (X + Y) - 157464000000000
   * (`polmodular.c:2039-2050`)
   */
  it('Phi_2', () => {
    const M = polmodular_ZM(2);
    const c = (i: number, j: number) => M[j]![i]!;
    expect(c(3, 0)).toBe(1n);
    expect(c(0, 3)).toBe(1n);
    expect(c(2, 2)).toBe(-1n);
    expect(c(2, 1)).toBe(1488n);
    expect(c(1, 2)).toBe(1488n);
    expect(c(2, 0)).toBe(-162000n);
    expect(c(0, 2)).toBe(-162000n);
    expect(c(1, 1)).toBe(40773375n);
    expect(c(1, 0)).toBe(8748000000n);
    expect(c(0, 1)).toBe(8748000000n);
    expect(c(0, 0)).toBe(-157464000000000n);
    /* every other coefficient is zero */
    for (let i = 0; i < 4; i++)
      for (let j = 0; j < 4; j++) if (i + j > 3 && !(i === 2 && j === 2)) expect(c(i, j)).toBe(0n);
  });

  /*
   * Phi_3: [4,0,1], [3,3,-1], [3,2,2232], [3,1,-1069956], [3,0,36864000],
   * [2,2,2587918086], [2,1,8900222976000], [2,0,452984832000000],
   * [1,1,-770845966336000000], [1,0,1855425871872000000000], [0,0,0]
   * (`polmodular.c:2081-2098`)
   */
  it('Phi_3', () => {
    const M = polmodular_ZM(3);
    const c = (i: number, j: number) => M[j]![i]!;
    const want: Array<[number, number, bigint]> = [
      [4, 0, 1n],
      [3, 3, -1n],
      [3, 2, 2232n],
      [3, 1, -1069956n],
      [3, 0, 36864000n],
      [2, 2, 2587918086n],
      [2, 1, 8900222976000n],
      [2, 0, 452984832000000n],
      [1, 1, -770845966336000000n],
      [1, 0, 1855425871872000000000n],
      [0, 0, 0n],
    ];
    for (const [i, j, v] of want) {
      expect(c(i, j)).toBe(v);
      expect(c(j, i)).toBe(v);
    }
  });
});

describe('structural properties of Phi_L', () => {
  const levels = [2, 3, 5, 7, 11, 13];
  for (const L of levels) {
    it(`Phi_${L} is symmetric, monic and of degree ${L + 1} in each variable`, () => {
      const M = polmodular_ZM(L);
      expect(M.length).toBe(L + 2);
      for (const col of M) expect(col.length).toBe(L + 2);
      /* symmetry */
      for (let i = 0; i < L + 2; i++)
        for (let j = 0; j < L + 2; j++) expect(M[j]![i]).toBe(M[i]![j]!);
      /* degree L+1 in X, monic in X */
      expect(M[0]![L + 1]).toBe(1n);
      expect(M[L + 1]![0]).toBe(1n);
      /* the coefficient of X^{L+1} Y^k vanishes for k > 0 */
      for (let j = 1; j < L + 2; j++) expect(M[j]![L + 1]).toBe(0n);
      /* the leading coefficient of the "diagonal" is -1 for L > 1 */
      expect(M[L]![L]).toBe(-1n);
    }, 120000);
  }

  /*
   * Kronecker's congruence: for prime L,
   *   Phi_L(X, Y) = (X^L - Y)(X - Y^L)  (mod L)
   *              = X^{L+1} - X^L Y^L - X Y + Y^{L+1}  (mod L).
   * This pins every coefficient mod L and is completely independent of the
   * algorithm used to produce them.
   */
  for (const L of levels) {
    it(`Phi_${L} satisfies Kronecker's congruence mod ${L}`, () => {
      const M = polmodular_ZM(L);
      const Lb = BigInt(L);
      const want = (i: number, j: number): bigint => {
        if ((i === L + 1 && j === 0) || (i === 0 && j === L + 1)) return 1n;
        if ((i === L && j === L) || (i === 1 && j === 1)) return Lb - 1n;
        return 0n;
      };
      for (let i = 0; i < L + 2; i++)
        for (let j = 0; j < L + 2; j++) expect(((M[j]![i]! % Lb) + Lb) % Lb).toBe(want(i, j));
    }, 120000);
  }
});

describe("PARI's own regression hashes (test/in/polmodular)", () => {
  /* `modpoly_hashes` indexed by `lvl_idx`, for inv = 0. */
  const modpoly_hashes: Array<[number, bigint]> = [
    [2, 953115400354185n],
    [3, 619732354788530567n],
    [5, 7671381920119322245n],
    [7, 1662362517513198972n],
    [11, 11499552816775494464n],
    [13, 10945716853871337038n],
    [17, 1858790070632847848n],
    [19, 16279119036202003022n],
    [23, 9091292905489559584n],
  ];

  for (const [L, want] of modpoly_hashes) {
    it(`hash_ZXX(polmodular(${L})) = ${want}`, () => {
      expect(hash_ZXX(polmodular_ZM(L), L + 2)).toBe(want);
    }, 300000);
  }

  /* `modfn_in` entries for the class invariants we support (INV_G2 = 5). */
  const modfn_in: Array<[number, number, bigint]> = [
    [2, INV_G2, 7818678061185n],
    [5, INV_G2, 10135583858468178383n],
    [7, INV_G2, 9634555674574853739n],
  ];
  for (const [L, inv, want] of modfn_in) {
    it(`hash_ZXX(polmodular(${L}, ${inv})) = ${want}`, () => {
      expect(hash_ZXX(polmodular_ZM(L, inv), L + 2)).toBe(want);
    }, 300000);
  }

  /* `check_eval_modpoly(L, j, p, [h(Phi), h(Phi'), h(Phi''])` */
  const evalCases: Array<[number, bigint, bigint, [bigint, bigint, bigint]]> = [
    [5, 7n, 151n, [8033941431460000n, 243641761686181n, 243612090562303n]],
    [19, 7n, 151n, [11844895572672018496n, 369501438945078285n, 13082720985735388448n]],
    [5, 7n, 1099511627791n, [3901199766181530739n, 4054334766401667256n, 16751141247645108349n]],
    [23, 7n, 1099511627791n, [2360118342899681926n, 2787294817779511277n, 18359991236545579908n]],
  ];
  for (const [L, j, p, want] of evalCases) {
    it(`check_eval_modpoly(${L}, ${j}, ${p})`, () => {
      const r = Fp_polmodular_evalx(L, INV_J, j, p, true) as [bigint[], bigint[], bigint[]];
      expect(r.map((v) => hash_ZX(toVec(v)))).toEqual(want);
      /* without derivatives we must get the same polynomial back */
      const r0 = Fp_polmodular_evalx(L, INV_J, j, p, false) as bigint[];
      expect(r0).toEqual(r[0]);
    }, 300000);
  }
});

describe('evaluation entry points', () => {
  it('polmodular(5, , Mod(7, 151)) matches PARI', () => {
    expect(polmodular(5, INV_J, { j: 7n, p: 151n })).toEqual([
      125n,
      106n,
      58n,
      19n,
      144n,
      140n,
      1n,
    ]);
  });

  it('Flm_Fl_polmodular_evalx agrees with the generic matrix-vector product', () => {
    const p = 107n;
    for (const L of [2, 3, 5]) {
      const M = polmodular_ZM(L).map((c) => c.map((x) => ((x % p) + p) % p));
      for (const j of [0n, 1n, 33n, 57n, 106n]) {
        const n = L + 2;
        const direct = new Array<bigint>(n).fill(0n);
        let jp = 1n;
        for (let k = 0; k < n; k++) {
          for (let i = 0; i < n; i++) direct[i] = (direct[i]! + M[k]![i]! * jp) % p;
          jp = (jp * j) % p;
        }
        expect(Flm_Fl_polmodular_evalx(M, L, j, p)).toEqual(direct);
      }
    }
  });
});

describe('the polmodular database (polmodular.c:941-1029)', () => {
  it('caches levels and reduces mod p on demand', () => {
    const db = polmodular_db_init(INV_J);
    expect(db.inv).toBeNull();
    polmodular_db_add_level(db, 3, INV_J);
    expect(db.j[3]).not.toBeNull();
    const p = 97n;
    const red = polmodular_db_getp(db.j, 3, p);
    expect(red.length).toBe(5);
    for (const col of red) for (const x of col) expect(x >= 0n && x < p).toBe(true);
    expect(() => polmodular_db_getp(db.j, 7, p)).toThrow(PariBugError);
  });

  it('grows the table beyond its initial length', () => {
    const db = polmodular_db_init(INV_J);
    expect(db.j.length).toBe(33);
    /* level 2 is in the internal tables, so this is cheap */
    polmodular_db_add_level(db, 2, INV_J);
    expect(db.j[2]).not.toBeNull();
  });
});

describe('sympol_to_ZM (polmodular.c:1810-1821)', () => {
  it('desymmetrises a coefficient vector', () => {
    /* L = 1 is not a legal level, but the routine is purely combinatorial:
     * (L+1)(L+2)/2 = 3 coefficients for L = 1. */
    const M = sympol_to_ZM([10n, 20n, 30n], 1);
    expect(M.length).toBe(3);
    expect(M[0]![0]).toBe(10n);
    expect(M[0]![1]).toBe(20n);
    expect(M[1]![0]).toBe(20n);
    expect(M[1]![1]).toBe(30n);
    expect(M[0]![2]).toBe(1n);
    expect(M[2]![0]).toBe(1n);
  });
});

describe('error behaviour (test/in/polmodular argument checks)', () => {
  it('level 1 is a domain error', () => {
    expect(() => polmodular(1)).toThrow(PariDomainError);
  });
  it('composite level is not implemented', () => {
    expect(() => polmodular(6)).toThrow(PariImplError);
  });
  it('incompatible level/invariant pair is a domain error', () => {
    expect(() => polmodular_ZM(2, INV_F)).toThrow(PariDomainError);
    /* but the same invariant is fine at a level coprime to its own */
    expect(() => polmodular_ZM(5, INV_F)).not.toThrow();
  });
  it('bad invariant is a domain error', () => {
    expect(() => polmodular(19, 7)).toThrow(PariDomainError);
  });
  it('bad variable priority is a priority error', () => {
    expect(() => polmodular_ZXX(7, INV_J, 1, 0)).toThrow(PariPriorityError);
    expect(() => polmodular_ZXX(7, INV_J, 1, 1)).toThrow(PariPriorityError);
  });
  it('derivatives are a flag error when no argument is given', () => {
    expect(() => polmodular(7, INV_J, null, 1, true)).toThrow(/incorrect flag/);
  });
  it('unported class invariants throw NotImplementedError, naming what is missing', () => {
    /* Weber f at a level that needs the CM algorithm */
    expect(() => polmodular_ZM(7, INV_F)).toThrow(NotImplementedError);
    try {
      polmodular_ZM(7, INV_F);
    } catch (e) {
      expect((e as Error).message).toContain('SAGE_NOT_IMPLEMENTED');
      expect((e as Error).message).toContain('INV_J and INV_G2');
    }
  });
});
