/**
 * Unit tests for sage/misc/randstate
 *
 * Reference: reference/sage/src/sage/misc/randstate.pyx
 */
import { describe, expect, test } from 'bun:test';
import { ValueError } from '../errors.js';
import {
  PythonRandom,
  RandState,
  SAGE_RAND_MAX,
  current_randstate,
  initial_seed,
  random,
  set_random_seed,
} from './randstate.js';

/**
 * Largest run of a repeating period-`p` pattern at the start of `seq`.
 */
function hasPeriod(seq: string[], p: number): boolean {
  if (seq.length <= p) return false;
  for (let i = p; i < seq.length; i++) {
    if (seq[i] !== seq[i - p]) return false;
  }
  return true;
}

describe('RandState', () => {
  describe('seed handling', () => {
    test('seed() returns the initial seed, not the evolving state', () => {
      // randstate.pyx:562-580
      const r = new RandState(314159);
      expect(r.seed()).toBe(314159n);
      r.c_rand_double();
      r.random_bits(64);
      r.randint(0n, 1000n);
      expect(r.seed()).toBe(314159n);
    });

    test('set_seed updates the reported seed and restarts the stream', () => {
      const r = new RandState(0);
      r.set_seed(12345);
      expect(r.seed()).toBe(12345n);
      const first = [r.random_bits(32), r.random_bits(32), r.random_bits(32)];
      r.set_seed(12345);
      const second = [r.random_bits(32), r.random_bits(32), r.random_bits(32)];
      expect(second).toEqual(first);
    });

    test('seed s and -s give the same stream', () => {
      // randstate.pyx:925-937 shows set_random_seed(-12345) and seed(12345)
      // producing the same value.  This is not because GMP takes an absolute
      // value: `randseed_mt` reduces the seed modulo 2^19937-20027 and adds 2
      // (randmts.c:126-130), which sends s and -s to +-(s+2) modulo the
      // powering modulus 2^19937-20023, and the exponent 1074888996 is even.
      const a = new RandState(12345);
      const b = new RandState(-12345);
      expect(b.random_bits(64)).toBe(a.random_bits(64));
    });

    test('accepts arbitrarily large bigint seeds', () => {
      const big = 305866218880103397618377824640007711767n;
      const r = new RandState(big);
      expect(r.seed()).toBe(big);
      const s = new RandState(big);
      expect(s.random_bits(64)).toBe(r.random_bits(64));
    });

    test('different seeds give different streams', () => {
      const a = new RandState(1).random_bits(64);
      const b = new RandState(2).random_bits(64);
      expect(a).not.toBe(b);
    });
  });

  describe('no short cycles (C1 regression)', () => {
    test('random_below(2) does not alternate', () => {
      const r = new RandState(7);
      const bits: string[] = [];
      for (let i = 0; i < 200; i++) bits.push(r.random_below(2n).toString());
      expect(hasPeriod(bits, 2)).toBe(false);
      expect(hasPeriod(bits, 1)).toBe(false);
      const ones = bits.filter((b) => b === '1').length;
      expect(ones).toBeGreaterThan(70);
      expect(ones).toBeLessThan(130);
    });

    test('random_below(4) has no period-4 cycle', () => {
      const r = new RandState(7);
      const vals: string[] = [];
      for (let i = 0; i < 200; i++) vals.push(r.random_below(4n).toString());
      for (let p = 1; p <= 8; p++) {
        expect(hasPeriod(vals, p)).toBe(false);
      }
    });

    test('randint(0, 9) has no period-10 cycle', () => {
      const r = new RandState(7);
      const vals: string[] = [];
      for (let i = 0; i < 500; i++) vals.push(r.randint(0n, 9n).toString());
      expect(hasPeriod(vals, 10)).toBe(false);
    });

    test('every low bit is balanced', () => {
      // For an LCG mod 2^64, bit k has period 2^(k+1); the low bits of the
      // Mersenne Twister are not structured this way.
      const r = new RandState(99);
      const n = 20000;
      const ones = new Array(8).fill(0);
      for (let i = 0; i < n; i++) {
        const v = r.random_bits(64);
        for (let k = 0; k < 8; k++) {
          if ((v >> BigInt(k)) & 1n) ones[k]!++;
        }
      }
      for (let k = 0; k < 8; k++) {
        expect(Math.abs(ones[k]! / n - 0.5)).toBeLessThan(0.02);
      }
    });

    test('random_below(10) is uniform (chi-squared)', () => {
      const r = new RandState(12345);
      const n = 200000;
      const counts = new Array(10).fill(0);
      for (let i = 0; i < n; i++) counts[Number(r.random_below(10n))]!++;
      const expected = n / 10;
      let chi2 = 0;
      for (const c of counts) chi2 += ((c - expected) * (c - expected)) / expected;
      // 99.9% quantile of chi^2 with 9 degrees of freedom is 27.88.
      expect(chi2).toBeLessThan(27.88);
    });
  });

  describe('random_bits', () => {
    test('returns values in [0, 2^bits)', () => {
      const r = new RandState(3);
      for (const bits of [1, 7, 31, 32, 33, 64, 65, 200]) {
        for (let i = 0; i < 50; i++) {
          const v = r.random_bits(bits);
          expect(v >= 0n).toBe(true);
          expect(v < 1n << BigInt(bits)).toBe(true);
        }
      }
    });

    test('random_bits(0) is 0', () => {
      expect(new RandState(3).random_bits(0)).toBe(0n);
    });

    test('rejects a negative bit count', () => {
      expect(() => new RandState(3).random_bits(-1)).toThrow(ValueError);
    });

    test('large requests use the full width', () => {
      const r = new RandState(5);
      let maxSeen = 0n;
      for (let i = 0; i < 200; i++) {
        const v = r.random_bits(128);
        if (v > maxSeen) maxSeen = v;
      }
      expect(maxSeen > 1n << 120n).toBe(true);
    });
  });

  describe('c_random and c_rand_double', () => {
    test('c_random returns a 31-bit integer', () => {
      const r = new RandState(1207);
      for (let i = 0; i < 1000; i++) {
        const v = r.c_random();
        expect(Number.isInteger(v)).toBe(true);
        expect(v).toBeGreaterThanOrEqual(0);
        expect(v).toBeLessThanOrEqual(SAGE_RAND_MAX);
      }
    });

    test('c_rand_double lies in [0, 1) and consumes two draws', () => {
      // randstate.pyx:883-895
      const r = new RandState(2718281828);
      for (let i = 0; i < 1000; i++) {
        const d = r.c_rand_double();
        expect(d).toBeGreaterThanOrEqual(0);
        expect(d).toBeLessThan(1);
      }

      const a = new RandState(11);
      const b = new RandState(11);
      a.c_rand_double();
      b.random_bits(25);
      b.random_bits(28);
      // Both states consumed the same amount of randomness.
      expect(a.random_bits(32)).toBe(b.random_bits(32));
    });

    test('c_rand_double is uniform on [0, 1)', () => {
      const r = new RandState(4242);
      const n = 100000;
      const buckets = new Array(10).fill(0);
      let sum = 0;
      for (let i = 0; i < n; i++) {
        const d = r.c_rand_double();
        sum += d;
        buckets[Math.min(9, Math.floor(d * 10))]!++;
      }
      expect(Math.abs(sum / n - 0.5)).toBeLessThan(0.01);
      for (const b of buckets) {
        expect(Math.abs(b / n - 0.1)).toBeLessThan(0.01);
      }
    });

    test('random() is an alias for c_rand_double', () => {
      const a = new RandState(77);
      const b = new RandState(77);
      expect(a.random()).toBe(b.c_rand_double());
    });
  });

  describe('random_below and randint', () => {
    test('random_below rejects a non-positive bound', () => {
      expect(() => new RandState(1).random_below(0n)).toThrow(ValueError);
      expect(() => new RandState(1).random_below(-5n)).toThrow(ValueError);
    });

    test('randint covers the whole inclusive range', () => {
      const r = new RandState(8);
      const seen = new Set<string>();
      for (let i = 0; i < 2000; i++) {
        const v = r.randint(-3n, 4n);
        expect(v >= -3n).toBe(true);
        expect(v <= 4n).toBe(true);
        seen.add(v.toString());
      }
      expect(seen.size).toBe(8);
    });

    test('randint rejects min > max', () => {
      expect(() => new RandState(1).randint(5n, 4n)).toThrow(ValueError);
    });
  });
});

describe('module-level functions', () => {
  test('set_random_seed makes the global stream reproducible', () => {
    set_random_seed(0);
    const first = [random(), random(), random()];
    set_random_seed(0);
    const second = [random(), random(), random()];
    expect(second).toEqual(first);
  });

  test('random() returns a 31-bit integer', () => {
    // randstate.pyx:1008-1020
    set_random_seed(31);
    for (let i = 0; i < 100; i++) {
      const v = random();
      expect(typeof v).toBe('bigint');
      expect(v >= 0n).toBe(true);
      expect(v < 1n << 31n).toBe(true);
    }
  });

  test('initial_seed returns the seed of the current randstate', () => {
    // randstate.pyx:1024-1042
    set_random_seed(42);
    expect(initial_seed()).toBe(42n);
    random();
    current_randstate().random_bits(64);
    expect(initial_seed()).toBe(42n);

    set_random_seed(5);
    expect(initial_seed()).toBe(5n);
  });
});

/**
 * Bit-for-bit equality with GMP and with SageMath.
 *
 * The `u32`/`b31`/`b64`/`urandomm` vectors were produced by calling the real
 * libgmp 6.3.0 (`gmp_randinit_default` + `gmp_randseed` + `gmp_urandomb_ui` /
 * `mpz_urandomm`) from C; the Sage values were produced by running SageMath
 * 10.3, and several are quoted from doctests in
 * `reference/sage/src/sage/misc/randstate.pyx`.
 */
describe('GMP parity', () => {
  test('gmp_urandomb_ui(state, 32) after gmp_randseed', () => {
    expect(Array.from({ length: 8 }, () => new RandState(42).random_bits(32))).toBeDefined();
    const r = new RandState(42);
    expect(Array.from({ length: 8 }, () => Number(r.random_bits(32)))).toEqual([
      1046334137, 3369841305, 4035822279, 2116953745, 1412997532, 1798563773, 1283018953,
      1273316211,
    ]);
    const s = new RandState(12345);
    expect(Array.from({ length: 8 }, () => Number(s.random_bits(32)))).toEqual([
      4239593463, 3179460376, 814094032, 2096529494, 2194027650, 2354695116, 618918688, 3640908660,
    ]);
    const big = new RandState(123456789012345678901234567890n);
    expect(Array.from({ length: 8 }, () => Number(big.random_bits(32)))).toEqual([
      3930944604, 1359561654, 3868168932, 3908873923, 1290190283, 1001594458, 2497613979, 649146301,
    ]);
  });

  test('the buffer is regenerated correctly across the 624-word boundary', () => {
    const r = new RandState(42);
    const all = Array.from({ length: 700 }, () => Number(r.random_bits(32)));
    expect(all.slice(620, 632)).toEqual([
      3828237464, 1393354641, 47043546, 63311296, 265568767, 1835891644, 580896834, 1363414749,
      2337062852, 755589586, 1565641461, 3290505273,
    ]);
  });

  test('partial-width requests consume one whole word each', () => {
    // GMP's randget_mt takes ceil(nbits/32) tempered outputs and masks.
    const b1 = new RandState(42);
    expect(Array.from({ length: 10 }, () => Number(b1.random_bits(1)))).toEqual([
      1, 1, 1, 1, 0, 1, 1, 1, 0, 1,
    ]);
    const b31 = new RandState(42);
    expect(Array.from({ length: 6 }, () => Number(b31.random_bits(31)))).toEqual([
      1046334137, 1222357657, 1888338631, 2116953745, 1412997532, 1798563773,
    ]);
    const b64 = new RandState(42);
    expect(Array.from({ length: 4 }, () => b64.random_bits(64))).toEqual([
      14473358198731295417n,
      9092247105955545799n,
      7724772586218365340n,
      5468851484994654409n,
    ]);
  });

  test('random_below is mpz_urandomm', () => {
    const r = new RandState(42);
    expect(Array.from({ length: 24 }, () => Number(r.random_below(11n)))).toEqual([
      9, 9, 7, 1, 9, 3, 2, 9, 4, 3, 6, 2, 10, 1, 4, 3, 3, 8, 7, 8, 7, 2, 10, 7,
    ]);
    const s = new RandState(42);
    expect(Array.from({ length: 6 }, () => Number(s.random_below(1000000n)))).toEqual([
      903865, 766617, 901831, 927377, 565660, 255933,
    ]);
    // A power of two consumes exactly log2(n) bits and never rejects.
    const t = new RandState(42);
    expect(Array.from({ length: 3 }, () => t.random_below(1n << 128n))).toEqual([
      167722355418488286110758738271573756601n,
      100882503720822822072470797230485840284n,
      328855610751865861188446067748149349250n,
    ]);
    // mpz_urandomm(rop, state, 1) returns 0 without drawing (urandomm.c:57-61).
    const u = new RandState(42);
    expect(u.random_below(1n)).toBe(0n);
    expect(Number(u.random_bits(32))).toBe(1046334137);
  });

  test("seed 0 keeps GMP's built-in default buffer", () => {
    // randstate.pyx:556-559: `if seed:` skips gmp_randseed for 0, so the
    // generator stays on randmt.c's `default_state` with mti = 2000 % 624.
    const r = new RandState(0);
    expect(Array.from({ length: 6 }, () => r.c_random())).toEqual([
      968665204, 1490187595, 43003183, 1984426947, 1745599944, 347717857,
    ]);
  });
});

describe('SageMath parity', () => {
  /** `ZZ.random_element(n)` (`integer_ring.pyx:798-812`). */
  const zzRandomElement = (r: RandState, n: bigint): bigint => {
    r.c_random(); // the unused `den` variable at integer_ring.pyx:801
    return r.random_below(n);
  };

  test('c_random and c_rand_double doctests', () => {
    // randstate.pyx:868-881 and :884-895
    expect(new RandState(1207).c_random()).toBe(2008037228);
    expect(new RandState(2718281828).c_rand_double()).toBe(0.22437207488974298);
  });

  test('ZZ.random_element(10^30) doctests', () => {
    // randstate.pyx:925-937
    const r = new RandState(-12345);
    expect(zzRandomElement(r, 10n ** 30n)).toBe(197130468050826967386035500824n);
    expect(zzRandomElement(r, 10n ** 30n)).toBe(601704412330400807050962541983n);
    const s = new RandState(12345);
    expect(zzRandomElement(s, 10n ** 30n)).toBe(197130468050826967386035500824n);
    // Verified against SageMath 10.3 with set_random_seed(0).
    const z = new RandState(0);
    expect([0, 1, 2].map(() => zzRandomElement(z, 10n ** 30n))).toEqual([
      670431516147804558529383265611n,
      772308321268490156498894882619n,
      551349305655019862415052218319n,
    ]);
  });

  test('python_random reproduces CPython seeded from the GMP state', () => {
    // randstate.pyx:562-580: set_random_seed(0); randstate(314159).python_random().random()
    set_random_seed(0);
    expect(new RandState(314159).python_random().random()).toBe(0.111439293741037);

    // randstate.pyx:603-610 (verified against SageMath 10.3)
    set_random_seed(5);
    const p = current_randstate().python_random();
    expect(p.random()).toBe(0.013558022446944151);
    expect(p.randrange(1000)).toBe(557n);

    set_random_seed(314159);
    expect(current_randstate().python_random().random()).toBe(0.29929142114291285);
  });

  test('python_random is cached per randstate', () => {
    set_random_seed(11);
    const a = current_randstate().python_random();
    const b = current_randstate().python_random();
    expect(a).toBe(b);
  });

  test('IntegerModRing(11).random_element() stream', () => {
    // integer_mod_ring.py:1543-1547 uses random.randint(0, order-1).
    set_random_seed(42);
    const p = current_randstate().python_random();
    expect(Array.from({ length: 24 }, () => Number(p.randint(0, 10)))).toEqual([
      5, 2, 3, 6, 3, 5, 1, 6, 2, 0, 5, 2, 6, 5, 8, 8, 6, 2, 9, 1, 6, 6, 2, 3,
    ]);
  });

  test('normalvariate matches sage.misc.prandom.normalvariate', () => {
    set_random_seed(42);
    const p = current_randstate().python_random();
    expect(Array.from({ length: 6 }, () => p.normalvariate(0, 1))).toEqual([
      -2.3229235471909804, 1.934962438857738, -0.3744876410111312, -2.377297348843676,
      0.11902611877180731, 1.1513070511415486,
    ]);
  });

  test('gen_lattice(m=10, seed=42) is reproducible from this stream', () => {
    // sage/crypto/lattice.py:238-243 builds the random block with
    // MatrixSpace(ZZ_q, m-n, n).random_element(), which for a matrix over
    // Z/qZ is `rstate.c_random() % q` per entry, row by row
    // (matrix_modn_dense_template.pxi randomize).
    const q = 11n;
    set_random_seed(42);
    const st = current_randstate();
    const A: bigint[][] = [];
    for (let i = 0; i < 6; i++) {
      A.push(Array.from({ length: 4 }, () => BigInt(st.c_random()) % q));
    }
    const abs = (x: bigint) => (x < 0n ? -x : x);
    const minrep = (a: bigint) => (abs(a - q) < abs(a) ? a - q : a);
    expect(A.map((r) => r.map((x) => Number(minrep(x))))).toEqual([
      [2, 4, 3, 5],
      [1, -5, -4, 2],
      [-4, 3, -1, 1],
      [-2, -3, -4, -1],
      [-5, -5, 3, 3],
      [-4, -3, 2, -5],
    ]);
  });
});

describe('PythonRandom', () => {
  test('matches CPython for an explicit seed', () => {
    // python3 -c "import random; r=random.Random(); r.seed(42); ..."
    const p = new PythonRandom();
    p.seed(42n);
    expect(p.random()).toBe(0.6394267984578837);
    expect(p.random()).toBe(0.025010755222666936);
    const q = new PythonRandom();
    q.seed(42n);
    expect(Array.from({ length: 5 }, () => Number(q.randint(1, 6)))).toEqual([6, 1, 1, 6, 3]);
    const t = new PythonRandom();
    t.seed(12345n);
    expect(t.getrandbits(64)).toBe(13515657874892102023n);
    expect(t.getrandbits(7)).toBe(1n);

    const u = new PythonRandom();
    u.seed(5n);
    expect(u.random()).toBe(0.6229016948897019);
    expect(u.randrange(1000)).toBe(759n);

    // getrandbits at every interesting width boundary.
    const v = new PythonRandom();
    v.seed(999n);
    expect([1, 3, 32, 33, 100].map((k) => v.getrandbits(k))).toEqual([
      1n,
      5n,
      343879192n,
      8145208016n,
      596881323810598836294328680313n,
    ]);

    const w = new PythonRandom();
    w.seed(7n);
    expect(Array.from({ length: 3 }, () => w.normalvariate(0, 1))).toEqual([
      -0.35590824951057143, 0.27915309343878736, 0.0970447726200016,
    ]);
  });

  test('getrandbits(0) is 0 and negative widths are rejected', () => {
    const p = new PythonRandom(1);
    expect(p.getrandbits(0)).toBe(0n);
    expect(() => p.getrandbits(-1)).toThrow(ValueError);
  });

  test('randrange rejects an empty range', () => {
    const p = new PythonRandom(1);
    expect(() => p.randrange(5, 5)).toThrow(ValueError);
    expect(() => p.randrange(0)).toThrow(ValueError);
  });
});
