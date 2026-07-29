/**
 * @module mpqs
 * @description Self-Initializing Multi-Polynomial Quadratic Sieve (SIMPQS),
 * ported from PARI/GP.
 *
 * Reference: reference/pari/src/basemath/mpqs.c (algorithm)
 *            reference/pari/src/basemath/mpqs.h (sizing tables, structures)
 *            reference/pari/src/basemath/F2v.c:397 (F2m_ker_sp)
 *            reference/pari/src/basemath/F2v.c:1063 (F2Ms_ker)
 *
 * Upstream's own sketch of the algorithm (mpqs.c:21-59):
 *
 *   Given an odd integer N > 1 to be factored, we throw in a small odd
 *   squarefree multiplier k so as to make kN = 1 mod 4 and to have many small
 *   primes over which X^2 - kN splits. We compute a factor base FB of such
 *   primes then look for values x0 such that Q0(x0) = x0^2 - kN can be
 *   decomposed over FB, up to a possible factor dividing k and a possible
 *   "large prime". Relations involving the latter can be combined into full
 *   relations which don't; full relations, by Gaussian elimination over F2 for
 *   the exponent vectors lead us to an expression X^2 - Y^2 divisible by N and
 *   hopefully to a nontrivial splitting when we compute gcd(X + Y, N).
 *
 * The Multi-Polynomial version sieves Q(x) = Q0(2Ax + B) = 4A(Ax^2 + Bx + C)
 * with A a product of omega_A factor base primes; "self-initializing" refers to
 * the fact that switching from one B to the next (there are 2^(omega_A-1) of
 * them for each A) is a cheap update of the sieve start positions.
 *
 * PORTING NOTES
 * -------------
 * - Upstream is *explicitly inexact* in exactly three places, and we reproduce
 *   its floating point there rather than replacing it: the Knuth-Schroeppel
 *   multiplier score (mpqs.c:221), the byte-sized scaled logarithms used by the
 *   sieve (mpqs.c:415), and the target size of A (mpqs.c:461). C `double` and
 *   JS `number` are both IEEE-754 binary64, and the `float`-typed `fbe_flogp`
 *   field is reproduced with a Float32Array. Everything else (the relations,
 *   the linear algebra, the gcds) is exact integer arithmetic.
 * - The sieve is a byte array with wrapping unsigned char arithmetic; a
 *   Uint8Array wraps identically.
 * - The unrolled sieve loops (mpqs_sieve_p/p1/p2, mpqs.c:914-964) are pure
 *   loop unrolling: they add `logp` to exactly the positions of the two
 *   arithmetic progressions which lie in the sieve array, in increasing order.
 *   We keep the plain loops.
 * - `mpqs_eval_sieve` scans the byte array `sizeof(mpqs_bit_array)` bytes at a
 *   time; we use the 8-byte (LONG_IS_64BIT, no SSE2) layout of mpqs.c:93.
 * - `Fl_sqrt` (arith1.c:847) picks *a* square root via a randomized generator
 *   search; which of the two roots comes back is not deterministic upstream
 *   either. We use deterministic Tonelli-Shanks. Both give a valid B.
 * - `F2Ms_ker` (F2v.c:1063) switches to block Lanczos above 640 rows; we always
 *   run the dense Gaussian elimination of `F2m_ker_sp` (F2v.c:397). Same kernel
 *   space, worse asymptotics; unreachable sizes are limited by relation
 *   collection long before the linear algebra becomes the bottleneck.
 * - `mpqs_class_init` / `mpqs_class_rels` (mpqs.c:1775-1865, the class group
 *   entry points used by buch2.c) are not exported: buch2.c is not ported. The
 *   MPQS_MODE_CLASSGROUP branches of the shared routines are transcribed
 *   nonetheless, so adding those entry points later is mechanical.
 */

import { Z_issquareall, isPrime, is_357_power } from './ifactor.js';

/* ------------------------------------------------------------------ */
/* mpqs.h: non-configurable sizing parameters                          */
/* ------------------------------------------------------------------ */

/** mpqs.h:27 - 'large primes' must be smaller than min(this, largest_FB_p) */
const MPQS_LP_BOUND = 12500000;
/** mpqs.h:33 - ~ -log2(0.9) */
const MPQS_A_FUDGE = 0.15;
/** mpqs.h:35 - max. this many candidates per polynomial */
const MPQS_CANDIDATE_ARRAY_SIZE = 2000;
/** mpqs.h:400 */
const MPQS_MAX_DIGIT_SIZE_KN = 107;
/** mpqs.h:108 - how many values for k we'll try */
const MPQS_POSSIBLE_MULTIPLIERS = 15;

/** mpqs.c:78-80 */
const REL_OFFSET = 20;
const REL_MASK = (1 << REL_OFFSET) - 1;

/** mpqs.h:86-92 - flag bits for fbe_flags */
const MPQS_FBE_CLEAR = 0x0;
const MPQS_FBE_DIVIDES_A = 0x1;
const MPQS_FBE_DIVIDES_N = 0x2;

/** mpqs.h:21-22 */
export const MPQS_MODE_FACTOR = 0;
export const MPQS_MODE_CLASSGROUP = 1;

/** mpqs.h:102-106 */
interface Multiplier {
  k: number;
  omega_k: number;
  kp: readonly number[];
}

/** mpqs.h:111-153 (cand_multipliers) */
const cand_multipliers: readonly Multiplier[] = [
  { k: 1, omega_k: 0, kp: [] },
  { k: 3, omega_k: 1, kp: [3] },
  { k: 5, omega_k: 1, kp: [5] },
  { k: 7, omega_k: 1, kp: [7] },
  { k: 11, omega_k: 1, kp: [11] },
  { k: 13, omega_k: 1, kp: [13] },
  { k: 15, omega_k: 2, kp: [3, 5] },
  { k: 17, omega_k: 1, kp: [17] },
  { k: 19, omega_k: 1, kp: [19] },
  { k: 21, omega_k: 2, kp: [3, 7] },
  { k: 23, omega_k: 1, kp: [23] },
  { k: 29, omega_k: 1, kp: [29] },
  { k: 31, omega_k: 1, kp: [31] },
  { k: 33, omega_k: 2, kp: [3, 11] },
  { k: 35, omega_k: 2, kp: [5, 7] },
  { k: 37, omega_k: 1, kp: [37] },
  { k: 39, omega_k: 2, kp: [3, 13] },
  { k: 41, omega_k: 1, kp: [41] },
  { k: 43, omega_k: 1, kp: [43] },
  { k: 47, omega_k: 1, kp: [47] },
  { k: 51, omega_k: 2, kp: [3, 17] },
  { k: 53, omega_k: 1, kp: [53] },
  { k: 55, omega_k: 2, kp: [5, 11] },
  { k: 57, omega_k: 2, kp: [3, 19] },
  { k: 59, omega_k: 1, kp: [59] },
  { k: 61, omega_k: 1, kp: [61] },
  { k: 65, omega_k: 2, kp: [5, 13] },
  { k: 67, omega_k: 1, kp: [67] },
  { k: 69, omega_k: 2, kp: [3, 23] },
  { k: 71, omega_k: 1, kp: [71] },
  { k: 73, omega_k: 1, kp: [73] },
  { k: 77, omega_k: 2, kp: [7, 11] },
  { k: 79, omega_k: 1, kp: [79] },
  { k: 83, omega_k: 1, kp: [83] },
  { k: 85, omega_k: 2, kp: [5, 17] },
  { k: 87, omega_k: 2, kp: [3, 29] },
  { k: 89, omega_k: 1, kp: [89] },
  { k: 91, omega_k: 2, kp: [7, 13] },
  { k: 93, omega_k: 2, kp: [3, 31] },
  { k: 95, omega_k: 2, kp: [5, 19] },
  { k: 97, omega_k: 1, kp: [97] },
];

/**
 * mpqs.h:292-398 (mpqs_parameters), indexed by the number of decimal digits of
 * kN, subscript 0 corresponding to 9 (or fewer) digits.
 * Columns: tolerance, lp_scale, M, size_of_FB, omega_A, pmin_index1.
 * `tolerance` is a C `float` in the table; Math.fround reproduces the value
 * that ends up in the (double) handle field.
 */
const mpqs_parameters: ReadonlyArray<readonly [number, number, number, number, number, number]> = (
  [
    /*  9 */ [0.8, 1, 350, 19, 3, 5],
    /* 10 */ [0.8, 1, 300, 23, 3, 5],
    /* 11 */ [0.8, 1, 1000, 27, 3, 5],
    /* 12 */ [0.8, 1, 1100, 27, 3, 5],
    /* 13 */ [0.8, 1, 1400, 31, 3, 5],
    /* 14 */ [0.8, 1, 2200, 33, 3, 5],
    /* 15 */ [0.8, 1, 2300, 39, 3, 5],
    /* 16 */ [0.8, 1, 2900, 43, 3, 5],
    /* 17 */ [0.8, 1, 3200, 51, 3, 5],
    /* 18 */ [0.8, 1, 2800, 55, 3, 5],
    /* 19 */ [0.8, 1, 3400, 65, 3, 5],
    /* 20 */ [0.8, 1, 3400, 71, 3, 5],
    /* 21 */ [0.8, 1, 5400, 90, 3, 5],
    /* 22 */ [0.8, 1, 5700, 95, 3, 5],
    /* 23 */ [0.8, 1, 5700, 110, 3, 5],
    /* 24 */ [0.8, 1, 6000, 130, 4, 7],
    /* 25 */ [0.8, 1, 6500, 140, 4, 7],
    /* 26 */ [0.9, 1, 9000, 160, 4, 7],
    /* 27 */ [1.12, 1, 10000, 160, 4, 7],
    /* 28 */ [1.17, 1, 13000, 180, 4, 11],
    /* 29 */ [1.22, 1, 14000, 220, 4, 11],
    /* 30 */ [1.3, 1, 13000, 240, 4, 11],
    /* 31 */ [1.33, 1, 11000, 240, 4, 13],
    /* 32 */ [1.36, 1, 14000, 300, 5, 13],
    /* 33 */ [1.4, 1, 15000, 340, 5, 13],
    /* 34 */ [1.43, 1, 15000, 380, 5, 17],
    /* 35 */ [1.48, 30, 15000, 380, 5, 17],
    /* 36 */ [1.53, 45, 16000, 440, 5, 17],
    /* 37 */ [1.6, 60, 15000, 420, 6, 19],
    /* 38 */ [1.66, 70, 15000, 520, 6, 19],
    /* 39 */ [1.69, 80, 16000, 540, 6, 23],
    /* 40 */ [1.69, 80, 16000, 600, 6, 23],
    /* 41 */ [1.69, 80, 16000, 700, 6, 23],
    /* 42 */ [1.69, 80, 24000, 900, 6, 29],
    /* 43 */ [1.69, 80, 26000, 1000, 6, 29],
    /* 44 */ [1.69, 80, 18000, 1100, 7, 31],
    /* 45 */ [1.69, 80, 20000, 1200, 7, 31],
    /* 46 */ [1.69, 80, 22000, 1300, 7, 37],
    /* 47 */ [1.69, 80, 24000, 1400, 7, 37],
    /* 48 */ [1.69, 80, 24000, 1600, 7, 37],
    /* 49 */ [1.72, 80, 28000, 1900, 7, 41],
    /* 50 */ [1.75, 80, 36000, 2100, 7, 41],
    /* 51 */ [1.8, 80, 32000, 2100, 7, 43],
    /* 52 */ [1.85, 80, 44000, 2300, 7, 43],
    /* 53 */ [1.9, 80, 44000, 2600, 7, 47],
    /* 54 */ [1.95, 80, 40000, 2700, 7, 47],
    /* 55 */ [1.95, 80, 48000, 3200, 7, 53],
    /* 56 */ [1.95, 80, 56000, 3400, 7, 53],
    /* 57 */ [2.0, 80, 40000, 3000, 8, 53],
    /* 58 */ [2.05, 80, 64000, 3400, 8, 59],
    /* 59 */ [2.1, 80, 64000, 3800, 8, 59],
    /* 60 */ [2.15, 80, 80000, 4300, 8, 61],
    /* 61 */ [2.2, 80, 80000, 4800, 8, 61],
    /* 62 */ [2.25, 80, 80000, 4600, 8, 67],
    /* 63 */ [2.39, 80, 80000, 4800, 8, 67],
    /* 64 */ [2.3, 80, 88000, 5400, 8, 67],
    /* 65 */ [2.31, 80, 120000, 6600, 8, 71],
    /* 66 */ [2.32, 80, 120000, 6800, 8, 71],
    /* 67 */ [2.33, 80, 144000, 7600, 8, 73],
    /* 68 */ [2.34, 80, 144000, 9000, 8, 73],
    /* 69 */ [2.35, 80, 160000, 9500, 8, 79],
    /* 70 */ [2.36, 80, 176000, 10500, 8, 79],
    /* 71 */ [2.37, 80, 240000, 11000, 9, 79],
    /* 72 */ [2.38, 80, 240000, 12500, 9, 83],
    /* 73 */ [2.41, 80, 240000, 13000, 9, 83],
    /* 74 */ [2.46, 80, 256000, 13250, 9, 83],
    /* 75 */ [2.51, 80, 256000, 14500, 9, 89],
    /* 76 */ [2.56, 80, 256000, 15250, 9, 89],
    /* 77 */ [2.58, 80, 320000, 17000, 9, 89],
    /* 78 */ [2.6, 80, 320000, 18000, 9, 89],
    /* 79 */ [2.63, 80, 320000, 19500, 9, 97],
    /* 80 */ [2.65, 80, 448000, 21000, 9, 97],
    /* 81 */ [2.72, 80, 448000, 22000, 9, 97],
    /* 82 */ [2.77, 80, 448000, 24000, 9, 101],
    /* 83 */ [2.82, 80, 480000, 23000, 10, 101],
    /* 84 */ [2.84, 80, 480000, 24000, 10, 103],
    /* 85 */ [2.86, 80, 512000, 28000, 10, 103],
    /* 86 */ [2.88, 80, 448000, 29000, 10, 107],
    /* 87 */ [2.9, 80, 512000, 32000, 10, 107],
    /* 88 */ [2.91, 80, 512000, 35000, 10, 109],
    /* 89 */ [2.92, 80, 512000, 38000, 10, 109],
    /* 90 */ [2.93, 80, 512000, 40000, 10, 113],
    /* 91 */ [2.94, 80, 770000, 32200, 10, 113],
    /* 92 */ [3.6, 90, 2000000, 35000, 9, 113],
    /* 93 */ [3.7, 90, 2000000, 37000, 9, 113],
    /* 94 */ [3.7, 90, 2000000, 39500, 9, 127],
    /* 95 */ [3.7, 90, 2500000, 41500, 9, 127],
    /* 96 */ [3.8, 90, 2500000, 45000, 10, 127],
    /* 97 */ [3.8, 90, 2500000, 47500, 10, 131],
    /* 98 */ [3.7, 90, 3000000, 51000, 10, 131],
    /* 99 */ [3.8, 90, 3000000, 53000, 10, 133],
    /*100 */ [3.8, 90, 875000, 50000, 10, 133],
    /*101 */ [3.8, 90, 3500000, 54000, 10, 139],
    /*102 */ [3.8, 90, 3500000, 57000, 10, 139],
    /*103 */ [3.9, 90, 4000000, 61000, 10, 139],
    /*104 */ [3.9, 90, 4000000, 66000, 10, 149],
    /*105 */ [3.9, 90, 4000000, 70000, 10, 149],
    /*106 */ [3.9, 90, 4000000, 75000, 10, 151],
    /*107 */ [3.9, 90, 4000000, 80000, 10, 151],
  ] as ReadonlyArray<readonly [number, number, number, number, number, number]>
).map(
  (r) =>
    [Math.fround(r[0]), r[1], r[2], r[3], r[4], r[5]] as readonly [
      number,
      number,
      number,
      number,
      number,
      number,
    ]
);

/* ------------------------------------------------------------------ */
/* small helpers (word-sized modular arithmetic, primes, gcd)          */
/* ------------------------------------------------------------------ */

/** gcd of two nonnegative bigints (PARI: gcdii). */
function gcdii(a: bigint, b: bigint): bigint {
  if (a < 0n) a = -a;
  if (b < 0n) b = -b;
  while (b) {
    const t = a % b;
    a = b;
    b = t;
  }
  return a;
}

/** number of decimal digits of |N| (mpqs.c:115 decimal_len = 1 + logint(N,10)) */
function decimal_len(N: bigint): number {
  return (N < 0n ? -N : N).toString().length;
}

/** largest e with B^e <= N, for N >= 1 (PARI: logint) */
function logint(N: bigint, B: bigint): number {
  let e = 0;
  let p = B;
  while (p <= N) {
    p *= B;
    e++;
  }
  return e;
}

function Fl_add(a: number, b: number, p: number): number {
  const s = a + b;
  return s >= p ? s - p : s;
}
function Fl_sub(a: number, b: number, p: number): number {
  const s = a - b;
  return s < 0 ? s + p : s;
}
/** valid for p < 2^26 (all FB primes are far below this) */
function Fl_mul(a: number, b: number, p: number): number {
  return (a * b) % p;
}
/** modular inverse mod p (p prime, a != 0 mod p) */
function Fl_inv(a: number, p: number): number {
  let t = 0;
  let newt = 1;
  let r = p;
  let newr = a % p;
  while (newr !== 0) {
    const q = Math.floor(r / newr);
    [t, newt] = [newt, t - q * newt];
    [r, newr] = [newr, r - q * newr];
  }
  if (r !== 1) throw new Error(`mpqs: Fl_inv: ${a} not invertible mod ${p}`);
  return t < 0 ? t + p : t;
}
function Fl_div(a: number, b: number, p: number): number {
  return Fl_mul(a, Fl_inv(b, p), p);
}
function Fl_powu(a: number, n: number, p: number): number {
  let r = 1;
  a %= p;
  while (n > 0) {
    if (n & 1) r = Fl_mul(r, a, p);
    a = Fl_mul(a, a, p);
    n >>>= 1;
  }
  return r;
}

/**
 * Square root of a mod p (p an odd prime, a a QR): deterministic
 * Tonelli-Shanks. PARI's Fl_sqrt (arith1.c:847) searches for a generator of the
 * 2-Sylow at random, so it is not deterministic either; both roots are valid.
 */
function Fl_sqrt(a: number, p: number): number {
  a %= p;
  if (a === 0) return 0;
  if (p === 2) return a;
  if (p % 4 === 3) return Fl_powu(a, (p + 1) / 4, p);
  let q = p - 1;
  let e = 0;
  while ((q & 1) === 0) {
    q >>= 1;
    e++;
  }
  /* smallest non-residue */
  let n = 2;
  while (Fl_powu(n, (p - 1) / 2, p) !== p - 1) n++;
  let y = Fl_powu(n, q, p);
  let r = e;
  let x = Fl_powu(a, (q - 1) / 2, p);
  let b = Fl_mul(Fl_mul(a, x, p), x, p);
  x = Fl_mul(a, x, p);
  while (b !== 1) {
    let m = 0;
    let t = b;
    while (t !== 1) {
      t = Fl_mul(t, t, p);
      m++;
    }
    let tt = y;
    for (let i = 0; i < r - m - 1; i++) tt = Fl_mul(tt, tt, p);
    y = Fl_mul(tt, tt, p);
    r = m;
    x = Fl_mul(x, tt, p);
    b = Fl_mul(b, y, p);
  }
  return x;
}

/** Kronecker symbol (x|y) for x, y >= 0 (PARI: krouu) */
function krouu(x: number, y: number): number {
  if (y === 0) return x === 1 ? 1 : 0;
  let s = 1;
  if ((y & 1) === 0) {
    if ((x & 1) === 0) return 0;
    let v = 0;
    while ((y & 1) === 0) {
      y >>= 1;
      v++;
    }
    if (v & 1) {
      const m = x & 7;
      if (m === 3 || m === 5) s = -s;
    }
  }
  /* Jacobi symbol (x|y), y odd > 0 */
  x %= y;
  while (x !== 0) {
    while ((x & 1) === 0) {
      x >>= 1;
      const m = y & 7;
      if (m === 3 || m === 5) s = -s;
    }
    const t = x;
    x = y;
    y = t;
    if ((x & 3) === 3 && (y & 3) === 3) s = -s;
    x %= y;
  }
  return y === 1 ? s : 0;
}

/** Kronecker symbol (N|p) for a bigint N and a word p (PARI: kroiu) */
function kroiu(N: bigint, p: number): number {
  if (p === 2) {
    if ((N & 1n) === 0n) return 0;
    const m = Number(N & 7n);
    return m === 1 || m === 7 ? 1 : -1;
  }
  return krouu(Number(N % BigInt(p)), p);
}

/**
 * x mod p as a nonnegative machine integer (PARI: umodiu, kernel/none/mp.c:500
 * returns `p - r` for negative x, not JS's truncated remainder).
 */
function umodiu(x: bigint, p: number): number {
  const r = Number(x % BigInt(p));
  return r < 0 ? r + p : r;
}

/** number of trailing zero bits of a nonzero uint32 (PARI: vals) */
function vals(x: number): number {
  return 31 - Math.clz32(x & -x);
}

/** simple sieve of Eratosthenes, all primes <= limit */
function simpleSieve(limit: number): number[] {
  const s = new Uint8Array(limit + 1);
  const out: number[] = [];
  for (let i = 2; i <= limit; i++) {
    if (!s[i]) {
      out.push(i);
      if (i <= limit / i) for (let j = i * i; j <= limit; j += i) s[j] = 1;
    }
  }
  return out;
}

/** unbounded prime iterator starting at `start` (PARI: u_forprime_init/next) */
function* u_forprime(start: number): Generator<number> {
  const SEG = 1 << 16;
  let lo = Math.max(2, start);
  let base: number[] = simpleSieve(1024);
  let baseLimit = 1024;
  for (;;) {
    const hi = lo + SEG - 1;
    const need = Math.floor(Math.sqrt(hi)) + 1;
    if (need > baseLimit) {
      baseLimit = Math.max(need, baseLimit * 2);
      base = simpleSieve(baseLimit);
    }
    const seg = new Uint8Array(SEG);
    for (const p of base) {
      if (p > hi / p) break;
      let s = Math.ceil(lo / p) * p;
      if (s < p * p) s = p * p;
      for (let j = s; j <= hi; j += p) seg[j - lo] = 1;
    }
    for (let i = lo; i <= hi; i++) if (!seg[i - lo] && i >= 2) yield i;
    lo = hi + 1;
  }
}

/* ------------------------------------------------------------------ */
/* the handle (mpqs.h:199-248)                                         */
/* ------------------------------------------------------------------ */

/**
 * Factor base, mpqs.h:65-83. One "struct of arrays" instead of an array of
 * 32-byte structs; the cache-line alignment of mpqs_FB_ctor (mpqs.c:164) is a
 * pure memory-layout optimization.
 *
 * Layout (mpqs.c:287-308): FB[0] unused, FB[1] stands for -1, FB[2] is 2,
 * FB[3..index0_FB-1] are the prime factors of k, the real odd FB primes start
 * at index0_FB, FB[size_of_FB+1] is the largest prime and FB[size_of_FB+2] is a
 * sentinel with p = 0.
 */
interface FBArrays {
  p: Int32Array;
  start1: Int32Array;
  start2: Int32Array;
  sqrt_kN: Int32Array;
  flogp: Float32Array;
  logval: Uint8Array;
  flags: Uint8Array;
}

interface Handle {
  N: bigint;
  kN: bigint;
  A: bigint;
  B: bigint;
  FB: FBArrays;
  /** (1/A) H[i] mod p_j, flat: index j * (omega_A - 1) + i (mpqs.h:190-194) */
  inv_A_H: Int32Array;
  /** per_A_pr[i]._H (mpqs.h:159-162) */
  per_A_H: bigint[];
  /** per_A_pr[i]._i */
  per_A_i: Int32Array;
  sieve_array: Uint8Array;
  sieve_words: Uint32Array;
  candidates: Int32Array;
  relaprimes: Int32Array;

  M: number;
  size_of_FB: number;
  index0_FB: number;
  index1_FB: number;
  index2_FB: number;
  index2_moved: number;
  sieve_threshold: number;
  two_is_norm: number;
  two_is_bad: number;
  omega_A: number;
  no_B: number;
  l2_target_A: number;
  bin_index: number;
  index_i: number;
  index_j: number;
  target_rels: number;
  largest_FB_p: number;
  pmin_index1: number;
  lp_scale: number;
  lp_bound: number;
  digit_size_kN: number;
  _k: Multiplier;
  tolerance: number;
  dkN: number;
  l2sqrtkN: number;
  l2M: number;
  /** not in the C handle: turns on upstream's -DMPQS_DEBUG relation checks */
  debug: boolean;
}

/* accessors of mpqs.h:168-173 */
const MPQS_I = (h: Handle, i: number) => h.per_A_i[i];
const MPQS_AP = (h: Handle, i: number) => h.FB.p[h.per_A_i[i]];
const MPQS_LP = (h: Handle, i: number) => h.FB.flogp[h.per_A_i[i]];
const MPQS_SQRT = (h: Handle, i: number) => h.FB.sqrt_kN[h.per_A_i[i]];

/* ------------------------------------------------------------------ */
/* INITIAL SIZING (mpqs.c:111-152)                                     */
/* ------------------------------------------------------------------ */

/**
 * mpqs.c:122 (mpqs_set_parameters). To be called after choosing k and putting
 * kN into the handle. Return false when kN is too large, true when we're ok.
 */
function mpqs_set_parameters(h: Handle): boolean {
  const D = decimal_len(h.kN < 0n ? -h.kN : h.kN);
  h.digit_size_kN = D;
  if (D > MPQS_MAX_DIGIT_SIZE_KN) return false;
  const P = mpqs_parameters[Math.max(0, D - 9)];
  h.tolerance = P[0];
  h.lp_scale = P[1];
  /* make room for prime factors of k if any */
  const s = P[3] + h._k.omega_k;
  h.size_of_FB = s;
  /* prime factors of k behave like real FB primes at the Gauss stage */
  h.target_rels = s >= 200 ? s + 10 : Math.trunc(s * 1.05);
  h.M = P[2];
  h.omega_A = P[4];
  h.no_B = 1 << (P[4] - 1);
  h.pmin_index1 = P[5];
  h.index0_FB = 3 + h._k.omega_k;
  return true;
}

/* ------------------------------------------------------------------ */
/* OBJECT HOUSEKEEPING (mpqs.c:154-208)                                */
/* ------------------------------------------------------------------ */

/** mpqs.c:164 (mpqs_FB_ctor) */
function mpqs_FB_ctor(h: Handle): FBArrays {
  const n = h.size_of_FB + 3; /* slots 0, 1 and a sentinel at the end */
  const FB: FBArrays = {
    p: new Int32Array(n),
    start1: new Int32Array(n),
    start2: new Int32Array(n),
    sqrt_kN: new Int32Array(n),
    flogp: new Float32Array(n),
    logval: new Uint8Array(n),
    flags: new Uint8Array(n),
  };
  h.FB = FB;
  h.inv_A_H = new Int32Array((h.size_of_FB + 3) * Math.max(1, h.omega_A - 1));
  return FB;
}

/** mpqs.c:178 (mpqs_sieve_array_ctor) */
function mpqs_sieve_array_ctor(h: Handle): void {
  const size = (h.M << 1) + 1;
  /* round up to a multiple of 8 so that the word scan of mpqs_eval_sieve
   * (which reads sizeof(mpqs_bit_array) bytes at a time) stays in bounds */
  const alloc = (size + 7) & ~7;
  const buf = new ArrayBuffer(alloc);
  h.sieve_array = new Uint8Array(buf);
  h.sieve_words = new Uint32Array(buf);
  h.sieve_array[h.M << 1] = 255; /* sentinel, mpqs.c:186 */
  h.candidates = new Int32Array(MPQS_CANDIDATE_ARRAY_SIZE + 16);
  /* mpqs.c:190 caps this at MAX_PE_PAIR = 60 pairs, which a candidate with
   * more than 60 distinct FB divisors would overrun; we size it so that it
   * cannot overrun. Same behaviour, no buffer overflow. */
  h.relaprimes = new Int32Array((h.size_of_FB + 2) << 1);
}

/** mpqs.c:196 (mpqs_poly_ctor) */
function mpqs_poly_ctor(h: Handle): void {
  h.per_A_H = new Array<bigint>(h.omega_A).fill(0n);
  h.per_A_i = new Int32Array(h.omega_A);
  h.A = 0n;
  h.B = 0n;
}

/* ------------------------------------------------------------------ */
/* FACTOR BASE SETUP (mpqs.c:210-383)                                  */
/* ------------------------------------------------------------------ */

/**
 * mpqs.c:221 (mpqs_find_k). Fill in the best-guess multiplier k for N, forcing
 * kN = 1 mod 4 (Knuth-Schroeppel, see Silverman, Math. Comp. 48 (1987)).
 * Returns a prime factor of N if one is found, else 0.
 */
function mpqs_find_k(h: Handle): number {
  const N_mod_8 = Number(h.N & 7n);
  const N_mod_4 = N_mod_8 & 3;
  const dl = decimal_len(h.N);
  const D = Math.max(0, Math.min(dl, MPQS_MAX_DIGIT_SIZE_KN) - 9);
  const MPQS_MULTIPLIER_SEARCH_DEPTH = mpqs_parameters[D][3];
  const MPQS_NB_MULTIPLIERS = dl < 40 ? 5 : MPQS_POSSIBLE_MULTIPLIERS;

  const cacheK: Multiplier[] = [];
  const cacheNp: number[] = [];
  const cacheVal: number[] = [];
  let nbk = 0;
  for (let i = 0; i < cand_multipliers.length; i++) {
    const cand_k = cand_multipliers[i];
    const k = cand_k.k;
    if ((k & 3) !== N_mod_4) continue; /* want kN = 1 (mod 4) */
    let v = -Math.log(k) / 2;
    if ((k & 7) === N_mod_8) v += Math.LN2; /* kN = 1 (mod 8) */
    cacheNp[nbk] = 0;
    cacheK[nbk] = cand_k;
    cacheVal[nbk] = v;
    if (++nbk === MPQS_NB_MULTIPLIERS) break;
  }
  if (nbk > MPQS_POSSIBLE_MULTIPLIERS) nbk = MPQS_POSSIBLE_MULTIPLIERS;

  for (const p of u_forprime(2)) {
    const kroNp = kroiu(h.N, p);
    let seen = 0;
    if (!kroNp) return p;
    for (let i = 0; i < nbk; i++) {
      if (cacheNp[i] > MPQS_MULTIPLIER_SEARCH_DEPTH) continue;
      seen++;
      const krokp = krouu(cacheK[i].k % p, p);
      if (krokp === kroNp) {
        cacheVal[i] += (2 * Math.log(p)) / p;
        cacheNp[i]++;
      } else if (krokp === 0) {
        cacheVal[i] += Math.log(p) / p;
        cacheNp[i]++;
      }
    }
    if (!seen) break; /* gone through SEARCH_DEPTH primes for all k */
  }
  let best_i = 0;
  let v = cacheVal[0];
  for (let i = 1; i < nbk; i++) {
    if (cacheVal[i] > v) {
      best_i = i;
      v = cacheVal[i];
    }
  }
  h._k = cacheK[best_i];
  return 0;
}

/**
 * mpqs.c:311 (mpqs_create_FB). Create a factor base of `size_of_FB` primes p_i
 * with legendre(kN, p_i) != -1. If a prime factor of N is found during the
 * construction it is returned, else 0.
 */
function mpqs_create_FB(h: Handle, wantFactor: boolean): number {
  const FB = mpqs_FB_ctor(h);
  const size = h.size_of_FB;
  const k = h._k.k;

  h.largest_FB_p = 0;
  FB.p[1] = -1;
  FB.p[2] = 2;
  FB.flags[2] = MPQS_FBE_CLEAR;
  let i: number;
  for (i = 3; i < h.index0_FB; i++) {
    /* executes omega_k = 0, 1 or 2 times */
    const kp = h._k.kp[i - 3];
    FB.p[i] = kp;
    FB.flags[i] = MPQS_FBE_CLEAR;
    FB.flogp[i] = Math.fround(Math.log2(kp));
    FB.sqrt_kN[i] = 0;
  }
  for (const p of u_forprime(3)) {
    if (i >= size + 2) break;
    if (p > k || k % p) {
      const kNp = umodiu(h.kN, p);
      const kr = krouu(kNp, p);
      if (kr >= 0) {
        FB.flags[i] = MPQS_FBE_CLEAR;
        if (kr === 0) {
          if (wantFactor) return p;
          /* classgroup mode: skip primes whose square divides kN */
          let t = h.kN;
          const P = BigInt(p);
          let val = 0;
          while (t % P === 0n) {
            t /= P;
            val++;
          }
          if (val > 1) continue;
          FB.flags[i] = MPQS_FBE_DIVIDES_N;
        }
        FB.p[i] = p;
        FB.flogp[i] = Math.fround(Math.log2(p));
        /* x such that x^2 = kN (mod p_i) */
        FB.sqrt_kN[i] = Fl_sqrt(kNp, p);
        i++;
      }
    }
  }

  FB.p[i] = 0; /* sentinel */
  h.largest_FB_p = FB.p[i - 1]; /* at subscript size_of_FB + 1 */

  /* locate the smallest prime that will be used for sieving */
  for (i = h.index0_FB; FB.p[i] !== 0; i++) if (FB.p[i] >= h.pmin_index1) break;
  h.index1_FB = i;
  return 0;
}

/* ------------------------------------------------------------------ */
/* MISC HELPER FUNCTIONS (mpqs.c:385-493)                              */
/* ------------------------------------------------------------------ */

/**
 * mpqs.c:415 (mpqs_set_sieve_threshold). Rescale log2(sqrt(kN)*M/lp^tolerance)
 * to 232 and fill in the byte-sized approximate scaled logarithms of the p_i.
 */
function mpqs_set_sieve_threshold(h: Handle): void {
  const FB = h.FB;
  h.l2sqrtkN = 0.5 * Math.log2(h.dkN);
  h.l2M = Math.log2(h.M);
  const log_maxval = h.l2sqrtkN + h.l2M - MPQS_A_FUDGE;
  const log_multiplier = 232.0 / log_maxval;
  /* (unsigned char) cast truncates toward 0 and wraps mod 256 */
  h.sieve_threshold =
    ((Math.trunc(log_multiplier * (log_maxval - h.tolerance * Math.log2(h.largest_FB_p))) & 255) +
      1) &
    255;
  if (h.sieve_threshold < 128) h.sieve_threshold = 128;
  for (let i = h.index0_FB; i < h.size_of_FB + 2; i++) {
    FB.logval[i] = Math.trunc(log_multiplier * FB.flogp[i]) & 255;
  }
}

/**
 * mpqs.c:461 (mpqs_locate_A_range). Find the optimum place in the FB to pick
 * the prime factors of A. Return true on success.
 */
function mpqs_locate_A_range(h: Handle): boolean {
  let i = h.index0_FB + 2 * h.omega_A - 4;
  const FB = h.FB;
  h.l2_target_A = h.l2sqrtkN - h.l2M - MPQS_A_FUDGE;
  const l2_target_pA = h.l2_target_A / h.omega_A;
  while (FB.p[i] && FB.flogp[i] <= l2_target_pA) i++;
  if (i > h.size_of_FB - 3) return false; /* should never happen */
  h.index2_FB = i - 1;
  return true;
}

/* ------------------------------------------------------------------ */
/* SELF-INITIALIZATION (mpqs.c:495-908)                                */
/* ------------------------------------------------------------------ */

/**
 * mpqs.c:528 (mpqs_increment). Increment x to the next larger value with the
 * same number of 1 bits, avoiding values obtained by moving a single 1 bit one
 * position to the left.
 */
function mpqs_increment(x: number): number {
  let r1_mask: number;
  let r01_mask: number;
  let slider = 1;
  let post = false;
  switch (x & 0x1f) {
    case 29:
      x = (x + 1) >>> 0;
      post = true;
      break;
    case 26:
      x = (x + 2) >>> 0;
      post = true;
      break;
    case 1:
    case 3:
    case 6:
    case 9:
    case 11:
    case 17:
    case 19:
    case 22:
    case 25:
    case 27:
      return (x + 3) >>> 0;
    case 20:
      x = (x + 4) >>> 0;
      post = true;
      break;
    case 5:
    case 12:
    case 14:
    case 21:
      return (x + 5) >>> 0;
    case 2:
    case 7:
    case 13:
    case 18:
    case 23:
      return (x + 6) >>> 0;
    case 10:
      return (x + 7) >>> 0;
    case 8:
      x = (x + 8) >>> 0;
      post = true;
      break;
    case 4:
    case 15:
      return (x + 12) >>> 0;
    default: {
      /* 0, 16, 24, 28, 30, 31 */
      r1_mask = ((((x ^ (x - 1)) >>> 0) + 1) >>> 1) >>> 0;
      r01_mask = ((((x ^ (x + r1_mask)) >>> 0) + r1_mask) >>> 2) >>> 0;
      if (r1_mask === r01_mask) {
        x = (x + r1_mask) >>> 0;
        post = true;
        break;
      }
      if (r1_mask === 1) {
        x = (x + r01_mask) >>> 0;
        post = true;
        break;
      }
      if (r1_mask === 2) return (x + (r01_mask >>> 1) + 1) >>> 0;
      while (r01_mask > r1_mask && slider < r1_mask) {
        r01_mask >>>= 1;
        slider <<= 1;
      }
      return (x + r01_mask + slider - 1) >>> 0;
    }
  }
  if (!post) return x >>> 0;
  /* post-process cases which couldn't be finalized above */
  r1_mask = ((((x ^ (x - 1)) >>> 0) + 1) >>> 1) >>> 0;
  r01_mask = ((((x ^ (x + r1_mask)) >>> 0) + r1_mask) >>> 2) >>> 0;
  if (r1_mask === r01_mask) return (x + r1_mask) >>> 0;
  if (r1_mask === 1) return (x + r01_mask) >>> 0;
  if (r1_mask === 2) return (x + (r01_mask >>> 1) + 1) >>> 0;
  slider = 1;
  while (r01_mask > r1_mask && slider < r1_mask) {
    r01_mask >>>= 1;
    slider <<= 1;
  }
  return (x + r01_mask + slider - 1) >>> 0;
}

/**
 * mpqs.c:599 (mpqs_si_choose_primes). Self-init (1): advance the bit pattern
 * and choose the primes for A. Returns 1 when all is fine, 0 when the caller
 * should retry (index2_FB was bumped, or a chosen prime divides N).
 */
function mpqs_si_choose_primes(h: Handle, missing_primes: number[] | null): number {
  const FB = h.FB;
  let l2_last_p = h.l2_target_A;
  const omega_A = h.omega_A;
  let i: number;
  let j: number;
  let v2: number;
  let prev_last_p_idx = 0;
  let room = h.index2_FB - h.index0_FB - omega_A + 4;
  let p: number;

  if (h.bin_index === 0) {
    h.bin_index = ((1 << (omega_A - 1)) - 1) >>> 0;
    prev_last_p_idx = 0;
  } else {
    for (i = 0; i < omega_A; i++) FB.flags[MPQS_I(h, i)] &= ~MPQS_FBE_DIVIDES_A;
    prev_last_p_idx = MPQS_I(h, omega_A - 1);

    if (room > 30) room = 30;
    const room_mask = ~((1 << room) - 1) >>> 0;

    h.bin_index = mpqs_increment(h.bin_index);
    if (h.index2_moved) {
      while ((h.bin_index & (room_mask | 0x3)) >>> 0 === 0)
        h.bin_index = mpqs_increment(h.bin_index);
    }
    if ((h.bin_index & room_mask) >>> 0 !== 0) {
      /* fell off the edge on the left */
      h.index2_FB += 2; /* caller to check this isn't too large */
      h.index2_moved = 1;
      h.bin_index = 0;
      return 0; /* back off - caller should retry */
    }
  }
  let bits = h.bin_index >>> 0;
  /* map bits to FB subscripts, counting downward with bit 0 corresponding to
   * index2_FB, and accumulate logarithms against l2_last_p */
  j = h.index2_FB;
  v2 = vals(bits);
  if (v2) {
    j -= v2;
    bits >>>= v2;
  }
  for (i = omega_A - 2; i >= 0; i--) {
    h.per_A_i[i] = j;
    l2_last_p -= MPQS_LP(h, i);
    if (FB.flags[MPQS_I(h, i)] & MPQS_FBE_DIVIDES_N) return 0; /* retry */
    FB.flags[MPQS_I(h, i)] |= MPQS_FBE_DIVIDES_A;
    bits &= ~1;
    if (!bits) break; /* i = 0 */
    v2 = vals(bits);
    bits >>>= v2;
    j -= v2;
  }
  if (missing_primes) {
    const lm = missing_primes.length;
    j = missing_primes[h.bin_index % lm];
    j += h.index0_FB - 1;
    if (h.two_is_norm) j--;
    if (FB.flags[j] & MPQS_FBE_DIVIDES_A) return 0; /* retry */
    p = FB.p[j];
  } else {
    /* choose the larger prime; index2_FB <= size_of_FB - 3.
     * mpqs.c:696: `for (j = h->index2_FB + 1; (p = FB[j].fbe_p); j++)`, i.e.
     * the loop also stops on the p = 0 sentinel, leaving p = 0. */
    for (j = h.index2_FB + 1; ; j++) {
      p = FB.p[j];
      if (!p) break;
      if (!(FB.flags[j] & MPQS_FBE_DIVIDES_N) && FB.flogp[j] > l2_last_p) break;
    }
    /* avoid re-using the same last prime as the previous A */
    if (p && j === prev_last_p_idx) {
      j++;
      p = FB.p[j];
    }
  }
  h.per_A_i[omega_A - 1] = p ? j : h.size_of_FB + 1;
  if (FB.flags[MPQS_I(h, omega_A - 1)] & MPQS_FBE_DIVIDES_N) return 0; /* retry */
  FB.flags[MPQS_I(h, omega_A - 1)] |= MPQS_FBE_DIVIDES_A;
  return 1;
}

/**
 * mpqs.c:746 (mpqs_self_init). Compute the coefficients A and B of the sieving
 * polynomial and the attached start positions. Returns false when we have run
 * out of primes for A.
 */
function mpqs_self_init(h: Handle, missing_primes: number[] | null): boolean {
  const size_of_FB = h.size_of_FB + 1;
  const FB = h.FB;
  const omega_A = h.omega_A;
  const w1 = Math.max(1, omega_A - 1);
  let i: number;
  let j: number;

  h.index_j = (h.index_j + 1) >>> 0;
  if (h.index_j === h.no_B) {
    /* all the B's have been used, choose new A */
    h.index_j = 0;
    h.index_i = (h.index_i + 1) >>> 0;
  }

  if (missing_primes || h.index_j === 0) {
    /* compute first polynomial with new A */
    while (mpqs_si_choose_primes(h, missing_primes) <= 0) {
      if (size_of_FB - h.index2_FB < 4) return false; /* fail */
    }
    /* compute A = product of omega_A primes given by bin_index */
    let a = 1n;
    for (i = 0; i < omega_A; i++) a *= BigInt(MPQS_AP(h, i));
    h.A = a;
    /* compute H[i] and the initial B = sum(H[i]) */
    let b = 0n;
    for (i = 0; i < omega_A; i++) {
      const p = MPQS_AP(h, i);
      const P = BigInt(p);
      const t0 = h.A / P;
      const t = (t0 * BigInt(Fl_inv(umodiu(t0, p), p) * MPQS_SQRT(h, i))) % h.A;
      h.per_A_H[i] = t;
      b += t;
    }
    if (h.kN & 1n) {
      /* ensure b = 1 mod 4 */
      if ((b & 1n) === 0n) b += (h.A % 4n) * h.A;
    } else {
      /* ensure b = 0 mod 2 */
      if (b & 1n) b += h.A;
    }
    h.B = b;

    const A2 = h.A << 1n;
    /* roots z1, z2 of Q(x) mod p_j, and the start positions.
     * Primes dividing A are skipped here (handled in the common part). */
    for (j = 3; j <= size_of_FB; j++) {
      if (FB.flags[j] & MPQS_FBE_DIVIDES_A) continue;
      const p = FB.p[j];
      const m = h.M % p;
      const iA2 = Fl_inv(umodiu(A2, p), p); /* 1/(2A) mod p_j */
      let iA = iA2 << 1;
      if (iA > p) iA -= p;
      let mb = umodiu(h.B, p);
      if (mb) mb = p - mb; /* -B mod p */
      const s = FB.sqrt_kN[j];
      const t = Fl_add(m, Fl_mul(Fl_sub(mb, s, p), iA2, p), p);
      FB.start1[j] = t;
      FB.start2[j] = Fl_add(t, Fl_mul(s, iA, p), p);
      for (i = 0; i < omega_A - 1; i++) {
        const hh = umodiu(h.per_A_H[i], p);
        h.inv_A_H[j * w1 + i] = Fl_mul(hh, iA, p); /* 1/A * H[i] mod p_j */
      }
    }
  } else {
    /* no "real" computation -- use the recursive (Gray code) formula */
    const v2 = vals(h.index_j);
    j = h.index_j >>> v2;
    const p1 = h.per_A_H[v2] << 1n;
    if (j & 2) {
      /* j = 3 mod 4 */
      for (j = 3; j <= size_of_FB; j++) {
        if (FB.flags[j] & MPQS_FBE_DIVIDES_A) continue;
        const p = FB.p[j];
        const d = h.inv_A_H[j * w1 + v2];
        FB.start1[j] = Fl_sub(FB.start1[j], d, p);
        FB.start2[j] = Fl_sub(FB.start2[j], d, p);
      }
      h.B = h.B + p1;
    } else {
      /* j = 1 mod 4 */
      for (j = 3; j <= size_of_FB; j++) {
        if (FB.flags[j] & MPQS_FBE_DIVIDES_A) continue;
        const p = FB.p[j];
        const d = h.inv_A_H[j * w1 + v2];
        FB.start1[j] = Fl_add(FB.start1[j], d, p);
        FB.start2[j] = Fl_add(FB.start2[j], d, p);
      }
      h.B = h.B - p1;
    }
  }

  /* p = 2 is a special case: start1[2], start2[2] are never looked at.
   * Compute the zeros of the polynomials that have only one zero mod p (p | A) */
  const mC = (h.kN - h.B * h.B) / (h.A << 2n); /* coefficient -C */
  for (i = 0; i < omega_A; i++) {
    const p = MPQS_AP(h, i);
    const s = h.M + Fl_div(umodiu(mC, p), umodiu(h.B, p), p);
    FB.start1[MPQS_I(h, i)] = FB.start2[MPQS_I(h, i)] = s % p;
  }
  return true;
}

/* ------------------------------------------------------------------ */
/* THE SIEVE (mpqs.c:910-1034)                                         */
/* ------------------------------------------------------------------ */

/**
 * mpqs.c:965 (mpqs_sieve). Add the scaled logarithms of the FB primes to the
 * sieve array at the positions of their two arithmetic progressions.
 *
 * The three loops of the original differ only in the amount of unrolling
 * (mpqs_sieve_p / _p1 / _p2, mpqs.c:914-964); all of them add `logp` to
 * exactly the positions < 2M of the two progressions.
 */
function mpqs_sieve(h: Handle): void {
  const FB = h.FB;
  const S = h.sieve_array;
  const size = h.M << 1;
  S.fill(0, 0, size);
  for (let l = h.index1_FB; ; l++) {
    const p = FB.p[l];
    if (!p) break;
    const logp = FB.logval[l];
    const s1 = FB.start1[l];
    const s2 = FB.start2[l];
    for (let x = s1; x < size; x += p) S[x] += logp;
    if (s2 !== s1) for (let x = s2; x < size; x += p) S[x] += logp;
  }
}

/**
 * mpqs.c:1007 (mpqs_eval_sieve). Collect the subscripts whose sieve value
 * reached the threshold. Returns the number of candidates.
 *
 * We use the 8-byte (LONG_IS_64BIT, no SSE2) bit array layout of mpqs.c:93.
 */
function mpqs_eval_sieve(h: Handle): number {
  let x = 0;
  let count = 0;
  const M2 = h.M << 1;
  const t = h.sieve_threshold;
  const S = h.sieve_array;
  const W = h.sieve_words;
  const cand = h.candidates;
  const sizemask = 8;

  while (count < MPQS_CANDIDATE_ARRAY_SIZE - 1) {
    /* the sentinel S[2M] = 255 stops this loop */
    while (!((W[2 * x] | W[2 * x + 1]) & 0x80808080)) x++;
    let y = x * sizemask;
    for (let j = 0; j < sizemask; j++, y++) {
      if (y >= M2) {
        cand[count] = 0;
        return count;
      }
      if (S[y] >= t) cand[count++] = y;
    }
    x++;
  }
  cand[count] = 0;
  return count;
}

/* ------------------------------------------------------------------ */
/* CONSTRUCTING RELATIONS (mpqs.c:1036-1386)                           */
/* ------------------------------------------------------------------ */

/** A relation: Y and the list of (exponent, FB subscript) pairs, packed as in
 * mpqs.c:1130 (mpqs_add_factor): `pi | (ei << REL_OFFSET)`. */
interface Rel {
  Y: bigint;
  relp: number[];
}

function mpqs_add_factor(relp: number[], ei: number, pi: number): void {
  /* PARI packs this into a long and unpacks with an *arithmetic* shift, so a
   * negative exponent (class group mode, rel_sub_ei) round-trips; we do the
   * same in 32 bits, which bounds |ei| by 2^11 instead of PARI's 2^43. */
  relp.push(pi | (ei << REL_OFFSET));
}

/** mpqs.c:1095 (rel_to_ei) */
function rel_to_ei(ei: Int32Array, relp: number[]): void {
  for (const r of relp) {
    const e = r >> REL_OFFSET;
    const i = r & REL_MASK;
    ei[i] += e;
  }
}

/** mpqs.c:1106 (rel_add_ei), class group mode only */
function rel_add_ei(h: Handle, ei: Int32Array, relp: number[], b: bigint): void {
  for (const r of relp) {
    const e = r >> REL_OFFSET;
    const i = r & REL_MASK;
    const p = h.FB.p[i];
    ei[i] += umodiu(b, p << 1) > p ? -e : e;
  }
}

/** mpqs.c:1118 (rel_sub_ei), class group mode only */
function rel_sub_ei(h: Handle, ei: Int32Array, relp: number[], b: bigint): void {
  for (const r of relp) {
    const e = r >> REL_OFFSET;
    const i = r & REL_MASK;
    const p = h.FB.p[i];
    ei[i] -= umodiu(b, p << 1) > p ? -e : e;
  }
}

/** mpqs.c:1134 (zv_is_even) */
function zv_is_even(ei: Int32Array, lei: number): boolean {
  for (let i = 1; i <= lei; i++) if (ei[i] & 1) return false;
  return true;
}

/** modular inverse; on failure returns the gcd (PARI: invmod) */
function invmod(a: bigint, N: bigint): { inv: bigint | null; gcd: bigint } {
  let [old_r, r] = [((a % N) + N) % N, N];
  let [old_s, s] = [1n, 0n];
  while (r !== 0n) {
    const q = old_r / r;
    [old_r, r] = [r, old_r - q * r];
    [old_s, s] = [s, old_s - q * s];
  }
  if (old_r !== 1n) return { inv: null, gcd: old_r };
  return { inv: ((old_s % N) + N) % N, gcd: 1n };
}

/**
 * mpqs.c:1143 (combine_large_primes). Combine two relations sharing the same
 * large prime q into a full relation. Returns the relation, or a factor of N
 * (as a bigint, "very unlikely"), or null.
 */
function combine_large_primes(
  h: Handle,
  q: number,
  rel1: Rel,
  rel2: Rel,
  mode: number
): Rel | bigint | null {
  const Y1 = rel1.Y;
  const Y2 = rel2.Y;
  const lei = h.size_of_FB + 1;
  const { inv: iq, gcd } = invmod(BigInt(q), h.N);
  if (iq === null) return gcd === h.N ? null : gcd; /* rare */
  const ei = new Int32Array(lei + 1);
  let new_Y: bigint;
  if (mode === MPQS_MODE_CLASSGROUP) {
    rel_add_ei(h, ei, rel1.relp, Y1);
    const Q = BigInt(q);
    if (Y1 % Q === Y2 % Q) rel_sub_ei(h, ei, rel2.relp, Y2);
    else rel_add_ei(h, ei, rel2.relp, Y2);
    new_Y = 0n;
  } else {
    rel_to_ei(ei, rel1.relp);
    rel_to_ei(ei, rel2.relp);
    if (zv_is_even(ei, lei)) return null;
    new_Y = (Y1 * Y2 * iq) % h.N;
    const new_Y1 = h.N - new_Y;
    if (new_Y1 < new_Y) new_Y = new_Y1;
  }
  const relp: number[] = [];
  if (ei[1] & 1) mpqs_add_factor(relp, 1, 1);
  for (let l = 2; l <= lei; l++) if (ei[l]) mpqs_add_factor(relp, ei[l], l);
  return { Y: new_Y, relp };
}

/**
 * mpqs.c:1056 (mpqs_factorback, MPQS_DEBUG only): the product of the FB primes
 * of a relation, with their exponents, mod N.
 */
function mpqs_factorback(h: Handle, relp: number[]): bigint {
  let Q = 1n;
  for (const r of relp) {
    const e = r >> REL_OFFSET;
    const i = r & REL_MASK;
    if (i === 1) Q = (h.N - Q) % h.N; /* special case -1 */
    else Q = (Q * modpow(BigInt(h.FB.p[i]), e, h.N)) % h.N;
  }
  return Q;
}

/**
 * mpqs.c:1069 (mpqs_check_rel, MPQS_DEBUG only): verify that Y^2 = q * prod
 * p_i^e_i mod N. Enabled by the `debug` option; upstream compiles it in with
 * -DMPQS_DEBUG.
 */
function mpqs_check_rel(h: Handle, c: Rel, q: number, mode: number): void {
  const Y = c.Y;
  const Qx_2 = (Y * Y) % h.N;
  if (mode === MPQS_MODE_CLASSGROUP) {
    if (Y === 0n || q !== 1) return;
    q = 4;
  }
  const rhs = (mpqs_factorback(h, c.relp) * BigInt(q)) % h.N;
  if (Qx_2 !== rhs) {
    /* mpqs.c:1087: the message depends on q, exactly as upstream */
    throw new Error(
      q ? 'MPQS: wrong large prime relation found' : 'MPQS: wrong full relation found'
    );
  }
}

/** key identifying a relation for the `frel` hash table (mpqs.c:103 frel_add) */
function relKey(r: Rel): string {
  return r.Y.toString(36) + '|' + r.relp.join(',');
}

/**
 * mpqs.c:1199 (mpqs_eval_cand). Evaluate `nc` candidates, adding full
 * relations to `frel` and large-prime relations to `lprel`. Returns a factor of
 * N when one turns up while combining, else null.
 */
function mpqs_eval_cand(
  h: Handle,
  nc: number,
  frel: Map<string, Rel>,
  lprel: Map<number, Rel>,
  mode: number
): bigint | null {
  const FB = h.FB;
  const A = h.A;
  const B = h.B;
  const relaprimes = h.relaprimes;
  const candidates = h.candidates;
  const two_bad = h.two_is_bad;

  for (let i = 0; i < nc; i++) {
    const relp: number[] = [];
    const x = candidates[i];
    let relaprpos = 0;
    let thr = h.sieve_array[x];
    /* Y = 2*A*x + B, Qx = Y^2/(4*A) = Q(x) */
    const Y = A * BigInt(2 * (x - h.M)) + B;
    let Qx = Y * Y - h.kN; /* != 0 since N is not a square and (N,k) = 1 */
    if (Qx < 0n) {
      Qx = -Qx;
      mpqs_add_factor(relp, 1, 1); /* i = 1: the -1 slot */
    }
    /* divide by powers of 2: we are really dealing with 4*A*Q(x) */
    let powers_of_2: number;
    if (two_bad) powers_of_2 = 2;
    else {
      powers_of_2 = 0;
      while ((Qx & 1n) === 0n) {
        Qx >>= 1n;
        powers_of_2++;
      }
    }
    if (two_bad) Qx >>= 2n;
    if (mode === MPQS_MODE_CLASSGROUP) {
      if (powers_of_2 !== 2) mpqs_add_factor(relp, powers_of_2 - 2, 2);
    } else mpqs_add_factor(relp, powers_of_2, 2);

    /* Pass 1 over the odd primes in FB: pick up all possible divisors of Qx */
    let Qx_part = A;
    let pi: number;
    for (pi = 3; pi < h.index1_FB; pi++) {
      const p = FB.p[pi];
      const xp = x % p;
      if (xp === FB.start1[pi] || xp === FB.start2[pi]) {
        /* p divides Q(x)/A and possibly A */
        const ei = FB.flags[pi] & MPQS_FBE_DIVIDES_A;
        relaprimes[relaprpos++] = pi;
        relaprimes[relaprpos++] = 1 + ei;
        Qx_part *= BigInt(p);
      }
    }
    for (; thr; pi++) {
      const p = FB.p[pi];
      if (!p) break;
      const xp = x % p;
      if (xp === FB.start1[pi] || xp === FB.start2[pi]) {
        const ei = FB.flags[pi] & MPQS_FBE_DIVIDES_A;
        relaprimes[relaprpos++] = pi;
        relaprimes[relaprpos++] = 1 + ei;
        Qx_part *= BigInt(p);
        thr = (thr - FB.logval[pi]) & 255;
      }
    }
    for (let k = 0; k < h.omega_A; k++) {
      const ppi = MPQS_I(h, k);
      const p = FB.p[ppi];
      const xp = x % p;
      if (!(xp === FB.start1[ppi] || xp === FB.start2[ppi])) {
        /* p divides A but does not divide Q(x)/A */
        relaprimes[relaprpos++] = ppi;
        relaprimes[relaprpos++] = 0;
      }
    }
    /* divide off what we know */
    Qx = Qx / Qx_part;

    /* Pass 2: deal with repeated factors and store the tentative relation */
    for (let pii = 0; pii < relaprpos; pii += 2) {
      let ei = relaprimes[pii + 1];
      const pj = relaprimes[pii];
      /* p | k (index below index0_FB) or p | A (ei = 0) */
      if (pj < h.index0_FB || ei === 0) {
        mpqs_add_factor(relp, 1, pj);
        continue;
      }
      const p = BigInt(FB.p[pj]);
      /* p might still divide the current adjusted Qx */
      if (Qx === p) {
        ei++;
        Qx = 1n;
      } else if (Qx > p) {
        while (Qx % p === 0n) {
          ei++;
          Qx /= p;
        }
      }
      mpqs_add_factor(relp, ei, pj);
    }

    if (Qx === 1n) {
      const rel: Rel = { Y: Y < 0n ? -Y : Y, relp };
      if (h.debug) mpqs_check_rel(h, rel, 1, mode);
      const key = relKey(rel);
      if (!frel.has(key)) frel.set(key, rel);
    } else if (Qx <= BigInt(h.lp_bound) && mode !== MPQS_MODE_CLASSGROUP) {
      const q = Number(Qx);
      const rel: Rel = { Y: Y < 0n ? -Y : Y, relp };
      if (h.debug) mpqs_check_rel(h, rel, q, mode);
      const col = lprel.get(q);
      if (!col) lprel.set(q, rel);
      else {
        const c = combine_large_primes(h, q, rel, col, mode);
        if (c !== null) {
          if (typeof c === 'bigint') return c; /* very unlikely */
          if (h.debug) mpqs_check_rel(h, c, 1, mode);
          const key = relKey(c);
          if (!frel.has(key)) frel.set(key, c);
        }
      }
    }
  }
  return null;
}

/* ------------------------------------------------------------------ */
/* FROM RELATIONS TO DIVISORS (mpqs.c:1388-1606)                       */
/* ------------------------------------------------------------------ */

/** mpqs.c:1424 (rels_to_F2Ms): create an F2Ms (list of columns) from relations */
function rels_to_F2Ms(rels: Rel[]): number[][] {
  const m: number[][] = [];
  for (const rel of rels) {
    const rel2: number[] = [];
    for (const r of rel.relp) if ((r >> REL_OFFSET) & 1) rel2.push(r & REL_MASK);
    m.push(rel2);
  }
  return m;
}

/**
 * F2v.c:397 (F2m_ker_sp, deplin = 0): kernel of a matrix over F2 given by its
 * columns, each column a bitset over rows 1..nbrow.
 *
 * Returns the kernel basis: an array of bitsets over the column indices 1..n.
 * (F2v.c:1063 F2Ms_ker switches to block Lanczos above 640 rows; see the
 * porting notes at the top of this file.)
 */
function F2Ms_ker(M: number[][], nbrow: number): Uint32Array[] {
  const n = M.length;
  const W = (nbrow >>> 5) + 1;
  const x = new Uint32Array(n * W);
  for (let k = 0; k < n; k++) {
    for (const row of M[k]) x[k * W + (row >>> 5)] |= 1 << (row & 31);
  }
  /* c[j] = 1 while row j has not been used as a pivot */
  const cW = (nbrow >>> 5) + 1;
  const c = new Uint32Array(cW);
  for (let j = 1; j <= nbrow; j++) c[j >>> 5] |= 1 << (j & 31);
  const d = new Int32Array(n + 2);
  let r = 0;
  for (let k = 1; k <= n; k++) {
    const off = (k - 1) * W;
    /* F2v_find_nonzero: smallest j with x[k][j] = 1 and c[j] = 1 */
    let j = nbrow + 1;
    for (let w = 0; w < W; w++) {
      const e = x[off + w] & c[w];
      if (e) {
        j = (w << 5) + vals(e);
        break;
      }
    }
    if (j > nbrow) {
      r++;
      d[k] = 0;
    } else {
      c[j >>> 5] &= ~(1 << (j & 31));
      d[k] = j;
      x[off + (j >>> 5)] &= ~(1 << (j & 31));
      for (let i = k + 1; i <= n; i++) {
        const offi = (i - 1) * W;
        if (x[offi + (j >>> 5)] & (1 << (j & 31))) {
          for (let w = 0; w < W; w++) x[offi + w] ^= x[off + w];
        }
      }
      x[off + (j >>> 5)] |= 1 << (j & 31);
    }
  }
  const yW = (n >>> 5) + 1;
  const y: Uint32Array[] = [];
  for (let j = 1, k = 1; j <= r; j++, k++) {
    const C = new Uint32Array(yW);
    while (d[k]) k++;
    const offk = (k - 1) * W;
    for (let i = 1; i < k; i++) {
      const di = d[i];
      if (di && x[offk + (di >>> 5)] & (1 << (di & 31))) C[i >>> 5] |= 1 << (i & 31);
    }
    C[k >>> 5] |= 1 << (k & 31);
    y.push(C);
  }
  return y;
}

/**
 * mpqs.c:1444 (split). Replace D by its root when D is a perfect power and set
 * the exponent; return true when D is a probable prime or a proper power.
 */
function split(D: bigint): { ok: boolean; D: bigint; e: number } {
  if (isPrime(D)) return { ok: true, D, e: 1 }; /* PARI: MR_Jaeschke */
  const s = Z_issquareall(D);
  if (s !== null) return { ok: true, D: s, e: 2 };
  const [flag, root] = is_357_power(D, 7);
  if (flag) return { ok: true, D: root, e: flag };
  return { ok: false, D, e: 0 }; /* known composite */
}

/**
 * mpqs.c:1468 (mpqs_solve_linear_system). Gaussian elimination over F2 on the
 * exponent vectors, then gcds. Returns [value, exponent] pairs whose product
 * is N (PARI's ifac format), or null when no factor was found.
 */
function mpqs_solve_linear_system(
  h: Handle,
  frel: Map<string, Rel>
): Array<[bigint, bigint]> | null {
  const FB = h.FB;
  const N = h.N;
  const rels = [...frel.values()];
  const M = rels_to_F2Ms(rels);
  const Ker = F2Ms_ker(M, h.size_of_FB + 1);
  const rank = Ker.length;
  if (!rank) return null; /* trivial kernel: main loop may look for more */

  const ei = new Int32Array(h.size_of_FB + 2);
  let rmax = logint(N, 3n);
  if (rank <= 62) rmax = Math.min(rmax, 2 ** rank);
  const r: bigint[] = new Array(rmax + 1).fill(0n);
  const c: number[] = new Array(rmax + 1).fill(0);
  let rnext = 1;
  let rlast = 1;
  const nrows = M.length;

  for (let i = 1; i <= rank; i++) {
    /* loop over the kernel basis */
    let X = 1n;
    let Y_prod = 1n;
    let done = 0;
    ei.fill(0);
    const K = Ker[i - 1];
    for (let j = 1; j <= nrows; j++) {
      if (K[j >>> 5] & (1 << (j & 31))) {
        const R = rels[j - 1];
        Y_prod = (Y_prod * R.Y) % N;
        rel_to_ei(ei, R.relp);
      }
    }
    for (let j = 2; j <= h.size_of_FB + 1; j++) {
      if (ei[j]) {
        const q = BigInt(FB.p[j]);
        if (ei[j] & 1) throw new Error('MPQS (relation is a nonsquare)');
        X = (X * modpow(q, ei[j] >> 1, N)) % N;
      }
    }
    /* mpqs.c:1525 (MPQS_DEBUGLEVEL >= 1): X^2 - Y^2 must be divisible by N */
    if (h.debug && (X * X - Y_prod * Y_prod) % N !== 0n) {
      throw new Error('MPQS: wrong relation found after Gauss');
    }
    /* gcd(X-Y,N) * gcd(X+Y,N) = N and X is coprime to N, so gcd(X+Y,N) alone */
    const X_plus_Y = X + Y_prod;
    if (rnext === 1) {
      /* we still haven't decomposed, and want both a gcd and its cofactor */
      const D = gcdii(X_plus_Y, N);
      if (D === 1n || D === N) continue;
      r[1] = N / D;
      r[2] = D;
      rlast = rnext = 3;
      const s1 = split(r[1]);
      r[1] = s1.D;
      c[1] = s1.e;
      if (s1.ok) done++;
      const s2 = split(r[2]);
      r[2] = s2.D;
      c[2] = s2.e;
      if (s2.ok) done++;
      if (done === 2 || rmax === 2) break;
    } else {
      /* we already have factors */
      for (let j = 1; j < rnext; j++) {
        /* loop over known-composite factors */
        if (c[j]) {
          done++;
          continue;
        }
        const D = gcdii(X_plus_Y, r[j]);
        if (D === 1n || D === r[j]) continue;
        r[j] = r[j] / D;
        r[rnext] = D;
        const sj = split(r[j]);
        r[j] = sj.D;
        c[j] = sj.e;
        if (sj.ok) done++;
        const sn = split(r[rnext]);
        r[rnext] = sn.D;
        c[rnext] = sn.e;
        if (++rnext > rmax) break;
      }
      if (rnext > rlast) rlast = rnext;
      if (rnext > rmax || done === rnext - 1) break;
    }
  }
  if (rnext === 1) return null; /* no factors found */

  rlast = rnext - 1; /* # of distinct factors found */
  const res: Array<[bigint, bigint]> = [];
  for (let i = 1; i <= rlast; i++) {
    const C = c[i];
    res.push([r[i], C <= 1 ? 1n : BigInt(C)]);
  }
  return res;
}

/** x^e mod m, e a nonnegative machine integer */
function modpow(x: bigint, e: number, m: bigint): bigint {
  let r = 1n;
  let b = x % m;
  while (e > 0) {
    if (e & 1) r = (r * b) % m;
    b = (b * b) % m;
    e >>>= 1;
  }
  return r;
}

/* ------------------------------------------------------------------ */
/* MAIN ENTRY POINT AND DRIVER ROUTINE (mpqs.c:1608-1773)              */
/* ------------------------------------------------------------------ */

function newHandle(): Handle {
  return {
    N: 0n,
    kN: 0n,
    A: 0n,
    B: 0n,
    FB: null as unknown as FBArrays,
    inv_A_H: new Int32Array(0),
    per_A_H: [],
    per_A_i: new Int32Array(0),
    sieve_array: new Uint8Array(0),
    sieve_words: new Uint32Array(0),
    candidates: new Int32Array(0),
    relaprimes: new Int32Array(0),
    M: 0,
    size_of_FB: 0,
    index0_FB: 0,
    index1_FB: 0,
    index2_FB: 0,
    index2_moved: 0,
    sieve_threshold: 0,
    two_is_norm: 0,
    two_is_bad: 0,
    omega_A: 0,
    no_B: 0,
    l2_target_A: 0,
    bin_index: 0,
    index_i: 0,
    index_j: 0,
    target_rels: 0,
    largest_FB_p: 0,
    pmin_index1: 0,
    lp_scale: 0,
    lp_bound: 0,
    digit_size_kN: 0,
    _k: cand_multipliers[0],
    tolerance: 0,
    dkN: 0,
    l2sqrtkN: 0,
    l2M: 0,
    debug: false,
  };
}

/** Options of our port; PARI has no equivalent (it never gives up early). */
export interface MpqsOptions {
  /**
   * Optional cap on the number of polynomials tried, so that a caller can put
   * a bound on the running time. PARI has no such cap: it stops only when it
   * runs out of primes for A, or when Gauss failed with 1.5*target_rels
   * relations (mpqs.c:1770). 0 (the default) means "no cap".
   */
  maxPolys?: number;
  /**
   * Turn on upstream's own consistency checks, which it compiles in with
   * -DMPQS_DEBUG: every relation is verified (mpqs.c:1069 mpqs_check_rel) and
   * X^2 = Y^2 mod N is checked after Gauss (mpqs.c:1525). Off by default, as
   * upstream.
   */
  debug?: boolean;
}

/**
 * mpqs.c:1639 (mpqs). Factor N using the self-initializing multipolynomial
 * quadratic sieve. N must be odd, composite, not a perfect power and have no
 * small prime factors (PARI only ever calls this from ifac_crack).
 *
 * @returns [value, exponent] pairs (PARI's ifac format) whose product is N, or
 * null when we can't seem to make any headway (PARI returns NULL there too).
 */
export function mpqs(N: bigint, options?: MpqsOptions): Array<[bigint, bigint]> | null {
  const size_N = decimal_len(N);
  const maxPolys = options?.maxPolys ?? 0;
  if (size_N > MPQS_MAX_DIGIT_SIZE_KN) return null; /* toolarge() */

  const H = newHandle();
  H.debug = options?.debug ?? false;
  H.N = N;
  H.two_is_norm = 0;
  H.two_is_bad = 0;
  H.bin_index = 0;
  H.index_i = 0;
  H.index_j = 0;
  H.index2_moved = 0;
  const p0 = mpqs_find_k(H);
  if (p0)
    return [
      [BigInt(p0), 1n],
      [N / BigInt(p0), 1n],
    ];
  H.kN = N * BigInt(H._k.k);
  if (!mpqs_set_parameters(H)) return null; /* toolarge() */

  const pf = mpqs_create_FB(H, true);
  if (pf)
    return [
      [BigInt(pf), 1n],
      [N / BigInt(pf), 1n],
    ];
  mpqs_sieve_array_ctor(H);
  mpqs_poly_ctor(H);

  H.lp_bound = Math.min(H.largest_FB_p, MPQS_LP_BOUND);
  /* don't allow large primes to have room for two factors both bigger than
   * what the FB contains */
  H.lp_bound *= Math.min(H.lp_scale, H.largest_FB_p - 1);
  H.dkN = Number(H.kN);
  mpqs_set_sieve_threshold(H);
  if (!mpqs_locate_A_range(H)) return null;

  /* Let (A, B_i) be the current pair of coeffs. If i == 0 a new A is generated */
  H.index_j = 0xffffffff; /* (mpqs_uint32_t)-1: increment below starts at 0 */

  const frel = new Map<string, Rel>();
  const lprel = new Map<number, Rel>();
  /* mpqs.c:1699: computed once, from the initial target_rels */
  const DEFEAT = H.target_rels * 1.5;
  let polys = 0;
  for (;;) {
    /* self initialization: compute polynomial and its zeros */
    if (!mpqs_self_init(H, null)) return null; /* ran out of primes for A */
    if (maxPolys && ++polys > maxPolys) return null;
    mpqs_sieve(H);
    const tc = mpqs_eval_sieve(H);
    if (tc) {
      const fact = mpqs_eval_cand(H, tc, frel, lprel, MPQS_MODE_FACTOR);
      if (fact)
        return [
          [fact, 1n],
          [N / fact, 1n],
        ]; /* factor found while combining */
    }
    if (frel.size < H.target_rels) continue; /* main loop */

    const fact = mpqs_solve_linear_system(H, frel);
    if (fact) return fact;
    if (frel.size >= DEFEAT) return null;
    H.target_rels += 10;
  }
}

/**
 * Internals exposed for the test suite only (they are `static` in mpqs.c and
 * are not part of the module's interface).
 */
export const mpqsInternals = {
  mpqs_increment,
  Fl_sqrt,
  Fl_inv,
  krouu,
  kroiu,
  F2Ms_ker,
  decimal_len,
  logint,
  mpqs_parameters,
  cand_multipliers,
};

/* mpqs.c:1775 (mpqs_class_init) and mpqs.c:1815 (mpqs_class_rels) are the
 * class-group entry points; they are not ported because buch2.c, their only
 * caller, is not ported. See the porting notes at the top of this file. */
