/**
 * @module ifactor
 * @description Integer factorization ported from PARI/GP
 *
 * Reference: reference/pari/src/basemath/ifactor1.c
 *            reference/pari/src/basemath/ispower.c (pure power detection)
 *
 * PARI's factorization (`ifactor_sign`, ifactor1.c:4279) proceeds as follows:
 *
 * 1. Trial division by primes up to `tridiv_bound(n)` (ifactor1.c:3158), using
 *    the gcd-with-primorial trick of `Z_oddprimedivisors_fast` (ifactor1.c:3306).
 * 2. BPSW primality test on what is left (`ifac_isprime`, ifactor1.c:2514).
 * 3. If still composite, `ifac_decomp` (ifactor1.c:3006) repeatedly cracks
 *    composite entries with `ifac_crack` (ifactor1.c:2786), which tries, in
 *    this order:
 *      a. pure powers (square, then 3/5/7-th, then p-th powers),
 *      b. Shanks' SQUFOF          (`squfof`,       ifactor1.c:1474),
 *      c. Pollard-Brent rho       (`pollardbrent`, ifactor1.c:1361),
 *      d. Lenstra-Montgomery ECM  (`ellfacteur`,   ifactor1.c:1038), not insisting,
 *      e. MPQS                    (`mpqs`, mpqs.c) -- NOT PORTED, see below,
 *      f. ECM again, insisting.
 *
 * Everything except MPQS is ported here. Deviations from PARI are marked with
 * `@see Deviation:` and listed at the top of each affected function:
 *
 * - MPQS (the multiple polynomial quadratic sieve) is not ported. Its slot in
 *   the chain is taken by a Pollard-Brent run with a much larger round budget,
 *   which is a complete (if exponentially slower) algorithm. If every stage
 *   fails, `Z_factor` throws `NotImplementedError` rather than declaring a
 *   composite prime (which is what PARI does only at `DEBUGLEVEL >= 2` with a
 *   warning, and only when factorization was explicitly bounded).
 * - ECM is run one curve at a time instead of `nbc` curves in parallel: PARI
 *   batches the modular inversions across curves with Montgomery's trick,
 *   which is a constant-factor speedup with no effect on which factors are
 *   found. Same curve family, same B1 schedule, same seeds.
 * - Where PARI uses floating point square/n-th roots to guess a root before
 *   verifying it exactly, we use exact integer Newton roots (project rule: no
 *   floating point). The result is identical.
 */

import { Fp_pow } from './ff.js';

/**
 * Thrown when factorization would require an algorithm we have not ported.
 * Mirrors sagemath-ts' `NotImplementedError`; parigp-ts cannot import from
 * sagemath-ts (dependency direction), so it is re-declared here.
 */
export class NotImplementedError extends Error {
  override name = 'NotImplementedError';

  constructor(message: string = 'not implemented') {
    super(message);
  }
}

export type Factorization = Array<[bigint, bigint]>;

// ============================================================================
// Small helpers on bigints
// ============================================================================

/** Integer square root (floor), exact for arbitrarily large input. */
export function isqrt(n: bigint): bigint {
  if (n < 0n) throw new Error('isqrt: negative argument');
  if (n < 2n) return n;
  let x = 1n << (BigInt(n.toString(2).length + 1) >> 1n);
  for (;;) {
    const x1 = (x + n / x) >> 1n;
    if (x1 >= x) return x;
    x = x1;
  }
}

/**
 * Exact integer k-th root: `[r, exact]` with `r = floor(x^(1/k))`.
 *
 * PARI computes an approximate root with `sqrtnr` (floating point) and then
 * verifies `y^k == x` exactly (ifactor1.c:2054); we compute the root itself
 * exactly by Newton iteration, which yields the same `y` and the same verdict.
 */
export function Z_iroot(x: bigint, k: number): [bigint, boolean] {
  if (k <= 0) throw new Error('Z_iroot: k must be positive');
  if (k === 1) return [x, true];
  if (x < 2n) return [x, true];
  const K = BigInt(k);
  const bits = BigInt(x.toString(2).length);
  let r = 1n << ((bits + K - 1n) / K);
  for (;;) {
    const next = ((K - 1n) * r + x / r ** (K - 1n)) / K;
    if (next >= r) break;
    r = next;
  }
  return [r, r ** K === x];
}

/** gcd of two nonnegative bigints. */
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

/** `[v, q]` with `n = p^v * q`, `p` not dividing `q` (PARI: Z_lvalrem). */
function lvalrem(n: bigint, p: bigint): [number, bigint] {
  let v = 0;
  while (n % p === 0n) {
    n /= p;
    v++;
  }
  return [v, n];
}

/**
 * Modular inverse with factor detection.
 * PARI's `invmod(a,N,&g)` returns 0 and sets `g = gcd(a,N)` on failure
 * (ifactor1.c uses this to spot factors of N during curve arithmetic).
 */
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

// ============================================================================
// Prime generation (PARI: forprime.c, a byte sieve; same here)
// ============================================================================

let basePrimesCache: number[] = [];
let basePrimesLimit = 0;

/** All primes <= limit, by a simple sieve of Eratosthenes (cached, growing). */
function primesUpTo(limit: number): number[] {
  if (limit <= basePrimesLimit) {
    if (limit === basePrimesLimit) return basePrimesCache;
    const out: number[] = [];
    for (const p of basePrimesCache) {
      if (p > limit) break;
      out.push(p);
    }
    return out;
  }
  const sieve = new Uint8Array(limit + 1);
  const primes: number[] = [];
  for (let i = 2; i <= limit; i++) {
    if (!sieve[i]) {
      primes.push(i);
      if (i <= limit / i) for (let j = i * i; j <= limit; j += i) sieve[j] = 1;
    }
  }
  basePrimesCache = primes;
  basePrimesLimit = limit;
  return primes;
}

const SIEVE_SEGMENT = 1 << 16;

/** Iterate over primes p with a <= p <= b (PARI: u_forprime_init/next). */
export function* forprime(a: number, b: number): Generator<number> {
  if (b < 2) return;
  if (a < 2) a = 2;
  const base = primesUpTo(Math.floor(Math.sqrt(b)) + 1);
  for (let lo = a; lo <= b; lo += SIEVE_SEGMENT) {
    const hi = Math.min(lo + SIEVE_SEGMENT - 1, b);
    const seg = new Uint8Array(hi - lo + 1);
    for (const p of base) {
      if (p * p > hi) break;
      const start = Math.max(p * p, Math.ceil(lo / p) * p);
      for (let j = start; j <= hi; j += p) seg[j - lo] = 1;
    }
    for (let i = lo; i <= hi; i++) {
      if (!seg[i - lo] && i >= 2) yield i;
    }
  }
}

/**
 * Stateful prime iterator, the analogue of PARI's `forprime_t` (used by
 * `is_pth_power`, ifactor1.c:2073).
 */
interface ForprimeT {
  next(): number | null;
}

function forprime_init(a: number, b: number): ForprimeT {
  const gen = forprime(a, b);
  return {
    next() {
      const r = gen.next();
      return r.done ? null : r.value;
    },
  };
}

// ============================================================================
// Primality testing -- BPSW (PARI: BPSW_psp, ifac_isprime ifactor1.c:2514)
// ============================================================================

/** First 168 primes; used for trial division inside the primality test. */
const SMALL_PRIMES: readonly bigint[] = primesUpTo(1000).map(BigInt);

/**
 * The first 26 primes, PARI's `tinyprimes` (ispower.c:903) truncated at 101
 * exactly as PARI's `Z_isanypower_aux`/`isprimepower_i` use it.
 */
const TINY_PRIMES: readonly bigint[] = [
  2n,
  3n,
  5n,
  7n,
  11n,
  13n,
  17n,
  19n,
  23n,
  29n,
  31n,
  37n,
  41n,
  43n,
  47n,
  53n,
  59n,
  61n,
  67n,
  71n,
  73n,
  79n,
  83n,
  89n,
  97n,
  101n,
];

/** Miller-Rabin test of `n` for a single witness `a` (n odd > 2). */
function millerRabinWitness(n: bigint, a: bigint): boolean {
  let d = n - 1n;
  let s = 0;
  while ((d & 1n) === 0n) {
    d >>= 1n;
    s++;
  }
  let x = Fp_pow(a, d, n);
  if (x === 1n || x === n - 1n) return true;
  for (let r = 1; r < s; r++) {
    x = (x * x) % n;
    if (x === n - 1n) return true;
    if (x === 1n) return false;
  }
  return false;
}

/** Jacobi symbol (a|n) by quadratic reciprocity; n > 0 odd. */
function jacobiSymbol(a: bigint, n: bigint): number {
  if (n <= 0n || (n & 1n) === 0n) {
    throw new Error('Jacobi symbol: n must be positive and odd');
  }
  a = ((a % n) + n) % n;
  let result = 1;
  while (a !== 0n) {
    while ((a & 1n) === 0n) {
      a >>= 1n;
      const nMod8 = Number(n & 7n);
      if (nMod8 === 3 || nMod8 === 5) result = -result;
    }
    [a, n] = [n, a];
    if ((a & 3n) === 3n && (n & 3n) === 3n) result = -result;
    a = a % n;
  }
  return n === 1n ? result : 0;
}

/** Strong Lucas probable prime test with Selfridge parameters. */
function strongLucas(n: bigint): boolean {
  let D = 5n;
  let sign = 1n;
  for (;;) {
    const jacobi = jacobiSymbol(D * sign, n);
    if (jacobi === 0) return n === D * sign || n === -(D * sign);
    if (jacobi === -1) {
      D = D * sign;
      break;
    }
    D += 2n;
    sign = -sign;
    if (D > 1000n) return true; /* n is a perfect square: caught elsewhere */
  }

  const P = 1n;
  const Q = (1n - D) / 4n;

  let d = n + 1n;
  let s = 0;
  while ((d & 1n) === 0n) {
    d >>= 1n;
    s++;
  }

  let U = 1n;
  let V = P;
  let Qk = Q;

  const bits: boolean[] = [];
  let temp = d;
  while (temp > 0n) {
    bits.push((temp & 1n) === 1n);
    temp >>= 1n;
  }

  for (let i = bits.length - 2; i >= 0; i--) {
    U = (U * V) % n;
    V = (V * V - 2n * Qk) % n;
    if (V < 0n) V += n;
    Qk = (Qk * Qk) % n;

    if (bits[i]) {
      const Unew = (P * U + V) % n;
      const Vnew = (D * U + P * V) % n;
      U = (Unew & 1n) === 0n ? Unew / 2n : (Unew + n) / 2n;
      V = (Vnew & 1n) === 0n ? Vnew / 2n : (Vnew + n) / 2n;
      U = ((U % n) + n) % n;
      V = ((V % n) + n) % n;
      Qk = (Qk * Q) % n;
      if (Qk < 0n) Qk += n;
    }
  }

  if (U === 0n) return true;
  for (let r = 0; r < s; r++) {
    if (V === 0n) return true;
    V = (V * V - 2n * Qk) % n;
    if (V < 0n) V += n;
    Qk = (Qk * Qk) % n;
  }
  return false;
}

/**
 * BPSW primality test (Baillie-Pomerance-Selfridge-Wagstaff).
 *
 * Reference: ifactor1.c:2514 (ifac_isprime -> BPSW_psp)
 */
export function isPrime(n: bigint): boolean {
  if (n < 2n) return false;
  if (n === 2n) return true;
  if ((n & 1n) === 0n) return false;

  for (const p of SMALL_PRIMES) {
    if (p * p > n) return true;
    if (n % p === 0n) return n === p;
  }
  if (!millerRabinWitness(n, 2n)) return false;
  return strongLucas(n);
}

// ============================================================================
// Pure power detection
// Reference: ifactor1.c:1957-2100 (is_357_power, is_kth_power, is_pth_power)
//            ispower.c:842 (Z_isanypower_101), ispower.c:1028 (isprimepower_i)
// ============================================================================

/** Squares mod 64: cheap filter used before extracting a square root. */
const SQUARES_MOD_64 = (() => {
  const t = new Uint8Array(64);
  for (let i = 0; i < 64; i++) t[(i * i) % 64] = 1;
  return t;
})();

/**
 * Is `x` a perfect square? Returns its square root or null.
 * Reference: ispower.c:Z_issquareall
 */
export function Z_issquareall(x: bigint): bigint | null {
  if (x < 0n) return null;
  if (x < 2n) return x;
  if (!SQUARES_MOD_64[Number(x & 63n)]) return null;
  const r = isqrt(x);
  return r * r === x ? r : null;
}

/**
 * Is `x` a 3rd, 5th or 7th power? `mask` has bit 0/1/2 set when a cube/5th/7th
 * power is still possible; the bit of a failed test is cleared.
 *
 * Returns `[exponent, root, newMask]`, exponent 0 on failure.
 *
 * Reference: ifactor1.c:1957 (is_357_power). PARI first runs a multistage
 * residue sieve mod 211*209*61*203 etc. and only then extracts a root; the
 * sieve is a pure optimization, we go straight to the exact root.
 */
export function is_357_power(x: bigint, mask: number): [number, bigint, number] {
  while (mask) {
    let b: number;
    let e: number;
    /* priority to higher powers -- ifactor1.c:1990 */
    if (mask & 4) {
      b = 4;
      e = 7;
    } else if (mask & 2) {
      b = 2;
      e = 5;
    } else {
      b = 1;
      e = 3;
    }
    const [y, exact] = Z_iroot(x, e);
    if (exact) return [e, y, mask];
    mask &= ~b;
  }
  return [0, 0n, mask];
}

/**
 * Is `x` an `n`-th power? Returns the `n`-th root or null.
 *
 * Reference: ifactor1.c:2010 (is_kth_power). Modular checks with primes
 * q = 1 (mod n) first, then an exact root.
 */
export function is_kth_power(x: bigint, n: number): bigint | null {
  /* number of modular checks, ifactor1.c:2027 */
  let j: number;
  if (n < 16) j = 5;
  else if (n < 32) j = 4;
  else if (n < 101) j = 3;
  else if (n < 1001) j = 2;
  else if (n < 17886697) j = 1;
  else j = 0;

  const N = BigInt(n);
  /* primes q = 1 mod n, starting at the smallest prime >= n */
  let q = n % 2 === 1 ? 2 * n + 1 : n + 1;
  const step = n % 2 === 1 ? 2 * n : n;
  for (; j > 0; j--) {
    while (q < Number.MAX_SAFE_INTEGER && !isPrime(BigInt(q))) q += step;
    if (q >= Number.MAX_SAFE_INTEGER) break;
    const Q = BigInt(q);
    const residue = x % Q;
    if (residue === 0n) {
      /* Z_lval(x,q) % n != 0 => not an n-th power */
      const [v] = lvalrem(x, Q);
      if (v % n) return null;
    } else if (Fp_pow(residue, (Q - 1n) / N, Q) !== 1n) {
      return null;
    }
    q += step;
  }

  const [y, exact] = Z_iroot(x, n);
  return exact ? y : null;
}

/**
 * Is `x` a p^i-th power for some prime p >= the current position of `T`?
 * Stops when `log2(x)/p < cutoffbits`, since such a base would have been found
 * by trial division.
 *
 * Returns `[v, root]` with `x = root^v` and `v` a power of p (v = 1: failure).
 *
 * Reference: ifactor1.c:2073 (is_pth_power)
 */
export function is_pth_power(x: bigint, T: ForprimeT, cutoffbits: number): [number, bigint] {
  let size = x.toString(2).length - 1; /* expi(x), not +1 */
  for (;;) {
    const p = T.next();
    if (p === null || size / p < cutoffbits) break;
    let v = 1;
    let y = x;
    for (;;) {
      const r = is_kth_power(y, p);
      if (r === null) break;
      v *= p;
      y = r;
      size = y.toString(2).length - 1;
    }
    if (v > 1) return [v, y];
  }
  return [1, x];
}

/**
 * Largest `k` with `x = y^k`, when every prime divisor of `x` is > 102.
 * Returns `[k, y]`.
 *
 * Reference: ispower.c:842 (Z_isanypower_101). PARI switches to floating point
 * logarithms for the large exponents; we keep exact integer roots.
 */
export function Z_isanypower_101(x: bigint): [number, bigint] {
  let k = 1;
  let y = x;
  for (;;) {
    const s = Z_issquareall(y);
    if (s === null) break;
    k <<= 1;
    y = s;
  }
  let mask = 7;
  for (;;) {
    const [e, r, m] = is_357_power(y, mask);
    mask = m;
    if (!e) break;
    k *= e;
    y = r;
  }
  /* every prime divisor is >= 103, so y = z^e forces 103^e <= y */
  const LOG2_103 = 6.6865; /* lower bound for log_2(103), ispower.c:845 */
  for (;;) {
    const e2 = Math.floor((y.toString(2).length - 1 + 1) / LOG2_103);
    if (e2 < 11) break;
    const T = forprime_init(11, e2);
    const [v, r] = is_pth_power(y, T, 1);
    if (v === 1) break;
    k *= v;
    y = r;
  }
  return [k, y];
}

/**
 * Largest `k >= 2` with `x = y^k`, or 0 when `x` is not a perfect power.
 * Returns `[k, y]` (k = 0 and y = x when x is not a perfect power).
 *
 * Reference: ispower.c:908 (Z_isanypower_aux) / ispower.c:987 (Z_isanypower)
 * Only |x| is considered; callers deal with the sign.
 */
export function Z_isanypower(x: bigint): [number, bigint] {
  if (x < 0n) x = -x;
  if (x < 2n) return [0, x];

  let k = 1;
  let e = 0; /* gcd of the valuations found by trial division */
  const P: bigint[] = [];
  const E: number[] = [];
  let rest = x;
  let done = false;
  for (const p of TINY_PRIMES) {
    const [v, q] = lvalrem(rest, p);
    if (v) {
      rest = q;
      P.push(p);
      E.push(v);
      e = e ? Number(gcdii(BigInt(e), BigInt(v))) : v;
      if (e === 1) {
        done = true;
        break;
      }
    }
    if (rest < p * p) {
      /* rest is 1 or prime */
      if (rest === 1n) k = e;
      done = true;
      break;
    }
  }

  if (!done) {
    if (e) {
      /* result divides e; strip powers of the primes dividing e */
      let e2 = e;
      let v2 = 0;
      while (e2 % 2 === 0) {
        e2 /= 2;
        v2++;
      }
      for (let i = 0; i < v2; i++) {
        const s = Z_issquareall(rest);
        if (s === null) break;
        k <<= 1;
        rest = s;
      }
      let mask = 0;
      const strip = (q: number, bit: number) => {
        let c = 0;
        while (e2 % q === 0) {
          e2 /= q;
          c++;
        }
        if (c) mask |= bit;
        return c;
      };
      let v3 = strip(3, 1);
      let v5 = strip(5, 2);
      let v7 = strip(7, 4);
      for (;;) {
        const [ex, r, m] = is_357_power(rest, mask);
        mask = m;
        if (!ex) break;
        rest = r;
        if (ex === 3) {
          k *= 3;
          if (--v3 === 0) mask &= ~1;
        } else if (ex === 5) {
          k *= 5;
          if (--v5 === 0) mask &= ~2;
        } else {
          k *= 7;
          if (--v7 === 0) mask &= ~4;
        }
      }
      /* split_exponent(e2, &rest): ispower.c:818 */
      if (e2 > 1) {
        for (const [q, mult] of factorSmallExponent(e2)) {
          for (let j = 0; j < mult; j++) {
            const r = is_kth_power(rest, q);
            if (r === null) break;
            k *= q;
            rest = r;
          }
        }
      }
    } else {
      const [k101, y] = Z_isanypower_101(rest);
      k = k101;
      rest = y;
    }
  }

  if (k === 1) return [0, x];
  /* add back the missing small factors -- ispower.c:966 */
  let y = rest;
  for (let i = 0; i < P.length; i++) {
    y *= P[i] ** BigInt(Math.floor(E[i] / k));
  }
  return [k, y];
}

/** Factor a small exponent (< 2^31) by trial division. */
function factorSmallExponent(e: number): Array<[number, number]> {
  const out: Array<[number, number]> = [];
  for (let p = 2; p * p <= e; p++) {
    let c = 0;
    while (e % p === 0) {
      e /= p;
      c++;
    }
    if (c) out.push([p, c]);
  }
  if (e > 1) out.push([e, 1]);
  return out;
}

/**
 * Is `n` a prime power `p^k` (k >= 1)? Returns `[p, k]` or null.
 *
 * NOTE: PARI keeps this in `ispower.c` (isprimepower_i, ispower.c:1028); we
 * host it here because parigp-ts has no `ispower.ts` yet.
 *
 * This never factors `n`: it strips perfect powers exactly and runs BPSW on
 * the base, which is what makes `is_prime_power(p*q)` cheap for huge p, q.
 */
export function isprimepower(n: bigint): [bigint, number] | null {
  if (n <= 0n) return null;
  if (n === 1n) return null;

  for (const p of TINY_PRIMES) {
    const [v, q] = lvalrem(n, p);
    if (v) return q === 1n ? [p, v] : null;
  }
  /* p | n => p >= 103 */
  const [k, base] = Z_isanypower_101(n);
  if (!isPrime(base)) return null;
  return [base, k];
}

// ============================================================================
// Trial division bound
// Reference: ifactor1.c:3137 (tridiv_boundu), ifactor1.c:3158 (tridiv_bound)
// ============================================================================

/**
 * Where to stop trial dividing. Must be >= 661.
 * Reference: ifactor1.c:3137-3168 (64-bit branch)
 */
export function tridiv_bound(n: bigint): number {
  const bits = n.toString(2).length;
  if (bits <= 64) {
    const e = bits - 1; /* expu(n) */
    if (e < 30) return 1 << 12;
    if (e < 34) return 1 << 13;
    if (e < 37) return 1 << 14;
    if (e < 42) return 1 << 15;
    if (e < 47) return 1 << 16;
    if (e < 56) return 1 << 17;
    /* ifactor1.c:3147 repeats the `e<56` test, so 1<<18 is unreachable there;
     * kept verbatim for fidelity. */
    if (e < 62) return 1 << 19;
    return 1 << 18;
  }
  if (bits <= 512) return (bits - 16) << 10;
  return 1 << 19; /* Rho is generally faster above this */
}

/** Product of the odd primes <= limit, cached (PARI: prodprimes()). */
const primorialCache = new Map<number, bigint>();
function oddPrimorial(limit: number): bigint {
  const cached = primorialCache.get(limit);
  if (cached !== undefined) return cached;
  const primes = primesUpTo(limit);
  /* product tree, to keep the cost quasi-linear */
  let level: bigint[] = [];
  for (const p of primes) if (p !== 2) level.push(BigInt(p));
  if (level.length === 0) level = [1n];
  while (level.length > 1) {
    const next: bigint[] = [];
    for (let i = 0; i < level.length; i += 2) {
      next.push(i + 1 < level.length ? level[i] * level[i + 1] : level[i]);
    }
    level = next;
  }
  primorialCache.set(limit, level[0]);
  return level[0];
}

// ============================================================================
// FACTORIZATION (Shanks' SQUFOF)
// Reference: ifactor1.c:1412 (squfof_ambig), ifactor1.c:1474 (squfof)
// ============================================================================

const SQUFOF_BLACKLIST_SZ = 64;

/** floor sqrt of a JS number < 2^53. */
function usqrt(x: number): number {
  let r = Math.floor(Math.sqrt(x));
  while (r * r > x) r--;
  while ((r + 1) * (r + 1) <= x) r++;
  return r;
}

function uissquare(x: number): number | null {
  if (x < 0) return null;
  const r = usqrt(x);
  return r * r === x ? r : null;
}

function ugcd(a: number, b: number): number {
  while (b) {
    const t = a % b;
    a = b;
    b = t;
  }
  return a;
}

/**
 * Walk back along the ambiguous cycle from a square form to the factor.
 * Reference: ifactor1.c:1412
 */
function squfof_ambig(a: number, B: number, dd: number, D: bigint): number {
  let q = Math.floor((dd + Math.floor(B / 2)) / a);
  const qa = q * a;
  let b = qa - B + qa;
  /* c = ((D - b^2)/4)/a -- needs exact arithmetic beyond 2^53 */
  let c = Number((D - BigInt(b) * BigInt(b)) / 4n / BigInt(a));
  const a0 = a;
  const b0 = b;
  let b1 = b;
  let cnt = 0;
  for (;;) {
    const c0 = c;
    let qcb: number;
    if (c0 > dd) q = 1;
    else q = Math.floor((dd + Math.floor(b / 2)) / c0);
    if (q === 1) {
      qcb = c0 - b;
      b = c0 + qcb;
      c = a - qcb;
    } else {
      const qc = q * c0;
      qcb = qc - b;
      b = qc + qcb;
      c = a - q * qcb;
    }
    a = c0;
    cnt++;
    if (b === b1) break;
    if (b === b0 && a === a0) return 0;
    b1 = b;
  }
  return a % 2 !== 0 ? a : a / 2;
}

/**
 * Shanks' square forms factorization. Assumes `(n,30) = 1`, `n` composite and
 * not a perfect square.
 *
 * Returns a list of factors whose product is `n`, or null on failure.
 *
 * Reference: ifactor1.c:1474
 *
 * @see Deviation: PARI's 64-bit build declines above 2^46 (ifactor1.c:1487)
 * because MPQS takes over there; MPQS is not ported, so we use the algorithm
 * up to its documented limit of 2^59 (ifactor1.c:1492). All arithmetic stays
 * exact (values are < 2^32 except the discriminant, which is a bigint).
 */
export function squfof(n: bigint): bigint[] | null {
  if (n >= 1n << 59n) return null;

  const nm4 = Number(n & 3n);
  let D1: bigint;
  let D2: bigint;
  let d2: number;
  let b1: number;
  let b2: number;
  let c1: number;
  let c2: number;
  let dd2: number;

  if (nm4 === 1) {
    /* n = 1 (mod 4): D1 = n, D2 = 5n */
    D1 = n;
    D2 = 5n * n;
    d2 = Number(isqrt(D2));
    dd2 = (d2 >> 1) + (d2 & 1);
    b2 = d2 % 2 === 1 ? d2 : d2 - 1; /* (d2-1)|1: largest odd <= d2 */
  } else {
    /* n = 3 (mod 4): D1 = 3n, D2 = 4n */
    D1 = 3n * n;
    D2 = 4n * n;
    dd2 = Number(isqrt(n));
    d2 = dd2 * 2;
    b2 = d2 - (d2 % 2); /* largest even <= d2 */
  }
  const d1 = Number(isqrt(D1));
  b1 = d1 % 2 === 1 ? d1 : d1 - 1; /* largest odd <= d1 */
  c1 = Number((D1 - BigInt(b1) * BigInt(b1)) / 4n);
  c2 = Number((D2 - BigInt(b2) * BigInt(b2)) / 4n);
  if (c1 === 0 || c2 === 0) return null; /* n or 3n/5n was a square */

  const L1 = usqrt(d1);
  const L2 = usqrt(d2);
  const dd1 = (d1 >> 1) + (d1 & 1);
  let a1 = 1;
  let a2 = 1;

  const blacklist1: number[] = [];
  const blacklist2: number[] = [];
  let act1 = true;
  let act2 = true;
  let cnt = 0;

  while (act1 || act2) {
    if (act1) {
      const c = c1;
      const q = c > dd1 ? 1 : Math.floor((dd1 + Math.floor(b1 / 2)) / c);
      if (q === 1) {
        const qcb = c - b1;
        b1 = c + qcb;
        c1 = a1 - qcb;
      } else {
        const qc = q * c;
        const qcb = qc - b1;
        b1 = qc + qcb;
        c1 = a1 - q * qcb;
      }
      a1 = c;
      if (a1 <= L1) {
        if (blacklist1.length >= SQUFOF_BLACKLIST_SZ) act1 = false;
        else blacklist1.push(a1);
      }
    }
    if (act2) {
      const c = c2;
      const q = c > dd2 ? 1 : Math.floor((dd2 + Math.floor(b2 / 2)) / c);
      if (q === 1) {
        const qcb = c - b2;
        b2 = c + qcb;
        c2 = a2 - qcb;
      } else {
        const qc = q * c;
        const qcb = qc - b2;
        b2 = qc + qcb;
        c2 = a2 - q * qcb;
      }
      a2 = c;
      if (a2 <= L2) {
        if (blacklist2.length >= SQUFOF_BLACKLIST_SZ) act2 = false;
        else blacklist2.push(a2);
      }
    }

    if (++cnt & 1) continue; /* odd iteration */

    if (act1 && a1 === 1) act1 = false; /* back to identity: drop it */
    if (act1) {
      const root = uissquare(a1);
      /* not blacklisted? (the root form must not be on the principal cycle) */
      if (root !== null && root > 0 && !(root <= L1 && blacklist1.includes(root))) {
        const g = ugcd(root, Math.abs(b1));
        if (g > 1) {
          /* imprimitive form: g^2 divides D1 hence n */
          const G = BigInt(g);
          if (n % (G * G) === 0n) return [G, G, n / (G * G)].filter((f) => f !== 1n);
        } else {
          let q = squfof_ambig(root, b1, dd1, D1);
          if (q > 3) {
            if (nm4 === 3 && q % 3 === 0) q /= 3;
            const Q = BigInt(q);
            if (Q > 1n && Q < n && n % Q === 0n) return [Q, n / Q];
          }
        }
      }
    }

    if (act2 && a2 === 1) act2 = false;
    if (act2) {
      const root = uissquare(a2);
      if (root !== null && root > 0 && !(root <= L2 && blacklist2.includes(root))) {
        const g = ugcd(root, Math.abs(b2));
        if (g > 1) {
          const G = BigInt(g);
          if (n % (G * G) === 0n) return [G, G, n / (G * G)].filter((f) => f !== 1n);
        } else {
          let q = squfof_ambig(root, b2, dd2, D2);
          if (q > 5) {
            if (nm4 === 1 && q % 5 === 0) q /= 5;
            const Q = BigInt(q);
            if (Q > 1n && Q < n && n % Q === 0n) return [Q, n / Q];
          }
        }
      }
    }
  }
  return null;
}

// ============================================================================
// FACTORIZATION (Pollard-Brent rho)
// Reference: ifactor1.c:1171 (one_iter), 1184 (pollardbrent_i), 1361 (pollardbrent)
// ============================================================================

/** x <- x^2 + delta mod n, P <- P * (x1 - x) mod n. Reference: ifactor1.c:1171 */
function one_iter(x: bigint, P: bigint, x1: bigint, n: bigint, delta: bigint): [bigint, bigint] {
  const xn = ((x * x) % n) + delta;
  let Pn = (P * (x1 - xn)) % n;
  if (Pn < 0n) Pn += n;
  return [xn, Pn];
}

/**
 * Pollard-Brent rho with Brent's cycle detection and batched gcds.
 *
 * Returns a list of factors whose product is `n` (2 or 3 entries, possibly
 * composite), or null when the round budget `c0` is exhausted.
 *
 * Reference: ifactor1.c:1184
 */
export function pollardbrent_i(
  n: bigint,
  size: number,
  c0: number,
  retries: number
): bigint[] | null {
  let c = c0 << 5; /* 2^5 iterations per round */

  for (;;) {
    /* 'random' choice of delta determined by n -- ifactor1.c:1206 */
    let delta: bigint;
    switch ((size + retries) & 7) {
      case 0:
        delta = 1n;
        break;
      case 1:
        delta = -1n;
        break;
      case 2:
        delta = 3n;
        break;
      case 3:
        delta = 5n;
        break;
      case 4:
        delta = -5n;
        break;
      case 5:
        delta = 7n;
        break;
      case 6:
        delta = 11n;
        break;
      default:
        delta = -11n;
        break;
    }

    let x = 2n;
    let P = 1n;
    let y = 2n;
    let x1 = 2n;
    let k = 1;
    let l = 1;
    let g = 1n;
    let gaveUp = false;

    for (;;) {
      [x, P] = one_iter(x, P, x1, n, delta);

      if ((--c & 0x1f) === 0) {
        g = gcdii(n, P);
        if (g !== 1n) break;
        if (c <= 0) {
          gaveUp = true;
          break;
        }
        P = 1n;
        y = x;
      }

      if (--k) continue;

      if (c & 0x1f) {
        g = gcdii(n, P);
        if (g !== 1n) break;
        P = 1n;
      }

      /* fast forward phase: l iterations without gcds */
      c -= l >> 1;
      if (c <= 0) {
        gaveUp = true;
        break;
      }
      c &= ~0x1f;

      x1 = x;
      k = l;
      l <<= 1;
      for (let k1 = k; k1; k1--) [x, P] = one_iter(x, P, x1, n, delta);
      y = x;
    }
    if (gaveUp) return null;

    /* fin: an accumulated gcd was > 1 -- ifactor1.c:1286 */
    let g1: bigint;
    if (g !== n) {
      if (isPrime(g)) return [g, n / g];
      g1 = g; /* known composite, work modulo g1 */
    } else {
      g1 = n;
    }

    x = y;
    for (;;) {
      /* backtrack until period recovered; must terminate */
      x = ((x * x) % g1) + delta;
      g = gcdii((((x1 - x) % g1) + g1) % g1, g1);
      if (g !== 1n) break;
    }

    if (g1 === n || g === g1) {
      if (g1 === n && g === g1) {
        /* out of luck: restart with another delta -- ifactor1.c:1324 */
        if (++retries >= 4) return null;
        continue;
      }
      /* half lucky: we split n; g may be composite */
      return [g, n / g];
    }
    /* g < g1 < n: we split g1 as well */
    return [g, g1 / g, n / g1];
  }
}

/**
 * Pollard-Brent rho with PARI's round budget.
 *
 * Reference: ifactor1.c:1361
 *
 * @see Deviation: PARI declines for n < 2^96 (ifactor1.c:1365) because MPQS
 * covers that range faster. MPQS is not ported, so we accept every size; the
 * factors returned are unaffected.
 */
export function pollardbrent(n: bigint): bigint[] | null {
  const tune = 14;
  const size = n.toString(2).length; /* expi(n) + 1 */
  let c0: number;
  if (size <= 301) {
    c0 = tune + size - 60 + ((size - 73) >> 1) * ((size - 70) >> 3) * ((size - 56) >> 4);
  } else {
    c0 = 49152; /* ECM is faster when it'd take longer */
  }
  /* PARI's formula goes negative below 60 bits, where it never calls rho
   * (MPQS covers that range); give those inputs the base budget instead. */
  if (c0 < tune) c0 = tune;
  return pollardbrent_i(n, size, c0, 0);
}

/**
 * `Z_pollardbrent`: rho with an explicit number of rounds and seed.
 * Reference: ifactor1.c:1379
 */
export function Z_pollardbrent(n: bigint, rounds: number, seed: number): bigint[] | null {
  return pollardbrent_i(n, n.toString(2).length, rounds, seed);
}

// ============================================================================
// FACTORIZATION (Lenstra-Montgomery ECM)
// Reference: ifactor1.c:379-1122 (curve arithmetic, ECM_loop, ellfacteur)
// ============================================================================

/* B1 schedules, ifactor1.c:302 (TB1) and ifactor1.c:322 (TB1_for_stage) */
const TB1: readonly number[] = [
  142, 172, 208, 252, 305, 370, 450, 545, 661, 801, 972, 1180, 1430, 1735, 2100, 2550, 3090, 3745,
  4540, 5505, 6675, 8090, 9810, 11900, 14420, 17490, 21200, 25700, 31160, 37780, 45810, 55550,
  67350, 81660, 99010, 120050, 145550, 176475, 213970, 259430, 314550, 381380, 462415, 560660,
  679780, 824220, 999340, 1211670, 1469110, 1781250, 2159700, 2618600, 3175000, 3849600, 4667500,
  5659200, 6861600, 8319500, 10087100, 12230300, 14828900, 17979600, 21799700, 26431500, 32047300,
  38856400, 47112200, 57122100, 69258800, 83974200, 101816200, 123449000, 149678200, 181480300,
  220039400, 266791100, 323476100, 392204900, 475536500, 576573500, 699077800, 847610500,
  1027701900, 1246057200, 1510806400, 1831806700, 2221009800, 2692906700, 3265067200, 3958794400,
  4799917500,
];
const TB1_for_stage: readonly number[] = [
  500, 520, 560, 620, 700, 800, 900, 1000, 1150, 1300, 1450, 1600, 1800, 2000, 2200, 2450, 2700,
  2950, 3250, 3600, 4000, 4400, 4850, 5300, 5800, 6400, 7100, 7850, 8700, 9600, 10600, 11700, 12900,
  14200, 15700, 17300, 19000, 21000, 23200, 25500, 28000, 31000, 34500, 38500, 43000, 48000, 53800,
  60400, 67750, 76000, 85300, 95700, 107400, 120500, 135400, 152000, 170800, 191800, 215400, 241800,
  271400, 304500, 341500, 383100, 429700, 481900, 540400, 606000, 679500, 761800, 854100, 957500,
  1073500,
];
const nbcmax = 64;

interface ECMPoint {
  x: bigint;
  y: bigint;
}

/** Result of a curve operation: a point, a factor of N, or 'infinity mod N'. */
type ECMOp = { t: 'p'; P: ECMPoint } | { t: 'g'; g: bigint } | { t: 'inf' };

/**
 * (Px,Py) + (Qx,Qy) on y^2 = x^3 + x + b over Z/NZ.
 * Reference: ifactor1.c:379 (FpE_add_i) + the inversion of ecm_elladd0
 */
function ecm_add(N: bigint, P: ECMPoint, Q: ECMPoint): ECMOp {
  const { inv, gcd } = invmod(P.x - Q.x, N);
  if (inv === null) return gcd === N ? { t: 'inf' } : { t: 'g', g: gcd };
  const slope = ((P.y - Q.y) * inv) % N;
  let x = (slope * slope - Q.x - P.x) % N;
  if (x < 0n) x += N;
  let y = (slope * (P.x - x) - P.y) % N;
  if (y < 0n) y += N;
  return { t: 'p', P: { x, y } };
}

/**
 * Doubling on y^2 = x^3 + x + b: L = (3x^2+1)/(2y).
 * Reference: ifactor1.c:516 (elldouble)
 */
function ecm_double(N: bigint, P: ECMPoint): ECMOp {
  const { inv, gcd } = invmod(P.y, N);
  if (inv === null) return gcd === N ? { t: 'inf' } : { t: 'g', g: gcd };
  let L = ((1n + 3n * ((P.x * P.x) % N)) * inv) % N;
  if (L !== 0n) L = (L & 1n) === 1n ? (L + N) / 2n : L / 2n; /* halve mod N */
  let x = (L * L - 2n * P.x) % N;
  if (x < 0n) x += N;
  let y = (L * (P.x - x) - P.y) % N;
  if (y < 0n) y += N;
  return { t: 'p', P: { x, y } };
}

/**
 * [k]P for k >= 2.
 *
 * @see Deviation: PARI uses Montgomery's PRAC addition chain (ifactor1.c:592
 * ellmult); we use the binary ladder. Both compute [k]P and both reveal a
 * factor whenever an intermediate denominator has a nontrivial gcd with N; the
 * particular intermediate points differ, so the curve on which a given factor
 * turns up may differ.
 */
function ecm_mul(N: bigint, P: ECMPoint, k: number): ECMOp {
  const bits = k.toString(2);
  let R = P;
  for (let i = 1; i < bits.length; i++) {
    const d = ecm_double(N, R);
    if (d.t !== 'p') return d;
    R = d.P;
    if (bits[i] === '1') {
      const a = ecm_add(N, R, P);
      if (a.t !== 'p') return a;
      R = a.P;
    }
  }
  return { t: 'p', P: R };
}

/**
 * One ECM round: `nbc` curves with the given B1.
 *
 * Reference: ifactor1.c:752 (ECM_loop)
 *
 * @see Deviation: PARI runs the nbc curves in parallel and batches their
 * modular inversions (Montgomery's trick); we run them one at a time. Stage 2
 * is PARI's "improved standard continuation" as well, but instead of
 * accumulating products of x-coordinate differences over a 210-helix with
 * baby/giant steps (ifactor1.c:885-1030) we step [p]Q additively through the
 * primes of (B1,B2] and let each addition's inversion report the factor.
 * Detection is equivalent: [p]Q vanishes mod a prime divisor exactly when the
 * corresponding denominator does; only the constant factor differs.
 *
 * A curve that reaches the point at infinity mod N (denominator divisible by
 * the whole of N) is abandoned; PARI, which is working on a whole batch, only
 * skips the current multiplier there (ifactor1.c:790) because the other curves
 * of the batch are still alive.
 */
function ECM_loop(N: bigint, nbc: number, seed: number, B1: number): bigint | null {
  const B2 = 110 * B1;
  const B2_rt = usqrt(B2);
  const nbc2 = nbc << 1;

  /* pick curves: X[i] = seed++ downwards, point i = (X[i], X[nbc+i]) */
  const coord: bigint[] = new Array(nbc2);
  for (let i = nbc2 - 1, s = seed; i >= 0; i--, s++) coord[i] = BigInt(s);

  for (let curve = 0; curve < nbc; curve++) {
    let Q: ECMPoint = { x: coord[curve] % N, y: coord[nbc + curve] % N };

    /* ---B1 PHASE--- (ifactor1.c:785) */
    let broke = false;
    for (let m = 1; m <= B2 / 2; m <<= 1) {
      const r = ecm_double(N, Q);
      if (r.t === 'g') return r.g;
      if (r.t === 'inf') {
        broke = true;
        break;
      }
      Q = r.P;
    }
    let p = 2;
    if (!broke) {
      /* p = 3,...,nextprime(B1); the loop conditions test the previous p, so
       * nextprime(B1) is included, exactly as in ifactor1.c:795-812 */
      const gen = forprime(3, B1 + 1000);
      const next = (): number | null => {
        const r = gen.next();
        return r.done ? null : r.value;
      };
      while (p < B1 && p <= B2_rt) {
        const np = next();
        if (np === null) break;
        p = np;
        const lim = B2 / p;
        for (let m = 1; m <= lim; m *= p) {
          const r = ecm_mul(N, Q, p);
          if (r.t === 'g') return r.g;
          if (r.t === 'inf') {
            broke = true;
            break;
          }
          Q = r.P;
        }
        if (broke) break;
      }
      /* primes larger than sqrt(B2) appear only to the 1st power */
      while (!broke && p < B1) {
        const np = next();
        if (np === null) break;
        p = np;
        const r = ecm_mul(N, Q, p);
        if (r.t === 'g') return r.g;
        if (r.t === 'inf') {
          broke = true;
          break;
        }
        Q = r.P;
      }
    }
    if (broke) continue; /* point at infinity mod N: this curve is spent */

    /* ---B2 PHASE--- (ifactor1.c:820) */
    /* precompute [2j]Q for the prime gaps we will need */
    const D: ECMPoint[] = [];
    const dbl = ecm_double(N, Q);
    if (dbl.t === 'g') return dbl.g;
    if (dbl.t === 'inf') continue;
    D[2] = dbl.P;
    let R: ECMPoint | null = null;
    let prev = 0;
    let failed = false;
    for (const q of forprime(p + 1, B2)) {
      if (R === null) {
        const r = ecm_mul(N, Q, q);
        if (r.t === 'g') return r.g;
        if (r.t === 'inf') {
          failed = true;
          break;
        }
        R = r.P;
      } else {
        const gap = q - prev;
        for (let g = 4; g <= gap; g += 2) {
          if (D[g] === undefined) {
            /* [4]Q = [2]([2]Q) is a doubling, not a generic addition */
            const s = g === 4 ? ecm_double(N, D[2]) : ecm_add(N, D[g - 2], D[2]);
            if (s.t === 'g') return s.g;
            if (s.t === 'inf') {
              failed = true;
              break;
            }
            D[g] = s.P;
          }
        }
        if (failed) break;
        const s = ecm_add(N, R, D[gap]);
        if (s.t === 'g') return s.g;
        if (s.t === 'inf') {
          failed = true;
          break;
        }
        R = s.P;
      }
      prev = q;
    }
    if (failed) continue;
  }
  return null;
}

/**
 * ECM driver.
 *
 * Reference: ifactor1.c:1038 (ellfacteur)
 *
 * @param insist - when false, decline small inputs and give up after `rep`
 *   rounds (PARI then calls MPQS); when true, keep escalating B1.
 * @param maxRounds - safety net for the insisting mode: PARI loops forever
 *   because MPQS has already dealt with everything ECM cannot reach, we have
 *   no MPQS so we must be able to report failure instead of hanging.
 */
export function ellfacteur(N: bigint, insist: boolean, maxRounds = 60): bigint | null {
  const size = N.toString(2).length; /* expi(N)+1 */
  let nbc: number;
  let dsn: number;
  let dsnmax: number;
  let rep = 0;
  let seed: number;

  if (insist) {
    const DSNMAX = TB1.length - 1;
    dsnmax = (size >> 2) - 10;
    if (dsnmax < 0) dsnmax = 0;
    else if (dsnmax > DSNMAX) dsnmax = DSNMAX;
    seed = 1 + (nbcmax << 7) * (size & 0xffff);
    dsn = (size >> 3) - 5;
    if (dsn < 0) dsn = 0;
    else if (dsn > 47) dsn = 47;
    nbc = dsn + (dsn >> 2) + 9;
    nbc &= ~3;
  } else {
    dsn = (size - 140) >> 3;
    if (dsn < 0) return null; /* too small, decline the task */
    if (dsn > 12) dsn = 12;
    rep = size <= 248 ? (size <= 176 ? (size - 124) >> 4 : (size - 148) >> 3) : (size - 224) >> 1;
    if (rep < 1) rep = 1;
    dsnmax = 72;
    seed = 1 + (nbcmax << 3) * (size & 0xf);
    nbc = ((size >> 3) << 2) - 80;
    if (nbc < 8) nbc = 8;
  }
  if (nbc > nbcmax) nbc = nbcmax;
  if (dsn > dsnmax) dsn = dsnmax;

  for (let round = 0; ; round++) {
    if (insist && round >= maxRounds) return null;
    const B1 = insist ? TB1[dsn] : TB1_for_stage[dsn];
    const g = ECM_loop(N, nbc, seed, B1);
    seed += nbc << 1;
    if (g && g > 1n && g < N) return g;
    if (dsn < dsnmax) {
      if (insist) dsn++;
      else {
        dsn += 2;
        if (dsn > dsnmax) dsn = dsnmax;
      }
    }
    if (!insist && !--rep) return null;
  }
}

// ============================================================================
// FACTORIZATION (master iteration)
// Reference: ifactor1.c:2786 (ifac_crack), 3006 (ifac_decomp)
// ============================================================================

/** Rounds given to Pollard-Brent in place of the missing MPQS stage. */
const MPQS_SUBSTITUTE_ROUNDS = 1 << 16;

/**
 * Options for `Z_factor`.
 *
 * `ecmRounds` bounds the final, "insisting" ECM stage. PARI loops there
 * forever (ifactor1.c:1131) because MPQS has already handled everything ECM
 * cannot reach; with no MPQS we must be able to stop and report failure
 * instead of hanging. Raise it to spend more time before giving up.
 */
export interface FactorOptions {
  ecmRounds?: number;
}

const DEFAULT_ECM_ROUNDS = 4;

/**
 * Split a composite `n` (no small prime factors, not a prime power) into a
 * list of factors whose product is `n`.
 *
 * Reference: ifactor1.c:2786 (ifac_crack). PARI's order is: pure powers,
 * SQUFOF, Pollard-Brent rho, ECM (non insisting), MPQS, ECM (insisting).
 *
 * @throws {NotImplementedError} when every ported stage fails; PARI would use
 * MPQS here (mpqs.c), which is not ported. We refuse to declare a composite
 * prime.
 */
function ifac_crack(n: bigint, ecmRounds: number): Array<[bigint, bigint]> {
  /* --- pure power stage (ifactor1.c:2810) ---
   * MPQS/rho cannot split p^k, so powers are peeled off first. */
  {
    let k = 1;
    let base = n;
    for (;;) {
      const s = Z_issquareall(base);
      if (s === null) break;
      k *= 2;
      base = s;
    }
    let mask = 7;
    for (;;) {
      const [e, r, m] = is_357_power(base, mask);
      mask = m;
      if (!e) break;
      k *= e;
      base = r;
    }
    /* cutoff at 15 bits: smaller bases were found by trial division */
    const T = forprime_init(11, 1 << 20);
    for (;;) {
      const [v, r] = is_pth_power(base, T, 15);
      if (v === 1) break;
      k *= v;
      base = r;
    }
    if (k > 1) {
      /* PARI records the common exponent in place (update_pow, ifactor1.c:2669)
       * instead of repeating the base k times, so the base is cracked once. */
      return [[base, BigInt(k)]];
    }
  }

  /* --- SQUFOF, then rho (ifactor1.c:2843) --- */
  let factors = squfof(n);
  if (!factors) factors = pollardbrent(n);
  /* --- first ECM stage, not insisting (ifactor1.c:2855) --- */
  if (!factors) {
    const g = ellfacteur(n, false);
    if (g) factors = [g, n / g];
  }
  /* --- MPQS stage (ifactor1.c:2860): not ported, extra rho rounds instead --- */
  if (!factors) {
    factors = Z_pollardbrent(n, MPQS_SUBSTITUTE_ROUNDS, 0);
  }
  /* --- final ECM stage, insisting (ifactor1.c:2865) --- */
  if (!factors) {
    const g = ellfacteur(n, true, ecmRounds);
    if (g) factors = [g, n / g];
  }
  if (!factors) {
    throw new NotImplementedError(
      'SAGE_NOT_IMPLEMENTED: mpqs (multiple polynomial quadratic sieve, ' +
        'reference/pari/src/basemath/mpqs.c) is not ported; trial division, ' +
        'SQUFOF, Pollard-Brent rho and ECM all failed to split the ' +
        `${n.toString(2).length}-bit composite ${n}`
    );
  }
  return factors.map((f) => [f, 1n] as [bigint, bigint]);
}

/**
 * Factor a composite `n` with no prime factor below the trial division bound.
 * Reference: ifactor1.c:3006 (ifac_decomp)
 */
function ifac_decomp(n: bigint, ecmRounds: number): Factorization {
  const found = new Map<string, bigint>(); /* prime -> exponent */
  const record = (p: bigint, e: bigint) => {
    const k = p.toString();
    found.set(k, (found.get(k) ?? 0n) + e);
  };

  const stack: Array<[bigint, bigint]> = [[n, 1n]];
  while (stack.length) {
    const [v, e] = stack.pop()!;
    if (v === 1n) continue;
    if (isPrime(v)) {
      record(v, e);
      continue;
    }
    for (const [f, m] of ifac_crack(v, ecmRounds)) {
      if (f > 1n) stack.push([f, e * m]);
    }
  }

  const out: Factorization = [];
  for (const [p, e] of found) out.push([BigInt(p), e]);
  out.sort((a, b) => (a[0] < b[0] ? -1 : a[0] > b[0] ? 1 : 0));
  return out;
}

/**
 * Z_factor - factor an integer completely.
 *
 * Reference: ifactor1.c:4514 (Z_factor) -> ifactor1.c:4279 (ifactor_sign)
 *
 * @param n - Integer to factor (must be nonzero)
 * @returns Array of `[prime, exponent]` pairs, sorted by prime; a leading
 *   `[-1n, 1n]` records a negative sign, exactly as PARI's `Z_factor`.
 * @throws {NotImplementedError} if the composite cannot be split by any ported
 *   algorithm (PARI would use MPQS).
 */
export function Z_factor(n: bigint, options?: FactorOptions): Factorization {
  const ecmRounds = options?.ecmRounds ?? DEFAULT_ECM_ROUNDS;
  if (n === 0n) {
    throw new Error('Z_factor: factorization of 0 is not defined');
  }

  const factors: Factorization = [];
  if (n < 0n) {
    factors.push([-1n, 1n]);
    n = -n;
  }
  if (n === 1n) return factors;

  /* trial division by 2 (ifactor1.c:4338) */
  if ((n & 1n) === 0n) {
    const [v, q] = lvalrem(n, 2n);
    factors.push([2n, BigInt(v)]);
    n = q;
    if (n === 1n) return factors;
  }

  /* PARI rounds the trial division bound down to the largest power of two for
   * which it has a cached primorial (ifactor1.c:3288 u_oddprimedivisors_gcd:
   * `b = expu(lim)-6`, `*pLIM = LIM[b]`); do the same. */
  const lim = 1 << (31 - Math.clz32(tridiv_bound(n)));
  {
    /* fast trial division: gcd with the primorial, then split the gcd
     * (ifactor1.c:3306 Z_oddprimedivisors_fast) */
    let g = gcdii(n, oddPrimorial(lim));
    if (g > 1n) {
      for (const p of forprime(3, lim)) {
        const P = BigInt(p);
        if (P * P > g) break;
        if (g % P === 0n) {
          g /= P;
          const [v, q] = lvalrem(n, P);
          factors.push([P, BigInt(v)]);
          n = q;
          if (g === 1n) break;
        }
      }
      if (g > 1n) {
        /* what is left of the gcd is prime */
        const [v, q] = lvalrem(n, g);
        factors.push([g, BigInt(v)]);
        n = q;
      }
    }
  }
  if (n === 1n) return factors;

  /* n has no prime factor <= lim: prime, or hand over to ifac */
  if (n <= BigInt(lim) * BigInt(lim) || isPrime(n)) {
    factors.push([n, 1n]);
    return factors;
  }
  for (const pe of ifac_decomp(n, ecmRounds)) factors.push(pe);

  /* keep the primes sorted (the sign, if any, stays in front) */
  const sign = factors.length && factors[0][0] === -1n ? factors.shift()! : null;
  factors.sort((a, b) => (a[0] < b[0] ? -1 : a[0] > b[0] ? 1 : 0));
  if (sign) factors.unshift(sign);
  return factors;
}

/**
 * factoru - factor a small unsigned integer.
 *
 * Reference: ifactor1.c:3522. PARI's version uses word arithmetic throughout;
 * the factorization returned is the same as `Z_factor`'s.
 */
export function factoru(n: bigint, options?: FactorOptions): Factorization {
  return Z_factor(n, options);
}

/**
 * Format a factorization as a string.
 *
 * @param f - Factorization to format
 * @returns String like "2^2 * 3 * 5^2"
 */
export function formatFactorization(f: Factorization): string {
  if (f.length === 0) return '1';
  return f.map(([p, e]) => (e === 1n ? `${p}` : `${p}^${e}`)).join(' * ');
}
