/**
 * GMP 6.3.0 integer factorial and prime-sieve kernels, using exact native BigInt.
 * Sources: mpz/fac_ui.c, mpz/oddfac_1.c, mpz/prodlimbs.c, primesieve.c,
 * gen-fac.c and gen-sieve.c, https://ftp.gnu.org/gnu/gmp/gmp-6.3.0.tar.xz.
 * Copyright 1991, 1993-1995, 2000-2002, 2005, 2009-2012, 2015-2017,
 * 2020-2022 Free Software Foundation, Inc.
 * Adapted under GPL-2.0-or-later, matching this package.
 *
 * @see Deviation: GMP factorial kernels and native resource bounds
 */
const WORD_MAX = (1n << 64n) - 1n;
// GMP 6.3.0 mpn/arm64/applem1/gmp-mparam.h: FAC_ODD=0, FAC_DSC=252, TOOM22=26.
const DSC_THRESHOLD = 252n;
const oddFac = [1n];
const oddDouble = [1n];
const fac = [1n];
for (let n = 1n; n <= 25n; n++) {
  let odd = n;
  while ((odd & 1n) === 0n) odd >>= 1n;
  oddFac.push(oddFac[oddFac.length - 1]! * odd);
  if (n <= 20n) fac.push(fac[fac.length - 1]! * n);
}
for (let n = 3n; n <= 33n; n += 2n) oddDouble.push(oddDouble[oddDouble.length - 1]! * n);
function mpz_prodlimbs(factors: readonly bigint[], start = 0, count = factors.length): bigint {
  if (count < 26) {
    let result = 1n;
    for (let i = start; i < start + count; i++) result *= factors[i]!;
    return result;
  }
  const left = Math.floor(count / 2);
  return mpz_prodlimbs(factors, start, left) * mpz_prodlimbs(factors, start + left, count - left);
}
const bit_to_n = (bit: number) => 3 * bit + 5 - (bit % 2);
const n_to_bit = (n: number) => Math.floor((n - 5 + (n % 2 === 1 ? 1 : 0)) / 3);
// Generated from GMP gen-sieve.c's 1792-bit initial sieve and 5/11, 7/13 masks.
const INITIAL: readonly number[] = [
  1762821248, 848611808, 3299549660, 2510511646, 3093902182, 1255657113, 1921893675, 1704310490,
  2276511454, 3933052807, 3442636201, 1062642164, 1957128923, 4248324347, 2716726959, 3686403537,
  3525810597, 3469209982, 3144777046, 3941341117, 1482358003, 990820275, 2682219599, 3848526070,
  2757661436, 4267419563, 1005886333, 361623151, 3991325978, 3193600964, 3397105325, 3613891391,
  535771113, 3287706519, 969495549, 1870576883, 3526745072, 3584421084, 3585498683, 3975838511,
  3365889969, 3532586489, 1037283151, 3414129786, 4285215436, 4005484237, 1590667644, 3585963000,
  3148695799, 570277455, 4005035495, 1580573621, 2816195785, 3656121683, 788406134, 4288601775,
];
const MASK_5_11 = 204251606677392680049439193302149n;
const MASK_7_13 = 479677307237979300934035444814585430017860845601588234n;
/** Composite-bit mask for candidates 6k +/- 1, beginning with 5. Internal n >= 5.
 * A pair of Uint32 words represents each original 64-bit limb, low word first.
 */
export function gmp_primesieve(n: bigint): Uint32Array {
  if (n < 5n || n > WORD_MAX) throw new RangeError('outside GMP prime-sieve domain');
  const last = ((n - 5n) | 1n) / 3n;
  const words = (last / 64n + 1n) * 2n;
  if (words > BigInt(Number.MAX_SAFE_INTEGER)) throw new RangeError('Invalid typed array length');
  const sieve = new Uint32Array(Number(words));
  sieve.set(INITIAL.slice(0, sieve.length));
  const has = (bit: number) => (sieve[Math.floor(bit / 32)]! & (1 << (bit % 32))) !== 0;
  const mark = (bit: number) => {
    const word = Math.floor(bit / 32);
    sieve[word] = sieve[word]! | (1 << (bit % 32));
  };
  function block_resieve(start: number, end: number): void {
    for (let bit = start; bit < end; bit += 32) {
      const offset1 = BigInt(bit % 110),
        offset2 = BigInt(bit % 182);
      const mask1 = (MASK_5_11 >> offset1) | (MASK_5_11 << (110n - offset1));
      const mask2 = (MASK_7_13 >> offset2) | (MASK_7_13 << (182n - offset2));
      sieve[bit / 32] = Number((mask1 | mask2) & 0xffffffffn);
    }
    for (let bit = 4; ; bit++) {
      if (has(bit)) continue;
      const id = bit + 1,
        prime = 3 * id + 1 + (id % 2),
        step = 2 * prime;
      const square = id * (prime + 1) - 1 + (id % 2) * (id + 1);
      if (square >= end) break;
      const next = id * (id * 3 + 6) + (id % 2);
      for (let index of [square, next]) {
        if (index < start) index += step * Math.ceil((start - index) / step);
        for (; index < end; index += step) mark(index);
      }
    }
  }
  const size64 = sieve.length / 2,
    block = 2048;
  if (size64 > 28) {
    let off = size64 > 2 * block ? block + (size64 % block) : size64;
    block_resieve(28 * 64, off * 64);
    for (; off < size64; off += block) block_resieve(off * 64, (off + block) * 64);
  }
  const count = Number(last + 1n);
  for (let bit = count; bit < sieve.length * 32; bit++) mark(bit);
  return sieve;
}
function mpz_2multiswing_1(original: bigint, sieve: Uint32Array): bigint {
  const n = original & ~1n,
    max = WORD_MAX / (n - 1n),
    factors: bigint[] = [];
  let prod = original % 2n ? original : 1n;
  const store = (p: bigint, limit: bigint) => {
    if (prod > limit) {
      factors.push(prod);
      prod = p;
    } else prod *= p;
  };
  const swingPrime = (p: bigint) => {
    if (prod > max) {
      factors.push(prod);
      prod = 1n;
    }
    let q = n;
    do {
      q /= p;
      if (q % 2n) prod *= p;
    } while (q >= p);
  };
  swingPrime(3n);
  const shift = BigInt(Math.floor(n.toString(2).length / 2));
  const approx = (1n << (shift - 1n)) + ((n >> 1n) >> shift);
  const has = (i: number) => (sieve[Math.floor(i / 32)]! & (1 << (i % 32))) !== 0;
  let i = 0;
  for (const limit of [Number(approx), Number(n / 3n)]) {
    const last = n_to_bit(limit);
    for (; i <= last; i++)
      if (!has(i)) {
        const p = BigInt(bit_to_n(i));
        if (limit === Number(approx)) swingPrime(p);
        else if ((n / p) % 2n) store(p, max * 3n);
      }
  }
  for (i = n_to_bit(Number(n >> 1n)) + 1; i <= n_to_bit(Number(n)); i++)
    if (!has(i)) store(BigInt(bit_to_n(i)), max);
  factors.push(prod);
  return mpz_prodlimbs(factors);
}
/** Odd part of n!. With flag=1, skip the final squaring (n must be >= 252). */
export function mpz_oddfac_1(n: bigint, flag: 0 | 1 = 0): bigint {
  if ((flag !== 0 && flag !== 1) || (flag === 1 && n < DSC_THRESHOLD))
    throw new RangeError('outside GMP odd-factorial flag domain');
  if (n < 0n || n > WORD_MAX) throw new RangeError('outside unsigned long domain');
  if (n <= 25n) return oddFac[Number(n)]!;
  if (n <= 34n) return oddDouble[Number((n - 1n) >> 1n)]! * oddFac[Number(n >> 1n)]!;
  let steps = 0n,
    tn = n;
  for (; tn >= DSC_THRESHOLD; steps++) tn >>= 1n;
  const factors: bigint[] = [];
  let prod = 1n,
    max = WORD_MAX / (DSC_THRESHOLD * DSC_THRESHOLD);
  const store = (p: bigint) => {
    if (prod > max) {
      factors.push(prod);
      prod = p;
    } else prod *= p;
  };
  do {
    factors.push(oddDouble[16]!);
    let diff = (tn - 33n) & ~1n;
    if (diff & 2n) {
      store(33n + diff);
      diff -= 2n;
    }
    if (diff) {
      let factor = 35n * (33n + diff);
      do {
        store(factor);
        diff -= 4n;
        factor += 2n * diff;
      } while (diff);
    }
    max <<= 2n;
    tn >>= 1n;
  } while (tn > 34n);
  factors.push(prod, oddDouble[Number((tn - 1n) >> 1n)]!, oddFac[Number(tn >> 1n)]!);
  let result = mpz_prodlimbs(factors);
  if (steps) {
    const sieve = gmp_primesieve(n - 1n);
    do {
      steps--;
      result =
        (flag === 1 && steps === 0n ? result : result * result) *
        mpz_2multiswing_1(n >> steps, sieve);
    } while (steps);
  }
  return result;
}
/** n! for an unsigned 64-bit argument; physical result sizes remain runtime-limited. */
export function mpz_fac_ui(n: bigint): bigint {
  if (n < 0n || n > WORD_MAX) throw new RangeError('outside unsigned long domain');
  if (n <= 20n) return fac[Number(n)]!;
  let count = 0n;
  for (let k = n; k; k >>= 1n) count += k & 1n;
  return mpz_oddfac_1(n) << (n - count);
}
