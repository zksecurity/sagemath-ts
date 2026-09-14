import { mantissa2nr } from '../gmp/mp.js';
import { real_0_bit, type MpReal } from '../../qfb.js';

/** Native 64-bit PARI dbltor, including normalization of IEEE subnormals. */
export function dbltor(x: number): MpReal {
  if (x === 0) return real_0_bit(-1023);
  const view = new DataView(new ArrayBuffer(8));
  view.setFloat64(0, x, false);
  const raw = view.getBigUint64(0, false),
    high = 1n << 63n;
  let e = Number((raw >> 52n) & 2047n) - 1023;
  if (e === 1024) throw new PariError('overflow in dbltor [NaN or Infinity]');
  let m = (raw << 11n) & (high - 1n);
  if (e === -1023) {
    const shift = 64 - m.toString(2).length;
    e -= shift - 1;
    m <<= BigInt(shift);
  } else m |= high;
  return { s: raw & high ? -1 : 1, e, m, p: 64 };
}

/** Native PARI rtodbl; preserves its rounding, exponent cutoff and subnormal behavior. */
export function rtodbl(x: MpReal<number | bigint>): number {
  if (!x.s || x.e < -1023) return 0;
  // Only exponents in the binary64 range can reach output encoding; conversion
  // of a wider native exponent is used solely to trigger the overflow cutoff.
  let e = Number(x.e);
  const high = 1n << 63n;
  const word = x.p >= 64 ? x.m >> BigInt(x.p - 64) : x.m << BigInt(64 - x.p);
  let a = (word & (high - 1n)) + 1024n;
  if (a & high) {
    e++;
    a = 0n;
  }
  if (e >= 1023) throw new PariError('overflow in t_REAL->double conversion');
  let raw = (BigInt(e + 1023) << 52n) | (a >> 11n);
  if (x.s < 0) raw |= high;
  const view = new DataView(new ArrayBuffer(8));
  view.setBigUint64(0, raw, false);
  return view.getFloat64(0, false);
}

import { addrr, divrr, negr, itor, real_1, rtor, subir } from '../../qfb.js';
import { quadratic_prec_mask } from '../../Zp.js';
import { PariError } from '../../errors.js';

/** Native reciprocal, including the precision-doubling Newton branch.
 * @see Deviation: Native PARI allocated zero records
 */
export function invr(b: MpReal): MpReal {
  // Large allocated zeros have unspecified native payload words; reject them.
  if (!b.s && (b.p === 0 || b.p > 4800)) throw new PariError('impossible inverse in invr: 0');
  const basecase = (y: MpReal) => rtor(divrr(real_1(y.p + 64), y), y.p);
  if (b.p <= 4800) return basecase(b);
  let mask = quadratic_prec_mask(b.p / 64),
    words = 1;
  for (let i = 0; i < 6; i++) {
    words *= 2;
    if (mask & 1n) words--;
    mask >>= 1n;
  }
  const a: MpReal = { ...b, s: 1, e: 0 };
  let x = rtor(basecase(rtor(a, words * 64)), b.p);
  while (mask > 1n) {
    words *= 2;
    if (mask & 1n) words--;
    mask >>= 1n;
    const p = words * 64,
      tail = b.p - p;
    const current: MpReal = { ...x, p, m: x.m >> BigInt(tail) };
    const factor: MpReal = { ...a, p, m: a.m >> BigInt(tail) };
    const correction = mulrr(current, subir(1n, mulrr(factor, current)));
    const next = rtor(addrr(current, correction), p);
    x = { ...next, p: b.p, m: (next.m << BigInt(tail)) | (x.m & ((1n << BigInt(tail)) - 1n)) };
  }
  return { ...x, s: b.s, e: x.e - b.e };
}

const B = 1n << 64n,
  H = 1n << 63n;
const bits = (v: bigint) => v.toString(2).length;
/** Native unsigned-word quotient with its normalization guard. */
export function divru(x: MpReal, y: bigint): MpReal {
  if (y <= 0n || y >= B) throw new RangeError('divru requires a nonzero unsigned word');
  if (!x.s) return real_0_bit(x.e - bits(y) + 1);
  if (!(y & (y - 1n))) return { ...x, e: x.e - bits(y) + 1 };
  const top = x.m >> BigInt(x.p - 64);
  let q: bigint,
    guard: bigint,
    e = x.e;
  if (y <= top) {
    q = x.m / y;
    guard = ((x.m % y) << 64n) / y;
  } else {
    q = (x.m << 64n) / y;
    guard = (x.m << 64n) % y;
    e -= 64;
  }
  const shift = 64 - bits(q >> BigInt(x.p - 64));
  let m = (q << BigInt(shift)) | (guard >> BigInt(64 - shift));
  e -= shift;
  if ((guard << BigInt(shift)) & H) {
    m++;
    if (m >> BigInt(x.p)) {
      m = 1n << BigInt(x.p - 1);
      e++;
    }
  }
  return { s: x.s, e, p: x.p, m };
}
/** Native integer/real quotient, including divur's reciprocal dispatch. */
export function divir(x: bigint, y: MpReal): MpReal {
  if (y.p === 0) throw new PariError('impossible inverse in divir: 0');
  if (!x) return real_0_bit(-y.p - y.e);
  const mag = x < 0n ? -x : x;
  if (mag < B && y.p > 4800) {
    const z = invr(y),
      v = mag === 1n ? z : mulir(mag, z);
    return x < 0n ? negr(v) : v;
  }
  return rtor(divrr(itor(x, y.p + 64), y), y.p);
}

/** Native truncated real products (mp_indep.c:139-449).
 * @see Deviation: Native PARI real multiplication kernels
 */
const M = B - 1n;
function rounded(m: bigint, e: number, p: number, s: -1 | 1): MpReal {
  if (m >> BigInt(p)) {
    m >>= 1n;
    e++;
  }
  return { s, e, m, p };
}
function mulrrz_end(high: bigint, guard: bigint, p: number, e: number, s: -1 | 1): MpReal {
  if (high >> BigInt(p - 1)) e++;
  else {
    high = (high << 1n) | (guard >> 63n);
    guard = (guard << 1n) & M;
  }
  if (guard & H) high++;
  return rounded(high, e, p, s);
}
function mulrrz_int(x: MpReal, y: MpReal, p: number, s: -1 | 1): MpReal {
  const product = x.m * y.m,
    shift = x.p + y.p - p;
  return mulrrz_end(product >> BigInt(shift), (product >> BigInt(shift - 64)) & M, p, x.e + y.e, s);
}
function mulrrz_i(x: MpReal, y: MpReal, flag: boolean, s: -1 | 1): MpReal {
  const lz = x.p / 64 + 2,
    z = Array<bigint>(lz).fill(0n);
  const word = (a: MpReal, i: number) => (a.m >> BigInt(a.p - 64 * (i - 1))) & M;
  let hi = 0n,
    overflow = 0n;
  const mul = (a: bigint, b: bigint) => {
    const t = a * b;
    hi = t >> 64n;
    return t & M;
  };
  const addmul = (a: bigint, b: bigint) => {
    const t = a * b + hi;
    hi = t >> 64n;
    return t & M;
  };
  const add = (a: bigint, b: bigint) => {
    const t = a + b;
    overflow = t >> 64n;
    return t & M;
  };
  if (lz === 3) {
    let guard: bigint;
    if (flag) {
      mul(word(x, 2), word(y, 3));
      guard = addmul(word(x, 2), word(y, 2));
    } else guard = mul(word(x, 2), word(y, 2));
    return mulrrz_end(hi, guard, x.p, x.e + y.e, s);
  }
  let guard = 0n;
  if (flag) {
    mul(word(x, 2), word(y, lz));
    guard = hi;
  }
  const lzz = lz - 1;
  let p1 = word(x, lzz);
  if (p1) {
    mul(p1, word(y, 3));
    guard = add(addmul(p1, word(y, 2)), guard);
    z[lzz] = hi + overflow;
  }
  for (let j = lz - 2; j >= 3; j--) {
    p1 = word(x, j);
    if (p1) {
      mul(p1, word(y, lz + 2 - j));
      guard = add(addmul(p1, word(y, lz + 1 - j)), guard);
      for (let i = lzz; i > j; i--) {
        hi += overflow;
        z[i] = add(addmul(p1, word(y, i + 1 - j)), z[i]!);
      }
      z[j] = hi + overflow;
    }
  }
  p1 = word(x, 2);
  guard = add(mul(p1, word(y, lz - 1)), guard);
  for (let i = lzz; i > 2; i--) {
    hi += overflow;
    z[i] = add(addmul(p1, word(y, i - 1)), z[i]!);
  }
  z[2] = hi + overflow;
  let high = 0n;
  for (let i = 2; i < lz; i++) high = (high << 64n) | z[i]!;
  return mulrrz_end(high, guard, x.p, x.e + y.e, s);
}
export function sqrr(x: MpReal): MpReal {
  if (!x.s) return real_0_bit(2 * x.e);
  return x.p > 768 ? mulrrz_int(x, x, x.p, 1) : mulrrz_i(x, x, false, 1);
}
export function mulrr(x: MpReal, y: MpReal): MpReal {
  if (x === y) return sqrr(x);
  if (!x.s || !y.s) return real_0_bit(x.e + y.e);
  if (x.p > y.p) [x, y] = [y, x];
  const s = x.s === y.s ? 1 : -1,
    flag = x.p !== y.p;
  if (x.p <= 3520) return mulrrz_i(x, y, flag, s);
  const p = x.p + (flag ? 64 : 0),
    yp = { ...y, p, m: y.m >> BigInt(y.p - p) };
  return mulrrz_int(x, yp, x.p, s);
}
export function mulir(x: bigint, y: MpReal): MpReal {
  if (!x) return real_0_bit(y.p ? y.e - y.p : y.e < 0 ? 2 * y.e : 0);
  const mag = x < 0n ? -x : x,
    s = (x < 0n ? -1 : 1) === y.s ? 1 : -1;
  if (!y.s) return real_0_bit(bits(mag) - 1 + y.e);
  if (mag < B) {
    if (mag === 1n) return { ...y, s };
    const product = mag * y.m,
      k = bits(product) - y.p;
    let m = product >> BigInt(k);
    if ((product >> BigInt(k - 1)) & 1n) m++;
    return rounded(m, y.e + k, y.p, s);
  }
  const xp = Math.ceil(bits(mag) / 64) * 64,
    lx = xp / 64 + 2,
    lz = y.p / 64 + 2;
  if (lx < Math.floor(lz / 2) || (lx < lz && y.p > 3520))
    return mulrrz_int(itor(mag, xp), y, y.p, s);
  return mulrr(itor(x, y.p), y);
}

/** PARI trunc2nr, used by gtrunc2n: unlike truncr, permits precision loss. */
export function trunc2nr(x: MpReal, n: number): bigint {
  const e = x.e + n;
  return !x.s || e < 0 ? 0n : mantissa2nr(x, e - x.p + 1);
}
