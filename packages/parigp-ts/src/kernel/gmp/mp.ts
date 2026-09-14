import { mpn_sqrtrem } from './sqrtrem.js';
/** Native real division: one-word, short Knuth and GMP quotient kernels.
 * Source: reference/pari/src/kernel/gmp/mp.c:636-844.
 * @see Deviation: Native PARI real division kernels
 */
import { type MpReal, negr, real_0_bit } from '../../qfb.js';
import { divru } from '../none/mp_indep.js';
const B = 1n << 64n,
  M = B - 1n,
  H = 1n << 63n;
const bits = (v: bigint) => v.toString(2).length;
function normalize(q: bigint, p: number, e: number, s: -1 | 1): MpReal {
  const high = q >> BigInt(p);
  let m = q & ((1n << BigInt(p)) - 1n);
  if (!high) e--;
  else if (high === 1n) m = (m >> 1n) | (1n << BigInt(p - 1));
  else {
    m = 1n << BigInt(p - 1);
    e++;
  }
  return { s, e, p, m };
}
function divrr_with_gmp(x: MpReal, y: MpReal): MpReal {
  const lx = x.p / 64,
    ly = y.p / 64,
    lw = Math.min(lx, ly),
    lly = Math.min(lw + 1, ly),
    lu = lw + lly,
    llx = Math.min(lu, lx);
  const den = y.m >> BigInt(y.p - lly * 64),
    num = (x.m >> BigInt(x.p - llx * 64)) << BigInt((lu - llx) * 64);
  let q = num / den;
  const r = num % den;
  if (r >> BigInt((lly - 1) * 64) > (den >> BigInt((lly - 1) * 64)) >> 1n) q++;
  return normalize(q, lw * 64, x.e - y.e, x.s === y.s ? 1 : -1);
}
function divrr_basecase(x: MpReal, y: MpReal): MpReal {
  const lx = x.p / 64 + 2,
    ly = y.p / 64 + 2,
    lr = Math.min(lx, ly),
    r = Array<bigint>(lr).fill(0n),
    yw = Array<bigint>(ly).fill(0n);
  const xw = (i: number) => (x.m >> BigInt(x.p - 64 * (i - 1))) & M;
  for (let i = 2; i < lr; i++) r[i - 1] = xw(i);
  r[lr - 1] = lx > ly ? xw(lr) : 0n;
  for (let i = 2; i < ly; i++) yw[i] = (y.m >> BigInt(y.p - 64 * (i - 1))) & M;
  const y0 = yw[2]!,
    y1 = yw[3]!;
  let carry = 0n,
    hi = 0n;
  const add = (a: bigint, b: bigint) => {
    const z = a + b;
    carry = z >> 64n;
    return z & M;
  };
  const addx = (a: bigint, b: bigint) => {
    const z = a + b + carry;
    carry = z >> 64n;
    return z & M;
  };
  const sub = (a: bigint, b: bigint) => {
    const z = a - b;
    carry = z < 0n ? 1n : 0n;
    return z & M;
  };
  const subx = (a: bigint, b: bigint) => {
    const z = a - b - carry;
    carry = z < 0n ? 1n : 0n;
    return z & M;
  };
  const mul = (a: bigint, b: bigint) => {
    const z = a * b;
    hi = z >> 64n;
    return z & M;
  };
  const addmul = (a: bigint, b: bigint) => {
    const z = a * b + hi;
    hi = z >> 64n;
    return z & M;
  };
  const increment = (end: number) => {
    let j = end;
    do {
      r[j] = (r[j]! + 1n) & M;
      if (r[j] || j === 0) break;
      j--;
    } while (true);
  };
  for (let i = 0; i < lr - 1; i++) {
    const at = (j: number) => r[i - 1 + j]!;
    const put = (j: number, z: bigint) => {
      r[i - 1 + j] = z;
    };
    let qp: bigint, k: bigint, j: number;
    if (at(1) === y0) {
      qp = M;
      k = add(y0, at(2));
    } else {
      if (at(1) > y0) {
        j = lr - i;
        put(j, sub(at(j), yw[j + 1]!));
        for (j--; j > 0; j--) put(j, subx(at(j), yw[j + 1]!));
        increment(i - 1);
      }
      hi = at(1);
      carry = 0n;
      const numerator = (hi << 64n) + at(2);
      qp = numerator / y0;
      hi = numerator % y0;
      k = hi;
    }
    j = lr - i + 1;
    if (!carry) {
      let k3 = mul(qp, y1),
        k4: bigint;
      if (j === 3) k4 = sub(hi, k);
      else {
        k3 = sub(k3, at(3));
        k4 = subx(hi, k);
      }
      while (!carry && k4) {
        qp = (qp - 1n) & M;
        k3 = sub(k3, y1);
        k4 = subx(k4, y0);
      }
    }
    if (j < ly) mul(qp, yw[j]!);
    else {
      hi = 0n;
      j = ly;
    }
    for (j--; j > 1; j--) {
      put(j, sub(at(j), addmul(qp, yw[j]!)));
      hi = (hi + carry) & M;
    }
    if (at(1) !== hi) {
      if (at(1) < hi) {
        qp = (qp - 1n) & M;
        j = lr - i - (lr - i >= ly ? 1 : 0);
        put(j, add(at(j), yw[j]!));
        for (j--; j > 1; j--) put(j, addx(at(j), yw[j]!));
      } else {
        put(1, (at(1) - hi) & M);
        while (at(1)) {
          qp = (qp + 1n) & M;
          if (!qp) increment(i - 1);
          j = lr - i - (lr - i >= ly ? 1 : 0);
          put(j, sub(at(j), yw[j]!));
          for (j--; j > 1; j--) put(j, subx(at(j), yw[j]!));
          put(1, (at(1) - carry) & M);
        }
      }
    }
    r[i] = qp;
  }
  if (r[lr - 1]! > y0 >> 1n) increment(lr - 2);
  let q = 0n;
  for (let i = 0; i < lr - 1; i++) q = (q << 64n) + r[i]!;
  return normalize(q, (lr - 2) * 64, x.e - y.e, x.s === y.s ? 1 : -1);
}
export function divrr(x: MpReal, y: MpReal): MpReal {
  if (!y.s) throw new RangeError('impossible inverse in divrr');
  if (!x.s) return real_0_bit(x.e - y.e);
  if (y.p === 64) {
    let e = x.e - y.e;
    let num = (x.m >> BigInt(Math.max(x.p - 128, 0))) << BigInt(Math.max(128 - x.p, 0));
    if (num >> 64n < y.m) e--;
    else num >>= 1n;
    let m = num / y.m;
    if (num % y.m > y.m >> 1n) m++;
    if (m === B) {
      m = H;
      e++;
    }
    return { s: x.s === y.s ? 1 : -1, e, p: 64, m };
  }
  return y.p >= 256 ? divrr_with_gmp(x, y) : divrr_basecase(x, y);
}
export function divri(x: MpReal, y: bigint): MpReal {
  if (y === 0n) throw new RangeError('impossible inverse in divri');
  const mag = y < 0n ? -y : y,
    by = bits(mag);
  if (!x.s) return real_0_bit(x.e - by + 1);
  if (by <= 63) {
    const z = divru(x, mag);
    return y < 0n ? negr(z) : z;
  }
  const ly = Math.ceil(by / 64),
    lly = Math.min(x.p / 64 + 1, ly),
    shift = ly * 64 - by;
  const den = (mag >> BigInt((ly - lly) * 64)) << BigInt(shift),
    num = x.m << BigInt(lly * 64);
  let q = num / den;
  const r = num % den;
  if (r >> BigInt((lly - 1) * 64) > (den >> BigInt((lly - 1) * 64)) >> 1n) q++;
  return normalize(q, x.p, x.e - by + 1, (y < 0n ? -1 : 1) === x.s ? 1 : -1);
}

/** Native PARI sqrtremi; the output pointer is represented by the second tuple entry. */
export function sqrtremi(a: bigint): [bigint, bigint] {
  return mpn_sqrtrem(a < 0n ? -a : a);
}

/** Native GMP real square root, including the even-exponent guard correction.
 * @see Deviation: Native PARI square-root results and GMP dependency
 */
export function sqrtr_abs(x: MpReal): MpReal {
  const p = x.p,
    e = Math.floor(x.e / 2),
    odd = x.e % 2 !== 0;
  const n = x.m << BigInt(p + (odd ? 0 : 127));
  const [root, remainder] = sqrtremi(n);
  let m: bigint;
  if (odd) m = root + (remainder > root ? 1n : 0n);
  else {
    const guard = root & M;
    m = root >> 64n;
    if (guard & H || (guard === H - 1n && remainder > m)) m++;
  }
  return { s: 1, e, m, p };
}

/** PARI mantissa2nr: signed mantissa times 2^n, truncated toward zero. */
export function mantissa2nr(x: MpReal, n: number): bigint {
  if (!x.s) return 0n;
  const magnitude = n >= 0 ? x.m << BigInt(n) : x.m >> BigInt(-n);
  return x.s < 0 ? -magnitude : magnitude;
}
