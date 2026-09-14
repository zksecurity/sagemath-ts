/** Integer-vector product from PARI basemath/ZV.c. */
import { gen_product } from './bb_group.js';

export function ZV_prod(values: readonly bigint[]): bigint {
  if (values.length <= 6) {
    let result = 1n;
    for (const value of values) result *= value;
    return result;
  }
  return gen_product(values, (a, b) => a * b);
}

/** PARI ZV.c integer matrix multiplication, with the original dispatch thresholds.
 * @see Deviation: PARI matrix product representation
 */
import { Fp_inv } from './ff.js';
import { nextprime } from './ifactor.js';
import { Flm_mul } from './FpV.js';
import {
  classical,
  fromColumns,
  toColumns,
  winograd,
  reduceMod,
  type Rows,
} from './_matrix_mul.js';

const integerLength = (x: bigint): number =>
  x === 0n ? 2 : 2 + Math.ceil((x < 0n ? -x : x).toString(2).length / 64);
const maxLength = (a: Rows, rows: number, cols: number): number =>
  a
    .slice(0, rows)
    .reduce((s, r) => r.slice(0, cols).reduce((t, x) => Math.max(t, integerLength(x)), s), 2);
const modularPrimes: bigint[] = [];
function primeAt(i: number): bigint {
  while (modularPrimes.length <= i)
    modularPrimes.push(
      nextprime(
        modularPrimes.length ? modularPrimes[modularPrimes.length - 1]! + 1n : (1n << 63n) + 1n
      )
    );
  return modularPrimes[i]!;
}

/** Native 2x2 kernel uses seven products only when every entry has at least 14 limbs. */
function mul2(a: Rows, b: Rows): Rows {
  if (
    [
      ...a.slice(0, 2).flatMap((r) => r.slice(0, 2)),
      ...b.slice(0, 2).flatMap((r) => r.slice(0, 2)),
    ].some((x) => integerLength(x) < 16)
  )
    return classical(a, b, 2, 2, 2);
  const [a11, a12] = a[0]!,
    [a21, a22] = a[1]!,
    [b11, b12] = b[0]!,
    [b21, b22] = b[1]!;
  const m1 = (a11! + a22!) * (b11! + b22!),
    m2 = (a21! + a22!) * b11!,
    m3 = a11! * (b12! - b22!),
    m4 = a22! * (b21! - b11!);
  const m5 = (a11! + a12!) * b22!,
    m6 = (a21! - a11!) * (b11! + b12!),
    m7 = (a12! - a22!) * (b21! + b22!);
  return [
    [m1 + m4 + m7 - m5, m3 + m5],
    [m2 + m4, m1 - m2 + m3 + m6],
  ];
}

/** Modular reconstruction uses native high-bit primes and a balanced CRT tree. */
function modularProduct(
  a: Rows,
  b: Rows,
  m: number,
  n: number,
  k: number,
  sx: number,
  sy: number
): Rows {
  if (sx === 2 || sy === 2) return Array.from({ length: m }, () => new Array<bigint>(k).fill(0n));
  const bound = 2 + (sx + sy - 4) * 64 + Math.floor(Math.log2(n));
  let modulus = 1n;
  const primes: bigint[] = [];
  while (modulus.toString(2).length - 1 < bound) {
    const count = Math.floor((bound - (modulus.toString(2).length - 1)) / 63) + 1;
    for (let i = 0; i < count; i++) {
      const p = primeAt(primes.length);
      primes.push(p);
      modulus *= p;
    }
  }
  type Residue = { value: Rows; modulus: bigint };
  const combine = (lo: number, hi: number): Residue => {
    if (hi - lo === 1) {
      const p = primes[lo]!,
        red = (v: Rows) => v.map((r) => r.map((z) => reduceMod(z, p)));
      return {
        value: fromColumns(Flm_mul(toColumns(red(a), n), toColumns(red(b), k), p)),
        modulus: p,
      };
    }
    const mid = Math.floor((lo + hi) / 2),
      left = combine(lo, mid),
      right = combine(mid, hi);
    const q = left.modulus,
      p = right.modulus,
      pq = p * q,
      inv = Fp_inv(q % p, p);
    return {
      value: left.value.map((r, i) =>
        r.map((v, j) => {
          const z = v + q * reduceMod((right.value[i]![j]! - v) * inv, p);
          return z > pq / 2n ? z - pq : z;
        })
      ),
      modulus: pq,
    };
  };
  const result = combine(0, primes.length);
  return result.value.map((r) => r.map((v) => (v > result.modulus / 2n ? v - result.modulus : v)));
}

function mul(a: Rows, b: Rows, m: number, n: number, k: number): Rows {
  if (!m || !n || !k) return classical(a, b, m, n, k);
  if (m === 2 && n === 2 && k === 2) return mul2(a, b);
  const sx = maxLength(a, m, n),
    sy = maxLength(b, n, k);
  if (Math.min(m + 1, n + 1, k + 1) > 70 && sx <= 10 * sy && sy <= 10 * sx)
    return modularProduct(a, b, m, n, k, sx, sy);
  const s = Math.min(sx, sy),
    cutoff = s > 60 ? 2 : s > 25 ? 4 : s > 15 ? 8 : s > 8 ? 16 : 32;
  return Math.min(m + 1, n + 1, k + 1) <= cutoff
    ? classical(a, b, m, n, k)
    : winograd(a, b, m, n, k, mul);
}

/** Columns/rows have unused slot zero, matching the existing buch ZM_mul API. */
export function ZM_mul(x: bigint[][], y: bigint[][]): bigint[][] {
  const a = fromColumns(x),
    b = fromColumns(y),
    k = Math.max(0, y.length - 1);
  return toColumns(mul(a, b, a.length, Math.max(0, x.length - 1), k), k);
}

/** PARI ZV.c:1540 assumes that every signed-word product fits.
 * Keep exact intermediates until the existing Number return. Reject overflow
 * before the native unchecked multiplication would leave its defined domain.
 * @see Deviation: PARI permutation and word-vector adapters
 */
export function zv_prod(v: readonly number[]): number {
  let product = 1n;
  for (const value of v) {
    product *= BigInt(value);
    if (product < -(1n << 63n) || product >= 1n << 63n)
      throw new RangeError('zv_prod requires a product within the signed-word range');
  }
  return Number(product);
}
