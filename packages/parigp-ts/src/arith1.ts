import { unsignedPower } from './_scalar_power.js';
import { ZV_prod } from './ZV.js';
import { PariError } from './errors.js';
/** PARI integer Hilbert symbols, arith1.c:hilbert/hilbertii/mphilbertoo. */
import { Z_pvalrem } from './gen2.js';
import { kronecker, Fp_pow } from './ff.js';
import { Z_factor } from './ifactor.js';

/**
 * Native integer domain; p=0 denotes the real place.
 * @see Deviation: Hilbert symbol dependency domain
 */
export function hilbert(x: bigint, y: bigint, p: bigint = 0n): number {
  if (p === 0n) return x === 0n || y === 0n ? 0 : x < 0n && y < 0n ? -1 : 1;
  if (p <= 1n) {
    const error = new Error(`not a prime number in hilbertii: ${p}`);
    error.name = 'PariError';
    throw error;
  }
  if (x === 0n || y === 0n) return 0;
  const [vx, ux] = Z_pvalrem(x, p),
    [vy, uy] = Z_pvalrem(y, p);
  const oddX = (vx & 1) !== 0,
    oddY = (vy & 1) !== 0;
  x = ux;
  y = uy;
  const mod4 = (n: bigint) => n & 3n;
  const omega = (n: bigint) => {
    const r = n & 7n;
    return r === 3n || r === 5n;
  };
  let z: number;
  if (p === 2n) {
    z = mod4(x) === 3n && mod4(y) === 3n ? -1 : 1;
    if (oddX && omega(y)) z = -z;
    if (oddY && omega(x)) z = -z;
  } else {
    z = oddX && oddY && mod4(p) === 3n ? -1 : 1;
    if (oddX && kronecker(y, p) < 0) z = -z;
    if (oddY && kronecker(x, p) < 0) z = -z;
  }
  return z;
}

/** Native unsigned interval product, pairing endpoints for the balanced branch. */
export function mulu_interval_step(a: bigint, b: bigint, step: bigint): bigint {
  if (a < 0n || b >= 1n << 64n || b < a || step <= 0n || step >= 1n << 64n)
    throw new RangeError(
      'mulu_interval_step requires an ordered unsigned interval and positive word step'
    );
  if (a === 0n) return 0n;
  const count = 1n + (b - a) / step;
  if (count < 61n) {
    let value = a;
    for (let k = a + step; k <= b; k += step) value *= k;
    return value;
  }
  b -= (b - a) % step;
  const sum = a + b,
    values: bigint[] = [];
  for (let k = a; ; k += step) {
    const other = sum - k;
    if (other < k) break;
    values.push(other === k ? k : k * other);
    if (other === k) break;
  }
  return ZV_prod(values);
}

/** PARI mpfact: odd interval products raised to their binary-layer powers. */
export function mpfact(n: bigint): bigint {
  if (n < 0n) throw new PariError('domain error in factorial: argument < 0');
  if (n >= 1n << 63n) throw new RangeError('mpfact requires a signed-word argument');
  if (n <= 12n)
    return [
      1n,
      1n,
      2n,
      6n,
      24n,
      120n,
      720n,
      5040n,
      40320n,
      362880n,
      3628800n,
      39916800n,
      479001600n,
    ][Number(n)]!;
  const layers: bigint[] = [];
  for (let k = 1n; ; k++) {
    const upper = n >> (k - 1n);
    if (upper <= 2n) break;
    const lower = (1n + (n >> k)) | 1n;
    layers.push(mulu_interval_step(lower, upper, 2n) ** k);
  }
  let value = layers.pop()!;
  while (layers.length) value *= layers.pop()!;
  const popcount = n.toString(2).replaceAll('0', '').length;
  return value << (n - BigInt(popcount));
}

/** PARI arith1.c:2432: unsigned-word powering, including small-power shortcuts.
 * Exponents enter through the native unsigned 64-bit cast.
 */
export function Fp_powu(a: bigint, exponent: bigint, modulus: bigint): bigint {
  return unsignedPower(a, BigInt.asUintN(64, exponent), modulus);
}

/** PARI arith1.c:pgener_Zp/pgener_Zl, Conrey prime-power normalization. */
export function pgener_Zp(p: bigint): bigint {
  if (p === 2n) throw new PariError('domain error in pgener_Zl: p = 2');
  if (p === 40487n) return 10n;
  if (p <= 19n) return p === 2n ? 1n : p === 7n || p === 17n ? 3n : 2n;
  let q = p / 2n;
  while (q % 2n === 0n) q /= 2n;
  const exponents = Z_factor(q)
    .map(([l]) => p / 2n / l)
    .reverse();
  for (let x = 2n; ; x++) {
    if (kronecker(x, p) >= 0) continue;
    if (
      exponents.some((e) => {
        const t = Fp_pow(x, e, p);
        return t === 1n || t === p - 1n;
      })
    )
      continue;
    if (p >= 1n << 32n && Fp_pow(x, p - 1n, p * p) === 1n) continue;
    return x;
  }
}

