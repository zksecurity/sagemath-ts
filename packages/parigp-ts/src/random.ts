/** PARI basemath/random.c, 64-bit XORGEN state and exact rejection sampling.
 * @see Deviation: PARI word random-state adapter
 */
import { PariError } from './errors.js';
const MASK = (1n << 64n) - 1n,
  WEYL = 0x61c8864680b583ebn;
const state = Array<bigint>(64).fill(0n);
let xorgen_w = 0n,
  xorgen_i = 63;
function block(): bigint {
  xorgen_i = (xorgen_i + 1) & 63;
  let t = state[xorgen_i]!,
    v = state[(xorgen_i + 11) & 63]!;
  t ^= (t << 33n) & MASK;
  t ^= t >> 26n;
  v ^= (v << 27n) & MASK;
  v ^= v >> 29n;
  return (state[xorgen_i] = t ^ v);
}
function init_xor4096i(v: bigint): void {
  const step = (v: bigint) => {
    v ^= (v << 10n) & MASK;
    v ^= v >> 15n;
    v ^= (v << 4n) & MASK;
    v ^= v >> 13n;
    return v;
  };
  for (let k = 64; k > 0; k--) v = step(v);
  xorgen_w = v;
  for (let k = 0; k < 64; k++) {
    v = step(v);
    xorgen_w = (xorgen_w + WEYL) & MASK;
    state[k] = (v + xorgen_w) & MASK;
  }
  xorgen_i = 63;
  for (let k = 256; k > 0; k--) block();
}
export function pari_init_rand(): void {
  init_xor4096i(1n);
}
pari_init_rand();
export function pari_rand(): bigint {
  const v = block();
  xorgen_w = (xorgen_w + WEYL) & MASK;
  return (v + (xorgen_w ^ (xorgen_w >> 27n))) & MASK;
}
export function setrand(seed: bigint): void {
  if (seed <= 0n) throw new PariError('domain error in setrand: n <= 0');
  if (seed <= MASK) {
    init_xor4096i(seed);
    return;
  }
  if (Math.ceil(seed.toString(2).length / 64) !== 66)
    throw new PariError('domain error in setrand: n != getrand()');
  for (let i = 0; i < 64; i++) {
    state[i] = seed & MASK;
    seed >>= 64n;
  }
  xorgen_w = seed & MASK;
  xorgen_i = Number((seed >> 64n) & 63n);
}
export function getrand(): bigint {
  let seed = BigInt(xorgen_i || 64);
  seed = (seed << 64n) | xorgen_w;
  for (let i = 63; i >= 0; i--) seed = (seed << 64n) | state[i]!;
  return seed;
}
/** Native long return type is signed, including the k=64 case.
 * @see Deviation: PARI word random-state adapter
 */
export function random_bits(k: number): bigint {
  if (!Number.isSafeInteger(k) || k < 1 || k > 64)
    throw new RangeError('bits must be between 1 and 64');
  return BigInt.asIntN(64, pari_rand() >> BigInt(64 - k));
}
export function random_Fl(n: bigint): bigint {
  if (n <= 0n || n > MASK) throw new RangeError('limit must be a positive unsigned word');
  if (n === 1n) return 0n;
  const bits = n.toString(2).length;
  if ((n & (n - 1n)) === 0n) return pari_rand() >> BigInt(65 - bits);
  for (;;) {
    const d = pari_rand() >> BigInt(64 - bits);
    if (d < n) return d;
  }
}
export function randomi(n: bigint): bigint {
  if (n <= 0n) throw new RangeError('limit must be positive');
  if (n <= MASK) return random_Fl(n);
  let bits = n.toString(2).length;
  if ((n & (n - 1n)) === 0n) bits--;
  const words = Math.ceil(bits / 64),
    shift = BigInt(words * 64 - bits);
  for (;;) {
    let value = 0n;
    for (let i = 0; i < words - 1; i++) value |= pari_rand() << BigInt(64 * i);
    value |= (pari_rand() >> shift) << BigInt(64 * (words - 1));
    if (value < n) return value;
  }
}
/** Packed binary coefficients, bit i being the coefficient of x^i. */
export function random_F2x(d: number): bigint {
  if (!Number.isSafeInteger(d) || d < 0) throw new RangeError('length must be nonnegative');
  let value = 0n;
  for (let i = 0; i < Math.ceil(d / 64); i++) value |= pari_rand() << BigInt(64 * i);
  return d ? value & ((1n << BigInt(d)) - 1n) : 0n;
}
/** Native t_VECSMALL entries are signed machine words. */
export function random_zv(n: number): bigint[] {
  if (!Number.isSafeInteger(n) || n < 0) throw new RangeError('length must be nonnegative');
  return Array.from({ length: n }, () => BigInt.asIntN(64, pari_rand()));
}
