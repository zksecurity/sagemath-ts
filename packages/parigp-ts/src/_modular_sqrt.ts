/** Native square-root kernels from PARI basemath/arith1.c:717–1289.
 * @see Deviation: PARI modular square-root adapters
 */
import { gen_pow_fold } from './bb_group.js';
import { PariError } from './errors.js';
import { Fp_pow, kronecker } from './ff.js';
import { forprime, nextprime } from './ifactor.js';
import { vali } from './types.js';

const WORD_LIMIT = 1n << 64n;
const SIGNED_LIMIT = 1n << 63n;
const mod = (a: bigint, p: bigint): bigint => ((a % p) + p) % p;

function primeError(name: string, p: bigint): never {
  throw new PariError(`not a prime number in ${name}: ${p}.`);
}

/** The native search enumerates primes, starting at 3, through ULONG_MAX. */
function* smallPrimes(): Generator<bigint> {
  for (let lo = 3; lo < 1 << 20; lo += 1 << 16)
    for (const q of forprime(lo, lo + (1 << 16) - 1)) yield BigInt(q);
  let q = nextprime(BigInt(3 + (1 << 20)));
  for (; q < WORD_LIMIT; q = nextprime(q + 1n)) yield q;
}

export function nonsquare1_Fl(p: bigint): bigint {
  if ((p & 7n) !== 1n) return 2n;
  for (const q of [3n, 5n, 7n]) {
    const r = p % q;
    if (q === 3n ? r === 2n : q === 5n ? r === 2n || r === 3n : r !== 4n && r >= 3n) return q;
    if (r === 0n) primeError('Fl_nonsquare', p);
  }
  for (const q of forprime(11, 1967)) if (kronecker(BigInt(q), p) < 0) return BigInt(q);
  return primeError('Fl_nonsquare', p);
}

function nonsquare_Fp(p: bigint): bigint {
  if ((p & 3n) === 3n) return -1n;
  if ((p & 7n) === 5n) return 2n;
  for (const q of smallPrimes()) if (kronecker(q, p) < 0) return q;
  return primeError('Fp_sqrt [modulus]', p);
}

function Fp_rootsof1(l: bigint, p: bigint): bigint {
  const exponent = (p - 1n) / l;
  for (const q of smallPrimes()) {
    const z = Fp_pow(q, exponent, p);
    if (z !== 1n) return z;
  }
  return primeError('Fp_sqrt [modulus]', p);
}

function Fp_gausssum(D: bigint, p: bigint): bigint {
  const l = D < 0n ? -D : D;
  const z = Fp_rootsof1(l, p);
  let s = z;
  let x = z;
  for (let i = 2n; i < l; i++) {
    x *= z;
    const k = kronecker(i, l);
    if (k === 1) s += x;
    else if (k === -1) s -= x;
  }
  return s;
}

/** Zero is the native fall-through sentinel; null means a nonresidue. */
function Fp_sqrts(a: bigint, p: bigint): bigint | null {
  const v = Math.floor(vali(a) / 2);
  a >>= BigInt(2 * v);
  let r: bigint;
  switch (a) {
    case 1n:
      r = 1n;
      break;
    case -1n:
      if ((p & 3n) !== 1n) return null;
      r = Fp_pow(nonsquare_Fp(p), p >> 2n, p);
      break;
    case 2n:
    case -2n:
      if ((p & 7n) === 1n) {
        const z = Fp_pow(nonsquare_Fp(p), p >> 3n, p);
        r = mod(z * (1n + (a === 2n ? -z * z : z * z)), p);
      } else if ((p & 7n) === (a === 2n ? 7n : 3n)) r = Fp_pow(a, (p + 1n) >> 2n, p);
      else return null;
      break;
    case -3n: {
      if (p % 3n !== 1n) return null;
      const z = Fp_rootsof1(3n, p);
      r = mod(z - z * z, p);
      break;
    }
    case 5n:
    case 13n:
    case 17n:
    case 21n:
    case 29n:
    case 33n:
    case -7n:
    case -11n:
    case -15n:
    case -19n:
    case -23n:
      if (p % (a < 0n ? -a : a) !== 1n) return 0n;
      r = Fp_gausssum(a, p);
      break;
    default:
      return 0n;
  }
  // remii preserves the sign of a Gauss sum, unlike modular reduction.
  return (r << BigInt(v)) % p;
}

function sqrt_Cipolla(a: bigint, p: bigint): bigint | null {
  if (kronecker(a, p) < 0) return null;
  const half = p >> 1n;
  if (a > half) a -= p;
  let t = 1n;
  while (kronecker(t * t - a, p) >= 0) t++;
  const n = t * t - a;
  const square = ([u, v]: [bigint, bigint]): [bigint, bigint] => {
    const u2 = u * u;
    const v2 = v * v;
    return [mod(u2 + v2 * n, p), mod((v + u) ** 2n - u2 - v2, p)];
  };
  const fused = ([u, v]: [bigint, bigint]): [bigint, bigint] => {
    const d = u + t * v;
    const d2 = d * d;
    const b = (a * v) % p;
    return [mod(t * d2 - b * (u + d), p), mod(d2 - b * v, p)];
  };
  const result = gen_pow_fold<[bigint, bigint]>([t, 1n], half, square, fused);
  return mod(result[1] * a, p);
}

function tonelli(a: bigint, y: bigint, p: bigint, e: number, word: boolean): bigint | null {
  const q = (p - 1n) >> BigInt(e);
  let first = word ? Fp_pow(a, q >> 1n, p) : 0n;
  if (word && first === 0n) return 0n;
  if (y === 0n) {
    if (word) y = Fp_pow(nonsquare1_Fl(p), q, p);
    else {
      let k = 2n;
      for (; ; k++) {
        const symbol = kronecker(k, p);
        if (symbol < 0) break;
        if (symbol === 0) primeError('Fp_sqrt [modulus]', p);
      }
      y = Fp_pow(k, q, p);
    }
  }
  if (!word) first = Fp_pow(a, q >> 1n, p);
  let v = mod(a * first, p);
  let w = mod(v * first, p);
  while (w !== 1n) {
    let power = mod(w * w, p);
    let k = 1;
    for (; power !== 1n && k < e; k++) power = mod(power * power, p);
    if (k === e) return null;
    power = y;
    for (let i = 1; i < e - k; i++) power = mod(power * power, p);
    y = mod(power * power, p);
    e = k;
    w = mod(y * w, p);
    v = mod(v * power, p);
  }
  return v;
}

/** Reduced unsigned-word input, with null adapting PARI's ULONG_MAX sentinel. */
export function wordSquareRoot(a: bigint, p: bigint): bigint | null {
  if (p <= 0n || p >= WORD_LIMIT) throw new RangeError('modulus must be a positive word integer');
  if (a < 0n || a >= p) throw new RangeError('word square-root argument must be reduced');
  if (a === 0n) return 0n;
  const e = vali(p - 1n);
  if (e === 0) {
    if (p !== 2n) primeError('Fl_sqrt [modulus]', p);
    return a & 1n;
  }
  let r: bigint | null;
  if (e === 1) {
    r = Fp_pow(a, (p + 1n) >> 2n, p);
    if (mod(r * r, p) !== a) return null;
  } else r = tonelli(a, 0n, p, e, true);
  return r === null ? null : r > p - r ? p - r : r;
}

/** Fp_sqrt_i keeps native word dispatch, signed-small shortcuts and cost test. */
export function modularSquareRoot(a: bigint, y: bigint | null, p: bigint): bigint | null {
  if (a === 0n) return 0n;
  if (p === 0n) throw new PariError('impossible inverse in dvmdii: 0.');
  const magnitude = p < 0n ? -p : p;
  a = mod(a, magnitude);
  if (magnitude < WORD_LIMIT) return wordSquareRoot(a, magnitude);
  if (a === 0n) return 0n;
  const e = vali(p - 1n);
  if (e === 0) primeError('Fp_sqrt [modulus]', p);
  let r: bigint | null = 0n;
  if (e === 1) {
    r = Fp_pow(a, ((p - 1n) >> 2n) + 1n, p);
    if (mod(r * r, p) !== a) return null;
  } else {
    let small = a < SIGNED_LIMIT ? a : 0n;
    if (small === 0n && a - p > -SIGNED_LIMIT && a - p < SIGNED_LIMIT) small = a - p;
    if (small !== 0n) r = Fp_sqrts(small, p);
    if (r === 0n) {
      if (e === 2) {
        const a2 = mod(2n * a, p);
        const v = Fp_pow(a2, (p - 1n) >> 3n, p);
        const I = mod(a2 * v * v, p);
        r = mod(a * v * (I - 1n), p);
        if (mod(r * r, p) !== a) return null;
      } else if (e * (e - 1) > 20 + 8 * (magnitude.toString(2).length - 1)) r = sqrt_Cipolla(a, p);
      else r = tonelli(a, y ?? 0n, p, e, false);
    }
  }
  return r === null ? null : r > p - r ? p - r : r;
}
