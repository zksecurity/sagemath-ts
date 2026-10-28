/** Named algorithms, not a claim to emulate any library's key generator. */
import { PythonRandom } from '../../../sagemath-ts/src/misc/randstate.js';

export const WIDTHS = {
  lrand48: 4,
  'drand48-state': 6,
  'glibc-rand': 4,
  'ansi-lcg31': 4,
  'python-mt19937': 4,
  'java-random': 4,
  'msvc-rand': 4,
  'cpp-mt19937': 4,
} as const;
export type Rng = keyof typeof WIDTHS;
export const RNGS = Object.keys(WIDTHS) as Rng[];
export type Word = [bigint, number];
export type WordStream = Iterator<Word>;

/** Exact integer words with explicit serialized widths; seed truncation is algorithm-specific. */
export function* words(name: Rng, seed: bigint): Generator<Word> {
  if (typeof seed !== 'bigint' || seed < 0n) throw new Error('Seed must be a nonnegative integer');
  if (name === 'python-mt19937') {
    const rng = new PythonRandom(seed);
    while (true) yield [rng.getrandbits(32), 4];
  } else if (name === 'glibc-rand') {
    const low = Number(seed & 0xffffffffn) || 1;
    const state = new Uint32Array(31);
    state[0] = low;
    let word = low | 0;
    for (let i = 1; i < 31; i++) {
      const hi = Math.trunc(word / 127773);
      word = 16807 * (word - hi * 127773) - 2836 * hi;
      if (word < 0) word += 2147483647;
      state[i] = word;
    }
    let front = 3;
    let rear = 0;
    for (let i = 0; ; i++) {
      state[front] = state[front]! + state[rear]!;
      const value = state[front]! >>> 1;
      front = (front + 1) % 31;
      rear = (rear + 1) % 31;
      if (i >= 310) yield [BigInt(value), 4];
    }
  } else if (name === 'cpp-mt19937') {
    const state = new Uint32Array(624);
    state[0] = Number(seed & 0xffffffffn);
    for (let i = 1; i < 624; i++) {
      const previous = state[i - 1]!;
      state[i] = Math.imul(1812433253, previous ^ (previous >>> 30)) + i;
    }
    while (true) {
      for (let i = 0; i < 624; i++) {
        const joined = (state[i]! & 0x80000000) | (state[(i + 1) % 624]! & 0x7fffffff);
        state[i] = state[(i + 397) % 624]! ^ (joined >>> 1) ^ (joined & 1 ? 0x9908b0df : 0);
      }
      for (let value of state) {
        value ^= value >>> 11;
        value ^= (value << 7) & 0x9d2c5680;
        value ^= (value << 15) & 0xefc60000;
        value ^= value >>> 18;
        yield [BigInt(value >>> 0), 4];
      }
    }
  } else {
    let state: bigint;
    let multiplier: bigint;
    let increment: bigint;
    let mask: bigint;
    if (name === 'lrand48' || name === 'drand48-state' || name === 'java-random') {
      state =
        name === 'java-random' ? seed ^ 0x5deece66dn : ((seed & 0xffffffffn) << 16n) | 0x330en;
      multiplier = 0x5deece66dn;
      increment = 11n;
      mask = (1n << 48n) - 1n;
    } else if (name === 'ansi-lcg31') {
      state = seed;
      multiplier = 1103515245n;
      increment = 12345n;
      mask = 0x7fffffffn;
    } else if (name === 'msvc-rand') {
      state = seed;
      multiplier = 214013n;
      increment = 2531011n;
      mask = 0xffffffffn;
    } else throw new Error(`Unknown RNG: ${name}`);
    state &= mask;
    while (true) {
      state = (multiplier * state + increment) & mask;
      const value =
        name === 'lrand48'
          ? state >> 17n
          : name === 'java-random'
            ? state >> 16n
            : name === 'msvc-rand'
              ? (state >> 16n) & 0x7fffn
              : state;
      yield [value, WIDTHS[name]];
    }
  }
}
