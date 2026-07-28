/**
 * @module sage/misc/randstate
 * @description Random state utilities (partial port)
 *
 * Port of: sage/misc/randstate.pyx (interface-level parity)
 * Reference: reference/sage/src/sage/misc/randstate.pyx
 *
 * Sage's `randstate` wraps a GMP `gmp_randstate_t` initialised with
 * `gmp_randinit_default` (`randstate.pyx:511`), which is GMP's Mersenne
 * Twister (`gmp_randinit_mt`).  All random bits are drawn through
 * `gmp_urandomb_ui(state, nbits)` (`randstate.pyx:881`).
 *
 * We therefore implement the reference MT19937 generator of Matsumoto and
 * Nishimura and expose the same `gmp_urandomb_ui` style primitive: a request
 * for `nbits` bits consumes `ceil(nbits/32)` 32-bit outputs, assembles them
 * little-endian (as GMP's `randget_mt` fills limbs) and masks off the excess
 * high bits.
 *
 * @see Deviation: the seeding step differs from GMP's `randseed_mt` (which is
 * not part of the vendored reference sources), so the concrete stream for a
 * given seed differs from Sage's even though the generator family, the
 * per-call bit consumption and the output distributions match.
 */

import { ValueError } from '../errors.js';

// Sage uses a 31-bit RAND_MAX for some distributions.
export const SAGE_RAND_MAX = 0x7fffffff;

/** Number of 32-bit words in the MT19937 state. */
const MT_N = 624;
/** MT19937 shift parameter. */
const MT_M = 397;
/** MT19937 twist constant. */
const MT_MATRIX_A = 0x9908b0df;
/** Most significant bit. */
const MT_UPPER_MASK = 0x80000000;
/** Least significant 31 bits. */
const MT_LOWER_MASK = 0x7fffffff;

let _currentRandState: RandState | null = null;

function defaultSeed(): bigint {
  const timeSeed = BigInt(Date.now());
  const perfSeed =
    typeof performance !== 'undefined' ? BigInt(Math.floor(performance.now() * 1_000_000)) : 0n;
  const hrSeed =
    typeof process !== 'undefined' && typeof process.hrtime?.bigint === 'function'
      ? process.hrtime.bigint()
      : 0n;

  if (typeof globalThis.crypto?.getRandomValues === 'function') {
    const buf = new Uint32Array(2);
    globalThis.crypto.getRandomValues(buf);
    const cryptoSeed = (BigInt(buf[0]) << 32n) | BigInt(buf[1]);
    return cryptoSeed ^ timeSeed ^ perfSeed ^ hrSeed;
  }

  return timeSeed ^ perfSeed ^ hrSeed;
}

/**
 * Split a nonnegative bigint into 32-bit words, least significant first.
 *
 * Mirrors GMP's `mpz_export(..., order = -1, size = 4, ...)`, which is how the
 * seed reaches the Mersenne Twister state.
 */
function seedWords(seed: bigint): number[] {
  // GMP's `mpz_export` writes the absolute value, which is why Sage's
  // `set_random_seed(-12345)` and `set_random_seed(12345)` agree
  // (`randstate.pyx:955-965`).
  let s = seed < 0n ? -seed : seed;
  const words: number[] = [];
  while (s > 0n) {
    words.push(Number(s & 0xffffffffn));
    s >>= 32n;
  }
  if (words.length === 0) {
    words.push(0);
  }
  return words;
}

/**
 * Random state container.
 *
 * This is not a cryptographic RNG. It is intended to match Sage's
 * centralized randstate semantics for reproducibility.
 */
export class RandState {
  /** The MT19937 state vector. */
  private readonly mt: Uint32Array;
  /** Index into the state vector; `MT_N` means "state exhausted". */
  private mti: number;
  /**
   * The seed this state was created (or reseeded) with.
   *
   * Sage keeps this separately from the generator state so that `seed()`
   * does not change as random numbers are drawn (`randstate.pyx:562-580`).
   */
  private _seed: bigint;

  constructor(seed?: bigint | number) {
    this.mt = new Uint32Array(MT_N);
    this.mti = MT_N + 1;
    this._seed = 0n;
    this.set_seed(seed === undefined ? defaultSeed() : BigInt(seed));
  }

  /**
   * Set the RNG seed.
   */
  set_seed(seed: bigint | number): void {
    const seedValue = BigInt(seed);
    this._seed = seedValue;
    this._init_by_array(seedWords(seedValue));
  }

  /**
   * Return the initial seed of this random state.
   *
   * This is *not* the current state: it does not change when random numbers
   * are drawn.
   *
   * Reference: `randstate.pyx:562` `randstate.seed`.
   */
  seed(): bigint {
    return this._seed;
  }

  /**
   * MT19937 `init_genrand`: seed the state from a single 32-bit word.
   */
  private _init_genrand(s: number): void {
    this.mt[0] = s >>> 0;
    for (let i = 1; i < MT_N; i++) {
      const prev = this.mt[i - 1]!;
      this.mt[i] = (Math.imul(1812433253, prev ^ (prev >>> 30)) + i) >>> 0;
    }
    this.mti = MT_N;
  }

  /**
   * MT19937 `init_by_array`: seed the state from an arbitrary-length key.
   *
   * This lets arbitrarily large integer seeds be used, as GMP does.
   */
  private _init_by_array(key: number[]): void {
    this._init_genrand(19650218);
    let i = 1;
    let j = 0;
    let k = Math.max(MT_N, key.length);
    for (; k > 0; k--) {
      const prev = this.mt[i - 1]!;
      this.mt[i] =
        (((this.mt[i]! ^ Math.imul(prev ^ (prev >>> 30), 1664525)) >>> 0) + key[j]! + j) >>> 0;
      i++;
      j++;
      if (i >= MT_N) {
        this.mt[0] = this.mt[MT_N - 1]!;
        i = 1;
      }
      if (j >= key.length) {
        j = 0;
      }
    }
    for (k = MT_N - 1; k > 0; k--) {
      const prev = this.mt[i - 1]!;
      this.mt[i] = (((this.mt[i]! ^ Math.imul(prev ^ (prev >>> 30), 1566083941)) >>> 0) - i) >>> 0;
      i++;
      if (i >= MT_N) {
        this.mt[0] = this.mt[MT_N - 1]!;
        i = 1;
      }
    }
    // MSB is 1, assuring a non-zero initial array.
    this.mt[0] = MT_UPPER_MASK;
    this.mti = MT_N;
  }

  /**
   * MT19937 `genrand_int32`: generate the next 32-bit output word.
   */
  private next_uint32(): number {
    if (this.mti >= MT_N) {
      // Generate MT_N words at a time.
      for (let kk = 0; kk < MT_N - MT_M; kk++) {
        const y = (this.mt[kk]! & MT_UPPER_MASK) | (this.mt[kk + 1]! & MT_LOWER_MASK);
        this.mt[kk] = (this.mt[kk + MT_M]! ^ (y >>> 1) ^ (y & 1 ? MT_MATRIX_A : 0)) >>> 0;
      }
      for (let kk = MT_N - MT_M; kk < MT_N - 1; kk++) {
        const y = (this.mt[kk]! & MT_UPPER_MASK) | (this.mt[kk + 1]! & MT_LOWER_MASK);
        this.mt[kk] = (this.mt[kk + (MT_M - MT_N)]! ^ (y >>> 1) ^ (y & 1 ? MT_MATRIX_A : 0)) >>> 0;
      }
      const y = (this.mt[MT_N - 1]! & MT_UPPER_MASK) | (this.mt[0]! & MT_LOWER_MASK);
      this.mt[MT_N - 1] = (this.mt[MT_M - 1]! ^ (y >>> 1) ^ (y & 1 ? MT_MATRIX_A : 0)) >>> 0;
      this.mti = 0;
    }

    let y = this.mt[this.mti++]!;

    // Tempering.
    y ^= y >>> 11;
    y = (y ^ ((y << 7) & 0x9d2c5680)) >>> 0;
    y = (y ^ ((y << 15) & 0xefc60000)) >>> 0;
    y ^= y >>> 18;

    return y >>> 0;
  }

  /**
   * Return a uniformly random integer in `[0, 2^nbits)`.
   *
   * This is GMP's `gmp_urandomb_ui`/`randget_mt`: fill the requested number of
   * bits from consecutive 32-bit generator outputs (least significant word
   * first) and mask off the excess high bits.
   */
  private urandomb(nbits: number): bigint {
    if (nbits === 0) {
      return 0n;
    }
    let result = 0n;
    let shift = 0n;
    for (let remaining = nbits; remaining > 0; remaining -= 32) {
      result |= BigInt(this.next_uint32()) << shift;
      shift += 32n;
    }
    return result & ((1n << BigInt(nbits)) - 1n);
  }

  /**
   * Return a 31-bit random number (0 <= x <= SAGE_RAND_MAX).
   *
   * Reference: `randstate.pyx:881` `randstate.c_random`.
   */
  c_random(): number {
    return Number(this.urandomb(31));
  }

  /**
   * Return a random floating-point number in [0, 1).
   *
   * Ported verbatim from `randstate.pyx:883-895`: it consumes *two* draws,
   * 25 bits scaled by 2^-25 plus 28 bits scaled by 2^-53.
   */
  c_rand_double(): number {
    const a = Number(this.urandomb(25)) * (1.0 / 33554432.0); // divide by 2^25
    const b = Number(this.urandomb(28)) * (1.0 / 9007199254740992.0); // divide by 2^53
    return a + b;
  }

  /**
   * Return a floating-point number in [0, 1).
   *
   * @deprecated Sage spells this `c_rand_double`; `random()` in
   * `sage.misc.randstate` is the module-level 31-bit integer generator. Kept
   * as an alias for existing callers.
   * @see Deviation: extra method not present in Sage.
   */
  random(): number {
    return this.c_rand_double();
  }

  /**
   * Return a uniformly random bigint in `[0, 2^bits)`.
   *
   * Mirrors GMP's `mpz_urandomb` / `gmp_urandomb_ui`.
   */
  random_bits(bits: number): bigint {
    if (bits < 0) {
      throw new ValueError('number of bits must be nonnegative');
    }
    return this.urandomb(bits);
  }

  /**
   * Return a random bigint in [0, n).
   */
  random_below(n: bigint): bigint {
    if (n <= 0n) {
      throw new ValueError('n must be positive');
    }
    const bits = n.toString(2).length;
    while (true) {
      const candidate = this.random_bits(bits);
      if (candidate < n) {
        return candidate;
      }
    }
  }

  /**
   * Return a random bigint in [min, max] inclusive.
   */
  randint(min: bigint, max: bigint): bigint {
    if (min > max) {
      throw new ValueError('min must be <= max');
    }
    const range = max - min + 1n;
    return min + this.random_below(range);
  }
}

/**
 * Return the current global random state.
 */
export function current_randstate(): RandState {
  if (_currentRandState === null) {
    _currentRandState = new RandState();
  }
  return _currentRandState;
}

/**
 * Set the global random seed.
 */
export function set_random_seed(seed?: number): void {
  if (_currentRandState === null) {
    _currentRandState = new RandState(seed);
  } else if (seed === undefined) {
    _currentRandState.set_seed(defaultSeed());
  } else {
    _currentRandState.set_seed(seed);
  }
}

/**
 * Return a 31-bit random number.
 *
 * Intended as a drop-in replacement for the libc `random()` function.
 *
 * Reference: `randstate.pyx:1008` `random`.
 */
export function random(): bigint {
  return BigInt(current_randstate().c_random());
}

/**
 * Return the initial seed used to create the current random state.
 *
 * Reference: `randstate.pyx:1024` `initial_seed`.
 */
export function initial_seed(): bigint {
  return current_randstate().seed();
}
