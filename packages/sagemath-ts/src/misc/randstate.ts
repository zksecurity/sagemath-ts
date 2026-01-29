/**
 * @module sage/misc/randstate
 * @description Random state utilities (partial port)
 *
 * Port of: sage/misc/randstate.pyx (interface-level parity)
 * Reference: reference/sage/src/sage/misc/randstate.pyx
 */

import { ValueError } from '../errors.js';

// Sage uses a 31-bit RAND_MAX for some distributions.
export const SAGE_RAND_MAX = 0x7fffffff;

const UINT64_MASK = (1n << 64n) - 1n;

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
 * Random state container.
 *
 * This is not a cryptographic RNG. It is intended to match Sage's
 * centralized randstate semantics for reproducibility.
 */
export class RandState {
  private state: bigint;

  constructor(seed?: bigint | number) {
    const seedValue = seed === undefined ? defaultSeed() : BigInt(seed);
    this.state = seedValue & UINT64_MASK;
  }

  /**
   * Set the RNG seed.
   */
  set_seed(seed: bigint | number): void {
    this.state = BigInt(seed) & UINT64_MASK;
  }

  /**
   * Return the current seed value.
   */
  seed(): bigint {
    return this.state;
  }

  /**
   * Generate next 64-bit value using an LCG.
   */
  private next_uint64(): bigint {
    // Constants from MMIX by Donald Knuth.
    const a = 6364136223846793005n;
    const c = 1442695040888963407n;
    this.state = (a * this.state + c) & UINT64_MASK;
    return this.state;
  }

  /**
   * Return a 31-bit random integer (0 <= x <= SAGE_RAND_MAX).
   */
  c_random(): number {
    return Number(this.next_uint64() & BigInt(SAGE_RAND_MAX));
  }

  /**
   * Return a floating-point number in [0, 1).
   */
  random(): number {
    const top53 = this.next_uint64() >> 11n; // 64 - 11 = 53 bits
    return Number(top53) / 2 ** 53;
  }

  /**
   * Return a random bigint with exactly `bits` bits (unless bits is 0).
   */
  random_bits(bits: number): bigint {
    if (bits < 0) {
      throw new ValueError('number of bits must be nonnegative');
    }
    if (bits === 0) {
      return 0n;
    }

    let remaining = bits;
    let result = 0n;
    while (remaining > 0) {
      const chunkBits = Math.min(remaining, 64);
      const mask = (1n << BigInt(chunkBits)) - 1n;
      const chunk = this.next_uint64() & mask;
      result = (result << BigInt(chunkBits)) | chunk;
      remaining -= chunkBits;
    }
    return result;
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
