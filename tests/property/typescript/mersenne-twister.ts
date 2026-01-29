/**
 * Mersenne Twister (MT19937) implementation matching Python's random module.
 *
 * This implementation produces identical sequences to Python's random.Random
 * for the same seed, allowing property tests to generate matching test cases.
 *
 * Reference: https://en.wikipedia.org/wiki/Mersenne_Twister
 * Python source: https://github.com/python/cpython/blob/main/Modules/_randommodule.c
 */

// MT19937 constants
const N = 624;
const M = 397;
const MATRIX_A = 0x9908b0df;
const UPPER_MASK = 0x80000000;
const LOWER_MASK = 0x7fffffff;

/**
 * Mersenne Twister PRNG matching Python's random module.
 */
export class MersenneTwister {
  private mt: Uint32Array;
  private mti: number;

  constructor(seed?: number | bigint) {
    this.mt = new Uint32Array(N);
    this.mti = N + 1;
    if (seed !== undefined) {
      this.seed(seed);
    }
  }

  /**
   * Initialize the generator with init_genrand (basic MT seeding).
   */
  private initGenrand(seed: number): void {
    seed = seed >>> 0;
    this.mt[0] = seed;
    for (let i = 1; i < N; i++) {
      const s = this.mt[i - 1]! ^ (this.mt[i - 1]! >>> 30);
      // mt[i] = 1812433253 * (mt[i-1] ^ (mt[i-1] >> 30)) + i
      // Use BigInt for intermediate to avoid JS number overflow
      const product = BigInt(s) * 1812433253n;
      this.mt[i] = Number((product + BigInt(i)) & 0xffffffffn);
    }
    this.mti = N;
  }

  /**
   * Initialize with an array of seeds (matching Python's init_by_array).
   */
  private initByArray(initKey: number[]): void {
    this.initGenrand(19650218);

    let i = 1;
    let j = 0;
    let k = N > initKey.length ? N : initKey.length;

    for (; k > 0; k--) {
      const s = this.mt[i - 1]! ^ (this.mt[i - 1]! >>> 30);
      // mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525)) + init_key[j] + j
      const product = BigInt(s) * 1664525n;
      this.mt[i] = Number(
        (BigInt(this.mt[i]!) ^ product) & 0xffffffffn
      );
      this.mt[i] = ((this.mt[i]! + initKey[j]!) >>> 0);
      this.mt[i] = ((this.mt[i]! + j) >>> 0);
      i++;
      j++;
      if (i >= N) {
        this.mt[0] = this.mt[N - 1]!;
        i = 1;
      }
      if (j >= initKey.length) {
        j = 0;
      }
    }

    for (k = N - 1; k > 0; k--) {
      const s = this.mt[i - 1]! ^ (this.mt[i - 1]! >>> 30);
      // mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941)) - i
      const product = BigInt(s) * 1566083941n;
      this.mt[i] = Number(
        (BigInt(this.mt[i]!) ^ product) & 0xffffffffn
      );
      this.mt[i] = (this.mt[i]! - i) >>> 0;
      i++;
      if (i >= N) {
        this.mt[0] = this.mt[N - 1]!;
        i = 1;
      }
    }

    // MSB is 1; assuring non-zero initial array
    this.mt[0] = 0x80000000;
  }

  /**
   * Initialize the generator with a seed (matching Python's seed()).
   * Python converts the integer to 32-bit words and uses init_by_array.
   */
  seed(seed: number | bigint): void {
    // Convert to bigint for consistent handling
    let n = typeof seed === 'number' ? BigInt(seed) : seed;

    // Handle negative seeds like Python does (use absolute value)
    if (n < 0n) {
      n = -n;
    }

    // Convert to array of 32-bit words (little-endian, like Python)
    const initKey: number[] = [];
    if (n === 0n) {
      initKey.push(0);
    } else {
      while (n > 0n) {
        initKey.push(Number(n & 0xffffffffn));
        n >>= 32n;
      }
    }

    this.initByArray(initKey);
  }

  /**
   * Generate a random 32-bit unsigned integer.
   */
  private genrandUint32(): number {
    let y: number;
    const mag01 = [0, MATRIX_A];

    if (this.mti >= N) {
      let kk: number;

      // Generate N words at once
      for (kk = 0; kk < N - M; kk++) {
        y = (this.mt[kk]! & UPPER_MASK) | (this.mt[kk + 1]! & LOWER_MASK);
        this.mt[kk] = this.mt[kk + M]! ^ (y >>> 1) ^ mag01[y & 1]!;
      }
      for (; kk < N - 1; kk++) {
        y = (this.mt[kk]! & UPPER_MASK) | (this.mt[kk + 1]! & LOWER_MASK);
        this.mt[kk] = this.mt[kk + (M - N)]! ^ (y >>> 1) ^ mag01[y & 1]!;
      }
      y = (this.mt[N - 1]! & UPPER_MASK) | (this.mt[0]! & LOWER_MASK);
      this.mt[N - 1] = this.mt[M - 1]! ^ (y >>> 1) ^ mag01[y & 1]!;

      this.mti = 0;
    }

    y = this.mt[this.mti++]!;

    // Tempering
    y ^= y >>> 11;
    y ^= (y << 7) & 0x9d2c5680;
    y ^= (y << 15) & 0xefc60000;
    y ^= y >>> 18;

    return y >>> 0;
  }

  /**
   * Generate k random bits as a bigint (matching Python's getrandbits).
   * Python takes high bits from each 32-bit word, not low bits.
   */
  getrandbits(k: number): bigint {
    if (k === 0) return 0n;

    let result = 0n;
    let bitsNeeded = k;

    while (bitsNeeded > 0) {
      const word = this.genrandUint32();
      const bitsFromWord = Math.min(32, bitsNeeded);

      // Python takes high bits: shift right to get top bitsFromWord bits
      const bits = word >>> (32 - bitsFromWord);

      result = (result << BigInt(bitsFromWord)) | BigInt(bits);
      bitsNeeded -= bitsFromWord;
    }

    return result;
  }

  /**
   * Generate a random integer in [0, n) (matching Python's _randbelow_with_getrandbits).
   */
  randbelow(n: bigint): bigint {
    if (n <= 0n) {
      throw new Error('n must be positive');
    }
    if (n === 1n) return 0n;

    // Calculate number of bits needed
    const k = bitLength(n);

    // Rejection sampling to ensure uniform distribution
    let r = this.getrandbits(k);
    while (r >= n) {
      r = this.getrandbits(k);
    }
    return r;
  }

  /**
   * Generate a random integer in [a, b] inclusive (matching Python's randint).
   */
  randint(a: bigint, b: bigint): bigint {
    if (a > b) {
      throw new Error(`empty range for randint(${a}, ${b})`);
    }
    return a + this.randbelow(b - a + 1n);
  }

  /**
   * Generate a random floating point number in [0, 1) (matching Python's random()).
   */
  random(): number {
    const a = this.genrandUint32() >>> 5; // 27 bits
    const b = this.genrandUint32() >>> 6; // 26 bits
    return (a * 67108864.0 + b) / 9007199254740992.0;
  }
}

/**
 * Calculate the bit length of a bigint (like Python's int.bit_length()).
 */
function bitLength(n: bigint): number {
  if (n === 0n) return 0;
  if (n < 0n) n = -n;

  let bits = 0;
  while (n > 0n) {
    bits++;
    n >>= 1n;
  }
  return bits;
}

/**
 * Default instance for convenience.
 */
export const mt = new MersenneTwister();
