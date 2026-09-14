/**
 * NTL arbitrary precision integer (ZZ).
 * @see Reference: ntl/src/ZZ.cpp
 *
 * ZZ represents arbitrary precision integers with support for
 * arithmetic operations, GCD computations, primality testing,
 * and modular arithmetic.
 */

/**
 * Arbitrary precision integer class.
 * TypeScript port of NTL's ZZ class.
 */
export class ZZ {
  private _value: bigint;

  /**
   * Creates a new ZZ instance.
   * @param value - Initial value (default 0)
   */
  constructor(value?: bigint | number | string) {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ');
  }

  /**
   * Returns the zero element.
   */
  static zero(): ZZ {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ.zero');
  }

  /**
   * Returns the one element.
   */
  static one(): ZZ {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ.one');
  }

  // ============================================
  // Arithmetic Operations
  // ============================================

  /**
   * Adds two ZZ values.
   * @param other - The value to add
   * @returns The sum
   */
  add(other: ZZ): ZZ {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ.add');
  }

  /**
   * Subtracts a ZZ value.
   * @param other - The value to subtract
   * @returns The difference
   */
  sub(other: ZZ): ZZ {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ.sub');
  }

  /**
   * Multiplies two ZZ values.
   * @param other - The value to multiply
   * @returns The product
   */
  mul(other: ZZ): ZZ {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ.mul');
  }

  /**
   * Divides by a ZZ value (integer division).
   * @param other - The divisor
   * @returns The quotient
   */
  div(other: ZZ): ZZ {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ.div');
  }

  /**
   * Computes remainder after division.
   * @param other - The divisor
   * @returns The remainder
   */
  rem(other: ZZ): ZZ {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ.rem');
  }

  /**
   * Computes quotient and remainder.
   * @param other - The divisor
   * @returns Tuple of [quotient, remainder]
   */
  divRem(other: ZZ): [ZZ, ZZ] {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ.divRem');
  }

  /**
   * Negates the value.
   * @returns The negation
   */
  negate(): ZZ {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ.negate');
  }

  /**
   * Computes absolute value.
   * @returns The absolute value
   */
  abs(): ZZ {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ.abs');
  }

  /**
   * Squares the value.
   * @returns The square
   */
  sqr(): ZZ {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ.sqr');
  }

  // ============================================
  // GCD and Related Functions
  // ============================================

  /**
   * Computes GCD of two ZZ values.
   * @param a - First value
   * @param b - Second value
   * @returns The GCD
   */
  static GCD(a: ZZ, b: ZZ): ZZ {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ.GCD');
  }

  /**
   * Computes LCM of two ZZ values.
   * @param a - First value
   * @param b - Second value
   * @returns The LCM
   */
  static LCM(a: ZZ, b: ZZ): ZZ {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ.LCM');
  }

  /**
   * Extended GCD computation.
   * Computes d, s, t such that d = GCD(a, b) = s*a + t*b
   * @param a - First value
   * @param b - Second value
   * @returns Tuple [d, s, t] where d = s*a + t*b
   */
  static XGCD(a: ZZ, b: ZZ): [ZZ, ZZ, ZZ] {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ.XGCD');
  }

  // ============================================
  // Modular Arithmetic
  // ============================================

  /**
   * Computes (a + b) mod n.
   * @param a - First operand
   * @param b - Second operand
   * @param n - Modulus
   * @returns (a + b) mod n
   */
  static AddMod(a: ZZ, b: ZZ, n: ZZ): ZZ {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ.AddMod');
  }

  /**
   * Computes (a - b) mod n.
   * @param a - First operand
   * @param b - Second operand
   * @param n - Modulus
   * @returns (a - b) mod n
   */
  static SubMod(a: ZZ, b: ZZ, n: ZZ): ZZ {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ.SubMod');
  }

  /**
   * Computes (a * b) mod n.
   * @param a - First operand
   * @param b - Second operand
   * @param n - Modulus
   * @returns (a * b) mod n
   */
  static MulMod(a: ZZ, b: ZZ, n: ZZ): ZZ {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ.MulMod');
  }

  /**
   * Computes a^2 mod n.
   * @param a - Base
   * @param n - Modulus
   * @returns a^2 mod n
   */
  static SqrMod(a: ZZ, n: ZZ): ZZ {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ.SqrMod');
  }

  /**
   * Computes modular inverse of a mod n.
   * @param a - Value to invert
   * @param n - Modulus
   * @returns a^(-1) mod n
   * @throws If inverse does not exist
   */
  static InvMod(a: ZZ, n: ZZ): ZZ {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ.InvMod');
  }

  /**
   * Computes modular inverse, returning status.
   * @param a - Value to invert
   * @param n - Modulus
   * @returns [status, result] where status is 0 if inverse exists
   */
  static InvModStatus(a: ZZ, n: ZZ): [number, ZZ] {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ.InvModStatus');
  }

  /**
   * Computes a^e mod n (modular exponentiation).
   * @param a - Base
   * @param e - Exponent
   * @param n - Modulus
   * @returns a^e mod n
   */
  static PowerMod(a: ZZ, e: ZZ, n: ZZ): ZZ {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ.PowerMod');
  }

  // ============================================
  // Primality Testing
  // ============================================

  /**
   * Miller-Rabin primality test.
   * @param n - Number to test
   * @param NumTrials - Number of trials (default 10)
   * @returns 1 if probably prime, 0 otherwise
   */
  static ProbPrime(n: ZZ, NumTrials?: number): number {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ.ProbPrime');
  }

  /**
   * Miller witness test.
   * @param n - Number to test
   * @param x - Witness
   * @returns 1 if x is a witness that n is composite
   */
  static MillerWitness(n: ZZ, x: ZZ): number {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ.MillerWitness');
  }

  /**
   * Generates a random prime of given bit length.
   * @param l - Bit length
   * @param NumTrials - Number of primality trials
   * @returns A random prime
   */
  static RandomPrime(l: number, NumTrials?: number): ZZ {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ.RandomPrime');
  }

  /**
   * Finds next prime >= m.
   * @param m - Starting point
   * @param NumTrials - Number of primality trials
   * @returns The next prime
   */
  static NextPrime(m: ZZ, NumTrials?: number): ZZ {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ.NextPrime');
  }

  // ============================================
  // Bit Operations
  // ============================================

  /**
   * Returns the number of bits in the representation.
   * @returns Number of bits
   */
  NumBits(): number {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ.NumBits');
  }

  /**
   * Returns bit at position k.
   * @param k - Bit position (0-indexed)
   * @returns 0 or 1
   */
  bit(k: number): number {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ.bit');
  }

  /**
   * Left shift by n bits.
   * @param n - Number of bits to shift
   * @returns Shifted value
   */
  LeftShift(n: number): ZZ {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ.LeftShift');
  }

  /**
   * Right shift by n bits.
   * @param n - Number of bits to shift
   * @returns Shifted value
   */
  RightShift(n: number): ZZ {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ.RightShift');
  }

  /**
   * Makes the number odd by dividing out powers of 2.
   * @returns The number of 2s divided out
   */
  MakeOdd(): number {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ.MakeOdd');
  }

  // ============================================
  // Comparison and Utility
  // ============================================

  /**
   * Returns the sign: -1, 0, or 1.
   * @returns Sign of the value
   */
  sign(): number {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ.sign');
  }

  /**
   * Checks if value is zero.
   * @returns True if zero
   */
  IsZero(): boolean {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ.IsZero');
  }

  /**
   * Checks if value is one.
   * @returns True if one
   */
  IsOne(): boolean {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ.IsOne');
  }

  /**
   * Checks if value is odd.
   * @returns True if odd
   */
  IsOdd(): boolean {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ.IsOdd');
  }

  /**
   * Compares with another ZZ.
   * @param other - Value to compare
   * @returns -1, 0, or 1
   */
  compare(other: ZZ): number {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ.compare');
  }

  /**
   * Checks equality.
   * @param other - Value to compare
   * @returns True if equal
   */
  equals(other: ZZ): boolean {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ.equals');
  }

  /**
   * Divisibility test: checks whether b divides a
   * (NTL: `long divide(const ZZ& a, const ZZ& b)`, "if b | a, returns 1;
   * otherwise returns 0", ntl/doc/ZZ.txt:234).
   * @param a - Dividend
   * @param b - Potential divisor
   * @returns True if b divides a
   */
  static divide(a: ZZ, b: ZZ): boolean {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ.divide');
  }

  // ============================================
  // Conversion
  // ============================================

  /**
   * Converts to bigint.
   * @returns The value as bigint
   */
  toBigInt(): bigint {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ.toBigInt');
  }

  /**
   * Converts to number (may lose precision).
   * @returns The value as number
   */
  toNumber(): number {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ.toNumber');
  }

  /**
   * Converts to string representation.
   * @returns String representation
   */
  toString(): string {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ.toString');
  }

  // ============================================
  // Random Number Generation
  // ============================================

  /**
   * Generates random integer in [0, n).
   * @param n - Upper bound (exclusive)
   * @returns Random integer
   */
  static RandomBnd(n: ZZ): ZZ {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ.RandomBnd');
  }

  /**
   * Generates random integer with exactly l bits.
   * @param l - Bit length
   * @returns Random integer
   */
  static RandomLen(l: number): ZZ {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ.RandomLen');
  }

  /**
   * Generates random integer with at most l bits.
   * @param l - Maximum bit length
   * @returns Random integer
   */
  static RandomBits(l: number): ZZ {
    throw new Error('NTL_NOT_IMPLEMENTED: ZZ.RandomBits');
  }
}

// ============================================
// Standalone Functions (NTL-style API)
// ============================================

/**
 * Computes GCD of two values.
 */
export function GCD(a: ZZ, b: ZZ): ZZ {
  return ZZ.GCD(a, b);
}

/**
 * Computes LCM of two values.
 */
export function LCM(a: ZZ, b: ZZ): ZZ {
  return ZZ.LCM(a, b);
}

/**
 * Extended GCD.
 */
export function XGCD(a: ZZ, b: ZZ): [ZZ, ZZ, ZZ] {
  return ZZ.XGCD(a, b);
}

/**
 * Modular exponentiation.
 */
export function PowerMod(a: ZZ, e: ZZ, n: ZZ): ZZ {
  return ZZ.PowerMod(a, e, n);
}

/**
 * Modular inverse.
 */
export function InvMod(a: ZZ, n: ZZ): ZZ {
  return ZZ.InvMod(a, n);
}

/**
 * Find next prime.
 */
export function NextPrime(m: ZZ, NumTrials?: number): ZZ {
  return ZZ.NextPrime(m, NumTrials);
}

// Native 60-bit scalar profile: NTL_PRIME_BND=(1<<14)-1.
const PRIME_BND = 16383;
const PRIME_LIMIT = BigInt((2 * PRIME_BND + 1) ** 2);
let lowsieve_storage: Uint8Array | undefined;
/** NTL segmented prime iterator, ending at the native small-prime bound.
 * @see Deviation: NTL prime sequence adapters
 */
export class PrimeSeq {
  private movesieve: Uint8Array | undefined;
  private movesieve_mem: Uint8Array | undefined;
  private pindex = -1;
  private pshift = -1;
  private exhausted = false;
  /** Return the next prime, or zero after the native sequence is exhausted. */
  next(): bigint {
    if (this.exhausted) return 0n;
    if (this.pshift < 0) {
      this.shift(0);
      return 2n;
    }
    for (;;) {
      let i = this.pindex;
      while (++i < PRIME_BND) {
        if (this.movesieve![i]) {
          this.pindex = i;
          return BigInt(this.pshift + 2 * i + 3);
        }
      }
      const newshift = this.pshift + 2 * PRIME_BND;
      if (newshift > 2 * PRIME_BND * (2 * PRIME_BND + 1)) {
        this.exhausted = true;
        return 0n;
      }
      this.shift(newshift);
    }
  }
  /** Reset so the next result is the first prime greater than or equal to b. */
  reset(b: bigint): void {
    if (b > PRIME_LIMIT) {
      this.exhausted = true;
      return;
    }
    if (b <= 2n) {
      this.shift(-1);
      return;
    }
    if (b % 2n === 0n) b++;
    const bound = Number(b);
    this.shift(Math.floor((bound - 3) / (2 * PRIME_BND)) * (2 * PRIME_BND));
    this.pindex = (bound - this.pshift - 3) / 2 - 1;
  }
  private shift(newshift: number): void {
    if (!lowsieve_storage) this.start();
    const low = lowsieve_storage!;
    if (newshift < 0) this.pshift = -1;
    else if (newshift === 0) {
      this.pshift = 0;
      this.movesieve = low;
    } else if (newshift !== this.pshift) {
      this.movesieve_mem ??= new Uint8Array(PRIME_BND);
      this.pshift = newshift;
      const p = this.movesieve_mem;
      this.movesieve = p;
      p.fill(1);
      const bound = this.pshift + 2 * PRIME_BND + 1;
      for (let i = 0, step = 3; step * step <= bound; i++, step += 2) {
        if (!low[i]) continue;
        let start = Math.floor((this.pshift + 2) / step) + 1;
        if (start % 2 === 0) start++;
        if (start <= step) start = step;
        start = (start * step - this.pshift - 3) / 2;
        for (let j = start; j < PRIME_BND; j += step) p[j] = 0;
      }
    }
    this.pindex = -1;
    this.exhausted = false;
  }
  private start(): void {
    if (lowsieve_storage) return;
    const p = new Uint8Array(PRIME_BND).fill(1);
    // floor(sqrt(2*PRIME_BND+1))=181 for this native word profile.
    const bound = Math.trunc((181 - 3) / 2);
    for (let i = 0, step = 1, start = -1; i <= bound; i++) {
      step += 2;
      start += 2 * (step - 1);
      if (p[i]) for (let j = start; j < PRIME_BND; j += step) p[j] = 0;
    }
    lowsieve_storage = p;
  }
}

// SHA-256 follows NTL's adaptation of Brad Conte's public-domain implementation.
const sha256_const = new Uint32Array([
  0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
  0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
  0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
  0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
  0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
  0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
  0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
  0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
]);
const rotateRight = (x: number, n: number): number => (x >>> n) | (x << (32 - n));
function sha256Parts(parts: readonly Uint8Array[]): Uint8Array {
  const state = new Uint32Array([
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19,
  ]);
  const words = new Uint32Array(64),
    buffer = new Uint8Array(64);
  let used = 0,
    length = 0n;
  const transform = (data: Uint8Array, offset: number): void => {
    for (let i = 0; i < 16; i++) {
      const j = offset + 4 * i;
      words[i] = (data[j]! << 24) | (data[j + 1]! << 16) | (data[j + 2]! << 8) | data[j + 3]!;
    }
    for (let i = 16; i < 64; i++) {
      const x = words[i - 15]!,
        y = words[i - 2]!;
      words[i] =
        (rotateRight(y, 17) ^ rotateRight(y, 19) ^ (y >>> 10)) +
        words[i - 7]! +
        (rotateRight(x, 7) ^ rotateRight(x, 18) ^ (x >>> 3)) +
        words[i - 16]!;
    }
    let [a, b, c, d, e, f, g, h] = Array.from(state) as [
      number,
      number,
      number,
      number,
      number,
      number,
      number,
      number,
    ];
    for (let i = 0; i < 64; i++) {
      const t1 =
        (h +
          (rotateRight(e, 6) ^ rotateRight(e, 11) ^ rotateRight(e, 25)) +
          ((e & f) ^ (~e & g)) +
          sha256_const[i]! +
          words[i]!) >>>
        0;
      const t2 =
        ((rotateRight(a, 2) ^ rotateRight(a, 13) ^ rotateRight(a, 22)) +
          ((a & b) ^ (a & c) ^ (b & c))) >>>
        0;
      [a, b, c, d, e, f, g, h] = [(t1 + t2) >>> 0, a, b, c, (d + t1) >>> 0, e, f, g];
    }
    const final = [a, b, c, d, e, f, g, h];
    for (let i = 0; i < 8; i++) state[i] = (state[i]! + final[i]!) >>> 0;
  };
  for (const part of parts) {
    length += BigInt(part.length);
    for (let offset = 0; offset < part.length; ) {
      const take = Math.min(64 - used, part.length - offset);
      buffer.set(part.subarray(offset, offset + take), used);
      used += take;
      offset += take;
      if (used === 64) {
        transform(buffer, 0);
        used = 0;
      }
    }
  }
  buffer[used++] = 0x80;
  if (used > 56) {
    buffer.fill(0, used);
    transform(buffer, 0);
    used = 0;
  }
  buffer.fill(0, used, 56);
  let bits = BigInt.asUintN(64, length * 8n);
  for (let i = 63; i >= 56; i--) {
    buffer[i] = Number(bits & 255n);
    bits >>= 8n;
  }
  transform(buffer, 0);
  const out = new Uint8Array(32);
  for (let i = 0; i < 8; i++) {
    const w = state[i]!;
    out[4 * i] = w >>> 24;
    out[4 * i + 1] = w >>> 16;
    out[4 * i + 2] = w >>> 8;
    out[4 * i + 3] = w;
  }
  return out;
}
/** NTL SHA-256 helper, returning the written digest prefix (at most 32 bytes).
 * @see Deviation: NTL deterministic byte stream adapters
 */
export function sha256(data: Uint8Array, length = 32): Uint8Array {
  return sha256Parts([data]).slice(0, Math.max(0, Math.min(32, length)));
}
/** NTL HMAC-SHA256 helper, with typed arrays replacing pointer/length pairs.
 * @see Deviation: NTL deterministic byte stream adapters
 */
export function hmac_sha256(key: Uint8Array, data: Uint8Array, length = 32): Uint8Array {
  const K = new Uint8Array(64);
  K.set(key.length > 64 ? sha256(key) : key);
  const inner = new Uint8Array(64),
    outer = new Uint8Array(64);
  for (let i = 0; i < 64; i++) {
    inner[i] = K[i]! ^ 0x36;
    outer[i] = K[i]! ^ 0x5c;
  }
  return sha256Parts([outer, sha256Parts([inner, data.subarray(0, data.length >>> 0)])]).slice(
    0,
    Math.max(0, Math.min(32, length))
  );
}
/** NTL HMAC counter-mode key derivation.
 * @see Deviation: NTL deterministic byte stream adapters
 */
export function DeriveKey(data: Uint8Array, length: number): Uint8Array {
  if (length < 0) throw new Error('DeriveKey: bad args');
  const key = hmac_sha256(new Uint8Array(), data),
    counter = new Uint8Array(8),
    out = new Uint8Array(length);
  for (let offset = 0; offset < length; offset += 32) {
    out.set(hmac_sha256(key, counter).subarray(0, Math.min(32, length - offset)), offset);
    for (let i = 0; i < 8; i++) {
      counter[i] = (counter[i]! + 1) & 255;
      if (counter[i] !== 0) break;
    }
  }
  return out;
}
/** Historical NTL name for initialization of its ChaCha20 state.
 * @see Deviation: NTL deterministic byte stream adapters
 */
export function salsa20_init(key: Uint8Array): Uint32Array {
  if (key.length < 32) throw new Error('RandomStream: key must contain at least 32 bytes');
  const state = new Uint32Array(16);
  state.set([0x61707865, 0x3320646e, 0x79622d32, 0x6b206574]);
  const view = new DataView(key.buffer, key.byteOffset, 32);
  for (let i = 0; i < 8; i++) state[4 + i] = view.getUint32(4 * i, true);
  return state;
}
/** Historical NTL name for the twenty ChaCha rounds, mutating the first 16 words.
 * @see Deviation: NTL deterministic byte stream adapters
 */
export function salsa20_core(x: Uint32Array): void {
  if (x.length < 16) throw new Error('salsa20_core: state must contain at least 16 words');
  const rotate = (a: number, n: number) => (a << n) | (a >>> (32 - n));
  const quarter = (a: number, b: number, c: number, d: number): void => {
    x[a] = (x[a]! + x[b]!) >>> 0;
    x[d] = rotate(x[d]! ^ x[a]!, 16);
    x[c] = (x[c]! + x[d]!) >>> 0;
    x[b] = rotate(x[b]! ^ x[c]!, 12);
    x[a] = (x[a]! + x[b]!) >>> 0;
    x[d] = rotate(x[d]! ^ x[a]!, 8);
    x[c] = (x[c]! + x[d]!) >>> 0;
    x[b] = rotate(x[b]! ^ x[c]!, 7);
  };
  for (let i = 0; i < 10; i++) {
    quarter(0, 4, 8, 12);
    quarter(1, 5, 9, 13);
    quarter(2, 6, 10, 14);
    quarter(3, 7, 11, 15);
    quarter(0, 5, 10, 15);
    quarter(1, 6, 11, 12);
    quarter(2, 7, 8, 13);
    quarter(3, 4, 9, 14);
  }
}
/** Historical NTL name for one ChaCha block; increments the input's 64-bit counter.
 * @see Deviation: NTL deterministic byte stream adapters
 */
export function salsa20_apply(state: Uint32Array): Uint32Array {
  if (state.length < 16) throw new Error('salsa20_apply: state must contain at least 16 words');
  const x = state.slice(0, 16);
  salsa20_core(x);
  for (let i = 0; i < 16; i++) x[i] = (x[i]! + state[i]!) >>> 0;
  state[12] = (state[12]! + 1) >>> 0;
  if (state[12] === 0) state[13] = (state[13]! + 1) >>> 0;
  return x;
}
/** NTL's generic ChaCha20 byte stream with its 64-bit counter and nonce.
 * @see Deviation: NTL deterministic byte stream adapters
 */
export class RandomStream {
  private readonly state: Uint32Array;
  private readonly buffer = new Uint8Array(64);
  private position = 64;
  constructor(key: Uint8Array | RandomStream) {
    if (key instanceof RandomStream) {
      this.state = key.state.slice();
      this.buffer.set(key.buffer);
      this.position = key.position;
      return;
    }
    this.state = salsa20_init(key);
  }
  /** Native copy assignment, preserving this object's identity. */
  assign(other: RandomStream): this {
    this.state.set(other.state);
    this.buffer.set(other.buffer);
    this.position = other.position;
    return this;
  }
  get(n: number): Uint8Array {
    if (n < 0) throw new Error('RandomStream::get: bad args');
    const out = new Uint8Array(n);
    for (let copied = 0; copied < n; ) {
      if (this.position === 64) {
        this.apply();
        this.position = 0;
      }
      const take = Math.min(n - copied, 64 - this.position);
      out.set(this.buffer.subarray(this.position, this.position + take), copied);
      copied += take;
      this.position += take;
    }
    return out;
  }
  set_nonce(nonce: bigint): void {
    const word = BigInt.asUintN(64, nonce);
    this.state[12] = 0;
    this.state[13] = 0;
    this.state[14] = Number(word & 0xffffffffn);
    this.state[15] = Number(word >> 32n);
    this.position = 64;
  }
  private apply(): void {
    const x = salsa20_apply(this.state);
    for (let i = 0; i < 16; i++) {
      const w = x[i]!;
      this.buffer[4 * i] = w;
      this.buffer[4 * i + 1] = w >>> 8;
      this.buffer[4 * i + 2] = w >>> 16;
      this.buffer[4 * i + 3] = w >>> 24;
    }
  }
}

function WordFromBytes(bytes: Uint8Array): bigint {
  let out = 0n;
  for (let i = bytes.length - 1; i >= 0; i--) out = (out << 8n) | BigInt(bytes[i]!);
  return out;
}
function integerFromRandomBytes(bytes: Uint8Array): bigint {
  return BigInt(
    '0x' +
      Array.from(bytes, (x) => x.toString(16).padStart(2, '0'))
        .reverse()
        .join('')
  );
}
/** NTL random unsigned word, using an explicit stream context.
 * @see Deviation: NTL integer sampling stream contexts
 */
export function RandomWord(stream: RandomStream): bigint {
  return WordFromBytes(stream.get(8));
}
/** NTL vector of random unsigned words, preserving successive byte consumption.
 * @see Deviation: NTL integer sampling stream contexts
 */
export function VectorRandomWord(k: number, stream: RandomStream): bigint[] {
  const out: bigint[] = [];
  for (let i = 0; i < k; i++) out.push(RandomWord(stream));
  return out;
}
/** NTL signed-word random bits.
 * @see Deviation: NTL integer sampling stream contexts
 */
export function RandomBits_long(l: number, stream: RandomStream): bigint {
  if (l <= 0) return 0n;
  if (l >= 64) throw new Error('RandomBits: length too big');
  return WordFromBytes(stream.get(Math.ceil(l / 8))) & ((1n << BigInt(l)) - 1n);
}
/** NTL unsigned-word random bits, including exactly 64 bits.
 * @see Deviation: NTL integer sampling stream contexts
 */
export function RandomBits_ulong(l: number, stream: RandomStream): bigint {
  if (l <= 0) return 0n;
  if (l > 64) throw new Error('RandomBits: length too big');
  return WordFromBytes(stream.get(Math.ceil(l / 8))) & ((1n << BigInt(l)) - 1n);
}
/** NTL exact signed-word bit length; consumes ceil((l-1)/8) bytes.
 * @see Deviation: NTL integer sampling stream contexts
 */
export function RandomLen_long(l: number, stream: RandomStream): bigint {
  if (l <= 0) return 0n;
  if (l === 1) return 1n;
  if (l >= 64) throw new Error('RandomLen: length too big');
  const mask = (1n << BigInt(l - 1)) - 1n;
  return (WordFromBytes(stream.get(Math.ceil((l - 1) / 8))) & mask) | (mask + 1n);
}
/** NTL arbitrary-integer random bits with the native resource guard.
 * @see Deviation: NTL integer sampling stream contexts
 */
export function RandomBits(l: number, stream: RandomStream): bigint {
  if (l <= 0) return 0n;
  if (l >= 2 ** 60) throw new Error('RandomBits: length too big');
  const nb = Math.ceil(l / 8),
    bytes = stream.get(nb);
  bytes[nb - 1] = bytes[nb - 1]! & ((1 << (l - (nb - 1) * 8)) - 1);
  return integerFromRandomBytes(bytes);
}
/** Native return-value alias for arbitrary-integer random bits.
 * @see Deviation: NTL integer sampling stream contexts
 */
export function RandomBits_ZZ(l: number, stream: RandomStream): bigint {
  return RandomBits(l, stream);
}
/** NTL exact arbitrary-integer bit length; consumes ceil(l/8) bytes.
 * @see Deviation: NTL integer sampling stream contexts
 */
export function RandomLen(l: number, stream: RandomStream): bigint {
  if (l <= 0) return 0n;
  if (l === 1) return 1n;
  if (l >= 2 ** 60) throw new Error('RandomLen: length too big');
  const nb = Math.ceil(l / 8),
    mask = (1 << (l - (nb - 1) * 8)) - 1,
    bytes = stream.get(nb);
  bytes[nb - 1] = (bytes[nb - 1]! & mask) | ((mask >> 1) + 1);
  return integerFromRandomBytes(bytes);
}
/** Native return-value alias for exact arbitrary-integer bit length.
 * @see Deviation: NTL integer sampling stream contexts
 */
export function RandomLen_ZZ(l: number, stream: RandomStream): bigint {
  return RandomLen(l, stream);
}
/** NTL bounded integer sampler; word=true selects the native long overload.
 * @see Deviation: NTL integer sampling stream contexts
 */
export function RandomBnd(
  bound: bigint,
  stream: RandomStream,
  options?: { word?: boolean }
): bigint {
  if (options?.word && (bound < -(1n << 63n) || bound >= 1n << 63n))
    throw new RangeError('RandomBnd: word bound must fit a signed 64-bit integer');
  if (bound <= 1n) return 0n;
  const l = (options?.word ? bound - 1n : bound).toString(2).length,
    nb = Math.ceil(l / 8);
  if (options?.word || nb <= 3) {
    const mask = (1n << BigInt(l)) - 1n;
    let out: bigint;
    do {
      out = WordFromBytes(stream.get(nb)) & mask;
    } while (out >= bound);
    return out;
  }
  const highBound = bound >> BigInt(8 * (nb - 2)),
    mask = (1n << BigInt(16 - nb * 8 + l)) - 1n,
    bytes = new Uint8Array(nb);
  for (;;) {
    const high = WordFromBytes(stream.get(2)) & mask;
    if (high > highBound) continue;
    bytes.set(stream.get(nb - 2));
    bytes[nb - 2] = Number(high & 255n);
    bytes[nb - 1] = Number(high >> 8n);
    const out = integerFromRandomBytes(bytes);
    if (high < highBound || out < bound) return out;
  }
}

/** NTL cached word-bound sampler, retaining its explicit stream reference.
 * @see Deviation: NTL cached word samplers and root products
 */
export class RandomBndGenerator {
  p = 0n;
  nb = 0;
  mask = 0n;
  private str: RandomStream;
  constructor(other: RandomBndGenerator);
  constructor(bound: bigint | null, stream: RandomStream);
  constructor(bound: bigint | null | RandomBndGenerator, stream?: RandomStream) {
    if (bound instanceof RandomBndGenerator) {
      this.p = bound.p;
      this.nb = bound.nb;
      this.mask = bound.mask;
      this.str = bound.str;
      return;
    }
    this.str = stream!;
    if (bound !== null) this.build(bound);
  }
  /** Native copy assignment; aliases observe the copied bound and stream reference. */
  assign(other: RandomBndGenerator): this {
    this.p = other.p;
    this.nb = other.nb;
    this.mask = other.mask;
    this.str = other.str;
    return this;
  }
  build(bound: bigint): void {
    if (bound < -(1n << 63n) || bound >= 1n << 63n)
      throw new RangeError('RandomBnd: word bound must fit a signed 64-bit integer');
    if (bound <= 1n) throw new Error('RandomBndGenerator::init: bad args');
    this.p = bound;
    const l = (bound - 1n).toString(2).length;
    this.nb = Math.ceil(l / 8);
    this.mask = (1n << BigInt(l)) - 1n;
  }
  next(): bigint {
    if (this.p === 0n) throw new Error('RandomBndGenerator::next: uninitialized generator');
    let out: bigint;
    do {
      out = WordFromBytes(this.str.get(this.nb)) & this.mask;
    } while (out >= this.p);
    return out;
  }
}
/** NTL vector of bounded signed-word samples, with a cached bound.
 * @see Deviation: NTL cached word samplers and root products
 */
export function VectorRandomBnd(k: number, n: bigint, stream: RandomStream): bigint[] {
  if (n < -(1n << 63n) || n >= 1n << 63n)
    throw new RangeError('RandomBnd: word bound must fit a signed 64-bit integer');
  if (k <= 0) return [];
  const out = Array<bigint>(k).fill(0n);
  if (n <= 1n) return out;
  const gen = new RandomBndGenerator(n, stream);
  for (let i = 0; i < k; i++) out[i] = gen.next();
  return out;
}

// Original ZZ.cpp403-791,983-998,3750-3867, with explicit stream contexts.
const primeBits = (x: bigint): number => (x === 0n ? 0 : (x < 0n ? -x : x).toString(2).length);
function primePower(a: bigint, e: bigint, n: bigint): bigint {
  let result = 1n;
  while (e) {
    if (e & 1n) result = (result * a) % n;
    a = (a * a) % n;
    e >>= 1n;
  }
  return result;
}
function primeWitness(n: bigint, x: bigint): number {
  if (x === 0n) return 0;
  let m = n - 1n,
    k = 0;
  while (!(m & 1n)) {
    m >>= 1n;
    k++;
  }
  let z = primePower(x, m, n);
  if (z === 1n) return 0;
  let y: bigint,
    j = 0;
  do {
    y = z;
    z = (y * y) % n;
    j++;
  } while (j !== k && z !== 1n);
  return Number(z !== 1n || y !== n - 1n);
}
/** Native trial-division tuning bound, with 64-bit GMP limbs and 60-bit word primes.
 * @see Deviation: NTL probable-prime stream contexts
 */
export function ComputePrimeBound(bn: number): bigint {
  const b = BigInt(bn),
    wn = (b + 63n) / 64n;
  const fn = wn <= 36n ? wn / 4n + 1n : BigInt(Math.trunc(1.67 * Math.sqrt(Number(wn))));
  return primeBits(b) + primeBits(fn) > 60 ? 1n << 60n : b * fn;
}
function wordProbPrime(n: bigint, trials: number, stream: RandomStream): number {
  if (n <= 1n) return 0;
  for (const p of [2n, 3n, 5n, 7n, 11n, 13n]) {
    if (n === p) return 1;
    if (n % p === 0n) return 0;
  }
  if (n >= 1n << 60n) return ProbPrime(n, stream, { NumTrials: trials });
  for (let i = 0; i < trials + 1; i++) {
    let x = 2n;
    if (i !== 0)
      do {
        x = RandomBnd(n, stream, { word: true });
      } while (x === 0n);
    if (primeWitness(n, x)) return 0;
  }
  return 1;
}
/** Native probable-prime test with explicit stream context.
 * @see Deviation: NTL probable-prime stream contexts
 */
export function ProbPrime(
  n: bigint,
  stream: RandomStream,
  options?: { NumTrials?: number; word?: boolean }
): number;
export function ProbPrime(n: ZZ, NumTrials?: number): number;
export function ProbPrime(
  n: bigint | ZZ,
  context?: RandomStream | number,
  options?: { NumTrials?: number; word?: boolean }
): number {
  if (typeof n !== 'bigint') return ZZ.ProbPrime(n, context as number | undefined);
  const stream = context as RandomStream;
  if (options?.word && (n < -(1n << 63n) || n >= 1n << 63n))
    throw new RangeError('ProbPrime: word input must fit a signed 64-bit integer');
  const trials = Math.max(0, options?.NumTrials ?? 10);
  if (options?.word) return wordProbPrime(n, trials, stream);
  if (n <= 1n) return 0;
  if (n < 1n << 60n) return wordProbPrime(n, trials, stream);
  const bound = ComputePrimeBound(primeBits(n)),
    s = new PrimeSeq();
  for (let p = s.next(); p && p < bound; p = s.next()) if (n % p === 0n) return 0;
  if (primeWitness(n, 2n)) return 0;
  for (let i = 0; i < trials; i++) {
    let x: bigint;
    do {
      x = RandomBnd(n, stream);
    } while (x === 0n);
    if (primeWitness(n, x)) return 0;
  }
  return 1;
}
function primeSeedStream(seed: bigint): RandomStream {
  const bytes: number[] = [];
  for (let s = seed < 0n ? -seed : seed; s; s >>= 8n) bytes.push(Number(s & 255n));
  return new RandomStream(DeriveKey(Uint8Array.from(bytes), 32));
}
// NTL's parallel candidate schedule, with the native oracle fixed to one thread.
// Candidate blocks and the winning counter do not depend on worker scheduling.
function multiThreadedRandomPrime(length: number, trials: number, stream: RandomStream): bigint {
  let seed = RandomBits(256, stream),
    initial = 0n;
  for (;;) {
    const local = primeSeedStream(seed);
    let counter = initial,
      candidate = 0n,
      winning = 0n,
      found = false;
    while (counter < 1n << 63n) {
      const current = counter++;
      local.set_nonce(current);
      for (let iter = 0; iter < 8; iter++) {
        let n = RandomLen(length, local);
        if (!(n & 1n)) n++;
        if (ProbPrime(n, local, { NumTrials: 0 })) {
          candidate = n;
          winning = current;
          found = true;
          break;
        }
      }
      if (found) break;
    }
    if (!found) {
      initial = 0n;
      seed = RandomBits(256, stream);
      continue;
    }
    const witnesses: bigint[] = [];
    for (let i = 0; i < trials; i++) {
      let x: bigint;
      do {
        x = RandomBnd(candidate, stream);
      } while (x === 0n);
      witnesses.push(x);
    }
    let passed = true;
    for (const x of witnesses)
      if (primeWitness(candidate, x)) {
        passed = false;
        break;
      }
    if (passed) return candidate;
    initial = winning + 1n;
  }
}
/** Native random prime, including its candidate-block schedule from 256 bits.
 * @see Deviation: NTL probable-prime stream contexts
 */
export function RandomPrime(
  length: number,
  stream: RandomStream,
  options?: { NumTrials?: number }
): bigint;
export function RandomPrime(length: number, NumTrials?: number): ZZ;
export function RandomPrime(
  length: number,
  context?: RandomStream | number,
  options?: { NumTrials?: number }
): bigint | ZZ {
  if (!(context instanceof RandomStream)) return ZZ.RandomPrime(length, context);
  const stream = context;
  const trials = Math.max(0, options?.NumTrials ?? 10);
  if (length >= 256) return multiThreadedRandomPrime(length, trials, stream);
  if (length <= 1) throw new Error('RandomPrime: l out of range');
  if (length === 2) return RandomBnd(2n, stream, { word: true }) ? 3n : 2n;
  for (;;) {
    let n = RandomLen(length, stream);
    if (!(n & 1n)) n++;
    if (ProbPrime(n, stream, { NumTrials: trials })) return n;
  }
}
/** Native historical sequential prime generator.
 * @see Deviation: NTL probable-prime stream contexts
 */
export function OldRandomPrime(
  length: number,
  stream: RandomStream,
  options?: { NumTrials?: number }
): bigint {
  if (length <= 1) throw new Error('RandomPrime: l out of range');
  if (length === 2) return RandomBnd(2n, stream, { word: true }) ? 3n : 2n;
  for (;;) {
    let n = RandomLen(length, stream);
    if (!(n & 1n)) n++;
    if (ProbPrime(n, stream, options)) return n;
  }
}
/** Native signed-word random prime; candidates are not forced odd.
 * @see Deviation: NTL probable-prime stream contexts
 */
export function RandomPrime_long(
  length: number,
  stream: RandomStream,
  options?: { NumTrials?: number }
): bigint {
  const trials = Math.max(0, options?.NumTrials ?? 10);
  if (length <= 1 || length >= 64) throw new Error('RandomPrime: length out of range');
  for (;;) {
    const n = RandomLen_long(length, stream);
    if (wordProbPrime(n, trials, stream)) return n;
  }
}
/** Native Damgard–Landrock–Pomerance bound selector, preserving all return codes.
 * @see Deviation: NTL probable-prime stream contexts
 */
export function ErrBoundTest(kk: number, tt: number, nn: number): number {
  const fudge = 1 + 1024 / 2 ** 52,
    log2 = (x: number) => Math.log(x) / Math.log(2);
  const k = kk,
    t = tt,
    n = nn;
  if (k < 3 || t < 1) return 0;
  if (n < 1) return 1;
  if (9 * t > 2 ** 52) throw new Error('ErrBoundTest: t too big');
  const lk = log2(k);
  if ((n + lk) * fudge <= 2 * t) return 1;
  if ((2 * lk + 4 + n) * fudge <= 2 * Math.sqrt(k)) return 2;
  if ((t === 2 && k >= 88) || (3 <= t && 9 * t <= k && k >= 21)) {
    if ((1.5 * lk + t + 4 + n) * fudge <= 0.5 * log2(t) + 2 * Math.sqrt(t * k)) return 3;
  }
  if (k <= 9 * t && 4 * t <= k && k >= 21) {
    if (
      (log2(3) + log2(7) + lk + n) * fudge <= log2(20) + 5 * t &&
      (log2(3) + (15 / 4) * lk + n) * fudge <= log2(7) + k / 2 + 2 * t &&
      (2 * log2(3) + 2 + lk + n) * fudge <= k / 4 + 3 * t
    )
      return 4;
  }
  if (4 * t >= k && k >= 21) if (((15 / 4) * lk + n) * fudge <= log2(7) + k / 2 + 2 * t) return 5;
  return 0;
}
/** Native error-bound-controlled integer prime generation.
 * @see Deviation: NTL probable-prime stream contexts
 */
export function GenPrime(length: number, stream: RandomStream, options?: { err?: number }): bigint {
  if (length <= 1) throw new Error('GenPrime: bad length');
  if (length > 2 ** 20) throw new Error('GenPrime: length too large');
  const err = Math.max(1, Math.min(512, options?.err ?? 80));
  if (length === 2) return RandomBnd(2n, stream, { word: true }) ? 3n : 2n;
  let trials = 1;
  while (!ErrBoundTest(length, trials, err)) trials++;
  return RandomPrime(length, stream, { NumTrials: trials });
}
/** Native error-bound-controlled signed-word prime generation.
 * @see Deviation: NTL probable-prime stream contexts
 */
export function GenPrime_long(
  length: number,
  stream: RandomStream,
  options?: { err?: number }
): bigint {
  if (length <= 1) throw new Error('GenPrime: bad length');
  if (length >= 64) throw new Error('GenPrime: length too large');
  const err = Math.max(1, Math.min(512, options?.err ?? 80));
  if (length === 2) return RandomBnd(2n, stream, { word: true }) ? 3n : 2n;
  let trials = 1;
  while (!ErrBoundTest(length, trials, err)) trials++;
  return RandomPrime_long(length, stream, { NumTrials: trials });
}

/** NTL's signed, positive-half CRT interval predicate.
 * @see Deviation: NTL integer CRT and modular determinant adapters
 */
export function CRTInRange(g: bigint, a: bigint): number {
  return a > 0n && -a < 2n * g && 2n * g <= a ? 1 : 0;
}
/** Incremental native CRT; returns [modified, residue, accumulated modulus].
 * The word overload uses the same mathematical schedule, delegating above 2^60.
 * @see Deviation: NTL integer CRT and modular determinant adapters
 */
export function CRT(
  gg: bigint,
  a: bigint,
  G: bigint,
  p: bigint,
  options: { word?: boolean } = {}
): [number, bigint, bigint] {
  if (options.word && (p < -(1n << 63n) || p >= 1n << 63n || G < -(1n << 63n) || G >= 1n << 63n))
    throw new RangeError('CRT: word inputs must fit signed 64-bit integers');
  if (a <= 0n || p <= 1n || G < 0n || G >= p)
    throw new RangeError('CRT: require a > 0, p > 1 and 0 <= G < p');
  const norm = (x: bigint, m: bigint) => ((x % m) + m) % m;
  let modified = 0,
    g = gg;
  if (!CRTInRange(g, a)) {
    modified = 1;
    g = norm(g, a);
    if (g > a / 2n) g -= a;
  }
  let u = norm(a, p),
    v = p,
    s = 1n,
    t = 0n;
  while (v) {
    const q = u / v;
    [u, v] = [v, u - q * v];
    [s, t] = [t, s - q * t];
  }
  if (u !== 1n) throw new Error('InvMod: inverse undefined');
  let h = norm((G - norm(g, p)) * norm(s, p), p);
  if (h > p / 2n) h -= p;
  if (h) {
    modified = 1;
    if (!(p & 1n) && g > 0n && h === p / 2n) g -= a * h;
    else g += a * h;
  }
  return [modified, g, a * p];
}
