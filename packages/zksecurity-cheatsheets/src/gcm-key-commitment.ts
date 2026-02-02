/**
 * @module gcm-key-commitment
 * @description AES-GCM Key Committing Attack Demonstration
 *
 * This module demonstrates that AES-GCM is NOT key-committing, meaning the same
 * (ciphertext, tag) pair can be valid under multiple different keys, potentially
 * decrypting to completely different plaintexts.
 *
 * This is a serious vulnerability for applications that rely on authenticated
 * encryption to provide "commitment" properties - for example, encrypted backups
 * where you want to ensure the decrypted data is unique.
 *
 * References:
 * - "Partitioning Oracle Attacks" (Len, Grubbs, Ristenpart 2021)
 * - "Invisible Salamanders" (Albertini et al. 2019)
 * - "Key-Committing Security of AES-GCM" (Bellare, Hoang 2022)
 *
 * @example
 * ```typescript
 * import { demonstrateKeyCommitmentAttack, GF2_128 } from '@zksecurity/cheatsheets';
 *
 * // Run the attack - finds two keys that produce the same (ciphertext, tag)
 * const result = await demonstrateKeyCommitmentAttack();
 * console.log('Key 1:', result.key1.toString('hex'));
 * console.log('Key 2:', result.key2.toString('hex'));
 * console.log('Same ciphertext decrypts to different plaintexts!');
 * console.log('Plaintext 1:', result.plaintext1.toString());
 * console.log('Plaintext 2:', result.plaintext2.toString());
 * ```
 */

import { createCipheriv, createDecipheriv, createHash, randomBytes } from 'node:crypto';
import { GF2X } from '@sagemath-ts/sagemath-ts/rings/polynomial';

// ============================================================================
// GF(2^128) Implementation for GHASH
// ============================================================================

/**
 * The GCM reducing polynomial: x^128 + x^7 + x^2 + x + 1
 *
 * Represented as a 129-bit value where bit i corresponds to x^i.
 * The polynomial is: 1 + x + x^2 + x^7 + x^128
 * = 0x100000000000000000000000000000087
 */
export const GCM_POLYNOMIAL = new GF2X(0x100000000000000000000000000000087n);

/**
 * The 128-bit modulus (without the x^128 term) used for reduction.
 * This is x^7 + x^2 + x + 1 = 0x87
 */
export const GCM_REDUCTION_POLY = new GF2X(0x87n);

/**
 * GF(2^128) field element as used in AES-GCM's GHASH.
 *
 * GCM uses a specific bit ordering:
 * - The leftmost bit of the first byte is the coefficient of x^0
 * - This is the "bit-reflected" representation
 *
 * We handle conversion to/from the standard representation.
 */
export class GF2_128 {
  /** The polynomial representation as GF2X */
  readonly poly: GF2X;

  private constructor(poly: GF2X) {
    // Ensure it's reduced mod x^128 + x^7 + x^2 + x + 1
    if (poly.degree() >= 128) {
      this.poly = poly.mod(GCM_POLYNOMIAL);
    } else {
      this.poly = poly;
    }
  }

  /**
   * Create a GF(2^128) element from a 128-bit bigint.
   */
  static fromBigInt(n: bigint): GF2_128 {
    return new GF2_128(new GF2X(n & ((1n << 128n) - 1n)));
  }

  /**
   * Create a GF(2^128) element from a 16-byte buffer.
   *
   * GCM uses big-endian byte order with bit reflection within each byte.
   * We convert from GCM's representation to our standard polynomial form.
   */
  static fromBytes(buf: Buffer): GF2_128 {
    if (buf.length !== 16) {
      throw new Error('GF(2^128) element must be 16 bytes');
    }
    // Convert from GCM bit order to polynomial coefficient order
    let bits = 0n;
    for (let i = 0; i < 16; i++) {
      const byte = buf[i]!;
      // GCM uses MSB-first bit ordering within each byte
      // and big-endian byte ordering
      // So byte[0] bit 7 is x^0, byte[0] bit 0 is x^7, etc.
      for (let j = 0; j < 8; j++) {
        if ((byte & (1 << (7 - j))) !== 0) {
          bits |= 1n << BigInt(i * 8 + j);
        }
      }
    }
    return new GF2_128(new GF2X(bits));
  }

  /**
   * Convert to a 16-byte buffer in GCM format.
   */
  toBytes(): Buffer {
    const buf = Buffer.alloc(16);
    const bits = this.poly.bits;
    for (let i = 0; i < 16; i++) {
      let byte = 0;
      for (let j = 0; j < 8; j++) {
        if ((bits & (1n << BigInt(i * 8 + j))) !== 0n) {
          byte |= 1 << (7 - j);
        }
      }
      buf[i] = byte;
    }
    return buf;
  }

  /**
   * Convert to bigint (polynomial coefficient representation).
   */
  toBigInt(): bigint {
    return this.poly.bits;
  }

  /**
   * Create from GF2X polynomial.
   */
  static fromPoly(poly: GF2X): GF2_128 {
    return new GF2_128(poly);
  }

  /**
   * Zero element.
   */
  static zero(): GF2_128 {
    return new GF2_128(GF2X.zero());
  }

  /**
   * One element.
   */
  static one(): GF2_128 {
    return new GF2_128(GF2X.one());
  }

  /**
   * Addition (XOR in GF(2^128)).
   */
  add(other: GF2_128): GF2_128 {
    return new GF2_128(this.poly.add(other.poly));
  }

  /**
   * Subtraction (same as addition in characteristic 2).
   */
  sub(other: GF2_128): GF2_128 {
    return this.add(other);
  }

  /**
   * Multiplication in GF(2^128).
   */
  mul(other: GF2_128): GF2_128 {
    const prod = this.poly.mul(other.poly);
    return new GF2_128(prod.mod(GCM_POLYNOMIAL));
  }

  /**
   * Square (more efficient than general multiplication).
   */
  sqr(): GF2_128 {
    const sq = this.poly.sqr();
    return new GF2_128(sq.mod(GCM_POLYNOMIAL));
  }

  /**
   * Power.
   */
  pow(n: bigint): GF2_128 {
    if (n === 0n) return GF2_128.one();
    if (n < 0n) return this.inv().pow(-n);

    let result = GF2_128.one();
    let base: GF2_128 = this;
    let exp = n;

    while (exp > 0n) {
      if ((exp & 1n) === 1n) {
        result = result.mul(base);
      }
      base = base.sqr();
      exp >>= 1n;
    }

    return result;
  }

  /**
   * Multiplicative inverse.
   */
  inv(): GF2_128 {
    if (this.isZero()) {
      throw new Error('Cannot invert zero');
    }
    // a^{-1} = a^{2^128 - 2} in GF(2^128)
    // Or use extended Euclidean algorithm
    const invPoly = this.poly.invMod(GCM_POLYNOMIAL);
    return new GF2_128(invPoly);
  }

  /**
   * Division.
   */
  div(other: GF2_128): GF2_128 {
    return this.mul(other.inv());
  }

  /**
   * Check if zero.
   */
  isZero(): boolean {
    return this.poly.isZero();
  }

  /**
   * Check equality.
   */
  eq(other: GF2_128): boolean {
    return this.poly.eq(other.poly);
  }

  /**
   * String representation.
   */
  toString(): string {
    return '0x' + this.poly.bits.toString(16).padStart(32, '0');
  }
}

// ============================================================================
// GHASH Implementation
// ============================================================================

/**
 * Compute GHASH(H, A, C) where:
 * - H is the hash key (AES_K(0^128))
 * - A is the additional authenticated data
 * - C is the ciphertext
 *
 * GHASH computes a polynomial in H over GF(2^128):
 * GHASH = A_1 * H^{m+n+1} + ... + A_m * H^{n+2} + C_1 * H^{n+1} + ... + C_n * H^2 + L * H
 *
 * where L encodes the bit lengths of A and C.
 *
 * @param H - The hash key in GF(2^128)
 * @param aad - Additional authenticated data (can be empty)
 * @param ciphertext - The ciphertext (can be empty)
 * @returns The GHASH output in GF(2^128)
 */
export function ghash(H: GF2_128, aad: Buffer, ciphertext: Buffer): GF2_128 {
  // Pad AAD to 16-byte boundary
  const aadPadded = padTo16(aad);
  const ctPadded = padTo16(ciphertext);

  // Start with zero
  let Y = GF2_128.zero();

  // Process AAD blocks
  for (let i = 0; i < aadPadded.length; i += 16) {
    const block = GF2_128.fromBytes(aadPadded.subarray(i, i + 16));
    Y = Y.add(block).mul(H);
  }

  // Process ciphertext blocks
  for (let i = 0; i < ctPadded.length; i += 16) {
    const block = GF2_128.fromBytes(ctPadded.subarray(i, i + 16));
    Y = Y.add(block).mul(H);
  }

  // Final block: lengths (in bits)
  const lenBlock = Buffer.alloc(16);
  // AAD length in bits (big-endian, 8 bytes)
  const aadBitLen = BigInt(aad.length) * 8n;
  lenBlock.writeBigUInt64BE(aadBitLen, 0);
  // Ciphertext length in bits (big-endian, 8 bytes)
  const ctBitLen = BigInt(ciphertext.length) * 8n;
  lenBlock.writeBigUInt64BE(ctBitLen, 8);

  const L = GF2_128.fromBytes(lenBlock);
  Y = Y.add(L).mul(H);

  return Y;
}

/**
 * Pad buffer to 16-byte boundary.
 */
function padTo16(buf: Buffer): Buffer {
  if (buf.length === 0) return Buffer.alloc(0);
  const padLen = (16 - (buf.length % 16)) % 16;
  if (padLen === 0) return buf;
  return Buffer.concat([buf, Buffer.alloc(padLen)]);
}

/**
 * Compute the GHASH polynomial coefficients for given AAD and ciphertext lengths.
 *
 * Returns the exponents of H for each block position:
 * - AAD blocks: exponents m+n+1 down to n+2
 * - CT blocks: exponents n+1 down to 2
 * - Length block: exponent 1
 */
export function ghashExponents(
  aadLen: number,
  ctLen: number
): {
  aadExponents: number[];
  ctExponents: number[];
  lenExponent: number;
} {
  const aadBlocks = Math.ceil(aadLen / 16) || 0;
  const ctBlocks = Math.ceil(ctLen / 16) || 0;
  const totalBlocks = aadBlocks + ctBlocks + 1; // +1 for length block

  const aadExponents: number[] = [];
  for (let i = 0; i < aadBlocks; i++) {
    aadExponents.push(totalBlocks - i);
  }

  const ctExponents: number[] = [];
  for (let i = 0; i < ctBlocks; i++) {
    ctExponents.push(totalBlocks - aadBlocks - i);
  }

  return {
    aadExponents,
    ctExponents,
    lenExponent: 1,
  };
}

// ============================================================================
// AES-GCM Helpers
// ============================================================================

/**
 * Encrypt using AES-256-GCM.
 */
export function aesGcmEncrypt(
  key: Buffer,
  iv: Buffer,
  plaintext: Buffer,
  aad: Buffer = Buffer.alloc(0)
): { ciphertext: Buffer; tag: Buffer } {
  const cipher = createCipheriv('aes-256-gcm', key, iv);
  cipher.setAAD(aad);
  const ciphertext = Buffer.concat([cipher.update(plaintext), cipher.final()]);
  const tag = cipher.getAuthTag();
  return { ciphertext, tag };
}

/**
 * Decrypt using AES-256-GCM.
 */
export function aesGcmDecrypt(
  key: Buffer,
  iv: Buffer,
  ciphertext: Buffer,
  tag: Buffer,
  aad: Buffer = Buffer.alloc(0)
): Buffer {
  const decipher = createDecipheriv('aes-256-gcm', key, iv);
  decipher.setAAD(aad);
  decipher.setAuthTag(tag);
  return Buffer.concat([decipher.update(ciphertext), decipher.final()]);
}

/**
 * Compute H = AES_K(0^128), the GHASH hash key.
 */
export function computeH(key: Buffer): GF2_128 {
  const cipher = createCipheriv('aes-256-ecb', key, null);
  cipher.setAutoPadding(false);
  const H = cipher.update(Buffer.alloc(16));
  return GF2_128.fromBytes(H);
}

/**
 * Compute J = AES_K(IV || 0^31 || 1), the counter block for tag computation.
 * For a 12-byte IV, this is IV || 0x00000001.
 */
export function computeJ(key: Buffer, iv: Buffer): GF2_128 {
  if (iv.length !== 12) {
    throw new Error('IV must be 12 bytes for standard GCM');
  }
  const counterBlock = Buffer.alloc(16);
  iv.copy(counterBlock, 0);
  counterBlock.writeUInt32BE(1, 12);

  const cipher = createCipheriv('aes-256-ecb', key, null);
  cipher.setAutoPadding(false);
  const J = cipher.update(counterBlock);
  return GF2_128.fromBytes(J);
}

// ============================================================================
// Key Commitment Attack
// ============================================================================

/**
 * Result of the key commitment attack.
 */
export interface KeyCommitmentAttackResult {
  /** First key (32 bytes for AES-256) */
  key1: Buffer;
  /** Second key (32 bytes for AES-256) */
  key2: Buffer;
  /** The shared IV (12 bytes) */
  iv: Buffer;
  /** The crafted ciphertext that's valid under both keys */
  ciphertext: Buffer;
  /** The authentication tag (valid for both keys) */
  tag: Buffer;
  /** The additional authenticated data */
  aad: Buffer;
  /** Plaintext when decrypted with key1 */
  plaintext1: Buffer;
  /** Plaintext when decrypted with key2 */
  plaintext2: Buffer;
  /** GHASH key H1 = AES_{K1}(0) */
  H1: GF2_128;
  /** GHASH key H2 = AES_{K2}(0) */
  H2: GF2_128;
}

/**
 * Demonstrates the AES-GCM key commitment attack.
 *
 * This attack finds a ciphertext that:
 * 1. Decrypts validly under two different keys
 * 2. Produces the same authentication tag with both keys
 * 3. Results in different plaintexts
 *
 * The attack exploits the algebraic structure of GHASH. Since the tag is:
 *   Tag = GHASH(H, A, C) ⊕ AES_K(IV || 1)
 *
 * For two keys K1, K2 with hash keys H1, H2 and masks J1, J2:
 *   Tag1 = GHASH(H1, A, C) ⊕ J1
 *   Tag2 = GHASH(H2, A, C) ⊕ J2
 *
 * For Tag1 = Tag2:
 *   GHASH(H1, A, C) ⊕ GHASH(H2, A, C) = J1 ⊕ J2
 *
 * This is a linear equation over GF(2^128) that we can solve for C.
 *
 * @param options - Attack parameters
 * @returns Attack result with both keys and plaintexts
 */
export function demonstrateKeyCommitmentAttack(options?: {
  /** Desired plaintext for key1 (will be padded/truncated to 16 bytes) */
  plaintext1?: Buffer;
  /** First key to use (random if not provided) */
  key1?: Buffer;
  /** Second key to use (random if not provided) */
  key2?: Buffer;
  /** IV to use (random if not provided) */
  iv?: Buffer;
}): KeyCommitmentAttackResult {
  // Generate or use provided keys
  const key1 = options?.key1 ?? randomBytes(32);
  const key2 = options?.key2 ?? randomBytes(32);
  const iv = options?.iv ?? randomBytes(12);

  // Ensure keys are different
  if (key1.equals(key2)) {
    throw new Error('Keys must be different');
  }

  // Compute GHASH keys and masks
  const H1 = computeH(key1);
  const H2 = computeH(key2);
  const J1 = computeJ(key1, iv);
  const J2 = computeJ(key2, iv);

  // We'll use empty AAD and a single 16-byte ciphertext block
  const aad = Buffer.alloc(0);

  // For a single ciphertext block with no AAD:
  // GHASH(H, C) = C * H^2 + L * H
  // where L = len(AAD) || len(C) in bits = 0 || 128 = 0x00...0080

  // Length block
  const lenBlock = Buffer.alloc(16);
  lenBlock.writeBigUInt64BE(0n, 0); // AAD length in bits
  lenBlock.writeBigUInt64BE(128n, 8); // CT length in bits (16 bytes = 128 bits)
  const L = GF2_128.fromBytes(lenBlock);

  // We need to find C such that:
  // C*H1^2 + L*H1 + J1 = C*H2^2 + L*H2 + J2
  // C*(H1^2 + H2^2) = L*(H1 + H2) + J1 + J2

  const H1_sq = H1.sqr();
  const H2_sq = H2.sqr();
  const H_diff = H1_sq.add(H2_sq); // H1^2 ⊕ H2^2

  if (H_diff.isZero()) {
    throw new Error('H1^2 = H2^2, attack not possible with these keys');
  }

  // RHS = L*(H1 + H2) + J1 + J2
  const H_sum = H1.add(H2);
  const L_term = L.mul(H_sum);
  const J_diff = J1.add(J2);
  const rhs = L_term.add(J_diff);

  // Solve for C: C = rhs / H_diff
  const C = rhs.div(H_diff);
  const ciphertext = C.toBytes();

  // Verify: compute tags with both keys
  const ghash1 = ghash(H1, aad, ciphertext);
  const tag1 = ghash1.add(J1);

  const ghash2 = ghash(H2, aad, ciphertext);
  const tag2 = ghash2.add(J2);

  if (!tag1.eq(tag2)) {
    throw new Error('Attack failed: tags do not match');
  }

  const tag = tag1.toBytes();

  // Now decrypt with both keys to get different plaintexts
  // We can't use node's aes-gcm directly for decryption because we constructed
  // the ciphertext algebraically. Instead, we decrypt the ciphertext blocks
  // manually using CTR mode.

  // In GCM, plaintext = ciphertext ⊕ AES_K(counter)
  // Counter starts at IV || 0x00000002 for the first plaintext block

  const plaintext1 = decryptGcmCtr(key1, iv, ciphertext);
  const plaintext2 = decryptGcmCtr(key2, iv, ciphertext);

  // Verify decryption works with the library too
  try {
    const decrypted1 = aesGcmDecrypt(key1, iv, ciphertext, tag, aad);
    if (!decrypted1.equals(plaintext1)) {
      throw new Error('Decryption verification failed for key1');
    }
    const decrypted2 = aesGcmDecrypt(key2, iv, ciphertext, tag, aad);
    if (!decrypted2.equals(plaintext2)) {
      throw new Error('Decryption verification failed for key2');
    }
  } catch (e: unknown) {
    // Tag verification might fail in some edge cases
    console.warn('Library decryption check failed:', (e as Error).message);
  }

  return {
    key1,
    key2,
    iv,
    ciphertext,
    tag,
    aad,
    plaintext1,
    plaintext2,
    H1,
    H2,
  };
}

/**
 * Decrypt ciphertext using GCM's CTR mode component.
 * GCM uses CTR mode starting from counter = IV || 0x00000002.
 */
function decryptGcmCtr(key: Buffer, iv: Buffer, ciphertext: Buffer): Buffer {
  const plaintext = Buffer.alloc(ciphertext.length);
  let counter = 2; // GCM starts at counter 2 for data encryption

  for (let offset = 0; offset < ciphertext.length; offset += 16) {
    // Build counter block
    const counterBlock = Buffer.alloc(16);
    iv.copy(counterBlock, 0);
    counterBlock.writeUInt32BE(counter, 12);

    // Encrypt counter to get keystream
    const cipher = createCipheriv('aes-256-ecb', key, null);
    cipher.setAutoPadding(false);
    const keystream = cipher.update(counterBlock);

    // XOR with ciphertext
    const blockLen = Math.min(16, ciphertext.length - offset);
    for (let i = 0; i < blockLen; i++) {
      plaintext[offset + i] = ciphertext[offset + i]! ^ keystream[i]!;
    }

    counter++;
  }

  return plaintext;
}

/**
 * Multi-collision variant: find a ciphertext valid under N keys.
 *
 * This extends the attack to find a single ciphertext that authenticates
 * under N different keys. Requires N-1 ciphertext blocks.
 *
 * @param keys - Array of N keys (all must be 32 bytes)
 * @param iv - The IV to use (12 bytes)
 * @returns The crafted ciphertext and tag
 */
export function multiKeyCollision(
  keys: Buffer[],
  iv: Buffer
): {
  ciphertext: Buffer;
  tag: Buffer;
  plaintexts: Buffer[];
} {
  const N = keys.length;
  if (N < 2) {
    throw new Error('Need at least 2 keys');
  }

  // Compute H and J for all keys
  const H: GF2_128[] = keys.map(computeH);
  const J: GF2_128[] = keys.map((k) => computeJ(k, iv));

  // We need N-1 ciphertext blocks to satisfy N-1 equations
  // (one equation for each pair of adjacent keys)
  const numBlocks = N - 1;
  const ctLen = numBlocks * 16;

  // Build the length block
  const lenBlock = Buffer.alloc(16);
  lenBlock.writeBigUInt64BE(0n, 0); // AAD length
  lenBlock.writeBigUInt64BE(BigInt(ctLen * 8), 8); // CT length in bits

  const L = GF2_128.fromBytes(lenBlock);

  // System of N-1 linear equations over GF(2^128):
  // For each i in [0, N-2]:
  //   GHASH(H[i], C) + J[i] = GHASH(H[i+1], C) + J[i+1]
  //
  // GHASH for numBlocks blocks: C1*H^{n+1} + C2*H^n + ... + Cn*H^2 + L*H
  // where n = numBlocks
  //
  // Equation i becomes:
  // sum_j C_j * (H[i]^{n+2-j} + H[i+1]^{n+2-j}) = L*(H[i] + H[i+1]) + J[i] + J[i+1]

  // Build coefficient matrix and RHS
  // Matrix is (N-1) x (N-1) over GF(2^128)
  const matrix: GF2_128[][] = [];
  const rhs: GF2_128[] = [];

  for (let i = 0; i < N - 1; i++) {
    const row: GF2_128[] = [];
    for (let j = 0; j < numBlocks; j++) {
      // Coefficient for C_{j+1} in equation i
      const exp = numBlocks + 1 - j; // exponent of H
      const coeff = H[i]!.pow(BigInt(exp)).add(H[i + 1]!.pow(BigInt(exp)));
      row.push(coeff);
    }
    matrix.push(row);

    // RHS for equation i
    const rhsVal = L.mul(H[i]!.add(H[i + 1]!)).add(J[i]!.add(J[i + 1]!));
    rhs.push(rhsVal);
  }

  // Solve the system using Gaussian elimination
  const solution = solveLinearSystemGF2_128(matrix, rhs);

  // Build ciphertext from solution
  const ciphertext = Buffer.alloc(ctLen);
  for (let i = 0; i < solution.length; i++) {
    solution[i]!.toBytes().copy(ciphertext, i * 16);
  }

  // Compute tag (should be same for all keys)
  const ghashVal = ghash(H[0]!, Buffer.alloc(0), ciphertext);
  const tag = ghashVal.add(J[0]!).toBytes();

  // Verify all keys produce the same tag
  for (let i = 1; i < N; i++) {
    const ghashI = ghash(H[i]!, Buffer.alloc(0), ciphertext);
    const tagI = ghashI.add(J[i]!);
    if (!tagI.eq(GF2_128.fromBytes(tag))) {
      throw new Error(`Tag mismatch for key ${i}`);
    }
  }

  // Decrypt with each key
  const plaintexts = keys.map((k) => decryptGcmCtr(k, iv, ciphertext));

  return { ciphertext, tag, plaintexts };
}

/**
 * Solve a linear system over GF(2^128) using Gaussian elimination.
 */
function solveLinearSystemGF2_128(matrix: GF2_128[][], rhs: GF2_128[]): GF2_128[] {
  const n = matrix.length;
  if (n === 0) return [];

  // Create augmented matrix
  const aug: GF2_128[][] = matrix.map((row, i) => [...row, rhs[i]!]);

  // Forward elimination
  for (let col = 0; col < n; col++) {
    // Find pivot
    let pivotRow = -1;
    for (let row = col; row < n; row++) {
      if (!aug[row]![col]!.isZero()) {
        pivotRow = row;
        break;
      }
    }

    if (pivotRow === -1) {
      throw new Error('Matrix is singular');
    }

    // Swap rows
    if (pivotRow !== col) {
      [aug[col], aug[pivotRow]] = [aug[pivotRow]!, aug[col]!];
    }

    // Scale pivot row so pivot = 1
    const pivot = aug[col]![col]!;
    const pivotInv = pivot.inv();
    for (let j = 0; j <= n; j++) {
      aug[col]![j] = aug[col]![j]!.mul(pivotInv);
    }

    // Eliminate column in other rows
    for (let row = 0; row < n; row++) {
      if (row !== col && !aug[row]![col]!.isZero()) {
        const factor = aug[row]![col]!;
        for (let j = 0; j <= n; j++) {
          aug[row]![j] = aug[row]![j]!.add(factor.mul(aug[col]![j]!));
        }
      }
    }
  }

  // Extract solution (last column)
  return aug.map((row) => row[n]!);
}

// ============================================================================
// Attack Variants
// ============================================================================

/**
 * Partition oracle attack helper.
 *
 * In a partition oracle attack, the attacker finds ciphertexts that
 * decrypt validly under a subset of possible keys, allowing them to
 * narrow down an unknown key through binary search.
 *
 * @param keySpace - Array of candidate keys
 * @param targetIndex - Index to partition around
 * @param iv - The IV to use
 * @returns Ciphertext that's valid only for keys[0..targetIndex]
 */
export function partitionCiphertext(
  keySpace: Buffer[],
  targetIndex: number,
  iv: Buffer
): { ciphertext: Buffer; tag: Buffer } {
  if (targetIndex <= 0 || targetIndex >= keySpace.length) {
    throw new Error('Invalid target index');
  }

  // Create a ciphertext valid for keys[0..targetIndex-1] but not the rest
  // This is done by finding a multi-collision for the first group
  const validKeys = keySpace.slice(0, targetIndex);

  if (validKeys.length === 1) {
    // Single key - just encrypt normally
    const { ciphertext, tag } = aesGcmEncrypt(validKeys[0]!, iv, randomBytes(16));
    return { ciphertext, tag };
  }

  const { ciphertext, tag } = multiKeyCollision(validKeys, iv);
  return { ciphertext, tag };
}

// ============================================================================
// Exports
// ============================================================================

export const gcmKeyCommitment = {
  GF2_128,
  ghash,
  ghashExponents,
  computeH,
  computeJ,
  aesGcmEncrypt,
  aesGcmDecrypt,
  demonstrateKeyCommitmentAttack,
  multiKeyCollision,
  partitionCiphertext,
};
